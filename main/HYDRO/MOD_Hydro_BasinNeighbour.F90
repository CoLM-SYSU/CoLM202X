#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_Hydro_BasinNeighbour
   !--------------------------------------------------------------------------------!
   ! DESCRIPTION:                                                                   !
   !                                                                                !
   !    Basin Neighbours : data and communication subroutines.                      !
   !                                                                                !
   ! Created by Shupeng Zhang, May 2023                                             !
   !--------------------------------------------------------------------------------!

   USE MOD_Precision
   USE MOD_DataType
   IMPLICIT NONE
   
   ! -- neighbour parameters --
   type basin_neighbour_type
      integer  :: nnb     ! number of neighbours
      real(r8) :: myarea  ! area of this basin [m^2]
      real(r8) :: myelva  ! elevation of this basin [m]
      integer , allocatable :: bindex  (:)   ! neighbour global index
      integer , allocatable :: addr    (:,:) ! data address: (1,:) refers to process, (2,:) refers to location
      real(r8), allocatable :: dist    (:)   ! distance between basin centers [m]
      real(r8), allocatable :: lenbdr  (:)   ! length of boundary line [m]
      real(r8), allocatable :: area    (:)   ! area of neighbours [m^2]
      real(r8), allocatable :: elva    (:)   ! elevation of neighbours [m]
      real(r8), allocatable :: slope   (:)   ! slope (positive) [-]
      logical , allocatable :: iswatb  (:)   ! whether a neighbour is water body
   END type basin_neighbour_type

   type(basin_neighbour_type), allocatable :: basinneighbour (:)

   ! -- neighbour variables --
   TYPE(pointer_real8_1d), allocatable :: theta_a_nb (:)  ! saturated volume content [-]
   TYPE(pointer_real8_1d), allocatable :: zwt_nb     (:)  ! water table depth [m]
   TYPE(pointer_real8_1d), allocatable :: Ks_nb      (:)  ! saturated hydraulic conductivity [m/s]
   TYPE(pointer_real8_1d), allocatable :: wdsrf_nb   (:)  ! depth of surface water [m]

   ! -- neighbour communication --
   TYPE neighbour_sendrecv_type
      INTEGER :: ndata
      INTEGER, allocatable :: bindex (:)
      INTEGER, allocatable :: ibasin (:)
   END TYPE neighbour_sendrecv_type

   TYPE(neighbour_sendrecv_type), allocatable :: recvaddr(:)
   TYPE(neighbour_sendrecv_type), allocatable :: sendaddr(:)

CONTAINS
   
   ! ----------
   SUBROUTINE basin_neighbour_init ()

      USE MOD_SPMD_Task
      USE MOD_Namelist
      USE MOD_NetCDFSerial
      USE MOD_Mesh
      USE MOD_LandElm
      USE MOD_Hydro_HillslopeNetwork
      USE MOD_Hydro_RiverLakeNetwork
      USE MOD_Utils
      IMPLICIT NONE

      ! Local Variables
      CHARACTER(len=256) :: neighbour_file

      INTEGER :: numbasin, ibasin
      INTEGER :: iwork, mesg(2), isrc, idest
      INTEGER :: nrecv, irecv
      INTEGER :: iloc, iloc1, iloc2
      INTEGER :: nnb, nnbinq, inb, ndata
   
      INTEGER :: maxnnb
      INTEGER , allocatable :: nnball   (:)
      INTEGER , allocatable :: idxnball (:,:)
      REAL(r8), allocatable :: lenbdall (:,:)

      INTEGER , allocatable :: addrbasin(:)

      INTEGER , allocatable :: bindex  (:)
      INTEGER , allocatable :: icache1 (:)
      INTEGER , allocatable :: icache2 (:,:)
      REAL(r8), allocatable :: rcache2 (:,:)

      INTEGER, allocatable :: basin_sorted(:), order(:)
      INTEGER, allocatable :: idxinq(:), addrinq(:)

      LOGICAL, allocatable :: mask(:)

      REAL(r8), allocatable :: rlon_b(:), rlat_b(:)
      TYPE(pointer_real8_1d), allocatable :: rlon_nb(:), rlat_nb(:)
   
      REAL(r8), allocatable :: area_b(:)
      REAL(r8), allocatable :: elva_b(:)
      real(r8), allocatable :: iswatb(:)

      TYPE(pointer_real8_1d), allocatable :: area_nb  (:)  ! m^2
      TYPE(pointer_real8_1d), allocatable :: elva_nb  (:)  ! m
      TYPE(pointer_real8_1d), allocatable :: iswat_nb (:)  

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      numbasin = numelm

      neighbour_file = DEF_CatchmentMesh_data 

      IF (p_is_master) THEN
         CALL ncio_read_serial (neighbour_file, 'basin_num_neighbour', nnball  )
         CALL ncio_read_serial (neighbour_file, 'basin_idx_neighbour', idxnball)
         CALL ncio_read_serial (neighbour_file, 'basin_len_border'   , lenbdall)

         maxnnb = size(idxnball,1)

         lenbdall = lenbdall * 1.e3 ! km to m
      ENDIF

#ifdef USEMPI

      CALL mpi_bcast (maxnnb, 1, MPI_INTEGER, p_root, p_comm_glb, p_err)

      IF (p_is_master) THEN

         allocate (addrbasin (size(nnball)))
         addrbasin(:) = -1

         DO iwork = 0, p_np_worker-1

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN
               
               allocate (bindex  (nrecv))
               allocate (icache1 (nrecv))
               allocate (icache2 (maxnnb,nrecv))
               allocate (rcache2 (maxnnb,nrecv))

               CALL mpi_recv (bindex, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               addrbasin(bindex) = isrc
               
               idest = isrc

               icache1 = nnball(bindex)
               CALL mpi_send (icache1, nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  icache2(:,irecv) = idxnball(:,bindex(irecv))
               ENDDO
               CALL mpi_send (icache2, maxnnb*nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache2(:,irecv) = lenbdall(:,bindex(irecv))
               ENDDO
               CALL mpi_send (rcache2, maxnnb*nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               deallocate (bindex )
               deallocate (icache1)
               deallocate (icache2)
               deallocate (rcache2)

            ENDIF
         ENDDO
      ENDIF
#endif

      IF (p_is_worker) THEN
               
         IF (numbasin > 0) THEN
            allocate (bindex (numbasin))
            bindex = landelm%eindex
         ENDIF 

#ifdef USEMPI
         mesg(1:2) = (/p_iam_glb, numbasin/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err) 

         IF (numbasin > 0) THEN
            CALL mpi_send (bindex, numbasin, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_err) 

            allocate (nnball (numbasin))
            CALL mpi_recv (nnball, numbasin, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (idxnball (maxnnb,numbasin))
            CALL mpi_recv (idxnball, maxnnb*numbasin, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (lenbdall (maxnnb,numbasin))
            CALL mpi_recv (lenbdall, maxnnb*numbasin, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
#else
         allocate (icache1 (numbasin))
         allocate (icache2 (maxnnb,numbasin))
         allocate (rcache2 (maxnnb,numbasin))

         icache1 = nnball  
         icache2 = idxnball
         rcache2 = lenbdall

         DO ibasin = 1, numbasin
            nnball   (ibasin)   = icache1 (bindex(ibasin))
            idxnball (:,ibasin) = icache2 (:,bindex(ibasin))
            lenbdall (:,ibasin) = rcache2 (:,bindex(ibasin))
         ENDDO

         deallocate (icache1, icache2, rcache2)
#endif

         IF (numbasin > 0) THEN
           
            allocate (basinneighbour (numbasin))
            
            DO ibasin = 1, numbasin
               nnb = nnball(ibasin)
               basinneighbour(ibasin)%nnb = nnb

               IF (nnb > 0) THEN
                  allocate (basinneighbour(ibasin)%bindex (nnb))
                  allocate (basinneighbour(ibasin)%lenbdr (nnb))
                  allocate (basinneighbour(ibasin)%addr (2,nnb))
                  
                  basinneighbour(ibasin)%bindex = idxnball(1:nnb,ibasin)
                  basinneighbour(ibasin)%lenbdr = lenbdall(1:nnb,ibasin)
                  basinneighbour(ibasin)%addr(1,:) = -9999
               ENDIF
            ENDDO
         ENDIF

      ENDIF 

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN
               allocate (bindex  (nrecv))
               allocate (icache1 (nrecv))
               
               CALL mpi_recv (bindex, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               icache1 = addrbasin(bindex)

               idest = isrc
               CALL mpi_send (icache1, nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               deallocate(bindex, icache1)
            ENDIF
         ENDDO
      ENDIF 
#endif

      IF (p_is_worker) THEN

         IF (numbasin > 0) THEN
            allocate (basin_sorted (numbasin))
            allocate (order (numbasin))

            basin_sorted = bindex
            order = (/(ibasin, ibasin = 1, numbasin)/)

            CALL quicksort (numbasin, basin_sorted, order)

#ifdef USEMPI
            allocate(idxinq (numbasin*maxnnb)) 
#endif

            nnbinq = 0
            DO ibasin = 1, numbasin
               DO inb = 1, basinneighbour(ibasin)%nnb

                  IF (basinneighbour(ibasin)%bindex(inb) <= 0) CYCLE ! skip ocean neighbour

                  iloc = find_in_sorted_list1 (basinneighbour(ibasin)%bindex(inb), numbasin, basin_sorted)
                  IF (iloc > 0) THEN
                     basinneighbour(ibasin)%addr(1,inb) = -1
                     basinneighbour(ibasin)%addr(2,inb) = order(iloc)
#ifdef USEMPI
                  ELSE
                     CALL insert_into_sorted_list1 (basinneighbour(ibasin)%bindex(inb), nnbinq, idxinq, iloc)
#endif
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            nnbinq = 0
         ENDIF 

#ifdef USEMPI
         mesg(1:2) = (/p_iam_glb, nnbinq/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err) 

         IF (nnbinq > 0) THEN
            
            CALL mpi_send (idxinq(1:nnbinq), nnbinq, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_err) 

            allocate (addrinq (nnbinq))
            CALL mpi_recv (addrinq, nnbinq, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

         ENDIF

         IF (nnbinq > 0) allocate(mask (nnbinq))

         allocate (recvaddr (0:p_np_worker-1))
         DO iwork = 0, p_np_worker-1
            IF (nnbinq > 0) THEN
               mask = (addrinq == p_address_worker(iwork))
               ndata = count(mask)
            ELSE
               ndata = 0
            ENDIF 

            recvaddr(iwork)%ndata = ndata
            IF (ndata > 0) THEN
               allocate (recvaddr(iwork)%bindex (ndata))
               recvaddr(iwork)%bindex = pack(idxinq(1:nnbinq), mask)
            ENDIF
         ENDDO

         IF (nnbinq > 0) deallocate(mask)
      
         DO ibasin = 1, numbasin
            DO inb = 1, basinneighbour(ibasin)%nnb
               IF ((basinneighbour(ibasin)%addr(1,inb) == -9999) &
                  .and. (basinneighbour(ibasin)%bindex(inb) > 0)) THEN ! skip ocean neighbour

                  iloc = find_in_sorted_list1 (basinneighbour(ibasin)%bindex(inb), nnbinq, idxinq(1:nnbinq))

                  iwork = p_itis_worker(addrinq(iloc))
                  iloc1 = find_in_sorted_list1 (basinneighbour(ibasin)%bindex(inb), &
                     recvaddr(iwork)%ndata, recvaddr(iwork)%bindex)
                  
                  basinneighbour(ibasin)%addr(1,inb) = iwork
                  basinneighbour(ibasin)%addr(2,inb) = iloc1
               ENDIF
            ENDDO
         ENDDO

         allocate (sendaddr (0:p_np_worker-1))
         DO iwork = 0, p_np_worker-1
            sendaddr(iwork)%ndata = 0
         ENDDO

         DO ibasin = 1, numbasin
            DO inb = 1, basinneighbour(ibasin)%nnb
               IF (basinneighbour(ibasin)%addr(1,inb) >= 0) THEN
                  iwork = basinneighbour(ibasin)%addr(1,inb)
                  sendaddr(iwork)%ndata = sendaddr(iwork)%ndata + 1
               ENDIF
            ENDDO
         ENDDO

         DO iwork = 0, p_np_worker-1
            IF (sendaddr(iwork)%ndata > 0) THEN 
               allocate (sendaddr(iwork)%bindex (sendaddr(iwork)%ndata))
               sendaddr(iwork)%ndata = 0
            ENDIF
         ENDDO

         DO ibasin = 1, numbasin
            DO inb = 1, basinneighbour(ibasin)%nnb
               IF (basinneighbour(ibasin)%addr(1,inb) >= 0) THEN
                  iwork = basinneighbour(ibasin)%addr(1,inb)
                  CALL insert_into_sorted_list1 (bindex(ibasin), &
                     sendaddr(iwork)%ndata, sendaddr(iwork)%bindex, iloc)
               ENDIF
            ENDDO
         ENDDO
         
         DO iwork = 0, p_np_worker-1
            IF (sendaddr(iwork)%ndata > 0) THEN 
               IF (sendaddr(iwork)%ndata < size(sendaddr(iwork)%bindex)) THEN
                  allocate (icache1 (sendaddr(iwork)%ndata))
                  icache1 = sendaddr(iwork)%bindex(1:sendaddr(iwork)%ndata)

                  deallocate (sendaddr(iwork)%bindex)
                  allocate (sendaddr(iwork)%bindex (sendaddr(iwork)%ndata))
                  sendaddr(iwork)%bindex = icache1

                  deallocate (icache1)
               ENDIF
            ENDIF
         ENDDO
         
         DO iwork = 0, p_np_worker-1
            IF (sendaddr(iwork)%ndata > 0) THEN
               allocate (sendaddr(iwork)%ibasin (sendaddr(iwork)%ndata))

               DO inb = 1, sendaddr(iwork)%ndata
                  iloc = find_in_sorted_list1 (sendaddr(iwork)%bindex(inb), numbasin, basin_sorted)
                  sendaddr(iwork)%ibasin(inb) = order(iloc)
               ENDDO
            ENDIF
         ENDDO
#endif
      ENDIF

      IF (allocated(nnball   )) deallocate(nnball   )
      IF (allocated(idxnball )) deallocate(idxnball )
      IF (allocated(lenbdall )) deallocate(lenbdall )
      IF (allocated(addrbasin)) deallocate(addrbasin)
      IF (allocated(bindex ))   deallocate(bindex )
      IF (allocated(icache1))   deallocate(icache1)
      IF (allocated(icache2))   deallocate(icache2)
      IF (allocated(rcache2))   deallocate(rcache2)

      IF (allocated(basin_sorted)) deallocate(basin_sorted)
      IF (allocated(order  )) deallocate(order  )
      IF (allocated(idxinq )) deallocate(idxinq )
      IF (allocated(addrinq)) deallocate(addrinq)
      IF (allocated(mask   )) deallocate(mask   )

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_worker) THEN

         DO ibasin = 1, numbasin
            nnb = basinneighbour(ibasin)%nnb
            IF (nnb > 0) THEN
               allocate (basinneighbour(ibasin)%dist  (nnb))
               allocate (basinneighbour(ibasin)%area  (nnb))
               allocate (basinneighbour(ibasin)%elva  (nnb))
               allocate (basinneighbour(ibasin)%slope (nnb))
               allocate (basinneighbour(ibasin)%iswatb(nnb))
            ENDIF
         ENDDO

         IF (numbasin > 0) THEN
            allocate (rlon_b(numbasin))
            allocate (rlat_b(numbasin))
            CALL landelm%get_lonlat_radian (rlon_b, rlat_b)
         ENDIF

         CALL allocate_neighbour_data (rlon_nb)
         CALL allocate_neighbour_data (rlat_nb)

         CALL retrieve_neighbour_data (rlon_b, rlon_nb)
         CALL retrieve_neighbour_data (rlat_b, rlat_nb)
         
         DO ibasin = 1, numbasin
            DO inb = 1, basinneighbour(ibasin)%nnb
               IF (basinneighbour(ibasin)%bindex(inb) > 0) THEN ! skip ocean neighbour
                  basinneighbour(ibasin)%dist(inb) = 1.0e3 * arclen ( &
                     rlat_b (ibasin),          rlon_b (ibasin), &
                     rlat_nb(ibasin)%val(inb), rlon_nb(ibasin)%val(inb))
               ENDIF
            ENDDO
         ENDDO
         
         IF (numbasin > 0) THEN
            allocate (area_b(numbasin))
            allocate (elva_b(numbasin))
            allocate (iswatb(numbasin))
            DO ibasin = 1, numbasin
               IF (lake_id(ibasin) <= 0) THEN
                  area_b(ibasin) = sum(hillslope_network(ibasin)%area)
                  IF ((hillslope_network(ibasin)%nhru == 1) .and. (hillslope_network(ibasin)%indx(1) == 0)) THEN
                     iswatb(ibasin) = 1.
                  ELSE
                     iswatb(ibasin) = 0.
                  ENDIF
               ELSE
                  area_b(ibasin) = sum(lakes(ibasin)%area0)
                  iswatb(ibasin) = 1.
               ENDIF
               
               elva_b(ibasin) = basinelv(ibasin)
               IF (lake_id(ibasin) > 0) THEN
                  elva_b(ibasin) = bedelv(ibasin)
               ENDIF
               
               basinneighbour(ibasin)%myarea = area_b(ibasin)
               basinneighbour(ibasin)%myelva = elva_b(ibasin)
            ENDDO
         ENDIF
         
         CALL allocate_neighbour_data (area_nb)
         CALL retrieve_neighbour_data (area_b, area_nb)
         
         CALL allocate_neighbour_data (elva_nb)
         CALL retrieve_neighbour_data (elva_b, elva_nb)
         
         CALL allocate_neighbour_data (iswat_nb)
         CALL retrieve_neighbour_data (iswatb, iswat_nb)
         
         DO ibasin = 1, numbasin
            DO inb = 1, basinneighbour(ibasin)%nnb
               IF (basinneighbour(ibasin)%bindex(inb) > 0) THEN ! skip ocean neighbour
                  basinneighbour(ibasin)%area (inb) = area_nb(ibasin)%val(inb)
                  basinneighbour(ibasin)%elva (inb) = elva_nb(ibasin)%val(inb)
                  basinneighbour(ibasin)%slope(inb) = &
                     abs(elva_nb(ibasin)%val(inb) - elva_b(ibasin)) / basinneighbour(ibasin)%dist(inb)

                  IF (iswat_nb(ibasin)%val(inb) > 0) THEN
                     basinneighbour(ibasin)%iswatb(inb) = .true.
                  ELSE
                     basinneighbour(ibasin)%iswatb(inb) = .false.
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         
         CALL allocate_neighbour_data (theta_a_nb)
         CALL allocate_neighbour_data (zwt_nb    )
         CALL allocate_neighbour_data (Ks_nb     )
         CALL allocate_neighbour_data (wdsrf_nb  )

         IF (allocated(rlon_b)) deallocate(rlon_b)
         IF (allocated(rlat_b)) deallocate(rlat_b)
         IF (allocated(elva_b)) deallocate(elva_b)
         IF (allocated(area_b)) deallocate(area_b)
         IF (allocated(iswatb)) deallocate(iswatb)
         
         IF (allocated(rlon_nb )) deallocate(rlon_nb )
         IF (allocated(rlat_nb )) deallocate(rlat_nb )
         IF (allocated(area_nb )) deallocate(area_nb )
         IF (allocated(elva_nb )) deallocate(elva_nb )
         IF (allocated(iswat_nb)) deallocate(iswat_nb)
         
      ENDIF

      If (p_is_worker) THEN

         DO ibasin = 1, numbasin
            IF (lake_id(ibasin) == 0) THEN
               IF ((to_lake(ibasin)) .or. (riverdown(ibasin) <= 0)) THEN
                  ! river to lake, ocean or inland depression
                  outletwth(ibasin) = riverwth(ibasin)
               ELSE
                  ! river to river
                  outletwth(ibasin) = (riverwth(ibasin) + riverwth_ds(ibasin)) * 0.5
               ENDIF
            ELSEIF (lake_id(ibasin) /= 0) THEN
               IF ((.not. to_lake(ibasin)) .and. (riverdown(ibasin) /= 0)) THEN
                  IF (riverdown(ibasin) > 0) THEN
                     ! lake to river
                     outletwth(ibasin) = riverwth_ds(ibasin)
                  ELSEIF (riverdown(ibasin) == -1) THEN
                     ! lake is inland depression
                     outletwth(ibasin) = 0
                  ENDIF
               ELSEIF (to_lake(ibasin) .or. (riverdown(ibasin) == 0)) THEN
                  ! lake to lake .or. lake catchment to lake .or. lake to ocean
                  IF (riverdown(ibasin) > 0) THEN
                     inb = findloc(basinneighbour(ibasin)%bindex, riverdown(ibasin), dim=1)
                  ELSE
                     inb = findloc(basinneighbour(ibasin)%bindex, -9, dim=1) ! -9 is ocean
                  ENDIF

                  IF (inb <= 0) THEN
                     outletwth(ibasin) = 0
                  ELSE
                     outletwth(ibasin) = basinneighbour(ibasin)%lenbdr(inb)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDIF

   END SUBROUTINE basin_neighbour_init

   ! ----------
   SUBROUTINE retrieve_neighbour_data (vec_in, nbdata)

      USE MOD_Precision
      USE MOD_SPMD_Task
      USE MOD_Mesh, only : numelm
      IMPLICIT NONE

      REAL(r8), intent(inout) :: vec_in (:)
      TYPE(pointer_real8_1d)  :: nbdata (:)

      ! Local Variables
      LOGICAL, allocatable :: smask(:), rmask(:)
      INTEGER, allocatable :: req_send(:), req_recv(:)
      TYPE(pointer_real8_1d), allocatable :: sbuff(:), rbuff(:)
      INTEGER :: numbasin, iwork, ibasin, inb, iloc

      IF (p_is_worker) THEN
      
         numbasin = numelm

         DO ibasin = 1, numbasin
            DO inb = 1, basinneighbour(ibasin)%nnb
               IF (basinneighbour(ibasin)%addr(1,inb)== -1) THEN
                  iloc = basinneighbour(ibasin)%addr(2,inb)
                  nbdata(ibasin)%val(inb) = vec_in(iloc)
               ENDIF
            ENDDO
         ENDDO

#ifdef USEMPI
         CALL mpi_barrier (p_comm_worker, p_err)

         allocate (smask    (0:p_np_worker-1))
         allocate (req_send (0:p_np_worker-1))
         allocate (sbuff    (0:p_np_worker-1))
         smask(:) = .false.

         DO iwork = 0, p_np_worker-1
            IF (sendaddr(iwork)%ndata > 0) THEN
               smask(iwork) = .true.

               allocate (sbuff(iwork)%val (sendaddr(iwork)%ndata))
               sbuff(iwork)%val = vec_in(sendaddr(iwork)%ibasin)

               CALL mpi_isend(sbuff(iwork)%val, sendaddr(iwork)%ndata, MPI_REAL8, &
                  p_address_worker(iwork), 101, p_comm_glb, req_send(iwork), p_err)
            ENDIF
         ENDDO
         
         allocate (rmask    (0:p_np_worker-1))
         allocate (req_recv (0:p_np_worker-1))
         allocate (rbuff    (0:p_np_worker-1))
         rmask(:) = .false.

         DO iwork = 0, p_np_worker-1
            IF (recvaddr(iwork)%ndata > 0) THEN
               rmask(iwork) = .true.

               allocate (rbuff(iwork)%val (recvaddr(iwork)%ndata))

               CALL mpi_irecv(rbuff(iwork)%val, recvaddr(iwork)%ndata, MPI_REAL8, &
                  p_address_worker(iwork), 101, p_comm_glb, req_recv(iwork), p_err)
            ENDIF
         ENDDO

         IF (any(rmask)) THEN

            CALL mpi_waitall(count(rmask), pack(req_recv,rmask), MPI_STATUSES_IGNORE, p_err)

            DO ibasin = 1, numbasin
               DO inb = 1, basinneighbour(ibasin)%nnb
                  IF (basinneighbour(ibasin)%addr(1,inb) >= 0) THEN
                     iwork = basinneighbour(ibasin)%addr(1,inb)
                     iloc  = basinneighbour(ibasin)%addr(2,inb)
                     nbdata(ibasin)%val(inb) = rbuff(iwork)%val(iloc)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         IF (any(smask)) THEN
            CALL mpi_waitall(count(smask), pack(req_send,smask), MPI_STATUSES_IGNORE, p_err)
         ENDIF

         IF (allocated(smask)) deallocate(smask)
         IF (allocated(rmask)) deallocate(rmask)
         
         IF (allocated(req_send)) deallocate(req_send)
         IF (allocated(req_recv)) deallocate(req_recv)
         
         IF (allocated(sbuff)) deallocate(sbuff)
         IF (allocated(rbuff)) deallocate(rbuff)

         CALL mpi_barrier (p_comm_worker, p_err)
#endif

      ENDIF

   END SUBROUTINE retrieve_neighbour_data

   ! ---
   SUBROUTINE allocate_neighbour_data (nbdata)
      
      USE MOD_Mesh, only : numelm
      IMPLICIT NONE

      TYPE(pointer_real8_1d), allocatable :: nbdata(:)
      INTEGER :: ibasin
            
      IF (numelm > 0) THEN
         allocate (nbdata(numelm))
         DO ibasin = 1, numelm
            IF (basinneighbour(ibasin)%nnb > 0) THEN
               allocate (nbdata(ibasin)%val (basinneighbour(ibasin)%nnb))
            ENDIF
         ENDDO
      ENDIF 

   END SUBROUTINE allocate_neighbour_data 

   ! ----------
   SUBROUTINE basin_neighbour_final ()

      IMPLICIT NONE
      INTEGER :: i

      IF (allocated(basinneighbour)) THEN
         DO i = 1, size(basinneighbour)
            IF (allocated(basinneighbour(i)%bindex)) deallocate(basinneighbour(i)%bindex)
            IF (allocated(basinneighbour(i)%addr  )) deallocate(basinneighbour(i)%addr  )
            IF (allocated(basinneighbour(i)%dist  )) deallocate(basinneighbour(i)%dist  )
            IF (allocated(basinneighbour(i)%lenbdr)) deallocate(basinneighbour(i)%lenbdr)
            IF (allocated(basinneighbour(i)%area  )) deallocate(basinneighbour(i)%area  )
            IF (allocated(basinneighbour(i)%elva  )) deallocate(basinneighbour(i)%elva  )
            IF (allocated(basinneighbour(i)%slope )) deallocate(basinneighbour(i)%slope )
            IF (allocated(basinneighbour(i)%iswatb)) deallocate(basinneighbour(i)%iswatb)
         ENDDO
         deallocate(basinneighbour)
      ENDIF
      
      IF (allocated(theta_a_nb)) deallocate(theta_a_nb)
      IF (allocated(zwt_nb    )) deallocate(zwt_nb    )
      IF (allocated(Ks_nb     )) deallocate(Ks_nb     )
      IF (allocated(wdsrf_nb  )) deallocate(wdsrf_nb  )
      
      DO i = lbound(recvaddr,1), ubound(recvaddr,1)
         IF (allocated(recvaddr(i)%bindex)) deallocate(recvaddr(i)%bindex)
         IF (allocated(recvaddr(i)%ibasin )) deallocate(recvaddr(i)%ibasin )
      ENDDO 

      DO i = lbound(sendaddr,1), ubound(sendaddr,1)
         IF (allocated(sendaddr(i)%bindex)) deallocate(sendaddr(i)%bindex)
         IF (allocated(sendaddr(i)%ibasin )) deallocate(sendaddr(i)%ibasin )
      ENDDO
         
      IF (allocated(recvaddr)) deallocate(recvaddr)
      IF (allocated(sendaddr)) deallocate(sendaddr)

   END SUBROUTINE basin_neighbour_final

END MODULE MOD_Hydro_BasinNeighbour
#endif
