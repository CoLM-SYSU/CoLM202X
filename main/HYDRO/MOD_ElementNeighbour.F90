#include <define.h>

MODULE MOD_ElementNeighbour
!--------------------------------------------------------------------------------!
! DESCRIPTION:                                                                   !
!                                                                                !
!    Element Neighbours : data and communication subroutines.                    !
!                                                                                !
! Created by Shupeng Zhang, May 2023                                             !
!--------------------------------------------------------------------------------!

   USE MOD_Precision
   USE MOD_DataType
   IMPLICIT NONE

   ! -- neighbour parameters --
   type element_neighbour_type

      integer  :: nnb    ! number of neighbours
      real(r8) :: myarea ! area of this element [m^2]
      real(r8) :: myelva ! elevation of this element [m]

      integer*8, allocatable :: glbindex (:) ! neighbour global index

      ! data address: (1,:) refers to process, (2,:) refers to location
      integer , allocatable :: addr (:,:)

      real(r8), allocatable :: dist   (:) ! distance between element centers [m]
      real(r8), allocatable :: lenbdr (:) ! length of boundary line [m]
      real(r8), allocatable :: area   (:) ! area of neighbours [m^2]
      real(r8), allocatable :: elva   (:) ! elevation of neighbours [m]
      real(r8), allocatable :: slope  (:) ! slope (positive) [-]

   END type element_neighbour_type

   type(element_neighbour_type), allocatable :: elementneighbour (:)

   ! -- neighbour communication --
   type neighbour_sendrecv_type
      integer :: ndata
      integer*8, allocatable :: glbindex (:)
      integer,   allocatable :: ielement (:)
   END type neighbour_sendrecv_type

   type(neighbour_sendrecv_type), allocatable :: recvaddr(:)
   type(neighbour_sendrecv_type), allocatable :: sendaddr(:)

   INTERFACE allocate_neighbour_data
      MODULE procedure allocate_neighbour_data_real8
      MODULE procedure allocate_neighbour_data_logic
   END INTERFACE allocate_neighbour_data

CONTAINS

   ! ----------
   SUBROUTINE element_neighbour_init (patcharea, lc_year)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_NetCDFVector
   USE MOD_Mesh
   USE MOD_Pixel
   USE MOD_LandElm
   USE MOD_LandPatch
   USE MOD_Utils
   IMPLICIT NONE

   integer,  intent(in) :: lc_year    ! which year of land cover data used
   real(r8), intent(in) :: patcharea (:)

   ! Local Variables
   character(len=256) :: neighbour_file

   integer :: ielm
   integer :: iwork, mesg(2), isrc, idest
   integer :: nrecv, irecv
   integer :: iloc, iloc1, iloc2
   integer :: nnb, nnbinq, inb, ndata

   integer :: maxnnb
   integer , allocatable :: nnball   (:)
   integer , allocatable :: idxnball (:,:)
   real(r8), allocatable :: lenbdall (:,:)

   integer , allocatable :: addrelement(:)

   integer*8, allocatable :: eindex  (:)
   integer,   allocatable :: icache1 (:)
   integer,   allocatable :: icache2 (:,:)
   real(r8),  allocatable :: rcache2 (:,:)

   integer*8, allocatable :: elm_sorted (:)
   integer,   allocatable :: order      (:)
   integer*8, allocatable :: idxinq     (:)
   integer,   allocatable :: addrinq    (:)

   logical, allocatable :: mask(:)

   real(r8), allocatable :: rlon_b(:), rlat_b(:)
   type(pointer_real8_1d), allocatable :: rlon_nb(:), rlat_nb(:)

   real(r8), allocatable :: area_b(:)
   real(r8), allocatable :: elva_b(:)

   character(len=256) :: lndname, cyear
   real(r8), allocatable :: elv_patches(:)

   type(pointer_real8_1d), allocatable :: area_nb (:)  ! m^2
   type(pointer_real8_1d), allocatable :: elva_nb (:)  ! m

   integer :: istt, iend

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      neighbour_file = DEF_ElementNeighbour_file

      IF (p_is_master) THEN
         CALL ncio_read_serial (neighbour_file, 'num_neighbour', nnball  )
         CALL ncio_read_serial (neighbour_file, 'idx_neighbour', idxnball)
         CALL ncio_read_serial (neighbour_file, 'len_border'   , lenbdall)

         maxnnb = size(idxnball,1)

         lenbdall = lenbdall * 1.e3 ! km to m
      ENDIF

#ifdef USEMPI

      CALL mpi_bcast (maxnnb, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)

      IF (p_is_master) THEN

         allocate (addrelement (size(nnball)))
         addrelement(:) = -1

         DO iwork = 0, p_np_worker-1

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN

               allocate (eindex  (nrecv))
               allocate (icache1 (nrecv))
               allocate (icache2 (maxnnb,nrecv))
               allocate (rcache2 (maxnnb,nrecv))

               CALL mpi_recv (eindex, nrecv, MPI_INTEGER8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               addrelement(eindex) = isrc

               idest = isrc

               icache1 = nnball(eindex)
               CALL mpi_send (icache1, nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               DO irecv = 1, nrecv
                  icache2(:,irecv) = idxnball(:,eindex(irecv))
               ENDDO
               CALL mpi_send (icache2, maxnnb*nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               DO irecv = 1, nrecv
                  rcache2(:,irecv) = lenbdall(:,eindex(irecv))
               ENDDO
               CALL mpi_send (rcache2, maxnnb*nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (eindex )
               deallocate (icache1)
               deallocate (icache2)
               deallocate (rcache2)

            ENDIF
         ENDDO
      ENDIF
#endif

      IF (p_is_worker) THEN

         IF (numelm > 0) THEN
            allocate (eindex (numelm))
            eindex = landelm%eindex
         ENDIF

#ifdef USEMPI
         mesg(1:2) = (/p_iam_glb, numelm/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (numelm > 0) THEN
            CALL mpi_send (eindex, numelm, MPI_INTEGER8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)

            allocate (nnball (numelm))
            CALL mpi_recv (nnball, numelm, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (idxnball (maxnnb,numelm))
            CALL mpi_recv (idxnball, maxnnb*numelm, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (lenbdall (maxnnb,numelm))
            CALL mpi_recv (lenbdall, maxnnb*numelm, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
#else
         allocate (icache1 (numelm))
         allocate (icache2 (maxnnb,numelm))
         allocate (rcache2 (maxnnb,numelm))

         icache1 = nnball
         icache2 = idxnball
         rcache2 = lenbdall

         DO ielm = 1, numelm
            nnball   (ielm)   = icache1 (eindex(ielm))
            idxnball (:,ielm) = icache2 (:,eindex(ielm))
            lenbdall (:,ielm) = rcache2 (:,eindex(ielm))
         ENDDO

         deallocate (icache1, icache2, rcache2)
#endif

         IF (numelm > 0) THEN

            allocate (elementneighbour (numelm))

            DO ielm = 1, numelm
               nnb = nnball(ielm)
               elementneighbour(ielm)%nnb = nnb

               IF (nnb > 0) THEN
                  allocate (elementneighbour(ielm)%glbindex (nnb))
                  allocate (elementneighbour(ielm)%lenbdr (nnb))
                  allocate (elementneighbour(ielm)%addr (2,nnb))

                  elementneighbour(ielm)%glbindex = idxnball(1:nnb,ielm)
                  elementneighbour(ielm)%lenbdr = lenbdall(1:nnb,ielm)
                  elementneighbour(ielm)%addr(1,:) = -9999
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
               allocate (eindex  (nrecv))
               allocate (icache1 (nrecv))

               CALL mpi_recv (eindex, nrecv, MPI_INTEGER8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               icache1 = addrelement(eindex)

               idest = isrc
               CALL mpi_send (icache1, nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate(eindex, icache1)
            ENDIF
         ENDDO
      ENDIF
#endif

      IF (p_is_worker) THEN

         IF (numelm > 0) THEN
            allocate (elm_sorted (numelm))
            allocate (order (numelm))

            elm_sorted = eindex
            order = (/(ielm, ielm = 1, numelm)/)

            CALL quicksort (numelm, elm_sorted, order)

#ifdef USEMPI
            allocate(idxinq (numelm*maxnnb))
#endif

            nnbinq = 0
            DO ielm = 1, numelm
               DO inb = 1, elementneighbour(ielm)%nnb

                  IF (elementneighbour(ielm)%glbindex(inb) <= 0) CYCLE ! skip ocean neighbour

                  iloc = find_in_sorted_list1 (elementneighbour(ielm)%glbindex(inb), numelm, elm_sorted)
                  IF (iloc > 0) THEN
                     elementneighbour(ielm)%addr(1,inb) = -1
                     elementneighbour(ielm)%addr(2,inb) = order(iloc)
#ifdef USEMPI
                  ELSE
                     CALL insert_into_sorted_list1 (elementneighbour(ielm)%glbindex(inb), nnbinq, idxinq, iloc)
#endif
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            nnbinq = 0
         ENDIF

#ifdef USEMPI
         mesg(1:2) = (/p_iam_glb, nnbinq/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (nnbinq > 0) THEN

            CALL mpi_send (idxinq(1:nnbinq), nnbinq, MPI_INTEGER8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)

            allocate (addrinq (nnbinq))
            CALL mpi_recv (addrinq, nnbinq, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

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
               allocate (recvaddr(iwork)%glbindex (ndata))
               recvaddr(iwork)%glbindex = pack(idxinq(1:nnbinq), mask)
            ENDIF
         ENDDO

         IF (nnbinq > 0) deallocate(mask)

         DO ielm = 1, numelm
            DO inb = 1, elementneighbour(ielm)%nnb
               IF ((elementneighbour(ielm)%addr(1,inb) == -9999) &
                  .and. (elementneighbour(ielm)%glbindex(inb) > 0)) THEN ! skip ocean neighbour

                  iloc = find_in_sorted_list1 (elementneighbour(ielm)%glbindex(inb), &
                     nnbinq, idxinq(1:nnbinq))

                  iwork = p_itis_worker(addrinq(iloc))
                  iloc1 = find_in_sorted_list1 (elementneighbour(ielm)%glbindex(inb), &
                     recvaddr(iwork)%ndata, recvaddr(iwork)%glbindex)

                  elementneighbour(ielm)%addr(1,inb) = iwork
                  elementneighbour(ielm)%addr(2,inb) = iloc1
               ENDIF
            ENDDO
         ENDDO

         allocate (sendaddr (0:p_np_worker-1))
         DO iwork = 0, p_np_worker-1
            sendaddr(iwork)%ndata = 0
         ENDDO

         DO ielm = 1, numelm
            DO inb = 1, elementneighbour(ielm)%nnb
               IF (elementneighbour(ielm)%addr(1,inb) >= 0) THEN
                  iwork = elementneighbour(ielm)%addr(1,inb)
                  sendaddr(iwork)%ndata = sendaddr(iwork)%ndata + 1
               ENDIF
            ENDDO
         ENDDO

         DO iwork = 0, p_np_worker-1
            IF (sendaddr(iwork)%ndata > 0) THEN
               allocate (sendaddr(iwork)%glbindex (sendaddr(iwork)%ndata))
               sendaddr(iwork)%ndata = 0
            ENDIF
         ENDDO

         DO ielm = 1, numelm
            DO inb = 1, elementneighbour(ielm)%nnb
               IF (elementneighbour(ielm)%addr(1,inb) >= 0) THEN
                  iwork = elementneighbour(ielm)%addr(1,inb)
                  CALL insert_into_sorted_list1 (eindex(ielm), &
                     sendaddr(iwork)%ndata, sendaddr(iwork)%glbindex, iloc)
               ENDIF
            ENDDO
         ENDDO

         DO iwork = 0, p_np_worker-1
            IF (sendaddr(iwork)%ndata > 0) THEN
               IF (sendaddr(iwork)%ndata < size(sendaddr(iwork)%glbindex)) THEN
                  allocate (icache1 (sendaddr(iwork)%ndata))
                  icache1 = sendaddr(iwork)%glbindex(1:sendaddr(iwork)%ndata)

                  deallocate (sendaddr(iwork)%glbindex)
                  allocate (sendaddr(iwork)%glbindex (sendaddr(iwork)%ndata))
                  sendaddr(iwork)%glbindex = icache1

                  deallocate (icache1)
               ENDIF
            ENDIF
         ENDDO

         DO iwork = 0, p_np_worker-1
            IF (sendaddr(iwork)%ndata > 0) THEN
               allocate (sendaddr(iwork)%ielement (sendaddr(iwork)%ndata))

               DO inb = 1, sendaddr(iwork)%ndata
                  iloc = find_in_sorted_list1 (sendaddr(iwork)%glbindex(inb), numelm, elm_sorted)
                  sendaddr(iwork)%ielement(inb) = order(iloc)
               ENDDO
            ENDIF
         ENDDO
#endif
      ENDIF

      IF (allocated(addrelement))   deallocate(addrelement)
      IF (allocated(elm_sorted ))   deallocate(elm_sorted )
      IF (allocated(nnball     ))   deallocate(nnball     )
      IF (allocated(idxnball   ))   deallocate(idxnball   )
      IF (allocated(lenbdall   ))   deallocate(lenbdall   )
      IF (allocated(eindex     ))   deallocate(eindex     )
      IF (allocated(icache1    ))   deallocate(icache1    )
      IF (allocated(icache2    ))   deallocate(icache2    )
      IF (allocated(rcache2    ))   deallocate(rcache2    )
      IF (allocated(order      ))   deallocate(order      )
      IF (allocated(idxinq     ))   deallocate(idxinq     )
      IF (allocated(addrinq    ))   deallocate(addrinq    )
      IF (allocated(mask       ))   deallocate(mask       )

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      write(cyear,'(i4.4)') lc_year
      lndname = trim(DEF_dir_landdata) // '/topography/'//trim(cyear)//'/elevation_patches.nc'
      CALL ncio_read_vector (lndname, 'elevation_patches', landpatch, elv_patches)

      IF (p_is_worker) THEN

         DO ielm = 1, numelm
            nnb = elementneighbour(ielm)%nnb
            IF (nnb > 0) THEN
               allocate (elementneighbour(ielm)%dist  (nnb))
               allocate (elementneighbour(ielm)%area  (nnb))
               allocate (elementneighbour(ielm)%elva  (nnb))
               allocate (elementneighbour(ielm)%slope (nnb))
            ENDIF
         ENDDO

         IF (numelm > 0) THEN
            allocate (rlon_b(numelm))
            allocate (rlat_b(numelm))
            CALL landelm%get_lonlat_radian (rlon_b, rlat_b)
         ENDIF

         CALL allocate_neighbour_data (rlon_nb)
         CALL allocate_neighbour_data (rlat_nb)

         CALL retrieve_neighbour_data (rlon_b, rlon_nb)
         CALL retrieve_neighbour_data (rlat_b, rlat_nb)

         DO ielm = 1, numelm
            DO inb = 1, elementneighbour(ielm)%nnb
               IF (elementneighbour(ielm)%glbindex(inb) > 0) THEN ! skip ocean neighbour
                  elementneighbour(ielm)%dist(inb) = 1.0e3 * arclen ( &
                     rlat_b (ielm),          rlon_b (ielm), &
                     rlat_nb(ielm)%val(inb), rlon_nb(ielm)%val(inb))
                  elementneighbour(ielm)%dist(inb) = max(elementneighbour(ielm)%dist(inb), 90.)
               ENDIF
            ENDDO
         ENDDO

         IF (numelm > 0) THEN
            allocate (area_b(numelm))
            allocate (elva_b(numelm))
            DO ielm = 1, numelm
               istt = elm_patch%substt(ielm)
               iend = elm_patch%subend(ielm)
               area_b(ielm) = sum(patcharea(istt:iend))
               elva_b(ielm) = sum(elv_patches(istt:iend) * elm_patch%subfrc(istt:iend))

               elementneighbour(ielm)%myarea = area_b(ielm)
               elementneighbour(ielm)%myelva = elva_b(ielm)
            ENDDO
         ENDIF

         CALL allocate_neighbour_data (area_nb)
         CALL retrieve_neighbour_data (area_b, area_nb)

         CALL allocate_neighbour_data (elva_nb)
         CALL retrieve_neighbour_data (elva_b, elva_nb)

         DO ielm = 1, numelm
            DO inb = 1, elementneighbour(ielm)%nnb
               IF (elementneighbour(ielm)%glbindex(inb) > 0) THEN ! skip ocean neighbour
                  elementneighbour(ielm)%area (inb) = area_nb(ielm)%val(inb)
                  elementneighbour(ielm)%elva (inb) = elva_nb(ielm)%val(inb)
                  elementneighbour(ielm)%slope(inb) = &
                     abs(elva_nb(ielm)%val(inb) - elva_b(ielm)) / elementneighbour(ielm)%dist(inb)
               ENDIF
            ENDDO
         ENDDO

         IF (allocated(rlon_b )) deallocate(rlon_b )
         IF (allocated(rlat_b )) deallocate(rlat_b )
         IF (allocated(elva_b )) deallocate(elva_b )
         IF (allocated(area_b )) deallocate(area_b )
         IF (allocated(rlon_nb)) deallocate(rlon_nb)
         IF (allocated(rlat_nb)) deallocate(rlat_nb)
         IF (allocated(area_nb)) deallocate(area_nb)
         IF (allocated(elva_nb)) deallocate(elva_nb)

      ENDIF

   END SUBROUTINE element_neighbour_init

   ! ----------
   SUBROUTINE retrieve_neighbour_data (vec_in, nbdata)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Mesh, only: numelm
   IMPLICIT NONE

   real(r8), intent(inout) :: vec_in (:)
   type(pointer_real8_1d)  :: nbdata (:)

   ! Local Variables
   logical, allocatable :: smask(:), rmask(:)
   integer, allocatable :: req_send(:), req_recv(:)
   type(pointer_real8_1d), allocatable :: sbuff(:), rbuff(:)
   integer :: iwork, ielm, inb, iloc

      IF (p_is_worker) THEN

         DO ielm = 1, numelm
            DO inb = 1, elementneighbour(ielm)%nnb
               IF (elementneighbour(ielm)%addr(1,inb)== -1) THEN
                  iloc = elementneighbour(ielm)%addr(2,inb)
                  nbdata(ielm)%val(inb) = vec_in(iloc)
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
               sbuff(iwork)%val = vec_in(sendaddr(iwork)%ielement)

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

            DO ielm = 1, numelm
               DO inb = 1, elementneighbour(ielm)%nnb
                  IF (elementneighbour(ielm)%addr(1,inb) >= 0) THEN
                     iwork = elementneighbour(ielm)%addr(1,inb)
                     iloc  = elementneighbour(ielm)%addr(2,inb)
                     nbdata(ielm)%val(inb) = rbuff(iwork)%val(iloc)
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
   SUBROUTINE allocate_neighbour_data_real8 (nbdata)

   USE MOD_Mesh, only: numelm
   IMPLICIT NONE

   type(pointer_real8_1d), allocatable :: nbdata(:)
   integer :: ielm

      IF (numelm > 0) THEN
         allocate (nbdata(numelm))
         DO ielm = 1, numelm
            IF (elementneighbour(ielm)%nnb > 0) THEN
               allocate (nbdata(ielm)%val (elementneighbour(ielm)%nnb))
            ENDIF
         ENDDO
      ENDIF

   END SUBROUTINE allocate_neighbour_data_real8

   ! ---
   SUBROUTINE allocate_neighbour_data_logic (nbdata)

   USE MOD_Mesh, only: numelm
   IMPLICIT NONE

   type(pointer_logic_1d), allocatable :: nbdata(:)
   integer :: ielm

      IF (numelm > 0) THEN
         allocate (nbdata(numelm))
         DO ielm = 1, numelm
            IF (elementneighbour(ielm)%nnb > 0) THEN
               allocate (nbdata(ielm)%val (elementneighbour(ielm)%nnb))
            ENDIF
         ENDDO
      ENDIF

   END SUBROUTINE allocate_neighbour_data_logic

   ! ----------
   SUBROUTINE element_neighbour_final ()

   IMPLICIT NONE
   integer :: i

      IF (allocated(elementneighbour)) THEN
         DO i = 1, size(elementneighbour)
            IF (allocated(elementneighbour(i)%glbindex)) deallocate(elementneighbour(i)%glbindex)
            IF (allocated(elementneighbour(i)%addr  )) deallocate(elementneighbour(i)%addr  )
            IF (allocated(elementneighbour(i)%dist  )) deallocate(elementneighbour(i)%dist  )
            IF (allocated(elementneighbour(i)%lenbdr)) deallocate(elementneighbour(i)%lenbdr)
            IF (allocated(elementneighbour(i)%area  )) deallocate(elementneighbour(i)%area  )
            IF (allocated(elementneighbour(i)%elva  )) deallocate(elementneighbour(i)%elva  )
            IF (allocated(elementneighbour(i)%slope )) deallocate(elementneighbour(i)%slope )
         ENDDO
         deallocate(elementneighbour)
      ENDIF

      IF (allocated(recvaddr)) THEN
         DO i = lbound(recvaddr,1), ubound(recvaddr,1)
            IF (allocated(recvaddr(i)%glbindex)) deallocate(recvaddr(i)%glbindex)
            IF (allocated(recvaddr(i)%ielement)) deallocate(recvaddr(i)%ielement)
         ENDDO
      ENDIF

      IF (allocated(sendaddr)) THEN
         DO i = lbound(sendaddr,1), ubound(sendaddr,1)
            IF (allocated(sendaddr(i)%glbindex)) deallocate(sendaddr(i)%glbindex)
            IF (allocated(sendaddr(i)%ielement)) deallocate(sendaddr(i)%ielement)
         ENDDO
      ENDIF

      IF (allocated(recvaddr)) deallocate(recvaddr)
      IF (allocated(sendaddr)) deallocate(sendaddr)

   END SUBROUTINE element_neighbour_final

END MODULE MOD_ElementNeighbour
