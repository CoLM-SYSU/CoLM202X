#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_RiverLakeNetwork
!--------------------------------------------------------------------------------
! DESCRIPTION:
! 
!    River networks: data and communication subroutines.
!
! Created by Shupeng Zhang, May 2023
!--------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global, only : spval
   IMPLICIT NONE

   ! -- river parameters --
   real(r8), allocatable :: riverlen  (:)
   real(r8), allocatable :: riverelv  (:)
   real(r8), allocatable :: riverwth  (:)
   real(r8), allocatable :: riverarea (:)
   real(r8), allocatable :: riverdpth (:)

   real(r8), allocatable :: basinelv  (:)
   real(r8), allocatable :: bedelv    (:)
   real(r8), allocatable :: handmin   (:)
   
   real(r8), allocatable :: wtsrfelv  (:)

   ! index of downstream river 
   ! > 0 : other catchment;   0 : river mouth; -1 : inland depression
   integer, allocatable :: riverdown  (:) 
   logical, allocatable :: to_lake (:)

   ! address of downstream river 
   ! > 0 : catchment on this process;   0 : catchment on other processes; 
   ! -1  : not found, including river mouth, out of domain, inland depression.
   integer, allocatable :: addrdown (:)

   real(r8), allocatable :: riverlen_ds  (:)
   real(r8), allocatable :: wtsrfelv_ds  (:)
   real(r8), allocatable :: riverwth_ds  (:)
   real(r8), allocatable :: bedelv_ds    (:)
   
   real(r8), allocatable :: outletwth (:)

   ! -- lake data type --
   type :: lake_info_type
      integer :: nsub
      real(r8), allocatable :: area0  (:) ! area data in HRU order
      real(r8), allocatable :: area   (:) ! area data in the order from deepest to shallowest HRU
      real(r8), allocatable :: depth0 (:) ! depth data in HRU order
      real(r8), allocatable :: depth  (:) ! depth data in the order from deepest to shallowest HRU
      ! a curve describing the relationship between depth of water from lake bottom and total water volume
      ! the i-th value corresponds to the volume when water depth is at i-th depth
      real(r8), allocatable :: dep_vol_curve (:)
   CONTAINS 
      procedure, PUBLIC :: surface => retrieve_lake_surface_from_volume
      procedure, PUBLIC :: volume  => retrieve_lake_volume_from_surface
   END type lake_info_type

   ! -- lake information --
   integer, allocatable :: lake_id (:)
   type(lake_info_type), allocatable :: lakes (:)

   ! -- communications --
   type :: river_sendrecv_type
      integer :: nproc
      integer, allocatable :: iproc (:)
      integer, allocatable :: wdsp  (:)
      integer, allocatable :: ndata (:)
      integer, allocatable :: ups   (:)
      integer, allocatable :: down  (:)
      integer, allocatable :: iloc  (:)
   CONTAINS
      final :: river_sendrecv_free_mem
   END type river_sendrecv_type

   type(river_sendrecv_type), target :: river_up
   type(river_sendrecv_type), target :: river_dn

   integer, parameter :: SEND_DATA_DOWN_TO_UP = 1
   integer, parameter :: SEND_DATA_UP_TO_DOWN = 2

CONTAINS
   
   ! ----------
   SUBROUTINE river_lake_network_init ()

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   USE MOD_Pixel
   USE MOD_LandElm
   USE MOD_LandPatch
   USE MOD_Catch_HillslopeNetwork
   USE MOD_ElementNeighbour
   USE MOD_DataType
   USE MOD_Utils
   USE MOD_UserDefFun
   USE MOD_Vars_TimeInvariants, only : lakedepth
   IMPLICIT NONE

   ! Local Variables
   character(len=256) :: river_file, rivdpt_file
   logical :: use_calc_rivdpt

   integer :: numbasin, ibasin, nbasin, inb
   integer :: iworker, mesg(4), isrc, idest, iproc
   integer :: irecv, ifrom, ito, iup, idn, idata
   integer :: nrecv, ndata, nup, ndn
   integer :: iloc, iloc1, iloc2
   
   integer , allocatable :: bindex (:)
   integer , allocatable :: icache (:)
   real(r8), allocatable :: rcache (:)
   logical , allocatable :: lcache (:)
   
   integer , allocatable :: addrbasin (:,:)
   integer , allocatable :: ndata_w (:)

   type(pointer_int32_2d), allocatable :: exchange_w (:)
   integer, allocatable :: exchange(:,:)
   integer, allocatable :: basin_sorted(:), order(:)

   ! for lakes
   integer :: ps, pe, nsublake, i, ipatch, ipxl

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      numbasin = numelm

      use_calc_rivdpt = DEF_USE_EstimatedRiverDepth
      river_file      = DEF_CatchmentMesh_data 

      IF (p_is_master) THEN
         
         CALL ncio_read_serial (river_file, 'lake_id'         , lake_id  )
         CALL ncio_read_serial (river_file, 'basin_downstream', riverdown)
         CALL ncio_read_serial (river_file, 'river_length'    , riverlen )
         CALL ncio_read_serial (river_file, 'river_elevation' , riverelv )
         CALL ncio_read_serial (river_file, 'basin_elevation' , basinelv )
         
         IF (.not. use_calc_rivdpt) THEN
            CALL ncio_read_serial (river_file, 'river_depth' , riverdpth)
         ENDIF

         riverlen = riverlen * 1.e3 ! km to m
         
         nbasin = size(riverdown)
         allocate (to_lake (nbasin))
         to_lake = .false.
         DO i = 1, nbasin
            IF (riverdown(i) > 0) THEN
               to_lake(i) = lake_id(riverdown(i)) > 0 
            ENDIF
         ENDDO

      ENDIF
         
      IF (use_calc_rivdpt) THEN
         ! Estimate river depth by using runoff data.
         CALL calc_riverdepth_from_runoff ()
      ENDIF

#ifdef USEMPI
      IF (p_is_master) THEN

         nbasin = size(riverdown)

         allocate (addrbasin (2,nbasin))
         addrbasin(:,:) = -1

         DO iworker = 1, p_np_worker

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN
               
               allocate (bindex (nrecv))
               allocate (icache (nrecv))
               allocate (rcache (nrecv))
               allocate (lcache (nrecv))

               CALL mpi_recv (bindex, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO irecv = 1, nrecv
                  addrbasin(1,bindex(irecv)) = isrc
               ENDDO
               
               idest = isrc

               DO irecv = 1, nrecv
                  icache(irecv) = lake_id(bindex(irecv))
               ENDDO
               CALL mpi_send (icache, nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  icache(irecv) = riverdown(bindex(irecv))
               ENDDO
               CALL mpi_send (icache, nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  lcache(irecv) = to_lake(bindex(irecv))
               ENDDO
               CALL mpi_send (lcache, nrecv, MPI_LOGICAL, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(irecv) = riverlen(bindex(irecv))
               ENDDO
               CALL mpi_send (rcache, nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(irecv) = riverelv(bindex(irecv))
               ENDDO
               CALL mpi_send (rcache, nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(irecv) = riverdpth(bindex(irecv))
               ENDDO
               CALL mpi_send (rcache, nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(irecv) = basinelv(bindex(irecv))
               ENDDO
               CALL mpi_send (rcache, nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               deallocate (bindex)
               deallocate (icache)
               deallocate (rcache)
               deallocate (lcache)

            ENDIF
         ENDDO

      ENDIF
#endif

      IF (p_is_worker) THEN
               
         IF (numbasin > 0) THEN
            allocate (bindex (numbasin))
            DO ibasin = 1, numbasin
               bindex(ibasin) = mesh(ibasin)%indx
            ENDDO
         ENDIF 

#ifdef USEMPI
         mesg(1:2) = (/p_iam_glb, numbasin/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err) 

         IF (numbasin > 0) THEN
            CALL mpi_send (bindex, numbasin, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err) 
            
            allocate (lake_id (numbasin))
            CALL mpi_recv (lake_id, numbasin, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (riverdown (numbasin))
            CALL mpi_recv (riverdown, numbasin, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (to_lake (numbasin))
            CALL mpi_recv (to_lake, numbasin, MPI_LOGICAL, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (riverlen (numbasin))
            CALL mpi_recv (riverlen, numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (riverelv (numbasin))
            CALL mpi_recv (riverelv, numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (riverdpth (numbasin))
            CALL mpi_recv (riverdpth, numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (basinelv (numbasin))
            CALL mpi_recv (basinelv, numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
#else
         IF (numbasin > 0) THEN

            lake_id   = lake_id  (bindex)
            riverdown = riverdown(bindex)
            to_lake   = to_lake  (bindex)
            riverlen  = riverlen (bindex)
            riverelv  = riverelv (bindex)
            riverdpth = riverdpth(bindex)
            basinelv  = basinelv (bindex)

         ENDIF
#endif

      ENDIF 

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
      IF (p_is_master) THEN

         allocate (ndata_w (0:p_np_worker-1))

         ndata_w(:) = 0
         DO ibasin = 1, nbasin
            IF (riverdown(ibasin) >= 1) THEN
               addrbasin(2,ibasin) = addrbasin(1,riverdown(ibasin))
               IF ((addrbasin(1,ibasin) /= -1) .and. (addrbasin(2,ibasin) /= -1) &
                  .and. (addrbasin(1,ibasin) /= addrbasin(2,ibasin))) THEN
                  ifrom = p_itis_worker(addrbasin(1,ibasin))
                  ito   = p_itis_worker(addrbasin(2,ibasin))
                  ndata_w(ifrom) = ndata_w(ifrom) + 1
                  ndata_w(ito)   = ndata_w(ito)   + 1
               ENDIF
            ENDIF
         ENDDO

         allocate (exchange_w (0:p_np_worker-1))
         DO iworker = 0, p_np_worker-1
            IF (ndata_w(iworker) > 0) THEN
               allocate (exchange_w(iworker)%val (4,ndata_w(iworker)))
            ENDIF
         ENDDO

         ndata_w(:) = 0
         DO ibasin = 1, nbasin
            IF ((addrbasin(1,ibasin) /= -1) .and. (addrbasin(2,ibasin) /= -1) &
               .and. (addrbasin(1,ibasin) /= addrbasin(2,ibasin))) THEN
               ifrom = p_itis_worker(addrbasin(1,ibasin))
               ito   = p_itis_worker(addrbasin(2,ibasin))
               ndata_w(ifrom) = ndata_w(ifrom) + 1
               ndata_w(ito)   = ndata_w(ito)   + 1

               exchange_w(ifrom)%val(:,ndata_w(ifrom)) = &
                  (/addrbasin(1,ibasin), ibasin, addrbasin(2,ibasin), riverdown(ibasin)/)
               exchange_w(ito)%val(:,ndata_w(ito)) = &
                  (/addrbasin(1,ibasin), ibasin, addrbasin(2,ibasin), riverdown(ibasin)/)
            ENDIF
         ENDDO

         DO iworker = 0, p_np_worker-1
            CALL mpi_send (ndata_w(iworker), 1, MPI_INTEGER, &
               p_address_worker(iworker), mpi_tag_size, p_comm_glb, p_err) 
            IF (ndata_w(iworker) > 0) THEN
               CALL mpi_send (exchange_w(iworker)%val, 4*ndata_w(iworker), MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err) 
            ENDIF
         ENDDO
      ENDIF
#endif

      IF (p_is_worker) THEN
#ifdef USEMPI
         CALL mpi_recv (ndata, 1, MPI_INTEGER, p_address_master, mpi_tag_size, p_comm_glb, p_stat, p_err)
         IF (ndata > 0) THEN
            allocate (exchange(4,ndata))
            CALL mpi_recv (exchange, 4*ndata, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
#endif

         IF (numbasin > 0) THEN
            
            allocate (basin_sorted (numbasin))
            allocate (order (numbasin))
            basin_sorted = bindex
            order = (/(ibasin, ibasin = 1, numbasin)/)

            CALL quicksort (numbasin, basin_sorted, order)

            allocate (addrdown (numbasin))
            addrdown(:) = -1

            DO ibasin = 1, numbasin
               IF (riverdown(ibasin) > 0) THEN
                  iloc = find_in_sorted_list1 (riverdown(ibasin), numbasin, basin_sorted)
                  IF (iloc > 0) THEN
                     addrdown(ibasin) = order(iloc)
                  ENDIF
               ENDIF
            ENDDO

#ifdef USEMPI
            IF (ndata > 0) THEN
               nup = count(exchange(3,:) == p_iam_glb)
               ndn = count(exchange(1,:) == p_iam_glb)

               IF (nup > 0) allocate (river_up%iproc (nup))
               IF (nup > 0) allocate (river_up%ups   (nup))
               IF (nup > 0) allocate (river_up%down  (nup))
               IF (nup > 0) allocate (river_up%iloc  (nup))
               
               IF (ndn > 0) allocate (river_dn%iproc (ndn))
               IF (ndn > 0) allocate (river_dn%ups   (ndn))
               IF (ndn > 0) allocate (river_dn%down  (ndn))
               IF (ndn > 0) allocate (river_dn%iloc  (ndn))
               
               iup = 0
               idn = 0
               DO idata = 1, ndata
                  IF (exchange(3,idata) == p_iam_glb) THEN
                     CALL insert_into_sorted_list2 (exchange(2,idata), exchange(1,idata), &
                        iup, river_up%ups, river_up%iproc, iloc)
                  ELSEIF (exchange(1,idata) == p_iam_glb) THEN
                     CALL insert_into_sorted_list2 (exchange(2,idata), exchange(3,idata), &
                        idn, river_dn%ups, river_dn%iproc, iloc)
                  ENDIF
               ENDDO

               DO idata = 1, ndata
                  IF (exchange(3,idata) == p_iam_glb) THEN
                     
                     iloc1 = find_in_sorted_list2 (exchange(2,idata), exchange(1,idata), &
                        nup, river_up%ups, river_up%iproc)

                     river_up%down(iloc1) = exchange(4,idata)

                     iloc2 = find_in_sorted_list1 (exchange(4,idata), numbasin, basin_sorted)
                     river_up%iloc(iloc1) = order(iloc2)

                  ELSEIF (exchange(1,idata) == p_iam_glb) THEN
                     
                     iloc1 = find_in_sorted_list2 (exchange(2,idata), exchange(3,idata), &
                        ndn, river_dn%ups, river_dn%iproc)
                     
                     river_dn%down(iloc1) = exchange(4,idata)
                     
                     iloc2 = find_in_sorted_list1 (exchange(2,idata), numbasin, basin_sorted)
                     river_dn%iloc(iloc1) = order(iloc2)
            
                  ENDIF
               ENDDO

               IF (nup > 0) THEN

                  river_up%nproc = 1
                  DO iup = 2, nup
                     IF (river_up%iproc(iup) /= river_up%iproc(iup-1)) THEN
                        river_up%nproc = river_up%nproc + 1
                     ENDIF
                  ENDDO

                  allocate (river_up%wdsp (river_up%nproc))
                  allocate (river_up%ndata(river_up%nproc))

                  river_up%ndata(:) = 0

                  iproc = 1
                  river_up%wdsp (1) = 0
                  river_up%ndata(1) = 1
                  DO iup = 2, nup
                     IF (river_up%iproc(iup) /= river_up%iproc(iup-1)) THEN
                        iproc = iproc + 1
                        river_up%wdsp (iproc) = iup - 1
                        river_up%ndata(iproc) = 1
                     ELSE
                        river_up%ndata(iproc) = river_up%ndata(iproc) + 1
                     ENDIF
                  ENDDO

               ELSE
                  river_up%nproc = 0
               ENDIF
               
               IF (ndn > 0) THEN

                  river_dn%nproc = 1
                  DO idn = 2, ndn
                     IF (river_dn%iproc(idn) /= river_dn%iproc(idn-1)) THEN
                        river_dn%nproc = river_dn%nproc + 1
                     ENDIF
                  ENDDO

                  allocate (river_dn%wdsp (river_dn%nproc))
                  allocate (river_dn%ndata(river_dn%nproc))

                  river_dn%ndata(:) = 0

                  iproc = 1
                  river_dn%wdsp (1) = 0
                  river_dn%ndata(1) = 1
                  DO idn = 2, ndn
                     IF (river_dn%iproc(idn) /= river_dn%iproc(idn-1)) THEN
                        iproc = iproc + 1
                        river_dn%wdsp (iproc) = idn - 1
                        river_dn%ndata(iproc) = 1
                     ELSE
                        river_dn%ndata(iproc) = river_dn%ndata(iproc) + 1
                     ENDIF
                  ENDDO
               
                  addrdown(river_dn%iloc)= 0

               ELSE
                  river_dn%nproc = 0
               ENDIF
                     
            ENDIF
#endif
         ENDIF
      ENDIF
      
      IF (allocated(bindex      )) deallocate(bindex      )
      IF (allocated(addrbasin   )) deallocate(addrbasin   )
      IF (allocated(ndata_w     )) deallocate(ndata_w     )
      IF (allocated(exchange_w  )) deallocate(exchange_w  )
      IF (allocated(exchange    )) deallocate(exchange    )
      IF (allocated(basin_sorted)) deallocate(basin_sorted)
      IF (allocated(order       )) deallocate(order       )

      IF (p_is_worker) THEN

         IF (numbasin > 0) THEN

            allocate (lakes       (numbasin))
            allocate (riverarea   (numbasin))
            allocate (riverwth    (numbasin))
            allocate (bedelv      (numbasin))
            allocate (handmin     (numbasin))
            allocate (wtsrfelv    (numbasin))
            allocate (riverlen_ds (numbasin))
            allocate (wtsrfelv_ds (numbasin))
            allocate (riverwth_ds (numbasin))
            allocate (bedelv_ds   (numbasin))
            allocate (outletwth   (numbasin))

            DO ibasin = 1, numbasin

               IF (lake_id(ibasin) == 0) THEN

                  riverarea(ibasin) = hillslope_network(ibasin)%area(1)
                  riverwth (ibasin) = riverarea(ibasin) / riverlen(ibasin)

                  ! modify height above nearest drainage data to consider river depth
                  IF (hillslope_network(ibasin)%nhru > 1) THEN
                     hillslope_network(ibasin)%hand(2:) = &
                        hillslope_network(ibasin)%hand(2:) + riverdpth(ibasin)
                  ENDIF

                  wtsrfelv(ibasin) = riverelv(ibasin)
                  bedelv  (ibasin) = riverelv(ibasin) - riverdpth(ibasin)

               ELSEIF (lake_id(ibasin) > 0) THEN
               
                  wtsrfelv(ibasin) = basinelv(ibasin)

                  ps = elm_patch%substt(ibasin)
                  pe = elm_patch%subend(ibasin)

                  bedelv(ibasin) = basinelv(ibasin) - maxval(lakedepth(ps:pe))

                  nsublake = pe - ps + 1
                  lakes(ibasin)%nsub = nsublake

                  allocate (lakes(ibasin)%area0  (nsublake))
                  allocate (lakes(ibasin)%area   (nsublake))
                  allocate (lakes(ibasin)%depth0 (nsublake))
                  allocate (lakes(ibasin)%depth  (nsublake))

                  DO i = 1, nsublake
                     ipatch = i + ps - 1
                     lakes(ibasin)%area(i) = 0
                     DO ipxl = landpatch%ipxstt(ipatch), landpatch%ipxend(ipatch)
                        lakes(ibasin)%area(i) = lakes(ibasin)%area(i) &
                           + 1.0e6 * areaquad ( &
                           pixel%lat_s(mesh(ibasin)%ilat(ipxl)), pixel%lat_n(mesh(ibasin)%ilat(ipxl)), &
                           pixel%lon_w(mesh(ibasin)%ilon(ipxl)), pixel%lon_e(mesh(ibasin)%ilon(ipxl)) )
                     ENDDO
                  ENDDO 

                  ! area data in HRU order
                  lakes(ibasin)%area0 = lakes(ibasin)%area

                  lakes(ibasin)%depth = lakedepth(ps:pe)
                  ! depth data in HRU order
                  lakes(ibasin)%depth0 = lakes(ibasin)%depth

                  allocate (order (1:nsublake))
                  order = (/(i, i = 1, nsublake)/)
            
                  CALL quicksort (nsublake, lakes(ibasin)%depth, order)

                  ! area data in depth order
                  lakes(ibasin)%area = lakes(ibasin)%area(order)
                  
                  ! adjust to be from deepest to shallowest
                  lakes(ibasin)%depth = lakes(ibasin)%depth(nsublake:1:-1)
                  lakes(ibasin)%area  = lakes(ibasin)%area (nsublake:1:-1)

                  allocate (lakes(ibasin)%dep_vol_curve (nsublake))

                  lakes(ibasin)%dep_vol_curve(1) = 0
                  DO i = 2, nsublake
                     lakes(ibasin)%dep_vol_curve(i) = lakes(ibasin)%dep_vol_curve(i-1) &
                        + sum(lakes(ibasin)%area(1:i-1)) * (lakes(ibasin)%depth(i-1)-lakes(ibasin)%depth(i))
                  ENDDO

                  riverlen(ibasin) = 0.

                  deallocate (order)

               ENDIF

               IF (lake_id(ibasin) <= 0) THEN
                  handmin(ibasin) = minval(hillslope_network(ibasin)%hand)
               ENDIF
            ENDDO
         ENDIF

         DO ibasin = 1, numbasin
            IF (addrdown(ibasin) > 0) THEN
               riverlen_ds (ibasin) = riverlen (addrdown(ibasin)) 
               wtsrfelv_ds (ibasin) = wtsrfelv (addrdown(ibasin)) 
               riverwth_ds (ibasin) = riverwth (addrdown(ibasin))
               bedelv_ds   (ibasin) = bedelv   (addrdown(ibasin))
            ELSE
               riverlen_ds (ibasin) = spval
               wtsrfelv_ds (ibasin) = spval
               riverwth_ds (ibasin) = spval
               bedelv_ds   (ibasin) = spval
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL river_data_exchange (SEND_DATA_DOWN_TO_UP, accum = .false., &
            vec_send1 = riverlen, vec_recv1 = riverlen_ds, &
            vec_send2 = wtsrfelv, vec_recv2 = wtsrfelv_ds, &
            vec_send3 = riverwth, vec_recv3 = riverwth_ds, &
            vec_send4 = bedelv  , vec_recv4 = bedelv_ds  )
#endif

         DO ibasin = 1, numbasin
            IF (lake_id(ibasin) < 0) THEN
               bedelv(ibasin) = wtsrfelv_ds(ibasin) + minval(hillslope_network(ibasin)%hand)
            ENDIF
         ENDDO

         DO ibasin = 1, numbasin
            IF (lake_id(ibasin) == 0) THEN
               IF ((to_lake(ibasin)) .or. (riverdown(ibasin) <= 0)) THEN
                  ! river to lake, ocean, inland depression or out of region
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
                     inb = findloc_ud(elementneighbour(ibasin)%glbindex == riverdown(ibasin))
                  ELSE
                     inb = findloc_ud(elementneighbour(ibasin)%glbindex == -9) ! -9 is ocean
                  ENDIF

                  IF (inb <= 0) THEN
                     outletwth(ibasin) = 0
                     IF (riverdown(ibasin) > 0) THEN
                        outletwth(ibasin) = 90.
                     ENDIF
                  ELSE
                     outletwth(ibasin) = elementneighbour(ibasin)%lenbdr(inb)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
      IF (p_is_master) write(*,'(A)') 'Read river network information done.'
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE river_lake_network_init

   ! ----- retrieve river depth from runoff -----
   SUBROUTINE calc_riverdepth_from_runoff ()
      
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_DataType
   USE MOD_NetCDFSerial
   USE MOD_NetCDFBlock
   USE MOD_Block
   USE MOD_Mesh
   USE MOD_Grid
   USE MOD_SpatialMapping
   USE MOD_LandElm
   USE MOD_ElmVector
   USE MOD_ElementNeighbour
   USE MOD_Hydro_IO
   IMPLICIT NONE

   ! Local Variables
   character(len=256) :: file_rnof, file_rivdpt
   type(grid_type)    :: grid_rnof
   type(block_data_real8_2d)  :: f_rnof
   type(spatial_mapping_type) :: mg2p_rnof

   real(r8), allocatable :: bsnrnof(:) , bsndis(:)
   integer,  allocatable :: nups_riv(:), iups_riv(:), b_up2down(:)

   integer :: i, j, ithis, ib, jb, iblkme
   integer :: iwork, mesg(2), isrc, ndata
   real(r8), allocatable :: rcache(:)
   
   real(r8), parameter :: cH_rivdpt   = 0.1
   real(r8), parameter :: pH_rivdpt   = 0.5
   real(r8), parameter :: B0_rivdpt   = 0.0
   real(r8), parameter :: Bmin_rivdpt = 1.0


      file_rnof = trim(DEF_dir_runtime) // '/runoff_clim.nc'

      CALL grid_rnof%define_from_file (file_rnof, 'lat', 'lon')

      CALL mg2p_rnof%build_arealweighted (grid_rnof, landelm)

      IF (p_is_io) THEN
         CALL allocate_block_data (grid_rnof, f_rnof)
         CALL ncio_read_block (file_rnof, 'ro', grid_rnof, f_rnof)

         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)
            DO j = 1, grid_rnof%ycnt(jb)
               DO i = 1, grid_rnof%xcnt(ib)
                  f_rnof%blk(ib,jb)%val(i,j) = max(f_rnof%blk(ib,jb)%val(i,j), 0.)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (p_is_worker) THEN
         IF (numelm > 0) allocate (bsnrnof (numelm))
      ENDIF

      CALL mg2p_rnof%grid2pset (f_rnof, bsnrnof)

      IF (p_is_worker) THEN
         IF (numelm > 0) THEN
            bsnrnof = bsnrnof /24.0/3600.0 ! from m/day to m/s
            DO i = 1, numelm
               ! total runoff in basin, from m/s to m3/s
               bsnrnof(i) = bsnrnof(i) * elementneighbour(i)%myarea
            ENDDO
         ENDIF
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         mesg = (/p_iam_glb, numelm/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err) 
         IF (numelm > 0) THEN
            CALL mpi_send (bsnrnof, numelm, MPI_REAL8, p_address_master, mpi_tag_data, p_comm_glb, p_err) 
         ENDIF
      ENDIF
      
      IF (p_is_master) THEN
         
         allocate (bsnrnof (totalnumelm))

         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndata))

               CALL mpi_recv (rcache, ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
               
               bsnrnof(elm_data_address(p_itis_worker(isrc))%val) = rcache

               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF
      
      CALL mpi_barrier (p_comm_glb, p_err)
#else
      bsnrnof(elm_data_address(0)%val) = bsnrnof
#endif


      IF (p_is_master) THEN

         allocate (nups_riv (totalnumelm))
         allocate (iups_riv (totalnumelm))
         allocate (b_up2down(totalnumelm))

         allocate (bsndis   (totalnumelm))

         nups_riv(:) = 0
         DO i = 1, totalnumelm
            j = riverdown(i)
            IF (j > 0) THEN
               nups_riv(j) = nups_riv(j) + 1
            ENDIF
         ENDDO

         iups_riv(:) = 0
         ithis = 0
         DO i = 1, totalnumelm
            IF (iups_riv(i) == nups_riv(i)) THEN
               ithis = ithis + 1
               b_up2down(ithis) = i

               j = riverdown(i)
               IF (j > 0) THEN
                  iups_riv(j) = iups_riv(j) + 1
                  DO WHILE (iups_riv(j) == nups_riv(j))
                     IF (j < i) THEN
                        ithis = ithis + 1
                        b_up2down(ithis) = j
                     ENDIF
                     j = riverdown(j)
                     IF (j > 0) THEN
                        iups_riv(j) = iups_riv(j) + 1
                     ELSE
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF
            ELSE
               CYCLE
            ENDIF
         ENDDO

         bsndis(:) = 0.
         DO i = 1, totalnumelm
            j = b_up2down(i)
            bsndis(j) = bsndis(j) + bsnrnof(j)
            IF (riverdown(j) > 0) THEN
               bsndis(riverdown(j)) = bsndis(riverdown(j)) + bsndis(j)
            ENDIF
         ENDDO

         allocate (riverdpth (totalnumelm))
         DO i = 1, totalnumelm
            riverdpth(i) = max(cH_rivdpt * (bsndis(i)**pH_rivdpt) + B0_rivdpt, Bmin_rivdpt)
         ENDDO

      ENDIF

      IF (allocated (bsnrnof  )) deallocate(bsnrnof  )
      IF (allocated (bsndis   )) deallocate(bsndis   )
      IF (allocated (nups_riv )) deallocate(nups_riv )
      IF (allocated (iups_riv )) deallocate(iups_riv )
      IF (allocated (b_up2down)) deallocate(b_up2down)

   END SUBROUTINE calc_riverdepth_from_runoff

   ! 
   FUNCTION retrieve_lake_surface_from_volume (this, volume) result(surface)

   IMPLICIT NONE

   class(lake_info_type) :: this
   real(r8), intent(in)  :: volume
   real(r8) :: surface
   
   ! Local Variables
   integer :: i

      IF (volume <= 0) THEN
         surface = 0
         RETURN
      ENDIF

      IF (this%nsub == 1) THEN
         surface = volume / this%area(1)
      ELSE
         i = 1
         DO WHILE (i < this%nsub)
            IF (volume >= this%dep_vol_curve(i+1)) THEN
               i = i + 1
            ELSE
               EXIT
            ENDIF
         ENDDO
         surface = this%depth(1) - this%depth(i) + &
            (volume - this%dep_vol_curve(i)) / sum(this%area(1:i)) 
      ENDIF

   END FUNCTION retrieve_lake_surface_from_volume

   ! 
   FUNCTION retrieve_lake_volume_from_surface (this, surface) result(volume)

   IMPLICIT NONE

   ! Local Variables
   integer :: i

   class(lake_info_type) :: this
   real(r8), intent(in)  :: surface 
   real(r8) :: volume 

      IF (surface <= 0) THEN
         volume = 0
         RETURN
      ENDIF
      
      IF (this%nsub == 1) THEN
         volume = this%area(1) * surface 
      ELSE
         i = 1
         DO WHILE (i < this%nsub)
            IF (surface >= this%depth(1)-this%depth(i+1)) THEN
               i = i + 1
            ELSE
               EXIT
            ENDIF
         ENDDO
         volume = this%dep_vol_curve(i) &
            + (surface - (this%depth(1) - this%depth(i))) * sum(this%area(1:i)) 
      ENDIF

   END FUNCTION retrieve_lake_volume_from_surface


   ! ----------
#ifdef USEMPI
   SUBROUTINE river_data_exchange (direction, accum, &
         vec_send1, vec_recv1, vec_send2, vec_recv2, &
         vec_send3, vec_recv3, vec_send4, vec_recv4 )

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   integer, intent(in) :: direction
   logical, intent(in) :: accum

   real(r8), intent(inout) :: vec_send1(:), vec_recv1(:)
   real(r8), intent(inout), optional :: vec_send2(:), vec_recv2(:)
   real(r8), intent(inout), optional :: vec_send3(:), vec_recv3(:)
   real(r8), intent(inout), optional :: vec_send4(:), vec_recv4(:)

   ! Local Variables
   type(river_sendrecv_type), pointer :: send_pointer
   integer :: nproc_send, ndatasend, idest
   integer,  allocatable :: req_send(:,:)
   real(r8), allocatable :: sendcache1(:)
   real(r8), allocatable :: sendcache2(:)
   real(r8), allocatable :: sendcache3(:)
   real(r8), allocatable :: sendcache4(:)

   type(river_sendrecv_type), pointer :: recv_pointer
   integer :: nproc_recv, ndatarecv, isrc
   integer,  allocatable :: req_recv(:,:)
   real(r8), allocatable :: recvcache1(:)
   real(r8), allocatable :: recvcache2(:)
   real(r8), allocatable :: recvcache3(:)
   real(r8), allocatable :: recvcache4(:)

   integer :: nvec, iproc, i, istt, iend, ndata

      IF (p_is_worker) THEN

         CALL mpi_barrier (p_comm_worker, p_err)

         IF (direction == SEND_DATA_DOWN_TO_UP) THEN
            send_pointer => river_up
            recv_pointer => river_dn
         ELSEIF (direction == SEND_DATA_UP_TO_DOWN) THEN
            send_pointer => river_dn
            recv_pointer => river_up
         ENDIF

         nproc_send = send_pointer%nproc
         IF (nproc_send > 0) THEN

            ndatasend = sum(send_pointer%ndata)
            
            nvec = 1
            allocate (sendcache1(ndatasend))
            DO i = 1, ndatasend
               sendcache1(i) = vec_send1(send_pointer%iloc(i))
            ENDDO

            IF (present(vec_send2) .and. present(vec_recv2)) THEN
               nvec = nvec + 1
               allocate (sendcache2(ndatasend))
               DO i = 1, ndatasend
                  sendcache2(i) = vec_send2(send_pointer%iloc(i))
               ENDDO

               IF (present(vec_send3) .and. present(vec_recv3)) THEN
                  nvec = nvec + 1
                  allocate (sendcache3(ndatasend))
                  DO i = 1, ndatasend
                     sendcache3(i) = vec_send3(send_pointer%iloc(i))
                  ENDDO

                  IF (present(vec_send4) .and. present(vec_recv4)) THEN
                     nvec = nvec + 1
                     allocate (sendcache4(ndatasend))
                     DO i = 1, ndatasend
                        sendcache4(i) = vec_send4(send_pointer%iloc(i))
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF

            allocate (req_send(nvec,nproc_send))
           
            DO iproc = 1, nproc_send
               ndata = send_pointer%ndata(iproc)
               istt  = send_pointer%wdsp (iproc) + 1
               iend  = send_pointer%wdsp (iproc) + ndata
               idest = send_pointer%iproc(istt)

               CALL mpi_isend(sendcache1(istt:iend), ndata, MPI_REAL8, &
                  idest, 101, p_comm_glb, req_send(1,iproc), p_err)
               IF (present(vec_send2) .and. present(vec_recv2)) THEN
                  CALL mpi_isend(sendcache2(istt:iend), ndata, MPI_REAL8, &
                     idest, 102, p_comm_glb, req_send(2,iproc), p_err)
                  IF (present(vec_send3) .and. present(vec_recv3)) THEN
                     CALL mpi_isend(sendcache3(istt:iend), ndata, MPI_REAL8, &
                        idest, 103, p_comm_glb, req_send(3,iproc), p_err)
                     IF (present(vec_send4) .and. present(vec_recv4)) THEN
                        CALL mpi_isend(sendcache4(istt:iend), ndata, MPI_REAL8, &
                           idest, 104, p_comm_glb, req_send(4,iproc), p_err)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

         ENDIF

         nproc_recv = recv_pointer%nproc
         IF (nproc_recv > 0) THEN

            ndatarecv = sum(recv_pointer%ndata)
            
            nvec = 1
            allocate (recvcache1(ndatarecv))
            IF (present(vec_send2) .and. present(vec_recv2)) THEN
               nvec = nvec + 1
               allocate (recvcache2(ndatarecv))
               IF (present(vec_send3) .and. present(vec_recv3)) THEN
                  nvec = nvec + 1
                  allocate (recvcache3(ndatarecv))
                  IF (present(vec_send4) .and. present(vec_recv4)) THEN
                     nvec = nvec + 1
                     allocate (recvcache4(ndatarecv))
                  ENDIF
               ENDIF
            ENDIF
            
            allocate (req_recv(nvec,nproc_recv))
            
            DO iproc = 1, nproc_recv
               ndata = recv_pointer%ndata(iproc)
               istt  = recv_pointer%wdsp(iproc) + 1
               iend  = recv_pointer%wdsp(iproc) + ndata
               isrc  = recv_pointer%iproc(istt)

               CALL mpi_irecv(recvcache1(istt:iend), ndata, MPI_REAL8, &
                  isrc, 101, p_comm_glb, req_recv(1,iproc), p_err)
               IF (present(vec_send2) .and. present(vec_recv2)) THEN
                  CALL mpi_irecv(recvcache2(istt:iend), ndata, MPI_REAL8, &
                     isrc, 102, p_comm_glb, req_recv(2,iproc), p_err)
                  IF (present(vec_send3) .and. present(vec_recv3)) THEN
                     CALL mpi_irecv(recvcache3(istt:iend), ndata, MPI_REAL8, &
                        isrc, 103, p_comm_glb, req_recv(3,iproc), p_err)
                     IF (present(vec_send4) .and. present(vec_recv4)) THEN
                        CALL mpi_irecv(recvcache4(istt:iend), ndata, MPI_REAL8, &
                           isrc, 104, p_comm_glb, req_recv(4,iproc), p_err)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

         ENDIF

         IF (nproc_recv > 0) THEN

            CALL mpi_waitall(nvec*nproc_recv, req_recv, MPI_STATUSES_IGNORE, p_err)
            ! write(*,*) 'p error', p_err

            IF (accum) THEN
               DO i = 1, ndatarecv
                  vec_recv1(recv_pointer%iloc(i)) = &
                     vec_recv1(recv_pointer%iloc(i)) + recvcache1(i)
               ENDDO

               IF (present(vec_send2) .and. present(vec_recv2)) THEN
                  DO i = 1, ndatarecv
                     vec_recv2(recv_pointer%iloc(i)) = &
                        vec_recv2(recv_pointer%iloc(i)) + recvcache2(i)
                  ENDDO

                  IF (present(vec_send3) .and. present(vec_recv3)) THEN
                     DO i = 1, ndatarecv
                        vec_recv3(recv_pointer%iloc(i)) = &
                           vec_recv3(recv_pointer%iloc(i)) + recvcache3(i)
                     ENDDO

                     IF (present(vec_send4) .and. present(vec_recv4)) THEN
                        DO i = 1, ndatarecv
                           vec_recv4(recv_pointer%iloc(i)) = &
                              vec_recv4(recv_pointer%iloc(i)) + recvcache4(i)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
            ELSE
               DO i = 1, ndatarecv
                  vec_recv1(recv_pointer%iloc(i)) = recvcache1(i)
               ENDDO

               IF (present(vec_send2) .and. present(vec_recv2)) THEN
                  DO i = 1, ndatarecv
                     vec_recv2(recv_pointer%iloc(i)) = recvcache2(i)
                  ENDDO

                  IF (present(vec_send3) .and. present(vec_recv3)) THEN
                     DO i = 1, ndatarecv
                        vec_recv3(recv_pointer%iloc(i)) = recvcache3(i)
                     ENDDO

                     IF (present(vec_send4) .and. present(vec_recv4)) THEN
                        DO i = 1, ndatarecv
                           vec_recv4(recv_pointer%iloc(i)) = recvcache4(i)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (nproc_send > 0) THEN
            CALL mpi_waitall(nvec*nproc_send, req_send, MPI_STATUSES_IGNORE, p_err)
         ENDIF

         IF (allocated(req_send  )) deallocate(req_send)
         IF (allocated(sendcache1)) deallocate(sendcache1)
         IF (allocated(sendcache2)) deallocate(sendcache2)
         IF (allocated(sendcache3)) deallocate(sendcache3)
         IF (allocated(sendcache4)) deallocate(sendcache4)
         
         IF (allocated(req_recv  )) deallocate(req_recv)
         IF (allocated(recvcache1)) deallocate(recvcache1)
         IF (allocated(recvcache2)) deallocate(recvcache2)
         IF (allocated(recvcache3)) deallocate(recvcache3)
         IF (allocated(recvcache4)) deallocate(recvcache4)

         CALL mpi_barrier (p_comm_worker, p_err)

      ENDIF

   END SUBROUTINE river_data_exchange
#endif

   ! ----------
   SUBROUTINE river_lake_network_final ()

   IMPLICIT NONE

   ! Local Variables
   integer :: ilake

      IF (allocated(lake_id  )) deallocate(lake_id  )
      IF (allocated(riverlen )) deallocate(riverlen )
      IF (allocated(riverelv )) deallocate(riverelv )
      IF (allocated(riverarea)) deallocate(riverarea)
      IF (allocated(riverwth )) deallocate(riverwth )
      IF (allocated(riverdpth)) deallocate(riverdpth)
      IF (allocated(basinelv )) deallocate(basinelv )
      IF (allocated(bedelv   )) deallocate(bedelv   )
      IF (allocated(handmin  )) deallocate(handmin  )
      IF (allocated(wtsrfelv )) deallocate(wtsrfelv )
      IF (allocated(riverdown)) deallocate(riverdown)
      IF (allocated(addrdown )) deallocate(addrdown )
      IF (allocated(to_lake  )) deallocate(to_lake  )

      IF (allocated(riverlen_ds))  deallocate(riverlen_ds)
      IF (allocated(wtsrfelv_ds))  deallocate(wtsrfelv_ds)
      IF (allocated(riverwth_ds))  deallocate(riverwth_ds)
      IF (allocated(bedelv_ds  ))  deallocate(bedelv_ds  )
      IF (allocated(outletwth  ))  deallocate(outletwth  )

      IF (allocated(lakes)) THEN
         DO ilake = 1, size(lakes)
            IF (allocated(lakes(ilake)%area0        )) deallocate(lakes(ilake)%area0        )
            IF (allocated(lakes(ilake)%area         )) deallocate(lakes(ilake)%area         )
            IF (allocated(lakes(ilake)%depth0       )) deallocate(lakes(ilake)%depth0       )
            IF (allocated(lakes(ilake)%depth        )) deallocate(lakes(ilake)%depth        )
            IF (allocated(lakes(ilake)%dep_vol_curve)) deallocate(lakes(ilake)%dep_vol_curve)
         ENDDO

         deallocate(lakes)
      ENDIF

   END SUBROUTINE river_lake_network_final

   ! ---------
   SUBROUTINE river_sendrecv_free_mem (this)
      
   IMPLICIT NONE
   type(river_sendrecv_type) :: this

      IF (allocated(this%iproc)) deallocate(this%iproc)
      IF (allocated(this%wdsp )) deallocate(this%wdsp )
      IF (allocated(this%ndata)) deallocate(this%ndata)
      IF (allocated(this%ups  )) deallocate(this%ups  )
      IF (allocated(this%down )) deallocate(this%down )
      IF (allocated(this%iloc )) deallocate(this%iloc )

   END SUBROUTINE river_sendrecv_free_mem

END MODULE MOD_Catch_RiverLakeNetwork
#endif
