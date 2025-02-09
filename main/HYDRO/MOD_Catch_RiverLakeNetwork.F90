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
   USE MOD_Pixelset
   USE MOD_Catch_BasinNetwork
   USE MOD_Catch_HillslopeNetwork
   USE MOD_Catch_Vars_TimeVariables
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
   logical, allocatable :: to_lake    (:)

   real(r8), allocatable :: riverlen_ds (:)
   real(r8), allocatable :: wtsrfelv_ds (:)
   real(r8), allocatable :: riverwth_ds (:)
   real(r8), allocatable :: bedelv_ds   (:)
   
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
      final :: lake_info_free_mem
   END type lake_info_type

   ! -- lake information --
   integer, allocatable :: lake_id (:)
   type(lake_info_type), allocatable :: lakeinfo (:)
   
   ! -- information of HRU in basin --
   type(hillslope_network_type), pointer :: hillslope_basin (:)

   type(basin_pushdata_type), target :: river_iam_dn
   type(basin_pushdata_type), target :: river_iam_up

CONTAINS
   
   ! ----------
   SUBROUTINE river_lake_network_init ()

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   USE MOD_Pixel
   USE MOD_LandElm
   USE MOD_LandHRU
   USE MOD_LandPatch
   USE MOD_ElementNeighbour
   USE MOD_DataType
   USE MOD_Utils
   USE MOD_UserDefFun
   USE MOD_Vars_TimeInvariants, only : lakedepth
   IMPLICIT NONE

   ! Local Variables
   character(len=256) :: river_file, rivdpt_file
   logical :: use_calc_rivdpt

   integer :: totalnumbasin, ibasin, inb
   integer :: iworker, mesg(4), isrc, idest, iproc
   integer :: nrecv, irecv, ifrom, ito, iup, idn
   integer :: ndata, idata, ip, nup, ndn
   integer :: iloc, iloc1, iloc2
   integer :: ielm, i, j, ithis, nave
   

   integer , allocatable :: icache (:)
   real(r8), allocatable :: rcache (:)
   logical , allocatable :: lcache (:)
   
   type(pointer_int32_2d), allocatable :: datapush_w (:)
   integer, allocatable :: datapush(:,:)
   integer, allocatable :: bindex(:), addrbasin(:), addrdown(:), nelm_wrk(:), paddr(:), ndata_w(:)
   integer, allocatable :: basin_sorted(:), basin_order(:), order (:)
   integer, allocatable :: river_up_ups(:), river_up_paddr(:), river_dn_ups(:), river_dn_paddr(:)

   ! for lakes
   integer :: ps, pe, nsublake, hs, he, ihru, ipxl
   integer,  allocatable :: lake_id_elm (:)
   integer , allocatable :: lakedown_id_elm(:), lakedown_id_bsn (:)
   real(r8), allocatable :: lakedepth_hru  (:), lakedepth_bsnhru(:)
   real(r8), allocatable :: lakeoutlet_elm (:), lakeoutlet_bsn  (:)
   real(r8), allocatable :: lakearea_hru   (:), lakearea_bsnhru (:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      use_calc_rivdpt = DEF_USE_EstimatedRiverDepth
      river_file      = DEF_CatchmentMesh_data 

      ! step 1: read in parameters from file.
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
         
         totalnumbasin = size(riverdown)
         allocate (to_lake (totalnumbasin))
         to_lake = .false.
         DO i = 1, totalnumbasin
            IF (riverdown(i) > 0) THEN
               to_lake(i) = lake_id(riverdown(i)) > 0 
            ENDIF
         ENDDO

      ENDIF
         
      ! step 2: Estimate river depth by using runoff data.
      IF (use_calc_rivdpt) THEN
         CALL calc_riverdepth_from_runoff ()
      ENDIF

#ifdef USEMPI
      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN
               
               allocate (bindex (nrecv))
               allocate (icache (nrecv))
               allocate (rcache (nrecv))
               allocate (lcache (nrecv))

               CALL mpi_recv (bindex, nrecv, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               idest = isrc

               icache = lake_id(bindex)
               CALL mpi_send (icache, nrecv, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err) 

               icache = riverdown(bindex)
               CALL mpi_send (icache, nrecv, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err) 

               lcache = to_lake(bindex)
               CALL mpi_send (lcache, nrecv, MPI_LOGICAL, idest, mpi_tag_data, p_comm_glb, p_err) 

               rcache = riverlen(bindex)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err) 

               rcache = riverelv(bindex)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err) 

               rcache = riverdpth(bindex)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err) 

               rcache = basinelv(bindex)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err) 

               deallocate (bindex)
               deallocate (icache)
               deallocate (rcache)
               deallocate (lcache)

            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN
               
         mesg(1:2) = (/p_iam_glb, numbasin/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err) 

         IF (numbasin > 0) THEN
            CALL mpi_send (basinindex, numbasin, MPI_INTEGER, &
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

      ENDIF 
      
      CALL mpi_barrier (p_comm_glb, p_err)
#else
      IF (numbasin > 0) THEN

         lake_id   = lake_id  (basinindex)
         riverdown = riverdown(basinindex)
         to_lake   = to_lake  (basinindex)
         riverlen  = riverlen (basinindex)
         riverelv  = riverelv (basinindex)
         riverdpth = riverdpth(basinindex)
         basinelv  = basinelv (basinindex)

      ENDIF
#endif

#ifdef USEMPI
      ! get address of basins
      IF (p_is_master) THEN
         
         allocate (addrbasin (totalnumbasin)); addrbasin(:) = -1

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN
               allocate (bindex (nrecv))

               CALL mpi_recv (bindex, nrecv, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               addrbasin(bindex) = isrc

               deallocate(bindex)
            ENDIF

         ENDDO

      ELSEIF (p_is_worker) THEN
      
         mesg(1:2) = (/p_iam_glb, numbasin/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err) 

         IF (numbasin > 0) THEN
            CALL mpi_send (basinindex, numbasin, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_err) 
         ENDIF

      ENDIF
      
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
      IF (p_is_master) THEN

         allocate (addrdown  (totalnumbasin))
         allocate (ndata_w (0:p_np_worker-1))

         ndata_w(:) = 0
         DO ibasin = 1, totalnumbasin
            IF (riverdown(ibasin) >= 1) THEN
               addrdown(ibasin) = addrbasin(riverdown(ibasin))
               IF (addrbasin(ibasin) /= addrdown(ibasin)) THEN
                  ifrom = p_itis_worker(addrbasin(ibasin))
                  ito   = p_itis_worker(addrdown(ibasin))
                  ndata_w(ifrom) = ndata_w(ifrom) + 1
                  ndata_w(ito)   = ndata_w(ito)   + 1
               ENDIF
            ENDIF
         ENDDO

         allocate (datapush_w (0:p_np_worker-1))
         DO iworker = 0, p_np_worker-1
            IF (ndata_w(iworker) > 0) THEN
               allocate (datapush_w(iworker)%val (4,ndata_w(iworker)))
            ENDIF
         ENDDO

         ndata_w(:) = 0
         DO ibasin = 1, totalnumbasin
            IF ((riverdown(ibasin) >= 1) .and. (addrbasin(ibasin) /= addrdown(ibasin))) THEN
               ifrom = p_itis_worker(addrbasin(ibasin))
               ito   = p_itis_worker(addrdown(ibasin))
               ndata_w(ifrom) = ndata_w(ifrom) + 1
               ndata_w(ito)   = ndata_w(ito)   + 1

               datapush_w(ifrom)%val(:,ndata_w(ifrom)) = &
                  (/addrbasin(ibasin), ibasin, addrdown(ibasin), riverdown(ibasin)/)
               datapush_w(ito)%val(:,ndata_w(ito)) = &
                  (/addrbasin(ibasin), ibasin, addrdown(ibasin), riverdown(ibasin)/)
            ENDIF
         ENDDO

         DO iworker = 0, p_np_worker-1
            CALL mpi_send (ndata_w(iworker), 1, MPI_INTEGER, &
               p_address_worker(iworker), mpi_tag_size, p_comm_glb, p_err) 
            IF (ndata_w(iworker) > 0) THEN
               CALL mpi_send (datapush_w(iworker)%val, 4*ndata_w(iworker), MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err) 
            ENDIF
         ENDDO

         deallocate (addrdown  )
         deallocate (ndata_w   )
         deallocate (datapush_w)

      ENDIF
#endif

      IF (p_is_worker) THEN
#ifdef USEMPI
         CALL mpi_recv (ndata, 1, MPI_INTEGER, p_address_master, mpi_tag_size, p_comm_glb, p_stat, p_err)
         IF (ndata > 0) THEN
            allocate (datapush(4,ndata))
            CALL mpi_recv (datapush, 4*ndata, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
#endif

         IF (numbasin > 0) THEN
            allocate (basin_sorted (numbasin))
            allocate (basin_order  (numbasin))
            basin_sorted = basinindex
            basin_order  = (/(ibasin, ibasin = 1, numbasin)/)

            CALL quicksort (numbasin, basin_sorted, basin_order)
         ENDIF

         river_iam_up%nself = 0
         river_iam_dn%nself = 0
         river_iam_up%nproc = 0
         river_iam_dn%nproc = 0

         IF (numbasin > 0) THEN
            
            DO ibasin = 1, numbasin
               IF (riverdown(ibasin) > 0) THEN
                  iloc = find_in_sorted_list1 (riverdown(ibasin), numbasin, basin_sorted)
                  IF (iloc > 0) THEN
                     river_iam_up%nself = river_iam_up%nself + 1
                  ENDIF
               ENDIF
            ENDDO
            
            IF (river_iam_up%nself > 0) THEN
               
               allocate (river_iam_up%iself (river_iam_up%nself))

               river_iam_dn%nself = river_iam_up%nself
               allocate (river_iam_dn%iself (river_iam_dn%nself))

               idata = 0
               DO ibasin = 1, numbasin
                  IF (riverdown(ibasin) > 0) THEN
                     iloc = find_in_sorted_list1 (riverdown(ibasin), numbasin, basin_sorted)
                     IF (iloc > 0) THEN
                        idata = idata + 1
                        river_iam_up%iself(idata) = ibasin
                        river_iam_dn%iself(idata) = basin_order(iloc)
                     ENDIF
                  ENDIF
               ENDDO

            ENDIF

#ifdef USEMPI
            IF (ndata > 0) THEN

               nup = count(datapush(3,:) == p_iam_glb)
               ndn = count(datapush(1,:) == p_iam_glb)

               IF (nup > 0) allocate (river_up_paddr (nup))
               IF (nup > 0) allocate (river_up_ups   (nup))
               IF (ndn > 0) allocate (river_dn_paddr (ndn))
               IF (ndn > 0) allocate (river_dn_ups   (ndn))
               
               iup = 0
               idn = 0
               DO idata = 1, ndata
                  IF (datapush(3,idata) == p_iam_glb) THEN
                     CALL insert_into_sorted_list2 (datapush(2,idata), datapush(1,idata), &
                        iup, river_up_ups, river_up_paddr, iloc)
                  ELSEIF (datapush(1,idata) == p_iam_glb) THEN
                     CALL insert_into_sorted_list2 (datapush(2,idata), datapush(3,idata), &
                        idn, river_dn_ups, river_dn_paddr, iloc)
                  ENDIF
               ENDDO
               
               IF (nup > 0) allocate (river_iam_dn%ipush (nup))
               IF (ndn > 0) allocate (river_iam_up%ipush (ndn))

               DO idata = 1, ndata
                  IF (datapush(3,idata) == p_iam_glb) THEN
                     
                     iloc1 = find_in_sorted_list2 (datapush(2,idata), datapush(1,idata), &
                        nup, river_up_ups, river_up_paddr)
                     iloc2 = find_in_sorted_list1 (datapush(4,idata), numbasin, basin_sorted)

                     river_iam_dn%ipush(iloc1) = basin_order(iloc2)

                  ELSEIF (datapush(1,idata) == p_iam_glb) THEN
                     
                     iloc1 = find_in_sorted_list2 (datapush(2,idata), datapush(3,idata), &
                        ndn, river_dn_ups, river_dn_paddr)
                     iloc2 = find_in_sorted_list1 (datapush(2,idata), numbasin, basin_sorted)
                     
                     river_iam_up%ipush(iloc1) = basin_order(iloc2)
            
                  ENDIF
               ENDDO
         

               IF (nup > 0) THEN

                  DO iup = 1, nup
                     IF (iup == 1) THEN
                        river_iam_dn%nproc = 1
                     ELSEIF (river_up_paddr(iup) /= river_up_paddr(iup-1)) THEN
                        river_iam_dn%nproc = river_iam_dn%nproc + 1
                     ENDIF
                  ENDDO
                  
                  allocate (river_iam_dn%paddr(river_iam_dn%nproc))
                  allocate (river_iam_dn%ndata(river_iam_dn%nproc))
                  
                  DO iup = 1, nup
                     IF (iup == 1) THEN
                        ip = 1
                        river_iam_dn%paddr(ip) = river_up_paddr(iup)
                        river_iam_dn%ndata(ip) = 1
                     ELSEIF (river_up_paddr(iup) /= river_up_paddr(iup-1)) THEN
                        ip = ip + 1
                        river_iam_dn%paddr(ip) = river_up_paddr(iup)
                        river_iam_dn%ndata(ip) = 1
                     ELSE
                        river_iam_dn%ndata(ip) = river_iam_dn%ndata(ip) + 1
                     ENDIF
                  ENDDO

               ENDIF
               

               IF (ndn > 0) THEN

                  DO idn = 1, ndn
                     IF (idn == 1) THEN
                        river_iam_up%nproc = 1
                     ELSEIF (river_dn_paddr(idn) /= river_dn_paddr(idn-1)) THEN
                        river_iam_up%nproc = river_iam_up%nproc + 1
                     ENDIF
                  ENDDO
                  
                  allocate (river_iam_up%paddr(river_iam_up%nproc))
                  allocate (river_iam_up%ndata(river_iam_up%nproc))
                  
                  DO idn = 1, ndn
                     IF (idn == 1) THEN
                        ip = 1
                        river_iam_up%paddr(ip) = river_dn_paddr(idn)
                        river_iam_up%ndata(ip) = 1
                     ELSEIF (river_dn_paddr(idn) /= river_dn_paddr(idn-1)) THEN
                        ip = ip + 1
                        river_iam_up%paddr(ip) = river_dn_paddr(idn)
                        river_iam_up%ndata(ip) = 1
                     ELSE
                        river_iam_up%ndata(ip) = river_iam_up%ndata(ip) + 1
                     ENDIF
                  ENDDO

               ENDIF
               
               IF (nup > 0) deallocate (river_up_paddr)
               IF (nup > 0) deallocate (river_up_ups  )
               IF (ndn > 0) deallocate (river_dn_paddr)
               IF (ndn > 0) deallocate (river_dn_ups  )
      
               deallocate(datapush)
                     
            ENDIF
#endif
         ENDIF
      ENDIF
      
      IF (allocated(basin_sorted  )) deallocate(basin_sorted  )
      IF (allocated(basin_order   )) deallocate(basin_order   )


      CALL hillslope_network_init (numbasin, basinindex, hillslope_basin)

      IF (p_is_worker) THEN

         IF (numelm > 0)    allocate (lake_id_elm     (numelm))
         IF (numhru > 0)    allocate (lakedepth_hru   (numhru))
         IF (numbsnhru > 0) allocate (lakedepth_bsnhru(numbsnhru))
         IF (numelm > 0)    allocate (lakedown_id_elm (numelm))
         IF (numbasin > 0)  allocate (lakedown_id_bsn (numbasin))
         IF (numelm > 0)    allocate (lakeoutlet_elm  (numelm))
         IF (numbasin > 0)  allocate (lakeoutlet_bsn  (numbasin))
         IF (numhru > 0)    allocate (lakearea_hru    (numhru))
         IF (numbsnhru > 0) allocate (lakearea_bsnhru (numbsnhru))

         DO ibasin = 1, numbasin

            lakedown_id_bsn(ibasin) = 0

            IF ((lake_id(ibasin) /= 0) .and. (to_lake(ibasin))) THEN
               ! lake to lake .or. lake catchment to lake
               lakedown_id_bsn(ibasin) = riverdown(ibasin)
            ENDIF
            IF ((lake_id(ibasin) > 0) .and. (riverdown(ibasin) == 0)) THEN
               ! lake to ocean
               lakedown_id_bsn(ibasin) = -9    ! -9 is ocean
            ENDIF
         ENDDO

      ENDIF

      CALL worker_push_data (iam_bsn, iam_elm, .false., lake_id, lake_id_elm)
      CALL worker_push_data (iam_bsn, iam_elm, .false., lakedown_id_bsn, lakedown_id_elm)

      IF (p_is_worker) THEN

         lakedepth_hru  = 0.
         lakeoutlet_elm = 0.
         lakearea_hru   = 0.

         DO ielm = 1, numelm
            IF (lake_id_elm(ielm) > 0) THEN
               hs = elm_hru%substt(ielm)
               he = elm_hru%subend(ielm)
               DO ihru = hs, he
                  ps = hru_patch%substt(ihru)
                  pe = hru_patch%subend(ihru)

                  lakedepth_hru(ihru) = maxval(lakedepth(ps:pe))

                  lakearea_hru(ihru) = 0.
                  DO ipxl = landhru%ipxstt(ihru), landhru%ipxend(ihru)
                     lakearea_hru(ihru) = lakearea_hru(ihru) &
                        + 1.0e6 * areaquad ( &
                        pixel%lat_s(mesh(ielm)%ilat(ipxl)), pixel%lat_n(mesh(ielm)%ilat(ipxl)), &
                        pixel%lon_w(mesh(ielm)%ilon(ipxl)), pixel%lon_e(mesh(ielm)%ilon(ipxl)) )
                  ENDDO
               ENDDO 
            ENDIF

            IF (lakedown_id_elm(ielm) /= 0) THEN
               inb = findloc_ud(elementneighbour(ielm)%glbindex == lakedown_id_elm(ielm))
               IF (inb > 0) lakeoutlet_elm(ielm) = elementneighbour(ielm)%lenbdr(inb)
            ENDIF
         ENDDO
      ENDIF
      
      CALL worker_push_data (iam_elm, iam_bsn, .false., lakeoutlet_elm, lakeoutlet_bsn)

      CALL worker_push_subset_data (iam_elm, iam_bsn, elm_hru, basin_hru, lakedepth_hru, lakedepth_bsnhru)
      CALL worker_push_subset_data (iam_elm, iam_bsn, elm_hru, basin_hru, lakearea_hru,  lakearea_bsnhru )

      IF (allocated (lake_id_elm    )) deallocate (lake_id_elm    )
      IF (allocated (lakedepth_hru  )) deallocate (lakedepth_hru  )
      IF (allocated (lakedown_id_elm)) deallocate (lakedown_id_elm)
      IF (allocated (lakedown_id_bsn)) deallocate (lakedown_id_bsn)
      IF (allocated (lakeoutlet_elm )) deallocate (lakeoutlet_elm )
      IF (allocated (lakearea_hru   )) deallocate (lakearea_hru   )

      IF (p_is_worker) THEN

         IF (numbasin > 0) THEN

            allocate (lakeinfo    (numbasin))
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

                  riverarea(ibasin) = hillslope_basin(ibasin)%area(1)
                  riverwth (ibasin) = riverarea(ibasin) / riverlen(ibasin)

                  ! modify height above nearest drainage data to consider river depth
                  IF (hillslope_basin(ibasin)%nhru > 1) THEN
                     hillslope_basin(ibasin)%hand(2:) = &
                        hillslope_basin(ibasin)%hand(2:) + riverdpth(ibasin)
                  ENDIF

                  wtsrfelv(ibasin) = riverelv(ibasin)
                  bedelv  (ibasin) = riverelv(ibasin) - riverdpth(ibasin)

               ELSEIF (lake_id(ibasin) > 0) THEN
               
                  hs = basin_hru%substt(ibasin)
                  he = basin_hru%subend(ibasin)

                  wtsrfelv(ibasin) = basinelv(ibasin)

                  bedelv(ibasin) = basinelv(ibasin) - minval(lakedepth_bsnhru(hs:he))

                  nsublake = he - hs + 1
                  lakeinfo(ibasin)%nsub = nsublake

                  allocate (lakeinfo(ibasin)%area0  (nsublake))
                  allocate (lakeinfo(ibasin)%area   (nsublake))
                  allocate (lakeinfo(ibasin)%depth0 (nsublake))
                  allocate (lakeinfo(ibasin)%depth  (nsublake))

                  lakeinfo(ibasin)%area  = lakearea_bsnhru (hs:he)
                  lakeinfo(ibasin)%depth = lakedepth_bsnhru(hs:he)

                  ! area data in HRU order
                  lakeinfo(ibasin)%area0 = lakeinfo(ibasin)%area

                  ! depth data in HRU order
                  lakeinfo(ibasin)%depth0 = lakeinfo(ibasin)%depth

                  allocate (order (1:nsublake))
                  order = (/(i, i = 1, nsublake)/)
            
                  CALL quicksort (nsublake, lakeinfo(ibasin)%depth, order)

                  ! area data in depth order
                  lakeinfo(ibasin)%area = lakeinfo(ibasin)%area(order)
                  
                  ! adjust to be from deepest to shallowest
                  lakeinfo(ibasin)%depth = lakeinfo(ibasin)%depth(nsublake:1:-1)
                  lakeinfo(ibasin)%area  = lakeinfo(ibasin)%area (nsublake:1:-1)

                  allocate (lakeinfo(ibasin)%dep_vol_curve (nsublake))

                  lakeinfo(ibasin)%dep_vol_curve(1) = 0
                  DO i = 2, nsublake
                     lakeinfo(ibasin)%dep_vol_curve(i) = &
                        lakeinfo(ibasin)%dep_vol_curve(i-1) &
                        + sum(lakeinfo(ibasin)%area(1:i-1)) &
                        * (lakeinfo(ibasin)%depth(i-1)-lakeinfo(ibasin)%depth(i))
                  ENDDO

                  riverlen(ibasin) = 0.

                  deallocate (order)

               ENDIF

               IF (lake_id(ibasin) <= 0) THEN
                  handmin(ibasin) = minval(hillslope_basin(ibasin)%hand)
               ENDIF
            ENDDO
         ENDIF

         CALL worker_push_data (river_iam_dn, river_iam_up, .false., riverlen, riverlen_ds)
         CALL worker_push_data (river_iam_dn, river_iam_up, .false., wtsrfelv, wtsrfelv_ds)
         CALL worker_push_data (river_iam_dn, river_iam_up, .false., riverwth, riverwth_ds)
         CALL worker_push_data (river_iam_dn, river_iam_up, .false., bedelv  , bedelv_ds  )

         DO ibasin = 1, numbasin
            IF (lake_id(ibasin) < 0) THEN
               bedelv(ibasin) = wtsrfelv_ds(ibasin) + minval(hillslope_basin(ibasin)%hand)
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
               ELSEIF (to_lake(ibasin)) THEN
                  ! lake to lake .or. lake catchment to lake
                  outletwth(ibasin) = lakeoutlet_bsn(ibasin)
               ELSEIF (riverdown(ibasin) == 0) THEN
                  ! lake to ocean
                  outletwth(ibasin) = lakeoutlet_bsn(ibasin)
               ENDIF
            ENDIF
         ENDDO

      ENDIF

      IF (allocated (lakedepth_bsnhru)) deallocate (lakedepth_bsnhru)
      IF (allocated (lakeoutlet_bsn  )) deallocate (lakeoutlet_bsn  )
      IF (allocated (lakearea_bsnhru )) deallocate (lakearea_bsnhru )

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
      IF (p_is_master) write(*,'(A)') 'Building river network information done.'
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE river_lake_network_init

   ! ----- retrieve river depth from runoff -----
   SUBROUTINE calc_riverdepth_from_runoff ()
      
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_DataType
   USE MOD_Utils
   USE MOD_NetCDFSerial
   USE MOD_NetCDFBlock
   USE MOD_Pixel
   USE MOD_Block
   USE MOD_Mesh
   USE MOD_Grid
   USE MOD_SpatialMapping
   USE MOD_LandElm
   USE MOD_ElmVector
   IMPLICIT NONE

   ! Local Variables
   character(len=256) :: file_rnof, file_rivdpt
   type(grid_type)    :: grid_rnof
   type(block_data_real8_2d)  :: f_rnof
   type(spatial_mapping_type) :: mg2p_rnof

   real(r8), allocatable :: bsnrnof(:) , bsndis(:)
   integer,  allocatable :: nups_riv(:), iups_riv(:), b_up2down(:)

   integer :: i, j, ithis, ib, jb, iblkme, ipxl
   integer :: iwork, mesg(2), isrc, ndata
   real(r8), allocatable :: rcache(:)
   real(r8) :: myarea
   
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
               myarea = 0.
               DO ipxl = landelm%ipxstt(i), landelm%ipxend(i)
                  myarea = myarea + 1.0e6 * areaquad ( &
                     pixel%lat_s(mesh(i)%ilat(ipxl)), pixel%lat_n(mesh(i)%ilat(ipxl)), &
                     pixel%lon_w(mesh(i)%ilon(ipxl)), pixel%lon_e(mesh(i)%ilon(ipxl)) )
               ENDDO

               ! total runoff in basin, from m/s to m3/s
               bsnrnof(i) = bsnrnof(i) * myarea
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
               DO WHILE (j > 0) 

                  iups_riv(j) = iups_riv(j) + 1

                  IF (iups_riv(j) == nups_riv(j)) THEN
                     ithis = ithis + 1
                     b_up2down(ithis) = j
                     j = riverdown(j)
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
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
      IF (allocated(to_lake  )) deallocate(to_lake  )

      IF (allocated(riverlen_ds))  deallocate(riverlen_ds)
      IF (allocated(wtsrfelv_ds))  deallocate(wtsrfelv_ds)
      IF (allocated(riverwth_ds))  deallocate(riverwth_ds)
      IF (allocated(bedelv_ds  ))  deallocate(bedelv_ds  )
      IF (allocated(outletwth  ))  deallocate(outletwth  )
      
      IF (allocated(lakeinfo)) deallocate(lakeinfo)

      IF (associated(hillslope_basin)) deallocate(hillslope_basin)

   END SUBROUTINE river_lake_network_final

   ! ---------
   SUBROUTINE lake_info_free_mem (this)

   IMPLICIT NONE
   type(lake_info_type) :: this

      IF (allocated(this%area0 )) deallocate (this%area0 )
      IF (allocated(this%area  )) deallocate (this%area  )
      IF (allocated(this%depth0)) deallocate (this%depth0)
      IF (allocated(this%depth )) deallocate (this%depth )
      IF (allocated(this%dep_vol_curve)) deallocate (this%dep_vol_curve)

   END SUBROUTINE lake_info_free_mem

END MODULE MOD_Catch_RiverLakeNetwork
#endif
