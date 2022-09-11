#include <define.h>

MODULE mod_landbasin

   USE precision
   IMPLICIT NONE

   ! ---- data types ----
   TYPE :: irregular_basin_type 
      
      INTEGER :: indx
      INTEGER :: xblk, yblk

      INTEGER :: npxl
      INTEGER, allocatable :: ilon(:)
      INTEGER, allocatable :: ilat(:)

   END TYPE irregular_basin_type

   ! ---- Instance ----
   INTEGER :: numbasin
   TYPE (irregular_basin_type), allocatable :: landbasin (:)

   INTEGER, allocatable :: nbasin_blk(:,:)
   
CONTAINS

   ! -------
   SUBROUTINE copy_basin (basin_from, basin_to)

      IMPLICIT NONE
      TYPE (irregular_basin_type), intent(in)  :: basin_from
      TYPE (irregular_basin_type), intent(out) :: basin_to

      basin_to%indx = basin_from%indx  
      basin_to%npxl = basin_from%npxl 
      basin_to%xblk = basin_from%xblk 
      basin_to%yblk = basin_from%yblk 

      IF (allocated(basin_to%ilat)) deallocate(basin_to%ilat)
      IF (allocated(basin_to%ilon)) deallocate(basin_to%ilon)

      allocate (basin_to%ilat (basin_to%npxl))
      allocate (basin_to%ilon (basin_to%npxl))
      basin_to%ilon = basin_from%ilon
      basin_to%ilat = basin_from%ilat 

   END SUBROUTINE copy_basin
   
   ! --------------------------------
   SUBROUTINE landbasin_build (gbasin)

      USE precision
      USE mod_namelist
      USE spmd_task
      USE ncio_block
      USE mod_block
      USE mod_pixel
      USE mod_grid
      USE mod_utils
      USE mod_data_type
      USE mod_hydro_data

      IMPLICIT NONE

      TYPE(grid_type), intent(in) :: gbasin

      ! Local Variables 
      CHARACTER(len=256) :: file_lbasin
      CHARACTER(len=256) :: varname
      TYPE(block_data_int32_2d) :: databasin

      INTEGER  :: iworker
      INTEGER  :: nbasin, iu, ju
      INTEGER  :: iblkme, iblk, jblk, xloc, yloc, xg, yg, ixloc, iyloc
      INTEGER  :: xp, yp, xblk, yblk, npxl, ipxl, ix, iy
      INTEGER  :: iloc, iloc_max(2)
      INTEGER  :: iproc, idest, isrc
      INTEGER  :: ylg, yug, ysp, ynp, nyp
      INTEGER  :: xlg, xug, xwp, xep, nxp
      REAL(r8) :: dlatp, dlonp
      LOGICAL  :: is_new
      INTEGER  :: nsend, nrecv, irecv
      INTEGER  :: smesg(5), rmesg(5)
      INTEGER  :: nave, nres
      
      INTEGER, allocatable :: nbasin_worker(:)
      TYPE(pointer_int32_1d), allocatable :: ulist_worker(:)

      INTEGER, allocatable :: ulist(:), iaddr(:)
      INTEGER, allocatable :: ulist2(:,:), xlist2(:,:), ylist2(:,:)
      INTEGER, allocatable :: sbuf(:), ipt2(:,:), rbuf(:)
      INTEGER, allocatable :: ulist_recv(:), xlist_recv(:), ylist_recv(:)
      INTEGER, allocatable :: npxl_blk(:,:)
      INTEGER, allocatable :: idest_all(:)
      LOGICAL, allocatable :: msk2(:,:), msk(:)
      INTEGER, allocatable :: xlist(:), ylist(:)
      TYPE(irregular_basin_type), allocatable :: lbasin (:)
      LOGICAL, allocatable :: work_done(:)
      INTEGER, allocatable :: blkdsp(:,:), blkcnt(:,:)
      INTEGER :: iblk_p, jblk_p
      INTEGER :: nbasin_glb

      INTEGER, allocatable :: basinindx(:), order(:)

#ifdef SinglePoint

      numbasin = 1
      allocate (landbasin(1))
      landbasin(1)%indx = 1
      
      landbasin(1)%npxl = 1

      allocate(landbasin(1)%ilat(1))
      landbasin(1)%ilat(1) = find_nearest_south (SITE_lat_location, pixel%nlat, pixel%lat_s)

      allocate(landbasin(1)%ilon(1))
      CALL normalize_longitude (SITE_lon_location)
      landbasin(1)%ilon(1) = find_nearest_west  (SITE_lon_location, pixel%nlon, pixel%lon_w)
      
      landbasin(1)%xblk = find_nearest_west  (pixel%lon_w(landbasin(1)%ilon(1)), gblock%nxblk, gblock%lon_w)
      landbasin(1)%yblk = find_nearest_south (pixel%lat_s(landbasin(1)%ilat(1)), gblock%nyblk, gblock%lat_s)

      allocate(nbasin_blk(gblock%nxblk,gblock%nyblk))
      nbasin_blk(:,:) = 0
      nbasin_blk(landbasin(1)%xblk,landbasin(1)%yblk) = 1

      xblkme(1) = landbasin(1)%xblk
      yblkme(1) = landbasin(1)%yblk
      
      RETURN
#endif

      IF (p_is_io) THEN 

         CALL allocate_block_data (gbasin, databasin)

#ifdef GRIDBASED 
         CALL ncio_read_block (DEF_file_landgrid, 'landmask', gbasin, databasin)
#endif
#ifdef CATCHMENT
         CALL hydro_data_read (DEF_dir_hydrodata, 'icat', gbasin, databasin, spv = -1)
#endif
      ENDIF

      ! Step 1: How many basins in each block?
      IF (p_is_io) THEN

         nbasin = 0

         allocate (nbasin_worker (0:p_np_worker-1))
         nbasin_worker(:) = 0

         allocate (ulist_worker (0:p_np_worker-1))
         DO iworker = 0, p_np_worker-1
            allocate (ulist_worker(iworker)%val (1000))
         ENDDO

         DO iblkme = 1, nblkme 
            iblk = xblkme(iblkme)
            jblk = yblkme(iblkme)
                  
            DO yloc = 1, gbasin%ycnt(jblk)
               DO xloc = 1, gbasin%xcnt(iblk)

#ifdef GRIDBASED 
                  IF (databasin%blk(iblk,jblk)%val(xloc,yloc) > 0) THEN
                     xg = gbasin%xdsp(iblk) + xloc
                     IF (xg > gbasin%nlon) xg = xg - gbasin%nlon

                     yg = gbasin%ydsp(jblk) + yloc

                     iu = gbasin%nlon * (yg-1) + xg
                  ELSE
                     iu = 0
                  ENDIF
#endif
#ifdef CATCHMENT
                  iu = databasin%blk(iblk,jblk)%val(xloc,yloc)
#endif

                  IF (iu > 0) THEN

                     iworker = mod(iu, p_np_worker)
                     CALL insert_into_sorted_list1 ( &
                        iu, nbasin_worker(iworker), ulist_worker(iworker)%val, iloc)

                     IF (nbasin_worker(iworker) == size(ulist_worker(iworker)%val)) THEN
                        CALL expand_list (ulist_worker(iworker)%val, 0.2_r8)
                     ENDIF

                  ENDIF

               ENDDO
            ENDDO
                        
#ifdef USEMPI
            DO iworker = 0, p_np_worker-1
               IF (nbasin_worker(iworker) > 0) then 
                  idest = p_address_worker(iworker)
                  smesg(1:2) = (/p_iam_glb, nbasin_worker(iworker)/) 
                  ! send(01)
                  CALL mpi_send (smesg(1:2), 2, MPI_INTEGER, &
                     idest, mpi_tag_size, p_comm_glb, p_err) 
               ENDIF
            ENDDO
#endif

            ! IF (sum(nbasin_worker) > 0) THEN
            !    write(*,*) 'Found ', sum(nbasin_worker), ' on block', iblk, jblk
            ! ENDIF

            nbasin = nbasin + sum(nbasin_worker)
            nbasin_worker(:) = 0
         ENDDO               
         
#ifdef USEMPI
         DO iworker = 0, p_np_worker-1
            idest = p_address_worker(iworker)
            ! send(02)
            smesg(1:2) = (/p_iam_glb, 0/)
            CALL mpi_send (smesg(1:2), 2, MPI_INTEGER, &
               idest, mpi_tag_size, p_comm_glb, p_err) 
         ENDDO
#endif

         deallocate (nbasin_worker)
         DO iworker = 0, p_np_worker-1
            deallocate (ulist_worker(iworker)%val)
         ENDDO
         deallocate (ulist_worker)

      ENDIF 

#ifdef USEMPI
      IF (p_is_worker) THEN
         nbasin = 0
         allocate(work_done(0:p_np_io-1))
         work_done(:) = .false.
         DO WHILE (.not. all(work_done))
            ! recv(01,02)
            CALL mpi_recv (rmesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_size, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            nrecv = rmesg(2)

            IF (nrecv > 0) THEN
               nbasin = nbasin + nrecv
            ELSE 
               work_done(p_itis_io(isrc)) = .true.
            ENDIF 
         ENDDO 

         deallocate(work_done)
      ENDIF
#endif

      ! Step 2: Build pixel list for each basin.
      IF (p_is_worker) THEN
         IF (nbasin > 0) THEN
            allocate (lbasin (nbasin))
            allocate (ulist (nbasin))
            allocate (iaddr (nbasin))
         ENDIF 
         nbasin = 0
      ENDIF 

      IF (p_is_io) THEN

         DO iblkme = 1, nblkme 
            iblk = xblkme(iblkme)
            jblk = yblkme(iblkme)
            IF (gbasin%xcnt(iblk) <= 0) cycle 
            IF (gbasin%ycnt(jblk) <= 0) cycle 

            ylg = gbasin%ydsp(jblk) + 1
            yug = gbasin%ydsp(jblk) + gbasin%ycnt(jblk)
            IF (gbasin%yinc == 1) THEN
               ysp = find_nearest_south (gbasin%lat_s(ylg), pixel%nlat, pixel%lat_s)
               ynp = find_nearest_north (gbasin%lat_n(yug), pixel%nlat, pixel%lat_n)
            ELSE
               ysp = find_nearest_south (gbasin%lat_s(yug), pixel%nlat, pixel%lat_s)
               ynp = find_nearest_north (gbasin%lat_n(ylg), pixel%nlat, pixel%lat_n)
            ENDIF

            nyp = ynp - ysp + 1

            xlg = gbasin%xdsp(iblk) + 1
            xug = gbasin%xdsp(iblk) + gbasin%xcnt(iblk)
            IF (xug > gbasin%nlon) xug = xug - gbasin%nlon

            xwp = find_nearest_west (gbasin%lon_w(xlg), pixel%nlon, pixel%lon_w)
            IF (.not. lon_between_floor(pixel%lon_w(xwp), gbasin%lon_w(xlg), gbasin%lon_e(xlg))) THEN
               xwp = mod(xwp,pixel%nlon) + 1
            ENDIF

            xep = find_nearest_east (gbasin%lon_e(xug), pixel%nlon, pixel%lon_e)
            IF (.not. lon_between_ceil(pixel%lon_e(xep), gbasin%lon_w(xug), gbasin%lon_e(xug))) THEN
               xep = xep - 1
               IF (xep == 0) xep = pixel%nlon
            ENDIF

            nxp = xep - xwp + 1
            IF (nxp <= 0) nxp = nxp + pixel%nlon 

            allocate (ulist2 (nxp,nyp))
            allocate (xlist2 (nxp,nyp))
            allocate (ylist2 (nxp,nyp))
            allocate (msk2   (nxp,nyp))

            DO iy = ysp, ynp
               yg = gbasin%ygrd(iy)
               yloc = gbasin%yloc(yg)

               iyloc = iy - ysp + 1
               dlatp = pixel%lat_n(iy) - pixel%lat_s(iy)
               IF (dlatp < 1.0e-7_r8) THEN
                  ulist2(:,iyloc) = 0
                  cycle
               ENDIF

               ix = xwp
               ixloc = 0
               DO while (.true.)
                  ixloc = ixloc + 1
                  dlonp = pixel%lon_e(ix) - pixel%lon_w(ix)
                  IF (dlonp < 0) dlonp = dlonp + 360.0_r8 

                  xg = gbasin%xgrd(ix)
                  xloc = gbasin%xloc(xg)

#ifdef GRIDBASED 
                  IF (databasin%blk(iblk,jblk)%val(xloc,yloc) > 0) THEN
                     iu = gbasin%nlon * (yg-1) + xg
                  ELSE
                     iu = 0
                  ENDIF
#endif
#ifdef CATCHMENT
                  iu = databasin%blk(iblk,jblk)%val(xloc,yloc)
#endif

                  xlist2(ixloc,iyloc) = ix
                  ylist2(ixloc,iyloc) = iy
                  ulist2(ixloc,iyloc) = iu

                  IF (dlonp < 1.0e-7_r8) THEN
                     ulist2(ixloc,iyloc) = 0
                  ENDIF

                  IF (ix == xep) exit
                  ix = mod(ix,pixel%nlon) + 1
               ENDDO
            ENDDO
                  
#ifdef USEMPI
            allocate (sbuf (nxp*nyp))
            allocate (ipt2 (nxp,nyp))
                  
            ipt2 = mod(ulist2, p_np_worker)
            DO iproc = 0, p_np_worker-1
               msk2  = (ipt2 == iproc) .and. (ulist2 > 0)
               nsend = count(msk2)
               IF (nsend > 0) THEN

                  idest = p_address_worker(iproc)

                  smesg(1:2) = (/p_iam_glb, nsend/)
                  ! send(03)
                  CALL mpi_send (smesg(1:2), 2, MPI_INTEGER, &
                     idest, mpi_tag_mesg, p_comm_glb, p_err) 

                  sbuf(1:nsend) = pack(ulist2, msk2)
                  ! send(04)
                  CALL mpi_send (sbuf(1:nsend), nsend, MPI_INTEGER, &
                     idest, mpi_tag_data, p_comm_glb, p_err) 

                  sbuf(1:nsend) = pack(xlist2, msk2)
                  ! send(05)
                  CALL mpi_send (sbuf(1:nsend), nsend, MPI_INTEGER, &
                     idest, mpi_tag_data, p_comm_glb, p_err) 

                  sbuf(1:nsend) = pack(ylist2, msk2)
                  ! send(06)
                  CALL mpi_send (sbuf(1:nsend), nsend, MPI_INTEGER, &
                     idest, mpi_tag_data, p_comm_glb, p_err) 

               ENDIF
            ENDDO

            deallocate (sbuf  )
            deallocate (ipt2  )
#else
                  
            DO iy = 1, nyp
               DO ix = 1, nxp

                  iu = ulist2(ix,iy)
                  IF (iu > 0) THEN 

                     CALL insert_into_sorted_list1 (iu, nbasin, ulist, iloc, is_new)

                     msk2 = (ulist2 == iu)
                     npxl = count(msk2)

                     IF (is_new) THEN
                        IF (iloc < nbasin) THEN
                           iaddr(iloc+1:nbasin) = iaddr(iloc:nbasin-1)
                        ENDIF
                        iaddr(iloc) = nbasin

                        lbasin(iaddr(iloc))%indx = iu
                        lbasin(iaddr(iloc))%npxl = npxl
                     ELSE
                        lbasin(iaddr(iloc))%npxl = lbasin(iaddr(iloc))%npxl + npxl
                     ENDIF

                     allocate (xlist(npxl))
                     allocate (ylist(npxl))
                     xlist = pack(xlist2, msk2)
                     ylist = pack(ylist2, msk2)

                     CALL append_to_list (lbasin(iaddr(iloc))%ilon, xlist)
                     CALL append_to_list (lbasin(iaddr(iloc))%ilat, ylist)

                     where(msk2) ulist2 = -1

                     deallocate (xlist)
                     deallocate (ylist)
                  ENDIF 

               ENDDO
            ENDDO
#endif

            deallocate (ulist2)
            deallocate (xlist2)
            deallocate (ylist2)
            deallocate (msk2  )

         ENDDO

#ifdef USEMPI
         DO iworker = 0, p_np_worker-1
            idest = p_address_worker(iworker)
            ! send(07)
            rmesg(1:2) = (/p_iam_glb, 0/)
            CALL mpi_send (rmesg(1:2), 2, MPI_INTEGER, &
               idest, mpi_tag_mesg, p_comm_glb, p_err) 
         ENDDO
#endif

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN

         allocate(work_done(0:p_np_io-1))
         work_done(:) = .false.
         DO WHILE (.not. all(work_done))
            ! recv(03,07)
            CALL mpi_recv (rmesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            nrecv = rmesg(2)
            IF (nrecv > 0) THEN

               allocate (ulist_recv (nrecv))
               ! recv(04)
               CALL mpi_recv (ulist_recv, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               
               allocate (xlist_recv (nrecv))
               ! recv(05)
               CALL mpi_recv (xlist_recv, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               
               allocate (ylist_recv (nrecv))
               ! recv(06)
               CALL mpi_recv (ylist_recv, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               allocate (msk(nrecv))

               DO irecv = 1, nrecv

                  iu = ulist_recv(irecv)

                  IF (iu > 0) THEN

                     CALL insert_into_sorted_list1 (iu, nbasin, ulist, iloc, is_new)

                     msk  = (ulist_recv == iu)
                     npxl = count(msk)

                     IF (is_new) THEN
                        IF (iloc < nbasin) THEN
                           iaddr(iloc+1:nbasin) = iaddr(iloc:nbasin-1)
                        ENDIF
                        iaddr(iloc) = nbasin

                        lbasin(iaddr(iloc))%indx = iu
                        lbasin(iaddr(iloc))%npxl = npxl
                     ELSE
                        lbasin(iaddr(iloc))%npxl = lbasin(iaddr(iloc))%npxl + npxl
                     ENDIF

                     allocate (xlist(npxl))
                     allocate (ylist(npxl))
                     xlist = pack(xlist_recv, msk)
                     ylist = pack(ylist_recv, msk)

                     CALL append_to_list (lbasin(iaddr(iloc))%ilon, xlist)
                     CALL append_to_list (lbasin(iaddr(iloc))%ilat, ylist)

                     where(msk) ulist_recv = -1
                     deallocate (xlist)
                     deallocate (ylist)
                  ENDIF 

               ENDDO 

               deallocate (msk)
               deallocate (ulist_recv)
               deallocate (xlist_recv)
               deallocate (ylist_recv)
            ELSE
               work_done(p_itis_io(isrc)) = .true.
            ENDIF
         ENDDO

      ENDIF 
#endif

      ! Step 3: Which block each basin locates at. 
      IF (p_is_worker) THEN 

         allocate (npxl_blk   (gblock%nxblk,gblock%nyblk))
         allocate (nbasin_blk (gblock%nxblk,gblock%nyblk))

         nbasin_blk(:,:) = 0

         DO iu = 1, nbasin

            npxl_blk (:,:) = 0
            
            DO ipxl = 1, lbasin(iu)%npxl
               xp = lbasin(iu)%ilon(ipxl)
               yp = lbasin(iu)%ilat(ipxl)

               xg = gbasin%xgrd(xp)
               yg = gbasin%ygrd(yp)

               xblk = gbasin%xblk(xg)
               yblk = gbasin%yblk(yg)

               npxl_blk(xblk,yblk) = npxl_blk(xblk,yblk) + 1
            ENDDO

            iloc_max = maxloc(npxl_blk)
            lbasin(iu)%xblk = iloc_max(1)
            lbasin(iu)%yblk = iloc_max(2)
            
            nbasin_blk(iloc_max(1), iloc_max(2)) = &
               nbasin_blk(iloc_max(1), iloc_max(2)) + 1 

         ENDDO

         deallocate (npxl_blk)
      ENDIF 

#ifdef USEMPI
      IF (.not. p_is_worker) THEN
         allocate (nbasin_blk (gblock%nxblk,gblock%nyblk))
         nbasin_blk(:,:) = 0
      ENDIF

      CALL mpi_allreduce (MPI_IN_PLACE, nbasin_blk, gblock%nxblk*gblock%nyblk, &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif 

      ! Step 4: IF MPI is used, sending basins from worker to their IO processes.
      
      IF (p_is_io) THEN 

         allocate (blkdsp (gblock%nxblk, gblock%nyblk))
         blkdsp(1,1) = 0
         DO jblk = 1, gblock%nyblk 
            DO iblk = 1, gblock%nxblk 
               IF ((iblk /= 1) .or. (jblk /= 1)) THEN
                  IF (iblk == 1) THEN
                     iblk_p = gblock%nxblk
                     jblk_p = jblk - 1
                  ELSE                       
                     iblk_p = iblk - 1
                     jblk_p = jblk
                  ENDIF 

                  IF (gblock%pio(iblk_p,jblk_p) == p_iam_glb) THEN
                     blkdsp(iblk,jblk) = blkdsp(iblk_p,jblk_p) + nbasin_blk(iblk_p,jblk_p) 
                  ELSE
                     blkdsp(iblk,jblk) = blkdsp(iblk_p,jblk_p) 
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN 

         allocate (idest_all(nbasin))
         DO iu = 1, nbasin
            idest_all(iu) = gblock%pio (lbasin(iu)%xblk, lbasin(iu)%yblk)
         ENDDO 
         
         DO iu = 1, nbasin
            idest = idest_all(iu)
            smesg(1) = p_iam_glb
            smesg(2:3) = (/lbasin(iu)%indx, lbasin(iu)%npxl/)
            smesg(4:5) = (/lbasin(iu)%xblk, lbasin(iu)%yblk/)
            ! send(09)
            CALL mpi_send (smesg(1:5), 5, MPI_INTEGER, &
               idest, mpi_tag_mesg, p_comm_glb, p_err) 
            ! send(10)
            CALL mpi_send (lbasin(iu)%ilon, lbasin(iu)%npxl, MPI_INTEGER, &
               idest, mpi_tag_data, p_comm_glb, p_err)
            ! send(11)
            CALL mpi_send (lbasin(iu)%ilat, lbasin(iu)%npxl, MPI_INTEGER, &
               idest, mpi_tag_data, p_comm_glb, p_err)
         ENDDO
            
         deallocate (idest_all)

      ENDIF 

      IF (p_is_io) THEN 

         numbasin = sum(nbasin_blk, mask = gblock%pio == p_iam_glb)

         IF (numbasin > 0) THEN

            allocate (landbasin (numbasin))

            allocate (blkcnt (gblock%nxblk, gblock%nyblk))
            blkcnt(:,:) = 0
            DO iu = 1, numbasin

               ! recv(09)
               CALL mpi_recv (rmesg, 5, MPI_INTEGER, &
                  MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

               isrc = rmesg(1)
               xblk = rmesg(4)
               yblk = rmesg(5)

               blkcnt(xblk,yblk) = blkcnt(xblk,yblk) + 1
               ju = blkdsp(xblk,yblk) + blkcnt(xblk,yblk)

               landbasin(ju)%indx = rmesg(2)
               landbasin(ju)%npxl = rmesg(3)
               landbasin(ju)%xblk = rmesg(4)
               landbasin(ju)%yblk = rmesg(5)

               allocate (landbasin(ju)%ilon (landbasin(ju)%npxl))
               allocate (landbasin(ju)%ilat (landbasin(ju)%npxl))

               ! recv(10)
               CALL mpi_recv (landbasin(ju)%ilon, landbasin(ju)%npxl, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               ! recv(11)
               CALL mpi_recv (landbasin(ju)%ilat, landbasin(ju)%npxl, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

            ENDDO

         ENDIF
      ENDIF 

#else
      numbasin = nbasin
      IF (numbasin > 0) THEN

         allocate (landbasin (numbasin))

         allocate (blkcnt (gblock%nxblk, gblock%nyblk))
         blkcnt(:,:) = 0
         DO iu = 1, numbasin

            xblk = lbasin(iu)%xblk
            yblk = lbasin(iu)%yblk

            blkcnt(xblk,yblk) = blkcnt(xblk,yblk) + 1
            ju = blkdsp(xblk,yblk) + blkcnt(xblk,yblk)

            CALL copy_basin (lbasin(iu), landbasin(ju))

         ENDDO

      ENDIF
#endif

      ! Step 4-2: sort basins.
      IF (p_is_io) THEN
         IF (allocated (lbasin)) THEN
            DO iu = 1, size(lbasin)
               IF (allocated(lbasin(iu)%ilon))  deallocate (lbasin(iu)%ilon)
               IF (allocated(lbasin(iu)%ilon))  deallocate (lbasin(iu)%ilat)
            ENDDO
            deallocate (lbasin)
         ENDIF

         IF (numbasin > 0) THEN
            allocate (lbasin (numbasin))
            DO iu = 1, numbasin
               CALL copy_basin(landbasin(iu), lbasin(iu))
            ENDDO

            DO iblkme = 1, nblkme 
               iblk = xblkme(iblkme)
               jblk = yblkme(iblkme)

               IF (blkcnt(iblk,jblk) > 0) THEN
                  allocate (basinindx (blkcnt(iblk,jblk)))
                  allocate (order     (blkcnt(iblk,jblk)))

                  DO iu = blkdsp(iblk,jblk)+1, blkdsp(iblk,jblk)+blkcnt(iblk,jblk)
                     basinindx(iu-blkdsp(iblk,jblk)) = landbasin(iu)%indx
                  ENDDO

                  order = (/ (iu, iu = 1, blkcnt(iblk,jblk)) /)
                  CALL quicksort (blkcnt(iblk,jblk), basinindx, order)

                  DO iu = 1, blkcnt(iblk,jblk)
                     CALL copy_basin (lbasin(blkdsp(iblk,jblk)+order(iu)), &
                        landbasin(blkdsp(iblk,jblk)+iu))
                  ENDDO

                  deallocate (basinindx)
                  deallocate (order    )
               ENDIF

            ENDDO
         ENDIF
      ENDIF

      IF (allocated(blkdsp)) deallocate(blkdsp)
      IF (allocated(blkcnt)) deallocate(blkcnt)
         
      IF (allocated (lbasin)) THEN
         DO iu = 1, size(lbasin)
            IF (allocated(lbasin(iu)%ilon))  deallocate (lbasin(iu)%ilon)
            IF (allocated(lbasin(iu)%ilon))  deallocate (lbasin(iu)%ilat)
         ENDDO

         deallocate (lbasin )
      ENDIF

      ! Step 5: IF MPI is used, scatter basins from IO to workers.
#ifdef USEMPI
      CALL scatter_landbasin_from_io_to_worker
#endif
         
      IF (allocated(ulist)) deallocate (ulist)
      IF (allocated(iaddr)) deallocate (iaddr)

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land basins :'
      ENDIF
      
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_io) THEN
         CALL mpi_reduce (numbasin, nbasin_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_io, p_err)
         IF (p_iam_io == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', nbasin_glb, ' basins on IO.'
         ENDIF
      ENDIF
      IF (p_is_worker) THEN
         CALL mpi_reduce (numbasin, nbasin_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', nbasin_glb, ' basins on worker.'
         ENDIF
      ENDIF
      
      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numbasin, ' basins.'
#endif

   END SUBROUTINE landbasin_build


#ifdef USEMPI
   ! --------------------------------
   SUBROUTINE scatter_landbasin_from_io_to_worker

      USE spmd_task
      USE mod_block
      IMPLICIT NONE

      ! Local variables
      INTEGER :: iblk, jblk, nave, nres, iproc, ndsp, nsend, idest, isrc, iu
      INTEGER :: smesg(4), rmesg(4)
      INTEGER, allocatable :: nbasin_worker(:)
      INTEGER :: iblkme
      
      IF (p_is_io) THEN

         allocate (nbasin_worker (1:p_np_group-1))
         nbasin_worker(:) = 0

         DO iblkme = 1, nblkme 
            iblk = xblkme(iblkme)
            jblk = yblkme(iblkme)
        
            nave = nbasin_blk(iblk,jblk) / (p_np_group-1)
            nres = mod(nbasin_blk(iblk,jblk), p_np_group-1)
            DO iproc = 1, p_np_group-1
               nbasin_worker(iproc) = nbasin_worker(iproc) + nave
               IF (iproc <= nres)  nbasin_worker(iproc) = nbasin_worker(iproc) + 1
            ENDDO
         ENDDO

         DO iproc = 1, p_np_group-1
            CALL mpi_send (nbasin_worker(iproc), 1, MPI_INTEGER, &
               iproc, mpi_tag_size, p_comm_group, p_err) 
         ENDDO
         deallocate (nbasin_worker)

         ndsp = 0
         DO iblkme = 1, nblkme 
            iblk = xblkme(iblkme)
            jblk = yblkme(iblkme)

            nave = nbasin_blk(iblk,jblk) / (p_np_group-1)
            nres = mod(nbasin_blk(iblk,jblk), p_np_group-1)
            DO iproc = 1, p_np_group-1
               nsend = nave
               IF (iproc <= nres)  nsend = nsend + 1

               DO iu = ndsp+1, ndsp+nsend
                  idest = iproc
                  smesg(1:2) = (/landbasin(iu)%indx, landbasin(iu)%npxl/)
                  smesg(3:4) = (/landbasin(iu)%xblk, landbasin(iu)%yblk/)
                  CALL mpi_send (smesg(1:4), 4, MPI_INTEGER, &
                     idest, mpi_tag_mesg, p_comm_group, p_err) 
                  CALL mpi_send (landbasin(iu)%ilon, landbasin(iu)%npxl, &
                     MPI_INTEGER, idest, mpi_tag_data, p_comm_group, p_err)
                  CALL mpi_send (landbasin(iu)%ilat, landbasin(iu)%npxl, &
                     MPI_INTEGER, idest, mpi_tag_data, p_comm_group, p_err)
               ENDDO
               ndsp = ndsp + nsend
            ENDDO
         ENDDO
            
      ENDIF

      IF (p_is_worker) THEN

         CALL mpi_recv (numbasin, 1, MPI_INTEGER, &
            p_root, mpi_tag_size, p_comm_group, p_stat, p_err)

         IF (numbasin > 0) THEN
            allocate (landbasin (numbasin))

            DO iu = 1, numbasin
               CALL mpi_recv (rmesg, 4, MPI_INTEGER, &
                  p_root, mpi_tag_mesg, p_comm_group, p_stat, p_err)

               landbasin(iu)%indx = rmesg(1)
               landbasin(iu)%npxl = rmesg(2)
               landbasin(iu)%xblk = rmesg(3)
               landbasin(iu)%yblk = rmesg(4)

               allocate (landbasin(iu)%ilon (landbasin(iu)%npxl))
               allocate (landbasin(iu)%ilat (landbasin(iu)%npxl))

               CALL mpi_recv (landbasin(iu)%ilon, landbasin(iu)%npxl, MPI_INTEGER, &
                  p_root, mpi_tag_data, p_comm_group, p_stat, p_err)
               CALL mpi_recv (landbasin(iu)%ilat, landbasin(iu)%npxl, MPI_INTEGER, &
                  p_root, mpi_tag_data, p_comm_group, p_stat, p_err)
            ENDDO
         ENDIF

      ENDIF

   END SUBROUTINE scatter_landbasin_from_io_to_worker

#endif

   ! --------------------------------
   SUBROUTINE landbasin_free_mem ()
      
      IMPLICIT NONE

      ! Local variables
      INTEGER :: iu

      IF (allocated(landbasin)) THEN
         DO iu = 1, numbasin
            deallocate (landbasin(iu)%ilon)
            deallocate (landbasin(iu)%ilat)
         ENDDO

         deallocate (landbasin)
      ENDIF

   END SUBROUTINE landbasin_free_mem
   

END MODULE mod_landbasin
