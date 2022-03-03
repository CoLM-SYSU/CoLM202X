#include <define.h>

MODULE mod_landunit

   USE precision
   IMPLICIT NONE

   ! ---- data types ----
   TYPE :: irregular_unit_type 
      
      INTEGER :: num
      INTEGER :: xblk, yblk

      INTEGER :: npxl
      INTEGER, allocatable :: ilon(:)
      INTEGER, allocatable :: ilat(:)

   END TYPE irregular_unit_type

   ! ---- Instance ----
   INTEGER :: numunit
   TYPE (irregular_unit_type), allocatable :: landunit (:)

   INTEGER, allocatable :: nunits_blk(:,:)
   
CONTAINS
   
   ! --------------------------------
   SUBROUTINE landunit_build (gunit)

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

      TYPE(grid_type), intent(in) :: gunit

      ! Local Variables 
      CHARACTER(len=256) :: file_lunit
      CHARACTER(len=256) :: varname
      TYPE(block_data_int32_2d) :: dataunit

      INTEGER  :: iworker
      INTEGER  :: nunit, iu, ju
      INTEGER  :: iblk, jblk, xloc, yloc, xg, yg, ixloc, iyloc
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
      
      INTEGER, allocatable :: nunit_worker(:)
      TYPE(pointer_int32_1d), allocatable :: ulist_worker(:)

      INTEGER, allocatable :: ulist(:), iaddr(:)
      INTEGER, allocatable :: ulist2(:,:), xlist2(:,:), ylist2(:,:)
      INTEGER, allocatable :: sbuf(:), ipt2(:,:), rbuf(:)
      INTEGER, allocatable :: ulist_recv(:), xlist_recv(:), ylist_recv(:)
      INTEGER, allocatable :: npxl_blk(:,:)
      INTEGER, allocatable :: idest_all(:)
      LOGICAL, allocatable :: msk2(:,:), msk(:)
      INTEGER, allocatable :: xlist(:), ylist(:)
      TYPE(irregular_unit_type), allocatable :: lunit (:)
      LOGICAL, allocatable :: work_done(:)
      INTEGER, allocatable :: blkdsp(:,:), blkcnt(:,:)
      INTEGER :: iblk_p, jblk_p

      IF (p_is_io) THEN 

         CALL allocate_block_data (gunit, dataunit)

#ifdef GRIDBASED 
         CALL ncio_read_block (DEF_file_landgrid, 'landmask', gunit, dataunit)
#endif
#ifdef CATCHMENT
         CALL hydro_data_read (DEF_dir_hydrodata, 'icat', gunit, dataunit, spv = -1)
#endif
      ENDIF

      ! Step 1: How many units in each block?
      IF (p_is_io) THEN

         nunit = 0

         allocate (nunit_worker (0:p_np_worker-1))
         nunit_worker(:) = 0

         allocate (ulist_worker (0:p_np_worker-1))
         DO iworker = 0, p_np_worker-1
            allocate (ulist_worker(iworker)%val (1000))
         ENDDO

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  DO yloc = 1, gunit%ycnt(jblk)
                     DO xloc = 1, gunit%xcnt(iblk)

#ifdef GRIDBASED 
                        IF (dataunit%blk(iblk,jblk)%val(xloc,yloc) > 0) THEN
                           xg = gunit%xdsp(iblk) + xloc
                           IF (xg > gunit%nlon) xg = xg - gunit%nlon

                           yg = gunit%ydsp(jblk) + yloc

                           iu = gunit%nlon * (yg-1) + xg
                        ELSE
                           iu = 0
                        ENDIF
#endif
#ifdef CATCHMENT
                        iu = dataunit%blk(iblk,jblk)%val(xloc,yloc)
#endif

                        IF (iu > 0) THEN

                           iworker = mod(iu, p_np_worker)
                           CALL insert_into_sorted_list1 ( &
                              iu, nunit_worker(iworker), ulist_worker(iworker)%val, iloc)
                        
                           IF (nunit_worker(iworker) == size(ulist_worker(iworker)%val)) THEN
                              CALL expand_list (ulist_worker(iworker)%val, 0.2_r8)
                           ENDIF
                           
                        ENDIF

                     ENDDO
                  ENDDO
                        
#ifdef USEMPI
                  DO iworker = 0, p_np_worker-1
                     IF (nunit_worker(iworker) > 0) then 
                        idest = p_address_worker(iworker)
                        smesg(1:2) = (/p_iam_glb, nunit_worker(iworker)/) 
                        ! send(01)
                        CALL mpi_send (smesg(1:2), 2, MPI_INTEGER, &
                           idest, mpi_tag_size, p_comm_glb, p_err) 
                     ENDIF
                  ENDDO
#endif

                  IF (sum(nunit_worker) > 0) THEN
                     write(*,*) 'Found ', sum(nunit_worker), ' on block', iblk, jblk
                  ENDIF

                  nunit = nunit + sum(nunit_worker)
                  nunit_worker(:) = 0
               ENDIF
            ENDDO
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

         deallocate (nunit_worker)
         DO iworker = 0, p_np_worker-1
            deallocate (ulist_worker(iworker)%val)
         ENDDO
         deallocate (ulist_worker)

      ENDIF 

#ifdef USEMPI
      IF (p_is_worker) THEN
         nunit = 0
         allocate(work_done(0:p_np_io-1))
         work_done(:) = .false.
         DO WHILE (.not. all(work_done))
            ! recv(01,02)
            CALL mpi_recv (rmesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_size, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            nrecv = rmesg(2)

            IF (nrecv > 0) THEN
               nunit = nunit + nrecv
            ELSE 
               work_done(p_itis_io(isrc)) = .true.
            ENDIF 
         ENDDO 

         deallocate(work_done)
      ENDIF
#endif

      ! Step 2: Build pixel list for each unit.
      IF (p_is_worker) THEN
         IF (nunit > 0) THEN
            allocate (lunit (nunit))
            allocate (ulist (nunit))
            allocate (iaddr (nunit))
         ENDIF 
         nunit = 0
      ENDIF 

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            IF (gunit%ycnt(jblk) <= 0) cycle 

            DO iblk = 1, gblock%nxblk
               IF (gunit%xcnt(iblk) <= 0) cycle 

               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
         
                  ylg = gunit%ydsp(jblk) + 1
                  yug = gunit%ydsp(jblk) + gunit%ycnt(jblk)
                  IF (gunit%yinc == 1) THEN
                     ysp = find_nearest_south (gunit%lat_s(ylg), pixel%nlat, pixel%lat_s)
                     ynp = find_nearest_north (gunit%lat_n(yug), pixel%nlat, pixel%lat_n)
                  ELSE
                     ysp = find_nearest_south (gunit%lat_s(yug), pixel%nlat, pixel%lat_s)
                     ynp = find_nearest_north (gunit%lat_n(ylg), pixel%nlat, pixel%lat_n)
                  ENDIF

                  nyp = ynp - ysp + 1

                  xlg = gunit%xdsp(iblk) + 1
                  xug = gunit%xdsp(iblk) + gunit%xcnt(iblk)
                  IF (xug > gunit%nlon) xug = xug - gunit%nlon
                  
                  xwp = find_nearest_west (gunit%lon_w(xlg), pixel%nlon, pixel%lon_w)
                  IF (.not. lon_between_floor(pixel%lon_w(xwp), gunit%lon_w(xlg), gunit%lon_e(xlg))) THEN
                     xwp = mod(xwp,pixel%nlon) + 1
                  ENDIF
                  
                  xep = find_nearest_east (gunit%lon_e(xug), pixel%nlon, pixel%lon_e)
                  IF (.not. lon_between_ceil(pixel%lon_e(xep), gunit%lon_w(xug), gunit%lon_e(xug))) THEN
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
                     yg = gunit%ygrd(iy)
                     yloc = gunit%yloc(yg)
                     
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

                        xg = gunit%xgrd(ix)
                        xloc = gunit%xloc(xg)

#ifdef GRIDBASED 
                        IF (dataunit%blk(iblk,jblk)%val(xloc,yloc) > 0) THEN
                           iu = gunit%nlon * (yg-1) + xg
                        ELSE
                           iu = 0
                        ENDIF
#endif
#ifdef CATCHMENT
                        iu = dataunit%blk(iblk,jblk)%val(xloc,yloc)
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
                           
                           CALL insert_into_sorted_list1 (iu, nunit, ulist, iloc, is_new)
                           
                           msk2 = (ulist2 == iu)
                           npxl = count(msk2)

                           IF (is_new) THEN
                              IF (iloc < nunit) THEN
                                 iaddr(iloc+1:nunit) = iaddr(iloc:nunit-1)
                              ENDIF
                              iaddr(iloc) = nunit

                              lunit(iaddr(iloc))%num  = iu
                              lunit(iaddr(iloc))%npxl = npxl
                           ELSE
                              lunit(iaddr(iloc))%npxl = lunit(iaddr(iloc))%npxl + npxl
                           ENDIF
                           
                           allocate (xlist(npxl))
                           allocate (ylist(npxl))
                           xlist = pack(xlist2, msk2)
                           ylist = pack(ylist2, msk2)
                        
                           CALL append_to_list (lunit(iaddr(iloc))%ilon, xlist)
                           CALL append_to_list (lunit(iaddr(iloc))%ilat, ylist)

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

               ENDIF
            ENDDO
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

                     CALL insert_into_sorted_list1 (iu, nunit, ulist, iloc, is_new)

                     msk  = (ulist_recv == iu)
                     npxl = count(msk)

                     IF (is_new) THEN
                        IF (iloc < nunit) THEN
                           iaddr(iloc+1:nunit) = iaddr(iloc:nunit-1)
                        ENDIF
                        iaddr(iloc) = nunit

                        lunit(iaddr(iloc))%num  = iu
                        lunit(iaddr(iloc))%npxl = npxl
                     ELSE
                        lunit(iaddr(iloc))%npxl = lunit(iaddr(iloc))%npxl + npxl
                     ENDIF

                     allocate (xlist(npxl))
                     allocate (ylist(npxl))
                     xlist = pack(xlist_recv, msk)
                     ylist = pack(ylist_recv, msk)

                     CALL append_to_list (lunit(iaddr(iloc))%ilon, xlist)
                     CALL append_to_list (lunit(iaddr(iloc))%ilat, ylist)

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

      ! Step 3: Which block each unit locates at. 
      IF (p_is_worker) THEN 

         allocate (npxl_blk   (gblock%nxblk,gblock%nyblk))
         allocate (nunits_blk (gblock%nxblk,gblock%nyblk))

         nunits_blk(:,:) = 0

         DO iu = 1, nunit

            npxl_blk (:,:) = 0
            
            DO ipxl = 1, lunit(iu)%npxl
               xp = lunit(iu)%ilon(ipxl)
               yp = lunit(iu)%ilat(ipxl)

               xg = gunit%xgrd(xp)
               yg = gunit%ygrd(yp)

               xblk = gunit%xblk(xg)
               yblk = gunit%yblk(yg)

               npxl_blk(xblk,yblk) = npxl_blk(xblk,yblk) + 1
            ENDDO

            iloc_max = maxloc(npxl_blk)
            lunit(iu)%xblk = iloc_max(1)
            lunit(iu)%yblk = iloc_max(2)
            
            nunits_blk(iloc_max(1), iloc_max(2)) = &
               nunits_blk(iloc_max(1), iloc_max(2)) + 1 

         ENDDO

         deallocate (npxl_blk)
      ENDIF 

#ifdef USEMPI
      IF (.not. p_is_worker) THEN
         allocate (nunits_blk (gblock%nxblk,gblock%nyblk))
         nunits_blk(:,:) = 0
      ENDIF

      CALL mpi_allreduce (MPI_IN_PLACE, nunits_blk, gblock%nxblk*gblock%nyblk, &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif 

      ! Step 4: IF MPI is used, sending units from worker to their IO processes.
      
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
                     blkdsp(iblk,jblk) = blkdsp(iblk_p,jblk_p) + nunits_blk(iblk_p,jblk_p) 
                  ELSE
                     blkdsp(iblk,jblk) = blkdsp(iblk_p,jblk_p) 
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN 

         allocate (idest_all(nunit))
         DO iu = 1, nunit
            idest_all(iu) = gblock%pio (lunit(iu)%xblk, lunit(iu)%yblk)
         ENDDO 
         
         DO iu = 1, nunit
            idest = idest_all(iu)
            smesg(1) = p_iam_glb
            smesg(2:3) = (/lunit(iu)%num,  lunit(iu)%npxl/)
            smesg(4:5) = (/lunit(iu)%xblk, lunit(iu)%yblk/)
            ! send(09)
            CALL mpi_send (smesg(1:5), 5, MPI_INTEGER, &
               idest, mpi_tag_mesg, p_comm_glb, p_err) 
            ! send(10)
            CALL mpi_send (lunit(iu)%ilon, lunit(iu)%npxl, MPI_INTEGER, &
               idest, mpi_tag_data, p_comm_glb, p_err)
            ! send(11)
            CALL mpi_send (lunit(iu)%ilat, lunit(iu)%npxl, MPI_INTEGER, &
               idest, mpi_tag_data, p_comm_glb, p_err)
         ENDDO
            
         deallocate (idest_all)

      ENDIF 

      IF (p_is_io) THEN 

         numunit = sum(nunits_blk, mask = gblock%pio == p_iam_glb)

         IF (numunit > 0) THEN

            allocate (landunit (numunit))

            allocate (blkcnt (gblock%nxblk, gblock%nyblk))
            blkcnt(:,:) = 0
            DO iu = 1, numunit

               ! recv(09)
               CALL mpi_recv (rmesg, 5, MPI_INTEGER, &
                  MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

               isrc = rmesg(1)
               xblk = rmesg(4)
               yblk = rmesg(5)

               blkcnt(xblk,yblk) = blkcnt(xblk,yblk) + 1
               ju = blkdsp(xblk,yblk) + blkcnt(xblk,yblk)

               landunit(ju)%num  = rmesg(2)
               landunit(ju)%npxl = rmesg(3)
               landunit(ju)%xblk = rmesg(4)
               landunit(ju)%yblk = rmesg(5)

               allocate (landunit(ju)%ilon (landunit(ju)%npxl))
               allocate (landunit(ju)%ilat (landunit(ju)%npxl))

               ! recv(10)
               CALL mpi_recv (landunit(ju)%ilon, landunit(ju)%npxl, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               ! recv(11)
               CALL mpi_recv (landunit(ju)%ilat, landunit(ju)%npxl, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

            ENDDO

         ENDIF
      ENDIF 

#else
      numunit = nunit
      IF (numunit > 0) THEN

         allocate (landunit (numunit))

         allocate (blkcnt (gblock%nxblk, gblock%nyblk))
         blkcnt(:,:) = 0
         DO iu = 1, numunit

            blkcnt(xblk,yblk) = blkcnt(xblk,yblk) + 1
            ju = blkdsp(xblk,yblk) + blkcnt(xblk,yblk)

            landunit(ju)%num  = lunit(iu)%num 
            landunit(ju)%npxl = lunit(iu)%npxl
            landunit(ju)%xblk = lunit(iu)%xblk
            landunit(ju)%yblk = lunit(iu)%yblk

            allocate (landunit(ju)%ilon (landunit(ju)%npxl))
            allocate (landunit(ju)%ilat (landunit(ju)%npxl))

            landunit(ju)%ilon = lunit(iu)%ilon
            landunit(ju)%ilat = lunit(iu)%ilat
         ENDDO

      ENDIF
#endif

      IF (allocated(blkdsp)) deallocate(blkdsp)
      IF (allocated(blkcnt)) deallocate(blkcnt)
         
      IF (allocated (lunit)) THEN
         DO iu = 1, nunit
            deallocate (lunit(iu)%ilon)
            deallocate (lunit(iu)%ilat)
         ENDDO

         deallocate (lunit )
      ENDIF

      ! Step 5: IF MPI is used, scatter units from IO to workers.
#ifdef USEMPI
      CALL scatter_landunit_from_io_to_worker
#endif
         
      IF (allocated(ulist)) deallocate (ulist)
      IF (allocated(iaddr)) deallocate (iaddr)

      IF (p_is_master) THEN
         write(*,*) 'Making land units :'
      ENDIF
      
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN
         write(*,'(I10,A12,I5)') numunit, ' units on IO', p_iam_glb
      ENDIF
      
      IF (p_is_worker) THEN
         write(*,'(I10,A16,I5)') numunit, ' units on worker', p_iam_glb
      ENDIF

   END SUBROUTINE landunit_build


#ifdef USEMPI
   ! --------------------------------
   SUBROUTINE scatter_landunit_from_io_to_worker

      USE spmd_task
      USE mod_block
      IMPLICIT NONE

      ! Local variables
      INTEGER :: iblk, jblk, nave, nres, iproc, ndsp, nsend, idest, isrc, iu
      INTEGER :: smesg(4), rmesg(4)
      INTEGER, allocatable :: nunit_worker(:)
      
      IF (p_is_io) THEN

         allocate (nunit_worker (1:p_np_group-1))
         nunit_worker(:) = 0

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
        
                  nave = nunits_blk(iblk,jblk) / (p_np_group-1)
                  nres = mod(nunits_blk(iblk,jblk), p_np_group-1)
                  DO iproc = 1, p_np_group-1
                     nunit_worker(iproc) = nunit_worker(iproc) + nave
                     IF (iproc <= nres)  nunit_worker(iproc) = nunit_worker(iproc) + 1
                  ENDDO

               ENDIF
            ENDDO
         ENDDO

         DO iproc = 1, p_np_group-1
            CALL mpi_send (nunit_worker(iproc), 1, MPI_INTEGER, &
               iproc, mpi_tag_size, p_comm_group, p_err) 
         ENDDO
         deallocate (nunit_worker)

         ndsp = 0
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                  nave = nunits_blk(iblk,jblk) / (p_np_group-1)
                  nres = mod(nunits_blk(iblk,jblk), p_np_group-1)
                  DO iproc = 1, p_np_group-1
                     nsend = nave
                     IF (iproc <= nres)  nsend = nsend + 1

                     DO iu = ndsp+1, ndsp+nsend
                        idest = iproc
                        smesg(1:2) = (/landunit(iu)%num,  landunit(iu)%npxl/)
                        smesg(3:4) = (/landunit(iu)%xblk, landunit(iu)%yblk/)
                        CALL mpi_send (smesg(1:4), 4, MPI_INTEGER, &
                           idest, mpi_tag_mesg, p_comm_group, p_err) 
                        CALL mpi_send (landunit(iu)%ilon, landunit(iu)%npxl, &
                           MPI_INTEGER, idest, mpi_tag_data, p_comm_group, p_err)
                        CALL mpi_send (landunit(iu)%ilat, landunit(iu)%npxl, &
                           MPI_INTEGER, idest, mpi_tag_data, p_comm_group, p_err)
                     ENDDO
                     ndsp = ndsp + nsend
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
            
      ENDIF

      IF (p_is_worker) THEN

         CALL mpi_recv (numunit, 1, MPI_INTEGER, &
            p_root, mpi_tag_size, p_comm_group, p_stat, p_err)

         IF (numunit > 0) THEN
            allocate (landunit (numunit))

            DO iu = 1, numunit
               CALL mpi_recv (rmesg, 4, MPI_INTEGER, &
                  p_root, mpi_tag_mesg, p_comm_group, p_stat, p_err)

               landunit(iu)%num  = rmesg(1)
               landunit(iu)%npxl = rmesg(2)
               landunit(iu)%xblk = rmesg(3)
               landunit(iu)%yblk = rmesg(4)

               allocate (landunit(iu)%ilon (landunit(iu)%npxl))
               allocate (landunit(iu)%ilat (landunit(iu)%npxl))

               CALL mpi_recv (landunit(iu)%ilon, landunit(iu)%npxl, MPI_INTEGER, &
                  p_root, mpi_tag_data, p_comm_group, p_stat, p_err)
               CALL mpi_recv (landunit(iu)%ilat, landunit(iu)%npxl, MPI_INTEGER, &
                  p_root, mpi_tag_data, p_comm_group, p_stat, p_err)
            ENDDO
         ENDIF

      ENDIF

   END SUBROUTINE scatter_landunit_from_io_to_worker

#endif

   ! --------------------------------
   SUBROUTINE landunit_free_mem ()
      
      IMPLICIT NONE

      ! Local variables
      INTEGER :: iu

      IF (allocated(landunit)) THEN
         DO iu = 1, numunit
            deallocate (landunit(iu)%ilon)
            deallocate (landunit(iu)%ilat)
         ENDDO

         deallocate (landunit)
      ENDIF

   END SUBROUTINE landunit_free_mem
   

END MODULE mod_landunit
