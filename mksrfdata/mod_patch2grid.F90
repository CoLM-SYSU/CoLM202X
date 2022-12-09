#include <define.h>

#ifdef MAP_PATCH_TO_GRID
MODULE mod_patch2grid

   USE precision
   USE spmd_task
   USE mod_block
   USE mod_data_type
   USE mod_namelist
   USE mod_grid
   USE mod_mapping_pset2grid
   USE GlobalVars, only : spval
   USE ncio_serial
   USE mod_landpatch
   
   INTEGER :: npatchtype
   INTEGER, allocatable :: patchindex(:) 
   
   TYPE(grid_type)              :: grid_patch2grid 
   TYPE(grid_concat_type)       :: concat_patch2grid
   TYPE(mapping_pset2grid_type) :: map_patch2grid
   TYPE(block_data_real8_3d)    :: wt_patch2grid

   CHARACTER(len=256) :: dir_patch2grid

CONTAINS 

   ! ----------
   SUBROUTINE patch2grid_init ()

      IMPLICIT NONE
      ! Local Variables
      INTEGER :: ipatchtype
      logical,  allocatable :: filter (:)
      REAL(r8), allocatable :: vectmp (:)
      TYPE(block_data_real8_2d) :: sumwt 
      INTEGER  :: iblkme, xblk, yblk, xloc, yloc 

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      
#if (defined IGBP_CLASSIFICATION || defined PFT_CLASSIFICATION || defined PC_CLASSIFICATION) 
      npatchtype = 18
#elif (defined USGS_CLASSIFICATION)
      npatchtype = 25
#endif

      allocate (patchindex(npatchtype))
      patchindex = (/ (ipatchtype, ipatchtype=0,npatchtype-1) /)

      call concat_patch2grid%set (grid_patch2grid)

#ifndef CROP
      call map_patch2grid%build (landpatch, grid_patch2grid)
#else
      call map_patch2grid%build (landpatch, grid_patch2grid, pctcrop)
#endif

      IF (p_is_master) THEN
         write(*,*) '  Map and Wrtie data on patches: patch fractions' 
      ENDIF

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (filter(numpatch))
            allocate (vectmp(numpatch))
         ENDIF
      ENDIF

      IF (p_is_io) THEN
         CALL allocate_block_data (grid_patch2grid, sumwt)
         call allocate_block_data (grid_patch2grid, wt_patch2grid, npatchtype, 0)  
      ENDIF
      
      DO ipatchtype = 0, npatchtype-1

         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = (landpatch%ltyp == ipatchtype)
               vectmp(:) = 1.
            end if
         end if

         call map_patch2grid%map (vectmp, sumwt, msk = filter)

         if (p_is_io) then
            DO iblkme = 1, gblock%nblkme 
               xblk = gblock%xblkme(iblkme)
               yblk = gblock%yblkme(iblkme)

               do yloc = 1, grid_patch2grid%ycnt(yblk) 
                  do xloc = 1, grid_patch2grid%xcnt(xblk) 

                     if (sumwt%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                        wt_patch2grid%blk(xblk,yblk)%val(ipatchtype,xloc,yloc) = &
                           sumwt%blk(xblk,yblk)%val(xloc,yloc)
                     else
                        wt_patch2grid%blk(xblk,yblk)%val(ipatchtype,xloc,yloc) = spval
                     end if

                  end do
               end do

            end do
         end if
      ENDDO 

      dir_patch2grid = trim(adjustl(DEF_dir_output)) &
         // '/' // trim(adjustl(DEF_CASE_NAME)) // '/landdata_gridded'

      IF (p_is_master) THEN
         CALL system('mkdir -p ' // trim(dir_patch2grid))
      ENDIF

      CALL write3d_patch2grid (wt_patch2grid, 'patchfrac', 'patchfrac.nc')

      IF (allocated(filter)) deallocate (filter)
      IF (allocated(vectmp)) deallocate (vectmp)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE patch2grid_init
   
   ! ----------
   SUBROUTINE map_patchdata_to_grid_and_write (patchdata, dataname, filename)

      USE spmd_task
      IMPLICIT NONE

      REAL(r8), intent(in) :: patchdata(:)
      CHARACTER(len=*), intent(in) :: dataname
      CHARACTER(len=*), intent(in) :: filename

      ! Local Variables
      TYPE(block_data_real8_2d) :: tmpgriddata2
      TYPE(block_data_real8_3d) :: tmpgriddata3
      logical,  allocatable :: filter (:)
      INTEGER  :: iblkme, xblk, yblk, xloc, yloc, ipatchtype
      REAL(r8) :: wt, tmp

      IF (p_is_master) THEN
         write(*,*) '  Map and Wrtie data on patches: ', dataname
      ENDIF

      IF (p_is_io) THEN 
         call allocate_block_data (grid_patch2grid, tmpgriddata2)  
         call allocate_block_data (grid_patch2grid, tmpgriddata3, npatchtype, 0)  
      ENDIF
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (filter(numpatch))
         ENDIF
      ENDIF
      
      DO ipatchtype = 0, npatchtype-1
         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = (landpatch%ltyp == ipatchtype)
            end if
         end if

         call map_patch2grid%map (patchdata, tmpgriddata2, spv = spval, msk = filter)   

         if (p_is_io) then
            DO iblkme = 1, gblock%nblkme 
               xblk = gblock%xblkme(iblkme)
               yblk = gblock%yblkme(iblkme)

               do yloc = 1, grid_patch2grid%ycnt(yblk) 
                  do xloc = 1, grid_patch2grid%xcnt(xblk) 

                     wt  = wt_patch2grid%blk(xblk,yblk)%val(ipatchtype,xloc,yloc)
                     tmp = tmpgriddata2%blk(xblk,yblk)%val(xloc,yloc)

                     if ((wt /= spval) .and. (tmp /= spval)) THEN
                        tmpgriddata3%blk(xblk,yblk)%val(ipatchtype,xloc,yloc) = tmp / wt
                     else
                        tmpgriddata3%blk(xblk,yblk)%val(ipatchtype,xloc,yloc) = spval
                     end if

                  end do
               end do

            end do
         end if
      ENDDO 

      CALL write3d_patch2grid (tmpgriddata3, dataname, filename)
      
      IF (allocated(filter)) deallocate (filter)

   END SUBROUTINE map_patchdata_to_grid_and_write

   ! -------------
   SUBROUTINE write3d_patch2grid (griddata, dataname, filename)

      IMPLICIT NONE
      
      TYPE(block_data_real8_3d), intent(in) :: griddata
      CHARACTER(len=*), intent(in) :: dataname
      CHARACTER(len=*), intent(in) :: filename
      
      ! Local Variables
      CHARACTER(len=256) :: fullfilename
      INTEGER :: idata, isrc, iblk, jblk
      INTEGER :: ixseg, iyseg, xgdsp, ygdsp, xcnt, ycnt, xbdsp, ybdsp
      REAL(r8), allocatable :: vdata(:,:,:), rbuf(:,:,:), sbuf(:,:,:)
      INTEGER :: smesg(3), rmesg(3)
      LOGICAL :: fexists
      INTEGER :: data_id = 2000

      fullfilename = trim(dir_patch2grid) // '/' // trim(filename)

      if (p_is_master) then
         inquire (file=fullfilename, exist=fexists)
         if (.not. fexists) then
            call ncio_create_file (trim(fullfilename))
            call ncio_define_dimension (fullfilename, 'lat' , concat_patch2grid%ginfo%nlat)
            call ncio_define_dimension (fullfilename, 'lon' , concat_patch2grid%ginfo%nlon)
            call ncio_define_dimension (fullfilename, 'patch', npatchtype) 

            call ncio_write_serial (fullfilename, 'lat',   concat_patch2grid%ginfo%lat_c, 'lat')
            call ncio_write_serial (fullfilename, 'lon',   concat_patch2grid%ginfo%lon_c, 'lon')
            call ncio_write_serial (fullfilename, 'patch', patchindex, 'patch')
         endif
      ENDIF

      data_id = mod(data_id+1,2000)

      if (p_is_master) then
               
         allocate (vdata (npatchtype, concat_patch2grid%ginfo%nlon, concat_patch2grid%ginfo%nlat))
         vdata(:,:,:) = spval
            
#ifdef USEMPI
         do idata = 1, concat_patch2grid%ndatablk

            call mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
               data_id, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            ixseg = rmesg(2)
            iyseg = rmesg(3)

            xgdsp = concat_patch2grid%xsegs(ixseg)%gdsp
            ygdsp = concat_patch2grid%ysegs(iyseg)%gdsp
            xcnt = concat_patch2grid%xsegs(ixseg)%cnt
            ycnt = concat_patch2grid%ysegs(iyseg)%cnt

            allocate (rbuf (npatchtype,xcnt,ycnt))

            call mpi_recv (rbuf, npatchtype * xcnt * ycnt, MPI_DOUBLE, &
               isrc, data_id, p_comm_glb, p_stat, p_err)

            vdata (:,xgdsp+1:xgdsp+xcnt,ygdsp+1:ygdsp+ycnt) = rbuf

            deallocate (rbuf)
         end do
#else

         do iyseg = 1, concat_patch2grid%nyseg
            do ixseg = 1, concat_patch2grid%nxseg
               iblk = concat_patch2grid%xsegs(ixseg)%blk
               jblk = concat_patch2grid%ysegs(iyseg)%blk
               if (gblock%pio(iblk,jblk) == p_iam_glb) then
                  xbdsp = concat_patch2grid%xsegs(ixseg)%bdsp
                  ybdsp = concat_patch2grid%ysegs(iyseg)%bdsp
                  xgdsp = concat_patch2grid%xsegs(ixseg)%gdsp
                  ygdsp = concat_patch2grid%ysegs(iyseg)%gdsp
                  xcnt = concat_patch2grid%xsegs(ixseg)%cnt
                  ycnt = concat_patch2grid%ysegs(iyseg)%cnt

                  vdata (:,xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                     griddata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
               ENDIF
            end do
         ENDDO
#endif

         call ncio_write_serial (fullfilename, dataname, vdata, 'patch', 'lon', 'lat', compress = 1)
         CALL ncio_put_attr (fullfilename, dataname, 'missing_value', spval)

         deallocate (vdata)
      ENDIF

#ifdef USEMPI
      if (p_is_io) then

         do iyseg = 1, concat_patch2grid%nyseg
            do ixseg = 1, concat_patch2grid%nxseg

               iblk = concat_patch2grid%xsegs(ixseg)%blk
               jblk = concat_patch2grid%ysegs(iyseg)%blk

               if (gblock%pio(iblk,jblk) == p_iam_glb) then

                  xbdsp = concat_patch2grid%xsegs(ixseg)%bdsp
                  ybdsp = concat_patch2grid%ysegs(iyseg)%bdsp
                  xcnt = concat_patch2grid%xsegs(ixseg)%cnt
                  ycnt = concat_patch2grid%ysegs(iyseg)%cnt

                  allocate (sbuf (npatchtype,xcnt,ycnt))
                  sbuf = griddata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

                  smesg = (/p_iam_glb, ixseg, iyseg/)
                  call mpi_send (smesg, 3, MPI_INTEGER, p_root, data_id, p_comm_glb, p_err) 
                  call mpi_send (sbuf, npatchtype*xcnt*ycnt, MPI_DOUBLE, &
                     p_root, data_id, p_comm_glb, p_err)

                  deallocate (sbuf)
               end if
            end do
         end do
      end if
#endif

   END SUBROUTINE write3d_patch2grid

END MODULE mod_patch2grid
#endif
