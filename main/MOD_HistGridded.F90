#include <define.h>

module MOD_HistGridded

   !----------------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !     Write out gridded model results to history files.
   !
   ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
   !
   ! REVISIONS:
   ! Shupeng Zhang, 05/2023: 1) porting codes to MPI parallel version
   !
   ! TODO...(need complement)
   !----------------------------------------------------------------------------

   use MOD_Precision
   use MOD_Grid
   use MOD_Mapping_Pset2Grid
   USE MOD_Namelist
   USE MOD_NetCDFSerial

   type(grid_type), target :: ghist
   type(mapping_pset2grid_type) :: mp2g_hist
   type(mapping_pset2grid_type) :: mp2g_hist_urb

   TYPE(grid_concat_type) :: hist_concat

   integer :: hist_data_id

!--------------------------------------------------------------------------
contains

   !---------------------------------------
   subroutine hist_gridded_init (dir_hist)

      USE MOD_Vars_Global
      USE MOD_Namelist
      use MOD_Grid
      USE MOD_LandPatch
#ifdef URBAN_MODEL
      USE MOD_LandUrban
#endif
      use MOD_Mapping_Pset2Grid
      use MOD_Vars_1DAccFluxes
      USE MOD_Forcing, only : gforc
#ifdef SinglePoint
      USE MOD_SingleSrfData
#endif
      implicit none

      character(len=*), intent(in) :: dir_hist

      IF (DEF_hist_grid_as_forcing) then
         CALL ghist%define_by_copy (gforc)
      ELSE
         call ghist%define_by_res (DEF_hist_lon_res, DEF_hist_lat_res)
      ENDIF

#ifndef CROP
      call mp2g_hist%build (landpatch, ghist)
#else
      call mp2g_hist%build (landpatch, ghist, pctcrop)
#endif

#ifdef URBAN_MODEL
      CALL mp2g_hist_urb%build (landurban, ghist)
#endif

      call hist_concat%set (ghist)
#ifdef SinglePoint
      hist_concat%ginfo%lat_c(:) = SITE_lat_location
      hist_concat%ginfo%lon_c(:) = SITE_lon_location
#endif

      if (trim(DEF_HIST_mode) == 'one') then
         hist_data_id = 1000
      end if
         
   end subroutine hist_gridded_init

   ! -------
   subroutine flux_map_and_write_2d ( &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units)

      use MOD_Precision
      use MOD_SPMD_Task
      use MOD_Namelist
      use MOD_DataType
      use MOD_Mapping_Pset2Grid
      use MOD_Block
      use MOD_Grid
      use MOD_Vars_1DAccFluxes,  only: nac
      use MOD_Vars_Global, only: spval
      implicit none

      real(r8), intent(inout) :: acc_vec(:)
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer,          intent(in) :: itime_in_file
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units

      type(block_data_real8_2d), intent(in) :: sumarea
      logical, intent(in) :: filter(:)

      ! Local variables
      type(block_data_real8_2d) :: flux_xy_2d
      integer :: iblkme, xblk, yblk, xloc, yloc
      integer :: compress

      if (p_is_worker)  where (acc_vec /= spval)  acc_vec = acc_vec / nac
      IF (p_is_io)      call allocate_block_data (ghist, flux_xy_2d)  

      call mp2g_hist%map (acc_vec, flux_xy_2d, spv = spval, msk = filter)

      if (p_is_io) then
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            do yloc = 1, ghist%ycnt(yblk)
               do xloc = 1, ghist%xcnt(xblk)

                  if (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                     IF (flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                        flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           = flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                     ENDIF
                  else
                     flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                  end if

               end do
            end do

         end do
      end if

      compress = DEF_HIST_COMPRESS_LEVEL
      call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy_2d, compress)

      IF (p_is_master .and. (itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_2d

   ! -------
   subroutine flux_map_and_write_urb_2d ( &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units)

      use MOD_Precision
      use MOD_SPMD_Task
      use MOD_Namelist
      use MOD_DataType
      use MOD_Mapping_Pset2Grid
      use MOD_Block
      use MOD_Grid
      use MOD_Vars_1DAccFluxes,  only: nac
      use MOD_Vars_Global, only: spval
      implicit none

      real(r8), intent(inout) :: acc_vec(:)
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer,          intent(in) :: itime_in_file
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units

      type(block_data_real8_2d), intent(in) :: sumarea
      logical, intent(in) :: filter(:)

      ! Local variables
      type(block_data_real8_2d) :: flux_xy_2d
      integer :: iblkme, xblk, yblk, xloc, yloc
      integer :: compress

      if (p_is_worker)  where (acc_vec /= spval)  acc_vec = acc_vec / nac
      IF (p_is_io)      call allocate_block_data (ghist, flux_xy_2d)  

      call mp2g_hist_urb%map (acc_vec, flux_xy_2d, spv = spval, msk = filter)

      if (p_is_io) then
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            do yloc = 1, ghist%ycnt(yblk)
               do xloc = 1, ghist%xcnt(xblk)

                  if (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                     IF (flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                        flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           = flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                     ENDIF
                  else
                     flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                  end if

               end do
            end do

         end do
      end if

      compress = DEF_HIST_COMPRESS_LEVEL
      call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy_2d, compress)

      IF (p_is_master .and. (itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_urb_2d

   ! -------
   subroutine flux_map_and_write_3d ( &
         acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, sumarea, filter, &
         longname, units)

      use MOD_Precision
      use MOD_SPMD_Task
      use MOD_Namelist
      use MOD_DataType
      use MOD_Mapping_Pset2Grid
      use MOD_Block
      use MOD_Grid
      use MOD_Vars_1DAccFluxes,  only: nac
      use MOD_Vars_Global, only: spval
      implicit none

      real(r8), intent(inout) :: acc_vec(:,:)
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer, intent(in) :: itime_in_file
      character(len=*), intent(in) :: dim1name
      integer, intent(in) :: lb1, ndim1

      type(block_data_real8_2d), intent(in) :: sumarea
      logical, intent(in) :: filter(:)
      character (len=*), intent(in) :: longname
      character (len=*), intent(in) :: units

      ! Local variables
      type(block_data_real8_3d) :: flux_xy_3d
      integer :: iblkme, xblk, yblk, xloc, yloc, i1
      integer :: compress

      if (p_is_worker)  THEN
         where (acc_vec /= spval)  acc_vec = acc_vec / nac
      ENDIF
      IF (p_is_io) then
         call allocate_block_data (ghist, flux_xy_3d, ndim1, lb1)  
      ENDIF

      call mp2g_hist%map (acc_vec, flux_xy_3d, spv = spval, msk = filter)

      if (p_is_io) then
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            do yloc = 1, ghist%ycnt(yblk)
               do xloc = 1, ghist%xcnt(xblk)

                  if (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                     DO i1 = flux_xy_3d%lb1, flux_xy_3d%ub1
                        IF (flux_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) /= spval) THEN
                           flux_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) &
                              = flux_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) &
                              / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                        ENDIF
                     ENDDO
                  else
                     flux_xy_3d%blk(xblk,yblk)%val(:,xloc,yloc) = spval
                  end if

               end do
            end do

         end do
      end if

      compress = DEF_HIST_COMPRESS_LEVEL
      call hist_write_var_real8_3d (file_hist, varname, dim1name, ghist, &
         itime_in_file, flux_xy_3d, compress)

      IF (p_is_master .and. (itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_3d

   ! -------
   subroutine flux_map_and_write_4d ( &
         acc_vec, file_hist, varname, itime_in_file, &
         dim1name, lb1, ndim1, dim2name, lb2, ndim2, &
         sumarea, filter, longname, units)

      use MOD_Precision
      use MOD_SPMD_Task
      use MOD_Namelist
      use MOD_DataType
      use MOD_Mapping_Pset2Grid
      use MOD_Block
      use MOD_Grid
      use MOD_Vars_1DAccFluxes,  only: nac
      use MOD_Vars_Global, only: spval
      implicit none

      real(r8), intent(inout) :: acc_vec(:,:,:)
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer, intent(in) :: itime_in_file
      character(len=*), intent(in) :: dim1name, dim2name
      integer, intent(in) :: lb1, ndim1, lb2, ndim2

      type(block_data_real8_2d), intent(in) :: sumarea
      logical, intent(in) :: filter(:)
      character (len=*), intent(in) :: longname
      character (len=*), intent(in) :: units

      ! Local variables
      type(block_data_real8_4d) :: flux_xy_4d
      integer :: iblkme, xblk, yblk, xloc, yloc, i1, i2
      integer :: compress

      if (p_is_worker) then
         where(acc_vec /= spval)  acc_vec = acc_vec / nac
      end if
      IF (p_is_io) then
         call allocate_block_data (ghist, flux_xy_4d, ndim1, ndim2, lb1 = lb1, lb2 = lb2)
      ENDIF

      call mp2g_hist%map (acc_vec, flux_xy_4d, spv = spval, msk = filter)

      if (p_is_io) then
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            do yloc = 1, ghist%ycnt(yblk)
               do xloc = 1, ghist%xcnt(xblk)

                  if (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                     DO i1 = flux_xy_4d%lb1, flux_xy_4d%ub1
                        DO i2 = flux_xy_4d%lb2, flux_xy_4d%ub2
                           IF (flux_xy_4d%blk(xblk,yblk)%val(i1,i2,xloc,yloc) /= spval) THEN
                              flux_xy_4d%blk(xblk,yblk)%val(i1,i2,xloc,yloc) &
                                 = flux_xy_4d%blk(xblk,yblk)%val(i1,i2,xloc,yloc) &
                                 / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                           ENDIF
                        ENDDO
                     ENDDO
                  else
                     flux_xy_4d%blk(xblk,yblk)%val(:,:,xloc,yloc) = spval
                  end if

               end do
            end do

         end do
      end if

      compress = DEF_HIST_COMPRESS_LEVEL
      call hist_write_var_real8_4d (file_hist, varname, dim1name, dim2name, &
         ghist, itime_in_file, flux_xy_4d, compress)

      IF (p_is_master .and. (itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_4d

   ! -------
   subroutine flux_map_and_write_ln ( &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units)

      use MOD_Precision
      use MOD_SPMD_Task
      use MOD_Namelist
      use MOD_DataType
      use MOD_Mapping_Pset2Grid
      use MOD_Block
      use MOD_Grid
      use MOD_Vars_1DAccFluxes,  only: nac_ln
      use MOD_Vars_Global, only: spval
      implicit none

      real(r8), intent(inout) :: acc_vec(:)
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer, intent(in) :: itime_in_file

      type(block_data_real8_2d), intent(in) :: sumarea
      logical,  intent(in) :: filter(:)
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units

      ! Local variables
      type(block_data_real8_2d) :: flux_xy_2d
      integer :: i, iblkme, xblk, yblk, xloc, yloc
      integer :: compress

      if (p_is_worker) then
         do i = lbound(acc_vec,1), ubound(acc_vec,1)
            if ((acc_vec(i) /= spval) .and. (nac_ln(i) > 0)) then
               acc_vec(i) = acc_vec(i) / nac_ln(i)
            end if
         end do
      end if

      IF (p_is_io) THEN
         call allocate_block_data (ghist, flux_xy_2d)
      ENDIF

      call mp2g_hist%map (acc_vec, flux_xy_2d, spv = spval, msk = filter)

      if (p_is_io) then
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            do yloc = 1, ghist%ycnt(yblk)
               do xloc = 1, ghist%xcnt(xblk)

                  if ((sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) &
                     .and. (flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval)) then
                     flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                        = flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                        / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                  else
                     flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                  end if

               end do
            end do

         end do
      end if

      compress = DEF_HIST_COMPRESS_LEVEL
      call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy_2d, &
         compress)

      IF (p_is_master .and. (itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_ln

   !------------------------------
   subroutine hist_gridded_write_time ( &
         filename, dataname, time, itime)

      use MOD_Namelist
      use MOD_Grid
      use MOD_Block
      use MOD_SPMD_Task
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname

      integer, intent(in)  :: time(3)
      integer, intent(out) :: itime

      ! Local variables
      character(len=256) :: fileblock
      integer :: iblkme, iblk, jblk
      logical :: fexists

      if (trim(DEF_HIST_mode) == 'one') then
         if (p_is_master) then
            inquire (file=filename, exist=fexists)
            if (.not. fexists) then
               call ncio_create_file (trim(filename))
               CALL ncio_define_dimension(filename, 'time', 0)
               call ncio_define_dimension(filename, 'lat' , hist_concat%ginfo%nlat)
               call ncio_define_dimension(filename, 'lon' , hist_concat%ginfo%nlon)

               call ncio_write_serial (filename, 'lat', hist_concat%ginfo%lat_c, 'lat')
               CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
               CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')

               call ncio_write_serial (filename, 'lon', hist_concat%ginfo%lon_c, 'lon')
               CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
               CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')

#ifndef SinglePoint
               call ncio_write_serial (filename, 'lat_s', hist_concat%ginfo%lat_s, 'lat')
               call ncio_write_serial (filename, 'lat_n', hist_concat%ginfo%lat_n, 'lat')
               call ncio_write_serial (filename, 'lon_w', hist_concat%ginfo%lon_w, 'lon')
               call ncio_write_serial (filename, 'lon_e', hist_concat%ginfo%lon_e, 'lon')
#endif
            endif
            call ncio_write_time (filename, dataname, time, itime, DEF_HIST_FREQ)

         ENDIF

      elseif (trim(DEF_HIST_mode) == 'block') then

         if (p_is_io) then

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)
               IF (ghist%ycnt(jblk) <= 0) cycle
               IF (ghist%xcnt(iblk) <= 0) cycle

               call get_filename_block (filename, iblk, jblk, fileblock)

               inquire (file=fileblock, exist=fexists)
               if (.not. fexists) then
                  call ncio_create_file (trim(fileblock))
                  CALL ncio_define_dimension (fileblock, 'time', 0)
                  call hist_write_grid_info  (fileblock, ghist, iblk, jblk)
               end if

               call ncio_write_time (fileblock, dataname, time, itime, DEF_HIST_FREQ)

            end do

         end if
      endif

   end subroutine hist_gridded_write_time

   !----------------------------------------------------------------------------
   subroutine hist_write_var_real8_2d ( &
         filename, dataname, grid, itime, wdata, compress)

      use MOD_Namelist
      use MOD_Block
      use MOD_Grid
      use MOD_DataType
      use MOD_SPMD_Task
      use MOD_Vars_Global, only: spval
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      type (grid_type),  intent(in) :: grid
      integer, intent(in) :: itime

      type (block_data_real8_2d), intent(in) :: wdata

      integer, intent(in) :: compress

      ! Local variables
      integer :: iblkme, iblk, jblk, idata, ixseg, iyseg
      integer :: xcnt, ycnt, xbdsp, ybdsp, xgdsp, ygdsp
      integer :: rmesg(3), smesg(3), isrc
      character(len=256) :: fileblock
      real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)

      if (trim(DEF_HIST_mode) == 'one') then

         if (p_is_master) then

            allocate (vdata (hist_concat%ginfo%nlon, hist_concat%ginfo%nlat))
            vdata(:,:) = spval

#ifdef USEMPI
            do idata = 1, hist_concat%ndatablk
               call mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
                  hist_data_id, p_comm_glb, p_stat, p_err)

               isrc  = rmesg(1)
               ixseg = rmesg(2)
               iyseg = rmesg(3)

               xgdsp = hist_concat%xsegs(ixseg)%gdsp
               ygdsp = hist_concat%ysegs(iyseg)%gdsp
               xcnt = hist_concat%xsegs(ixseg)%cnt
               ycnt = hist_concat%ysegs(iyseg)%cnt

               allocate (rbuf(xcnt,ycnt))

               call mpi_recv (rbuf, xcnt*ycnt, MPI_DOUBLE, &
                  isrc, hist_data_id, p_comm_glb, p_stat, p_err)

               vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = rbuf
               deallocate (rbuf)

            end do
#else
            do iyseg = 1, hist_concat%nyseg
               do ixseg = 1, hist_concat%nxseg
                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk
                  if (gblock%pio(iblk,jblk) == p_iam_glb) then
                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xgdsp = hist_concat%xsegs(ixseg)%gdsp
                     ygdsp = hist_concat%ysegs(iyseg)%gdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt

                     vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                        wdata%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                  ENDIF
               ENDDO
            end do
#endif

            call ncio_write_serial_time (filename, dataname, itime, vdata, &
               'lon', 'lat', 'time', compress)

            deallocate (vdata)
         ENDIF

#ifdef USEMPI
         if (p_is_io) then
            do iyseg = 1, hist_concat%nyseg
               do ixseg = 1, hist_concat%nxseg

                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk

                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt

                     allocate (sbuf (xcnt,ycnt))
                     sbuf = wdata%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

                     smesg = (/p_iam_glb, ixseg, iyseg/)
                     call mpi_send (smesg, 3, MPI_INTEGER, &
                        p_root, hist_data_id, p_comm_glb, p_err)
                     call mpi_send (sbuf, xcnt*ycnt, MPI_DOUBLE, &
                        p_root, hist_data_id, p_comm_glb, p_err)

                     deallocate (sbuf)

                  end if
               end do
            end do
         end if
#endif

         hist_data_id = hist_data_id + 1

      elseif (trim(DEF_HIST_mode) == 'block') then

         if (p_is_io) then

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               if ((grid%xcnt(iblk) == 0) .or. (grid%ycnt(jblk) == 0)) cycle

               call get_filename_block (filename, iblk, jblk, fileblock)

               call ncio_write_serial_time (fileblock, dataname, itime, &
                  wdata%blk(iblk,jblk)%val, 'lon', 'lat', 'time', compress)

            end do

         end if
      end if

   end subroutine hist_write_var_real8_2d

   !----------------------------------------------------------------------------
   subroutine hist_write_var_real8_3d ( &
         filename, dataname, dim1name, grid, itime, wdata, compress)

      use MOD_Namelist
      use MOD_Block
      use MOD_Grid
      use MOD_DataType
      use MOD_SPMD_Task
      use MOD_Vars_Global, only: spval
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      character (len=*), intent(in) :: dim1name
      type (grid_type),  intent(in) :: grid
      integer, intent(in) :: itime

      type (block_data_real8_3d), intent(in) :: wdata

      integer, intent(in) :: compress

      ! Local variables
      integer :: iblkme, iblk, jblk, idata, ixseg, iyseg
      integer :: xcnt, ycnt, ndim1, xbdsp, ybdsp, xgdsp, ygdsp
      integer :: rmesg(4), smesg(4), isrc
      character(len=256) :: fileblock
      real(r8), allocatable :: rbuf(:,:,:), sbuf(:,:,:), vdata(:,:,:)

      if (trim(DEF_HIST_mode) == 'one') then

         if (p_is_master) then

#ifdef USEMPI
            do idata = 1, hist_concat%ndatablk

               call mpi_recv (rmesg, 4, MPI_INTEGER, MPI_ANY_SOURCE, &
                  hist_data_id, p_comm_glb, p_stat, p_err)

               isrc  = rmesg(1)
               ixseg = rmesg(2)
               iyseg = rmesg(3)
               ndim1 = rmesg(4)

               xgdsp = hist_concat%xsegs(ixseg)%gdsp
               ygdsp = hist_concat%ysegs(iyseg)%gdsp
               xcnt = hist_concat%xsegs(ixseg)%cnt
               ycnt = hist_concat%ysegs(iyseg)%cnt

               allocate (rbuf (ndim1,xcnt,ycnt))

               call mpi_recv (rbuf, ndim1 * xcnt * ycnt, MPI_DOUBLE, &
                  isrc, hist_data_id, p_comm_glb, p_stat, p_err)

               IF (idata == 1) THEN
                  allocate (vdata (ndim1, hist_concat%ginfo%nlon, hist_concat%ginfo%nlat))
                  vdata(:,:,:) = spval
               ENDIF

               vdata (:,xgdsp+1:xgdsp+xcnt,ygdsp+1:ygdsp+ycnt) = rbuf

               deallocate (rbuf)
            end do
#else
            ndim1 = wdata%ub1 - wdata%lb1 + 1
            allocate (vdata (ndim1, hist_concat%ginfo%nlon, hist_concat%ginfo%nlat))
            vdata(:,:,:) = spval

            do iyseg = 1, hist_concat%nyseg
               do ixseg = 1, hist_concat%nxseg
                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk
                  if (gblock%pio(iblk,jblk) == p_iam_glb) then
                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xgdsp = hist_concat%xsegs(ixseg)%gdsp
                     ygdsp = hist_concat%ysegs(iyseg)%gdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt

                     vdata (:,xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                        wdata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                  ENDIF
               end do
            ENDDO
#endif

            call ncio_define_dimension (filename, dim1name, ndim1)

            call ncio_write_serial_time (filename, dataname, itime, &
               vdata, dim1name, 'lon', 'lat', 'time', compress)

            deallocate (vdata)
         ENDIF

#ifdef USEMPI
         if (p_is_io) then

            do iyseg = 1, hist_concat%nyseg
               do ixseg = 1, hist_concat%nxseg

                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk

                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt
                     ndim1 = size(wdata%blk(iblk,jblk)%val,1)

                     allocate (sbuf (ndim1,xcnt,ycnt))
                     sbuf = wdata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

                     smesg = (/p_iam_glb, ixseg, iyseg, ndim1/)
                     call mpi_send (smesg, 4, MPI_INTEGER, &
                        p_root, hist_data_id, p_comm_glb, p_err)
                     call mpi_send (sbuf, ndim1*xcnt*ycnt, MPI_DOUBLE, &
                        p_root, hist_data_id, p_comm_glb, p_err)

                     deallocate (sbuf)
                  end if
               end do
            end do
         end if
#endif

         hist_data_id = hist_data_id + 1

      elseif (trim(DEF_HIST_mode) == 'block') then

         if (p_is_io) then

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               if ((grid%xcnt(iblk) == 0) .or. (grid%ycnt(jblk) == 0)) cycle

               call get_filename_block (filename, iblk, jblk, fileblock)

               call ncio_define_dimension (fileblock, dim1name, wdata%ub1-wdata%lb1+1)

               call ncio_write_serial_time (fileblock, dataname, itime, &
                  wdata%blk(iblk,jblk)%val, dim1name, 'lon', 'lat', 'time', compress)

            end do

         end if
      end if

   end subroutine hist_write_var_real8_3d

   !----------------------------------------------------------------------------
   subroutine hist_write_var_real8_4d ( &
         filename, dataname, dim1name, dim2name, grid, itime, wdata, compress)

      use MOD_Namelist
      use MOD_Block
      use MOD_Grid
      use MOD_DataType
      use MOD_SPMD_Task
      use MOD_Vars_Global, only: spval
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      character (len=*), intent(in) :: dim1name, dim2name
      type (grid_type),  intent(in) :: grid
      integer, intent(in) :: itime

      type (block_data_real8_4d), intent(in) :: wdata

      integer, intent(in) :: compress

      ! Local variables
      integer :: iblkme, iblk, jblk, idata, ixseg, iyseg
      integer :: xcnt, ycnt, ndim1, ndim2, xbdsp, ybdsp, xgdsp, ygdsp
      integer :: rmesg(5), smesg(5), isrc
      character(len=256) :: fileblock
      real(r8), allocatable :: rbuf(:,:,:,:), sbuf(:,:,:,:), vdata(:,:,:,:)

      if (trim(DEF_HIST_mode) == 'one') then

         if (p_is_master) then

#ifdef USEMPI
            do idata = 1, hist_concat%ndatablk

               call mpi_recv (rmesg, 5, MPI_INTEGER, MPI_ANY_SOURCE, &
                  hist_data_id, p_comm_glb, p_stat, p_err)

               isrc  = rmesg(1)
               ixseg = rmesg(2)
               iyseg = rmesg(3)
               ndim1 = rmesg(4)
               ndim2 = rmesg(4)

               xgdsp = hist_concat%xsegs(ixseg)%gdsp
               ygdsp = hist_concat%ysegs(iyseg)%gdsp
               xcnt = hist_concat%xsegs(ixseg)%cnt
               ycnt = hist_concat%ysegs(iyseg)%cnt

               allocate (rbuf (ndim1,ndim2,xcnt,ycnt))

               call mpi_recv (rbuf, ndim1*ndim2*xcnt*ycnt, MPI_DOUBLE, &
                  isrc, hist_data_id, p_comm_glb, p_stat, p_err)

               IF (idata == 1) THEN
                  allocate (vdata (ndim1,ndim2,hist_concat%ginfo%nlon,hist_concat%ginfo%nlat))
                  vdata(:,:,:,:) = spval
               ENDIF

               vdata (:,:,xgdsp+1:xgdsp+xcnt,ygdsp+1:ygdsp+ycnt) = rbuf

               deallocate (rbuf)
            end do
#else
            ndim1 = wdata%ub1 - wdata%lb1 + 1
            ndim2 = wdata%ub2 - wdata%lb2 + 1
            allocate (vdata (ndim1,ndim2,hist_concat%ginfo%nlon,hist_concat%ginfo%nlat))
            vdata(:,:,:,:) = spval

            do iyseg = 1, hist_concat%nyseg
               do ixseg = 1, hist_concat%nxseg
                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk
                  if (gblock%pio(iblk,jblk) == p_iam_glb) then
                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xgdsp = hist_concat%xsegs(ixseg)%gdsp
                     ygdsp = hist_concat%ysegs(iyseg)%gdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt

                     vdata (:,:,xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                        wdata%blk(iblk,jblk)%val(:,:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                  ENDIF
               ENDDO
            ENDDO

#endif

            call ncio_define_dimension (filename, dim1name, ndim1)
            call ncio_define_dimension (filename, dim2name, ndim2)

            call ncio_write_serial_time (filename, dataname, itime, vdata, dim1name, dim2name, &
                  'lon', 'lat', 'time', compress)

            deallocate (vdata)
         ENDIF

#ifdef USEMPI
         if (p_is_io) then

            do iyseg = 1, hist_concat%nyseg
               do ixseg = 1, hist_concat%nxseg

                  iblk = hist_concat%xsegs(ixseg)%blk
                  jblk = hist_concat%ysegs(iyseg)%blk

                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     xbdsp = hist_concat%xsegs(ixseg)%bdsp
                     ybdsp = hist_concat%ysegs(iyseg)%bdsp
                     xcnt = hist_concat%xsegs(ixseg)%cnt
                     ycnt = hist_concat%ysegs(iyseg)%cnt

                     ndim1 = size(wdata%blk(iblk,jblk)%val,1)
                     ndim2 = size(wdata%blk(iblk,jblk)%val,2)
                     allocate (sbuf (ndim1,ndim2,xcnt,ycnt))
                     sbuf = wdata%blk(iblk,jblk)%val(:,:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

                     smesg = (/p_iam_glb, ixseg, iyseg, ndim1, ndim2/)
                     call mpi_send (smesg, 5, MPI_INTEGER, &
                        p_root, hist_data_id, p_comm_glb, p_err)
                     call mpi_send (sbuf, ndim1*ndim2*xcnt*ycnt, MPI_DOUBLE, &
                        p_root, hist_data_id, p_comm_glb, p_err)

                     deallocate (sbuf)
                  end if
               end do
            end do
         end if
#endif

         hist_data_id = hist_data_id + 1

      elseif (trim(DEF_HIST_mode) == 'block') then
         if (p_is_io) then

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               if ((grid%xcnt(iblk) == 0) .or. (grid%ycnt(jblk) == 0)) cycle

               call get_filename_block (filename, iblk, jblk, fileblock)

               call ncio_define_dimension (fileblock, dim1name, wdata%ub1-wdata%lb1+1)
               call ncio_define_dimension (fileblock, dim2name, wdata%ub2-wdata%lb2+1)

               call ncio_write_serial_time (fileblock, dataname, itime, &
                  wdata%blk(iblk,jblk)%val, dim1name, dim2name, 'lon', 'lat', 'time', compress)

            end do

         end if
      end if

   end subroutine hist_write_var_real8_4d

   !------------------
   subroutine hist_write_grid_info (fileblock, grid, iblk, jblk)

      use MOD_Block
      use MOD_Grid
      implicit none

      character(len=*), intent(in) :: fileblock
      type (grid_type), intent(in) :: grid
      integer, intent(in) :: iblk, jblk

      ! Local variable
      integer :: yl, yu, xl, xu, nx
      real(r8), allocatable :: lat_s(:), lat_n(:), lon_w(:), lon_e(:)

      allocate (lon_w (grid%xcnt(iblk)))
      allocate (lon_e (grid%xcnt(iblk)))
      allocate (lat_s (grid%ycnt(jblk)))
      allocate (lat_n (grid%ycnt(jblk)))

      yl = grid%ydsp(jblk) + 1
      yu = grid%ydsp(jblk) + grid%ycnt(jblk)

      lat_s = grid%lat_s(yl:yu)
      lat_n = grid%lat_n(yl:yu)

      if (grid%xdsp(iblk) + grid%xcnt(iblk) > grid%nlon) then
         xl = grid%xdsp(iblk) + 1
         xu = grid%nlon
         nx = grid%nlon - grid%xdsp(iblk)
         lon_w(1:nx) = grid%lon_w(xl:xu)
         lon_e(1:nx) = grid%lon_e(xl:xu)

         xl = 1
         xu = grid%xcnt(iblk) - nx
         lon_w(nx+1:grid%xcnt(iblk)) = grid%lon_w(xl:xu)
         lon_e(nx+1:grid%xcnt(iblk)) = grid%lon_e(xl:xu)
      else
         xl = grid%xdsp(iblk) + 1
         xu = grid%xdsp(iblk) + grid%xcnt(iblk)
         lon_w = grid%lon_w(xl:xu)
         lon_e = grid%lon_e(xl:xu)
      end if

      CALL ncio_define_dimension (fileblock, 'lat', grid%ycnt(jblk))
      CALL ncio_define_dimension (fileblock, 'lon', grid%xcnt(iblk))
      call ncio_write_serial (fileblock, 'lat_s', lat_s, 'lat')
      call ncio_write_serial (fileblock, 'lat_n', lat_n, 'lat')
      call ncio_write_serial (fileblock, 'lon_w', lon_w, 'lon')
      call ncio_write_serial (fileblock, 'lon_e', lon_e, 'lon')

   end subroutine hist_write_grid_info

end module MOD_HistGridded
