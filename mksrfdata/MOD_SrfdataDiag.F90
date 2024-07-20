#include <define.h>

#ifdef SrfdataDiag
MODULE MOD_SrfdataDiag
!-----------------------------------------------------------------------------------------
! DESCRIPTION:
!
!    This module includes subroutines for checking the results of making surface data.
!
!    The surface data in vector form is mapped to gridded data with last
!    three dimensions of [type,longitude,latitude], which can be viewed by other softwares.
!
!    In GRIDBASED, the grid of gridded data is just the grid of the mesh.
!    In UNSTRUCTURED or CATCHMENT, the grid is user defined and the mapping uses area
!    weighted scheme.
!
! Created by Shupeng Zhang, May 2023
!
! Revisions:
! TODO
!-----------------------------------------------------------------------------------------

   USE MOD_Grid
   USE MOD_SpatialMapping

   IMPLICIT NONE

   ! PUBLIC variables and subroutines
   type(grid_type) :: gdiag
   
   type(spatial_mapping_type) :: m_elm2diag

   type(spatial_mapping_type) :: m_patch2diag
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   type(spatial_mapping_type) :: m_pft2diag
#endif
#ifdef URBAN_MODEL
   type(spatial_mapping_type) :: m_urb2diag
#endif

   PUBLIC :: srfdata_diag_init
   PUBLIC :: srfdata_map_and_write

   ! PRIVATE
   type(grid_concat_type), PRIVATE :: srf_concat
   integer, PRIVATE :: srf_data_id

CONTAINS

   ! ------ SUBROUTINE ------
   SUBROUTINE srfdata_diag_init (dir_landdata)

   USE MOD_SPMD_Task
   USE MOD_LandElm
   USE MOD_LandPatch
#ifdef CROP
   USE MOD_LandCrop
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
#endif
#ifdef URBAN_MODEL
   USE MOD_LandUrban
#endif

   IMPLICIT NONE

   character(len=256), intent(in) :: dir_landdata

   ! Local Variables
   character(len=256) :: landdir, landname
   integer :: ityp
   integer :: typindex(N_land_classification+1)
   real(r8), allocatable :: elmid_r8(:)

      landdir = trim(dir_landdata) // '/diag/'
      IF (p_is_master) THEN
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF

      CALL srf_concat%set (gdiag)
      
      CALL m_elm2diag%build_arealweighted (gdiag, landelm)

      CALL m_patch2diag%build_arealweighted (gdiag, landpatch)

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL m_pft2diag%build_arealweighted (gdiag, landpft)
#endif

#ifdef URBAN_MODEL
      CALL m_urb2diag%build_arealweighted (gdiag, landurban)
#endif

      srf_data_id = 666

      IF (p_is_worker) THEN
         allocate (elmid_r8 (landelm%nset)); elmid_r8 = real(landelm%eindex, r8)
      ENDIF

      landname = trim(dir_landdata)//'/diag/element.nc'
      CALL srfdata_map_and_write (elmid_r8, landelm%settyp, (/0/), m_elm2diag, &
         -1.0e36_r8, landname, 'element', compress = 1, write_mode = 'one')
      
      IF (p_is_worker) deallocate (elmid_r8)

      typindex = (/(ityp, ityp = 0, N_land_classification)/)
      landname = trim(dir_landdata)//'/diag/patchfrac_elm.nc'
      CALL srfdata_map_and_write (elm_patch%subfrc, landpatch%settyp, typindex, m_patch2diag, &
         -1.0e36_r8, landname, 'patchfrac_elm', compress = 1, write_mode = 'one')

#ifdef CATCHMENT
      typindex = (/(ityp, ityp = 0, N_land_classification)/)
      landname = trim(dir_landdata)//'/diag/patchfrac_hru.nc'
      CALL srfdata_map_and_write (hru_patch%subfrc, landpatch%settyp, typindex, m_patch2diag, &
         -1.0e36_r8, landname, 'patchfrac_hru', compress = 1, write_mode = 'one')
#endif

   END SUBROUTINE srfdata_diag_init

   ! ------ SUBROUTINE ------
   SUBROUTINE srfdata_map_and_write ( &
         vsrfdata, settyp, typindex, m_srf, spv, filename, dataname, &
         compress, write_mode, lastdimname, lastdimvalue)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   real(r8), intent(in) :: vsrfdata (:)
   integer , intent(in) :: settyp   (:)
   integer , intent(in) :: typindex (:)

   type(spatial_mapping_type), intent(in) :: m_srf

   real(r8), intent(in) :: spv

   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: dataname
   integer, intent(in) :: compress

   character (len=*), intent(in), optional :: write_mode

   character (len=*), intent(in), optional :: lastdimname
   integer, intent(in), optional :: lastdimvalue

   ! Local variables
   type(block_data_real8_3d) :: wdata, sumwt
   real(r8), allocatable :: vecone (:)

   character(len=10) :: wmode
   integer :: iblkme, ib, jb, iblk, jblk, idata, ixseg, iyseg
   integer :: ntyps, xcnt, ycnt, xbdsp, ybdsp, xgdsp, ygdsp
   integer :: rmesg(3), smesg(3), isrc
   character(len=256) :: fileblock
   real(r8), allocatable :: rbuf(:,:,:), sbuf(:,:,:), vdata(:,:,:)
   logical :: fexists
   integer :: ilastdim

      IF (present(write_mode)) THEN
         wmode = trim(write_mode)
      ELSE
         wmode = 'one'
      ENDIF

      ntyps = size(typindex)

      IF (p_is_io) THEN
         CALL allocate_block_data (gdiag, sumwt, ntyps)
         CALL allocate_block_data (gdiag, wdata, ntyps)
      ENDIF

      IF (p_is_worker) THEN
         IF (size(vsrfdata) > 0) THEN
            allocate (vecone (size(vsrfdata)))
            vecone(:) = 1.0
         ENDIF
      ENDIF

      CALL m_srf%pset2grid_split (vecone  , settyp, typindex, sumwt, spv)
      CALL m_srf%pset2grid_split (vsrfdata, settyp, typindex, wdata, spv)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)

            WHERE ((sumwt%blk(ib,jb)%val > 0.) .and. (wdata%blk(ib,jb)%val /= spv))
               wdata%blk(ib,jb)%val = wdata%blk(ib,jb)%val / sumwt%blk(ib,jb)%val
            elsewhere
               wdata%blk(ib,jb)%val = spv
            END WHERE
         ENDDO
      ENDIF

      IF (trim(wmode) == 'one') THEN

         IF (p_is_master) THEN

            allocate (vdata (ntyps, srf_concat%ginfo%nlon, srf_concat%ginfo%nlat))
            vdata(:,:,:) = spv

#ifdef USEMPI
            DO idata = 1, srf_concat%ndatablk

               CALL mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
                  srf_data_id, p_comm_glb, p_stat, p_err)

               isrc  = rmesg(1)
               ixseg = rmesg(2)
               iyseg = rmesg(3)

               xgdsp = srf_concat%xsegs(ixseg)%gdsp
               ygdsp = srf_concat%ysegs(iyseg)%gdsp
               xcnt  = srf_concat%xsegs(ixseg)%cnt
               ycnt  = srf_concat%ysegs(iyseg)%cnt

               allocate (rbuf (ntyps,xcnt,ycnt))

               CALL mpi_recv (rbuf, ntyps * xcnt * ycnt, MPI_DOUBLE, &
                  isrc, srf_data_id, p_comm_glb, p_stat, p_err)

               vdata (:,xgdsp+1:xgdsp+xcnt,ygdsp+1:ygdsp+ycnt) = rbuf

               deallocate (rbuf)
            ENDDO
#else
            DO iyseg = 1, srf_concat%nyseg
               DO ixseg = 1, srf_concat%nxseg
                  iblk  = srf_concat%xsegs(ixseg)%blk
                  jblk  = srf_concat%ysegs(iyseg)%blk
                  xbdsp = srf_concat%xsegs(ixseg)%bdsp
                  ybdsp = srf_concat%ysegs(iyseg)%bdsp
                  xgdsp = srf_concat%xsegs(ixseg)%gdsp
                  ygdsp = srf_concat%ysegs(iyseg)%gdsp
                  xcnt  = srf_concat%xsegs(ixseg)%cnt
                  ycnt  = srf_concat%ysegs(iyseg)%cnt

                  vdata (:,xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                     wdata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
               ENDDO
            ENDDO
#endif

            write(*,*) 'Please check gridded data < ', trim(dataname), ' > in ', trim(filename)

            inquire (file=trim(filename), exist=fexists)
            IF (.not. fexists) THEN

               CALL ncio_create_file (filename)

               CALL ncio_define_dimension (filename, 'TypeIndex', ntyps)
               CALL ncio_define_dimension (filename, 'lon' , srf_concat%ginfo%nlon)
               CALL ncio_define_dimension (filename, 'lat' , srf_concat%ginfo%nlat)

               CALL ncio_write_serial (filename, 'lat', srf_concat%ginfo%lat_c, 'lat')
               CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
               CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')

               CALL ncio_write_serial (filename, 'lon', srf_concat%ginfo%lon_c, 'lon')
               CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
               CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')

               CALL ncio_write_serial (filename, 'lat_s', srf_concat%ginfo%lat_s, 'lat')
               CALL ncio_put_attr (filename, 'lat_s', 'long_name', 'southern latitude boundary')
               CALL ncio_put_attr (filename, 'lat_s', 'units', 'degrees_north')

               CALL ncio_write_serial (filename, 'lat_n', srf_concat%ginfo%lat_n, 'lat')
               CALL ncio_put_attr (filename, 'lat_n', 'long_name', 'northern latitude boundary')
               CALL ncio_put_attr (filename, 'lat_n', 'units', 'degrees_north')

               CALL ncio_write_serial (filename, 'lon_w', srf_concat%ginfo%lon_w, 'lon')
               CALL ncio_put_attr (filename, 'lon_w', 'long_name', 'western longitude boundary')
               CALL ncio_put_attr (filename, 'lon_w', 'units', 'degrees_east')

               CALL ncio_write_serial (filename, 'lon_e', srf_concat%ginfo%lon_e, 'lon')
               CALL ncio_put_attr (filename, 'lon_e', 'long_name', 'eastern longitude boundary')
               CALL ncio_put_attr (filename, 'lon_e', 'units', 'degrees_east')

               CALL ncio_write_serial (filename, 'TypeIndex', typindex, 'TypeIndex')

            ENDIF

            IF (present(lastdimname) .and. present(lastdimvalue)) THEN
               CALL ncio_write_lastdim (filename, lastdimname, lastdimvalue, ilastdim)
               CALL ncio_write_serial_time (filename, dataname, ilastdim, vdata, &
                  'TypeIndex', 'lon', 'lat', trim(lastdimname), compress)
            ELSE
               CALL ncio_write_serial (filename, dataname, vdata, 'TypeIndex', 'lon', 'lat', compress)
            ENDIF

            CALL ncio_put_attr (filename, dataname, 'missing_value', spv)

            deallocate (vdata)

         ENDIF

#ifdef USEMPI
         IF (p_is_io) THEN

            DO iyseg = 1, srf_concat%nyseg
               DO ixseg = 1, srf_concat%nxseg

                  iblk = srf_concat%xsegs(ixseg)%blk
                  jblk = srf_concat%ysegs(iyseg)%blk

                  IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                     xbdsp = srf_concat%xsegs(ixseg)%bdsp
                     ybdsp = srf_concat%ysegs(iyseg)%bdsp
                     xcnt  = srf_concat%xsegs(ixseg)%cnt
                     ycnt  = srf_concat%ysegs(iyseg)%cnt

                     allocate (sbuf (ntyps,xcnt,ycnt))
                     sbuf = wdata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

                     smesg = (/p_iam_glb, ixseg, iyseg/)
                     CALL mpi_send (smesg, 3, MPI_INTEGER, &
                        p_root, srf_data_id, p_comm_glb, p_err)
                     CALL mpi_send (sbuf, ntyps*xcnt*ycnt, MPI_DOUBLE, &
                        p_root, srf_data_id, p_comm_glb, p_err)

                     deallocate (sbuf)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
#endif

         srf_data_id = srf_data_id + 1

      ELSEIF (trim(wmode) == 'block') THEN

         IF (p_is_io) THEN

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               IF ((gdiag%xcnt(iblk) == 0) .or. (gdiag%ycnt(jblk) == 0)) CYCLE

               CALL get_filename_block (filename, iblk, jblk, fileblock)

               inquire (file=trim(filename), exist=fexists)
               IF (.not. fexists) THEN
                  CALL ncio_create_file (fileblock)
                  CALL ncio_define_dimension (fileblock, 'TypeIndex', ntyps)
                  CALL srf_write_grid_info   (fileblock, gdiag, iblk, jblk)
               ENDIF

               IF (present(lastdimname) .and. present(lastdimvalue)) THEN
                  CALL ncio_write_lastdim (fileblock, lastdimname, lastdimvalue, ilastdim)
                  CALL ncio_write_serial_time (fileblock, dataname, ilastdim, wdata%blk(iblk,jblk)%val, &
                     'TypeIndex', 'lon', 'lat', trim(lastdimname), compress)
               ELSE
                  CALL ncio_write_serial (fileblock, dataname, &
                     wdata%blk(iblk,jblk)%val, 'TypeIndex', 'lon', 'lat', compress)
               ENDIF

               CALL ncio_put_attr (fileblock, dataname, 'missing_value', spv)

            ENDDO

         ENDIF
      ENDIF

      IF (allocated(vecone)) deallocate(vecone)

   END SUBROUTINE srfdata_map_and_write

   !------------------
   SUBROUTINE srf_write_grid_info (fileblock, grid, iblk, jblk)

   USE MOD_Block
   USE MOD_Grid
   USE MOD_NetCDFSerial
   IMPLICIT NONE

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

      IF (grid%xdsp(iblk) + grid%xcnt(iblk) > grid%nlon) THEN
         xl = grid%xdsp(iblk) + 1
         xu = grid%nlon
         nx = grid%nlon - grid%xdsp(iblk)
         lon_w(1:nx) = grid%lon_w(xl:xu)
         lon_e(1:nx) = grid%lon_e(xl:xu)

         xl = 1
         xu = grid%xcnt(iblk) - nx
         lon_w(nx+1:grid%xcnt(iblk)) = grid%lon_w(xl:xu)
         lon_e(nx+1:grid%xcnt(iblk)) = grid%lon_e(xl:xu)
      ELSE
         xl = grid%xdsp(iblk) + 1
         xu = grid%xdsp(iblk) + grid%xcnt(iblk)
         lon_w = grid%lon_w(xl:xu)
         lon_e = grid%lon_e(xl:xu)
      ENDIF

      CALL ncio_define_dimension (fileblock, 'lat', grid%ycnt(jblk))
      CALL ncio_define_dimension (fileblock, 'lon', grid%xcnt(iblk))
      CALL ncio_write_serial (fileblock, 'lat_s', lat_s, 'lat')
      CALL ncio_write_serial (fileblock, 'lat_n', lat_n, 'lat')
      CALL ncio_write_serial (fileblock, 'lon_w', lon_w, 'lon')
      CALL ncio_write_serial (fileblock, 'lon_e', lon_e, 'lon')

   END SUBROUTINE srf_write_grid_info

END MODULE MOD_SrfdataDiag
#endif
