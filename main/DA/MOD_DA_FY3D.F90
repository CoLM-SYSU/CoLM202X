#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_FY3D
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Data assimilation of brightness temperature from FY3D satellite
!
! AUTHOR:
!   Lu Li, 07/2025: Initial version
!-----------------------------------------------------------------------------
   USE MOD_DataType
   USE MOD_SpatialMapping
   USE MOD_DA_ObsOperator
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Vars_Global, only: pi, nl_soil
   USE MOD_LandPatch
   USE MOD_Block
   USE MOD_Namelist
   IMPLICIT NONE
   SAVE

! public functions
   PUBLIC :: init_DA_FY3D
   PUBLIC :: run_DA_FY3D
   PUBLIC :: end_DA_FY3D

   PRIVATE

! local variables
   ! file path
   character(len=256) :: file_fy3d                           ! FY3D file path
   character(len=256) :: file_grid                           ! FY3D world grid file path

   ! grid
   type(grid_type) :: grid_fy3d                              ! FY3D world grid
   type(spatial_mapping_type) :: mg2p_fy3d                   ! mapping between world grid to patch

   ! time (UTC) at current step
   integer :: month, mday, hour                              ! month, day, hour of current step
   character(len=256) :: yearstr, monthstr, daystr, hourstr  ! string of year, month, day, hour
   integer :: idate_b, idate_e                               ! begin & end seconds since begin of current day (UTC)

   ! time variables used to determine whether has obs
   real(r8), allocatable :: fy3d_time(:)                     ! seconds of all obs since begin of current day (UTC)
   real(r8), allocatable :: dt_b(:)                          ! delta time between obs and begin seconds of current day
   real(r8), allocatable :: dt_e(:)                          ! delta time between obs and end seconds of current day

   ! logical variables
   logical :: has_file                                       ! whether has file of FY3D at target day
   logical :: has_obs                                        ! whether has obs at current step
   logical :: has_DA                                         ! whether has data assimilation 

   ! observations (dimensions changes with time)
   integer :: num_obs                                        ! number of all obs in current file
   real(r8), allocatable :: fy3d_lat(:)                      ! latitude of all obs
   real(r8), allocatable :: fy3d_lon(:)                      ! longitude of all obs
   real(r8), allocatable :: fy3d_tb_h(:)                     ! H- polarized brightness temperature of all obs ([K])
   real(r8), allocatable :: fy3d_tb_v(:)                     ! V- polarized brightness temperature of all obs ([K])
   integer,  allocatable :: fy3d_ii(:)                       ! i-th lat grid in world map of all obs
   integer,  allocatable :: fy3d_jj(:)                       ! j-th lon grid in world map of all obs
   
   ! setting of FY3D satellite
   real(r8), parameter :: fy3d_theta = 53.0*pi/180           ! incidence angle (rad)
   real(r8), parameter :: fy3d_fghz = 10.65                  ! frequency (GHz), X-band

   ! ensemble predicted observations (at patch)
   real(r8), allocatable :: pred_tb_h_pset_ens(:, :)         ! predicted H- polarized temp on patch
   real(r8), allocatable :: pred_tb_v_pset_ens(:, :)         ! predicted V- polarized temp on patch
   real(r8), allocatable :: area_pset(:)                     ! area of each patch across world grid

   ! ensemble predicted observations (at world grid)
   type(block_data_real8_3d) :: pred_tb_h_wgrid_ens          ! predicted H- polarized temp on world grid
   type(block_data_real8_3d) :: pred_tb_v_wgrid_ens          ! predicted V- polarized temp on world grid
   type(block_data_real8_2d) :: area_wgrid                   ! area of each patch across world grid

   ! ensemble predicted observations (at obs grid)
   real(r8), allocatable :: pred_tb_h_ogrid_ens(:, :)        ! predicted H- polarized temp on obs grid
   real(r8), allocatable :: pred_tb_v_ogrid_ens(:, :)        ! predicted V- polarized temp on obs grid
   real(r8) :: area_wgrid_obs                                ! area of each patch across world grid for each obs

   ! observations around patch (dimensions changes with patch)
   integer :: num_obs_p
   logical, allocatable :: index_p(:)                        ! index of obs around each patch
   real(r8), allocatable :: fy3d_lat_p(:)                    ! latitude of obs around each patch
   real(r8), allocatable :: fy3d_lon_p(:)                    ! longitude of obs around each patch
   real(r8), allocatable :: fy3d_tb_h_p(:)                   ! H- polarized brightness temperature of obs around each patch ([K])
   real(r8), allocatable :: fy3d_tb_v_p(:)                   ! V- polarized brightness temperature of obs around each patch ([K])
   real(r8), allocatable :: pred_tb_h_p_ens(:, :)            ! predicted H- polarized temp around patch
   real(r8), allocatable :: pred_tb_v_p_ens(:, :)            ! predicted V- polarized temp around patch
   real(r8), allocatable :: d_p(:)                           ! distance between obs and patch center

   ! parameters of LETKF
   real(r8), allocatable :: obs_err(:)                       ! observation error
   real(r8), parameter   :: dres = 0.4                       ! search localization radius (deg)
   real(r8), parameter   :: loc_r = 1.0                      ! localization radius
   real(r8), parameter   :: infl = 1.2                       ! inflation factor
   real(r8), parameter   :: static_obs_err = 0.1             ! static observation error

   ! temporary variables for data assimilation
   real(r8), allocatable :: trans(:,:)                       ! transformation matrix on each patch
   real(r8), allocatable :: wice_soi_ens(:,:)                ! soil ice content
   real(r8), allocatable :: wice_soi_ens_da(:,:)             ! soil ice content after data assimilation
   real(r8), allocatable :: wliq_soi_ens(:,:)                ! soil liquid water content
   real(r8), allocatable :: wliq_soi_ens_da(:,:)             ! soil liquid water content after data assimilation
   logical, allocatable ::  filter(:)                        ! to mask the water
   real(r8) :: eff_porsl                                     ! effective porosity of soil layer

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_FY3D()

!-----------------------------------------------------------------------------
      USE MOD_Spmd_Task
      USE MOD_Namelist, only: DEF_DA_obsdir
      USE MOD_Grid
      USE MOD_NetCDFSerial
      USE MOD_LandPatch
      USE MOD_Pixelset
      USE MOD_RangeCheck
      IMPLICIT NONE

!-----------------------------------------------------------------------

      ! grid file path of FY3D, 0.25 degree world grid
      file_grid = trim(DEF_DA_obsdir)//'/FY3D_LatLon_M25km_pre.nc'

#ifndef SinglePoint
      CALL grid_fy3d%define_from_file(file_grid, 'latitude', 'longitude')

      ! map FY3D grid to patch
      CALL mg2p_fy3d%build_arealweighted(grid_fy3d, landpatch)
#endif

   END SUBROUTINE init_DA_FY3D

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA_FY3D(idate, deltim)

!-----------------------------------------------------------------------------
      USE MOD_Spmd_task
      USE MOD_TimeManager
      USE MOD_NetCDFBlock
      USE MOD_Mesh
      USE MOD_LandElm
      USE MOD_LandPatch
      USE MOD_Vars_1DFluxes
      USE MOD_Vars_1DForcing
      USE MOD_Vars_TimeVariables
      USE MOD_Vars_TimeInvariants
      USE MOD_DA_Vars_TimeVariables
      USE MOD_RangeCheck
      USE MOD_UserDefFun
      USE MOD_DA_EnKF
      USE MOD_Const_Physical, only: denice, denh2o
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer, intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
      real(r8) :: lat_p_n, lat_p_s, lon_p_w, lon_p_e
      integer  :: ib, jb, il, jl, ip, iens, iobs, np, i, n, ndata, isrc
      integer  :: sdate(3), smesg(2), rmesg(2)
      integer, allocatable  :: iloc(:)               
      integer, allocatable  :: itemp(:)
      real(r8), allocatable :: dtemp(:,:)

!-----------------------------------------------------------------------------

!#############################################################################
! Identify if there are observations at this time step 
!#############################################################################
      ! Do not perform DA, only calcuate predict BRT for diagnostic
      IF (DEF_DA_ENS == 1) THEN
         has_file = .false.
         has_obs = .false.
      ELSE
         ! covert local time to UTC for single point
         sdate = idate
         CALL adj2begin(sdate)
#ifdef SinglePoint
         IF (.not. DEF_simulation_time%greenwich) THEN
            CALL localtime2gmt(sdate)
         ENDIF
#endif
         ! calculate year/month/day/hour of current step
         CALL julian2monthday(sdate(1), sdate(2), month, mday)
         hour = int(sdate(3)/3600)

         ! whether has file of FY3D at target day
         write (yearstr, '(I4.4)') sdate(1)
         write (monthstr, '(I2.2)') month
         write (daystr, '(I2.2)') mday
         write (hourstr, '(I2.2)') hour
      file_fy3d = trim(DEF_DA_obsdir)//'/pre/'//'/FY3D_L1_TB_'// &
         trim(yearstr)//'_'//trim(monthstr)//'_'//trim(daystr)//'_'//trim(hourstr)//'.nc'
      inquire (file=trim(file_fy3d), exist=has_file)

         ! whether have obs at this time interval
         has_obs = .false.
         IF (has_file) THEN
            CALL ncio_read_bcast_serial(file_fy3d, 'time', fy3d_time)
            num_obs = size(fy3d_time)
            idate_b = sdate(3)
            idate_e = sdate(3) + deltim
            allocate (dt_b(num_obs))
            allocate (dt_e(num_obs))
            dt_b = fy3d_time - idate_b
            dt_e = fy3d_time - idate_e
            IF (any(dt_b >= 0 .and. dt_e <= 0)) has_obs = .true.
            deallocate (dt_b)
            deallocate (dt_e)
         ELSE
            has_obs = .false.
         ENDIF
      ENDIF


!#############################################################################
! Allocate memory for variables 
!#############################################################################
      ! allocate memory for ensemble predicted observations at patch (for DA or no DA)
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            IF (allocated(pred_tb_h_pset_ens)) deallocate (pred_tb_h_pset_ens)
            IF (allocated(pred_tb_v_pset_ens)) deallocate (pred_tb_v_pset_ens)
            allocate (pred_tb_h_pset_ens(DEF_DA_ENS, numpatch))
            allocate (pred_tb_v_pset_ens(DEF_DA_ENS, numpatch))
         ENDIF
      ENDIF

      ! allocate memory for ensemble predicted observations in vectors (only for DA)
      IF (has_obs) THEN
#ifndef SinglePoint 
         IF (allocated(pred_tb_h_ogrid_ens)) deallocate (pred_tb_h_ogrid_ens)
         IF (allocated(pred_tb_v_ogrid_ens)) deallocate (pred_tb_v_ogrid_ens)
         IF (allocated(area_pset)) deallocate (area_pset)
         allocate (pred_tb_h_ogrid_ens(num_obs, DEF_DA_ENS))
         allocate (pred_tb_v_ogrid_ens(num_obs, DEF_DA_ENS))
         allocate (area_pset(numpatch))
#endif
      ENDIF


!#############################################################################
! Calculate predicted observations & mapping to world grid & read observations
!#############################################################################
      ! forward model (for DA or no DA)
      IF (p_is_worker) THEN
         DO iens = 1, DEF_DA_ENS
            DO np = 1, numpatch
               CALL forward( &
                  patchtype(np), patchclass(np), dz_sno_ens(:, iens, np), &
                  forc_topo(np), htop(np), &
                  tref_ens(iens, np), t_soisno_ens(:, iens, np), tleaf_ens(iens, np), &
                  wliq_soisno_ens(:, iens, np), wice_soisno_ens(:, iens, np), h2osoi_ens(:, iens, np), &
                  snowdp_ens(iens, np), lai_ens(iens, np), sai_ens(iens, np), &
                  vf_clay(:, np), vf_sand(:, np), BD_all(:, np), porsl(:, np), &
                  fy3d_theta, fy3d_fghz, &
                  pred_tb_h_pset_ens(iens, np), pred_tb_v_pset_ens(iens, np))
            ENDDO
         ENDDO
         t_brt_ens(1,:,:) = pred_tb_h_pset_ens
         t_brt_ens(2,:,:) = pred_tb_v_pset_ens
      ENDIF

      ! read observations data (only for DA)
      IF (has_obs) THEN
         CALL ncio_read_bcast_serial(file_fy3d, 'lat' , fy3d_lat )
         CALL ncio_read_bcast_serial(file_fy3d, 'lon' , fy3d_lon )
         CALL ncio_read_bcast_serial(file_fy3d, 'ii'  , fy3d_ii  )
         CALL ncio_read_bcast_serial(file_fy3d, 'jj'  , fy3d_jj  )
         CALL ncio_read_bcast_serial(file_fy3d, 'tb_h', fy3d_tb_h)
         CALL ncio_read_bcast_serial(file_fy3d, 'tb_v', fy3d_tb_v)
      ENDIF

      ! map predicted observation from patch to world grid & crop useful data (only for DA)
      IF (has_obs) THEN
#ifndef SinglePoint
         ! mapping predicted observations from patch to world grid
         IF (p_is_io) THEN
            CALL allocate_block_data(grid_fy3d, pred_tb_h_wgrid_ens, DEF_DA_ENS)
            CALL allocate_block_data(grid_fy3d, pred_tb_v_wgrid_ens, DEF_DA_ENS)
         ENDIF
         CALL mg2p_fy3d%pset2grid(pred_tb_h_pset_ens, pred_tb_h_wgrid_ens, spval)
         CALL mg2p_fy3d%pset2grid(pred_tb_v_pset_ens, pred_tb_v_wgrid_ens, spval)

         ! calculate area of each patch across world grid
         area_pset(:) = 1
         allocate (filter(numpatch))
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               filter(:) = patchtype <= 2 
            ENDIF
         ENDIF
         IF (p_is_io) CALL allocate_block_data(grid_fy3d, area_wgrid)
         CALL mg2p_fy3d%pset2grid(area_pset, area_wgrid, msk=filter)
         deallocate (filter)

         ! crop the predicted observations 
         IF (p_is_io) THEN
            allocate (iloc(num_obs))
            pred_tb_h_ogrid_ens = -9999.0
            pred_tb_v_ogrid_ens = -9999.0

            ! crop obs from world grid to all obs grids
            ndata = 0
            DO i = 1, num_obs
               !(Notes: Python index starts from 0)
               ib = grid_fy3d%xblk(fy3d_jj(i) + 1)
               jb = grid_fy3d%yblk(fy3d_ii(i) + 1)
               il = grid_fy3d%xloc(fy3d_jj(i) + 1)
               jl = grid_fy3d%yloc(fy3d_ii(i) + 1)
               IF (ib /= 0 .and. jb /= 0) THEN
                  IF (gblock%pio(ib, jb) == p_iam_glb) THEN
                     area_wgrid_obs = area_wgrid%blk(ib, jb)%val(il, jl)
                     IF (area_wgrid_obs /= 0) THEN
                        ndata = ndata + 1
                        iloc(ndata) = i
                        pred_tb_h_ogrid_ens(ndata, :) = (pred_tb_h_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_wgrid_obs
                        pred_tb_v_ogrid_ens(ndata, :) = (pred_tb_v_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_wgrid_obs
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

            ! send data from io to masters
#ifdef USEMPI
            smesg = (/p_iam_glb, ndata/)
            CALL mpi_send(smesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

            IF (ndata > 0) THEN
               CALL mpi_send(iloc(1:ndata), ndata, MPI_INTEGER, p_address_master, mpi_tag_data, p_comm_glb, p_err)

               allocate (dtemp(ndata, DEF_DA_ENS))
               dtemp = pred_tb_h_ogrid_ens(1:ndata, :)
               CALL mpi_send(dtemp, ndata*DEF_DA_ENS, MPI_DOUBLE, p_address_master, mpi_tag_data, p_comm_glb, p_err)
               dtemp = pred_tb_v_ogrid_ens(1:ndata, :)
               CALL mpi_send(dtemp, ndata*DEF_DA_ENS, MPI_DOUBLE, p_address_master, mpi_tag_data + 1, p_comm_glb, p_err)
               deallocate (dtemp)
            ENDIF
#endif
            deallocate (iloc)
         ENDIF

         ! broadcast from master to all workers
#ifdef USEMPI
         IF (p_is_master) THEN
            pred_tb_h_ogrid_ens = -9999.0
            pred_tb_v_ogrid_ens = -9999.0
            DO ip = 0, p_np_io - 1
               CALL mpi_recv(rmesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

               ndata = rmesg(2)
               IF (ndata > 0) THEN
                  allocate (itemp(ndata))
                  allocate (dtemp(ndata, DEF_DA_ENS))

                  isrc = rmesg(1)
                  CALL mpi_recv(itemp, ndata, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  CALL mpi_recv(dtemp, ndata*DEF_DA_ENS, MPI_DOUBLE, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                  pred_tb_h_ogrid_ens(itemp, :) = dtemp
                  CALL mpi_recv(dtemp, ndata*DEF_DA_ENS, MPI_DOUBLE, isrc, mpi_tag_data + 1, p_comm_glb, p_stat, p_err)
                  pred_tb_v_ogrid_ens(itemp, :) = dtemp

                  deallocate (itemp)
                  deallocate (dtemp)
               ENDIF
            ENDDO
         ENDIF
         CALL mpi_bcast(pred_tb_h_ogrid_ens, num_obs*DEF_DA_ENS, MPI_DOUBLE, p_address_master, p_comm_glb, p_err)
         CALL mpi_bcast(pred_tb_v_ogrid_ens, num_obs*DEF_DA_ENS, MPI_DOUBLE, p_address_master, p_comm_glb, p_err)
#endif

#endif
      ENDIF


!#############################################################################
! Run data assimilation (for only DA)
!#############################################################################
      has_DA = .false.
      IF (has_obs) THEN

         ! for single point data assimilation
#ifdef SinglePoint
         IF (p_is_worker) THEN

            ! regions info around target patch
            lat_p_n = patchlatr(1)*180/pi + dres
            lat_p_s = patchlatr(1)*180/pi - dres
            lon_p_w = patchlonr(1)*180/pi - dres
            lon_p_e = patchlonr(1)*180/pi + dres

            ! find observations around target patch
            num_obs_p = count(fy3d_lat(:) < lat_p_n .and. fy3d_lat(:) > lat_p_s .and. &
               fy3d_lon(:) > lon_p_w .and. fy3d_lon(:) < lon_p_e .and. &
               fy3d_time(:) - idate_b >= 0 .and. fy3d_time(:) - idate_e <= 0)

            ! perform data assimilation if have observations
            IF (num_obs_p > 0) THEN
               allocate (d_p            (num_obs_p         ))
               allocate (obs_err        (num_obs_p         ))
               allocate (index_p        (num_obs_p         ))
               allocate (fy3d_lat_p     (num_obs_p         ))
               allocate (fy3d_lon_p     (num_obs_p         ))
               allocate (fy3d_tb_h_p    (num_obs_p         ))
               allocate (fy3d_tb_v_p    (num_obs_p         ))
               allocate (pred_tb_h_p_ens(num_obs_p, DEF_DA_ENS))
               allocate (pred_tb_v_p_ens(num_obs_p, DEF_DA_ENS))
               allocate (trans          (DEF_DA_ENS,DEF_DA_ENS))
               allocate (wice_soi_ens   (nl_soil,   DEF_DA_ENS))
               allocate (wice_soi_ens_da(nl_soil,   DEF_DA_ENS))
               allocate (wliq_soi_ens   (nl_soil,   DEF_DA_ENS))
               allocate (wliq_soi_ens_da(nl_soil,   DEF_DA_ENS))

               ! index of observations around target patch
               index_p = (fy3d_lat(:) < lat_p_n .and. fy3d_lat(:) > lat_p_s .and. &
                  fy3d_lon(:) > lon_p_w .and. fy3d_lon(:) < lon_p_e .and. &
                  fy3d_time(:) - idate_b >= 0 .and. fy3d_time(:) - idate_e <= 0)

               ! crop observations around target patch
               fy3d_lat_p  = pack(fy3d_lat , index_p)
               fy3d_lon_p  = pack(fy3d_lon , index_p)
               fy3d_tb_h_p = pack(fy3d_tb_h, index_p)
               fy3d_tb_v_p = pack(fy3d_tb_v, index_p)

               ! predicted observations around target patch
               DO i = 1, num_obs_p
                  pred_tb_h_p_ens(i, :) = pred_tb_h_pset_ens(:, 1) ! only one patch
                  pred_tb_v_p_ens(i, :) = pred_tb_v_pset_ens(:, 1)
               ENDDO

               ! calculate distance between observations and target patch
               d_p = 2*6.3781e3*asin(sqrt(sin((fy3d_lat_p*pi/180 - patchlatr(1))/2.0)**2 + &
                  cos(fy3d_lat_p*pi/180)*cos(patchlatr(1))*sin((fy3d_lon_p*pi/180 - patchlonr(1))/2.0)**2))
               
               ! setting observation error matrix
               obs_err(:) = static_obs_err 

               ! calculate transformation matrix
               CALL letkf(DEF_DA_ENS, num_obs_p, &
                  pred_tb_h_p_ens, fy3d_tb_h_p, obs_err, d_p, loc_r, infl, &
                  trans)

               ! calculate analysis value
               IF (wliq_soisno_ens(1, 1, 1) /= spval) THEN
                  has_DA = .true.

                  ! soil layer
                  wliq_soi_ens = wliq_soisno_ens(1:, :, 1)
                  wice_soi_ens = wice_soisno_ens(1:, :, 1)

                  ! analysis
                  CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS, DEF_DA_ENS, 1.0_8, wliq_soi_ens, &
                     nl_soil, trans, DEF_DA_ENS, 0.0_8, wliq_soi_ens_da, nl_soil)
                  CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS, DEF_DA_ENS, 1.0_8, wice_soi_ens, &
                     nl_soil, trans, DEF_DA_ENS, 0.0_8, wice_soi_ens_da, nl_soil)
                  wliq_soisno_ens(1:, :, 1) = wliq_soi_ens_da
                  wice_soisno_ens(1:, :, 1) = wice_soi_ens_da

                  ! limit the soil liquid and ice water in a reasonable range
                  DO i = 1, nl_soil
                     DO iens = 1, DEF_DA_ENS
                        ! lower bound
                        wliq_soisno_ens(i, iens, 1) = max(0.0d0, wliq_soisno_ens(i, iens, 1))
                        wice_soisno_ens(i, iens, 1) = max(0.0d0, wice_soisno_ens(i, iens, 1))

                        ! upper bound
                        wice_soisno_ens(i, iens, 1) = min(porsl(i, 1)*(dz_soi(i)*denice), wice_soisno_ens(i, iens, 1))
                        eff_porsl = max(0.0d0, porsl(i, 1) - wice_soisno_ens(i, iens, 1)/(dz_soi(i)*denice))
                        wliq_soisno_ens(i, iens, 1) = min(eff_porsl*(dz_soi(i)*denh2o), wliq_soisno_ens(i, iens, 1))
                     ENDDO
                  ENDDO

                  ! move residual water to water table
                  wa_ens(:, 1) = wa_ens(:, 1) - sum(wliq_soisno_ens(1:, :, 1) + wice_soisno_ens(1:, :, 1) - wliq_soi_ens - wice_soi_ens, dim=1)

                  ! update volumetric water content for diagnostic
                  DO iens = 1, DEF_DA_ENS
                     h2osoi_ens(:, iens, 1) = wliq_soisno_ens(1:, iens, 1)/(dz_soi(:)*denh2o) + wice_soisno_ens(1:, iens, 1)/(dz_soi(:)*denice)
                     h2osoi_ens(:, iens, 1) = min(1.0d0, h2osoi_ens(:, iens, 1))
                     h2osoi_ens(:, iens, 1) = max(0.0d0, h2osoi_ens(:, iens, 1))
                  ENDDO
               ENDIF
            ENDIF

            ! deallocate memory (cuz dimensions changes with patch)
            IF (allocated(index_p         )) deallocate (index_p          )
            IF (allocated(fy3d_lat_p      )) deallocate (fy3d_lat_p       )
            IF (allocated(fy3d_lon_p      )) deallocate (fy3d_lon_p       )
            IF (allocated(fy3d_tb_h_p     )) deallocate (fy3d_tb_h_p      )
            IF (allocated(fy3d_tb_v_p     )) deallocate (fy3d_tb_v_p      )
            IF (allocated(pred_tb_h_p_ens )) deallocate (pred_tb_h_p_ens  )
            IF (allocated(pred_tb_v_p_ens )) deallocate (pred_tb_v_p_ens  )
            IF (allocated(d_p             )) deallocate (d_p              )
            IF (allocated(obs_err         )) deallocate (obs_err          )
            IF (allocated(trans           )) deallocate (trans            )
            IF (allocated(wice_soi_ens    )) deallocate (wice_soi_ens     )
            IF (allocated(wice_soi_ens_da )) deallocate (wice_soi_ens_da  )
            IF (allocated(wliq_soi_ens    )) deallocate (wliq_soi_ens     )
            IF (allocated(wliq_soi_ens_da )) deallocate (wliq_soi_ens_da  )
         ENDIF
#else
         ! for grid data assimilation
         IF (p_is_worker) THEN
            DO np = 1, numpatch

               ! regions info around target patch
               lat_p_n = patchlatr(np)*180/pi + dres
               lat_p_s = patchlatr(np)*180/pi - dres
               lon_p_w = patchlonr(np)*180/pi - dres
               lon_p_e = patchlonr(np)*180/pi + dres

               ! find observations around each patch
               num_obs_p = count(fy3d_lat(:) < lat_p_n .and. fy3d_lat(:) > lat_p_s .and. &
                  fy3d_lon(:) > lon_p_w .and. fy3d_lon(:) < lon_p_e .and. &
                  fy3d_time(:) - idate_b >= 0 .and. fy3d_time(:) - idate_e <= 0 .and. &
                  pred_tb_h_ogrid_ens(:, 1) > 0)

               ! perform data assimilation if have observations
               IF (num_obs_p > 0) THEN
                  allocate (index_p        (num_obs_p            ))
                  allocate (fy3d_lat_p     (num_obs_p            ))
                  allocate (fy3d_lon_p     (num_obs_p            ))
                  allocate (fy3d_tb_h_p    (num_obs_p            ))
                  allocate (fy3d_tb_v_p    (num_obs_p            ))
                  allocate (pred_tb_h_p_ens(num_obs_p, DEF_DA_ENS))
                  allocate (pred_tb_v_p_ens(num_obs_p, DEF_DA_ENS))
                  allocate (d_p            (num_obs_p            ))
                  allocate (obs_err        (num_obs_p            ))
                  allocate (trans          (DEF_DA_ENS,DEF_DA_ENS))
                  allocate (wice_soi_ens   (nl_soil,   DEF_DA_ENS))
                  allocate (wice_soi_ens_da(nl_soil,   DEF_DA_ENS))
                  allocate (wliq_soi_ens   (nl_soil,   DEF_DA_ENS))
                  allocate (wliq_soi_ens_da(nl_soil,   DEF_DA_ENS))

                  ! index of observations around target patch
                  index_p = (fy3d_lat(:) < lat_p_n .and. fy3d_lat(:) > lat_p_s .and. &
                     fy3d_lon(:) > lon_p_w .and. fy3d_lon(:) < lon_p_e .and. &
                     fy3d_time(:) - idate_b >= 0 .and. fy3d_time(:) - idate_e <= 0 .and. &
                     pred_tb_h_ogrid_ens(:, 1) > 0)

                  ! crop observations around target patch
                  fy3d_lat_p  = pack(fy3d_lat , index_p)
                  fy3d_lon_p  = pack(fy3d_lon , index_p)
                  fy3d_tb_h_p = pack(fy3d_tb_h, index_p)
                  fy3d_tb_v_p = pack(fy3d_tb_v, index_p)

                  ! predicted observations around target patch
                  DO i = 1, DEF_DA_ENS
                     pred_tb_h_p_ens(:, i) = pack(pred_tb_h_ogrid_ens(:, i), index_p)
                     pred_tb_v_p_ens(:, i) = pack(pred_tb_v_ogrid_ens(:, i), index_p)
                  ENDDO

                  ! calculate distance between observations and target patch
                  d_p = 2*6.3781e3*asin(sqrt(sin((fy3d_lat_p*pi/180 - patchlatr(np))/2.0)**2 + &
                     cos(fy3d_lat_p*pi/180)*cos(patchlatr(np))*sin((fy3d_lon_p*pi/180 - patchlonr(np))/2.0)**2))
                  
                  ! setting observation error matrix
                  obs_err(:) = static_obs_err

                  ! calculate transformation matrix
                  CALL letkf(DEF_DA_ENS, num_obs_p, &
                     pred_tb_h_p_ens, fy3d_tb_h_p, obs_err, d_p, loc_r, infl, &
                     trans)

                  ! calculate analysis value
                  IF (wliq_soisno_ens(1, 1, np) /= spval) THEN
                     has_DA = .TRUE.

                     ! soil layer
                     DO iens = 1, DEF_DA_ENS
                        wliq_soi_ens(:, iens) = wliq_soisno_ens(1:, iens, np)
                        wice_soi_ens(:, iens) = wice_soisno_ens(1:, iens, np)
                     ENDDO

                     ! analysis
                     CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS, DEF_DA_ENS, 1.0_8, wliq_soi_ens, &
                        nl_soil, trans, DEF_DA_ENS, 0.0_8, wliq_soi_ens_da, nl_soil)
                     CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS, DEF_DA_ENS, 1.0_8, wice_soi_ens, &
                        nl_soil, trans, DEF_DA_ENS, 0.0_8, wice_soi_ens_da, nl_soil)

                     DO iens = 1, DEF_DA_ENS
                        wliq_soisno_ens(1:, iens, np) = wliq_soi_ens_da(1:, iens)
                        wice_soisno_ens(1:, iens, np) = wice_soi_ens_da(1:, iens)
                     ENDDO

                     ! limit the soil liquid and ice water in a reasonable range
                     DO i = 1, nl_soil
                        DO iens = 1, DEF_DA_ENS
                           ! lower bound
                           wliq_soisno_ens(i, iens, np) = max(0.0d0, wliq_soisno_ens(i, iens, np))
                           wice_soisno_ens(i, iens, np) = max(0.0d0, wice_soisno_ens(i, iens, np))
                           IF (wliq_soisno_ens(i, iens, np) == 0.0 .and. wice_soisno_ens(i, iens, np) == 0.0) THEN
                              IF (t_soisno_ens(i, iens, np) < -5.0) THEN
                                 wice_soisno_ens(i, iens, np) = 1e-10
                              ELSE
                                 wliq_soisno_ens(i, iens, np) = 1e-10
                              ENDIF
                           ENDIF

                           ! upper bound
                           wice_soisno_ens(i, iens, np) = min(porsl(i, np)*(dz_soi(i)*denice), wice_soisno_ens(i, iens, np))
                           eff_porsl = max(0.0d0, porsl(i, np) - wice_soisno_ens(i, iens, np)/(dz_soi(i)*denice))
                           wliq_soisno_ens(i, iens, np) = min(eff_porsl*(dz_soi(i)*denh2o), wliq_soisno_ens(i, iens, np))
                        ENDDO
                     ENDDO

                     ! move residual water to water table
                     wa_ens(:, np) = wa_ens(:, np) - sum(wliq_soisno_ens(1:, :, np) + wice_soisno_ens(1:, :, np) - wliq_soi_ens - wice_soi_ens, dim=1)

                     ! update volumetric water content for diagnostic
                     DO iens = 1, DEF_DA_ENS
                        h2osoi_ens(:, iens, np) = wliq_soisno_ens(1:, iens, np)/(dz_soi(:)*denh2o) + wice_soisno_ens(1:, iens, np)/(dz_soi(:)*denice)
                        h2osoi_ens(:, iens, np) = min(1.0d0, h2osoi_ens(:, iens, np))
                        h2osoi_ens(:, iens, np) = max(0.0d0, h2osoi_ens(:, iens, np))
                     ENDDO
                  ENDIF
               ENDIF

               ! deallocate memory (cuz dimensions changes with patch)
               IF (allocated(index_p         )) deallocate (index_p          )
               IF (allocated(fy3d_lat_p      )) deallocate (fy3d_lat_p       )
               IF (allocated(fy3d_lon_p      )) deallocate (fy3d_lon_p       )
               IF (allocated(fy3d_tb_h_p     )) deallocate (fy3d_tb_h_p      )
               IF (allocated(fy3d_tb_v_p     )) deallocate (fy3d_tb_v_p      )
               IF (allocated(pred_tb_h_p_ens )) deallocate (pred_tb_h_p_ens  )
               IF (allocated(pred_tb_v_p_ens )) deallocate (pred_tb_v_p_ens  )
               IF (allocated(d_p             )) deallocate (d_p              )
               IF (allocated(obs_err         )) deallocate (obs_err          )
               IF (allocated(trans           )) deallocate (trans            )
               IF (allocated(wice_soi_ens    )) deallocate (wice_soi_ens     )
               IF (allocated(wice_soi_ens_da )) deallocate (wice_soi_ens_da  )
               IF (allocated(wliq_soi_ens    )) deallocate (wliq_soi_ens     )
               IF (allocated(wliq_soi_ens_da )) deallocate (wliq_soi_ens_da  )
            ENDDO
         ENDIF
#endif
      ENDIF


!#############################################################################
! Calculate ensemble brightness temperature after DA for diagnostic
!#############################################################################
      IF (has_obs) THEN
         IF (p_is_worker) THEN
            DO iens = 1, DEF_DA_ENS
               DO np = 1, numpatch
                  CALL forward( &
                     patchtype(np), patchclass(np), dz_sno_ens(:, iens, np), &
                     forc_topo(np), htop(np), &
                     tref_ens(iens, np), t_soisno_ens(:, iens, np), tleaf_ens(iens, np), &
                     wliq_soisno_ens(:, iens, np), wice_soisno_ens(:, iens, np), h2osoi_ens(:, iens, np), &
                     snowdp_ens(iens, np), lai_ens(iens, np), sai_ens(iens, np), &
                     vf_clay(:, np), vf_sand(:, np), BD_all(:, np), porsl(:, np), &
                     fy3d_theta, fy3d_fghz, &
                     t_brt_ens(1, iens, np), t_brt_ens(2, iens, np))
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      ! deallocate memory (cuz dimensions changes with patch)
      IF (allocated(index_p        )) deallocate (index_p        )
      IF (allocated(fy3d_lat_p     )) deallocate (fy3d_lat_p     )
      IF (allocated(fy3d_lon_p     )) deallocate (fy3d_lon_p     )
      IF (allocated(fy3d_tb_h_p    )) deallocate (fy3d_tb_h_p    )
      IF (allocated(fy3d_tb_v_p    )) deallocate (fy3d_tb_v_p    )
      IF (allocated(pred_tb_h_p_ens)) deallocate (pred_tb_h_p_ens)
      IF (allocated(pred_tb_v_p_ens)) deallocate (pred_tb_v_p_ens)
      IF (allocated(d_p            )) deallocate (d_p            )
      IF (allocated(obs_err        )) deallocate (obs_err        )

      IF (allocated(fy3d_time      )) deallocate (fy3d_time      )
      IF (allocated(fy3d_lat       )) deallocate (fy3d_lat       )
      IF (allocated(fy3d_lon       )) deallocate (fy3d_lon       )
      IF (allocated(fy3d_tb_h      )) deallocate (fy3d_tb_h      )
      IF (allocated(fy3d_tb_v      )) deallocate (fy3d_tb_v      )
      IF (allocated(fy3d_ii        )) deallocate (fy3d_ii        )
      IF (allocated(fy3d_jj        )) deallocate (fy3d_jj        )

      IF (allocated(trans           )) deallocate (trans            )
      IF (allocated(wice_soi_ens    )) deallocate (wice_soi_ens     )
      IF (allocated(wice_soi_ens_da )) deallocate (wice_soi_ens_da  )
      IF (allocated(wliq_soi_ens    )) deallocate (wliq_soi_ens     )
      IF (allocated(wliq_soi_ens_da )) deallocate (wliq_soi_ens_da  )

   END SUBROUTINE run_DA_FY3D

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_FY3D()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (allocated(fy3d_time              )) deallocate (fy3d_time              )
      IF (allocated(dt_b                   )) deallocate (dt_b                   )
      IF (allocated(dt_e                   )) deallocate (dt_e                   )
      IF (allocated(fy3d_lat               )) deallocate (fy3d_lat               )
      IF (allocated(fy3d_lon               )) deallocate (fy3d_lon               )
      IF (allocated(fy3d_tb_h              )) deallocate (fy3d_tb_h              )
      IF (allocated(fy3d_tb_v              )) deallocate (fy3d_tb_v              )
      IF (allocated(fy3d_ii                )) deallocate (fy3d_ii                )
      IF (allocated(fy3d_jj                )) deallocate (fy3d_jj                )
      IF (allocated(pred_tb_h_pset_ens     )) deallocate (pred_tb_h_pset_ens     )
      IF (allocated(pred_tb_v_pset_ens     )) deallocate (pred_tb_v_pset_ens     )
      IF (allocated(area_pset              )) deallocate (area_pset              )
      IF (allocated(pred_tb_h_ogrid_ens    )) deallocate (pred_tb_h_ogrid_ens    )
      IF (allocated(pred_tb_v_ogrid_ens    )) deallocate (pred_tb_v_ogrid_ens    )
      IF (allocated(obs_err                )) deallocate (obs_err                )
      IF (allocated(index_p                )) deallocate (index_p                )
      IF (allocated(fy3d_lat_p             )) deallocate (fy3d_lat_p             )
      IF (allocated(fy3d_lon_p             )) deallocate (fy3d_lon_p             )
      IF (allocated(fy3d_tb_h_p            )) deallocate (fy3d_tb_h_p            )
      IF (allocated(fy3d_tb_v_p            )) deallocate (fy3d_tb_v_p            )
      IF (allocated(pred_tb_h_p_ens        )) deallocate (pred_tb_h_p_ens        )
      IF (allocated(pred_tb_v_p_ens        )) deallocate (pred_tb_v_p_ens        )
      IF (allocated(d_p                    )) deallocate (d_p                    )
      IF (allocated(trans                  )) deallocate (trans                  )
      IF (allocated(wice_soi_ens           )) deallocate (wice_soi_ens           )
      IF (allocated(wice_soi_ens_da        )) deallocate (wice_soi_ens_da        )
      IF (allocated(wliq_soi_ens           )) deallocate (wliq_soi_ens           )
      IF (allocated(wliq_soi_ens_da        )) deallocate (wliq_soi_ens_da        )
   END SUBROUTINE end_DA_FY3D

!-----------------------------------------------------------------------------
END MODULE MOD_DA_FY3D
#endif