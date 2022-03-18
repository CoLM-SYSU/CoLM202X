#include <define.h>

module mod_hist

   use precision
   use mod_grid
   use mod_mapping_pset2grid
   USE mod_namelist
   USE GlobalVars, only : spval

#if(defined CaMa_Flood)
   USE YOS_CMF_INPUT,      ONLY: NX, NY
   USE YOS_CMF_MAP,        ONLY: D1LON, D1LAT
   USE MOD_1D_cama_Fluxes
   USE CMF_CALC_DIAG_MOD,  ONLY: CMF_DIAG_AVERAGE, CMF_DIAG_RESET
#endif

   type(grid_type), target :: ghist
   type(mapping_pset2grid_type) :: mp2g_hist

   public :: hist_init
   public :: hist_out
   public :: hist_final
#if(defined CaMa_Flood)
   public :: var_out
#endif

!------------------------------------------------------------------

!>>>>>add by zhongwang wei 
type :: segment_type
   integer :: blk
   integer :: cnt
   integer :: bdsp
   integer :: gdsp
end type segment_type

type :: grid_info_type
   integer :: nlat, nlon
   real(r8), allocatable :: lat_s(:)
   real(r8), allocatable :: lat_n(:)
   real(r8), allocatable :: lon_w(:)
   real(r8), allocatable :: lon_e(:)
   real(r8), allocatable :: lon_c(:) !grid center
   real(r8), allocatable :: lat_c(:) !grid center
end type grid_info_type

TYPE :: block_info_type
   integer :: ndatablk
   integer :: nxseg, nyseg
   type(segment_type), allocatable :: xsegs(:), ysegs(:)
   type(grid_info_type) :: ginfo
CONTAINS
   final :: block_info_free_mem
END TYPE block_info_type

TYPE(block_info_type) :: hist_block_info

integer :: hist_data_id
!<<<<<add by zhongwang wei
!--------------------------------------------------------------------------
contains

   !---------------------------------------
   subroutine hist_init (dir_hist, lon_points, lat_points)

      USE GlobalVars
      use spmd_task
      use mod_grid
      USE mod_landpatch
      use mod_mapping_pset2grid
      use MOD_1D_Acc_Fluxes
#if(defined CaMa_Flood)
      use MOD_1D_Acc_cama_Fluxes
#endif
      implicit none

      character(len=*), intent(in) :: dir_hist
      integer, intent(in) :: lon_points
      integer, intent(in) :: lat_points

      call ghist%define_by_ndims (lon_points, lat_points)
#ifndef CROP
      call mp2g_hist%build (landpatch, ghist)
#else
      call mp2g_hist%build (landpatch, ghist, pctcrop)
#endif

      call allocate_acc_fluxes ()
      call FLUSH_acc_fluxes ()
#if(defined CaMa_Flood)
      call allocate_acc_cama_fluxes ()
      call FLUSH_acc_cama_fluxes ()
#endif

      !>>>>>add by zhongwang wei 
      if (trim(DEF_HIST_mode) == 'one') then
         call set_hist_block_info (ghist, hist_block_info)
         hist_data_id = 1000
      end if
      !<<<<<add by zhongwang wei


   end subroutine hist_init 

   !--------------------------------------
   subroutine hist_final ()

      use MOD_1D_Acc_Fluxes
      use MOD_1D_Acc_cama_Fluxes

      implicit none
      
      call deallocate_acc_fluxes ()
#if(defined CaMa_Flood)
      call deallocate_acc_cama_fluxes()
#endif

   end subroutine hist_final

   !---------------------------------------
   SUBROUTINE hist_out (idate, deltim, itstamp, ptstamp, &
         dir_hist, site)

      !=======================================================================
      ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
      !=======================================================================

      use precision
      use mod_namelist
      use timemanager
      use spmd_task
      use ncio_serial
      use mod_2d_fluxes
      use MOD_1D_Acc_Fluxes
      use mod_block
      use mod_data_type
      use mod_landpatch
      use mod_mapping_pset2grid
      use MOD_2D_Fluxes
      use mod_colm_debug
      use GlobalVars, only : spval
      USE MOD_TimeInvariants, only : patchtype
#if(defined CaMa_Flood)
      USE YOS_CMF_DIAG,       ONLY: D2OUTFLW_AVG
      USE YOS_CMF_PROG,       ONLY: D2RIVSTO,     D2FLDSTO,     D2GDWSTO, &
      & d2damsto !!! added
      USE YOS_CMF_DIAG,       ONLY: D2RIVDPH,     D2FLDDPH,     D2FLDFRC,     D2FLDARE,     D2SFCELV,     D2STORGE, &
      & D2OUTFLW_AVG, D2RIVOUT_AVG, D2FLDOUT_AVG, D2PTHOUT_AVG, D1PTHFLW_AVG, &
      & D2RIVVEL_AVG, D2GDWRTN_AVG, D2RUNOFF_AVG, D2ROFSUB_AVG,               &
      & D2OUTFLW_MAX, D2STORGE_MAX, D2RIVDPH_MAX, &
      & d2daminf_avg   !!! added
#endif
      IMPLICIT NONE

      integer,  INTENT(in) :: idate(3)
      real(r8), INTENT(in) :: deltim
      type(timestamp), intent(in) :: itstamp
      type(timestamp), intent(in) :: ptstamp

      character(LEN=*), intent(in) :: dir_hist
      character(LEN=*), intent(in) :: site

      ! Local variables
      logical :: lwrite
      character(LEN=256) :: file_hist
      integer :: itime_in_file
      integer :: month, day
      integer :: days_month(1:12)
      character(len=10) :: cdate
      character(len=256) :: groupby
      
      type(block_data_real8_2d) :: sumwt
      real(r8), allocatable ::  vectmp(:)  
      logical,  allocatable ::  filter(:)

      
      if (itstamp <= ptstamp) then
         call FLUSH_acc_fluxes ()
         return 
      else

         call accumulate_fluxes 

      end if

      select case (trim(DEF_HIST_FREQ))
      case ('HOURLY')
         lwrite = isendofhour (idate, deltim)
      case ('DAILY')
         lwrite = isendofday(idate, deltim)
      case ('MONTHLY')
         lwrite = isendofmonth(idate, deltim)       
      case ('YEARLY')
         lwrite = isendofyear(idate, deltim)
      end select

      if (lwrite)  then

         call julian2monthday(idate(1), idate(2), month, day)

         days_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)
         if (isleapyear(idate(1))) days_month(2) = 29

         groupby = DEF_HIST_groupby

         if ( trim(groupby) == 'YEAR' ) then
            write(cdate,'(i4.4)') idate(1)
         elseif ( trim(groupby) == 'MONTH' ) then
            write(cdate,'(i4.4,"-",i2.2)') idate(1), month
         elseif ( trim(groupby) == 'DAY' ) then
            write(cdate,'(i4.4,"-",i2.2,"-",i2.2)') idate(1), month, day
         end if

         file_hist = trim(dir_hist) // '/' // trim(site) //'_hist_'//trim(cdate)//'.nc'

         call hist_write_time (file_hist, 'time', ghist, idate, itime_in_file)  

         if (p_is_worker) then
            if (numpatch > 0) then
               allocate (filter (numpatch))
               allocate (vectmp (numpatch))
            end if
         end if

         if (p_is_io) then
            call allocate_block_data (ghist, sumwt)
         end if

         ! ---------------------------------------------------
         ! Meteorological forcing
         ! ---------------------------------------------------
         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype < 99
               vectmp(:) = 1.
            end if
         end if

         call mp2g_hist%map (vectmp, sumwt, spv = spval, msk = filter)

         ! wind in eastward direction [m/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_us, &
            a_us, f_xy_us, file_hist, 'f_xy_us', itime_in_file, sumwt, filter, &
            'wind in eastward direction', 'm/s')

         ! wind in northward direction [m/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_vs, &
            a_vs, f_xy_vs, file_hist, 'f_xy_vs', itime_in_file, sumwt, filter, &
            'wind in northward direction','m/s')

         ! temperature at reference height [kelvin]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_t, &
            a_t, f_xy_t, file_hist, 'f_xy_t', itime_in_file, sumwt, filter, &
            'temperature at reference height','kelvin')

         ! specific humidity at reference height [kg/kg]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_q, &
            a_q, f_xy_q, file_hist, 'f_xy_q', itime_in_file, sumwt, filter, &
            'specific humidity at reference height','kg/kg')

         ! convective precipitation [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_prc, &
            a_prc, f_xy_prc, file_hist, 'f_xy_prc', itime_in_file, sumwt, filter, &
            'convective precipitation','mm/s')

         ! large scale precipitation [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_prl, &
            a_prl, f_xy_prl, file_hist, 'f_xy_prl', itime_in_file, sumwt, filter, &
            'large scale precipitation','mm/s')

         ! atmospheric pressure at the surface [pa]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_pbot, &
            a_pbot, f_xy_pbot, file_hist, 'f_xy_pbot', itime_in_file, sumwt, filter, &
            'atmospheric pressure at the surface','pa')

         ! atmospheric infrared (longwave) radiation [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_frl, &
            a_frl, f_xy_frl, file_hist, 'f_xy_frl', itime_in_file, sumwt, filter, &
            'atmospheric infrared (longwave) radiation','W/m2')
         
         ! downward solar radiation at surface [W/m2]       
         call flux_map_and_write_2d ( DEF_hist_vars%xy_solarin, &
            a_solarin, f_xy_solarin, file_hist, 'f_xy_solarin', itime_in_file, sumwt, filter, &
            'downward solar radiation at surface','W/m2')

         ! rain [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_rain, &
            a_rain, f_xy_rain, file_hist, 'f_xy_rain', itime_in_file, sumwt, filter, &
            'rain','mm/s')
         
         ! snow [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_snow, &
            a_snow, f_xy_snow, file_hist, 'f_xy_snow', itime_in_file, sumwt, filter, &
            'snow','mm/s')

         ! ------------------------------------------------------------------------------------------
         ! Mapping the fluxes and state variables at patch [numpatch] to grid 
         ! ------------------------------------------------------------------------------------------
         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype < 99
               vectmp (:) = 1.
            end if
         end if

         call mp2g_hist%map (vectmp, sumwt, spv = spval, msk = filter)

         ! wind stress: E-W [kg/m/s2]                                   
         call flux_map_and_write_2d ( DEF_hist_vars%taux, &
            a_taux, f_taux, file_hist, 'f_taux', itime_in_file, sumwt, filter, &
            'wind stress: E-W','kg/m/s2')

         ! wind stress: N-S [kg/m/s2]
         call flux_map_and_write_2d ( DEF_hist_vars%tauy, &
            a_tauy, f_tauy, file_hist, 'f_tauy', itime_in_file, sumwt, filter, &
            'wind stress: N-S','kg/m/s2')

         ! sensible heat from canopy height to atmosphere [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%fsena, &
            a_fsena, f_fsena, file_hist, 'f_fsena', itime_in_file, sumwt, filter, &
            'sensible heat from canopy height to atmosphere','W/m2')

         ! latent heat flux from canopy height to atmosphere [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%lfevpa, &
            a_lfevpa, f_lfevpa, file_hist, 'f_lfevpa', itime_in_file, sumwt, filter, &
            'latent heat flux from canopy height to atmosphere','W/m2')

         ! evapotranspiration from canopy to atmosphere [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%fevpa, &
            a_fevpa, f_fevpa, file_hist, 'f_fevpa', itime_in_file, sumwt, filter, &
            'evapotranspiration from canopy height to atmosphere','mm/s')

         ! sensible heat from leaves [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%fsenl, &
            a_fsenl, f_fsenl, file_hist, 'f_fsenl', itime_in_file, sumwt, filter, &
            'sensible heat from leaves','W/m2')

         ! evaporation+transpiration from leaves [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%fevpl, &
            a_fevpl, f_fevpl, file_hist, 'f_fevpl', itime_in_file, sumwt, filter, &
            'evaporation+transpiration from leaves','mm/s')

         ! transpiration rate [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%etr, &
            a_etr, f_etr, file_hist, 'f_etr', itime_in_file, sumwt, filter, &
            'transpiration rate','mm/s')

         ! sensible heat flux from ground [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%fseng, &
            a_fseng, f_fseng, file_hist, 'f_fseng', itime_in_file, sumwt, filter, &
            'sensible heat flux from ground','W/m2')

         ! evaporation heat flux from ground [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%fevpg, &
            a_fevpg, f_fevpg, file_hist, 'f_fevpg', itime_in_file, sumwt, filter, &
            'evaporation heat flux from ground','mm/s')

         ! ground heat flux [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%fgrnd, &
            a_fgrnd, f_fgrnd, file_hist, 'f_fgrnd', itime_in_file, sumwt, filter, &
            'ground heat flux','W/m2')

         ! solar absorbed by sunlit canopy [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%sabvsun, &
            a_sabvsun, f_sabvsun, file_hist, 'f_sabvsun', itime_in_file, sumwt, filter, &
            'solar absorbed by sunlit canopy','W/m2')

         ! solar absorbed by shaded [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%sabvsha, &
            a_sabvsha, f_sabvsha, file_hist, 'f_sabvsha', itime_in_file, sumwt, filter, &
            'solar absorbed by shaded','W/m2')

         ! solar absorbed by ground  [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%sabg, &
            a_sabg, f_sabg, file_hist, 'f_sabg', itime_in_file, sumwt, filter, &
            'solar absorbed by ground','W/m2')

         ! outgoing long-wave radiation from ground+canopy [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%olrg, &
            a_olrg, f_olrg, file_hist, 'f_olrg', itime_in_file, sumwt, filter, &
            'outgoing long-wave radiation from ground+canopy','W/m2')

         ! net radiation [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%rnet, &
            a_rnet, f_rnet, file_hist, 'f_rnet', itime_in_file, sumwt, filter, &
            'net radiation','W/m2')

         ! the error of water banace [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xerr, &
            a_xerr, f_xerr, file_hist, 'f_xerr', itime_in_file, sumwt, filter, &
            'the error of water banace','mm/s')

         ! the error of energy balance [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%zerr, &
            a_zerr, f_zerr, file_hist, 'f_zerr', itime_in_file, sumwt, filter, &
            'the error of energy balance','W/m2')

         ! surface runoff [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%rsur, &
            a_rsur, f_rsur, file_hist, 'f_rsur', itime_in_file, sumwt, filter, &
            'surface runoff','mm/s')

         ! total runoff [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%rnof, &
            a_rnof, f_rnof, file_hist, 'f_rnof', itime_in_file, sumwt, filter, &
            'total runoff','mm/s')

         ! interception [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%qintr, &
            a_qintr, f_qintr, file_hist, 'f_qintr', itime_in_file, sumwt, filter, &
            'interception','mm/s')

         ! inflitraton [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%qinfl, &
            a_qinfl, f_qinfl, file_hist, 'f_qinfl', itime_in_file, sumwt, filter, &
            'f_qinfl','mm/s')

         ! throughfall [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%qdrip, &
            a_qdrip, f_qdrip, file_hist, 'f_qdrip', itime_in_file, sumwt, filter, &
            'total throughfall','mm/s')

         ! total water storage [mm]
         call flux_map_and_write_2d ( DEF_hist_vars%wat, &
            a_wat, f_wat, file_hist, 'f_wat', itime_in_file, sumwt, filter, &
            'total water storage','mm')

         ! canopy assimilation rate [mol m-2 s-1]
         call flux_map_and_write_2d ( DEF_hist_vars%assim, &
            a_assim, f_assim, file_hist, 'f_assim', itime_in_file, sumwt, filter, &
            'canopy assimilation rate','umol m-2 s-1')

         ! respiration (plant+soil) [mol m-2 s-1]
         call flux_map_and_write_2d ( DEF_hist_vars%respc, &
            a_respc, f_respc, file_hist, 'f_respc', itime_in_file, sumwt, filter, &
            'respiration (plant+soil)','mol m-2 s-1')

         ! groundwater recharge rate [mm/s]                            
         call flux_map_and_write_2d ( DEF_hist_vars%qcharge, &
            a_qcharge, f_qcharge, file_hist, 'f_qcharge', itime_in_file, sumwt, filter, &
            'groundwater recharge rate','mm/s')

         ! ground surface temperature [K]                        
         call flux_map_and_write_2d ( DEF_hist_vars%t_grnd, &
            a_t_grnd, f_t_grnd, file_hist, 'f_t_grnd', itime_in_file, sumwt, filter, &
            'ground surface temperature','K')

         ! leaf temperature [K]
         call flux_map_and_write_2d ( DEF_hist_vars%tleaf, &
            a_tleaf, f_tleaf, file_hist, 'f_tleaf', itime_in_file, sumwt, filter, &
            'leaf temperature','K')

         ! depth of water on foliage [mm]
         call flux_map_and_write_2d ( DEF_hist_vars%ldew, &
            a_ldew, f_ldew, file_hist, 'f_ldew', itime_in_file, sumwt, filter, &
            'depth of water on foliage','mm')

         ! snow cover, water equivalent [mm]
         call flux_map_and_write_2d ( DEF_hist_vars%scv, &
            a_scv, f_scv, file_hist, 'f_scv', itime_in_file, sumwt, filter, &
            'snow cover, water equivalent','mm')

         ! snow depth [meter]
         call flux_map_and_write_2d ( DEF_hist_vars%snowdp, &
            a_snowdp, f_snowdp, file_hist, 'f_snowdp', itime_in_file, sumwt, filter, &
            'snow depth','meter')

         ! fraction of snow cover on ground
         call flux_map_and_write_2d ( DEF_hist_vars%fsno, &
            a_fsno, f_fsno, file_hist, 'f_fsno', itime_in_file, sumwt, filter, &
            'fraction of snow cover on ground','-')

         ! fraction of veg cover, excluding snow-covered veg [-]
         call flux_map_and_write_2d ( DEF_hist_vars%sigf, &
            a_sigf, f_sigf, file_hist, 'f_sigf', itime_in_file, sumwt, filter, &
            'fraction of veg cover, excluding snow-covered veg','-')

         ! leaf greenness
         call flux_map_and_write_2d ( DEF_hist_vars%green, &
            a_green, f_green, file_hist, 'f_green', itime_in_file, sumwt, filter, &
            'leaf greenness','-')

         ! leaf area index
         call flux_map_and_write_2d ( DEF_hist_vars%lai, &
            a_lai, f_lai, file_hist, 'f_lai', itime_in_file, sumwt, filter, &
            'leaf area index','m2/m2')

         ! leaf area index
         call flux_map_and_write_2d ( DEF_hist_vars%laisun, &
            a_laisun, f_laisun, file_hist, 'f_laisun', itime_in_file, sumwt, filter, &
            'sunlit leaf area index','m2/m2')

         ! leaf area index
         call flux_map_and_write_2d ( DEF_hist_vars%laisha, &
            a_laisha, f_laisha, file_hist, 'f_laisha', itime_in_file, sumwt, filter, &
            'shaded leaf area index','m2/m2')

         ! stem area index
         call flux_map_and_write_2d ( DEF_hist_vars%sai, &
            a_sai, f_sai, file_hist, 'f_sai', itime_in_file, sumwt, filter, &
            'stem area index','m2/m2')

         ! averaged albedo [visible, direct; direct, diffuse] 
         call flux_map_and_write_4d ( DEF_hist_vars%alb, &
            a_alb, f_alb, file_hist, 'f_alb', 'band', 'wetdry', itime_in_file, sumwt, filter, &
            'averaged albedo direct','%')

         ! averaged bulk surface emissivity
         call flux_map_and_write_2d ( DEF_hist_vars%emis, &
            a_emis, f_emis, file_hist, 'f_emis', itime_in_file, sumwt, filter, &
            'averaged bulk surface emissivity','-')

         ! effective roughness [m]
         call flux_map_and_write_2d ( DEF_hist_vars%z0m, &
            a_z0m, f_z0m, file_hist, 'f_z0m', itime_in_file, sumwt, filter, &
            'effective roughness','m')

         ! radiative temperature of surface [K]
         call flux_map_and_write_2d ( DEF_hist_vars%trad, &
            a_trad, f_trad, file_hist, 'f_trad', itime_in_file, sumwt, filter, &
            'radiative temperature of surface','kelvin')

         ! 2 m height air temperature [kelvin]
         call flux_map_and_write_2d ( DEF_hist_vars%tref, &
            a_tref, f_tref, file_hist, 'f_tref', itime_in_file, sumwt, filter, & 
            '2 m height air temperature','kelvin')

         ! 2 m height air specific humidity [kg/kg]
         call flux_map_and_write_2d ( DEF_hist_vars%qref, &
            a_qref, f_qref, file_hist, 'f_qref', itime_in_file, sumwt, filter, &
            '2 m height air specific humidity','kg/kg')

         ! --------------------------------------------------------------------
         ! Temperature and water (excluding land water bodies and ocean patches)
         ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3; land water bodies => 4; ocean => 99]
         ! --------------------------------------------------------------------

         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype <= 3
               vectmp (:) = 1.
            end if
         end if

         call mp2g_hist%map (vectmp, sumwt, spv = spval, msk = filter)

         ! soil temperature [K]
         call flux_map_and_write_3d ( DEF_hist_vars%t_soisno, &
            a_t_soisno, f_t_soisno, file_hist, 'f_t_soisno', 'soilsnow', itime_in_file, sumwt, filter, &
            'soil temperature','K')

         ! liquid water in soil layers [kg/m2]
         call flux_map_and_write_3d ( DEF_hist_vars%wliq_soisno, &
            a_wliq_soisno, f_wliq_soisno, file_hist, 'f_wliq_soisno', 'soilsnow', &
            itime_in_file, sumwt, filter,'liquid water in soil layers','kg/m2')

         ! ice lens in soil layers [kg/m2]
         call flux_map_and_write_3d ( DEF_hist_vars%wice_soisno, &
            a_wice_soisno, f_wice_soisno, file_hist, 'f_wice_soisno', 'soilsnow', &
            itime_in_file, sumwt, filter,'ice lens in soil layers','kg/m2')

         ! --------------------------------------------------------------------
         ! additial diagnostic variables for output (vegetated land only <=2)
         ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3; land water bodies => 4; ocean => 99]
         ! --------------------------------------------------------------------

         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype <= 2
               vectmp(:) = 1.
            end if
         end if

         call mp2g_hist%map (vectmp, sumwt, spv = spval, msk = filter)
         
         ! volumetric soil water in layers [m3/m3]
         call flux_map_and_write_3d ( DEF_hist_vars%h2osoi, &
            a_h2osoi, f_h2osoi, file_hist, 'f_h2osoi', 'soil', itime_in_file, sumwt, filter, &
            'volumetric water in soil layers','m3/m3')

         ! fraction of root water uptake from each soil layer, all layers add to 1, when PHS is not defined
         ! water exchange between soil layers and root. Positive: soil->root [mm h2o/s], when PHS is defined
         call flux_map_and_write_3d ( DEF_hist_vars%rootr, &
            a_rootr, f_rootr, file_hist, 'f_rootr', 'soil', itime_in_file, sumwt, filter)

#ifdef PLANT_HYDRAULIC_STRESS
         ! vegetation water potential [mm]
         call flux_map_and_write_3d ( DEF_hist_vars%vegwp, &
            a_vegwp, f_vegwp, file_hist, 'f_vegwp', 'vegnodes', itime_in_file, sumwt, filter)
#endif

         ! water depth [m]
         call flux_map_and_write_2d ( DEF_hist_vars%zwt, &
            a_zwt, f_zwt, file_hist, 'f_zwt', itime_in_file, sumwt, filter, &
            'the depth to water table','m')

         ! water storage in aquifer [mm]
         call flux_map_and_write_2d ( DEF_hist_vars%wa, &
            a_wa, f_wa, file_hist, 'f_wa', itime_in_file, sumwt, filter, &
            'water storage in aquifer','mm')

         ! -----------------------------------------------
         ! Land water bodies' ice fraction and temperature
         ! -----------------------------------------------

         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype == 4
               vectmp(:) = 1.
            end if
         end if

         call mp2g_hist%map (vectmp, sumwt, spv = spval, msk = filter)
         
         ! lake temperature [K]
         call flux_map_and_write_3d ( DEF_hist_vars%t_lake, &
            a_t_lake, f_t_lake, file_hist, 'f_t_lake', 'lake', itime_in_file, sumwt, filter, &
            'lake temperature','K')

         ! lake ice fraction cover [0-1]
         call flux_map_and_write_3d ( DEF_hist_vars%lake_icefrac, &
            a_lake_icefrac, f_lake_icefrac, file_hist, 'f_lake_icefrac', &
            'lake', itime_in_file, sumwt, filter,'lake ice fraction cover','0-1')

         ! --------------------------------
         ! Retrieve through averaged fluxes
         ! --------------------------------
         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype < 99
               vectmp(:) = 1.
            end if
         end if

         call mp2g_hist%map (vectmp, sumwt, spv = spval, msk = filter)
      
         ! u* in similarity theory [m/s]
         call flux_map_and_write_2d ( DEF_hist_vars%ustar, &
            a_ustar, f_ustar, file_hist, 'f_ustar', itime_in_file, sumwt, filter, &
            'u* in similarity theory','m/s')

         ! t* in similarity theory [kg/kg]
         call flux_map_and_write_2d ( DEF_hist_vars%tstar, &
            a_tstar, f_tstar, file_hist, 'f_tstar', itime_in_file, sumwt, filter, &
            't* in similarity theory','kg/kg')

         ! q* in similarity theory [kg/kg]
         call flux_map_and_write_2d ( DEF_hist_vars%qstar, &
            a_qstar, f_qstar, file_hist, 'f_qstar', itime_in_file, sumwt, filter, &
            'q* in similarity theory')

         ! dimensionless height (z/L) used in Monin-Obukhov theory
         call flux_map_and_write_2d ( DEF_hist_vars%zol, &
            a_zol, f_zol, file_hist, 'f_zol', itime_in_file, sumwt, filter, &
            'dimensionless height (z/L) used in Monin-Obukhov theory','-')

         ! bulk Richardson number in surface layer
         call flux_map_and_write_2d ( DEF_hist_vars%rib, &
            a_rib, f_rib, file_hist, 'f_rib', itime_in_file, sumwt, filter, &
            'bulk Richardson number in surface layer','-')

         ! integral of profile function for momentum
         call flux_map_and_write_2d ( DEF_hist_vars%fm, &
            a_fm, f_fm, file_hist, 'f_fm', itime_in_file, sumwt, filter, &
            'integral of profile function for momentum','-')

         ! integral of profile function for heat
         call flux_map_and_write_2d ( DEF_hist_vars%fh, &
            a_fh, f_fh, file_hist, 'f_fh', itime_in_file, sumwt, filter, &
            'integral of profile function for heat','-')

         ! integral of profile function for moisture
         call flux_map_and_write_2d ( DEF_hist_vars%fq, &
            a_fq, f_fq, file_hist, 'f_fq', itime_in_file, sumwt, filter, &
            'integral of profile function for moisture','-')

         ! 10m u-velocity [m/s]
         call flux_map_and_write_2d ( DEF_hist_vars%us10m, &
            a_us10m, f_us10m, file_hist, 'f_us10m', itime_in_file, sumwt, filter, &
            '10m u-velocity','m/s')

         ! 10m v-velocity [m/s]
         call flux_map_and_write_2d ( DEF_hist_vars%vs10m, &
            a_vs10m, f_vs10m, file_hist, 'f_vs10m', itime_in_file, sumwt, filter, &
            '10m v-velocity','m/s')

         ! integral of profile function for momentum at 10m [-]
         call flux_map_and_write_2d ( DEF_hist_vars%fm10m, &
            a_fm10m, f_fm10m, file_hist, 'f_fm10m', itime_in_file, sumwt, filter, &
            'integral of profile function for momentum at 10m','-')
                 
         ! total reflected solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%sr, &
            a_sr, f_sr, file_hist, 'f_sr', itime_in_file, sumwt, filter, &
            'reflected solar radiation at surface [W/m2]','W/m2')

         ! incident direct beam vis solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%solvd, &
            a_solvd, f_solvd, file_hist, 'f_solvd', itime_in_file, sumwt, filter, &
            'incident direct beam vis solar radiation (W/m2)','W/m2')

         ! incident diffuse beam vis solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%solvi, &
            a_solvi, f_solvi, file_hist, 'f_solvi', itime_in_file, sumwt, filter, &
            'incident diffuse beam vis solar radiation (W/m2)','W/m2')

         ! incident direct beam nir solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%solnd, &
            a_solnd, f_solnd, file_hist, 'f_solnd', itime_in_file, sumwt, filter, &
            'incident direct beam nir solar radiation (W/m2)','W/m2')

         ! incident diffuse beam nir solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%solni, &
            a_solni, f_solni, file_hist, 'f_solni', itime_in_file, sumwt, filter, &
            'incident diffuse beam nir solar radiation (W/m2)','W/m2')

         ! reflected direct beam vis solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%srvd, &
            a_srvd, f_srvd, file_hist, 'f_srvd', itime_in_file, sumwt, filter, &
            'reflected direct beam vis solar radiation (W/m2)','W/m2')

         ! reflected diffuse beam vis solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%srvi, &
            a_srvi, f_srvi, file_hist, 'f_srvi', itime_in_file, sumwt, filter, &
            'reflected diffuse beam vis solar radiation (W/m2)','W/m2')

         ! reflected direct beam nir solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%srnd, &
            a_srnd, f_srnd, file_hist, 'f_srnd', itime_in_file, sumwt, filter, &
            'reflected direct beam nir solar radiation (W/m2)','W/m2')

         ! reflected diffuse beam nir solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%srni, &
            a_srni, f_srni, file_hist, 'f_srni', itime_in_file, sumwt, filter, &
            'reflected diffuse beam nir solar radiation (W/m2)','W/m2')

         ! local noon fluxes 
         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = nac_ln > 0
               vectmp(:) = 1.
            end if
         end if

         call mp2g_hist%map (vectmp, sumwt, msk = filter)

         ! incident direct beam vis solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%solvdln, &
            a_solvdln, f_solvdln, file_hist, 'f_solvdln', itime_in_file, sumwt, filter, &
            'incident direct beam vis solar radiation at local noon(W/m2)','W/m2')

         ! incident diffuse beam vis solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%solviln, &
            a_solviln, f_solviln, file_hist, 'f_solviln', itime_in_file, sumwt, filter, &
            'incident diffuse beam vis solar radiation at local noon(W/m2)','W/m2')

         ! incident direct beam nir solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%solndln, &
            a_solndln, f_solndln, file_hist, 'f_solndln', itime_in_file, sumwt, filter, &
            'incident direct beam nir solar radiation at local noon(W/m2)','W/m2')

         ! incident diffuse beam nir solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%solniln, &
            a_solniln, f_solniln, file_hist, 'f_solniln', itime_in_file, sumwt, filter, &
            'incident diffuse beam nir solar radiation at local noon(W/m2)','W/m2')

         ! reflected direct beam vis solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%srvdln, &
            a_srvdln, f_srvdln, file_hist, 'f_srvdln', itime_in_file, sumwt, filter, &
            'reflected direct beam vis solar radiation at local noon(W/m2)','W/m2')

         ! reflected diffuse beam vis solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%srviln, &
            a_srviln, f_srviln, file_hist, 'f_srviln', itime_in_file, sumwt, filter, &
            'reflected diffuse beam vis solar radiation at local noon(W/m2)','W/m2')

         ! reflected direct beam nir solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%srndln, &
            a_srndln, f_srndln, file_hist, 'f_srndln', itime_in_file, sumwt, filter, &
            'reflected direct beam nir solar radiation at local noon(W/m2)','W/m2')

         ! reflected diffuse beam nir solar radiation at local noon(W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%srniln, &
            a_srniln, f_srniln, file_hist, 'f_srniln', itime_in_file, sumwt, filter, &
            'reflected diffuse beam nir solar radiation at local noon(W/m2)','W/m2')

#if(defined CaMa_Flood)
if (p_is_master) then
         !*** average variable
         CALL CMF_DIAG_AVERAGE
         !*** write output data
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%rivout,D2RIVOUT_AVG, file_hist, 'rivout', itime_in_file,'river discharge','m3/s')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%rivsto,D2RIVSTO, file_hist, 'rivsto', itime_in_file,'river storage','m3')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%rivdph,D2RIVDPH, file_hist, 'rivdph', itime_in_file,'river depth','m')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%rivvel,D2RIVVEL_AVG, file_hist, 'rivvel', itime_in_file,'river velocity','m/s')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldout,D2FLDOUT_AVG, file_hist, 'fldout', itime_in_file,'floodplain discharge','m3/s')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldsto, D2FLDSTO, file_hist, 'fldsto', itime_in_file,'floodplain storage','m3')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%flddph, D2FLDDPH, file_hist, 'flddph', itime_in_file,'floodplain depth','m')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldfrc, D2FLDFRC, file_hist, 'fldfrc', itime_in_file,'flooded fraction','0-1')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldare, D2FLDARE, file_hist, 'fldare', itime_in_file,'flooded area','m2')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%sfcelv, D2SFCELV, file_hist, 'sfcelv', itime_in_file,'water surface elevation','m')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%totout,D2OUTFLW_AVG, file_hist, 'totout', itime_in_file,'discharge (river+floodplain)','m3/s')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%outflw,D2OUTFLW_AVG, file_hist, 'outflw', itime_in_file,'discharge (river+floodplain)','m3/s')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%totsto,D2STORGE, file_hist, 'totsto', itime_in_file,'total storage (river+floodplain)','m3')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%storge,D2STORGE, file_hist, 'storge', itime_in_file,'total storage (river+floodplain)','m3')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%pthout,D2PTHOUT_AVG, file_hist, 'pthout', itime_in_file,'net bifurcation discharge','m3/s')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%gdwsto,D2GDWSTO, file_hist, 'gdwsto', itime_in_file,'ground water storage','m3')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%gwsto,D2GDWSTO, file_hist, 'gwsto', itime_in_file,'ground water storage','m3')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%gwout,D2GDWRTN_AVG, file_hist, 'gwout', itime_in_file,'ground water discharge','m3/s')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%runoff,D2RUNOFF_AVG, file_hist, 'runoff', itime_in_file,'Surface runoff','m3/s')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%runoffsub,D2ROFSUB_AVG  , file_hist, 'runoffsub', itime_in_file,'sub-surface runoff','m3/s')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%maxsto,D2STORGE_MAX, file_hist, 'maxsto', itime_in_file,'daily maximum storage','m3')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%maxflw,D2OUTFLW_MAX, file_hist, 'maxflw', itime_in_file,'daily maximum discharge','m3/s')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%maxdph,D2RIVDPH_MAX, file_hist, 'maxdph', itime_in_file,'daily maximum river depth','m')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%damsto,d2damsto, file_hist, 'damsto', itime_in_file,'reservoir storage','m3')
         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%daminf,d2daminf_avg, file_hist, 'daminf', itime_in_file,'reservoir inflow','m3/s')
         !*** reset variable
         CALL CMF_DIAG_RESET
endif
#endif

         if (allocated(filter)) deallocate(filter)
         if (allocated(vectmp)) deallocate(vectmp)

         call FLUSH_acc_fluxes ()

      end if

   END SUBROUTINE hist_out
   
   ! -------
#if(defined CaMa_Flood)
   subroutine flux_map_and_write_2d_cama (is_hist,var_in, file_hist, varname, itime_in_file,longname,units)
   USE YOS_CMF_MAP,        ONLY: NSEQALL
   USE PARKIND1,            ONLY: JPRM
   USE CMF_UTILS_MOD,           ONLY: VEC2MAPD
   use ncio_serial

   IMPLICIT NONE
   logical, intent(in) :: is_hist
   real(r8), INTENT(in) ::  var_in (NSEQALL, 1) 
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character (len=*), intent(in),optional :: longname
   character (len=*), intent(in),optional :: units

   REAL(r8)      :: R2OUT(NX,NY)
   integer :: compress

   if (.not. is_hist) return
   CALL VEC2MAPD(var_in,R2OUT)
   compress = DEF_HIST_COMPRESS_LEVEL 
   call ncio_write_serial_time (file_hist, varname, itime_in_file, R2OUT, 'lon_cama', 'lat_cama', 'time',compress,longname,units)
end subroutine flux_map_and_write_2d_cama


   SUBROUTINE var_out (runoff_2d)

   !=======================================================================
   ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
   !=======================================================================

   use precision
   use mod_namelist
   use timemanager
   use spmd_task
   use mod_2d_cama_fluxes
   use MOD_1D_cama_Fluxes
   use mod_block
   use mod_data_type
   use mod_landpatch
   use mod_mapping_pset2grid
   use mod_colm_debug
   USE MOD_TimeInvariants, only : patchtype
   USE MOD_1D_Acc_cama_Fluxes

   !use GlobalVars, only : spval
   IMPLICIT NONE
   real(r8), INTENT(out) ::  runoff_2d (DEF_nlon_hist, DEF_nlat_hist)

   type(block_data_real8_2d) :: sumwt
   real(r8), allocatable     ::  vectmp(:)  
   logical,  allocatable     ::  filter(:)
   integer :: xblk, yblk, xloc, yloc
   integer :: iblk, jblk, idata, ixseg, iyseg
   integer :: rmesg(3), smesg(3), isrc
   real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)
   integer :: xdsp, ydsp, xcnt, ycnt
   character(len=256) :: fileblock

   if(p_is_master)then 
      runoff_2d(:,:) = spval
   endif 

   call FLUSH_acc_cama_fluxes ()

   call accumulate_cama_fluxes 

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
   if (p_is_worker) then
      if (numpatch > 0) then
         allocate (filter (numpatch))
         allocate (vectmp (numpatch))
      end if
   end if
   if (p_is_io) then
      call allocate_block_data (ghist, sumwt)
   end if
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   ! ------------------------------------------------------------------------------------------
   ! Mapping the fluxes and state variables at patch [numpatch] to grid 
   ! ------------------------------------------------------------------------------------------
   if (p_is_worker) then
      if (numpatch > 0) then
         filter(:) = patchtype < 99
         vectmp (:) = 1.
      end if
   end if

   call mp2g_hist%map (vectmp, sumwt, spv = spval, msk = filter)
   
   call mp2g_hist%map (a_rnof_cama, f_rnof_cama, spv = spval, msk = filter)   

   if (p_is_io) then
      do yblk = 1, gblock%nyblk
         do xblk = 1, gblock%nxblk
            if (gblock%pio(xblk,yblk) == p_iam_glb) then
               do yloc = 1, ghist%ycnt(yblk) 
                  do xloc = 1, ghist%xcnt(xblk) 

                     if (sumwt%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                        IF (f_rnof_cama%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                           f_rnof_cama%blk(xblk,yblk)%val(xloc,yloc) &
                              = f_rnof_cama%blk(xblk,yblk)%val(xloc,yloc) &
                              / sumwt%blk(xblk,yblk)%val(xloc,yloc)
                        ENDIF
                     else
                        f_rnof_cama%blk(xblk,yblk)%val(xloc,yloc) = spval
                     end if

                  end do
               end do

            end if
         end do
      end do
   end if

#ifdef USEMPI
        CALL mpi_barrier (p_comm_glb, p_err)
#endif
        if (p_is_master) then
            do idata = 1, hist_block_info%ndatablk
                call mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, 10011, p_comm_glb, p_stat, p_err)
                isrc  = rmesg(1)
                ixseg = rmesg(2)
                iyseg = rmesg(3)

                xdsp = hist_block_info%xsegs(ixseg)%gdsp
                ydsp = hist_block_info%ysegs(iyseg)%gdsp
                xcnt = hist_block_info%xsegs(ixseg)%cnt
                ycnt = hist_block_info%ysegs(iyseg)%cnt

                allocate (rbuf(xcnt,ycnt))
                call mpi_recv (rbuf, xcnt * ycnt, MPI_DOUBLE, &
                        isrc, 10011, p_comm_glb, p_stat, p_err)
                runoff_2d (xdsp+1:xdsp+xcnt,ydsp+1:ydsp+ycnt) = rbuf
                deallocate (rbuf)
            end do

        elseif (p_is_io) then
            do iyseg = 1, hist_block_info%nyseg
                do ixseg = 1, hist_block_info%nxseg

                    iblk = hist_block_info%xsegs(ixseg)%blk
                    jblk = hist_block_info%ysegs(iyseg)%blk

                    if (gblock%pio(iblk,jblk) == p_iam_glb) then
                       xdsp = hist_block_info%xsegs(ixseg)%bdsp
                       ydsp = hist_block_info%ysegs(iyseg)%bdsp
                       xcnt = hist_block_info%xsegs(ixseg)%cnt
                       ycnt = hist_block_info%ysegs(iyseg)%cnt

                        allocate (sbuf (xcnt,ycnt))
                        sbuf = f_rnof_cama%blk(iblk,jblk)%val(xdsp+1:xdsp+xcnt,ydsp+1:ydsp+ycnt)

                        smesg = (/p_iam_glb, ixseg, iyseg/)
                        call mpi_send (smesg, 3, MPI_INTEGER, &
                        p_root, 10011, p_comm_glb, p_err) 
                        call mpi_send (sbuf, xcnt*ycnt, MPI_DOUBLE, &
                        p_root, 10011, p_comm_glb, p_err)

                        deallocate (sbuf)

                    end if
                end do
            end do
        end if


   if (allocated(filter)) deallocate(filter)
   if (allocated(vectmp)) deallocate(vectmp)

   END SUBROUTINE var_out

#endif

   ! -------
   subroutine flux_map_and_write_2d ( is_hist, &
         acc_vec, flux_xy, file_hist, varname, itime_in_file, sumwt, filter, longname, units)

      use precision
      use spmd_task
      use mod_namelist
      use mod_data_type
      use mod_mapping_pset2grid
      use mod_block
      use mod_grid
      use MOD_1D_Acc_Fluxes,  only: nac
      use GlobalVars, only: spval
      implicit none

      logical, intent(in) :: is_hist

      real(r8), intent(inout) :: acc_vec(:)
      type(block_data_real8_2d), intent(inout) :: flux_xy
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer, intent(in) :: itime_in_file
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units
 
      type(block_data_real8_2d), intent(in) :: sumwt
      logical, intent(in) :: filter(:)
      
      ! Local variables
      integer :: xblk, yblk, xloc, yloc
      integer :: compress

      if (.not. is_hist) return

      if (p_is_worker) then
         where (acc_vec /= spval)  acc_vec = acc_vec / nac
      end if
      
      call mp2g_hist%map (acc_vec, flux_xy, spv = spval, msk = filter)   

      if (p_is_io) then
         do yblk = 1, gblock%nyblk
            do xblk = 1, gblock%nxblk
               if (gblock%pio(xblk,yblk) == p_iam_glb) then

                  do yloc = 1, ghist%ycnt(yblk) 
                     do xloc = 1, ghist%xcnt(xblk) 

                        if (sumwt%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                           IF (flux_xy%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                              flux_xy%blk(xblk,yblk)%val(xloc,yloc) &
                                 = flux_xy%blk(xblk,yblk)%val(xloc,yloc) &
                                 / sumwt%blk(xblk,yblk)%val(xloc,yloc)
                           ENDIF
                        else
                           flux_xy%blk(xblk,yblk)%val(xloc,yloc) = spval
                        end if

                     end do
                  end do

               end if
            end do
         end do
      end if
      
      compress = DEF_HIST_COMPRESS_LEVEL 
      IF (present(longname) .and. present(units)) THEN
         call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy, &
            compress, longname, units) 
      ELSE
         call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy, compress) 
      ENDIF

   end subroutine flux_map_and_write_2d

   ! -------
   subroutine flux_map_and_write_3d ( is_hist, &
         acc_vec, flux_xy, file_hist, varname, dim1name, &
         itime_in_file, sumwt, filter, longname, units)

      use precision
      use spmd_task
      use mod_namelist
      use mod_data_type
      use mod_mapping_pset2grid
      use mod_block
      use mod_grid
      use MOD_1D_Acc_Fluxes,  only: nac
      use GlobalVars, only: spval
      implicit none

      logical, intent(in) :: is_hist

      real(r8), intent(inout) :: acc_vec(:,:)
      type(block_data_real8_3d), intent(inout) :: flux_xy
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: dim1name
      integer, intent(in) :: itime_in_file
      
      type(block_data_real8_2d), intent(in) :: sumwt
      logical, intent(in) :: filter(:)
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units

      ! Local variables
      integer :: xblk, yblk, xloc, yloc, i1
      integer :: compress

      if (.not. is_hist) return

      if (p_is_worker) then
         where (acc_vec /= spval)  acc_vec = acc_vec / nac
      end if
      
      call mp2g_hist%map (acc_vec, flux_xy, spv = spval, msk = filter)   

      if (p_is_io) then
         do yblk = 1, gblock%nyblk
            do xblk = 1, gblock%nxblk
               if (gblock%pio(xblk,yblk) == p_iam_glb) then

                  do yloc = 1, ghist%ycnt(yblk) 
                     do xloc = 1, ghist%xcnt(xblk) 

                        if (sumwt%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                           DO i1 = flux_xy%lb1, flux_xy%ub1
                              IF (flux_xy%blk(xblk,yblk)%val(i1,xloc,yloc) /= spval) THEN
                                 flux_xy%blk(xblk,yblk)%val(i1,xloc,yloc) &
                                    = flux_xy%blk(xblk,yblk)%val(i1,xloc,yloc) &
                                    / sumwt%blk(xblk,yblk)%val(xloc,yloc)
                              ENDIF
                           ENDDO 
                        else
                           flux_xy%blk(xblk,yblk)%val(:,xloc,yloc) = spval
                        end if

                     end do
                  end do

               end if
            end do
         end do
      end if
      
      compress = DEF_HIST_COMPRESS_LEVEL 
      IF (present(longname) .and. present(units)) THEN
         call hist_write_var_real8_3d (file_hist, varname, dim1name, ghist, &
            itime_in_file, flux_xy, compress, longname, units) 
      ELSE
         call hist_write_var_real8_3d (file_hist, varname, dim1name, ghist, &
            itime_in_file, flux_xy, compress) 
      ENDIF

   end subroutine flux_map_and_write_3d

   ! -------
   subroutine flux_map_and_write_4d ( is_hist, &
         acc_vec, flux_xy, file_hist, varname, dim1name, dim2name, itime_in_file, sumwt, filter, &
         longname, units)

      use precision
      use spmd_task
      use mod_namelist
      use mod_data_type
      use mod_mapping_pset2grid
      use mod_block
      use mod_grid
      use MOD_1D_Acc_Fluxes,  only: nac
      use GlobalVars, only: spval
      implicit none

      logical, intent(in) :: is_hist

      real(r8), intent(inout) :: acc_vec(:,:,:)
      type(block_data_real8_4d), intent(inout) :: flux_xy
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: dim1name, dim2name
      integer, intent(in) :: itime_in_file
      
      type(block_data_real8_2d), intent(in) :: sumwt
      logical, intent(in) :: filter(:)
      character (len=*), intent(in),optional :: longname
      character (len=*), intent(in),optional :: units

      ! Local variables
      integer :: xblk, yblk, xloc, yloc, i1, i2
      integer :: compress

      if (.not. is_hist) return

      if (p_is_worker) then
         where(acc_vec /= spval)  acc_vec = acc_vec / nac
      end if
      
      call mp2g_hist%map (acc_vec, flux_xy, spv = spval, msk = filter)   

      if (p_is_io) then
         do yblk = 1, gblock%nyblk
            do xblk = 1, gblock%nxblk
               if (gblock%pio(xblk,yblk) == p_iam_glb) then

                  do yloc = 1, ghist%ycnt(yblk) 
                     do xloc = 1, ghist%xcnt(xblk) 

                        if (sumwt%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                           DO i1 = flux_xy%lb1, flux_xy%ub1
                              DO i2 = flux_xy%lb2, flux_xy%ub2
                                 IF (flux_xy%blk(xblk,yblk)%val(i1,i2,xloc,yloc) /= spval) THEN
                                    flux_xy%blk(xblk,yblk)%val(i1,i2,xloc,yloc) &
                                       = flux_xy%blk(xblk,yblk)%val(i1,i2,xloc,yloc) &
                                       / sumwt%blk(xblk,yblk)%val(xloc,yloc)
                                 ENDIF
                              ENDDO
                           ENDDO
                        else
                           flux_xy%blk(xblk,yblk)%val(:,:,xloc,yloc) = spval
                        end if

                     end do
                  end do

               end if
            end do
         end do
      end if
      
      compress = DEF_HIST_COMPRESS_LEVEL 
      IF (present(longname) .and. present(units)) THEN
         call hist_write_var_real8_4d (file_hist, varname, dim1name, dim2name, &
            ghist, itime_in_file, flux_xy, compress, longname, units) 
      ELSE
         call hist_write_var_real8_4d (file_hist, varname, dim1name, dim2name, &
            ghist, itime_in_file, flux_xy, compress) 
      ENDIF

   end subroutine flux_map_and_write_4d

   ! -------
   subroutine flux_map_and_write_ln ( is_hist, &
         acc_vec, flux_xy, file_hist, varname, itime_in_file, sumwt, filter, longname, units)

      use precision
      use spmd_task
      use mod_namelist
      use mod_data_type
      use mod_mapping_pset2grid
      use mod_block
      use mod_grid
      use MOD_1D_Acc_Fluxes,  only: nac_ln
      use GlobalVars, only: spval
      implicit none

      logical, intent(in) :: is_hist

      real(r8), intent(inout) :: acc_vec(:)
      type(block_data_real8_2d), intent(inout) :: flux_xy
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer, intent(in) :: itime_in_file
      
      type(block_data_real8_2d), intent(in) :: sumwt
      logical,  intent(in) :: filter(:)
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units

      ! Local variables
      integer :: i, xblk, yblk, xloc, yloc
      integer :: compress

      if (.not. is_hist) return

      if (p_is_worker) then
         do i = lbound(acc_vec,1), ubound(acc_vec,1)
            if ((acc_vec(i) /= spval) .and. (nac_ln(i) > 0)) then
               acc_vec(i) = acc_vec(i) / nac_ln(i)
            end if
         end do
      end if
      
      call mp2g_hist%map (acc_vec, flux_xy, spv = spval, msk = filter)   

      if (p_is_io) then
         do yblk = 1, gblock%nyblk
            do xblk = 1, gblock%nxblk
               if (gblock%pio(xblk,yblk) == p_iam_glb) then

                  do yloc = 1, ghist%ycnt(yblk) 
                     do xloc = 1, ghist%xcnt(xblk) 

                        if ((sumwt%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) &
                           .and. (flux_xy%blk(xblk,yblk)%val(xloc,yloc) /= spval)) then
                              flux_xy%blk(xblk,yblk)%val(xloc,yloc) &
                                 = flux_xy%blk(xblk,yblk)%val(xloc,yloc) &
                                 / sumwt%blk(xblk,yblk)%val(xloc,yloc)
                        else
                           flux_xy%blk(xblk,yblk)%val(xloc,yloc) = spval
                        end if

                     end do
                  end do

               end if
            end do
         end do
      end if
      
      compress = DEF_HIST_COMPRESS_LEVEL 
      IF (present(longname) .and. present(units)) THEN
         call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy, &
            compress, longname, units) 
      ELSE
         call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy, &
            compress) 
      ENDIF

   end subroutine flux_map_and_write_ln

   subroutine set_hist_block_info (ghist, block_info)

      use mod_grid
      use mod_block
      USE mod_utils
      implicit none

      type(grid_type), intent(in) :: ghist
      TYPE(block_info_type), intent(out) :: block_info

      ! Local variables
      integer :: ilat_l, ilat_u, ilat, ilatloc, jblk, iyseg
      integer :: ilon_w, ilon_e, ilon, ilonloc, iblk, ixseg

      ilat_l = findloc(ghist%yblk /= 0, .true., dim=1)
      ilat_u = findloc(ghist%yblk /= 0, .true., dim=1, back=.true.)

      block_info%ginfo%nlat = ilat_u - ilat_l + 1
      allocate (block_info%ginfo%lat_s (block_info%ginfo%nlat))
      allocate (block_info%ginfo%lat_n (block_info%ginfo%nlat))
      allocate (block_info%ginfo%lat_c (block_info%ginfo%nlat))

      block_info%nyseg = 0
      jblk  = 0
      ilatloc = 0
      do ilat = ilat_l, ilat_u
         if (ghist%yblk(ilat) /= jblk) then
            block_info%nyseg = block_info%nyseg + 1
            jblk  = ghist%yblk(ilat)
         end if

         ilatloc = ilatloc + 1
         block_info%ginfo%lat_s(ilatloc) = ghist%lat_s(ilat)
         block_info%ginfo%lat_n(ilatloc) = ghist%lat_n(ilat)
         block_info%ginfo%lat_c(ilatloc) = (ghist%lat_s(ilat)+ghist%lat_n(ilat)) * 0.5
      end do

      allocate (block_info%ysegs (block_info%nyseg))

      iyseg = 0
      jblk  = 0
      do ilat = ilat_l, ilat_u
         if (ghist%yblk(ilat) /= jblk) then
            iyseg = iyseg + 1
            jblk  = ghist%yblk(ilat)
            block_info%ysegs(iyseg)%blk  = jblk
            block_info%ysegs(iyseg)%bdsp = ghist%yloc(ilat) - 1
            block_info%ysegs(iyseg)%gdsp = ilat - ilat_l
            block_info%ysegs(iyseg)%cnt  = 1
         else
            block_info%ysegs(iyseg)%cnt  = block_info%ysegs(iyseg)%cnt + 1
         end if
      end do

      if (all(ghist%xblk > 0)) then
         ilon_w = 1
         ilon_e = ghist%nlon
      else
         ilon_w = findloc(ghist%xblk /= 0, .true., dim=1)
         do while (.true.)
            ilon = ilon_w - 1
            if (ilon == 0) ilon = ghist%nlon

            if (ghist%xblk(ilon) /= 0) then
               ilon_w = ilon
            else
               exit
            end if
         end do

         ilon_e = ilon_w
         do while (.true.)
            ilon = mod(ilon_e,ghist%nlon) + 1

            if (ghist%xblk(ilon) /= 0) then
               ilon_e = ilon
            else
               exit
            end if
         end do
      end if

      block_info%ginfo%nlon = ilon_e - ilon_w + 1
      if (block_info%ginfo%nlon <= 0) THEN
         block_info%ginfo%nlon = block_info%ginfo%nlon + ghist%nlon
      ENDIF

      allocate (block_info%ginfo%lon_w (block_info%ginfo%nlon))
      allocate (block_info%ginfo%lon_e (block_info%ginfo%nlon))
      allocate (block_info%ginfo%lon_c (block_info%ginfo%nlon))

      block_info%nxseg = 0
      ilon = ilon_w - 1
      iblk = 0
      ilonloc = 0
      do while (.true.) 
         ilon = mod(ilon,ghist%nlon) + 1
         if (ghist%xblk(ilon) /= iblk) then
            block_info%nxseg = block_info%nxseg + 1
            iblk = ghist%xblk(ilon)
         end if

         ilonloc = ilonloc + 1
         block_info%ginfo%lon_w(ilonloc) = ghist%lon_w(ilon)
         block_info%ginfo%lon_e(ilonloc) = ghist%lon_e(ilon)

         block_info%ginfo%lon_c(ilonloc) = (ghist%lon_w(ilon) + ghist%lon_e(ilon)) * 0.5
         IF (ghist%lon_w(ilon) > ghist%lon_e(ilon)) THEN
            block_info%ginfo%lon_c(ilonloc) = block_info%ginfo%lon_c(ilonloc) + 180.0
            CALL normalize_longitude (block_info%ginfo%lon_c(ilonloc))
         ENDIF

         if (ilon == ilon_e) exit
      end do

      allocate (block_info%xsegs (block_info%nxseg))

      ixseg = 0
      iblk = 0
      ilon = ilon_w - 1
      ilonloc = 0
      do while (.true.) 
         ilon = mod(ilon,ghist%nlon) + 1
         ilonloc = ilonloc + 1
         if (ghist%xblk(ilon) /= iblk) then
            ixseg = ixseg + 1
            iblk = ghist%xblk(ilon)
            block_info%xsegs(ixseg)%blk  = iblk
            block_info%xsegs(ixseg)%bdsp = ghist%xloc(ilon) - 1
            block_info%xsegs(ixseg)%gdsp = ilonloc - 1
            block_info%xsegs(ixseg)%cnt = 1
         else
            block_info%xsegs(ixseg)%cnt = block_info%xsegs(ixseg)%cnt + 1
         end if

         if (ilon == ilon_e) exit
      end do

      block_info%ndatablk = 0

      do iyseg = 1, block_info%nyseg
         do ixseg = 1, block_info%nxseg
            iblk = block_info%xsegs(ixseg)%blk
            jblk = block_info%ysegs(iyseg)%blk
            if (gblock%pio(iblk,jblk) >= 0) then
               block_info%ndatablk = block_info%ndatablk + 1
            end if
         end do
      end do

   end subroutine set_hist_block_info



   !------------------------------
   subroutine hist_write_time ( &
         filename, dataname, grid, time, itime)

      use mod_namelist
      use mod_grid
      use mod_block
      use spmd_task
      use ncio_serial
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      type (grid_type),  intent(in) :: grid

      integer, intent(in)  :: time(3)
      integer, intent(out) :: itime

      ! Local variables
      character(len=256) :: fileblock
      integer :: iblk, jblk
      logical :: fexists

      if (trim(DEF_HIST_mode) == 'one') then
         if (p_is_master) then
            inquire (file=filename, exist=fexists)
            if (.not. fexists) then
               call ncio_create_file (trim(filename))
               CALL ncio_define_dimension(filename, 'time', 0)
               call ncio_define_dimension(filename, 'lat' , hist_block_info%ginfo%nlat)
               call ncio_define_dimension(filename, 'lon' , hist_block_info%ginfo%nlon)
#if(defined CaMa_Flood)
               call ncio_define_dimension(filename,'lat_cama', NY)
               call ncio_define_dimension(filename,'lon_cama', NX)
#endif
               call ncio_write_serial (filename, 'lat_s', hist_block_info%ginfo%lat_s, 'lat')
               call ncio_write_serial (filename, 'lat_n', hist_block_info%ginfo%lat_n, 'lat')
               call ncio_write_serial (filename, 'lon_w', hist_block_info%ginfo%lon_w, 'lon')
               call ncio_write_serial (filename, 'lon_e', hist_block_info%ginfo%lon_e, 'lon')
               call ncio_write_serial (filename, 'lat',   hist_block_info%ginfo%lat_c, 'lat')
               call ncio_write_serial (filename, 'lon',   hist_block_info%ginfo%lon_c, 'lon')
#if(defined CaMa_Flood)
               call ncio_write_serial (filename, 'lat_cama', D1LAT,'lat_cama')
               call ncio_write_serial (filename, 'lon_cama', D1LON,'lon_cama')
#endif
            endif
      
            call ncio_write_time (filename, dataname, time, itime)

         ENDIF

      elseif (trim(DEF_HIST_mode) == 'block') then

         if (p_is_io) then

            do jblk = 1, gblock%nyblk
               IF (ghist%ycnt(jblk) <= 0) cycle
               do iblk = 1, gblock%nxblk
                  IF (ghist%xcnt(iblk) <= 0) cycle
                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     call get_filename_block (filename, iblk, jblk, fileblock)

                     inquire (file=fileblock, exist=fexists)
                     if (.not. fexists) then
                        call ncio_create_file (trim(fileblock))
                        CALL ncio_define_dimension (fileblock, 'time', 0)
                        call hist_write_grid_info  (fileblock, grid, iblk, jblk)
                     end if

                     call ncio_write_time (fileblock, dataname, time, itime)

                  end if
               end do
            end do

         end if
      endif

   end subroutine hist_write_time

   !----------------------------------------------------------------------------
   subroutine hist_write_var_real8_2d ( &
         filename, dataname, grid, itime, wdata, compress, longname, units)

      use mod_namelist
      use mod_block
      use mod_grid
      use mod_data_type
      use spmd_task
      use ncio_serial
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      type (grid_type),  intent(in) :: grid
      integer, intent(in) :: itime

      type (block_data_real8_2d), intent(in) :: wdata

      integer, intent(in) :: compress
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units

      ! Local variables
      integer :: iblk, jblk, idata, ixseg, iyseg
      integer :: xcnt, ycnt, xbdsp, ybdsp, xgdsp, ygdsp
      integer :: rmesg(3), smesg(3), isrc
      character(len=256) :: fileblock
      real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)

      if (trim(DEF_HIST_mode) == 'one') then

         if (p_is_master) then

            allocate (vdata (hist_block_info%ginfo%nlon, hist_block_info%ginfo%nlat))
            vdata(:,:) = spval

#ifdef USEMPI
            do idata = 1, hist_block_info%ndatablk
               call mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
                  hist_data_id, p_comm_glb, p_stat, p_err)

               isrc  = rmesg(1)
               ixseg = rmesg(2)
               iyseg = rmesg(3)
                     
               xgdsp = hist_block_info%xsegs(ixseg)%gdsp
               ygdsp = hist_block_info%ysegs(iyseg)%gdsp
               xcnt = hist_block_info%xsegs(ixseg)%cnt
               ycnt = hist_block_info%ysegs(iyseg)%cnt

               allocate (rbuf(xcnt,ycnt))

               call mpi_recv (rbuf, xcnt*ycnt, MPI_DOUBLE, &
                  isrc, hist_data_id, p_comm_glb, p_stat, p_err)

               vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = rbuf
               deallocate (rbuf)

            end do
#else
            do iyseg = 1, hist_block_info%nyseg
               do ixseg = 1, hist_block_info%nxseg
                  iblk = hist_block_info%xsegs(ixseg)%blk
                  jblk = hist_block_info%ysegs(iyseg)%blk
                  if (gblock%pio(iblk,jblk) == p_iam_glb) then
                     xbdsp = hist_block_info%xsegs(ixseg)%bdsp
                     ybdsp = hist_block_info%ysegs(iyseg)%bdsp
                     xgdsp = hist_block_info%xsegs(ixseg)%gdsp
                     ygdsp = hist_block_info%ysegs(iyseg)%gdsp
                     xcnt = hist_block_info%xsegs(ixseg)%cnt
                     ycnt = hist_block_info%ysegs(iyseg)%cnt

                     vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                        wdata%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                  ENDIF
               ENDDO
            end do
#endif

            IF (present(longname) .and. present(units)) THEN
               call ncio_write_serial_time (filename, dataname, itime, vdata, &
                  'lon', 'lat', 'time',compress, longname, units)
            ELSE
               call ncio_write_serial_time (filename, dataname, itime, vdata, &
                  'lon', 'lat', 'time',compress)
            ENDIF

            deallocate (vdata)
         ENDIF

#ifdef USEMPI
         if (p_is_io) then
            do iyseg = 1, hist_block_info%nyseg
               do ixseg = 1, hist_block_info%nxseg

                  iblk = hist_block_info%xsegs(ixseg)%blk
                  jblk = hist_block_info%ysegs(iyseg)%blk

                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     xbdsp = hist_block_info%xsegs(ixseg)%bdsp
                     ybdsp = hist_block_info%ysegs(iyseg)%bdsp
                     xcnt = hist_block_info%xsegs(ixseg)%cnt
                     ycnt = hist_block_info%ysegs(iyseg)%cnt

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

            do jblk = 1, gblock%nyblk
               do iblk = 1, gblock%nxblk
                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     if ((grid%xcnt(iblk) == 0) .or. (grid%ycnt(jblk) == 0)) cycle

                     call get_filename_block (filename, iblk, jblk, fileblock)

                     call ncio_write_serial_time (fileblock, dataname, itime, &
                        wdata%blk(iblk,jblk)%val, 'lon', 'lat', 'time', compress)

                  end if
               end do
            end do

         end if
      end if

   end subroutine hist_write_var_real8_2d

   !----------------------------------------------------------------------------
   subroutine hist_write_var_real8_3d ( &
         filename, dataname, dim1name, grid, itime, wdata, compress, longname, units)

      use mod_namelist
      use mod_block
      use mod_grid
      use mod_data_type
      use spmd_task
      use ncio_serial
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      character (len=*), intent(in) :: dim1name
      type (grid_type),  intent(in) :: grid
      integer, intent(in) :: itime

      type (block_data_real8_3d), intent(in) :: wdata

      integer, intent(in) :: compress
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units

      ! Local variables
      integer :: iblk, jblk, idata, ixseg, iyseg
      integer :: xcnt, ycnt, ndim1, xbdsp, ybdsp, xgdsp, ygdsp 
      integer :: rmesg(4), smesg(4), isrc
      character(len=256) :: fileblock
      real(r8), allocatable :: rbuf(:,:,:), sbuf(:,:,:), vdata(:,:,:)

      if (trim(DEF_HIST_mode) == 'one') then

         if (p_is_master) then
            
#ifdef USEMPI
            do idata = 1, hist_block_info%ndatablk
            
               call mpi_recv (rmesg, 4, MPI_INTEGER, MPI_ANY_SOURCE, &
                  hist_data_id, p_comm_glb, p_stat, p_err)

               isrc  = rmesg(1)
               ixseg = rmesg(2)
               iyseg = rmesg(3)
               ndim1 = rmesg(4)
               
               xgdsp = hist_block_info%xsegs(ixseg)%gdsp
               ygdsp = hist_block_info%ysegs(iyseg)%gdsp
               xcnt = hist_block_info%xsegs(ixseg)%cnt
               ycnt = hist_block_info%ysegs(iyseg)%cnt

               allocate (rbuf (ndim1,xcnt,ycnt))

               call mpi_recv (rbuf, ndim1 * xcnt * ycnt, MPI_DOUBLE, &
                  isrc, hist_data_id, p_comm_glb, p_stat, p_err)

               IF (idata == 1) THEN
                  allocate (vdata (ndim1, hist_block_info%ginfo%nlon, hist_block_info%ginfo%nlat))
                  vdata(:,:,:) = spval
               ENDIF

               vdata (:,xgdsp+1:xgdsp+xcnt,ygdsp+1:ygdsp+ycnt) = rbuf

               deallocate (rbuf)
            end do
#else
            ndim1 = wdata%ub1 - wdata%lb1 + 1
            allocate (vdata (ndim1, hist_block_info%ginfo%nlon, hist_block_info%ginfo%nlat))
            vdata(:,:,:) = spval
                  
            do iyseg = 1, hist_block_info%nyseg
               do ixseg = 1, hist_block_info%nxseg
                  iblk = hist_block_info%xsegs(ixseg)%blk
                  jblk = hist_block_info%ysegs(iyseg)%blk
                  if (gblock%pio(iblk,jblk) == p_iam_glb) then
                     xbdsp = hist_block_info%xsegs(ixseg)%bdsp
                     ybdsp = hist_block_info%ysegs(iyseg)%bdsp
                     xgdsp = hist_block_info%xsegs(ixseg)%gdsp
                     ygdsp = hist_block_info%ysegs(iyseg)%gdsp
                     xcnt = hist_block_info%xsegs(ixseg)%cnt
                     ycnt = hist_block_info%ysegs(iyseg)%cnt

                     vdata (:,xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                        wdata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                  ENDIF
               end do
            ENDDO
#endif

            call ncio_define_dimension (filename, dim1name, ndim1) 

            IF (present(longname) .and. present(units)) THEN
               call ncio_write_serial_time (filename, dataname, itime, &
                  vdata, dim1name, 'lon', 'lat', 'time', compress, longname, units)
            ELSE
               call ncio_write_serial_time (filename, dataname, itime, &
                  vdata, dim1name, 'lon', 'lat', 'time', compress)
            ENDIF

            deallocate (vdata)
         ENDIF

#ifdef USEMPI
         if (p_is_io) then

            do iyseg = 1, hist_block_info%nyseg
               do ixseg = 1, hist_block_info%nxseg

                  iblk = hist_block_info%xsegs(ixseg)%blk
                  jblk = hist_block_info%ysegs(iyseg)%blk

                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     xbdsp = hist_block_info%xsegs(ixseg)%bdsp
                     ybdsp = hist_block_info%ysegs(iyseg)%bdsp
                     xcnt = hist_block_info%xsegs(ixseg)%cnt
                     ycnt = hist_block_info%ysegs(iyseg)%cnt
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

            do jblk = 1, gblock%nyblk
               do iblk = 1, gblock%nxblk
                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     if ((grid%xcnt(iblk) == 0) .or. (grid%ycnt(jblk) == 0)) cycle

                     call get_filename_block (filename, iblk, jblk, fileblock)

                     call ncio_define_dimension (fileblock, dim1name, wdata%ub1-wdata%lb1+1)

                     call ncio_write_serial_time (fileblock, dataname, itime, &
                        wdata%blk(iblk,jblk)%val, dim1name, 'lon', 'lat', 'time', compress)

                  end if
               end do
            end do

         end if
      end if

   end subroutine hist_write_var_real8_3d

   !----------------------------------------------------------------------------
   subroutine hist_write_var_real8_4d ( &
         filename, dataname, dim1name, dim2name, grid, itime, wdata, compress, longname, units)

      use mod_namelist
      use mod_block
      use mod_grid
      use mod_data_type
      use spmd_task
      use ncio_serial
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      character (len=*), intent(in) :: dim1name, dim2name
      type (grid_type),  intent(in) :: grid
      integer, intent(in) :: itime

      type (block_data_real8_4d), intent(in) :: wdata

      integer, intent(in) :: compress
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units
      ! Local variables
      integer :: iblk, jblk, idata, ixseg, iyseg
      integer :: xcnt, ycnt, ndim1, ndim2, xbdsp, ybdsp, xgdsp, ygdsp
      integer :: rmesg(5), smesg(5), isrc
      character(len=256) :: fileblock
      real(r8), allocatable :: rbuf(:,:,:,:), sbuf(:,:,:,:), vdata(:,:,:,:)

      if (trim(DEF_HIST_mode) == 'one') then

         if (p_is_master) then
               
#ifdef USEMPI
            do idata = 1, hist_block_info%ndatablk

               call mpi_recv (rmesg, 5, MPI_INTEGER, MPI_ANY_SOURCE, &
                  hist_data_id, p_comm_glb, p_stat, p_err)

               isrc  = rmesg(1)
               ixseg = rmesg(2)
               iyseg = rmesg(3)
               ndim1 = rmesg(4)
               ndim2 = rmesg(4)

               xgdsp = hist_block_info%xsegs(ixseg)%gdsp
               ygdsp = hist_block_info%ysegs(iyseg)%gdsp
               xcnt = hist_block_info%xsegs(ixseg)%cnt
               ycnt = hist_block_info%ysegs(iyseg)%cnt

               allocate (rbuf (ndim1,ndim2,xcnt,ycnt))

               call mpi_recv (rbuf, ndim1*ndim2*xcnt*ycnt, MPI_DOUBLE, &
                  isrc, hist_data_id, p_comm_glb, p_stat, p_err)
            
               IF (idata == 1) THEN
                  allocate (vdata (ndim1,ndim2,hist_block_info%ginfo%nlon,hist_block_info%ginfo%nlat))
                  vdata(:,:,:,:) = spval
               ENDIF 

               vdata (:,:,xgdsp+1:xgdsp+xcnt,ygdsp+1:ygdsp+ycnt) = rbuf

               deallocate (rbuf)
            end do
#else
            ndim1 = wdata%ub1 - wdata%lb1 + 1
            ndim2 = wdata%ub2 - wdata%lb2 + 1
            allocate (vdata (ndim1,ndim2,hist_block_info%ginfo%nlon,hist_block_info%ginfo%nlat))
            vdata(:,:,:,:) = spval

            do iyseg = 1, hist_block_info%nyseg
               do ixseg = 1, hist_block_info%nxseg
                  iblk = hist_block_info%xsegs(ixseg)%blk
                  jblk = hist_block_info%ysegs(iyseg)%blk
                  if (gblock%pio(iblk,jblk) == p_iam_glb) then
                     xbdsp = hist_block_info%xsegs(ixseg)%bdsp
                     ybdsp = hist_block_info%ysegs(iyseg)%bdsp
                     xgdsp = hist_block_info%xsegs(ixseg)%gdsp
                     ygdsp = hist_block_info%ysegs(iyseg)%gdsp
                     xcnt = hist_block_info%xsegs(ixseg)%cnt
                     ycnt = hist_block_info%ysegs(iyseg)%cnt

                     vdata (:,:,xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                        wdata%blk(iblk,jblk)%val(:,:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
                  ENDIF
               ENDDO
            ENDDO

#endif
            
            call ncio_define_dimension (filename, dim1name, ndim1) 
            call ncio_define_dimension (filename, dim2name, ndim2) 

            IF (present(longname) .and. present(units)) THEN
               call ncio_write_serial_time (filename, dataname, itime, vdata, dim1name, dim2name, &
                  'lon', 'lat', 'time', compress, longname, units)
            ELSE
               call ncio_write_serial_time (filename, dataname, itime, vdata, dim1name, dim2name, &
                  'lon', 'lat', 'time', compress)
            ENDIF

            deallocate (vdata)
         ENDIF

#ifdef USEMPI
         if (p_is_io) then

            do iyseg = 1, hist_block_info%nyseg
               do ixseg = 1, hist_block_info%nxseg

                  iblk = hist_block_info%xsegs(ixseg)%blk
                  jblk = hist_block_info%ysegs(iyseg)%blk

                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     xbdsp = hist_block_info%xsegs(ixseg)%bdsp
                     ybdsp = hist_block_info%ysegs(iyseg)%bdsp
                     xcnt = hist_block_info%xsegs(ixseg)%cnt
                     ycnt = hist_block_info%ysegs(iyseg)%cnt

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

            do jblk = 1, gblock%nyblk
               do iblk = 1, gblock%nxblk
                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     if ((grid%xcnt(iblk) == 0) .or. (grid%ycnt(jblk) == 0)) cycle

                     call get_filename_block (filename, iblk, jblk, fileblock)

                     call ncio_define_dimension (fileblock, dim1name, wdata%ub1-wdata%lb1+1)
                     call ncio_define_dimension (fileblock, dim2name, wdata%ub2-wdata%lb2+1)

                     IF (present(longname) .and. present(units)) THEN
                        call ncio_write_serial_time (fileblock, dataname, itime, &
                           wdata%blk(iblk,jblk)%val, dim1name, dim2name, 'lon', 'lat', 'time', compress, &
                           longname, units)
                     else
                        call ncio_write_serial_time (fileblock, dataname, itime, &
                           wdata%blk(iblk,jblk)%val, dim1name, dim2name, 'lon', 'lat', 'time', compress)
                     end if

                  end if
               end do
            end do

         end if
      end if

   end subroutine hist_write_var_real8_4d

   !------------------
   subroutine hist_write_grid_info (fileblock, grid, iblk, jblk)

      use mod_block
      use mod_grid
      use ncio_serial
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

   SUBROUTINE block_info_free_mem (this)

      IMPLICIT NONE

      TYPE(block_info_type) :: this

      IF (allocated(this%xsegs)) deallocate(this%xsegs)
      IF (allocated(this%ysegs)) deallocate(this%ysegs)

      IF (allocated(this%ginfo%lat_s)) deallocate(this%ginfo%lat_s)
      IF (allocated(this%ginfo%lat_n)) deallocate(this%ginfo%lat_n)
      IF (allocated(this%ginfo%lat_c)) deallocate(this%ginfo%lat_c)
      IF (allocated(this%ginfo%lon_w)) deallocate(this%ginfo%lon_w)
      IF (allocated(this%ginfo%lon_e)) deallocate(this%ginfo%lon_e)
      IF (allocated(this%ginfo%lon_c)) deallocate(this%ginfo%lon_c)

   END SUBROUTINE block_info_free_mem

end module mod_hist
