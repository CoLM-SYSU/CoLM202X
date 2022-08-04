#include <define.h>

module mod_hist

   use precision
   use mod_grid
   use mod_mapping_pset2grid
   USE mod_namelist
   USE GlobalVars, only : spval
   USE ncio_serial

   type(grid_type), target :: ghist
   type(mapping_pset2grid_type) :: mp2g_hist

   public :: hist_init
   public :: hist_out
   public :: hist_final

   TYPE(grid_concat_type) :: hist_concat

   integer :: hist_data_id

!--------------------------------------------------------------------------
contains

   !---------------------------------------
   subroutine hist_init (dir_hist, lon_res, lat_res)

      USE GlobalVars
      use spmd_task
      use mod_grid
      USE mod_landpatch
      use mod_mapping_pset2grid
      use MOD_1D_Acc_Fluxes
      implicit none

      character(len=*), intent(in) :: dir_hist
      REAL(r8), intent(in) :: lon_res
      REAL(r8), intent(in) :: lat_res

      ! Local Variables
      INTEGER :: lon_points, lat_points

      lon_points = nint(360.0/lon_res)
      lat_points = nint(180.0/lat_res)

      call ghist%define_by_ndims (lon_points, lat_points)
#ifndef CROP
      call mp2g_hist%build (landpatch, ghist)
#else
      call mp2g_hist%build (landpatch, ghist, pctcrop)
#endif

      call allocate_acc_fluxes ()
      call FLUSH_acc_fluxes ()

      !>>>>>add by zhongwang wei 
      call hist_concat%set (ghist)
#ifdef SinglePoint
      hist_concat%ginfo%lat_c(:) = SITE_lat_location
      hist_concat%ginfo%lon_c(:) = SITE_lon_location
#endif

      if (trim(DEF_HIST_mode) == 'one') then
         hist_data_id = 1000
      end if
      !<<<<<add by zhongwang wei


   end subroutine hist_init 

   !--------------------------------------
   subroutine hist_final ()

      use MOD_1D_Acc_Fluxes
      implicit none
      
      call deallocate_acc_fluxes ()

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
#ifdef USE_DEPTH_TO_BEDROCK
      USE MOD_TimeInvariants, only : ibedrock
#endif
#if(defined CaMa_Flood)
      use MOD_CaMa_Variables
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
#if(defined CaMa_Flood)
      character(LEN=256) :: file_hist_cama
      integer :: itime_in_file_cama
#endif
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

#if(defined CaMa_Flood)
         file_hist_cama = trim(dir_hist) // '/' // trim(site) //'_hist_cama_'//trim(cdate)//'.nc'
         call hist_write_cama_time (file_hist_cama, 'time', idate, itime_in_file_cama)
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
!#ifdef CLM5_INTERCEPTION
!         call flux_map_and_write_2d ( DEF_hist_vars%ldew_rain, &
!         a_ldew, f_ldew_rain, file_hist, 'f_ldew_rain', itime_in_file, sumwt, filter, &
!         'depth of rain on foliage','mm')
!         call flux_map_and_write_2d ( DEF_hist_vars%ldew_snow, &
!         a_ldew, f_ldew_snow, file_hist, 'f_ldew_snow', itime_in_file, sumwt, filter, &
!         'depth of snow on foliage','mm')
!#endif
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
            a_alb, f_alb, file_hist, 'f_alb', 'band', 'rtyp', itime_in_file, sumwt, filter, &
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

#ifdef BGC
         ! leaf carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafc, &
             a_leafc, f_leafc, file_hist, 'f_leafc', itime_in_file, sumwt, filter, &
             'leaf carbon display pool','gC/m2')

         ! leaf carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_storage, &
             a_leafc_storage, f_leafc_storage, file_hist, 'f_leafc_storage', itime_in_file, sumwt, filter, &
             'leaf carbon storage pool','gC/m2')

         ! leaf carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_xfer, &
             a_leafc_xfer, f_leafc_xfer, file_hist, 'f_leafc_xfer', itime_in_file, sumwt, filter, &
             'leaf carbon transfer pool','gC/m2')

         ! fine root carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootc, &
             a_frootc, f_frootc, file_hist, 'f_frootc', itime_in_file, sumwt, filter, &
             'fine root carbon display pool','gC/m2')

         ! fine root carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootc_storage, &
             a_frootc_storage, f_frootc_storage, file_hist, 'f_frootc_storage', itime_in_file, sumwt, filter, &
             'fine root carbon storage pool','gC/m2')

         ! fine root carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootc_xfer, &
             a_frootc_xfer, f_frootc_xfer, file_hist, 'f_frootc_xfer', itime_in_file, sumwt, filter, &
             'fine root carbon transfer pool','gC/m2')

         ! live stem carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemc, &
             a_livestemc, f_livestemc, file_hist, 'f_livestemc', itime_in_file, sumwt, filter, &
             'live stem carbon display pool','gC/m2')

         ! live stem carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemc_storage, &
             a_livestemc_storage, f_livestemc_storage, file_hist, 'f_livestemc_storage', itime_in_file, sumwt, filter, &
             'live stem carbon storage pool','gC/m2')

         ! live stem carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemc_xfer, &
             a_livestemc_xfer, f_livestemc_xfer, file_hist, 'f_livestemc_xfer', itime_in_file, sumwt, filter, &
             'live stem carbon transfer pool','gC/m2')

         ! dead stem carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemc, &
             a_deadstemc, f_deadstemc, file_hist, 'f_deadstemc', itime_in_file, sumwt, filter, &
             'dead stem carbon display pool','gC/m2')

         ! dead stem carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemc_storage, &
             a_deadstemc_storage, f_deadstemc_storage, file_hist, 'f_deadstemc_storage', itime_in_file, sumwt, filter, &
             'dead stem carbon storage pool','gC/m2')

         ! dead stem carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemc_xfer, &
             a_deadstemc_xfer, f_deadstemc_xfer, file_hist, 'f_deadstemc_xfer', itime_in_file, sumwt, filter, &
             'dead stem carbon transfer pool','gC/m2')

         ! live coarse root carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootc, &
             a_livecrootc, f_livecrootc, file_hist, 'f_livecrootc', itime_in_file, sumwt, filter, &
             'live coarse root carbon display pool','gC/m2')

         ! live coarse root carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootc_storage, &
             a_livecrootc_storage, f_livecrootc_storage, file_hist, 'f_livecrootc_storage', itime_in_file, sumwt, filter, &
             'live coarse root carbon storage pool','gC/m2')

         ! live coarse root carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootc_xfer, &
             a_livecrootc_xfer, f_livecrootc_xfer, file_hist, 'f_livecrootc_xfer', itime_in_file, sumwt, filter, &
             'live coarse root carbon transfer pool','gC/m2')

         ! dead coarse root carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootc, &
             a_deadcrootc, f_deadcrootc, file_hist, 'f_deadcrootc', itime_in_file, sumwt, filter, &
             'dead coarse root carbon display pool','gC/m2')

         ! dead coarse root carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootc_storage, &
             a_deadcrootc_storage, f_deadcrootc_storage, file_hist, 'f_deadcrootc_storage', itime_in_file, sumwt, filter, &
             'dead coarse root carbon storage pool','gC/m2')

         ! dead coarse root carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootc_xfer, &
             a_deadcrootc_xfer, f_deadcrootc_xfer, file_hist, 'f_deadcrootc_xfer', itime_in_file, sumwt, filter, &
             'dead coarse root carbon transfer pool','gC/m2')

         ! grain carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainc, &
             a_grainc, f_grainc, file_hist, 'f_grainc', itime_in_file, sumwt, filter, &
             'grain carbon display pool','gC/m2')

         ! grain carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainc_storage, &
             a_grainc_storage, f_grainc_storage, file_hist, 'f_grainc_storage', itime_in_file, sumwt, filter, &
             'grain carbon storage pool','gC/m2')

         ! grain carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainc_xfer, &
             a_grainc_xfer, f_grainc_xfer, file_hist, 'f_grainc_xfer', itime_in_file, sumwt, filter, &
             'grain carbon transfer pool','gC/m2')

         ! leaf nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafn, &
             a_leafn, f_leafn, file_hist, 'f_leafn', itime_in_file, sumwt, filter, &
             'leaf nitrogen display pool','gN/m2')

         ! leaf nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafn_storage, &
             a_leafn_storage, f_leafn_storage, file_hist, 'f_leafn_storage', itime_in_file, sumwt, filter, &
             'leaf nitrogen storage pool','gN/m2')

         ! leaf nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafn_xfer, &
             a_leafn_xfer, f_leafn_xfer, file_hist, 'f_leafn_xfer', itime_in_file, sumwt, filter, &
             'leaf nitrogen transfer pool','gN/m2')

         ! fine root nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootn, &
             a_frootn, f_frootn, file_hist, 'f_frootn', itime_in_file, sumwt, filter, &
             'fine root nitrogen display pool','gN/m2')

         ! fine root nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootn_storage, &
             a_frootn_storage, f_frootn_storage, file_hist, 'f_frootn_storage', itime_in_file, sumwt, filter, &
             'fine root nitrogen storage pool','gN/m2')

         ! fine root nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootn_xfer, &
             a_frootn_xfer, f_frootn_xfer, file_hist, 'f_frootn_xfer', itime_in_file, sumwt, filter, &
             'fine root nitrogen transfer pool','gN/m2')

         ! live stem nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemn, &
             a_livestemn, f_livestemn, file_hist, 'f_livestemn', itime_in_file, sumwt, filter, &
             'live stem nitrogen display pool','gN/m2')

         ! live stem nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemn_storage, &
             a_livestemn_storage, f_livestemn_storage, file_hist, 'f_livestemn_storage', itime_in_file, sumwt, filter, &
             'live stem nitrogen storage pool','gN/m2')

         ! live stem nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemn_xfer, &
             a_livestemn_xfer, f_livestemn_xfer, file_hist, 'f_livestemn_xfer', itime_in_file, sumwt, filter, &
             'live stem nitrogen transfer pool','gN/m2')

         ! dead stem nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemn, &
             a_deadstemn, f_deadstemn, file_hist, 'f_deadstemn', itime_in_file, sumwt, filter, &
             'dead stem nitrogen display pool','gN/m2')

         ! dead stem nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemn_storage, &
             a_deadstemn_storage, f_deadstemn_storage, file_hist, 'f_deadstemn_storage', itime_in_file, sumwt, filter, &
             'dead stem nitrogen storage pool','gN/m2')

         ! dead stem nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemn_xfer, &
             a_deadstemn_xfer, f_deadstemn_xfer, file_hist, 'f_deadstemn_xfer', itime_in_file, sumwt, filter, &
             'dead stem nitrogen transfer pool','gN/m2')

         ! live coarse root nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootn, &
             a_livecrootn, f_livecrootn, file_hist, 'f_livecrootn', itime_in_file, sumwt, filter, &
             'live coarse root nitrogen display pool','gN/m2')

         ! live coarse root nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootn_storage, &
             a_livecrootn_storage, f_livecrootn_storage, file_hist, 'f_livecrootn_storage', itime_in_file, sumwt, filter, &
             'live coarse root nitrogen storage pool','gN/m2')

         ! live coarse root nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootn_xfer, &
             a_livecrootn_xfer, f_livecrootn_xfer, file_hist, 'f_livecrootn_xfer', itime_in_file, sumwt, filter, &
             'live coarse root nitrogen transfer pool','gN/m2')

         ! dead coarse root nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootn, &
             a_deadcrootn, f_deadcrootn, file_hist, 'f_deadcrootn', itime_in_file, sumwt, filter, &
             'dead coarse root nitrogen display pool','gN/m2')

         ! dead coarse root nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootn_storage, &
             a_deadcrootn_storage, f_deadcrootn_storage, file_hist, 'f_deadcrootn_storage', itime_in_file, sumwt, filter, &
             'dead coarse root nitrogen storage pool','gN/m2')

         ! dead coarse root nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootn_xfer, &
             a_deadcrootn_xfer, f_deadcrootn_xfer, file_hist, 'f_deadcrootn_xfer', itime_in_file, sumwt, filter, &
             'dead coarse root nitrogen transfer pool','gN/m2')

         ! grain nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainn, &
             a_grainn, f_grainn, file_hist, 'f_grainn', itime_in_file, sumwt, filter, &
             'grain nitrogen display pool','gN/m2')

         ! grain nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainn_storage, &
             a_grainn_storage, f_grainn_storage, file_hist, 'f_grainn_storage', itime_in_file, sumwt, filter, &
             'grain nitrogen storage pool','gN/m2')

         ! grain nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainn_xfer, &
             a_grainn_xfer, f_grainn_xfer, file_hist, 'f_grainn_xfer', itime_in_file, sumwt, filter, &
             'grain nitrogen transfer pool','gN/m2')

         ! retranslocation nitrogen pool
         call flux_map_and_write_2d ( DEF_hist_vars%retrasn, &
             a_retransn, f_retransn, file_hist, 'f_retrasn', itime_in_file, sumwt, filter, &
             'retranslocation nitrogen pool','gN/m2')

         ! gross primary productivity
         call flux_map_and_write_2d ( DEF_hist_vars%gpp, &
             a_gpp, f_gpp, file_hist, 'f_gpp', itime_in_file, sumwt, filter, &
             'gross primary productivity','gC/m2/s')

         ! gross primary productivity
         call flux_map_and_write_2d ( DEF_hist_vars%downreg, &
             a_downreg, f_downreg, file_hist, 'f_downreg', itime_in_file, sumwt, filter, &
             'gpp downregulation due to N limitation','unitless')

         ! autotrophic respiration
         call flux_map_and_write_2d ( DEF_hist_vars%ar , &
             a_ar , f_ar , file_hist, 'f_ar', itime_in_file, sumwt, filter, &
             'autotrophic respiration','gC/m2/s')

#ifdef CROP
         ! crop phase
         call flux_map_and_write_2d ( DEF_hist_vars%cphase, &
             a_cphase, f_cphase, file_hist, 'f_cphase', itime_in_file, sumwt, filter, &
             'crop phase','unitless')

         ! 1-yr crop production carbon
         call flux_map_and_write_2d ( DEF_hist_vars%cropprod1c, &
             a_cropprod1c, f_cropprod1c, file_hist, 'f_cropprod1c', itime_in_file, sumwt, filter, &
             '1-yr crop production carbon','gC/m2')

         ! loss rate of 1-yr crop production carbon
         call flux_map_and_write_2d ( DEF_hist_vars%cropprod1c_loss, &
             a_cropprod1c_loss, f_cropprod1c_loss, file_hist, 'f_cropprod1c_loss', itime_in_file, sumwt, filter, &
             'loss rate of 1-yr crop production carbon','gC/m2/s')

         ! crop seed deficit
         call flux_map_and_write_2d ( DEF_hist_vars%cropseedc_deficit, &
             a_cropseedc_deficit, f_cropseedc_deficit, file_hist, 'f_cropseedc_deficit', itime_in_file, sumwt, filter, &
             'crop seed deficit','gC/m2/s')

         ! grain to crop production carbon
         call flux_map_and_write_2d ( DEF_hist_vars%grainc_to_cropprodc, &
             a_grainc_to_cropprodc, f_grainc_to_cropprodc, file_hist, 'f_grainc_to_cropprodc', itime_in_file, sumwt, filter, &
             'grain to crop production carbon','gC/m2/s')

         ! grain to crop seed carbon
         call flux_map_and_write_2d ( DEF_hist_vars%grainc_to_seed, &
             a_grainc_to_seed, f_grainc_to_seed, file_hist, 'f_grainc_to_seed', itime_in_file, sumwt, filter, &
             'grain to crop seed carbon','gC/m2/s')
#endif

         ! litter 1 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr1c_vr, &
            a_litr1c_vr, f_litr1c_vr, file_hist, 'f_litr1c_vr', 'soil', &
            itime_in_file, sumwt, filter,'litter 1 carbon density in soil layers','gC/m3')

         ! litter 2 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr2c_vr, &
            a_litr2c_vr, f_litr2c_vr, file_hist, 'f_litr2c_vr', 'soil', &
            itime_in_file, sumwt, filter,'litter 2 carbon density in soil layers','gC/m3')

         ! litter 3 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr3c_vr, &
            a_litr3c_vr, f_litr3c_vr, file_hist, 'f_litr3c_vr', 'soil', &
            itime_in_file, sumwt, filter,'litter 3 carbon density in soil layers','gC/m3')

         ! soil 1 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil1c_vr, &
            a_soil1c_vr, f_soil1c_vr, file_hist, 'f_soil1c_vr', 'soil', &
            itime_in_file, sumwt, filter,'soil 1 carbon density in soil layers','gC/m3')

         ! soil 2 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil2c_vr, &
            a_soil2c_vr, f_soil2c_vr, file_hist, 'f_soil2c_vr', 'soil', &
            itime_in_file, sumwt, filter,'soil 2 carbon density in soil layers','gC/m3')

         ! soil 3 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil3c_vr, &
            a_soil3c_vr, f_soil3c_vr, file_hist, 'f_soil3c_vr', 'soil', &
            itime_in_file, sumwt, filter,'soil 3 carbon density in soil layers','gC/m3')

         ! coarse woody debris carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%cwdc_vr, &
            a_cwdc_vr, f_cwdc_vr, file_hist, 'f_cwdc_vr', 'soil', &
            itime_in_file, sumwt, filter,'coarse woody debris carbon density in soil layers','gC/m3')

         ! litter 1 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr1n_vr, &
            a_litr1n_vr, f_litr1n_vr, file_hist, 'f_litr1n_vr', 'soil', &
            itime_in_file, sumwt, filter,'litter 1 nitrogen density in soil layers','gN/m3')

         ! litter 2 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr2n_vr, &
            a_litr2n_vr, f_litr2n_vr, file_hist, 'f_litr2n_vr', 'soil', &
            itime_in_file, sumwt, filter,'litter 2 nitrogen density in soil layers','gN/m3')

         ! litter 3 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr3n_vr, &
            a_litr3n_vr, f_litr3n_vr, file_hist, 'f_litr3n_vr', 'soil', &
            itime_in_file, sumwt, filter,'litter 3 nitrogen density in soil layers','gN/m3')

         ! soil 1 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil1n_vr, &
            a_soil1n_vr, f_soil1n_vr, file_hist, 'f_soil1n_vr', 'soil', &
            itime_in_file, sumwt, filter,'soil 1 nitrogen density in soil layers','gN/m3')

         ! soil 2 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil2n_vr, &
            a_soil2n_vr, f_soil2n_vr, file_hist, 'f_soil2n_vr', 'soil', &
            itime_in_file, sumwt, filter,'soil 2 nitrogen density in soil layers','gN/m3')

         ! soil 3 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil3n_vr, &
            a_soil3n_vr, f_soil3n_vr, file_hist, 'f_soil3n_vr', 'soil', &
            itime_in_file, sumwt, filter,'soil 3 nitrogen density in soil layers','gN/m3')

         ! coarse woody debris nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%cwdn_vr, &
            a_cwdn_vr, f_cwdn_vr, file_hist, 'f_cwdn_vr', 'soil', &
            itime_in_file, sumwt, filter,'coarse woody debris nitrogen density in soil layers','gN/m3')

         ! mineral nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%sminn_vr, &
            a_sminn_vr, f_sminn_vr, file_hist, 'f_sminn_vr', 'soil', &
            itime_in_file, sumwt, filter,'mineral nitrogen density in soil layers','gN/m3')

#endif
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
            a_rootr, f_rootr, file_hist, 'f_rootr', 'soil', itime_in_file, sumwt, filter, &
            'root water uptake', 'mm h2o/s')

#ifdef PLANT_HYDRAULIC_STRESS
         ! vegetation water potential [mm]
         call flux_map_and_write_3d ( DEF_hist_vars%vegwp, &
            a_vegwp, f_vegwp, file_hist, 'f_vegwp', 'vegnodes', itime_in_file, sumwt, filter, &
            'vegetation water potential', 'mm')
#endif

#ifdef VARIABLY_SATURATED_FLOW
         ! depth of ponding water [m]
         call flux_map_and_write_2d ( DEF_hist_vars%dpond, &
            a_dpond, f_dpond, file_hist, 'f_dpond', itime_in_file, sumwt, filter, &
            'depth of ponding water','mm')
#endif

#ifndef USE_DEPTH_TO_BEDROCK
         ! water table depth [m]
         call flux_map_and_write_2d ( DEF_hist_vars%zwt, &
            a_zwt, f_zwt, file_hist, 'f_zwt', itime_in_file, sumwt, filter, &
            'the depth to water table','m')

         ! water storage in aquifer [mm]
         call flux_map_and_write_2d ( DEF_hist_vars%wa, &
            a_wa, f_wa, file_hist, 'f_wa', itime_in_file, sumwt, filter, &
            'water storage in aquifer','mm')
#else
         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = (patchtype <= 2) .and. (ibedrock <= nl_soil)
               vectmp(:) = 1.
            end if
         end if

         call mp2g_hist%map (vectmp, sumwt, spv = spval, msk = filter)
         
         ! water table depth [m]
         call flux_map_and_write_2d ( DEF_hist_vars%dwatsub, &
            a_dwatsub, f_dwatsub, file_hist, 'f_dwatsub', itime_in_file, sumwt, filter, &
            'depth of saturated subsurface water above bedrock','m')

         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = (patchtype <= 2) .and. (ibedrock > nl_soil)
               vectmp(:) = 1.
            end if
         end if

         call mp2g_hist%map (vectmp, sumwt, spv = spval, msk = filter)
         
         ! water table depth [m]
         call flux_map_and_write_2d ( DEF_hist_vars%zwt, &
            a_zwt, f_zwt, file_hist, 'f_zwt', itime_in_file, sumwt, filter, &
            'the depth to water table','m')

         ! water storage in aquifer [mm]
         call flux_map_and_write_2d ( DEF_hist_vars%wa, &
            a_wa, f_wa, file_hist, 'f_wa', itime_in_file, sumwt, filter, &
            'water storage in aquifer','mm')
#endif

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
            'q* in similarity theory', 'kg/kg')

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
         CALL hist_out_cama (file_hist_cama, itime_in_file_cama)
#endif

         if (allocated(filter)) deallocate(filter)
         if (allocated(vectmp)) deallocate(vectmp)

         call FLUSH_acc_fluxes ()

      end if

   END SUBROUTINE hist_out
   
   ! -------
   subroutine flux_map_and_write_2d ( is_hist, &
         acc_vec, flux_xy, file_hist, varname, itime_in_file, sumwt, filter, &
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

      real(r8), intent(inout) :: acc_vec(:)
      type(block_data_real8_2d), intent(inout) :: flux_xy
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer,          intent(in) :: itime_in_file
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
 
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
      call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy, compress) 

      IF ((itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
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
      character (len=*), intent(in) :: longname
      character (len=*), intent(in) :: units

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
      call hist_write_var_real8_3d (file_hist, varname, dim1name, ghist, &
         itime_in_file, flux_xy, compress) 

      IF ((itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
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
      character (len=*), intent(in) :: longname
      character (len=*), intent(in) :: units

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
      call hist_write_var_real8_4d (file_hist, varname, dim1name, dim2name, &
         ghist, itime_in_file, flux_xy, compress) 

      IF ((itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_4d

   ! -------
   subroutine flux_map_and_write_ln ( is_hist, &
         acc_vec, flux_xy, file_hist, varname, itime_in_file, sumwt, filter, &
         longname, units)

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
      call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy, &
         compress) 

      IF ((itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_ln

   !------------------------------
   subroutine hist_write_time ( &
         filename, dataname, grid, time, itime)

      use mod_namelist
      use mod_grid
      use mod_block
      use spmd_task
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
               call ncio_define_dimension(filename, 'lat' , hist_concat%ginfo%nlat)
               call ncio_define_dimension(filename, 'lon' , hist_concat%ginfo%nlon)

               call ncio_write_serial (filename, 'lat',   hist_concat%ginfo%lat_c, 'lat')
               call ncio_write_serial (filename, 'lon',   hist_concat%ginfo%lon_c, 'lon')
#ifndef SinglePoint
               call ncio_write_serial (filename, 'lat_s', hist_concat%ginfo%lat_s, 'lat')
               call ncio_write_serial (filename, 'lat_n', hist_concat%ginfo%lat_n, 'lat')
               call ncio_write_serial (filename, 'lon_w', hist_concat%ginfo%lon_w, 'lon')
               call ncio_write_serial (filename, 'lon_e', hist_concat%ginfo%lon_e, 'lon')
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
         filename, dataname, grid, itime, wdata, compress)

      use mod_namelist
      use mod_block
      use mod_grid
      use mod_data_type
      use spmd_task
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      type (grid_type),  intent(in) :: grid
      integer, intent(in) :: itime

      type (block_data_real8_2d), intent(in) :: wdata

      integer, intent(in) :: compress

      ! Local variables
      integer :: iblk, jblk, idata, ixseg, iyseg
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
         filename, dataname, dim1name, grid, itime, wdata, compress)

      use mod_namelist
      use mod_block
      use mod_grid
      use mod_data_type
      use spmd_task
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      character (len=*), intent(in) :: dim1name
      type (grid_type),  intent(in) :: grid
      integer, intent(in) :: itime

      type (block_data_real8_3d), intent(in) :: wdata

      integer, intent(in) :: compress

      ! Local variables
      integer :: iblk, jblk, idata, ixseg, iyseg
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
         filename, dataname, dim1name, dim2name, grid, itime, wdata, compress)

      use mod_namelist
      use mod_block
      use mod_grid
      use mod_data_type
      use spmd_task
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      character (len=*), intent(in) :: dim1name, dim2name
      type (grid_type),  intent(in) :: grid
      integer, intent(in) :: itime

      type (block_data_real8_4d), intent(in) :: wdata

      integer, intent(in) :: compress

      ! Local variables
      integer :: iblk, jblk, idata, ixseg, iyseg
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

            do jblk = 1, gblock%nyblk
               do iblk = 1, gblock%nxblk
                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     if ((grid%xcnt(iblk) == 0) .or. (grid%ycnt(jblk) == 0)) cycle

                     call get_filename_block (filename, iblk, jblk, fileblock)

                     call ncio_define_dimension (fileblock, dim1name, wdata%ub1-wdata%lb1+1)
                     call ncio_define_dimension (fileblock, dim2name, wdata%ub2-wdata%lb2+1)

                     call ncio_write_serial_time (fileblock, dataname, itime, &
                        wdata%blk(iblk,jblk)%val, dim1name, dim2name, 'lon', 'lat', 'time', compress)

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

end module mod_hist
