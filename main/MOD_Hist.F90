#include <define.h>

module MOD_Hist

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
#ifdef URBAN_MODEL
   USE MOD_LandUrban
#endif
#ifdef LULC_IGBP_PFT
   USE MOD_Vars_PFTimeInvariants, only: pftclass
   USE MOD_LandPFT, only : patch_pft_s
#endif
   USE MOD_Vars_Global, only : spval
   USE MOD_NetCDFSerial
#if (defined UNSTRUCTURED || defined CATCHMENT)
   USE MOD_HistVector
#endif

   type(grid_type), target :: ghist
   TYPE(grid_type), target :: ghist_urb
   type(mapping_pset2grid_type) :: mp2g_hist
   type(mapping_pset2grid_type) :: mp2g_hist_urb

   public :: hist_init
   public :: hist_out
   public :: hist_final

   TYPE(grid_concat_type) :: hist_concat

   integer :: hist_data_id

!--------------------------------------------------------------------------
contains

   !---------------------------------------
   subroutine hist_init (dir_hist, lon_res, lat_res)

      USE MOD_Vars_Global
      use MOD_SPMD_Task
      use MOD_Grid
      USE MOD_LandPatch
      use MOD_Mapping_Pset2Grid
      use MOD_Vars_1DAccFluxes
#ifdef LATERAL_FLOW
      USE MOD_Hydro_Hist
#endif
      USE MOD_Forcing, only : gforc
      implicit none

      character(len=*), intent(in) :: dir_hist
      REAL(r8), intent(in) :: lon_res
      REAL(r8), intent(in) :: lat_res

      call allocate_acc_fluxes ()
      call FLUSH_acc_fluxes ()

#ifdef LATERAL_FLOW
      CALL hist_basin_init ()
#endif

#if (defined UNSTRUCTURED || defined CATCHMENT)
      IF (DEF_HISTORY_IN_VECTOR) THEN
         RETURN
      ENDIF
#endif

      IF (DEF_hist_grid_as_forcing) then
         CALL ghist%define_by_copy (gforc)
      ELSE
         call ghist%define_by_res (lon_res, lat_res)
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


   end subroutine hist_init

   !--------------------------------------
   subroutine hist_final ()

      use MOD_Vars_1DAccFluxes
#ifdef LATERAL_FLOW
      USE MOD_Hydro_Hist
#endif
      implicit none

      call deallocate_acc_fluxes ()

#ifdef LATERAL_FLOW
      CALL hist_basin_final ()
#endif

   end subroutine hist_final

   !---------------------------------------
   SUBROUTINE hist_out (idate, deltim, itstamp, ptstamp, &
         dir_hist, site)

      !=======================================================================
      ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
      !=======================================================================

      use MOD_Precision
      use MOD_Namelist
      use MOD_TimeManager
      use MOD_SPMD_Task
      use MOD_Vars_2DFluxes
      use MOD_Vars_1DAccFluxes
      use MOD_Block
      use MOD_DataType
      use MOD_LandPatch
      use MOD_Mapping_Pset2Grid
      use MOD_Vars_2DFluxes
      use MOD_CoLMDebug
      use MOD_Vars_Global, only : spval
      USE MOD_Vars_TimeInvariants, only : patchtype, patchclass
#ifdef URBAN_MODEL
      USE MOD_LandUrban
#endif
#if(defined CaMa_Flood)
      use MOD_CaMa_Vars !defination of CaMa variables
#endif
#ifdef LATERAL_FLOW
      USE MOD_Hydro_Hist
#endif
      USE MOD_Forcing, only : forcmask
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

      type(block_data_real8_2d) :: sumarea
      type(block_data_real8_2d) :: sumarea_urb
      real(r8), allocatable ::  vectmp(:)
      real(r8), allocatable ::  vecacc(:)
      logical,  allocatable ::  filter(:)

      integer i, u
#ifdef URBAN_MODEL
      real(r8), allocatable ::  vectmp_urb(:)
      logical,  allocatable ::  filter_urb(:)
#endif

      if (itstamp <= ptstamp) then
         call FLUSH_acc_fluxes ()
         return
      else

         call accumulate_fluxes

      end if

      select case (trim(adjustl(DEF_HIST_FREQ)))
      case ('TIMESTEP')
         lwrite = .true.
      case ('HOURLY')
         lwrite = isendofhour (idate, deltim)
      case ('DAILY')
         lwrite = isendofday(idate, deltim)
      case ('MONTHLY')
         lwrite = isendofmonth(idate, deltim)
      case ('YEARLY')
         lwrite = isendofyear(idate, deltim)
      case default
         write(*,*) 'Warning : Please use one of TIMESTEP/HOURLY/DAILY/MONTHLY/YEARLY for history frequency.'
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
         else
            write(*,*) 'Warning : Please use one of DAY/MONTH/YEAR for history group.'
         end if

#if(defined CaMa_Flood)
         !add variables to write cama-flood output.
         file_hist_cama = trim(dir_hist) // '/' // trim(site) //'_hist_cama_'//trim(cdate)//'.nc' !file name of cama-flood output
         call hist_write_cama_time (file_hist_cama, 'time', idate, itime_in_file_cama)         ! write CaMa-Flood output
#endif

         file_hist = trim(dir_hist) // '/' // trim(site) //'_hist_'//trim(cdate)//'.nc'

#ifdef LATERAL_FLOW
         CALL hist_basin_out (file_hist, idate)
#endif

#if (defined UNSTRUCTURED || defined CATCHMENT)
         IF (DEF_HISTORY_IN_VECTOR) THEN
            CALL hist_vector_out (file_hist, idate)
            RETURN
         ENDIF
#endif

         call hist_write_time (file_hist, 'time', ghist, idate, itime_in_file)

         if (p_is_worker) then
            if (numpatch > 0) then
               allocate (filter (numpatch))
               allocate (vectmp (numpatch))
               allocate (vecacc (numpatch))
            end if
#ifdef URBAN_MODEL
            IF (numurban > 0) THEN
               allocate (filter_urb (numurban))
               allocate (vectmp_urb (numurban))
            ENDIF
#endif
         end if

         if (p_is_io) then
            call allocate_block_data (ghist, sumarea)
#ifdef URBAN_MODEL
            call allocate_block_data (ghist, sumarea_urb)
#endif
         end if

         ! ---------------------------------------------------
         ! Meteorological forcing
         ! ---------------------------------------------------
         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype < 99
               vectmp(:) = 1.
               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask
               ENDIF
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! wind in eastward direction [m/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_us, &
            a_us, f_xy_us, file_hist, 'f_xy_us', itime_in_file, sumarea, filter, &
            'wind in eastward direction', 'm/s')

         ! wind in northward direction [m/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_vs, &
            a_vs, f_xy_vs, file_hist, 'f_xy_vs', itime_in_file, sumarea, filter, &
            'wind in northward direction','m/s')

         ! temperature at reference height [kelvin]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_t, &
            a_t, f_xy_t, file_hist, 'f_xy_t', itime_in_file, sumarea, filter, &
            'temperature at reference height','kelvin')

         ! specific humidity at reference height [kg/kg]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_q, &
            a_q, f_xy_q, file_hist, 'f_xy_q', itime_in_file, sumarea, filter, &
            'specific humidity at reference height','kg/kg')

         ! convective precipitation [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_prc, &
            a_prc, f_xy_prc, file_hist, 'f_xy_prc', itime_in_file, sumarea, filter, &
            'convective precipitation','mm/s')

         ! large scale precipitation [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_prl, &
            a_prl, f_xy_prl, file_hist, 'f_xy_prl', itime_in_file, sumarea, filter, &
            'large scale precipitation','mm/s')

         ! atmospheric pressure at the surface [pa]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_pbot, &
            a_pbot, f_xy_pbot, file_hist, 'f_xy_pbot', itime_in_file, sumarea, filter, &
            'atmospheric pressure at the surface','pa')

         ! atmospheric infrared (longwave) radiation [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_frl, &
            a_frl, f_xy_frl, file_hist, 'f_xy_frl', itime_in_file, sumarea, filter, &
            'atmospheric infrared (longwave) radiation','W/m2')

         ! downward solar radiation at surface [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_solarin, &
            a_solarin, f_xy_solarin, file_hist, 'f_xy_solarin', itime_in_file, sumarea, filter, &
            'downward solar radiation at surface','W/m2')

         ! rain [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_rain, &
            a_rain, f_xy_rain, file_hist, 'f_xy_rain', itime_in_file, sumarea, filter, &
            'rain','mm/s')

         ! snow [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xy_snow, &
            a_snow, f_xy_snow, file_hist, 'f_xy_snow', itime_in_file, sumarea, filter, &
            'snow','mm/s')

		 if (DEF_USE_CBL_HEIGHT) then
         ! atmospheric boundary layer height [m]
           call flux_map_and_write_2d ( DEF_hist_vars%xy_hpbl, &
              a_hpbl, f_xy_hpbl, file_hist, 'f_xy_hpbl', itime_in_file, sumarea, filter, &
              'boundary layer height','m')
		 endif

         ! ------------------------------------------------------------------------------------------
         ! Mapping the fluxes and state variables at patch [numpatch] to grid
         ! ------------------------------------------------------------------------------------------
         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype < 99
               vectmp (:) = 1.
               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask
               ENDIF
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! wind stress: E-W [kg/m/s2]
         call flux_map_and_write_2d ( DEF_hist_vars%taux, &
            a_taux, f_taux, file_hist, 'f_taux', itime_in_file, sumarea, filter, &
            'wind stress: E-W','kg/m/s2')

         ! wind stress: N-S [kg/m/s2]
         call flux_map_and_write_2d ( DEF_hist_vars%tauy, &
            a_tauy, f_tauy, file_hist, 'f_tauy', itime_in_file, sumarea, filter, &
            'wind stress: N-S','kg/m/s2')

         ! sensible heat from canopy height to atmosphere [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%fsena, &
            a_fsena, f_fsena, file_hist, 'f_fsena', itime_in_file, sumarea, filter, &
            'sensible heat from canopy height to atmosphere','W/m2')

         ! latent heat flux from canopy height to atmosphere [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%lfevpa, &
            a_lfevpa, f_lfevpa, file_hist, 'f_lfevpa', itime_in_file, sumarea, filter, &
            'latent heat flux from canopy height to atmosphere','W/m2')

         ! evapotranspiration from canopy to atmosphere [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%fevpa, &
            a_fevpa, f_fevpa, file_hist, 'f_fevpa', itime_in_file, sumarea, filter, &
            'evapotranspiration from canopy height to atmosphere','mm/s')

         ! sensible heat from leaves [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%fsenl, &
            a_fsenl, f_fsenl, file_hist, 'f_fsenl', itime_in_file, sumarea, filter, &
            'sensible heat from leaves','W/m2')

         ! evaporation+transpiration from leaves [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%fevpl, &
            a_fevpl, f_fevpl, file_hist, 'f_fevpl', itime_in_file, sumarea, filter, &
            'evaporation+transpiration from leaves','mm/s')

         ! transpiration rate [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%etr, &
            a_etr, f_etr, file_hist, 'f_etr', itime_in_file, sumarea, filter, &
            'transpiration rate','mm/s')

         ! sensible heat flux from ground [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%fseng, &
            a_fseng, f_fseng, file_hist, 'f_fseng', itime_in_file, sumarea, filter, &
            'sensible heat flux from ground','W/m2')

         ! evaporation heat flux from ground [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%fevpg, &
            a_fevpg, f_fevpg, file_hist, 'f_fevpg', itime_in_file, sumarea, filter, &
            'evaporation heat flux from ground','mm/s')

         ! ground heat flux [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%fgrnd, &
            a_fgrnd, f_fgrnd, file_hist, 'f_fgrnd', itime_in_file, sumarea, filter, &
            'ground heat flux','W/m2')

         ! solar absorbed by sunlit canopy [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%sabvsun, &
            a_sabvsun, f_sabvsun, file_hist, 'f_sabvsun', itime_in_file, sumarea, filter, &
            'solar absorbed by sunlit canopy','W/m2')

         ! solar absorbed by shaded [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%sabvsha, &
            a_sabvsha, f_sabvsha, file_hist, 'f_sabvsha', itime_in_file, sumarea, filter, &
            'solar absorbed by shaded','W/m2')

         ! solar absorbed by ground  [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%sabg, &
            a_sabg, f_sabg, file_hist, 'f_sabg', itime_in_file, sumarea, filter, &
            'solar absorbed by ground','W/m2')

         ! outgoing long-wave radiation from ground+canopy [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%olrg, &
            a_olrg, f_olrg, file_hist, 'f_olrg', itime_in_file, sumarea, filter, &
            'outgoing long-wave radiation from ground+canopy','W/m2')

         ! net radiation [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%rnet, &
            a_rnet, f_rnet, file_hist, 'f_rnet', itime_in_file, sumarea, filter, &
            'net radiation','W/m2')

         ! the error of water banace [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%xerr, &
            a_xerr, f_xerr, file_hist, 'f_xerr', itime_in_file, sumarea, filter, &
            'the error of water banace','mm/s')

         ! the error of energy balance [W/m2]
         call flux_map_and_write_2d ( DEF_hist_vars%zerr, &
            a_zerr, f_zerr, file_hist, 'f_zerr', itime_in_file, sumarea, filter, &
            'the error of energy balance','W/m2')

         ! surface runoff [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%rsur, &
            a_rsur, f_rsur, file_hist, 'f_rsur', itime_in_file, sumarea, filter, &
            'surface runoff / surface water change by lateral flow)','mm/s')

         ! subsurface runoff [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%rsub, &
            a_rsub, f_rsub, file_hist, 'f_rsub', itime_in_file, sumarea, filter, &
            'subsurface runoff / groundwater change by lateral flow','mm/s')

         ! total runoff [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%rnof, &
            a_rnof, f_rnof, file_hist, 'f_rnof', itime_in_file, sumarea, filter, &
            'total runoff / total change of surface water and groundwater by lateral flow','mm/s')

         ! interception [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%qintr, &
            a_qintr, f_qintr, file_hist, 'f_qintr', itime_in_file, sumarea, filter, &
            'interception','mm/s')

         ! inflitraton [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%qinfl, &
            a_qinfl, f_qinfl, file_hist, 'f_qinfl', itime_in_file, sumarea, filter, &
            'f_qinfl','mm/s')

         ! throughfall [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%qdrip, &
            a_qdrip, f_qdrip, file_hist, 'f_qdrip', itime_in_file, sumarea, filter, &
            'total throughfall','mm/s')

         ! total water storage [mm]
         call flux_map_and_write_2d ( DEF_hist_vars%wat, &
            a_wat, f_wat, file_hist, 'f_wat', itime_in_file, sumarea, filter, &
            'total water storage','mm')

         ! canopy assimilation rate [mol m-2 s-1]
         call flux_map_and_write_2d ( DEF_hist_vars%assim, &
            a_assim, f_assim, file_hist, 'f_assim', itime_in_file, sumarea, filter, &
            'canopy assimilation rate','mol m-2 s-1')

         ! respiration (plant+soil) [mol m-2 s-1]
         call flux_map_and_write_2d ( DEF_hist_vars%respc, &
            a_respc, f_respc, file_hist, 'f_respc', itime_in_file, sumarea, filter, &
            'respiration (plant+soil)','mol m-2 s-1')

         ! groundwater recharge rate [mm/s]
         call flux_map_and_write_2d ( DEF_hist_vars%qcharge, &
            a_qcharge, f_qcharge, file_hist, 'f_qcharge', itime_in_file, sumarea, filter, &
            'groundwater recharge rate','mm/s')

         ! ground surface temperature [K]
         call flux_map_and_write_2d ( DEF_hist_vars%t_grnd, &
            a_t_grnd, f_t_grnd, file_hist, 'f_t_grnd', itime_in_file, sumarea, filter, &
            'ground surface temperature','K')

         ! leaf temperature [K]
         call flux_map_and_write_2d ( DEF_hist_vars%tleaf, &
            a_tleaf, f_tleaf, file_hist, 'f_tleaf', itime_in_file, sumarea, filter, &
            'leaf temperature','K')

         ! depth of water on foliage [mm]
         call flux_map_and_write_2d ( DEF_hist_vars%ldew, &
            a_ldew, f_ldew, file_hist, 'f_ldew', itime_in_file, sumarea, filter, &
            'depth of water on foliage','mm')

         ! snow cover, water equivalent [mm]
         call flux_map_and_write_2d ( DEF_hist_vars%scv, &
            a_scv, f_scv, file_hist, 'f_scv', itime_in_file, sumarea, filter, &
            'snow cover, water equivalent','mm')

         ! snow depth [meter]
         call flux_map_and_write_2d ( DEF_hist_vars%snowdp, &
            a_snowdp, f_snowdp, file_hist, 'f_snowdp', itime_in_file, sumarea, filter, &
            'snow depth','meter')

         ! fraction of snow cover on ground
         call flux_map_and_write_2d ( DEF_hist_vars%fsno, &
            a_fsno, f_fsno, file_hist, 'f_fsno', itime_in_file, sumarea, filter, &
            'fraction of snow cover on ground','-')

         ! fraction of veg cover, excluding snow-covered veg [-]
         call flux_map_and_write_2d ( DEF_hist_vars%sigf, &
            a_sigf, f_sigf, file_hist, 'f_sigf', itime_in_file, sumarea, filter, &
            'fraction of veg cover, excluding snow-covered veg','-')

         ! leaf greenness
         call flux_map_and_write_2d ( DEF_hist_vars%green, &
            a_green, f_green, file_hist, 'f_green', itime_in_file, sumarea, filter, &
            'leaf greenness','-')

         ! leaf area index
         call flux_map_and_write_2d ( DEF_hist_vars%lai, &
            a_lai, f_lai, file_hist, 'f_lai', itime_in_file, sumarea, filter, &
            'leaf area index','m2/m2')

         ! leaf area index
         call flux_map_and_write_2d ( DEF_hist_vars%laisun, &
            a_laisun, f_laisun, file_hist, 'f_laisun', itime_in_file, sumarea, filter, &
            'sunlit leaf area index','m2/m2')

         ! leaf area index
         call flux_map_and_write_2d ( DEF_hist_vars%laisha, &
            a_laisha, f_laisha, file_hist, 'f_laisha', itime_in_file, sumarea, filter, &
            'shaded leaf area index','m2/m2')

         ! stem area index
         call flux_map_and_write_2d ( DEF_hist_vars%sai, &
            a_sai, f_sai, file_hist, 'f_sai', itime_in_file, sumarea, filter, &
            'stem area index','m2/m2')

         ! averaged albedo [visible, direct; direct, diffuse]
         call flux_map_and_write_4d ( DEF_hist_vars%alb, &
            a_alb, f_alb, file_hist, 'f_alb', 'band', 'rtyp', itime_in_file, sumarea, filter, &
            'averaged albedo direct','%')

         ! averaged bulk surface emissivity
         call flux_map_and_write_2d ( DEF_hist_vars%emis, &
            a_emis, f_emis, file_hist, 'f_emis', itime_in_file, sumarea, filter, &
            'averaged bulk surface emissivity','-')

         ! effective roughness [m]
         call flux_map_and_write_2d ( DEF_hist_vars%z0m, &
            a_z0m, f_z0m, file_hist, 'f_z0m', itime_in_file, sumarea, filter, &
            'effective roughness','m')

         ! radiative temperature of surface [K]
         call flux_map_and_write_2d ( DEF_hist_vars%trad, &
            a_trad, f_trad, file_hist, 'f_trad', itime_in_file, sumarea, filter, &
            'radiative temperature of surface','kelvin')

         ! soil resistance [m/s]
         call flux_map_and_write_2d ( DEF_hist_vars%rss, &
            a_rss, f_rss, file_hist, 'f_rss', itime_in_file, sumarea, filter, &
            'soil resistance','m/s')

         ! 2 m height air temperature [kelvin]
         call flux_map_and_write_2d ( DEF_hist_vars%tref, &
            a_tref, f_tref, file_hist, 'f_tref', itime_in_file, sumarea, filter, &
            '2 m height air temperature','kelvin')

         ! 2 m height air specific humidity [kg/kg]
         call flux_map_and_write_2d ( DEF_hist_vars%qref, &
            a_qref, f_qref, file_hist, 'f_qref', itime_in_file, sumarea, filter, &
            '2 m height air specific humidity','kg/kg')

         ! ------------------------------------------------------------------------------------------
         ! Mapping the urban variables at patch [numurban] to grid
         ! ------------------------------------------------------------------------------------------

#ifdef URBAN_MODEL
         if (p_is_worker) then
            if (numpatch > 0) then
               DO i = 1, numpatch
                  IF (patchtype(i) == 1) THEN
                     u = patch2urban(i)

                     filter_urb(u) = .true.
                     vectmp_urb(u) = 1.

                     IF (DEF_forcing%has_missing_value) THEN
                        filter_urb(u) = filter_urb(u) .and. forcmask(i)
                     ENDIF
                  ENDIF
               ENDDO
            end if
         end if
         
         call mp2g_hist_urb%map (vectmp_urb, sumarea_urb, spv = spval, msk = filter_urb)
         
         ! sensible heat from building roof [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%fsen_roof, &
            a_senroof, f_senroof, file_hist, 'f_fsenroof', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban roof [W/m2]','W/m2')

         ! sensible heat from building sunlit wall [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%fsen_wsun, &
            a_senwsun, f_senwsun, file_hist, 'f_fsenwsun', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban sunlit wall [W/m2]','W/m2')

         ! sensible heat from building shaded wall [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%fsen_wsha, &
            a_senwsha, f_senwsha, file_hist, 'f_fsenwsha', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban shaded wall [W/m2]','W/m2')

         ! sensible heat from impervious ground [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%fsen_gimp, &
            a_sengimp, f_sengimp, file_hist, 'f_fsengimp', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban impervious ground [W/m2]','W/m2')

         ! sensible heat from pervious ground [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%fsen_gper, &
            a_sengper, f_sengper, file_hist, 'f_fsengper', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban pervious ground [W/m2]','W/m2')

         ! sensible heat from urban tree [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%fsen_urbl, &
            a_senurbl, f_senurbl, file_hist, 'f_fsenurbl', itime_in_file, sumarea_urb, filter_urb, &
            'sensible heat from urban tree [W/m2]','W/m2')

         ! latent heat flux from building roof [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%lfevp_roof, &
            a_lfevproof, f_lfevproof, file_hist, 'f_lfevproof', itime_in_file, sumarea_urb, filter_urb, &
            'latent heat from urban roof [W/m2]','W/m2')

         ! latent heat flux from impervious ground [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%lfevp_gimp, &
            a_lfevpgimp, f_lfevpgimp, file_hist, 'f_lfevpgimp', itime_in_file, sumarea_urb, filter_urb, &
            'latent heat from urban impervious ground [W/m2]','W/m2')

         ! latent heat flux from pervious ground [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%lfevp_gper, &
            a_lfevpgper, f_lfevpgper, file_hist, 'f_lfevpgper', itime_in_file, sumarea_urb, filter_urb, &
            'latent heat from urban pervious ground [W/m2]','W/m2')

         ! latent heat flux from urban tree [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%lfevp_urbl, &
            a_lfevpurbl, f_lfevpurbl, file_hist, 'f_lfevpurbl', itime_in_file, sumarea_urb, filter_urb, &
            'latent heat from urban tree [W/m2]','W/m2')

         ! sensible flux from heat or cool AC [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%fhac, &
            a_fhac, f_fhac, file_hist, 'f_fhac', itime_in_file, sumarea_urb, filter_urb, &
            'sensible flux from heat or cool AC [W/m2]','W/m2')

         ! waste heat flux from heat or cool AC [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%fwst, &
            a_fwst, f_fwst, file_hist, 'f_fwst', itime_in_file, sumarea_urb, filter_urb, &
            'waste heat flux from heat or cool AC [W/m2]','W/m2')

         ! flux from inner and outter air exchange [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%fach, &
            a_fach, f_fach, file_hist, 'f_fach', itime_in_file, sumarea_urb, filter_urb, &
            'flux from inner and outter air exchange [W/m2]','W/m2')

         ! flux from total heating/cooling [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%fhah, &
            a_fhah, f_fhah, file_hist, 'f_fhah', itime_in_file, sumarea_urb, filter_urb, &
            'flux from heating/cooling [W/m2]','W/m2')

         ! flux from metabolism [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%meta, &
            a_meta, f_fmeta, file_hist, 'f_fmeta', itime_in_file, sumarea_urb, filter_urb, &
            'flux from human metabolism [W/m2]','W/m2')

         ! flux from vehicle [W/m2]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%vehc, &
            a_vehc, f_fvehc, file_hist, 'f_fvehc', itime_in_file, sumarea_urb, filter_urb, &
            'flux from traffic [W/m2]','W/m2')

         ! temperature of inner building [K]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%t_room, &
            a_t_room, f_t_room, file_hist, 'f_t_room', itime_in_file, sumarea_urb, filter_urb, &
            'temperature of inner building [K]','kelvin')

         ! temperature of outer building [K]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%tafu, &
            a_tafu, f_tafu, file_hist, 'f_tafu', itime_in_file, sumarea_urb, filter_urb, &
            'temperature of outer building [K]','kelvin')

         ! temperature of building roof [K]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%t_roof, &
            a_troof, f_troof, file_hist, 'f_t_roof', itime_in_file, sumarea_urb, filter_urb, &
            'temperature of urban roof [K]','kelvin')

         ! temperature of building wall [K]
         call flux_map_and_write_urb_2d ( DEF_hist_vars%t_wall, &
            a_twall, f_twall, file_hist, 'f_t_wall', itime_in_file, sumarea_urb, filter_urb, &
            'temperature of urban wall [K]','kelvin')
#endif



         ! 1: assimsun enf temperate
         call flux_map_and_write_2d ( DEF_hist_vars%assimsun, &
             a_assimsun, f_assimsun, file_hist, 'f_assimsun', itime_in_file, sumarea, filter, &
             'Photosynthetic assimilation rate of sunlit leaf for needleleaf evergreen temperate tree','mol m-2 s-1')

         ! 1: assimsha enf temperate
         call flux_map_and_write_2d ( DEF_hist_vars%assimsha, &
             a_assimsha, f_assimsha, file_hist, 'f_assimsha', itime_in_file, sumarea, filter, &
             'Photosynthetic assimilation rate of shaded leaf for needleleaf evergreen temperate tree','mol m-2 s-1')

         ! 1: etrsun enf temperate
         call flux_map_and_write_2d ( DEF_hist_vars%etrsun, &
             a_etrsun, f_etrsun, file_hist, 'f_etrsun', itime_in_file, sumarea, filter, &
             'Transpiration rate of sunlit leaf for needleleaf evergreen temperate tree','mm s-1')

         ! 1: etrsha enf temperate
         call flux_map_and_write_2d ( DEF_hist_vars%etrsha, &
             a_etrsha, f_etrsha, file_hist, 'f_etrsha', itime_in_file, sumarea, filter, &
             'Transpiration rate of shaded leaf for needleleaf evergreen temperate tree','mm s-1')

         ! rstfacsun
         call flux_map_and_write_2d ( DEF_hist_vars%rstfacsun, &
             a_rstfacsun, f_rstfacsun, file_hist, 'f_rstfacsun', itime_in_file, sumarea, filter, &
             'Ecosystem level Water stress factor on sunlit canopy','unitless')

         ! rstfacsha
         call flux_map_and_write_2d ( DEF_hist_vars%rstfacsha, &
             a_rstfacsha, f_rstfacsha, file_hist, 'f_rstfacsha', itime_in_file, sumarea, filter, &
             'Ecosystem level Water stress factor on shaded canopy','unitless')

         ! gssun
         call flux_map_and_write_2d ( DEF_hist_vars%gssun, &
             a_gssun, f_gssun, file_hist, 'f_gssun', itime_in_file, sumarea, filter, &
             'Ecosystem level canopy conductance on sunlit canopy','mol m-2 s-1')

         ! gssha
         call flux_map_and_write_2d ( DEF_hist_vars%gssha, &
             a_gssha, f_gssha, file_hist, 'f_gssha', itime_in_file, sumarea, filter, &
             'Ecosystem level canopy conductance on shaded canopy','mol m-2 s-1')

#ifdef BGC
         ! leaf carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafc, &
             a_leafc, f_leafc, file_hist, 'f_leafc', itime_in_file, sumarea, filter, &
             'leaf carbon display pool','gC/m2')

         ! leaf carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_storage, &
             a_leafc_storage, f_leafc_storage, file_hist, 'f_leafc_storage', itime_in_file, sumarea, filter, &
             'leaf carbon storage pool','gC/m2')

         ! leaf carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_xfer, &
             a_leafc_xfer, f_leafc_xfer, file_hist, 'f_leafc_xfer', itime_in_file, sumarea, filter, &
             'leaf carbon transfer pool','gC/m2')

         ! fine root carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootc, &
             a_frootc, f_frootc, file_hist, 'f_frootc', itime_in_file, sumarea, filter, &
             'fine root carbon display pool','gC/m2')

         ! fine root carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootc_storage, &
             a_frootc_storage, f_frootc_storage, file_hist, 'f_frootc_storage', itime_in_file, sumarea, filter, &
             'fine root carbon storage pool','gC/m2')

         ! fine root carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootc_xfer, &
             a_frootc_xfer, f_frootc_xfer, file_hist, 'f_frootc_xfer', itime_in_file, sumarea, filter, &
             'fine root carbon transfer pool','gC/m2')

         ! live stem carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemc, &
             a_livestemc, f_livestemc, file_hist, 'f_livestemc', itime_in_file, sumarea, filter, &
             'live stem carbon display pool','gC/m2')

         ! live stem carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemc_storage, &
             a_livestemc_storage, f_livestemc_storage, file_hist, 'f_livestemc_storage', itime_in_file, sumarea, filter, &
             'live stem carbon storage pool','gC/m2')

         ! live stem carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemc_xfer, &
             a_livestemc_xfer, f_livestemc_xfer, file_hist, 'f_livestemc_xfer', itime_in_file, sumarea, filter, &
             'live stem carbon transfer pool','gC/m2')

         ! dead stem carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemc, &
             a_deadstemc, f_deadstemc, file_hist, 'f_deadstemc', itime_in_file, sumarea, filter, &
             'dead stem carbon display pool','gC/m2')

         ! dead stem carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemc_storage, &
             a_deadstemc_storage, f_deadstemc_storage, file_hist, 'f_deadstemc_storage', itime_in_file, sumarea, filter, &
             'dead stem carbon storage pool','gC/m2')

         ! dead stem carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemc_xfer, &
             a_deadstemc_xfer, f_deadstemc_xfer, file_hist, 'f_deadstemc_xfer', itime_in_file, sumarea, filter, &
             'dead stem carbon transfer pool','gC/m2')

         ! live coarse root carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootc, &
             a_livecrootc, f_livecrootc, file_hist, 'f_livecrootc', itime_in_file, sumarea, filter, &
             'live coarse root carbon display pool','gC/m2')

         ! live coarse root carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootc_storage, &
             a_livecrootc_storage, f_livecrootc_storage, file_hist, 'f_livecrootc_storage', itime_in_file, sumarea, filter, &
             'live coarse root carbon storage pool','gC/m2')

         ! live coarse root carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootc_xfer, &
             a_livecrootc_xfer, f_livecrootc_xfer, file_hist, 'f_livecrootc_xfer', itime_in_file, sumarea, filter, &
             'live coarse root carbon transfer pool','gC/m2')

         ! dead coarse root carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootc, &
             a_deadcrootc, f_deadcrootc, file_hist, 'f_deadcrootc', itime_in_file, sumarea, filter, &
             'dead coarse root carbon display pool','gC/m2')

         ! dead coarse root carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootc_storage, &
             a_deadcrootc_storage, f_deadcrootc_storage, file_hist, 'f_deadcrootc_storage', itime_in_file, sumarea, filter, &
             'dead coarse root carbon storage pool','gC/m2')

         ! dead coarse root carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootc_xfer, &
             a_deadcrootc_xfer, f_deadcrootc_xfer, file_hist, 'f_deadcrootc_xfer', itime_in_file, sumarea, filter, &
             'dead coarse root carbon transfer pool','gC/m2')

#ifdef CROP
         ! grain carbon display pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainc, &
             a_grainc, f_grainc, file_hist, 'f_grainc', itime_in_file, sumarea, filter, &
             'grain carbon display pool','gC/m2')

         ! grain carbon storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainc_storage, &
             a_grainc_storage, f_grainc_storage, file_hist, 'f_grainc_storage', itime_in_file, sumarea, filter, &
             'grain carbon storage pool','gC/m2')

         ! grain carbon transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainc_xfer, &
             a_grainc_xfer, f_grainc_xfer, file_hist, 'f_grainc_xfer', itime_in_file, sumarea, filter, &
             'grain carbon transfer pool','gC/m2')
#endif

         ! leaf nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafn, &
             a_leafn, f_leafn, file_hist, 'f_leafn', itime_in_file, sumarea, filter, &
             'leaf nitrogen display pool','gN/m2')

         ! leaf nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafn_storage, &
             a_leafn_storage, f_leafn_storage, file_hist, 'f_leafn_storage', itime_in_file, sumarea, filter, &
             'leaf nitrogen storage pool','gN/m2')

         ! leaf nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%leafn_xfer, &
             a_leafn_xfer, f_leafn_xfer, file_hist, 'f_leafn_xfer', itime_in_file, sumarea, filter, &
             'leaf nitrogen transfer pool','gN/m2')

         ! fine root nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootn, &
             a_frootn, f_frootn, file_hist, 'f_frootn', itime_in_file, sumarea, filter, &
             'fine root nitrogen display pool','gN/m2')

         ! fine root nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootn_storage, &
             a_frootn_storage, f_frootn_storage, file_hist, 'f_frootn_storage', itime_in_file, sumarea, filter, &
             'fine root nitrogen storage pool','gN/m2')

         ! fine root nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%frootn_xfer, &
             a_frootn_xfer, f_frootn_xfer, file_hist, 'f_frootn_xfer', itime_in_file, sumarea, filter, &
             'fine root nitrogen transfer pool','gN/m2')

         ! live stem nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemn, &
             a_livestemn, f_livestemn, file_hist, 'f_livestemn', itime_in_file, sumarea, filter, &
             'live stem nitrogen display pool','gN/m2')

         ! live stem nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemn_storage, &
             a_livestemn_storage, f_livestemn_storage, file_hist, 'f_livestemn_storage', itime_in_file, sumarea, filter, &
             'live stem nitrogen storage pool','gN/m2')

         ! live stem nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%livestemn_xfer, &
             a_livestemn_xfer, f_livestemn_xfer, file_hist, 'f_livestemn_xfer', itime_in_file, sumarea, filter, &
             'live stem nitrogen transfer pool','gN/m2')

         ! dead stem nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemn, &
             a_deadstemn, f_deadstemn, file_hist, 'f_deadstemn', itime_in_file, sumarea, filter, &
             'dead stem nitrogen display pool','gN/m2')

         ! dead stem nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemn_storage, &
             a_deadstemn_storage, f_deadstemn_storage, file_hist, 'f_deadstemn_storage', itime_in_file, sumarea, filter, &
             'dead stem nitrogen storage pool','gN/m2')

         ! dead stem nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadstemn_xfer, &
             a_deadstemn_xfer, f_deadstemn_xfer, file_hist, 'f_deadstemn_xfer', itime_in_file, sumarea, filter, &
             'dead stem nitrogen transfer pool','gN/m2')

         ! live coarse root nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootn, &
             a_livecrootn, f_livecrootn, file_hist, 'f_livecrootn', itime_in_file, sumarea, filter, &
             'live coarse root nitrogen display pool','gN/m2')

         ! live coarse root nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootn_storage, &
             a_livecrootn_storage, f_livecrootn_storage, file_hist, 'f_livecrootn_storage', itime_in_file, sumarea, filter, &
             'live coarse root nitrogen storage pool','gN/m2')

         ! live coarse root nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%livecrootn_xfer, &
             a_livecrootn_xfer, f_livecrootn_xfer, file_hist, 'f_livecrootn_xfer', itime_in_file, sumarea, filter, &
             'live coarse root nitrogen transfer pool','gN/m2')

         ! dead coarse root nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootn, &
             a_deadcrootn, f_deadcrootn, file_hist, 'f_deadcrootn', itime_in_file, sumarea, filter, &
             'dead coarse root nitrogen display pool','gN/m2')

         ! dead coarse root nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootn_storage, &
             a_deadcrootn_storage, f_deadcrootn_storage, file_hist, 'f_deadcrootn_storage', itime_in_file, sumarea, filter, &
             'dead coarse root nitrogen storage pool','gN/m2')

         ! dead coarse root nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%deadcrootn_xfer, &
             a_deadcrootn_xfer, f_deadcrootn_xfer, file_hist, 'f_deadcrootn_xfer', itime_in_file, sumarea, filter, &
             'dead coarse root nitrogen transfer pool','gN/m2')

#ifdef CROP
         ! grain nitrogen display pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainn, &
             a_grainn, f_grainn, file_hist, 'f_grainn', itime_in_file, sumarea, filter, &
             'grain nitrogen display pool','gN/m2')

         ! grain nitrogen storage pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainn_storage, &
             a_grainn_storage, f_grainn_storage, file_hist, 'f_grainn_storage', itime_in_file, sumarea, filter, &
             'grain nitrogen storage pool','gN/m2')

         ! grain nitrogen transfer pool
         call flux_map_and_write_2d ( DEF_hist_vars%grainn_xfer, &
             a_grainn_xfer, f_grainn_xfer, file_hist, 'f_grainn_xfer', itime_in_file, sumarea, filter, &
             'grain nitrogen transfer pool','gN/m2')
#endif

         ! retranslocation nitrogen pool
         call flux_map_and_write_2d ( DEF_hist_vars%retrasn, &
             a_retransn, f_retransn, file_hist, 'f_retrasn', itime_in_file, sumarea, filter, &
             'retranslocation nitrogen pool','gN/m2')

         ! gross primary productivity
         call flux_map_and_write_2d ( DEF_hist_vars%gpp, &
             a_gpp, f_gpp, file_hist, 'f_gpp', itime_in_file, sumarea, filter, &
             'gross primary productivity','gC/m2/s')

         ! gross primary productivity
         call flux_map_and_write_2d ( DEF_hist_vars%downreg, &
             a_downreg, f_downreg, file_hist, 'f_downreg', itime_in_file, sumarea, filter, &
             'gpp downregulation due to N limitation','unitless')

         call flux_map_and_write_2d ( DEF_hist_vars%fpg, &
             a_fpg, f_fpg, file_hist, 'f_fpg', itime_in_file, sumarea, filter, &
             'fraction of gpp potential','unitless')

         call flux_map_and_write_2d ( DEF_hist_vars%fpi, &
             a_fpi, f_fpi, file_hist, 'f_fpi', itime_in_file, sumarea, filter, &
             'fraction of immobalization','unitless')

         ! autotrophic respiration
         call flux_map_and_write_2d ( DEF_hist_vars%ar , &
             a_ar , f_ar , file_hist, 'f_ar', itime_in_file, sumarea, filter, &
             'autotrophic respiration','gC/m2/s')

         ! CWD production
         call flux_map_and_write_2d ( DEF_hist_vars%cwdprod , &
             a_cwdprod , f_cwdprod , file_hist, 'f_cwdprod', itime_in_file, sumarea, filter, &
             'CWD production','gC/m2/s')

         ! CWD decomposition
         call flux_map_and_write_2d ( DEF_hist_vars%cwddecomp , &
             a_cwddecomp , f_cwddecomp , file_hist, 'f_cwddecomp', itime_in_file, sumarea, filter, &
             'CWD decomposition','gC/m2/s')

         ! heterotrophic respiration
         call flux_map_and_write_2d ( DEF_hist_vars%hr , &
             a_hr , f_hr , file_hist, 'f_hr', itime_in_file, sumarea, filter, &
             'heterotrophic respiration','gC/m2/s')

#ifdef CROP
         ! crop phase
         call flux_map_and_write_2d ( DEF_hist_vars%cphase, &
             a_cphase, f_cphase, file_hist, 'f_cphase', itime_in_file, sumarea, filter, &
             'crop phase','unitless')

        ! heat unit index
        if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_hui (:)
            end if
         end if

        call flux_map_and_write_2d ( DEF_hist_vars%hui, &
             vecacc, f_hui, file_hist, 'f_hui', itime_in_file, sumarea, filter, &
             'heat unit index','unitless')

         ! gdd needed to harvest
        call flux_map_and_write_2d ( DEF_hist_vars%gddmaturity, &
             a_gddmaturity, f_gddmaturity, file_hist, 'f_gddmaturity', itime_in_file, sumarea, filter, &
             'gdd needed to harvest','ddays')

        ! gdd past planting date for crop
        call flux_map_and_write_2d ( DEF_hist_vars%gddplant, &
             a_gddplant, f_gddplant, file_hist, 'f_gddplant', itime_in_file, sumarea, filter, &
             'gdd past planting date for crop','ddays')

        ! vernalization response
        call flux_map_and_write_2d ( DEF_hist_vars%vf, &
             a_vf, f_vf, file_hist, 'f_vf', itime_in_file, sumarea, filter, &
             'vernalization response', 'unitless')

         ! 1-yr crop production carbon
         call flux_map_and_write_2d ( DEF_hist_vars%cropprod1c, &
             a_cropprod1c, f_cropprod1c, file_hist, 'f_cropprod1c', itime_in_file, sumarea, filter, &
             '1-yr crop production carbon','gC/m2')

         ! loss rate of 1-yr crop production carbon
         call flux_map_and_write_2d ( DEF_hist_vars%cropprod1c_loss, &
             a_cropprod1c_loss, f_cropprod1c_loss, file_hist, 'f_cropprod1c_loss', itime_in_file, sumarea, filter, &
             'loss rate of 1-yr crop production carbon','gC/m2/s')

         ! crop seed deficit
         call flux_map_and_write_2d ( DEF_hist_vars%cropseedc_deficit, &
             a_cropseedc_deficit, f_cropseedc_deficit, file_hist, 'f_cropseedc_deficit', itime_in_file, sumarea, filter, &
             'crop seed deficit','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         ! grain to crop production carbon
         call flux_map_and_write_2d ( DEF_hist_vars%grainc_to_cropprodc, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_grainc_to_cropprodc', itime_in_file, sumarea, filter, &
             'grain to crop production carbon','gC/m2/s')

         ! grain to crop seed carbon
         call flux_map_and_write_2d ( DEF_hist_vars%grainc_to_seed, &
             a_grainc_to_seed, f_grainc_to_seed, file_hist, 'f_grainc_to_seed', itime_in_file, sumarea, filter, &
             'grain to crop seed carbon','gC/m2/s')

         ! grain to crop seed carbon
         call flux_map_and_write_2d ( DEF_hist_vars%fert_to_sminn, &
             a_fert_to_sminn, f_fert_to_sminn, file_hist, 'f_fert_to_sminn', itime_in_file, sumarea, filter, &
             'fertilization','gN/m2/s')

#endif

         ! grain to crop seed carbon
         call flux_map_and_write_2d ( DEF_hist_vars%ndep_to_sminn, &
             a_ndep_to_sminn, f_ndep_to_sminn, file_hist, 'f_ndep_to_sminn', itime_in_file, sumarea, filter, &
             'nitrogen deposition','gN/m2/s')

         IF(DEF_USE_OZONESTRESS)THEN
         ! ozone concentration
            call flux_map_and_write_2d ( DEF_hist_vars%xy_ozone, &
               a_ozone, f_xy_ozone, file_hist, 'f_xy_ozone', itime_in_file, sumarea, filter, &
               'Ozone concentration','mol/mol')
         ENDIF

         ! litter 1 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr1c_vr, &
            a_litr1c_vr, f_litr1c_vr, file_hist, 'f_litr1c_vr', 'soil', &
            itime_in_file, sumarea, filter,'litter 1 carbon density in soil layers','gC/m3')

         ! litter 2 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr2c_vr, &
            a_litr2c_vr, f_litr2c_vr, file_hist, 'f_litr2c_vr', 'soil', &
            itime_in_file, sumarea, filter,'litter 2 carbon density in soil layers','gC/m3')

         ! litter 3 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr3c_vr, &
            a_litr3c_vr, f_litr3c_vr, file_hist, 'f_litr3c_vr', 'soil', &
            itime_in_file, sumarea, filter,'litter 3 carbon density in soil layers','gC/m3')

         ! soil 1 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil1c_vr, &
            a_soil1c_vr, f_soil1c_vr, file_hist, 'f_soil1c_vr', 'soil', &
            itime_in_file, sumarea, filter,'soil 1 carbon density in soil layers','gC/m3')

         ! soil 2 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil2c_vr, &
            a_soil2c_vr, f_soil2c_vr, file_hist, 'f_soil2c_vr', 'soil', &
            itime_in_file, sumarea, filter,'soil 2 carbon density in soil layers','gC/m3')

         ! soil 3 carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil3c_vr, &
            a_soil3c_vr, f_soil3c_vr, file_hist, 'f_soil3c_vr', 'soil', &
            itime_in_file, sumarea, filter,'soil 3 carbon density in soil layers','gC/m3')

         ! coarse woody debris carbon density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%cwdc_vr, &
            a_cwdc_vr, f_cwdc_vr, file_hist, 'f_cwdc_vr', 'soil', &
            itime_in_file, sumarea, filter,'coarse woody debris carbon density in soil layers','gC/m3')

         ! litter 1 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr1n_vr, &
            a_litr1n_vr, f_litr1n_vr, file_hist, 'f_litr1n_vr', 'soil', &
            itime_in_file, sumarea, filter,'litter 1 nitrogen density in soil layers','gN/m3')

         ! litter 2 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr2n_vr, &
            a_litr2n_vr, f_litr2n_vr, file_hist, 'f_litr2n_vr', 'soil', &
            itime_in_file, sumarea, filter,'litter 2 nitrogen density in soil layers','gN/m3')

         ! litter 3 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%litr3n_vr, &
            a_litr3n_vr, f_litr3n_vr, file_hist, 'f_litr3n_vr', 'soil', &
            itime_in_file, sumarea, filter,'litter 3 nitrogen density in soil layers','gN/m3')

         ! soil 1 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil1n_vr, &
            a_soil1n_vr, f_soil1n_vr, file_hist, 'f_soil1n_vr', 'soil', &
            itime_in_file, sumarea, filter,'soil 1 nitrogen density in soil layers','gN/m3')

         ! soil 2 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil2n_vr, &
            a_soil2n_vr, f_soil2n_vr, file_hist, 'f_soil2n_vr', 'soil', &
            itime_in_file, sumarea, filter,'soil 2 nitrogen density in soil layers','gN/m3')

         ! soil 3 nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%soil3n_vr, &
            a_soil3n_vr, f_soil3n_vr, file_hist, 'f_soil3n_vr', 'soil', &
            itime_in_file, sumarea, filter,'soil 3 nitrogen density in soil layers','gN/m3')

         ! coarse woody debris nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%cwdn_vr, &
            a_cwdn_vr, f_cwdn_vr, file_hist, 'f_cwdn_vr', 'soil', &
            itime_in_file, sumarea, filter,'coarse woody debris nitrogen density in soil layers','gN/m3')

         ! mineral nitrogen density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%sminn_vr, &
            a_sminn_vr, f_sminn_vr, file_hist, 'f_sminn_vr', 'soil', &
            itime_in_file, sumarea, filter,'mineral nitrogen density in soil layers','gN/m3')

         ! bulk density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%BD_all, &
            a_BD_all, f_BD_all, file_hist, 'f_BD_all', 'soil', &
            itime_in_file, sumarea, filter,'bulk density in soil layers','kg/m3')

         ! field capacity in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%wfc, &
            a_wfc, f_wfc, file_hist, 'f_wfc', 'soil', &
            itime_in_file, sumarea, filter,'field capacity in soil layers','m3/m3')

         ! organic matter density in soil layers
         call flux_map_and_write_3d ( DEF_hist_vars%OM_density, &
            a_OM_density, f_OM_density, file_hist, 'f_OM_density', 'soil', &
            itime_in_file, sumarea, filter,'organic matter density in soil layers','kg/m3')

         if(DEF_USE_NITRIF)then
            ! O2 soil Concentration for non-inundated area
            call flux_map_and_write_3d ( DEF_hist_vars%CONC_O2_UNSAT, &
               a_conc_o2_unsat, f_conc_o2_unsat, file_hist, 'f_CONC_O2_UNSAT', 'soil', &
               itime_in_file, sumarea, filter,'O2 soil Concentration for non-inundated area','mol/m3')

            ! O2 consumption from HR and AR for non-inundated area
            call flux_map_and_write_3d ( DEF_hist_vars%O2_DECOMP_DEPTH_UNSAT, &
               a_o2_decomp_depth_unsat, f_o2_decomp_depth_unsat, file_hist, 'f_O2_DECOMP_DEPTH_UNSAT', 'soil', &
               itime_in_file, sumarea, filter,'O2 consumption from HR and AR for non-inundated area','mol/m3/s')
         end if

         if(DEF_USE_FIRE)then
            call flux_map_and_write_2d ( DEF_hist_vars%abm, &
                 vecacc, f_abm, file_hist, 'f_abm', itime_in_file, sumarea, filter, &
                 'peak crop fire month','unitless')

            call flux_map_and_write_2d ( DEF_hist_vars%gdp, &
                 vecacc, f_gdp, file_hist, 'f_gdp', itime_in_file, sumarea, filter, &
                 'gdp','unitless')

            call flux_map_and_write_2d ( DEF_hist_vars%peatf, &
                 vecacc, f_peatf, file_hist, 'f_peatf', itime_in_file, sumarea, filter, &
                 'peatf','unitless')

            call flux_map_and_write_2d ( DEF_hist_vars%hdm, &
                 vecacc, f_hdm, file_hist, 'f_hdm', itime_in_file, sumarea, filter, &
                 'hdm','unitless')

            call flux_map_and_write_2d ( DEF_hist_vars%lnfm, &
                 vecacc, f_lnfm, file_hist, 'f_lnfm', itime_in_file, sumarea, filter, &
                 'lnfm','unitless')
         end if


         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) .ne. 12 .and. patchtype(i) .eq. 0)then
                     filter(i) = .true.
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! 1: gpp enf temperate
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_enftemp, &
             a_gpp_enftemp, f_gpp_enftemp, file_hist, 'f_gpp_enftemp', itime_in_file, sumarea, filter, &
             'gross primary productivity for needleleaf evergreen temperate tree','gC/m2/s')

         ! 1: leaf carbon display pool enf temperate
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_enftemp, &
             a_leafc_enftemp, f_leafc_enftemp, file_hist, 'f_leafc_enftemp', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for needleleaf evergreen temperate tree','gC/m2')

         ! 2: gpp enf boreal
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_enfboreal, &
             a_gpp_enfboreal, f_gpp_enfboreal, file_hist, 'f_gpp_enfboreal', itime_in_file, sumarea, filter, &
             'gross primary productivity for needleleaf evergreen boreal tree','gC/m2/s')

         ! 2: leaf carbon display pool enf boreal
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_enfboreal, &
             a_leafc_enfboreal, f_leafc_enfboreal, file_hist, 'f_leafc_enfboreal', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for needleleaf evergreen boreal tree','gC/m2')

         ! 3: gpp dnf boreal
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_dnfboreal, &
             a_gpp_dnfboreal, f_gpp_dnfboreal, file_hist, 'f_gpp_dnfboreal', itime_in_file, sumarea, filter, &
             'gross primary productivity for needleleaf deciduous boreal tree','gC/m2/s')

         ! 3: leaf carbon display pool dnf boreal
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_dnfboreal, &
             a_leafc_dnfboreal, f_leafc_dnfboreal, file_hist, 'f_leafc_dnfboreal', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for needleleaf deciduous boreal tree','gC/m2')

         ! 4: gpp ebf trop
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_ebftrop, &
             a_gpp_ebftrop, f_gpp_ebftrop, file_hist, 'f_gpp_ebftrop', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf evergreen tropical tree','gC/m2/s')

         ! 4: leaf carbon display pool ebf trop
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_ebftrop, &
             a_leafc_ebftrop, f_leafc_ebftrop, file_hist, 'f_leafc_ebftrop', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf evergreen tropical tree','gC/m2')

         ! 5: gpp ebf temp
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_ebftemp, &
             a_gpp_ebftemp, f_gpp_ebftemp, file_hist, 'f_gpp_ebftemp', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf evergreen temperate tree','gC/m2/s')

         ! 5: leaf carbon display pool ebf temp
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_ebftemp, &
             a_leafc_ebftemp, f_leafc_ebftemp, file_hist, 'f_leafc_ebftemp', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf evergreen temperate tree','gC/m2')

         ! 6: gpp dbf trop
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_dbftrop, &
             a_gpp_dbftrop, f_gpp_dbftrop, file_hist, 'f_gpp_dbftrop', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf deciduous tropical tree','gC/m2/s')

         ! 6: leaf carbon display pool dbf trop
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_dbftrop, &
             a_leafc_dbftrop, f_leafc_dbftrop, file_hist, 'f_leafc_dbftrop', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf deciduous tropical tree','gC/m2')

         ! 7: gpp dbf temp
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_dbftemp, &
             a_gpp_dbftemp, f_gpp_dbftemp, file_hist, 'f_gpp_dbftemp', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf deciduous temperate tree','gC/m2/s')

         ! 7: leaf carbon display pool dbf temp
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_dbftemp, &
             a_leafc_dbftemp, f_leafc_dbftemp, file_hist, 'f_leafc_dbftemp', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf deciduous temperate tree','gC/m2')

         ! 8: gpp dbf boreal
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_dbfboreal, &
             a_gpp_dbfboreal, f_gpp_dbfboreal, file_hist, 'f_gpp_dbfboreal', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf deciduous boreal tree','gC/m2/s')

         ! 8: leaf carbon display pool dbf boreal
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_dbfboreal, &
             a_leafc_dbfboreal, f_leafc_dbfboreal, file_hist, 'f_leafc_dbfboreal', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf deciduous boreal tree','gC/m2')

         ! 9: gpp ebs temp
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_ebstemp, &
             a_gpp_ebstemp, f_gpp_ebstemp, file_hist, 'f_gpp_ebstemp', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf evergreen temperate shrub','gC/m2/s')

         ! 9: leaf carbon display pool ebs temp
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_ebstemp, &
             a_leafc_ebstemp, f_leafc_ebstemp, file_hist, 'f_leafc_ebstemp', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf evergreen temperate shrub','gC/m2')

         ! 10: gpp dbs temp
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_dbstemp, &
             a_gpp_dbstemp, f_gpp_dbstemp, file_hist, 'f_gpp_dbstemp', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf deciduous temperate shrub','gC/m2/s')

         ! 10: leaf carbon display pool dbs temp
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_dbstemp, &
             a_leafc_dbstemp, f_leafc_dbstemp, file_hist, 'f_leafc_dbstemp', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf deciduous temperate shrub','gC/m2')

         ! 11: gpp dbs boreal
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_dbsboreal, &
             a_gpp_dbsboreal, f_gpp_dbsboreal, file_hist, 'f_gpp_dbsboreal', itime_in_file, sumarea, filter, &
             'gross primary productivity for broadleaf deciduous boreal shrub','gC/m2/s')

         ! 11: leaf carbon display pool dbs boreal
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_dbsboreal, &
             a_leafc_dbsboreal, f_leafc_dbsboreal, file_hist, 'f_leafc_dbsboreal', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for broadleaf deciduous boreal shrub','gC/m2')

         ! 12: gpp arctic c3 grass
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_c3arcgrass, &
             a_gpp_c3arcgrass, f_gpp_c3arcgrass, file_hist, 'f_gpp_c3arcgrass', itime_in_file, sumarea, filter, &
             'gross primary productivity for c3 arctic grass','gC/m2/s')

         ! 12: leaf carbon display pool c3 grass
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_c3grass, &
             a_leafc_c3grass, f_leafc_c3grass, file_hist, 'f_leafc_c3grass', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for c3 grass','gC/m2')

         ! 13: gpp c3 grass
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_c3grass, &
             a_gpp_c3grass, f_gpp_c3grass, file_hist, 'f_gpp_c3grass', itime_in_file, sumarea, filter, &
             'gross primary productivity for c3 grass','gC/m2/s')

         ! 13: leaf carbon display pool arctic c3 grass
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_c3grass, &
             a_leafc_c3grass, f_leafc_c3grass, file_hist, 'f_leafc_c3grass', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for c3 arctic grass','gC/m2')

         ! 14: gpp c4 grass
         call flux_map_and_write_2d ( DEF_hist_vars%gpp_c4grass, &
             a_gpp_c4grass, f_gpp_c4grass, file_hist, 'f_gpp_c4grass', itime_in_file, sumarea, filter, &
             'gross primary productivity for c4 grass','gC/m2/s')

         ! 14: leaf carbon display pool arctic c4 grass
         call flux_map_and_write_2d ( DEF_hist_vars%leafc_c4grass, &
             a_leafc_c4grass, f_leafc_c4grass, file_hist, 'f_leafc_c4grass', itime_in_file, sumarea, filter, &
             'leaf carbon display pool for c4 arctic grass','gC/m2')

#ifdef CROP
!*****************************************
         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                      if(pftclass(patch_pft_s(i)) .eq. 19 .or. pftclass(patch_pft_s(i)) .eq. 20)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if
         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_hui (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%huiswheat, &
             vecacc, f_hui, file_hist, 'f_huiswheat', itime_in_file, sumarea, filter, &
             'heat unit index  (rainfed spring wheat)','unitless')

!************************************************************
         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 17 .or. pftclass(patch_pft_s(i)) .eq. 18 &
                   .or. pftclass(patch_pft_s(i)) .eq. 75 .or. pftclass(patch_pft_s(i)) .eq. 76)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%pdcorn, &
            a_pdcorn, f_pdcorn, file_hist, 'f_pdcorn', &
            itime_in_file, sumarea, filter,'planting date of corn','day')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 19 .or. pftclass(patch_pft_s(i)) .eq. 20)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%pdswheat, &
            a_pdswheat, f_pdswheat, file_hist, 'f_pdswheat', &
            itime_in_file, sumarea, filter,'planting date of spring wheat','day')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 21 .or. pftclass(patch_pft_s(i)) .eq. 22)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%pdwwheat, &
            a_pdwwheat, f_pdwwheat, file_hist, 'f_pdwwheat', &
            itime_in_file, sumarea, filter,'planting date of winter wheat','day')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 23 .or. pftclass(patch_pft_s(i)) .eq. 24 &
                   .or. pftclass(patch_pft_s(i)) .eq. 77 .or. pftclass(patch_pft_s(i)) .eq. 78)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%pdsoybean, &
            a_pdsoybean, f_pdsoybean, file_hist, 'f_pdsoybean', &
            itime_in_file, sumarea, filter,'planting date of soybean','day')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 41 .or. pftclass(patch_pft_s(i)) .eq. 42)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%pdcotton, &
            a_pdcotton, f_pdcotton, file_hist, 'f_pdcotton', &
            itime_in_file, sumarea, filter,'planting date of cotton','day')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 61 .or. pftclass(patch_pft_s(i)) .eq. 62)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%pdrice1, &
            a_pdrice1, f_pdrice1, file_hist, 'f_pdrice1', &
            itime_in_file, sumarea, filter,'planting date of rice1','day')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 67 .or. pftclass(patch_pft_s(i)) .eq. 68)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%pdrice2, &
            a_pdrice2, f_pdrice2, file_hist, 'f_pdrice2', &
            itime_in_file, sumarea, filter,'planting date of rice2','day')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 67 .or. pftclass(patch_pft_s(i)) .eq. 68)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%pdsugarcane, &
            a_pdsugarcane, f_pdsugarcane, file_hist, 'f_pdsugarcane', &
            itime_in_file, sumarea, filter,'planting date of sugarcane','day')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 17 .or. pftclass(patch_pft_s(i)) .eq. 18 &
                   .or. pftclass(patch_pft_s(i)) .eq. 75 .or. pftclass(patch_pft_s(i)) .eq. 76)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%fertnitro_corn, &
            a_fertnitro_corn, f_fertnitro_corn, file_hist, 'f_fertnitro_corn', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for corn','gN/m2/yr')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 19 .or. pftclass(patch_pft_s(i)) .eq. 20)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%fertnitro_swheat, &
            a_fertnitro_swheat, f_fertnitro_swheat, file_hist, 'f_fertnitro_swheat', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for spring wheat','gN/m2/yr')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 21 .or. pftclass(patch_pft_s(i)) .eq. 22)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%fertnitro_wwheat, &
            a_fertnitro_wwheat, f_fertnitro_wwheat, file_hist, 'f_fertnitro_wwheat', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for winter wheat','gN/m2/yr')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 23 .or. pftclass(patch_pft_s(i)) .eq. 24 &
                   .or. pftclass(patch_pft_s(i)) .eq. 77 .or. pftclass(patch_pft_s(i)) .eq. 78)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%fertnitro_soybean, &
            a_fertnitro_soybean, f_fertnitro_soybean, file_hist, 'f_fertnitro_soybean', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for soybean','gN/m2/yr')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 41 .or. pftclass(patch_pft_s(i)) .eq. 42)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%fertnitro_cotton, &
            a_fertnitro_cotton, f_fertnitro_cotton, file_hist, 'f_fertnitro_cotton', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for cotton','gN/m2/yr')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 61 .or. pftclass(patch_pft_s(i)) .eq. 62)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%fertnitro_rice1, &
            a_fertnitro_rice1, f_fertnitro_rice1, file_hist, 'f_fertnitro_rice1', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for rice1','gN/m2/yr')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 61 .or. pftclass(patch_pft_s(i)) .eq. 62)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%fertnitro_rice2, &
            a_fertnitro_rice2, f_fertnitro_rice2, file_hist, 'f_fertnitro_rice2', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for rice2','gN/m2/yr')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 67 .or. pftclass(patch_pft_s(i)) .eq. 68)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         call flux_map_and_write_2d ( DEF_hist_vars%fertnitro_sugarcane, &
            a_fertnitro_sugarcane, f_fertnitro_sugarcane, file_hist, 'f_fertnitro_sugarcane', &
            itime_in_file, sumarea, filter,'nitrogen fertilizer for sugarcane','gN/m2/yr')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 17)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to corn production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_rainfed_temp_corn, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_rainfed_temp_corn', itime_in_file, sumarea, filter, &
             'Crop production (rainfed temperate corn)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 18)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to corn production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_irrigated_temp_corn, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_irrigated_temp_corn', itime_in_file, sumarea, filter, &
             'Crop production (irrigated temperate corn)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 19)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if


         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to spring wheat production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_rainfed_spwheat, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_rainfed_spwheat', itime_in_file, sumarea, filter, &
             'Crop production (rainfed spring wheat)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 20)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if


         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to spring wheat production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_irrigated_spwheat, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_irrigated_spwheat', itime_in_file, sumarea, filter, &
             'Crop production (irrigated spring wheat)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 21)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to winter wheat production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_rainfed_wtwheat, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_rainfed_wtwheat', itime_in_file, sumarea, filter, &
             'Crop production (rainfed winter wheat)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 22)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to winter wheat production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_irrigated_wtwheat, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_irrigated_wtwheat', itime_in_file, sumarea, filter, &
             'Crop production (irrigated winter wheat)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 23)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if
         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to soybean production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_rainfed_temp_soybean, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_rainfed_temp_soybean', itime_in_file, sumarea, filter, &
             'Crop production (rainfed temperate soybean)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 24)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if
         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to soybean production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_irrigated_temp_soybean, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_irrigated_temp_soybean', itime_in_file, sumarea, filter, &
             'Crop production (irrigated temperate soybean)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 41)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to cotton production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_rainfed_cotton, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_rainfed_cotton', itime_in_file, sumarea, filter, &
             'Crop production (rainfed cotton)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 42)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to cotton production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_irrigated_cotton, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_irrigated_cotton', itime_in_file, sumarea, filter, &
             'Crop production (irrigated cotton)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 61)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to rice production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_rainfed_rice, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_rainfed_rice', itime_in_file, sumarea, filter, &
             'Crop production (rainfed rice)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 62)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to rice production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_irrigated_rice, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_irrigated_rice', itime_in_file, sumarea, filter, &
             'Crop production (irrigated rice)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 67)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_rainfed_sugarcane, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_rainfed_sugarcane', itime_in_file, sumarea, filter, &
             'Crop production (rainfed sugarcane)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 68)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_irrigated_sugarcane, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_irrigated_sugarcane', itime_in_file, sumarea, filter, &
             'Crop production (irrigated sugarcane)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 75)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_rainfed_trop_corn, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_rainfed_trop_corn', itime_in_file, sumarea, filter, &
             'Crop production (rainfed_trop_corn)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 76)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_irrigated_trop_corn, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_irrigated_trop_corn', itime_in_file, sumarea, filter, &
             'Crop production (irrigated_trop_corn)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 77)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_rainfed_trop_soybean, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_rainfed_trop_soybean', itime_in_file, sumarea, filter, &
             'Crop production (rainfed trop soybean)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 78)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_irrigated_trop_soybean, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_irrigated_trop_soybean', itime_in_file, sumarea, filter, &
             'Crop production (irrigated trop soybean)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 15)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if


         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to unmanaged crop production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_plantdate (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%plantdate_unmanagedcrop, &
             vecacc, f_plantdate, file_hist, 'f_plantdate_unmanagedcrop', itime_in_file, sumarea, filter, &
             'Crop production (unmanaged crop production)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 17)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to corn production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_rainfed_temp_corn, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_rainfed_temp_corn', itime_in_file, sumarea, filter, &
             'Crop production (rainfed temperate corn)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 18)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to corn production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_irrigated_temp_corn, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_irrigated_temp_corn', itime_in_file, sumarea, filter, &
             'Crop production (irrigated temperate corn)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 19)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if


         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to spring wheat production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_rainfed_spwheat, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_rainfed_spwheat', itime_in_file, sumarea, filter, &
             'Crop production (rainfed spring wheat)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 20)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if


         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to spring wheat production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_irrigated_spwheat, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_irrigated_spwheat', itime_in_file, sumarea, filter, &
             'Crop production (irrigated spring wheat)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 21)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to winter wheat production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_rainfed_wtwheat, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_rainfed_wtwheat', itime_in_file, sumarea, filter, &
             'Crop production (rainfed winter wheat)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 22)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to winter wheat production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_irrigated_wtwheat, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_irrigated_wtwheat', itime_in_file, sumarea, filter, &
             'Crop production (irrigated winter wheat)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 23)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if
         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to soybean production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_rainfed_temp_soybean, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_rainfed_temp_soybean', itime_in_file, sumarea, filter, &
             'Crop production (rainfed temperate soybean)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 24)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if
         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to soybean production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_irrigated_temp_soybean, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_irrigated_temp_soybean', itime_in_file, sumarea, filter, &
             'Crop production (irrigated temperate soybean)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 41)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to cotton production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_rainfed_cotton, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_rainfed_cotton', itime_in_file, sumarea, filter, &
             'Crop production (rainfed cotton)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 42)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to cotton production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_irrigated_cotton, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_irrigated_cotton', itime_in_file, sumarea, filter, &
             'Crop production (irrigated cotton)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 61)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to rice production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_rainfed_rice, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_rainfed_rice', itime_in_file, sumarea, filter, &
             'Crop production (rainfed rice)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 62)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to rice production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_irrigated_rice, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_irrigated_rice', itime_in_file, sumarea, filter, &
             'Crop production (irrigated rice)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 67)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_rainfed_sugarcane, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_rainfed_sugarcane', itime_in_file, sumarea, filter, &
             'Crop production (rainfed sugarcane)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 68)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_irrigated_sugarcane, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_irrigated_sugarcane', itime_in_file, sumarea, filter, &
             'Crop production (irrigated sugarcane)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 75)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_rainfed_trop_corn, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_rainfed_trop_corn', itime_in_file, sumarea, filter, &
             'Crop production (rainfed_trop_corn)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 76)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_irrigated_trop_corn, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_irrigated_trop_corn', itime_in_file, sumarea, filter, &
             'Crop production (irrigated_trop_corn)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 77)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_rainfed_trop_soybean, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_rainfed_trop_soybean', itime_in_file, sumarea, filter, &
             'Crop production (rainfed trop soybean)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 78)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to sugarcane production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_irrigated_trop_soybean, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_irrigated_trop_soybean', itime_in_file, sumarea, filter, &
             'Crop production (irrigated trop soybean)','gC/m2/s')

         if (p_is_worker) then
            if (numpatch > 0) then
               do i=1,numpatch
                  if(patchclass(i) == 12)then
                     if(pftclass(patch_pft_s(i)) .eq. 15)then
                        filter(i) = .true.
                     else
                        filter(i) = .false.
                     end if
                  else
                     filter(i) = .false.
                  end if
                  vectmp(i) = 1.
               end do
            end if
         end if


         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! grain to unmanaged crop production carbon
         if (p_is_worker) then
            if (numpatch > 0) then
               vecacc (:) = a_grainc_to_cropprodc (:)
            end if
         end if
         call flux_map_and_write_2d ( DEF_hist_vars%cropprodc_unmanagedcrop, &
             vecacc, f_grainc_to_cropprodc, file_hist, 'f_cropprodc_unmanagedcrop', itime_in_file, sumarea, filter, &
             'Crop production (unmanaged crop production)','gC/m2/s')
#endif
#endif
         ! --------------------------------------------------------------------
         ! Temperature and water (excluding land water bodies and ocean patches)
         ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3; land water bodies => 4; ocean => 99]
         ! --------------------------------------------------------------------

         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype <= 3
               vectmp (:) = 1.
               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask
               ENDIF
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! soil temperature [K]
         call flux_map_and_write_3d ( DEF_hist_vars%t_soisno, &
            a_t_soisno, f_t_soisno, file_hist, 'f_t_soisno', 'soilsnow', itime_in_file, sumarea, filter, &
            'soil temperature','K')

         ! liquid water in soil layers [kg/m2]
         call flux_map_and_write_3d ( DEF_hist_vars%wliq_soisno, &
            a_wliq_soisno, f_wliq_soisno, file_hist, 'f_wliq_soisno', 'soilsnow', &
            itime_in_file, sumarea, filter,'liquid water in soil layers','kg/m2')

         ! ice lens in soil layers [kg/m2]
         call flux_map_and_write_3d ( DEF_hist_vars%wice_soisno, &
            a_wice_soisno, f_wice_soisno, file_hist, 'f_wice_soisno', 'soilsnow', &
            itime_in_file, sumarea, filter,'ice lens in soil layers','kg/m2')

         ! --------------------------------------------------------------------
         ! additial diagnostic variables for output (vegetated land only <=2)
         ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3; land water bodies => 4; ocean => 99]
         ! --------------------------------------------------------------------

         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype <= 2
               vectmp(:) = 1.
               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask
               ENDIF
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! volumetric soil water in layers [m3/m3]
         call flux_map_and_write_3d ( DEF_hist_vars%h2osoi, &
            a_h2osoi, f_h2osoi, file_hist, 'f_h2osoi', 'soil', itime_in_file, sumarea, filter, &
            'volumetric water in soil layers','m3/m3')

         ! fraction of root water uptake from each soil layer, all layers add to 1, when PHS is not defined
         ! water exchange between soil layers and root. Positive: soil->root [mm h2o/s], when PHS is defined
         call flux_map_and_write_3d ( DEF_hist_vars%rootr, &
            a_rootr, f_rootr, file_hist, 'f_rootr', 'soil', itime_in_file, sumarea, filter, &
            'root water uptake', 'mm h2o/s')

         if(DEF_USE_PLANTHYDRAULICS)then
         ! vegetation water potential [mm]
            call flux_map_and_write_3d ( DEF_hist_vars%vegwp, &
               a_vegwp, f_vegwp, file_hist, 'f_vegwp', 'vegnodes', itime_in_file, sumarea, filter, &
               'vegetation water potential', 'mm')
         end if

         ! water table depth [m]
         call flux_map_and_write_2d ( DEF_hist_vars%zwt, &
            a_zwt, f_zwt, file_hist, 'f_zwt', itime_in_file, sumarea, filter, &
            'the depth to water table','m')

         ! water storage in aquifer [mm]
         call flux_map_and_write_2d ( DEF_hist_vars%wa, &
            a_wa, f_wa, file_hist, 'f_wa', itime_in_file, sumarea, filter, &
            'water storage in aquifer','mm')

         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = (patchtype <= 2) .or. (patchtype == 4)
               vectmp(:) = 1.
               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask
               ENDIF
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! depth of surface water [m]
         call flux_map_and_write_2d ( DEF_hist_vars%wdsrf, &
            a_wdsrf, f_wdsrf, file_hist, 'f_wdsrf', itime_in_file, sumarea, filter, &
            'depth of surface water','mm')

         ! -----------------------------------------------
         ! Land water bodies' ice fraction and temperature
         ! -----------------------------------------------

         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype == 4
               vectmp(:) = 1.
               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask
               ENDIF
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! lake temperature [K]
         call flux_map_and_write_3d ( DEF_hist_vars%t_lake, &
            a_t_lake, f_t_lake, file_hist, 'f_t_lake', 'lake', itime_in_file, sumarea, filter, &
            'lake temperature','K')

         ! lake ice fraction cover [0-1]
         call flux_map_and_write_3d ( DEF_hist_vars%lake_icefrac, &
            a_lake_icefrac, f_lake_icefrac, file_hist, 'f_lake_icefrac', &
            'lake', itime_in_file, sumarea, filter,'lake ice fraction cover','0-1')

         ! --------------------------------
         ! Retrieve through averaged fluxes
         ! --------------------------------
         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = patchtype < 99
               vectmp(:) = 1.
               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask
               ENDIF
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, spv = spval, msk = filter)

         ! u* in similarity theory [m/s]
         call flux_map_and_write_2d ( DEF_hist_vars%ustar, &
            a_ustar, f_ustar, file_hist, 'f_ustar', itime_in_file, sumarea, filter, &
            'u* in similarity theory','m/s')

         ! t* in similarity theory [kg/kg]
         call flux_map_and_write_2d ( DEF_hist_vars%tstar, &
            a_tstar, f_tstar, file_hist, 'f_tstar', itime_in_file, sumarea, filter, &
            't* in similarity theory','kg/kg')

         ! q* in similarity theory [kg/kg]
         call flux_map_and_write_2d ( DEF_hist_vars%qstar, &
            a_qstar, f_qstar, file_hist, 'f_qstar', itime_in_file, sumarea, filter, &
            'q* in similarity theory', 'kg/kg')

         ! dimensionless height (z/L) used in Monin-Obukhov theory
         call flux_map_and_write_2d ( DEF_hist_vars%zol, &
            a_zol, f_zol, file_hist, 'f_zol', itime_in_file, sumarea, filter, &
            'dimensionless height (z/L) used in Monin-Obukhov theory','-')

         ! bulk Richardson number in surface layer
         call flux_map_and_write_2d ( DEF_hist_vars%rib, &
            a_rib, f_rib, file_hist, 'f_rib', itime_in_file, sumarea, filter, &
            'bulk Richardson number in surface layer','-')

         ! integral of profile function for momentum
         call flux_map_and_write_2d ( DEF_hist_vars%fm, &
            a_fm, f_fm, file_hist, 'f_fm', itime_in_file, sumarea, filter, &
            'integral of profile function for momentum','-')

         ! integral of profile function for heat
         call flux_map_and_write_2d ( DEF_hist_vars%fh, &
            a_fh, f_fh, file_hist, 'f_fh', itime_in_file, sumarea, filter, &
            'integral of profile function for heat','-')

         ! integral of profile function for moisture
         call flux_map_and_write_2d ( DEF_hist_vars%fq, &
            a_fq, f_fq, file_hist, 'f_fq', itime_in_file, sumarea, filter, &
            'integral of profile function for moisture','-')

         ! 10m u-velocity [m/s]
         call flux_map_and_write_2d ( DEF_hist_vars%us10m, &
            a_us10m, f_us10m, file_hist, 'f_us10m', itime_in_file, sumarea, filter, &
            '10m u-velocity','m/s')

         ! 10m v-velocity [m/s]
         call flux_map_and_write_2d ( DEF_hist_vars%vs10m, &
            a_vs10m, f_vs10m, file_hist, 'f_vs10m', itime_in_file, sumarea, filter, &
            '10m v-velocity','m/s')

         ! integral of profile function for momentum at 10m [-]
         call flux_map_and_write_2d ( DEF_hist_vars%fm10m, &
            a_fm10m, f_fm10m, file_hist, 'f_fm10m', itime_in_file, sumarea, filter, &
            'integral of profile function for momentum at 10m','-')

         ! total reflected solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%sr, &
            a_sr, f_sr, file_hist, 'f_sr', itime_in_file, sumarea, filter, &
            'reflected solar radiation at surface [W/m2]','W/m2')

         ! incident direct beam vis solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%solvd, &
            a_solvd, f_solvd, file_hist, 'f_solvd', itime_in_file, sumarea, filter, &
            'incident direct beam vis solar radiation (W/m2)','W/m2')

         ! incident diffuse beam vis solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%solvi, &
            a_solvi, f_solvi, file_hist, 'f_solvi', itime_in_file, sumarea, filter, &
            'incident diffuse beam vis solar radiation (W/m2)','W/m2')

         ! incident direct beam nir solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%solnd, &
            a_solnd, f_solnd, file_hist, 'f_solnd', itime_in_file, sumarea, filter, &
            'incident direct beam nir solar radiation (W/m2)','W/m2')

         ! incident diffuse beam nir solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%solni, &
            a_solni, f_solni, file_hist, 'f_solni', itime_in_file, sumarea, filter, &
            'incident diffuse beam nir solar radiation (W/m2)','W/m2')

         ! reflected direct beam vis solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%srvd, &
            a_srvd, f_srvd, file_hist, 'f_srvd', itime_in_file, sumarea, filter, &
            'reflected direct beam vis solar radiation (W/m2)','W/m2')

         ! reflected diffuse beam vis solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%srvi, &
            a_srvi, f_srvi, file_hist, 'f_srvi', itime_in_file, sumarea, filter, &
            'reflected diffuse beam vis solar radiation (W/m2)','W/m2')

         ! reflected direct beam nir solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%srnd, &
            a_srnd, f_srnd, file_hist, 'f_srnd', itime_in_file, sumarea, filter, &
            'reflected direct beam nir solar radiation (W/m2)','W/m2')

         ! reflected diffuse beam nir solar radiation (W/m2)
         call flux_map_and_write_2d ( DEF_hist_vars%srni, &
            a_srni, f_srni, file_hist, 'f_srni', itime_in_file, sumarea, filter, &
            'reflected diffuse beam nir solar radiation (W/m2)','W/m2')

         ! local noon fluxes
         if (p_is_worker) then
            if (numpatch > 0) then
               filter(:) = nac_ln > 0
               vectmp(:) = 1.
               IF (DEF_forcing%has_missing_value) THEN
                  filter = filter .and. forcmask
               ENDIF
            end if
         end if

         call mp2g_hist%map (vectmp, sumarea, msk = filter)

         ! incident direct beam vis solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%solvdln, &
            a_solvdln, f_solvdln, file_hist, 'f_solvdln', itime_in_file, sumarea, filter, &
            'incident direct beam vis solar radiation at local noon(W/m2)','W/m2')

         ! incident diffuse beam vis solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%solviln, &
            a_solviln, f_solviln, file_hist, 'f_solviln', itime_in_file, sumarea, filter, &
            'incident diffuse beam vis solar radiation at local noon(W/m2)','W/m2')

         ! incident direct beam nir solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%solndln, &
            a_solndln, f_solndln, file_hist, 'f_solndln', itime_in_file, sumarea, filter, &
            'incident direct beam nir solar radiation at local noon(W/m2)','W/m2')

         ! incident diffuse beam nir solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%solniln, &
            a_solniln, f_solniln, file_hist, 'f_solniln', itime_in_file, sumarea, filter, &
            'incident diffuse beam nir solar radiation at local noon(W/m2)','W/m2')

         ! reflected direct beam vis solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%srvdln, &
            a_srvdln, f_srvdln, file_hist, 'f_srvdln', itime_in_file, sumarea, filter, &
            'reflected direct beam vis solar radiation at local noon(W/m2)','W/m2')

         ! reflected diffuse beam vis solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%srviln, &
            a_srviln, f_srviln, file_hist, 'f_srviln', itime_in_file, sumarea, filter, &
            'reflected diffuse beam vis solar radiation at local noon(W/m2)','W/m2')

         ! reflected direct beam nir solar radiation at local noon (W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%srndln, &
            a_srndln, f_srndln, file_hist, 'f_srndln', itime_in_file, sumarea, filter, &
            'reflected direct beam nir solar radiation at local noon(W/m2)','W/m2')

         ! reflected diffuse beam nir solar radiation at local noon(W/m2)
         call flux_map_and_write_ln ( DEF_hist_vars%srniln, &
            a_srniln, f_srniln, file_hist, 'f_srniln', itime_in_file, sumarea, filter, &
            'reflected diffuse beam nir solar radiation at local noon(W/m2)','W/m2')

#if(defined CaMa_Flood)
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
if (p_is_master) then
         CALL hist_out_cama (file_hist_cama, itime_in_file_cama)
ENDIF
#endif

         if (allocated(filter)) deallocate(filter)
         if (allocated(vectmp)) deallocate(vectmp)
#ifdef URBAN_MODEL
         if (allocated(filter_urb)) deallocate(filter_urb)
         if (allocated(vectmp_urb)) deallocate(vectmp_urb)
#endif
         call FLUSH_acc_fluxes ()

      end if

   END SUBROUTINE hist_out

   ! -------
   subroutine flux_map_and_write_2d ( is_hist, &
         acc_vec, flux_xy, file_hist, varname, itime_in_file, sumarea, filter, &
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

      logical, intent(in) :: is_hist

      real(r8), intent(inout) :: acc_vec(:)
      type(block_data_real8_2d), intent(inout) :: flux_xy
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer,          intent(in) :: itime_in_file
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units

      type(block_data_real8_2d), intent(in) :: sumarea
      logical, intent(in) :: filter(:)

      ! Local variables
      integer :: iblkme, xblk, yblk, xloc, yloc
      integer :: compress

      if (.not. is_hist) return

      if (p_is_worker) then
         where (acc_vec /= spval)  acc_vec = acc_vec / nac
      end if

      call mp2g_hist%map (acc_vec, flux_xy, spv = spval, msk = filter)

      if (p_is_io) then
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            do yloc = 1, ghist%ycnt(yblk)
               do xloc = 1, ghist%xcnt(xblk)

                  if (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                     IF (flux_xy%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                        flux_xy%blk(xblk,yblk)%val(xloc,yloc) &
                           = flux_xy%blk(xblk,yblk)%val(xloc,yloc) &
                           / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                     ENDIF
                  else
                     flux_xy%blk(xblk,yblk)%val(xloc,yloc) = spval
                  end if

               end do
            end do

         end do
      end if

      compress = DEF_HIST_COMPRESS_LEVEL
      call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy, compress)

      IF (p_is_master .and. (itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_2d

   ! -------
   subroutine flux_map_and_write_urb_2d ( is_hist, &
         acc_vec, flux_xy, file_hist, varname, itime_in_file, sumarea, filter, &
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

      logical, intent(in) :: is_hist

      real(r8), intent(inout) :: acc_vec(:)
      type(block_data_real8_2d), intent(inout) :: flux_xy
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer,          intent(in) :: itime_in_file
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units

      type(block_data_real8_2d), intent(in) :: sumarea
      logical, intent(in) :: filter(:)

      ! Local variables
      integer :: iblkme, xblk, yblk, xloc, yloc
      integer :: compress

      if (.not. is_hist) return

      if (p_is_worker) then
         where (acc_vec /= spval)  acc_vec = acc_vec / nac
      end if

      call mp2g_hist_urb%map (acc_vec, flux_xy, spv = spval, msk = filter)

      if (p_is_io) then
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            do yloc = 1, ghist%ycnt(yblk)
               do xloc = 1, ghist%xcnt(xblk)

                  if (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                     IF (flux_xy%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                        flux_xy%blk(xblk,yblk)%val(xloc,yloc) &
                           = flux_xy%blk(xblk,yblk)%val(xloc,yloc) &
                           / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                     ENDIF
                  else
                     flux_xy%blk(xblk,yblk)%val(xloc,yloc) = spval
                  end if

               end do
            end do

         end do
      end if

      compress = DEF_HIST_COMPRESS_LEVEL
      call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy, compress)

      IF (p_is_master .and. (itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_urb_2d

   ! -------
   subroutine flux_map_and_write_3d ( is_hist, &
         acc_vec, flux_xy, file_hist, varname, dim1name, &
         itime_in_file, sumarea, filter, longname, units)

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

      logical, intent(in) :: is_hist

      real(r8), intent(inout) :: acc_vec(:,:)
      type(block_data_real8_3d), intent(inout) :: flux_xy
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: dim1name
      integer, intent(in) :: itime_in_file

      type(block_data_real8_2d), intent(in) :: sumarea
      logical, intent(in) :: filter(:)
      character (len=*), intent(in) :: longname
      character (len=*), intent(in) :: units

      ! Local variables
      integer :: iblkme, xblk, yblk, xloc, yloc, i1
      integer :: compress

      if (.not. is_hist) return

      if (p_is_worker) then
         where (acc_vec /= spval)  acc_vec = acc_vec / nac
      end if

      call mp2g_hist%map (acc_vec, flux_xy, spv = spval, msk = filter)

      if (p_is_io) then
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            do yloc = 1, ghist%ycnt(yblk)
               do xloc = 1, ghist%xcnt(xblk)

                  if (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                     DO i1 = flux_xy%lb1, flux_xy%ub1
                        IF (flux_xy%blk(xblk,yblk)%val(i1,xloc,yloc) /= spval) THEN
                           flux_xy%blk(xblk,yblk)%val(i1,xloc,yloc) &
                              = flux_xy%blk(xblk,yblk)%val(i1,xloc,yloc) &
                              / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                        ENDIF
                     ENDDO
                  else
                     flux_xy%blk(xblk,yblk)%val(:,xloc,yloc) = spval
                  end if

               end do
            end do

         end do
      end if

      compress = DEF_HIST_COMPRESS_LEVEL
      call hist_write_var_real8_3d (file_hist, varname, dim1name, ghist, &
         itime_in_file, flux_xy, compress)

      IF (p_is_master .and. (itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_3d

   ! -------
   subroutine flux_map_and_write_4d ( is_hist, &
         acc_vec, flux_xy, file_hist, varname, dim1name, dim2name, itime_in_file, sumarea, filter, &
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

      logical, intent(in) :: is_hist

      real(r8), intent(inout) :: acc_vec(:,:,:)
      type(block_data_real8_4d), intent(inout) :: flux_xy
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: dim1name, dim2name
      integer, intent(in) :: itime_in_file

      type(block_data_real8_2d), intent(in) :: sumarea
      logical, intent(in) :: filter(:)
      character (len=*), intent(in) :: longname
      character (len=*), intent(in) :: units

      ! Local variables
      integer :: iblkme, xblk, yblk, xloc, yloc, i1, i2
      integer :: compress

      if (.not. is_hist) return

      if (p_is_worker) then
         where(acc_vec /= spval)  acc_vec = acc_vec / nac
      end if

      call mp2g_hist%map (acc_vec, flux_xy, spv = spval, msk = filter)

      if (p_is_io) then
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            do yloc = 1, ghist%ycnt(yblk)
               do xloc = 1, ghist%xcnt(xblk)

                  if (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                     DO i1 = flux_xy%lb1, flux_xy%ub1
                        DO i2 = flux_xy%lb2, flux_xy%ub2
                           IF (flux_xy%blk(xblk,yblk)%val(i1,i2,xloc,yloc) /= spval) THEN
                              flux_xy%blk(xblk,yblk)%val(i1,i2,xloc,yloc) &
                                 = flux_xy%blk(xblk,yblk)%val(i1,i2,xloc,yloc) &
                                 / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                           ENDIF
                        ENDDO
                     ENDDO
                  else
                     flux_xy%blk(xblk,yblk)%val(:,:,xloc,yloc) = spval
                  end if

               end do
            end do

         end do
      end if

      compress = DEF_HIST_COMPRESS_LEVEL
      call hist_write_var_real8_4d (file_hist, varname, dim1name, dim2name, &
         ghist, itime_in_file, flux_xy, compress)

      IF (p_is_master .and. (itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_4d

   ! -------
   subroutine flux_map_and_write_ln ( is_hist, &
         acc_vec, flux_xy, file_hist, varname, itime_in_file, sumarea, filter, &
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

      logical, intent(in) :: is_hist

      real(r8), intent(inout) :: acc_vec(:)
      type(block_data_real8_2d), intent(inout) :: flux_xy
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer, intent(in) :: itime_in_file

      type(block_data_real8_2d), intent(in) :: sumarea
      logical,  intent(in) :: filter(:)
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units

      ! Local variables
      integer :: i, iblkme, xblk, yblk, xloc, yloc
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
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            do yloc = 1, ghist%ycnt(yblk)
               do xloc = 1, ghist%xcnt(xblk)

                  if ((sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) &
                     .and. (flux_xy%blk(xblk,yblk)%val(xloc,yloc) /= spval)) then
                     flux_xy%blk(xblk,yblk)%val(xloc,yloc) &
                        = flux_xy%blk(xblk,yblk)%val(xloc,yloc) &
                        / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                  else
                     flux_xy%blk(xblk,yblk)%val(xloc,yloc) = spval
                  end if

               end do
            end do

         end do
      end if

      compress = DEF_HIST_COMPRESS_LEVEL
      call hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, flux_xy, &
         compress)

      IF (p_is_master .and. (itime_in_file == 1) .and. (trim(DEF_HIST_mode) == 'one')) then
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
      ENDIF

   end subroutine flux_map_and_write_ln

   !------------------------------
   subroutine hist_write_time ( &
         filename, dataname, grid, time, itime)

      use MOD_Namelist
      use MOD_Grid
      use MOD_Block
      use MOD_SPMD_Task
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      type (grid_type),  intent(in) :: grid

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
            call ncio_write_time (filename, dataname, time, itime)

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
                  call hist_write_grid_info  (fileblock, grid, iblk, jblk)
               end if

               call ncio_write_time (fileblock, dataname, time, itime)

            end do

         end if
      endif

   end subroutine hist_write_time

   !----------------------------------------------------------------------------
   subroutine hist_write_var_real8_2d ( &
         filename, dataname, grid, itime, wdata, compress)

      use MOD_Namelist
      use MOD_Block
      use MOD_Grid
      use MOD_DataType
      use MOD_SPMD_Task
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

end module MOD_Hist
