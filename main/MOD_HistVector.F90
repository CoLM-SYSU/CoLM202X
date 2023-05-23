#include <define.h>

#if (defined UNSTRUCTURED || defined CATCHMENT)
module MOD_HistVector

   use precision
   USE spmd_task
   USE mod_namelist
   USE MOD_Vars_Global, only : spval
   USE mod_mesh
   USE mod_landelm
#ifdef CATCHMENT
   USE mod_landhru
#endif
   USE mod_pixelset
   USE ncio_serial
#ifdef CATCHMENT
   USE mod_hru_vector
#else
   USE mod_elm_vector
#endif

   public :: hist_vector_out

contains
! ----- subroutines ------

   SUBROUTINE hist_vector_out (file_hist, idate)

      use precision
      use mod_namelist
      use timemanager
      use spmd_task
      use MOD_Vars_1DAccFluxes
      use mod_landpatch
      use mod_colm_debug
      use MOD_Vars_Global, only: spval
      USE MOD_Vars_TimeInvariants, only : patchtype
      USE MOD_Forcing, only : forcmask
      IMPLICIT NONE

      character(LEN=*), intent(in) :: file_hist
      integer,  INTENT(in) :: idate(3)

      ! Local variables
      LOGICAL :: fexists
      integer :: itime_in_file
      logical,  allocatable ::  filter(:)

      if (p_is_master) then
         inquire (file=file_hist, exist=fexists)
         if (.not. fexists) then
            call ncio_create_file (trim(file_hist))
            CALL ncio_define_dimension(file_hist, 'time', 0)

#ifdef CATCHMENT
            call ncio_define_dimension(file_hist, 'hydrounit', totalnumhru)

            call ncio_write_serial (file_hist, 'bsn_hru', eindx_hru, 'hydrounit')
            CALL ncio_put_attr (file_hist, 'bsn_hru', 'long_name', &
               'basin index of hydrological units in mesh')

            call ncio_write_serial (file_hist, 'typ_hru' , htype_hru, 'hydrounit')
            CALL ncio_put_attr (file_hist, 'typ_hru' , 'long_name', &
               'index of hydrological units inside basin')
#else
            call ncio_define_dimension(file_hist, 'element', totalnumelm)
            call ncio_write_serial (file_hist, 'elmindex', eindex_glb, 'element')
            CALL ncio_put_attr (file_hist, 'elmindex', 'long_name', &
               'element index in mesh')

#endif
         endif

         call ncio_write_time (file_hist, 'time', idate, itime_in_file)

      ENDIF

      ! ---------------------------------------------------
      ! Meteorological forcing
      ! ---------------------------------------------------
      if (p_is_worker) then
         if (numpatch > 0) then
            allocate (filter (numpatch))
            filter(:) = patchtype < 99
            IF (DEF_forcing%has_missing_value) THEN
               filter = filter .and. forcmask
            ENDIF
         end if
      end if

      ! wind in eastward direction [m/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xy_us, &
         a_us,  file_hist, 'f_xy_us', itime_in_file, filter, &
         'wind in eastward direction', 'm/s')

      ! wind in northward direction [m/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xy_vs, &
         a_vs, file_hist, 'f_xy_vs', itime_in_file, filter, &
         'wind in northward direction','m/s')

      ! temperature at reference height [kelvin]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xy_t, &
         a_t, file_hist, 'f_xy_t', itime_in_file, filter, &
         'temperature at reference height','kelvin')

      ! specific humidity at reference height [kg/kg]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xy_q, &
         a_q, file_hist, 'f_xy_q', itime_in_file, filter, &
         'specific humidity at reference height','kg/kg')

      ! convective precipitation [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xy_prc, &
         a_prc, file_hist, 'f_xy_prc', itime_in_file, filter, &
         'convective precipitation','mm/s')

      ! large scale precipitation [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xy_prl, &
         a_prl, file_hist, 'f_xy_prl', itime_in_file, filter, &
         'large scale precipitation','mm/s')

      ! atmospheric pressure at the surface [pa]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xy_pbot, &
         a_pbot, file_hist, 'f_xy_pbot', itime_in_file, filter, &
         'atmospheric pressure at the surface','pa')

      ! atmospheric infrared (longwave) radiation [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xy_frl, &
         a_frl, file_hist, 'f_xy_frl', itime_in_file, filter, &
         'atmospheric infrared (longwave) radiation','W/m2')

      ! downward solar radiation at surface [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xy_solarin, &
         a_solarin, file_hist, 'f_xy_solarin', itime_in_file, filter, &
         'downward solar radiation at surface','W/m2')

      ! rain [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xy_rain, &
         a_rain, file_hist, 'f_xy_rain', itime_in_file, filter, &
         'rain','mm/s')

      ! snow [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xy_snow, &
         a_snow, file_hist, 'f_xy_snow', itime_in_file, filter, &
         'snow','mm/s')

      ! ------------------------------------------------------------------------------------------
      ! Mapping the fluxes and state variables at patch [numpatch] to grid
      ! ------------------------------------------------------------------------------------------

      ! wind stress: E-W [kg/m/s2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%taux, &
         a_taux, file_hist, 'f_taux', itime_in_file, filter, &
         'wind stress: E-W','kg/m/s2')

      ! wind stress: N-S [kg/m/s2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%tauy, &
         a_tauy, file_hist, 'f_tauy', itime_in_file, filter, &
         'wind stress: N-S','kg/m/s2')

      ! sensible heat from canopy height to atmosphere [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fsena, &
         a_fsena, file_hist, 'f_fsena', itime_in_file, filter, &
         'sensible heat from canopy height to atmosphere','W/m2')

      ! latent heat flux from canopy height to atmosphere [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%lfevpa, &
         a_lfevpa, file_hist, 'f_lfevpa', itime_in_file, filter, &
         'latent heat flux from canopy height to atmosphere','W/m2')

      ! evapotranspiration from canopy to atmosphere [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fevpa, &
         a_fevpa, file_hist, 'f_fevpa', itime_in_file, filter, &
         'evapotranspiration from canopy height to atmosphere','mm/s')

      ! sensible heat from leaves [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fsenl, &
         a_fsenl, file_hist, 'f_fsenl', itime_in_file, filter, &
         'sensible heat from leaves','W/m2')

      ! evaporation+transpiration from leaves [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fevpl, &
         a_fevpl, file_hist, 'f_fevpl', itime_in_file, filter, &
         'evaporation+transpiration from leaves','mm/s')

      ! transpiration rate [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%etr, &
         a_etr, file_hist, 'f_etr', itime_in_file, filter, &
         'transpiration rate','mm/s')

      ! sensible heat flux from ground [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fseng, &
         a_fseng, file_hist, 'f_fseng', itime_in_file, filter, &
         'sensible heat flux from ground','W/m2')

      ! evaporation heat flux from ground [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fevpg, &
         a_fevpg, file_hist, 'f_fevpg', itime_in_file, filter, &
         'evaporation heat flux from ground','mm/s')

      ! ground heat flux [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fgrnd, &
         a_fgrnd, file_hist, 'f_fgrnd', itime_in_file, filter, &
         'ground heat flux','W/m2')

      ! solar absorbed by sunlit canopy [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%sabvsun, &
         a_sabvsun, file_hist, 'f_sabvsun', itime_in_file, filter, &
         'solar absorbed by sunlit canopy','W/m2')

      ! solar absorbed by shaded [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%sabvsha, &
         a_sabvsha, file_hist, 'f_sabvsha', itime_in_file, filter, &
         'solar absorbed by shaded','W/m2')

      ! solar absorbed by ground  [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%sabg, &
         a_sabg, file_hist, 'f_sabg', itime_in_file, filter, &
         'solar absorbed by ground','W/m2')

      ! outgoing long-wave radiation from ground+canopy [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%olrg, &
         a_olrg, file_hist, 'f_olrg', itime_in_file, filter, &
         'outgoing long-wave radiation from ground+canopy','W/m2')

      ! net radiation [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%rnet, &
         a_rnet, file_hist, 'f_rnet', itime_in_file, filter, &
         'net radiation','W/m2')

      ! the error of water banace [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%xerr, &
         a_xerr, file_hist, 'f_xerr', itime_in_file, filter, &
         'the error of water banace','mm/s')

      ! the error of energy balance [W/m2]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%zerr, &
         a_zerr, file_hist, 'f_zerr', itime_in_file, filter, &
         'the error of energy balance','W/m2')

      ! surface runoff [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%rsur, &
         a_rsur, file_hist, 'f_rsur', itime_in_file, filter, &
         'surface runoff','mm/s')

      ! subsurface runoff [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%rsub, &
         a_rsub, file_hist, 'f_rsub', itime_in_file, filter, &
         'subsurface runoff','mm/s')

      ! total runoff [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%rnof, &
         a_rnof, file_hist, 'f_rnof', itime_in_file, filter, &
         'total runoff','mm/s')

      ! interception [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%qintr, &
         a_qintr, file_hist, 'f_qintr', itime_in_file, filter, &
         'interception','mm/s')

      ! inflitraton [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%qinfl, &
         a_qinfl, file_hist, 'f_qinfl', itime_in_file, filter, &
         'f_qinfl','mm/s')

      ! throughfall [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%qdrip, &
         a_qdrip, file_hist, 'f_qdrip', itime_in_file, filter, &
         'total throughfall','mm/s')

      ! total water storage [mm]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%wat, &
         a_wat, file_hist, 'f_wat', itime_in_file, filter, &
         'total water storage','mm')

      ! canopy assimilation rate [mol m-2 s-1]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%assim, &
         a_assim, file_hist, 'f_assim', itime_in_file, filter, &
         'canopy assimilation rate','umol m-2 s-1')

      ! respiration (plant+soil) [mol m-2 s-1]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%respc, &
         a_respc, file_hist, 'f_respc', itime_in_file, filter, &
         'respiration (plant+soil)','mol m-2 s-1')

      ! groundwater recharge rate [mm/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%qcharge, &
         a_qcharge, file_hist, 'f_qcharge', itime_in_file, filter, &
         'groundwater recharge rate','mm/s')

      ! ground surface temperature [K]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%t_grnd, &
         a_t_grnd, file_hist, 'f_t_grnd', itime_in_file, filter, &
         'ground surface temperature','K')

      ! leaf temperature [K]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%tleaf, &
         a_tleaf, file_hist, 'f_tleaf', itime_in_file, filter, &
         'leaf temperature','K')

      ! depth of water on foliage [mm]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%ldew, &
         a_ldew, file_hist, 'f_ldew', itime_in_file, filter, &
         'depth of water on foliage','mm')

      ! snow cover, water equivalent [mm]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%scv, &
         a_scv, file_hist, 'f_scv', itime_in_file, filter, &
         'snow cover, water equivalent','mm')

      ! snow depth [meter]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%snowdp, &
         a_snowdp, file_hist, 'f_snowdp', itime_in_file, filter, &
         'snow depth','meter')

      ! fraction of snow cover on ground
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fsno, &
         a_fsno, file_hist, 'f_fsno', itime_in_file, filter, &
         'fraction of snow cover on ground','-')

      ! fraction of veg cover, excluding snow-covered veg [-]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%sigf, &
         a_sigf, file_hist, 'f_sigf', itime_in_file, filter, &
         'fraction of veg cover, excluding snow-covered veg','-')

      ! leaf greenness
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%green, &
         a_green, file_hist, 'f_green', itime_in_file, filter, &
         'leaf greenness','-')

      ! leaf area index
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%lai, &
         a_lai, file_hist, 'f_lai', itime_in_file, filter, &
         'leaf area index','m2/m2')

      ! leaf area index
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%laisun, &
         a_laisun, file_hist, 'f_laisun', itime_in_file, filter, &
         'sunlit leaf area index','m2/m2')

      ! leaf area index
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%laisha, &
         a_laisha, file_hist, 'f_laisha', itime_in_file, filter, &
         'shaded leaf area index','m2/m2')

      ! stem area index
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%sai, &
         a_sai, file_hist, 'f_sai', itime_in_file, filter, &
         'stem area index','m2/m2')

      ! averaged albedo [visible, direct; direct, diffuse]
      call aggregate_to_vector_and_write_4d ( DEF_hist_vars%alb, &
         a_alb, file_hist, 'f_alb', itime_in_file, filter, &
         'band', 2, 'rtyp', 2, &
         'averaged albedo direct','%')

      ! averaged bulk surface emissivity
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%emis, &
         a_emis, file_hist, 'f_emis', itime_in_file, filter, &
         'averaged bulk surface emissivity','-')

      ! effective roughness [m]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%z0m, &
         a_z0m, file_hist, 'f_z0m', itime_in_file, filter, &
         'effective roughness','m')

      ! radiative temperature of surface [K]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%trad, &
         a_trad, file_hist, 'f_trad', itime_in_file, filter, &
         'radiative temperature of surface','kelvin')

      ! 2 m height air temperature [kelvin]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%tref, &
         a_tref, file_hist, 'f_tref', itime_in_file, filter, &
         '2 m height air temperature','kelvin')

      ! 2 m height air specific humidity [kg/kg]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%qref, &
         a_qref, file_hist, 'f_qref', itime_in_file, filter, &
         '2 m height air specific humidity','kg/kg')

#ifdef BGC
      ! leaf carbon display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%leafc, &
         a_leafc, file_hist, 'f_leafc', itime_in_file, filter, &
         'leaf carbon display pool','gC/m2')

      ! leaf carbon storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%leafc_storage, &
         a_leafc_storage, file_hist, 'f_leafc_storage', itime_in_file, filter, &
         'leaf carbon storage pool','gC/m2')

      ! leaf carbon transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%leafc_xfer, &
         a_leafc_xfer, file_hist, 'f_leafc_xfer', itime_in_file, filter, &
         'leaf carbon transfer pool','gC/m2')

      ! fine root carbon display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%frootc, &
         a_frootc, file_hist, 'f_frootc', itime_in_file, filter, &
         'fine root carbon display pool','gC/m2')

      ! fine root carbon storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%frootc_storage, &
         a_frootc_storage, file_hist, 'f_frootc_storage', itime_in_file, filter, &
         'fine root carbon storage pool','gC/m2')

      ! fine root carbon transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%frootc_xfer, &
         a_frootc_xfer, file_hist, 'f_frootc_xfer', itime_in_file, filter, &
         'fine root carbon transfer pool','gC/m2')

      ! live stem carbon display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livestemc, &
         a_livestemc, file_hist, 'f_livestemc', itime_in_file, filter, &
         'live stem carbon display pool','gC/m2')

      ! live stem carbon storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livestemc_storage, &
         a_livestemc_storage, file_hist, 'f_livestemc_storage', itime_in_file, filter, &
         'live stem carbon storage pool','gC/m2')

      ! live stem carbon transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livestemc_xfer, &
         a_livestemc_xfer, file_hist, 'f_livestemc_xfer', itime_in_file, filter, &
         'live stem carbon transfer pool','gC/m2')

      ! dead stem carbon display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadstemc, &
         a_deadstemc, file_hist, 'f_deadstemc', itime_in_file, filter, &
         'dead stem carbon display pool','gC/m2')

      ! dead stem carbon storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadstemc_storage, &
         a_deadstemc_storage, file_hist, 'f_deadstemc_storage', itime_in_file, filter, &
         'dead stem carbon storage pool','gC/m2')

      ! dead stem carbon transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadstemc_xfer, &
         a_deadstemc_xfer, file_hist, 'f_deadstemc_xfer', itime_in_file, filter, &
         'dead stem carbon transfer pool','gC/m2')

      ! live coarse root carbon display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livecrootc, &
         a_livecrootc, file_hist, 'f_livecrootc', itime_in_file, filter, &
         'live coarse root carbon display pool','gC/m2')

      ! live coarse root carbon storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livecrootc_storage, &
         a_livecrootc_storage, file_hist, 'f_livecrootc_storage', itime_in_file, filter, &
         'live coarse root carbon storage pool','gC/m2')

      ! live coarse root carbon transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livecrootc_xfer, &
         a_livecrootc_xfer, file_hist, 'f_livecrootc_xfer', itime_in_file, filter, &
         'live coarse root carbon transfer pool','gC/m2')

      ! dead coarse root carbon display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadcrootc, &
         a_deadcrootc, file_hist, 'f_deadcrootc', itime_in_file, filter, &
         'dead coarse root carbon display pool','gC/m2')

      ! dead coarse root carbon storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadcrootc_storage, &
         a_deadcrootc_storage, file_hist, 'f_deadcrootc_storage', itime_in_file, filter, &
         'dead coarse root carbon storage pool','gC/m2')

      ! dead coarse root carbon transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadcrootc_xfer, &
         a_deadcrootc_xfer, file_hist, 'f_deadcrootc_xfer', itime_in_file, filter, &
         'dead coarse root carbon transfer pool','gC/m2')

#ifdef CROP
      ! grain carbon display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%grainc, &
         a_grainc, file_hist, 'f_grainc', itime_in_file, filter, &
         'grain carbon display pool','gC/m2')

      ! grain carbon storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%grainc_storage, &
         a_grainc_storage, file_hist, 'f_grainc_storage', itime_in_file, filter, &
         'grain carbon storage pool','gC/m2')

      ! grain carbon transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%grainc_xfer, &
         a_grainc_xfer, file_hist, 'f_grainc_xfer', itime_in_file, filter, &
         'grain carbon transfer pool','gC/m2')
#endif

      ! leaf nitrogen display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%leafn, &
         a_leafn, file_hist, 'f_leafn', itime_in_file, filter, &
         'leaf nitrogen display pool','gN/m2')

      ! leaf nitrogen storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%leafn_storage, &
         a_leafn_storage, file_hist, 'f_leafn_storage', itime_in_file, filter, &
         'leaf nitrogen storage pool','gN/m2')

      ! leaf nitrogen transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%leafn_xfer, &
         a_leafn_xfer, file_hist, 'f_leafn_xfer', itime_in_file, filter, &
         'leaf nitrogen transfer pool','gN/m2')

      ! fine root nitrogen display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%frootn, &
         a_frootn, file_hist, 'f_frootn', itime_in_file, filter, &
         'fine root nitrogen display pool','gN/m2')

      ! fine root nitrogen storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%frootn_storage, &
         a_frootn_storage, file_hist, 'f_frootn_storage', itime_in_file, filter, &
         'fine root nitrogen storage pool','gN/m2')

      ! fine root nitrogen transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%frootn_xfer, &
         a_frootn_xfer, file_hist, 'f_frootn_xfer', itime_in_file, filter, &
         'fine root nitrogen transfer pool','gN/m2')

      ! live stem nitrogen display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livestemn, &
         a_livestemn, file_hist, 'f_livestemn', itime_in_file, filter, &
         'live stem nitrogen display pool','gN/m2')

      ! live stem nitrogen storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livestemn_storage, &
         a_livestemn_storage, file_hist, 'f_livestemn_storage', itime_in_file, filter, &
         'live stem nitrogen storage pool','gN/m2')

      ! live stem nitrogen transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livestemn_xfer, &
         a_livestemn_xfer, file_hist, 'f_livestemn_xfer', itime_in_file, filter, &
         'live stem nitrogen transfer pool','gN/m2')

      ! dead stem nitrogen display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadstemn, &
         a_deadstemn, file_hist, 'f_deadstemn', itime_in_file, filter, &
         'dead stem nitrogen display pool','gN/m2')

      ! dead stem nitrogen storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadstemn_storage, &
         a_deadstemn_storage, file_hist, 'f_deadstemn_storage', itime_in_file, filter, &
         'dead stem nitrogen storage pool','gN/m2')

      ! dead stem nitrogen transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadstemn_xfer, &
         a_deadstemn_xfer, file_hist, 'f_deadstemn_xfer', itime_in_file, filter, &
         'dead stem nitrogen transfer pool','gN/m2')

      ! live coarse root nitrogen display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livecrootn, &
         a_livecrootn, file_hist, 'f_livecrootn', itime_in_file, filter, &
         'live coarse root nitrogen display pool','gN/m2')

      ! live coarse root nitrogen storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livecrootn_storage, &
         a_livecrootn_storage, file_hist, 'f_livecrootn_storage', itime_in_file, filter, &
         'live coarse root nitrogen storage pool','gN/m2')

      ! live coarse root nitrogen transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%livecrootn_xfer, &
         a_livecrootn_xfer, file_hist, 'f_livecrootn_xfer', itime_in_file, filter, &
         'live coarse root nitrogen transfer pool','gN/m2')

      ! dead coarse root nitrogen display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadcrootn, &
         a_deadcrootn, file_hist, 'f_deadcrootn', itime_in_file, filter, &
         'dead coarse root nitrogen display pool','gN/m2')

      ! dead coarse root nitrogen storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadcrootn_storage, &
         a_deadcrootn_storage, file_hist, 'f_deadcrootn_storage', itime_in_file, filter, &
         'dead coarse root nitrogen storage pool','gN/m2')

      ! dead coarse root nitrogen transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%deadcrootn_xfer, &
         a_deadcrootn_xfer, file_hist, 'f_deadcrootn_xfer', itime_in_file, filter, &
         'dead coarse root nitrogen transfer pool','gN/m2')

#ifdef CROP
      ! grain nitrogen display pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%grainn, &
         a_grainn, file_hist, 'f_grainn', itime_in_file, filter, &
         'grain nitrogen display pool','gN/m2')

      ! grain nitrogen storage pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%grainn_storage, &
         a_grainn_storage, file_hist, 'f_grainn_storage', itime_in_file, filter, &
         'grain nitrogen storage pool','gN/m2')

      ! grain nitrogen transfer pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%grainn_xfer, &
         a_grainn_xfer, file_hist, 'f_grainn_xfer', itime_in_file, filter, &
         'grain nitrogen transfer pool','gN/m2')
#endif

      ! retranslocation nitrogen pool
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%retrasn, &
         a_retransn, file_hist, 'f_retrasn', itime_in_file, filter, &
         'retranslocation nitrogen pool','gN/m2')

      ! gross primary productivity
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%gpp, &
         a_gpp, file_hist, 'f_gpp', itime_in_file, filter, &
         'gross primary productivity','gC/m2/s')

      ! gross primary productivity
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%downreg, &
         a_downreg, file_hist, 'f_downreg', itime_in_file, filter, &
         'gpp downregulation due to N limitation','unitless')

      ! autotrophic respiration
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%ar , &
         a_ar , file_hist, 'f_ar', itime_in_file, filter, &
         'autotrophic respiration','gC/m2/s')

#ifdef CROP
      ! crop phase
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%cphase, &
         a_cphase, file_hist, 'f_cphase', itime_in_file, filter, &
         'crop phase','unitless')

      ! 1-yr crop production carbon
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%cropprod1c, &
         a_cropprod1c, file_hist, 'f_cropprod1c', itime_in_file, filter, &
         '1-yr crop production carbon','gC/m2')

      ! loss rate of 1-yr crop production carbon
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%cropprod1c_loss, &
         a_cropprod1c_loss, file_hist, 'f_cropprod1c_loss', itime_in_file, filter, &
         'loss rate of 1-yr crop production carbon','gC/m2/s')

      ! crop seed deficit
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%cropseedc_deficit, &
         a_cropseedc_deficit, file_hist, 'f_cropseedc_deficit', itime_in_file, filter, &
         'crop seed deficit','gC/m2/s')

      ! grain to crop production carbon
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%grainc_to_cropprodc, &
         a_grainc_to_cropprodc, file_hist, 'f_grainc_to_cropprodc', itime_in_file, filter, &
         'grain to crop production carbon','gC/m2/s')

      ! grain to crop seed carbon
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%grainc_to_seed, &
         a_grainc_to_seed, file_hist, 'f_grainc_to_seed', itime_in_file, filter, &
         'grain to crop seed carbon','gC/m2/s')
#endif

      ! litter 1 carbon density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%litr1c_vr, &
         a_litr1c_vr, file_hist, 'f_litr1c_vr', itime_in_file, filter, &
         'soil', nl_soil, 'litter 1 carbon density in soil layers','gC/m3')

      ! litter 2 carbon density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%litr2c_vr, &
         a_litr2c_vr, file_hist, 'f_litr2c_vr', itime_in_file, filter, &
         'soil', nl_soil, 'litter 2 carbon density in soil layers','gC/m3')

      ! litter 3 carbon density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%litr3c_vr, &
         a_litr3c_vr, file_hist, 'f_litr3c_vr', itime_in_file, filter, &
         'soil', nl_soil, 'litter 3 carbon density in soil layers','gC/m3')

      ! soil 1 carbon density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%soil1c_vr, &
         a_soil1c_vr, file_hist, 'f_soil1c_vr', itime_in_file, filter, &
         'soil', nl_soil, 'soil 1 carbon density in soil layers','gC/m3')

      ! soil 2 carbon density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%soil2c_vr, &
         a_soil2c_vr, file_hist, 'f_soil2c_vr', itime_in_file, filter, &
         'soil', nl_soil, 'soil 2 carbon density in soil layers','gC/m3')

      ! soil 3 carbon density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%soil3c_vr, &
         a_soil3c_vr, file_hist, 'f_soil3c_vr', itime_in_file, filter, &
         'soil', nl_soil, 'soil 3 carbon density in soil layers','gC/m3')

      ! coarse woody debris carbon density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%cwdc_vr, &
         a_cwdc_vr, file_hist, 'f_cwdc_vr', itime_in_file, filter, &
         'soil', nl_soil, 'coarse woody debris carbon density in soil layers','gC/m3')

      ! litter 1 nitrogen density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%litr1n_vr, &
         a_litr1n_vr, file_hist, 'f_litr1n_vr', itime_in_file, filter, &
         'soil', nl_soil, 'litter 1 nitrogen density in soil layers','gN/m3')

      ! litter 2 nitrogen density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%litr2n_vr, &
         a_litr2n_vr, file_hist, 'f_litr2n_vr', itime_in_file, filter, &
         'soil', nl_soil, 'litter 2 nitrogen density in soil layers','gN/m3')

      ! litter 3 nitrogen density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%litr3n_vr, &
         a_litr3n_vr, file_hist, 'f_litr3n_vr', itime_in_file, filter, &
         'soil', nl_soil, 'litter 3 nitrogen density in soil layers','gN/m3')

      ! soil 1 nitrogen density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%soil1n_vr, &
         a_soil1n_vr, file_hist, 'f_soil1n_vr', itime_in_file, filter, &
         'soil', nl_soil, 'soil 1 nitrogen density in soil layers','gN/m3')

      ! soil 2 nitrogen density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%soil2n_vr, &
         a_soil2n_vr, file_hist, 'f_soil2n_vr', itime_in_file, filter, &
         'soil', nl_soil, 'soil 2 nitrogen density in soil layers','gN/m3')

      ! soil 3 nitrogen density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%soil3n_vr, &
         a_soil3n_vr, file_hist, 'f_soil3n_vr', itime_in_file, filter, &
         'soil', nl_soil, 'soil 3 nitrogen density in soil layers','gN/m3')

      ! coarse woody debris nitrogen density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%cwdn_vr, &
         a_cwdn_vr, file_hist, 'f_cwdn_vr', itime_in_file, filter, &
         'soil', nl_soil, 'coarse woody debris nitrogen density in soil layers','gN/m3')

      ! mineral nitrogen density in soil layers
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%sminn_vr, &
         a_sminn_vr, file_hist, 'f_sminn_vr', itime_in_file, filter, &
         'soil', nl_soil, 'mineral nitrogen density in soil layers','gN/m3')

#endif

      ! --------------------------------------------------------------------
      ! Temperature and water (excluding land water bodies and ocean patches)
      ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3; land water bodies => 4; ocean => 99]
      ! --------------------------------------------------------------------

      if (p_is_worker) then
         if (numpatch > 0) then
            filter(:) = patchtype <= 3
            IF (DEF_forcing%has_missing_value) THEN
               filter = filter .and. forcmask
            ENDIF
         end if
      end if

      ! soil temperature [K]
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%t_soisno, &
         a_t_soisno, file_hist, 'f_t_soisno', itime_in_file, filter, &
         'soilsnow', nl_soil-maxsnl, 'soil temperature','K')

      ! liquid water in soil layers [kg/m2]
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%wliq_soisno, &
         a_wliq_soisno, file_hist, 'f_wliq_soisno', itime_in_file, filter,&
         'soilsnow', nl_soil-maxsnl, 'liquid water in soil layers','kg/m2')

      ! ice lens in soil layers [kg/m2]
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%wliq_soisno, &
         a_wice_soisno, file_hist, 'f_wice_soisno', itime_in_file, filter,&
         'soilsnow', nl_soil-maxsnl, 'ice lens in soil layers','kg/m2')

      ! --------------------------------------------------------------------
      ! additial diagnostic variables for output (vegetated land only <=2)
      ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3; land water bodies => 4; ocean => 99]
      ! --------------------------------------------------------------------

      if (p_is_worker) then
         if (numpatch > 0) then
            filter(:) = patchtype <= 2
            IF (DEF_forcing%has_missing_value) THEN
               filter = filter .and. forcmask
            ENDIF
         end if
      end if

      ! volumetric soil water in layers [m3/m3]
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%h2osoi, &
         a_h2osoi, file_hist, 'f_h2osoi', itime_in_file, filter, &
         'soil', nl_soil, 'volumetric water in soil layers','m3/m3')

      ! fraction of root water uptake from each soil layer, all layers add to 1, when PHS is not defined
      ! water exchange between soil layers and root. Positive: soil->root [mm h2o/s], when PHS is defined
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%rootr, &
         a_rootr, file_hist, 'f_rootr', itime_in_file, filter, &
         'soil', nl_soil, 'root water uptake', 'mm h2o/s')

#ifdef PLANT_HYDRAULIC_STRESS
      ! vegetation water potential [mm]
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%vegwp, &
         a_vegwp, file_hist, 'f_vegwp', itime_in_file, filter, &
         'vegnodes', nvegwcs, 'vegetation water potential', 'mm')
#endif

      ! water table depth [m]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%zwt, &
         a_zwt, file_hist, 'f_zwt', itime_in_file, filter, &
         'the depth to water table','m')

      ! water storage in aquifer [mm]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%wa, &
         a_wa, file_hist, 'f_wa', itime_in_file, filter, &
         'water storage in aquifer','mm')

      if (p_is_worker) then
         if (numpatch > 0) then
            filter(:) = .true.
            IF (DEF_forcing%has_missing_value) THEN
               filter = filter .and. forcmask
            ENDIF
         end if
      end if
      ! depth of ponding water [m]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%dpond, &
         a_dpond, file_hist, 'f_dpond', itime_in_file, filter, &
         'depth of ponding water','mm')

      ! -----------------------------------------------
      ! Land water bodies' ice fraction and temperature
      ! -----------------------------------------------

      if (p_is_worker) then
         if (numpatch > 0) then
            filter(:) = patchtype == 4
            IF (DEF_forcing%has_missing_value) THEN
               filter = filter .and. forcmask
            ENDIF
         end if
      end if

      ! lake temperature [K]
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%t_lake, &
         a_t_lake, file_hist, 'f_t_lake', itime_in_file, filter, &
         'lake', nl_lake, 'lake temperature','K')

      ! lake ice fraction cover [0-1]
      call aggregate_to_vector_and_write_3d ( DEF_hist_vars%lake_icefrac, &
         a_lake_icefrac, file_hist, 'f_lake_icefrac', itime_in_file, filter, &
         'lake', nl_lake, 'lake ice fraction cover','0-1')

      ! --------------------------------
      ! Retrieve through averaged fluxes
      ! --------------------------------
      if (p_is_worker) then
         if (numpatch > 0) then
            filter(:) = patchtype < 99
            IF (DEF_forcing%has_missing_value) THEN
               filter = filter .and. forcmask
            ENDIF
         end if
      end if

      ! u* in similarity theory [m/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%ustar, &
         a_ustar, file_hist, 'f_ustar', itime_in_file, filter, &
         'u* in similarity theory','m/s')

      ! t* in similarity theory [kg/kg]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%tstar, &
         a_tstar, file_hist, 'f_tstar', itime_in_file, filter, &
         't* in similarity theory','kg/kg')

      ! q* in similarity theory [kg/kg]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%qstar, &
         a_qstar, file_hist, 'f_qstar', itime_in_file, filter, &
         'q* in similarity theory', 'kg/kg')

      ! dimensionless height (z/L) used in Monin-Obukhov theory
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%zol, &
         a_zol, file_hist, 'f_zol', itime_in_file, filter, &
         'dimensionless height (z/L) used in Monin-Obukhov theory','-')

      ! bulk Richardson number in surface layer
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%rib, &
         a_rib, file_hist, 'f_rib', itime_in_file, filter, &
         'bulk Richardson number in surface layer','-')

      ! integral of profile function for momentum
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fm, &
         a_fm, file_hist, 'f_fm', itime_in_file, filter, &
         'integral of profile function for momentum','-')

      ! integral of profile function for heat
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fh, &
         a_fh, file_hist, 'f_fh', itime_in_file, filter, &
         'integral of profile function for heat','-')

      ! integral of profile function for moisture
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fq, &
         a_fq, file_hist, 'f_fq', itime_in_file, filter, &
         'integral of profile function for moisture','-')

      ! 10m u-velocity [m/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%us10m, &
         a_us10m, file_hist, 'f_us10m', itime_in_file, filter, &
         '10m u-velocity','m/s')

      ! 10m v-velocity [m/s]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%vs10m, &
         a_vs10m, file_hist, 'f_vs10m', itime_in_file, filter, &
         '10m v-velocity','m/s')

      ! integral of profile function for momentum at 10m [-]
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%fm10m, &
         a_fm10m, file_hist, 'f_fm10m', itime_in_file, filter, &
         'integral of profile function for momentum at 10m','-')

      ! total reflected solar radiation (W/m2)
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%sr, &
         a_sr, file_hist, 'f_sr', itime_in_file, filter, &
         'reflected solar radiation at surface [W/m2]','W/m2')

      ! incident direct beam vis solar radiation (W/m2)
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%solvd, &
         a_solvd, file_hist, 'f_solvd', itime_in_file, filter, &
         'incident direct beam vis solar radiation (W/m2)','W/m2')

      ! incident diffuse beam vis solar radiation (W/m2)
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%solvi, &
         a_solvi, file_hist, 'f_solvi', itime_in_file, filter, &
         'incident diffuse beam vis solar radiation (W/m2)','W/m2')

      ! incident direct beam nir solar radiation (W/m2)
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%solnd, &
         a_solnd, file_hist, 'f_solnd', itime_in_file, filter, &
         'incident direct beam nir solar radiation (W/m2)','W/m2')

      ! incident diffuse beam nir solar radiation (W/m2)
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%solni, &
         a_solni, file_hist, 'f_solni', itime_in_file, filter, &
         'incident diffuse beam nir solar radiation (W/m2)','W/m2')

      ! reflected direct beam vis solar radiation (W/m2)
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%srvd, &
         a_srvd, file_hist, 'f_srvd', itime_in_file, filter, &
         'reflected direct beam vis solar radiation (W/m2)','W/m2')

      ! reflected diffuse beam vis solar radiation (W/m2)
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%srvi, &
         a_srvi, file_hist, 'f_srvi', itime_in_file, filter, &
         'reflected diffuse beam vis solar radiation (W/m2)','W/m2')

      ! reflected direct beam nir solar radiation (W/m2)
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%srnd, &
         a_srnd, file_hist, 'f_srnd', itime_in_file, filter, &
         'reflected direct beam nir solar radiation (W/m2)','W/m2')

      ! reflected diffuse beam nir solar radiation (W/m2)
      call aggregate_to_vector_and_write_2d ( DEF_hist_vars%srni, &
         a_srni, file_hist, 'f_srni', itime_in_file, filter, &
         'reflected diffuse beam nir solar radiation (W/m2)','W/m2')

      ! local noon fluxes
      if (p_is_worker) then
         if (numpatch > 0) then
            filter(:) = nac_ln > 0
            IF (DEF_forcing%has_missing_value) THEN
               filter = filter .and. forcmask
            ENDIF
         end if
      end if

      ! incident direct beam vis solar radiation at local noon (W/m2)
      call aggregate_to_vector_and_write_ln ( DEF_hist_vars%solvdln, &
         a_solvdln, file_hist, 'f_solvdln', itime_in_file, filter, &
         'incident direct beam vis solar radiation at local noon(W/m2)','W/m2')

      ! incident diffuse beam vis solar radiation at local noon (W/m2)
      call aggregate_to_vector_and_write_ln ( DEF_hist_vars%solviln, &
         a_solviln, file_hist, 'f_solviln', itime_in_file, filter, &
         'incident diffuse beam vis solar radiation at local noon(W/m2)','W/m2')

      ! incident direct beam nir solar radiation at local noon (W/m2)
      call aggregate_to_vector_and_write_ln ( DEF_hist_vars%solndln, &
         a_solndln, file_hist, 'f_solndln', itime_in_file, filter, &
         'incident direct beam nir solar radiation at local noon(W/m2)','W/m2')

      ! incident diffuse beam nir solar radiation at local noon (W/m2)
      call aggregate_to_vector_and_write_ln ( DEF_hist_vars%solniln, &
         a_solniln, file_hist, 'f_solniln', itime_in_file, filter, &
         'incident diffuse beam nir solar radiation at local noon(W/m2)','W/m2')

      ! reflected direct beam vis solar radiation at local noon (W/m2)
      call aggregate_to_vector_and_write_ln ( DEF_hist_vars%srvdln, &
         a_srvdln, file_hist, 'f_srvdln', itime_in_file, filter, &
         'reflected direct beam vis solar radiation at local noon(W/m2)','W/m2')

      ! reflected diffuse beam vis solar radiation at local noon (W/m2)
      call aggregate_to_vector_and_write_ln ( DEF_hist_vars%srviln, &
         a_srviln, file_hist, 'f_srviln', itime_in_file, filter, &
         'reflected diffuse beam vis solar radiation at local noon(W/m2)','W/m2')

      ! reflected direct beam nir solar radiation at local noon (W/m2)
      call aggregate_to_vector_and_write_ln ( DEF_hist_vars%srndln, &
         a_srndln, file_hist, 'f_srndln', itime_in_file, filter, &
         'reflected direct beam nir solar radiation at local noon(W/m2)','W/m2')

      ! reflected diffuse beam nir solar radiation at local noon(W/m2)
      call aggregate_to_vector_and_write_ln ( DEF_hist_vars%srniln, &
         a_srniln, file_hist, 'f_srniln', itime_in_file, filter, &
         'reflected diffuse beam nir solar radiation at local noon(W/m2)','W/m2')

      if (allocated(filter)) deallocate(filter)

      call FLUSH_acc_fluxes ()

   END SUBROUTINE hist_vector_out

   ! -------
   subroutine aggregate_to_vector_and_write_2d ( is_hist, &
         acc_vec_patch, file_hist, varname, itime_in_file, filter, &
         longname, units)

      use precision
      use spmd_task
      use mod_namelist
      USE mod_landpatch
      use MOD_Vars_1DAccFluxes,  only: nac
      use MOD_Vars_Global, only: spval
      implicit none

      logical, intent(in) :: is_hist

      real(r8), intent(in) :: acc_vec_patch (:)
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer,          intent(in) :: itime_in_file
      logical, intent(in) :: filter(:)

      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units

      ! Local variables
      INTEGER :: numset, totalnumset, iset, istt, iend, iwork, mesg(2), isrc, ndata, compress
      LOGICAL,  allocatable :: mask(:)
      REAL(r8), allocatable :: frac(:)
      REAL(r8), allocatable :: acc_vec(:), rcache(:)
      REAL(r8) :: sumwt

      if (.not. is_hist) return

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      if (p_is_worker) then
#ifdef CATCHMENT
         numset = numhru
#else
         numset = numelm
#endif

         IF (numset > 0) THEN

            allocate (acc_vec (numset))
            acc_vec(:) = spval

            DO iset = 1, numset
#ifdef CATCHMENT
               istt = hru_patch%substt(iset)
               iend = hru_patch%subend(iset)
#else
               istt = elm_patch%substt(iset)
               iend = elm_patch%subend(iset)
#endif

               IF ((istt > 0) .and. (iend >= istt)) THEN
                  allocate (mask(istt:iend))
                  allocate (frac(istt:iend))
                  mask = (acc_vec_patch(istt:iend) /= spval) .and. filter(istt:iend)
                  IF (any(mask)) THEN
#ifdef CATCHMENT
                     frac = hru_patch%subfrc(istt:iend)
#else
                     frac = elm_patch%subfrc(istt:iend)
#endif
                     sumwt = sum(frac, mask = mask)
                     acc_vec(iset) = sum(frac * acc_vec_patch(istt:iend), mask = mask)
                     acc_vec(iset) = acc_vec(iset) / sumwt / nac
                  ENDIF
                  deallocate(mask)
                  deallocate(frac)
               ENDIF
            ENDDO
         end if

#ifdef USEMPI
         mesg = (/p_iam_glb, numset/)
         call mpi_send (mesg, 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numset > 0) THEN
            call mpi_send (acc_vec, numset, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN

#ifdef CATCHMENT
         totalnumset = totalnumhru
#else
         totalnumset = totalnumelm
#endif

         IF (.not. allocated(acc_vec)) THEN
            allocate (acc_vec (totalnumset))
         ENDIF

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            call mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndata))
               call mpi_recv (rcache, ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

#ifdef CATCHMENT
               acc_vec(hru_data_address(p_itis_worker(isrc))%val) = rcache
#else
               acc_vec(elm_data_address(p_itis_worker(isrc))%val) = rcache
#endif

               deallocate (rcache)
            ENDIF
         ENDDO
#else
#ifdef CATCHMENT
         acc_vec(hru_data_address(0)%val) = acc_vec
#else
         acc_vec(elm_data_address(0)%val) = acc_vec
#endif
#endif
      ENDIF

      IF (p_is_master) THEN

         compress = DEF_HIST_COMPRESS_LEVEL
#ifdef CATCHMENT
         call ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
            'hydrounit', 'time', compress)
#else
         call ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
            'element', 'time', compress)
#endif

         IF (itime_in_file == 1) then
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF

      ENDIF

      IF (allocated(acc_vec)) deallocate (acc_vec)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   end subroutine aggregate_to_vector_and_write_2d

   ! -------
   subroutine aggregate_to_vector_and_write_3d ( is_hist, &
         acc_vec_patch, file_hist, varname, itime_in_file, filter, &
         dim1name, ndim1, &
         longname, units)

      use precision
      use spmd_task
      use mod_namelist
      USE mod_landpatch
      use MOD_Vars_1DAccFluxes,  only: nac
      use MOD_Vars_Global, only: spval
      implicit none

      logical, intent(in) :: is_hist

      real(r8), intent(in) :: acc_vec_patch (:,:)
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer,          intent(in) :: itime_in_file
      logical, intent(in) :: filter(:)

      character(len=*), intent(in) :: dim1name
      integer,          intent(in) :: ndim1

      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units

      ! Local variables
      INTEGER :: numset, totalnumset, iset, istt, iend, iwork, mesg(2), isrc, ndata, compress
      INTEGER :: lb1, ub1, i1
      LOGICAL,  allocatable :: mask(:)
      REAL(r8), allocatable :: frac(:)
      REAL(r8), allocatable :: acc_vec(:,:), rcache(:,:)
      REAL(r8) :: sumwt

      if (.not. is_hist) return

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      if (p_is_worker) then
#ifdef CATCHMENT
         numset = numhru
#else
         numset = numelm
#endif

         IF (numset > 0) THEN

            IF (numpatch > 0) THEN
               lb1 = lbound(acc_vec_patch,1)
               ub1 = ubound(acc_vec_patch,1)
            ELSE
               lb1 = 1
               ub1 = ndim1
            ENDIF

            allocate (acc_vec (lb1:ub1,numset))

            acc_vec(:,:) = spval

            DO iset = 1, numset
#ifdef CATCHMENT
               istt = hru_patch%substt(iset)
               iend = hru_patch%subend(iset)
#else
               istt = elm_patch%substt(iset)
               iend = elm_patch%subend(iset)
#endif

               IF ((istt > 0) .and. (iend >= istt)) THEN
                  allocate (mask(istt:iend))
                  allocate (frac(istt:iend))
                  DO i1 = lb1, ub1
                     mask = (acc_vec_patch(i1,istt:iend) /= spval) .and. filter(istt:iend)
                     IF (any(mask)) THEN
#ifdef CATCHMENT
                        frac = hru_patch%subfrc(istt:iend)
#else
                        frac = elm_patch%subfrc(istt:iend)
#endif
                        sumwt = sum(frac, mask = mask)
                        acc_vec(i1,iset) = sum(frac * acc_vec_patch(i1,istt:iend), mask = mask)
                        acc_vec(i1,iset) = acc_vec(i1,iset) / sumwt / nac
                     ENDIF
                  ENDDO
                  deallocate(mask)
                  deallocate(frac)
               ENDIF
            ENDDO
         end if

#ifdef USEMPI
         mesg = (/p_iam_glb, numset/)
         call mpi_send (mesg, 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numset > 0) THEN
            call mpi_send (acc_vec, ndim1 * numset, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN

#ifdef CATCHMENT
         totalnumset = totalnumhru
#else
         totalnumset = totalnumelm
#endif

         IF (.not. allocated(acc_vec)) THEN
            allocate (acc_vec (ndim1,totalnumset))
         ENDIF

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            call mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndim1,ndata))
               call mpi_recv (rcache, ndim1*ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO i1 = 1, ndim1
#ifdef CATCHMENT
                  acc_vec(i1,hru_data_address(p_itis_worker(isrc))%val) = rcache(i1,:)
#else
                  acc_vec(i1,elm_data_address(p_itis_worker(isrc))%val) = rcache(i1,:)
#endif

               ENDDO

               deallocate (rcache)
            ENDIF
         ENDDO
#else
         DO i1 = lb1, ub1
#ifdef CATCHMENT
            acc_vec(i1,hru_data_address(0)%val) = acc_vec(i1,:)
#else
            acc_vec(i1,elm_data_address(0)%val) = acc_vec(i1,:)
#endif
         ENDDO
#endif
      ENDIF

      IF (p_is_master) THEN

         call ncio_define_dimension (file_hist, dim1name, ndim1)

         compress = DEF_HIST_COMPRESS_LEVEL
#ifdef CATCHMENT
         call ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
            dim1name, 'hydrounit', 'time', compress)
#else
         call ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
            dim1name, 'element', 'time', compress)
#endif

         IF (itime_in_file == 1) then
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF

      ENDIF

      IF (allocated(acc_vec)) deallocate (acc_vec)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   end subroutine aggregate_to_vector_and_write_3d

   ! -------
   subroutine aggregate_to_vector_and_write_4d ( is_hist, &
         acc_vec_patch, file_hist, varname, itime_in_file, filter, &
         dim1name, ndim1, dim2name, ndim2, &
         longname, units)

      use precision
      use spmd_task
      use mod_namelist
      USE mod_landpatch
      use MOD_Vars_1DAccFluxes,  only: nac
      use MOD_Vars_Global, only: spval
      implicit none

      logical, intent(in) :: is_hist

      real(r8), intent(in) :: acc_vec_patch (:,:,:)
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer,          intent(in) :: itime_in_file
      logical, intent(in) :: filter(:)

      character(len=*), intent(in) :: dim1name, dim2name
      integer,          intent(in) :: ndim1, ndim2

      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units

      ! Local variables
      INTEGER :: numset, totalnumset, iset, istt, iend, iwork, mesg(2), isrc, ndata, compress
      INTEGER :: lb1, ub1, i1, lb2, ub2, i2
      LOGICAL,  allocatable :: mask(:)
      REAL(r8), allocatable :: frac(:)
      REAL(r8), allocatable :: acc_vec(:,:,:), rcache(:,:,:)
      REAL(r8) :: sumwt

      if (.not. is_hist) return

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      if (p_is_worker) then
#ifdef CATCHMENT
         numset = numhru
#else
         numset = numelm
#endif

         IF (numset > 0) THEN

            IF (numpatch > 0) THEN
               lb1 = lbound(acc_vec_patch,1)
               ub1 = ubound(acc_vec_patch,1)
               lb2 = lbound(acc_vec_patch,2)
               ub2 = ubound(acc_vec_patch,2)
            ELSE
               lb1 = 1
               ub1 = ndim1
               lb2 = 1
               ub2 = ndim2
            ENDIF

            allocate (acc_vec (lb1:ub1,lb2:ub2,numset))

            acc_vec(:,:,:) = spval

            DO iset = 1, numset
#ifdef CATCHMENT
               istt = hru_patch%substt(iset)
               iend = hru_patch%subend(iset)
#else
               istt = elm_patch%substt(iset)
               iend = elm_patch%subend(iset)
#endif

               IF ((istt > 0) .and. (iend >= istt)) THEN
                  allocate (mask(istt:iend))
                  allocate (frac(istt:iend))
                  DO i1 = lb1, ub1
                     DO i2 = lb2, ub2
                        mask = (acc_vec_patch(i1,i2,istt:iend) /= spval) .and. filter(istt:iend)
                        IF (any(mask)) THEN
#ifdef CATCHMENT
                           frac = hru_patch%subfrc(istt:iend)
#else
                           frac = elm_patch%subfrc(istt:iend)
#endif
                           sumwt = sum(frac, mask = mask)
                           acc_vec(i1,i2,iset) = sum(frac * acc_vec_patch(i1,i2,istt:iend), mask = mask)
                           acc_vec(i1,i2,iset) = acc_vec(i1,i2,iset) / sumwt / nac
                        ENDIF
                     ENDDO
                  ENDDO
                  deallocate(mask)
                  deallocate(frac)
               ENDIF
            ENDDO
         end if

#ifdef USEMPI
         mesg = (/p_iam_glb, numset/)
         call mpi_send (mesg, 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numset > 0) THEN
            call mpi_send (acc_vec, ndim1 * ndim2 * numset, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN

#ifdef CATCHMENT
         totalnumset = totalnumhru
#else
         totalnumset = totalnumelm
#endif

         IF (.not. allocated(acc_vec)) THEN
            allocate (acc_vec (ndim1,ndim2,totalnumset))
         ENDIF

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            call mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndim1,ndim2,ndata))
               call mpi_recv (rcache, ndim1 * ndim2 * ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO i1 = 1, ndim1
                  DO i2 = 1, ndim2
#ifdef CATCHMENT
                     acc_vec(i1,i2,hru_data_address(p_itis_worker(isrc))%val) = rcache(i1,i2,:)
#else
                     acc_vec(i1,i2,elm_data_address(p_itis_worker(isrc))%val) = rcache(i1,i2,:)
#endif
                  ENDDO
               ENDDO

               deallocate (rcache)
            ENDIF
         ENDDO
#else
         DO i1 = lb1, ub1
            DO i2 = lb2, ub2
#ifdef CATCHMENT
               acc_vec(i1,i2,hru_data_address(0)%val) = acc_vec(i1,i2,:)
#else
               acc_vec(i1,i2,elm_data_address(0)%val) = acc_vec(i1,i2,:)
#endif
            ENDDO
         ENDDO
#endif
      ENDIF

      IF (p_is_master) THEN

         call ncio_define_dimension (file_hist, dim1name, ndim1)
         call ncio_define_dimension (file_hist, dim2name, ndim2)

         compress = DEF_HIST_COMPRESS_LEVEL
#ifdef CATCHMENT
         call ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
            dim1name, dim2name, 'hydrounit', 'time', compress)
#else
         call ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
            dim1name, dim2name, 'element', 'time', compress)
#endif

         IF (itime_in_file == 1) then
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF

      ENDIF

      IF (allocated(acc_vec)) deallocate (acc_vec)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   end subroutine aggregate_to_vector_and_write_4d

   ! -------
   subroutine aggregate_to_vector_and_write_ln ( is_hist, &
         acc_vec_patch, file_hist, varname, itime_in_file, filter, &
         longname, units)

      use precision
      use spmd_task
      use mod_namelist
      USE mod_landpatch
      use MOD_Vars_1DAccFluxes,  only: nac_ln
      use MOD_Vars_Global, only: spval
      implicit none

      logical, intent(in) :: is_hist

      real(r8), intent(in) :: acc_vec_patch (:)
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer,          intent(in) :: itime_in_file
      logical, intent(in) :: filter(:)

      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units

      ! Local variables
      INTEGER :: numset, totalnumset, iset, istt, iend, iwork, mesg(2), isrc, ndata, compress
      LOGICAL,  allocatable :: mask(:)
      REAL(r8), allocatable :: frac(:)
      REAL(r8), allocatable :: acc_vec(:), rcache(:)
      REAL(r8) :: sumwt

      if (.not. is_hist) return

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      if (p_is_worker) then
#ifdef CATCHMENT
         numset = numhru
#else
         numset = numelm
#endif

         IF (numset > 0) THEN

            allocate (acc_vec (numset))
            acc_vec(:) = spval

            DO iset = 1, numset
#ifdef CATCHMENT
               istt = hru_patch%substt(iset)
               iend = hru_patch%subend(iset)
#else
               istt = elm_patch%substt(iset)
               iend = elm_patch%subend(iset)
#endif

               IF ((istt > 0) .and. (iend >= istt)) THEN
                  allocate (mask(istt:iend))
                  allocate (frac(istt:iend))
                  mask = (acc_vec_patch(istt:iend) /= spval) .and. filter(istt:iend) .and. (nac_ln > 0)
                  IF (any(mask)) THEN
#ifdef CATCHMENT
                     frac = hru_patch%subfrc(istt:iend)
#else
                     frac = elm_patch%subfrc(istt:iend)
#endif
                     sumwt = sum(frac, mask = mask)
                     acc_vec(iset) = sum(frac * acc_vec_patch(istt:iend) / nac_ln(istt:iend), mask = mask)
                     acc_vec(iset) = acc_vec(iset) / sumwt
                  ENDIF
                  deallocate(mask)
                  deallocate(frac)
               ENDIF
            ENDDO
         end if

#ifdef USEMPI
         mesg = (/p_iam_glb, numset/)
         call mpi_send (mesg, 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numset > 0) THEN
            call mpi_send (acc_vec, numset, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN

#ifdef CATCHMENT
         totalnumset = totalnumhru
#else
         totalnumset = totalnumelm
#endif

         IF (.not. allocated(acc_vec)) THEN
            allocate (acc_vec (totalnumset))
         ENDIF

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            call mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndata))
               call mpi_recv (rcache, ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

#ifdef CATCHMENT
               acc_vec(hru_data_address(p_itis_worker(isrc))%val) = rcache
#else
               acc_vec(elm_data_address(p_itis_worker(isrc))%val) = rcache
#endif

               deallocate (rcache)
            ENDIF
         ENDDO
#else
#ifdef CATCHMENT
         acc_vec(hru_data_address(0)%val) = acc_vec
#else
         acc_vec(elm_data_address(0)%val) = acc_vec
#endif
#endif
      ENDIF

      IF (p_is_master) THEN

         compress = DEF_HIST_COMPRESS_LEVEL
#ifdef CATCHMENT
         call ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
            'hydrounit', 'time', compress)
#else
         call ncio_write_serial_time (file_hist, varname, itime_in_file, acc_vec, &
            'element', 'time', compress)
#endif

         IF (itime_in_file == 1) then
            CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
            CALL ncio_put_attr (file_hist, varname, 'units', units)
            CALL ncio_put_attr (file_hist, varname, 'missing_value', spval)
         ENDIF

      ENDIF

      IF (allocated(acc_vec)) deallocate (acc_vec)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   end subroutine aggregate_to_vector_and_write_ln

end module MOD_HistVector
#endif
