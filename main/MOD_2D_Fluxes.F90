#include <define.h>

MODULE MOD_2D_Fluxes
! ----------------------------------------------------------------------
! perfrom the grid average mapping: average a subgrid input 1d vector 
! of length numpatch to a output 2d array of length [lon_points,lat_points]
!
! Created by Yongjiu Dai, 03/2014
!---------------------------------------------------------------------

   use mod_data_type
   USE GlobalVars
#ifdef BGC
   USE MOD_2D_BGCFluxes
#endif

   IMPLICIT NONE
   SAVE

   type(block_data_real8_2d) :: f_taux    ! wind stress: E-W [kg/m/s2]
   type(block_data_real8_2d) :: f_tauy    ! wind stress: N-S [kg/m/s2]
   type(block_data_real8_2d) :: f_fsena   ! sensible heat from canopy height to atmosphere [W/m2]
   type(block_data_real8_2d) :: f_lfevpa  ! latent heat flux from canopy height to atmosphere [W/m2]
   type(block_data_real8_2d) :: f_fevpa   ! evapotranspiration from canopy to atmosphere [mm/s]
   type(block_data_real8_2d) :: f_fsenl   ! sensible heat from leaves [W/m2]
   type(block_data_real8_2d) :: f_fevpl   ! evaporation+transpiration from leaves [mm/s]
   type(block_data_real8_2d) :: f_etr     ! transpiration rate [mm/s]
   type(block_data_real8_2d) :: f_fseng   ! sensible heat flux from ground [W/m2]
   type(block_data_real8_2d) :: f_fevpg   ! evaporation heat flux from ground [mm/s]
   type(block_data_real8_2d) :: f_fgrnd   ! ground heat flux [W/m2]
   type(block_data_real8_2d) :: f_sabvsun ! solar absorbed by sunlit canopy [W/m2]
   type(block_data_real8_2d) :: f_sabvsha ! solar absorbed by shaded [W/m2]
   type(block_data_real8_2d) :: f_sabg    ! solar absorbed by ground  [W/m2]
   type(block_data_real8_2d) :: f_sr      ! total reflected solar radiation (W/m2)
   type(block_data_real8_2d) :: f_solvd   ! incident direct beam vis solar radiation (W/m2)
   type(block_data_real8_2d) :: f_solvi   ! incident diffuse beam vis solar radiation (W/m2)
   type(block_data_real8_2d) :: f_solnd   ! incident direct beam nir solar radiation (W/m2)
   type(block_data_real8_2d) :: f_solni   ! incident diffuse beam nir solar radiation (W/m2)
   type(block_data_real8_2d) :: f_srvd    ! reflected direct beam vis solar radiation (W/m2)
   type(block_data_real8_2d) :: f_srvi    ! reflected diffuse beam vis solar radiation (W/m2)
   type(block_data_real8_2d) :: f_srnd    ! reflected direct beam nir solar radiation (W/m2)
   type(block_data_real8_2d) :: f_srni    ! reflected diffuse beam nir solar radiation (W/m2)
   type(block_data_real8_2d) :: f_solvdln ! incident direct beam vis solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_solviln ! incident diffuse beam vis solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_solndln ! incident direct beam nir solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_solniln ! incident diffuse beam nir solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_srvdln  ! reflected direct beam vis solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_srviln  ! reflected diffuse beam vis solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_srndln  ! reflected direct beam nir solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_srniln  ! reflected diffuse beam nir solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_olrg    ! outgoing long-wave radiation from ground+canopy [W/m2]
   type(block_data_real8_2d) :: f_rnet    ! net radiation [W/m2]
   type(block_data_real8_2d) :: f_xerr    ! the error of water banace [mm/s]
   type(block_data_real8_2d) :: f_zerr    ! the error of energy balance [W/m2]
   type(block_data_real8_2d) :: f_rsur    ! surface runoff [mm/s]
   type(block_data_real8_2d) :: f_rsub    ! subsurface runoff [mm/s]
   type(block_data_real8_2d) :: f_rnof    ! total runoff [mm/s]
   type(block_data_real8_2d) :: f_qintr   ! interception [mm/s]
   type(block_data_real8_2d) :: f_qinfl   ! inflitration [mm/s]
   type(block_data_real8_2d) :: f_qdrip   ! throughfall [mm/s]
   type(block_data_real8_2d) :: f_assim   ! canopy assimilation rate [mol m-2 s-1]
   type(block_data_real8_2d) :: f_respc   ! respiration (plant+soil) [mol m-2 s-1]
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
   type(block_data_real8_2d) :: f_assim_RuBP_sun        
   type(block_data_real8_2d) :: f_assim_RuBP_sha        
   type(block_data_real8_2d) :: f_assim_Rubisco_sun        
   type(block_data_real8_2d) :: f_assim_Rubisco_sha        
   type(block_data_real8_2d) :: f_assimsun        
   type(block_data_real8_2d) :: f_assimsha        
   type(block_data_real8_2d) :: f_etrsun        
   type(block_data_real8_2d) :: f_etrsha        
   type(block_data_real8_2d) :: f_cisun        
   type(block_data_real8_2d) :: f_cisha        
   type(block_data_real8_2d) :: f_Dsun        
   type(block_data_real8_2d) :: f_Dsha        
   type(block_data_real8_2d) :: f_gammasun        
   type(block_data_real8_2d) :: f_gammasha        
   type(block_data_real8_2d) :: f_lambdasun      
   type(block_data_real8_2d) :: f_lambdasha      
   type(block_data_real8_2d) :: f_lambda                 
#endif
#endif
   type(block_data_real8_2d) :: f_qcharge ! groundwater recharge rate [mm/s] 

   !---------------------------------------------------------------------
   type(block_data_real8_2d) :: f_t_grnd   ! ground surface temperature [K]
   type(block_data_real8_2d) :: f_tleaf    ! sunlit leaf temperature [K]
   type(block_data_real8_2d) :: f_ldew     ! depth of water on foliage [mm]
!#ifdef CLM5_INTERCEPTION
   type(block_data_real8_2d) :: f_ldew_rain     ! depth of rain on foliage [mm]
   type(block_data_real8_2d) :: f_ldew_snow     ! depth of snow on foliage [mm]
!#endif
   type(block_data_real8_2d) :: f_scv      ! snow cover, water equivalent [mm]
   type(block_data_real8_2d) :: f_snowdp   ! snow depth [meter]
   type(block_data_real8_2d) :: f_fsno     ! fraction of snow cover on ground
   type(block_data_real8_2d) :: f_sigf     ! fraction of veg cover, excluding snow-covered veg [-]
   type(block_data_real8_2d) :: f_green    ! leaf greenness
   type(block_data_real8_2d) :: f_lai      ! leaf area index
   type(block_data_real8_2d) :: f_laisun   ! sunlit leaf area index
   type(block_data_real8_2d) :: f_laisha   ! shaded leaf area index
   type(block_data_real8_2d) :: f_sai      ! stem area index
   type(block_data_real8_4d) :: f_alb      ! averaged albedo [visible, direct; direct, diffuse]
   type(block_data_real8_2d) :: f_emis     ! averaged bulk surface emissivity
   type(block_data_real8_2d) :: f_z0m      ! effective roughness [m]
   type(block_data_real8_2d) :: f_trad     ! radiative temperature of surface [K]
   type(block_data_real8_2d) :: f_tref     ! 2 m height air temperature [kelvin]
   type(block_data_real8_2d) :: f_qref     ! 2 m height air specific humidity [kg/kg]

   !---------------------------------------------------------------------
   type(block_data_real8_3d) :: f_t_soisno     ! soil temperature [K]
   type(block_data_real8_3d) :: f_wliq_soisno  ! liquid water in soil layers [kg/m2]
   type(block_data_real8_3d) :: f_wice_soisno  ! ice lens in soil layers [kg/m2]
   type(block_data_real8_3d) :: f_h2osoi       ! volumetric soil water in layers [m3/m3]
   type(block_data_real8_3d) :: f_rootr        ! water exchange between soil layers and root
   type(block_data_real8_3d) :: f_BD_all       ! bulk density in soil layers [kg/m3]
   type(block_data_real8_3d) :: f_wfc          ! water field capacity [m3/m3]
   type(block_data_real8_3d) :: f_OM_density   ! soil organic matter density [kg/m3]
#ifdef PLANT_HYDRAULIC_STRESS
   type(block_data_real8_3d) :: f_vegwp        ! vegetation water potential [mm]
#endif
   type(block_data_real8_2d) :: f_rstfacsun    ! factor of soil water stress
   type(block_data_real8_2d) :: f_rstfacsha    ! factor of soil water stress
   type(block_data_real8_2d) :: f_gssun        ! factor of soil water stress
   type(block_data_real8_2d) :: f_gssha        ! factor of soil water stress
   type(block_data_real8_2d) :: f_dpond        ! depth of ponding water [mm]
   type(block_data_real8_2d) :: f_zwt          ! the depth to water table [m]
   type(block_data_real8_2d) :: f_wa           ! water storage in aquifer [mm]
   type(block_data_real8_2d) :: f_wat          ! total water storage [mm]

   type(block_data_real8_3d) :: f_t_lake       ! lake temperature [K]
   type(block_data_real8_3d) :: f_lake_icefrac ! lake ice fraction cover [0-1]

   !---------------------------------------------------------------------
   type(block_data_real8_2d) :: f_ustar   ! u* in similarity theory [m/s]
   type(block_data_real8_2d) :: f_tstar   ! t* in similarity theory [kg/kg]
   type(block_data_real8_2d) :: f_qstar   ! q* in similarity theory [kg/kg]
   type(block_data_real8_2d) :: f_zol     ! dimensionless height (z/L) used in Monin-Obukhov theory
   type(block_data_real8_2d) :: f_rib     ! bulk Richardson number in surface layer
   type(block_data_real8_2d) :: f_fm      ! integral of profile function for momentum
   type(block_data_real8_2d) :: f_fh      ! integral of profile function for heat
   type(block_data_real8_2d) :: f_fq      ! integral of profile function for moisture
   type(block_data_real8_2d) :: f_us10m   ! 10m u-velocity [m/s]
   type(block_data_real8_2d) :: f_vs10m   ! 10m v-velocity [m/s]
   type(block_data_real8_2d) :: f_fm10m   ! integral of profile function for momentum at 10m [-]

   !---------------------------------------------------------------------
   type(block_data_real8_2d) :: f_xy_us      ! wind in eastward direction [m/s]
   type(block_data_real8_2d) :: f_xy_vs      ! wind in northward direction [m/s]
   type(block_data_real8_2d) :: f_xy_t       ! temperature at reference height [kelvin]
   type(block_data_real8_2d) :: f_xy_q       ! specific humidity at reference height [kg/kg]
   type(block_data_real8_2d) :: f_xy_prc     ! convective precipitation [mm/s]
   type(block_data_real8_2d) :: f_xy_prl     ! large scale precipitation [mm/s]
   type(block_data_real8_2d) :: f_xy_pbot    ! atmospheric pressure at the surface [pa]
   type(block_data_real8_2d) :: f_xy_frl     ! atmospheric infrared (longwave) radiation [W/m2]
   type(block_data_real8_2d) :: f_xy_solarin ! downward solar radiation at surface [W/m2]
   type(block_data_real8_2d) :: f_xy_rain    ! rain [mm/s]
   type(block_data_real8_2d) :: f_xy_snow    ! snow [mm/s]
   type(block_data_real8_2d) :: f_xy_ozone   ! ozone concentration [mol/mol]
   
   ! PUBLIC MEMBER FUNCTIONS:
   public :: allocate_2D_Fluxes

CONTAINS

   SUBROUTINE allocate_2D_Fluxes (grid)
      ! --------------------------------------------------------------------
      ! Allocates memory for CLM 2d [lon_points,lat_points] variables
      ! --------------------------------------------------------------------

      use spmd_task
      use mod_grid
      use mod_data_type
      USE mod_namelist
      implicit none

      type(grid_type), intent(in) :: grid

#if (defined UNSTRUCTURED || defined CATCHMENT) 
      IF (DEF_HISTORY_IN_VECTOR) THEN
         RETURN
      ENDIF 
#endif

      if (p_is_io) then
         
         call allocate_block_data (grid, f_taux   )  ! wind stress: E-W [kg/m/s2]
         call allocate_block_data (grid, f_tauy   )  ! wind stress: N-S [kg/m/s2]
         call allocate_block_data (grid, f_fsena  )  ! sensible heat from canopy height to atmosphere [W/m2]
         call allocate_block_data (grid, f_lfevpa )  ! latent heat flux from canopy height to atmosphere [W/m2]
         call allocate_block_data (grid, f_fevpa  )  ! evapotranspiration from canopy to atmosphere [mm/s]
         call allocate_block_data (grid, f_fsenl  )  ! sensible heat from leaves [W/m2]
         call allocate_block_data (grid, f_fevpl  )  ! evaporation+transpiration from leaves [mm/s]
         call allocate_block_data (grid, f_etr    )  ! transpiration rate [mm/s]
         call allocate_block_data (grid, f_fseng  )  ! sensible heat flux from ground [W/m2]
         call allocate_block_data (grid, f_fevpg  )  ! evaporation heat flux from ground [mm/s]
         call allocate_block_data (grid, f_fgrnd  )  ! ground heat flux [W/m2]
         call allocate_block_data (grid, f_sabvsun)  ! solar absorbed by sunlit canopy [W/m2]
         call allocate_block_data (grid, f_sabvsha)  ! solar absorbed by shaded [W/m2]
         call allocate_block_data (grid, f_sabg   )  ! solar absorbed by ground  [W/m2]
         call allocate_block_data (grid, f_sr     )  ! total reflected solar radiation (W/m2)
         call allocate_block_data (grid, f_solvd  )  ! incident direct beam vis solar radiation (W/m2)
         call allocate_block_data (grid, f_solvi  )  ! incident diffuse beam vis solar radiation (W/m2)
         call allocate_block_data (grid, f_solnd  )  ! incident direct beam nir solar radiation (W/m2)
         call allocate_block_data (grid, f_solni  )  ! incident diffuse beam nir solar radiation (W/m2)
         call allocate_block_data (grid, f_srvd   )  ! reflected direct beam vis solar radiation (W/m2)
         call allocate_block_data (grid, f_srvi   )  ! reflected diffuse beam vis solar radiation (W/m2)
         call allocate_block_data (grid, f_srnd   )  ! reflected direct beam nir solar radiation (W/m2)
         call allocate_block_data (grid, f_srni   )  ! reflected diffuse beam nir solar radiation (W/m2)
         call allocate_block_data (grid, f_solvdln)  ! incident direct beam vis solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_solviln)  ! incident diffuse beam vis solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_solndln)  ! incident direct beam nir solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_solniln)  ! incident diffuse beam nir solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_srvdln )  ! reflected direct beam vis solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_srviln )  ! reflected diffuse beam vis solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_srndln )  ! reflected direct beam nir solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_srniln )  ! reflected diffuse beam nir solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_olrg   )  ! outgoing long-wave radiation from ground+canopy [W/m2]
         call allocate_block_data (grid, f_rnet   )  ! net radiation [W/m2]
         call allocate_block_data (grid, f_xerr   )  ! the error of water banace [mm/s]
         call allocate_block_data (grid, f_zerr   )  ! the error of energy balance [W/m2]
         call allocate_block_data (grid, f_rsur   )  ! surface runoff [mm/s]
         call allocate_block_data (grid, f_rsub   )  ! surface runoff [mm/s]
         call allocate_block_data (grid, f_rnof   )  ! total runoff [mm/s]
         call allocate_block_data (grid, f_qintr  )  ! interception [mm/s]
         call allocate_block_data (grid, f_qinfl  )  ! inflitration [mm/s]
         call allocate_block_data (grid, f_qdrip  )  ! throughfall [mm/s]
         call allocate_block_data (grid, f_assim  )  ! canopy assimilation rate [mol m-2 s-1]
         call allocate_block_data (grid, f_respc  )  ! respiration (plant+soil) [mol m-2 s-1]
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
         call allocate_block_data (grid, f_assim_RuBP_sun        )
         call allocate_block_data (grid, f_assim_RuBP_sha        )
         call allocate_block_data (grid, f_assim_Rubisco_sun        )
         call allocate_block_data (grid, f_assim_Rubisco_sha        )
         call allocate_block_data (grid, f_assimsun        )
         call allocate_block_data (grid, f_assimsha        )
         call allocate_block_data (grid, f_etrsun        )
         call allocate_block_data (grid, f_etrsha        )
         call allocate_block_data (grid, f_cisun        )
         call allocate_block_data (grid, f_cisha        )
         call allocate_block_data (grid, f_Dsun        )
         call allocate_block_data (grid, f_Dsha        )
         call allocate_block_data (grid, f_gammasun        )
         call allocate_block_data (grid, f_gammasha        )
         call allocate_block_data (grid, f_lambdasun        )
         call allocate_block_data (grid, f_lambdasha        )
         call allocate_block_data (grid, f_lambda                   )
#endif
#endif
         call allocate_block_data (grid, f_qcharge)  ! groundwater recharge rate [mm/s] 

         !---------------------------------------------------------------------
         call allocate_block_data (grid, f_t_grnd  )  ! ground surface temperature [K]
         call allocate_block_data (grid, f_tleaf   )  ! sunlit leaf temperature [K]
         call allocate_block_data (grid, f_ldew    )  ! depth of water on foliage [mm]
         call allocate_block_data (grid, f_ldew_rain    )  ! depth of rain on foliage [mm]
         call allocate_block_data (grid, f_ldew_snow    )  ! depth of snow on foliage [mm]
         call allocate_block_data (grid, f_scv     )  ! snow cover, water equivalent [mm]
         call allocate_block_data (grid, f_snowdp  )  ! snow depth [meter]
         call allocate_block_data (grid, f_fsno    )  ! fraction of snow cover on ground
         call allocate_block_data (grid, f_sigf    )  ! fraction of veg cover, excluding snow-covered veg [-]
         call allocate_block_data (grid, f_green   )  ! leaf greenness
         call allocate_block_data (grid, f_lai     )  ! leaf area index
         call allocate_block_data (grid, f_laisun  )  ! sunlit leaf area index
         call allocate_block_data (grid, f_laisha  )  ! shaded leaf area index
         call allocate_block_data (grid, f_sai     )  ! stem area index
         call allocate_block_data (grid, f_alb, 2, 2) ! averaged albedo [visible, direct; direct, diffuse]
         call allocate_block_data (grid, f_emis    )  ! averaged bulk surface emissivity
         call allocate_block_data (grid, f_z0m     )  ! effective roughness [m]
         call allocate_block_data (grid, f_trad    )  ! radiative temperature of surface [K]
         call allocate_block_data (grid, f_tref    )  ! 2 m height air temperature [kelvin]
         call allocate_block_data (grid, f_qref    )  ! 2 m height air specific humidity [kg/kg]

         !---------------------------------------------------------------------
         call allocate_block_data (grid, f_t_soisno   , nl_soil-maxsnl, maxsnl+1)  ! soil temperature [K]
         call allocate_block_data (grid, f_wliq_soisno, nl_soil-maxsnl, maxsnl+1)  ! liquid water in soil layers [kg/m2]
         call allocate_block_data (grid, f_wice_soisno, nl_soil-maxsnl, maxsnl+1)  ! ice lens in soil layers [kg/m2]
         call allocate_block_data (grid, f_h2osoi     , nl_soil)  ! volumetric soil water in layers [m3/m3]
         call allocate_block_data (grid, f_rootr      , nl_soil)  ! water exchange between soil layers and root
         call allocate_block_data (grid, f_BD_all     , nl_soil)
         call allocate_block_data (grid, f_wfc        , nl_soil)
         call allocate_block_data (grid, f_OM_density , nl_soil)
#ifdef PLANT_HYDRAULIC_STRESS
         call allocate_block_data (grid, f_vegwp      , nvegwcs)  ! vegetation water potential [mm]
#endif
         call allocate_block_data (grid, f_rstfacsun)  ! factor of soil water stress
         call allocate_block_data (grid, f_rstfacsha)  ! factor of soil water stress
         call allocate_block_data (grid, f_gssun)  ! factor of soil water stress
         call allocate_block_data (grid, f_gssha)  ! factor of soil water stress

         call allocate_block_data (grid, f_dpond  )  ! depth of ponding water [m]

         call allocate_block_data (grid, f_zwt   )  ! the depth to water table [m]
         call allocate_block_data (grid, f_wa    )  ! water storage in aquifer [mm]
         call allocate_block_data (grid, f_wat   )  ! total water storage [mm]

         call allocate_block_data (grid, f_t_lake      , nl_lake)  ! lake temperature [K]
         call allocate_block_data (grid, f_lake_icefrac, nl_lake)  ! lake ice fraction cover [0-1]

         !---------------------------------------------------------------------
         call allocate_block_data (grid, f_ustar)  ! u* in similarity theory [m/s]
         call allocate_block_data (grid, f_tstar)  ! t* in similarity theory [kg/kg]
         call allocate_block_data (grid, f_qstar)  ! q* in similarity theory [kg/kg]
         call allocate_block_data (grid, f_zol  )  ! dimensionless height (z/L) used in Monin-Obukhov theory
         call allocate_block_data (grid, f_rib  )  ! bulk Richardson number in surface layer
         call allocate_block_data (grid, f_fm   )  ! integral of profile function for momentum
         call allocate_block_data (grid, f_fh   )  ! integral of profile function for heat
         call allocate_block_data (grid, f_fq   )  ! integral of profile function for moisture
         call allocate_block_data (grid, f_us10m)  ! 10m u-velocity [m/s]
         call allocate_block_data (grid, f_vs10m)  ! 10m v-velocity [m/s]
         call allocate_block_data (grid, f_fm10m)  ! integral of profile function for momentum at 10m [-]

         !---------------------------------------------------------------------
         call allocate_block_data (grid, f_xy_us     )  ! wind in eastward direction [m/s]
         call allocate_block_data (grid, f_xy_vs     )  ! wind in northward direction [m/s]
         call allocate_block_data (grid, f_xy_t      )  ! temperature at reference height [kelvin]
         call allocate_block_data (grid, f_xy_q      )  ! specific humidity at reference height [kg/kg]
         call allocate_block_data (grid, f_xy_prc    )  ! convective precipitation [mm/s]
         call allocate_block_data (grid, f_xy_prl    )  ! large scale precipitation [mm/s]
         call allocate_block_data (grid, f_xy_pbot   )  ! atmospheric pressure at the surface [pa]
         call allocate_block_data (grid, f_xy_frl    )  ! atmospheric infrared (longwave) radiation [W/m2]
         call allocate_block_data (grid, f_xy_solarin)  ! downward solar radiation at surface [W/m2]
         call allocate_block_data (grid, f_xy_rain   )  ! rain [mm/s]
         call allocate_block_data (grid, f_xy_snow   )  ! snow [mm/s]
         call allocate_block_data (grid, f_xy_ozone  )  ! ozone concentration [mol/mol]

      end if

#ifdef BGC
      CALL allocate_2D_BGCFluxes (grid)
#endif

   END SUBROUTINE allocate_2D_Fluxes

END MODULE MOD_2D_Fluxes
