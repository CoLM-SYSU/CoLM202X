#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_TimeVariables
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Process time-varying state variables for data assimilation
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!   Lu Li, 07/2025: Remove unused variables and clean codes
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_TimeManager
   USE MOD_Namelist, only: DEF_DA_ENS_NUM
   IMPLICIT NONE
   SAVE

   ! public functions
   PUBLIC :: allocate_DATimeVariables
   PUBLIC :: deallocate_DATimeVariables
   PUBLIC :: READ_DATimeVariables
   PUBLIC :: WRITE_DATimeVariables
#ifdef RangeCheck
   PUBLIC :: check_DATimeVariables
#endif

   ! define variables
   ! Time-varying state variables which required by restart run, used to store tha noda trajectory state variable
   real(r8), allocatable :: z_sno_noda       (:,:) ! node depth [m]
   real(r8), allocatable :: dz_sno_noda      (:,:) ! interface depth [m]
   real(r8), allocatable :: t_soisno_noda    (:,:) ! soil temperature [K]
   real(r8), allocatable :: wliq_soisno_noda (:,:) ! liquid water in layers [kg/m2]
   real(r8), allocatable :: wice_soisno_noda (:,:) ! ice lens in layers [kg/m2]
   real(r8), allocatable :: smp_noda         (:,:) ! soil matrix potential [mm]
   real(r8), allocatable :: hk_noda          (:,:) ! hydraulic conductivity [mm h2o/s]
   real(r8), allocatable :: t_grnd_noda        (:) ! ground surface temperature [K]
   real(r8), allocatable :: tleaf_noda         (:) ! leaf temperature [K]
   real(r8), allocatable :: ldew_noda          (:) ! depth of water on foliage [mm]
   real(r8), allocatable :: ldew_rain_noda     (:) ! depth of rain on foliage [mm]
   real(r8), allocatable :: ldew_snow_noda     (:) ! depth of rain on foliage [mm]
   real(r8), allocatable :: fwet_snow_noda     (:) ! vegetation snow fractional cover [-]
   real(r8), allocatable :: sag_noda           (:) ! non dimensional snow age [-]
   real(r8), allocatable :: scv_noda           (:) ! snow cover, water equivalent [mm]
   real(r8), allocatable :: snowdp_noda        (:) ! snow depth [meter]
   real(r8), allocatable :: fveg_noda          (:) ! fraction of vegetation cover
   real(r8), allocatable :: fsno_noda          (:) ! fraction of snow cover on ground
   real(r8), allocatable :: sigf_noda          (:) ! fraction of veg cover, excluding snow-covered veg [-]
   real(r8), allocatable :: green_noda         (:) ! leaf greenness
   real(r8), allocatable :: tlai_noda          (:) ! leaf area index
   real(r8), allocatable :: lai_noda           (:) ! leaf area index
   real(r8), allocatable :: tsai_noda          (:) ! stem area index
   real(r8), allocatable :: sai_noda           (:) ! stem area index
   real(r8), allocatable :: alb_noda       (:,:,:) ! averaged albedo [-]
   real(r8), allocatable :: ssun_noda      (:,:,:) ! sunlit canopy absorption for solar radiation (0-1)
   real(r8), allocatable :: ssha_noda      (:,:,:) ! shaded canopy absorption for solar radiation (0-1)
   real(r8), allocatable :: ssoi_noda      (:,:,:) ! soil absorption for solar radiation (0-1)
   real(r8), allocatable :: ssno_noda      (:,:,:) ! snow absorption for solar radiation (0-1)
   real(r8), allocatable :: thermk_noda        (:) ! canopy gap fraction for tir radiation
   real(r8), allocatable :: extkb_noda         (:) ! (k, g(mu)/mu) direct solar extinction coefficient
   real(r8), allocatable :: extkd_noda         (:) ! diffuse and scattered diffuse PAR extinction coefficient
   real(r8), allocatable :: zwt_noda           (:) ! the depth to water table [m]
   real(r8), allocatable :: wdsrf_noda         (:) ! depth of surface water [mm]
   real(r8), allocatable :: wa_noda            (:) ! water storage in aquifer [mm]
   real(r8), allocatable :: wetwat_noda        (:) ! water storage in wetland [mm]
   real(r8), allocatable :: t_lake_noda      (:,:) ! lake layer teperature [K]
   real(r8), allocatable :: lake_icefrac_noda(:,:) ! lake mass fraction of lake layer that is frozen
   real(r8), allocatable :: savedtke1_noda     (:) ! top level eddy conductivity (W/m K)

   ! diagnostic variables for RTM forward operator
   real(r8), allocatable :: tref_noda          (:) ! 2 m height air temperature [kelvin]
   real(r8), allocatable :: h2osoi_noda      (:,:) ! volumetric soil water in layers [m3/m3]

   ! Time-varying state variables which required by restart run
   real(r8), allocatable :: z_sno_ens       (:,:,:) ! node depth [m]
   real(r8), allocatable :: dz_sno_ens      (:,:,:) ! interface depth [m]
   real(r8), allocatable :: t_soisno_ens    (:,:,:) ! soil temperature [K]
   real(r8), allocatable :: wliq_soisno_ens (:,:,:) ! liquid water in layers [kg/m2]
   real(r8), allocatable :: wice_soisno_ens (:,:,:) ! ice lens in layers [kg/m2]
   real(r8), allocatable :: smp_ens         (:,:,:) ! soil matrix potential [mm]
   real(r8), allocatable :: hk_ens          (:,:,:) ! hydraulic conductivity [mm h2o/s]
   real(r8), allocatable :: t_grnd_ens        (:,:) ! ground surface temperature [K]
   real(r8), allocatable :: tleaf_ens         (:,:) ! leaf temperature [K]
   real(r8), allocatable :: ldew_ens          (:,:) ! depth of water on foliage [mm]
   real(r8), allocatable :: ldew_rain_ens     (:,:) ! depth of rain on foliage [mm]
   real(r8), allocatable :: ldew_snow_ens     (:,:) ! depth of rain on foliage [mm]
   real(r8), allocatable :: fwet_snow_ens     (:,:) ! vegetation snow fractional cover [-]
   real(r8), allocatable :: sag_ens           (:,:) ! non dimensional snow age [-]
   real(r8), allocatable :: scv_ens           (:,:) ! snow cover, water equivalent [mm]
   real(r8), allocatable :: snowdp_ens        (:,:) ! snow depth [meter]
   real(r8), allocatable :: fveg_ens          (:,:) ! fraction of vegetation cover
   real(r8), allocatable :: fsno_ens          (:,:) ! fraction of snow cover on ground
   real(r8), allocatable :: sigf_ens          (:,:) ! fraction of veg cover, excluding snow-covered veg [-]
   real(r8), allocatable :: green_ens         (:,:) ! leaf greenness
   real(r8), allocatable :: tlai_ens          (:,:) ! leaf area index
   real(r8), allocatable :: lai_ens           (:,:) ! leaf area index
   real(r8), allocatable :: tsai_ens          (:,:) ! stem area index
   real(r8), allocatable :: sai_ens           (:,:) ! stem area index
   real(r8), allocatable :: alb_ens       (:,:,:,:) ! averaged albedo [-]
   real(r8), allocatable :: ssun_ens      (:,:,:,:) ! sunlit canopy absorption for solar radiation (0-1)
   real(r8), allocatable :: ssha_ens      (:,:,:,:) ! shaded canopy absorption for solar radiation (0-1)
   real(r8), allocatable :: ssoi_ens      (:,:,:,:) ! soil absorption for solar radiation (0-1)
   real(r8), allocatable :: ssno_ens      (:,:,:,:) ! snow absorption for solar radiation (0-1)
   real(r8), allocatable :: thermk_ens        (:,:) ! canopy gap fraction for tir radiation
   real(r8), allocatable :: extkb_ens         (:,:) ! (k, g(mu)/mu) direct solar extinction coefficient
   real(r8), allocatable :: extkd_ens         (:,:) ! diffuse and scattered diffuse PAR extinction coefficient
   real(r8), allocatable :: zwt_ens           (:,:) ! the depth to water table [m]
   real(r8), allocatable :: wdsrf_ens         (:,:) ! depth of surface water [mm]
   real(r8), allocatable :: wa_ens            (:,:) ! water storage in aquifer [mm]
   real(r8), allocatable :: wetwat_ens        (:,:) ! water storage in wetland [mm]
   real(r8), allocatable :: t_lake_ens      (:,:,:) ! lake layer teperature [K]
   real(r8), allocatable :: lake_icefrac_ens(:,:,:) ! lake mass fraction of lake layer that is frozen
   real(r8), allocatable :: savedtke1_ens     (:,:) ! top level eddy conductivity (W/m K)

   ! diagnostic variables for DA
   real(r8), allocatable :: h2osoi_ens      (:,:,:) ! volumetric soil water in layers [m3/m3]
   real(r8), allocatable :: t_brt_smap_ens  (:,:,:) ! brightness temperature for radiance calculation [K]
   real(r8), allocatable :: t_brt_fy3d_ens  (:,:,:) ! brightness temperature for radiance calculation [K]
   real(r8), allocatable :: t_brt_smap        (:,:) ! brightness temperature for radiance calculation [K]
   real(r8), allocatable :: t_brt_fy3d        (:,:) ! brightness temperature for radiance calculation [K]
   real(r8), allocatable :: trad_ens          (:,:) ! radiative temperature of surface [K]
   real(r8), allocatable :: tref_ens          (:,:) ! 2 m height air temperature [kelvin]
   real(r8), allocatable :: qref_ens          (:,:) ! 2 m height air specific humidity
   real(r8), allocatable :: rhref_ens         (:,:) ! 2 m height air relative humidity
   real(r8), allocatable :: ustar_ens         (:,:) ! u* in similarity theory [m/s]
   real(r8), allocatable :: qstar_ens         (:,:) ! q* in similarity theory [kg/kg]
   real(r8), allocatable :: tstar_ens         (:,:) ! t* in similarity theory [K]
   real(r8), allocatable :: fm_ens            (:,:) ! integral of profile FUNCTION for momentum
   real(r8), allocatable :: fh_ens            (:,:) ! integral of profile FUNCTION for heat
   real(r8), allocatable :: fq_ens            (:,:) ! integral of profile FUNCTION for moisture

   ! ensemble forcing variables used for ensemble DA
   real(r8), allocatable :: forc_t_ens        (:,:) ! temperature [K]
   real(r8), allocatable :: forc_frl_ens      (:,:) ! atmospheric infrared (longwave) radiation [W/m2]
   real(r8), allocatable :: forc_prc_ens      (:,:) ! convective precipitation [mm/s]
   real(r8), allocatable :: forc_prl_ens      (:,:) ! large scale precipitation [mm/s]
   real(r8), allocatable :: forc_sols_ens     (:,:) ! atm vis direct beam solar rad onto srf [W/m2]
   real(r8), allocatable :: forc_soll_ens     (:,:) ! atm nir direct beam solar rad onto srf [W/m2]
   real(r8), allocatable :: forc_solsd_ens    (:,:) ! atm vis diffuse solar rad onto srf [W/m2]
   real(r8), allocatable :: forc_solld_ens    (:,:) ! atm nir diffuse solar rad onto srf [W/m2]

   ! save for analysis increment
   real(r8), allocatable :: t_soisno_a        (:,:) !
   real(r8), allocatable :: wliq_soisno_a     (:,:) !
   real(r8), allocatable :: wice_soisno_a     (:,:) !
   real(r8), allocatable :: t_grnd_a            (:) !
   real(r8), allocatable :: tleaf_a             (:) !
   real(r8), allocatable :: snowdp_a            (:) !
   real(r8), allocatable :: h2osoi_a          (:,:) !
   real(r8), allocatable :: t_brt_smap_a      (:,:) !
   real(r8), allocatable :: t_brt_fy3d_a      (:,:) !
   real(r8), allocatable :: trad_a              (:) !
   real(r8), allocatable :: tref_a              (:) !
   real(r8), allocatable :: qref_a              (:) !
   real(r8), allocatable :: ustar_a             (:) !
   real(r8), allocatable :: qstar_a             (:) !
   real(r8), allocatable :: tstar_a             (:) !
   real(r8), allocatable :: fm_a                (:) !
   real(r8), allocatable :: fh_a                (:) !
   real(r8), allocatable :: fq_a                (:) !

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE allocate_DATimeVariables()

!-----------------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Vars_Global
      USE MOD_SPMD_Task
      USE MOD_LandPatch, only: numpatch
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (z_sno_noda      (maxsnl+1:0,      numpatch)); z_sno_noda       (:,:) = spval
            allocate (dz_sno_noda     (maxsnl+1:0,      numpatch)); dz_sno_noda      (:,:) = spval
            allocate (t_soisno_noda   (maxsnl+1:nl_soil,numpatch)); t_soisno_noda    (:,:) = spval
            allocate (wliq_soisno_noda(maxsnl+1:nl_soil,numpatch)); wliq_soisno_noda (:,:) = spval
            allocate (wice_soisno_noda(maxsnl+1:nl_soil,numpatch)); wice_soisno_noda (:,:) = spval
            allocate (smp_noda               (1:nl_soil,numpatch)); smp_noda         (:,:) = spval
            allocate (hk_noda                (1:nl_soil,numpatch)); hk_noda          (:,:) = spval
            allocate (t_grnd_noda                      (numpatch)); t_grnd_noda        (:) = spval
            allocate (tleaf_noda                       (numpatch)); tleaf_noda         (:) = spval
            allocate (ldew_noda                        (numpatch)); ldew_noda          (:) = spval
            allocate (ldew_rain_noda                   (numpatch)); ldew_rain_noda     (:) = spval
            allocate (ldew_snow_noda                   (numpatch)); ldew_snow_noda     (:) = spval
            allocate (fwet_snow_noda                   (numpatch)); fwet_snow_noda     (:) = spval
            allocate (sag_noda                         (numpatch)); sag_noda           (:) = spval
            allocate (scv_noda                         (numpatch)); scv_noda           (:) = spval
            allocate (snowdp_noda                      (numpatch)); snowdp_noda        (:) = spval
            allocate (fveg_noda                        (numpatch)); fveg_noda          (:) = spval
            allocate (fsno_noda                        (numpatch)); fsno_noda          (:) = spval
            allocate (sigf_noda                        (numpatch)); sigf_noda          (:) = spval
            allocate (green_noda                       (numpatch)); green_noda         (:) = spval
            allocate (tlai_noda                        (numpatch)); tlai_noda          (:) = spval
            allocate (lai_noda                         (numpatch)); lai_noda           (:) = spval
            allocate (tsai_noda                        (numpatch)); tsai_noda          (:) = spval
            allocate (sai_noda                         (numpatch)); sai_noda           (:) = spval
            allocate (alb_noda                     (2,2,numpatch)); alb_noda       (:,:,:) = spval
            allocate (ssun_noda                    (2,2,numpatch)); ssun_noda      (:,:,:) = spval
            allocate (ssha_noda                    (2,2,numpatch)); ssha_noda      (:,:,:) = spval
            allocate (ssoi_noda                    (2,2,numpatch)); ssoi_noda      (:,:,:) = spval
            allocate (ssno_noda                    (2,2,numpatch)); ssno_noda      (:,:,:) = spval
            allocate (thermk_noda                      (numpatch)); thermk_noda        (:) = spval
            allocate (extkb_noda                       (numpatch)); extkb_noda         (:) = spval
            allocate (extkd_noda                       (numpatch)); extkd_noda         (:) = spval
            allocate (zwt_noda                         (numpatch)); zwt_noda           (:) = spval
            allocate (wdsrf_noda                       (numpatch)); wdsrf_noda         (:) = spval
            allocate (wa_noda                          (numpatch)); wa_noda            (:) = spval
            allocate (wetwat_noda                      (numpatch)); wetwat_noda        (:) = spval
            allocate (t_lake_noda              (nl_lake,numpatch)); t_lake_noda      (:,:) = spval
            allocate (lake_icefrac_noda        (nl_lake,numpatch)); lake_icefrac_noda(:,:) = spval
            allocate (savedtke1_noda                   (numpatch)); savedtke1_noda     (:) = spval

            ! diagnostic variables for DA
            allocate (h2osoi_noda            (1:nl_soil,numpatch)); h2osoi_noda      (:,:) = spval
            allocate (tref_noda                        (numpatch)); tref_noda          (:) = spval

            ! allocate all time-varying state variables
            allocate (z_sno_ens      (maxsnl+1:0,      DEF_DA_ENS_NUM,numpatch)); z_sno_ens       (:,:,:) = spval
            allocate (dz_sno_ens     (maxsnl+1:0,      DEF_DA_ENS_NUM,numpatch)); dz_sno_ens      (:,:,:) = spval
            allocate (t_soisno_ens   (maxsnl+1:nl_soil,DEF_DA_ENS_NUM,numpatch)); t_soisno_ens    (:,:,:) = spval
            allocate (wliq_soisno_ens(maxsnl+1:nl_soil,DEF_DA_ENS_NUM,numpatch)); wliq_soisno_ens (:,:,:) = spval
            allocate (wice_soisno_ens(maxsnl+1:nl_soil,DEF_DA_ENS_NUM,numpatch)); wice_soisno_ens (:,:,:) = spval
            allocate (smp_ens               (1:nl_soil,DEF_DA_ENS_NUM,numpatch)); smp_ens         (:,:,:) = spval
            allocate (hk_ens                (1:nl_soil,DEF_DA_ENS_NUM,numpatch)); hk_ens          (:,:,:) = spval
            allocate (t_grnd_ens                      (DEF_DA_ENS_NUM,numpatch)); t_grnd_ens        (:,:) = spval
            allocate (tleaf_ens                       (DEF_DA_ENS_NUM,numpatch)); tleaf_ens         (:,:) = spval
            allocate (ldew_ens                        (DEF_DA_ENS_NUM,numpatch)); ldew_ens          (:,:) = spval
            allocate (ldew_rain_ens                   (DEF_DA_ENS_NUM,numpatch)); ldew_rain_ens     (:,:) = spval
            allocate (ldew_snow_ens                   (DEF_DA_ENS_NUM,numpatch)); ldew_snow_ens     (:,:) = spval
            allocate (fwet_snow_ens                   (DEF_DA_ENS_NUM,numpatch)); fwet_snow_ens     (:,:) = spval
            allocate (sag_ens                         (DEF_DA_ENS_NUM,numpatch)); sag_ens           (:,:) = spval
            allocate (scv_ens                         (DEF_DA_ENS_NUM,numpatch)); scv_ens           (:,:) = spval
            allocate (snowdp_ens                      (DEF_DA_ENS_NUM,numpatch)); snowdp_ens        (:,:) = spval
            allocate (fveg_ens                        (DEF_DA_ENS_NUM,numpatch)); fveg_ens          (:,:) = spval
            allocate (fsno_ens                        (DEF_DA_ENS_NUM,numpatch)); fsno_ens          (:,:) = spval
            allocate (sigf_ens                        (DEF_DA_ENS_NUM,numpatch)); sigf_ens          (:,:) = spval
            allocate (green_ens                       (DEF_DA_ENS_NUM,numpatch)); green_ens         (:,:) = spval
            allocate (tlai_ens                        (DEF_DA_ENS_NUM,numpatch)); tlai_ens          (:,:) = spval
            allocate (lai_ens                         (DEF_DA_ENS_NUM,numpatch)); lai_ens           (:,:) = spval
            allocate (tsai_ens                        (DEF_DA_ENS_NUM,numpatch)); tsai_ens          (:,:) = spval
            allocate (sai_ens                         (DEF_DA_ENS_NUM,numpatch)); sai_ens           (:,:) = spval
            allocate (alb_ens                     (2,2,DEF_DA_ENS_NUM,numpatch)); alb_ens       (:,:,:,:) = spval
            allocate (ssun_ens                    (2,2,DEF_DA_ENS_NUM,numpatch)); ssun_ens      (:,:,:,:) = spval
            allocate (ssha_ens                    (2,2,DEF_DA_ENS_NUM,numpatch)); ssha_ens      (:,:,:,:) = spval
            allocate (ssoi_ens                    (2,2,DEF_DA_ENS_NUM,numpatch)); ssoi_ens      (:,:,:,:) = spval
            allocate (ssno_ens                    (2,2,DEF_DA_ENS_NUM,numpatch)); ssno_ens      (:,:,:,:) = spval
            allocate (thermk_ens                      (DEF_DA_ENS_NUM,numpatch)); thermk_ens        (:,:) = spval
            allocate (extkb_ens                       (DEF_DA_ENS_NUM,numpatch)); extkb_ens         (:,:) = spval
            allocate (extkd_ens                       (DEF_DA_ENS_NUM,numpatch)); extkd_ens         (:,:) = spval
            allocate (zwt_ens                         (DEF_DA_ENS_NUM,numpatch)); zwt_ens           (:,:) = spval
            allocate (wdsrf_ens                       (DEF_DA_ENS_NUM,numpatch)); wdsrf_ens         (:,:) = spval
            allocate (wa_ens                          (DEF_DA_ENS_NUM,numpatch)); wa_ens            (:,:) = spval
            allocate (wetwat_ens                      (DEF_DA_ENS_NUM,numpatch)); wetwat_ens        (:,:) = spval
            allocate (t_lake_ens              (nl_lake,DEF_DA_ENS_NUM,numpatch)); t_lake_ens      (:,:,:) = spval
            allocate (lake_icefrac_ens        (nl_lake,DEF_DA_ENS_NUM,numpatch)); lake_icefrac_ens(:,:,:) = spval
            allocate (savedtke1_ens                   (DEF_DA_ENS_NUM,numpatch)); savedtke1_ens     (:,:) = spval

            ! diagnostic variables for DA
            allocate (h2osoi_ens            (1:nl_soil,DEF_DA_ENS_NUM,numpatch)); h2osoi_ens      (:,:,:) = spval
            allocate (t_brt_smap_ens                (2,DEF_DA_ENS_NUM,numpatch)); t_brt_smap_ens  (:,:,:) = spval
            allocate (t_brt_fy3d_ens                (2,DEF_DA_ENS_NUM,numpatch)); t_brt_fy3d_ens  (:,:,:) = spval
            allocate (t_brt_smap                                   (2,numpatch)); t_brt_smap        (:,:) = spval
            allocate (t_brt_fy3d                                   (2,numpatch)); t_brt_fy3d        (:,:) = spval
            allocate (trad_ens                        (DEF_DA_ENS_NUM,numpatch)); trad_ens          (:,:) = spval
            allocate (tref_ens                        (DEF_DA_ENS_NUM,numpatch)); tref_ens          (:,:) = spval
            allocate (qref_ens                        (DEF_DA_ENS_NUM,numpatch)); qref_ens          (:,:) = spval
            allocate (rhref_ens                       (DEF_DA_ENS_NUM,numpatch)); rhref_ens         (:,:) = spval
            allocate (ustar_ens                       (DEF_DA_ENS_NUM,numpatch)); ustar_ens         (:,:) = spval
            allocate (qstar_ens                       (DEF_DA_ENS_NUM,numpatch)); qstar_ens         (:,:) = spval
            allocate (tstar_ens                       (DEF_DA_ENS_NUM,numpatch)); tstar_ens         (:,:) = spval
            allocate (fm_ens                          (DEF_DA_ENS_NUM,numpatch)); fm_ens            (:,:) = spval
            allocate (fh_ens                          (DEF_DA_ENS_NUM,numpatch)); fh_ens            (:,:) = spval
            allocate (fq_ens                          (DEF_DA_ENS_NUM,numpatch)); fq_ens            (:,:) = spval

            ! ensemble forcing variables used for ensemble DA
            allocate (forc_t_ens                      (DEF_DA_ENS_NUM,numpatch)); forc_t_ens        (:,:) = spval
            allocate (forc_frl_ens                    (DEF_DA_ENS_NUM,numpatch)); forc_frl_ens      (:,:) = spval
            allocate (forc_prc_ens                    (DEF_DA_ENS_NUM,numpatch)); forc_prc_ens      (:,:) = spval
            allocate (forc_prl_ens                    (DEF_DA_ENS_NUM,numpatch)); forc_prl_ens      (:,:) = spval
            allocate (forc_sols_ens                   (DEF_DA_ENS_NUM,numpatch)); forc_sols_ens     (:,:) = spval
            allocate (forc_soll_ens                   (DEF_DA_ENS_NUM,numpatch)); forc_soll_ens     (:,:) = spval
            allocate (forc_solsd_ens                  (DEF_DA_ENS_NUM,numpatch)); forc_solsd_ens    (:,:) = spval
            allocate (forc_solld_ens                  (DEF_DA_ENS_NUM,numpatch)); forc_solld_ens    (:,:) = spval

            ! allocate variables for analysis increment
            allocate (t_soisno_a                (maxsnl+1:nl_soil,numpatch)); t_soisno_a        (:,:) = spval
            allocate (wliq_soisno_a             (maxsnl+1:nl_soil,numpatch)); wliq_soisno_a     (:,:) = spval
            allocate (wice_soisno_a             (maxsnl+1:nl_soil,numpatch)); wice_soisno_a     (:,:) = spval
            allocate (t_grnd_a                                   (numpatch)); t_grnd_a            (:) = spval
            allocate (tleaf_a                                    (numpatch)); tleaf_a             (:) = spval
            allocate (snowdp_a                                   (numpatch)); snowdp_a            (:) = spval
            allocate (h2osoi_a                         (1:nl_soil,numpatch)); h2osoi_a          (:,:) = spval
            allocate (t_brt_smap_a                             (2,numpatch)); t_brt_smap_a      (:,:) = spval
            allocate (t_brt_fy3d_a                             (2,numpatch)); t_brt_fy3d_a      (:,:) = spval

            allocate (trad_a                                     (numpatch)); trad_a              (:) = spval
            allocate (tref_a                                     (numpatch)); tref_a              (:) = spval
            allocate (qref_a                                     (numpatch)); qref_a              (:) = spval
            allocate (ustar_a                                    (numpatch)); ustar_a             (:) = spval
            allocate (qstar_a                                    (numpatch)); qstar_a             (:) = spval
            allocate (tstar_a                                    (numpatch)); tstar_a             (:) = spval
            allocate (fm_a                                       (numpatch)); fm_a                (:) = spval
            allocate (fh_a                                       (numpatch)); fh_a                (:) = spval
            allocate (fq_a                                       (numpatch)); fq_a                (:) = spval
         ENDIF
      ENDIF

   END SUBROUTINE allocate_DATimeVariables

!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_DATimeVariables()

!-----------------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_LandPatch, only: numpatch
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            deallocate (z_sno_noda                 )
            deallocate (dz_sno_noda                )
            deallocate (t_soisno_noda              )
            deallocate (wliq_soisno_noda           )
            deallocate (wice_soisno_noda           )
            deallocate (smp_noda                   )
            deallocate (hk_noda                    )
            deallocate (t_grnd_noda                )
            deallocate (tleaf_noda                 )
            deallocate (ldew_noda                  )
            deallocate (ldew_rain_noda             )
            deallocate (ldew_snow_noda             )
            deallocate (fwet_snow_noda             )
            deallocate (sag_noda                   )
            deallocate (scv_noda                   )
            deallocate (snowdp_noda                )
            deallocate (fveg_noda                  )
            deallocate (fsno_noda                  )
            deallocate (sigf_noda                  )
            deallocate (green_noda                 )
            deallocate (tlai_noda                  )
            deallocate (lai_noda                   )
            deallocate (tsai_noda                  )
            deallocate (sai_noda                   )
            deallocate (alb_noda                   )
            deallocate (ssun_noda                  )
            deallocate (ssha_noda                  )
            deallocate (ssoi_noda                  )
            deallocate (ssno_noda                  )
            deallocate (thermk_noda                )
            deallocate (extkb_noda                 )
            deallocate (extkd_noda                 )
            deallocate (zwt_noda                   )
            deallocate (wdsrf_noda                 )
            deallocate (wa_noda                    )
            deallocate (wetwat_noda                )
            deallocate (t_lake_noda                )
            deallocate (lake_icefrac_noda          )
            deallocate (savedtke1_noda             )

            deallocate (tref_noda                  )
            deallocate (h2osoi_noda                )

            deallocate (z_sno_ens                  )
            deallocate (dz_sno_ens                 )
            deallocate (t_soisno_ens               )
            deallocate (wliq_soisno_ens            )
            deallocate (wice_soisno_ens            )
            deallocate (smp_ens                    )
            deallocate (hk_ens                     )
            deallocate (t_grnd_ens                 )
            deallocate (tleaf_ens                  )
            deallocate (ldew_ens                   )
            deallocate (ldew_rain_ens              )
            deallocate (ldew_snow_ens              )
            deallocate (fwet_snow_ens              )
            deallocate (sag_ens                    )
            deallocate (scv_ens                    )
            deallocate (snowdp_ens                 )
            deallocate (fveg_ens                   )
            deallocate (fsno_ens                   )
            deallocate (sigf_ens                   )
            deallocate (green_ens                  )
            deallocate (tlai_ens                   )
            deallocate (lai_ens                    )
            deallocate (tsai_ens                   )
            deallocate (sai_ens                    )
            deallocate (alb_ens                    )
            deallocate (ssun_ens                   )
            deallocate (ssha_ens                   )
            deallocate (ssoi_ens                   )
            deallocate (ssno_ens                   )
            deallocate (thermk_ens                 )
            deallocate (extkb_ens                  )
            deallocate (extkd_ens                  )
            deallocate (zwt_ens                    )
            deallocate (wdsrf_ens                  )
            deallocate (wa_ens                     )
            deallocate (wetwat_ens                 )
            deallocate (t_lake_ens                 )
            deallocate (lake_icefrac_ens           )
            deallocate (savedtke1_ens              )

            deallocate (h2osoi_ens                 )
            deallocate (t_brt_smap_ens             )
            deallocate (t_brt_fy3d_ens             )
            deallocate (t_brt_smap                 )
            deallocate (t_brt_fy3d                 )
            deallocate (trad_ens                   )
            deallocate (tref_ens                   )
            deallocate (qref_ens                   )
            deallocate (rhref_ens                  )
            deallocate (ustar_ens                  )
            deallocate (qstar_ens                  )
            deallocate (tstar_ens                  )
            deallocate (fm_ens                     )
            deallocate (fh_ens                     )
            deallocate (fq_ens                     )

            deallocate (forc_t_ens                 )
            deallocate (forc_frl_ens               )
            deallocate (forc_prc_ens               )
            deallocate (forc_prl_ens               )
            deallocate (forc_sols_ens              )
            deallocate (forc_soll_ens              )
            deallocate (forc_solsd_ens             )
            deallocate (forc_solld_ens             )

            deallocate (t_soisno_a                 )
            deallocate (wliq_soisno_a              )
            deallocate (wice_soisno_a              )
            deallocate (t_grnd_a                   )
            deallocate (tleaf_a                    )
            deallocate (snowdp_a                   )
            deallocate (h2osoi_a                   )
            deallocate (t_brt_smap_a               )
            deallocate (t_brt_fy3d_a               )
            deallocate (trad_a                     )
            deallocate (tref_a                     )
            deallocate (qref_a                     )
            deallocate (ustar_a                    )
            deallocate (qstar_a                    )
            deallocate (tstar_a                    )
            deallocate (fm_a                       )
            deallocate (fh_a                       )
            deallocate (fq_a                       )
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_DATimeVariables

!-----------------------------------------------------------------------------

   SUBROUTINE WRITE_DATimeVariables (idate, lc_year, site, dir_restart)

!-----------------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_Namelist, only : DEF_REST_CompressLevel, DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS, &
         DEF_USE_IRRIGATION, DEF_USE_Dynamic_Lake
      USE MOD_LandPatch
      USE MOD_NetCDFVector
      USE MOD_Vars_Global
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer, intent(in) :: idate(3)
      integer, intent(in) :: lc_year      ! year of land cover type data
      character(len=*), intent(in) :: site
      character(len=*), intent(in) :: dir_restart

!------------------------ Local Variables ------------------------------------
      character(len=256) :: file_restart
      character(len=14)  :: cdate
      character(len=256) :: cyear         ! character for lc_year
      integer :: compress
      integer :: i

!-----------------------------------------------------------------------------
      compress = DEF_REST_CompressLevel

      ! land cover type year
      write(cyear,'(i4.4)') lc_year
      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)

      IF (p_is_master) THEN
         CALL system('mkdir -p ' // trim(dir_restart)//'/'//trim(cdate))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      file_restart = trim(dir_restart)// '/'//trim(cdate)//'/' // trim(site) //'_restart_DA_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

      CALL ncio_create_file_vector      (file_restart, landpatch)

      CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'snow',     -maxsnl       )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'snowp1',   -maxsnl+1     )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soilsnow', nl_soil-maxsnl)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil',     nl_soil)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'lake',     nl_lake)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'band', 2)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'rtyp', 2)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'ens', DEF_DA_ENS_NUM)

      ! Time-varying state variables which reaquired by restart run
      CALL ncio_write_vector (file_restart, 'z_sno   '   , 'snow',     -maxsnl,        'ens', DEF_DA_ENS_NUM, 'patch', landpatch, z_sno_ens,       compress) ! node depth [m]
      CALL ncio_write_vector (file_restart, 'dz_sno  '   , 'snow',     -maxsnl,        'ens', DEF_DA_ENS_NUM, 'patch', landpatch, dz_sno_ens,      compress) ! interface depth [m]
      CALL ncio_write_vector (file_restart, 't_soisno'   , 'soilsnow', nl_soil-maxsnl, 'ens', DEF_DA_ENS_NUM, 'patch', landpatch, t_soisno_ens,    compress) ! soil temperature [K]
      CALL ncio_write_vector (file_restart, 'wliq_soisno', 'soilsnow', nl_soil-maxsnl, 'ens', DEF_DA_ENS_NUM, 'patch', landpatch, wliq_soisno_ens, compress) ! liquid water in layers [kg/m2]
      CALL ncio_write_vector (file_restart, 'wice_soisno', 'soilsnow', nl_soil-maxsnl, 'ens', DEF_DA_ENS_NUM, 'patch', landpatch, wice_soisno_ens, compress) ! ice lens in layers [kg/m2]
      CALL ncio_write_vector (file_restart, 'smp',         'soil',     nl_soil,        'ens', DEF_DA_ENS_NUM, 'patch', landpatch, smp_ens,         compress) ! soil matrix potential [mm]
      CALL ncio_write_vector (file_restart, 'hk',          'soil',     nl_soil,        'ens', DEF_DA_ENS_NUM, 'patch', landpatch, hk_ens,          compress) ! hydraulic conductivity [mm h2o/s]
      CALL ncio_write_vector (file_restart, 't_grnd',      'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, t_grnd_ens,    compress) ! ground surface temperature [K]
      CALL ncio_write_vector (file_restart, 'tleaf',       'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, tleaf_ens,     compress) ! leaf temperature [K]
      CALL ncio_write_vector (file_restart, 'ldew',        'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, ldew_ens,      compress) ! depth of water on foliage [mm]
      CALL ncio_write_vector (file_restart, 'ldew_rain',   'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, ldew_rain_ens, compress) ! depth of water on foliage [mm]
      CALL ncio_write_vector (file_restart, 'ldew_snow',   'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, ldew_snow_ens, compress) ! depth of water on foliage [mm]
      CALL ncio_write_vector (file_restart, 'fwet_snow',   'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, fwet_snow_ens, compress) ! vegetation snow fractional cover [-]
      CALL ncio_write_vector (file_restart, 'sag',         'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, sag_ens,       compress) ! non dimensional snow age [-]
      CALL ncio_write_vector (file_restart, 'scv',         'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, scv_ens,       compress) ! snow cover, water equivalent [mm]
      CALL ncio_write_vector (file_restart, 'snowdp',      'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, snowdp_ens,    compress) ! snow depth [meter]
      CALL ncio_write_vector (file_restart, 'fveg',        'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, fveg_ens,      compress) ! fraction of vegetation cover
      CALL ncio_write_vector (file_restart, 'fsno',        'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, fsno_ens,      compress) ! fraction of snow cover on ground
      CALL ncio_write_vector (file_restart, 'sigf',        'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, sigf_ens,      compress) ! fraction of veg cover, excluding snow-covered veg [-]
      CALL ncio_write_vector (file_restart, 'green',       'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, green_ens,     compress) ! leaf greenness
      CALL ncio_write_vector (file_restart, 'tlai',        'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, tlai_ens,      compress) ! leaf area index
      CALL ncio_write_vector (file_restart, 'lai',         'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, lai_ens,       compress) ! leaf area index
      CALL ncio_write_vector (file_restart, 'tsai',        'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, tsai_ens,      compress) ! stem area index
      CALL ncio_write_vector (file_restart, 'sai',         'ens',      DEF_DA_ENS_NUM, 'patch', landpatch, sai_ens,       compress) ! stem area index
      CALL ncio_write_vector (file_restart, 'alb',         'band',     2, 'rtyp', 2, 'ens', DEF_DA_ENS_NUM, 'patch', landpatch, alb_ens,  compress) ! averaged albedo [-]
      CALL ncio_write_vector (file_restart, 'ssun',        'band',     2, 'rtyp', 2, 'ens', DEF_DA_ENS_NUM, 'patch', landpatch, ssun_ens, compress) ! sunlit canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'ssha',        'band',     2, 'rtyp', 2, 'ens', DEF_DA_ENS_NUM, 'patch', landpatch, ssha_ens, compress) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'ssoi',        'band',     2, 'rtyp', 2, 'ens', DEF_DA_ENS_NUM, 'patch', landpatch, ssoi_ens, compress) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'ssno',        'band',     2, 'rtyp', 2, 'ens', DEF_DA_ENS_NUM, 'patch', landpatch, ssno_ens, compress) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'thermk',      'ens',      DEF_DA_ENS_NUM,     'patch', landpatch, thermk_ens,    compress) ! canopy gap fraction for tir radiation
      CALL ncio_write_vector (file_restart, 'extkb',       'ens',      DEF_DA_ENS_NUM,     'patch', landpatch, extkb_ens,     compress) ! (k, g(mu)/mu) direct solar extinction coefficient
      CALL ncio_write_vector (file_restart, 'extkd',       'ens',      DEF_DA_ENS_NUM,     'patch', landpatch, extkd_ens,     compress) ! diffuse and scattered diffuse PAR extinction coefficient
      CALL ncio_write_vector (file_restart, 'zwt',         'ens',      DEF_DA_ENS_NUM,     'patch', landpatch, zwt_ens,       compress) ! the depth to water table [m]
      CALL ncio_write_vector (file_restart, 'wdsrf',       'ens',      DEF_DA_ENS_NUM,     'patch', landpatch, wdsrf_ens,     compress) ! depth of surface water [mm]
      CALL ncio_write_vector (file_restart, 'wa',          'ens',      DEF_DA_ENS_NUM,     'patch', landpatch, wa_ens,        compress) ! water storage in aquifer [mm]
      CALL ncio_write_vector (file_restart, 'wetwat',      'ens',      DEF_DA_ENS_NUM,     'patch', landpatch, wetwat_ens,    compress) ! water storage in wetland [mm]
      CALL ncio_write_vector (file_restart, 't_lake',      'lake',     nl_lake,        'ens',   DEF_DA_ENS_NUM, 'patch', landpatch, t_lake_ens,       compress)
      CALL ncio_write_vector (file_restart, 'lake_icefrc', 'lake',     nl_lake,        'ens',   DEF_DA_ENS_NUM, 'patch', landpatch, lake_icefrac_ens, compress)
      CALL ncio_write_vector (file_restart, 'savedtke1',   'ens',      DEF_DA_ENS_NUM,     'patch', landpatch, savedtke1_ens, compress)

   END SUBROUTINE WRITE_DATimeVariables

!-----------------------------------------------------------------------------

   SUBROUTINE READ_DATimeVariables (idate, lc_year, site, dir_restart)

!-----------------------------------------------------------------------------
      USE MOD_Namelist
      USE MOD_SPMD_Task
      USE MOD_NetCDFVector
#ifdef RangeCheck
      USE MOD_RangeCheck
#endif
      USE MOD_LandPatch
      USE MOD_Vars_Global
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer, intent(in) :: idate(3)
      integer, intent(in) :: lc_year      !year of land cover type data
      character(len=*), intent(in) :: site
      character(len=*), intent(in) :: dir_restart

!------------------------ Local Variables ------------------------------------
      character(len=256) :: file_restart
      character(len=14)  :: cdate, cyear

!-----------------------------------------------------------------------------
#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         write (*, *) 'Loading DA Time Variables ...'
      END IF

      ! land cover type year
      write (cyear, '(i4.4)') lc_year

      write (cdate, '(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
      file_restart = trim(dir_restart)//'/'//trim(cdate)//'/'//trim(site)//'_restart_DA_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

      ! Time-varying state variables which reaquired by restart run
      CALL ncio_read_vector(file_restart, 'z_sno   ', -maxsnl, DEF_DA_ENS_NUM, landpatch, z_sno_ens)             ! node depth [m]
      CALL ncio_read_vector(file_restart, 'dz_sno  ', -maxsnl, DEF_DA_ENS_NUM, landpatch, dz_sno_ens)            ! interface depth [m]
      CALL ncio_read_vector(file_restart, 't_soisno', nl_soil - maxsnl, DEF_DA_ENS_NUM, landpatch, t_soisno_ens)   ! soil temperature [K]
      CALL ncio_read_vector(file_restart, 'wliq_soisno', nl_soil - maxsnl, DEF_DA_ENS_NUM, landpatch, wliq_soisno_ens)! liquid water in layers [kg/m2]
      CALL ncio_read_vector(file_restart, 'wice_soisno', nl_soil - maxsnl, DEF_DA_ENS_NUM, landpatch, wice_soisno_ens)! ice lens in layers [kg/m2]
      CALL ncio_read_vector(file_restart, 'smp', nl_soil, DEF_DA_ENS_NUM, landpatch, smp_ens)        ! soil matrix potential [mm]
      CALL ncio_read_vector(file_restart, 'hk', nl_soil, DEF_DA_ENS_NUM, landpatch, hk_ens)         ! hydraulic conductivity [mm h2o/s]
      CALL ncio_read_vector(file_restart, 't_grnd  ', DEF_DA_ENS_NUM, landpatch, t_grnd_ens) ! ground surface temperature [K]
      CALL ncio_read_vector(file_restart, 'tleaf   ', DEF_DA_ENS_NUM, landpatch, tleaf_ens) ! leaf temperature [K]
      CALL ncio_read_vector(file_restart, 'ldew    ', DEF_DA_ENS_NUM, landpatch, ldew_ens) ! depth of water on foliage [mm]
      CALL ncio_read_vector(file_restart, 'ldew_rain', DEF_DA_ENS_NUM, landpatch, ldew_rain_ens) ! depth of water on foliage [mm]
      CALL ncio_read_vector(file_restart, 'ldew_snow', DEF_DA_ENS_NUM, landpatch, ldew_snow_ens) ! depth of water on foliage [mm]
      CALL ncio_read_vector(file_restart, 'fwet_snow', DEF_DA_ENS_NUM, landpatch, fwet_snow_ens) ! vegetation snow fractional cover [-]
      CALL ncio_read_vector(file_restart, 'sag     ', DEF_DA_ENS_NUM, landpatch, sag_ens) ! non dimensional snow age [-]
      CALL ncio_read_vector(file_restart, 'scv     ', DEF_DA_ENS_NUM, landpatch, scv_ens) ! snow cover, water equivalent [mm]
      CALL ncio_read_vector(file_restart, 'snowdp  ', DEF_DA_ENS_NUM, landpatch, snowdp_ens) ! snow depth [meter]
      CALL ncio_read_vector(file_restart, 'fveg    ', DEF_DA_ENS_NUM, landpatch, fveg_ens) ! fraction of vegetation cover
      CALL ncio_read_vector(file_restart, 'fsno    ', DEF_DA_ENS_NUM, landpatch, fsno_ens) ! fraction of snow cover on ground
      CALL ncio_read_vector(file_restart, 'sigf    ', DEF_DA_ENS_NUM, landpatch, sigf_ens) ! fraction of veg cover, excluding snow-covered veg [-]
      CALL ncio_read_vector(file_restart, 'green   ', DEF_DA_ENS_NUM, landpatch, green_ens) ! leaf greenness
      CALL ncio_read_vector(file_restart, 'lai     ', DEF_DA_ENS_NUM, landpatch, lai_ens) ! leaf area index
      CALL ncio_read_vector(file_restart, 'tlai    ', DEF_DA_ENS_NUM, landpatch, tlai_ens) ! leaf area index
      CALL ncio_read_vector(file_restart, 'sai     ', DEF_DA_ENS_NUM, landpatch, sai_ens) ! stem area index
      CALL ncio_read_vector(file_restart, 'tsai    ', DEF_DA_ENS_NUM, landpatch, tsai_ens) ! stem area index
      CALL ncio_read_vector(file_restart, 'alb     ', 2, 2, DEF_DA_ENS_NUM, landpatch, alb_ens) ! averaged albedo [-]
      CALL ncio_read_vector(file_restart, 'ssun    ', 2, 2, DEF_DA_ENS_NUM, landpatch, ssun_ens) ! sunlit canopy absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'ssha    ', 2, 2, DEF_DA_ENS_NUM, landpatch, ssha_ens) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'ssoi    ', 2, 2, DEF_DA_ENS_NUM, landpatch, ssoi_ens) ! soil absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'ssno    ', 2, 2, DEF_DA_ENS_NUM, landpatch, ssno_ens) ! snow absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'thermk  ', DEF_DA_ENS_NUM, landpatch, thermk_ens) ! canopy gap fraction for tir radiation
      CALL ncio_read_vector(file_restart, 'extkb   ', DEF_DA_ENS_NUM, landpatch, extkb_ens) ! (k, g(mu)/mu) direct solar extinction coefficient
      CALL ncio_read_vector(file_restart, 'extkd   ', DEF_DA_ENS_NUM, landpatch, extkd_ens) ! diffuse and scattered diffuse PAR extinction coefficient
      CALL ncio_read_vector(file_restart, 'zwt     ', DEF_DA_ENS_NUM, landpatch, zwt_ens) ! the depth to water table [m]
      CALL ncio_read_vector(file_restart, 'wdsrf   ', DEF_DA_ENS_NUM, landpatch, wdsrf_ens) ! depth of surface water [mm]
      CALL ncio_read_vector(file_restart, 'wa      ', DEF_DA_ENS_NUM, landpatch, wa_ens) ! water storage in aquifer [mm]
      CALL ncio_read_vector(file_restart, 'wetwat  ', DEF_DA_ENS_NUM, landpatch, wetwat_ens) ! water storage in wetland [mm]
      CALL ncio_read_vector(file_restart, 't_lake  ', nl_lake, DEF_DA_ENS_NUM, landpatch, t_lake_ens) ! lake temperature [K]
      CALL ncio_read_vector(file_restart, 'lake_icefrc', nl_lake, DEF_DA_ENS_NUM, landpatch, lake_icefrac_ens) ! lake ice fraction [-]
      CALL ncio_read_vector(file_restart, 'savedtke1  ', DEF_DA_ENS_NUM, landpatch, savedtke1_ens) ! saved tke1 [m2/s2]

#ifdef RangeCheck
      CALL check_DATimeVariables
#endif

      IF (p_is_master) THEN
         write (*, *) 'Loading DA Time Variables done.'
      END IF

   END SUBROUTINE READ_DATimeVariables

#ifdef RangeCheck
!-----------------------------------------------------------------------

   SUBROUTINE check_DATimeVariables ()

!-----------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_RangeCheck
      USE MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS, DEF_USE_IRRIGATION, &
         DEF_USE_SNICAR, DEF_USE_Dynamic_Lake
      IMPLICIT NONE

#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write (*, *) 'Checking DA Time Variables ...'
      END IF

      CALL check_vector_data ('z_sno       [m]    ', z_sno_ens       ) ! node depth [m]
      CALL check_vector_data ('dz_sno      [m]    ', dz_sno_ens      ) ! interface depth [m]
      CALL check_vector_data ('t_soisno    [K]    ', t_soisno_ens    ) ! soil temperature [K]
      CALL check_vector_data ('wliq_soisno [kg/m2]', wliq_soisno_ens ) ! liquid water in layers [kg/m2]
      CALL check_vector_data ('wice_soisno [kg/m2]', wice_soisno_ens ) ! ice lens in layers [kg/m2]
      CALL check_vector_data ('smp         [mm]   ', smp_ens         ) ! soil matrix potential [mm]
      CALL check_vector_data ('hk          [mm/s] ', hk_ens          ) ! hydraulic conductivity [mm h2o/s]
      CALL check_vector_data ('t_grnd      [K]    ', t_grnd_ens      ) ! ground surface temperature [K]
      CALL check_vector_data ('tleaf       [K]    ', tleaf_ens       ) ! leaf temperature [K]
      CALL check_vector_data ('ldew        [mm]   ', ldew_ens        ) ! depth of water on foliage [mm]
      CALL check_vector_data ('ldew_rain   [mm]   ', ldew_rain_ens   ) ! depth of rain on foliage [mm]
      CALL check_vector_data ('ldew_snow   [mm]   ', ldew_snow_ens   ) ! depth of snow on foliage [mm]
      CALL check_vector_data ('fwet_snow   [mm]   ', fwet_snow_ens   ) ! vegetation snow fractional cover [-]
      CALL check_vector_data ('sag         [-]    ', sag_ens         ) ! non dimensional snow age [-]
      CALL check_vector_data ('scv         [mm]   ', scv_ens         ) ! snow cover, water equivalent [mm]
      CALL check_vector_data ('snowdp      [m]    ', snowdp_ens      ) ! snow depth [meter]
      CALL check_vector_data ('fveg        [-]    ', fveg_ens        ) ! fraction of vegetation cover
      CALL check_vector_data ('fsno        [-]    ', fsno_ens        ) ! fraction of snow cover on ground
      CALL check_vector_data ('sigf        [-]    ', sigf_ens        ) ! fraction of veg cover, excluding snow-covered veg [-]
      CALL check_vector_data ('green       [-]    ', green_ens       ) ! leaf greenness
      CALL check_vector_data ('lai         [-]    ', lai_ens         ) ! leaf area index
      CALL check_vector_data ('tlai        [-]    ', tlai_ens        ) ! leaf area index
      CALL check_vector_data ('sai         [-]    ', sai_ens         ) ! stem area index
      CALL check_vector_data ('tsai        [-]    ', tsai_ens        ) ! stem area index
      CALL check_vector_data ('alb         [-]    ', alb_ens         ) ! averaged albedo [-]
      CALL check_vector_data ('ssun        [-]    ', ssun_ens        ) ! sunlit canopy absorption for solar radiation (0-1)
      CALL check_vector_data ('ssha        [-]    ', ssha_ens        ) ! shaded canopy absorption for solar radiation (0-1)
      CALL check_vector_data ('ssoi        [-]    ', ssoi_ens        ) ! soil absorption for solar radiation (0-1)
      CALL check_vector_data ('ssno        [-]    ', ssno_ens        ) ! snow absorption for solar radiation (0-1)
      CALL check_vector_data ('thermk      [-]    ', thermk_ens      ) ! canopy gap fraction for tir radiation
      CALL check_vector_data ('extkb       [-]    ', extkb_ens       ) ! (k, g(mu)/mu) direct solar extinction coefficient
      CALL check_vector_data ('extkd       [-]    ', extkd_ens       ) ! diffuse and scattered diffuse PAR extinction coefficient
      CALL check_vector_data ('zwt         [m]    ', zwt_ens         ) ! the depth to water table [m]
      CALL check_vector_data ('wdsrf       [mm]   ', wdsrf_ens       ) ! depth of surface water [mm]
      CALL check_vector_data ('wa          [mm]   ', wa_ens          ) ! water storage in aquifer [mm]
      CALL check_vector_data ('wetwat      [mm]   ', wetwat_ens      ) ! water storage in wetland [mm]
      CALL check_vector_data ('t_lake      [K]    ', t_lake_ens      ) ! lake temperature [K]
      CALL check_vector_data ('lake_icefrc [-]    ', lake_icefrac_ens) ! lake ice fraction [-]
      CALL check_vector_data ('savedtke1   [W/m K]', savedtke1_ens   ) ! saved tke1 [m2/s2]

#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif

   END SUBROUTINE check_DATimeVariables
!-----------------------------------------------------------------------------
#endif

END MODULE MOD_DA_Vars_TimeVariables
#endif
