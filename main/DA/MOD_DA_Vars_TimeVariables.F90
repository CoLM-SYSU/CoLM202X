#include <define.h>

MODULE MOD_DA_Vars_TimeVariables
   !//TODO: Lu Li: Only support default vars and IGBP now, need
   !               to extend to support other land cover types

   USE MOD_Precision
   USE MOD_TimeManager

   IMPLICIT NONE
   SAVE
! -----------------------------------------------------------------
   ! Time-varying state variables which reaquired by restart run
   real(r8), allocatable :: z_sno_ens       (:,:,:) ! node depth [m]
   real(r8), allocatable :: dz_sno_ens      (:,:,:) ! interface depth [m]
   real(r8), allocatable :: t_soisno_ens    (:,:,:) ! soil temperature [K]
   real(r8), allocatable :: wliq_soisno_ens (:,:,:) ! liquid water in layers [kg/m2]
   real(r8), allocatable :: wice_soisno_ens (:,:,:) ! ice lens in layers [kg/m2]
   real(r8), allocatable :: h2osoi_ens      (:,:,:) ! volumetric soil water in layers [m3/m3]

   real(r8), allocatable :: smp_ens         (:,:,:) ! soil matrix potential [mm]
   real(r8), allocatable :: hk_ens          (:,:,:) ! hydraulic conductivity [mm h2o/s]
   real(r8), allocatable :: rootr_ens       (:,:,:) ! transpiration contribution fraction from different layers
   real(r8), allocatable :: rootflux_ens    (:,:,:) ! water exchange between soil and root. Positive: soil->root [?]

   real(r8), allocatable :: vegwp_ens       (:,:,:) ! vegetation water potential [mm]
   real(r8), allocatable :: gs0sun_ens        (:,:) ! working copy of sunlit stomata conductance
   real(r8), allocatable :: gs0sha_ens        (:,:) ! working copy of shalit stomata conductance

   real(r8), allocatable :: lai_old_ens       (:,:) ! lai in last time step
   real(r8), allocatable :: o3uptakesun_ens   (:,:) ! Ozone does, sunlit leaf (mmol O3/m^2)
   real(r8), allocatable :: o3uptakesha_ens   (:,:) ! Ozone does, shaded leaf (mmol O3/m^2)

   real(r8), allocatable :: rstfacsun_out_ens (:,:) ! factor of soil water stress on sunlit leaf
   real(r8), allocatable :: rstfacsha_out_ens (:,:) ! factor of soil water stress on shaded leaf
   real(r8), allocatable :: gssun_out_ens     (:,:) ! stomata conductance on sunlit leaf
   real(r8), allocatable :: gssha_out_ens     (:,:) ! stomata conductance on shaded leaf
   real(r8), allocatable :: assimsun_out_ens  (:,:) ! diagnostic sunlit leaf assim value for output
   real(r8), allocatable :: assimsha_out_ens  (:,:) ! diagnostic sunlit leaf etr value for output
   real(r8), allocatable :: etrsun_out_ens    (:,:) ! diagnostic shaded leaf assim for output
   real(r8), allocatable :: etrsha_out_ens    (:,:) ! diagnostic shaded leaf etr for output

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
   real(r8), allocatable :: tlai_ens          (:,:) ! leaf area index !TODO: do not in driver?
   real(r8), allocatable :: lai_ens           (:,:) ! leaf area index
   real(r8), allocatable :: laisun_ens        (:,:) ! leaf area index for sunlit leaf
   real(r8), allocatable :: laisha_ens        (:,:) ! leaf area index for shaded leaf
   real(r8), allocatable :: tsai_ens          (:,:) ! stem area index !TODO: do not have
   real(r8), allocatable :: sai_ens           (:,:) ! stem area index
   real(r8), allocatable :: coszen_ens        (:,:) ! cosine of solar zenith angle
   real(r8), allocatable :: alb_ens       (:,:,:,:) ! averaged albedo [-]
   real(r8), allocatable :: ssun_ens      (:,:,:,:) ! sunlit canopy absorption for solar radiation (0-1)
   real(r8), allocatable :: ssha_ens      (:,:,:,:) ! shaded canopy absorption for solar radiation (0-1)
   real(r8), allocatable :: ssoi_ens      (:,:,:,:) ! soil absorption for solar radiation (0-1)
   real(r8), allocatable :: ssno_ens      (:,:,:,:) ! snow absorption for solar radiation (0-1)
   real(r8), allocatable :: thermk_ens        (:,:) ! canopy gap fraction for tir radiation
   real(r8), allocatable :: extkb_ens         (:,:) ! (k, g(mu)/mu) direct solar extinction coefficient
   real(r8), allocatable :: extkd_ens         (:,:) ! diffuse and scattered diffuse PAR extinction coefficient
   real(r8), allocatable :: zwt_ens           (:,:) ! the depth to water table [m]
   real(r8), allocatable :: wa_ens            (:,:) ! water storage in aquifer [mm]
   real(r8), allocatable :: wetwat_ens        (:,:) ! water storage in wetland [mm]
   real(r8), allocatable :: wat_ens           (:,:) ! total water storage [mm]
   real(r8), allocatable :: wdsrf_ens         (:,:) ! depth of surface water [mm]
   real(r8), allocatable :: rss_ens           (:,:) ! soil surface resistance [s/m]

   real(r8), allocatable :: t_lake_ens      (:,:,:) ! lake layer teperature [K]
   real(r8), allocatable :: lake_icefrac_ens(:,:,:) ! lake mass fraction of lake layer that is frozen
   real(r8), allocatable :: savedtke1_ens     (:,:) ! top level eddy conductivity (W/m K)
   real(r8), allocatable :: dz_lake_ens     (:,:,:) ! lake depth [m]

   real(r8), allocatable :: snw_rds_ens     (:,:,:) ! effective grain radius (col,lyr) [microns, m-6]
   real(r8), allocatable :: mss_bcpho_ens   (:,:,:) ! mass of hydrophobic BC in snow  (col,lyr) [kg]
   real(r8), allocatable :: mss_bcphi_ens   (:,:,:) ! mass of hydrophillic BC in snow (col,lyr) [kg]
   real(r8), allocatable :: mss_ocpho_ens   (:,:,:) ! mass of hydrophobic OC in snow  (col,lyr) [kg]
   real(r8), allocatable :: mss_ocphi_ens   (:,:,:) ! mass of hydrophillic OC in snow (col,lyr) [kg]
   real(r8), allocatable :: mss_dst1_ens    (:,:,:) ! mass of dust species 1 in snow  (col,lyr) [kg]
   real(r8), allocatable :: mss_dst2_ens    (:,:,:) ! mass of dust species 2 in snow  (col,lyr) [kg]
   real(r8), allocatable :: mss_dst3_ens    (:,:,:) ! mass of dust species 3 in snow  (col,lyr) [kg]
   real(r8), allocatable :: mss_dst4_ens    (:,:,:) ! mass of dust species 4 in snow  (col,lyr) [kg]
   real(r8), allocatable :: ssno_lyr_ens(:,:,:,:,:) ! snow layer absorption [-]

   real(r8), allocatable :: trad_ens          (:,:) ! radiative temperature of surface [K]
   real(r8), allocatable :: tref_ens          (:,:) ! 2 m height air temperature [kelvin]
   real(r8), allocatable :: qref_ens          (:,:) ! 2 m height air specific humidity
   real(r8), allocatable :: rst_ens           (:,:) ! canopy stomatal resistance (s/m)
   real(r8), allocatable :: emis_ens          (:,:) ! averaged bulk surface emissivity
   real(r8), allocatable :: z0m_ens           (:,:) ! effective roughness [m]
   real(r8), allocatable :: displa_ens        (:,:) ! zero displacement height [m] !TODO: do not have
   real(r8), allocatable :: zol_ens           (:,:) ! dimensionless height (z/L) used in Monin-Obukhov theory
   real(r8), allocatable :: rib_ens           (:,:) ! bulk Richardson number in surface layer
   real(r8), allocatable :: ustar_ens         (:,:) ! u* in similarity theory [m/s]
   real(r8), allocatable :: qstar_ens         (:,:) ! q* in similarity theory [kg/kg]
   real(r8), allocatable :: tstar_ens         (:,:) ! t* in similarity theory [K]
   real(r8), allocatable :: fm_ens            (:,:) ! integral of profile FUNCTION for momentum
   real(r8), allocatable :: fh_ens            (:,:) ! integral of profile FUNCTION for heat
   real(r8), allocatable :: fq_ens            (:,:) ! integral of profile FUNCTION for moisture

   real(r8), allocatable :: brt_temp_ens    (:,:,:) ! brightness temperature for radiance calculation [K]

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_TimeVariables_ens
   PUBLIC :: deallocate_TimeVariables_ens
   PUBLIC :: READ_TimeVariables_ens
   PUBLIC :: WRITE_TimeVariables_ens
#ifdef RangeCheck
   PUBLIC :: check_TimeVariables_ens
#endif

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_TimeVariables_ens

!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Vars_Global
      USE MOD_SPMD_Task
      USE MOD_LandPatch, only: numpatch
      IMPLICIT NONE

!-----------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (z_sno_ens      (maxsnl+1:0,      num_ens,numpatch)); z_sno_ens       (:,:,:) = spval
            allocate (dz_sno_ens     (maxsnl+1:0,      num_ens,numpatch)); dz_sno_ens      (:,:,:) = spval
            allocate (t_soisno_ens   (maxsnl+1:nl_soil,num_ens,numpatch)); t_soisno_ens    (:,:,:) = spval
            allocate (wliq_soisno_ens(maxsnl+1:nl_soil,num_ens,numpatch)); wliq_soisno_ens (:,:,:) = spval
            allocate (wice_soisno_ens(maxsnl+1:nl_soil,num_ens,numpatch)); wice_soisno_ens (:,:,:) = spval
            allocate (h2osoi_ens            (1:nl_soil,num_ens,numpatch)); h2osoi_ens      (:,:,:) = spval

            allocate (smp_ens               (1:nl_soil,num_ens,numpatch)); smp_ens         (:,:,:) = spval
            allocate (hk_ens                (1:nl_soil,num_ens,numpatch)); hk_ens          (:,:,:) = spval
            allocate (rootr_ens             (1:nl_soil,num_ens,numpatch)); rootr_ens       (:,:,:) = spval
            allocate (rootflux_ens          (1:nl_soil,num_ens,numpatch)); rootflux_ens    (:,:,:) = spval

            allocate (vegwp_ens             (1:nvegwcs,num_ens,numpatch)); vegwp_ens       (:,:,:) = spval
            allocate (gs0sun_ens                      (num_ens,numpatch)); gs0sun_ens        (:,:) = spval
            allocate (gs0sha_ens                      (num_ens,numpatch)); gs0sha_ens        (:,:) = spval

            allocate (lai_old_ens                     (num_ens,numpatch)); lai_old_ens       (:,:) = spval
            allocate (o3uptakesun_ens                 (num_ens,numpatch)); o3uptakesun_ens   (:,:) = spval
            allocate (o3uptakesha_ens                 (num_ens,numpatch)); o3uptakesha_ens   (:,:) = spval

            allocate (rstfacsun_out_ens               (num_ens,numpatch)); rstfacsun_out_ens (:,:) = spval
            allocate (rstfacsha_out_ens               (num_ens,numpatch)); rstfacsha_out_ens (:,:) = spval
            allocate (gssun_out_ens                   (num_ens,numpatch)); gssun_out_ens     (:,:) = spval
            allocate (gssha_out_ens                   (num_ens,numpatch)); gssha_out_ens     (:,:) = spval
            allocate (assimsun_out_ens                (num_ens,numpatch)); assimsun_out_ens  (:,:) = spval
            allocate (assimsha_out_ens                (num_ens,numpatch)); assimsha_out_ens  (:,:) = spval
            allocate (etrsun_out_ens                  (num_ens,numpatch)); etrsun_out_ens    (:,:) = spval
            allocate (etrsha_out_ens                  (num_ens,numpatch)); etrsha_out_ens    (:,:) = spval

            allocate (t_grnd_ens                      (num_ens,numpatch)); t_grnd_ens        (:,:) = spval
            allocate (tleaf_ens                       (num_ens,numpatch)); tleaf_ens         (:,:) = spval
            allocate (ldew_ens                        (num_ens,numpatch)); ldew_ens          (:,:) = spval
            allocate (ldew_rain_ens                   (num_ens,numpatch)); ldew_rain_ens     (:,:) = spval
            allocate (ldew_snow_ens                   (num_ens,numpatch)); ldew_snow_ens     (:,:) = spval
            allocate (fwet_snow_ens                   (num_ens,numpatch)); fwet_snow_ens     (:,:) = spval
            allocate (sag_ens                         (num_ens,numpatch)); sag_ens           (:,:) = spval
            allocate (scv_ens                         (num_ens,numpatch)); scv_ens           (:,:) = spval
            allocate (snowdp_ens                      (num_ens,numpatch)); snowdp_ens        (:,:) = spval
            allocate (fveg_ens                        (num_ens,numpatch)); fveg_ens          (:,:) = spval
            allocate (fsno_ens                        (num_ens,numpatch)); fsno_ens          (:,:) = spval
            allocate (sigf_ens                        (num_ens,numpatch)); sigf_ens          (:,:) = spval
            allocate (green_ens                       (num_ens,numpatch)); green_ens         (:,:) = spval
            allocate (tlai_ens                        (num_ens,numpatch)); tlai_ens          (:,:) = spval
            allocate (lai_ens                         (num_ens,numpatch)); lai_ens           (:,:) = spval
            allocate (laisun_ens                      (num_ens,numpatch)); laisun_ens        (:,:) = spval
            allocate (laisha_ens                      (num_ens,numpatch)); laisha_ens        (:,:) = spval
            allocate (tsai_ens                        (num_ens,numpatch)); tsai_ens          (:,:) = spval
            allocate (sai_ens                         (num_ens,numpatch)); sai_ens           (:,:) = spval
            allocate (coszen_ens                      (num_ens,numpatch)); coszen_ens        (:,:) = spval
            allocate (alb_ens                     (2,2,num_ens,numpatch)); alb_ens       (:,:,:,:) = spval
            allocate (ssun_ens                    (2,2,num_ens,numpatch)); ssun_ens      (:,:,:,:) = spval
            allocate (ssha_ens                    (2,2,num_ens,numpatch)); ssha_ens      (:,:,:,:) = spval
            allocate (ssoi_ens                    (2,2,num_ens,numpatch)); ssoi_ens      (:,:,:,:) = spval
            allocate (ssno_ens                    (2,2,num_ens,numpatch)); ssno_ens      (:,:,:,:) = spval
            allocate (thermk_ens                      (num_ens,numpatch)); thermk_ens        (:,:) = spval
            allocate (extkb_ens                       (num_ens,numpatch)); extkb_ens         (:,:) = spval
            allocate (extkd_ens                       (num_ens,numpatch)); extkd_ens         (:,:) = spval
            allocate (zwt_ens                         (num_ens,numpatch)); zwt_ens           (:,:) = spval
            allocate (wa_ens                          (num_ens,numpatch)); wa_ens            (:,:) = spval
            allocate (wetwat_ens                      (num_ens,numpatch)); wetwat_ens        (:,:) = spval
            allocate (wat_ens                         (num_ens,numpatch)); wat_ens           (:,:) = spval
            allocate (wdsrf_ens                       (num_ens,numpatch)); wdsrf_ens         (:,:) = spval
            allocate (rss_ens                         (num_ens,numpatch)); rss_ens           (:,:) = spval

            allocate (t_lake_ens              (nl_lake,num_ens,numpatch)); t_lake_ens      (:,:,:) = spval
            allocate (lake_icefrac_ens        (nl_lake,num_ens,numpatch)); lake_icefrac_ens(:,:,:) = spval
            allocate (savedtke1_ens                   (num_ens,numpatch)); savedtke1_ens     (:,:) = spval
            allocate (dz_lake_ens             (nl_lake,num_ens,numpatch)); t_lake_ens      (:,:,:) = spval

            allocate (snw_rds_ens          (maxsnl+1:0,num_ens,numpatch)); snw_rds_ens       (:,:,:) = spval
            allocate (mss_bcpho_ens        (maxsnl+1:0,num_ens,numpatch)); mss_bcpho_ens     (:,:,:) = spval
            allocate (mss_bcphi_ens        (maxsnl+1:0,num_ens,numpatch)); mss_bcphi_ens     (:,:,:) = spval
            allocate (mss_ocpho_ens        (maxsnl+1:0,num_ens,numpatch)); mss_ocpho_ens     (:,:,:) = spval
            allocate (mss_ocphi_ens        (maxsnl+1:0,num_ens,numpatch)); mss_ocphi_ens     (:,:,:) = spval
            allocate (mss_dst1_ens         (maxsnl+1:0,num_ens,numpatch)); mss_dst1_ens      (:,:,:) = spval
            allocate (mss_dst2_ens         (maxsnl+1:0,num_ens,numpatch)); mss_dst2_ens      (:,:,:) = spval
            allocate (mss_dst3_ens         (maxsnl+1:0,num_ens,numpatch)); mss_dst3_ens      (:,:,:) = spval
            allocate (mss_dst4_ens         (maxsnl+1:0,num_ens,numpatch)); mss_dst4_ens      (:,:,:) = spval
            allocate (ssno_lyr_ens     (2,2,maxsnl+1:1,num_ens,numpatch)); ssno_lyr_ens  (:,:,:,:,:) = spval

            allocate (trad_ens                        (num_ens,numpatch)); trad_ens            (:,:) = spval
            allocate (tref_ens                        (num_ens,numpatch)); tref_ens            (:,:) = spval
            allocate (qref_ens                        (num_ens,numpatch)); qref_ens            (:,:) = spval
            allocate (rst_ens                         (num_ens,numpatch)); rst_ens             (:,:) = spval
            allocate (emis_ens                        (num_ens,numpatch)); emis_ens            (:,:) = spval
            allocate (z0m_ens                         (num_ens,numpatch)); z0m_ens             (:,:) = spval
            allocate (displa_ens                      (num_ens,numpatch)); displa_ens          (:,:) = spval
            allocate (zol_ens                         (num_ens,numpatch)); zol_ens             (:,:) = spval
            allocate (rib_ens                         (num_ens,numpatch)); rib_ens             (:,:) = spval
            allocate (ustar_ens                       (num_ens,numpatch)); ustar_ens           (:,:) = spval
            allocate (qstar_ens                       (num_ens,numpatch)); qstar_ens           (:,:) = spval
            allocate (tstar_ens                       (num_ens,numpatch)); tstar_ens           (:,:) = spval
            allocate (fm_ens                          (num_ens,numpatch)); fm_ens              (:,:) = spval
            allocate (fh_ens                          (num_ens,numpatch)); fh_ens              (:,:) = spval
            allocate (fq_ens                          (num_ens,numpatch)); fq_ens              (:,:) = spval

            allocate (brt_temp_ens                  (2,num_ens,numpatch)); brt_temp_ens      (:,:,:) = spval

         ENDIF
      ENDIF

   END SUBROUTINE allocate_TimeVariables_ens

!-----------------------------------------------------------------------

   SUBROUTINE deallocate_TimeVariables_ens()

!-----------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_LandPatch, only: numpatch
      IMPLICIT NONE

!-----------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            deallocate (z_sno_ens                  )
            deallocate (dz_sno_ens                 )
            deallocate (t_soisno_ens               )
            deallocate (wliq_soisno_ens            )
            deallocate (wice_soisno_ens            )
            deallocate (h2osoi_ens                 )

            deallocate (smp_ens                    )
            deallocate (hk_ens                     )
            deallocate (rootr_ens                  )
            deallocate (rootflux_ens               )

            deallocate (vegwp_ens                  )
            deallocate (gs0sun_ens                 )
            deallocate (gs0sha_ens                 )

            deallocate (lai_old_ens                )
            deallocate (o3uptakesun_ens            )
            deallocate (o3uptakesha_ens            )

            deallocate (rstfacsun_out_ens          )
            deallocate (rstfacsha_out_ens          )
            deallocate (gssun_out_ens              )
            deallocate (gssha_out_ens              )
            deallocate (assimsun_out_ens           )
            deallocate (assimsha_out_ens           )
            deallocate (etrsun_out_ens             )
            deallocate (etrsha_out_ens             )

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
            deallocate (laisun_ens                 )
            deallocate (laisha_ens                 )
            deallocate (tsai_ens                   )
            deallocate (sai_ens                    )
            deallocate (coszen_ens                 )
            deallocate (alb_ens                    )
            deallocate (ssun_ens                   )
            deallocate (ssha_ens                   )
            deallocate (ssoi_ens                   )
            deallocate (ssno_ens                   )
            deallocate (thermk_ens                 )
            deallocate (extkb_ens                  )
            deallocate (extkd_ens                  )
            deallocate (zwt_ens                    )
            deallocate (wa_ens                     )
            deallocate (wetwat_ens                 )
            deallocate (wat_ens                    )
            deallocate (wdsrf_ens                  )
            deallocate (rss_ens                    )

            deallocate (t_lake_ens                 )
            deallocate (lake_icefrac_ens           )
            deallocate (savedtke1_ens              )
            deallocate (dz_lake_ens                )

            deallocate (snw_rds_ens                )
            deallocate (mss_bcpho_ens              )
            deallocate (mss_bcphi_ens              )
            deallocate (mss_ocpho_ens              )
            deallocate (mss_ocphi_ens              )
            deallocate (mss_dst1_ens               )
            deallocate (mss_dst2_ens               )
            deallocate (mss_dst3_ens               )
            deallocate (mss_dst4_ens               )
            deallocate (ssno_lyr_ens               )

            deallocate (trad_ens                   )
            deallocate (tref_ens                   )
            deallocate (qref_ens                   )
            deallocate (rst_ens                    )
            deallocate (emis_ens                   )
            deallocate (z0m_ens                    )
            deallocate (displa_ens                 )
            deallocate (zol_ens                    )
            deallocate (rib_ens                    )
            deallocate (ustar_ens                  )
            deallocate (qstar_ens                  )
            deallocate (tstar_ens                  )
            deallocate (fm_ens                     )
            deallocate (fh_ens                     )
            deallocate (fq_ens                     )

            deallocate (brt_temp_ens               )
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_TimeVariables_ens

!-----------------------------------------------------------------------

   SUBROUTINE WRITE_TimeVariables_ens (idate, lc_year, site, dir_restart)

!-----------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_Namelist, only : DEF_REST_CompressLevel, DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS, &
         DEF_USE_IRRIGATION, DEF_USE_Dynamic_Lake
      USE MOD_LandPatch
      USE MOD_NetCDFVector
      USE MOD_Vars_Global
      USE MOD_Vars_TimeInvariants, only : dz_lake
      IMPLICIT NONE

!-----------------------------------------------------------------------
      integer, intent(in) :: idate(3)
      integer, intent(in) :: lc_year      !year of land cover type data
      character(len=*), intent(in) :: site
      character(len=*), intent(in) :: dir_restart

!-----------------------------------------------------------------------
      ! Local variables
      character(len=256) :: file_restart
      character(len=14)  :: cdate
      character(len=256) :: cyear         !character for lc_year
      integer :: compress
      integer :: i

!-----------------------------------------------------------------------
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

      file_restart = trim(dir_restart)// '/'//trim(cdate)//'/' // trim(site) //'_restart_'//trim(cdate)//'_lc'//trim(cyear)//'_ens'//'.nc'

      CALL ncio_create_file_vector      (file_restart, landpatch)

      CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'snow',     -maxsnl       )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'snowp1',   -maxsnl+1     )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soilsnow', nl_soil-maxsnl)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil',     nl_soil)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'lake',     nl_lake)
      IF(DEF_USE_PLANTHYDRAULICS)THEN
         CALL ncio_define_dimension_vector (file_restart, landpatch, 'vegnodes', nvegwcs)
      ENDIF
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'band', 2)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'rtyp', 2)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'ens', num_ens)

      ! Time-varying state variables which reaquired by restart run
      CALL ncio_write_vector (file_restart, 'z_sno   '   , 'snow',     -maxsnl,        'ens', num_ens, 'patch', landpatch, z_sno_ens,       compress) ! node depth [m]
      CALL ncio_write_vector (file_restart, 'dz_sno  '   , 'snow',     -maxsnl,        'ens', num_ens, 'patch', landpatch, dz_sno_ens,      compress) ! interface depth [m]
      CALL ncio_write_vector (file_restart, 't_soisno'   , 'soilsnow', nl_soil-maxsnl, 'ens', num_ens, 'patch', landpatch, t_soisno_ens,    compress) ! soil temperature [K]
      CALL ncio_write_vector (file_restart, 'wliq_soisno', 'soilsnow', nl_soil-maxsnl, 'ens', num_ens, 'patch', landpatch, wliq_soisno_ens, compress) ! liquid water in layers [kg/m2]
      CALL ncio_write_vector (file_restart, 'wice_soisno', 'soilsnow', nl_soil-maxsnl, 'ens', num_ens, 'patch', landpatch, wice_soisno_ens, compress) ! ice lens in layers [kg/m2]
      CALL ncio_write_vector (file_restart, 'smp',         'soil',     nl_soil,        'ens', num_ens, 'patch', landpatch, smp_ens,         compress) ! soil matrix potential [mm]
      CALL ncio_write_vector (file_restart, 'hk',          'soil',     nl_soil,        'ens', num_ens, 'patch', landpatch, hk_ens,          compress) ! hydraulic conductivity [mm h2o/s]

      IF (DEF_USE_PLANTHYDRAULICS) THEN
         CALL ncio_write_vector (file_restart, 'vegwp',      'vegnodes',  nvegwcs,        'ens', num_ens, 'patch', landpatch, vegwp_ens,  compress) ! vegetation water potential [mm]
         CALL ncio_write_vector (file_restart, 'gs0sun',                                  'ens', num_ens, 'patch', landpatch, gs0sun_ens, compress) ! working copy of sunlit stomata conductance
         CALL ncio_write_vector (file_restart, 'gs0sha',                                  'ens', num_ens, 'patch', landpatch, gs0sha_ens, compress) ! working copy of shaded stomata conductance
      ENDIF

      IF (DEF_USE_OZONESTRESS) THEN
         CALL ncio_write_vector (file_restart, 'lai_old',                                 'ens', num_ens, 'patch', landpatch, lai_old_ens,     compress)
         CALL ncio_write_vector (file_restart, 'o3uptakesun',                             'ens', num_ens, 'patch', landpatch, o3uptakesun_ens, compress)
         CALL ncio_write_vector (file_restart, 'o3uptakesha',                             'ens', num_ens, 'patch', landpatch, o3uptakesha_ens, compress)
      ENDIF

      CALL ncio_write_vector (file_restart, 't_grnd',    'ens', num_ens, 'patch', landpatch, t_grnd_ens,    compress) ! ground surface temperature [K]
      CALL ncio_write_vector (file_restart, 'tleaf',     'ens', num_ens, 'patch', landpatch, tleaf_ens,     compress) ! leaf temperature [K]
      CALL ncio_write_vector (file_restart, 'ldew',      'ens', num_ens, 'patch', landpatch, ldew_ens,      compress) ! depth of water on foliage [mm]
      CALL ncio_write_vector (file_restart, 'ldew_rain', 'ens', num_ens, 'patch', landpatch, ldew_rain_ens, compress) ! depth of water on foliage [mm]
      CALL ncio_write_vector (file_restart, 'ldew_snow', 'ens', num_ens, 'patch', landpatch, ldew_snow_ens, compress) ! depth of water on foliage [mm]
      CALL ncio_write_vector (file_restart, 'fwet_snow', 'ens', num_ens, 'patch', landpatch, fwet_snow_ens, compress) ! vegetation snow fractional cover [-]
      CALL ncio_write_vector (file_restart, 'sag',       'ens', num_ens, 'patch', landpatch, sag_ens,       compress) ! non dimensional snow age [-]
      CALL ncio_write_vector (file_restart, 'scv',       'ens', num_ens, 'patch', landpatch, scv_ens,       compress) ! snow cover, water equivalent [mm]
      CALL ncio_write_vector (file_restart, 'snowdp',    'ens', num_ens, 'patch', landpatch, snowdp_ens,    compress) ! snow depth [meter]
      CALL ncio_write_vector (file_restart, 'fveg',      'ens', num_ens, 'patch', landpatch, fveg_ens,      compress) ! fraction of vegetation cover
      CALL ncio_write_vector (file_restart, 'fsno',      'ens', num_ens, 'patch', landpatch, fsno_ens,      compress) ! fraction of snow cover on ground
      CALL ncio_write_vector (file_restart, 'sigf',      'ens', num_ens, 'patch', landpatch, sigf_ens,      compress) ! fraction of veg cover, excluding snow-covered veg [-]
      CALL ncio_write_vector (file_restart, 'green',     'ens', num_ens, 'patch', landpatch, green_ens,     compress) ! leaf greenness
      CALL ncio_write_vector (file_restart, 'lai',       'ens', num_ens, 'patch', landpatch, lai_ens,       compress) ! leaf area index
      CALL ncio_write_vector (file_restart, 'tlai',      'ens', num_ens, 'patch', landpatch, tlai_ens,      compress) ! leaf area index
      CALL ncio_write_vector (file_restart, 'sai',       'ens', num_ens, 'patch', landpatch, sai_ens,       compress) ! stem area index
      CALL ncio_write_vector (file_restart, 'tsai',      'ens', num_ens, 'patch', landpatch, tsai_ens,      compress) ! stem area index
      CALL ncio_write_vector (file_restart, 'coszen',    'ens', num_ens, 'patch', landpatch, coszen_ens,    compress) ! cosine of solar zenith angle
      CALL ncio_write_vector (file_restart, 'alb     '   , 'band', 2, 'rtyp', 2, 'ens', num_ens, 'patch', landpatch, alb_ens,  compress) ! averaged albedo [-]
      CALL ncio_write_vector (file_restart, 'ssun    '   , 'band', 2, 'rtyp', 2, 'ens', num_ens, 'patch', landpatch, ssun_ens, compress) ! sunlit canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'ssha    '   , 'band', 2, 'rtyp', 2, 'ens', num_ens, 'patch', landpatch, ssha_ens, compress) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'ssoi    '   , 'band', 2, 'rtyp', 2, 'ens', num_ens, 'patch', landpatch, ssoi_ens, compress) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'ssno    '   , 'band', 2, 'rtyp', 2, 'ens', num_ens, 'patch', landpatch, ssno_ens, compress) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'thermk  ',  'ens', num_ens, 'patch', landpatch, thermk_ens,    compress) ! canopy gap fraction for tir radiation
      CALL ncio_write_vector (file_restart, 'extkb',     'ens', num_ens, 'patch', landpatch, extkb_ens,     compress) ! (k, g(mu)/mu) direct solar extinction coefficient
      CALL ncio_write_vector (file_restart, 'extkd',     'ens', num_ens, 'patch', landpatch, extkd_ens,     compress) ! diffuse and scattered diffuse PAR extinction coefficient
      CALL ncio_write_vector (file_restart, 'zwt',       'ens', num_ens, 'patch', landpatch, zwt_ens,       compress) ! the depth to water table [m]
      CALL ncio_write_vector (file_restart, 'wa',        'ens', num_ens, 'patch', landpatch, wa_ens,        compress) ! water storage in aquifer [mm]
      CALL ncio_write_vector (file_restart, 'wetwat',    'ens', num_ens, 'patch', landpatch, wetwat_ens,    compress) ! water storage in wetland [mm]
      CALL ncio_write_vector (file_restart, 'wdsrf',     'ens', num_ens, 'patch', landpatch, wdsrf_ens,     compress) ! depth of surface water [mm]
      CALL ncio_write_vector (file_restart, 'rss',       'ens', num_ens, 'patch', landpatch, rss_ens,       compress) ! soil surface resistance [s/m]

      IF (DEF_USE_Dynamic_Lake) THEN
         DO i = 1, num_ens
            dz_lake_ens(:,i,:) = dz_lake
         ENDDO
         CALL ncio_write_vector (file_restart, 'dz_lake'    , 'lake', nl_lake, 'ens', num_ens, 'patch', landpatch,  dz_lake_ens, compress)
      ENDIF
      CALL ncio_write_vector (file_restart, 't_lake  '   , 'lake', nl_lake, 'ens', num_ens, 'patch', landpatch, t_lake_ens,       compress)
      CALL ncio_write_vector (file_restart, 'lake_icefrc', 'lake', nl_lake, 'ens', num_ens, 'patch', landpatch, lake_icefrac_ens, compress)
      CALL ncio_write_vector (file_restart, 'savedtke1  ', 'ens', num_ens, 'patch', landpatch, savedtke1_ens, compress)

      CALL ncio_write_vector (file_restart, 'snw_rds  ', 'snow', -maxsnl, 'ens', num_ens, 'patch', landpatch, snw_rds_ens,   compress)
      CALL ncio_write_vector (file_restart, 'mss_bcpho', 'snow', -maxsnl, 'ens', num_ens, 'patch', landpatch, mss_bcpho_ens, compress)
      CALL ncio_write_vector (file_restart, 'mss_bcphi', 'snow', -maxsnl, 'ens', num_ens, 'patch', landpatch, mss_bcphi_ens, compress)
      CALL ncio_write_vector (file_restart, 'mss_ocpho', 'snow', -maxsnl, 'ens', num_ens, 'patch', landpatch, mss_ocpho_ens, compress)
      CALL ncio_write_vector (file_restart, 'mss_ocphi', 'snow', -maxsnl, 'ens', num_ens, 'patch', landpatch, mss_ocphi_ens, compress)
      CALL ncio_write_vector (file_restart, 'mss_dst1 ', 'snow', -maxsnl, 'ens', num_ens, 'patch', landpatch, mss_dst1_ens,  compress)
      CALL ncio_write_vector (file_restart, 'mss_dst2 ', 'snow', -maxsnl, 'ens', num_ens, 'patch', landpatch, mss_dst2_ens,  compress)
      CALL ncio_write_vector (file_restart, 'mss_dst3 ', 'snow', -maxsnl, 'ens', num_ens, 'patch', landpatch, mss_dst3_ens,  compress)
      CALL ncio_write_vector (file_restart, 'mss_dst4 ', 'snow', -maxsnl, 'ens', num_ens, 'patch', landpatch, mss_dst4_ens,  compress)
      CALL ncio_write_vector (file_restart, 'ssno_lyr', 'band', 2, 'rtyp', 2, 'snowp1', -maxsnl+1, 'ens', num_ens, 'patch', landpatch, ssno_lyr_ens, compress) ! TODO: Add write vector for 5d

      CALL ncio_write_vector (file_restart, 'trad ', 'ens', num_ens, 'patch', landpatch, trad_ens, compress) ! radiative temperature of surface [K]
      CALL ncio_write_vector (file_restart, 'tref ', 'ens', num_ens, 'patch', landpatch, tref_ens, compress) ! 2 m height air temperature [kelvin]
      CALL ncio_write_vector (file_restart, 'qref ', 'ens', num_ens, 'patch', landpatch, qref_ens, compress) ! 2 m height air specific humidity
      CALL ncio_write_vector (file_restart, 'rst  ', 'ens', num_ens, 'patch', landpatch, rst_ens,  compress) ! canopy stomatal resistance (s/m)
      CALL ncio_write_vector (file_restart, 'emis ', 'ens', num_ens, 'patch', landpatch, emis_ens, compress) ! averaged bulk surface emissivity
      CALL ncio_write_vector (file_restart, 'z0m  ', 'ens', num_ens, 'patch', landpatch, z0m_ens,  compress) ! effective roughness [m]
      CALL ncio_write_vector (file_restart, 'zol  ', 'ens', num_ens, 'patch', landpatch, zol_ens,  compress) ! dimensionless height (z/L) used in Monin-Obukhov theory
      CALL ncio_write_vector (file_restart, 'rib  ', 'ens', num_ens, 'patch', landpatch, rib_ens,  compress) ! bulk Richardson number in surface layer
      CALL ncio_write_vector (file_restart, 'ustar', 'ens', num_ens, 'patch', landpatch, ustar_ens,compress) ! u* in similarity theory [m/s]
      CALL ncio_write_vector (file_restart, 'qstar', 'ens', num_ens, 'patch', landpatch, qstar_ens,compress) ! q* in similarity theory [kg/kg]
      CALL ncio_write_vector (file_restart, 'tstar', 'ens', num_ens, 'patch', landpatch, tstar_ens,compress) ! t* in similarity theory [K]
      CALL ncio_write_vector (file_restart, 'fm   ', 'ens', num_ens, 'patch', landpatch, fm_ens,   compress) ! integral of profile FUNCTION for momentum
      CALL ncio_write_vector (file_restart, 'fh   ', 'ens', num_ens, 'patch', landpatch, fh_ens,   compress) ! integral of profile FUNCTION for heat
      CALL ncio_write_vector (file_restart, 'fq   ', 'ens', num_ens, 'patch', landpatch, fq_ens,   compress) ! integral of profile FUNCTION for moisture

   END SUBROUTINE WRITE_TimeVariables_ens

!-----------------------------------------------------------------------

   SUBROUTINE READ_TimeVariables_ens(idate, lc_year, site, dir_restart)

!-----------------------------------------------------------------------
      USE MOD_Namelist
      USE MOD_SPMD_Task
      USE MOD_NetCDFVector
#ifdef RangeCheck
      USE MOD_RangeCheck
#endif
      USE MOD_LandPatch
      USE MOD_Vars_Global
      USE MOD_Vars_TimeInvariants, only: dz_lake
      IMPLICIT NONE

!-----------------------------------------------------------------------
      integer, intent(in) :: idate(3)
      integer, intent(in) :: lc_year      !year of land cover type data
      character(len=*), intent(in) :: site
      character(len=*), intent(in) :: dir_restart

!-----------------------------------------------------------------------
      character(len=256) :: file_restart
      character(len=14)  :: cdate, cyear

!-----------------------------------------------------------------------
#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         write (*, *) 'Loading Ensemble Time Variables ...'
      END IF

      ! land cover type year
      write (cyear, '(i4.4)') lc_year

      write (cdate, '(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
      file_restart = trim(dir_restart)//'/'//trim(cdate)//'/'//trim(site)//'_restart_'//trim(cdate)//'_lc'//trim(cyear)//'_ens'//'.nc'

      ! Time-varying state variables which reaquired by restart run
      CALL ncio_read_vector(file_restart, 'z_sno   ', -maxsnl, num_ens, landpatch, z_sno_ens)             ! node depth [m]
      CALL ncio_read_vector(file_restart, 'dz_sno  ', -maxsnl, num_ens, landpatch, dz_sno_ens)            ! interface depth [m]
      CALL ncio_read_vector(file_restart, 't_soisno', nl_soil - maxsnl, num_ens, landpatch, t_soisno_ens)   ! soil temperature [K]
      CALL ncio_read_vector(file_restart, 'wliq_soisno', nl_soil - maxsnl, num_ens, landpatch, wliq_soisno_ens)! liquid water in layers [kg/m2]
      CALL ncio_read_vector(file_restart, 'wice_soisno', nl_soil - maxsnl, num_ens, landpatch, wice_soisno_ens)! ice lens in layers [kg/m2]
      CALL ncio_read_vector(file_restart, 'smp', nl_soil, num_ens, landpatch, smp_ens)        ! soil matrix potential [mm]
      CALL ncio_read_vector(file_restart, 'hk', nl_soil, num_ens, landpatch, hk_ens)         ! hydraulic conductivity [mm h2o/s]

      IF (DEF_USE_PLANTHYDRAULICS) THEN
         CALL ncio_read_vector(file_restart, 'vegwp', nvegwcs, num_ens, landpatch, vegwp_ens) ! vegetation water potential [mm]
         CALL ncio_read_vector(file_restart, 'gs0sun  ', num_ens, landpatch, gs0sun_ens) ! working copy of sunlit stomata conductance
         CALL ncio_read_vector(file_restart, 'gs0sha  ', num_ens, landpatch, gs0sha_ens) ! working copy of shalit stomata conductance
      END IF

      IF (DEF_USE_OZONESTRESS) THEN
         CALL ncio_read_vector(file_restart, 'lai_old    ', num_ens, landpatch, lai_old_ens)
         CALL ncio_read_vector(file_restart, 'o3uptakesun', num_ens, landpatch, o3uptakesun_ens)
         CALL ncio_read_vector(file_restart, 'o3uptakesha', num_ens, landpatch, o3uptakesha_ens)
      END IF

      CALL ncio_read_vector(file_restart, 't_grnd  ', num_ens, landpatch, t_grnd_ens) ! ground surface temperature [K]
      CALL ncio_read_vector(file_restart, 'tleaf   ', num_ens, landpatch, tleaf_ens) ! leaf temperature [K]
      CALL ncio_read_vector(file_restart, 'ldew    ', num_ens, landpatch, ldew_ens) ! depth of water on foliage [mm]
      CALL ncio_read_vector(file_restart, 'ldew_rain', num_ens, landpatch, ldew_rain_ens) ! depth of water on foliage [mm]
      CALL ncio_read_vector(file_restart, 'ldew_snow', num_ens, landpatch, ldew_snow_ens) ! depth of water on foliage [mm]
      CALL ncio_read_vector(file_restart, 'fwet_snow', num_ens, landpatch, fwet_snow_ens) ! vegetation snow fractional cover [-]
      CALL ncio_read_vector(file_restart, 'sag     ', num_ens, landpatch, sag_ens) ! non dimensional snow age [-]
      CALL ncio_read_vector(file_restart, 'scv     ', num_ens, landpatch, scv_ens) ! snow cover, water equivalent [mm]
      CALL ncio_read_vector(file_restart, 'snowdp  ', num_ens, landpatch, snowdp_ens) ! snow depth [meter]
      CALL ncio_read_vector(file_restart, 'fveg    ', num_ens, landpatch, fveg_ens) ! fraction of vegetation cover
      CALL ncio_read_vector(file_restart, 'fsno    ', num_ens, landpatch, fsno_ens) ! fraction of snow cover on ground
      CALL ncio_read_vector(file_restart, 'sigf    ', num_ens, landpatch, sigf_ens) ! fraction of veg cover, excluding snow-covered veg [-]
      CALL ncio_read_vector(file_restart, 'green   ', num_ens, landpatch, green_ens) ! leaf greenness
      CALL ncio_read_vector(file_restart, 'lai     ', num_ens, landpatch, lai_ens) ! leaf area index
      CALL ncio_read_vector(file_restart, 'tlai    ', num_ens, landpatch, tlai_ens) ! leaf area index
      CALL ncio_read_vector(file_restart, 'sai     ', num_ens, landpatch, sai_ens) ! stem area index
      CALL ncio_read_vector(file_restart, 'tsai    ', num_ens, landpatch, tsai_ens) ! stem area index
      CALL ncio_read_vector(file_restart, 'coszen  ', num_ens, landpatch, coszen_ens) ! cosine of solar zenith angle
      CALL ncio_read_vector(file_restart, 'alb     ', 2, 2, num_ens, landpatch, alb_ens) ! averaged albedo [-]
      CALL ncio_read_vector(file_restart, 'ssun    ', 2, 2, num_ens, landpatch, ssun_ens) ! sunlit canopy absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'ssha    ', 2, 2, num_ens, landpatch, ssha_ens) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'ssoi    ', 2, 2, num_ens, landpatch, ssoi_ens) ! soil absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'ssno    ', 2, 2, num_ens, landpatch, ssno_ens) ! snow absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'thermk  ', num_ens, landpatch, thermk_ens) ! canopy gap fraction for tir radiation
      CALL ncio_read_vector(file_restart, 'extkb   ', num_ens, landpatch, extkb_ens) ! (k, g(mu)/mu) direct solar extinction coefficient
      CALL ncio_read_vector(file_restart, 'extkd   ', num_ens, landpatch, extkd_ens) ! diffuse and scattered diffuse PAR extinction coefficient
      CALL ncio_read_vector(file_restart, 'zwt     ', num_ens, landpatch, zwt_ens) ! the depth to water table [m]
      CALL ncio_read_vector(file_restart, 'wa      ', num_ens, landpatch, wa_ens) ! water storage in aquifer [mm]
      CALL ncio_read_vector(file_restart, 'wetwat  ', num_ens, landpatch, wetwat_ens) ! water storage in wetland [mm]
      CALL ncio_read_vector(file_restart, 'wdsrf   ', num_ens, landpatch, wdsrf_ens) ! depth of surface water [mm]
      CALL ncio_read_vector(file_restart, 'rss     ', num_ens, landpatch, rss_ens) ! soil surface resistance [s/m]

      IF (DEF_USE_Dynamic_Lake) THEN
         CALL ncio_read_vector(file_restart, 'dz_lake', nl_lake, num_ens, landpatch, dz_lake_ens)
      END IF
      CALL ncio_read_vector(file_restart, 't_lake  ', nl_lake, num_ens, landpatch, t_lake_ens) ! lake temperature [K]
      CALL ncio_read_vector(file_restart, 'lake_icefrc', nl_lake, num_ens, landpatch, lake_icefrac_ens) ! lake ice fraction [-]
      CALL ncio_read_vector(file_restart, 'savedtke1  ', num_ens, landpatch, savedtke1_ens) ! saved tke1 [m2/s2]

      CALL ncio_read_vector(file_restart, 'snw_rds  ', -maxsnl, num_ens, landpatch, snw_rds_ens) ! snow radius [m]
      CALL ncio_read_vector(file_restart, 'mss_bcpho', -maxsnl, num_ens, landpatch, mss_bcpho_ens) ! mass of BC in snow [kg/m2]
      CALL ncio_read_vector(file_restart, 'mss_bcphi', -maxsnl, num_ens, landpatch, mss_bcphi_ens) ! mass of BC in snow [kg/m2]
      CALL ncio_read_vector(file_restart, 'mss_ocpho', -maxsnl, num_ens, landpatch, mss_ocpho_ens) ! mass of OC in snow [kg/m2]
      CALL ncio_read_vector(file_restart, 'mss_ocphi', -maxsnl, num_ens, landpatch, mss_ocphi_ens) ! mass of OC in snow [kg/m2]
      CALL ncio_read_vector(file_restart, 'mss_dst1 ', -maxsnl, num_ens, landpatch, mss_dst1_ens) ! mass of DST1 in snow [kg/m2]
      CALL ncio_read_vector(file_restart, 'mss_dst2 ', -maxsnl, num_ens, landpatch, mss_dst2_ens) ! mass of DST2 in snow [kg/m2]
      CALL ncio_read_vector(file_restart, 'mss_dst3 ', -maxsnl, num_ens, landpatch, mss_dst3_ens) ! mass of DST3 in snow [kg/m2]
      CALL ncio_read_vector(file_restart, 'mss_dst4 ', -maxsnl, num_ens, landpatch, mss_dst4_ens) ! mass of DST4 in snow [kg/m2]
      CALL ncio_read_vector(file_restart, 'ssno_lyr', 2, 2, -maxsnl + 1, num_ens, landpatch, ssno_lyr_ens) ! snow layer thickness [m]

      CALL ncio_read_vector(file_restart, 'trad ', num_ens, landpatch, trad_ens) ! radiative temperature of surface [K]
      CALL ncio_read_vector(file_restart, 'tref ', num_ens, landpatch, tref_ens) ! 2 m height air temperature [kelvin]
      CALL ncio_read_vector(file_restart, 'qref ', num_ens, landpatch, qref_ens) ! 2 m height air specific humidity
      CALL ncio_read_vector(file_restart, 'rst  ', num_ens, landpatch, rst_ens) ! canopy stomatal resistance (s/m)
      CALL ncio_read_vector(file_restart, 'emis ', num_ens, landpatch, emis_ens) ! averaged bulk surface emissivity
      CALL ncio_read_vector(file_restart, 'z0m  ', num_ens, landpatch, z0m_ens) ! effective roughness [m]
      CALL ncio_read_vector(file_restart, 'zol  ', num_ens, landpatch, zol_ens) ! dimensionless height (z/L) used in Monin-Obukhov theory
      CALL ncio_read_vector(file_restart, 'rib  ', num_ens, landpatch, rib_ens) ! bulk Richardson number in surface layer
      CALL ncio_read_vector(file_restart, 'ustar', num_ens, landpatch, ustar_ens) ! u* in similarity theory [m/s]
      CALL ncio_read_vector(file_restart, 'qstar', num_ens, landpatch, qstar_ens) ! q* in similarity theory [kg/kg]
      CALL ncio_read_vector(file_restart, 'tstar', num_ens, landpatch, tstar_ens) ! t* in similarity theory [K]
      CALL ncio_read_vector(file_restart, 'fm   ', num_ens, landpatch, fm_ens) ! integral of profile FUNCTION for momentum
      CALL ncio_read_vector(file_restart, 'fh   ', num_ens, landpatch, fh_ens) ! integral of profile FUNCTION for heat
      CALL ncio_read_vector(file_restart, 'fq   ', num_ens, landpatch, fq_ens) ! integral of profile FUNCTION for moisture

#ifdef RangeCheck
      CALL check_TimeVariables_ens
#endif

      IF (p_is_master) THEN
         write (*, *) 'Loading Ensemble Time Variables done.'
      END IF

   END SUBROUTINE READ_TimeVariables_ens

#ifdef RangeCheck
!-----------------------------------------------------------------------

   SUBROUTINE check_TimeVariables_ens()

!-----------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_RangeCheck
      USE MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS, DEF_USE_IRRIGATION, &
         DEF_USE_SNICAR, DEF_USE_Dynamic_Lake
      USE MOD_Vars_TimeInvariants, only: dz_lake
      IMPLICIT NONE

#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write (*, *) 'Checking Ensemble Time Variables ...'
      END IF

      CALL check_vector_data ('z_sno       [m]    ', z_sno_ens      ) ! node depth [m]
      CALL check_vector_data ('dz_sno      [m]    ', dz_sno_ens     ) ! interface depth [m]
      CALL check_vector_data ('t_soisno    [K]    ', t_soisno_ens   ) ! soil temperature [K]
      CALL check_vector_data ('wliq_soisno [kg/m2]', wliq_soisno_ens) ! liquid water in layers [kg/m2]
      CALL check_vector_data ('wice_soisno [kg/m2]', wice_soisno_ens) ! ice lens in layers [kg/m2]
      CALL check_vector_data ('smp         [mm]   ', smp_ens        ) ! soil matrix potential [mm]
      CALL check_vector_data ('hk          [mm/s] ', hk_ens         ) ! hydraulic conductivity [mm h2o/s]

      IF (DEF_USE_PLANTHYDRAULICS) THEN
         CALL check_vector_data ('vegwp       [m]    ', vegwp_ens      ) ! vegetation water potential [mm]
         CALL check_vector_data ('gs0sun      []     ', gs0sun_ens     ) ! working copy of sunlit stomata conductance
         CALL check_vector_data ('gs0sha      []     ', gs0sha_ens     ) ! working copy of shalit stomata conductance
      ENDIF

      IF (DEF_USE_OZONESTRESS) THEN
         CALL check_vector_data ('lai_old     []     ', lai_old_ens    )
         CALL check_vector_data ('o3uptakesun []     ', o3uptakesun_ens)
         CALL check_vector_data ('o3uptakesha []     ', o3uptakesha_ens)
      ENDIF

      CALL check_vector_data ('t_grnd      [K]    ', t_grnd_ens     ) ! ground surface temperature [K]
      CALL check_vector_data ('tleaf       [K]    ', tleaf_ens      ) ! leaf temperature [K]
      CALL check_vector_data ('ldew        [mm]   ', ldew_ens       ) ! depth of water on foliage [mm]
      CALL check_vector_data ('ldew_rain   [mm]   ', ldew_rain_ens  ) ! depth of rain on foliage [mm]
      CALL check_vector_data ('ldew_snow   [mm]   ', ldew_snow_ens  ) ! depth of snow on foliage [mm]
      CALL check_vector_data ('fwet_snow   [mm]   ', fwet_snow_ens  ) ! vegetation snow fractional cover [-]
      CALL check_vector_data ('sag         [-]    ', sag_ens        ) ! non dimensional snow age [-]
      CALL check_vector_data ('scv         [mm]   ', scv_ens        ) ! snow cover, water equivalent [mm]
      CALL check_vector_data ('snowdp      [m]    ', snowdp_ens     ) ! snow depth [meter]
      CALL check_vector_data ('fveg        [-]    ', fveg_ens       ) ! fraction of vegetation cover
      CALL check_vector_data ('fsno        [-]    ', fsno_ens       ) ! fraction of snow cover on ground
      CALL check_vector_data ('sigf        [-]    ', sigf_ens       ) ! fraction of veg cover, excluding snow-covered veg [-]
      CALL check_vector_data ('green       [-]    ', green_ens      ) ! leaf greenness
      CALL check_vector_data ('lai         [-]    ', lai_ens        ) ! leaf area index
      CALL check_vector_data ('tlai        [-]    ', tlai_ens       ) ! leaf area index
      CALL check_vector_data ('sai         [-]    ', sai_ens        ) ! stem area index
      CALL check_vector_data ('tsai        [-]    ', tsai_ens       ) ! stem area index
      CALL check_vector_data ('coszen      [-]    ', coszen_ens     ) ! cosine of solar zenith angle
      CALL check_vector_data ('alb         [-]    ', alb_ens        ) ! averaged albedo [-]
      CALL check_vector_data ('ssun        [-]    ', ssun_ens       ) ! sunlit canopy absorption for solar radiation (0-1)
      CALL check_vector_data ('ssha        [-]    ', ssha_ens       ) ! shaded canopy absorption for solar radiation (0-1)
      CALL check_vector_data ('ssoi        [-]    ', ssoi_ens       ) ! soil absorption for solar radiation (0-1)
      CALL check_vector_data ('ssno        [-]    ', ssno_ens       ) ! snow absorption for solar radiation (0-1)
      CALL check_vector_data ('thermk      [-]    ', thermk_ens     ) ! canopy gap fraction for tir radiation
      CALL check_vector_data ('extkb       [-]    ', extkb_ens      ) ! (k, g(mu)/mu) direct solar extinction coefficient
      CALL check_vector_data ('extkd       [-]    ', extkd_ens      ) ! diffuse and scattered diffuse PAR extinction coefficient
      CALL check_vector_data ('zwt         [m]    ', zwt_ens        ) ! the depth to water table [m]
      CALL check_vector_data ('wa          [mm]   ', wa_ens         ) ! water storage in aquifer [mm]
      CALL check_vector_data ('wetwat      [mm]   ', wetwat_ens     ) ! water storage in wetland [mm]
      CALL check_vector_data ('wdsrf       [mm]   ', wdsrf_ens      ) ! depth of surface water [mm]
      CALL check_vector_data ('rss         [s/m]  ', rss_ens        ) ! soil surface resistance [s/m]

      IF (DEF_USE_Dynamic_Lake) THEN
         CALL check_vector_data ('dz_lake     [m]    ', dz_lake_ens    )!
      ENDIF
      CALL check_vector_data ('t_lake      [K]    ', t_lake_ens     ) ! lake temperature [K]
      CALL check_vector_data ('lake_icefrc [-]    ', lake_icefrac_ens) ! lake ice fraction [-]
      CALL check_vector_data ('savedtke1   [W/m K]', savedtke1_ens   ) ! saved tke1 [m2/s2]

      IF (DEF_USE_SNICAR) THEN
         CALL check_vector_data ('snw_rds     [m-6]  ',  snw_rds_ens   )
         CALL check_vector_data ('mss_bcpho   [Kg]   ',  mss_bcpho_ens )
         CALL check_vector_data ('mss_bcphi   [Kg]   ',  mss_bcphi_ens )
         CALL check_vector_data ('mss_ocpho   [Kg]   ',  mss_ocpho_ens )
         CALL check_vector_data ('mss_ocphi   [Kg]   ',  mss_ocphi_ens )
         CALL check_vector_data ('mss_dst1    [Kg]   ',  mss_dst1_ens  )
         CALL check_vector_data ('mss_dst2    [Kg]   ',  mss_dst2_ens  )
         CALL check_vector_data ('mss_dst3    [Kg]   ',  mss_dst3_ens  )
         CALL check_vector_data ('mss_dst4    [Kg]   ',  mss_dst4_ens  )
         CALL check_vector_data ('ssno_lyr    [-]    ',  ssno_lyr_ens  )
      ENDIF

      CALL check_vector_data ('brt_temp [K]', brt_temp_ens)

#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif

   END SUBROUTINE check_TimeVariables_ens
#endif

END MODULE MOD_DA_Vars_TimeVariables
! ---------- EOP ------------
