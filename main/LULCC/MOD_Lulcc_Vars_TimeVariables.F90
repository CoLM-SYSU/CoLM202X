#include <define.h>

MODULE MOD_Lulcc_Vars_TimeVariables

!-----------------------------------------------------------------------
!
!  Created by Hua Yuan, 04/2022
!
! !REVISIONS:
!  07/2023, Wenzong Dong: porting to MPI version
!  08/2023, Hua Yuan: unified PFT and PC process
!  10/2023, Wanyi Lin: check with MOD_Vars_TimeVariables.F90, add
!           variables, and remove unnecessary variables
!
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   IMPLICIT NONE
   SAVE
! ----------------------------------------------------------------------
! Time-varying state variables which required by restart run
   !TODO: need to check with MOD_Vars_TimeVariables.F90 whether
   !      there are any variables missing. - DONE
   real(r8), allocatable :: z_sno_         (:,:)  !node depth [m]
   real(r8), allocatable :: dz_sno_        (:,:)  !interface depth [m]
   real(r8), allocatable :: t_soisno_      (:,:)  !soil temperature [K]
   real(r8), allocatable :: wliq_soisno_   (:,:)  !liquid water in layers [kg/m2]
   real(r8), allocatable :: wice_soisno_   (:,:)  !ice lens in layers [kg/m2]
   real(r8), allocatable :: smp_           (:,:)  !soil matrix potential [mm]
   real(r8), allocatable :: hk_            (:,:)  !hydraulic conductivity [mm h2o/s]
   real(r8), allocatable :: t_grnd_          (:)  !ground surface temperature [K]

   real(r8), allocatable :: tleaf_           (:)  !leaf temperature [K]
   real(r8), allocatable :: ldew_            (:)  !depth of water on foliage [mm]
   real(r8), allocatable :: ldew_rain_       (:)  !depth of rain on foliage [mm]
   real(r8), allocatable :: ldew_snow_       (:)  !depth of rain on foliage [mm]
   real(r8), allocatable :: fwet_snow_       (:)  !vegetation snow fractional cover [-]
   real(r8), allocatable :: sag_             (:)  !non dimensional snow age [-]
   real(r8), allocatable :: scv_             (:)  !snow cover, water equivalent [mm]
   real(r8), allocatable :: snowdp_          (:)  !snow depth [meter]
   real(r8), allocatable :: fsno_            (:)  !frac of snow cover on ground
   real(r8), allocatable :: sigf_            (:)  !frac of veg cover, excluding snow-covered veg [-]
   real(r8), allocatable :: zwt_             (:)  !the depth to water table [m]
   real(r8), allocatable :: wa_              (:)  !water storage in aquifer [mm]
   real(r8), allocatable :: wdsrf_           (:)  !depth of surface water [mm]
   real(r8), allocatable :: rss_             (:)  !soil surface resistance [s/m]

   real(r8), allocatable :: t_lake_        (:,:)  !lake layer temperature [K]
   real(r8), allocatable :: lake_icefrac_  (:,:)  !lake mass fraction of lake layer that is frozen
   real(r8), allocatable :: savedtke1_       (:)  !top level eddy conductivity (W/m K)

   !Plant Hydraulic variables
   real(r8), allocatable :: vegwp_         (:,:)  !vegetation water potential [mm]
   real(r8), allocatable :: gs0sun_          (:)  !working copy of sunlit stomata conductance
   real(r8), allocatable :: gs0sha_          (:)  !working copy of shaded stomata conductance
   !END plant hydraulic variables

   !Ozone stress variables
   real(r8), allocatable :: lai_old_         (:)  !lai in last time step
   real(r8), allocatable :: o3uptakesun_     (:)  !Ozone does, sunlit leaf (mmol O3/m^2)
   real(r8), allocatable :: o3uptakesha_     (:)  !Ozone does, shaded leaf (mmol O3/m^2)
   !End ozone stress variables

   real(r8), allocatable :: snw_rds_       (:,:)  !effective grain radius (col,lyr) [microns, m-6]
   real(r8), allocatable :: mss_bcpho_     (:,:)  !mass of hydrophobic BC in snow  (col,lyr) [kg]
   real(r8), allocatable :: mss_bcphi_     (:,:)  !mass of hydrophillic BC in snow (col,lyr) [kg]
   real(r8), allocatable :: mss_ocpho_     (:,:)  !mass of hydrophobic OC in snow  (col,lyr) [kg]
   real(r8), allocatable :: mss_ocphi_     (:,:)  !mass of hydrophillic OC in snow (col,lyr) [kg]
   real(r8), allocatable :: mss_dst1_      (:,:)  !mass of dust species 1 in snow  (col,lyr) [kg]
   real(r8), allocatable :: mss_dst2_      (:,:)  !mass of dust species 2 in snow  (col,lyr) [kg]
   real(r8), allocatable :: mss_dst3_      (:,:)  !mass of dust species 3 in snow  (col,lyr) [kg]
   real(r8), allocatable :: mss_dst4_      (:,:)  !mass of dust species 4 in snow  (col,lyr) [kg]
   real(r8), allocatable :: ssno_lyr_  (:,:,:,:)  !snow layer absorption [-]

   ! Additional variables required by regional model (such as WRF ) RSM)
   real(r8), allocatable :: trad_            (:)  !radiative temperature of surface [K]
   real(r8), allocatable :: tref_            (:)  !2 m height air temperature [kelvin]
   real(r8), allocatable :: qref_            (:)  !2 m height air specific humidity
   real(r8), allocatable :: rst_             (:)  !canopy stomatal resistance (s/m)
   real(r8), allocatable :: emis_            (:)  !averaged bulk surface emissivity
   real(r8), allocatable :: z0m_             (:)  !effective roughness [m]
   real(r8), allocatable :: displa_          (:)  !zero displacement height [m]
   real(r8), allocatable :: zol_             (:)  !dimensionless height (z/L) used in M-O theory
   real(r8), allocatable :: rib_             (:)  !bulk Richardson number in surface layer
   real(r8), allocatable :: ustar_           (:)  !u* in similarity theory [m/s]
   real(r8), allocatable :: qstar_           (:)  !q* in similarity theory [kg/kg]
   real(r8), allocatable :: tstar_           (:)  !t* in similarity theory [K]
   real(r8), allocatable :: fm_              (:)  !integral of profile function for momentum
   real(r8), allocatable :: fh_              (:)  !integral of profile function for heat
   real(r8), allocatable :: fq_              (:)  !integral of profile function for moisture

   real(r8), allocatable :: sum_irrig_       (:)  !total irrigation amount [kg/m2]
   real(r8), allocatable :: sum_irrig_count_ (:)  !total irrigation counts [-]

   ! for LULC_IGBP_PFT and LULC_IGBP_PC
   real(r8), allocatable :: tleaf_p_         (:)  !shaded leaf temperature [K]
   real(r8), allocatable :: ldew_p_          (:)  !depth of water on foliage [mm]
   real(r8), allocatable :: ldew_rain_p_     (:)  !depth of rain on foliage [mm]
   real(r8), allocatable :: ldew_snow_p_     (:)  !depth of snow on foliage [mm]
   real(r8), allocatable :: fwet_snow_p_     (:)  !vegetation snow fractional cover [-]
   real(r8), allocatable :: sigf_p_          (:)  !frac of veg cover, excluding snow-covered veg [-]

   !TODO@yuan: to check the below for PC whether they are needed
   real(r8), allocatable :: tref_p_          (:)  !2 m height air temperature [kelvin]
   real(r8), allocatable :: qref_p_          (:)  !2 m height air specific humidity
   real(r8), allocatable :: rst_p_           (:)  !canopy stomatal resistance (s/m)
   real(r8), allocatable :: z0m_p_           (:)  !effective roughness [m]

   ! Plant Hydraulic variables
   real(r8), allocatable :: vegwp_p_       (:,:)  !vegetation water potential [mm]
   real(r8), allocatable :: gs0sun_p_        (:)  !working copy of sunlit stomata conductance
   real(r8), allocatable :: gs0sha_p_        (:)  !working copy of shaded stomata conductance
   ! end plant hydraulic variables

   ! Ozone Stress Variables
   real(r8), allocatable :: lai_old_p_       (:)  !lai in last time step
   real(r8), allocatable :: o3uptakesun_p_   (:)  !Ozone does, sunlit leaf (mmol O3/m^2)
   real(r8), allocatable :: o3uptakesha_p_   (:)  !Ozone does, shaded leaf (mmol O3/m^2)
   ! End Ozone Stress Variables

   ! for URBAN_MODEL
   real(r8), allocatable :: fwsun_           (:)  !sunlit fraction of walls [-]
   real(r8), allocatable :: dfwsun_          (:)  !change of sunlit fraction of walls [-]

   ! shortwave absorption
   real(r8), allocatable :: sroof_       (:,:,:)  !roof absorption [-]
   real(r8), allocatable :: swsun_       (:,:,:)  !sunlit wall absorption [-]
   real(r8), allocatable :: swsha_       (:,:,:)  !shaded wall absorption [-]
   real(r8), allocatable :: sgimp_       (:,:,:)  !impervious absorption [-]
   real(r8), allocatable :: sgper_       (:,:,:)  !pervious absorption [-]
   real(r8), allocatable :: slake_       (:,:,:)  !urban lake absorption [-]

   ! net longwave radiation for last time temperature change
   real(r8), allocatable :: lwsun_           (:)  !net longwave of sunlit wall [W/m2]
   real(r8), allocatable :: lwsha_           (:)  !net longwave of shaded wall [W/m2]
   real(r8), allocatable :: lgimp_           (:)  !net longwave of impervious  [W/m2]
   real(r8), allocatable :: lgper_           (:)  !net longwave of pervious [W/m2]
   real(r8), allocatable :: lveg_            (:)  !net longwave of vegetation [W/m2]

   real(r8), allocatable :: z_sno_roof_    (:,:)  !node depth of roof [m]
   real(r8), allocatable :: z_sno_gimp_    (:,:)  !node depth of impervious [m]
   real(r8), allocatable :: z_sno_gper_    (:,:)  !node depth pervious [m]
   real(r8), allocatable :: z_sno_lake_    (:,:)  !node depth lake [m]

   real(r8), allocatable :: dz_sno_roof_   (:,:)  !interface depth of roof [m]
   real(r8), allocatable :: dz_sno_gimp_   (:,:)  !interface depth of impervious [m]
   real(r8), allocatable :: dz_sno_gper_   (:,:)  !interface depth pervious [m]
   real(r8), allocatable :: dz_sno_lake_   (:,:)  !interface depth lake [m]

   real(r8), allocatable :: troof_inner_     (:)  !temperature of roof [K]
   real(r8), allocatable :: twsun_inner_     (:)  !temperature of sunlit wall [K]
   real(r8), allocatable :: twsha_inner_     (:)  !temperature of shaded wall [K]

   real(r8), allocatable :: t_roofsno_     (:,:)  !temperature of roof [K]
   real(r8), allocatable :: t_wallsun_     (:,:)  !temperature of sunlit wall [K]
   real(r8), allocatable :: t_wallsha_     (:,:)  !temperature of shaded wall [K]
   real(r8), allocatable :: t_gimpsno_     (:,:)  !temperature of impervious [K]
   real(r8), allocatable :: t_gpersno_     (:,:)  !temperature of pervious [K]
   real(r8), allocatable :: t_lakesno_     (:,:)  !temperature of pervious [K]

   real(r8), allocatable :: wliq_roofsno_  (:,:)  !liquid water in layers [kg/m2]
   real(r8), allocatable :: wliq_gimpsno_  (:,:)  !liquid water in layers [kg/m2]
   real(r8), allocatable :: wliq_gpersno_  (:,:)  !liquid water in layers [kg/m2]
   real(r8), allocatable :: wliq_lakesno_  (:,:)  !liquid water in layers [kg/m2]
   real(r8), allocatable :: wice_roofsno_  (:,:)  !ice lens in layers [kg/m2]
   real(r8), allocatable :: wice_gimpsno_  (:,:)  !ice lens in layers [kg/m2]
   real(r8), allocatable :: wice_gpersno_  (:,:)  !ice lens in layers [kg/m2]
   real(r8), allocatable :: wice_lakesno_  (:,:)  !ice lens in layers [kg/m2]

   real(r8), allocatable :: sag_roof_        (:)  !roof snow age [-]
   real(r8), allocatable :: sag_gimp_        (:)  !impervious ground snow age [-]
   real(r8), allocatable :: sag_gper_        (:)  !pervious ground snow age [-]
   real(r8), allocatable :: sag_lake_        (:)  !urban lake snow age [-]

   real(r8), allocatable :: scv_roof_        (:)  !roof snow cover [-]
   real(r8), allocatable :: scv_gimp_        (:)  !impervious ground snow cover [-]
   real(r8), allocatable :: scv_gper_        (:)  !pervious ground snow cover [-]
   real(r8), allocatable :: scv_lake_        (:)  !urban lake snow cover [-]

   real(r8), allocatable :: fsno_roof_       (:)  !roof snow fraction [-]
   real(r8), allocatable :: fsno_gimp_       (:)  !impervious ground snow fraction [-]
   real(r8), allocatable :: fsno_gper_       (:)  !pervious ground snow fraction [-]
   real(r8), allocatable :: fsno_lake_       (:)  !urban lake snow fraction [-]

   real(r8), allocatable :: snowdp_roof_     (:)  !roof snow depth [m]
   real(r8), allocatable :: snowdp_gimp_     (:)  !impervious ground snow depth [m]
   real(r8), allocatable :: snowdp_gper_     (:)  !pervious ground snow depth [m]
   real(r8), allocatable :: snowdp_lake_     (:)  !urban lake snow depth [m]

   !TODO: consider renaming the below variables
   real(r8), allocatable :: Fhac_            (:)  !sensible flux from heat or cool AC [W/m2]
   real(r8), allocatable :: Fwst_            (:)  !waste heat flux from heat or cool AC [W/m2]
   real(r8), allocatable :: Fach_            (:)  !flux from inner and outer air exchange [W/m2]
   real(r8), allocatable :: Fahe_            (:)  !flux from metabolism and vehicle [W/m2]
   real(r8), allocatable :: Fhah_            (:)  !sensible heat flux from heating [W/m2]
   real(r8), allocatable :: vehc_            (:)  !flux from vehicle [W/m2]
   real(r8), allocatable :: meta_            (:)  !flux from metabolism [W/m2]

   real(r8), allocatable :: t_room_          (:)  !temperature of inner building [K]
   real(r8), allocatable :: t_roof_          (:)  !temperature of roof [K]
   real(r8), allocatable :: t_wall_          (:)  !temperature of wall [K]
   real(r8), allocatable :: tafu_            (:)  !temperature of outer building [K]

   real(r8), allocatable :: urb_green_       (:)  !fractional of green leaf in urban patch [-]

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_LulccTimeVariables
   PUBLIC :: deallocate_LulccTimeVariables
   PUBLIC :: SAVE_LulccTimeVariables
   PUBLIC :: REST_LulccTimeVariables

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_LulccTimeVariables
   ! --------------------------------------------------------------------
   ! Allocates memory for Lulcc time variant variables
   ! --------------------------------------------------------------------

   USE MOD_SPMD_Task
   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_LandPatch
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_Vars_PFTimeVariables
   USE MOD_LandPFT
#endif
#ifdef URBAN_MODEL
   USE MOD_Urban_Vars_TimeVariables
   USE MOD_LandUrban
#endif

   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (z_sno_             (maxsnl+1:0,numpatch))
            allocate (dz_sno_            (maxsnl+1:0,numpatch))
            allocate (t_soisno_    (maxsnl+1:nl_soil,numpatch))
            allocate (wliq_soisno_ (maxsnl+1:nl_soil,numpatch))
            allocate (wice_soisno_ (maxsnl+1:nl_soil,numpatch))
            allocate (smp_                (1:nl_soil,numpatch))
            allocate (hk_                 (1:nl_soil,numpatch))
            allocate (t_grnd_                       (numpatch))
            allocate (tleaf_                        (numpatch))
            allocate (ldew_                         (numpatch))
            allocate (ldew_rain_                    (numpatch))
            allocate (ldew_snow_                    (numpatch))
            allocate (fwet_snow_                    (numpatch))
            allocate (sag_                          (numpatch))
            allocate (scv_                          (numpatch))
            allocate (snowdp_                       (numpatch))
            allocate (fsno_                         (numpatch))
            allocate (sigf_                         (numpatch))
            allocate (zwt_                          (numpatch))
            allocate (wa_                           (numpatch))
            allocate (wdsrf_                        (numpatch))
            allocate (rss_                          (numpatch))

            allocate (t_lake_               (nl_lake,numpatch))
            allocate (lake_icefrac_         (nl_lake,numpatch))
            allocate (savedtke1_                    (numpatch))

            !Plant Hydraulic variables
            allocate (vegwp_              (1:nvegwcs,numpatch))
            allocate (gs0sun_                       (numpatch))
            allocate (gs0sha_                       (numpatch))
            !END plant hydraulic variables

            !Ozone Stress variables
            allocate (lai_old_                      (numpatch))
            allocate (o3uptakesun_                  (numpatch))
            allocate (o3uptakesha_                  (numpatch))
            !End ozone stress variables

            allocate (snw_rds_          (maxsnl+1:0,numpatch))
            allocate (mss_bcpho_        (maxsnl+1:0,numpatch))
            allocate (mss_bcphi_        (maxsnl+1:0,numpatch))
            allocate (mss_ocpho_        (maxsnl+1:0,numpatch))
            allocate (mss_ocphi_        (maxsnl+1:0,numpatch))
            allocate (mss_dst1_         (maxsnl+1:0,numpatch))
            allocate (mss_dst2_         (maxsnl+1:0,numpatch))
            allocate (mss_dst3_         (maxsnl+1:0,numpatch))
            allocate (mss_dst4_         (maxsnl+1:0,numpatch))
            allocate (ssno_lyr_     (2,2,maxsnl+1:1,numpatch))

            allocate (trad_                        (numpatch))
            allocate (tref_                        (numpatch))
            allocate (qref_                        (numpatch))
            allocate (rst_                         (numpatch))
            allocate (emis_                        (numpatch))
            allocate (z0m_                         (numpatch))
            allocate (zol_                         (numpatch))
            allocate (rib_                         (numpatch))
            allocate (ustar_                       (numpatch))
            allocate (qstar_                       (numpatch))
            allocate (tstar_                       (numpatch))
            allocate (fm_                          (numpatch))
            allocate (fh_                          (numpatch))
            allocate (fq_                          (numpatch))

            allocate (sum_irrig_                   (numpatch))
            allocate (sum_irrig_count_             (numpatch))
         ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         IF (numpft > 0) THEN
            allocate (tleaf_p_                        (numpft))
            allocate (ldew_p_                         (numpft))
            allocate (ldew_rain_p_                    (numpft))
            allocate (ldew_snow_p_                    (numpft))
            allocate (fwet_snow_p_                    (numpft))
            allocate (sigf_p_                         (numpft))
            allocate (tref_p_                         (numpft))
            allocate (qref_p_                         (numpft))
            allocate (rst_p_                          (numpft))
            allocate (z0m_p_                          (numpft))

            ! Plant Hydraulic variables
            allocate (vegwp_p_              (1:nvegwcs,numpft))
            allocate (gs0sun_p_                       (numpft))
            allocate (gs0sha_p_                       (numpft))
            ! end plant hydraulic variables

            ! Allocate Ozone Stress Variables
            allocate (lai_old_p_                      (numpft))
            allocate (o3uptakesun_p_                  (numpft))
            allocate (o3uptakesha_p_                  (numpft))
            ! End allocate Ozone Stress Variables
         ENDIF
#endif

#ifdef URBAN_MODEL
         IF (numurban > 0) THEN
            allocate (fwsun_                        (numurban))
            allocate (dfwsun_                       (numurban))

            allocate (sroof_                    (2,2,numurban))
            allocate (swsun_                    (2,2,numurban))
            allocate (swsha_                    (2,2,numurban))
            allocate (sgimp_                    (2,2,numurban))
            allocate (sgper_                    (2,2,numurban))
            allocate (slake_                    (2,2,numurban))

            allocate (lwsun_                        (numurban))
            allocate (lwsha_                        (numurban))
            allocate (lgimp_                        (numurban))
            allocate (lgper_                        (numurban))
            allocate (lveg_                         (numurban))

            allocate (z_sno_roof_        (maxsnl+1:0,numurban))
            allocate (z_sno_gimp_        (maxsnl+1:0,numurban))
            allocate (z_sno_gper_        (maxsnl+1:0,numurban))
            allocate (z_sno_lake_        (maxsnl+1:0,numurban))

            allocate (dz_sno_roof_       (maxsnl+1:0,numurban))
            allocate (dz_sno_gimp_       (maxsnl+1:0,numurban))
            allocate (dz_sno_gper_       (maxsnl+1:0,numurban))
            allocate (dz_sno_lake_       (maxsnl+1:0,numurban))

            allocate (t_roofsno_   (maxsnl+1:nl_roof,numurban))
            allocate (t_wallsun_   (maxsnl+1:nl_wall,numurban))
            allocate (t_wallsha_   (maxsnl+1:nl_wall,numurban))
            allocate (t_gimpsno_   (maxsnl+1:nl_soil,numurban))
            allocate (t_gpersno_   (maxsnl+1:nl_soil,numurban))
            allocate (t_lakesno_   (maxsnl+1:nl_soil,numurban))

            allocate (troof_inner_                  (numurban))
            allocate (twsun_inner_                  (numurban))
            allocate (twsha_inner_                  (numurban))

            allocate (wliq_roofsno_(maxsnl+1:nl_roof,numurban))
            allocate (wice_roofsno_(maxsnl+1:nl_roof,numurban))
            allocate (wliq_gimpsno_(maxsnl+1:nl_soil,numurban))
            allocate (wice_gimpsno_(maxsnl+1:nl_soil,numurban))
            allocate (wliq_gpersno_(maxsnl+1:nl_soil,numurban))
            allocate (wice_gpersno_(maxsnl+1:nl_soil,numurban))
            allocate (wliq_lakesno_(maxsnl+1:nl_soil,numurban))
            allocate (wice_lakesno_(maxsnl+1:nl_soil,numurban))

            allocate (sag_roof_                     (numurban))
            allocate (sag_gimp_                     (numurban))
            allocate (sag_gper_                     (numurban))
            allocate (sag_lake_                     (numurban))
            allocate (scv_roof_                     (numurban))
            allocate (scv_gimp_                     (numurban))
            allocate (scv_gper_                     (numurban))
            allocate (scv_lake_                     (numurban))
            allocate (fsno_roof_                    (numurban))
            allocate (fsno_gimp_                    (numurban))
            allocate (fsno_gper_                    (numurban))
            allocate (fsno_lake_                    (numurban))
            allocate (snowdp_roof_                  (numurban))
            allocate (snowdp_gimp_                  (numurban))
            allocate (snowdp_gper_                  (numurban))
            allocate (snowdp_lake_                  (numurban))

            allocate (Fhac_                         (numurban))
            allocate (Fwst_                         (numurban))
            allocate (Fach_                         (numurban))
            allocate (Fahe_                         (numurban))
            allocate (Fhah_                         (numurban))
            allocate (vehc_                         (numurban))
            allocate (meta_                         (numurban))
            allocate (t_room_                       (numurban))
            allocate (t_roof_                       (numurban))
            allocate (t_wall_                       (numurban))
            allocate (tafu_                         (numurban))
            allocate (urb_green_                    (numurban))
         ENDIF
#endif
      ENDIF
   END SUBROUTINE allocate_LulccTimeVariables


   SUBROUTINE SAVE_LulccTimeVariables

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Vars_TimeVariables
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_Vars_PFTimeVariables
#endif
#ifdef URBAN_MODEL
   USE MOD_Urban_Vars_TimeVariables
#endif

   IMPLICIT NONE

      IF (p_is_worker) THEN
         z_sno_        = z_sno
         dz_sno_       = dz_sno
         t_soisno_     = t_soisno
         wliq_soisno_  = wliq_soisno
         wice_soisno_  = wice_soisno
         smp_          = smp
         hk_           = hk
         t_grnd_       = t_grnd
         tleaf_        = tleaf
         ldew_         = ldew
         ldew_rain_    = ldew_rain
         ldew_snow_    = ldew_snow
         fwet_snow_    = fwet_snow
         sag_          = sag
         scv_          = scv
         snowdp_       = snowdp
         fsno_         = fsno
         sigf_         = sigf
         zwt_          = zwt
         wa_           = wa
         wdsrf_        = wdsrf
         rss_          = rss

         t_lake_       = t_lake
         lake_icefrac_ = lake_icefrac
         savedtke1_    = savedtke1

IF(DEF_USE_PLANTHYDRAULICS)THEN
         vegwp_        = vegwp
         gs0sun_       = gs0sun
         gs0sha_       = gs0sha
ENDIF

IF(DEF_USE_OZONESTRESS)THEN
         lai_old_      = lai_old
         o3uptakesun_  = o3uptakesun
         o3uptakesha_  = o3uptakesha
ENDIF
         snw_rds_      = snw_rds
         mss_bcpho_    = mss_bcpho
         mss_bcphi_    = mss_bcphi
         mss_ocpho_    = mss_ocpho
         mss_ocphi_    = mss_ocphi
         mss_dst1_     = mss_dst1
         mss_dst2_     = mss_dst2
         mss_dst3_     = mss_dst3
         mss_dst4_     = mss_dst4
         ssno_lyr_     = ssno_lyr

         trad_         = trad
         tref_         = tref
         qref_         = qref
         rst_          = rst
         emis_         = emis
         z0m_          = z0m
         displa_       = displa
         zol_          = zol
         rib_          = rib
         ustar_        = ustar
         qstar_        = qstar
         tstar_        = tstar
         fm_           = fm
         fh_           = fh
         fq_           = fq

IF (DEF_USE_IRRIGATION) THEN
         sum_irrig_                    = sum_irrig
         sum_irrig_count_              = sum_irrig_count
ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         tleaf_p_      = tleaf_p
         ldew_p_       = ldew_p
         ldew_rain_p_  = ldew_rain_p
         ldew_snow_p_  = ldew_snow_p
         fwet_snow_p_  = fwet_snow_p
         sigf_p_       = sigf_p

         tref_p_       = tref_p
         qref_p_       = qref_p
         rst_p_        = rst_p
         z0m_p_        = z0m_p
IF(DEF_USE_PLANTHYDRAULICS)THEN
         ! Plant Hydraulic variables
         vegwp_p_      = vegwp_p
         gs0sun_p_     = gs0sun_p
         gs0sha_p_     = gs0sha_p
         ! end plant hydraulic variables
ENDIF
IF(DEF_USE_OZONESTRESS)THEN
         ! Ozone Stress Variables
         lai_old_p_      = lai_old_p
         o3uptakesun_p_  = o3uptakesun_p
         o3uptakesha_p_  = o3uptakesha_p
         ! End allocate Ozone Stress Variables
ENDIF
#endif

#ifdef URBAN_MODEL
         fwsun_        = fwsun
         dfwsun_       = dfwsun

         sroof_        = sroof
         swsun_        = swsun
         swsha_        = swsha
         sgimp_        = sgimp
         sgper_        = sgper
         slake_        = slake

         lwsun_        = lwsun
         lwsha_        = lwsha
         lgimp_        = lgimp
         lgper_        = lgper
         lveg_         = lveg

         z_sno_roof_   = z_sno_roof
         z_sno_gimp_   = z_sno_gimp
         z_sno_gper_   = z_sno_gper
         z_sno_lake_   = z_sno_lake

         dz_sno_roof_  = dz_sno_roof
         dz_sno_gimp_  = dz_sno_gimp
         dz_sno_gper_  = dz_sno_gper
         dz_sno_lake_  = dz_sno_lake

         t_roofsno_    = t_roofsno
         t_wallsun_    = t_wallsun
         t_wallsha_    = t_wallsha
         t_gimpsno_    = t_gimpsno
         t_gpersno_    = t_gpersno
         t_lakesno_    = t_lakesno

         troof_inner_  = troof_inner
         twsun_inner_  = twsun_inner
         twsha_inner_  = twsha_inner

         wliq_roofsno_ = wliq_roofsno
         wice_roofsno_ = wice_roofsno
         wliq_gimpsno_ = wliq_gimpsno
         wice_gimpsno_ = wice_gimpsno
         wliq_gpersno_ = wliq_gpersno
         wice_gpersno_ = wice_gpersno
         wliq_lakesno_ = wliq_lakesno
         wice_lakesno_ = wice_lakesno

         sag_roof_     = sag_roof
         sag_gimp_     = sag_gimp
         sag_gper_     = sag_gper
         sag_lake_     = sag_lake
         scv_roof_     = scv_roof
         scv_gimp_     = scv_gimp
         scv_gper_     = scv_gper
         scv_lake_     = scv_lake
         fsno_roof_    = fsno_roof
         fsno_gimp_    = fsno_gimp
         fsno_gper_    = fsno_gper
         fsno_lake_    = fsno_lake
         snowdp_roof_  = snowdp_roof
         snowdp_gimp_  = snowdp_gimp
         snowdp_gper_  = snowdp_gper
         snowdp_lake_  = snowdp_lake

         Fhac_         = Fhac
         Fwst_         = Fwst
         Fach_         = Fach
         Fahe_         = Fahe
         Fhah_         = Fhah
         vehc_         = vehc
         meta_         = meta
         t_room_       = t_room
         t_roof_       = t_roof
         t_wall_       = t_wall
         tafu_         = tafu
         urb_green_    = urb_green
#endif
      ENDIF

   END SUBROUTINE SAVE_LulccTimeVariables


   SUBROUTINE REST_LulccTimeVariables

   USE MOD_SPMD_Task
   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_LandElm
   USE MOD_Mesh
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables
   USE MOD_Lulcc_Vars_TimeInvariants
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_Vars_PFTimeInvariants
   USE MOD_Vars_PFTimeVariables
   USE MOD_LandPFT
#endif
#ifdef URBAN_MODEL
   USE MOD_Urban_Vars_TimeVariables
   USE MOD_LandUrban
#endif

   IMPLICIT NONE

   real(r8), allocatable, dimension(:) :: grid_patch_s , grid_patch_e
   real(r8), allocatable, dimension(:) :: grid_patch_s_, grid_patch_e_
   integer , allocatable, dimension(:) :: locpxl
   integer i, j, np, np_, ip, ip_, pc, pc_, u, u_
   integer ps, ps_, pe, pe_
   integer numpxl, ipxl

      IF (p_is_worker) THEN
         ! allocate with numelm
         allocate(grid_patch_s (numelm ))
         allocate(grid_patch_e (numelm ))
         allocate(grid_patch_s_(numelm_))
         allocate(grid_patch_e_(numelm_))

         grid_patch_e (:) = -1.
         grid_patch_s (:) = -1.
         grid_patch_e_(:) = -1.
         grid_patch_s_(:) = -1.

         ! loop for numelm of next year, patches at the beginning and end of
         ! the element were recorded landpatch%eindex is arranged in order,
         ! and the not land element is skipped so, IF element is missing, the
         ! recorder is -1.
         DO i=1, numelm
            ! how many patches in ith element in this worker
            numpxl = count(landpatch%eindex==landelm%eindex(i))

            IF (allocated(locpxl)) deallocate(locpxl)
            allocate(locpxl(numpxl))

            ! get all patches' index that eindex is equal the i element
            locpxl = pack([(ipxl, ipxl=1, numpatch)], &
                          landpatch%eindex==landelm%eindex(i))
            ! the min index is the start of patch's index
            grid_patch_s(i) = minval(locpxl)
            ! the max index is the end of patch's index
            grid_patch_e(i) = maxval(locpxl)
         ENDDO

         ! same as above, loop for numelm of previous year
         ! patches at the beginning and end of the element were recorded
         DO i=1, numelm_
            numpxl = count(landpatch_%eindex==landelm_%eindex(i))

            IF (allocated(locpxl)) deallocate(locpxl)
            allocate(locpxl(numpxl))

            locpxl = pack([(ipxl, ipxl=1, numpatch_)], &
                          landpatch_%eindex==landelm_%eindex(i))

            grid_patch_s_(i) = minval(locpxl)
            grid_patch_e_(i) = maxval(locpxl)
         ENDDO

         ! loop for element
         ! print*, 'minelm is', minelm, 'maxelm is', maxelm
         DO i=1, numelm
            DO j=1,numelm_
               IF (landelm%eindex(i) == landelm_%eindex(j)) THEN
                  np = grid_patch_s (i)
                  np_= grid_patch_s_(j)

                  IF (np.le.0 .or. np_.le.0) CYCLE

                  ! IF element is still present, loop for patches in same element
                  DO WHILE (np.le.grid_patch_e(i) .and. np_.le.grid_patch_e_(j))

                     ! IF a patch is missing, CYCLE
                     IF (patchclass(np) > patchclass_(np_)) THEN
                        np_= np_+ 1
                        CYCLE
                     ENDIF

                     ! IF a patch is added, CYCLE
                     IF (patchclass(np) < patchclass_(np_)) THEN
                        np = np + 1
                        CYCLE
                     ENDIF

#ifdef URBAN_MODEL
                     IF (numurban > 0) THEN
                        u = patch2urban (np )
                        u_= patch2urban_(np_)

                        ! vars assignment needs same urb class for urban patch
                        IF (patchclass(np) == URBAN) THEN
                           ! IF a Urban type is missing, CYCLE
                           IF (landurban%settyp(u) > urbclass_(u_)) THEN
                              np_= np_+ 1
                              CYCLE
                           ENDIF

                           ! IF a urban type is added, CYCLE
                           IF (landurban%settyp(u) < urbclass_(u_)) THEN
                              np = np + 1
                              CYCLE
                           ENDIF
                        ENDIF
                     ENDIF
#endif
                     ! otherwise, set patch value
                     ! only for the same patch type
                     z_sno       (:,np) = z_sno_       (:,np_)
                     dz_sno      (:,np) = dz_sno_      (:,np_)
                     t_soisno    (:,np) = t_soisno_    (:,np_)
                     wliq_soisno (:,np) = wliq_soisno_ (:,np_)
                     wice_soisno (:,np) = wice_soisno_ (:,np_)
                     scv           (np) = scv_           (np_)
                     smp         (:,np) = smp_         (:,np_)
                     hk          (:,np) = hk_          (:,np_)
                     t_grnd        (np) = t_grnd_        (np_)
                     tleaf         (np) = tleaf_         (np_)
                     ldew          (np) = ldew_          (np_)
                     ldew_rain     (np) = ldew_rain_     (np_)
                     ldew_snow     (np) = ldew_snow_     (np_)
                     fwet_snow     (np) = fwet_snow_     (np_)
                     sag           (np) = sag_           (np_)
                     snowdp        (np) = snowdp_        (np_)
                     fsno          (np) = fsno_          (np_)
                     sigf          (np) = sigf_          (np_)
                     ! In case lai+sai come into existence this year, set sigf to 1
                     IF ( (sigf(np) .eq. 0) .and. ((lai(np) + sai(np)) .gt. 0) ) THEN
                        sigf(np) = 1
                     ENDIF
                     zwt           (np) = zwt_           (np_)
                     wa            (np) = wa_            (np_)
                     wdsrf         (np) = wdsrf_         (np_)
                     rss           (np) = rss_           (np_)

                     t_lake      (:,np) = t_lake_      (:,np_)
                     lake_icefrac(:,np) = lake_icefrac_(:,np_)
                     savedtke1     (np) = savedtke1_     (np_)
IF(DEF_USE_PLANTHYDRAULICS)THEN
                     !Plant Hydraulic variables
                     vegwp       (:,np) = vegwp_       (:,np_)
                     gs0sun        (np) = gs0sun_        (np_)
                     gs0sha        (np) = gs0sha_        (np_)
                     !END plant hydraulic variables
ENDIF
IF(DEF_USE_OZONESTRESS)THEN
                     !Ozone Stress variables
                     lai_old       (np) = lai_old_       (np_)
                     o3uptakesun   (np) = o3uptakesun_   (np_)
                     o3uptakesha   (np) = o3uptakesha_   (np_)
                     !End ozone stress variables
ENDIF
                     snw_rds     (:,np) = snw_rds_     (:,np_)
                     mss_bcpho   (:,np) = mss_bcpho_   (:,np_)
                     mss_bcphi   (:,np) = mss_bcphi_   (:,np_)
                     mss_ocpho   (:,np) = mss_ocpho_   (:,np_)
                     mss_ocphi   (:,np) = mss_ocphi_   (:,np_)
                     mss_dst1    (:,np) = mss_dst1_    (:,np_)
                     mss_dst2    (:,np) = mss_dst2_    (:,np_)
                     mss_dst3    (:,np) = mss_dst3_    (:,np_)
                     mss_dst4    (:,np) = mss_dst4_    (:,np_)
                     ssno_lyr(2,2,:,np) = ssno_lyr_(2,2,:,np_)

                     trad          (np) = trad_          (np_)
                     tref          (np) = tref_          (np_)
                     qref          (np) = qref_          (np_)
                     rst           (np) = rst_           (np_)
                     emis          (np) = emis_          (np_)
                     z0m           (np) = z0m_           (np_)
                     zol           (np) = zol_           (np_)
                     rib           (np) = rib_           (np_)
                     ustar         (np) = ustar_         (np_)
                     qstar         (np) = qstar_         (np_)
                     tstar         (np) = tstar_         (np_)
                     fm            (np) = fm_            (np_)
                     fh            (np) = fh_            (np_)
                     fq            (np) = fq_            (np_)

IF(DEF_USE_IRRIGATION)THEN
                     sum_irrig       (np) = sum_irrig_       (np_)
                     sum_irrig_count (np) = sum_irrig_count_ (np_)
ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
IF (patchtype(np)==0 .and. patchtype_(np_)==0) THEN
                     ip = patch_pft_s (np )
                     ip_= patch_pft_s_(np_)

                     IF (ip.le.0 .or. ip_.le.0) THEN
                        print *, "Error in REST_LulccTimeVariables LULC_IGBP_PFT|LULC_IGBP_PC!"
                        CALL CoLM_stop ()
                     ENDIF

                     DO WHILE (ip.le.patch_pft_e(np) .and. ip_.le.patch_pft_e_(np_))

                        ! IF a PFT is missing, CYCLE
                        IF (pftclass(ip) > pftclass_(ip_)) THEN
                           ip_= ip_+ 1
                           CYCLE
                        ENDIF

                        ! IF a PFT is added, CYCLE
                        IF (pftclass(ip) < pftclass_(ip_)) THEN
                           ip = ip + 1
                           CYCLE
                        ENDIF

                        ! for the same PFT, set PFT value
                        tleaf_p    (ip) = tleaf_p_    (ip_)
                        ldew_p     (ip) = ldew_p_     (ip_)
                        ldew_rain_p(ip) = ldew_rain_p_(ip_)
                        ldew_snow_p(ip) = ldew_snow_p_(ip_)
                        fwet_snow_p(ip) = fwet_snow_p_(ip_)
                        sigf_p     (ip) = sigf_p_     (ip_)

                        tref_p     (ip) = tref_p_     (ip_)
                        qref_p     (ip) = qref_p_     (ip_)
                        rst_p      (ip) = rst_p_      (ip_)
                        z0m_p      (ip) = z0m_p_      (ip_)

IF(DEF_USE_PLANTHYDRAULICS)THEN
                        ! Plant Hydraulic variables
                        vegwp_p  (:,ip) = vegwp_p_  (:,ip_)
                        gs0sun_p   (ip) = gs0sun_p_   (ip_)
                        gs0sha_p   (ip) = gs0sha_p_   (ip_)
                        ! end plant hydraulic variables
ENDIF
IF(DEF_USE_OZONESTRESS)THEN
                        ! Ozone Stress Variables
                        lai_old_p     (ip) = lai_old_p_     (ip_)
                        o3uptakesun_p (ip) = o3uptakesun_p_ (ip_)
                        o3uptakesha_p (ip) = o3uptakesha_p_ (ip_)
                        ! End allocate Ozone Stress Variables
ENDIF
                        ip = ip + 1
                        ip_= ip_+ 1
                     ENDDO
                     ps       = patch_pft_s(np)
                     pe       = patch_pft_e(np)
                     ldew(np) = sum( ldew_p(ps:pe)*pftfrac(ps:pe) )
ENDIF
#endif

#ifdef URBAN_MODEL
IF (patchclass(np)==URBAN .and. patchclass_(np_)==URBAN) THEN

                     ! u = patch2urban (np )
                     ! u_= patch2urban_(np_)

                     IF (u.le.0 .or. u_.le.0) THEN
                        print *, "Error in REST_LulccTimeVariables URBAN_MODEL!"
                        CALL CoLM_stop ()
                     ENDIF

                     ! ! IF a Urban type is missing, CYCLE
                     ! IF (landurban%settyp(u) > urbclass_(u_)) THEN
                     !    np_= np_+ 1
                     !    CYCLE
                     ! ENDIF

                     ! ! IF a urban type is added, CYCLE
                     ! IF (landurban%settyp(u) < urbclass_(u_)) THEN
                     !    np = np + 1
                     !    CYCLE
                     ! ENDIF

                     ! otherwise, set urban value
                     ! include added urban and the same urban type
                     fwsun          (u) = fwsun_          (u_)
                     dfwsun         (u) = dfwsun_         (u_)

                     sroof      (:,:,u) = sroof_      (:,:,u_)
                     swsun      (:,:,u) = swsun_      (:,:,u_)
                     swsha      (:,:,u) = swsha_      (:,:,u_)
                     sgimp      (:,:,u) = sgimp_      (:,:,u_)
                     sgper      (:,:,u) = sgper_      (:,:,u_)
                     slake      (:,:,u) = slake_      (:,:,u_)

                     lwsun          (u) = lwsun_          (u_)
                     lwsha          (u) = lwsha_          (u_)
                     lgimp          (u) = lgimp_          (u_)
                     lgper          (u) = lgper_          (u_)
                     lveg           (u) = lveg_           (u_)

                     z_sno_roof   (:,u) = z_sno_roof_   (:,u_)
                     z_sno_gimp   (:,u) = z_sno_gimp_   (:,u_)
                     z_sno_gper   (:,u) = z_sno_gper_   (:,u_)
                     z_sno_lake   (:,u) = z_sno_lake_   (:,u_)

                     dz_sno_roof  (:,u) = dz_sno_roof_  (:,u_)
                     dz_sno_gimp  (:,u) = dz_sno_gimp_  (:,u_)
                     dz_sno_gper  (:,u) = dz_sno_gper_  (:,u_)
                     dz_sno_lake  (:,u) = dz_sno_lake_  (:,u_)

                     t_roofsno    (:,u) = t_roofsno_    (:,u_)
                     t_wallsun    (:,u) = t_wallsun_    (:,u_)
                     t_wallsha    (:,u) = t_wallsha_    (:,u_)
                     t_gimpsno    (:,u) = t_gimpsno_    (:,u_)
                     t_gpersno    (:,u) = t_gpersno_    (:,u_)
                     t_lakesno    (:,u) = t_lakesno_    (:,u_)

                     troof_inner    (u) = troof_inner_    (u_)
                     twsun_inner    (u) = twsun_inner_    (u_)
                     twsha_inner    (u) = twsha_inner_    (u_)

                     wliq_roofsno (:,u) = wliq_roofsno_ (:,u_)
                     wice_roofsno (:,u) = wice_roofsno_ (:,u_)
                     wliq_gimpsno (:,u) = wliq_gimpsno_ (:,u_)
                     wice_gimpsno (:,u) = wice_gimpsno_ (:,u_)
                     wliq_gpersno (:,u) = wliq_gpersno_ (:,u_)
                     wice_gpersno (:,u) = wice_gpersno_ (:,u_)
                     wliq_lakesno (:,u) = wliq_lakesno_ (:,u_)
                     wice_lakesno (:,u) = wice_lakesno_ (:,u_)

                     sag_roof       (u) = sag_roof_       (u_)
                     sag_gimp       (u) = sag_gimp_       (u_)
                     sag_gper       (u) = sag_gper_       (u_)
                     sag_lake       (u) = sag_lake_       (u_)
                     scv_roof       (u) = scv_roof_       (u_)
                     scv_gimp       (u) = scv_gimp_       (u_)
                     scv_gper       (u) = scv_gper_       (u_)
                     scv_lake       (u) = scv_lake_       (u_)
                     fsno_roof      (u) = fsno_roof_      (u_)
                     fsno_gimp      (u) = fsno_gimp_      (u_)
                     fsno_gper      (u) = fsno_gper_      (u_)
                     fsno_lake      (u) = fsno_lake_      (u_)
                     snowdp_roof    (u) = snowdp_roof_    (u_)
                     snowdp_gimp    (u) = snowdp_gimp_    (u_)
                     snowdp_gper    (u) = snowdp_gper_    (u_)
                     snowdp_lake    (u) = snowdp_lake_    (u_)

                     Fhac           (u) = Fhac_           (u_)
                     Fwst           (u) = Fwst_           (u_)
                     Fach           (u) = Fach_           (u_)
                     Fahe           (u) = Fahe_           (u_)
                     Fhah           (u) = Fhah_           (u_)
                     vehc           (u) = vehc_           (u_)
                     meta           (u) = meta_           (u_)
                     t_room         (u) = t_room_         (u_)
                     t_roof         (u) = t_roof_         (u_)
                     t_wall         (u) = t_wall_         (u_)
                     tafu           (u) = tafu_           (u_)
                     urb_green      (u) = urb_green_      (u_)

                     wliq_soisno(: ,np) = 0.
                     wliq_soisno(:1,np) = wliq_roofsno(:1,u )*froof(u)
                     wliq_soisno(: ,np) = wliq_soisno (: ,np) &
                                        + wliq_gpersno(: ,u )*(1-froof(u))*fgper(u)
                     wliq_soisno(:1,np) = wliq_soisno (:1,np) &
                                        + wliq_gimpsno(:1,u )*(1-froof(u))*(1-fgper(u))

                     wice_soisno(: ,np) = 0.
                     wice_soisno(:1,np) = wice_roofsno(:1,u )*froof(u)
                     wice_soisno(: ,np) = wice_soisno (: ,np) &
                                        + wice_gpersno(: ,u )*(1-froof(u))*fgper(u)
                     wice_soisno(:1,np) = wice_soisno (:1,np) &
                                        + wice_gimpsno(:1,u )*(1-froof(u))*(1-fgper(u))

                     scv(np) = scv_roof(u)*froof(u) + scv_gper(u)*(1-froof(u))*fgper(u) &
                             + scv_gimp(u)*(1-froof(u))*(1-fgper(u))
ENDIF
#endif
                     np = np + 1
                     np_= np_+ 1
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (p_is_worker) THEN
         IF (allocated(grid_patch_s )) deallocate (grid_patch_s  )
         IF (allocated(grid_patch_e )) deallocate (grid_patch_e  )
         IF (allocated(grid_patch_s_)) deallocate (grid_patch_s_ )
         IF (allocated(grid_patch_e_)) deallocate (grid_patch_e_ )
         IF (allocated(locpxl       )) deallocate (locpxl        )
      ENDIF
   END SUBROUTINE REST_LulccTimeVariables


   SUBROUTINE deallocate_LulccTimeVariables

   USE MOD_SPMD_Task
   USE MOD_Lulcc_Vars_TimeInvariants, only: numpatch_, numpft_, numpc_, numurban_

! --------------------------------------------------
! Deallocates memory for Lulcc time variant variables
! --------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch_ > 0) THEN
            deallocate (z_sno_           )
            deallocate (dz_sno_          )
            deallocate (t_soisno_        )
            deallocate (wliq_soisno_     )
            deallocate (wice_soisno_     )
            deallocate (smp_             )
            deallocate (hk_              )
            deallocate (t_grnd_          )
            deallocate (tleaf_           )
            deallocate (ldew_            )
            deallocate (ldew_rain_       )
            deallocate (ldew_snow_       )
            deallocate (fwet_snow_       )
            deallocate (sag_             )
            deallocate (scv_             )
            deallocate (snowdp_          )
            deallocate (fsno_            )
            deallocate (sigf_            )
            deallocate (zwt_             )
            deallocate (wa_              )
            deallocate (wdsrf_           )
            deallocate (rss_             )

            deallocate (t_lake_          )
            deallocate (lake_icefrac_    )
            deallocate (savedtke1_       )

            !Plant Hydraulic variables
            deallocate (vegwp_           )
            deallocate (gs0sun_          )
            deallocate (gs0sha_          )
            !END plant hydraulic variables

            !Ozone Stress variables
            deallocate (lai_old_         )
            deallocate (o3uptakesun_     )
            deallocate (o3uptakesha_     )
            !End ozone stress variables

            deallocate (snw_rds_         )
            deallocate (mss_bcpho_       )
            deallocate (mss_bcphi_       )
            deallocate (mss_ocpho_       )
            deallocate (mss_ocphi_       )
            deallocate (mss_dst1_        )
            deallocate (mss_dst2_        )
            deallocate (mss_dst3_        )
            deallocate (mss_dst4_        )
            deallocate (ssno_lyr_        )

            deallocate (trad_            )
            deallocate (tref_            )
            deallocate (qref_            )
            deallocate (rst_             )
            deallocate (emis_            )
            deallocate (z0m_             )
            deallocate (zol_             )
            deallocate (rib_             )
            deallocate (ustar_           )
            deallocate (qstar_           )
            deallocate (tstar_           )
            deallocate (fm_              )
            deallocate (fh_              )
            deallocate (fq_              )

            deallocate (sum_irrig_       )
            deallocate (sum_irrig_count_ )

         ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         IF (numpft_ > 0) THEN
            deallocate (tleaf_p_         )
            deallocate (ldew_p_          )
            deallocate (ldew_rain_p_     )
            deallocate (ldew_snow_p_     )
            deallocate (fwet_snow_p_     )
            deallocate (sigf_p_          )
            deallocate (tref_p_          )
            deallocate (qref_p_          )
            deallocate (rst_p_           )
            deallocate (z0m_p_           )

            ! Plant Hydraulic variables
            deallocate (vegwp_p_         )
            deallocate (gs0sun_p_        )
            deallocate (gs0sha_p_        )
            ! end plant hydraulic variables

            ! Allocate Ozone Stress Variables
            deallocate (lai_old_p_       )
            deallocate (o3uptakesun_p_   )
            deallocate (o3uptakesha_p_   )
            ! End allocate Ozone Stress Variables
         ENDIF
#endif

#ifdef URBAN_MODEL
         IF (numurban_ > 0) THEN
            deallocate (fwsun_           )
            deallocate (dfwsun_          )

            deallocate (sroof_           )
            deallocate (swsun_           )
            deallocate (swsha_           )
            deallocate (sgimp_           )
            deallocate (sgper_           )
            deallocate (slake_           )

            deallocate (lwsun_           )
            deallocate (lwsha_           )
            deallocate (lgimp_           )
            deallocate (lgper_           )
            deallocate (lveg_            )

            deallocate (z_sno_roof_      )
            deallocate (z_sno_gimp_      )
            deallocate (z_sno_gper_      )
            deallocate (z_sno_lake_      )

            deallocate (dz_sno_roof_     )
            deallocate (dz_sno_gimp_     )
            deallocate (dz_sno_gper_     )
            deallocate (dz_sno_lake_     )

            deallocate (t_roofsno_       )
            deallocate (t_wallsun_       )
            deallocate (t_wallsha_       )
            deallocate (t_gimpsno_       )
            deallocate (t_gpersno_       )
            deallocate (t_lakesno_       )

            deallocate (troof_inner_     )
            deallocate (twsun_inner_     )
            deallocate (twsha_inner_     )

            deallocate (wliq_roofsno_    )
            deallocate (wice_roofsno_    )
            deallocate (wliq_gimpsno_    )
            deallocate (wice_gimpsno_    )
            deallocate (wliq_gpersno_    )
            deallocate (wice_gpersno_    )
            deallocate (wliq_lakesno_    )
            deallocate (wice_lakesno_    )

            deallocate (sag_roof_        )
            deallocate (sag_gimp_        )
            deallocate (sag_gper_        )
            deallocate (sag_lake_        )
            deallocate (scv_roof_        )
            deallocate (scv_gimp_        )
            deallocate (scv_gper_        )
            deallocate (scv_lake_        )
            deallocate (fsno_roof_       )
            deallocate (fsno_gimp_       )
            deallocate (fsno_gper_       )
            deallocate (fsno_lake_       )
            deallocate (snowdp_roof_     )
            deallocate (snowdp_gimp_     )
            deallocate (snowdp_gper_     )
            deallocate (snowdp_lake_     )

            deallocate (Fhac_            )
            deallocate (Fwst_            )
            deallocate (Fach_            )
            deallocate (Fahe_            )
            deallocate (Fhah_            )
            deallocate (vehc_            )
            deallocate (meta_            )
            deallocate (t_room_          )
            deallocate (t_roof_          )
            deallocate (t_wall_          )
            deallocate (tafu_            )
            deallocate (urb_green_       )
         ENDIF
#endif
      ENDIF

   END SUBROUTINE deallocate_LulccTimeVariables

END MODULE MOD_Lulcc_Vars_TimeVariables
! ---------- EOP ------------
