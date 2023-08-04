#include <define.h>

MODULE MOD_LuLcc_Vars_TimeVariables
! -------------------------------
! Created by Hua Yuan, 04/2022
! -------------------------------

  USE MOD_Precision
  USE MOD_Vars_Global
  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run
  REAL(r8), allocatable :: z_sno_       (:,:)  !node depth [m]
  REAL(r8), allocatable :: dz_sno_      (:,:)  !interface depth [m]
  REAL(r8), allocatable :: t_soisno_    (:,:)  !soil temperature [K]
  REAL(r8), allocatable :: wliq_soisno_ (:,:)  !liquid water in layers [kg/m2]
  REAL(r8), allocatable :: wice_soisno_ (:,:)  !ice lens in layers [kg/m2]
  REAL(r8), allocatable :: t_grnd_        (:)  !ground surface temperature [K]

  REAL(r8), allocatable :: tleaf_         (:)  !leaf temperature [K]
  REAL(r8), allocatable :: ldew_          (:)  !depth of water on foliage [mm]
  REAL(r8), allocatable :: sag_           (:)  !non dimensional snow age [-]
  REAL(r8), allocatable :: scv_           (:)  !snow cover, water equivalent [mm]
  REAL(r8), allocatable :: snowdp_        (:)  !snow depth [meter]
  REAL(r8), allocatable :: fveg_          (:)  !fraction of vegetation cover
  REAL(r8), allocatable :: fsno_          (:)  !fraction of snow cover on ground
  REAL(r8), allocatable :: sigf_          (:)  !fraction of veg cover, excluding snow-covered veg [-]
  REAL(r8), allocatable :: green_         (:)  !leaf greenness
  REAL(r8), allocatable :: lai_           (:)  !leaf area index
  REAL(r8), allocatable :: sai_           (:)  !stem area index
  REAL(r8), allocatable :: coszen_        (:)  !cosine of solar zenith angle
  REAL(r8), allocatable :: alb_       (:,:,:)  !averaged albedo [-]
  REAL(r8), allocatable :: ssun_      (:,:,:)  !sunlit canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: ssha_      (:,:,:)  !shaded canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: thermk_        (:)  !canopy gap fraction for tir radiation
  REAL(r8), allocatable :: extkb_         (:)  !(k, g(mu)/mu) direct solar extinction coefficient
  REAL(r8), allocatable :: extkd_         (:)  !diffuse and scattered diffuse PAR extinction coefficient
  REAL(r8), allocatable :: zwt_           (:)  !the depth to water table [m]
  REAL(r8), allocatable :: wa_            (:)  !water storage in aquifer [mm]

  REAL(r8), allocatable :: t_lake_      (:,:)  !lake layer teperature [K]
  REAL(r8), allocatable :: lake_icefrac_(:,:)  !lake mass fraction of lake layer that is frozen

  ! for LULC_IGBP_PFT and LULC_IGBP_PC
  REAL(r8), allocatable :: tleaf_p_       (:)  !shaded leaf temperature [K]
  REAL(r8), allocatable :: ldew_p_        (:)  !depth of water on foliage [mm]
  REAL(r8), allocatable :: sigf_p_        (:)  !fraction of veg cover, excluding snow-covered veg [-]
  REAL(r8), allocatable :: lai_p_         (:)  !leaf area index
  REAL(r8), allocatable :: sai_p_         (:)  !stem area index
  REAL(r8), allocatable :: ssun_p_    (:,:,:)  !sunlit canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: ssha_p_    (:,:,:)  !shaded canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: thermk_p_      (:)  !canopy gap fraction for tir radiation
  REAL(r8), allocatable :: fshade_p_      (:)  !canopy shade fraction for tir radiation
  REAL(r8), allocatable :: extkb_p_       (:)  !(k, g(mu)/mu) direct solar extinction coefficient
  REAL(r8), allocatable :: extkd_p_       (:)  !diffuse and scattered diffuse PAR extinction coefficient

  ! for URBAN_MODEL
  REAL(r8), allocatable :: fwsun_         (:)  !sunlit fraction of walls [-]
  REAL(r8), allocatable :: dfwsun_        (:)  !change of sunlit fraction of walls [-]

  ! shortwave absorption
  REAL(r8), allocatable :: sroof_     (:,:,:)  !roof aborption [-]
  REAL(r8), allocatable :: swsun_     (:,:,:)  !sunlit wall absorption [-]
  REAL(r8), allocatable :: swsha_     (:,:,:)  !shaded wall absorption [-]
  REAL(r8), allocatable :: sgimp_     (:,:,:)  !impervious absorptioin [-]
  REAL(r8), allocatable :: sgper_     (:,:,:)  !pervious absorptioin [-]
  REAL(r8), allocatable :: slake_     (:,:,:)  !urban lake absorptioin [-]

  ! net longwave radiation for last time temperature change
  REAL(r8), allocatable :: lwsun_         (:)  !net longwave of sunlit wall [W/m2]
  REAL(r8), allocatable :: lwsha_         (:)  !net longwave of shaded wall [W/m2]
  REAL(r8), allocatable :: lgimp_         (:)  !net longwave of impervious  [W/m2]
  REAL(r8), allocatable :: lgper_         (:)  !net longwave of pervious [W/m2]
  REAL(r8), allocatable :: lveg_          (:)  !net longwave of vegetation [W/m2]

  REAL(r8), allocatable :: z_sno_roof_  (:,:)  !node depth of roof [m]
  REAL(r8), allocatable :: z_sno_gimp_  (:,:)  !node depth of impervious [m]
  REAL(r8), allocatable :: z_sno_gper_  (:,:)  !node depth pervious [m]
  REAL(r8), allocatable :: z_sno_lake_  (:,:)  !node depth lake [m]

  REAL(r8), allocatable :: dz_sno_roof_ (:,:)  !interface depth of roof [m]
  REAL(r8), allocatable :: dz_sno_gimp_ (:,:)  !interface depth of impervious [m]
  REAL(r8), allocatable :: dz_sno_gper_ (:,:)  !interface depth pervious [m]
  REAL(r8), allocatable :: dz_sno_lake_ (:,:)  !interface depth lake [m]

  REAL(r8), allocatable :: troof_inner_   (:)  !temperature of roof [K]
  REAL(r8), allocatable :: twsun_inner_   (:)  !temperature of sunlit wall [K]
  REAL(r8), allocatable :: twsha_inner_   (:)  !temperature of shaded wall [K]

  REAL(r8), allocatable :: t_roofsno_   (:,:)  !temperature of roof [K]
  REAL(r8), allocatable :: t_wallsun_   (:,:)  !temperature of sunlit wall [K]
  REAL(r8), allocatable :: t_wallsha_   (:,:)  !temperature of shaded wall [K]
  REAL(r8), allocatable :: t_gimpsno_   (:,:)  !temperature of impervious [K]
  REAL(r8), allocatable :: t_gpersno_   (:,:)  !temperature of pervious [K]
  REAL(r8), allocatable :: t_lakesno_   (:,:)  !temperature of pervious [K]

  REAL(r8), allocatable :: wliq_roofsno_(:,:)  !liquid water in layers [kg/m2]
  REAL(r8), allocatable :: wliq_gimpsno_(:,:)  !liquid water in layers [kg/m2]
  REAL(r8), allocatable :: wliq_gpersno_(:,:)  !liquid water in layers [kg/m2]
  REAL(r8), allocatable :: wliq_lakesno_(:,:)  !liquid water in layers [kg/m2]
  REAL(r8), allocatable :: wice_roofsno_(:,:)  !ice lens in layers [kg/m2]
  REAL(r8), allocatable :: wice_gimpsno_(:,:)  !ice lens in layers [kg/m2]
  REAL(r8), allocatable :: wice_gpersno_(:,:)  !ice lens in layers [kg/m2]
  REAL(r8), allocatable :: wice_lakesno_(:,:)  !ice lens in layers [kg/m2]

  REAL(r8), allocatable :: sag_roof_      (:)  !roof snow age [-]
  REAL(r8), allocatable :: sag_gimp_      (:)  !impervious ground snow age [-]
  REAL(r8), allocatable :: sag_gper_      (:)  !pervious ground snow age [-]
  REAL(r8), allocatable :: sag_lake_      (:)  !urban lake snow age [-]

  REAL(r8), allocatable :: scv_roof_      (:)  !roof snow cover [-]
  REAL(r8), allocatable :: scv_gimp_      (:)  !impervious ground snow cover [-]
  REAL(r8), allocatable :: scv_gper_      (:)  !pervious ground snow cover [-]
  REAL(r8), allocatable :: scv_lake_      (:)  !urban lake snow cover [-]

  REAL(r8), allocatable :: fsno_roof_     (:)  !roof snow fraction [-]
  REAL(r8), allocatable :: fsno_gimp_     (:)  !impervious ground snow fraction [-]
  REAL(r8), allocatable :: fsno_gper_     (:)  !pervious ground snow fraction [-]
  REAL(r8), allocatable :: fsno_lake_     (:)  !urban lake snow fraction [-]

  REAL(r8), allocatable :: snowdp_roof_   (:)  !roof snow depth [m]
  REAL(r8), allocatable :: snowdp_gimp_   (:)  !impervious ground snow depth [m]
  REAL(r8), allocatable :: snowdp_gper_   (:)  !pervious ground snow depth [m]
  REAL(r8), allocatable :: snowdp_lake_   (:)  !urban lake snow depth [m]

  REAL(r8), allocatable :: t_room_        (:)  !temperature of inner building [K]
  REAL(r8), allocatable :: tafu_          (:)  !temperature of outer building [K]
  REAL(r8), allocatable :: Fhac_          (:)  !sensible flux from heat or cool AC [W/m2]
  REAL(r8), allocatable :: Fwst_          (:)  !waste heat flux from heat or cool AC [W/m2]
  REAL(r8), allocatable :: Fach_          (:)  !flux from inner and outter air exchange [W/m2]


! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_LuLccTimeVariables
  PUBLIC :: deallocate_LuLccTimeVariables
  PUBLIC :: SAVE_LuLccTimeVariables
  PUBLIC :: REST_LuLccTimeVariables

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_LuLccTimeVariables
  ! --------------------------------------------------------------------
  ! Allocates memory for LuLcc time variant variables
  ! --------------------------------------------------------------------

     use MOD_SPMD_Task
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
           allocate (t_grnd_                       (numpatch))
           allocate (tleaf_                        (numpatch))
           allocate (ldew_                         (numpatch))
           allocate (sag_                          (numpatch))
           allocate (scv_                          (numpatch))
           allocate (snowdp_                       (numpatch))
           allocate (fveg_                         (numpatch))
           allocate (fsno_                         (numpatch))
           allocate (sigf_                         (numpatch))
           allocate (green_                        (numpatch))
           allocate (lai_                          (numpatch))
           allocate (sai_                          (numpatch))
           allocate (coszen_                       (numpatch))
           allocate (alb_                      (2,2,numpatch))
           allocate (ssun_                     (2,2,numpatch))
           allocate (ssha_                     (2,2,numpatch))
           allocate (thermk_                       (numpatch))
           allocate (extkb_                        (numpatch))
           allocate (extkd_                        (numpatch))
           allocate (zwt_                          (numpatch))
           allocate (wa_                           (numpatch))

           allocate (t_lake_               (nl_lake,numpatch))
           allocate (lake_icefrac_         (nl_lake,numpatch))
        ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
        IF (numpft > 0) THEN
           allocate (tleaf_p_                        (numpft))
           allocate (ldew_p_                         (numpft))
           allocate (sigf_p_                         (numpft))
           allocate (lai_p_                          (numpft))
           allocate (sai_p_                          (numpft))
           allocate (ssun_p_                     (2,2,numpft))
           allocate (ssha_p_                     (2,2,numpft))
           allocate (thermk_p_                       (numpft))
           allocate (fshade_p_                       (numpft))
           allocate (extkb_p_                        (numpft))
           allocate (extkd_p_                        (numpft))
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

           allocate (t_room_                       (numurban))
           allocate (tafu_                         (numurban))
           allocate (Fhac_                         (numurban))
           allocate (Fwst_                         (numurban))
           allocate (Fach_                         (numurban))
        ENDIF
#endif
     ENDIF
  END SUBROUTINE allocate_LuLccTimeVariables


  SUBROUTINE SAVE_LuLccTimeVariables

     USE MOD_Precision
     use MOD_SPMD_Task
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
         t_grnd_       = t_grnd
         tleaf_        = tleaf
         ldew_         = ldew
         sag_          = sag
         scv_          = scv
         snowdp_       = snowdp
         fveg_         = fveg
         fsno_         = fsno
         sigf_         = sigf
         green_        = green
         lai_          = lai
         sai_          = sai
         coszen_       = coszen
         alb_          = alb
         ssun_         = ssun
         ssha_         = ssha
         thermk_       = thermk
         extkb_        = extkb
         extkd_        = extkd
         zwt_          = zwt
         wa_           = wa

         t_lake_       = t_lake
         lake_icefrac_ = lake_icefrac

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         tleaf_p_      = tleaf_p
         ldew_p_       = ldew_p
         sigf_p_       = sigf_p
         lai_p_        = lai_p
         sai_p_        = sai_p
         ssun_p_       = ssun_p
         ssha_p_       = ssha_p
         thermk_p_     = thermk_p
         fshade_p_     = fshade_p
         extkb_p_      = extkb_p
         extkd_p_      = extkd_p
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

         t_room_       = t_room
         tafu_         = tafu
         Fhac_         = Fhac
         Fwst_         = Fwst
         Fach_         = Fach
#endif
     ENDIF

  END SUBROUTINE SAVE_LuLccTimeVariables


  SUBROUTINE REST_LuLccTimeVariables

     use MOD_SPMD_Task
     USE MOD_Precision
     USE MOD_Vars_Global
     USE MOD_LandPatch
     USE MOD_LandElm
     USE MOD_Mesh
     USE MOD_Vars_TimeInvariants
     USE MOD_Vars_TimeVariables
     USE MOD_LuLcc_Vars_TimeInvariants
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

     REAL(r8), allocatable, dimension(:) :: grid_patch_s , grid_patch_e
     REAL(r8), allocatable, dimension(:) :: grid_patch_s_, grid_patch_e_
     INTEGER , allocatable, dimension(:) :: locpxl
     INTEGER i, j, np, np_, ip, ip_, pc, pc_, u, u_
     INTEGER ps, ps_, pe, pe_
     INTEGER numpxl, ipxl

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
        ! and the not land element is skipped so, if element is missing, the
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

                 ! if element is still present, loop for patches in same element
                 DO WHILE (np.le.grid_patch_e(i) .and. np_.le.grid_patch_e_(j))

                    ! if a patch is missing, CYCLE
                    IF (patchclass(np) > patchclass_(np_)) THEN
                       np_= np_+ 1
                       CYCLE
                    ENDIF

                    ! if a patch is added, CYCLE
                    IF (patchclass(np) < patchclass_(np_)) THEN
                       np = np + 1
                       CYCLE
                    ENDIF

                    ! otherwise, set patch value
                    ! only for the same patch TYPE
                    z_sno       (:,np) = z_sno_       (:,np_)
                    dz_sno      (:,np) = dz_sno_      (:,np_)
                    t_soisno    (:,np) = t_soisno_    (:,np_)
                    wliq_soisno (:,np) = wliq_soisno_ (:,np_)
                    wice_soisno (:,np) = wice_soisno_ (:,np_)
                    t_grnd        (np) = t_grnd_        (np_)
                    tleaf         (np) = tleaf_         (np_)
                    ldew          (np) = ldew_          (np_)
                    sag           (np) = sag_           (np_)
                    scv           (np) = scv_           (np_)
                    snowdp        (np) = snowdp_        (np_)
                    fveg          (np) = fveg_          (np_)
                    fsno          (np) = fsno_          (np_)
                    sigf          (np) = sigf_          (np_)
                    green         (np) = green_         (np_)
                    lai           (np) = lai_           (np_)
                    sai           (np) = sai_           (np_)
                    coszen        (np) = coszen_        (np_)
                    alb       (:,:,np) = alb_       (:,:,np_)
                    ssun      (:,:,np) = ssun_      (:,:,np_)
                    ssha      (:,:,np) = ssha_      (:,:,np_)
                    thermk        (np) = thermk_        (np_)
                    extkb         (np) = extkb_         (np_)
                    extkd         (np) = extkd_         (np_)
                    zwt           (np) = zwt_           (np_)
                    wa            (np) = wa_            (np_)

                    t_lake      (:,np) = t_lake_      (:,np_)
                    lake_icefrac(:,np) = lake_icefrac_(:,np_)

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
IF (patchtype(np)==0 .and. patchtype_(np_)==0) THEN
                    ip = patch_pft_s (np )
                    ip_= patch_pft_s_(np_)

                    IF (ip.le.0 .or. ip_.le.0) THEN
                       print *, "Error in REST_LuLccTimeVariables LULC_IGBP_PFT|LULC_IGBP_PC!"
                       STOP
                    ENDIF

                    DO WHILE (ip.le.patch_pft_e(np) .and. ip_.le.patch_pft_e_(np_))

                       ! if a PFT is missing, CYCLE
                       IF (pftclass(ip) > pftclass_(ip_)) THEN
                          ip_= ip_+ 1
                          CYCLE
                       ENDIF

                       ! if a PFT is added, CYCLE
                       IF (pftclass(ip) < pftclass_(ip_)) THEN
                          ip = ip + 1
                          CYCLE
                       ENDIF

                       ! for the same PFT, set PFT value
                       tleaf_p    (ip) = tleaf_p_    (ip_)
                       ldew_p     (ip) = ldew_p_     (ip_)
                       sigf_p     (ip) = sigf_p_     (ip_)
                       lai_p      (ip) = lai_p_      (ip_)
                       sai_p      (ip) = sai_p_      (ip_)
                       ssun_p (:,:,ip) = ssun_p_ (:,:,ip_)
                       ssha_p (:,:,ip) = ssha_p_ (:,:,ip_)
                       thermk_p   (ip) = thermk_p_   (ip_)
                       fshade_p   (ip) = fshade_p_   (ip_)
                       extkb_p    (ip) = extkb_p_    (ip_)
                       extkd_p    (ip) = extkd_p_    (ip_)

                       ip = ip + 1
                       ip_= ip_+ 1
                    ENDDO
ENDIF
#endif

#ifdef URBAN_MODEL
IF (patchclass(np)==URBAN .and. patchclass_(np_)==URBAN) THEN
                    u = patch2urban (np )
                    u_= patch2urban_(np_)

                    IF (u.le.0 .or. u_.le.0) THEN
                       print *, "Error in REST_LuLccTimeVariables URBAN_MODEL!"
                       STOP
                    ENDIF

                    ! if a Urban TYPE is missing, CYCLE
                    IF (landurban%settyp(u) > urbclass_(u_)) THEN
                       np_= np_+ 1
                       CYCLE
                    ENDIF

                    ! if a urban TYPE is added, CYCLE
                    IF (landurban%settyp(u) < urbclass_(u_)) THEN
                       np = np + 1
                       CYCLE
                    ENDIF

                    ! otherwise, set urban value
                    ! include added urban and the same urban TYPE
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

                    t_room         (u) = t_room_         (u_)
                    tafu           (u) = tafu_           (u_)
                    Fhac           (u) = Fhac_           (u_)
                    Fwst           (u) = Fwst_           (u_)
                    Fach           (u) = Fach_           (u_)
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
        IF (allocated(grid_patch_s )) deallocate(grid_patch_s )
        IF (allocated(grid_patch_e )) deallocate(grid_patch_e )
        IF (allocated(grid_patch_s_)) deallocate(grid_patch_s_)
        IF (allocated(grid_patch_e_)) deallocate(grid_patch_e_)
        IF (allocated(locpxl       )) deallocate(locpxl       )
     ENDIF
  END SUBROUTINE REST_LuLccTimeVariables


  SUBROUTINE deallocate_LuLccTimeVariables
     use MOD_SPMD_Task
     USE MOD_LuLcc_Vars_TimeInvariants, only: numpatch_, numpft_, numpc_, numurban_

! --------------------------------------------------
! Deallocates memory for LuLcc time variant variables
! --------------------------------------------------
     IF (p_is_worker) THEN
        IF (numpatch_ > 0) THEN
           deallocate (z_sno_        )
           deallocate (dz_sno_       )
           deallocate (t_soisno_     )
           deallocate (wliq_soisno_  )
           deallocate (wice_soisno_  )
           deallocate (t_grnd_       )
           deallocate (tleaf_        )
           deallocate (ldew_         )
           deallocate (sag_          )
           deallocate (scv_          )
           deallocate (snowdp_       )
           deallocate (fveg_         )
           deallocate (fsno_         )
           deallocate (sigf_         )
           deallocate (green_        )
           deallocate (lai_          )
           deallocate (sai_          )
           deallocate (coszen_       )
           deallocate (alb_          )
           deallocate (ssun_         )
           deallocate (ssha_         )
           deallocate (thermk_       )
           deallocate (extkb_        )
           deallocate (extkd_        )
           deallocate (zwt_          )
           deallocate (wa_           )

           deallocate (t_lake_       )
           deallocate (lake_icefrac_ )
        ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
        IF (numpft_ > 0) THEN
           deallocate (tleaf_p_      )
           deallocate (ldew_p_       )
           deallocate (sigf_p_       )
           deallocate (lai_p_        )
           deallocate (sai_p_        )
           deallocate (ssun_p_       )
           deallocate (ssha_p_       )
           deallocate (thermk_p_     )
           deallocate (fshade_p_     )
           deallocate (extkb_p_      )
           deallocate (extkd_p_      )
        ENDIF
#endif

#ifdef URBAN_MODEL
        IF (numurban_ > 0) THEN
           deallocate (fwsun_        )
           deallocate (dfwsun_       )

           deallocate (sroof_        )
           deallocate (swsun_        )
           deallocate (swsha_        )
           deallocate (sgimp_        )
           deallocate (sgper_        )
           deallocate (slake_        )

           deallocate (lwsun_        )
           deallocate (lwsha_        )
           deallocate (lgimp_        )
           deallocate (lgper_        )
           deallocate (lveg_         )

           deallocate (z_sno_roof_   )
           deallocate (z_sno_gimp_   )
           deallocate (z_sno_gper_   )
           deallocate (z_sno_lake_   )

           deallocate (dz_sno_roof_  )
           deallocate (dz_sno_gimp_  )
           deallocate (dz_sno_gper_  )
           deallocate (dz_sno_lake_  )

           deallocate (t_roofsno_    )
           deallocate (t_wallsun_    )
           deallocate (t_wallsha_    )
           deallocate (t_gimpsno_    )
           deallocate (t_gpersno_    )
           deallocate (t_lakesno_    )

           deallocate (troof_inner_  )
           deallocate (twsun_inner_  )
           deallocate (twsha_inner_  )

           deallocate (wliq_roofsno_ )
           deallocate (wice_roofsno_ )
           deallocate (wliq_gimpsno_ )
           deallocate (wice_gimpsno_ )
           deallocate (wliq_gpersno_ )
           deallocate (wice_gpersno_ )
           deallocate (wliq_lakesno_ )
           deallocate (wice_lakesno_ )

           deallocate (sag_roof_     )
           deallocate (sag_gimp_     )
           deallocate (sag_gper_     )
           deallocate (sag_lake_     )
           deallocate (scv_roof_     )
           deallocate (scv_gimp_     )
           deallocate (scv_gper_     )
           deallocate (scv_lake_     )
           deallocate (fsno_roof_    )
           deallocate (fsno_gimp_    )
           deallocate (fsno_gper_    )
           deallocate (fsno_lake_    )
           deallocate (snowdp_roof_  )
           deallocate (snowdp_gimp_  )
           deallocate (snowdp_gper_  )
           deallocate (snowdp_lake_  )

           deallocate (t_room_       )
           deallocate (tafu_         )
           deallocate (Fhac_         )
           deallocate (Fwst_         )
           deallocate (Fach_         )
        ENDIF
#endif
     ENDIF

  END SUBROUTINE deallocate_LuLccTimeVariables

END MODULE MOD_LuLcc_Vars_TimeVariables
! ---------- EOP ------------
