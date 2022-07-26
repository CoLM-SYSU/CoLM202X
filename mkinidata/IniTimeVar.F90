#include <define.h>

SUBROUTINE IniTimeVar(ipatch, patchtype&
                     ,porsl,psi0,hksati,soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb&
                     ,z0m,zlnd,chil,rho,tau,z_soisno,dz_soisno&
                     ,t_soisno,wliq_soisno,wice_soisno,smp,hk,zwt,wa&
#ifdef PLANT_HYDRAULIC_STRESS
                     ,vegwp,gs0sun,gs0sha&
#endif
                     ,t_grnd,tleaf,ldew,ldew_rain,ldew_snow,sag,scv&   
                     ,snowdp,fveg,fsno,sigf,green,lai,sai,coszen&
                     ,alb,ssun,ssha,thermk,extkb,extkd&
                     ,trad,tref,qref,rst,emis,zol,rib&
                     ,ustar,qstar,tstar,fm,fh,fq&
#if(defined BGC)
                     ,totlitc, totsomc, totcwdc, decomp_cpools, decomp_cpools_vr, ctrunc_veg, ctrunc_soil, ctrunc_vr &
                     ,totlitn, totsomn, totcwdn, decomp_npools, decomp_npools_vr, ntrunc_veg, ntrunc_soil, ntrunc_vr &
                     ,totvegc, totvegn, totcolc, totcoln, col_endcb, col_begcb, col_endnb, col_begnb &
                     ,col_vegendcb, col_vegbegcb, col_soilendcb, col_soilbegcb &
                     ,col_vegendnb, col_vegbegnb, col_soilendnb, col_soilbegnb &
                     ,col_sminnendnb, col_sminnbegnb &
                     ,altmax, altmax_lastyear, altmax_lastyear_indx &
                     ,sminn_vr, sminn, smin_no3_vr, smin_nh4_vr &
                     ,prec10, prec60, prec365, prec_today, prec_daily, tsoi17, rh30, accumnstep, skip_balance_check &
#ifdef SASU
!------------------------SASU variables-----------------
                     ,decomp0_cpools_vr          , decomp0_npools_vr           &
                     ,I_met_c_vr_acc             , I_cel_c_vr_acc             , I_lig_c_vr_acc             , I_cwd_c_vr_acc              &
                     ,AKX_met_to_soil1_c_vr_acc  , AKX_cel_to_soil1_c_vr_acc  , AKX_lig_to_soil2_c_vr_acc  , AKX_soil1_to_soil2_c_vr_acc &
                     ,AKX_cwd_to_cel_c_vr_acc    , AKX_cwd_to_lig_c_vr_acc    , AKX_soil1_to_soil3_c_vr_acc, AKX_soil2_to_soil1_c_vr_acc &
                     ,AKX_soil2_to_soil3_c_vr_acc, AKX_soil3_to_soil1_c_vr_acc &
                     ,AKX_met_exit_c_vr_acc      , AKX_cel_exit_c_vr_acc      , AKX_lig_exit_c_vr_acc      , AKX_cwd_exit_c_vr_acc       &
                     ,AKX_soil1_exit_c_vr_acc    , AKX_soil2_exit_c_vr_acc    , AKX_soil3_exit_c_vr_acc     &
                     ,diagVX_c_vr_acc            , upperVX_c_vr_acc           , lowerVX_c_vr_acc            &
                     ,I_met_n_vr_acc             , I_cel_n_vr_acc             , I_lig_n_vr_acc             , I_cwd_n_vr_acc              &
                     ,AKX_met_to_soil1_n_vr_acc  , AKX_cel_to_soil1_n_vr_acc  , AKX_lig_to_soil2_n_vr_acc  , AKX_soil1_to_soil2_n_vr_acc &
                     ,AKX_cwd_to_cel_n_vr_acc    , AKX_cwd_to_lig_n_vr_acc    , AKX_soil1_to_soil3_n_vr_acc, AKX_soil2_to_soil1_n_vr_acc &
                     ,AKX_soil2_to_soil3_n_vr_acc, AKX_soil3_to_soil1_n_vr_acc &
                     ,AKX_met_exit_n_vr_acc      , AKX_cel_exit_n_vr_acc      , AKX_lig_exit_n_vr_acc      , AKX_cwd_exit_n_vr_acc       &
                     ,AKX_soil1_exit_n_vr_acc    , AKX_soil2_exit_n_vr_acc    , AKX_soil3_exit_n_vr_acc     &
                     ,diagVX_n_vr_acc            , upperVX_n_vr_acc           , lowerVX_n_vr_acc           &
#endif
!------------------------------------------------------------
#endif
#if(defined SOILINI)
                     ,nl_soil_ini,soil_z,soil_t,soil_w,snow_d)
#else
                     )
#endif

!=======================================================================
! Created by Yongjiu Dai, 09/15/1999
! Revised by Yongjiu Dai, 08/30/2002
!                         03/2014
!=======================================================================

  USE precision
  USE PhysicalConstants, only: tfrz
  USE MOD_TimeVariables, only: tlai, tsai
  USE PFT_Const, only: isevg, woody, leafcn, deadwdcn
#ifdef USE_DEPTH_TO_BEDROCK
  USE MOD_TimeInvariants, only : ibedrock, dbedrock
#endif
#if(defined PFT_CLASSIFICATION)
  USE mod_landpft, only : patch_pft_s, patch_pft_e
  USE MOD_PFTimeInvars
  USE MOD_PFTimeVars
#endif
#if(defined PC_CLASSIFICATION)
  USE mod_landpc
  USE MOD_PCTimeInvars
  USE MOD_PCTimeVars
#endif
  USE GlobalVars
  USE ALBEDO
#ifdef VARIABLY_SATURATED_FLOW
  USE MOD_TimeVariables, only: dpond
#endif

  IMPLICIT NONE 

  INTEGER, intent(in) ::        &! 
        ipatch,                 &! patch index
        patchtype                ! index for land cover TYPE [-]

  REAL(r8), intent(in) ::       &!
        fveg,                   &! fraction of vegetation cover
        green,                  &! leaf greenness
        coszen,                 &! cosine of solar zenith angle
        soil_s_v_alb,           &! albedo of visible of the saturated soil
        soil_d_v_alb,           &! albedo of visible of the dry soil
        soil_s_n_alb,           &! albedo of near infrared of the saturated soil
        soil_d_n_alb,           &! albedo of near infrared of the dry soil
        z0m,                    &! aerodynamic roughness length [m]
        zlnd,                   &! aerodynamic roughness length over soil surface [m] 
        chil,                   &! leaf angle distribution factor
        rho(2,2),               &! leaf reflectance (iw=iband, il=life and dead)
        tau(2,2),               &! leaf transmittance (iw=iband, il=life and dead)
        porsl(1:nl_soil),       &! porosity of soil
        psi0 (1:nl_soil),       &! saturated soil suction (mm) (NEGATIVE)
        hksati(1:nl_soil)        ! hydraulic conductivity at saturation [mm h2o/s]

#if(defined SOILINI)
  INTEGER, intent(in)  :: nl_soil_ini
  REAL(r8), intent(in) ::       &!
        soil_z(nl_soil_ini),    &! soil layer depth for initial (m)
        soil_t(nl_soil_ini),    &! soil temperature from initial file (K)
        soil_w(nl_soil_ini),    &! soil wetness from initial file (-)
        snow_d                   ! snow depth (m)
#endif

  REAL(r8), intent(inout) ::    &!
        z_soisno (maxsnl+1:nl_soil),   &! node depth [m]
        dz_soisno(maxsnl+1:nl_soil)     ! layer thickness [m]

  REAL(r8), intent(out) ::      &!
        t_soisno (maxsnl+1:nl_soil),   &! soil temperature [K]
        wliq_soisno(maxsnl+1:nl_soil), &! liquid water in layers [kg/m2]
        wice_soisno(maxsnl+1:nl_soil), &! ice lens in layers [kg/m2]
        smp        (1:nl_soil)       , &! soil matrix potential
        hk         (1:nl_soil)       , &! soil hydraulic conductance
#ifdef PLANT_HYDRAULIC_STRESS
        vegwp(1:nvegwcs),       &! vegetation water potential
        gs0sun,                 &! working copy of sunlit stomata conductance
        gs0sha,                 &! working copy of shalit stomata conductance
#endif
        t_grnd,                 &! ground surface temperature [K]
        tleaf,                  &! sunlit leaf temperature [K]
!#ifdef CLM5_INTERCEPTION
       ldew_rain,               &! depth of rain on foliage [mm]
       ldew_snow,               &! depth of snow on foliage [mm]   
!#endif
        ldew,                   &! depth of water on foliage [mm]
        sag,                    &! non dimensional snow age [-]
        scv,                    &! snow cover, water equivalent [mm]
        snowdp,                 &! snow depth [meter]
        fsno,                   &! fraction of snow cover on ground
        sigf,                   &! fraction of veg cover, excluding snow-covered veg [-]
        lai,                    &! leaf area index
        sai,                    &! stem area index

        alb (2,2),              &! averaged albedo [-]
        ssun(2,2),              &! sunlit canopy absorption for solar radiation
        ssha(2,2),              &! shaded canopy absorption for solar radiation
        thermk,                 &! canopy gap fraction for tir radiation
        extkb,                  &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,                  &! diffuse and scattered diffuse PAR extinction coefficient
        wa,                     &! water storage in aquifer [mm]
        zwt,                    &! the depth to water table [m]

                    ! Additional variables required by reginal model (WRF & RSM) 
                    ! ---------------------------------------------------------
        trad,                   &! radiative temperature of surface [K]
        tref,                   &! 2 m height air temperature [kelvin]
        qref,                   &! 2 m height air specific humidity
        rst,                    &! canopy stomatal resistance (s/m)
        emis,                   &! averaged bulk surface emissivity
        zol,                    &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,                    &! bulk Richardson number in surface layer
        ustar,                  &! u* in similarity theory [m/s]
        qstar,                  &! q* in similarity theory [kg/kg]
        tstar,                  &! t* in similarity theory [K]
        fm,                     &! integral of profile function for momentum
        fh,                     &! integral of profile function for heat
        fq                       ! integral of profile function for moisture

#ifdef BGC
   REAL(r8),intent(out) ::      &
        totlitc               , &
        totsomc               , &
        totcwdc               , &
        totvegc               , &
        totcolc               , &
        totlitn               , &
        totsomn               , &
        totcwdn               , &
        totvegn               , &
        totcoln               , &
        col_endcb             , &
        col_begcb             , &
        col_vegendcb          , &
        col_vegbegcb          , &
        col_soilendcb         , &
        col_soilbegcb         , &
        col_endnb             , &
        col_begnb             , &
        col_vegendnb          , &
        col_vegbegnb          , &
        col_soilendnb         , &
        col_soilbegnb         , &
        col_sminnendnb        , &
        col_sminnbegnb        , &
        decomp_cpools_vr          (nl_soil_full,ndecomp_pools), &
        decomp_cpools             (ndecomp_pools)             , &
        ctrunc_vr                 (nl_soil_full)              , &
        ctrunc_veg            , &    
        ctrunc_soil           , &    
        altmax                                                , &
        altmax_lastyear                                       
   INTEGER, intent(out) :: altmax_lastyear_indx     
   REAL(r8),intent(out) ::      &
        decomp_npools_vr          (nl_soil_full,ndecomp_pools), &
        decomp_npools             (ndecomp_pools)             , &
        ntrunc_vr                 (nl_soil_full)              , &
        ntrunc_veg            , &    
        ntrunc_soil           , &    
        sminn_vr                  (nl_soil)                   , &
        sminn                 , &
        smin_no3_vr               (nl_soil)                   , &
        smin_nh4_vr               (nl_soil)                   , &
        prec10                                                , &
        prec60                                                , &
        prec365                                               , &
        prec_today                                            , &
        prec_daily                                            , &
        tsoi17                                                , &
        rh30                                                  , &
        accumnstep                                            
    
#ifdef SASU
 !---------------SASU variables-----------------------
   REAL(r8),intent(out) ::      &
        decomp0_cpools_vr          (nl_soil,ndecomp_pools)  , &
        decomp0_npools_vr          (nl_soil,ndecomp_pools)  , &

        I_met_c_vr_acc             (nl_soil)  , &
        I_cel_c_vr_acc             (nl_soil)  , &
        I_lig_c_vr_acc             (nl_soil)  , &
        I_cwd_c_vr_acc             (nl_soil)  , &
        AKX_met_to_soil1_c_vr_acc  (nl_soil)  , &
        AKX_cel_to_soil1_c_vr_acc  (nl_soil)  , &
        AKX_lig_to_soil2_c_vr_acc  (nl_soil)  , &
        AKX_soil1_to_soil2_c_vr_acc(nl_soil)  , &
        AKX_cwd_to_cel_c_vr_acc    (nl_soil)  , &
        AKX_cwd_to_lig_c_vr_acc    (nl_soil)  , &
        AKX_soil1_to_soil3_c_vr_acc(nl_soil)  , &
        AKX_soil2_to_soil1_c_vr_acc(nl_soil)  , &
        AKX_soil2_to_soil3_c_vr_acc(nl_soil)  , &
        AKX_soil3_to_soil1_c_vr_acc(nl_soil)  , &
        AKX_met_exit_c_vr_acc      (nl_soil)  , &
        AKX_cel_exit_c_vr_acc      (nl_soil)  , &
        AKX_lig_exit_c_vr_acc      (nl_soil)  , &
        AKX_cwd_exit_c_vr_acc      (nl_soil)  , &
        AKX_soil1_exit_c_vr_acc    (nl_soil)  , &
        AKX_soil2_exit_c_vr_acc    (nl_soil)  , &
        AKX_soil3_exit_c_vr_acc    (nl_soil)  , &

        I_met_n_vr_acc             (nl_soil)  , &
        I_cel_n_vr_acc             (nl_soil)  , &
        I_lig_n_vr_acc             (nl_soil)  , &
        I_cwd_n_vr_acc             (nl_soil)  , &
        AKX_met_to_soil1_n_vr_acc  (nl_soil)  , &
        AKX_cel_to_soil1_n_vr_acc  (nl_soil)  , &
        AKX_lig_to_soil2_n_vr_acc  (nl_soil)  , &
        AKX_soil1_to_soil2_n_vr_acc(nl_soil)  , &
        AKX_cwd_to_cel_n_vr_acc    (nl_soil)  , &
        AKX_cwd_to_lig_n_vr_acc    (nl_soil)  , &
        AKX_soil1_to_soil3_n_vr_acc(nl_soil)  , &
        AKX_soil2_to_soil1_n_vr_acc(nl_soil)  , &
        AKX_soil2_to_soil3_n_vr_acc(nl_soil)  , &
        AKX_soil3_to_soil1_n_vr_acc(nl_soil)  , &
        AKX_met_exit_n_vr_acc      (nl_soil)  , &
        AKX_cel_exit_n_vr_acc      (nl_soil)  , &
        AKX_lig_exit_n_vr_acc      (nl_soil)  , &
        AKX_cwd_exit_n_vr_acc      (nl_soil)  , &
        AKX_soil1_exit_n_vr_acc    (nl_soil)  , &
        AKX_soil2_exit_n_vr_acc    (nl_soil)  , &
        AKX_soil3_exit_n_vr_acc    (nl_soil)  , &

        diagVX_c_vr_acc            (nl_soil,ndecomp_pools)  , &
        upperVX_c_vr_acc           (nl_soil,ndecomp_pools)  , &
        lowerVX_c_vr_acc           (nl_soil,ndecomp_pools)  , &
        diagVX_n_vr_acc            (nl_soil,ndecomp_pools)  , &
        upperVX_n_vr_acc           (nl_soil,ndecomp_pools)  , &
        lowerVX_n_vr_acc           (nl_soil,ndecomp_pools)  

#endif
 !----------------------------------------------------
   LOGICAL, intent(out) :: &
        skip_balance_check

#endif
        
        INTEGER j, snl, m, ivt                      
        REAL(r8) wet(nl_soil), wt, ssw, oro, rhosno_ini, a
        real(r8) alpha       (1:nl_soil) ! used in calculating hk
        real(r8) zmm         (1:nl_soil) ! z in mm
        real(r8) den         (1:nl_soil) !

        INTEGER ps, pe, pc
!-----------------------------------------------------------------------

  IF(patchtype <= 5)THEN ! land grid
     rhosno_ini = 250.
#if(defined SOILINI)
     DO j = 1, nl_soil
        CALL polint(soil_z,soil_t,nl_soil_ini,z_soisno(j),t_soisno(j))
        CALL polint(soil_z,soil_w,nl_soil_ini,z_soisno(j),wet(j))
        a = min(soil_t(1),soil_t(2),soil_t(3))-5.
        t_soisno(j) = max(t_soisno(j), a)
        a = max(soil_t(1),soil_t(2),soil_t(3))+5.
        t_soisno(j) = min(t_soisno(j), a)

        a = min(soil_w(1),soil_w(2),soil_w(3))
        wet(j) = max(wet(j), a, 0.1)
        a = max(soil_w(1),soil_w(2),soil_w(3))
        wet(j) = min(wet(j), a, 0.5)

        IF(t_soisno(j).ge.tfrz)THEN
           wliq_soisno(j) = wet(j)*dz_soisno(j)*1000.
!          wliq_soisno(j) = porsl(j)*wet(j)*dz_soisno(j)*1000.
           wice_soisno(j) = 0.
        ELSE
           wliq_soisno(j) = 0.
           wice_soisno(j) = wet(j)*dz_soisno(j)*1000.
!          wliq_soisno(j) = porsl(j)*wet(j)*dz_soisno(j)*1000.
        ENDIF
     ENDDO

     snowdp = snow_d
     sag    = 0.
     scv    = snowdp*rhosno_ini

! yuan, 08/02/2019: TODO: need to be changed in future
! for PFT_CLASSIFICATION or PC_CLASSIFICATION
! have done but not for SOILINI right now
     CALL snowfraction (lai,sai,z0m,zlnd,scv,snowdp,wt,sigf,fsno)
     CALL snow_ini (patchtype,maxsnl,snowdp,snl,z_soisno,dz_soisno)

     IF(snl.lt.0)THEN
        DO j = snl+1, 0
           t_soisno(j) = min(tfrz-1., t_soisno(1))
           wliq_soisno(j) = 0.
           wice_soisno(j) = dz_soisno(j)*rhosno_ini         !m * kg m-3 = kg m-2
        ENDDO
     ENDIF

     IF(snl>maxsnl)THEN
        t_soisno (maxsnl+1:snl) = -999.
        wice_soisno(maxsnl+1:snl) = 0.
        wliq_soisno(maxsnl+1:snl) = 0.
        z_soisno   (maxsnl+1:snl) = 0.
        dz_soisno  (maxsnl+1:snl) = 0.
     ENDIF
!#ifdef CLM5_INTERCEPTION
   ldew_rain  = 0.
   ldew_snow  = 0.
!#endif
     ldew  = 0.
     tleaf = t_soisno(1)
     t_grnd = t_soisno(1)

#ifdef PLANT_HYDRAULIC_STRESS
     vegwp(1:nvegwcs) = -2.5e4
     gs0sun = 1.0e4
     gs0sha = 1.0e4
#endif

#else
! soil temperature, water content and matrix potential
     DO j = 1, nl_soil
        IF(patchtype==3)THEN !land ice 
           t_soisno(j) = 253.
           wliq_soisno(j) = 0.
           wice_soisno(j) = dz_soisno(j)*1000.
           smp(j) = 1.e3 * 0.3336e6/9.80616*(t_soisno(j)-tfrz)/t_soisno(j)
        ELSE
           t_soisno(j) = 283.
           wliq_soisno(j) = dz_soisno(j)*porsl(j)*1000.
           wice_soisno(j) = 0.
           smp(j) = psi0(j)
        ENDIF
     ENDDO

! soil hydraulic conductivity
     zmm(1:) = z_soisno(1:)*1000.
     DO j = 1, nl_soil
        IF(patchtype==3)THEN !land ice 
           hk(j) = 0.
        ELSE
           if(j<nl_soil)then
              den(j)   = (zmm(j+1)-zmm(j))
              alpha(j) = (smp(j+1)-smp(j))/den(j) - 1._r8
           else
              alpha(j) = 0._r8
           end if
           if(alpha(j) <= 0.)then
              hk(j) = hksati(j)
           else
              hk(j) = hksati(j+1) 
           end if
        ENDIF
     ENDDO

#ifdef USE_DEPTH_TO_BEDROCK
     IF (patchtype <= 2) THEN
        IF (ibedrock(ipatch) <= nl_soil) THEN
           j = ibedrock(ipatch)
           IF (j == 1) THEN
              wliq_soisno(j) = dbedrock(ipatch) *porsl(j)*1000.
           else
              wliq_soisno(j) = (dbedrock(ipatch) - zi_soi(j-1)) *porsl(j)*1000.
           ENDIF 

           DO j = ibedrock(ipatch)+1, nl_soil
              wliq_soisno(j) = 0.
           ENDDO
        ENDIF 
     ENDIF
#endif


! water table depth (initially at 1.0 m below the model bottom; wa when zwt
!                    is below the model bottom zi(nl_soil)

     wa  = 4800.                             !assuming aquifer capacity is 5000 mm
     zwt = (25. + z_soisno(nl_soil))+dz_soisno(nl_soil)/2. - wa/1000./0.2 !to result in zwt = zi(nl_soil) + 1.0 m
#ifdef VARIABLY_SATURATED_FLOW
     wa = 0.
     zwt = zi_soi(nl_soil)
     dpond = 0.
#endif

! snow temperature and water content
     t_soisno(maxsnl+1:0) = -999.
     wice_soisno(maxsnl+1:0) = 0.
     wliq_soisno(maxsnl+1:0) = 0.
     z_soisno (maxsnl+1:0) = 0.
     dz_soisno(maxsnl+1:0) = 0.

IF (patchtype == 0) THEN

#if(defined USGS_CLASSIFICATION || defined IGBP_CLASSIFICATION)
     sigf   = fveg
!#ifdef CLM5_INTERCEPTION
     ldew_rain  = 0.
     ldew_snow  = 0.
!#endif
     ldew  = 0.
     tleaf  = t_soisno(1)
     lai    = tlai(ipatch)
     sai    = tsai(ipatch) * sigf
#ifdef PLANT_HYDRAULIC_STRESS
     vegwp(1:nvegwcs) = -2.5e4
     gs0sun = 1.0e4
     gs0sha = 1.0e4
#endif
#endif

#ifdef PFT_CLASSIFICATION
     ps = patch_pft_s(ipatch)      
     pe = patch_pft_e(ipatch)
     sigf_p(ps:pe)   = 1.
!#ifdef CLM5_INTERCEPTION
     ldew_p_rain(ps:pe)  = 0.
     ldew_p_snow(ps:pe)  = 0.
!#endif
     ldew_p(ps:pe)   = 0.
     tleaf_p(ps:pe)  = t_soisno(1)
#ifdef PLANT_HYDRAULIC_STRESS
     vegwp_p(1:nvegwcs,ps:pe) = -2.5e4
     gs0sun_p(ps:pe) = 1.0e4
     gs0sha_p(ps:pe) = 1.0e4
#endif
     lai_p(ps:pe)    = tlai_p(ps:pe)
     sai_p(ps:pe)    = tsai_p(ps:pe) * sigf_p(ps:pe)

     sigf  = 1.
!#ifdef CLM5_INTERCEPTION
     ldew_rain  = 0.
     ldew_snow  = 0.
!#endif
     ldew  = 0.
     tleaf = t_soisno(1)
#ifdef PLANT_HYDRAULIC_STRESS
     vegwp(1:nvegwcs) = -2.5e4
     gs0sun = 1.0e4
     gs0sha = 1.0e4
#endif
     lai   = tlai(ipatch)
     sai   = sum(sai_p(ps:pe) * pftfrac(ps:pe))
#endif

#ifdef PC_CLASSIFICATION
     pc = patch2pc(ipatch)
     sigf_c(:,pc)   = 1.
!#ifdef CLM5_INTERCEPTION
     ldew_c_rain(:,pc)  = 0.
     ldew_c_snow(:,pc)  = 0. 
!#endif
     ldew_c(:,pc)   = 0.
     tleaf_c(:,pc)  = t_soisno(1)
#ifdef PLANT_HYDRAULIC_STRESS
     vegwp_c(1:nvegwcs,:,pc) = -2.5e4
     gs0sun_c(:,pc) = 1.0e4
     gs0sha_c(:,pc) = 1.0e4
#endif
     lai_c(:,pc)    = tlai_c(:,pc)
     sai_c(:,pc)    = tsai_c(:,pc) * sigf_c(:,pc)

     sigf  = 1.
!#ifdef CLM5_INTERCEPTION
     ldew_rain  = 0.
     ldew_snow  = 0.
!#endif
     ldew  = 0.
     tleaf = t_soisno(1)
#ifdef PLANT_HYDRAULIC_STRESS
     vegwp(1:nvegwcs) = -2.5e4
     gs0sun = 1.0e4
     gs0sha = 1.0e4
#endif
     lai   = tlai(ipatch)
     sai   = sum(sai_c(:,pc)*pcfrac(:,pc))
#endif

ELSE
     sigf   = fveg
!#ifdef CLM5_INTERCEPTION
     ldew_rain  = 0.
     ldew_snow  = 0. 
!#endif
     ldew  = 0.
     tleaf  = t_soisno(1)
#ifdef PLANT_HYDRAULIC_STRESS
     vegwp(1:nvegwcs) = -2.5e4
     gs0sun = 1.0e4
     gs0sha = 1.0e4
#endif
     lai    = tlai(ipatch)
     sai    = tsai(ipatch) * sigf
ENDIF
     
     fsno   = 0.
     scv    = 0.
     sag    = 0.
     snowdp = 0.

     wt     = 0.
     t_grnd = t_soisno(1)
#endif

! surface albedo
     ssw = min(1.,1.e-3*wliq_soisno(1)/dz_soisno(1))
     CALL albland (ipatch,patchtype,soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb,&
                   chil,rho,tau,fveg,green,lai,sai,coszen,wt,fsno,scv,sag,ssw,t_grnd,&
                   alb,ssun,ssha,thermk,extkb,extkd)

  ELSE                 !ocean grid
     t_soisno(:) = 300.
     wice_soisno(:) = 0.
     wliq_soisno(:) = 1000.
     z_soisno (maxsnl+1:0) = 0.
     dz_soisno(maxsnl+1:0) = 0.
     sigf   = 0.
     fsno   = 0.
!#ifdef CLM5_INTERCEPTION
     ldew_rain  = 0.
     ldew_snow  = 0.
!#endif
     ldew  = 0.
     scv    = 0.
     sag    = 0.
     snowdp = 0.
     tleaf  = 300.
#ifdef PLANT_HYDRAULIC_STRESS
     vegwp(1:nvegwcs) = -2.5e4
     gs0sun = 1.0e4
     gs0sha = 1.0e4
#endif
     t_grnd = 300.

     oro = 0
     CALL albocean (oro,scv,coszen,alb)
     ssun(:,:) = 0.0
     ssha(:,:) = 0.0
     thermk = 0.0
     extkb = 0.0
     extkd = 0.0
  ENDIF

! Additional variables required by reginal model (WRF & RSM)
! totally arbitrarily assigned here
  trad  = t_grnd      
  tref  = t_grnd      
  qref  = 0.3     
  rst   = 1.e36   
  emis  = 1.0     
  zol   = -1.0    
  rib   = -0.1    
  ustar = 0.25    
  qstar = 0.001   
  tstar = -1.5    
  fm    = alog(30.)  
  fh    = alog(30.)  
  fq    = alog(30.)  

#ifdef BGC
    totlitc                         = 0.0
    totsomc                         = 0.0
    totcwdc                         = 0.0
    totvegc                         = 0.0
    totcolc                         = 0.0
    totlitn                         = 0.0
    totsomn                         = 0.0
    totcwdn                         = 0.0
    totvegn                         = 0.0
    totcoln                         = 0.0
    col_endcb                       = 0.0
    col_begcb                       = 0.0
    col_vegendcb                    = 0.0
    col_vegbegcb                    = 0.0
    col_soilendcb                   = 0.0
    col_soilbegcb                   = 0.0
    col_endnb                       = 0.0
    col_begnb                       = 0.0
    col_vegendnb                    = 0.0
    col_vegbegnb                    = 0.0
    col_soilendnb                   = 0.0
    col_soilbegnb                   = 0.0
    col_sminnendnb                  = 0.0
    col_sminnbegnb                  = 0.0
    decomp_cpools_vr          (:,:) = 0.0
    decomp_cpools             (:)   = 0.0
    ctrunc_vr                 (:)   = 0.0
    ctrunc_veg                      = 0.0
    ctrunc_soil                     = 0.0
    altmax                          = 10.0
    altmax_lastyear                 = 10.0
    altmax_lastyear_indx            = 10
    decomp_npools_vr          (:,:) = 0.0
    decomp_npools             (:)   = 0.0
    ntrunc_vr                 (:)   = 0.0
    ntrunc_veg                      = 0.0
    ntrunc_soil                     = 0.0
    sminn_vr                  (:)   = 0.0
    sminn                           = 0.0
    smin_no3_vr               (:)   = 0.0
    smin_nh4_vr               (:)   = 0.0
    sminn_vr                  (:)   = 0.0
    prec10                          = 0._r8
    prec60                          = 0._r8
    prec365                         = 0._r8
    prec_today                      = 0._r8
    prec_daily                      = 0._r8
    tsoi17                          = 273.15_r8
    rh30                            = 0._r8
    accumnstep                      = 0._r8
#ifdef SASU
 !---------------SASU variables-----------------------
    decomp0_cpools_vr         (:,:) = 0.0
    I_met_c_vr_acc              (:) = 0.0
    I_cel_c_vr_acc              (:) = 0.0
    I_lig_c_vr_acc              (:) = 0.0
    I_cwd_c_vr_acc              (:) = 0.0
    AKX_met_to_soil1_c_vr_acc   (:) = 0.0
    AKX_cel_to_soil1_c_vr_acc   (:) = 0.0
    AKX_lig_to_soil2_c_vr_acc   (:) = 0.0
    AKX_soil1_to_soil2_c_vr_acc (:) = 0.0
    AKX_cwd_to_cel_c_vr_acc     (:) = 0.0
    AKX_cwd_to_lig_c_vr_acc     (:) = 0.0
    AKX_soil1_to_soil3_c_vr_acc (:) = 0.0
    AKX_soil2_to_soil1_c_vr_acc (:) = 0.0
    AKX_soil2_to_soil3_c_vr_acc (:) = 0.0
    AKX_soil3_to_soil1_c_vr_acc (:) = 0.0
    AKX_met_exit_c_vr_acc       (:) = 0.0
    AKX_cel_exit_c_vr_acc       (:) = 0.0
    AKX_lig_exit_c_vr_acc       (:) = 0.0
    AKX_cwd_exit_c_vr_acc       (:) = 0.0
    AKX_soil1_exit_c_vr_acc     (:) = 0.0
    AKX_soil2_exit_c_vr_acc     (:) = 0.0
    AKX_soil3_exit_c_vr_acc     (:) = 0.0

    decomp0_npools_vr         (:,:) = 0.0
    I_met_n_vr_acc              (:) = 0.0
    I_cel_n_vr_acc              (:) = 0.0
    I_lig_n_vr_acc              (:) = 0.0
    I_cwd_n_vr_acc              (:) = 0.0
    AKX_met_to_soil1_n_vr_acc   (:) = 0.0
    AKX_cel_to_soil1_n_vr_acc   (:) = 0.0
    AKX_lig_to_soil2_n_vr_acc   (:) = 0.0
    AKX_soil1_to_soil2_n_vr_acc (:) = 0.0
    AKX_cwd_to_cel_n_vr_acc     (:) = 0.0
    AKX_cwd_to_lig_n_vr_acc     (:) = 0.0
    AKX_soil1_to_soil3_n_vr_acc (:) = 0.0
    AKX_soil2_to_soil1_n_vr_acc (:) = 0.0
    AKX_soil2_to_soil3_n_vr_acc (:) = 0.0
    AKX_soil3_to_soil1_n_vr_acc (:) = 0.0
    AKX_met_exit_n_vr_acc       (:) = 0.0
    AKX_cel_exit_n_vr_acc       (:) = 0.0
    AKX_lig_exit_n_vr_acc       (:) = 0.0
    AKX_cwd_exit_n_vr_acc       (:) = 0.0
    AKX_soil1_exit_n_vr_acc     (:) = 0.0
    AKX_soil2_exit_n_vr_acc     (:) = 0.0
    AKX_soil3_exit_n_vr_acc     (:) = 0.0

    diagVX_c_vr_acc           (:,:) = 0.0
    upperVX_c_vr_acc          (:,:) = 0.0
    lowerVX_c_vr_acc          (:,:) = 0.0
    diagVX_n_vr_acc           (:,:) = 0.0
    upperVX_n_vr_acc          (:,:) = 0.0
    lowerVX_n_vr_acc          (:,:) = 0.0

 !----------------------------------------------------
#endif
    skip_balance_check              = .false.

#if(defined PFT_CLASSIFICATION)
    IF (patchtype == 0) THEN
       do m = ps, pe
          ivt = pftclass(m)
          if(ivt .eq.  0)then  !no vegetation
             leafc_p                  (m) = 0.0
             leafc_storage_p          (m) = 0.0
             leafn_p                  (m) = 0.0
             leafn_storage_p          (m) = 0.0
          else
             if(isevg(ivt))then
                leafc_p               (m) = 100.0
                leafc_storage_p       (m) = 0.0
             else if(ivt >= npcropmin) then
                leafc_p               (m) = 0.0
                leafc_storage_p       (m) = 0.0
             else
                leafc_p               (m) = 0.0
                leafc_storage_p       (m) = 100.0
             end if
             leafn_p                  (m) = leafc_p        (m) / leafcn  (ivt)
             leafn_storage_p          (m) = leafc_storage_p(m) / leafcn  (ivt)
          end if
          if(woody(ivt) .eq. 1)then
             deadstemc_p              (m) = 0.1
             deadstemn_p              (m) = deadstemc_p    (m) / deadwdcn(ivt)
          else
             deadstemc_p              (m) = 0.0
             deadstemn_p              (m) = 0.0
          end if
          totcolc = totcolc + (leafc_p(m) + leafc_storage_p(m) + deadstemc_p(m))* pftfrac(m)
          totcoln = totcoln + (leafn_p(m) + leafn_storage_p(m) + deadstemn_p(m))* pftfrac(m)
       end do
       leafc_xfer_p             (ps:pe) = 0.0
       frootc_p                 (ps:pe) = 0.0
       frootc_storage_p         (ps:pe) = 0.0
       frootc_xfer_p            (ps:pe) = 0.0
       livestemc_p              (ps:pe) = 0.0
       livestemc_storage_p      (ps:pe) = 0.0
       livestemc_xfer_p         (ps:pe) = 0.0
       deadstemc_storage_p      (ps:pe) = 0.0
       deadstemc_xfer_p         (ps:pe) = 0.0
       livecrootc_p             (ps:pe) = 0.0
       livecrootc_storage_p     (ps:pe) = 0.0
       livecrootc_xfer_p        (ps:pe) = 0.0
       deadcrootc_p             (ps:pe) = 0.0
       deadcrootc_storage_p     (ps:pe) = 0.0
       deadcrootc_xfer_p        (ps:pe) = 0.0
       grainc_p                 (ps:pe) = 0.0
       grainc_storage_p         (ps:pe) = 0.0
       grainc_xfer_p            (ps:pe) = 0.0
       cropseedc_deficit_p      (ps:pe) = 0.0
       xsmrpool_p               (ps:pe) = 0.0
       gresp_storage_p          (ps:pe) = 0.0
       gresp_xfer_p             (ps:pe) = 0.0
       cpool_p                  (ps:pe) = 0.0
       cropprod1c_p             (ps:pe) = 0.0

       leafn_xfer_p             (ps:pe) = 0.0
       frootn_p                 (ps:pe) = 0.0
       frootn_storage_p         (ps:pe) = 0.0
       frootn_xfer_p            (ps:pe) = 0.0
       livestemn_p              (ps:pe) = 0.0
       livestemn_storage_p      (ps:pe) = 0.0
       livestemn_xfer_p         (ps:pe) = 0.0
       deadstemn_storage_p      (ps:pe) = 0.0
       deadstemn_xfer_p         (ps:pe) = 0.0
       livecrootn_p             (ps:pe) = 0.0
       livecrootn_storage_p     (ps:pe) = 0.0
       livecrootn_xfer_p        (ps:pe) = 0.0
       deadcrootn_p             (ps:pe) = 0.0
       deadcrootn_storage_p     (ps:pe) = 0.0
       deadcrootn_xfer_p        (ps:pe) = 0.0
       grainn_p                 (ps:pe) = 0.0
       grainn_storage_p         (ps:pe) = 0.0
       grainn_xfer_p            (ps:pe) = 0.0
       cropseedn_deficit_p      (ps:pe) = 0.0
       retransn_p               (ps:pe) = 0.0
       
       harvdate_p               (ps:pe) = 99999999
   
       tempsum_potential_gpp_p  (ps:pe) = 0.0
       tempmax_retransn_p       (ps:pe) = 0.0
       tempavg_tref_p           (ps:pe) = 0.0
       tempsum_npp_p            (ps:pe) = 0.0
       tempsum_litfall_p        (ps:pe) = 0.0
       annsum_potential_gpp_p   (ps:pe) = 0.0
       annmax_retransn_p        (ps:pe) = 0.0
       annavg_tref_p            (ps:pe) = 280.0
       annsum_npp_p             (ps:pe) = 0.0
       annsum_litfall_p         (ps:pe) = 0.0
       
       bglfr_p                  (ps:pe) = 0.0
       bgtr_p                   (ps:pe) = 0.0
       lgsf_p                   (ps:pe) = 0.0
       gdd0_p                   (ps:pe) = 0.0
       gdd8_p                   (ps:pe) = 0.0
       gdd10_p                  (ps:pe) = 0.0
       gdd020_p                 (ps:pe) = 0.0
       gdd820_p                 (ps:pe) = 0.0
       gdd1020_p                (ps:pe) = 0.0
       nyrs_crop_active_p       (ps:pe) = 0
       
       offset_flag_p            (ps:pe) = 0.0
       offset_counter_p         (ps:pe) = 0.0
       onset_flag_p             (ps:pe) = 0.0
       onset_counter_p          (ps:pe) = 0.0
       onset_gddflag_p          (ps:pe) = 0.0
       onset_gdd_p              (ps:pe) = 0.0
       onset_fdd_p              (ps:pe) = 0.0
       onset_swi_p              (ps:pe) = 0.0
       offset_fdd_p             (ps:pe) = 0.0
       offset_swi_p             (ps:pe) = 0.0
       dormant_flag_p           (ps:pe) = 1.0
       prev_leafc_to_litter_p   (ps:pe) = 0.0
       prev_frootc_to_litter_p  (ps:pe) = 0.0
       days_active_p            (ps:pe) = 0.0
    
       burndate_p               (ps:pe) = 10000
       grain_flag_p             (ps:pe) = 0.0

#ifdef CROP
! crop variables
       croplive_p               (ps:pe) = .false.
       gddtsoi_p                (ps:pe) =  spval
       huileaf_p                (ps:pe) =  spval
       gddplant_p               (ps:pe) =  spval
       huigrain_p               (ps:pe) =  0.0_r8
       peaklai_p                (ps:pe) =  0
       aroot_p                  (ps:pe) =  spval
       astem_p                  (ps:pe) =  spval
       arepr_p                  (ps:pe) =  spval
       aleaf_p                  (ps:pe) =  spval
       astemi_p                 (ps:pe) =  spval
       aleafi_p                 (ps:pe) =  spval
       gddmaturity_p            (ps:pe) =  spval

       cropplant_p              (ps:pe) = .false.
       idop_p                   (ps:pe) = 99999999
       a5tmin_p                 (ps:pe) = spval
       a10tmin_p                (ps:pe) = spval
       t10_p                    (ps:pe) = spval
       cumvd_p                  (ps:pe) = spval
       hdidx_p                  (ps:pe) = spval
       vf_p                     (ps:pe) = 0._r8
       cphase_p                 (ps:pe) = 4._r8
       fert_counter_p           (ps:pe) = 0._r8
       fert_p                   (ps:pe) = 0._r8
       tref_min_p               (ps:pe) = 273.15_r8
       tref_max_p               (ps:pe) = 273.15_r8
       tref_min_inst_p          (ps:pe) = spval
       tref_max_inst_p          (ps:pe) = spval
       fertnitro_p              (ps:pe) = spval
       latbaset_p               (ps:pe) = spval
#endif

#ifdef SASU
! SASU varaibles
       leafc0_p                 (ps:pe) = 0.0
       leafc0_storage_p         (ps:pe) = 0.0
       leafc0_xfer_p            (ps:pe) = 0.0
       frootc0_p                (ps:pe) = 0.0
       frootc0_storage_p        (ps:pe) = 0.0
       frootc0_xfer_p           (ps:pe) = 0.0
       livestemc0_p             (ps:pe) = 0.0
       livestemc0_storage_p     (ps:pe) = 0.0
       livestemc0_xfer_p        (ps:pe) = 0.0
       deadstemc0_p             (ps:pe) = 0.0
       deadstemc0_storage_p     (ps:pe) = 0.0
       deadstemc0_xfer_p        (ps:pe) = 0.0
       livecrootc0_p            (ps:pe) = 0.0
       livecrootc0_storage_p    (ps:pe) = 0.0
       livecrootc0_xfer_p       (ps:pe) = 0.0
       deadcrootc0_p            (ps:pe) = 0.0
       deadcrootc0_storage_p    (ps:pe) = 0.0
       deadcrootc0_xfer_p       (ps:pe) = 0.0
       grainc0_p                (ps:pe) = 0.0
       grainc0_storage_p        (ps:pe) = 0.0
       grainc0_xfer_p           (ps:pe) = 0.0

       leafn0_p                 (ps:pe) = 0.0
       leafn0_storage_p         (ps:pe) = 0.0
       leafn0_xfer_p            (ps:pe) = 0.0
       frootn0_p                (ps:pe) = 0.0
       frootn0_storage_p        (ps:pe) = 0.0
       frootn0_xfer_p           (ps:pe) = 0.0
       livestemn0_p             (ps:pe) = 0.0
       livestemn0_storage_p     (ps:pe) = 0.0
       livestemn0_xfer_p        (ps:pe) = 0.0
       deadstemn0_p             (ps:pe) = 0.0
       deadstemn0_storage_p     (ps:pe) = 0.0
       deadstemn0_xfer_p        (ps:pe) = 0.0
       livecrootn0_p            (ps:pe) = 0.0
       livecrootn0_storage_p    (ps:pe) = 0.0
       livecrootn0_xfer_p       (ps:pe) = 0.0
       deadcrootn0_p            (ps:pe) = 0.0
       deadcrootn0_storage_p    (ps:pe) = 0.0
       deadcrootn0_xfer_p       (ps:pe) = 0.0
       grainn0_p                (ps:pe) = 0.0
       grainn0_storage_p        (ps:pe) = 0.0
       grainn0_xfer_p           (ps:pe) = 0.0
       retransn0_p              (ps:pe) = 0.0

       I_leafc_p_acc            (ps:pe) = 0._r8
       I_leafc_st_p_acc         (ps:pe) = 0._r8
       I_frootc_p_acc           (ps:pe) = 0._r8
       I_frootc_st_p_acc        (ps:pe) = 0._r8
       I_livestemc_p_acc        (ps:pe) = 0._r8
       I_livestemc_st_p_acc     (ps:pe) = 0._r8
       I_deadstemc_p_acc        (ps:pe) = 0._r8
       I_deadstemc_st_p_acc     (ps:pe) = 0._r8
       I_livecrootc_p_acc       (ps:pe) = 0._r8
       I_livecrootc_st_p_acc    (ps:pe) = 0._r8
       I_deadcrootc_p_acc       (ps:pe) = 0._r8
       I_deadcrootc_st_p_acc    (ps:pe) = 0._r8
       I_grainc_p_acc           (ps:pe) = 0._r8
       I_grainc_st_p_acc        (ps:pe) = 0._r8
       I_leafn_p_acc            (ps:pe) = 0._r8
       I_leafn_st_p_acc         (ps:pe) = 0._r8
       I_frootn_p_acc           (ps:pe) = 0._r8
       I_frootn_st_p_acc        (ps:pe) = 0._r8
       I_livestemn_p_acc        (ps:pe) = 0._r8
       I_livestemn_st_p_acc     (ps:pe) = 0._r8
       I_deadstemn_p_acc        (ps:pe) = 0._r8
       I_deadstemn_st_p_acc     (ps:pe) = 0._r8
       I_livecrootn_p_acc       (ps:pe) = 0._r8
       I_livecrootn_st_p_acc    (ps:pe) = 0._r8
       I_deadcrootn_p_acc       (ps:pe) = 0._r8
       I_deadcrootn_st_p_acc    (ps:pe) = 0._r8
       I_grainn_p_acc           (ps:pe) = 0._r8
       I_grainn_st_p_acc        (ps:pe) = 0._r8

       AKX_leafc_xf_to_leafc_p_acc                 (ps:pe) = 0._r8
       AKX_frootc_xf_to_frootc_p_acc               (ps:pe) = 0._r8
       AKX_livestemc_xf_to_livestemc_p_acc         (ps:pe) = 0._r8
       AKX_deadstemc_xf_to_deadstemc_p_acc         (ps:pe) = 0._r8
       AKX_livecrootc_xf_to_livecrootc_p_acc       (ps:pe) = 0._r8
       AKX_deadcrootc_xf_to_deadcrootc_p_acc       (ps:pe) = 0._r8
       AKX_grainc_xf_to_grainc_p_acc               (ps:pe) = 0._r8
       AKX_livestemc_to_deadstemc_p_acc            (ps:pe) = 0._r8
       AKX_livecrootc_to_deadcrootc_p_acc          (ps:pe) = 0._r8
      
       AKX_leafc_st_to_leafc_xf_p_acc              (ps:pe) = 0._r8
       AKX_frootc_st_to_frootc_xf_p_acc            (ps:pe) = 0._r8
       AKX_livestemc_st_to_livestemc_xf_p_acc      (ps:pe) = 0._r8
       AKX_deadstemc_st_to_deadstemc_xf_p_acc      (ps:pe) = 0._r8
       AKX_livecrootc_st_to_livecrootc_xf_p_acc    (ps:pe) = 0._r8
       AKX_deadcrootc_st_to_deadcrootc_xf_p_acc    (ps:pe) = 0._r8
       AKX_grainc_st_to_grainc_xf_p_acc            (ps:pe) = 0._r8

       AKX_leafc_exit_p_acc                        (ps:pe) = 0._r8
       AKX_frootc_exit_p_acc                       (ps:pe) = 0._r8
       AKX_livestemc_exit_p_acc                    (ps:pe) = 0._r8
       AKX_deadstemc_exit_p_acc                    (ps:pe) = 0._r8
       AKX_livecrootc_exit_p_acc                   (ps:pe) = 0._r8
       AKX_deadcrootc_exit_p_acc                   (ps:pe) = 0._r8
       AKX_grainc_exit_p_acc                       (ps:pe) = 0._r8

       AKX_leafc_st_exit_p_acc                     (ps:pe) = 0._r8
       AKX_frootc_st_exit_p_acc                    (ps:pe) = 0._r8
       AKX_livestemc_st_exit_p_acc                 (ps:pe) = 0._r8
       AKX_deadstemc_st_exit_p_acc                 (ps:pe) = 0._r8
       AKX_livecrootc_st_exit_p_acc                (ps:pe) = 0._r8
       AKX_deadcrootc_st_exit_p_acc                (ps:pe) = 0._r8
       AKX_grainc_st_exit_p_acc                    (ps:pe) = 0._r8

       AKX_leafc_xf_exit_p_acc                     (ps:pe) = 0._r8
       AKX_frootc_xf_exit_p_acc                    (ps:pe) = 0._r8
       AKX_livestemc_xf_exit_p_acc                 (ps:pe) = 0._r8
       AKX_deadstemc_xf_exit_p_acc                 (ps:pe) = 0._r8
       AKX_livecrootc_xf_exit_p_acc                (ps:pe) = 0._r8
       AKX_deadcrootc_xf_exit_p_acc                (ps:pe) = 0._r8
       AKX_grainc_xf_exit_p_acc                    (ps:pe) = 0._r8
      
       AKX_leafn_xf_to_leafn_p_acc                 (ps:pe) = 0._r8        
       AKX_frootn_xf_to_frootn_p_acc               (ps:pe) = 0._r8
       AKX_livestemn_xf_to_livestemn_p_acc         (ps:pe) = 0._r8
       AKX_deadstemn_xf_to_deadstemn_p_acc         (ps:pe) = 0._r8
       AKX_livecrootn_xf_to_livecrootn_p_acc       (ps:pe) = 0._r8
       AKX_deadcrootn_xf_to_deadcrootn_p_acc       (ps:pe) = 0._r8
       AKX_grainn_xf_to_grainn_p_acc               (ps:pe) = 0._r8
       AKX_livestemn_to_deadstemn_p_acc            (ps:pe) = 0._r8
       AKX_livecrootn_to_deadcrootn_p_acc          (ps:pe) = 0._r8

       AKX_leafn_st_to_leafn_xf_p_acc              (ps:pe) = 0._r8
       AKX_frootn_st_to_frootn_xf_p_acc            (ps:pe) = 0._r8
       AKX_livestemn_st_to_livestemn_xf_p_acc      (ps:pe) = 0._r8
       AKX_deadstemn_st_to_deadstemn_xf_p_acc      (ps:pe) = 0._r8
       AKX_livecrootn_st_to_livecrootn_xf_p_acc    (ps:pe) = 0._r8
       AKX_deadcrootn_st_to_deadcrootn_xf_p_acc    (ps:pe) = 0._r8
       AKX_grainn_st_to_grainn_xf_p_acc            (ps:pe) = 0._r8

       AKX_leafn_to_retransn_p_acc                 (ps:pe) = 0._r8
       AKX_frootn_to_retransn_p_acc                (ps:pe) = 0._r8
       AKX_livestemn_to_retransn_p_acc             (ps:pe) = 0._r8
       AKX_livecrootn_to_retransn_p_acc            (ps:pe) = 0._r8

       AKX_retransn_to_leafn_p_acc                 (ps:pe) = 0._r8
       AKX_retransn_to_frootn_p_acc                (ps:pe) = 0._r8
       AKX_retransn_to_livestemn_p_acc             (ps:pe) = 0._r8
       AKX_retransn_to_deadstemn_p_acc             (ps:pe) = 0._r8
       AKX_retransn_to_livecrootn_p_acc            (ps:pe) = 0._r8
       AKX_retransn_to_deadcrootn_p_acc            (ps:pe) = 0._r8
       AKX_retransn_to_grainn_p_acc                (ps:pe) = 0._r8

       AKX_retransn_to_leafn_st_p_acc              (ps:pe) = 0._r8
       AKX_retransn_to_frootn_st_p_acc             (ps:pe) = 0._r8
       AKX_retransn_to_livestemn_st_p_acc          (ps:pe) = 0._r8
       AKX_retransn_to_deadstemn_st_p_acc          (ps:pe) = 0._r8
       AKX_retransn_to_livecrootn_st_p_acc         (ps:pe) = 0._r8
       AKX_retransn_to_deadcrootn_st_p_acc         (ps:pe) = 0._r8
       AKX_retransn_to_grainn_st_p_acc             (ps:pe) = 0._r8

       AKX_leafn_exit_p_acc                        (ps:pe) = 0._r8
       AKX_frootn_exit_p_acc                       (ps:pe) = 0._r8
       AKX_livestemn_exit_p_acc                    (ps:pe) = 0._r8
       AKX_deadstemn_exit_p_acc                    (ps:pe) = 0._r8
       AKX_livecrootn_exit_p_acc                   (ps:pe) = 0._r8
       AKX_deadcrootn_exit_p_acc                   (ps:pe) = 0._r8
       AKX_grainn_exit_p_acc                       (ps:pe) = 0._r8
       AKX_retransn_exit_p_acc                     (ps:pe) = 0._r8

       AKX_leafn_st_exit_p_acc                     (ps:pe) = 0._r8
       AKX_frootn_st_exit_p_acc                    (ps:pe) = 0._r8
       AKX_livestemn_st_exit_p_acc                 (ps:pe) = 0._r8
       AKX_deadstemn_st_exit_p_acc                 (ps:pe) = 0._r8
       AKX_livecrootn_st_exit_p_acc                (ps:pe) = 0._r8
       AKX_deadcrootn_st_exit_p_acc                (ps:pe) = 0._r8
       AKX_grainn_st_exit_p_acc                    (ps:pe) = 0._r8

       AKX_leafn_xf_exit_p_acc                     (ps:pe) = 0._r8
       AKX_frootn_xf_exit_p_acc                    (ps:pe) = 0._r8
       AKX_livestemn_xf_exit_p_acc                 (ps:pe) = 0._r8
       AKX_deadstemn_xf_exit_p_acc                 (ps:pe) = 0._r8
       AKX_livecrootn_xf_exit_p_acc                (ps:pe) = 0._r8
       AKX_deadcrootn_xf_exit_p_acc                (ps:pe) = 0._r8
       AKX_grainn_xf_exit_p_acc                    (ps:pe) = 0._r8
#endif
       
    end if
#endif
#endif
END SUBROUTINE IniTimeVar
!-----------------------------------------------------------------------
! EOP


SUBROUTINE snow_ini(patchtype,maxsnl,snowdp,snl,z_soisno,dz_soisno)

! Snow spatial discretization initially

  USE precision
  IMPLICIT NONE

  INTEGER,  intent(in) :: maxsnl    !maximum of snow layers
  INTEGER,  intent(in) :: patchtype !index for land cover TYPE [-]
  REAL(r8), intent(in) :: snowdp    !snow depth [m]
  REAL(r8), intent(out) :: z_soisno (maxsnl+1:0) !node depth [m]
  REAL(r8), intent(out) :: dz_soisno(maxsnl+1:0) !layer thickness [m]
  INTEGER,  intent(out) :: snl                   !number of snow layer
  REAL(r8) zi
  INTEGER i
!-----------------------------------------------------------------------

  dz_soisno(:0) = 0.
  z_soisno(:0) = 0.
  snl = 0
  IF(patchtype.le.3)THEN !non water bodies

     IF(snowdp.lt.0.01)THEN
        snl = 0
     ELSE
        IF(snowdp>=0.01 .and. snowdp<=0.03)THEN
           snl = -1
           dz_soisno(0)  = snowdp
        ELSE IF(snowdp>0.03 .and. snowdp<=0.04)THEN
           snl = -2
           dz_soisno(-1) = snowdp/2.
           dz_soisno( 0) = dz_soisno(-1)
        ELSE IF(snowdp>0.04 .and. snowdp<=0.07)THEN
           snl = -2
           dz_soisno(-1) = 0.02
           dz_soisno( 0) = snowdp - dz_soisno(-1)
        ELSE IF(snowdp>0.07 .and. snowdp<=0.12)THEN
           snl = -3
           dz_soisno(-2) = 0.02
           dz_soisno(-1) = (snowdp - 0.02)/2.
           dz_soisno( 0) = dz_soisno(-1)
        ELSE IF(snowdp>0.12 .and. snowdp<=0.18)THEN
           snl = -3
           dz_soisno(-2) = 0.02
           dz_soisno(-1) = 0.05
           dz_soisno( 0) = snowdp - dz_soisno(-2) - dz_soisno(-1)
        ELSE IF(snowdp>0.18 .and. snowdp<=0.29)THEN
           snl = -4
           dz_soisno(-3) = 0.02
           dz_soisno(-2) = 0.05
           dz_soisno(-1) = (snowdp - dz_soisno(-3) - dz_soisno(-2))/2.
           dz_soisno( 0) = dz_soisno(-1)
        ELSE IF(snowdp>0.29 .and. snowdp<=0.41)THEN
           snl = -4
           dz_soisno(-3) = 0.02
           dz_soisno(-2) = 0.05
           dz_soisno(-1) = 0.11
           dz_soisno( 0) = snowdp - dz_soisno(-3) - dz_soisno(-2) - dz_soisno(-1)
        ELSE IF(snowdp>0.41 .and. snowdp<=0.64)THEN
           snl = -5
           dz_soisno(-4) = 0.02
           dz_soisno(-3) = 0.05
           dz_soisno(-2) = 0.11
           dz_soisno(-1) = (snowdp - dz_soisno(-4) - dz_soisno(-3) - dz_soisno(-2))/2.
           dz_soisno( 0) = dz_soisno(-1)
        ELSE IF(snowdp>0.64)THEN
           snl = -5
           dz_soisno(-4) = 0.02
           dz_soisno(-3) = 0.05
           dz_soisno(-2) = 0.11
           dz_soisno(-1) = 0.23
           dz_soisno( 0) = snowdp - dz_soisno(-4) - dz_soisno(-3) - dz_soisno(-2) - dz_soisno(-1)
        ENDIF

        zi = 0.
        DO i = 0, snl+1, -1
           z_soisno(i) = zi - dz_soisno(i)/2.
           zi = -zi-dz_soisno(i)
        ENDDO
     ENDIF

  ENDIF

END SUBROUTINE snow_ini
!-----------------------------------------------------------------------
! EOP


SUBROUTINE polint(xa,ya,n,x,y)

! Given arrays xa and ya, each of length n, and gi
! value y, and an error estimate dy. If P (x) is the p
! P (xa(i)) = ya(i), i = 1, . . . , n, then the returned value
! (from: "Numerical Recipes")

  USE precision
  IMPLICIT NONE
  INTEGER n,NMAX
  REAL(r8) dy,x,y,xa(n),ya(n)
  parameter (NMAX=10)      !Largest anticipated val
  INTEGER i,m,ns
  REAL(r8) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

  ns=1
  dif=abs(x-xa(1))

  DO i=1,n       !Here we find the index ns of the closest table entry,
     dift=abs(x-xa(i))
     IF(dift.lt.dif) THEN
        ns=i
        dif=dift
     ENDIF
     c(i)=ya(i)  !and initialize the tableau of c's and d's.
     d(i)=ya(i)
  ENDDO

  y=ya(ns)       !This is the initial approximation to y.
  ns=ns-1

  DO m=1,n-1  !For each column of the tableau,
     DO i=1,n-m   !we loop over the current c's and d's and update them.
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        IF(den.eq.0.) print*, 'failure in polint'  !two input xa's are identical.
        den=w/den
        d(i)=hp*den                                !here the c's and d's are updated.
        c(i)=ho*den
     ENDDO
     IF(2*ns.lt.n-m)THEN  !After each column in the tableau is completed, we decide
        dy=c(ns+1)        !which correction, c or d, we want to add to our accumulating
     ELSE                 !value of y, i.e., which path to take through
        dy=d(ns)          !the tableau-forking up or down. We DO this in such a
        ns=ns-1           !way as to take the most "straight line" route through the
     ENDIF                !tableau to its apex, updating ns accordingly to keep track
     y=y+dy               !of where we are. This route keeps the partial approximations
  ENDDO                   !centered (insofar as possible) on the target x. T he
                          !last dy added is thus the error indication.
END SUBROUTINE polint
!-----------------------------------------------------------------------
! EOP
