#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Main
!-----------------------------------------------------------------------------
! DESCRIPTION:
!     Main procedures for data assimilation
!
! AUTHOR:
!     Lu Li, 12/2024
!     Zhilong Fan, Lu Li, 03/2024
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Spmd_task
   USE MOD_Namelist
   USE MOD_TimeManager
   USE MOD_LandPatch
   USE MOD_DA_TWS
   USE MOD_DA_SM
   USE MOD_DA_Ensemble
   USE MOD_Vars_1DFluxes
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_1DForcing
   IMPLICIT NONE
   SAVE

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!-----------------------------------------------------------------------------

!#############################################################################
! Init data assimilation for different products
!#############################################################################
   IF (DEF_DA_TWS) THEN
      IF (DEF_DA_TWS_GRACE) THEN
         CALL init_DA_GRACE ()
      ENDIF
   ENDIF

   IF (DEF_DA_SM) THEN
      IF (p_is_master) THEN
         print *, '[CoLM-DA] initialize surface soil moisture & temperature data assimilation.'
      ENDIF
      CALL init_DA_SM ()
   ENDIF

   END SUBROUTINE init_DA

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA (idate, deltim, dolai, doalb, dosst, oro)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim
   logical,  intent(in) :: dolai    ! true if time for time-varying vegetation parameter
   logical,  intent(in) :: doalb    ! true if time for surface albedo calculation
   logical,  intent(in) :: dosst    ! true if time for update sst/ice/snow

   real(r8), intent(inout) :: oro(numpatch)  ! ocean(0)/seaice(2)/ flag

!------------------------ Local Variables ------------------------------------
   integer :: np, i, maxday
   integer :: sdate(3) ! backward date from the end to the begin of the time step

!-----------------------------------------------------------------------------

      ! use OI for terrestrial water storage
      IF (DEF_DA_TWS) THEN
         IF (DEF_DA_TWS_GRACE) THEN
            CALL run_DA_GRACE (idate, deltim)
         ENDIF
      ENDIF

      ! use ensemble DA for soil moisture
      IF ((DEF_DA_SM) .and. (DEF_DA_ENS_NUM > 1))THEN
!#############################################################################
! Generate ensemble members
!#############################################################################
         IF (p_is_worker) THEN

            !  store the non-DA trajectory
            z_sno_noda        = z_sno
            dz_sno_noda       = dz_sno
            t_soisno_noda     = t_soisno
            wliq_soisno_noda  = wliq_soisno
            wice_soisno_noda  = wice_soisno
            smp_noda          = smp
            hk_noda           = hk
            t_grnd_noda       = t_grnd
            tleaf_noda        = tleaf
            ldew_noda         = ldew
            ldew_rain_noda    = ldew_rain
            ldew_snow_noda    = ldew_snow
            fwet_snow_noda    = fwet_snow
            sag_noda          = sag
            scv_noda          = scv
            snowdp_noda       = snowdp
            fveg_noda         = fveg
            fsno_noda         = fsno
            sigf_noda         = sigf
            green_noda        = green
            tlai_noda         = tlai
            lai_noda          = lai
            tsai_noda         = tsai
            sai_noda          = sai
            alb_noda          = alb
            ssun_noda         = ssun
            ssha_noda         = ssha
            ssoi_noda         = ssoi
            ssno_noda         = ssno
            thermk_noda       = thermk
            extkb_noda        = extkb
            extkd_noda        = extkd
            zwt_noda          = zwt
            wdsrf_noda        = wdsrf
            wa_noda           = wa
            wetwat_noda       = wetwat
            t_lake_noda       = t_lake
            lake_icefrac_noda = lake_icefrac
            savedtke1_noda    = savedtke1

            tref_noda         = tref
            h2osoi_noda       = h2osoi

            ! Generate ensemble members
            CALL ensemble (deltim)

            ! loop over ensemble members
            DO i = 1, DEF_DA_ENS_NUM
               ! set ensemble forcing variables
               forc_t     = forc_t_ens(i,:)
               forc_prc   = forc_prc_ens(i,:)
               forc_prl   = forc_prl_ens(i,:)
               forc_sols  = forc_sols_ens(i,:)
               forc_soll  = forc_soll_ens(i,:)
               forc_solsd = forc_solsd_ens(i,:)
               forc_solld = forc_solld_ens(i,:)
               forc_frl   = forc_frl_ens(i,:)

               ! give the i-th trajectory state value to state variables
               z_sno        = z_sno_ens       (  :,i,:)
               dz_sno       = dz_sno_ens      (  :,i,:)
               t_soisno     = t_soisno_ens    (  :,i,:)
               wliq_soisno  = wliq_soisno_ens (  :,i,:)
               wice_soisno  = wice_soisno_ens (  :,i,:)
               smp          = smp_ens         (  :,i,:)
               hk           = hk_ens          (  :,i,:)
               t_grnd       = t_grnd_ens      (    i,:)
               tleaf        = tleaf_ens       (    i,:)
               ldew         = ldew_ens        (    i,:)
               ldew_rain    = ldew_rain_ens   (    i,:)
               ldew_snow    = ldew_snow_ens   (    i,:)
               fwet_snow    = fwet_snow_ens   (    i,:)
               sag          = sag_ens         (    i,:)
               scv          = scv_ens         (    i,:)
               snowdp       = snowdp_ens      (    i,:)
               fveg         = fveg_ens        (    i,:)
               fsno         = fsno_ens        (    i,:)
               sigf         = sigf_ens        (    i,:)
               green        = green_ens       (    i,:)
               tlai         = tlai_ens        (    i,:)
               lai          = lai_ens         (    i,:)
               tsai         = tsai_ens        (    i,:)
               sai          = sai_ens         (    i,:)
               alb          = alb_ens         (:,:,i,:)
               ssun         = ssun_ens        (:,:,i,:)
               ssha         = ssha_ens        (:,:,i,:)
               ssoi         = ssoi_ens        (:,:,i,:)
               ssno         = ssno_ens        (:,:,i,:)
               thermk       = thermk_ens      (    i,:)
               extkb        = extkb_ens       (    i,:)
               extkd        = extkd_ens       (    i,:)
               zwt          = zwt_ens         (    i,:)
               wdsrf        = wdsrf_ens       (    i,:)
               wa           = wa_ens          (    i,:)
               wetwat       = wetwat_ens      (    i,:)
               t_lake       = t_lake_ens      (  :,i,:)
               lake_icefrac = lake_icefrac_ens(  :,i,:)
               savedtke1    = savedtke1_ens   (    i,:)

               ! run colm
               CALL CoLMDRIVER (idate, deltim, dolai, doalb, dosst, oro)

               ! output ensemble members
               z_sno_ens       (  :,i,:) = z_sno
               dz_sno_ens      (  :,i,:) = dz_sno
               t_soisno_ens    (  :,i,:) = t_soisno
               wliq_soisno_ens (  :,i,:) = wliq_soisno
               wice_soisno_ens (  :,i,:) = wice_soisno
               smp_ens         (  :,i,:) = smp
               hk_ens          (  :,i,:) = hk
               t_grnd_ens      (    i,:) = t_grnd
               tleaf_ens       (    i,:) = tleaf
               ldew_ens        (    i,:) = ldew
               ldew_rain_ens   (    i,:) = ldew_rain
               ldew_snow_ens   (    i,:) = ldew_snow
               fwet_snow_ens   (    i,:) = fwet_snow
               sag_ens         (    i,:) = sag
               scv_ens         (    i,:) = scv
               snowdp_ens      (    i,:) = snowdp
               fveg_ens        (    i,:) = fveg
               fsno_ens        (    i,:) = fsno
               sigf_ens        (    i,:) = sigf
               green_ens       (    i,:) = green
               tlai_ens        (    i,:) = tlai
               lai_ens         (    i,:) = lai
               tsai_ens        (    i,:) = tsai
               sai_ens         (    i,:) = sai
               alb_ens         (:,:,i,:) = alb
               ssun_ens        (:,:,i,:) = ssun
               ssha_ens        (:,:,i,:) = ssha
               ssoi_ens        (:,:,i,:) = ssoi
               ssno_ens        (:,:,i,:) = ssno
               thermk_ens      (    i,:) = thermk
               extkb_ens       (    i,:) = extkb
               extkd_ens       (    i,:) = extkd
               zwt_ens         (    i,:) = zwt
               wdsrf_ens       (    i,:) = wdsrf
               wa_ens          (    i,:) = wa
               wetwat_ens      (    i,:) = wetwat
               t_lake_ens      (  :,i,:) = t_lake
               lake_icefrac_ens(  :,i,:) = lake_icefrac
               savedtke1_ens   (    i,:) = savedtke1

               h2osoi_ens(:,i,:) = h2osoi
               trad_ens  (i,:)   = trad
               tref_ens  (i,:)   = tref
               qref_ens  (i,:)   = qref
               ustar_ens (i,:)   = ustar
               qstar_ens (i,:)   = qstar
               tstar_ens (i,:)   = tstar
               fm_ens    (i,:)   = fm
               fh_ens    (i,:)   = fh
               fq_ens    (i,:)   = fq

               fsena_ens (i,:)   = fsena
               lfevpa_ens(i,:)   = lfevpa
               fevpa_ens (i,:)   = fevpa
               rsur_ens  (i,:)   = rsur
            ENDDO

            !  recover the no-DA trajectory
            z_sno        = z_sno_noda
            dz_sno       = dz_sno_noda
            t_soisno     = t_soisno_noda
            wliq_soisno  = wliq_soisno_noda
            wice_soisno  = wice_soisno_noda
            smp          = smp_noda
            hk           = hk_noda
            t_grnd       = t_grnd_noda
            tleaf        = tleaf_noda
            ldew         = ldew_noda
            ldew_rain    = ldew_rain_noda
            ldew_snow    = ldew_snow_noda
            fwet_snow    = fwet_snow_noda
            sag          = sag_noda
            scv          = scv_noda
            snowdp       = snowdp_noda
            fveg         = fveg_noda
            fsno         = fsno_noda
            sigf         = sigf_noda
            green        = green_noda
            tlai         = tlai_noda
            lai          = lai_noda
            tsai         = tsai_noda
            sai          = sai_noda
            alb          = alb_noda
            ssun         = ssun_noda
            ssha         = ssha_noda
            ssoi         = ssoi_noda
            ssno         = ssno_noda
            thermk       = thermk_noda
            extkb        = extkb_noda
            extkd        = extkd_noda
            zwt          = zwt_noda
            wdsrf        = wdsrf_noda
            wa           = wa_noda
            wetwat       = wetwat_noda
            t_lake       = t_lake_noda
            lake_icefrac = lake_icefrac_noda
            savedtke1    = savedtke1_noda

            tref         = tref_noda
            h2osoi       = h2osoi_noda

         ENDIF
      ENDIF

!#############################################################################
! Perform data assimilation for different satellite products
!#############################################################################
      IF (DEF_DA_SM) THEN
         ! backward date from the end to the begin of the time step
         sdate = idate
         sdate(3) = sdate(3) - nint(deltim)
         IF (sdate(3) < 0) THEN
            sdate(2) = sdate(2) - 1
            sdate(3) = sdate(3) + 86400

            IF (sdate(2) < 1) THEN
               sdate(1) = sdate(1) - 1
               IF (isleapyear(sdate(1))) THEN
                  maxday = 366
               ELSE
                  maxday = 365
               ENDIF
               sdate(2) = maxday
            ENDIF
         ENDIF

         ! data assimilation
         IF (DEF_DA_SM) THEN
            IF (p_is_master) THEN
                print *, '[CoLM-DA] Start surface soil moisture & temperature data assimilation.'
            ENDIF
            CALL mpi_barrier (p_comm_glb, p_err)
            CALL run_DA_SM  (sdate, deltim)
         ENDIF
      ENDIF

!#############################################################################
! Use ensemble mean for important outputs
!#############################################################################
      IF ((DEF_DA_SM) .and. (DEF_DA_ENS_NUM > 1))THEN
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               DO np = 1, numpatch
                  ! important state variables
                  t_soisno_a(:,np) = sum(t_soisno_ens(:,:,np), dim=2) / DEF_DA_ENS_NUM
                  wliq_soisno_a(:,np) = sum(wliq_soisno_ens(:,:,np), dim=2) / DEF_DA_ENS_NUM
                  wice_soisno_a(:,np) = sum(wice_soisno_ens(:,:,np), dim=2) / DEF_DA_ENS_NUM
                  t_grnd_a(np) = sum(t_grnd_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  tleaf_a(np) = sum(tleaf_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  snowdp_a(np) = sum(snowdp_ens(:,np), dim=1) / DEF_DA_ENS_NUM

                  ! diagnostic variables
                  h2osoi_a(:,np) = sum(h2osoi_ens(:,:,np), dim=2) / DEF_DA_ENS_NUM
                  t_brt_smap_a(:,np) = sum(t_brt_smap_ens(:,:,np), dim=2) / DEF_DA_ENS_NUM
                  t_brt_fy3d_a(:,np) = sum(t_brt_fy3d_ens(:,:,np), dim=2) / DEF_DA_ENS_NUM
                  trad_a(np) = sum(trad_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  tref_a(np) = sum(tref_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  qref_a(np) = sum(qref_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  ustar_a(np) = sum(ustar_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  qstar_a(np) = sum(qstar_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  tstar_a(np) = sum(tstar_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  fm_a(np) = sum(fm_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  fh_a(np) = sum(fh_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  fq_a(np) = sum(fq_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  fsena_a(np) = sum(fsena_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  fevpa_a(np) = sum(fevpa_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  lfevpa_a(np) = sum(lfevpa_ens(:,np), dim=1) / DEF_DA_ENS_NUM
                  rsur_a(np) = sum(rsur_ens(:,np), dim=1) / DEF_DA_ENS_NUM
               ENDDO
            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE run_DA

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!-----------------------------------------------------------------------------

   IF (DEF_DA_TWS) THEN
      IF (DEF_DA_TWS_GRACE) THEN
         CALL end_DA_GRACE ()
      ENDIF
   ENDIF

   IF (DEF_DA_SM) THEN
      CALL end_DA_SM ()
   ENDIF

   END SUBROUTINE end_DA

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Main
#endif
