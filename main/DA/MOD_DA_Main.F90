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
   USE MOD_LandPatch
   USE MOD_DA_GRACE
   USE MOD_DA_SMAP
   USE MOD_DA_FY3D
   USE MOD_Vars_1DFluxes
   USE MOD_Vars_TimeVariables
   USE MOD_DA_Vars_1DFluxes
   USE MOD_DA_Vars_TimeVariables
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
! Init data assimilation for different satellite products
!#############################################################################
   IF (DEF_DA_GRACE) THEN
      CALL init_DA_GRACE ()
   ENDIF

   IF (DEF_DA_SMAP) THEN
      CALL init_DA_SMAP  ()
   ENDIF

   IF (DEF_DA_FY3D) THEN
      CALL init_DA_FY3D  ()
   ENDIF

   END SUBROUTINE init_DA

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA (idate, deltim)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
   integer :: np

!-----------------------------------------------------------------------------

!#############################################################################
! Perform data assimilation for different satellite products
!#############################################################################
      IF (DEF_DA_GRACE) THEN
         CALL run_DA_GRACE (idate, deltim)
      ENDIF

      IF (DEF_DA_SMAP) THEN
         CALL run_DA_SMAP  (idate, deltim)
      ENDIF

      IF (DEF_DA_FY3D) THEN
         CALL run_DA_FY3D  (idate, deltim)
      ENDIF

!#############################################################################
! Use ensemble mean for important outputs
!#############################################################################
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            DO np = 1, numpatch
               ! important state variables
               t_soisno(:,np) = sum(t_soisno_ens(:,:,np), dim=2) / DEF_DA_ENS
               wliq_soisno(:,np) = sum(wliq_soisno_ens(:,:,np), dim=2) / DEF_DA_ENS
               wice_soisno(:,np) = sum(wice_soisno_ens(:,:,np), dim=2) / DEF_DA_ENS
               t_grnd(np) = sum(t_grnd_ens(:,np), dim=1) / DEF_DA_ENS
               tleaf(np) = sum(tleaf_ens(:,np), dim=1) / DEF_DA_ENS
               snowdp(np) = sum(snowdp_ens(:,np), dim=1) / DEF_DA_ENS

               ! diagnostic variables
               h2osoi(:,np) = sum(h2osoi_ens(:,:,np), dim=2) / DEF_DA_ENS
               t_brt(:,np) = sum(t_brt_ens(:,:,np), dim=2) / DEF_DA_ENS
               trad(np) = sum(trad_ens(:,np), dim=1) / DEF_DA_ENS
               tref(np) = sum(tref_ens(:,np), dim=1) / DEF_DA_ENS
               qref(np) = sum(qref_ens(:,np), dim=1) / DEF_DA_ENS
               rsur(np) = sum(rsur_ens(:,np), dim=1) / DEF_DA_ENS
               ustar(np) = sum(ustar_ens(:,np), dim=1) / DEF_DA_ENS
               qstar(np) = sum(qstar_ens(:,np), dim=1) / DEF_DA_ENS
               tstar(np) = sum(tstar_ens(:,np), dim=1) / DEF_DA_ENS
               fm(np) = sum(fm_ens(:,np), dim=1) / DEF_DA_ENS
               fh(np) = sum(fh_ens(:,np), dim=1) / DEF_DA_ENS
               fq(np) = sum(fq_ens(:,np), dim=1) / DEF_DA_ENS
               fsena(np) = sum(fsena_ens(:,np), dim=1) / DEF_DA_ENS
               fevpa(np) = sum(fevpa_ens(:,np), dim=1) / DEF_DA_ENS
               lfevpa(np) = sum(lfevpa_ens(:,np), dim=1) / DEF_DA_ENS   
               rsur(np) = sum(rsur_ens(:,np), dim=1) / DEF_DA_ENS
            ENDDO
         ENDIF
      ENDIF

   END SUBROUTINE run_DA

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!-----------------------------------------------------------------------------

   IF (DEF_DA_GRACE) THEN
      CALL end_DA_GRACE ()
   ENDIF

   IF (DEF_DA_SMAP) THEN
      CALL end_DA_SMAP  ()
   ENDIF

   IF (DEF_DA_FY3D) THEN
      CALL end_DA_FY3D  ()
   ENDIF

   END SUBROUTINE end_DA

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Main
#endif
