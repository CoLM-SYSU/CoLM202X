#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Ensemble
!-----------------------------------------------------------------------
! DESCRIPTION:
!    Provide functions to generate ensemble samples for data assimilation
!
! REFERENCES:
!    [1] SMOS brightness temperature assimilation into the Community Land Model
!        Hydrol. Earth Syst. Sci., 21, 5929–5951, 2017,
!        https://doi.org/10.5194/hess-21-5929-2017
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Vars_TimeVariables
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Vars_1DForcing
   USE MOD_LandPatch
   IMPLICIT NONE
   SAVE

   ! public functions
   PUBLIC :: ensemble

   ! local parameters
   real(r8) :: std_sr = 0.3 ! multiplicative noise
   real(r8) :: std_ta = 1.0 ! additive noise
   real(r8) :: std_lw  = 20.0 ! multiplicative noise, log-normal noise
   real(r8) :: std_p  = 0.5 ! multiplicative noise, log-normal noise

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE ensemble()

!-----------------------------------------------------------------------
      IMPLICIT NONE

!------------------------ Local Variables ------------------------------
      integer  ::  np

!-----------------------------------------------------------------------

!#############################################################################
! Generate ensemble samples for forcing variables
!#############################################################################

      DO np = 1, numpatch
         ! generate ensemble samples for forcing variables
         CALL disturb_ens(DEF_DA_ENS, 2, forc_prc(np),   0.0, std_p,  forc_prc_ens(:,np))
         CALL disturb_ens(DEF_DA_ENS, 2, forc_prl(np),   0.0, std_p,  forc_prl_ens(:,np))
         CALL disturb_ens(DEF_DA_ENS, 0, forc_t(np),     0.0, std_ta, forc_t_ens(:,np))
         CALL disturb_ens(DEF_DA_ENS, 0, forc_frl(np),   0.0, std_lw, forc_frl_ens(:,np))
         CALL disturb_ens(DEF_DA_ENS, 1, forc_sols(np),  0.0, std_sr, forc_sols_ens(:,np))
         CALL disturb_ens(DEF_DA_ENS, 1, forc_soll(np),  0.0, std_sr, forc_soll_ens(:,np))
         CALL disturb_ens(DEF_DA_ENS, 1, forc_solsd(np), 0.0, std_sr, forc_solsd_ens(:,np))
         CALL disturb_ens(DEF_DA_ENS, 1, forc_solld(np), 0.0, std_sr, forc_solld_ens(:,np))
      ENDDO

   END SUBROUTINE ensemble

!-----------------------------------------------------------------------

   SUBROUTINE disturb_ens(num_ens, type, A, mu, sigma, A_ens)

!-----------------------------------------------------------------------
      USE MOD_Vars_Global, only: pi
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)  :: num_ens
      integer, intent(in)  :: type
      real(r8), intent(in)  :: A
      real(r8), intent(in)  :: mu
      real(r8), intent(in)  :: sigma
      real(r8), intent(out) :: A_ens(num_ens)

!------------------------ Local Variables ------------------------------
      real(r8) :: u1, u2
      real(r8) :: z(num_ens)
      real(r8) :: mean_z, max_z
      integer  :: i

!-----------------------------------------------------------------------

      ! generate disturbance ensemble samples ~ N(mu, sigma)
      DO i = 1, num_ens
         CALL random_number(u1)
         CALL random_number(u2)
         u1 = max(u1, 1e-10)
         u2 = max(u2, 1e-10)
         z(i) = sqrt(-2.0*log(u1))*cos(2.0*pi*u2)
      ENDDO
      z = sigma*z + mu

      ! normalize the disturbance ensemble samples to mean 0
      mean_z = sum(z)/num_ens
      DO i = 1, num_ens
         z(i) = z(i) - mean_z
      ENDDO

      ! generate ensemble samples according different types (0: additive, 1: multiplicative)
      IF (type == 0) THEN ! additive normal
         DO i = 1, num_ens
            A_ens(i) = A + z(i)
         ENDDO
      ELSEIF (type == 1) THEN ! multiplicative normal
         DO i = 1, num_ens
            A_ens(i) = A*(1.0 + z(i))
         ENDDO
      ELSEIF (type == 2) THEN ! multiplicative log-normal
         max_z = maxval(abs(z))
         DO i = 1, num_ens
            A_ens(i) = A*(1 + z(i)/max_z)
         ENDDO
      ENDIF

   END SUBROUTINE disturb_ens

!-----------------------------------------------------------------------
END MODULE MOD_DA_Ensemble
#endif
