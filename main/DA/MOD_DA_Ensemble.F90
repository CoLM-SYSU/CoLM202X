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
   IMPLICIT NONE
   SAVE

   ! public functions
   PUBLIC :: disturb_forc_ens

   ! local parameters
   real(r8) :: std_sr = 0.3 ! multiplicative noise
   real(r8) :: std_ta = 2.5 ! additive noise
   real(r8) :: std_p  = 0.3 ! multiplicative noise, log-normal noise

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE disturb_ens(num_ens, type, rand_seed, A, mu, sigma, A_ens)

!-----------------------------------------------------------------------
      USE MOD_Vars_Global, only: pi
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)  :: num_ens
      integer, intent(in)  :: rand_seed
      integer, intent(in)  :: type
      real(r8), intent(in)  :: A
      real(r8), intent(in)  :: mu
      real(r8), intent(in)  :: sigma
      real(r8), intent(out) :: A_ens(num_ens)

!------------------------ Local Variables ------------------------------
      real(r8) :: u1, u2
      real(r8) :: z(num_ens)
      real(r8) :: mean_z
      integer  :: i
      integer  :: seed_size
      integer, allocatable :: seed(:)

!-----------------------------------------------------------------------

      ! initialize random number generator
      CALL random_seed(size=seed_size)
      allocate (seed(seed_size))
      seed = rand_seed
      CALL random_seed(put=seed)

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
         DO i = 1, num_ens
            A_ens(i) = A*exp(z(i))
         ENDDO
      ENDIF

      deallocate (seed)

   END SUBROUTINE disturb_ens

!-----------------------------------------------------------------------

   SUBROUTINE disturb_forc_ens( &
      idate, &
      num_ens, &
      forc_t, &
      forc_prc, forc_prl, &
      forc_sols, forc_soll, forc_solsd, forc_solld, &
      forc_t_ens, &
      forc_prc_ens, forc_prl_ens, &
      forc_sols_ens, forc_soll_ens, forc_solsd_ens, forc_solld_ens)

!-----------------------------------------------------------------------
      USE MOD_TimeManager, only: minutes_since_1900
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in) :: num_ens
      integer, intent(in) :: idate(3)  
      real(r8), intent(in) :: &
         forc_t, &
         forc_prc, &
         forc_prl, &
         forc_sols, &
         forc_soll, &
         forc_solsd, &
         forc_solld

      real(r8), intent(out) :: &
         forc_t_ens(num_ens), &
         forc_prc_ens(num_ens), &
         forc_prl_ens(num_ens), &
         forc_sols_ens(num_ens), &
         forc_soll_ens(num_ens), &
         forc_solsd_ens(num_ens), &
         forc_solld_ens(num_ens)

!------------------------ Local Variables ------------------------------
      integer  ::  rand_seed

!-----------------------------------------------------------------------
      ! to guarantee reproducibility, all random processes are initialized with a fixed seed.
      rand_seed = minutes_since_1900(idate(1), idate(2), idate(3))

      ! generate ensemble samples for forcing variables
      CALL disturb_ens(num_ens, 0, rand_seed + 0, forc_t,     0.0, std_ta, forc_t_ens)
      CALL disturb_ens(num_ens, 2, rand_seed + 1, forc_prc,   0.0, std_p,  forc_prc_ens)
      CALL disturb_ens(num_ens, 2, rand_seed + 2, forc_prl,   0.0, std_p,  forc_prl_ens)
      CALL disturb_ens(num_ens, 1, rand_seed + 3, forc_sols,  0.0, std_sr, forc_sols_ens)
      CALL disturb_ens(num_ens, 1, rand_seed + 4, forc_soll,  0.0, std_sr, forc_soll_ens)
      CALL disturb_ens(num_ens, 1, rand_seed + 5, forc_solsd, 0.0, std_sr, forc_solsd_ens)
      CALL disturb_ens(num_ens, 1, rand_seed + 6, forc_solld, 0.0, std_sr, forc_solld_ens)

   END SUBROUTINE disturb_forc_ens

!-----------------------------------------------------------------------
END MODULE MOD_DA_Ensemble
#endif