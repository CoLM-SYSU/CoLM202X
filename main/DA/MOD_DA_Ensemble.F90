#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Ensemble
!-----------------------------------------------------------------------
! DESCRIPTION:
!    Provide functions to generate ensemble samples for data assimilation
!
! REFERENCES:
!    [1] Algorithm Theoretical Basis Document Level 4 Surface and Root
!        Zone Soil Moisture (L4_SM) Data Product
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!   Lu Li, 10/2025: Consider correlation & AR(1) process
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
   ! forcing [parameters used here is consistent with SMAP L4 (Table 4 in [1])]
   integer,  parameter :: nvar = 4                                   ! number of pertutated forcing variables
   real(r8), parameter :: tau_ar = 24.0                              ! correlation time scale in hours

   real(r8) :: dt                                                    ! time step in hours
   real(r8) :: phi                                                   ! AR(1) autocorrelation coefficient consider time scale
   real(r8) :: sigma_eps                                             ! standard deviation of noise in AR(1) process
   real(r8), dimension(nvar) :: sigma = (/0.5, 0.3, 20.0, 1.0/)      ! standard deviation of perturbed forcing variables (prcp, sw, lw, t)
   real(r8), dimension(nvar, nvar) :: C = reshape([ &
         1.0, -0.8,  0.5, 0.0, &
        -0.8,  1.0, -0.5, 0.4, &
         0.5, -0.5,  1.0, 0.4, &
         0.0,  0.4,  0.4, 1.0], shape=[nvar, nvar])                  ! cross-correlation matrix between perturbed forcing variables
   real(r8), allocatable :: r_prev(:,:,:)                            ! previous perturbation (numpatch, nvar, DEF_DA_ENS_NUM)
   real(r8), allocatable :: r_curr(:,:,:)                            ! current perturbation (numpatch, nvar, DEF_DA_ENS_NUM)
   logical :: initialized = .false.                                  ! flag to indicate if is initialized

   ! soil moisture [default set 0.002 m3/m3 disterbulance]
   integer,  parameter :: nvar_sm = 2                                ! number of pertutated soil moisture layers
   real(r8), parameter :: tau_sm = 3.0                               ! correlation time scale in hours
   real(r8) :: phi_sm                                                ! AR(1) autocorrelation coefficient consider time scale
   real(r8) :: sigma_eps_sm                                          ! standard deviation of noise in AR(1) process
   real(r8), dimension(nvar_sm) :: sigma_sm = (/0.035, 0.0552/)      ! standard deviation of perturbed soil moisture (equal to 0.002 m3/m3)
   real(r8), allocatable :: r_prev_sm(:,:,:)                         ! previous perturbation (numpatch, nvar_sm, DEF_DA_ENS_NUM)
   real(r8), allocatable :: r_curr_sm(:,:,:)                         ! current perturbation (numpatch, nvar_sm, DEF_DA_ENS_NUM)

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE ensemble(deltim)

!-----------------------------------------------------------------------
      IMPLICIT NONE

      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------
      integer  ::  np, i, j

      real(r8) :: cov_matrix(nvar, nvar)                      ! covariance matrix between perturbed forcing variables
      real(r8) :: L(nvar, nvar)                               ! Cholesky decomposition of correlation matrix
      integer  :: info                                        ! info flag for Cholesky decomposition
      real(r8) :: u1(DEF_DA_ENS_NUM/2), u2(DEF_DA_ENS_NUM/2)  ! uniform random variables
      real(r8) :: z(nvar, DEF_DA_ENS_NUM)                     ! standard normal random variables
      real(r8) :: mean_z(nvar)                                ! mean of perturbation (nvar)
      real(r8) :: std_z(nvar)                                 ! std of perturbation (nvar)
      real(r8) :: zxL(nvar, DEF_DA_ENS_NUM)                   ! correlated random variables (nvar, DEF_DA_ENS_NUM)

      real(r8) :: z_sm(DEF_DA_ENS_NUM)                        ! standard normal random variables
      real(r8) :: mean_z_sm                                   ! mean of perturbation
      real(r8) :: std_z_sm                                    ! std of perturbation
      real(r8) :: mean_r_sm(nvar_sm)                          ! mean of perturbation for soil moisture (nvar_sm)
      real(r8) :: std_r_sm(nvar_sm)                           ! std of perturbation for soil moisture (nvar_sm)
      real(r8) :: a1(DEF_DA_ENS_NUM)                          ! temporary disturbed variable for soil moisture layer 1
      real(r8) :: a2(DEF_DA_ENS_NUM)                          ! temporary disturbed variable for soil moisture layer 2

!-----------------------------------------------------------------------

      ! initialize persistent variables
      IF (.not. initialized) THEN
         allocate(r_prev(numpatch, nvar, DEF_DA_ENS_NUM))
         allocate(r_curr(numpatch, nvar, DEF_DA_ENS_NUM))
         allocate(r_prev_sm(numpatch, nvar_sm, DEF_DA_ENS_NUM))
         allocate(r_curr_sm(numpatch, nvar_sm, DEF_DA_ENS_NUM))
         r_prev = 0.0_r8
         r_curr = 0.0_r8
         r_prev_sm = 0.0_r8
         r_curr_sm = 0.0_r8
         initialized = .true.
      ENDIF

      ! calculate AR(1) parameters
      dt = deltim/3600
      phi = exp(-dt/tau_ar)
      sigma_eps = sqrt(1.0 - phi**2)
      phi_sm = exp(-dt/tau_sm)
      sigma_eps_sm = sqrt(1.0 - phi_sm**2)

      ! calculate covariance matrix by cross-correlation matrix and standard deviation
      cov_matrix = 0.0_r8
      DO i = 1, nvar
         DO j = 1, nvar
            cov_matrix(i,j) = C(i,j) * sigma(i) * sigma(j)
         ENDDO
      ENDDO

      ! perform Cholesky decomposition of covariance matrix
      L = cov_matrix
      CALL dpotrf('L', nvar, L, nvar, info)
      DO i = 1, nvar
         DO j = i+1, nvar
            L(i,j) = 0.0_r8
         ENDDO
      ENDDO
      IF (info /= 0) THEN
         print *, 'Error: Cholesky decomposition failed'
         stop
      ENDIF

      ! Generate ensemble samples for forcing variables
      DO np = 1, numpatch
         ! generate disturbance ensemble samples ~ N(0, I)
         CALL random_seed()
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)  ! ensure u1 is not zero
            z(1,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z(1,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO
         CALL random_seed()
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)
            z(2,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z(2,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO
         CALL random_seed()
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)
            z(3,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z(3,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO
         CALL random_seed()
         CALL random_number(u1)
         CALL random_number(u2)
         DO i = 1, DEF_DA_ENS_NUM/2
            u1(i) = max(u1(i), 1e-10)
            z(4,i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
            z(4,i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
         ENDDO

         ! normalize z to mean 0 and std 1
         DO i = 1, nvar
            mean_z(i) = sum(z(i, :))/DEF_DA_ENS_NUM
            std_z(i) = sqrt(sum((z(i, :) - mean_z(i))**2)/(DEF_DA_ENS_NUM - 1))
         ENDDO
         DO i = 1, nvar
            z(i,:) = (z(i,:)-mean_z(i))/std_z(i)
         ENDDO

         ! multiply by Cholesky factor to introduce correlation (z*L)
         CALL dgemm('N', 'N', nvar, DEF_DA_ENS_NUM, nvar, 1.0_r8, L, nvar, z, nvar, 0.0_r8, zxL, nvar)

         ! introduce correlation using AR(1) process
         DO i = 1, nvar
            DO j = 1, DEF_DA_ENS_NUM
               r_curr(np,i,j) = phi * r_prev(np,i,j) + sigma_eps * zxL(i,j)
            ENDDO
         ENDDO
         ! no AR(1) process, directly use correlated random variables
         ! r_curr(np,:,:) = zxL

         ! normalize the disturbance ensemble samples to mean 0
         mean_z = sum(r_curr(np,:,:), dim=2)/DEF_DA_ENS_NUM
         DO i = 1, nvar
            DO j = 1, DEF_DA_ENS_NUM
               r_curr(np,i,j) = r_curr(np,i,j) - mean_z(i)
            ENDDO
         ENDDO

         ! save current perturbation as previous perturbation for next time step
         r_prev = r_curr

         ! generate ensemble samples according different types
         DO j = 1, DEF_DA_ENS_NUM
            forc_prc_ens(j,np) = forc_prc(np) * exp(r_curr(np,1,j) - 0.5 * sigma(1)**2)
            forc_prl_ens(j,np) = forc_prl(np) * exp(r_curr(np,1,j) - 0.5 * sigma(1)**2)
            forc_sols_ens(j,np) = forc_sols(np) * exp(r_curr(np,2,j) - 0.5 * sigma(2)**2)
            forc_soll_ens(j,np) = forc_soll(np) * exp(r_curr(np,2,j) - 0.5 * sigma(2)**2)
            forc_solsd_ens(j,np) = forc_solsd(np) * exp(r_curr(np,2,j) - 0.5 * sigma(2)**2)
            forc_solld_ens(j,np) = forc_solld(np) * exp(r_curr(np,2,j) - 0.5 * sigma(2)**2)
            forc_frl_ens(j,np) = forc_frl(np) + r_curr(np,3,j)
            forc_t_ens(j,np) = forc_t(np) + r_curr(np,4,j)
         ENDDO

         IF (DEF_DA_ENS_SM) THEN
            ! generate ensemble samples (0, I) for soil moisture
            CALL random_seed()
            CALL random_number(u1)
            CALL random_number(u2)
            DO i = 1, DEF_DA_ENS_NUM/2
               u1(i) = max(u1(i), 1e-10)
               z_sm(i*2-1) = sqrt(-2.0*log(u1(i))) * cos(2.0*pi*u2(i))
               z_sm(i*2)   = sqrt(-2.0*log(u1(i))) * sin(2.0*pi*u2(i))
            ENDDO
            mean_z_sm = sum(z_sm)/DEF_DA_ENS_NUM
            std_z_sm = sqrt(sum((z_sm - mean_z_sm)**2)/(DEF_DA_ENS_NUM - 1))
            z_sm = (z_sm - mean_z_sm)/std_z_sm

            ! introduce correlation using AR(1) process
            DO i = 1, nvar_sm
               DO j = 1, DEF_DA_ENS_NUM
                  r_curr_sm(np,i,j) = phi_sm * r_prev_sm(np,i,j) + sigma_eps_sm * sigma_sm(i) * z_sm(j)
               ENDDO
            ENDDO

            ! normalize the disturbance ensemble samples to mean 0
            mean_r_sm = sum(r_curr_sm(np,:,:), dim=2)/DEF_DA_ENS_NUM
            DO i = 1, nvar_sm
               DO j = 1, DEF_DA_ENS_NUM
                  r_curr_sm(np,i,j) = r_curr_sm(np,i,j) - mean_r_sm(i)
               ENDDO
            ENDDO

            ! save current perturbation as previous perturbation for next time step
            r_prev_sm = r_curr_sm

            ! generate ensemble samples according different types
            DO j = 1, DEF_DA_ENS_NUM
               a1(j) = wliq_soisno_ens(1,j,np) + r_curr_sm(np,1,j)
               a2(j) = wliq_soisno_ens(2,j,np) + r_curr_sm(np,2,j)
            ENDDO
            DO j = 1, DEF_DA_ENS_NUM
               a1(j) = max(1e-10, a1(j))
               a2(j) = max(1e-10, a2(j))
            ENDDO

            ! move residual water to water table
            DO j = 1, DEF_DA_ENS_NUM
               wa_ens(j, np) = wa_ens(j, np) - (sum(wliq_soisno_ens(1:2, j, np)) - a1(j) - a2(j))
               wliq_soisno_ens(1, j, np) = a1(j)
               wliq_soisno_ens(2, j, np) = a2(j)
            ENDDO
         ENDIF
      ENDDO

   END SUBROUTINE ensemble

!-----------------------------------------------------------------------
END MODULE MOD_DA_Ensemble
#endif
