#include <define.h>

#ifdef LULC_IGBP_PC

MODULE MOD_Vars_1DPCFluxes
! -----------------------------------------------------------------
! !DESCRIPTION:
! Define Plant Community flux variables
!
! Created by Hua Yuan, 08/2019
! -----------------------------------------------------------------

  USE MOD_Precision
  IMPLICIT NONE
  SAVE

! -----------------------------------------------------------------
! Fluxes
! -----------------------------------------------------------------
  REAL(r8), allocatable :: fsenl_c   (:,:) !sensible heat from leaves [W/m2]
  REAL(r8), allocatable :: fevpl_c   (:,:) !evaporation+transpiration from leaves [mm/s]
  REAL(r8), allocatable :: etr_c     (:,:) !transpiration rate [mm/s]
  REAL(r8), allocatable :: fseng_c   (:,:) !sensible heat flux from ground [W/m2]
  REAL(r8), allocatable :: fevpg_c   (:,:) !evaporation heat flux from ground [mm/s]
  REAL(r8), allocatable :: parsun_c  (:,:) !solar absorbed by sunlit vegetation [W/m2]
  REAL(r8), allocatable :: parsha_c  (:,:) !solar absorbed by shaded vegetation [W/m2]
  REAL(r8), allocatable :: sabvsun_c (:,:) !solar absorbed by sunlit vegetation [W/m2]
  REAL(r8), allocatable :: sabvsha_c (:,:) !solar absorbed by shaded vegetation [W/m2]
  REAL(r8), allocatable :: qintr_c   (:,:) !interception (mm h2o/s)
  REAL(r8), allocatable :: qintr_rain_c(:,:)!rainfall interception (mm h2o/s)
  REAL(r8), allocatable :: qintr_snow_c(:,:)!snowfall interception (mm h2o/s)
  REAL(r8), allocatable :: assim_c   (:,:) !canopy assimilation rate (mol m-2 s-1)
  REAL(r8), allocatable :: respc_c   (:,:) !canopy respiration (mol m-2 s-1)

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_1D_PCFluxes
  PUBLIC :: deallocate_1D_PCFluxes
  PUBLIC :: set_1D_PCFluxes

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_1D_PCFluxes
  ! --------------------------------------------------------------------
  ! Allocates memory for CoLM 1d [numpc] variables
  ! --------------------------------------------------------------------

     USE MOD_SPMD_Task
     USE MOD_LandPC
     USE MOD_Precision
     USE MOD_Vars_Global
     IMPLICIT NONE

     IF (p_is_worker) THEN
        IF (numpc > 0) THEN

           allocate (fsenl_c   (0:N_PFT-1,numpc))   ; fsenl_c   (:,:) = spval   ! sensible heat from leaves [W/m2]
           allocate (fevpl_c   (0:N_PFT-1,numpc))   ; fevpl_c   (:,:) = spval   ! evaporation+transpiration from leaves [mm/s]
           allocate (etr_c     (0:N_PFT-1,numpc))   ; etr_c     (:,:) = spval   ! transpiration rate [mm/s]
           allocate (fseng_c   (0:N_PFT-1,numpc))   ; fseng_c   (:,:) = spval   ! sensible heat flux from ground [W/m2]
           allocate (fevpg_c   (0:N_PFT-1,numpc))   ; fevpg_c   (:,:) = spval   ! evaporation heat flux from ground [mm/s]
           allocate (parsun_c  (0:N_PFT-1,numpc))   ; parsun_c  (:,:) = spval   ! solar absorbed by sunlit vegetation [W/m2]
           allocate (parsha_c  (0:N_PFT-1,numpc))   ; parsha_c  (:,:) = spval   ! solar absorbed by shaded vegetation [W/m2]
           allocate (sabvsun_c (0:N_PFT-1,numpc))   ; sabvsun_c (:,:) = spval   ! solar absorbed by sunlit vegetation [W/m2]
           allocate (sabvsha_c (0:N_PFT-1,numpc))   ; sabvsha_c (:,:) = spval   ! solar absorbed by shaded vegetation [W/m2]
           allocate (qintr_c   (0:N_PFT-1,numpc))   ; qintr_c   (:,:) = spval   ! interception (mm h2o/s)
           allocate (qintr_rain_c(0:N_PFT-1,numpc)) ; qintr_rain_c(:,:) = spval ! rainfall interception (mm h2o/s)
           allocate (qintr_snow_c(0:N_PFT-1,numpc)) ; qintr_snow_c(:,:) = spval ! snowfall interception (mm h2o/s)
           allocate (assim_c   (0:N_PFT-1,numpc))   ; assim_c   (:,:) = spval   ! canopy assimilation rate (mol m-2 s-1)
           allocate (respc_c   (0:N_PFT-1,numpc))   ; respc_c   (:,:) = spval   ! canopy respiration (mol m-2 s-1)

        ENDIF
     ENDIF

  END SUBROUTINE allocate_1D_PCFluxes

  SUBROUTINE deallocate_1D_PCFluxes
  ! --------------------------------------------------------------------
  ! deallocates memory for CoLM 1d [numpc] variables
  ! --------------------------------------------------------------------
     USE MOD_SPMD_Task
     USE MOD_LandPC

     IF (p_is_worker) THEN
        IF (numpc > 0) THEN

           deallocate (fsenl_c   )
           deallocate (fevpl_c   )
           deallocate (etr_c     )
           deallocate (fseng_c   )
           deallocate (fevpg_c   )
           deallocate (parsun_c  )
           deallocate (parsha_c  )
           deallocate (sabvsun_c )
           deallocate (sabvsha_c )
           deallocate (qintr_c   )
           deallocate (qintr_rain_c)
           deallocate (qintr_snow_c)
           deallocate (assim_c   )
           deallocate (respc_c   )

        ENDIF
     ENDIF

  END SUBROUTINE deallocate_1D_PCFluxes

  SUBROUTINE set_1D_PCFluxes (Values, Nan)
  ! --------------------------------------------------------------------
  ! Allocates memory for CoLM 1d [numpc] variables
  ! --------------------------------------------------------------------

     USE MOD_SPMD_Task
     USE MOD_LandPC
     USE MOD_Precision
     USE MOD_Vars_Global
     IMPLICIT NONE
     REAL(r8),intent(in) :: Values
     REAL(r8),intent(in) :: Nan

     IF (p_is_worker) THEN
        IF (numpc > 0) THEN

           fsenl_c     (:,:)  = Values ! sensible heat from leaves [W/m2]
           fevpl_c     (:,:)  = Values ! evaporation+transpiration from leaves [mm/s]
           etr_c       (:,:)  = Values ! transpiration rate [mm/s]
           fseng_c     (:,:)  = Values ! sensible heat flux from ground [W/m2]
           fevpg_c     (:,:)  = Values ! evaporation heat flux from ground [mm/s]
           parsun_c    (:,:)  = Values ! solar absorbed by sunlit vegetation [W/m2]
           parsha_c    (:,:)  = Values ! solar absorbed by shaded vegetation [W/m2]
           sabvsun_c   (:,:)  = Values ! solar absorbed by sunlit vegetation [W/m2]
           sabvsha_c   (:,:)  = Values ! solar absorbed by shaded vegetation [W/m2]
           qintr_c     (:,:)  = Values ! interception (mm h2o/s)
           qintr_rain_c(:,:)  = Values ! rainfall interception (mm h2o/s)
           qintr_snow_c(:,:)  = Values ! snowfall interception (mm h2o/s)
           assim_c     (:,:)  = Values ! canopy assimilation rate (mol m-2 s-1)
           respc_c     (:,:)  = Values ! canopy respiration (mol m-2 s-1)

        ENDIF
     ENDIF

  END SUBROUTINE set_1D_PCFluxes

END MODULE MOD_Vars_1DPCFluxes

#endif
! ---------- EOP ------------
