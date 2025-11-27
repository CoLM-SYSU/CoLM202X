#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_1DFluxes
!-----------------------------------------------------------------------------
! DESCRIPTION:
!   Process fluxes variables for diagnostic for data assimilation
!
! AUTHOR:
!   Lu Li, 07/2025: Initial version
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist, only: DEF_DA_ENS_NUM
   IMPLICIT NONE
   SAVE

   ! public functions
   PUBLIC :: allocate_1D_DAFluxes
   PUBLIC :: deallocate_1D_DAFluxes

   ! define variables
   real(r8), allocatable :: fsena_ens  (:,:) ! sensible heat from canopy height to atmosphere [W/m2]
   real(r8), allocatable :: lfevpa_ens (:,:) ! latent heat flux from canopy height to atmosphere [W/m2]
   real(r8), allocatable :: fevpa_ens  (:,:) ! evapotranspiration from canopy to atmosphere [mm/s]
   real(r8), allocatable :: rsur_ens   (:,:) ! surface runoff (mm h2o/s)

   ! save for analysis increment
   real(r8), allocatable :: fsena_a      (:) ! 
   real(r8), allocatable :: fevpa_a      (:) !
   real(r8), allocatable :: lfevpa_a     (:) ! 
   real(r8), allocatable :: rsur_a       (:) !

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE allocate_1D_DAFluxes ()

!-----------------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Vars_Global
      USE MOD_SPMD_Task
      USE MOD_LandPatch
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate ( fsena_ens  (DEF_DA_ENS_NUM, numpatch) )  ; fsena_ens  (:,:) = spval 
            allocate ( lfevpa_ens (DEF_DA_ENS_NUM, numpatch) )  ; lfevpa_ens (:,:) = spval 
            allocate ( fevpa_ens  (DEF_DA_ENS_NUM, numpatch) )  ; fevpa_ens  (:,:) = spval 
            allocate ( rsur_ens   (DEF_DA_ENS_NUM, numpatch) )  ; rsur_ens   (:,:) = spval 

            allocate ( fsena_a                (numpatch) )  ; fsena_a      (:) = spval
            allocate ( fevpa_a                (numpatch) )  ; fevpa_a      (:) = spval
            allocate ( lfevpa_a               (numpatch) )  ; lfevpa_a     (:) = spval
            allocate ( rsur_a                 (numpatch) )  ; rsur_a       (:) = spval
         ENDIF
      ENDIF

   END SUBROUTINE allocate_1D_DAFluxes

!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_1D_DAFluxes()

!-----------------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_LandPatch

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            deallocate ( fsena_ens  )
            deallocate ( lfevpa_ens )
            deallocate ( fevpa_ens  )
            deallocate ( rsur_ens   )

            deallocate ( fsena_a    )
            deallocate ( fevpa_a    )
            deallocate ( lfevpa_a   )
            deallocate ( rsur_a     )
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_1D_DAFluxes

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_1DFluxes
#endif