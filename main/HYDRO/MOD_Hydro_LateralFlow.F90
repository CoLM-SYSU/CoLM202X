#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_Hydro_LateralFlow
   !-------------------------------------------------------------------------------------
   ! DESCRIPTION:
   !   
   !   Lateral flow.
   !
   !   Lateral flows in CoLM include
   !   1. Surface flow over hillslopes; 
   !   2. Routing flow in rivers;
   !   3. Groundwater (subsurface) lateral flow. 
   !
   !   Water exchanges between
   !   1. surface flow and rivers;
   !   2. subsurface flow and rivers.
   !
   ! Created by Shupeng Zhang, May 2023
   !-------------------------------------------------------------------------------------
   
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Hydro_RiverNetwork
   USE MOD_Hydro_SubsurfaceNetwork
   USE MOD_Hydro_SurfaceNetwork
   USE MOD_Hydro_SurfaceFlow
   USE MOD_Hydro_SubsurfaceFlow
   USE MOD_Hydro_RiverFlow
   IMPLICIT NONE 

   INTEGER, parameter :: nsubstep = 20

CONTAINS

   ! ----------
   SUBROUTINE lateral_flow_init ()

      IMPLICIT NONE

      CALL surface_network_init    ()
      CALL river_network_init      ()
      CALL subsurface_network_init ()

   END SUBROUTINE lateral_flow_init

   ! ----------
   SUBROUTINE lateral_flow (deltime)

      USE MOD_Mesh,      only : numelm
      USE MOD_LandHRU,   only : numhru
      USE MOD_LandPatch, only : numpatch

      USE MOD_Vars_1DFluxes,      only : rsur
      USE MOD_Vars_TimeVariables, only : wdsrf
      USE MOD_Hydro_Vars_1DFluxes
      USE MOD_Hydro_Vars_TimeVariables

      USE MOD_RangeCheck
      IMPLICIT NONE

      REAL(r8), intent(in) :: deltime

      ! Local Variables
      INTEGER :: nriver
      INTEGER :: istep
      real(r8), allocatable :: wdsrf_p (:)

      IF (p_is_worker) THEN

         nriver = numelm

         IF (nriver > 0) THEN
            riverheight_ta(:) = 0
            rivermomtem_ta(:) = 0
         ENDIF

         IF (numhru > 0) THEN
            wdsrf_hru_ta(:) = 0
            momtm_hru_ta(:) = 0
         ENDIF

         IF (numpatch > 0) THEN
            allocate (wdsrf_p (numpatch))
            wdsrf_p = wdsrf
         ENDIF

         DO istep = 1, nsubstep
            ! (1) Surface flow over hillslopes.
            CALL surface_flow (deltime/nsubstep)
            ! (2) River flow.
            CALL river_flow   (deltime/nsubstep)
         ENDDO

         ! (3) Subsurface lateral flow.
         CALL subsurface_flow (deltime)

         IF (nriver > 0) THEN
            riverheight_ta(:) = riverheight_ta(:) / deltime
            rivermomtem_ta(:) = rivermomtem_ta(:) / deltime

            where (riverheight_ta > 0)
               riverveloct_ta = rivermomtem_ta / riverheight_ta
            ELSE where
               riverveloct_ta = 0
            END where
         ENDIF

         IF (numhru > 0) THEN
            wdsrf_hru_ta(:) = wdsrf_hru_ta(:) / deltime
            momtm_hru_ta(:) = momtm_hru_ta(:) / deltime

            where (wdsrf_hru_ta > 0)
               veloc_hru_ta = momtm_hru_ta / wdsrf_hru_ta
            ELSE where
               veloc_hru_ta = 0.
            END where
         ENDIF

         IF (numpatch > 0) THEN
            rsur(:) = (wdsrf_p(:) - wdsrf(:)) / deltime
         ENDIF

         IF (allocated(wdsrf_p)) deallocate(wdsrf_p)

      ENDIF

#ifdef RangeCheck
      if (p_is_worker .and. (p_iam_worker == 0)) then
         write(*,'(/,A)') 'Checking Lateral Flow Variables ...'
      end if
      CALL check_vector_data ('River Height          ', riverheight)
      CALL check_vector_data ('River Velocity        ', riverveloct)
      CALL check_vector_data ('Surface Water Depth   ', wdsrf_hru)
      CALL check_vector_data ('Surface Water Velocity', veloc_hru)
#endif

   END SUBROUTINE lateral_flow

   ! ----------
   SUBROUTINE lateral_flow_final ()

      IMPLICIT NONE

      CALL surface_network_final    ()
      CALL river_network_final      ()
      CALL subsurface_network_final ()

   END SUBROUTINE lateral_flow_final

END MODULE MOD_Hydro_LateralFlow
#endif
