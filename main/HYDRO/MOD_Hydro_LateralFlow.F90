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
   USE MOD_Hydro_Vars_TimeVariables
   USE MOD_Hydro_RiverLakeNetwork
   USE MOD_Hydro_BasinNeighbour
   USE MOD_Hydro_HillslopeNetwork
   USE MOD_Hydro_HillslopeFlow
   USE MOD_Hydro_SubsurfaceFlow
   USE MOD_Hydro_RiverLakeFlow
   IMPLICIT NONE 

   INTEGER, parameter :: nsubstep = 20

CONTAINS

   ! ----------
   SUBROUTINE lateral_flow_init ()

      IMPLICIT NONE

      CALL hillslope_network_init  ()
      CALL river_lake_network_init ()
      CALL basin_neighbour_init    ()

      IF (p_is_worker) THEN
         wdsrf_bsn_prev(:) = wdsrf_bsn(:)
         wdsrf_hru_prev(:) = wdsrf_hru(:)
      ENDIF

   END SUBROUTINE lateral_flow_init

   ! ----------
   SUBROUTINE lateral_flow (deltime)

      USE MOD_Mesh,      only : numelm
      USE MOD_LandHRU,   only : landhru,  numhru,    basin_hru
      USE MOD_LandPatch, only : numpatch, elm_patch, hru_patch

      USE MOD_Vars_1DFluxes,       only : rsur, rsub, rnof
      USE MOD_Vars_TimeVariables,  only : wdsrf
      USE MOD_Vars_TimeInvariants, only : lakedepth
      USE MOD_Hydro_Vars_1DFluxes
      USE MOD_Hydro_Vars_TimeVariables

      USE MOD_RangeCheck
      IMPLICIT NONE

      REAL(r8), intent(in) :: deltime

      ! Local Variables
      INTEGER  :: nbasin, ibasin, ihru, i, j, istt, iend, istep
      real(r8), allocatable :: wdsrf_p (:)

      IF (p_is_worker) THEN

         nbasin = numelm

         ! a) The smallest unit in surface lateral flow (including hillslope flow and river-lake flow)
         !    is HRU and the main prognostic variable is "wdsrf_hru" (surface water depth).
         ! b) "wdsrf_hru" is updated by aggregating water depths in patches.
         ! c) Water surface in a basin ("wdsrf_bsn", defined as the lowest surface water in the basin) 
         ! is derived from "wdsrf_hru".
         DO i = 1, numhru
            istt = hru_patch%substt(i)
            iend = hru_patch%subend(i)
            wdsrf_hru(i) = sum(wdsrf(istt:iend) * hru_patch%subfrc(istt:iend))
            wdsrf_hru(i) = wdsrf_hru(i) / 1.0e3 ! mm to m
         ENDDO

         wdsrf_hru_ta(:) = 0
         momen_hru_ta(:) = 0
         wdsrf_bsn_ta(:) = 0
         momen_riv_ta(:) = 0

         IF (numpatch > 0) THEN
            allocate (wdsrf_p (numpatch))
            wdsrf_p = wdsrf
         ENDIF

         DO istep = 1, nsubstep

            ! (1) Surface flow over hillslopes.
            CALL hillslope_flow (deltime/nsubstep)
         
            ! (2) River and Lake flow.
            CALL river_lake_flow (deltime/nsubstep)
         
         ENDDO

         IF (nbasin > 0) THEN
            wdsrf_bsn_ta(:) = wdsrf_bsn_ta(:) / deltime
            momen_riv_ta(:) = momen_riv_ta(:) / deltime

            where (wdsrf_bsn_ta > 0)
               veloc_riv_ta = momen_riv_ta / wdsrf_bsn_ta
            ELSE where
               veloc_riv_ta = 0
            END where
         ENDIF

         IF (numhru > 0) THEN
            wdsrf_hru_ta(:) = wdsrf_hru_ta(:) / deltime
            momen_hru_ta(:) = momen_hru_ta(:) / deltime

            where (wdsrf_hru_ta > 0)
               veloc_hru_ta = momen_hru_ta / wdsrf_hru_ta
            ELSE where
               veloc_hru_ta = 0.
            END where
         ENDIF

         ! update surface water depth on patches
         DO i = 1, numhru
            istt = hru_patch%substt(i)
            iend = hru_patch%subend(i)
            wdsrf(istt:iend) = wdsrf_hru(i) * 1.0e3 ! m to mm
         ENDDO
            
         IF (numpatch > 0) THEN
            rsur(:) = (wdsrf_p(:) - wdsrf(:)) / deltime
         ENDIF

         IF (allocated(wdsrf_p)) deallocate(wdsrf_p)

         ! (3) Subsurface lateral flow.
         CALL subsurface_flow (deltime)
         
         IF (numpatch > 0) THEN
            rnof(:) = rsur(:) + rsub(:)
         ENDIF

      ENDIF

#ifdef RangeCheck
      if (p_is_worker .and. (p_iam_worker == 0)) then
         write(*,'(/,A)') 'Checking Lateral Flow Variables ...'
      end if
      CALL check_vector_data ('Basin Water Depth   [m]  ', wdsrf_bsn)
      CALL check_vector_data ('River Velocity      [m/s]', veloc_riv)
      CALL check_vector_data ('HRU Water Depth     [m]  ', wdsrf_hru)
      CALL check_vector_data ('HRU Water Velocity  [m/s]', veloc_hru)
      CALL check_vector_data ('Subsurface bt basin [m/s]', rsubs_bsn)
      CALL check_vector_data ('Subsurface bt HRU   [m/s]', rsubs_hru)
      CALL check_vector_data ('Subsurface bt patch [m/s]', rsubs_pch)
#endif

   END SUBROUTINE lateral_flow

   ! ----------
   SUBROUTINE lateral_flow_final ()

      IMPLICIT NONE

      CALL hillslope_network_final  ()
      CALL river_lake_network_final ()
      CALL basin_neighbour_final    ()

   END SUBROUTINE lateral_flow_final

END MODULE MOD_Hydro_LateralFlow
#endif
