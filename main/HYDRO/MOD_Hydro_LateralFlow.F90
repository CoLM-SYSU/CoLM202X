#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_Hydro_LateralFlow
   
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

      CALL surface_network_init ()
      CALL river_init ()
      CALL ssrf_init ()

   END SUBROUTINE lateral_flow_init

   ! ----------
   SUBROUTINE lateral_flow (deltime)

      USE MOD_Mesh,      only : numelm
      USE MOD_LandHRU,   only : numhru
      USE MOD_LandPatch, only : numpatch
      USE MOD_Vars_1DFluxes
      USE MOD_Hydro_Vars_1DFluxes
      USE MOD_Hydro_Vars_Timevars
      USE MOD_CoLMDebug
      IMPLICIT NONE

      REAL(r8), intent(in) :: deltime

      ! Local Variables
      INTEGER :: nriver
      INTEGER :: istep

      IF (p_is_worker) THEN

         nriver = numelm

         IF (nriver > 0) THEN
            riverheight_ta(:) = 0
            rivermomtem_ta(:) = 0
         ENDIF

         IF (numhru > 0) THEN
            dpond_hru_ta(:) = 0
            momtm_hru_ta(:) = 0
            rsurf_hru(:) = 0
         ENDIF

         IF (numpatch > 0) THEN
            rsur(:) = 0
         ENDIF

         DO istep = 1, nsubstep

            CALL surface_runoff (deltime/nsubstep)

            CALL river_flow     (deltime/nsubstep)

         ENDDO

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
            dpond_hru_ta(:) = dpond_hru_ta(:) / deltime
            momtm_hru_ta(:) = momtm_hru_ta(:) / deltime

            where (dpond_hru_ta > 0)
               veloc_hru_ta = momtm_hru_ta / dpond_hru_ta
            ELSE where
               veloc_hru_ta = 0.
            END where

            rsurf_hru(:) = rsurf_hru(:) / deltime
         ENDIF

         IF (numpatch > 0) THEN
            rsur(:) = rsur(:) / deltime
         ENDIF

         CALL subsurface_runoff (deltime)

      ENDIF

#ifdef CoLMDEBUG
      if (p_is_worker .and. (p_iam_worker == 0)) then
         write(*,'(/,A)') 'Checking Lateral Flow Variables ...'
      end if
      CALL check_vector_data ('River Height          ', riverheight)
      CALL check_vector_data ('River Height time ave ', riverheight_ta)
      CALL check_vector_data ('River Velocity        ', riverveloct)
      CALL check_vector_data ('Surface Water Depth   ', dpond_hru)
      CALL check_vector_data ('Surface Water Velocity', veloc_hru)
      CALL check_vector_data ('Surface Runoff        ', rsur)
      CALL check_vector_data ('Depth to Water Table  ', zwt_hru)
      CALL check_vector_data ('Subsurface runoff     ', rsubs_pch)

#endif

   END SUBROUTINE lateral_flow

   ! ----------
   SUBROUTINE lateral_flow_final ()

      IMPLICIT NONE

      CALL surface_network_final ()
      CALL river_final ()
      CALL ssrf_final  ()

   END SUBROUTINE lateral_flow_final

   ! ----------
   SUBROUTINE check_catchment_water_inout ()
      
      USE MOD_Mesh
      USE MOD_LandPatch
      USE MOD_Vars_Global, only : nl_soil
      USE MOD_Vars_1DForcing,    only : forc_prc,    forc_prl
      USE MOD_Vars_1DFluxes,     only : fevpa,       rnof,        rsur
      USE MOD_Vars_TimeVariables, only : wliq_soisno, wice_soisno, ldew, scv, wa
      IMPLICIT NONE

      INTEGER , pointer :: ihru_p  (:)
      REAL(r8), pointer :: area_p  (:)

      REAL(r8) :: area_all , area_cat
      REAL(r8) :: precp_all, precp_cat
      REAL(r8) :: ldew_all , ldew_cat
      REAL(r8) :: scv_all  , scv_cat
      REAL(r8) :: evpa_all , evpa_cat
      REAL(r8) :: rsur_all , rsur_cat
      REAL(r8) :: rnof_all , rnof_cat
      REAL(r8) :: wsoi_all , wsoi_cat
      REAL(r8) :: wa_all   , wa_cat

      REAL(r8) :: wsoi_l (nl_soil)

      INTEGER :: numbasin, nhru, ibasin, ih, istt, iend, ilev

      area_all  = 0.
      precp_all = 0.
      ldew_all  = 0.
      scv_all   = 0.
      evpa_all  = 0.
      rsur_all  = 0.
      rnof_all  = 0.
      wsoi_all  = 0.
      wa_all    = 0.

      IF (p_is_worker) THEN

         numbasin = numelm

         DO ibasin = 1, numbasin

            nhru = surface_network(ibasin)%nhru
            ihru_p => surface_network(ibasin)%ihru
            area_p => surface_network(ibasin)%area

            area_cat  = 0.
            precp_cat = 0.
            ldew_cat  = 0.
            scv_cat   = 0.
            evpa_cat  = 0.
            rsur_cat  = 0.
            rnof_cat  = 0.
            wsoi_cat  = 0.
            wa_cat    = 0.

            DO ih = 1, nhru
               istt = hru_patch%substt(ihru_p(ih))
               iend = hru_patch%subend(ihru_p(ih))

               area_cat  = area_cat + area_p(ih)

               precp_cat = precp_cat + sum(forc_prc(istt:iend) * hru_patch%subfrc(istt:iend)) * area_p(ih) &
                  + sum(forc_prl(istt:iend) * hru_patch%subfrc(istt:iend)) * area_p(ih)
               ldew_cat = ldew_cat + sum(ldew (istt:iend) * hru_patch%subfrc(istt:iend)) * area_p(ih)
               scv_cat  = scv_cat  + sum(scv  (istt:iend) * hru_patch%subfrc(istt:iend)) * area_p(ih)
               evpa_cat = evpa_cat + sum(fevpa(istt:iend) * hru_patch%subfrc(istt:iend)) * area_p(ih)
               rsur_cat = rsur_cat + sum(rsur (istt:iend) * hru_patch%subfrc(istt:iend)) * area_p(ih)
               rnof_cat = rnof_cat + sum(rnof (istt:iend) * hru_patch%subfrc(istt:iend)) * area_p(ih)
               wa_cat   = wa_cat   + sum(wa   (istt:iend) * hru_patch%subfrc(istt:iend)) * area_p(ih)

               DO ilev = 1, nl_soil
                  wsoi_l(ilev) = sum(wliq_soisno(ilev,istt:iend) * hru_patch%subfrc(istt:iend)) &
                     + sum(wice_soisno(ilev,istt:iend) * hru_patch%subfrc(istt:iend))
               ENDDO
               wsoi_cat = wsoi_cat + sum(wsoi_l) * area_p(ih)
            ENDDO

            ! write(*,100) landbasin(ibasin)%indx, precp_cat, evpa_cat, rsur_cat, rnof_cat-rsur_cat

         ENDDO
      ENDIF

      100 format('CAT ', I8.8, ', Precp ', E12.4, ', Evpa ', E12.4, ', Rsurf ', E12.4, ', Rsub ', E12.4)

   END SUBROUTINE check_catchment_water_inout

END MODULE MOD_Hydro_LateralFlow
#endif
