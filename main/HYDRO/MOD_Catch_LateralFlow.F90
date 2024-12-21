#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_LateralFlow
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
   USE MOD_ElementNeighbour
   USE MOD_Catch_RiverLakeNetwork
   USE MOD_Catch_HillslopeNetwork
   USE MOD_Catch_HillslopeFlow
   USE MOD_Catch_SubsurfaceFlow
   USE MOD_Catch_RiverLakeFlow
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_Global,    only : dz_soi
   USE MOD_Const_Physical, only : denice, denh2o
   IMPLICIT NONE 

   integer, parameter :: nsubstep = 20
   real(r8) :: dt_average

#ifdef CoLMDEBUG
   real(r8) :: landarea
   real(r8), allocatable :: patcharea (:)  ! m^2
#endif

CONTAINS

   ! ----------
   SUBROUTINE lateral_flow_init (lc_year)

#ifdef CoLMDEBUG
   USE MOD_SPMD_Task
   USE MOD_Mesh
   USE MOD_Pixel
   USE MOD_LandPatch
   USE MOD_Utils
#endif
   IMPLICIT NONE
      
   integer, intent(in) :: lc_year    ! which year of land cover data used

#ifdef CoLMDEBUG
   integer :: ip ,ie, ipxl
#endif

      CALL element_neighbour_init  (lc_year)

      CALL hillslope_network_init  ()
      CALL river_lake_network_init ()
      CALL basin_neighbour_init    ()

#ifdef CoLMDEBUG
      IF (p_is_worker) THEN
         allocate (patcharea (numpatch))
         DO ip = 1, numpatch
            patcharea(ip) = 0.
            ie = landpatch%ielm(ip)
            DO ipxl = landpatch%ipxstt(ip), landpatch%ipxend(ip)
               patcharea(ip) = patcharea(ip) + 1.0e6 * areaquad ( &
                  pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                  pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
            ENDDO
         ENDDO

         landarea = 0.
         IF (numpatch > 0) landarea = sum(patcharea)
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, landarea, 1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
#endif
      ENDIF
#endif

   END SUBROUTINE lateral_flow_init

   ! ----------
   SUBROUTINE lateral_flow (deltime)

   USE MOD_Namelist,  only : DEF_USE_Dynamic_Lake
   USE MOD_Mesh,      only : numelm
   USE MOD_LandHRU,   only : landhru,  numhru,    basin_hru
   USE MOD_LandPatch, only : numpatch, elm_patch, hru_patch

   USE MOD_Vars_Global,         only : nl_lake
   USE MOD_Const_Physical,      only : tfrz
   USE MOD_Vars_1DFluxes,       only : rsur, rsub, rnof
   USE MOD_Vars_TimeVariables,  only : wdsrf, t_lake, lake_icefrac, t_soisno
   USE MOD_Vars_TimeInvariants, only : lakedepth, dz_lake
   USE MOD_Hydro_Vars_1DFluxes
   USE MOD_Hydro_Vars_TimeVariables

   USE MOD_Lake, only : adjust_lake_layer

   USE MOD_RangeCheck
   IMPLICIT NONE

   real(r8), intent(in) :: deltime

   ! Local Variables
   integer  :: numbasin, ibasin, ihru, i, j, ps, pe, istep
   real(r8), allocatable :: wdsrf_p (:)
#ifdef CoLMDEBUG
   real(r8) :: dtolw, toldis
#endif

      IF (p_is_worker) THEN

         numbasin = numelm

         ! a) The smallest unit in surface lateral flow (including hillslope flow and river-lake flow)
         !    is HRU and the main prognostic variable is "wdsrf_hru" (surface water depth).
         ! b) "wdsrf_hru" is updated by aggregating water depths in patches.
         ! c) Water surface in a basin ("wdsrf_bsn", defined as the lowest surface water in the basin) 
         ! is derived from "wdsrf_hru".
         DO i = 1, numhru
            ps = hru_patch%substt(i)
            pe = hru_patch%subend(i)
            wdsrf_hru(i) = sum(wdsrf(ps:pe) * hru_patch%subfrc(ps:pe))
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

         dt_average = 0.

         IF (numpatch > 0) rsur     (:) = 0.
         IF (numbasin > 0) discharge(:) = 0.

         DO istep = 1, nsubstep

            ! (1) Surface flow over hillslopes.
            CALL hillslope_flow (deltime/nsubstep)
         
            ! (2) River and Lake flow.
            CALL river_lake_flow (deltime/nsubstep)
      
            dt_average = dt_average + deltime/nsubstep/ntimestep_riverlake
         
         ENDDO
         
         IF (numpatch > 0) rsur = rsur / deltime
         IF (numbasin > 0) discharge = discharge / deltime

         IF (numbasin > 0) THEN
            wdsrf_bsn_ta(:) = wdsrf_bsn_ta(:) / deltime
            momen_riv_ta(:) = momen_riv_ta(:) / deltime

            WHERE (wdsrf_bsn_ta > 0)
               veloc_riv_ta = momen_riv_ta / wdsrf_bsn_ta
            ELSE WHERE
               veloc_riv_ta = 0
            END WHERE
         ENDIF

         IF (numhru > 0) THEN
            wdsrf_hru_ta(:) = wdsrf_hru_ta(:) / deltime
            momen_hru_ta(:) = momen_hru_ta(:) / deltime

            WHERE (wdsrf_hru_ta > 0)
               veloc_hru_ta = momen_hru_ta / wdsrf_hru_ta
            ELSE WHERE
               veloc_hru_ta = 0.
            END WHERE
         ENDIF

         ! update surface water depth on patches
         DO i = 1, numhru
            ps = hru_patch%substt(i)
            pe = hru_patch%subend(i)
            wdsrf(ps:pe) = wdsrf_hru(i) * 1.0e3 ! m to mm
         ENDDO
            
         IF (numpatch > 0) THEN
            xwsur(:) = (wdsrf_p(:) - wdsrf(:)) / deltime
         ENDIF

         ! (3) Subsurface lateral flow.
         CALL subsurface_flow (deltime)
         
         IF (numpatch > 0) THEN
            rnof(:) = rsur(:) + rsub(:)
         ENDIF
       
         DO i = 1, numpatch
            h2osoi(:,i) = wliq_soisno(1:,i)/(dz_soi(1:)*denh2o) + wice_soisno(1:,i)/(dz_soi(1:)*denice)
            wat(i)      = sum(wice_soisno(1:,i)+wliq_soisno(1:,i)) + ldew(i) + scv(i) + wetwat(i)
         ENDDO

         ! (4) vertical layers adjustment
         IF (DEF_USE_Dynamic_Lake) THEN
            DO i = 1, numpatch
               IF (wdsrf_p(i) >= 100.) THEN
                  ! wet previously
                  dz_lake(:,i) = dz_lake(:,i) * wdsrf(i)*1.e-3/sum(dz_lake(:,i))
               ELSE
                  ! dry previously
                  dz_lake(:,i) = wdsrf(i)*1.e-3/nl_lake
                  t_lake (:,i) = t_soisno(1,i)
                  IF (t_soisno(1,i) >= tfrz) THEN
                     lake_icefrac(:,i) = 0.
                  ELSE
                     lake_icefrac(:,i) = 1.
                  ENDIF
               ENDIF
                  
               IF (wdsrf(i) >= 100.) THEN
                  CALL adjust_lake_layer (nl_lake, dz_lake(:,i), t_lake(:,i), lake_icefrac(:,i))
               ENDIF
            ENDDO
         ENDIF

         IF (allocated(wdsrf_p)) deallocate(wdsrf_p)

      ENDIF

#ifdef RangeCheck
      IF (p_is_worker .and. (p_iam_worker == 0)) THEN
         write(*,'(/,A)') 'Checking Lateral Flow Variables ...'
         write(*,'(A,F12.5,A)') 'River Lake Flow average timestep: ', &
               dt_average/nsubstep, ' seconds'
      ENDIF

      CALL check_vector_data ('Basin Water Depth   [m]  ', wdsrf_bsn)
      CALL check_vector_data ('River Velocity      [m/s]', veloc_riv)
      CALL check_vector_data ('HRU Water Depth     [m]  ', wdsrf_hru)
      CALL check_vector_data ('HRU Water Velocity  [m/s]', veloc_hru)
      CALL check_vector_data ('Subsurface bt basin [m/s]', xsubs_bsn)
      CALL check_vector_data ('Subsurface bt HRU   [m/s]', xsubs_hru)
      CALL check_vector_data ('Subsurface bt patch [m/s]', xsubs_pch)

#ifdef CoLMDEBUG
      IF (p_is_worker) THEN
         
         dtolw  = 0
         toldis = 0
        
         IF (numpatch > 0) THEN
            dtolw = sum(patcharea * xwsur) / 1.e3 * deltime
         ENDIF
         IF (numelm > 0) THEN
            toldis = sum(discharge*deltime, mask = (riverdown == 0) .or. (riverdown == -3))
            dtolw  = dtolw - toldis
         ENDIF

#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, dtolw,  1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, toldis, 1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
#endif
         IF (p_iam_worker == 0) THEN
            write(*,'(A,F10.2,A,ES10.3,A,ES10.3,A)') 'Total surface water error: ', dtolw, &
               '(m^3) in area ', landarea, '(m^2), discharge ', toldis, '(m^3)' 
         ENDIF

         dtolw = 0
         IF (numpatch > 0) dtolw = sum(patcharea * xwsub) / 1.e3 * deltime
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, dtolw,  1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
#endif
         IF (p_iam_worker == 0) THEN
            write(*,'(A,F10.2,A,ES10.3,A)') 'Total ground  water error: ', dtolw, &
               '(m^3) in area ', landarea, '(m^2)'
         ENDIF
      ENDIF
#endif
#endif

   END SUBROUTINE lateral_flow

   ! ----------
   SUBROUTINE lateral_flow_final ()

   IMPLICIT NONE

      CALL hillslope_network_final  ()
      CALL river_lake_network_final ()
      CALL basin_neighbour_final    ()

#ifdef CoLMDEBUG
      IF (allocated(patcharea)) deallocate(patcharea)
#endif

   END SUBROUTINE lateral_flow_final

END MODULE MOD_Catch_LateralFlow
#endif
