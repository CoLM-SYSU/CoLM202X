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
   USE MOD_Catch_Vars_TimeVariables
   USE MOD_ElementNeighbour
   USE MOD_Catch_BasinNetwork
   USE MOD_Catch_RiverLakeNetwork
   USE MOD_Catch_HillslopeNetwork
   USE MOD_Catch_HillslopeFlow
   USE MOD_Catch_SubsurfaceFlow
   USE MOD_Catch_RiverLakeFlow
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_Global,    only: dz_soi
   USE MOD_Const_Physical, only: denice, denh2o
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
   USE MOD_Catch_WriteParameters
   IMPLICIT NONE

   integer, intent(in) :: lc_year    ! which year of land cover data used

#ifdef CoLMDEBUG
   integer :: ip ,ie, ipxl
#endif

      CALL element_neighbour_init  (lc_year)
      CALL river_lake_network_init ()
      CALL subsurface_network_init ()

      CALL write_catch_parameters ()

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

   USE MOD_Namelist,  only: DEF_USE_Dynamic_Lake
   USE MOD_Mesh,      only: numelm
   USE MOD_LandHRU,   only: landhru,  numhru,    elm_hru
   USE MOD_LandPatch, only: numpatch, elm_patch, hru_patch

   USE MOD_Vars_Global,         only: nl_lake
   USE MOD_Const_Physical,      only: tfrz
   USE MOD_Vars_TimeVariables,  only: wdsrf, t_lake, lake_icefrac, t_soisno
   USE MOD_Vars_TimeInvariants, only: lakedepth, dz_lake
   USE MOD_Vars_1DFluxes,       only: rsur, rsub, rnof
   USE MOD_Catch_Vars_1DFluxes
   USE MOD_Catch_Vars_TimeVariables
   USE MOD_Catch_RiverLakeNetwork

   USE MOD_Lake, only: adjust_lake_layer

   USE MOD_UserDefFun, only : findloc_ud
   USE MOD_RangeCheck
   IMPLICIT NONE

   real(r8), intent(in) :: deltime

   ! Local Variables
   integer  :: i, j, j0, h, ps, pe, istep, s
   real(r8) :: rnofsrf, sumarea
   real(r8), allocatable :: wdsrf_p (:), wdsrf_hru_p (:)
#ifdef CoLMDEBUG
   real(r8) :: dtolw, toldis
#endif

      IF (p_is_worker) THEN

         ! a) The smallest unit in surface lateral flow (including hillslope flow and river-lake flow)
         !    is HRU and the main prognostic variable is "wdsrf_bsnhru" (surface water depth).
         ! b) "wdsrf_bsnhru" is updated by aggregating water depths in patches.
         ! c) Water surface in a basin ("wdsrf_bsn", defined as the lowest surface water in the basin)
         !    is derived from "wdsrf_bsnhru".
         DO i = 1, numhru
            ps = hru_patch%substt(i)
            pe = hru_patch%subend(i)
            wdsrf_hru(i) = sum(wdsrf(ps:pe) * hru_patch%subfrc(ps:pe))
            wdsrf_hru(i) = wdsrf_hru(i) / 1.0e3 ! mm to m
         ENDDO

         IF (numhru > 0) THEN
            allocate (wdsrf_hru_p (numhru))
            wdsrf_hru_p = wdsrf_hru
         ENDIF

         IF (numpatch > 0) THEN
            allocate (wdsrf_p (numpatch))
            wdsrf_p = wdsrf
         ENDIF

         dt_average = 0.

         IF (numbasin > 0)  wdsrf_bsn_ta    (:) = 0.
         IF (numbasin > 0)  momen_riv_ta    (:) = 0.
         IF (numbasin > 0)  discharge_ta    (:) = 0.
         IF (numbsnhru > 0) wdsrf_bsnhru_ta (:) = 0.
         IF (numbsnhru > 0) momen_bsnhru_ta (:) = 0.

         CALL worker_push_subset_data (iam_elm, iam_bsn, elm_hru, basin_hru, wdsrf_hru, wdsrf_bsnhru)

         DO istep = 1, nsubstep

            ! (1) ------------------- Surface flow over hillslopes. -------------------
            CALL hillslope_flow (deltime/nsubstep)

            ! (2) ----------------------- River and Lake flow. ------------------------
            CALL river_lake_flow (deltime/nsubstep)

            dt_average = dt_average + deltime/nsubstep/ntimestep_riverlake

         ENDDO

         IF (numbasin > 0) THEN
            wdsrf_bsn_ta(:) = wdsrf_bsn_ta(:) / deltime
            momen_riv_ta(:) = momen_riv_ta(:) / deltime

            WHERE (wdsrf_bsn_ta > 0)
               veloc_riv_ta = momen_riv_ta / wdsrf_bsn_ta
            ELSE WHERE
               veloc_riv_ta = 0
            END WHERE

            discharge_ta = discharge_ta / deltime
         ENDIF

         IF (numbsnhru > 0) THEN
            wdsrf_bsnhru_ta(:) = wdsrf_bsnhru_ta(:) / deltime
            momen_bsnhru_ta(:) = momen_bsnhru_ta(:) / deltime

            WHERE (wdsrf_bsnhru_ta > 0)
               veloc_bsnhru_ta = momen_bsnhru_ta / wdsrf_bsnhru_ta
            ELSE WHERE
               veloc_bsnhru_ta = 0.
            END WHERE
         ENDIF

         ! update surface water depth on patches
         CALL worker_push_subset_data (iam_bsn, iam_elm, basin_hru, elm_hru, wdsrf_bsnhru, wdsrf_hru)
         DO i = 1, numhru
            wdsrf_hru(i) = max(0., wdsrf_hru(i))
            ps = hru_patch%substt(i)
            pe = hru_patch%subend(i)
            wdsrf(ps:pe) = wdsrf_hru(i) * 1.0e3 ! m to mm
         ENDDO

         IF (numpatch > 0) THEN
            xwsur(:) = (wdsrf_p(:) - wdsrf(:)) / deltime
         ENDIF

         ! update surface runoff from hillslope to river
         IF (numpatch > 0) rsur(:) = 0.
         DO i = 1, numelm
            IF (lake_id_elm(i) <= 0) THEN

               rnofsrf = 0.
               sumarea = 0.

               IF (lake_id_elm(i) == 0) THEN
                  j0 = 2 ! regular catchment with river
               ELSE
                  j0 = 1 ! catchment directly to lake
               ENDIF

               DO j = j0, hillslope_element(i)%nhru
                  h = hillslope_element(i)%ihru(j)
                  rnofsrf = rnofsrf + (wdsrf_hru_p(h) - wdsrf_hru(h)) * hillslope_element(i)%area(j)
                  sumarea = sumarea + hillslope_element(i)%area(j)
               ENDDO

               IF (sumarea > 0.) THEN
                  rnofsrf = rnofsrf / sumarea * 1.e3 / deltime ! unit: mm/s
                  DO j = j0, hillslope_element(i)%nhru
                     h  = hillslope_element(i)%ihru(j)
                     ps = hru_patch%substt(h)
                     pe = hru_patch%subend(h)
                     rsur(ps:pe) = rnofsrf
                  ENDDO
               ENDIF
            ENDIF
         ENDDO

         ! update fraction of flooded area
         DO i = 1, numelm
            IF (lake_id_elm(i) <= 0) THEN
               DO j = 1, hillslope_element(i)%nhru
                  h  = hillslope_element(i)%ihru(j)
                  ps = hru_patch%substt(h)
                  pe = hru_patch%subend(h)

                  IF (hillslope_element(i)%indx(j) == 0) THEN
                     fldarea(ps:pe) = 1.0 ! river
                  ELSE
                     s = findloc_ud(hillslope_element(i)%fldprof(:,j) <= wdsrf_hru(h), back = .true.)
                     IF (s == nfldstep) THEN
                        fldarea(ps:pe) = 1.0
                     ELSEIF (s == 0) THEN
                        fldarea(ps:pe) = sqrt(wdsrf_hru(h)/hillslope_element(i)%fldprof(1,j) * 1./nfldstep**2)
                     ELSE
                        fldarea(ps:pe) = sqrt( (wdsrf_hru(h)-hillslope_element(i)%fldprof(s,j))            &
                           / (hillslope_element(i)%fldprof(s+1,j)-hillslope_element(i)%fldprof(s,j)) &
                           * real(2*s+1) / nfldstep**2 + real(s**2)/nfldstep**2 )
                     ENDIF
                  ENDIF
               ENDDO
            ELSE
               ps = elm_patch%substt(i)
               pe = elm_patch%subend(i)
               fldarea(ps:pe) = 1.0 ! lake
            ENDIF
         ENDDO

         ! (3) ------------------- Subsurface lateral flow -------------------
         CALL subsurface_flow (deltime)

         DO i = 1, numpatch
            h2osoi(:,i) = wliq_soisno(1:,i)/(dz_soi(1:)*denh2o) + wice_soisno(1:,i)/(dz_soi(1:)*denice)
            wat(i)      = sum(wice_soisno(1:,i)+wliq_soisno(1:,i)) + ldew(i) + scv(i) + wetwat(i)
         ENDDO

         IF (numpatch > 0) rnof = rsur + rsub

         ! (4) ---------------- vertical layers adjustment ---------------------
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

         IF (allocated(wdsrf_p    )) deallocate(wdsrf_p    )
         IF (allocated(wdsrf_hru_p)) deallocate(wdsrf_hru_p)

      ENDIF

#ifdef RangeCheck
      IF (p_is_worker .and. (p_iam_worker == 0)) THEN
         write(*,'(/,A)') 'Checking Lateral Flow Variables ...'
         write(*,'(A,F12.5,A)') 'River Lake Flow minimum average timestep: ', &
               dt_average/nsubstep, ' seconds'
      ENDIF

      CALL check_vector_data ('Basin Water Depth   [m]  ', wdsrf_bsn)
      CALL check_vector_data ('River Velocity      [m/s]', veloc_riv)
      CALL check_vector_data ('HRU Water Depth     [m]  ', wdsrf_bsnhru)
      CALL check_vector_data ('HRU Water Velocity  [m/s]', veloc_bsnhru)
      CALL check_vector_data ('Subsurface bt basin [m/s]', xsubs_elm)
      CALL check_vector_data ('Subsurface bt HRU   [m/s]', xsubs_hru)
      CALL check_vector_data ('Subsurface bt patch [m/s]', xsubs_pch)

#ifdef CoLMDEBUG
      IF (p_is_worker) THEN

         dtolw  = 0
         toldis = 0

         IF (numpatch > 0) THEN
            dtolw = sum(patcharea * xwsur) / 1.e3 * deltime
         ENDIF
         IF (numbasin > 0) THEN
            toldis = sum(discharge_ta*deltime, mask = (riverdown == 0) .or. (riverdown == -3))
            dtolw  = dtolw - toldis
         ENDIF

#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, dtolw,  1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, toldis, 1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
#endif
         IF (p_iam_worker == 0) THEN
            write(*,'(A,F12.2,A,ES8.1,A,ES10.3,A)') 'Total surface water error: ', dtolw, &
               '(m^3) in area ', landarea, '(m^2), discharge ', toldis, '(m^3)'
         ENDIF

         dtolw = 0
         IF (numpatch > 0) dtolw = sum(patcharea * xwsub) / 1.e3 * deltime
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, dtolw,  1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
#endif
         IF (p_iam_worker == 0) THEN
            write(*,'(A,F12.2,A,ES8.1,A)') 'Total ground  water error: ', dtolw, &
               '(m^3) in area ', landarea, '(m^2)'
         ENDIF
      ENDIF
#endif
#endif

   END SUBROUTINE lateral_flow

   ! ----------
   SUBROUTINE lateral_flow_final ()

   IMPLICIT NONE

      CALL river_lake_network_final ()
      CALL subsurface_network_final ()
      CALL basin_network_final      ()

#ifdef CoLMDEBUG
      IF (allocated(patcharea)) deallocate(patcharea)
#endif

   END SUBROUTINE lateral_flow_final

END MODULE MOD_Catch_LateralFlow
#endif
