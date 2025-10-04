#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_RiverLakeFlow
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Shallow water equation solver in rivers.
!
!   References
!   [1] Toro EF. Shock-capturing methods for free-surface shallow flows.
!      Chichester: John Wiley & Sons; 2001.
!   [2] Liang, Q., Borthwick, A. G. L. (2009). Adaptive quadtree simulation of shallow
!      flows with wet-dry fronts over complex topography.
!      Computers and Fluids, 38(2), 221-234.
!   [3] Audusse, E., Bouchut, F., Bristeau, M.-O., Klein, R., Perthame, B. (2004).
!      A Fast and Stable Well-Balanced Scheme with Hydrostatic Reconstruction for
!      Shallow Water Flows. SIAM Journal on Scientific Computing, 25(6), 2050-2065.
!
! Created by Shupeng Zhang, May 2023
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE

   real(r8), parameter :: nmanning_riv = 0.03

   real(r8), parameter :: RIVERMIN  = 1.e-5_r8
   real(r8), parameter :: VOLUMEMIN = 1.e-5_r8

   integer :: ntimestep_riverlake

CONTAINS

   ! ---------
   SUBROUTINE river_lake_flow (year, dt)

   USE MOD_SPMD_Task
   USE MOD_Namelist, only : DEF_Reservoir_Method
   USE MOD_Utils
   USE MOD_Catch_BasinNetwork
   USE MOD_Catch_HillslopeNetwork
   USE MOD_Catch_RiverLakeNetwork
   USE MOD_Catch_Reservoir
   USE MOD_Catch_Vars_TimeVariables
   USE MOD_Catch_Vars_1DFluxes
   USE MOD_Const_Physical, only: grav
   IMPLICIT NONE

   integer,  intent(in) :: year
   real(r8), intent(in) :: dt

   ! Local Variables
   integer  :: hs, he, i, j, irsv
   real(r8) :: dt_this

   real(r8),  allocatable :: wdsrf_bsn_ds(:)
   real(r8),  allocatable :: veloc_riv_ds(:)
   real(r8),  allocatable :: momen_riv_ds(:)

   real(r16), allocatable :: hflux_fc(:)
   real(r16), allocatable :: mflux_fc(:)
   real(r16), allocatable :: zgrad_dn(:)

   real(r16), allocatable :: hflux_resv(:)
   real(r16), allocatable :: mflux_resv(:)

   real(r16), allocatable :: sum_hflux_riv(:)
   real(r16), allocatable :: sum_mflux_riv(:)
   real(r16), allocatable :: sum_zgrad_riv(:)

   real(r8) :: veloct_fc, height_fc, momen_fc, zsurf_fc
   real(r8) :: bedelv_fc, height_up, height_dn
   real(r8) :: vwave_up, vwave_dn, hflux_up, hflux_dn, mflux_up, mflux_dn
   real(r8) :: totalvolume, friction, dvol, nextl, nexta, nextv, ddep
   real(r8),  allocatable :: dt_res(:), dt_all(:), dt_tmp(:)
   logical,   allocatable :: hmask(:), bsnfilter(:)


      IF (p_is_worker) THEN

         ! update water depth in basin by aggregating water depths in patches
         DO i = 1, numbasin
            hs = basin_hru%substt(i)
            he = basin_hru%subend(i)

            IF (lake_id(i) <= 0) THEN
               ! river or lake catchment
               ! Water surface in a basin is defined as the lowest surface water in the basin
               wdsrf_bsn(i) = minval(hillslope_basin(i)%hand + wdsrf_bsnhru(hs:he)) - handmin(i)

            ELSEIF (lake_id(i) > 0) THEN
               ! lake
               totalvolume  = sum(wdsrf_bsnhru(hs:he) * lakeinfo(i)%area0)
               wdsrf_bsn(i) = lakeinfo(i)%surface(totalvolume)
            ENDIF


            IF (lake_id(i) == 0) THEN
               ! river momentum is less or equal than the momentum at last time step.
               IF (wdsrf_bsn_prev(i) < wdsrf_bsn(i)) THEN
                  momen_riv(i) = wdsrf_bsn_prev(i) * veloc_riv(i)
                  veloc_riv(i) = momen_riv(i) / wdsrf_bsn(i)
               ELSE
                  momen_riv(i) = wdsrf_bsn(i) * veloc_riv(i)
               ENDIF
            ELSE
               ! water in lake or lake catchment is assumued to be stationary.
               ! TODO: lake dynamics
               momen_riv(i) = 0
               veloc_riv(i) = 0
            ENDIF

         ENDDO

         IF (numbasin > 0) THEN

            allocate (wdsrf_bsn_ds  (numbasin))
            allocate (veloc_riv_ds  (numbasin))
            allocate (momen_riv_ds  (numbasin))
            allocate (hflux_fc      (numbasin))
            allocate (mflux_fc      (numbasin))
            allocate (zgrad_dn      (numbasin))
            allocate (sum_hflux_riv (numbasin))
            allocate (sum_mflux_riv (numbasin))
            allocate (sum_zgrad_riv (numbasin))
            allocate (bsnfilter     (numbasin))

            IF (DEF_Reservoir_Method > 0) THEN
               allocate (hflux_resv (numbasin))
               allocate (mflux_resv (numbasin))
            ENDIF

         ENDIF

         allocate (dt_res (numrivsys))
         allocate (dt_all (numrivsys))

         ntimestep_riverlake = 0

         dt_res(:) = dt

         DO WHILE (any(dt_res > 0))

            ntimestep_riverlake = ntimestep_riverlake + 1

            DO i = 1, numbasin
               bsnfilter(i) = dt_res(irivsys(i)) > 0
               IF (bsnfilter(i)) THEN
                  sum_hflux_riv(i) = 0.
                  sum_mflux_riv(i) = 0.
                  sum_zgrad_riv(i) = 0.
               ENDIF
            ENDDO

            WHERE (bsnfilter) ntacc_bsn = ntacc_bsn + 1

            CALL pull_from_downstream (wdsrf_bsn, wdsrf_bsn_ds, bsnfilter)
            CALL pull_from_downstream (veloc_riv, veloc_riv_ds, bsnfilter)
            CALL pull_from_downstream (momen_riv, momen_riv_ds, bsnfilter)

            ! velocity in ocean or inland depression is assumed to be 0.
            IF (numbasin > 0) THEN
               WHERE (riverdown <= 0)
                  veloc_riv_ds = 0.
               END WHERE
            ENDIF

            dt_all(:) = dt_res(:)

            DO i = 1, numbasin

               IF (.not. bsnfilter(i)) CYCLE

               IF (riverdown(i) >= 0) THEN

                  IF (riverdown(i) > 0) THEN
                     ! both elements are dry.
                     IF ((wdsrf_bsn(i) < RIVERMIN) .and. (wdsrf_bsn_ds(i) < RIVERMIN)) THEN
                        hflux_fc(i) = 0
                        mflux_fc(i) = 0
                        zgrad_dn(i) = 0
                        CYCLE
                     ENDIF
                  ENDIF

                  ! reconstruction of height of water near interface
                  IF (riverdown(i) > 0) THEN
                     bedelv_fc = max(bedelv(i), bedelv_ds(i))
                     IF ((lake_type(i) == 2) .or. (lake_type(i) == 3)) THEN
                        ! for reservoir (type=2) or controlled lake (type=3)
                        IF ((DEF_Reservoir_Method > 0) &
                           .and. (year >= dam_build_year(bsn2resv(i)))) THEN
                           bedelv_fc = max(bedelv_fc, dam_elv(bsn2resv(i)))
                        ENDIF
                     ENDIF
                     height_up = max(0., wdsrf_bsn(i)   +bedelv(i)   -bedelv_fc)
                     height_dn = max(0., wdsrf_bsn_ds(i)+bedelv_ds(i)-bedelv_fc)
                  ELSEIF (riverdown(i) == 0) THEN ! for river mouth
                     bedelv_fc = bedelv(i)
                     height_up = wdsrf_bsn(i)
                     ! sea level is assumed to be 0. and sea bed is assumed to be negative infinity.
                     height_dn = max(0., - bedelv_fc)
                  ENDIF

                  ! velocity at river downstream face (middle region in Riemann problem)
                  veloct_fc = 0.5 * (veloc_riv(i) + veloc_riv_ds(i)) &
                     + sqrt(grav * height_up) - sqrt(grav * height_dn)

                  ! height of water at downstream face (middle region in Riemann problem)
                  height_fc = 1/grav * (0.5*(sqrt(grav*height_up) + sqrt(grav*height_dn)) &
                     + 0.25 * (veloc_riv(i) - veloc_riv_ds(i))) ** 2

                  IF (height_up > 0) THEN
                     vwave_up = min(veloc_riv(i)-sqrt(grav*height_up), veloct_fc-sqrt(grav*height_fc))
                  ELSE
                     vwave_up = veloc_riv_ds(i) - 2.0 * sqrt(grav*height_dn)
                  ENDIF

                  IF (height_dn > 0) THEN
                     vwave_dn = max(veloc_riv_ds(i)+sqrt(grav*height_dn), veloct_fc+sqrt(grav*height_fc))
                  ELSE
                     vwave_dn = veloc_riv(i) + 2.0 * sqrt(grav*height_up)
                  ENDIF

                  hflux_up = veloc_riv(i)    * height_up
                  hflux_dn = veloc_riv_ds(i) * height_dn
                  mflux_up = veloc_riv(i)**2    * height_up + 0.5*grav * height_up**2
                  mflux_dn = veloc_riv_ds(i)**2 * height_dn + 0.5*grav * height_dn**2

                  IF (vwave_up >= 0.) THEN
                     hflux_fc(i) = outletwth(i) * hflux_up
                     mflux_fc(i) = outletwth(i) * mflux_up
                  ELSEIF (vwave_dn <= 0.) THEN
                     hflux_fc(i) = outletwth(i) * hflux_dn
                     mflux_fc(i) = outletwth(i) * mflux_dn
                  ELSE
                     hflux_fc(i) = outletwth(i) * (vwave_dn*hflux_up - vwave_up*hflux_dn &
                        + vwave_up*vwave_dn*(height_dn-height_up)) / (vwave_dn-vwave_up)
                     mflux_fc(i) = outletwth(i) * (vwave_dn*mflux_up - vwave_up*mflux_dn &
                        + vwave_up*vwave_dn*(hflux_dn-hflux_up)) / (vwave_dn-vwave_up)
                  ENDIF

                  sum_zgrad_riv(i) = sum_zgrad_riv(i) + outletwth(i) * 0.5*grav * height_up**2

                  zgrad_dn(i) = outletwth(i) * 0.5*grav * height_dn**2

               ELSEIF (riverdown(i) == -3) THEN
                  ! downstream is not in model region.
                  ! assume: 1. downstream river bed is equal to this river bed.
                  !         2. downstream water surface is equal to this river depth.
                  !         3. downstream water velocity is equal to this velocity.

                  veloc_riv(i) = max(veloc_riv(i), 0.)

                  IF (wdsrf_bsn(i) > riverdpth(i)) THEN

                     ! reconstruction of height of water near interface
                     height_up = wdsrf_bsn(i)
                     height_dn = riverdpth(i)

                     veloct_fc = veloc_riv(i) + sqrt(grav * height_up) - sqrt(grav * height_dn)
                     height_fc = 1/grav * (0.5*(sqrt(grav*height_up) + sqrt(grav*height_dn))) ** 2

                     vwave_up = min(veloc_riv(i)-sqrt(grav*height_up), veloct_fc-sqrt(grav*height_fc))
                     vwave_dn = max(veloc_riv(i)+sqrt(grav*height_dn), veloct_fc+sqrt(grav*height_fc))

                     hflux_up = veloc_riv(i) * height_up
                     hflux_dn = veloc_riv(i) * height_dn
                     mflux_up = veloc_riv(i)**2 * height_up + 0.5*grav * height_up**2
                     mflux_dn = veloc_riv(i)**2 * height_dn + 0.5*grav * height_dn**2

                     IF (vwave_up >= 0.) THEN
                        hflux_fc(i) = outletwth(i) * hflux_up
                        mflux_fc(i) = outletwth(i) * mflux_up
                     ELSEIF (vwave_dn <= 0.) THEN
                        hflux_fc(i) = outletwth(i) * hflux_dn
                        mflux_fc(i) = outletwth(i) * mflux_dn
                     ELSE
                        hflux_fc(i) = outletwth(i) * (vwave_dn*hflux_up - vwave_up*hflux_dn &
                           + vwave_up*vwave_dn*(height_dn-height_up)) / (vwave_dn-vwave_up)
                        mflux_fc(i) = outletwth(i) * (vwave_dn*mflux_up - vwave_up*mflux_dn &
                           + vwave_up*vwave_dn*(hflux_dn-hflux_up)) / (vwave_dn-vwave_up)
                     ENDIF

                     sum_zgrad_riv(i) = sum_zgrad_riv(i) + outletwth(i) * 0.5*grav * height_up**2

                  ELSE
                     hflux_fc(i) = 0
                     mflux_fc(i) = 0
                  ENDIF

               ELSEIF (riverdown(i) == -1) THEN ! inland depression
                  hflux_fc(i) = 0
                  mflux_fc(i) = 0
               ENDIF

               IF ((lake_id(i) < 0) .and. (hflux_fc(i) < 0)) THEN
                  dt_this = dt_all(irivsys(i))
                  hflux_fc(i) = &
                     max(hflux_fc(i), (height_up-height_dn) / dt_this * sum(hillslope_basin(i)%area, &
                     mask = hillslope_basin(i)%hand <= wdsrf_bsn(i) + handmin(i)))
               ENDIF

               ! reservoir operation.
               IF (DEF_Reservoir_Method > 0) THEN
                  IF (lake_type(i) == 2) THEN
                     hflux_fc(i) = 0.
                     mflux_fc(i) = 0.
                     zgrad_dn(i) = 0.
                  ENDIF
               ENDIF

               sum_hflux_riv(i) = sum_hflux_riv(i) + hflux_fc(i)
               sum_mflux_riv(i) = sum_mflux_riv(i) + mflux_fc(i)

            ENDDO

            IF (numbasin > 0) THEN
               hflux_fc = - hflux_fc;  mflux_fc = - mflux_fc;  zgrad_dn = - zgrad_dn
            ENDIF

            CALL push_to_downstream (hflux_fc, sum_hflux_riv, bsnfilter)
            CALL push_to_downstream (mflux_fc, sum_mflux_riv, bsnfilter)
            CALL push_to_downstream (zgrad_dn, sum_zgrad_riv, bsnfilter)

            IF (numbasin > 0) THEN
               hflux_fc = - hflux_fc;  mflux_fc = - mflux_fc;  zgrad_dn = - zgrad_dn
            ENDIF

            ! reservoir operation.
            IF (DEF_Reservoir_Method > 0) THEN

               DO i = 1, numbasin
                  IF (.not. bsnfilter(i)) CYCLE
                  IF (lake_type(i) == 2) THEN
                     IF (year >= dam_build_year(bsn2resv(i))) THEN

                        irsv = bsn2resv(i)
                        qresv_in(irsv) = - sum_hflux_riv(i)
                        volresv (irsv) = lakeinfo(i)%volume( wdsrf_bsn(i) )

                        IF (volresv(irsv) > 1.e-4 * volresv_total(irsv)) THEN
                           CALL reservoir_operation (DEF_Reservoir_Method, &
                              irsv, qresv_in(irsv), volresv(irsv), qresv_out(irsv))
                        ELSE
                           qresv_out (irsv) = 0.
                        ENDIF

                        hflux_fc(i) = qresv_out(irsv)
                        mflux_fc(i) = qresv_out(irsv) * sqrt(2*grav*wdsrf_bsn(i))

                        sum_hflux_riv(i) = sum_hflux_riv(i) + hflux_fc(i)
                        sum_mflux_riv(i) = sum_mflux_riv(i) + mflux_fc(i)

                        hflux_resv(i) = - hflux_fc(i)
                        mflux_resv(i) = - mflux_fc(i)

                     ELSE
                        hflux_resv(i) = 0.
                        mflux_resv(i) = 0.
                     ENDIF
                  ELSE
                     hflux_resv(i) = 0.
                     mflux_resv(i) = 0.
                  ENDIF
               ENDDO

               CALL push_to_downstream (hflux_resv, sum_hflux_riv, bsnfilter)
               CALL push_to_downstream (mflux_resv, sum_mflux_riv, bsnfilter)

            ENDIF

            DO i = 1, numbasin

               IF (.not. bsnfilter(i)) CYCLE

               dt_this = dt_all(irivsys(i))

               ! constraint 1: CFL condition (only for rivers)
               IF (lake_id(i) == 0) THEN
                  IF ((veloc_riv(i) /= 0.) .or. (wdsrf_bsn(i) > 0.)) THEN
                     dt_this = min(dt_this, riverlen(i)/(abs(veloc_riv(i))+sqrt(grav*wdsrf_bsn(i)))*0.8)
                  ENDIF
               ENDIF

               ! constraint 2: Avoid negative values of water
               IF (sum_hflux_riv(i) > 0) THEN
                  IF (lake_id(i) <= 0) THEN
                     ! for river or lake catchment
                     totalvolume = sum((wdsrf_bsn(i) + handmin(i) - hillslope_basin(i)%hand) &
                        * hillslope_basin(i)%area, &
                        mask = wdsrf_bsn(i) + handmin(i) >= hillslope_basin(i)%hand)
                  ELSEIF (lake_id(i) > 0) THEN
                     ! for lake
                     totalvolume = lakeinfo(i)%volume(wdsrf_bsn(i))
                  ENDIF

                  dt_this = min(dt_this, totalvolume / sum_hflux_riv(i))

               ENDIF

               ! constraint 3: Avoid change of flow direction (only for rivers)
               IF (lake_id(i) == 0) THEN
                  IF ((abs(veloc_riv(i)) > 0.1) &
                     .and. (veloc_riv(i) * (sum_mflux_riv(i)-sum_zgrad_riv(i)) > 0)) THEN
                     dt_this = min(dt_this, &
                        abs(momen_riv(i) * riverarea(i) / (sum_mflux_riv(i)-sum_zgrad_riv(i))))
                  ENDIF
               ENDIF

               dt_all(irivsys(i)) = min(dt_this, dt_all(irivsys(i)))

            ENDDO

#ifdef USEMPI
            IF (riversystem /= -1) THEN
               CALL mpi_allreduce (MPI_IN_PLACE, dt_all, 1, MPI_REAL8, MPI_MIN, p_comm_rivsys, p_err)
            ENDIF
#endif

            DO i = 1, numbasin

               IF (.not. bsnfilter(i)) CYCLE

               IF (lake_id(i) <= 0) THEN
                  ! rivers or lake catchments
                  hs = basin_hru%substt(i)
                  he = basin_hru%subend(i)
                  allocate (hmask (hillslope_basin(i)%nhru))

                  totalvolume = sum((wdsrf_bsn(i) + handmin(i) - hillslope_basin(i)%hand) &
                     * hillslope_basin(i)%area, &
                     mask = wdsrf_bsn(i) + handmin(i) >= hillslope_basin(i)%hand)

                  totalvolume = totalvolume - sum_hflux_riv(i) * dt_all(irivsys(i))

                  IF (totalvolume < VOLUMEMIN) THEN
                     DO j = 1, hillslope_basin(i)%nhru
                        IF (hillslope_basin(i)%hand(j) <= wdsrf_bsn(i) + handmin(i)) THEN
                           wdsrf_bsnhru(j+hs-1) = wdsrf_bsnhru(j+hs-1) &
                              - (wdsrf_bsn(i) + handmin(i) - hillslope_basin(i)%hand(j))
                        ENDIF
                     ENDDO
                     wdsrf_bsn(i) = 0
                  ELSE

                     dvol = sum_hflux_riv(i) * dt_all(irivsys(i))
                     IF (dvol > VOLUMEMIN) THEN
                        DO WHILE (dvol > VOLUMEMIN)
                           hmask = hillslope_basin(i)%hand < wdsrf_bsn(i) + handmin(i)
                           nextl = maxval(hillslope_basin(i)%hand, mask = hmask)
                           nexta = sum   (hillslope_basin(i)%area, mask = hmask)
                           nextv = nexta * (wdsrf_bsn(i)+handmin(i)-nextl)
                           IF (nextv > dvol) THEN
                              ddep = dvol/nexta
                              dvol = 0.
                           ELSE
                              ddep = wdsrf_bsn(i)+handmin(i) - nextl
                              dvol = dvol - nextv
                           ENDIF

                           wdsrf_bsn(i) = wdsrf_bsn(i) - ddep

                           DO j = 1, hillslope_basin(i)%nhru
                              IF (hmask(j)) THEN
                                 wdsrf_bsnhru(j+hs-1) = wdsrf_bsnhru(j+hs-1) - ddep
                              ENDIF
                           ENDDO
                        ENDDO
                     ELSEIF (dvol < -VOLUMEMIN) THEN
                        hmask = .true.
                        nexta = 0.
                        DO WHILE (dvol < -VOLUMEMIN)
                           IF (any(hmask)) THEN
                              j = minloc(hillslope_basin(i)%hand + wdsrf_bsnhru(hs:he), 1, mask = hmask)
                              nexta = nexta + hillslope_basin(i)%area(j)
                              hmask(j) = .false.
                           ENDIF
                           IF (any(hmask)) THEN
                              nextl = minval(hillslope_basin(i)%hand + wdsrf_bsnhru(hs:he), mask = hmask)
                              nextv = nexta*(nextl-(wdsrf_bsn(i)+handmin(i)))
                              IF ((-dvol) > nextv) THEN
                                 ddep = nextl - (wdsrf_bsn(i)+handmin(i))
                                 dvol = dvol + nextv
                              ELSE
                                 ddep = (-dvol)/nexta
                                 dvol = 0.
                              ENDIF
                           ELSE
                              ddep = (-dvol)/nexta
                              dvol = 0.
                           ENDIF

                           wdsrf_bsn(i) = wdsrf_bsn(i) + ddep

                           DO j = 1, hillslope_basin(i)%nhru
                              IF (.not. hmask(j)) THEN
                                 wdsrf_bsnhru(j+hs-1) = wdsrf_bsnhru(j+hs-1) + ddep
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDIF

                  ENDIF
                  deallocate(hmask)

               ELSE
                  totalvolume  = lakeinfo(i)%volume(wdsrf_bsn(i))
                  totalvolume  = totalvolume - sum_hflux_riv(i) * dt_all(irivsys(i))
                  wdsrf_bsn(i) = lakeinfo(i)%surface(totalvolume)
               ENDIF

               IF ((lake_id(i) /= 0) .or. (wdsrf_bsn(i) < RIVERMIN)) THEN
                  momen_riv(i) = 0
                  veloc_riv(i) = 0
               ELSE
                  friction = grav * nmanning_riv**2 / wdsrf_bsn(i)**(7.0/3.0) * abs(momen_riv(i))
                  momen_riv(i) = (momen_riv(i) &
                     - (sum_mflux_riv(i) - sum_zgrad_riv(i)) / riverarea(i) * dt_all(irivsys(i))) &
                     / (1 + friction * dt_all(irivsys(i)))
                  veloc_riv(i) = momen_riv(i) / wdsrf_bsn(i)
               ENDIF

               ! inland depression river
               IF ((lake_id(i) == 0) .and. (riverdown(i) == -1)) THEN
                  momen_riv(i) = min(0., momen_riv(i))
                  veloc_riv(i) = min(0., veloc_riv(i))
               ENDIF

               veloc_riv(i) = min(veloc_riv(i),  20.)
               veloc_riv(i) = max(veloc_riv(i), -20.)

            ENDDO

            DO i = 1, numbasin
               IF (bsnfilter(i)) THEN

                  wdsrf_bsn_ta (i) = wdsrf_bsn_ta (i) + wdsrf_bsn(i) * dt_all(irivsys(i))
                  momen_riv_ta (i) = momen_riv_ta (i) + momen_riv(i) * dt_all(irivsys(i))
                  discharge_ta (i) = discharge_ta (i) + hflux_fc (i) * dt_all(irivsys(i))

                  IF (lake_id(i) > 0) THEN ! for lakes
                     hs = basin_hru%substt(i)
                     he = basin_hru%subend(i)
                     DO j = hs, he
                        wdsrf_bsnhru(j) = &
                           max(wdsrf_bsn(i) - (lakeinfo(i)%depth(1) - lakeinfo(i)%depth0(j-hs+1)), 0.)
                        wdsrf_bsnhru_ta(j) = wdsrf_bsnhru_ta(j) + wdsrf_bsnhru(j) * dt_all(irivsys(i))
                     ENDDO
                  ENDIF

                  IF (DEF_Reservoir_Method > 0) THEN
                     IF (lake_type(i) == 2) THEN
                        IF (year >= dam_build_year(bsn2resv(i))) THEN
                           irsv = bsn2resv(i)
                           volresv_ta  (irsv) = volresv_ta  (irsv) + volresv  (irsv) * dt_all(irivsys(i))
                           qresv_in_ta (irsv) = qresv_in_ta (irsv) + qresv_in (irsv) * dt_all(irivsys(i))
                           qresv_out_ta(irsv) = qresv_out_ta(irsv) + qresv_out(irsv) * dt_all(irivsys(i))
                        ENDIF
                     ENDIF
                  ENDIF

               ENDIF
            ENDDO

            dt_res = dt_res - dt_all

         ENDDO

         IF (numbasin > 0) wdsrf_bsn_prev(:) = wdsrf_bsn(:)

         CALL mpi_allreduce (MPI_IN_PLACE, ntimestep_riverlake, 1, MPI_INTEGER, MPI_MAX, p_comm_worker, p_err)

         IF (allocated(wdsrf_bsn_ds )) deallocate(wdsrf_bsn_ds )
         IF (allocated(veloc_riv_ds )) deallocate(veloc_riv_ds )
         IF (allocated(momen_riv_ds )) deallocate(momen_riv_ds )
         IF (allocated(hflux_fc     )) deallocate(hflux_fc     )
         IF (allocated(mflux_fc     )) deallocate(mflux_fc     )
         IF (allocated(zgrad_dn     )) deallocate(zgrad_dn     )
         IF (allocated(hflux_resv   )) deallocate(hflux_resv   )
         IF (allocated(mflux_resv   )) deallocate(mflux_resv   )
         IF (allocated(sum_hflux_riv)) deallocate(sum_hflux_riv)
         IF (allocated(sum_mflux_riv)) deallocate(sum_mflux_riv)
         IF (allocated(sum_zgrad_riv)) deallocate(sum_zgrad_riv)
         IF (allocated(bsnfilter    )) deallocate(bsnfilter    )

      ENDIF

   END SUBROUTINE river_lake_flow

END MODULE MOD_Catch_RiverLakeFlow
#endif
