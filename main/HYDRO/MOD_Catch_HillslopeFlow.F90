#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_HillslopeFlow
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!   
!   Shallow water equation solver over hillslopes.
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

   ! -- Parameters --
   real(r8), parameter :: PONDMIN  = 1.e-4 
   real(r8), parameter :: nmanning_hslp = 0.3
   
CONTAINS
   
   ! ----------
   SUBROUTINE hillslope_flow (dt)

   USE MOD_SPMD_Task
   USE MOD_Catch_BasinNetwork
   USE MOD_Catch_HillslopeNetwork
   USE MOD_Catch_RiverLakeNetwork
   USE MOD_Catch_Vars_TimeVariables
   USE MOD_Catch_Vars_1DFluxes
   USE MOD_Const_Physical, only: grav

   IMPLICIT NONE
   
   real(r8), intent(in) :: dt

   ! Local Variables
   integer :: nhru, hs, he, ibasin, i, j

   type(hillslope_network_type), pointer :: hillslope

   real(r8), allocatable :: wdsrf_h (:) ! [m]
   real(r8), allocatable :: momen_h (:) ! [m^2/s]
   real(r8), allocatable :: veloc_h (:) ! [m/s]

   real(r8), allocatable :: sum_hflux_h (:) 
   real(r8), allocatable :: sum_mflux_h (:) 
   real(r8), allocatable :: sum_zgrad_h (:) 
   
   real(r8) :: hand_fc,  wdsrf_fc, veloc_fc, hflux_fc, mflux_fc
   real(r8) :: wdsrf_up, wdsrf_dn, vwave_up, vwave_dn
   real(r8) :: hflux_up, hflux_dn, mflux_up, mflux_dn
   
   real(r8), allocatable :: xsurf_h (:) ! [m/s]

   real(r8) :: friction
   real(r8) :: dt_res, dt_this

   logical, allocatable :: mask(:)
   real(r8) :: srfbsn, dvol, nextl, nexta, nextv, ddep

      IF (p_is_worker) THEN

         DO ibasin = 1, numbasin

            hs = basin_hru%substt(ibasin)
            he = basin_hru%subend(ibasin)

            IF (lake_id(ibasin) > 0) THEN
               veloc_bsnhru(hs:he) = 0
               momen_bsnhru(hs:he) = 0
               CYCLE ! skip lakes
            ELSE
               DO i = hs, he
                  ! momentum is less or equal than the momentum at last time step.
                  momen_bsnhru(i) = min(wdsrf_bsnhru_prev(i), wdsrf_bsnhru(i)) * veloc_bsnhru(i)
               ENDDO
            ENDIF

            hillslope => hillslope_basin(ibasin)

            nhru = hillslope%nhru

            allocate (wdsrf_h (nhru))
            allocate (veloc_h (nhru))
            allocate (momen_h (nhru))

            allocate (sum_hflux_h (nhru))
            allocate (sum_mflux_h (nhru))
            allocate (sum_zgrad_h (nhru))
               
            allocate (xsurf_h (nhru))

            DO i = 1, nhru
               wdsrf_h(i) = wdsrf_bsnhru(hillslope%ihru(i)) 
               momen_h(i) = momen_bsnhru(hillslope%ihru(i)) 
               IF (wdsrf_h(i) > 0.) THEN
                  veloc_h(i) = momen_h(i) / wdsrf_h(i)
               ELSE
                  veloc_h(i) = 0.
               ENDIF
            ENDDO

            dt_res = dt 
            DO WHILE (dt_res > 0.)

               DO i = 1, nhru
                  sum_hflux_h(i) = 0.
                  sum_mflux_h(i) = 0.
                  sum_zgrad_h(i) = 0.
               ENDDO 

               dt_this = dt_res

               DO i = 1, nhru

                  j = hillslope%inext(i)

                  IF (j <= 0) CYCLE ! lowest HRUs

                  ! dry HRU
                  IF ((wdsrf_h(i) < PONDMIN) .and. (wdsrf_h(j) < PONDMIN)) THEN
                     CYCLE
                  ENDIF

                  ! reconstruction of height of water near interface
                  hand_fc  = max(hillslope%hand(i), hillslope%hand(j))
                  wdsrf_up = max(0., hillslope%hand(i)+wdsrf_h(i) - hand_fc)
                  wdsrf_dn = max(0., hillslope%hand(j)+wdsrf_h(j) - hand_fc)

                  ! velocity at hydrounit downstream face
                  veloc_fc = 0.5 * (veloc_h(i) + veloc_h(j)) &
                     + sqrt(grav * wdsrf_up) - sqrt(grav * wdsrf_dn)

                  ! depth of water at downstream face
                  wdsrf_fc = 1/grav * (0.5*(sqrt(grav*wdsrf_up) + sqrt(grav*wdsrf_dn)) &
                     + 0.25 * (veloc_h(i) - veloc_h(j)))**2.0

                  IF (wdsrf_up > 0) THEN
                     vwave_up = min(veloc_h(i)-sqrt(grav*wdsrf_up), veloc_fc-sqrt(grav*wdsrf_fc))
                  ELSE
                     vwave_up = veloc_h(j) - 2.0 * sqrt(grav*wdsrf_dn)
                  ENDIF

                  IF (wdsrf_dn > 0) THEN
                     vwave_dn = max(veloc_h(j)+sqrt(grav*wdsrf_dn), veloc_fc+sqrt(grav*wdsrf_fc))
                  ELSE
                     vwave_dn = veloc_h(i) + 2.0 * sqrt(grav*wdsrf_up)
                  ENDIF

                  hflux_up = veloc_h(i) * wdsrf_up
                  hflux_dn = veloc_h(j) * wdsrf_dn
                  mflux_up = veloc_h(i)**2 * wdsrf_up + 0.5*grav * wdsrf_up**2
                  mflux_dn = veloc_h(j)**2 * wdsrf_dn + 0.5*grav * wdsrf_dn**2

                  IF (vwave_up >= 0.) THEN
                     hflux_fc = hillslope%flen(i) * hflux_up
                     mflux_fc = hillslope%flen(i) * mflux_up
                  ELSEIF (vwave_dn <= 0.) THEN
                     hflux_fc = hillslope%flen(i) * hflux_dn
                     mflux_fc = hillslope%flen(i) * mflux_dn
                  ELSE
                     hflux_fc = hillslope%flen(i) * (vwave_dn*hflux_up - vwave_up*hflux_dn &
                        + vwave_up*vwave_dn*(wdsrf_dn-wdsrf_up)) / (vwave_dn-vwave_up)
                     mflux_fc = hillslope%flen(i) * (vwave_dn*mflux_up - vwave_up*mflux_dn &
                        + vwave_up*vwave_dn*(hflux_dn-hflux_up)) / (vwave_dn-vwave_up)
                  ENDIF

                  sum_hflux_h(i) = sum_hflux_h(i) + hflux_fc
                  sum_hflux_h(j) = sum_hflux_h(j) - hflux_fc

                  sum_mflux_h(i) = sum_mflux_h(i) + mflux_fc
                  sum_mflux_h(j) = sum_mflux_h(j) - mflux_fc
                  
                  sum_zgrad_h(i) = sum_zgrad_h(i) + hillslope%flen(i) * 0.5*grav * wdsrf_up**2
                  sum_zgrad_h(j) = sum_zgrad_h(j) - hillslope%flen(i) * 0.5*grav * wdsrf_dn**2 

               ENDDO

               DO i = 1, nhru
                  ! constraint 1: CFL condition
                  IF (hillslope%inext(i) > 0) THEN
                     IF ((veloc_h(i) /= 0.) .or. (wdsrf_h(i) > 0.)) THEN
                        dt_this = min(dt_this, hillslope%plen(i)/(abs(veloc_h(i)) + sqrt(grav*wdsrf_h(i)))*0.8)
                     ENDIF
                  ENDIF

                  ! constraint 2: Avoid negative values of water
                  xsurf_h(i) = sum_hflux_h(i) / hillslope%area(i)
                  IF (xsurf_h(i) > 0) THEN
                     dt_this = min(dt_this, wdsrf_h(i) / xsurf_h(i))
                  ENDIF
                     
                  ! constraint 3: Avoid change of flow direction
                  IF ((abs(veloc_h(i)) > 0.1) &
                     .and. (veloc_h(i) * (sum_mflux_h(i) - sum_zgrad_h(i)) > 0)) THEN
                     dt_this = min(dt_this, &
                        abs(momen_h(i) * hillslope%area(i) / (sum_mflux_h(i) - sum_zgrad_h(i))))
                  ENDIF
               ENDDO 

               DO i = 1, nhru

                  wdsrf_h(i) = max(0., wdsrf_h(i) - xsurf_h(i) * dt_this)

                  IF (wdsrf_h(i) < PONDMIN) THEN
                     momen_h(i) = 0
                  ELSE
                     friction = grav * nmanning_hslp**2 * abs(momen_h(i)) / wdsrf_h(i)**(7.0/3.0) 
                     momen_h(i) = (momen_h(i) - &
                        (sum_mflux_h(i) - sum_zgrad_h(i)) / hillslope%area(i) * dt_this) &
                        / (1 + friction * dt_this) 

                     IF (hillslope%inext(i) <= 0) THEN
                        momen_h(i) = min(momen_h(i), 0.)
                     ENDIF

                     IF (all(hillslope%inext /= i)) THEN
                        momen_h(i) = max(momen_h(i), 0.)
                     ENDIF
                  ENDIF

               ENDDO

               IF (hillslope%indx(1) == 0) THEN
                  srfbsn = minval(hillslope%hand + wdsrf_h)
                  IF (srfbsn < wdsrf_h(1)) THEN
                     allocate (mask (hillslope%nhru))
                     dvol = (wdsrf_h(1) - srfbsn) * hillslope%area(1)
                     momen_h(1) = srfbsn/wdsrf_h(1) * momen_h(1)
                     wdsrf_h(1) = srfbsn
                     DO WHILE (dvol > 0)
                        mask  = hillslope%hand + wdsrf_h > srfbsn
                        nexta = sum(hillslope%area, mask = (.not. mask)) 
                        IF (any(mask)) THEN
                           nextl = minval(hillslope%hand + wdsrf_h, mask = mask)
                           nextv = nexta*(nextl-srfbsn)
                           IF (dvol > nextv) THEN
                              ddep = nextl - srfbsn
                              dvol = dvol - nextv
                           ELSE
                              ddep = dvol/nexta
                              dvol = 0.
                           ENDIF
                        ELSE
                           ddep = dvol/nexta
                           dvol = 0.
                        ENDIF

                        srfbsn = srfbsn + ddep

                        WHERE (.not. mask)
                           wdsrf_h = wdsrf_h + ddep
                        END WHERE 
                     ENDDO
                     deallocate(mask)
                  ENDIF
               ENDIF

               DO i = 1, nhru
                  IF (wdsrf_h(i) < PONDMIN) THEN
                     veloc_h(i) = 0
                  ELSE
                     veloc_h(i) = momen_h(i) / wdsrf_h(i)
                  ENDIF

                  wdsrf_bsnhru_ta(hillslope%ihru(i)) = wdsrf_bsnhru_ta(hillslope%ihru(i)) + wdsrf_h(i) * dt_this
                  momen_bsnhru_ta(hillslope%ihru(i)) = momen_bsnhru_ta(hillslope%ihru(i)) + momen_h(i) * dt_this
               ENDDO

               dt_res = dt_res - dt_this

            ENDDO
            
            ! SAVE depth of surface water
            DO i = 1, nhru
               wdsrf_bsnhru(hillslope%ihru(i)) = wdsrf_h(i)
               veloc_bsnhru(hillslope%ihru(i)) = veloc_h(i)
            ENDDO

            deallocate (wdsrf_h)
            deallocate (veloc_h)
            deallocate (momen_h)

            deallocate (sum_hflux_h)
            deallocate (sum_mflux_h)
            deallocate (sum_zgrad_h)

            deallocate (xsurf_h)

         ENDDO

         IF (numbsnhru > 0) wdsrf_bsnhru_prev(:) = wdsrf_bsnhru(:)

      ENDIF

   END SUBROUTINE hillslope_flow

END MODULE MOD_Catch_HillslopeFlow
#endif
