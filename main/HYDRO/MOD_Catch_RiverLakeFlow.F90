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
!      Computers and Fluids, 38(2), 221–234.
!   [3] Audusse, E., Bouchut, F., Bristeau, M.-O., Klein, R., Perthame, B. (2004). 
!      A Fast and Stable Well-Balanced Scheme with Hydrostatic Reconstruction for 
!      Shallow Water Flows. SIAM Journal on Scientific Computing, 25(6), 2050–2065.
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
   SUBROUTINE river_lake_flow (dt)

   USE MOD_SPMD_Task
   USE MOD_Mesh
   USE MOD_LandHRU
   USE MOD_LandPatch
   USE MOD_Vars_TimeVariables
   USE MOD_Hydro_Vars_1DFluxes
   USE MOD_Catch_HillslopeNetwork
   USE MOD_Catch_RiverLakeNetwork
   USE MOD_Const_Physical, only : grav
   IMPLICIT NONE

   real(r8), intent(in) :: dt

   ! Local Variables
   integer :: nbasin
   integer :: hs, he, i, j
   
   real(r8), allocatable :: wdsrf_bsn_ds(:)
   real(r8), allocatable :: veloc_riv_ds(:)
   real(r8), allocatable :: momen_riv_ds(:)

   real(r8), allocatable :: hflux_fc(:)
   real(r8), allocatable :: mflux_fc(:)
   real(r8), allocatable :: zgrad_dn(:)
   
   real(r8), allocatable :: sum_hflux_riv(:)
   real(r8), allocatable :: sum_mflux_riv(:)
   real(r8), allocatable :: sum_zgrad_riv(:)
   
   real(r8) :: veloct_fc, height_fc, momen_fc, zsurf_fc
   real(r8) :: bedelv_fc, height_up, height_dn
   real(r8) :: vwave_up, vwave_dn, hflux_up, hflux_dn, mflux_up, mflux_dn
   real(r8) :: totalvolume, loss, friction, dvol, nextl, nexta, nextv, ddep
   real(r8) :: dt_res, dt_this
   logical, allocatable :: mask(:)

      
      IF (p_is_worker) THEN
      
         nbasin = numelm

         ! update water depth in basin by aggregating water depths in patches
         DO i = 1, nbasin
            hs = basin_hru%substt(i)
            he = basin_hru%subend(i)

            IF (lake_id(i) <= 0) THEN
               ! river or lake catchment
               ! Water surface in a basin is defined as the lowest surface water in the basin
               wdsrf_bsn(i) = minval(hillslope_network(i)%hand + wdsrf_hru(hs:he)) - handmin(i)

            ELSEIF (lake_id(i) > 0) THEN
               ! lake 
               totalvolume  = sum(wdsrf_hru(hs:he) * lakes(i)%area0)
               wdsrf_bsn(i) = lakes(i)%surface(totalvolume)
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

         IF (nbasin > 0) THEN
            allocate (wdsrf_bsn_ds  (nbasin))
            allocate (veloc_riv_ds  (nbasin))
            allocate (momen_riv_ds  (nbasin))
            allocate (hflux_fc      (nbasin))
            allocate (mflux_fc      (nbasin))
            allocate (zgrad_dn      (nbasin))
            allocate (sum_hflux_riv (nbasin))
            allocate (sum_mflux_riv (nbasin))
            allocate (sum_zgrad_riv (nbasin))
         ENDIF

         ntimestep_riverlake = 0
         dt_res = dt
         DO WHILE (dt_res > 0)

            ntimestep_riverlake = ntimestep_riverlake + 1
            
            DO i = 1, nbasin
               sum_hflux_riv(i) = 0.
               sum_mflux_riv(i) = 0.
               sum_zgrad_riv(i) = 0.
               
               IF (addrdown(i) > 0) THEN
                  wdsrf_bsn_ds(i) = wdsrf_bsn(addrdown(i))
                  veloc_riv_ds(i) = veloc_riv(addrdown(i))
                  momen_riv_ds(i) = momen_riv(addrdown(i))
               ELSE
                  wdsrf_bsn_ds(i) = 0
                  veloc_riv_ds(i) = 0
                  momen_riv_ds(i) = 0
               ENDIF
            ENDDO 
#ifdef USEMPI
            CALL river_data_exchange (SEND_DATA_DOWN_TO_UP, accum = .false., &
               vec_send1 = wdsrf_bsn, vec_recv1 = wdsrf_bsn_ds, &
               vec_send2 = veloc_riv, vec_recv2 = veloc_riv_ds, &
               vec_send3 = momen_riv, vec_recv3 = momen_riv_ds )
#endif
            ! velocity in ocean or inland depression is assumed to be 0.
            WHERE (riverdown <= 0)
               veloc_riv_ds = 0.
            END WHERE

            dt_this = dt_res

            DO i = 1, nbasin
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
                  hflux_fc(i) = &
                     max(hflux_fc(i), (height_up-height_dn) / dt_this * sum(hillslope_network(i)%area, &
                     mask = hillslope_network(i)%hand <= wdsrf_bsn(i) + handmin(i)))
               ENDIF

               sum_hflux_riv(i) = sum_hflux_riv(i) + hflux_fc(i)
               sum_mflux_riv(i) = sum_mflux_riv(i) + mflux_fc(i)

               IF (addrdown(i) > 0) THEN
                  j = addrdown(i)
                  sum_hflux_riv(j) = sum_hflux_riv(j) - hflux_fc(i)
                  sum_mflux_riv(j) = sum_mflux_riv(j) - mflux_fc(i)
                  sum_zgrad_riv(j) = sum_zgrad_riv(j) - zgrad_dn(i)
               ENDIF
            
            ENDDO
            
#ifdef USEMPI
            hflux_fc = - hflux_fc;  mflux_fc = - mflux_fc;  zgrad_dn = - zgrad_dn

            CALL river_data_exchange (SEND_DATA_UP_TO_DOWN, accum = .true., &
               vec_send1 = hflux_fc, vec_recv1 = sum_hflux_riv, &
               vec_send2 = mflux_fc, vec_recv2 = sum_mflux_riv, &
               vec_send3 = zgrad_dn, vec_recv3 = sum_zgrad_riv)
            
            hflux_fc = - hflux_fc;  mflux_fc = - mflux_fc;  zgrad_dn = - zgrad_dn
#endif

            DO i = 1, nbasin
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
                     totalvolume = sum((wdsrf_bsn(i) + handmin(i) - hillslope_network(i)%hand) &
                        * hillslope_network(i)%area, &
                        mask = wdsrf_bsn(i) + handmin(i) >= hillslope_network(i)%hand) 
                  ELSEIF (lake_id(i) > 0) THEN
                     ! for lake
                     totalvolume = lakes(i)%volume(wdsrf_bsn(i))
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
            ENDDO 

#ifdef USEMPI
            CALL mpi_allreduce (MPI_IN_PLACE, dt_this, 1, MPI_REAL8, MPI_MIN, p_comm_worker, p_err)
#endif

            DO i = 1, nbasin

               IF (lake_id(i) <= 0) THEN
                  ! rivers or lake catchments
                  hs = basin_hru%substt(i)
                  he = basin_hru%subend(i)
                  allocate (mask (hillslope_network(i)%nhru))
                  
                  totalvolume = sum((wdsrf_bsn(i) + handmin(i) - hillslope_network(i)%hand) &
                     * hillslope_network(i)%area, &
                     mask = wdsrf_bsn(i) + handmin(i) >= hillslope_network(i)%hand) 

                  totalvolume = totalvolume - sum_hflux_riv(i) * dt_this
                           
                  IF (totalvolume < VOLUMEMIN) THEN
                     DO j = 1, hillslope_network(i)%nhru
                        IF (hillslope_network(i)%hand(j) <= wdsrf_bsn(i) + handmin(i)) THEN
                           wdsrf_hru(j+hs-1) = wdsrf_hru(j+hs-1) &
                              - (wdsrf_bsn(i) + handmin(i) - hillslope_network(i)%hand(j))
                        ENDIF
                     ENDDO
                     wdsrf_bsn(i) = 0
                  ELSE

                     dvol = sum_hflux_riv(i) * dt_this
                     IF (dvol > VOLUMEMIN) THEN
                        DO WHILE (dvol > VOLUMEMIN)
                           mask  = hillslope_network(i)%hand < wdsrf_bsn(i) + handmin(i)
                           nextl = maxval(hillslope_network(i)%hand, mask = mask)
                           nexta = sum   (hillslope_network(i)%area, mask = mask) 
                           nextv = nexta * (wdsrf_bsn(i)+handmin(i)-nextl)
                           IF (nextv > dvol) THEN
                              ddep = dvol/nexta
                              dvol = 0.
                           ELSE
                              ddep = wdsrf_bsn(i)+handmin(i) - nextl
                              dvol = dvol - nextv
                           ENDIF

                           wdsrf_bsn(i) = wdsrf_bsn(i) - ddep

                           DO j = 1, hillslope_network(i)%nhru
                              IF (mask(j)) THEN
                                 wdsrf_hru(j+hs-1) = wdsrf_hru(j+hs-1) - ddep
                              ENDIF
                           ENDDO
                        ENDDO
                     ELSEIF (dvol < -VOLUMEMIN) THEN
                        mask  = .true.
                        nexta = 0.
                        DO WHILE (dvol < -VOLUMEMIN)
                           IF (any(mask)) THEN
                              j = minloc(hillslope_network(i)%hand + wdsrf_hru(hs:he), 1, mask = mask)
                              nexta = nexta + hillslope_network(i)%area(j)
                              mask(j) = .false.
                           ENDIF
                           IF (any(mask)) THEN
                              nextl = minval(hillslope_network(i)%hand + wdsrf_hru(hs:he), mask = mask)
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

                           DO j = 1, hillslope_network(i)%nhru
                              IF (.not. mask(j)) THEN
                                 wdsrf_hru(j+hs-1) = wdsrf_hru(j+hs-1) + ddep
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDIF
                     
                  ENDIF
                  deallocate(mask)

               ELSE
                  totalvolume  = lakes(i)%volume(wdsrf_bsn(i))
                  totalvolume  = totalvolume - sum_hflux_riv(i) * dt_this
                  wdsrf_bsn(i) = lakes(i)%surface(totalvolume) 
               ENDIF

               IF ((lake_id(i) /= 0) .or. (wdsrf_bsn(i) < RIVERMIN)) THEN
                  momen_riv(i) = 0
                  veloc_riv(i) = 0
               ELSE
                  friction = grav * nmanning_riv**2 / wdsrf_bsn(i)**(7.0/3.0) * abs(momen_riv(i))
                  momen_riv(i) = (momen_riv(i) &
                     - (sum_mflux_riv(i) - sum_zgrad_riv(i)) / riverarea(i) * dt_this) &
                     / (1 + friction * dt_this) 
                  veloc_riv(i) = momen_riv(i) / wdsrf_bsn(i)
               ENDIF
            
               ! inland depression river
               IF ((lake_id(i) == 0) .and. (riverdown(i) == -1)) THEN
                  momen_riv(i) = min(0., momen_riv(i))
                  veloc_riv(i) = min(0., veloc_riv(i))
               ENDIF

            ENDDO

            IF (nbasin > 0) THEN
               wdsrf_bsn_ta(:) = wdsrf_bsn_ta(:) + wdsrf_bsn(:) * dt_this
               momen_riv_ta(:) = momen_riv_ta(:) + momen_riv(:) * dt_this
               discharge   (:) = discharge   (:) + hflux_fc (:) * dt_this
            ENDIF
         
            DO i = 1, nbasin
               IF (lake_id(i) > 0) THEN ! for lakes
                  hs = basin_hru%substt(i)
                  he = basin_hru%subend(i)
                  DO j = hs, he
                     wdsrf_hru(j) = max(wdsrf_bsn(i) - (lakes(i)%depth(1) - lakes(i)%depth0(j-hs+1)), 0.)
                     wdsrf_hru_ta(j) = wdsrf_hru_ta(j) + wdsrf_hru(j) * dt_this
                  ENDDO
               ENDIF
            ENDDO

            dt_res = dt_res - dt_this

         ENDDO

         wdsrf_bsn_prev(:) = wdsrf_bsn(:)

         IF (allocated(wdsrf_bsn_ds )) deallocate(wdsrf_bsn_ds )
         IF (allocated(veloc_riv_ds )) deallocate(veloc_riv_ds )
         IF (allocated(momen_riv_ds )) deallocate(momen_riv_ds )
         IF (allocated(hflux_fc     )) deallocate(hflux_fc     )
         IF (allocated(mflux_fc     )) deallocate(mflux_fc     )
         IF (allocated(zgrad_dn     )) deallocate(zgrad_dn     )
         IF (allocated(sum_hflux_riv)) deallocate(sum_hflux_riv)
         IF (allocated(sum_mflux_riv)) deallocate(sum_mflux_riv)
         IF (allocated(sum_zgrad_riv)) deallocate(sum_zgrad_riv)

      ENDIF

   END SUBROUTINE river_lake_flow

END MODULE MOD_Catch_RiverLakeFlow
#endif
