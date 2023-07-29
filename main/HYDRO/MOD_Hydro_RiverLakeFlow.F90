#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_Hydro_RiverLakeFlow
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
   
   REAL(r8), parameter :: RIVERMIN  = 1.e-4 
   REAL(r8), parameter :: nmanning_riv = 0.03
   
CONTAINS
   
   ! ---------
   SUBROUTINE river_lake_flow (dt)

      USE MOD_SPMD_Task
      USE MOD_Mesh
      USE MOD_LandHRU
      USE MOD_LandPatch
      USE MOD_Vars_TimeVariables
      USE MOD_Hydro_Vars_1DFluxes
      USE MOD_Hydro_RiverLakeNetwork
      USE MOD_Const_Physical, only : grav
      IMPLICIT NONE

      REAL(r8), intent(in) :: dt

      ! Local Variables
      INTEGER :: nriver
      INTEGER :: istt, iend, i, j
      
      REAL(r8), allocatable :: wdsrf_bsn_ds(:)
      REAL(r8), allocatable :: veloc_riv_ds (:)
      REAL(r8), allocatable :: momen_riv_ds (:)

      REAL(r8), allocatable :: hflux_fc(:)
      REAL(r8), allocatable :: mflux_fc(:)
      REAL(r8), allocatable :: zgrad_dn(:)
      
      REAL(r8), allocatable :: sum_hflux_riv(:)
      REAL(r8), allocatable :: sum_mflux_riv(:)
      REAL(r8), allocatable :: sum_zgrad_riv(:)
      
      REAL(r8) :: veloct_fc, height_fc, momen_fc, zsurf_fc
      REAL(r8) :: bedelv_fc, height_up, height_dn
      REAL(r8) :: vwave_up, vwave_dn, hflux_up, hflux_dn, mflux_up, mflux_dn
      REAL(r8) :: loss, totalvolume, friction
      REAL(r8) :: dt_res, dt_this
      CHARACTER(len=50) :: fmtt


      IF (p_is_worker) THEN
      
         nriver = numelm
         nbasin = numelm
         
         ! update water depth in basin by aggregating water depths in patches
         DO i = 1, nbasin
            IF (lake_id(i) <= 0) THEN
               ! river or lake catchment
               istt = basin_hru%substt(i)
               iend = basin_hru%subend(i)

               wdsrf_this = minval(surface_network%hand + wdsrf_hru(istt:iend))

               totalvolume = sum(wdsrf_this - surface_network%hand, mask = surface_network%hand <= wdsrf_this)
               IF (totalvolume < 1.0e-4) THEN
                  wdsrf_this = 0
               ENDIF
            ELSEIF (lake_id(i) > 0) THEN
               ! lake : the whole lake is an element 
               istt = elm_patch%substt(i)
               iend = elm_patch%subend(i)
               totalvolume = sum(wdsrf(istt:iend)/1.0e3 * lakes(i)%area0) ! mm to m
               wdsrf_this  = lakes(i)%surface(totalvolume)
            ENDIF

            ! river momentum is less or equal than the momentum at last time step.
            IF (lake_id(i) == 0) THEN
               momen_riv(i) = min(wdsrf_bsn(i), wdsrf_this) * veloc_riv(i)
            ELSE
               momen_riv(i) = 0
            ENDIF

            wdsrf_bsn(i) = wdsrf_this

         ENDDO

         IF (nriver > 0) THEN
            allocate (wdsrf_bsn_ds  (nriver))
            allocate (veloc_riv_ds  (nriver))
            allocate (momen_riv_ds  (nriver))
            allocate (hflux_fc      (nriver))
            allocate (mflux_fc      (nriver))
            allocate (zgrad_dn      (nriver))
            allocate (sum_hflux_riv (nriver))
            allocate (sum_mflux_riv (nriver))
            allocate (sum_zgrad_riv (nriver))
         ENDIF

         dt_res = dt
         DO WHILE (dt_res > 0)
            
            DO i = 1, nriver
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
               
            dt_this = dt_res

            DO i = 1, nriver
               IF (riverdown(i) >= 0) THEN
                  
                  ! velocity in lakes is assumed to be 0.
                  IF (lake_id(i) > 0) THEN
                     veloc_riv(i) = 0.
                  ENDIF
                  ! velocity in lakes or ocean is assumed to be 0.
                  IF (to_lake(i) .or. (riverdown(i) == 0)) THEN
                     veloc_riv_ds(i) = 0.
                  ENDIF

                  IF (riverdown(i) > 0) THEN
                     ! both elements are dry.
                     IF ((wdsrf_bsn(i) < RIVERMIN) .and. (wdsrf_bsn_ds(i) < RIVERMIN)) THEN
                        hflux_fc(i) = 0
                        mflux_fc(i) = 0
                        zgrad_dn(i) = 0
                        cycle
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
                  
               ELSE ! inland depression
                  hflux_fc(i) = 0
                  mflux_fc(i) = 0
                  zgrad_dn(i) = 0
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
            hflux_fc = - hflux_fc
            mflux_fc = - mflux_fc
            zgrad_dn = - zgrad_dn
            CALL river_data_exchange (SEND_DATA_UP_TO_DOWN, accum = .true., &
               vec_send1 = hflux_fc, vec_recv1 = sum_hflux_riv, &
               vec_send2 = mflux_fc, vec_recv2 = sum_mflux_riv, &
               vec_send3 = zgrad_dn, vec_recv3 = sum_zgrad_riv)
#endif

            DO i = 1, nriver
               ! constraint 1: CFL condition
               IF (lake_id(i) <= 0) THEN
                  IF ((veloc_riv(i) /= 0.) .or. (wdsrf_bsn(i) > 0.)) THEN
                     dt_this = min(dt_this, riverlen(i)/(abs(veloc_riv(i))+sqrt(grav*wdsrf_bsn(i)))*0.8)
                  ENDIF
               ENDIF

               ! constraint 2: Avoid negative values of water
               IF (sum_hflux_riv(i) > 0) THEN
                  IF (lake_id(i) <= 0) THEN
                     ! for river or lake catchment
                     istt = basin_hru%substt(i)
                     iend = basin_hru%subend(i)
                     allocate (mask (surface_network(i)%nhru))

                     mask = surface_network(i)%hand <= wdsrf_bsn(i)
                     totalvolume = sum(wdsrf_bsn(i) - surface_network(i)%hand, mask = mask)

                     deallocate(mask)
                  ELSEIF (lake_id(i) > 0) THEN
                     ! for lake
                     totalvolume = lakes(i)%volume(wdsrf_bsn(i))
                  ENDIF
               
                  dt_this = min(dt_this, totalvolume / sum_hflux_riv(i))
                  
               ENDIF

               ! constraint 3: Avoid change of flow direction
               IF (lake_id(i) <= 0) THEN
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

            DO i = 1, nriver

               IF (lake_id(i) <= 0) THEN
                  istt = basin_hru%substt(i)
                  iend = basin_hru%subend(i)
                  
                  allocate (mask (surface_network(i)%nhru))
                  dvol = sum_hflux_riv(i) * dt_this
                  IF (dvol > 0) THEN
                     DO WHILE (dvol > 0)
                        mask  = surface_network(i)%hand < wdsrf_bsn(i)
                        nextl = maxval(surface_network(i)%hand, mask = mask)
                        nexta = sum   (surface_network(i)%area, mask = mask) 
                        nextv = nexta * (wdsrf_bsn(i)-nextl)
                        IF (nextv > dvol) THEN
                           ddep = dvol/nexta
                           dvol = 0
                        ELSE
                           ddep = wdsrf_bsn(i) - nextl
                           dvol = dvol - (wdsrf_bsn(i)-nextl)*nexta
                        ENDIF
                           
                        wdsrf_bsn(i) = wdsrf_bsn(i) - ddep

                        DO j = 1, surface_network(i)%nhru
                           IF (mask(j)) THEN
                              wdsrf_hru(j+istt-1) = wdsrf_hru(j+istt-1) - ddep
                           ENDIF
                        ENDDO
                     ENDDO
                  ELSEIF (dvol < 0) THEN
                     DO WHILE (dvol < 0)
                        mask  = surface_network(i)%hand + wdsrf_hru(istt:iend) > wdsrf_bsn(i)
                        nexta = sum(surface_network(i)%area, mask = (.not. mask)) 
                        IF (any(mask)) THEN
                           nextl = minval(surface_network(i)%hand + wdsrf_hru(istt:iend), mask = mask)
                           nextv = nexta*(nextl-wdsrf_bsn(i))
                           IF ((-dvol) > nextv) THEN
                              ddep = nextl - wdsrf_bsn(i)
                              dvol = dvol + nextv
                           ELSE
                              ddep = (-dvol)/nexta
                              dvol = 0
                           ENDIF
                        ELSE
                           ddep = (-dvol)/nexta
                           dvol = 0
                        ENDIF

                        wdsrf_bsn(i) = wdsrf_bsn(i) + ddep

                        DO j = 1, surface_network(i)%nhru
                           IF (.not. mask(j)) THEN
                              wdsrf_hru(j+istt-1) = wdsrf_hru(j+istt-1) + ddep
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
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
            
               ! inland depression
               IF (riverdown(i) < 0) THEN
                  momen_riv(i) = min(0., momen_riv(i))
                  veloc_riv(i) = min(0., veloc_riv(i))
               ENDIF

            ENDDO

            IF (nriver > 0) THEN
               wdsrf_bsn_ta(:) = wdsrf_bsn_ta(:) + wdsrf_bsn(:) * dt_this
               momen_riv_ta(:) = momen_riv_ta(:) + momen_riv(:) * dt_this
            ENDIF

            dt_res = dt_res - dt_this

         ENDDO

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
   
END MODULE MOD_Hydro_RiverLakeFlow
#endif
