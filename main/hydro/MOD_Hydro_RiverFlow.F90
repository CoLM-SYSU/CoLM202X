#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_Hydro_RiverFlow

   USE MOD_Precision
   USE MOD_Hydro_RiverNetwork
   IMPLICIT NONE
   
   REAL(r8), parameter :: RIVERMIN  = 1.e-4 
   REAL(r8), parameter :: nmanning_riv = 0.03
   
CONTAINS
   
   ! ---------
   SUBROUTINE river_flow (dt)

      USE MOD_SPMD_Task
      USE MOD_Mesh
      USE mod_landhru
      USE mod_landpatch
      USE MOD_Vars_TimeVariables
      USE MOD_Hydro_Vars_1DFluxes
      USE PhysicalConstants, only : grav
      USE MOD_CoLMDebug
      IMPLICIT NONE

      REAL(r8), intent(in) :: dt

      ! Local Variables
      INTEGER :: nriver
      INTEGER :: istt, iend, i, j
      REAL(r8), allocatable :: height_riv(:)
      REAL(r8), allocatable :: veloct_riv(:)
      REAL(r8), allocatable :: momtem_riv(:)
      
      REAL(r8), allocatable :: height_riv_ds(:)
      REAL(r8), allocatable :: veloct_riv_ds(:)
      REAL(r8), allocatable :: momtem_riv_ds(:)

      REAL(r8), allocatable :: hflux_fc(:)
      REAL(r8), allocatable :: mflux_fc(:)
      REAL(r8), allocatable :: zgrad_fc(:)
      
      REAL(r8), allocatable :: sum_hflux_riv(:)
      REAL(r8), allocatable :: sum_mflux_riv(:)
      REAL(r8), allocatable :: sum_zgrad_riv(:)
      
      REAL(r8) :: veloct_fc, height_fc, momtem_fc, zsurf_fc
      REAL(r8) :: riverelv_fc, height_up, height_dn
      REAL(r8) :: vwave_up, vwave_dn, hflux_up, hflux_dn, mflux_up, mflux_dn
      REAL(r8) :: loss, friction
      REAL(r8) :: dt_res, dt_this
      CHARACTER(len=50) :: fmtt


      IF (p_is_worker) THEN
      
         nriver = numelm

         IF (nriver > 0) THEN
            allocate (height_riv (nriver))
            allocate (veloct_riv (nriver))
            allocate (momtem_riv (nriver))
            
            allocate (height_riv_ds (nriver))
            allocate (veloct_riv_ds (nriver))
            allocate (momtem_riv_ds (nriver))
            
            allocate (hflux_fc (nriver))
            allocate (mflux_fc (nriver))
            allocate (zgrad_fc (nriver))
            
            allocate (sum_hflux_riv (nriver))
            allocate (sum_mflux_riv (nriver))
            allocate (sum_zgrad_riv (nriver))
         ENDIF

         DO i = 1, nriver
            istt = hru_patch%substt(basin_hru%substt(i))
            iend = hru_patch%subend(basin_hru%substt(i))
            height_riv(i) = sum(dpond(istt:iend) * hru_patch%subfrc(istt:iend))
            height_riv(i) = height_riv(i) / 1.0e3 ! mm to m
         ENDDO

         DO i = 1, nriver
            veloct_riv(i) = riverveloct(i)
            ! IF water is increased, momentum is not increased.
            ! IF water is decreased, momentum is also decreased.
            momtem_riv(i) = min(riverheight(i), height_riv(i)) * veloct_riv(i)
         ENDDO

         dt_res = dt
         DO WHILE (dt_res > 0)
            
            DO i = 1, nriver
               sum_hflux_riv(i) = 0.
               sum_mflux_riv(i) = 0.
               sum_zgrad_riv(i) = 0.
               
               IF (addrdown(i) > 0) THEN
                  height_riv_ds(i) = height_riv(addrdown(i))
                  veloct_riv_ds(i) = veloct_riv(addrdown(i))
                  momtem_riv_ds(i) = momtem_riv(addrdown(i))
               ELSE
                  height_riv_ds(i) = 0
                  veloct_riv_ds(i) = 0
                  momtem_riv_ds(i) = 0
               ENDIF
            ENDDO 
#ifdef USEMPI
            CALL river_data_exchange (SEND_DATA_DOWN_TO_UP, accum = .false., &
               vec_send1 = height_riv, vec_recv1 = height_riv_ds, &
               vec_send2 = veloct_riv, vec_recv2 = veloct_riv_ds, &
               vec_send3 = momtem_riv, vec_recv3 = momtem_riv_ds)
#endif
               
            dt_this = dt_res

            DO i = 1, nriver
               IF (addrdown(i) >= 0) THEN

                  IF ((height_riv(i) < RIVERMIN) .and. (height_riv_ds(i) < RIVERMIN)) THEN
                     hflux_fc(i) = 0
                     mflux_fc(i) = 0
                     zgrad_fc(i) = 0
                     cycle
                  ENDIF
                  
                  ! reconstruction of height of water near interface
                  riverelv_fc = max(riverelv(i), riverelv_ds(i))
                  height_up = max(0., height_riv(i)   +riverelv(i)   -riverelv_fc) 
                  height_dn = max(0., height_riv_ds(i)+riverelv_ds(i)-riverelv_fc) 

                  ! velocity at river downstream face
                  veloct_fc = 0.5 * (veloct_riv(i) + veloct_riv_ds(i)) & 
                     + sqrt(grav * height_up) - sqrt(grav * height_dn)

                  ! height of water at downstream face
                  height_fc = 1/grav * (0.5*(sqrt(grav*height_up) + sqrt(grav*height_dn)) &
                     + 0.25 * (veloct_riv(i) - veloct_riv_ds(i))) ** 2

                  IF (height_up > 0) THEN
                     vwave_up = min(veloct_riv(i)-sqrt(grav*height_up), veloct_fc-sqrt(grav*height_fc))
                  ELSE
                     vwave_up = veloct_riv_ds(i) - 2.0 * sqrt(grav*height_dn)
                  ENDIF

                  IF (height_dn > 0) THEN
                     vwave_dn = max(veloct_riv_ds(i)+sqrt(grav*height_dn), veloct_fc+sqrt(grav*height_fc))
                  ELSE
                     vwave_dn = veloct_riv(i) + 2.0 * sqrt(grav*height_up)
                  ENDIF

                  hflux_up = veloct_riv(i)    * height_up
                  hflux_dn = veloct_riv_ds(i) * height_dn
                  mflux_up = veloct_riv(i)**2    * height_up + 0.5*grav * height_up**2 
                  mflux_dn = veloct_riv_ds(i)**2 * height_dn + 0.5*grav * height_dn**2 

                  IF (vwave_up >= 0.) THEN
                     hflux_fc(i) = riverfac_ds(i) * hflux_up
                     mflux_fc(i) = riverfac_ds(i) * mflux_up
                  ELSEIF (vwave_dn <= 0.) THEN
                     hflux_fc(i) = riverfac_ds(i) * hflux_dn
                     mflux_fc(i) = riverfac_ds(i) * mflux_dn
                  ELSE
                     hflux_fc(i) = riverfac_ds(i) * (vwave_dn*hflux_up - vwave_up*hflux_dn &
                        + vwave_up*vwave_dn*(height_dn-height_up)) / (vwave_dn-vwave_up)
                     mflux_fc(i) = riverfac_ds(i) * (vwave_dn*mflux_up - vwave_up*mflux_dn &
                        + vwave_up*vwave_dn*(hflux_dn-hflux_up)) / (vwave_dn-vwave_up)
                  ENDIF
               
                  sum_zgrad_riv(i) = sum_zgrad_riv(i) - riverfac_ds(i) * 0.5*grav * height_up**2

                  zgrad_fc(i) = riverfac_ds(i) * 0.5*grav * height_dn**2
                  
               ELSE
                  IF (riverdown(i) == 0) THEN
                     ! river mouth
                     veloct_fc = max(0., veloct_riv(i))
                     hflux_fc(i) = height_riv(i) * veloct_fc * riverwdth(i)
                     mflux_fc(i) = momtem_riv(i) * veloct_fc * riverwdth(i)
                     zgrad_fc(i) = 0
                  ELSE 
                     ! inland depression
                     hflux_fc(i) = 0
                     mflux_fc(i) = 0
                     zgrad_fc(i) = 0
                  ENDIF
               ENDIF

               sum_hflux_riv(i) = sum_hflux_riv(i) + hflux_fc(i)
               sum_mflux_riv(i) = sum_mflux_riv(i) + mflux_fc(i)

               IF (addrdown(i) > 0) THEN
                  j = addrdown(i)
                  sum_hflux_riv(j) = sum_hflux_riv(j) - hflux_fc(i)
                  sum_mflux_riv(j) = sum_mflux_riv(j) - mflux_fc(i)
                  sum_zgrad_riv(j) = sum_zgrad_riv(j) + zgrad_fc(i)
               ENDIF
            
            ENDDO
            
#ifdef USEMPI
            hflux_fc = - hflux_fc
            mflux_fc = - mflux_fc
            CALL river_data_exchange (SEND_DATA_UP_TO_DOWN, accum = .true., &
               vec_send1 = hflux_fc, vec_recv1 = sum_hflux_riv, &
               vec_send2 = mflux_fc, vec_recv2 = sum_mflux_riv, &
               vec_send3 = zgrad_fc, vec_recv3 = sum_zgrad_riv)
#endif

            DO i = 1, nriver
               ! CFL condition
               IF ((veloct_riv(i) /= 0.) .or. (height_riv(i) > 0.)) THEN
                  dt_this = min(dt_this, riverlen(i)/(abs(veloct_riv(i))+sqrt(grav*height_riv(i)))*0.8)
               ENDIF

               loss = sum_hflux_riv(i) / riverarea(i)
               IF (loss > 0) THEN
                  dt_this = min(dt_this, height_riv(i) / loss)
               ENDIF

               IF ((abs(veloct_riv(i)) > 0.1) &
                  .and. (veloct_riv(i) * (sum_mflux_riv(i)+sum_zgrad_riv(i)) > 0)) THEN
                  dt_this = min(dt_this, &
                     abs(momtem_riv(i) * riverarea(i) / (sum_mflux_riv(i)+sum_zgrad_riv(i))))
               ENDIF
            ENDDO 

#ifdef USEMPI
            CALL mpi_allreduce (MPI_IN_PLACE, dt_this, 1, MPI_REAL8, MPI_MIN, p_comm_worker, p_err)
#endif

            DO i = 1, nriver

               height_riv(i) = height_riv(i) - sum_hflux_riv(i) / riverarea(i) * dt_this

               IF (height_riv(i) < RIVERMIN) THEN
                  momtem_riv(i) = 0
                  veloct_riv(i) = 0
               ELSE
                  friction = grav * nmanning_riv**2 / height_riv(i)**(7.0/3.0) * abs(momtem_riv(i))
                  momtem_riv(i) = (momtem_riv(i) &
                     - (sum_mflux_riv(i) + sum_zgrad_riv(i)) / riverarea(i) * dt_this) &
                     / (1 + friction * dt_this) 
                  veloct_riv(i) = momtem_riv(i) / height_riv(i)
               ENDIF
            
               ! inland depression
               IF (riverdown(i) < 0) THEN
                  momtem_riv(i) = min(0., momtem_riv(i))
                  veloct_riv(i) = min(0., veloct_riv(i))
               ENDIF

            ENDDO

            IF (nriver > 0) THEN
               riverheight_ta(:) = riverheight_ta(:) + height_riv(:) * dt_this
               rivermomtem_ta(:) = rivermomtem_ta(:) + momtem_riv(:) * dt_this
            ENDIF

            dt_res = dt_res - dt_this

         ENDDO

         DO i = 1, nriver
         
            riverheight(i) = height_riv(i)
            riverveloct(i) = veloct_riv(i)

            istt = hru_patch%substt(basin_hru%substt(i))
            iend = hru_patch%subend(basin_hru%substt(i))
            dpond(istt:iend) = height_riv(i) * 1.0e3 ! m to mm

         ENDDO
            
         IF (allocated(height_riv)) deallocate(height_riv)
         IF (allocated(veloct_riv)) deallocate(veloct_riv)
         IF (allocated(momtem_riv)) deallocate(momtem_riv)

         IF (allocated(height_riv_ds)) deallocate(height_riv_ds)
         IF (allocated(veloct_riv_ds)) deallocate(veloct_riv_ds)
         IF (allocated(momtem_riv_ds)) deallocate(momtem_riv_ds)

         IF (allocated(hflux_fc)) deallocate(hflux_fc)
         IF (allocated(mflux_fc)) deallocate(mflux_fc)
         IF (allocated(zgrad_fc)) deallocate(zgrad_fc)

         IF (allocated(sum_hflux_riv)) deallocate(sum_hflux_riv)
         IF (allocated(sum_mflux_riv)) deallocate(sum_mflux_riv)
         IF (allocated(sum_zgrad_riv)) deallocate(sum_zgrad_riv)

      ENDIF

   END SUBROUTINE river_flow
   
END MODULE MOD_Hydro_RiverFlow
#endif
