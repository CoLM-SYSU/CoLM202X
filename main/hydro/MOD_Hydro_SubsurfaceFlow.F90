#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_Hydro_SubsurfaceFlow

   USE MOD_Precision
   IMPLICIT NONE
    
   REAL(r8), parameter :: e_ice  = 6.0   ! soil ice impedance factor
   REAL(r8), parameter :: raniso = 1.   ! anisotropy ratio, unitless
   
CONTAINS
   
   ! ---------
   SUBROUTINE subsurface_runoff (deltime)
      
      USE MOD_SPMD_Task
      USE MOD_Mesh
      USE MOD_LandElm
      USE MOD_LandPatch
      USE MOD_Vars_TimeVariables
      USE MOD_Vars_TimeInvariants
      USE MOD_Vars_1DFluxes
      USE MOD_Hydro_DrainageNetwork
      USE MOD_Hydro_RiverNetwork
      USE MOD_Hydro_SubsurfaceNetwork
      USE MOD_Const_Physical, only : denice, denh2o
      USE MOD_Vars_Global, only : pi

      IMPLICIT NONE
      
      REAL(r8), intent(in) :: deltime

      ! Local Variables
      INTEGER :: numbasin, nhru, ibasin, i, j, ihru, ipatch, ps, pe, ilev

      TYPE(drainage_network_info_type), pointer :: hrus

      REAL(r8), allocatable :: theta_a_h (:) 
      REAL(r8), allocatable :: zwt_h     (:) 
      REAL(r8), allocatable :: Ks_h      (:) ! [m/s]
      REAL(r8), allocatable :: rsubs_h   (:) ! [m/s]
      REAL(r8), allocatable :: rsubs_fc  (:) ! [m/s]
      REAL(r8) :: rsubs_bsn

      REAL(r8) :: theta_s_h, air_h, icefrac, imped, delp
      REAL(r8) :: zsubs_h_up, zsubs_h_dn
      REAL(r8) :: slope, bdamp, Ks_fc, Ks_in
      REAL(r8) :: ca, cb
      REAL(r8) :: alp
      
      REAL(r8), allocatable :: theta_a_bsn (:) 
      REAL(r8), allocatable :: zwt_bsn     (:) 
      REAL(r8), allocatable :: Ks_bsn      (:) ! [m/s]

      INTEGER  :: jnb
      REAL(r8) :: zsubs_up, zwt_up, Ks_up, theta_a_up, area_up
      REAL(r8) :: zsubs_dn, zwt_dn, Ks_dn, theta_a_dn, area_dn
      REAL(r8) :: lenbdr, rsubs_nb

      IF (p_is_worker) THEN

         numbasin = numelm
            
         rsubs_pch(:) = 0.  ! subsurface runoff in each patch
         rsub     (:) = 0.  ! subsurface runoff into river network
         rsubs_hru(:) = 0.  ! subsurface runoff in each hydro unit

         bdamp = 4.8

         IF (numbasin > 0) THEN
            allocate (theta_a_bsn (numbasin))
            allocate (zwt_bsn     (numbasin))
            allocate (Ks_bsn      (numbasin))
         ENDIF

         DO ibasin = 1, numbasin

            hrus => drainagenetwork(ibasin)
            
            theta_a_bsn (ibasin) = 0.
            zwt_bsn     (ibasin) = 0.
            Ks_bsn      (ibasin) = 0.
            
            nhru = hrus%nhru
            IF (nhru <= 1) THEN
               ps = hru_patch%substt(hrus%ihru(1))
               pe = hru_patch%subend(hrus%ihru(1))
               zwt_hru(hrus%ihru(1)) = sum(zwt(ps:pe) * hru_patch%subfrc(ps:pe))
               cycle
            ENDIF
            
            allocate (theta_a_h (nhru))
            allocate (zwt_h     (nhru))
            allocate (Ks_h      (nhru))

            theta_a_h = 0.
            zwt_h     = 0.
            Ks_h      = 0.

            DO i = 2, nhru
               ps = hru_patch%substt(hrus%ihru(i))
               pe = hru_patch%subend(hrus%ihru(i))

               theta_s_h = 0
               DO ipatch = ps, pe
                  theta_s_h = theta_s_h + hru_patch%subfrc(ipatch) &
                     * sum(porsl(1:nl_soil,ipatch) * dz_soi(1:nl_soil) &
                      - wice_soisno(1:nl_soil,ipatch)/denice) / sum(dz_soi(1:nl_soil)) 
               ENDDO

               IF (theta_s_h > 0.) THEN
                  
                  zwt_h(i) = sum(zwt(ps:pe) * hru_patch%subfrc(ps:pe))

                  air_h = 0.
                  DO ipatch = ps, pe
                     air_h = air_h + hru_patch%subfrc(ipatch) &
                        * (sum( porsl(1:nl_soil,ipatch) * dz_soi(1:nl_soil) &
                        - wliq_soisno(1:nl_soil,ipatch)/denh2o &
                        - wice_soisno(1:nl_soil,ipatch)/denice ) - wa(ipatch)/1.0e3)
                     air_h = max(0., air_h)
                  ENDDO

                  IF ((air_h <= 0.) .or. (zwt_h(i) <= 0.)) THEN
                     theta_a_h(i) = theta_s_h
                     zwt_h(i) = 0.
                  ELSE
                     theta_a_h(i) = air_h / zwt_h(i)
                     IF (theta_a_h(i) > theta_s_h) THEN
                        theta_a_h(i) = theta_s_h
                        zwt_h(i) = air_h / theta_a_h(i)
                     ENDIF
                  ENDIF

                  Ks_h(i) = 0.
                  DO ipatch = ps, pe
                     DO ilev = 1, nl_soil
                        icefrac = min(1., wice_soisno(ilev,ipatch)/denice/dz_soi(ilev)/porsl(ilev,ipatch))
                        imped   = 10.**(-e_ice*icefrac)
                        Ks_h(i) = Ks_h(i) + hru_patch%subfrc(ipatch) &
                           * hksati(ilev,ipatch)/1.0e3 * imped * dz_soi(ilev)/zi_soi(nl_soil) 
                     ENDDO
                  ENDDO
               ELSE
                  ! Frozen soil.
                  Ks_h(i) = 0.
               ENDIF

            ENDDO
            
            allocate (rsubs_h  (nhru))
            allocate (rsubs_fc (nhru))

            rsubs_h (:) = 0.
            rsubs_fc(:) = 0.

            DO i = 2, nhru
                  
               j = hrus%inext(i)

               IF (j == 1) THEN
                  IF (Ks_h(i) == 0.) cycle
               ELSE
                  IF ((Ks_h(i) == 0.) .or. (Ks_h(j) == 0.)) THEN
                     cycle
                  ENDIF
               ENDIF
                  
               zsubs_h_up = hrus%elva(i) - zwt_h(i)

               IF (j > 1) THEN
                  zsubs_h_dn = hrus%elva(j) - zwt_h(j)
               ELSE
                  zsubs_h_dn = hrus%elva(1) - riverdpth(ibasin) + dpond_hru(hrus%ihru(1)) 
               ENDIF

               IF (j > 1) THEN
                  delp = hrus%plen(i) + hrus%plen(j)
               ELSE
                  delp = hrus%plen(i)
               ENDIF

               slope = abs(hrus%elva(i)-hrus%elva(j))/delp
               ! from Fan et al., JGR 112(D10125)
               IF (slope > 0.16) THEN
                  bdamp = 4.8
               ELSE
                  bdamp = 120./(1+150.*slope)
               ENDIF

               ! Upstream scheme for hydraulic conductivity
               IF ((zsubs_h_up > zsubs_h_dn) .or. (j == 1)) THEN
                  IF (zwt_h(i) > 1.5) THEN
                     ! from Fan et al., JGR 112(D10125)
                     Ks_fc = raniso * Ks_h(i) * bdamp * exp(-(zwt_h(i)-1.5)/bdamp)
                  ELSE
                     Ks_fc = raniso * Ks_h(i) * ((1.5-zwt_h(i)) + bdamp)
                  ENDIF
               ELSE
                  IF (zwt_h(j) > 1.5) THEN
                     Ks_fc = raniso * Ks_h(j) * bdamp * exp(-(zwt_h(j)-1.5)/bdamp)
                  ELSE
                     Ks_fc = raniso * Ks_h(j) * ((1.5-zwt_h(j)) + bdamp)
                  ENDIF
               ENDIF

               ca = hrus%flen(i) * Ks_fc / theta_a_h(i) / delp / hrus%area(i) * deltime

               IF (j > 1) THEN
                  cb = hrus%flen(i) * Ks_fc / theta_a_h(j) / delp / hrus%area(j) * deltime
               ELSE
                  cb = hrus%flen(i) * Ks_fc / delp / hrus%area(j) * deltime
               ENDIF
               
               rsubs_fc(i) = (zsubs_h_up - zsubs_h_dn) * hrus%flen(i) * Ks_fc / (1+ca+cb) / delp

               rsubs_h(i) = rsubs_h(i) + rsubs_fc(i) / hrus%area(i)
               rsubs_h(j) = rsubs_h(j) - rsubs_fc(i) / hrus%area(j)

            ENDDO
            
            rsubs_bsn = - rsubs_h(1) * hrus%area(1)/sum(hrus%area) * 1.0e3 ! (positive = out of soil column) 

            IF (rsubs_h(1)*deltime > riverheight(ibasin)*riverarea(ibasin)) THEN 
               alp = riverheight(ibasin)*riverarea(ibasin) / (rsubs_h(1)*deltime)
               rsubs_bsn  = rsubs_bsn  * alp
               rsubs_h(1) = rsubs_h(1) * alp
               DO i = 2, nhru
                  j = hrus%inext(i)
                  IF (j == 1) THEN
                     rsubs_h(i) = rsubs_h(i) - (1.0-alp)*rsubs_fc(i)/hrus%area(i)
                  ENDIF
               ENDDO
            ENDIF
            
            DO i = 1, nhru
               ps = hru_patch%substt(hrus%ihru(i))
               pe = hru_patch%subend(hrus%ihru(i))
               
               ! Between hydrological units
               rsubs_pch(ps:pe) = rsubs_h(i) ! (positive = out of soil column) 

               ! Inside hydrological units
               IF (i > 1) THEN
                  IF (zwt_h(i) > 1.5) THEN
                     ! from Fan et al., JGR 112(D10125)
                     Ks_in = raniso * Ks_h(i) * bdamp * exp(-(zwt_h(i)-1.5)/bdamp)
                  ELSE
                     Ks_in = raniso * Ks_h(i) * ((1.5-zwt_h(i)) + bdamp)
                  ENDIF

                  rsubs_pch(ps:pe) = rsubs_pch(ps:pe) &
                     - Ks_in * (zwt(ps:pe) - sum(zwt(ps:pe)*hru_patch%subfrc(ps:pe))) *6.0*pi/hrus%area(i)
               ENDIF

               rsubs_pch(ps:pe) = rsubs_pch(ps:pe) * 1.e3 ! m/s to mm/s

               zwt_hru(hrus%ihru(i)) = zwt_h(i)
               rsubs_hru(hrus%ihru(i)) = rsubs_h(i)

               rsub(ps:pe) = rsubs_bsn 
            ENDDO
            
            theta_a_bsn (ibasin) = sum(theta_a_h(2:) * hrus%area(2:)) / sum(hrus%area(2:))
            zwt_bsn     (ibasin) = sum(zwt_h    (2:) * hrus%area(2:)) / sum(hrus%area(2:))
            Ks_bsn      (ibasin) = sum(Ks_h     (2:) * hrus%area(2:)) / sum(hrus%area(2:))

            deallocate (theta_a_h)
            deallocate (zwt_h    )
            deallocate (Ks_h     )
            deallocate (rsubs_h  )
            deallocate (rsubs_fc )

         ENDDO

         CALL retrieve_neighbour_data (theta_a_bsn, theta_a_nb)
         CALL retrieve_neighbour_data (zwt_bsn    , zwt_nb    )
         CALL retrieve_neighbour_data (Ks_bsn     , Ks_nb     )

         DO ibasin = 1, numbasin
            
            hrus => drainagenetwork(ibasin)
            
            DO jnb = 1, num_nb(ibasin)

               Ks_up = Ks_bsn(ibasin)
               Ks_dn = Ks_nb(ibasin)%val(jnb)

               IF ((Ks_up == 0.) .or. (Ks_dn == 0.)) THEN
                  cycle
               ENDIF

               zwt_up = zwt_bsn(ibasin)
               zwt_dn = zwt_nb(ibasin)%val(jnb)

               theta_a_up = theta_a_bsn(ibasin)
               theta_a_dn = theta_a_nb(ibasin)%val(jnb)

               zsubs_up = elva_b(ibasin) - zwt_up 
               zsubs_dn = elva_nb(ibasin)%val(jnb) - zwt_dn

               delp = dist_nb(ibasin)%val(jnb)

               area_up = area_b(ibasin)
               area_dn = area_nb(ibasin)%val(jnb)

               lenbdr = lenbdr_nb(ibasin)%val(jnb)

               ! from Fan et al., JGR 112(D10125)
               slope = slope_nb(ibasin)%val(jnb)
               IF (slope > 0.16) THEN
                  bdamp = 4.8
               ELSE
                  bdamp = 120./(1+150.*slope)
               ENDIF

               ! Upstream scheme for hydraulic conductivity
               IF (zsubs_up > zsubs_dn) THEN
                  IF (zwt_up > 1.5) THEN
                     ! from Fan et al., JGR 112(D10125)
                     Ks_fc = raniso * Ks_up * bdamp * exp(-(zwt_up-1.5)/bdamp)
                  ELSE
                     Ks_fc = raniso * Ks_up * ((1.5-zwt_up) + bdamp)
                  ENDIF
               ELSE
                  IF (zwt_dn > 1.5) THEN
                     Ks_fc = raniso * Ks_dn * bdamp * exp(-(zwt_dn-1.5)/bdamp)
                  ELSE
                     Ks_fc = raniso * Ks_dn * ((1.5-zwt_dn) + bdamp)
                  ENDIF
               ENDIF

               ca = lenbdr * Ks_fc / theta_a_up / delp / area_up * deltime
               cb = lenbdr * Ks_fc / theta_a_dn / delp / area_dn * deltime

               rsubs_nb = (zsubs_up - zsubs_dn) * lenbdr * Ks_fc / (1+ca+cb) / delp
               rsubs_nb = rsubs_nb / sum(hrus%area(2:))

               DO i = 2, hrus%nhru
                  ps = hru_patch%substt(hrus%ihru(i))
                  pe = hru_patch%subend(hrus%ihru(i))
               
                  rsubs_pch(ps:pe) = rsubs_pch(ps:pe) + rsubs_nb * 1.e3 ! m/s to mm/s
                  rsubs_hru(hrus%ihru(i)) = rsubs_hru(hrus%ihru(i)) + rsubs_nb
               ENDDO

            ENDDO
         ENDDO

         IF (allocated(theta_a_bsn)) deallocate(theta_a_bsn)
         IF (allocated(zwt_bsn    )) deallocate(zwt_bsn    )
         IF (allocated(Ks_bsn     )) deallocate(Ks_bsn     )

      ENDIF


   END SUBROUTINE subsurface_runoff

END MODULE MOD_Hydro_SubsurfaceFlow
#endif
