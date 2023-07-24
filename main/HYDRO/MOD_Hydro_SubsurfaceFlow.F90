#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_Hydro_SubsurfaceFlow
   !-------------------------------------------------------------------------------------
   ! DESCRIPTION:
   !   
   !   Ground water lateral flow.
   !
   !   Ground water fluxes are calculated
   !   1. between basins
   !   2. between hydrological response units
   !   3. between patches inside one HRU  
   !
   ! Created by Shupeng Zhang, May 2023
   !-------------------------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE
    
   REAL(r8), parameter :: e_ice  = 6.0   ! soil ice impedance factor
   REAL(r8), parameter :: raniso = 1.    ! anisotropy ratio, unitless
   
CONTAINS
   
   ! ---------
   SUBROUTINE subsurface_flow (deltime)
      
      USE MOD_SPMD_Task
      USE MOD_Mesh
      USE MOD_LandElm
      USE MOD_LandPatch
      USE MOD_Vars_TimeVariables
      USE MOD_Vars_TimeInvariants
      USE MOD_Vars_1DFluxes
      USE MOD_Hydro_SurfaceNetwork
      USE MOD_Hydro_RiverNetwork
      USE MOD_Hydro_SubsurfaceNetwork
      USE MOD_Const_Physical,  only : denice, denh2o
      USE MOD_Vars_Global,     only : pi, nl_soil, zi_soi
      USE MOD_Hydro_SoilWater, only : soilwater_aquifer_exchange

      IMPLICIT NONE
      
      REAL(r8), intent(in) :: deltime

      ! Local Variables
      INTEGER :: numbasin, nhru, ibasin, i, j, ihru, ipatch, ps, pe, ilev

      TYPE(surface_network_info_type), pointer :: hrus

      REAL(r8), allocatable :: theta_a_h (:) 
      REAL(r8), allocatable :: zwt_h     (:) 
      REAL(r8), allocatable :: Ks_h      (:) ! [m/s]
      REAL(r8), allocatable :: rsubs_h   (:) ! [m/s]
      REAL(r8), allocatable :: rsubs_fc  (:) ! [m/s]
      REAL(r8) :: rsubs_riv

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

      ! for water exchange
      integer  :: izwt
      real(r8) :: exwater
      REAL(r8) :: sp_zi(0:nl_soil), sp_dz(1:nl_soil), zwtmm ! [mm]
      real(r8) :: vl_r (1:nl_soil)
#ifdef Campbell_SOIL_MODEL
      INTEGER, parameter :: nprms = 1
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      INTEGER, parameter :: nprms = 5
#endif
      REAL(r8) :: prms  (nprms,1:nl_soil)
      real(r8) :: vol_ice     (1:nl_soil)
      real(r8) :: vol_liq     (1:nl_soil)
      real(r8) :: eff_porosity(1:nl_soil)
      logical  :: is_permeable(1:nl_soil)
      real(r8) :: wresi       (1:nl_soil)
      real(r8) :: w_sum_before, w_sum_after, errblc


      IF (p_is_worker) THEN

         numbasin = numelm
            
         rsubs_bsn(:) = 0.  ! subsurface lateral flow between basins                     
         rsubs_hru(:) = 0.  ! subsurface lateral flow between hydrological response units
         rsubs_pch(:) = 0.  ! subsurface lateral flow between patches inside one HRU     

         rsub(:) = 0. ! total recharge/discharge from subsurface lateral flow

         bdamp = 4.8

         IF (numbasin > 0) THEN
            allocate (theta_a_bsn (numbasin))
            allocate (zwt_bsn     (numbasin))
            allocate (Ks_bsn      (numbasin))
         ENDIF

         DO ibasin = 1, numbasin

            hrus => surface_network(ibasin)
            
            theta_a_bsn (ibasin) = 0.
            zwt_bsn     (ibasin) = 0.
            Ks_bsn      (ibasin) = 0.
            
            nhru = hrus%nhru

            IF (nhru <= 1) THEN
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
                  zsubs_h_dn = hrus%elva(1) - riverdpth(ibasin) + wdsrf_hru(hrus%ihru(1)) 
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
            
            rsubs_riv = - rsubs_h(1) * hrus%area(1)/sum(hrus%area) * 1.0e3 ! (positive = out of soil column) 

            IF (rsubs_h(1)*deltime > riverheight(ibasin)*riverarea(ibasin)) THEN 
               alp = riverheight(ibasin)*riverarea(ibasin) / (rsubs_h(1)*deltime)
               rsubs_riv  = rsubs_riv  * alp
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
               
               ! Update total subsurface lateral flow (1): Between hydrological units
               rsub(ps:pe) = rsub(ps:pe) + rsubs_h(i) * 1.e3 ! (positive = out of soil column) 

               ! Inside hydrological units
               IF (i > 1) THEN
                  IF (zwt_h(i) > 1.5) THEN
                     ! from Fan et al., JGR 112(D10125)
                     Ks_in = raniso * Ks_h(i) * bdamp * exp(-(zwt_h(i)-1.5)/bdamp)
                  ELSE
                     Ks_in = raniso * Ks_h(i) * ((1.5-zwt_h(i)) + bdamp)
                  ENDIF

                  rsubs_pch(ps:pe) = &
                     - Ks_in * (zwt(ps:pe) - sum(zwt(ps:pe)*hru_patch%subfrc(ps:pe))) *6.0*pi/hrus%area(i)
               
                  ! Update total subsurface lateral flow (2): Between patches
                  rsub(ps:pe) = rsub(ps:pe) + rsubs_pch(ps:pe) * 1.e3 ! m/s to mm/s

               ENDIF

               rsubs_hru(hrus%ihru(i)) = rsubs_h(i)

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
            
            hrus => surface_network(ibasin)
            
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
               
               rsubs_bsn(ibasin) = rsubs_bsn(ibasin) + rsubs_nb

            ENDDO

            DO i = 2, hrus%nhru
               ps = hru_patch%substt(hrus%ihru(i))
               pe = hru_patch%subend(hrus%ihru(i))

               ! Update total subsurface lateral flow (3): Between basins
               rsub(ps:pe) = rsub(ps:pe) + rsubs_bsn(ibasin) * 1.e3 ! m/s to mm/s
            ENDDO

         ENDDO

         IF (allocated(theta_a_bsn)) deallocate(theta_a_bsn)
         IF (allocated(zwt_bsn    )) deallocate(zwt_bsn    )
         IF (allocated(Ks_bsn     )) deallocate(Ks_bsn     )

      ENDIF

      ! Exchange between soil water and aquifer.
      IF (p_is_worker) THEN
      
         sp_zi(0) = 0.
         sp_zi(1:nl_soil) = zi_soi(1:nl_soil) * 1000.0   ! from meter to mm
         sp_dz(1:nl_soil) = sp_zi(1:nl_soil) - sp_zi(0:nl_soil-1)
         
         DO ipatch = 1, numpatch

#if(defined CoLMDEBUG)
            ! For water balance check, the sum of water in soil column before the calcultion
            w_sum_before = sum(wliq_soisno(1:nl_soil,ipatch)) + sum(wice_soisno(1:nl_soil,ipatch)) &
               + wa(ipatch) + wdsrf(ipatch)
#endif

            exwater = rsub(ipatch) * deltime

#ifdef Campbell_SOIL_MODEL
            vl_r(1:nl_soil) = 0._r8
            prms(1,:) = bsw(1:nl_soil,ipatch)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
            vl_r  (1:nl_soil) = theta_r  (1:nl_soil,ipatch)
            prms(1,1:nl_soil) = alpha_vgm(1:nl_soil,ipatch)
            prms(2,1:nl_soil) = n_vgm    (1:nl_soil,ipatch)
            prms(3,1:nl_soil) = L_vgm    (1:nl_soil,ipatch)
            prms(4,1:nl_soil) = sc_vgm   (1:nl_soil,ipatch)
            prms(5,1:nl_soil) = fc_vgm   (1:nl_soil,ipatch)
#endif

            DO ilev = 1, nl_soil
               vol_ice(ilev) = wice_soisno(ilev,ipatch)/denice*1000. / sp_dz(ilev)
               vol_ice(ilev) = min(vol_ice(ilev), porsl(ilev,ipatch))

               eff_porosity(ilev) = max(wimp, porsl(ilev,ipatch)-vol_ice(ilev))
               is_permeable(ilev) = eff_porosity(ilev) > max(wimp, vl_r(ilev))
               IF (is_permeable(ilev)) THEN
                  vol_liq(ilev) = wliq_soisno(ilev,ipatch)/denh2o*1000. / sp_dz(ilev)
                  vol_liq(ilev) = min(eff_porosity(ilev), max(0., vol_liq(ilev)))
                  wresi(ilev) = wliq_soisno(ilev,ipatch) - sp_dz(ilev)*vol_liq(ilev)/1000. * denh2o
               ENDIF
            ENDDO

            zwtmm = zwt(ipatch) * 1000. ! m -> mm
            izwt  = findloc(zwtmm >= sp_zi, .true., dim=1, back=.true.)
            
            IF (izwt <= nl_soil) THEN
               IF (is_permeable(izwt)) THEN
                  vol_liq(izwt) = (wliq_soisno(izwt,ipatch)/denh2o*1000.0 &
                     - eff_porosity(izwt)*(sp_zi(izwt)-zwtmm)) / (zwtmm - sp_zi(izwt-1))

                  IF (vol_liq(izwt) < 0.) THEN
                     zwtmm = sp_zi(izwt)
                     vol_liq(izwt) = wliq_soisno(izwt,ipatch)/denh2o*1000.0 / (sp_zi(izwt)-sp_zi(izwt-1))
                  ENDIF

                  vol_liq(izwt) = max(0., min(eff_porosity(izwt), vol_liq(izwt)))
                  wresi(izwt) = wliq_soisno(izwt,ipatch) - (eff_porosity(izwt)*(sp_zi(izwt)-zwtmm) &
                     + vol_liq(izwt)*(zwtmm-sp_zi(izwt-1))) /1000. * denh2o
               ENDIF
            ENDIF

            CALL soilwater_aquifer_exchange ( &
               nl_soil, exwater, sp_zi, is_permeable, porsl(:,ipatch), vl_r, psi0(:,ipatch), &
               hksati(:,ipatch), nprms, prms,   porsl(nl_soil,ipatch), wdsrf(ipatch), &
               vol_liq, zwtmm, wa(ipatch), izwt)
      
            ! update the mass of liquid water
            DO ilev = nl_soil, 1, -1
               IF (is_permeable(ilev)) THEN
                  IF (zwtmm < sp_zi(ilev)) THEN
                     IF (zwtmm >= sp_zi(ilev-1)) THEN
                        wliq_soisno(ilev,ipatch) = ((eff_porosity(ilev)*(sp_zi(ilev)-zwtmm))  &
                           + vol_liq(ilev)*(zwtmm-sp_zi(ilev-1)))/1000.0 * denh2o
                     ELSE
                        wliq_soisno(ilev,ipatch) = denh2o * eff_porosity(ilev)*sp_dz(ilev)/1000.0
                     ENDIF
                  ELSE
                     wliq_soisno(ilev,ipatch) = denh2o * vol_liq(ilev)*sp_dz(ilev)/1000.0
                  ENDIF

                  wliq_soisno(ilev,ipatch) = wliq_soisno(ilev,ipatch) + wresi(ilev)
               ENDIF
            ENDDO

            zwt(ipatch) = zwtmm/1000.0

#if(defined CoLMDEBUG)
            ! For water balance check, the sum of water in soil column after the calcultion
            w_sum_after = sum(wliq_soisno(1:nl_soil,ipatch)) + sum(wice_soisno(1:nl_soil,ipatch)) &
               + wa(ipatch) + wdsrf(ipatch)
            errblc = w_sum_after - w_sum_before + exwater

            if(abs(errblc) > 1.e-3)then
               write(6,'(A,E20.5)') 'Warning (Subsurface Runoff): water balance violation', errblc
            endif
#endif

         ENDDO
      ENDIF

   END SUBROUTINE subsurface_flow

END MODULE MOD_Hydro_SubsurfaceFlow
#endif
