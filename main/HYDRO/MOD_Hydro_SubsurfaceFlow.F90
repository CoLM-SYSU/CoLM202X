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
      USE MOD_Hydro_HillslopeNetwork
      USE MOD_Hydro_RiverLakeNetwork
      USE MOD_Hydro_BasinNeighbour
      USE MOD_Const_Physical,  only : denice, denh2o
      USE MOD_Vars_Global,     only : pi, nl_soil, zi_soi
      USE MOD_Hydro_SoilWater, only : soilwater_aquifer_exchange

      IMPLICIT NONE
      
      REAL(r8), intent(in) :: deltime

      ! Local Variables
      INTEGER :: numbasin, nhru, ibasin, i, i0, j, ihru, ipatch, ps, pe, ilev

      TYPE(hillslope_network_info_type), pointer :: hrus

      REAL(r8), allocatable :: theta_a_h (:) 
      REAL(r8), allocatable :: zwt_h     (:) 
      REAL(r8), allocatable :: Ks_h      (:) ! [m/s]
      REAL(r8), allocatable :: rsubs_h   (:) ! [m/s]
      REAL(r8), allocatable :: rsubs_fc  (:) ! [m/s]
      REAL(r8) :: rsubs_riv

      logical  :: j_is_river
      REAL(r8) :: theta_s_h, air_h, icefrac, imped, delp
      real(r8) :: sumwt, sumarea, zwt_mean
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
      logical  :: iam_watb, nb_is_watb, has_river

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

            hrus => hillslope_network(ibasin)
            
            theta_a_bsn (ibasin) = 0.
            zwt_bsn     (ibasin) = 0.
            Ks_bsn      (ibasin) = 0.

            nhru = hrus%nhru

            IF (lake_id(ibasin) > 0)               CYCLE  ! lake
            IF ((hrus%indx(1)==0) .and. (nhru==1)) CYCLE  ! only river in a catchment
            
            allocate (theta_a_h (nhru))
            allocate (zwt_h     (nhru))
            allocate (Ks_h      (nhru))

            theta_a_h = 0.
            zwt_h     = 0.
            Ks_h      = 0.

            DO i = 1, nhru
               
               IF (hrus%indx(i) == 0) CYCLE ! river
               IF (hrus%awat(i) == 0) CYCLE ! only land ice or water body in this HRU

               ps = hru_patch%substt(hrus%ihru(i))
               pe = hru_patch%subend(hrus%ihru(i))

               theta_s_h = 0
               DO ipatch = ps, pe
                  IF (patchtype(ipatch) <= 2) THEN
                     theta_s_h = theta_s_h + hru_patch%subfrc(ipatch) &
                        * sum(porsl(1:nl_soil,ipatch) * dz_soi(1:nl_soil) &
                        - wice_soisno(1:nl_soil,ipatch)/denice) / sum(dz_soi(1:nl_soil)) 
                  ENDIF 
               ENDDO

               IF (theta_s_h > 0.) THEN
                  
                  air_h    = 0.
                  zwt_h(i) = 0.
                  DO ipatch = ps, pe
                     IF (patchtype(ipatch) <= 2) THEN
                        air_h = air_h + hru_patch%subfrc(ipatch) &
                           * (sum( porsl(1:nl_soil,ipatch) * dz_soi(1:nl_soil) &
                           - wliq_soisno(1:nl_soil,ipatch)/denh2o &
                           - wice_soisno(1:nl_soil,ipatch)/denice ) - wa(ipatch)/1.0e3)
                        air_h = max(0., air_h)
                        
                        zwt_h(i) = zwt_h(i) + zwt(ipatch) * hru_patch%subfrc(ipatch)
                     ENDIF
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
                     IF (patchtype(ipatch) <= 2) THEN
                        DO ilev = 1, nl_soil
                           icefrac = min(1., wice_soisno(ilev,ipatch)/denice/dz_soi(ilev)/porsl(ilev,ipatch))
                           imped   = 10.**(-e_ice*icefrac)
                           Ks_h(i) = Ks_h(i) + hru_patch%subfrc(ipatch) &
                              * hksati(ilev,ipatch)/1.0e3 * imped * dz_soi(ilev)/zi_soi(nl_soil) 
                        ENDDO
                     ENDIF
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

            DO i = 1, nhru
                  
               j = hrus%inext(i)

               IF (j <= 0)        CYCLE ! downstream is out of catchment
               IF (Ks_h(i) == 0.) CYCLE ! this HRU is frozen
               
               j_is_river = (hrus%indx(j) == 0)

               IF ((.not. j_is_river) .and. (Ks_h(j) == 0.)) CYCLE ! non-river downstream HRU is frozen
                  
               zsubs_h_up = hrus%elva(i) - zwt_h(i)

               IF (.not. j_is_river) THEN
                  zsubs_h_dn = hrus%elva(j) - zwt_h(j)
               ELSE
                  zsubs_h_dn = hrus%elva(1) - riverdpth(ibasin) + wdsrf_hru(hrus%ihru(1)) 
               ENDIF

               IF (.not. j_is_river) THEN
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
               IF ((zsubs_h_up > zsubs_h_dn) .or. j_is_river) THEN
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

               ca = hrus%flen(i) * Ks_fc / theta_a_h(i) / delp / hrus%awat(i) * deltime

               IF (.not. j_is_river) THEN
                  cb = hrus%flen(i) * Ks_fc / theta_a_h(j) / delp / hrus%awat(j) * deltime
               ELSE
                  cb = hrus%flen(i) * Ks_fc / delp / hrus%awat(j) * deltime
               ENDIF
               
               rsubs_fc(i) = (zsubs_h_up - zsubs_h_dn) * hrus%flen(i) * Ks_fc / (1+ca+cb) / delp

               rsubs_h(i) = rsubs_h(i) + rsubs_fc(i) / hrus%awat(i)
               rsubs_h(j) = rsubs_h(j) - rsubs_fc(i) / hrus%awat(j)
               
            ENDDO
            
            IF (hrus%indx(1) == 0) THEN
               ! rsubs_h(1) is positive = out of soil column
               IF (rsubs_h(1)*deltime > wdsrf_bsn(ibasin)*riverarea(ibasin)) THEN 
                  alp = wdsrf_bsn(ibasin)*riverarea(ibasin) / (rsubs_h(1)*deltime)
                  rsubs_h(1) = rsubs_h(1) * alp
                  DO i = 2, nhru
                     IF (hrus%inext(i) == 1) THEN
                        rsubs_h(i) = rsubs_h(i) - (1.0-alp)*rsubs_fc(i)/hrus%awat(i)
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
               
            ! Update total subsurface lateral flow (1): Between hydrological units
            DO i = 1, nhru
               rsubs_hru(hrus%ihru(i)) = rsubs_h(i)
               
               ps = hru_patch%substt(hrus%ihru(i))
               pe = hru_patch%subend(hrus%ihru(i))
               DO ipatch = ps, pe
                  IF (patchtype(ipatch) <= 2) THEN 
                     rsub(ipatch) = rsub(ipatch) + rsubs_h(i) * 1.e3 ! (positive = out of soil column) 
                  ENDIF
               ENDDO
            ENDDO
            
            DO i = 1, nhru
               ! Inside hydrological units
               IF (hrus%indx(i) /= 0) THEN
                  ps = hru_patch%substt(hrus%ihru(i))
                  pe = hru_patch%subend(hrus%ihru(i))
               
                  IF (zwt_h(i) > 1.5) THEN
                     ! from Fan et al., JGR 112(D10125)
                     Ks_in = raniso * Ks_h(i) * bdamp * exp(-(zwt_h(i)-1.5)/bdamp)
                  ELSE
                     Ks_in = raniso * Ks_h(i) * ((1.5-zwt_h(i)) + bdamp)
                  ENDIF

                  sumwt = sum(hru_patch%subfrc(ps:pe), mask = patchtype(ps:pe) <= 2)
                  IF (sumwt > 0) THEN
                     zwt_mean = sum(zwt(ps:pe)*hru_patch%subfrc(ps:pe), mask = patchtype(ps:pe) <= 2) / sumwt

                     DO ipatch = ps, pe
                        IF (patchtype(ipatch) <= 2) THEN
                           rsubs_pch(ipatch) = - Ks_in * (zwt(ipatch) - zwt_mean) *6.0*pi/hrus%awat(i)
                           ! Update total subsurface lateral flow (2): Between patches
                           rsub(ipatch) = rsub(ipatch) + rsubs_pch(ipatch) * 1.e3 ! m/s to mm/s
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO

            sumarea = sum(hrus%awat)
            IF (sumarea > 0) THEN
               theta_a_bsn (ibasin) = sum(theta_a_h * hrus%awat) / sumarea
               zwt_bsn     (ibasin) = sum(zwt_h     * hrus%awat) / sumarea
               Ks_bsn      (ibasin) = sum(Ks_h      * hrus%awat) / sumarea
            ENDIF

            deallocate (theta_a_h)
            deallocate (zwt_h    )
            deallocate (Ks_h     )
            deallocate (rsubs_h  )
            deallocate (rsubs_fc )

         ENDDO

         CALL retrieve_neighbour_data (theta_a_bsn, theta_a_nb)
         CALL retrieve_neighbour_data (zwt_bsn    , zwt_nb    )
         CALL retrieve_neighbour_data (Ks_bsn     , Ks_nb     )
         CALL retrieve_neighbour_data (wdsrf_bsn  , wdsrf_nb  )

         DO ibasin = 1, numbasin
            
            hrus => hillslope_network(ibasin)
               
            iam_watb = .false.
            IF (lake_id(ibasin) > 0) THEN
               iam_watb = .true.
            ELSEIF ((hrus%nhru == 1) .and. (hrus%indx(1) == 0)) THEN
               iam_watb = .true.
            ENDIF
            
            DO jnb = 1, basinneighbour(ibasin)%nnb

               IF (basinneighbour(ibasin)%bindex(jnb) == -9) CYCLE ! skip ocean neighbour

               nb_is_watb = basinneighbour(ibasin)%iswatb(jnb)

               IF (iam_watb .and. nb_is_watb) then
                  CYCLE
               ENDIF

               IF (.not. iam_watb)   Ks_up = Ks_bsn(ibasin)
               IF (.not. nb_is_watb) Ks_dn = Ks_nb(ibasin)%val(jnb)

               IF ((Ks_up == 0.) .or. (Ks_dn == 0.)) THEN
                  cycle
               ENDIF

               IF (.not. iam_watb)   zwt_up = zwt_bsn(ibasin)
               IF (.not. nb_is_watb) zwt_dn = zwt_nb(ibasin)%val(jnb)

               IF (.not. iam_watb) then
                  theta_a_up = theta_a_bsn(ibasin)
               ELSE
                  theta_a_up = 1.
               ENDIF 

               IF (.not. nb_is_watb) THEN
                  theta_a_dn = theta_a_nb(ibasin)%val(jnb)
               ELSE
                  theta_a_dn = 1.
               ENDIF

               IF (iam_watb) then
                  zsubs_up = basinneighbour(ibasin)%myelva + wdsrf_bsn(ibasin) 
                  IF ((zsubs_up > zsubs_dn) .and. (wdsrf_bsn(ibasin) == 0.)) THEN
                     CYCLE
                  ENDIF
                  delp = basinneighbour(ibasin)%area(jnb) / basinneighbour(ibasin)%lenbdr(jnb) * 0.5
               ELSE
                  zsubs_up = basinneighbour(ibasin)%myelva - zwt_up 
                  delp = basinneighbour(ibasin)%dist(jnb)
               ENDIF

               IF (nb_is_watb) THEN
                  zsubs_dn = basinneighbour(ibasin)%elva(jnb) + wdsrf_nb(ibasin)%val(jnb)
                  IF ((zsubs_up < zsubs_dn) .and. (wdsrf_nb(ibasin)%val(jnb) == 0.)) THEN
                     CYCLE
                  ENDIF
                  delp = basinneighbour(ibasin)%myarea / basinneighbour(ibasin)%lenbdr(jnb) * 0.5
               ELSE
                  zsubs_dn = basinneighbour(ibasin)%elva(jnb) - zwt_dn
                  delp = basinneighbour(ibasin)%dist(jnb)
               ENDIF 
                  
               area_up = basinneighbour(ibasin)%myarea
               area_dn = basinneighbour(ibasin)%area(jnb)

               lenbdr = basinneighbour(ibasin)%lenbdr(jnb)

               ! from Fan et al., JGR 112(D10125)
               slope = abs(basinneighbour(ibasin)%slope(jnb))
               IF (slope > 0.16) THEN
                  bdamp = 4.8
               ELSE
                  bdamp = 120./(1+150.*slope)
               ENDIF

               ! Upstream scheme for hydraulic conductivity
               IF (nb_is_watb .or. ((.not. iam_watb) .and. (zsubs_up > zsubs_dn))) THEN
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

               IF (.not. iam_watb) THEN
                  IF (hrus%indx(1) == 0) THEN
                     rsubs_nb = rsubs_nb / sum(hrus%area(2:))
                  ELSE
                     rsubs_nb = rsubs_nb / basinneighbour(ibasin)%myarea
                  ENDIF
               ELSE
                  rsubs_nb = rsubs_nb / basinneighbour(ibasin)%myarea
               ENDIF
               
               rsubs_bsn(ibasin) = rsubs_bsn(ibasin) + rsubs_nb

            ENDDO

            IF (iam_watb) THEN
               ps = elm_patch%substt(ibasin)
               pe = elm_patch%subend(ibasin)
               ! Update total subsurface lateral flow (3): Between basins
               rsub(ps:pe) = rsub(ps:pe) + rsubs_bsn(ibasin) * 1.e3 ! m/s to mm/s
            ELSE
               IF (hrus%indx(1) == 0) THEN
                  i0 = 2 ! excluding river HRU
               ELSE
                  i0 = 1
               ENDIF

               DO i = i0, hrus%nhru
                  ps = hru_patch%substt(hrus%ihru(i))
                  pe = hru_patch%subend(hrus%ihru(i))
                  ! Update total subsurface lateral flow (3): Between basins
                  rsub(ps:pe) = rsub(ps:pe) + rsubs_bsn(ibasin) * 1.e3 ! m/s to mm/s
               ENDDO
            ENDIF

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

            IF (patchtype(ipatch) <= 2) THEN 
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
                  ELSE
                     vol_liq(ilev) = eff_porosity(ilev)
                     wresi(ilev) = 0.
                  ENDIF
               ENDDO
      
               zwtmm = zwt(ipatch) * 1000. ! m -> mm

               ! check consistancy between water table location and liquid water content
               DO ilev = 1, nl_soil
                  IF ((vol_liq(ilev) < eff_porosity(ilev)-1.e-6) .and. (zwtmm <= sp_zi(ilev-1))) THEN
                     zwtmm = sp_zi(ilev)
                  ENDIF
               ENDDO

               izwt = findloc(zwtmm >= sp_zi, .true., dim=1, back=.true.)

               IF (izwt <= nl_soil) THEN
                  IF (is_permeable(izwt) .and. (zwtmm > sp_zi(izwt-1))) THEN
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
                  nl_soil, exwater, sp_zi, is_permeable, eff_porosity, vl_r, psi0(:,ipatch), &
                  hksati(:,ipatch), nprms, prms, porsl(nl_soil,ipatch), wdsrf(ipatch), &
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
                  write(6,'(A,I0,4E20.5)') 'Warning (Subsurface Runoff): water balance violation ', &
                     ipatch, errblc, exwater, zwtmm
                  STOP
               endif
#endif
            ELSEIF (patchtype(ipatch) == 4) THEN ! land water bodies
               IF (wa(ipatch) < 0) THEN
                  wa(ipatch) = wa(ipatch) - rsub(ipatch)*deltime
                  IF (wa(ipatch) > 0) THEN
                     wdsrf(ipatch) = wa(ipatch) 
                     wa(ipatch) = 0
                  ENDIF
               ELSE
                  wdsrf(ipatch) = wdsrf(ipatch) - rsub(ipatch)*deltime
                  IF (wdsrf(ipatch) < 0) THEN
                     wa(ipatch) = wdsrf(ipatch) 
                     wdsrf(ipatch) = 0
                  ENDIF
               ENDIF
            ENDIF

         ENDDO
      ENDIF

   END SUBROUTINE subsurface_flow

END MODULE MOD_Hydro_SubsurfaceFlow
#endif
