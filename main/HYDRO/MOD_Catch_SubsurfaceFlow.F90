#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_SubsurfaceFlow
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Ground water lateral flow.
!
!   Ground water fluxes are calculated
!   1. between elements
!   2. between hydrological response units
!   3. between patches inside one HRU
!
! Created by Shupeng Zhang, May 2023
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_DataType
   USE MOD_Catch_HillslopeNetwork
   IMPLICIT NONE

   ! --- information of HRU on hillslope ---
   type(hillslope_network_type), pointer :: hillslope_element (:)

   integer,  allocatable :: lake_id_elm  (:)
   real(r8), allocatable :: lakedepth_elm(:)
   real(r8), allocatable :: riverdpth_elm(:)
   real(r8), allocatable :: wdsrf_elm    (:)

   real(r8), parameter :: e_ice  = 6.0   ! soil ice impedance factor

   ! anisotropy ratio of lateral/vertical hydraulic conductivity (unitless)
   ! for USDA soil texture class:
   ! 0: undefined
   ! 1: clay;  2: silty clay;  3: sandy clay;   4: clay loam;   5: silty clay loam;   6: sandy clay loam; &
   ! 7: loam;  8: silty loam;  9: sandy loam;  10: silt;       11: loamy sand;       12: sand
   real(r8), parameter :: raniso(0:12) = (/ 1., &
                                            48., 40., 28., 24., 20., 14., 12., 10., 4., 2., 3., 2. /)

   ! -- neighbour variables --
   type(pointer_real8_1d), allocatable :: agwt_nb    (:)  ! ground water area (for patchtype <= 2) of neighbours [m^2]
   type(pointer_real8_1d), allocatable :: theta_a_nb (:)  ! saturated volume content [-]
   type(pointer_real8_1d), allocatable :: zwt_nb     (:)  ! water table depth [m]
   type(pointer_real8_1d), allocatable :: Kl_nb      (:)  ! lateral hydraulic conductivity [m/s]
   type(pointer_real8_1d), allocatable :: wdsrf_nb   (:)  ! depth of surface water [m]
   type(pointer_logic_1d), allocatable :: islake_nb  (:)  ! whether a neighbour is water body
   type(pointer_real8_1d), allocatable :: lakedp_nb  (:)  ! lake depth of neighbour [m]

CONTAINS

   ! ----------
   SUBROUTINE subsurface_network_init (patcharea)

   USE MOD_SPMD_Task
   USE MOD_Utils
   USE MOD_Mesh
   USE MOD_Pixel
   USE MOD_LandElm
   USE MOD_LandPatch
   USE MOD_ElementNeighbour
   USE MOD_WorkerPushData,         only: worker_push_data
   USE MOD_Catch_BasinNetwork,     only: push_bsn2elm
   USE MOD_Catch_RiverLakeNetwork, only: lake_id, riverdpth
   USE MOD_Vars_TimeInvariants,    only: patchtype, lakedepth
   IMPLICIT NONE

   real(r8), intent(in) :: patcharea (:)

   integer :: ielm, inb, i, ihru, ps, pe, ipatch

   real(r8), allocatable :: agwt_b(:)
   real(r8), allocatable :: islake(:)
   type(pointer_real8_1d), allocatable :: iswat_nb (:)

   integer, allocatable :: eindex(:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_worker) THEN
         IF (numelm > 0) THEN
            allocate (eindex (numelm))
            eindex = landelm%eindex
         ENDIF
      ENDIF

      CALL hillslope_network_init (numelm, eindex, hillslope_element)

      IF (allocated(eindex)) deallocate (eindex)

      IF (p_is_worker) THEN

         IF (numelm > 0) allocate (lake_id_elm  (numelm))
         IF (numelm > 0) allocate (riverdpth_elm(numelm))
         IF (numelm > 0) allocate (lakedepth_elm(numelm))
         IF (numelm > 0) allocate (wdsrf_elm    (numelm))

         CALL worker_push_data (push_bsn2elm, lake_id,   lake_id_elm,   -9999)
         CALL worker_push_data (push_bsn2elm, riverdpth, riverdpth_elm, spval)

         DO ielm = 1, numelm
            IF (lake_id_elm(ielm) <= 0) THEN
               DO i = 1, hillslope_element(ielm)%nhru

                  hillslope_element(ielm)%agwt(i) = 0
                  hillslope_element(ielm)%area(i) = 0

                  ihru = hillslope_element(ielm)%ihru(i)
                  ps = hru_patch%substt(ihru)
                  pe = hru_patch%subend(ihru)
                  DO ipatch = ps, pe
                     hillslope_element(ielm)%area(i) = hillslope_element(ielm)%area(i) + patcharea(ipatch)
                     IF (patchtype(ipatch) <= 2) THEN
                        hillslope_element(ielm)%agwt(i) = hillslope_element(ielm)%agwt(i) + patcharea(ipatch)
                     ENDIF
                  ENDDO

               ENDDO
            ENDIF
         ENDDO

         lakedepth_elm(:) = 0.
         DO ielm = 1, numelm
            IF (lake_id_elm(ielm) > 0) THEN
               ps = elm_patch%substt(ielm)
               pe = elm_patch%subend(ielm)
               lakedepth_elm(ielm) = sum(lakedepth(ps:pe) * elm_patch%subfrc(ps:pe))
            ENDIF
         ENDDO

         CALL allocate_neighbour_data (agwt_nb   )
         CALL allocate_neighbour_data (theta_a_nb)
         CALL allocate_neighbour_data (zwt_nb    )
         CALL allocate_neighbour_data (Kl_nb     )
         CALL allocate_neighbour_data (wdsrf_nb  )
         CALL allocate_neighbour_data (islake_nb )
         CALL allocate_neighbour_data (lakedp_nb )
         CALL allocate_neighbour_data (iswat_nb  )

         IF (numelm > 0) THEN
            allocate (agwt_b(numelm))
            allocate (islake(numelm))
            DO ielm = 1, numelm
               IF (lake_id_elm(ielm) <= 0) THEN
                  agwt_b(ielm) = sum(hillslope_element(ielm)%agwt)
                  islake(ielm) = 0.
               ELSE
                  agwt_b(ielm) = 0.
                  islake(ielm) = 1.
               ENDIF
            ENDDO
         ENDIF

         CALL retrieve_neighbour_data (lakedepth_elm, lakedp_nb)

         CALL retrieve_neighbour_data (agwt_b, agwt_nb )
         CALL retrieve_neighbour_data (islake, iswat_nb)

         DO ielm = 1, numelm
            DO inb = 1, elementneighbour(ielm)%nnb
               IF (elementneighbour(ielm)%glbindex(inb) > 0) THEN ! skip ocean neighbour
                  islake_nb(ielm)%val(inb) = (iswat_nb(ielm)%val(inb) > 0)
               ENDIF
            ENDDO
         ENDDO

         IF (allocated(agwt_b  )) deallocate(agwt_b  )
         IF (allocated(islake  )) deallocate(islake  )
         IF (allocated(iswat_nb)) deallocate(iswat_nb)

      ENDIF

   END SUBROUTINE subsurface_network_init

   ! ---------
   SUBROUTINE subsurface_flow (deltime)

   USE MOD_SPMD_Task
   USE MOD_UserDefFun
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_LandHRU
   USE MOD_LandPatch
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_1DFluxes
   USE MOD_Catch_HillslopeNetwork
   USE MOD_ElementNeighbour
   USE MOD_Const_Physical,  only: denice, denh2o
   USE MOD_Vars_Global,     only: pi, nl_soil, zi_soi
   USE MOD_Hydro_SoilWater, only: soilwater_aquifer_exchange

   IMPLICIT NONE

   real(r8), intent(in) :: deltime

   ! Local Variables
   integer :: nhru, ielm, i, i0, j, ihru, ipatch, ps, pe, hs, he, ilev

   type(hillslope_network_type), pointer :: hrus

   real(r8), allocatable :: theta_a_h (:)
   real(r8), allocatable :: zwt_h     (:)
   real(r8), allocatable :: Kl_h      (:) ! [m/s]
   real(r8), allocatable :: xsubs_h   (:) ! [m/s]
   real(r8), allocatable :: xsubs_fc  (:) ! [m/s]

   logical  :: j_is_river
   real(r8) :: theta_s_h, air_h, icefrac, imped, delp
   real(r8) :: sumwt, sumarea, zwt_mean
   real(r8) :: zsubs_h_up, zsubs_h_dn
   real(r8) :: slope, bdamp, Kl_fc, Kl_in
   real(r8) :: ca, cb
   real(r8) :: alp

   real(r8), allocatable :: theta_a_elm (:)
   real(r8), allocatable :: zwt_elm     (:)
   real(r8), allocatable :: Kl_elm      (:) ! [m/s]

   integer  :: jnb
   real(r8) :: zsubs_up, zwt_up, Kl_up, theta_a_up, area_up
   real(r8) :: zsubs_dn, zwt_dn, Kl_dn, theta_a_dn, area_dn
   real(r8) :: lenbdr, xsubs_nb
   logical  :: iam_lake, nb_is_lake, has_river

   ! for water exchange
   logical  :: is_dry_lake
   integer  :: izwt
   real(r8) :: exwater
   real(r8) :: sp_zi(0:nl_soil), sp_dz(1:nl_soil), zwtmm ! [mm]
   real(r8) :: vl_r (1:nl_soil)
#ifdef Campbell_SOIL_MODEL
   integer, parameter :: nprms = 1
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   integer, parameter :: nprms = 5
#endif
   real(r8) :: prms  (nprms,1:nl_soil)
   real(r8) :: vol_ice     (1:nl_soil)
   real(r8) :: vol_liq     (1:nl_soil)
   real(r8) :: eff_porosity(1:nl_soil)
   logical  :: is_permeable(1:nl_soil)
   real(r8) :: wresi       (1:nl_soil)
   real(r8) :: w_sum_before, w_sum_after, errblc


      IF (p_is_worker) THEN

         xsubs_elm(:) = 0.  ! subsurface lateral flow between element basins
         xsubs_hru(:) = 0.  ! subsurface lateral flow between hydrological response units
         xsubs_pch(:) = 0.  ! subsurface lateral flow between patches inside one HRU

         xwsub(:) = 0. ! total recharge/discharge from subsurface lateral flow

         IF (numpatch > 0) rsub(:) = 0.

         IF (numelm > 0) THEN
            allocate (theta_a_elm (numelm));  theta_a_elm = 0.
            allocate (zwt_elm     (numelm));  zwt_elm     = 0.
            allocate (Kl_elm      (numelm));  Kl_elm      = 0.
         ENDIF

         DO ielm = 1, numelm

            hrus => hillslope_element(ielm)

            nhru = hrus%nhru

            IF (lake_id_elm(ielm) > 0) CYCLE  ! lake
            IF (sum(hrus%agwt) <= 0)   CYCLE  ! no area of soil, urban or wetland

            allocate (theta_a_h (nhru));  theta_a_h = 0.
            allocate (zwt_h     (nhru));  zwt_h     = 0.
            allocate (Kl_h      (nhru));  Kl_h      = 0.

            DO i = 1, nhru

               IF (hrus%indx(i) == 0) CYCLE ! river
               IF (hrus%agwt(i) == 0) CYCLE ! no area of soil, urban or wetland

               ps = hru_patch%substt(hrus%ihru(i))
               pe = hru_patch%subend(hrus%ihru(i))

               theta_s_h = 0
               sumwt = 0
               DO ipatch = ps, pe
                  IF (patchtype(ipatch) <= 2) THEN
                     theta_s_h = theta_s_h + hru_patch%subfrc(ipatch) &
                        * sum(porsl(1:nl_soil,ipatch) * dz_soi(1:nl_soil) &
                        - wice_soisno(1:nl_soil,ipatch)/denice) / sum(dz_soi(1:nl_soil))
                     sumwt = sumwt + hru_patch%subfrc(ipatch)
                  ENDIF
               ENDDO
               IF (sumwt > 0) theta_s_h = theta_s_h / sumwt

               IF (theta_s_h > 0.) THEN

                  air_h    = 0.
                  zwt_h(i) = 0.
                  sumwt    = 0.
                  DO ipatch = ps, pe
                     IF (patchtype(ipatch) <= 2) THEN
                        air_h = air_h + hru_patch%subfrc(ipatch) &
                           * (sum( porsl(1:nl_soil,ipatch) * dz_soi(1:nl_soil) &
                           - wliq_soisno(1:nl_soil,ipatch)/denh2o &
                           - wice_soisno(1:nl_soil,ipatch)/denice ) - wa(ipatch)/1.0e3)
                        air_h = max(0., air_h)

                        zwt_h(i) = zwt_h(i) + zwt(ipatch) * hru_patch%subfrc(ipatch)

                        sumwt = sumwt + hru_patch%subfrc(ipatch)
                     ENDIF
                  ENDDO
                  IF (sumwt > 0) air_h    = air_h / sumwt
                  IF (sumwt > 0) zwt_h(i) = zwt_h(i) / sumwt

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

                  Kl_h(i) = 0.
                  sumwt   = 0.
                  DO ipatch = ps, pe
                     IF (patchtype(ipatch) <= 2) THEN
                        DO ilev = 1, nl_soil
                           icefrac = min(1., wice_soisno(ilev,ipatch)/denice/dz_soi(ilev)/porsl(ilev,ipatch))
                           imped   = 10.**(-e_ice*icefrac)
                           Kl_h(i) = Kl_h(i) + hru_patch%subfrc(ipatch) * raniso(soiltext(ipatch)) &
                              * hksati(ilev,ipatch)/1.0e3 * imped * dz_soi(ilev)/zi_soi(nl_soil)
                        ENDDO
                        sumwt = sumwt + hru_patch%subfrc(ipatch)
                     ENDIF
                  ENDDO
                  IF (sumwt > 0) Kl_h(i) = Kl_h(i) / sumwt
               ELSE
                  ! Frozen soil.
                  Kl_h(i) = 0.
               ENDIF

            ENDDO

            allocate (xsubs_h  (nhru))
            allocate (xsubs_fc (nhru))

            xsubs_h (:) = 0.
            xsubs_fc(:) = 0.

            DO i = 1, nhru

               j = hrus%inext(i)

               IF (j <= 0)        CYCLE ! downstream is out of catchment
               IF (Kl_h(i) == 0.) CYCLE ! this HRU is frozen

               j_is_river = (hrus%indx(j) == 0)

               IF ((.not. j_is_river) .and. (Kl_h(j) == 0.)) CYCLE ! non-river downstream HRU is frozen

               zsubs_h_up = hrus%elva(i) - zwt_h(i)

               IF (.not. j_is_river) THEN
                  zsubs_h_dn = hrus%elva(j) - zwt_h(j)
               ELSE
                  zsubs_h_dn = hrus%elva(1) - riverdpth_elm(ielm) + wdsrf_hru(hrus%ihru(1))
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
                     Kl_fc = Kl_h(i) * bdamp * exp(-(zwt_h(i)-1.5)/bdamp)
                  ELSE
                     Kl_fc = Kl_h(i) * ((1.5-zwt_h(i)) + bdamp)
                  ENDIF
               ELSE
                  IF (zwt_h(j) > 1.5) THEN
                     Kl_fc = Kl_h(j) * bdamp * exp(-(zwt_h(j)-1.5)/bdamp)
                  ELSE
                     Kl_fc = Kl_h(j) * ((1.5-zwt_h(j)) + bdamp)
                  ENDIF
               ENDIF

               ca = hrus%flen(i) * Kl_fc / theta_a_h(i) / delp / hrus%agwt(i) * deltime

               IF (.not. j_is_river) THEN
                  cb = hrus%flen(i) * Kl_fc / theta_a_h(j) / delp / hrus%agwt(j) * deltime
               ELSE
                  cb = hrus%flen(i) * Kl_fc / delp / hrus%area(j) * deltime
               ENDIF

               xsubs_fc(i) = (zsubs_h_up - zsubs_h_dn) * hrus%flen(i) * Kl_fc / (1+ca+cb) / delp

               xsubs_h(i) = xsubs_h(i) + xsubs_fc(i) / hrus%agwt(i)

               IF (j_is_river) THEN
                  xsubs_h(j) = xsubs_h(j) - xsubs_fc(i) / hrus%area(j)
               ELSE
                  xsubs_h(j) = xsubs_h(j) - xsubs_fc(i) / hrus%agwt(j)
               ENDIF

            ENDDO

            IF (hrus%indx(1) == 0) THEN
               ! xsubs_h(1) is positive = out of soil column
               IF (xsubs_h(1)*deltime > wdsrf_hru(hrus%ihru(1))) THEN
                  alp = wdsrf_hru(hrus%ihru(1)) / (xsubs_h(1)*deltime)
                  xsubs_h(1) = xsubs_h(1) * alp
                  DO i = 2, nhru
                     IF ((hrus%inext(i) == 1) .and. (hrus%agwt(i) > 0.)) THEN
                        xsubs_h(i) = xsubs_h(i) - (1.0-alp)*xsubs_fc(i)/hrus%agwt(i)
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            ! Update total subsurface lateral flow (1): Between hydrological units
            ! for soil, urban, wetland or river patches
            DO i = 1, nhru
               xsubs_hru(hrus%ihru(i)) = xsubs_h(i)

               ps = hru_patch%substt(hrus%ihru(i))
               pe = hru_patch%subend(hrus%ihru(i))
               DO ipatch = ps, pe
                  IF ((patchtype(ipatch) <= 2) .or. (hrus%indx(i) == 0)) THEN
                     xwsub(ipatch) = xwsub(ipatch) + xsubs_h(i) * 1.e3 ! (positive = out of soil column)
                  ENDIF
               ENDDO

               IF (hrus%indx(1) == 0) THEN
                  DO ipatch = ps, pe
                     IF (patchtype(ipatch) <= 2) THEN
                        rsub(ipatch) = - xsubs_h(1) * hrus%area(1) / sum(hrus%agwt) * 1.0e3 ! m/s to mm/s
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO

            DO i = 1, nhru
               ! Inside hydrological units
               IF (hrus%agwt(i) > 0) THEN

                  bdamp = 4.8

                  IF (zwt_h(i) > 1.5) THEN
                     ! from Fan et al., JGR 112(D10125)
                     Kl_in = Kl_h(i) * bdamp * exp(-(zwt_h(i)-1.5)/bdamp)
                  ELSE
                     Kl_in = Kl_h(i) * ((1.5-zwt_h(i)) + bdamp)
                  ENDIF

                  ps = hru_patch%substt(hrus%ihru(i))
                  pe = hru_patch%subend(hrus%ihru(i))
                  sumwt = sum(hru_patch%subfrc(ps:pe), mask = patchtype(ps:pe) <= 2)
                  IF (sumwt > 0) THEN
                     zwt_mean = sum(zwt(ps:pe)*hru_patch%subfrc(ps:pe), mask = patchtype(ps:pe) <= 2) / sumwt

                     DO ipatch = ps, pe
                        IF (patchtype(ipatch) <= 2) THEN
                           xsubs_pch(ipatch) = - Kl_in * (zwt(ipatch) - zwt_mean) *6.0*pi/hrus%agwt(i)
                           ! Update total subsurface lateral flow (2): Between patches
                           xwsub(ipatch) = xwsub(ipatch) + xsubs_pch(ipatch) * 1.e3 ! m/s to mm/s
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO

            sumarea = sum(hrus%agwt)
            IF (sumarea > 0) THEN
               theta_a_elm (ielm) = sum(theta_a_h * hrus%agwt) / sumarea
               zwt_elm     (ielm) = sum(zwt_h     * hrus%agwt) / sumarea
               Kl_elm      (ielm) = sum(Kl_h      * hrus%agwt) / sumarea
            ENDIF

            deallocate (theta_a_h)
            deallocate (zwt_h    )
            deallocate (Kl_h     )
            deallocate (xsubs_h  )
            deallocate (xsubs_fc )

         ENDDO

         DO ielm = 1, numelm
            hs = elm_hru%substt(ielm)
            he = elm_hru%subend(ielm)
            wdsrf_elm(ielm) = sum(wdsrf_hru(hs:he) * elm_hru%subfrc(hs:he))
         ENDDO

         CALL retrieve_neighbour_data (theta_a_elm, theta_a_nb)
         CALL retrieve_neighbour_data (zwt_elm    , zwt_nb    )
         CALL retrieve_neighbour_data (Kl_elm     , Kl_nb     )
         CALL retrieve_neighbour_data (wdsrf_elm  , wdsrf_nb  )

         DO ielm = 1, numelm

            hrus => hillslope_element(ielm)

            iam_lake = (lake_id_elm(ielm) > 0)

            DO jnb = 1, elementneighbour(ielm)%nnb

               IF (elementneighbour(ielm)%glbindex(jnb) == -9) CYCLE ! skip ocean neighbour

               nb_is_lake = islake_nb(ielm)%val(jnb)

               IF (iam_lake .and. nb_is_lake) THEN
                  CYCLE
               ENDIF

               IF (.not. iam_lake) THEN
                  Kl_up      = Kl_elm    (ielm)
                  zwt_up     = zwt_elm    (ielm)
                  theta_a_up = theta_a_elm(ielm)
                  zsubs_up   = elementneighbour(ielm)%myelva - zwt_up
                  area_up    = sum(hrus%agwt)
               ELSE
                  theta_a_up = 1.
                  zsubs_up   = elementneighbour(ielm)%myelva - lakedepth_elm(ielm) + wdsrf_elm(ielm)
                  area_up    = elementneighbour(ielm)%myarea
               ENDIF

               IF (.not. nb_is_lake) THEN
                  Kl_dn      = Kl_nb(ielm)%val(jnb)
                  zwt_dn     = zwt_nb(ielm)%val(jnb)
                  theta_a_dn = theta_a_nb(ielm)%val(jnb)
                  zsubs_dn   = elementneighbour(ielm)%elva(jnb) - zwt_dn
                  area_dn    = agwt_nb(ielm)%val(jnb)
               ELSE
                  theta_a_dn = 1.
                  zsubs_dn   = elementneighbour(ielm)%elva(jnb) - lakedp_nb(ielm)%val(jnb) + wdsrf_nb(ielm)%val(jnb)
                  area_dn    = elementneighbour(ielm)%area(jnb)
               ENDIF

               IF ((.not. iam_lake)   .and. (area_up <= 0)) CYCLE
               IF ((.not. nb_is_lake) .and. (area_dn <= 0)) CYCLE
               IF ((.not. iam_lake)   .and. (Kl_up == 0. )) CYCLE
               IF ((.not. nb_is_lake) .and. (Kl_dn == 0. )) CYCLE

               ! water body is dry.
               IF (iam_lake .and. (zsubs_up > zsubs_dn) .and. (wdsrf_elm(ielm) == 0.)) THEN
                  CYCLE
               ENDIF
               IF (nb_is_lake .and. (zsubs_up < zsubs_dn) .and. (wdsrf_nb(ielm)%val(jnb) == 0.)) THEN
                  CYCLE
               ENDIF

               lenbdr = elementneighbour(ielm)%lenbdr(jnb)

               delp = elementneighbour(ielm)%dist(jnb)
               IF (iam_lake) THEN
                  delp = elementneighbour(ielm)%area(jnb) / lenbdr * 0.5
               ENDIF
               IF (nb_is_lake) THEN
                  delp = elementneighbour(ielm)%myarea / lenbdr * 0.5
               ENDIF

               ! from Fan et al., JGR 112(D10125)
               slope = abs(elementneighbour(ielm)%slope(jnb))
               IF (slope > 0.16) THEN
                  bdamp = 4.8
               ELSE
                  bdamp = 120./(1+150.*slope)
               ENDIF

               ! Upstream scheme for hydraulic conductivity
               IF (nb_is_lake .or. ((.not. iam_lake) .and. (zsubs_up > zsubs_dn))) THEN
                  IF (zwt_up > 1.5) THEN
                     ! from Fan et al., JGR 112(D10125)
                     Kl_fc = Kl_up * bdamp * exp(-(zwt_up-1.5)/bdamp)
                  ELSE
                     Kl_fc = Kl_up * ((1.5-zwt_up) + bdamp)
                  ENDIF
               ELSE
                  IF (zwt_dn > 1.5) THEN
                     Kl_fc = Kl_dn * bdamp * exp(-(zwt_dn-1.5)/bdamp)
                  ELSE
                     Kl_fc = Kl_dn * ((1.5-zwt_dn) + bdamp)
                  ENDIF
               ENDIF

               ca = lenbdr * Kl_fc / theta_a_up / delp / area_up * deltime
               cb = lenbdr * Kl_fc / theta_a_dn / delp / area_dn * deltime

               xsubs_nb = (zsubs_up - zsubs_dn) * lenbdr * Kl_fc / (1+ca+cb) / delp

               IF (.not. iam_lake) THEN
                  xsubs_nb = xsubs_nb / sum(hrus%agwt)
               ELSE
                  xsubs_nb = xsubs_nb / elementneighbour(ielm)%myarea
               ENDIF

               xsubs_elm(ielm) = xsubs_elm(ielm) + xsubs_nb

               IF (nb_is_lake) THEN
                  ps = elm_patch%substt(ielm)
                  pe = elm_patch%subend(ielm)
                  DO ipatch = ps, pe
                     IF (patchtype(ipatch) <= 2) THEN
                        rsub(ipatch) = rsub(ipatch) + xsubs_nb * 1.e3
                     ENDIF
                  ENDDO
               ENDIF

            ENDDO

            ! Update total subsurface lateral flow (3): Between basins
            ps = elm_patch%substt(ielm)
            pe = elm_patch%subend(ielm)
            DO ipatch = ps, pe
               IF (iam_lake .or. (patchtype(ipatch) <= 2)) THEN
                  xwsub(ipatch) = xwsub(ipatch) + xsubs_elm(ielm) * 1.e3 ! m/s to mm/s
               ENDIF
            ENDDO

         ENDDO

         IF (allocated(theta_a_elm)) deallocate(theta_a_elm)
         IF (allocated(zwt_elm    )) deallocate(zwt_elm    )
         IF (allocated(Kl_elm     )) deallocate(Kl_elm     )

      ENDIF

      ! Exchange between soil water and aquifer.
      IF (p_is_worker) THEN

         sp_zi(0) = 0.
         sp_zi(1:nl_soil) = zi_soi(1:nl_soil) * 1000.0   ! from meter to mm
         sp_dz(1:nl_soil) = sp_zi(1:nl_soil) - sp_zi(0:nl_soil-1)

         DO ipatch = 1, numpatch

#if (defined CoLMDEBUG)
            ! For water balance check, the sum of water in soil column before the calcultion
            w_sum_before = sum(wliq_soisno(1:nl_soil,ipatch)) + sum(wice_soisno(1:nl_soil,ipatch)) &
               + wa(ipatch) + wdsrf(ipatch) + wetwat(ipatch)
#endif

            IF (DEF_USE_Dynamic_Lake) THEN
               is_dry_lake = (patchtype(ipatch) == 4) .and. (zwt(ipatch) > 0.)
            ELSE
               is_dry_lake = .false.
            ENDIF

            IF ((patchtype(ipatch) <= 1) .or. is_dry_lake) THEN

               exwater = xwsub(ipatch) * deltime

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
                  IF ((vol_liq(ilev) < eff_porosity(ilev)-1.e-8) .and. (zwtmm <= sp_zi(ilev-1))) THEN
                     zwtmm = sp_zi(ilev)
                  ENDIF
               ENDDO

               izwt = findloc_ud(zwtmm >= sp_zi, back=.true.)

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

            ELSEIF (patchtype(ipatch) == 2) THEN ! wetland

               wetwat(ipatch) = wdsrf(ipatch) + wa(ipatch) + wetwat(ipatch)  - xwsub(ipatch)*deltime

               IF (wetwat(ipatch) > wetwatmax) THEN
                  wdsrf (ipatch) = wetwat(ipatch) - wetwatmax
                  wetwat(ipatch) = wetwatmax
                  wa    (ipatch) = 0.
               ELSEIF (wetwat(ipatch) < 0) THEN
                  wa    (ipatch) = wetwat(ipatch)
                  wdsrf (ipatch) = 0.
                  wetwat(ipatch) = 0.
               ELSE
                  wdsrf(ipatch)  = 0.
                  wa   (ipatch)  = 0.
               ENDIF

            ELSEIF (patchtype(ipatch) == 4) THEN ! land water bodies

               wdsrf(ipatch) = wa(ipatch) + wdsrf(ipatch) - xwsub(ipatch)*deltime

               IF (wdsrf(ipatch) < 0) THEN
                  wa   (ipatch) = wdsrf(ipatch)
                  wdsrf(ipatch) = 0
               ELSE
                  wa(ipatch) = 0
               ENDIF

            ENDIF

#if (defined CoLMDEBUG)
            ! For water balance check, the sum of water in soil column after the calcultion
            w_sum_after = sum(wliq_soisno(1:nl_soil,ipatch)) + sum(wice_soisno(1:nl_soil,ipatch)) &
               + wa(ipatch) + wdsrf(ipatch) + wetwat(ipatch)
            errblc = w_sum_after - w_sum_before + xwsub(ipatch)*deltime

            IF(abs(errblc) > 1.e-3)THEN
               write(6,'(A,I0,4E20.5)') 'Warning (Subsurface Runoff): water balance violation ', &
                  ipatch, errblc, xwsub(ipatch), zwtmm
               write(*,*) patchtype(ipatch)
               CALL CoLM_stop ()
            ENDIF
#endif
         ENDDO
      ENDIF

   END SUBROUTINE subsurface_flow

   ! ----------
   SUBROUTINE subsurface_network_final ()

   IMPLICIT NONE

      IF (allocated(lake_id_elm  )) deallocate(lake_id_elm  )
      IF (allocated(riverdpth_elm)) deallocate(riverdpth_elm)
      IF (allocated(lakedepth_elm)) deallocate(lakedepth_elm)
      IF (allocated(wdsrf_elm    )) deallocate(wdsrf_elm    )

      IF (allocated(theta_a_nb)) deallocate(theta_a_nb)
      IF (allocated(zwt_nb    )) deallocate(zwt_nb    )
      IF (allocated(Kl_nb     )) deallocate(Kl_nb     )
      IF (allocated(wdsrf_nb  )) deallocate(wdsrf_nb  )
      IF (allocated(agwt_nb   )) deallocate(agwt_nb   )
      IF (allocated(islake_nb )) deallocate(islake_nb )
      IF (allocated(lakedp_nb )) deallocate(lakedp_nb )

      IF (associated(hillslope_element)) deallocate(hillslope_element)

   END SUBROUTINE subsurface_network_final

END MODULE MOD_Catch_SubsurfaceFlow
#endif
