#include <define.h>

MODULE MOD_Aerosol

  !-----------------------------------------------------------------------
  use MOD_Precision
  IMPLICIT NONE
  SAVE
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: AerosolMasses
  public :: AerosolFluxes
  !
  ! !PUBLIC DATA MEMBERS:
  !-----------------------------------------------------------------------
    !integer, parameter :: begc = 1      !  beginning column index
    !integer, parameter :: endc = 1      !  beginning and ending column index
    real(r8) :: dtime = 1800.0_r8       !  land model time step (sec)

    logical, parameter :: use_extrasnowlayers = .false.
    integer, parameter :: nlevsno = 5   !  maximum number of snow layers

    real(r8), public, parameter :: snw_rds_min = 54.526_r8  ! minimum allowed snow effective radius (also "fresh snow" value) [microns

contains

  !-----------------------------------------------------------------------
  subroutine AerosolMasses( snl           ,do_capsnow     ,&
             h2osno_ice    ,h2osno_liq    ,qflx_snwcp_ice ,snw_rds       ,&

             mss_bcpho     ,mss_bcphi     ,mss_ocpho      ,mss_ocphi     ,&
             mss_dst1      ,mss_dst2      ,mss_dst3       ,mss_dst4      ,&

             mss_cnc_bcphi ,mss_cnc_bcpho ,mss_cnc_ocphi  ,mss_cnc_ocpho ,&
             mss_cnc_dst1  ,mss_cnc_dst2  ,mss_cnc_dst3   ,mss_cnc_dst4  )

    !
    ! !DESCRIPTION:
    ! Calculate column-integrated aerosol masses, and
    ! mass concentrations for radiative calculations and output
    ! (based on new snow level state, after SnowFilter is rebuilt.
    ! NEEDS TO BE AFTER SnowFiler is rebuilt in Hydrology2, otherwise there
    ! can be zero snow layers but an active column in filter)

    IMPLICIT NONE

    ! !ARGUMENTS:
    !
    integer, intent(in)     ::  snl              !  number of snow layers

    logical,  intent(in)    ::  do_capsnow       !  true => do snow capping
    real(r8), intent(in)    ::  h2osno_ice    ( -nlevsno+1:0 ) !  ice lens (kg/m2)
    real(r8), intent(in)    ::  h2osno_liq    ( -nlevsno+1:0 ) !  liquid water (kg/m2)
    real(r8), intent(in)    ::  qflx_snwcp_ice   !  excess snowfall due to snow capping (mm H2O /s) [+]

    real(r8), intent(inout) ::  snw_rds       ( -nlevsno+1:0 ) !  effective snow grain radius (col,lyr) [microns, m^-6]

    real(r8), intent(inout) ::  mss_bcpho     ( -nlevsno+1:0 ) !  mass of hydrophobic BC in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_bcphi     ( -nlevsno+1:0 ) !  mass of hydrophillic BC in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_ocpho     ( -nlevsno+1:0 ) !  mass of hydrophobic OC in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_ocphi     ( -nlevsno+1:0 ) !  mass of hydrophillic OC in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_dst1      ( -nlevsno+1:0 ) !  mass of dust species 1 in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_dst2      ( -nlevsno+1:0 ) !  mass of dust species 2 in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_dst3      ( -nlevsno+1:0 ) !  mass of dust species 3 in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_dst4      ( -nlevsno+1:0 ) !  mass of dust species 4 in snow (col,lyr) [kg]

    real(r8), intent(out)   ::  mss_cnc_bcphi ( -nlevsno+1:0 ) !  mass concentration of BC species 1 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_bcpho ( -nlevsno+1:0 ) !  mass concentration of BC species 2 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_ocphi ( -nlevsno+1:0 ) !  mass concentration of OC species 1 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_ocpho ( -nlevsno+1:0 ) !  mass concentration of OC species 2 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_dst1  ( -nlevsno+1:0 ) !  mass concentration of dust species 1 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_dst2  ( -nlevsno+1:0 ) !  mass concentration of dust species 2 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_dst3  ( -nlevsno+1:0 ) !  mass concentration of dust species 3 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_dst4  ( -nlevsno+1:0 ) !  mass concentration of dust species 4 (col,lyr) [kg/kg]

    ! !LOCAL VARIABLES:
    integer  :: c,j             ! indices
    real(r8) :: snowmass        ! liquid+ice snow mass in a layer [kg/m2]
    real(r8) :: snowcap_scl_fct ! temporary factor used to correct for snow capping

    !-----------------------------------------------------------------------

      !do c = begc, endc
         do j = -nlevsno+1, 0

            ! layer mass of snow:
            snowmass = h2osno_ice(j) + h2osno_liq(j)

            if (.not. use_extrasnowlayers) then
               ! Correct the top layer aerosol mass to account for snow capping.
               ! This approach conserves the aerosol mass concentration
               ! (but not the aerosol amss) when snow-capping is invoked

               if (j == snl+1) then
                  if (do_capsnow) then

                     snowcap_scl_fct = snowmass / (snowmass + (qflx_snwcp_ice*dtime))

                     mss_bcpho(j) = mss_bcpho(j)*snowcap_scl_fct
                     mss_bcphi(j) = mss_bcphi(j)*snowcap_scl_fct
                     mss_ocpho(j) = mss_ocpho(j)*snowcap_scl_fct
                     mss_ocphi(j) = mss_ocphi(j)*snowcap_scl_fct

                     mss_dst1(j)  = mss_dst1(j)*snowcap_scl_fct
                     mss_dst2(j)  = mss_dst2(j)*snowcap_scl_fct
                     mss_dst3(j)  = mss_dst3(j)*snowcap_scl_fct
                     mss_dst4(j)  = mss_dst4(j)*snowcap_scl_fct
                  endif
               endif
            endif

            if (j >= snl+1) then

               mss_cnc_bcphi(j) = mss_bcphi(j) / snowmass
               mss_cnc_bcpho(j) = mss_bcpho(j) / snowmass

               mss_cnc_ocphi(j) = mss_ocphi(j) / snowmass
               mss_cnc_ocpho(j) = mss_ocpho(j) / snowmass

               mss_cnc_dst1(j)  = mss_dst1(j)  / snowmass
               mss_cnc_dst2(j)  = mss_dst2(j)  / snowmass
               mss_cnc_dst3(j)  = mss_dst3(j)  / snowmass
               mss_cnc_dst4(j)  = mss_dst4(j)  / snowmass

            else
               ! 01/10/2023, yuan: set empty snow layers to snw_rds_min
               !snw_rds(j)       = 0._r8
               snw_rds(j)       = snw_rds_min

               mss_bcpho(j)     = 0._r8
               mss_bcphi(j)     = 0._r8
               mss_cnc_bcphi(j) = 0._r8
               mss_cnc_bcpho(j) = 0._r8

               mss_ocpho(j)     = 0._r8
               mss_ocphi(j)     = 0._r8
               mss_cnc_ocphi(j) = 0._r8
               mss_cnc_ocpho(j) = 0._r8

               mss_dst1(j)      = 0._r8
               mss_dst2(j)      = 0._r8
               mss_dst3(j)      = 0._r8
               mss_dst4(j)      = 0._r8
               mss_cnc_dst1(j)  = 0._r8
               mss_cnc_dst2(j)  = 0._r8
               mss_cnc_dst3(j)  = 0._r8
               mss_cnc_dst4(j)  = 0._r8
            endif
         enddo

      !enddo

  end subroutine AerosolMasses



  !-----------------------------------------------------------------------
  subroutine AerosolFluxes( snl, forc_aer, &
                            mss_bcphi  ,mss_bcpho  ,mss_ocphi  ,mss_ocpho ,&
                            mss_dst1   ,mss_dst2   ,mss_dst3   ,mss_dst4   )
    !
    ! !DESCRIPTION:
    ! Compute aerosol fluxes through snowpack and aerosol deposition fluxes into top layere
    !
    IMPLICIT NONE
    !
    !-----------------------------------------------------------------------
    ! !ARGUMENTS:
    integer, intent(in) :: snl    ! number of snow layers

    real(r8), intent(in) :: forc_aer (14 )  ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]

    real(r8), intent(inout) :: mss_bcphi  (-nlevsno+1:0 )  ! hydrophillic BC mass in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_bcpho  (-nlevsno+1:0 )  ! hydrophobic  BC mass in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_ocphi  (-nlevsno+1:0 )  ! hydrophillic OC mass in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_ocpho  (-nlevsno+1:0 )  ! hydrophobic  OC mass in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_dst1   (-nlevsno+1:0 )  ! mass of dust species 1 in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_dst2   (-nlevsno+1:0 )  ! mass of dust species 2 in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_dst3   (-nlevsno+1:0 )  ! mass of dust species 3 in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_dst4   (-nlevsno+1:0 )  ! mass of dust species 4 in snow (col,lyr) [kg]

    ! !LOCAL VARIABLES:
    real(r8) :: flx_bc_dep          ! total BC deposition (col) [kg m-2 s-1]
    real(r8) :: flx_bc_dep_phi      ! hydrophillic BC deposition (col) [kg m-1 s-1]
    real(r8) :: flx_bc_dep_pho      ! hydrophobic BC deposition (col) [kg m-1 s-1]
    real(r8) :: flx_oc_dep          ! total OC deposition (col) [kg m-2 s-1]
    real(r8) :: flx_oc_dep_phi      ! hydrophillic OC deposition (col) [kg m-1 s-1]
    real(r8) :: flx_oc_dep_pho      ! hydrophobic OC deposition (col) [kg m-1 s-1]
    real(r8) :: flx_dst_dep         ! total dust deposition (col) [kg m-2 s-1]

    real(r8) :: flx_dst_dep_wet1    ! wet dust (species 1) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_dry1    ! dry dust (species 1) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_wet2    ! wet dust (species 2) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_dry2    ! dry dust (species 2) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_wet3    ! wet dust (species 3) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_dry3    ! dry dust (species 3) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_wet4    ! wet dust (species 4) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_dry4    ! dry dust (species 4) deposition (col) [kg m-2 s-1]

    integer  :: c

    !-----------------------------------------------------------------------
    ! set aerosol deposition fluxes from forcing array
    ! The forcing array is either set from an external file
    ! or from fluxes received from the atmosphere model
#ifdef MODAL_AER
    ! Mapping for modal aerosol scheme where within-hydrometeor and
    ! interstitial aerosol fluxes are differentiated. Here, "phi"
    ! flavors of BC and OC correspond to within-hydrometeor
    ! (cloud-borne) aerosol, and "pho" flavors are interstitial
    ! aerosol. "wet" and "dry" fluxes of BC and OC specified here are
    ! purely diagnostic
    !do c = begc, endc

       flx_bc_dep_phi   = forc_aer(3)
       flx_bc_dep_pho   = forc_aer(1) + forc_aer(2)
       flx_bc_dep       = forc_aer(1) + forc_aer(2) + forc_aer(3)

       flx_oc_dep_phi   = forc_aer(6)
       flx_oc_dep_pho   = forc_aer(4) + forc_aer(5)
       flx_oc_dep       = forc_aer(4) + forc_aer(5) + forc_aer(6)

       flx_dst_dep_wet1 = forc_aer(7)
       flx_dst_dep_dry1 = forc_aer(8)
       flx_dst_dep_wet2 = forc_aer(9)
       flx_dst_dep_dry2 = forc_aer(10)
       flx_dst_dep_wet3 = forc_aer(11)
       flx_dst_dep_dry3 = forc_aer(12)
       flx_dst_dep_wet4 = forc_aer(13)
       flx_dst_dep_dry4 = forc_aer(14)
       flx_dst_dep      = forc_aer(7) + forc_aer(8) + forc_aer(9) + &
                             forc_aer(10) + forc_aer(11) + forc_aer(12) + &
                             forc_aer(13) + forc_aer(14)
    !end do

#else

    ! Original mapping for bulk aerosol deposition. phi and pho BC/OC
    ! species are distinguished in model, other fluxes (e.g., dry and
    ! wet BC/OC) are purely diagnostic.

      !do c = begc,endc

         flx_bc_dep_phi   = forc_aer(1) + forc_aer(3)
         flx_bc_dep_pho   = forc_aer(2)
         flx_bc_dep       = forc_aer(1) + forc_aer(2) + forc_aer(3)

         flx_oc_dep_phi   = forc_aer(4) + forc_aer(6)
         flx_oc_dep_pho   = forc_aer(5)
         flx_oc_dep       = forc_aer(4) + forc_aer(5) + forc_aer(6)

         flx_dst_dep_wet1 = forc_aer(7)
         flx_dst_dep_dry1 = forc_aer(8)
         flx_dst_dep_wet2 = forc_aer(9)
         flx_dst_dep_dry2 = forc_aer(10)
         flx_dst_dep_wet3 = forc_aer(11)
         flx_dst_dep_dry3 = forc_aer(12)
         flx_dst_dep_wet4 = forc_aer(13)
         flx_dst_dep_dry4 = forc_aer(14)
         flx_dst_dep      = forc_aer(7) + forc_aer(8) + forc_aer(9) + &
                               forc_aer(10) + forc_aer(11) + forc_aer(12) + &
                               forc_aer(13) + forc_aer(14)
      !end do
#endif

      ! aerosol deposition fluxes into top layer
      ! This is done after the inter-layer fluxes so that some aerosol
      ! is in the top layer after deposition, and is not immediately
      ! washed out before radiative calculations are done


      !do c = begc,endc
         mss_bcphi(snl+1) = mss_bcphi(snl+1) + (flx_bc_dep_phi*dtime)
         mss_bcpho(snl+1) = mss_bcpho(snl+1) + (flx_bc_dep_pho*dtime)
         mss_ocphi(snl+1) = mss_ocphi(snl+1) + (flx_oc_dep_phi*dtime)
         mss_ocpho(snl+1) = mss_ocpho(snl+1) + (flx_oc_dep_pho*dtime)

         mss_dst1(snl+1) = mss_dst1(snl+1) + (flx_dst_dep_dry1 + flx_dst_dep_wet1)*dtime
         mss_dst2(snl+1) = mss_dst2(snl+1) + (flx_dst_dep_dry2 + flx_dst_dep_wet2)*dtime
         mss_dst3(snl+1) = mss_dst3(snl+1) + (flx_dst_dep_dry3 + flx_dst_dep_wet3)*dtime
         mss_dst4(snl+1) = mss_dst4(snl+1) + (flx_dst_dep_dry4 + flx_dst_dep_wet4)*dtime
      !end do

  end subroutine AerosolFluxes


END MODULE MOD_Aerosol

