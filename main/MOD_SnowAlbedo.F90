MODULE MOD_SnowAlbedo

!-----------------------------------------------------------------------
   USE precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: SnowAlbedo


!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------


   subroutine SnowAlbedo( use_snicar_frc,use_snicar_ad ,coszen_col    ,&
                          albsod        ,albsoi        ,snl           ,frac_sno      ,&
                          h2osno        ,h2osno_liq    ,h2osno_ice    ,snw_rds       ,&

                          mss_cnc_bcphi ,mss_cnc_bcpho ,mss_cnc_ocphi ,mss_cnc_ocpho ,&
                          mss_cnc_dst1  ,mss_cnc_dst2  ,mss_cnc_dst3  ,mss_cnc_dst4  ,&

                          albgrd        ,albgri        ,albgrd_pur    ,albgri_pur    ,&
                          albgrd_bc     ,albgri_bc     ,albgrd_oc     ,albgri_oc     ,&
                          albgrd_dst    ,albgri_dst    ,flx_absdv     ,flx_absdn     ,&
                          flx_absiv     ,flx_absin      )

   ! !DESCRIPTION:
   ! The calling sequence is:
   ! -> SNICAR_RT:   snow albedos: direct beam (SNICAR)
   !    or
   !    SNICAR_AD_RT: snow albedos: direct beam (SNICAR-AD)
   ! -> SNICAR_RT:   snow albedos: diffuse (SNICAR)
   !    or
   !    SNICAR_AD_RT:   snow albedos: diffuse (SNICAR-AD)
   !
   ! ORIGINAL:
   ! 1) The Community Land Model version5.0 (CLM5.0)
   ! 2) Energy Exascale Earth System Model version 2.0 (E3SM v2.0) Land Model (ELM v2.0)
   !
   ! REFERENCES:
   ! 1) Flanner et al, 2021, SNICAR-ADv3: a community tool for modeling spectral snow albedo.
   ! Geosci. Model Dev., 14, 7673–7704, https://doi.org/10.5194/gmd-14-7673-2021
   ! 2) Hao et al., 2023, Improving snow albedo modeling in the E3SM land model (version 2.0)
   ! and assessing its impacts on snow and surface fluxes over the Tibetan Plateau.
   ! Geosci. Model Dev., 16, 75–94, https://doi.org/10.5194/gmd-16-75-2023
   !
   ! REVISIONS:
   ! Yongjiu Dai, and Hua Yuan, December, 2022 : ASSEMBLING and FITTING

   !-----------------------------------------------------------------------
   ! !USES:
     use MOD_SnowSnicar , only : SNICAR_RT, SNICAR_AD_RT

   ! and the evolution of snow effective radius
   !
   ! DAI, Dec. 28, 2022

    IMPLICIT NONE

!-------------------------------------------------------------------------
! temporay setting

    integer, parameter :: numrad  = 2            !  number of solar radiation bands: vis, nir
    integer, parameter :: nlevsno = 5            !  maximum number of snow layers

    integer, parameter :: sno_nbr_aer = 8        !  number of aerosol species in snowpack
    logical, parameter :: DO_SNO_OC   = .true.   !  parameter to include organic carbon (OC)
    logical, parameter :: DO_SNO_AER  = .true.   !  parameter to include aerosols in snowpack radiative calculations
    integer, parameter :: subgridflag = 1        !  = 0 use subgrid fluxes, = 1 not use subgrid fluxes
    !
    ! !ARGUMENTS:
    !
    logical , INTENT(in) :: use_snicar_frc       !  true: if radiative forcing is being calculated, first estimate clean-snow albedo
    logical , INTENT(in) :: use_snicar_ad        !  true: use SNICAR_AD_RT, false: use SNICAR_RT

    real(r8), INTENT(in) :: coszen_col      !  cosine of solar zenith angle
    real(r8), INTENT(in) :: albsod        ( numrad )  !  direct-beam soil albedo (col,bnd) [frc]
    real(r8), INTENT(in) :: albsoi        ( numrad )  !  diffuse soil albedo (col,bnd) [frc]

    integer , INTENT(in) :: snl             !  negative number of snow layers (col) [nbr]
    real(r8), INTENT(in) :: frac_sno        ! fraction of ground covered by snow (0 to 1)
    real(r8), INTENT(in) :: h2osno          ! snow water equivalent (mm H2O)
    real(r8), INTENT(in) :: h2osno_liq    ( -nlevsno+1:0 )  ! liquid water content (col,lyr) [kg/m2]
    real(r8), INTENT(in) :: h2osno_ice    ( -nlevsno+1:0 )  ! ice lens content (col,lyr) [kg/m2]
    real(r8), INTENT(in) :: snw_rds       ( -nlevsno+1:0 )  ! snow grain radius (col,lyr) [microns]

    real(r8), INTENT(in) :: mss_cnc_bcphi ( -nlevsno+1:0 )  !  mass concentration of hydrophilic BC (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_bcpho ( -nlevsno+1:0 )  !  mass concentration of hydrophobic BC (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_ocphi ( -nlevsno+1:0 )  !  mass concentration of hydrophilic OC (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_ocpho ( -nlevsno+1:0 )  !  mass concentration of hydrophobic OC (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_dst1  ( -nlevsno+1:0 )  !  mass concentration of dust aerosol species 1 (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_dst2  ( -nlevsno+1:0 )  !  mass concentration of dust aerosol species 2 (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_dst3  ( -nlevsno+1:0 )  !  mass concentration of dust aerosol species 3 (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_dst4  ( -nlevsno+1:0 )  !  mass concentration of dust aerosol species 4 (col,lyr) [kg/kg]

    real(r8), INTENT(out) :: albgrd       ( numrad )  !  ground albedo (direct)
    real(r8), INTENT(out) :: albgri       ( numrad )  !  ground albedo (diffuse)
    real(r8), INTENT(out) :: albgrd_pur   ( numrad )  !  pure snow ground albedo (direct)
    real(r8), INTENT(out) :: albgri_pur   ( numrad )  !  pure snow ground albedo (diffuse)
    real(r8), INTENT(out) :: albgrd_bc    ( numrad )  !  ground albedo without BC (direct)
    real(r8), INTENT(out) :: albgri_bc    ( numrad )  !  ground albedo without BC (diffuse)
    real(r8), INTENT(out) :: albgrd_oc    ( numrad )  !  ground albedo without OC (direct)
    real(r8), INTENT(out) :: albgri_oc    ( numrad )  !  ground albedo without OC (diffuse)
    real(r8), INTENT(out) :: albgrd_dst   ( numrad )  !  ground albedo without dust (direct)
    real(r8), INTENT(out) :: albgri_dst   ( numrad )  !  ground albedo without dust (diffuse)
    real(r8), INTENT(out) :: flx_absdv    ( -nlevsno+1:1 )  !  direct flux absorption factor (col,lyr): VIS [frc]
    real(r8), INTENT(out) :: flx_absdn    ( -nlevsno+1:1 )  !  direct flux absorption factor (col,lyr): NIR [frc]
    real(r8), INTENT(out) :: flx_absiv    ( -nlevsno+1:1 )  !  diffuse flux absorption factor (col,lyr): VIS [frc]
    real(r8), INTENT(out) :: flx_absin    ( -nlevsno+1:1 )  !  diffuse flux absorption factor (col,lyr): NIR [frc]

  !-----------------------------------------------------------------------
    !
    ! !LOCAL VARIABLES:
    integer  :: i            ! index for layers [idx]
    integer  :: aer          ! index for sno_nbr_aer
    integer  :: ib           ! band index
    integer  :: ic           ! 0=unit incoming direct; 1=unit incoming diffuse
    integer  :: flg_slr      ! flag for SNICAR (=1 if direct, =2 if diffuse)
    integer  :: flg_snw_ice  ! flag for SNICAR (=1 when called from ELM, =2 when called from sea-ice)

    real(r8) :: mss_cnc_aer_in_frc_pur (-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for forcing calculation (zero) (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_bc  (-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for BC forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_oc  (-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for OC forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_dst (-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for dust forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_fdb     (-nlevsno+1:0,sno_nbr_aer) ! mass concentration of all aerosol species for feedback calculation (col,lyr,aer) [kg kg-1]

    real(r8) :: albsfc       (numrad)               ! albedo of surface underneath snow (col,bnd)
    real(r8) :: albsnd       (numrad)               ! snow albedo (direct)
    real(r8) :: albsni       (numrad)               ! snow albedo (diffuse)
    real(r8) :: albsnd_pur   (numrad)               ! direct pure snow albedo (radiative forcing)
    real(r8) :: albsni_pur   (numrad)               ! diffuse pure snow albedo (radiative forcing)
    real(r8) :: albsnd_bc    (numrad)               ! direct snow albedo without BC (radiative forcing)
    real(r8) :: albsni_bc    (numrad)               ! diffuse snow albedo without BC (radiative forcing)
    real(r8) :: albsnd_oc    (numrad)               ! direct snow albedo without OC (radiative forcing)
    real(r8) :: albsni_oc    (numrad)               ! diffuse snow albedo without OC (radiative forcing)
    real(r8) :: albsnd_dst   (numrad)               ! direct snow albedo without dust (radiative forcing)
    real(r8) :: albsni_dst   (numrad)               ! diffuse snow albedo without dust (radiative forcing)
    real(r8) :: flx_absd_snw (-nlevsno+1:1,numrad)  ! flux absorption factor for just snow (direct) [frc]
    real(r8) :: flx_absi_snw (-nlevsno+1:1,numrad)  ! flux absorption factor for just snow (diffuse) [frc]
    real(r8) :: foo_snw      (-nlevsno+1:1,numrad)  ! dummy array for forcing calls

    integer  :: snw_rds_in   (-nlevsno+1:0)         ! snow grain size sent to SNICAR (col,lyr) [microns]

    integer , parameter :: nband =numrad   ! number of solar radiation waveband classes

  !-----------------------------------------------------------------------

    ! Initialize output because solar radiation only done if coszen > 0

    do ib = 1, numrad
       albgrd(ib)     = 0._r8
       albgri(ib)     = 0._r8
       albgrd_pur(ib) = 0._r8
       albgri_pur(ib) = 0._r8
       albgrd_bc(ib)  = 0._r8
       albgri_bc(ib)  = 0._r8
       albgrd_oc(ib)  = 0._r8
       albgri_oc(ib)  = 0._r8
       albgrd_dst(ib) = 0._r8
       albgri_dst(ib) = 0._r8
       do i=-nlevsno+1,1,1
          flx_absdv(i) = 0._r8
          flx_absdn(i) = 0._r8
          flx_absiv(i) = 0._r8
          flx_absin(i) = 0._r8
       enddo
    end do  ! end of numrad loop

    ! set variables to pass to SNICAR.

    flg_snw_ice = 1
    albsfc(:)     = albsoi(:)
    snw_rds_in(:) = nint(snw_rds(:))

    ! zero aerosol input arrays
    do aer = 1, sno_nbr_aer
       do i = -nlevsno+1, 0
          mss_cnc_aer_in_frc_pur(i,aer) = 0._r8
          mss_cnc_aer_in_frc_bc(i,aer)  = 0._r8
          mss_cnc_aer_in_frc_oc(i,aer)  = 0._r8
          mss_cnc_aer_in_frc_dst(i,aer) = 0._r8
          mss_cnc_aer_in_fdb(i,aer)     = 0._r8
       end do
    end do

    ! If radiative forcing is being calculated, first estimate clean-snow albedo

    if (use_snicar_frc) then

       ! 1. PURE SNOW ALBEDO CALCULATIONS
          flg_slr = 1  ! direct-beam
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_pur(:, :), &
                             albsfc(:), &
                             albsnd_pur(:), &
                             foo_snw(:, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_pur(:, :), &
                             albsfc(:), &
                             albsnd_pur(:), &
                             foo_snw(:, :) )
          endif ! end if use_snicar_ad

          flg_slr = 2  ! diffuse
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_pur(:, :), &
                             albsfc(:), &
                             albsni_pur(:), &
                             foo_snw(:, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_pur(:, :), &
                             albsfc(:), &
                             albsni_pur(:), &
                             foo_snw(:, :) )
          endif ! end if use_snicar_ad

       ! 2. BC input array:
       !  set dust and (optionally) OC concentrations, so BC_FRC=[(BC+OC+dust)-(OC+dust)]
       if (DO_SNO_OC) then
          mss_cnc_aer_in_frc_bc(:,3) = mss_cnc_ocphi(:)
          mss_cnc_aer_in_frc_bc(:,4) = mss_cnc_ocpho(:)
       endif
       mss_cnc_aer_in_frc_bc(:,5) = mss_cnc_dst1(:)
       mss_cnc_aer_in_frc_bc(:,6) = mss_cnc_dst2(:)
       mss_cnc_aer_in_frc_bc(:,7) = mss_cnc_dst3(:)
       mss_cnc_aer_in_frc_bc(:,8) = mss_cnc_dst4(:)

       ! BC FORCING CALCULATIONS
       flg_slr = 1  ! direct-beam
       if (use_snicar_ad) then
           call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_bc(:, :), &
                             albsfc(:), &
                             albsnd_bc(:), &
                             foo_snw(:, :) )
       else
           call SNICAR_RT   (flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_bc(:, :), &
                             albsfc(:), &
                             albsnd_bc(:), &
                             foo_snw(:, :) )
       endif ! end if use_snicar_ad

       flg_slr = 2  ! diffuse
       if (use_snicar_ad) then
           call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_bc(:, :), &
                             albsfc(:), &
                             albsni_bc(:), &
                             foo_snw(:, :) )
       else
           call SNICAR_RT   (flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_bc(:, :), &
                             albsfc(:), &
                             albsni_bc(:), &
                             foo_snw(:, :) )
       endif ! end if use_snicar_ad

       ! 3. OC input array:
       !  set BC and dust concentrations, so OC_FRC=[(BC+OC+dust)-(BC+dust)]
       if (DO_SNO_OC) then
          mss_cnc_aer_in_frc_oc(:,1) = mss_cnc_bcphi(:)
          mss_cnc_aer_in_frc_oc(:,2) = mss_cnc_bcpho(:)

          mss_cnc_aer_in_frc_oc(:,5) = mss_cnc_dst1(:)
          mss_cnc_aer_in_frc_oc(:,6) = mss_cnc_dst2(:)
          mss_cnc_aer_in_frc_oc(:,7) = mss_cnc_dst3(:)
          mss_cnc_aer_in_frc_oc(:,8) = mss_cnc_dst4(:)

       ! OC FORCING CALCULATIONS
          flg_slr = 1  ! direct-beam
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_oc(:, :), &
                             albsfc(:), &
                             albsnd_oc(:), &
                             foo_snw(:, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_oc(:, :), &
                             albsfc(:), &
                             albsnd_oc(:), &
                             foo_snw(:, :) )
          endif ! end if use_snicar_ad

          flg_slr = 2  ! diffuse
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_oc(:, :), &
                             albsfc(:), &
                             albsni_oc(:), &
                             foo_snw(:, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_oc(:, :), &
                             albsfc(:), &
                             albsni_oc(:), &
                             foo_snw(:, :) )
          endif ! end if use_snicar_ad
       endif  ! end if (DO_SNO_OC)

       ! 4. DUST FORCING CALCULATIONS
          ! DUST input array:
          ! set BC and OC concentrations, so DST_FRC=[(BC+OC+dust)-(BC+OC)]
          mss_cnc_aer_in_frc_dst(:,1) = mss_cnc_bcphi(:)
          mss_cnc_aer_in_frc_dst(:,2) = mss_cnc_bcpho(:)

          if (DO_SNO_OC) then
              mss_cnc_aer_in_frc_dst(:,3) = mss_cnc_ocphi(:)
              mss_cnc_aer_in_frc_dst(:,4) = mss_cnc_ocpho(:)
          endif

          flg_slr = 1  ! direct-beam
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_dst(:, :), &
                             albsfc(:), &
                             albsnd_dst(:), &
                             foo_snw(:, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_dst(:, :), &
                             albsfc(:), &
                             albsnd_dst(:), &
                             foo_snw(:, :) )
          endif ! end if use_snicar_ad

          flg_slr = 2  ! diffuse
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_dst(:, :), &
                             albsfc(:), &
                             albsni_dst(:), &
                             foo_snw(:, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col, &
                             snl, &
                             h2osno, &
                             frac_sno, &
                             h2osno_liq(:), &
                             h2osno_ice(:), &
                             snw_rds_in(:), &
                             mss_cnc_aer_in_frc_dst(:, :), &
                             albsfc(:), &
                             albsni_dst(:), &
                             foo_snw(:, :)  )
          endif ! end if use_snicar_ad

    end if !end if use_snicar_frc


    ! --------------------------------------------
    ! CLIMATE FEEDBACK CALCULATIONS, ALL AEROSOLS:
    ! --------------------------------------------
    ! Set aerosol input arrays
    ! feedback input arrays have been zeroed
    ! set soot and dust aerosol concentrations:
    if (DO_SNO_AER) then
        mss_cnc_aer_in_fdb(:,1) = mss_cnc_bcphi(:)
        mss_cnc_aer_in_fdb(:,2) = mss_cnc_bcpho(:)

        ! DO_SNO_OC is set in SNICAR_varpar. Default case is to ignore OC concentrations because:
        !  1) Knowledge of their optical properties is primitive
        !  2) When 'water-soluble' OPAC optical properties are applied to OC in snow,
        !     it has a negligible darkening effect.
        if (DO_SNO_OC) then
           mss_cnc_aer_in_fdb(:,3) = mss_cnc_ocphi(:)
           mss_cnc_aer_in_fdb(:,4) = mss_cnc_ocpho(:)
        endif

        mss_cnc_aer_in_fdb(:,5) = mss_cnc_dst1(:)
        mss_cnc_aer_in_fdb(:,6) = mss_cnc_dst2(:)
        mss_cnc_aer_in_fdb(:,7) = mss_cnc_dst3(:)
        mss_cnc_aer_in_fdb(:,8) = mss_cnc_dst4(:)
    endif

    flg_slr = 1  ! direct-beam
    if (use_snicar_ad) then
        call SNICAR_AD_RT(flg_snw_ice, &
                          flg_slr, &
                          coszen_col, &
                          snl, &
                          h2osno, &
                          frac_sno, &
                          h2osno_liq(:), &
                          h2osno_ice(:), &
                          snw_rds_in(:), &
                          mss_cnc_aer_in_fdb(:, :), &
                          albsfc(:), &
                          albsnd(:), &
                          flx_absd_snw(:, :) )
    else
        call SNICAR_RT   (flg_snw_ice, &
                          flg_slr, &
                          coszen_col, &
                          snl, &
                          h2osno, &
                          frac_sno, &
                          h2osno_liq(:), &
                          h2osno_ice(:), &
                          snw_rds_in(:), &
                          mss_cnc_aer_in_fdb(:, :), &
                          albsfc(:), &
                          albsnd(:), &
                          flx_absd_snw(:, :) )
    endif ! end if use_snicar_ad

    flg_slr = 2  ! diffuse
    if (use_snicar_ad) then
        call SNICAR_AD_RT(flg_snw_ice, &
                          flg_slr, &
                          coszen_col, &
                          snl, &
                          h2osno, &
                          frac_sno, &
                          h2osno_liq(:), &
                          h2osno_ice(:), &
                          snw_rds_in(:), &
                          mss_cnc_aer_in_fdb(:, :), &
                          albsfc(:), &
                          albsni(:), &
                          flx_absi_snw(:, :) )
    else
        call SNICAR_RT   (flg_snw_ice, &
                          flg_slr, &
                          coszen_col, &
                          snl, &
                          h2osno, &
                          frac_sno, &
                          h2osno_liq(:), &
                          h2osno_ice(:), &
                          snw_rds_in(:), &
                          mss_cnc_aer_in_fdb(:, :), &
                          albsfc(:), &
                          albsni(:), &
                          flx_absi_snw(:, :) )
    endif ! end if use_snicar_ad


    ! ground albedos and snow-fraction weighting of snow absorption factors
    do ib = 1, nband
       if (coszen_col > 0._r8) then
          ! ground albedo was originally computed in SoilAlbedo, but is now computed here
          ! because the order of SoilAlbedo and SNICAR_RT/SNICAR_AD_RT was switched for SNICAR/SNICAR_AD_RT.
          albgrd(ib) = albsod(ib)*(1._r8-frac_sno) + albsnd(ib)*frac_sno
          albgri(ib) = albsoi(ib)*(1._r8-frac_sno) + albsni(ib)*frac_sno

          ! albedos for radiative forcing calculations:
          if (use_snicar_frc) then
             ! pure snow albedo for all-aerosol radiative forcing
             albgrd_pur(ib) = albsod(ib)*(1.-frac_sno) + albsnd_pur(ib)*frac_sno
             albgri_pur(ib) = albsoi(ib)*(1.-frac_sno) + albsni_pur(ib)*frac_sno

             ! BC forcing albedo
             albgrd_bc(ib) = albsod(ib)*(1.-frac_sno) + albsnd_bc(ib)*frac_sno
             albgri_bc(ib) = albsoi(ib)*(1.-frac_sno) + albsni_bc(ib)*frac_sno

             if (DO_SNO_OC) then
                ! OC forcing albedo
                albgrd_oc(ib) = albsod(ib)*(1.-frac_sno) + albsnd_oc(ib)*frac_sno
                albgri_oc(ib) = albsoi(ib)*(1.-frac_sno) + albsni_oc(ib)*frac_sno
             endif

             ! dust forcing albedo
             albgrd_dst(ib) = albsod(ib)*(1.-frac_sno) + albsnd_dst(ib)*frac_sno
             albgri_dst(ib) = albsoi(ib)*(1.-frac_sno) + albsni_dst(ib)*frac_sno
          end if

          ! also in this loop (but optionally in a different loop for vectorized code)
          !  weight snow layer radiative absorption factors based on snow fraction and soil albedo
          !  (NEEDED FOR ENERGY CONSERVATION)
          do i = -nlevsno+1,1,1
             if (subgridflag == 0 ) then
                if (ib == 1) then
                   flx_absdv(i) = flx_absd_snw(i,ib)*frac_sno + &
                        ((1.-frac_sno)*(1-albsod(ib))*(flx_absd_snw(i,ib)/(1.-albsnd(ib))))
                   flx_absiv(i) = flx_absi_snw(i,ib)*frac_sno + &
                        ((1.-frac_sno)*(1-albsoi(ib))*(flx_absi_snw(i,ib)/(1.-albsni(ib))))
                elseif (ib == 2) then
                   flx_absdn(i) = flx_absd_snw(i,ib)*frac_sno + &
                        ((1.-frac_sno)*(1-albsod(ib))*(flx_absd_snw(i,ib)/(1.-albsnd(ib))))
                   flx_absin(i) = flx_absi_snw(i,ib)*frac_sno + &
                        ((1.-frac_sno)*(1-albsoi(ib))*(flx_absi_snw(i,ib)/(1.-albsni(ib))))
                endif
             else
                if (ib == 1) then
                   flx_absdv(i) = flx_absd_snw(i,ib)*(1.-albsnd(ib))
                   flx_absiv(i) = flx_absi_snw(i,ib)*(1.-albsni(ib))
                elseif (ib == 2) then
                   flx_absdn(i) = flx_absd_snw(i,ib)*(1.-albsnd(ib))
                   flx_absin(i) = flx_absi_snw(i,ib)*(1.-albsni(ib))
                endif
             endif
          enddo
       endif
    enddo

    end subroutine SnowAlbedo

END MODULE MOD_SnowAlbedo
