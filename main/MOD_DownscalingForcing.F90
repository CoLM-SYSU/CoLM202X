#include <define.h>

module MOD_DownscalingForcing
  !-----------------------------------------------------------------------
  ! DESCRIPTION:
  ! Downscaling of gridcell-level meteorological forcings into column-level forcings
  ! (not included wind and solar radiation)
  !
  ! INITIAL:
  ! The Community Land Model version 5.0 (CLM5.0)
  !
  ! REVISIONS:
  ! Updated by Yongjiu Dai, January 16, 2023
  !

  USE precision
  use MOD_Qsadv
  IMPLICIT NONE

  real(r8), parameter :: SHR_CONST_MWDAIR = 28.966_r8       ! molecular weight dry air [kg/kmole]
  real(r8), parameter :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
  real(r8), parameter :: SHR_CONST_AVOGAD = 6.02214e26_r8   ! Avogadro's number [molecules/kmole]
  real(r8), parameter :: SHR_CONST_BOLTZ  = 1.38065e-23_r8  ! Boltzmann's constant [J/K/molecule]
  real(r8), parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ  ! Universal gas constant [J/K/kmole]
  real(r8), parameter :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR   ! Dry air gas constant [J/K/kg]
  real(R8), parameter :: SHR_CONST_TKFRZ  = 273.15_r8       ! freezing T of fresh water [K]

  real(r8), parameter :: cpair  = 1.00464e3_r8  ! specific heat of dry air [J/kg/K]
  real(r8), parameter :: grav   = 9.80616_r8    ! acceleration of gravity [m/s^2]
  real(r8), parameter :: denh2o = 1.000e3_r8    ! density of liquid water [kg/m3]
  real(r8), parameter :: hfus   = 3.337e5_r8    ! latent heat of fusion for ice [J/kg]
  real(r8) :: rair = SHR_CONST_RDAIR  ! Dry air gas constant [J/K/kg]
  real(r8) :: tfrz = SHR_CONST_TKFRZ  ! freezing T of fresh water [K]

  ! On the windward side of the range, annual mean lapse rates of 3.9-5.2 (deg km-1),
  ! substantially smaller than the often-assumed 6.5 (deg km-1).
  ! The data sets show similar seasonal and diurnal variability, with lapse rates smallest
  ! (2.5-3.5 deg km-1) in late-summer minimum temperatures, and largest (6.5-7.5 deg km-1)
  ! in spring maximum temperatures. Geographic (windward versus lee side) differences in
  ! lapse rates are substantial. [Minder et al., 2010, Surface temperature lapse rates over complex terrain:
  ! Lessons from the Cascade Mountains. JGR, 115, doi:10.1029/2009JD013493]
  !
  ! Kunkel, K. E., 1989: Simple procedures for extrapolation of humidity variables in the mountainous western United States.
  ! J. Climate, 2, 656-669. lapse_rate = /Jan - Dec/4.4,5.9,7.1,7.8,8.1,8.2,8.1,8.1,7.7,6.8,5.5,4.7/ (deg km-1)

  real(r8), parameter :: lapse_rate = 0.006_r8 ! surface temperature lapse rate (deg m-1)

  SAVE
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: downscale_forcings   ! Downscale atm forcing fields from gridcell to column

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: rhos                ! calculate atmospheric density
  private :: downscale_longwave  ! Downscale longwave radiation from gridcell to column

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine downscale_forcings(num_gridcells,num_columns,begc,endc,glaciers,wt_column,&

             forc_topo_g ,forc_t_g   ,forc_th_g  ,forc_q_g     ,forc_pbot_g ,&
             forc_rho_g  ,forc_prc_g ,forc_prl_g ,forc_lwrad_g ,forc_hgt_grc,&

             forc_topo_c ,forc_t_c   ,forc_th_c  ,forc_q_c     ,forc_pbot_c ,&
             forc_rho_c  ,forc_prc_c ,forc_prl_c ,forc_lwrad_c )
    !
    ! !DESCRIPTION:
    ! Downscale atmospheric forcing fields from gridcell to column.
    !
    ! Downscaling is done based on the difference between each land model column's elevation and
    ! the atmosphere's surface elevation (which is the elevation at which the atmospheric
    ! forcings are valid).
    !
    ! Note that the downscaling procedure can result in changes in grid cell mean values
    ! compared to what was provided by the atmosphere. We conserve fluxes of mass and
    ! energy, but allow states such as temperature to differ.
    !
    ! !USES:
    !
    IMPLICIT NONE

    ! !ARGUMENTS:
    integer,  INTENT(in) :: num_gridcells  ! number of gridcell
    integer,  INTENT(in) :: num_columns    ! number of column
    integer,  INTENT(in) :: begc (1:num_gridcells) ! beginning column index
    integer,  INTENT(in) :: endc (1:num_gridcells) ! ending column index
    logical,  INTENT(in) :: glaciers (1:num_columns) ! true: glacier column (itypwat = 3)
    real(r8), INTENT(in) :: wt_column(1:num_columns) ! weight of the column relative to the grid cell

    ! Gridcell-level non-downscaled fields:
    real(r8), INTENT(in) :: forc_topo_g  (1:num_gridcells) ! atmospheric surface height [m]
    real(r8), INTENT(in) :: forc_t_g     (1:num_gridcells) ! atmospheric temperature [Kelvin]
    real(r8), INTENT(in) :: forc_th_g    (1:num_gridcells) ! atmospheric potential temperature [Kelvin]
    real(r8), INTENT(in) :: forc_q_g     (1:num_gridcells) ! atmospheric specific humidity [kg/kg]
    real(r8), INTENT(in) :: forc_pbot_g  (1:num_gridcells) ! atmospheric pressure [Pa]
    real(r8), INTENT(in) :: forc_rho_g   (1:num_gridcells) ! atmospheric density [kg/m**3]
    real(r8), INTENT(in) :: forc_prc_g   (1:num_gridcells) ! convective precipitation in grid [mm/s]
    real(r8), INTENT(in) :: forc_prl_g   (1:num_gridcells) ! large-scale precipitation in grid [mm/s]
    real(r8), INTENT(in) :: forc_lwrad_g (1:num_gridcells) ! grid downward longwave [W/m**2]
    real(r8), INTENT(in) :: forc_hgt_grc (1:num_gridcells) ! atmospheric reference height [m]

    ! Column-level downscaled fields:
    real(r8), INTENT(in)  :: forc_topo_c (1:num_columns) ! column surface height [m]
    real(r8), INTENT(out) :: forc_t_c    (1:num_columns) ! atmospheric temperature [Kelvin]
    real(r8), INTENT(out) :: forc_th_c   (1:num_columns) ! atmospheric potential temperature [Kelvin]
    real(r8), INTENT(out) :: forc_q_c    (1:num_columns) ! atmospheric specific humidity [kg/kg]
    real(r8), INTENT(out) :: forc_pbot_c (1:num_columns) ! atmospheric pressure [Pa]
    real(r8), INTENT(out) :: forc_rho_c  (1:num_columns) ! atmospheric density [kg/m**3]
    real(r8), INTENT(out) :: forc_prc_c  (1:num_columns) ! column convective precipitation [mm/s]
    real(r8), INTENT(out) :: forc_prl_c  (1:num_columns) ! column large-scale precipitation [mm/s]
    real(r8), INTENT(out) :: forc_lwrad_c(1:num_columns) ! column downward longwave [W/m**2]
    !
    ! !LOCAL VARIABLES:
    integer :: g,c   ! indices

    ! !Local varaiables for topo downscaling:
    real(r8) :: hsurf_g, hsurf_c
    real(r8) :: Hbot, zbot
    real(r8) :: tbot_g, pbot_g, thbot_g, qbot_g, qs_g, es_g, rhos_g
    real(r8) :: tbot_c, pbot_c, thbot_c, qbot_c, qs_c, es_c, rhos_c
    real(r8) :: rhos_c_estimate, rhos_g_estimate
    real(r8) :: dum1, dum2

    real(r8) :: max_elev_c    ! the maximum column level elevation value within the grid
    real(r8) :: delta_prc_c   ! deviation of the column convective precipitation from the grid level precipitation
    real(r8) :: delta_prl_c   ! deviation of the column large-scale precipitation from the grid level precipitation

    !-----------------------------------------------------------------------

    ! Initialize column forcing (needs to be done for ALL active columns)
    do g = 1, num_gridcells
       do c = begc(g), endc(g)
          forc_t_c    (c) = forc_t_g    (g)
          forc_th_c   (c) = forc_th_g   (g)
          forc_q_c    (c) = forc_q_g    (g)
          forc_pbot_c (c) = forc_pbot_g (g)
          forc_rho_c  (c) = forc_rho_g  (g)

          forc_prc_c  (c) = forc_prc_g  (g)
          forc_prl_c  (c) = forc_prl_g  (g)
          forc_lwrad_c(c) = forc_lwrad_g(g)
       end do
    end do

    ! Downscale forc_t, forc_th, forc_q, forc_pbot, and forc_rho to columns.
    do g = 1, num_gridcells

       hsurf_g = forc_topo_g(g)             ! gridcell sfc elevation
       tbot_g  = forc_t_g(g)                ! atm sfc temp
       thbot_g = forc_th_g(g)               ! atm sfc pot temp
       qbot_g  = forc_q_g(g)                ! atm sfc spec humid
       pbot_g  = forc_pbot_g(g)             ! atm sfc pressure
       rhos_g  = forc_rho_g(g)              ! atm density
       zbot    = forc_hgt_grc(g)            ! atm ref height

       max_elev_c = maxval(forc_topo_c(begc(g):endc(g))) ! maximum column level elevation value within the grid

       do c = begc(g), endc(g)
          ! This is a simple downscaling procedure
          ! Note that forc_hgt, forc_u, forc_v and solar radiation are not downscaled.

          hsurf_c = forc_topo_c(c)                        ! column sfc elevation
          tbot_c  = tbot_g-lapse_rate*(hsurf_c-hsurf_g)   ! sfc temp for column
          Hbot    = rair*0.5_r8*(tbot_g+tbot_c)/grav      ! scale ht at avg temp
          pbot_c  = pbot_g*exp(-(hsurf_c-hsurf_g)/Hbot)   ! column sfc press

          ! Derivation of potential temperature calculation:
          !
          ! The textbook definition would be:
          ! thbot_c = tbot_c * (p0/pbot_c)^(rair/cpair)
          !
          ! Note that pressure is related to scale height as:
          ! pbot_c = p0 * exp(-zbot/Hbot)
          !
          ! Plugging this in to the textbook definition, then manipulating, we get:
          ! thbot_c = tbot_c * (p0/(p0*exp(-zbot/Hbot)))^(rair/cpair)
          !         = tbot_c * (1/exp(-zbot/Hbot))^(rair/cpair)
          !         = tbot_c * (exp(zbot/Hbot))^(rair/cpair)
          !         = tbot_c * exp((zbot/Hbot) * (rair/cpair))

          ! But we want everything expressed in delta form, resulting in:
          thbot_c = thbot_g + (tbot_c - tbot_g)*exp((zbot/Hbot)*(rair/cpair))

          call Qsadv(tbot_g,pbot_g,es_g,dum1,qs_g,dum2)
          call Qsadv(tbot_c,pbot_c,es_c,dum1,qs_c,dum2)

          qbot_c = qbot_g*(qs_c/qs_g)
          rhos_c_estimate = rhos(qbot=qbot_c, pbot=pbot_c, tbot=tbot_c)
          rhos_g_estimate = rhos(qbot=qbot_g, pbot=pbot_g, tbot=tbot_g)
          rhos_c = rhos_g * (rhos_c_estimate / rhos_g_estimate)

          forc_t_c(c)    = tbot_c
          forc_th_c(c)   = thbot_c
          forc_q_c(c)    = qbot_c
          forc_pbot_c(c) = pbot_c
          forc_rho_c(c)  = rhos_c

#if(defined option_precipitation_adjust_I)
          ! Tesfa et al, 2020: Exploring Topography-Based Methods for Downscaling
          ! Subgrid Precipitation for Use in Earth System Models. Equation (5)
          ! https://doi.org/ 10.1029/2019JD031456

          delta_prc_c = forc_prc_g(g) * (forc_topo_c(c) - forc_topo_g(g)) / max_elev_c
          forc_prc_c(c) = forc_prc_g(g) + delta_prc_c   ! convective precipitation [mm/s]

          delta_prl_c = forc_prl_g(g) * (forc_topo_c(c) - forc_topo_g(g)) / max_elev_c
          forc_prl_c(c) = forc_prl_g(g) + delta_prl_c   ! large scale precipitation [mm/s]

#elif(defined option_precipitation_adjust_II)
          ! Liston, G. E. and Elder, K.: A meteorological distribution system
          ! for high-resolution terrestrial modeling (MicroMet), J. Hydrometeorol., 7, 217-234, 2006.
          ! Equation (33) and Table 1: chi range from January to December:
          ! [0.35,0.35,0.35,0.30,0.25,0.20,0.20,0.20,0.20,0.25,0.30,0.35] (1/m)

          delta_prc_c = forc_prc_g(g) * 2.0*0.27e-3*(forc_topo_c(c) - forc_topo_g(g)) &
                                    /(1.0 - 0.27e-3*(forc_topo_c(c) - forc_topo_g(g)))
          forc_prc_c(c) = forc_prc_g(g) + delta_prc_c   ! large scale precipitation [mm/s]

          delta_prl_c = forc_prl_g(g) * 2.0*0.27e-3*(forc_topo_c(c) - forc_topo_g(g)) &
                                    /(1.0 - 0.27e-3*(forc_topo_c(c) - forc_topo_g(g)))
          forc_prl_c(c) = forc_prl_g(g) + delta_prl_c   ! large scale precipitation [mm/s]
#else

#endif

          IF (forc_prl_c(c) < 0) THEN
             write(*,*) 'negative prl', forc_prl_g(g), max_elev_c, forc_topo_c(c), forc_topo_g(g)
          ENDIF

          IF (forc_prc_c(c) < 0) THEN
             write(*,*) 'negative prc', forc_prc_g(g), max_elev_c, forc_topo_c(c), forc_topo_g(g)
          ENDIF

          forc_prc_c(c) = max(forc_prc_c(c), 0.)
          forc_prl_c(c) = max(forc_prl_c(c), 0.)

       end do
    end do

    call downscale_longwave(num_gridcells, num_columns, begc, endc, glaciers, wt_column, &
                   forc_topo_g, forc_t_g, forc_q_g, forc_pbot_g, forc_lwrad_g, &
                   forc_topo_c, forc_t_c, forc_q_c, forc_pbot_c, forc_lwrad_c)

  end subroutine downscale_forcings


  !-----------------------------------------------------------------------
  pure function rhos(qbot, pbot, tbot)
    !
    ! !DESCRIPTION:
    ! Compute atmospheric density (kg/m**3)
    !
    IMPLICIT NONE
    !
    ! !ARGUMENTS:
    real(r8) :: rhos  ! function result: atmospheric density (kg/m**3)
    real(r8), intent(in) :: qbot  ! atmospheric specific humidity (kg/kg)
    real(r8), intent(in) :: pbot  ! atmospheric pressure (Pa)
    real(r8), intent(in) :: tbot  ! atmospheric temperature (K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: egcm
    real(r8) :: wv_to_dair_weight_ratio  ! ratio of molecular weight of water vapor to that of dry air [-]

    !-----------------------------------------------------------------------
    wv_to_dair_weight_ratio = SHR_CONST_MWWV/SHR_CONST_MWDAIR

    egcm = qbot*pbot / (wv_to_dair_weight_ratio + (1._r8 - wv_to_dair_weight_ratio)*qbot)
    rhos = (pbot - (1._r8 - wv_to_dair_weight_ratio)*egcm) / (rair*tbot)

  end function rhos


  ! -----------------------------------------------------------------------
  subroutine downscale_longwave(num_gridcells, num_columns, begc, endc, glaciers, wt_column, &
                       forc_topo_g, forc_t_g, forc_q_g, forc_pbot_g, forc_lwrad_g, &
                       forc_topo_c, forc_t_c, forc_q_c, forc_pbot_c, forc_lwrad_c)
    !
    ! !DESCRIPTION:
    ! Downscale longwave radiation from gridcell to column
    ! Must be done AFTER temperature downscaling

    IMPLICIT NONE

    ! !ARGUMENTS:
    integer,  INTENT(in) :: num_gridcells  ! number of gridcell
    integer,  INTENT(in) :: num_columns    ! number of column
    integer,  INTENT(in) :: begc (1:num_gridcells) ! beginning column index
    integer,  INTENT(in) :: endc (1:num_gridcells) ! ending column index
    logical,  INTENT(in) :: glaciers  (1:num_columns) ! true: glacier column
    real(r8), INTENT(in) :: wt_column (1:num_columns) ! weight of the column relative to the grid cell

    real(r8), INTENT(in) :: forc_topo_g  (1:num_gridcells)  ! atmospheric surface height (m)
    real(r8), INTENT(in) :: forc_t_g     (1:num_gridcells)  ! atmospheric temperature [Kelvin]
    real(r8), INTENT(in) :: forc_q_g     (1:num_gridcells) ! atmospheric specific humidity [kg/kg]
    real(r8), INTENT(in) :: forc_pbot_g  (1:num_gridcells) ! atmospheric pressure [Pa]
    real(r8), INTENT(in) :: forc_lwrad_g (1:num_gridcells)  ! downward longwave (W/m**2)
    real(r8), INTENT(in) :: forc_topo_c  (1:num_columns) ! column surface height (m)
    real(r8), INTENT(in) :: forc_t_c     (1:num_columns) ! atmospheric temperature [Kelvin]
    real(r8), INTENT(in) :: forc_q_c     (1:num_columns) ! atmospheric specific humidity [kg/kg]
    real(r8), INTENT(in) :: forc_pbot_c  (1:num_columns) ! atmospheric pressure [Pa]
    real(r8), INTENT(out) :: forc_lwrad_c(1:num_columns) ! downward longwave (W/m**2)


    ! !LOCAL VARIABLES:
    integer  :: c,g      ! indices
    real(r8) :: hsurf_c  ! column-level elevation (m)
    real(r8) :: hsurf_g  ! gridcell-level elevation (m)
    real(r8) :: sum_lwrad_g    (1:num_gridcells) ! weighted sum of column-level lwrad
    real(r8) :: sum_wts_g      (1:num_gridcells) ! sum of weights that contribute to sum_lwrad_g
    real(r8) :: lwrad_norm_g   (1:num_gridcells) ! normalization factors
    real(r8) :: newsum_lwrad_g (1:num_gridcells) ! weighted sum of column-level lwrad after normalization

    real(r8) :: pv_g ! the water vapor pressure at grid cell (hPa)
    real(r8) :: pv_c ! the water vapor pressure at column (hPa)
    real(r8) :: emissivity_clearsky_g ! clear-sky emissivity at grid cell
    real(r8) :: emissivity_clearsky_c ! clear-sky emissivity at grid column
    real(r8) :: emissivity_allsky_g   ! all-sky emissivity at grid cell
    real(r8) :: es_g, es_c, dum1, dum2, dum3

    real(r8), parameter :: lapse_rate_longwave = 0.032_r8  ! longwave radiation lapse rate (W m-2 m-1)
    real(r8), parameter :: longwave_downscaling_limit = 0.5_r8  ! relative limit for how much longwave downscaling can be done (unitless)

    !-----------------------------------------------------------------------

    ! Initialize column forcing (needs to be done for ALL active columns)
    do g = 1, num_gridcells
       do c = begc(g), endc(g)
          forc_lwrad_c(c) = forc_lwrad_g(g)
       end do
    end do

    ! Downscale the longwave radiation, conserving energy

    ! Initialize variables related to normalization
    do g = 1, num_gridcells
       sum_lwrad_g(g) = 0._r8
       sum_wts_g(g) = 0._r8
       newsum_lwrad_g(g) = 0._r8
    end do

    ! Do the downscaling
    do g = 1, num_gridcells
       do c = begc(g), endc(g)

          hsurf_g = forc_topo_g(g)
          hsurf_c = forc_topo_c(c)

#if(defined option_longwave_adjust_I)
          ! Fiddes and Gruber, 2014, TopoSCALE v.1.0: downscaling gridded climate data in
          ! complex terrain. Geosci. Model Dev., 7, 387-405. doi:10.5194/gmd-7-387-2014.
          ! Equation (1) (2) (3); here, the empirical parameters x1 and x2 are different from
          ! Konzelmann et al. (1994) where x1 = 0.443 and x2 = 8 (optimal for measurements on the Greenland ice sheet)

          call Qsadv(forc_t_g(g)  ,forc_pbot_g(g)  ,es_g,dum1,dum2,dum3)
          call Qsadv(forc_t_c(c),forc_pbot_c(c),es_c,dum1,dum2,dum3)
          pv_g = forc_q_g(g)  *es_g/100._r8  ! (hPa)
          pv_c = forc_q_c(c)*es_c/100._r8  ! (hPa)

          emissivity_clearsky_g = 0.23_r8 + 0.43_r8*(pv_g/forc_t_g(g))**(1._r8/5.7_r8)
          emissivity_clearsky_c = 0.23_r8 + 0.43_r8*(pv_c/forc_t_c(c))**(1._r8/5.7_r8)
          emissivity_allsky_g = forc_lwrad_g(g) / (5.67e-8_r8*forc_t_g(g)**4)

          forc_lwrad_c(c) = (emissivity_clearsky_c + (emissivity_allsky_g - emissivity_clearsky_g)) &
                            * 5.67e-8_r8*forc_t_c(c)**4
#else
          ! Longwave radiation is downscaled by assuming a linear decrease in downwelling longwave radiation
          ! with increasing elevation (0.032 W m-2 m-1, limited to 0.5 - 1.5 times the gridcell mean value,
          ! then normalized to conserve gridcell total energy) (Van Tricht et al., 2016, TC) Figure 6,
          ! doi:10.5194/tc-10-2379-2016

          if (glaciers(c)) then
              forc_lwrad_c(c) = forc_lwrad_g(g) - lapse_rate_longwave * (hsurf_c-hsurf_g)

          ! Here we assume that deltaLW = (dLW/dT)*(dT/dz)*deltaz
          ! We get dLW/dT = 4*eps*sigma*T^3 = 4*LW/T from the Stefan-Boltzmann law,
          ! evaluated at the mean temp. We assume the same temperature lapse rate as above.

          else
              forc_lwrad_c(c) = forc_lwrad_g(g) &
                       - 4.0_r8 * forc_lwrad_g(g)/(0.5_r8*(forc_t_c(c)+forc_t_g(g))) &
                       * lapse_rate * (hsurf_c - hsurf_g)
          end if
#endif

          ! But ensure that we don't depart too far from the atmospheric forcing value:
          ! negative values of lwrad are certainly bad, but small positive values might
          ! also be bad. We can especially run into trouble due to the normalization: a
          ! small lwrad value in one column can lead to a big normalization factor,
          ! leading to huge lwrad values in other columns.

          forc_lwrad_c(c) = min(forc_lwrad_c(c), &
                                  forc_lwrad_g(g) * (1._r8 + longwave_downscaling_limit))
          forc_lwrad_c(c) = max(forc_lwrad_c(c), &
                                  forc_lwrad_g(g) * (1._r8 - longwave_downscaling_limit))

          ! Keep track of the gridcell-level weighted sum for later normalization.
          ! This gridcell-level weighted sum just includes points for which we do the
          ! downscaling (e.g., glc_mec points). Thus the contributing weights
          ! generally do not add to 1. So to do the normalization properly, we also
          ! need to keep track of the weights that have contributed to this sum.

          sum_lwrad_g(g) = sum_lwrad_g(g) + wt_column(c)*forc_lwrad_c(c)
          sum_wts_g(g) = sum_wts_g(g) + wt_column(c)
       end do

       ! Normalize forc_lwrad_c(c) to conserve energy
       if (sum_wts_g(g) == 0._r8) then
          lwrad_norm_g(g) = 1.0_r8
       else if (sum_lwrad_g(g) == 0._r8) then
          lwrad_norm_g(g) = 1.0_r8
       else  ! The standard case
          lwrad_norm_g(g) = forc_lwrad_g(g) / (sum_lwrad_g(g) / sum_wts_g(g))
       end if

       do c = begc(g), endc(g)
          forc_lwrad_c(c) = forc_lwrad_c(c) * lwrad_norm_g(g)
          newsum_lwrad_g(g) = newsum_lwrad_g(g) + wt_column(c)*forc_lwrad_c(c)
       end do

    end do

    ! Make sure that, after normalization, the grid cell mean is conserved
    do g = 1, num_gridcells
       if (sum_wts_g(g) > 0._r8) then
          if (abs((newsum_lwrad_g(g) / sum_wts_g(g)) - forc_lwrad_g(g)) > 1.e-8_r8) then
             write(6,*) 'g, newsum_lwrad_g, sum_wts_g, forc_lwrad_g: ', &
                  g, newsum_lwrad_g(g), sum_wts_g(g), forc_lwrad_g(g)
             call abort
          end if
       end if
    end do

  end subroutine downscale_longwave


!-----------------------------------------------------------------------
end module MOD_DownscalingForcing
