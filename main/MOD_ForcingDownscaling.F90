#include <define.h>

MODULE MOD_ForcingDownscaling

!-----------------------------------------------------------------------------
! DESCRIPTION:
! Downscaling of gridcell-level meteorological forcings into column-level forcings
! (not included wind and solar radiation)
!
! INITIAL:
! The Community Land Model version 5.0 (CLM5.0)
!
! REVISIONS:
! Updated by Yongjiu Dai, January 16, 2023
! Revised by Lu Li, January 24, 2024
!

   USE MOD_Precision
   USE MOD_Qsadv
   USE MOD_Namelist
   USE MOD_Const_Physical
   IMPLICIT NONE

   real(r8), parameter :: SHR_CONST_MWDAIR = 28.966_r8       ! molecular weight dry air [kg/kmole]
   real(r8), parameter :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
   real(r8), parameter :: SHR_CONST_AVOGAD = 6.02214e26_r8   ! Avogadro's number [molecules/kmole]
   real(r8), parameter :: SHR_CONST_BOLTZ  = 1.38065e-23_r8  ! Boltzmann's constant [J/K/molecule]
   real(r8), parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ  ! Universal gas constant [J/K/kmole]
   real(r8), parameter :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR   ! Dry air gas constant [J/K/kg]

   real(r8) :: rair = SHR_CONST_RDAIR  ! Dry air gas constant [J/K/kg]

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

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: downscale_forcings     ! Downscale atm forcing fields from gridcell to column
   PUBLIC :: downscale_forcings_1c  ! Downscale atm forcing fields from gridcell to column

! PRIVATE MEMBER FUNCTIONS:
   PRIVATE :: rhos                  ! calculate atmospheric density
   PRIVATE :: downscale_longwave    ! Downscale longwave radiation from gridcell to column
   PRIVATE :: downscale_longwave_1c ! Downscale longwave radiation from gridcell to column


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------
   
   SUBROUTINE downscale_forcings(&
                   num_gridcells,num_columns,begc,endc,glaciers,wt_column,&

                   !slp_c, asp_c, cur_c, svf_c, sf_c,&

                   forc_topo_g ,forc_t_g   ,forc_th_g  ,forc_q_g     ,forc_pbot_g ,&
                   forc_rho_g  ,forc_prc_g ,forc_prl_g ,forc_lwrad_g ,forc_hgt_grc,&
                   !forc_us_g   ,forc_vs_g  ,forc_swrad_g,&

                   forc_topo_c ,forc_t_c   ,forc_th_c  ,forc_q_c     ,forc_pbot_c ,&
                   forc_rho_c  ,forc_prc_c ,forc_prl_c ,forc_lwrad_c) 
                   !forc_swrad_c,forc_us_c   ,forc_vs_c)

!-----------------------------------------------------------------------------
! DESCRIPTION:
! Downscale atmospheric forcing fields from gridcell to column.
!
! Downscaling is done based on the difference between each land model column's elevation and
! the atmosphere's surface elevation (which is the elevation at which the atmospheric
! forcings are valid).
!
! Note that the downscaling procedure can result in changes in grid cell mean values
! compared to what was provided by the atmosphere. We conserve fluxes of mass and
! energy, but allow states such as temperature to differ.
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   ! ARGUMENTS:
   integer,  intent(in) :: num_gridcells            ! number of gridcell (element)
   integer,  intent(in) :: num_columns              ! number of column (patch)
   integer,  intent(in) :: begc (1:num_gridcells)   ! beginning column index
   integer,  intent(in) :: endc (1:num_gridcells)   ! ending column index
   logical,  intent(in) :: glaciers (1:num_columns) ! true: glacier column (itypwat = 3)
   real(r8), intent(in) :: wt_column(1:num_columns) ! weight of the column relative to the grid cell

   ! topography factor
   !real(r8), intent(in) :: forc_slp_c(1:num_columns) ! slope factor
   !real(r8), intent(in) :: forc_asp_c(1:num_columns) ! aspect factor
   !real(r8), intent(in) :: forc_cur_c(1:num_columns) ! curvature factor
   !real(r8), intent(in) :: forc_svf_c(1:num_columns) ! sky view factor
   !real(r8), intent(in) :: forc_sf_c(1:num_columns)  ! shadow factor 

   ! Gridcell-level non-downscaled fields:
   real(r8), intent(in) :: forc_topo_g  (1:num_gridcells) ! atmospheric surface height [m]
   real(r8), intent(in) :: forc_t_g     (1:num_gridcells) ! atmospheric temperature [Kelvin]
   real(r8), intent(in) :: forc_th_g    (1:num_gridcells) ! atmospheric potential temperature [Kelvin]
   real(r8), intent(in) :: forc_q_g     (1:num_gridcells) ! atmospheric specific humidity [kg/kg]
   real(r8), intent(in) :: forc_pbot_g  (1:num_gridcells) ! atmospheric pressure [Pa]
   real(r8), intent(in) :: forc_rho_g   (1:num_gridcells) ! atmospheric density [kg/m**3]
   real(r8), intent(in) :: forc_prc_g   (1:num_gridcells) ! convective precipitation in grid [mm/s]
   real(r8), intent(in) :: forc_prl_g   (1:num_gridcells) ! large-scale precipitation in grid [mm/s]
   real(r8), intent(in) :: forc_lwrad_g (1:num_gridcells) ! grid downward longwave [W/m**2]
   real(r8), intent(in) :: forc_hgt_grc (1:num_gridcells) ! atmospheric reference height [m]
   !real(r8), intent(in) :: forc_us_g    (1:num_gridcells) ! atmospheric u wind [m/s]
   !real(r8), intent(in) :: forc_vs_g    (1:num_gridcells) ! atmospheric v wind [m/s]
   !real(r8), intent(in) :: forc_swrad_g (1:num_gridcells) ! grid downward shortwave [W/m**2]

   ! Column-level downscaled fields:
   real(r8), intent(in)  :: forc_topo_c (1:num_columns) ! column surface height [m]
   real(r8), intent(out) :: forc_t_c    (1:num_columns) ! atmospheric temperature [Kelvin]
   real(r8), intent(out) :: forc_th_c   (1:num_columns) ! atmospheric potential temperature [Kelvin]
   real(r8), intent(out) :: forc_q_c    (1:num_columns) ! atmospheric specific humidity [kg/kg]
   real(r8), intent(out) :: forc_pbot_c (1:num_columns) ! atmospheric pressure [Pa]
   real(r8), intent(out) :: forc_rho_c  (1:num_columns) ! atmospheric density [kg/m**3]
   real(r8), intent(out) :: forc_prc_c  (1:num_columns) ! column convective precipitation [mm/s]
   real(r8), intent(out) :: forc_prl_c  (1:num_columns) ! column large-scale precipitation [mm/s]
   real(r8), intent(out) :: forc_lwrad_c(1:num_columns) ! column downward longwave [W/m**2]
   !real(r8), intent(out) :: forc_swrad_c(1:num_columns) ! column downward shortwave [W/m**2]
   !real(r8), intent(out) :: forc_us_c   (1:num_columns) ! column u wind [m/s]
   !real(r8), intent(out) :: forc_vs_c   (1:num_columns) ! column v wind [m/s]
    
   ! Local variables for topo downscaling:
   integer :: g,c   ! indices
   
   real(r8) :: hsurf_g, hsurf_c
   real(r8) :: Hbot, zbot
   real(r8) :: tbot_g, pbot_g, thbot_g, qbot_g, qs_g, es_g, rhos_g
   real(r8) :: tbot_c, pbot_c, thbot_c, qbot_c, qs_c, es_c, rhos_c
   real(r8) :: rhos_c_estimate, rhos_g_estimate
   real(r8) :: dum1, dum2

   real(r8) :: max_elev_c    ! the maximum column level elevation value within the grid
   real(r8) :: delta_prc_c   ! deviation of the column convective precipitation from the grid level precipitation
   real(r8) :: delta_prl_c   ! deviation of the column large-scale precipitation from the grid level precipitation

   !--------------------------------------------------------------------------

      ! Initialize column forcing (needs to be done for ALL active columns)
      DO g = 1, num_gridcells

         DO c = begc(g), endc(g)
            forc_t_c    (c) = forc_t_g    (g)
            forc_th_c   (c) = forc_th_g   (g)
            forc_q_c    (c) = forc_q_g    (g)
            forc_pbot_c (c) = forc_pbot_g (g)
            forc_rho_c  (c) = forc_rho_g  (g)
            forc_prc_c  (c) = forc_prc_g  (g)
            forc_prl_c  (c) = forc_prl_g  (g)
            forc_lwrad_c(c) = forc_lwrad_g(g)
         END DO
      END DO

      ! Downscale forc_t, forc_th, forc_q, forc_pbot, and forc_rho to columns.
      DO g = 1, num_gridcells
         hsurf_g = forc_topo_g(g)             ! gridcell sfc elevation
         tbot_g  = forc_t_g(g)                ! atm sfc temp
         thbot_g = forc_th_g(g)               ! atm sfc pot temp
         qbot_g  = forc_q_g(g)                ! atm sfc spec humid
         pbot_g  = forc_pbot_g(g)             ! atm sfc pressure
         rhos_g  = forc_rho_g(g)              ! atm density
         zbot    = forc_hgt_grc(g)            ! atm ref height

         max_elev_c = maxval(forc_topo_c(begc(g):endc(g))) ! maximum column level elevation value within the grid
         
         !real(r8) :: wd(1:num_columns)            ! wind direction (rad)
         !real(r8) :: slp_wd(1:num_columns)        ! the slope in the direction of wind
         !real(r8) :: norm_slp_wd   ! normalize the slope in the direction of wind
         !real(r8) :: norm_cur_c    ! normalize curvature factor
         
         ! wind adjust factor
         !wd = atan(forc_vs_g(g)/forc_us_g(g)) ! cal wind direction (rad)         
         !slp_wd = forc_slp_c(begc(g):endc(g))*cos(wd-forc_asp_c(begc(g):endc(g))) ! the slope in the direction of wind
         !norm_slp_wd = slp_wd/(2*maxval(slp_wd)) ! normalize the slope in the direction of wind
         !norm_cur_c = forc_cur_c(begc(g):endc(g))/(2*maxval(forc_cur_c(begc(g):endc(g)))) ! normalize curvature factor

         DO c = begc(g), endc(g)
            ! This is a simple downscaling procedure
            ! Note that forc_hgt, forc_u, forc_v and solar radiation are not downscaled.

            !asp_c = forc_asp_c(c)
            !cur_c = forc_cur_c(c)

            hsurf_c = forc_topo_c(c)                        ! column sfc elevation
            tbot_c  = tbot_g-lapse_rate*(hsurf_c-hsurf_g)   ! adjust temp for column
            Hbot    = rair*0.5_r8*(tbot_g+tbot_c)/grav      ! scale ht at avg temp
            pbot_c  = pbot_g*exp(-(hsurf_c-hsurf_g)/Hbot)   ! adjust press for column

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
            thbot_c = thbot_g + (tbot_c - tbot_g)*exp((zbot/Hbot)*(rair/cpair)) ! adjust pot temp for column

            CALL Qsadv(tbot_g,pbot_g,es_g,dum1,qs_g,dum2) ! es, qs for gridcell
            CALL Qsadv(tbot_c,pbot_c,es_c,dum1,qs_c,dum2) ! es, qs for column
            qbot_c = qbot_g*(qs_c/qs_g) ! adjust q for column

            rhos_c_estimate = rhos(qbot=qbot_c, pbot=pbot_c, tbot=tbot_c)
            rhos_g_estimate = rhos(qbot=qbot_g, pbot=pbot_g, tbot=tbot_g)
            rhos_c = rhos_g * (rhos_c_estimate / rhos_g_estimate) ! adjust density for column

            forc_t_c(c)    = tbot_c
            forc_th_c(c)   = thbot_c
            forc_q_c(c)    = qbot_c
            forc_pbot_c(c) = pbot_c
            forc_rho_c(c)  = rhos_c

            ! adjust precipitation 
            IF (trim(DEF_DS_precipitation_adjust_scheme) == 'I') THEN
               ! Tesfa et al, 2020: Exploring Topography-Based Methods for Downscaling
               ! Subgrid Precipitation for Use in Earth System Models. Equation (5)
               ! https://doi.org/ 10.1029/2019JD031456

               delta_prc_c = forc_prc_g(g) * (forc_topo_c(c) - forc_topo_g(g)) / max_elev_c
               forc_prc_c(c) = forc_prc_g(g) + delta_prc_c   ! convective precipitation [mm/s]

               delta_prl_c = forc_prl_g(g) * (forc_topo_c(c) - forc_topo_g(g)) / max_elev_c
               forc_prl_c(c) = forc_prl_g(g) + delta_prl_c   ! large scale precipitation [mm/s]

            ELSEIF (trim(DEF_DS_precipitation_adjust_scheme) == 'II') THEN
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
             
            ELSEIF (trim(DEF_DS_precipitation_adjust_scheme) == 'III') THEN 
               ! Mei, Y., Maggioni, V., Houser, P., Xue, Y., & Rouf, T. (2020). A nonparametric statistical 
               ! technique for spatial downscaling of precipitation over High Mountain Asia. Water Resources Research, 
               ! 56, e2020WR027472. https://doi.org/ 10.1029/2020WR027472
               ! Change Random forest model to AutoML model.
               !TODO: Lu Li; Need to done after all other forcings are downscaled
            END IF

            IF (forc_prl_c(c) < 0) THEN
               write(*,*) 'negative prl', forc_prl_g(g), max_elev_c, forc_topo_c(c), forc_topo_g(g)
            END IF

            IF (forc_prc_c(c) < 0) THEN
               write(*,*) 'negative prc', forc_prc_g(g), max_elev_c, forc_topo_c(c), forc_topo_g(g)
            END IF

            forc_prc_c(c) = max(forc_prc_c(c), 0.)
            forc_prl_c(c) = max(forc_prl_c(c), 0.)

         END DO
      END DO

      CALL downscale_longwave(num_gridcells, num_columns, begc, endc, glaciers, wt_column, &
                   forc_topo_g, forc_t_g, forc_q_g, forc_pbot_g, forc_lwrad_g, &
                   forc_topo_c, forc_t_c, forc_q_c, forc_pbot_c, forc_lwrad_c)

   END SUBROUTINE downscale_forcings



!-----------------------------------------------------------------------------

   PURE FUNCTION rhos(qbot, pbot, tbot)
   
!-----------------------------------------------------------------------------
! DESCRIPTION:
! Compute atmospheric density (kg/m**3)
!-----------------------------------------------------------------------------

   IMPLICIT NONE
   
   ! ARGUMENTS:
   real(r8) :: rhos  ! function result: atmospheric density (kg/m**3)
   real(r8), intent(in) :: qbot  ! atmospheric specific humidity (kg/kg)
   real(r8), intent(in) :: pbot  ! atmospheric pressure (Pa)
   real(r8), intent(in) :: tbot  ! atmospheric temperature (K)
   
   ! LOCAL VARIABLES:
   real(r8) :: egcm
   real(r8) :: wv_to_dair_weight_ratio  ! ratio of molecular weight of water vapor to that of dry air [-]

   !--------------------------------------------------------------------------
      wv_to_dair_weight_ratio = SHR_CONST_MWWV/SHR_CONST_MWDAIR

      egcm = qbot*pbot / (wv_to_dair_weight_ratio + (1._r8 - wv_to_dair_weight_ratio)*qbot)
      rhos = (pbot - (1._r8 - wv_to_dair_weight_ratio)*egcm) / (rair*tbot)

   END FUNCTION rhos



!-----------------------------------------------------------------------------

   SUBROUTINE downscale_longwave(&
      num_gridcells, num_columns, begc, endc, glaciers, wt_column, &
      forc_topo_g, forc_t_g, forc_q_g, forc_pbot_g, forc_lwrad_g, &
      forc_topo_c, forc_t_c, forc_q_c, forc_pbot_c, forc_lwrad_c)

!-----------------------------------------------------------------------------
! DESCRIPTION:
! Downscale longwave radiation from gridcell to column
! Must be done AFTER temperature downscaling
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   ! ARGUMENTS:
   integer,  intent(in) :: num_gridcells  ! number of gridcell
   integer,  intent(in) :: num_columns    ! number of column
   integer,  intent(in) :: begc (1:num_gridcells) ! beginning column index
   integer,  intent(in) :: endc (1:num_gridcells) ! ending column index
   logical,  intent(in) :: glaciers  (1:num_columns) ! true: glacier column
   real(r8), intent(in) :: wt_column (1:num_columns) ! weight of the column relative to the grid cell

   real(r8), intent(in) :: forc_topo_g  (1:num_gridcells) ! atmospheric surface height (m)
   real(r8), intent(in) :: forc_t_g     (1:num_gridcells) ! atmospheric temperature [Kelvin]
   real(r8), intent(in) :: forc_q_g     (1:num_gridcells) ! atmospheric specific humidity [kg/kg]
   real(r8), intent(in) :: forc_pbot_g  (1:num_gridcells) ! atmospheric pressure [Pa]
   real(r8), intent(in) :: forc_lwrad_g (1:num_gridcells) ! downward longwave (W/m**2)
   real(r8), intent(in) :: forc_topo_c  (1:num_columns) ! column surface height (m)
   real(r8), intent(in) :: forc_t_c     (1:num_columns) ! atmospheric temperature [Kelvin]
   real(r8), intent(in) :: forc_q_c     (1:num_columns) ! atmospheric specific humidity [kg/kg]
   real(r8), intent(in) :: forc_pbot_c  (1:num_columns) ! atmospheric pressure [Pa]
   real(r8), intent(out) :: forc_lwrad_c(1:num_columns) ! downward longwave (W/m**2)

   ! LOCAL VARIABLES:
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

   !--------------------------------------------------------------------------

      ! Initialize column forcing (needs to be done for ALL active columns)
      DO g = 1, num_gridcells
         DO c = begc(g), endc(g)
            forc_lwrad_c(c) = forc_lwrad_g(g)
         END DO
      END DO

      ! Downscale the longwave radiation, conserving energy

      ! Initialize variables related to normalization
      DO g = 1, num_gridcells
         sum_lwrad_g(g) = 0._r8
         sum_wts_g(g) = 0._r8
         newsum_lwrad_g(g) = 0._r8
      END DO

      ! Do the downscaling
      DO g = 1, num_gridcells
         DO c = begc(g), endc(g)

            hsurf_g = forc_topo_g(g)
            hsurf_c = forc_topo_c(c)

            IF (trim(DEF_DS_longwave_adjust_scheme) == 'I') THEN
               ! Fiddes and Gruber, 2014, TopoSCALE v.1.0: downscaling gridded climate data in
               ! complex terrain. Geosci. Model Dev., 7, 387-405. doi:10.5194/gmd-7-387-2014.
               ! Equation (1) (2) (3); here, the empirical parameters x1 and x2 are different from
               ! Konzelmann et al. (1994) where x1 = 0.443 and x2 = 8 (optimal for measurements on the Greenland ice sheet)

               CALL Qsadv(forc_t_g(g)  ,forc_pbot_g(g)  ,es_g,dum1,dum2,dum3)
               CALL Qsadv(forc_t_c(c),forc_pbot_c(c),es_c,dum1,dum2,dum3)
               pv_g = forc_q_g(g)  *es_g/100._r8  ! (hPa)
               pv_c = forc_q_c(c)*es_c/100._r8  ! (hPa)

               emissivity_clearsky_g = 0.23_r8 + 0.43_r8*(pv_g/forc_t_g(g))**(1._r8/5.7_r8)
               emissivity_clearsky_c = 0.23_r8 + 0.43_r8*(pv_c/forc_t_c(c))**(1._r8/5.7_r8)
               emissivity_allsky_g = forc_lwrad_g(g) / (5.67e-8_r8*forc_t_g(g)**4)

               forc_lwrad_c(c) = (emissivity_clearsky_c + (emissivity_allsky_g - emissivity_clearsky_g)) &
                  * 5.67e-8_r8*forc_t_c(c)**4
            ELSE
               ! Longwave radiation is downscaled by assuming a linear decrease in downwelling longwave radiation
               ! with increasing elevation (0.032 W m-2 m-1, limited to 0.5 - 1.5 times the gridcell mean value,
               ! then normalized to conserve gridcell total energy) (Van Tricht et al., 2016, TC) Figure 6,
               ! doi:10.5194/tc-10-2379-2016

               IF (glaciers(c)) THEN
                  forc_lwrad_c(c) = forc_lwrad_g(g) - lapse_rate_longwave * (hsurf_c-hsurf_g)

                  ! Here we assume that deltaLW = (dLW/dT)*(dT/dz)*deltaz
                  ! We get dLW/dT = 4*eps*sigma*T^3 = 4*LW/T from the Stefan-Boltzmann law,
                  ! evaluated at the mean temp. We assume the same temperature lapse rate as above.

               ELSE
                  forc_lwrad_c(c) = forc_lwrad_g(g) &
                     - 4.0_r8 * forc_lwrad_g(g)/(0.5_r8*(forc_t_c(c)+forc_t_g(g))) &
                     * lapse_rate * (hsurf_c - hsurf_g)
               END IF
            END IF

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
         END DO

         ! Normalize forc_lwrad_c(c) to conserve energy
         IF (sum_wts_g(g) == 0._r8) THEN
            lwrad_norm_g(g) = 1.0_r8
         ELSE IF (sum_lwrad_g(g) == 0._r8) THEN
            lwrad_norm_g(g) = 1.0_r8
         ELSE  ! The standard case
            lwrad_norm_g(g) = forc_lwrad_g(g) / (sum_lwrad_g(g) / sum_wts_g(g))
         END IF

         DO c = begc(g), endc(g)
            forc_lwrad_c(c) = forc_lwrad_c(c) * lwrad_norm_g(g)
            newsum_lwrad_g(g) = newsum_lwrad_g(g) + wt_column(c)*forc_lwrad_c(c)
         END DO

      END DO

      ! Make sure that, after normalization, the grid cell mean is conserved
      DO g = 1, num_gridcells
         IF (sum_wts_g(g) > 0._r8) THEN
            IF (abs((newsum_lwrad_g(g) / sum_wts_g(g)) - forc_lwrad_g(g)) > 1.e-8_r8) THEN
               write(6,*) 'g, newsum_lwrad_g, sum_wts_g, forc_lwrad_g: ', &
                     g, newsum_lwrad_g(g), sum_wts_g(g), forc_lwrad_g(g)
               CALL abort
            END IF
         END IF
      END DO

   END SUBROUTINE downscale_longwave



   !-----------------------------------------------------------------------------
   
   SUBROUTINE downscale_forcings_1c (&
                   glaciers, &

                   !slp_c, asp_c, cur_c, svf_c, sf_c,&

                   forc_topo_g ,forc_maxelv_g ,forc_t_g   ,forc_th_g  ,forc_q_g     ,&
                   forc_pbot_g ,forc_rho_g    ,forc_prc_g ,forc_prl_g ,forc_lwrad_g ,&
                   forc_hgt_grc,&
                   !forc_us_g   ,forc_vs_g  ,forc_swrad_g,&

                   forc_topo_c ,forc_t_c   ,forc_th_c  ,forc_q_c     ,forc_pbot_c ,&
                   forc_rho_c  ,forc_prc_c ,forc_prl_c ,forc_lwrad_c) 
                   !forc_swrad_c,forc_us_c   ,forc_vs_c)

!-----------------------------------------------------------------------------
! DESCRIPTION:
! Downscale atmospheric forcing fields from gridcell to column.
!
! Downscaling is done based on the difference between each land model column's elevation and
! the atmosphere's surface elevation (which is the elevation at which the atmospheric
! forcings are valid).
!
! Note that the downscaling procedure can result in changes in grid cell mean values
! compared to what was provided by the atmosphere. We conserve fluxes of mass and
! energy, but allow states such as temperature to differ.
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   ! ARGUMENTS:
   logical,  intent(in) :: glaciers ! true: glacier column (itypwat = 3)

   ! Gridcell-level non-downscaled fields:
   real(r8), intent(in) :: forc_topo_g   ! atmospheric surface height [m]
   real(r8), intent(in) :: forc_maxelv_g ! max atmospheric surface height [m]
   real(r8), intent(in) :: forc_t_g      ! atmospheric temperature [Kelvin]
   real(r8), intent(in) :: forc_th_g     ! atmospheric potential temperature [Kelvin]
   real(r8), intent(in) :: forc_q_g      ! atmospheric specific humidity [kg/kg]
   real(r8), intent(in) :: forc_pbot_g   ! atmospheric pressure [Pa]
   real(r8), intent(in) :: forc_rho_g    ! atmospheric density [kg/m**3]
   real(r8), intent(in) :: forc_prc_g    ! convective precipitation in grid [mm/s]
   real(r8), intent(in) :: forc_prl_g    ! large-scale precipitation in grid [mm/s]
   real(r8), intent(in) :: forc_lwrad_g  ! grid downward longwave [W/m**2]
   real(r8), intent(in) :: forc_hgt_grc  ! atmospheric reference height [m]

   ! Column-level downscaled fields:
   real(r8), intent(in)  :: forc_topo_c  ! column surface height [m]
   real(r8), intent(out) :: forc_t_c     ! atmospheric temperature [Kelvin]
   real(r8), intent(out) :: forc_th_c    ! atmospheric potential temperature [Kelvin]
   real(r8), intent(out) :: forc_q_c     ! atmospheric specific humidity [kg/kg]
   real(r8), intent(out) :: forc_pbot_c  ! atmospheric pressure [Pa]
   real(r8), intent(out) :: forc_rho_c   ! atmospheric density [kg/m**3]
   real(r8), intent(out) :: forc_prc_c   ! column convective precipitation [mm/s]
   real(r8), intent(out) :: forc_prl_c   ! column large-scale precipitation [mm/s]
   real(r8), intent(out) :: forc_lwrad_c ! column downward longwave [W/m**2]
    
   ! Local variables for topo downscaling:
   
   real(r8) :: hsurf_g, hsurf_c
   real(r8) :: Hbot, zbot
   real(r8) :: tbot_g, pbot_g, thbot_g, qbot_g, qs_g, es_g, rhos_g
   real(r8) :: tbot_c, pbot_c, thbot_c, qbot_c, qs_c, es_c, rhos_c
   real(r8) :: rhos_c_estimate, rhos_g_estimate
   real(r8) :: dum1, dum2

   real(r8) :: max_elev_c    ! the maximum column level elevation value within the grid
   real(r8) :: delta_prc_c   ! deviation of the column convective precipitation from the grid level precipitation
   real(r8) :: delta_prl_c   ! deviation of the column large-scale precipitation from the grid level precipitation


      ! ! Downscale forc_t, forc_th, forc_q, forc_pbot, and forc_rho to columns.
      hsurf_g = forc_topo_g             ! gridcell sfc elevation
      tbot_g  = forc_t_g                ! atm sfc temp
      thbot_g = forc_th_g               ! atm sfc pot temp
      qbot_g  = forc_q_g                ! atm sfc spec humid
      pbot_g  = forc_pbot_g             ! atm sfc pressure
      rhos_g  = forc_rho_g              ! atm density
      zbot    = forc_hgt_grc            ! atm ref height

      ! This is a simple downscaling procedure
      ! Note that forc_hgt, forc_u, forc_v and solar radiation are not downscaled.

      !asp_c = forc_asp_c(c)
      !cur_c = forc_cur_c(c)

      hsurf_c = forc_topo_c                        ! column sfc elevation
      tbot_c  = tbot_g-lapse_rate*(hsurf_c-hsurf_g)   ! adjust temp for column
      Hbot    = rair*0.5_r8*(tbot_g+tbot_c)/grav      ! scale ht at avg temp
      pbot_c  = pbot_g*exp(-(hsurf_c-hsurf_g)/Hbot)   ! adjust press for column

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
      thbot_c = thbot_g + (tbot_c - tbot_g)*exp((zbot/Hbot)*(rair/cpair)) ! adjust pot temp for column

      CALL Qsadv(tbot_g,pbot_g,es_g,dum1,qs_g,dum2) ! es, qs for gridcell
      CALL Qsadv(tbot_c,pbot_c,es_c,dum1,qs_c,dum2) ! es, qs for column
      qbot_c = qbot_g*(qs_c/qs_g) ! adjust q for column

      rhos_c_estimate = rhos(qbot=qbot_c, pbot=pbot_c, tbot=tbot_c)
      rhos_g_estimate = rhos(qbot=qbot_g, pbot=pbot_g, tbot=tbot_g)
      rhos_c = rhos_g * (rhos_c_estimate / rhos_g_estimate) ! adjust density for column

      forc_t_c    = tbot_c
      forc_th_c   = thbot_c
      forc_q_c    = qbot_c
      forc_pbot_c = pbot_c
      forc_rho_c  = rhos_c

      ! adjust precipitation 
      IF (trim(DEF_DS_precipitation_adjust_scheme) == 'I') THEN
         ! Tesfa et al, 2020: Exploring Topography-Based Methods for Downscaling
         ! Subgrid Precipitation for Use in Earth System Models. Equation (5)
         ! https://doi.org/ 10.1029/2019JD031456

         delta_prc_c = forc_prc_g * (forc_topo_c - forc_topo_g) / forc_maxelv_g
         forc_prc_c  = forc_prc_g + delta_prc_c   ! convective precipitation [mm/s]

         delta_prl_c = forc_prl_g * (forc_topo_c - forc_topo_g) / forc_maxelv_g
         forc_prl_c  = forc_prl_g + delta_prl_c   ! large scale precipitation [mm/s]

      ELSEIF (trim(DEF_DS_precipitation_adjust_scheme) == 'II') THEN
         ! Liston, G. E. and Elder, K.: A meteorological distribution system
         ! for high-resolution terrestrial modeling (MicroMet), J. Hydrometeorol., 7, 217-234, 2006.
         ! Equation (33) and Table 1: chi range from January to December:
         ! [0.35,0.35,0.35,0.30,0.25,0.20,0.20,0.20,0.20,0.25,0.30,0.35] (1/m)

         delta_prc_c = forc_prc_g * 2.0*0.27e-3*(forc_topo_c - forc_topo_g) &
            /(1.0 - 0.27e-3*(forc_topo_c - forc_topo_g))
         forc_prc_c = forc_prc_g + delta_prc_c   ! large scale precipitation [mm/s]

         delta_prl_c = forc_prl_g * 2.0*0.27e-3*(forc_topo_c - forc_topo_g) &
            /(1.0 - 0.27e-3*(forc_topo_c - forc_topo_g))
         forc_prl_c = forc_prl_g + delta_prl_c   ! large scale precipitation [mm/s]

      ELSEIF (trim(DEF_DS_precipitation_adjust_scheme) == 'III') THEN 
         ! Mei, Y., Maggioni, V., Houser, P., Xue, Y., & Rouf, T. (2020). A nonparametric statistical 
         ! technique for spatial downscaling of precipitation over High Mountain Asia. Water Resources Research, 
         ! 56, e2020WR027472. https://doi.org/ 10.1029/2020WR027472
         ! Change Random forest model to AutoML model.
         !TODO: Lu Li; Need to done after all other forcings are downscaled
      END IF

      IF (forc_prl_c < 0) THEN
         write(*,*) 'negative prl', forc_prl_g, forc_maxelv_g, forc_topo_c, forc_topo_g
         forc_prl_c = 0.
      END IF

      IF (forc_prc_c < 0) THEN
         write(*,*) 'negative prc', forc_prc_g, forc_maxelv_g, forc_topo_c, forc_topo_g
         forc_prc_c = 0.
      END IF


      CALL downscale_longwave_1c (glaciers, &
                   forc_topo_g, forc_t_g, forc_q_g, forc_pbot_g, forc_lwrad_g, &
                   forc_topo_c, forc_t_c, forc_q_c, forc_pbot_c, forc_lwrad_c)

   END SUBROUTINE downscale_forcings_1c



!-----------------------------------------------------------------------------

   SUBROUTINE downscale_longwave_1c (glaciers, &
      forc_topo_g, forc_t_g, forc_q_g, forc_pbot_g, forc_lwrad_g, &
      forc_topo_c, forc_t_c, forc_q_c, forc_pbot_c, forc_lwrad_c)

!-----------------------------------------------------------------------------
! DESCRIPTION:
! Downscale longwave radiation from gridcell to column
! Must be done AFTER temperature downscaling
!-----------------------------------------------------------------------------

   IMPLICIT NONE

   ! ARGUMENTS:
   logical,  intent(in) :: glaciers  ! true: glacier column

   real(r8), intent(in) :: forc_topo_g   ! atmospheric surface height (m)
   real(r8), intent(in) :: forc_t_g      ! atmospheric temperature [Kelvin]
   real(r8), intent(in) :: forc_q_g      ! atmospheric specific humidity [kg/kg]
   real(r8), intent(in) :: forc_pbot_g   ! atmospheric pressure [Pa]
   real(r8), intent(in) :: forc_lwrad_g  ! downward longwave (W/m**2)
   real(r8), intent(in) :: forc_topo_c   ! column surface height (m)
   real(r8), intent(in) :: forc_t_c      ! atmospheric temperature [Kelvin]
   real(r8), intent(in) :: forc_q_c      ! atmospheric specific humidity [kg/kg]
   real(r8), intent(in) :: forc_pbot_c   ! atmospheric pressure [Pa]
   real(r8), intent(out) :: forc_lwrad_c ! downward longwave (W/m**2)

   ! LOCAL VARIABLES:
   real(r8) :: hsurf_c  ! column-level elevation (m)
   real(r8) :: hsurf_g  ! gridcell-level elevation (m)

   real(r8) :: pv_g ! the water vapor pressure at grid cell (hPa)
   real(r8) :: pv_c ! the water vapor pressure at column (hPa)
   real(r8) :: emissivity_clearsky_g ! clear-sky emissivity at grid cell
   real(r8) :: emissivity_clearsky_c ! clear-sky emissivity at grid column
   real(r8) :: emissivity_allsky_g   ! all-sky emissivity at grid cell
   real(r8) :: es_g, es_c, dum1, dum2, dum3

   real(r8), parameter :: lapse_rate_longwave = 0.032_r8  ! longwave radiation lapse rate (W m-2 m-1)
   real(r8), parameter :: longwave_downscaling_limit = 0.5_r8  ! relative limit for how much longwave downscaling can be done (unitless)

   !--------------------------------------------------------------------------

      ! Initialize column forcing (needs to be done for ALL active columns)
      forc_lwrad_c = forc_lwrad_g

      ! Do the downscaling

      hsurf_g = forc_topo_g
      hsurf_c = forc_topo_c

      IF (trim(DEF_DS_longwave_adjust_scheme) == 'I') THEN
         ! Fiddes and Gruber, 2014, TopoSCALE v.1.0: downscaling gridded climate data in
         ! complex terrain. Geosci. Model Dev., 7, 387-405. doi:10.5194/gmd-7-387-2014.
         ! Equation (1) (2) (3); here, the empirical parameters x1 and x2 are different from
         ! Konzelmann et al. (1994) where x1 = 0.443 and x2 = 8 (optimal for measurements on the Greenland ice sheet)

         CALL Qsadv(forc_t_g, forc_pbot_g, es_g,dum1,dum2,dum3)
         CALL Qsadv(forc_t_c, forc_pbot_c, es_c,dum1,dum2,dum3)
         pv_g = forc_q_g*es_g/100._r8  ! (hPa)
         pv_c = forc_q_c*es_c/100._r8  ! (hPa)

         emissivity_clearsky_g = 0.23_r8 + 0.43_r8*(pv_g/forc_t_g)**(1._r8/5.7_r8)
         emissivity_clearsky_c = 0.23_r8 + 0.43_r8*(pv_c/forc_t_c)**(1._r8/5.7_r8)
         emissivity_allsky_g = forc_lwrad_g / (5.67e-8_r8*forc_t_g**4)

         forc_lwrad_c = &
            (emissivity_clearsky_c + (emissivity_allsky_g - emissivity_clearsky_g)) &
            * 5.67e-8_r8*forc_t_c**4
      ELSE
         ! Longwave radiation is downscaled by assuming a linear decrease in downwelling longwave radiation
         ! with increasing elevation (0.032 W m-2 m-1, limited to 0.5 - 1.5 times the gridcell mean value,
         ! then normalized to conserve gridcell total energy) (Van Tricht et al., 2016, TC) Figure 6,
         ! doi:10.5194/tc-10-2379-2016

         IF (glaciers) THEN
            forc_lwrad_c = forc_lwrad_g - lapse_rate_longwave * (hsurf_c-hsurf_g)

            ! Here we assume that deltaLW = (dLW/dT)*(dT/dz)*deltaz
            ! We get dLW/dT = 4*eps*sigma*T^3 = 4*LW/T from the Stefan-Boltzmann law,
            ! evaluated at the mean temp. We assume the same temperature lapse rate as above.

         ELSE
            forc_lwrad_c = forc_lwrad_g &
               - 4.0_r8 * forc_lwrad_g/(0.5_r8*(forc_t_c+forc_t_g)) &
               * lapse_rate * (hsurf_c - hsurf_g)
         END IF
      END IF

      ! But ensure that we don't depart too far from the atmospheric forcing value:
      ! negative values of lwrad are certainly bad, but small positive values might
      ! also be bad. We can especially run into trouble due to the normalization: a
      ! small lwrad value in one column can lead to a big normalization factor,
      ! leading to huge lwrad values in other columns.

      forc_lwrad_c = min(forc_lwrad_c, forc_lwrad_g * (1._r8 + longwave_downscaling_limit))
      forc_lwrad_c = max(forc_lwrad_c, forc_lwrad_g * (1._r8 - longwave_downscaling_limit))


   END SUBROUTINE downscale_longwave_1c

END MODULE MOD_ForcingDownscaling
