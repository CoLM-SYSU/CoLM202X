#include <define.h>
MODULE  MOD_Lake_Const

   USE MOD_Precision

   IMPLICIT NONE

!--------------------------------------------------------------
PUBLIC
!--------------------------------------------------------------

   ! lake model 
   integer,  parameter         :: &
      COLML        =  1          ,&! CoLM-Lake
      FLAKE        =  2          ,&! FLake
      SIMSTRAT     =  3          ,&! Simstrat
      XOML         =  4            ! XOML

   ! partition between shallow and deep lake
   real(r8), parameter         :: &
      ShallowDeepPartition = 20.   ! partition between shallow and deep lake

   ! lake type
   integer , parameter :: &
      SHALLOW      =  1          ,&! shallow lake
      DEEP         =  2            ! deep lake

   ! lake volume scheme
   integer , parameter         :: &
      TRAP_PRISM =  1            ,&! trapezoidal prism (GLM scheme)
      TRUNC_CONE =  2              ! truncated cone 


!-----------------  Layer Scheme  -----------------
!****************************************************************************************************
   integer,  parameter         :: &
      COLMLLayer   =  1          ,&! CoLM-Lake Layer scheme
      CSSPLLayer   =  2          ,&! CSSP-Lake Layer scheme
      EqualLayer   =  3          ,&! Equal Layer scheme
      XMlogLayer   =  4          ,&! MIN XU, log layer scheme
      WMlogLayer   =  5            ! WMEJ, log Layer scheme

   real(r8), parameter, DIMENSION(11,1:2) :: & !- CSSP divides lakes into 10 layers
   DZ_CSSPL= reshape(          &! layer thickness for (shallow,deep) lake [m]
            (/ 0.250, 0.500, 0.750, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
               1.000, 2.000, 3.000, 4.000, 5.000, 7.000, 7.000, 7.000, 7.000, 7.000, 7.000/)&
                     , (/11,2/))

   real(r8), parameter, DIMENSION(10,1:2) :: & !- CoLM divides lakes into 10 layers
      DZ_COLML= reshape(       &! layer thickness for (shallow,deep) lake [m]
            (/ 0.100, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
               0.100, 1.000, 2.000, 3.000, 4.000, 5.000, 7.000, 7.000, 10.450, 10.450/)&
                     , (/10,2/))

!****************************************************************************************************


!********************************************************
!               --- CoLML Parameter ---               
!********************************************************
   ! constants parameters
   real(r8), parameter :: emisi             = 0.97             ! surface emissivity for ice or snow [-]
   real(r8), parameter :: emisw             = 0.97             ! surface emissivity for water [-]
   real(r8), parameter :: emiss             = 0.96             ! surface emissivity for soil  [-]
   real(r8), parameter :: beta_wc           = 1.0              ! coefficient of convective velocity [-]
   real(r8), parameter :: betavis           = 0.0              ! The fraction of the visible (e.g. vis not nir from atm) sunlight
                                                               ! absorbed in ~1 m of water (the surface layer za_lake).
                                                               ! This is roughly the fraction over 700 nm but may depend on the details
                                                               ! of atmospheric radiative transfer.
                                                               ! As long as NIR = 700 nm and up, this can be zero.
   real(r8), parameter :: kva0              = 1.51e-5          ! kinematic viscosity of air (m^2/s) at 20C and 1.013e5 Pa
   real(r8), parameter :: prn               = 0.713            ! Prandtl # for air at neutral stability
   real(r8), parameter :: sch               = 0.66             ! Schmidt # for water in air at neutral stability
   real(r8), parameter :: fcrit             = 22.              ! critical dimensionless fetch for Charnock parameter (Vickers & Mahrt 1997)
                                                               ! but converted to USE u instead of u* (Subin et al. 2011)
   real(r8), parameter :: z0sice            = 0.0400           ! roughness length for sea ice [m]

   ! For calculating prognostic roughness length
   real(r8), parameter :: cur0              = 0.01             ! min. Charnock parameter
   real(r8), parameter :: cus               = 0.1              ! empirical constant for roughness under smooth flow
   real(r8), parameter :: curm              = 0.1              ! maximum Charnock parameter
   real(r8), parameter :: beta1             = 1.               ! coefficient of convective velocity (in computing W_*) [-]
   real(r8), parameter :: thk_bedrock       = 3.               ! xum refer to CLM4 user's guide

   real(r8), parameter :: one3rd            = 1./3.            ! one third
   real(r8), parameter :: PI                = 4*atan(1.)       ! pi value
   real(r8), parameter :: d2r               = 1.745329251994330e-2_r8         ! conversion from degree to radian
   real(r8), parameter :: densnow           = 250.             ! initial snow density [kg/m3]
   real(r8), parameter :: denice            = 917.             ! density of ice [kg/m3]
   real(r8), parameter :: denh2o            = 1000.            ! density of liquid water [kg/m3]
   real(r8), parameter :: cpliq             = 4188.            ! Specific heat of water [J/kg-K]
   real(r8), parameter :: cpice             = 2117.27          ! Specific heat of ice [J/kg-K]
   real(r8), parameter :: cpair             = 1004.64          ! e.g: 1005.70, specific heat of dry air [J/kg/K] 
   real(r8), parameter :: rair              = 287.04           ! gas constant for dry air [J/kg/K]
   real(r8), parameter :: rovcp             = rair/cpair       !
   real(r8), parameter :: hfus              = 0.3336e6         ! latent heat of fusion for ice [J/kg]
   real(r8), parameter :: hvap              = 2.5104e6         ! latent heat of evap for water [J/kg]
   real(r8), parameter :: hsub              = 2.8440e6         ! latent heat of sublimation [J/kg]
   real(r8), parameter :: tkair             = 0.023            ! thermal conductivity of air [W/m/k]
   real(r8), parameter :: tkice             = 2.290            ! thermal conductivity of ice [W/m/k]
   real(r8), parameter :: tkwat             = 0.6              ! thermal conductivity of water [W/m/k]
   real(r8), parameter :: tfrz              = 273.16           ! freezing temperature [K]
   real(r8), parameter :: rgas              = 287.04           ! gas constant for dry air [J/kg/K]
   real(r8), parameter :: roverg            = 4.71047e4        ! rw/g=(8.3144/0.018)/(9.80616)*1000. mm/K
   real(r8), parameter :: rwat              = 461.296          ! gas constant for water vapor [J/(kg K)]
   real(r8), parameter :: grav              = 9.80616          ! gravity constant [m/s2]
   real(r8), parameter :: vonkar            = 0.4              ! von Karman constant [-]
   real(r8), parameter :: stefnc            = 5.67e-8          ! Stefan-Boltzmann constant  [W/m2/K4]
   real(r8), parameter :: rd                = 287.0423         ! gas constant for dry air[J/kg/K]
   real(r8), parameter :: ssi               = 0.033            ! Irreducible water saturation of snow
   real(r8), parameter :: wimp              = 0.05             ! Water impremeable IF porosity less than wimp
   real(r8), parameter :: twmax             = 3.98             ! temperature for maximum water denisty [0C]
   real(r8), parameter :: tdmax             = twmax + tfrz     ! kevin 

   real(r8), parameter :: depthcrit         = 25.              ! (m) Depth beneath which is subject to enhance mixing.
                                                               ! See discussion in Subin et al. 2011
   real(r8), parameter :: n2min             = 7.5e-5           ! (s^-2) (yields diffusivity about 6 times km) ! Fang & Stefan 1996
   real(r8), parameter :: p0                = 1.               ! neutral value of turbulent prandtl number
   real(r8), parameter :: cnfac             = 0.5              ! Crank Nicholson factor between 0 and 1
   real(r8), parameter :: lsadz             = 0.03             ! thickness avoid numeric unstabilty [m]
   real(r8), parameter :: Visc              = 0.00000114       ! molecular diffusivity of momentum [m2/s]

!********************************************************
!               --- FLake Parameter ---               
!********************************************************
   ! constants parameters
   real(r8), parameter :: omega             = 7.292e-5         ! earth self rotation velocity [s-1]
   real(r8), parameter :: c_cbl_1           = 0.17             ! Constant in the CBL entrainment equation
   real(r8), parameter :: c_cbl_2           = 1.               ! Constant in the CBL entrainment equation
   real(r8), parameter :: c_sbl_ZM_n        = 0.5              ! Constant in the ZM1996 equation for the equilibrium SBL depth
   real(r8), parameter :: c_sbl_ZM_s        = 10.              ! Constant in the ZM1996 equation for the equilibrium SBL depth
   real(r8), parameter :: c_sbl_ZM_i        = 20.              ! Constant in the ZM1996 equation for the equilibrium SBL depth
   real(r8), parameter :: c_relax_h         = 0.030            ! Constant in the relaxation equation for the SBL depth
   real(r8), parameter :: c_relax_C         = 0.0030           ! Constant in the relaxation equation for the shape factor
                                                               ! with respect to the temperature profile in the thermocline
   real(r8), parameter :: h_Snow_min_flk    = 1.0E-5           ! Minimum snow thickness [m]
   real(r8), parameter :: h_Ice_min_flk     = 1.0E-9           ! Minimum ice thickness [m]
   real(r8), parameter :: h_ML_min_flk      = 1.0E-2           ! Minimum mixed-layer depth [m]
   real(r8), parameter :: h_ML_max_flk      = 1.0E+3           ! Maximum mixed-layer depth [m]
   real(r8), parameter :: H_B1_min_flk      = 1.0E-3           ! Minimum thickness of the upper layer of bottom sediments [m]
   real(r8), parameter :: u_star_min_flk    = 1.0E-6           ! Minimum value of the surface friction velocity [m s^{-1}   ]
   real(r8), parameter :: tpsf_C_StefBoltz  = 5.67E-08         ! The Stefan-Boltzmann constant [W m^{-2} K^{-4}]
   real(r8), parameter :: tpsf_R_dryair     = 2.8705E+02       ! Gas constant for dry air [J kg^{-1} K^{-1}]
   real(r8), parameter :: tpsf_R_watvap     = 4.6151E+02       ! Gas constant for water vapour [J kg^{-1} K^{-1}]
   real(r8), parameter :: tpsf_c_a_p        = 1.005E+03        ! Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
   real(r8), parameter :: tpsf_L_evap       = 2.501E+06        ! Specific heat of evaporation [J kg^{-1}]
   real(r8), parameter :: tpsf_nu_u_a       = 1.50E-05         ! Kinematic molecular viscosity of air [m^{2} s^{-1}]
   real(r8), parameter :: tpsf_kappa_t_a    = 2.20E-05         ! Molecular temperature conductivity of air [m^{2} s^{-1}]
   real(r8), parameter :: tpsf_kappa_q_a    = 2.40E-05         ! Molecular diffusivity of air for water vapour [m^{2} s^{-1}]
   real(r8), parameter :: c_lwrad_emis      = 0.99             ! Surface emissivity with respect to the long-wave radiation
   real(r8), parameter :: c_free_conv       = 0.14             ! Constant in the expressions for fluxes in free convection
   real(r8), parameter :: c_small_flk       = 1.0E-10          ! A small number
   real(r8), parameter :: tpl_T_r           = 277.13           ! Temperature of maximum density of fresh water [K]
   real(r8), parameter :: tpl_T_f           = 273.15           ! Fresh water freezing point [K]
   real(r8), parameter :: tpl_a_T           = 1.6509E-05       ! Constant in the fresh-water equation of state [K^{-2}]
   real(r8), parameter :: tpl_rho_w_r       = 1.0E+03          ! Maximum density of fresh water [kg m^{-3}]
   real(r8), parameter :: tpl_rho_I         = 9.1E+02          ! Density of ice [kg m^{-3}]
   real(r8), parameter :: tpl_rho_S_min     = 1.0E+02          ! Minimum snow density [kg m^{-3}]
   real(r8), parameter :: tpl_rho_S_max     = 4.0E+02          ! Maximum snow density [kg m^{-3}]
   real(r8), parameter :: tpl_Gamma_rho_S   = 2.0E+02          ! Empirical parameter [kg m^{-4}] in the expression for the snow density 
   real(r8), parameter :: tpl_L_f           = 3.3E+05          ! Latent heat of fusion [J kg^{-1}]
   real(r8), parameter :: tpl_c_w           = 4.2E+03          ! Specific heat of water [J kg^{-1} K^{-1}]
   real(r8), parameter :: tpl_c_I           = 2.1E+03          ! Specific heat of ice [J kg^{-1} K^{-1}]
   real(r8), parameter :: tpl_c_S           = 2.1E+03          ! Specific heat of snow [J kg^{-1} K^{-1}]
   real(r8), parameter :: tpl_kappa_w       = 5.46E-01         ! Molecular heat conductivity of water [J m^{-1} s^{-1} K^{-1}]
   real(r8), parameter :: tpl_kappa_I       = 2.29             ! Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]
   real(r8), parameter :: tpl_kappa_S_min   = 0.2              ! Minimum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
   real(r8), parameter :: tpl_kappa_S_max   = 1.5              ! Maximum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
   real(r8), parameter :: tpl_Gamma_kappa_S = 1.3              ! Empirical parameter [J m^{-2} s^{-1} K^{-1}] in the expression for the snow heat conductivity 
   real(r8), parameter :: u_wind_min_sf     = 1.0E-02          ! Minimum wind speed [m s^{-1}]
   real(r8), parameter :: u_star_min_sf     = 1.0E-04          ! Minimum value of friction velocity [m s^{-1}]
   real(r8), parameter :: c_accur_sf        = 1.0E-07          ! A small number (accuracy)
   real(r8), parameter :: c_small_sf        = 1.0E-04          ! A small number (used to compute fluxes)
   real(r8), parameter :: C_T_min           = 0.5              ! Minimum value of the shape factor C_T (thermocline)
   real(r8), parameter :: C_T_max           = 0.8              ! Maximum value of the shape factor C_T (thermocline)
   real(r8), parameter :: Phi_T_pr0_1       = 40./3.           ! Constant in the expression for the T shape-function derivative 
   real(r8), parameter :: Phi_T_pr0_2       = 20./3.           ! Constant in the expression for the T shape-function derivative 
   real(r8), parameter :: C_TT_1            = 11./18.          ! Constant in the expression for C_TT (thermocline)
   real(r8), parameter :: C_TT_2            = 7./45.           ! Constant in the expression for C_TT (thermocline)
   real(r8), parameter :: C_B1              = 2./3.            ! Shape factor (upper layer of bottom sediments)
   real(r8), parameter :: C_B2              = 3./5.            ! Shape factor (lower layer of bottom sediments)
   real(r8), parameter :: Phi_B1_pr0        = 2.               ! B1 shape-function derivative 
   real(r8), parameter :: C_S_lin           = 0.5              ! Shape factor (linear temperature profile in the snow layer)
   real(r8), parameter :: Phi_S_pr0_lin     = 1.               ! S shape-function derivative (linear profile) 
   real(r8), parameter :: C_I_lin           = 0.5              ! Shape factor (linear temperature profile in the ice layer)
   real(r8), parameter :: Phi_I_pr0_lin     = 1.               ! I shape-function derivative (linear profile) 
   real(r8), parameter :: Phi_I_pr1_lin     = 1.               ! I shape-function derivative (linear profile) 
   real(r8), parameter :: Phi_I_ast_MR      = 2.               ! Constant in the MR2004 expression for I shape factor
   real(r8), parameter :: C_I_MR            = 1./12.           ! Constant in the MR2004 expression for I shape factor
   real(r8), parameter :: H_Ice_max         = 3.               ! Maximum ice tickness in the Mironov and Ritter (2004, MR2004) ice model [m] 
   real(r8), parameter :: Pr_neutral        = 1.0              ! Turbulent Prandtl number at neutral static stability
   real(r8), parameter :: Sc_neutral        = 1.0              ! Turbulent Schmidt number at neutral static stability
   real(r8), parameter :: c_MO_u_stab       = 5.0              ! Constant of the MO theory (wind, stable stratification)
   real(r8), parameter :: c_MO_t_stab       = 5.0              ! Constant of the MO theory (temperature, stable stratification)
   real(r8), parameter :: c_MO_q_stab       = 5.0              ! Constant of the MO theory (humidity, stable stratification)
   real(r8), parameter :: c_MO_u_conv       = 15.0             ! Constant of the MO theory (wind, convection)
   real(r8), parameter :: c_MO_t_conv       = 15.0             ! Constant of the MO theory (temperature, convection)
   real(r8), parameter :: c_MO_q_conv       = 15.0             ! Constant of the MO theory (humidity, convection)
   real(r8), parameter :: c_MO_u_exp        = 0.25             ! Constant of the MO theory (wind, exponent)
   real(r8), parameter :: c_MO_t_exp        = 0.5              ! Constant of the MO theory (temperature, exponent)
   real(r8), parameter :: c_MO_q_exp        = 0.5              ! Constant of the MO theory (humidity, exponent)
   real(r8), parameter :: z0u_ice_rough     = 1.0E-03          ! Aerodynamic roughness of the ice surface [m] (rough flow)
   real(r8), parameter :: c_z0u_smooth      = 0.1              ! Constant in the expression for z0u (smooth flow) 
   real(r8), parameter :: c_z0u_rough       = 1.23E-02         ! The Charnock constant in the expression for z0u (rough flow)
   real(r8), parameter :: c_z0u_rough_L     = 1.00E-01         ! An increased Charnock constant (used as the upper limit)
   real(r8), parameter :: c_z0u_ftch_f      = 0.70             ! Factor in the expression for fetch-dependent Charnock parameter
   real(r8), parameter :: c_z0u_ftch_ex     = 1./3.            ! Exponent in the expression for fetch-dependent Charnock parameter
   real(r8), parameter :: c_z0t_rough_1     = 4.0              ! Constant in the expression for z0t (factor) 
   real(r8), parameter :: c_z0t_rough_2     = 3.2              ! Constant in the expression for z0t (factor)
   real(r8), parameter :: c_z0t_rough_3     = 0.5              ! Constant in the expression for z0t (exponent) 
   real(r8), parameter :: c_z0q_rough_1     = 4.0              ! Constant in the expression for z0q (factor)
   real(r8), parameter :: c_z0q_rough_2     = 4.2              ! Constant in the expression for z0q (factor)
   real(r8), parameter :: c_z0q_rough_3     = 0.5              ! Constant in the expression for z0q (exponent)
   real(r8), parameter :: c_z0t_ice_b0s     = 1.250            ! Constant in the expression for z0t over ice
   real(r8), parameter :: c_z0t_ice_b0t     = 0.149            ! Constant in the expression for z0t over ice
   real(r8), parameter :: c_z0t_ice_b1t     = -0.550           ! Constant in the expression for z0t over ice
   real(r8), parameter :: c_z0t_ice_b0r     = 0.317            ! Constant in the expression for z0t over ice
   real(r8), parameter :: c_z0t_ice_b1r     = -0.565           ! Constant in the expression for z0t over ice
   real(r8), parameter :: c_z0t_ice_b2r     = -0.183           ! Constant in the expression for z0t over ice
   real(r8), parameter :: c_z0q_ice_b0s     = 1.610            ! Constant in the expression for z0q over ice
   real(r8), parameter :: c_z0q_ice_b0t     = 0.351            ! Constant in the expression for z0q over ice
   real(r8), parameter :: c_z0q_ice_b1t     = -0.628           ! Constant in the expression for z0q over ice
   real(r8), parameter :: c_z0q_ice_b0r     = 0.396            ! Constant in the expression for z0q over ice
   real(r8), parameter :: c_z0q_ice_b1r     = -0.512           ! Constant in the expression for z0q over ice
   real(r8), parameter :: c_z0q_ice_b2r     = -0.180           ! Constant in the expression for z0q over ice
   real(r8), parameter :: c_albice_MR       = 95.6             ! Constant in the interpolation formula for the ice albedo (Mironov and Ritter 2004)
   real(r8), parameter :: Re_z0s_ice_t      = 2.5              ! Threshold value of the surface Reynolds number 
                                                               !  used to compute z0t and z0q over ice (Andreas 2002)
   real(r8), parameter :: Re_z0u_thresh     = 0.1              ! Threshold value of the roughness Reynolds number 
                                                               ! [value from Zilitinkevich, Grachev, and Fairall (200), currently not used] 

   ! Derived thermodynamic parameters
   real(r8), parameter :: tpsf_Rd_o_Rv = tpsf_R_dryair/tpsf_R_watvap  ! Ratio of gas constants (Rd/Rv)
   real(r8), parameter :: tpsf_alpha_q = (1.-tpsf_Rd_o_Rv)/tpsf_Rd_o_Rv! Diemsnionless ratio 

   logical, parameter  :: lflk_botsed_use   = .TRUE.           ! .TRUE. indicates that the bottom-sediment scheme is used
                                                               ! to compute the depth penetrated by the thermal wave, 
                                                               ! the temperature at this depth and the bottom heat flux.
                                                               ! Otherwise, the heat flux at the water-bottom sediment interface
                                                               ! is set to zero, the depth penetrated by the thermal wave 
                                                               ! is set to a reference value defined below,
                                                               ! and the temperature at this depth is set to 
                                                               ! the temperature of maximum density of the fresh water.

   real(r8), parameter :: rflk_depth_bs_ref = 10.0             ! Reference value of the depth of the thermally active
                                                               ! layer of bottom sediments [m].
                                                               ! This value is used to (formally) define
                                                               ! the depth penetrated by the thermal wave
                                                               ! in case the bottom-sediment scheme is not used.

   ! opticpar parameters
   integer , parameter :: nband_optic_max   = 10               ! Maximum value of the wave-length bands 
                                                               ! in the exponential decay law for the radiation flux.
                                                               ! A storage for a ten-band approximation is allocated,
                                                               ! although a smaller number of bands is actually used.
   type opticpar_medium
      integer  :: nband_optic                                  ! Number of wave-length bands
      real(r8) :: frac_optic(nband_optic_max)                  ! Fractions of total radiation flux
      real(r8) :: extincoef_optic(nband_optic_max)             ! Extinction coefficients
    END type opticpar_medium

    ! Optical characteristics for water, ice and snow.
    ! The simplest one-band approximation is used as a reference.
    integer, PRIVATE :: i
    TYPE (opticpar_medium), PARAMETER ::                   & 
        opticpar_water_trans  = opticpar_medium(2,         &   ! Transparent Water (two-band)
            (/0.10, 0.90, (0.,i=3,nband_optic_max)/),      &
            (/2.0,  0.20, (1.E+10,i=3,nband_optic_max)/)) ,&
        !_nu  opticpar_water_tr = opticpar_medium(1,         & ! Transparent Water (one-band)
        !_nu    (/1., (0.,i=2,nband_optic_max)/),            &
        !_nu    (/0.30, (1.E+10,i=2,nband_optic_max)/))    , &
        
        opticpar_water_ref    = opticpar_medium(1,         & ! Water (reference)
            (/1.,      (0.    , i=2, nband_optic_max)/),   &
            (/3.,      (1.E+10, i=2, nband_optic_max)/))  ,&
        
        opticpar_whiteice_ref = opticpar_medium(1,         & ! White ice
            (/1.,      (0.    , i=2, nband_optic_max)/),   &   
            (/17.1,    (1.E+10, i=2, nband_optic_max)/))  ,&
        
        opticpar_blueice_ref  = opticpar_medium(1,         & ! Blue ice
            (/1.,      (0.    , i=2, nband_optic_max)/),   &
            (/8.4,     (1.E+10, i=2, nband_optic_max)/))  ,&
        
        opticpar_drysnow_ref  = opticpar_medium(1,         & ! Dry snow 
            (/1.,      (0.    , i=2, nband_optic_max)/),   &
            (/25.0,    (1.E+10, i=2, nband_optic_max)/))  ,&
        
        opticpar_meltingsnow_ = opticpar_medium(1,         & ! Melting snow 
            (/1.,      (0.    , i=2, nband_optic_max)/),   &
            (/15.0,    (1.E+10, i=2, nband_optic_max)/))  ,&
        
        opticpar_ice_opaque   = opticpar_medium(1,         & ! Opaque ice
            (/1.,      (0.    , i=2, nband_optic_max)/),   &
            (/1.0E+07, (1.E+10, i=2, nband_optic_max)/))  ,&

        opticpar_snow_opaque  = opticpar_medium(1,         & ! Opaque snow
            (/1.,      (0.    , i=2, nband_optic_max)/),   &
            (/1.0E+07, (1.E+10, i=2, nband_optic_max)/)) 


!********************************************************
!               --- Simstrat Parameter ---               
!********************************************************
    ! Control parameters of the scheme, modify carefully
    logical, parameter ::  CoupleAED2           = .false.      ! Switch to turn on/off the biochemical model AED2
    logical, parameter ::  SplitSeicheParameter = .false.      ! True: USE a_seiche IF N2 exceeds strat_sumr and a_seiche_w otherwise, False: Always USE a_seiche
    logical, parameter ::  UseFilteredWind      = .false.      !-WMEJ Switch to activate filtered wind, current version is not available!!! 
    logical, parameter ::  has_advection        = .false.      !-WMEJ Switch to turn on/off advection, current version is not available!!!
    
    integer, parameter ::  TurbulenceModel      = 1            ! 1: k-epsilon, 2: Mellor-Yamada
    integer, parameter ::  StabilityFunction    = 2            ! 1: constant, 2: quasi-equilibrium
    integer, parameter ::  FluxCondition        = 1            ! 0: Dirichlet condition, 1: no-flux
    integer, parameter ::  ForcingMode          = 5            ! 1: Wind + Temp + SolRad,
                                                               ! 2: 1 + VapP,
                                                               ! 3: 2 + Cloud,
                                                               ! 4: Wind + HeatFlux + SolRad,
                                                               ! 5: 2 + Incoming Long Wave
    integer, parameter ::  SeicheNormalization  = 2            ! 1: max N2, 2: integral              
    integer, parameter ::  WindDragModel        = 3            ! 1: constant, 2: ocean (increasing), 3: lake (Wüest and Lor8e, 2003)
    integer, parameter ::  InflowPlacement      = 0            ! 0: inflow depths are chosen manually, 1: inflow is density-driven
    integer, parameter ::  PressureGradients    = 0            ! 0: off, 1: Svensson 1978, 2: Bottom drag
    integer, parameter ::  IceModel             = 1            ! Switch to turn on/off the ice model
    integer, parameter ::  SnowModel            = 0            ! Switch to turn on/off the snow model (needs ice model and an additional column in forcing)

    !-WMEJ These parameters are all adjustable parameters
    real(r8), parameter :: a_seiche             = 0.002        ! Fraction of wind energy which goes into seiche energy [-]
    real(r8), parameter :: a_seiche_w           = 0.           ! Fraction of wind energy which goes into seiche energy in winter [-] 
                                                               ! (only used IF SplitSeicheParameter=true)
    real(r8), parameter :: strat_sumr           = 0.           ! IF maximum N2 (Brundt-Vaisala frequency) is below this threshold,
                                                               ! a_seiche_w is used instead of a_seiche (only used IF SplitSeicheParameter=true)
    real(r8), parameter :: q_nn                 = 0.75         ! Fit parameter for distribution of seiche energy [-]
    real(r8), parameter :: f_wind               = 1.6423       ! Ratio of forcing wind to wind speed at 10 m above lake level [-]
    real(r8), parameter :: c10                  = 1.           ! Wind drag coefficient, a physical constant around 0.001 IF WindDragModel = 1
                                                               ! and a calibration parameter around 1 IF WindDragModel = 2 or 3
    real(r8), parameter :: cd                   = 0.002        ! Bottom drag coefficient [-]
    real(r8), parameter :: hgeo                 = 0.1          ! Geothermal heat flux [W/m2]
    real(r8), parameter :: p_sw                 = 1.           ! Fit parameter for absorption of short wave radiation from sky [-]
    real(r8), parameter :: p_lw                 = 1.1978       ! Fit parameter for absorption of IR radiation from sky [-]
    real(r8), parameter :: p_windf              = 1.           ! Fit parameter for convective and latent heat fluxes [-]
    real(r8), parameter :: beta_sol             = 0.3          ! Fraction of short-wave radiation directly absorbed as heat by water [-]
    real(r8), parameter :: p_albedo             = 1            ! Fit parameter for albedo of ice, snow-ice and snow [-] (only used IF IceModel = 1)
    real(r8), parameter :: freez_temp           = 0.00         ! Freezing temperature of water [°C] (only used IF IceModel = 1)
    real(r8), parameter :: snow_temp            = 1.           ! Temperature of snow [°C] (only used IF IceModel = 1 and SnowModel = 1)
    
    ! *** General constants ***
    real(r8), parameter :: rho_air          = 1.2              ! density of air [kg/m3]
    real(r8), parameter :: g                = 9.81             ! Earth gravitational acceleration [m/s2]
    real(r8), parameter :: kappa            = 0.41             ! Von Karman constant [-]
    real(r8), parameter :: K_s              = 0.05             ! Bottom roughness [m]
    real(r8), parameter :: z0               = 0.5              ! Surface roughness [m]

    ! *** Constants for freshwater ***
    real(r8), parameter :: rho_0            = 1000             ! Mean freshwater density (seawater: 1023) [kg/m3]
    real(r8), parameter :: cp               = 4182             ! Mean freshwater heat capacity (seawater: 3992) [J/kg/K]
    real(r8), parameter :: cp_air           = 1005             ! Mean air heat capacity [J/kg/K]

    ! *** Further parameters controlling water dynamic ***
    real(r8), parameter :: Prndtl           = 0.8              ! Prandtl number for air??
    real(r8), parameter :: k_min            = 1.0e-12          ! Min value allowed for turbulent kinetic energy
    real(r8), parameter :: eps_min          = 1.0e-30          ! Min value allowed for TKE dissipation
    real(r8), parameter :: avh_min          = 1.0e-8           ! Min value allowed for turbulent viscosity at boundaries

    ! *** Ice/Snow Parameters ***
    real(r8), parameter :: k_ice            = 2.22             ! Thermal conductivity ice at 0°C WaterTemp (W K-1 m-1)
    real(r8), parameter :: k_snow           = 0.2              ! Thermal conductivity snow at 0°C WaterTemp (W K-1 m-1)
    real(r8), parameter :: l_h              = 3.34e+5          ! Latent heat of melting, [J kg-1]
    real(r8), parameter :: l_e              = 0 !2.265e+6      ! Latent heat of evaporation, [J kg-1], l_e 0 => non-sublimation 
                                                               ! or l_e ~ 0 -> sublimation (i.e. solid to gas)  
    real(r8), parameter :: ice_dens         = 916.2            ! ice density kg m-3
    real(r8), parameter :: snowice_dens     = 875              ! snow ice density  [kg m-3] Saloranta 2000
    real(r8), parameter :: rho_s_0          = 250              ! New snowfall density [kg/m3]
    real(r8), parameter :: rho_s_max        = 450              ! maximum (wet) snow dens [kg/m3]
    real(r8), parameter :: cp_s             = 2090             ! Mean snow heat capacity (-5°C) [J/kg/K]
    real(r8), parameter :: emiss_water      = 0.97             ! Emissivity water
    real(r8), parameter :: emiss_ice        = 0.98             ! Emissivity ice
    ! Emissivity of snow ranges from 0.8 to 0.9 depending on snow density in strat_forcing
    real(r8), parameter :: C01              = 5.8/ 3600.       ! snow compresion in [m^-1 sec^-1] to match model timestep, initaly [m/h],
                                                               ! from Yen (1981) page 5, C1 range [2.6 to 9.0] 
    real(r8), parameter :: C02              = 21.0/ 1000.      ! snow compresion in [m^3/kg] initaly [m^3/Mg], from Yen (1981) page 5
    real(r8), parameter :: Ha_a             = 0.68             ! Longwave emision parameter, Matti Leppäranta 2009
    real(r8), parameter :: Ha_b             = 0.036            ! Longwave emision parameter, Matti Leppäranta 2009 [mbar^-1/2]
    real(r8), parameter :: Ha_c             = 0.18             ! Longwave emision parameter, Matti Leppäranta 2009
    real(r8), parameter :: Hk_CH            = 1.5e-3           ! convectiv bulk exchange coefficient, Matti Leppäranta 2009 and Gill 1982
    real(r8), parameter :: Hv_CE            = 1.5e-3           ! latent bulk exchange coefficient, Matti Leppäranta 2009 and Gill 1982
    real(r8), parameter :: lambda_snow      = 24               ! lambda (light absorption) snow [m-1]
    real(r8), parameter :: lambda_snowice   = 3                ! lambda (light absorption) white ice [m-1]
    real(r8), parameter :: lambda_ice       = 1                ! lambda (light absorption) black ice [m-1]

    ! *** Parameters for k-eps model ***
    real(r8), parameter :: ce1              = 1.44   
    real(r8), parameter :: ce2              = 1.92   
    real(r8), parameter :: ce3              = -0.4   
    real(r8), parameter :: sig_k            = 1.0   
    real(r8), save      :: sig_e            = 1.3   
    real(r8), parameter :: cmue             = 0.09   
    real(r8), parameter :: r_a              = 0.03             ! Ratio of reflected to total long-wave iradiance
    real(r8), parameter :: B0               = 0.61             ! Bowen constant

    ! *** Parameters for MY model ***
    real(r8), parameter :: a1               = 0.92   
    real(r8), parameter :: a2               = 0.74   
    real(r8), parameter :: b1               = 16.6   
    real(r8), parameter :: b2               = 10.1   
    real(r8), parameter :: c1               = 0.08   
    real(r8), parameter :: e1               = 1.8   
    real(r8), parameter :: e2               = 1.33   
    real(r8), parameter :: sl               = 0.2   


!********************************************************
!         --- Reference Parameter(NOT USED) ---               
!********************************************************
    ! Albedo for water, ice and snow. not used
    real(r8), parameter :: albedo_water_ref       = 0.07       ! Water
    real(r8), parameter :: albedo_whiteice_ref    = 0.60       ! White ice
    real(r8), parameter :: albedo_blueice_ref     = 0.10       ! Blue ice
    real(r8), parameter :: albedo_drysnow_ref     = 0.60       ! Dry snow 
    real(r8), parameter :: albedo_meltingsnow_ref = 0.10       ! Melting snow 
    real(r8), parameter :: ice_albedo             = 0.25       ! Albedo black ice
    real(r8), parameter :: snowice_albedo         = 0.35       ! Albedo white ice
    real(r8), parameter :: snow_albedo            = 0.7        ! Albedo snow


!*********************************************************
!         --- Reference Parameters For Inflows/Outflow ---               
!*********************************************************
    ! Inflow/outflow parameters
    real(r8), parameter :: strm_hf_angle = 65.0                ! Inflow angle [degree]
    real(r8), parameter :: strm_slope = 2.0                    ! Inflow slope [m/m]
    real(r8), parameter :: strm_drag_coeff = 0.0160            ! Inflow drag coefficient [-]
    real(r8), parameter :: InflowFactor = 1.0                  ! Inflow factor [-]
    real(r8), parameter :: OutflowFactor = 1.0                 ! Outflow factor [-]
    real(r8), parameter :: CrestFactor = 1.0                   ! Crest factor [-]
    real(r8), parameter :: cwnsq2 = 0.5
    real(r8), parameter :: zero = 0.0_r8


END MODULE  MOD_Lake_Const
