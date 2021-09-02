MODULE PhysicalConstants

!=======================================================================
! physical constants 
!=======================================================================

  USE precision
  IMPLICIT NONE

  PUBLIC
  REAL(r8), parameter :: denice = 917.      ! density of ice [kg/m3]
  REAL(r8), parameter :: denh2o = 1000.     ! density of liquid water [kg/m3]
  REAL(r8), parameter :: cpliq  = 4188.     ! Specific heat of water [J/kg-K]
  REAL(r8), parameter :: cpice  = 2117.27   ! Specific heat of ice [J/kg-K]
  REAL(r8), parameter :: cpair  = 1004.64   ! specific heat of dry air [J/kg/K]
  REAL(r8), parameter :: hfus   = 0.3336e6  ! latent heat of fusion for ice [J/kg]
  REAL(r8), parameter :: hvap   = 2.5104e6  ! latent heat of evap for water [J/kg]
  REAL(r8), parameter :: hsub   = 2.8440e6  ! latent heat of sublimation [J/kg]
  REAL(r8), parameter :: tkair  = 0.023     ! thermal conductivity of air [W/m/k]
  REAL(r8), parameter :: tkice  = 2.290     ! thermal conductivity of ice [W/m/k]
  REAL(r8), parameter :: tkwat  = 0.6       ! thermal conductivity of water [W/m/k]
  REAL(r8), parameter :: tfrz   = 273.16    ! freezing temperature [K]
  REAL(r8), parameter :: rgas   = 287.04    ! gas constant for dry air [J/kg/K]
  REAL(r8), parameter :: roverg = 4.71047e4 ! rw/g = (8.3144/0.018)/(9.80616)*1000. mm/K
  REAL(r8), parameter :: rwat   = 461.296   ! gas constant for water vapor [J/(kg K)]
  REAL(r8), parameter :: grav   = 9.80616   ! gravity constant [m/s2]
  REAL(r8), parameter :: vonkar = 0.4       ! von Karman constant [-]
  REAL(r8), parameter :: stefnc = 5.67e-8   ! Stefan-Boltzmann constant  [W/m2/K4]

END MODULE PhysicalConstants
