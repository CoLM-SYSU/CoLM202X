#include <define.h> 

MODULE MOD_2D_Forcing
! -------------------------------
! Meteorogical Forcing
!
! Created by Yongjiu Dai, 03/2014
! -------------------------------

use precision
IMPLICIT NONE
SAVE

! -----------------------------------------------------------------
  real(r8), allocatable :: forc_xy_pco2m (:,:) ! CO2 concentration in atmos. (pascals)
  real(r8), allocatable :: forc_xy_po2m  (:,:) ! O2 concentration in atmos. (pascals)
  real(r8), allocatable :: forc_xy_us    (:,:) ! wind in eastward direction [m/s]
  real(r8), allocatable :: forc_xy_vs    (:,:) ! wind in northward direction [m/s]
  real(r8), allocatable :: forc_xy_t     (:,:) ! temperature at reference height [kelvin]
  real(r8), allocatable :: forc_xy_q     (:,:) ! specific humidity at reference height [kg/kg]
  real(r8), allocatable :: forc_xy_prc   (:,:) ! convective precipitation [mm/s]
  real(r8), allocatable :: forc_xy_prl   (:,:) ! large scale precipitation [mm/s]
  real(r8), allocatable :: forc_xy_psrf  (:,:) ! atmospheric pressure at the surface [pa]
  real(r8), allocatable :: forc_xy_pbot  (:,:) ! atm bottom level pressure (or reference height) (pa)
  real(r8), allocatable :: forc_xy_sols  (:,:) ! atm vis direct beam solar rad onto srf [W/m2]
  real(r8), allocatable :: forc_xy_soll  (:,:) ! atm nir direct beam solar rad onto srf [W/m2]
  real(r8), allocatable :: forc_xy_solsd (:,:) ! atm vis diffuse solar rad onto srf [W/m2]
  real(r8), allocatable :: forc_xy_solld (:,:) ! atm nir diffuse solar rad onto srf [W/m2]
  real(r8), allocatable :: forc_xy_frl   (:,:) ! atmospheric infrared (longwave) radiation [W/m2]
  real(r8), allocatable :: forc_xy_hgt_u (:,:) ! observational height of wind [m]
  real(r8), allocatable :: forc_xy_hgt_t (:,:) ! observational height of temperature [m]
  real(r8), allocatable :: forc_xy_hgt_q (:,:) ! observational height of humidity [m]
  real(r8), allocatable :: forc_xy_rhoair(:,:) ! air density [kg/m3]

! PUBLIC MEMBER FUNCTIONS:
      public :: allocate_2D_Forcing
      public :: deallocate_2D_Forcing

! PRIVATE MEMBER FUNCTIONS:


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_2D_Forcing (lon_points,lat_points)
! ------------------------------------------------
! Allocates memory for CLM 2d [lon_points,lat_points] variables
! ------------------------------------------------
  use precision
  IMPLICIT NONE
  integer, INTENT(in) :: lon_points
  integer, INTENT(in) :: lat_points

      allocate ( forc_xy_pco2m  (lon_points,lat_points) ) ! CO2 concentration in atmos. (pascals)
      allocate ( forc_xy_po2m   (lon_points,lat_points) ) ! O2 concentration in atmos. (pascals)
      allocate ( forc_xy_us     (lon_points,lat_points) ) ! wind in eastward direction [m/s]
      allocate ( forc_xy_vs     (lon_points,lat_points) ) ! wind in northward direction [m/s]
      allocate ( forc_xy_t      (lon_points,lat_points) ) ! temperature at reference height [kelvin]
      allocate ( forc_xy_q      (lon_points,lat_points) ) ! specific humidity at reference height [kg/kg]
      allocate ( forc_xy_prc    (lon_points,lat_points) ) ! convective precipitation [mm/s]
      allocate ( forc_xy_prl    (lon_points,lat_points) ) ! large scale precipitation [mm/s]
      allocate ( forc_xy_psrf   (lon_points,lat_points) ) ! atmospheric pressure at the surface [pa]
      allocate ( forc_xy_pbot   (lon_points,lat_points) ) ! atm bottom level pressure (or reference height) (pa)
      allocate ( forc_xy_sols   (lon_points,lat_points) ) ! atm vis direct beam solar rad onto srf [W/m2]
      allocate ( forc_xy_soll   (lon_points,lat_points) ) ! atm nir direct beam solar rad onto srf [W/m2]
      allocate ( forc_xy_solsd  (lon_points,lat_points) ) ! atm vis diffuse solar rad onto srf [W/m2]
      allocate ( forc_xy_solld  (lon_points,lat_points) ) ! atm nir diffuse solar rad onto srf [W/m2]
      allocate ( forc_xy_frl    (lon_points,lat_points) ) ! atmospheric infrared (longwave) radiation [W/m2]
      allocate ( forc_xy_hgt_u  (lon_points,lat_points) ) ! observational height of wind [m]
      allocate ( forc_xy_hgt_t  (lon_points,lat_points) ) ! observational height of temperature [m]
      allocate ( forc_xy_hgt_q  (lon_points,lat_points) ) ! observational height of humidity [m]
      allocate ( forc_xy_rhoair (lon_points,lat_points) ) ! air density [kg/m3]

  END SUBROUTINE allocate_2D_Forcing


  SUBROUTINE deallocate_2D_Forcing
      deallocate ( forc_xy_pco2m  ) ! CO2 concentration in atmos. (pascals)
      deallocate ( forc_xy_po2m   ) ! O2 concentration in atmos. (pascals)
      deallocate ( forc_xy_us     ) ! wind in eastward direction [m/s]
      deallocate ( forc_xy_vs     ) ! wind in northward direction [m/s]
      deallocate ( forc_xy_t      ) ! temperature at reference height [kelvin]
      deallocate ( forc_xy_q      ) ! specific humidity at reference height [kg/kg]
      deallocate ( forc_xy_prc    ) ! convective precipitation [mm/s]
      deallocate ( forc_xy_prl    ) ! large scale precipitation [mm/s]
      deallocate ( forc_xy_psrf   ) ! atmospheric pressure at the surface [pa]
      deallocate ( forc_xy_pbot   ) ! atm bottom level pressure (or reference height) (pa)
      deallocate ( forc_xy_sols   ) ! atm vis direct beam solar rad onto srf [W/m2]
      deallocate ( forc_xy_soll   ) ! atm nir direct beam solar rad onto srf [W/m2]
      deallocate ( forc_xy_solsd  ) ! atm vis diffuse solar rad onto srf [W/m2]
      deallocate ( forc_xy_solld  ) ! atm nir diffuse solar rad onto srf [W/m2]
      deallocate ( forc_xy_frl    ) ! atmospheric infrared (longwave) radiation [W/m2]
      deallocate ( forc_xy_hgt_u  ) ! observational height of wind [m]
      deallocate ( forc_xy_hgt_t  ) ! observational height of temperature [m]
      deallocate ( forc_xy_hgt_q  ) ! observational height of humidity [m]
      deallocate ( forc_xy_rhoair ) ! air density [kg/m3]
  END SUBROUTINE deallocate_2D_Forcing

END MODULE MOD_2D_Forcing
! ------ EOP --------
