#include <define.h>

#ifdef OzoneStress
Module OzoneMod

  use precision
  use MOD_1D_Forcing, only: forc_ozone, forc_psrf
  USE PhysicalConstants, only: rgas
  use PFT_const, only: isevg, leaf_long, woody
  IMPLICIT NONE
  SAVE

  public :: CalcOzoneStress

  CONTAINS
  
  subroutine CalcOzoneStress (o3coefv,o3coefg,forc_ozone, forc_psrf, th, ram, &
                              rs, rb, lai, lai_old, ivt, o3uptake, deltim)
     ! convert o3 from mol/mol to nmol m^-3
     real(r8), intent(out)   :: o3coefv
     real(r8), intent(out)   :: o3coefg
     real(r8), intent(inout) :: forc_ozone 
     real(r8), intent(in)    :: forc_psrf
     real(r8), intent(in)    :: th
     real(r8), intent(in)    :: ram
     real(r8), intent(in)    :: rs
     real(r8), intent(in)    :: rb
     real(r8), intent(in)    :: lai
     real(r8), intent(in)    :: lai_old
     integer , intent(in)    :: ivt
     real(r8), intent(inout) :: o3uptake
     real(r8), intent(in)    :: deltim

     real(r8) :: o3concnmolm3   ! o3 concentration (nmol/m^3)
     real(r8) :: o3flux         ! instantaneous o3 flux (nmol m^-2 s^-1)
     real(r8) :: o3fluxcrit     ! instantaneous o3 flux beyond threshold (nmol m^-2 s^-1)
     real(r8) :: o3fluxperdt    ! o3 flux per timestep (mmol m^-2)
     real(r8) :: heal           ! o3uptake healing rate based on % of new leaves growing (mmol m^-2)
     real(r8) :: leafturn       ! leaf turnover time / mortality rate (per hour)
     real(r8) :: decay          ! o3uptake decay rate based on leaf lifetime (mmol m^-2)
     real(r8) :: photoInt       ! intercept for photosynthesis
     real(r8) :: photoSlope     ! slope for photosynthesis
     real(r8) :: condInt        ! intercept for conductance
     real(r8) :: condSlope      ! slope for conductance

     real(r8), parameter :: ko3 = 1.67_r8
 
    ! LAI threshold for LAIs that asymptote and don't reach 0
     real(r8), parameter :: lai_thresh = 0.5_r8
  
     ! threshold below which o3flux is set to 0 (nmol m^-2 s^-1)
     real(r8), parameter :: o3_flux_threshold = 0.8_r8
  
     ! o3 intercepts and slopes for photosynthesis
     real(r8), parameter :: needleleafPhotoInt   = 0.8390_r8  ! units = unitless
     real(r8), parameter :: needleleafPhotoSlope = 0._r8      ! units = per mmol m^-2
     real(r8), parameter :: broadleafPhotoInt    = 0.8752_r8  ! units = unitless
     real(r8), parameter :: broadleafPhotoSlope  = 0._r8      ! units = per mmol m^-2
     real(r8), parameter :: nonwoodyPhotoInt     = 0.8021_r8  ! units = unitless
     real(r8), parameter :: nonwoodyPhotoSlope   = -0.0009_r8 ! units = per mmol m^-2
  
     ! o3 intercepts and slopes for conductance
     real(r8), parameter :: needleleafCondInt    = 0.7823_r8  ! units = unitless
     real(r8), parameter :: needleleafCondSlope  = 0.0048_r8  ! units = per mmol m^-2
     real(r8), parameter :: broadleafCondInt     = 0.9125_r8  ! units = unitless
     real(r8), parameter :: broadleafCondSlope   = 0._r8      ! units = per mmol m^-2
     real(r8), parameter :: nonwoodyCondInt      = 0.7511_r8  ! units = unitless
     real(r8), parameter :: nonwoodyCondSlope    = 0._r8      ! units = per mmol m^-2

     forc_ozone = 100._r8 * 1.e-9_r8 ! ozone partial pressure [mol/mol]

     o3concnmolm3 = forc_ozone * 1.e9_r8 * (forc_psrf/(th*rgas*0.001_r8 ))

     ! calculate instantaneous flux
     o3flux = o3concnmolm3/ (ko3*rs+ rb + ram)

     ! apply o3 flux threshold
     if (o3flux < o3_flux_threshold) then
        o3fluxcrit = 0._r8
     else
        o3fluxcrit = o3flux - o3_flux_threshold
     endif

     ! calculate o3 flux per timestep
     o3fluxperdt = o3fluxcrit * deltim * 0.000001_r8

     if (lai > lai_thresh) then
        ! checking if new leaf area was added
        if (lai - lai_old > 0) then
           ! minimizing o3 damage to new leaves
           heal = max(0._r8,(((lai-lai_old)/lai)*o3fluxperdt))
        else
           heal = 0._r8
        endif

        if (isevg(ivt)) then
           leafturn = 1._r8/(leaf_long(ivt)*365._r8*24._r8)
        else
           leafturn = 0._r8
        endif

        ! o3 uptake decay based on leaf lifetime for evergreen plants
        decay = o3uptake * leafturn * deltim/3600._r8
        !cumulative uptake (mmol m^-2)
        o3uptake = max(0._r8, o3uptake + o3fluxperdt - decay - heal)

     else
        o3uptake = 0._r8
     end if

    if (o3uptake == 0._r8) then
        ! No o3 damage if no o3 uptake
        o3coefv = 1._r8
        o3coefg = 1._r8
     else
        ! Determine parameter values for this pft
        ! TODO(wjs, 2014-10-01) Once these parameters are moved into the params file,     this
        ! logic can be removed.
        if (ivt>3) then
           if (woody(ivt)==0) then
              photoInt   = nonwoodyPhotoInt
              photoSlope = nonwoodyPhotoSlope
              condInt    = nonwoodyCondInt
              condSlope  = nonwoodyCondSlope
           else
              photoInt   = broadleafPhotoInt
              photoSlope = broadleafPhotoSlope
              condInt    = broadleafCondInt
              condSlope  = broadleafCondSlope
           end if
        else
           photoInt   = needleleafPhotoInt
           photoSlope = needleleafPhotoSlope
           condInt    = needleleafCondInt
           condSlope  = needleleafCondSlope
        end if

        ! Apply parameter values to compute o3 coefficients
        o3coefv = max(0._r8, min(1._r8, photoInt + photoSlope * o3uptake))
        o3coefg = max(0._r8, min(1._r8, condInt  + condSlope  * o3uptake))

     end if

  end subroutine CalcOzoneStress
  end module OzoneMod
#endif
