#include <define.h>

MODULE MOD_2D_Fluxes
! ----------------------------------------------------------------------
! perfrom the grid average mapping: average a subgrid input 1d vector 
! of length numpatch to a output 2d array of length [lon_points,lat_points]
!
! Created by Yongjiu Dai, 03/2014
!---------------------------------------------------------------------

use precision
USE GlobalVars

IMPLICIT NONE
SAVE

integer,  allocatable :: mask     (:,:)  ! grid mask flag
real(r8), allocatable :: frac     (:,:)  ! grid total fraction
real(r8), allocatable :: area     (:,:)  ! grid cell area

real(r8), allocatable :: f_taux   (:,:)  ! wind stress: E-W [kg/m/s2]
real(r8), allocatable :: f_tauy   (:,:)  ! wind stress: N-S [kg/m/s2]
real(r8), allocatable :: f_fsena  (:,:)  ! sensible heat from canopy height to atmosphere [W/m2]
real(r8), allocatable :: f_lfevpa (:,:)  ! latent heat flux from canopy height to atmosphere [W/m2]
real(r8), allocatable :: f_fevpa  (:,:)  ! evapotranspiration from canopy to atmosphere [mm/s]
real(r8), allocatable :: f_fsenl  (:,:)  ! sensible heat from leaves [W/m2]
real(r8), allocatable :: f_fevpl  (:,:)  ! evaporation+transpiration from leaves [mm/s]
real(r8), allocatable :: f_etr    (:,:)  ! transpiration rate [mm/s]
real(r8), allocatable :: f_fseng  (:,:)  ! sensible heat flux from ground [W/m2]
real(r8), allocatable :: f_fevpg  (:,:)  ! evaporation heat flux from ground [mm/s]
real(r8), allocatable :: f_fgrnd  (:,:)  ! ground heat flux [W/m2]
real(r8), allocatable :: f_sabvsun(:,:)  ! solar absorbed by sunlit canopy [W/m2]
real(r8), allocatable :: f_sabvsha(:,:)  ! solar absorbed by shaded [W/m2]
real(r8), allocatable :: f_sabg   (:,:)  ! solar absorbed by ground  [W/m2]
real(r8), allocatable :: f_sr     (:,:)  ! total reflected solar radiation (W/m2)
real(r8), allocatable :: f_solvd  (:,:)  ! incident direct beam vis solar radiation (W/m2)
real(r8), allocatable :: f_solvi  (:,:)  ! incident diffuse beam vis solar radiation (W/m2)
real(r8), allocatable :: f_solnd  (:,:)  ! incident direct beam nir solar radiation (W/m2)
real(r8), allocatable :: f_solni  (:,:)  ! incident diffuse beam nir solar radiation (W/m2)
real(r8), allocatable :: f_srvd   (:,:)  ! reflected direct beam vis solar radiation (W/m2)
real(r8), allocatable :: f_srvi   (:,:)  ! reflected diffuse beam vis solar radiation (W/m2)
real(r8), allocatable :: f_srnd   (:,:)  ! reflected direct beam nir solar radiation (W/m2)
real(r8), allocatable :: f_srni   (:,:)  ! reflected diffuse beam nir solar radiation (W/m2)
real(r8), allocatable :: f_solvdln(:,:)  ! incident direct beam vis solar radiation at local noon (W/m2)
real(r8), allocatable :: f_solviln(:,:)  ! incident diffuse beam vis solar radiation at local noon (W/m2)
real(r8), allocatable :: f_solndln(:,:)  ! incident direct beam nir solar radiation at local noon (W/m2)
real(r8), allocatable :: f_solniln(:,:)  ! incident diffuse beam nir solar radiation at local noon (W/m2)
real(r8), allocatable :: f_srvdln (:,:)  ! reflected direct beam vis solar radiation at local noon (W/m2)
real(r8), allocatable :: f_srviln (:,:)  ! reflected diffuse beam vis solar radiation at local noon (W/m2)
real(r8), allocatable :: f_srndln (:,:)  ! reflected direct beam nir solar radiation at local noon (W/m2)
real(r8), allocatable :: f_srniln (:,:)  ! reflected diffuse beam nir solar radiation at local noon (W/m2)
real(r8), allocatable :: f_olrg   (:,:)  ! outgoing long-wave radiation from ground+canopy [W/m2]
real(r8), allocatable :: f_rnet   (:,:)  ! net radiation [W/m2]
real(r8), allocatable :: f_xerr   (:,:)  ! the error of water banace [mm/s]
real(r8), allocatable :: f_zerr   (:,:)  ! the error of energy balance [W/m2]
real(r8), allocatable :: f_rsur   (:,:)  ! surface runoff [mm/s]
real(r8), allocatable :: f_rnof   (:,:)  ! total runoff [mm/s]
real(r8), allocatable :: f_qintr  (:,:)  ! interception [mm/s]
real(r8), allocatable :: f_qinfl  (:,:)  ! inflitration [mm/s]
real(r8), allocatable :: f_qdrip  (:,:)  ! throughfall [mm/s]
real(r8), allocatable :: f_assim  (:,:)  ! canopy assimilation rate [mol m-2 s-1]
real(r8), allocatable :: f_respc  (:,:)  ! respiration (plant+soil) [mol m-2 s-1]
real(r8), allocatable :: f_qcharge(:,:)  ! groundwater recharge rate [mm/s] 

!---------------------------------------------------------------------
real(r8), allocatable :: f_t_grnd (:,:)  ! ground surface temperature [K]
real(r8), allocatable :: f_tleaf  (:,:)  ! sunlit leaf temperature [K]
real(r8), allocatable :: f_ldew   (:,:)  ! depth of water on foliage [mm]
real(r8), allocatable :: f_scv    (:,:)  ! snow cover, water equivalent [mm]
real(r8), allocatable :: f_snowdp (:,:)  ! snow depth [meter]
real(r8), allocatable :: f_fsno   (:,:)  ! fraction of snow cover on ground
real(r8), allocatable :: f_sigf   (:,:)  ! fraction of veg cover, excluding snow-covered veg [-]
real(r8), allocatable :: f_green  (:,:)  ! leaf greenness
real(r8), allocatable :: f_lai    (:,:)  ! leaf area index
real(r8), allocatable :: f_laisun (:,:)  ! sunlit leaf area index
real(r8), allocatable :: f_laisha (:,:)  ! shaded leaf area index
real(r8), allocatable :: f_sai    (:,:)  ! stem area index
real(r8), allocatable :: f_alb(:,:,:,:)  ! averaged albedo [visible, direct; direct, diffuse]
real(r8), allocatable :: f_emis   (:,:)  ! averaged bulk surface emissivity
real(r8), allocatable :: f_z0m    (:,:)  ! effective roughness [m]
real(r8), allocatable :: f_trad   (:,:)  ! radiative temperature of surface [K]
real(r8), allocatable :: f_tref   (:,:)  ! 2 m height air temperature [kelvin]
real(r8), allocatable :: f_qref   (:,:)  ! 2 m height air specific humidity [kg/kg]

!---------------------------------------------------------------------
real(r8), allocatable :: f_t_soisno   (:,:,:)  ! soil temperature [K]
real(r8), allocatable :: f_wliq_soisno(:,:,:)  ! liquid water in soil layers [kg/m2]
real(r8), allocatable :: f_wice_soisno(:,:,:)  ! ice lens in soil layers [kg/m2]
real(r8), allocatable :: f_h2osoi     (:,:,:)  ! volumetric soil water in layers [m3/m3]
real(r8), allocatable :: f_rstfac     (:,:)    ! factor of soil water stress 
real(r8), allocatable :: f_zwt        (:,:)    ! the depth to water table [m]
real(r8), allocatable :: f_wa         (:,:)    ! water storage in aquifer [mm]
real(r8), allocatable :: f_wat        (:,:)    ! total water storage [mm]

real(r8), allocatable :: f_t_lake     (:,:,:) ! lake temperature [K]
real(r8), allocatable :: f_lake_icefrac (:,:,:) ! lake ice fraction cover [0-1]

!---------------------------------------------------------------------
real(r8), allocatable :: f_ustar  (:,:)  ! u* in similarity theory [m/s]
real(r8), allocatable :: f_tstar  (:,:)  ! t* in similarity theory [kg/kg]
real(r8), allocatable :: f_qstar  (:,:)  ! q* in similarity theory [kg/kg]
real(r8), allocatable :: f_zol    (:,:)  ! dimensionless height (z/L) used in Monin-Obukhov theory
real(r8), allocatable :: f_rib    (:,:)  ! bulk Richardson number in surface layer
real(r8), allocatable :: f_fm     (:,:)  ! integral of profile function for momentum
real(r8), allocatable :: f_fh     (:,:)  ! integral of profile function for heat
real(r8), allocatable :: f_fq     (:,:)  ! integral of profile function for moisture
real(r8), allocatable :: f_us10m  (:,:)  ! 10m u-velocity [m/s]
real(r8), allocatable :: f_vs10m  (:,:)  ! 10m v-velocity [m/s]
real(r8), allocatable :: f_fm10m  (:,:)  ! integral of profile function for momentum at 10m [-]

!---------------------------------------------------------------------
real(r8), allocatable :: f_xy_us  (:,:)  ! wind in eastward direction [m/s]
real(r8), allocatable :: f_xy_vs  (:,:)  ! wind in northward direction [m/s]
real(r8), allocatable :: f_xy_t   (:,:)  ! temperature at reference height [kelvin]
real(r8), allocatable :: f_xy_q   (:,:)  ! specific humidity at reference height [kg/kg]
real(r8), allocatable :: f_xy_prc (:,:)  ! convective precipitation [mm/s]
real(r8), allocatable :: f_xy_prl (:,:)  ! large scale precipitation [mm/s]
real(r8), allocatable :: f_xy_pbot(:,:)  ! atmospheric pressure at the surface [pa]
real(r8), allocatable :: f_xy_frl (:,:)  ! atmospheric infrared (longwave) radiation [W/m2]
real(r8), allocatable :: f_xy_solarin(:,:)  ! downward solar radiation at surface [W/m2]
real(r8), allocatable :: f_xy_rain(:,:)  ! rain [mm/s]
real(r8), allocatable :: f_xy_snow(:,:)  ! snow [mm/s]


! PUBLIC MEMBER FUNCTIONS:
      public :: allocate_2D_Fluxes
      public :: deallocate_2D_Fluxes
      public :: FLUSH_2D_Fluxes

! PRIVATE MEMBER FUNCTIONS:


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

SUBROUTINE allocate_2D_Fluxes(lon_points,lat_points)
! --------------------------------------------------------------------
! Allocates memory for CLM 2d [lon_points,lat_points] variables
! --------------------------------------------------------------------

use precision
IMPLICIT NONE

Integer, INTENT(in) :: lon_points
Integer, INTENT(in) :: lat_points

allocate ( mask     (lon_points,lat_points) )  ! grid mask
allocate ( frac     (lon_points,lat_points) )  ! grid total fraction
allocate ( area     (lon_points,lat_points) )  ! grid cell area
allocate ( f_taux   (lon_points,lat_points) )  ! wind stress: E-W [kg/m/s2]
allocate ( f_tauy   (lon_points,lat_points) )  ! wind stress: N-S [kg/m/s2]
allocate ( f_fsena  (lon_points,lat_points) )  ! sensible heat from canopy height to atmosphere [W/m2]
allocate ( f_lfevpa (lon_points,lat_points) )  ! latent heat flux from canopy height to atmosphere [W/m2]
allocate ( f_fevpa  (lon_points,lat_points) )  ! evapotranspiration from canopy to atmosphere [mm/s]
allocate ( f_fsenl  (lon_points,lat_points) )  ! sensible heat from leaves [W/m2]
allocate ( f_fevpl  (lon_points,lat_points) )  ! evaporation+transpiration from leaves [mm/s]
allocate ( f_etr    (lon_points,lat_points) )  ! transpiration rate [mm/s]
allocate ( f_fseng  (lon_points,lat_points) )  ! sensible heat flux from ground [W/m2]
allocate ( f_fevpg  (lon_points,lat_points) )  ! evaporation heat flux from ground [mm/s]
allocate ( f_fgrnd  (lon_points,lat_points) )  ! ground heat flux [W/m2]
allocate ( f_sabvsun(lon_points,lat_points) )  ! solar absorbed by sunlit canopy [W/m2]
allocate ( f_sabvsha(lon_points,lat_points) )  ! solar absorbed by shaded [W/m2]
allocate ( f_sabg   (lon_points,lat_points) )  ! solar absorbed by ground  [W/m2]
allocate ( f_sr     (lon_points,lat_points) )  ! total reflected solar radiation (W/m2)
allocate ( f_solvd  (lon_points,lat_points) )  ! incident direct beam vis solar radiation (W/m2)
allocate ( f_solvi  (lon_points,lat_points) )  ! incident diffuse beam vis solar radiation (W/m2)
allocate ( f_solnd  (lon_points,lat_points) )  ! incident direct beam nir solar radiation (W/m2)
allocate ( f_solni  (lon_points,lat_points) )  ! incident diffuse beam nir solar radiation (W/m2)
allocate ( f_srvd   (lon_points,lat_points) )  ! reflected direct beam vis solar radiation (W/m2)
allocate ( f_srvi   (lon_points,lat_points) )  ! reflected diffuse beam vis solar radiation (W/m2)
allocate ( f_srnd   (lon_points,lat_points) )  ! reflected direct beam nir solar radiation (W/m2)
allocate ( f_srni   (lon_points,lat_points) )  ! reflected diffuse beam nir solar radiation (W/m2)
allocate ( f_solvdln(lon_points,lat_points) )  ! incident direct beam vis solar radiation at local noon(W/m2)
allocate ( f_solviln(lon_points,lat_points) )  ! incident diffuse beam vis solar radiation at local noon(W/m2)
allocate ( f_solndln(lon_points,lat_points) )  ! incident direct beam nir solar radiation at local noon(W/m2)
allocate ( f_solniln(lon_points,lat_points) )  ! incident diffuse beam nir solar radiation at local noon(W/m2)
allocate ( f_srvdln (lon_points,lat_points) )  ! reflected direct beam vis solar radiation at local noon(W/m2)
allocate ( f_srviln (lon_points,lat_points) )  ! reflected diffuse beam vis solar radiation at local noon(W/m2)
allocate ( f_srndln (lon_points,lat_points) )  ! reflected direct beam nir solar radiation at local noon(W/m2)
allocate ( f_srniln (lon_points,lat_points) )  ! reflected diffuse beam nir solar radiation at local noon(W/m2)
allocate ( f_olrg   (lon_points,lat_points) )  ! outgoing long-wave radiation from ground+canopy [W/m2]
allocate ( f_rnet   (lon_points,lat_points) )  ! net radiation [W/m2]
allocate ( f_xerr   (lon_points,lat_points) )  ! the error of water banace [mm/s]
allocate ( f_zerr   (lon_points,lat_points) )  ! the error of energy balance [W/m2]
allocate ( f_rsur   (lon_points,lat_points) )  ! surface runoff [mm/s]
allocate ( f_rnof   (lon_points,lat_points) )  ! total runoff [mm/s]
allocate ( f_qintr  (lon_points,lat_points) )  ! interception [mm/s]
allocate ( f_qinfl  (lon_points,lat_points) )  ! inflitration [mm/s]
allocate ( f_qdrip  (lon_points,lat_points) )  ! throughfall [mm/s]
allocate ( f_assim  (lon_points,lat_points) )  ! canopy assimilation rate [mol m-2 s-1]
allocate ( f_respc  (lon_points,lat_points) )  ! respiration (plant+soil) [mol m-2 s-1]
allocate ( f_qcharge(lon_points,lat_points) )  ! groundwater recharge rate [mm/s] 

!---------------------------------------------------------------------
allocate ( f_t_grnd (lon_points,lat_points) )  ! ground surface temperature [K]
allocate ( f_tleaf  (lon_points,lat_points) )  ! sunlit leaf temperature [K]
allocate ( f_ldew   (lon_points,lat_points) )  ! depth of water on foliage [mm]
allocate ( f_scv    (lon_points,lat_points) )  ! snow cover, water equivalent [mm]
allocate ( f_snowdp (lon_points,lat_points) )  ! snow depth [meter]
allocate ( f_fsno   (lon_points,lat_points) )  ! fraction of snow cover on ground
allocate ( f_sigf   (lon_points,lat_points) )  ! fraction of veg cover, excluding snow-covered veg [-]
allocate ( f_green  (lon_points,lat_points) )  ! leaf greenness
allocate ( f_lai    (lon_points,lat_points) )  ! leaf area index
allocate ( f_laisun (lon_points,lat_points) )  ! sunlit leaf area index
allocate ( f_laisha (lon_points,lat_points) )  ! shaded leaf area index
allocate ( f_sai    (lon_points,lat_points) )  ! stem area index
allocate ( f_alb(2,2,lon_points,lat_points) )  ! averaged albedo [visible, direct; direct, diffuse]
allocate ( f_emis   (lon_points,lat_points) )  ! averaged bulk surface emissivity
allocate ( f_z0m    (lon_points,lat_points) )  ! effective roughness [m]
allocate ( f_trad   (lon_points,lat_points) )  ! radiative temperature of surface [K]
allocate ( f_tref   (lon_points,lat_points) )  ! 2 m height air temperature [kelvin]
allocate ( f_qref   (lon_points,lat_points) )  ! 2 m height air specific humidity [kg/kg]

!---------------------------------------------------------------------
allocate ( f_t_soisno   (maxsnl+1:nl_soil,lon_points,lat_points) )  ! soil temperature [K]
allocate ( f_wliq_soisno(maxsnl+1:nl_soil,lon_points,lat_points) )  ! liquid water in soil layers [kg/m2]
allocate ( f_wice_soisno(maxsnl+1:nl_soil,lon_points,lat_points) )  ! ice lens in soil layers [kg/m2]
allocate ( f_h2osoi            (1:nl_soil,lon_points,lat_points) )  ! volumetric soil water in layers [m3/m3]
allocate ( f_rstfac                      (lon_points,lat_points) )  ! factor of soil water stress 
allocate ( f_zwt                         (lon_points,lat_points) )  ! the depth to water table [m]
allocate ( f_wa                          (lon_points,lat_points) )  ! water storage in aquifer [mm]
allocate ( f_wat                         (lon_points,lat_points) )  ! total water storage [mm]

allocate ( f_t_lake      (nl_lake,lon_points,lat_points) )  ! lake temperature [K]
allocate ( f_lake_icefrac(nl_lake,lon_points,lat_points) )  ! lake ice fraction cover [0-1]

!---------------------------------------------------------------------
allocate ( f_ustar  (lon_points,lat_points) )  ! u* in similarity theory [m/s]
allocate ( f_tstar  (lon_points,lat_points) )  ! t* in similarity theory [kg/kg]
allocate ( f_qstar  (lon_points,lat_points) )  ! q* in similarity theory [kg/kg]
allocate ( f_zol    (lon_points,lat_points) )  ! dimensionless height (z/L) used in Monin-Obukhov theory
allocate ( f_rib    (lon_points,lat_points) )  ! bulk Richardson number in surface layer
allocate ( f_fm     (lon_points,lat_points) )  ! integral of profile function for momentum
allocate ( f_fh     (lon_points,lat_points) )  ! integral of profile function for heat
allocate ( f_fq     (lon_points,lat_points) )  ! integral of profile function for moisture
allocate ( f_us10m  (lon_points,lat_points) )  ! 10m u-velocity [m/s]
allocate ( f_vs10m  (lon_points,lat_points) )  ! 10m v-velocity [m/s]
allocate ( f_fm10m  (lon_points,lat_points) )  ! integral of profile function for momentum at 10m [-]

!---------------------------------------------------------------------
allocate ( f_xy_us  (lon_points,lat_points) )  ! wind in eastward direction [m/s]
allocate ( f_xy_vs  (lon_points,lat_points) )  ! wind in northward direction [m/s]
allocate ( f_xy_t   (lon_points,lat_points) )  ! temperature at reference height [kelvin]
allocate ( f_xy_q   (lon_points,lat_points) )  ! specific humidity at reference height [kg/kg]
allocate ( f_xy_prc (lon_points,lat_points) )  ! convective precipitation [mm/s]
allocate ( f_xy_prl (lon_points,lat_points) )  ! large scale precipitation [mm/s]
allocate ( f_xy_pbot(lon_points,lat_points) )  ! atmospheric pressure at the surface [pa]
allocate ( f_xy_frl (lon_points,lat_points) )  ! atmospheric infrared (longwave) radiation [W/m2]
allocate ( f_xy_solarin(lon_points,lat_points) )  ! downward solar radiation at surface [W/m2]
allocate ( f_xy_rain(lon_points,lat_points) )  ! rain [mm/s]
allocate ( f_xy_snow(lon_points,lat_points) )  ! snow [mm/s]


END SUBROUTINE allocate_2D_Fluxes



SUBROUTINE FLUSH_2D_Fluxes

! flush the 2D_Fluxes for accumulation
f_xy_us     (:,:) = spval
f_xy_vs     (:,:) = spval
f_xy_t      (:,:) = spval
f_xy_q      (:,:) = spval
f_xy_prc    (:,:) = spval
f_xy_prl    (:,:) = spval
f_xy_pbot   (:,:) = spval
f_xy_frl    (:,:) = spval
f_xy_solarin(:,:) = spval
f_xy_rain   (:,:) = spval
f_xy_snow   (:,:) = spval

mask        (:,:) = 0
frac        (:,:) = spval
area        (:,:) = spval

f_taux      (:,:) = spval
f_tauy      (:,:) = spval
f_fsena     (:,:) = spval
f_lfevpa    (:,:) = spval
f_fevpa     (:,:) = spval
f_fsenl     (:,:) = spval
f_fevpl     (:,:) = spval
f_etr       (:,:) = spval
f_fseng     (:,:) = spval
f_fevpg     (:,:) = spval
f_fgrnd     (:,:) = spval
f_sabvsun   (:,:) = spval
f_sabvsha   (:,:) = spval
f_sabg      (:,:) = spval
f_sr        (:,:) = spval
f_solvd     (:,:) = spval
f_solvi     (:,:) = spval
f_solnd     (:,:) = spval
f_solni     (:,:) = spval
f_srvd      (:,:) = spval
f_srvi      (:,:) = spval
f_srnd      (:,:) = spval
f_srni      (:,:) = spval
f_solvdln   (:,:) = spval
f_solviln   (:,:) = spval
f_solndln   (:,:) = spval
f_solniln   (:,:) = spval
f_srvdln    (:,:) = spval
f_srviln    (:,:) = spval
f_srndln    (:,:) = spval
f_srniln    (:,:) = spval
f_olrg      (:,:) = spval
f_rnet      (:,:) = spval
f_xerr      (:,:) = spval
f_zerr      (:,:) = spval
f_rsur      (:,:) = spval
f_rnof      (:,:) = spval
f_qintr     (:,:) = spval
f_qinfl     (:,:) = spval
f_qdrip     (:,:) = spval
f_assim     (:,:) = spval
f_respc     (:,:) = spval

f_qcharge   (:,:) = spval

f_t_grnd    (:,:) = spval
f_tleaf     (:,:) = spval
f_ldew      (:,:) = spval
f_scv       (:,:) = spval
f_snowdp    (:,:) = spval
f_fsno      (:,:) = spval
f_sigf      (:,:) = spval
f_green     (:,:) = spval
f_lai       (:,:) = spval
f_laisun    (:,:) = spval
f_laisha    (:,:) = spval
f_sai       (:,:) = spval
f_alb       (:,:,:,:) = spval
f_emis      (:,:) = spval
f_z0m       (:,:) = spval
f_trad      (:,:) = spval
f_tref      (:,:) = spval
f_qref      (:,:) = spval

f_t_soisno    (:,:,:) = spval
f_wliq_soisno (:,:,:) = spval
f_wice_soisno (:,:,:) = spval
f_h2osoi      (:,:,:) = spval
f_rstfac        (:,:) = spval
f_zwt           (:,:) = spval
f_wa            (:,:) = spval
f_wat           (:,:) = spval
f_t_lake      (:,:,:) = spval
f_lake_icefrac(:,:,:) = spval

f_ustar     (:,:) = spval
f_tstar     (:,:) = spval
f_qstar     (:,:) = spval
f_zol       (:,:) = spval
f_rib       (:,:) = spval
f_fm        (:,:) = spval
f_fh        (:,:) = spval
f_fq        (:,:) = spval

f_us10m     (:,:) = spval
f_vs10m     (:,:) = spval
f_fm10m     (:,:) = spval


END SUBROUTINE FLUSH_2D_Fluxes



SUBROUTINE deallocate_2D_Fluxes
! --------------------------------------------------------------------
! Deallocates memory for CLM 2d [:,:] variables
! --------------------------------------------------------------------

deallocate ( mask     )  ! grid mask
deallocate ( frac     )  ! grid total fraction
deallocate ( area     )  ! grid cell area
deallocate ( f_taux   )  ! wind stress: E-W [kg/m/s2]
deallocate ( f_tauy   )  ! wind stress: N-S [kg/m/s2]
deallocate ( f_fsena  )  ! sensible heat from canopy height to atmosphere [W/m2]
deallocate ( f_lfevpa )  ! latent heat flux from canopy height to atmosphere [W/m2]
deallocate ( f_fevpa  )  ! evapotranspiration from canopy to atmosphere [mm/s]
deallocate ( f_fsenl  )  ! sensible heat from leaves [W/m2]
deallocate ( f_fevpl  )  ! evaporation+transpiration from leaves [mm/s]
deallocate ( f_etr    )  ! transpiration rate [mm/s]
deallocate ( f_fseng  )  ! sensible heat flux from ground [W/m2]
deallocate ( f_fevpg  )  ! evaporation heat flux from ground [mm/s]
deallocate ( f_fgrnd  )  ! ground heat flux [W/m2]
deallocate ( f_sabvsun)  ! solar absorbed by sunlit canopy [W/m2]
deallocate ( f_sabvsha)  ! solar absorbed by shaded [W/m2]
deallocate ( f_sabg   )  ! solar absorbed by ground  [W/m2]
deallocate ( f_sr     )  ! total reflected solar radiation (W/m2)
deallocate ( f_solvd  )  ! incident direct beam vis solar radiation (W/m2)
deallocate ( f_solvi  )  ! incident diffuse beam vis solar radiation (W/m2)
deallocate ( f_solnd  )  ! incident direct beam nir solar radiation (W/m2)
deallocate ( f_solni  )  ! incident diffuse beam nir solar radiation (W/m2)
deallocate ( f_srvd   )  ! reflected direct beam vis solar radiation (W/m2)
deallocate ( f_srvi   )  ! reflected diffuse beam vis solar radiation (W/m2)
deallocate ( f_srnd   )  ! reflected direct beam nir solar radiation (W/m2)
deallocate ( f_srni   )  ! reflected diffuse beam nir solar radiation (W/m2)
deallocate ( f_solvdln)  ! incident direct beam vis solar radiation at local noon(W/m2)
deallocate ( f_solviln)  ! incident diffuse beam vis solar radiation at local noon(W/m2)
deallocate ( f_solndln)  ! incident direct beam nir solar radiation at local noon(W/m2)
deallocate ( f_solniln)  ! incident diffuse beam nir solar radiation at local noon(W/m2)
deallocate ( f_srvdln )  ! reflected direct beam vis solar radiation at local noon(W/m2)
deallocate ( f_srviln )  ! reflected diffuse beam vis solar radiation at local noon(W/m2)
deallocate ( f_srndln )  ! reflected direct beam nir solar radiation at local noon(W/m2)
deallocate ( f_srniln )  ! reflected diffuse beam nir solar radiation at local noon(W/m2)
deallocate ( f_olrg   )  ! outgoing long-wave radiation from ground+canopy [W/m2]
deallocate ( f_rnet   )  ! net radiation [W/m2]
deallocate ( f_xerr   )  ! the error of water banace [mm/s]
deallocate ( f_zerr   )  ! the error of energy balance [W/m2]
deallocate ( f_rsur   )  ! surface runoff [mm/s]
deallocate ( f_rnof   )  ! total runoff [mm/s]
deallocate ( f_qintr  )  ! interception [mm/s]
deallocate ( f_qinfl  )  ! inflitration [mm/s]
deallocate ( f_qdrip  )  ! throughfall [mm/s]
deallocate ( f_assim  )  ! canopy assimilation rate [mol m-2 s-1]
deallocate ( f_respc  )  ! respiration (plant+soil) [mol m-2 s-1]
deallocate ( f_qcharge)  ! groundwater recharge rate [mm/s] 

!---------------------------------------------------------------------
deallocate ( f_t_grnd )  ! ground surface temperature [K]
deallocate ( f_tleaf  )  ! sunlit leaf temperature [K]
deallocate ( f_ldew   )  ! depth of water on foliage [mm]
deallocate ( f_scv    )  ! snow cover, water equivalent [mm]
deallocate ( f_snowdp )  ! snow depth [meter]
deallocate ( f_fsno   )  ! fraction of snow cover on ground
deallocate ( f_sigf   )  ! fraction of veg cover, excluding snow-covered veg [-]
deallocate ( f_green  )  ! leaf greenness
deallocate ( f_lai    )  ! leaf area index
deallocate ( f_laisun )  ! sunlit leaf area index
deallocate ( f_laisha )  ! shaded leaf area index
deallocate ( f_sai    )  ! stem area index
deallocate ( f_alb    )  ! averaged albedo [visible, direct; direct, diffuse]
deallocate ( f_emis   )  ! averaged bulk surface emissivity
deallocate ( f_z0m    )  ! effective roughness [m]
deallocate ( f_trad   )  ! radiative temperature of surface [K]
deallocate ( f_tref   )  ! 2 m height air temperature [kelvin]
deallocate ( f_qref   )  ! 2 m height air specific humidity [kg/kg]

!---------------------------------------------------------------------
deallocate ( f_t_soisno    )  ! soil temperature [K]
deallocate ( f_wliq_soisno )  ! liquid water in soil layers [kg/m2]
deallocate ( f_wice_soisno )  ! ice lens in soil layers [kg/m2]
deallocate ( f_h2osoi      )  ! ice lens in soil layers [kg/m2]
deallocate ( f_rstfac      )  ! factor of soil water stress 
deallocate ( f_zwt         )  ! the depth to water table [m]
deallocate ( f_wa          )  ! water storage in aquifer [mm]
deallocate ( f_wat         )  ! total water storage [mm]

deallocate ( f_t_lake       ) ! lake temperature [K]
deallocate ( f_lake_icefrac ) ! lake ice fraction cover [0-1]

!---------------------------------------------------------------------
deallocate ( f_ustar  )  ! u* in similarity theory [m/s]
deallocate ( f_tstar  )  ! t* in similarity theory [kg/kg]
deallocate ( f_qstar  )  ! q* in similarity theory [kg/kg]
deallocate ( f_zol    )  ! dimensionless height (z/L) used in Monin-Obukhov theory
deallocate ( f_rib    )  ! bulk Richardson number in surface layer
deallocate ( f_fm     )  ! integral of profile function for momentum
deallocate ( f_fh     )  ! integral of profile function for heat
deallocate ( f_fq     )  ! integral of profile function for moisture
deallocate ( f_us10m  )  ! 10m u-velocity [m/s]
deallocate ( f_vs10m  )  ! 10m v-velocity [m/s]
deallocate ( f_fm10m  )  ! integral of profile function for momentum at 10m [-]

!---------------------------------------------------------------------
deallocate ( f_xy_us  )  ! wind in eastward direction [m/s]
deallocate ( f_xy_vs  )  ! wind in northward direction [m/s]
deallocate ( f_xy_t   )  ! temperature at reference height [kelvin]
deallocate ( f_xy_q   )  ! specific humidity at reference height [kg/kg]
deallocate ( f_xy_prc )  ! convective precipitation [mm/s]
deallocate ( f_xy_prl )  ! large scale precipitation [mm/s]
deallocate ( f_xy_pbot)  ! atmospheric pressure at the surface [pa]
deallocate ( f_xy_frl )  ! atmospheric infrared (longwave) radiation [W/m2]
deallocate ( f_xy_solarin)  ! downward solar radiation at surface [W/m2]
deallocate ( f_xy_rain)  ! rain [mm/s]
deallocate ( f_xy_snow)  ! snow [mm/s]

END SUBROUTINE deallocate_2D_Fluxes

END MODULE MOD_2D_Fluxes
