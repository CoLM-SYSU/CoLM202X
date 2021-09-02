#include <define.h>

 SUBROUTINE flxwrite (idate,nac,nac_ln,lon_points,lat_points,&
                      dir_output,casename)

!=======================================================================
! Original version: Yongjiu Dai, September 15, 1999, 03/2014
!=======================================================================

  use precision
  USE GlobalVars
  use MOD_2D_Fluxes
  use MOD_TimeInvariants, only: gridlond, gridlatd 
  use timemanager
  IMPLICIT NONE

  integer, INTENT(in) :: idate(3)
  integer, INTENT(in) :: nac
  integer, INTENT(in) :: nac_ln(lon_points,lat_points)
  integer, INTENT(in) :: lon_points
  integer, INTENT(in) :: lat_points
  
  character(LEN=256) :: dir_output
  character(LEN=256) :: casename
 
  integer luout, month, day, i, j, l
  character(LEN=256) fout
  character(LEN=256) cdate
  real(r8) a

! ----------------------------------------------------------------------
! Open for model time varying data (model state variables) and history filed

     luout = 100
#if(defined WO_MONTHLY)
     call julian2monthday(idate(1), idate(2), month, day)
     write(cdate,'(i4.4,"-",i2.2)') idate(1), month
#else
     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1),idate(2),idate(3)
#endif
     fout = trim(dir_output)//trim(casename)//'_'//'2D_Fluxes'//'_'//trim(cdate)
     print*,trim(fout)
     OPEN(unit=luout,file=fout,access='sequential',form='unformatted',&
                     status='unknown',action='write')

     a = float(nac)

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,l)
#endif
     do j = 1, lat_points
        do i = 1, lon_points
           
           if (f_taux   (i,j) /= spval) f_taux   (i,j) = f_taux   (i,j) / a  ! wind stress: E-W [kg/m/s2]
           if (f_tauy   (i,j) /= spval) f_tauy   (i,j) = f_tauy   (i,j) / a  ! wind stress: N-S [kg/m/s2]
           if (f_fsena  (i,j) /= spval) f_fsena  (i,j) = f_fsena  (i,j) / a  ! sensible heat from canopy height to atmosphere [W/m2]
           if (f_lfevpa (i,j) /= spval) f_lfevpa (i,j) = f_lfevpa (i,j) / a  ! latent heat flux from canopy height to atmosphere [W/m2]
           if (f_fevpa  (i,j) /= spval) f_fevpa  (i,j) = f_fevpa  (i,j) / a  ! evapotranspiration from canopy to atmosphere [mm/s]
           if (f_fsenl  (i,j) /= spval) f_fsenl  (i,j) = f_fsenl  (i,j) / a  ! sensible heat from leaves [W/m2]
           if (f_fevpl  (i,j) /= spval) f_fevpl  (i,j) = f_fevpl  (i,j) / a  ! evaporation+transpiration from leaves [mm/s]
           if (f_etr    (i,j) /= spval) f_etr    (i,j) = f_etr    (i,j) / a  ! transpiration rate [mm/s]
           if (f_fseng  (i,j) /= spval) f_fseng  (i,j) = f_fseng  (i,j) / a  ! sensible heat flux from ground [W/m2]
           if (f_fevpg  (i,j) /= spval) f_fevpg  (i,j) = f_fevpg  (i,j) / a  ! evaporation heat flux from ground [mm/s]
           if (f_fgrnd  (i,j) /= spval) f_fgrnd  (i,j) = f_fgrnd  (i,j) / a  ! ground heat flux [W/m2]
           if (f_sabvsun(i,j) /= spval) f_sabvsun(i,j) = f_sabvsun(i,j) / a  ! solar absorbed by sunlit canopy [W/m2]
           if (f_sabvsha(i,j) /= spval) f_sabvsha(i,j) = f_sabvsha(i,j) / a  ! solar absorbed by shaded [W/m2]
           if (f_sabg   (i,j) /= spval) f_sabg   (i,j) = f_sabg   (i,j) / a  ! solar absorbed by ground  [W/m2]
           if (f_olrg   (i,j) /= spval) f_olrg   (i,j) = f_olrg   (i,j) / a  ! outgoing long-wave radiation from ground+canopy [W/m2]
           if (f_rnet   (i,j) /= spval) f_rnet   (i,j) = f_rnet   (i,j) / a  ! net radiation [W/m2]
           if (f_xerr   (i,j) /= spval) f_xerr   (i,j) = f_xerr   (i,j) / a  ! the error of water banace [mm/s]
           if (f_zerr   (i,j) /= spval) f_zerr   (i,j) = f_zerr   (i,j) / a  ! the error of energy balance [W/m2]
           if (f_rsur   (i,j) /= spval) f_rsur   (i,j) = f_rsur   (i,j) / a  ! surface runoff [mm/s]
           if (f_rnof   (i,j) /= spval) f_rnof   (i,j) = f_rnof   (i,j) / a  ! total runoff [mm/s]
           if (f_qintr  (i,j) /= spval) f_qintr  (i,j) = f_qintr  (i,j) / a  ! interception [mm/s]
           if (f_qinfl  (i,j) /= spval) f_qinfl  (i,j) = f_qinfl  (i,j) / a  ! inflitraton [mm/s]
           if (f_qdrip  (i,j) /= spval) f_qdrip  (i,j) = f_qdrip  (i,j) / a  ! throughfall [mm/s]
           if (f_rstfac (i,j) /= spval) f_rstfac (i,j) = f_rstfac (i,j) / a  ! factor of soil water stress 
           if (f_zwt    (i,j) /= spval) f_zwt    (i,j) = f_zwt    (i,j) / a  ! water depth [m]
           if (f_wa     (i,j) /= spval) f_wa     (i,j) = f_wa     (i,j) / a  ! water storage in aquifer [mm]
           if (f_wat    (i,j) /= spval) f_wat    (i,j) = f_wat    (i,j) / a  ! total water storage [mm]
           if (f_assim  (i,j) /= spval) f_assim  (i,j) = f_assim  (i,j) / a  ! canopy assimilation rate [mol m-2 s-1]
           if (f_respc  (i,j) /= spval) f_respc  (i,j) = f_respc  (i,j) / a  ! respiration (plant+soil) [mol m-2 s-1]
           if (f_qcharge(i,j) /= spval) f_qcharge(i,j) = f_qcharge(i,j) / a  ! groundwater recharge rate [mm/s] 

!---------------------------------------------------------------------
           if (f_t_grnd (i,j) /= spval) f_t_grnd (i,j) = f_t_grnd (i,j) / a  ! ground surface temperature [K]
           if (f_tleaf  (i,j) /= spval) f_tleaf  (i,j) = f_tleaf  (i,j) / a  ! sunlit leaf temperature [K]
           if (f_ldew   (i,j) /= spval) f_ldew   (i,j) = f_ldew   (i,j) / a  ! depth of water on foliage [mm]
           if (f_scv    (i,j) /= spval) f_scv    (i,j) = f_scv    (i,j) / a  ! snow cover, water equivalent [mm]
           if (f_snowdp (i,j) /= spval) f_snowdp (i,j) = f_snowdp (i,j) / a  ! snow depth [meter]
           if (f_fsno   (i,j) /= spval) f_fsno   (i,j) = f_fsno   (i,j) / a  ! fraction of snow cover on ground
           if (f_sigf   (i,j) /= spval) f_sigf   (i,j) = f_sigf   (i,j) / a  ! fraction of veg cover, excluding snow-covered veg [-]
           if (f_green  (i,j) /= spval) f_green  (i,j) = f_green  (i,j) / a  ! leaf greenness
           if (f_lai    (i,j) /= spval) f_lai    (i,j) = f_lai    (i,j) / a  ! leaf area index
           if (f_laisun (i,j) /= spval) f_laisun (i,j) = f_laisun (i,j) / a  ! leaf area index
           if (f_laisha (i,j) /= spval) f_laisha (i,j) = f_laisha (i,j) / a  ! leaf area index
           if (f_sai    (i,j) /= spval) f_sai    (i,j) = f_sai    (i,j) / a  ! stem area index
           if (f_alb(1,1,i,j) /= spval) f_alb(1,1,i,j) = f_alb(1,1,i,j) / a  ! averaged albedo [visible, direct; direct, diffuse]
           if (f_alb(2,1,i,j) /= spval) f_alb(2,1,i,j) = f_alb(2,1,i,j) / a  ! averaged albedo [visible, direct; direct, diffuse]
           if (f_alb(1,2,i,j) /= spval) f_alb(1,2,i,j) = f_alb(1,2,i,j) / a  ! averaged albedo [visible, direct; direct, diffuse]
           if (f_alb(2,2,i,j) /= spval) f_alb(2,2,i,j) = f_alb(2,2,i,j) / a  ! averaged albedo [visible, direct; direct, diffuse]
           if (f_emis   (i,j) /= spval) f_emis   (i,j) = f_emis   (i,j) / a  ! averaged bulk surface emissivity
           if (f_z0m    (i,j) /= spval) f_z0m    (i,j) = f_z0m    (i,j) / a  ! effective roughness [m]
           if (f_trad   (i,j) /= spval) f_trad   (i,j) = f_trad   (i,j) / a  ! radiative temperature of surface [K]
           if (f_tref   (i,j) /= spval) f_tref   (i,j) = f_tref   (i,j) / a  ! 2 m height air temperature [kelvin]
           if (f_qref   (i,j) /= spval) f_qref   (i,j) = f_qref   (i,j) / a  ! 2 m height air specific humidity [kg/kg]
           if (f_xy_rain(i,j) /= spval) f_xy_rain(i,j) = f_xy_rain(i,j) / a  ! rain [mm/s]
           if (f_xy_snow(i,j) /= spval) f_xy_snow(i,j) = f_xy_snow(i,j) / a  ! snow [mm/s]

!---------------------------------------------------------------------
           do l = maxsnl+1, nl_soil
              if (f_t_soisno   (l,i,j) /= spval) f_t_soisno   (l,i,j) = f_t_soisno   (l,i,j) / a  ! soil temperature [K]
              if (f_wliq_soisno(l,i,j) /= spval) f_wliq_soisno(l,i,j) = f_wliq_soisno(l,i,j) / a  ! liquid water in soil layers [kg/m2]
              if (f_wice_soisno(l,i,j) /= spval) f_wice_soisno(l,i,j) = f_wice_soisno(l,i,j) / a  ! ice lens in soil layers [kg/m2]
           end do

           do l = 1, nl_soil
              if (f_h2osoi     (l,i,j) /= spval) f_h2osoi     (l,i,j) = f_h2osoi     (l,i,j) / a  ! volumetric soil water in layers [m3/m3]
           end do

           do l = 1, nl_lake
              if (f_t_lake      (l,i,j) /= spval) f_t_lake      (l,i,j) = f_t_lake      (l,i,j) / a  ! lake temperature [K]
              if (f_lake_icefrac(l,i,j) /= spval) f_lake_icefrac(l,i,j) = f_lake_icefrac(l,i,j) / a  ! lake ice fraction cover [0-1]
           end do

           if (f_ustar  (i,j) /= spval) f_ustar  (i,j) = f_ustar  (i,j) / a  ! u* in similarity theory [m/s]
           if (f_tstar  (i,j) /= spval) f_tstar  (i,j) = f_tstar  (i,j) / a  ! t* in similarity theory [kg/kg]
           if (f_qstar  (i,j) /= spval) f_qstar  (i,j) = f_qstar  (i,j) / a  ! q* in similarity theory [kg/kg]
           if (f_zol    (i,j) /= spval) f_zol    (i,j) = f_zol    (i,j) / a  ! dimensionless height (z/L) used in Monin-Obukhov theory
           if (f_rib    (i,j) /= spval) f_rib    (i,j) = f_rib    (i,j) / a  ! bulk Richardson number in surface layer
           if (f_fm     (i,j) /= spval) f_fm     (i,j) = f_fm     (i,j) / a  ! integral of profile function for momentum
           if (f_fh     (i,j) /= spval) f_fh     (i,j) = f_fh     (i,j) / a  ! integral of profile function for heat
           if (f_fq     (i,j) /= spval) f_fq     (i,j) = f_fq     (i,j) / a  ! integral of profile function for moisture
           if (f_us10m  (i,j) /= spval) f_us10m  (i,j) = f_us10m  (i,j) / a  ! 10m u-velocity [m/s]
           if (f_vs10m  (i,j) /= spval) f_vs10m  (i,j) = f_vs10m  (i,j) / a  ! 10m v-velocity [m/s]
           if (f_fm10m  (i,j) /= spval) f_fm10m  (i,j) = f_fm10m  (i,j) / a  ! integral of profile function for momentum at 10m [-]

           if (f_xy_us  (i,j) /= spval) f_xy_us  (i,j) = f_xy_us  (i,j) / a  ! wind in eastward direction [m/s]
           if (f_xy_vs  (i,j) /= spval) f_xy_vs  (i,j) = f_xy_vs  (i,j) / a  ! wind in northward direction [m/s]
           if (f_xy_t   (i,j) /= spval) f_xy_t   (i,j) = f_xy_t   (i,j) / a  ! temperature at reference height [kelvin]
           if (f_xy_q   (i,j) /= spval) f_xy_q   (i,j) = f_xy_q   (i,j) / a  ! specific humidity at reference height [kg/kg]
           if (f_xy_prc (i,j) /= spval) f_xy_prc (i,j) = f_xy_prc (i,j) / a  ! convective precipitation [mm/s]
           if (f_xy_prl (i,j) /= spval) f_xy_prl (i,j) = f_xy_prl (i,j) / a  ! large scale precipitation [mm/s]
           if (f_xy_pbot(i,j) /= spval) f_xy_pbot(i,j) = f_xy_pbot(i,j) / a  ! atmospheric pressure at the surface [pa]
           if (f_xy_frl (i,j) /= spval) f_xy_frl (i,j) = f_xy_frl (i,j) / a  ! atmospheric infrared (longwave) radiation [W/m2]
           if (f_xy_solarin(i,j) /= spval) f_xy_solarin(i,j) = f_xy_solarin(i,j) / a  ! downward solar radiation at surface [W/m2]

           if (f_sr     (i,j) /= spval) f_sr     (i,j) = f_sr     (i,j) / a  ! total reflected solar radiation (W/m2)
           if (f_solvd  (i,j) /= spval) f_solvd  (i,j) = f_solvd  (i,j) / a  ! incident direct beam vis solar radiation (W/m2)
           if (f_solvi  (i,j) /= spval) f_solvi  (i,j) = f_solvi  (i,j) / a  ! incident diffuse beam vis solar radiation (W/m2)
           if (f_solnd  (i,j) /= spval) f_solnd  (i,j) = f_solnd  (i,j) / a  ! incident direct beam nir solar radiation (W/m2)
           if (f_solni  (i,j) /= spval) f_solni  (i,j) = f_solni  (i,j) / a  ! incident diffuse beam nir solar radiation (W/m2)
           if (f_srvd   (i,j) /= spval) f_srvd   (i,j) = f_srvd   (i,j) / a  ! reflected direct beam vis solar radiation (W/m2)
           if (f_srvi   (i,j) /= spval) f_srvi   (i,j) = f_srvi   (i,j) / a  ! reflected diffuse beam vis solar radiation (W/m2)
           if (f_srnd   (i,j) /= spval) f_srnd   (i,j) = f_srnd   (i,j) / a  ! reflected direct beam nir solar radiation (W/m2)
           if (f_srni   (i,j) /= spval) f_srni   (i,j) = f_srni   (i,j) / a  ! reflected diffuse beam nir solar radiation (W/m2)
           if (f_solvdln(i,j) /= spval) f_solvdln(i,j) = f_solvdln(i,j) / nac_ln(i,j) ! incident direct beam vis solar radiation at local noon (W/m2)
           if (f_solviln(i,j) /= spval) f_solviln(i,j) = f_solviln(i,j) / nac_ln(i,j) ! incident diffuse beam vis solar radiation at local noon (W/m2)
           if (f_solndln(i,j) /= spval) f_solndln(i,j) = f_solndln(i,j) / nac_ln(i,j) ! incident direct beam nir solar radiation at local noon (W/m2)
           if (f_solniln(i,j) /= spval) f_solniln(i,j) = f_solniln(i,j) / nac_ln(i,j) ! incident diffuse beam nir solar radiation at local noon (W/m2)
           if (f_srvdln (i,j) /= spval) f_srvdln (i,j) = f_srvdln (i,j) / nac_ln(i,j) ! reflected direct beam vis solar radiation at local noon (W/m2)
           if (f_srviln (i,j) /= spval) f_srviln (i,j) = f_srviln (i,j) / nac_ln(i,j) ! reflected diffuse beam vis solar radiation at local noon (W/m2)
           if (f_srndln (i,j) /= spval) f_srndln (i,j) = f_srndln (i,j) / nac_ln(i,j) ! reflected direct beam nir solar radiation at local noon (W/m2)
           if (f_srniln (i,j) /= spval) f_srniln (i,j) = f_srniln (i,j) / nac_ln(i,j) ! reflected diffuse beam nir solar radiation at local noon(W/m2)

        end do
     end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

     write(luout) gridlond (:)    ! longitude in degree
     write(luout) gridlatd (:)    ! latitude in degree
     write(luout) mask     (:,:)  ! grid mask
     write(luout) frac     (:,:)  ! grid total fraction
     write(luout) area     (:,:)  ! grid cell area

     write(luout) f_taux   (:,:)  ! wind stress: E-W [kg/m/s2]
     write(luout) f_tauy   (:,:)  ! wind stress: N-S [kg/m/s2]
     write(luout) f_fsena  (:,:)  ! sensible heat from canopy height to atmosphere [W/m2]
     write(luout) f_lfevpa (:,:)  ! latent heat flux from canopy height to atmosphere [W/m2]
     write(luout) f_fevpa  (:,:)  ! evapotranspiration from canopy to atmosphere [mm/s]
     write(luout) f_fsenl  (:,:)  ! sensible heat from leaves [W/m2]
     write(luout) f_fevpl  (:,:)  ! evaporation+transpiration from leaves [mm/s]
     write(luout) f_etr    (:,:)  ! transpiration rate [mm/s]
     write(luout) f_fseng  (:,:)  ! sensible heat flux from ground [W/m2]
     write(luout) f_fevpg  (:,:)  ! evaporation heat flux from ground [mm/s]
     write(luout) f_fgrnd  (:,:)  ! ground heat flux [W/m2]
     write(luout) f_sabvsun(:,:)  ! solar absorbed by sunlit canopy [W/m2]
     write(luout) f_sabvsha(:,:)  ! solar absorbed by shaded [W/m2]
     write(luout) f_sabg   (:,:)  ! solar absorbed by ground  [W/m2]
     write(luout) f_olrg   (:,:)  ! outgoing long-wave radiation from ground+canopy [W/m2]
     write(luout) f_rnet   (:,:)  ! net radiation [W/m2]
     write(luout) f_xerr   (:,:)  ! the error of water banace [mm/s]
     write(luout) f_zerr   (:,:)  ! the error of energy balance [W/m2]
     write(luout) f_rsur   (:,:)  ! surface runoff [mm/s]
     write(luout) f_rnof   (:,:)  ! total runoff [mm/s]
     write(luout) f_qintr  (:,:)  ! interception [mm/s]
     write(luout) f_qinfl  (:,:)  ! inflitration [mm/s]
     write(luout) f_qdrip  (:,:)  ! throughfall [mm/s]
     write(luout) f_assim  (:,:)  ! canopy assimilation rate [mol m-2 s-1]
     write(luout) f_respc  (:,:)  ! respiration (plant+soil) [mol m-2 s-1]
     write(luout) f_qcharge(:,:)  ! groundwater recharge rate [mm/s] 

!---------------------------------------------------------------------
     write(luout) f_t_grnd (:,:)  ! ground surface temperature [K]
     write(luout) f_tleaf  (:,:)  ! sunlit leaf temperature [K]
     write(luout) f_ldew   (:,:)  ! depth of water on foliage [mm]
     write(luout) f_scv    (:,:)  ! snow cover, water equivalent [mm]
     write(luout) f_snowdp (:,:)  ! snow depth [meter]
     write(luout) f_fsno   (:,:)  ! fraction of snow cover on ground
     write(luout) f_sigf   (:,:)  ! fraction of veg cover, excluding snow-covered veg [-]
     write(luout) f_green  (:,:)  ! leaf greenness
     write(luout) f_lai    (:,:)  ! leaf area index
     write(luout) f_laisun (:,:)  ! sunlit leaf area index
     write(luout) f_laisha (:,:)  ! shaded leaf area index
     write(luout) f_sai    (:,:)  ! stem area index
     write(luout) f_alb(:,:,:,:)  ! averaged albedo [visible, direct; direct, diffuse]
     write(luout) f_emis   (:,:)  ! averaged bulk surface emissivity
     write(luout) f_z0m    (:,:)  ! effective roughness [m]
     write(luout) f_trad   (:,:)  ! radiative temperature of surface [K]
     write(luout) f_tref   (:,:)  ! 2 m height air temperature [kelvin]
     write(luout) f_qref   (:,:)  ! 2 m height air specific humidity [kg/kg]
     write(luout) f_xy_rain(:,:)  ! rain [mm/s]
     write(luout) f_xy_snow(:,:)  ! snow [mm/s]
 
!---------------------------------------------------------------------
     write(luout) f_t_soisno   (:,:,:)  ! soil temperature [K]
     write(luout) f_wliq_soisno(:,:,:)  ! liquid water in soil layers [kg/m2]
     write(luout) f_wice_soisno(:,:,:)  ! ice lens in soil layers [kg/m2]
     write(luout) f_h2osoi     (:,:,:)  ! volumetric soil water in layers [m3/m3]
     write(luout) f_rstfac     (:,:)    ! factor of soil water stress 
     write(luout) f_zwt        (:,:)    ! the depth to water table [m]
     write(luout) f_wa         (:,:)    ! water storage in aquifer [mm]
     write(luout) f_wat        (:,:)    ! total water storage [mm]
 
     write(luout) f_t_lake      (:,:,:) ! lake temperature [K]
     write(luout) f_lake_icefrac(:,:,:) ! lake ice fraction cover [0-1]

     write(luout) f_ustar  (:,:)   ! u* in similarity theory [m/s]
     write(luout) f_tstar  (:,:)   ! t* in similarity theory [kg/kg]
     write(luout) f_qstar  (:,:)   ! q* in similarity theory [kg/kg]
     write(luout) f_zol    (:,:)   ! dimensionless height (z/L) used in Monin-Obukhov theory
     write(luout) f_rib    (:,:)   ! bulk Richardson number in surface layer
     write(luout) f_fm     (:,:)   ! integral of profile function for momentum
     write(luout) f_fh     (:,:)   ! integral of profile function for heat
     write(luout) f_fq     (:,:)   ! integral of profile function for moisture
     write(luout) f_us10m  (:,:)   ! 10m u-velocity [m/s]
     write(luout) f_vs10m  (:,:)   ! 10m v-velocity [m/s]
     write(luout) f_fm10m  (:,:)   ! integral of profile function for momentum at 10m [-]

     write(luout) f_xy_us  (:,:)   ! wind in eastward direction [m/s]
     write(luout) f_xy_vs  (:,:)   ! wind in northward direction [m/s]
     write(luout) f_xy_t   (:,:)   ! temperature at reference height [kelvin]
     write(luout) f_xy_q   (:,:)   ! specific humidity at reference height [kg/kg]
     write(luout) f_xy_prc (:,:)   ! convective precipitation [mm/s]
     write(luout) f_xy_prl (:,:)   ! large scale precipitation [mm/s]
     write(luout) f_xy_pbot(:,:)   ! atmospheric pressure at the surface [pa]
     write(luout) f_xy_frl (:,:)   ! atmospheric infrared (longwave) radiation [W/m2]
     write(luout) f_xy_solarin(:,:)   ! downward solar radiation at surface [W/m2]

     write(luout) f_sr     (:,:)  ! total reflected solar radiation (W/m2)
     write(luout) f_solvd  (:,:)  ! incident direct beam vis solar radiation (W/m2)
     write(luout) f_solvi  (:,:)  ! incident diffuse beam vis solar radiation (W/m2)
     write(luout) f_solnd  (:,:)  ! incident direct beam nir solar radiation (W/m2)
     write(luout) f_solni  (:,:)  ! incident diffuse beam nir solar radiation (W/m2)
     write(luout) f_srvd   (:,:)  ! reflected direct beam vis solar radiation (W/m2)
     write(luout) f_srvi   (:,:)  ! reflected diffuse beam vis solar radiation (W/m2)
     write(luout) f_srnd   (:,:)  ! reflected direct beam nir solar radiation (W/m2)
     write(luout) f_srni   (:,:)  ! reflected diffuse beam nir solar radiation (W/m2)
     write(luout) f_solvdln(:,:)  ! incident direct beam vis solar radiation at local noon(W/m2)
     write(luout) f_solviln(:,:)  ! incident diffuse beam vis solar radiation at local noon(W/m2)
     write(luout) f_solndln(:,:)  ! incident direct beam nir solar radiation at local noon(W/m2)
     write(luout) f_solniln(:,:)  ! incident diffuse beam nir solar radiation at local noon(W/m2)
     write(luout) f_srvdln (:,:)  ! reflected direct beam vis solar radiation at local noon(W/m2)
     write(luout) f_srviln (:,:)  ! reflected diffuse beam vis solar radiation at local noon(W/m2)
     write(luout) f_srndln (:,:)  ! reflected direct beam nir solar radiation at local noon(W/m2)
     write(luout) f_srniln (:,:)  ! reflected diffuse beam nir solar radiation at local noon(W/m2)

     CLOSE (luout)

  END SUBROUTINE flxwrite

! ----------------------------------------------------------------------
! EOP
