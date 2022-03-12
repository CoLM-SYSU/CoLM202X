#include <define.h>

program hist_concatenate

   USE spmd_task
   USE GlobalVars
   use mod_namelist
   use mod_block
   use mod_grid
   use mod_concatenate
   USE ncio_serial

   implicit none

   ! Local variables 
   TYPE(grid_type) :: ghist
   character(len=256) :: nlfile, filename
   character(len=256) :: filehist
   INTEGER :: timelen
   integer, parameter :: compress = 1

#ifdef USEMPI
   CALL spmd_init
#endif

   IF (p_is_master) THEN
      call getarg (1, nlfile)
      call getarg (2, filehist)
   ENDIF

   call read_namelist (nlfile)

   IF (p_is_master) THEN
      ! Load block information.
      filename = trim(DEF_dir_landdata) // '/block.nc'
      CALL ncio_read_serial (filename, 'lat_s', gblock%lat_s)
      CALL ncio_read_serial (filename, 'lat_n', gblock%lat_n)
      CALL ncio_read_serial (filename, 'lon_w', gblock%lon_w)
      CALL ncio_read_serial (filename, 'lon_e', gblock%lon_e)
      gblock%nxblk = size(gblock%lon_w)
      gblock%nyblk = size(gblock%lat_s)

      ! Define grid of history data.
      call ghist%define_by_ndims (DEF_nlon_hist, DEF_nlat_hist)
         
      call set_hist_block_info (ghist, hist_block_info)

      CALL ncio_create_file (filehist)

      call hist_concatenate_time      (filehist, 'time', timelen)         
      CALL hist_concatenate_grid_info (filehist)

      CALL ncio_define_dimension (filehist, 'soil', nl_soil)
      CALL ncio_define_dimension (filehist, 'lake', nl_lake)
      CALL ncio_define_dimension (filehist, 'soilsnow', nl_soil-maxsnl)
      CALL ncio_define_dimension (filehist, 'band', 2)
      CALL ncio_define_dimension (filehist, 'wetdry', 2)
#ifdef PLANT_HYDRAULIC_STRESS
      CALL ncio_define_dimension (filehist, 'vegnodes', nvegwcs)
#endif

   ENDIF

   call hist_concatenate_var_2d (filehist, 'f_taux   ', timelen, compress, &  
      'wind stress: E-W','kg/m/s2')

   call hist_concatenate_var_2d (filehist, 'f_tauy   ', timelen, compress, &  
      'wind stress: N-S','kg/m/s2')

   call hist_concatenate_var_2d (filehist, 'f_fsena  ', timelen, compress, &  
      'sensible heat from canopy height to atmosphere','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_lfevpa ', timelen, compress, &  
      'latent heat flux from canopy height to atmosphere','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_fevpa  ', timelen, compress, &  
      'evapotranspiration from canopy height to atmosphere','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_fsenl  ', timelen, compress, &  
      'sensible heat from leaves','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_fevpl  ', timelen, compress, &  
      'evaporation+transpiration from leaves','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_etr    ', timelen, compress, &  
      'transpiration rate','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_fseng  ', timelen, compress, &  
      'sensible heat flux from ground','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_fevpg  ', timelen, compress, &  
      'evaporation heat flux from ground','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_fgrnd  ', timelen, compress, &  
      'ground heat flux','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_sabvsun', timelen, compress, &  
      'solar absorbed by sunlit canopy','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_sabvsha', timelen, compress, &  
      'solar absorbed by shaded','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_sabg   ', timelen, compress, &  
      'solar absorbed by ground','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_olrg   ', timelen, compress, &  
      'outgoing long-wave radiation from ground+canopy','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_rnet   ', timelen, compress, &  
      'net radiation','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_xerr   ', timelen, compress, &  
      'the error of water banace','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_zerr   ', timelen, compress, &  
      'the error of energy balance','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_rsur   ', timelen, compress, &  
      'surface runoff','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_rnof   ', timelen, compress, &  
      'total runoff','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_qintr  ', timelen, compress, &  
      'interception','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_qinfl  ', timelen, compress, &  
      'f_qinfl','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_qdrip  ', timelen, compress, &  
      'total throughfall','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_assim  ', timelen, compress, &  
      'canopy assimilation rate','umol m-2 s-1')

   call hist_concatenate_var_2d (filehist, 'f_respc  ', timelen, compress, &  
      'respiration (plant+soil)','mol m-2 s-1')

   call hist_concatenate_var_2d (filehist, 'f_qcharge', timelen, compress, &
      'groundwater recharge rate','mm/s')

   !-------------------------------------------------------
   call hist_concatenate_var_2d (filehist, 'f_t_grnd ', timelen, compress, & 
      'ground surface temperature','K')

   call hist_concatenate_var_2d (filehist, 'f_tleaf  ', timelen, compress, & 
      'leaf temperature','K')

   call hist_concatenate_var_2d (filehist, 'f_ldew   ', timelen, compress, & 
      'depth of water on foliage','mm')

   call hist_concatenate_var_2d (filehist, 'f_scv    ', timelen, compress, & 
      'snow cover, water equivalent','mm')

   call hist_concatenate_var_2d (filehist, 'f_snowdp ', timelen, compress, &
      'snow depth','meter')

   call hist_concatenate_var_2d (filehist, 'f_fsno   ', timelen, compress, &
      'fraction of snow cover on ground','-')

   call hist_concatenate_var_2d (filehist, 'f_sigf   ', timelen, compress, & 
      'fraction of veg cover, excluding snow-covered veg','-')

   call hist_concatenate_var_2d (filehist, 'f_green  ', timelen, compress, & 
      'leaf greenness','-')

   call hist_concatenate_var_2d (filehist, 'f_lai    ', timelen, compress, & 
      'leaf area index','m2/m2')

   call hist_concatenate_var_2d (filehist, 'f_laisun ', timelen, compress, & 
      'sunlit leaf area index','m2/m2')

   call hist_concatenate_var_2d (filehist, 'f_laisha ', timelen, compress, & 
      'shaded leaf area index','m2/m2')

   call hist_concatenate_var_2d (filehist, 'f_sai    ', timelen, compress, & 
      'stem area index','m2/m2')

   call hist_concatenate_var_4d (filehist, 'f_alb    ', timelen, &
      'band', 'wetdry', 2, 2, compress, 'averaged albedo direct','%')

   call hist_concatenate_var_2d (filehist, 'f_emis   ', timelen, compress, &  
      'averaged bulk surface emissivity','-')

   call hist_concatenate_var_2d (filehist, 'f_z0m    ', timelen, compress, &  
      'effective roughness','m')

   call hist_concatenate_var_2d (filehist, 'f_trad   ', timelen, compress, &  
      'radiative temperature of surface','kelvin')

   call hist_concatenate_var_2d (filehist, 'f_tref   ', timelen, compress, &  
      '2 m height air temperature','kelvin')

   call hist_concatenate_var_2d (filehist, 'f_qref   ', timelen, compress, &  
      '2 m height air specific humidity','kg/kg')

   call hist_concatenate_var_2d (filehist, 'f_xy_rain', timelen, compress, &  
      'rain','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_xy_snow', timelen, compress, &  
      'snow','mm/s')

   !--------------------------------------------------------------
   call hist_concatenate_var_3d (filehist, 'f_t_soisno   ', timelen, 'soilsnow', &
      nl_soil-maxsnl, compress, 'soil temperature','K')

   call hist_concatenate_var_3d (filehist, 'f_wliq_soisno', timelen, 'soilsnow', &
      nl_soil-maxsnl, compress, 'liquid water in soil layers','kg/m2')

   call hist_concatenate_var_3d (filehist, 'f_wice_soisno', timelen, 'soilsnow', &
      nl_soil-maxsnl, compress, 'ice lens in soil layers','kg/m2')

   call hist_concatenate_var_3d (filehist, 'f_h2osoi     ', timelen, 'soil', &
      nl_soil, compress, 'volumetric water in soil layers','m3/m3')

   ! call hist_concatenate_var_2d (filehist, 'f_rstfacsun  ', timelen, compress)  ! factor of soil water stress to transpiration on sunlit leaf
   ! call hist_concatenate_var_2d (filehist, 'f_rstfacsha  ', timelen, compress)  ! factor of soil water stress to transpiration on shaded leaf
   
   ! fraction of root water uptake from each soil layer (PHS undefined)
   ! water exchange between soil layers and root. Positive: soil->root [mm h2o/s] (PHS defined)
   call hist_concatenate_var_3d (filehist, 'f_rootr      ', timelen, 'soil', nl_soil, compress)  

#ifdef PLANT_HYDRAULIC_STRESS
   call hist_concatenate_var_3d (filehist, 'f_vegwp      ', timelen, 'vegnodes', nvegwcs, compress)  ! vegetation water potential [mm]
#endif
   
   call hist_concatenate_var_2d (filehist, 'f_zwt        ', timelen, compress, & 
      'the depth to water table','m')

   call hist_concatenate_var_2d (filehist, 'f_wa         ', timelen, compress, &
      'water storage in aquifer','mm')

   call hist_concatenate_var_2d (filehist, 'f_wat        ', timelen, compress, &
      'total water storage','mm')

   call hist_concatenate_var_3d (filehist, 'f_t_lake      ', timelen, 'lake', &
      nl_lake, compress, 'lake temperature','K')

   call hist_concatenate_var_3d (filehist, 'f_lake_icefrac', timelen, 'lake', &
      nl_lake, compress, 'lake ice fraction cover','0-1')

   call hist_concatenate_var_2d (filehist, 'f_ustar  ', timelen, compress, & 
      'u* in similarity theory','m/s')

   call hist_concatenate_var_2d (filehist, 'f_tstar  ', timelen, compress, & 
      't* in similarity theory','kg/kg')

   call hist_concatenate_var_2d (filehist, 'f_qstar  ', timelen, compress, & 
      'q* in similarity theory')

   call hist_concatenate_var_2d (filehist, 'f_zol    ', timelen, compress, & 
      'dimensionless height (z/L) used in Monin-Obukhov theory','-')

   call hist_concatenate_var_2d (filehist, 'f_rib    ', timelen, compress, & 
      'bulk Richardson number in surface layer','-')

   call hist_concatenate_var_2d (filehist, 'f_fm     ', timelen, compress, & 
      'integral of profile function for momentum','-')

   call hist_concatenate_var_2d (filehist, 'f_fh     ', timelen, compress, & 
      'integral of profile function for heat','-')

   call hist_concatenate_var_2d (filehist, 'f_fq     ', timelen, compress, & 
      'integral of profile function for moisture','-')

   call hist_concatenate_var_2d (filehist, 'f_us10m  ', timelen, compress, & 
      '10m u-velocity','m/s')

   call hist_concatenate_var_2d (filehist, 'f_vs10m  ', timelen, compress, & 
      '10m v-velocity','m/s')

   call hist_concatenate_var_2d (filehist, 'f_fm10m  ', timelen, compress, &  
      'integral of profile function for momentum at 10m','-')

   call hist_concatenate_var_2d (filehist, 'f_xy_us  ', timelen, compress, &  
      'wind in eastward direction', 'm/s')

   call hist_concatenate_var_2d (filehist, 'f_xy_vs  ', timelen, compress, &
      'wind in northward direction','m/s')

   call hist_concatenate_var_2d (filehist, 'f_xy_t   ', timelen, compress, &
      'temperature at reference height','kelvin')

   call hist_concatenate_var_2d (filehist, 'f_xy_q   ', timelen, compress, & 
      'specific humidity at reference height','kg/kg')

   call hist_concatenate_var_2d (filehist, 'f_xy_prc ', timelen, compress, &  
      'convective precipitation','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_xy_prl ', timelen, compress, &  
      'large scale precipitation','mm/s')

   call hist_concatenate_var_2d (filehist, 'f_xy_pbot', timelen, compress, &  
      'atmospheric pressure at the surface','pa')

   call hist_concatenate_var_2d (filehist, 'f_xy_frl ', timelen, compress, &  
      'atmospheric infrared (longwave) radiation','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_xy_solarin', timelen, compress, & 
      'downward solar radiation at surface','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_sr     ', timelen, compress, & 
      'reflected solar radiation at surface [W/m2]','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_solvd  ', timelen, compress, & 
      'incident direct beam vis solar radiation (W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_solvi  ', timelen, compress, & 
      'incident diffuse beam vis solar radiation (W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_solnd  ', timelen, compress, & 
      'incident direct beam nir solar radiation (W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_solni  ', timelen, compress, & 
      'incident diffuse beam nir solar radiation (W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_srvd   ', timelen, compress, & 
      'reflected direct beam vis solar radiation (W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_srvi   ', timelen, compress, & 
      'reflected diffuse beam vis solar radiation (W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_srnd   ', timelen, compress, & 
      'reflected direct beam nir solar radiation (W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_srni   ', timelen, compress, & 
      'reflected diffuse beam nir solar radiation (W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_solvdln', timelen, compress, &
      'incident direct beam vis solar radiation at local noon(W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_solviln', timelen, compress, & 
      'incident diffuse beam vis solar radiation at local noon(W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_solndln', timelen, compress, & 
      'incident direct beam nir solar radiation at local noon(W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_solniln', timelen, compress, & 
      'incident diffuse beam nir solar radiation at local noon(W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_srvdln ', timelen, compress, & 
      'reflected direct beam vis solar radiation at local noon(W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_srviln ', timelen, compress, & 
      'reflected diffuse beam vis solar radiation at local noon(W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_srndln ', timelen, compress, & 
      'reflected direct beam nir solar radiation at local noon(W/m2)','W/m2')

   call hist_concatenate_var_2d (filehist, 'f_srniln ', timelen, compress, &  
      'reflected diffuse beam nir solar radiation at local noon(W/m2)','W/m2')

end program hist_concatenate
! ----------------------------------------------------------------------
! EOP
