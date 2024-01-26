#include <define.h>

PROGRAM hist_concatenate

   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_Grid
   USE MOD_Concatenate
   USE MOD_NetCDFSerial

   IMPLICIT NONE

   ! Local variables
   type(grid_type) :: ghist
   integer :: lon_points, lat_points
   character(len=256) :: nlfile, filename
   character(len=256) :: filehist
   integer :: timelen
   integer, parameter :: compress = 1

#ifdef USEMPI
      CALL spmd_init
#endif

      IF (p_is_master) THEN
         write(*,*) 'Usage : PATH/hist_concatenate colm_namelist_file filehist'
         write(*,*) '      : < filehist is the file of hist without block information. >'
      ENDIF
   
      IF (p_is_master) THEN
         CALL getarg (1, nlfile)
         CALL getarg (2, filehist)
      ENDIF
   
      CALL read_namelist (nlfile)
   
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
         lon_points = nint(360.0/DEF_hist_lon_res)
         lat_points = nint(180.0/DEF_hist_lat_res)
   
         CALL ghist%define_by_ndims (lon_points, lat_points)
   
         CALL set_hist_block_info (ghist, hist_block_info)
   
         CALL ncio_create_file (filehist)
   
         CALL hist_concatenate_time      (filehist, 'time', timelen)
         CALL hist_concatenate_grid_info (filehist)
   
         CALL ncio_define_dimension (filehist, 'soil', nl_soil)
         CALL ncio_define_dimension (filehist, 'lake', nl_lake)
         CALL ncio_define_dimension (filehist, 'soilsnow', nl_soil-maxsnl)
         CALL ncio_define_dimension (filehist, 'band', 2)
         CALL ncio_define_dimension (filehist, 'rtyp', 2)
         IF(DEF_USE_PLANTHYDRAULICS)THEN
            CALL ncio_define_dimension (filehist, 'vegnodes', nvegwcs)
         ENDIF
   
      ENDIF
   
      CALL hist_concatenate_var_2d (filehist, 'f_taux   ', timelen, compress, &
         'wind stress: E-W','kg/m/s2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_tauy   ', timelen, compress, &
         'wind stress: N-S','kg/m/s2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fsena  ', timelen, compress, &
         'sensible heat from canopy height to atmosphere','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_lfevpa ', timelen, compress, &
         'latent heat flux from canopy height to atmosphere','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fevpa  ', timelen, compress, &
         'evapotranspiration from canopy height to atmosphere','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fevpg  ', timelen, compress, &
         'evaporation heat flux from ground [mm/s]')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fsenl  ', timelen, compress, &
         'sensible heat from leaves','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fevpl  ', timelen, compress, &
         'evaporation+transpiration from leaves','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_etr    ', timelen, compress, &
         'transpiration rate','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fseng  ', timelen, compress, &
         'sensible heat flux from ground','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fevpg  ', timelen, compress, &
         'evaporation heat flux from ground','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fgrnd  ', timelen, compress, &
         'ground heat flux','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_sabvsun', timelen, compress, &
         'solar absorbed by sunlit canopy','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_sabvsha', timelen, compress, &
         'solar absorbed by shaded','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_sabg   ', timelen, compress, &
         'solar absorbed by ground','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_olrg   ', timelen, compress, &
         'outgoing long-wave radiation from ground+canopy','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_rnet   ', timelen, compress, &
         'net radiation','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xerr   ', timelen, compress, &
         'the error of water banace','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_zerr   ', timelen, compress, &
         'the error of energy balance','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_rsur   ', timelen, compress, &
         'surface runoff','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_rnof   ', timelen, compress, &
         'total runoff','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_qintr  ', timelen, compress, &
         'interception','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_qinfl  ', timelen, compress, &
         'f_qinfl','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_qdrip  ', timelen, compress, &
         'total throughfall','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_assim  ', timelen, compress, &
         'canopy assimilation rate','umol m-2 s-1')
   
      CALL hist_concatenate_var_2d (filehist, 'f_respc  ', timelen, compress, &
         'respiration (plant+soil)','mol m-2 s-1')
   
      CALL hist_concatenate_var_2d (filehist, 'f_qcharge', timelen, compress, &
         'groundwater recharge rate','mm/s')
   
      !-------------------------------------------------------
      CALL hist_concatenate_var_2d (filehist, 'f_t_grnd ', timelen, compress, &
         'ground surface temperature','K')
   
      CALL hist_concatenate_var_2d (filehist, 'f_tleaf  ', timelen, compress, &
         'leaf temperature','K')
   
      CALL hist_concatenate_var_2d (filehist, 'f_ldew   ', timelen, compress, &
         'depth of water on foliage','mm')
   
      CALL hist_concatenate_var_2d (filehist, 'f_scv    ', timelen, compress, &
         'snow cover, water equivalent','mm')
   
      CALL hist_concatenate_var_2d (filehist, 'f_snowdp ', timelen, compress, &
         'snow depth','meter')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fsno   ', timelen, compress, &
         'fraction of snow cover on ground','-')
   
      CALL hist_concatenate_var_2d (filehist, 'f_sigf   ', timelen, compress, &
         'fraction of veg cover, excluding snow-covered veg','-')
   
      CALL hist_concatenate_var_2d (filehist, 'f_green  ', timelen, compress, &
         'leaf greenness','-')
   
      CALL hist_concatenate_var_2d (filehist, 'f_lai    ', timelen, compress, &
         'leaf area index','m2/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_laisun ', timelen, compress, &
         'sunlit leaf area index','m2/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_laisha ', timelen, compress, &
         'shaded leaf area index','m2/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_sai    ', timelen, compress, &
         'stem area index','m2/m2')
   
      CALL hist_concatenate_var_4d (filehist, 'f_alb    ', timelen, &
         'band', 'rtyp', 2, 2, compress, 'averaged albedo direct','%')
   
      CALL hist_concatenate_var_2d (filehist, 'f_emis   ', timelen, compress, &
         'averaged bulk surface emissivity','-')
   
      CALL hist_concatenate_var_2d (filehist, 'f_z0m    ', timelen, compress, &
         'effective roughness','m')
   
      CALL hist_concatenate_var_2d (filehist, 'f_trad   ', timelen, compress, &
         'radiative temperature of surface','kelvin')
   
      CALL hist_concatenate_var_2d (filehist, 'f_tref   ', timelen, compress, &
         '2 m height air temperature','kelvin')
   
      CALL hist_concatenate_var_2d (filehist, 'f_qref   ', timelen, compress, &
         '2 m height air specific humidity','kg/kg')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xy_rain', timelen, compress, &
         'rain','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xy_snow', timelen, compress, &
         'snow','mm/s')
   
#ifdef BGC
            ! leaf carbon display pool
      CALL hist_concatenate_var_2d (filehist, 'f_leafc', timelen, compress, &
                'leaf carbon display pool','gC/m2')
   
            ! leaf carbon storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_leafc_storage', timelen, compress, &
                'leaf carbon storage pool','gC/m2')
   
            ! leaf carbon transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_leafc_xfer', timelen, compress, &
                'leaf carbon transfer pool','gC/m2')
   
            ! fine root carbon display pool
      CALL hist_concatenate_var_2d (filehist, 'f_frootc', timelen, compress, &
                'fine root carbon display pool','gC/m2')
   
            ! fine root carbon storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_frootc_storage', timelen, compress, &
                'fine root carbon storage pool','gC/m2')
   
            ! fine root carbon transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_frootc_xfer', timelen, compress, &
                'fine root carbon transfer pool','gC/m2')
   
            ! live stem carbon display pool
      CALL hist_concatenate_var_2d (filehist, 'f_livestemc', timelen, compress, &
                'live stem carbon display pool','gC/m2')
   
            ! live stem carbon storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_livestemc_storage', timelen, compress, &
                'live stem carbon storage pool','gC/m2')
   
            ! live stem carbon transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_livestemc_xfer', timelen, compress, &
                'live stem carbon transfer pool','gC/m2')
   
            ! dead stem carbon display pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadstemc', timelen, compress, &
                'dead stem carbon display pool','gC/m2')
   
            ! dead stem carbon storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadstemc_storage', timelen, compress, &
                'dead stem carbon storage pool','gC/m2')
   
            ! dead stem carbon transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadstemc_xfer', timelen, compress, &
                'dead stem carbon transfer pool','gC/m2')
   
            ! live coarse root carbon display pool
      CALL hist_concatenate_var_2d (filehist, 'f_livecrootc', timelen, compress, &
                'live coarse root carbon display pool','gC/m2')
   
            ! live coarse root carbon storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_livecrootc_storage', timelen, compress, &
                'live coarse root carbon storage pool','gC/m2')
   
            ! live coarse root carbon transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_livecrootc_xfer', timelen, compress, &
                'live coarse root carbon transfer pool','gC/m2')
   
            ! dead coarse root carbon display pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadcrootc', timelen, compress, &
                'dead coarse root carbon display pool','gC/m2')
   
            ! dead coarse root carbon storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadcrootc_storage', timelen, compress, &
                'dead coarse root carbon storage pool','gC/m2')
   
            ! dead coarse root carbon transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadcrootc_xfer', timelen, compress, &
                'dead coarse root carbon transfer pool','gC/m2')
   
            ! grain carbon display pool
      CALL hist_concatenate_var_2d (filehist, 'f_grainc', timelen, compress, &
                'grain carbon display pool','gC/m2')
   
            ! grain carbon storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_grainc_storage', timelen, compress, &
                'grain carbon storage pool','gC/m2')
   
            ! grain carbon transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_grainc_xfer', timelen, compress, &
                'grain carbon transfer pool','gC/m2')
   
            ! leaf nitrogen display pool
      CALL hist_concatenate_var_2d (filehist, 'f_leafn', timelen, compress, &
                'leaf nitrogen display pool','gN/m2')
   
            ! leaf nitrogen storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_leafn_storage', timelen, compress, &
                'leaf nitrogen storage pool','gN/m2')
   
            ! leaf nitrogen transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_leafn_xfer', timelen, compress, &
                'leaf nitrogen transfer pool','gN/m2')
   
            ! fine root nitrogen display pool
      CALL hist_concatenate_var_2d (filehist, 'f_frootn', timelen, compress, &
                'fine root nitrogen display pool','gN/m2')
   
            ! fine root nitrogen storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_frootn_storage', timelen, compress, &
                'fine root nitrogen storage pool','gN/m2')
   
            ! fine root nitrogen transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_frootn_xfer', timelen, compress, &
                'fine root nitrogen transfer pool','gN/m2')
   
            ! live stem nitrogen display pool
      CALL hist_concatenate_var_2d (filehist, 'f_livestemn', timelen, compress, &
                'live stem nitrogen display pool','gN/m2')
   
            ! live stem nitrogen storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_livestemn_storage', timelen, compress, &
                'live stem nitrogen storage pool','gN/m2')
   
            ! live stem nitrogen transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_livestemn_xfer', timelen, compress, &
                'live stem nitrogen transfer pool','gN/m2')
   
            ! dead stem nitrogen display pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadstemn', timelen, compress, &
                'dead stem nitrogen display pool','gN/m2')
   
            ! dead stem nitrogen storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadstemn_storage', timelen, compress, &
                'dead stem nitrogen storage pool','gN/m2')
   
            ! dead stem nitrogen transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadstemn_xfer', timelen, compress, &
                'dead stem nitrogen transfer pool','gN/m2')
   
            ! live coarse root nitrogen display pool
      CALL hist_concatenate_var_2d (filehist, 'f_livecrootn', timelen, compress, &
                'live coarse root nitrogen display pool','gN/m2')
   
            ! live coarse root nitrogen storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_livecrootn_storage', timelen, compress, &
                'live coarse root nitrogen storage pool','gN/m2')
   
            ! live coarse root nitrogen transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_livecrootn_xfer', timelen, compress, &
                'live coarse root nitrogen transfer pool','gN/m2')
   
            ! dead coarse root nitrogen display pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadcrootn', timelen, compress, &
                'dead coarse root nitrogen display pool','gN/m2')
   
            ! dead coarse root nitrogen storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadcrootn_storage', timelen, compress, &
                'dead coarse root nitrogen storage pool','gN/m2')
   
            ! dead coarse root nitrogen transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_deadcrootn_xfer', timelen, compress, &
                'dead coarse root nitrogen transfer pool','gN/m2')
   
            ! grain nitrogen display pool
      CALL hist_concatenate_var_2d (filehist, 'f_grainn', timelen, compress, &
                'grain nitrogen display pool','gN/m2')
   
            ! grain nitrogen storage pool
      CALL hist_concatenate_var_2d (filehist, 'f_grainn_storage', timelen, compress, &
                'grain nitrogen storage pool','gN/m2')
   
            ! grain nitrogen transfer pool
      CALL hist_concatenate_var_2d (filehist, 'f_grainn_xfer', timelen, compress, &
                'grain nitrogen transfer pool','gN/m2')
   
            ! retranslocation nitrogen pool
      CALL hist_concatenate_var_2d (filehist, 'f_retrasn', timelen, compress, &
                'retranslocation nitrogen pool','gN/m2')
   
            ! gross primary productivity
      CALL hist_concatenate_var_2d (filehist, 'f_gpp', timelen, compress, &
                'gross primary productivity','gC/m2/s')
   
            ! gross primary productivity
      CALL hist_concatenate_var_2d (filehist, 'f_downreg', timelen, compress, &
                'gpp downregulation due to N limitation','unitless')
   
            ! autotrophic respiration
      CALL hist_concatenate_var_2d (filehist, 'f_ar ', timelen, compress, &
                'autotrophic respiration','gC/m2/s')
   
#ifdef CROP
            ! crop phase
      CALL hist_concatenate_var_2d (filehist, 'f_cphase', timelen, compress, &
                'crop phase','unitless')
   
            ! 1-yr crop production carbon
      CALL hist_concatenate_var_2d (filehist, 'f_cropprod1c', timelen, compress, &
                '1-yr crop production carbon','gC/m2')
   
            ! loss rate of 1-yr crop production carbon
      CALL hist_concatenate_var_2d (filehist, 'f_cropprod1c_loss', timelen, compress, &
                'loss rate of 1-yr crop production carbon','gC/m2/s')
   
            ! crop seed deficit
      CALL hist_concatenate_var_2d (filehist, 'f_cropseedc_deficit', timelen, compress, &
                'crop seed deficit','gC/m2/s')
   
            ! grain to crop production carbon
      CALL hist_concatenate_var_2d (filehist, 'f_grainc_to_cropprodc', timelen, compress, &
                'grain to crop production carbon','gC/m2/s')
   
            ! grain to crop seed carbon
      CALL hist_concatenate_var_2d (filehist, 'f_grainc_to_seed', timelen, compress, &
                'grain to crop seed carbon','gC/m2/s')
#endif
#endif
      !--------------------------------------------------------------
      CALL hist_concatenate_var_3d (filehist, 'f_t_soisno   ', timelen, 'soilsnow', &
         nl_soil-maxsnl, compress, 'soil temperature','K')
   
      CALL hist_concatenate_var_3d (filehist, 'f_wliq_soisno', timelen, 'soilsnow', &
         nl_soil-maxsnl, compress, 'liquid water in soil layers','kg/m2')
   
      CALL hist_concatenate_var_3d (filehist, 'f_wice_soisno', timelen, 'soilsnow', &
         nl_soil-maxsnl, compress, 'ice lens in soil layers','kg/m2')
   
      CALL hist_concatenate_var_3d (filehist, 'f_h2osoi     ', timelen, 'soil', &
         nl_soil, compress, 'volumetric water in soil layers','m3/m3')
   
      ! CALL hist_concatenate_var_2d (filehist, 'f_rstfacsun  ', timelen, compress)  ! factor of soil water stress to transpiration on sunlit leaf
      ! CALL hist_concatenate_var_2d (filehist, 'f_rstfacsha  ', timelen, compress)  ! factor of soil water stress to transpiration on shaded leaf
   
      ! fraction of root water uptake from each soil layer (PHS undefined)
      ! water exchange between soil layers and root. Positive: soil->root [mm h2o/s] (PHS defined)
      CALL hist_concatenate_var_3d (filehist, 'f_rootr      ', timelen, 'soil', nl_soil, compress)
   
      IF(DEF_USE_PLANTHYDRAULICS)THEN
         CALL hist_concatenate_var_3d (filehist, 'f_vegwp      ', timelen, 'vegnodes', nvegwcs, compress)  ! vegetation water potential [mm]
      ENDIF
   
      CALL hist_concatenate_var_2d (filehist, 'f_zwt        ', timelen, compress, &
         'the depth to water table','m')
   
      CALL hist_concatenate_var_2d (filehist, 'f_wa         ', timelen, compress, &
         'water storage in aquifer','mm')
   
      CALL hist_concatenate_var_2d (filehist, 'f_wat        ', timelen, compress, &
         'total water storage','mm')
   
      CALL hist_concatenate_var_3d (filehist, 'f_t_lake      ', timelen, 'lake', &
         nl_lake, compress, 'lake temperature','K')
   
      CALL hist_concatenate_var_3d (filehist, 'f_lake_icefrac', timelen, 'lake', &
         nl_lake, compress, 'lake ice fraction cover','0-1')
   
#ifdef BGC
            ! litter 1 carbon density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_litr1c_vr', timelen, 'soil', &
         nl_soil, compress, 'litter 1 carbon density in soil layers','gC/m3')
   
            ! litter 2 carbon density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_litr2c_vr', timelen, 'soil', &
         nl_soil, compress, 'litter 2 carbon density in soil layers','gC/m3')
   
            ! litter 3 carbon density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_litr3c_vr', timelen, 'soil', &
         nl_soil, compress, 'litter 3 carbon density in soil layers','gC/m3')
   
            ! soil 1 carbon density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_soil1c_vr', timelen, 'soil', &
         nl_soil, compress, 'soil 1 carbon density in soil layers','gC/m3')
   
            ! soil 2 carbon density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_soil2c_vr', timelen, 'soil', &
         nl_soil, compress, 'soil 2 carbon density in soil layers','gC/m3')
   
            ! soil 3 carbon density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_soil3c_vr', timelen, 'soil', &
         nl_soil, compress, 'soil 3 carbon density in soil layers','gC/m3')
   
            ! coarse woody debris carbon density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_cwdc_vr', timelen, 'soil', &
         nl_soil, compress, 'coarse woody debris carbon density in soil layers','gC/m3')
   
            ! litter 1 nitrogen density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_litr1n_vr', timelen, 'soil', &
         nl_soil, compress, 'litter 1 nitrogen density in soil layers','gN/m3')
   
            ! litter 2 nitrogen density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_litr2n_vr', timelen, 'soil', &
         nl_soil, compress, 'litter 2 nitrogen density in soil layers','gN/m3')
   
            ! litter 3 nitrogen density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_litr3n_vr', timelen, 'soil', &
         nl_soil, compress, 'litter 3 nitrogen density in soil layers','gN/m3')
   
            ! soil 1 nitrogen density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_soil1n_vr', timelen, 'soil', &
         nl_soil, compress, 'soil 1 nitrogen density in soil layers','gN/m3')
   
            ! soil 2 nitrogen density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_soil2n_vr', timelen, 'soil', &
         nl_soil, compress, 'soil 2 nitrogen density in soil layers','gN/m3')
   
            ! soil 3 nitrogen density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_soil3n_vr', timelen, 'soil', &
         nl_soil, compress, 'soil 3 nitrogen density in soil layers','gN/m3')
   
            ! coarse woody debris nitrogen density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_cwdn_vr', timelen, 'soil', &
         nl_soil, compress, 'coarse woody debris nitrogen density in soil layers','gN/m3')
   
            ! mineral nitrogen density in soil layers
      CALL hist_concatenate_var_3d (filehist, 'f_sminn_vr', timelen, 'soil', &
         nl_soil, compress, 'mineral nitrogen density in soil layers','gN/m3')
   
#endif
      CALL hist_concatenate_var_2d (filehist, 'f_ustar  ', timelen, compress, &
         'u* in similarity theory','m/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_tstar  ', timelen, compress, &
         't* in similarity theory','K')
   
      CALL hist_concatenate_var_2d (filehist, 'f_qstar  ', timelen, compress, &
         'q* in similarity theory')
   
      CALL hist_concatenate_var_2d (filehist, 'f_zol    ', timelen, compress, &
         'dimensionless height (z/L) used in Monin-Obukhov theory','-')
   
      CALL hist_concatenate_var_2d (filehist, 'f_rib    ', timelen, compress, &
         'bulk Richardson number in surface layer','-')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fm     ', timelen, compress, &
         'integral of profile function for momentum','-')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fh     ', timelen, compress, &
         'integral of profile function for heat','-')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fq     ', timelen, compress, &
         'integral of profile function for moisture','-')
   
      CALL hist_concatenate_var_2d (filehist, 'f_us10m  ', timelen, compress, &
         '10m u-velocity','m/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_vs10m  ', timelen, compress, &
         '10m v-velocity','m/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_fm10m  ', timelen, compress, &
         'integral of profile function for momentum at 10m','-')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xy_us  ', timelen, compress, &
         'wind in eastward direction', 'm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xy_vs  ', timelen, compress, &
         'wind in northward direction','m/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xy_t   ', timelen, compress, &
         'temperature at reference height','kelvin')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xy_q   ', timelen, compress, &
         'specific humidity at reference height','kg/kg')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xy_prc ', timelen, compress, &
         'convective precipitation','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xy_prl ', timelen, compress, &
         'large scale precipitation','mm/s')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xy_pbot', timelen, compress, &
         'atmospheric pressure at the surface','pa')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xy_frl ', timelen, compress, &
         'atmospheric infrared (longwave) radiation','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_xy_solarin', timelen, compress, &
         'downward solar radiation at surface','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_sr     ', timelen, compress, &
         'reflected solar radiation at surface [W/m2]','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_solvd  ', timelen, compress, &
         'incident direct beam vis solar radiation (W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_solvi  ', timelen, compress, &
         'incident diffuse beam vis solar radiation (W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_solnd  ', timelen, compress, &
         'incident direct beam nir solar radiation (W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_solni  ', timelen, compress, &
         'incident diffuse beam nir solar radiation (W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_srvd   ', timelen, compress, &
         'reflected direct beam vis solar radiation (W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_srvi   ', timelen, compress, &
         'reflected diffuse beam vis solar radiation (W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_srnd   ', timelen, compress, &
         'reflected direct beam nir solar radiation (W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_srni   ', timelen, compress, &
         'reflected diffuse beam nir solar radiation (W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_solvdln', timelen, compress, &
         'incident direct beam vis solar radiation at local noon(W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_solviln', timelen, compress, &
         'incident diffuse beam vis solar radiation at local noon(W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_solndln', timelen, compress, &
         'incident direct beam nir solar radiation at local noon(W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_solniln', timelen, compress, &
         'incident diffuse beam nir solar radiation at local noon(W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_srvdln ', timelen, compress, &
         'reflected direct beam vis solar radiation at local noon(W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_srviln ', timelen, compress, &
         'reflected diffuse beam vis solar radiation at local noon(W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_srndln ', timelen, compress, &
         'reflected direct beam nir solar radiation at local noon(W/m2)','W/m2')
   
      CALL hist_concatenate_var_2d (filehist, 'f_srniln ', timelen, compress, &
         'reflected diffuse beam nir solar radiation at local noon(W/m2)','W/m2')

END PROGRAM hist_concatenate
! ----------------------------------------------------------------------
! EOP
