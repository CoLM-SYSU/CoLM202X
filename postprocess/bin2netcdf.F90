
#define lon_points 720
#define lat_points 360

PROGRAM bin2netcdf
! ================================================================
!  PURPOSE:
!     post process for CoLM output file. convert binary data to NetCDF
!  format. 
!
!  USAGE: 
!     1) for CoLM output format
!        ./bin2netcdf [casename]_2D_Fluxes_YYYY-MM
!
!     2) for NCAR CLM output format (can be used as input for Land 
!        Model Diagnostics Package)
!        ./bin2netcdf [casename]_2D_Fluxes_YYYY-MM ncar
! ================================================================

!TODO: 根据main修改的代码修改输出变量
   use precision
   USE GlobalVars
   use PhysicalConstants, only: hvap
   use MOD_2D_Fluxes
   use netcdf
   
   implicit none

   !integer,  parameter :: nl_soil    = 10      ! number of soil layers
   !integer,  parameter :: nl_lake    = 10      ! number of lake layers
   !integer,  parameter :: maxsnl     = -5      ! max number of snow layers

   real(r4), parameter :: spval_r4   = -1.e36_r4 ! a special value for missing and filling use

   real(r8) :: lons_r8(lon_points)
   real(r8) :: lats_r8(lat_points)
   
   character(len=256)  :: filename
   character(len=256)  :: outtype

 ! read in the file name
   CALL getarg(1, filename)
  
 ! allocate memory
   print *, 'allocate memory...'
   CALL allocate_2D_Fluxes(lon_points,lat_points)

 ! read in flux data
   print *, 'read in flux data...'
   CALL readfluxes

 ! write to NetCDF format
   print *, 'write to NetCDF format...'
   CALL getarg(2, outtype)
   if (trim(outtype) == 'ncar') then
      CALL writenetcdf_ncar
   else
      CALL writenetcdf
   end if

 ! free memory
   print *, 'free memory. done.'
   CALL deallocate_2D_Fluxes     
   

CONTAINS

   SUBROUTINE readfluxes
      
      implicit none
      
      open(unit=11, file=trim(filename), access='sequential', form='unformatted', status='old')
      print *, 'open file: ', trim(filename)

      read(11) lons_r8  (:)     ! longitudes in degree
      read(11) lats_r8  (:)     ! latitudes in degree
      read(11) mask     (:,:)   ! grid mask [1: target (land), 0: none]
      read(11) frac     (:,:)   ! grid total fraction [fraction]
      read(11) area     (:,:)   ! grid cell area [km2]
      read(11) f_taux   (:,:)   ! wind stress: E-W [kg/m/s2]
      read(11) f_tauy   (:,:)   ! wind stress: N-S [kg/m/s2]
      read(11) f_fsena  (:,:)   ! sensible heat from canopy height to atmosphere [W/m2]
      read(11) f_lfevpa (:,:)   ! latent heat flux from canopy height to atmosphere [W/m2]
      read(11) f_fevpa  (:,:)   ! evapotranspiration from canopy to atmosphere [mm/s]
      read(11) f_fsenl  (:,:)   ! sensible heat from leaves [W/m2]
      read(11) f_fevpl  (:,:)   ! evaporation+transpiration from leaves [mm/s]
      read(11) f_etr    (:,:)   ! transpiration rate [mm/s]
      read(11) f_fseng  (:,:)   ! sensible heat flux from ground [W/m2]
      read(11) f_fevpg  (:,:)   ! evaporation heat flux from ground [mm/s]
      read(11) f_fgrnd  (:,:)   ! ground heat flux [W/m2]
      read(11) f_sabvsun(:,:)   ! solar absorbed by sunlit canopy [W/m2]
      read(11) f_sabvsha(:,:)   ! solar absorbed by shaded [W/m2]
      read(11) f_sabg   (:,:)   ! solar absorbed by ground  [W/m2]
      read(11) f_olrg   (:,:)   ! outgoing long-wave radiation from ground+canopy [W/m2]
      read(11) f_rnet   (:,:)   ! net radiation [W/m2]
      read(11) f_xerr   (:,:)   ! the error of water banace [mm/s]
      read(11) f_zerr   (:,:)   ! the error of energy balance [W/m2]
      read(11) f_rsur   (:,:)   ! surface runoff [mm/s]
      read(11) f_rnof   (:,:)   ! total runoff [mm/s]
      read(11) f_qintr  (:,:)   ! interception [mm/s]
      read(11) f_qinfl  (:,:)   ! inflitration [mm/s]
      read(11) f_qdrip  (:,:)   ! throughfall [mm/s]
      read(11) f_assim  (:,:)   ! canopy assimilation rate [mol m-2 s-1]
      read(11) f_respc  (:,:)   ! respiration (plant+soil) [mol m-2 s-1]
      read(11) f_qcharge(:,:)   ! groundwater recharge rate [mm/s] 

!---------------------------------------------------------------------
      read(11) f_t_grnd (:,:)   ! ground surface temperature [K]
      read(11) f_tleaf  (:,:)   ! sunlit leaf temperature [K]
      !read(11) f_tlsha  (:,:)   ! shaded leaf temperature [K]
      read(11) f_ldew   (:,:)   ! depth of water on foliage [mm]
      read(11) f_scv    (:,:)   ! snow cover, water equivalent [mm]
      read(11) f_snowdp (:,:)   ! snow depth [meter]
      read(11) f_fsno   (:,:)   ! fraction of snow cover on ground
      read(11) f_sigf   (:,:)   ! fraction of veg cover, excluding snow-covered veg [-]
      read(11) f_green  (:,:)   ! leaf greenness
      read(11) f_lai    (:,:)   ! leaf area index
      read(11) f_laisun (:,:)   ! sunlit leaf area index
      read(11) f_laisha (:,:)   ! shaded leaf area index
      read(11) f_sai    (:,:)   ! stem area index
      read(11) f_alb(:,:,:,:)   ! averaged albedo [visible, direct; direct, diffuse]
      read(11) f_emis   (:,:)   ! averaged bulk surface emissivity
      read(11) f_z0m    (:,:)   ! effective roughness [m]
      read(11) f_trad   (:,:)   ! radiative temperature of surface [K]
      read(11) f_tref   (:,:)   ! 2 m height air temperature [kelvin]
      read(11) f_qref   (:,:)   ! 2 m height air specific humidity [kg/kg]
      read(11) f_xy_rain(:,:)   ! rain [mm/s]
      read(11) f_xy_snow(:,:)   ! snow[mm/s]

!---------------------------------------------------------------------
      read(11) f_t_soisno   (:,:,:)   ! soil temperature [K]
      read(11) f_wliq_soisno(:,:,:)   ! liquid water in soil layers [kg/m2]
      read(11) f_wice_soisno(:,:,:)   ! ice lens in soil layers [kg/m2]
      read(11) f_h2osoi     (:,:,:)  ! volumetric soil water in layers [m3/m3]
      read(11) f_rstfac     (:,:)    ! factor of soil water stress 
      read(11) f_zwt        (:,:)    ! the depth to water table [m]
      read(11) f_wa         (:,:)    ! water storage in aquifer [mm]
      read(11) f_wat        (:,:)    ! total water storage [mm]

      read(11) f_t_lake     (:,:,:)   ! lake temperature [K]
      read(11) f_lake_icefrac (:,:,:)   ! lake ice fraction cover [0-1]

      read(11) f_ustar  (:,:)   ! u* in similarity theory [m/s]
      read(11) f_tstar  (:,:)   ! t* in similarity theory [kg/kg]
      read(11) f_qstar  (:,:)   ! q* in similarity theory [kg/kg]
      read(11) f_zol    (:,:)   ! dimensionless height (z/L) used in Monin-Obukhov theory
      read(11) f_rib    (:,:)   ! bulk Richardson number in surface layer
      read(11) f_fm     (:,:)   ! integral of profile function for momentum
      read(11) f_fh     (:,:)   ! integral of profile function for heat
      read(11) f_fq     (:,:)   ! integral of profile function for moisture
      read(11) f_us10m  (:,:)   ! 10m u-velocity [m/s]
      read(11) f_vs10m  (:,:)   ! 10m v-velocity [m/s]
      read(11) f_fm10m  (:,:)   ! integral of profile function for momentum at 10m [-]

      read(11) f_xy_us  (:,:)   ! wind in eastward direction [m/s]
      read(11) f_xy_vs  (:,:)   ! wind in northward direction [m/s]
      read(11) f_xy_t   (:,:)   ! temperature at reference height [kelvin]
      read(11) f_xy_q   (:,:)   ! specific humidity at reference height [kg/kg]
      read(11) f_xy_prc (:,:)   ! convective precipitation [mm/s]
      read(11) f_xy_prl (:,:)   ! large scale precipitation [mm/s]
      read(11) f_xy_pbot(:,:)   ! atmospheric pressure at the surface [pa]
      read(11) f_xy_frl (:,:)   ! atmospheric infrared (longwave) radiation [W/m2]
      read(11) f_xy_solarin(:,:)! downward solar radiation at surface [W/m2]

      read(11) f_sr     (:,:)   ! reflected solar radiation at surface [W/m2]
      read(11) f_solvd  (:,:)   ! incident direct beam vis solar radiation (W/m2)
      read(11) f_solvi  (:,:)   ! incident diffuse beam vis solar radiation (W/m2)
      read(11) f_solnd  (:,:)   ! incident direct beam nir solar radiation (W/m2)
      read(11) f_solni  (:,:)   ! incident diffuse beam nir solar radiation (W/m2)
      read(11) f_srvd   (:,:)   ! reflected direct beam vis solar radiation (W/m2)
      read(11) f_srvi   (:,:)   ! reflected diffuse beam vis solar radiation (W/m2)
      read(11) f_srnd   (:,:)   ! reflected direct beam nir solar radiation (W/m2)
      read(11) f_srni   (:,:)   ! reflected diffuse beam nir solar radiation (W/m2)
      read(11) f_solvdln(:,:)   ! incident direct beam vis solar radiation at local noon(W/m2)
      read(11) f_solviln(:,:)   ! incident diffuse beam vis solar radiation at local noon(W/m2)
      read(11) f_solndln(:,:)   ! incident direct beam nir solar radiation at local noon(W/m2)
      read(11) f_solniln(:,:)   ! incident diffuse beam nir solar radiation at local noon(W/m2)
      read(11) f_srvdln (:,:)   ! reflected direct beam vis solar radiation at local noon(W/m2)
      read(11) f_srviln (:,:)   ! reflected diffuse beam vis solar radiation at local noon(W/m2)
      read(11) f_srndln (:,:)   ! reflected direct beam nir solar radiation at local noon(W/m2)
      read(11) f_srniln (:,:)   ! reflected diffuse beam nir solar radiation at local noon(W/m2)

      close(11)

   END SUBROUTINE readfluxes
   
   SUBROUTINE writenetcdf
   
      implicit none

      integer  :: ix,iy,ilev,i
      integer  :: ncid, xid, yid, varid, sslevid, lakelevid, bandid
      integer  :: lake_lev(nl_lake), sslev(nl_soil-maxsnl), band(2)
      real(r4) :: lons(lon_points), lats(lat_points)
      real(r8) :: tmp(lon_points, lat_points, 2)
      real(r8) :: tmp1(lon_points, lat_points, maxsnl+1:nl_soil)
      real(r8) :: tmp2(lon_points, lat_points, nl_lake)
      
    ! create netcdf and define dimensions
      call sanity( nf90_create(trim(filename)//'.nc', nf90_64bit_offset, ncid) )

    ! dimensions
      call sanity( nf90_def_dim(ncid, 'lon',           lon_points,     xid) )
      call sanity( nf90_def_dim(ncid, 'lat',           lat_points,     yid) )
      call sanity( nf90_def_dim(ncid, 'soil_snow_lev', nl_soil-maxsnl, sslevid) )
      call sanity( nf90_def_dim(ncid, 'lake_lev',      nl_lake,        lakelevid) )
      call sanity( nf90_def_dim(ncid, 'band',          2,              bandid) )

    ! global attr
      call sanity( nf90_put_att(ncid, NF90_GLOBAL, 'title','CLM SIMULATED SURFACE FLUXES') )
      
    ! variables
    ! dimension variables
      call sanity( nf90_def_var(ncid, 'lon', nf90_float, (/xid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','longitude') )
      call sanity( nf90_put_att(ncid, varid, 'units','degrees_east') )

      call sanity( nf90_def_var(ncid, 'lat', nf90_float, (/yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','latitude') )
      call sanity( nf90_put_att(ncid, varid, 'units','degrees_north') )

      call sanity( nf90_def_var(ncid, 'soil_snow_lev', nf90_int, (/sslevid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name',"soil plus snow level") )
      
      call sanity( nf90_def_var(ncid, 'lake_lev', nf90_int, (/lakelevid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name',"lake level") )
      
      call sanity( nf90_def_var(ncid, 'band', nf90_int, (/bandid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name',"band (vis/nir)") )
 
      ! grid mask
      call sanity( nf90_def_var(ncid, 'landmask', nf90_int, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','grid mask') )
      call sanity( nf90_put_att(ncid, varid, 'units','none') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', -1) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', -1) )

      ! grid total fraction
      call sanity( nf90_def_var(ncid, 'landfrac', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','grid total fraction') )
      call sanity( nf90_put_att(ncid, varid, 'units','fraction') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
       
      ! grid cell area [km2]
      call sanity( nf90_def_var(ncid, 'area', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','grid cell area [km2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','fraction') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
     
      ! wind stress: E-W [kg/m/s2]
      call sanity( nf90_def_var(ncid, 'f_taux', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','wind stress: E-W [kg/m/s2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/m/s2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! wind stress: N-S [kg/m/s2]
      call sanity( nf90_def_var(ncid, 'f_tauy', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','wind stress: N-S [kg/m/s2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/m/s2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! sensible heat from canopy height to atmosphere [W/m2]
      call sanity( nf90_def_var(ncid, 'f_fsena', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','sensible heat from canopy height to atmosphere [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! latent heat flux from canopy height to atmosphere [W/m2]
      call sanity( nf90_def_var(ncid, 'f_lfevpa', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','latent heat flux from canopy height to atmosphere [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! evapotranspiration from canopy to atmosphere [mm/s]
      call sanity( nf90_def_var(ncid, 'f_fevpa', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','evapotranspiration from canopy to atmosphere [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! sensible heat from leaves [W/m2]
      call sanity( nf90_def_var(ncid, 'f_fsenl', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','sensible heat from leaves [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! evaporation+transpiration from leaves [mm/s]
      call sanity( nf90_def_var(ncid, 'f_fevpl', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','evaporation+transpiration from leaves [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! transpiration rate [mm/s]
      call sanity( nf90_def_var(ncid, 'f_etr', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','transpiration rate [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! sensible heat flux from ground [W/m2]
      call sanity( nf90_def_var(ncid, 'f_fseng', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','sensible heat flux from ground [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! evaporation heat flux from ground [mm/s]
      call sanity( nf90_def_var(ncid, 'f_fevpg', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','evaporation heat flux from ground [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! ground heat flux [W/m2]
      call sanity( nf90_def_var(ncid, 'f_fgrnd', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','ground heat flux [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! solar absorbed by sunlit canopy [W/m2]
      call sanity( nf90_def_var(ncid, 'f_sabvsun', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','solar absorbed by sunlit canopy [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! solar absorbed by shaded [W/m2]
      call sanity( nf90_def_var(ncid, 'f_sabvsha', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','solar absorbed by shaded [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! solar absorbed by ground  [W/m2]
      call sanity( nf90_def_var(ncid, 'f_sabg', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','solar absorbed by ground  [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! outgoing long-wave radiation from ground+canopy [W/m2]
      call sanity( nf90_def_var(ncid, 'f_olrg', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','outgoing long-wave radiation from ground+canopy [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! net radiation [W/m2]
      call sanity( nf90_def_var(ncid, 'f_rnet', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','net radiation [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! the error of water banace [mm/s]
      call sanity( nf90_def_var(ncid, 'f_xerr', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','the error of water banace [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! the error of energy balance [W/m2]
      call sanity( nf90_def_var(ncid, 'f_zerr', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','the error of energy balance [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! surface runoff [mm/s]
      call sanity( nf90_def_var(ncid, 'f_rsur', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','surface runoff [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! total runoff [mm/s]
      call sanity( nf90_def_var(ncid, 'f_rnof', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','total runoff [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! canopy assimilation rate [mol m-2 s-1]
      call sanity( nf90_def_var(ncid, 'f_assim', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','canopy assimilation rate [mol m-2 s-1]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mol m-2 s-1') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! respiration (plant+soil) [mol m-2 s-1]
      call sanity( nf90_def_var(ncid, 'f_respc', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','respiration (plant+soil) [mol m-2 s-1]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mol m-2 s-1') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! groundwater recharge rate [mm/s] 
      call sanity( nf90_def_var(ncid, 'f_qcharge', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','groundwater recharge rate [mm/s] ') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      !---------------------------------------------------------------------
      ! ground surface temperature [K]
      call sanity( nf90_def_var(ncid, 'f_t_grnd', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','ground surface temperature [K]') )
      call sanity( nf90_put_att(ncid, varid, 'units','K') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! sunlit leaf temperature [K]
      call sanity( nf90_def_var(ncid, 'f_tleaf', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','sunlit leaf temperature [K]') )
      call sanity( nf90_put_att(ncid, varid, 'units','K') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! shaded leaf temperature [K]
      !call sanity( nf90_def_var(ncid, 'f_tleafsha', nf90_double, (/xid,yid/), varid) )
      !call sanity( nf90_put_att(ncid, varid, 'long_name','shaded leaf temperature [K]') )
      !call sanity( nf90_put_att(ncid, varid, 'units','K') )
      !call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      !call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! depth of water on foliage [mm]
      call sanity( nf90_def_var(ncid, 'f_ldew', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','depth of water on foliage [mm]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! snow cover, water equivalent [mm]
      call sanity( nf90_def_var(ncid, 'f_scv', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','snow cover, water equivalent [mm]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! snow depth [meter]
      call sanity( nf90_def_var(ncid, 'f_snowdp', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','snow depth [meter]') )
      call sanity( nf90_put_att(ncid, varid, 'units','meter') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! fraction of snow cover on ground [-]
      call sanity( nf90_def_var(ncid, 'f_fsno', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','fraction of snow cover on ground [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! fraction of veg cover, excluding snow-covered veg [-]
      call sanity( nf90_def_var(ncid, 'f_sigf', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','fraction of veg cover, excluding snow-covered veg [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! leaf greenness [fraction]
      call sanity( nf90_def_var(ncid, 'f_green', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','leaf greenness [fraction]') )
      call sanity( nf90_put_att(ncid, varid, 'units','fraction') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! leaf area index [m2/m2]
      call sanity( nf90_def_var(ncid, 'f_lai', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','leaf area index [m2/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m2/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! stem area index [m2/m2]
      call sanity( nf90_def_var(ncid, 'f_sai', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','stem area index [m2/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m2/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! averaged albedo direct [%]
      call sanity( nf90_def_var(ncid, 'f_albd', nf90_double, (/xid,yid,bandid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','averaged albedo direct [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! averaged albedo diffuse [%]
      call sanity( nf90_def_var(ncid, 'f_albi', nf90_double, (/xid,yid,bandid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','averaged albedo diffuse [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! averaged bulk surface emissivity [-]
      call sanity( nf90_def_var(ncid, 'f_emis', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','averaged bulk surface emissivity [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! effective roughness [m]
      call sanity( nf90_def_var(ncid, 'f_z0m', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','effective roughness [m]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! radiative temperature of surface [K]
      call sanity( nf90_def_var(ncid, 'f_trad', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','radiative temperature of surface [K]') )
      call sanity( nf90_put_att(ncid, varid, 'units','K') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! 2 m height air temperature [kelvin]
      call sanity( nf90_def_var(ncid, 'f_tref', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','2 m height air temperature [kelvin]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kelvin') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! 2 m height air specific humidity [kg/kg]
      call sanity( nf90_def_var(ncid, 'f_qref', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','2 m height air specific humidity [kg/kg]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/kg') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
       
      ! rain [mm/s]
      call sanity( nf90_def_var(ncid, 'f_xy_rain', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','rain [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
 
      ! snow [mm/s]
      call sanity( nf90_def_var(ncid, 'f_xy_snow', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','snow [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
 
      !---------------------------------------------------------------------
      ! soil temperature [K]
      call sanity( nf90_def_var(ncid, 'f_t_soisno', nf90_double, (/xid,yid,sslevid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','soil temperature [K]') )
      call sanity( nf90_put_att(ncid, varid, 'units','K') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! liquid water in soil layers [kg/m2]
      call sanity( nf90_def_var(ncid, 'f_wliq_soisno', nf90_double, (/xid,yid,sslevid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','liquid water in soil layers [kg/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! ice lens in soil layers [kg/m2]
      call sanity( nf90_def_var(ncid, 'f_wice_soisno', nf90_double, (/xid,yid,sslevid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','ice lens in soil layers [kg/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! lake temperature [K]
      call sanity( nf90_def_var(ncid, 'f_t_lake', nf90_double, (/xid,yid,lakelevid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','lake temperature [K]') )
      call sanity( nf90_put_att(ncid, varid, 'units','K') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! lake ice fraction cover [0-1]
      call sanity( nf90_def_var(ncid, 'f_lake_icefrac', nf90_double, (/xid,yid,lakelevid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','lake ice fraction cover [0-1]') )
      call sanity( nf90_put_att(ncid, varid, 'units','0-1') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! u* in similarity theory [m/s]
      call sanity( nf90_def_var(ncid, 'f_ustar', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','u* in similarity theory [m/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! t* in similarity theory [kg/kg]
      call sanity( nf90_def_var(ncid, 'f_tstar', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','t* in similarity theory [kg/kg]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/kg') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! q* in similarity theory [kg/kg]
      call sanity( nf90_def_var(ncid, 'f_qstar', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','q* in similarity theory [kg/kg]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/kg') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! dimensionless height (z/L) used in Monin-Obukhov theory [-]
      call sanity( nf90_def_var(ncid, 'f_zol', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','dimensionless height (z/L) used in Monin-Obukhov theory [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! bulk Richardson number in surface layer [-]
      call sanity( nf90_def_var(ncid, 'f_rib', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','bulk Richardson number in surface layer [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! integral of profile function for momentum [-]
      call sanity( nf90_def_var(ncid, 'f_fm', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','integral of profile function for momentum [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! integral of profile function for heat [-]
      call sanity( nf90_def_var(ncid, 'f_fh', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','integral of profile function for heat [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! integral of profile function for moisture [-]
      call sanity( nf90_def_var(ncid, 'f_fq', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','integral of profile function for moisture [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! 10m u-velocity [m/s]
      call sanity( nf90_def_var(ncid, 'f_us10m', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','10m u-velocity [m/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! 10m v-velocity [m/s]
      call sanity( nf90_def_var(ncid, 'f_vs10m', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','10m v-velocity [m/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! integral of profile function for momentum at 10m [-]
      call sanity( nf90_def_var(ncid, 'f_fm10m', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','integral of profile function for momentum at 10m [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! wind in eastward direction [m/s]
      call sanity( nf90_def_var(ncid, 'f_xy_us', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','wind in eastward direction [m/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! wind in northward direction [m/s]
      call sanity( nf90_def_var(ncid, 'f_xy_vs', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','wind in northward direction [m/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! temperature at reference height [kelvin]
      call sanity( nf90_def_var(ncid, 'f_xy_t', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','temperature at reference height [kelvin]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kelvin') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! specific humidity at reference height [kg/kg]
      call sanity( nf90_def_var(ncid, 'f_xy_q', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','specific humidity at reference height [kg/kg]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/kg') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! convective precipitation [mm/s]
      call sanity( nf90_def_var(ncid, 'f_xy_prc', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','convective precipitation [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! large scale precipitation [mm/s]
      call sanity( nf90_def_var(ncid, 'f_xy_prl', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','large scale precipitation [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! atmospheric pressure at the surface [pa]
      call sanity( nf90_def_var(ncid, 'f_xy_pbot', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','atmospheric pressure at the surface [pa]') )
      call sanity( nf90_put_att(ncid, varid, 'units','pa') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! atmospheric infrared (longwave) radiation [W/m2]
      call sanity( nf90_def_var(ncid, 'f_xy_frl', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','atmospheric infrared (longwave) radiation [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )

      ! downward solar radiation at surface [W/m2]
      call sanity( nf90_def_var(ncid, 'f_xy_solarin', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','downward solar radiation at surface [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! reflected solar radiation at surface [W/m2]
      call sanity( nf90_def_var(ncid, 'f_sr', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected solar radiation at surface [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! incident direct beam vis solar radiation (W/m2)
      call sanity( nf90_def_var(ncid, 'f_solvd', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident direct beam vis solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! incident diffuse beam vis solar radiation (W/m2)
      call sanity( nf90_def_var(ncid, 'f_solvi', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam vis solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! incident direct beam nir solar radiation (W/m2)
      call sanity( nf90_def_var(ncid, 'f_solnd', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident direct beam nir solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! incident diffuse beam nir solar radiation (W/m2)
      call sanity( nf90_def_var(ncid, 'f_solni', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam nir solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! reflected direct beam vis solar radiation (W/m2)
      call sanity( nf90_def_var(ncid, 'f_srvd', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected direct beam vis solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! reflected diffuse beam vis solar radiation (W/m2)
      call sanity( nf90_def_var(ncid, 'f_srvi', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam vis solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! reflected direct beam nir solar radiation (W/m2)
      call sanity( nf90_def_var(ncid, 'f_srnd', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected direct beam nir solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! reflected diffuse beam nir solar radiation (W/m2)
      call sanity( nf90_def_var(ncid, 'f_srni', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam nir solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! incident direct beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_def_var(ncid, 'f_solvdln', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident direct beam vis solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! incident diffuse beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_def_var(ncid, 'f_solviln', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam vis solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! incident direct beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_def_var(ncid, 'f_solndln', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident direct beam nir solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! incident diffuse beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_def_var(ncid, 'f_solniln', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam nir solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! reflected direct beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_def_var(ncid, 'f_srvdln', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected direct beam vis solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! reflected diffuse beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_def_var(ncid, 'f_srviln', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam vis solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! reflected direct beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_def_var(ncid, 'f_srndln', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected direct beam nir solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      
      ! reflected diffuse beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_def_var(ncid, 'f_srniln', nf90_double, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam nir solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval) )
      

    ! end defination
      call sanity( nf90_enddef(ncid) )

    ! write data
    ! ------------------------------------------------------------

    ! dimension data
      lons = lons_r8
      call sanity( nf90_inq_varid(ncid,'lon',varid) )
      call sanity( nf90_put_var(ncid,varid,lons) )

      lats = lats_r8
      call sanity( nf90_inq_varid(ncid,'lat',varid) )
      call sanity( nf90_put_var(ncid,varid,lats) )

      do ilev = 1, nl_soil-maxsnl
         sslev(ilev) = maxsnl+ilev
      end do
      call sanity( nf90_inq_varid(ncid,'soil_snow_lev',varid) )
      call sanity( nf90_put_var(ncid,varid,sslev) )

      do ilev = 1, nl_lake
         lake_lev(ilev) = ilev
      end do
      call sanity( nf90_inq_varid(ncid,'lake_lev',varid) )
      call sanity( nf90_put_var(ncid,varid,lake_lev) )

      do ilev = 1, 2
         band(ilev) = ilev
      end do
      call sanity( nf90_inq_varid(ncid,'band',varid) )
      call sanity( nf90_put_var(ncid,varid,band) )
      
      ! grid mask
      call sanity( nf90_inq_varid(ncid,'landmask',varid) )
      call sanity( nf90_put_var(ncid,varid,mask) )
      
      ! grid total fraction
      call sanity( nf90_inq_varid(ncid,'landfrac',varid) )
      call sanity( nf90_put_var(ncid,varid,frac) )
 
      ! grid cell area [km2]
      call sanity( nf90_inq_varid(ncid,'area',varid) )
      call sanity( nf90_put_var(ncid,varid,area) )
 
      ! wind stress: E-W [kg/m/s2]
      call sanity( nf90_inq_varid(ncid,'f_taux',varid) )
      call sanity( nf90_put_var(ncid,varid,f_taux) )
      
      ! wind stress: N-S [kg/m/s2]
      call sanity( nf90_inq_varid(ncid,'f_tauy',varid) )
      call sanity( nf90_put_var(ncid,varid,f_tauy) )
      
      ! sensible heat from canopy height to atmosphere [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_fsena',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fsena) )
      
      ! latent heat flux from canopy height to atmosphere [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_lfevpa',varid) )
      call sanity( nf90_put_var(ncid,varid,f_lfevpa) )
     
      ! evapotranspiration from canopy to atmosphere [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_fevpa',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fevpa) )
      
      ! sensible heat from leaves [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_fsenl',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fsenl) )
      
      ! evaporation+transpiration from leaves [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_fevpl',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fevpl) )
      
      ! transpiration rate [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_etr',varid) )
      call sanity( nf90_put_var(ncid,varid,f_etr) )
      
      ! sensible heat flux from ground [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_fseng',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fseng) )
      
      ! evaporation heat flux from ground [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_fevpg',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fevpg) )
      
      ! ground heat flux [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_fgrnd',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fgrnd) )
      
      ! solar absorbed by sunlit canopy [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_sabvsun',varid) )
      call sanity( nf90_put_var(ncid,varid,f_sabvsun) )
      
      ! solar absorbed by shaded [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_sabvsha',varid) )
      call sanity( nf90_put_var(ncid,varid,f_sabvsha) )
      
      ! solar absorbed by ground  [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_sabg',varid) )
      call sanity( nf90_put_var(ncid,varid,f_sabg) )
      
      ! outgoing long-wave radiation from ground+canopy [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_olrg',varid) )
      call sanity( nf90_put_var(ncid,varid,f_olrg) )
      
      ! net radiation [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_rnet',varid) )
      call sanity( nf90_put_var(ncid,varid,f_rnet) )
      
      ! the error of water banace [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_xerr',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xerr) )
      
      ! the error of energy balance [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_zerr',varid) )
      call sanity( nf90_put_var(ncid,varid,f_zerr) )
      
      ! surface runoff [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_rsur',varid) )
      call sanity( nf90_put_var(ncid,varid,f_rsur) )
      
      ! total runoff [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_rnof',varid) )
      call sanity( nf90_put_var(ncid,varid,f_rnof) )
      
      ! canopy assimilation rate [mol m-2 s-1]
      call sanity( nf90_inq_varid(ncid,'f_assim',varid) )
      call sanity( nf90_put_var(ncid,varid,f_assim) )
      
      ! respiration (plant+soil) [mol m-2 s-1]
      call sanity( nf90_inq_varid(ncid,'f_respc',varid) )
      call sanity( nf90_put_var(ncid,varid,f_respc) )
      
      ! groundwater recharge rate [mm/s] 
      call sanity( nf90_inq_varid(ncid,'f_qcharge',varid) )
      call sanity( nf90_put_var(ncid,varid,f_qcharge) )
      
!---------------------------------------------------------------------
      ! ground surface temperature [K]
      call sanity( nf90_inq_varid(ncid,'f_t_grnd',varid) )
      call sanity( nf90_put_var(ncid,varid,f_t_grnd) )
      
      ! sunlit leaf temperature [K]
      call sanity( nf90_inq_varid(ncid,'f_tleaf',varid) )
      call sanity( nf90_put_var(ncid,varid,f_tleaf) )
      
      ! shaded leaf temperature [K]
      !call sanity( nf90_inq_varid(ncid,'f_tlsha',varid) )
      !call sanity( nf90_put_var(ncid,varid,f_tlsha) )
      
      ! depth of water on foliage [mm]
      call sanity( nf90_inq_varid(ncid,'f_ldew',varid) )
      call sanity( nf90_put_var(ncid,varid,f_ldew) )
      
      ! snow cover, water equivalent [mm]
      call sanity( nf90_inq_varid(ncid,'f_scv',varid) )
      call sanity( nf90_put_var(ncid,varid,f_scv) )
      
      ! snow depth [meter]
      call sanity( nf90_inq_varid(ncid,'f_snowdp',varid) )
      call sanity( nf90_put_var(ncid,varid,f_snowdp) )
      
      ! fraction of snow cover on ground
      call sanity( nf90_inq_varid(ncid,'f_fsno',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fsno) )
      
      ! fraction of veg cover, excluding snow-covered veg [-]
      call sanity( nf90_inq_varid(ncid,'f_sigf',varid) )
      call sanity( nf90_put_var(ncid,varid,f_sigf) )
      
      ! leaf greenness
      call sanity( nf90_inq_varid(ncid,'f_green',varid) )
      call sanity( nf90_put_var(ncid,varid,f_green) )
      
      ! leaf area index
      call sanity( nf90_inq_varid(ncid,'f_lai',varid) )
      call sanity( nf90_put_var(ncid,varid,f_lai) )
      
      ! stem area index
      call sanity( nf90_inq_varid(ncid,'f_sai',varid) )
      call sanity( nf90_put_var(ncid,varid,f_sai) )
      
      ! averaged albedo direct 
      tmp(:,:,1) = f_alb(1,1,:,:) 
      tmp(:,:,2) = f_alb(2,1,:,:) 
      call sanity( nf90_inq_varid(ncid,'f_albd',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp) )
      
      ! averaged albedo diffuse 
      tmp(:,:,1) = f_alb(1,2,:,:) 
      tmp(:,:,2) = f_alb(2,2,:,:) 
      call sanity( nf90_inq_varid(ncid,'f_albi',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp) )
      
      ! averaged bulk surface emissivity
      call sanity( nf90_inq_varid(ncid,'f_emis',varid) )
      call sanity( nf90_put_var(ncid,varid,f_emis) )
      
      ! effective roughness [m]
      call sanity( nf90_inq_varid(ncid,'f_z0m',varid) )
      call sanity( nf90_put_var(ncid,varid,f_z0m) )
      
      ! radiative temperature of surface [K]
      call sanity( nf90_inq_varid(ncid,'f_trad',varid) )
      call sanity( nf90_put_var(ncid,varid,f_trad) )
      
      ! 2 m height air temperature [kelvin]
      call sanity( nf90_inq_varid(ncid,'f_tref',varid) )
      call sanity( nf90_put_var(ncid,varid,f_tref) )
      
      ! 2 m height air specific humidity [kg/kg]
      call sanity( nf90_inq_varid(ncid,'f_qref',varid) )
      call sanity( nf90_put_var(ncid,varid,f_qref) )
      
      ! rain [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_xy_rain',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xy_rain) )
 
      ! snow [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_xy_snow',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xy_snow) )
 
 
!---------------------------------------------------------------------
      ! soil temperature [K]
      do i = maxsnl+1, nl_soil
         tmp1(:,:,i) = f_t_soisno(i,:,:)
      end do
      call sanity( nf90_inq_varid(ncid,'f_t_soisno',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp1) )
      
      ! liquid water in soil layers [kg/m2]
      do i = maxsnl+1, nl_soil
         tmp1(:,:,i) = f_wliq_soisno(i,:,:)
      end do
      call sanity( nf90_inq_varid(ncid,'f_wliq_soisno',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp1) )
      
      ! ice lens in soil layers [kg/m2]
      do i = maxsnl+1, nl_soil
         tmp1(:,:,i) = f_wice_soisno(i,:,:)
      end do
      call sanity( nf90_inq_varid(ncid,'f_wice_soisno',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp1) )
 
      ! lake temperature [K]
      do i = 1, nl_lake
         tmp2(:,:,i) = f_t_lake(i,:,:)
      end do
      call sanity( nf90_inq_varid(ncid,'f_t_lake',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp2) )
      
      ! lake ice fraction cover [0-1]
      do i = 1, nl_lake
         tmp2(:,:,i) = f_lake_icefrac(i,:,:)
      end do
      call sanity( nf90_inq_varid(ncid,'f_lake_icefrac',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp2) )

      ! u* in similarity theory [m/s]
      call sanity( nf90_inq_varid(ncid,'f_ustar',varid) )
      call sanity( nf90_put_var(ncid,varid,f_ustar) )
      
      ! t* in similarity theory [kg/kg]
      call sanity( nf90_inq_varid(ncid,'f_tstar',varid) )
      call sanity( nf90_put_var(ncid,varid,f_tstar) )
      
      ! q* in similarity theory [kg/kg]
      call sanity( nf90_inq_varid(ncid,'f_qstar',varid) )
      call sanity( nf90_put_var(ncid,varid,f_qstar) )
      
      ! dimensionless height (z/L) used in Monin-Obukhov theory
      call sanity( nf90_inq_varid(ncid,'f_zol',varid) )
      call sanity( nf90_put_var(ncid,varid,f_zol) )
      
      ! bulk Richardson number in surface layer
      call sanity( nf90_inq_varid(ncid,'f_rib',varid) )
      call sanity( nf90_put_var(ncid,varid,f_rib) )
      
      ! integral of profile function for momentum
      call sanity( nf90_inq_varid(ncid,'f_fm',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fm) )
      
      ! integral of profile function for heat
      call sanity( nf90_inq_varid(ncid,'f_fh',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fh) )
      
      ! integral of profile function for moisture
      call sanity( nf90_inq_varid(ncid,'f_fq',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fq) )
      
      ! 10m u-velocity [m/s]
      call sanity( nf90_inq_varid(ncid,'f_us10m',varid) )
      call sanity( nf90_put_var(ncid,varid,f_us10m) )
      
      ! 10m v-velocity [m/s]
      call sanity( nf90_inq_varid(ncid,'f_vs10m',varid) )
      call sanity( nf90_put_var(ncid,varid,f_vs10m) )
      
      ! integral of profile function for momentum at 10m [-]
      call sanity( nf90_inq_varid(ncid,'f_fm10m',varid) )
      call sanity( nf90_put_var(ncid,varid,f_fm10m) )

      ! wind in eastward direction [m/s]
      call sanity( nf90_inq_varid(ncid,'f_xy_us',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xy_us) )
      
      ! wind in northward direction [m/s]
      call sanity( nf90_inq_varid(ncid,'f_xy_vs',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xy_vs) )
      
      ! temperature at reference height [kelvin]
      call sanity( nf90_inq_varid(ncid,'f_xy_t',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xy_t) )
      
      ! specific humidity at reference height [kg/kg]
      call sanity( nf90_inq_varid(ncid,'f_xy_q',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xy_q) )
      
      ! convective precipitation [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_xy_prc',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xy_prc) )
      
      ! large scale precipitation [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_xy_prl',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xy_prl) )
      
      ! atmospheric pressure at the surface [pa]
      call sanity( nf90_inq_varid(ncid,'f_xy_pbot',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xy_pbot) )
      
      ! atmospheric infrared (longwave) radiation [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_xy_frl',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xy_frl) )
      
      ! downward solar radiation at surface [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_xy_solarin',varid) )
      call sanity( nf90_put_var(ncid,varid,f_xy_solarin) )
      
      ! total reflected solar radiation at surface [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_sr',varid) )
      call sanity( nf90_put_var(ncid,varid,f_sr) )
      
      ! incident direct beam vis solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'f_solvd',varid) )
      call sanity( nf90_put_var(ncid,varid,f_solvd) )
      
      ! incident diffuse beam vis solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'f_solvi',varid) )
      call sanity( nf90_put_var(ncid,varid,f_solvi) )
      
      ! incident direct beam nir solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'f_solnd',varid) )
      call sanity( nf90_put_var(ncid,varid,f_solnd) )
      
      ! incident diffuse beam nir solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'f_solni',varid) )
      call sanity( nf90_put_var(ncid,varid,f_solni) )
      
      ! reflected direct beam vis solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'f_srvd',varid) )
      call sanity( nf90_put_var(ncid,varid,f_srvd) )
      
      ! reflected diffuse beam vis solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'f_srvi',varid) )
      call sanity( nf90_put_var(ncid,varid,f_srvi) )
      
      ! reflected direct beam nir solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'f_srnd',varid) )
      call sanity( nf90_put_var(ncid,varid,f_srnd) )
      
      ! reflected diffuse beam nir solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'f_srni',varid) )
      call sanity( nf90_put_var(ncid,varid,f_srni) )
      
      ! incident direct beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'f_solvdln',varid) )
      call sanity( nf90_put_var(ncid,varid,f_solvdln) )
      
      ! incident diffuse beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'f_solviln',varid) )
      call sanity( nf90_put_var(ncid,varid,f_solviln) )
      
      ! incident direct beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'f_solndln',varid) )
      call sanity( nf90_put_var(ncid,varid,f_solndln) )
      
      ! incident diffuse beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'f_solniln',varid) )
      call sanity( nf90_put_var(ncid,varid,f_solniln) )
      
      ! reflected direct beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'f_srvdln',varid) )
      call sanity( nf90_put_var(ncid,varid,f_srvdln) )
      
      ! reflected diffuse beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'f_srviln',varid) )
      call sanity( nf90_put_var(ncid,varid,f_srviln) )
      
      ! reflected direct beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'f_srndln',varid) )
      call sanity( nf90_put_var(ncid,varid,f_srndln) )
      
      ! reflected diffuse beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'f_srniln',varid) )
      call sanity( nf90_put_var(ncid,varid,f_srniln) )
 
      call sanity( nf90_close(ncid) )

   END SUBROUTINE writenetcdf

   SUBROUTINE writenetcdf_ncar
   
      implicit none

      integer  :: months(0:12)
      integer  :: ix,iy,ilev,i,j,iyear,iday,isec,imonth
      integer  :: ncid, xid, yid, varid, sslevid, lakelevid, bandid, timeid
      integer  :: lake_lev(nl_lake), sslev(nl_soil), band(2)
      real(r4) :: lons(lon_points), lats(lat_points)
      real(r4) :: vars(lon_points, lat_points)
      real(r4) :: tmp0(lon_points, lat_points)
      real(r4) :: tmp(lon_points, lat_points, 2)
      real(r4) :: tmp1(lon_points, lat_points, nl_soil)
      real(r4) :: tmp2(lon_points, lat_points, nl_lake)
      real(r4) :: calyear
   
      character(len=256) :: casename
      character(len=256) :: year, day, sec, month
      
    ! create netcdf and define dimensions
    ! need to rename the file (split GLOBAL_2D_Fluxes_1990-01/1990-001-84600)
      i = index(filename, '_2D_Fluxes_')
      j = index(filename, '-', back=.true.)
     
      read(filename(    :i-1 ), *) casename
      read(filename(i+11:i+14), *) year
      read(filename(i+11:i+14), "(I4)") iyear
      
      if( (mod(iyear,4)==0 .AND. mod(iyear,100)/=0) .OR. mod(iyear,400)==0 ) then
         months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
      else
         months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      end if

      !call sanity( nf90_create(trim(filename)//'.nc', nf90_64bit_offset, ncid) )
      if ( (j-i) == 15) then
         read(filename(j+1:j+2), *) month
         read(filename(j+1:j+2), "(I2)") imonth
         calyear = iyear + months(imonth)*1./months(12)
         call sanity( nf90_create(trim(casename)//'.clm2.h0.'//trim(year)//'-'//trim(month)//'.nc', &
                      nf90_64bit_offset, ncid) )
      else
         read(filename(j-3:j-1), *) day
         read(filename(j-3:j-1), "(I3)") iday
         read(filename(j+1:j+5), *) sec
         read(filename(j+1:j+5), "(I5)") isec
         calyear = iyear + (iday - 1. + isec/86400.)/months(12)
         call sanity( nf90_create(trim(casename)//'.clm2.h0.'//trim(year)//'-'//trim(day)//'-'//trim(sec)//'.nc', &
                      nf90_64bit_offset, ncid) )
      end if
      

    ! dimensions
      call sanity( nf90_def_dim(ncid, 'lon',           lon_points,      xid) )        ! lon
      call sanity( nf90_def_dim(ncid, 'lat',           lat_points,      yid) )        ! lat
      call sanity( nf90_def_dim(ncid, 'levgrnd',       nl_soil,         sslevid) )    ! levgrnd
      call sanity( nf90_def_dim(ncid, 'levlak',        nl_lake,         lakelevid) )  ! levlak
      call sanity( nf90_def_dim(ncid, 'numrad',        2,               bandid) )     ! numrad
      call sanity( nf90_def_dim(ncid, 'time',          nf90_unlimited,  timeid) )     ! time

    ! global attr
      call sanity( nf90_put_att(ncid, NF90_GLOBAL, 'title','CLM SIMULATED SURFACE FLUXES') )
      
    ! variables
    ! dimension variables
      call sanity( nf90_def_var(ncid, 'lon', nf90_float, (/xid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','longitude') )
      call sanity( nf90_put_att(ncid, varid, 'units','degrees_east') )

      call sanity( nf90_def_var(ncid, 'lat', nf90_float, (/yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','latitude') )
      call sanity( nf90_put_att(ncid, varid, 'units','degrees_north') )

      call sanity( nf90_def_var(ncid, 'levgrnd', nf90_int, (/sslevid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name',"soil plus snow level") )
      
      call sanity( nf90_def_var(ncid, 'levlak', nf90_int, (/lakelevid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name',"lake level") )
      
      call sanity( nf90_def_var(ncid, 'numrad', nf90_int, (/bandid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name',"band (vis/nir)") )
      
      call sanity( nf90_def_var(ncid, 'time', nf90_float, (/timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name',"time") )

      ! grid mask
      call sanity( nf90_def_var(ncid, 'landmask', nf90_int, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','grid mask') )
      call sanity( nf90_put_att(ncid, varid, 'units','none') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', -1) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', -1) )

      ! grid total fraction
      call sanity( nf90_def_var(ncid, 'landfrac', nf90_float, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','grid total fraction') )
      call sanity( nf90_put_att(ncid, varid, 'units','fraction') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! grid cell area [km2]
      call sanity( nf90_def_var(ncid, 'area', nf90_float, (/xid,yid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','grid cell area [km2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','fraction') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! wind stress: E-W [kg/m/s2], TAUX
      call sanity( nf90_def_var(ncid, 'TAUX', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','wind stress: E-W [kg/m/s2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/m/s2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! wind stress: N-S [kg/m/s2], TAUY
      call sanity( nf90_def_var(ncid, 'TAUY', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','wind stress: N-S [kg/m/s2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/m/s2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! sensible heat from canopy height to atmosphere [W/m2], FSH
      call sanity( nf90_def_var(ncid, 'FSH', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','sensible heat from canopy height to atmosphere [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! latent heat flux from canopy height to atmosphere [W/m2], FCTR + FCEV + FGEV
      call sanity( nf90_def_var(ncid, 'f_lfevpa', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','latent heat flux from canopy height to atmosphere [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! canopy transpiration
      call sanity( nf90_def_var(ncid, 'FCTR', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','canopy tranpiration [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! canopy evaporation
      call sanity( nf90_def_var(ncid, 'FCEV', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','canopy evaporation [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! ground evaporation
      call sanity( nf90_def_var(ncid, 'FGEV', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','ground evaporation [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! evapotranspiration from canopy height to atmosphere [mm/s], QVEGT + QVEGE + QSOIL
      call sanity( nf90_def_var(ncid, 'f_fevpa', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','evapotranspiration from canopy height to atmosphere [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! sensible heat from leaves [W/m2], FSH_V
      call sanity( nf90_def_var(ncid, 'FSH_V', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','sensible heat from leaves [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! evaporation+transpiration from leaves [mm/s], QVEGT + QVEGE
      call sanity( nf90_def_var(ncid, 'f_fevpl', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','evaporation+transpiration from leaves [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! evaporation from leaves [mm/s], QVEGE
      call sanity( nf90_def_var(ncid, 'QVEGE', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','evaporation from leaves [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! transpiration rate [mm/s], QVEGT
      call sanity( nf90_def_var(ncid, 'QVEGT', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','transpiration rate [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! sensible heat flux from ground [W/m2], FSH_G
      call sanity( nf90_def_var(ncid, 'FSH_G', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','sensible heat flux from ground [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! evaporation heat flux from ground [mm/s], QSOIL
      call sanity( nf90_def_var(ncid, 'QSOIL', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','evaporation heat flux from ground [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! ground heat flux [W/m2], FGR
      call sanity( nf90_def_var(ncid, 'FGR', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','ground heat flux [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! solar absorbed by sunlit canopy [W/m2], SUN_ATOT
      call sanity( nf90_def_var(ncid, 'SUN_ATOT', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','solar absorbed by sunlit canopy [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! solar absorbed by shaded [W/m2], SHA_ATOT
      call sanity( nf90_def_var(ncid, 'SHA_ATOT', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','solar absorbed by shaded [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! solar absorbed by vegetation  [W/m2], SABV
      call sanity( nf90_def_var(ncid, 'SABV', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','solar absorbed by vegetaion [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! solar absorbed by ground  [W/m2], SABG
      call sanity( nf90_def_var(ncid, 'SABG', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','solar absorbed by ground  [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! outgoing long-wave radiation from ground+canopy [W/m2], LWup/FIRE
      call sanity( nf90_def_var(ncid, 'FIRE', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','outgoing long-wave radiation from ground+canopy [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! net long-wave radiation [W/m2], FIRA
      call sanity( nf90_def_var(ncid, 'FIRA', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','net long-wave radiation [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! net radiation [W/m2], Rnet
      call sanity( nf90_def_var(ncid, 'Rnet', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','net radiation [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! the error of water banace [mm/s], ERRH2O
      call sanity( nf90_def_var(ncid, 'ERRH2O', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','the error of water banace [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! the error of energy balance [W/m2], ERRSEB
      call sanity( nf90_def_var(ncid, 'ERRSEB', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','the error of energy balance [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! surface runoff [mm/s], QOVER
      call sanity( nf90_def_var(ncid, 'QOVER', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','surface runoff [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! total runoff [mm/s], QOVER + QDRAI + QRGWL
      call sanity( nf90_def_var(ncid, 'f_rnof', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','total runoff [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! sub-surface drainage [mm/s], QDRAI
      call sanity( nf90_def_var(ncid, 'QDRAI', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','sub-surface drainage [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! surface runoff at glaciers, wetlands, lakes [mm/s], QRGWL
      call sanity( nf90_def_var(ncid, 'QRGWL', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','surface runoff at glaciers, wetlands, lakes (included in QDRAI) [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! interception [mm/s], QINTR
      call sanity( nf90_def_var(ncid, 'QINTR', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','interception [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! inflitration [mm/s], QINFL
      call sanity( nf90_def_var(ncid, 'QINFL', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','inflitration [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! total throughfall [mm/s], QDRIP
      call sanity( nf90_def_var(ncid, 'QDRIP', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','total throughfall [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! canopy assimilation rate [mol m-2 s-1], FPSN, mol -> umol
      call sanity( nf90_def_var(ncid, 'FPSN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','canopy assimilation rate [umol m-2 s-1]') )
      call sanity( nf90_put_att(ncid, varid, 'units','umol m-2 s-1') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! respiration (plant+soil) [mol m-2 s-1], NONE
      call sanity( nf90_def_var(ncid, 'f_respc', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','respiration (plant+soil) [mol m-2 s-1]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mol m-2 s-1') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! groundwater recharge rate [mm/s], QCHARGE
      call sanity( nf90_def_var(ncid, 'QCHARGE', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','groundwater recharge rate [mm/s] ') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      !---------------------------------------------------------------------
      ! ground surface temperature [K], TG
      call sanity( nf90_def_var(ncid, 'TG', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','ground surface temperature [K]') )
      call sanity( nf90_put_att(ncid, varid, 'units','K') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! sunlit leaf temperature [K], TV
      call sanity( nf90_def_var(ncid, 'TV', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','sunlit leaf temperature [K]') )
      call sanity( nf90_put_att(ncid, varid, 'units','K') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! shaded leaf temperature [K], TV
      !call sanity( nf90_def_var(ncid, 'f_tlsha', nf90_float, (/xid,yid,timeid/), varid) )
      !call sanity( nf90_put_att(ncid, varid, 'long_name','shaded leaf temperature [K]') )
      !call sanity( nf90_put_att(ncid, varid, 'units','K') )
      !call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      !call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! depth of water on foliage [mm], H2OCAN
      call sanity( nf90_def_var(ncid, 'H2OCAN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','depth of water on foliage [mm]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! snow cover, water equivalent [mm], H2OSNO
      call sanity( nf90_def_var(ncid, 'H2OSNO', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','snow cover, water equivalent [mm]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! snow liquid water [kg/m2], SNOWLIQ
      call sanity( nf90_def_var(ncid, 'SNOWLIQ', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','snow liquid water [kg/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! snow ice [kg/m2], SNOWICE
      call sanity( nf90_def_var(ncid, 'SNOWICE', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','snow ice [kg/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! snow depth [meter], SNOWDP
      call sanity( nf90_def_var(ncid, 'SNOWDP', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','snow depth [m]') )
      call sanity( nf90_put_att(ncid, varid, 'units','meter') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! fraction of snow cover on ground [-], FSNO
      call sanity( nf90_def_var(ncid, 'FSNO', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','fraction of snow cover on ground [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! fraction of veg cover, excluding snow-covered veg [-], NONE
      call sanity( nf90_def_var(ncid, 'f_sigf', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','fraction of veg cover, excluding snow-covered veg [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! leaf greenness [fraction], NONE
      call sanity( nf90_def_var(ncid, 'f_green', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','leaf greenness [fraction]') )
      call sanity( nf90_put_att(ncid, varid, 'units','fraction') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! leaf area index [m2/m2], ELAI
      call sanity( nf90_def_var(ncid, 'ELAI', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','leaf area index [m2/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m2/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! stem area index [m2/m2], ESAI
      call sanity( nf90_def_var(ncid, 'ESAI', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','stem area index [m2/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m2/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! leaf area index [m2/m2], TLAI
      call sanity( nf90_def_var(ncid, 'TLAI', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','leaf area index [m2/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m2/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! stem area index [m2/m2], TSAI
      call sanity( nf90_def_var(ncid, 'TSAI', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','stem area index [m2/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m2/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! sunlit leaf area index [m2/m2], LAISUN
      call sanity( nf90_def_var(ncid, 'LAISUN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','sunlit leaf area index [m2/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m2/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! shaded leaf area index [m2/m2], LAISUN
      call sanity( nf90_def_var(ncid, 'LAISHA', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','shaded leaf area index [m2/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m2/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! averaged albedo direct [%], ALBD
      call sanity( nf90_def_var(ncid, 'ALBD', nf90_float, (/xid,yid,bandid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','averaged albedo direct [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! averaged albedo diffuse [%], ALBI
      call sanity( nf90_def_var(ncid, 'ALBI', nf90_float, (/xid,yid,bandid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','averaged albedo diffuse [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! averaged bulk surface emissivity [-], NONE
      call sanity( nf90_def_var(ncid, 'f_emis', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','averaged bulk surface emissivity [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! effective roughness [m], Z0MG
      call sanity( nf90_def_var(ncid, 'Z0MG', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','effective roughness [m]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! radiative temperature of surface [K], NONE
      call sanity( nf90_def_var(ncid, 'f_trad', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','radiative temperature of surface [K]') )
      call sanity( nf90_put_att(ncid, varid, 'units','K') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! 2 m height air temperature [kelvin], TSA
      call sanity( nf90_def_var(ncid, 'TSA', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','2 m height air temperature [kelvin]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kelvin') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! 2 m height air specific humidity [kg/kg], Q2M
      call sanity( nf90_def_var(ncid, 'Q2M', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','2 m height air specific humidity [kg/kg]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/kg') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
 
      ! rain [mm/s], RAIN
      call sanity( nf90_def_var(ncid, 'RAIN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','rain [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
 
      ! snow [mm/s], SNOW
      call sanity( nf90_def_var(ncid, 'SNOW', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','snow [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
 
      !---------------------------------------------------------------------
      ! soil temperature [K], TSOI
      call sanity( nf90_def_var(ncid, 'TSOI', nf90_float, (/xid,yid,sslevid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','soil temperature [K]') )
      call sanity( nf90_put_att(ncid, varid, 'units','K') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! liquid water in soil layers [kg/m2], SOILLIQ
      call sanity( nf90_def_var(ncid, 'SOILLIQ', nf90_float, (/xid,yid,sslevid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','liquid water in soil layers [kg/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! ice lens in soil layers [kg/m2], SOILICE
      call sanity( nf90_def_var(ncid, 'SOILICE', nf90_float, (/xid,yid,sslevid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','ice lens in soil layers [kg/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! volumetric water in soil layers [m3/m3], H2OSOI
      call sanity( nf90_def_var(ncid, 'H2OSOI', nf90_float, (/xid,yid,sslevid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','volumetric water in soil layers [m3/m3]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m3/m3') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! factor of soil water stress, BTRAN
      call sanity( nf90_def_var(ncid, 'BTRAN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','factor of soil water stress') )
      call sanity( nf90_put_att(ncid, varid, 'units','0-1') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! the depth to water table [m], ZWT
      call sanity( nf90_def_var(ncid, 'ZWT', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','the depth to water table [m]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! water storage in aquifer [mm], WA
      call sanity( nf90_def_var(ncid, 'WA', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','water storage in aquifer [mm]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! total water storage [mm], WT
      call sanity( nf90_def_var(ncid, 'WT', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','total water storage [mm]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! lake temperature [K], TLAKE
      call sanity( nf90_def_var(ncid, 'TLAKE', nf90_float, (/xid,yid,lakelevid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','lake temperature [K]') )
      call sanity( nf90_put_att(ncid, varid, 'units','K') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! lake ice fraction cover [0-1], LAKEICEFRAC
      call sanity( nf90_def_var(ncid, 'LAKEICEFRAC', nf90_float, (/xid,yid,lakelevid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','lake ice fraction cover [0-1]') )
      call sanity( nf90_put_att(ncid, varid, 'units','0-1') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! u* in similarity theory [m/s], NONE
      call sanity( nf90_def_var(ncid, 'f_ustar', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','u* in similarity theory [m/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! t* in similarity theory [kg/kg], NONE
      call sanity( nf90_def_var(ncid, 'f_tstar', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','t* in similarity theory [kg/kg]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/kg') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! q* in similarity theory [kg/kg], NONE
      call sanity( nf90_def_var(ncid, 'f_qstar', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','q* in similarity theory [kg/kg]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/kg') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! dimensionless height (z/L) used in Monin-Obukhov theory [-], NONE
      call sanity( nf90_def_var(ncid, 'f_zol', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','dimensionless height (z/L) used in Monin-Obukhov theory [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! bulk Richardson number in surface layer [-], NONE
      call sanity( nf90_def_var(ncid, 'f_rib', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','bulk Richardson number in surface layer [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! integral of profile function for momentum [-], NONE
      call sanity( nf90_def_var(ncid, 'f_fm', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','integral of profile function for momentum [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! integral of profile function for heat [-], NONE
      call sanity( nf90_def_var(ncid, 'f_fh', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','integral of profile function for heat [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! integral of profile function for moisture [-], NONE
      call sanity( nf90_def_var(ncid, 'f_fq', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','integral of profile function for moisture [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! 10m u-velocity [m/s], U10
      call sanity( nf90_def_var(ncid, 'U10', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','10m u-velocity [m/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! 10m v-velocity [m/s], NONE
      call sanity( nf90_def_var(ncid, 'f_vs10m', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','10m v-velocity [m/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! integral of profile function for momentum at 10m [-], NONE
      call sanity( nf90_def_var(ncid, 'f_fm10m', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','integral of profile function for momentum at 10m [-]') )
      call sanity( nf90_put_att(ncid, varid, 'units','-') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! wind in eastward direction [m/s], WIND
      call sanity( nf90_def_var(ncid, 'WIND', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','wind in eastward direction [m/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! wind in northward direction [m/s], NONE
      call sanity( nf90_def_var(ncid, 'f_xy_vs', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','wind in northward direction [m/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','m/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! temperature at reference height [kelvin], TBOT
      call sanity( nf90_def_var(ncid, 'TBOT', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','temperature at reference height [kelvin]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kelvin') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! specific humidity at reference height [kg/kg], QBOT
      call sanity( nf90_def_var(ncid, 'QBOT', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','specific humidity at reference height [kg/kg]') )
      call sanity( nf90_put_att(ncid, varid, 'units','kg/kg') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! convective precipitation [mm/s], NONE
      call sanity( nf90_def_var(ncid, 'f_xy_prc', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','convective precipitation [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! large scale precipitation [mm/s], NONE
      call sanity( nf90_def_var(ncid, 'f_xy_prl', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','large scale precipitation [mm/s]') )
      call sanity( nf90_put_att(ncid, varid, 'units','mm/s') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! atmospheric pressure at the surface [pa], PSurf
      call sanity( nf90_def_var(ncid, 'PSurf', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','atmospheric pressure at the surface [pa]') )
      call sanity( nf90_put_att(ncid, varid, 'units','pa') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! atmospheric infrared (longwave) radiation [W/m2], LWdown/FLDS
      call sanity( nf90_def_var(ncid, 'FLDS', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','atmospheric infrared (longwave) radiation [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! downward solar radiation at surface [W/m2], FSDS
      call sanity( nf90_def_var(ncid, 'FSDS', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','downward solar radiation at surface [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! reflected solar radiation at surface [W/m2], FSR
      call sanity( nf90_def_var(ncid, 'FSR', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected solar radiation at surface [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! absorbed solar radiation at surface [W/m2], FSA
      call sanity( nf90_def_var(ncid, 'FSA', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','absorbed solar radiation at surface [W/m2]') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )

      ! incident direct beam vis solar radiation (W/m2), FSDSVD
      call sanity( nf90_def_var(ncid, 'FSDSVD', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident direct beam vis solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! incident diffuse beam vis solar radiation (W/m2), FSDSVI
      call sanity( nf90_def_var(ncid, 'FSDSVI', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam vis solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! incident direct beam nir solar radiation (W/m2), FSDSND
      call sanity( nf90_def_var(ncid, 'FSDSND', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident direct beam nir solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! incident diffuse beam nir solar radiation (W/m2), FSDSNI
      call sanity( nf90_def_var(ncid, 'FSDSNI', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam nir solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! reflected direct beam vis solar radiation (W/m2), FSRVD
      call sanity( nf90_def_var(ncid, 'FSRVD', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected direct beam vis solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! reflected diffuse beam vis solar radiation (W/m2), FSRVI
      call sanity( nf90_def_var(ncid, 'FSRVI', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam vis solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! reflected direct beam nir solar radiation (W/m2), FSRND
      call sanity( nf90_def_var(ncid, 'FSRND', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected direct beam nir solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! reflected diffuse beam nir solar radiation (W/m2), FSRNI
      call sanity( nf90_def_var(ncid, 'FSRNI', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam nir solar radiation (W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! incident direct beam vis solar radiation at local noon(W/m2), FSDSVDLN
      call sanity( nf90_def_var(ncid, 'FSDSVDLN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident direct beam vis solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! incident diffuse beam vis solar radiation at local noon(W/m2), FSDSVILN
      call sanity( nf90_def_var(ncid, 'FSDSVILN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam vis solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! incident direct beam nir solar radiation at local noon(W/m2), FSDSNDLN
      call sanity( nf90_def_var(ncid, 'FSDSNDLN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident direct beam nir solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! incident diffuse beam nir solar radiation at local noon(W/m2), FSDSNILN
      call sanity( nf90_def_var(ncid, 'FSDSNILN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam nir solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! reflected direct beam vis solar radiation at local noon(W/m2), FSRVDLN
      call sanity( nf90_def_var(ncid, 'FSRVDLN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected direct beam vis solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! reflected diffuse beam vis solar radiation at local noon(W/m2), FSRVILN
      call sanity( nf90_def_var(ncid, 'FSRVILN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam vis solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! reflected direct beam nir solar radiation at local noon(W/m2), FSRNDLN
      call sanity( nf90_def_var(ncid, 'FSRNDLN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected direct beam nir solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
      ! reflected diffuse beam nir solar radiation at local noon(W/m2), FSRNILN
      call sanity( nf90_def_var(ncid, 'FSRNILN', nf90_float, (/xid,yid,timeid/), varid) )
      call sanity( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam nir solar radiation at local noon(W/m2)') )
      call sanity( nf90_put_att(ncid, varid, 'units','W/m2') )
      call sanity( nf90_put_att(ncid, varid, 'missing_value', spval_r4) )
      call sanity( nf90_put_att(ncid, varid, '_FillValue', spval_r4) )
      
    ! end defination
      call sanity( nf90_enddef(ncid) )

    ! write data
    ! ------------------------------------------------------------

    ! dimension data
      lons = lons_r8 + 180.
      call sanity( nf90_inq_varid(ncid,'lon',varid) )
      call sanity( nf90_put_var(ncid,varid,lons) )

      lats = lats_r8(lat_points:1:-1)
      call sanity( nf90_inq_varid(ncid,'lat',varid) )
      call sanity( nf90_put_var(ncid,varid,lats) )

      do ilev = 1, nl_soil
         sslev(ilev) = ilev
      end do
      call sanity( nf90_inq_varid(ncid,'levgrnd',varid) )
      call sanity( nf90_put_var(ncid,varid,sslev) )

      do ilev = 1, nl_lake
         lake_lev(ilev) = ilev
      end do
      call sanity( nf90_inq_varid(ncid,'levlak',varid) )
      call sanity( nf90_put_var(ncid,varid,lake_lev) )

      do ilev = 1, 2
         band(ilev) = ilev
      end do
      call sanity( nf90_inq_varid(ncid,'numrad',varid) )
      call sanity( nf90_put_var(ncid,varid,band) )
      
      call sanity( nf90_inq_varid(ncid,'time',varid) )
      call sanity( nf90_put_var(ncid,varid,calyear) )

! added by yuan, 07/10/2016
! bug found, coordinates adjust missing.
      vars(1:max(lon_points/2,1),:) = mask((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = mask(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! grid mask
      call sanity( nf90_inq_varid(ncid,'landmask',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = frac((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = frac(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! grid total fraction
      call sanity( nf90_inq_varid(ncid,'landfrac',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = area((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = area(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! grid cell area [km2]
      call sanity( nf90_inq_varid(ncid,'area',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_taux((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_taux(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! wind stress: E-W [kg/m/s2]
      call sanity( nf90_inq_varid(ncid,'TAUX',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_tauy((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_tauy(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! wind stress: N-S [kg/m/s2]
      call sanity( nf90_inq_varid(ncid,'TAUY',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_fsena((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fsena(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! sensible heat from canopy height to atmosphere [W/m2]
      call sanity( nf90_inq_varid(ncid,'FSH',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_lfevpa((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_lfevpa(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! latent heat flux from canopy height to atmosphere [W/m2]
      call sanity( nf90_inq_varid(ncid,'f_lfevpa',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
     
      vars(1:max(lon_points/2,1),:) = f_fevpa((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fevpa(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! evapotranspiration from canopy to atmosphere [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_fevpa',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_fsenl((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fsenl(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! sensible heat from leaves [W/m2]
      call sanity( nf90_inq_varid(ncid,'FSH_V',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_fevpl((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fevpl(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! evaporation+transpiration from leaves [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_fevpl',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_etr((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_etr(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! transpiration rate [mm/s]
      call sanity( nf90_inq_varid(ncid,'QVEGT',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      where (vars /= spval_r4) vars = hvap * vars
      ! canopy transpiration [W/m2]
      call sanity( nf90_inq_varid(ncid,'FCTR',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      tmp0 = spval_r4
      where (f_fevpl /= spval) tmp0 = f_fevpl - f_etr

      vars(1:max(lon_points/2,1),:) = tmp0((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = tmp0(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! evaporation from leaves [mm/s]
      call sanity( nf90_inq_varid(ncid,'QVEGE',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      where (vars /= spval_r4) vars = hvap * vars
      ! canopy evaporation [W/m2]
      call sanity( nf90_inq_varid(ncid,'FCEV',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_fseng((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fseng(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! sensible heat flux from ground [W/m2]
      call sanity( nf90_inq_varid(ncid,'FSH_G',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_fevpg((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fevpg(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! evaporation heat flux from ground [mm/s]
      call sanity( nf90_inq_varid(ncid,'QSOIL',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      where (vars /= spval_r4) vars = hvap * vars
      ! ground evaporation [W/m2]
      call sanity( nf90_inq_varid(ncid,'FGEV',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_fgrnd((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fgrnd(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! ground heat flux [W/m2]
      call sanity( nf90_inq_varid(ncid,'FGR',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_sabvsun((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_sabvsun(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! solar absorbed by sunlit canopy [W/m2]
      call sanity( nf90_inq_varid(ncid,'SUN_ATOT',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_sabvsha((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_sabvsha(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! solar absorbed by shaded [W/m2]
      call sanity( nf90_inq_varid(ncid,'SHA_ATOT',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      tmp0 = spval_r4
      where (f_sabvsun /= spval) tmp0 = f_sabvsun + f_sabvsha
      
      vars(1:max(lon_points/2,1),:) = tmp0((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = tmp0(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! solar absorbed by vegetation  [W/m2]
      call sanity( nf90_inq_varid(ncid,'SABV',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
         
      vars(1:max(lon_points/2,1),:) = f_sabg((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_sabg(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! solar absorbed by ground  [W/m2]
      call sanity( nf90_inq_varid(ncid,'SABG',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_olrg((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_olrg(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! outgoing long-wave radiation from ground+canopy [W/m2]
      call sanity( nf90_inq_varid(ncid,'FIRE',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      tmp0 = spval_r4
      where (f_olrg /= spval) tmp0 = f_olrg - f_xy_frl
      
      vars(1:max(lon_points/2,1),:) = tmp0((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = tmp0(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! net long-wave radiation [W/m2]
      call sanity( nf90_inq_varid(ncid,'FIRA',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
         
      vars(1:max(lon_points/2,1),:) = f_rnet((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_rnet(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! net radiation [W/m2]
      call sanity( nf90_inq_varid(ncid,'Rnet',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_xerr((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xerr(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! the error of water banace [mm/s]
      call sanity( nf90_inq_varid(ncid,'ERRH2O',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_zerr((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_zerr(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! the error of energy balance [W/m2]
      call sanity( nf90_inq_varid(ncid,'ERRSEB',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_rsur((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_rsur(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! surface runoff [mm/s]
      call sanity( nf90_inq_varid(ncid,'QOVER',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_rnof((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_rnof(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! total runoff [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_rnof',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      tmp0 = spval_r4
      where (f_rnof /= spval) tmp0 = f_rnof - f_rsur

      vars(1:max(lon_points/2,1),:) = tmp0((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = tmp0(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! sub-surface drainage  [mm/s]
      call sanity( nf90_inq_varid(ncid,'QDRAI',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      where (vars /= spval) vars = 0.
      ! surface runoff at glaciers, wetlands, lakes  [mm/s]
      call sanity( nf90_inq_varid(ncid,'QRGWL',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_qintr((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_qintr(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! inflitration [mm/s]
      call sanity( nf90_inq_varid(ncid,'QINTR',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_qinfl((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_qinfl(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! inflitration [mm/s]
      call sanity( nf90_inq_varid(ncid,'QINFL',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_qdrip((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_qdrip(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! throughfall [mm/s]
      call sanity( nf90_inq_varid(ncid,'QDRIP',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      tmp0 = spval_r4
      where (f_assim /= spval) tmp0 = f_assim * 1.e6_r8

      vars(1:max(lon_points/2,1),:) = tmp0((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = tmp0(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! canopy assimilation rate [mol m-2 s-1]
      call sanity( nf90_inq_varid(ncid,'FPSN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_respc((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_respc(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! respiration (plant+soil) [mol m-2 s-1]
      call sanity( nf90_inq_varid(ncid,'f_respc',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_qcharge((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_qcharge(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! groundwater recharge rate [mm/s] 
      call sanity( nf90_inq_varid(ncid,'QCHARGE',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
!---------------------------------------------------------------------
      vars(1:max(lon_points/2,1),:) = f_t_grnd((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_t_grnd(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! ground surface temperature [K]
      call sanity( nf90_inq_varid(ncid,'TG',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_tleaf((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_tleaf(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! sunlit leaf temperature [K]
      call sanity( nf90_inq_varid(ncid,'TV',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      !vars(1:max(lon_points/2,1),:) = f_tlsha((lon_points/2+1):lon_points,:)
      !vars((lon_points/2+1):lon_points,:) = f_tlsha(1:max(lon_points/2,1),:)
      !vars = vars(:,lat_points:1:-1)
      !! shaded leaf temperature [K]
      !call sanity( nf90_inq_varid(ncid,'f_tlsha',varid) )
      !call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_ldew((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_ldew(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! depth of water on foliage [mm]
      call sanity( nf90_inq_varid(ncid,'H2OCAN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_scv((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_scv(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! snow cover, water equivalent [mm]
      call sanity( nf90_inq_varid(ncid,'H2OSNO',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      tmp0 = f_wliq_soisno(0,:,:)
      do i = maxsnl+1, -1
         where (f_wliq_soisno(i,:,:) /= spval) tmp0 = tmp0 + f_wliq_soisno(i,:,:)
      end do
      
      vars(1:max(lon_points/2,1),:) = tmp0((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = tmp0(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! snow liquid water [kg/m2]
      call sanity( nf90_inq_varid(ncid,'SNOWLIQ',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      tmp0 = f_wice_soisno(0,:,:)
      do i = maxsnl+1, -1
         where (f_wice_soisno(i,:,:) /= spval) tmp0 = tmp0 + f_wice_soisno(i,:,:)
      end do
      
      vars(1:max(lon_points/2,1),:) = tmp0((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = tmp0(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! snow ice [kg/m2]
      call sanity( nf90_inq_varid(ncid,'SNOWICE',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_snowdp((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_snowdp(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! snow depth [meter]
      call sanity( nf90_inq_varid(ncid,'SNOWDP',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_fsno((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fsno(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! fraction of snow cover on ground
      call sanity( nf90_inq_varid(ncid,'FSNO',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_sigf((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_sigf(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! fraction of veg cover, excluding snow-covered veg [-]
      call sanity( nf90_inq_varid(ncid,'f_sigf',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_green((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_green(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! leaf greenness
      call sanity( nf90_inq_varid(ncid,'f_green',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_lai((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_lai(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! leaf area index
      call sanity( nf90_inq_varid(ncid,'ELAI',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      ! leaf area index
      call sanity( nf90_inq_varid(ncid,'TLAI',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_sai((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_sai(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! stem area index
      call sanity( nf90_inq_varid(ncid,'ESAI',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      ! stem area index
      call sanity( nf90_inq_varid(ncid,'TSAI',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_laisun((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_laisun(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! sunlit leaf area index
      call sanity( nf90_inq_varid(ncid,'LAISUN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_laisha((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_laisha(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! shaded leaf area index
      call sanity( nf90_inq_varid(ncid,'LAISHA',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      ! averaged albedo direct 
      tmp(1:max(lon_points/2,1),:,1) = f_alb(1,1,(lon_points/2+1):lon_points,:)
      tmp((lon_points/2+1):lon_points,:,1) = f_alb(1,1,1:max(lon_points/2,1),:)
      tmp(1:max(lon_points/2,1),:,2) = f_alb(2,1,(lon_points/2+1):lon_points,:)
      tmp((lon_points/2+1):lon_points,:,2) = f_alb(2,1,1:max(lon_points/2,1),:)
      tmp = tmp(:,lat_points:1:-1,:)
      call sanity( nf90_inq_varid(ncid,'ALBD',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp) )
      
      ! averaged albedo diffuse 
      tmp(1:max(lon_points/2,1),:,1) = f_alb(1,2,(lon_points/2+1):lon_points,:)
      ! 08/22/2021, yuan: bug
      !tmp((lon_points/2+1):lon_points,:,1) = f_alb(1,1,1:max(lon_points/2,1),:)
      tmp((lon_points/2+1):lon_points,:,1) = f_alb(1,2,1:max(lon_points/2,1),:)
      tmp(1:max(lon_points/2,1),:,2) = f_alb(2,2,(lon_points/2+1):lon_points,:)
      ! 08/22/2021, yuan: bug
      !tmp((lon_points/2+1):lon_points,:,2) = f_alb(2,1,1:max(lon_points/2,1),:)
      tmp((lon_points/2+1):lon_points,:,2) = f_alb(2,2,1:max(lon_points/2,1),:)
      tmp = tmp(:,lat_points:1:-1,:)
      call sanity( nf90_inq_varid(ncid,'ALBI',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp) )
      
      vars(1:max(lon_points/2,1),:) = f_emis((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_emis(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! averaged bulk surface emissivity
      call sanity( nf90_inq_varid(ncid,'f_emis',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_z0m((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_z0m(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! effective roughness [m]
      call sanity( nf90_inq_varid(ncid,'Z0MG',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_trad((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_trad(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! radiative temperature of surface [K]
      call sanity( nf90_inq_varid(ncid,'f_trad',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_tref((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_tref(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! 2 m height air temperature [kelvin]
      call sanity( nf90_inq_varid(ncid,'TSA',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_qref((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_qref(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! 2 m height air specific humidity [kg/kg]
      call sanity( nf90_inq_varid(ncid,'Q2M',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
 
      vars(1:max(lon_points/2,1),:) = f_xy_rain((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xy_rain(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! rain [mm/s]
      call sanity( nf90_inq_varid(ncid,'RAIN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
 
      vars(1:max(lon_points/2,1),:) = f_xy_snow((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xy_snow(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! snow [mm/s]
      call sanity( nf90_inq_varid(ncid,'SNOW',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
 
!---------------------------------------------------------------------
      ! soil temperature [K]
      do i = 1, nl_soil
         tmp1(1:max(lon_points/2,1),:,i) = f_t_soisno(i,(lon_points/2+1):lon_points,:)
         tmp1((lon_points/2+1):lon_points,:,i) = f_t_soisno(i,1:max(lon_points/2,1),:)
      end do
      tmp1 = tmp1(:,lat_points:1:-1,:)
      call sanity( nf90_inq_varid(ncid,'TSOI',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp1) )
      
      ! liquid water in soil layers [kg/m2]
      do i = 1, nl_soil
         tmp1(1:max(lon_points/2,1),:,i) = f_wliq_soisno(i,(lon_points/2+1):lon_points,:)
         tmp1((lon_points/2+1):lon_points,:,i) = f_wliq_soisno(i,1:max(lon_points/2,1),:)
      end do
      tmp1 = tmp1(:,lat_points:1:-1,:)
      call sanity( nf90_inq_varid(ncid,'SOILLIQ',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp1) )
      
      ! ice lens in soil layers [kg/m2]
      do i = 1, nl_soil
         tmp1(1:max(lon_points/2,1),:,i) = f_wice_soisno(i,(lon_points/2+1):lon_points,:)
         tmp1((lon_points/2+1):lon_points,:,i) = f_wice_soisno(i,1:max(lon_points/2,1),:)
      end do
      tmp1 = tmp1(:,lat_points:1:-1,:)
      call sanity( nf90_inq_varid(ncid,'SOILICE',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp1) )
 
      ! volumetric soil in soil layers [m3/m3]
      do i = 1, nl_soil
         tmp1(1:max(lon_points/2,1),:,i) = f_h2osoi(i,(lon_points/2+1):lon_points,:)
         tmp1((lon_points/2+1):lon_points,:,i) = f_h2osoi(i,1:max(lon_points/2,1),:)
      end do
      tmp1 = tmp1(:,lat_points:1:-1,:)
      call sanity( nf90_inq_varid(ncid,'H2OSOI',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp1) )
 
      vars(1:max(lon_points/2,1),:) = f_rstfac((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_rstfac(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! factor of soil water stress [m/s]
      call sanity( nf90_inq_varid(ncid,'BTRAN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_zwt((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_zwt(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! soil water depth [m/s]
      call sanity( nf90_inq_varid(ncid,'ZWT',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_wa((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_wa(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! water storage in aquifer [m/s]
      call sanity( nf90_inq_varid(ncid,'WA',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_wat((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_wat(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! total water storage [mm]
      call sanity( nf90_inq_varid(ncid,'WT',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      ! lake temperature [K]
      do i = 1, nl_lake
         tmp2(1:max(lon_points/2,1),:,i) = f_t_lake(i,(lon_points/2+1):lon_points,:)
         tmp2((lon_points/2+1):lon_points,:,i) = f_t_lake(i,1:max(lon_points/2,1),:)
      end do
      tmp2 = tmp2(:,lat_points:1:-1,:)
      call sanity( nf90_inq_varid(ncid,'TLAKE',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp2) )
      
      ! lake ice fraction cover [0-1]
      do i = 1, nl_lake
         tmp2(1:max(lon_points/2,1),:,i) = f_lake_icefrac(i,(lon_points/2+1):lon_points,:)
         tmp2((lon_points/2+1):lon_points,:,i) = f_lake_icefrac(i,1:max(lon_points/2,1),:)
      end do
      tmp2 = tmp2(:,lat_points:1:-1,:)
      call sanity( nf90_inq_varid(ncid,'LAKEICEFRAC',varid) )
      call sanity( nf90_put_var(ncid,varid,tmp2) )
      
      vars(1:max(lon_points/2,1),:) = f_ustar((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_ustar(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! u* in similarity theory [m/s]
      call sanity( nf90_inq_varid(ncid,'f_ustar',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_tstar((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_tstar(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! t* in similarity theory [kg/kg]
      call sanity( nf90_inq_varid(ncid,'f_tstar',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_qstar((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_qstar(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! q* in similarity theory [kg/kg]
      call sanity( nf90_inq_varid(ncid,'f_qstar',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_zol((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_zol(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! dimensionless height (z/L) used in Monin-Obukhov theory
      call sanity( nf90_inq_varid(ncid,'f_zol',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_rib((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_rib(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! bulk Richardson number in surface layer
      call sanity( nf90_inq_varid(ncid,'f_rib',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_fm((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fm(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! integral of profile function for momentum
      call sanity( nf90_inq_varid(ncid,'f_fm',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_fh((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fh(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! integral of profile function for heat
      call sanity( nf90_inq_varid(ncid,'f_fh',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_fq((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fq(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! integral of profile function for moisture
      call sanity( nf90_inq_varid(ncid,'f_fq',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_us10m((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_us10m(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! 10m u-velocity [m/s]
      call sanity( nf90_inq_varid(ncid,'U10',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_vs10m((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_vs10m(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! 10m v-velocity [m/s]
      call sanity( nf90_inq_varid(ncid,'f_vs10m',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_fm10m((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_fm10m(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! integral of profile function for momentum at 10m [-]
      call sanity( nf90_inq_varid(ncid,'f_fm10m',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )

      vars(1:max(lon_points/2,1),:) = f_xy_us((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xy_us(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! wind in eastward direction [m/s]
      call sanity( nf90_inq_varid(ncid,'WIND',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_xy_vs((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xy_vs(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! wind in northward direction [m/s]
      call sanity( nf90_inq_varid(ncid,'f_xy_vs',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_xy_t((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xy_t(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! temperature at reference height [kelvin]
      call sanity( nf90_inq_varid(ncid,'TBOT',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_xy_q((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xy_q(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! specific humidity at reference height [kg/kg]
      call sanity( nf90_inq_varid(ncid,'QBOT',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_xy_prc((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xy_prc(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! convective precipitation [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_xy_prc',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_xy_prl((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xy_prl(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! large scale precipitation [mm/s]
      call sanity( nf90_inq_varid(ncid,'f_xy_prl',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_xy_pbot((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xy_pbot(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! atmospheric pressure at the surface [pa]
      call sanity( nf90_inq_varid(ncid,'PSurf',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_xy_frl((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xy_frl(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! atmospheric infrared (longwave) radiation [W/m2]
      call sanity( nf90_inq_varid(ncid,'FLDS',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_xy_solarin((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_xy_solarin(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! downward solar radiation at surface [W/m2]
      call sanity( nf90_inq_varid(ncid,'FSDS',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_sr((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_sr(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! total reflected solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'FSR',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      tmp0 = spval_r4
      where (f_sr /= spval) tmp0 = f_xy_solarin - f_sr

      vars(1:max(lon_points/2,1),:) = tmp0((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = tmp0(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! total reflected solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'FSA',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_solvd((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_solvd(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! incident direct beam vis solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'FSDSVD',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_solvi((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_solvi(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! incident diffuse beam vis solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'FSDSVI',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_solnd((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_solnd(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! incident direct beam nir solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'FSDSND',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_solni((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_solni(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! incident diffuse beam nir solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'FSDSNI',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_srvd((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_srvd(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! reflected direct beam vis solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'FSRVD',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_srvi((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_srvi(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! reflected diffuse beam vis solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'FSRVI',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_srnd((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_srnd(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! reflected direct beam nir solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'FSRND',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_srni((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_srni(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! reflected diffuse beam nir solar radiation (W/m2)
      call sanity( nf90_inq_varid(ncid,'FSRNI',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_solvdln((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_solvdln(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! incident direct beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'FSDSVDLN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_solviln((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_solviln(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! incident diffuse beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'FSDSVILN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_solndln((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_solndln(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! incident direct beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'FSDSNDLN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_solniln((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_solniln(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! incident diffuse beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'FSDSNILN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_srvdln((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_srvdln(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! reflected direct beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'FSRVDLN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_srviln((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_srviln(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! reflected diffuse beam vis solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'FSRVILN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_srndln((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_srndln(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! reflected direct beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'FSRNDLN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      vars(1:max(lon_points/2,1),:) = f_srniln((lon_points/2+1):lon_points,:)
      vars((lon_points/2+1):lon_points,:) = f_srniln(1:max(lon_points/2,1),:)
      vars = vars(:,lat_points:1:-1)
      ! reflected diffuse beam nir solar radiation at local noon(W/m2)
      call sanity( nf90_inq_varid(ncid,'FSRNILN',varid) )
      call sanity( nf90_put_var(ncid,varid,vars) )
      
      call sanity( nf90_close(ncid) )

   END SUBROUTINE writenetcdf_ncar


 ! nc operation check
 ! ------------------------------------------------------------
   SUBROUTINE sanity(ret)

      implicit none
      integer, intent(in) :: ret

      if (ret .ne. nf90_noerr) then
         write(6, *) trim(nf90_strerror(ret)); stop
      end if

   END SUBROUTINE sanity

END PROGRAM bin2netcdf
