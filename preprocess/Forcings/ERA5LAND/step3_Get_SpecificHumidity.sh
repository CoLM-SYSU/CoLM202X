#!/bin/zsh
# code for calculating Specific humidity from dew temperature and pressure
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

Years="2005"
Months="01 02 03 04 05 06 07 08 09 10 11 12"

for Year in ${Years}; do

for Month in ${Months}; do

cat <<EOF >q.ncl
;-write netcdf 
begin
  outfile = "ERA5LAND_${Year}_${Month}_specific_humidity.nc"
  if (isfilepresent(outfile)) then
    system("rm -rf "+outfile)          ;-- make sure that file does not exist
  end if
;-- open data file
  fin  = addfile("../../data/2m_dewpoint_temperature/ERA5LAND_${Year}_${Month}_2m_dewpoint_temperature.nc4","r")  ;-- open data file
  fin1 = addfile("../../data/surface_pressure/ERA5LAND_${Year}_${Month}_surface_pressure.nc4","r")  ;-- open data file
;-- get variable t and its dimensions and dimension sizes
  tdk    =  short2flt(fin->d2m)                     ;-- get variable
  p    =  short2flt(fin1->sp)                   ;-- get variable
  time   =  fin->time                             ;-- get dimension time
  latitude    =  fin->latitude                    ;-- get dimension lat
  longitude    =  fin->longitude                  ;-- get dimension lon
  ntim   =  dimsizes(time)                        ;-- get dimension sizes of time
  nlat   =  dimsizes(latitude)                    ;-- get dimension sizes of lat
  nlon   =  dimsizes(longitude)                   ;-- get dimension sizes of lon
  SP = mixhum_ptd (p, tdk, 2)
  SP@units     = "kg/kg"             ;-- define new units



;-- create new netCDF file
  fout = addfile(outfile,"c")

;-- begin output file settings
  setfileoption(fout,"DefineMode",True) ;-- explicitly declare file definition mode

;-- create global attributes of the file
  fAtt                  =  True        ;-- assign file attributes
  fAtt@title            = "Specific humidity"  
  fAtt@source_file      = "temp_${Year}_${Month}.nc"
  fAtt@Conventions      = "CF"   
  fAtt@creation_date    =  systemfunc ("date")        
  fAtt@history          =  "NCL script: ex_netcdf_2.ncl" 
  fAtt@comment          = "calculate Specific humidity from dew temperature and pressure usine NCL"       
  fileattdef(fout,fAtt)                ;-- copy file attributes    

;-- predefine the coordinate variables and their dimensionality
  dimNames = (/"time", "latitude", "longitude"/)  
  dimSizes = (/ -1   , nlat,  nlon /) 
  dimUnlim = (/ True , False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

;-- predefine the the dimensionality of the variables to be written out
  filevardef(fout, "time" ,typeof(time),getvardims(time)) 
  filevardef(fout, "latitude"  ,typeof(latitude), getvardims(latitude))                          
  filevardef(fout, "longitude"  ,typeof(longitude), getvardims(longitude))                          
  filevardef(fout, "SP" ,typeof(p),  getvardims(p))

;-- copy attributes associated with each variable to the file
  filevarattdef(fout,"time" ,time)       ;-- copy time attributes
  filevarattdef(fout,"latitude",latitude)        ;-- copy lat attributes
  filevarattdef(fout,"longitude"  ,longitude)        ;-- copy lon attributes
  filevarattdef(fout,"SP",   SP)     ;-- copy wdir attributes

;-- explicitly exit file definition mode (not required)
  setfileoption(fout,"DefineMode",False)

;-- output only the data values since the dimensionality and such have been predefined.
;-- locations on the file.
  fout->time   =  (/time/)               ;-- write time to new netCDF file
  fout->latitude    =  (/latitude/)                ;-- write lat to new netCDF file
  fout->longitude    =  (/longitude/)                ;-- write lon to new netCDF file
  fout->SP     =  (/SP/)                 ;-- write variable to new netCDF file
end
EOF
ncl q.ncl
rm q.ncl
cdo -b F32 -f nc4c -z zip_3 chname,SP,Q ERA5LAND_${Year}_${Month}_specific_humidity.nc ERA5LAND_${Year}_${Month}_specific_humidity.nc4
rm ERA5LAND_${Year}_${Month}_specific_humidity.nc
done
done
