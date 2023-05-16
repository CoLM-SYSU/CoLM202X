#!/bin/zsh

# code for unzip WFDE5 data
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

DIR=./
SYear=1979
EYear=2019
varnames="rainfall_flux snowfall_flux"
#varnames="mean_surface_downward_long_wave_radiation_flux"
#varnames="near_surface_wind_speed surface_air_pressure surface_downwelling_longwave_radiation surface_downwelling_shortwave_radiation  near_surface_specific_humidity grid_point_altitude near_surface_air_temperature "
Month1="01 02 03 04 05 06 07 08 09 10 11 12"
#Month1="10"
while [ ${SYear} -le ${EYear} ] ; do
#mkdir -p ${varnames}
for Month in ${Month1};do
for varname in ${varnames}; do
tar -xvzf ${DIR}/${varname}_CRU/WFDE5_${SYear}_${Month}_${varname}.tar.gz -C ${DIR}/${varname}_CRU/
tar -xvzf ${DIR}/${varname}_CRU_GPCP/WFDE5_${SYear}_${Month}_${varname}.tar.gz -C ${DIR}/${varname}_CRU_GPCP/

#mv 100m_uv_component_of_wind_${Year}_${Month}.nc ${Year}
done
done
SYear=`expr $SYear + 1`
done
