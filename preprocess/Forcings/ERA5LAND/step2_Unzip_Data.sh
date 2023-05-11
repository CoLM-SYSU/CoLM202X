#!/bin/zsh
# code for unzip ERA5LAND data
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

DIR=./
SYear=1979
EYear=1989
varnames="volumetric_soil_water_layer_4"
Month1="01 02 03 04 05 06 07 08 09 10 11 12"
while [ ${SYear} -le ${EYear} ] ; do
for Month in ${Month1};do
for varname in ${varnames}; do
unzip ERA5LAND_${SYear}_${Month}_${varname}.nc.zip
mv data.nc ERA5LAND_${SYear}_${Month}_${varname}.nc
done
done
SYear=`expr $SYear + 1`
done

