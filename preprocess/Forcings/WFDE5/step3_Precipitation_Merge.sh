#!/bin/zsh
DIR=./
SYear=1979
EYear=2019
Month1="01 02 03 04 05 06 07 08 09 10 11 12"
while [ ${SYear} -le ${EYear} ] ; do
for Month in ${Month1};do
 cdo merge rainfall_flux_CRU_GPCP/Rainf_WFDE5_CRU+GPCC_${SYear}${Month}_v2.0.nc snowfall_flux_CRU_GPCP/Snowf_WFDE5_CRU+GPCC_${SYear}${Month}_v2.0.nc  temp.nc
 cdo expr,'precipitation=Rainf+Snowf' temp.nc rainfall_plus_snowfall_flux_CRU_GPCP/precipitation_CRU_GPCP_${SYear}${Month}_v2.0.nc
 rm -rf temp.nc
done
SYear=`expr $SYear + 1`
done
