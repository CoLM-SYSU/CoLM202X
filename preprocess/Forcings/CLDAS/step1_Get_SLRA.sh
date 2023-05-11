#!/bin/zsh
# description: code for getting SLRA data from ERA5LAND, since the cldas longwave thermal radiation data is not available now
# it is also calculated by an empirical expression, which is not accurate
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com
# remove the offset and reduce the size of the data by converting the data to float32

# PLEASE CHANGE THE PATH TO YOUR OWN PATH
INPATH=/home/zhongwang/CLDAS_NRT_ASI_0P0625_HOR/SSRA
OUTPATH=/home/zhongwang/CLDAS_NRT_ASI_0P0625_HOR/SLRA
SYear=2008
EYear=2020
Months="01 02 03 04 05 06 07 08 09 10 11 12"
    while [ ${SYear} -le ${EYear} ] ; do
        Year=${SYear}
        for Month in ${Months};do
		ln -sf SSRA/CLDAS_NRT_ASI_0P0625_HOR-SSRA-${Year}${Month}.nc ${INPATH}/sample.nc
		cdo remapbil,sample.nc ${INPATH}/ERA5LAND/surface_thermal_radiation_downwards_w_m2/ERA5LAND_${Year}_${Month}_surface_thermal_radiation_downwards_w_m2.nc ${OUTPATH}/tstr/CLDAS_NRT_ASI_0P0625_HOR-tstr-${Year}${Month}.nc
        done
        SYear=`expr $SYear + 1`
    done




