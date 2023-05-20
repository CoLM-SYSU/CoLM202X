#!/bin/zsh
# description: code for remove spikes, since the cldas SHU data has lots of spikes.
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

# PLEASE CHANGE THE PATH TO YOUR OWN PATH
INPATH=./
OUTPATH=./

SYear=2008
EYear=2020
Months="01 02 03 04 05 06 07 08 09 10 11 12"
    while [ ${SYear} -le ${EYear} ] ; do
        Year=${SYear}
        for Month in ${Months};do
		cdo setrtoc,-100.0,0.0,0.001 ${INPATH}/SHU-bk/CLDAS_NRT_ASI_0P0625_HOR-SHU-${SYear}${Month}.nc  ${OUTPATH}/SHU/CLDAS_NRT_ASI_0P0625_HOR-SHU-${SYear}${Month}.nc
        done
        SYear=`expr $SYear + 1`
    done




