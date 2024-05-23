#!/bin/bash
# code for merging dlwrf data from 1979 to 2013
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

SYear=2003
EYear=2023
#varnames="spfh"
Months="01 02 03 04 05 06 07 08 09 10 11 12"

while [ ${SYear} -le ${EYear} ] ; do
    Year=${SYear}
    for Month in ${Months};do
        cdo  --reduce_dim -b F32 -f nc selname,tprate1have-sfc-fc-gauss tprate/jra3q.fcst_phy2m.0_1_52.tprate1have-sfc-fc-gauss.${Year}${Month}*.nc tprate/tprate_${Year}_${Month}.nc
    done
    SYear=`expr $SYear + 1`
done
