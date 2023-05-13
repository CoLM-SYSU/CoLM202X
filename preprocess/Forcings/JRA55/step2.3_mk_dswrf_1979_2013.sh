#!/bin/bash
# code for merging dswrf data from 1979 to 2013
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com
SYear=1979
EYear=2013
#varnames="spfh"
while [ ${SYear} -le ${EYear} ] ; do
        Year=${SYear}
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.204_dswrf.reg_tl319.${Year}010100_${Year}033121 ${Year}_1.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.204_dswrf.reg_tl319.${Year}040100_${Year}063021 ${Year}_2.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.204_dswrf.reg_tl319.${Year}070100_${Year}093021 ${Year}_3.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.204_dswrf.reg_tl319.${Year}100100_${Year}123121 ${Year}_4.nc
        cdo mergetime ${Year}_*.nc dswrf_${Year}.nc4
        cdo selname,var204 dswrf_${Year}.nc4 dswrf_${Year}.nc
    SYear=`expr $SYear + 1`
done
