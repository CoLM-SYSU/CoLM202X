#!/bin/bash
# code for merging tprat data from 2014 to 2022  
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

SYear=2014
EYear=2022
#varnames="spfh"
while [ ${SYear} -le ${EYear} ] ; do
        Year=${SYear}
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}01*  ${Year}_1.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}02*  ${Year}_2.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}03*  ${Year}_3.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}04*  ${Year}_4.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}05*  ${Year}_5.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}06*  ${Year}_6.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}07*  ${Year}_7.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}08*  ${Year}_8.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}09*  ${Year}_9.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}10*  ${Year}_10.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}11*  ${Year}_11.nc
        cdo  --reduce_dim -b F32 -f nc copy fcst_phy2m.061_tprat.reg_tl319.${Year}12*  ${Year}_12.nc

        cdo mergetime ${Year}_*.nc tprat_${Year}.nc
        rm ${Year}_*.nc

    SYear=`expr $SYear + 1`
done
