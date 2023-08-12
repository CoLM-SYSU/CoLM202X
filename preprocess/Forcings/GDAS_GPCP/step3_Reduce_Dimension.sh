#!/bin/zsh

# code for remove the dimension of z to reduce the size of data and make it easy to read in colm
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

#!/bin/bash
SYear=2002
EYear=2022
Months="01 02 03 04 05 06 07 08 09 10 11 12"
while [ ${SYear} -le ${EYear} ] ; do
    Year=${SYear}
    for Month in ${Months};do
        cdo --reduce_dim -f nc4c copy  GDAS_GPCP/GLDAS_GDAS_3H_tot_prcip.${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_tot_prcip.${Year}${Month}.nc4
        cdo --reduce_dim -f nc4c copy  GDAS_GPCP/GLDAS_GDAS_3H_Psurf.${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_Psurf.${Year}${Month}.nc4
        cdo --reduce_dim -f nc4c copy GDAS_GPCP/GLDAS_GDAS_3H_Qair.${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_Qair.${Year}${Month}.nc4
        cdo --reduce_dim -f nc4c copy GDAS_GPCP/GLDAS_GDAS_3H_LWdown.${Year}${Month}.nc  GDAS_GPCP/GLDAS_GDAS_3H_LWdown.${Year}${Month}.nc4
        cdo --reduce_dim -f nc4c copy GDAS_GPCP/GLDAS_GDAS_3H_SWdown.${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_SWdown.${Year}${Month}.nc4
        cdo --reduce_dim -f nc4c copy GDAS_GPCP/GLDAS_GDAS_3H_Tair.${Year}${Month}.nc GDAS_GPCP/GLDAS_GDAS_3H_Tair.${Year}${Month}.nc4
        cdo --reduce_dim -f nc4c copy  GDAS_GPCP/GLDAS_GDAS_3H_Wind.${Year}${Month}.nc  GDAS_GPCP/GLDAS_GDAS_3H_Wind.${Year}${Month}.nc4
    done 
    SYear=`expr $SYear + 1`
done
