#!/bin/zsh

# code for remove the dimension of z to reduce the size of data and make it easy to read in colm
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

SYear=1948
EYear=2006
Months="01 02 03 04 05 06 07 08 09 10 11 12"
varnames="dswrf dlwrf prcp pres shum tas wind"
while [ ${SYear} -le ${EYear} ] ; do
    for varname in ${varnames}; do
        Year=${SYear}
        echo $varname ${Year}
        cdo --reduce_dim copy /tera06/zhwei/CoLM_Forcing/princeton/${varname}/${varname}_3hourly_${Year}-${Year}.nc ${varname}/${varname}_3hourly_${Year}-${Year}.nc	
    done
    SYear=`expr $SYear + 1`
done
