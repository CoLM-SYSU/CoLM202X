#!/bin/zsh

# code for unzip WFDE5 data
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

DIR=/Volumes/work/Dataset/LSM_Forcings/crujra/

varnames="dlwrf dswrf pre pres spfh tmax tmin tmp ugrd vgrd"
for varname in ${varnames}; do
SYear=1901
EYear=2020
while [ ${SYear} -le ${EYear} ] ; do
gunzip ${DIR}/${varname}/crujra.v2.3.5d.${varname}.${SYear}.365d.noc.nc.gz 
SYear=`expr $SYear + 1`
done

done
