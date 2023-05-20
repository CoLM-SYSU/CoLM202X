#!/bin/zsh
# code for remove offset and reduce the size of the data
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com
# remove the offset and reduce the size of the data by converting the data to float32
# PLEASE CHANGE THE PATH TO YOUR OWN PATH
DIR=./
SYear=2000
EYear=2021
varnames="boundary_layer_height"
Month1="01 02 03 04 05 06 07 08 09 10 11 12"
while [ ${SYear} -le ${EYear} ] ; do
for Month in ${Month1};do
for varname in ${varnames}; do
 cdo -f nc4c -b F32 copy ERA5_${SYear}_${Month}_${varname}.nc ERA5_${SYear}_${Month}_${varname}.nc4
done
done
SYear=`expr $SYear + 1`
done

