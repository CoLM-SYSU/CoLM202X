#!/bin/zsh
# code for merge the data
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com
# remove the offset and reduce the size of the data by converting the data to float32

# PLEASE CHANGE THE PATH TO YOUR OWN PATH
INPATH=/Users/wei/Desktop/TPMFD
OUTPATH=/Volumes/SourceData_1/TPMFD
SYear=1979
EYear=2020
Month1="01 02 03 04 05 06 07 08 09 10 11 12"
Vars="shum srad temp wind lrad prcp pres"
Vars="srad temp wind lrad prcp pres"
for var in ${Vars};do
mkdir $OUTPATH/${var}
SYear=1979
EYear=2020
while [ ${SYear} -le ${EYear} ] ; do
for Month in ${Month1};do
   cdo -f nc4c -b F32 -z zip_1 mergetime $INPATH/${var}/hourly/${SYear}/tpmfd_${var}_h_${SYear}${Month}*_00_23.nc    $OUTPATH/${var}/${var}_${SYear}${Month}.nc
done
SYear=`expr $SYear + 1`
done
done
