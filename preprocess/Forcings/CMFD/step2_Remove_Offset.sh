#!/bin/zsh
# code for remove offset and reduce the size of the data
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com
# remove the offset and reduce the size of the data by converting the data to float32

# PLEASE CHANGE THE PATH TO YOUR OWN PATH
SYear=1979
EYear=2022
varnames="prec pres shum srad temp wind lrad"
Months="01 02 03 04 05 06 07 08 09 10 11 12"
INPATH=./
OUTPATH=./
for varname in ${varnames}; do
    while [ ${SYear} -le ${EYear} ] ; do
        Year=${SYear}
        for Month in ${Months};do
           cdo -f nc4c -b F32 copy ${INPATH}/${varname}/${varname}_CMFD_V0106_B-01_03hr_010deg_${Year}${Month}.nc ${OUTPATH}/${varname}/${varname}_CMFD_V0106_B-01_03hr_010deg_${Year}${Month}.nc4
	done 
        SYear=`expr $SYear + 1`
     done
done




