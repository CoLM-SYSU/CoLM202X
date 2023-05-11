#!/bin/zsh
# code for offset time axis of ERA5LAND specific_humidity data (1hour)
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

ifort times.f newtime.f -o newtime
DIR=./scratch
source=../../data/
Out=../../processing/
mkdir -p $DIR $Out
SYear=1989
EYear=1989
varnames="specific_humidity"
Months="06"

for varname in ${varnames}; do
    mkdir -p ${DIR}/${varname}
    while [ ${SYear} -le ${EYear} ] ; do
        for Month in ${Months};do
            rm temp.nc
            DAY=1
            EDDAY=`cal ${Month} ${SYear} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
            echo "Last Day in ${SYear} ${Month} is $EDDAY"
            
            Ndate=${SYear}${Month}${EDDAY}_23


            date_temp=`./newtime ${Ndate} +100`
            echo ${date_temp}
            nyr=`echo ${date_temp:0:4}`
            nmn=`echo ${date_temp:4:2}`
            ndy=`echo ${date_temp:6:2}`
            nhr=`echo ${date_temp:9:2}`

            cdo -b F32 -f nc4c mergetime ${source}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}.nc4 ${source}/${varname}/ERA5LAND_${nyr}_${nmn}_${varname}.nc4 temp.nc
            cdo shifttime,-1hour  temp.nc temp0.nc
            cdo selmonth,${Month} temp0.nc temp1.nc
            cdo -f nc4c -z zip_6 settaxis,${SYear}-${Month}-01,00:30:00,1hour temp1.nc $Out/${varname}/ERA5LAND_${SYear}_${Month}_${varname}.nc
            rm temp*.nc

        done
        SYear=`expr $SYear + 1`
    done
done
#rm -rf $DIR/${varname}
