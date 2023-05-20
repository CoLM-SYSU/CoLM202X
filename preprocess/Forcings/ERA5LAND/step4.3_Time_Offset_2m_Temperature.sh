#!/bin/zsh
# code for offset time axis of ERA5LAND 2m temperature data (1hour)
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

ifort times.f newtime.f -o newtime
DIR=./scratch
source=/Volumes/原始数据/ERA5LAND/data/
Out=../../processing/
mkdir -p $DIR $Out
SYear=2000
EYear=2000
varnames="2m_temperature"
Months="12"

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
        
            #cdo -b F32 -f nc4c -shifttime,-1hour ${source}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}.nc4  temp.nc
            #cdo selmonth,${Month} temp.nc $DIR/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_1.nc4
            #cdo selmonth,${pmn}   temp.nc $DIR/${varname}/ERA5LAND_${pyr}_${pmn}_${varname}_last.nc4
            #cdo mergetime $DIR/${varname}/ERA5LAND_${pyr}_${pmn}_${varname}_*.nc4  ERA5LAND_${pyr}_${pmn}_${varname}_temp.nc
            #cdo -f nc4c -z zip_6 settaxis,${pyr}-${pmn}-01,00:30:00,1hour ERA5LAND_${pyr}_${pmn}_${varname}_temp.nc $Out/${varname}/ERA5LAND_${pyr}_${pmn}_${varname}.nc

            #cdo settaxis,${pyr}-${pmn}-01,00:30:00,1hour ERA5LAND_${pyr}_${pmn}_${varname}_temp.nc ERA5LAND_${pyr}_${pmn}_${varname}.nc
            #cdo setattribute,slhf@units="Kg/Kg" ERA5LAND_${pyr}_${pmn}_${varname}.nc ERA5LAND_${pyr}_${pmn}_${varname}_temp.nc
            #cdo -f nc4c -z zip_6 expr,'slhf=slhf/3600.0' ERA5LAND_${pyr}_${pmn}_${varname}_temp.nc ERA5LAND_${pyr}_${pmn}_${varname}_w_m2.nc
            #rm ERA5LAND_${pyr}_${pmn}_${varname}_temp.nc RA5LAND_${pyr}_${pmn}_${varname}.nc  ERA5LAND_${pyr}_${pmn}_${varname}_temp.nc
        done
        SYear=`expr $SYear + 1`
    done
done
#rm -rf $DIR/${varname}
