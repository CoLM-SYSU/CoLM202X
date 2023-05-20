#!/bin/zsh
# code for offseting time axis of ERA5LAND surface_solar_radiation_downwards data (1hour), and convert it to hourly areverage data from accumulation data.
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

ifort times.f newtime.f -o newtime
DIR=./scratch
source=../../original/surface_solar_radiation_downwards/
SYear=1969
EYear=1969
varnames="surface_solar_radiation_downwards"
Months="09 10 11"

for varname in ${varnames}; do
    mkdir -p ${DIR}/${varname}
    while [ ${SYear} -le ${EYear} ] ; do
        for Month in ${Months};do
            DAY=1
            EDDAY=`cal ${Month} ${SYear} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
            echo "Last Day in ${SYear} ${Month} is $EDDAY"
            Cdate=${SYear}${Month}01_00
            date_temp=`./newtime ${Cdate} -1`
            echo ${date_temp}
            pyr=`echo ${date_temp:0:4}`
            pmn=`echo ${date_temp:4:2}`
            pdy=`echo ${date_temp:6:2}`
            phr=`echo ${date_temp:9:2}`

            Cdate=${SYear}${Month}${EDDAY}_23
            date_temp=`./newtime ${Cdate} 1`
            echo ${date_temp}
            nyr=`echo ${date_temp:0:4}`
            nmn=`echo ${date_temp:4:2}`
            ndy=`echo ${date_temp:6:2}`
            nhr=`echo ${date_temp:9:2}`
            #echo $nyr $nmn $ndy $nhr
            cdo  -b F32 seldate,${nyr}-${nmn}-${ndy}T${nhr}:00:00 $source/ERA5LAND_${nyr}_${nmn}_${varname}.nc4 ${DIR}/${varname}/NEXT.nc
            cdo  -b F32 seldate,${pyr}-${pmn}-${pdy}T${phr}:00:00 $source/ERA5LAND_${pyr}_${pmn}_${varname}.nc4 ${DIR}/${varname}/PREV.nc
            cdo  -b F32 mergetime ${DIR}/${varname}/PREV.nc ${DIR}/${varname}/NEXT.nc $source/ERA5LAND_${SYear}_${Month}_${varname}.nc4 ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_withNEXTPREV.nc
            Hours="00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23"
            #Hours="00 23"
            for Hour in ${Hours}; do
                cdo -b F32 selhour,${Hour} ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_withNEXTPREV.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_${Hour}_temp.nc
            done                
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_00_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_23_temp.nc ${DIR}/${varname}/temp_00.nc
            cdo -b F32 copy ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_01_temp.nc  ${DIR}/${varname}/temp_01.nc
            cdo delete,timestep=1  ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_23_temp.nc  ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_23_tempr.nc
            rm ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_23_temp.nc 
            mv ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_23_tempr.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_23_temp.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_02_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_01_temp.nc ${DIR}/${varname}/temp_02.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_03_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_02_temp.nc ${DIR}/${varname}/temp_03.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_04_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_03_temp.nc ${DIR}/${varname}/temp_04.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_05_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_04_temp.nc ${DIR}/${varname}/temp_05.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_06_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_05_temp.nc ${DIR}/${varname}/temp_06.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_07_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_06_temp.nc ${DIR}/${varname}/temp_07.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_08_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_07_temp.nc ${DIR}/${varname}/temp_08.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_09_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_08_temp.nc ${DIR}/${varname}/temp_09.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_10_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_09_temp.nc ${DIR}/${varname}/temp_10.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_11_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_10_temp.nc ${DIR}/${varname}/temp_11.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_12_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_11_temp.nc ${DIR}/${varname}/temp_12.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_13_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_12_temp.nc ${DIR}/${varname}/temp_13.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_14_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_13_temp.nc ${DIR}/${varname}/temp_14.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_15_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_14_temp.nc ${DIR}/${varname}/temp_15.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_16_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_15_temp.nc ${DIR}/${varname}/temp_16.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_17_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_16_temp.nc ${DIR}/${varname}/temp_17.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_18_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_17_temp.nc ${DIR}/${varname}/temp_18.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_19_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_18_temp.nc ${DIR}/${varname}/temp_19.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_20_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_19_temp.nc ${DIR}/${varname}/temp_20.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_21_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_20_temp.nc ${DIR}/${varname}/temp_21.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_22_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_21_temp.nc ${DIR}/${varname}/temp_22.nc
            cdo -b F32 sub ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_23_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_22_temp.nc ${DIR}/${varname}/temp_23.nc

            cdo mergetime ${DIR}/${varname}/temp_*.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_w_m2_temp.nc
            cdo delete,timestep=1 ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_w_m2_temp.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_w_m2_temp1.nc
            #cdo setattribute,tp@units="m/hr" ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_m_hr_temp1.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_m_hr_temp2.nc
            #ncap2 -s 'where(tp<0.) tp=0.;' ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_m_hr_temp2.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_m_hr_temp3.nc
            cdo setattribute,ssrd@units="w/m2" ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_w_m2_temp1.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_w_m2_temp2.nc
            cdo expr,'ssrd=ssrd/3600' ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_w_m2_temp2.nc ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_w_m2_temp3.nc
            #cdo -f nc4c -z zip_6 settaxis,${SYear}-${Month}-01,00:30:00,1hour ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_m_hr_temp3.nc  ERA5LAND_${SYear}_${Month}_${varname}_m_hr.nc
            #rm ${DIR}/${varname}/*.nc            
            cdo -f nc4c -z zip_6 settaxis,${SYear}-${Month}-01,00:30:00,1hour ${DIR}/${varname}/ERA5LAND_${SYear}_${Month}_${varname}_w_m2_temp3.nc  ERA5LAND_${SYear}_${Month}_${varname}_w_m2.nc
            rm ${DIR}/${varname}/*.nc

        done
        SYear=`expr $SYear + 1`
    done
done




