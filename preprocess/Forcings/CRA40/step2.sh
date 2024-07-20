#!/bin/zsh
# code for merging CRA40 data 
# prepared by: zhongwang Wei @ SYSU 2024-04-10, Zhongwang007@gmail.com


#varnames="spfh"
INDIR=/Users/wei/Desktop/CRA40-old/19800102
OUTDIR=/Users/wei/Desktop/

SYear=1980
EYear=1980

for file in ${INDIR}/*.grib2; do
    timestamp=$(echo $file | cut -d'_' -f3)  # Extract timestamp like 1980010106
    year=${timestamp:0:4}
    month=${timestamp:4:2}
    day=${timestamp:6:2}
    hour=${timestamp:8:2}
    filename="${file##*/}"  # Get just the filename
    filename="${filename%.*}"

    echo $year $month $day $hour $filename


    ncl_convert2nc ${file} #${OUTDIR}${filename}.nc
    #use cdo to set time axis
    cdo settaxis,${year}-${month}-${day},${hour}:00:00,6hour ${filename}.nc ${OUTDIR}${filename}.nc
done


for file in ${INDIR}/*.grib; do
    timestamp=$(echo $file | cut -d'_' -f3)  # Extract timestamp like 1980010106
    year=${timestamp:0:4}
    month=${timestamp:4:2}
    day=${timestamp:6:2}
    hour=${timestamp:8:2}
    filename="${file##*/}"  # Get just the filename
    filename="${filename%.*}"

    echo $year $month $day $hour $filename

    ncl_convert2nc ${file} #${OUTDIR}${filename}.nc
    #use cdo to set time axis
    cdo settaxis,${year}-${month}-${day},${hour}:00:00,6hour ${filename}.nc ${OUTDIR}${filename}.nc

done

# -chname,lat_0,lat -chname,lon_0,lon 

while [ ${SYear} -le ${EYear} ] ; do
    Year=${SYear}
    cdo -b F32 -f nc4c -selname,DLWRF_P8_L1_GLL0_avg,DSWRF_P8_L1_GLL0_avg,PRATE_P8_L1_GLL0_avg -mergetime ${OUTDIR}/CRA40_SINGLEF_${Year}*_GLB_0P25_HOUR_V1_1_2.nc ${OUTDIR}/CRA40_Radiation_precip_${Year}.nc
    cdo -b F32 -f nc4c  -mergetime ${OUTDIR}/CRA40LAND_SURFACE_${Year}*_GLB_0P25_HOUR_V1_0_0.nc ${OUTDIR}/CRA40_sp_t_u_v_${Year}.nc
    cdo -b F32 -f nc4c -selname,PRES_P0_L1_GLL0 -mergetime ${OUTDIR}/CRA40_SINGLE_${Year}*_GLB_0P25_HOUR_V1_0_0.nc ${OUTDIR}/CRA40_pres_${Year}.nc
    cdo remapbil,${OUTDIR}/CRA40_sp_t_u_v_${Year}.nc ${OUTDIR}/CRA40_pres_${Year}.nc ${OUTDIR}/CRA40_pres_remap_${Year}.nc
    cdo remapbil,${OUTDIR}/CRA40_sp_t_u_v_${Year}.nc  ${OUTDIR}/CRA40_Radiation_precip_${Year}.nc  ${OUTDIR}/CRA40_Radiation_precip_remap_${Year}.nc

    SYear=`expr $SYear + 1`
done
