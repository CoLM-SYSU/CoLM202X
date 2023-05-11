#!/bin/zsh
SSYear=1979
EYear=2022
varnames="volumetric_soil_water_layer_4"
Month1="01 02 03 04 05 06 07 08 09 10 11 12"
for varname in ${varnames}; do
SYear=$SSYear
INDIR="./"
OUTDIR="./"
while [ ${SYear} -le ${EYear} ] ; do
for Month in ${Month1};do
zip $OUTDIR/ERA5LAND_${SYear}_${Month}_${varname}.nc.zip $INDIR/ERA5LAND_${SYear}_${Month}_${varname}.nc
#rm $INDIR/ERA5LAND_${SYear}_${Month}_${varname}.nc
done
SYear=`expr $SYear + 1`
done
done


