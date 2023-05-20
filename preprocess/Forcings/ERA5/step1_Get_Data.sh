#!/bin/zsh
# code for downloading ERA5LAND data
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com
# Please set up cds api key first:
# --->https://cds.climate.copernicus.eu/api-how-to
# --->https://cds.climate.copernicus.eu/api-how-to/installation-of-the-client

# PLEASE CHANGE THE PATH TO YOUR OWN PATH
DIR=./
SYear=2000
EYear=2021
#varnames="surface_solar_radiation_downwards"
varnames="volumetric_soil_water_layer_2"
Month1="01 02 03 04 05 06 07 08 09 10 11 12"
while [ ${SYear} -le ${EYear} ] ; do
mkdir -p ${varnames}
for Month in ${Month1};do
for varname in ${varnames}; do
if [ ! -s ${DIR}/${varnames}/ERA5_${SYear}_${Month}_${varname}.nc ] ; then
cat <<EOF >ERA.py
#!/usr/bin/env python
import cdsapi
c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
         'variable':[
           '${varname}'
        ],
        'year':'${SYear}',
        'month':'${Month}',
        'day':[
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31','01','02',
            '03','04','05',
            '06','07','08',
            '09'
        ],
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ],
        'format':'netcdf'
    },
    '${DIR}/${varnames}/ERA5_${SYear}_${Month}_${varname}.nc')
EOF
python ERA.py
else
  echo "ERA5_${SYear}_${Month}_${varname}.nc exist"
fi
#mv 100m_uv_component_of_wind_${Year}_${Month}.nc ${Year}
done
done
SYear=`expr $SYear + 1`
done
