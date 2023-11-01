
#!/bin/zsh

# code for remove the dimension of z to reduce the size of data and make it easy to read in colm
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com

SYear=2002
EYear=2022
INPath=MSWX
OUTPath=MSWX
Vars="LWd P Pres RelHum SWd SpecHum Temp Tmax Tmin Wind"
Months="01 02 03 04 05 06 07 08 09 10 11 12"

for Var in ${Vars};do
cd ${INPath}/${Var}/3hourly/
echo "${INPath}/${Var}/3hourly/"
while [ ${SYear} -le ${EYear} ] ; do
    Year=${SYear}
    for file in ${Year}*.nc; do
        echo "${file}"
        year=$(echo $file | cut -c1-4)
        day=$(echo $file | cut -c5-7)
        day=$(expr $day + 0)
        hour=$(echo $file | cut -c9-10)
        date=$(gdate -d "$year-01-01 +$((day-1)) days" "+%Y%m%d")
        cp "$file" "${date}.${hour}.nc4"
    done
    SYear=`expr $SYear + 1`
done
cd ../../
done

for Var in ${Vars};do
cd ${INPath}/${Var}
while [ ${SYear} -le ${EYear} ] ; do
    Year=${SYear}
        for Month in ${Months};do
           cdo -b F32 -f nc mergetime 3hourly/${Year}${Month}*.nc4 var_temp.nc
           cdo  invertlat var_temp.nc ${Var}_${Year}${Month}.nc
        done 
    SYear=`expr $SYear + 1`
done
cd ${INPath}
done

