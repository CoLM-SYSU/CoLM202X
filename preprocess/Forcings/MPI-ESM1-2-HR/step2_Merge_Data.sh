
#!/bin/zsh

# code for merge the dataset into one file
# prepared by: zhongwang Wei @ SYSU 2021-10-20, Zhongwang007@gmail.com


INPath=./
OUTPath=./
SSPs="585"
models="MPI-ESM1-2-HR"
Vars="huss  pr  ps  rlds  rsds  rsus  tas  uas  vas"
Months="01 02 03 04 05 06 07 08 09 10 11 12"

for Var in ${Vars};do
for model in ${models};do
for SSP in ${SSPs};do
cdo -b F32 -f nc mergetime ${INPath}/${Var}/${Var}_3hr_${model}_ssp${SSP}_r2i1p1f1_gn_*.nc ${INPath}/${Var}/${Var}_3hr_${model}_ssp${SSP}_r2i1p1f1_gn.nc
done
done
done



