
#!/bin/zsh

# code for remap the data to 1d grid
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
 cdo splityear ${INPath}/${Var}/${Var}_3hr_${model}_ssp${SSP}_r2i1p1f1_gn_1d.nc ${INPath}/${Var}/
done
done
done



