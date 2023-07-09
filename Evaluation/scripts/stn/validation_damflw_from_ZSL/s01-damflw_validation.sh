#!/bin/sh

#### initial setting ========================================
## project name
out_tag='src_v411_15min_h06'
map_tag='glb_15min'

## specify validation period
SYEAR=1979
SMON=1
SDAY=1
EYEAR=1981
EMON=12
EDAY=31

para_flag=1
num_cores=72

# ===========================================================
# Link input data in validate/ dir
# ===========================================================
PWDD=`pwd`
BASE="/tera02/zhangsl/cama-flood/cama_v4.10"

# obs data
OBSDIR="/tera02/zhangsl/cama-flood/obs/ResOpsUS/time_series_all"

# sim data
SIMDIR=${BASE}/out/${out_tag}

# dam data
DAMLIST=${BASE}/map/dam_para/dam_${map_tag}/dam_params_${map_tag}.csv

# save path
SAVEDIR=${PWDD}/validation/${out_tag}

echo "### DAMLIST = " $DAMLIST
echo "### SIMDIR = " $SIMDIR
echo "### SAVEDIR = " $SAVEDIR
# #===========================================================

mkdir -p ${SAVEDIR}/inp

cd /${SAVEDIR}


ln -fs $OBSDIR ./inp/obs

ln -fs $SIMDIR ./inp/sim

ln -fs $DAMLIST ./inp/damlist.csv


echo "#==========================================================="
echo "### dam validation"
echo "### checkpwdd" $PWDD
echo "#==========================================================="
python ${PWDD}/src/p01_damflw_validation.py $SYEAR $SMON $SDAY $EYEAR $EMON $EDAY $out_tag $para_flag $num_cores


exit 0