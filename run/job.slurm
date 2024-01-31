#!/bin/bash

#SBATCH -N 1 -n 40

export LD_LIBRARY_PATH=/opt/netcdf-c-4.9.2-fortran-4.6.0-gnu/lib:/opt/hdf5-1.14.0-gnu/lib:$LD_LIBRARY_PATH

# mpirun -np 40 run/mksrfdata.x run/GreaterBay_Grid_10km_IGBP_VG.nml
mpirun -np 40 run/mkinidata.x run/GreaterBay_Grid_10km_IGBP_VG.nml
mpirun -np 40 run/colm.x      run/GreaterBay_Grid_10km_IGBP_VG.nml
