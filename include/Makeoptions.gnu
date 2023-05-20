# =======================================================
# mpif90 - ifort 
# 


 FF = /opt/mpich-3.4.2-intel/bin/mpif90

 NETCDF_LIB = /home/zhwei/software/NETCDF/c-4.7.3-f4.5.2/lib
 NETCDF_INC = /home/zhwei/software/NETCDF/c-4.7.3-f4.5.2/include

MOD_CMD = -module 

#<<<<<<< HEAD
#FOPTS = -qopenmp -g -traceback -r8 -free -check uninit 
       # -r8 -free -O0 -check uninit -check bounds -check pointers \
       # -traceback  -assume byterecl -pthread -heap-arrays #-nogen-interface

#INCLUDE_DIR = -I../include -I../share -I../mksrfdata \
#               -I../mkinidata -I../main -I../hydro -I${NETCDF_INC}
#LDFLAGS = -L${NETCDF_LIB} -lnetcdff


 FOPTS = -g -qopenmp -traceback -r8 #-free -check uninit #-check bounds

 LDFLAGS = -L${NETCDF_LIB} -lnetcdff -llapack -lblas

# =======================================================
# mpif90 - gfortran 
# 

#  FF = /opt/mpich/bin/mpif90
#   
#  MPI_LIB = /opt/mpich/lib
#  NETCDF_LIB = /usr/lib/x86_64-linux-gnu
#  NETCDF_INC = /usr/include
#   
#  MOD_CMD = -J
# 
#  FOPTS = -fdefault-real-8 -ffree-form -C -g -u -xcheck=stkovf \
#           -ffpe-trap=invalid,zero,overflow -fbounds-check \
#           -mcmodel=medium -fbacktrace -fdump-core -cpp
#  
#  INCLUDE_DIR = -I../include -I../share -I../mksrfdata -I../mkinidata -I../main -I$(NETCDF_INC)
#  LDFLAGS = -L${MPI_LIB} -L$(NETCDF_LIB) -lnetcdff -lnetcdf
