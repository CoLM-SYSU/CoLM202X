# =======================================================
# mpif90 - gfortran 
# 

  FF = /usr/bin/mpif90 
   
  NETCDF_LIB = /usr/lib/x86_64-linux-gnu
  NETCDF_INC = /usr/include
   
  MOD_CMD = -J
 
  FOPTS = -fdefault-real-8 -ffree-form -C -g -u -xcheck=stkovf \
           -ffpe-trap=invalid,zero,overflow -fbounds-check \
           -mcmodel=medium -fbacktrace -fdump-core -cpp \
           -ffree-line-length-0
  
  INCLUDE_DIR = -I../include -I../share -I../mksrfdata -I../mkinidata -I../main -I$(NETCDF_INC)
  LDFLAGS = -L$(NETCDF_LIB) -lnetcdff -lnetcdf
