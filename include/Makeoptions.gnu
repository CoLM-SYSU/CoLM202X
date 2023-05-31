# =======================================================
# mpif90 - gfortran 
# 

  FF = /usr/bin/mpif90 -fopenmp
   
  NETCDF_LIB = /usr/lib/x86_64-linux-gnu
  NETCDF_INC = /usr/include
   
  MOD_CMD = -J
 
  FOPTS = -fdefault-real-8 -ffree-form -C -g -u -xcheck=stkovf \
           -ffpe-trap=invalid,zero,overflow -fbounds-check \
           -mcmodel=medium -fbacktrace -fdump-core -cpp \
           -ffree-line-length-0  
  
  INCLUDE_DIR = -I../include -I../share -I../mksrfdata -I../mkinidata -I../main -I$(NETCDF_INC)
  LDFLAGS = -L$(NETCDF_LIB) -lnetcdff -lnetcdf -llapack -lblas 



#============================================================
# CaMa-Flood Mkinclude (for Linux, gfortran)

RM = /bin/rm -f
CP = /bin/cp
#----
# Pre-Prosessing options
# DMPI=-DUseMPI: activate when MPI parallelization is used
# DCDF=-DUseCDF: activate when using netCDF, comment out when not needed
# DATM=-DNoAtom: activate when OMP ATOMIC calculation should be avoided (bit identical simulation)
#----
#DMPI=-DUseMPI
DCDF=-DUseCDF
#DATM=-DNoAtom
CFLAGS=$(DMPI) $(DCDF) $(DATM) 
#----
FCMP = /usr/bin/gfortran -fopenmp
FC   = /usr/bin/gfortran 

LFLAGS =
FFLAGS = -O3 -Wall -cpp -free -fimplicit-none -fbounds-check -fbacktrace 
