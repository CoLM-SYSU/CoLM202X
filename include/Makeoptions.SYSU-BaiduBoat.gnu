# ==========================================================
# mpif90 - gnu 

#   please "source /share/home/dq089/soft/gnu-env" first.

FF = mpif90

NETCDF_LIB = /share/home/dq089/soft/netcdf-fortran-4.6.1-gnu/lib         
NETCDF_INC = /share/home/dq089/soft/netcdf-fortran-4.6.1-gnu/include     

MOD_CMD = -J

FOPTS = -fdefault-real-8 -ffree-form -C -g -u -xcheck=stkovf \
        -ffpe-trap=invalid,zero,overflow -fbounds-check \
        -mcmodel=medium -fbacktrace -fdump-core -cpp \
        -ffree-line-length-0 -fopenmp

INCLUDE_DIR = -I../include -I../share -I../mksrfdata -I../mkinidata -I../main -I$(NETCDF_INC)
LDFLAGS = -fopenmp -L${NETCDF_LIB} -lnetcdff -llapack -L/share/home/dq089/soft/lib -lblas.gnu

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
DSINGLE=-DSinglePrec_CMF

CFLAGS=$(DMPI) $(DCDF) $(DATM) $(DSINGLE)
#----
FCMP = gfortran -fopenmp
FC   = gfortran -fopenmp
FFLAGS = -O3 -Wall -cpp -ffree-line-length-none -fimplicit-none -ftree-vectorize 
