# Makefile for CLM main program

####################################################################

include ../include/Makeoptions
 
INCLUDE_DIR = -I../include -I${NETCDF_INC}

HEADER = ../include/define.h

OBJS =                          \
	precision.o             \
	GlobalVars.o            \
	LC_Const.o              \
	MathConstants.o         \
	PhysicalConstants.o     \
	spmd_task.o             \
	mod_namelist.o          \
	mod_utils.o             \
	timemanager.o           \
	ncio_serial.o           \
	mod_block.o             \
	mod_grid.o              \
	mod_pixel.o             \
	mod_data_type.o         \
	ncio_block.o            \
	mod_catchment_data.o    \
	mod_5x5_data.o          \
	mod_mesh.o              \
	mod_pixelset.o          \
	ncio_vector.o           \
	mod_colm_debug.o        \
	mod_mapping_pset2grid.o \
	mod_mapping_grid2pset.o


####################################################################

main : ${OBJS} ${HEADER}
	@echo 'Compiling modules in share directory completed!'

$(OBJS) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

####################################################################

clean :
	${RM} -f *.o *.mod
