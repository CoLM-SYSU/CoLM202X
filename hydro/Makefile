# Makefile for CLM surface data processing

include ../include/Makeoptions

####################################################################
.DEFAULT :

OBJS_SHARED =  ../share/precision.o \
					../share/spmd_task.o \
					../share/hdf5_io_serial.o \
					../share/mod_block.o \
					../share/mod_grid.o  \
					../share/mod_data_type.o

OBJS = mod_hydro_data.o

####################################################################

main : ${OBJS}
	@echo 'Compiling modules in hydro completed!'

$(OBJS):%.o : %.F90 ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

clean :
	${RM} -f *.o *.mod
