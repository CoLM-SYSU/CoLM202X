# Makefile for CLM main program

####################################################################

include ../include/Makeoptions

HEADER = ../include/define.h

OBJS =  \
	precision.o         \
	MathConstants.o     \
	PhysicalConstants.o \
	spmd_task.o         \
	mod_namelist.o      \
	ncio_serial.o       \
	mod_utils.o         \
	mod_block.o         \
	mod_grid.o          \
	mod_pixel.o         \
	mod_data_type.o     \
	ncio_block.o        \
	mod_hydro_data.o    \
	mod_landunit.o      \
	mod_pixelset.o      \
	ncio_vector.o       \
	timemanager.o       \
	mod_colm_debug.o     


####################################################################

main : ${OBJS} ${HEADER}
	@echo 'Compiling modules in share directory completed!'

$(OBJS) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

####################################################################

clean :
	${RM} -f *.o *.mod
