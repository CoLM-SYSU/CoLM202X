# Makefile for post processing

include ../include/Makeoptions

# name of executable
EXE1 = bin2netcdf

####################################################################
.DEFAULT :

OBJS =  precision.o \
	GlobalVars.o \
	PhysicalConstants.o \
	MOD_2D_Fluxes.o \
	bin2netcdf.o

####################################################################

${EXE1} : ${OBJS} 
	${FF} ${FOPTS} ${OBJS} -o $@ ${LDFLAGS} 
	@echo 'making bin2netcdf completed!'
	@echo ''
	@echo '!!NOTE: make sure to define lat/lon points in '
	@echo 'the beginning of bin2netcdf.F90 as you expect.'

$(OBJS):%.o:%.F90
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

clean :
	${RM} -f *.o *.mod 
