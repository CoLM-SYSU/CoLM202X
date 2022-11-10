# Makefile for CLM surface data processing

include ../../include/Makeoptions
 
INCLUDE_DIR = -I../../include -I../../share 

VPATH  = ../../share 
HEADER = ../../include/define.h

####################################################################
.DEFAULT :

OBJS_SHARED = ../../share/precision.o ../../share/PhysicalConstants.o

OBJS = mod_soil_function.o mod_soil_water.o

####################################################################

main : ${OBJS}
	@echo 'Compiling modules in hydro completed!'

$(OBJS):%.o : %.F90 ${HEADER} ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

$(OBJS_SHARED) : %.o : %.F90 ${HEADER} 
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../../share

clean :
	${RM} -f *.o *.mod
