# Makefile for post processing

include ../include/Makeoptions

# Name of executable
EXECUTABLE = hist_concatenate.x

####################################################################
.DEFAULT :

HEADER = ../include/define.h

OBJS_SHARED = \
				  ../share/precision.o \
				  ../share/MathConstants.o \
				  ../share/spmd_task.o \
				  ../share/mod_utils.o \
				  ../share/GlobalVars.o \
				  ../share/mod_namelist.o \
				  ../share/ncio_serial.o \
				  ../share/mod_block.o \
				  ../share/mod_grid.o 

OBJS = mod_hist_nc.o \
		 hist_concatenate.o

####################################################################

# Compile CoLM decoupled with river routing scheme (CaMa-Flood)

${EXECUTABLE} : ${HEADER} ${OBJS_SHARED} ${OBJS}
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS} -o $@ ${LDFLAGS}

$(OBJS):%.o:%.F90
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

$(OBJS_SHARED) : %.o : %.F90 ${HEADER} 
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../share

####################################################################

clean :
	${RM} -f *.o *.mod ${EXECUTABLE} 
