# Makefile for post processing

include ../include/Makeoptions

INCLUDE_DIR = -I../include -I../share -I../main -I${NETCDF_INC}

VPATH  = ../share
HEADER = ../include/define.h

# Name of executable
EXECUTABLE = hist_concatenate.x

####################################################################
.DEFAULT :


OBJS_SHARED = \
				  ../share/precision.o \
				  ../share/MathConstants.o \
				  ../share/spmd_task.o \
				  ../share/mod_utils.o \
				  ../share/timemanager.o \
				  ../share/GlobalVars.o \
				  ../share/mod_namelist.o \
				  ../share/ncio_serial.o \
				  ../share/mod_block.o \
				  ../share/mod_grid.o 

OBJS = mod_concatenate.o hist_concatenate.o

####################################################################

# Compile CoLM decoupled with river routing scheme (CaMa-Flood)

${EXECUTABLE} : ${HEADER} ${OBJS_SHARED} ${OBJS_MAIN} ${OBJS}
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS} -o $@ ${LDFLAGS}

$(OBJS):%.o:%.F90 ${HEADER} ${OBJS_MAIN}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

$(OBJS_SHARED) : %.o : %.F90 ${HEADER} 
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../share
####################################################################

clean :
	${RM} -f *.o *.mod ${EXECUTABLE} 
