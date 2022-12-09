# Makefile for post processing

include ../include/Makeoptions

INCLUDE_DIR = -I../include -I../share -I../main -I${NETCDF_INC}

VPATH  = ../share
HEADER = ../include/define.h

DEF = $(shell grep -i CATCHMENT ../include/define.h)
vector2grid = $(word 1, ${DEF})

ifneq (${vector2grid},\#define)
DEF = $(shell grep -i UNSTRUCTURED ../include/define.h)
vector2grid = $(word 1, ${DEF})
endif

# Name of executable
EXECUTABLE = hist_concatenate.x srfdata_concatenate.x post_vector2grid.x

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

OBJS1 = mod_concatenate.o hist_concatenate.o
OBJS2 = mod_vector2grid.o post_vector2grid.o
OBJS3 = srfdata_concatenate.o

####################################################################

ifneq (${vector2grid},\#define)
all : hist_concatenate.x srfdata_concatenate.x
else
all : hist_concatenate.x srfdata_concatenate.x post_vector2grid.x
endif

hist_concatenate.x : ${HEADER} ${OBJS_SHARED} ${OBJS1}
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS1} -o $@ ${LDFLAGS}

post_vector2grid.x : ${HEADER} ${OBJS_SHARED} ${OBJS2}
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS2} -o $@ ${LDFLAGS}

srfdata_concatenate.x : ${HEADER} ${OBJS_SHARED} ${OBJS3}
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS3} -o $@ ${LDFLAGS}

$(OBJS1):%.o:%.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

$(OBJS2):%.o:%.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

$(OBJS3):%.o:%.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

$(OBJS_SHARED) : %.o : %.F90 ${HEADER} 
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../share
####################################################################

clean :
	${RM} -f *.o *.mod ${EXECUTABLE} 
