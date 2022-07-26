# Makefile for preprocess programs

include ../include/Makeoptions

all : rawdata_to_nc
.PHONY : all

OBJS_SHARED =  ../share/precision.o \
					../io/spmd_task.o  \
					../io/ncio_serial.o 

OBJS = rawdata_to_nc.o

rawdata_to_nc : $(OBJS) ${OBJS_SHARED}
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS} -o $@ ${LDFLAGS}

$(OBJS) : %.o : %.F90 ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

$(OBJS_SHARED) : %.o : %.F90
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

.PHONY : clean
clean :
	${RM} -f *.o *.mod
