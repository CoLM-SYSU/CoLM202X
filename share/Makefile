# Makefile for CoLM main program

####################################################################

include ../include/Makeoptions
 
INCLUDE_DIR = -I../include -I${NETCDF_INC}

HEADER = ../include/define.h

OBJS =                       \
	MOD_Precision.o           \
	GlobalVars.o              \
	LC_Const.o                \
	MathConstants.o           \
	PhysicalConstants.o       \
	MOD_SPMD_Task.o           \
	MOD_Namelist.o            \
	MOD_Utils.o               \
	MOD_TimeManager.o         \
	MOD_NetCDFSerial.o        \
	MOD_Block.o               \
	MOD_Grid.o                \
	MOD_Pixel.o               \
	MOD_DataType.o            \
	MOD_NetCDFBlock.o         \
	MOD_CatchmentDataReadin.o \
	MOD_5x5DataReadin.o       \
	MOD_Mesh.o                \
	MOD_Pixelset.o            \
	MOD_NetCDFVector.o        \
	MOD_CoLMDebug.o           \
	MOD_Mapping_Pset2Grid.o   \
	MOD_Mapping_Grid2Pset.o


####################################################################

main : ${OBJS} ${HEADER}
	@echo 'Compiling modules in share directory completed!'

$(OBJS) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

####################################################################

clean :
	${RM} -f *.o *.mod
