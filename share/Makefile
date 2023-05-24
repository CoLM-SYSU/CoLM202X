# Makefile for CoLM main program

####################################################################

include ../include/Makeoptions
 
INCLUDE_DIR = -I../include -I${NETCDF_INC}

HEADER = ../include/define.h

<<<<<<< HEAD
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
=======
OBJS =                              \
	  precision.o               \
	  ../main/MOD_Vars_Global.o \
	  MathConstants.o           \
	  spmd_task.o               \
	  mod_namelist.o            \
	  mod_utils.o               \
	  timemanager.o             \
	  ncio_serial.o             \
	  mod_block.o               \
	  mod_grid.o                \
	  mod_pixel.o               \
	  mod_data_type.o           \
	  ncio_block.o              \
	  mod_catchment_data.o      \
	  mod_5x5_data.o            \
	  mod_mesh.o                \
	  mod_pixelset.o            \
	  ncio_vector.o             \
	  mod_colm_debug.o          \
	  mod_mapping_pset2grid.o   \
	  mod_mapping_grid2pset.o
>>>>>>> 3e1016d9a4bd5972c19f496b221102686701cc1d


####################################################################

main : ${OBJS} ${HEADER}
	@echo 'Compiling modules in share directory completed!'

$(OBJS) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

####################################################################

clean :
	${RM} -f *.o *.mod
