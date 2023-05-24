# Makefile for CoLM surface data processing

include ../include/Makeoptions

INCLUDE_DIR = -I../include -I../share -I${NETCDF_INC}

VPATH = ../share

# name of executable
EXECUTABLE = mksrfdata.x

####################################################################
.DEFAULT :

HEADER = ../include/define.h

<<<<<<< HEAD
OBJS_SHARED =                         \
	../share/MOD_Precision.o           \
	../share/MathConstants.o           \
	../share/MOD_SPMD_Task.o           \
	../share/GlobalVars.o              \
	../share/LC_Const.o                \
	../share/UrbanLCZ_Const.o          \
	../share/MOD_Namelist.o            \
	../share/MOD_TimeManager.o         \
	../share/MOD_Utils.o               \
	../share/MOD_NetCDFSerial.o        \
	../share/MOD_Block.o               \
	../share/MOD_Grid.o                \
	../share/MOD_Pixel.o               \
	../share/MOD_DataType.o            \
	../share/MOD_NetCDFBlock.o         \
	../share/MOD_CatchmentDataReadin.o \
	../share/MOD_5x5DataReadin.o       \
	../share/MOD_Mesh.o                \
	../share/MOD_Pixelset.o            \
	../share/MOD_NetCDFVector.o        \
	../share/MOD_CoLMDebug.o           \
	../share/MOD_Mapping_Pset2Grid.o
=======
OBJS_SHARED =                              \
	  ../share/precision.o             \
	  ../main/MOD_Vars_Global.o        \
	  ../main/MOD_Const_LC.o           \
	  ../share/MathConstants.o         \
	  ../share/spmd_task.o             \
	  ../share/mod_namelist.o          \
	  ../share/timemanager.o           \
	  ../share/mod_utils.o             \
	  ../share/ncio_serial.o           \
	  ../share/mod_block.o             \
	  ../share/mod_grid.o              \
	  ../share/mod_pixel.o             \
	  ../share/mod_data_type.o         \
	  ../share/ncio_block.o            \
	  ../share/mod_catchment_data.o    \
	  ../share/mod_5x5_data.o          \
	  ../share/mod_mesh.o              \
	  ../share/mod_pixelset.o          \
	  ../share/ncio_vector.o           \
	  ../share/mod_colm_debug.o        \
	  ../share/mod_mapping_pset2grid.o
>>>>>>> 3e1016d9a4bd5972c19f496b221102686701cc1d

OBJS =                                     \
	 MOD_AggregationRequestData.o      \
	 MOD_SingleSrfdata.o               \
	 MOD_SrfdataRestart.o              \
	 MOD_PixelsetShadow.o              \
	 MOD_LandElm.o                     \
	 MOD_LandHRU.o                     \
	 MOD_LandPatch.o                   \
	 MOD_LandUrban.o                   \
	 MOD_LandPFT.o                     \
	 MOD_LandPC.o                      \
	 MOD_MeshFilter.o                  \
	 MOD_SrfdataDiag.o                 \
	 Aggregation_PercentagesPFT.o      \
	 Aggregation_LAI.o                 \
	 Aggregation_NDeposition.o         \
	 Aggregation_Fire.o                \
	 Aggregation_CropParameters.o      \
	 Aggregation_NitrifParameters.o    \
	 Aggregation_SoilBrightness.o      \
	 Aggregation_LakeDepth.o           \
	 Aggregation_ForestHeight.o        \
	 Aggregation_SoilParameters.o      \
	 Aggregation_DBedrock.o            \
	 Aggregation_Topography.o          \
	 Aggregation_Urban.o               \
	 MOD_RegionClip.o                  \
	 MKSRFDATA.o

####################################################################

${EXECUTABLE} : ${HEADER} ${OBJS_SHARED} ${OBJS}
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS} -o $@ ${LDFLAGS}
	@echo 'making CoLM surface data completed!'

$(OBJS) : %.o : %.F90 ${HEADER} ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

$(OBJS_SHARED) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../share

clean :
	${RM} -f *.o *.mod ${EXECUTABLE}
