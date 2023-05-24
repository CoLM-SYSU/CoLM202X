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
	  mod_aggregation.o                \
	  mod_single_srfdata.o             \
	  mod_srfdata_restart.o            \
	  mod_pixelsetshadow.o             \
	  mod_landelm.o                    \
	  mod_landhru.o                    \
	  mod_landpatch.o                  \
	  mod_landurban.o                  \
	  mod_landpft.o                    \
	  mod_landpc.o                     \
	  mod_mesh_filter.o                \
	  mod_srfdata_diag.o               \
	  aggregation_percentages.o        \
	  aggregation_LAI.o                \
	  aggregation_NDEP.o               \
	  aggregation_fire.o               \
	  aggregation_crop_parameters.o    \
	  aggregation_nitrif_parameters.o  \
	  aggregation_soil_brightness.o    \
	  aggregation_lakedepth.o          \
	  aggregation_forest_height.o      \
	  aggregation_soil_parameters.o    \
	  aggregation_dbedrock.o           \
	  aggregation_topography.o         \
	  aggregation_urban.o              \
	  mod_region_clip.o                \
	  mksrfdata.o

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
