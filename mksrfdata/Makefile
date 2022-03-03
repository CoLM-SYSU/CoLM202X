# Makefile for CLM surface data processing

include ../include/Makeoptions

VPATH = ../share

# name of executable
EXECUTABLE = mksrfdata.x

####################################################################
.DEFAULT :

HEADER = ../include/define.h

OBJS_SHARED =  ../share/precision.o \
					../share/MathConstants.o \
					../share/spmd_task.o \
					../share/mod_namelist.o \
					../share/mod_utils.o \
					../share/ncio_serial.o \
					../share/mod_block.o \
					../share/mod_grid.o  \
					../share/mod_pixel.o \
					../share/mod_data_type.o \
					../share/ncio_block.o \
					../share/mod_hydro_data.o \
					../share/mod_landunit.o \
					../share/mod_pixelset.o \
					../share/ncio_vector.o \
					../share/mod_colm_debug.o \
					../share/GlobalVars.o \
					../share/LC_Const.o

OBJS = \
		 mod_srfdata_restart.o \
		 mod_modis_data.o \
		 mod_landcell.o \
		 mod_landpatch.o \
		 mod_aggregation_LC.o \
		 mod_aggregation_PFT.o \
		 mod_landpft.o \
		 mod_landpc.o \
		 aggregation_landtype.o \
		 aggregation_LAI.o \
		 aggregation_soil_brightness.o \
		 aggregation_lakedepth.o \
		 aggregation_forest_height.o \
		 aggregation_soil_parameters.o \
		 aggregation_dbedrock.o \
		 mksrfdata.o

####################################################################

${EXECUTABLE} : ${HEADER} ${OBJS_SHARED} ${OBJS}  
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS} -o $@ ${LDFLAGS}
	@echo 'making CLM surface data completed!'

$(OBJS) : %.o : %.F90 ${HEADER} ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

$(OBJS_SHARED) : %.o : %.F90 ${HEADER} 
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../share

clean :
	${RM} -f *.o *.mod
