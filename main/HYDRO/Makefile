# Makefile for CoLM surface data processing

include ../../include/Makeoptions
 
INCLUDE_DIR = -I../../include -I../../share -I../../mksrfdata -I.. \
				  -I${NETCDF_INC}

VPATH  = ../../share  ../../mksrfdata ..
HEADER = ../../include/define.h

####################################################################
.DEFAULT :

OBJS_SHARED = \
				  ../../share/precision.o \
				  ../../share/MathConstants.o \
				  ../../share/PhysicalConstants.o \
				  ../../share/spmd_task.o \
				  ../../share/mod_utils.o \
				  ../../share/timemanager.o \
				  ../../share/ncio_serial.o \
				  ../../share/mod_namelist.o \
				  ../../share/mod_block.o \
				  ../../share/mod_grid.o           \
				  ../../share/mod_pixel.o          \
				  ../../share/mod_data_type.o      \
				  ../../share/ncio_block.o \
				  ../../share/mod_catchment_data.o \
				  ../../share/mod_landbasin.o \
				  ../../share/mod_pixelset.o \
				  ../../share/GlobalVars.o \
				  ../../share/LC_Const.o \
				  ../../share/ncio_vector.o \
				  ../../share/mod_colm_debug.o

OBJS_MKSRFDATA = \
				  ../../mksrfdata/mod_modis_data.o \
				  ../../mksrfdata/mod_hydrounit.o \
				  ../../mksrfdata/mod_landpatch.o \
				  ../../mksrfdata/mod_aggregation_PFT.o \
				  ../../mksrfdata/mod_landpft.o \
				  ../../mksrfdata/mod_landpc.o

BASIC =  \
			../mod_unstructured_mesh.o \
		  ../MOD_Vars_PFTimeInvars.o \
		  ../MOD_Vars_PCTimeInvars.o \
		  ../MOD_Vars_TimeInvariants.o \
		  ../MOD_Vars_PFTimeVars.o \
		  ../MOD_Vars_PCTimeVars.o \
		  ../MOD_Vars_TimeVariables.o \
		  ../MOD_Vars_1DPFTFluxes.o \
		  ../MOD_Vars_1DPCFluxes.o \
		  ../MOD_Vars_1DFluxes.o

OBJS = \
		 mod_io_basin.o \
		 MOD_HydroTimeVars.o \
		 MOD_1D_HydroFluxes.o \
		 mod_soil_function.o \
		 mod_soil_water.o    \
		 mod_drainage_network.o \
		 mod_river_network.o \
		 mod_surface_runoff.o \
		 mod_subsurface_runoff.o \
		 mod_river_flow.o \
		 mod_hist_basin.o \
		 mod_lateral_flow.o

####################################################################

main : ${OBJS} ${OBJS_SHARED} ${OBJS_MKSRFDATA} ${BASIC}
	@echo 'Compiling modules in hydro completed!'

$(OBJS):%.o : %.F90 ${HEADER} ${OBJS_SHARED} ${OBJS_MKSRFDATA} ${BASIC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

$(OBJS_SHARED) : %.o : %.F90 ${HEADER} 
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../../share

$(OBJS_MKSRFDATA) : %.o : %.F90 ${HEADER} ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../../mksrfdata

$(BASIC) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../

clean :
	${RM} -f *.o *.mod
