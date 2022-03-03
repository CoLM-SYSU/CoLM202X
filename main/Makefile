# Makefile for CLM main program

include ../include/Makeoptions

VPATH = ../share : ../mksrfdata : ../main

HEADER = ../include/define.h

# Name of executable
EXECUTABLE = clm.x

####################################################################
.DEFAULT :

OBJS_SHARED = \
				  ../share/precision.o \
				  ../share/GlobalVars.o \
				  ../share/LC_Const.o \
				  ../share/spmd_task.o          \
				  ../share/mod_namelist.o \
				  ../share/PhysicalConstants.o  \
				  ../share/MathConstants.o  \
				  ../share/mod_utils.o          \
				  ../share/timemanager.o        \
				  ../share/ncio_serial.o \
				  ../share/mod_block.o          \
				  ../share/mod_grid.o           \
				  ../share/mod_pixel.o          \
				  ../share/mod_pixelset.o          \
				  ../share/mod_data_type.o      \
				  ../share/ncio_block.o  \
				  ../share/mod_hydro_data.o \
				  ../share/mod_landunit.o   \
				  ../share/ncio_vector.o \
				  ../share/mod_colm_debug.o

OBJS_MKSRFDATA = \
				  ../mksrfdata/mod_srfdata_restart.o \
				  ../mksrfdata/mod_modis_data.o \
				  ../mksrfdata/mod_aggregation_PFT.o \
				  ../mksrfdata/mod_landcell.o \
				  ../mksrfdata/mod_landpatch.o \
				  ../mksrfdata/mod_landpft.o \
				  ../mksrfdata/mod_landpc.o

OBJS =  \
		  PFT_Const.o \
		  MOD_PFTimeVars.o \
		  MOD_PFTimeInvars.o \
		  MOD_PCTimeVars.o \
		  MOD_PCTimeInvars.o \
		  MOD_TimeInvariants.o        \
		  MOD_TimeVariables.o         \
		  mod_mapping_grid2pset.o \
		  mod_mapping_pset2grid.o \
		  MOD_1D_PFTFluxes.o          \
		  MOD_1D_PCFluxes.o           \
		  MOD_1D_Fluxes.o             \
		  MOD_1D_Forcing.o            \
		  MOD_2D_Fluxes.o             \
		  MOD_2D_Forcing.o            \
		  user_specified_forcing.o    \
		  mod_forcing.o               \
		  ThreeDCanopy.o              \
		  ALBEDO.o                    \
		  ASSIM_STOMATA_conductance.o \
		  FRICTION_VELOCITY.o         \
		  LEAF_temperature.o          \
		  LEAF_temperature_PC.o       \
		  SOIL_thermal_parameters.o   \
		  SOIL_SNOW_hydrology.o       \
		  SNOW_Layers_CombineDivide.o \
		  GLACIER.o                   \
		  LAKE.o                      \
		  SIMPLE_OCEAN.o              \
		  eroot.o                     \
		  groundfluxes.o              \
		  groundtem.o                 \
		  LAI_empirical.o             \
		  LAI_readin.o                \
		  LEAF_interception.o         \
		  meltf.o                     \
		  netsolar.o                  \
		  rain_snow_temp.o            \
		  newsnow.o                   \
		  orb_coszen.o                \
		  qsadv.o                     \
		  snowage.o                   \
		  snowfraction.o              \
		  THERMAL.o                   \
		  tridia.o                    \
		  wetbulb.o                   \
		  MOD_1D_Acc_Fluxes.o         \
		  mod_hist.o                  \
		  CLMDRIVER.o                 \
		  CLMMAIN.o                   \
		  CLM.o                       

####################################################################

# Compile CoLM decoupled with river routing scheme (CaMa-Flood)

${EXECUTABLE} : ${HEADER} ${OBJS_SHARED} ${OBJS_MKSRFDATA} ${OBJS} 
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS_MKSRFDATA} ${OBJS} -o $@ ${LDFLAGS}
	@echo 'making CLM completed!'

$(OBJS) : %.o : %.F90 ${HEADER} 
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

$(OBJS_SHARED) : %.o : %.F90 ${HEADER} 
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../share

$(OBJS_MKSRFDATA) : %.o : %.F90 ${HEADER} 
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../mksrfdata

####################################################################

clean :
	${RM} -f *.o *.mod
