# Makefile for CLM initial data processing

include ../include/Makeoptions

INCLUDE_DIR = -I../include -I../share -I../mksrfdata -I../main \
	  -I../main/hydro -I../main/bgc -I../main/urban -I${NETCDF_INC}
VPATH = ../share : ../mksrfdata : ../main : ../main/hydro : ../main/bgc : ../main/urban
HEADER = ../include/define.h

# name of executable
EXECUTABLE = initial.x

####################################################################

OBJS_SHARED =                                 \
	  ../share/precision.o                \
	  ../share/spmd_task.o                \
	  ../share/mod_namelist.o             \
	  ../share/PhysicalConstants.o        \
	  ../share/MathConstants.o            \
	  ../share/timemanager.o              \
	  ../share/mod_utils.o                \
	  ../share/ncio_serial.o              \
	  ../share/mod_block.o                \
	  ../share/mod_grid.o                 \
	  ../share/mod_pixel.o                \
	  ../share/mod_data_type.o            \
	  ../share/ncio_block.o               \
	  ../share/mod_catchment_data.o       \
	  ../share/mod_5x5_data.o             \
	  ../share/mod_mesh.o                 \
	  ../share/mod_pixelset.o             \
	  ../share/ncio_vector.o              \
	  ../share/mod_colm_debug.o           \
	  ../share/GlobalVars.o               \
	  ../share/LC_Const.o                 \
	  ../share/UrbanLCZ_Const.o           \
	  ../share/mod_mapping_grid2pset.o

OBJS_MKSRFDATA =                              \
	  ../mksrfdata/mod_aggregation.o      \
	  ../mksrfdata/mod_srfdata_restart.o  \
	  ../mksrfdata/mod_landelm.o          \
	  ../mksrfdata/mod_landhru.o          \
	  ../mksrfdata/mod_pixelsetshadow.o   \
	  ../mksrfdata/mod_landpatch.o        \
	  ../mksrfdata/mod_landpft.o          \
	  ../mksrfdata/mod_landpc.o           \
	  ../mksrfdata/mod_single_srfdata.o   \
	  ../mksrfdata/mod_elm_vector.o       \
	  ../mksrfdata/mod_hru_vector.o       \
	  ../mksrfdata/mod_landurban.o

BASIC_HYDRO =                                 \
	  ../main/hydro/mod_io_basin.o        \
	  ../main/hydro/MOD_HydroTimeVars.o   \
	  ../main/hydro/MOD_1D_HydroFluxes.o  \
	  ../main/hydro/mod_soil_function.o   \
	  ../main/hydro/mod_soil_water.o

BASIC_BGC =                                   \
	  ../main/bgc/MOD_1D_BGCFluxes.o      \
	  ../main/bgc/MOD_1D_BGCPFTFluxes.o   \
	  ../main/bgc/MOD_2D_BGCFluxes.o      \
	  ../main/bgc/MOD_BGCPFTimeVars.o     \
	  ../main/bgc/MOD_BGCTimeInvars.o     \
	  ../main/bgc/MOD_BGCTimeVars.o

BASIC_URBAN =                                 \
	  ../main/urban/MOD_1D_UrbanFluxes.o  \
	  ../main/urban/MOD_UrbanTimeVars.o   \
	  ../main/urban/MOD_UrbanTimeInvars.o \

BASIC =                                       \
	  ../main/PFT_Const.o                 \
	  ../main/MOD_PFTimeVars.o            \
	  ../main/MOD_PFTimeInvars.o          \
	  ../main/MOD_PCTimeVars.o            \
	  ../main/MOD_PCTimeInvars.o          \
	  ../main/MOD_TimeInvariants.o        \
	  ../main/MOD_TimeVariables.o         \
	  ../main/MOD_1D_PFTFluxes.o          \
	  ../main/MOD_1D_PCFluxes.o           \
	  ../main/MOD_1D_Fluxes.o             \
	  ../main/MOD_1D_Forcing.o

OBJS_MAIN =                                   \
	  ../main/LAI_empirical.o             \
	  ../main/LAI_readin.o                \
	  ../main/orb_coszen.o                \
	  ../main/ThreeDCanopy.o              \
	  ../main/SnowSnicarMod.o             \
	  ../main/snowage.o                   \
	  ../main/SnowAlbedo.o                \
	  ../main/AerosolMod.o                \
	  ../main/ALBEDO.o                    \
	  ../main/snowfraction.o              \
	  ../main/urban/UrbanLAI_readin.o     \
	  ../main/urban/UrbanShortwave.o      \
	  ../main/urban/UrbanALBEDO.o

OBJS =                                        \
	  pct_readin.o                        \
	  lakedepth_readin.o                  \
	  dbedrock_readin.o                   \
	  soil_color_refl.o                   \
	  soil_parameters_readin.o            \
	  HTOP_readin.o                       \
	  Urban_readin.o                      \
	  ../main/CROP_readin.o               \
	  ../main/NITRIF_readin.o             \
	  ../main/NDEP_readin.o               \
	  ../main/Fire_readin.o               \
	  IniTimeVar.o                        \
	  UrbanIniTimeVar.o                   \
	  initialize.o                        \
	  CLMINI.o

####################################################################

all: ${EXECUTABLE}

${EXECUTABLE} : ${HEADER} ${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC} ${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${OBJS_MAIN} ${OBJS}
	\
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC} ${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${OBJS_MAIN} ${OBJS} -o $@ ${LDFLAGS}
	@echo 'making CLM initial data completed!'

$(OBJS) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

$(OBJS_SHARED) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../share

$(OBJS_MKSRFDATA) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../mksrfdata

$(BASIC_BGC) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../main/bgc

$(BASIC_HYDRO) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../main/hydro

$(BASIC_URBAN) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../main/urban

$(BASIC) : %.o : %.F90 ${HEADER} ${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../main

$(OBJS_MAIN) : %.o : %.F90 ${HEADER} ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../main

clean :
	${RM} -f *.o *.mod ${EXECUTABLE}
