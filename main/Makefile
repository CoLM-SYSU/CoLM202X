# Makefile for CLM main program

include ../include/Makeoptions
INCLUDE_DIR = -I../include -I../share -I../mksrfdata \
               -I../mkinidata -I./hydro -I./bgc -I${NETCDF_INC}
VPATH = ../share : ../mksrfdata : ../main : hydro : bgc
HEADER = ../include/define.h

DEF  = $(shell grep -i cama_flood ../include/define.h)
CaMa = $(word 1, ${DEF})

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
				  ../share/mod_data_type.o      \
				  ../share/ncio_block.o  \
				  ../share/mod_catchment_data.o \
				  ../share/mod_modis_data.o \
				  ../share/mod_mesh.o   \
				  ../share/mod_pixelset.o          \
				  ../share/ncio_vector.o \
				  ../share/mod_colm_debug.o \
				  ../share/mod_mapping_grid2pset.o \
				  ../share/mod_mapping_pset2grid.o

OBJS_MKSRFDATA = \
				  ../mksrfdata/mod_aggregation_generic.o \
				  ../mksrfdata/mod_srfdata_restart.o \
				  ../mksrfdata/mod_pixelsetshadow.o \
				  ../mksrfdata/mod_landelm.o \
				  ../mksrfdata/mod_landhru.o \
				  ../mksrfdata/mod_landpatch.o \
				  ../mksrfdata/mod_aggregation_PFT.o \
				  ../mksrfdata/mod_landpft.o \
				  ../mksrfdata/mod_landpc.o  \
				  ../mksrfdata/mod_single_srfdata.o \
				  ../mksrfdata/mod_elm_vector.o \
				  ../mksrfdata/mod_hru_vector.o

OBJS_HYDRO = \
				  hydro/mod_soil_function.o \
				  hydro/mod_soil_water.o

OBJS_BGC = \
			  bgc/bgc_CNCStateUpdate1Mod.o    \
			  bgc/bgc_CNCStateUpdate2Mod.o    \
			  bgc/bgc_CNCStateUpdate3Mod.o    \
			  bgc/bgc_CNNStateUpdate1Mod.o    \
			  bgc/bgc_CNNStateUpdate2Mod.o    \
			  bgc/bgc_CNNStateUpdate3Mod.o    \
			  bgc/bgc_SoilBiogeochemNStateUpdate1Mod.o \
			  bgc/bgc_soil_SoilBiogeochemNitrifDenitrifMod.o \
			  bgc/bgc_soil_SoilBiogeochemCompetitionMod.o  \
			  bgc/bgc_soil_SoilBiogeochemDecompCascadeBGCMod.o  \
			  bgc/bgc_soil_SoilBiogeochemDecompMod.o  \
			  bgc/bgc_soil_SoilBiogeochemLittVertTranspMod.o  \
			  bgc/bgc_soil_SoilBiogeochemNLeachingMod.o  \
			  bgc/bgc_soil_SoilBiogeochemPotentialMod.o  \
			  bgc/bgc_soil_SoilBiogeochemVerticalProfileMod.o  \
			  bgc/bgc_veg_CNGapMortalityMod.o  \
			  bgc/bgc_veg_CNGRespMod.o  \
			  bgc/bgc_veg_CNMRespMod.o  \
			  bgc/bgc_DaylengthMod.o \
			  bgc/bgc_veg_CNPhenologyMod.o  \
			  bgc/bgc_veg_NutrientCompetitionMod.o  \
			  bgc/bgc_veg_CNVegStructUpdateMod.o \
			  bgc/bgc_CNSummaryMod.o  \
			  bgc/bgc_CNAnnualUpdateMod.o \
			  bgc/bgc_CNZeroFluxesMod.o \
			  bgc/bgc_CNBalanceCheckMod.o \
			  bgc/bgc_CNSASUMod.o \
			  bgc/bgc_veg_CNNDynamicsMod.o \
			  bgc/bgc_veg_CNFireBaseMod.o \
                          bgc/bgc_veg_CNFireLi2016Mod.o \
			  bgc/bgc_driver.o

BASIC =  \
		  PFT_Const.o \
		  MOD_PFTimeVars.o \
		  MOD_PFTimeInvars.o \
		  MOD_PCTimeVars.o \
		  MOD_PCTimeInvars.o \
		  MOD_TimeInvariants.o        \
		  MOD_TimeVariables.o    \
		  MOD_1D_PFTFluxes.o          \
		  MOD_1D_PCFluxes.o           \
		  MOD_1D_Fluxes.o \
		  MOD_1D_Forcing.o \
		  monthly_in_situ_co2_mlo.o

BASIC_BGC = \
			  bgc/MOD_1D_BGCFluxes.o  \
			  bgc/MOD_1D_BGCPFTFluxes.o  \
			  bgc/MOD_2D_BGCFluxes.o  \
			  bgc/MOD_BGCPFTimeVars.o  \
			  bgc/MOD_BGCTimeInvars.o  \
			  bgc/MOD_BGCTimeVars.o 

OBJS =  \
		  MOD_2D_Fluxes.o             \
		  MOD_2D_Forcing.o            \
		  user_specified_forcing.o    \
		  DownscalingForcingMod.o     \
		  mod_forcing.o               \
		  ThreeDCanopy.o              \
		  AerosolMod.o                \
		  SnowSnicarMod.o             \
		  SnowAlbedo.o                \
		  ALBEDO.o                    \
		  ASSIM_STOMATA_conductance.o \
		  PlantHydraulic.o            \
		  FRICTION_VELOCITY.o         \
		  OzoneMod.o                  \
		  LEAF_temperature.o          \
		  LEAF_temperature_PC.o       \
		  soil_hcap_cond.o            \
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
		  CROP_readin.o               \
		  NITRIF_readin.o             \
		  NDEP_readin.o               \
		  Fire_readin.o               \
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
		  MOD_CaMa_Variables.o        \
		  mod_hist_vector.o           \
		  mod_hist.o                  \
		  mod_lightning_data.o        \
		  mod_ozone_data.o            \
		  colm_CaMaMod.o              \
		  CLMDRIVER.o                 \
		  CLMMAIN.o                   \
		  CLM.o                       

####################################################################
ifneq (${CaMa},\#define)
# Compile CoLM decoupled without river routing scheme (CaMa-Flood)

${EXECUTABLE} : ${HEADER} \
	${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC_BGC} ${BASIC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS} 
	\
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC_BGC} ${BASIC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS} \
	-o $@ ${LDFLAGS}
	@echo 'making CLM completed!'

$(OBJS) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

else
####################################################################
# The global river model CaMa-Flood (version 4.0.1)
CaMa_DIR = ../CaMa

# CaMa Flood Model modules directories
CaMa_MODS = -I$(CaMa_DIR)/src

# CaMa Flood Model libs (static) directories
CaMa_LIBS = $(CaMa_DIR)/src/libcama.a
#--------------------------------------------

${EXECUTABLE} : ${HEADER} \
	${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC_BGC} ${BASIC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS} mk_CaMa
	\
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC_BGC} ${BASIC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS} \
	${CaMa_LIBS} -o $@ ${LDFLAGS}
	@echo 'making CLM with CaMa Flood Model completed!'

$(OBJS) : %.o : %.F90 ${HEADER} mk_CaMa
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${CaMa_MODS}

mk_CaMa : 
	cd ../CaMa/src && make

endif

$(OBJS_HYDRO) : %.o : %.F90 ${HEADER} ${BASIC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}hydro

$(BASIC_BGC) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}bgc

$(OBJS_BGC) : %.o : %.F90 ${HEADER} ${BASIC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}bgc

$(BASIC) : %.o : %.F90 ${HEADER} ${BASIC_BGC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

$(OBJS_SHARED) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../share

$(OBJS_MKSRFDATA) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../mksrfdata

.PHONY: clean
clean :
	${RM} -f *.o *.mod ${EXECUTABLE}
	cd hydro && ${RM} -f *.o *.mod
	cd bgc   && make clean
