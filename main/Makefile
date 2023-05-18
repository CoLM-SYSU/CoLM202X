# Makefile for CLM main program

include ../include/Makeoptions
INCLUDE_DIR = -I../include -I../share -I../mksrfdata \
               -I../mkinidata -I./hydro -I./bgc -I./urban \
	       -I./lulcc -I${NETCDF_INC}
VPATH = ../share : ../mksrfdata : ../main : hydro : bgc : urban : lulcc
HEADER = ../include/define.h

DEF  = $(shell grep -i cama_flood ../include/define.h)
CaMa = $(word 1, ${DEF})

# Name of executable
EXECUTABLE = clm.x

####################################################################
.DEFAULT :

OBJS_SHARED =                                              \
	  ../share/precision.o                             \
	  ../share/GlobalVars.o                            \
	  ../share/LC_Const.o                              \
	  ../share/UrbanLCZ_Const.o                        \
	  ../share/spmd_task.o                             \
	  ../share/mod_namelist.o                          \
	  ../share/PhysicalConstants.o                     \
	  ../share/MathConstants.o                         \
	  ../share/mod_utils.o                             \
	  ../share/timemanager.o                           \
	  ../share/ncio_serial.o                           \
	  ../share/mod_block.o                             \
	  ../share/mod_grid.o                              \
	  ../share/mod_pixel.o                             \
	  ../share/mod_data_type.o                         \
	  ../share/ncio_block.o                            \
	  ../share/mod_catchment_data.o                    \
	  ../share/mod_5x5_data.o                          \
	  ../share/mod_mesh.o                              \
	  ../share/mod_pixelset.o                          \
	  ../share/ncio_vector.o                           \
	  ../share/mod_colm_debug.o                        \
	  ../share/mod_mapping_grid2pset.o                 \
	  ../share/mod_mapping_pset2grid.o

OBJS_MKSRFDATA =                                           \
	  ../mksrfdata/mod_aggregation.o                   \
	  ../mksrfdata/mod_srfdata_restart.o               \
	  ../mksrfdata/mod_pixelsetshadow.o                \
	  ../mksrfdata/mod_landelm.o                       \
	  ../mksrfdata/mod_landhru.o                       \
	  ../mksrfdata/mod_landpatch.o                     \
	  ../mksrfdata/mod_landpft.o                       \
	  ../mksrfdata/mod_landpc.o                        \
	  ../mksrfdata/mod_landurban.o                     \
	  ../mksrfdata/mod_single_srfdata.o                \
	  ../mksrfdata/mod_elm_vector.o                    \
	  ../mksrfdata/mod_hru_vector.o

OBJS_MKINIDATA =                                           \
		../mkinidata/pct_readin.o                  \
		../mkinidata/HTOP_readin.o                 \
		../mkinidata/lakedepth_readin.o            \
		../mkinidata/soil_parameters_readin.o      \
		../mkinidata/Urban_readin.o                \
		../mkinidata/UrbanIniTimeVar.o             \
		../mkinidata/IniTimeVar.o                  

BASIC_HYDRO =                                              \
	  hydro/mod_io_basin.o                             \
	  hydro/MOD_HydroTimeVars.o                        \
	  hydro/MOD_1D_HydroFluxes.o

OBJS_HYDRO =                                               \
	  hydro/mod_soil_function.o                        \
	  hydro/mod_soil_water.o                           \
	  hydro/mod_drainage_network.o                     \
	  hydro/mod_river_network.o                        \
	  hydro/mod_ssrf_network.o                         \
	  hydro/mod_surface_runoff.o                       \
	  hydro/mod_subsurface_runoff.o                    \
	  hydro/mod_river_flow.o                           \
	  hydro/mod_hist_basin.o                           \
	  hydro/mod_lateral_flow.o

OBJS_BGC = \
			  bgc/MOD_BGC_CNCStateUpdate1.o    \
			  bgc/MOD_BGC_CNCStateUpdate2.o    \
			  bgc/MOD_BGC_CNCStateUpdate3.o    \
			  bgc/MOD_BGC_CNNStateUpdate1.o    \
			  bgc/MOD_BGC_CNNStateUpdate2.o    \
			  bgc/MOD_BGC_CNNStateUpdate3.o    \
			  bgc/MOD_BGC_Soil_BiogeochemNStateUpdate1.o \
			  bgc/MOD_BGC_Soil_BiogeochemNitrifDenitrif.o \
			  bgc/MOD_BGC_Soil_BiogeochemCompetition.o  \
			  bgc/MOD_BGC_Soil_BiogeochemDecompCascadeBGC.o  \
			  bgc/MOD_BGC_Soil_BiogeochemDecomp.o  \
			  bgc/MOD_BGC_Soil_BiogeochemLittVertTransp.o  \
			  bgc/MOD_BGC_Soil_BiogeochemNLeaching.o  \
			  bgc/MOD_BGC_Soil_BiogeochemPotential.o  \
			  bgc/MOD_BGC_Soil_BiogeochemVerticalProfile.o  \
			  bgc/MOD_BGC_Veg_CNGapMortality.o  \
			  bgc/MOD_BGC_Veg_CNGResp.o  \
			  bgc/MOD_BGC_Veg_CNMResp.o  \
			  bgc/MOD_BGC_Daylength.o \
			  bgc/MOD_BGC_Veg_CNPhenology.o  \
			  bgc/MOD_BGC_Veg_NutrientCompetition.o  \
			  bgc/MOD_BGC_Veg_CNVegStructUpdate.o \
			  bgc/MOD_BGC_CNSummary.o  \
			  bgc/MOD_BGC_CNAnnualUpdate.o \
			  bgc/MOD_BGC_CNZeroFluxes.o \
			  bgc/MOD_BGC_CNBalanceCheck.o \
			  bgc/MOD_BGC_CNSASU.o \
			  bgc/MOD_BGC_Veg_CNNDynamics.o \
			  bgc/MOD_BGC_Veg_CNFireBase.o \
                          bgc/MOD_BGC_Veg_CNFireLi2016.o \
			  bgc/MOD_BGC_driver.o

OBJS_URBAN =                                               \
	  urban/UrbanShortwave.o                           \
	  urban/UrbanALBEDO.o                              \
	  urban/UrbanLAI_readin.o                          
OBJS_LULCC =                                               \
	  lulcc/LuLccInitialize.o                          \
          lulcc/LuLccDRIVER.o	  
BASIC =                                                    \
	  PFT_Const.o                                      \
	  MOD_PFTimeVars.o                                 \
	  MOD_PFTimeInvars.o                               \
	  MOD_PCTimeVars.o                                 \
	  MOD_PCTimeInvars.o                               \
	  MOD_TimeInvariants.o                             \
	  MOD_TimeVariables.o                              \
	  MOD_1D_PFTFluxes.o                               \
	  MOD_1D_PCFluxes.o                                \
	  MOD_1D_Fluxes.o                                  \
	  MOD_1D_Forcing.o                                 \
	  monthly_in_situ_co2_mlo.o

BASIC_BGC =                                                \
	  bgc/MOD_BGC_Vars_1DFluxes.o                           \
	  bgc/MOD_BGC_Vars_1DPFTFluxes.o                        \
	  bgc/MOD_BGC_Vars_2DFluxes.o                           \
	  bgc/MOD_BGC_Vars_PFTimeVars.o                          \
	  bgc/MOD_BGC_Vars_TimeInvars.o                          \
	  bgc/MOD_BGC_Vars_TimeVars.o

BASIC_URBAN =                                              \
	  urban/MOD_UrbanTimeInvars.o                      \
	  urban/MOD_UrbanTimeVars.o
BASIC_LULCC =                                              \
	  lulcc/MOD_LuLccTimeInvars.o                      \
	  lulcc/MOD_LuLccTimeVars.o
OBJS =                                                     \
	  MOD_2D_Fluxes.o                                  \
	  MOD_2D_Forcing.o                                 \
	  user_specified_forcing.o                         \
	  DownscalingForcingMod.o                          \
	  mod_forcing.o                                    \
	  ThreeDCanopy.o                                   \
	  AerosolMod.o                                     \
	  SnowSnicarMod.o                                  \
	  SnowAlbedo.o                                     \
	  ALBEDO.o                                         \
	  ASSIM_STOMATA_conductance.o                      \
	  PlantHydraulic.o                                 \
	  FRICTION_VELOCITY.o                              \
	  OzoneMod.o                                       \
	  LEAF_temperature.o                               \
	  LEAF_temperature_PC.o                            \
	  soil_hcap_cond.o                                 \
	  SOIL_SNOW_hydrology.o                            \
	  SNOW_Layers_CombineDivide.o                      \
	  GLACIER.o                                        \
	  LAKE.o                                           \
	  SIMPLE_OCEAN.o                                   \
	  eroot.o                                          \
	  groundfluxes.o                                   \
	  groundtem.o                                      \
	  LAI_empirical.o                                  \
	  LAI_readin.o                                     \
	  CROP_readin.o                                    \
	  NITRIF_readin.o                                  \
	  NDEP_readin.o                                    \
	  Fire_readin.o                                    \
	  LEAF_interception.o                              \
	  meltf.o                                          \
	  netsolar.o                                       \
	  rain_snow_temp.o                                 \
	  newsnow.o                                        \
	  orb_coszen.o                                     \
	  qsadv.o                                          \
	  snowage.o                                        \
	  snowfraction.o                                   \
	  THERMAL.o                                        \
	  tridia.o                                         \
	  wetbulb.o                                        \
	  MOD_1D_Acc_Fluxes.o                              \
	  MOD_CaMa_Variables.o                             \
	  mod_hist_vector.o                                \
	  mod_hist.o                                       \
	  mod_lightning_data.o                             \
	  mod_ozone_data.o                                 \
	  colm_CaMaMod.o                                   \
	  CLMDRIVER.o                                      \
	  CLMMAIN.o                                        \
	  CLM.o

####################################################################
ifneq (${CaMa},\#define)
# Compile CoLM decoupled without river routing scheme (CaMa-Flood)

${EXECUTABLE} : ${HEADER} \
	${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${BASIC_LULCC} ${BASIC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS_URBAN} ${OBJS_MKINIDATA} ${OBJS_LULCC} ${OBJS} 
	\
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${BASIC_LULCC} ${BASIC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS_URBAN} ${OBJS_MKINIDATA} ${OBJS_LULCC} ${OBJS} \
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
	${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${BASIC_LULCC} ${BASIC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS_URBAN} ${OBJS_LULCC} ${OBJS} mk_CaMa
	\
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${BASIC} ${BASIC_LULCC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS_URBAN} ${OBJS_LULCC} ${OBJS} \
	${CaMa_LIBS} -o $@ ${LDFLAGS}
	@echo 'making CLM with CaMa Flood Model completed!'

$(OBJS) : %.o : %.F90 ${HEADER} mk_CaMa
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${CaMa_MODS}

mk_CaMa :
	cd ../CaMa/src && make

endif

$(BASIC_HYDRO) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}hydro

$(OBJS_HYDRO) : %.o : %.F90 ${HEADER} ${BASIC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}hydro

$(BASIC_BGC) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}bgc

$(OBJS_BGC) : %.o : %.F90 ${HEADER} ${BASIC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}bgc

$(BASIC_URBAN) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}urban

$(OBJS_URBAN) : %.o : %.F90 ${HEADER} ${BASIC} ${OBJS}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}urban

$(BASIC_LULCC) : %.o : %.F90 ${HEADER} ${BASIC} $(BASIC_URBAN)
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}lulcc

$(OBJS_MKINIDATA) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../mkinidata

$(OBJS_LULCC) : %.o : %.F90 ${HEADER} ${BASIC} $(BASIC_URBAN) $(OBJS_MKINIDATA)
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}lulcc
	
$(BASIC) : %.o : %.F90 ${HEADER} ${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

$(OBJS_SHARED) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../share

$(OBJS_MKSRFDATA) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../mksrfdata

.PHONY: clean
clean :
	${RM} -f *.o *.mod ${EXECUTABLE}
	cd hydro && ${RM} -f *.o *.mod
	cd urban && ${RM} -f *.o *.mod
	cd bgc   && make clean
