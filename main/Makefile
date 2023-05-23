# Makefile for CoLM main program

include ../include/Makeoptions
INCLUDE_DIR = -I../include -I../share -I../mksrfdata \
               -I../mkinidata -I./hydro -I./bgc -I./urban -I${NETCDF_INC}
VPATH = ../share : ../mksrfdata : ../main : hydro : bgc : urban
HEADER = ../include/define.h

DEF  = $(shell grep -i cama_flood ../include/define.h)
CaMa = $(word 1, ${DEF})

# Name of executable
EXECUTABLE = colm.x

####################################################################
.DEFAULT :

OBJS_SHARED =                                           \
	  ../share/precision.o                          \
	  MOD_Vars_Global.o                             \
	  MOD_Vars_LCConst.o                            \
	  MOD_Vars_PhysicalConst.o                      \
	  ../share/spmd_task.o                          \
	  ../share/mod_namelist.o                       \
	  ../share/MathConstants.o                      \
	  ../share/mod_utils.o                          \
	  ../share/timemanager.o                        \
	  ../share/ncio_serial.o                        \
	  ../share/mod_block.o                          \
	  ../share/mod_grid.o                           \
	  ../share/mod_pixel.o                          \
	  ../share/mod_data_type.o                      \
	  ../share/ncio_block.o                         \
	  ../share/mod_catchment_data.o                 \
	  ../share/mod_5x5_data.o                       \
	  ../share/mod_mesh.o                           \
	  ../share/mod_pixelset.o                       \
	  ../share/ncio_vector.o                        \
	  ../share/mod_colm_debug.o                     \
	  ../share/mod_mapping_grid2pset.o              \
	  ../share/mod_mapping_pset2grid.o              \
	  MOD_Eroot.o                                   \
	  MOD_Qsadv.o

OBJS_MKSRFDATA =                                        \
	  ../mksrfdata/mod_aggregation.o                \
	  ../mksrfdata/mod_srfdata_restart.o            \
	  ../mksrfdata/mod_pixelsetshadow.o             \
	  ../mksrfdata/mod_landelm.o                    \
	  ../mksrfdata/mod_landhru.o                    \
	  ../mksrfdata/mod_landpatch.o                  \
	  ../mksrfdata/mod_landpft.o                    \
	  ../mksrfdata/mod_landpc.o                     \
	  ../mksrfdata/mod_landurban.o                  \
	  ../mksrfdata/mod_single_srfdata.o             \
	  ../mksrfdata/mod_elm_vector.o                 \
	  ../mksrfdata/mod_hru_vector.o

BASIC_HYDRO =                                           \
	  hydro/mod_io_basin.o                          \
	  hydro/MOD_HydroTimeVars.o                     \
	  hydro/MOD_1D_HydroFluxes.o

OBJS_HYDRO =                                            \
	  hydro/mod_soil_function.o                     \
	  hydro/mod_soil_water.o                        \
	  hydro/mod_drainage_network.o                  \
	  hydro/mod_river_network.o                     \
	  hydro/mod_ssrf_network.o                      \
	  hydro/mod_surface_runoff.o                    \
	  hydro/mod_subsurface_runoff.o                 \
	  hydro/mod_river_flow.o                        \
	  hydro/mod_hist_basin.o                        \
	  hydro/mod_lateral_flow.o

OBJS_BGC =                                              \
	  bgc/MOD_BGC_CNCStateUpdate1.o                 \
	  bgc/MOD_BGC_CNCStateUpdate2.o                 \
	  bgc/MOD_BGC_CNCStateUpdate3.o                 \
	  bgc/MOD_BGC_CNNStateUpdate1.o                 \
	  bgc/MOD_BGC_CNNStateUpdate2.o                 \
	  bgc/MOD_BGC_CNNStateUpdate3.o                 \
	  bgc/MOD_BGC_Soil_BiogeochemNStateUpdate1.o    \
	  bgc/MOD_BGC_Soil_BiogeochemNitrifDenitrif.o   \
	  bgc/MOD_BGC_Soil_BiogeochemCompetition.o      \
	  bgc/MOD_BGC_Soil_BiogeochemDecompCascadeBGC.o \
	  bgc/MOD_BGC_Soil_BiogeochemDecomp.o           \
	  bgc/MOD_BGC_Soil_BiogeochemLittVertTransp.o   \
	  bgc/MOD_BGC_Soil_BiogeochemNLeaching.o        \
	  bgc/MOD_BGC_Soil_BiogeochemPotential.o        \
	  bgc/MOD_BGC_Soil_BiogeochemVerticalProfile.o  \
	  bgc/MOD_BGC_Veg_CNGapMortality.o              \
	  bgc/MOD_BGC_Veg_CNGResp.o                     \
	  bgc/MOD_BGC_Veg_CNMResp.o                     \
	  bgc/MOD_BGC_Daylength.o                       \
	  bgc/MOD_BGC_Veg_CNPhenology.o                 \
	  bgc/MOD_BGC_Veg_NutrientCompetition.o         \
	  bgc/MOD_BGC_Veg_CNVegStructUpdate.o           \
	  bgc/MOD_BGC_CNSummary.o                       \
	  bgc/MOD_BGC_CNAnnualUpdate.o                  \
	  bgc/MOD_BGC_CNZeroFluxes.o                    \
	  bgc/MOD_BGC_CNBalanceCheck.o                  \
	  bgc/MOD_BGC_CNSASU.o                          \
	  bgc/MOD_BGC_Veg_CNNDynamics.o                 \
	  bgc/MOD_BGC_Veg_CNFireBase.o                  \
	  bgc/MOD_BGC_Veg_CNFireLi2016.o                \
	  bgc/MOD_BGC_driver.o

OBJS_URBAN =                                            \
	  urban/MOD_Urban_Shortwave.o                   \
	  urban/MOD_Urban_Longwave.o                    \
	  urban/MOD_Urban_Albedo.o                      \
	  urban/MOD_Urban_NetSolar.o                    \
	  urban/MOD_Urban_Flux.o                        \
	  urban/MOD_Urban_GroundFlux.o                  \
	  urban/MOD_Urban_RoofFlux.o                    \
	  urban/MOD_Urban_RoofTem.o                     \
	  urban/MOD_Urban_WallTem.o                     \
	  urban/MOD_Urban_PerviousTem.o                 \
	  urban/MOD_Urban_ImperviousTem.o               \
	  urban/MOD_Urban_Hydrology.o                   \
	  urban/MOD_Urban_BEM.o                         \
	  urban/MOD_Urban_LUCY.o                        \
	  urban/MOD_Urban_Thermal.o                     \
	  urban/Urban_CoLMMAIN.o

BASIC =                                                 \
	  MOD_Vars_PFTConst.o                           \
	  MOD_Vars_PFTimeVars.o                         \
	  MOD_Vars_PFTimeInvars.o                       \
	  MOD_Vars_PCTimeVars.o                         \
	  MOD_Vars_PCTimeInvars.o                       \
	  MOD_Vars_TimeInvariants.o                     \
	  MOD_Vars_TimeVariables.o                      \
	  MOD_Vars_1DPFTFluxes.o                        \
	  MOD_Vars_1DPCFluxes.o                         \
	  MOD_Vars_1DFluxes.o                           \
	  MOD_Vars_1DForcing.o                          \
	  MOD_MonthlyinSituCO2mlo.o

BASIC_BGC =                                             \
	  bgc/MOD_BGC_Vars_1DFluxes.o                   \
	  bgc/MOD_BGC_Vars_1DPFTFluxes.o                \
	  bgc/MOD_BGC_Vars_2DFluxes.o                   \
	  bgc/MOD_BGC_Vars_PFTimeVars.o                 \
	  bgc/MOD_BGC_Vars_TimeInvars.o                 \
	  bgc/MOD_BGC_Vars_TimeVars.o

BASIC_URBAN =                                           \
	  urban/MOD_Urban_Vars_LCZConst.o               \
	  urban/MOD_Urban_Vars_1DFluxes.o               \
	  urban/MOD_Urban_Vars_TimeInvars.o             \
	  urban/MOD_Urban_Vars_TimeVars.o

OBJS =                                                  \
	  MOD_Vars_2DFluxes.o                           \
	  MOD_Vars_2DForcing.o                          \
	  MOD_OrbCoszen.o                               \
	  MOD_UserSpecifiedForcing.o                    \
	  MOD_DownscalingForcing.o                      \
	  MOD_Forcing.o                                 \
	  MOD_ThreeDCanopy.o                            \
	  MOD_Aerosol.o                                 \
	  MOD_SnowSnicar.o                              \
	  MOD_Albedo.o                                  \
	  MOD_AssimStomataConductance.o                 \
	  MOD_PlantHydraulic.o                          \
	  MOD_FrictionVelocity.o                        \
	  MOD_Ozone.o                                   \
	  MOD_LeafTemperature.o                         \
	  MOD_LeafTemperaturePC.o                       \
	  MOD_SoilHcapCond.o                            \
	  MOD_SoilThermalParameters.o                   \
	  MOD_SoilSnowHydrology.o                       \
	  MOD_SnowLayersCombineDivide.o                 \
	  MOD_Meltf.o                                   \
	  MOD_Glacier.o                                 \
	  MOD_Lake.o                                    \
	  MOD_SimpleOcean.o                             \
	  MOD_GroundFluxes.o                            \
	  MOD_GroundTem.o                               \
	  MOD_LAIEmpirical.o                            \
	  MOD_LAIReadin.o                               \
	  MOD_CropReadin.o                              \
	  MOD_NitrifReadin.o                            \
	  MOD_NdepReadin.o                              \
	  MOD_FireReadin.o                              \
	  urban/MOD_Urban_LAIReadin.o                   \
	  MOD_LeafInterception.o                        \
	  MOD_NetSolar.o                                \
	  MOD_WetBulb.o                                 \
	  MOD_RainSnowTemp.o                            \
	  MOD_NewSnow.o                                 \
	  MOD_SnowFraction.o                            \
	  MOD_Thermal.o                                 \
	  MOD_Vars_1DAccFluxes.o                        \
	  ../CaMa/src/MOD_CaMa_Vars.o                   \
	  MOD_HistVector.o                              \
	  MOD_Hist.o                                    \
	  MOD_LightningData.o                           \
	  MOD_OzoneData.o                               \
	  ../CaMa/src/MOD_CaMa_colmCaMa.o               \
	  CoLMDRIVER.o                                  \
	  CoLMMAIN.o                                    \
	  CoLM.o

####################################################################
ifneq (${CaMa},\#define)
# Compile CoLM decoupled without river routing scheme (CaMa-Flood)

${EXECUTABLE} : ${HEADER} \
	${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${BASIC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS_URBAN} ${OBJS}
	\
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${BASIC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS_URBAN} ${OBJS} \
	-o $@ ${LDFLAGS}
	@echo 'making CoLM completed!'

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
	${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${BASIC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS_URBAN} ${OBJS} mk_CaMa
	\
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${BASIC} \
	${OBJS_HYDRO} ${OBJS_BGC} ${OBJS_URBAN} ${OBJS} \
	${CaMa_LIBS} -o $@ ${LDFLAGS}
	@echo 'making CoLM with CaMa Flood Model completed!'

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
