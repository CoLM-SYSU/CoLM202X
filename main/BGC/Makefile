# Makefile for CoLM surface data processing

include ../../include/Makeoptions
 
INCLUDE_DIR = -I../../include -I../../share -I../../mksrfdata \
				  -I.. -I${NETCDF_INC}

VPATH  = ../../share ../../mksrfdata ..
HEADER = ../../include/define.h

####################################################################
.DEFAULT :

OBJS_SHARED = \
				  ../../share/precision.o \
				  ../../share/GlobalVars.o \
				  ../../share/LC_Const.o \
				  ../../share/spmd_task.o          \
				  ../../share/mod_namelist.o \
				  ../../share/PhysicalConstants.o  \
				  ../../share/MathConstants.o  \
				  ../../share/mod_utils.o          \
				  ../../share/timemanager.o        \
				  ../../share/ncio_serial.o \
				  ../../share/mod_block.o          \
				  ../../share/mod_grid.o           \
				  ../../share/mod_pixel.o          \
				  ../../share/mod_data_type.o      \
				  ../../share/ncio_block.o  \
				  ../../share/mod_hydro_data.o \
				  ../../share/mod_landunit.o   \
				  ../../share/mod_pixelset.o          \
				  ../../share/ncio_vector.o \
				  ../../share/mod_colm_debug.o

OBJS_MKSRFDATA = \
				  ../../mksrfdata/mod_modis_data.o \
				  ../../mksrfdata/MOD_PixelsetShadow.o \
				  ../../mksrfdata/mod_landcell.o \
				  ../../mksrfdata/MOD_LandPatch.o \
				  ../../mksrfdata/MOD_LandPFT.o

OBJS_MAIN =  \
		  ../MOD_Const_PFT.o \
		  ../MOD_Vars_PFTimeVars.o \
		  ../MOD_Vars_PFTimeInvars.o \
		  ../MOD_Vars_PCTimeVars.o \
		  ../MOD_Vars_PCTimeInvars.o \
		  ../MOD_Vars_TimeInvariants.o        \
		  ../MOD_Vars_TimeVariables.o         \
		  ../MOD_Vars_1DPFTFluxes.o          \
		  ../MOD_Vars_1DPCFluxes.o           \
		  ../MOD_Vars_1DFluxes.o             \
		  ../MOD_Vars_1DForcing.o            \
		  ../MOD_Vars_2DFluxes.o             \
		  ../MOD_Vars_2DForcing.o

OBJS = \
		 MOD_1D_BGCFluxes.o  \
		 MOD_1D_BGCPFTFluxes.o  \
		 MOD_2D_BGCFluxes.o  \
		 MOD_BGCPFTimeVars.o  \
		 MOD_BGCTimeInvars.o  \
		 MOD_BGCTimeVars.o    \
		 bgc_CNCStateUpdate1Mod.o    \
		 bgc_CNCStateUpdate2Mod.o    \
		 bgc_CNCStateUpdate3Mod.o    \
		 bgc_CNNStateUpdate1Mod.o    \
		 bgc_CNNStateUpdate2Mod.o    \
		 bgc_CNNStateUpdate3Mod.o    \
		 bgc_SoilBiogeochemNStateUpdate1Mod.o \
		 bgc_soil_SoilBiogeochemCompetitionMod.o  \
		 bgc_soil_SoilBiogeochemDecompCascadeBGCMod.o  \
		 bgc_soil_SoilBiogeochemDecompMod.o  \
		 bgc_soil_SoilBiogeochemLittVertTranspMod.o  \
		 bgc_soil_SoilBiogeochemNLeachingMod.o  \
		 bgc_soil_SoilBiogeochemPotentialMod.o  \
		 bgc_soil_SoilBiogeochemVerticalProfileMod.o  \
		 bgc_veg_CNGapMortalityMod.o  \
		 bgc_veg_CNGRespMod.o  \
		 bgc_veg_CNMRespMod.o  \
		 bgc_DaylengthMod.o \
		 bgc_veg_CNPhenologyMod.o  \
		 bgc_veg_NutrientCompetitionMod.o  \
		 bgc_veg_CNVegStructUpdateMod.o \
		 bgc_CNSummaryMod.o  \
		 bgc_CNAnnualUpdateMod.o \
		 bgc_CNZeroFluxesMod.o \
		 bgc_CNBalanceCheckMod.o \
		 bgc_CNSASUMod.o \
		 bgc_driver.o \

####################################################################

main : ${OBJS}
	@echo 'Compiling modules in BGC completed!'

$(OBJS):%.o : %.F90 ${HEADER} ${OBJS_SHARED} ${OBJS_MKSRFDATA} ${OBJS_MAIN}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

$(OBJS_SHARED) : %.o : %.F90 ${HEADER} 
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../../share

$(OBJS_MKSRFDATA) : %.o : %.F90 ${HEADER} ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../../mksrfdata

$(OBJS_MAIN) : %.o : %.F90 ${HEADER} ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}..

clean :
	${RM} -f *.o *.mod
