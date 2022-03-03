# Makefile for CLM main program

include ../include/Makeoptions

DEF  = $(shell grep -i cama_flood ../include/define.h)
CaMa = $(word 1, ${DEF})

# Name of executable
EXECUTABLE = clm.x

####################################################################
.DEFAULT :

OBJS =  precision.o                 \
	GlobalVars.o                \
	PhysicalConstants.o         \
	LC_Const.o                  \
	ncio.o                      \
	MOD_PFTimeInvars.o          \
	MOD_PCTimeInvars.o          \
	MOD_TimeInvariants.o        \
	timemanager.o               \
	PFT_Const.o                 \
	user_specified_forcing.o    \
	metdata.o                   \
	MOD_PFTimeVars.o            \
	MOD_PCTimeVars.o            \
	MOD_TimeVariables.o         \
	GETMET.o                    \
	MOD_1D_PFTFluxes.o          \
	MOD_1D_PCFluxes.o           \
	MOD_1D_Fluxes.o             \
	MOD_1D_Forcing.o            \
	MOD_2D_Fluxes.o             \
	MOD_2D_Forcing.o            \
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
	CLMDRIVER.o                 \
	CLM.o                       \
	CLMMAIN.o                   \
	eroot.o                     \
	flxwrite.o                  \
	groundfluxes.o              \
	groundtem.o                 \
	LAI_empirical.o             \
	LAI_readin.o                \
	LAI_readin_nc.o             \
	LEAF_interception.o         \
	lpwrite.o                   \
	meltf.o                     \
	netsolar.o                  \
	rain_snow_temp.o            \
	newsnow.o                   \
	orb_coszen.o                \
	qsadv.o                     \
	rd_forcing.o                \
	snowage.o                   \
	snowfraction.o              \
	THERMAL.o                   \
	tridia.o                    \
	vec2xy.o                    \
	wetbulb.o                   \
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
        bgc_NDEP_readin.o


####################################################################

ifneq (${CaMa},\#define)
# Compile CoLM decoupled with river routing scheme (CaMa-Flood)

${EXECUTABLE} : ${OBJS}
	${FF} ${FOPTS} ${OBJS} -o $@ ${LDFLAGS}
	@echo 'making CLM completed!'

$(OBJS):%.o:%.F90
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

else
####################################################################

# The global river model CaMa-Flood (version 3.4.4)
CaMa_DIR = ../CaMa

# CaMa Flood Model modules directories
CaMa_MODS = -I$(CaMa_DIR)/mod \
            -I$(CaMa_DIR)/lib \
            -I$(CaMa_DIR)/src

# CaMa Flood Model libs (static) directories
CaMa_LIBS = $(CaMa_DIR)/src/srcs.a \
            $(CaMa_DIR)/mod/mods.a \
            $(CaMa_DIR)/lib/libs.a

#--------------------------------------------

${EXECUTABLE} : ${OBJS} 
	${FF} ${FOPTS} ${OBJS} ${CaMa_LIBS} -o $@ ${LDFLAGS}
	@echo 'making CLM with CaMa Flood Model completed!'

$(OBJS):%.o:%.F90
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${CaMa_MODS}

endif
####################################################################

clean :
	${RM} -f *.o *.mod
