# Makefile for CoLM initial data processing

include ../include/Makeoptions

INCLUDE_DIR = -I../include -I../share -I../mksrfdata -I../main \
	  -I../main/hydro -I../main/bgc -I../main/urban -I${NETCDF_INC}
VPATH = ../share : ../mksrfdata : ../main : ../main/hydro : ../main/bgc : ../main/urban
HEADER = ../include/define.h

# name of executable
EXECUTABLE = initial.x

####################################################################

OBJS_SHARED = \
				  ../share/precision.o \
				  ../share/spmd_task.o \
				  ../share/mod_namelist.o \
				  ../share/PhysicalConstants.o \
				  ../share/MathConstants.o \
				  ../share/timemanager.o \
				  ../share/mod_utils.o \
				  ../share/ncio_serial.o \
				  ../share/mod_block.o \
				  ../share/mod_grid.o \
				  ../share/mod_pixel.o \
				  ../share/mod_data_type.o \
				  ../share/ncio_block.o \
				  ../share/mod_catchment_data.o \
				  ../share/mod_5x5_data.o \
				  ../share/mod_mesh.o \
				  ../share/mod_pixelset.o \
				  ../share/ncio_vector.o \
				  ../share/mod_colm_debug.o \
				  ../share/GlobalVars.o \
				  ../share/LC_Const.o \
				  ../share/UrbanLCZ_Const.o              \
				  ../share/mod_mapping_grid2pset.o

OBJS_MKSRFDATA = \
				  ../mksrfdata/MOD_AggregationRequestData.o \
				  ../mksrfdata/MOD_SrfdataRestart.o \
				  ../mksrfdata/MOD_LandElm.o \
				  ../mksrfdata/MOD_LandHRU.o \
		 		  ../mksrfdata/MOD_PixelsetShadow.o \
				  ../mksrfdata/MOD_LandPatch.o \
				  ../mksrfdata/MOD_LandPFT.o \
				  ../mksrfdata/MOD_LandPC.o \
				  ../mksrfdata/MOD_SingleSrfdata.o \
				  ../mksrfdata/MOD_ElmVector.o \
				  ../mksrfdata/MOD_HRUVector.o \
				  ../mksrfdata/MOD_LandUrban.o

BASIC_HYDRO = \
				  ../main/hydro/mod_io_basin.o \
				  ../main/hydro/MOD_HydroTimeVars.o \
				  ../main/hydro/MOD_1D_HydroFluxes.o \
				  ../main/hydro/mod_soil_function.o \
				  ../main/hydro/mod_soil_water.o

BASIC_BGC = \
			  ../main/bgc/MOD_BGC_Vars_1DFluxes.o  \
			  ../main/bgc/MOD_BGC_Vars_1DPFTFluxes.o  \
			  ../main/bgc/MOD_BGC_Vars_2DFluxes.o  \
			  ../main/bgc/MOD_BGC_Vars_PFTimeVars.o  \
			  ../main/bgc/MOD_BGC_Vars_TimeInvars.o  \
			  ../main/bgc/MOD_BGC_Vars_TimeVars.o 

BASIC_URBAN = \
	                  ../main/urban/MOD_1D_UrbanFluxes.o     \
                 	  ../main/urban/MOD_UrbanTimeVars.o      \
                 	  ../main/urban/MOD_UrbanTimeInvars.o    

BASIC =  \
		  ../main/MOD_Vars_PFTConst.o \
		  ../main/MOD_Vars_PFTimeVars.o \
		  ../main/MOD_Vars_PFTimeInvars.o \
		  ../main/MOD_Vars_PCTimeVars.o \
		  ../main/MOD_Vars_PCTimeInvars.o \
		  ../main/MOD_Vars_TimeInvariants.o        \
		  ../main/MOD_Vars_TimeVariables.o    \
		  ../main/MOD_Vars_1DPFTFluxes.o          \
		  ../main/MOD_Vars_1DPCFluxes.o           \
		  ../main/MOD_Vars_1DFluxes.o \
		  ../main/MOD_Vars_1DForcing.o 

OBJS_MAIN = \
				  ../main/MOD_LAIEmpirical.o  \
				  ../main/MOD_LAIReadin.o \
				  ../main/MOD_OrbCoszen.o \
				  ../main/MOD_ThreeDCanopy.o \
				  ../main/MOD_SnowSnicar.o \
				  ../main/MOD_SnowAge.o \
				  ../main/MOD_SnowAlbedo.o \
				  ../main/MOD_Aerosol.o    \
				  ../main/MOD_Albedo.o    \
				  ../main/MOD_SnowFraction.o \
				  ../main/urban/UrbanLAI_readin.o        \
  				  ../main/urban/UrbanShortwave.o         \
				  ../main/urban/UrbanALBEDO.o

OBJS = \
		 MOD_PercentagesPFTReadin.o       \
		 MOD_LakeDepthReadin.o       \
		 MOD_DBedrockReadin.o       \
		 MOD_SoilColorRefl.o        \
		 MOD_SoilParametersReadin.o \
		 MOD_HtopReadin.o \
		 MOD_UrbanReadin.o \
		 ../main/MOD_CropReadin.o \
		 ../main/MOD_NitrifReadin.o \
	         ../main/MOD_NdepReadin.o   \
		 ../main/MOD_FireReadin.o   \
		 MOD_IniTimeVar.o             \
		 MOD_IniTimeVarUrban.o        \
		 MOD_Initialize.o             \
		 CoLMINI.o 

####################################################################

all: ${EXECUTABLE}

${EXECUTABLE} : ${HEADER} ${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC} ${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${OBJS_MAIN} ${OBJS}
	\
	${FF} ${FOPTS} ${OBJS_SHARED} ${OBJS_MKSRFDATA} \
	${BASIC} ${BASIC_BGC} ${BASIC_HYDRO} ${BASIC_URBAN} ${OBJS_MAIN} ${OBJS} -o $@ ${LDFLAGS}
	@echo 'making CoLM initial data completed!'

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
