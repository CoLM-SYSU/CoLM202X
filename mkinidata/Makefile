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
				  ../share/MOD_Precision.o           \
				  ../share/MOD_SPMD_Task.o           \
				  ../share/MOD_Namelist.o            \
				  ../share/PhysicalConstants.o       \
				  ../share/MathConstants.o           \
				  ../share/MOD_TimeManager.o         \
				  ../share/MOD_Utils.o               \
				  ../share/MOD_NetCDFSerial.o        \
				  ../share/MOD_Block.o               \
				  ../share/MOD_Grid.o                \
				  ../share/MOD_Pixel.o               \
				  ../share/MOD_DataType.o            \
				  ../share/MOD_NetCDFBlock.o         \
				  ../share/MOD_CatchmentDataReadin.o \
				  ../share/MOD_5x5DataReadin.o       \
				  ../share/MOD_Mesh.o                \
				  ../share/MOD_Pixelset.o            \
				  ../share/MOD_NetCDFVector.o        \
				  ../share/MOD_CoLMDebug.o           \
				  ../share/GlobalVars.o              \
				  ../share/LC_Const.o                \
				  ../share/UrbanLCZ_Const.o          \
				  ../share/MOD_Mapping_Grid2Pset.o

OBJS_MKSRFDATA = \
				  ../mksrfdata/mod_aggregation.o \
				  ../mksrfdata/mod_srfdata_restart.o \
				  ../mksrfdata/mod_landelm.o \
				  ../mksrfdata/mod_landhru.o \
		 		  ../mksrfdata/mod_pixelsetshadow.o \
				  ../mksrfdata/mod_landpatch.o \
				  ../mksrfdata/mod_landpft.o \
				  ../mksrfdata/mod_landpc.o \
				  ../mksrfdata/mod_single_srfdata.o \
				  ../mksrfdata/mod_elm_vector.o \
				  ../mksrfdata/mod_hru_vector.o \
				  ../mksrfdata/mod_landurban.o

BASIC_HYDRO = \
				  ../main/hydro/MOD_Hydro_IO.o \
				  ../main/hydro/MOD_Hydro_Vars_Timevars.o \
				  ../main/hydro/MOD_Hydro_Vars_1DFluxes.o \
				  ../main/hydro/MOD_Hydro_SoilFunction.o \
				  ../main/hydro/MOD_Hydro_SoilWater.o

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
		 pct_readin.o       \
		 lakedepth_readin.o       \
		 dbedrock_readin.o       \
		 soil_color_refl.o        \
		 soil_parameters_readin.o \
		 HTOP_readin.o \
		 Urban_readin.o \
		 ../main/MOD_CropReadin.o \
		 ../main/MOD_NitrifReadin.o \
	         ../main/MOD_NdepReadin.o   \
		 ../main/MOD_FireReadin.o   \
		 IniTimeVar.o             \
		 UrbanIniTimeVar.o        \
		 initialize.o             \
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
