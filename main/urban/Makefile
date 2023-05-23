# Makefile for CLM surface data processing

include ../../include/Makeoptions

INCLUDE_DIR = -I../../include -I../../share -I../../mksrfdata -I.. \
				  -I${NETCDF_INC}

VPATH  = ../../share  ../../mksrfdata ..
HEADER = ../../include/define.h

####################################################################
.DEFAULT :

OBJS_SHARED =                              \
	  ../../share/precision.o          \
	  ../../share/MathConstants.o      \
	  ../MOD_Vars_PhysicalConst.o      \
	  ../MOD_Vars_Global.o             \
	  ../MOD_Vars_LCConst.o            \
	  ../../share/spmd_task.o          \
	  ../../share/mod_utils.o          \
	  ../../share/timemanager.o        \
	  ../../share/ncio_serial.o        \
	  ../../share/mod_namelist.o       \
	  ../../share/mod_block.o          \
	  ../../share/mod_grid.o           \
	  ../../share/mod_pixel.o          \
	  ../../share/mod_data_type.o      \
	  ../../share/ncio_block.o         \
	  ../../share/mod_catchment_data.o \
	  ../../share/mod_pixelset.o       \
	  ../../share/ncio_vector.o        \
	  ../../share/mod_colm_debug.o

BASIC =                                    \
	  ../MOD_Vars_TimeInvariants.o     \
	  ../MOD_Vars_TimeVariables.o      \
	  ../MOD_Vars_1DPFTFluxes.o        \
	  ../MOD_Vars_1DPCFluxes.o         \
	  ../MOD_Vars_1DFluxes.o

OBJS =                                     \
	  MOD_Urban_Vars_LCZConst.o        \
	  MOD_Urban_Vars_1DFluxes.o        \
	  MOD_Urban_Vars_TimeInvars.o      \
	  MOD_Urban_Vars_TimeVars.o        \
	  MOD_Urban_Shortwave.o            \
	  MOD_Urban_Longwave.o             \
	  MOD_Urban_Albedo.o               \
	  MOD_Urban_NetSolar.o             \
	  MOD_Urban_Flux.o                 \
	  MOD_Urban_GroundFlux.o           \
	  MOD_Urban_RoofFlux.o             \
	  MOD_Urban_RoofTem.o              \
	  MOD_Urban_WallTem.o              \
	  MOD_Urban_PerviousTem.o          \
	  MOD_Urban_ImperviousTem.o        \
	  MOD_Urban_Hydrology.o            \
	  MOD_Urban_BEM.o                  \
	  MOD_Urban_LUCY.o                 \
	  MOD_Urban_Thermal.o              \
	  Urban_CoLMMAIN.o


####################################################################

main : ${OBJS} ${OBJS_SHARED} ${BASIC}
	@echo 'Compiling modules in urban completed!'

$(OBJS):%.o : %.F90 ${HEADER} ${OBJS_SHARED} ${BASIC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $<

$(OBJS_SHARED) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../../share

$(BASIC) : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< ${MOD_CMD}../

clean :
	${RM} -f *.o *.mod
