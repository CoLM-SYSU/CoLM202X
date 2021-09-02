# Makefile for CLM initial data processing

include ../include/Makeoptions

# name of executable
EXECUTABLE = initial.x

####################################################################

OBJS =  precision.o              \
	GlobalVars.o             \
	LC_Const.o               \
	PFT_Const.o              \
	ncio.o                   \
	MOD_PFTimeInvars.o       \
	MOD_PCTimeInvars.o       \
	MOD_TimeInvariants.o     \
	timemanager.o            \
	PhysicalConstants.o      \
	MOD_PFTimeVars.o         \
	MOD_PCTimeVars.o         \
	MOD_TimeVariables.o      \
	CLMINI.o                 \
	initialize.o             \
	ThreeDCanopy.o           \
	ALBEDO.o                 \
	IniTimeVar.o             \
	LAI_empirical.o          \
	LAI_readin.o             \
	LAI_readin_nc.o          \
	HTOP_readin_nc.o         \
	lakedepth_readin.o       \
	orb_coszen.o             \
	snowfraction.o           \
	soil_color_refl.o        \
	soil_parameters_readin.o

####################################################################

all: ${EXECUTABLE}

${EXECUTABLE} : ${OBJS} 
	${FF} ${FOPTS} ${OBJS} -o $@ ${LDFLAGS}
	@echo 'making CLM initial data completed!'

$(OBJS):%.o:%.F90
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

clean :
	${RM} -f *.o *.mod
