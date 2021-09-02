# Makefile for CLM surface data processing

include ../include/Makeoptions

# name of executable
EXECUTABLE = srf.x

####################################################################
.DEFAULT :

OBJS =  precision.o \
	GlobalVars.o \
	mksrfdata.o \
	rdgrid.o \
	crgrid.o \
	cellarea.o \
	info_gridcell.o \
	rd_land_types.o \
	rd_soil_properties.o \
	aggregation_forest_height.o \
	aggregation_LAI.o \
	aggregation_lakedepth.o \
	aggregation_landtypes.o \
	aggregation_soil_parameters.o \
	aggregation_soil_brightness.o \
	aggregation_landwaters.o \
	aggregation_wetland.o \
	soil_hydraulic_parameters.o \
	soil_thermal_parameters.o \
	array_sorting_mod.o \
	median.o

####################################################################

${EXECUTABLE} : ${OBJS} 
	${FF} ${FOPTS} ${OBJS} -o $@ ${LDFLAGS}
	@echo 'making CLM surface data completed!'

$(OBJS):%.o:%.F90
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

clean :
	${RM} -f *.o *.mod
