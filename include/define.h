#undef CATCHMENT
#define GRIDBASED
#define	USGS_CLASSIFICATION       
#define LANDONLY

#undef  USE_DEPTH_TO_BEDROCK
#undef  VARIABLY_SATURATED_FLOW
#define Campbell_SOIL_MODEL
#undef  vanGenuchten_Mualem_SOIL_MODEL

#undef  DYN_PHENOLOGY             
#undef	SOILINI                   
#undef	SOIL_REFL_GUESSED         
#define	SOIL_REFL_READ            

#define	CLMDEBUG                  
#define USEMPI
#undef  CaMa_Flood
