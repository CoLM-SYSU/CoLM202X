#include <define.h>

MODULE MOD_2D_Fluxes
! ----------------------------------------------------------------------
! perfrom the grid average mapping: average a subgrid input 1d vector 
! of length numpatch to a output 2d array of length [lon_points,lat_points]
!
! Created by Yongjiu Dai, 03/2014
!---------------------------------------------------------------------

   use mod_data_type
   USE GlobalVars
#ifdef BGC
   USE MOD_2D_BGCFluxes
#endif

   IMPLICIT NONE
   SAVE

   type(block_data_real8_2d) :: f_taux    ! wind stress: E-W [kg/m/s2]
   type(block_data_real8_2d) :: f_tauy    ! wind stress: N-S [kg/m/s2]
   type(block_data_real8_2d) :: f_fsena   ! sensible heat from canopy height to atmosphere [W/m2]
   type(block_data_real8_2d) :: f_lfevpa  ! latent heat flux from canopy height to atmosphere [W/m2]
   type(block_data_real8_2d) :: f_fevpa   ! evapotranspiration from canopy to atmosphere [mm/s]
   type(block_data_real8_2d) :: f_fsenl   ! sensible heat from leaves [W/m2]
   type(block_data_real8_2d) :: f_fevpl   ! evaporation+transpiration from leaves [mm/s]
   type(block_data_real8_2d) :: f_etr     ! transpiration rate [mm/s]
   type(block_data_real8_2d) :: f_fseng   ! sensible heat flux from ground [W/m2]
   type(block_data_real8_2d) :: f_fevpg   ! evaporation heat flux from ground [mm/s]
   type(block_data_real8_2d) :: f_fgrnd   ! ground heat flux [W/m2]
   type(block_data_real8_2d) :: f_sabvsun ! solar absorbed by sunlit canopy [W/m2]
   type(block_data_real8_2d) :: f_sabvsha ! solar absorbed by shaded [W/m2]
   type(block_data_real8_2d) :: f_sabg    ! solar absorbed by ground  [W/m2]
   type(block_data_real8_2d) :: f_sr      ! total reflected solar radiation (W/m2)
   type(block_data_real8_2d) :: f_solvd   ! incident direct beam vis solar radiation (W/m2)
   type(block_data_real8_2d) :: f_solvi   ! incident diffuse beam vis solar radiation (W/m2)
   type(block_data_real8_2d) :: f_solnd   ! incident direct beam nir solar radiation (W/m2)
   type(block_data_real8_2d) :: f_solni   ! incident diffuse beam nir solar radiation (W/m2)
   type(block_data_real8_2d) :: f_srvd    ! reflected direct beam vis solar radiation (W/m2)
   type(block_data_real8_2d) :: f_srvi    ! reflected diffuse beam vis solar radiation (W/m2)
   type(block_data_real8_2d) :: f_srnd    ! reflected direct beam nir solar radiation (W/m2)
   type(block_data_real8_2d) :: f_srni    ! reflected diffuse beam nir solar radiation (W/m2)
   type(block_data_real8_2d) :: f_solvdln ! incident direct beam vis solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_solviln ! incident diffuse beam vis solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_solndln ! incident direct beam nir solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_solniln ! incident diffuse beam nir solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_srvdln  ! reflected direct beam vis solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_srviln  ! reflected diffuse beam vis solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_srndln  ! reflected direct beam nir solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_srniln  ! reflected diffuse beam nir solar radiation at local noon (W/m2)
   type(block_data_real8_2d) :: f_olrg    ! outgoing long-wave radiation from ground+canopy [W/m2]
   type(block_data_real8_2d) :: f_rnet    ! net radiation [W/m2]
   type(block_data_real8_2d) :: f_xerr    ! the error of water banace [mm/s]
   type(block_data_real8_2d) :: f_zerr    ! the error of energy balance [W/m2]
   type(block_data_real8_2d) :: f_rsur    ! surface runoff [mm/s]
   type(block_data_real8_2d) :: f_rnof    ! total runoff [mm/s]
   type(block_data_real8_2d) :: f_qintr   ! interception [mm/s]
   type(block_data_real8_2d) :: f_qinfl   ! inflitration [mm/s]
   type(block_data_real8_2d) :: f_qdrip   ! throughfall [mm/s]
   type(block_data_real8_2d) :: f_assim   ! canopy assimilation rate [mol m-2 s-1]
   type(block_data_real8_2d) :: f_respc   ! respiration (plant+soil) [mol m-2 s-1]
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
   type(block_data_real8_2d) :: f_assim_RuBP_sun_enftemp        
   type(block_data_real8_2d) :: f_assim_RuBP_sun_enfboreal      
   type(block_data_real8_2d) :: f_assim_RuBP_sun_dnfboreal      
   type(block_data_real8_2d) :: f_assim_RuBP_sun_ebftrop        
   type(block_data_real8_2d) :: f_assim_RuBP_sun_ebftemp        
   type(block_data_real8_2d) :: f_assim_RuBP_sun_dbftrop        
   type(block_data_real8_2d) :: f_assim_RuBP_sun_dbftemp        
   type(block_data_real8_2d) :: f_assim_RuBP_sun_dbfboreal      
   type(block_data_real8_2d) :: f_assim_RuBP_sun_ebstemp        
   type(block_data_real8_2d) :: f_assim_RuBP_sun_dbstemp        
   type(block_data_real8_2d) :: f_assim_RuBP_sun_dbsboreal      
   type(block_data_real8_2d) :: f_assim_RuBP_sun_c3arcgrass     
   type(block_data_real8_2d) :: f_assim_RuBP_sun_c3grass        
   type(block_data_real8_2d) :: f_assim_RuBP_sun_c4grass        
   type(block_data_real8_2d) :: f_assim_RuBP_sha_enftemp        
   type(block_data_real8_2d) :: f_assim_RuBP_sha_enfboreal      
   type(block_data_real8_2d) :: f_assim_RuBP_sha_dnfboreal      
   type(block_data_real8_2d) :: f_assim_RuBP_sha_ebftrop        
   type(block_data_real8_2d) :: f_assim_RuBP_sha_ebftemp        
   type(block_data_real8_2d) :: f_assim_RuBP_sha_dbftrop        
   type(block_data_real8_2d) :: f_assim_RuBP_sha_dbftemp        
   type(block_data_real8_2d) :: f_assim_RuBP_sha_dbfboreal      
   type(block_data_real8_2d) :: f_assim_RuBP_sha_ebstemp        
   type(block_data_real8_2d) :: f_assim_RuBP_sha_dbstemp        
   type(block_data_real8_2d) :: f_assim_RuBP_sha_dbsboreal      
   type(block_data_real8_2d) :: f_assim_RuBP_sha_c3arcgrass     
   type(block_data_real8_2d) :: f_assim_RuBP_sha_c3grass        
   type(block_data_real8_2d) :: f_assim_RuBP_sha_c4grass        
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_enftemp        
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_enfboreal      
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_dnfboreal      
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_ebftrop        
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_ebftemp        
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_dbftrop        
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_dbftemp        
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_dbfboreal      
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_ebstemp        
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_dbstemp        
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_dbsboreal      
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_c3arcgrass     
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_c3grass        
   type(block_data_real8_2d) :: f_assim_Rubisco_sun_c4grass        
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_enftemp        
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_enfboreal      
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_dnfboreal      
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_ebftrop        
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_ebftemp        
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_dbftrop        
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_dbftemp        
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_dbfboreal      
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_ebstemp        
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_dbstemp        
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_dbsboreal      
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_c3arcgrass     
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_c3grass        
   type(block_data_real8_2d) :: f_assim_Rubisco_sha_c4grass        
   type(block_data_real8_2d) :: f_assimsun_enftemp        
   type(block_data_real8_2d) :: f_assimsun_enfboreal      
   type(block_data_real8_2d) :: f_assimsun_dnfboreal      
   type(block_data_real8_2d) :: f_assimsun_ebftrop        
   type(block_data_real8_2d) :: f_assimsun_ebftemp        
   type(block_data_real8_2d) :: f_assimsun_dbftrop        
   type(block_data_real8_2d) :: f_assimsun_dbftemp        
   type(block_data_real8_2d) :: f_assimsun_dbfboreal      
   type(block_data_real8_2d) :: f_assimsun_ebstemp        
   type(block_data_real8_2d) :: f_assimsun_dbstemp        
   type(block_data_real8_2d) :: f_assimsun_dbsboreal      
   type(block_data_real8_2d) :: f_assimsun_c3arcgrass     
   type(block_data_real8_2d) :: f_assimsun_c3grass        
   type(block_data_real8_2d) :: f_assimsun_c4grass        
   type(block_data_real8_2d) :: f_assimsha_enftemp        
   type(block_data_real8_2d) :: f_assimsha_enfboreal      
   type(block_data_real8_2d) :: f_assimsha_dnfboreal      
   type(block_data_real8_2d) :: f_assimsha_ebftrop        
   type(block_data_real8_2d) :: f_assimsha_ebftemp        
   type(block_data_real8_2d) :: f_assimsha_dbftrop        
   type(block_data_real8_2d) :: f_assimsha_dbftemp        
   type(block_data_real8_2d) :: f_assimsha_dbfboreal      
   type(block_data_real8_2d) :: f_assimsha_ebstemp        
   type(block_data_real8_2d) :: f_assimsha_dbstemp        
   type(block_data_real8_2d) :: f_assimsha_dbsboreal      
   type(block_data_real8_2d) :: f_assimsha_c3arcgrass     
   type(block_data_real8_2d) :: f_assimsha_c3grass        
   type(block_data_real8_2d) :: f_assimsha_c4grass        
   type(block_data_real8_2d) :: f_etrsun_enftemp        
   type(block_data_real8_2d) :: f_etrsun_enfboreal      
   type(block_data_real8_2d) :: f_etrsun_dnfboreal      
   type(block_data_real8_2d) :: f_etrsun_ebftrop        
   type(block_data_real8_2d) :: f_etrsun_ebftemp        
   type(block_data_real8_2d) :: f_etrsun_dbftrop        
   type(block_data_real8_2d) :: f_etrsun_dbftemp        
   type(block_data_real8_2d) :: f_etrsun_dbfboreal      
   type(block_data_real8_2d) :: f_etrsun_ebstemp        
   type(block_data_real8_2d) :: f_etrsun_dbstemp        
   type(block_data_real8_2d) :: f_etrsun_dbsboreal      
   type(block_data_real8_2d) :: f_etrsun_c3arcgrass     
   type(block_data_real8_2d) :: f_etrsun_c3grass        
   type(block_data_real8_2d) :: f_etrsun_c4grass        
   type(block_data_real8_2d) :: f_etrsha_enftemp        
   type(block_data_real8_2d) :: f_etrsha_enfboreal      
   type(block_data_real8_2d) :: f_etrsha_dnfboreal      
   type(block_data_real8_2d) :: f_etrsha_ebftrop        
   type(block_data_real8_2d) :: f_etrsha_ebftemp        
   type(block_data_real8_2d) :: f_etrsha_dbftrop        
   type(block_data_real8_2d) :: f_etrsha_dbftemp        
   type(block_data_real8_2d) :: f_etrsha_dbfboreal      
   type(block_data_real8_2d) :: f_etrsha_ebstemp        
   type(block_data_real8_2d) :: f_etrsha_dbstemp        
   type(block_data_real8_2d) :: f_etrsha_dbsboreal      
   type(block_data_real8_2d) :: f_etrsha_c3arcgrass     
   type(block_data_real8_2d) :: f_etrsha_c3grass        
   type(block_data_real8_2d) :: f_etrsha_c4grass        
   type(block_data_real8_2d) :: f_cisun_enftemp        
   type(block_data_real8_2d) :: f_cisun_enfboreal      
   type(block_data_real8_2d) :: f_cisun_dnfboreal      
   type(block_data_real8_2d) :: f_cisun_ebftrop        
   type(block_data_real8_2d) :: f_cisun_ebftemp        
   type(block_data_real8_2d) :: f_cisun_dbftrop        
   type(block_data_real8_2d) :: f_cisun_dbftemp        
   type(block_data_real8_2d) :: f_cisun_dbfboreal      
   type(block_data_real8_2d) :: f_cisun_ebstemp        
   type(block_data_real8_2d) :: f_cisun_dbstemp        
   type(block_data_real8_2d) :: f_cisun_dbsboreal      
   type(block_data_real8_2d) :: f_cisun_c3arcgrass     
   type(block_data_real8_2d) :: f_cisun_c3grass        
   type(block_data_real8_2d) :: f_cisun_c4grass        
   type(block_data_real8_2d) :: f_cisha_enftemp        
   type(block_data_real8_2d) :: f_cisha_enfboreal      
   type(block_data_real8_2d) :: f_cisha_dnfboreal      
   type(block_data_real8_2d) :: f_cisha_ebftrop        
   type(block_data_real8_2d) :: f_cisha_ebftemp        
   type(block_data_real8_2d) :: f_cisha_dbftrop        
   type(block_data_real8_2d) :: f_cisha_dbftemp        
   type(block_data_real8_2d) :: f_cisha_dbfboreal      
   type(block_data_real8_2d) :: f_cisha_ebstemp        
   type(block_data_real8_2d) :: f_cisha_dbstemp        
   type(block_data_real8_2d) :: f_cisha_dbsboreal      
   type(block_data_real8_2d) :: f_cisha_c3arcgrass     
   type(block_data_real8_2d) :: f_cisha_c3grass        
   type(block_data_real8_2d) :: f_cisha_c4grass        
   type(block_data_real8_2d) :: f_essun_enftemp        
   type(block_data_real8_2d) :: f_essun_enfboreal      
   type(block_data_real8_2d) :: f_essun_dnfboreal      
   type(block_data_real8_2d) :: f_essun_ebftrop        
   type(block_data_real8_2d) :: f_essun_ebftemp        
   type(block_data_real8_2d) :: f_essun_dbftrop        
   type(block_data_real8_2d) :: f_essun_dbftemp        
   type(block_data_real8_2d) :: f_essun_dbfboreal      
   type(block_data_real8_2d) :: f_essun_ebstemp        
   type(block_data_real8_2d) :: f_essun_dbstemp        
   type(block_data_real8_2d) :: f_essun_dbsboreal      
   type(block_data_real8_2d) :: f_essun_c3arcgrass     
   type(block_data_real8_2d) :: f_essun_c3grass        
   type(block_data_real8_2d) :: f_essun_c4grass        
   type(block_data_real8_2d) :: f_essha_enftemp        
   type(block_data_real8_2d) :: f_essha_enfboreal      
   type(block_data_real8_2d) :: f_essha_dnfboreal      
   type(block_data_real8_2d) :: f_essha_ebftrop        
   type(block_data_real8_2d) :: f_essha_ebftemp        
   type(block_data_real8_2d) :: f_essha_dbftrop        
   type(block_data_real8_2d) :: f_essha_dbftemp        
   type(block_data_real8_2d) :: f_essha_dbfboreal      
   type(block_data_real8_2d) :: f_essha_ebstemp        
   type(block_data_real8_2d) :: f_essha_dbstemp        
   type(block_data_real8_2d) :: f_essha_dbsboreal      
   type(block_data_real8_2d) :: f_essha_c3arcgrass     
   type(block_data_real8_2d) :: f_essha_c3grass        
   type(block_data_real8_2d) :: f_essha_c4grass        
   type(block_data_real8_2d) :: f_gssun_enftemp        
   type(block_data_real8_2d) :: f_gssun_enfboreal      
   type(block_data_real8_2d) :: f_gssun_dnfboreal      
   type(block_data_real8_2d) :: f_gssun_ebftrop        
   type(block_data_real8_2d) :: f_gssun_ebftemp        
   type(block_data_real8_2d) :: f_gssun_dbftrop        
   type(block_data_real8_2d) :: f_gssun_dbftemp        
   type(block_data_real8_2d) :: f_gssun_dbfboreal      
   type(block_data_real8_2d) :: f_gssun_ebstemp        
   type(block_data_real8_2d) :: f_gssun_dbstemp        
   type(block_data_real8_2d) :: f_gssun_dbsboreal      
   type(block_data_real8_2d) :: f_gssun_c3arcgrass     
   type(block_data_real8_2d) :: f_gssun_c3grass        
   type(block_data_real8_2d) :: f_gssun_c4grass        
   type(block_data_real8_2d) :: f_gssha_enftemp        
   type(block_data_real8_2d) :: f_gssha_enfboreal      
   type(block_data_real8_2d) :: f_gssha_dnfboreal      
   type(block_data_real8_2d) :: f_gssha_ebftrop        
   type(block_data_real8_2d) :: f_gssha_ebftemp        
   type(block_data_real8_2d) :: f_gssha_dbftrop        
   type(block_data_real8_2d) :: f_gssha_dbftemp        
   type(block_data_real8_2d) :: f_gssha_dbfboreal      
   type(block_data_real8_2d) :: f_gssha_ebstemp        
   type(block_data_real8_2d) :: f_gssha_dbstemp        
   type(block_data_real8_2d) :: f_gssha_dbsboreal      
   type(block_data_real8_2d) :: f_gssha_c3arcgrass     
   type(block_data_real8_2d) :: f_gssha_c3grass        
   type(block_data_real8_2d) :: f_gssha_c4grass        
   type(block_data_real8_2d) :: f_gammasun_enftemp        
   type(block_data_real8_2d) :: f_gammasun_enfboreal      
   type(block_data_real8_2d) :: f_gammasun_dnfboreal      
   type(block_data_real8_2d) :: f_gammasun_ebftrop        
   type(block_data_real8_2d) :: f_gammasun_ebftemp        
   type(block_data_real8_2d) :: f_gammasun_dbftrop        
   type(block_data_real8_2d) :: f_gammasun_dbftemp        
   type(block_data_real8_2d) :: f_gammasun_dbfboreal      
   type(block_data_real8_2d) :: f_gammasun_ebstemp        
   type(block_data_real8_2d) :: f_gammasun_dbstemp        
   type(block_data_real8_2d) :: f_gammasun_dbsboreal      
   type(block_data_real8_2d) :: f_gammasun_c3arcgrass     
   type(block_data_real8_2d) :: f_gammasun_c3grass        
   type(block_data_real8_2d) :: f_gammasun_c4grass        
   type(block_data_real8_2d) :: f_gammasha_enftemp        
   type(block_data_real8_2d) :: f_gammasha_enfboreal      
   type(block_data_real8_2d) :: f_gammasha_dnfboreal      
   type(block_data_real8_2d) :: f_gammasha_ebftrop        
   type(block_data_real8_2d) :: f_gammasha_ebftemp        
   type(block_data_real8_2d) :: f_gammasha_dbftrop        
   type(block_data_real8_2d) :: f_gammasha_dbftemp        
   type(block_data_real8_2d) :: f_gammasha_dbfboreal      
   type(block_data_real8_2d) :: f_gammasha_ebstemp        
   type(block_data_real8_2d) :: f_gammasha_dbstemp        
   type(block_data_real8_2d) :: f_gammasha_dbsboreal      
   type(block_data_real8_2d) :: f_gammasha_c3arcgrass     
   type(block_data_real8_2d) :: f_gammasha_c3grass        
   type(block_data_real8_2d) :: f_gammasha_c4grass        
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_enftemp      
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_enfboreal   
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_dnfboreal  
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_ebftrop   
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_ebftemp    
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_dbftrop     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_dbftemp     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_dbfboreal   
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_ebstemp     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_dbstemp     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_dbsboreal   
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_c3arcgrass  
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_c3grass     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sun_c4grass     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_enftemp     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_enfboreal   
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_dnfboreal   
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_ebftrop     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_ebftemp     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_dbftrop     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_dbftemp     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_dbfboreal   
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_ebstemp     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_dbstemp     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_dbsboreal   
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_c3arcgrass  
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_c3grass     
   type(block_data_real8_2d) :: f_RuBPlimitfrac_sha_c4grass     
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_enftemp  
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_enfboreal 
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_dnfboreal 
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_ebftrop   
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_ebftemp   
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_dbftrop   
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_dbftemp   
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_dbfboreal 
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_ebstemp   
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_dbstemp   
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_dbsboreal 
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_c3arcgrass 
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_c3grass    
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sun_c4grass    
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_enftemp    
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_enfboreal  
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_dnfboreal  
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_ebftrop    
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_ebftemp    
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_dbftrop    
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_dbftemp    
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_dbfboreal  
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_ebstemp    
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_dbstemp    
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_dbsboreal  
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_c3arcgrass 
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_c3grass    
   type(block_data_real8_2d) :: f_Rubiscolimitfrac_sha_c4grass    
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_enftemp       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_enfboreal     
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_dnfboreal     
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_ebftrop       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_ebftemp       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_dbftrop       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_dbftemp       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_dbfboreal     
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_ebstemp       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_dbstemp       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_dbsboreal     
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_c3arcgrass    
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_c3grass       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sun_c4grass       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_enftemp       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_enfboreal     
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_dnfboreal     
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_ebftrop       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_ebftemp       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_dbftrop       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_dbftemp       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_dbfboreal     
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_ebstemp       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_dbstemp       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_dbsboreal     
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_c3arcgrass    
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_c3grass       
   type(block_data_real8_2d) :: f_Sinklimitfrac_sha_c4grass       
   type(block_data_real8_2d) :: f_rstfacsun_enftemp      
   type(block_data_real8_2d) :: f_rstfacsun_enfboreal    
   type(block_data_real8_2d) :: f_rstfacsun_dnfboreal    
   type(block_data_real8_2d) :: f_rstfacsun_ebftrop      
   type(block_data_real8_2d) :: f_rstfacsun_ebftemp      
   type(block_data_real8_2d) :: f_rstfacsun_dbftrop      
   type(block_data_real8_2d) :: f_rstfacsun_dbftemp      
   type(block_data_real8_2d) :: f_rstfacsun_dbfboreal    
   type(block_data_real8_2d) :: f_rstfacsun_ebstemp      
   type(block_data_real8_2d) :: f_rstfacsun_dbstemp      
   type(block_data_real8_2d) :: f_rstfacsun_dbsboreal    
   type(block_data_real8_2d) :: f_rstfacsun_c3arcgrass   
   type(block_data_real8_2d) :: f_rstfacsun_c3grass      
   type(block_data_real8_2d) :: f_rstfacsun_c4grass      
   type(block_data_real8_2d) :: f_rstfacsha_enftemp      
   type(block_data_real8_2d) :: f_rstfacsha_enfboreal    
   type(block_data_real8_2d) :: f_rstfacsha_dnfboreal    
   type(block_data_real8_2d) :: f_rstfacsha_ebftrop      
   type(block_data_real8_2d) :: f_rstfacsha_ebftemp      
   type(block_data_real8_2d) :: f_rstfacsha_dbftrop      
   type(block_data_real8_2d) :: f_rstfacsha_dbftemp      
   type(block_data_real8_2d) :: f_rstfacsha_dbfboreal    
   type(block_data_real8_2d) :: f_rstfacsha_ebstemp      
   type(block_data_real8_2d) :: f_rstfacsha_dbstemp      
   type(block_data_real8_2d) :: f_rstfacsha_dbsboreal    
   type(block_data_real8_2d) :: f_rstfacsha_c3arcgrass   
   type(block_data_real8_2d) :: f_rstfacsha_c3grass      
   type(block_data_real8_2d) :: f_rstfacsha_c4grass      
   type(block_data_real8_2d) :: f_lambdasun_enftemp      
   type(block_data_real8_2d) :: f_lambdasun_enfboreal    
   type(block_data_real8_2d) :: f_lambdasun_dnfboreal    
   type(block_data_real8_2d) :: f_lambdasun_ebftrop      
   type(block_data_real8_2d) :: f_lambdasun_ebftemp      
   type(block_data_real8_2d) :: f_lambdasun_dbftrop      
   type(block_data_real8_2d) :: f_lambdasun_dbftemp      
   type(block_data_real8_2d) :: f_lambdasun_dbfboreal    
   type(block_data_real8_2d) :: f_lambdasun_ebstemp      
   type(block_data_real8_2d) :: f_lambdasun_dbstemp      
   type(block_data_real8_2d) :: f_lambdasun_dbsboreal    
   type(block_data_real8_2d) :: f_lambdasun_c3arcgrass   
   type(block_data_real8_2d) :: f_lambdasun_c3grass      
   type(block_data_real8_2d) :: f_lambdasun_c4grass      
   type(block_data_real8_2d) :: f_lambdasha_enftemp      
   type(block_data_real8_2d) :: f_lambdasha_enfboreal    
   type(block_data_real8_2d) :: f_lambdasha_dnfboreal    
   type(block_data_real8_2d) :: f_lambdasha_ebftrop      
   type(block_data_real8_2d) :: f_lambdasha_ebftemp      
   type(block_data_real8_2d) :: f_lambdasha_dbftrop      
   type(block_data_real8_2d) :: f_lambdasha_dbftemp      
   type(block_data_real8_2d) :: f_lambdasha_dbfboreal    
   type(block_data_real8_2d) :: f_lambdasha_ebstemp      
   type(block_data_real8_2d) :: f_lambdasha_dbstemp      
   type(block_data_real8_2d) :: f_lambdasha_dbsboreal    
   type(block_data_real8_2d) :: f_lambdasha_c3arcgrass   
   type(block_data_real8_2d) :: f_lambdasha_c3grass      
   type(block_data_real8_2d) :: f_lambdasha_c4grass      
   type(block_data_real8_2d) :: f_lambda                 
#endif
#endif
   type(block_data_real8_2d) :: f_qcharge ! groundwater recharge rate [mm/s] 

   !---------------------------------------------------------------------
   type(block_data_real8_2d) :: f_t_grnd   ! ground surface temperature [K]
   type(block_data_real8_2d) :: f_tleaf    ! sunlit leaf temperature [K]
   type(block_data_real8_2d) :: f_ldew     ! depth of water on foliage [mm]
!#ifdef CLM5_INTERCEPTION
   type(block_data_real8_2d) :: f_ldew_rain     ! depth of rain on foliage [mm]
   type(block_data_real8_2d) :: f_ldew_snow     ! depth of snow on foliage [mm]
!#endif
   type(block_data_real8_2d) :: f_scv      ! snow cover, water equivalent [mm]
   type(block_data_real8_2d) :: f_snowdp   ! snow depth [meter]
   type(block_data_real8_2d) :: f_fsno     ! fraction of snow cover on ground
   type(block_data_real8_2d) :: f_sigf     ! fraction of veg cover, excluding snow-covered veg [-]
   type(block_data_real8_2d) :: f_green    ! leaf greenness
   type(block_data_real8_2d) :: f_lai      ! leaf area index
   type(block_data_real8_2d) :: f_laisun   ! sunlit leaf area index
   type(block_data_real8_2d) :: f_laisha   ! shaded leaf area index
   type(block_data_real8_2d) :: f_sai      ! stem area index
   type(block_data_real8_4d) :: f_alb      ! averaged albedo [visible, direct; direct, diffuse]
   type(block_data_real8_2d) :: f_emis     ! averaged bulk surface emissivity
   type(block_data_real8_2d) :: f_z0m      ! effective roughness [m]
   type(block_data_real8_2d) :: f_trad     ! radiative temperature of surface [K]
   type(block_data_real8_2d) :: f_tref     ! 2 m height air temperature [kelvin]
   type(block_data_real8_2d) :: f_qref     ! 2 m height air specific humidity [kg/kg]

   !---------------------------------------------------------------------
   type(block_data_real8_3d) :: f_t_soisno     ! soil temperature [K]
   type(block_data_real8_3d) :: f_wliq_soisno  ! liquid water in soil layers [kg/m2]
   type(block_data_real8_3d) :: f_wice_soisno  ! ice lens in soil layers [kg/m2]
   type(block_data_real8_3d) :: f_h2osoi       ! volumetric soil water in layers [m3/m3]
   type(block_data_real8_3d) :: f_rootr        ! water exchange between soil layers and root
   type(block_data_real8_3d) :: f_BD_all       ! bulk density in soil layers [kg/m3]
   type(block_data_real8_3d) :: f_wfc          ! water field capacity [m3/m3]
   type(block_data_real8_3d) :: f_OM_density   ! soil organic matter density [kg/m3]
#ifdef PLANT_HYDRAULIC_STRESS
   type(block_data_real8_3d) :: f_vegwp        ! vegetation water potential [mm]
#endif
   type(block_data_real8_2d) :: f_rstfacsun    ! factor of soil water stress
   type(block_data_real8_2d) :: f_rstfacsha    ! factor of soil water stress
   type(block_data_real8_2d) :: f_gs_sun    ! factor of soil water stress
   type(block_data_real8_2d) :: f_gs_sha    ! factor of soil water stress
   type(block_data_real8_2d) :: f_dpond        ! depth of ponding water [mm]
   type(block_data_real8_2d) :: f_zwt          ! the depth to water table [m]
   type(block_data_real8_2d) :: f_wa           ! water storage in aquifer [mm]
   type(block_data_real8_2d) :: f_wat          ! total water storage [mm]

   type(block_data_real8_3d) :: f_t_lake       ! lake temperature [K]
   type(block_data_real8_3d) :: f_lake_icefrac ! lake ice fraction cover [0-1]

   !---------------------------------------------------------------------
   type(block_data_real8_2d) :: f_ustar   ! u* in similarity theory [m/s]
   type(block_data_real8_2d) :: f_tstar   ! t* in similarity theory [kg/kg]
   type(block_data_real8_2d) :: f_qstar   ! q* in similarity theory [kg/kg]
   type(block_data_real8_2d) :: f_zol     ! dimensionless height (z/L) used in Monin-Obukhov theory
   type(block_data_real8_2d) :: f_rib     ! bulk Richardson number in surface layer
   type(block_data_real8_2d) :: f_fm      ! integral of profile function for momentum
   type(block_data_real8_2d) :: f_fh      ! integral of profile function for heat
   type(block_data_real8_2d) :: f_fq      ! integral of profile function for moisture
   type(block_data_real8_2d) :: f_us10m   ! 10m u-velocity [m/s]
   type(block_data_real8_2d) :: f_vs10m   ! 10m v-velocity [m/s]
   type(block_data_real8_2d) :: f_fm10m   ! integral of profile function for momentum at 10m [-]

   !---------------------------------------------------------------------
   type(block_data_real8_2d) :: f_xy_us      ! wind in eastward direction [m/s]
   type(block_data_real8_2d) :: f_xy_vs      ! wind in northward direction [m/s]
   type(block_data_real8_2d) :: f_xy_t       ! temperature at reference height [kelvin]
   type(block_data_real8_2d) :: f_xy_q       ! specific humidity at reference height [kg/kg]
   type(block_data_real8_2d) :: f_xy_prc     ! convective precipitation [mm/s]
   type(block_data_real8_2d) :: f_xy_prl     ! large scale precipitation [mm/s]
   type(block_data_real8_2d) :: f_xy_pbot    ! atmospheric pressure at the surface [pa]
   type(block_data_real8_2d) :: f_xy_frl     ! atmospheric infrared (longwave) radiation [W/m2]
   type(block_data_real8_2d) :: f_xy_solarin ! downward solar radiation at surface [W/m2]
   type(block_data_real8_2d) :: f_xy_rain    ! rain [mm/s]
   type(block_data_real8_2d) :: f_xy_snow    ! snow [mm/s]
   type(block_data_real8_2d) :: f_xy_ozone   ! ozone concentration [mol/mol]
   
   ! PUBLIC MEMBER FUNCTIONS:
   public :: allocate_2D_Fluxes

CONTAINS

   SUBROUTINE allocate_2D_Fluxes (grid)
      ! --------------------------------------------------------------------
      ! Allocates memory for CLM 2d [lon_points,lat_points] variables
      ! --------------------------------------------------------------------

      use spmd_task
      use mod_grid
      use mod_data_type
      USE mod_namelist
      implicit none

      type(grid_type), intent(in) :: grid

#if (defined UNSTRUCTURED || defined CATCHMENT) 
      IF (DEF_HISTORY_IN_VECTOR) THEN
         RETURN
      ENDIF 
#endif

      if (p_is_io) then
         
         call allocate_block_data (grid, f_taux   )  ! wind stress: E-W [kg/m/s2]
         call allocate_block_data (grid, f_tauy   )  ! wind stress: N-S [kg/m/s2]
         call allocate_block_data (grid, f_fsena  )  ! sensible heat from canopy height to atmosphere [W/m2]
         call allocate_block_data (grid, f_lfevpa )  ! latent heat flux from canopy height to atmosphere [W/m2]
         call allocate_block_data (grid, f_fevpa  )  ! evapotranspiration from canopy to atmosphere [mm/s]
         call allocate_block_data (grid, f_fsenl  )  ! sensible heat from leaves [W/m2]
         call allocate_block_data (grid, f_fevpl  )  ! evaporation+transpiration from leaves [mm/s]
         call allocate_block_data (grid, f_etr    )  ! transpiration rate [mm/s]
         call allocate_block_data (grid, f_fseng  )  ! sensible heat flux from ground [W/m2]
         call allocate_block_data (grid, f_fevpg  )  ! evaporation heat flux from ground [mm/s]
         call allocate_block_data (grid, f_fgrnd  )  ! ground heat flux [W/m2]
         call allocate_block_data (grid, f_sabvsun)  ! solar absorbed by sunlit canopy [W/m2]
         call allocate_block_data (grid, f_sabvsha)  ! solar absorbed by shaded [W/m2]
         call allocate_block_data (grid, f_sabg   )  ! solar absorbed by ground  [W/m2]
         call allocate_block_data (grid, f_sr     )  ! total reflected solar radiation (W/m2)
         call allocate_block_data (grid, f_solvd  )  ! incident direct beam vis solar radiation (W/m2)
         call allocate_block_data (grid, f_solvi  )  ! incident diffuse beam vis solar radiation (W/m2)
         call allocate_block_data (grid, f_solnd  )  ! incident direct beam nir solar radiation (W/m2)
         call allocate_block_data (grid, f_solni  )  ! incident diffuse beam nir solar radiation (W/m2)
         call allocate_block_data (grid, f_srvd   )  ! reflected direct beam vis solar radiation (W/m2)
         call allocate_block_data (grid, f_srvi   )  ! reflected diffuse beam vis solar radiation (W/m2)
         call allocate_block_data (grid, f_srnd   )  ! reflected direct beam nir solar radiation (W/m2)
         call allocate_block_data (grid, f_srni   )  ! reflected diffuse beam nir solar radiation (W/m2)
         call allocate_block_data (grid, f_solvdln)  ! incident direct beam vis solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_solviln)  ! incident diffuse beam vis solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_solndln)  ! incident direct beam nir solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_solniln)  ! incident diffuse beam nir solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_srvdln )  ! reflected direct beam vis solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_srviln )  ! reflected diffuse beam vis solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_srndln )  ! reflected direct beam nir solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_srniln )  ! reflected diffuse beam nir solar radiation at local noon(W/m2)
         call allocate_block_data (grid, f_olrg   )  ! outgoing long-wave radiation from ground+canopy [W/m2]
         call allocate_block_data (grid, f_rnet   )  ! net radiation [W/m2]
         call allocate_block_data (grid, f_xerr   )  ! the error of water banace [mm/s]
         call allocate_block_data (grid, f_zerr   )  ! the error of energy balance [W/m2]
         call allocate_block_data (grid, f_rsur   )  ! surface runoff [mm/s]
         call allocate_block_data (grid, f_rnof   )  ! total runoff [mm/s]
         call allocate_block_data (grid, f_qintr  )  ! interception [mm/s]
         call allocate_block_data (grid, f_qinfl  )  ! inflitration [mm/s]
         call allocate_block_data (grid, f_qdrip  )  ! throughfall [mm/s]
         call allocate_block_data (grid, f_assim  )  ! canopy assimilation rate [mol m-2 s-1]
         call allocate_block_data (grid, f_respc  )  ! respiration (plant+soil) [mol m-2 s-1]
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
         call allocate_block_data (grid, f_assim_RuBP_sun_enftemp        )
         call allocate_block_data (grid, f_assim_RuBP_sun_enfboreal      )
         call allocate_block_data (grid, f_assim_RuBP_sun_dnfboreal      )
         call allocate_block_data (grid, f_assim_RuBP_sun_ebftrop        )
         call allocate_block_data (grid, f_assim_RuBP_sun_ebftemp        )
         call allocate_block_data (grid, f_assim_RuBP_sun_dbftrop        )
         call allocate_block_data (grid, f_assim_RuBP_sun_dbftemp        )
         call allocate_block_data (grid, f_assim_RuBP_sun_dbfboreal      )
         call allocate_block_data (grid, f_assim_RuBP_sun_ebstemp        )
         call allocate_block_data (grid, f_assim_RuBP_sun_dbstemp        )
         call allocate_block_data (grid, f_assim_RuBP_sun_dbsboreal      )
         call allocate_block_data (grid, f_assim_RuBP_sun_c3arcgrass     )
         call allocate_block_data (grid, f_assim_RuBP_sun_c3grass        )
         call allocate_block_data (grid, f_assim_RuBP_sun_c4grass        )
         call allocate_block_data (grid, f_assim_RuBP_sha_enftemp        )
         call allocate_block_data (grid, f_assim_RuBP_sha_enfboreal      )
         call allocate_block_data (grid, f_assim_RuBP_sha_dnfboreal      )
         call allocate_block_data (grid, f_assim_RuBP_sha_ebftrop        )
         call allocate_block_data (grid, f_assim_RuBP_sha_ebftemp        )
         call allocate_block_data (grid, f_assim_RuBP_sha_dbftrop        )
         call allocate_block_data (grid, f_assim_RuBP_sha_dbftemp        )
         call allocate_block_data (grid, f_assim_RuBP_sha_dbfboreal      )
         call allocate_block_data (grid, f_assim_RuBP_sha_ebstemp        )
         call allocate_block_data (grid, f_assim_RuBP_sha_dbstemp        )
         call allocate_block_data (grid, f_assim_RuBP_sha_dbsboreal      )
         call allocate_block_data (grid, f_assim_RuBP_sha_c3arcgrass     )
         call allocate_block_data (grid, f_assim_RuBP_sha_c3grass        )
         call allocate_block_data (grid, f_assim_RuBP_sha_c4grass        )
         call allocate_block_data (grid, f_assim_Rubisco_sun_enftemp        )
         call allocate_block_data (grid, f_assim_Rubisco_sun_enfboreal      )
         call allocate_block_data (grid, f_assim_Rubisco_sun_dnfboreal      )
         call allocate_block_data (grid, f_assim_Rubisco_sun_ebftrop        )
         call allocate_block_data (grid, f_assim_Rubisco_sun_ebftemp        )
         call allocate_block_data (grid, f_assim_Rubisco_sun_dbftrop        )
         call allocate_block_data (grid, f_assim_Rubisco_sun_dbftemp        )
         call allocate_block_data (grid, f_assim_Rubisco_sun_dbfboreal      )
         call allocate_block_data (grid, f_assim_Rubisco_sun_ebstemp        )
         call allocate_block_data (grid, f_assim_Rubisco_sun_dbstemp        )
         call allocate_block_data (grid, f_assim_Rubisco_sun_dbsboreal      )
         call allocate_block_data (grid, f_assim_Rubisco_sun_c3arcgrass     )
         call allocate_block_data (grid, f_assim_Rubisco_sun_c3grass        )
         call allocate_block_data (grid, f_assim_Rubisco_sun_c4grass        )
         call allocate_block_data (grid, f_assim_Rubisco_sha_enftemp        )
         call allocate_block_data (grid, f_assim_Rubisco_sha_enfboreal      )
         call allocate_block_data (grid, f_assim_Rubisco_sha_dnfboreal      )
         call allocate_block_data (grid, f_assim_Rubisco_sha_ebftrop        )
         call allocate_block_data (grid, f_assim_Rubisco_sha_ebftemp        )
         call allocate_block_data (grid, f_assim_Rubisco_sha_dbftrop        )
         call allocate_block_data (grid, f_assim_Rubisco_sha_dbftemp        )
         call allocate_block_data (grid, f_assim_Rubisco_sha_dbfboreal      )
         call allocate_block_data (grid, f_assim_Rubisco_sha_ebstemp        )
         call allocate_block_data (grid, f_assim_Rubisco_sha_dbstemp        )
         call allocate_block_data (grid, f_assim_Rubisco_sha_dbsboreal      )
         call allocate_block_data (grid, f_assim_Rubisco_sha_c3arcgrass     )
         call allocate_block_data (grid, f_assim_Rubisco_sha_c3grass        )
         call allocate_block_data (grid, f_assim_Rubisco_sha_c4grass        )
         call allocate_block_data (grid, f_assimsun_enftemp        )
         call allocate_block_data (grid, f_assimsun_enfboreal      )
         call allocate_block_data (grid, f_assimsun_dnfboreal      )
         call allocate_block_data (grid, f_assimsun_ebftrop        )
         call allocate_block_data (grid, f_assimsun_ebftemp        )
         call allocate_block_data (grid, f_assimsun_dbftrop        )
         call allocate_block_data (grid, f_assimsun_dbftemp        )
         call allocate_block_data (grid, f_assimsun_dbfboreal      )
         call allocate_block_data (grid, f_assimsun_ebstemp        )
         call allocate_block_data (grid, f_assimsun_dbstemp        )
         call allocate_block_data (grid, f_assimsun_dbsboreal      )
         call allocate_block_data (grid, f_assimsun_c3arcgrass     )
         call allocate_block_data (grid, f_assimsun_c3grass        )
         call allocate_block_data (grid, f_assimsun_c4grass        )
         call allocate_block_data (grid, f_assimsha_enftemp        )
         call allocate_block_data (grid, f_assimsha_enfboreal      )
         call allocate_block_data (grid, f_assimsha_dnfboreal      )
         call allocate_block_data (grid, f_assimsha_ebftrop        )
         call allocate_block_data (grid, f_assimsha_ebftemp        )
         call allocate_block_data (grid, f_assimsha_dbftrop        )
         call allocate_block_data (grid, f_assimsha_dbftemp        )
         call allocate_block_data (grid, f_assimsha_dbfboreal      )
         call allocate_block_data (grid, f_assimsha_ebstemp        )
         call allocate_block_data (grid, f_assimsha_dbstemp        )
         call allocate_block_data (grid, f_assimsha_dbsboreal      )
         call allocate_block_data (grid, f_assimsha_c3arcgrass     )
         call allocate_block_data (grid, f_assimsha_c3grass        )
         call allocate_block_data (grid, f_assimsha_c4grass        )
         call allocate_block_data (grid, f_etrsun_enftemp        )
         call allocate_block_data (grid, f_etrsun_enfboreal      )
         call allocate_block_data (grid, f_etrsun_dnfboreal      )
         call allocate_block_data (grid, f_etrsun_ebftrop        )
         call allocate_block_data (grid, f_etrsun_ebftemp        )
         call allocate_block_data (grid, f_etrsun_dbftrop        )
         call allocate_block_data (grid, f_etrsun_dbftemp        )
         call allocate_block_data (grid, f_etrsun_dbfboreal      )
         call allocate_block_data (grid, f_etrsun_ebstemp        )
         call allocate_block_data (grid, f_etrsun_dbstemp        )
         call allocate_block_data (grid, f_etrsun_dbsboreal      )
         call allocate_block_data (grid, f_etrsun_c3arcgrass     )
         call allocate_block_data (grid, f_etrsun_c3grass        )
         call allocate_block_data (grid, f_etrsun_c4grass        )
         call allocate_block_data (grid, f_etrsha_enftemp        )
         call allocate_block_data (grid, f_etrsha_enfboreal      )
         call allocate_block_data (grid, f_etrsha_dnfboreal      )
         call allocate_block_data (grid, f_etrsha_ebftrop        )
         call allocate_block_data (grid, f_etrsha_ebftemp        )
         call allocate_block_data (grid, f_etrsha_dbftrop        )
         call allocate_block_data (grid, f_etrsha_dbftemp        )
         call allocate_block_data (grid, f_etrsha_dbfboreal      )
         call allocate_block_data (grid, f_etrsha_ebstemp        )
         call allocate_block_data (grid, f_etrsha_dbstemp        )
         call allocate_block_data (grid, f_etrsha_dbsboreal      )
         call allocate_block_data (grid, f_etrsha_c3arcgrass     )
         call allocate_block_data (grid, f_etrsha_c3grass        )
         call allocate_block_data (grid, f_etrsha_c4grass        )
         call allocate_block_data (grid, f_cisun_enftemp        )
         call allocate_block_data (grid, f_cisun_enfboreal      )
         call allocate_block_data (grid, f_cisun_dnfboreal      )
         call allocate_block_data (grid, f_cisun_ebftrop        )
         call allocate_block_data (grid, f_cisun_ebftemp        )
         call allocate_block_data (grid, f_cisun_dbftrop        )
         call allocate_block_data (grid, f_cisun_dbftemp        )
         call allocate_block_data (grid, f_cisun_dbfboreal      )
         call allocate_block_data (grid, f_cisun_ebstemp        )
         call allocate_block_data (grid, f_cisun_dbstemp        )
         call allocate_block_data (grid, f_cisun_dbsboreal      )
         call allocate_block_data (grid, f_cisun_c3arcgrass     )
         call allocate_block_data (grid, f_cisun_c3grass        )
         call allocate_block_data (grid, f_cisun_c4grass        )
         call allocate_block_data (grid, f_cisha_enftemp        )
         call allocate_block_data (grid, f_cisha_enfboreal      )
         call allocate_block_data (grid, f_cisha_dnfboreal      )
         call allocate_block_data (grid, f_cisha_ebftrop        )
         call allocate_block_data (grid, f_cisha_ebftemp        )
         call allocate_block_data (grid, f_cisha_dbftrop        )
         call allocate_block_data (grid, f_cisha_dbftemp        )
         call allocate_block_data (grid, f_cisha_dbfboreal      )
         call allocate_block_data (grid, f_cisha_ebstemp        )
         call allocate_block_data (grid, f_cisha_dbstemp        )
         call allocate_block_data (grid, f_cisha_dbsboreal      )
         call allocate_block_data (grid, f_cisha_c3arcgrass     )
         call allocate_block_data (grid, f_cisha_c3grass        )
         call allocate_block_data (grid, f_cisha_c4grass        )
         call allocate_block_data (grid, f_essun_enftemp        )
         call allocate_block_data (grid, f_essun_enfboreal      )
         call allocate_block_data (grid, f_essun_dnfboreal      )
         call allocate_block_data (grid, f_essun_ebftrop        )
         call allocate_block_data (grid, f_essun_ebftemp        )
         call allocate_block_data (grid, f_essun_dbftrop        )
         call allocate_block_data (grid, f_essun_dbftemp        )
         call allocate_block_data (grid, f_essun_dbfboreal      )
         call allocate_block_data (grid, f_essun_ebstemp        )
         call allocate_block_data (grid, f_essun_dbstemp        )
         call allocate_block_data (grid, f_essun_dbsboreal      )
         call allocate_block_data (grid, f_essun_c3arcgrass     )
         call allocate_block_data (grid, f_essun_c3grass        )
         call allocate_block_data (grid, f_essun_c4grass        )
         call allocate_block_data (grid, f_essha_enftemp        )
         call allocate_block_data (grid, f_essha_enfboreal      )
         call allocate_block_data (grid, f_essha_dnfboreal      )
         call allocate_block_data (grid, f_essha_ebftrop        )
         call allocate_block_data (grid, f_essha_ebftemp        )
         call allocate_block_data (grid, f_essha_dbftrop        )
         call allocate_block_data (grid, f_essha_dbftemp        )
         call allocate_block_data (grid, f_essha_dbfboreal      )
         call allocate_block_data (grid, f_essha_ebstemp        )
         call allocate_block_data (grid, f_essha_dbstemp        )
         call allocate_block_data (grid, f_essha_dbsboreal      )
         call allocate_block_data (grid, f_essha_c3arcgrass     )
         call allocate_block_data (grid, f_essha_c3grass        )
         call allocate_block_data (grid, f_essha_c4grass        )
         call allocate_block_data (grid, f_gssun_enftemp        )
         call allocate_block_data (grid, f_gssun_enfboreal      )
         call allocate_block_data (grid, f_gssun_dnfboreal      )
         call allocate_block_data (grid, f_gssun_ebftrop        )
         call allocate_block_data (grid, f_gssun_ebftemp        )
         call allocate_block_data (grid, f_gssun_dbftrop        )
         call allocate_block_data (grid, f_gssun_dbftemp        )
         call allocate_block_data (grid, f_gssun_dbfboreal      )
         call allocate_block_data (grid, f_gssun_ebstemp        )
         call allocate_block_data (grid, f_gssun_dbstemp        )
         call allocate_block_data (grid, f_gssun_dbsboreal      )
         call allocate_block_data (grid, f_gssun_c3arcgrass     )
         call allocate_block_data (grid, f_gssun_c3grass        )
         call allocate_block_data (grid, f_gssun_c4grass        )
         call allocate_block_data (grid, f_gssha_enftemp        )
         call allocate_block_data (grid, f_gssha_enfboreal      )
         call allocate_block_data (grid, f_gssha_dnfboreal      )
         call allocate_block_data (grid, f_gssha_ebftrop        )
         call allocate_block_data (grid, f_gssha_ebftemp        )
         call allocate_block_data (grid, f_gssha_dbftrop        )
         call allocate_block_data (grid, f_gssha_dbftemp        )
         call allocate_block_data (grid, f_gssha_dbfboreal      )
         call allocate_block_data (grid, f_gssha_ebstemp        )
         call allocate_block_data (grid, f_gssha_dbstemp        )
         call allocate_block_data (grid, f_gssha_dbsboreal      )
         call allocate_block_data (grid, f_gssha_c3arcgrass     )
         call allocate_block_data (grid, f_gssha_c3grass        )
         call allocate_block_data (grid, f_gssha_c4grass        )
         call allocate_block_data (grid, f_gammasun_enftemp        )
         call allocate_block_data (grid, f_gammasun_enfboreal      )
         call allocate_block_data (grid, f_gammasun_dnfboreal      )
         call allocate_block_data (grid, f_gammasun_ebftrop        )
         call allocate_block_data (grid, f_gammasun_ebftemp        )
         call allocate_block_data (grid, f_gammasun_dbftrop        )
         call allocate_block_data (grid, f_gammasun_dbftemp        )
         call allocate_block_data (grid, f_gammasun_dbfboreal      )
         call allocate_block_data (grid, f_gammasun_ebstemp        )
         call allocate_block_data (grid, f_gammasun_dbstemp        )
         call allocate_block_data (grid, f_gammasun_dbsboreal      )
         call allocate_block_data (grid, f_gammasun_c3arcgrass     )
         call allocate_block_data (grid, f_gammasun_c3grass        )
         call allocate_block_data (grid, f_gammasun_c4grass        )
         call allocate_block_data (grid, f_gammasha_enftemp        )
         call allocate_block_data (grid, f_gammasha_enfboreal      )
         call allocate_block_data (grid, f_gammasha_dnfboreal      )
         call allocate_block_data (grid, f_gammasha_ebftrop        )
         call allocate_block_data (grid, f_gammasha_ebftemp        )
         call allocate_block_data (grid, f_gammasha_dbftrop        )
         call allocate_block_data (grid, f_gammasha_dbftemp        )
         call allocate_block_data (grid, f_gammasha_dbfboreal      )
         call allocate_block_data (grid, f_gammasha_ebstemp        )
         call allocate_block_data (grid, f_gammasha_dbstemp        )
         call allocate_block_data (grid, f_gammasha_dbsboreal      )
         call allocate_block_data (grid, f_gammasha_c3arcgrass     )
         call allocate_block_data (grid, f_gammasha_c3grass        )
         call allocate_block_data (grid, f_gammasha_c4grass        )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_enftemp          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_enfboreal        )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_dnfboreal        )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_ebftrop          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_ebftemp          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_dbftrop          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_dbftemp          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_dbfboreal        )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_ebstemp          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_dbstemp          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_dbsboreal        )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_c3arcgrass       )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_c3grass          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sun_c4grass          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_enftemp          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_enfboreal        )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_dnfboreal        )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_ebftrop          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_ebftemp          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_dbftrop          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_dbftemp          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_dbfboreal        )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_ebstemp          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_dbstemp          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_dbsboreal        )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_c3arcgrass       )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_c3grass          )
         call allocate_block_data (grid, f_RuBPlimitfrac_sha_c4grass          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_enftemp          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_enfboreal        )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_dnfboreal        )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_ebftrop          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_ebftemp          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_dbftrop          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_dbftemp          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_dbfboreal        )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_ebstemp          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_dbstemp          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_dbsboreal        )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_c3arcgrass       )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_c3grass          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sun_c4grass          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_enftemp          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_enfboreal        )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_dnfboreal        )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_ebftrop          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_ebftemp          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_dbftrop          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_dbftemp          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_dbfboreal        )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_ebstemp          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_dbstemp          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_dbsboreal        )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_c3arcgrass       )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_c3grass          )
         call allocate_block_data (grid, f_Rubiscolimitfrac_sha_c4grass          )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_enftemp          )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_enfboreal        )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_dnfboreal        )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_ebftrop          )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_ebftemp          )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_dbftrop          )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_dbftemp          )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_dbfboreal        )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_ebstemp          )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_dbstemp          )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_dbsboreal        )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_c3arcgrass       )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_c3grass          )
         call allocate_block_data (grid, f_Sinklimitfrac_sun_c4grass          )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_enftemp          )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_enfboreal        )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_dnfboreal        )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_ebftrop          )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_ebftemp          )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_dbftrop          )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_dbftemp          )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_dbfboreal        )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_ebstemp          )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_dbstemp          )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_dbsboreal        )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_c3arcgrass       )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_c3grass          )
         call allocate_block_data (grid, f_Sinklimitfrac_sha_c4grass          )
         call allocate_block_data (grid, f_rstfacsun_enftemp        )
         call allocate_block_data (grid, f_rstfacsun_enfboreal      )
         call allocate_block_data (grid, f_rstfacsun_dnfboreal      )
         call allocate_block_data (grid, f_rstfacsun_ebftrop        )
         call allocate_block_data (grid, f_rstfacsun_ebftemp        )
         call allocate_block_data (grid, f_rstfacsun_dbftrop        )
         call allocate_block_data (grid, f_rstfacsun_dbftemp        )
         call allocate_block_data (grid, f_rstfacsun_dbfboreal      )
         call allocate_block_data (grid, f_rstfacsun_ebstemp        )
         call allocate_block_data (grid, f_rstfacsun_dbstemp        )
         call allocate_block_data (grid, f_rstfacsun_dbsboreal      )
         call allocate_block_data (grid, f_rstfacsun_c3arcgrass     )
         call allocate_block_data (grid, f_rstfacsun_c3grass        )
         call allocate_block_data (grid, f_rstfacsun_c4grass        )
         call allocate_block_data (grid, f_rstfacsha_enftemp        )
         call allocate_block_data (grid, f_rstfacsha_enfboreal      )
         call allocate_block_data (grid, f_rstfacsha_dnfboreal      )
         call allocate_block_data (grid, f_rstfacsha_ebftrop        )
         call allocate_block_data (grid, f_rstfacsha_ebftemp        )
         call allocate_block_data (grid, f_rstfacsha_dbftrop        )
         call allocate_block_data (grid, f_rstfacsha_dbftemp        )
         call allocate_block_data (grid, f_rstfacsha_dbfboreal      )
         call allocate_block_data (grid, f_rstfacsha_ebstemp        )
         call allocate_block_data (grid, f_rstfacsha_dbstemp        )
         call allocate_block_data (grid, f_rstfacsha_dbsboreal      )
         call allocate_block_data (grid, f_rstfacsha_c3arcgrass     )
         call allocate_block_data (grid, f_rstfacsha_c3grass        )
         call allocate_block_data (grid, f_rstfacsha_c4grass        )
         call allocate_block_data (grid, f_lambdasun_enftemp        )
         call allocate_block_data (grid, f_lambdasun_enfboreal      )
         call allocate_block_data (grid, f_lambdasun_dnfboreal      )
         call allocate_block_data (grid, f_lambdasun_ebftrop        )
         call allocate_block_data (grid, f_lambdasun_ebftemp        )
         call allocate_block_data (grid, f_lambdasun_dbftrop        )
         call allocate_block_data (grid, f_lambdasun_dbftemp        )
         call allocate_block_data (grid, f_lambdasun_dbfboreal      )
         call allocate_block_data (grid, f_lambdasun_ebstemp        )
         call allocate_block_data (grid, f_lambdasun_dbstemp        )
         call allocate_block_data (grid, f_lambdasun_dbsboreal      )
         call allocate_block_data (grid, f_lambdasun_c3arcgrass     )
         call allocate_block_data (grid, f_lambdasun_c3grass        )
         call allocate_block_data (grid, f_lambdasun_c4grass        )
         call allocate_block_data (grid, f_lambdasha_enftemp        )
         call allocate_block_data (grid, f_lambdasha_enfboreal      )
         call allocate_block_data (grid, f_lambdasha_dnfboreal      )
         call allocate_block_data (grid, f_lambdasha_ebftrop        )
         call allocate_block_data (grid, f_lambdasha_ebftemp        )
         call allocate_block_data (grid, f_lambdasha_dbftrop        )
         call allocate_block_data (grid, f_lambdasha_dbftemp        )
         call allocate_block_data (grid, f_lambdasha_dbfboreal      )
         call allocate_block_data (grid, f_lambdasha_ebstemp        )
         call allocate_block_data (grid, f_lambdasha_dbstemp        )
         call allocate_block_data (grid, f_lambdasha_dbsboreal      )
         call allocate_block_data (grid, f_lambdasha_c3arcgrass     )
         call allocate_block_data (grid, f_lambdasha_c3grass        )
         call allocate_block_data (grid, f_lambdasha_c4grass        )
         call allocate_block_data (grid, f_lambda                   )
#endif
#endif
         call allocate_block_data (grid, f_qcharge)  ! groundwater recharge rate [mm/s] 

         !---------------------------------------------------------------------
         call allocate_block_data (grid, f_t_grnd  )  ! ground surface temperature [K]
         call allocate_block_data (grid, f_tleaf   )  ! sunlit leaf temperature [K]
         call allocate_block_data (grid, f_ldew    )  ! depth of water on foliage [mm]
         call allocate_block_data (grid, f_ldew_rain    )  ! depth of rain on foliage [mm]
         call allocate_block_data (grid, f_ldew_snow    )  ! depth of snow on foliage [mm]
         call allocate_block_data (grid, f_scv     )  ! snow cover, water equivalent [mm]
         call allocate_block_data (grid, f_snowdp  )  ! snow depth [meter]
         call allocate_block_data (grid, f_fsno    )  ! fraction of snow cover on ground
         call allocate_block_data (grid, f_sigf    )  ! fraction of veg cover, excluding snow-covered veg [-]
         call allocate_block_data (grid, f_green   )  ! leaf greenness
         call allocate_block_data (grid, f_lai     )  ! leaf area index
         call allocate_block_data (grid, f_laisun  )  ! sunlit leaf area index
         call allocate_block_data (grid, f_laisha  )  ! shaded leaf area index
         call allocate_block_data (grid, f_sai     )  ! stem area index
         call allocate_block_data (grid, f_alb, 2, 2) ! averaged albedo [visible, direct; direct, diffuse]
         call allocate_block_data (grid, f_emis    )  ! averaged bulk surface emissivity
         call allocate_block_data (grid, f_z0m     )  ! effective roughness [m]
         call allocate_block_data (grid, f_trad    )  ! radiative temperature of surface [K]
         call allocate_block_data (grid, f_tref    )  ! 2 m height air temperature [kelvin]
         call allocate_block_data (grid, f_qref    )  ! 2 m height air specific humidity [kg/kg]

         !---------------------------------------------------------------------
         call allocate_block_data (grid, f_t_soisno   , nl_soil-maxsnl, maxsnl+1)  ! soil temperature [K]
         call allocate_block_data (grid, f_wliq_soisno, nl_soil-maxsnl, maxsnl+1)  ! liquid water in soil layers [kg/m2]
         call allocate_block_data (grid, f_wice_soisno, nl_soil-maxsnl, maxsnl+1)  ! ice lens in soil layers [kg/m2]
         call allocate_block_data (grid, f_h2osoi     , nl_soil)  ! volumetric soil water in layers [m3/m3]
         call allocate_block_data (grid, f_rootr      , nl_soil)  ! water exchange between soil layers and root
         call allocate_block_data (grid, f_BD_all     , nl_soil)
         call allocate_block_data (grid, f_wfc        , nl_soil)
         call allocate_block_data (grid, f_OM_density , nl_soil)
#ifdef PLANT_HYDRAULIC_STRESS
         call allocate_block_data (grid, f_vegwp      , nvegwcs)  ! vegetation water potential [mm]
#endif
         call allocate_block_data (grid, f_rstfacsun)  ! factor of soil water stress
         call allocate_block_data (grid, f_rstfacsha)  ! factor of soil water stress
         call allocate_block_data (grid, f_gs_sun)  ! factor of soil water stress
         call allocate_block_data (grid, f_gs_sha)  ! factor of soil water stress

         call allocate_block_data (grid, f_dpond  )  ! depth of ponding water [m]

         call allocate_block_data (grid, f_zwt   )  ! the depth to water table [m]
         call allocate_block_data (grid, f_wa    )  ! water storage in aquifer [mm]
         call allocate_block_data (grid, f_wat   )  ! total water storage [mm]

         call allocate_block_data (grid, f_t_lake      , nl_lake)  ! lake temperature [K]
         call allocate_block_data (grid, f_lake_icefrac, nl_lake)  ! lake ice fraction cover [0-1]

         !---------------------------------------------------------------------
         call allocate_block_data (grid, f_ustar)  ! u* in similarity theory [m/s]
         call allocate_block_data (grid, f_tstar)  ! t* in similarity theory [kg/kg]
         call allocate_block_data (grid, f_qstar)  ! q* in similarity theory [kg/kg]
         call allocate_block_data (grid, f_zol  )  ! dimensionless height (z/L) used in Monin-Obukhov theory
         call allocate_block_data (grid, f_rib  )  ! bulk Richardson number in surface layer
         call allocate_block_data (grid, f_fm   )  ! integral of profile function for momentum
         call allocate_block_data (grid, f_fh   )  ! integral of profile function for heat
         call allocate_block_data (grid, f_fq   )  ! integral of profile function for moisture
         call allocate_block_data (grid, f_us10m)  ! 10m u-velocity [m/s]
         call allocate_block_data (grid, f_vs10m)  ! 10m v-velocity [m/s]
         call allocate_block_data (grid, f_fm10m)  ! integral of profile function for momentum at 10m [-]

         !---------------------------------------------------------------------
         call allocate_block_data (grid, f_xy_us     )  ! wind in eastward direction [m/s]
         call allocate_block_data (grid, f_xy_vs     )  ! wind in northward direction [m/s]
         call allocate_block_data (grid, f_xy_t      )  ! temperature at reference height [kelvin]
         call allocate_block_data (grid, f_xy_q      )  ! specific humidity at reference height [kg/kg]
         call allocate_block_data (grid, f_xy_prc    )  ! convective precipitation [mm/s]
         call allocate_block_data (grid, f_xy_prl    )  ! large scale precipitation [mm/s]
         call allocate_block_data (grid, f_xy_pbot   )  ! atmospheric pressure at the surface [pa]
         call allocate_block_data (grid, f_xy_frl    )  ! atmospheric infrared (longwave) radiation [W/m2]
         call allocate_block_data (grid, f_xy_solarin)  ! downward solar radiation at surface [W/m2]
         call allocate_block_data (grid, f_xy_rain   )  ! rain [mm/s]
         call allocate_block_data (grid, f_xy_snow   )  ! snow [mm/s]
         call allocate_block_data (grid, f_xy_ozone  )  ! ozone concentration [mol/mol]

      end if

#ifdef BGC
      CALL allocate_2D_BGCFluxes (grid)
#endif

   END SUBROUTINE allocate_2D_Fluxes

END MODULE MOD_2D_Fluxes
