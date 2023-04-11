#include <define.h>

module MOD_1D_Acc_Fluxes

   use precision

   real(r8) :: nac              ! number of accumulation
   real(r8), allocatable :: nac_ln   (:)

   real(r8), allocatable :: a_us     (:)  
   real(r8), allocatable :: a_vs     (:)
   real(r8), allocatable :: a_t      (:)
   real(r8), allocatable :: a_q      (:)
   real(r8), allocatable :: a_prc    (:)
   real(r8), allocatable :: a_prl    (:)
   real(r8), allocatable :: a_pbot   (:)
   real(r8), allocatable :: a_frl    (:)
   real(r8), allocatable :: a_solarin(:)

   real(r8), allocatable :: a_taux   (:)
   real(r8), allocatable :: a_tauy   (:)
   real(r8), allocatable :: a_fsena  (:)
   real(r8), allocatable :: a_lfevpa (:)
   real(r8), allocatable :: a_fevpa  (:)
   real(r8), allocatable :: a_fsenl  (:)
   real(r8), allocatable :: a_fevpl  (:)
   real(r8), allocatable :: a_etr    (:)
   real(r8), allocatable :: a_fseng  (:)
   real(r8), allocatable :: a_fevpg  (:)
   real(r8), allocatable :: a_fgrnd  (:)
   real(r8), allocatable :: a_sabvsun(:)  
   real(r8), allocatable :: a_sabvsha(:) 
   real(r8), allocatable :: a_sabg   (:)
   real(r8), allocatable :: a_olrg   (:)
   real(r8), allocatable :: a_rnet   (:)
   real(r8), allocatable :: a_xerr   (:)
   real(r8), allocatable :: a_zerr   (:)
   real(r8), allocatable :: a_rsur   (:)
   real(r8), allocatable :: a_rnof   (:)
   real(r8), allocatable :: a_qintr  (:)
   real(r8), allocatable :: a_qinfl  (:)
   real(r8), allocatable :: a_qdrip  (:)
   real(r8), allocatable :: a_rstfacsun (:)
   real(r8), allocatable :: a_rstfacsha (:)
   real(r8), allocatable :: a_gs_sun (:)
   real(r8), allocatable :: a_gs_sha (:)
   real(r8), allocatable :: a_dpond  (:)
   real(r8), allocatable :: a_zwt    (:)
   real(r8), allocatable :: a_wa     (:)
   real(r8), allocatable :: a_wat    (:)
   real(r8), allocatable :: a_assim  (:)
   real(r8), allocatable :: a_respc  (:)

   real(r8), allocatable :: a_qcharge(:)

   real(r8), allocatable :: a_t_grnd(:) 
   real(r8), allocatable :: a_tleaf (:) 
   real(r8), allocatable :: a_ldew  (:) 
   real(r8), allocatable :: a_ldew_rain  (:) 
   real(r8), allocatable :: a_ldew_snow  (:) 
   real(r8), allocatable :: a_scv   (:) 
   real(r8), allocatable :: a_snowdp(:) 
   real(r8), allocatable :: a_fsno  (:) 
   real(r8), allocatable :: a_sigf  (:) 
   real(r8), allocatable :: a_green (:) 
   real(r8), allocatable :: a_lai   (:) 
   real(r8), allocatable :: a_laisun(:) 
   real(r8), allocatable :: a_laisha(:) 
   real(r8), allocatable :: a_sai   (:) 

   real(r8), allocatable :: a_alb(:,:,:)    

   real(r8), allocatable :: a_emis (:)
   real(r8), allocatable :: a_z0m  (:)
   real(r8), allocatable :: a_trad (:)
   real(r8), allocatable :: a_tref (:)
   real(r8), allocatable :: a_qref (:)
   real(r8), allocatable :: a_rain (:)
   real(r8), allocatable :: a_snow (:)  

#ifdef BGC
   real(r8), allocatable :: a_leafc              (:)
   real(r8), allocatable :: a_leafc_storage      (:)
   real(r8), allocatable :: a_leafc_xfer         (:)
   real(r8), allocatable :: a_frootc             (:)
   real(r8), allocatable :: a_frootc_storage     (:)
   real(r8), allocatable :: a_frootc_xfer        (:)
   real(r8), allocatable :: a_livestemc          (:)
   real(r8), allocatable :: a_livestemc_storage  (:)
   real(r8), allocatable :: a_livestemc_xfer     (:)
   real(r8), allocatable :: a_deadstemc          (:)
   real(r8), allocatable :: a_deadstemc_storage  (:)
   real(r8), allocatable :: a_deadstemc_xfer     (:)
   real(r8), allocatable :: a_livecrootc         (:)
   real(r8), allocatable :: a_livecrootc_storage (:)
   real(r8), allocatable :: a_livecrootc_xfer    (:)
   real(r8), allocatable :: a_deadcrootc         (:)
   real(r8), allocatable :: a_deadcrootc_storage (:)
   real(r8), allocatable :: a_deadcrootc_xfer    (:)
   real(r8), allocatable :: a_grainc             (:)
   real(r8), allocatable :: a_grainc_storage     (:)
   real(r8), allocatable :: a_grainc_xfer        (:)
   real(r8), allocatable :: a_leafn              (:)
   real(r8), allocatable :: a_leafn_storage      (:)
   real(r8), allocatable :: a_leafn_xfer         (:)
   real(r8), allocatable :: a_frootn             (:)
   real(r8), allocatable :: a_frootn_storage     (:)
   real(r8), allocatable :: a_frootn_xfer        (:)
   real(r8), allocatable :: a_livestemn          (:)
   real(r8), allocatable :: a_livestemn_storage  (:)
   real(r8), allocatable :: a_livestemn_xfer     (:)
   real(r8), allocatable :: a_deadstemn          (:)
   real(r8), allocatable :: a_deadstemn_storage  (:)
   real(r8), allocatable :: a_deadstemn_xfer     (:)
   real(r8), allocatable :: a_livecrootn         (:)
   real(r8), allocatable :: a_livecrootn_storage (:)
   real(r8), allocatable :: a_livecrootn_xfer    (:)
   real(r8), allocatable :: a_deadcrootn         (:)
   real(r8), allocatable :: a_deadcrootn_storage (:)
   real(r8), allocatable :: a_deadcrootn_xfer    (:)
   real(r8), allocatable :: a_grainn             (:)
   real(r8), allocatable :: a_grainn_storage     (:)
   real(r8), allocatable :: a_grainn_xfer        (:)
   real(r8), allocatable :: a_retransn           (:)
   real(r8), allocatable :: a_gpp                (:)
   real(r8), allocatable :: a_downreg            (:)
   real(r8), allocatable :: a_ar                 (:)
   real(r8), allocatable :: a_cwdprod            (:)
   real(r8), allocatable :: a_cwddecomp          (:)
   real(r8), allocatable :: a_hr                 (:)
   real(r8), allocatable :: a_fpg                (:)
   real(r8), allocatable :: a_fpi                (:)
   real(r8), allocatable :: a_gpp_enftemp        (:) !1
   real(r8), allocatable :: a_gpp_enfboreal      (:) !2
   real(r8), allocatable :: a_gpp_dnfboreal      (:) !3
   real(r8), allocatable :: a_gpp_ebftrop        (:) !4
   real(r8), allocatable :: a_gpp_ebftemp        (:) !5
   real(r8), allocatable :: a_gpp_dbftrop        (:) !6
   real(r8), allocatable :: a_gpp_dbftemp        (:) !7
   real(r8), allocatable :: a_gpp_dbfboreal      (:) !8
   real(r8), allocatable :: a_gpp_ebstemp        (:) !9
   real(r8), allocatable :: a_gpp_dbstemp        (:) !10
   real(r8), allocatable :: a_gpp_dbsboreal      (:) !11
   real(r8), allocatable :: a_gpp_c3arcgrass     (:) !12
   real(r8), allocatable :: a_gpp_c3grass        (:) !13
   real(r8), allocatable :: a_gpp_c4grass        (:) !14
   real(r8), allocatable :: a_leafc_enftemp        (:) !1
   real(r8), allocatable :: a_leafc_enfboreal    (:) !2
   real(r8), allocatable :: a_leafc_dnfboreal    (:) !3
   real(r8), allocatable :: a_leafc_ebftrop      (:) !4
   real(r8), allocatable :: a_leafc_ebftemp      (:) !5
   real(r8), allocatable :: a_leafc_dbftrop      (:) !6
   real(r8), allocatable :: a_leafc_dbftemp      (:) !7
   real(r8), allocatable :: a_leafc_dbfboreal    (:) !8
   real(r8), allocatable :: a_leafc_ebstemp      (:) !9
   real(r8), allocatable :: a_leafc_dbstemp      (:) !10
   real(r8), allocatable :: a_leafc_dbsboreal    (:) !11
   real(r8), allocatable :: a_leafc_c3arcgrass   (:) !12
   real(r8), allocatable :: a_leafc_c3grass      (:) !13
   real(r8), allocatable :: a_leafc_c4grass      (:) !14
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
   real(r8), allocatable :: a_assim_RuBP_sun_enftemp        (:) !1
   real(r8), allocatable :: a_assim_RuBP_sun_enfboreal      (:) !2
   real(r8), allocatable :: a_assim_RuBP_sun_dnfboreal      (:) !3
   real(r8), allocatable :: a_assim_RuBP_sun_ebftrop        (:) !4
   real(r8), allocatable :: a_assim_RuBP_sun_ebftemp        (:) !5
   real(r8), allocatable :: a_assim_RuBP_sun_dbftrop        (:) !6
   real(r8), allocatable :: a_assim_RuBP_sun_dbftemp        (:) !7
   real(r8), allocatable :: a_assim_RuBP_sun_dbfboreal      (:) !8
   real(r8), allocatable :: a_assim_RuBP_sun_ebstemp        (:) !9
   real(r8), allocatable :: a_assim_RuBP_sun_dbstemp        (:) !10
   real(r8), allocatable :: a_assim_RuBP_sun_dbsboreal      (:) !11
   real(r8), allocatable :: a_assim_RuBP_sun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_assim_RuBP_sun_c3grass        (:) !13
   real(r8), allocatable :: a_assim_RuBP_sun_c4grass        (:) !14
   real(r8), allocatable :: a_assim_RuBP_sha_enftemp        (:) !1
   real(r8), allocatable :: a_assim_RuBP_sha_enfboreal      (:) !2
   real(r8), allocatable :: a_assim_RuBP_sha_dnfboreal      (:) !3
   real(r8), allocatable :: a_assim_RuBP_sha_ebftrop        (:) !4
   real(r8), allocatable :: a_assim_RuBP_sha_ebftemp        (:) !5
   real(r8), allocatable :: a_assim_RuBP_sha_dbftrop        (:) !6
   real(r8), allocatable :: a_assim_RuBP_sha_dbftemp        (:) !7
   real(r8), allocatable :: a_assim_RuBP_sha_dbfboreal      (:) !8
   real(r8), allocatable :: a_assim_RuBP_sha_ebstemp        (:) !9
   real(r8), allocatable :: a_assim_RuBP_sha_dbstemp        (:) !10
   real(r8), allocatable :: a_assim_RuBP_sha_dbsboreal      (:) !11
   real(r8), allocatable :: a_assim_RuBP_sha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_assim_RuBP_sha_c3grass        (:) !13
   real(r8), allocatable :: a_assim_RuBP_sha_c4grass        (:) !14
   real(r8), allocatable :: a_assim_Rubisco_sun_enftemp        (:) !1
   real(r8), allocatable :: a_assim_Rubisco_sun_enfboreal      (:) !2
   real(r8), allocatable :: a_assim_Rubisco_sun_dnfboreal      (:) !3
   real(r8), allocatable :: a_assim_Rubisco_sun_ebftrop        (:) !4
   real(r8), allocatable :: a_assim_Rubisco_sun_ebftemp        (:) !5
   real(r8), allocatable :: a_assim_Rubisco_sun_dbftrop        (:) !6
   real(r8), allocatable :: a_assim_Rubisco_sun_dbftemp        (:) !7
   real(r8), allocatable :: a_assim_Rubisco_sun_dbfboreal      (:) !8
   real(r8), allocatable :: a_assim_Rubisco_sun_ebstemp        (:) !9
   real(r8), allocatable :: a_assim_Rubisco_sun_dbstemp        (:) !10
   real(r8), allocatable :: a_assim_Rubisco_sun_dbsboreal      (:) !11
   real(r8), allocatable :: a_assim_Rubisco_sun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_assim_Rubisco_sun_c3grass        (:) !13
   real(r8), allocatable :: a_assim_Rubisco_sun_c4grass        (:) !14
   real(r8), allocatable :: a_assim_Rubisco_sha_enftemp        (:) !1
   real(r8), allocatable :: a_assim_Rubisco_sha_enfboreal      (:) !2
   real(r8), allocatable :: a_assim_Rubisco_sha_dnfboreal      (:) !3
   real(r8), allocatable :: a_assim_Rubisco_sha_ebftrop        (:) !4
   real(r8), allocatable :: a_assim_Rubisco_sha_ebftemp        (:) !5
   real(r8), allocatable :: a_assim_Rubisco_sha_dbftrop        (:) !6
   real(r8), allocatable :: a_assim_Rubisco_sha_dbftemp        (:) !7
   real(r8), allocatable :: a_assim_Rubisco_sha_dbfboreal      (:) !8
   real(r8), allocatable :: a_assim_Rubisco_sha_ebstemp        (:) !9
   real(r8), allocatable :: a_assim_Rubisco_sha_dbstemp        (:) !10
   real(r8), allocatable :: a_assim_Rubisco_sha_dbsboreal      (:) !11
   real(r8), allocatable :: a_assim_Rubisco_sha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_assim_Rubisco_sha_c3grass        (:) !13
   real(r8), allocatable :: a_assim_Rubisco_sha_c4grass        (:) !14
   real(r8), allocatable :: a_assimsun_enftemp        (:) !1
   real(r8), allocatable :: a_assimsun_enfboreal      (:) !2
   real(r8), allocatable :: a_assimsun_dnfboreal      (:) !3
   real(r8), allocatable :: a_assimsun_ebftrop        (:) !4
   real(r8), allocatable :: a_assimsun_ebftemp        (:) !5
   real(r8), allocatable :: a_assimsun_dbftrop        (:) !6
   real(r8), allocatable :: a_assimsun_dbftemp        (:) !7
   real(r8), allocatable :: a_assimsun_dbfboreal      (:) !8
   real(r8), allocatable :: a_assimsun_ebstemp        (:) !9
   real(r8), allocatable :: a_assimsun_dbstemp        (:) !10
   real(r8), allocatable :: a_assimsun_dbsboreal      (:) !11
   real(r8), allocatable :: a_assimsun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_assimsun_c3grass        (:) !13
   real(r8), allocatable :: a_assimsun_c4grass        (:) !14
   real(r8), allocatable :: a_assimsha_enftemp        (:) !1
   real(r8), allocatable :: a_assimsha_enfboreal      (:) !2
   real(r8), allocatable :: a_assimsha_dnfboreal      (:) !3
   real(r8), allocatable :: a_assimsha_ebftrop        (:) !4
   real(r8), allocatable :: a_assimsha_ebftemp        (:) !5
   real(r8), allocatable :: a_assimsha_dbftrop        (:) !6
   real(r8), allocatable :: a_assimsha_dbftemp        (:) !7
   real(r8), allocatable :: a_assimsha_dbfboreal      (:) !8
   real(r8), allocatable :: a_assimsha_ebstemp        (:) !9
   real(r8), allocatable :: a_assimsha_dbstemp        (:) !10
   real(r8), allocatable :: a_assimsha_dbsboreal      (:) !11
   real(r8), allocatable :: a_assimsha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_assimsha_c3grass        (:) !13
   real(r8), allocatable :: a_assimsha_c4grass        (:) !14
   real(r8), allocatable :: a_etrsun_enftemp        (:) !1
   real(r8), allocatable :: a_etrsun_enfboreal      (:) !2
   real(r8), allocatable :: a_etrsun_dnfboreal      (:) !3
   real(r8), allocatable :: a_etrsun_ebftrop        (:) !4
   real(r8), allocatable :: a_etrsun_ebftemp        (:) !5
   real(r8), allocatable :: a_etrsun_dbftrop        (:) !6
   real(r8), allocatable :: a_etrsun_dbftemp        (:) !7
   real(r8), allocatable :: a_etrsun_dbfboreal      (:) !8
   real(r8), allocatable :: a_etrsun_ebstemp        (:) !9
   real(r8), allocatable :: a_etrsun_dbstemp        (:) !10
   real(r8), allocatable :: a_etrsun_dbsboreal      (:) !11
   real(r8), allocatable :: a_etrsun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_etrsun_c3grass        (:) !13
   real(r8), allocatable :: a_etrsun_c4grass        (:) !14
   real(r8), allocatable :: a_etrsha_enftemp        (:) !1
   real(r8), allocatable :: a_etrsha_enfboreal      (:) !2
   real(r8), allocatable :: a_etrsha_dnfboreal      (:) !3
   real(r8), allocatable :: a_etrsha_ebftrop        (:) !4
   real(r8), allocatable :: a_etrsha_ebftemp        (:) !5
   real(r8), allocatable :: a_etrsha_dbftrop        (:) !6
   real(r8), allocatable :: a_etrsha_dbftemp        (:) !7
   real(r8), allocatable :: a_etrsha_dbfboreal      (:) !8
   real(r8), allocatable :: a_etrsha_ebstemp        (:) !9
   real(r8), allocatable :: a_etrsha_dbstemp        (:) !10
   real(r8), allocatable :: a_etrsha_dbsboreal      (:) !11
   real(r8), allocatable :: a_etrsha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_etrsha_c3grass        (:) !13
   real(r8), allocatable :: a_etrsha_c4grass        (:) !14
   real(r8), allocatable :: a_cisun_enftemp        (:) !1
   real(r8), allocatable :: a_cisun_enfboreal      (:) !2
   real(r8), allocatable :: a_cisun_dnfboreal      (:) !3
   real(r8), allocatable :: a_cisun_ebftrop        (:) !4
   real(r8), allocatable :: a_cisun_ebftemp        (:) !5
   real(r8), allocatable :: a_cisun_dbftrop        (:) !6
   real(r8), allocatable :: a_cisun_dbftemp        (:) !7
   real(r8), allocatable :: a_cisun_dbfboreal      (:) !8
   real(r8), allocatable :: a_cisun_ebstemp        (:) !9
   real(r8), allocatable :: a_cisun_dbstemp        (:) !10
   real(r8), allocatable :: a_cisun_dbsboreal      (:) !11
   real(r8), allocatable :: a_cisun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_cisun_c3grass        (:) !13
   real(r8), allocatable :: a_cisun_c4grass        (:) !14
   real(r8), allocatable :: a_cisha_enftemp        (:) !1
   real(r8), allocatable :: a_cisha_enfboreal      (:) !2
   real(r8), allocatable :: a_cisha_dnfboreal      (:) !3
   real(r8), allocatable :: a_cisha_ebftrop        (:) !4
   real(r8), allocatable :: a_cisha_ebftemp        (:) !5
   real(r8), allocatable :: a_cisha_dbftrop        (:) !6
   real(r8), allocatable :: a_cisha_dbftemp        (:) !7
   real(r8), allocatable :: a_cisha_dbfboreal      (:) !8
   real(r8), allocatable :: a_cisha_ebstemp        (:) !9
   real(r8), allocatable :: a_cisha_dbstemp        (:) !10
   real(r8), allocatable :: a_cisha_dbsboreal      (:) !11
   real(r8), allocatable :: a_cisha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_cisha_c3grass        (:) !13
   real(r8), allocatable :: a_cisha_c4grass        (:) !14
   real(r8), allocatable :: a_essun_enftemp        (:) !1
   real(r8), allocatable :: a_essun_enfboreal      (:) !2
   real(r8), allocatable :: a_essun_dnfboreal      (:) !3
   real(r8), allocatable :: a_essun_ebftrop        (:) !4
   real(r8), allocatable :: a_essun_ebftemp        (:) !5
   real(r8), allocatable :: a_essun_dbftrop        (:) !6
   real(r8), allocatable :: a_essun_dbftemp        (:) !7
   real(r8), allocatable :: a_essun_dbfboreal      (:) !8
   real(r8), allocatable :: a_essun_ebstemp        (:) !9
   real(r8), allocatable :: a_essun_dbstemp        (:) !10
   real(r8), allocatable :: a_essun_dbsboreal      (:) !11
   real(r8), allocatable :: a_essun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_essun_c3grass        (:) !13
   real(r8), allocatable :: a_essun_c4grass        (:) !14
   real(r8), allocatable :: a_essha_enftemp        (:) !1
   real(r8), allocatable :: a_essha_enfboreal      (:) !2
   real(r8), allocatable :: a_essha_dnfboreal      (:) !3
   real(r8), allocatable :: a_essha_ebftrop        (:) !4
   real(r8), allocatable :: a_essha_ebftemp        (:) !5
   real(r8), allocatable :: a_essha_dbftrop        (:) !6
   real(r8), allocatable :: a_essha_dbftemp        (:) !7
   real(r8), allocatable :: a_essha_dbfboreal      (:) !8
   real(r8), allocatable :: a_essha_ebstemp        (:) !9
   real(r8), allocatable :: a_essha_dbstemp        (:) !10
   real(r8), allocatable :: a_essha_dbsboreal      (:) !11
   real(r8), allocatable :: a_essha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_essha_c3grass        (:) !13
   real(r8), allocatable :: a_essha_c4grass        (:) !14
   real(r8), allocatable :: a_gssun_enftemp        (:) !1
   real(r8), allocatable :: a_gssun_enfboreal      (:) !2
   real(r8), allocatable :: a_gssun_dnfboreal      (:) !3
   real(r8), allocatable :: a_gssun_ebftrop        (:) !4
   real(r8), allocatable :: a_gssun_ebftemp        (:) !5
   real(r8), allocatable :: a_gssun_dbftrop        (:) !6
   real(r8), allocatable :: a_gssun_dbftemp        (:) !7
   real(r8), allocatable :: a_gssun_dbfboreal      (:) !8
   real(r8), allocatable :: a_gssun_ebstemp        (:) !9
   real(r8), allocatable :: a_gssun_dbstemp        (:) !10
   real(r8), allocatable :: a_gssun_dbsboreal      (:) !11
   real(r8), allocatable :: a_gssun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_gssun_c3grass        (:) !13
   real(r8), allocatable :: a_gssun_c4grass        (:) !14
   real(r8), allocatable :: a_gssha_enftemp        (:) !1
   real(r8), allocatable :: a_gssha_enfboreal      (:) !2
   real(r8), allocatable :: a_gssha_dnfboreal      (:) !3
   real(r8), allocatable :: a_gssha_ebftrop        (:) !4
   real(r8), allocatable :: a_gssha_ebftemp        (:) !5
   real(r8), allocatable :: a_gssha_dbftrop        (:) !6
   real(r8), allocatable :: a_gssha_dbftemp        (:) !7
   real(r8), allocatable :: a_gssha_dbfboreal      (:) !8
   real(r8), allocatable :: a_gssha_ebstemp        (:) !9
   real(r8), allocatable :: a_gssha_dbstemp        (:) !10
   real(r8), allocatable :: a_gssha_dbsboreal      (:) !11
   real(r8), allocatable :: a_gssha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_gssha_c3grass        (:) !13
   real(r8), allocatable :: a_gssha_c4grass        (:) !14
   real(r8), allocatable :: a_gammasun_enftemp        (:) !1
   real(r8), allocatable :: a_gammasun_enfboreal      (:) !2
   real(r8), allocatable :: a_gammasun_dnfboreal      (:) !3
   real(r8), allocatable :: a_gammasun_ebftrop        (:) !4
   real(r8), allocatable :: a_gammasun_ebftemp        (:) !5
   real(r8), allocatable :: a_gammasun_dbftrop        (:) !6
   real(r8), allocatable :: a_gammasun_dbftemp        (:) !7
   real(r8), allocatable :: a_gammasun_dbfboreal      (:) !8
   real(r8), allocatable :: a_gammasun_ebstemp        (:) !9
   real(r8), allocatable :: a_gammasun_dbstemp        (:) !10
   real(r8), allocatable :: a_gammasun_dbsboreal      (:) !11
   real(r8), allocatable :: a_gammasun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_gammasun_c3grass        (:) !13
   real(r8), allocatable :: a_gammasun_c4grass        (:) !14
   real(r8), allocatable :: a_gammasha_enftemp        (:) !1
   real(r8), allocatable :: a_gammasha_enfboreal      (:) !2
   real(r8), allocatable :: a_gammasha_dnfboreal      (:) !3
   real(r8), allocatable :: a_gammasha_ebftrop        (:) !4
   real(r8), allocatable :: a_gammasha_ebftemp        (:) !5
   real(r8), allocatable :: a_gammasha_dbftrop        (:) !6
   real(r8), allocatable :: a_gammasha_dbftemp        (:) !7
   real(r8), allocatable :: a_gammasha_dbfboreal      (:) !8
   real(r8), allocatable :: a_gammasha_ebstemp        (:) !9
   real(r8), allocatable :: a_gammasha_dbstemp        (:) !10
   real(r8), allocatable :: a_gammasha_dbsboreal      (:) !11
   real(r8), allocatable :: a_gammasha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_gammasha_c3grass        (:) !13
   real(r8), allocatable :: a_gammasha_c4grass        (:) !14
   real(r8), allocatable :: a_RuBPlimitfrac_sun_enftemp        (:) !1
   real(r8), allocatable :: a_RuBPlimitfrac_sun_enfboreal      (:) !2
   real(r8), allocatable :: a_RuBPlimitfrac_sun_dnfboreal      (:) !3
   real(r8), allocatable :: a_RuBPlimitfrac_sun_ebftrop        (:) !4
   real(r8), allocatable :: a_RuBPlimitfrac_sun_ebftemp        (:) !5
   real(r8), allocatable :: a_RuBPlimitfrac_sun_dbftrop        (:) !6
   real(r8), allocatable :: a_RuBPlimitfrac_sun_dbftemp        (:) !7
   real(r8), allocatable :: a_RuBPlimitfrac_sun_dbfboreal      (:) !8
   real(r8), allocatable :: a_RuBPlimitfrac_sun_ebstemp        (:) !9
   real(r8), allocatable :: a_RuBPlimitfrac_sun_dbstemp        (:) !10
   real(r8), allocatable :: a_RuBPlimitfrac_sun_dbsboreal      (:) !11
   real(r8), allocatable :: a_RuBPlimitfrac_sun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_RuBPlimitfrac_sun_c3grass        (:) !13
   real(r8), allocatable :: a_RuBPlimitfrac_sun_c4grass        (:) !14
   real(r8), allocatable :: a_RuBPlimitfrac_sha_enftemp        (:) !1
   real(r8), allocatable :: a_RuBPlimitfrac_sha_enfboreal      (:) !2
   real(r8), allocatable :: a_RuBPlimitfrac_sha_dnfboreal      (:) !3
   real(r8), allocatable :: a_RuBPlimitfrac_sha_ebftrop        (:) !4
   real(r8), allocatable :: a_RuBPlimitfrac_sha_ebftemp        (:) !5
   real(r8), allocatable :: a_RuBPlimitfrac_sha_dbftrop        (:) !6
   real(r8), allocatable :: a_RuBPlimitfrac_sha_dbftemp        (:) !7
   real(r8), allocatable :: a_RuBPlimitfrac_sha_dbfboreal      (:) !8
   real(r8), allocatable :: a_RuBPlimitfrac_sha_ebstemp        (:) !9
   real(r8), allocatable :: a_RuBPlimitfrac_sha_dbstemp        (:) !10
   real(r8), allocatable :: a_RuBPlimitfrac_sha_dbsboreal      (:) !11
   real(r8), allocatable :: a_RuBPlimitfrac_sha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_RuBPlimitfrac_sha_c3grass        (:) !13
   real(r8), allocatable :: a_RuBPlimitfrac_sha_c4grass        (:) !14
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_enftemp        (:) !1
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_enfboreal      (:) !2
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_dnfboreal      (:) !3
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_ebftrop        (:) !4
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_ebftemp        (:) !5
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_dbftrop        (:) !6
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_dbftemp        (:) !7
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_dbfboreal      (:) !8
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_ebstemp        (:) !9
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_dbstemp        (:) !10
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_dbsboreal      (:) !11
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_c3grass        (:) !13
   real(r8), allocatable :: a_Rubiscolimitfrac_sun_c4grass        (:) !14
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_enftemp        (:) !1
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_enfboreal      (:) !2
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_dnfboreal      (:) !3
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_ebftrop        (:) !4
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_ebftemp        (:) !5
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_dbftrop        (:) !6
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_dbftemp        (:) !7
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_dbfboreal      (:) !8
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_ebstemp        (:) !9
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_dbstemp        (:) !10
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_dbsboreal      (:) !11
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_c3grass        (:) !13
   real(r8), allocatable :: a_Rubiscolimitfrac_sha_c4grass        (:) !14
   real(r8), allocatable :: a_Sinklimitfrac_sun_enftemp        (:) !1
   real(r8), allocatable :: a_Sinklimitfrac_sun_enfboreal      (:) !2
   real(r8), allocatable :: a_Sinklimitfrac_sun_dnfboreal      (:) !3
   real(r8), allocatable :: a_Sinklimitfrac_sun_ebftrop        (:) !4
   real(r8), allocatable :: a_Sinklimitfrac_sun_ebftemp        (:) !5
   real(r8), allocatable :: a_Sinklimitfrac_sun_dbftrop        (:) !6
   real(r8), allocatable :: a_Sinklimitfrac_sun_dbftemp        (:) !7
   real(r8), allocatable :: a_Sinklimitfrac_sun_dbfboreal      (:) !8
   real(r8), allocatable :: a_Sinklimitfrac_sun_ebstemp        (:) !9
   real(r8), allocatable :: a_Sinklimitfrac_sun_dbstemp        (:) !10
   real(r8), allocatable :: a_Sinklimitfrac_sun_dbsboreal      (:) !11
   real(r8), allocatable :: a_Sinklimitfrac_sun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_Sinklimitfrac_sun_c3grass        (:) !13
   real(r8), allocatable :: a_Sinklimitfrac_sun_c4grass        (:) !14
   real(r8), allocatable :: a_Sinklimitfrac_sha_enftemp        (:) !1
   real(r8), allocatable :: a_Sinklimitfrac_sha_enfboreal      (:) !2
   real(r8), allocatable :: a_Sinklimitfrac_sha_dnfboreal      (:) !3
   real(r8), allocatable :: a_Sinklimitfrac_sha_ebftrop        (:) !4
   real(r8), allocatable :: a_Sinklimitfrac_sha_ebftemp        (:) !5
   real(r8), allocatable :: a_Sinklimitfrac_sha_dbftrop        (:) !6
   real(r8), allocatable :: a_Sinklimitfrac_sha_dbftemp        (:) !7
   real(r8), allocatable :: a_Sinklimitfrac_sha_dbfboreal      (:) !8
   real(r8), allocatable :: a_Sinklimitfrac_sha_ebstemp        (:) !9
   real(r8), allocatable :: a_Sinklimitfrac_sha_dbstemp        (:) !10
   real(r8), allocatable :: a_Sinklimitfrac_sha_dbsboreal      (:) !11
   real(r8), allocatable :: a_Sinklimitfrac_sha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_Sinklimitfrac_sha_c3grass        (:) !13
   real(r8), allocatable :: a_Sinklimitfrac_sha_c4grass        (:) !14
   real(r8), allocatable :: a_rstfacsun_enftemp        (:) !1
   real(r8), allocatable :: a_rstfacsun_enfboreal      (:) !2
   real(r8), allocatable :: a_rstfacsun_dnfboreal      (:) !3
   real(r8), allocatable :: a_rstfacsun_ebftrop        (:) !4
   real(r8), allocatable :: a_rstfacsun_ebftemp        (:) !5
   real(r8), allocatable :: a_rstfacsun_dbftrop        (:) !6
   real(r8), allocatable :: a_rstfacsun_dbftemp        (:) !7
   real(r8), allocatable :: a_rstfacsun_dbfboreal      (:) !8
   real(r8), allocatable :: a_rstfacsun_ebstemp        (:) !9
   real(r8), allocatable :: a_rstfacsun_dbstemp        (:) !10
   real(r8), allocatable :: a_rstfacsun_dbsboreal      (:) !11
   real(r8), allocatable :: a_rstfacsun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_rstfacsun_c3grass        (:) !13
   real(r8), allocatable :: a_rstfacsun_c4grass        (:) !14
   real(r8), allocatable :: a_rstfacsha_enftemp        (:) !1
   real(r8), allocatable :: a_rstfacsha_enfboreal      (:) !2
   real(r8), allocatable :: a_rstfacsha_dnfboreal      (:) !3
   real(r8), allocatable :: a_rstfacsha_ebftrop        (:) !4
   real(r8), allocatable :: a_rstfacsha_ebftemp        (:) !5
   real(r8), allocatable :: a_rstfacsha_dbftrop        (:) !6
   real(r8), allocatable :: a_rstfacsha_dbftemp        (:) !7
   real(r8), allocatable :: a_rstfacsha_dbfboreal      (:) !8
   real(r8), allocatable :: a_rstfacsha_ebstemp        (:) !9
   real(r8), allocatable :: a_rstfacsha_dbstemp        (:) !10
   real(r8), allocatable :: a_rstfacsha_dbsboreal      (:) !11
   real(r8), allocatable :: a_rstfacsha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_rstfacsha_c3grass        (:) !13
   real(r8), allocatable :: a_rstfacsha_c4grass        (:) !14
   real(r8), allocatable :: a_lambdasun_enftemp        (:) !1
   real(r8), allocatable :: a_lambdasun_enfboreal      (:) !2
   real(r8), allocatable :: a_lambdasun_dnfboreal      (:) !3
   real(r8), allocatable :: a_lambdasun_ebftrop        (:) !4
   real(r8), allocatable :: a_lambdasun_ebftemp        (:) !5
   real(r8), allocatable :: a_lambdasun_dbftrop        (:) !6
   real(r8), allocatable :: a_lambdasun_dbftemp        (:) !7
   real(r8), allocatable :: a_lambdasun_dbfboreal      (:) !8
   real(r8), allocatable :: a_lambdasun_ebstemp        (:) !9
   real(r8), allocatable :: a_lambdasun_dbstemp        (:) !10
   real(r8), allocatable :: a_lambdasun_dbsboreal      (:) !11
   real(r8), allocatable :: a_lambdasun_c3arcgrass     (:) !12
   real(r8), allocatable :: a_lambdasun_c3grass        (:) !13
   real(r8), allocatable :: a_lambdasun_c4grass        (:) !14
   real(r8), allocatable :: a_lambdasha_enftemp        (:) !1
   real(r8), allocatable :: a_lambdasha_enfboreal      (:) !2
   real(r8), allocatable :: a_lambdasha_dnfboreal      (:) !3
   real(r8), allocatable :: a_lambdasha_ebftrop        (:) !4
   real(r8), allocatable :: a_lambdasha_ebftemp        (:) !5
   real(r8), allocatable :: a_lambdasha_dbftrop        (:) !6
   real(r8), allocatable :: a_lambdasha_dbftemp        (:) !7
   real(r8), allocatable :: a_lambdasha_dbfboreal      (:) !8
   real(r8), allocatable :: a_lambdasha_ebstemp        (:) !9
   real(r8), allocatable :: a_lambdasha_dbstemp        (:) !10
   real(r8), allocatable :: a_lambdasha_dbsboreal      (:) !11
   real(r8), allocatable :: a_lambdasha_c3arcgrass     (:) !12
   real(r8), allocatable :: a_lambdasha_c3grass        (:) !13
   real(r8), allocatable :: a_lambdasha_c4grass        (:) !14
   real(r8), allocatable :: a_lambda                   (:) !14
#endif
#endif
#ifdef NITRIF
   real(r8), allocatable :: a_O2_DECOMP_DEPTH_UNSAT (:,:)
   real(r8), allocatable :: a_CONC_O2_UNSAT         (:,:)
#endif
#ifdef CROP
   real(r8), allocatable :: a_pdcorn                (:)
   real(r8), allocatable :: a_pdswheat              (:)
   real(r8), allocatable :: a_pdwwheat              (:)
   real(r8), allocatable :: a_pdsoybean             (:)
   real(r8), allocatable :: a_pdcotton              (:)
   real(r8), allocatable :: a_pdrice1               (:)
   real(r8), allocatable :: a_pdrice2               (:)
   real(r8), allocatable :: a_pdsugarcane           (:)
   real(r8), allocatable :: a_plantdate             (:)
   real(r8), allocatable :: a_fertnitro_corn        (:)
   real(r8), allocatable :: a_fertnitro_swheat      (:)
   real(r8), allocatable :: a_fertnitro_wwheat      (:)
   real(r8), allocatable :: a_fertnitro_soybean     (:)
   real(r8), allocatable :: a_fertnitro_cotton      (:)
   real(r8), allocatable :: a_fertnitro_rice1       (:)
   real(r8), allocatable :: a_fertnitro_rice2       (:)
   real(r8), allocatable :: a_fertnitro_sugarcane   (:)
   real(r8), allocatable :: a_cphase             (:)
   real(r8), allocatable :: a_gddplant           (:)
   real(r8), allocatable :: a_gddmaturity        (:)
   real(r8), allocatable :: a_vf                 (:)
   real(r8), allocatable :: a_hui                (:)
   real(r8), allocatable :: a_cropprod1c         (:)
   real(r8), allocatable :: a_cropprod1c_loss    (:)
   real(r8), allocatable :: a_cropseedc_deficit  (:)
   real(r8), allocatable :: a_grainc_to_cropprodc(:)
   real(r8), allocatable :: a_grainc_to_seed     (:)
   real(r8), allocatable :: a_fert_to_sminn      (:)
#endif
   real(r8), allocatable :: a_ndep_to_sminn      (:)
#ifdef Fire
   real(r8), allocatable :: a_abm                (:)
   real(r8), allocatable :: a_gdp                (:)
   real(r8), allocatable :: a_peatf              (:)
   real(r8), allocatable :: a_hdm                (:)
   real(r8), allocatable :: a_lnfm               (:)
#endif
#endif
#ifdef OzoneStress
   real(r8), allocatable :: a_ozone              (:)
#endif

   real(r8), allocatable :: a_t_soisno    (:,:)    
   real(r8), allocatable :: a_wliq_soisno (:,:)
   real(r8), allocatable :: a_wice_soisno (:,:)
   real(r8), allocatable :: a_h2osoi      (:,:)
   real(r8), allocatable :: a_rootr       (:,:)
   real(r8), allocatable :: a_BD_all      (:,:)
   real(r8), allocatable :: a_wfc         (:,:)
   real(r8), allocatable :: a_OM_density  (:,:)
#ifdef PLANT_HYDRAULIC_STRESS
   real(r8), allocatable :: a_vegwp       (:,:)
#endif
   real(r8), allocatable :: a_t_lake      (:,:) 
   real(r8), allocatable :: a_lake_icefrac(:,:) 

#ifdef BGC
   real(r8), allocatable :: a_litr1c_vr   (:,:)
   real(r8), allocatable :: a_litr2c_vr   (:,:)
   real(r8), allocatable :: a_litr3c_vr   (:,:)
   real(r8), allocatable :: a_soil1c_vr   (:,:)
   real(r8), allocatable :: a_soil2c_vr   (:,:)
   real(r8), allocatable :: a_soil3c_vr   (:,:)
   real(r8), allocatable :: a_cwdc_vr     (:,:)
   real(r8), allocatable :: a_litr1n_vr   (:,:)
   real(r8), allocatable :: a_litr2n_vr   (:,:)
   real(r8), allocatable :: a_litr3n_vr   (:,:)
   real(r8), allocatable :: a_soil1n_vr   (:,:)
   real(r8), allocatable :: a_soil2n_vr   (:,:)
   real(r8), allocatable :: a_soil3n_vr   (:,:)
   real(r8), allocatable :: a_cwdn_vr     (:,:)
   real(r8), allocatable :: a_sminn_vr    (:,:)
   real(r8), allocatable :: decomp_vr_tmp (:,:)
#endif

   real(r8), allocatable :: a_ustar(:) 
   real(r8), allocatable :: a_tstar(:)
   real(r8), allocatable :: a_qstar(:)
   real(r8), allocatable :: a_zol  (:)
   real(r8), allocatable :: a_rib  (:)
   real(r8), allocatable :: a_fm   (:)
   real(r8), allocatable :: a_fh   (:)
   real(r8), allocatable :: a_fq   (:)

   real(r8), allocatable :: a_us10m(:) 
   real(r8), allocatable :: a_vs10m(:) 
   real(r8), allocatable :: a_fm10m(:) 

   real(r8), allocatable :: a_sr     (:)
   real(r8), allocatable :: a_solvd  (:)
   real(r8), allocatable :: a_solvi  (:)
   real(r8), allocatable :: a_solnd  (:)
   real(r8), allocatable :: a_solni  (:)
   real(r8), allocatable :: a_srvd   (:)
   real(r8), allocatable :: a_srvi   (:)
   real(r8), allocatable :: a_srnd   (:)
   real(r8), allocatable :: a_srni   (:)
   real(r8), allocatable :: a_solvdln(:)
   real(r8), allocatable :: a_solviln(:)
   real(r8), allocatable :: a_solndln(:)
   real(r8), allocatable :: a_solniln(:)
   real(r8), allocatable :: a_srvdln (:)
   real(r8), allocatable :: a_srviln (:)
   real(r8), allocatable :: a_srndln (:)
   real(r8), allocatable :: a_srniln (:)

   public :: allocate_acc_fluxes
   public :: deallocate_acc_fluxes
   public :: flush_acc_fluxes
   public :: accumulate_fluxes

contains

   subroutine allocate_acc_fluxes 

      use spmd_task
      USE GlobalVars
      use mod_landpatch, only : numpatch
      implicit none

      if (p_is_worker) then
         if (numpatch > 0) then

            allocate (a_us     (numpatch))  
            allocate (a_vs     (numpatch))
            allocate (a_t      (numpatch))
            allocate (a_q      (numpatch))
            allocate (a_prc    (numpatch))
            allocate (a_prl    (numpatch))
            allocate (a_pbot   (numpatch))
            allocate (a_frl    (numpatch))
            allocate (a_solarin(numpatch))

            allocate (a_taux      (numpatch))
            allocate (a_tauy      (numpatch))
            allocate (a_fsena     (numpatch))
            allocate (a_lfevpa    (numpatch))
            allocate (a_fevpa     (numpatch))
            allocate (a_fsenl     (numpatch))
            allocate (a_fevpl     (numpatch))
            allocate (a_etr       (numpatch))
            allocate (a_fseng     (numpatch))
            allocate (a_fevpg     (numpatch))
            allocate (a_fgrnd     (numpatch))
            allocate (a_sabvsun   (numpatch))  
            allocate (a_sabvsha   (numpatch)) 
            allocate (a_sabg      (numpatch))
            allocate (a_olrg      (numpatch))
            allocate (a_rnet      (numpatch))
            allocate (a_xerr      (numpatch))
            allocate (a_zerr      (numpatch))
            allocate (a_rsur      (numpatch))
            allocate (a_rnof      (numpatch))
            allocate (a_qintr     (numpatch))
            allocate (a_qinfl     (numpatch))
            allocate (a_qdrip     (numpatch))
            allocate (a_rstfacsun (numpatch))
            allocate (a_rstfacsha (numpatch))
            allocate (a_gs_sun     (numpatch))
            allocate (a_gs_sha     (numpatch))
            allocate (a_dpond     (numpatch))

            allocate (a_zwt       (numpatch))
            allocate (a_wa        (numpatch))
            allocate (a_wat       (numpatch))
            allocate (a_assim     (numpatch))
            allocate (a_respc     (numpatch))

#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
            allocate (a_assim_RuBP_sun_enftemp        (numpatch)) !1
            allocate (a_assim_RuBP_sun_enfboreal      (numpatch)) !2
            allocate (a_assim_RuBP_sun_dnfboreal      (numpatch)) !3
            allocate (a_assim_RuBP_sun_ebftrop        (numpatch)) !4
            allocate (a_assim_RuBP_sun_ebftemp        (numpatch)) !5
            allocate (a_assim_RuBP_sun_dbftrop        (numpatch)) !6
            allocate (a_assim_RuBP_sun_dbftemp        (numpatch)) !7
            allocate (a_assim_RuBP_sun_dbfboreal      (numpatch)) !8
            allocate (a_assim_RuBP_sun_ebstemp        (numpatch)) !9
            allocate (a_assim_RuBP_sun_dbstemp        (numpatch)) !10
            allocate (a_assim_RuBP_sun_dbsboreal      (numpatch)) !11
            allocate (a_assim_RuBP_sun_c3arcgrass     (numpatch)) !12
            allocate (a_assim_RuBP_sun_c3grass        (numpatch)) !13
            allocate (a_assim_RuBP_sun_c4grass        (numpatch)) !14
            allocate (a_assim_RuBP_sha_enftemp        (numpatch)) !1
            allocate (a_assim_RuBP_sha_enfboreal      (numpatch)) !2
            allocate (a_assim_RuBP_sha_dnfboreal      (numpatch)) !3
            allocate (a_assim_RuBP_sha_ebftrop        (numpatch)) !4
            allocate (a_assim_RuBP_sha_ebftemp        (numpatch)) !5
            allocate (a_assim_RuBP_sha_dbftrop        (numpatch)) !6
            allocate (a_assim_RuBP_sha_dbftemp        (numpatch)) !7
            allocate (a_assim_RuBP_sha_dbfboreal      (numpatch)) !8
            allocate (a_assim_RuBP_sha_ebstemp        (numpatch)) !9
            allocate (a_assim_RuBP_sha_dbstemp        (numpatch)) !10
            allocate (a_assim_RuBP_sha_dbsboreal      (numpatch)) !11
            allocate (a_assim_RuBP_sha_c3arcgrass     (numpatch)) !12
            allocate (a_assim_RuBP_sha_c3grass        (numpatch)) !13
            allocate (a_assim_RuBP_sha_c4grass        (numpatch)) !14
            allocate (a_assim_Rubisco_sun_enftemp        (numpatch)) !1
            allocate (a_assim_Rubisco_sun_enfboreal      (numpatch)) !2
            allocate (a_assim_Rubisco_sun_dnfboreal      (numpatch)) !3
            allocate (a_assim_Rubisco_sun_ebftrop        (numpatch)) !4
            allocate (a_assim_Rubisco_sun_ebftemp        (numpatch)) !5
            allocate (a_assim_Rubisco_sun_dbftrop        (numpatch)) !6
            allocate (a_assim_Rubisco_sun_dbftemp        (numpatch)) !7
            allocate (a_assim_Rubisco_sun_dbfboreal      (numpatch)) !8
            allocate (a_assim_Rubisco_sun_ebstemp        (numpatch)) !9
            allocate (a_assim_Rubisco_sun_dbstemp        (numpatch)) !10
            allocate (a_assim_Rubisco_sun_dbsboreal      (numpatch)) !11
            allocate (a_assim_Rubisco_sun_c3arcgrass     (numpatch)) !12
            allocate (a_assim_Rubisco_sun_c3grass        (numpatch)) !13
            allocate (a_assim_Rubisco_sun_c4grass        (numpatch)) !14
            allocate (a_assim_Rubisco_sha_enftemp        (numpatch)) !1
            allocate (a_assim_Rubisco_sha_enfboreal      (numpatch)) !2
            allocate (a_assim_Rubisco_sha_dnfboreal      (numpatch)) !3
            allocate (a_assim_Rubisco_sha_ebftrop        (numpatch)) !4
            allocate (a_assim_Rubisco_sha_ebftemp        (numpatch)) !5
            allocate (a_assim_Rubisco_sha_dbftrop        (numpatch)) !6
            allocate (a_assim_Rubisco_sha_dbftemp        (numpatch)) !7
            allocate (a_assim_Rubisco_sha_dbfboreal      (numpatch)) !8
            allocate (a_assim_Rubisco_sha_ebstemp        (numpatch)) !9
            allocate (a_assim_Rubisco_sha_dbstemp        (numpatch)) !10
            allocate (a_assim_Rubisco_sha_dbsboreal      (numpatch)) !11
            allocate (a_assim_Rubisco_sha_c3arcgrass     (numpatch)) !12
            allocate (a_assim_Rubisco_sha_c3grass        (numpatch)) !13
            allocate (a_assim_Rubisco_sha_c4grass        (numpatch)) !14
            allocate (a_assimsun_enftemp        (numpatch)) !1
            allocate (a_assimsun_enfboreal      (numpatch)) !2
            allocate (a_assimsun_dnfboreal      (numpatch)) !3
            allocate (a_assimsun_ebftrop        (numpatch)) !4
            allocate (a_assimsun_ebftemp        (numpatch)) !5
            allocate (a_assimsun_dbftrop        (numpatch)) !6
            allocate (a_assimsun_dbftemp        (numpatch)) !7
            allocate (a_assimsun_dbfboreal      (numpatch)) !8
            allocate (a_assimsun_ebstemp        (numpatch)) !9
            allocate (a_assimsun_dbstemp        (numpatch)) !10
            allocate (a_assimsun_dbsboreal      (numpatch)) !11
            allocate (a_assimsun_c3arcgrass     (numpatch)) !12
            allocate (a_assimsun_c3grass        (numpatch)) !13
            allocate (a_assimsun_c4grass        (numpatch)) !14
            allocate (a_assimsha_enftemp        (numpatch)) !1
            allocate (a_assimsha_enfboreal      (numpatch)) !2
            allocate (a_assimsha_dnfboreal      (numpatch)) !3
            allocate (a_assimsha_ebftrop        (numpatch)) !4
            allocate (a_assimsha_ebftemp        (numpatch)) !5
            allocate (a_assimsha_dbftrop        (numpatch)) !6
            allocate (a_assimsha_dbftemp        (numpatch)) !7
            allocate (a_assimsha_dbfboreal      (numpatch)) !8
            allocate (a_assimsha_ebstemp        (numpatch)) !9
            allocate (a_assimsha_dbstemp        (numpatch)) !10
            allocate (a_assimsha_dbsboreal      (numpatch)) !11
            allocate (a_assimsha_c3arcgrass     (numpatch)) !12
            allocate (a_assimsha_c3grass        (numpatch)) !13
            allocate (a_assimsha_c4grass        (numpatch)) !14
            allocate (a_etrsun_enftemp        (numpatch)) !1
            allocate (a_etrsun_enfboreal      (numpatch)) !2
            allocate (a_etrsun_dnfboreal      (numpatch)) !3
            allocate (a_etrsun_ebftrop        (numpatch)) !4
            allocate (a_etrsun_ebftemp        (numpatch)) !5
            allocate (a_etrsun_dbftrop        (numpatch)) !6
            allocate (a_etrsun_dbftemp        (numpatch)) !7
            allocate (a_etrsun_dbfboreal      (numpatch)) !8
            allocate (a_etrsun_ebstemp        (numpatch)) !9
            allocate (a_etrsun_dbstemp        (numpatch)) !10
            allocate (a_etrsun_dbsboreal      (numpatch)) !11
            allocate (a_etrsun_c3arcgrass     (numpatch)) !12
            allocate (a_etrsun_c3grass        (numpatch)) !13
            allocate (a_etrsun_c4grass        (numpatch)) !14
            allocate (a_etrsha_enftemp        (numpatch)) !1
            allocate (a_etrsha_enfboreal      (numpatch)) !2
            allocate (a_etrsha_dnfboreal      (numpatch)) !3
            allocate (a_etrsha_ebftrop        (numpatch)) !4
            allocate (a_etrsha_ebftemp        (numpatch)) !5
            allocate (a_etrsha_dbftrop        (numpatch)) !6
            allocate (a_etrsha_dbftemp        (numpatch)) !7
            allocate (a_etrsha_dbfboreal      (numpatch)) !8
            allocate (a_etrsha_ebstemp        (numpatch)) !9
            allocate (a_etrsha_dbstemp        (numpatch)) !10
            allocate (a_etrsha_dbsboreal      (numpatch)) !11
            allocate (a_etrsha_c3arcgrass     (numpatch)) !12
            allocate (a_etrsha_c3grass        (numpatch)) !13
            allocate (a_etrsha_c4grass        (numpatch)) !14
            allocate (a_cisun_enftemp        (numpatch)) !1
            allocate (a_cisun_enfboreal      (numpatch)) !2
            allocate (a_cisun_dnfboreal      (numpatch)) !3
            allocate (a_cisun_ebftrop        (numpatch)) !4
            allocate (a_cisun_ebftemp        (numpatch)) !5
            allocate (a_cisun_dbftrop        (numpatch)) !6
            allocate (a_cisun_dbftemp        (numpatch)) !7
            allocate (a_cisun_dbfboreal      (numpatch)) !8
            allocate (a_cisun_ebstemp        (numpatch)) !9
            allocate (a_cisun_dbstemp        (numpatch)) !10
            allocate (a_cisun_dbsboreal      (numpatch)) !11
            allocate (a_cisun_c3arcgrass     (numpatch)) !12
            allocate (a_cisun_c3grass        (numpatch)) !13
            allocate (a_cisun_c4grass        (numpatch)) !14
            allocate (a_cisha_enftemp        (numpatch)) !1
            allocate (a_cisha_enfboreal      (numpatch)) !2
            allocate (a_cisha_dnfboreal      (numpatch)) !3
            allocate (a_cisha_ebftrop        (numpatch)) !4
            allocate (a_cisha_ebftemp        (numpatch)) !5
            allocate (a_cisha_dbftrop        (numpatch)) !6
            allocate (a_cisha_dbftemp        (numpatch)) !7
            allocate (a_cisha_dbfboreal      (numpatch)) !8
            allocate (a_cisha_ebstemp        (numpatch)) !9
            allocate (a_cisha_dbstemp        (numpatch)) !10
            allocate (a_cisha_dbsboreal      (numpatch)) !11
            allocate (a_cisha_c3arcgrass     (numpatch)) !12
            allocate (a_cisha_c3grass        (numpatch)) !13
            allocate (a_cisha_c4grass        (numpatch)) !14
            allocate (a_essun_enftemp        (numpatch)) !1
            allocate (a_essun_enfboreal      (numpatch)) !2
            allocate (a_essun_dnfboreal      (numpatch)) !3
            allocate (a_essun_ebftrop        (numpatch)) !4
            allocate (a_essun_ebftemp        (numpatch)) !5
            allocate (a_essun_dbftrop        (numpatch)) !6
            allocate (a_essun_dbftemp        (numpatch)) !7
            allocate (a_essun_dbfboreal      (numpatch)) !8
            allocate (a_essun_ebstemp        (numpatch)) !9
            allocate (a_essun_dbstemp        (numpatch)) !10
            allocate (a_essun_dbsboreal      (numpatch)) !11
            allocate (a_essun_c3arcgrass     (numpatch)) !12
            allocate (a_essun_c3grass        (numpatch)) !13
            allocate (a_essun_c4grass        (numpatch)) !14
            allocate (a_essha_enftemp        (numpatch)) !1
            allocate (a_essha_enfboreal      (numpatch)) !2
            allocate (a_essha_dnfboreal      (numpatch)) !3
            allocate (a_essha_ebftrop        (numpatch)) !4
            allocate (a_essha_ebftemp        (numpatch)) !5
            allocate (a_essha_dbftrop        (numpatch)) !6
            allocate (a_essha_dbftemp        (numpatch)) !7
            allocate (a_essha_dbfboreal      (numpatch)) !8
            allocate (a_essha_ebstemp        (numpatch)) !9
            allocate (a_essha_dbstemp        (numpatch)) !10
            allocate (a_essha_dbsboreal      (numpatch)) !11
            allocate (a_essha_c3arcgrass     (numpatch)) !12
            allocate (a_essha_c3grass        (numpatch)) !13
            allocate (a_essha_c4grass        (numpatch)) !14
            allocate (a_gssun_enftemp        (numpatch)) !1
            allocate (a_gssun_enfboreal      (numpatch)) !2
            allocate (a_gssun_dnfboreal      (numpatch)) !3
            allocate (a_gssun_ebftrop        (numpatch)) !4
            allocate (a_gssun_ebftemp        (numpatch)) !5
            allocate (a_gssun_dbftrop        (numpatch)) !6
            allocate (a_gssun_dbftemp        (numpatch)) !7
            allocate (a_gssun_dbfboreal      (numpatch)) !8
            allocate (a_gssun_ebstemp        (numpatch)) !9
            allocate (a_gssun_dbstemp        (numpatch)) !10
            allocate (a_gssun_dbsboreal      (numpatch)) !11
            allocate (a_gssun_c3arcgrass     (numpatch)) !12
            allocate (a_gssun_c3grass        (numpatch)) !13
            allocate (a_gssun_c4grass        (numpatch)) !14
            allocate (a_gssha_enftemp        (numpatch)) !1
            allocate (a_gssha_enfboreal      (numpatch)) !2
            allocate (a_gssha_dnfboreal      (numpatch)) !3
            allocate (a_gssha_ebftrop        (numpatch)) !4
            allocate (a_gssha_ebftemp        (numpatch)) !5
            allocate (a_gssha_dbftrop        (numpatch)) !6
            allocate (a_gssha_dbftemp        (numpatch)) !7
            allocate (a_gssha_dbfboreal      (numpatch)) !8
            allocate (a_gssha_ebstemp        (numpatch)) !9
            allocate (a_gssha_dbstemp        (numpatch)) !10
            allocate (a_gssha_dbsboreal      (numpatch)) !11
            allocate (a_gssha_c3arcgrass     (numpatch)) !12
            allocate (a_gssha_c3grass        (numpatch)) !13
            allocate (a_gssha_c4grass        (numpatch)) !14
            allocate (a_gammasun_enftemp        (numpatch)) !1
            allocate (a_gammasun_enfboreal      (numpatch)) !2
            allocate (a_gammasun_dnfboreal      (numpatch)) !3
            allocate (a_gammasun_ebftrop        (numpatch)) !4
            allocate (a_gammasun_ebftemp        (numpatch)) !5
            allocate (a_gammasun_dbftrop        (numpatch)) !6
            allocate (a_gammasun_dbftemp        (numpatch)) !7
            allocate (a_gammasun_dbfboreal      (numpatch)) !8
            allocate (a_gammasun_ebstemp        (numpatch)) !9
            allocate (a_gammasun_dbstemp        (numpatch)) !10
            allocate (a_gammasun_dbsboreal      (numpatch)) !11
            allocate (a_gammasun_c3arcgrass     (numpatch)) !12
            allocate (a_gammasun_c3grass        (numpatch)) !13
            allocate (a_gammasun_c4grass        (numpatch)) !14
            allocate (a_gammasha_enftemp        (numpatch)) !1
            allocate (a_gammasha_enfboreal      (numpatch)) !2
            allocate (a_gammasha_dnfboreal      (numpatch)) !3
            allocate (a_gammasha_ebftrop        (numpatch)) !4
            allocate (a_gammasha_ebftemp        (numpatch)) !5
            allocate (a_gammasha_dbftrop        (numpatch)) !6
            allocate (a_gammasha_dbftemp        (numpatch)) !7
            allocate (a_gammasha_dbfboreal      (numpatch)) !8
            allocate (a_gammasha_ebstemp        (numpatch)) !9
            allocate (a_gammasha_dbstemp        (numpatch)) !10
            allocate (a_gammasha_dbsboreal      (numpatch)) !11
            allocate (a_gammasha_c3arcgrass     (numpatch)) !12
            allocate (a_gammasha_c3grass        (numpatch)) !13
            allocate (a_gammasha_c4grass        (numpatch)) !14
            allocate (a_RuBPlimitfrac_sun_enftemp          (numpatch)) !1
            allocate (a_RuBPlimitfrac_sun_enfboreal        (numpatch)) !2
            allocate (a_RuBPlimitfrac_sun_dnfboreal        (numpatch)) !3
            allocate (a_RuBPlimitfrac_sun_ebftrop          (numpatch)) !4
            allocate (a_RuBPlimitfrac_sun_ebftemp          (numpatch)) !5
            allocate (a_RuBPlimitfrac_sun_dbftrop          (numpatch)) !6
            allocate (a_RuBPlimitfrac_sun_dbftemp          (numpatch)) !7
            allocate (a_RuBPlimitfrac_sun_dbfboreal        (numpatch)) !8
            allocate (a_RuBPlimitfrac_sun_ebstemp          (numpatch)) !9
            allocate (a_RuBPlimitfrac_sun_dbstemp          (numpatch)) !10
            allocate (a_RuBPlimitfrac_sun_dbsboreal        (numpatch)) !11
            allocate (a_RuBPlimitfrac_sun_c3arcgrass       (numpatch)) !12
            allocate (a_RuBPlimitfrac_sun_c3grass          (numpatch)) !13
            allocate (a_RuBPlimitfrac_sun_c4grass          (numpatch)) !14
            allocate (a_RuBPlimitfrac_sha_enftemp          (numpatch)) !1
            allocate (a_RuBPlimitfrac_sha_enfboreal        (numpatch)) !2
            allocate (a_RuBPlimitfrac_sha_dnfboreal        (numpatch)) !3
            allocate (a_RuBPlimitfrac_sha_ebftrop          (numpatch)) !4
            allocate (a_RuBPlimitfrac_sha_ebftemp          (numpatch)) !5
            allocate (a_RuBPlimitfrac_sha_dbftrop          (numpatch)) !6
            allocate (a_RuBPlimitfrac_sha_dbftemp          (numpatch)) !7
            allocate (a_RuBPlimitfrac_sha_dbfboreal        (numpatch)) !8
            allocate (a_RuBPlimitfrac_sha_ebstemp          (numpatch)) !9
            allocate (a_RuBPlimitfrac_sha_dbstemp          (numpatch)) !10
            allocate (a_RuBPlimitfrac_sha_dbsboreal        (numpatch)) !11
            allocate (a_RuBPlimitfrac_sha_c3arcgrass       (numpatch)) !12
            allocate (a_RuBPlimitfrac_sha_c3grass          (numpatch)) !13
            allocate (a_RuBPlimitfrac_sha_c4grass          (numpatch)) !14
            allocate (a_Rubiscolimitfrac_sun_enftemp          (numpatch)) !1
            allocate (a_Rubiscolimitfrac_sun_enfboreal        (numpatch)) !2
            allocate (a_Rubiscolimitfrac_sun_dnfboreal        (numpatch)) !3
            allocate (a_Rubiscolimitfrac_sun_ebftrop          (numpatch)) !4
            allocate (a_Rubiscolimitfrac_sun_ebftemp          (numpatch)) !5
            allocate (a_Rubiscolimitfrac_sun_dbftrop          (numpatch)) !6
            allocate (a_Rubiscolimitfrac_sun_dbftemp          (numpatch)) !7
            allocate (a_Rubiscolimitfrac_sun_dbfboreal        (numpatch)) !8
            allocate (a_Rubiscolimitfrac_sun_ebstemp          (numpatch)) !9
            allocate (a_Rubiscolimitfrac_sun_dbstemp          (numpatch)) !10
            allocate (a_Rubiscolimitfrac_sun_dbsboreal        (numpatch)) !11
            allocate (a_Rubiscolimitfrac_sun_c3arcgrass       (numpatch)) !12
            allocate (a_Rubiscolimitfrac_sun_c3grass          (numpatch)) !13
            allocate (a_Rubiscolimitfrac_sun_c4grass          (numpatch)) !14
            allocate (a_Rubiscolimitfrac_sha_enftemp          (numpatch)) !1
            allocate (a_Rubiscolimitfrac_sha_enfboreal        (numpatch)) !2
            allocate (a_Rubiscolimitfrac_sha_dnfboreal        (numpatch)) !3
            allocate (a_Rubiscolimitfrac_sha_ebftrop          (numpatch)) !4
            allocate (a_Rubiscolimitfrac_sha_ebftemp          (numpatch)) !5
            allocate (a_Rubiscolimitfrac_sha_dbftrop          (numpatch)) !6
            allocate (a_Rubiscolimitfrac_sha_dbftemp          (numpatch)) !7
            allocate (a_Rubiscolimitfrac_sha_dbfboreal        (numpatch)) !8
            allocate (a_Rubiscolimitfrac_sha_ebstemp          (numpatch)) !9
            allocate (a_Rubiscolimitfrac_sha_dbstemp          (numpatch)) !10
            allocate (a_Rubiscolimitfrac_sha_dbsboreal        (numpatch)) !11
            allocate (a_Rubiscolimitfrac_sha_c3arcgrass       (numpatch)) !12
            allocate (a_Rubiscolimitfrac_sha_c3grass          (numpatch)) !13
            allocate (a_Rubiscolimitfrac_sha_c4grass          (numpatch)) !14
            allocate (a_Sinklimitfrac_sun_enftemp          (numpatch)) !1
            allocate (a_Sinklimitfrac_sun_enfboreal        (numpatch)) !2
            allocate (a_Sinklimitfrac_sun_dnfboreal        (numpatch)) !3
            allocate (a_Sinklimitfrac_sun_ebftrop          (numpatch)) !4
            allocate (a_Sinklimitfrac_sun_ebftemp          (numpatch)) !5
            allocate (a_Sinklimitfrac_sun_dbftrop          (numpatch)) !6
            allocate (a_Sinklimitfrac_sun_dbftemp          (numpatch)) !7
            allocate (a_Sinklimitfrac_sun_dbfboreal        (numpatch)) !8
            allocate (a_Sinklimitfrac_sun_ebstemp          (numpatch)) !9
            allocate (a_Sinklimitfrac_sun_dbstemp          (numpatch)) !10
            allocate (a_Sinklimitfrac_sun_dbsboreal        (numpatch)) !11
            allocate (a_Sinklimitfrac_sun_c3arcgrass       (numpatch)) !12
            allocate (a_Sinklimitfrac_sun_c3grass          (numpatch)) !13
            allocate (a_Sinklimitfrac_sun_c4grass          (numpatch)) !14
            allocate (a_Sinklimitfrac_sha_enftemp          (numpatch)) !1
            allocate (a_Sinklimitfrac_sha_enfboreal        (numpatch)) !2
            allocate (a_Sinklimitfrac_sha_dnfboreal        (numpatch)) !3
            allocate (a_Sinklimitfrac_sha_ebftrop          (numpatch)) !4
            allocate (a_Sinklimitfrac_sha_ebftemp          (numpatch)) !5
            allocate (a_Sinklimitfrac_sha_dbftrop          (numpatch)) !6
            allocate (a_Sinklimitfrac_sha_dbftemp          (numpatch)) !7
            allocate (a_Sinklimitfrac_sha_dbfboreal        (numpatch)) !8
            allocate (a_Sinklimitfrac_sha_ebstemp          (numpatch)) !9
            allocate (a_Sinklimitfrac_sha_dbstemp          (numpatch)) !10
            allocate (a_Sinklimitfrac_sha_dbsboreal        (numpatch)) !11
            allocate (a_Sinklimitfrac_sha_c3arcgrass       (numpatch)) !12
            allocate (a_Sinklimitfrac_sha_c3grass          (numpatch)) !13
            allocate (a_Sinklimitfrac_sha_c4grass          (numpatch)) !14
            allocate (a_rstfacsun_enftemp        (numpatch)) !1
            allocate (a_rstfacsun_enfboreal      (numpatch)) !2
            allocate (a_rstfacsun_dnfboreal      (numpatch)) !3
            allocate (a_rstfacsun_ebftrop        (numpatch)) !4
            allocate (a_rstfacsun_ebftemp        (numpatch)) !5
            allocate (a_rstfacsun_dbftrop        (numpatch)) !6
            allocate (a_rstfacsun_dbftemp        (numpatch)) !7
            allocate (a_rstfacsun_dbfboreal      (numpatch)) !8
            allocate (a_rstfacsun_ebstemp        (numpatch)) !9
            allocate (a_rstfacsun_dbstemp        (numpatch)) !10
            allocate (a_rstfacsun_dbsboreal      (numpatch)) !11
            allocate (a_rstfacsun_c3arcgrass     (numpatch)) !12
            allocate (a_rstfacsun_c3grass        (numpatch)) !13
            allocate (a_rstfacsun_c4grass        (numpatch)) !14
            allocate (a_rstfacsha_enftemp        (numpatch)) !1
            allocate (a_rstfacsha_enfboreal      (numpatch)) !2
            allocate (a_rstfacsha_dnfboreal      (numpatch)) !3
            allocate (a_rstfacsha_ebftrop        (numpatch)) !4
            allocate (a_rstfacsha_ebftemp        (numpatch)) !5
            allocate (a_rstfacsha_dbftrop        (numpatch)) !6
            allocate (a_rstfacsha_dbftemp        (numpatch)) !7
            allocate (a_rstfacsha_dbfboreal      (numpatch)) !8
            allocate (a_rstfacsha_ebstemp        (numpatch)) !9
            allocate (a_rstfacsha_dbstemp        (numpatch)) !10
            allocate (a_rstfacsha_dbsboreal      (numpatch)) !11
            allocate (a_rstfacsha_c3arcgrass     (numpatch)) !12
            allocate (a_rstfacsha_c3grass        (numpatch)) !13
            allocate (a_rstfacsha_c4grass        (numpatch)) !14
            allocate (a_lambdasun_enftemp        (numpatch)) !1
            allocate (a_lambdasun_enfboreal      (numpatch)) !2
            allocate (a_lambdasun_dnfboreal      (numpatch)) !3
            allocate (a_lambdasun_ebftrop        (numpatch)) !4
            allocate (a_lambdasun_ebftemp        (numpatch)) !5
            allocate (a_lambdasun_dbftrop        (numpatch)) !6
            allocate (a_lambdasun_dbftemp        (numpatch)) !7
            allocate (a_lambdasun_dbfboreal      (numpatch)) !8
            allocate (a_lambdasun_ebstemp        (numpatch)) !9
            allocate (a_lambdasun_dbstemp        (numpatch)) !10
            allocate (a_lambdasun_dbsboreal      (numpatch)) !11
            allocate (a_lambdasun_c3arcgrass     (numpatch)) !12
            allocate (a_lambdasun_c3grass        (numpatch)) !13
            allocate (a_lambdasun_c4grass        (numpatch)) !14
            allocate (a_lambdasha_enftemp        (numpatch)) !1
            allocate (a_lambdasha_enfboreal      (numpatch)) !2
            allocate (a_lambdasha_dnfboreal      (numpatch)) !3
            allocate (a_lambdasha_ebftrop        (numpatch)) !4
            allocate (a_lambdasha_ebftemp        (numpatch)) !5
            allocate (a_lambdasha_dbftrop        (numpatch)) !6
            allocate (a_lambdasha_dbftemp        (numpatch)) !7
            allocate (a_lambdasha_dbfboreal      (numpatch)) !8
            allocate (a_lambdasha_ebstemp        (numpatch)) !9
            allocate (a_lambdasha_dbstemp        (numpatch)) !10
            allocate (a_lambdasha_dbsboreal      (numpatch)) !11
            allocate (a_lambdasha_c3arcgrass     (numpatch)) !12
            allocate (a_lambdasha_c3grass        (numpatch)) !13
            allocate (a_lambdasha_c4grass        (numpatch)) !14
            allocate (a_lambda                   (numpatch)) !1
#endif 
#endif

            allocate (a_qcharge   (numpatch))

            allocate (a_t_grnd    (numpatch)) 
            allocate (a_tleaf     (numpatch)) 
            allocate (a_ldew_rain (numpatch)) 
            allocate (a_ldew_snow (numpatch)) 
            allocate (a_ldew      (numpatch)) 
            allocate (a_scv       (numpatch)) 
            allocate (a_snowdp    (numpatch)) 
            allocate (a_fsno      (numpatch)) 
            allocate (a_sigf      (numpatch)) 
            allocate (a_green     (numpatch)) 
            allocate (a_lai       (numpatch)) 
            allocate (a_laisun    (numpatch)) 
            allocate (a_laisha    (numpatch)) 
            allocate (a_sai       (numpatch)) 

            allocate (a_alb   (2,2,numpatch))    

            allocate (a_emis      (numpatch))
            allocate (a_z0m       (numpatch))
            allocate (a_trad      (numpatch))
            allocate (a_tref      (numpatch))
            allocate (a_qref      (numpatch))
            allocate (a_rain      (numpatch))
            allocate (a_snow      (numpatch))  

#ifdef BGC
            allocate (a_leafc              (numpatch))
            allocate (a_leafc_storage      (numpatch))
            allocate (a_leafc_xfer         (numpatch))
            allocate (a_frootc             (numpatch))
            allocate (a_frootc_storage     (numpatch))
            allocate (a_frootc_xfer        (numpatch))
            allocate (a_livestemc          (numpatch))
            allocate (a_livestemc_storage  (numpatch))
            allocate (a_livestemc_xfer     (numpatch))
            allocate (a_deadstemc          (numpatch))
            allocate (a_deadstemc_storage  (numpatch))
            allocate (a_deadstemc_xfer     (numpatch))
            allocate (a_livecrootc         (numpatch))
            allocate (a_livecrootc_storage (numpatch))
            allocate (a_livecrootc_xfer    (numpatch))
            allocate (a_deadcrootc         (numpatch))
            allocate (a_deadcrootc_storage (numpatch))
            allocate (a_deadcrootc_xfer    (numpatch))
            allocate (a_grainc             (numpatch))
            allocate (a_grainc_storage     (numpatch))
            allocate (a_grainc_xfer        (numpatch))
            allocate (a_leafn              (numpatch))
            allocate (a_leafn_storage      (numpatch))
            allocate (a_leafn_xfer         (numpatch))
            allocate (a_frootn             (numpatch))
            allocate (a_frootn_storage     (numpatch))
            allocate (a_frootn_xfer        (numpatch))
            allocate (a_livestemn          (numpatch))
            allocate (a_livestemn_storage  (numpatch))
            allocate (a_livestemn_xfer     (numpatch))
            allocate (a_deadstemn          (numpatch))
            allocate (a_deadstemn_storage  (numpatch))
            allocate (a_deadstemn_xfer     (numpatch))
            allocate (a_livecrootn         (numpatch))
            allocate (a_livecrootn_storage (numpatch))
            allocate (a_livecrootn_xfer    (numpatch))
            allocate (a_deadcrootn         (numpatch))
            allocate (a_deadcrootn_storage (numpatch))
            allocate (a_deadcrootn_xfer    (numpatch))
            allocate (a_grainn             (numpatch))
            allocate (a_grainn_storage     (numpatch))
            allocate (a_grainn_xfer        (numpatch))
            allocate (a_retransn           (numpatch))
            allocate (a_gpp                (numpatch))
            allocate (a_downreg            (numpatch))
            allocate (a_ar                 (numpatch))
            allocate (a_cwdprod            (numpatch))
            allocate (a_cwddecomp          (numpatch))
            allocate (a_hr                 (numpatch))
            allocate (a_fpg                (numpatch))
            allocate (a_fpi                (numpatch))
            allocate (a_gpp_enftemp        (numpatch)) !1
            allocate (a_gpp_enfboreal      (numpatch)) !2
            allocate (a_gpp_dnfboreal      (numpatch)) !3
            allocate (a_gpp_ebftrop        (numpatch)) !4
            allocate (a_gpp_ebftemp        (numpatch)) !5
            allocate (a_gpp_dbftrop        (numpatch)) !6
            allocate (a_gpp_dbftemp        (numpatch)) !7
            allocate (a_gpp_dbfboreal      (numpatch)) !8
            allocate (a_gpp_ebstemp        (numpatch)) !9
            allocate (a_gpp_dbstemp        (numpatch)) !10
            allocate (a_gpp_dbsboreal      (numpatch)) !11
            allocate (a_gpp_c3arcgrass     (numpatch)) !12
            allocate (a_gpp_c3grass        (numpatch)) !13
            allocate (a_gpp_c4grass        (numpatch)) !14
            allocate (a_leafc_enftemp      (numpatch)) !1
            allocate (a_leafc_enfboreal    (numpatch)) !2
            allocate (a_leafc_dnfboreal    (numpatch)) !3
            allocate (a_leafc_ebftrop      (numpatch)) !4
            allocate (a_leafc_ebftemp      (numpatch)) !5
            allocate (a_leafc_dbftrop      (numpatch)) !6
            allocate (a_leafc_dbftemp      (numpatch)) !7
            allocate (a_leafc_dbfboreal    (numpatch)) !8
            allocate (a_leafc_ebstemp      (numpatch)) !9
            allocate (a_leafc_dbstemp      (numpatch)) !10
            allocate (a_leafc_dbsboreal    (numpatch)) !11
            allocate (a_leafc_c3arcgrass   (numpatch)) !12
            allocate (a_leafc_c3grass      (numpatch)) !13
            allocate (a_leafc_c4grass      (numpatch)) !14
#ifdef NITRIF
            allocate (a_O2_DECOMP_DEPTH_UNSAT (1:nl_soil,numpatch))
            allocate (a_CONC_O2_UNSAT         (1:nl_soil,numpatch))
#endif
#ifdef CROP
            allocate (a_pdcorn             (numpatch))
            allocate (a_pdswheat           (numpatch))
            allocate (a_pdwwheat           (numpatch))
            allocate (a_pdsoybean          (numpatch))
            allocate (a_pdcotton           (numpatch))
            allocate (a_pdrice1            (numpatch))
            allocate (a_pdrice2            (numpatch))
            allocate (a_pdsugarcane        (numpatch))
            allocate (a_plantdate          (numpatch))
            allocate (a_fertnitro_corn     (numpatch))
            allocate (a_fertnitro_swheat   (numpatch))
            allocate (a_fertnitro_wwheat   (numpatch))
            allocate (a_fertnitro_soybean  (numpatch))
            allocate (a_fertnitro_cotton   (numpatch))
            allocate (a_fertnitro_rice1    (numpatch))
            allocate (a_fertnitro_rice2    (numpatch))
            allocate (a_fertnitro_sugarcane(numpatch))
            allocate (a_cphase             (numpatch))
            allocate (a_hui                (numpatch))
            allocate (a_gddmaturity        (numpatch))
            allocate (a_gddplant           (numpatch))
            allocate (a_vf                 (numpatch))
            allocate (a_cropprod1c         (numpatch))
            allocate (a_cropprod1c_loss    (numpatch))
            allocate (a_cropseedc_deficit  (numpatch))
            allocate (a_grainc_to_cropprodc(numpatch))
            allocate (a_grainc_to_seed     (numpatch))
            allocate (a_fert_to_sminn      (numpatch))
#endif
            allocate (a_ndep_to_sminn      (numpatch))
#ifdef Fire
            allocate (a_abm                (numpatch))
            allocate (a_gdp                (numpatch))
            allocate (a_peatf              (numpatch))
            allocate (a_hdm                (numpatch))
            allocate (a_lnfm               (numpatch))
#endif
#endif
#ifdef OzoneStress
            allocate (a_ozone              (numpatch))
#endif
            allocate (a_t_soisno    (maxsnl+1:nl_soil,numpatch))    
            allocate (a_wliq_soisno (maxsnl+1:nl_soil,numpatch))
            allocate (a_wice_soisno (maxsnl+1:nl_soil,numpatch))
            allocate (a_h2osoi      (1:nl_soil,       numpatch))
            allocate (a_rootr       (1:nl_soil,       numpatch))
            allocate (a_BD_all      (1:nl_soil,       numpatch))
            allocate (a_wfc         (1:nl_soil,       numpatch))
            allocate (a_OM_density  (1:nl_soil,       numpatch))
#ifdef PLANT_HYDRAULIC_STRESS
            allocate (a_vegwp       (1:nvegwcs,       numpatch))
#endif
            allocate (a_t_lake      (nl_lake,numpatch)) 
            allocate (a_lake_icefrac(nl_lake,numpatch)) 

#ifdef BGC
            allocate (a_litr1c_vr   (1:nl_soil,       numpatch))
            allocate (a_litr2c_vr   (1:nl_soil,       numpatch))
            allocate (a_litr3c_vr   (1:nl_soil,       numpatch))
            allocate (a_soil1c_vr   (1:nl_soil,       numpatch))
            allocate (a_soil2c_vr   (1:nl_soil,       numpatch))
            allocate (a_soil3c_vr   (1:nl_soil,       numpatch))
            allocate (a_cwdc_vr     (1:nl_soil,       numpatch))
            allocate (a_litr1n_vr   (1:nl_soil,       numpatch))
            allocate (a_litr2n_vr   (1:nl_soil,       numpatch))
            allocate (a_litr3n_vr   (1:nl_soil,       numpatch))
            allocate (a_soil1n_vr   (1:nl_soil,       numpatch))
            allocate (a_soil2n_vr   (1:nl_soil,       numpatch))
            allocate (a_soil3n_vr   (1:nl_soil,       numpatch))
            allocate (a_cwdn_vr     (1:nl_soil,       numpatch))
            allocate (a_sminn_vr    (1:nl_soil,       numpatch))
            allocate (decomp_vr_tmp (1:nl_soil,       numpatch))
#endif

            allocate (a_ustar     (numpatch)) 
            allocate (a_tstar     (numpatch))
            allocate (a_qstar     (numpatch))
            allocate (a_zol       (numpatch))
            allocate (a_rib       (numpatch))
            allocate (a_fm        (numpatch))
            allocate (a_fh        (numpatch))
            allocate (a_fq        (numpatch))

            allocate (a_us10m     (numpatch)) 
            allocate (a_vs10m     (numpatch)) 
            allocate (a_fm10m     (numpatch)) 

            allocate (a_sr        (numpatch))
            allocate (a_solvd     (numpatch))
            allocate (a_solvi     (numpatch))
            allocate (a_solnd     (numpatch))
            allocate (a_solni     (numpatch))
            allocate (a_srvd      (numpatch))
            allocate (a_srvi      (numpatch))
            allocate (a_srnd      (numpatch))
            allocate (a_srni      (numpatch))
            allocate (a_solvdln   (numpatch))
            allocate (a_solviln   (numpatch))
            allocate (a_solndln   (numpatch))
            allocate (a_solniln   (numpatch))
            allocate (a_srvdln    (numpatch))
            allocate (a_srviln    (numpatch))
            allocate (a_srndln    (numpatch))
            allocate (a_srniln    (numpatch))

            allocate (nac_ln      (numpatch))

         end if
      end if

   end subroutine allocate_acc_fluxes

   subroutine deallocate_acc_fluxes ()

      use spmd_task
      use mod_landpatch, only : numpatch
      implicit none

      if (p_is_worker) then
         if (numpatch > 0) then

            deallocate (a_us     )  
            deallocate (a_vs     )
            deallocate (a_t      )
            deallocate (a_q      )
            deallocate (a_prc    )
            deallocate (a_prl    )
            deallocate (a_pbot   )
            deallocate (a_frl    )
            deallocate (a_solarin)

            deallocate (a_taux      )
            deallocate (a_tauy      )
            deallocate (a_fsena     )
            deallocate (a_lfevpa    )
            deallocate (a_fevpa     )
            deallocate (a_fsenl     )
            deallocate (a_fevpl     )
            deallocate (a_etr       )
            deallocate (a_fseng     )
            deallocate (a_fevpg     )
            deallocate (a_fgrnd     )
            deallocate (a_sabvsun   )  
            deallocate (a_sabvsha   ) 
            deallocate (a_sabg      )
            deallocate (a_olrg      )
            deallocate (a_rnet      )
            deallocate (a_xerr      )
            deallocate (a_zerr      )
            deallocate (a_rsur      )
            deallocate (a_rnof      )
            deallocate (a_qintr     )
            deallocate (a_qinfl     )
            deallocate (a_qdrip     )
            deallocate (a_rstfacsun )
            deallocate (a_rstfacsha )
            deallocate (a_gs_sun )
            deallocate (a_gs_sha )
            deallocate (a_dpond     )

            deallocate (a_zwt       )
            deallocate (a_wa        )
            deallocate (a_wat       )
            deallocate (a_assim     )
            deallocate (a_respc     )

            deallocate (a_qcharge   )

            deallocate (a_t_grnd    ) 
            deallocate (a_tleaf     )
            deallocate (a_ldew_rain ) 
            deallocate (a_ldew_snow ) 
            deallocate (a_ldew      ) 
            deallocate (a_scv       ) 
            deallocate (a_snowdp    ) 
            deallocate (a_fsno      ) 
            deallocate (a_sigf      ) 
            deallocate (a_green     ) 
            deallocate (a_lai       ) 
            deallocate (a_laisun    ) 
            deallocate (a_laisha    ) 
            deallocate (a_sai       ) 

            deallocate (a_alb  ) 

            deallocate (a_emis      )
            deallocate (a_z0m       )
            deallocate (a_trad      )
            deallocate (a_tref      )
            deallocate (a_qref      )
            deallocate (a_rain      )
            deallocate (a_snow      )  
#ifdef BGC
            deallocate (a_leafc              )
            deallocate (a_leafc_storage      )
            deallocate (a_leafc_xfer         )
            deallocate (a_frootc             )
            deallocate (a_frootc_storage     )
            deallocate (a_frootc_xfer        )
            deallocate (a_livestemc          )
            deallocate (a_livestemc_storage  )
            deallocate (a_livestemc_xfer     )
            deallocate (a_deadstemc          )
            deallocate (a_deadstemc_storage  )
            deallocate (a_deadstemc_xfer     )
            deallocate (a_livecrootc         )
            deallocate (a_livecrootc_storage )
            deallocate (a_livecrootc_xfer    )
            deallocate (a_deadcrootc         )
            deallocate (a_deadcrootc_storage )
            deallocate (a_deadcrootc_xfer    )
            deallocate (a_grainc             )
            deallocate (a_grainc_storage     )
            deallocate (a_grainc_xfer        )
            deallocate (a_leafn              )
            deallocate (a_leafn_storage      )
            deallocate (a_leafn_xfer         )
            deallocate (a_frootn             )
            deallocate (a_frootn_storage     )
            deallocate (a_frootn_xfer        )
            deallocate (a_livestemn          )
            deallocate (a_livestemn_storage  )
            deallocate (a_livestemn_xfer     )
            deallocate (a_deadstemn          )
            deallocate (a_deadstemn_storage  )
            deallocate (a_deadstemn_xfer     )
            deallocate (a_livecrootn         )
            deallocate (a_livecrootn_storage )
            deallocate (a_livecrootn_xfer    )
            deallocate (a_deadcrootn         )
            deallocate (a_deadcrootn_storage )
            deallocate (a_deadcrootn_xfer    )
            deallocate (a_grainn             )
            deallocate (a_grainn_storage     )
            deallocate (a_grainn_xfer        )
            deallocate (a_retransn           )
            deallocate (a_gpp                )
            deallocate (a_downreg            )
            deallocate (a_ar                 )
            deallocate (a_cwdprod            )
            deallocate (a_cwddecomp          )
            deallocate (a_hr                 )
            deallocate (a_fpg                )
            deallocate (a_fpi                )
            deallocate (a_gpp_enftemp        ) !1
            deallocate (a_gpp_enfboreal      ) !2
            deallocate (a_gpp_dnfboreal      ) !3
            deallocate (a_gpp_ebftrop        ) !4
            deallocate (a_gpp_ebftemp        ) !5
            deallocate (a_gpp_dbftrop        ) !6
            deallocate (a_gpp_dbftemp        ) !7
            deallocate (a_gpp_dbfboreal      ) !8
            deallocate (a_gpp_ebstemp        ) !9
            deallocate (a_gpp_dbstemp        ) !10
            deallocate (a_gpp_dbsboreal      ) !11
            deallocate (a_gpp_c3arcgrass     ) !12
            deallocate (a_gpp_c3grass        ) !13
            deallocate (a_gpp_c4grass        ) !14
            deallocate (a_leafc_enftemp      ) !1
            deallocate (a_leafc_enfboreal    ) !2
            deallocate (a_leafc_dnfboreal    ) !3
            deallocate (a_leafc_ebftrop      ) !4
            deallocate (a_leafc_ebftemp      ) !5
            deallocate (a_leafc_dbftrop      ) !6
            deallocate (a_leafc_dbftemp      ) !7
            deallocate (a_leafc_dbfboreal    ) !8
            deallocate (a_leafc_ebstemp      ) !9
            deallocate (a_leafc_dbstemp      ) !10
            deallocate (a_leafc_dbsboreal    ) !11
            deallocate (a_leafc_c3arcgrass   ) !12
            deallocate (a_leafc_c3grass      ) !13
            deallocate (a_leafc_c4grass      ) !14
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
            deallocate (a_assim_RuBP_sun_enftemp        ) !1
            deallocate (a_assim_RuBP_sun_enfboreal      ) !2
            deallocate (a_assim_RuBP_sun_dnfboreal      ) !3
            deallocate (a_assim_RuBP_sun_ebftrop        ) !4
            deallocate (a_assim_RuBP_sun_ebftemp        ) !5
            deallocate (a_assim_RuBP_sun_dbftrop        ) !6
            deallocate (a_assim_RuBP_sun_dbftemp        ) !7
            deallocate (a_assim_RuBP_sun_dbfboreal      ) !8
            deallocate (a_assim_RuBP_sun_ebstemp        ) !9
            deallocate (a_assim_RuBP_sun_dbstemp        ) !10
            deallocate (a_assim_RuBP_sun_dbsboreal      ) !11
            deallocate (a_assim_RuBP_sun_c3arcgrass     ) !12
            deallocate (a_assim_RuBP_sun_c3grass        ) !13
            deallocate (a_assim_RuBP_sun_c4grass        ) !14
            deallocate (a_assim_RuBP_sha_enftemp        ) !1
            deallocate (a_assim_RuBP_sha_enfboreal      ) !2
            deallocate (a_assim_RuBP_sha_dnfboreal      ) !3
            deallocate (a_assim_RuBP_sha_ebftrop        ) !4
            deallocate (a_assim_RuBP_sha_ebftemp        ) !5
            deallocate (a_assim_RuBP_sha_dbftrop        ) !6
            deallocate (a_assim_RuBP_sha_dbftemp        ) !7
            deallocate (a_assim_RuBP_sha_dbfboreal      ) !8
            deallocate (a_assim_RuBP_sha_ebstemp        ) !9
            deallocate (a_assim_RuBP_sha_dbstemp        ) !10
            deallocate (a_assim_RuBP_sha_dbsboreal      ) !11
            deallocate (a_assim_RuBP_sha_c3arcgrass     ) !12
            deallocate (a_assim_RuBP_sha_c3grass        ) !13
            deallocate (a_assim_RuBP_sha_c4grass        ) !14
            deallocate (a_assim_Rubisco_sun_enftemp        ) !1
            deallocate (a_assim_Rubisco_sun_enfboreal      ) !2
            deallocate (a_assim_Rubisco_sun_dnfboreal      ) !3
            deallocate (a_assim_Rubisco_sun_ebftrop        ) !4
            deallocate (a_assim_Rubisco_sun_ebftemp        ) !5
            deallocate (a_assim_Rubisco_sun_dbftrop        ) !6
            deallocate (a_assim_Rubisco_sun_dbftemp        ) !7
            deallocate (a_assim_Rubisco_sun_dbfboreal      ) !8
            deallocate (a_assim_Rubisco_sun_ebstemp        ) !9
            deallocate (a_assim_Rubisco_sun_dbstemp        ) !10
            deallocate (a_assim_Rubisco_sun_dbsboreal      ) !11
            deallocate (a_assim_Rubisco_sun_c3arcgrass     ) !12
            deallocate (a_assim_Rubisco_sun_c3grass        ) !13
            deallocate (a_assim_Rubisco_sun_c4grass        ) !14
            deallocate (a_assim_Rubisco_sha_enftemp        ) !1
            deallocate (a_assim_Rubisco_sha_enfboreal      ) !2
            deallocate (a_assim_Rubisco_sha_dnfboreal      ) !3
            deallocate (a_assim_Rubisco_sha_ebftrop        ) !4
            deallocate (a_assim_Rubisco_sha_ebftemp        ) !5
            deallocate (a_assim_Rubisco_sha_dbftrop        ) !6
            deallocate (a_assim_Rubisco_sha_dbftemp        ) !7
            deallocate (a_assim_Rubisco_sha_dbfboreal      ) !8
            deallocate (a_assim_Rubisco_sha_ebstemp        ) !9
            deallocate (a_assim_Rubisco_sha_dbstemp        ) !10
            deallocate (a_assim_Rubisco_sha_dbsboreal      ) !11
            deallocate (a_assim_Rubisco_sha_c3arcgrass     ) !12
            deallocate (a_assim_Rubisco_sha_c3grass        ) !13
            deallocate (a_assim_Rubisco_sha_c4grass        ) !14
            deallocate (a_assimsun_enftemp        ) !1
            deallocate (a_assimsun_enfboreal      ) !2
            deallocate (a_assimsun_dnfboreal      ) !3
            deallocate (a_assimsun_ebftrop        ) !4
            deallocate (a_assimsun_ebftemp        ) !5
            deallocate (a_assimsun_dbftrop        ) !6
            deallocate (a_assimsun_dbftemp        ) !7
            deallocate (a_assimsun_dbfboreal      ) !8
            deallocate (a_assimsun_ebstemp        ) !9
            deallocate (a_assimsun_dbstemp        ) !10
            deallocate (a_assimsun_dbsboreal      ) !11
            deallocate (a_assimsun_c3arcgrass     ) !12
            deallocate (a_assimsun_c3grass        ) !13
            deallocate (a_assimsun_c4grass        ) !14
            deallocate (a_assimsha_enftemp        ) !1
            deallocate (a_assimsha_enfboreal      ) !2
            deallocate (a_assimsha_dnfboreal      ) !3
            deallocate (a_assimsha_ebftrop        ) !4
            deallocate (a_assimsha_ebftemp        ) !5
            deallocate (a_assimsha_dbftrop        ) !6
            deallocate (a_assimsha_dbftemp        ) !7
            deallocate (a_assimsha_dbfboreal      ) !8
            deallocate (a_assimsha_ebstemp        ) !9
            deallocate (a_assimsha_dbstemp        ) !10
            deallocate (a_assimsha_dbsboreal      ) !11
            deallocate (a_assimsha_c3arcgrass     ) !12
            deallocate (a_assimsha_c3grass        ) !13
            deallocate (a_assimsha_c4grass        ) !14
            deallocate (a_etrsun_enftemp        ) !1
            deallocate (a_etrsun_enfboreal      ) !2
            deallocate (a_etrsun_dnfboreal      ) !3
            deallocate (a_etrsun_ebftrop        ) !4
            deallocate (a_etrsun_ebftemp        ) !5
            deallocate (a_etrsun_dbftrop        ) !6
            deallocate (a_etrsun_dbftemp        ) !7
            deallocate (a_etrsun_dbfboreal      ) !8
            deallocate (a_etrsun_ebstemp        ) !9
            deallocate (a_etrsun_dbstemp        ) !10
            deallocate (a_etrsun_dbsboreal      ) !11
            deallocate (a_etrsun_c3arcgrass     ) !12
            deallocate (a_etrsun_c3grass        ) !13
            deallocate (a_etrsun_c4grass        ) !14
            deallocate (a_etrsha_enftemp        ) !1
            deallocate (a_etrsha_enfboreal      ) !2
            deallocate (a_etrsha_dnfboreal      ) !3
            deallocate (a_etrsha_ebftrop        ) !4
            deallocate (a_etrsha_ebftemp        ) !5
            deallocate (a_etrsha_dbftrop        ) !6
            deallocate (a_etrsha_dbftemp        ) !7
            deallocate (a_etrsha_dbfboreal      ) !8
            deallocate (a_etrsha_ebstemp        ) !9
            deallocate (a_etrsha_dbstemp        ) !10
            deallocate (a_etrsha_dbsboreal      ) !11
            deallocate (a_etrsha_c3arcgrass     ) !12
            deallocate (a_etrsha_c3grass        ) !13
            deallocate (a_etrsha_c4grass        ) !14
            deallocate (a_cisun_enftemp        ) !1
            deallocate (a_cisun_enfboreal      ) !2
            deallocate (a_cisun_dnfboreal      ) !3
            deallocate (a_cisun_ebftrop        ) !4
            deallocate (a_cisun_ebftemp        ) !5
            deallocate (a_cisun_dbftrop        ) !6
            deallocate (a_cisun_dbftemp        ) !7
            deallocate (a_cisun_dbfboreal      ) !8
            deallocate (a_cisun_ebstemp        ) !9
            deallocate (a_cisun_dbstemp        ) !10
            deallocate (a_cisun_dbsboreal      ) !11
            deallocate (a_cisun_c3arcgrass     ) !12
            deallocate (a_cisun_c3grass        ) !13
            deallocate (a_cisun_c4grass        ) !14
            deallocate (a_cisha_enftemp        ) !1
            deallocate (a_cisha_enfboreal      ) !2
            deallocate (a_cisha_dnfboreal      ) !3
            deallocate (a_cisha_ebftrop        ) !4
            deallocate (a_cisha_ebftemp        ) !5
            deallocate (a_cisha_dbftrop        ) !6
            deallocate (a_cisha_dbftemp        ) !7
            deallocate (a_cisha_dbfboreal      ) !8
            deallocate (a_cisha_ebstemp        ) !9
            deallocate (a_cisha_dbstemp        ) !10
            deallocate (a_cisha_dbsboreal      ) !11
            deallocate (a_cisha_c3arcgrass     ) !12
            deallocate (a_cisha_c3grass        ) !13
            deallocate (a_cisha_c4grass        ) !14
            deallocate (a_essun_enftemp        ) !1
            deallocate (a_essun_enfboreal      ) !2
            deallocate (a_essun_dnfboreal      ) !3
            deallocate (a_essun_ebftrop        ) !4
            deallocate (a_essun_ebftemp        ) !5
            deallocate (a_essun_dbftrop        ) !6
            deallocate (a_essun_dbftemp        ) !7
            deallocate (a_essun_dbfboreal      ) !8
            deallocate (a_essun_ebstemp        ) !9
            deallocate (a_essun_dbstemp        ) !10
            deallocate (a_essun_dbsboreal      ) !11
            deallocate (a_essun_c3arcgrass     ) !12
            deallocate (a_essun_c3grass        ) !13
            deallocate (a_essun_c4grass        ) !14
            deallocate (a_essha_enftemp        ) !1
            deallocate (a_essha_enfboreal      ) !2
            deallocate (a_essha_dnfboreal      ) !3
            deallocate (a_essha_ebftrop        ) !4
            deallocate (a_essha_ebftemp        ) !5
            deallocate (a_essha_dbftrop        ) !6
            deallocate (a_essha_dbftemp        ) !7
            deallocate (a_essha_dbfboreal      ) !8
            deallocate (a_essha_ebstemp        ) !9
            deallocate (a_essha_dbstemp        ) !10
            deallocate (a_essha_dbsboreal      ) !11
            deallocate (a_essha_c3arcgrass     ) !12
            deallocate (a_essha_c3grass        ) !13
            deallocate (a_essha_c4grass        ) !14
            deallocate (a_gssun_enftemp        ) !1
            deallocate (a_gssun_enfboreal      ) !2
            deallocate (a_gssun_dnfboreal      ) !3
            deallocate (a_gssun_ebftrop        ) !4
            deallocate (a_gssun_ebftemp        ) !5
            deallocate (a_gssun_dbftrop        ) !6
            deallocate (a_gssun_dbftemp        ) !7
            deallocate (a_gssun_dbfboreal      ) !8
            deallocate (a_gssun_ebstemp        ) !9
            deallocate (a_gssun_dbstemp        ) !10
            deallocate (a_gssun_dbsboreal      ) !11
            deallocate (a_gssun_c3arcgrass     ) !12
            deallocate (a_gssun_c3grass        ) !13
            deallocate (a_gssun_c4grass        ) !14
            deallocate (a_gssha_enftemp        ) !1
            deallocate (a_gssha_enfboreal      ) !2
            deallocate (a_gssha_dnfboreal      ) !3
            deallocate (a_gssha_ebftrop        ) !4
            deallocate (a_gssha_ebftemp        ) !5
            deallocate (a_gssha_dbftrop        ) !6
            deallocate (a_gssha_dbftemp        ) !7
            deallocate (a_gssha_dbfboreal      ) !8
            deallocate (a_gssha_ebstemp        ) !9
            deallocate (a_gssha_dbstemp        ) !10
            deallocate (a_gssha_dbsboreal      ) !11
            deallocate (a_gssha_c3arcgrass     ) !12
            deallocate (a_gssha_c3grass        ) !13
            deallocate (a_gssha_c4grass        ) !14
            deallocate (a_gammasun_enftemp        ) !1
            deallocate (a_gammasun_enfboreal      ) !2
            deallocate (a_gammasun_dnfboreal      ) !3
            deallocate (a_gammasun_ebftrop        ) !4
            deallocate (a_gammasun_ebftemp        ) !5
            deallocate (a_gammasun_dbftrop        ) !6
            deallocate (a_gammasun_dbftemp        ) !7
            deallocate (a_gammasun_dbfboreal      ) !8
            deallocate (a_gammasun_ebstemp        ) !9
            deallocate (a_gammasun_dbstemp        ) !10
            deallocate (a_gammasun_dbsboreal      ) !11
            deallocate (a_gammasun_c3arcgrass     ) !12
            deallocate (a_gammasun_c3grass        ) !13
            deallocate (a_gammasun_c4grass        ) !14
            deallocate (a_gammasha_enftemp        ) !1
            deallocate (a_gammasha_enfboreal      ) !2
            deallocate (a_gammasha_dnfboreal      ) !3
            deallocate (a_gammasha_ebftrop        ) !4
            deallocate (a_gammasha_ebftemp        ) !5
            deallocate (a_gammasha_dbftrop        ) !6
            deallocate (a_gammasha_dbftemp        ) !7
            deallocate (a_gammasha_dbfboreal      ) !8
            deallocate (a_gammasha_ebstemp        ) !9
            deallocate (a_gammasha_dbstemp        ) !10
            deallocate (a_gammasha_dbsboreal      ) !11
            deallocate (a_gammasha_c3arcgrass     ) !12
            deallocate (a_gammasha_c3grass        ) !13
            deallocate (a_gammasha_c4grass        ) !14
            deallocate (a_RuBPlimitfrac_sun_enftemp          )
            deallocate (a_RuBPlimitfrac_sun_enfboreal        )
            deallocate (a_RuBPlimitfrac_sun_dnfboreal        )
            deallocate (a_RuBPlimitfrac_sun_ebftrop          )
            deallocate (a_RuBPlimitfrac_sun_ebftemp          )
            deallocate (a_RuBPlimitfrac_sun_dbftrop          )
            deallocate (a_RuBPlimitfrac_sun_dbftemp          )
            deallocate (a_RuBPlimitfrac_sun_dbfboreal        )
            deallocate (a_RuBPlimitfrac_sun_ebstemp          )
            deallocate (a_RuBPlimitfrac_sun_dbstemp          )
            deallocate (a_RuBPlimitfrac_sun_dbsboreal        )
            deallocate (a_RuBPlimitfrac_sun_c3arcgrass       )
            deallocate (a_RuBPlimitfrac_sun_c3grass          )
            deallocate (a_RuBPlimitfrac_sun_c4grass          )
            deallocate (a_RuBPlimitfrac_sha_enftemp          )
            deallocate (a_RuBPlimitfrac_sha_enfboreal        )
            deallocate (a_RuBPlimitfrac_sha_dnfboreal        )
            deallocate (a_RuBPlimitfrac_sha_ebftrop          )
            deallocate (a_RuBPlimitfrac_sha_ebftemp          )
            deallocate (a_RuBPlimitfrac_sha_dbftrop          )
            deallocate (a_RuBPlimitfrac_sha_dbftemp          )
            deallocate (a_RuBPlimitfrac_sha_dbfboreal        )
            deallocate (a_RuBPlimitfrac_sha_ebstemp          )
            deallocate (a_RuBPlimitfrac_sha_dbstemp          )
            deallocate (a_RuBPlimitfrac_sha_dbsboreal        )
            deallocate (a_RuBPlimitfrac_sha_c3arcgrass       )
            deallocate (a_RuBPlimitfrac_sha_c3grass          )
            deallocate (a_RuBPlimitfrac_sha_c4grass          )
            deallocate (a_Rubiscolimitfrac_sun_enftemp          )
            deallocate (a_Rubiscolimitfrac_sun_enfboreal        )
            deallocate (a_Rubiscolimitfrac_sun_dnfboreal        )
            deallocate (a_Rubiscolimitfrac_sun_ebftrop          )
            deallocate (a_Rubiscolimitfrac_sun_ebftemp          )
            deallocate (a_Rubiscolimitfrac_sun_dbftrop          )
            deallocate (a_Rubiscolimitfrac_sun_dbftemp          )
            deallocate (a_Rubiscolimitfrac_sun_dbfboreal        )
            deallocate (a_Rubiscolimitfrac_sun_ebstemp          )
            deallocate (a_Rubiscolimitfrac_sun_dbstemp          )
            deallocate (a_Rubiscolimitfrac_sun_dbsboreal        )
            deallocate (a_Rubiscolimitfrac_sun_c3arcgrass       )
            deallocate (a_Rubiscolimitfrac_sun_c3grass          )
            deallocate (a_Rubiscolimitfrac_sun_c4grass          )
            deallocate (a_Rubiscolimitfrac_sha_enftemp          )
            deallocate (a_Rubiscolimitfrac_sha_enfboreal        )
            deallocate (a_Rubiscolimitfrac_sha_dnfboreal        )
            deallocate (a_Rubiscolimitfrac_sha_ebftrop          )
            deallocate (a_Rubiscolimitfrac_sha_ebftemp          )
            deallocate (a_Rubiscolimitfrac_sha_dbftrop          )
            deallocate (a_Rubiscolimitfrac_sha_dbftemp          )
            deallocate (a_Rubiscolimitfrac_sha_dbfboreal        )
            deallocate (a_Rubiscolimitfrac_sha_ebstemp          )
            deallocate (a_Rubiscolimitfrac_sha_dbstemp          )
            deallocate (a_Rubiscolimitfrac_sha_dbsboreal        )
            deallocate (a_Rubiscolimitfrac_sha_c3arcgrass       )
            deallocate (a_Rubiscolimitfrac_sha_c3grass          )
            deallocate (a_Rubiscolimitfrac_sha_c4grass          )
            deallocate (a_Sinklimitfrac_sun_enftemp          )
            deallocate (a_Sinklimitfrac_sun_enfboreal        )
            deallocate (a_Sinklimitfrac_sun_dnfboreal        )
            deallocate (a_Sinklimitfrac_sun_ebftrop          )
            deallocate (a_Sinklimitfrac_sun_ebftemp          )
            deallocate (a_Sinklimitfrac_sun_dbftrop          )
            deallocate (a_Sinklimitfrac_sun_dbftemp          )
            deallocate (a_Sinklimitfrac_sun_dbfboreal        )
            deallocate (a_Sinklimitfrac_sun_ebstemp          )
            deallocate (a_Sinklimitfrac_sun_dbstemp          )
            deallocate (a_Sinklimitfrac_sun_dbsboreal        )
            deallocate (a_Sinklimitfrac_sun_c3arcgrass       )
            deallocate (a_Sinklimitfrac_sun_c3grass          )
            deallocate (a_Sinklimitfrac_sun_c4grass          )
            deallocate (a_Sinklimitfrac_sha_enftemp          )
            deallocate (a_Sinklimitfrac_sha_enfboreal        )
            deallocate (a_Sinklimitfrac_sha_dnfboreal        )
            deallocate (a_Sinklimitfrac_sha_ebftrop          )
            deallocate (a_Sinklimitfrac_sha_ebftemp          )
            deallocate (a_Sinklimitfrac_sha_dbftrop          )
            deallocate (a_Sinklimitfrac_sha_dbftemp          )
            deallocate (a_Sinklimitfrac_sha_dbfboreal        )
            deallocate (a_Sinklimitfrac_sha_ebstemp          )
            deallocate (a_Sinklimitfrac_sha_dbstemp          )
            deallocate (a_Sinklimitfrac_sha_dbsboreal        )
            deallocate (a_Sinklimitfrac_sha_c3arcgrass       )
            deallocate (a_Sinklimitfrac_sha_c3grass          )
            deallocate (a_Sinklimitfrac_sha_c4grass          )
            deallocate (a_rstfacsun_enftemp        ) !1
            deallocate (a_rstfacsun_enfboreal      ) !2
            deallocate (a_rstfacsun_dnfboreal      ) !3
            deallocate (a_rstfacsun_ebftrop        ) !4
            deallocate (a_rstfacsun_ebftemp        ) !5
            deallocate (a_rstfacsun_dbftrop        ) !6
            deallocate (a_rstfacsun_dbftemp        ) !7
            deallocate (a_rstfacsun_dbfboreal      ) !8
            deallocate (a_rstfacsun_ebstemp        ) !9
            deallocate (a_rstfacsun_dbstemp        ) !10
            deallocate (a_rstfacsun_dbsboreal      ) !11
            deallocate (a_rstfacsun_c3arcgrass     ) !12
            deallocate (a_rstfacsun_c3grass        ) !13
            deallocate (a_rstfacsun_c4grass        ) !14
            deallocate (a_rstfacsha_enftemp        ) !1
            deallocate (a_rstfacsha_enfboreal      ) !2
            deallocate (a_rstfacsha_dnfboreal      ) !3
            deallocate (a_rstfacsha_ebftrop        ) !4
            deallocate (a_rstfacsha_ebftemp        ) !5
            deallocate (a_rstfacsha_dbftrop        ) !6
            deallocate (a_rstfacsha_dbftemp        ) !7
            deallocate (a_rstfacsha_dbfboreal      ) !8
            deallocate (a_rstfacsha_ebstemp        ) !9
            deallocate (a_rstfacsha_dbstemp        ) !10
            deallocate (a_rstfacsha_dbsboreal      ) !11
            deallocate (a_rstfacsha_c3arcgrass     ) !12
            deallocate (a_rstfacsha_c3grass        ) !13
            deallocate (a_rstfacsha_c4grass        ) !14
            deallocate (a_lambdasun_enftemp        ) !1
            deallocate (a_lambdasun_enfboreal      ) !2
            deallocate (a_lambdasun_dnfboreal      ) !3
            deallocate (a_lambdasun_ebftrop        ) !4
            deallocate (a_lambdasun_ebftemp        ) !5
            deallocate (a_lambdasun_dbftrop        ) !6
            deallocate (a_lambdasun_dbftemp        ) !7
            deallocate (a_lambdasun_dbfboreal      ) !8
            deallocate (a_lambdasun_ebstemp        ) !9
            deallocate (a_lambdasun_dbstemp        ) !10
            deallocate (a_lambdasun_dbsboreal      ) !11
            deallocate (a_lambdasun_c3arcgrass     ) !12
            deallocate (a_lambdasun_c3grass        ) !13
            deallocate (a_lambdasun_c4grass        ) !14
            deallocate (a_lambdasha_enftemp        ) !1
            deallocate (a_lambdasha_enfboreal      ) !2
            deallocate (a_lambdasha_dnfboreal      ) !3
            deallocate (a_lambdasha_ebftrop        ) !4
            deallocate (a_lambdasha_ebftemp        ) !5
            deallocate (a_lambdasha_dbftrop        ) !6
            deallocate (a_lambdasha_dbftemp        ) !7
            deallocate (a_lambdasha_dbfboreal      ) !8
            deallocate (a_lambdasha_ebstemp        ) !9
            deallocate (a_lambdasha_dbstemp        ) !10
            deallocate (a_lambdasha_dbsboreal      ) !11
            deallocate (a_lambdasha_c3arcgrass     ) !12
            deallocate (a_lambdasha_c3grass        ) !13
            deallocate (a_lambdasha_c4grass        ) !14
            deallocate (a_lambda                   ) !1
#endif
#endif
#ifdef NITRIF
            deallocate (a_O2_DECOMP_DEPTH_UNSAT )
            deallocate (a_CONC_O2_UNSAT         )
#endif
#ifdef CROP
            deallocate (a_pdcorn             )
            deallocate (a_pdswheat           )
            deallocate (a_pdwwheat           )
            deallocate (a_pdsoybean          )
            deallocate (a_pdcotton           )
            deallocate (a_pdrice1            )
            deallocate (a_pdrice2            )
            deallocate (a_pdsugarcane        )
            deallocate (a_plantdate          )
            deallocate (a_fertnitro_corn     )
            deallocate (a_fertnitro_swheat   )
            deallocate (a_fertnitro_wwheat   )
            deallocate (a_fertnitro_soybean  )
            deallocate (a_fertnitro_cotton   )
            deallocate (a_fertnitro_rice1    )
            deallocate (a_fertnitro_rice2    )
            deallocate (a_fertnitro_sugarcane)
            deallocate (a_cphase             )
            deallocate (a_hui                )
            deallocate (a_vf                 )
            deallocate (a_gddmaturity        )
            deallocate (a_gddplant           )
            deallocate (a_cropprod1c         )
            deallocate (a_cropprod1c_loss    )
            deallocate (a_cropseedc_deficit  )
            deallocate (a_grainc_to_cropprodc)
            deallocate (a_grainc_to_seed     )
            deallocate (a_fert_to_sminn      )
#endif
            deallocate (a_ndep_to_sminn      )
#ifdef Fire
            deallocate (a_abm                )
            deallocate (a_gdp                )
            deallocate (a_peatf              )
            deallocate (a_hdm                )
            deallocate (a_lnfm               )
#endif
#endif
#ifdef OzoneStress
            deallocate (a_ozone              )
#endif

            deallocate (a_t_soisno    )    
            deallocate (a_wliq_soisno )
            deallocate (a_wice_soisno )
            deallocate (a_h2osoi      ) 
            deallocate (a_rootr       )
            deallocate (a_BD_all      )
            deallocate (a_wfc         )
            deallocate (a_OM_density  )
#ifdef PLANT_HYDRAULIC_STRESS
            deallocate (a_vegwp       )
#endif
            deallocate (a_t_lake      ) 
            deallocate (a_lake_icefrac) 
#ifdef BGC
            deallocate (a_litr1c_vr   )
            deallocate (a_litr2c_vr   )
            deallocate (a_litr3c_vr   )
            deallocate (a_soil1c_vr   )
            deallocate (a_soil2c_vr   )
            deallocate (a_soil3c_vr   )
            deallocate (a_cwdc_vr     )
            deallocate (a_litr1n_vr   )
            deallocate (a_litr2n_vr   )
            deallocate (a_litr3n_vr   )
            deallocate (a_soil1n_vr   )
            deallocate (a_soil2n_vr   )
            deallocate (a_soil3n_vr   )
            deallocate (a_cwdn_vr     )
            deallocate (a_sminn_vr    )
            deallocate (decomp_vr_tmp )
#endif

            deallocate (a_ustar     ) 
            deallocate (a_tstar     )
            deallocate (a_qstar     )
            deallocate (a_zol       )
            deallocate (a_rib       )
            deallocate (a_fm        )
            deallocate (a_fh        )
            deallocate (a_fq        )

            deallocate (a_us10m     ) 
            deallocate (a_vs10m     ) 
            deallocate (a_fm10m     ) 

            deallocate (a_sr        )
            deallocate (a_solvd     )
            deallocate (a_solvi     )
            deallocate (a_solnd     )
            deallocate (a_solni     )
            deallocate (a_srvd      )
            deallocate (a_srvi      )
            deallocate (a_srnd      )
            deallocate (a_srni      )
            deallocate (a_solvdln   )
            deallocate (a_solviln   )
            deallocate (a_solndln   )
            deallocate (a_solniln   )
            deallocate (a_srvdln    )
            deallocate (a_srviln    )
            deallocate (a_srndln    )
            deallocate (a_srniln    )

            deallocate (nac_ln      )

         end if
      end if

   end subroutine deallocate_acc_fluxes

   !-----------------------
   SUBROUTINE FLUSH_acc_fluxes ()

      use spmd_task
      use mod_landpatch, only : numpatch
      use GlobalVars,    only : spval 
      implicit none

      if (p_is_worker) then

         nac = 0

         if (numpatch > 0) then

            ! flush the Fluxes for accumulation
            a_us     (:) = spval
            a_vs     (:) = spval
            a_t      (:) = spval
            a_q      (:) = spval
            a_prc    (:) = spval
            a_prl    (:) = spval
            a_pbot   (:) = spval
            a_frl    (:) = spval
            a_solarin(:) = spval

            a_taux    (:) = spval
            a_tauy    (:) = spval
            a_fsena   (:) = spval
            a_lfevpa  (:) = spval
            a_fevpa   (:) = spval
            a_fsenl   (:) = spval
            a_fevpl   (:) = spval
            a_etr     (:) = spval
            a_fseng   (:) = spval
            a_fevpg   (:) = spval
            a_fgrnd   (:) = spval
            a_sabvsun (:) = spval
            a_sabvsha (:) = spval
            a_sabg    (:) = spval
            a_olrg    (:) = spval
            a_rnet    (:) = spval
            a_xerr    (:) = spval
            a_zerr    (:) = spval
            a_rsur    (:) = spval
            a_rnof    (:) = spval
            a_qintr   (:) = spval
            a_qinfl   (:) = spval
            a_qdrip   (:) = spval
            a_rstfacsun(:) = spval
            a_rstfacsha(:) = spval
            a_gs_sun   (:) = spval
            a_gs_sha   (:) = spval

            a_dpond   (:) = spval
            a_zwt     (:) = spval
            a_wa      (:) = spval
            a_wat     (:) = spval
            a_assim   (:) = spval
            a_respc   (:) = spval

            a_qcharge (:) = spval

            a_t_grnd  (:) = spval
            a_tleaf   (:) = spval
            a_ldew_rain(:) = spval
            a_ldew_snow(:) = spval
            a_ldew    (:) = spval
            a_scv     (:) = spval
            a_snowdp  (:) = spval
            a_fsno    (:) = spval
            a_sigf    (:) = spval
            a_green   (:) = spval
            a_lai     (:) = spval
            a_laisun  (:) = spval
            a_laisha  (:) = spval
            a_sai     (:) = spval
            
            a_alb   (:,:,:) = spval

            a_emis      (:) = spval
            a_z0m       (:) = spval
            a_trad      (:) = spval
            a_tref      (:) = spval
            a_qref      (:) = spval 
            a_rain      (:) = spval
            a_snow      (:) = spval
#ifdef BGC
            a_leafc              (:) = spval
            a_leafc_storage      (:) = spval
            a_leafc_xfer         (:) = spval
            a_frootc             (:) = spval
            a_frootc_storage     (:) = spval
            a_frootc_xfer        (:) = spval
            a_livestemc          (:) = spval
            a_livestemc_storage  (:) = spval
            a_livestemc_xfer     (:) = spval
            a_deadstemc          (:) = spval
            a_deadstemc_storage  (:) = spval
            a_deadstemc_xfer     (:) = spval
            a_livecrootc         (:) = spval
            a_livecrootc_storage (:) = spval
            a_livecrootc_xfer    (:) = spval
            a_deadcrootc         (:) = spval
            a_deadcrootc_storage (:) = spval
            a_deadcrootc_xfer    (:) = spval
            a_grainc             (:) = spval
            a_grainc_storage     (:) = spval
            a_grainc_xfer        (:) = spval
            a_leafn              (:) = spval
            a_leafn_storage      (:) = spval
            a_leafn_xfer         (:) = spval
            a_frootn             (:) = spval
            a_frootn_storage     (:) = spval
            a_frootn_xfer        (:) = spval
            a_livestemn          (:) = spval
            a_livestemn_storage  (:) = spval
            a_livestemn_xfer     (:) = spval
            a_deadstemn          (:) = spval
            a_deadstemn_storage  (:) = spval
            a_deadstemn_xfer     (:) = spval
            a_livecrootn         (:) = spval
            a_livecrootn_storage (:) = spval
            a_livecrootn_xfer    (:) = spval
            a_deadcrootn         (:) = spval
            a_deadcrootn_storage (:) = spval
            a_deadcrootn_xfer    (:) = spval
            a_grainn             (:) = spval
            a_grainn_storage     (:) = spval
            a_grainn_xfer        (:) = spval
            a_retransn           (:) = spval
            a_gpp                (:) = spval
            a_downreg            (:) = spval
            a_ar                 (:) = spval
            a_cwdprod            (:) = spval
            a_cwddecomp          (:) = spval
            a_hr                 (:) = spval
            a_fpg                (:) = spval
            a_fpi                (:) = spval
            a_gpp_enftemp        (:) = spval
            a_gpp_enfboreal      (:) = spval
            a_gpp_dnfboreal      (:) = spval
            a_gpp_ebftrop        (:) = spval
            a_gpp_ebftemp        (:) = spval
            a_gpp_dbftrop        (:) = spval
            a_gpp_dbftemp        (:) = spval
            a_gpp_dbfboreal      (:) = spval
            a_gpp_ebstemp        (:) = spval
            a_gpp_dbstemp        (:) = spval
            a_gpp_dbsboreal      (:) = spval
            a_gpp_c3arcgrass     (:) = spval
            a_gpp_c3grass        (:) = spval
            a_gpp_c4grass        (:) = spval
            a_leafc_enftemp      (:) = spval
            a_leafc_enfboreal    (:) = spval
            a_leafc_dnfboreal    (:) = spval
            a_leafc_ebftrop      (:) = spval
            a_leafc_ebftemp      (:) = spval
            a_leafc_dbftrop      (:) = spval
            a_leafc_dbftemp      (:) = spval
            a_leafc_dbfboreal    (:) = spval
            a_leafc_ebstemp      (:) = spval
            a_leafc_dbstemp      (:) = spval
            a_leafc_dbsboreal    (:) = spval
            a_leafc_c3arcgrass   (:) = spval
            a_leafc_c3grass      (:) = spval
            a_leafc_c4grass      (:) = spval
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
            a_assim_RuBP_sun_enftemp        (:) = spval !1
            a_assim_RuBP_sun_enfboreal      (:) = spval !2
            a_assim_RuBP_sun_dnfboreal      (:) = spval !3
            a_assim_RuBP_sun_ebftrop        (:) = spval !4
            a_assim_RuBP_sun_ebftemp        (:) = spval !5
            a_assim_RuBP_sun_dbftrop        (:) = spval !6
            a_assim_RuBP_sun_dbftemp        (:) = spval !7
            a_assim_RuBP_sun_dbfboreal      (:) = spval !8
            a_assim_RuBP_sun_ebstemp        (:) = spval !9
            a_assim_RuBP_sun_dbstemp        (:) = spval !10
            a_assim_RuBP_sun_dbsboreal      (:) = spval !11
            a_assim_RuBP_sun_c3arcgrass     (:) = spval !12
            a_assim_RuBP_sun_c3grass        (:) = spval !13
            a_assim_RuBP_sun_c4grass        (:) = spval !14
            a_assim_RuBP_sha_enftemp        (:) = spval !1
            a_assim_RuBP_sha_enfboreal      (:) = spval !2
            a_assim_RuBP_sha_dnfboreal      (:) = spval !3
            a_assim_RuBP_sha_ebftrop        (:) = spval !4
            a_assim_RuBP_sha_ebftemp        (:) = spval !5
            a_assim_RuBP_sha_dbftrop        (:) = spval !6
            a_assim_RuBP_sha_dbftemp        (:) = spval !7
            a_assim_RuBP_sha_dbfboreal      (:) = spval !8
            a_assim_RuBP_sha_ebstemp        (:) = spval !9
            a_assim_RuBP_sha_dbstemp        (:) = spval !10
            a_assim_RuBP_sha_dbsboreal      (:) = spval !11
            a_assim_RuBP_sha_c3arcgrass     (:) = spval !12
            a_assim_RuBP_sha_c3grass        (:) = spval !13
            a_assim_RuBP_sha_c4grass        (:) = spval !14
            a_assim_Rubisco_sun_enftemp        (:) = spval !1
            a_assim_Rubisco_sun_enfboreal      (:) = spval !2
            a_assim_Rubisco_sun_dnfboreal      (:) = spval !3
            a_assim_Rubisco_sun_ebftrop        (:) = spval !4
            a_assim_Rubisco_sun_ebftemp        (:) = spval !5
            a_assim_Rubisco_sun_dbftrop        (:) = spval !6
            a_assim_Rubisco_sun_dbftemp        (:) = spval !7
            a_assim_Rubisco_sun_dbfboreal      (:) = spval !8
            a_assim_Rubisco_sun_ebstemp        (:) = spval !9
            a_assim_Rubisco_sun_dbstemp        (:) = spval !10
            a_assim_Rubisco_sun_dbsboreal      (:) = spval !11
            a_assim_Rubisco_sun_c3arcgrass     (:) = spval !12
            a_assim_Rubisco_sun_c3grass        (:) = spval !13
            a_assim_Rubisco_sun_c4grass        (:) = spval !14
            a_assim_Rubisco_sha_enftemp        (:) = spval !1
            a_assim_Rubisco_sha_enfboreal      (:) = spval !2
            a_assim_Rubisco_sha_dnfboreal      (:) = spval !3
            a_assim_Rubisco_sha_ebftrop        (:) = spval !4
            a_assim_Rubisco_sha_ebftemp        (:) = spval !5
            a_assim_Rubisco_sha_dbftrop        (:) = spval !6
            a_assim_Rubisco_sha_dbftemp        (:) = spval !7
            a_assim_Rubisco_sha_dbfboreal      (:) = spval !8
            a_assim_Rubisco_sha_ebstemp        (:) = spval !9
            a_assim_Rubisco_sha_dbstemp        (:) = spval !10
            a_assim_Rubisco_sha_dbsboreal      (:) = spval !11
            a_assim_Rubisco_sha_c3arcgrass     (:) = spval !12
            a_assim_Rubisco_sha_c3grass        (:) = spval !13
            a_assim_Rubisco_sha_c4grass        (:) = spval !14
            a_assimsun_enftemp        (:) = spval !1
            a_assimsun_enfboreal      (:) = spval !2
            a_assimsun_dnfboreal      (:) = spval !3
            a_assimsun_ebftrop        (:) = spval !4
            a_assimsun_ebftemp        (:) = spval !5
            a_assimsun_dbftrop        (:) = spval !6
            a_assimsun_dbftemp        (:) = spval !7
            a_assimsun_dbfboreal      (:) = spval !8
            a_assimsun_ebstemp        (:) = spval !9
            a_assimsun_dbstemp        (:) = spval !10
            a_assimsun_dbsboreal      (:) = spval !11
            a_assimsun_c3arcgrass     (:) = spval !12
            a_assimsun_c3grass        (:) = spval !13
            a_assimsun_c4grass        (:) = spval !14
            a_assimsha_enftemp        (:) = spval !1
            a_assimsha_enfboreal      (:) = spval !2
            a_assimsha_dnfboreal      (:) = spval !3
            a_assimsha_ebftrop        (:) = spval !4
            a_assimsha_ebftemp        (:) = spval !5
            a_assimsha_dbftrop        (:) = spval !6
            a_assimsha_dbftemp        (:) = spval !7
            a_assimsha_dbfboreal      (:) = spval !8
            a_assimsha_ebstemp        (:) = spval !9
            a_assimsha_dbstemp        (:) = spval !10
            a_assimsha_dbsboreal      (:) = spval !11
            a_assimsha_c3arcgrass     (:) = spval !12
            a_assimsha_c3grass        (:) = spval !13
            a_assimsha_c4grass        (:) = spval !14
            a_etrsun_enftemp        (:) = spval !1
            a_etrsun_enfboreal      (:) = spval !2
            a_etrsun_dnfboreal      (:) = spval !3
            a_etrsun_ebftrop        (:) = spval !4
            a_etrsun_ebftemp        (:) = spval !5
            a_etrsun_dbftrop        (:) = spval !6
            a_etrsun_dbftemp        (:) = spval !7
            a_etrsun_dbfboreal      (:) = spval !8
            a_etrsun_ebstemp        (:) = spval !9
            a_etrsun_dbstemp        (:) = spval !10
            a_etrsun_dbsboreal      (:) = spval !11
            a_etrsun_c3arcgrass     (:) = spval !12
            a_etrsun_c3grass        (:) = spval !13
            a_etrsun_c4grass        (:) = spval !14
            a_etrsha_enftemp        (:) = spval !1
            a_etrsha_enfboreal      (:) = spval !2
            a_etrsha_dnfboreal      (:) = spval !3
            a_etrsha_ebftrop        (:) = spval !4
            a_etrsha_ebftemp        (:) = spval !5
            a_etrsha_dbftrop        (:) = spval !6
            a_etrsha_dbftemp        (:) = spval !7
            a_etrsha_dbfboreal      (:) = spval !8
            a_etrsha_ebstemp        (:) = spval !9
            a_etrsha_dbstemp        (:) = spval !10
            a_etrsha_dbsboreal      (:) = spval !11
            a_etrsha_c3arcgrass     (:) = spval !12
            a_etrsha_c3grass        (:) = spval !13
            a_etrsha_c4grass        (:) = spval !14
            a_cisun_enftemp        (:) = spval !1
            a_cisun_enfboreal      (:) = spval !2
            a_cisun_dnfboreal      (:) = spval !3
            a_cisun_ebftrop        (:) = spval !4
            a_cisun_ebftemp        (:) = spval !5
            a_cisun_dbftrop        (:) = spval !6
            a_cisun_dbftemp        (:) = spval !7
            a_cisun_dbfboreal      (:) = spval !8
            a_cisun_ebstemp        (:) = spval !9
            a_cisun_dbstemp        (:) = spval !10
            a_cisun_dbsboreal      (:) = spval !11
            a_cisun_c3arcgrass     (:) = spval !12
            a_cisun_c3grass        (:) = spval !13
            a_cisun_c4grass        (:) = spval !14
            a_cisha_enftemp        (:) = spval !1
            a_cisha_enfboreal      (:) = spval !2
            a_cisha_dnfboreal      (:) = spval !3
            a_cisha_ebftrop        (:) = spval !4
            a_cisha_ebftemp        (:) = spval !5
            a_cisha_dbftrop        (:) = spval !6
            a_cisha_dbftemp        (:) = spval !7
            a_cisha_dbfboreal      (:) = spval !8
            a_cisha_ebstemp        (:) = spval !9
            a_cisha_dbstemp        (:) = spval !10
            a_cisha_dbsboreal      (:) = spval !11
            a_cisha_c3arcgrass     (:) = spval !12
            a_cisha_c3grass        (:) = spval !13
            a_cisha_c4grass        (:) = spval !14
            a_essun_enftemp        (:) = spval !1
            a_essun_enfboreal      (:) = spval !2
            a_essun_dnfboreal      (:) = spval !3
            a_essun_ebftrop        (:) = spval !4
            a_essun_ebftemp        (:) = spval !5
            a_essun_dbftrop        (:) = spval !6
            a_essun_dbftemp        (:) = spval !7
            a_essun_dbfboreal      (:) = spval !8
            a_essun_ebstemp        (:) = spval !9
            a_essun_dbstemp        (:) = spval !10
            a_essun_dbsboreal      (:) = spval !11
            a_essun_c3arcgrass     (:) = spval !12
            a_essun_c3grass        (:) = spval !13
            a_essun_c4grass        (:) = spval !14
            a_essha_enftemp        (:) = spval !1
            a_essha_enfboreal      (:) = spval !2
            a_essha_dnfboreal      (:) = spval !3
            a_essha_ebftrop        (:) = spval !4
            a_essha_ebftemp        (:) = spval !5
            a_essha_dbftrop        (:) = spval !6
            a_essha_dbftemp        (:) = spval !7
            a_essha_dbfboreal      (:) = spval !8
            a_essha_ebstemp        (:) = spval !9
            a_essha_dbstemp        (:) = spval !10
            a_essha_dbsboreal      (:) = spval !11
            a_essha_c3arcgrass     (:) = spval !12
            a_essha_c3grass        (:) = spval !13
            a_essha_c4grass        (:) = spval !14
            a_gssun_enftemp        (:) = spval !1
            a_gssun_enfboreal      (:) = spval !2
            a_gssun_dnfboreal      (:) = spval !3
            a_gssun_ebftrop        (:) = spval !4
            a_gssun_ebftemp        (:) = spval !5
            a_gssun_dbftrop        (:) = spval !6
            a_gssun_dbftemp        (:) = spval !7
            a_gssun_dbfboreal      (:) = spval !8
            a_gssun_ebstemp        (:) = spval !9
            a_gssun_dbstemp        (:) = spval !10
            a_gssun_dbsboreal      (:) = spval !11
            a_gssun_c3arcgrass     (:) = spval !12
            a_gssun_c3grass        (:) = spval !13
            a_gssun_c4grass        (:) = spval !14
            a_gssha_enftemp        (:) = spval !1
            a_gssha_enfboreal      (:) = spval !2
            a_gssha_dnfboreal      (:) = spval !3
            a_gssha_ebftrop        (:) = spval !4
            a_gssha_ebftemp        (:) = spval !5
            a_gssha_dbftrop        (:) = spval !6
            a_gssha_dbftemp        (:) = spval !7
            a_gssha_dbfboreal      (:) = spval !8
            a_gssha_ebstemp        (:) = spval !9
            a_gssha_dbstemp        (:) = spval !10
            a_gssha_dbsboreal      (:) = spval !11
            a_gssha_c3arcgrass     (:) = spval !12
            a_gssha_c3grass        (:) = spval !13
            a_gssha_c4grass        (:) = spval !14
            a_gammasun_enftemp        (:) = spval !1
            a_gammasun_enfboreal      (:) = spval !2
            a_gammasun_dnfboreal      (:) = spval !3
            a_gammasun_ebftrop        (:) = spval !4
            a_gammasun_ebftemp        (:) = spval !5
            a_gammasun_dbftrop        (:) = spval !6
            a_gammasun_dbftemp        (:) = spval !7
            a_gammasun_dbfboreal      (:) = spval !8
            a_gammasun_ebstemp        (:) = spval !9
            a_gammasun_dbstemp        (:) = spval !10
            a_gammasun_dbsboreal      (:) = spval !11
            a_gammasun_c3arcgrass     (:) = spval !12
            a_gammasun_c3grass        (:) = spval !13
            a_gammasun_c4grass        (:) = spval !14
            a_gammasha_enftemp        (:) = spval !1
            a_gammasha_enfboreal      (:) = spval !2
            a_gammasha_dnfboreal      (:) = spval !3
            a_gammasha_ebftrop        (:) = spval !4
            a_gammasha_ebftemp        (:) = spval !5
            a_gammasha_dbftrop        (:) = spval !6
            a_gammasha_dbftemp        (:) = spval !7
            a_gammasha_dbfboreal      (:) = spval !8
            a_gammasha_ebstemp        (:) = spval !9
            a_gammasha_dbstemp        (:) = spval !10
            a_gammasha_dbsboreal      (:) = spval !11
            a_gammasha_c3arcgrass     (:) = spval !12
            a_gammasha_c3grass        (:) = spval !13
            a_gammasha_c4grass        (:) = spval !14
            a_RuBPlimitfrac_sun_enftemp          (:) = spval
            a_RuBPlimitfrac_sun_enfboreal        (:) = spval
            a_RuBPlimitfrac_sun_dnfboreal        (:) = spval
            a_RuBPlimitfrac_sun_ebftrop          (:) = spval
            a_RuBPlimitfrac_sun_ebftemp          (:) = spval
            a_RuBPlimitfrac_sun_dbftrop          (:) = spval
            a_RuBPlimitfrac_sun_dbftemp          (:) = spval
            a_RuBPlimitfrac_sun_dbfboreal        (:) = spval
            a_RuBPlimitfrac_sun_ebstemp          (:) = spval
            a_RuBPlimitfrac_sun_dbstemp          (:) = spval
            a_RuBPlimitfrac_sun_dbsboreal        (:) = spval
            a_RuBPlimitfrac_sun_c3arcgrass       (:) = spval
            a_RuBPlimitfrac_sun_c3grass          (:) = spval
            a_RuBPlimitfrac_sun_c4grass          (:) = spval
            a_RuBPlimitfrac_sha_enftemp          (:) = spval
            a_RuBPlimitfrac_sha_enfboreal        (:) = spval
            a_RuBPlimitfrac_sha_dnfboreal        (:) = spval
            a_RuBPlimitfrac_sha_ebftrop          (:) = spval
            a_RuBPlimitfrac_sha_ebftemp          (:) = spval
            a_RuBPlimitfrac_sha_dbftrop          (:) = spval
            a_RuBPlimitfrac_sha_dbftemp          (:) = spval
            a_RuBPlimitfrac_sha_dbfboreal        (:) = spval
            a_RuBPlimitfrac_sha_ebstemp          (:) = spval
            a_RuBPlimitfrac_sha_dbstemp          (:) = spval
            a_RuBPlimitfrac_sha_dbsboreal        (:) = spval
            a_RuBPlimitfrac_sha_c3arcgrass       (:) = spval
            a_RuBPlimitfrac_sha_c3grass          (:) = spval
            a_RuBPlimitfrac_sha_c4grass          (:) = spval
            a_Rubiscolimitfrac_sun_enftemp          (:) = spval
            a_Rubiscolimitfrac_sun_enfboreal        (:) = spval
            a_Rubiscolimitfrac_sun_dnfboreal        (:) = spval
            a_Rubiscolimitfrac_sun_ebftrop          (:) = spval
            a_Rubiscolimitfrac_sun_ebftemp          (:) = spval
            a_Rubiscolimitfrac_sun_dbftrop          (:) = spval
            a_Rubiscolimitfrac_sun_dbftemp          (:) = spval
            a_Rubiscolimitfrac_sun_dbfboreal        (:) = spval
            a_Rubiscolimitfrac_sun_ebstemp          (:) = spval
            a_Rubiscolimitfrac_sun_dbstemp          (:) = spval
            a_Rubiscolimitfrac_sun_dbsboreal        (:) = spval
            a_Rubiscolimitfrac_sun_c3arcgrass       (:) = spval
            a_Rubiscolimitfrac_sun_c3grass          (:) = spval
            a_Rubiscolimitfrac_sun_c4grass          (:) = spval
            a_Rubiscolimitfrac_sha_enftemp          (:) = spval
            a_Rubiscolimitfrac_sha_enfboreal        (:) = spval
            a_Rubiscolimitfrac_sha_dnfboreal        (:) = spval
            a_Rubiscolimitfrac_sha_ebftrop          (:) = spval
            a_Rubiscolimitfrac_sha_ebftemp          (:) = spval
            a_Rubiscolimitfrac_sha_dbftrop          (:) = spval
            a_Rubiscolimitfrac_sha_dbftemp          (:) = spval
            a_Rubiscolimitfrac_sha_dbfboreal        (:) = spval
            a_Rubiscolimitfrac_sha_ebstemp          (:) = spval
            a_Rubiscolimitfrac_sha_dbstemp          (:) = spval
            a_Rubiscolimitfrac_sha_dbsboreal        (:) = spval
            a_Rubiscolimitfrac_sha_c3arcgrass       (:) = spval
            a_Rubiscolimitfrac_sha_c3grass          (:) = spval
            a_Rubiscolimitfrac_sha_c4grass          (:) = spval
            a_Sinklimitfrac_sun_enftemp          (:) = spval
            a_Sinklimitfrac_sun_enfboreal        (:) = spval
            a_Sinklimitfrac_sun_dnfboreal        (:) = spval
            a_Sinklimitfrac_sun_ebftrop          (:) = spval
            a_Sinklimitfrac_sun_ebftemp          (:) = spval
            a_Sinklimitfrac_sun_dbftrop          (:) = spval
            a_Sinklimitfrac_sun_dbftemp          (:) = spval
            a_Sinklimitfrac_sun_dbfboreal        (:) = spval
            a_Sinklimitfrac_sun_ebstemp          (:) = spval
            a_Sinklimitfrac_sun_dbstemp          (:) = spval
            a_Sinklimitfrac_sun_dbsboreal        (:) = spval
            a_Sinklimitfrac_sun_c3arcgrass       (:) = spval
            a_Sinklimitfrac_sun_c3grass          (:) = spval
            a_Sinklimitfrac_sun_c4grass          (:) = spval
            a_Sinklimitfrac_sha_enftemp          (:) = spval
            a_Sinklimitfrac_sha_enfboreal        (:) = spval
            a_Sinklimitfrac_sha_dnfboreal        (:) = spval
            a_Sinklimitfrac_sha_ebftrop          (:) = spval
            a_Sinklimitfrac_sha_ebftemp          (:) = spval
            a_Sinklimitfrac_sha_dbftrop          (:) = spval
            a_Sinklimitfrac_sha_dbftemp          (:) = spval
            a_Sinklimitfrac_sha_dbfboreal        (:) = spval
            a_Sinklimitfrac_sha_ebstemp          (:) = spval
            a_Sinklimitfrac_sha_dbstemp          (:) = spval
            a_Sinklimitfrac_sha_dbsboreal        (:) = spval
            a_Sinklimitfrac_sha_c3arcgrass       (:) = spval
            a_Sinklimitfrac_sha_c3grass          (:) = spval
            a_Sinklimitfrac_sha_c4grass          (:) = spval
            a_rstfacsun_enftemp        (:) = spval  !1
            a_rstfacsun_enfboreal      (:) = spval  !2
            a_rstfacsun_dnfboreal      (:) = spval  !3
            a_rstfacsun_ebftrop        (:) = spval  !4
            a_rstfacsun_ebftemp        (:) = spval  !5
            a_rstfacsun_dbftrop        (:) = spval  !6
            a_rstfacsun_dbftemp        (:) = spval  !7
            a_rstfacsun_dbfboreal      (:) = spval  !8
            a_rstfacsun_ebstemp        (:) = spval  !9
            a_rstfacsun_dbstemp        (:) = spval  !10
            a_rstfacsun_dbsboreal      (:) = spval  !11
            a_rstfacsun_c3arcgrass     (:) = spval  !12
            a_rstfacsun_c3grass        (:) = spval  !13
            a_rstfacsun_c4grass        (:) = spval  !14
            a_rstfacsha_enftemp        (:) = spval  !1
            a_rstfacsha_enfboreal      (:) = spval  !2
            a_rstfacsha_dnfboreal      (:) = spval  !3
            a_rstfacsha_ebftrop        (:) = spval  !4
            a_rstfacsha_ebftemp        (:) = spval  !5
            a_rstfacsha_dbftrop        (:) = spval  !6
            a_rstfacsha_dbftemp        (:) = spval  !7
            a_rstfacsha_dbfboreal      (:) = spval  !8
            a_rstfacsha_ebstemp        (:) = spval  !9
            a_rstfacsha_dbstemp        (:) = spval  !10
            a_rstfacsha_dbsboreal      (:) = spval  !11
            a_rstfacsha_c3arcgrass     (:) = spval  !12
            a_rstfacsha_c3grass        (:) = spval  !13
            a_rstfacsha_c4grass        (:) = spval  !14
            a_lambdasun_enftemp        (:) = spval  !1
            a_lambdasun_enfboreal      (:) = spval  !2
            a_lambdasun_dnfboreal      (:) = spval  !3
            a_lambdasun_ebftrop        (:) = spval  !4
            a_lambdasun_ebftemp        (:) = spval  !5
            a_lambdasun_dbftrop        (:) = spval  !6
            a_lambdasun_dbftemp        (:) = spval  !7
            a_lambdasun_dbfboreal      (:) = spval  !8
            a_lambdasun_ebstemp        (:) = spval  !9
            a_lambdasun_dbstemp        (:) = spval  !10
            a_lambdasun_dbsboreal      (:) = spval  !11
            a_lambdasun_c3arcgrass     (:) = spval  !12
            a_lambdasun_c3grass        (:) = spval  !13
            a_lambdasun_c4grass        (:) = spval  !14
            a_lambdasha_enftemp        (:) = spval  !1
            a_lambdasha_enfboreal      (:) = spval  !2
            a_lambdasha_dnfboreal      (:) = spval  !3
            a_lambdasha_ebftrop        (:) = spval  !4
            a_lambdasha_ebftemp        (:) = spval  !5
            a_lambdasha_dbftrop        (:) = spval  !6
            a_lambdasha_dbftemp        (:) = spval  !7
            a_lambdasha_dbfboreal      (:) = spval  !8
            a_lambdasha_ebstemp        (:) = spval  !9
            a_lambdasha_dbstemp        (:) = spval  !10
            a_lambdasha_dbsboreal      (:) = spval  !11
            a_lambdasha_c3arcgrass     (:) = spval  !12
            a_lambdasha_c3grass        (:) = spval  !13
            a_lambdasha_c4grass        (:) = spval  !14
            a_lambda                   (:) = spval  !1
#endif
#endif
#ifdef NITRIF
            a_O2_DECOMP_DEPTH_UNSAT (:,:) = spval
            a_CONC_O2_UNSAT         (:,:) = spval
#endif
#ifdef CROP
            a_pdcorn             (:) = spval
            a_pdswheat           (:) = spval
            a_pdwwheat           (:) = spval
            a_pdsoybean          (:) = spval
            a_pdcotton           (:) = spval
            a_pdrice1            (:) = spval
            a_pdrice2            (:) = spval
            a_pdsugarcane        (:) = spval
            a_plantdate          (:) = spval
            a_fertnitro_corn     (:) = spval
            a_fertnitro_swheat   (:) = spval
            a_fertnitro_wwheat   (:) = spval
            a_fertnitro_soybean  (:) = spval
            a_fertnitro_cotton   (:) = spval
            a_fertnitro_rice1    (:) = spval
            a_fertnitro_rice2    (:) = spval
            a_fertnitro_sugarcane(:) = spval
            a_cphase             (:) = spval
            a_vf                 (:) = spval
            a_gddmaturity        (:) = spval
            a_gddplant           (:) = spval
            a_hui                (:) = spval
            a_cropprod1c         (:) = spval
            a_cropprod1c_loss    (:) = spval
            a_cropseedc_deficit  (:) = spval
            a_grainc_to_cropprodc(:) = spval
            a_grainc_to_seed     (:) = spval
            a_fert_to_sminn      (:) = spval
#endif
            a_ndep_to_sminn      (:) = spval
#ifdef Fire
            a_abm                (:) = spval
            a_gdp                (:) = spval
            a_peatf              (:) = spval
            a_hdm                (:) = spval
            a_lnfm               (:) = spval
#endif            
#endif
#ifdef OzoneStress
            a_ozone              (:) = spval
#endif

            a_t_soisno     (:,:) = spval 
            a_wliq_soisno  (:,:) = spval
            a_wice_soisno  (:,:) = spval
            a_h2osoi       (:,:) = spval
            a_rootr        (:,:) = spval
            a_BD_all       (:,:) = spval
            a_wfc          (:,:) = spval
            a_OM_density   (:,:) = spval
#ifdef PLANT_HYDRAULIC_STRESS
            a_vegwp        (:,:) = spval
#endif
            a_t_lake       (:,:) = spval
            a_lake_icefrac (:,:) = spval
#ifdef BGC
            a_litr1c_vr    (:,:) = spval
            a_litr2c_vr    (:,:) = spval
            a_litr3c_vr    (:,:) = spval
            a_soil1c_vr    (:,:) = spval
            a_soil2c_vr    (:,:) = spval
            a_soil3c_vr    (:,:) = spval
            a_cwdc_vr      (:,:) = spval
            a_litr1n_vr    (:,:) = spval
            a_litr2n_vr    (:,:) = spval
            a_litr3n_vr    (:,:) = spval
            a_soil1n_vr    (:,:) = spval
            a_soil2n_vr    (:,:) = spval
            a_soil3n_vr    (:,:) = spval
            a_cwdn_vr      (:,:) = spval
            a_sminn_vr     (:,:) = spval
#endif

            a_ustar (:) = spval 
            a_tstar (:) = spval
            a_qstar (:) = spval
            a_zol   (:) = spval
            a_rib   (:) = spval
            a_fm    (:) = spval
            a_fh    (:) = spval
            a_fq    (:) = spval

            a_us10m (:) = spval
            a_vs10m (:) = spval
            a_fm10m (:) = spval

            a_sr       (:) = spval
            a_solvd    (:) = spval
            a_solvi    (:) = spval
            a_solnd    (:) = spval
            a_solni    (:) = spval
            a_srvd     (:) = spval
            a_srvi     (:) = spval
            a_srnd     (:) = spval
            a_srni     (:) = spval
            a_solvdln  (:) = spval
            a_solviln  (:) = spval
            a_solndln  (:) = spval
            a_solniln  (:) = spval
            a_srvdln   (:) = spval
            a_srviln   (:) = spval
            a_srndln   (:) = spval
            a_srniln   (:) = spval

            nac_ln  (:) = 0

         end if
      end if

   END SUBROUTINE FLUSH_acc_fluxes

   SUBROUTINE accumulate_fluxes 
      ! ----------------------------------------------------------------------
      ! perfrom the grid average mapping: average a subgrid input 1d vector 
      ! of length numpatch to a output 2d array of length [ghist%xcnt,ghist%ycnt]
      !
      ! Created by Yongjiu Dai, 03/2014
      !---------------------------------------------------------------------

      use precision
      use spmd_task
      use mod_landpatch,     only : numpatch
      use PhysicalConstants, only : vonkar, stefnc, cpair, rgas, grav
      use MOD_TimeInvariants
      use MOD_TimeVariables
      use MOD_1D_Forcing
      use MOD_1D_Fluxes
      use FRICTION_VELOCITY
      use mod_colm_debug
      use GlobalVars
      

      IMPLICIT NONE

      ! Local Variables

      real(r8), allocatable :: r_trad  (:)

      real(r8), allocatable :: r_ustar (:)
      real(r8), allocatable :: r_tstar (:)
      real(r8), allocatable :: r_qstar (:)
      real(r8), allocatable :: r_zol   (:)
      real(r8), allocatable :: r_rib   (:)
      real(r8), allocatable :: r_fm    (:)
      real(r8), allocatable :: r_fh    (:)
      real(r8), allocatable :: r_fq    (:)

      real(r8), allocatable :: r_us10m (:)
      real(r8), allocatable :: r_vs10m (:)
      real(r8), allocatable :: r_fm10m (:)

      !---------------------------------------------------------------------
      integer  ib, jb, i, j
      real(r8) rhoair,thm,th,thv,ur,displa_av,zldis,hgt_u,hgt_t,hgt_q
      real(r8) z0m_av,z0h_av,z0q_av,us,vs,tm,qm,psrf
      real(r8) obu,fh2m,fq2m
      real(r8) um,thvstar,beta,zii,wc,wc2

      if (p_is_worker) then
         if (numpatch > 0) then

            nac = nac + 1

            call acc1d (forc_us  , a_us  )
            call acc1d (forc_vs  , a_vs  )
            call acc1d (forc_t   , a_t   )
            call acc1d (forc_q   , a_q   )
            call acc1d (forc_prc , a_prc )
            call acc1d (forc_prl , a_prl )
            call acc1d (forc_pbot, a_pbot)
            call acc1d (forc_frl , a_frl )

            call acc1d (forc_sols,  a_solarin)
            call acc1d (forc_soll,  a_solarin)
            call acc1d (forc_solsd, a_solarin)
            call acc1d (forc_solld, a_solarin)

            call acc1d (taux    , a_taux   )
            call acc1d (tauy    , a_tauy   )
            call acc1d (fsena   , a_fsena  )
            call acc1d (lfevpa  , a_lfevpa )
            call acc1d (fevpa   , a_fevpa  )
            call acc1d (fsenl   , a_fsenl  )
            call acc1d (fevpl   , a_fevpl  )
            call acc1d (etr     , a_etr    )
            call acc1d (fseng   , a_fseng  )
            call acc1d (fevpg   , a_fevpg  )
            call acc1d (fgrnd   , a_fgrnd  )
            call acc1d (sabvsun , a_sabvsun)
            call acc1d (sabvsha , a_sabvsha)
            call acc1d (sabg    , a_sabg   )
            call acc1d (olrg    , a_olrg   )

            rnet = sabg + sabvsun + sabvsha - olrg + forc_frl
            call acc1d (rnet    , a_rnet   )

            call acc1d (xerr   , a_xerr   )
            call acc1d (zerr   , a_zerr   )
            call acc1d (rsur   , a_rsur   )
            call acc1d (rnof   , a_rnof   )
            call acc1d (qintr  , a_qintr  )
            call acc1d (qinfl  , a_qinfl  )
            call acc1d (qdrip  , a_qdrip  )
            call acc1d (rstfacsun , a_rstfacsun )
            call acc1d (rstfacsha , a_rstfacsha )
            call acc1d (gs_sun     , a_gs_sun )
            call acc1d (gs_sha     , a_gs_sha )

            call acc1d (dpond  , a_dpond  )
            call acc1d (zwt    , a_zwt    )
            call acc1d (wa     , a_wa     )
            call acc1d (wat    , a_wat    )
            call acc1d (assim  , a_assim  )
            call acc1d (respc  , a_respc  )

            call acc1d (qcharge, a_qcharge)

            call acc1d (t_grnd , a_t_grnd )
            call acc1d (tleaf  , a_tleaf  )
            call acc1d (ldew_rain, a_ldew_rain)
            call acc1d (ldew_snow, a_ldew_snow)
            call acc1d (ldew   , a_ldew   )
            call acc1d (scv    , a_scv    )
            call acc1d (snowdp , a_snowdp )
            call acc1d (fsno   , a_fsno   )
            call acc1d (sigf   , a_sigf   )
            call acc1d (green  , a_green  )
            call acc1d (lai    , a_lai    )
            call acc1d (laisun , a_laisun )
            call acc1d (laisha , a_laisha )
            call acc1d (sai    , a_sai    )

            call acc3d (alb    , a_alb    )

            call acc1d (emis   , a_emis   )
            call acc1d (z0m    , a_z0m    )

            allocate (r_trad (numpatch))
            do i = 1, numpatch
               r_trad(i) = (olrg(i)/stefnc)**0.25
            end do
            call acc1d (r_trad , a_trad   )
            deallocate (r_trad )

            call acc1d (tref   , a_tref   )
            call acc1d (qref   , a_qref   )

            call acc1d (forc_rain, a_rain )
            call acc1d (forc_snow, a_snow )

#ifdef BGC
            call acc1d (leafc              , a_leafc               )
            call acc1d (leafc_storage      , a_leafc_storage       )
            call acc1d (leafc_xfer         , a_leafc_xfer          )
            call acc1d (frootc             , a_frootc              )
            call acc1d (frootc_storage     , a_frootc_storage      )
            call acc1d (frootc_xfer        , a_frootc_xfer         )
            call acc1d (livestemc          , a_livestemc           )
            call acc1d (livestemc_storage  , a_livestemc_storage   )
            call acc1d (livestemc_xfer     , a_livestemc_xfer      )
            call acc1d (deadstemc          , a_deadstemc           )
            call acc1d (deadstemc_storage  , a_deadstemc_storage   )
            call acc1d (deadstemc_xfer     , a_deadstemc_xfer      )
            call acc1d (livecrootc         , a_livecrootc          )
            call acc1d (livecrootc_storage , a_livecrootc_storage  )
            call acc1d (livecrootc_xfer    , a_livecrootc_xfer     )
            call acc1d (deadcrootc         , a_deadcrootc          )
            call acc1d (deadcrootc_storage , a_deadcrootc_storage  )
            call acc1d (deadcrootc_xfer    , a_deadcrootc_xfer     )
            call acc1d (grainc             , a_grainc              )
            call acc1d (grainc_storage     , a_grainc_storage      )
            call acc1d (grainc_xfer        , a_grainc_xfer         )
            call acc1d (leafn              , a_leafn               )
            call acc1d (leafn_storage      , a_leafn_storage       )
            call acc1d (leafn_xfer         , a_leafn_xfer          )
            call acc1d (frootn             , a_frootn              )
            call acc1d (frootn_storage     , a_frootn_storage      )
            call acc1d (frootn_xfer        , a_frootn_xfer         )
            call acc1d (livestemn          , a_livestemn           )
            call acc1d (livestemn_storage  , a_livestemn_storage   )
            call acc1d (livestemn_xfer     , a_livestemn_xfer      )
            call acc1d (deadstemn          , a_deadstemn           )
            call acc1d (deadstemn_storage  , a_deadstemn_storage   )
            call acc1d (deadstemn_xfer     , a_deadstemn_xfer      )
            call acc1d (livecrootn         , a_livecrootn          )
            call acc1d (livecrootn_storage , a_livecrootn_storage  )
            call acc1d (livecrootn_xfer    , a_livecrootn_xfer     )
            call acc1d (deadcrootn         , a_deadcrootn          )
            call acc1d (deadcrootn_storage , a_deadcrootn_storage  )
            call acc1d (deadcrootn_xfer    , a_deadcrootn_xfer     )
            call acc1d (grainn             , a_grainn              )
            call acc1d (grainn_storage     , a_grainn_storage      )
            call acc1d (grainn_xfer        , a_grainn_xfer         )
            call acc1d (retransn           , a_retransn            )
            call acc1d (gpp                , a_gpp                 )
            call acc1d (downreg            , a_downreg             )
            call acc1d (ar                 , a_ar                  )
            call acc1d (cwdprod            , a_cwdprod             )
            call acc1d (cwddecomp          , a_cwddecomp           )
            call acc1d (decomp_hr          , a_hr                  )
            call acc1d (fpg                , a_fpg                 )
            call acc1d (fpi                , a_fpi                 )
            call acc1d (gpp_enftemp        , a_gpp_enftemp         )
            call acc1d (gpp_enfboreal      , a_gpp_enfboreal       )
            call acc1d (gpp_dnfboreal      , a_gpp_dnfboreal       )
            call acc1d (gpp_ebftrop        , a_gpp_ebftrop         )
            call acc1d (gpp_ebftemp        , a_gpp_ebftemp         )
            call acc1d (gpp_dbftrop        , a_gpp_dbftrop         )
            call acc1d (gpp_dbftemp        , a_gpp_dbftemp         )
            call acc1d (gpp_dbfboreal      , a_gpp_dbfboreal       )
            call acc1d (gpp_ebstemp        , a_gpp_ebstemp         )
            call acc1d (gpp_dbstemp        , a_gpp_dbstemp         )
            call acc1d (gpp_dbsboreal      , a_gpp_dbsboreal       )
            call acc1d (gpp_c3arcgrass     , a_gpp_c3arcgrass      )
            call acc1d (gpp_c3grass        , a_gpp_c3grass         )
            call acc1d (gpp_c4grass        , a_gpp_c4grass         )
            call acc1d (leafc_enftemp      , a_leafc_enftemp       )
            call acc1d (leafc_enfboreal    , a_leafc_enfboreal     )
            call acc1d (leafc_dnfboreal    , a_leafc_dnfboreal     )
            call acc1d (leafc_ebftrop      , a_leafc_ebftrop       )
            call acc1d (leafc_ebftemp      , a_leafc_ebftemp       )
            call acc1d (leafc_dbftrop      , a_leafc_dbftrop       )
            call acc1d (leafc_dbftemp      , a_leafc_dbftemp       )
            call acc1d (leafc_dbfboreal    , a_leafc_dbfboreal     )
            call acc1d (leafc_ebstemp      , a_leafc_ebstemp       )
            call acc1d (leafc_dbstemp      , a_leafc_dbstemp       )
            call acc1d (leafc_dbsboreal    , a_leafc_dbsboreal     )
            call acc1d (leafc_c3arcgrass   , a_leafc_c3arcgrass    )
            call acc1d (leafc_c3grass      , a_leafc_c3grass       )
            call acc1d (leafc_c4grass      , a_leafc_c4grass       )
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
            call acc1d (assim_RuBP_sun_enftemp        , a_assim_RuBP_sun_enftemp        )
            call acc1d (assim_RuBP_sun_enfboreal      , a_assim_RuBP_sun_enfboreal      )
            call acc1d (assim_RuBP_sun_dnfboreal      , a_assim_RuBP_sun_dnfboreal      )
            call acc1d (assim_RuBP_sun_ebftrop        , a_assim_RuBP_sun_ebftrop        )
            call acc1d (assim_RuBP_sun_ebftemp        , a_assim_RuBP_sun_ebftemp        )
            call acc1d (assim_RuBP_sun_dbftrop        , a_assim_RuBP_sun_dbftrop        )
            call acc1d (assim_RuBP_sun_dbftemp        , a_assim_RuBP_sun_dbftemp        )
            call acc1d (assim_RuBP_sun_dbfboreal      , a_assim_RuBP_sun_dbfboreal      )
            call acc1d (assim_RuBP_sun_ebstemp        , a_assim_RuBP_sun_ebstemp        )
            call acc1d (assim_RuBP_sun_dbstemp        ,  a_assim_RuBP_sun_dbstemp        )
            call acc1d (assim_RuBP_sun_dbsboreal      ,  a_assim_RuBP_sun_dbsboreal      )
            call acc1d (assim_RuBP_sun_c3arcgrass     ,  a_assim_RuBP_sun_c3arcgrass     )
            call acc1d (assim_RuBP_sun_c3grass        ,  a_assim_RuBP_sun_c3grass        )
            call acc1d (assim_RuBP_sun_c4grass        ,  a_assim_RuBP_sun_c4grass        )
            call acc1d (assim_RuBP_sha_enftemp        , a_assim_RuBP_sha_enftemp        )
            call acc1d (assim_RuBP_sha_enfboreal      , a_assim_RuBP_sha_enfboreal      )
            call acc1d (assim_RuBP_sha_dnfboreal      , a_assim_RuBP_sha_dnfboreal      )
            call acc1d (assim_RuBP_sha_ebftrop        , a_assim_RuBP_sha_ebftrop        )
            call acc1d (assim_RuBP_sha_ebftemp        , a_assim_RuBP_sha_ebftemp        )
            call acc1d (assim_RuBP_sha_dbftrop        , a_assim_RuBP_sha_dbftrop        )
            call acc1d (assim_RuBP_sha_dbftemp        , a_assim_RuBP_sha_dbftemp        )
            call acc1d (assim_RuBP_sha_dbfboreal      , a_assim_RuBP_sha_dbfboreal      )
            call acc1d (assim_RuBP_sha_ebstemp        , a_assim_RuBP_sha_ebstemp        )
            call acc1d (assim_RuBP_sha_dbstemp        ,  a_assim_RuBP_sha_dbstemp        )
            call acc1d (assim_RuBP_sha_dbsboreal      ,  a_assim_RuBP_sha_dbsboreal      )
            call acc1d (assim_RuBP_sha_c3arcgrass     ,  a_assim_RuBP_sha_c3arcgrass     )
            call acc1d (assim_RuBP_sha_c3grass        ,  a_assim_RuBP_sha_c3grass        )
            call acc1d (assim_RuBP_sha_c4grass        ,  a_assim_RuBP_sha_c4grass        )
            call acc1d (assim_Rubisco_sun_enftemp        , a_assim_Rubisco_sun_enftemp        )
            call acc1d (assim_Rubisco_sun_enfboreal      , a_assim_Rubisco_sun_enfboreal      )
            call acc1d (assim_Rubisco_sun_dnfboreal      , a_assim_Rubisco_sun_dnfboreal      )
            call acc1d (assim_Rubisco_sun_ebftrop        , a_assim_Rubisco_sun_ebftrop        )
            call acc1d (assim_Rubisco_sun_ebftemp        , a_assim_Rubisco_sun_ebftemp        )
            call acc1d (assim_Rubisco_sun_dbftrop        , a_assim_Rubisco_sun_dbftrop        )
            call acc1d (assim_Rubisco_sun_dbftemp        , a_assim_Rubisco_sun_dbftemp        )
            call acc1d (assim_Rubisco_sun_dbfboreal      , a_assim_Rubisco_sun_dbfboreal      )
            call acc1d (assim_Rubisco_sun_ebstemp        , a_assim_Rubisco_sun_ebstemp        )
            call acc1d (assim_Rubisco_sun_dbstemp        ,  a_assim_Rubisco_sun_dbstemp        )
            call acc1d (assim_Rubisco_sun_dbsboreal      ,  a_assim_Rubisco_sun_dbsboreal      )
            call acc1d (assim_Rubisco_sun_c3arcgrass     ,  a_assim_Rubisco_sun_c3arcgrass     )
            call acc1d (assim_Rubisco_sun_c3grass        ,  a_assim_Rubisco_sun_c3grass        )
            call acc1d (assim_Rubisco_sun_c4grass        ,  a_assim_Rubisco_sun_c4grass        )
            call acc1d (assim_Rubisco_sha_enftemp        , a_assim_Rubisco_sha_enftemp        )
            call acc1d (assim_Rubisco_sha_enfboreal      , a_assim_Rubisco_sha_enfboreal      )
            call acc1d (assim_Rubisco_sha_dnfboreal      , a_assim_Rubisco_sha_dnfboreal      )
            call acc1d (assim_Rubisco_sha_ebftrop        , a_assim_Rubisco_sha_ebftrop        )
            call acc1d (assim_Rubisco_sha_ebftemp        , a_assim_Rubisco_sha_ebftemp        )
            call acc1d (assim_Rubisco_sha_dbftrop        , a_assim_Rubisco_sha_dbftrop        )
            call acc1d (assim_Rubisco_sha_dbftemp        , a_assim_Rubisco_sha_dbftemp        )
            call acc1d (assim_Rubisco_sha_dbfboreal      , a_assim_Rubisco_sha_dbfboreal      )
            call acc1d (assim_Rubisco_sha_ebstemp        , a_assim_Rubisco_sha_ebstemp        )
            call acc1d (assim_Rubisco_sha_dbstemp        ,  a_assim_Rubisco_sha_dbstemp        )
            call acc1d (assim_Rubisco_sha_dbsboreal      ,  a_assim_Rubisco_sha_dbsboreal      )
            call acc1d (assim_Rubisco_sha_c3arcgrass     ,  a_assim_Rubisco_sha_c3arcgrass     )
            call acc1d (assim_Rubisco_sha_c3grass        ,  a_assim_Rubisco_sha_c3grass        )
            call acc1d (assim_Rubisco_sha_c4grass        ,  a_assim_Rubisco_sha_c4grass        )
            call acc1d (assimsun_enftemp        , a_assimsun_enftemp        )
            call acc1d (assimsun_enfboreal      , a_assimsun_enfboreal      )
            call acc1d (assimsun_dnfboreal      , a_assimsun_dnfboreal      )
            call acc1d (assimsun_ebftrop        , a_assimsun_ebftrop        )
            call acc1d (assimsun_ebftemp        , a_assimsun_ebftemp        )
            call acc1d (assimsun_dbftrop        , a_assimsun_dbftrop        )
            call acc1d (assimsun_dbftemp        , a_assimsun_dbftemp        )
            call acc1d (assimsun_dbfboreal      , a_assimsun_dbfboreal      )
            call acc1d (assimsun_ebstemp        , a_assimsun_ebstemp        )
            call acc1d (assimsun_dbstemp        ,  a_assimsun_dbstemp        )
            call acc1d (assimsun_dbsboreal      ,  a_assimsun_dbsboreal      )
            call acc1d (assimsun_c3arcgrass     ,  a_assimsun_c3arcgrass     )
            call acc1d (assimsun_c3grass        ,  a_assimsun_c3grass        )
            call acc1d (assimsun_c4grass        ,  a_assimsun_c4grass        )
            call acc1d (assimsha_enftemp        , a_assimsha_enftemp        )
            call acc1d (assimsha_enfboreal      , a_assimsha_enfboreal      )
            call acc1d (assimsha_dnfboreal      , a_assimsha_dnfboreal      )
            call acc1d (assimsha_ebftrop        , a_assimsha_ebftrop        )
            call acc1d (assimsha_ebftemp        , a_assimsha_ebftemp        )
            call acc1d (assimsha_dbftrop        , a_assimsha_dbftrop        )
            call acc1d (assimsha_dbftemp        , a_assimsha_dbftemp        )
            call acc1d (assimsha_dbfboreal      , a_assimsha_dbfboreal      )
            call acc1d (assimsha_ebstemp        , a_assimsha_ebstemp        )
            call acc1d (assimsha_dbstemp        ,  a_assimsha_dbstemp        )
            call acc1d (assimsha_dbsboreal      ,  a_assimsha_dbsboreal      )
            call acc1d (assimsha_c3arcgrass     ,  a_assimsha_c3arcgrass     )
            call acc1d (assimsha_c3grass        ,  a_assimsha_c3grass        )
            call acc1d (assimsha_c4grass        ,  a_assimsha_c4grass        )
            call acc1d (etrsun_enftemp        , a_etrsun_enftemp        )
            call acc1d (etrsun_enfboreal      , a_etrsun_enfboreal      )
            call acc1d (etrsun_dnfboreal      , a_etrsun_dnfboreal      )
            call acc1d (etrsun_ebftrop        , a_etrsun_ebftrop        )
            call acc1d (etrsun_ebftemp        , a_etrsun_ebftemp        )
            call acc1d (etrsun_dbftrop        , a_etrsun_dbftrop        )
            call acc1d (etrsun_dbftemp        , a_etrsun_dbftemp        )
            call acc1d (etrsun_dbfboreal      , a_etrsun_dbfboreal      )
            call acc1d (etrsun_ebstemp        , a_etrsun_ebstemp        )
            call acc1d (etrsun_dbstemp        ,  a_etrsun_dbstemp        )
            call acc1d (etrsun_dbsboreal      ,  a_etrsun_dbsboreal      )
            call acc1d (etrsun_c3arcgrass     ,  a_etrsun_c3arcgrass     )
            call acc1d (etrsun_c3grass        ,  a_etrsun_c3grass        )
            call acc1d (etrsun_c4grass        ,  a_etrsun_c4grass        )
            call acc1d (etrsha_enftemp        , a_etrsha_enftemp        )
            call acc1d (etrsha_enfboreal      , a_etrsha_enfboreal      )
            call acc1d (etrsha_dnfboreal      , a_etrsha_dnfboreal      )
            call acc1d (etrsha_ebftrop        , a_etrsha_ebftrop        )
            call acc1d (etrsha_ebftemp        , a_etrsha_ebftemp        )
            call acc1d (etrsha_dbftrop        , a_etrsha_dbftrop        )
            call acc1d (etrsha_dbftemp        , a_etrsha_dbftemp        )
            call acc1d (etrsha_dbfboreal      , a_etrsha_dbfboreal      )
            call acc1d (etrsha_ebstemp        , a_etrsha_ebstemp        )
            call acc1d (etrsha_dbstemp        ,  a_etrsha_dbstemp        )
            call acc1d (etrsha_dbsboreal      ,  a_etrsha_dbsboreal      )
            call acc1d (etrsha_c3arcgrass     ,  a_etrsha_c3arcgrass     )
            call acc1d (etrsha_c3grass        ,  a_etrsha_c3grass        )
            call acc1d (etrsha_c4grass        ,  a_etrsha_c4grass        )
            call acc1d (cisun_enftemp        , a_cisun_enftemp        )
            call acc1d (cisun_enfboreal      , a_cisun_enfboreal      )
            call acc1d (cisun_dnfboreal      , a_cisun_dnfboreal      )
            call acc1d (cisun_ebftrop        , a_cisun_ebftrop        )
            call acc1d (cisun_ebftemp        , a_cisun_ebftemp        )
            call acc1d (cisun_dbftrop        , a_cisun_dbftrop        )
            call acc1d (cisun_dbftemp        , a_cisun_dbftemp        )
            call acc1d (cisun_dbfboreal      , a_cisun_dbfboreal      )
            call acc1d (cisun_ebstemp        , a_cisun_ebstemp        )
            call acc1d (cisun_dbstemp        ,  a_cisun_dbstemp        )
            call acc1d (cisun_dbsboreal      ,  a_cisun_dbsboreal      )
            call acc1d (cisun_c3arcgrass     ,  a_cisun_c3arcgrass     )
            call acc1d (cisun_c3grass        ,  a_cisun_c3grass        )
            call acc1d (cisun_c4grass        ,  a_cisun_c4grass        )
            call acc1d (cisha_enftemp        , a_cisha_enftemp        )
            call acc1d (cisha_enfboreal      , a_cisha_enfboreal      )
            call acc1d (cisha_dnfboreal      , a_cisha_dnfboreal      )
            call acc1d (cisha_ebftrop        , a_cisha_ebftrop        )
            call acc1d (cisha_ebftemp        , a_cisha_ebftemp        )
            call acc1d (cisha_dbftrop        , a_cisha_dbftrop        )
            call acc1d (cisha_dbftemp        , a_cisha_dbftemp        )
            call acc1d (cisha_dbfboreal      , a_cisha_dbfboreal      )
            call acc1d (cisha_ebstemp        , a_cisha_ebstemp        )
            call acc1d (cisha_dbstemp        ,  a_cisha_dbstemp        )
            call acc1d (cisha_dbsboreal      ,  a_cisha_dbsboreal      )
            call acc1d (cisha_c3arcgrass     ,  a_cisha_c3arcgrass     )
            call acc1d (cisha_c3grass        ,  a_cisha_c3grass        )
            call acc1d (cisha_c4grass        ,  a_cisha_c4grass        )
            call acc1d (essun_enftemp        , a_essun_enftemp        )
            call acc1d (essun_enfboreal      , a_essun_enfboreal      )
            call acc1d (essun_dnfboreal      , a_essun_dnfboreal      )
            call acc1d (essun_ebftrop        , a_essun_ebftrop        )
            call acc1d (essun_ebftemp        , a_essun_ebftemp        )
            call acc1d (essun_dbftrop        , a_essun_dbftrop        )
            call acc1d (essun_dbftemp        , a_essun_dbftemp        )
            call acc1d (essun_dbfboreal      , a_essun_dbfboreal      )
            call acc1d (essun_ebstemp        , a_essun_ebstemp        )
            call acc1d (essun_dbstemp        ,  a_essun_dbstemp        )
            call acc1d (essun_dbsboreal      ,  a_essun_dbsboreal      )
            call acc1d (essun_c3arcgrass     ,  a_essun_c3arcgrass     )
            call acc1d (essun_c3grass        ,  a_essun_c3grass        )
            call acc1d (essun_c4grass        ,  a_essun_c4grass        )
            call acc1d (essha_enftemp        , a_essha_enftemp        )
            call acc1d (essha_enfboreal      , a_essha_enfboreal      )
            call acc1d (essha_dnfboreal      , a_essha_dnfboreal      )
            call acc1d (essha_ebftrop        , a_essha_ebftrop        )
            call acc1d (essha_ebftemp        , a_essha_ebftemp        )
            call acc1d (essha_dbftrop        , a_essha_dbftrop        )
            call acc1d (essha_dbftemp        , a_essha_dbftemp        )
            call acc1d (essha_dbfboreal      , a_essha_dbfboreal      )
            call acc1d (essha_ebstemp        , a_essha_ebstemp        )
            call acc1d (essha_dbstemp        ,  a_essha_dbstemp        )
            call acc1d (essha_dbsboreal      ,  a_essha_dbsboreal      )
            call acc1d (essha_c3arcgrass     ,  a_essha_c3arcgrass     )
            call acc1d (essha_c3grass        ,  a_essha_c3grass        )
            call acc1d (essha_c4grass        ,  a_essha_c4grass        )
            call acc1d (gssun_enftemp        , a_gssun_enftemp        )
            call acc1d (gssun_enfboreal      , a_gssun_enfboreal      )
            call acc1d (gssun_dnfboreal      , a_gssun_dnfboreal      )
            call acc1d (gssun_ebftrop        , a_gssun_ebftrop        )
            call acc1d (gssun_ebftemp        , a_gssun_ebftemp        )
            call acc1d (gssun_dbftrop        , a_gssun_dbftrop        )
            call acc1d (gssun_dbftemp        , a_gssun_dbftemp        )
            call acc1d (gssun_dbfboreal      , a_gssun_dbfboreal      )
            call acc1d (gssun_ebstemp        , a_gssun_ebstemp        )
            call acc1d (gssun_dbstemp        ,  a_gssun_dbstemp        )
            call acc1d (gssun_dbsboreal      ,  a_gssun_dbsboreal      )
            call acc1d (gssun_c3arcgrass     ,  a_gssun_c3arcgrass     )
            call acc1d (gssun_c3grass        ,  a_gssun_c3grass        )
            call acc1d (gssun_c4grass        ,  a_gssun_c4grass        )
            call acc1d (gssha_enftemp        , a_gssha_enftemp        )
            call acc1d (gssha_enfboreal      , a_gssha_enfboreal      )
            call acc1d (gssha_dnfboreal      , a_gssha_dnfboreal      )
            call acc1d (gssha_ebftrop        , a_gssha_ebftrop        )
            call acc1d (gssha_ebftemp        , a_gssha_ebftemp        )
            call acc1d (gssha_dbftrop        , a_gssha_dbftrop        )
            call acc1d (gssha_dbftemp        , a_gssha_dbftemp        )
            call acc1d (gssha_dbfboreal      , a_gssha_dbfboreal      )
            call acc1d (gssha_ebstemp        , a_gssha_ebstemp        )
            call acc1d (gssha_dbstemp        ,  a_gssha_dbstemp        )
            call acc1d (gssha_dbsboreal      ,  a_gssha_dbsboreal      )
            call acc1d (gssha_c3arcgrass     ,  a_gssha_c3arcgrass     )
            call acc1d (gssha_c3grass        ,  a_gssha_c3grass        )
            call acc1d (gssha_c4grass        ,  a_gssha_c4grass        )
            call acc1d (gammasun_enftemp        , a_gammasun_enftemp        )
            call acc1d (gammasun_enfboreal      , a_gammasun_enfboreal      )
            call acc1d (gammasun_dnfboreal      , a_gammasun_dnfboreal      )
            call acc1d (gammasun_ebftrop        , a_gammasun_ebftrop        )
            call acc1d (gammasun_ebftemp        , a_gammasun_ebftemp        )
            call acc1d (gammasun_dbftrop        , a_gammasun_dbftrop        )
            call acc1d (gammasun_dbftemp        , a_gammasun_dbftemp        )
            call acc1d (gammasun_dbfboreal      , a_gammasun_dbfboreal      )
            call acc1d (gammasun_ebstemp        , a_gammasun_ebstemp        )
            call acc1d (gammasun_dbstemp        ,  a_gammasun_dbstemp        )
            call acc1d (gammasun_dbsboreal      ,  a_gammasun_dbsboreal      )
            call acc1d (gammasun_c3arcgrass     ,  a_gammasun_c3arcgrass     )
            call acc1d (gammasun_c3grass        ,  a_gammasun_c3grass        )
            call acc1d (gammasun_c4grass        ,  a_gammasun_c4grass        )
            call acc1d (gammasha_enftemp        , a_gammasha_enftemp        )
            call acc1d (gammasha_enfboreal      , a_gammasha_enfboreal      )
            call acc1d (gammasha_dnfboreal      , a_gammasha_dnfboreal      )
            call acc1d (gammasha_ebftrop        , a_gammasha_ebftrop        )
            call acc1d (gammasha_ebftemp        , a_gammasha_ebftemp        )
            call acc1d (gammasha_dbftrop        , a_gammasha_dbftrop        )
            call acc1d (gammasha_dbftemp        , a_gammasha_dbftemp        )
            call acc1d (gammasha_dbfboreal      , a_gammasha_dbfboreal      )
            call acc1d (gammasha_ebstemp        , a_gammasha_ebstemp        )
            call acc1d (gammasha_dbstemp        ,  a_gammasha_dbstemp        )
            call acc1d (gammasha_dbsboreal      ,  a_gammasha_dbsboreal      )
            call acc1d (gammasha_c3arcgrass     ,  a_gammasha_c3arcgrass     )
            call acc1d (gammasha_c3grass        ,  a_gammasha_c3grass        )
            call acc1d (gammasha_c4grass        ,  a_gammasha_c4grass        )
            call acc1d (RuBPlimitfrac_sun_enftemp        , a_RuBPlimitfrac_sun_enftemp        )
            call acc1d (RuBPlimitfrac_sun_enfboreal      , a_RuBPlimitfrac_sun_enfboreal      )
            call acc1d (RuBPlimitfrac_sun_dnfboreal      , a_RuBPlimitfrac_sun_dnfboreal      )
            call acc1d (RuBPlimitfrac_sun_ebftrop        , a_RuBPlimitfrac_sun_ebftrop        )
            call acc1d (RuBPlimitfrac_sun_ebftemp        , a_RuBPlimitfrac_sun_ebftemp        )
            call acc1d (RuBPlimitfrac_sun_dbftrop        , a_RuBPlimitfrac_sun_dbftrop        )
            call acc1d (RuBPlimitfrac_sun_dbftemp        , a_RuBPlimitfrac_sun_dbftemp        )
            call acc1d (RuBPlimitfrac_sun_dbfboreal      , a_RuBPlimitfrac_sun_dbfboreal      )
            call acc1d (RuBPlimitfrac_sun_ebstemp        , a_RuBPlimitfrac_sun_ebstemp        )
            call acc1d (RuBPlimitfrac_sun_dbstemp        , a_RuBPlimitfrac_sun_dbstemp        )
            call acc1d (RuBPlimitfrac_sun_dbsboreal      , a_RuBPlimitfrac_sun_dbsboreal      )
            call acc1d (RuBPlimitfrac_sun_c3arcgrass     , a_RuBPlimitfrac_sun_c3arcgrass     )
            call acc1d (RuBPlimitfrac_sun_c3grass        , a_RuBPlimitfrac_sun_c3grass        )
            call acc1d (RuBPlimitfrac_sun_c4grass        , a_RuBPlimitfrac_sun_c4grass        )
            call acc1d (RuBPlimitfrac_sha_enftemp        , a_RuBPlimitfrac_sha_enftemp        )
            call acc1d (RuBPlimitfrac_sha_enfboreal      , a_RuBPlimitfrac_sha_enfboreal      )
            call acc1d (RuBPlimitfrac_sha_dnfboreal      , a_RuBPlimitfrac_sha_dnfboreal      )
            call acc1d (RuBPlimitfrac_sha_ebftrop        , a_RuBPlimitfrac_sha_ebftrop        )
            call acc1d (RuBPlimitfrac_sha_ebftemp        , a_RuBPlimitfrac_sha_ebftemp        )
            call acc1d (RuBPlimitfrac_sha_dbftrop        , a_RuBPlimitfrac_sha_dbftrop        )
            call acc1d (RuBPlimitfrac_sha_dbftemp        , a_RuBPlimitfrac_sha_dbftemp        )
            call acc1d (RuBPlimitfrac_sha_dbfboreal      , a_RuBPlimitfrac_sha_dbfboreal      )
            call acc1d (RuBPlimitfrac_sha_ebstemp        , a_RuBPlimitfrac_sha_ebstemp        )
            call acc1d (RuBPlimitfrac_sha_dbstemp        , a_RuBPlimitfrac_sha_dbstemp        )
            call acc1d (RuBPlimitfrac_sha_dbsboreal      , a_RuBPlimitfrac_sha_dbsboreal      )
            call acc1d (RuBPlimitfrac_sha_c3arcgrass     , a_RuBPlimitfrac_sha_c3arcgrass     )
            call acc1d (RuBPlimitfrac_sha_c3grass        , a_RuBPlimitfrac_sha_c3grass        )
            call acc1d (RuBPlimitfrac_sha_c4grass        , a_RuBPlimitfrac_sha_c4grass        )
            call acc1d (Rubiscolimitfrac_sun_enftemp     , a_Rubiscolimitfrac_sun_enftemp     )
            call acc1d (Rubiscolimitfrac_sun_enfboreal   , a_Rubiscolimitfrac_sun_enfboreal   )
            call acc1d (Rubiscolimitfrac_sun_dnfboreal   , a_Rubiscolimitfrac_sun_dnfboreal   )
            call acc1d (Rubiscolimitfrac_sun_ebftrop     , a_Rubiscolimitfrac_sun_ebftrop     )
            call acc1d (Rubiscolimitfrac_sun_ebftemp     , a_Rubiscolimitfrac_sun_ebftemp     )
            call acc1d (Rubiscolimitfrac_sun_dbftrop     , a_Rubiscolimitfrac_sun_dbftrop     )
            call acc1d (Rubiscolimitfrac_sun_dbftemp     , a_Rubiscolimitfrac_sun_dbftemp     )
            call acc1d (Rubiscolimitfrac_sun_dbfboreal   , a_Rubiscolimitfrac_sun_dbfboreal   )
            call acc1d (Rubiscolimitfrac_sun_ebstemp     , a_Rubiscolimitfrac_sun_ebstemp     )
            call acc1d (Rubiscolimitfrac_sun_dbstemp     , a_Rubiscolimitfrac_sun_dbstemp     )
            call acc1d (Rubiscolimitfrac_sun_dbsboreal   , a_Rubiscolimitfrac_sun_dbsboreal   )
            call acc1d (Rubiscolimitfrac_sun_c3arcgrass  , a_Rubiscolimitfrac_sun_c3arcgrass  )
            call acc1d (Rubiscolimitfrac_sun_c3grass     , a_Rubiscolimitfrac_sun_c3grass     )
            call acc1d (Rubiscolimitfrac_sun_c4grass     , a_Rubiscolimitfrac_sun_c4grass     )
            call acc1d (Rubiscolimitfrac_sha_enftemp     , a_Rubiscolimitfrac_sha_enftemp     )
            call acc1d (Rubiscolimitfrac_sha_enfboreal   , a_Rubiscolimitfrac_sha_enfboreal   )
            call acc1d (Rubiscolimitfrac_sha_dnfboreal   , a_Rubiscolimitfrac_sha_dnfboreal   )
            call acc1d (Rubiscolimitfrac_sha_ebftrop     , a_Rubiscolimitfrac_sha_ebftrop     )
            call acc1d (Rubiscolimitfrac_sha_ebftemp     , a_Rubiscolimitfrac_sha_ebftemp     )
            call acc1d (Rubiscolimitfrac_sha_dbftrop     , a_Rubiscolimitfrac_sha_dbftrop     )
            call acc1d (Rubiscolimitfrac_sha_dbftemp     , a_Rubiscolimitfrac_sha_dbftemp     )
            call acc1d (Rubiscolimitfrac_sha_dbfboreal   , a_Rubiscolimitfrac_sha_dbfboreal   )
            call acc1d (Rubiscolimitfrac_sha_ebstemp     , a_Rubiscolimitfrac_sha_ebstemp     )
            call acc1d (Rubiscolimitfrac_sha_dbstemp     , a_Rubiscolimitfrac_sha_dbstemp     )
            call acc1d (Rubiscolimitfrac_sha_dbsboreal   , a_Rubiscolimitfrac_sha_dbsboreal   )
            call acc1d (Rubiscolimitfrac_sha_c3arcgrass  , a_Rubiscolimitfrac_sha_c3arcgrass  )
            call acc1d (Rubiscolimitfrac_sha_c3grass     , a_Rubiscolimitfrac_sha_c3grass     )
            call acc1d (Rubiscolimitfrac_sha_c4grass     , a_Rubiscolimitfrac_sha_c4grass     )
            call acc1d (Sinklimitfrac_sun_enftemp        , a_Sinklimitfrac_sun_enftemp        )
            call acc1d (Sinklimitfrac_sun_enfboreal      , a_Sinklimitfrac_sun_enfboreal      )
            call acc1d (Sinklimitfrac_sun_dnfboreal      , a_Sinklimitfrac_sun_dnfboreal      )
            call acc1d (Sinklimitfrac_sun_ebftrop        , a_Sinklimitfrac_sun_ebftrop        )
            call acc1d (Sinklimitfrac_sun_ebftemp        , a_Sinklimitfrac_sun_ebftemp        )
            call acc1d (Sinklimitfrac_sun_dbftrop        , a_Sinklimitfrac_sun_dbftrop        )
            call acc1d (Sinklimitfrac_sun_dbftemp        , a_Sinklimitfrac_sun_dbftemp        )
            call acc1d (Sinklimitfrac_sun_dbfboreal      , a_Sinklimitfrac_sun_dbfboreal      )
            call acc1d (Sinklimitfrac_sun_ebstemp        , a_Sinklimitfrac_sun_ebstemp        )
            call acc1d (Sinklimitfrac_sun_dbstemp        , a_Sinklimitfrac_sun_dbstemp        )
            call acc1d (Sinklimitfrac_sun_dbsboreal      , a_Sinklimitfrac_sun_dbsboreal      )
            call acc1d (Sinklimitfrac_sun_c3arcgrass     , a_Sinklimitfrac_sun_c3arcgrass     )
            call acc1d (Sinklimitfrac_sun_c3grass        , a_Sinklimitfrac_sun_c3grass        )
            call acc1d (Sinklimitfrac_sun_c4grass        , a_Sinklimitfrac_sun_c4grass        )
            call acc1d (Sinklimitfrac_sha_enftemp        , a_Sinklimitfrac_sha_enftemp        )
            call acc1d (Sinklimitfrac_sha_enfboreal      , a_Sinklimitfrac_sha_enfboreal      )
            call acc1d (Sinklimitfrac_sha_dnfboreal      , a_Sinklimitfrac_sha_dnfboreal      )
            call acc1d (Sinklimitfrac_sha_ebftrop        , a_Sinklimitfrac_sha_ebftrop        )
            call acc1d (Sinklimitfrac_sha_ebftemp        , a_Sinklimitfrac_sha_ebftemp        )
            call acc1d (Sinklimitfrac_sha_dbftrop        , a_Sinklimitfrac_sha_dbftrop        )
            call acc1d (Sinklimitfrac_sha_dbftemp        , a_Sinklimitfrac_sha_dbftemp        )
            call acc1d (Sinklimitfrac_sha_dbfboreal      , a_Sinklimitfrac_sha_dbfboreal      )
            call acc1d (Sinklimitfrac_sha_ebstemp        , a_Sinklimitfrac_sha_ebstemp        )
            call acc1d (Sinklimitfrac_sha_dbstemp        , a_Sinklimitfrac_sha_dbstemp        )
            call acc1d (Sinklimitfrac_sha_dbsboreal      , a_Sinklimitfrac_sha_dbsboreal      )
            call acc1d (Sinklimitfrac_sha_c3arcgrass     , a_Sinklimitfrac_sha_c3arcgrass     )
            call acc1d (Sinklimitfrac_sha_c3grass        , a_Sinklimitfrac_sha_c3grass        )
            call acc1d (Sinklimitfrac_sha_c4grass        , a_Sinklimitfrac_sha_c4grass        )
            call acc1d (rstfacsun_enftemp        , a_rstfacsun_enftemp        )
            call acc1d (rstfacsun_enfboreal      , a_rstfacsun_enfboreal      )
            call acc1d (rstfacsun_dnfboreal      , a_rstfacsun_dnfboreal      )
            call acc1d (rstfacsun_ebftrop        , a_rstfacsun_ebftrop        )
            call acc1d (rstfacsun_ebftemp        , a_rstfacsun_ebftemp        )
            call acc1d (rstfacsun_dbftrop        , a_rstfacsun_dbftrop        )
            call acc1d (rstfacsun_dbftemp        , a_rstfacsun_dbftemp        )
            call acc1d (rstfacsun_dbfboreal      , a_rstfacsun_dbfboreal      )
            call acc1d (rstfacsun_ebstemp        , a_rstfacsun_ebstemp        )
            call acc1d (rstfacsun_dbstemp        , a_rstfacsun_dbstemp        )
            call acc1d (rstfacsun_dbsboreal      , a_rstfacsun_dbsboreal      )
            call acc1d (rstfacsun_c3arcgrass     , a_rstfacsun_c3arcgrass     )
            call acc1d (rstfacsun_c3grass        , a_rstfacsun_c3grass        )
            call acc1d (rstfacsun_c4grass        , a_rstfacsun_c4grass        )
            call acc1d (rstfacsha_enftemp        , a_rstfacsha_enftemp        )
            call acc1d (rstfacsha_enfboreal      , a_rstfacsha_enfboreal      )
            call acc1d (rstfacsha_dnfboreal      , a_rstfacsha_dnfboreal      )
            call acc1d (rstfacsha_ebftrop        , a_rstfacsha_ebftrop        )
            call acc1d (rstfacsha_ebftemp        , a_rstfacsha_ebftemp        )
            call acc1d (rstfacsha_dbftrop        , a_rstfacsha_dbftrop        )
            call acc1d (rstfacsha_dbftemp        , a_rstfacsha_dbftemp        )
            call acc1d (rstfacsha_dbfboreal      , a_rstfacsha_dbfboreal      )
            call acc1d (rstfacsha_ebstemp        , a_rstfacsha_ebstemp        )
            call acc1d (rstfacsha_dbstemp        , a_rstfacsha_dbstemp        )
            call acc1d (rstfacsha_dbsboreal      , a_rstfacsha_dbsboreal      )
            call acc1d (rstfacsha_c3arcgrass     , a_rstfacsha_c3arcgrass     )
            call acc1d (rstfacsha_c3grass        , a_rstfacsha_c3grass        )
            call acc1d (rstfacsha_c4grass        , a_rstfacsha_c4grass        )
            call acc1d (lambdasun_enftemp        , a_lambdasun_enftemp        )
            call acc1d (lambdasun_enfboreal      , a_lambdasun_enfboreal      )
            call acc1d (lambdasun_dnfboreal      , a_lambdasun_dnfboreal      )
            call acc1d (lambdasun_ebftrop        , a_lambdasun_ebftrop        )
            call acc1d (lambdasun_ebftemp        , a_lambdasun_ebftemp        )
            call acc1d (lambdasun_dbftrop        , a_lambdasun_dbftrop        )
            call acc1d (lambdasun_dbftemp        , a_lambdasun_dbftemp        )
            call acc1d (lambdasun_dbfboreal      , a_lambdasun_dbfboreal      )
            call acc1d (lambdasun_ebstemp        , a_lambdasun_ebstemp        )
            call acc1d (lambdasun_dbstemp        , a_lambdasun_dbstemp        )
            call acc1d (lambdasun_dbsboreal      , a_lambdasun_dbsboreal      )
            call acc1d (lambdasun_c3arcgrass     , a_lambdasun_c3arcgrass     )
            call acc1d (lambdasun_c3grass        , a_lambdasun_c3grass        )
            call acc1d (lambdasun_c4grass        , a_lambdasun_c4grass        )
            call acc1d (lambdasha_enftemp        , a_lambdasha_enftemp        )
            call acc1d (lambdasha_enfboreal      , a_lambdasha_enfboreal      )
            call acc1d (lambdasha_dnfboreal      , a_lambdasha_dnfboreal      )
            call acc1d (lambdasha_ebftrop        , a_lambdasha_ebftrop        )
            call acc1d (lambdasha_ebftemp        , a_lambdasha_ebftemp        )
            call acc1d (lambdasha_dbftrop        , a_lambdasha_dbftrop        )
            call acc1d (lambdasha_dbftemp        , a_lambdasha_dbftemp        )
            call acc1d (lambdasha_dbfboreal      , a_lambdasha_dbfboreal      )
            call acc1d (lambdasha_ebstemp        , a_lambdasha_ebstemp        )
            call acc1d (lambdasha_dbstemp        , a_lambdasha_dbstemp        )
            call acc1d (lambdasha_dbsboreal      , a_lambdasha_dbsboreal      )
            call acc1d (lambdasha_c3arcgrass     , a_lambdasha_c3arcgrass     )
            call acc1d (lambdasha_c3grass        , a_lambdasha_c3grass        )
            call acc1d (lambdasha_c4grass        , a_lambdasha_c4grass        )
            call acc1d (lambda                   , a_lambda                   )
#endif
#endif
#ifdef NITRIF
            call acc2d (to2_decomp_depth_unsat, a_O2_DECOMP_DEPTH_UNSAT)
            call acc2d (tconc_o2_unsat        , a_CONC_O2_UNSAT        )
#endif
#ifdef CROP
            call acc1d (pdcorn             ,   a_pdcorn             )
            call acc1d (pdswheat           ,   a_pdswheat           )
            call acc1d (pdwwheat           ,   a_pdwwheat           )
            call acc1d (pdsoybean          ,   a_pdsoybean          )
            call acc1d (pdcotton           ,   a_pdcotton           )
            call acc1d (pdrice1            ,   a_pdrice1            )
            call acc1d (pdrice2            ,   a_pdrice2            )
            call acc1d (pdsugarcane        ,   a_pdsugarcane        )
            call acc1d (plantdate          ,   a_plantdate          )
            call acc1d (fertnitro_corn     ,   a_fertnitro_corn     )
            call acc1d (fertnitro_swheat   ,   a_fertnitro_swheat   )
            call acc1d (fertnitro_wwheat   ,   a_fertnitro_wwheat   )
            call acc1d (fertnitro_soybean  ,   a_fertnitro_soybean  )
            call acc1d (fertnitro_cotton   ,   a_fertnitro_cotton   )
            call acc1d (fertnitro_rice1    ,   a_fertnitro_rice1    )
            call acc1d (fertnitro_rice2    ,   a_fertnitro_rice2    )
            call acc1d (fertnitro_sugarcane, a_fertnitro_sugarcane  )
            call acc1d (cphase             ,   a_cphase             )
            call acc1d (hui                ,   a_hui                )
            call acc1d (vf                 ,   a_vf                 )
            call acc1d (gddmaturity        ,   a_gddmaturity        )
            call acc1d (gddplant           ,   a_gddplant           )           
            call acc1d (cropprod1c         ,   a_cropprod1c         )
            call acc1d (cropprod1c_loss    ,   a_cropprod1c_loss    )
            call acc1d (cropseedc_deficit  ,   a_cropseedc_deficit  )
            call acc1d (grainc_to_cropprodc,   a_grainc_to_cropprodc)
            call acc1d (grainc_to_seed     ,   a_grainc_to_seed     )
            call acc1d (fert_to_sminn      ,   a_fert_to_sminn      )
#endif
            call acc1d (ndep_to_sminn      ,   a_ndep_to_sminn      )
#ifdef Fire
            call acc1d (abm_lf             ,   a_abm                )
            call acc1d (gdp_lf             ,   a_gdp                )
            call acc1d (peatf_lf           ,   a_peatf              )
            call acc1d (hdm_lf             ,   a_hdm                )
            call acc1d (lnfm               ,   a_lnfm               )
#endif
#endif
#ifdef OzoneStress
            call acc1d (forc_ozone         ,   a_ozone              )
#endif

            call acc2d (t_soisno   , a_t_soisno   )
            call acc2d (wliq_soisno, a_wliq_soisno)
            call acc2d (wice_soisno, a_wice_soisno)

            call acc2d (h2osoi     , a_h2osoi     )
            call acc2d (rootr      , a_rootr      )
            call acc2d (BD_all     , a_BD_all      )
            call acc2d (wfc        , a_wfc         )
            call acc2d (OM_density , a_OM_density  )
#ifdef PLANT_HYDRAULIC_STRESS
            call acc2d (vegwp      , a_vegwp      )
#endif
            call acc2d (t_lake      , a_t_lake      )
            call acc2d (lake_icefrac, a_lake_icefrac)
#ifdef BGC
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_met_lit,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_litr1c_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_cel_lit,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_litr2c_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_lig_lit,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_litr3c_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_soil1,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_soil1c_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_soil2,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_soil2c_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_soil3,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_soil3c_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_cpools_vr(j,i_cwd,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_cwdc_vr     )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_met_lit,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_litr1n_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_cel_lit,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_litr2n_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_lig_lit,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_litr3n_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_soil1,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_soil1n_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_soil2,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_soil2n_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_soil3,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_soil3n_vr   )
            do i = 1, numpatch
               do j = 1, nl_soil
                  decomp_vr_tmp(j,i)  = decomp_npools_vr(j,i_cwd,i)
               end do
            end do
            call acc2d (decomp_vr_tmp, a_cwdn_vr     )
            call acc2d (sminn_vr     , a_sminn_vr    )
#endif
            allocate (r_ustar (numpatch))
            allocate (r_tstar (numpatch))
            allocate (r_qstar (numpatch))
            allocate (r_zol   (numpatch))
            allocate (r_rib   (numpatch))
            allocate (r_fm    (numpatch))
            allocate (r_fh    (numpatch))
            allocate (r_fq    (numpatch))

            allocate (r_us10m (numpatch))
            allocate (r_vs10m (numpatch))
            allocate (r_fm10m (numpatch))

            do i = 1, numpatch

               z0m_av = z0m(i) 
               z0h_av = z0m(i)
               z0q_av = z0m(i)

               displa_av = 2./3.*z0m_av/0.07

               hgt_u = max(forc_hgt_u(i), 5.+displa_av)
               hgt_t = max(forc_hgt_t(i), 5.+displa_av)
               hgt_q = max(forc_hgt_q(i), 5.+displa_av)
               zldis = hgt_u-displa_av

               us = forc_us(i)
               vs = forc_vs(i)
               tm = forc_t (i)
               qm = forc_q (i)
               psrf = forc_psrf(i)
               rhoair = (psrf - 0.378*qm*psrf/(0.622+0.378*qm)) / (rgas*tm)

               r_ustar(i) = sqrt(max(1.e-6,sqrt(taux(i)**2+tauy(i)**2))/rhoair)
               r_tstar(i) = -fsena(i)/(rhoair*r_ustar(i))/cpair
               r_qstar(i) = -fevpa(i)/(rhoair*r_ustar(i))

               thm = tm + 0.0098*hgt_t
               th  = tm*(100000./psrf)**(rgas/cpair)
               thv = th*(1.+0.61*qm)

               r_zol(i) = zldis*vonkar*grav * (r_tstar(i)+0.61*th*r_qstar(i)) &
                  / (r_ustar(i)**2*thv)

               if(r_zol(i) >= 0.)then   !stable
                  r_zol(i) = min(2.,max(r_zol(i),1.e-6))
               else                       !unstable
                  r_zol(i) = max(-100.,min(r_zol(i),-1.e-6))
               endif

               beta = 1.
               zii = 1000.
               thvstar=r_tstar(i)+0.61*th*r_qstar(i)
               ur = sqrt(us*us+vs*vs)
               if(r_zol(i) >= 0.)then
                  um = max(ur,0.1)
               else
                  wc = (-grav*r_ustar(i)*thvstar*zii/thv)**(1./3.)
                  wc2 = beta*beta*(wc*wc)
                  um = max(0.1,sqrt(ur*ur+wc2))
               endif

               obu = zldis/r_zol(i)
               call moninobuk(hgt_u,hgt_t,hgt_q,displa_av,z0m_av,z0h_av,z0q_av,&
                  obu,um,r_ustar(i),fh2m,fq2m,r_fm10m(i),r_fm(i),r_fh(i),r_fq(i))

               ! bug found by chen qiying 2013/07/01 
               r_rib(i) = r_zol(i) /vonkar * r_ustar(i)**2 / (vonkar/r_fh(i)*um**2)
               r_rib(i) = min(5.,r_rib(i))

               r_us10m(i) = us/um * r_ustar(i) /vonkar * r_fm10m(i)
               r_vs10m(i) = vs/um * r_ustar(i) /vonkar * r_fm10m(i)

            end do

            call acc1d (r_ustar, a_ustar)
            call acc1d (r_tstar, a_tstar)
            call acc1d (r_qstar, a_qstar)
            call acc1d (r_zol  , a_zol  )
            call acc1d (r_rib  , a_rib  )
            call acc1d (r_fm   , a_fm   )
            call acc1d (r_fh   , a_fh   )
            call acc1d (r_fq   , a_fq   )

            call acc1d (r_us10m, a_us10m)
            call acc1d (r_vs10m, a_vs10m)
            call acc1d (r_fm10m, a_fm10m)

            deallocate (r_ustar )
            deallocate (r_tstar )
            deallocate (r_qstar )
            deallocate (r_zol   )
            deallocate (r_rib   )
            deallocate (r_fm    )
            deallocate (r_fh    )
            deallocate (r_fq    )

            deallocate (r_us10m )
            deallocate (r_vs10m )
            deallocate (r_fm10m )

            call acc1d (sr     , a_sr     )
            call acc1d (solvd  , a_solvd  )
            call acc1d (solvi  , a_solvi  )
            call acc1d (solnd  , a_solnd  )
            call acc1d (solni  , a_solni  )
            call acc1d (srvd   , a_srvd   )
            call acc1d (srvi   , a_srvi   )
            call acc1d (srnd   , a_srnd   )
            call acc1d (srni   , a_srni   )
            call acc1d (solvdln, a_solvdln)
            call acc1d (solviln, a_solviln)
            call acc1d (solndln, a_solndln)
            call acc1d (solniln, a_solniln)
            call acc1d (srvdln , a_srvdln )
            call acc1d (srviln , a_srviln )
            call acc1d (srndln , a_srndln )
            call acc1d (srniln , a_srniln )

            do i = 1, numpatch
               if (solvdln(i) /= spval) then
                  nac_ln(i) = nac_ln(i) + 1
               end if
            end do

         end if
      end if

   END SUBROUTINE accumulate_fluxes


   !------
   SUBROUTINE acc1d (var, s)

      use precision
      use GlobalVars, only: spval

      IMPLICIT NONE

      real(r8), intent(in)    :: var(:)
      real(r8), intent(inout) :: s  (:)
      ! Local variables
      integer :: i

      do i = lbound(var,1), ubound(var,1)
         if (var(i) /= spval) then
            if (s(i) /= spval) then
               s(i) = s(i) + var(i)
            else
               s(i) = var(i)
            end if
         end if
      end do
      
   END SUBROUTINE acc1d

   !------
   SUBROUTINE acc2d (var, s)

      use precision
      use GlobalVars, only: spval

      IMPLICIT NONE

      real(r8), intent(in)    :: var(:,:)
      real(r8), intent(inout) :: s  (:,:)
      ! Local variables
      integer :: i1, i2

      do i2 = lbound(var,2), ubound(var,2)
         do i1 = lbound(var,1), ubound(var,1)
            if (var(i1,i2) /= spval) then
               if (s(i1,i2) /= spval) then
                  s(i1,i2) = s(i1,i2) + var(i1,i2)
               else
                  s(i1,i2) = var(i1,i2)
               end if
            end if
         end do
      end do

   END SUBROUTINE acc2d

   !------
   SUBROUTINE acc3d (var, s)

      use precision
      use GlobalVars, only: spval

      IMPLICIT NONE

      real(r8), intent(in)    :: var(:,:,:)
      real(r8), intent(inout) :: s  (:,:,:)
      ! Local variables
      integer :: i1, i2, i3

      do i3 = lbound(var,3), ubound(var,3)
         do i2 = lbound(var,2), ubound(var,2)
            do i1 = lbound(var,1), ubound(var,1)
               if (var(i1,i2,i3) /= spval) then
                  if (s(i1,i2,i3) /= spval) then
                     s(i1,i2,i3) = s(i1,i2,i3) + var(i1,i2,i3)
                  else
                     s(i1,i2,i3) = var(i1,i2,i3)
                  end if
               end if
            end do
         end do
      end do

   END SUBROUTINE acc3d

end module MOD_1D_Acc_Fluxes
! ----- EOP ---------
