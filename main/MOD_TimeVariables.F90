#include <define.h>

MODULE MOD_TimeVariables
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

use precision
use timemanager
#ifdef PFT_CLASSIFICATION
USE MOD_PFTimeVars
#endif
#ifdef PC_CLASSIFICATION
USE MOD_PCTimeVars
#endif
#ifdef BGC
USE MOD_BGCTimeVars
#endif
IMPLICIT NONE
SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run
      real(r8), allocatable :: z_sno    (:,:)   ! node depth [m]
      real(r8), allocatable :: dz_sno   (:,:)   ! interface depth [m]
      real(r8), allocatable :: t_soisno (:,:)   ! soil temperature [K]
      real(r8), allocatable :: wliq_soisno(:,:) ! liquid water in layers [kg/m2]
      real(r8), allocatable :: wice_soisno(:,:) ! ice lens in layers [kg/m2]
      real(r8), allocatable :: h2osoi (:,:)     ! volumetric soil water in layers [m3/m3]
      real(r8), allocatable :: smp(:,:)         ! soil matrix potential [mm]
      real(r8), allocatable :: hk (:,:)         ! hydraulic conductivity [mm h2o/s]
      real(r8), allocatable :: rootr(:,:)       ! water exchange between soil and root. Positive: soil->root [?]
#ifdef PLANT_HYDRAULIC_STRESS
      real(r8), allocatable :: vegwp(:,:)       ! vegetation water potential [mm]
      real(r8), allocatable :: gs0sun   (:)     ! working copy of sunlit stomata conductance
      real(r8), allocatable :: gs0sha   (:)     ! working copy of shalit stomata conductance
#endif
#ifdef OzoneStress
      real(r8), allocatable :: o3coefv_sun(:) ! Ozone stress factor for photosynthesis on sunlit leaf
      real(r8), allocatable :: o3coefv_sha(:) ! Ozone stress factor for photosynthesis on shaded leaf
      real(r8), allocatable :: o3coefg_sun(:) ! Ozone stress factor for stomata on sunlit leaf
      real(r8), allocatable :: o3coefg_sha(:) ! Ozone stress factor for stomata on shaded leaf
      real(r8), allocatable :: lai_old    (:) ! lai in last time step
      real(r8), allocatable :: o3uptakesun(:) ! Ozone does, sunlit leaf (mmol O3/m^2)
      real(r8), allocatable :: o3uptakesha(:) ! Ozone does, shaded leaf (mmol O3/m^2)
#endif
      real(r8), allocatable :: rstfacsun(:)     ! factor of soil water stress on sunlit leaf
      real(r8), allocatable :: rstfacsha(:)     ! factor of soil water stress on shaded leaf
      real(r8), allocatable :: gs_sun(:)      ! stomata conductance on sunlit leaf
      real(r8), allocatable :: gs_sha(:)      ! stomata conductance on shaded leaf
      real(r8), allocatable :: t_grnd   (:)     ! ground surface temperature [K]

#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
      real(r8), allocatable :: assim_RuBP_sun_enftemp        (:) !1
      real(r8), allocatable :: assim_RuBP_sun_enfboreal      (:) !2
      real(r8), allocatable :: assim_RuBP_sun_dnfboreal      (:) !3
      real(r8), allocatable :: assim_RuBP_sun_ebftrop        (:) !4
      real(r8), allocatable :: assim_RuBP_sun_ebftemp        (:) !5
      real(r8), allocatable :: assim_RuBP_sun_dbftrop        (:) !6
      real(r8), allocatable :: assim_RuBP_sun_dbftemp        (:) !7
      real(r8), allocatable :: assim_RuBP_sun_dbfboreal      (:) !8
      real(r8), allocatable :: assim_RuBP_sun_ebstemp        (:) !9
      real(r8), allocatable :: assim_RuBP_sun_dbstemp        (:) !10
      real(r8), allocatable :: assim_RuBP_sun_dbsboreal      (:) !11
      real(r8), allocatable :: assim_RuBP_sun_c3arcgrass     (:) !12
      real(r8), allocatable :: assim_RuBP_sun_c3grass        (:) !13
      real(r8), allocatable :: assim_RuBP_sun_c4grass        (:) !14
      real(r8), allocatable :: assim_RuBP_sha_enftemp        (:) !1
      real(r8), allocatable :: assim_RuBP_sha_enfboreal      (:) !2
      real(r8), allocatable :: assim_RuBP_sha_dnfboreal      (:) !3
      real(r8), allocatable :: assim_RuBP_sha_ebftrop        (:) !4
      real(r8), allocatable :: assim_RuBP_sha_ebftemp        (:) !5
      real(r8), allocatable :: assim_RuBP_sha_dbftrop        (:) !6
      real(r8), allocatable :: assim_RuBP_sha_dbftemp        (:) !7
      real(r8), allocatable :: assim_RuBP_sha_dbfboreal      (:) !8
      real(r8), allocatable :: assim_RuBP_sha_ebstemp        (:) !9
      real(r8), allocatable :: assim_RuBP_sha_dbstemp        (:) !10
      real(r8), allocatable :: assim_RuBP_sha_dbsboreal      (:) !11
      real(r8), allocatable :: assim_RuBP_sha_c3arcgrass     (:) !12
      real(r8), allocatable :: assim_RuBP_sha_c3grass        (:) !13
      real(r8), allocatable :: assim_RuBP_sha_c4grass        (:) !14
      real(r8), allocatable :: assim_Rubisco_sun_enftemp        (:) !1
      real(r8), allocatable :: assim_Rubisco_sun_enfboreal      (:) !2
      real(r8), allocatable :: assim_Rubisco_sun_dnfboreal      (:) !3
      real(r8), allocatable :: assim_Rubisco_sun_ebftrop        (:) !4
      real(r8), allocatable :: assim_Rubisco_sun_ebftemp        (:) !5
      real(r8), allocatable :: assim_Rubisco_sun_dbftrop        (:) !6
      real(r8), allocatable :: assim_Rubisco_sun_dbftemp        (:) !7
      real(r8), allocatable :: assim_Rubisco_sun_dbfboreal      (:) !8
      real(r8), allocatable :: assim_Rubisco_sun_ebstemp        (:) !9
      real(r8), allocatable :: assim_Rubisco_sun_dbstemp        (:) !10
      real(r8), allocatable :: assim_Rubisco_sun_dbsboreal      (:) !11
      real(r8), allocatable :: assim_Rubisco_sun_c3arcgrass     (:) !12
      real(r8), allocatable :: assim_Rubisco_sun_c3grass        (:) !13
      real(r8), allocatable :: assim_Rubisco_sun_c4grass        (:) !14
      real(r8), allocatable :: assim_Rubisco_sha_enftemp        (:) !1
      real(r8), allocatable :: assim_Rubisco_sha_enfboreal      (:) !2
      real(r8), allocatable :: assim_Rubisco_sha_dnfboreal      (:) !3
      real(r8), allocatable :: assim_Rubisco_sha_ebftrop        (:) !4
      real(r8), allocatable :: assim_Rubisco_sha_ebftemp        (:) !5
      real(r8), allocatable :: assim_Rubisco_sha_dbftrop        (:) !6
      real(r8), allocatable :: assim_Rubisco_sha_dbftemp        (:) !7
      real(r8), allocatable :: assim_Rubisco_sha_dbfboreal      (:) !8
      real(r8), allocatable :: assim_Rubisco_sha_ebstemp        (:) !9
      real(r8), allocatable :: assim_Rubisco_sha_dbstemp        (:) !10
      real(r8), allocatable :: assim_Rubisco_sha_dbsboreal      (:) !11
      real(r8), allocatable :: assim_Rubisco_sha_c3arcgrass     (:) !12
      real(r8), allocatable :: assim_Rubisco_sha_c3grass        (:) !13
      real(r8), allocatable :: assim_Rubisco_sha_c4grass        (:) !14
      real(r8), allocatable :: assimsun_enftemp        (:) !1
      real(r8), allocatable :: assimsun_enfboreal      (:) !2
      real(r8), allocatable :: assimsun_dnfboreal      (:) !3
      real(r8), allocatable :: assimsun_ebftrop        (:) !4
      real(r8), allocatable :: assimsun_ebftemp        (:) !5
      real(r8), allocatable :: assimsun_dbftrop        (:) !6
      real(r8), allocatable :: assimsun_dbftemp        (:) !7
      real(r8), allocatable :: assimsun_dbfboreal      (:) !8
      real(r8), allocatable :: assimsun_ebstemp        (:) !9
      real(r8), allocatable :: assimsun_dbstemp        (:) !10
      real(r8), allocatable :: assimsun_dbsboreal      (:) !11
      real(r8), allocatable :: assimsun_c3arcgrass     (:) !12
      real(r8), allocatable :: assimsun_c3grass        (:) !13
      real(r8), allocatable :: assimsun_c4grass        (:) !14
      real(r8), allocatable :: assimsha_enftemp        (:) !1
      real(r8), allocatable :: assimsha_enfboreal      (:) !2
      real(r8), allocatable :: assimsha_dnfboreal      (:) !3
      real(r8), allocatable :: assimsha_ebftrop        (:) !4
      real(r8), allocatable :: assimsha_ebftemp        (:) !5
      real(r8), allocatable :: assimsha_dbftrop        (:) !6
      real(r8), allocatable :: assimsha_dbftemp        (:) !7
      real(r8), allocatable :: assimsha_dbfboreal      (:) !8
      real(r8), allocatable :: assimsha_ebstemp        (:) !9
      real(r8), allocatable :: assimsha_dbstemp        (:) !10
      real(r8), allocatable :: assimsha_dbsboreal      (:) !11
      real(r8), allocatable :: assimsha_c3arcgrass     (:) !12
      real(r8), allocatable :: assimsha_c3grass        (:) !13
      real(r8), allocatable :: assimsha_c4grass        (:) !14
      real(r8), allocatable :: etrsun_enftemp        (:) !1
      real(r8), allocatable :: etrsun_enfboreal      (:) !2
      real(r8), allocatable :: etrsun_dnfboreal      (:) !3
      real(r8), allocatable :: etrsun_ebftrop        (:) !4
      real(r8), allocatable :: etrsun_ebftemp        (:) !5
      real(r8), allocatable :: etrsun_dbftrop        (:) !6
      real(r8), allocatable :: etrsun_dbftemp        (:) !7
      real(r8), allocatable :: etrsun_dbfboreal      (:) !8
      real(r8), allocatable :: etrsun_ebstemp        (:) !9
      real(r8), allocatable :: etrsun_dbstemp        (:) !10
      real(r8), allocatable :: etrsun_dbsboreal      (:) !11
      real(r8), allocatable :: etrsun_c3arcgrass     (:) !12
      real(r8), allocatable :: etrsun_c3grass        (:) !13
      real(r8), allocatable :: etrsun_c4grass        (:) !14
      real(r8), allocatable :: etrsha_enftemp        (:) !1
      real(r8), allocatable :: etrsha_enfboreal      (:) !2
      real(r8), allocatable :: etrsha_dnfboreal      (:) !3
      real(r8), allocatable :: etrsha_ebftrop        (:) !4
      real(r8), allocatable :: etrsha_ebftemp        (:) !5
      real(r8), allocatable :: etrsha_dbftrop        (:) !6
      real(r8), allocatable :: etrsha_dbftemp        (:) !7
      real(r8), allocatable :: etrsha_dbfboreal      (:) !8
      real(r8), allocatable :: etrsha_ebstemp        (:) !9
      real(r8), allocatable :: etrsha_dbstemp        (:) !10
      real(r8), allocatable :: etrsha_dbsboreal      (:) !11
      real(r8), allocatable :: etrsha_c3arcgrass     (:) !12
      real(r8), allocatable :: etrsha_c3grass        (:) !13
      real(r8), allocatable :: etrsha_c4grass        (:) !14
      real(r8), allocatable :: cisun_enftemp        (:) !1
      real(r8), allocatable :: cisun_enfboreal      (:) !2
      real(r8), allocatable :: cisun_dnfboreal      (:) !3
      real(r8), allocatable :: cisun_ebftrop        (:) !4
      real(r8), allocatable :: cisun_ebftemp        (:) !5
      real(r8), allocatable :: cisun_dbftrop        (:) !6
      real(r8), allocatable :: cisun_dbftemp        (:) !7
      real(r8), allocatable :: cisun_dbfboreal      (:) !8
      real(r8), allocatable :: cisun_ebstemp        (:) !9
      real(r8), allocatable :: cisun_dbstemp        (:) !10
      real(r8), allocatable :: cisun_dbsboreal      (:) !11
      real(r8), allocatable :: cisun_c3arcgrass     (:) !12
      real(r8), allocatable :: cisun_c3grass        (:) !13
      real(r8), allocatable :: cisun_c4grass        (:) !14
      real(r8), allocatable :: cisha_enftemp        (:) !1
      real(r8), allocatable :: cisha_enfboreal      (:) !2
      real(r8), allocatable :: cisha_dnfboreal      (:) !3
      real(r8), allocatable :: cisha_ebftrop        (:) !4
      real(r8), allocatable :: cisha_ebftemp        (:) !5
      real(r8), allocatable :: cisha_dbftrop        (:) !6
      real(r8), allocatable :: cisha_dbftemp        (:) !7
      real(r8), allocatable :: cisha_dbfboreal      (:) !8
      real(r8), allocatable :: cisha_ebstemp        (:) !9
      real(r8), allocatable :: cisha_dbstemp        (:) !10
      real(r8), allocatable :: cisha_dbsboreal      (:) !11
      real(r8), allocatable :: cisha_c3arcgrass     (:) !12
      real(r8), allocatable :: cisha_c3grass        (:) !13
      real(r8), allocatable :: cisha_c4grass        (:) !14
      real(r8), allocatable :: essun_enftemp        (:) !1
      real(r8), allocatable :: essun_enfboreal      (:) !2
      real(r8), allocatable :: essun_dnfboreal      (:) !3
      real(r8), allocatable :: essun_ebftrop        (:) !4
      real(r8), allocatable :: essun_ebftemp        (:) !5
      real(r8), allocatable :: essun_dbftrop        (:) !6
      real(r8), allocatable :: essun_dbftemp        (:) !7
      real(r8), allocatable :: essun_dbfboreal      (:) !8
      real(r8), allocatable :: essun_ebstemp        (:) !9
      real(r8), allocatable :: essun_dbstemp        (:) !10
      real(r8), allocatable :: essun_dbsboreal      (:) !11
      real(r8), allocatable :: essun_c3arcgrass     (:) !12
      real(r8), allocatable :: essun_c3grass        (:) !13
      real(r8), allocatable :: essun_c4grass        (:) !14
      real(r8), allocatable :: essha_enftemp        (:) !1
      real(r8), allocatable :: essha_enfboreal      (:) !2
      real(r8), allocatable :: essha_dnfboreal      (:) !3
      real(r8), allocatable :: essha_ebftrop        (:) !4
      real(r8), allocatable :: essha_ebftemp        (:) !5
      real(r8), allocatable :: essha_dbftrop        (:) !6
      real(r8), allocatable :: essha_dbftemp        (:) !7
      real(r8), allocatable :: essha_dbfboreal      (:) !8
      real(r8), allocatable :: essha_ebstemp        (:) !9
      real(r8), allocatable :: essha_dbstemp        (:) !10
      real(r8), allocatable :: essha_dbsboreal      (:) !11
      real(r8), allocatable :: essha_c3arcgrass     (:) !12
      real(r8), allocatable :: essha_c3grass        (:) !13
      real(r8), allocatable :: essha_c4grass        (:) !14
      real(r8), allocatable :: gssun_enftemp        (:) !1
      real(r8), allocatable :: gssun_enfboreal      (:) !2
      real(r8), allocatable :: gssun_dnfboreal      (:) !3
      real(r8), allocatable :: gssun_ebftrop        (:) !4
      real(r8), allocatable :: gssun_ebftemp        (:) !5
      real(r8), allocatable :: gssun_dbftrop        (:) !6
      real(r8), allocatable :: gssun_dbftemp        (:) !7
      real(r8), allocatable :: gssun_dbfboreal      (:) !8
      real(r8), allocatable :: gssun_ebstemp        (:) !9
      real(r8), allocatable :: gssun_dbstemp        (:) !10
      real(r8), allocatable :: gssun_dbsboreal      (:) !11
      real(r8), allocatable :: gssun_c3arcgrass     (:) !12
      real(r8), allocatable :: gssun_c3grass        (:) !13
      real(r8), allocatable :: gssun_c4grass        (:) !14
      real(r8), allocatable :: gssha_enftemp        (:) !1
      real(r8), allocatable :: gssha_enfboreal      (:) !2
      real(r8), allocatable :: gssha_dnfboreal      (:) !3
      real(r8), allocatable :: gssha_ebftrop        (:) !4
      real(r8), allocatable :: gssha_ebftemp        (:) !5
      real(r8), allocatable :: gssha_dbftrop        (:) !6
      real(r8), allocatable :: gssha_dbftemp        (:) !7
      real(r8), allocatable :: gssha_dbfboreal      (:) !8
      real(r8), allocatable :: gssha_ebstemp        (:) !9
      real(r8), allocatable :: gssha_dbstemp        (:) !10
      real(r8), allocatable :: gssha_dbsboreal      (:) !11
      real(r8), allocatable :: gssha_c3arcgrass     (:) !12
      real(r8), allocatable :: gssha_c3grass        (:) !13
      real(r8), allocatable :: gssha_c4grass        (:) !14
      real(r8), allocatable :: gammasun_enftemp        (:) !1
      real(r8), allocatable :: gammasun_enfboreal      (:) !2
      real(r8), allocatable :: gammasun_dnfboreal      (:) !3
      real(r8), allocatable :: gammasun_ebftrop        (:) !4
      real(r8), allocatable :: gammasun_ebftemp        (:) !5
      real(r8), allocatable :: gammasun_dbftrop        (:) !6
      real(r8), allocatable :: gammasun_dbftemp        (:) !7
      real(r8), allocatable :: gammasun_dbfboreal      (:) !8
      real(r8), allocatable :: gammasun_ebstemp        (:) !9
      real(r8), allocatable :: gammasun_dbstemp        (:) !10
      real(r8), allocatable :: gammasun_dbsboreal      (:) !11
      real(r8), allocatable :: gammasun_c3arcgrass     (:) !12
      real(r8), allocatable :: gammasun_c3grass        (:) !13
      real(r8), allocatable :: gammasun_c4grass        (:) !14
      real(r8), allocatable :: gammasha_enftemp        (:) !1
      real(r8), allocatable :: gammasha_enfboreal      (:) !2
      real(r8), allocatable :: gammasha_dnfboreal      (:) !3
      real(r8), allocatable :: gammasha_ebftrop        (:) !4
      real(r8), allocatable :: gammasha_ebftemp        (:) !5
      real(r8), allocatable :: gammasha_dbftrop        (:) !6
      real(r8), allocatable :: gammasha_dbftemp        (:) !7
      real(r8), allocatable :: gammasha_dbfboreal      (:) !8
      real(r8), allocatable :: gammasha_ebstemp        (:) !9
      real(r8), allocatable :: gammasha_dbstemp        (:) !10
      real(r8), allocatable :: gammasha_dbsboreal      (:) !11
      real(r8), allocatable :: gammasha_c3arcgrass     (:) !12
      real(r8), allocatable :: gammasha_c3grass        (:) !13
      real(r8), allocatable :: gammasha_c4grass        (:) !14
      real(r8), allocatable :: RuBPlimitfrac_sun_enftemp          (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_enfboreal        (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_dnfboreal        (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_ebftrop          (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_ebftemp          (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_dbftrop          (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_dbftemp          (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_dbfboreal        (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_ebstemp          (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_dbstemp          (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_dbsboreal        (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_c3arcgrass       (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_c3grass          (:)
      real(r8), allocatable :: RuBPlimitfrac_sun_c4grass          (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_enftemp          (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_enfboreal        (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_dnfboreal        (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_ebftrop          (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_ebftemp          (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_dbftrop          (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_dbftemp          (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_dbfboreal        (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_ebstemp          (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_dbstemp          (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_dbsboreal        (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_c3arcgrass       (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_c3grass          (:)
      real(r8), allocatable :: RuBPlimitfrac_sha_c4grass          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_enftemp          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_enfboreal        (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_dnfboreal        (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_ebftrop          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_ebftemp          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_dbftrop          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_dbftemp          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_dbfboreal        (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_ebstemp          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_dbstemp          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_dbsboreal        (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_c3arcgrass       (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_c3grass          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sun_c4grass          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_enftemp          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_enfboreal        (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_dnfboreal        (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_ebftrop          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_ebftemp          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_dbftrop          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_dbftemp          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_dbfboreal        (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_ebstemp          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_dbstemp          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_dbsboreal        (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_c3arcgrass       (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_c3grass          (:)
      real(r8), allocatable :: Rubiscolimitfrac_sha_c4grass          (:)
      real(r8), allocatable :: Sinklimitfrac_sun_enftemp          (:)
      real(r8), allocatable :: Sinklimitfrac_sun_enfboreal        (:)
      real(r8), allocatable :: Sinklimitfrac_sun_dnfboreal        (:)
      real(r8), allocatable :: Sinklimitfrac_sun_ebftrop          (:)
      real(r8), allocatable :: Sinklimitfrac_sun_ebftemp          (:)
      real(r8), allocatable :: Sinklimitfrac_sun_dbftrop          (:)
      real(r8), allocatable :: Sinklimitfrac_sun_dbftemp          (:)
      real(r8), allocatable :: Sinklimitfrac_sun_dbfboreal        (:)
      real(r8), allocatable :: Sinklimitfrac_sun_ebstemp          (:)
      real(r8), allocatable :: Sinklimitfrac_sun_dbstemp          (:)
      real(r8), allocatable :: Sinklimitfrac_sun_dbsboreal        (:)
      real(r8), allocatable :: Sinklimitfrac_sun_c3arcgrass       (:)
      real(r8), allocatable :: Sinklimitfrac_sun_c3grass          (:)
      real(r8), allocatable :: Sinklimitfrac_sun_c4grass          (:)
      real(r8), allocatable :: Sinklimitfrac_sha_enftemp          (:)
      real(r8), allocatable :: Sinklimitfrac_sha_enfboreal        (:)
      real(r8), allocatable :: Sinklimitfrac_sha_dnfboreal        (:)
      real(r8), allocatable :: Sinklimitfrac_sha_ebftrop          (:)
      real(r8), allocatable :: Sinklimitfrac_sha_ebftemp          (:)
      real(r8), allocatable :: Sinklimitfrac_sha_dbftrop          (:)
      real(r8), allocatable :: Sinklimitfrac_sha_dbftemp          (:)
      real(r8), allocatable :: Sinklimitfrac_sha_dbfboreal        (:)
      real(r8), allocatable :: Sinklimitfrac_sha_ebstemp          (:)
      real(r8), allocatable :: Sinklimitfrac_sha_dbstemp          (:)
      real(r8), allocatable :: Sinklimitfrac_sha_dbsboreal        (:)
      real(r8), allocatable :: Sinklimitfrac_sha_c3arcgrass       (:)
      real(r8), allocatable :: Sinklimitfrac_sha_c3grass          (:)
      real(r8), allocatable :: Sinklimitfrac_sha_c4grass          (:)
      real(r8), allocatable :: rstfacsun_enftemp        (:)
      real(r8), allocatable :: rstfacsun_enfboreal      (:)
      real(r8), allocatable :: rstfacsun_dnfboreal      (:)
      real(r8), allocatable :: rstfacsun_ebftrop        (:)
      real(r8), allocatable :: rstfacsun_ebftemp        (:)
      real(r8), allocatable :: rstfacsun_dbftrop        (:)
      real(r8), allocatable :: rstfacsun_dbftemp        (:)
      real(r8), allocatable :: rstfacsun_dbfboreal      (:)
      real(r8), allocatable :: rstfacsun_ebstemp        (:)
      real(r8), allocatable :: rstfacsun_dbstemp        (:)
      real(r8), allocatable :: rstfacsun_dbsboreal      (:)
      real(r8), allocatable :: rstfacsun_c3arcgrass     (:)
      real(r8), allocatable :: rstfacsun_c3grass        (:)
      real(r8), allocatable :: rstfacsun_c4grass        (:)
      real(r8), allocatable :: rstfacsha_enftemp        (:)
      real(r8), allocatable :: rstfacsha_enfboreal      (:)
      real(r8), allocatable :: rstfacsha_dnfboreal      (:)
      real(r8), allocatable :: rstfacsha_ebftrop        (:)
      real(r8), allocatable :: rstfacsha_ebftemp        (:)
      real(r8), allocatable :: rstfacsha_dbftrop        (:)
      real(r8), allocatable :: rstfacsha_dbftemp        (:)
      real(r8), allocatable :: rstfacsha_dbfboreal      (:)
      real(r8), allocatable :: rstfacsha_ebstemp        (:)
      real(r8), allocatable :: rstfacsha_dbstemp        (:)
      real(r8), allocatable :: rstfacsha_dbsboreal      (:)
      real(r8), allocatable :: rstfacsha_c3arcgrass     (:)
      real(r8), allocatable :: rstfacsha_c3grass        (:)
      real(r8), allocatable :: rstfacsha_c4grass        (:)
      real(r8), allocatable :: lambdasun_enftemp        (:)
      real(r8), allocatable :: lambdasun_enfboreal      (:)
      real(r8), allocatable :: lambdasun_dnfboreal      (:)
      real(r8), allocatable :: lambdasun_ebftrop        (:)
      real(r8), allocatable :: lambdasun_ebftemp        (:)
      real(r8), allocatable :: lambdasun_dbftrop        (:)
      real(r8), allocatable :: lambdasun_dbftemp        (:)
      real(r8), allocatable :: lambdasun_dbfboreal      (:)
      real(r8), allocatable :: lambdasun_ebstemp        (:)
      real(r8), allocatable :: lambdasun_dbstemp        (:)
      real(r8), allocatable :: lambdasun_dbsboreal      (:)
      real(r8), allocatable :: lambdasun_c3arcgrass     (:)
      real(r8), allocatable :: lambdasun_c3grass        (:)
      real(r8), allocatable :: lambdasun_c4grass        (:)
      real(r8), allocatable :: lambdasha_enftemp        (:)
      real(r8), allocatable :: lambdasha_enfboreal      (:)
      real(r8), allocatable :: lambdasha_dnfboreal      (:)
      real(r8), allocatable :: lambdasha_ebftrop        (:)
      real(r8), allocatable :: lambdasha_ebftemp        (:)
      real(r8), allocatable :: lambdasha_dbftrop        (:)
      real(r8), allocatable :: lambdasha_dbftemp        (:)
      real(r8), allocatable :: lambdasha_dbfboreal      (:)
      real(r8), allocatable :: lambdasha_ebstemp        (:)
      real(r8), allocatable :: lambdasha_dbstemp        (:)
      real(r8), allocatable :: lambdasha_dbsboreal      (:)
      real(r8), allocatable :: lambdasha_c3arcgrass     (:)
      real(r8), allocatable :: lambdasha_c3grass        (:)
      real(r8), allocatable :: lambdasha_c4grass        (:)
      real(r8), allocatable :: lambda                   (:)
#endif
#endif

      real(r8), allocatable :: tleaf    (:)     ! leaf temperature [K]
      real(r8), allocatable :: ldew     (:)     ! depth of water on foliage [mm]
      real(r8), allocatable :: ldew_rain(:)     ! depth of rain on foliage [mm]
      real(r8), allocatable :: ldew_snow(:)     ! depth of rain on foliage [mm]
      real(r8), allocatable :: sag      (:)     ! non dimensional snow age [-]
      real(r8), allocatable :: scv      (:)     ! snow cover, water equivalent [mm]
      real(r8), allocatable :: snowdp   (:)     ! snow depth [meter]
      real(r8), allocatable :: fveg     (:)     ! fraction of vegetation cover
      real(r8), allocatable :: fsno     (:)     ! fraction of snow cover on ground
      real(r8), allocatable :: sigf     (:)     ! fraction of veg cover, excluding snow-covered veg [-]
      real(r8), allocatable :: green    (:)     ! leaf greenness
      real(r8), allocatable :: tlai     (:)     ! leaf area index
      real(r8), allocatable :: lai      (:)     ! leaf area index
      real(r8), allocatable :: laisun   (:)     ! leaf area index for sunlit leaf
      real(r8), allocatable :: laisha   (:)     ! leaf area index for shaded leaf
      real(r8), allocatable :: tsai     (:)     ! stem area index
      real(r8), allocatable :: sai      (:)     ! stem area index
      real(r8), allocatable :: coszen   (:)     ! cosine of solar zenith angle
      real(r8), allocatable :: alb  (:,:,:)     ! averaged albedo [-]
      real(r8), allocatable :: ssun (:,:,:)     ! sunlit canopy absorption for solar radiation (0-1)
      real(r8), allocatable :: ssha (:,:,:)     ! shaded canopy absorption for solar radiation (0-1)
      real(r8), allocatable :: thermk   (:)     ! canopy gap fraction for tir radiation
      real(r8), allocatable :: extkb    (:)     ! (k, g(mu)/mu) direct solar extinction coefficient
      real(r8), allocatable :: extkd    (:)     ! diffuse and scattered diffuse PAR extinction coefficient
      real(r8), allocatable :: zwt      (:)     ! the depth to water table [m]
      real(r8), allocatable :: wa       (:)     ! water storage in aquifer [mm]
      real(r8), allocatable :: wat      (:)     ! total water storage [mm]
      real(r8), allocatable :: dpond    (:)     ! depth of ponding water [mm]

      real(r8), allocatable :: t_lake(:,:)      ! lake layer teperature [K]
      real(r8), allocatable :: lake_icefrac(:,:)! lake mass fraction of lake layer that is frozen
      real(r8), allocatable :: savedtke1(:)     ! top level eddy conductivity (W/m K) 

      REAL(r8), allocatable :: snw_rds     (:,:) !effective grain radius (col,lyr) [microns, m-6]
      REAL(r8), allocatable :: mss_bcpho   (:,:) !mass of hydrophobic BC in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: mss_bcphi   (:,:) !mass of hydrophillic BC in snow (col,lyr) [kg]
      REAL(r8), allocatable :: mss_ocpho   (:,:) !mass of hydrophobic OC in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: mss_ocphi   (:,:) !mass of hydrophillic OC in snow (col,lyr) [kg]
      REAL(r8), allocatable :: mss_dst1    (:,:) !mass of dust species 1 in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: mss_dst2    (:,:) !mass of dust species 2 in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: mss_dst3    (:,:) !mass of dust species 3 in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: mss_dst4    (:,:) !mass of dust species 4 in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: ssno    (:,:,:,:) !snow layer absorption [-]

      real(r8), allocatable :: trad     (:) ! radiative temperature of surface [K]
      real(r8), allocatable :: tref     (:) ! 2 m height air temperature [kelvin]
      real(r8), allocatable :: qref     (:) ! 2 m height air specific humidity
      real(r8), allocatable :: rst      (:) ! canopy stomatal resistance (s/m)
      real(r8), allocatable :: emis     (:) ! averaged bulk surface emissivity
      real(r8), allocatable :: z0m      (:) ! effective roughness [m]
      real(r8), allocatable :: displa   (:) ! zero displacement height [m]
      real(r8), allocatable :: zol      (:) ! dimensionless height (z/L) used in Monin-Obukhov theory
      real(r8), allocatable :: rib      (:) ! bulk Richardson number in surface layer
      real(r8), allocatable :: ustar    (:) ! u* in similarity theory [m/s]
      real(r8), allocatable :: qstar    (:) ! q* in similarity theory [kg/kg]
      real(r8), allocatable :: tstar    (:) ! t* in similarity theory [K]
      real(r8), allocatable :: fm       (:) ! integral of profile function for momentum
      real(r8), allocatable :: fh       (:) ! integral of profile function for heat
      real(r8), allocatable :: fq       (:) ! integral of profile function for moisture

      ! PUBLIC MEMBER FUNCTIONS:
      public :: allocate_TimeVariables
      public :: deallocate_TimeVariables
      public :: READ_TimeVariables
      public :: WRITE_TimeVariables
#ifdef CLMDEBUG
      public :: check_TimeVariables
#endif


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_TimeVariables 
! --------------------------------------------------------------------
! Allocates memory for CLM 1d [numpatch] variables
! ------------------------------------------------------

  use precision
  USE GlobalVars
  use spmd_task
  use mod_landpatch, only : numpatch
  IMPLICIT NONE


  if (p_is_worker) then

     if (numpatch > 0) then

        allocate (z_sno      (maxsnl+1:0,      numpatch))
        allocate (dz_sno     (maxsnl+1:0,      numpatch))
        allocate (t_soisno   (maxsnl+1:nl_soil,numpatch))
        allocate (wliq_soisno(maxsnl+1:nl_soil,numpatch))
        allocate (wice_soisno(maxsnl+1:nl_soil,numpatch))
        allocate (smp        (1:nl_soil,numpatch))
        allocate (hk         (1:nl_soil,numpatch))
        allocate (h2osoi     (1:nl_soil,numpatch))
        allocate (rootr      (1:nl_soil,numpatch))
#ifdef PLANT_HYDRAULIC_STRESS
        allocate (vegwp      (1:nvegwcs,numpatch))
        allocate (gs0sun               (numpatch))
        allocate (gs0sha               (numpatch))
#endif
#ifdef OzoneStress
        allocate (o3coefv_sun          (numpatch)) ! Ozone stress factor for photosynthesis on sunlit leaf
        allocate (o3coefv_sha          (numpatch)) ! Ozone stress factor for photosynthesis on shaded leaf
        allocate (o3coefg_sun          (numpatch)) ! Ozone stress factor for stomata on sunlit leaf
        allocate (o3coefg_sha          (numpatch)) ! Ozone stress factor for stomata on shaded leaf
        allocate (lai_old              (numpatch)) ! lai in last time step
        allocate (o3uptakesun          (numpatch)) ! Ozone does, sunlit leaf (mmol O3/m^2)
        allocate (o3uptakesha          (numpatch)) ! Ozone does, shaded leaf (mmol O3/m^2)
#endif
        allocate (rstfacsun            (numpatch))
        allocate (rstfacsha            (numpatch))
        allocate (gs_sun            (numpatch))
        allocate (gs_sha            (numpatch))
        allocate (t_grnd               (numpatch))
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
            allocate ( assim_RuBP_sun_enftemp        (numpatch) )
            allocate ( assim_RuBP_sun_enfboreal      (numpatch) )
            allocate ( assim_RuBP_sun_dnfboreal      (numpatch) )
            allocate ( assim_RuBP_sun_ebftrop        (numpatch) )
            allocate ( assim_RuBP_sun_ebftemp        (numpatch) )
            allocate ( assim_RuBP_sun_dbftrop        (numpatch) )
            allocate ( assim_RuBP_sun_dbftemp        (numpatch) )
            allocate ( assim_RuBP_sun_dbfboreal      (numpatch) )
            allocate ( assim_RuBP_sun_ebstemp        (numpatch) )
            allocate ( assim_RuBP_sun_dbstemp        (numpatch) )
            allocate ( assim_RuBP_sun_dbsboreal      (numpatch) )
            allocate ( assim_RuBP_sun_c3arcgrass     (numpatch) )
            allocate ( assim_RuBP_sun_c3grass        (numpatch) )
            allocate ( assim_RuBP_sun_c4grass        (numpatch) )
            allocate ( assim_RuBP_sha_enftemp        (numpatch) )
            allocate ( assim_RuBP_sha_enfboreal      (numpatch) )
            allocate ( assim_RuBP_sha_dnfboreal      (numpatch) )
            allocate ( assim_RuBP_sha_ebftrop        (numpatch) )
            allocate ( assim_RuBP_sha_ebftemp        (numpatch) )
            allocate ( assim_RuBP_sha_dbftrop        (numpatch) )
            allocate ( assim_RuBP_sha_dbftemp        (numpatch) )
            allocate ( assim_RuBP_sha_dbfboreal      (numpatch) )
            allocate ( assim_RuBP_sha_ebstemp        (numpatch) )
            allocate ( assim_RuBP_sha_dbstemp        (numpatch) )
            allocate ( assim_RuBP_sha_dbsboreal      (numpatch) )
            allocate ( assim_RuBP_sha_c3arcgrass     (numpatch) )
            allocate ( assim_RuBP_sha_c3grass        (numpatch) )
            allocate ( assim_RuBP_sha_c4grass        (numpatch) )
            allocate ( assim_Rubisco_sun_enftemp        (numpatch) )
            allocate ( assim_Rubisco_sun_enfboreal      (numpatch) )
            allocate ( assim_Rubisco_sun_dnfboreal      (numpatch) )
            allocate ( assim_Rubisco_sun_ebftrop        (numpatch) )
            allocate ( assim_Rubisco_sun_ebftemp        (numpatch) )
            allocate ( assim_Rubisco_sun_dbftrop        (numpatch) )
            allocate ( assim_Rubisco_sun_dbftemp        (numpatch) )
            allocate ( assim_Rubisco_sun_dbfboreal      (numpatch) )
            allocate ( assim_Rubisco_sun_ebstemp        (numpatch) )
            allocate ( assim_Rubisco_sun_dbstemp        (numpatch) )
            allocate ( assim_Rubisco_sun_dbsboreal      (numpatch) )
            allocate ( assim_Rubisco_sun_c3arcgrass     (numpatch) )
            allocate ( assim_Rubisco_sun_c3grass        (numpatch) )
            allocate ( assim_Rubisco_sun_c4grass        (numpatch) )
            allocate ( assim_Rubisco_sha_enftemp        (numpatch) )
            allocate ( assim_Rubisco_sha_enfboreal      (numpatch) )
            allocate ( assim_Rubisco_sha_dnfboreal      (numpatch) )
            allocate ( assim_Rubisco_sha_ebftrop        (numpatch) )
            allocate ( assim_Rubisco_sha_ebftemp        (numpatch) )
            allocate ( assim_Rubisco_sha_dbftrop        (numpatch) )
            allocate ( assim_Rubisco_sha_dbftemp        (numpatch) )
            allocate ( assim_Rubisco_sha_dbfboreal      (numpatch) )
            allocate ( assim_Rubisco_sha_ebstemp        (numpatch) )
            allocate ( assim_Rubisco_sha_dbstemp        (numpatch) )
            allocate ( assim_Rubisco_sha_dbsboreal      (numpatch) )
            allocate ( assim_Rubisco_sha_c3arcgrass     (numpatch) )
            allocate ( assim_Rubisco_sha_c3grass        (numpatch) )
            allocate ( assim_Rubisco_sha_c4grass        (numpatch) )
            allocate ( assimsun_enftemp        (numpatch) )
            allocate ( assimsun_enfboreal      (numpatch) )
            allocate ( assimsun_dnfboreal      (numpatch) )
            allocate ( assimsun_ebftrop        (numpatch) )
            allocate ( assimsun_ebftemp        (numpatch) )
            allocate ( assimsun_dbftrop        (numpatch) )
            allocate ( assimsun_dbftemp        (numpatch) )
            allocate ( assimsun_dbfboreal      (numpatch) )
            allocate ( assimsun_ebstemp        (numpatch) )
            allocate ( assimsun_dbstemp        (numpatch) )
            allocate ( assimsun_dbsboreal      (numpatch) )
            allocate ( assimsun_c3arcgrass     (numpatch) )
            allocate ( assimsun_c3grass        (numpatch) )
            allocate ( assimsun_c4grass        (numpatch) )
            allocate ( assimsha_enftemp        (numpatch) )
            allocate ( assimsha_enfboreal      (numpatch) )
            allocate ( assimsha_dnfboreal      (numpatch) )
            allocate ( assimsha_ebftrop        (numpatch) )
            allocate ( assimsha_ebftemp        (numpatch) )
            allocate ( assimsha_dbftrop        (numpatch) )
            allocate ( assimsha_dbftemp        (numpatch) )
            allocate ( assimsha_dbfboreal      (numpatch) )
            allocate ( assimsha_ebstemp        (numpatch) )
            allocate ( assimsha_dbstemp        (numpatch) )
            allocate ( assimsha_dbsboreal      (numpatch) )
            allocate ( assimsha_c3arcgrass     (numpatch) )
            allocate ( assimsha_c3grass        (numpatch) )
            allocate ( assimsha_c4grass        (numpatch) )
            allocate ( etrsun_enftemp        (numpatch) )
            allocate ( etrsun_enfboreal      (numpatch) )
            allocate ( etrsun_dnfboreal      (numpatch) )
            allocate ( etrsun_ebftrop        (numpatch) )
            allocate ( etrsun_ebftemp        (numpatch) )
            allocate ( etrsun_dbftrop        (numpatch) )
            allocate ( etrsun_dbftemp        (numpatch) )
            allocate ( etrsun_dbfboreal      (numpatch) )
            allocate ( etrsun_ebstemp        (numpatch) )
            allocate ( etrsun_dbstemp        (numpatch) )
            allocate ( etrsun_dbsboreal      (numpatch) )
            allocate ( etrsun_c3arcgrass     (numpatch) )
            allocate ( etrsun_c3grass        (numpatch) )
            allocate ( etrsun_c4grass        (numpatch) )
            allocate ( etrsha_enftemp        (numpatch) )
            allocate ( etrsha_enfboreal      (numpatch) )
            allocate ( etrsha_dnfboreal      (numpatch) )
            allocate ( etrsha_ebftrop        (numpatch) )
            allocate ( etrsha_ebftemp        (numpatch) )
            allocate ( etrsha_dbftrop        (numpatch) )
            allocate ( etrsha_dbftemp        (numpatch) )
            allocate ( etrsha_dbfboreal      (numpatch) )
            allocate ( etrsha_ebstemp        (numpatch) )
            allocate ( etrsha_dbstemp        (numpatch) )
            allocate ( etrsha_dbsboreal      (numpatch) )
            allocate ( etrsha_c3arcgrass     (numpatch) )
            allocate ( etrsha_c3grass        (numpatch) )
            allocate ( etrsha_c4grass        (numpatch) )
            allocate ( cisun_enftemp        (numpatch) )
            allocate ( cisun_enfboreal      (numpatch) )
            allocate ( cisun_dnfboreal      (numpatch) )
            allocate ( cisun_ebftrop        (numpatch) )
            allocate ( cisun_ebftemp        (numpatch) )
            allocate ( cisun_dbftrop        (numpatch) )
            allocate ( cisun_dbftemp        (numpatch) )
            allocate ( cisun_dbfboreal      (numpatch) )
            allocate ( cisun_ebstemp        (numpatch) )
            allocate ( cisun_dbstemp        (numpatch) )
            allocate ( cisun_dbsboreal      (numpatch) )
            allocate ( cisun_c3arcgrass     (numpatch) )
            allocate ( cisun_c3grass        (numpatch) )
            allocate ( cisun_c4grass        (numpatch) )
            allocate ( cisha_enftemp        (numpatch) )
            allocate ( cisha_enfboreal      (numpatch) )
            allocate ( cisha_dnfboreal      (numpatch) )
            allocate ( cisha_ebftrop        (numpatch) )
            allocate ( cisha_ebftemp        (numpatch) )
            allocate ( cisha_dbftrop        (numpatch) )
            allocate ( cisha_dbftemp        (numpatch) )
            allocate ( cisha_dbfboreal      (numpatch) )
            allocate ( cisha_ebstemp        (numpatch) )
            allocate ( cisha_dbstemp        (numpatch) )
            allocate ( cisha_dbsboreal      (numpatch) )
            allocate ( cisha_c3arcgrass     (numpatch) )
            allocate ( cisha_c3grass        (numpatch) )
            allocate ( cisha_c4grass        (numpatch) )
            allocate ( essun_enftemp        (numpatch) )
            allocate ( essun_enfboreal      (numpatch) )
            allocate ( essun_dnfboreal      (numpatch) )
            allocate ( essun_ebftrop        (numpatch) )
            allocate ( essun_ebftemp        (numpatch) )
            allocate ( essun_dbftrop        (numpatch) )
            allocate ( essun_dbftemp        (numpatch) )
            allocate ( essun_dbfboreal      (numpatch) )
            allocate ( essun_ebstemp        (numpatch) )
            allocate ( essun_dbstemp        (numpatch) )
            allocate ( essun_dbsboreal      (numpatch) )
            allocate ( essun_c3arcgrass     (numpatch) )
            allocate ( essun_c3grass        (numpatch) )
            allocate ( essun_c4grass        (numpatch) )
            allocate ( essha_enftemp        (numpatch) )
            allocate ( essha_enfboreal      (numpatch) )
            allocate ( essha_dnfboreal      (numpatch) )
            allocate ( essha_ebftrop        (numpatch) )
            allocate ( essha_ebftemp        (numpatch) )
            allocate ( essha_dbftrop        (numpatch) )
            allocate ( essha_dbftemp        (numpatch) )
            allocate ( essha_dbfboreal      (numpatch) )
            allocate ( essha_ebstemp        (numpatch) )
            allocate ( essha_dbstemp        (numpatch) )
            allocate ( essha_dbsboreal      (numpatch) )
            allocate ( essha_c3arcgrass     (numpatch) )
            allocate ( essha_c3grass        (numpatch) )
            allocate ( essha_c4grass        (numpatch) )
            allocate ( gssun_enftemp        (numpatch) )
            allocate ( gssun_enfboreal      (numpatch) )
            allocate ( gssun_dnfboreal      (numpatch) )
            allocate ( gssun_ebftrop        (numpatch) )
            allocate ( gssun_ebftemp        (numpatch) )
            allocate ( gssun_dbftrop        (numpatch) )
            allocate ( gssun_dbftemp        (numpatch) )
            allocate ( gssun_dbfboreal      (numpatch) )
            allocate ( gssun_ebstemp        (numpatch) )
            allocate ( gssun_dbstemp        (numpatch) )
            allocate ( gssun_dbsboreal      (numpatch) )
            allocate ( gssun_c3arcgrass     (numpatch) )
            allocate ( gssun_c3grass        (numpatch) )
            allocate ( gssun_c4grass        (numpatch) )
            allocate ( gssha_enftemp        (numpatch) )
            allocate ( gssha_enfboreal      (numpatch) )
            allocate ( gssha_dnfboreal      (numpatch) )
            allocate ( gssha_ebftrop        (numpatch) )
            allocate ( gssha_ebftemp        (numpatch) )
            allocate ( gssha_dbftrop        (numpatch) )
            allocate ( gssha_dbftemp        (numpatch) )
            allocate ( gssha_dbfboreal      (numpatch) )
            allocate ( gssha_ebstemp        (numpatch) )
            allocate ( gssha_dbstemp        (numpatch) )
            allocate ( gssha_dbsboreal      (numpatch) )
            allocate ( gssha_c3arcgrass     (numpatch) )
            allocate ( gssha_c3grass        (numpatch) )
            allocate ( gssha_c4grass        (numpatch) )
            allocate ( gammasun_enftemp        (numpatch) )
            allocate ( gammasun_enfboreal      (numpatch) )
            allocate ( gammasun_dnfboreal      (numpatch) )
            allocate ( gammasun_ebftrop        (numpatch) )
            allocate ( gammasun_ebftemp        (numpatch) )
            allocate ( gammasun_dbftrop        (numpatch) )
            allocate ( gammasun_dbftemp        (numpatch) )
            allocate ( gammasun_dbfboreal      (numpatch) )
            allocate ( gammasun_ebstemp        (numpatch) )
            allocate ( gammasun_dbstemp        (numpatch) )
            allocate ( gammasun_dbsboreal      (numpatch) )
            allocate ( gammasun_c3arcgrass     (numpatch) )
            allocate ( gammasun_c3grass        (numpatch) )
            allocate ( gammasun_c4grass        (numpatch) )
            allocate ( gammasha_enftemp        (numpatch) )
            allocate ( gammasha_enfboreal      (numpatch) )
            allocate ( gammasha_dnfboreal      (numpatch) )
            allocate ( gammasha_ebftrop        (numpatch) )
            allocate ( gammasha_ebftemp        (numpatch) )
            allocate ( gammasha_dbftrop        (numpatch) )
            allocate ( gammasha_dbftemp        (numpatch) )
            allocate ( gammasha_dbfboreal      (numpatch) )
            allocate ( gammasha_ebstemp        (numpatch) )
            allocate ( gammasha_dbstemp        (numpatch) )
            allocate ( gammasha_dbsboreal      (numpatch) )
            allocate ( gammasha_c3arcgrass     (numpatch) )
            allocate ( gammasha_c3grass        (numpatch) )
            allocate ( gammasha_c4grass        (numpatch) )
            allocate ( RuBPlimitfrac_sun_enftemp          (numpatch) )
            allocate ( RuBPlimitfrac_sun_enfboreal        (numpatch) )
            allocate ( RuBPlimitfrac_sun_dnfboreal        (numpatch) )
            allocate ( RuBPlimitfrac_sun_ebftrop          (numpatch) )
            allocate ( RuBPlimitfrac_sun_ebftemp          (numpatch) )
            allocate ( RuBPlimitfrac_sun_dbftrop          (numpatch) )
            allocate ( RuBPlimitfrac_sun_dbftemp          (numpatch) )
            allocate ( RuBPlimitfrac_sun_dbfboreal        (numpatch) )
            allocate ( RuBPlimitfrac_sun_ebstemp          (numpatch) )
            allocate ( RuBPlimitfrac_sun_dbstemp          (numpatch) )
            allocate ( RuBPlimitfrac_sun_dbsboreal        (numpatch) )
            allocate ( RuBPlimitfrac_sun_c3arcgrass       (numpatch) )
            allocate ( RuBPlimitfrac_sun_c3grass          (numpatch) )
            allocate ( RuBPlimitfrac_sun_c4grass          (numpatch) )
            allocate ( RuBPlimitfrac_sha_enftemp          (numpatch) )
            allocate ( RuBPlimitfrac_sha_enfboreal        (numpatch) )
            allocate ( RuBPlimitfrac_sha_dnfboreal        (numpatch) )
            allocate ( RuBPlimitfrac_sha_ebftrop          (numpatch) )
            allocate ( RuBPlimitfrac_sha_ebftemp          (numpatch) )
            allocate ( RuBPlimitfrac_sha_dbftrop          (numpatch) )
            allocate ( RuBPlimitfrac_sha_dbftemp          (numpatch) )
            allocate ( RuBPlimitfrac_sha_dbfboreal        (numpatch) )
            allocate ( RuBPlimitfrac_sha_ebstemp          (numpatch) )
            allocate ( RuBPlimitfrac_sha_dbstemp          (numpatch) )
            allocate ( RuBPlimitfrac_sha_dbsboreal        (numpatch) )
            allocate ( RuBPlimitfrac_sha_c3arcgrass       (numpatch) )
            allocate ( RuBPlimitfrac_sha_c3grass          (numpatch) )
            allocate ( RuBPlimitfrac_sha_c4grass          (numpatch) )
            allocate ( Rubiscolimitfrac_sun_enftemp          (numpatch) )
            allocate ( Rubiscolimitfrac_sun_enfboreal        (numpatch) )
            allocate ( Rubiscolimitfrac_sun_dnfboreal        (numpatch) )
            allocate ( Rubiscolimitfrac_sun_ebftrop          (numpatch) )
            allocate ( Rubiscolimitfrac_sun_ebftemp          (numpatch) )
            allocate ( Rubiscolimitfrac_sun_dbftrop          (numpatch) )
            allocate ( Rubiscolimitfrac_sun_dbftemp          (numpatch) )
            allocate ( Rubiscolimitfrac_sun_dbfboreal        (numpatch) )
            allocate ( Rubiscolimitfrac_sun_ebstemp          (numpatch) )
            allocate ( Rubiscolimitfrac_sun_dbstemp          (numpatch) )
            allocate ( Rubiscolimitfrac_sun_dbsboreal        (numpatch) )
            allocate ( Rubiscolimitfrac_sun_c3arcgrass       (numpatch) )
            allocate ( Rubiscolimitfrac_sun_c3grass          (numpatch) )
            allocate ( Rubiscolimitfrac_sun_c4grass          (numpatch) )
            allocate ( Rubiscolimitfrac_sha_enftemp          (numpatch) )
            allocate ( Rubiscolimitfrac_sha_enfboreal        (numpatch) )
            allocate ( Rubiscolimitfrac_sha_dnfboreal        (numpatch) )
            allocate ( Rubiscolimitfrac_sha_ebftrop          (numpatch) )
            allocate ( Rubiscolimitfrac_sha_ebftemp          (numpatch) )
            allocate ( Rubiscolimitfrac_sha_dbftrop          (numpatch) )
            allocate ( Rubiscolimitfrac_sha_dbftemp          (numpatch) )
            allocate ( Rubiscolimitfrac_sha_dbfboreal        (numpatch) )
            allocate ( Rubiscolimitfrac_sha_ebstemp          (numpatch) )
            allocate ( Rubiscolimitfrac_sha_dbstemp          (numpatch) )
            allocate ( Rubiscolimitfrac_sha_dbsboreal        (numpatch) )
            allocate ( Rubiscolimitfrac_sha_c3arcgrass       (numpatch) )
            allocate ( Rubiscolimitfrac_sha_c3grass          (numpatch) )
            allocate ( Rubiscolimitfrac_sha_c4grass          (numpatch) )
            allocate ( Sinklimitfrac_sun_enftemp          (numpatch) )
            allocate ( Sinklimitfrac_sun_enfboreal        (numpatch) )
            allocate ( Sinklimitfrac_sun_dnfboreal        (numpatch) )
            allocate ( Sinklimitfrac_sun_ebftrop          (numpatch) )
            allocate ( Sinklimitfrac_sun_ebftemp          (numpatch) )
            allocate ( Sinklimitfrac_sun_dbftrop          (numpatch) )
            allocate ( Sinklimitfrac_sun_dbftemp          (numpatch) )
            allocate ( Sinklimitfrac_sun_dbfboreal        (numpatch) )
            allocate ( Sinklimitfrac_sun_ebstemp          (numpatch) )
            allocate ( Sinklimitfrac_sun_dbstemp          (numpatch) )
            allocate ( Sinklimitfrac_sun_dbsboreal        (numpatch) )
            allocate ( Sinklimitfrac_sun_c3arcgrass       (numpatch) )
            allocate ( Sinklimitfrac_sun_c3grass          (numpatch) )
            allocate ( Sinklimitfrac_sun_c4grass          (numpatch) )
            allocate ( Sinklimitfrac_sha_enftemp          (numpatch) )
            allocate ( Sinklimitfrac_sha_enfboreal        (numpatch) )
            allocate ( Sinklimitfrac_sha_dnfboreal        (numpatch) )
            allocate ( Sinklimitfrac_sha_ebftrop          (numpatch) )
            allocate ( Sinklimitfrac_sha_ebftemp          (numpatch) )
            allocate ( Sinklimitfrac_sha_dbftrop          (numpatch) )
            allocate ( Sinklimitfrac_sha_dbftemp          (numpatch) )
            allocate ( Sinklimitfrac_sha_dbfboreal        (numpatch) )
            allocate ( Sinklimitfrac_sha_ebstemp          (numpatch) )
            allocate ( Sinklimitfrac_sha_dbstemp          (numpatch) )
            allocate ( Sinklimitfrac_sha_dbsboreal        (numpatch) )
            allocate ( Sinklimitfrac_sha_c3arcgrass       (numpatch) )
            allocate ( Sinklimitfrac_sha_c3grass          (numpatch) )
            allocate ( Sinklimitfrac_sha_c4grass          (numpatch) )
            allocate ( rstfacsun_enftemp        (numpatch) )
            allocate ( rstfacsun_enfboreal      (numpatch) )
            allocate ( rstfacsun_dnfboreal      (numpatch) )
            allocate ( rstfacsun_ebftrop        (numpatch) )
            allocate ( rstfacsun_ebftemp        (numpatch) )
            allocate ( rstfacsun_dbftrop        (numpatch) )
            allocate ( rstfacsun_dbftemp        (numpatch) )
            allocate ( rstfacsun_dbfboreal      (numpatch) )
            allocate ( rstfacsun_ebstemp        (numpatch) )
            allocate ( rstfacsun_dbstemp        (numpatch) )
            allocate ( rstfacsun_dbsboreal      (numpatch) )
            allocate ( rstfacsun_c3arcgrass     (numpatch) )
            allocate ( rstfacsun_c3grass        (numpatch) )
            allocate ( rstfacsun_c4grass        (numpatch) )
            allocate ( rstfacsha_enftemp        (numpatch) )
            allocate ( rstfacsha_enfboreal      (numpatch) )
            allocate ( rstfacsha_dnfboreal      (numpatch) )
            allocate ( rstfacsha_ebftrop        (numpatch) )
            allocate ( rstfacsha_ebftemp        (numpatch) )
            allocate ( rstfacsha_dbftrop        (numpatch) )
            allocate ( rstfacsha_dbftemp        (numpatch) )
            allocate ( rstfacsha_dbfboreal      (numpatch) )
            allocate ( rstfacsha_ebstemp        (numpatch) )
            allocate ( rstfacsha_dbstemp        (numpatch) )
            allocate ( rstfacsha_dbsboreal      (numpatch) )
            allocate ( rstfacsha_c3arcgrass     (numpatch) )
            allocate ( rstfacsha_c3grass        (numpatch) )
            allocate ( rstfacsha_c4grass        (numpatch) )
            allocate ( lambdasun_enftemp        (numpatch) )
            allocate ( lambdasun_enfboreal      (numpatch) )
            allocate ( lambdasun_dnfboreal      (numpatch) )
            allocate ( lambdasun_ebftrop        (numpatch) )
            allocate ( lambdasun_ebftemp        (numpatch) )
            allocate ( lambdasun_dbftrop        (numpatch) )
            allocate ( lambdasun_dbftemp        (numpatch) )
            allocate ( lambdasun_dbfboreal      (numpatch) )
            allocate ( lambdasun_ebstemp        (numpatch) )
            allocate ( lambdasun_dbstemp        (numpatch) )
            allocate ( lambdasun_dbsboreal      (numpatch) )
            allocate ( lambdasun_c3arcgrass     (numpatch) )
            allocate ( lambdasun_c3grass        (numpatch) )
            allocate ( lambdasun_c4grass        (numpatch) )
            allocate ( lambdasha_enftemp        (numpatch) )
            allocate ( lambdasha_enfboreal      (numpatch) )
            allocate ( lambdasha_dnfboreal      (numpatch) )
            allocate ( lambdasha_ebftrop        (numpatch) )
            allocate ( lambdasha_ebftemp        (numpatch) )
            allocate ( lambdasha_dbftrop        (numpatch) )
            allocate ( lambdasha_dbftemp        (numpatch) )
            allocate ( lambdasha_dbfboreal      (numpatch) )
            allocate ( lambdasha_ebstemp        (numpatch) )
            allocate ( lambdasha_dbstemp        (numpatch) )
            allocate ( lambdasha_dbsboreal      (numpatch) )
            allocate ( lambdasha_c3arcgrass     (numpatch) )
            allocate ( lambdasha_c3grass        (numpatch) )
            allocate ( lambdasha_c4grass        (numpatch) )
            allocate ( lambda                   (numpatch) )
#endif
#endif
        allocate (tleaf                (numpatch))
        allocate (ldew                 (numpatch))
        allocate (ldew_rain            (numpatch))
        allocate (ldew_snow            (numpatch))
        allocate (sag                  (numpatch))
        allocate (scv                  (numpatch))
        allocate (snowdp               (numpatch))
        allocate (fveg                 (numpatch))
        allocate (fsno                 (numpatch))
        allocate (sigf                 (numpatch))
        allocate (green                (numpatch))
        allocate (tlai                 (numpatch))
        allocate (lai                  (numpatch))
        allocate (laisun               (numpatch))
        allocate (laisha               (numpatch))
        allocate (tsai                 (numpatch))
        allocate (sai                  (numpatch))
        allocate (coszen               (numpatch))
        allocate (alb              (2,2,numpatch))
        allocate (ssun             (2,2,numpatch))
        allocate (ssha             (2,2,numpatch))
        allocate (thermk               (numpatch))
        allocate (extkb                (numpatch))
        allocate (extkd                (numpatch))
        allocate (zwt                  (numpatch))
        allocate (wa                   (numpatch))
        allocate (wat                  (numpatch))
        allocate (dpond                (numpatch))
        
        allocate (t_lake       (nl_lake,numpatch))    !new lake scheme
        allocate (lake_icefrac (nl_lake,numpatch))    !new lake scheme
        allocate (savedtke1            (numpatch))    !new lake scheme

        allocate (snw_rds           (maxsnl+1:0,numpatch))
        allocate (mss_bcpho         (maxsnl+1:0,numpatch))
        allocate (mss_bcphi         (maxsnl+1:0,numpatch))
        allocate (mss_ocpho         (maxsnl+1:0,numpatch))
        allocate (mss_ocphi         (maxsnl+1:0,numpatch))
        allocate (mss_dst1          (maxsnl+1:0,numpatch))
        allocate (mss_dst2          (maxsnl+1:0,numpatch))
        allocate (mss_dst3          (maxsnl+1:0,numpatch))
        allocate (mss_dst4          (maxsnl+1:0,numpatch))
        allocate (ssno          (2,2,maxsnl+1:1,numpatch))

        allocate (trad                 (numpatch))
        allocate (tref                 (numpatch))
        allocate (qref                 (numpatch))
        allocate (rst                  (numpatch))
        allocate (emis                 (numpatch))
        allocate (z0m                  (numpatch))
        allocate (displa               (numpatch))
        allocate (zol                  (numpatch))
        allocate (rib                  (numpatch))
        allocate (ustar                (numpatch))
        allocate (qstar                (numpatch))
        allocate (tstar                (numpatch))
        allocate (fm                   (numpatch))
        allocate (fh                   (numpatch))
        allocate (fq                   (numpatch))

     end if
  end if

#ifdef PFT_CLASSIFICATION
     CALL allocate_PFTimeVars
#endif

#ifdef PC_CLASSIFICATION
     CALL allocate_PCTimeVars
#endif

#ifdef BGC 
     CALL allocate_BGCTimeVars
#endif

  END SUBROUTINE allocate_TimeVariables



  SUBROUTINE deallocate_TimeVariables ()

     use spmd_task
     use mod_landpatch, only : numpatch
     implicit none

     ! --------------------------------------------------
     ! Deallocates memory for CLM 1d [numpatch] variables
     ! --------------------------------------------------

     if (p_is_worker) then
        
        if (numpatch > 0) then

           deallocate (z_sno    )
           deallocate (dz_sno   )
           deallocate (t_soisno    )
           deallocate (wliq_soisno )
           deallocate (wice_soisno )
           deallocate (smp )
           deallocate (hk  )
           deallocate (h2osoi )
           deallocate (rootr  )
           deallocate (rstfacsun )
           deallocate (rstfacsha )
           deallocate (gs_sun )
           deallocate (gs_sha )
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
           deallocate ( assim_RuBP_sun_enftemp        )
           deallocate ( assim_RuBP_sun_enfboreal      )
           deallocate ( assim_RuBP_sun_dnfboreal      )
           deallocate ( assim_RuBP_sun_ebftrop        )
           deallocate ( assim_RuBP_sun_ebftemp        )
           deallocate ( assim_RuBP_sun_dbftrop        )
           deallocate ( assim_RuBP_sun_dbftemp        )
           deallocate ( assim_RuBP_sun_dbfboreal      )
           deallocate ( assim_RuBP_sun_ebstemp        )
           deallocate ( assim_RuBP_sun_dbstemp        )
           deallocate ( assim_RuBP_sun_dbsboreal      )
           deallocate ( assim_RuBP_sun_c3arcgrass     )
           deallocate ( assim_RuBP_sun_c3grass        )
           deallocate ( assim_RuBP_sun_c4grass        )
           deallocate ( assim_RuBP_sha_enftemp        )
           deallocate ( assim_RuBP_sha_enfboreal      )
           deallocate ( assim_RuBP_sha_dnfboreal      )
           deallocate ( assim_RuBP_sha_ebftrop        )
           deallocate ( assim_RuBP_sha_ebftemp        )
           deallocate ( assim_RuBP_sha_dbftrop        )
           deallocate ( assim_RuBP_sha_dbftemp        )
           deallocate ( assim_RuBP_sha_dbfboreal      )
           deallocate ( assim_RuBP_sha_ebstemp        )
           deallocate ( assim_RuBP_sha_dbstemp        )
           deallocate ( assim_RuBP_sha_dbsboreal      )
           deallocate ( assim_RuBP_sha_c3arcgrass     )
           deallocate ( assim_RuBP_sha_c3grass        )
           deallocate ( assim_RuBP_sha_c4grass        )
           deallocate ( assim_Rubisco_sun_enftemp        )
           deallocate ( assim_Rubisco_sun_enfboreal      )
           deallocate ( assim_Rubisco_sun_dnfboreal      )
           deallocate ( assim_Rubisco_sun_ebftrop        )
           deallocate ( assim_Rubisco_sun_ebftemp        )
           deallocate ( assim_Rubisco_sun_dbftrop        )
           deallocate ( assim_Rubisco_sun_dbftemp        )
           deallocate ( assim_Rubisco_sun_dbfboreal      )
           deallocate ( assim_Rubisco_sun_ebstemp        )
           deallocate ( assim_Rubisco_sun_dbstemp        )
           deallocate ( assim_Rubisco_sun_dbsboreal      )
           deallocate ( assim_Rubisco_sun_c3arcgrass     )
           deallocate ( assim_Rubisco_sun_c3grass        )
           deallocate ( assim_Rubisco_sun_c4grass        )
           deallocate ( assim_Rubisco_sha_enftemp        )
           deallocate ( assim_Rubisco_sha_enfboreal      )
           deallocate ( assim_Rubisco_sha_dnfboreal      )
           deallocate ( assim_Rubisco_sha_ebftrop        )
           deallocate ( assim_Rubisco_sha_ebftemp        )
           deallocate ( assim_Rubisco_sha_dbftrop        )
           deallocate ( assim_Rubisco_sha_dbftemp        )
           deallocate ( assim_Rubisco_sha_dbfboreal      )
           deallocate ( assim_Rubisco_sha_ebstemp        )
           deallocate ( assim_Rubisco_sha_dbstemp        )
           deallocate ( assim_Rubisco_sha_dbsboreal      )
           deallocate ( assim_Rubisco_sha_c3arcgrass     )
           deallocate ( assim_Rubisco_sha_c3grass        )
           deallocate ( assim_Rubisco_sha_c4grass        )
           deallocate ( assimsun_enftemp        )
           deallocate ( assimsun_enfboreal      )
           deallocate ( assimsun_dnfboreal      )
           deallocate ( assimsun_ebftrop        )
           deallocate ( assimsun_ebftemp        )
           deallocate ( assimsun_dbftrop        )
           deallocate ( assimsun_dbftemp        )
           deallocate ( assimsun_dbfboreal      )
           deallocate ( assimsun_ebstemp        )
           deallocate ( assimsun_dbstemp        )
           deallocate ( assimsun_dbsboreal      )
           deallocate ( assimsun_c3arcgrass     )
           deallocate ( assimsun_c3grass        )
           deallocate ( assimsun_c4grass        )
           deallocate ( assimsha_enftemp        )
           deallocate ( assimsha_enfboreal      )
           deallocate ( assimsha_dnfboreal      )
           deallocate ( assimsha_ebftrop        )
           deallocate ( assimsha_ebftemp        )
           deallocate ( assimsha_dbftrop        )
           deallocate ( assimsha_dbftemp        )
           deallocate ( assimsha_dbfboreal      )
           deallocate ( assimsha_ebstemp        )
           deallocate ( assimsha_dbstemp        )
           deallocate ( assimsha_dbsboreal      )
           deallocate ( assimsha_c3arcgrass     )
           deallocate ( assimsha_c3grass        )
           deallocate ( assimsha_c4grass        )
           deallocate ( etrsun_enftemp        )
           deallocate ( etrsun_enfboreal      )
           deallocate ( etrsun_dnfboreal      )
           deallocate ( etrsun_ebftrop        )
           deallocate ( etrsun_ebftemp        )
           deallocate ( etrsun_dbftrop        )
           deallocate ( etrsun_dbftemp        )
           deallocate ( etrsun_dbfboreal      )
           deallocate ( etrsun_ebstemp        )
           deallocate ( etrsun_dbstemp        )
           deallocate ( etrsun_dbsboreal      )
           deallocate ( etrsun_c3arcgrass     )
           deallocate ( etrsun_c3grass        )
           deallocate ( etrsun_c4grass        )
           deallocate ( etrsha_enftemp        )
           deallocate ( etrsha_enfboreal      )
           deallocate ( etrsha_dnfboreal      )
           deallocate ( etrsha_ebftrop        )
           deallocate ( etrsha_ebftemp        )
           deallocate ( etrsha_dbftrop        )
           deallocate ( etrsha_dbftemp        )
           deallocate ( etrsha_dbfboreal      )
           deallocate ( etrsha_ebstemp        )
           deallocate ( etrsha_dbstemp        )
           deallocate ( etrsha_dbsboreal      )
           deallocate ( etrsha_c3arcgrass     )
           deallocate ( etrsha_c3grass        )
           deallocate ( etrsha_c4grass        )
           deallocate ( cisun_enftemp        )
           deallocate ( cisun_enfboreal      )
           deallocate ( cisun_dnfboreal      )
           deallocate ( cisun_ebftrop        )
           deallocate ( cisun_ebftemp        )
           deallocate ( cisun_dbftrop        )
           deallocate ( cisun_dbftemp        )
           deallocate ( cisun_dbfboreal      )
           deallocate ( cisun_ebstemp        )
           deallocate ( cisun_dbstemp        )
           deallocate ( cisun_dbsboreal      )
           deallocate ( cisun_c3arcgrass     )
           deallocate ( cisun_c3grass        )
           deallocate ( cisun_c4grass        )
           deallocate ( cisha_enftemp        )
           deallocate ( cisha_enfboreal      )
           deallocate ( cisha_dnfboreal      )
           deallocate ( cisha_ebftrop        )
           deallocate ( cisha_ebftemp        )
           deallocate ( cisha_dbftrop        )
           deallocate ( cisha_dbftemp        )
           deallocate ( cisha_dbfboreal      )
           deallocate ( cisha_ebstemp        )
           deallocate ( cisha_dbstemp        )
           deallocate ( cisha_dbsboreal      )
           deallocate ( cisha_c3arcgrass     )
           deallocate ( cisha_c3grass        )
           deallocate ( cisha_c4grass        )
           deallocate ( essun_enftemp        )
           deallocate ( essun_enfboreal      )
           deallocate ( essun_dnfboreal      )
           deallocate ( essun_ebftrop        )
           deallocate ( essun_ebftemp        )
           deallocate ( essun_dbftrop        )
           deallocate ( essun_dbftemp        )
           deallocate ( essun_dbfboreal      )
           deallocate ( essun_ebstemp        )
           deallocate ( essun_dbstemp        )
           deallocate ( essun_dbsboreal      )
           deallocate ( essun_c3arcgrass     )
           deallocate ( essun_c3grass        )
           deallocate ( essun_c4grass        )
           deallocate ( essha_enftemp        )
           deallocate ( essha_enfboreal      )
           deallocate ( essha_dnfboreal      )
           deallocate ( essha_ebftrop        )
           deallocate ( essha_ebftemp        )
           deallocate ( essha_dbftrop        )
           deallocate ( essha_dbftemp        )
           deallocate ( essha_dbfboreal      )
           deallocate ( essha_ebstemp        )
           deallocate ( essha_dbstemp        )
           deallocate ( essha_dbsboreal      )
           deallocate ( essha_c3arcgrass     )
           deallocate ( essha_c3grass        )
           deallocate ( essha_c4grass        )
           deallocate ( gssun_enftemp        )
           deallocate ( gssun_enfboreal      )
           deallocate ( gssun_dnfboreal      )
           deallocate ( gssun_ebftrop        )
           deallocate ( gssun_ebftemp        )
           deallocate ( gssun_dbftrop        )
           deallocate ( gssun_dbftemp        )
           deallocate ( gssun_dbfboreal      )
           deallocate ( gssun_ebstemp        )
           deallocate ( gssun_dbstemp        )
           deallocate ( gssun_dbsboreal      )
           deallocate ( gssun_c3arcgrass     )
           deallocate ( gssun_c3grass        )
           deallocate ( gssun_c4grass        )
           deallocate ( gssha_enftemp        )
           deallocate ( gssha_enfboreal      )
           deallocate ( gssha_dnfboreal      )
           deallocate ( gssha_ebftrop        )
           deallocate ( gssha_ebftemp        )
           deallocate ( gssha_dbftrop        )
           deallocate ( gssha_dbftemp        )
           deallocate ( gssha_dbfboreal      )
           deallocate ( gssha_ebstemp        )
           deallocate ( gssha_dbstemp        )
           deallocate ( gssha_dbsboreal      )
           deallocate ( gssha_c3arcgrass     )
           deallocate ( gssha_c3grass        )
           deallocate ( gssha_c4grass        )
           deallocate ( gammasun_enftemp        )
           deallocate ( gammasun_enfboreal      )
           deallocate ( gammasun_dnfboreal      )
           deallocate ( gammasun_ebftrop        )
           deallocate ( gammasun_ebftemp        )
           deallocate ( gammasun_dbftrop        )
           deallocate ( gammasun_dbftemp        )
           deallocate ( gammasun_dbfboreal      )
           deallocate ( gammasun_ebstemp        )
           deallocate ( gammasun_dbstemp        )
           deallocate ( gammasun_dbsboreal      )
           deallocate ( gammasun_c3arcgrass     )
           deallocate ( gammasun_c3grass        )
           deallocate ( gammasun_c4grass        )
           deallocate ( gammasha_enftemp        )
           deallocate ( gammasha_enfboreal      )
           deallocate ( gammasha_dnfboreal      )
           deallocate ( gammasha_ebftrop        )
           deallocate ( gammasha_ebftemp        )
           deallocate ( gammasha_dbftrop        )
           deallocate ( gammasha_dbftemp        )
           deallocate ( gammasha_dbfboreal      )
           deallocate ( gammasha_ebstemp        )
           deallocate ( gammasha_dbstemp        )
           deallocate ( gammasha_dbsboreal      )
           deallocate ( gammasha_c3arcgrass     )
           deallocate ( gammasha_c3grass        )
           deallocate ( gammasha_c4grass        )
           deallocate ( RuBPlimitfrac_sun_enftemp          )
           deallocate ( RuBPlimitfrac_sun_enfboreal        )
           deallocate ( RuBPlimitfrac_sun_dnfboreal        )
           deallocate ( RuBPlimitfrac_sun_ebftrop          )
           deallocate ( RuBPlimitfrac_sun_ebftemp          )
           deallocate ( RuBPlimitfrac_sun_dbftrop          )
           deallocate ( RuBPlimitfrac_sun_dbftemp          )
           deallocate ( RuBPlimitfrac_sun_dbfboreal        )
           deallocate ( RuBPlimitfrac_sun_ebstemp          )
           deallocate ( RuBPlimitfrac_sun_dbstemp          )
           deallocate ( RuBPlimitfrac_sun_dbsboreal        )
           deallocate ( RuBPlimitfrac_sun_c3arcgrass       )
           deallocate ( RuBPlimitfrac_sun_c3grass          )
           deallocate ( RuBPlimitfrac_sun_c4grass          )
           deallocate ( RuBPlimitfrac_sha_enftemp          )
           deallocate ( RuBPlimitfrac_sha_enfboreal        )
           deallocate ( RuBPlimitfrac_sha_dnfboreal        )
           deallocate ( RuBPlimitfrac_sha_ebftrop          )
           deallocate ( RuBPlimitfrac_sha_ebftemp          )
           deallocate ( RuBPlimitfrac_sha_dbftrop          )
           deallocate ( RuBPlimitfrac_sha_dbftemp          )
           deallocate ( RuBPlimitfrac_sha_dbfboreal        )
           deallocate ( RuBPlimitfrac_sha_ebstemp          )
           deallocate ( RuBPlimitfrac_sha_dbstemp          )
           deallocate ( RuBPlimitfrac_sha_dbsboreal        )
           deallocate ( RuBPlimitfrac_sha_c3arcgrass       )
           deallocate ( RuBPlimitfrac_sha_c3grass          )
           deallocate ( RuBPlimitfrac_sha_c4grass          )
           deallocate ( Rubiscolimitfrac_sun_enftemp          )
           deallocate ( Rubiscolimitfrac_sun_enfboreal        )
           deallocate ( Rubiscolimitfrac_sun_dnfboreal        )
           deallocate ( Rubiscolimitfrac_sun_ebftrop          )
           deallocate ( Rubiscolimitfrac_sun_ebftemp          )
           deallocate ( Rubiscolimitfrac_sun_dbftrop          )
           deallocate ( Rubiscolimitfrac_sun_dbftemp          )
           deallocate ( Rubiscolimitfrac_sun_dbfboreal        )
           deallocate ( Rubiscolimitfrac_sun_ebstemp          )
           deallocate ( Rubiscolimitfrac_sun_dbstemp          )
           deallocate ( Rubiscolimitfrac_sun_dbsboreal        )
           deallocate ( Rubiscolimitfrac_sun_c3arcgrass       )
           deallocate ( Rubiscolimitfrac_sun_c3grass          )
           deallocate ( Rubiscolimitfrac_sun_c4grass          )
           deallocate ( Rubiscolimitfrac_sha_enftemp          )
           deallocate ( Rubiscolimitfrac_sha_enfboreal        )
           deallocate ( Rubiscolimitfrac_sha_dnfboreal        )
           deallocate ( Rubiscolimitfrac_sha_ebftrop          )
           deallocate ( Rubiscolimitfrac_sha_ebftemp          )
           deallocate ( Rubiscolimitfrac_sha_dbftrop          )
           deallocate ( Rubiscolimitfrac_sha_dbftemp          )
           deallocate ( Rubiscolimitfrac_sha_dbfboreal        )
           deallocate ( Rubiscolimitfrac_sha_ebstemp          )
           deallocate ( Rubiscolimitfrac_sha_dbstemp          )
           deallocate ( Rubiscolimitfrac_sha_dbsboreal        )
           deallocate ( Rubiscolimitfrac_sha_c3arcgrass       )
           deallocate ( Rubiscolimitfrac_sha_c3grass          )
           deallocate ( Rubiscolimitfrac_sha_c4grass          )
           deallocate ( Sinklimitfrac_sun_enftemp          )
           deallocate ( Sinklimitfrac_sun_enfboreal        )
           deallocate ( Sinklimitfrac_sun_dnfboreal        )
           deallocate ( Sinklimitfrac_sun_ebftrop          )
           deallocate ( Sinklimitfrac_sun_ebftemp          )
           deallocate ( Sinklimitfrac_sun_dbftrop          )
           deallocate ( Sinklimitfrac_sun_dbftemp          )
           deallocate ( Sinklimitfrac_sun_dbfboreal        )
           deallocate ( Sinklimitfrac_sun_ebstemp          )
           deallocate ( Sinklimitfrac_sun_dbstemp          )
           deallocate ( Sinklimitfrac_sun_dbsboreal        )
           deallocate ( Sinklimitfrac_sun_c3arcgrass       )
           deallocate ( Sinklimitfrac_sun_c3grass          )
           deallocate ( Sinklimitfrac_sun_c4grass          )
           deallocate ( Sinklimitfrac_sha_enftemp          )
           deallocate ( Sinklimitfrac_sha_enfboreal        )
           deallocate ( Sinklimitfrac_sha_dnfboreal        )
           deallocate ( Sinklimitfrac_sha_ebftrop          )
           deallocate ( Sinklimitfrac_sha_ebftemp          )
           deallocate ( Sinklimitfrac_sha_dbftrop          )
           deallocate ( Sinklimitfrac_sha_dbftemp          )
           deallocate ( Sinklimitfrac_sha_dbfboreal        )
           deallocate ( Sinklimitfrac_sha_ebstemp          )
           deallocate ( Sinklimitfrac_sha_dbstemp          )
           deallocate ( Sinklimitfrac_sha_dbsboreal        )
           deallocate ( Sinklimitfrac_sha_c3arcgrass       )
           deallocate ( Sinklimitfrac_sha_c3grass          )
           deallocate ( Sinklimitfrac_sha_c4grass          )
           deallocate ( rstfacsun_enftemp        )
           deallocate ( rstfacsun_enfboreal      )
           deallocate ( rstfacsun_dnfboreal      )
           deallocate ( rstfacsun_ebftrop        )
           deallocate ( rstfacsun_ebftemp        )
           deallocate ( rstfacsun_dbftrop        )
           deallocate ( rstfacsun_dbftemp        )
           deallocate ( rstfacsun_dbfboreal      )
           deallocate ( rstfacsun_ebstemp        )
           deallocate ( rstfacsun_dbstemp        )
           deallocate ( rstfacsun_dbsboreal      )
           deallocate ( rstfacsun_c3arcgrass     )
           deallocate ( rstfacsun_c3grass        )
           deallocate ( rstfacsun_c4grass        )
           deallocate ( rstfacsha_enftemp        )
           deallocate ( rstfacsha_enfboreal      )
           deallocate ( rstfacsha_dnfboreal      )
           deallocate ( rstfacsha_ebftrop        )
           deallocate ( rstfacsha_ebftemp        )
           deallocate ( rstfacsha_dbftrop        )
           deallocate ( rstfacsha_dbftemp        )
           deallocate ( rstfacsha_dbfboreal      )
           deallocate ( rstfacsha_ebstemp        )
           deallocate ( rstfacsha_dbstemp        )
           deallocate ( rstfacsha_dbsboreal      )
           deallocate ( rstfacsha_c3arcgrass     )
           deallocate ( rstfacsha_c3grass        )
           deallocate ( rstfacsha_c4grass        )
           deallocate ( lambdasun_enftemp        )
           deallocate ( lambdasun_enfboreal      )
           deallocate ( lambdasun_dnfboreal      )
           deallocate ( lambdasun_ebftrop        )
           deallocate ( lambdasun_ebftemp        )
           deallocate ( lambdasun_dbftrop        )
           deallocate ( lambdasun_dbftemp        )
           deallocate ( lambdasun_dbfboreal      )
           deallocate ( lambdasun_ebstemp        )
           deallocate ( lambdasun_dbstemp        )
           deallocate ( lambdasun_dbsboreal      )
           deallocate ( lambdasun_c3arcgrass     )
           deallocate ( lambdasun_c3grass        )
           deallocate ( lambdasun_c4grass        )
           deallocate ( lambdasha_enftemp        )
           deallocate ( lambdasha_enfboreal      )
           deallocate ( lambdasha_dnfboreal      )
           deallocate ( lambdasha_ebftrop        )
           deallocate ( lambdasha_ebftemp        )
           deallocate ( lambdasha_dbftrop        )
           deallocate ( lambdasha_dbftemp        )
           deallocate ( lambdasha_dbfboreal      )
           deallocate ( lambdasha_ebstemp        )
           deallocate ( lambdasha_dbstemp        )
           deallocate ( lambdasha_dbsboreal      )
           deallocate ( lambdasha_c3arcgrass     )
           deallocate ( lambdasha_c3grass        )
           deallocate ( lambdasha_c4grass        )
           deallocate ( lambda                   )
#endif
#endif
#ifdef PLANT_HYDRAULIC_STRESS 
           deallocate (vegwp  )
           deallocate (gs0sun )
           deallocate (gs0sha )
#endif
#ifdef OzoneStress
           deallocate (o3coefv_sun) ! Ozone stress factor for photosynthesis on sunlit leaf
           deallocate (o3coefv_sha) ! Ozone stress factor for photosynthesis on shaded leaf
           deallocate (o3coefg_sun) ! Ozone stress factor for stomata on sunlit leaf
           deallocate (o3coefg_sha) ! Ozone stress factor for stomata on shaded leaf
           deallocate (lai_old    ) ! lai in last time step
           deallocate (o3uptakesun) ! Ozone does, sunlit leaf (mmol O3/m^2)
           deallocate (o3uptakesha) ! Ozone does, shaded leaf (mmol O3/m^2)
#endif
           deallocate (t_grnd )
           deallocate (tleaf  )
           deallocate (ldew   )
           deallocate (ldew_rain)
           deallocate (ldew_snow)   
           deallocate (sag    )
           deallocate (scv    )
           deallocate (snowdp )
           deallocate (fveg   )
           deallocate (fsno   )
           deallocate (sigf   )
           deallocate (green  )
           deallocate (tlai   )
           deallocate (lai    )
           deallocate (laisun )
           deallocate (laisha )
           deallocate (tsai   )
           deallocate (sai    )
           deallocate (coszen )
           deallocate (alb    )
           deallocate (ssun   )
           deallocate (ssha   )
           deallocate (thermk )
           deallocate (extkb  )
           deallocate (extkd  )
           deallocate (zwt    )
           deallocate (wa     )
           deallocate (wat    )
           deallocate (dpond  )

           deallocate (t_lake )      ! new lake scheme
           deallocate (lake_icefrac) ! new lake scheme
           deallocate (savedtke1)    ! new lake scheme

           deallocate (snw_rds  )
           deallocate (mss_bcpho)
           deallocate (mss_bcphi)
           deallocate (mss_ocpho)
           deallocate (mss_ocphi)
           deallocate (mss_dst1 )
           deallocate (mss_dst2 )
           deallocate (mss_dst3 )
           deallocate (mss_dst4 )
           deallocate (ssno     )

           deallocate (trad   )
           deallocate (tref   )
           deallocate (qref   )
           deallocate (rst    )
           deallocate (emis   )
           deallocate (z0m    )
           deallocate (displa )
           deallocate (zol    )
           deallocate (rib    )
           deallocate (ustar  )
           deallocate (qstar  )
           deallocate (tstar  )
           deallocate (fm     )
           deallocate (fh     )
           deallocate (fq     )

        end if
     end if

#if (defined PFT_CLASSIFICATION)
     CALL deallocate_PFTimeVars
#endif

#if (defined PC_CLASSIFICATION)
     CALL deallocate_PCTimeVars
#endif

#if (defined BGC)
     CALL deallocate_BGCTimeVars
#endif

  END SUBROUTINE deallocate_TimeVariables


  !---------------------------------------
  function save_to_restart (idate, deltim, itstamp, ptstamp) result(rwrite)

     use mod_namelist
     implicit none

     logical :: rwrite

     integer,  intent(in) :: idate(3)
     real(r8), intent(in) :: deltim
     type(timestamp), intent(in) :: itstamp, ptstamp


     ! added by yuan, 08/31/2014
     select case (trim(DEF_WRST_FREQ))
     case ('TIMESTEP')
        rwrite = .true.
     case ('HOURLY')
        rwrite = isendofhour (idate, deltim)
     case ('DAILY')
        rwrite = isendofday(idate, deltim)
     case ('MONTHLY')
        rwrite = isendofmonth(idate, deltim)       
     case ('YEARLY')
        rwrite = isendofyear(idate, deltim)
     end select

     if (rwrite) then
        rwrite = (ptstamp < itstamp)
     end if

  end function save_to_restart

  !---------------------------------------
  SUBROUTINE WRITE_TimeVariables (idate, site, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use mod_namelist, only : DEF_REST_COMPRESS_LEVEL 
     USE mod_landpatch
     use ncio_vector
     USE GlobalVars
     IMPLICIT NONE

     integer, INTENT(in) :: idate(3)
     character(LEN=*), intent(in) :: site
     character(LEN=*), intent(in) :: dir_restart
     
     ! Local variables
     character(LEN=256) :: file_restart
     character(len=14)  :: cdate
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL 

     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_'//trim(cdate)//'.nc'

     call ncio_create_file_vector (file_restart, landpatch)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')
     
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'snow',     -maxsnl       )
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'snowp1',   -maxsnl+1     )
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'soilsnow', nl_soil-maxsnl)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil',     nl_soil)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'lake',     nl_lake)

#ifdef PLANT_HYDRAULIC_STRESS
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'vegnodes', nvegwcs)
#endif
     
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'band', 2)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'rtyp', 2)

     ! Time-varying state variables which reaquired by restart run
     call ncio_write_vector (file_restart, 'z_sno   '   , 'snow', -maxsnl, 'patch', landpatch, z_sno , compress) !  node depth [m]
     call ncio_write_vector (file_restart, 'dz_sno  '   , 'snow', -maxsnl, 'patch', landpatch, dz_sno, compress) !  interface depth [m]
     call ncio_write_vector (file_restart, 't_soisno'   , 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, t_soisno   , compress) !  soil temperature [K]
     call ncio_write_vector (file_restart, 'wliq_soisno', 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, wliq_soisno, compress) !  liquid water in layers [kg/m2]
     call ncio_write_vector (file_restart, 'wice_soisno', 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, wice_soisno, compress) !  ice lens in layers [kg/m2]
     call ncio_write_vector (file_restart, 'smp',         'soil', nl_soil, 'patch', landpatch, smp, compress) !  soil matrix potential [mm]
     call ncio_write_vector (file_restart, 'hk',          'soil', nl_soil, 'patch', landpatch, hk, compress) !  hydraulic conductivity [mm h2o/s]
#ifdef PLANT_HYDRAULIC_STRESS
     call ncio_write_vector (file_restart, 'vegwp',   'vegnodes', nvegwcs, 'patch', landpatch, vegwp, compress) !  vegetation water potential [mm]
     call ncio_write_vector (file_restart, 'gs0sun',    'patch', landpatch, gs0sun, compress) !  working copy of sunlit stomata conductance
     call ncio_write_vector (file_restart, 'gs0sha',    'patch', landpatch, gs0sha, compress) !  working copy of shalit stomata conductance
#endif
#ifdef OzoneStress
     call ncio_write_vector (file_restart, 'lai_old    ', 'patch', landpatch, lai_old    , compress)
     call ncio_write_vector (file_restart, 'o3uptakesun', 'patch', landpatch, o3uptakesun, compress)
     call ncio_write_vector (file_restart, 'o3uptakesha', 'patch', landpatch, o3uptakesha, compress)
#endif
     call ncio_write_vector (file_restart, 't_grnd  '   , 'patch', landpatch, t_grnd    , compress) !  ground surface temperature [K]
     call ncio_write_vector (file_restart, 'tleaf   '   , 'patch', landpatch, tleaf     , compress) !  leaf temperature [K]
     call ncio_write_vector (file_restart, 'ldew    '   , 'patch', landpatch, ldew      , compress) !  depth of water on foliage [mm]
     call ncio_write_vector (file_restart, 'ldew_rain'  , 'patch', landpatch, ldew_rain , compress) !  depth of water on foliage [mm]
     call ncio_write_vector (file_restart, 'ldew_snow'  , 'patch', landpatch, ldew_snow , compress) !  depth of water on foliage [mm]
     call ncio_write_vector (file_restart, 'sag     '   , 'patch', landpatch, sag       , compress) !  non dimensional snow age [-]
     call ncio_write_vector (file_restart, 'scv     '   , 'patch', landpatch, scv       , compress) !  snow cover, water equivalent [mm]
     call ncio_write_vector (file_restart, 'snowdp  '   , 'patch', landpatch, snowdp    , compress) !  snow depth [meter]
     call ncio_write_vector (file_restart, 'fveg    '   , 'patch', landpatch, fveg      , compress) !  fraction of vegetation cover
     call ncio_write_vector (file_restart, 'fsno    '   , 'patch', landpatch, fsno      , compress) !  fraction of snow cover on ground
     call ncio_write_vector (file_restart, 'sigf    '   , 'patch', landpatch, sigf      , compress) !  fraction of veg cover, excluding snow-covered veg [-]
     call ncio_write_vector (file_restart, 'green   '   , 'patch', landpatch, green     , compress) !  leaf greenness
     call ncio_write_vector (file_restart, 'lai     '   , 'patch', landpatch, lai       , compress) !  leaf area index
     call ncio_write_vector (file_restart, 'tlai    '   , 'patch', landpatch, tlai      , compress) !  leaf area index
     call ncio_write_vector (file_restart, 'sai     '   , 'patch', landpatch, sai       , compress) !  stem area index
     call ncio_write_vector (file_restart, 'tsai    '   , 'patch', landpatch, tsai      , compress) !  stem area index
     call ncio_write_vector (file_restart, 'coszen  '   , 'patch', landpatch, coszen    , compress) !  cosine of solar zenith angle
     call ncio_write_vector (file_restart, 'alb     '   , 'band', 2, 'rtyp', 2, 'patch', landpatch, alb , compress) !  averaged albedo [-]
     call ncio_write_vector (file_restart, 'ssun    '   , 'band', 2, 'rtyp', 2, 'patch', landpatch, ssun, compress) !  sunlit canopy absorption for solar radiation (0-1)
     call ncio_write_vector (file_restart, 'ssha    '   , 'band', 2, 'rtyp', 2, 'patch', landpatch, ssha, compress) !  shaded canopy absorption for solar radiation (0-1)
     call ncio_write_vector (file_restart, 'thermk  '   , 'patch', landpatch, thermk    , compress) !  canopy gap fraction for tir radiation
     call ncio_write_vector (file_restart, 'extkb   '   , 'patch', landpatch, extkb     , compress) !  (k, g(mu)/mu) direct solar extinction coefficient
     call ncio_write_vector (file_restart, 'extkd   '   , 'patch', landpatch, extkd     , compress) !  diffuse and scattered diffuse PAR extinction coefficient
     call ncio_write_vector (file_restart, 'zwt     '   , 'patch', landpatch, zwt       , compress) !  the depth to water table [m]
     call ncio_write_vector (file_restart, 'wa      '   , 'patch', landpatch, wa        , compress) !  water storage in aquifer [mm]
     call ncio_write_vector (file_restart, 'dpond   '   , 'patch', landpatch, dpond     , compress) ! depth of ponding water

     call ncio_write_vector (file_restart, 't_lake  '   , 'lake', nl_lake, 'patch', landpatch, t_lake      , compress) !
     call ncio_write_vector (file_restart, 'lake_icefrc', 'lake', nl_lake, 'patch', landpatch, lake_icefrac, compress) !
     call ncio_write_vector (file_restart, 'savedtke1  ', 'patch', landpatch, savedtke1   , compress) !
     call ncio_write_vector (file_restart, 'snw_rds  ', 'snow', -maxsnl, 'patch', landpatch, snw_rds  , compress) 
     call ncio_write_vector (file_restart, 'mss_bcpho', 'snow', -maxsnl, 'patch', landpatch, mss_bcpho, compress) 
     call ncio_write_vector (file_restart, 'mss_bcphi', 'snow', -maxsnl, 'patch', landpatch, mss_bcphi, compress) 
     call ncio_write_vector (file_restart, 'mss_ocpho', 'snow', -maxsnl, 'patch', landpatch, mss_ocpho, compress) 
     call ncio_write_vector (file_restart, 'mss_ocphi', 'snow', -maxsnl, 'patch', landpatch, mss_ocphi, compress) 
     call ncio_write_vector (file_restart, 'mss_dst1 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst1 , compress) 
     call ncio_write_vector (file_restart, 'mss_dst2 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst2 , compress) 
     call ncio_write_vector (file_restart, 'mss_dst3 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst3 , compress) 
     call ncio_write_vector (file_restart, 'mss_dst4 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst4 , compress) 
     call ncio_write_vector (file_restart, 'ssno', 'band', 2, 'rtyp', 2, 'snowp1', -maxsnl+1, 'patch', landpatch, ssno, compress) 

     ! Additional va_vectorriables required by reginal model (such as WRF ) RSM) 
     call ncio_write_vector (file_restart, 'trad ', 'patch', landpatch, trad , compress) !     radiative temperature of surface [K]
     call ncio_write_vector (file_restart, 'tref ', 'patch', landpatch, tref , compress) !     2 m height air temperature [kelvin]
     call ncio_write_vector (file_restart, 'qref ', 'patch', landpatch, qref , compress) !     2 m height air specific humidity
     call ncio_write_vector (file_restart, 'rst  ', 'patch', landpatch, rst  , compress) !     canopy stomatal resistance (s/m)
     call ncio_write_vector (file_restart, 'emis ', 'patch', landpatch, emis , compress) !     averaged bulk surface emissivity
     call ncio_write_vector (file_restart, 'z0m  ', 'patch', landpatch, z0m  , compress) !     effective roughness [m]
     call ncio_write_vector (file_restart, 'zol  ', 'patch', landpatch, zol  , compress) !     dimensionless height (z/L) used in Monin-Obukhov theory
     call ncio_write_vector (file_restart, 'rib  ', 'patch', landpatch, rib  , compress) !     bulk Richardson number in surface layer
     call ncio_write_vector (file_restart, 'ustar', 'patch', landpatch, ustar, compress) !     u* in similarity theory [m/s]
     call ncio_write_vector (file_restart, 'qstar', 'patch', landpatch, qstar, compress) !     q* in similarity theory [kg/kg]
     call ncio_write_vector (file_restart, 'tstar', 'patch', landpatch, tstar, compress) !     t* in similarity theory [K]
     call ncio_write_vector (file_restart, 'fm   ', 'patch', landpatch, fm   , compress) !     integral of profile function for momentum
     call ncio_write_vector (file_restart, 'fh   ', 'patch', landpatch, fh   , compress) !     integral of profile function for heat
     call ncio_write_vector (file_restart, 'fq   ', 'patch', landpatch, fq   , compress) !     integral of profile function for moisture

#if (defined PFT_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pft_'//trim(cdate)//'.nc'
     CALL WRITE_PFTimeVars (file_restart)
#endif

#if (defined PC_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pc_'//trim(cdate)//'.nc'
     CALL WRITE_PCTimeVars (file_restart)
#endif

#if (defined BGC)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_bgc_'//trim(cdate)//'.nc'
     CALL WRITE_BGCTimeVars (file_restart)
#endif

  end subroutine WRITE_TimeVariables

  !---------------------------------------
  SUBROUTINE READ_TimeVariables (idate, site, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use mod_namelist
     use spmd_task
     use ncio_vector
#ifdef CLMDEBUG 
     USE mod_colm_debug
#endif
     USE mod_landpatch
     USE GlobalVars

     IMPLICIT NONE

     integer, INTENT(in) :: idate(3)
     character(LEN=*), intent(in) :: site
     character(LEN=*), intent(in) :: dir_restart

     ! Local variables
     character(LEN=256) :: file_restart
     character(len=14)  :: cdate
     
#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

     if (p_is_master) then
        write(*,'(/,A26)') 'Loading Time Variables ...'
     end if

     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
     file_restart = trim(dir_restart) // '/' // trim(site) //'_restart_'//trim(cdate)//'.nc'

     ! Time-varying state variables which reaquired by restart run
     call ncio_read_vector (file_restart, 'z_sno   '   , -maxsnl, landpatch, z_sno ) !  node depth [m]
     call ncio_read_vector (file_restart, 'dz_sno  '   , -maxsnl, landpatch, dz_sno) !  interface depth [m]
     call ncio_read_vector (file_restart, 't_soisno'   , nl_soil-maxsnl, landpatch, t_soisno   ) !  soil temperature [K]
     call ncio_read_vector (file_restart, 'wliq_soisno', nl_soil-maxsnl, landpatch, wliq_soisno) !  liquid water in layers [kg/m2]
     call ncio_read_vector (file_restart, 'wice_soisno', nl_soil-maxsnl, landpatch, wice_soisno) !  ice lens in layers [kg/m2]
     call ncio_read_vector (file_restart, 'smp',         nl_soil,        landpatch, smp        ) !  soil matrix potential [mm]
     call ncio_read_vector (file_restart, 'hk',          nl_soil,        landpatch, hk         ) !  hydraulic conductivity [mm h2o/s]
#ifdef PLANT_HYDRAULIC_STRESS
     call ncio_read_vector (file_restart, 'vegwp',       nvegwcs,        landpatch, vegwp      ) !  vegetation water potential [mm]
     call ncio_read_vector (file_restart, 'gs0sun  ',    landpatch, gs0sun     ) !  working copy of sunlit stomata conductance
     call ncio_read_vector (file_restart, 'gs0sha  ',    landpatch, gs0sha     ) !  working copy of shalit stomata conductance
#endif
#ifdef OzoneStress
     call ncio_read_vector (file_restart, 'lai_old    ', landpatch, lai_old    )
     call ncio_read_vector (file_restart, 'o3uptakesun', landpatch, o3uptakesun)
     call ncio_read_vector (file_restart, 'o3uptakesha', landpatch, o3uptakesha)
#endif
     call ncio_read_vector (file_restart, 't_grnd  '   , landpatch, t_grnd     ) !  ground surface temperature [K]
     call ncio_read_vector (file_restart, 'tleaf   '   , landpatch, tleaf      ) !  leaf temperature [K]
     call ncio_read_vector (file_restart, 'ldew    '   , landpatch, ldew       ) !  depth of water on foliage [mm]
     call ncio_read_vector (file_restart, 'ldew_rain    '   , landpatch, ldew_rain       ) !  depth of water on foliage [mm]
     call ncio_read_vector (file_restart, 'ldew_snow    '   , landpatch, ldew_snow       ) !  depth of water on foliage [mm]
     call ncio_read_vector (file_restart, 'sag     '   , landpatch, sag        ) !  non dimensional snow age [-]
     call ncio_read_vector (file_restart, 'scv     '   , landpatch, scv        ) !  snow cover, water equivalent [mm]
     call ncio_read_vector (file_restart, 'snowdp  '   , landpatch, snowdp     ) !  snow depth [meter]
     call ncio_read_vector (file_restart, 'fveg    '   , landpatch, fveg       ) !  fraction of vegetation cover
     call ncio_read_vector (file_restart, 'fsno    '   , landpatch, fsno       ) !  fraction of snow cover on ground
     call ncio_read_vector (file_restart, 'sigf    '   , landpatch, sigf       ) !  fraction of veg cover, excluding snow-covered veg [-]
     call ncio_read_vector (file_restart, 'green   '   , landpatch, green      ) !  leaf greenness
     call ncio_read_vector (file_restart, 'lai     '   , landpatch, lai        ) !  leaf area index
     call ncio_read_vector (file_restart, 'tlai    '   , landpatch, tlai       ) !  leaf area index
     call ncio_read_vector (file_restart, 'sai     '   , landpatch, sai        ) !  stem area index
     call ncio_read_vector (file_restart, 'tsai    '   , landpatch, tsai       ) !  stem area index
     call ncio_read_vector (file_restart, 'coszen  '   , landpatch, coszen     ) !  cosine of solar zenith angle
     call ncio_read_vector (file_restart, 'alb     '   , 2, 2, landpatch, alb  ) !  averaged albedo [-]
     call ncio_read_vector (file_restart, 'ssun    '   , 2, 2, landpatch, ssun ) !  sunlit canopy absorption for solar radiation (0-1)
     call ncio_read_vector (file_restart, 'ssha    '   , 2, 2, landpatch, ssha ) !  shaded canopy absorption for solar radiation (0-1)
     call ncio_read_vector (file_restart, 'thermk  '   , landpatch, thermk     ) !  canopy gap fraction for tir radiation
     call ncio_read_vector (file_restart, 'extkb   '   , landpatch, extkb      ) !  (k, g(mu)/mu) direct solar extinction coefficient
     call ncio_read_vector (file_restart, 'extkd   '   , landpatch, extkd      ) !  diffuse and scattered diffuse PAR extinction coefficient
     call ncio_read_vector (file_restart, 'zwt     '   , landpatch, zwt        ) !  the depth to water table [m]
     call ncio_read_vector (file_restart, 'wa      '   , landpatch, wa         ) !  water storage in aquifer [mm]
     call ncio_read_vector (file_restart, 'dpond   '   , landpatch, dpond      ) ! depth of ponding water

     call ncio_read_vector (file_restart, 't_lake  '   , nl_lake, landpatch, t_lake      ) !
     call ncio_read_vector (file_restart, 'lake_icefrc', nl_lake, landpatch, lake_icefrac) !
     call ncio_read_vector (file_restart, 'savedtke1', landpatch, savedtke1) !
     
     call ncio_read_vector (file_restart, 'snw_rds  ', -maxsnl, landpatch, snw_rds  ) !
     call ncio_read_vector (file_restart, 'mss_bcpho', -maxsnl, landpatch, mss_bcpho) !
     call ncio_read_vector (file_restart, 'mss_bcphi', -maxsnl, landpatch, mss_bcphi) !
     call ncio_read_vector (file_restart, 'mss_ocpho', -maxsnl, landpatch, mss_ocpho) !
     call ncio_read_vector (file_restart, 'mss_ocphi', -maxsnl, landpatch, mss_ocphi) !
     call ncio_read_vector (file_restart, 'mss_dst1 ', -maxsnl, landpatch, mss_dst1 ) !
     call ncio_read_vector (file_restart, 'mss_dst2 ', -maxsnl, landpatch, mss_dst2 ) !
     call ncio_read_vector (file_restart, 'mss_dst3 ', -maxsnl, landpatch, mss_dst3 ) !
     call ncio_read_vector (file_restart, 'mss_dst4 ', -maxsnl, landpatch, mss_dst4 ) !
     call ncio_read_vector (file_restart, 'ssno', 2,2, -maxsnl+1, landpatch, ssno) !

     ! Additional variables required by reginal model (such as WRF ) RSM) 
     call ncio_read_vector (file_restart, 'trad ', landpatch, trad ) !     radiative temperature of surface [K]
     call ncio_read_vector (file_restart, 'tref ', landpatch, tref ) !     2 m height air temperature [kelvin]
     call ncio_read_vector (file_restart, 'qref ', landpatch, qref ) !     2 m height air specific humidity
     call ncio_read_vector (file_restart, 'rst  ', landpatch, rst  ) !     canopy stomatal resistance (s/m)
     call ncio_read_vector (file_restart, 'emis ', landpatch, emis ) !     averaged bulk surface emissivity
     call ncio_read_vector (file_restart, 'z0m  ', landpatch, z0m  ) !     effective roughness [m]
     call ncio_read_vector (file_restart, 'zol  ', landpatch, zol  ) !     dimensionless height (z/L) used in Monin-Obukhov theory
     call ncio_read_vector (file_restart, 'rib  ', landpatch, rib  ) !     bulk Richardson number in surface layer
     call ncio_read_vector (file_restart, 'ustar', landpatch, ustar) !     u* in similarity theory [m/s]
     call ncio_read_vector (file_restart, 'qstar', landpatch, qstar) !     q* in similarity theory [kg/kg]
     call ncio_read_vector (file_restart, 'tstar', landpatch, tstar) !     t* in similarity theory [K]
     call ncio_read_vector (file_restart, 'fm   ', landpatch, fm   ) !     integral of profile function for momentum
     call ncio_read_vector (file_restart, 'fh   ', landpatch, fh   ) !     integral of profile function for heat
     call ncio_read_vector (file_restart, 'fq   ', landpatch, fq   ) !     integral of profile function for moisture

#if (defined PFT_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pft_'//trim(cdate)//'.nc'
     CALL READ_PFTimeVars (file_restart)
#endif

#if (defined PC_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pc_'//trim(cdate)//'.nc'
     CALL READ_PCTimeVars (file_restart)
#endif
     
#if (defined BGC)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_bgc_'//trim(cdate)//'.nc'
     CALL READ_BGCTimeVars (file_restart)
#endif

#ifdef CLMDEBUG
     call check_TimeVariables
#endif
     
     if (p_is_master) then
        write(*,*) 'Loading Time Variables done.'
     end if

  end subroutine READ_TimeVariables

  !---------------------------------------
#ifdef CLMDEBUG
  SUBROUTINE check_TimeVariables ()

     use spmd_task
     use mod_colm_debug

     IMPLICIT NONE

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif
     if (p_is_master) then
        write(*,'(/,A27)') 'Checking Time Variables ...'
     end if

     call check_vector_data ('z_sno       ', z_sno )      !  node depth [m]
     call check_vector_data ('dz_sno      ', dz_sno)      !  interface depth [m]
     call check_vector_data ('t_soisno    ', t_soisno   ) !  soil temperature [K]
     call check_vector_data ('wliq_soisno ', wliq_soisno) !  liquid water in layers [kg/m2]
     call check_vector_data ('wice_soisno ', wice_soisno) !  ice lens in layers [kg/m2]
     call check_vector_data ('smp         ', smp        ) !  soil matrix potential [mm]
     call check_vector_data ('hk          ', hk         ) !  hydraulic conductivity [mm h2o/s]
#ifdef PLANT_HYDRAULIC_STRESS
     call check_vector_data ('vegwp       ', vegwp      ) !  vegetation water potential [mm]
     call check_vector_data ('gs0sun      ', gs0sun     ) !  working copy of sunlit stomata conductance
     call check_vector_data ('gs0sha      ', gs0sha     ) !  working copy of shalit stomata conductance
#endif
#ifdef OzoneStress
     call check_vector_data ('o3coefv_sun', o3coefv_sun)
     call check_vector_data ('o3coefv_sha', o3coefv_sha)
     call check_vector_data ('o3coefg_sun', o3coefg_sun)
     call check_vector_data ('o3coefg_sha', o3coefg_sha)
     call check_vector_data ('lai_old    ', lai_old    )
     call check_vector_data ('o3uptakesun', o3uptakesun)
     call check_vector_data ('o3uptakesha', o3uptakesha)
#endif
     call check_vector_data ('t_grnd      ', t_grnd     ) !  ground surface temperature [K]
     call check_vector_data ('tleaf       ', tleaf      ) !  leaf temperature [K]
     call check_vector_data ('ldew        ', ldew       ) !  depth of water on foliage [mm]
     call check_vector_data ('ldew_rain        ', ldew_rain       ) !  depth of water on foliage [mm]
     call check_vector_data ('ldew_snow        ', ldew_snow       ) !  depth of water on foliage [mm]
     call check_vector_data ('sag         ', sag        ) !  non dimensional snow age [-]
     call check_vector_data ('scv         ', scv        ) !  snow cover, water equivalent [mm]
     call check_vector_data ('snowdp      ', snowdp     ) !  snow depth [meter]
     call check_vector_data ('fveg        ', fveg       ) !  fraction of vegetation cover
     call check_vector_data ('fsno        ', fsno       ) !  fraction of snow cover on ground
     call check_vector_data ('sigf        ', sigf       ) !  fraction of veg cover, excluding snow-covered veg [-]
     call check_vector_data ('green       ', green      ) !  leaf greenness
     call check_vector_data ('lai         ', lai        ) !  leaf area index
     call check_vector_data ('tlai        ', tlai       ) !  leaf area index
     call check_vector_data ('sai         ', sai        ) !  stem area index
     call check_vector_data ('tsai        ', tsai       ) !  stem area index
     call check_vector_data ('coszen      ', coszen     ) !  cosine of solar zenith angle
     call check_vector_data ('alb         ', alb  ) !  averaged albedo [-]
     call check_vector_data ('ssun        ', ssun ) !  sunlit canopy absorption for solar radiation (0-1)
     call check_vector_data ('ssha        ', ssha ) !  shaded canopy absorption for solar radiation (0-1)
     call check_vector_data ('thermk      ', thermk     ) !  canopy gap fraction for tir radiation
     call check_vector_data ('extkb       ', extkb      ) !  (k, g(mu)/mu) direct solar extinction coefficient
     call check_vector_data ('extkd       ', extkd      ) !  diffuse and scattered diffuse PAR extinction coefficient
     call check_vector_data ('zwt         ', zwt        ) !  the depth to water table [m]
     call check_vector_data ('wa          ', wa         ) !  water storage in aquifer [mm]
     call check_vector_data ('dpond       ', dpond      ) !  depth of ponding water

     call check_vector_data ('t_lake      ', t_lake      ) !
     call check_vector_data ('lake_icefrc ', lake_icefrac) !
     call check_vector_data ('savedtke1   ', savedtke1   ) !

#if (defined PFT_CLASSIFICATION)
     CALL check_PFTimeVars
#endif

#if (defined PC_CLASSIFICATION)
     CALL check_PCTimeVars
#endif

#if (defined BGC)
     CALL check_BGCTimeVars
#endif

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif
     
  end subroutine check_TimeVariables
#endif


END MODULE MOD_TimeVariables
! ------ EOP --------------
