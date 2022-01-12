#include <define.h> 

MODULE MOD_TimeInvariants 
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

  USE precision
  USE GlobalVars, only: ndecomp_transitions, ndecomp_pools
  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
! surface classification and soil information
  INTEGER,  allocatable :: patch2lon      (:)  !patch longitude index
  INTEGER,  allocatable :: patch2lat      (:)  !patch latitude index
  INTEGER,  allocatable :: patchclass     (:)  !index of land cover type of the patches at the fraction > 0
  INTEGER,  allocatable :: patchtype      (:)  !land water type
  INTEGER,  allocatable :: grid_patch_s (:,:)  !start patch number of grid
  INTEGER,  allocatable :: grid_patch_e (:,:)  !end patch number of grid
  REAL(r8), allocatable :: gridarea     (:,:)  !area of grid cell the patch located in
  REAL(r8), allocatable :: patchfrac      (:)  !patch weight

  REAL(r8), allocatable :: patchlatr      (:)  !latitude in radians
  REAL(r8), allocatable :: patchlonr      (:)  !longitude in radians
  REAL(r8), allocatable :: gridlatd       (:)  !latitude in degrees
  REAL(r8), allocatable :: gridlond       (:)  !longitude in degrees

  REAL(r8), allocatable :: lakedepth      (:)  !lake depth
  REAL(r8), allocatable :: dz_lake      (:,:)  !new lake scheme

  REAL(r8), allocatable :: soil_s_v_alb   (:)  !albedo of visible of the saturated soil
  REAL(r8), allocatable :: soil_d_v_alb   (:)  !albedo of visible of the dry soil
  REAL(r8), allocatable :: soil_s_n_alb   (:)  !albedo of near infrared of the saturated soil
  REAL(r8), allocatable :: soil_d_n_alb   (:)  !albedo of near infrared of the dry soil
  REAL(r8), allocatable :: porsl        (:,:)  !fraction of soil that is voids [-]
  REAL(r8), allocatable :: psi0         (:,:)  !minimum soil suction [mm] (NOTE: "-" valued)
  REAL(r8), allocatable :: bsw          (:,:)  !clapp and hornbereger "b" parameter [-]
  REAL(r8), allocatable :: hksati       (:,:)  !hydraulic conductivity at saturation [mm h2o/s]
  REAL(r8), allocatable :: csol         (:,:)  !heat capacity of soil solids [J/(m3 K)]
  REAL(r8), allocatable :: dksatu       (:,:)  !thermal conductivity of saturated soil [W/m-K]
  REAL(r8), allocatable :: dkdry        (:,:)  !thermal conductivity for dry soil  [W/(m-K)]

  REAL(r8), allocatable :: htop           (:)  !canopy top height [m]
  REAL(r8), allocatable :: hbot           (:)  !canopy bottom height [m]


  REAL(r8) :: zlnd         !roughness length for soil [m]
  REAL(r8) :: zsno         !roughness length for snow [m]
  REAL(r8) :: csoilc       !drag coefficient for soil under canopy [-]
  REAL(r8) :: dewmx        !maximum dew
  REAL(r8) :: wtfact       !fraction of model area with high water table
  REAL(r8) :: capr         !tuning factor to turn first layer T into surface T
  REAL(r8) :: cnfac        !Crank Nicholson factor between 0 and 1
  REAL(r8) :: ssi          !irreducible water saturation of snow
  REAL(r8) :: wimp         !water impremeable if porosity less than wimp
  REAL(r8) :: pondmx       !ponding depth (mm)
  REAL(r8) :: smpmax       !wilting point potential in mm
  REAL(r8) :: smpmin       !restriction for min of soil poten. (mm)
  REAL(r8) :: trsmx0       !max transpiration for moist soil+100% veg.  [mm/s]
  REAL(r8) :: tcrit        !critical temp. to determine rain or snow

! bgc constant
  INTEGER  :: donor_pool       (1:ndecomp_transitions)
  INTEGER  :: receiver_pool    (1:ndecomp_transitions)
  REAL(r8) :: am
  LOGICAL  :: floating_cn_ratio(1:ndecomp_pools)
  REAL(r8) :: initial_cn_ratio (1:ndecomp_pools)
  REAL(r8), allocatable :: rf_decomp        (:,:,:)
  REAL(r8), allocatable :: pathfrac_decomp  (:,:,:)

  INTEGER  :: i_met_lit
  INTEGER  :: i_cel_lit
  INTEGER  :: i_lig_lit
  INTEGER  :: i_cwd
  INTEGER  :: i_soil1
  INTEGER  :: i_soil2
  INTEGER  :: i_soil3
  INTEGER  :: i_atm

  LOGICAL  :: is_cwd    (1:ndecomp_pools) ! True => is a coarse woody debris pool
  LOGICAL  :: is_litter (1:ndecomp_pools) ! True => is a litter pool
  LOGICAL  :: is_soil   (1:ndecomp_pools) ! True => is a soil pool

  REAL(r8), allocatable :: gdp_lf (:) !
  REAL(r8), allocatable :: abm_lf (:) !
  REAL(r8), allocatable :: peatf_lf(:)!
  REAL(r8) :: cmb_cmplt_fact(1:2)

  REAL(r8) :: nitrif_n2o_loss_frac ! fraction of N lost as N2O in nitrification (Li et al., 2000)
  REAL(r8) :: dnp     ! denitrification proportion
  REAL(r8) :: bdnr    ! bulk denitrification rate (1/day)
  REAL(r8) :: Q10
  REAL(r8) :: froz_q10 
  REAL(r8) :: tau_l1
  REAL(r8) :: tau_l2_l3
  REAL(r8) :: tau_s1   
  REAL(r8) :: tau_s2  
  REAL(r8) :: tau_s3 
  REAL(r8) :: tau_cwd
  REAL(r8) :: lwtop

  REAL(r8) :: som_adv_flux
  REAL(r8) :: som_diffus
  REAL(r8) :: cryoturb_diffusion_k
  REAL(r8) :: max_altdepth_cryoturbation
  REAL(r8) :: max_depth_cryoturb

  REAL(r8) :: br         ! basal respiration rate for aboveground biomass
  REAL(r8) :: br_root    ! basal respiration rate for belowground biomass
  REAL(r8) :: fstor2tran
  REAL(r8) :: ndays_on
  REAL(r8) :: ndays_off
  REAL(r8) :: crit_dayl
  REAL(r8) :: crit_onset_fdd
  REAL(r8) :: crit_onset_swi
  REAL(r8) :: crit_offset_fdd
  REAL(r8) :: crit_offset_swi
  REAL(r8) :: soilpsi_on
  REAL(r8) :: soilpsi_off

  REAL(r8) :: occur_hi_gdp_tree
  REAL(r8) :: lfuel
  REAL(r8) :: ufuel
  REAL(r8) :: cropfire_a1
  REAL(r8) :: borealat
  REAL(r8) :: troplat
  REAL(r8) :: non_boreal_peatfire_c
  REAL(r8) :: boreal_peatfire_c
  REAL(r8) :: rh_low
  REAL(r8) :: rh_hgh
  REAL(r8) :: bt_min
  REAL(r8) :: bt_max
  REAL(r8) :: pot_hmn_ign_counts_alpha
  REAL(r8) :: g0

  REAL(r8) :: sf
  REAL(r8) :: sf_no3

!----

          
! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_TimeInvariants
  PUBLIC :: deallocate_TimeInvariants
  PUBLIC :: WRITE_TimeInvariants
  PUBLIC :: READ_TimeInvariants

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_TimeInvariants (lon_points, lat_points)
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM 1d [numpatch] variables
  ! --------------------------------------------------------------------

     USE precision
     USE GlobalVars
     USE MOD_PFTimeInvars
     USE MOD_PCTimeInvars
     IMPLICIT NONE

     INTEGER, intent(in) :: lon_points
     INTEGER, intent(in) :: lat_points

     allocate (patch2lon                 (numpatch))
     allocate (patch2lat                 (numpatch))
     allocate (patchclass                (numpatch))
     allocate (patchtype                 (numpatch))
     allocate (grid_patch_s (lon_points,lat_points))
     allocate (grid_patch_e (lon_points,lat_points))
     allocate (gridarea     (lon_points,lat_points))
     allocate (patchfrac                 (numpatch))

     allocate (patchlatr                 (numpatch))
     allocate (patchlonr                 (numpatch))
     allocate (gridlatd                (lat_points))
     allocate (gridlond                (lon_points))

     allocate (lakedepth                 (numpatch))
     allocate (dz_lake           (nl_lake,numpatch))

     allocate (soil_s_v_alb              (numpatch))
     allocate (soil_d_v_alb              (numpatch))
     allocate (soil_s_n_alb              (numpatch))
     allocate (soil_d_n_alb              (numpatch))
     allocate (porsl             (nl_soil,numpatch))
     allocate (psi0              (nl_soil,numpatch))
     allocate (bsw               (nl_soil,numpatch))
     allocate (hksati            (nl_soil,numpatch))
     allocate (csol              (nl_soil,numpatch))
     allocate (dksatu            (nl_soil,numpatch))
     allocate (dkdry             (nl_soil,numpatch))

     allocate (htop                      (numpatch))
     allocate (hbot                      (numpatch))

! bgc varaibles
     allocate (rf_decomp         (nl_soil,ndecomp_transitions,numpatch))
     allocate (pathfrac_decomp   (nl_soil,ndecomp_transitions,numpatch))
     allocate (gdp_lf            (numpatch))
     allocate (abm_lf            (numpatch))
     allocate (peatf_lf          (numpatch))

! end bgc variables

#ifdef PFT_CLASSIFICATION
     CALL allocate_PFTimeInvars
#endif

#ifdef PC_CLASSIFICATION
     CALL allocate_PCTimeInvars
#endif

  END SUBROUTINE allocate_TimeInvariants


  SUBROUTINE READ_TimeInvariants(dir_restart_hist,casename)
! --------------------------------------------------------------------
! Write out as a restart file [histTimeConst]
! ...............................................
     USE precision
     USE GlobalVars
     USE MOD_PFTimeInvars
     USE MOD_PCTimeInvars
     IMPLICIT NONE

     CHARACTER(LEN=256), intent(in) :: casename           !casename name
     CHARACTER(LEN=256), intent(in) :: dir_restart_hist

     CHARACTER(LEN=256) :: fhistTimeConst
     INTEGER :: lhistTimeConst

     lhistTimeConst = 100
     fhistTimeConst = trim(dir_restart_hist)//trim(casename)//'-'//'rstTimeConst'
     open(unit=lhistTimeConst,file=trim(fhistTimeConst),status='unknown',&
                              form='unformatted',action='read')

     read (lhistTimeConst)  &!
           patch2lon,       &! longitude index for each patch point
           patch2lat,       &! latitude index for each patch point
           patchclass,      &! index of land cover type of the patches at the fraction > 0
           patchtype,       &! land water TYPE
           grid_patch_s,    &! start patch number of grid
           grid_patch_e,    &! end patch number of grid
           gridarea,        &! area of grid the patch located in
           patchfrac         ! subgrid weight for each patch point

     read (lhistTimeConst)  &!
           patchlatr,       &! latitude in radians
           patchlonr,       &! longitude in radians
           gridlatd,        &! latitude in degrees
           gridlond,        &! longitude in degrees

     ! Soil and plant parameters OF CLM
           lakedepth,       &! lake depth
           dz_lake,         &! new lake scheme

           soil_s_v_alb,    &! albedo of visible of the saturated soil
           soil_d_v_alb,    &! albedo of visible of the dry soil
           soil_s_n_alb,    &! albedo of near infrared of the saturated soil
           soil_d_n_alb,    &! albedo of near infrared of the dry soil
           porsl,           &! fraction of soil that is voids [-]
           psi0,            &! minimum soil suction [mm] (NOTE: "-" valued)
           bsw,             &! clapp and hornbereger "b" parameter [-]
           hksati,          &! hydraulic conductivity at saturation [mm h2o/s]
           csol,            &! heat capacity of soil solids [J/(m3 K)]
           dksatu,          &! thermal conductivity of saturated soil [W/m-K]
           dkdry,           &! thermal conductivity for dry soil  [W/(m-K)]

           htop,            &! canopy top height [m]
           hbot,            &! canopy bottom height [m]

     ! CLM TUNABLE constants
           zlnd,            &! roughness length for soil [m]
           zsno,            &! roughness length for snow [m]
           csoilc,          &! drag coefficient for soil under canopy [-]
           dewmx,           &! maximum dew
           wtfact,          &! fraction of model area with high water table
           capr,            &! tuning factor to turn first layer T into surface T
           cnfac,           &! Crank Nicholson factor between 0 and 1
           ssi,             &! irreducible water saturation of snow
           wimp,            &! water impremeable if porosity less than wimp
           pondmx,          &! ponding depth (mm)
           smpmax,          &! wilting point potential in mm
           smpmin,          &! restriction for min of soil poten. (mm)
           trsmx0,          &! max transpiration for moist soil+100% veg.  [mm/s]
           tcrit,           &! critical temp. to determine rain or snow

! bgc constants
           donor_pool     , &
           receiver_pool  , &
           floating_cn_ratio,& ! TRUE => pool has fixed C:N ratio
           initial_cn_ratio, & ! c:n ratio for initialization of pools
           rf_decomp       , & ! (frac) respired fraction in decomposition step
           pathfrac_decomp , & ! what fraction of C leaving a given pool passes through a given transition (frac)

           i_met_lit      , &
           i_cel_lit      , &
           i_lig_lit      , &
           i_cwd          , &
           i_soil1        , &
           i_soil2        , &
           i_soil3        , &
           i_atm          , &
           is_cwd         , &
           is_litter      , &
           is_soil        , &

           gdp_lf         , &
           abm_lf         , &
           peatf_lf       , &
           cmb_cmplt_fact , &

           nitrif_n2o_loss_frac, &! fraction of N lost as N2O in nitrification (Li et al., 2000)
           dnp                 , &!
           bdnr                , &!
           Q10                 , &!
           froz_q10            , &!
           tau_l1              , &!
           tau_l2_l3           , &!
           tau_s1              , &!
           tau_s2              , &!
           tau_s3              , &!
           tau_cwd             , &
           lwtop               , &

           som_adv_flux        , &
           som_diffus          , &
           cryoturb_diffusion_k, &
           max_altdepth_cryoturbation, &
           max_depth_cryoturb  , &

           am                  , &
           br                  , &
           br_root             , &
           fstor2tran          , &
           ndays_on            , &
           ndays_off           , &
           crit_dayl           , &
           crit_onset_fdd      , &
           crit_onset_swi      , &
           crit_offset_fdd     , &
           crit_offset_swi     , &
           soilpsi_on          , &
           soilpsi_off         , &

           occur_hi_gdp_tree   , &
           lfuel               , &
           ufuel               , &
           cropfire_a1         , &
           borealat            , &
           troplat             , &
           non_boreal_peatfire_c, &
           boreal_peatfire_c   , &
           rh_low              , &
           rh_hgh              , &
           bt_min              , &
           bt_max              , &
           pot_hmn_ign_counts_alpha, &
           g0, &

           sf, &
           sf_no3

     
     ! PFT/PC time invariants
#ifdef PFT_CLASSIFICATION
     read (lhistTimeConst)  &!
           pftclass,        &! PFT type
           pftfrac,         &! PFT fractional cover
           patch_pft_s,     &! patch start index of PFT
           patch_pft_e,     &! patch end index of PFT
           pft2patch,       &! projection from PFT to patch
           htop_p,          &! canopy top height [m]
           hbot_p            ! canopy bottom height [m]
#endif

#ifdef PC_CLASSIFICATION
     read (lhistTimeConst)  &!
           patch2pc,        &! projection from patch to PC
           pc2patch,        &! projection from PC to patch
           pcfrac,          &! PC fractional cover
           htop_c,          &! canopy top height [m]
           hbot_c            ! canopy bottom height [m]   
#endif

     close(lhistTimeConst)

  END SUBROUTINE READ_TimeInvariants


  SUBROUTINE WRITE_TimeInvariants(dir_restart_hist,casename)
! --------------------------------------------------------------------
! Write out as a restart file [histTimeConst]
! ...............................................
     USE precision
     USE GlobalVars
     USE MOD_PFTimeInvars
     USE MOD_PCTimeInvars
     IMPLICIT NONE

     CHARACTER(LEN=256), intent(in) :: casename           !casename name
     CHARACTER(LEN=256), intent(in) :: dir_restart_hist

     CHARACTER(LEN=256) :: fhistTimeConst
     INTEGER :: lhistTimeConst

     lhistTimeConst = 100
     fhistTimeConst = trim(dir_restart_hist)//trim(casename)//'-'//'rstTimeConst'
     open(unit=lhistTimeConst,file=trim(fhistTimeConst),status='unknown',&
                              form='unformatted',action='write')


     write(lhistTimeConst)  &!
           patch2lon,       &! longitude index for each patch point
           patch2lat,       &! latitude index for each patch point
           patchclass,      &! index of land cover type of the patches at the fraction > 0
           patchtype,       &! land water TYPE
           grid_patch_s,    &! start patch number of grid
           grid_patch_e,    &! end patch number of grid
           gridarea,        &! area of grid the patch located in
           patchfrac         ! subgrid weight for each patch point

     write(lhistTimeConst)  &!
           patchlatr,       &! patch latitude in radians
           patchlonr,       &! patch longitude in radians
           gridlatd,        &! grid latitude in degrees
           gridlond,        &! grid longitude in degrees

     ! Soil and plant parameters OF CLM
           lakedepth,       &! lake depth
           dz_lake,         &! new lake scheme

           soil_s_v_alb,    &! albedo of visible of the saturated soil
           soil_d_v_alb,    &! albedo of visible of the dry soil
           soil_s_n_alb,    &! albedo of near infrared of the saturated soil
           soil_d_n_alb,    &! albedo of near infrared of the dry soil
           porsl,           &! fraction of soil that is voids [-]
           psi0,            &! minimum soil suction [mm] (NOTE: "-" valued)
           bsw,             &! clapp and hornbereger "b" parameter [-]
           hksati,          &! hydraulic conductivity at saturation [mm h2o/s]
           csol,            &! heat capacity of soil solids [J/(m3 K)]
           dksatu,          &! thermal conductivity of saturated soil [W/m-K]
           dkdry,           &! thermal conductivity for dry soil  [W/(m-K)]

           htop,            &! canopy top height [m]
           hbot,            &! canopy bottom height [m]

     ! CLM TUNABLE constants
           zlnd,            &! roughness length for soil [m]
           zsno,            &! roughness length for snow [m]
           csoilc,          &! drag coefficient for soil under canopy [-]
           dewmx,           &! maximum dew
           wtfact,          &! fraction of model area with high water table
           capr,            &! tuning factor to turn first layer T into surface T
           cnfac,           &! Crank Nicholson factor between 0 and 1
           ssi,             &! irreducible water saturation of snow
           wimp,            &! water impremeable if porosity less than wimp
           pondmx,          &! ponding depth (mm)
           smpmax,          &! wilting point potential in mm
           smpmin,          &! restriction for min of soil poten. (mm)
           trsmx0,          &! max transpiration for moist soil+100% veg.  [mm/s]
           tcrit,           &! critical temp. to determine rain or snow

! bgc variables
           donor_pool     , &
           receiver_pool  , &
           floating_cn_ratio,& ! TRUE => pool has fixed C:N ratio
           initial_cn_ratio, & ! c:n ratio for initialization of pools
           rf_decomp       , & ! (frac) respired fraction in decomposition step
           pathfrac_decomp , & ! what fraction of C leaving a given pool passes through a given transition (frac)

           i_met_lit      , &
           i_cel_lit      , &
           i_lig_lit      , &
           i_cwd          , &
           i_soil1        , &
           i_soil2        , &
           i_soil3        , &
           i_atm          , &
           is_cwd         , &
           is_litter      , &
           is_soil        , &

           gdp_lf         , &
           abm_lf         , &
           peatf_lf       , &
           cmb_cmplt_fact , &

           nitrif_n2o_loss_frac, &! fraction of N lost as N2O in nitrification (Li et al., 2000)
           dnp                 , &!
           bdnr                , &!
           Q10                 , &!
           froz_q10            , &!
           tau_l1              , &!
           tau_l2_l3           , &!
           tau_s1              , &!
           tau_s2              , &!
           tau_s3              , &!
           tau_cwd             , &
           lwtop               , &

           som_adv_flux        , &
           som_diffus          , &
           cryoturb_diffusion_k, &
           max_altdepth_cryoturbation, &
           max_depth_cryoturb  , &

           am                  , &
           br                  , &
           br_root             , &
           fstor2tran          , &
           ndays_on            , &
           ndays_off           , &
           crit_dayl           , &
           crit_onset_fdd      , &
           crit_onset_swi      , &
           crit_offset_fdd     , &
           crit_offset_swi     , &
           soilpsi_on          , &
           soilpsi_off         , &

           occur_hi_gdp_tree   , &
           lfuel               , &
           ufuel               , &
           cropfire_a1         , &
           borealat            , &
           troplat             , &
           non_boreal_peatfire_c, &
           boreal_peatfire_c   , &
           rh_low              , &
           rh_hgh              , &
           bt_min              , &
           bt_max              , &
           pot_hmn_ign_counts_alpha, &
           g0, &

           sf, &
           sf_no3

     ! PFT/PC time invariants
#ifdef PFT_CLASSIFICATION
     write(lhistTimeConst)  &!
           pftclass,        &! PFT type
           pftfrac,         &! PFT fractional cover
           patch_pft_s,     &! patch start index of PFT
           patch_pft_e,     &! patch end index of PFT
           pft2patch,       &! projection from PFT to patch
           htop_p,          &! canopy top height [m]
           hbot_p            ! canopy bottom height [m]
#endif

#ifdef PC_CLASSIFICATION
     write(lhistTimeConst)  &!
           patch2pc,        &! projection from patch to PC
           pc2patch,        &! projection from PC to patch
           pcfrac,          &! PC fractional cover
           htop_c,          &! canopy top height [m]
           hbot_c            ! canopy bottom height [m]   
#endif
     
     close(lhistTimeConst)

  END SUBROUTINE WRITE_TimeInvariants


  SUBROUTINE deallocate_TimeInvariants
! --------------------------------------------------
! Deallocates memory for CLM 1d [numpatch] variables
! --------------------------------------------------
     USE MOD_PFTimeInvars
     USE MOD_PCTimeInvars

     deallocate (patch2lon    )
     deallocate (patch2lat    )
     deallocate (patchclass   )
     deallocate (patchtype    )
     deallocate (grid_patch_s )
     deallocate (grid_patch_e )
     deallocate (gridarea     )
     deallocate (patchfrac    )

     deallocate (patchlatr    )
     deallocate (patchlonr    )
     deallocate (gridlatd     )
     deallocate (gridlond     )

     deallocate (lakedepth    )
     deallocate (dz_lake      )

     deallocate (soil_s_v_alb )
     deallocate (soil_d_v_alb )
     deallocate (soil_s_n_alb )
     deallocate (soil_d_n_alb )
     deallocate (porsl        )
     deallocate (psi0         )
     deallocate (bsw          )
     deallocate (hksati       )
     deallocate (csol         )
     deallocate (dksatu       )
     deallocate (dkdry        )

     deallocate (htop         )
     deallocate (hbot         )

! bgc variables
     deallocate (rf_decomp      )
     deallocate (pathfrac_decomp)
     deallocate (gdp_lf         )
     deallocate (abm_lf         )
     deallocate (peatf_lf       )
  
#ifdef PFT_CLASSIFICATION
     CALL deallocate_PFTimeInvars
#endif

#ifdef PC_CLASSIFICATION
     CALL deallocate_PCTimeInvars
#endif

  END SUBROUTINE deallocate_TimeInvariants

END MODULE MOD_TimeInvariants
! ---------- EOP ------------
