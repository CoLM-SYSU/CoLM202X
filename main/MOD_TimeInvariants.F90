#include <define.h>

MODULE MOD_TimeInvariants 
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

use precision
IMPLICIT NONE
SAVE
! -----------------------------------------------------------------
! surface classification and soil information
  INTEGER,  allocatable :: patchclass     (:)  !index of land cover type of the patches at the fraction > 0
  INTEGER,  allocatable :: patchtype      (:)  !land water type

  REAL(r8), allocatable :: patchlatr      (:)  !latitude in radians
  REAL(r8), allocatable :: patchlonr      (:)  !longitude in radians

  REAL(r8), allocatable :: lakedepth      (:)  !lake depth
  REAL(r8), allocatable :: dz_lake      (:,:)  !new lake scheme

  REAL(r8), allocatable :: soil_s_v_alb   (:)  !albedo of visible of the saturated soil
  REAL(r8), allocatable :: soil_d_v_alb   (:)  !albedo of visible of the dry soil
  REAL(r8), allocatable :: soil_s_n_alb   (:)  !albedo of near infrared of the saturated soil
  REAL(r8), allocatable :: soil_d_n_alb   (:)  !albedo of near infrared of the dry soil

  REAL(r8), allocatable :: vf_quartz (:,:)     ! volumetric fraction of quartz within mineral soil
  REAL(r8), allocatable :: vf_gravels(:,:)     ! volumetric fraction of gravels
  REAL(r8), allocatable :: vf_om     (:,:)     ! volumetric fraction of organic matter
  REAL(r8), allocatable :: vf_sand   (:,:)     ! volumetric fraction of sand
  REAL(r8), allocatable :: wf_gravels(:,:)     ! gravimetric fraction of gravels
  REAL(r8), allocatable :: wf_sand   (:,:)     ! gravimetric fraction of sand

  REAL(r8), allocatable :: porsl        (:,:)  !fraction of soil that is voids [-]
  REAL(r8), allocatable :: psi0         (:,:)  !minimum soil suction [mm] (NOTE: "-" valued)
#ifdef Campbell_SOIL_MODEL
  REAL(r8), allocatable :: bsw          (:,:)  !clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
  REAL(r8), allocatable :: theta_r      (:,:)
  REAL(r8), allocatable :: alpha_vgm    (:,:)
  REAL(r8), allocatable :: L_vgm        (:,:)
  REAL(r8), allocatable :: n_vgm        (:,:)
  REAL(r8), allocatable :: sc_vgm       (:,:)
  REAL(r8), allocatable :: fc_vgm       (:,:)
#endif
  REAL(r8), allocatable :: hksati       (:,:)  !hydraulic conductivity at saturation [mm h2o/s]
  REAL(r8), allocatable :: csol         (:,:)  !heat capacity of soil solids [J/(m3 K)]
  REAL(r8), allocatable :: k_solids     (:,:)  !thermal conductivity of soil solids [W/m-K]
  REAL(r8), allocatable :: dksatu       (:,:)  !thermal conductivity of saturated soil [W/m-K]
  real(r8), allocatable :: dksatf       (:,:)  !thermal conductivity of saturated frozen soil [W/m-K]
  REAL(r8), allocatable :: dkdry        (:,:)  !thermal conductivity for dry soil  [W/(m-K)]
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
  REAL(r8), allocatable :: BA_alpha     (:,:)  !alpha in Balland and Arp(2005) thermal conductivity scheme
  REAL(r8), allocatable :: BA_beta      (:,:)  !beta in Balland and Arp(2005) thermal conductivity scheme
#endif
  REAL(r8), allocatable :: htop           (:)  !canopy top height [m]
  REAL(r8), allocatable :: hbot           (:)  !canopy bottom height [m]

#ifdef USE_DEPTH_TO_BEDROCK
  real(r8), allocatable :: dbedrock       (:)  ! depth to bedrock
  integer , allocatable :: ibedrock       (:)  ! bedrock level
#endif


  REAL(r8) :: zlnd         ! roughness length for soil [m]
  REAL(r8) :: zsno         ! roughness length for snow [m]
  REAL(r8) :: csoilc       ! drag coefficient for soil under canopy [-]
  REAL(r8) :: dewmx        ! maximum dew
  REAL(r8) :: wtfact       ! fraction of model area with high water table
  REAL(r8) :: capr         ! tuning factor to turn first layer T into surface T
  REAL(r8) :: cnfac        ! Crank Nicholson factor between 0 and 1
  REAL(r8) :: ssi          ! irreducible water saturation of snow
  REAL(r8) :: wimp         ! water impremeable if porosity less than wimp
  REAL(r8) :: pondmx       ! ponding depth (mm)
  REAL(r8) :: smpmax       ! wilting point potential in mm
  REAL(r8) :: smpmin       ! restriction for min of soil poten. (mm)
  REAL(r8) :: trsmx0       ! max transpiration for moist soil+100% veg.  [mm/s]
  REAL(r8) :: tcrit        ! critical temp. to determine rain or snow

!------------------------- BGC constant --------------------------------------
  INTEGER , allocatable  :: donor_pool       (:)
  INTEGER , allocatable  :: receiver_pool    (:)
  REAL(r8) :: am
  LOGICAL , allocatable  :: floating_cn_ratio(:)
  REAL(r8), allocatable :: initial_cn_ratio (:)
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

  LOGICAL , allocatable :: is_cwd    (:)!(1:ndecomp_pools) ! True => is a coarse woody debris pool
  LOGICAL , allocatable :: is_litter (:)!(1:ndecomp_pools) ! True => is a litter pool
  LOGICAL , allocatable :: is_soil   (:)!(1:ndecomp_pools) ! True => is a soil pool

  REAL(r8), allocatable :: gdp_lf (:) !
  REAL(r8), allocatable :: abm_lf (:) !
  REAL(r8), allocatable :: peatf_lf(:)!
  REAL(r8), allocatable :: cmb_cmplt_fact(:)

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

!----------------------------------- end BGC constants -----------------

          
! PUBLIC MEMBER FUNCTIONS:
  public :: allocate_TimeInvariants
  public :: deallocate_TimeInvariants
  public :: READ_TimeInvariants
  public :: WRITE_TimeInvariants

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_TimeInvariants ()
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM 1d [numpatch] variables
  ! --------------------------------------------------------------------

     use precision
     use GlobalVars, only: nl_soil, ndecomp_transitions, ndecomp_pools
     use spmd_task
     use mod_landpatch, only : numpatch
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
#endif
     IMPLICIT NONE

  if (p_is_worker) then

     if (numpatch > 0) then

        allocate (patchclass           (numpatch))
        allocate (patchtype            (numpatch))
        
        allocate (patchlonr            (numpatch))
        allocate (patchlatr            (numpatch))

        allocate (lakedepth            (numpatch))
        allocate (dz_lake      (nl_lake,numpatch))

        allocate (soil_s_v_alb         (numpatch))
        allocate (soil_d_v_alb         (numpatch))
        allocate (soil_s_n_alb         (numpatch))
        allocate (soil_d_n_alb         (numpatch))

        allocate (vf_quartz    (nl_soil,numpatch))
        allocate (vf_gravels   (nl_soil,numpatch))
        allocate (vf_om        (nl_soil,numpatch))
        allocate (vf_sand      (nl_soil,numpatch))
        allocate (wf_gravels   (nl_soil,numpatch))
        allocate (wf_sand      (nl_soil,numpatch))
        allocate (porsl        (nl_soil,numpatch))
        allocate (psi0         (nl_soil,numpatch))
#ifdef Campbell_SOIL_MODEL
        allocate (bsw          (nl_soil,numpatch))
#endif

#ifdef vanGenuchten_Mualem_SOIL_MODEL
        allocate (theta_r      (nl_soil,numpatch))
        allocate (alpha_vgm    (nl_soil,numpatch))
        allocate (L_vgm        (nl_soil,numpatch))
        allocate (n_vgm        (nl_soil,numpatch))
        allocate (sc_vgm       (nl_soil,numpatch))
        allocate (fc_vgm       (nl_soil,numpatch))
#endif
        allocate (hksati       (nl_soil,numpatch))
        allocate (csol         (nl_soil,numpatch))
        allocate (k_solids     (nl_soil,numpatch))
        allocate (dksatu       (nl_soil,numpatch))
        allocate (dksatf       (nl_soil,numpatch))
        allocate (dkdry        (nl_soil,numpatch))
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4     
        allocate (BA_alpha     (nl_soil,numpatch))
        allocate (BA_beta      (nl_soil,numpatch))
#endif
        allocate (htop                 (numpatch))
        allocate (hbot                 (numpatch))

#ifdef USE_DEPTH_TO_BEDROCK
        allocate (dbedrock             (numpatch))
        allocate (ibedrock             (numpatch))
#endif
     end if

#ifdef BGC
! bgc varaibles
     allocate (donor_pool        (ndecomp_transitions))
     allocate (receiver_pool     (ndecomp_transitions))
     allocate (floating_cn_ratio (ndecomp_pools))
     allocate (initial_cn_ratio  (ndecomp_pools))
     allocate (rf_decomp         (nl_soil,ndecomp_transitions,numpatch))
     allocate (pathfrac_decomp   (nl_soil,ndecomp_transitions,numpatch))
     allocate (is_cwd            (ndecomp_pools)) ! True => is a coarse woody debris pool
     allocate (is_litter         (ndecomp_pools)) ! True => is a litter pool
     allocate (is_soil           (ndecomp_pools)) ! True => is a soil pool
     allocate (gdp_lf            (numpatch))
     allocate (abm_lf            (numpatch))
     allocate (peatf_lf          (numpatch))
     allocate (cmb_cmplt_fact    (2))

! end bgc variables
#endif

#ifdef PFT_CLASSIFICATION
     CALL allocate_PFTimeInvars
#endif

#ifdef PC_CLASSIFICATION
     CALL allocate_PCTimeInvars
#endif

  end if

  END SUBROUTINE allocate_TimeInvariants

  !---------------------------------------
  SUBROUTINE READ_TimeInvariants (casename, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use mod_namelist
     use spmd_task
     use ncio_vector
     use ncio_serial
#ifdef CLMDEBUG 
     USE mod_colm_debug
#endif
     USE mod_landpatch
     USE GlobalVars
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
#endif

     IMPLICIT NONE

     character(LEN=*), intent(in) :: casename
     character(LEN=*), intent(in) :: dir_restart
  
     ! Local variables
     character(LEN=256) :: file_restart

     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_const.nc'

     call ncio_read_vector (file_restart, 'patchclass',   landpatch, patchclass) !
     call ncio_read_vector (file_restart, 'patchtype' ,   landpatch, patchtype ) !
     
     call ncio_read_vector (file_restart, 'patchlonr' ,   landpatch, patchlonr ) !
     call ncio_read_vector (file_restart, 'patchlatr' ,   landpatch, patchlatr ) !
     
     call ncio_read_vector (file_restart, 'lakedepth',    landpatch, lakedepth) !
     call ncio_read_vector (file_restart, 'dz_lake' ,     nl_lake, landpatch, dz_lake) !
     
     call ncio_read_vector (file_restart, 'soil_s_v_alb', landpatch, soil_s_v_alb) ! albedo of visible of the saturated soil
     call ncio_read_vector (file_restart, 'soil_d_v_alb', landpatch, soil_d_v_alb) ! albedo of visible of the dry soil
     call ncio_read_vector (file_restart, 'soil_s_n_alb', landpatch, soil_s_n_alb) ! albedo of near infrared of the saturated soil
     call ncio_read_vector (file_restart, 'soil_d_n_alb', landpatch, soil_d_n_alb) ! albedo of near infrared of the dry soil

     call ncio_read_vector (file_restart, 'vf_quartz ',   nl_soil, landpatch, vf_quartz ) ! volumetric fraction of quartz within mineral soil
     call ncio_read_vector (file_restart, 'vf_gravels',   nl_soil, landpatch, vf_gravels) ! volumetric fraction of gravels
     call ncio_read_vector (file_restart, 'vf_om     ',   nl_soil, landpatch, vf_om     ) ! volumetric fraction of organic matter
     call ncio_read_vector (file_restart, 'vf_sand   ',   nl_soil, landpatch, vf_sand   ) ! volumetric fraction of sand
     call ncio_read_vector (file_restart, 'wf_gravels',   nl_soil, landpatch, wf_gravels) ! gravimetric fraction of gravels
     call ncio_read_vector (file_restart, 'wf_sand   ',   nl_soil, landpatch, wf_sand   ) ! gravimetric fraction of sand
     call ncio_read_vector (file_restart, 'porsl  ' ,     nl_soil, landpatch, porsl     ) ! fraction of soil that is voids [-]
     call ncio_read_vector (file_restart, 'psi0   ' ,     nl_soil, landpatch, psi0      ) ! minimum soil suction [mm] (NOTE: "-" valued)
#ifdef Campbell_SOIL_MODEL                                         
     call ncio_read_vector (file_restart, 'bsw    ' ,     nl_soil, landpatch, bsw       ) ! clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
     call ncio_read_vector (file_restart, 'theta_r  ' ,   nl_soil, landpatch, theta_r   ) 
     call ncio_read_vector (file_restart, 'alpha_vgm' ,   nl_soil, landpatch, alpha_vgm ) 
     call ncio_read_vector (file_restart, 'L_vgm    ' ,   nl_soil, landpatch, L_vgm     ) 
     call ncio_read_vector (file_restart, 'n_vgm    ' ,   nl_soil, landpatch, n_vgm     ) 
     call ncio_read_vector (file_restart, 'sc_vgm   ' ,   nl_soil, landpatch, sc_vgm    ) 
     call ncio_read_vector (file_restart, 'fc_vgm   ' ,   nl_soil, landpatch, fc_vgm    ) 
#endif                                                             
     call ncio_read_vector (file_restart, 'hksati ' ,     nl_soil, landpatch, hksati ) ! hydraulic conductivity at saturation [mm h2o/s]
     call ncio_read_vector (file_restart, 'csol   ' ,     nl_soil, landpatch, csol   ) ! heat capacity of soil solids [J/(m3 K)]
     call ncio_read_vector (file_restart, 'k_solids',     nl_soil, landpatch, k_solids)! thermal conductivity of soil solids [W/m-K] 
     call ncio_read_vector (file_restart, 'dksatu ' ,     nl_soil, landpatch, dksatu ) ! thermal conductivity of unfrozen saturated soil [W/m-K]
     call ncio_read_vector (file_restart, 'dksatf ' ,     nl_soil, landpatch, dksatf ) ! thermal conductivity of frozen saturated soil [W/m-K]
     call ncio_read_vector (file_restart, 'dkdry  ' ,     nl_soil, landpatch, dkdry  ) ! thermal conductivity for dry soil  [W/(m-K)]
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4     
     call ncio_read_vector (file_restart, 'BA_alpha',     nl_soil, landpatch, BA_alpha)! alpha in Balland and Arp(2005) thermal conductivity scheme
     call ncio_read_vector (file_restart, 'BA_beta' ,     nl_soil, landpatch, BA_beta )! beta in Balland and Arp(2005) thermal conductivity scheme
#endif
     call ncio_read_vector (file_restart, 'htop' ,    landpatch, htop) !
     call ncio_read_vector (file_restart, 'hbot' ,    landpatch, hbot) !

#ifdef USE_DEPTH_TO_BEDROCK
     call ncio_read_vector (file_restart, 'debdrock' ,    landpatch, dbedrock) !
     call ncio_read_vector (file_restart, 'ibedrock' ,    landpatch, ibedrock) !
#endif

     call ncio_read_bcast_serial (file_restart, 'zlnd  ', zlnd  ) ! roughness length for soil [m]
     call ncio_read_bcast_serial (file_restart, 'zsno  ', zsno  ) ! roughness length for snow [m]
     call ncio_read_bcast_serial (file_restart, 'csoilc', csoilc) ! drag coefficient for soil under canopy [-]
     call ncio_read_bcast_serial (file_restart, 'dewmx ', dewmx ) ! maximum dew
     call ncio_read_bcast_serial (file_restart, 'wtfact', wtfact) ! fraction of model area with high water table
     call ncio_read_bcast_serial (file_restart, 'capr  ', capr  ) ! tuning factor to turn first layer T into surface T
     call ncio_read_bcast_serial (file_restart, 'cnfac ', cnfac ) ! Crank Nicholson factor between 0 and 1
     call ncio_read_bcast_serial (file_restart, 'ssi   ', ssi   ) ! irreducible water saturation of snow
     call ncio_read_bcast_serial (file_restart, 'wimp  ', wimp  ) ! water impremeable if porosity less than wimp
     call ncio_read_bcast_serial (file_restart, 'pondmx', pondmx) ! ponding depth (mm)
     call ncio_read_bcast_serial (file_restart, 'smpmax', smpmax) ! wilting point potential in mm
     call ncio_read_bcast_serial (file_restart, 'smpmin', smpmin) ! restriction for min of soil poten. (mm)
     call ncio_read_bcast_serial (file_restart, 'trsmx0', trsmx0) ! max transpiration for moist soil+100% veg.  [mm/s]
     call ncio_read_bcast_serial (file_restart, 'tcrit ', tcrit ) ! critical temp. to determine rain or snow

#ifdef BGC
! bgc constants
     call ncio_read_bcast_serial (file_restart, 'donor_pool     ', donor_pool     )
     call ncio_read_bcast_serial (file_restart, 'receiver_pool  ', receiver_pool  )
     call ncio_read_bcast_serial (file_restart, 'floating_cn_ratio', floating_cn_ratio)
     call ncio_read_bcast_serial (file_restart, 'initial_cn_ratio' , initial_cn_ratio) 
     call ncio_read_vector       (file_restart, 'rf_decomp       ', nl_soil,ndecomp_transitions,landpatch, rf_decomp      )
     call ncio_read_vector       (file_restart, 'pathfrac_decomp ', nl_soil,ndecomp_transitions,landpatch,pathfrac_decomp )

     call ncio_read_bcast_serial (file_restart, 'i_met_lit      ', i_met_lit      )
     call ncio_read_bcast_serial (file_restart, 'i_cel_lit      ', i_cel_lit      )
     call ncio_read_bcast_serial (file_restart, 'i_lig_lit      ', i_lig_lit      )
     call ncio_read_bcast_serial (file_restart, 'i_cwd          ', i_cwd          )
     call ncio_read_bcast_serial (file_restart, 'i_soil1        ', i_soil1        )
     call ncio_read_bcast_serial (file_restart, 'i_soil2        ', i_soil2        )
     call ncio_read_bcast_serial (file_restart, 'i_soil3        ', i_soil3        )
     call ncio_read_bcast_serial (file_restart, 'i_atm          ', i_atm          )
     call ncio_read_bcast_serial (file_restart, 'is_cwd         ', is_cwd         )
     call ncio_read_bcast_serial (file_restart, 'is_litter      ', is_litter      )
     call ncio_read_bcast_serial (file_restart, 'is_soil        ', is_soil        )

     call ncio_read_vector       (file_restart, 'gdp_lf         ', landpatch, gdp_lf        )
     call ncio_read_vector       (file_restart, 'abm_lf         ', landpatch, abm_lf        )
     call ncio_read_vector       (file_restart, 'peatf_lf       ', landpatch, peatf_lf      )
     call ncio_read_bcast_serial (file_restart, 'cmb_cmplt_fact ', cmb_cmplt_fact )

     call ncio_read_bcast_serial (file_restart, 'nitrif_n2o_loss_frac', nitrif_n2o_loss_frac)
     call ncio_read_bcast_serial (file_restart, 'dnp                 ', dnp                 )!
     call ncio_read_bcast_serial (file_restart, 'bdnr                ', bdnr                )!
     call ncio_read_bcast_serial (file_restart, 'Q10                 ', Q10                 )!
     call ncio_read_bcast_serial (file_restart, 'froz_q10            ', froz_q10            )!
     call ncio_read_bcast_serial (file_restart, 'tau_l1              ', tau_l1              )!
     call ncio_read_bcast_serial (file_restart, 'tau_l2_l3           ', tau_l2_l3           )!
     call ncio_read_bcast_serial (file_restart, 'tau_s1              ', tau_s1              )!
     call ncio_read_bcast_serial (file_restart, 'tau_s2              ', tau_s2              )!
     call ncio_read_bcast_serial (file_restart, 'tau_s3              ', tau_s3              )!
     call ncio_read_bcast_serial (file_restart, 'tau_cwd             ', tau_cwd             )
     call ncio_read_bcast_serial (file_restart, 'lwtop               ', lwtop               )

     call ncio_read_bcast_serial (file_restart, 'som_adv_flux        ', som_adv_flux        )
     call ncio_read_bcast_serial (file_restart, 'som_diffus          ', som_diffus          )
     call ncio_read_bcast_serial (file_restart, 'cryoturb_diffusion_k', cryoturb_diffusion_k)
     call ncio_read_bcast_serial (file_restart, 'max_altdepth_cryoturbation', max_altdepth_cryoturbation)
     call ncio_read_bcast_serial (file_restart, 'max_depth_cryoturb  ', max_depth_cryoturb  )

     call ncio_read_bcast_serial (file_restart, 'am                  ', am                  )
     call ncio_read_bcast_serial (file_restart, 'br                  ', br                  )
     call ncio_read_bcast_serial (file_restart, 'br_root             ', br_root             )
     call ncio_read_bcast_serial (file_restart, 'fstor2tran          ', fstor2tran          )
     call ncio_read_bcast_serial (file_restart, 'ndays_on            ', ndays_on            )
     call ncio_read_bcast_serial (file_restart, 'ndays_off           ', ndays_off           )
     call ncio_read_bcast_serial (file_restart, 'crit_dayl           ', crit_dayl           )
     call ncio_read_bcast_serial (file_restart, 'crit_onset_fdd      ', crit_onset_fdd      )
     call ncio_read_bcast_serial (file_restart, 'crit_onset_swi      ', crit_onset_swi      )
     call ncio_read_bcast_serial (file_restart, 'crit_offset_fdd     ', crit_offset_fdd     )
     call ncio_read_bcast_serial (file_restart, 'crit_offset_swi     ', crit_offset_swi     )
     call ncio_read_bcast_serial (file_restart, 'soilpsi_on          ', soilpsi_on          )
     call ncio_read_bcast_serial (file_restart, 'soilpsi_off         ', soilpsi_off         )

     call ncio_read_bcast_serial (file_restart, 'occur_hi_gdp_tree   ', occur_hi_gdp_tree   )
     call ncio_read_bcast_serial (file_restart, 'lfuel               ', lfuel               )
     call ncio_read_bcast_serial (file_restart, 'ufuel               ', ufuel               )
     call ncio_read_bcast_serial (file_restart, 'cropfire_a1         ', cropfire_a1         )
     call ncio_read_bcast_serial (file_restart, 'borealat            ', borealat            )
     call ncio_read_bcast_serial (file_restart, 'troplat             ', troplat             )
     call ncio_read_bcast_serial (file_restart, 'non_boreal_peatfire_c', non_boreal_peatfire_c)
     call ncio_read_bcast_serial (file_restart, 'boreal_peatfire_c   ', boreal_peatfire_c   )
     call ncio_read_bcast_serial (file_restart, 'rh_low              ', rh_low              )
     call ncio_read_bcast_serial (file_restart, 'rh_hgh              ', rh_hgh              )
     call ncio_read_bcast_serial (file_restart, 'bt_min              ', bt_min              )
     call ncio_read_bcast_serial (file_restart, 'bt_max              ', bt_max              )
     call ncio_read_bcast_serial (file_restart, 'pot_hmn_ign_counts_alpha', pot_hmn_ign_counts_alpha)
     call ncio_read_bcast_serial (file_restart, 'g0', g0)

     call ncio_read_bcast_serial (file_restart, 'sf', sf)
     call ncio_read_bcast_serial (file_restart, 'sf_no3', sf_no3)
#endif

#if (defined PFT_CLASSIFICATION)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_pft_const.nc'
     CALL READ_PFTimeInvars (file_restart)
#endif 

#if (defined PC_CLASSIFICATION)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_pc_const.nc'
     CALL READ_PCTimeInvars (file_restart)
#endif 

#ifdef CLMDEBUG 
     call check_TimeInvariants ()
#endif

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

     if (p_is_master) then
        write(*,'(A29)') 'Loading Time Invariants done.'
     end if

  end subroutine READ_TimeInvariants

  !---------------------------------------
  SUBROUTINE WRITE_TimeInvariants (casename, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use mod_namelist, only : DEF_REST_COMPRESS_LEVEL 
     use spmd_task
     use ncio_serial
     use ncio_vector
     use mod_landpatch
     USE GlobalVars
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
#endif

     IMPLICIT NONE

     character(len=*), intent(in) :: casename
     character(len=*), intent(in) :: dir_restart
     
     ! Local Variables
     character(len=256) :: file_restart
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL 

     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_const.nc'

     call ncio_create_file_vector (file_restart, landpatch)

     CALL ncio_define_pixelset_dimension (file_restart, landpatch)
     CALL ncio_define_dimension_vector (file_restart, 'soil', nl_soil)
     CALL ncio_define_dimension_vector (file_restart, 'lake', nl_lake)
     CALL ncio_define_dimension_vector (file_restart, 'band',   2)
     CALL ncio_define_dimension_vector (file_restart, 'wetdry', 2)
     CALL ncio_define_dimension_vector (file_restart, 'ndecomp_transitions',ndecomp_transitions)
     CALL ncio_define_dimension_vector (file_restart, 'ndecomp_pools'      ,ndecomp_pools)
     
     call ncio_write_vector (file_restart, 'patchclass', 'vector', landpatch, patchclass) !
     call ncio_write_vector (file_restart, 'patchtype' , 'vector', landpatch, patchtype ) !
     
     call ncio_write_vector (file_restart, 'patchlonr' , 'vector', landpatch, patchlonr ) !
     call ncio_write_vector (file_restart, 'patchlatr' , 'vector', landpatch, patchlatr ) !
     
     call ncio_write_vector (file_restart, 'lakedepth', 'vector', landpatch, lakedepth ,      compress) !
     call ncio_write_vector (file_restart, 'dz_lake' ,  'lake', nl_lake, 'vector', landpatch, dz_lake, compress) !
     
     call ncio_write_vector (file_restart, 'soil_s_v_alb', 'vector', landpatch, soil_s_v_alb, compress) ! albedo of visible of the saturated soil
     call ncio_write_vector (file_restart, 'soil_d_v_alb', 'vector', landpatch, soil_d_v_alb, compress) ! albedo of visible of the dry soil
     call ncio_write_vector (file_restart, 'soil_s_n_alb', 'vector', landpatch, soil_s_n_alb, compress) ! albedo of near infrared of the saturated soil
     call ncio_write_vector (file_restart, 'soil_d_n_alb', 'vector', landpatch, soil_d_n_alb, compress) ! albedo of near infrared of the dry soil

     call ncio_write_vector (file_restart, 'vf_quartz ', 'soil', nl_soil, 'vector', landpatch, vf_quartz , compress) ! volumetric fraction of quartz within mineral soil
     call ncio_write_vector (file_restart, 'vf_gravels', 'soil', nl_soil, 'vector', landpatch, vf_gravels, compress) ! volumetric fraction of gravels
     call ncio_write_vector (file_restart, 'vf_om     ', 'soil', nl_soil, 'vector', landpatch, vf_om     , compress) ! volumetric fraction of organic matter
     call ncio_write_vector (file_restart, 'vf_sand   ', 'soil', nl_soil, 'vector', landpatch, vf_sand   , compress) ! volumetric fraction of sand
     call ncio_write_vector (file_restart, 'wf_gravels', 'soil', nl_soil, 'vector', landpatch, wf_gravels, compress) ! gravimetric fraction of gravels
     call ncio_write_vector (file_restart, 'wf_sand   ', 'soil', nl_soil, 'vector', landpatch, wf_sand   , compress) ! gravimetric fraction of sand
     call ncio_write_vector (file_restart, 'porsl     ', 'soil', nl_soil, 'vector', landpatch, porsl     , compress) ! fraction of soil that is voids [-]
     call ncio_write_vector (file_restart, 'psi0      ', 'soil', nl_soil, 'vector', landpatch, psi0      , compress) ! minimum soil suction [mm] (NOTE: "-" valued)
#ifdef Campbell_SOIL_MODEL
     call ncio_write_vector (file_restart, 'bsw       ', 'soil', nl_soil, 'vector', landpatch, bsw       , compress) ! clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
     call ncio_write_vector (file_restart, 'theta_r  ' , 'soil', nl_soil, 'vector', landpatch, theta_r   , compress) 
     call ncio_write_vector (file_restart, 'alpha_vgm' , 'soil', nl_soil, 'vector', landpatch, alpha_vgm , compress) 
     call ncio_write_vector (file_restart, 'L_vgm    ' , 'soil', nl_soil, 'vector', landpatch, L_vgm     , compress) 
     call ncio_write_vector (file_restart, 'n_vgm    ' , 'soil', nl_soil, 'vector', landpatch, n_vgm     , compress) 
     call ncio_write_vector (file_restart, 'sc_vgm   ' , 'soil', nl_soil, 'vector', landpatch, sc_vgm    , compress) 
     call ncio_write_vector (file_restart, 'fc_vgm   ' , 'soil', nl_soil, 'vector', landpatch, fc_vgm    , compress) 
#endif
     call ncio_write_vector (file_restart, 'hksati   ' , 'soil', nl_soil, 'vector', landpatch, hksati    , compress) ! hydraulic conductivity at saturation [mm h2o/s]
     call ncio_write_vector (file_restart, 'csol     ' , 'soil', nl_soil, 'vector', landpatch, csol      , compress) ! heat capacity of soil solids [J/(m3 K)]
     call ncio_write_vector (file_restart, 'k_solids ' , 'soil', nl_soil, 'vector', landpatch, k_solids  , compress) ! thermal conductivity of soil solids [W/m-K]
     call ncio_write_vector (file_restart, 'dksatu   ' , 'soil', nl_soil, 'vector', landpatch, dksatu    , compress) ! thermal conductivity of saturated soil [W/m-K]
     call ncio_write_vector (file_restart, 'dksatf   ' , 'soil', nl_soil, 'vector', landpatch, dksatf    , compress) ! thermal conductivity of saturated soil [W/m-K]
     call ncio_write_vector (file_restart, 'dkdry    ' , 'soil', nl_soil, 'vector', landpatch, dkdry     , compress) ! thermal conductivity for dry soil  [W/(m-K)]
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4     
     call ncio_write_vector (file_restart, 'BA_alpha ' , 'soil', nl_soil, 'vector', landpatch, BA_alpha  , compress) ! alpha in Balland and Arp(2005) thermal conductivity scheme
     call ncio_write_vector (file_restart, 'BA_beta  ' , 'soil', nl_soil, 'vector', landpatch, BA_beta   , compress) ! beta in Balland and Arp(2005) thermal conductivity scheme
#endif

     call ncio_write_vector (file_restart, 'htop' , 'vector', landpatch, htop) !
     call ncio_write_vector (file_restart, 'hbot' , 'vector', landpatch, hbot) !

#ifdef USE_DEPTH_TO_BEDROCK
     call ncio_write_vector (file_restart, 'debdrock' , 'vector', landpatch, dbedrock) !
     call ncio_write_vector (file_restart, 'ibedrock' , 'vector', landpatch, ibedrock) !
#endif

#ifdef BGC
     call ncio_write_vector       (file_restart, 'rf_decomp      ', 'soil'   , nl_soil  , &
                          'ndecomp_transitions', ndecomp_transitions,'vector', landpatch, rf_decomp      , compress)
     call ncio_write_vector       (file_restart, 'pathfrac_decomp', 'soil'   , nl_soil  , &
                          'ndecomp_transitions', ndecomp_transitions,'vector', landpatch, pathfrac_decomp, compress)
     call ncio_write_vector       (file_restart, 'gdp_lf         ',  'vector', landpatch, gdp_lf         , compress)
     call ncio_write_vector       (file_restart, 'abm_lf         ',  'vector', landpatch, abm_lf         , compress)
     call ncio_write_vector       (file_restart, 'peatf_lf       ',  'vector', landpatch, peatf_lf       , compress)
#endif
     if (p_is_master) then

        call ncio_create_file (file_restart)
        call ncio_define_dimension(file_restart, 'ndecomp_transitions',ndecomp_transitions)
        call ncio_define_dimension(file_restart, 'ndecomp_pools'      ,ndecomp_pools)
        call ncio_define_dimension(file_restart, 'nlitter_fire'       ,2            )

        call ncio_write_serial (file_restart, 'zlnd  ', zlnd  ) ! roughness length for soil [m]
        call ncio_write_serial (file_restart, 'zsno  ', zsno  ) ! roughness length for snow [m]
        call ncio_write_serial (file_restart, 'csoilc', csoilc) ! drag coefficient for soil under canopy [-]
        call ncio_write_serial (file_restart, 'dewmx ', dewmx ) ! maximum dew
        call ncio_write_serial (file_restart, 'wtfact', wtfact) ! fraction of model area with high water table
        call ncio_write_serial (file_restart, 'capr  ', capr  ) ! tuning factor to turn first layer T into surface T
        call ncio_write_serial (file_restart, 'cnfac ', cnfac ) ! Crank Nicholson factor between 0 and 1
        call ncio_write_serial (file_restart, 'ssi   ', ssi   ) ! irreducible water saturation of snow
        call ncio_write_serial (file_restart, 'wimp  ', wimp  ) ! water impremeable if porosity less than wimp
        call ncio_write_serial (file_restart, 'pondmx', pondmx) ! ponding depth (mm)
        call ncio_write_serial (file_restart, 'smpmax', smpmax) ! wilting point potential in mm
        call ncio_write_serial (file_restart, 'smpmin', smpmin) ! restriction for min of soil poten. (mm)
        call ncio_write_serial (file_restart, 'trsmx0', trsmx0) ! max transpiration for moist soil+100% veg.  [mm/s]
        call ncio_write_serial (file_restart, 'tcrit ', tcrit ) ! critical temp. to determine rain or snow

#ifdef BGC
! bgc constants
        call ncio_write_serial (file_restart, 'donor_pool     '  , donor_pool       , 'ndecomp_transitions')
        call ncio_write_serial (file_restart, 'receiver_pool  '  , receiver_pool    , 'ndecomp_transitions')
        call ncio_write_serial (file_restart, 'floating_cn_ratio', floating_cn_ratio, 'ndecomp_pools')
        call ncio_write_serial (file_restart, 'initial_cn_ratio' , initial_cn_ratio , 'ndecomp_pools') 
        call ncio_write_serial (file_restart, 'is_cwd         '  , is_cwd           , 'ndecomp_pools')
        call ncio_write_serial (file_restart, 'is_litter      '  , is_litter        , 'ndecomp_pools')
        call ncio_write_serial (file_restart, 'is_soil        '  , is_soil          , 'ndecomp_pools')
        call ncio_write_serial (file_restart, 'cmb_cmplt_fact '  , cmb_cmplt_fact   , 'nlitter_fire' )

        call ncio_write_serial (file_restart, 'i_met_lit      ', i_met_lit      )
        call ncio_write_serial (file_restart, 'i_cel_lit      ', i_cel_lit      )
        call ncio_write_serial (file_restart, 'i_lig_lit      ', i_lig_lit      )
        call ncio_write_serial (file_restart, 'i_cwd          ', i_cwd          )
        call ncio_write_serial (file_restart, 'i_soil1        ', i_soil1        )
        call ncio_write_serial (file_restart, 'i_soil2        ', i_soil2        )
        call ncio_write_serial (file_restart, 'i_soil3        ', i_soil3        )
        call ncio_write_serial (file_restart, 'i_atm          ', i_atm          )


        call ncio_write_serial (file_restart, 'nitrif_n2o_loss_frac', nitrif_n2o_loss_frac)
        call ncio_write_serial (file_restart, 'dnp                 ', dnp                 )!
        call ncio_write_serial (file_restart, 'bdnr                ', bdnr                )!
        call ncio_write_serial (file_restart, 'Q10                 ', Q10                 )!
        call ncio_write_serial (file_restart, 'froz_q10            ', froz_q10            )!
        call ncio_write_serial (file_restart, 'tau_l1              ', tau_l1              )!
        call ncio_write_serial (file_restart, 'tau_l2_l3           ', tau_l2_l3           )!
        call ncio_write_serial (file_restart, 'tau_s1              ', tau_s1              )!
        call ncio_write_serial (file_restart, 'tau_s2              ', tau_s2              )!
        call ncio_write_serial (file_restart, 'tau_s3              ', tau_s3              )!
        call ncio_write_serial (file_restart, 'tau_cwd             ', tau_cwd             )
        call ncio_write_serial (file_restart, 'lwtop               ', lwtop               )
   
        call ncio_write_serial (file_restart, 'som_adv_flux        ', som_adv_flux        )
        call ncio_write_serial (file_restart, 'som_diffus          ', som_diffus          )
        call ncio_write_serial (file_restart, 'cryoturb_diffusion_k', cryoturb_diffusion_k)
        call ncio_write_serial (file_restart, 'max_altdepth_cryoturbation', max_altdepth_cryoturbation)
        call ncio_write_serial (file_restart, 'max_depth_cryoturb  ', max_depth_cryoturb  )

        call ncio_write_serial (file_restart, 'am                  ', am                  )
        call ncio_write_serial (file_restart, 'br                  ', br                  )
        call ncio_write_serial (file_restart, 'br_root             ', br_root             )
        call ncio_write_serial (file_restart, 'fstor2tran          ', fstor2tran          )
        call ncio_write_serial (file_restart, 'ndays_on            ', ndays_on            )
        call ncio_write_serial (file_restart, 'ndays_off           ', ndays_off           )
        call ncio_write_serial (file_restart, 'crit_dayl           ', crit_dayl           )
        call ncio_write_serial (file_restart, 'crit_onset_fdd      ', crit_onset_fdd      )
        call ncio_write_serial (file_restart, 'crit_onset_swi      ', crit_onset_swi      )
        call ncio_write_serial (file_restart, 'crit_offset_fdd     ', crit_offset_fdd     )
        call ncio_write_serial (file_restart, 'crit_offset_swi     ', crit_offset_swi     )
        call ncio_write_serial (file_restart, 'soilpsi_on          ', soilpsi_on          )
        call ncio_write_serial (file_restart, 'soilpsi_off         ', soilpsi_off         )

        call ncio_write_serial (file_restart, 'occur_hi_gdp_tree   ', occur_hi_gdp_tree   )
        call ncio_write_serial (file_restart, 'lfuel               ', lfuel               )
        call ncio_write_serial (file_restart, 'ufuel               ', ufuel               )
        call ncio_write_serial (file_restart, 'cropfire_a1         ', cropfire_a1         )
        call ncio_write_serial (file_restart, 'borealat            ', borealat            )
        call ncio_write_serial (file_restart, 'troplat             ', troplat             )
        call ncio_write_serial (file_restart, 'non_boreal_peatfire_c', non_boreal_peatfire_c)
        call ncio_write_serial (file_restart, 'boreal_peatfire_c   ', boreal_peatfire_c   )
        call ncio_write_serial (file_restart, 'rh_low              ', rh_low              )
        call ncio_write_serial (file_restart, 'rh_hgh              ', rh_hgh              )
        call ncio_write_serial (file_restart, 'bt_min              ', bt_min              )
        call ncio_write_serial (file_restart, 'bt_max              ', bt_max              )
        call ncio_write_serial (file_restart, 'pot_hmn_ign_counts_alpha', pot_hmn_ign_counts_alpha)
        call ncio_write_serial (file_restart, 'g0', g0)

        call ncio_write_serial (file_restart, 'sf', sf)
        call ncio_write_serial (file_restart, 'sf_no3', sf_no3)
#endif

     end if

#if (defined PFT_CLASSIFICATION)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_pft_const.nc'
     CALL WRITE_PFTimeInvars (file_restart)
#endif 

#if (defined PC_CLASSIFICATION)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_pc_const.nc'
     CALL WRITE_PCTimeInvars (file_restart)
#endif 

   end subroutine WRITE_TimeInvariants 

  SUBROUTINE deallocate_TimeInvariants ()

     use spmd_task
     use mod_landpatch, only : numpatch
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
#endif
     implicit none

     ! --------------------------------------------------
     ! Deallocates memory for CLM 1d [numpatch] variables
     ! --------------------------------------------------

     if (p_is_worker) then

        if (numpatch > 0) then

           deallocate (patchclass)
           deallocate (patchtype )
           
           deallocate (patchlonr )
           deallocate (patchlatr )

           deallocate (lakedepth)
           deallocate (dz_lake  )

           deallocate (soil_s_v_alb)
           deallocate (soil_d_v_alb)
           deallocate (soil_s_n_alb)
           deallocate (soil_d_n_alb)

           deallocate (vf_quartz )
           deallocate (vf_gravels)
           deallocate (vf_om     )
           deallocate (vf_sand   )
           deallocate (wf_gravels)
           deallocate (wf_sand   )
           deallocate (porsl  )
           deallocate (psi0   )
#ifdef Campbell_SOIL_MODEL
           deallocate (bsw    )
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
           deallocate (theta_r  )
           deallocate (alpha_vgm)
           deallocate (L_vgm    )
           deallocate (n_vgm    )
           deallocate (sc_vgm   )
           deallocate (fc_vgm   )
#endif
           deallocate (hksati )
           deallocate (csol   )
           deallocate (k_solids)
           deallocate (dksatu )
           deallocate (dksatf )
           deallocate (dkdry  )
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
           deallocate (BA_alpha )
           deallocate (BA_beta  )
#endif

           deallocate (htop)
           deallocate (hbot)
          
#ifdef USE_DEPTH_TO_BEDROCK
           deallocate (dbedrock)
           deallocate (ibedrock)
#endif

#ifdef BGC
! bgc variables
     deallocate (donor_pool       )
     deallocate (receiver_pool    )
     deallocate (floating_cn_ratio)
     deallocate (initial_cn_ratio )
     deallocate (rf_decomp      )
     deallocate (pathfrac_decomp)
     deallocate (is_cwd         )
     deallocate (is_litter      )
     deallocate (is_soil        )
     deallocate (gdp_lf         )
     deallocate (abm_lf         )
     deallocate (peatf_lf       )
#endif
        end if
     end if

#ifdef PFT_CLASSIFICATION
     CALL deallocate_PFTimeInvars
#endif

#ifdef PC_CLASSIFICATION
     CALL deallocate_PCTimeInvars
#endif

  END SUBROUTINE deallocate_TimeInvariants

#ifdef CLMDEBUG
   !---------------------------------------
   SUBROUTINE check_TimeInvariants ()

      use spmd_task
      use mod_colm_debug
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
#endif

      IMPLICIT NONE
      
      if (p_is_master) then
         write(*,'(/,A29)') 'Checking Time Invariantes ...'
      end if

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

      call check_vector_data ('lakedepth   ', lakedepth   ) !
      call check_vector_data ('dz_lake     ', dz_lake     ) ! new lake scheme
      
      call check_vector_data ('soil_s_v_alb', soil_s_v_alb) ! albedo of visible of the saturated soil
      call check_vector_data ('soil_d_v_alb', soil_d_v_alb) ! albedo of visible of the dry soil
      call check_vector_data ('soil_s_n_alb', soil_s_n_alb) ! albedo of near infrared of the saturated soil
      call check_vector_data ('soil_d_n_alb', soil_d_n_alb) ! albedo of near infrared of the dry soil
      call check_vector_data ('vf_quartz   ', vf_quartz   ) ! volumetric fraction of quartz within mineral soil
      call check_vector_data ('vf_gravels  ', vf_gravels  ) ! volumetric fraction of gravels
      call check_vector_data ('vf_om       ', vf_om       ) ! volumetric fraction of organic matter
      call check_vector_data ('vf_sand     ', vf_sand     ) ! volumetric fraction of sand
      call check_vector_data ('wf_gravels  ', wf_gravels  ) ! gravimetric fraction of gravels
      call check_vector_data ('wf_sand     ', wf_sand     ) ! gravimetric fraction of sand
      call check_vector_data ('porsl       ', porsl       ) ! fraction of soil that is voids [-]
      call check_vector_data ('psi0        ', psi0        ) ! minimum soil suction [mm] (NOTE: "-" valued)
#ifdef Campbell_SOIL_MODEL
      call check_vector_data ('bsw         ', bsw         ) ! clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      call check_vector_data ('theta_r     ', theta_r     ) 
      call check_vector_data ('alpha_vgm   ', alpha_vgm   ) 
      call check_vector_data ('L_vgm       ', L_vgm       ) 
      call check_vector_data ('n_vgm       ', n_vgm       ) 
      call check_vector_data ('sc_vgm      ', sc_vgm      ) 
      call check_vector_data ('fc_vgm      ', fc_vgm      ) 
#endif
      call check_vector_data ('hksati      ', hksati      ) ! hydraulic conductivity at saturation [mm h2o/s]
      call check_vector_data ('csol        ', csol        ) ! heat capacity of soil solids [J/(m3 K)]
      call check_vector_data ('k_solids    ', k_solids    ) ! thermal conductivity of soil solids [W/m-K]
      call check_vector_data ('dksatu      ', dksatu      ) ! thermal conductivity of unfrozen saturated soil [W/m-K]
      call check_vector_data ('dksatf      ', dksatf      ) ! thermal conductivity of frozen saturated soil [W/m-K]
      call check_vector_data ('dkdry       ', dkdry       ) ! thermal conductivity for dry soil  [W/(m-K)]
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
      call check_vector_data ('BA_alpha    ', BA_alpha    ) ! alpha in Balland and Arp(2005) thermal conductivity scheme
      call check_vector_data ('BA_beta     ', BA_beta     ) ! beta in Balland and Arp(2005) thermal conductivity scheme
#endif

      call check_vector_data ('htop        ', htop        ) 
      call check_vector_data ('hbot        ', hbot        ) 

#ifdef USE_DEPTH_TO_BEDROCK
      call check_vector_data ('dbedrock    ', dbedrock    ) !
#endif

#ifdef BGC
      call check_vector_data ('rf_decomp      ',  rf_decomp      )
      call check_vector_data ('pathfrac_decomp',  pathfrac_decomp)
      call check_vector_data ('gdp_lf         ',  gdp_lf         )
      call check_vector_data ('abm_lf         ',  abm_lf         )
      call check_vector_data ('peatf_lf       ',  peatf_lf       )
#endif

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif
     
      if (p_is_master) then
         write(*,'(A7,E20.10)') 'zlnd  ', zlnd   ! roughness length for soil [m]
         write(*,'(A7,E20.10)') 'zsno  ', zsno   ! roughness length for snow [m]
         write(*,'(A7,E20.10)') 'csoilc', csoilc ! drag coefficient for soil under canopy [-]
         write(*,'(A7,E20.10)') 'dewmx ', dewmx  ! maximum dew
         write(*,'(A7,E20.10)') 'wtfact', wtfact ! fraction of model area with high water table
         write(*,'(A7,E20.10)') 'capr  ', capr   ! tuning factor to turn first layer T into surface T
         write(*,'(A7,E20.10)') 'cnfac ', cnfac  ! Crank Nicholson factor between 0 and 1
         write(*,'(A7,E20.10)') 'ssi   ', ssi    ! irreducible water saturation of snow
         write(*,'(A7,E20.10)') 'wimp  ', wimp   ! water impremeable if porosity less than wimp
         write(*,'(A7,E20.10)') 'pondmx', pondmx ! ponding depth (mm)
         write(*,'(A7,E20.10)') 'smpmax', smpmax ! wilting point potential in mm
         write(*,'(A7,E20.10)') 'smpmin', smpmin ! restriction for min of soil poten. (mm)
         write(*,'(A7,E20.10)') 'trsmx0', trsmx0 ! max transpiration for moist soil+100% veg.  [mm/s]
         write(*,'(A7,E20.10)') 'tcrit ', tcrit  ! critical temp. to determine rain or snow
      end if

#ifdef PFT_CLASSIFICATION
     CALL check_PFTimeInvars
#endif

#ifdef PC_CLASSIFICATION
     CALL check_PCTimeInvars
#endif

   end subroutine check_TimeInvariants 
#endif

END MODULE MOD_TimeInvariants
