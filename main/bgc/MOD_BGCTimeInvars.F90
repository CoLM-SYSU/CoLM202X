#include <define.h>

MODULE MOD_BGCTimeInvars 
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------
#ifdef BGC

use precision
IMPLICIT NONE
SAVE
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
  INTEGER , allocatable :: rice2pdt(:)

  REAL(r8) :: nitrif_n2o_loss_frac ! fraction of N lost as N2O in nitrification (Li et al., 2000)
  REAL(r8) :: dnp     ! denitrification proportion
  REAL(r8) :: bdnr    ! bulk denitrification rate (1/day)
  REAL(r8) :: compet_plant_no3
  REAL(r8) :: compet_plant_nh4
  REAL(r8) :: compet_decomp_no3
  REAL(r8) :: compet_decomp_nh4
  REAL(r8) :: compet_denit
  REAL(r8) :: compet_nit
  REAL(r8) :: surface_tension_water
  REAL(r8) :: rij_kro_a
  REAL(r8) :: rij_kro_alpha
  REAL(r8) :: rij_kro_beta
  REAL(r8) :: rij_kro_gamma
  REAL(r8) :: rij_kro_delta
  REAL(r8) :: nfix_timeconst
  REAL(r8) :: organic_max
  REAL(r8) :: d_con_g21
  REAL(r8) :: d_con_g22
  REAL(r8) :: d_con_w21
  REAL(r8) :: d_con_w22
  REAL(r8) :: d_con_w23
  REAL(r8) :: denit_resp_coef
  REAL(r8) :: denit_resp_exp
  REAL(r8) :: denit_nitrate_coef
  REAL(r8) :: denit_nitrate_exp
  REAL(r8) :: k_nitr_max
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
  public :: allocate_BGCTimeInvars
  public :: deallocate_BGCTimeInvars
  public :: READ_BGCTimeInvars
  public :: WRITE_BGCTimeInvars

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_BGCTimeInvars ()
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM 1d [numpatch] variables
  ! --------------------------------------------------------------------

     use precision
     use GlobalVars, only: nl_soil, ndecomp_transitions, ndecomp_pools
     use spmd_task
     use mod_landpatch, only : numpatch
     IMPLICIT NONE

  if (p_is_worker) then

     if (numpatch > 0) then
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
     allocate (rice2pdt          (numpatch))

! end bgc variables
  end if
  ENDIF

  END SUBROUTINE allocate_BGCTimeInvars

  !---------------------------------------
  SUBROUTINE READ_BGCTimeInvars (file_restart)

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

     IMPLICIT NONE

     character(LEN=*), intent(in) :: file_restart
  
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
     call ncio_read_vector       (file_restart, 'rice2pdt       ', landpatch, rice2pdt      )

     call ncio_read_bcast_serial (file_restart, 'nitrif_n2o_loss_frac', nitrif_n2o_loss_frac)
     call ncio_read_bcast_serial (file_restart, 'dnp                 ', dnp                 )!
     call ncio_read_bcast_serial (file_restart, 'bdnr                ', bdnr                )!
     call ncio_read_bcast_serial (file_restart, 'compet_plant_no3    ', compet_plant_no3    )!
     call ncio_read_bcast_serial (file_restart, 'compet_plant_nh4    ', compet_plant_nh4    )!
     call ncio_read_bcast_serial (file_restart, 'compet_decomp_no3   ', compet_decomp_no3   )!
     call ncio_read_bcast_serial (file_restart, 'compet_decomp_nh4   ', compet_decomp_nh4   )!
     call ncio_read_bcast_serial (file_restart, 'compet_denit        ', compet_denit        )!
     call ncio_read_bcast_serial (file_restart, 'compet_nit          ', compet_nit          )!
     call ncio_read_bcast_serial (file_restart, 'surface_tension_water',surface_tension_water)
     call ncio_read_bcast_serial (file_restart, 'rij_kro_a           ', rij_kro_a           )
     call ncio_read_bcast_serial (file_restart, 'rij_kro_alpha       ', rij_kro_alpha       )
     call ncio_read_bcast_serial (file_restart, 'rij_kro_beta        ', rij_kro_beta        )
     call ncio_read_bcast_serial (file_restart, 'rij_kro_gamma       ', rij_kro_gamma       )
     call ncio_read_bcast_serial (file_restart, 'rij_kro_delta       ', rij_kro_delta       )
     call ncio_read_bcast_serial (file_restart, 'nfix_timeconst      ', nfix_timeconst      )
     call ncio_read_bcast_serial (file_restart, 'organic_max         ', organic_max         )
     call ncio_read_bcast_serial (file_restart, 'd_con_g21           ', d_con_g21           )
     call ncio_read_bcast_serial (file_restart, 'd_con_g22           ', d_con_g22           )
     call ncio_read_bcast_serial (file_restart, 'd_con_w21           ', d_con_w21           )
     call ncio_read_bcast_serial (file_restart, 'd_con_w22           ', d_con_w22           )
     call ncio_read_bcast_serial (file_restart, 'd_con_w23           ', d_con_w23           )
     call ncio_read_bcast_serial (file_restart, 'denit_resp_coef     ', denit_resp_coef     )
     call ncio_read_bcast_serial (file_restart, 'denit_resp_exp      ', denit_resp_exp      )
     call ncio_read_bcast_serial (file_restart, 'denit_nitrate_coef  ', denit_nitrate_coef  )
     call ncio_read_bcast_serial (file_restart, 'denit_nitrate_exp   ', denit_nitrate_exp   )!
     call ncio_read_bcast_serial (file_restart, 'k_nitr_max          ', k_nitr_max          )!
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

#ifdef CLMDEBUG 
     call check_BGCTimeInvars ()
#endif

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

  end subroutine READ_BGCTimeInvars

  !---------------------------------------
  SUBROUTINE WRITE_BGCTimeInvars (file_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use mod_namelist, only : DEF_REST_COMPRESS_LEVEL 
     use spmd_task
     use ncio_serial
     use ncio_vector
     use mod_landpatch
     USE GlobalVars

     IMPLICIT NONE

     character(len=*), intent(in) :: file_restart
     
     ! Local Variables
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL 

     call ncio_create_file_vector (file_restart, landpatch)

     CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil', nl_soil)
     call ncio_define_dimension_vector (file_restart, landpatch, 'ndecomp_transitions',ndecomp_transitions)
     
     call ncio_write_vector       (file_restart, 'rf_decomp      ', 'soil'   , nl_soil  , &
                          'ndecomp_transitions', ndecomp_transitions,'patch', landpatch, rf_decomp      , compress)
     call ncio_write_vector       (file_restart, 'pathfrac_decomp', 'soil'   , nl_soil  , &
                          'ndecomp_transitions', ndecomp_transitions,'patch', landpatch, pathfrac_decomp, compress)
     call ncio_write_vector       (file_restart, 'gdp_lf         ',  'patch', landpatch, gdp_lf         , compress)
     call ncio_write_vector       (file_restart, 'abm_lf         ',  'patch', landpatch, abm_lf         , compress)
     call ncio_write_vector       (file_restart, 'peatf_lf       ',  'patch', landpatch, peatf_lf       , compress)
     call ncio_write_vector       (file_restart, 'rice2pdt       ',  'patch', landpatch, rice2pdt       , compress)

     if (p_is_master) then

        call ncio_create_file (file_restart)
        call ncio_define_dimension(file_restart, 'ndecomp_transitions',ndecomp_transitions)
        call ncio_define_dimension(file_restart, 'ndecomp_pools'      ,ndecomp_pools)
        call ncio_define_dimension(file_restart, 'nlitter_fire'       ,2            )

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
        call ncio_write_serial (file_restart, 'compet_plant_no3    ', compet_plant_no3    )!
        call ncio_write_serial (file_restart, 'compet_plant_nh4    ', compet_plant_nh4    )!
        call ncio_write_serial (file_restart, 'compet_decomp_no3   ', compet_decomp_no3   )!
        call ncio_write_serial (file_restart, 'compet_decomp_nh4   ', compet_decomp_nh4   )!
        call ncio_write_serial (file_restart, 'compet_denit        ', compet_denit        )!
        call ncio_write_serial (file_restart, 'compet_nit          ', compet_nit          )!
        call ncio_write_serial (file_restart, 'surface_tension_water',surface_tension_water)
        call ncio_write_serial (file_restart, 'rij_kro_a           ', rij_kro_a           )
        call ncio_write_serial (file_restart, 'rij_kro_alpha       ', rij_kro_alpha       )
        call ncio_write_serial (file_restart, 'rij_kro_beta        ', rij_kro_beta        )
        call ncio_write_serial (file_restart, 'rij_kro_gamma       ', rij_kro_gamma       )
        call ncio_write_serial (file_restart, 'rij_kro_delta       ', rij_kro_delta       )
        call ncio_write_serial (file_restart, 'nfix_timeconst      ', nfix_timeconst      )
        call ncio_write_serial (file_restart, 'organic_max         ', organic_max         )
        call ncio_write_serial (file_restart, 'd_con_g21           ', d_con_g21           )
        call ncio_write_serial (file_restart, 'd_con_g22           ', d_con_g22           )
        call ncio_write_serial (file_restart, 'd_con_w21           ', d_con_w21           )
        call ncio_write_serial (file_restart, 'd_con_w22           ', d_con_w22           )
        call ncio_write_serial (file_restart, 'd_con_w23           ', d_con_w23           )
        call ncio_write_serial (file_restart, 'denit_resp_coef     ', denit_resp_coef     )
        call ncio_write_serial (file_restart, 'denit_resp_exp      ', denit_resp_exp      )
        call ncio_write_serial (file_restart, 'denit_nitrate_coef  ', denit_nitrate_coef  )
        call ncio_write_serial (file_restart, 'denit_nitrate_exp   ', denit_nitrate_exp   )!
        call ncio_write_serial (file_restart, 'k_nitr_max          ', k_nitr_max          )!
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

     end if

   end subroutine WRITE_BGCTimeInvars 

  SUBROUTINE deallocate_BGCTimeInvars ()

     use spmd_task
     use mod_landpatch, only : numpatch
     implicit none

     ! --------------------------------------------------
     ! Deallocates memory for CLM 1d [numpatch] variables
     ! --------------------------------------------------

     if (p_is_worker) then

        if (numpatch > 0) then

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
     deallocate (rice2pdt       )

        end if
     end if

  END SUBROUTINE deallocate_BGCTimeInvars

#ifdef CLMDEBUG
   !---------------------------------------
   SUBROUTINE check_BGCTimeInvars ()

      use spmd_task
      use mod_colm_debug

      IMPLICIT NONE
      
      call check_vector_data ('rf_decomp      ',  rf_decomp      )
      call check_vector_data ('pathfrac_decomp',  pathfrac_decomp)
      call check_vector_data ('gdp_lf         ',  gdp_lf         )
      call check_vector_data ('abm_lf         ',  abm_lf         )
      call check_vector_data ('peatf_lf       ',  peatf_lf       )
      call check_vector_data ('rice2pdt       ',  float(rice2pdt))
     
   end subroutine check_BGCTimeInvars 
#endif

#endif
END MODULE MOD_BGCTimeInvars
