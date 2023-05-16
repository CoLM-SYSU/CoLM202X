#include <define.h>
#ifdef BGC
module bgc_veg_CNFireBaseMod

use precision
use PFT_Const, only: &
    cc_leaf   , cc_lstem   , cc_dstem , cc_other, fm_leaf, fm_lstem, fm_lroot, fm_root, fm_droot, fm_other, &
    fr_fcel   , fr_flig    , fr_flab  , lf_fcel , lf_flig, lf_flab
use MOD_TimeInvariants, only: &
    cmb_cmplt_fact, patchlatr, borealat, is_cwd, is_litter

use MOD_TimeVariables, only: &
    ! decomposition pools & fluxes variables (inout)
    decomp_cpools_vr, decomp_npools_vr,  cropf, farea_burned, baf_crop, baf_peatf, totsomc

use MOD_1D_Fluxes, only: &
    m_decomp_cpools_to_fire_vr, m_decomp_npools_to_fire_vr, &
    fire_mortality_to_met_c, fire_mortality_to_cel_c, fire_mortality_to_lig_c, fire_mortality_to_cwdc, &
    fire_mortality_to_met_n, fire_mortality_to_cel_n, fire_mortality_to_lig_n, fire_mortality_to_cwdn, &
    somc_fire

use MOD_PFTimeVars, only: &
    leafc_p     , leafc_storage_p     , leafc_xfer_p     , frootc_p    , frootc_storage_p    , frootc_xfer_p    , &
    livestemc_p , livestemc_storage_p , livestemc_xfer_p , deadstemc_p , deadstemc_storage_p , deadstemc_xfer_p , &
    livecrootc_p, livecrootc_storage_p, livecrootc_xfer_p, deadcrootc_p, deadcrootc_storage_p, deadcrootc_xfer_p, &
    leafn_p     , leafn_storage_p     , leafn_xfer_p     , frootn_p    , frootn_storage_p    , frootn_xfer_p    , &
    livestemn_p , livestemn_storage_p , livestemn_xfer_p , deadstemn_p , deadstemn_storage_p , deadstemn_xfer_p , &
    livecrootn_p, livecrootn_storage_p, livecrootn_xfer_p, deadcrootn_p, deadcrootn_storage_p, deadcrootn_xfer_p, &
    livecrootn_p, livecrootn_storage_p, livecrootn_xfer_p, deadcrootn_p, deadcrootn_storage_p, deadcrootn_xfer_p, &
    gresp_xfer_p, gresp_storage_p     , retransn_p       , &
    leaf_prof_p , froot_prof_p        , croot_prof_p     , stem_prof_p

use MOD_1D_PFTFluxes, only: &
    m_leafc_to_fire_p       , m_leafc_storage_to_fire_p       , m_leafc_xfer_to_fire_p     , &
    m_frootc_to_fire_p      , m_frootc_storage_to_fire_p      , m_frootc_xfer_to_fire_p    , &
    m_livestemc_to_fire_p   , m_livestemc_storage_to_fire_p   , m_livestemc_xfer_to_fire_p , &
    m_deadstemc_to_fire_p   , m_deadstemc_storage_to_fire_p   , m_deadstemc_xfer_to_fire_p , &
    m_livecrootc_to_fire_p  , m_livecrootc_storage_to_fire_p  , m_livecrootc_xfer_to_fire_p, &
    m_deadcrootc_to_fire_p  , m_deadcrootc_storage_to_fire_p  , m_deadcrootc_xfer_to_fire_p, &
    m_livestemc_to_deadstemc_fire_p, m_livecrootc_to_deadcrootc_fire_p, &
    m_gresp_xfer_to_fire_p  , m_gresp_storage_to_fire_p       , m_retransn_to_fire_p       , &
    m_leafn_to_fire_p       , m_leafn_storage_to_fire_p       , m_leafn_xfer_to_fire_p     , &
    m_frootn_to_fire_p      , m_frootn_storage_to_fire_p      , m_frootn_xfer_to_fire_p    , &
    m_livestemn_to_fire_p   , m_livestemn_storage_to_fire_p   , m_livestemn_xfer_to_fire_p , &
    m_deadstemn_to_fire_p   , m_deadstemn_storage_to_fire_p   , m_deadstemn_xfer_to_fire_p , &
    m_livecrootn_to_fire_p  , m_livecrootn_storage_to_fire_p  , m_livecrootn_xfer_to_fire_p, &
    m_deadcrootn_to_fire_p  , m_deadcrootn_storage_to_fire_p  , m_deadcrootn_xfer_to_fire_p, &
    m_livestemn_to_deadstemn_fire_p, m_livecrootn_to_deadcrootn_fire_p, &
      
    m_leafc_to_litter_fire_p     , m_leafc_storage_to_litter_fire_p     , m_leafc_xfer_to_litter_fire_p     , &
    m_frootc_to_litter_fire_p    , m_frootc_storage_to_litter_fire_p    , m_frootc_xfer_to_litter_fire_p    , &
    m_livestemc_to_litter_fire_p , m_livestemc_storage_to_litter_fire_p , m_livestemc_xfer_to_litter_fire_p , &
    m_deadstemc_to_litter_fire_p , m_deadstemc_storage_to_litter_fire_p , m_deadstemc_xfer_to_litter_fire_p , &
    m_livecrootc_to_litter_fire_p, m_livecrootc_storage_to_litter_fire_p, m_livecrootc_xfer_to_litter_fire_p, &
    m_deadcrootc_to_litter_fire_p, m_deadcrootc_storage_to_litter_fire_p, m_deadcrootc_xfer_to_litter_fire_p, &
    m_gresp_xfer_to_litter_fire_p, m_gresp_storage_to_litter_fire_p     , m_retransn_to_litter_fire_p       , &
    m_leafn_to_litter_fire_p     , m_leafn_storage_to_litter_fire_p     , m_leafn_xfer_to_litter_fire_p     , &
    m_frootn_to_litter_fire_p    , m_frootn_storage_to_litter_fire_p    , m_frootn_xfer_to_litter_fire_p    , &
    m_livestemn_to_litter_fire_p , m_livestemn_storage_to_litter_fire_p , m_livestemn_xfer_to_litter_fire_p , &
    m_deadstemn_to_litter_fire_p , m_deadstemn_storage_to_litter_fire_p , m_deadstemn_xfer_to_litter_fire_p , &
    m_livecrootn_to_litter_fire_p, m_livecrootn_storage_to_litter_fire_p, m_livecrootn_xfer_to_litter_fire_p, &
    m_deadcrootn_to_litter_fire_p, m_deadcrootn_storage_to_litter_fire_p, m_deadcrootn_xfer_to_litter_fire_p
   
use MOD_PFTimeInvars, only: pftfrac

implicit none

public CNFireFluxes

contains

subroutine CNFireFluxes(i,ps,pe,dlat,nl_soil,ndecomp_pools)

    integer ,intent(in) :: i
    integer ,intent(in) :: ps
    integer ,intent(in) :: pe
    real(r8),intent(in) :: dlat
    integer ,intent(in) :: nl_soil
    integer ,intent(in) :: ndecomp_pools

    ! !LOCAL VARIABLES:
    integer :: j,l   ! indices
    real(r8):: f
    real(r8):: mort
    INTEGER :: ivt, m

    integer, parameter :: lit_fp = 1   ! Pool for liter
    integer, parameter :: cwd_fp = 2   ! Pool for CWD Course woody debris


     !num_actfirep = 0
!     do fp = 1,num_soilp
!        p = filter_soilp(fp)
!        c = patch%column(p)

     do m = ps, pe
        if(cropf(i) < 1.0_r8)then
           ! For non-crop (bare-soil and natural vegetation)
!           if (transient_landcover) then
!              f = (fbac(c)-baf_crop(c))/(1.0_r8-cropf_col(c))
!           else
              f = (farea_burned(i)-baf_crop(i))/(1.0_r8-cropf(i))
!           end if
        else
           ! For crops
           if(cropf(i) > 0._r8)then
             f = baf_crop(i) /cropf(i)
           else
             f = 0._r8
           end if
        end if

        ! apply this rate to the patch state variables to get flux rates
        ! biomass burning
        ! carbon fluxes
        mort = 1._r8
!        if (spinup_state == 2) then
!           m = 10._r8
!        end if

!        if(f /= 0)then
!           num_actfirep = num_actfirep + 1
!           filter_actfirep(num_actfirep) = p
!        end if
        m_leafc_to_fire_p(m)               =  leafc_p(m)              * f * cc_leaf(ivt)
        m_leafc_storage_to_fire_p(m)       =  leafc_storage_p(m)      * f * cc_other(ivt)
        m_leafc_xfer_to_fire_p(m)          =  leafc_xfer_p(m)         * f * cc_other(ivt)
        m_livestemc_to_fire_p(m)           =  livestemc_p(m)          * f * cc_lstem(ivt)
        m_livestemc_storage_to_fire_p(m)   =  livestemc_storage_p(m)  * f * cc_other(ivt)
        m_livestemc_xfer_to_fire_p(m)      =  livestemc_xfer_p(m)     * f * cc_other(ivt)
        m_deadstemc_to_fire_p(m)           =  deadstemc_p(m)          * f * cc_dstem(ivt)
        m_deadstemc_storage_to_fire_p(m)   =  deadstemc_storage_p(m)  * f * cc_other(ivt)
        m_deadstemc_xfer_to_fire_p(m)      =  deadstemc_xfer_p(m)     * f * cc_other(ivt)
        m_frootc_to_fire_p(m)              =  frootc_p(m)             * f * 0._r8
        m_frootc_storage_to_fire_p(m)      =  frootc_storage_p(m)     * f * cc_other(ivt) 
        m_frootc_xfer_to_fire_p(m)         =  frootc_xfer_p(m)        * f * cc_other(ivt)
        m_livecrootc_to_fire_p(m)          =  livecrootc_p(m)         * f * 0._r8
        m_livecrootc_storage_to_fire_p(m)  =  livecrootc_storage_p(m) * f * cc_other(ivt) 
        m_livecrootc_xfer_to_fire_p(m)     =  livecrootc_xfer_p(m)    * f * cc_other(ivt) 
        m_deadcrootc_to_fire_p(m)          =  deadcrootc_p(m)         * f * 0._r8
        m_deadcrootc_storage_to_fire_p(m)  =  deadcrootc_storage_p(m) * f*  cc_other(ivt) 
        m_deadcrootc_xfer_to_fire_p(m)     =  deadcrootc_xfer_p(m)    * f * cc_other(ivt) 
        m_gresp_storage_to_fire_p(m)       =  gresp_storage_p(m)      * f * cc_other(ivt)
        m_gresp_xfer_to_fire_p(m)          =  gresp_xfer_p(m)         * f * cc_other(ivt)


        ! nitrogen fluxes
        m_leafn_to_fire_p(m)               =  leafn_p(m)              * f * cc_leaf(ivt)
        m_leafn_storage_to_fire_p(m)       =  leafn_storage_p(m)      * f * cc_other(ivt)
        m_leafn_xfer_to_fire_p(m)          =  leafn_xfer_p(m)         * f * cc_other(ivt)
        m_livestemn_to_fire_p(m)           =  livestemn_p(m)          * f * cc_lstem(ivt)
        m_livestemn_storage_to_fire_p(m)   =  livestemn_storage_p(m)  * f * cc_other(ivt)
        m_livestemn_xfer_to_fire_p(m)      =  livestemn_xfer_p(m)     * f * cc_other(ivt)
        m_deadstemn_to_fire_p(m)           =  deadstemn_p(m)          * f * cc_dstem(ivt) 
        m_deadstemn_storage_to_fire_p(m)   =  deadstemn_storage_p(m)  * f * cc_other(ivt)
        m_deadstemn_xfer_to_fire_p(m)      =  deadstemn_xfer_p(m)     * f * cc_other(ivt)
        m_frootn_to_fire_p(m)              =  frootn_p(m)             * f * 0._r8
        m_frootn_storage_to_fire_p(m)      =  frootn_storage_p(m)     * f * cc_other(ivt)
        m_frootn_xfer_to_fire_p(m)         =  frootn_xfer_p(m)        * f * cc_other(ivt)
        m_livecrootn_to_fire_p(m)          =  livecrootn_p(m)         * f * 0._r8 
        m_livecrootn_storage_to_fire_p(m)  =  livecrootn_storage_p(m) * f * cc_other(ivt) 
        m_livecrootn_xfer_to_fire_p(m)     =  livecrootn_xfer_p(m)    * f * cc_other(ivt)
        m_deadcrootn_to_fire_p(m)          =  deadcrootn_p(m)         * f * 0._r8
        m_deadcrootn_xfer_to_fire_p(m)     =  deadcrootn_xfer_p(m)    * f * cc_other(ivt) 
        m_deadcrootn_storage_to_fire_p(m)  =  deadcrootn_storage_p(m) * f * cc_other(ivt)
        m_retransn_to_fire_p(m)            =  retransn_p(m)           * f * cc_other(ivt)

!        if(use_matrixcn)then
!           matrix_fitransfer(p,ileaf_to_iout_fic)        = matrix_fitransfer(p,ileaf_to_iout_fic)          + f * cc_leaf(patch%itype(p))
!           matrix_fitransfer(p,ileafst_to_iout_fic)      = matrix_fitransfer(p,ileafst_to_iout_fic)        + f * cc_other(patch%itype(p))
!           matrix_fitransfer(p,ileafxf_to_iout_fic)      = matrix_fitransfer(p,ileafxf_to_iout_fic)        + f * cc_other(patch%itype(p))
!           matrix_fitransfer(p,ilivestem_to_iout_fic)    = matrix_fitransfer(p,ilivestem_to_iout_fic)      + f * cc_lstem(patch%itype(p))
!           matrix_fitransfer(p,ilivestemst_to_iout_fic)  = matrix_fitransfer(p,ilivestemst_to_iout_fic)    + f * cc_other(patch%itype(p))
!           matrix_fitransfer(p,ilivestemxf_to_iout_fic)  = matrix_fitransfer(p,ilivestemxf_to_iout_fic)    + f * cc_other(patch%itype(p))
!           matrix_fitransfer(p,ideadstem_to_iout_fic)    = matrix_fitransfer(p,ideadstem_to_iout_fic)      + f * cc_dstem(patch%itype(p))*m
!           matrix_fitransfer(p,ideadstemst_to_iout_fic)  = matrix_fitransfer(p,ideadstemst_to_iout_fic)    + f * cc_other(patch%itype(p))
!           matrix_fitransfer(p,ideadstemxf_to_iout_fic)  = matrix_fitransfer(p,ideadstemxf_to_iout_fic)    + f * cc_other(patch%itype(p))
!           matrix_fitransfer(p,ifroot_to_iout_fic)       = matrix_fitransfer(p,ifroot_to_iout_fic)         + f * 0._r8
!           matrix_fitransfer(p,ifrootst_to_iout_fic)     = matrix_fitransfer(p,ifrootst_to_iout_fic)       + f * cc_other(patch%itype(p))
!           matrix_fitransfer(p,ifrootxf_to_iout_fic)     = matrix_fitransfer(p,ifrootxf_to_iout_fic)       + f * cc_other(patch%itype(p))
!           matrix_fitransfer(p,ilivecroot_to_iout_fic)   = matrix_fitransfer(p,ilivecroot_to_iout_fic)     + f * 0._r8
!           matrix_fitransfer(p,ilivecrootst_to_iout_fic) = matrix_fitransfer(p,ilivecrootst_to_iout_fic)   + f * cc_other(patch%itype(p))
!           matrix_fitransfer(p,ilivecrootxf_to_iout_fic) = matrix_fitransfer(p,ilivecrootxf_to_iout_fic)   + f * cc_other(patch%itype(p))
!           matrix_fitransfer(p,ideadcroot_to_iout_fic)   = matrix_fitransfer(p,ideadcroot_to_iout_fic)     + f * 0._r8
!           matrix_fitransfer(p,ideadcrootst_to_iout_fic) = matrix_fitransfer(p,ideadcrootst_to_iout_fic)   + f * cc_other(patch%itype(p))
!           matrix_fitransfer(p,ideadcrootxf_to_iout_fic) = matrix_fitransfer(p,ideadcrootxf_to_iout_fic)   + f * cc_other(patch%itype(p))
!
!           matrix_nfitransfer(p,ileaf_to_iout_fin)        = matrix_nfitransfer(p,ileaf_to_iout_fin)        + f * cc_leaf(patch%itype(p))
!           matrix_nfitransfer(p,ileafst_to_iout_fin)      = matrix_nfitransfer(p,ileafst_to_iout_fin)      + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,ileafxf_to_iout_fin)      = matrix_nfitransfer(p,ileafxf_to_iout_fin)      + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,ilivestem_to_iout_fin)    = matrix_nfitransfer(p,ilivestem_to_iout_fin)    + f * cc_lstem(patch%itype(p))
!           matrix_nfitransfer(p,ilivestemst_to_iout_fin)  = matrix_nfitransfer(p,ilivestemst_to_iout_fin)  + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,ilivestemxf_to_iout_fin)  = matrix_nfitransfer(p,ilivestemxf_to_iout_fin)  + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,ideadstem_to_iout_fin)    = matrix_nfitransfer(p,ideadstem_to_iout_fin)    + f * cc_dstem(patch%itype(p))*m
!           matrix_nfitransfer(p,ideadstemst_to_iout_fin)  = matrix_nfitransfer(p,ideadstemst_to_iout_fin)  + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,ideadstemxf_to_iout_fin)  = matrix_nfitransfer(p,ideadstemxf_to_iout_fin)  + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,ifroot_to_iout_fin)       = matrix_nfitransfer(p,ifroot_to_iout_fin)       + f * 0._r8
!           matrix_nfitransfer(p,ifrootst_to_iout_fin)     = matrix_nfitransfer(p,ifrootst_to_iout_fin)     + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,ifrootxf_to_iout_fin)     = matrix_nfitransfer(p,ifrootxf_to_iout_fin)     + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,ilivecroot_to_iout_fin)   = matrix_nfitransfer(p,ilivecroot_to_iout_fin)   + f * 0._r8
!           matrix_nfitransfer(p,ilivecrootst_to_iout_fin) = matrix_nfitransfer(p,ilivecrootst_to_iout_fin) + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,ilivecrootxf_to_iout_fin) = matrix_nfitransfer(p,ilivecrootxf_to_iout_fin) + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,ideadcroot_to_iout_fin)   = matrix_nfitransfer(p,ideadcroot_to_iout_fin)   + f * 0._r8
!           matrix_nfitransfer(p,ideadcrootst_to_iout_fin) = matrix_nfitransfer(p,ideadcrootst_to_iout_fin) + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,ideadcrootxf_to_iout_fin) = matrix_nfitransfer(p,ideadcrootxf_to_iout_fin) + f * cc_other(patch%itype(p))
!           matrix_nfitransfer(p,iretransn_to_iout_fin)    = matrix_nfitransfer(p,iretransn_to_iout_fin)    + f * cc_other(patch%itype(p))
!        end if
        ! mortality due to fire
        ! carbon pools
        m_leafc_to_litter_fire_p(m)                   =  leafc_p(m) * f * &
             (1._r8 - cc_leaf(ivt)) * &
             fm_leaf(ivt)
        m_leafc_storage_to_litter_fire_p(m)           =  leafc_storage_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        m_leafc_xfer_to_litter_fire_p(m)              =  leafc_xfer_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from livestem/livecroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_livestemc_to_litter_fire_p(m)               =  livestemc_p(m) * f * &
             (1._r8 - cc_lstem(ivt)) * &
             fm_droot(ivt)    
        m_livestemc_storage_to_litter_fire_p(m)       =  livestemc_storage_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        m_livestemc_xfer_to_litter_fire_p(m)          =  livestemc_xfer_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt) 
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent the fraction of plant-tissue mortality for deadstem/deadcroot
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516
        m_livestemc_to_deadstemc_fire_p(m)            =  livestemc_p(m) * f * &
             (1._r8 - cc_lstem(ivt)) * &
             (fm_lstem(ivt)-fm_droot(ivt))
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from deadstem/deadcroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_deadstemc_to_litter_fire_p(m)               =  deadstemc_p(m) * f * m * &
             (1._r8 - cc_dstem(ivt)) * &
             fm_droot(ivt)    
        m_deadstemc_storage_to_litter_fire_p(m)       =  deadstemc_storage_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        m_deadstemc_xfer_to_litter_fire_p(m)          =  deadstemc_xfer_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        m_frootc_to_litter_fire_p(m)                  =  frootc_p(m)             * f * &
             fm_root(ivt)
        m_frootc_storage_to_litter_fire_p(m)          =  frootc_storage_p(m)     * f * &
             (1._r8- cc_other(ivt)) * &
             fm_other(ivt)
        m_frootc_xfer_to_litter_fire_p(m)             =  frootc_xfer_p(m)        * f * &
             (1._r8- cc_other(ivt)) * &
             fm_other(ivt)
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from livestem/livecroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_livecrootc_to_litter_fire_p(m)              =  livecrootc_p(m)         * f * &
             fm_droot(ivt)
        m_livecrootc_storage_to_litter_fire_p(m)      =  livecrootc_storage_p(m) * f * &
             (1._r8- cc_other(ivt)) * &
             fm_other(ivt) 
        m_livecrootc_xfer_to_litter_fire_p(m)         =  livecrootc_xfer_p(m)    * f * &
             (1._r8- cc_other(ivt)) * &
             fm_other(ivt) 
        m_livecrootc_to_deadcrootc_fire_p(m)          =  livecrootc_p(m)         * f * &
             (fm_lroot(ivt)-fm_droot(ivt))
        m_deadcrootc_to_litter_fire_p(m)              =  deadcrootc_p(m)         * f * m * &
             fm_droot(ivt)
        m_deadcrootc_storage_to_litter_fire_p(m)      =  deadcrootc_storage_p(m) * f * &
             (1._r8- cc_other(ivt)) * &
             fm_other(ivt)
        m_deadcrootc_xfer_to_litter_fire_p(m)         =  deadcrootc_xfer_p(m)    * f * &
             (1._r8- cc_other(ivt)) * &
             fm_other(ivt)      
        m_gresp_storage_to_litter_fire_p(m)           =  gresp_storage_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)  
        m_gresp_xfer_to_litter_fire_p(m)              =  gresp_xfer_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt) 


        ! nitrogen pools    
        m_leafn_to_litter_fire_p(m)                  =  leafn_p(m) * f * &
             (1._r8 - cc_leaf(ivt)) * &
             fm_leaf(ivt)
        m_leafn_storage_to_litter_fire_p(m)          =  leafn_storage_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)  
        m_leafn_xfer_to_litter_fire_p(m)             =  leafn_xfer_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from livestem/livecroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_livestemn_to_litter_fire_p(m)              =  livestemn_p(m) * f * &
             (1._r8 - cc_lstem(ivt)) * &
             fm_droot(ivt)
        m_livestemn_storage_to_litter_fire_p(m)      =  livestemn_storage_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)   
        m_livestemn_xfer_to_litter_fire_p(m)         =  livestemn_xfer_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent the fraction of plant-tissue mortality for deadstem/deadcroot
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516
        m_livestemn_to_deadstemn_fire_p(m)           =  livestemn_p(m) * f * &
             (1._r8 - cc_lstem(ivt)) * &
             (fm_lstem(ivt)-fm_droot(ivt))
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from deadstem/deadcroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_deadstemn_to_litter_fire_p(m)              =  deadstemn_p(m) * f * m * &
             (1._r8 - cc_dstem(ivt)) * &
             fm_droot(ivt)    
        m_deadstemn_storage_to_litter_fire_p(m)      =  deadstemn_storage_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        m_deadstemn_xfer_to_litter_fire_p(m)         =  deadstemn_xfer_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        m_frootn_to_litter_fire_p(m)                 =  frootn_p(m)             * f * &
             fm_root(ivt)
        m_frootn_storage_to_litter_fire_p(m)         =  frootn_storage_p(m)     * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        m_frootn_xfer_to_litter_fire_p(m)            =  frootn_xfer_p(m)        * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        ! NOTE: It looks incorrect to use fm_droot here, but it's used to represent fraction of transport from livestem/livecroot to litter
        ! EBK Oct/06/2017 see bug 2516 http://bugs.cgd.ucar.edu/show_bug.cgi?id=2516 (stem and root live or dead assumed to have the same transport)
        m_livecrootn_to_litter_fire_p(m)             =  livecrootn_p(m)         * f * &
             fm_droot(ivt)
        m_livecrootn_storage_to_litter_fire_p(m)     =  livecrootn_storage_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        m_livecrootn_xfer_to_litter_fire_p(m)        =  livecrootn_xfer_p(m)    * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt) 
        m_livecrootn_to_deadcrootn_fire_p(m)         =  livecrootn_p(m)         * f * &
             (fm_lroot(ivt)-fm_droot(ivt))
        m_deadcrootn_to_litter_fire_p(m)             =  deadcrootn_p(m)         * f * &
             fm_droot(ivt)
        m_deadcrootn_storage_to_litter_fire_p(m)     =  deadcrootn_storage_p(m) * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        m_deadcrootn_xfer_to_litter_fire_p(m)        =  deadcrootn_xfer_p(m)    * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt)
        m_retransn_to_litter_fire_p(m)               =  retransn_p(m)           * f * &
             (1._r8 - cc_other(ivt)) * &
             fm_other(ivt) 

!        if(use_matrixcn)then
!           matrix_fitransfer(p,ileaf_to_iout_fic)             = matrix_fitransfer(p,ileaf_to_iout_fic) &
!             + f * (1._r8 - cc_leaf(ivt))    * fm_leaf(ivt)
!           matrix_fitransfer(p,ileafst_to_iout_fic)           = matrix_fitransfer(p,ileafst_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt)
!           matrix_fitransfer(p,ileafxf_to_iout_fic)           = matrix_fitransfer(p,ileafxf_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt)
!           matrix_fitransfer(p,ilivestem_to_iout_fic)         = matrix_fitransfer(p,ilivestem_to_iout_fic) &
!             + f * (1._r8 - cc_lstem(ivt))   * fm_droot(ivt)
!           matrix_fitransfer(p,ilivestemst_to_iout_fic)       = matrix_fitransfer(p,ilivestemst_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt)
!           matrix_fitransfer(p,ilivestemxf_to_iout_fic)       = matrix_fitransfer(p,ilivestemxf_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt)
!           matrix_fitransfer(p,ilivestem_to_ideadstem_fic)    = matrix_fitransfer(p,ilivestem_to_ideadstem_fic) &
!             + f * (1._r8 - cc_lstem(ivt))   * (fm_lstem(ivt)-fm_droot(ivt))
!           matrix_fitransfer(p,ideadstem_to_iout_fic)         = matrix_fitransfer(p,ideadstem_to_iout_fic) &
!             + f * m*(1._r8 - cc_dstem(ivt)) * fm_droot(ivt)
!           matrix_fitransfer(p,ideadstemst_to_iout_fic)       = matrix_fitransfer(p,ideadstemst_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt)
!           matrix_fitransfer(p,ideadstemxf_to_iout_fic)       = matrix_fitransfer(p,ideadstemxf_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt)
!           matrix_fitransfer(p,ifroot_to_iout_fic)            = matrix_fitransfer(p,ifroot_to_iout_fic) &
!             + f * fm_root(ivt)
!           matrix_fitransfer(p,ifrootst_to_iout_fic)          = matrix_fitransfer(p,ifrootst_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt)
!           matrix_fitransfer(p,ifrootxf_to_iout_fic)          = matrix_fitransfer(p,ifrootxf_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt)
!           matrix_fitransfer(p,ilivecroot_to_iout_fic)        = matrix_fitransfer(p,ilivecroot_to_iout_fic) &
!             + f * fm_droot(ivt)
!           matrix_fitransfer(p,ilivecrootst_to_iout_fic)      = matrix_fitransfer(p,ilivecrootst_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt) 
!           matrix_fitransfer(p,ilivecrootxf_to_iout_fic)      = matrix_fitransfer(p,ilivecrootxf_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt) 
!           matrix_fitransfer(p,ilivecroot_to_ideadcroot_fic)  = matrix_fitransfer(p,ilivecroot_to_ideadcroot_fic) &
!             + f * (fm_lroot(ivt)-fm_droot(ivt))
!           matrix_fitransfer(p,ideadcroot_to_iout_fic)        = matrix_fitransfer(p,ideadcroot_to_iout_fic) &
!             + f * m * fm_droot(ivt)
!           matrix_fitransfer(p,ideadcrootst_to_iout_fic)      = matrix_fitransfer(p,ideadcrootst_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt) 
!           matrix_fitransfer(p,ideadcrootxf_to_iout_fic)      = matrix_fitransfer(p,ideadcrootxf_to_iout_fic) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt) 
!
!           matrix_nfitransfer(p,ileaf_to_iout_fin)            = matrix_nfitransfer(p,ileaf_to_iout_fin) &
!             + f * (1._r8 - cc_leaf(ivt))    * fm_leaf(ivt)
!           matrix_nfitransfer(p,ileafst_to_iout_fin)          = matrix_nfitransfer(p,ileafst_to_iout_fin) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt)
!           matrix_nfitransfer(p,ileafxf_to_iout_fin)          = matrix_nfitransfer(p,ileafxf_to_iout_fin) &
!             + f * (1._r8 - cc_other(ivt))   * fm_other(ivt)
!           matrix_nfitransfer(p,ilivestem_to_iout_fin)        = matrix_nfitransfer(p,ilivestem_to_iout_fin) &
!             + f * (1._r8 - cc_lstem(ivt))   * fm_droot(ivt)
!           matrix_nfitransfer(p,ilivestemst_to_iout_fin)      = matrix_nfitransfer(p,ilivestemst_to_iout_fin) &
!             + f * (1._r8 - cc_other(patch%itype(p)))   * fm_other(patch%itype(p))
!           matrix_nfitransfer(p,ilivestemxf_to_iout_fin)      = matrix_nfitransfer(p,ilivestemxf_to_iout_fin) &
!             + f * (1._r8 - cc_other(patch%itype(p)))   * fm_other(patch%itype(p))
!           matrix_nfitransfer(p,ilivestem_to_ideadstem_fin)   = matrix_nfitransfer(p,ilivestem_to_ideadstem_fin) &
!             + f * (1._r8 - cc_lstem(patch%itype(p)))   * (fm_lstem(patch%itype(p))-fm_droot(patch%itype(p)))
!           matrix_nfitransfer(p,ideadstem_to_iout_fin)        = matrix_nfitransfer(p,ideadstem_to_iout_fin) &
!             + f * m*(1._r8 - cc_dstem(patch%itype(p))) * fm_droot(patch%itype(p))
!           matrix_nfitransfer(p,ideadstemst_to_iout_fin)      = matrix_nfitransfer(p,ideadstemst_to_iout_fin) &
!             + f * (1._r8 - cc_other(patch%itype(p)))   * fm_other(patch%itype(p))
!           matrix_nfitransfer(p,ideadstemxf_to_iout_fin)      = matrix_nfitransfer(p,ideadstemxf_to_iout_fin) &
!             + f * (1._r8 - cc_other(patch%itype(p)))   * fm_other(patch%itype(p))
!           matrix_nfitransfer(p,ifroot_to_iout_fin)           = matrix_nfitransfer(p,ifroot_to_iout_fin) &
!             + f * fm_root(patch%itype(p))
!           matrix_nfitransfer(p,ifrootst_to_iout_fin)         = matrix_nfitransfer(p,ifrootst_to_iout_fin) &
!             + f * (1._r8 - cc_other(patch%itype(p)))   * fm_other(patch%itype(p))
!           matrix_nfitransfer(p,ifrootxf_to_iout_fin)         = matrix_nfitransfer(p,ifrootxf_to_iout_fin) &
!             + f * (1._r8 - cc_other(patch%itype(p)))   * fm_other(patch%itype(p))
!           matrix_nfitransfer(p,ilivecroot_to_iout_fin)       = matrix_nfitransfer(p,ilivecroot_to_iout_fin) &
!             + f * fm_droot(patch%itype(p))
!           matrix_nfitransfer(p,ilivecrootst_to_iout_fin)     = matrix_nfitransfer(p,ilivecrootst_to_iout_fin) &
!             + f * (1._r8 - cc_other(patch%itype(p)))   * fm_other(patch%itype(p)) 
!           matrix_nfitransfer(p,ilivecrootxf_to_iout_fin)     = matrix_nfitransfer(p,ilivecrootxf_to_iout_fin) &
!             + f * (1._r8 - cc_other(patch%itype(p)))   * fm_other(patch%itype(p)) 
!           matrix_nfitransfer(p,ilivecroot_to_ideadcroot_fin) = matrix_nfitransfer(p,ilivecroot_to_ideadcroot_fin) &
!             + f * (fm_lroot(patch%itype(p))-fm_droot(patch%itype(p)))
!           matrix_nfitransfer(p,ideadcroot_to_iout_fin)       = matrix_nfitransfer(p,ideadcroot_to_iout_fin) &
!             + f * m * fm_droot(patch%itype(p))
!           matrix_nfitransfer(p,ideadcrootst_to_iout_fin)     = matrix_nfitransfer(p,ideadcrootst_to_iout_fin) &
!             + f * (1._r8 - cc_other(patch%itype(p)))   * fm_other(patch%itype(p)) 
!           matrix_nfitransfer(p,ideadcrootxf_to_iout_fin)     = matrix_nfitransfer(p,ideadcrootxf_to_iout_fin) &
!             + f * (1._r8 - cc_other(patch%itype(p)))   * fm_other(patch%itype(p)) 
!           matrix_nfitransfer(p,iretransn_to_iout_fin)        = matrix_nfitransfer(p,iretransn_to_iout_fin) &
!             + f * (1._r8 - cc_other(patch%itype(p)))   * fm_other(patch%itype(p)) 
!        end if

!        if (use_cndv) then
!           if ( woody(patch%itype(p)) == 1._r8 )then
!              if ( livestemc(p)+deadstemc(p) > 0._r8 )then
!                 nind(p) = nind(p)*(1._r8-1._r8*fm_droot(patch%itype(p))*f) 
!              else
!                 nind(p) = 0._r8
!              end if
!           end if
!           leafcmax(p) = max(leafc(p)-m_leafc_to_fire(p)*dt, leafcmax(p))
!           if (patch%itype(p) == noveg) leafcmax(p) = 0._r8
!        end if

     end do  ! end of patches loop  

     ! fire-induced transfer of carbon and nitrogen pools to litter and cwd

     do j = 1,nl_soil
        fire_mortality_to_cwdc (j,i) = 0._r8
        fire_mortality_to_cwdn (j,i) = 0._r8
        fire_mortality_to_met_c(j,i) = 0._r8
!        do pi = 1,max_patch_per_col
!           do fc = 1,num_soilc
!              c = filter_soilc(fc)
!              if (pi <=  col%npatches(c)) then
!                 p = col%patchi(c) + pi - 1
!                 if ( patch%active(p) ) then
        do m = ps, pe
           fire_mortality_to_cwdc(j,i) = fire_mortality_to_cwdc(j,i) + &
                         m_deadstemc_to_litter_fire_p(m) * stem_prof_p(j,m) * pftfrac(m)
           fire_mortality_to_cwdc(j,i) = fire_mortality_to_cwdc(j,i) + &
                         m_deadcrootc_to_litter_fire_p(m) * croot_prof_p(j,m) * pftfrac(m)
           fire_mortality_to_cwdn(j,i) = fire_mortality_to_cwdn(j,i) + &
                         m_deadstemn_to_litter_fire_p(m) * stem_prof_p(j,m) * pftfrac(m)
           fire_mortality_to_cwdn(j,i) = fire_mortality_to_cwdn(j,i) + &
                         m_deadcrootn_to_litter_fire_p(m) * croot_prof_p(j,m) * pftfrac(m)


           fire_mortality_to_cwdc(j,i) = fire_mortality_to_cwdc(j,i) + &
                         m_livestemc_to_litter_fire_p(m) * stem_prof_p(j,m) * pftfrac(m)
           fire_mortality_to_cwdc(j,i) = fire_mortality_to_cwdc(j,i) + &
                         m_livecrootc_to_litter_fire_p(m) * croot_prof_p(j,m) * pftfrac(m)
           fire_mortality_to_cwdn(j,i) = fire_mortality_to_cwdn(j,i) + &
                         m_livestemn_to_litter_fire_p(m) * stem_prof_p(j,m) * pftfrac(m)
           fire_mortality_to_cwdn(j,i) = fire_mortality_to_cwdn(j,i) + &
                         m_livecrootn_to_litter_fire_p(m) * croot_prof_p(j,m) * pftfrac(m)


!        m_c_to_litr_met_fire(j,i)=m_c_to_litr_met_fire(j,i) + &
           fire_mortality_to_met_c(j,i)=fire_mortality_to_met_c(j,i) &
                         +((m_leafc_to_litter_fire_p(m)*lf_flab(ivt) &
                         +  m_leafc_storage_to_litter_fire_p(m) &
                         +  m_leafc_xfer_to_litter_fire_p(m) &
                         +  m_gresp_storage_to_litter_fire_p(m) &
                         +  m_gresp_xfer_to_litter_fire_p(m)) * leaf_prof_p(j,m) &
                         + (m_frootc_to_litter_fire_p(m)*fr_flab(ivt) &
                         +  m_frootc_storage_to_litter_fire_p(m) &
                         +  m_frootc_xfer_to_litter_fire_p(m)) * froot_prof_p(j,m) &
                         + (m_livestemc_storage_to_litter_fire_p(m) &
                         +  m_livestemc_xfer_to_litter_fire_p(m) &
                         +  m_deadstemc_storage_to_litter_fire_p(m) &
                         +  m_deadstemc_xfer_to_litter_fire_p(m)) * stem_prof_p(j,m) &
                         + (m_livecrootc_storage_to_litter_fire_p(m) &
                         +  m_livecrootc_xfer_to_litter_fire_p(m) &
                         +  m_deadcrootc_storage_to_litter_fire_p(m) &
                         +  m_deadcrootc_xfer_to_litter_fire_p(m)) * croot_prof_p(j,m)) * pftfrac(m)
!        m_c_to_litr_cel_fire(j,i)=m_c_to_litr_cel_fire(j,i) + &
           fire_mortality_to_cel_c(j,i)=fire_mortality_to_cel_c(j,i) &
                         + (m_leafc_to_litter_fire_p(m)*lf_fcel(ivt)*leaf_prof_p(j,m) &
                         +  m_frootc_to_litter_fire_p(m)*fr_fcel(ivt)*froot_prof_p(j,m)) * pftfrac(m)
!        m_c_to_litr_lig_fire(j,i)=m_c_to_litr_lig_fire(j,i) + &
           fire_mortality_to_lig_c(j,i)=fire_mortality_to_lig_c(j,i) &
                         + (m_leafc_to_litter_fire_p(m)*lf_flig(ivt)*leaf_prof_p(j,m) &
                         + m_frootc_to_litter_fire_p(m)*fr_flig(ivt)*froot_prof_p(j,m)) * pftfrac(m)

!        m_n_to_litr_met_fire(j,i)=m_n_to_litr_met_fire(j,i) + &
           fire_mortality_to_met_n(j,i)=fire_mortality_to_met_n(j,i) &
                         + ((m_leafn_to_litter_fire_p(m)*lf_flab(ivt) &
                         +   m_leafn_storage_to_litter_fire_p(m) &
                         +   m_leafn_xfer_to_litter_fire_p(m) &
                         +   m_retransn_to_litter_fire_p(m)) *leaf_prof_p(j,m) &
                         +  (m_frootn_to_litter_fire_p(m)*fr_flab(ivt) &
                         +   m_frootn_storage_to_litter_fire_p(m) &
                         +   m_frootn_xfer_to_litter_fire_p(m))*froot_prof_p(j,m) &
                         +  (m_livestemn_storage_to_litter_fire_p(m) &
                         +   m_livestemn_xfer_to_litter_fire_p(m) &
                         +   m_deadstemn_storage_to_litter_fire_p(m) &
                         +   m_deadstemn_xfer_to_litter_fire_p(m))* stem_prof_p(j,m)&
                         +  (m_livecrootn_storage_to_litter_fire_p(m) &
                         +   m_livecrootn_xfer_to_litter_fire_p(m) &
                         +   m_deadcrootn_storage_to_litter_fire_p(m) &
                         +   m_deadcrootn_xfer_to_litter_fire_p(m)) * croot_prof_p(j,m)) * pftfrac(m)
!        m_n_to_litr_cel_fire(j,i)=m_n_to_litr_cel_fire(j,i) + &
           fire_mortality_to_cel_n(j,i)=fire_mortality_to_cel_n(j,i) &
                         +  (m_leafn_to_litter_fire_p(m)*lf_fcel(i)*leaf_prof_p(j,m) &
                         +   m_frootn_to_litter_fire_p(m)*fr_fcel(i)*froot_prof_p(j,m)) * pftfrac(m)
!        m_n_to_litr_lig_fire(j,i)=m_n_to_litr_lig_fire(j,i) + &
           fire_mortality_to_lig_n(j,i)=fire_mortality_to_lig_n(j,i) &
                         +  (m_leafn_to_litter_fire_p(m)*lf_flig(i)*leaf_prof_p(j,m) &
                         +   m_frootn_to_litter_fire_p(m)*fr_flig(i)*froot_prof_p(j,m)) * pftfrac(m)
!                 end if
!              end if
!           end do
        end do
     end do
     !
     ! vertically-resolved decomposing C/N fire loss   
     ! column loop
     !
!     num_actfirec = 0
!     do fc = 1,num_soilc
!        c = filter_soilc(fc)
!
!        f = farea_burned(c) 

!        if(f /= 0 .or. f /= baf_crop(c))then
!           num_actfirec = num_actfirec + 1
!           filter_actfirec(num_actfirec) = c
!        end if
     do j = 1, nl_soil
           ! carbon fluxes
        do l = 1, ndecomp_pools
           if ( is_litter(l) ) then
              m_decomp_cpools_to_fire_vr(j,l,i) = decomp_cpools_vr(j,l,i) * f * &
                      cmb_cmplt_fact(lit_fp)
!                 if(use_soil_matrixcn)then! matrix is the same for C and N in the fire.
!                    matrix_decomp_fire_k(c,j+nlevdecomp*(l-1)) = matrix_decomp_fire_k(c,j+nlevdecomp*(l-1)) &
!                     - f * cmb_cmplt_fact(lit_fp) * dt
!                 end if
           end if
           if ( is_cwd(l) ) then
              m_decomp_cpools_to_fire_vr(j,l,i) = decomp_cpools_vr(j,l,i) * &
                      (f-baf_crop(i)) * cmb_cmplt_fact(cwd_fp)
!                 if(use_soil_matrixcn)then
!                    matrix_decomp_fire_k(c,j+nlevdecomp*(l-1)) = matrix_decomp_fire_k(c,j+nlevdecomp*(l-1)) &
!                     - (f-baf_crop(c)) * cmb_cmplt_fact(cwd_fp) * dt
!                 end if
           end if
        end do

           ! nitrogen fluxes
        do l = 1, ndecomp_pools
           if ( is_litter(l) ) then
              m_decomp_npools_to_fire_vr(j,l,i) = decomp_npools_vr(j,l,i) * f * &
                      cmb_cmplt_fact(lit_fp)
           end if
           if ( is_cwd(l) ) then
              m_decomp_npools_to_fire_vr(j,l,i) = decomp_npools_vr(j,l,i) * &
                     (f-baf_crop(i)) * cmb_cmplt_fact(cwd_fp)
           end if
        end do

     end do
!     end do  ! end of column loop

     ! carbon loss due to deforestation fires

!     if (transient_landcover) then
!        call get_curr_date (kyr, kmo, kda, mcsec)
!        do fc = 1,num_soilc
!           c = filter_soilc(fc)
!           lfc2(c)=0._r8
!           if( .not. (kmo == 1 .and. kda == 1 .and. mcsec == 0) )then
!              if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 .and. dtrotr_col(c) > 0._r8 .and. &
!                   lfc(c) > 0._r8 .and. fbac1(c) == 0._r8) then
!                 lfc2(c) = max(0._r8, min(lfc(c), (farea_burned(c)-baf_crop(c) - &
!                      baf_peatf(c))/2.0*dt))/(dtrotr_col(c)*dayspyr*secspday/dt)/dt
!                 lfc(c)  = lfc(c) - max(0._r8, min(lfc(c), (farea_burned(c)-baf_crop(c) - &
!                      baf_peatf(c))*dt/2.0_r8))
!              end if
!           end if
!        end do
!     end if
     !
     ! Carbon loss due to peat fires
     !
     ! somc_fire is not connected to clm45 soil carbon pool, ie does not decrease
     ! soil carbon b/c clm45 soil carbon was very low in several peatland grids
     !
!     do fc = 1,num_soilc
!        c = filter_soilc(fc)
!        g = col%gridcell(c)
        if( patchlatr(i)  <  borealat)then
           somc_fire(i)= totsomc(i)*baf_peatf(i)*6.0_r8/33.9_r8
        else
           somc_fire(i)= baf_peatf(i)*2.2e3_r8
        end if
!     end do

     ! Fang Li has not added aerosol and trace gas emissions due to fire, yet
     ! They will be added here in proportion to the carbon emission
     ! Emission factors differ for various fire types


  end subroutine CNFireFluxes

end module bgc_veg_CNFireBaseMod
#endif
