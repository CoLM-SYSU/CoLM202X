#include <define.h>
#ifdef BGC
module bgc_CNSASUMod

!#include "shr_assert.h"
  !-----------------------------------------------------------------------
  ! The matrix model of CLM5.0 was developed by Yiqi Luo EcoLab members,
  ! Drs. Xingjie Lu, Yuanyuan Huang and Zhengguang Du, at Northern Arizona University
  !----------------------------------------------------------------------------------
  ! 
  ! DESCRIPTION:
  ! Module for CLM5.0BGC matrices
  ! The matrix equation 
  ! Xn+1 = Xn + I*dt + (A*K(ksi) - Kfire - tri/dz)*Xn*dt
  ! Or
  ! Xn+1 = Xn + I*dt + (A*K(ksi) - Kfire - V)*Xn*dt
  
  ! !USES:
  use precision
  use MOD_TimeInvariants, only: &
      i_met_lit, i_cel_lit, i_lig_lit, i_cwd, i_soil1, i_soil2, i_soil3, floating_cn_ratio

  use MOD_TimeVariables, only: &
      decomp_cpools_vr           , decomp_npools_vr           , decomp0_cpools_vr          , decomp0_npools_vr          , &
      I_met_c_vr_acc             , I_cel_c_vr_acc             , I_lig_c_vr_acc             , I_cwd_c_vr_acc             , &
      AKX_met_to_soil1_c_vr_acc  , AKX_cel_to_soil1_c_vr_acc  , AKX_lig_to_soil2_c_vr_acc  , AKX_soil1_to_soil2_c_vr_acc, &
      AKX_cwd_to_cel_c_vr_acc    , AKX_cwd_to_lig_c_vr_acc    , AKX_soil1_to_soil3_c_vr_acc, AKX_soil2_to_soil1_c_vr_acc, &
      AKX_soil2_to_soil3_c_vr_acc, AKX_soil3_to_soil1_c_vr_acc, &
      AKX_met_exit_c_vr_acc      , AKX_cel_exit_c_vr_acc      , AKX_lig_exit_c_vr_acc      , AKX_cwd_exit_c_vr_acc      , &
      AKX_soil1_exit_c_vr_acc    , AKX_soil2_exit_c_vr_acc    , AKX_soil3_exit_c_vr_acc    , &
      diagVX_c_vr_acc            , upperVX_c_vr_acc           , lowerVX_c_vr_acc           , &
      I_met_n_vr_acc             , I_cel_n_vr_acc             , I_lig_n_vr_acc             , I_cwd_n_vr_acc             , &
      AKX_met_to_soil1_n_vr_acc  , AKX_cel_to_soil1_n_vr_acc  , AKX_lig_to_soil2_n_vr_acc  , AKX_soil1_to_soil2_n_vr_acc, &
      AKX_cwd_to_cel_n_vr_acc    , AKX_cwd_to_lig_n_vr_acc    , AKX_soil1_to_soil3_n_vr_acc, AKX_soil2_to_soil1_n_vr_acc, &
      AKX_soil2_to_soil3_n_vr_acc, AKX_soil3_to_soil1_n_vr_acc, &
      AKX_met_exit_n_vr_acc      , AKX_cel_exit_n_vr_acc      , AKX_lig_exit_n_vr_acc      , AKX_cwd_exit_n_vr_acc      , &
      AKX_soil1_exit_n_vr_acc    , AKX_soil2_exit_n_vr_acc    , AKX_soil3_exit_n_vr_acc    , &
      diagVX_n_vr_acc            , upperVX_n_vr_acc           , lowerVX_n_vr_acc           , skip_balance_check         , &
      cn_decomp_pools 

use MOD_PFTimeInvars, only: pftclass
use MOD_PFTimeVars, only: &
      leafc_p           , leafc_storage_p      , leafc_xfer_p      , leafc0_p             , leafc0_storage_p     , leafc0_xfer_p     , &
      frootc_p          , frootc_storage_p     , frootc_xfer_p     , frootc0_p            , frootc0_storage_p    , frootc0_xfer_p    , &
      livestemc_p       , livestemc_storage_p  , livestemc_xfer_p  , livestemc0_p         , livestemc0_storage_p , livestemc0_xfer_p , &
      deadstemc_p       , deadstemc_storage_p  , deadstemc_xfer_p  , deadstemc0_p         , deadstemc0_storage_p , deadstemc0_xfer_p , &
      livecrootc_p      , livecrootc_storage_p , livecrootc_xfer_p , livecrootc0_p        , livecrootc0_storage_p, livecrootc0_xfer_p, &
      deadcrootc_p      , deadcrootc_storage_p , deadcrootc_xfer_p , deadcrootc0_p        , deadcrootc0_storage_p, deadcrootc0_xfer_p, &
      grainc_p          , grainc_storage_p     , grainc_xfer_p     , grainc0_p            , grainc0_storage_p    , grainc0_xfer_p    , &

      leafn_p           , leafn_storage_p      , leafn_xfer_p      , leafn0_p             , leafn0_storage_p     , leafn0_xfer_p     , &
      frootn_p          , frootn_storage_p     , frootn_xfer_p     , frootn0_p            , frootn0_storage_p    , frootn0_xfer_p    , &
      livestemn_p       , livestemn_storage_p  , livestemn_xfer_p  , livestemn0_p         , livestemn0_storage_p , livestemn0_xfer_p , &
      deadstemn_p       , deadstemn_storage_p  , deadstemn_xfer_p  , deadstemn0_p         , deadstemn0_storage_p , deadstemn0_xfer_p , &
      livecrootn_p      , livecrootn_storage_p , livecrootn_xfer_p , livecrootn0_p        , livecrootn0_storage_p, livecrootn0_xfer_p, &
      deadcrootn_p      , deadcrootn_storage_p , deadcrootn_xfer_p , deadcrootn0_p        , deadcrootn0_storage_p, deadcrootn0_xfer_p, &
      grainn_p          , grainn_storage_p     , grainn_xfer_p     , grainn0_p            , grainn0_storage_p    , grainn0_xfer_p    , &
      retransn_p        , retransn0_p          , &

      I_leafc_p_acc     , I_leafc_st_p_acc     , I_frootc_p_acc    , I_frootc_st_p_acc    , &
      I_livestemc_p_acc , I_livestemc_st_p_acc , I_deadstemc_p_acc , I_deadstemc_st_p_acc , &
      I_livecrootc_p_acc, I_livecrootc_st_p_acc, I_deadcrootc_p_acc, I_deadcrootc_st_p_acc, &
      I_grainc_p_acc    , I_grainc_st_p_acc    , &

      I_leafn_p_acc     , I_leafn_st_p_acc     , I_frootn_p_acc    , I_frootn_st_p_acc    , &
      I_livestemn_p_acc , I_livestemn_st_p_acc , I_deadstemn_p_acc , I_deadstemn_st_p_acc , &
      I_livecrootn_p_acc, I_livecrootn_st_p_acc, I_deadcrootn_p_acc, I_deadcrootn_st_p_acc, &
      I_grainn_p_acc    , I_grainn_st_p_acc    , &

      AKX_leafc_xf_to_leafc_p_acc           , AKX_frootc_xf_to_frootc_p_acc           , AKX_livestemc_xf_to_livestemc_p_acc     , &
      AKX_deadstemc_xf_to_deadstemc_p_acc   , AKX_livecrootc_xf_to_livecrootc_p_acc   , AKX_deadcrootc_xf_to_deadcrootc_p_acc   , &
      AKX_grainc_xf_to_grainc_p_acc         , AKX_livestemc_to_deadstemc_p_acc        , AKX_livecrootc_to_deadcrootc_p_acc      , &
           
      AKX_leafc_st_to_leafc_xf_p_acc        , AKX_frootc_st_to_frootc_xf_p_acc        , AKX_livestemc_st_to_livestemc_xf_p_acc  , &
      AKX_deadstemc_st_to_deadstemc_xf_p_acc, AKX_livecrootc_st_to_livecrootc_xf_p_acc, AKX_deadcrootc_st_to_deadcrootc_xf_p_acc, &
      AKX_grainc_st_to_grainc_xf_p_acc      , &

      AKX_leafc_exit_p_acc                  , AKX_frootc_exit_p_acc                   , AKX_livestemc_exit_p_acc                , &
      AKX_deadstemc_exit_p_acc              , AKX_livecrootc_exit_p_acc               , AKX_deadcrootc_exit_p_acc               , &
      AKX_grainc_exit_p_acc                 , &

      AKX_leafc_st_exit_p_acc               , AKX_frootc_st_exit_p_acc                , AKX_livestemc_st_exit_p_acc             , &
      AKX_deadstemc_st_exit_p_acc           , AKX_livecrootc_st_exit_p_acc            , AKX_deadcrootc_st_exit_p_acc            , &
      AKX_grainc_st_exit_p_acc              , &

      AKX_leafc_xf_exit_p_acc               , AKX_frootc_xf_exit_p_acc                , AKX_livestemc_xf_exit_p_acc             , &
      AKX_deadstemc_xf_exit_p_acc           , AKX_livecrootc_xf_exit_p_acc            , AKX_deadcrootc_xf_exit_p_acc            , &
      AKX_grainc_xf_exit_p_acc              , &
           
      AKX_leafn_xf_to_leafn_p_acc           , AKX_frootn_xf_to_frootn_p_acc           , AKX_livestemn_xf_to_livestemn_p_acc     , &
      AKX_deadstemn_xf_to_deadstemn_p_acc   , AKX_livecrootn_xf_to_livecrootn_p_acc   , AKX_deadcrootn_xf_to_deadcrootn_p_acc   , &
      AKX_grainn_xf_to_grainn_p_acc         , AKX_livestemn_to_deadstemn_p_acc        , AKX_livecrootn_to_deadcrootn_p_acc      , &

      AKX_leafn_st_to_leafn_xf_p_acc        , AKX_frootn_st_to_frootn_xf_p_acc        , AKX_livestemn_st_to_livestemn_xf_p_acc  , &
      AKX_deadstemn_st_to_deadstemn_xf_p_acc, AKX_livecrootn_st_to_livecrootn_xf_p_acc, AKX_deadcrootn_st_to_deadcrootn_xf_p_acc, &
      AKX_grainn_st_to_grainn_xf_p_acc      , &

      AKX_leafn_to_retransn_p_acc           , AKX_frootn_to_retransn_p_acc            , AKX_livestemn_to_retransn_p_acc         , &
      AKX_livecrootn_to_retransn_p_acc      , &

      AKX_retransn_to_leafn_p_acc           , AKX_retransn_to_frootn_p_acc            , AKX_retransn_to_livestemn_p_acc         , &
      AKX_retransn_to_deadstemn_p_acc       , AKX_retransn_to_livecrootn_p_acc        , AKX_retransn_to_deadcrootn_p_acc        , &
      AKX_retransn_to_grainn_p_acc          , &

      AKX_retransn_to_leafn_st_p_acc        , AKX_retransn_to_frootn_st_p_acc         , AKX_retransn_to_livestemn_st_p_acc      , &
      AKX_retransn_to_deadstemn_st_p_acc    , AKX_retransn_to_livecrootn_st_p_acc     , AKX_retransn_to_deadcrootn_st_p_acc     , &
      AKX_retransn_to_grainn_st_p_acc       , &

      AKX_leafn_exit_p_acc                  , AKX_frootn_exit_p_acc                   , AKX_livestemn_exit_p_acc                , &
      AKX_deadstemn_exit_p_acc              , AKX_livecrootn_exit_p_acc               , AKX_deadcrootn_exit_p_acc               , &
      AKX_grainn_exit_p_acc                 , AKX_retransn_exit_p_acc                 , &

      AKX_leafn_st_exit_p_acc               , AKX_frootn_st_exit_p_acc                , AKX_livestemn_st_exit_p_acc             , &
      AKX_deadstemn_st_exit_p_acc           , AKX_livecrootn_st_exit_p_acc            , AKX_deadcrootn_st_exit_p_acc            , &
      AKX_grainn_st_exit_p_acc              , &

      AKX_leafn_xf_exit_p_acc               , AKX_frootn_xf_exit_p_acc                , AKX_livestemn_xf_exit_p_acc             , &
      AKX_deadstemn_xf_exit_p_acc           , AKX_livecrootn_xf_exit_p_acc            , AKX_deadcrootn_xf_exit_p_acc            , &
      AKX_grainn_xf_exit_p_acc
!
  implicit none

  public :: CNSASU
  public :: inverse
  !-----------------------------------------------------------------------

contains


  !-----------------------------------------------------------------------
  subroutine CNSASU(i,ps,pe,deltim,idate,nl_soil,ndecomp_transitions, ndecomp_pools, ndecomp_pools_vr)
    ! !DESCRIPTION:
    ! !ARGUMENTS:

    integer, intent(in) :: i
    integer, intent(in) :: ps
    integer, intent(in) :: pe
    real(r8),intent(in) :: deltim
    integer ,intent(in) :: idate(3)
    integer, intent(in) :: nl_soil
    integer, intent(in) :: ndecomp_transitions
    integer, intent(in) :: ndecomp_pools
    integer, intent(in) :: ndecomp_pools_vr

    !-----------------------------------------------------------------------

    integer :: k, m, j
     ! set time steps
    real(r8),parameter :: epsi          = 1.e-8_r8 
    integer ,parameter :: nvegc         = 21
    integer ,parameter :: nvegn         = 22
    integer ,parameter :: ileaf         = 1
    integer ,parameter :: ileaf_st      = 2
    integer ,parameter :: ileaf_xf      = 3
    integer ,parameter :: ifroot        = 4
    integer ,parameter :: ifroot_st     = 5
    integer ,parameter :: ifroot_xf     = 6
    integer ,parameter :: ilivestem     = 7
    integer ,parameter :: ilivestem_st  = 8
    integer ,parameter :: ilivestem_xf  = 9
    integer ,parameter :: ideadstem     = 10
    integer ,parameter :: ideadstem_st  = 11
    integer ,parameter :: ideadstem_xf  = 12
    integer ,parameter :: ilivecroot    = 13
    integer ,parameter :: ilivecroot_st = 14
    integer ,parameter :: ilivecroot_xf = 15
    integer ,parameter :: ideadcroot    = 16
    integer ,parameter :: ideadcroot_st = 17
    integer ,parameter :: ideadcroot_xf = 18
    integer ,parameter :: igrain        = 19
    integer ,parameter :: igrain_st     = 20
    integer ,parameter :: igrain_xf     = 21
    integer ,parameter :: iretrans      = 22
 
    
    real(r8),dimension(1:nvegc,1:nvegc)                       :: AK_veg_acc
    real(r8),dimension(1:nvegn,1:nvegn)                       :: AK_veg_nacc
    real(r8),dimension(1:nvegc)                               :: I_veg_acc
    real(r8),dimension(1:nvegn)                               :: I_veg_nacc
    real(r8),dimension(1:ndecomp_pools_vr,1:ndecomp_pools_vr)       :: AK_soil_acc
    real(r8),dimension(1:ndecomp_pools_vr,1:ndecomp_pools_vr)       :: AK_soil_nacc
    real(r8),dimension(1:ndecomp_pools_vr)                          :: I_soil_acc
    real(r8),dimension(1:ndecomp_pools_vr)                          :: I_soil_nacc
    real(r8),dimension(1:nvegc,1:nvegc)                       :: AKinv_veg
    real(r8),dimension(1:nvegn,1:nvegn)                       :: AKinvn_veg
    real(r8),dimension(1:ndecomp_pools_vr,1:ndecomp_pools_vr)       :: AKinv_soil
    real(r8),dimension(1:ndecomp_pools_vr,1:ndecomp_pools_vr)       :: AKinvn_soil
    real(r8),dimension(1:nvegc,1)                             :: vegmatrixc_cap
    real(r8),dimension(1:nvegn,1)                             :: vegmatrixn_cap
    real(r8),dimension(1:ndecomp_pools_vr,1)                        :: soilmatrixc_cap
    real(r8),dimension(1:ndecomp_pools_vr,1)                        :: soilmatrixn_cap

      ! calculate c:n ratios of applicable pools
!      do l = 1, ndecomp_pools
!         if ( floating_cn_ratio(l)) then
!            do j = 1, nl_soil
!               if ( decomp_npools_vr(j,l,i) > 0._r8 ) then
!                  cn_decomp_pools(j,l,i) = decomp_cpools_vr(j,l,i) / decomp_npools_vr(j,l,i)
!               else
!                  cn_decomp_pools(j,l,i) = initial_cn_ratio(l)
!               end if
!            end do
!         else
!            do j = 1,nl_soil
!               cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
!            end do
!         end if 
!      end do

!      do j = 1, nl_soil
!         if(pmnf_decomp(j,1,i) .gt. 0)then
!            AKX_met_to_soil1_c_vr_acc(j,i)   = AKX_met_to_soil1_c_vr_acc  (j,i) + (1._r8 - rf_decomp(j,1,i)) * pathfrac_decomp(j,1,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(1),i) * decomp_k (j,donor_pool(1),i) * fpi_vr(j,i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(1)))then
!               AKX_met_to_soil1_n_vr_acc(j,i)   = AKX_met_to_soil1_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,1,i)) * pathfrac_decomp(j,1,i) &
!                                                * (cn_decomp_pools(j,donor_pool(1),i) / cn_decomp_pools(j,receiver_pool(1),i)) &
!                                                * decomp_last_npools_vr(j,donor_pool(1),i) * decomp_k (j,donor_pool(1),i) * fpi_vr(j,i) * deltim
!            else
!               AKX_met_to_soil1_n_vr_acc(j,i)   = AKX_met_to_soil1_n_vr_acc  (j,i) + pathfrac_decomp(j,1,i) &
!                                                * decomp_last_npools_vr(j,donor_pool(1),i) * decomp_k (j,donor_pool(1),i) * fpi_vr(j,i) * deltim
!            end if
!            AKX_met_exit_c_vr_acc(j,i) = AKX_met_exit_c_vr_acc(j,i) + pathfrac_decomp(j,1,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(1),i) * decomp_k (j,donor_pool(1),i) * fpi_vr(j,i) * deltim   
!            AKX_met_exit_n_vr_acc(j,i) = AKX_met_exit_n_vr_acc(j,i) + pathfrac_decomp(j,1,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(1),i) * decomp_k (j,donor_pool(1),i) * fpi_vr(j,i) * deltim   
!         else
!            AKX_met_to_soil1_c_vr_acc(j,i)   = AKX_met_to_soil1_c_vr_acc  (j,i) + (1._r8 - rf_decomp(j,1,i)) * pathfrac_decomp(j,1,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(1),i) * decomp_k (j,donor_pool(1),i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(1)))then
!               AKX_met_to_soil1_n_vr_acc(j,i)  = AKX_met_to_soil1_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,1,i)) * pathfrac_decomp(j,1,i) &
!                                               * (cn_decomp_pools(j,donor_pool(1),i) / cn_decomp_pools(j,receiver_pool(1),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(1),i) * decomp_k (j,donor_pool(1),i) * deltim
!            else
!               AKX_met_to_soil1_n_vr_acc(j,i)  = AKX_met_to_soil1_n_vr_acc  (j,i) + pathfrac_decomp(j,1,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(1),i) * decomp_k (j,donor_pool(1),i) * deltim
!            end if
!            AKX_met_exit_c_vr_acc(j,i) = AKX_met_exit_c_vr_acc(j,i) + pathfrac_decomp(j,1,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(1),i) * decomp_k (j,donor_pool(1),i) * deltim   
!            AKX_met_exit_n_vr_acc(j,i) = AKX_met_exit_n_vr_acc(j,i) + pathfrac_decomp(j,1,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(1),i) * decomp_k (j,donor_pool(1),i) * deltim   
!         end if
!         if(pmnf_decomp(j,2,i) .gt. 0)then
!            AKX_cel_to_soil1_c_vr_acc(j,i)   = AKX_cel_to_soil1_c_vr_acc  (j,i) + (1._r8 - rf_decomp(j,2,i)) * pathfrac_decomp(j,2,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(2),i) * decomp_k (j,donor_pool(2),i) * fpi_vr(j,i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(2)))then
!               AKX_cel_to_soil1_n_vr_acc(j,i)  = AKX_cel_to_soil1_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,2,i)) * pathfrac_decomp(j,2,i) &
!                                               * (cn_decomp_pools(j,donor_pool(2),i) / cn_decomp_pools(j,receiver_pool(2),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(2),i) * decomp_k (j,donor_pool(2),i) * fpi_vr(j,i) * deltim
!            else
!               AKX_cel_to_soil1_n_vr_acc(j,i)  = AKX_cel_to_soil1_n_vr_acc  (j,i) + pathfrac_decomp(j,2,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(2),i) * decomp_k (j,donor_pool(2),i) * fpi_vr(j,i) * deltim
!            end if
!            AKX_cel_exit_c_vr_acc(j,i) = AKX_cel_exit_c_vr_acc(j,i) + pathfrac_decomp(j,2,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(2),i) * decomp_k (j,donor_pool(2),i) * fpi_vr(j,i) * deltim   
!            AKX_cel_exit_n_vr_acc(j,i) = AKX_cel_exit_n_vr_acc(j,i) + pathfrac_decomp(j,2,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(2),i) * decomp_k (j,donor_pool(2),i) * fpi_vr(j,i) * deltim   
!         else
!            AKX_cel_to_soil1_c_vr_acc(j,i)   = AKX_cel_to_soil1_c_vr_acc  (j,i) + (1._r8 - rf_decomp(j,2,i)) * pathfrac_decomp(j,2,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(2),i) * decomp_k (j,donor_pool(2),i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(2)))then
!               AKX_cel_to_soil1_n_vr_acc(j,i)  = AKX_cel_to_soil1_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,2,i)) * pathfrac_decomp(j,2,i) &
!                                               * (cn_decomp_pools(j,donor_pool(2),i) / cn_decomp_pools(j,receiver_pool(2),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(2),i) * decomp_k (j,donor_pool(2),i) * deltim
!            else
!               AKX_cel_to_soil1_n_vr_acc(j,i)  = AKX_cel_to_soil1_n_vr_acc  (j,i) + pathfrac_decomp(j,2,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(2),i) * decomp_k (j,donor_pool(2),i) * deltim
!            end if
!            AKX_cel_exit_c_vr_acc(j,i) = AKX_cel_exit_c_vr_acc(j,i) + pathfrac_decomp(j,2,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(2),i) * decomp_k (j,donor_pool(2),i) * deltim   
!            AKX_cel_exit_n_vr_acc(j,i) = AKX_cel_exit_n_vr_acc(j,i) + pathfrac_decomp(j,2,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(2),i) * decomp_k (j,donor_pool(2),i) * deltim   
!         end if
!         if(pmnf_decomp(j,3,i) .gt. 0)then
!            AKX_lig_to_soil2_c_vr_acc(j,i)   = AKX_lig_to_soil2_c_vr_acc  (j,i) + (1._r8 - rf_decomp(j,3,i)) * pathfrac_decomp(j,3,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(3),i) * decomp_k (j,donor_pool(3),i) * fpi_vr(j,i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(3)))then
!               AKX_lig_to_soil2_n_vr_acc(j,i)  = AKX_lig_to_soil2_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,3,i)) * pathfrac_decomp(j,3,i) &
!                                               * (cn_decomp_pools(j,donor_pool(3),i) / cn_decomp_pools(j,receiver_pool(3),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(3),i) * decomp_k (j,donor_pool(3),i) * fpi_vr(j,i) * deltim
!            else
!               AKX_lig_to_soil2_n_vr_acc(j,i)  = AKX_lig_to_soil2_n_vr_acc  (j,i) + pathfrac_decomp(j,3,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(3),i) * decomp_k (j,donor_pool(3),i) * fpi_vr(j,i) * deltim
!            end if
!            AKX_lig_exit_c_vr_acc(j,i) = AKX_lig_exit_c_vr_acc(j,i) + pathfrac_decomp(j,3,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(3),i) * decomp_k (j,donor_pool(3),i) * fpi_vr(j,i) * deltim   
!            AKX_lig_exit_n_vr_acc(j,i) = AKX_lig_exit_n_vr_acc(j,i) + pathfrac_decomp(j,3,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(3),i) * decomp_k (j,donor_pool(3),i) * fpi_vr(j,i) * deltim   
!         else
!            AKX_lig_to_soil2_c_vr_acc(j,i)   = AKX_lig_to_soil2_c_vr_acc  (j,i) + (1._r8 - rf_decomp(j,3,i)) * pathfrac_decomp(j,3,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(3),i) * decomp_k (j,donor_pool(3),i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(3)))then
!               AKX_lig_to_soil2_n_vr_acc(j,i)  = AKX_lig_to_soil2_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,3,i)) * pathfrac_decomp(j,3,i) &
!                                               * (cn_decomp_pools(j,donor_pool(3),i) / cn_decomp_pools(j,receiver_pool(3),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(3),i) * decomp_k (j,donor_pool(3),i) * deltim
!            else
!               AKX_lig_to_soil2_n_vr_acc(j,i)  = AKX_lig_to_soil2_n_vr_acc  (j,i) + pathfrac_decomp(j,3,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(3),i) * decomp_k (j,donor_pool(3),i) * deltim
!            end if
!            AKX_lig_exit_c_vr_acc(j,i) = AKX_lig_exit_c_vr_acc(j,i) + pathfrac_decomp(j,3,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(3),i) * decomp_k (j,donor_pool(3),i) * deltim   
!            AKX_lig_exit_n_vr_acc(j,i) = AKX_lig_exit_n_vr_acc(j,i) + pathfrac_decomp(j,3,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(3),i) * decomp_k (j,donor_pool(3),i) * deltim   
!         end if
!         if(pmnf_decomp(j,4,i) .gt. 0)then
!            AKX_soil1_to_soil2_c_vr_acc(j,i) = AKX_soil1_to_soil2_c_vr_acc(j,i) + (1._r8 - rf_decomp(j,4,i)) * pathfrac_decomp(j,4,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(4),i) * decomp_k (j,donor_pool(4),i) * fpi_vr(j,i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(4)))then
!               AKX_soil1_to_soil2_n_vr_acc(j,i) = AKX_soil1_to_soil2_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,4,i)) * pathfrac_decomp(j,4,i) &
!                                                * (cn_decomp_pools(j,donor_pool(4),i) / cn_decomp_pools(j,receiver_pool(4),i)) &
!                                                * decomp_last_npools_vr(j,donor_pool(4),i) * decomp_k (j,donor_pool(4),i) * fpi_vr(j,i) * deltim
!            else
!               AKX_soil1_to_soil2_n_vr_acc(j,i) = AKX_soil1_to_soil2_n_vr_acc  (j,i) + pathfrac_decomp(j,4,i) &
!                                                * decomp_last_npools_vr(j,donor_pool(4),i) * decomp_k (j,donor_pool(4),i) * fpi_vr(j,i) * deltim
!            end if
!            AKX_soil1_exit_c_vr_acc(j,i) = AKX_soil1_exit_c_vr_acc(j,i) + pathfrac_decomp(j,4,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(4),i) * decomp_k (j,donor_pool(4),i) * fpi_vr(j,i) * deltim   
!            AKX_soil1_exit_n_vr_acc(j,i) = AKX_soil1_exit_n_vr_acc(j,i) + pathfrac_decomp(j,4,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(4),i) * decomp_k (j,donor_pool(4),i) * fpi_vr(j,i) * deltim   
!         else
!            AKX_soil1_to_soil2_c_vr_acc(j,i)    = AKX_soil1_to_soil2_c_vr_acc(j,i) + (1._r8 - rf_decomp(j,4,i)) * pathfrac_decomp(j,4,i) &
!                                                * decomp_last_cpools_vr(j,donor_pool(4),i) * decomp_k (j,donor_pool(4),i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(4)))then
!               AKX_soil1_to_soil2_n_vr_acc(j,i) = AKX_soil1_to_soil2_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,4,i)) * pathfrac_decomp(j,4,i) &
!                                                * (cn_decomp_pools(j,donor_pool(4),i) / cn_decomp_pools(j,receiver_pool(4),i)) &
!                                                * decomp_last_npools_vr(j,donor_pool(4),i) * decomp_k (j,donor_pool(4),i) * deltim
!            else
!               AKX_soil1_to_soil2_n_vr_acc(j,i)  = AKX_soil1_to_soil2_n_vr_acc  (j,i) + pathfrac_decomp(j,4,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(4),i) * decomp_k (j,donor_pool(4),i) * deltim
!            end if
!            AKX_soil1_exit_c_vr_acc(j,i) = AKX_soil1_exit_c_vr_acc(j,i) + pathfrac_decomp(j,4,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(4),i) * decomp_k (j,donor_pool(4),i) * deltim   
!            AKX_soil1_exit_n_vr_acc(j,i) = AKX_soil1_exit_n_vr_acc(j,i) + pathfrac_decomp(j,4,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(4),i) * decomp_k (j,donor_pool(4),i) * deltim   
!         end if
!         if(pmnf_decomp(j,5,i) .gt. 0)then
!            AKX_cwd_to_cel_c_vr_acc(j,i)     = AKX_cwd_to_cel_c_vr_acc    (j,i) + (1._r8 - rf_decomp(j,5,i)) * pathfrac_decomp(j,5,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(5),i) * decomp_k (j,donor_pool(5),i) * fpi_vr(j,i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(5)))then
!               AKX_cwd_to_cel_n_vr_acc(j,i)  = AKX_cwd_to_cel_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,5,i)) * pathfrac_decomp(j,5,i) &
!                                               * (cn_decomp_pools(j,donor_pool(5),i) / cn_decomp_pools(j,receiver_pool(5),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(5),i) * decomp_k (j,donor_pool(5),i) * fpi_vr(j,i) * deltim
!            else
!               AKX_cwd_to_cel_n_vr_acc(j,i)  = AKX_cwd_to_cel_n_vr_acc  (j,i) + pathfrac_decomp(j,5,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(5),i) * decomp_k (j,donor_pool(5),i) * fpi_vr(j,i) * deltim
!            end if
!            AKX_cwd_exit_c_vr_acc(j,i) = AKX_cwd_exit_c_vr_acc(j,i) + pathfrac_decomp(j,5,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(5),i) * decomp_k (j,donor_pool(5),i) * fpi_vr(j,i) * deltim   
!            AKX_cwd_exit_n_vr_acc(j,i) = AKX_cwd_exit_n_vr_acc(j,i) + pathfrac_decomp(j,5,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(5),i) * decomp_k (j,donor_pool(5),i) * fpi_vr(j,i) * deltim   
!         else
!            AKX_cwd_to_cel_c_vr_acc(j,i)     = AKX_cwd_to_cel_c_vr_acc    (j,i) + (1._r8 - rf_decomp(j,5,i)) * pathfrac_decomp(j,5,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(5),i) * decomp_k (j,donor_pool(5),i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(5)))then
!               AKX_cwd_to_cel_n_vr_acc(j,i)  = AKX_cwd_to_cel_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,5,i)) * pathfrac_decomp(j,5,i) &
!                                               * (cn_decomp_pools(j,donor_pool(5),i) / cn_decomp_pools(j,receiver_pool(5),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(5),i) * decomp_k (j,donor_pool(5),i) * deltim
!            else
!               AKX_cwd_to_cel_n_vr_acc(j,i)  = AKX_cwd_to_cel_n_vr_acc  (j,i) + pathfrac_decomp(j,5,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(5),i) * decomp_k (j,donor_pool(5),i) * deltim
!            end if
!            AKX_cwd_exit_c_vr_acc(j,i) = AKX_cwd_exit_c_vr_acc(j,i) + pathfrac_decomp(j,5,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(5),i) * decomp_k (j,donor_pool(5),i) * deltim   
!            AKX_cwd_exit_n_vr_acc(j,i) = AKX_cwd_exit_n_vr_acc(j,i) + pathfrac_decomp(j,5,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(5),i) * decomp_k (j,donor_pool(5),i) * deltim   
!         end if
!         if(pmnf_decomp(j,6,i) .gt. 0)then
!            AKX_cwd_to_lig_c_vr_acc(j,i)     = AKX_cwd_to_lig_c_vr_acc    (j,i) + (1._r8 - rf_decomp(j,6,i)) * pathfrac_decomp(j,6,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(6),i) * decomp_k (j,donor_pool(6),i) * fpi_vr(j,i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(6)))then
!               AKX_cwd_to_lig_n_vr_acc(j,i)  = AKX_cwd_to_lig_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,6,i)) * pathfrac_decomp(j,6,i) &
!                                               * (cn_decomp_pools(j,donor_pool(6),i) / cn_decomp_pools(j,receiver_pool(6),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(6),i) * decomp_k (j,donor_pool(6),i) * fpi_vr(j,i) * deltim
!            else
!               AKX_cwd_to_lig_n_vr_acc(j,i)  = AKX_cwd_to_lig_n_vr_acc  (j,i) + pathfrac_decomp(j,6,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(6),i) * decomp_k (j,donor_pool(6),i) * fpi_vr(j,i) * deltim
!            end if
!            AKX_cwd_exit_c_vr_acc(j,i) = AKX_cwd_exit_c_vr_acc(j,i) + pathfrac_decomp(j,6,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(6),i) * decomp_k (j,donor_pool(6),i) * fpi_vr(j,i) * deltim   
!            AKX_cwd_exit_n_vr_acc(j,i) = AKX_cwd_exit_n_vr_acc(j,i) + pathfrac_decomp(j,6,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(6),i) * decomp_k (j,donor_pool(6),i) * fpi_vr(j,i) * deltim   
!         else
!            AKX_cwd_to_lig_c_vr_acc(j,i)     = AKX_cwd_to_lig_c_vr_acc    (j,i) + (1._r8 - rf_decomp(j,6,i)) * pathfrac_decomp(j,6,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(6),i) * decomp_k (j,donor_pool(6),i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(6)))then
!               AKX_cwd_to_lig_n_vr_acc(j,i)  = AKX_cwd_to_lig_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,6,i)) * pathfrac_decomp(j,6,i) &
!                                               * (cn_decomp_pools(j,donor_pool(6),i) / cn_decomp_pools(j,receiver_pool(6),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(6),i) * decomp_k (j,donor_pool(6),i) * deltim
!            else
!               AKX_cwd_to_lig_n_vr_acc(j,i)  = AKX_cwd_to_lig_n_vr_acc  (j,i) + pathfrac_decomp(j,6,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(6),i) * decomp_k (j,donor_pool(6),i) * deltim
!            end if
!            AKX_cwd_exit_c_vr_acc(j,i) = AKX_cwd_exit_c_vr_acc(j,i) + pathfrac_decomp(j,6,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(6),i) * decomp_k (j,donor_pool(6),i) * deltim   
!            AKX_cwd_exit_n_vr_acc(j,i) = AKX_cwd_exit_n_vr_acc(j,i) + pathfrac_decomp(j,6,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(6),i) * decomp_k (j,donor_pool(6),i) * deltim   
!         end if
!         if(pmnf_decomp(j,7,i) .gt. 0)then
!            AKX_soil1_to_soil3_c_vr_acc(j,i) = AKX_soil1_to_soil3_c_vr_acc(j,i) + (1._r8 - rf_decomp(j,7,i)) * pathfrac_decomp(j,7,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(7),i) * decomp_k (j,donor_pool(7),i) * fpi_vr(j,i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(7)))then
!               AKX_soil1_to_soil3_n_vr_acc(j,i)  = AKX_soil1_to_soil3_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,7,i)) * pathfrac_decomp(j,7,i) &
!                                               * (cn_decomp_pools(j,donor_pool(7),i) / cn_decomp_pools(j,receiver_pool(7),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(7),i) * decomp_k (j,donor_pool(7),i) * fpi_vr(j,i) * deltim
!            else
!               AKX_soil1_to_soil3_n_vr_acc(j,i)  = AKX_soil1_to_soil3_n_vr_acc  (j,i) + pathfrac_decomp(j,7,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(7),i) * decomp_k (j,donor_pool(7),i) * fpi_vr(j,i) * deltim
!            end if
!            AKX_soil1_exit_c_vr_acc(j,i) = AKX_soil1_exit_c_vr_acc(j,i) + pathfrac_decomp(j,7,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(7),i) * decomp_k (j,donor_pool(7),i) * fpi_vr(j,i) * deltim   
!            AKX_soil1_exit_n_vr_acc(j,i) = AKX_soil1_exit_n_vr_acc(j,i) + pathfrac_decomp(j,7,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(7),i) * decomp_k (j,donor_pool(7),i) * fpi_vr(j,i) * deltim   
!         else
!            AKX_soil1_to_soil3_c_vr_acc(j,i) = AKX_soil1_to_soil3_c_vr_acc(j,i) + (1._r8 - rf_decomp(j,7,i)) * pathfrac_decomp(j,7,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(7),i) * decomp_k (j,donor_pool(7),i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(7)))then
!               AKX_soil1_to_soil3_n_vr_acc(j,i)  = AKX_soil1_to_soil3_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,7,i)) * pathfrac_decomp(j,7,i) &
!                                               * (cn_decomp_pools(j,donor_pool(7),i) / cn_decomp_pools(j,receiver_pool(7),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(7),i) * decomp_k (j,donor_pool(7),i) * deltim
!            else
!               AKX_soil1_to_soil3_n_vr_acc(j,i)  = AKX_soil1_to_soil3_n_vr_acc  (j,i) + pathfrac_decomp(j,7,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(7),i) * decomp_k (j,donor_pool(7),i) * deltim
!            end if
!            AKX_soil1_exit_c_vr_acc(j,i) = AKX_soil1_exit_c_vr_acc(j,i) + pathfrac_decomp(j,7,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(7),i) * decomp_k (j,donor_pool(7),i) * deltim   
!            AKX_soil1_exit_n_vr_acc(j,i) = AKX_soil1_exit_n_vr_acc(j,i) + pathfrac_decomp(j,7,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(7),i) * decomp_k (j,donor_pool(7),i) * deltim   
!         end if
!         if(pmnf_decomp(j,8,i) .gt. 0)then
!            AKX_soil2_to_soil1_c_vr_acc(j,i) = AKX_soil2_to_soil1_c_vr_acc(j,i) + (1._r8 - rf_decomp(j,8,i)) * pathfrac_decomp(j,8,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(8),i) * decomp_k (j,donor_pool(8),i) * fpi_vr(j,i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(8)))then
!               AKX_soil2_to_soil1_n_vr_acc(j,i) = AKX_soil2_to_soil1_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,8,i)) * pathfrac_decomp(j,8,i) &
!                                                * (cn_decomp_pools(j,donor_pool(8),i) / cn_decomp_pools(j,receiver_pool(8),i)) &
!                                                * decomp_last_npools_vr(j,donor_pool(8),i) * decomp_k (j,donor_pool(8),i) * fpi_vr(j,i) * deltim
!            else
!               AKX_soil2_to_soil1_n_vr_acc(j,i) = AKX_soil2_to_soil1_n_vr_acc  (j,i) + pathfrac_decomp(j,8,i) &
!                                                * decomp_last_npools_vr(j,donor_pool(8),i) * decomp_k (j,donor_pool(8),i) * fpi_vr(j,i) * deltim
!            end if
!            AKX_soil2_exit_c_vr_acc(j,i) = AKX_soil2_exit_c_vr_acc(j,i) + pathfrac_decomp(j,8,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(8),i) * decomp_k (j,donor_pool(8),i) * fpi_vr(j,i) * deltim   
!            AKX_soil2_exit_n_vr_acc(j,i) = AKX_soil2_exit_n_vr_acc(j,i) + pathfrac_decomp(j,8,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(8),i) * decomp_k (j,donor_pool(8),i) * fpi_vr(j,i) * deltim   
!         else
!            AKX_soil2_to_soil1_c_vr_acc(j,i) = AKX_soil2_to_soil1_c_vr_acc(j,i) + (1._r8 - rf_decomp(j,8,i)) * pathfrac_decomp(j,8,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(8),i) * decomp_k (j,donor_pool(8),i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(8)))then
!               AKX_soil2_to_soil1_n_vr_acc(j,i) = AKX_soil2_to_soil1_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,8,i)) * pathfrac_decomp(j,8,i) &
!                                                * (cn_decomp_pools(j,donor_pool(8),i) / cn_decomp_pools(j,receiver_pool(8),i)) &
!                                                * decomp_last_npools_vr(j,donor_pool(8),i) * decomp_k (j,donor_pool(8),i) * deltim
!            else
!               AKX_soil2_to_soil1_n_vr_acc(j,i) = AKX_soil2_to_soil1_n_vr_acc  (j,i) + pathfrac_decomp(j,8,i) &
!                                                * decomp_last_npools_vr(j,donor_pool(8),i) * decomp_k (j,donor_pool(8),i) * deltim
!            end if
!            AKX_soil2_exit_c_vr_acc(j,i) = AKX_soil2_exit_c_vr_acc(j,i) + pathfrac_decomp(j,8,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(8),i) * decomp_k (j,donor_pool(8),i) * deltim   
!            AKX_soil2_exit_n_vr_acc(j,i) = AKX_soil2_exit_n_vr_acc(j,i) + pathfrac_decomp(j,8,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(8),i) * decomp_k (j,donor_pool(8),i) * deltim   
!         end if
!         if(pmnf_decomp(j,9,i) .gt. 0)then
!            AKX_soil2_to_soil3_c_vr_acc(j,i) = AKX_soil2_to_soil3_c_vr_acc(j,i) + (1._r8 - rf_decomp(j,9,i)) * pathfrac_decomp(j,9,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(9),i) * decomp_k (j,donor_pool(9),i) * fpi_vr(j,i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(9)))then
!               AKX_soil2_to_soil3_n_vr_acc(j,i)  = AKX_soil2_to_soil3_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,9,i)) * pathfrac_decomp(j,9,i) &
!                                               * (cn_decomp_pools(j,donor_pool(9),i) / cn_decomp_pools(j,receiver_pool(9),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(9),i) * decomp_k (j,donor_pool(9),i) * fpi_vr(j,i) * deltim
!            else
!               AKX_soil2_to_soil3_n_vr_acc(j,i)  = AKX_soil2_to_soil3_n_vr_acc  (j,i) + pathfrac_decomp(j,9,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(9),i) * decomp_k (j,donor_pool(9),i) * fpi_vr(j,i) * deltim
!            end if
!            AKX_soil2_exit_c_vr_acc(j,i) = AKX_soil2_exit_c_vr_acc(j,i) + pathfrac_decomp(j,9,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(9),i) * decomp_k (j,donor_pool(9),i) * fpi_vr(j,i) * deltim   
!            AKX_soil2_exit_n_vr_acc(j,i) = AKX_soil2_exit_n_vr_acc(j,i) + pathfrac_decomp(j,9,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(9),i) * decomp_k (j,donor_pool(9),i) * fpi_vr(j,i) * deltim   
!         else
!            AKX_soil2_to_soil3_c_vr_acc(j,i) = AKX_soil2_to_soil3_c_vr_acc(j,i) + (1._r8 - rf_decomp(j,9,i)) * pathfrac_decomp(j,9,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(9),i) * decomp_k (j,donor_pool(9),i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(9)))then
!               AKX_soil2_to_soil3_n_vr_acc(j,i)  = AKX_soil2_to_soil3_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,9,i)) * pathfrac_decomp(j,9,i) &
!                                               * (cn_decomp_pools(j,donor_pool(9),i) / cn_decomp_pools(j,receiver_pool(9),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(9),i) * decomp_k (j,donor_pool(9),i) * deltim
!            else
!               AKX_soil2_to_soil3_n_vr_acc(j,i)  = AKX_soil2_to_soil3_n_vr_acc  (j,i) + pathfrac_decomp(j,9,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(9),i) * decomp_k (j,donor_pool(9),i) * deltim
!            end if
!            AKX_soil2_exit_c_vr_acc(j,i) = AKX_soil2_exit_c_vr_acc(j,i) + pathfrac_decomp(j,9,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(9),i) * decomp_k (j,donor_pool(9),i) * deltim   
!            AKX_soil2_exit_n_vr_acc(j,i) = AKX_soil2_exit_n_vr_acc(j,i) + pathfrac_decomp(j,9,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(9),i) * decomp_k (j,donor_pool(9),i) * deltim   
!         end if
!         if(pmnf_decomp(j,10,i) .gt. 0)then
!            AKX_soil3_to_soil1_c_vr_acc(j,i) = AKX_soil3_to_soil1_c_vr_acc(j,i) + (1._r8 - rf_decomp(j,10,i)) * pathfrac_decomp(j,10,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(10),i) * decomp_k (j,donor_pool(10),i) * fpi_vr(j,i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(10)))then
!               AKX_soil3_to_soil1_n_vr_acc(j,i)  = AKX_soil3_to_soil1_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,10,i)) * pathfrac_decomp(j,10,i) &
!                                               * (cn_decomp_pools(j,donor_pool(10),i) / cn_decomp_pools(j,receiver_pool(10),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(10),i) * decomp_k (j,donor_pool(10),i) * fpi_vr(j,i) * deltim
!            else
!               AKX_soil3_to_soil1_n_vr_acc(j,i)  = AKX_soil3_to_soil1_n_vr_acc  (j,i) + pathfrac_decomp(j,10,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(10),i) * decomp_k (j,donor_pool(10),i) * fpi_vr(j,i) * deltim
!            end if
!            AKX_soil3_exit_c_vr_acc(j,i) = AKX_soil3_exit_c_vr_acc(j,i) + pathfrac_decomp(j,10,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(10),i) * decomp_k (j,donor_pool(10),i) * fpi_vr(j,i) * deltim   
!            AKX_soil3_exit_n_vr_acc(j,i) = AKX_soil3_exit_n_vr_acc(j,i) + pathfrac_decomp(j,10,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(10),i) * decomp_k (j,donor_pool(10),i) * fpi_vr(j,i) * deltim   
!         else
!            AKX_soil3_to_soil1_c_vr_acc(j,i) = AKX_soil3_to_soil1_c_vr_acc(j,i) + (1._r8 - rf_decomp(j,10,i)) * pathfrac_decomp(j,10,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(10),i) * decomp_k (j,donor_pool(10),i) * deltim
!            if(.not. floating_cn_ratio(receiver_pool(10)))then
!               AKX_soil3_to_soil1_n_vr_acc(j,i)  = AKX_soil3_to_soil1_n_vr_acc  (j,i) + (1._r8 - rf_decomp(j,10,i)) * pathfrac_decomp(j,10,i) &
!                                               * (cn_decomp_pools(j,donor_pool(10),i) / cn_decomp_pools(j,receiver_pool(10),i)) &
!                                               * decomp_last_npools_vr(j,donor_pool(10),i) * decomp_k (j,donor_pool(10),i) * deltim
!            else
!               AKX_soil3_to_soil1_n_vr_acc(j,i)  = AKX_soil3_to_soil1_n_vr_acc  (j,i) + pathfrac_decomp(j,10,i) &
!                                               * decomp_last_npools_vr(j,donor_pool(10),i) * decomp_k (j,donor_pool(10),i) * deltim
!            end if
!            AKX_soil3_exit_c_vr_acc(j,i) = AKX_soil3_exit_c_vr_acc(j,i) + pathfrac_decomp(j,10,i) &
!                                             * decomp_last_cpools_vr(j,donor_pool(10),i) * decomp_k (j,donor_pool(10),i) * deltim   
!            AKX_soil3_exit_n_vr_acc(j,i) = AKX_soil3_exit_n_vr_acc(j,i) + pathfrac_decomp(j,10,i) &
!                                             * decomp_last_npools_vr(j,donor_pool(10),i) * decomp_k (j,donor_pool(10),i) * deltim   
!         end if
!
!      end do
      
                                
      ! Save the C and N pool size at begin of each year, which are used to calculate C and N capacity at end of each year.
      if (idate(2) .eq. 1 .and. idate(3) .eq. 1800)then  
         do m = ps, pe
            leafc0_p             (m) = max(leafc_p             (m),epsi)
            leafc0_storage_p     (m) = max(leafc_storage_p     (m),epsi)
            leafc0_xfer_p        (m) = max(leafc_xfer_p        (m),epsi)
            frootc0_p            (m) = max(frootc_p            (m),epsi)
            frootc0_storage_p    (m) = max(frootc_storage_p    (m),epsi)
            frootc0_xfer_p       (m) = max(frootc_xfer_p       (m),epsi)
            livestemc0_p         (m) = max(livestemc_p         (m),epsi)
            livestemc0_storage_p (m) = max(livestemc_storage_p (m),epsi)
            livestemc0_xfer_p    (m) = max(livestemc_xfer_p    (m),epsi)
            deadstemc0_p         (m) = max(deadstemc_p         (m),epsi)
            deadstemc0_storage_p (m) = max(deadstemc_storage_p (m),epsi)
            deadstemc0_xfer_p    (m) = max(deadstemc_xfer_p    (m),epsi)
            livecrootc0_p        (m) = max(livecrootc_p        (m),epsi)
            livecrootc0_storage_p(m) = max(livecrootc_storage_p(m),epsi)
            livecrootc0_xfer_p   (m) = max(livecrootc_xfer_p   (m),epsi)
            deadcrootc0_p        (m) = max(deadcrootc_p        (m),epsi)
            deadcrootc0_storage_p(m) = max(deadcrootc_storage_p(m),epsi)
            deadcrootc0_xfer_p   (m) = max(deadcrootc_xfer_p   (m),epsi)
            grainc0_p            (m) = max(grainc_p            (m),epsi)
            grainc0_storage_p    (m) = max(grainc_storage_p    (m),epsi)
            grainc0_xfer_p       (m) = max(grainc_xfer_p       (m),epsi)
            leafn0_p             (m) = max(leafn_p             (m),epsi)
            leafn0_storage_p     (m) = max(leafn_storage_p     (m),epsi)
            leafn0_xfer_p        (m) = max(leafn_xfer_p        (m),epsi)
            frootn0_p            (m) = max(frootn_p            (m),epsi)
            frootn0_storage_p    (m) = max(frootn_storage_p    (m),epsi)
            frootn0_xfer_p       (m) = max(frootn_xfer_p       (m),epsi)
            livestemn0_p         (m) = max(livestemn_p         (m),epsi)
            livestemn0_storage_p (m) = max(livestemn_storage_p (m),epsi)
            livestemn0_xfer_p    (m) = max(livestemn_xfer_p    (m),epsi)
            deadstemn0_p         (m) = max(deadstemn_p         (m),epsi)
            deadstemn0_storage_p (m) = max(deadstemn_storage_p (m),epsi)
            deadstemn0_xfer_p    (m) = max(deadstemn_xfer_p    (m),epsi)
            livecrootn0_p        (m) = max(livecrootn_p        (m),epsi)
            livecrootn0_storage_p(m) = max(livecrootn_storage_p(m),epsi)
            livecrootn0_xfer_p   (m) = max(livecrootn_xfer_p   (m),epsi)
            deadcrootn0_p        (m) = max(deadcrootn_p        (m),epsi)
            deadcrootn0_storage_p(m) = max(deadcrootn_storage_p(m),epsi)
            deadcrootn0_xfer_p   (m) = max(deadcrootn_xfer_p   (m),epsi)
            grainn0_p            (m) = max(grainn_p            (m),epsi)
            grainn0_storage_p    (m) = max(grainn_storage_p    (m),epsi)
            grainn0_xfer_p       (m) = max(grainn_xfer_p       (m),epsi)
            retransn0_p          (m) = max(retransn_p          (m),epsi)
         end do
         do k = 1, ndecomp_pools
            do j = 1, nl_soil
               decomp0_cpools_vr(j,k,i)=max(decomp_cpools_vr(j,k,i),epsi)
               decomp0_npools_vr(j,k,i)=max(decomp_npools_vr(j,k,i),epsi)
            end do
         end do
      end if

      if(idate(2) .eq. 365 .and. idate(3) .eq. 84600)then
         ! Copy C transfers from sparse matrix to 2D temporary variables tran_acc and tran_nacc
         ! Calculate the C and N transfer rate by dividing CN transfer by base value saved at begin of each year.
         
         print*,'SASU:update pool size',i
         do m = ps, pe
            AK_veg_acc  (1:nvegc,1:nvegc)                       = 0._r8
            AK_veg_nacc (1:nvegn,1:nvegn)                       = 0._r8
            I_veg_acc   (1:nvegc)                               = 0._r8
            I_veg_nacc  (1:nvegn)                               = 0._r8

            AK_veg_acc  (        ileaf,     ileaf_xf) = AKX_leafc_xf_to_leafc_p_acc             (m) / leafc0_xfer_p        (m)
            AK_veg_acc  (       ifroot,    ifroot_xf) = AKX_frootc_xf_to_frootc_p_acc           (m) / frootc0_xfer_p       (m)
            AK_veg_acc  (    ilivestem, ilivestem_xf) = AKX_livestemc_xf_to_livestemc_p_acc     (m) / livestemc0_xfer_p    (m)
            AK_veg_acc  (    ideadstem, ideadstem_xf) = AKX_deadstemc_xf_to_deadstemc_p_acc     (m) / deadstemc0_xfer_p    (m)
            AK_veg_acc  (   ilivecroot,ilivecroot_xf) = AKX_livecrootc_xf_to_livecrootc_p_acc   (m) / livecrootc0_xfer_p   (m)
            AK_veg_acc  (   ideadcroot,ideadcroot_xf) = AKX_deadcrootc_xf_to_deadcrootc_p_acc   (m) / deadcrootc0_xfer_p   (m)
            AK_veg_acc  (       igrain,    igrain_xf) = AKX_grainc_xf_to_grainc_p_acc           (m) / grainc0_xfer_p       (m)
            AK_veg_acc  (    ideadstem,    ilivestem) = AKX_livestemc_to_deadstemc_p_acc        (m) / livestemc0_p         (m)
            AK_veg_acc  (   ideadcroot,   ilivecroot) = AKX_livecrootc_to_deadcrootc_p_acc      (m) / livecrootc0_p        (m)

            AK_veg_acc  (     ileaf_xf,     ileaf_st) = AKX_leafc_st_to_leafc_xf_p_acc          (m) / leafc0_storage_p     (m)
            AK_veg_acc  (    ifroot_xf,    ifroot_st) = AKX_frootc_st_to_frootc_xf_p_acc        (m) / frootc0_storage_p    (m)
            AK_veg_acc  ( ilivestem_xf, ilivestem_st) = AKX_livestemc_st_to_livestemc_xf_p_acc  (m) / livestemc0_storage_p (m)
            AK_veg_acc  ( ideadstem_xf, ideadstem_st) = AKX_deadstemc_st_to_deadstemc_xf_p_acc  (m) / deadstemc0_storage_p (m)
            AK_veg_acc  (ilivecroot_xf,ilivecroot_st) = AKX_livecrootc_st_to_livecrootc_xf_p_acc(m) / livecrootc0_storage_p(m)
            AK_veg_acc  (ideadcroot_xf,ideadcroot_st) = AKX_deadcrootc_st_to_deadcrootc_xf_p_acc(m) / deadcrootc0_storage_p(m)
            AK_veg_acc  (    igrain_xf,    igrain_st) = AKX_grainc_st_to_grainc_xf_p_acc        (m) / grainc0_storage_p    (m)

            AK_veg_acc  (        ileaf,        ileaf) = - AKX_leafc_exit_p_acc                  (m) / leafc0_p             (m)
            AK_veg_acc  (     ileaf_st,     ileaf_st) = - AKX_leafc_st_exit_p_acc               (m) / leafc0_storage_p     (m)
            AK_veg_acc  (     ileaf_xf,     ileaf_xf) = - AKX_leafc_xf_exit_p_acc               (m) / leafc0_xfer_p        (m)
            AK_veg_acc  (       ifroot,       ifroot) = - AKX_frootc_exit_p_acc                 (m) / frootc0_p            (m)
            AK_veg_acc  (    ifroot_st,    ifroot_st) = - AKX_frootc_st_exit_p_acc              (m) / frootc0_storage_p    (m)
            AK_veg_acc  (    ifroot_xf,    ifroot_xf) = - AKX_frootc_xf_exit_p_acc              (m) / frootc0_xfer_p       (m)
            AK_veg_acc  (    ilivestem,    ilivestem) = - AKX_livestemc_exit_p_acc              (m) / livestemc0_p         (m)
            AK_veg_acc  ( ilivestem_st, ilivestem_st) = - AKX_livestemc_st_exit_p_acc           (m) / livestemc0_storage_p (m)
            AK_veg_acc  ( ilivestem_xf, ilivestem_xf) = - AKX_livestemc_xf_exit_p_acc           (m) / livestemc0_xfer_p    (m)
            AK_veg_acc  (    ideadstem,    ideadstem) = - AKX_deadstemc_exit_p_acc              (m) / deadstemc0_p         (m)
            AK_veg_acc  ( ideadstem_st, ideadstem_st) = - AKX_deadstemc_st_exit_p_acc           (m) / deadstemc0_storage_p (m)
            AK_veg_acc  ( ideadstem_xf, ideadstem_xf) = - AKX_deadstemc_xf_exit_p_acc           (m) / deadstemc0_xfer_p    (m)
            AK_veg_acc  (   ilivecroot,   ilivecroot) = - AKX_livecrootc_exit_p_acc             (m) / livecrootc0_p        (m)
            AK_veg_acc  (ilivecroot_st,ilivecroot_st) = - AKX_livecrootc_st_exit_p_acc          (m) / livecrootc0_storage_p(m)
            AK_veg_acc  (ilivecroot_xf,ilivecroot_xf) = - AKX_livecrootc_xf_exit_p_acc          (m) / livecrootc0_xfer_p   (m)
            AK_veg_acc  (   ideadcroot,   ideadcroot) = - AKX_deadcrootc_exit_p_acc             (m) / deadcrootc0_p        (m)
            AK_veg_acc  (ideadcroot_st,ideadcroot_st) = - AKX_deadcrootc_st_exit_p_acc          (m) / deadcrootc0_storage_p(m)
            AK_veg_acc  (ideadcroot_xf,ideadcroot_xf) = - AKX_deadcrootc_xf_exit_p_acc          (m) / deadcrootc0_xfer_p   (m)
            AK_veg_acc  (       igrain,       igrain) = - AKX_grainc_exit_p_acc                 (m) / grainc0_p            (m)
            AK_veg_acc  (    igrain_st,    igrain_st) = - AKX_grainc_st_exit_p_acc              (m) / grainc0_storage_p    (m)
            AK_veg_acc  (    igrain_xf,    igrain_xf) = - AKX_grainc_xf_exit_p_acc              (m) / grainc0_xfer_p       (m)

            I_veg_acc   (        ileaf) = I_leafc_p_acc         (m)
            I_veg_acc   (     ileaf_st) = I_leafc_st_p_acc      (m)
            I_veg_acc   (       ifroot) = I_frootc_p_acc        (m)
            I_veg_acc   (    ifroot_st) = I_frootc_st_p_acc     (m)
            I_veg_acc   (    ilivestem) = I_livestemc_p_acc     (m)
            I_veg_acc   ( ilivestem_st) = I_livestemc_st_p_acc  (m)
            I_veg_acc   (    ideadstem) = I_deadstemc_p_acc     (m)
            I_veg_acc   ( ideadstem_st) = I_deadstemc_st_p_acc  (m)
            I_veg_acc   (   ilivecroot) = I_livecrootc_p_acc    (m)
            I_veg_acc   (ilivecroot_st) = I_livecrootc_st_p_acc (m)
            I_veg_acc   (   ideadcroot) = I_deadcrootc_p_acc    (m)
            I_veg_acc   (ideadcroot_st) = I_deadcrootc_st_p_acc (m)
            I_veg_acc   (       igrain) = I_grainc_p_acc        (m)
            I_veg_acc   (    igrain_st) = I_grainc_st_p_acc     (m)

            AK_veg_nacc (        ileaf,     ileaf_xf) = AKX_leafn_xf_to_leafn_p_acc             (m) / leafn0_xfer_p        (m)
            AK_veg_nacc (       ifroot,    ifroot_xf) = AKX_frootn_xf_to_frootn_p_acc           (m) / frootn0_xfer_p       (m)
            AK_veg_nacc (    ilivestem, ilivestem_xf) = AKX_livestemn_xf_to_livestemn_p_acc     (m) / livestemn0_xfer_p    (m)
            AK_veg_nacc (    ideadstem, ideadstem_xf) = AKX_deadstemn_xf_to_deadstemn_p_acc     (m) / deadstemn0_xfer_p    (m)
            AK_veg_nacc (   ilivecroot,ilivecroot_xf) = AKX_livecrootn_xf_to_livecrootn_p_acc   (m) / livecrootn0_xfer_p   (m)
            AK_veg_nacc (   ideadcroot,ideadcroot_xf) = AKX_deadcrootn_xf_to_deadcrootn_p_acc   (m) / deadcrootn0_xfer_p   (m)
            AK_veg_nacc (       igrain,    igrain_xf) = AKX_grainn_xf_to_grainn_p_acc           (m) / grainn0_xfer_p       (m)
            AK_veg_nacc (    ideadstem,    ilivestem) = AKX_livestemn_to_deadstemn_p_acc        (m) / livestemn0_p         (m)
            AK_veg_nacc (   ideadcroot,   ilivecroot) = AKX_livecrootn_to_deadcrootn_p_acc      (m) / livecrootn0_p        (m)

            AK_veg_nacc (     ileaf_xf,     ileaf_st) = AKX_leafn_st_to_leafn_xf_p_acc          (m) / leafn0_storage_p     (m)
            AK_veg_nacc (    ifroot_xf,    ifroot_st) = AKX_frootn_st_to_frootn_xf_p_acc        (m) / frootn0_storage_p    (m)
            AK_veg_nacc ( ilivestem_xf, ilivestem_st) = AKX_livestemn_st_to_livestemn_xf_p_acc  (m) / livestemn0_storage_p (m)
            AK_veg_nacc ( ideadstem_xf, ideadstem_st) = AKX_deadstemn_st_to_deadstemn_xf_p_acc  (m) / deadstemn0_storage_p (m)
            AK_veg_nacc (ilivecroot_xf,ilivecroot_st) = AKX_livecrootn_st_to_livecrootn_xf_p_acc(m) / livecrootn0_storage_p(m)
            AK_veg_nacc (ideadcroot_xf,ideadcroot_st) = AKX_deadcrootn_st_to_deadcrootn_xf_p_acc(m) / deadcrootn0_storage_p(m)
            AK_veg_nacc (    igrain_xf,    igrain_st) = AKX_grainn_st_to_grainn_xf_p_acc        (m) / grainn0_storage_p    (m)

            AK_veg_nacc (     iretrans,        ileaf) = AKX_leafn_to_retransn_p_acc             (m) / leafn0_p             (m)
            AK_veg_nacc (     iretrans,       ifroot) = AKX_frootn_to_retransn_p_acc            (m) / frootn0_p            (m)
            AK_veg_nacc (     iretrans,    ilivestem) = AKX_livestemn_to_retransn_p_acc         (m) / livestemn0_p         (m)
            AK_veg_nacc (     iretrans,   ilivecroot) = AKX_livecrootn_to_retransn_p_acc        (m) / livecrootn0_p        (m)

            AK_veg_nacc (        ileaf,     iretrans) = AKX_retransn_to_leafn_p_acc             (m) / retransn0_p          (m)
            AK_veg_nacc (       ifroot,     iretrans) = AKX_retransn_to_frootn_p_acc            (m) / retransn0_p          (m)
            AK_veg_nacc (    ilivestem,     iretrans) = AKX_retransn_to_livestemn_p_acc         (m) / retransn0_p          (m)
            AK_veg_nacc (    ideadstem,     iretrans) = AKX_retransn_to_deadstemn_p_acc         (m) / retransn0_p          (m)
            AK_veg_nacc (   ilivecroot,     iretrans) = AKX_retransn_to_livecrootn_p_acc        (m) / retransn0_p          (m)
            AK_veg_nacc (   ideadcroot,     iretrans) = AKX_retransn_to_deadcrootn_p_acc        (m) / retransn0_p          (m)
            AK_veg_nacc (       igrain,     iretrans) = AKX_retransn_to_grainn_p_acc            (m) / retransn0_p          (m)

            AK_veg_nacc (     ileaf_st,     iretrans) = AKX_retransn_to_leafn_st_p_acc          (m) / retransn0_p          (m)
            AK_veg_nacc (    ifroot_st,     iretrans) = AKX_retransn_to_frootn_st_p_acc         (m) / retransn0_p          (m)
            AK_veg_nacc ( ilivestem_st,     iretrans) = AKX_retransn_to_livestemn_st_p_acc      (m) / retransn0_p          (m)
            AK_veg_nacc ( ideadstem_st,     iretrans) = AKX_retransn_to_deadstemn_st_p_acc      (m) / retransn0_p          (m)
            AK_veg_nacc (ilivecroot_st,     iretrans) = AKX_retransn_to_livecrootn_st_p_acc     (m) / retransn0_p          (m)
            AK_veg_nacc (ideadcroot_st,     iretrans) = AKX_retransn_to_deadcrootn_st_p_acc     (m) / retransn0_p          (m)
            AK_veg_nacc (    igrain_st,     iretrans) = AKX_retransn_to_grainn_st_p_acc         (m) / retransn0_p          (m)

            AK_veg_nacc (        ileaf,        ileaf) = - AKX_leafn_exit_p_acc                  (m) / leafn0_p             (m)
            AK_veg_nacc (     ileaf_st,     ileaf_st) = - AKX_leafn_st_exit_p_acc               (m) / leafn0_storage_p     (m)
            AK_veg_nacc (     ileaf_xf,     ileaf_xf) = - AKX_leafn_xf_exit_p_acc               (m) / leafn0_xfer_p        (m)
            AK_veg_nacc (       ifroot,       ifroot) = - AKX_frootn_exit_p_acc                 (m) / frootn0_p            (m)
            AK_veg_nacc (    ifroot_st,    ifroot_st) = - AKX_frootn_st_exit_p_acc              (m) / frootn0_storage_p    (m)
            AK_veg_nacc (    ifroot_xf,    ifroot_xf) = - AKX_frootn_xf_exit_p_acc              (m) / frootn0_xfer_p       (m)
            AK_veg_nacc (    ilivestem,    ilivestem) = - AKX_livestemn_exit_p_acc              (m) / livestemn0_p         (m)
            AK_veg_nacc ( ilivestem_st, ilivestem_st) = - AKX_livestemn_st_exit_p_acc           (m) / livestemn0_storage_p (m)
            AK_veg_nacc ( ilivestem_xf, ilivestem_xf) = - AKX_livestemn_xf_exit_p_acc           (m) / livestemn0_xfer_p    (m)
            AK_veg_nacc (    ideadstem,    ideadstem) = - AKX_deadstemn_exit_p_acc              (m) / deadstemn0_p         (m)
            AK_veg_nacc ( ideadstem_st, ideadstem_st) = - AKX_deadstemn_st_exit_p_acc           (m) / deadstemn0_storage_p (m)
            AK_veg_nacc ( ideadstem_xf, ideadstem_xf) = - AKX_deadstemn_xf_exit_p_acc           (m) / deadstemn0_xfer_p    (m)
            AK_veg_nacc (   ilivecroot,   ilivecroot) = - AKX_livecrootn_exit_p_acc             (m) / livecrootn0_p        (m)
            AK_veg_nacc (ilivecroot_st,ilivecroot_st) = - AKX_livecrootn_st_exit_p_acc          (m) / livecrootn0_storage_p(m)
            AK_veg_nacc (ilivecroot_xf,ilivecroot_xf) = - AKX_livecrootn_xf_exit_p_acc          (m) / livecrootn0_xfer_p   (m)
            AK_veg_nacc (   ideadcroot,   ideadcroot) = - AKX_deadcrootn_exit_p_acc             (m) / deadcrootn0_p        (m)
            AK_veg_nacc (ideadcroot_st,ideadcroot_st) = - AKX_deadcrootn_st_exit_p_acc          (m) / deadcrootn0_storage_p(m)
            AK_veg_nacc (ideadcroot_xf,ideadcroot_xf) = - AKX_deadcrootn_xf_exit_p_acc          (m) / deadcrootn0_xfer_p   (m)
            AK_veg_nacc (       igrain,       igrain) = - AKX_grainn_exit_p_acc                 (m) / grainn0_p            (m)
            AK_veg_nacc (    igrain_st,    igrain_st) = - AKX_grainn_st_exit_p_acc              (m) / grainn0_storage_p    (m)
            AK_veg_nacc (    igrain_xf,    igrain_xf) = - AKX_grainn_xf_exit_p_acc              (m) / grainn0_xfer_p       (m)
            AK_veg_nacc (     iretrans,     iretrans) = - AKX_retransn_exit_p_acc               (m) / retransn0_p          (m)

            I_veg_nacc  (        ileaf) = I_leafn_p_acc         (m)
            I_veg_nacc  (     ileaf_st) = I_leafn_st_p_acc      (m)
            I_veg_nacc  (       ifroot) = I_frootn_p_acc        (m)
            I_veg_nacc  (    ifroot_st) = I_frootn_st_p_acc     (m)
            I_veg_nacc  (    ilivestem) = I_livestemn_p_acc     (m)
            I_veg_nacc  ( ilivestem_st) = I_livestemn_st_p_acc  (m)
            I_veg_nacc  (    ideadstem) = I_deadstemn_p_acc     (m)
            I_veg_nacc  ( ideadstem_st) = I_deadstemn_st_p_acc  (m)
            I_veg_nacc  (   ilivecroot) = I_livecrootn_p_acc    (m)
            I_veg_nacc  (ilivecroot_st) = I_livecrootn_st_p_acc (m)
            I_veg_nacc  (   ideadcroot) = I_deadcrootn_p_acc    (m)
            I_veg_nacc  (ideadcroot_st) = I_deadcrootn_st_p_acc (m)
            I_veg_nacc  (       igrain) = I_grainn_p_acc        (m)
            I_veg_nacc  (    igrain_st) = I_grainn_st_p_acc     (m)

            do j = 1, nvegc
               if(AK_veg_acc(j,j) .eq. 0)then
                  AK_veg_acc(j,j) = - 1.e+36
               end if
            end do
            do j = 1, nvegn
               if(AK_veg_nacc(j,j) .eq. 0)then
                  AK_veg_nacc(j,j) = - 1.e+36
               end if
            end do

            ! Calculate capacity
            call inverse(AK_veg_acc (1:nvegc,1:nvegc),AKinv_veg (1:nvegc,1:nvegc),nvegc)
            call inverse(AK_veg_nacc(1:nvegn,1:nvegn),AKinvn_veg(1:nvegn,1:nvegn),nvegn)
            vegmatrixc_cap(:,1) = -matmul(AKinv_veg (1:nvegc,1:nvegc),I_veg_acc (1:nvegc))
            vegmatrixn_cap(:,1) = -matmul(AKinvn_veg(1:nvegn,1:nvegn),I_veg_nacc(1:nvegn))
  
            do k = 1, nvegc
               if(vegmatrixc_cap(k,1) .lt. 0)then
                  vegmatrixc_cap(k,1) = epsi
               endif
            end do
            do k = 1, nvegn
               if(vegmatrixn_cap(k,1) .lt. 0)then
                  vegmatrixn_cap(k,1) = epsi
               endif
            end do
            deadstemc_p         (m) = vegmatrixc_cap(ideadstem,1)
            deadstemc_storage_p (m) = vegmatrixc_cap(ideadstem_st,1)
            deadcrootc_p        (m) = vegmatrixc_cap(ideadcroot,1)
            deadcrootc_storage_p(m) = vegmatrixc_cap(ideadcroot_st,1)
            deadstemn_p         (m) = vegmatrixn_cap(ideadstem,1)
            deadstemn_storage_p (m) = vegmatrixn_cap(ideadstem_st,1)
            deadcrootn_p        (m) = vegmatrixn_cap(ideadcroot,1)
            deadcrootn_storage_p(m) = vegmatrixn_cap(ideadcroot_st,1)
         end do

         AK_soil_acc (1:ndecomp_pools_vr,1:ndecomp_pools_vr) = 0._r8
         AK_soil_nacc(1:ndecomp_pools_vr,1:ndecomp_pools_vr) = 0._r8
         I_soil_acc  (1:ndecomp_pools_vr)                    = 0._r8
         I_soil_nacc (1:ndecomp_pools_vr)                    = 0._r8
         do j=1, nl_soil
            ! C exit rate
            AK_soil_acc ((i_met_lit-1)*nl_soil+j,(i_met_lit-1)*nl_soil+j) &
                = - (AKX_met_exit_c_vr_acc(j,i) + diagVX_c_vr_acc(j,i_met_lit,i)) / decomp0_cpools_vr(j,i_met_lit,i)
            AK_soil_acc ((i_cel_lit-1)*nl_soil+j,(i_cel_lit-1)*nl_soil+j) &
                = - (AKX_cel_exit_c_vr_acc(j,i) + diagVX_c_vr_acc(j,i_cel_lit,i)) / decomp0_cpools_vr(j,i_cel_lit,i)
            AK_soil_acc ((i_lig_lit-1)*nl_soil+j,(i_lig_lit-1)*nl_soil+j) &
                = - (AKX_lig_exit_c_vr_acc(j,i) + diagVX_c_vr_acc(j,i_lig_lit,i)) / decomp0_cpools_vr(j,i_lig_lit,i)
            AK_soil_acc ((i_cwd    -1)*nl_soil+j,(i_cwd    -1)*nl_soil+j) &
                = -  AKX_cwd_exit_c_vr_acc(j,i)                                   / decomp0_cpools_vr(j,i_cwd    ,i)
            AK_soil_acc ((i_soil1  -1)*nl_soil+j,(i_soil1  -1)*nl_soil+j) &
                = - (AKX_soil1_exit_c_vr_acc(j,i) + diagVX_c_vr_acc(j,i_soil1,i)) / decomp0_cpools_vr(j,i_soil1  ,i)
            AK_soil_acc ((i_soil2  -1)*nl_soil+j,(i_soil2  -1)*nl_soil+j) &
                = - (AKX_soil2_exit_c_vr_acc(j,i) + diagVX_c_vr_acc(j,i_soil2,i)) / decomp0_cpools_vr(j,i_soil2  ,i)
            AK_soil_acc ((i_soil3  -1)*nl_soil+j,(i_soil3  -1)*nl_soil+j) &
                = - (AKX_soil3_exit_c_vr_acc(j,i) + diagVX_c_vr_acc(j,i_soil3,i)) / decomp0_cpools_vr(j,i_soil3  ,i)

            ! C transfer
            AK_soil_acc ((i_soil1  -1)*nl_soil+j,(i_met_lit-1)*nl_soil+j) &
                = AKX_met_to_soil1_c_vr_acc  (j,i) / decomp0_cpools_vr(j,i_met_lit,i)
            AK_soil_acc ((i_soil1  -1)*nl_soil+j,(i_cel_lit-1)*nl_soil+j) &
                = AKX_cel_to_soil1_c_vr_acc  (j,i) / decomp0_cpools_vr(j,i_cel_lit,i)
            AK_soil_acc ((i_soil2  -1)*nl_soil+j,(i_lig_lit-1)*nl_soil+j) &
                = AKX_lig_to_soil2_c_vr_acc  (j,i) / decomp0_cpools_vr(j,i_lig_lit,i)
            AK_soil_acc ((i_soil2  -1)*nl_soil+j,(i_soil1  -1)*nl_soil+j) &
                = AKX_soil1_to_soil2_c_vr_acc(j,i) / decomp0_cpools_vr(j,i_soil1  ,i)
            AK_soil_acc ((i_cel_lit-1)*nl_soil+j,(i_cwd    -1)*nl_soil+j) &
                = AKX_cwd_to_cel_c_vr_acc    (j,i) / decomp0_cpools_vr(j,i_cwd    ,i)
            AK_soil_acc ((i_lig_lit-1)*nl_soil+j,(i_cwd    -1)*nl_soil+j) &
                = AKX_cwd_to_lig_c_vr_acc    (j,i) / decomp0_cpools_vr(j,i_cwd    ,i)
            AK_soil_acc ((i_soil3  -1)*nl_soil+j,(i_soil1  -1)*nl_soil+j) &
                = AKX_soil1_to_soil3_c_vr_acc(j,i) / decomp0_cpools_vr(j,i_soil1  ,i)
            AK_soil_acc ((i_soil1  -1)*nl_soil+j,(i_soil2  -1)*nl_soil+j) &
                = AKX_soil2_to_soil1_c_vr_acc(j,i) / decomp0_cpools_vr(j,i_soil2  ,i)
            AK_soil_acc ((i_soil3  -1)*nl_soil+j,(i_soil2  -1)*nl_soil+j) &
                = AKX_soil2_to_soil3_c_vr_acc(j,i) / decomp0_cpools_vr(j,i_soil2  ,i)
            AK_soil_acc ((i_soil1  -1)*nl_soil+j,(i_soil3  -1)*nl_soil+j) &
                = AKX_soil3_to_soil1_c_vr_acc(j,i) / decomp0_cpools_vr(j,i_soil3  ,i)

            ! C input
            I_soil_acc((i_met_lit-1)*nl_soil+j) = I_met_c_vr_acc(j,i)
            I_soil_acc((i_cel_lit-1)*nl_soil+j) = I_cel_c_vr_acc(j,i)
            I_soil_acc((i_lig_lit-1)*nl_soil+j) = I_lig_c_vr_acc(j,i)
            I_soil_acc((i_cwd    -1)*nl_soil+j) = I_cwd_c_vr_acc(j,i)

            ! N exit rate
            AK_soil_nacc((i_met_lit-1)*nl_soil+j,(i_met_lit-1)*nl_soil+j) &
                = - (AKX_met_exit_n_vr_acc(j,i) + diagVX_n_vr_acc(j,i_met_lit,i)) / decomp0_npools_vr(j,i_met_lit,i)
            AK_soil_nacc((i_cel_lit-1)*nl_soil+j,(i_cel_lit-1)*nl_soil+j) &
                = - (AKX_cel_exit_n_vr_acc(j,i) + diagVX_n_vr_acc(j,i_cel_lit,i)) / decomp0_npools_vr(j,i_cel_lit,i)
            AK_soil_nacc((i_lig_lit-1)*nl_soil+j,(i_lig_lit-1)*nl_soil+j) &
                = - (AKX_lig_exit_n_vr_acc(j,i) + diagVX_n_vr_acc(j,i_lig_lit,i)) / decomp0_npools_vr(j,i_lig_lit,i)
            AK_soil_nacc((i_cwd    -1)*nl_soil+j,(i_cwd    -1)*nl_soil+j) &
                = -  AKX_cwd_exit_n_vr_acc(j,i)                                   / decomp0_npools_vr(j,i_cwd    ,i)
            AK_soil_nacc((i_soil1  -1)*nl_soil+j,(i_soil1  -1)*nl_soil+j) &
                = - (AKX_soil1_exit_n_vr_acc(j,i) + diagVX_n_vr_acc(j,i_soil1,i)) / decomp0_npools_vr(j,i_soil1  ,i)
            AK_soil_nacc((i_soil2  -1)*nl_soil+j,(i_soil2  -1)*nl_soil+j) &
                = - (AKX_soil2_exit_n_vr_acc(j,i) + diagVX_n_vr_acc(j,i_soil2,i)) / decomp0_npools_vr(j,i_soil2  ,i)
            AK_soil_nacc((i_soil3  -1)*nl_soil+j,(i_soil3  -1)*nl_soil+j) &
                = - (AKX_soil3_exit_n_vr_acc(j,i) + diagVX_n_vr_acc(j,i_soil3,i)) / decomp0_npools_vr(j,i_soil3  ,i)

            ! N transfer
            AK_soil_nacc((i_soil1  -1)*nl_soil+j,(i_met_lit-1)*nl_soil+j) &
                = AKX_met_to_soil1_n_vr_acc  (j,i) / decomp0_npools_vr(j,i_met_lit,i)
            AK_soil_nacc((i_soil1  -1)*nl_soil+j,(i_cel_lit-1)*nl_soil+j) &
                = AKX_cel_to_soil1_n_vr_acc  (j,i) / decomp0_npools_vr(j,i_cel_lit,i)
            AK_soil_nacc((i_soil2  -1)*nl_soil+j,(i_lig_lit-1)*nl_soil+j) &
                = AKX_lig_to_soil2_n_vr_acc  (j,i) / decomp0_npools_vr(j,i_lig_lit,i)
            AK_soil_nacc((i_soil2  -1)*nl_soil+j,(i_soil1  -1)*nl_soil+j) &
                = AKX_soil1_to_soil2_n_vr_acc(j,i) / decomp0_npools_vr(j,i_soil1  ,i)
            AK_soil_nacc((i_cel_lit-1)*nl_soil+j,(i_cwd    -1)*nl_soil+j) &
                = AKX_cwd_to_cel_n_vr_acc    (j,i) / decomp0_npools_vr(j,i_cwd    ,i)
            AK_soil_nacc((i_lig_lit-1)*nl_soil+j,(i_cwd    -1)*nl_soil+j) &
                = AKX_cwd_to_lig_n_vr_acc    (j,i) / decomp0_npools_vr(j,i_cwd    ,i)
            AK_soil_nacc((i_soil3  -1)*nl_soil+j,(i_soil1  -1)*nl_soil+j) &
                = AKX_soil1_to_soil3_n_vr_acc(j,i) / decomp0_npools_vr(j,i_soil1  ,i)
            AK_soil_nacc((i_soil1  -1)*nl_soil+j,(i_soil2  -1)*nl_soil+j) &
                = AKX_soil2_to_soil1_n_vr_acc(j,i) / decomp0_npools_vr(j,i_soil2  ,i)
            AK_soil_nacc((i_soil3  -1)*nl_soil+j,(i_soil2  -1)*nl_soil+j) &
                = AKX_soil2_to_soil3_n_vr_acc(j,i) / decomp0_npools_vr(j,i_soil2  ,i)
            AK_soil_nacc((i_soil1  -1)*nl_soil+j,(i_soil3  -1)*nl_soil+j) &
                = AKX_soil3_to_soil1_n_vr_acc(j,i) / decomp0_npools_vr(j,i_soil3  ,i)

         end do

         do j=1,nl_soil-1
            ! upper triadiagnonal entries for C
            AK_soil_acc ((i_met_lit-1)*nl_soil+j,(i_met_lit-1)*nl_soil+j+1) &
                = upperVX_c_vr_acc(j,i_met_lit,i) / decomp0_cpools_vr(j+1,i_met_lit,i)
            AK_soil_acc ((i_cel_lit-1)*nl_soil+j,(i_cel_lit-1)*nl_soil+j+1) &
                = upperVX_c_vr_acc(j,i_cel_lit,i) / decomp0_cpools_vr(j+1,i_cel_lit,i)
            AK_soil_acc ((i_lig_lit-1)*nl_soil+j,(i_lig_lit-1)*nl_soil+j+1) &
                = upperVX_c_vr_acc(j,i_lig_lit,i) / decomp0_cpools_vr(j+1,i_lig_lit,i)
            AK_soil_acc ((i_soil1  -1)*nl_soil+j,(i_soil1  -1)*nl_soil+j+1) &
                = upperVX_c_vr_acc(j,i_soil1  ,i) / decomp0_cpools_vr(j+1,i_soil1  ,i)
            AK_soil_acc ((i_soil2  -1)*nl_soil+j,(i_soil2  -1)*nl_soil+j+1) &
                = upperVX_c_vr_acc(j,i_soil2  ,i) / decomp0_cpools_vr(j+1,i_soil2  ,i)
            AK_soil_acc ((i_soil3  -1)*nl_soil+j,(i_soil3  -1)*nl_soil+j+1) &
                = upperVX_c_vr_acc(j,i_soil3  ,i) / decomp0_cpools_vr(j+1,i_soil3  ,i)

            ! lower triadiagnonal entries for C
            AK_soil_acc ((i_met_lit-1)*nl_soil+j+1,(i_met_lit-1)*nl_soil+j) &
                = lowerVX_c_vr_acc(j+1,i_met_lit,i) / decomp0_cpools_vr(j,i_met_lit,i)
            AK_soil_acc ((i_cel_lit-1)*nl_soil+j+1,(i_cel_lit-1)*nl_soil+j) &
                = lowerVX_c_vr_acc(j+1,i_cel_lit,i) / decomp0_cpools_vr(j,i_cel_lit,i)
            AK_soil_acc ((i_lig_lit-1)*nl_soil+j+1,(i_lig_lit-1)*nl_soil+j) &
                = lowerVX_c_vr_acc(j+1,i_lig_lit,i) / decomp0_cpools_vr(j,i_lig_lit,i)
            AK_soil_acc ((i_soil1  -1)*nl_soil+j+1,(i_soil1  -1)*nl_soil+j) &
                = lowerVX_c_vr_acc(j+1,i_soil1  ,i) / decomp0_cpools_vr(j,i_soil1  ,i)
            AK_soil_acc ((i_soil2  -1)*nl_soil+j+1,(i_soil2  -1)*nl_soil+j) &
                = lowerVX_c_vr_acc(j+1,i_soil2  ,i) / decomp0_cpools_vr(j,i_soil2  ,i)
            AK_soil_acc ((i_soil3  -1)*nl_soil+j+1,(i_soil3  -1)*nl_soil+j) &
                = lowerVX_c_vr_acc(j+1,i_soil3  ,i) / decomp0_cpools_vr(j,i_soil3  ,i)


            ! upper triadiagnonal entries for N
            AK_soil_nacc((i_met_lit-1)*nl_soil+j,(i_met_lit-1)*nl_soil+j+1) &
                = upperVX_n_vr_acc(j,i_met_lit,i) / decomp0_cpools_vr(j+1,i_met_lit,i)
            AK_soil_nacc((i_cel_lit-1)*nl_soil+j,(i_cel_lit-1)*nl_soil+j+1) &
                = upperVX_n_vr_acc(j,i_cel_lit,i) / decomp0_cpools_vr(j+1,i_cel_lit,i)
            AK_soil_nacc((i_lig_lit-1)*nl_soil+j,(i_lig_lit-1)*nl_soil+j+1) &
                = upperVX_n_vr_acc(j,i_lig_lit,i) / decomp0_cpools_vr(j+1,i_lig_lit,i)
            AK_soil_nacc((i_soil1  -1)*nl_soil+j,(i_soil1  -1)*nl_soil+j+1) &
                = upperVX_n_vr_acc(j,i_soil1  ,i) / decomp0_cpools_vr(j+1,i_soil1  ,i)
            AK_soil_nacc((i_soil2  -1)*nl_soil+j,(i_soil2  -1)*nl_soil+j+1) &
                = upperVX_n_vr_acc(j,i_soil2  ,i) / decomp0_cpools_vr(j+1,i_soil2  ,i)
            AK_soil_nacc((i_soil3  -1)*nl_soil+j,(i_soil3  -1)*nl_soil+j+1) &
                = upperVX_n_vr_acc(j,i_soil3  ,i) / decomp0_cpools_vr(j+1,i_soil3  ,i)

            ! lower triadiagnonal entries for N
            AK_soil_nacc((i_met_lit-1)*nl_soil+j+1,(i_met_lit-1)*nl_soil+j) &
                = lowerVX_n_vr_acc(j+1,i_met_lit,i) / decomp0_npools_vr(j,i_met_lit,i)
            AK_soil_nacc((i_cel_lit-1)*nl_soil+j+1,(i_cel_lit-1)*nl_soil+j) &
                = lowerVX_n_vr_acc(j+1,i_cel_lit,i) / decomp0_npools_vr(j,i_cel_lit,i)
            AK_soil_nacc((i_lig_lit-1)*nl_soil+j+1,(i_lig_lit-1)*nl_soil+j) &
                = lowerVX_n_vr_acc(j+1,i_lig_lit,i) / decomp0_npools_vr(j,i_lig_lit,i)
            AK_soil_nacc((i_soil1  -1)*nl_soil+j+1,(i_soil1  -1)*nl_soil+j) &
                = lowerVX_n_vr_acc(j+1,i_soil1  ,i) / decomp0_npools_vr(j,i_soil1  ,i)
            AK_soil_nacc((i_soil2  -1)*nl_soil+j+1,(i_soil2  -1)*nl_soil+j) &
                = lowerVX_n_vr_acc(j+1,i_soil2  ,i) / decomp0_npools_vr(j,i_soil2  ,i)
            AK_soil_nacc((i_soil3  -1)*nl_soil+j+1,(i_soil3  -1)*nl_soil+j) &
                = lowerVX_n_vr_acc(j+1,i_soil3  ,i) / decomp0_npools_vr(j,i_soil3  ,i)

            ! N input
            I_soil_nacc((i_met_lit-1)*nl_soil+j) = I_met_n_vr_acc(j,i)
            I_soil_nacc((i_cel_lit-1)*nl_soil+j) = I_cel_n_vr_acc(j,i)
            I_soil_nacc((i_lig_lit-1)*nl_soil+j) = I_lig_n_vr_acc(j,i)
            I_soil_nacc((i_cwd    -1)*nl_soil+j) = I_cwd_n_vr_acc(j,i)

         end do

         do k=1,ndecomp_pools_vr
            if (abs(AK_soil_acc(k,k)) .le. epsi)then !avoid inversion nan
               AK_soil_acc(k,k) = - 1.e+36_r8
            end if 
         end do

         do k=1,ndecomp_pools_vr
            if (abs(AK_soil_nacc(k,k)) .le. epsi)then
               AK_soil_nacc(k,k) = - 1.e+36_r8
            end if 
         end do

         ! Calculate capacity 
         call inverse(AK_soil_acc (1:ndecomp_pools_vr,1:ndecomp_pools_vr),AKinv_soil (1:ndecomp_pools_vr,1:ndecomp_pools_vr),ndecomp_pools_vr)
         call inverse(AK_soil_nacc(1:ndecomp_pools_vr,1:ndecomp_pools_vr),AKinvn_soil(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ndecomp_pools_vr)
         soilmatrixc_cap(:,1) = -matmul(AKinv_soil(1:ndecomp_pools_vr,1:ndecomp_pools_vr), I_soil_acc (1:ndecomp_pools_vr))
         soilmatrixn_cap(:,1) = -matmul(AKinvn_soil(1:ndecomp_pools_vr,1:ndecomp_pools_vr),I_soil_nacc(1:ndecomp_pools_vr))

         do k = 1, ndecomp_pools
            do j = 1, nl_soil   
               if(soilmatrixc_cap(j+(k-1)*nl_soil,1) .lt. 0)then
                  soilmatrixc_cap(j+(k-1)*nl_soil,1) = 0._r8
               endif
               if(soilmatrixn_cap(j+(k-1)*nl_soil,1) .lt. 0)then
                  soilmatrixn_cap(j+(k-1)*nl_soil,1) = 0._r8
               endif
            end do
         end do
                     
         do k = 1, ndecomp_pools
            do j = 1, nl_soil
               if((soilmatrixc_cap(j+(k-1)*nl_soil,1)/decomp0_cpools_vr(j,k,i) .gt. 100 .and. soilmatrixc_cap(j+(k-1)*nl_soil,1) .gt. 1.e+5_r8  &
              .or. soilmatrixn_cap(j+(k-1)*nl_soil,1)/decomp0_npools_vr(j,k,i) .gt. 100 .and. soilmatrixn_cap(j+(k-1)*nl_soil,1) .gt. 1.e+3_r8) &
              .or. i .eq. i_cwd .and. (soilmatrixc_cap(j+(k-1)*nl_soil,1)/decomp0_cpools_vr(j,k,i) .gt. 100 .and. soilmatrixc_cap(j+(k-1)*nl_soil,1) .gt. 1.e+5_r8  &
              .or. soilmatrixn_cap(j+(k-1)*nl_soil,1)/decomp0_npools_vr(j,k,i) .gt. 100 .and. soilmatrixn_cap(j+(k-1)*nl_soil,1) .gt. 1.e+3_r8) )then
                  soilmatrixc_cap(j+(k-1)*nl_soil,1) = decomp_cpools_vr(j,k,i)
                  soilmatrixn_cap(j+(k-1)*nl_soil,1) = decomp_npools_vr(j,k,i)
               end if
            end do
         end do

         if(any(soilmatrixc_cap(:,1) .gt. 1.e+8_r8) .or. any(soilmatrixn_cap(:,1) .gt. 1.e+8_r8))then
            do k = 1, ndecomp_pools
               do j = 1, nl_soil
                  soilmatrixc_cap(j+(k-1)*nl_soil,1) = decomp_cpools_vr(j,k,i)
                  soilmatrixn_cap(j+(k-1)*nl_soil,1) = decomp_npools_vr(j,k,i)
               end do
            end do
         end if

         ! If spin up is on, the capacity replaces the pool size with capacity.
         ! Copy the capacity into a 3D variable, and be ready to write to history files.

         do k = 1, ndecomp_pools 
            do j = 1, nl_soil
               decomp_cpools_vr(j,k,i) = soilmatrixc_cap(j+(k-1)*nl_soil,1)
               if(floating_cn_ratio(k))then
                  decomp_npools_vr(j,k,i) = soilmatrixn_cap(j+(k-1)*nl_soil,1)
               else
                  decomp_npools_vr(j,k,i) = decomp_cpools_vr(j,k,i) / cn_decomp_pools(j,k,i)
               end if
            end do
         end do
         
         skip_balance_check(i) = .true.

         ! Reset to accumulation variables to 0 at end of each year
         do m=ps, pe
            I_leafc_p_acc        (m) = 0._r8
            I_leafc_st_p_acc     (m) = 0._r8
            I_frootc_p_acc       (m) = 0._r8
            I_frootc_st_p_acc    (m) = 0._r8
            I_livestemc_p_acc    (m) = 0._r8
            I_livestemc_st_p_acc (m) = 0._r8
            I_deadstemc_p_acc    (m) = 0._r8
            I_deadstemc_st_p_acc (m) = 0._r8
            I_livecrootc_p_acc   (m) = 0._r8
            I_livecrootc_st_p_acc(m) = 0._r8
            I_deadcrootc_p_acc   (m) = 0._r8
            I_deadcrootc_st_p_acc(m) = 0._r8
            I_grainc_p_acc       (m) = 0._r8
            I_grainc_st_p_acc    (m) = 0._r8
            I_leafn_p_acc        (m) = 0._r8
            I_leafn_st_p_acc     (m) = 0._r8
            I_frootn_p_acc       (m) = 0._r8
            I_frootn_st_p_acc    (m) = 0._r8
            I_livestemn_p_acc    (m) = 0._r8
            I_livestemn_st_p_acc (m) = 0._r8
            I_deadstemn_p_acc    (m) = 0._r8
            I_deadstemn_st_p_acc (m) = 0._r8
            I_livecrootn_p_acc   (m) = 0._r8
            I_livecrootn_st_p_acc(m) = 0._r8
            I_deadcrootn_p_acc   (m) = 0._r8
            I_deadcrootn_st_p_acc(m) = 0._r8
            I_grainn_p_acc       (m) = 0._r8
            I_grainn_st_p_acc    (m) = 0._r8

            AKX_leafc_xf_to_leafc_p_acc             (m) = 0._r8
            AKX_frootc_xf_to_frootc_p_acc           (m) = 0._r8
            AKX_livestemc_xf_to_livestemc_p_acc     (m) = 0._r8
            AKX_deadstemc_xf_to_deadstemc_p_acc     (m) = 0._r8
            AKX_livecrootc_xf_to_livecrootc_p_acc   (m) = 0._r8
            AKX_deadcrootc_xf_to_deadcrootc_p_acc   (m) = 0._r8
            AKX_grainc_xf_to_grainc_p_acc           (m) = 0._r8
            AKX_livestemc_to_deadstemc_p_acc        (m) = 0._r8
            AKX_livecrootc_to_deadcrootc_p_acc      (m) = 0._r8
           
            AKX_leafc_st_to_leafc_xf_p_acc          (m) = 0._r8
            AKX_frootc_st_to_frootc_xf_p_acc        (m) = 0._r8
            AKX_livestemc_st_to_livestemc_xf_p_acc  (m) = 0._r8
            AKX_deadstemc_st_to_deadstemc_xf_p_acc  (m) = 0._r8
            AKX_livecrootc_st_to_livecrootc_xf_p_acc(m) = 0._r8
            AKX_deadcrootc_st_to_deadcrootc_xf_p_acc(m) = 0._r8
            AKX_grainc_st_to_grainc_xf_p_acc        (m) = 0._r8

            AKX_leafc_exit_p_acc                    (m) = 0._r8
            AKX_frootc_exit_p_acc                   (m) = 0._r8
            AKX_livestemc_exit_p_acc                (m) = 0._r8
            AKX_deadstemc_exit_p_acc                (m) = 0._r8
            AKX_livecrootc_exit_p_acc               (m) = 0._r8
            AKX_deadcrootc_exit_p_acc               (m) = 0._r8
            AKX_grainc_exit_p_acc                   (m) = 0._r8

            AKX_leafc_st_exit_p_acc                 (m) = 0._r8
            AKX_frootc_st_exit_p_acc                (m) = 0._r8
            AKX_livestemc_st_exit_p_acc             (m) = 0._r8
            AKX_deadstemc_st_exit_p_acc             (m) = 0._r8
            AKX_livecrootc_st_exit_p_acc            (m) = 0._r8
            AKX_deadcrootc_st_exit_p_acc            (m) = 0._r8
            AKX_grainc_st_exit_p_acc                (m) = 0._r8

            AKX_leafc_xf_exit_p_acc                 (m) = 0._r8
            AKX_frootc_xf_exit_p_acc                (m) = 0._r8
            AKX_livestemc_xf_exit_p_acc             (m) = 0._r8
            AKX_deadstemc_xf_exit_p_acc             (m) = 0._r8
            AKX_livecrootc_xf_exit_p_acc            (m) = 0._r8
            AKX_deadcrootc_xf_exit_p_acc            (m) = 0._r8
            AKX_grainc_xf_exit_p_acc                (m) = 0._r8
           
            AKX_leafn_xf_to_leafn_p_acc             (m) = 0._r8        
            AKX_frootn_xf_to_frootn_p_acc           (m) = 0._r8
            AKX_livestemn_xf_to_livestemn_p_acc     (m) = 0._r8
            AKX_deadstemn_xf_to_deadstemn_p_acc     (m) = 0._r8
            AKX_livecrootn_xf_to_livecrootn_p_acc   (m) = 0._r8
            AKX_deadcrootn_xf_to_deadcrootn_p_acc   (m) = 0._r8
            AKX_grainn_xf_to_grainn_p_acc           (m) = 0._r8
            AKX_livestemn_to_deadstemn_p_acc        (m) = 0._r8
            AKX_livecrootn_to_deadcrootn_p_acc      (m) = 0._r8

            AKX_leafn_st_to_leafn_xf_p_acc          (m) = 0._r8
            AKX_frootn_st_to_frootn_xf_p_acc        (m) = 0._r8
            AKX_livestemn_st_to_livestemn_xf_p_acc  (m) = 0._r8
            AKX_deadstemn_st_to_deadstemn_xf_p_acc  (m) = 0._r8
            AKX_livecrootn_st_to_livecrootn_xf_p_acc(m) = 0._r8
            AKX_deadcrootn_st_to_deadcrootn_xf_p_acc(m) = 0._r8
            AKX_grainn_st_to_grainn_xf_p_acc        (m) = 0._r8

            AKX_leafn_to_retransn_p_acc             (m) = 0._r8
            AKX_frootn_to_retransn_p_acc            (m) = 0._r8
            AKX_livestemn_to_retransn_p_acc         (m) = 0._r8
            AKX_livecrootn_to_retransn_p_acc        (m) = 0._r8

            AKX_retransn_to_leafn_p_acc             (m) = 0._r8
            AKX_retransn_to_frootn_p_acc            (m) = 0._r8
            AKX_retransn_to_livestemn_p_acc         (m) = 0._r8
            AKX_retransn_to_deadstemn_p_acc         (m) = 0._r8
            AKX_retransn_to_livecrootn_p_acc        (m) = 0._r8
            AKX_retransn_to_deadcrootn_p_acc        (m) = 0._r8
            AKX_retransn_to_grainn_p_acc            (m) = 0._r8

            AKX_retransn_to_leafn_st_p_acc          (m) = 0._r8
            AKX_retransn_to_frootn_st_p_acc         (m) = 0._r8
            AKX_retransn_to_livestemn_st_p_acc      (m) = 0._r8
            AKX_retransn_to_deadstemn_st_p_acc      (m) = 0._r8
            AKX_retransn_to_livecrootn_st_p_acc     (m) = 0._r8
            AKX_retransn_to_deadcrootn_st_p_acc     (m) = 0._r8
            AKX_retransn_to_grainn_st_p_acc         (m) = 0._r8

            AKX_leafn_exit_p_acc                    (m) = 0._r8
            AKX_frootn_exit_p_acc                   (m) = 0._r8
            AKX_livestemn_exit_p_acc                (m) = 0._r8
            AKX_deadstemn_exit_p_acc                (m) = 0._r8
            AKX_livecrootn_exit_p_acc               (m) = 0._r8
            AKX_deadcrootn_exit_p_acc               (m) = 0._r8
            AKX_grainn_exit_p_acc                   (m) = 0._r8
            AKX_retransn_exit_p_acc                 (m) = 0._r8

            AKX_leafn_st_exit_p_acc                 (m) = 0._r8
            AKX_frootn_st_exit_p_acc                (m) = 0._r8
            AKX_livestemn_st_exit_p_acc             (m) = 0._r8
            AKX_deadstemn_st_exit_p_acc             (m) = 0._r8
            AKX_livecrootn_st_exit_p_acc            (m) = 0._r8
            AKX_deadcrootn_st_exit_p_acc            (m) = 0._r8
            AKX_grainn_st_exit_p_acc                (m) = 0._r8

            AKX_leafn_xf_exit_p_acc                 (m) = 0._r8
            AKX_frootn_xf_exit_p_acc                (m) = 0._r8
            AKX_livestemn_xf_exit_p_acc             (m) = 0._r8
            AKX_deadstemn_xf_exit_p_acc             (m) = 0._r8
            AKX_livecrootn_xf_exit_p_acc            (m) = 0._r8
            AKX_deadcrootn_xf_exit_p_acc            (m) = 0._r8
            AKX_grainn_xf_exit_p_acc                (m) = 0._r8
         end do

         do j=1,nl_soil
            AKX_met_exit_c_vr_acc      (j,i) = 0._r8   
            AKX_cel_exit_c_vr_acc      (j,i) = 0._r8   
            AKX_lig_exit_c_vr_acc      (j,i) = 0._r8   
            AKX_cwd_exit_c_vr_acc      (j,i) = 0._r8   
            AKX_soil1_exit_c_vr_acc    (j,i) = 0._r8   
            AKX_soil2_exit_c_vr_acc    (j,i) = 0._r8   
            AKX_soil3_exit_c_vr_acc    (j,i) = 0._r8   

            AKX_met_to_soil1_c_vr_acc  (j,i) = 0._r8
            AKX_cel_to_soil1_c_vr_acc  (j,i) = 0._r8
            AKX_lig_to_soil2_c_vr_acc  (j,i) = 0._r8
            AKX_soil1_to_soil2_c_vr_acc(j,i) = 0._r8
            AKX_cwd_to_cel_c_vr_acc    (j,i) = 0._r8
            AKX_cwd_to_lig_c_vr_acc    (j,i) = 0._r8
            AKX_soil1_to_soil3_c_vr_acc(j,i) = 0._r8
            AKX_soil2_to_soil1_c_vr_acc(j,i) = 0._r8
            AKX_soil2_to_soil3_c_vr_acc(j,i) = 0._r8
            AKX_soil3_to_soil1_c_vr_acc(j,i) = 0._r8

            diagVX_c_vr_acc  (j,i_met_lit,i) = 0._r8   
            diagVX_c_vr_acc  (j,i_cel_lit,i) = 0._r8   
            diagVX_c_vr_acc  (j,i_lig_lit,i) = 0._r8   
            diagVX_c_vr_acc  (j,i_cwd    ,i) = 0._r8   
            diagVX_c_vr_acc  (j,i_soil1  ,i) = 0._r8   
            diagVX_c_vr_acc  (j,i_soil2  ,i) = 0._r8   
            diagVX_c_vr_acc  (j,i_soil3  ,i) = 0._r8   
            
            upperVX_c_vr_acc (j,i_met_lit,i) = 0._r8   
            upperVX_c_vr_acc (j,i_cel_lit,i) = 0._r8   
            upperVX_c_vr_acc (j,i_lig_lit,i) = 0._r8   
            upperVX_c_vr_acc (j,i_cwd    ,i) = 0._r8   
            upperVX_c_vr_acc (j,i_soil1  ,i) = 0._r8   
            upperVX_c_vr_acc (j,i_soil2  ,i) = 0._r8   
            upperVX_c_vr_acc (j,i_soil3  ,i) = 0._r8   
            
            lowerVX_c_vr_acc (j,i_met_lit,i) = 0._r8   
            lowerVX_c_vr_acc (j,i_cel_lit,i) = 0._r8   
            lowerVX_c_vr_acc (j,i_lig_lit,i) = 0._r8   
            lowerVX_c_vr_acc (j,i_cwd    ,i) = 0._r8   
            lowerVX_c_vr_acc (j,i_soil1  ,i) = 0._r8   
            lowerVX_c_vr_acc (j,i_soil2  ,i) = 0._r8   
            lowerVX_c_vr_acc (j,i_soil3  ,i) = 0._r8   
            
            AKX_met_exit_n_vr_acc      (j,i) = 0._r8   
            AKX_cel_exit_n_vr_acc      (j,i) = 0._r8   
            AKX_lig_exit_n_vr_acc      (j,i) = 0._r8   
            AKX_cwd_exit_n_vr_acc      (j,i) = 0._r8   
            AKX_soil1_exit_n_vr_acc    (j,i) = 0._r8   
            AKX_soil2_exit_n_vr_acc    (j,i) = 0._r8   
            AKX_soil3_exit_n_vr_acc    (j,i) = 0._r8   

            AKX_met_to_soil1_n_vr_acc  (j,i) = 0._r8
            AKX_cel_to_soil1_n_vr_acc  (j,i) = 0._r8
            AKX_lig_to_soil2_n_vr_acc  (j,i) = 0._r8
            AKX_soil1_to_soil2_n_vr_acc(j,i) = 0._r8
            AKX_cwd_to_cel_n_vr_acc    (j,i) = 0._r8
            AKX_cwd_to_lig_n_vr_acc    (j,i) = 0._r8
            AKX_soil1_to_soil3_n_vr_acc(j,i) = 0._r8
            AKX_soil2_to_soil1_n_vr_acc(j,i) = 0._r8
            AKX_soil2_to_soil3_n_vr_acc(j,i) = 0._r8
            AKX_soil3_to_soil1_n_vr_acc(j,i) = 0._r8

            diagVX_n_vr_acc  (j,i_met_lit,i) = 0._r8   
            diagVX_n_vr_acc  (j,i_cel_lit,i) = 0._r8   
            diagVX_n_vr_acc  (j,i_lig_lit,i) = 0._r8   
            diagVX_n_vr_acc  (j,i_cwd    ,i) = 0._r8   
            diagVX_n_vr_acc  (j,i_soil1  ,i) = 0._r8   
            diagVX_n_vr_acc  (j,i_soil2  ,i) = 0._r8   
            diagVX_n_vr_acc  (j,i_soil3  ,i) = 0._r8   
            
            upperVX_n_vr_acc (j,i_met_lit,i) = 0._r8   
            upperVX_n_vr_acc (j,i_cel_lit,i) = 0._r8   
            upperVX_n_vr_acc (j,i_lig_lit,i) = 0._r8   
            upperVX_n_vr_acc (j,i_cwd    ,i) = 0._r8   
            upperVX_n_vr_acc (j,i_soil1  ,i) = 0._r8   
            upperVX_n_vr_acc (j,i_soil2  ,i) = 0._r8   
            upperVX_n_vr_acc (j,i_soil3  ,i) = 0._r8   
            
            lowerVX_n_vr_acc (j,i_met_lit,i) = 0._r8   
            lowerVX_n_vr_acc (j,i_cel_lit,i) = 0._r8   
            lowerVX_n_vr_acc (j,i_lig_lit,i) = 0._r8   
            lowerVX_n_vr_acc (j,i_cwd    ,i) = 0._r8   
            lowerVX_n_vr_acc (j,i_soil1  ,i) = 0._r8   
            lowerVX_n_vr_acc (j,i_soil2  ,i) = 0._r8   
            lowerVX_n_vr_acc (j,i_soil3  ,i) = 0._r8   
            
            I_met_c_vr_acc(j,i)              = 0._r8
            I_cel_c_vr_acc(j,i)              = 0._r8
            I_lig_c_vr_acc(j,i)              = 0._r8
            I_cwd_c_vr_acc(j,i)              = 0._r8
            
            I_met_n_vr_acc(j,i)              = 0._r8
            I_cel_n_vr_acc(j,i)              = 0._r8
            I_lig_n_vr_acc(j,i)              = 0._r8
            I_cwd_n_vr_acc(j,i)              = 0._r8
         end do
      end if

   end subroutine CNSASU
 
   subroutine inverse(a,c,n)
!============================================================
 ! Inverse matrix
 ! Method: Based on Doolittle LU factorization for Ax=b
 ! Alex G. December 2009
 !-----------------------------------------------------------
 ! input ...
 ! a(n,n) - array of coefficients for matrix A
 ! n      - dimension
 ! output ...
 ! c(n,n) - inverse matrix of A
 ! comments ...
 ! the original matrix a(n,n) will be destroyed
 ! during the calculation
 !===========================================================
   implicit none
   ! Arguments
   integer,intent(in) :: n           ! Size of matrix
   real(r8),intent(in)  :: a(:,:)    ! Input matrix to fine the inverse of
   real(r8),intent(out) :: c(:,:)    ! Output inverse
   ! Local variables
   real(r8) :: L(n,n)   ! matrix of the elimination coefficient
   real(r8) :: U(n,n)   ! Upper triangular part of input matrix A
   real(r8) :: aa(n,n)  ! Temporary equal to input matrix a
   real(r8) :: b(n)     ! Temporary vector
   real(r8) :: d(n)     ! Temporary vector (solution of L*d)
   real(r8) :: x(n)     ! Temporary vector (U*x = d)
   real(r8) :: coeff    ! coefficient
   integer i, j, k      ! Indices
   character(len=*), parameter :: subname = 'inverse'

   do k=1,n
      if ( a(k,k) == 0.0_r8 )then
         call abort
      end if
   end do
   !
   ! step 0: initialization for matrices L and U and b
   ! Fortran 90/95 aloows such operations on matrices
   !
   L=0.0
   U=0.0
   b=0.0

   aa=a
   !
   ! Step 1: forward elimination
   !
   do k=1, n-1
      do i=k+1,n
         ! Already verifieid that divisor isn't zero
         coeff=aa(i,k)/aa(k,k)
         L(i,k) = coeff
         do j=k+1,n
            aa(i,j) = aa(i,j)-coeff*aa(k,j)
         end do
      end do
   end do

   !
   ! Step 2: prepare L and U matrices
   ! L matrix is a matrix of the elimination coefficient
   ! + the diagonal elements are 1.0
   !
   do i=1,n
     L(i,i) = 1.0
   end do
   !
   ! U matrix is the upper triangular part of A
   !
   do j=1,n
     do i=1,j
       U(i,j) = aa(i,j)
     end do
   end do
   !
   ! Step 3: compute columns of the inverse matrix C
   !
   do k=1,n
     b(k)=1.0
     d(1) = b(1)
     ! Step 3a: Solve Ld=b using the forward substitution
     do i=2,n
       d(i)=b(i)
       do j=1,i-1
         d(i) = d(i) - L(i,j)*d(j)
       end do
     end do
     ! Step 3b: Solve Ux=d using the back substitution
     x(n)=d(n)/U(n,n)
     do i = n-1,1,-1
       x(i) = d(i)
       do j=n,i+1,-1
         x(i)=x(i)-U(i,j)*x(j)
       end do
       ! Already verifieid that divisor isn't zero
       x(i) = x(i)/u(i,i)
     end do
     ! Step 3c: fill the solutions x(n) into column k of C
     do i=1,n
       c(i,k) = x(i)
     end do
     b(k)=0.0
   end do

 end subroutine inverse

end module bgc_CNSASUMod

#endif
