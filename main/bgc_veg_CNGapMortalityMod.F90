module bgc_veg_CNGapMortalityMod

use precision
use PFT_Const, only: lf_flab, lf_fcel, lf_flig, fr_flab, fr_fcel, fr_flig
use MOD_TimeInvariants, only: &
  ! bgc constants
    am
use MOD_PFTimeInvars, only: pftclass, pftfrac

use MOD_1D_Fluxes, only: &
    ! decomposition carbon flux varables (in)
           gap_mortality_to_met_c, gap_mortality_to_cel_c , &
           gap_mortality_to_lig_c, gap_mortality_to_cwdc  , &
  
    ! decompositionn nitrogen fluxes variables (inout)
           gap_mortality_to_met_n, gap_mortality_to_cel_n , &
           gap_mortality_to_lig_n, gap_mortality_to_cwdn 

use MOD_1D_PFTFluxes, only: &
    ! vegetation carbon flux variables
           m_leafc_to_litter_p        , m_leafc_storage_to_litter_p     , m_leafc_xfer_to_litter_p     , &
           m_frootc_to_litter_p       , m_frootc_storage_to_litter_p    , m_frootc_xfer_to_litter_p    , &
           m_livestemc_to_litter_p    , m_livestemc_storage_to_litter_p , m_livestemc_xfer_to_litter_p , &
           m_deadstemc_to_litter_p    , m_deadstemc_storage_to_litter_p , m_deadstemc_xfer_to_litter_p , &
           m_livecrootc_to_litter_p   , m_livecrootc_storage_to_litter_p, m_livecrootc_xfer_to_litter_p, &
           m_deadcrootc_to_litter_p   , m_deadcrootc_storage_to_litter_p, m_deadcrootc_xfer_to_litter_p, &
           m_gresp_storage_to_litter_p, m_gresp_xfer_to_litter_p        , &

    ! vegetation nitrogen flux variables
           m_leafn_to_litter_p        , m_leafn_storage_to_litter_p     , m_leafn_xfer_to_litter_p     , &
           m_frootn_to_litter_p       , m_frootn_storage_to_litter_p    , m_frootn_xfer_to_litter_p    , &
           m_livestemn_to_litter_p    , m_livestemn_storage_to_litter_p , m_livestemn_xfer_to_litter_p , &
           m_deadstemn_to_litter_p    , m_deadstemn_storage_to_litter_p , m_deadstemn_xfer_to_litter_p , &
           m_livecrootn_to_litter_p   , m_livecrootn_storage_to_litter_p, m_livecrootn_xfer_to_litter_p, &
           m_deadcrootn_to_litter_p   , m_deadcrootn_storage_to_litter_p, m_deadcrootn_xfer_to_litter_p, &
           m_retransn_to_litter_p

use MOD_PFTimeVars, only: &
    ! vegetation carbon state variables (inout)
           leafc_p            , leafc_storage_p     , leafc_xfer_p     , &
           frootc_p           , frootc_storage_p    , frootc_xfer_p    , &
           livestemc_p        , livestemc_storage_p , livestemc_xfer_p , &
           deadstemc_p        , deadstemc_storage_p , deadstemc_xfer_p , &
           livecrootc_p       , livecrootc_storage_p, livecrootc_xfer_p, &
           deadcrootc_p       , deadcrootc_storage_p, deadcrootc_xfer_p, &
           gresp_storage_p    , gresp_xfer_p        , &

    ! vegetation nitrogen state variables (inout)
           leafn_p            , leafn_storage_p     , leafn_xfer_p     , &
           frootn_p           , frootn_storage_p    , frootn_xfer_p    , &
           livestemn_p        , livestemn_storage_p , livestemn_xfer_p , &
           deadstemn_p        , deadstemn_storage_p , deadstemn_xfer_p , &
           livecrootn_p       , livecrootn_storage_p, livecrootn_xfer_p, &
           deadcrootn_p       , deadcrootn_storage_p, deadcrootn_xfer_p, &
           retransn_p         , &

    ! profiles
           leaf_prof_p,       stem_prof_p,        froot_prof_p,        croot_prof_p

implicit none

public CNGapMortality

private CNGap_VegToLitter

contains

subroutine CNGapMortality(i, ps, pe, nl_soil, npcropmin)

integer ,intent(in) :: i
integer ,intent(in) :: ps
integer ,intent(in) :: pe
integer ,intent(in) :: nl_soil
integer ,intent(in) :: npcropmin

real(r8):: mort             ! rate for fractional mortality (1/s)
integer :: ivt, m
!    real(r8):: mort_max      ! asymptotic max mortality rate (/yr)

!         if (use_cndv) then
!            ! Stress mortality from lpj's subr Mortality.
!
!            if (woody(ivt(p)) == 1._r8) then
!
!               if (ivt(p) == 8) then
!                  mort_max = 0.03_r8 ! BDT boreal
!               else
!                  mort_max = 0.01_r8 ! original value for all patches
!               end if
!
!               ! heatstress and greffic calculated in Establishment once/yr
!
!               ! Mortality rate inversely related to growth efficiency
!               ! (Prentice et al 1993)
!               am = mort_max / (1._r8 + k_mort * greffic(p))
!
!               ! Mortality rate inversely related to growth efficiency
!               ! (Prentice et al 1993)
!               am = mort_max / (1._r8 + k_mort * greffic(p))
!
!               am = min(1._r8, am + heatstress(p))
!            else ! lpj didn't set this for grasses; cn does
               ! set the mortality rate based on annual rate
!               am = params_inst%am
!            end if
!
!         end if
      do m = ps , pe
         ivt = pftclass(m)

         mort  = am/(365._r8 * 86400._r8)

         !------------------------------------------------------
         ! patch-level gap mortality carbon fluxes
         !------------------------------------------------------

         ! displayed pools
         m_leafc_to_litter_p(m)               = leafc_p(m)               * mort
         m_frootc_to_litter_p(m)              = frootc_p(m)              * mort
         m_livestemc_to_litter_p(m)           = livestemc_p(m)           * mort
         m_livecrootc_to_litter_p(m)          = livecrootc_p(m)          * mort
!         if (spinup_state == 2 .and. .not. use_cndv) then   !accelerate mortality of dead woody pools
!           m_deadstemc_to_litter(i)         = deadstemc(i)  * m * 10._r8
!           m_deadcrootc_to_litter(i)        = deadcrootc(i) * m * 10._r8
!           if (use_matrixcn) then
!              matrix_gmtransfer_patch(p,ideadstem_to_iout_gmc)     = matrix_gmtransfer_patch(p,ideadstem_to_iout_gmc)     + m * 10._r8
!              matrix_gmtransfer_patch(p,ideadcroot_to_iout_gmc)    = matrix_gmtransfer_patch(p,ideadcroot_to_iout_gmc)    + m * 10._r8
!           end if !use_matrixcn
!         else
           m_deadstemc_to_litter_p(m)         = deadstemc_p(m)           * mort
           m_deadcrootc_to_litter_p(m)        = deadcrootc_p(m)          * mort
!           if (use_matrixcn) then
!              matrix_gmtransfer_patch(p,ideadstem_to_iout_gmc)     = matrix_gmtransfer_patch(p,ideadstem_to_iout_gmc)     + m
!              matrix_gmtransfer_patch(p,ideadcroot_to_iout_gmc)    = matrix_gmtransfer_patch(p,ideadcroot_to_iout_gmc)    + m
!           end if !use_matrixcn
!         end if

         ! storage pools
         m_leafc_storage_to_litter_p(m)       = leafc_storage_p(m)       * mort
         m_frootc_storage_to_litter_p(m)      = frootc_storage_p(m)      * mort
         m_livestemc_storage_to_litter_p(m)   = livestemc_storage_p(m)   * mort
         m_deadstemc_storage_to_litter_p(m)   = deadstemc_storage_p(m)   * mort
         m_livecrootc_storage_to_litter_p(m)  = livecrootc_storage_p(m)  * mort
         m_deadcrootc_storage_to_litter_p(m)  = deadcrootc_storage_p(m)  * mort
         m_gresp_storage_to_litter_p(m)       = gresp_storage_p(m)       * mort

         ! transfer pools
         m_leafc_xfer_to_litter_p(m)          = leafc_xfer_p(m)          * mort
         m_frootc_xfer_to_litter_p(m)         = frootc_xfer_p(m)         * mort
         m_livestemc_xfer_to_litter_p(m)      = livestemc_xfer_p(m)      * mort
         m_deadstemc_xfer_to_litter_p(m)      = deadstemc_xfer_p(m)      * mort
         m_livecrootc_xfer_to_litter_p(m)     = livecrootc_xfer_p(m)     * mort
         m_deadcrootc_xfer_to_litter_p(m)     = deadcrootc_xfer_p(m)     * mort
         m_gresp_xfer_to_litter_p(m)          = gresp_xfer_p(m)          * mort
!         if (use_matrixcn) then
!         ! displayed pools
!            matrix_gmtransfer_patch(p,ileaf_to_iout_gmc)         = m
!            matrix_gmtransfer_patch(p,ifroot_to_iout_gmc)        = m
!            matrix_gmtransfer_patch(p,ilivestem_to_iout_gmc)     = m
!            matrix_gmtransfer_patch(p,ilivecroot_to_iout_gmc)    = m
!
!         ! storage pools
!            matrix_gmtransfer_patch(p,ileafst_to_iout_gmc)       = m
!            matrix_gmtransfer_patch(p,ifrootst_to_iout_gmc)      = m
!            matrix_gmtransfer_patch(p,ilivestemst_to_iout_gmc)   = m
!            matrix_gmtransfer_patch(p,ilivecrootst_to_iout_gmc)  = m
!            matrix_gmtransfer_patch(p,ideadstemst_to_iout_gmc)   = m
!            matrix_gmtransfer_patch(p,ideadcrootst_to_iout_gmc)  = m
!
!         ! transfer pools
!            matrix_gmtransfer_patch(p,ileafxf_to_iout_gmc)       = m
!            matrix_gmtransfer_patch(p,ifrootxf_to_iout_gmc)      = m
!            matrix_gmtransfer_patch(p,ilivestemxf_to_iout_gmc)   = m
!            matrix_gmtransfer_patch(p,ilivecrootxf_to_iout_gmc)  = m
!            matrix_gmtransfer_patch(p,ideadstemxf_to_iout_gmc)   = m
!            matrix_gmtransfer_patch(p,ideadcrootxf_to_iout_gmc)  = m
!         end if !use_matrixcn

         !------------------------------------------------------
         ! patch-level gap mortality nitrogen fluxes
         !------------------------------------------------------

         ! displayed pools
         m_leafn_to_litter_p(m)            = leafn_p(m)               * mort
         m_frootn_to_litter_p(m)           = frootn_p(m)              * mort
         m_livestemn_to_litter_p(m)        = livestemn_p(m)           * mort
         m_livecrootn_to_litter_p(m)       = livecrootn_p(m)          * mort

!         if (spinup_state == 2 .and. .not. use_cndv) then   !accelerate mortality of dead woody pools
!           m_deadstemn_to_litter_patch(p)      = deadstemn_patch(p)  * m * 10._r8
!           m_deadcrootn_to_litter_patch(p)     = deadcrootn_patch(p) * m * 10._r8
!           if (use_matrixcn) then
!             matrix_ngmtransfer_patch(p,ideadstem_to_iout_gmn)     = matrix_ngmtransfer_patch(p,ideadstem_to_iout_gmn)     + m* 10._r8
!             matrix_ngmtransfer_patch(p,ideadcroot_to_iout_gmn)    = matrix_ngmtransfer_patch(p,ideadcroot_to_iout_gmn)    + m*10._r8
!           end if !use_matrixcn
!         else
           m_deadstemn_to_litter_p(m)      = deadstemn_p(m)           * mort
           m_deadcrootn_to_litter_p(m)     = deadcrootn_p(m)          * mort
!           if (use_matrixcn) then
!             matrix_ngmtransfer_patch(p,ideadstem_to_iout_gmn)     = matrix_ngmtransfer_patch(p,ideadstem_to_iout_gmn)     + m
!             matrix_ngmtransfer_patch(p,ideadcroot_to_iout_gmn)    = matrix_ngmtransfer_patch(p,ideadcroot_to_iout_gmn)    + m
!           end if !use_matrixcn
!        end if

         if (ivt < npcropmin) then
            m_retransn_to_litter_p(m) = retransn_p(m) * mort
         end if

         ! storage pools
         m_leafn_storage_to_litter_p(m)       = leafn_storage_p(m)      * mort
         m_frootn_storage_to_litter_p(m)      = frootn_storage_p(m)     * mort
         m_livestemn_storage_to_litter_p(m)   = livestemn_storage_p(m)  * mort
         m_deadstemn_storage_to_litter_p(m)   = deadstemn_storage_p(m)  * mort
         m_livecrootn_storage_to_litter_p(m)  = livecrootn_storage_p(m) * mort
         m_deadcrootn_storage_to_litter_p(m)  = deadcrootn_storage_p(m) * mort

         ! transfer pools
         m_leafn_xfer_to_litter_p(m)          = leafn_xfer_p(m)         * mort
         m_frootn_xfer_to_litter_p(m)         = frootn_xfer_p(m)        * mort
         m_livestemn_xfer_to_litter_p(m)      = livestemn_xfer_p(m)     * mort
         m_deadstemn_xfer_to_litter_p(m)      = deadstemn_xfer_p(m)     * mort
         m_livecrootn_xfer_to_litter_p(m)     = livecrootn_xfer_p(m)    * mort
         m_deadcrootn_xfer_to_litter_p(m)     = deadcrootn_xfer_p(m)    * mort

!         if (use_matrixcn) then
         ! displayed pools
!            matrix_ngmtransfer_patch(p,ileaf_to_iout_gmn)        = m
!            matrix_ngmtransfer_patch(p,ifroot_to_iout_gmn)       = m
!            matrix_ngmtransfer_patch(p,ilivestem_to_iout_gmn)    = m
!            matrix_ngmtransfer_patch(p,ilivecroot_to_iout_gmn)   = m

         ! storage pools
!            matrix_ngmtransfer_patch(p,ileafst_to_iout_gmn)      = m
!            matrix_ngmtransfer_patch(p,ifrootst_to_iout_gmn)     = m
!            matrix_ngmtransfer_patch(p,ilivestemst_to_iout_gmn)  = m
!            matrix_ngmtransfer_patch(p,ilivecrootst_to_iout_gmn) = m
!            matrix_ngmtransfer_patch(p,ideadstemst_to_iout_gmn)  = m
!            matrix_ngmtransfer_patch(p,ideadcrootst_to_iout_gmn) = m

         ! transfer pools
!            matrix_ngmtransfer_patch(p,ileafxf_to_iout_gmn)      = m
!            matrix_ngmtransfer_patch(p,ifrootxf_to_iout_gmn)     = m
!            matrix_ngmtransfer_patch(p,ilivestemxf_to_iout_gmn)  = m
!            matrix_ngmtransfer_patch(p,ilivecrootxf_to_iout_gmn) = m
!            matrix_ngmtransfer_patch(p,ideadstemxf_to_iout_gmn)  = m
!            matrix_ngmtransfer_patch(p,ideadcrootxf_to_iout_gmn) = m

!            if (ivt(p) < npcropmin) then
!               matrix_ngmtransfer_patch(p,iretransn_to_iout_gmn) = m
!            end if
!         end if !use_matrixcn
       end do

!      if(use_matrixcn)then
!         do fp = 1,num_soilp
!            p = filter_soilp(fp)
!            if(matrix_phtransfer_patch(p,ifroot_to_iout_ph) .ge. 1._r8 / get_step_size_real())then
!               m_frootc_to_litter_patch(p) = 0._r8
!               matrix_gmtransfer_patch(p,ifroot_to_iout_gmc) = 0._r8
!            end if
!         end do
!      end if

       call CNGap_VegToLitter(i, ps, pe, nl_soil)

end subroutine CNGapMortality

subroutine CNGap_VegToLitter(i, ps, pe, nl_soil)

integer ,intent(in) :: i
integer ,intent(in) :: ps
integer ,intent(in) :: pe
integer ,intent(in) :: nl_soil

integer j,m,ivt,wtcol

   do j = 1,nl_soil
      do m = ps, pe
         ivt = pftclass(m)
         wtcol = pftfrac(m)

         ! leaf gap mortality carbon fluxes
         gap_mortality_to_met_c(j,i) = gap_mortality_to_met_c(j,i) + &
              m_leafc_to_litter_p(m) * lf_flab(ivt) * wtcol * leaf_prof_p(j,m)
         gap_mortality_to_cel_c(j,i) = gap_mortality_to_cel_c(j,i) + &
              m_leafc_to_litter_p(m) * lf_fcel(ivt) * wtcol * leaf_prof_p(j,m)
         gap_mortality_to_lig_c(j,i) = gap_mortality_to_lig_c(j,i) + &
              m_leafc_to_litter_p(m) * lf_flig(ivt) * wtcol * leaf_prof_p(j,m)

         ! fine root gap mortality carbon fluxes
         gap_mortality_to_met_c(j,i) = gap_mortality_to_met_c(j,i) + &
              m_frootc_to_litter_p(m) * fr_flab(ivt) * wtcol * froot_prof_p(j,m)
         gap_mortality_to_cel_c(j,i) = gap_mortality_to_cel_c(j,i) + &
              m_frootc_to_litter_p(m) * fr_fcel(ivt) * wtcol * froot_prof_p(j,m)
         gap_mortality_to_lig_c(j,i) = gap_mortality_to_lig_c(j,i) + &
              m_frootc_to_litter_p(m) * fr_flig(ivt) * wtcol * froot_prof_p(j,m)

         ! wood gap mortality carbon fluxes
         gap_mortality_to_cwdc(j,i)  = gap_mortality_to_cwdc(j,i)  + &
              (m_livestemc_to_litter_p(m) + m_deadstemc_to_litter_p(m))  * wtcol * stem_prof_p(j,m)
         gap_mortality_to_cwdc(j,i) = gap_mortality_to_cwdc(j,i) + &
              (m_livecrootc_to_litter_p(m) + m_deadcrootc_to_litter_p(m)) * wtcol * croot_prof_p(j,m)

         ! storage gap mortality carbon fluxes
         gap_mortality_to_met_c(j,i) = gap_mortality_to_met_c(j,i) + &
              (m_leafc_storage_to_litter_p(m) + m_gresp_storage_to_litter_p(m)) * wtcol * leaf_prof_p(j,m)
         gap_mortality_to_met_c(j,i) = gap_mortality_to_met_c(j,i) + &
              m_frootc_storage_to_litter_p(m) * wtcol * froot_prof_p(j,m)
         gap_mortality_to_met_c(j,i) = gap_mortality_to_met_c(j,i)  + &
              (m_livestemc_storage_to_litter_p(m) + m_deadstemc_storage_to_litter_p(m)) * wtcol * stem_prof_p(j,m)
         gap_mortality_to_met_c(j,i) = gap_mortality_to_met_c(j,i) + &
              (m_livecrootc_storage_to_litter_p(m) + m_deadcrootc_storage_to_litter_p(m)) * wtcol * croot_prof_p(j,m)

         ! transfer gap mortality carbon fluxes
         gap_mortality_to_met_c(j,i) = gap_mortality_to_met_c(j,i) + &
              (m_leafc_xfer_to_litter_p(m) + m_gresp_xfer_to_litter_p(m)) * wtcol * leaf_prof_p(j,m)
         gap_mortality_to_met_c(j,i) = gap_mortality_to_met_c(j,i) + &
              m_frootc_xfer_to_litter_p(m) * wtcol * froot_prof_p(j,m)
         gap_mortality_to_met_c(j,i)  = gap_mortality_to_met_c(j,i)  + &
              (m_livestemc_xfer_to_litter_p(m) + m_deadstemc_xfer_to_litter_p(m))  * wtcol * stem_prof_p(j,m)
         gap_mortality_to_met_c(j,i) = gap_mortality_to_met_c(j,i) + &
              (m_livecrootc_xfer_to_litter_p(m) + m_deadcrootc_xfer_to_litter_p(m)) * wtcol * croot_prof_p(j,m)

         ! leaf gap mortality nitrogen fluxes
         gap_mortality_to_met_n(j,i) = gap_mortality_to_met_n(j,i) + &
              m_leafn_to_litter_p(m) * lf_flab(ivt) * wtcol * leaf_prof_p(j,m)
         gap_mortality_to_cel_n(j,i) = gap_mortality_to_cel_n(j,i) + &
              m_leafn_to_litter_p(m) * lf_fcel(ivt) * wtcol * leaf_prof_p(j,m)
         gap_mortality_to_lig_n(j,i) = gap_mortality_to_lig_n(j,i) + &
              m_leafn_to_litter_p(m) * lf_flig(ivt) * wtcol * leaf_prof_p(j,m)

         ! fine root litter nitrogen fluxes
         gap_mortality_to_met_n(j,i) = gap_mortality_to_met_n(j,i) + &
              m_frootn_to_litter_p(m) * fr_flab(ivt) * wtcol * froot_prof_p(j,m)
         gap_mortality_to_cel_n(j,i) = gap_mortality_to_cel_n(j,i) + &
              m_frootn_to_litter_p(m) * fr_fcel(ivt) * wtcol * froot_prof_p(j,m)
         gap_mortality_to_lig_n(j,i) = gap_mortality_to_lig_n(j,i) + &
              m_frootn_to_litter_p(m) * fr_flig(ivt) * wtcol * froot_prof_p(j,m)

         ! wood gap mortality nitrogen fluxes
         gap_mortality_to_cwdn(j,i) = gap_mortality_to_cwdn(j,i)  + &
              (m_livestemn_to_litter_p(m) + m_deadstemn_to_litter_p(m))  * wtcol * stem_prof_p(j,m)
         gap_mortality_to_cwdn(j,i) = gap_mortality_to_cwdn(j,i) + &
              (m_livecrootn_to_litter_p(m) + m_deadcrootn_to_litter_p(m)) * wtcol * croot_prof_p(j,m)

         ! retranslocated N pool gap mortality fluxes
         gap_mortality_to_met_n(j,i) = gap_mortality_to_met_n(j,i) + &
              m_retransn_to_litter_p(m) * wtcol * leaf_prof_p(j,m)

         ! storage gap mortality nitrogen fluxes
         gap_mortality_to_met_n(j,i) = gap_mortality_to_met_n(j,i) + &
              m_leafn_storage_to_litter_p(m) * wtcol * leaf_prof_p(j,m)
         gap_mortality_to_met_n(j,i) = gap_mortality_to_met_n(j,i) + &
              m_frootn_storage_to_litter_p(m) * wtcol * froot_prof_p(j,m)
         gap_mortality_to_met_n(j,i)  = gap_mortality_to_met_n(j,i) + &
              (m_livestemn_storage_to_litter_p(m) + m_deadstemn_storage_to_litter_p(m))  * wtcol * stem_prof_p(j,m)
         gap_mortality_to_met_n(j,i) = gap_mortality_to_met_n(j,i) + &
              (m_livecrootn_storage_to_litter_p(m) + m_deadcrootn_storage_to_litter_p(m)) * wtcol * croot_prof_p(j,m)

         ! transfer gap mortality nitrogen fluxes
         gap_mortality_to_met_n(j,i) = gap_mortality_to_met_n(j,i) + &
              m_leafn_xfer_to_litter_p(m) * wtcol * leaf_prof_p(j,m)
         gap_mortality_to_met_n(j,i) = gap_mortality_to_met_n(j,i) + &
              m_frootn_xfer_to_litter_p(m) * wtcol * froot_prof_p(j,m)
         gap_mortality_to_met_n(j,i) = gap_mortality_to_met_n(j,i) + &
              (m_livestemn_xfer_to_litter_p(m) + m_deadstemn_xfer_to_litter_p(m))  * wtcol * stem_prof_p(j,m)
         gap_mortality_to_met_n(j,i) = gap_mortality_to_met_n(j,i) + &
              (m_livecrootn_xfer_to_litter_p(m) + m_deadcrootn_xfer_to_litter_p(m)) * wtcol * croot_prof_p(j,m)


      end do
   end do

end subroutine CNGap_VegToLitter

end module bgc_veg_CNGapMortalityMod
