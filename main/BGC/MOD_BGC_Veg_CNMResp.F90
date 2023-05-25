#include <define.h>
#ifdef BGC
module MOD_BGC_Veg_CNMResp

  !-----------------------------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module calculates plant maintenance respiration 
  !
  ! !REFERENCE:
  ! Atkin OK, Bloomfield KJ, Reich PB, Tjoelker MG, Asner GP, Bonal D et al (2015) Global variability in leaf respiration
  ! in relation to climate, plant functional types and leaf traits. New Phytologist 206:614–636
  !
  ! !ORIGINAL:
  ! The Community Land Model version 5.0 (CLM5)
  !
  ! !REVISION:
  ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.

  use precision
  use MOD_BGC_Vars_TimeInvars, only: &
      Q10,br, br_root
  use MOD_Vars_PFTimeInvars, only: pftclass
  use MOD_Vars_TimeVariables, only: &
      t_soisno, tref
  use MOD_Vars_PFTimeVars, only: &
      laisun_p, laisha_p, sigf_p
  use MOD_BGC_Vars_PFTimeVars, only: &
      frootn_p, livestemn_p, livecrootn_p, grainn_p
  use MOD_Vars_1DPFTFluxes, only: &
      respc_p
  use MOD_BGC_Vars_1DPFTFluxes, only: &
      leaf_mr_p, froot_mr_p, livestem_mr_p, livecroot_mr_p, grain_mr_p
  use MOD_Const_PFT, only: &
      woody, rootfr_p

  implicit none

  public CNMResp

contains

  subroutine CNMResp(i, ps, pe, nl_soil, npcropmin)

    integer ,intent(in) :: i         ! patch index
    integer ,intent(in) :: ps        ! start pft index
    integer ,intent(in) :: pe        ! end pft index
    integer ,intent(in) :: nl_soil   ! number of total soil layers
    integer ,intent(in) :: npcropmin ! first crop pft index

    ! !LOCAL VARIABLES:
    integer :: j   ! indices
    integer :: ivt, m

    real(r8):: tc      ! temperature correction, 2m air temp (unitless)
    real(r8):: tcsoi(nl_soil) ! temperature correction by soil layer (unitless)

    ! base rate for maintenance respiration is from:
    ! M. Ryan, 1991. Effects of climate change on plant respiration.
    ! Ecological Applications, 1(2), 157-167.
    ! Original expression is br = 0.0106 molC/(molN h)
    ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
    ! set constants

    ! Peter Thornton: 3/13/09
    ! Q10 was originally set to 2.0, an arbitrary choice, but reduced to 1.5 as part of the tuning
    ! to improve seasonal cycle of atmospheric CO2 concentration in global
    ! simulatoins

    ! column loop to calculate temperature factors in each soil layer
    do j=1,nl_soil

       ! calculate temperature corrections for each soil layer, for use in
       ! estimating fine root maintenance respiration with depth
       tcsoi(j) = Q10**((t_soisno(j,i) - 273.15_r8 - 20.0_r8)/10.0_r8)
    end do

    ! calculate maintenance respiration fluxes in
    ! gC/m2/s for each of the live plant tissues.
    ! Leaf and live wood MR

    tc = Q10**((tref(i) - 273.15_r8 - 20.0_r8)/10.0_r8)

    !RF: acclimation of root and stem respiration fluxes
    ! n.b. we do not yet know if this is defensible scientifically (awaiting data analysis)
    ! turning this on will increase R and decrease productivity in boreal forests, A LOT. :)

    do m = ps, pe
       ivt = pftclass(m)
       if (sigf_p(m) == 1) then
          leaf_mr_p(m) = respc_p(m) * 12.011_r8
       else !nosno
          leaf_mr_p(m) = 0._r8
       end if

       if (woody(ivt) == 1) then
          livestem_mr_p (m) = livestemn_p (m)*br*tc
          livecroot_mr_p(m) = livecrootn_p(m)*br_root*tc
       else if (ivt >= npcropmin) then
          livestem_mr_p (m) = livestemn_p (m)*br*tc
          grain_mr_p    (m) = grainn_p    (m)*br*tc
       end if
       ! soil and patch loop for fine root

       do j = 1,nl_soil

       ! Fine root MR
       ! crootfr(j) sums to 1.0 over all soil layers, and
       ! describes the fraction of root mass for carbon that is in each
       ! layer.  This is used with the layer temperature correction
       ! to estimate the total fine root maintenance respiration as a
       ! function of temperature and N content.
          froot_mr_p(m) = froot_mr_p(m) + frootn_p(m)*br_root*tcsoi(j)*rootfr_p(j,ivt)
       end do
    end do

  end subroutine CNMResp

end module MOD_BGC_Veg_CNMResp
#endif
