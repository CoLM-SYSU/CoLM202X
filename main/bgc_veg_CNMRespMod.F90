module bgc_veg_CNMRespMod

use precision
use MOD_TimeInvariants, only: &
    Q10,br, br_root
use MOD_PFTimeInvars, only: pftclass
use MOD_TimeVariables, only: &
    t_soisno, tref
use MOD_PFTimeVars, only: &
    laisun_p, laisha_p, frootn_p, livestemn_p, livecrootn_p, grainn_p, sigf_p
use MOD_1D_PFTFluxes, only: &
    leaf_mr_p, froot_mr_p, livestem_mr_p, livecroot_mr_p, grain_mr_p, respc_p
use PFT_Const, only: &
    woody, rootfr_p

implicit none

public CNMResp

contains

  subroutine CNMResp(i, ps, pe, nl_soil, npcropmin)

    integer ,intent(in) :: i
    integer ,intent(in) :: ps
    integer ,intent(in) :: pe
    integer ,intent(in) :: nl_soil
    integer ,intent(in) :: npcropmin

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
!      br      = params_inst%br
!      br_root = params_inst%br_root

      ! Peter Thornton: 3/13/09
      ! Q10 was originally set to 2.0, an arbitrary choice, but reduced to 1.5 as part of the tuning
      ! to improve seasonal cycle of atmospheric CO2 concentration in global
      ! simulatoins
!      Q10 = CNParamsShareInst%Q10
      ! column loop to calculate temperature factors in each soil layer
      do j=1,nl_soil

         ! calculate temperature corrections for each soil layer, for use in
         ! estimating fine root maintenance respiration with depth
         tcsoi(j) = Q10**((t_soisno(j,i) - 273.15_r8 - 20.0_r8)/10.0_r8)
      end do

      ! patch loop for leaves and live wood

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

end module bgc_veg_CNMRespMod
