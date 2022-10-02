#include <define.h>
#ifdef BGC
module bgc_CNSummaryMod

use precision
use MOD_PFTimeInvars, only: pftclass
use MOD_BGCTimeVars, only: &
    totlitc, totsomc, totcwdc, decomp_cpools, decomp_cpools_vr, ctrunc_soil,ctrunc_veg, ctrunc_vr, &
    totlitn, totsomn, totcwdn, decomp_npools, decomp_npools_vr, ntrunc_soil,ntrunc_veg, ntrunc_vr, &
    totvegc, totvegn, totcolc, totcoln, sminn, sminn_vr, &
    leafc, frootc, livestemc, deadstemc, livecrootc, deadcrootc, leafc_storage, frootc_storage, livestemc_storage, &
    deadstemc_storage, livecrootc_storage, deadcrootc_storage, leafc_xfer, frootc_xfer, livestemc_xfer, &
    deadstemc_xfer, livecrootc_xfer, deadcrootc_xfer, xsmrpool, grainc, grainc_storage, grainc_xfer, &
    cropseedc_deficit, cropprod1c, cphase, & 
    fertnitro_corn, fertnitro_swheat, fertnitro_wwheat, fertnitro_soybean, &
    fertnitro_cotton, fertnitro_rice1, fertnitro_rice2, fertnitro_sugarcane, &
    leafn, frootn, livestemn, deadstemn, livecrootn, deadcrootn, leafn_storage, frootn_storage, livestemn_storage, &
    deadstemn_storage, livecrootn_storage, deadcrootn_storage, leafn_xfer, frootn_xfer, livestemn_xfer, &
    deadstemn_xfer, livecrootn_xfer, deadcrootn_xfer, retransn, grainn, grainn_storage, grainn_xfer, downreg
use MOD_BGCTimeInvars, only: &
    is_litter, is_soil, is_cwd
use MOD_BGCPFTimeVars, only: &
    leafc_p, frootc_p, livestemc_p, deadstemc_p, livecrootc_p, deadcrootc_p, &
    leafc_storage_p, frootc_storage_p, livestemc_storage_p, &
    deadstemc_storage_p, livecrootc_storage_p, deadcrootc_storage_p, gresp_storage_p, &
    leafc_xfer_p, frootc_xfer_p, livestemc_xfer_p, &
    deadstemc_xfer_p, livecrootc_xfer_p, deadcrootc_xfer_p, gresp_xfer_p, xsmrpool_p, &
    grainc_p, grainc_storage_p, grainc_xfer_p, ctrunc_p, totvegc_p, &
    cropseedc_deficit_p, cropprod1c_p, cpool_p, &
#ifdef CROP
    cphase_p, fertnitro_p, &
#endif
    leafn_p, frootn_p, livestemn_p, deadstemn_p, livecrootn_p, deadcrootn_p, &
    leafn_storage_p, frootn_storage_p, livestemn_storage_p, &
    deadstemn_storage_p, livecrootn_storage_p, deadcrootn_storage_p, &
    leafn_xfer_p, frootn_xfer_p, livestemn_xfer_p, &
    deadstemn_xfer_p, livecrootn_xfer_p, deadcrootn_xfer_p, retransn_p, npool_p, &
    grainn_p, grainn_storage_p, grainn_xfer_p, ntrunc_p, totvegn_p, cropseedn_deficit_p, downreg_p
use MOD_PFTimeInvars,  only: pftfrac
use MOD_1D_BGCFluxes, only: &
    decomp_hr, decomp_hr_vr, gpp, ar, er, supplement_to_sminn, supplement_to_sminn_vr, &
    cropprod1c_loss, grainc_to_cropprodc, grainc_to_seed, &
    sminn_leached, sminn_leached_vr, smin_no3_leached, smin_no3_leached_vr, smin_no3_runoff, smin_no3_runoff_vr, &
    f_n2o_nit, f_n2o_nit_vr, decomp_cpools_transport_tendency, decomp_npools_transport_tendency, &
    denit, f_denit_vr, fire_closs, hrv_xsmrpool_to_atm, som_c_leached, som_n_leached, sminn_to_denit_excess_vr, &
    sminn_to_denit_decomp_vr
use MOD_1D_BGCPFTFluxes, only: &
    psn_to_cpool_p, leaf_mr_p, froot_mr_p, livestem_mr_p, livecroot_mr_p, &
    cpool_leaf_gr_p, cpool_froot_gr_p, cpool_livestem_gr_p, cpool_deadstem_gr_p, &
    cpool_livecroot_gr_p, cpool_deadcroot_gr_p, transfer_leaf_gr_p, transfer_froot_gr_p, &
    transfer_livestem_gr_p, transfer_deadstem_gr_p, &
    transfer_livecroot_gr_p, transfer_deadcroot_gr_p, &
    cpool_leaf_storage_gr_p, cpool_froot_storage_gr_p, &
    cpool_livestem_storage_gr_p, cpool_deadstem_storage_gr_p, &
    cpool_livecroot_storage_gr_p, cpool_deadcroot_storage_gr_p, &
    grain_mr_p, xsmrpool_to_atm_p, cpool_grain_gr_p, &
    transfer_grain_gr_p, cpool_grain_storage_gr_p, soil_change_p, &
    fire_closs_p, hrv_xsmrpool_to_atm_p, &
    m_leafc_to_fire_p, m_leafc_storage_to_fire_p, m_leafc_xfer_to_fire_p, &
    m_frootc_to_fire_p, m_frootc_storage_to_fire_p, m_frootc_xfer_to_fire_p, &
    m_livestemc_to_fire_p, m_livestemc_storage_to_fire_p, m_livestemc_xfer_to_fire_p, &
    m_deadstemc_to_fire_p, m_deadstemc_storage_to_fire_p, m_deadstemc_xfer_to_fire_p, &
    m_livecrootc_to_fire_p, m_livecrootc_storage_to_fire_p, m_livecrootc_xfer_to_fire_p, &
    m_deadcrootc_to_fire_p, m_deadcrootc_storage_to_fire_p, m_deadcrootc_xfer_to_fire_p, &
    m_gresp_storage_to_fire_p, m_gresp_xfer_to_fire_p, &
    cropprod1c_loss_p, grainc_to_seed_p, grainc_to_food_p
    
implicit none

public CNDriverSummarizeStates
public CNDriverSummarizeFluxes

private soilbiogeochem_carbonstate_summary
private soilbiogeochem_nitrogenstate_summary
private cnveg_carbonstate_summary
private cnveg_nitrogenstate_summary
private soilbiogeochem_carbonflux_summary
private soilbiogeochem_nitrogenflux_summary
private cnveg_carbonflux_summary
private cnveg_nitrogenflux_summary

contains

subroutine CNDriverSummarizeStates(i,ps,pe,nl_soil,dz_soi,ndecomp_pools)

integer, intent(in) :: i
integer, intent(in) :: ps
integer, intent(in) :: pe
integer, intent(in) :: nl_soil
real(r8),intent(in) :: dz_soi(1:nl_soil)
integer, intent(in) :: ndecomp_pools

call soilbiogeochem_carbonstate_summary(i,nl_soil,dz_soi,ndecomp_pools)
call soilbiogeochem_nitrogenstate_summary(i,nl_soil,dz_soi,ndecomp_pools)

call cnveg_carbonstate_summary(i,ps,pe)
call cnveg_nitrogenstate_summary(i,ps,pe)

end subroutine CNDriverSummarizeStates

subroutine CNDriverSummarizeFluxes(i,ps,pe,nl_soil,dz_soi,ndecomp_transitions,ndecomp_pools)

integer, intent(IN) :: i
integer, intent(IN) :: ps
integer, intent(IN) :: pe
integer, intent(IN) :: nl_soil
real(r8),intent(IN) :: dz_soi(1:nl_soil)
integer, intent(IN) :: ndecomp_transitions
integer, intent(IN) :: ndecomp_pools

call soilbiogeochem_carbonflux_summary(i,nl_soil,dz_soi,ndecomp_transitions,ndecomp_pools)

call soilbiogeochem_nitrogenflux_summary(i,nl_soil,dz_soi,ndecomp_transitions,ndecomp_pools)

call cnveg_carbonflux_summary(i,ps,pe)

call cnveg_nitrogenflux_summary()

end subroutine CNDriverSummarizeFluxes

subroutine soilbiogeochem_carbonstate_summary(i,nl_soil,dz_soi, ndecomp_pools)

integer, intent(in) :: i
integer, intent(in) :: nl_soil
real(r8),intent(in) :: dz_soi(1:nl_soil)
integer, intent(in) :: ndecomp_pools

integer :: l,j

totsomc(i) = 0._r8
totlitc(i) = 0._r8
totcwdc(i) = 0._r8
ctrunc_soil(i) = 0._r8

do l = 1, ndecomp_pools
   decomp_cpools(l,i) = 0._r8
   do j = 1, nl_soil
      decomp_cpools(l,i) = decomp_cpools(l,i) + decomp_cpools_vr(j,l,i) * dz_soi(j)
   end do
end do
   
do l = 1, ndecomp_pools
   if(is_litter(l))then
      totlitc(i) = totlitc(i) + decomp_cpools(l,i)
   end if
   if(is_soil(l))then
      totsomc(i) = totsomc(i) + decomp_cpools(l,i) 
   end if
   if(is_cwd(l))then
      totcwdc(i) = totcwdc(i) + decomp_cpools(l,i)
   end if
end do

do j = 1, nl_soil
   ctrunc_soil(i) = ctrunc_soil(i) + ctrunc_vr(j,i) * dz_soi(j)
end do

end subroutine soilbiogeochem_carbonstate_summary

subroutine soilbiogeochem_nitrogenstate_summary(i,nl_soil,dz_soi,ndecomp_pools)

integer, intent(in) :: i
integer, intent(in) :: nl_soil
real(r8),intent(in) :: dz_soi(1:nl_soil)
integer, intent(in) :: ndecomp_pools

integer :: l,j

totsomn(i) = 0._r8
totlitn(i) = 0._r8
totcwdn(i) = 0._r8
sminn(i)   = 0._r8
ntrunc_soil(i)  = 0._r8

do l = 1, ndecomp_pools
   decomp_npools(l,i) = 0._r8
   do j = 1, nl_soil
      decomp_npools(l,i) = decomp_npools(l,i) + decomp_npools_vr(j,l,i) * dz_soi(j)
   end do
end do
   
do j = 1, nl_soil
   sminn(i) = sminn(i) + sminn_vr(j,i) * dz_soi(j)
end do

do l = 1, ndecomp_pools
   if(is_litter(l))then
      totlitn(i) = totlitn(i) + decomp_npools(l,i)
   end if
   if(is_soil(l))then
      totsomn(i) = totsomn(i) + decomp_npools(l,i) 
   end if
   if(is_cwd(l))then
      totcwdn(i) = totcwdn(i) + decomp_npools(l,i)
   end if
end do

do j = 1, nl_soil
   ntrunc_soil(i) = ntrunc_soil(i) + ntrunc_vr(j,i) * dz_soi(j)
end do

end subroutine soilbiogeochem_nitrogenstate_summary

subroutine cnveg_carbonstate_summary(i,ps,pe)

integer, intent(in) :: i
integer, intent(in) :: ps
integer, intent(in) :: pe

integer m

leafc(i)              = sum(leafc_p(ps:pe)              * pftfrac(ps:pe))
leafc_storage(i)      = sum(leafc_storage_p(ps:pe)      * pftfrac(ps:pe))
leafc_xfer(i)         = sum(leafc_xfer_p(ps:pe)         * pftfrac(ps:pe))
frootc(i)             = sum(frootc_p(ps:pe)             * pftfrac(ps:pe))
frootc_storage(i)     = sum(frootc_storage_p(ps:pe)     * pftfrac(ps:pe))
frootc_xfer(i)        = sum(frootc_xfer_p(ps:pe)        * pftfrac(ps:pe))
livestemc(i)          = sum(livestemc_p(ps:pe)          * pftfrac(ps:pe))
livestemc_storage(i)  = sum(livestemc_storage_p(ps:pe)  * pftfrac(ps:pe))
livestemc_xfer(i)     = sum(livestemc_xfer_p(ps:pe)     * pftfrac(ps:pe))
deadstemc(i)          = sum(deadstemc_p(ps:pe)          * pftfrac(ps:pe))
deadstemc_storage(i)  = sum(deadstemc_storage_p(ps:pe)  * pftfrac(ps:pe))
deadstemc_xfer(i)     = sum(deadstemc_xfer_p(ps:pe)     * pftfrac(ps:pe))
livecrootc(i)         = sum(livecrootc_p(ps:pe)         * pftfrac(ps:pe))
livecrootc_storage(i) = sum(livecrootc_storage_p(ps:pe) * pftfrac(ps:pe))
livecrootc_xfer(i)    = sum(livecrootc_xfer_p(ps:pe)    * pftfrac(ps:pe))
deadcrootc(i)         = sum(deadcrootc_p(ps:pe)         * pftfrac(ps:pe))
deadcrootc_storage(i) = sum(deadcrootc_storage_p(ps:pe) * pftfrac(ps:pe))
deadcrootc_xfer(i)    = sum(deadcrootc_xfer_p(ps:pe)    * pftfrac(ps:pe))
grainc(i)             = sum(grainc_p(ps:pe)             * pftfrac(ps:pe))
grainc_storage(i)     = sum(grainc_storage_p(ps:pe)     * pftfrac(ps:pe))
grainc_xfer(i)        = sum(grainc_xfer_p(ps:pe)        * pftfrac(ps:pe))
xsmrpool(i)           = sum(xsmrpool_p(ps:pe)           * pftfrac(ps:pe))
cropseedc_deficit(i)  = sum(cropseedc_deficit_p(ps:pe)  * pftfrac(ps:pe))
cropprod1c(i)         = sum(cropprod1c_p(ps:pe)         * pftfrac(ps:pe))
#ifdef CROP
cphase(i)             = sum(cphase_p(ps:pe)             * pftfrac(ps:pe))
!prodc10(i)            = sum(prodc10_p(ps:pe)            * pftfrac(ps:pe))
#endif

fertnitro_corn(i) = 0._r8
fertnitro_swheat(i) = 0._r8
fertnitro_wwheat(i) = 0._r8
fertnitro_soybean(i) = 0._r8
fertnitro_cotton(i) = 0._r8
fertnitro_rice1(i) = 0._r8
fertnitro_rice2(i) = 0._r8
fertnitro_sugarcane(i) = 0._r8
do m = ps, pe
   totvegc_p(m) = leafc_p(m)             + frootc_p(m)             + livestemc_p(m) &
                + deadstemc_p(m)         + livecrootc_p(m)         + deadcrootc_p(m) &
                + leafc_storage_p(m)     + frootc_storage_p(m)     + livestemc_storage_p(m) &
                + deadstemc_storage_p(m) + livecrootc_storage_p(m) + deadcrootc_storage_p(m) &
                + leafc_xfer_p(m)        + frootc_xfer_p(m)        + livestemc_xfer_p(m) &
                + deadstemc_xfer_p(m)    + livecrootc_xfer_p(m)    + deadcrootc_xfer_p(m) &
                + gresp_storage_p(m)     + gresp_xfer_p(m)         + xsmrpool_p(m) &
                + grainc_p(m)            + grainc_storage_p(m)     + grainc_xfer_p(m) &
                + cropseedc_deficit_p(m) + cpool_p(m)

   if( pftclass(m) .eq. 17 .or. pftclass(m) .eq. 18 .or. pftclass(m) .eq. 63 .or. pftclass(m) .eq. 64)then
      fertnitro_corn(i) = fertnitro_p(m) 
   else
      if( pftclass(m) .eq. 19 .or. pftclass(m) .eq. 20)then
         fertnitro_swheat(i) = fertnitro_p(m)
      else
         if(pftclass(m) .eq. 21 .or. pftclass(m) .eq. 22)then
            fertnitro_wwheat(i) = fertnitro_p(m)
         else
            if(pftclass(m) .eq. 23 .or. pftclass(m) .eq. 24 .or. pftclass(m) .eq. 77 .or. pftclass(m) .eq. 78)then
               fertnitro_soybean(i) = fertnitro_p(m)
            else
               if(pftclass(m) .eq. 41 .or. pftclass(m) .eq. 42)then
                  fertnitro_cotton(i) = fertnitro_p(m)
               else
                  if(pftclass(m) .eq. 61 .or. pftclass(m) .eq. 62)then
                     fertnitro_rice1(i) = fertnitro_p(m)
                     fertnitro_rice2(i) = fertnitro_p(m)
                  else
                     if(pftclass(m) .eq. 67 .or. pftclass(m) .eq. 68)then
                        fertnitro_sugarcane(i) = fertnitro_p(m)
                     end if
                  end if
               end if
            end if
         end if
      end if
   end if
end do

totvegc(i) = sum(totvegc_p(ps:pe)*pftfrac(ps:pe))
ctrunc_veg(i) = sum(ctrunc_p(ps:pe) *pftfrac(ps:pe))
totcolc(i) = totvegc(i) + totcwdc(i) + totlitc(i) + totsomc(i) + ctrunc_veg(i) +ctrunc_soil(i)


end subroutine cnveg_carbonstate_summary

subroutine cnveg_nitrogenstate_summary(i,ps,pe)

integer, intent(in) :: i
integer, intent(in) :: ps
integer, intent(in) :: pe

integer m

leafn(i)              = sum(leafn_p(ps:pe)              * pftfrac(ps:pe))
leafn_storage(i)      = sum(leafn_storage_p(ps:pe)      * pftfrac(ps:pe))
leafn_xfer(i)         = sum(leafn_xfer_p(ps:pe)         * pftfrac(ps:pe))
frootn(i)             = sum(frootn_p(ps:pe)             * pftfrac(ps:pe))
frootn_storage(i)     = sum(frootn_storage_p(ps:pe)     * pftfrac(ps:pe))
frootn_xfer(i)        = sum(frootn_xfer_p(ps:pe)        * pftfrac(ps:pe))
livestemn(i)          = sum(livestemn_p(ps:pe)          * pftfrac(ps:pe))
livestemn_storage(i)  = sum(livestemn_storage_p(ps:pe)  * pftfrac(ps:pe))
livestemn_xfer(i)     = sum(livestemn_xfer_p(ps:pe)     * pftfrac(ps:pe))
deadstemn(i)          = sum(deadstemn_p(ps:pe)          * pftfrac(ps:pe))
deadstemn_storage(i)  = sum(deadstemn_storage_p(ps:pe)  * pftfrac(ps:pe))
deadstemn_xfer(i)     = sum(deadstemn_xfer_p(ps:pe)     * pftfrac(ps:pe))
livecrootn(i)         = sum(livecrootn_p(ps:pe)         * pftfrac(ps:pe))
livecrootn_storage(i) = sum(livecrootn_storage_p(ps:pe) * pftfrac(ps:pe))
livecrootn_xfer(i)    = sum(livecrootn_xfer_p(ps:pe)    * pftfrac(ps:pe))
deadcrootn(i)         = sum(deadcrootn_p(ps:pe)         * pftfrac(ps:pe))
deadcrootn_storage(i) = sum(deadcrootn_storage_p(ps:pe) * pftfrac(ps:pe))
deadcrootn_xfer(i)    = sum(deadcrootn_xfer_p(ps:pe)    * pftfrac(ps:pe))
grainn(i)             = sum(grainn_p(ps:pe)             * pftfrac(ps:pe))
grainn_storage(i)     = sum(grainn_storage_p(ps:pe)     * pftfrac(ps:pe))
grainn_xfer(i)        = sum(grainn_xfer_p(ps:pe)        * pftfrac(ps:pe))
retransn(i)           = sum(retransn_p(ps:pe)           * pftfrac(ps:pe))
do m = ps, pe
   totvegn_p(m) = leafn_p(m)             + frootn_p(m)             + livestemn_p(m) &
                + deadstemn_p(m)         + livecrootn_p(m)         + deadcrootn_p(m) &
                + leafn_storage_p(m)     + frootn_storage_p(m)     + livestemn_storage_p(m) &
                + deadstemn_storage_p(m) + livecrootn_storage_p(m) + deadcrootn_storage_p(m) &
                + leafn_xfer_p(m)        + frootn_xfer_p(m)        + livestemn_xfer_p(m) &
                + deadstemn_xfer_p(m)    + livecrootn_xfer_p(m)    + deadcrootn_xfer_p(m) &
                + npool_p(m)             + retransn_p(m) &
                + grainn_p(m)            + grainn_storage_p(m)     + grainn_xfer_p(m) &
                + cropseedn_deficit_p(m)
end do

totvegn(i) = sum(totvegn_p(ps:pe)*pftfrac(ps:pe))
ntrunc_veg(i)  = sum(ntrunc_p(ps:pe) *pftfrac(ps:pe))
totcoln(i) = totvegn(i) + totcwdn(i) + totlitn(i) + totsomn(i) + sminn(i) + ntrunc_veg(i) + ntrunc_soil(i)

end subroutine cnveg_nitrogenstate_summary

subroutine soilbiogeochem_carbonflux_summary(i,nl_soil,dz_soi,ndecomp_transitions,ndecomp_pools)

integer, intent(in) :: i
integer, intent(in) :: nl_soil
real(r8),intent(in) :: dz_soi(1:nl_soil)
integer, intent(in) :: ndecomp_transitions
integer, intent(in) :: ndecomp_pools

integer k,j,l

do k = 1, ndecomp_transitions
   do j = 1, nl_soil
      decomp_hr(i) = decomp_hr(i) &
                   + decomp_hr_vr(j,k,i) * dz_soi(j)
   end do
end do

do l = 1, ndecomp_pools
   do j = 1, nl_soil
      som_c_leached(i) = som_c_leached(i) + decomp_cpools_transport_tendency(j,l,i) * dz_soi(j)
   end do
end do


end subroutine soilbiogeochem_carbonflux_summary

subroutine soilbiogeochem_nitrogenflux_summary(i,nl_soil,dz_soi,ndecomp_transitions,ndecomp_pools)

integer, intent(in) :: i
integer, intent(in) :: nl_soil
real(r8),intent(in) :: dz_soi(1:nl_soil)
integer, intent(in) :: ndecomp_transitions
integer, intent(in) :: ndecomp_pools

integer j,l,k

do l = 1, ndecomp_pools
   do j = 1, nl_soil
      som_n_leached(i) = som_n_leached(i) + decomp_npools_transport_tendency(j,l,i) * dz_soi(j)
   end do
end do

do j = 1, nl_soil
   supplement_to_sminn(i) = supplement_to_sminn(i) + supplement_to_sminn_vr(j,i) * dz_soi(j)
   smin_no3_leached(i)    = smin_no3_leached(i)    + smin_no3_leached_vr(j,i) * dz_soi(j)
   smin_no3_runoff(i)     = smin_no3_runoff(i)     + smin_no3_runoff_vr(j,i) * dz_soi(j)
   sminn_leached(i)       = sminn_leached(i)       + sminn_leached_vr(j,i) * dz_soi(j)
   f_n2o_nit(i)           = f_n2o_nit(i)           + f_n2o_nit_vr(j,i) * dz_soi(j)
   denit(i)               = denit(i)               &
#ifdef NITRIF
                          + f_denit_vr(j,i) * dz_soi(j)
#else
                          + sminn_to_denit_excess_vr(j,i) * dz_soi(j)
#endif

#ifndef NITRIF
   do k = 1, ndecomp_transitions
      denit(i) = denit(i) + sminn_to_denit_decomp_vr(j,k,i) * dz_soi(j)
   end do
#endif
end do
end subroutine soilbiogeochem_nitrogenflux_summary

subroutine cnveg_carbonflux_summary(i,ps,pe)

integer, intent(in) :: i
integer, intent(in) :: ps
integer, intent(in) :: pe

integer m

gpp(i) = sum(psn_to_cpool_p(ps:pe) * pftfrac(ps:pe))
downreg(i) = sum(downreg_p(ps:pe) * pftfrac(ps:pe))
ar (i) = sum((leaf_mr_p(ps:pe)                    + froot_mr_p(ps:pe) &
            + livestem_mr_p(ps:pe)                + livecroot_mr_p(ps:pe) &
            + cpool_leaf_gr_p(ps:pe)              + cpool_froot_gr_p(ps:pe) &
            + cpool_livestem_gr_p(ps:pe)          + cpool_deadstem_gr_p(ps:pe) &
            + cpool_livecroot_gr_p(ps:pe)         + cpool_deadcroot_gr_p(ps:pe) &
            + transfer_leaf_gr_p(ps:pe)           + transfer_froot_gr_p(ps:pe) &
            + transfer_livestem_gr_p(ps:pe)       + transfer_deadstem_gr_p(ps:pe) &
            + transfer_livecroot_gr_p(ps:pe)      + transfer_deadcroot_gr_p(ps:pe) &
            + cpool_leaf_storage_gr_p(ps:pe)      + cpool_froot_storage_gr_p(ps:pe) &  
            + cpool_livestem_storage_gr_p(ps:pe)  + cpool_deadstem_storage_gr_p(ps:pe) &
            + cpool_livecroot_storage_gr_p(ps:pe) + cpool_deadcroot_storage_gr_p(ps:pe) &
            + grain_mr_p(ps:pe)                   + xsmrpool_to_atm_p(ps:pe) &
            + cpool_grain_gr_p(ps:pe)             + transfer_grain_gr_p(ps:pe) &
            + cpool_grain_storage_gr_p(ps:pe))    * pftfrac(ps:pe))
#ifdef FUN
ar(i) = ar(i) + sum(soil_change_p(ps:pe) * pftfrac(ps:pe))
#endif
er(i) = ar(i) + decomp_hr(i)
#ifdef CROP
cropprod1c_loss(i) = sum(cropprod1c_loss_p(ps:pe) * pftfrac(ps:pe))
grainc_to_cropprodc (i) = sum(grainc_to_food_p(ps:pe)  * pftfrac(ps:pe))
grainc_to_seed (i) = sum(grainc_to_seed_p(ps:pe)  * pftfrac(ps:pe))
#endif

do m = ps, pe
   fire_closs_p(m) = m_leafc_to_fire_p(m) &
                   + m_leafc_storage_to_fire_p(m) &
                   + m_leafc_xfer_to_fire_p(m) &
                   + m_frootc_to_fire_p(m) &
                   + m_frootc_storage_to_fire_p(m) &
                   + m_frootc_xfer_to_fire_p(m) &
                   + m_livestemc_to_fire_p(m) &
                   + m_livestemc_storage_to_fire_p(m) &
                   + m_livestemc_xfer_to_fire_p(m) &
                   + m_deadstemc_to_fire_p(m) &
                   + m_deadstemc_storage_to_fire_p(m) &
                   + m_deadstemc_xfer_to_fire_p(m) &
                   + m_livecrootc_to_fire_p(m) &
                   + m_livecrootc_storage_to_fire_p(m) &
                   + m_livecrootc_xfer_to_fire_p(m) &
                   + m_deadcrootc_to_fire_p(m) &
                   + m_deadcrootc_storage_to_fire_p(m) &
                   + m_deadcrootc_xfer_to_fire_p(m) &
                   + m_gresp_storage_to_fire_p(m) &
                   + m_gresp_xfer_to_fire_p(m)
end do

fire_closs(i)          = sum(fire_closs_p(ps:pe)          * pftfrac(ps:pe))
hrv_xsmrpool_to_atm(i) = sum(hrv_xsmrpool_to_atm_p(ps:pe) * pftfrac(ps:pe))


end subroutine cnveg_carbonflux_summary

subroutine cnveg_nitrogenflux_summary()

end subroutine cnveg_nitrogenflux_summary

end module bgc_CNSummaryMod


#endif
