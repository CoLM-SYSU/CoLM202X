#include <define.h>
#ifdef BGC

module MOD_BGC_Soil_BiogeochemNLeaching

  !----------------------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module calculates the soil mineral N loss due to leaching. The leaching flux is a function of
  ! dissolved N concentration and sub-surface drainage flux.
  !
  ! !ORIGINAL:
  ! The Community Land Model version 5.0 (CLM5.0)
  !
  ! !REVISION:
  ! Xingjie Lu, 2021, revised original CLM5 code to be compatible with CoLM code structure.

  use MOD_Precision
  use MOD_BGC_Vars_TimeInvars, only: sf, sf_no3
  use MOD_Vars_TimeVariables, only: wliq_soisno
  use MOD_BGC_Vars_TimeVars,   only: sminn_vr, smin_no3_vr
  use MOD_Vars_1DFluxes,     only: rnof, rsur
  use MOD_BGC_Vars_1DFluxes, only: &
      sminn_leached_vr, smin_no3_leached_vr, smin_no3_runoff_vr

  implicit none

  public SoilBiogeochemNLeaching

contains

  subroutine SoilBiogeochemNLeaching(i,deltim,nl_soil,zi_soi,dz_soi)

    integer ,intent(in) :: i                 ! patch index
    real(r8),intent(in) :: deltim            ! time step in seconds
    integer ,intent(in) :: nl_soil           ! number of total soil layers
    real(r8),intent(in) :: zi_soi(0:nl_soil) ! interface level below a zsoi level (m)
    real(r8),intent(in) :: dz_soi(1:nl_soil) ! thicknesses of each soil layer (m)

    integer  :: j                                      ! indices
    real(r8) :: disn_conc                              ! dissolved mineral N concentration (gN/kg water)
    real(r8) :: tot_water                              ! total column liquid water (kg water/m2)
    real(r8) :: surface_water                          ! liquid water to shallow surface depth (kg water/m2)
    real(r8) :: drain_tot                              ! total drainage flux (mm H2O /s)
    real(r8), parameter :: depth_runoff_Nloss = 0.05   ! (m) depth over which runoff mixes with soil water for N loss to runoff

    ! calculate the total soil water
    tot_water = 0._r8
    do j = 1,nl_soil
       tot_water = tot_water + wliq_soisno(j,i)
    end do

    ! for runoff calculation; calculate total water to a given depth
    surface_water = 0._r8
    do j = 1,nl_soil
       if ( zi_soi(j) <= depth_runoff_Nloss)  then
          surface_water = surface_water + wliq_soisno(j,i)
       elseif ( zi_soi(j-1) < depth_runoff_Nloss)  then
          surface_water = surface_water + wliq_soisno(j,i) * ( (depth_runoff_Nloss - zi_soi(j-1)) / dz_soi(j))
       endif
    end do

    ! Loop through columns
    drain_tot = rnof(i) - rsur(i)


#ifndef NITRIF
    !----------------------------------------
    ! --------- NITRIF_NITRIF OFF------------
    !----------------------------------------
    do j = 1,nl_soil
       ! calculate the dissolved mineral N concentration (gN/kg water)
       ! assumes that 10% of mineral nitrogen is soluble
       disn_conc = 0._r8
       if (wliq_soisno(j,i) > 0._r8) then
          disn_conc = (sf * sminn_vr(j,i) * dz_soi(j) )/(wliq_soisno(j,i) )
       end if

       ! calculate the N leaching flux as a function of the dissolved
       ! concentration and the sub-surface drainage flux
       if(tot_water > 0._r8)then
          sminn_leached_vr(j,i) = disn_conc * drain_tot * wliq_soisno(j,i) / ( tot_water * dz_soi(j) )
       else
          sminn_leached_vr(j,i) = 0._r8
       end if

          ! limit the flux based on current sminn state
          ! only let at most the assumed soluble fraction
          ! of sminn be leached on any given timestep
       sminn_leached_vr(j,i) = min(sminn_leached_vr(j,i), (sf * sminn_vr(j,i))/deltim)

          ! limit the flux to a positive value
       sminn_leached_vr(j,i) = max(sminn_leached_vr(j,i), 0._r8)

    end do

#else

    !----------------------------------------
    ! --------- NITRIF_NITRIF ON-------------
    !----------------------------------------

    do j = 1,nl_soil
       ! calculate the dissolved mineral N concentration (gN/kg water)
       ! assumes that 10% of mineral nitrogen is soluble
       disn_conc = 0._r8
       if (wliq_soisno(j,i) > 0._r8) then
          disn_conc = (sf_no3 * smin_no3_vr(j,i) * dz_soi(j) )/(wliq_soisno(j,i) )
       end if
       !
       ! calculate the N leaching flux as a function of the dissolved
       ! concentration and the sub-surface drainage flux
       smin_no3_leached_vr(j,i) = disn_conc * drain_tot * wliq_soisno(j,i) / ( tot_water * dz_soi(j) )
       !
       ! ensure that leaching rate isn't larger than soil N pool
       smin_no3_leached_vr(j,i) = min(smin_no3_leached_vr(j,i), smin_no3_vr(j,i) / deltim )
       !
       ! limit the leaching flux to a positive value
       smin_no3_leached_vr(j,i) = max(smin_no3_leached_vr(j,i), 0._r8)
       !
       !
       ! calculate the N loss from surface runoff, assuming a shallow mixing of surface waters into soil and removal based on runoff
       if ( zi_soi(j) <= depth_runoff_Nloss )  then
          smin_no3_runoff_vr(j,i) = disn_conc * rsur(i) * &
               wliq_soisno(j,i) / ( surface_water * dz_soi(j) )
       elseif ( zi_soi(j-1) < depth_runoff_Nloss )  then
          smin_no3_runoff_vr(j,i) = disn_conc * rsur(i) * &
               wliq_soisno(j,i) * ((depth_runoff_Nloss - zi_soi(j-1)) / &
               dz_soi(j)) / ( surface_water * (depth_runoff_Nloss-zi_soi(j-1) ))
       else
          smin_no3_runoff_vr(j,i) = 0._r8
       endif
       !
       ! ensure that runoff rate isn't larger than soil N pool
       smin_no3_runoff_vr(j,i) = min(smin_no3_runoff_vr(j,i), smin_no3_vr(j,i) / deltim - smin_no3_leached_vr(j,i))
       !
       ! limit the flux to a positive value
       smin_no3_runoff_vr(j,i) = max(smin_no3_runoff_vr(j,i), 0._r8)


       ! limit the flux based on current smin_no3 state
       ! only let at most the assumed soluble fraction
       ! of smin_no3 be leached on any given timestep
       smin_no3_leached_vr(j,i) = min(smin_no3_leached_vr(j,i), (sf_no3 * smin_no3_vr(j,i))/deltim)

       ! limit the flux to a positive value
       smin_no3_leached_vr(j,i) = max(smin_no3_leached_vr(j,i), 0._r8)

    end do
#endif

  end subroutine SoilBiogeochemNLeaching

end module MOD_BGC_Soil_BiogeochemNLeaching
#endif
