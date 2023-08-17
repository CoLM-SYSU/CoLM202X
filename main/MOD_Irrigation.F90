#include <define.h>
#ifdef CROP
module MOD_Irrigation

!   DESCRIPTION:
!       This module has all irrigation related subroutines for irrigated crop at either IGBP/USGS or PFT Land type classification and even in the C and N cycle.
    use MOD_Precision
    USE MOD_TimeManager
    USE MOD_Namelist, only: DEF_simulation_time
    ! ,DEF_IRRIGATION_METHOD
    use MOD_Const_Physical, only: tfrz
    use MOD_Const_PFT, only: irrig_crop
    use MOD_Vars_Global, only: irrig_start_time, irrig_max_depth, irrig_threshold_fraction, irrig_min_cphase, irrig_max_cphase, irrig_time_per_day    
    use MOD_Qsadv, only: qsadv
    use MOD_Vars_TimeInvariants, only: &
#ifdef vanGenuchten_Mualem_SOIL_MODEL
        theta_r, alpha_vgm, n_vgm, &
#endif
        porsl, psi0, bsw
    use MOD_Vars_TimeVariables, only : tref, t_soisno, wliq_soisno, irrig_rate, deficit_irrig, sum_irrig, sum_irrig_count, n_irrig_steps_left, &
        tairday, usday, vsday, pairday, rnetday, fgrndday, potential_evapotranspiration
    use MOD_Vars_PFTimeInvariants, only: pftclass
    use MOD_Vars_PFTimeVariables, only: irrig_method_p
    use MOD_BGC_Vars_PFTimeVariables, only: cphase_p
    use MOD_Vars_1DForcing, only: forc_t, forc_frl, forc_psrf, forc_us, forc_vs
    use MOD_Vars_1DFluxes, only: sabg, sabvsun, sabvsha, olrg, fgrnd

    implicit none

    public :: CalIrrigationNeeded
    public :: CalIrrigationApplicationFluxes

    !   local variable
    integer :: irrig_method_drip = 1
    integer :: irrig_method_sprinkler = 2
    integer :: irrig_method_flood = 3
    integer :: irrig_method_paddy = 4

contains
            
    subroutine CalIrrigationNeeded(i,ps,pe,idate,nl_soil,nbedrock,z_soi,dz_soi,deltim,dlon,npcropmin)

        !   DESCRIPTION:
        !   This subroutine is used to calculate how much irrigation needed in each irrigated crop patch
        integer , intent(in) :: i
        integer , intent(in) :: ps, pe
        integer , intent(in) :: idate(3)
        integer , intent(in) :: nl_soil
        integer , intent(in) :: nbedrock
        real(r8), intent(in) :: z_soi(1:nl_soil)
        real(r8), intent(in) :: dz_soi(1:nl_soil)
        real(r8), intent(in) :: deltim
        real(r8), intent(in) :: dlon
        integer , intent(in) :: npcropmin

        ! local 
        integer :: m
        integer :: irrig_nsteps_per_day
        logical :: check_for_irrig 

        ! !   calculate last day potential evapotranspiration 
        ! call CalPotentialEvapotranspiration(i,idate,dlon,deltim)

        !   calculate whether irrigation needed
        call PointNeedsCheckForIrrig(i,ps,pe,idate,deltim,dlon,npcropmin,check_for_irrig)

        !   calculate irrigation needed
        if (check_for_irrig) then
            call CalIrrigationPotentialNeeded(i,ps,pe,nl_soil,nbedrock,z_soi,dz_soi)
            ! call CalIrrigationLimitedNeeded(i,ps,pe)
        end if

        !   calculate irrigation rate kg/m2->mm/s
        if ((check_for_irrig) .and. (deficit_irrig(i) > 0)) then
            irrig_nsteps_per_day = nint(irrig_time_per_day/deltim)
            irrig_rate(i) = deficit_irrig(i)/deltim/irrig_nsteps_per_day
            n_irrig_steps_left(i) = irrig_nsteps_per_day
            sum_irrig(i) = sum_irrig(i) + deficit_irrig(i)
            sum_irrig_count(i) = sum_irrig_count(i) + 1._r8
        end if
        
        ! !   zero irrigation at the end of growing season 
        ! do m = ps, pe
        !     if (cphase_p(m) >= 4._r8) then
        !         sum_irrig(i) = 0._r8
        !         sum_irrig_count(i) = 0._r8
        !     end if
        ! end do
    end subroutine CalIrrigationNeeded


    subroutine CalIrrigationPotentialNeeded(i,ps,pe,nl_soil,nbedrock,z_soi,dz_soi)

        !   DESCRIPTION:
        !   This subroutine is used to calculate how much irrigation needed in each irrigated crop patch without water supply restriction
        integer , intent(in) :: i
        integer , intent(in) :: ps, pe
        integer , intent(in) :: nbedrock
        integer , intent(in) :: nl_soil
        real(r8), intent(in) :: z_soi(1:nl_soil)
        real(r8), intent(in) :: dz_soi(1:nl_soil)

        !   local variables
        integer  :: j
        integer  :: m
        logical  :: reached_max_depth
        real(r8) :: h2osoi_liq_tot
        real(r8) :: h2osoi_liq_target_tot
        real(r8) :: h2osoi_liq_wilting_point_tot
        real(r8) :: h2osoi_liq_saturation_capacity_tot
        real(r8) :: smpswc,smpsfc
        real(r8) :: h2osoi_liq_wilting_point(1:nl_soil)
        real(r8) :: h2osoi_liq_field_capacity(1:nl_soil)
        real(r8) :: h2osoi_liq_saturation_capacity(1:nl_soil)
        real(r8) :: h2osoi_liq_at_threshold

        real(r8) :: smpswc = -1.5e5
        real(r8) :: smpsfc = -3.3e3  

        !   initialize local variables
        reached_max_depth = .false.
        h2osoi_liq_tot = 0._r8
        h2osoi_liq_target_tot = 0._r8
        h2osoi_liq_wilting_point_tot = 0._r8
        h2osoi_liq_saturation_capacity_tot = 0._r8

        ! !   single site initialization
        ! do m = ps, pe
        !     irrig_method_p(m) = DEF_IRRIGATION_METHOD
        ! enddo

!   calculate wilting point and field capacity
        do j = 1, nl_soil
            if (t_soisno(j,i) > tfrz .and. porsl(j,i) >= 1.e-6) then           
#ifdef Campbell_SOIL_MODEL
                h2osoi_liq_wilting_point(j) = 1000.*dz_soi(j)*porsl(j,i)*((smpswc/psi0(j,i))**(-1/bsw(j,i)))
                h2osoi_liq_field_capacity(j) = 1000.*dz_soi(j)*porsl(j,i)*((smpsfc/psi0(j,i))**(-1/bsw(j,i)))
                h2osoi_liq_saturation_capacity(j) = 1000.*dz_soi(j)*porsl(j,i)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
                h2osoi_liq_wilting_point(j) = soil_vliq_from_psi(smpswc, porsl(j,i), theta_r(j,i), psi0(j,i), 5, &
                  (/alpha_vgm(j,i), n_vgm(j,i), L_vgm(j,i), sc_vgm(j,i), fc_vgm(j,i)/))
                h2osoi_liq_wilting_point(j) = 1000.*dz_soi(j)*h2osoi_liq_wilting_point(j)
                h2osoi_liq_field_capacity(j) = soil_vliq_from_psi(smpsfc, porsl(j,i), theta_r(j,i), psi0(j,i), 5, &
                  (/alpha_vgm(j,i), n_vgm(j,i), L_vgm(j,i), sc_vgm(j,i), fc_vgm(j,i)/))
                h2osoi_liq_field_capacity(j) = 1000.*dz_soi(j)*h2osoi_liq_field_capacity(j)
                h2osoi_liq_saturation_capacity(j) = 1000.*dz_soi(j)*porsl(j,i)      
#endif
            end if
        end do 

        !   calculate total irrigation needed in all soil layers
        do m = ps, pe
            do j = 1, nl_soil
                if (.not. reached_max_depth) then
                    if (z_soi(j) > irrig_max_depth) then
                        reached_max_depth = .true.
                    else if (j > nbedrock) then
                        reached_max_depth = .true.
                    else if (t_soisno(j,i) <= tfrz) then
                        reached_max_depth = .true.
                    else 
                        h2osoi_liq_tot = h2osoi_liq_tot + wliq_soisno(j,i)
                        h2osoi_liq_wilting_point_tot = h2osoi_liq_wilting_point_tot + h2osoi_liq_wilting_point(j)
                        if (irrig_method_p(m) == irrig_method_drip .or. irrig_method_p(m) == irrig_method_sprinkler) then
                            h2osoi_liq_target_tot = h2osoi_liq_target_tot + h2osoi_liq_field_capacity(j)
                        !   irrigation threshold at field capacity, but irrigation amount at saturation capacity
                        else if (irrig_method_p(m) == irrig_method_flood) then
                            h2osoi_liq_target_tot = h2osoi_liq_target_tot + h2osoi_liq_field_capacity(j)
                            h2osoi_liq_saturation_capacity_tot = h2osoi_liq_saturation_capacity_tot + h2osoi_liq_saturation_capacity(j)
                        else if (irrig_method_p(m) == irrig_method_paddy) then
                            h2osoi_liq_target_tot = h2osoi_liq_target_tot + h2osoi_liq_saturation_capacity(j)
                        else
                            h2osoi_liq_target_tot = h2osoi_liq_target_tot + h2osoi_liq_field_capacity(j)
                        end if
                    end if 
                end if
            end do 
        end do

        !   calculate irrigation threshold
        deficit_irrig(i) = 0._r8
        h2osoi_liq_at_threshold = h2osoi_liq_wilting_point_tot + irrig_threshold_fraction * (h2osoi_liq_target_tot - h2osoi_liq_wilting_point_tot)

        !   calculate total irrigation
        do m = ps, pe
            if (h2osoi_liq_tot < h2osoi_liq_at_threshold) then
                if (irrig_method_p(m) == irrig_method_sprinkler) then 
                    deficit_irrig(i) = h2osoi_liq_target_tot - h2osoi_liq_tot
                    ! deficit_irrig(i) = h2osoi_liq_target_tot - h2osoi_liq_tot + potential_evapotranspiration(i)
                else if (irrig_method_p(m) == irrig_method_flood) then
                    deficit_irrig(i) = h2osoi_liq_saturation_capacity_tot - h2osoi_liq_tot
                else
                    deficit_irrig(i) = h2osoi_liq_at_threshold - h2osoi_liq_tot
                end if
            else
                deficit_irrig(i) = 0
            end if
        end do

    end subroutine CalIrrigationPotentialNeeded

    subroutine CalIrrigationApplicationFluxes(i,ps,pe,deltim,qflx_irrig_drip,qflx_irrig_sprinkler,qflx_irrig_flood,qflx_irrig_paddy,irrig_flag)
        !   DESCRIPTION:
        !   This subroutine is used to calculate irrigation application fluxes for each irrigated crop patch
        integer , intent(in) :: i
        integer , intent(in) :: ps, pe
        real(r8), intent(in) :: deltim
        integer , intent(in) :: irrig_flag  ! 1 if sprinker, 2 if others 
        real(r8), intent(out):: qflx_irrig_drip,qflx_irrig_sprinkler,qflx_irrig_flood,qflx_irrig_paddy

        integer :: m 

        qflx_irrig_drip = 0._r8
        qflx_irrig_sprinkler = 0._r8
        qflx_irrig_flood = 0._r8
        qflx_irrig_paddy = 0._r8

        ! !   single site initialization
        ! do m = ps, pe
        !     irrig_method_p(m) = DEF_IRRIGATION_METHOD
        ! enddo

        !   add irrigation fluxes to precipitation or land surface
        do m = ps, pe
            if (n_irrig_steps_left(i) > 0) then
                if ((irrig_flag == 1) .and. (irrig_method_p(m) == irrig_method_sprinkler)) then 
                    qflx_irrig_sprinkler = irrig_rate(i)
                    n_irrig_steps_left(i) = n_irrig_steps_left(i) -1
                    deficit_irrig(i) = deficit_irrig(i) - irrig_rate(i)*deltim
                else if (irrig_flag == 2) then
                    if (irrig_method_p(m) == irrig_method_drip) then
                        qflx_irrig_drip = irrig_rate(i)
                    else if (irrig_method_p(m) == irrig_method_flood) then
                        qflx_irrig_flood = irrig_rate(i)
                    else if (irrig_method_p(m) == irrig_method_paddy) then
                        qflx_irrig_paddy = irrig_rate(i)
                    else if ((irrig_method_p(m) /= irrig_method_drip) .and. (irrig_method_p(m) /= irrig_method_sprinkler) &
                        .and. (irrig_method_p(m) /= irrig_method_flood) .and. (irrig_method_p(m) /= irrig_method_paddy)) then
                        qflx_irrig_drip = irrig_rate(i)
                    end if
                    n_irrig_steps_left(i) = n_irrig_steps_left(i) -1
                    deficit_irrig(i) = deficit_irrig(i) - irrig_rate(i)*deltim
                end if
                if (deficit_irrig(i) < 0._r8) then
                    deficit_irrig(i) = 0._r8
                end if
            else 
                irrig_rate(i) = 0._r8
            end if
        end do
    end subroutine CalIrrigationApplicationFluxes
    
    subroutine PointNeedsCheckForIrrig(i,ps,pe,idate,deltim,dlon,npcropmin,check_for_irrig)
        !   DESCRIPTION:
        !   This subroutine is used to calculate whether irrigation needed in each patch
        integer , intent(in) :: i
        integer , intent(in) :: ps, pe
        integer , intent(in) :: idate(3)
        real(r8), intent(in) :: deltim
        real(r8), intent(in) :: dlon
        integer , intent(in) :: npcropmin
        logical , intent(out):: check_for_irrig

        !   local variable
        integer :: m, ivt
        real(r8):: ldate(3)
        real(r8):: seconds_since_irrig_start_time

        do m = ps, pe 
            ivt = pftclass(m)
            if ((ivt >= npcropmin) .and. (irrig_crop(ivt)) .and. &
                (cphase_p(m) >= irrig_min_cphase) .and. (cphase_p(m)<irrig_max_cphase)) then
                if (DEF_simulation_time%greenwich) then
                    call gmt2local(idate, dlon, ldate)
                    seconds_since_irrig_start_time = ldate(3) - irrig_start_time + deltim
                else
                    seconds_since_irrig_start_time = idate(3) - irrig_start_time + deltim
                end if
                if ((seconds_since_irrig_start_time >= 0._r8) .and. (seconds_since_irrig_start_time < deltim)) then
                    check_for_irrig = .true.
                else
                    check_for_irrig = .false.
                end if
            else
                check_for_irrig = .false.
            end if
        end do

    end subroutine PointNeedsCheckForIrrig

    ! subroutine CalPotentialEvapotranspiration(i,idate,dlon,deltim)
    !     !   DESCRIPTION:
    !     !   This subroutine is used to calculate daily potential evapotranspiration
    !     integer , intent(in) :: i
    !     integer , intent(in) :: idate(3)
    !     real(r8), intent(in) :: dlon
    !     real(r8), intent(in) :: deltim

    !     !   local variable
    !     real(r8):: ldate(3)
    !     real(r8):: seconds_since_irrig_start_time
    !     real(r8) :: es,esdT,qs,qsdT     ! saturation vapour pressure
    !     real(r8) :: evsat               ! vapour pressure
    !     real(r8) :: ur                  ! wind speed
    !     real(r8) :: delta               ! slope of saturation vapour pressure curve 
    !     real(r8) :: gamma               ! Psychrometric constant

    !     if (DEF_simulation_time%greenwich) then
    !         call gmt2local(idate, dlon, ldate)
    !         seconds_since_irrig_start_time = ldate(3) - irrig_start_time + deltim
    !     else
    !         seconds_since_irrig_start_time = idate(3) - irrig_start_time + deltim
    !     end if

    !     if (((seconds_since_irrig_start_time-deltim) >= 0) .and. ((seconds_since_irrig_start_time-deltim) < deltim)) then
    !         tairday(i) = (forc_t(i)-tfrz)*deltim/86400
    !         usday(i) = forc_us(i)*deltim/86400
    !         vsday(i) = forc_vs(i)*deltim/86400
    !         pairday(i) = forc_psrf(i)*deltim/86400/1000
    !         rnetday(i) = (sabg(i)+sabvsun(i)+sabvsha(i)-olrg(i)+forc_frl(i))*deltim/1000000
    !         fgrndday(i) = fgrnd(i)*deltim/1000000
    !     else
    !         tairday(i) = tairday(i) + (forc_t(i)-tfrz)*deltim/86400
    !         usday(i) = usday(i) + forc_us(i)*deltim/86400
    !         vsday(i) = vsday(i) + forc_vs(i)*deltim/86400
    !         pairday(i) = pairday(i) + forc_psrf(i)*deltim/86400/1000
    !         rnetday(i) = rnetday(i) + (sabg(i)+sabvsun(i)+sabvsha(i)-olrg(i)+forc_frl(i))*deltim/1000000
    !         fgrndday(i) = fgrndday(i) + fgrnd(i)*deltim/1000000
    !     endif

    !     if ((seconds_since_irrig_start_time >= 0) .and. (seconds_since_irrig_start_time < deltim)) then
    !         call qsadv(tairday(i),pairday(i),es,esdT,qs,qsdT)
    !         if (tairday(i) > 0)then
    !             evsat = 0.611*EXP(17.27*tairday(i)/(tairday(i)+237.3))
    !         else
    !             evsat = 0.611*EXP(21.87*tairday(i)/(tairday(i)+265.5))
    !         endif
    !         ur = max(0.1,sqrt(usday(i)*usday(i)+vsday(i)*vsday(i)))
    !         delta = 4098*evsat/((tairday(i)+237.3)*(tairday(i)+237.3))
    !         gamma = 0.665*0.001*pairday(i)
    !         potential_evapotranspiration(i) = (0.408*delta*(rnetday(i)-fgrndday(i))+gamma*(900/(tairday(i)+273))*ur* &
    !             (evsat-es))/(delta+(gamma*(1+0.34*ur)))
    !     end if
    ! end subroutine CalPotentialEvapotranspiration

end module MOD_Irrigation
#endif
