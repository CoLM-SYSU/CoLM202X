#include <define.h>
#ifdef CROP
module MOD_Irrigation

!   DESCRIPTION:
!       This module has all irrigation related subroutines for irrigated crop at either IGBP/USGS or PFT Land type classification and even in the C and N cycle.
    
    use MOD_Precision
    USE MOD_TimeManager
    USE MOD_Namelist, only: DEF_simulation_time,DEF_IRRIGATION_METHOD
    use MOD_Const_Physical, only: tfrz
    use MOD_Const_PFT, only: irrig_crop
    use MOD_Vars_Global, only: irrig_start_time, irrig_max_depth, irrig_threshold_fraction, irrig_min_cphase, irrig_max_cphase, irrig_time_per_day    
    use MOD_Vars_TimeInvariants, only: &
#ifdef vanGenuchten_Mualem_SOIL_MODEL
        theta_r, alpha_vgm, n_vgm, &
#endif
        porsl, psi0, bsw
    use MOD_Vars_PFTimeInvariants, only: pftclass
    use MOD_Vars_TimeVariables, only : t_soisno, wliq_soisno, irrig_rate, deficit_irrig, sum_irrig, sum_irrig_count, n_irrig_steps_left, irrig_method_p
    use MOD_BGC_Vars_PFTimeVariables, only: cphase_p

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
        integer :: irrig_nsteps_per_day
        logical :: check_for_irrig 

        ! !   calculate last day potential evapotranspiration 
        ! call CalPotentialEvapotranspiration()

        !   calculate whether irrigation needed
        call PointNeedsCheckForIrrig(i,ps,pe,idate,deltim,dlon,npcropmin,check_for_irrig)
        !   calculate irrigation needed
        if (check_for_irrig) then
            call CalIrrigationPotentialNeeded(i,ps,pe,nl_soil,nbedrock,z_soi,dz_soi)
            ! call CalIrrigationLimitedNeeded(i,ps,pe)
        end if

        !   calculate irrigation rate kg/m2->mm/s
        if (check_for_irrig) then
            irrig_nsteps_per_day = nint(irrig_time_per_day/deltim)
            irrig_rate(i) = deficit_irrig(i)/deltim/irrig_time_per_day
            n_irrig_steps_left(i) = irrig_nsteps_per_day
            sum_irrig(i) = sum_irrig(i) + deficit_irrig(i)
            sum_irrig_count(i) = sum_irrig_count(i) + 1._r8
        end if

    end subroutine CalIrrigationNeeded

!   需要修改，调试潜在蒸散发计算
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

        !   initialize local variables
        reached_max_depth = .false.
        h2osoi_liq_tot = 0._r8
        h2osoi_liq_target_tot = 0._r8
        h2osoi_liq_wilting_point_tot = 0._r8
        h2osoi_liq_saturation_capacity_tot = 0._r8

        !   single site initialization
        do m = ps, pe
            irrig_method_p(m) = DEF_IRRIGATION_METHOD
        enddo

!   calculate wilting point and field capacity
        do j = 1, nl_soil
            if (t_soisno(i,j)>tfrz .and. porsl(i,j)>=1.e-6) then         
#ifdef Campbell_SOIL_MODEL
                smpswc = -1.5e5
                smpsfc = -3.3e3
                h2osoi_liq_wilting_point(j) = 1000.*dz_soi(j)*porsl(j,i)*((smpswc/psi0(j,i))**(-1/bsw(j,i)))
                h2osoi_liq_field_capacity(j) = 1000.*dz_soi(j)*porsl(j,i)*((smpsfc/psi0(j,i))**(-1/bsw(j,i)))
                h2osoi_liq_saturation_capacity(j) = 1000.*dz_soi(j)*porsl(j,i)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
                smpswc = 1.5e5
                smpsfc = 3.3e3
                h2osoi_liq_wilting_point(j) = theta_r(j,i)+(1000.*dz_soi(j)*porsl(j,i)-1000.*dz_soi(j)* &
                    theta_r(j,i))*((1+(alpha_vgm(j,i)*smpswc)**n_vgm(j,i))**(1./n_vgm(j,i)-1.))
                h2osoi_liq_field_capacity(j) = theta_r(j,i)+(1000.*dz_soi(j)*porsl(j,i)-1000.*dz_soi(j)* &
                    theta_r(j,i))*((1+(alpha_vgm(j,i)*smpsfc)**n_vgm(j,i))**(1./n_vgm(j,i)-1.))
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



!   需要调整方式在PFT上
!   需要修改水量输出方式，增加新的水量变量
    subroutine CalIrrigationApplicationFluxes(i,ps,pe,deltim,qflx_irrig_drip,qflx_irrig_sprinkler,qflx_irrig_flood,qflx_irrig_paddy)
        !   DESCRIPTION:
        !   This subroutine is used to calculate irrigation application fluxes for each irrigated crop patch
        integer , intent(in) :: i
        integer , intent(in) :: ps, pe
        real(r8), intent(in) :: deltim
        real(r8), intent(out):: qflx_irrig_drip,qflx_irrig_sprinkler,qflx_irrig_flood,qflx_irrig_paddy

        integer :: m 

        qflx_irrig_drip = 0._r8
        qflx_irrig_sprinkler = 0._r8
        qflx_irrig_flood = 0._r8
        qflx_irrig_paddy = 0._r8

        !   single site initialization
        do m = ps, pe
            irrig_method_p(m) = DEF_IRRIGATION_METHOD
        enddo

        !   add irrigation fluxes to precipitation or land surface
        do m = ps, pe
            if (n_irrig_steps_left(i) > 0) then
                if (irrig_method_p(m) == irrig_method_drip) then
                    qflx_irrig_drip = irrig_rate(i)
                else if (irrig_method_p(m) == irrig_method_sprinkler) then
                    qflx_irrig_drip = irrig_rate(i)
                else if (irrig_method_p(m) == irrig_method_flood) then
                    qflx_irrig_flood = irrig_rate(i)
                else if (irrig_method_p(m) == irrig_method_paddy) then
                    qflx_irrig_paddy = irrig_rate(i)
                end if
                n_irrig_steps_left(i) = n_irrig_steps_left(i) -1
                deficit_irrig(i) = deficit_irrig(i) - irrig_rate(i)*deltim
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

    ! subroutine CalPotentialEvapotranspiration(i,idate,dlon)
    !     integer , intent(in) :: i
    !     integer , intent(in) :: ps, pe
    !     integer , intent(in) :: idate(3)
    !     if (DEF_simulation_time%greenwich) then
    !         call gmt2local(idate, dlon, ldate)
    !         seconds_since_irrig_start_time = ldate(3) - irrig_start_time
    !     else
    !         seconds_since_irrig_start_time = idate(3) - irrig_start_time + deltim
    !     end if
    !     if (seconds_since_irrig_start_time >= 0) .and. (seconds_since_irrig_start_time < deltim) then
    !         potential_evapotranspiration(i) = 0.408*delta*(rnet-fgrnd)+gamma*(900/(t+273))*ur* &
    !             (es-ea)/(delta+(gamma)*(1+0.34*u))*deltim/86400
    !     else
    !         potential_evapotranspiration(i) = potential_evapotranspiration(i) + 0.408*delta*(rnet-fgrnd)* &
    !             deltim+gamma*(900/(tm+273))*ur*(es-ea)/(delta+(gamma)*(1+0.34*ur))*deltim/86400
    !     end if
    ! end subroutine CalPotentialEvapotranspiration

end module MOD_Irrigation
#endif
