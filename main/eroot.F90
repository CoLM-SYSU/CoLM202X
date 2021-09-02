
 subroutine eroot (nl_soil,trsmx0,porsl,bsw,psi0,rootfr,&
                   dz_soisno,t_soisno,wliq_soisno,rootr,etrc,rstfac)
                   
!=======================================================================
! effective root fraction and maximum possible transpiration rate
! Original author : Yongjiu Dai, 08/30/2002
!=======================================================================

  use precision
  use PhysicalConstants, only : tfrz
  implicit none
    
!-----------------------Argument-----------------------------------------

  integer, INTENT(in) :: nl_soil            ! upper bound of array

  real(r8), INTENT(in) :: trsmx0            ! max transpiration for moist soil+100% veg.[mm/s]
  real(r8), INTENT(in) :: porsl(1:nl_soil)  ! soil porosity [-]
  real(r8), INTENT(in) :: bsw(1:nl_soil)    ! Clapp-Hornberger "B"
  real(r8), INTENT(in) :: psi0(1:nl_soil)   ! saturated soil suction (mm) (NEGATIVE)
  real(r8), INTENT(in) :: rootfr(1:nl_soil) ! fraction of roots in a layer, 
  real(r8), INTENT(in) :: dz_soisno(1:nl_soil)   ! layer thickness (m)
  real(r8), INTENT(in) :: t_soisno(1:nl_soil)    ! soil/snow skin temperature (K)
  real(r8), INTENT(in) :: wliq_soisno(1:nl_soil) ! liquid water (kg/m2)

  real(r8), INTENT(out) :: rootr(1:nl_soil) ! root resistance of a layer, all layers add to 1
  real(r8), INTENT(out) :: etrc             ! maximum possible transpiration rate (mm h2o/s)
  real(r8), INTENT(out) :: rstfac           ! factor of soil water stress for photosynthesis
                   
!-----------------------Local Variables------------------------------
                   
  real(r8) roota             ! accumulates root resistance factors
  real(r8) rresis(1:nl_soil) ! soil water contribution to root resistance
  real(r8) s_node            ! vol_liq/porosity
  real(r8) smpmax            ! wilting point potential in mm
  real(r8) smp_node          ! matrix potential

  integer i                  ! loop counter

!-----------------------End Variables list---------------------------

      ! transpiration potential(etrc) and root resistance factors (rstfac)

      roota = 1.e-10         ! must be non-zero to begin
      do i = 1, nl_soil

        if(t_soisno(i)>tfrz .and. porsl(i)>=1.e-6)then
           smpmax = -1.5e5
           s_node = max(wliq_soisno(i)/(1000.*dz_soisno(i)*porsl(i)),0.001)
           s_node = min(1., s_node)
           smp_node = max(smpmax, psi0(i)*s_node**(-bsw(i))) 
           rresis(i) =(1.-smp_node/smpmax)/(1.-psi0(i)/smpmax) 
           rootr(i) = rootfr(i)*rresis(i)
           roota = roota + rootr(i)
        else
           rootr(i) = 0.
        endif

      end do

      ! normalize root resistances to get layer contribution to ET
      rootr(:) = rootr(:)/roota
        
      ! determine maximum possible transpiration rate 
      etrc = trsmx0*roota
      rstfac = roota

 end subroutine eroot
