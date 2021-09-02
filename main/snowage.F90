
 subroutine snowage ( deltim,tg,scv,scvold,sag )

!=======================================================================
! Original version: Robert Dickinson
! Update snow cover and snow age, based on BATS code
!=======================================================================

  use precision
  use PhysicalConstants, only : tfrz
  implicit none

!-------------------------- Dummy Argument -----------------------------

  real(r8), INTENT(in) :: deltim ! seconds in a time step [second]
  real(r8), INTENT(in) :: tg     ! temperature of soil at surface [K]
  real(r8), INTENT(in) :: scv    ! snow cover, water equivalent [mm]
  real(r8), INTENT(in) :: scvold ! snow cover for previous time step [mm]
  real(r8), INTENT(inout) :: sag ! non dimensional snow age [-]

!-------------------------- Local variables ----------------------------

  real(r8) :: age1   ! snow aging factor due to crystal growth [-]
  real(r8) :: age2   ! snow aging factor due to surface growth [-]
  real(r8) :: age3   ! snow aging factor due to accum of other particles [-]
  real(r8) :: arg    ! temporary variable used in snow age calculation [-]
  real(r8) :: arg2   ! temporary variable used in snow age calculation [-]
  real(r8) :: dela   ! temporary variable used in snow age calculation [-]
  real(r8) :: dels   ! temporary variable used in snow age calculation [-]
  real(r8) :: sge    ! temporary variable used in snow age calculation [-]

!-----------------------------------------------------------------------
      if(scv <= 0.) then
         sag = 0.
!
! Over antarctica
!
      else if (scv > 800.) then
         sag = 0.
!
! Away from antarctica
!
      else
         age3  = 0.3
         arg   = 5.e3*(1./tfrz-1./tg)
         arg2  = min(0.,10.*arg)
         age2  = exp(arg2)
         age1  = exp(arg)
         dela  = 1.e-6*deltim*(age1+age2+age3)
         dels  = 0.1*max(0.0,scv-scvold)
         sge   = (sag+dela)*(1.0-dels)
         sag   = max(0.0,sge)
      end if

 end subroutine snowage
