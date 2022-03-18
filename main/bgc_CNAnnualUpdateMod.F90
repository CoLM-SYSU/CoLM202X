module bgc_CNAnnualUpdateMod

   use MOD_PFTimeVars, only: &
       tempsum_potential_gpp_p, tempmax_retransn_p, tempavg_tref_p, tempsum_npp_p, tempsum_litfall_p, &
       annsum_potential_gpp_p , annmax_retransn_p , annavg_tref_p , annsum_npp_p , annsum_litfall_p 

   use timemanager, only: isendofyear
   use precision

implicit none

public CNAnnualUpdate

contains

 subroutine CNAnnualUpdate(i,ps,pe,deltim,idate)

  !
  ! !DESCRIPTION:
  ! On the radiation time step, update annual summation variables
  !
  ! !USES:
  !
  ! !ARGUMENTS:
  !
  integer ,intent(in) :: i
  integer ,intent(in) :: ps
  integer ,intent(in) :: pe
  real(r8),intent(in) :: deltim
  integer ,intent(in) :: idate(3)

  ! !LOCAL VARIABLES:
  integer m
!  real(r8):: secspyear
!  logical :: end_of_year ! whether each column has reached the end of the year, according to its own annsum_counter
  !-----------------------------------------------------------------------

!  if(isleapyear(idate(1)))then
!     secspyear = 86400._r8 * 366
!  else
!     secspyear = 86400._r8 * 365
!  end if

! annsum_counter(i) = annsum_counter(i) + deltim
! if (annsum_counter(i) >= secspyear) then
!     end_of_year = .true.
!     annsum_counter(i) = 0._r8
!  else
!     end_of_year = .false.
!  end if

  if (isendofyear(idate,deltim)) then

     do m = ps, pe
     ! update annual plant ndemand accumulator
        annsum_potential_gpp_p(m)  = tempsum_potential_gpp_p(m)
        tempsum_potential_gpp_p(m) = 0._r8

     ! update annual total N retranslocation accumulator
        annmax_retransn_p(m)  = tempmax_retransn_p(m)
        tempmax_retransn_p(m) = 0._r8

     ! update annual average 2m air temperature accumulator
        annavg_tref_p(m)  = tempavg_tref_p(m)
        tempavg_tref_p(m) = 0._r8

     ! update annual NPP accumulator, convert to annual total
        annsum_npp_p(m)  = tempsum_npp_p(m) * deltim
        tempsum_npp_p(m) = 0._r8

     ! update annual litfall accumulator, convert to annual total
        annsum_litfall_p(m)  = tempsum_litfall_p(m) * deltim
        tempsum_litfall_p(m) = 0._r8
     end do

  end if

 end subroutine CNAnnualUpdate

end module bgc_CNAnnualUpdateMod
