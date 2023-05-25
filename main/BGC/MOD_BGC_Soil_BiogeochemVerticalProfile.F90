#include <define.h>
#ifdef BGC
module MOD_BGC_Soil_BiogeochemVerticalProfile

  !------------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module calculate soil vertical profile of different C and N inputs, including:
  ! nitrogen fixation, nitrogen deposition, fine root litter, coarse root litter, 
  ! leaf litter and stem litter.
  !
  ! !ORIGINAL:
  ! The Community Land Model version 5.0 (CLM5)
  !
  ! !REVISION:
  ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.


  use precision
  use MOD_Const_PFT, only: rootfr_p
  use MOD_BGC_Vars_TimeVars, only: &
      nfixation_prof, ndep_prof, altmax_lastyear_indx
  use MOD_BGC_Vars_PFTimeVars, only: &
      leaf_prof_p, froot_prof_p, croot_prof_p, stem_prof_p, cinput_rootfr_p
  use MOD_Vars_PFTimeInvars, only: &
      pftclass, pftfrac
  implicit none

  public SoilBiogeochemVerticalProfile

  real(r8), public :: surfprof_exp  = 10. ! how steep profile is for surface components (1/ e_folding depth) (1/m)

contains

  subroutine SoilBiogeochemVerticalProfile(i,ps,pe,nl_soil,nl_soil_full,nbedrock,zmin_bedrock,z_soi,dz_soi)

    integer ,intent(in) :: i                      ! patch index
    integer ,intent(in) :: ps                     ! start pft index
    integer ,intent(in) :: pe                     ! end pft index
    integer ,intent(in) :: nl_soil                ! number of total soil layers
    integer ,intent(in) :: nl_soil_full           ! number of total soil layers plus bedrock layers
    integer ,intent(in) :: nbedrock               ! where bedrock layer starts (ith soil layer)
    real(r8),intent(in) :: zmin_bedrock           ! depth where bedrock layer starts (m)
    real(r8),intent(in) :: z_soi (1:nl_soil_full) ! depth of each soil layer (m)
    real(r8),intent(in) :: dz_soi(1:nl_soil_full) ! thicknesses of each soil layer (m)

    ! !LOCAL VARIABLES:
    real(r8) :: surface_prof(1:nl_soil)
    real(r8) :: surface_prof_tot
    real(r8) :: rootfr_tot
    real(r8) :: col_cinput_rootfr(1:nl_soil_full)
    integer  :: ivt, m
    integer  :: j 
    ! debugging temp variables
    real(r8) :: froot_prof_sum
    real(r8) :: croot_prof_sum
    real(r8) :: leaf_prof_sum
    real(r8) :: stem_prof_sum
    real(r8) :: ndep_prof_sum
    real(r8) :: nfixation_prof_sum
    real(r8) :: delta = 1.e-10

    surface_prof(:) = 0._r8
    do j = 1, nl_soil
       surface_prof(j) = exp(-surfprof_exp * z_soi(j)) / dz_soi(j)
       if (z_soi(j) > zmin_bedrock) then
          surface_prof(j) = 0._r8
       end if
    end do

    ! initialize profiles to zero
    col_cinput_rootfr(:)   = 0._r8
    nfixation_prof   (:,i) = 0._r8
    ndep_prof        (:,i) = 0._r8
    do m = ps , pe
       ivt = pftclass(m)
       leaf_prof_p (:,m)     = 0._r8
       froot_prof_p(:,m)     = 0._r8
       croot_prof_p(:,m)     = 0._r8
       stem_prof_p (:,m)     = 0._r8

       cinput_rootfr_p(:,m)    = 0._r8

       if (ivt /= 0) then
          do j = 1, nl_soil
             cinput_rootfr_p(j,m) = rootfr_p(j,ivt) / dz_soi(j)
          end do

       else
          cinput_rootfr_p(1,m) = 0.
       endif
    end do

    do m = ps , pe
    ! integrate rootfr over active layer of soil column
       rootfr_tot = 0._r8
       surface_prof_tot = 0._r8
       do j = 1, min(max(altmax_lastyear_indx(i), 1), nl_soil)
          rootfr_tot = rootfr_tot + cinput_rootfr_p(j,m) * dz_soi(j)
          surface_prof_tot = surface_prof_tot + surface_prof(j)  * dz_soi(j)
       end do
       if ( (altmax_lastyear_indx(i) > 0) .and. (rootfr_tot > 0._r8) .and. (surface_prof_tot > 0._r8) ) then
       ! where there is not permafrost extending to the surface, integrate the profiles over the active layer
       ! this is equivalnet to integrating over all soil layers outside of permafrost regions
          do j = 1, min(max(altmax_lastyear_indx(i), 1), nl_soil)
             froot_prof_p(j,m) = cinput_rootfr_p(j,m) / rootfr_tot
             croot_prof_p(j,m) = cinput_rootfr_p(j,m) / rootfr_tot

             if (j > nbedrock .and. cinput_rootfr_p(j,m) > 0._r8) then
                write(*,*) 'ERROR: cinput_rootfr_p > 0 in bedrock'
             end if
          ! set all surface processes to shallower profile
             leaf_prof_p(j,m) = surface_prof(j)/ surface_prof_tot
             stem_prof_p(j,m) = surface_prof(j)/ surface_prof_tot
          end do
       else
          ! if fully frozen, or no roots, put everything in the top layer
          froot_prof_p(1,m) = 1./dz_soi(1)
          croot_prof_p(1,m) = 1./dz_soi(1)
          leaf_prof_p(1,m)  = 1./dz_soi(1)
          stem_prof_p(1,m)  = 1./dz_soi(1)
       endif
    end do


    !! aggregate root profile to column
    do m = ps , pe
       do j = 1,nl_soil
          col_cinput_rootfr(j) = col_cinput_rootfr(j) + cinput_rootfr_p(j,m) * pftfrac(m)
       end do
    end do

    ! repeat for column-native profiles: Ndep and Nfix
    rootfr_tot = 0._r8
    surface_prof_tot = 0._r8
    ! redo column ntegration over active layer for column-native profiles
    do j = 1, min(max(altmax_lastyear_indx(i), 1), nl_soil)
       rootfr_tot = rootfr_tot + col_cinput_rootfr(j) * dz_soi(j)
       surface_prof_tot = surface_prof_tot + surface_prof(j) * dz_soi(j)
    end do
    if ( (altmax_lastyear_indx(i) > 0) .and. (rootfr_tot > 0._r8) .and. (surface_prof_tot > 0._r8) ) then
       do j = 1,  min(max(altmax_lastyear_indx(i), 1), nl_soil)
          nfixation_prof(j,i) = col_cinput_rootfr(j) / rootfr_tot
          ndep_prof(j,i) = surface_prof(j)/ surface_prof_tot
       end do
    else
       nfixation_prof(1,i) = 1./dz_soi(1)
       ndep_prof(1,i) = 1./dz_soi(1)
    endif

  ! check to make sure integral of all profiles = 1.
    ndep_prof_sum = 0.
    nfixation_prof_sum = 0.
    do j = 1, nl_soil
       ndep_prof_sum = ndep_prof_sum + ndep_prof(j,i) *  dz_soi(j)
       nfixation_prof_sum = nfixation_prof_sum + nfixation_prof(j,i) *  dz_soi(j)
    end do
    if ( ( abs(ndep_prof_sum - 1._r8) > delta ) .or.  ( abs(nfixation_prof_sum - 1._r8) > delta ) ) then
       print*,'i',i,delta
       write(*,*) 'profile sums:',ndep_prof_sum-1._r8,nfixation_prof_sum-1._r8
       write(*,*) 'altmax_lastyear_indx: ', altmax_lastyear_indx(i)
       write(*,*) 'nfixation_prof: ', nfixation_prof(:,i)
       write(*,*) 'ndep_prof: ', ndep_prof(:,i)
       write(*,*) 'cinput_rootfr: ', cinput_rootfr_p(:,ps:pe)
       write(*,*) 'dz_soi: ', dz_soi(:)
       write(*,*) 'surface_prof: ', surface_prof(:)
       write(*,*) 'p, itype(p) : ', i, pftclass(ps:pe)
       write(*,*) 'cinput_rootfr(p,:): ', cinput_rootfr_p(:,ps:pe)
       write(*,*) 'ERROR: _prof_sum-1>delta'
       call abort()
    endif

    do m = ps , pe
       froot_prof_sum = 0.
       croot_prof_sum = 0.
       leaf_prof_sum = 0.
       stem_prof_sum = 0.
       do j = 1, nl_soil
          froot_prof_sum = froot_prof_sum + froot_prof_p(j,m) *  dz_soi(j)
          croot_prof_sum = croot_prof_sum + croot_prof_p(j,m) *  dz_soi(j)
          leaf_prof_sum = leaf_prof_sum + leaf_prof_p(j,m) *  dz_soi(j)
          stem_prof_sum = stem_prof_sum + stem_prof_p(j,m) *  dz_soi(j)
       end do
       if ( ( abs(froot_prof_sum - 1._r8) > delta ) .or.  ( abs(croot_prof_sum - 1._r8) > delta ) .or. &
            ( abs(stem_prof_sum - 1._r8) > delta ) .or.  ( abs(leaf_prof_sum - 1._r8) > delta ) ) then
          write(*,*) 'profile sums: ', froot_prof_sum, croot_prof_sum, leaf_prof_sum, stem_prof_sum
          write(*,*) ' ERROR: sum-1 > delta'
          call abort
       endif
    end do

  end subroutine SoilBiogeochemVerticalProfile

end module MOD_BGC_Soil_BiogeochemVerticalProfile
#endif
