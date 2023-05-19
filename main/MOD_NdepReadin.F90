#include <define.h>
#ifdef BGC

SUBROUTINE NDEP_readin (year, dir_landdata, isread, iswrite)
! ===========================================================
!
! !DESCRIPTION:
! Read in the Nitrogen deposition data from CLM5.
!
! !REFERENCE:
! Galloway, J.N., et al. 2004. Nitrogen cycles: past, present, and future. Biogeochem. 70:153-226.
!
! !ORIGINAL:
! Created by Xingjie Lu and Shupeng Zhang, 2022
! ===========================================================

      use precision
      use spmd_task
      use ncio_vector
      use mod_landpatch
      use MOD_BGC_Vars_TimeVars,  only: ndep
      use MOD_BGC_Vars_1DFluxes, only: ndep_to_sminn
      use MOD_TimeInvariants

      IMPLICIT NONE

      integer, INTENT(in) ::  year
      character(LEN=256), INTENT(in) :: dir_landdata
      logical, INTENT(in) :: isread
      logical, INTENT(in) :: iswrite

      character(LEN=256) :: cyear
      character(LEN=256) :: landdir, lndname

      integer npatch, m

! added by Xingjie Lu, 11/10/2022

      if(isread)then
         write(cyear,'(i4.4)') year
         landdir = trim(dir_landdata) // '/NDEP/'
         lndname = trim(landdir) // '/NDEP_'//trim(cyear)//'_patches.nc'
         call ncio_read_vector (lndname, 'NDEP_patches',  landpatch, ndep)
      end if
      if (p_is_worker .and. iswrite) then
         if (numpatch > 0) then
            do npatch = 1, numpatch
               m = patchclass(npatch)
               if(m == 0)then
                  ndep_to_sminn(npatch) = 0.
               else
                  ndep_to_sminn(npatch)  = ndep(npatch) / 3600. / 365. / 24. 
               end if
            end do

         ENDIF
      ENDIF 

END SUBROUTINE NDEP_readin
#endif
