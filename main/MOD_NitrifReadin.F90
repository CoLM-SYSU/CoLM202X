#include <define.h>

MODULE MOD_NitrifReadin

!-----------------------------------------------------------------------
   USE precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: NITRIF_readin


!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------


   SUBROUTINE NITRIF_readin (time, dir_landdata)
      ! ==========================================================================================================
      ! !DESCRIPTION:
      ! Read in the O2 data for nitrification module. The climatological data of soil O2 was from CLM5 simulation
      !
      ! !ORIGINAL:
      ! Created by Xingjie Lu and Shupeng Zhang, 2022
      ! ==========================================================================================================

      use precision
      use mod_namelist
      use spmd_task
      use ncio_vector
      use mod_landpatch
      use MOD_Vars_TimeInvariants
      use MOD_Vars_TimeVariables
#ifdef CoLMDEBUG
      use mod_colm_debug
#endif

      USE MOD_Vars_Global
      USE MOD_Vars_LCConst
#ifdef PFT_CLASSIFICATION
      USE mod_landpft
      USE MOD_Vars_PFTimeVars
#endif
#ifdef PC_CLASSIFICATION
      USE mod_landpc
      USE MOD_Vars_PCTimeVars
#endif
#ifdef SinglePoint
      USE mod_single_srfdata
#endif

      IMPLICIT NONE

      integer, INTENT(in) ::  time
      character(LEN=256), INTENT(in) :: dir_landdata

      ! Local variables
      integer :: iyear, itime
      character(LEN=256) :: cyear, ctime
      character(LEN=256) :: landdir, lndname
      integer :: m, npatch,nsl
      CHARACTER(LEN=4) ::cx
      REAL(r8),allocatable :: tCONC_O2_UNSAT_tmp(:)
      REAL(r8),allocatable :: tO2_DECOMP_DEPTH_UNSAT_tmp(:)
      ! READ in nitrif

#ifdef NITRIF
      allocate(tCONC_O2_UNSAT_tmp(numpatch))
      allocate(tO2_DECOMP_DEPTH_UNSAT_tmp(numpatch))

      landdir = trim(dir_landdata) // '/nitrif'
      DO nsl = 1, nl_soil
         write(cx,'(i2.2)') nsl
         write(ctime,'(i2.2)') time
         lndname = trim(landdir) // '/CONC_O2_UNSAT_patches_l' // trim(cx)//'_'// trim(ctime) // '.nc'
         call ncio_read_vector (lndname, 'CONC_O2_UNSAT_patches',  landpatch, tCONC_O2_UNSAT_tmp)
         if (p_is_worker) then
            if (numpatch > 0) then
               do npatch = 1, numpatch
                  m = patchclass(npatch)
                  if( m == 0 )then
                     tCONC_O2_UNSAT(nsl,npatch)  = 0.
                  else
                     tCONC_O2_UNSAT(nsl,npatch)  = tCONC_O2_UNSAT_tmp(npatch)
                  endif
                  if (tCONC_O2_UNSAT(nsl,npatch) < 1E-10) then
                     tCONC_O2_UNSAT(nsl,npatch)=0.0
                  endif
               end do

            ENDIF
         ENDIF
      END do

      DO nsl = 1, nl_soil
         write(cx,'(i2.2)') nsl
         write(ctime,'(i2.2)') time
         lndname = trim(landdir) // '/O2_DECOMP_DEPTH_UNSAT_patches_l' // trim(cx)//'_'// trim(ctime) // '.nc'
         call ncio_read_vector (lndname, 'O2_DECOMP_DEPTH_UNSAT_patches',  landpatch, tO2_DECOMP_DEPTH_UNSAT_tmp)
         if (p_is_worker) then
            if (numpatch > 0) then
               do npatch = 1, numpatch
                  m = patchclass(npatch)
                  if( m == 0 )then
                     tO2_DECOMP_DEPTH_UNSAT(nsl,npatch)  = 0.
                  else
                     tO2_DECOMP_DEPTH_UNSAT(nsl,npatch)  = tO2_DECOMP_DEPTH_UNSAT_tmp(npatch)
                  endif
                  if (tO2_DECOMP_DEPTH_UNSAT(nsl,npatch) < 1E-10) then
                     tO2_DECOMP_DEPTH_UNSAT(nsl,npatch)=0.0
                  endif
               end do

            ENDIF
         ENDIF
      END do
      deallocate(tCONC_O2_UNSAT_tmp)
      deallocate(tO2_DECOMP_DEPTH_UNSAT_tmp)
#endif
   END SUBROUTINE NITRIF_readin

END MODULE MOD_NitrifReadin
