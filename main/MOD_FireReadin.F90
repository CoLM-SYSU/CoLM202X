#include <define.h>
#ifdef BGC

MODULE MOD_FireReadin

!-----------------------------------------------------------------------
   USE precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Fire_readin


!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE Fire_readin (year, dir_landdata)
      ! ======================================================================================================
      !
      ! !DESCRIPTION:
      ! Read in the Fire data from CLM5 dataset (month when crop fire peak (abm), GDP, peatland fraction (peatf),
      ! and population density
      !
      ! !ORIGINAL: Xingjie Lu and Shupeng Zhang, 2022
      ! ======================================================================================================

      use precision
      use mod_namelist
      use spmd_task
      use ncio_vector
      use mod_landpatch
      use MOD_Vars_TimeInvariants, only: abm_lf, gdp_lf, peatf_lf
      use MOD_Vars_TimeVariables,  only: hdm_lf
#ifdef CoLMDEBUG
      use mod_colm_debug
#endif

      USE MOD_Vars_Global
      USE MOD_Const_LC
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

      integer, INTENT(in) :: year
      character(LEN=*), INTENT(in) :: dir_landdata
      real(r8),allocatable :: abm_tmp   (:)
      real(r8),allocatable :: gdp_tmp   (:)
      real(r8),allocatable :: peatf_tmp (:)
      real(r8),allocatable :: hdm_tmp   (:)

      ! Local variables
      integer :: iyear, itime
      character(LEN=256) :: cyear, ctime
      character(LEN=256) :: landdir, lndname
      integer :: m, npatch,ns, ipft
      CHARACTER(LEN=4) ::cx
      ! READ in crops

#ifdef Fire
      allocate(abm_tmp     (numpatch))
      allocate(gdp_tmp     (numpatch))
      allocate(peatf_tmp   (numpatch))
      allocate(hdm_tmp     (numpatch))

      landdir = trim(dir_landdata)//'/FIRE'
      lndname = trim(landdir)//'/abm_patches.nc'
      call ncio_read_vector (lndname, 'abm_patches'    , landpatch, abm_tmp)
      IF (p_is_worker) then
         IF (numpatch > 0) then
            DO npatch = 1, numpatch
               abm_lf    (npatch) = int(abm_tmp(npatch))
            ENDDO
         ENDIF
      ENDIF

      landdir = trim(dir_landdata) // '/FIRE'
      lndname = trim(landdir) // '/gdp_patches.nc'
      call ncio_read_vector (lndname, 'gdp_patches'    , landpatch, gdp_tmp)
      IF (p_is_worker) then
         IF (numpatch > 0) then
            DO npatch = 1, numpatch
               gdp_lf    (npatch) = gdp_tmp(npatch)
            ENDDO
         ENDIF
      ENDIF

      landdir = trim(dir_landdata) // '/FIRE'
      lndname = trim(landdir) // '/peatf_patches.nc'
      call ncio_read_vector (lndname, 'peatf_patches'    , landpatch, peatf_tmp)
      IF (p_is_worker) then
         IF (numpatch > 0) then
            DO npatch = 1, numpatch
               peatf_lf    (npatch) = peatf_tmp(npatch)
            ENDDO
         ENDIF
      ENDIF

      landdir = trim(dir_landdata) // '/FIRE'
      write(cyear,'(i4.4)') year
      lndname = trim(landdir) // '/hdm_'//trim(cyear)//'_patches.nc'
      call ncio_read_vector (lndname, 'hdm_patches'    , landpatch, hdm_tmp)
      IF (p_is_worker) then
         IF (numpatch > 0) then
            DO npatch = 1, numpatch
               hdm_lf    (npatch) = hdm_tmp(npatch)
            ENDDO
         ENDIF
      ENDIF

      deallocate(abm_tmp   )
      deallocate(gdp_tmp   )
      deallocate(peatf_tmp )
      deallocate(hdm_tmp   )

#endif
   END SUBROUTINE Fire_readin

END MODULE MOD_FireReadin
#endif
