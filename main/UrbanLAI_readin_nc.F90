#include <define.h>

SUBROUTINE UrbanLAI_readin_nc (year, time, dir_landdata)!(year, month, dir_srfdata, nam_urbdata)

! ===========================================================
! Read in urban LAI, SAI and urban tree cover data
! ===========================================================

      USE precision
      USE mod_namelist
      USE spmd_task
      USE GlobalVars
      USE LC_Const
      USE mod_landurban
      USE MOD_TimeVariables
      USE MOD_TimeInvariants
      USE MOD_UrbanTimeInvars
      USE ncio_vector

      IMPLICIT NONE

      INTEGER, intent(in) :: year
      INTEGER, intent(in) :: time
      ! INTEGER, intent(in) :: month
      ! CHARACTER(LEN=256), intent(in) :: dir_srfdata
      ! CHARACTER(LEN=256), intent(in) :: nam_urbdata
      CHARACTER(LEN=256), intent(in) :: dir_landdata

      CHARACTER(LEN=256) :: lndname
      CHARACTER(len=256) :: cyear, ctime
      INTEGER :: ncid
      INTEGER :: urbantreelai_vid, urbantreesai_vid
      INTEGER :: i, j, t, u, npatch

! READ in Leaf area index and stem area index
      write(ctime,'(i2.2)') time

      !TODO: parameter input for time year
      lndname = trim(dir_landdata)//'/urban/2005/urban_LAI_'//trim(ctime)//'.nc'
      call ncio_read_vector (lndname, 'TREE_LAI',  landurban, urb_lai)

      lndname = trim(dir_landdata)//'/urban/2005/urban_SAI_'//trim(ctime)//'.nc'
      call ncio_read_vector (lndname, 'TREE_SAI',  landurban, urb_sai)

      !TODO: usage?
      IF (p_is_worker) THEN
         DO u = 1, numurban
            ! urb_lai  (u) = urbantreelai(u) !leaf area index
            ! urb_sai  (u) = urbantreesai(u) !stem are index
            urb_green(u) = 1.              !fraction of green leaf
         ENDDO
      ENDIF


END SUBROUTINE UrbanLAI_readin_nc
