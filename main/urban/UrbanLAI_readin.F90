#include <define.h>

SUBROUTINE UrbanLAI_readin (year, time, dir_landdata)

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
      USE MOD_UrbanTimeVars
      USE MOD_UrbanTimeInvars
      USE ncio_vector

      IMPLICIT NONE

      INTEGER, intent(in) :: year
      INTEGER, intent(in) :: time
      CHARACTER(LEN=256), intent(in) :: dir_landdata

      CHARACTER(LEN=256) :: lndname
      CHARACTER(len=256) :: cyear, ctime
      INTEGER :: u, npatch

      ! READ in Leaf area index and stem area index
      write(ctime,'(i2.2)') time
      write(cyear,'(i4.4)') year

      !TODO: parameter input for time year
      ! done 
      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/urban_LAI_'//trim(ctime)//'.nc'
      call ncio_read_vector (lndname, 'TREE_LAI',  landurban, urb_lai)

      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/urban_SAI_'//trim(ctime)//'.nc'
      call ncio_read_vector (lndname, 'TREE_SAI',  landurban, urb_sai)

      !TODO: usage?
      ! loop for urban atch to assign fraction of green leaf
      IF (p_is_worker) THEN
         DO u = 1, numurban
            npatch = urban2patch(u)
            urb_green(u) = 1.              !fraction of green leaf
            green(npatch)= 1.              !fraction of green leaf
         ENDDO
      ENDIF


END SUBROUTINE UrbanLAI_readin
