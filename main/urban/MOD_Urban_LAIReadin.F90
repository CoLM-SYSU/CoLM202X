#include <define.h>

#ifdef URBAN_MODEL
MODULE MOD_Urban_LAIReadin

  USE precision
  IMPLICIT NONE
  SAVE

  PUBLIC :: UrbanLAI_readin

CONTAINS

 SUBROUTINE UrbanLAI_readin (year, time, dir_landdata)

! ===========================================================
! Read in urban LAI, SAI and urban tree cover data
! ===========================================================

      USE precision
      USE mod_namelist
      USE spmd_task
      USE MOD_Vars_Global
      USE MOD_Vars_LCConst
      USE mod_landurban
      USE MOD_Vars_TimeVariables
      USE MOD_Vars_TimeInvariants
      USE MOD_Urban_Vars_TimeInvars
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

      !TODO-done: parameter input for time year
      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/urban_LAI_'//trim(ctime)//'.nc'
      call ncio_read_vector (lndname, 'TREE_LAI',  landurban, urb_lai)

      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/urban_SAI_'//trim(ctime)//'.nc'
      call ncio_read_vector (lndname, 'TREE_SAI',  landurban, urb_sai)

      ! loop for urban atch to assign fraction of green leaf
      IF (p_is_worker) THEN
         DO u = 1, numurban
            npatch = urban2patch(u)
            tlai(npatch) = urb_lai(u)
            tsai(npatch) = urb_sai(u)
            urb_green(u) = 1.              !TODO: usage? fraction of green leaf
            green(npatch)= 1.              !fraction of green leaf
         ENDDO
      ENDIF

 END SUBROUTINE UrbanLAI_readin

END MODULE MOD_Urban_LAIReadin
#endif
