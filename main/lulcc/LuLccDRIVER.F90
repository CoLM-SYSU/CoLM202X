#include <define.h>

SUBROUTINE LuLccDRIVER (casename,dir_srfdata,dir_restart,&
                        nam_srfdata,nam_urbdata,idate,greenwich)

!=======================================================================
! PURPOSE:
!   the main subroutine for Land use and land cover change simulation
!
! Created by Hua Yuan, 04/08/2022
!=======================================================================

   USE precision
   USE MOD_LuLccTimeInvars
   USE MOD_LuLccTimeVars
   USE MOD_LuLccTMatrix

   IMPLICIT NONE

   CHARACTER(LEN=256), intent(in) :: casename      !casename name
   CHARACTER(LEN=256), intent(in) :: dir_srfdata   !surface data directory
   CHARACTER(LEN=256), intent(in) :: dir_restart   !case restart data directory
   CHARACTER(LEN=256), intent(in) :: nam_srfdata   !surface data filename
   CHARACTER(LEN=256), intent(in) :: nam_urbdata   !urban data filename

   LOGICAL, intent(in)    :: greenwich   !true: greenwich time, false: local time
   INTEGER, intent(inout) :: idate(3)    !year, julian day, seconds of the starting time

   ! allocate LuLcc memory
   CALL allocate_LuLccTimeInvars
   CALL allocate_LuLccTimeVars

   ! SAVE variables
   CALL SAVE_LuLccTimeInvars
   CALL SAVE_LuLccTimeVars

   ! cold start for LuLcc
   print *, ">>> LULCC: initializing..."
   CALL LuLccInitialize (casename,dir_srfdata,dir_restart,&
                         nam_srfdata,nam_urbdata,idate,greenwich)

   ! simple method for variable recovery
   print *, ">>> LULCC: simple method for variable recovery..."
   CALL REST_LuLccTimeVars

   ! conserved method for variable revocery
   print *, ">>> LULCC: Mass&Energy conserve for variable recovery..."
   !CALL READ_LuLccTMatrix()
   !CALL LuLccEnergyConserve()
   !CALL LuLccWaterConserve()

   ! deallocate LuLcc memory
   CALL deallocate_LuLccTimeInvars()
   CALL deallocate_LuLccTimeVars()

END SUBROUTINE LuLccDRIVER
