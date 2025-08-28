#include <define.h>

MODULE MOD_DBedrockReadin

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: dbedrock_readin

CONTAINS

   SUBROUTINE dbedrock_readin (dir_landdata, lc_year)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_UserDefFun
   USE MOD_LandPatch
   USE MOD_NetCDFVector
   USE MOD_Vars_Global, only: nl_soil, dz_soi
   USE MOD_Vars_TimeInvariants, only: dbedrock, ibedrock
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year    ! which year of data used
   character(len=256), intent(in) :: dir_landdata

   ! Local Variables
   character(len=256) :: lndname, cyear
   integer  :: ipatch, L

#ifdef SinglePoint
      dbedrock(:) = SITE_dbedrock
#else
      write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/dbedrock/'//trim(cyear)//'/dbedrock_patches.nc'
      CALL ncio_read_vector (lndname, 'dbedrock_patches', landpatch, dbedrock)
#endif

      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch

            L = landpatch%settyp(ipatch)

            IF (L == 0) THEN
               ibedrock(ipatch) = 0
            ELSE

               dbedrock(ipatch) = dbedrock(ipatch) / 100.0 ! from cm to meter
               dbedrock(ipatch) = max(dbedrock(ipatch), dz_soi(1))

               IF (dbedrock(ipatch) > zi_soi(1)) THEN
                  ibedrock(ipatch) = findloc_ud(dbedrock(ipatch)>zi_soi, back=.true.) + 1
               ELSE
                  ibedrock(ipatch) = 1
               ENDIF

            ENDIF

         ENDDO

      ENDIF

   END SUBROUTINE dbedrock_readin

END MODULE MOD_DBedrockReadin
