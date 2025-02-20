#include <define.h>

MODULE MOD_SoilTextureReadin

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: soiltext_readin

CONTAINS

   SUBROUTINE soiltext_readin (dir_landdata, lc_year)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_UserDefFun
   USE MOD_LandPatch
   USE MOD_NetCDFVector
   USE MOD_Vars_TimeInvariants, only: soiltext
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

   IMPLICIT NONE

   character(len=256), intent(in) :: dir_landdata
   integer, intent(in) :: lc_year    ! which year of land cover data used

   ! Local Variables
   character(len=256) :: lndname
   character(len=4)   :: cyear
   integer  :: ipatch

#ifdef SinglePoint
      soiltext(:) = SITE_soil_texture
#else
      write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)// '/soil/' // trim(cyear) // '/soiltexture_patches.nc'
      CALL ncio_read_vector (lndname, 'soiltext_patches', landpatch, soiltext)
#endif

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            WHERE (soiltext < 0 )  soiltext = 0
            WHERE (soiltext > 12)  soiltext = 0
         ENDIF
      ENDIF

   END SUBROUTINE soiltext_readin

END MODULE MOD_SoilTextureReadin

