#include <define.h>

MODULE MOD_PercentagesPFTReadin

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: pct_readin


!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------


   SUBROUTINE pct_readin (dir_landdata, lc_year)

      use MOD_Precision
      USE MOD_Vars_Global
      use MOD_SPMD_Task
      USE MOD_NetCDFVector
      USE MOD_LandPatch
#ifdef RangeCheck
      USE MOD_RangeCheck
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      use MOD_LandPFT
      use MOD_Vars_PFTimeInvariants, only : pftfrac
#endif
#ifdef SinglePoint
      USE MOD_SingleSrfdata
#endif
      IMPLICIT NONE

      INTEGER, intent(in) :: lc_year
      character(LEN=256), INTENT(in) :: dir_landdata
      ! Local Variables
      character(len=256) :: lndname, cyear
      REAL(r8), allocatable :: sumpct (:)
      INTEGER :: npatch, ipatch

      write(cyear,'(i4.4)') lc_year
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
#ifndef SinglePoint
      lndname = trim(dir_landdata)//'/pctpft/'//trim(cyear)//'/pct_pfts.nc'
      call ncio_read_vector (lndname, 'pct_pfts', landpft, pftfrac)
#else
      pftfrac = pack(SITE_pctpfts, SITE_pctpfts > 0.)
#endif

#if (defined CROP)
#ifndef SinglePoint
      lndname = trim(dir_landdata)//'/pctpft/'//trim(cyear)//'/pct_crops.nc'
      call ncio_read_vector (lndname, 'pct_crops', landpatch, pctcrop)
#else
      allocate (pctcrop (numpatch))
      IF (SITE_landtype == CROPLAND) THEN
         pctcrop = pack(SITE_pctcrop, SITE_pctcrop > 0.)
      ELSE
         pctcrop = 0.
      ENDIF
#endif
#endif

#ifdef RangeCheck
      IF (p_is_worker) THEN
         npatch = count(patchtypes(landpatch%settyp) == 0)
         allocate (sumpct (npatch))

         npatch = 0
         DO ipatch = 1, numpatch
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
               npatch = npatch + 1
               sumpct(npatch) = sum(pftfrac(patch_pft_s(ipatch):patch_pft_e(ipatch)))
            ENDIF
         ENDDO

      ENDIF

      CALL check_vector_data ('Sum PFT pct', sumpct)
#if (defined CROP)
      CALL check_vector_data ('CROP pct', pctcrop)
#endif

#endif
#endif

      IF (allocated(sumpct)) deallocate(sumpct)

   END SUBROUTINE pct_readin

END MODULE MOD_PercentagesPFTReadin
