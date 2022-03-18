#include <define.h>

SUBROUTINE pct_readin (dir_landdata)

   use precision
   USE GlobalVars
   use spmd_task
   use ncio_vector
   USE mod_landpatch
#ifdef CLMDEBUG 
   USE mod_colm_debug
#endif
#ifdef PFT_CLASSIFICATION
   use mod_landpft
   use MOD_PFTimeInvars, only : pftfrac
#endif
#ifdef PC_CLASSIFICATION
   use mod_landpc
   use MOD_PCTimeInvars, only : pcfrac
#endif
   IMPLICIT NONE

   character(LEN=256), INTENT(in) :: dir_landdata
   ! Local Variables
   character(len=256) :: lndname
   REAL(r8), allocatable :: sumpct (:)
   INTEGER :: npatch, ipatch

#ifdef PFT_CLASSIFICATION
   lndname = trim(dir_landdata)//'/pct_pfts.nc'
   call ncio_read_vector (lndname, 'pct_pfts', landpft, pftfrac) 

#if (defined CROP) 
   lndname = trim(dir_landdata)//'/pct_crops.nc'
   call ncio_read_vector (lndname, 'pct_crops', landpatch, pctcrop) 
#endif

#ifdef CLMDEBUG
   IF (p_is_worker) THEN 
      npatch = count(landpatch%ltyp == 1)
      allocate (sumpct (npatch))

      npatch = 0
      DO ipatch = 1, numpatch
         IF (landpatch%ltyp(ipatch) == 1) THEN
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
      
#ifdef PC_CLASSIFICATION
   lndname = trim(dir_landdata)//'/pct_pcs.nc'
   CALL ncio_read_vector (lndname, 'pct_pcs', N_PFT, landpc, pcfrac)

#ifdef CLMDEBUG
   IF (p_is_worker) THEN 
      allocate (sumpct (numpc))
      sumpct = sum(pcfrac,dim=1)
   ENDIF
   CALL check_vector_data ('Sum PFT pct', sumpct)
#endif
#endif

END SUBROUTINE pct_readin
