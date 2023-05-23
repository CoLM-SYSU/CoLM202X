#include <define.h>

SUBROUTINE pct_readin (dir_landdata)

   use MOD_Precision
   USE GlobalVars
   use MOD_SPMD_Task
   use MOD_NetCDFVector
   USE mod_landpatch
#ifdef CoLMDEBUG 
   USE MOD_CoLMDebug
#endif
#ifdef PFT_CLASSIFICATION
   use mod_landpft
   use MOD_Vars_PFTimeInvars, only : pftfrac
#endif
#ifdef PC_CLASSIFICATION
   use mod_landpc
   use MOD_Vars_PCTimeInvars, only : pcfrac
#endif
   IMPLICIT NONE

   character(LEN=256), INTENT(in) :: dir_landdata
   ! Local Variables
   character(len=256) :: lndname
   REAL(r8), allocatable :: sumpct (:)
   INTEGER :: npatch, ipatch

#ifdef PFT_CLASSIFICATION
#ifndef SinglePoint
   lndname = trim(dir_landdata)//'/pctpft/pct_pfts.nc'
   call ncio_read_vector (lndname, 'pct_pfts', landpft, pftfrac) 
#else
   pftfrac = pack(SITE_pctpfts, SITE_pctpfts > 0.)
#endif

#if (defined CROP) 
#ifndef SinglePoint
   lndname = trim(dir_landdata)//'/pctpft/pct_crops.nc'
   call ncio_read_vector (lndname, 'pct_crops', landpatch, pctcrop) 
#else
   allocate (pctcrop (numpatch))
   pctcrop = pack(SITE_pctcrop, SITE_pctcrop > 0.)
#endif
#endif

#ifdef CoLMDEBUG
   IF (p_is_worker) THEN 
      npatch = count(landpatch%settyp == 1)
      allocate (sumpct (npatch))

      npatch = 0
      DO ipatch = 1, numpatch
         IF (landpatch%settyp(ipatch) == 1) THEN
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
#ifndef SinglePoint
   lndname = trim(dir_landdata)//'/pctpft/pct_pcs.nc'
   CALL ncio_read_vector (lndname, 'pct_pcs', N_PFT, landpc, pcfrac)
#else
   pcfrac(:,1) = SITE_pctpfts 
#endif

#ifdef CoLMDEBUG
   IF (p_is_worker) THEN 
      allocate (sumpct (numpc))
      sumpct = sum(pcfrac,dim=1)
   ENDIF
   CALL check_vector_data ('Sum PFT pct', sumpct)
#endif
#endif

   IF (allocated(sumpct)) deallocate(sumpct)

END SUBROUTINE pct_readin
