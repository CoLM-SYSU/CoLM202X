#include <define.h>

MODULE MOD_LuLccTimeInvars
! -------------------------------
! Created by Hua Yuan, 04/2022
! -------------------------------

  USE precision
  USE GlobalVars
  USE mod_pixelset
  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
  ! for patch time invariant information
  TYPE(pixelset_type)  :: landpatch_
  INTEGER              :: numpatch_
  INTEGER              :: numelm_
  INTEGER, allocatable :: patchclass_    (:)  !index of land cover type
  INTEGER, allocatable :: patchtype_     (:)  !land water type

  ! for PFT_CLASSIFICATION
  INTEGER, allocatable :: pftclass_      (:)  !PFT type
  INTEGER, allocatable :: patch_pft_s_   (:)  !start PFT index of a patch
  INTEGER, allocatable :: patch_pft_e_   (:)  !end PFT index of a patch

  ! for PC_CLASSIFICATION
  INTEGER, allocatable :: patch2pc_      (:)  !projection from patch to PC

  ! for Urban model
  INTEGER, allocatable :: urbclass_      (:)  !urban TYPE
  INTEGER, allocatable :: patch2urban_   (:)  !projection from patch to Urban

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_LuLccTimeInvars
  PUBLIC :: deallocate_LuLccTimeInvars
  PUBLIC :: SAVE_LuLccTimeInvars

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_LuLccTimeInvars
  ! --------------------------------------------------------------------
  ! Allocates memory for LuLcc time invariant variables
  ! --------------------------------------------------------------------

     USE precision
     USE GlobalVars
     USE mod_landpatch
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
     USE MOD_PFTimeVars
     USE mod_landpft
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
     USE MOD_PCTimeVars
     USE mod_landpc
#endif
#ifdef URBAN_MODEL
     USE MOD_UrbanTimeInvars
     USE MOD_UrbanTimeVars
     USE mod_landurban
#endif
     IMPLICIT NONE

     
     allocate (landpatch_%eindex          (numpatch))
     allocate (landpatch_%ipxstt          (numpatch))
     allocate (landpatch_%ipxend          (numpatch))
     allocate (landpatch_%settyp          (numpatch))
     allocate (landpatch_%ielm            (numpatch))
     allocate (landpatch_%xblkgrp         (numpatch))
     allocate (landpatch_%yblkgrp         (numpatch))

     allocate (patchclass_                (numpatch))
     allocate (patchtype_                 (numpatch))

#ifdef PFT_CLASSIFICATION
     allocate (pftclass_                    (numpft))
     allocate (patch_pft_s_               (numpatch))
     allocate (patch_pft_e_               (numpatch))
#endif

#ifdef PC_CLASSIFICATION
     allocate (patch2pc_                  (numpatch))
#endif

#ifdef URBAN_MODEL
     allocate (urbclass_                  (numurban))
     allocate (patch2urban_               (numpatch))
#endif

  END SUBROUTINE allocate_LuLccTimeInvars


  SUBROUTINE SAVE_LuLccTimeInvars

     USE precision
     USE GlobalVars
     USE mod_pixelset
     USE MOD_TimeInvariants
     USE mod_mesh
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
     USE MOD_PFTimeVars
     USE mod_landpft
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
     USE MOD_PCTimeVars
     USE mod_landpc
#endif
#ifdef URBAN_MODEL
     USE MOD_UrbanTimeInvars
     USE MOD_UrbanTimeVars
     USE mod_landurban
#endif
     USE mod_landpatch
     IMPLICIT NONE

     CALL copy_pixelset(landpatch, landpatch_)
     numpatch_          = numpatch
     numelm_            = numelm
     patchclass_    (:) = patchclass    (:)
     patchtype_     (:) = patchtype     (:)

#ifdef PFT_CLASSIFICATION
     pftclass_      (:) = pftclass      (:)
     patch_pft_s_   (:) = patch_pft_s   (:)
     patch_pft_e_   (:) = patch_pft_e   (:)
#endif

#ifdef PC_CLASSIFICATION
     patch2pc_      (:) = patch2pc      (:)
#endif

#ifdef URBAN_MODEL
     urbclass_      (:) = landurban%settyp (:)
     patch2urban_   (:) = patch2urban      (:)
#endif

  END SUBROUTINE SAVE_LuLccTimeInvars


  SUBROUTINE deallocate_LuLccTimeInvars
      USE mod_pixelset
! --------------------------------------------------
! Deallocates memory for LuLcc time invariant variables
! --------------------------------------------------

     CALL landpatch_%forc_free_mem
     deallocate (patchclass_   )
     deallocate (patchtype_    )

#ifdef PFT_CLASSIFICATION
     deallocate (pftclass_     )
     deallocate (patch_pft_s_  )
     deallocate (patch_pft_e_  )
#endif

#ifdef PC_CLASSIFICATION
     deallocate (patch2pc_     )
#endif

#ifdef URBAN_MODEL
     deallocate (urbclass_     )
     deallocate (patch2urban_  )
#endif

  END SUBROUTINE deallocate_LuLccTimeInvars

END MODULE MOD_LuLccTimeInvars
! ---------- EOP ------------
