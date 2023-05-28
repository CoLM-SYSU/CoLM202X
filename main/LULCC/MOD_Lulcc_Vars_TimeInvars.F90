#include <define.h>

MODULE MOD_LuLcc_Vars_TimeInvars
! -------------------------------
! Created by Hua Yuan, 04/2022
! -------------------------------

  USE MOD_Precision
  USE MOD_Vars_Global
  USE MOD_PixelSet

  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
  ! for patch time invariant information
  TYPE(pixelset_type)  :: landpatch_
  TYPE(pixelset_type)  :: landelm_
  INTEGER              :: numpatch_
  INTEGER              :: numelm_
  INTEGER              :: numpft_
  INTEGER              :: numpc_
  INTEGER              :: numurban_
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

     use MOD_SPMD_Task
     USE MOD_Precision
     USE MOD_Vars_Global
     USE MOD_LandPatch
     USE MOD_Mesh
#ifdef PFT_CLASSIFICATION
     USE MOD_LandPFT
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_LandPC
#endif
#ifdef URBAN_MODEL
     USE MOD_LandUrban
#endif

     IMPLICIT NONE

     IF (p_is_worker) THEN
        IF (numpatch > 0) THEN
           allocate (landpatch_%eindex          (numpatch))
           allocate (landpatch_%ipxstt          (numpatch))
           allocate (landpatch_%ipxend          (numpatch))
           allocate (landpatch_%settyp          (numpatch))
           allocate (landpatch_%ielm            (numpatch))
           allocate (landpatch_%xblkgrp         (numpatch))
           allocate (landpatch_%yblkgrp         (numpatch))

           allocate (landelm_%eindex              (numelm))
           allocate (landelm_%ipxstt              (numelm))
           allocate (landelm_%ipxend              (numelm))
           allocate (landelm_%settyp              (numelm))

           allocate (patchclass_                (numpatch))
           allocate (patchtype_                 (numpatch))

#ifdef PFT_CLASSIFICATION
           IF (numpft > 0) THEN
              allocate (pftclass_                 (numpft))
              allocate (patch_pft_s_            (numpatch))
              allocate (patch_pft_e_            (numpatch))
           ENDIF
#endif

#ifdef PC_CLASSIFICATION
           IF (numpc > 0) THEN
              allocate (patch2pc_               (numpatch))
           ENDIF
#endif

#ifdef URBAN_MODEL
           IF (numurban > 0) THEN
               allocate (urbclass_              (numurban))
               allocate (patch2urban_           (numpatch))
           ENDIF
#endif
        ENDIF
     ENDIF
  END SUBROUTINE allocate_LuLccTimeInvars


  SUBROUTINE SAVE_LuLccTimeInvars

     USE MOD_Precision
     USE MOD_Vars_Global
     use MOD_SPMD_Task
     USE MOD_Pixelset
     USE MOD_Vars_TimeInvariants
     USE MOD_Landpatch
     USE MOD_Landelm
     USE MOD_Mesh
#ifdef PFT_CLASSIFICATION
     USE MOD_Vars_PFTimeInvars
     USE MOD_LandPFT
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_LandPC
#endif
#ifdef URBAN_MODEL
     USE MOD_LandUrban
#endif

     IMPLICIT NONE

     IF (p_is_worker) THEN
        IF (numpatch > 0) THEN
           CALL copy_pixelset(landpatch, landpatch_)
           CALL copy_pixelset(landelm  , landelm_  )
           numpatch_             = numpatch
           numelm_               = numelm
           patchclass_       (:) = patchclass       (:)
           patchtype_        (:) = patchtype        (:)

#ifdef PFT_CLASSIFICATION
           IF (numpft > 0) THEN
              numpft_            = numpft
              pftclass_      (:) = pftclass         (:)
              patch_pft_s_   (:) = patch_pft_s      (:)
              patch_pft_e_   (:) = patch_pft_e      (:)
           ENDIF
#endif

#ifdef PC_CLASSIFICATION
           IF (numpc > 0) THEN
              numpc_             = numpc
              patch2pc_      (:) = patch2pc         (:)
           ENDIF
#endif

#ifdef URBAN_MODEL
           IF (numurban > 0) THEN
              numurban_          = numurban
              urbclass_      (:) = landurban%settyp (:)
              patch2urban_   (:) = patch2urban      (:)
           ENDIF
#endif
        ENDIF
     ENDIF
  END SUBROUTINE SAVE_LuLccTimeInvars


  SUBROUTINE deallocate_LuLccTimeInvars
      use MOD_SPMD_Task
      USE MOD_PixelSet
! --------------------------------------------------
! Deallocates memory for LuLcc time invariant variables
! --------------------------------------------------
     IF (p_is_worker) THEN
        IF (numpatch_ > 0) THEN
           CALL landpatch_%forc_free_mem
           CALL landelm_%forc_free_mem
           deallocate    (patchclass_   )
           deallocate    (patchtype_    )

#ifdef PFT_CLASSIFICATION
           IF (numpft_ > 0) THEN
              deallocate (pftclass_     )
              deallocate (patch_pft_s_  )
              deallocate (patch_pft_e_  )
           ENDIF
#endif

#ifdef PC_CLASSIFICATION
           IF (numpc_ > 0) THEN
              deallocate (patch2pc_     )
           ENDIF
#endif

#ifdef URBAN_MODEL
           IF (numurban_ > 0) THEN
              deallocate (urbclass_     )
              deallocate (patch2urban_  )
           ENDIF
#endif
        ENDIF
     ENDIF

  END SUBROUTINE deallocate_LuLccTimeInvars

END MODULE MOD_LuLcc_Vars_TimeInvars
! ---------- EOP ------------
