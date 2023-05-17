#include <define.h>

MODULE MOD_LuLccTimeInvars
! -------------------------------
! Created by Hua Yuan, 04/2022
! -------------------------------

  USE precision
  USE GlobalVars
  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
  ! for patch time invariant information
  INTEGER, allocatable :: patchclass_    (:)  !index of land cover type
  INTEGER, allocatable :: patchtype_     (:)  !land water type
  INTEGER, allocatable :: grid_patch_s_(:,:)  !start patch number of grid
  INTEGER, allocatable :: grid_patch_e_(:,:)  !end patch number of grid

  ! for PFT_CLASSIFICATION
  INTEGER, allocatable :: pftclass_      (:)  !PFT type
  INTEGER, allocatable :: patch_pft_s_   (:)  !start PFT index of a patch
  INTEGER, allocatable :: patch_pft_e_   (:)  !end PFT index of a patch

  ! for PC_CLASSIFICATION
  INTEGER, allocatable :: patch2pc_      (:)  !projection from patch to PC

  ! for Urban model
  INTEGER, allocatable :: urbclass_      (:)  !urban TYPE
  INTEGER, allocatable :: patch2urb_     (:)  !projection from patch to Urban

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
     IMPLICIT NONE

     allocate (patchclass_                (numpatch))
     allocate (patchtype_                 (numpatch))
     allocate (grid_patch_s_ (lon_points,lat_points))
     allocate (grid_patch_e_ (lon_points,lat_points))

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
     allocate (patch2urb_                 (numpatch))
#endif

  END SUBROUTINE allocate_LuLccTimeInvars


  SUBROUTINE SAVE_LuLccTimeInvars

     USE precision
     USE GlobalVars
     USE MOD_TimeInvariants
     USE MOD_PFTimeInvars
     USE MOD_PCTimeInvars
     USE MOD_UrbanTimeInvars
     IMPLICIT NONE

     patchclass_    (:) = patchclass    (:)
     patchtype_     (:) = patchtype     (:)
     grid_patch_s_(:,:) = grid_patch_s(:,:)
     grid_patch_e_(:,:) = grid_patch_e(:,:)

#ifdef PFT_CLASSIFICATION
     pftclass_      (:) = pftclass      (:)
     patch_pft_s_   (:) = patch_pft_s   (:)
     patch_pft_e_   (:) = patch_pft_e   (:)
#endif

#ifdef PC_CLASSIFICATION
     patch2pc_      (:) = patch2pc      (:)
#endif

#ifdef URBAN_MODEL
     urbclass_      (:) = urbclass      (:)
     patch2urb_     (:) = patch2urb     (:)
#endif

  END SUBROUTINE SAVE_LuLccTimeInvars


  SUBROUTINE deallocate_LuLccTimeInvars
! --------------------------------------------------
! Deallocates memory for LuLcc time invariant variables
! --------------------------------------------------

     deallocate (patchclass_   )
     deallocate (patchtype_    )
     deallocate (grid_patch_s_ )
     deallocate (grid_patch_e_ )

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
     deallocate (patch2urb_    )
#endif

  END SUBROUTINE deallocate_LuLccTimeInvars

END MODULE MOD_LuLccTimeInvars
! ---------- EOP ------------
