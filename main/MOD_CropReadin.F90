#include <define.h>

MODULE MOD_CropReadin

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: CROP_readin


!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------


   SUBROUTINE CROP_readin (dir_landdata)
      ! ===========================================================
      ! ! DESCRIPTION:
      ! Read in crop planting date from data, and fertilization from data.
      ! Save these data in patch vector.
      !
      ! Original: Shupeng Zhang, Zhongwang Wei, and Xingjie Lu, 2022
      ! ===========================================================

      use MOD_Precision
      use MOD_Namelist
      use MOD_SPMD_Task
      use MOD_NetCDFVector
      use MOD_LandPatch
      use MOD_Vars_TimeInvariants
      use MOD_Vars_TimeVariables
#ifdef CoLMDEBUG
      use MOD_CoLMDebug
#endif

      USE MOD_Vars_Global
      USE MOD_Const_LC
#ifdef LULC_IGBP_PFT
      USE MOD_LandPFT
      USE MOD_Vars_PFTimeVariables
#endif
#ifdef LULC_IGBP_PC
      USE MOD_LandPC
      USE MOD_Vars_PCTimeVariables
#endif
#ifdef SinglePoint
      USE MOD_SingleSrfdata
#endif

      IMPLICIT NONE

      character(LEN=*), INTENT(in) :: dir_landdata
      real(r8),allocatable :: pdrice2_tmp   (:)
      real(r8),allocatable :: plantdate_tmp (:)
      real(r8),allocatable :: fertnitro_tmp (:)

      ! Local variables
      integer :: iyear, itime
      character(LEN=256) :: cyear, ctime
      character(LEN=256) :: landdir, lndname
      integer :: m, npatch,ns, ipft
      CHARACTER(LEN=4) ::cx
      ! READ in crops

#ifdef CROP
      allocate(pdrice2_tmp     (numpatch))
      allocate(plantdate_tmp   (numpft))
      allocate(fertnitro_tmp   (numpft))

      landdir = trim(dir_landdata) // '/crop'
      !print*,'landdir cropreadin',trim(landdir)
      lndname = trim(landdir) // '/plantdate_patches.nc'
      call ncio_read_vector (lndname, 'plantdate_rice2_patches'    , landpatch, pdrice2_tmp)
      IF (p_is_worker) then
         IF (numpatch > 0) then
            DO npatch = 1, numpatch
               pdrice2    (npatch) = int(pdrice2_tmp    (npatch))
            ENDDO

         ENDIF
      ENDIF

      lndname = trim(landdir) // '/plantdate_pfts.nc'
      call ncio_read_vector (lndname, 'plantdate_pfts',  landpft, plantdate_tmp)
      if (p_is_worker) then
         if (numpatch > 0) then
            do npatch = 1, numpatch
               m = patchclass(npatch)
               if( m == 12 )then
                   do ipft = patch_pft_s(npatch), patch_pft_e(npatch)
                      plantdate_p(ipft) = plantdate_tmp(ipft)
                      if(plantdate_p(ipft) <= 0._r8) then
                         plantdate_p(ipft) = -99999999._r8
                      end if
                   end do
               else
                   if(m == 1)then
                      do ipft = patch_pft_s(npatch), patch_pft_e(npatch)
                         plantdate_p(ipft) = -99999999._r8
                      end do
                   end if
               end if
            end do
         end if
      end if

      lndname = trim(landdir) // '/fertnitro_pfts.nc'
      call ncio_read_vector (lndname, 'fertnitro_pfts',  landpft, fertnitro_tmp)
      if (p_is_worker) then
         if (numpatch > 0) then
            do npatch = 1, numpatch
               m = patchclass(npatch)
               if( m == 12 )then
                  do ipft = patch_pft_s(npatch), patch_pft_e(npatch)
                     fertnitro_p(ipft)  = fertnitro_tmp(ipft)
                     if (fertnitro_p(ipft) < 0._r8) then
                        fertnitro_p(ipft)  = -99999999._r8
                     endif
                  end do
               else
                  if(m == 1)then
                     do ipft = patch_pft_s(npatch), patch_pft_e(npatch)
                        fertnitro_p(ipft)  = -99999999._r8
                        if (fertnitro_p(ipft) < 0._r8) then
                           fertnitro_p(ipft)  = -99999999._r8
                        endif
                     end do
                  end if
               endif
            end do
         ENDIF
      ENDIF
      deallocate(pdrice2_tmp     )
      deallocate(plantdate_tmp   )
      deallocate(fertnitro_tmp   )

#endif
   END SUBROUTINE CROP_readin

END MODULE MOD_CropReadin
