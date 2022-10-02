#include <define.h>

SUBROUTINE CROP_readin (dir_landdata)
   ! ===========================================================
   ! Read in the LAI, the LAI dataset was created by Yuan et al. (2011)
   ! http://globalchange.bnu.edu.cn
   !
   ! Created by Yongjiu Dai, March, 2014
   ! ===========================================================

   use precision
   use mod_namelist
   use spmd_task
   use ncio_vector
   use mod_landpatch
   use MOD_TimeInvariants
   use MOD_TimeVariables
#ifdef CLMDEBUG 
   use mod_colm_debug
#endif
   
   USE GlobalVars
   USE LC_Const
#ifdef PFT_CLASSIFICATION
   USE mod_landpft
   USE MOD_PFTimeVars
#endif
#ifdef PC_CLASSIFICATION
   USE mod_landpc
   USE MOD_PCTimeVars
#endif
#ifdef SinglePoint
   USE mod_single_srfdata
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
!            m = patchclass(npatch)
!            IF( m == 12 )then
               pdrice2    (npatch) = int(pdrice2_tmp    (npatch))
!            ELSE
!               pdcorn     (npatch) = -9999
!               pdswheat   (npatch) = -9999
!               pdwwheat   (npatch) = -9999
!               pdsoybean  (npatch) = -9999
!               pdcotton   (npatch) = -9999
!               pdrice1    (npatch) = -9999
!               pdrice2    (npatch) = -9999
!               pdsugarcane(npatch) = -9999
!            ENDIF
!            IF (pdcorn_tmp(npatch)      < 1E-10) THEN
!               pdcorn(npatch)      = -9999
!            ENDIF
!            IF (pdswheat_tmp(npatch)    < 1E-10) THEN
!               pdswheat(npatch)    = -9999
!            ENDIF
!            IF (pdwwheat_tmp(npatch)    < 1E-10) THEN
!               pdwwheat(npatch)    = -9999
!            ENDIF
!            IF (pdsoybean_tmp(npatch)   < 1E-10) THEN
!               pdsoybean(npatch)   = -9999
!            ENDIF
!            IF (pdcotton_tmp(npatch)    < 1E-10) THEN
!               pdcotton(npatch)    = -9999
!            ENDIF
!            IF (pdrice1_tmp(npatch)     < 1E-10) THEN
!               pdrice1(npatch)     = -9999
!            ENDIF
!            IF (pdrice2_tmp(npatch)     < 1E-10) THEN
!               pdrice2(npatch)     = -9999
!            ENDIF
!            IF (pdsugarcane_tmp(npatch) < 1E-10) THEN
!               pdsugarcane(npatch) = -9999
!            ENDIF
         ENDDO

      ENDIF
   ENDIF 

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

   !print*,'landdir cropreadin',trim(landdir)
   lndname = trim(landdir) // '/fertnitro_pfts.nc'
   !print*,'in cropreadin',lndname
   call ncio_read_vector (lndname, 'fertnitro_pfts',  landpft, fertnitro_tmp)
!   print*,'fertnitro_tmp',p_iam_glb,fertnitro_tmp
   if (p_is_worker) then
      if (numpatch > 0) then
         do npatch = 1, numpatch
            m = patchclass(npatch)
   !         print*,'read fertnitro_pfts',npatch,p_iam_glb,m,patch_pft_s(npatch), patch_pft_e(npatch)
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
!            if(m .eq. 1 .or. m .eq. 12)print*,'fertnitro_p',p_iam_glb,npatch,patch_pft_s(npatch),patch_pft_e(npatch),fertnitro_p(patch_pft_s(npatch):patch_pft_e(npatch))
!            if(m .ne. 1 .and. m .ne. 12 .and. patch_pft_s(npatch) .ne. -1 .and. patch_pft_e(npatch) .ne. -1)print*,'missing pft'
         end do
      ENDIF
   ENDIF 
   deallocate(pdrice2_tmp     )
   deallocate(plantdate_tmp   )
   deallocate(fertnitro_tmp   )

#endif
END SUBROUTINE CROP_readin
