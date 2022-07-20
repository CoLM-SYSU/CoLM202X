#include <define.h>

SUBROUTINE LAI_readin (time, dir_landdata)
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

   IMPLICIT NONE

   integer, INTENT(in) :: time 
   character(LEN=256), INTENT(in) :: dir_landdata

   ! Local variables
   character(LEN=256) :: c
   character(LEN=256) :: lndname
   integer :: m, npatch

#ifdef USGS_CLASSIFICATION
   real(r8), dimension(24), parameter :: &   ! Maximum fractional cover of vegetation [-]
      vegc=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, &
      1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0 /)
#endif

   ! READ in Leaf area index and stem area index



#ifdef USGS_CLASSIFICATION
   
   write(c,'(i3.3)') time
   lndname = trim(dir_landdata)//'/LAI_patches'//trim(c)//'.nc'
   call ncio_read_vector (lndname, 'LAI_patches',  landpatch, tlai)

   if (p_is_worker) then
      if (numpatch > 0) then

         do npatch = 1, numpatch
            m = patchclass(npatch)
            if( m == 0 )then
               fveg(npatch)  = 0.
               tlai(npatch)  = 0.
               tsai(npatch)  = 0.
               green(npatch) = 0.
            else
               fveg(npatch)  = vegc(m)    !fraction of veg. cover
               IF (vegc(m) > 0) THEN 
                  tlai(npatch)  = tlai(npatch)/vegc(m) !leaf area index
                  tsai(npatch)  = sai0(m) !stem are index
                  green(npatch) = 1.      !fraction of green leaf
               ELSE 
                  tlai(npatch)  = 0.  
                  tsai(npatch)  = 0.   
                  green(npatch) = 0.    
               ENDIF 
            ENDIF
         end do

      ENDIF
   ENDIF 
         
#endif

#ifdef IGBP_CLASSIFICATION

   write(c,'(i2.2)') time
   lndname = trim(dir_landdata)//'/LAI_patches'//trim(c)//'.nc'
   call ncio_read_vector (lndname, 'LAI_patches',  landpatch, tlai )

   write(c,'(i2.2)') time
   lndname = trim(dir_landdata)//'/SAI_patches'//trim(c)//'.nc'
   call ncio_read_vector (lndname, 'SAI_patches',  landpatch, tsai )

   if (p_is_worker) then
      if (numpatch > 0) then

         do npatch = 1, numpatch
            m = patchclass(npatch)
            if( m == 0 )then
               fveg(npatch)  = 0.
               tlai(npatch)  = 0.
               tsai(npatch)  = 0.
               green(npatch) = 0.
            else
               fveg(npatch)  = fveg0(m)           !fraction of veg. cover
               IF (fveg0(m) > 0) THEN
                  tlai(npatch)  = tlai(npatch)/fveg0(m) !leaf area index
                  tsai(npatch)  = tsai(npatch)/fveg0(m) !stem are index
                  green(npatch) = 1.                    !fraction of green leaf
               ELSE 
                  tlai(npatch)  = 0.       !leaf area index
                  tsai(npatch)  = 0.       !stem are index
                  green(npatch) = 0.       !fraction of green leaf
               ENDIF 
            endif
         end do

      ENDIF
   ENDIF 

#endif

#ifdef PFT_CLASSIFICATION

   write(c,'(i2.2)') time
#ifndef LAIfdbk
   lndname = trim(dir_landdata)//'/LAI_patches'//trim(c)//'.nc'
   call ncio_read_vector (lndname, 'LAI_patches',  landpatch, tlai )
#endif
   
   lndname = trim(dir_landdata)//'/SAI_patches'//trim(c)//'.nc'
   call ncio_read_vector (lndname, 'SAI_patches',  landpatch, tsai )
   
#ifndef LAIfdbk
   lndname = trim(dir_landdata)//'/LAI_pfts'//trim(c)//'.nc'
   call ncio_read_vector (lndname, 'LAI_pfts', landpft, tlai_p )
#endif
   
   lndname = trim(dir_landdata)//'/SAI_pfts'//trim(c)//'.nc'
   call ncio_read_vector (lndname, 'SAI_pfts', landpft, tsai_p )

   if (p_is_worker) then
      if (numpatch > 0) then
         do npatch = 1, numpatch
            m = patchclass(npatch)

            green(npatch) = 1.                
            fveg (npatch)  = fveg0(m)

         end do
      ENDIF
   ENDIF 

#endif

#ifdef PC_CLASSIFICATION

   write(c,'(i2.2)') time
   lndname = trim(dir_landdata)//'/LAI_patches'//trim(c)//'.nc'
   call ncio_read_vector (lndname, 'LAI_patches',  landpatch, tlai )
   
   lndname = trim(dir_landdata)//'/SAI_patches'//trim(c)//'.nc'
   call ncio_read_vector (lndname, 'SAI_patches',  landpatch, tsai )

   lndname = trim(dir_landdata)//'/LAI_pcs'//trim(c)//'.nc'
   call ncio_read_vector (lndname, 'LAI_pcs', N_PFT, landpc, tlai_c )

   lndname = trim(dir_landdata)//'/SAI_pcs'//trim(c)//'.nc'
   call ncio_read_vector (lndname, 'SAI_pcs', N_PFT, landpc, tsai_c )

   if (p_is_worker) then
      if (numpatch > 0) then
         do npatch = 1, numpatch
            m = patchclass(npatch)
            fveg (npatch)  = fveg0(m)
            green(npatch) = 1.                 
         end do
      ENDIF 
   ENDIF 

#endif

END SUBROUTINE LAI_readin
