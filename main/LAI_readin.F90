#include <define.h>

SUBROUTINE LAI_readin (year, time, dir_landdata)
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

   integer, INTENT(in) :: year, time
   character(LEN=256), INTENT(in) :: dir_landdata

   ! Local variables
   integer :: iyear, itime
   character(LEN=256) :: cyear, ctime
   character(LEN=256) :: landdir, lndname
   integer :: m, npatch

#ifdef USGS_CLASSIFICATION
   real(r8), dimension(24), parameter :: &   ! Maximum fractional cover of vegetation [-]
      vegc=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, &
      1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0 /)
#endif

   ! READ in Leaf area index and stem area index

   landdir = trim(dir_landdata) // '/LAI'

#ifdef SinglePoint
   IF (.not. DEF_LAI_CLIM) THEN
      iyear = findloc(SITE_LAI_year, year, dim=1)
      itime = (time-1)/8 + 1
   ENDIF
#endif

#if (defined USGS_CLASSIFICATION || defined IGBP_CLASSIFICATION)

#ifdef SinglePoint
   IF (DEF_LAI_CLIM) THEN
      tlai(:) = SITE_LAI_clim(time)
      tsai(:) = SITE_SAI_clim(time)
   ELSE
      tlai(:) = SITE_LAI_modis(itime,iyear)
   ENDIF
#else
   IF (DEF_LAI_CLIM) THEN
      write(cyear,'(i4.4)') year
      write(ctime,'(i2.2)') time

      lndname = trim(landdir)//'/'//trim(cyear)//'/LAI_patches'//trim(ctime)//'.nc'
      call ncio_read_vector (lndname, 'LAI_patches',  landpatch, tlai)

      lndname = trim(landdir)//'/'//trim(cyear)//'/SAI_patches'//trim(ctime)//'.nc'
      call ncio_read_vector (lndname, 'SAI_patches',  landpatch, tsai)
   ELSE
      write(cyear,'(i4.4)') year
      write(ctime,'(i3.3)') time
      lndname = trim(landdir)// '/' // trim(cyear) //'/LAI_patches'//trim(ctime)//'.nc'
      call ncio_read_vector (lndname, 'LAI_patches',  landpatch, tlai)
   ENDIF
#endif

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
                  IF (DEF_LAI_CLIM) THEN
                     tsai(npatch)  = tsai(npatch)/fveg0(m) !stem are index
                  ELSE
                     tsai(npatch)  = sai0(m) !stem are index
                  ENDIF
                  green(npatch) = 1.      !fraction of green leaf
               ELSE
                  tlai(npatch)  = 0.
                  tsai(npatch)  = 0.
                  green(npatch) = 0.
               ENDIF
            endif
         end do

      ENDIF
   ENDIF

#endif

#ifdef PFT_CLASSIFICATION

#ifdef SinglePoint
   !TODO: how to add time parameter in single point case
   IF (DEF_LAI_CLIM) THEN
      tlai_p(:) = pack(SITE_LAI_pfts_clim(:,time), SITE_pctpfts > 0.)
      tsai_p(:) = pack(SITE_SAI_pfts_clim(:,time), SITE_pctpfts > 0.)
      tlai(:)   = sum (SITE_LAI_pfts_clim(:,time) * SITE_pctpfts)
      tsai(:)   = sum (SITE_SAI_pfts_clim(:,time) * SITE_pctpfts)
   ENDIF
#else

   write(cyear,'(i4.4)') year
   write(ctime,'(i2.2)') time
#ifndef LAIfdbk
   lndname = trim(landdir)//'/'//trim(cyear)//'/LAI_patches'//trim(ctime)//'.nc'
   call ncio_read_vector (lndname, 'LAI_patches',  landpatch, tlai )
#endif
   lndname = trim(landdir)//'/'//trim(cyear)//'/SAI_patches'//trim(ctime)//'.nc'
   call ncio_read_vector (lndname, 'SAI_patches',  landpatch, tsai )
#ifndef LAIfdbk
   lndname = trim(landdir)//'/'//trim(cyear)//'/LAI_pfts'//trim(ctime)//'.nc'
   call ncio_read_vector (lndname, 'LAI_pfts', landpft, tlai_p )
#endif
   lndname = trim(landdir)//'/'//trim(cyear)//'/SAI_pfts'//trim(ctime)//'.nc'
   call ncio_read_vector (lndname, 'SAI_pfts', landpft, tsai_p )

#endif

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

#ifdef SinglePoint
   IF (DEF_LAI_CLIM) THEN
      tlai(:)   = sum(SITE_LAI_pfts_clim(:,time) * SITE_pctpfts)
      tsai(:)   = sum(SITE_SAI_pfts_clim(:,time) * SITE_pctpfts)
      tlai_c(:,1) = SITE_LAI_pfts_clim(:,time)
      tsai_c(:,1) = SITE_SAI_pfts_clim(:,time)
   ENDIF
#else

   write(cyear,'(i4.4)') year
   write(ctime,'(i2.2)') time
   lndname = trim(landdir)//'/'//trim(cyear)//'/LAI_patches'//trim(ctime)//'.nc'
   call ncio_read_vector (lndname, 'LAI_patches',  landpatch, tlai )

   lndname = trim(landdir)//'/'//trim(cyear)//'/SAI_patches'//trim(ctime)//'.nc'
   call ncio_read_vector (lndname, 'SAI_patches',  landpatch, tsai )

   lndname = trim(landdir)//'/'//trim(cyear)//'/LAI_pcs'//trim(ctime)//'.nc'
   call ncio_read_vector (lndname, 'LAI_pcs', N_PFT, landpc, tlai_c )

   lndname = trim(landdir)//'/'//trim(cyear)//'/SAI_pcs'//trim(ctime)//'.nc'
   call ncio_read_vector (lndname, 'SAI_pcs', N_PFT, landpc, tsai_c )

#endif

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
