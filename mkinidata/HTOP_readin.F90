#include <define.h>

SUBROUTINE HTOP_readin (dir_landdata)

! ===========================================================
! Read in the canopy tree top height
! ===========================================================

      USE precision
      USE spmd_task
      USE GlobalVars
      USE LC_Const
      USE PFT_Const
      USE MOD_TimeInvariants
      USE mod_landpatch
#ifdef PFT_CLASSIFICATION
      USE mod_landpft
      USE MOD_PFTimeInvars
      USE MOD_PFTimeVars
#endif
#ifdef PC_CLASSIFICATION
      USE mod_landpc
      USE MOD_PCTimeInvars
      USE MOD_PCTimeVars
#endif
      USE ncio_vector
#ifdef SinglePoint
      USE mod_single_srfdata
#endif

      IMPLICIT NONE

      character(LEN=256), INTENT(in) :: dir_landdata

      ! Local Variables
      character(LEN=256) :: c
      character(LEN=256) :: landdir, lndname
      integer :: i,j,t,p,ps,pe,m,n,npatch

      REAL(r8), allocatable :: htoplc  (:)
      REAL(r8), allocatable :: htoppft (:)

      landdir = trim(dir_landdata) // '/htop'


#ifdef USGS_CLASSIFICATION

      IF (p_is_worker) THEN
         do npatch = 1, numpatch
            m = patchclass(npatch)

            htop(npatch) = htop0(m)
            hbot(npatch) = hbot0(m)

         end do
      ENDIF

#endif

#ifdef IGBP_CLASSIFICATION
#ifdef SinglePoint
      allocate (htoplc (numpatch))
      htoplc(:) = SITE_htop
#else
      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_read_vector (lndname, 'htop_patches', landpatch, htoplc)
#endif

      IF (p_is_worker) THEN
         do npatch = 1, numpatch
            m = patchclass(npatch)

            htop(npatch) = htop0(m)
            hbot(npatch) = hbot0(m)

            ! trees or woody savannas
            IF ( m<6 .or. m==8) THEN
               ! 01/06/2020, yuan: adjust htop reading
               IF (htoplc(npatch) > 2.) THEN
                  htop(npatch) = htoplc(npatch)
                  hbot(npatch) = htoplc(npatch)*hbot0(m)/htop0(m)
                  hbot(npatch) = max(1., hbot(npatch))
                  !htop(npatch) = max(htop(npatch), hbot0(m)*1.2)
               ENDIF
            ENDIF

         end do
      ENDIF

      IF (allocated(htoplc))   deallocate ( htoplc )
#endif


#ifdef PFT_CLASSIFICATION
#ifdef SinglePoint
      allocate(htoppft(numpft))
      htoppft = pack(SITE_htop_pfts, SITE_pctpfts > 0.)
#else
      lndname = trim(landdir)//'/htop_pfts.nc'
      CALL ncio_read_vector (lndname, 'htop_pfts', landpft,   htoppft)
#endif

      IF (p_is_worker) THEN
         do npatch = 1, numpatch
            t = patchtype(npatch)
            m = patchclass(npatch)

            IF (t == 0) THEN
               ps = patch_pft_s(npatch)
               pe = patch_pft_e(npatch)

               DO p = ps, pe
                  n = pftclass(p)

                  htop_p(p) = htop0_p(n)
                  hbot_p(p) = hbot0_p(n)

                  ! for trees
                  ! 01/06/2020, yuan: adjust htop reading
                  IF ( n>0 .and. n<9 .and. htoppft(p)>2.) THEN
                     htop_p(p) = htoppft(p)
                     hbot_p(p) = htoppft(p)*hbot0_p(n)/htop0_p(n)
                     hbot_p(p) = max(1., hbot_p(p)) !weinan
                  ENDIF
               ENDDO

               htop(npatch) = sum(htop_p(ps:pe)*pftfrac(ps:pe))
               hbot(npatch) = sum(hbot_p(ps:pe)*pftfrac(ps:pe))

            ELSE
               htop(npatch) = htop0(m)
               hbot(npatch) = hbot0(m)
            ENDIF

         ENDDO
      ENDIF

      IF (allocated(htoppft)) deallocate(htoppft)
#endif

#ifdef PC_CLASSIFICATION
#ifdef SinglePoint
      allocate(htoplc(1))
      htoplc(:) = sum(SITE_htop_pfts * SITE_pctpfts)
#else
      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_read_vector (lndname, 'htop_patches', landpatch, htoplc )
#endif

      IF (p_is_worker) THEN
         do npatch = 1, numpatch
            t = patchtype(npatch)
            m = patchclass(npatch)
            IF (t == 0) THEN
               p = patch2pc(npatch)
               htop_c(:,p) = htop0_p(:)
               hbot_c(:,p) = hbot0_p(:)

               DO n = 1, N_PFT-1
                  ! 01/06/2020, yuan: adjust htop reading
                  IF (n < 9 .and. htoplc(npatch)>2.) THEN
                     htop_c(n,p) = htoplc(npatch)
                  ENDIF
               ENDDO
               htop(npatch) = sum(htop_c(:,p)*pcfrac(:,p))
               hbot(npatch) = sum(hbot_c(:,p)*pcfrac(:,p))
            ELSE
               htop(npatch) = htop0(m)
               hbot(npatch) = hbot0(m)
            ENDIF
         end do
      ENDIF

      IF (allocated(htoplc)) deallocate(htoplc)
#endif

END SUBROUTINE HTOP_readin
