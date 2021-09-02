#include <define.h>

SUBROUTINE HTOP_readin_nc (lon_points,lat_points,dir_model_landdata)

! ===========================================================
! Read in the canopy tree top height
! ===========================================================

      USE precision
      USE GlobalVars
      USE LC_Const
      USE PFT_Const
      USE MOD_TimeInvariants
      USE MOD_PFTimeInvars
      USE MOD_PCTimeInvars
      USE MOD_PFTimeVars
      USE MOD_PCTimeVars
      USE ncio
      USE omp_lib

      IMPLICIT NONE

      integer, INTENT(in) :: lon_points
      integer, INTENT(in) :: lat_points
      character(LEN=256), INTENT(in) :: dir_model_landdata

      character(LEN=256) :: c
      character(LEN=256) :: lndname
      integer :: ncid
      INTEGER :: htoplc_vid, htoppft_vid
      integer :: i,j,t,p,ps,pe,m,n,npatch

      REAL(r8), allocatable :: htoplc(:,:,:)
      REAL(r8), allocatable :: htoppft(:,:,:)

      lndname = trim(dir_model_landdata)//'global_0.5x0.5.MOD2005_V4.5.nc'
      print*,trim(lndname)
      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )

#ifdef IGBP_CLASSIFICATION
      allocate ( htoplc(1:lon_points,1:lat_points,1:N_land_classification) )
      CALL nccheck( nf90_inq_varid(ncid, "HTOP_LC", htoplc_vid ) )
      CALL nccheck( nf90_get_var(ncid, htoplc_vid, htoplc ) )

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,m)
#endif
      do npatch = 1, numpatch
         i = patch2lon(npatch)
         j = patch2lat(npatch)
         m = patchclass(npatch)

         htop(npatch) = htop0(m)
         hbot(npatch) = hbot0(m)
         
         ! trees or woody savannas
         IF ( m<6 .or. m==8) THEN
! yuan, 01/06/2020: adjust htop reading
            IF (htoplc(i,j,m) > 2.) THEN
               htop(npatch) = htoplc(i,j,m)
               hbot(npatch) = htoplc(i,j,m)*hbot0(m)/htop0(m)
               hbot(npatch) = max(1., hbot(npatch))
               !htop(npatch) = max(htop(npatch), hbot0(m)*1.2)
            ENDIF
         ENDIF
         
      end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
      deallocate ( htoplc )
#endif

#ifdef PFT_CLASSIFICATION
      allocate ( htoppft(1:lon_points,1:lat_points,0:N_PFT-1) )
      CALL nccheck( nf90_inq_varid(ncid, "HTOP_PFT", htoppft_vid ) )
      CALL nccheck( nf90_get_var(ncid, htoppft_vid, htoppft ) )

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,t,p,ps,pe,m,n)
#endif
      do npatch = 1, numpatch
         i = patch2lon(npatch)
         j = patch2lat(npatch)
         t = patchtype(npatch)
         m = patchclass(npatch)
         
         IF (t == 0) THEN
            ps = patch_pft_s(npatch)
            pe = patch_pft_e(npatch)

            DO p = ps, pe
               n = pftclass(p)
               
               htop_p(p) = htop0_p(n)
               hbot_p(p) = hbot0_p(n)
               
               ! trees
! yuan, 01/06/2020: adjust htop reading
               IF ( n>0 .and. n<9 .and. htoppft(i,j,n)>2.) THEN
                  htop_p(p) = htoppft(i,j,n)
                  hbot_p(p) = htoppft(i,j,n)*hbot0_p(n)/htop0_p(n)
                  hbot_p(p) = max(1., hbot_p(n))
               ENDIF
            ENDDO

            htop(npatch) = sum(htop_p(ps:pe)*pftfrac(ps:pe))
            hbot(npatch) = sum(hbot_p(ps:pe)*pftfrac(ps:pe))

         ELSE 
            htop(npatch) = htop0(m)
            hbot(npatch) = hbot0(m)
         ENDIF 
 
      ENDDO 
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
      deallocate ( htoppft  )
#endif

#ifdef PC_CLASSIFICATION
      allocate ( htoplc(1:lon_points,1:lat_points,1:N_land_classification) )
      CALL nccheck( nf90_inq_varid(ncid, "HTOP_LC", htoplc_vid ) )
      CALL nccheck( nf90_get_var(ncid, htoplc_vid, htoplc ) )

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,t,p,m,n)
#endif
      do npatch = 1, numpatch
         i = patch2lon(npatch)
         j = patch2lat(npatch)
         t = patchtype(npatch)
         m = patchclass(npatch)
         IF (t == 0) THEN
            p = patch2pc(npatch)
            htop_c(:,p) = htop0_p(:)
            hbot_c(:,p) = hbot0_p(:)

            DO n = 1, N_PFT-1
! yuan, 01/06/2020: adjust htop reading
               IF (n < 9 .and. htoplc(i,j,m)>2.) THEN
                  htop_c(n,p) = htoplc(i,j,m)
               ENDIF
            ENDDO
            htop(npatch) = sum(htop_c(:,p)*pcfrac(:,p))
            hbot(npatch) = sum(hbot_c(:,p)*pcfrac(:,p))
         ELSE
            htop(npatch) = htop0(m)
            hbot(npatch) = hbot0(m)
         ENDIF
      end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
      deallocate ( htoplc )
#endif

      CALL nccheck( nf90_close(ncid) )

END SUBROUTINE HTOP_readin_nc
