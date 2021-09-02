#include <define.h>

SUBROUTINE LAI_readin_nc (lon_points,lat_points,&
           month, dir_model_landdata)

! ===========================================================
! Read in the LAI, the LAI dataset was created by Yuan et al. (2019)
! http://globalchange.bnu.edu.cn
! ===========================================================

      use precision
      USE GlobalVars
      USE LC_Const
      USE PFT_Const
      use MOD_TimeInvariants
      use MOD_TimeVariables
      USE MOD_PFTimeVars
      USE MOD_PCTimeVars
      USE MOD_PFTimeInvars
      USE MOD_PCTimeInvars
      USE ncio
      use omp_lib

      IMPLICIT NONE

      integer, INTENT(in) :: lon_points
      integer, INTENT(in) :: lat_points
      integer, INTENT(in) :: month
      character(LEN=256), INTENT(in) :: dir_model_landdata

      character(LEN=256) :: c
      character(LEN=256) :: lndname
      integer :: ncid
      INTEGER :: lclai_vid, lcsai_vid, pftlai_vid, pftsai_vid
      INTEGER :: pclai_vid, pcsai_vid, pctpc_vid
      integer :: i, j, m, n, t, p, ps, pe, ep, npatch

      REAL(r8), allocatable :: lclai(:,:,:)
      REAL(r8), allocatable :: lcsai(:,:,:)
      REAL(r8), allocatable :: pftlai(:,:,:)
      REAL(r8), allocatable :: pftsai(:,:,:)
      REAL(r8), allocatable :: pclai(:,:,:,:)
      REAL(r8), allocatable :: pcsai(:,:,:,:)
      REAL(r8), allocatable :: pctpc(:,:,:,:)

! READ in Leaf area index and stem area index

      lndname = trim(dir_model_landdata)//'global_0.5x0.5.MOD2005_V4.5.nc'
      print*,trim(lndname)
      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )

#ifdef IGBP_CLASSIFICATION
      allocate ( lclai(1:lon_points,1:lat_points,1:N_land_classification) )
      allocate ( lcsai(1:lon_points,1:lat_points,1:N_land_classification) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_LC_LAI", lclai_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_LC_SAI", lcsai_vid ) )
      CALL nccheck( nf90_get_var(ncid, lclai_vid, lclai, &
           start=(/1,1,1,month/), &
           count=(/lon_points,lat_points,N_land_classification,1/)) )
      CALL nccheck( nf90_get_var(ncid, lcsai_vid, lcsai, &
           start=(/1,1,1,month/), &
           count=(/lon_points,lat_points,N_land_classification,1/)) )

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,m)
#endif
      do npatch = 1, numpatch
         i = patch2lon(npatch)
         j = patch2lat(npatch)
         m = patchclass(npatch)
         if( m == 0 )then
             fveg(npatch)  = 0.
             tlai(npatch)  = 0.
             tsai(npatch)  = 0.
             green(npatch) = 0.
         else
             fveg(npatch)  = fveg0(m)           !fraction of veg. cover
             IF (fveg0(m) > 0) THEN
                tlai(npatch)  = lclai(i,j,m)/fveg0(m) !leaf area index
                tsai(npatch)  = lcsai(i,j,m)/fveg0(m) !stem are index
                green(npatch) = 1.                    !fraction of green leaf
             ELSE 
                tlai(npatch)  = 0.       !leaf area index
                tsai(npatch)  = 0.       !stem are index
                green(npatch) = 0.       !fraction of green leaf
             ENDIF 
         endif
      end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( lclai )
      deallocate ( lcsai )

#endif

#ifdef PFT_CLASSIFICATION
      allocate ( pftlai(1:lon_points,1:lat_points,0:N_PFT-1) )
      allocate ( pftsai(1:lon_points,1:lat_points,0:N_PFT-1) )
      allocate ( pclai(1:lon_points,1:lat_points, &
                         0:N_PFT-1,1:N_land_classification) )
      allocate ( pcsai(1:lon_points,1:lat_points, &
                         0:N_PFT-1,1:N_land_classification) )
      allocate ( pctpc(1:lon_points,1:lat_points, &
                         0:N_PFT-1,1:N_land_classification) )

      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_PFT_LAI",  pftlai_vid  ) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_PFT_SAI",  pftsai_vid  ) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_ePFT_LAI",   pclai_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_ePFT_SAI",   pcsai_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_ePFT",           pctpc_vid ) )

      CALL nccheck( nf90_get_var(ncid, pftlai_vid, pftlai, &
                    start=(/1,1,1,month/), &
                    count=(/lon_points,lat_points,N_PFT,1/)) )
      CALL nccheck( nf90_get_var(ncid, pftsai_vid, pftsai, &
                    start=(/1,1,1,month/), & 
                    count=(/lon_points,lat_points,N_PFT,1/)) )
      CALL nccheck( nf90_get_var(ncid, pclai_vid, pclai, &
                    start=(/1,1,1,1,month/), &
                    count=(/lon_points,lat_points,N_PFT,N_land_classification,1/)) )
      CALL nccheck( nf90_get_var(ncid, pcsai_vid, pcsai, &
                    start=(/1,1,1,1,month/), &
                    count=(/lon_points,lat_points,N_PFT,N_land_classification,1/)) )
      CALL nccheck( nf90_get_var(ncid, pctpc_vid, pctpc ) )

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,m,n,t,p,ps,pe)
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
               tlai_p(p) = pftlai(i,j,n)
               tsai_p(p) = pftsai(i,j,n)
            ENDDO

            tlai(npatch) = sum(tlai_p(ps:pe) * pftfrac(ps:pe))
            tsai(npatch) = sum(tsai_p(ps:pe) * pftfrac(ps:pe))

         ELSE
! 12/28/2019, yuan: Bug
            ! pctpc from 1-100% -> 0-1
            tlai(npatch)  = sum(pclai(i,j,:,m) * pctpc(i,j,:,m)/100.) 
            tsai(npatch)  = sum(pcsai(i,j,:,m) * pctpc(i,j,:,m)/100.) 
         ENDIF
         
         green(npatch) = 1.                
         fveg(npatch)  = fveg0(m)

      end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( pftlai  )
      deallocate ( pftsai  )
      deallocate ( pclai   )
      deallocate ( pcsai   )
      deallocate ( pctpc   )

#endif

#ifdef PC_CLASSIFICATION
      allocate ( pclai(1:lon_points,1:lat_points, &
                          0:N_PFT-1,1:N_land_classification) )
      allocate ( pcsai(1:lon_points,1:lat_points, &
                          0:N_PFT-1,1:N_land_classification) )
      allocate ( pctpc(1:lon_points,1:lat_points, &
                          0:N_PFT-1,1:N_land_classification) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_ePFT_LAI", pclai_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_ePFT_SAI", pcsai_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_ePFT",         pctpc_vid ) )
      
      CALL nccheck( nf90_get_var(ncid, pclai_vid, pclai, &
                    start=(/1,1,1,1,month/), &
                    count=(/lon_points,lat_points,N_PFT,N_land_classification,1/)) )
      CALL nccheck( nf90_get_var(ncid, pcsai_vid, pcsai, &
                    start=(/1,1,1,1,month/), &
                    count=(/lon_points,lat_points,N_PFT,N_land_classification,1/)) )
      CALL nccheck( nf90_get_var(ncid, pctpc_vid, pctpc ) )

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,t,m,ep)
#endif
      do npatch = 1, numpatch
         i = patch2lon(npatch)
         j = patch2lat(npatch)
         t = patchtype(npatch)
         m = patchclass(npatch)
         IF (t == 0) THEN
            ep = patch2pc(npatch)
            tlai_c(:,ep) = pclai(i,j,:,m)
            tsai_c(:,ep) = pcsai(i,j,:,m)
         ENDIF 

! 12/28/2019, yuan: Bug
         ! pctpc from 1-100% -> 0-1
         tlai(npatch)  = sum(pclai(i,j,:,m)*pctpc(i,j,:,m)/100.) 
         tsai(npatch)  = sum(pcsai(i,j,:,m)*pctpc(i,j,:,m)/100.) 
         fveg(npatch)  = fveg0(m)
         green(npatch) = 1.                 

      end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( pclai )
      deallocate ( pcsai )
      deallocate ( pctpc )

#endif
      
      CALL nccheck( nf90_close(ncid) )

END SUBROUTINE LAI_readin_nc
