#include <define.h>

SUBROUTINE initialize (casename,dir_model_landdata,dir_restart_hist,&
                       idate,greenwich,lon_points,lat_points,&
                       tg_xy,albvb_xy,albvd_xy,albnb_xy,albnd_xy,&
                       trad_xy,rib_xy,fm_xy,fh_xy,fq_xy)

! ======================================================================
! initialization routine for land surface model.
!
! Created by Yongjiu Dai, 09/15/1999
! Revised by Yongjiu Dai, 08/30/2002
! Revised by Yongjiu Dai, 03/2014
!             
! ======================================================================
   USE precision
   USE GlobalVars
   USE PhysicalConstants
   USE MOD_TimeInvariants
   USE MOD_TimeVariables
   USE MOD_PFTimeInvars
   USE MOD_PCTimeInvars
   USE MOD_PFTimeVars
   USE MOD_PCTimeVars
   USE LC_Const
   USE PFT_Const
   USE timemanager
   USE ncio
   USE netcdf
   USE omp_lib

   IMPLICIT NONE

! ----------------------------------------------------------------------
   CHARACTER(LEN=256), INTENT(in) :: casename           !casename name
   CHARACTER(LEN=256), INTENT(in) :: dir_model_landdata !
   CHARACTER(LEN=256), INTENT(in) :: dir_restart_hist   !
   LOGICAL, INTENT(in)    :: greenwich   !true: greenwich time, false: local time
   INTEGER, INTENT(in)    :: lon_points  !number of longitude points on model grid
   INTEGER, INTENT(in)    :: lat_points  !number of latitude points on model grid
   INTEGER, INTENT(inout) :: idate(3)    !year, julian day, seconds of the starting time

! required by atmospheric models initialization (such as GRAPES, RSM, ...)
   REAL(r8), intent(out) :: tg_xy   (lon_points,lat_points) ! 
   REAL(r8), intent(out) :: albvb_xy(lon_points,lat_points) ! 
   REAL(r8), intent(out) :: albvd_xy(lon_points,lat_points) ! 
   REAL(r8), intent(out) :: albnb_xy(lon_points,lat_points) ! 
   REAL(r8), intent(out) :: albnd_xy(lon_points,lat_points) ! 
   REAL(r8), intent(out) :: trad_xy (lon_points,lat_points) ! 
   REAL(r8), intent(out) :: rib_xy  (lon_points,lat_points) ! 
   REAL(r8), intent(out) :: fm_xy   (lon_points,lat_points) ! 
   REAL(r8), intent(out) :: fh_xy   (lon_points,lat_points) ! 
   REAL(r8), intent(out) :: fq_xy   (lon_points,lat_points) ! 

! ------------------------ local variables -----------------------------
! surface classification and soil information

  REAL(r8) latixy  (lon_points,lat_points)          !latitude in radians
  REAL(r8) longxy  (lon_points,lat_points)          !longitude in radians
  REAL(r8) latdeg  (lat_points)                     !latitude in degree
  REAL(r8) londeg  (lon_points)                     !longitude in degree
  REAL(r8) area_gridcells (lon_points,lat_points)   !area of gridcells (km^2)

  REAL(r8), allocatable :: landfrac(:,:)            ! 
  REAL(r8), allocatable :: pctlc(:,:,:)             ! 
  REAL(r8), allocatable :: pctelc(:,:,:)            ! 
  REAL(r8), allocatable :: pctpft(:,:,:)            ! 
  REAL(r8), allocatable :: pcturban(:,:)            ! 
  REAL(r8), allocatable :: pctwater(:,:)            ! 
  REAL(r8), allocatable :: pctwetland(:,:)          ! 
  REAL(r8), allocatable :: pctglacier(:,:)          ! 
  REAL(r8), allocatable :: pctpc(:,:,:,:)           ! 
  REAL(r8), allocatable :: fraction_patches(:,:,:)  !fraction of the patch of landtypes in gridcells

#if(defined SOILINI)
  INTEGER :: lusoil
  INTEGER :: nl_soil_ini
  REAL(r8), allocatable :: snow_d_grid(:,:)
  REAL(r8), allocatable :: snow_d(:)

  REAL(r8), allocatable :: soil_z_grid(:)
  REAL(r8), allocatable :: soil_t_grid(:,:,:)
  REAL(r8), allocatable :: soil_w_grid(:,:,:)
  REAL(r8), allocatable :: soil_z(:)
  REAL(r8), allocatable :: soil_t(:,:)
  REAL(r8), allocatable :: soil_w(:,:)
#endif

  REAL(r8), allocatable :: z_soisno (:,:)
  REAL(r8), allocatable :: dz_soisno(:,:)

  REAL(r8) :: calday                   !Julian cal day (1.xx to 365.xx)
  INTEGER  :: year, jday, msec         !Julian day and seconds
  INTEGER  :: month, mday              !month and day of month
  INTEGER  :: i,j,k,l,m,npatch,np,nsl  !indices
  INTEGER  :: numpatch_lat(lat_points) !number of patches of grids at lon. strip

  CHARACTER(LEN=256) :: c
  CHARACTER(LEN=255) :: cdate          !CHARACTER for date
  CHARACTER(len=256) :: lndname
  INTEGER iunit
  INTEGER Julian_8day

  INTEGER npft, npc
  INTEGER ncid, landfrac_vid, pctlc_vid, pctelc_vid
  INTEGER pctpft_vid, pctpc_vid
  INTEGER pcturban_vid, pctwater_vid, pctwetland_vid, pctglacier_vid

  REAL(r8) sumpctpft
  REAL(r8), external :: orb_coszen     !cosine of the solar zenith angle


! ----------------------------------------------------------------------
! [1] READ IN LAND INFORMATION
! read time-invariant boundary data on [lon_points] x [lat_points] grid.
! ----------------------------------------------------------------------
! Read in the coordinate of the center of the model grids and area of grid cells
      iunit = 100
      lndname = trim(dir_model_landdata)//'model_lonlat_gridcell.bin'
      print*,trim(lndname)
      OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
      READ(iunit) latixy
      READ(iunit) longxy
      READ(iunit) area_gridcells
      close(iunit)

      ! get grid latitudes and longitudes
      latdeg = latixy(1,:)
      londeg = longxy(:,1)
      
      ! convert latitudes and longitudes from degress to radians
      latixy(:,:) = latixy(:,:)*PI/180. 
      longxy(:,:) = longxy(:,:)*PI/180. 


! ----------------------------------------------------------------------
! [2] MAPPING and ALLOCATE
! Build 1d subgrid patch <-> 2d grid mapping indices and weights
! 
! Build mapping indices and weights: [lon_points]x[lat_points] 2d grid <->
! <-> [numpatch] vector of subgrid patches. 
! The land surface model works by gathering all the land points on a
! [lon_points]x[lat_points] grid into a vector, and then expanded into 
! a vector of [numpatch] subgrid patches, allowing
! for up to [maxpatch=N_land_classification + 1] subgrid patches per land point. 
! [ixy], [jxy], [patch], and [land] are indices for the mapping: 
! [lon_points]x[lat_points] grid <-> [numpatch] vector of subgrid points. 
!
!-----------------------------------------------------------------------
! Find total number of patches [numpatch] allowing for multiple subgrid 
! patches in a grid cell.
! --------------------------------------------------------------------

#ifdef USGS_CLASSIFICATION

! Read in the patch fraction of the lantypes of the gridcells
      allocate (fraction_patches(0:N_land_classification,1:lon_points,1:lat_points))
      lndname = trim(dir_model_landdata)//'model_landtypes.bin'
      print*,trim(lndname)
      OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
      READ(iunit,err=100) fraction_patches
      close(iunit)
      print*,'fraction   =', minval(fraction_patches, mask = fraction_patches .gt. -1.0e30), &
                             maxval(fraction_patches ,mask = fraction_patches .gt. -1.0e30)

      npatch = 0
      numpatch_lat(:) = 0

      DO j = 1, lat_points
         DO i = 1, lon_points
#if(defined LANDONLY)
            DO np = 1, N_land_classification
               IF(fraction_patches(np,i,j)> 0.)THEN
                  npatch = npatch+1 !subgrid patch number
                  numpatch_lat(j) = numpatch_lat(j) + 1
               ENDIF
            ENDDO
#elif(defined LAND_SEA || defined USE_POINT_DATA)
            DO np = 0, N_land_classification
               IF(fraction_patches(np,i,j)> 0.)THEN
                  npatch = npatch+1 !subgrid patch number
                  numpatch_lat(j) = numpatch_lat(j) + 1
               ENDIF
            ENDDO
#endif
         ENDDO
      ENDDO
      numpatch = npatch
      IF(numpatch.ne.sum(numpatch_lat))THEN
         write(6,*) 'Total number of patches NOT as the summation of numpatch_lat'
         CALL abort
      ENDIF
      write(6,*) 'Total land patches = ', numpatch
#endif

! yuan, 07/31/2019: read land grid info from nc file
#ifdef IGBP_CLASSIFICATION
      
      allocate (landfrac(1:lon_points,1:lat_points))
      allocate (pctlc(1:lon_points,1:lat_points,1:N_land_classification))
      allocate (fraction_patches(1:lon_points,1:lat_points,1:N_land_classification))
!TODO: here for year 2005, need a input PARAMETER for a perscribed year
      lndname = trim(dir_model_landdata)//'global_0.5x0.5.MOD2005_V4.5.nc'
      print*,trim(lndname)

      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, "LANDFRAC", landfrac_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_LC", pctlc_vid ) )
      CALL nccheck( nf90_get_var(ncid, landfrac_vid, landfrac) )
      CALL nccheck( nf90_get_var(ncid, pctlc_vid, pctlc) )

      landfrac = landfrac / 100.
      pctlc    = pctlc / 100.

      DO np = 1, N_land_classification
         ! sum(pctlc) = 100%, landfrac: land%
         fraction_patches(:,:,np) = pctlc(:,:,np) * landfrac(:,:)
      ENDDO

      print*,'fraction   =', minval(fraction_patches, mask = fraction_patches .gt. -1.0e30), &
                             maxval(fraction_patches ,mask = fraction_patches .gt. -1.0e30)

      npatch = 0
      numpatch_lat(:) = 0

      ! NOTE: support for land ONLY right now
      DO j = 1, lat_points
         DO i = 1, lon_points
            DO np = 1, N_land_classification
               IF(fraction_patches(i,j,np)> 0.)THEN
                  npatch = npatch+1 !subgrid patch number
                  numpatch_lat(j) = numpatch_lat(j) + 1
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      numpatch = npatch
      IF(numpatch.ne.sum(numpatch_lat))THEN
         write(6,*) 'Total number of patches NOT as the summation of numpatch_lat'
         CALL abort
      ENDIF
      write(6,*) 'Total land patches = ', numpatch

#endif

#ifdef PFT_CLASSIFICATION
      
      allocate (landfrac(1:lon_points,1:lat_points))
      allocate (pctpft(1:lon_points,1:lat_points,0:N_PFT-1))
      allocate (pcturban(1:lon_points,1:lat_points))
      allocate (pctwater(1:lon_points,1:lat_points))
      allocate (pctwetland(1:lon_points,1:lat_points))
      allocate (pctglacier(1:lon_points,1:lat_points))
      lndname = trim(dir_model_landdata)//'global_0.5x0.5.MOD2005_V4.5.nc'
      print*,trim(lndname)

      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, "LANDFRAC", landfrac_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_PFT",  pctpft_vid   ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_URBAN", pcturban_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_WATER", pctwater_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_WETLAND", pctwetland_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_GLACIER", pctglacier_vid ) )

      CALL nccheck( nf90_get_var(ncid, landfrac_vid, landfrac) )
      CALL nccheck( nf90_get_var(ncid, pctpft_vid, pctpft) )
      CALL nccheck( nf90_get_var(ncid, pcturban_vid, pcturban) )
      CALL nccheck( nf90_get_var(ncid, pctwater_vid, pctwater) )
      CALL nccheck( nf90_get_var(ncid, pctwetland_vid, pctwetland) )
      CALL nccheck( nf90_get_var(ncid, pctglacier_vid, pctglacier) )

      landfrac   = landfrac   / 100.
      pctpft     = pctpft     / 100.
      pcturban   = pcturban   / 100.
      pctwater   = pctwater   / 100.
      pctwetland = pctwetland / 100.
      pctglacier = pctglacier / 100.

      npatch = 0
      numpatch_lat(:) = 0
      npft = 0

      ! NOTE: support for land ONLY right now
      DO j = 1, lat_points
         DO i = 1, lon_points

            sumpctpft = sum(pctpft(i,j,:)) 
            IF (sumpctpft > 0.) THEN
               npatch = npatch + 1 !subgrid patch number
               numpatch_lat(j) = numpatch_lat(j) + 1

               DO np = 0, N_PFT-1
                  IF (pctpft(i,j,np) > 0.) THEN
                     npft = npft + 1
                  ENDIF
               ENDDO
            ENDIF

            IF (pcturban(i,j) > 0.) THEN
               npatch = npatch + 1 
               numpatch_lat(j) = numpatch_lat(j) + 1
            ENDIF 
            IF (pctwetland(i,j) > 0.) THEN
               npatch = npatch + 1 
               numpatch_lat(j) = numpatch_lat(j) + 1
            ENDIF 
            IF (pctglacier(i,j) > 0.) THEN
               npatch = npatch + 1 
               numpatch_lat(j) = numpatch_lat(j) + 1
            ENDIF 
            IF (pctwater(i,j) > 0.) THEN
               npatch = npatch + 1 
               numpatch_lat(j) = numpatch_lat(j) + 1
            ENDIF 
         ENDDO
      ENDDO

      numpatch = npatch
      numpft   = npft
      IF(numpatch.ne.sum(numpatch_lat))THEN
         write(6,*) 'Total number of patches NOT as the summation of numpatch_lat'
         CALL abort
      ENDIF
      write(6,*) 'Total land patches = ', numpatch
      write(6,*) 'Total PFT  patches = ', numpft

#endif

#ifdef PC_CLASSIFICATION
      
      allocate (landfrac(1:lon_points,1:lat_points))
      allocate (pctelc(1:lon_points,1:lat_points,1:N_land_classification))
      allocate (pctpc(1:lon_points,1:lat_points,0:N_PFT-1,1:N_land_classification))
      allocate (fraction_patches(1:lon_points,1:lat_points,1:N_land_classification))
      lndname = trim(dir_model_landdata)//'global_0.5x0.5.MOD2005_V4.5.nc'
      print*,trim(lndname)

      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, "LANDFRAC", landfrac_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_eLC",  pctelc_vid  ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_ePFT", pctpc_vid   ) )
      CALL nccheck( nf90_get_var(ncid, landfrac_vid, landfrac) )
      CALL nccheck( nf90_get_var(ncid, pctelc_vid,   pctelc  ) )
      CALL nccheck( nf90_get_var(ncid, pctpc_vid,    pctpc   ) )

      landfrac = landfrac / 100.
      pctelc   = pctelc   / 100.
      pctpc    = pctpc    / 100.

      DO np = 1, N_land_classification
         fraction_patches(:,:,np) = pctelc(:,:,np) * landfrac(:,:)
      ENDDO

      print*,'fraction   =', minval(fraction_patches, mask = fraction_patches .gt. -1.0e30), &
                             maxval(fraction_patches ,mask = fraction_patches .gt. -1.0e30)

      npatch = 0
      numpatch_lat(:) = 0
      npc = 0
      
      ! NOTE: support for land ONLY right now
      DO j = 1, lat_points
         DO i = 1, lon_points
            DO np = 1, N_land_classification
               IF (fraction_patches(i,j,np) > 0.) THEN
                  npatch = npatch + 1 !subgrid patch number
                  numpatch_lat(j) = numpatch_lat(j) + 1
                  IF (patchtypes(np) == 0) THEN
                     npc = npc + 1
                  ENDIF 
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      numpatch = npatch
      numpc = npc
      IF(numpatch.ne.sum(numpatch_lat))THEN
         write(6,*) 'Total number of patches NOT as the summation of numpatch_lat'
         CALL abort
      ENDIF
      write(6,*) 'Total land patches = ', numpatch
      write(6,*) 'Total PC patches = ', numpc

#endif

! --------------------------------------------------------------------
! Allocates memory for CLM 1d [numpatch] variables
! --------------------------------------------------------------------

      CALL allocate_TimeInvariants (lon_points,lat_points)
      CALL allocate_TimeVariables 


! set the lat/lon values
      gridlatd(:) = latdeg(:)
      gridlond(:) = londeg(:)

! --------------------------------------------------------------------
! Build 1d land vector and 1d patch vector mapping components
! --------------------------------------------------------------------

! Determine land vector and patch vector mapping components

#ifdef USGS_CLASSIFICATION
      npatch = 0
      patchfrac(:) = 0.
      l = 0; m = 0
      grid_patch_s(:,:) = -1
      grid_patch_e(:,:) = -1

      DO j = 1, lat_points
         DO i = 1, lon_points

            gridarea(i,j) = area_gridcells(i,j)     !grid cell area

#if(defined LANDONLY)
            DO np = 1, N_land_classification
#elif(defined LAND_SEA || defined USE_POINT_DATA)
            DO np = 0, N_land_classification
#endif
               IF(fraction_patches(np,i,j)>0.)THEN                

                  npatch             = npatch+1                      
                  patch2lon(npatch)  = i  !patch longitude index
                  patch2lat(npatch)  = j  !patch latitude index
                  patchclass(npatch) = np !index of land cover type 
                  patchlatr(npatch)  = latixy(i,j) !latitude in radians
                  patchlonr(npatch)  = longxy(i,j) !longitude in radians

                  patchfrac(npatch)  = fraction_patches(np,i,j) !patch weight
                  patchtype(npatch)  = patchtypes(np)
                  grid_patch_e(i,j)  = npatch

                  IF (l.ne.i .OR. m.ne.j) THEN
                     l = i; m = j; grid_patch_s(i,j) = npatch
                  ENDIF

               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF(numpatch.ne.npatch)THEN
         write(6,*) 'the number of patches is not identical ', numpatch, npatch
         CALL abort
      ENDIF
      
      deallocate (fraction_patches)

#endif

#ifdef IGBP_CLASSIFICATION
      npatch = 0
      patchfrac(:) = 0.
      l = 0; m = 0
      grid_patch_s(:,:) = -1
      grid_patch_e(:,:) = -1

      DO j = 1, lat_points
         DO i = 1, lon_points

            gridarea(i,j) = area_gridcells(i,j) !grid cell area
            
            DO np = 1, N_land_classification
               IF(fraction_patches(i,j,np)>0.)THEN                

                  npatch             = npatch+1                      
                  patch2lon(npatch)  = i  !patch longitude index
                  patch2lat(npatch)  = j  !patch latitude index
                  patchclass(npatch) = np !index of land cover type 
                  patchlatr(npatch)  = latixy(i,j) !latitude in radians
                  patchlonr(npatch)  = longxy(i,j) !longitude in radians

                  patchfrac(npatch)  = fraction_patches(i,j,np) !patch weight
                  patchtype(npatch)  = patchtypes(np)
                  grid_patch_e(i,j)  = npatch

                  IF (l.ne.i .OR. m.ne.j) THEN
                     l = i; m = j; grid_patch_s(i,j) = npatch
                  ENDIF

               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF(numpatch.ne.npatch)THEN
         write(6,*) 'the number of patches is not identical ', numpatch, npatch
         CALL abort
      ENDIF

      deallocate (landfrac)
      deallocate (pctlc)
      deallocate (fraction_patches)

#endif
   
#ifdef PFT_CLASSIFICATION
      npatch = 0
      npft   = 0
      patchfrac(:) = 0.
      l = 0; m = 0
      grid_patch_s(:,:) = -1
      grid_patch_e(:,:) = -1
      patch_pft_s(:)    = -1
      patch_pft_e(:)    = -1

      DO j = 1, lat_points
         DO i = 1, lon_points

            gridarea(i,j) = area_gridcells(i,j) !grid cell area

            sumpctpft = sum(pctpft(i,j,:)) 
            IF (sumpctpft > 0.) THEN

               npatch             = npatch + 1 !subgrid patch number
               patch2lon(npatch)  = i !patch longitude index
               patch2lat(npatch)  = j !patch latitude index
               patchclass(npatch) = 1 !no meaning here
               patchlatr(npatch)  = latixy(i,j) !latitude in radians
               patchlonr(npatch)  = longxy(i,j) !longitude in radians

               patchfrac(npatch)  = sumpctpft*landfrac(i,j) !patch weight
               patchtype(npatch)  = 0                       !soil patch
               grid_patch_e(i,j)  = npatch
               patch_pft_s(npatch)= -1

               IF (l.ne.i .OR. m.ne.j) THEN
                  l = i; m = j; grid_patch_s(i,j) = npatch
               ENDIF

               DO np = 0, N_PFT-1
                  IF (pctpft(i,j,np) > 0.) THEN
                     npft            = npft + 1
                     pftclass(npft)  = np
                     pftfrac(npft)   = pctpft(i,j,np) / sumpctpft
                     pft2patch(npft) = npatch

                     IF (patch_pft_s(npatch) == -1) THEN
                        patch_pft_s(npatch) = npft
                     ENDIF 

                     patch_pft_e(npatch)= npft
                  ENDIF
               ENDDO
            ENDIF
            
            IF (pcturban(i,j) > 0.) THEN
               npatch             = npatch + 1 
               patch2lon(npatch)  = i  !patch longitude index
               patch2lat(npatch)  = j  !patch latitude index
               patchclass(npatch) = 13 !index of land cover type 
               patchlatr(npatch)  = latixy(i,j) !latitude in radians
               patchlonr(npatch)  = longxy(i,j) !longitude in radians

               patchfrac(npatch)  = pcturban(i,j)*landfrac(i,j) !patch weight
               patchtype(npatch)  = 1                           !urban patch
               grid_patch_e(i,j)  = npatch

               IF (l.ne.i .OR. m.ne.j) THEN
                  l = i; m = j; grid_patch_s(i,j) = npatch
               ENDIF
            ENDIF 

            IF (pctwetland(i,j) > 0.) THEN
               npatch             = npatch + 1 
               patch2lon(npatch)  = i  !patch longitude index
               patch2lat(npatch)  = j  !patch latitude index
               patchclass(npatch) = 11 !index of land cover type
               patchlatr(npatch)  = latixy(i,j) !latitude in radians
               patchlonr(npatch)  = longxy(i,j) !longitude in radians

               patchfrac(npatch)  = pctwetland(i,j)*landfrac(i,j) !patch weight
               patchtype(npatch)  = 2                             !wetlant patch
               grid_patch_e(i,j)  = npatch

               IF (l.ne.i .OR. m.ne.j) THEN
                  l = i; m = j; grid_patch_s(i,j) = npatch
               ENDIF
            ENDIF 

            IF (pctglacier(i,j) > 0.) THEN
               npatch             = npatch + 1 
               patch2lon(npatch)  = i  !patch longitude index
               patch2lat(npatch)  = j  !patch latitude index
               patchclass(npatch) = 15 !index of land cover type
               patchlatr(npatch)  = latixy(i,j) !latitude in radians
               patchlonr(npatch)  = longxy(i,j) !longitude in radians

               patchfrac(npatch)  = pctglacier(i,j)*landfrac(i,j) !patch weight
               patchtype(npatch)  = 3                             !glacier patch
               grid_patch_e(i,j)  = npatch

               IF (l.ne.i .OR. m.ne.j) THEN
                  l = i; m = j; grid_patch_s(i,j) = npatch
               ENDIF
            ENDIF 
            
            IF (pctwater(i,j) > 0.) THEN
               npatch             = npatch + 1 
               patch2lon(npatch)  = i  !patch longitude index
               patch2lat(npatch)  = j  !patch latitude index
               patchclass(npatch) = 17 !index of land cover type
               patchlatr(npatch)  = latixy(i,j) !latitude in radians
               patchlonr(npatch)  = longxy(i,j) !longitude in radians

               patchfrac(npatch)  = pctwater(i,j)*landfrac(i,j) !patch weight
               patchtype(npatch)  = 4                           !water patch
               grid_patch_e(i,j)  = npatch

               IF (l.ne.i .OR. m.ne.j) THEN
                  l = i; m = j; grid_patch_s(i,j) = npatch
               ENDIF
            ENDIF 

         ENDDO
      ENDDO

      IF(numpatch.ne.npatch)THEN
         write(6,*) 'the number of patches is not identical ', numpatch, npatch
         CALL abort
      ENDIF
      
      deallocate (landfrac  )
      deallocate (pctpft    )
      deallocate (pcturban  )
      deallocate (pctwater  )
      deallocate (pctwetland)
      deallocate (pctglacier)

#endif

#ifdef PC_CLASSIFICATION
      npatch = 0
      npc = 0
      patchfrac(:) = 0.
      l = 0; m = 0
      grid_patch_s(:,:) = -1
      grid_patch_e(:,:) = -1

      DO j = 1, lat_points
         DO i = 1, lon_points

            gridarea(i,j) = area_gridcells(i,j) !grid cell area

            DO np = 1, N_land_classification
               IF(fraction_patches(i,j,np)>0.)THEN                

                  npatch             = npatch+1                      
                  patch2lon(npatch)  = i  !patch longitude index
                  patch2lat(npatch)  = j  !patch latitude index
                  patchclass(npatch) = np !index of land cover type
                  patchlatr(npatch)  = latixy(i,j) !latitude in radians
                  patchlonr(npatch)  = longxy(i,j) !longitude in radians

                  patchfrac(npatch)  = fraction_patches(i,j,np) !patch weight
                  patchtype(npatch)  = patchtypes(np)
                  grid_patch_e(i,j)  = npatch

                  IF (l.ne.i .OR. m.ne.j) THEN
                     l = i; m = j; grid_patch_s(i,j) = npatch
                  ENDIF

                  ! for soil patches
                  IF (patchtypes(np) == 0) THEN
                     npc              = npc + 1
                     pcfrac(:,npc)    = pctpc(i,j,:,np)
                     patch2pc(npatch) = npc
                     pc2patch(npc)    = npatch
                  ENDIF 

               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF(numpatch.ne.npatch)THEN
         write(6,*) 'the number of patches is not identical ', numpatch, npatch
         CALL abort
      ENDIF

      deallocate (landfrac        )
      deallocate (pctelc          )
      deallocate (pctpc           )
      deallocate (fraction_patches)

#endif

! ---------------------------------------------------------------
! [3] INITIALIZE TIME INVARIANT VARIABLES
! ---------------------------------------------------------------
! 3.1 Define the soil and lake layers's thickness 
! ...............................................................
      allocate ( z_soisno (maxsnl+1:nl_soil,numpatch) )
      allocate ( dz_soisno(maxsnl+1:nl_soil,numpatch) )

    ! ------------------------------------------
    ! Lake depth and layers' thickness
    ! ------------------------------------------
      CALL lakedepth_readin (lon_points,lat_points,dir_model_landdata)

! ...............................................................
! 3.2 Read in the soil parameters of the patches of the gridcells
! ...............................................................

      CALL soil_parameters_readin (lon_points,lat_points,dir_model_landdata)

! ...............................................................
! 3.3 Plant time-invariant variables (based on the look-up tables or global map) 
! ...............................................................

#ifdef USGS_CLASSIFICATION
      ! set canopy top height with look-up tables
      DO i = 1, numpatch
         htop(i) = htop0(patchclass(i))
         hbot(i) = hbot0(patchclass(i))
      ENDDO 
#else
      ! read global tree top height from nc file
      CALL HTOP_readin_nc (lon_points, lat_points, dir_model_landdata)
#endif

! ................................
! 3.4 Initialize TUNABLE constants
! ................................
      zlnd   = 0.01    !Roughness length for soil [m]
      zsno   = 0.0024  !Roughness length for snow [m]
      csoilc = 0.004   !Drag coefficient for soil under canopy [-]
      dewmx  = 0.1     !maximum dew
      wtfact = 0.38    !Maximum saturated fraction (global mean; see Niu et al., 2005)
      capr   = 0.34    !Tuning factor to turn first layer T into surface T
      cnfac  = 0.5     !Crank Nicholson factor between 0 and 1
      ssi    = 0.033   !Irreducible water saturation of snow
      wimp   = 0.05    !Water impremeable if porosity less than wimp
      pondmx = 10.0    !Ponding depth (mm)
      smpmax = -1.5e5  !Wilting point potential in mm
      smpmin = -1.e8   !Restriction for min of soil poten. (mm)
      trsmx0 = 2.e-4   !Max transpiration for moist soil+100% veg. [mm/s]
      tcrit  = 2.5     !critical temp. to determine rain or snow

! ...............................................
! 3.5 Write out as a restart file [histTimeConst]
! ...............................................

      CALL WRITE_TimeInvariants (dir_restart_hist,casename)

      write (6,*)
      write (6,*) ('Successfully to Initialize the Land Time-Invariants')

! ----------------------------------------------------------------------
! [4] INITIALIZE TIME-VARYING VARIABLES 
! as subgrid vectors of length [numpatch]
! initial run: create the time-varying variables based on :
!              i) observation (NOT CODING CURRENTLY), or
!             ii) some already-known information (NO CODING CURRENTLY), or
!            iii) arbitrarily 
! continuation run: time-varying data read in from restart file 
! ----------------------------------------------------------------------

! 4.1 current time of model run
! ............................

      CALL initimetype (greenwich)

      IF(.not. greenwich)THEN
         print *, ".........greenwich false", longxy(1,1)
      ENDIF

      year = idate(1)
      jday = idate(2)
      msec = idate(3)

! ................................
! 4.2 cosine of solar zenith angle 
! ................................
      calday = calendarday(idate, gridlond(1))
      DO i = 1, numpatch
         coszen(i) = orb_coszen(calday,patchlonr(i),patchlatr(i))
      ENDDO

#if(defined SOILINI)
! ...........................................
!4.3 READ in or GUSSES land state information
! ...........................................
!!! PLEASE CHANGE
!!! PLEASE CHANGE
!!! PLEASE CHANGE the root of directory when the soil T and W are ready !!!
     lusoil = 100
     fsoildat_name = trim(fsoildat)//'-'//trim(cdate)
     OPEN(lusoil,file=fsoildat_name,status='old',form='unformatted',action='read')

     read(lusoil) nl_soil_ini

     allocate (snow_d_grid(lon_points,lat_points))
     allocate (snow_d(numpatch))

     allocate (soil_z_grid(nl_soil_ini))
     allocate (soil_t_grid(lon_points,lat_points,nl_soil_ini))
     allocate (soil_w_grid(lon_points,lat_points,nl_soil_ini))

     allocate (soil_z(nl_soil_ini))
     allocate (soil_t(nl_soil_ini,numpatch))
     allocate (soil_w(nl_soil_ini,numpatch))

     read(lusoil) soil_z_grid  !soil layer node depth (m)
     read(lusoil) soil_t_grid  !soil layer temperature (K)
     read(lusoil) soil_w_grid  !soil layer wetness (-)
     read(lusoil) snow_d_grid  !snow depth (m)

     close (lusoil)

     soil_z(:) = soil_z_grid(:) 
     DO npatch = 1, numpatch
        i = patch2lon(npatch)
        j = patch2lat(npatch)

        snow_d(npatch) = snow_d_grid(i,j)
        DO l = 1, nl_soil_ini
           soil_t(l,npatch) = soil_t_grid(i,j,l)
           soil_w(l,npatch) = soil_w_grid(i,j,l)
        ENDDO
     ENDDO
#endif

! ...................
! 4.4 LEAF area index
! ...................
#if(defined DYN_PHENOLOGY)
    ! CREAT fraction of vegetation cover, greenness, leaf area index, stem index
      lai(:)=0.0; sai(:)=0.0; green(:)=0.0; fveg(:)=0.0
      DO i = 1, numpatch
#if(defined SOILINI)
         DO l = 1, nl_soil
            t_soisno(l,i) = soil_t(nl_soil_ini,i)
         ENDDO
#else
         t_soisno(1:,i) = 283.
#endif
       ! Call Ecological Model() 
         m = patchclass(i)
         IF(m>0)THEN
            CALL lai_empirical (m,nl_soil,rootfr(1:,m),&
                                t_soisno(1:,i),lai(i),sai(i),fveg(i),green(i))

         ENDIF
      ENDDO
#else

#ifdef USGS_CLASSIFICATION
    ! READ in Leaf area index and stem area index
      Julian_8day = int(calendarday(idate)-1)/8*8 + 1
      CALL LAI_readin (lon_points,lat_points,&
                       Julian_8day,numpatch,dir_model_landdata)
#else
! yuan, 08/03/2019: read global LAI/SAI data
      CALL julian2monthday (year, jday, month, mday)
      CALL LAI_readin_nc (lon_points, lat_points, month, dir_model_landdata)
#endif

#endif
! ..............................................................................
! 4.5 initialize time-varying variables, as subgrid vectors of length [numpatch]
! ..............................................................................
      DO i = 1, numpatch
         z_soisno (1:nl_soil ,i) = z_soi (1:nl_soil)
         dz_soisno(1:nl_soil ,i) = dz_soi(1:nl_soil)
      ENDDO

! 03/09/2020, yuan: TODO, may change the below from serial->parallel
#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE (m) &
!$OMP SCHEDULE(STATIC, 1)
#endif
      DO i = 1, numpatch
         m = patchclass(i)
      CALL iniTimeVar (i, patchtype(i)&
          ,porsl(1:,i)&
          ,soil_s_v_alb(i),soil_d_v_alb(i),soil_s_n_alb(i),soil_d_n_alb(i)&
          ,z0m(i),zlnd,chil(m),rho(1:,1:,m),tau(1:,1:,m)&
          ,z_soisno(maxsnl+1:,i),dz_soisno(maxsnl+1:,i)&
          ,t_soisno(maxsnl+1:,i),wliq_soisno(maxsnl+1:,i),wice_soisno(maxsnl+1:,i)&
          ,zwt(i),wa(i)&
          ,t_grnd(i),tleaf(i),ldew(i),sag(i),scv(i)&
          ,snowdp(i),fveg(i),fsno(i),sigf(i),green(i),lai(i),sai(i),coszen(i)&
          ,alb(1:,1:,i),ssun(1:,1:,i),ssha(1:,1:,i)&
          ,thermk(i),extkb(i),extkd(i)&
          ,trad(i),tref(i),qref(i),rst(i),emis(i),zol(i),rib(i)&
          ,ustar(i),qstar(i),tstar(i),fm(i),fh(i),fq(i)&
#if(defined SOILINI)
          ,nl_soil_ini,soil_z,soil_t(1:,i),soil_w(1:,i),snow_d(i))
#else
          )
#endif
      ENDDO

! 03/13/2020, yuan: 
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      DO i = 1, numpatch
         z_sno (maxsnl+1:0,i) = z_soisno (maxsnl+1:0,i)
         dz_sno(maxsnl+1:0,i) = dz_soisno(maxsnl+1:0,i)
      ENDDO

    ! ------------------------------------------
    ! PLEASE  
    ! PLEASE UPDATE
    ! PLEASE UPDATE when have the observed lake status
      t_lake      (:,:) = 285.
      lake_icefrac(:,:) = 0.
    ! ------------------------------------------

! ...............................................................
! 4.6 Write out the model variables for restart run [histTimeVar]
! ...............................................................

      CALL WRITE_TimeVariables (idate,dir_restart_hist,casename)

      write (6,*)
      write (6,*) ('Successfully to Initialize the Land Time-Vraying Variables')


! ----------------------------------------------------------------------
! average subgrid albedos, srf temperature, etc. for atmospheric model
! ----------------------------------------------------------------------
      tg_xy   (:,:) = 0.0
      albvb_xy(:,:) = 0.0
      albvd_xy(:,:) = 0.0
      albnb_xy(:,:) = 0.0
      albnd_xy(:,:) = 0.0
      trad_xy (:,:) = 0.0
      rib_xy  (:,:) = 0.0
      fm_xy   (:,:) = 0.0
      fh_xy   (:,:) = 0.0
      fq_xy   (:,:) = 0.0

      DO npatch = 1, numpatch
         i = patch2lon(npatch)
         j = patch2lat(npatch)
            tg_xy(i,j) =    tg_xy(i,j) + patchfrac(npatch)*t_grnd (npatch)   !
         albvb_xy(i,j) = albvb_xy(i,j) + patchfrac(npatch)*alb(1,1,npatch)   !(v=visble, b=direct beam)
         albvd_xy(i,j) = albvd_xy(i,j) + patchfrac(npatch)*alb(1,2,npatch)   !(v=visble, d=diffuse)
         albnb_xy(i,j) = albnb_xy(i,j) + patchfrac(npatch)*alb(2,1,npatch)   !(n=near infrared, b=direct beam)
         albnd_xy(i,j) = albnd_xy(i,j) + patchfrac(npatch)*alb(2,2,npatch)   !(n=near infrared, d=diffuse)
         trad_xy (i,j) =  trad_xy(i,j) + patchfrac(npatch)*t_grnd (npatch)   !

         rib_xy(i,j) = -0.1       !
         fm_xy (i,j) = alog(30.)  !
         fh_xy (i,j) = alog(30.)  !
         fq_xy (i,j) = alog(30.)  !
      ENDDO

    ! --------------------------------------------------
    ! Deallocates memory for CLM 1d [numpatch] variables
    ! --------------------------------------------------

      CALL deallocate_TimeInvariants
      CALL deallocate_TimeVariables

      deallocate (z_soisno )
      deallocate (dz_soisno)

#if(defined SOILINI)
      deallocate (snow_d_grid)
      deallocate (snow_d)
      deallocate (soil_z_grid)
      deallocate (soil_t_grid)
      deallocate (soil_w_grid)
      deallocate (soil_z)
      deallocate (soil_t)
      deallocate (soil_w)
#endif


      go to 1000
100   print 101,lndname
101   format(' error occured on file: ',a50)
1000  continue


END SUBROUTINE initialize
! --------------------------------------------------
! EOP
