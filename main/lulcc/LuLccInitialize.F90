#include <define.h>

SUBROUTINE LuLccInitialize (casename,dir_srfdata,dir_restart,&
                            nam_srfdata,nam_urbdata,idate,greenwich)

! ======================================================================
!
! Initialization routine for Land-use-Land-cover-change (LuLcc) case
!
! ======================================================================
   USE precision
   USE GlobalVars
   USE PhysicalConstants
   USE MOD_TimeInvariants
   USE MOD_TimeVariables
   USE MOD_PFTimeInvars
   USE MOD_PCTimeInvars
   USE MOD_UrbanTimeInvars
   USE MOD_PFTimeVars
   USE MOD_PCTimeVars
   USE MOD_UrbanTimeVars
   USE MOD_LuLccTMatrix
   USE LC_Const
   USE PFT_Const
   USE timemanager
   USE ncio
   USE netcdf
   USE omp_lib

   IMPLICIT NONE

! ----------------------------------------------------------------------
   CHARACTER(LEN=256), intent(in) :: casename      !casename name
   CHARACTER(LEN=256), intent(in) :: dir_srfdata   !surface data directory
   CHARACTER(LEN=256), intent(in) :: dir_restart   !case restart data directory
   CHARACTER(LEN=256), intent(in) :: nam_srfdata   !surface data filename
   CHARACTER(LEN=256), intent(in) :: nam_urbdata   !urban data filename
   LOGICAL, intent(in)    :: greenwich   !true: greenwich time, false: local time
   INTEGER, intent(inout) :: idate(3)    !year, julian day, seconds of the starting time

! ------------------------ local variables -----------------------------
! surface classification and soil information

  REAL(r8) latixy  (lon_points,lat_points)          !latitude in radians
  REAL(r8) longxy  (lon_points,lat_points)          !longitude in radians
  REAL(r8) latdeg  (lat_points)                     !latitude in degree
  REAL(r8) londeg  (lon_points)                     !longitude in degree
  REAL(r8) area_gridcells (lon_points,lat_points)   !area of gridcells (km^2)

  REAL(r8), allocatable :: landfrac(:,:)            !land fractional cover
  REAL(r8), allocatable :: pctlc(:,:,:)             !percent each land cover TYPE
  REAL(r8), allocatable :: pctpft(:,:,:)            !percent PFT
  REAL(r8), allocatable :: pcturban(:,:)            !percent urban
  REAL(r8), allocatable :: pctwater(:,:)            !percent water body
  REAL(r8), allocatable :: pctwetland(:,:)          !percent wetland
  REAL(r8), allocatable :: pctglacier(:,:)          !percent glacier
  REAL(r8), allocatable :: pctpc(:,:,:,:)           !percent each component in PC
  REAL(r8), allocatable :: fraction_patches(:,:,:)  !fraction of the patch of landtypes in gridcells

! 定义城市读取数据变量
#if(defined URBAN_MODEL)
  INTEGER urbanpct_vid
  REAL(r8), allocatable :: urbanpct(:,:,:) !percent urban TYPE (density)
#endif

  REAL(r8), allocatable :: z_soisno (:,:)
  REAL(r8), allocatable :: dz_soisno(:,:)

  REAL(r8) :: calday                   !Julian cal day (1.xx to 365.xx)
  INTEGER  :: year, jday, msec         !Julian day and seconds
  INTEGER  :: month, mday              !month and day of month
  INTEGER  :: i,j,l,m,u,t,npatch,np    !indices
  INTEGER  :: numpatch_lat(lat_points) !number of patches of grids at lon. strip

  CHARACTER(LEN=255) :: cdate          !CHARACTER for date
  CHARACTER(len=256) :: cyear          !character for year
  CHARACTER(len=256) :: lndname        !land surface data file name
  CHARACTER(LEN=256) :: finfolist      !file name of run information
  INTEGER iunit

  INTEGER npft, npc, nurb
  INTEGER ncid, landfrac_vid, pctlc_vid
  INTEGER pctpft_vid, pctpc_vid
  INTEGER pcturban_vid, pctwater_vid, pctwetland_vid, pctglacier_vid

  REAL(r8) sumpctpft
  REAL(r8), external :: orb_coszen     !cosine of the solar zenith angle


! initial time of model run
! ............................
      CALL adj2begin(idate)

      year = idate(1)
      jday = idate(2)
      msec = idate(3)

#ifdef USGS_CLASSIFICATION
      cyear = ''
#else
      write(cyear,'(i4.4)') year
#endif

! ----------------------------------------------------------------------
! [1] READ IN LAND INFORMATION
! read time-invariant boundary data on [lon_points] x [lat_points] grid.
! ----------------------------------------------------------------------
! Read in the coordinate of the center of the model grids and area of grid cells
      iunit = 100
      lndname = trim(dir_srfdata)//trim(cyear)//'/model_lonlat_gridcell.bin'
      print*,trim(lndname)
      open(iunit,file=trim(lndname),form='unformatted',status='old')
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
      lndname = trim(dir_srfdata)//'model_landtypes.bin'
      print*,trim(lndname)
      open(iunit,file=trim(lndname),form='unformatted',status='old')
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
               IF (fraction_patches(np,i,j) > 0.) THEN
                  npatch = npatch+1 !subgrid patch number
                  numpatch_lat(j) = numpatch_lat(j) + 1
               ENDIF
            ENDDO
#else
            DO np = 0, N_land_classification
               IF (fraction_patches(np,i,j) > 0.) THEN
                  npatch = npatch+1 !subgrid patch number
                  numpatch_lat(j) = numpatch_lat(j) + 1
               ENDIF
            ENDDO
#endif
         ENDDO
      ENDDO
      numpatch = npatch
      IF (numpatch .ne. sum(numpatch_lat)) THEN
         write(6,*) 'Total number of patches NOT as the summation of numpatch_lat'
         CALL abort
      ENDIF
      write(6,*) 'Total land patches = ', numpatch
#endif

! 添加城市数据读取，目前仅支持MODIS IGBP数据
#if(!defined USGS_CLASSIFICATION && defined URBAN_MODEL)
      allocate (urbanpct(1:lon_points,1:lat_points,N_URB))
      lndname = trim(dir_srfdata)//trim(cyear)//'/'//trim(nam_urbdata)
      print*,trim(lndname)

      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, "URBAN_PCT", urbanpct_vid ) )
      CALL nccheck( nf90_get_var(ncid, urbanpct_vid, urbanpct) )
      CALL nccheck( nf90_close(ncid) )
      urbanpct = urbanpct / 100.
#endif

! yuan, 07/31/2019: read land grid info from nc file
#ifdef IGBP_CLASSIFICATION

      allocate (landfrac(1:lon_points,1:lat_points))
      allocate (pctlc   (1:lon_points,1:lat_points,1:N_land_classification))
      allocate (fraction_patches(1:lon_points,1:lat_points,1:N_land_classification))
      lndname = trim(dir_srfdata)//trim(cyear)//'/'//trim(nam_srfdata)
      print*,trim(lndname)

      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, "LANDFRAC", landfrac_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_LC",   pctlc_vid   ) )
      CALL nccheck( nf90_get_var(ncid, landfrac_vid, landfrac) )
      CALL nccheck( nf90_get_var(ncid, pctlc_vid,    pctlc   ) )
      CALL nccheck( nf90_close(ncid) )

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
      nurb = 0

      ! NOTE: support for land ONLY right now
      DO j = 1, lat_points
         DO i = 1, lon_points
            DO np = 1, N_land_classification
               IF (fraction_patches(i,j,np) > 0.) THEN
#ifndef URBAN_MODEL
                  npatch = npatch+1 !subgrid patch number
                  numpatch_lat(j) = numpatch_lat(j) + 1
#else
                  IF (np == URBAN) THEN
                     ! 进行细分，对不同城市类型计数
                     DO t = 1, N_URB
                        IF (urbanpct(i,j,t) > 0.) THEN
                           npatch = npatch+1 !subgrid patch number
                           numpatch_lat(j) = numpatch_lat(j) + 1
                           nurb = nurb + 1
                        ENDIF
                     ENDDO
                  ELSE
                     npatch = npatch+1 !subgrid patch number
                     numpatch_lat(j) = numpatch_lat(j) + 1
                  ENDIF
#endif
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      numpatch = npatch
      numurban = nurb
      IF (numpatch .ne. sum(numpatch_lat)) THEN
         write(6,*) 'Total number of patches NOT as the summation of numpatch_lat'
         CALL abort
      ENDIF
      write(6,*) 'Total land patches = ', numpatch

#endif

#ifdef PFT_CLASSIFICATION

      allocate (landfrac  (1:lon_points,1:lat_points))
      allocate (pcturban  (1:lon_points,1:lat_points))
      allocate (pctwater  (1:lon_points,1:lat_points))
      allocate (pctwetland(1:lon_points,1:lat_points))
      allocate (pctglacier(1:lon_points,1:lat_points))
      allocate (pctpft    (1:lon_points,1:lat_points,0:N_PFT-1))
      lndname = trim(dir_srfdata)//trim(cyear)//'/'//trim(nam_srfdata)
      print*,trim(lndname)

      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, "LANDFRAC",    landfrac_vid  ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_PFT",     pctpft_vid    ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_URBAN",   pcturban_vid  ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_WATER",   pctwater_vid  ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_WETLAND", pctwetland_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_GLACIER", pctglacier_vid) )

      CALL nccheck( nf90_get_var(ncid, landfrac_vid,   landfrac  ) )
      CALL nccheck( nf90_get_var(ncid, pctpft_vid,     pctpft    ) )
      CALL nccheck( nf90_get_var(ncid, pcturban_vid,   pcturban  ) )
      CALL nccheck( nf90_get_var(ncid, pctwater_vid,   pctwater  ) )
      CALL nccheck( nf90_get_var(ncid, pctwetland_vid, pctwetland) )
      CALL nccheck( nf90_get_var(ncid, pctglacier_vid, pctglacier) )
      CALL nccheck( nf90_close(ncid) )

      landfrac   = landfrac   / 100.
      pctpft     = pctpft     / 100.
      pcturban   = pcturban   / 100.
      pctwater   = pctwater   / 100.
      pctwetland = pctwetland / 100.
      pctglacier = pctglacier / 100.

      npatch = 0
      numpatch_lat(:) = 0
      npft = 0
      nurb = 0

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
#ifndef URBAN_MODEL
               npatch = npatch + 1
               numpatch_lat(j) = numpatch_lat(j) + 1
#else
               ! 进行细分，对不同城市类型计数
               DO t = 1, N_URB
                  IF (urbanpct(i,j,t) > 0.) THEN
                     npatch = npatch + 1
                     numpatch_lat(j) = numpatch_lat(j) + 1
                     nurb = nurb + 1
                  ENDIF
               ENDDO
#endif
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
      numurban = nurb
      IF (numpatch .ne. sum(numpatch_lat)) THEN
         write(6,*) 'Total number of patches NOT as the summation of numpatch_lat'
         CALL abort
      ENDIF
      write(6,*) 'Total land patches = ', numpatch
      write(6,*) 'Total PFT  patches = ', numpft

#endif

#ifdef PC_CLASSIFICATION

      allocate (landfrac(1:lon_points,1:lat_points))
      allocate (pctlc   (1:lon_points,1:lat_points,1:N_land_classification))
      allocate (pctpc   (1:lon_points,1:lat_points,0:N_PFT-1,1:N_land_classification))
      allocate (fraction_patches(1:lon_points,1:lat_points,1:N_land_classification))
      lndname = trim(dir_srfdata)//trim(cyear)//'/'//trim(nam_srfdata)
      print*,trim(lndname)

      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, "LANDFRAC", landfrac_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_LC",   pctlc_vid   ) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_PC",   pctpc_vid   ) )
      CALL nccheck( nf90_get_var(ncid, landfrac_vid, landfrac) )
      CALL nccheck( nf90_get_var(ncid, pctlc_vid,    pctlc   ) )
      CALL nccheck( nf90_get_var(ncid, pctpc_vid,    pctpc   ) )
      CALL nccheck( nf90_close(ncid) )

      landfrac = landfrac / 100.
      pctlc    = pctlc    / 100.
      pctpc    = pctpc    / 100.

      DO np = 1, N_land_classification
         fraction_patches(:,:,np) = pctlc(:,:,np) * landfrac(:,:)
      ENDDO

      print*,'fraction   =', minval(fraction_patches, mask = fraction_patches .gt. -1.0e30), &
                             maxval(fraction_patches ,mask = fraction_patches .gt. -1.0e30)

      npatch = 0
      numpatch_lat(:) = 0
      npc = 0
      nurb = 0

      ! NOTE: support for land ONLY right now
      DO j = 1, lat_points
         DO i = 1, lon_points
            DO np = 1, N_land_classification
               IF (fraction_patches(i,j,np) > 0.) THEN
                  IF (patchtypes(np) == 0) npc = npc + 1
#ifndef URBAN_MODEL
                  npatch = npatch + 1 !subgrid patch number
                  numpatch_lat(j) = numpatch_lat(j) + 1
#else
                  IF (np == URBAN) THEN
                     ! 需要进行细分，对不同城市类型计数
                     DO t = 1, N_URB
                        IF (urbanpct(i,j,t) > 0.) THEN
                           npatch = npatch + 1 !subgrid patch number
                           numpatch_lat(j) = numpatch_lat(j) + 1
                           nurb = nurb + 1
                        ENDIF
                     ENDDO
                  ELSE
                     npatch = npatch + 1 !subgrid patch number
                     numpatch_lat(j) = numpatch_lat(j) + 1
                  ENDIF
#endif
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      numpatch = npatch
      numpc    = npc
      numurban = nurb
      IF (numpatch .ne. sum(numpatch_lat)) THEN
         write(6,*) 'Total number of patches NOT as the summation of numpatch_lat'
         CALL abort
      ENDIF
      write(6,*) 'Total land patches = ', numpatch
      write(6,*) 'Total PC patches = ', numpc

#endif

! --------------------------------------------------------------------
! Deallocates memory for CLM 1d [numpatch] variables
! --------------------------------------------------------------------

      CALL deallocate_TimeInvariants
      CALL deallocate_TimeVariables

! --------------------------------------------------------------------
! Allocates memory for CLM 1d [numpatch] variables
! --------------------------------------------------------------------

      CALL allocate_TimeInvariants
      CALL allocate_TimeVariables

! set the lat/lon values
      gridlatd(:) = latdeg(:)
      gridlond(:) = londeg(:)

#if(defined URBAN_MODEL)
      patch2urb(:) = -1
#endif

! --------------------------------------------------------------------
! Build 1d land vector and 1d patch vector mapping components
! --------------------------------------------------------------------

      grid_patch_s(:,:) = -1
      grid_patch_e(:,:) = -1

! Determine land vector and patch vector mapping components

#ifdef USGS_CLASSIFICATION
      npatch = 0
      patchfrac(:) = 0.
      l = 0; m = 0

      DO j = 1, lat_points
         DO i = 1, lon_points

            gridarea(i,j) = area_gridcells(i,j)     !grid cell area

#if(defined LANDONLY)
            DO np = 1, N_land_classification
#else
            DO np = 0, N_land_classification
#endif
               IF (fraction_patches(np,i,j) > 0.) THEN

                  npatch             = npatch+1
                  patch2lon(npatch)  = i  !patch longitude index
                  patch2lat(npatch)  = j  !patch latitude index
                  patchclass(npatch) = np !index of land cover type
                  patchlatr(npatch)  = latixy(i,j) !latitude in radians
                  patchlonr(npatch)  = longxy(i,j) !longitude in radians

                  patchfrac(npatch)  = fraction_patches(np,i,j) !patch weight
                  patchtype(npatch)  = patchtypes(np)
                  grid_patch_e(i,j)  = npatch

                  IF (l.ne.i .or. m.ne.j) THEN
                     l = i; m = j; grid_patch_s(i,j) = npatch
                  ENDIF

               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF (numpatch .ne. npatch) THEN
         write(6,*) 'the number of patches is not identical ', numpatch, npatch
         CALL abort
      ENDIF

      deallocate (fraction_patches)

#endif

#ifdef IGBP_CLASSIFICATION
      npatch = 0
      nurb = 0
      patchfrac(:) = 0.
      l = 0; m = 0

      DO j = 1, lat_points
         DO i = 1, lon_points

            gridarea(i,j) = area_gridcells(i,j) !grid cell area

            DO np = 1, N_land_classification
               IF (fraction_patches(i,j,np) > 0.) THEN
#ifndef URBAN_MODEL
                  npatch             = npatch+1
                  patch2lon(npatch)  = i  !patch longitude index
                  patch2lat(npatch)  = j  !patch latitude index
                  patchclass(npatch) = np !index of land cover type
                  patchlatr(npatch)  = latixy(i,j) !latitude in radians
                  patchlonr(npatch)  = longxy(i,j) !longitude in radians

                  patchfrac(npatch)  = fraction_patches(i,j,np) !patch weight
                  patchtype(npatch)  = patchtypes(np)
                  grid_patch_e(i,j)  = npatch

                  IF (l.ne.i .or. m.ne.j) THEN
                     l = i; m = j; grid_patch_s(i,j) = npatch
                  ENDIF
#else
                  IF (np == URBAN) THEN
                     DO t = 1, N_URB
                        IF (urbanpct(i,j,t) > 0.) THEN
                           npatch             = npatch+1
                           patch2lon(npatch)  = i  !patch longitude index
                           patch2lat(npatch)  = j  !patch latitude index
                           patchclass(npatch) = np !index of land cover type
                           patchlatr(npatch)  = latixy(i,j) !latitude in radians
                           patchlonr(npatch)  = longxy(i,j) !longitude in radians

                           patchfrac(npatch)  = fraction_patches(i,j,np)*urbanpct(i,j,t) !patch weight
                           patchtype(npatch)  = patchtypes(np)
                           grid_patch_e(i,j)  = npatch

                           IF (l.ne.i .or. m.ne.j) THEN
                              l = i; m = j; grid_patch_s(i,j) = npatch
                           ENDIF

                           nurb = nurb + 1
                           urbclass(nurb)    = t
                           patch2urb(npatch) = nurb
                           urb2patch(nurb)   = npatch
                        ENDIF
                     ENDDO
                  ELSE
                     npatch             = npatch+1
                     patch2lon(npatch)  = i  !patch longitude index
                     patch2lat(npatch)  = j  !patch latitude index
                     patchclass(npatch) = np !index of land cover type
                     patchlatr(npatch)  = latixy(i,j) !latitude in radians
                     patchlonr(npatch)  = longxy(i,j) !longitude in radians

                     patchfrac(npatch)  = fraction_patches(i,j,np) !patch weight
                     patchtype(npatch)  = patchtypes(np)
                     grid_patch_e(i,j)  = npatch

                     IF (l.ne.i .or. m.ne.j) THEN
                        l = i; m = j; grid_patch_s(i,j) = npatch
                     ENDIF
                  ENDIF
#endif
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF (numpatch .ne. npatch) THEN
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
      nurb   = 0
      patchfrac(:) = 0.
      l = 0; m = 0
      patch_pft_s(:) = -1
      patch_pft_e(:) = -1

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

               IF (l.ne.i .or. m.ne.j) THEN
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
#ifndef URBAN_MODEL
               npatch             = npatch + 1
               patch2lon(npatch)  = i  !patch longitude index
               patch2lat(npatch)  = j  !patch latitude index
               patchclass(npatch) = URBAN !index of land cover type
               patchlatr(npatch)  = latixy(i,j) !latitude in radians
               patchlonr(npatch)  = longxy(i,j) !longitude in radians

               patchfrac(npatch)  = pcturban(i,j)*landfrac(i,j) !patch weight
               patchtype(npatch)  = 1                           !urban patch
               grid_patch_e(i,j)  = npatch

               IF (l.ne.i .or. m.ne.j) THEN
                  l = i; m = j; grid_patch_s(i,j) = npatch
               ENDIF
#else
               DO t = 1, N_URB
                  IF (urbanpct(i,j,t) > 0.) THEN
                     npatch             = npatch + 1
                     patch2lon(npatch)  = i  !patch longitude index
                     patch2lat(npatch)  = j  !patch latitude index
                     patchclass(npatch) = URBAN !index of land cover type
                     patchlatr(npatch)  = latixy(i,j) !latitude in radians
                     patchlonr(npatch)  = longxy(i,j) !longitude in radians

                     patchfrac(npatch)  = pcturban(i,j)*landfrac(i,j)*urbanpct(i,j,t) !patch weight
                     patchtype(npatch)  = 1                           !urban patch
                     grid_patch_e(i,j)  = npatch

                     IF (l.ne.i .or. m.ne.j) THEN
                        l = i; m = j; grid_patch_s(i,j) = npatch
                     ENDIF

                     nurb = nurb + 1
                     urbclass(nurb)    = t
                     patch2urb(npatch) = nurb
                     urb2patch(nurb)   = npatch
                  ENDIF
               ENDDO
#endif
            ENDIF

            IF (pctwetland(i,j) > 0.) THEN
               npatch             = npatch + 1
               patch2lon(npatch)  = i  !patch longitude index
               patch2lat(npatch)  = j  !patch latitude index
               patchclass(npatch) = WETLAND !index of land cover type
               patchlatr(npatch)  = latixy(i,j) !latitude in radians
               patchlonr(npatch)  = longxy(i,j) !longitude in radians

               patchfrac(npatch)  = pctwetland(i,j)*landfrac(i,j) !patch weight
               patchtype(npatch)  = 2                             !wetlant patch
               grid_patch_e(i,j)  = npatch

               IF (l.ne.i .or. m.ne.j) THEN
                  l = i; m = j; grid_patch_s(i,j) = npatch
               ENDIF
            ENDIF

            IF (pctglacier(i,j) > 0.) THEN
               npatch             = npatch + 1
               patch2lon(npatch)  = i  !patch longitude index
               patch2lat(npatch)  = j  !patch latitude index
               patchclass(npatch) = GLACIER !index of land cover type
               patchlatr(npatch)  = latixy(i,j) !latitude in radians
               patchlonr(npatch)  = longxy(i,j) !longitude in radians

               patchfrac(npatch)  = pctglacier(i,j)*landfrac(i,j) !patch weight
               patchtype(npatch)  = 3                             !glacier patch
               grid_patch_e(i,j)  = npatch

               IF (l.ne.i .or. m.ne.j) THEN
                  l = i; m = j; grid_patch_s(i,j) = npatch
               ENDIF
            ENDIF

            IF (pctwater(i,j) > 0.) THEN
               npatch             = npatch + 1
               patch2lon(npatch)  = i  !patch longitude index
               patch2lat(npatch)  = j  !patch latitude index
               patchclass(npatch) = WATERBODY !index of land cover type
               patchlatr(npatch)  = latixy(i,j) !latitude in radians
               patchlonr(npatch)  = longxy(i,j) !longitude in radians

               patchfrac(npatch)  = pctwater(i,j)*landfrac(i,j) !patch weight
               patchtype(npatch)  = 4                           !water patch
               grid_patch_e(i,j)  = npatch

               IF (l.ne.i .or. m.ne.j) THEN
                  l = i; m = j; grid_patch_s(i,j) = npatch
               ENDIF
            ENDIF

         ENDDO
      ENDDO

      IF (numpatch .ne. npatch) THEN
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
      nurb = 0
      patchfrac(:) = 0.
      patch2pc(:)  = -1
      l = 0; m = 0

      DO j = 1, lat_points
         DO i = 1, lon_points

            gridarea(i,j) = area_gridcells(i,j) !grid cell area

            DO np = 1, N_land_classification
               IF (fraction_patches(i,j,np) > 0.) THEN
#ifndef URBAN_MODEL
                  npatch             = npatch+1
                  patch2lon(npatch)  = i  !patch longitude index
                  patch2lat(npatch)  = j  !patch latitude index
                  patchclass(npatch) = np !index of land cover type
                  patchlatr(npatch)  = latixy(i,j) !latitude in radians
                  patchlonr(npatch)  = longxy(i,j) !longitude in radians

                  patchfrac(npatch)  = fraction_patches(i,j,np) !patch weight
                  patchtype(npatch)  = patchtypes(np)
                  grid_patch_e(i,j)  = npatch

                  IF (l.ne.i .or. m.ne.j) THEN
                     l = i; m = j; grid_patch_s(i,j) = npatch
                  ENDIF

                  ! for soil patches
                  IF (patchtypes(np) == 0) THEN
                     npc              = npc + 1
                     pcfrac(:,npc)    = pctpc(i,j,:,np)
                     patch2pc(npatch) = npc
                     pc2patch(npc)    = npatch
                  ENDIF
#else
                  ! for urban patches
                  IF (np == URBAN) THEN
                     DO t = 1, N_URB
                        IF (urbanpct(i,j,t) > 0.) THEN
                           npatch             = npatch+1
                           patch2lon(npatch)  = i  !patch longitude index
                           patch2lat(npatch)  = j  !patch latitude index
                           patchclass(npatch) = np !index of land cover type
                           patchlatr(npatch)  = latixy(i,j) !latitude in radians
                           patchlonr(npatch)  = longxy(i,j) !longitude in radians

                           patchfrac(npatch)  = fraction_patches(i,j,np)*urbanpct(i,j,t) !patch weight
                           patchtype(npatch)  = patchtypes(np)
                           grid_patch_e(i,j)  = npatch

                           IF (l.ne.i .or. m.ne.j) THEN
                              l = i; m = j; grid_patch_s(i,j) = npatch
                           ENDIF

                           nurb = nurb + 1
                           urbclass(nurb)    = t
                           patch2urb(npatch) = nurb
                           urb2patch(nurb)   = npatch
                        ENDIF
                     ENDDO
                  ELSE
                     npatch             = npatch+1
                     patch2lon(npatch)  = i  !patch longitude index
                     patch2lat(npatch)  = j  !patch latitude index
                     patchclass(npatch) = np !index of land cover type
                     patchlatr(npatch)  = latixy(i,j) !latitude in radians
                     patchlonr(npatch)  = longxy(i,j) !longitude in radians

                     patchfrac(npatch)  = fraction_patches(i,j,np) !patch weight
                     patchtype(npatch)  = patchtypes(np)
                     grid_patch_e(i,j)  = npatch

                     IF (l.ne.i .or. m.ne.j) THEN
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
#endif
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF (numpatch .ne. npatch) THEN
         write(6,*) 'the number of patches is not identical ', numpatch, npatch
         CALL abort
      ENDIF

      deallocate (landfrac        )
      deallocate (pctlc           )
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
      CALL lakedepth_readin (dir_srfdata,year)

! ...............................................................
! 3.2 Read in the soil parameters of the patches of the gridcells
! ...............................................................

      CALL soil_parameters_readin (dir_srfdata,year)

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
      CALL HTOP_readin_nc (dir_srfdata, nam_srfdata, year)
#endif

! ...............................................................
! 3.3+ Urban time-invariant variables (based on the look-up tables or global map)
! CALL Urban_readin_nc() !包含htop, fveg. TODO: 设置lakedepth为常数1m?
! ...............................................................

#ifdef URBAN_MODEL
      CALL Urban_readin_nc (dir_srfdata, nam_urbdata, year)
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

      CALL WRITE_TimeInvariants (year,dir_restart,casename)

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

      IF (.not. greenwich) THEN
         print *, ".........greenwich false", longxy(1,1)
      ENDIF

! ................................
! 4.2 cosine of solar zenith angle
! ................................
      calday = calendarday(idate, gridlond(1))
      DO i = 1, numpatch
         coszen(i) = orb_coszen(calday,patchlonr(i),patchlatr(i))
      ENDDO

! ...................
! 4.4 LEAF area index
! ...................
#if(defined DYN_PHENOLOGY)
    ! CREAT fraction of vegetation cover, greenness, leaf area index, stem index
      lai(:)=0.0; sai(:)=0.0; green(:)=0.0; fveg(:)=0.0
      DO i = 1, numpatch
         t_soisno(1:,i) = 283.

       ! Call Ecological Model()
         m = patchclass(i)
         IF (m > 0) THEN
            CALL lai_empirical (m,nl_soil,rootfr(1:,m),&
                                t_soisno(1:,i),lai(i),sai(i),fveg(i),green(i))

         ENDIF
      ENDDO
#else

! yuan, 08/03/2019: read global LAI/SAI data
      CALL julian2monthday (year, jday, month, mday)
      CALL LAI_readin_nc (year,month,dir_srfdata,nam_srfdata)

#ifdef URBAN_MODEL
      ! 读取城市LAI/SAI数据
      CALL UrbanLAI_readin_nc (year,month,dir_srfdata,nam_urbdata)
#endif

#endif

! ..............................................................................
! 4.5 initialize time-varying variables, as subgrid vectors of length [numpatch]
! ..............................................................................
      DO i = 1, numpatch
         z_soisno (1:nl_soil ,i) = z_soi (1:nl_soil)
         dz_soisno(1:nl_soil ,i) = dz_soi(1:nl_soil)
      ENDDO

    ! ------------------------------------------
    ! PLEASE
    ! PLEASE UPDATE
    ! PLEASE UPDATE when have the observed lake status
      t_lake      (:,:) = 285.
      lake_icefrac(:,:) = 0.
    ! ------------------------------------------

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE (i, m, u) &
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
          )

#if(defined URBAN_MODEL)
      IF (m == URBAN) THEN

         u = patch2urb(i)
         !print *, "patch:", i, "urban:", u, "coszen:", coszen(i)
         lwsun         (u) = 0.   !net longwave radiation of sunlit wall
         lwsha         (u) = 0.   !net longwave radiation of shaded wall
         lgimp         (u) = 0.   !net longwave radiation of impervious road
         lgper         (u) = 0.   !net longwave radiation of pervious road
         lveg          (u) = 0.   !net longwave radiation of vegetation [W/m2]

         t_roofsno   (:,u) = 283. !temperatures of roof layers
         t_wallsun   (:,u) = 283. !temperatures of sunlit wall layers
         t_wallsha   (:,u) = 283. !temperatures of shaded wall layers
         t_gimpsno   (:,u) = 283. !temperatures of impervious road layers
         t_gpersno   (:,u) = 283. !soil temperature [K]
         t_lakesno   (:,u) = 283. !lake soil temperature [K]

         wice_roofsno(:,u) = 0.   !ice lens [kg/m2]
         wice_gimpsno(:,u) = 0.   !ice lens [kg/m2]
         wice_gpersno(:,u) = 0.   !ice lens [kg/m2]
         wice_lakesno(:,u) = 0.   !ice lens [kg/m2]
         wliq_roofsno(:,u) = 0.   !liqui water [kg/m2]
         wliq_gimpsno(:,u) = 0.   !liqui water [kg/m2]
         wliq_gpersno(:,u) = wliq_soisno(:,i) !liqui water [kg/m2]
         wliq_lakesno(:,u) = wliq_soisno(:,i) !liqui water [kg/m2]

         wliq_soisno (:,i) = 0.
         wliq_soisno (:,i) = wliq_roofsno(:,u)*froof(u)
         wliq_soisno (:,i) = wliq_soisno(:,i) + wliq_gpersno(:,u)*(1-froof(u))*fgper(u)
         wliq_soisno (:,i) = wliq_soisno(:,i) + wliq_gimpsno(:,u)*(1-froof(u))*(1-fgper(u))

         snowdp_roof   (u) = 0.   !snow depth [m]
         snowdp_gimp   (u) = 0.   !snow depth [m]
         snowdp_gper   (u) = 0.   !snow depth [m]

         z_sno_roof  (:,u) = 0.   !node depth of roof [m]
         z_sno_gimp  (:,u) = 0.   !node depth of impervious [m]
         z_sno_gper  (:,u) = 0.   !node depth pervious [m]
         z_sno_lake  (:,u) = 0.   !node depth lake [m]

         dz_sno_roof (:,u) = 0.   !interface depth of roof [m]
         dz_sno_gimp (:,u) = 0.   !interface depth of impervious [m]
         dz_sno_gper (:,u) = 0.   !interface depth pervious [m]
         dz_sno_lake (:,u) = 0.   !interface depth lake [m]

         t_room        (u) = 283. !temperature of inner building [K]
         troof_inner   (u) = 283. !temperature of inner roof [K]
         twsun_inner   (u) = 283. !temperature of inner sunlit wall [K]
         twsha_inner   (u) = 283. !temperature of inner shaded wall [K]
         Fhac          (u) = 0.   !sensible flux from heat or cool AC [W/m2]
         Fwst          (u) = 0.   !waste heat flux from heat or cool AC [W/m2]
         Fach          (u) = 0.   !flux from inner and outter air exchange [W/m2]

         CALL UrbanIniTimeVar(i,froof(u),fgper(u),flake(u),hwr(u),hroof(u),&
            alb_roof(:,:,u),alb_wall(:,:,u),alb_gimp(:,:,u),alb_gper(:,:,u),&
            rho(:,:,m),tau(:,:,m),fveg(i),htop(i),hbot(i),lai(i),sai(i),coszen(i),&
            fsno_roof(u),fsno_gimp(u),fsno_gper(u),fsno_lake(u),&
            scv_roof(u),scv_gimp(u),scv_gper(u),scv_lake(u),&
            sag_roof(u),sag_gimp(u),sag_gper(u),sag_lake(u),t_lake(1,i),&
            fwsun(u),dfwsun(u),alb(:,:,i),ssun(:,:,i),ssha(:,:,i),sroof(:,:,u),&
            swsun(:,:,u),swsha(:,:,u),sgimp(:,:,u),sgper(:,:,u),slake(:,:,u))
      ENDIF
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

! ...............................................................
! 4.6 Write out the model variables for restart run [histTimeVar]
! ...............................................................

      write (6,*)
      write (6,*) ('Successfully to Initialize the Land Time-Vraying Variables')

      finfolist = trim(dir_restart)//'clmini.infolist'//'.lc'//trim(cyear)
      open(100,file=trim(finfolist),form='formatted')
      write(100,*) 'numpatch  = ', numpatch   !1.1
      write(100,*) 'numpft    = ', numpft     !1.2
      write(100,*) 'numpc     = ', numpc      !1.3
      write(100,*) 'numurban  = ', numurban   !1.4
      IF ( greenwich ) THEN
      write(100,*) 'greenwich =       .true.' !2
      ELSE
      write(100,*) 'greenwich =      .false.' !2
      ENDIF
      write(100,*) 's_year    = ', year       !3
      write(100,*) 's_month   = ', month      !4
      write(100,*) 's_day     = ', mday       !5
      write(100,*) 's_seconds = ', msec       !6
      write(100,*) '/'
      close(100)

      deallocate (z_soisno )
      deallocate (dz_soisno)

#if(defined URBAN_MODEL)
      deallocate (urbanpct)
#endif

      go to 1000
100   print 101,lndname
101   format(' error occured on file: ',a50)
1000  continue


END SUBROUTINE LuLccInitialize
! --------------------------------------------------
! EOP
