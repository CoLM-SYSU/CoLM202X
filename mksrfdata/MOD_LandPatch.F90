#include <define.h>

MODULE MOD_LandPatch

   !------------------------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !    Build pixelset "landpatch".
   !
   !    In CoLM, the global/regional area is divided into a hierarchical structure:
   !    1. If GRIDBASED or UNSTRUCTURED is defined, it is
   !       ELEMENT >>> PATCH
   !    2. If CATCHMENT is defined, it is
   !       ELEMENT >>> HRU >>> PATCH
   !    If Plant Function Type classification is used, PATCH is further divided into PFT.
   !    If Plant Community classification is used,     PATCH is further divided into PC.
   !
   !    "landpatch" refers to pixelset PATCH.
   !
   ! Created by Shupeng Zhang, May 2023
   !    porting codes from Hua Yuan's OpenMP version to MPI parallel version.
   !------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_Vars_Global
   USE MOD_Const_LC
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif
   IMPLICIT NONE

   ! ---- Instance ----
   INTEGER :: numpatch
   TYPE(grid_type)     :: gpatch
   TYPE(pixelset_type) :: landpatch

#if (defined CROP)
   TYPE(grid_type) :: gcrop
   REAL(r8), allocatable :: pctcrop   (:)
   INTEGER,  allocatable :: cropclass (:)
#endif

   TYPE(subset_type)   :: elm_patch
   TYPE(superset_type) :: patch2elm

#ifdef CATCHMENT
   TYPE(subset_type)   :: hru_patch
   TYPE(superset_type) :: patch2hru
#endif


CONTAINS

   ! -------------------------------
   SUBROUTINE landpatch_build (lc_year)

      USE MOD_Precision
      USE MOD_SPMD_Task
      USE MOD_Utils
      USE MOD_Grid
      USE MOD_DataType
      USE MOD_Mesh
      USE MOD_LandElm
#ifdef CATCHMENT
      USE MOD_LandHRU
#endif
      USE MOD_Namelist
      USE MOD_NetCDFBlock
#if (defined CROP)
      USE MOD_PixelsetShadow
#endif
      USE MOD_AggregationRequestData

      IMPLICIT NONE

      INTEGER, intent(in) :: lc_year
      ! Local Variables
      CHARACTER(len=256) :: file_patch
      CHARACTER(len=255) :: cyear
      TYPE (block_data_int32_2d) :: patchdata
      INTEGER :: iloc, npxl, ipxl, numset
      INTEGER :: ie, iset, ipxstt, ipxend
      INTEGER, allocatable :: types(:), order(:), ibuff(:)
      INTEGER, allocatable :: eindex_tmp(:), settyp_tmp(:), ipxstt_tmp(:), ipxend_tmp(:), ielm_tmp(:)
      LOGICAL, allocatable :: msk(:)
      INTEGER :: npatch_glb
#if (defined CROP)
      TYPE(block_data_real8_3d) :: cropdata
      INTEGER :: cropfilter(1)
#endif
      INTEGER :: dominant_type
      INTEGER, allocatable :: npxl_types (:)

      write(cyear,'(i4.4)') lc_year
      IF (p_is_master) THEN
         write(*,'(A)') 'Making land patches :'
      ENDIF

#if (defined SinglePoint && defined LULC_IGBP_PFT && defined CROP)
      IF ((SITE_landtype == CROPLAND) .and. (USE_SITE_pctcrop)) THEN

         numpatch = count(SITE_pctcrop > 0.)

         allocate (pctcrop  (numpatch))
         allocate (cropclass(numpatch))
         cropclass = pack(SITE_croptyp, SITE_pctcrop > 0.)
         pctcrop   = pack(SITE_pctcrop, SITE_pctcrop > 0.)

         pctcrop = pctcrop / sum(pctcrop)

         allocate (landpatch%eindex (numpatch))
         allocate (landpatch%ipxstt (numpatch))
         allocate (landpatch%ipxend (numpatch))
         allocate (landpatch%settyp (numpatch))
         allocate (landpatch%ielm   (numpatch))

         landpatch%eindex(:) = 1
         landpatch%ielm  (:) = 1
         landpatch%ipxstt(:) = 1
         landpatch%ipxend(:) = 1
         landpatch%settyp(:) = CROPLAND

         landpatch%nset = numpatch
         CALL landpatch%set_vecgs

         RETURN
      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifndef SinglePoint
      IF (p_is_io) THEN
         CALL allocate_block_data (gpatch, patchdata)

#ifndef LULC_USGS
         ! add parameter input for time year
         file_patch = trim(DEF_dir_rawdata)//'landtypes/landtype-igbp-modis-'//trim(cyear)//'.nc'
#else
         !TODO: need usgs land cover TYPE data
         file_patch = trim(DEF_dir_rawdata) //'/landtypes/landtype_usgs_update.nc'
#endif
         CALL ncio_read_block (file_patch, 'landtype', gpatch, patchdata)

#ifdef USEMPI
         CALL aggregation_data_daemon (gpatch, data_i4_2d_in1 = patchdata)
#endif
      ENDIF
#endif

      IF (p_is_worker) THEN

#ifdef CATCHMENT
         numset = numhru
#else
         numset = numelm
#endif

         IF (numset > 0) THEN
            allocate (eindex_tmp (numset*N_land_classification))
            allocate (settyp_tmp (numset*N_land_classification))
            allocate (ipxstt_tmp (numset*N_land_classification))
            allocate (ipxend_tmp (numset*N_land_classification))
            allocate (ielm_tmp   (numset*N_land_classification))
         ENDIF

         numpatch = 0

         DO iset = 1, numset
#ifdef CATCHMENT
            ie     = landhru%ielm  (iset)
            ipxstt = landhru%ipxstt(iset)
            ipxend = landhru%ipxend(iset)
#else
            ie     = landelm%ielm  (iset)
            ipxstt = landelm%ipxstt(iset)
            ipxend = landelm%ipxend(iset)
#endif

            npxl   = ipxend - ipxstt + 1

            allocate (types (ipxstt:ipxend))

#ifndef SinglePoint
#ifdef CATCHMENT
            CALL aggregation_request_data (landhru, iset, gpatch, &
#else
            CALL aggregation_request_data (landelm, iset, gpatch, &
#endif
               data_i4_2d_in1 = patchdata, data_i4_2d_out1 = ibuff)


            types(:) = ibuff
            deallocate (ibuff)
#else
            types(:) = SITE_landtype
#endif

#ifdef CATCHMENT
            IF (landhru%settyp(iset) <= 0) THEN
               types(ipxstt:ipxend) = WATERBODY
            ENDIF
#endif

IF ((DEF_USE_PFT .and. .not.DEF_SOLO_PFT) .or. DEF_FAST_PC) THEN
            ! For classification of plant function types or fast PC,
            ! merge all land types with soil ground
            DO ipxl = ipxstt, ipxend
               IF (types(ipxl) > 0) THEN
                  IF (patchtypes(types(ipxl)) == 0) THEN
#if (defined CROP)
                     !12  Croplands
                     !14  Cropland/Natural Vegetation Mosaics  ?
                     IF (types(ipxl) /= CROPLAND) THEN
                        types(ipxl) = 1
                     ENDIF
#else
                     types(ipxl) = 1
#endif
                  ENDIF
               ENDIF
            ENDDO
ENDIF

            allocate (order (ipxstt:ipxend))
            order = (/ (ipxl, ipxl = ipxstt, ipxend) /)

            CALL quicksort (npxl, types, order)

            mesh(ie)%ilon(ipxstt:ipxend) = mesh(ie)%ilon(order)
            mesh(ie)%ilat(ipxstt:ipxend) = mesh(ie)%ilat(order)

            IF (DEF_USE_DOMINANT_PATCHTYPE) THEN
               allocate (npxl_types (0:maxval(types)))
               npxl_types(:) = 0
               DO ipxl = ipxstt, ipxend
                  npxl_types(types(ipxl)) = npxl_types(types(ipxl)) + 1
               ENDDO

               IF (any(types > 0)) THEN
                  iloc = findloc(types > 0, .true., dim=1) + ipxstt - 1
                  dominant_type = maxloc(npxl_types(1:), dim=1)
                  types(iloc:ipxend) = dominant_type
               ENDIF

               deallocate(npxl_types)
            ENDIF

            DO ipxl = ipxstt, ipxend
               IF (ipxl == ipxstt) THEN
                  numpatch = numpatch + 1
                  eindex_tmp(numpatch) = mesh(ie)%indx
                  settyp_tmp(numpatch) = types(ipxl)
                  ipxstt_tmp(numpatch) = ipxl
                  ielm_tmp  (numpatch) = ie
               ELSEIF (types(ipxl) /= types(ipxl-1)) THEN
                  ipxend_tmp(numpatch) = ipxl - 1

                  numpatch = numpatch + 1
                  eindex_tmp(numpatch) = mesh(ie)%indx
                  settyp_tmp(numpatch) = types(ipxl)
                  ipxstt_tmp(numpatch) = ipxl
                  ielm_tmp  (numpatch) = ie
               ENDIF
            ENDDO
            ipxend_tmp(numpatch) = ipxend

            deallocate (types)
            deallocate (order)

         ENDDO

         IF (numpatch > 0) THEN
            allocate (landpatch%eindex (numpatch))
            allocate (landpatch%settyp (numpatch))
            allocate (landpatch%ipxstt (numpatch))
            allocate (landpatch%ipxend (numpatch))
            allocate (landpatch%ielm   (numpatch))

            landpatch%eindex = eindex_tmp(1:numpatch)
            landpatch%ipxstt = ipxstt_tmp(1:numpatch)
            landpatch%ipxend = ipxend_tmp(1:numpatch)
            landpatch%settyp = settyp_tmp(1:numpatch)
            landpatch%ielm   = ielm_tmp  (1:numpatch)
         ENDIF

         IF (numset > 0) THEN
            deallocate (eindex_tmp)
            deallocate (ipxstt_tmp)
            deallocate (ipxend_tmp)
            deallocate (settyp_tmp)
            deallocate (ielm_tmp  )
         ENDIF

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif

      ENDIF

      landpatch%nset = numpatch

      CALL landpatch%set_vecgs

      IF (DEF_LANDONLY) THEN
         IF ((p_is_worker) .and. (numpatch > 0)) THEN
            allocate(msk(numpatch))
            msk = (landpatch%settyp /= 0)
         ENDIF

         CALL landpatch%pset_pack (msk, numpatch)

         IF (allocated(msk)) deallocate(msk)
      ENDIF

#ifdef URBAN_MODEL
      continue
#else
#if (defined CROP)
      IF (p_is_io) THEN
!         file_patch = trim(DEF_dir_rawdata) // '/global_0.5x0.5.MOD2005_V4.5_CFT_mergetoclmpft.nc'
         file_patch = trim(DEF_dir_rawdata) // '/global_0.5x0.5.MOD2005_V4.5_CFT_lf-merged-20220930.nc'
         CALL allocate_block_data (gcrop, cropdata, N_CFT)
         CALL ncio_read_block (file_patch, 'PCT_CFT', gcrop, N_CFT, cropdata)
      ENDIF

      cropfilter = (/ CROPLAND /)

      CALL pixelsetshadow_build (landpatch, gcrop, cropdata, N_CFT, cropfilter, &
         pctcrop, cropclass)

      numpatch = landpatch%nset
#endif

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numpatch, npatch_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npatch_glb, ' patches.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numpatch, ' patches.'
#endif

#if (defined CROP)
      CALL elm_patch%build (landelm, landpatch, use_frac = .true., shadowfrac = pctcrop)
#else
      CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
#endif

#ifdef CATCHMENT
#if (defined CROP)
      CALL hru_patch%build (landhru, landpatch, use_frac = .true., shadowfrac = pctcrop)
#else
      CALL hru_patch%build (landhru, landpatch, use_frac = .true.)
#endif
#endif

      CALL write_patchfrac (DEF_dir_landdata, lc_year)
#endif
   END SUBROUTINE landpatch_build

   ! -----
   SUBROUTINE write_patchfrac (dir_landdata, lc_year)

      USE MOD_Namelist
      USE MOD_NetCDFVector
      IMPLICIT NONE

      INTEGER, intent(in) :: lc_year
      CHARACTER(LEN=*), intent(in) :: dir_landdata
      CHARACTER(len=256) :: lndname, cyear

      write(cyear,'(i4.4)') lc_year
      CALL system('mkdir -p ' // trim(dir_landdata) // '/landpatch/' // trim(cyear))

      lndname = trim(dir_landdata)//'/landpatch/'//trim(cyear)//'/patchfrac_elm.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'patchfrac_elm', 'patch', landpatch, elm_patch%subfrc, 1)

#ifdef CATCHMENT
      lndname = trim(dir_landdata)//'/landpatch/'//trim(cyear)//'/patchfrac_hru.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'patchfrac_hru', 'patch', landpatch, hru_patch%subfrc, 1)
#endif

   END SUBROUTINE write_patchfrac

END MODULE MOD_LandPatch
