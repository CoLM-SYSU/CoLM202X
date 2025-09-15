#include <define.h>

#ifdef CROP
MODULE MOD_CropReadin

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: CROP_readin

CONTAINS

   SUBROUTINE CROP_readin ()
!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Read in crop planting date from data, and fertilization from data.
!  Save these data in patch vector.
!
!  Original: Shupeng Zhang, Zhongwang Wei, and Xingjie Lu, 2022
!
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_LandPatch
   USE MOD_NetCDFSerial
   USE MOD_NetCDFBlock
   USE MOD_SpatialMapping
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables

   USE MOD_Vars_Global
   USE MOD_LandPFT
   USE MOD_Vars_PFTimeVariables
   USE MOD_RangeCheck
   USE MOD_Block

   IMPLICIT NONE

   character(len=256) :: file_crop
   type(grid_type)    :: grid_crop
   type(block_data_real8_2d)  :: f_xy_crop
   type(spatial_mapping_type) :: mg2patch_crop
   type(spatial_mapping_type) :: mg2pft_crop
   character(len=256) :: file_irrig
   type(grid_type)    :: grid_irrig
   type(block_data_int32_2d)  :: f_xy_irrig
   type(spatial_mapping_type) :: mg2pft_irrig

   character(len=256) :: file_irrigalloc
   type(grid_type)    :: grid_irrigalloc
   type(block_data_real8_2d)  :: f_xy_irrigalloc
   type(spatial_mapping_type) :: mg2p_irrigalloc

   character(len=256) :: file_fert
   type(grid_type)    :: grid_fert
   type(block_data_real8_2d) :: f_xy_fert
   type(spatial_mapping_type) :: mg2pft_fert

   real(r8),allocatable :: pdrice2_tmp      (:)
   real(r8),allocatable :: plantdate_tmp    (:)
   real(r8),allocatable :: fertnitro_tmp    (:)
   integer ,allocatable :: irrig_method_tmp (:)

   real(r8),allocatable :: fertilizer_tmp   (:)
   real(r8),allocatable :: manure_tmp       (:)

   ! Local variables
   real(r8), allocatable :: lat(:), lon(:)
   real(r8) :: missing_value
   integer  :: cft, npatch, ipft
   character(len=2) :: cx
   integer  :: iblkme, iblk, jblk
   integer  :: maxvalue, minvalue

      ! READ in crops

      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'

      CALL ncio_read_bcast_serial (file_crop, 'lat', lat)
      CALL ncio_read_bcast_serial (file_crop, 'lon', lon)

      CALL grid_crop%define_by_center (lat, lon)

      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_crop, f_xy_crop)
      ENDIF

      ! missing value
      IF (p_is_master) THEN
         CALL ncio_get_attr (file_crop, 'pdrice2', 'missing_value', missing_value)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (missing_value, 1, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif
      IF (p_is_io) THEN
         CALL ncio_read_block (file_crop, 'pdrice2', grid_crop, f_xy_crop)
      ENDIF

      CALL mg2patch_crop%build_arealweighted (grid_crop, landpatch)
      CALL mg2patch_crop%set_missing_value   (f_xy_crop, missing_value)

      CALL mg2pft_crop%build_arealweighted   (grid_crop, landpft)
      CALL mg2pft_crop%set_missing_value     (f_xy_crop, missing_value)

      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      IF (p_is_worker) THEN
         IF (numpatch > 0)  allocate(pdrice2_tmp    (numpatch))
         IF (numpft   > 0)  allocate(plantdate_tmp    (numpft))
         IF (numpft   > 0)  allocate(fertnitro_tmp    (numpft))
         IF (numpft   > 0)  allocate(irrig_method_tmp (numpft))
      ENDIF

      ! (1) Read in plant date for rice2.
      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'
      IF (p_is_io) THEN
         CALL ncio_read_block (file_crop, 'pdrice2', grid_crop, f_xy_crop)
      ENDIF

      CALL mg2patch_crop%grid2pset (f_xy_crop, pdrice2_tmp)

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            IF (pdrice2_tmp(npatch) /= spval) THEN
               pdrice2 (npatch) = int(pdrice2_tmp (npatch))
            ELSE
               pdrice2 (npatch) = 0
            ENDIF
         ENDDO
      ENDIF

#ifdef RangeCheck
      CALL check_vector_data ('plant date value for rice2 ', pdrice2)
#endif

      ! (2) Read in plant date.
      IF (p_is_worker) THEN
         IF (numpft > 0) plantdate_p(:) = -99999999._r8
      ENDIF

      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'
      DO cft = 15, 78
         write(cx, '(i2.2)') cft
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_crop, &
               'PLANTDATE_CFT_'//trim(cx), grid_crop, 1, f_xy_crop)
         ENDIF

         CALL mg2pft_crop%grid2pset (f_xy_crop, plantdate_tmp)

         IF (p_is_worker) THEN
            DO ipft = 1, numpft
               IF(landpft%settyp(ipft) .eq. cft)THEN
                  plantdate_p(ipft) = plantdate_tmp(ipft)
                  IF(plantdate_p(ipft) <= 0._r8) THEN
                     plantdate_p(ipft) = -99999999._r8
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

#ifdef RangeCheck
      CALL check_vector_data ('plantdate_pfts value ', plantdate_p)
#endif

      ! (3) Read in fertlization
      IF (DEF_FERT_SOURCE == 1) THEN
         IF (p_is_worker) THEN
            IF (numpft > 0) fertnitro_p(:) = -99999999._r8
         ENDIF

         file_fert = trim(DEF_dir_runtime) // '/crop/fertnitro_fillcoast.nc'
         DO cft = 15, 78
            write(cx, '(i2.2)') cft
            IF (p_is_io) THEN
               CALL ncio_read_block_time (file_fert, &
                  'CONST_FERTNITRO_CFT_'//trim(cx), grid_crop, 1, f_xy_crop)
            ENDIF

            CALL mg2pft_crop%grid2pset (f_xy_crop, fertnitro_tmp)

            IF (p_is_worker) THEN
               DO ipft = 1, numpft
                  IF(landpft%settyp(ipft) .eq. cft)THEN
                     fertnitro_p(ipft) = fertnitro_tmp(ipft)
                     IF(fertnitro_p(ipft) <= 0._r8) THEN
                        fertnitro_p(ipft) = 0._r8
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

      ELSEIF (DEF_FERT_SOURCE == 2)THEN
         IF (p_is_worker) THEN
            IF (numpft > 0) THEN
               allocate(fertilizer_tmp(numpft))
               allocate(manure_tmp(numpft))
            ENDIF
         ENDIF

         file_fert = trim(DEF_dir_runtime) // '/crop/fertilizer_2015soc.nc'

         CALL ncio_read_bcast_serial (file_fert, 'lat', lat)
         CALL ncio_read_bcast_serial (file_fert, 'lon', lon)
         CALL grid_fert%define_by_center (lat, lon)
         
         IF (p_is_io) THEN
            CALL allocate_block_data (grid_fert, f_xy_fert)
         ENDIF

         CALL mg2pft_fert%build_arealweighted (grid_fert, landpft)

         IF (allocated(lon)) deallocate(lon)
         IF (allocated(lat)) deallocate(lat)

         IF (p_is_worker) THEN
            fertnitro_p(:) = -99999999._r8
            manunitro_p(:) = -99999999._r8
         ENDIF
         
         ! read manure
         IF (p_is_io) THEN
            CALL ncio_read_block (file_fert, 'manure', grid_fert, f_xy_fert)
         ENDIF

         CALL mg2pft_fert%grid2pset (f_xy_fert, manure_tmp)

         IF (p_is_worker) THEN
            DO ipft = 1, numpft
               IF (landpft%settyp(ipft) .ge. 15) THEN
                  manunitro_p(ipft) = manure_tmp(ipft)
                  IF (manunitro_p(ipft) < 0._r8) THEN
                     manunitro_p(ipft) = 0._r8
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

         ! read fertilizer
         DO cft = 1, 64
            IF (p_is_io) THEN
               CALL ncio_read_block_time (file_fert, 'fertilizer', grid_fert, cft, f_xy_fert)
            ENDIF

            CALL mg2pft_fert%grid2pset (f_xy_fert, fertilizer_tmp)

            IF (p_is_worker) THEN
               DO ipft = 1, numpft
                  IF (landpft%settyp(ipft) .eq. cft+14) THEN
                     fertnitro_p(ipft) = fertilizer_tmp(ipft)
                     IF (fertnitro_p(ipft) < 0._r8) THEN
                        fertnitro_p(ipft) = 0._r8
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         IF (allocated (fertilizer_tmp)) deallocate(fertilizer_tmp)
         IF (allocated (manure_tmp))     deallocate(manure_tmp)
      ENDIF

#ifdef RangeCheck
      CALL check_vector_data ('fert nitro value ', fertnitro_p)
#endif

      ! (4) Read in irrigation method
      !file_irrig = trim(DEF_dir_runtime) // '/crop/surfdata_irrigation_method.nc'
      file_irrig = trim(DEF_dir_runtime) // '/crop/surfdata_irrigation_method_96x144.nc'

      CALL ncio_read_bcast_serial (file_irrig, 'lat', lat)
      CALL ncio_read_bcast_serial (file_irrig, 'lon', lon)

      CALL grid_irrig%define_by_center (lat, lon)

      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_irrig, f_xy_irrig)
      ENDIF

      CALL mg2pft_irrig%build_arealweighted (grid_irrig, landpft)

      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      IF (p_is_worker) THEN
         IF (numpft > 0) irrig_method_p(:) = -99999999
      ENDIF

      DO cft = 1, N_CFT
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_irrig, 'irrigation_method', grid_irrig, cft, f_xy_irrig)
         ENDIF

         CALL mg2pft_irrig%grid2pset_dominant (f_xy_irrig, irrig_method_tmp)

         IF (p_is_worker) THEN
            DO ipft = 1, numpft

               IF(landpft%settyp(ipft) .eq. cft + 14)THEN
                  irrig_method_p(ipft) = irrig_method_tmp(ipft)
                  IF(irrig_method_p(ipft) < 0) THEN
                     irrig_method_p(ipft) = -99999999
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

#ifdef RangeCheck
      CALL check_vector_data ('irrigation method ', irrig_method_p)
#endif

      IF (allocated (pdrice2_tmp  ))    deallocate(pdrice2_tmp  )
      IF (allocated (plantdate_tmp))    deallocate(plantdate_tmp)
      IF (allocated (fertnitro_tmp))    deallocate(fertnitro_tmp)
      IF (allocated (irrig_method_tmp)) deallocate(irrig_method_tmp)

      IF (DEF_IRRIGATION_ALLOCATION == 3) THEN 
         ! (5) Read in irrigation allocated to groundwater
         file_irrigalloc = trim(DEF_dir_runtime) // '/crop/surfdata_irrigation_allocation.nc'

         CALL ncio_read_bcast_serial (file_irrigalloc, 'lat', lat)
         CALL ncio_read_bcast_serial (file_irrigalloc, 'lon', lon)

         CALL grid_irrigalloc%define_by_center (lat, lon)

         IF (p_is_io) THEN
            CALL allocate_block_data (grid_irrigalloc, f_xy_irrigalloc)
         ENDIF

         call mg2p_irrigalloc%build_arealweighted (grid_irrigalloc, landpatch)

         IF (allocated(lon)) deallocate(lon)
         IF (allocated(lat)) deallocate(lat)

         IF (p_is_worker) THEN
            irrig_gw_alloc(:) = 0._r8
         ENDIF

         IF (p_is_io) THEN
            CALL ncio_read_block (file_irrigalloc, 'irrig_gw_alloc', grid_irrigalloc, f_xy_irrigalloc)
         ENDIF

         call mg2p_irrigalloc%grid2pset (f_xy_irrigalloc, irrig_gw_alloc)

#ifdef RangeCheck
         CALL check_vector_data ('irrigation goundwater allocation ', irrig_gw_alloc)
#endif

         ! (6) Read in irrigation allocated to surfacewater
         IF (p_is_worker) THEN
            irrig_sw_alloc(:) = 0._r8
         ENDIF

         IF (p_is_io) THEN
            CALL ncio_read_block (file_irrigalloc, 'irrig_sw_alloc', grid_irrigalloc, f_xy_irrigalloc)
         ENDIF

         call mg2p_irrigalloc%grid2pset (f_xy_irrigalloc, irrig_sw_alloc)

#ifdef RangeCheck
         CALL check_vector_data ('irrigation surfacewater allocation ', irrig_sw_alloc)
#endif
      ENDIF

   END SUBROUTINE CROP_readin

END MODULE MOD_CropReadin
#endif
