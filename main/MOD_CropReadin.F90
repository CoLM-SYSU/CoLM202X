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
      ! ===========================================================
      ! ! DESCRIPTION:
      ! Read in crop planting date from data, and fertilization from data.
      ! Save these data in patch vector.
      !
      ! Original: Shupeng Zhang, Zhongwang Wei, and Xingjie Lu, 2022
      ! ===========================================================

      use MOD_Precision
      use MOD_Namelist
      use MOD_SPMD_Task
      use MOD_LandPatch
      USE MOD_NetCDFSerial
      USE MOD_NetCDFBlock
      USE MOD_Mapping_Grid2Pset
      use MOD_Vars_TimeInvariants
      use MOD_Vars_TimeVariables

      USE MOD_Vars_Global
      USE MOD_LandPFT
      USE MOD_Vars_PFTimeVariables
      USE MOD_RangeCheck

      IMPLICIT NONE

      CHARACTER(len=256) :: file_crop
      TYPE(grid_type) :: grid_crop
      TYPE(block_data_real8_2d)    :: f_xy_crop
      type(mapping_grid2pset_type) :: mg2patch_crop
      type(mapping_grid2pset_type) :: mg2pft_crop

      real(r8),allocatable :: pdrice2_tmp   (:)
      real(r8),allocatable :: plantdate_tmp (:)
      real(r8),allocatable :: fertnitro_tmp (:)

      ! Local variables
      REAL(r8), allocatable :: lat(:), lon(:)
      integer :: cft, npatch, ipft
      CHARACTER(LEN=2) :: cx
      ! READ in crops
      
      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'

      CALL ncio_read_bcast_serial (file_crop, 'lat', lat)
      CALL ncio_read_bcast_serial (file_crop, 'lon', lon)

      CALL grid_crop%define_by_center (lat, lon)

      call mg2patch_crop%build (grid_crop, landpatch)
      call mg2pft_crop%build   (grid_crop, landpft)

      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      
      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_crop, f_xy_crop)
      ENDIF
      
      IF (p_is_worker) THEN
         IF (numpatch > 0)  allocate(pdrice2_tmp   (numpatch))
         IF (numpft   > 0)  allocate(plantdate_tmp (numpft))
         IF (numpft   > 0)  allocate(fertnitro_tmp (numpft))
      ENDIF

      ! (1) Read in plant date for rice2.
      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'
      IF (p_is_io) THEN
         CALL ncio_read_block (file_crop, 'pdrice2', grid_crop, f_xy_crop)
      ENDIF
      
      call mg2patch_crop%map_aweighted (f_xy_crop, pdrice2_tmp)

      IF (p_is_worker) then
         DO npatch = 1, numpatch
            pdrice2 (npatch) = int(pdrice2_tmp (npatch))
         ENDDO
      ENDIF

#ifdef RangeCheck
      CALL check_vector_data ('plant date value for rice2 ', pdrice2)
#endif

      ! (2) Read in plant date.
      IF (p_is_worker) THEN
         plantdate_p(:) = -99999999._r8
      ENDIF

      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'
      DO cft = 15, 78
         write(cx, '(i2.2)') cft
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_crop, 'PLANTDATE_CFT_'//trim(cx), grid_crop, 1, f_xy_crop)
         ENDIF
      
         call mg2pft_crop%map_aweighted (f_xy_crop, plantdate_tmp)
      
         if (p_is_worker) then
            do ipft = 1, numpft
               IF(landpft%settyp(ipft) .eq. cft)THEN
                  plantdate_p(ipft) = plantdate_tmp(ipft)
                  if(plantdate_p(ipft) <= 0._r8) then
                     plantdate_p(ipft) = -99999999._r8
                  end if
               endif
            end do
         ENDIF
      ENDDO

#ifdef RangeCheck
      CALL check_vector_data ('plantdate_pfts value ', plantdate_p)
#endif

      IF (p_is_worker) THEN
         fertnitro_p(:) = -99999999._r8
      ENDIF

      file_crop = trim(DEF_dir_runtime) // '/crop/fertnitro_fillcoast.nc'
      DO cft = 15, 78
         write(cx, '(i2.2)') cft
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_crop, 'CONST_FERTNITRO_CFT_'//trim(cx), grid_crop, 1, f_xy_crop)
         ENDIF
      
         call mg2pft_crop%map_aweighted (f_xy_crop, fertnitro_tmp)
      
         if (p_is_worker) then
            do ipft = 1, numpft
               IF(landpft%settyp(ipft) .eq. cft)THEN
                  fertnitro_p(ipft) = fertnitro_tmp(ipft)
                  if(fertnitro_p(ipft) <= 0._r8) then
                     fertnitro_p(ipft) = -99999999._r8
                  end if
               endif
            end do
         ENDIF
      ENDDO
      
#ifdef RangeCheck
      CALL check_vector_data ('fert nitro value ', fertnitro_p)
#endif

      IF (allocated (pdrice2_tmp  )) deallocate(pdrice2_tmp  )
      IF (allocated (plantdate_tmp)) deallocate(plantdate_tmp)
      IF (allocated (fertnitro_tmp)) deallocate(fertnitro_tmp)

   END SUBROUTINE CROP_readin

END MODULE MOD_CropReadin
#endif
