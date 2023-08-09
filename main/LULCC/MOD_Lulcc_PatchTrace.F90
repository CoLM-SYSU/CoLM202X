#include <define.h>

MODULE MOD_Lulcc_PatchTrace
! -------------------------------
! Created by Wanyi Lin and Hua Yuan, 07/2023
! may be renamed as MOD_Lulcc_PatchTrace
! -------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   IMPLICIT NONE
   SAVE
! -----------------------------------------------------------------

   real(r8), allocatable :: lccpct_patches(:,:) ! Percent area of source patches

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_LulccPatchTrace
   PUBLIC :: deallocate_LulccPatchTrace
   PUBLIC :: READ_LulccPatchTrace

   ! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_LulccPatchTrace
   ! --------------------------------------------------------------------
   ! Allocates memory for Lulcc time invariant variables
   ! --------------------------------------------------------------------

      USE MOD_Precision
      USE MOD_Vars_Global
      USE MOD_LandPatch
      USE MOD_SPMD_Task

      IMPLICIT NONE

      integer :: nlc = N_land_classification

      IF (p_is_worker) THEN
         allocate (lccpct_patches (numpatch, nlc))
         lccpct_patches (:,:) = 0
      ENDIF

   END SUBROUTINE allocate_LulccPatchTrace

   SUBROUTINE READ_LulccPatchTrace (lc_year)

      USE MOD_Precision
      USE MOD_Namelist
      USE MOD_SPMD_Task
      USE MOD_Grid
      USE MOD_LandPatch
      USE MOD_NetCDFVector
      USE MOD_NetCDFBlock
      USE MOD_5x5DataReadin
      USE MOD_Namelist, only: DEF_dir_rawdata
#ifdef CoLMDEBUG
      USE MOD_RangeCheck
#endif
      USE MOD_AggregationRequestData
      USE MOD_Utils
      USE MOD_Mesh
      USE MOD_MeshFilter
      USE MOD_LandElm
      USE MOD_DataType
      USE MOD_Block
      USE MOD_Pixel
      USE MOD_5x5DataReadin
      USE MOD_RegionClip
      USE MOD_Utils

      IMPLICIT NONE

      integer, intent(in) :: lc_year

      ! local variables:
      ! ---------------------------------------------------------------
      character(len=256) :: dir_5x5, suffix, cyear
      integer :: ipatch,ipxl,ipxstt, ipxend

      type (block_data_int32_2d) :: lcdatafr ! land cover data of last year
      integer,  allocatable :: lcdatafr_one(:), ibuff(:)
      real(r8), allocatable :: area_one(:)    , areabuff(:)
      real(r8) :: sum_areabuff

      write(cyear,'(i4.4)') lc_year-1

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      CALL gpatch%define_by_name ('colm_500m')
      CALL pixel%assimilate_grid (gpatch)
      CALL pixel%map_to_grid (gpatch)

      IF (p_is_io) THEN
         CALL allocate_block_data (gpatch, lcdatafr)
         dir_5x5 = trim(DEF_dir_rawdata) // '/plant_15s'
         suffix  = 'MOD'//trim(cyear)
         CALL read_5x5_data (dir_5x5, suffix, gpatch, 'LC', lcdatafr)

#ifdef USEMPI
         CALL aggregation_data_daemon (gpatch, data_i4_2d_in1 = lcdatafr)
#endif
      ENDIF

      ! -----------------------------------------------------------------
      ! extract the land cover type of pixels of last year for each patch
      ! -----------------------------------------------------------------

      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch

            CALL aggregation_request_data (landpatch, ipatch, gpatch, zip = .true., area = area_one, &
                  data_i4_2d_in1 = lcdatafr, data_i4_2d_out1 = lcdatafr_one)

            ipxstt = landpatch%ipxstt(ipatch)
            ipxend = landpatch%ipxend(ipatch)

            IF (allocated(ibuff)) deallocate(ibuff)
            allocate(ibuff(ipxstt:ipxend))
            ibuff(:) = lcdatafr_one(:)

            IF (allocated(areabuff)) deallocate(areabuff)
            allocate(areabuff(ipxstt:ipxend))
            areabuff(:) = area_one(:)

            sum_areabuff = sum(areabuff)
            DO ipxl = ipxstt, ipxend
               lccpct_patches(ipatch, ibuff(ipxl)) = lccpct_patches(ipatch, ibuff(ipxl)) + areabuff(ipxl) / sum_areabuff
            ENDDO

         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif

         IF (allocated(area_one))    deallocate (area_one)

      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CoLMDEBUG
   CALL check_vector_data ('lccpct_patches ', lccpct_patches)
#endif

   END SUBROUTINE READ_LulccPatchTrace

   SUBROUTINE deallocate_LulccPatchTrace
      ! --------------------------------------------------
      ! Deallocates memory for Lulcc time invariant variables
      ! --------------------------------------------------
      USE MOD_SPMD_Task

      IF (p_is_worker) THEN
         IF (allocated(lccpct_patches)) deallocate (lccpct_patches)
      ENDIF

   END SUBROUTINE deallocate_LulccPatchTrace

END MODULE MOD_Lulcc_PatchTrace
! ---------- EOP ------------
