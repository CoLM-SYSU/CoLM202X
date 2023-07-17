#include <define.h>

MODULE MOD_Lulcc_TMatrix
! -------------------------------
! Created by Hua Yuan, 04/2022
! may be renamed as MOD_Lulcc_Ptrace
! -------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   IMPLICIT NONE
   SAVE
! -----------------------------------------------------------------
   ! TODO: need coding below...
   real(r8), allocatable :: lccpct_patches(:,:) ! Percent area of source patches 

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_LulccTMatrix
   PUBLIC :: deallocate_LulccTMatrix
   PUBLIC :: READ_LulccTMatrix

   ! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_LulccTMatrix
   ! --------------------------------------------------------------------
   ! Allocates memory for Lulcc time invariant variables
   ! --------------------------------------------------------------------

      USE MOD_Precision
      USE MOD_Vars_Global
      USE MOD_LandPatch
      use MOD_SPMD_Task

      IMPLICIT NONE
      !TODO: need coding below...
      INTEGER :: nlc = N_land_classification

      IF (p_is_worker) THEN
         allocate (lccpct_patches (numpatch, nlc))
         lccpct_patches (:,:) = 0
      ENDIF
      
   END SUBROUTINE allocate_LulccTMatrix

   SUBROUTINE READ_LulccTMatrix (lc_year)

      USE MOD_Precision
      USE MOD_Namelist
      USE MOD_SPMD_Task
      USE MOD_Grid
      USE MOD_LandPatch
      USE MOD_NetCDFVector
      USE MOD_NetCDFBlock
      USE MOD_5x5DataReadin
      USE MOD_Namelist, only : DEF_dir_rawdata
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

      INTEGER, intent(in) :: lc_year

      ! local variables:
      ! ---------------------------------------------------------------
      CHARACTER(len=256) :: dir_5x5, suffix, cyear
      INTEGER :: ipatch,ipxl,ipxstt, ipxend

      TYPE (block_data_int32_2d) :: lcdatafr ! land cover data of last year
      REAL(r8), allocatable :: area_one(:)    , areabuff(:)
      INTEGER,  allocatable :: lcdatafr_one(:), ibuff(:)

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

      !   -----------------------------------------------------------------
      !   extract the land cover type of pixels of last year for each patch
      !   -----------------------------------------------------------------

      IF (p_is_worker) THEN
   
         DO ipatch = 1, numpatch

            CALL aggregation_request_data (landpatch, ipatch, gpatch, area = area_one, &
                  data_i4_2d_in1 = lcdatafr, data_i4_2d_out1 = lcdatafr_one)
            
            ipxstt = landpatch%ipxstt(ipatch)
            ipxend = landpatch%ipxend(ipatch)
            
            if (allocated(ibuff)) deallocate(ibuff)
            allocate(ibuff(ipxstt:ipxend))
            ibuff(:) = lcdatafr_one(:)
            
            if (allocated(areabuff)) deallocate(areabuff)
            allocate(areabuff(ipxstt:ipxend))
            areabuff(:) = area_one(:)

            DO ipxl = ipxstt, ipxend
               lccpct_patches(ipatch, ibuff(ipxl)) = lccpct_patches(ipatch, ibuff(ipxl)) + areabuff(ipxl) / sum(areabuff)
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

   END SUBROUTINE READ_LulccTMatrix

   SUBROUTINE deallocate_LulccTMatrix
      ! --------------------------------------------------
      ! Deallocates memory for Lulcc time invariant variables
      ! --------------------------------------------------
      use MOD_SPMD_Task
      !TODO: need coding below...
      IF (p_is_worker) THEN
         IF (allocated(lccpct_patches)) deallocate (lccpct_patches)
      ENDIF

   END SUBROUTINE deallocate_LulccTMatrix

END MODULE MOD_Lulcc_TMatrix
! ---------- EOP ------------
