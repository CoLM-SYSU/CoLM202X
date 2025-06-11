#include <define.h>

MODULE MOD_Lulcc_TransferTraceReadin

!------------------------------------------------------------------------
!
! !DESCRIPTION:
!  The transfer matrix and patch tracing vector were created using the
!  land cover type data of the adjacent two years. Based on next year's
!  patch, the pixels within the patch and last years' land cover type of
!  these pixels were obtained. Then the percent of source land cover
!  type of each patch was derived.
!
!  Created by Wanyi Lin, Shupeng Zhang and Hua Yuan, 07/2023
!
! !HISTORY:
!  05/2025, Wanyi Lin and Hua Yuan: change to reading the lulcc transfer
!           tracing vector from data files which are produced in the making
!           surface data stage. See mksrfdata/MOD_Lulcc_TransferTrace.F90.
!
!------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   IMPLICIT NONE
   SAVE

   !frac of source patches in a patch
   real(r8), allocatable, dimension(:,:) :: lccpct_patches(:,:)

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_LulccTransferTrace
   PUBLIC :: deallocate_LulccTransferTrace
   PUBLIC :: LulccTransferTraceReadin

   ! PRIVATE MEMBER FUNCTIONS:


CONTAINS


   SUBROUTINE allocate_LulccTransferTrace
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
        allocate (lccpct_patches (numpatch, 0:nlc))
        lccpct_patches (:,:) = 0
      ENDIF

   END SUBROUTINE allocate_LulccTransferTrace


   SUBROUTINE LulccTransferTraceReadin (lc_year)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
   USE MOD_AggregationRequestData
   USE MOD_Mesh
   USE MOD_MeshFilter
   USE MOD_LandElm
   USE MOD_DataType
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_5x5DataReadin
   USE MOD_RegionClip
   USE MOD_Utils
#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_Vars_TimeInvariants

   IMPLICIT NONE

   integer, intent(in) :: lc_year

   character(len=256) :: dir_landdata, lndname, thisyr
   character(len=4)   :: c2
   integer :: ilc
   real(r8), allocatable :: tmpvec(:)


      IF (p_is_worker) THEN
         IF (allocated(tmpvec)) deallocate (tmpvec)
         allocate (tmpvec(numpatch))
      ENDIF

      write(thisyr,'(i4.4)') lc_year
      dir_landdata = DEF_dir_landdata
      DO ilc = 0, N_land_classification
         write(c2, '(i2.2)') ilc
         lndname = trim(dir_landdata)//'/lulcc/'//trim(thisyr)//&
            '/lccpct_patches_lc'//trim(c2)//'.nc'
         CALL ncio_read_vector(lndname, 'lccpct_patches', landpatch, tmpvec)

         IF (p_is_worker) THEN
            lccpct_patches(:,ilc) = tmpvec
         ENDIF
      ENDDO


   END SUBROUTINE LulccTransferTraceReadin


   SUBROUTINE deallocate_LulccTransferTrace
   ! --------------------------------------------------
   ! Deallocates memory for Lulcc time invariant variables
   ! --------------------------------------------------
   USE MOD_SPMD_Task

      IF (p_is_worker) THEN
         IF (allocated(lccpct_patches)) deallocate (lccpct_patches)
      ENDIF

   END SUBROUTINE deallocate_LulccTransferTrace

END MODULE MOD_Lulcc_TransferTraceReadin
! ---------- EOP ------------
