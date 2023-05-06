#include <define.h>

MODULE mod_data_type

   !-----------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !    Definations of data types used in CoLM.
   !   
   ! Created by Shupeng Zhang, May 2023

   USE precision

   ! ---- data types ----
   !-------
   TYPE :: pointer_real8_1d
      REAL(r8), allocatable :: val(:)
   CONTAINS 
      final :: pointer_real8_1d_free_mem
   END TYPE pointer_real8_1d

   !-------
   TYPE :: pointer_int32_1d
      INTEGER, allocatable :: val(:)
   CONTAINS 
      final :: pointer_int32_1d_free_mem
   END TYPE pointer_int32_1d

   !-------
   TYPE :: pointer_int64_1d
      INTEGER*8, allocatable :: val(:)
   CONTAINS 
      final :: pointer_int64_1d_free_mem
   END TYPE pointer_int64_1d

   !-------
   TYPE :: pointer_int32_2d
      INTEGER, allocatable :: val (:,:)
   CONTAINS 
      final :: pointer_int32_2d_free_mem
   END TYPE pointer_int32_2d

   TYPE :: block_data_int32_2d
      TYPE(pointer_int32_2d), allocatable :: blk (:,:)
   CONTAINS
      final :: block_data_int32_2d_free_mem
   END TYPE block_data_int32_2d

   !-------
   TYPE :: pointer_real8_2d
      REAL(r8), allocatable :: val (:,:)
   CONTAINS 
      final :: pointer_real8_2d_free_mem
   END TYPE pointer_real8_2d

   TYPE :: block_data_real8_2d
      TYPE(pointer_real8_2d), allocatable :: blk (:,:)
   CONTAINS
      final :: block_data_real8_2d_free_mem
   END TYPE block_data_real8_2d

   !-------
   TYPE :: pointer_real8_3d
      REAL(r8), allocatable :: val (:,:,:)
   CONTAINS 
      final :: pointer_real8_3d_free_mem
   END TYPE pointer_real8_3d

   TYPE :: block_data_real8_3d
      INTEGER :: lb1, ub1
      TYPE(pointer_real8_3d), allocatable :: blk (:,:)
   CONTAINS
      final :: block_data_real8_3d_free_mem
   END TYPE block_data_real8_3d

   !-------
   TYPE :: pointer_real8_4d
      REAL(r8), allocatable :: val (:,:,:,:)
   CONTAINS 
      final :: pointer_real8_4d_free_mem
   END TYPE pointer_real8_4d

   TYPE :: block_data_real8_4d
      INTEGER :: lb1, ub1, lb2, ub2
      TYPE(pointer_real8_4d), allocatable :: blk (:,:)
   CONTAINS
      final :: block_data_real8_4d_free_mem
   END TYPE block_data_real8_4d

   ! ---- PUBLIC subroutines ----
   !------
   interface allocate_block_data
      MODULE procedure allocate_block_data_int32_2d
      MODULE procedure allocate_block_data_real8_2d
      MODULE procedure allocate_block_data_real8_3d
      MODULE procedure allocate_block_data_real8_4d
   END interface allocate_block_data

   !------
   interface flush_block_data
      MODULE procedure flush_block_data_int32_2d
      MODULE procedure flush_block_data_real8_2d
      MODULE procedure flush_block_data_real8_3d
      MODULE procedure flush_block_data_real8_4d
   END interface flush_block_data

   !-----
   PUBLIC :: block_data_linear_transform
   PUBLIC :: block_data_copy
   PUBLIC :: block_data_linear_interp

CONTAINS

   !------------------
   SUBROUTINE pointer_real8_1d_free_mem (this)

      IMPLICIT NONE

      TYPE(pointer_real8_1d) :: this

      IF (allocated(this%val)) THEN
         deallocate(this%val)
      ENDIF

   END SUBROUTINE pointer_real8_1d_free_mem 

   !------------------
   SUBROUTINE pointer_int32_1d_free_mem (this)

      IMPLICIT NONE

      TYPE(pointer_int32_1d) :: this

      IF (allocated(this%val)) THEN
         deallocate(this%val)
      ENDIF

   END SUBROUTINE pointer_int32_1d_free_mem 

   !------------------
   SUBROUTINE pointer_int64_1d_free_mem (this)

      IMPLICIT NONE

      TYPE(pointer_int64_1d) :: this

      IF (allocated(this%val)) THEN
         deallocate(this%val)
      ENDIF

   END SUBROUTINE pointer_int64_1d_free_mem 

   !------------------
   SUBROUTINE pointer_int32_2d_free_mem (this)

      IMPLICIT NONE

      TYPE(pointer_int32_2d) :: this

      IF (allocated(this%val)) THEN
         deallocate(this%val)
      ENDIF

   END SUBROUTINE pointer_int32_2d_free_mem 

   !------------------
   SUBROUTINE allocate_block_data_int32_2d (grid, gdata)

      USE mod_grid
      USE mod_block
      USE spmd_task
      IMPLICIT NONE

      TYPE(grid_type), intent(in) :: grid
      TYPE(block_data_int32_2d), intent(out) :: gdata

      ! Local variables
      INTEGER :: iblkme, iblk, jblk

      allocate (gdata%blk (gblock%nxblk,gblock%nyblk))

      DO iblkme = 1, gblock%nblkme 
         iblk = gblock%xblkme(iblkme)
         jblk = gblock%yblkme(iblkme)
         allocate (gdata%blk(iblk,jblk)%val (grid%xcnt(iblk), grid%ycnt(jblk)))
      ENDDO

   END SUBROUTINE allocate_block_data_int32_2d

   !------------------
   SUBROUTINE block_data_int32_2d_free_mem (this)

      USE mod_block
      IMPLICIT NONE

      TYPE(block_data_int32_2d) :: this

      ! Local variables
      INTEGER :: iblk, jblk

      IF (allocated (this%blk)) THEN
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (allocated (this%blk(iblk,jblk)%val)) THEN
                  deallocate (this%blk(iblk,jblk)%val)
               ENDIF
            ENDDO
         ENDDO

         deallocate (this%blk)
      ENDIF

   END SUBROUTINE block_data_int32_2d_free_mem 

   !------------------
   SUBROUTINE pointer_real8_2d_free_mem (this)

      IMPLICIT NONE

      TYPE(pointer_real8_2d) :: this

      IF (allocated(this%val)) THEN
         deallocate(this%val)
      ENDIF

   END SUBROUTINE pointer_real8_2d_free_mem 

   !------------------
   SUBROUTINE allocate_block_data_real8_2d (grid, gdata)

      USE mod_grid
      USE mod_block
      USE spmd_task
      IMPLICIT NONE

      TYPE(grid_type), intent(in) :: grid
      TYPE(block_data_real8_2d), intent(out) :: gdata

      ! Local variables
      INTEGER :: iblkme, iblk, jblk

      allocate (gdata%blk (gblock%nxblk,gblock%nyblk))

      DO iblkme = 1, gblock%nblkme 
         iblk = gblock%xblkme(iblkme)
         jblk = gblock%yblkme(iblkme)
         allocate (gdata%blk(iblk,jblk)%val (grid%xcnt(iblk), grid%ycnt(jblk)))
      ENDDO

   END SUBROUTINE allocate_block_data_real8_2d

   !------------------
   SUBROUTINE block_data_real8_2d_free_mem (this)

      USE mod_block
      IMPLICIT NONE

      TYPE(block_data_real8_2d) :: this

      ! Local variables
      INTEGER :: iblk, jblk
      
      IF (allocated (this%blk)) THEN
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (allocated (this%blk(iblk,jblk)%val)) THEN
                  deallocate (this%blk(iblk,jblk)%val)
               ENDIF
            ENDDO
         ENDDO

         deallocate (this%blk)
      ENDIF
      
   END SUBROUTINE block_data_real8_2d_free_mem 

   !------------------
   SUBROUTINE pointer_real8_3d_free_mem (this)

      IMPLICIT NONE

      TYPE(pointer_real8_3d) :: this

      IF (allocated(this%val)) THEN
         deallocate(this%val)
      ENDIF

   END SUBROUTINE pointer_real8_3d_free_mem 

   !------------------
   SUBROUTINE allocate_block_data_real8_3d (grid, gdata, ndim1, lb1)

      USE mod_grid
      USE mod_block
      USE spmd_task
      IMPLICIT NONE

      TYPE(grid_type), intent(in) :: grid
      TYPE(block_data_real8_3d), intent(out) :: gdata
      INTEGER, intent(in) :: ndim1
      INTEGER, intent(in), optional :: lb1

      ! Local variables
      INTEGER :: iblkme, iblk, jblk

      allocate (gdata%blk (gblock%nxblk,gblock%nyblk))
      
      IF (present(lb1)) THEN
         gdata%lb1 = lb1
      ELSE
         gdata%lb1 = 1
      ENDIF

      gdata%ub1 = gdata%lb1-1+ndim1

      DO iblkme = 1, gblock%nblkme 
         iblk = gblock%xblkme(iblkme)
         jblk = gblock%yblkme(iblkme)
         allocate (gdata%blk(iblk,jblk)%val (gdata%lb1:gdata%ub1, grid%xcnt(iblk), grid%ycnt(jblk)))
      ENDDO

   END SUBROUTINE allocate_block_data_real8_3d

   !------------------
   SUBROUTINE block_data_real8_3d_free_mem (this)

      USE mod_block
      IMPLICIT NONE

      TYPE(block_data_real8_3d) :: this

      ! Local variables
      INTEGER :: iblk, jblk

      IF (allocated (this%blk)) THEN
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (allocated (this%blk(iblk,jblk)%val)) THEN
                  deallocate (this%blk(iblk,jblk)%val)
               ENDIF
            ENDDO
         ENDDO

         deallocate (this%blk)
      ENDIF

   END SUBROUTINE block_data_real8_3d_free_mem 

   !------------------
   SUBROUTINE pointer_real8_4d_free_mem (this)

      IMPLICIT NONE

      TYPE(pointer_real8_4d) :: this

      IF (allocated(this%val)) THEN
         deallocate(this%val)
      ENDIF

   END SUBROUTINE pointer_real8_4d_free_mem 

   !------------------
   SUBROUTINE allocate_block_data_real8_4d (grid, gdata, ndim1, ndim2, lb1, lb2)

      USE mod_grid
      USE mod_block
      USE spmd_task
      IMPLICIT NONE

      TYPE(grid_type), intent(in) :: grid
      TYPE(block_data_real8_4d), intent(out) :: gdata
      INTEGER, intent(in) :: ndim1, ndim2
      INTEGER, intent(in), optional :: lb1, lb2

      ! Local variables
      INTEGER :: iblkme, iblk, jblk

      allocate (gdata%blk (gblock%nxblk,gblock%nyblk))

      IF (present(lb1)) THEN
         gdata%lb1 = lb1
      ELSE
         gdata%lb1 = 1
      ENDIF

      gdata%ub1 = gdata%lb1-1+ndim1

      IF (present(lb2)) THEN
         gdata%lb2 = lb2
      ELSE
         gdata%lb2 = 1
      ENDIF

      gdata%ub2 = gdata%lb2-1+ndim2

      DO iblkme = 1, gblock%nblkme 
         iblk = gblock%xblkme(iblkme)
         jblk = gblock%yblkme(iblkme)
         allocate (gdata%blk(iblk,jblk)%val ( &
            gdata%lb1:gdata%ub1, gdata%lb2:gdata%ub2, grid%xcnt(iblk), grid%ycnt(jblk)))
      ENDDO

   END SUBROUTINE allocate_block_data_real8_4d

   !------------------
   SUBROUTINE block_data_real8_4d_free_mem (this)

      USE mod_block
      IMPLICIT NONE

      TYPE(block_data_real8_4d) :: this

      ! Local variables
      INTEGER :: iblk, jblk

      IF (allocated (this%blk)) THEN
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (allocated (this%blk(iblk,jblk)%val)) THEN
                  deallocate (this%blk(iblk,jblk)%val)
               ENDIF
            ENDDO
         ENDDO

         deallocate (this%blk)
      ENDIF

   END SUBROUTINE block_data_real8_4d_free_mem 

   !------------------
   SUBROUTINE flush_block_data_real8_2d (gdata, spval)

      USE precision
      USE mod_block
      USE spmd_task
      IMPLICIT NONE

      TYPE(block_data_real8_2d), intent(inout) :: gdata
      REAL(r8), intent(in) :: spval

      ! Local variables
      INTEGER :: iblkme, iblk, jblk

      DO iblkme = 1, gblock%nblkme 
         iblk = gblock%xblkme(iblkme)
         jblk = gblock%yblkme(iblkme)
         gdata%blk(iblk,jblk)%val = spval
      ENDDO

   END SUBROUTINE flush_block_data_real8_2d

   !------------------
   SUBROUTINE flush_block_data_int32_2d (gdata, spval)

      USE precision
      USE mod_block
      USE spmd_task
      IMPLICIT NONE

      TYPE(block_data_int32_2d), intent(inout) :: gdata
      INTEGER, intent(in) :: spval

      ! Local variables
      INTEGER :: iblkme, iblk, jblk

      DO iblkme = 1, gblock%nblkme 
         iblk = gblock%xblkme(iblkme)
         jblk = gblock%yblkme(iblkme)
         gdata%blk(iblk,jblk)%val = spval
      ENDDO

   END SUBROUTINE flush_block_data_int32_2d

   !------------------
   SUBROUTINE flush_block_data_real8_3d (gdata, spval)

      USE precision
      USE mod_block
      USE spmd_task
      IMPLICIT NONE

      TYPE(block_data_real8_3d), intent(inout) :: gdata
      REAL(r8), intent(in) :: spval

      ! Local variables
      INTEGER :: iblkme, iblk, jblk

      DO iblkme = 1, gblock%nblkme 
         iblk = gblock%xblkme(iblkme)
         jblk = gblock%yblkme(iblkme)
         gdata%blk(iblk,jblk)%val = spval
      ENDDO

   END SUBROUTINE flush_block_data_real8_3d

   !------------------
   SUBROUTINE flush_block_data_real8_4d (gdata, spval)

      USE precision
      USE mod_block
      USE spmd_task
      IMPLICIT NONE

      TYPE(block_data_real8_4d), intent(inout) :: gdata
      REAL(r8), intent(in) :: spval

      ! Local variables
      INTEGER :: iblkme, iblk, jblk

      DO iblkme = 1, gblock%nblkme 
         iblk = gblock%xblkme(iblkme)
         jblk = gblock%yblkme(iblkme)
         gdata%blk(iblk,jblk)%val = spval
      ENDDO

   END SUBROUTINE flush_block_data_real8_4d

   !------------------
   SUBROUTINE block_data_linear_transform (gdata, scl, dsp)

      USE precision
      USE mod_block
      USE spmd_task
      IMPLICIT NONE

      TYPE(block_data_real8_2d), intent(inout) :: gdata
      REAL(r8), intent(in), optional :: scl
      REAL(r8), intent(in), optional :: dsp

      ! Local variables
      INTEGER :: iblkme, iblk, jblk

      IF (present(scl)) THEN
         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
            gdata%blk(iblk,jblk)%val = gdata%blk(iblk,jblk)%val * scl
         ENDDO
      ENDIF

      IF (present(dsp)) THEN
         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
            gdata%blk(iblk,jblk)%val = gdata%blk(iblk,jblk)%val + dsp
         ENDDO
      ENDIF

   END SUBROUTINE block_data_linear_transform

   !------------------
   SUBROUTINE block_data_copy (gdata_from, gdata_to, sca)

      USE precision
      USE mod_block
      USE spmd_task
      IMPLICIT NONE

      TYPE(block_data_real8_2d), intent(in)    :: gdata_from
      TYPE(block_data_real8_2d), intent(inout) :: gdata_to
      REAL(r8), intent(in), optional :: sca

      ! Local variables
      INTEGER :: iblkme, iblk, jblk

      DO iblkme = 1, gblock%nblkme 
         iblk = gblock%xblkme(iblkme)
         jblk = gblock%yblkme(iblkme)
         IF (present(sca)) THEN
            gdata_to%blk(iblk,jblk)%val = gdata_from%blk(iblk,jblk)%val * sca
         ELSE
            gdata_to%blk(iblk,jblk)%val = gdata_from%blk(iblk,jblk)%val
         ENDIF
      ENDDO

   END SUBROUTINE block_data_copy

   !------------------
   SUBROUTINE block_data_linear_interp ( &
         gdata_from1, alp1, gdata_from2, alp2, gdata_to)

      USE precision
      USE mod_block
      USE spmd_task
      IMPLICIT NONE

      TYPE(block_data_real8_2d), intent(in)    :: gdata_from1, gdata_from2
      REAL(r8), intent(in) :: alp1, alp2 
      TYPE(block_data_real8_2d), intent(inout) :: gdata_to

      ! Local variables
      INTEGER :: iblkme, iblk, jblk

      DO iblkme = 1, gblock%nblkme 
         iblk = gblock%xblkme(iblkme)
         jblk = gblock%yblkme(iblkme)
         gdata_to%blk(iblk,jblk)%val = &
            gdata_from1%blk(iblk,jblk)%val * alp1 &
            + gdata_from2%blk(iblk,jblk)%val * alp2
      ENDDO

   END SUBROUTINE block_data_linear_interp

      
END MODULE mod_data_type
