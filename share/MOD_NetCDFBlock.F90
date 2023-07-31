#include <define.h>

MODULE MOD_NetCDFBlock

   !----------------------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !    High-level Subroutines to read and write variables in files with netCDF format.
   !
   !    CoLM read and write netCDF files mainly in three ways:
   !    1. Serial: read and write data by a single process;
   !    2. Vector: 1) read vector data by IO and scatter from IO to workers
   !               2) gather from workers to IO and write vectors by IO
   !               Notice: each file contains vector data in one block.
   !    3. Block : read blocked data by IO
   !               Notice: input file is a single file.
   !    
   !    This module contains subroutines of "3. Block".
   !
   ! Created by Shupeng Zhang, May 2023
   !----------------------------------------------------------------------------------

   USE netcdf 
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   ! PUBLIC subroutines
   interface ncio_read_block
      MODULE procedure ncio_read_block_int32_2d 
      MODULE procedure ncio_read_block_real8_2d 
      MODULE procedure ncio_read_block_real8_3d 
   END interface ncio_read_block
   
   interface ncio_read_block_time
      MODULE procedure ncio_read_block_int32_2d_time 
      MODULE procedure ncio_read_block_real8_2d_time 
   END interface ncio_read_block_time

   PUBLIC :: ncio_read_site_time

CONTAINS

   ! ----
   SUBROUTINE ncio_read_block_int32_2d (filename, dataname, grid, rdata)
     
      USE netcdf
      USE MOD_Block
      USE MOD_Grid
      USE MOD_DataType
      USE MOD_SPMD_Task
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      TYPE (grid_type),  intent(in) :: grid

      TYPE (block_data_int32_2d), intent(inout) :: rdata

      ! Local variables
      INTEGER :: iblk, jblk, ndims(2), start2(2), count2(2), start_mem
      INTEGER :: ncid, varid
      INTEGER :: iblkme

      IF (p_is_io) THEN
         
         CALL check_ncfile_exist (filename)
         CALL nccheck (nf90_open(trim(filename), NF90_NOWRITE, ncid) )
         CALL nccheck (nf90_inq_varid(ncid, trim(dataname), varid) )
         
         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

            ndims = (/grid%xcnt(iblk), grid%ycnt(jblk)/)
            IF (any(ndims == 0)) cycle

            start2 = (/grid%xdsp(iblk)+1, grid%ydsp(jblk)+1/)
            count2(1) = min(grid%xcnt(iblk), grid%nlon-grid%xdsp(iblk))
            count2(2) = grid%ycnt(jblk)

            IF (count2(1) == grid%xcnt(iblk)) THEN
               CALL nccheck (nf90_get_var(ncid, varid, rdata%blk(iblk,jblk)%val, &
                  start2, count2) )
            ELSE
               CALL nccheck (nf90_get_var(ncid, varid, &
                  rdata%blk(iblk,jblk)%val(1:count2(1),:), start2, count2) )

               start2(1) = 1
               start_mem = count2(1) + 1
               count2(1) = grid%xdsp(iblk) + grid%xcnt(iblk) - grid%nlon
               CALL nccheck (nf90_get_var(ncid, varid, &
                  rdata%blk(iblk,jblk)%val(start_mem:ndims(1),:), start2, count2) )
            ENDIF

         ENDDO
                
         CALL nccheck( nf90_close(ncid) )

      ENDIF

   END SUBROUTINE ncio_read_block_int32_2d

   ! ----
   SUBROUTINE ncio_read_block_real8_2d (filename, dataname, grid, rdata)
     
      USE netcdf
      USE MOD_Block
      USE MOD_Grid
      USE MOD_DataType
      USE MOD_SPMD_Task
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      TYPE (grid_type),  intent(in) :: grid

      TYPE (block_data_real8_2d), intent(inout) :: rdata

      ! Local variables
      INTEGER :: iblk, jblk, ndims(2), start2(2), count2(2), start_mem
      INTEGER :: ncid, varid
      INTEGER :: iblkme

      IF (p_is_io) THEN
                     
         CALL check_ncfile_exist (filename)
         CALL nccheck (nf90_open(trim(filename), NF90_NOWRITE, ncid) )
         CALL nccheck (nf90_inq_varid(ncid, trim(dataname), varid) )
         
         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

            ndims = (/grid%xcnt(iblk), grid%ycnt(jblk)/)
            IF (any(ndims == 0)) cycle

            start2 = (/grid%xdsp(iblk)+1, grid%ydsp(jblk)+1/)
            count2(1) = min(grid%xcnt(iblk), grid%nlon-grid%xdsp(iblk))
            count2(2) = grid%ycnt(jblk)

            IF (count2(1) == grid%xcnt(iblk)) THEN
               CALL nccheck (nf90_get_var(ncid, varid, rdata%blk(iblk,jblk)%val, &
                  start2, count2) )
            ELSE
               CALL nccheck (nf90_get_var(ncid, varid, &
                  rdata%blk(iblk,jblk)%val(1:count2(1),:), start2, count2) )

               start2(1) = 1
               start_mem = count2(1) + 1
               count2(1) = grid%xdsp(iblk) + grid%xcnt(iblk) - grid%nlon
               CALL nccheck (nf90_get_var(ncid, varid, &
                  rdata%blk(iblk,jblk)%val(start_mem:ndims(1),:), start2, count2) )
            ENDIF

         ENDDO
                
         CALL nccheck( nf90_close(ncid) )

      ENDIF

   END SUBROUTINE ncio_read_block_real8_2d

   ! ----
   SUBROUTINE ncio_read_block_real8_3d (filename, dataname, grid, ndim1, rdata)
     
      USE netcdf
      USE MOD_Block
      USE MOD_Grid
      USE MOD_DataType
      USE MOD_SPMD_Task
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      TYPE (grid_type),  intent(in) :: grid
      INTEGER, intent(in) :: ndim1

      TYPE (block_data_real8_3d), intent(inout) :: rdata
      INTEGER :: ncid, varid

      ! Local variables
      INTEGER :: iblk, jblk, ndims(3), start3(3), count3(3), start_mem
      INTEGER :: iblkme

      IF (p_is_io) THEN
                     
         CALL check_ncfile_exist (filename)
         CALL nccheck (nf90_open(trim(filename), NF90_NOWRITE, ncid) )
         CALL nccheck (nf90_inq_varid(ncid, trim(dataname), varid) )
         
         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

            ndims = (/ndim1, grid%xcnt(iblk), grid%ycnt(jblk)/)
            IF (any(ndims == 0)) cycle

            start3 = (/1, grid%xdsp(iblk)+1, grid%ydsp(jblk)+1/)
            count3(1) = ndim1
            count3(2) = min(grid%xcnt(iblk), grid%nlon-grid%xdsp(iblk))
            count3(3) = grid%ycnt(jblk)

            IF (count3(2) == grid%xcnt(iblk)) THEN
               CALL nccheck (nf90_get_var(ncid, varid, rdata%blk(iblk,jblk)%val, &
                  start3, count3) )
            ELSE
               CALL nccheck (nf90_get_var(ncid, varid, &
                  rdata%blk(iblk,jblk)%val(:,1:count3(1),:), start3, count3) )

               start3(2) = 1
               start_mem = count3(2) + 1
               count3(2) = grid%xdsp(iblk) + grid%xcnt(iblk) - grid%nlon
               CALL nccheck (nf90_get_var(ncid, varid, &
                  rdata%blk(iblk,jblk)%val(:,start_mem:ndims(2),:), start3, count3) )
            ENDIF

         ENDDO
                
         CALL nccheck( nf90_close(ncid) )

      ENDIF

   END SUBROUTINE ncio_read_block_real8_3d

   ! ----
   SUBROUTINE ncio_read_block_int32_2d_time (filename, dataname, grid, itime, rdata)
     
      USE netcdf
      USE MOD_Block
      USE MOD_Grid
      USE MOD_DataType
      USE MOD_SPMD_Task
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      TYPE (grid_type),  intent(in) :: grid
      INTEGER, intent(in) :: itime

      TYPE (block_data_int32_2d), intent(inout) :: rdata

      ! Local variables
      INTEGER :: iblk, jblk, ndims(2), start3(3), count3(3), start_mem
      INTEGER :: ncid, varid
      INTEGER :: iblkme

      IF (p_is_io) THEN
                     
         CALL check_ncfile_exist (filename)
         CALL nccheck (nf90_open(trim(filename), NF90_NOWRITE, ncid) )
         CALL nccheck (nf90_inq_varid(ncid, trim(dataname), varid) )
         
         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

            ndims = (/grid%xcnt(iblk), grid%ycnt(jblk)/)
            IF (any(ndims == 0)) cycle

            start3 = (/grid%xdsp(iblk)+1, grid%ydsp(jblk)+1, itime/)
            count3(1) = min(grid%xcnt(iblk), grid%nlon-grid%xdsp(iblk))
            count3(2) = grid%ycnt(jblk)
            count3(3) = 1

            IF (count3(1) == grid%xcnt(iblk)) THEN
               CALL nccheck (nf90_get_var(ncid, varid, rdata%blk(iblk,jblk)%val, &
                  start3, count3) )
            ELSE
               CALL nccheck (nf90_get_var(ncid, varid, &
                  rdata%blk(iblk,jblk)%val(1:count3(1),:), start3, count3) )

               start3(1) = 1
               start_mem = count3(1) + 1
               count3(1) = grid%xdsp(iblk) + grid%xcnt(iblk) - grid%nlon
               CALL nccheck (nf90_get_var(ncid, varid, &
                  rdata%blk(iblk,jblk)%val(start_mem:ndims(1),:), start3, count3) )
            ENDIF

         ENDDO
                
         CALL nccheck( nf90_close(ncid) )

      ENDIF

   END SUBROUTINE ncio_read_block_int32_2d_time

   ! ----
   SUBROUTINE ncio_read_block_real8_2d_time (filename, dataname, grid, itime, rdata)
     
      USE netcdf
      USE MOD_Block
      USE MOD_Grid
      USE MOD_DataType
      USE MOD_SPMD_Task
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      TYPE (grid_type),  intent(in) :: grid
      INTEGER, intent(in) :: itime

      TYPE (block_data_real8_2d), intent(inout) :: rdata

      ! Local variables
      INTEGER :: iblk, jblk, ndims(2), start3(3), count3(3), start_mem
      INTEGER :: ncid, varid
      INTEGER :: iblkme

      IF (p_is_io) THEN
                     
         CALL check_ncfile_exist (filename)
         CALL nccheck (nf90_open(trim(filename), NF90_NOWRITE, ncid) )
         CALL nccheck (nf90_inq_varid(ncid, trim(dataname), varid) )
         
         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

            ndims = (/grid%xcnt(iblk), grid%ycnt(jblk)/)
            IF (any(ndims == 0)) cycle

            start3 = (/grid%xdsp(iblk)+1, grid%ydsp(jblk)+1, itime/)
            count3(1) = min(grid%xcnt(iblk), grid%nlon-grid%xdsp(iblk))
            count3(2) = grid%ycnt(jblk)
            count3(3) = 1
            IF (count3(1) == grid%xcnt(iblk)) THEN
               CALL nccheck (nf90_get_var(ncid, varid, rdata%blk(iblk,jblk)%val, &
                  start3, count3) )
            ELSE
               CALL nccheck (nf90_get_var(ncid, varid, &
                  rdata%blk(iblk,jblk)%val(1:count3(1),:), start3, count3) )

               start3(1) = 1
               start_mem = count3(1) + 1
               count3(1) = grid%xdsp(iblk) + grid%xcnt(iblk) - grid%nlon
               CALL nccheck (nf90_get_var(ncid, varid, &
                  rdata%blk(iblk,jblk)%val(start_mem:ndims(1),:), start3, count3) )
            ENDIF

         ENDDO
                
         CALL nccheck( nf90_close(ncid) )

      ENDIF

   END SUBROUTINE ncio_read_block_real8_2d_time

   ! ----
   SUBROUTINE ncio_read_site_time (filename, dataname, itime, rdata)
     
      USE netcdf
      USE MOD_Block
      USE MOD_DataType
      USE MOD_SPMD_Task
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      INTEGER, intent(in) :: itime

      TYPE (block_data_real8_2d), intent(inout) :: rdata

      ! Local variables
      INTEGER :: start3(3), count3(3)
      INTEGER :: ncid, varid

      IF (p_is_io) THEN
                     
         CALL check_ncfile_exist (filename)
         CALL nccheck (nf90_open(trim(filename), NF90_NOWRITE, ncid) ,trace=trim(filename))
         CALL nccheck (nf90_inq_varid(ncid, trim(dataname), varid) ,trace=trim(dataname))
         
         start3 = (/1, 1, itime/)
         count3 = (/1, 1, 1/)
         CALL nccheck (nf90_get_var(ncid, varid, &
            rdata%blk(gblock%xblkme(1),gblock%yblkme(1))%val, start3, count3) )
                
         CALL nccheck( nf90_close(ncid) )

      ENDIF

   END SUBROUTINE ncio_read_site_time

END MODULE MOD_NetCDFBlock
