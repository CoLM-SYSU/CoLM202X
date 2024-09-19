#include <define.h>

!----------------------------------------------------------------------------------
! DESCRIPTION:
!
!    High-level Subroutines to read and write variables in files with netCDF format.
!
!    CoLM read and write netCDF files mainly in three ways:
!    1. Serial: read and write data by a single process;
!    2. Vector: 1) read vector data by IO and scatter from IO to workers
!               2) gather from workers to IO and write vectors by IO
!    3. Block : read blocked data by IO
!               Notice: input file is a single file.
!    
!    This MODULE CONTAINS subroutines of "2. Vector".
!    
!    Two implementations can be used,
!    1) "MOD_NetCDFVectorBlk.F90": 
!       A vector is saved in separated files, each associated with a block. 
!       READ/WRITE are fast in this way and compression can be used.
!       However, there may be too many files, especially when blocks are small. 
!       CHOOSE this implementation by "#undef VectorInOneFile" in include/define.h
!    2) "MOD_NetCDFVectorOne.F90": 
!       A vector is saved in one file. 
!       READ/WRITE may be slow in this way.
!       CHOOSE this implementation by "#define VectorInOneFileS" in include/define.h
!
! Created by Shupeng Zhang, May 2023
!----------------------------------------------------------------------------------

! Put vector in one file.
#ifdef VectorInOneFileS

MODULE MOD_NetCDFVector

   USE netcdf
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_Pixelset
   USE MOD_NetCDFSerial, only : nccheck
   IMPLICIT NONE

   ! PUBLIC subroutines

   PUBLIC :: ncio_create_file_vector 
   PUBLIC :: ncio_define_dimension_vector 

   INTERFACE ncio_read_vector
      MODULE procedure ncio_read_vector_logical_1d 
      MODULE procedure ncio_read_vector_int32_1d 
      MODULE procedure ncio_read_vector_int64_1d 
      MODULE procedure ncio_read_vector_real8_1d 
      MODULE procedure ncio_read_vector_real8_2d 
      MODULE procedure ncio_read_vector_real8_3d 
      MODULE procedure ncio_read_vector_real8_4d 
   END INTERFACE ncio_read_vector

   INTERFACE ncio_write_vector
      MODULE procedure ncio_write_vector_logical_1d
      MODULE procedure ncio_write_vector_int32_1d 
      MODULE procedure ncio_write_vector_int64_1d 
      MODULE procedure ncio_write_vector_real8_1d 
      MODULE procedure ncio_write_vector_real8_2d 
      MODULE procedure ncio_write_vector_real8_3d 
      MODULE procedure ncio_write_vector_real8_4d 
   END INTERFACE ncio_write_vector

CONTAINS

   ! -----
   SUBROUTINE ncio_open_vector (filename, dataname, exit_on_err, ncid, grpid, vecname, noerr)
      
   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   logical, intent(in) :: exit_on_err

   integer, intent(out) :: ncid, grpid
   logical, intent(out) :: noerr
   character(len=*), intent(out) :: vecname
      
      noerr = (nf90_open(trim(filename), NF90_NOWRITE, ncid) == NF90_NOERR)
      IF (.not. noerr) write(*,*) 'Warning: '//trim(filename)//' not found.'

      IF (noerr) noerr = (nf90_inq_ncid(ncid, trim(dataname), grpid) == NF90_NOERR)
      IF (.not. noerr) write(*,*) 'Warning: '//trim(dataname)//' in '//trim(filename)//' not found.'

      IF (noerr) noerr = (nf90_get_att(grpid, NF90_GLOBAL, 'vector_name', vecname) == NF90_NOERR)
      IF (.not. noerr) write(*,*) 'Warning: '//trim(vecname)//' in '//trim(filename)//' not found.'

      IF ((.not. noerr) .and. (exit_on_err)) THEN
         write(*,'(A)') 'Netcdf error in reading ' // trim(dataname) // ' from ' // trim(filename) 
         CALL CoLM_Stop ()
      ENDIF

   END SUBROUTINE ncio_open_vector
   
   !---------------------------------------------------------
   SUBROUTINE ncio_inquire_length_grp (filename, dataname, blkname, length)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: blkname
   integer, intent(out) :: length

   ! Local variables
   integer :: ncid, varid, grpid, ndims
   integer, allocatable :: dimids(:)
   logical :: noerr

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_ncid(ncid, trim(dataname), grpid) )

      noerr = (nf90_inq_varid(grpid, trim(blkname), varid) == NF90_NOERR)

      IF (noerr) THEN
         CALL nccheck( nf90_inquire_variable(grpid, varid, ndims = ndims) )
         allocate (dimids(ndims))
         CALL nccheck( nf90_inquire_variable(grpid, varid, dimids = dimids) )
         CALL nccheck( nf90_inquire_dimension(grpid, dimids(ndims), len = length) )
         deallocate (dimids)
      ELSE
         length = 0
      ENDIF

      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_inquire_length_grp

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_grp_int64_1d (filename, dataname, blkname, rdata)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: blkname
   integer*8, allocatable, intent(out) :: rdata (:)

   ! Local variables
   integer :: ncid, grpid, varid, varlen
   integer, allocatable :: varsize(:)

      CALL ncio_inquire_length_grp (filename, dataname, blkname, varlen)

      IF (varlen > 0) THEN
         allocate (rdata (varlen) )
         CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
         CALL nccheck( nf90_inq_ncid(ncid, trim(dataname), grpid) )
         CALL nccheck( nf90_inq_varid(grpid, trim(blkname), varid) )
         CALL nccheck( nf90_get_var(grpid, varid, rdata) )
         CALL nccheck( nf90_close(ncid) )
      ENDIF

   END SUBROUTINE ncio_read_serial_grp_int64_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_int32_1d ( &
         filename, dataname, pixelset, rdata, defval)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   type(pixelset_type), intent(in) :: pixelset

   integer, allocatable, intent(inout) :: rdata (:)
   integer, intent(in), optional :: defval

   ! Local variables
   integer :: ncid, grpid, varid, iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: blockname, vecname, varname
   integer, allocatable :: sbuff(:), rbuff(:)
   logical :: noerr, ok

      IF (p_is_worker) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (pixelset%nset))
         ENDIF
      ENDIF

      IF (p_is_io) THEN

         CALL ncio_open_vector (filename, dataname, .not. present(defval), &
            ncid, grpid, vecname, noerr)

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))

            IF (noerr) THEN
               CALL get_blockname (iblk, jblk, blockname)
               varname = trim(vecname)//'_'//trim(blockname)
               ok = (nf90_inq_varid(grpid, trim(varname), varid) == NF90_NOERR)
               IF (ok) ok = (nf90_get_var(grpid, varid, sbuff) == NF90_NOERR)
            ELSE
               ok = .false.
            ENDIF

            IF (.not. ok) THEN
               IF (.not. present(defval)) THEN 
                  write(*,'(A)') 'Netcdf error in reading ' &
                     // trim(varname) // ' from ' // trim(filename) 
                  CALL CoLM_Stop ()
               ELSE
                  sbuff = defval
               ENDIF
            ENDIF

#ifdef USEMPI
            CALL mpi_scatterv ( &
               sbuff, pixelset%vecgs%vcnt(:,iblk,jblk), &
               pixelset%vecgs%vdsp(:,iblk,jblk), MPI_INTEGER, &
               MPI_IN_PLACE, 0, MPI_INTEGER, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rdata(istt:iend) = sbuff
#endif

            deallocate (sbuff)

         ENDDO
      
         IF (noerr) CALL nccheck( nf90_close(ncid), trim(filename) // ' close failed' )

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)
                     
            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
            ELSE
               allocate (rbuff(1))
            ENDIF

            CALL mpi_scatterv ( &
               MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER, & ! insignificant on workers
               rbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_INTEGER, &
               p_root, p_comm_group, p_err)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               rdata(istt:iend) = rbuff
            ENDIF

            IF (allocated(rbuff)) deallocate (rbuff)

         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_int32_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_int64_1d ( &
         filename, dataname, pixelset, rdata, defval)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   type(pixelset_type), intent(in) :: pixelset

   integer*8, allocatable, intent(inout) :: rdata (:)
   integer, intent(in), optional :: defval

   ! Local variables
   integer :: ncid, grpid, varid, iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: blockname, varname, vecname
   integer*8, allocatable :: sbuff(:), rbuff(:)
   logical :: noerr, ok

      IF (p_is_worker) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (pixelset%nset))
         ENDIF
      ENDIF

      IF (p_is_io) THEN

         CALL ncio_open_vector (filename, dataname, .not. present(defval), &
            ncid, grpid, vecname, noerr)

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))

            IF (noerr) THEN
               CALL get_blockname (iblk, jblk, blockname)
               varname = trim(vecname)//'_'//trim(blockname)
               ok = (nf90_inq_varid(grpid, trim(varname), varid) == NF90_NOERR)
               IF (ok) ok = (nf90_get_var(grpid, varid, sbuff) == NF90_NOERR)
            ELSE
               ok = .false.
            ENDIF

            IF (.not. ok) THEN
               IF (.not. present(defval)) THEN 
                  write(*,'(A)') 'Netcdf error in reading ' &
                     // trim(varname) // ' from ' // trim(filename) 
                  CALL CoLM_Stop ()
               ELSE
                  sbuff = defval
               ENDIF
            ENDIF

#ifdef USEMPI
            CALL mpi_scatterv ( &
               sbuff, pixelset%vecgs%vcnt(:,iblk,jblk), &
               pixelset%vecgs%vdsp(:,iblk,jblk), MPI_INTEGER8, &
               MPI_IN_PLACE, 0, MPI_INTEGER8, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rdata(istt:iend) = sbuff
#endif

            deallocate (sbuff)

         ENDDO

         IF (noerr) CALL nccheck( nf90_close(ncid), trim(filename) // ' close failed' )

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)
                     
            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
            ELSE
               allocate (rbuff(1))
            ENDIF

            CALL mpi_scatterv ( &
               MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER8, & ! insignificant on workers
               rbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_INTEGER8, &
               p_root, p_comm_group, p_err)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               rdata(istt:iend) = rbuff
            ENDIF

            IF (allocated(rbuff)) deallocate (rbuff)

         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_int64_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_logical_1d (filename, dataname, pixelset, rdata, &
         defval)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   type(pixelset_type), intent(in) :: pixelset

   logical, allocatable, intent(inout) :: rdata (:)
   logical, intent(in), optional :: defval

   ! Local variables
   integer :: ncid, grpid, varid, iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: blockname, varname, vecname
   integer(1), allocatable :: sbuff(:), rbuff(:)
   logical :: noerr, ok

      IF (p_is_worker) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (pixelset%nset))
         ENDIF
      ENDIF

      IF (p_is_io) THEN

         CALL ncio_open_vector (filename, dataname, .not. present(defval), &
            ncid, grpid, vecname, noerr)

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))

            IF (noerr) THEN
               CALL get_blockname (iblk, jblk, blockname)
               varname = trim(vecname)//'_'//trim(blockname)
               ok = (nf90_inq_varid(grpid, trim(varname), varid) == NF90_NOERR)
               IF (ok) ok = (nf90_get_var(grpid, varid, sbuff) == NF90_NOERR)
            ELSE
               ok = .false.
            ENDIF

            IF (.not. ok) THEN
               IF (.not. present(defval)) THEN 
                  write(*,'(A)') 'Netcdf error in reading ' &
                     // trim(varname) // ' from ' // trim(filename) 
                  CALL CoLM_Stop ()
               ELSE
                  IF (defval) THEN
                     sbuff = 1
                  ELSE
                     sbuff = 0
                  ENDIF
               ENDIF
            ENDIF

#ifdef USEMPI
            CALL mpi_scatterv ( &
               sbuff, pixelset%vecgs%vcnt(:,iblk,jblk), &
               pixelset%vecgs%vdsp(:,iblk,jblk), MPI_INTEGER1, &
               MPI_IN_PLACE, 0, MPI_INTEGER1, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rdata(istt:iend) = (sbuff == 1)
#endif

            deallocate (sbuff)

         ENDDO

         IF (noerr) CALL nccheck( nf90_close(ncid), trim(filename) // ' close failed' )
      
      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
            ELSE
               allocate (rbuff(1))
            ENDIF

            CALL mpi_scatterv ( &
               MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER1, & ! insignificant on workers
               rbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_INTEGER1, &
               p_root, p_comm_group, p_err)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               rdata(istt:iend) = (rbuff == 1)
            ENDIF

            IF (allocated(rbuff)) deallocate (rbuff)

         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_logical_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_1d (filename, dataname, pixelset, rdata, &
         defval)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   type(pixelset_type), intent(in) :: pixelset

   real(r8), allocatable, intent(inout) :: rdata (:)
   real(r8), intent(in), optional :: defval

   ! Local variables
   integer :: ncid, grpid, varid, iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: blockname, varname, vecname
   real(r8), allocatable :: sbuff(:), rbuff(:)
   logical :: noerr, ok
         
      IF (p_is_worker) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (pixelset%nset))
         ENDIF
      ENDIF
      
      IF (p_is_io) THEN

         CALL ncio_open_vector (filename, dataname, .not. present(defval), &
            ncid, grpid, vecname, noerr)

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))

            IF (noerr) THEN
               CALL get_blockname (iblk, jblk, blockname)
               varname = trim(vecname)//'_'//trim(blockname)
               ok = (nf90_inq_varid(grpid, trim(varname), varid) == NF90_NOERR)
               IF (ok) ok = (nf90_get_var(grpid, varid, sbuff) == NF90_NOERR)
            ELSE
               ok = .false.
            ENDIF

            IF (.not. ok) THEN
               IF (.not. present(defval)) THEN 
                  write(*,'(A)') 'Netcdf error in reading ' &
                     // trim(varname) // ' from ' // trim(filename) 
                  CALL CoLM_Stop ()
               ELSE
                  sbuff = defval
               ENDIF
            ENDIF

#ifdef USEMPI
            CALL mpi_scatterv ( &
               sbuff, pixelset%vecgs%vcnt(:,iblk,jblk), &
               pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
               MPI_IN_PLACE, 0, MPI_REAL8, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rdata(istt:iend) = sbuff
#endif

            deallocate (sbuff)

         ENDDO

         IF (noerr) CALL nccheck( nf90_close(ncid), trim(filename) // ' close failed' )
      
      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
            ELSE
               allocate (rbuff(1))
            ENDIF

            CALL mpi_scatterv ( &
               MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
               rbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
               p_root, p_comm_group, p_err)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               rdata(istt:iend) = rbuff
            ENDIF

            IF (allocated(rbuff)) deallocate (rbuff)

         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_real8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_2d ( &
         filename, dataname, ndim1, pixelset, rdata, defval)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   integer, intent(in) :: ndim1
   type(pixelset_type), intent(in) :: pixelset

   real(r8), allocatable, intent(inout) :: rdata (:,:)
   real(r8), intent(in), optional :: defval

   ! Local variables
   integer :: ncid, grpid, varid, iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: blockname, varname, vecname
   real(r8), allocatable :: sbuff(:,:), rbuff(:,:)
   logical :: noerr, ok

      IF (p_is_worker) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (ndim1, pixelset%nset))
         ENDIF
      ENDIF

      IF (p_is_io) THEN

         CALL ncio_open_vector (filename, dataname, .not. present(defval), &
            ncid, grpid, vecname, noerr)

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (sbuff (ndim1, pixelset%vecgs%vlen(iblk,jblk)))
            
            IF (noerr) THEN
               CALL get_blockname (iblk, jblk, blockname)
               varname = trim(vecname)//'_'//trim(blockname)
               ok = (nf90_inq_varid(grpid, trim(varname), varid) == NF90_NOERR)
               IF (ok) ok = (nf90_get_var(grpid, varid, sbuff) == NF90_NOERR)
            ELSE
               ok = .false.
            ENDIF

            IF (.not. ok) THEN
               IF (.not. present(defval)) THEN 
                  write(*,'(A)') 'Netcdf error in reading ' &
                     // trim(varname) // ' from ' // trim(filename) 
                  CALL CoLM_Stop ()
               ELSE
                  sbuff = defval
               ENDIF
            ENDIF

#ifdef USEMPI
            CALL mpi_scatterv ( &
               sbuff, ndim1 * pixelset%vecgs%vcnt(:,iblk,jblk), &
               ndim1 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
               MPI_IN_PLACE, 0, MPI_REAL8, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rdata(:,istt:iend) = sbuff
#endif

            deallocate (sbuff)

         ENDDO

         IF (noerr) CALL nccheck( nf90_close(ncid), trim(filename) // ' close failed' )

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (rbuff (ndim1, pixelset%vecgs%vlen(iblk,jblk)))
            ELSE
               allocate (rbuff(1,1))
            ENDIF

            CALL mpi_scatterv ( &
               MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
               rbuff, ndim1 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
               p_root, p_comm_group, p_err)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               rdata(:,istt:iend) = rbuff
            ENDIF

            IF (allocated(rbuff)) deallocate (rbuff)

         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_real8_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_3d ( &
         filename, dataname, ndim1, ndim2, pixelset, rdata, defval)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   integer, intent(in) :: ndim1, ndim2
   type(pixelset_type), intent(in) :: pixelset

   real(r8), allocatable, intent(inout) :: rdata (:,:,:)
   real(r8), intent(in), optional :: defval

   ! Local variables
   integer :: ncid, grpid, varid, iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: blockname, varname, vecname
   real(r8), allocatable :: sbuff(:,:,:), rbuff(:,:,:)
   logical :: noerr, ok

      IF (p_is_worker) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (ndim1,ndim2, pixelset%nset))
         ENDIF
      ENDIF

      IF (p_is_io) THEN

         CALL ncio_open_vector (filename, dataname, .not. present(defval), &
            ncid, grpid, vecname, noerr)

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (sbuff (ndim1,ndim2, pixelset%vecgs%vlen(iblk,jblk)))

            IF (noerr) THEN
               CALL get_blockname (iblk, jblk, blockname)
               varname = trim(vecname)//'_'//trim(blockname)
               ok = (nf90_inq_varid(grpid, trim(varname), varid) == NF90_NOERR)
               IF (ok) ok = (nf90_get_var(grpid, varid, sbuff) == NF90_NOERR)
            ELSE
               ok = .false.
            ENDIF

            IF (.not. ok) THEN
               IF (.not. present(defval)) THEN 
                  write(*,'(A)') 'Netcdf error in reading ' &
                     // trim(varname) // ' from ' // trim(filename) 
                  CALL CoLM_Stop ()
               ELSE
                  sbuff = defval
               ENDIF
            ENDIF

#ifdef USEMPI
            CALL mpi_scatterv ( &
               sbuff, ndim1 * ndim2 * pixelset%vecgs%vcnt(:,iblk,jblk), &
               ndim1 * ndim2 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
               MPI_IN_PLACE, 0, MPI_REAL8, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rdata(:,:,istt:iend) = sbuff
#endif

            deallocate (sbuff)

         ENDDO
         
         IF (noerr) CALL nccheck( nf90_close(ncid), trim(filename) // ' close failed' )

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (rbuff (ndim1,ndim2, pixelset%vecgs%vlen(iblk,jblk)))
            ELSE
               allocate (rbuff(1,1,1))
            ENDIF

            CALL mpi_scatterv ( &
               MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
               rbuff, ndim1 * ndim2 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
               p_root, p_comm_group, p_err)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               rdata(:,:,istt:iend) = rbuff
            ENDIF

            IF (allocated(rbuff)) deallocate (rbuff)

         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_real8_3d
   
   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_4d ( &
         filename, dataname, ndim1, ndim2, ndim3, pixelset, rdata, defval)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   integer, intent(in) :: ndim1, ndim2, ndim3
   type(pixelset_type), intent(in) :: pixelset

   real(r8), allocatable, intent(inout) :: rdata (:,:,:,:)
   real(r8), intent(in), optional :: defval

   ! Local variables
   integer :: ncid, grpid, varid, iblkgrp, iblk, jblk, istt, iend
   character(len=256) :: blockname, varname, vecname
   real(r8), allocatable :: sbuff(:,:,:,:), rbuff(:,:,:,:)
   logical :: noerr, ok

      IF (p_is_worker) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (ndim1,ndim2,ndim3, pixelset%nset))
         ENDIF
      ENDIF

      IF (p_is_io) THEN

         CALL ncio_open_vector (filename, dataname, .not. present(defval), &
            ncid, grpid, vecname, noerr)

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (sbuff (ndim1,ndim2,ndim3, pixelset%vecgs%vlen(iblk,jblk)))

            IF (noerr) THEN
               CALL get_blockname (iblk, jblk, blockname)
               varname = trim(vecname)//'_'//trim(blockname)
               ok = (nf90_inq_varid(grpid, trim(varname), varid) == NF90_NOERR)
               IF (ok) ok = (nf90_get_var(grpid, varid, sbuff) == NF90_NOERR)
            ELSE
               ok = .false.
            ENDIF

            IF (.not. ok) THEN
               IF (.not. present(defval)) THEN 
                  write(*,'(A)') 'Netcdf error in reading ' &
                     // trim(varname) // ' from ' // trim(filename) 
                  CALL CoLM_Stop ()
               ELSE
                  sbuff = defval
               ENDIF
            ENDIF

#ifdef USEMPI
            CALL mpi_scatterv ( &
               sbuff, ndim1 * ndim2 * ndim3 * pixelset%vecgs%vcnt(:,iblk,jblk), &
               ndim1 * ndim2 * ndim3 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
               MPI_IN_PLACE, 0, MPI_REAL8, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rdata(:,:,:,istt:iend) = sbuff
#endif

            deallocate (sbuff)

         ENDDO

         IF (noerr) CALL nccheck( nf90_close(ncid), trim(filename) // ' close failed' )

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         
         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (rbuff (ndim1,ndim2,ndim3, pixelset%vecgs%vlen(iblk,jblk)))
            ELSE
               allocate (rbuff(1,1,1,1))
            ENDIF

            CALL mpi_scatterv ( &
               MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
               rbuff, ndim1 * ndim2 * ndim3 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
               p_root, p_comm_group, p_err)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               rdata(:,:,:,istt:iend) = rbuff
            ENDIF

            IF (allocated(rbuff)) deallocate (rbuff)

         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_real8_4d
   
   !---------------------------------------------------------
   SUBROUTINE ncio_create_file_vector (filename, pixelset)

   USE MOD_NetCDFSerial, only : get_time_now
   IMPLICIT NONE

   character(len=*),    intent(in) :: filename
   type(pixelset_type), intent(in) :: pixelset

   ! Local Variables
   integer :: ncid, mode

      IF (p_is_master) THEN
      
         mode = ior(NF90_NETCDF4,NF90_CLOBBER)
         CALL nccheck( nf90_create(trim(filename), mode, ncid) )

         CALL nccheck (nf90_put_att(ncid, NF90_GLOBAL, 'create_time', get_time_now()))
         CALL nccheck (nf90_enddef (ncid))

         CALL nccheck( nf90_close(ncid) )
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE ncio_create_file_vector

   !---------------------------------------------------------
   SUBROUTINE ncio_define_dimension_vector (filename, pixelset, dimname, dimlen)

      IMPLICIT NONE

      character(len=*),    intent(in) :: filename
      type(pixelset_type), intent(in) :: pixelset
      character(len=*),    intent(in) :: dimname
      integer, optional,   intent(in) :: dimlen

      ! Local variables
      integer :: ncid, dimid, iblkall, iblk, jblk, err
      character(len=8) :: blockname

      IF (p_is_io) THEN
         
         IF (p_iam_io == 0) THEN 

            CALL nccheck( nf90_open (trim(filename), ior(NF90_WRITE,NF90_NETCDF4), ncid) )

            CALL nccheck (nf90_redef(ncid))

            IF (present(dimlen)) THEN
               err = nf90_inq_dimid(ncid, trim(dimname), dimid)
               IF (err /= NF90_NOERR) THEN
                  CALL nccheck( nf90_def_dim(ncid, trim(dimname), dimlen, dimid) )
               ENDIF
            ELSE

               DO iblkall = 1, pixelset%nblkall

                  iblk = pixelset%xblkall(iblkall)
                  jblk = pixelset%yblkall(iblkall)
                  CALL get_blockname (iblk, jblk, blockname)

                  err = nf90_inq_dimid(ncid, trim(dimname)//'_'//trim(blockname), dimid)
                  IF (err /= NF90_NOERR) THEN
                     CALL nccheck( nf90_def_dim(ncid, trim(dimname)//'_'//trim(blockname), &
                        pixelset%vlenall(iblk,jblk), dimid) )
                  ENDIF

               ENDDO
            ENDIF

            CALL nccheck (nf90_enddef(ncid))
            CALL nccheck (nf90_close (ncid))

         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_io, p_err)
#endif
      ENDIF
               
   END SUBROUTINE ncio_define_dimension_vector

   !---------------------------------------------------------
   SUBROUTINE ncio_define_variable_vector ( &
         filename, pixelset, vecname, dataname, datatype, &
         dim1name, dim2name, dim3name, compress)

   IMPLICIT NONE

   character(len=*),    intent(in) :: filename
   integer, intent(in) :: datatype
   type(pixelset_type), intent(in) :: pixelset
   character(len=*),    intent(in) :: vecname
   character(len=*),    intent(in) :: dataname

   character(len=*), optional, intent(in) :: dim1name, dim2name, dim3name
   integer, optional, intent(in) :: compress

   ! Local variables
   integer :: ncid, ndims, idim, iblkall, grpid, varid
   character(len=256) :: varname, blockname
   integer, allocatable :: dimids(:)

      IF (p_iam_io == 0) THEN

         CALL nccheck( nf90_open(trim(filename), ior(NF90_WRITE,NF90_NETCDF4), ncid) )

         ndims = 1
         IF (present(dim1name)) ndims = ndims + 1
         IF (present(dim2name)) ndims = ndims + 1
         IF (present(dim3name)) ndims = ndims + 1

         allocate (dimids(ndims))

         idim = 1
         IF (present(dim1name)) THEN
            CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimids(idim)))
            idim = idim + 1
         ENDIF
         IF (present(dim2name)) THEN
            CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimids(idim)))
            idim = idim + 1
         ENDIF
         IF (present(dim3name)) THEN
            CALL nccheck (nf90_inq_dimid(ncid, trim(dim3name), dimids(idim)))
            idim = idim + 1
         ENDIF

         CALL nccheck( nf90_redef(ncid))

         CALL nccheck( nf90_def_grp(ncid, trim(dataname), grpid) )

         CALL nccheck( nf90_put_att(grpid, NF90_GLOBAL, 'vector_name', trim(vecname)))

         DO iblkall = 1, pixelset%nblkall
            CALL get_blockname (pixelset%xblkall(iblkall), pixelset%yblkall(iblkall), blockname)
            varname = trim(vecname)//'_'//trim(blockname)
            CALL nccheck (nf90_inq_dimid(ncid, trim(varname), dimids(ndims)))
            IF (present(compress)) THEN
               CALL nccheck (nf90_def_var(grpid, trim(varname), datatype, dimids, varid, &
                  deflate_level = compress) )
            ELSE
               CALL nccheck (nf90_def_var(grpid, trim(varname), datatype, dimids, varid) )
            ENDIF
         ENDDO

         CALL nccheck (nf90_enddef(ncid))
         CALL nccheck (nf90_close (ncid))

         deallocate (dimids)

      ENDIF 

#ifdef USEMPI
      CALL mpi_barrier (p_comm_io, p_err)
#endif

   END SUBROUTINE ncio_define_variable_vector

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_int32_1d ( &
         filename, dataname, vecname, pixelset, wdata, compress_level)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: vecname
   type(pixelset_type), intent(in) :: pixelset
   integer, intent(in) :: wdata (:)

   integer, intent(in), optional :: compress_level

   ! Local variables
   integer :: ncid, grpid, varid, iblkall, iblkgrp, iblk, jblk, istt, iend, i, isrc, lock = 0
   character(len=8) :: blockname
   integer, allocatable :: sbuff(:)
   type(pointer_int32_1d), allocatable :: rbuff(:)

#ifdef USEMPI
      IF (p_is_master) THEN
         DO i = 0, p_np_io-1
            CALL mpi_recv (isrc, 1, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, &
               p_comm_glb, p_stat, p_err)
            CALL mpi_send (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_err)
            CALL mpi_recv (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
         ENDDO
      ENDIF
#endif

      IF (p_is_io) THEN

         IF (present(compress_level)) THEN
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_INT, &
               compress = compress_level)
         ELSE
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_INT)
         ENDIF
         
         allocate(rbuff(pixelset%nblkgrp))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff(iblkgrp)%val (pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
            CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
               rbuff(iblkgrp)%val, pixelset%vecgs%vcnt(:,iblk,jblk), &
               pixelset%vecgs%vdsp(:,iblk,jblk), MPI_INTEGER, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff(iblkgrp)%val = wdata(istt:iend)
#endif
         ENDDO

#ifdef USEMPI
         CALL mpi_send (p_iam_glb, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         CALL mpi_recv (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
#endif
         CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
         CALL nccheck( nf90_inq_ncid (ncid,  trim(dataname), grpid) )

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            CALL get_blockname (iblk, jblk, blockname)
            CALL nccheck( nf90_inq_varid(grpid, trim(vecname)//'_'//trim(blockname), varid) ) 
            CALL nccheck( nf90_put_var  (grpid, varid, rbuff(iblkgrp)%val) )

            deallocate (rbuff(iblkgrp)%val)
         ENDDO

         deallocate (rbuff)
         CALL nccheck( nf90_sync (ncid) )
         CALL nccheck( nf90_close(ncid) )
#ifdef USEMPI
         CALL mpi_send (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
#endif
      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               sbuff = wdata(istt:iend)
            ELSE
               allocate (sbuff (1))
            ENDIF

            CALL mpi_gatherv ( &
               sbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_INTEGER, &
               MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER, & ! insignificant on workers
               p_root, p_comm_group, p_err)

            IF (allocated(sbuff)) deallocate (sbuff)

         ENDDO

      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE ncio_write_vector_int32_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_logical_1d ( &
         filename, dataname, vecname, pixelset, wdata, compress_level)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: vecname
   type(pixelset_type), intent(in) :: pixelset
   logical, intent(in) :: wdata (:)

   integer, intent(in), optional :: compress_level

   ! Local variables
   integer :: ncid, grpid, varid, iblkall, iblkgrp, iblk, jblk, istt, iend, i, isrc, lock = 0
   character(len=8) :: blockname
   integer(1), allocatable :: sbuff(:)
   type(pointer_int8_1d), allocatable :: rbuff(:)

#ifdef USEMPI
      IF (p_is_master) THEN
         DO i = 0, p_np_io-1
            CALL mpi_recv (isrc, 1, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, &
               p_comm_glb, p_stat, p_err)
            CALL mpi_send (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_err)
            CALL mpi_recv (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
         ENDDO
      ENDIF
#endif

      IF (p_is_io) THEN

         IF (present(compress_level)) THEN
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_BYTE, &
               compress = compress_level)
         ELSE
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_BYTE)
         ENDIF
         
         allocate(rbuff(pixelset%nblkgrp))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff(iblkgrp)%val (pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
            CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER1, &
               rbuff(iblkgrp)%val, pixelset%vecgs%vcnt(:,iblk,jblk), &
               pixelset%vecgs%vdsp(:,iblk,jblk), MPI_INTEGER1, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            DO i = istt, iend
               IF(wdata(i))THEN
                  rbuff(iblkgrp)%val(i-istt+1) = 1
               ELSE
                  rbuff(iblkgrp)%val(i-istt+1) = 0
               ENDIF
            ENDDO
#endif
         ENDDO

#ifdef USEMPI
         CALL mpi_send (p_iam_glb, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         CALL mpi_recv (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
#endif
         CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
         CALL nccheck( nf90_inq_ncid (ncid,  trim(dataname), grpid) )

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)
         
            CALL get_blockname (iblk, jblk, blockname)
            CALL nccheck( nf90_inq_varid(grpid, trim(vecname)//'_'//trim(blockname), varid))
            CALL nccheck( nf90_put_var  (grpid, varid, rbuff(iblkgrp)%val) )
            
            deallocate (rbuff(iblkgrp)%val)
         ENDDO
            
         deallocate (rbuff)
         CALL nccheck( nf90_sync (ncid) )
         CALL nccheck( nf90_close(ncid) )
#ifdef USEMPI
         CALL mpi_send (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
#endif
      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               DO i = istt, iend
                  IF(wdata(i))THEN
                     sbuff(i-istt+1) = 1
                  ELSE
                     sbuff(i-istt+1) = 0
                  ENDIF
               ENDDO
            ELSE
               allocate (sbuff (1))
            ENDIF

            CALL mpi_gatherv ( &
               sbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_INTEGER1, &
               MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER1, & ! insignificant on workers
               p_root, p_comm_group, p_err)

            IF (allocated(sbuff)) deallocate (sbuff)

         ENDDO

      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE ncio_write_vector_logical_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_int64_1d ( &
         filename, dataname, vecname, pixelset, wdata, compress_level)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: vecname
   type(pixelset_type), intent(in) :: pixelset
   integer*8, intent(in) :: wdata (:)

   integer, intent(in), optional :: compress_level

   ! Local variables
   integer :: ncid, grpid, varid, iblkall, iblkgrp, iblk, jblk, istt, iend, i, isrc, lock = 0
   character(len=8) :: blockname
   integer*8, allocatable :: sbuff(:)
   type(pointer_int64_1d), allocatable :: rbuff(:)

#ifdef USEMPI
      IF (p_is_master) THEN
         DO i = 0, p_np_io-1
            CALL mpi_recv (isrc, 1, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, &
               p_comm_glb, p_stat, p_err)
            CALL mpi_send (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_err)
            CALL mpi_recv (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
         ENDDO
      ENDIF
#endif

      IF (p_is_io) THEN

         IF (present(compress_level)) THEN
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_INT64, &
               compress = compress_level)
         ELSE
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_INT64)
         ENDIF
         
         allocate(rbuff(pixelset%nblkgrp))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff(iblkgrp)%val (pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
            CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER8, &
               rbuff(iblkgrp)%val, pixelset%vecgs%vcnt(:,iblk,jblk), &
               pixelset%vecgs%vdsp(:,iblk,jblk), MPI_INTEGER8, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff(iblkgrp)%val = wdata(istt:iend)
#endif
         ENDDO

#ifdef USEMPI
         CALL mpi_send (p_iam_glb, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         CALL mpi_recv (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
#endif
         CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
         CALL nccheck( nf90_inq_ncid (ncid,  trim(dataname), grpid) )

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            CALL get_blockname (iblk, jblk, blockname)
            CALL nccheck( nf90_inq_varid(grpid, trim(vecname)//'_'//trim(blockname), varid))
            CALL nccheck( nf90_put_var  (grpid, varid, rbuff(iblkgrp)%val) )
            
            deallocate (rbuff(iblkgrp)%val)
         ENDDO

         deallocate (rbuff)
         CALL nccheck( nf90_sync (ncid) )
         CALL nccheck( nf90_close(ncid) )
#ifdef USEMPI
         CALL mpi_send (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
#endif
      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               sbuff = wdata(istt:iend)
            ELSE
               allocate (sbuff (1))
            ENDIF

            CALL mpi_gatherv ( &
               sbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_INTEGER8, &
               MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER8, & ! insignificant on workers
               p_root, p_comm_group, p_err)

            IF (allocated(sbuff)) deallocate (sbuff)

         ENDDO

      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE ncio_write_vector_int64_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_1d ( &
         filename, dataname, vecname, pixelset, wdata, compress_level)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: vecname
   type(pixelset_type), intent(in) :: pixelset
   real(r8), intent(in) :: wdata (:)
   
   integer, intent(in), optional :: compress_level

   ! Local variables
   integer :: ncid, grpid, varid, iblkall, iblkgrp, iblk, jblk, istt, iend, i, isrc, lock = 0
   character(len=8) :: blockname
   real(r8), allocatable :: sbuff(:)
   type(pointer_real8_1d), allocatable :: rbuff(:)

#ifdef USEMPI
      IF (p_is_master) THEN
         DO i = 0, p_np_io-1
            CALL mpi_recv (isrc, 1, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, &
               p_comm_glb, p_stat, p_err)
            CALL mpi_send (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_err)
            CALL mpi_recv (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
         ENDDO
      ENDIF
#endif

      IF (p_is_io) THEN
         
         IF (present(compress_level)) THEN
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_DOUBLE, &
               compress = compress_level)
         ELSE
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_DOUBLE)
         ENDIF

         allocate(rbuff(pixelset%nblkgrp))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff(iblkgrp)%val (pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
            CALL mpi_gatherv ( MPI_IN_PLACE, 0, MPI_REAL8, &
               rbuff(iblkgrp)%val, pixelset%vecgs%vcnt(:,iblk,jblk), &
               pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff(iblkgrp)%val = wdata(istt:iend)
#endif
         ENDDO

#ifdef USEMPI
         CALL mpi_send (p_iam_glb, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         CALL mpi_recv (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
#endif
         CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
         CALL nccheck( nf90_inq_ncid (ncid,  trim(dataname), grpid) )

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            CALL get_blockname (iblk, jblk, blockname)
            CALL nccheck( nf90_inq_varid(grpid, trim(vecname)//'_'//trim(blockname), varid))
            CALL nccheck( nf90_put_var  (grpid, varid, rbuff(iblkgrp)%val) )
            
            deallocate (rbuff(iblkgrp)%val)
         ENDDO

         deallocate (rbuff)
         CALL nccheck( nf90_sync (ncid) )
         CALL nccheck( nf90_close(ncid) )
#ifdef USEMPI
         CALL mpi_send (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
#endif
      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               sbuff = wdata(istt:iend)
            ELSE
               allocate (sbuff (1))
            ENDIF

            CALL mpi_gatherv ( &
               sbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
               MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
               p_root, p_comm_group, p_err)

            IF (allocated(sbuff)) deallocate (sbuff)

         ENDDO

      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE ncio_write_vector_real8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_2d ( &
         filename, dataname, dim1name, ndim1, &
         vecname, pixelset, wdata, compress_level)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dim1name, vecname
   integer,  intent(in) :: ndim1
   type(pixelset_type), intent(in) :: pixelset
   real(r8), intent(in) :: wdata (:,:)

   integer,  intent(in), optional :: compress_level

   ! Local variables
   integer :: ncid, grpid, varid, iblkall, iblkgrp, iblk, jblk, istt, iend, i, isrc, lock = 0
   character(len=8) :: blockname
   real(r8), allocatable :: sbuff(:,:)
   type(pointer_real8_2d), allocatable :: rbuff(:)

#ifdef USEMPI
      IF (p_is_master) THEN
         DO i = 0, p_np_io-1
            CALL mpi_recv (isrc, 1, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, &
               p_comm_glb, p_stat, p_err)
            CALL mpi_send (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_err)
            CALL mpi_recv (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
         ENDDO
      ENDIF
#endif

      IF (p_is_io) THEN

         IF (present(compress_level)) THEN
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_DOUBLE, &
               dim1name = dim1name, compress = compress_level)
         ELSE
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_DOUBLE, &
               dim1name = dim1name)
         ENDIF
         
         allocate(rbuff(pixelset%nblkgrp))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff(iblkgrp)%val (ndim1, pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
            CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_REAL8, &
               rbuff(iblkgrp)%val, ndim1 * pixelset%vecgs%vcnt(:,iblk,jblk), &
               ndim1 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff(iblkgrp)%val = wdata(:,istt:iend)
#endif
         ENDDO

#ifdef USEMPI
         CALL mpi_send (p_iam_glb, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         CALL mpi_recv (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
#endif
         CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
         CALL nccheck( nf90_inq_ncid (ncid,  trim(dataname), grpid) )

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            CALL get_blockname (iblk, jblk, blockname)
            CALL nccheck( nf90_inq_varid(grpid, trim(vecname)//'_'//trim(blockname), varid))
            CALL nccheck( nf90_put_var  (grpid, varid, rbuff(iblkgrp)%val) )
            
            deallocate (rbuff(iblkgrp)%val)
         ENDDO

         deallocate (rbuff)
         CALL nccheck( nf90_sync (ncid) )
         CALL nccheck( nf90_close(ncid) )
#ifdef USEMPI
         CALL mpi_send (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
#endif
      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (sbuff (ndim1,pixelset%vecgs%vlen(iblk,jblk)))
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               sbuff = wdata(:,istt:iend)
            ELSE
               allocate (sbuff (1,1))
            ENDIF

            CALL mpi_gatherv ( &
               sbuff, ndim1 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
               MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
               p_root, p_comm_group, p_err)

            IF (allocated(sbuff)) deallocate (sbuff)

         ENDDO

      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE ncio_write_vector_real8_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_3d ( &
         filename, dataname, dim1name, ndim1, dim2name, ndim2, &
         vecname, pixelset, wdata, compress_level)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dim1name, dim2name, vecname
   type(pixelset_type), intent(in) :: pixelset
   integer,  intent(in) :: ndim1, ndim2
   real(r8), intent(in) :: wdata (:,:,:)
   
   integer,  intent(in), optional :: compress_level

   ! Local variables
   integer :: ncid, grpid, varid, iblkall, iblkgrp, iblk, jblk, istt, iend, i, isrc, lock = 0
   character(len=8) :: blockname
   real(r8), allocatable :: sbuff(:,:,:)
   type(pointer_real8_3d), allocatable :: rbuff(:)

#ifdef USEMPI
      IF (p_is_master) THEN
         DO i = 0, p_np_io-1
            CALL mpi_recv (isrc, 1, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, &
               p_comm_glb, p_stat, p_err)
            CALL mpi_send (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_err)
            CALL mpi_recv (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
         ENDDO
      ENDIF
#endif

      IF (p_is_io) THEN

         IF (present(compress_level)) THEN
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_DOUBLE, &
               dim1name = dim1name, dim2name = dim2name, compress = compress_level)
         ELSE
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_DOUBLE, &
               dim1name = dim1name, dim2name = dim2name)
         ENDIF
         
         allocate(rbuff(pixelset%nblkgrp))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff(iblkgrp)%val (ndim1, ndim2, pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
            CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_REAL8, &
               rbuff(iblkgrp)%val, ndim1 * ndim2 * pixelset%vecgs%vcnt(:,iblk,jblk), &
               ndim1 * ndim2 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff(iblkgrp)%val = wdata(:,:,istt:iend)
#endif
         ENDDO

#ifdef USEMPI
         CALL mpi_send (p_iam_glb, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         CALL mpi_recv (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
#endif
         CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
         CALL nccheck( nf90_inq_ncid (ncid,  trim(dataname), grpid) )

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            CALL get_blockname (iblk, jblk, blockname)
            CALL nccheck( nf90_inq_varid(grpid, trim(vecname)//'_'//trim(blockname), varid))
            CALL nccheck( nf90_put_var  (grpid, varid, rbuff(iblkgrp)%val) )
            
            deallocate (rbuff(iblkgrp)%val)
         ENDDO
         
         deallocate (rbuff)
         CALL nccheck( nf90_sync (ncid) )
         CALL nccheck( nf90_close(ncid) )
#ifdef USEMPI
         CALL mpi_send (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
#endif
      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (sbuff (ndim1,ndim2,pixelset%vecgs%vlen(iblk,jblk)))
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               sbuff = wdata(:,:,istt:iend)
            ELSE
               allocate (sbuff (1,1,1))
            ENDIF

            CALL mpi_gatherv ( sbuff, &
               ndim1 * ndim2 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
               MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
               p_root, p_comm_group, p_err)

            IF (allocated(sbuff)) deallocate (sbuff)

         ENDDO

      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE ncio_write_vector_real8_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_4d ( &
         filename, dataname, dim1name, ndim1, dim2name, ndim2, dim3name, ndim3, &
         vecname, pixelset, wdata, compress_level)

   IMPLICIT NONE

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: dataname
   character(len=*), intent(in) :: dim1name, dim2name, dim3name, vecname
   integer,  intent(in) :: ndim1, ndim2, ndim3
   type(pixelset_type), intent(in) :: pixelset
   real(r8), intent(in) :: wdata (:,:,:,:)
   
   integer,  intent(in), optional :: compress_level

   ! Local variables
   integer :: ncid, grpid, varid, iblkall, iblkgrp, iblk, jblk, istt, iend, i, isrc, lock = 0
   character(len=8) :: blockname
   real(r8), allocatable :: sbuff(:,:,:,:)
   type(pointer_real8_4d), allocatable :: rbuff(:)

#ifdef USEMPI
      IF (p_is_master) THEN
         DO i = 0, p_np_io-1
            CALL mpi_recv (isrc, 1, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, &
               p_comm_glb, p_stat, p_err)
            CALL mpi_send (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_err)
            CALL mpi_recv (lock, 1, MPI_INTEGER, isrc, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
         ENDDO
      ENDIF
#endif

      IF (p_is_io) THEN

         IF (present(compress_level)) THEN
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_DOUBLE, &
               dim1name = dim1name, dim2name = dim2name, dim3name = dim3name, &
               compress = compress_level)
         ELSE
            CALL ncio_define_variable_vector (filename, pixelset, vecname, dataname, NF90_DOUBLE, &
               dim1name = dim1name, dim2name = dim2name, dim3name = dim3name)
         ENDIF
         
         allocate(rbuff(pixelset%nblkgrp))

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            allocate (rbuff(iblkgrp)%val (ndim1, ndim2, ndim3, pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
            CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_REAL8, &
               rbuff(iblkgrp)%val, ndim1 * ndim2 * ndim3 * pixelset%vecgs%vcnt(:,iblk,jblk), &
               ndim1 * ndim2 * ndim3 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
               p_root, p_comm_group, p_err)
#else
            istt = pixelset%vecgs%vstt(iblk,jblk)
            iend = pixelset%vecgs%vend(iblk,jblk)
            rbuff(iblkgrp)%val = wdata(:,:,:,istt:iend)
#endif
         ENDDO

#ifdef USEMPI
         CALL mpi_send (p_iam_glb, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         CALL mpi_recv (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
#endif
         CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
         CALL nccheck( nf90_inq_ncid (ncid,  trim(dataname), grpid) )

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            CALL get_blockname (iblk, jblk, blockname)
            CALL nccheck( nf90_inq_varid(grpid, trim(vecname)//'_'//trim(blockname), varid))
            CALL nccheck( nf90_put_var  (grpid, varid, rbuff(iblkgrp)%val) )
            
            deallocate (rbuff(iblkgrp)%val)
         ENDDO

         deallocate (rbuff)
         CALL nccheck( nf90_sync (ncid) )
         CALL nccheck( nf90_close(ncid) )
#ifdef USEMPI
         CALL mpi_send (lock, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
#endif
      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iblkgrp = 1, pixelset%nblkgrp
            iblk = pixelset%xblkgrp(iblkgrp)
            jblk = pixelset%yblkgrp(iblkgrp)

            IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
               allocate (sbuff (ndim1,ndim2,ndim3,pixelset%vecgs%vlen(iblk,jblk)))
               istt = pixelset%vecgs%vstt(iblk,jblk)
               iend = pixelset%vecgs%vend(iblk,jblk)
               sbuff = wdata(:,:,:,istt:iend)
            ELSE
               allocate (sbuff (1,1,1,1))
            ENDIF

            CALL mpi_gatherv ( sbuff, &
               ndim1 * ndim2 * ndim3 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
               MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
               p_root, p_comm_group, p_err)

            IF (allocated(sbuff)) deallocate (sbuff)

         ENDDO

      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE ncio_write_vector_real8_4d

END MODULE MOD_NetCDFVector

#endif
