#include <define.h>

MODULE ncio_serial

   USE netcdf
   USE precision
   IMPLICIT NONE

   ! PUBLIC subroutines

   PUBLIC :: ncio_create_file
   PUBLIC :: check_ncfile_exist

   INTERFACE ncio_put_attr
      MODULE procedure ncio_put_attr_str
      MODULE procedure ncio_put_attr_real8
   END INTERFACE ncio_put_attr

   INTERFACE ncio_get_attr
      MODULE procedure ncio_get_attr_str
      MODULE procedure ncio_get_attr_real8
   END INTERFACE ncio_get_attr

   PUBLIC :: ncio_var_exist
   PUBLIC :: ncio_inquire_varsize
   PUBLIC :: ncio_inquire_length

   INTERFACE ncio_read_serial
      MODULE procedure ncio_read_serial_int32_0d 
      MODULE procedure ncio_read_serial_real8_0d 
      MODULE procedure ncio_read_serial_int8_1d 
      MODULE procedure ncio_read_serial_int32_1d 
      MODULE procedure ncio_read_serial_real8_1d 
      MODULE procedure ncio_read_serial_int8_2d 
      MODULE procedure ncio_read_serial_int16_2d 
      MODULE procedure ncio_read_serial_int32_2d 
      MODULE procedure ncio_read_serial_real4_2d 
      MODULE procedure ncio_read_serial_real8_2d 
      MODULE procedure ncio_read_serial_int32_3d 
      MODULE procedure ncio_read_serial_real8_3d 
      MODULE procedure ncio_read_serial_real8_4d 
      MODULE procedure ncio_read_serial_real8_5d 
   END INTERFACE ncio_read_serial

   INTERFACE ncio_read_bcast_serial
      MODULE procedure ncio_read_bcast_serial_int32_0d 
      MODULE procedure ncio_read_bcast_serial_real8_0d 
      MODULE procedure ncio_read_bcast_serial_int32_1d 
      MODULE procedure ncio_read_bcast_serial_int32_2d 
      MODULE procedure ncio_read_bcast_serial_real8_1d 
      MODULE procedure ncio_read_bcast_serial_real8_2d 
      MODULE procedure ncio_read_bcast_serial_real8_3d
      MODULE procedure ncio_read_bcast_serial_real8_4d
      MODULE procedure ncio_read_bcast_serial_real8_5d
      MODULE procedure ncio_read_bcast_serial_logical_1d
   END INTERFACE ncio_read_bcast_serial

   interface ncio_read_part_serial
      MODULE procedure ncio_read_part_serial_int32_2d
   END interface ncio_read_part_serial
   

   interface ncio_define_dimension
      MODULE procedure ncio_define_dimension_int32
      MODULE procedure ncio_define_dimension_int64
   END interface ncio_define_dimension

   INTERFACE ncio_write_serial
      MODULE procedure ncio_write_serial_int32_0d 
      MODULE procedure ncio_write_serial_real8_0d 
      MODULE procedure ncio_write_serial_int8_1d 
      MODULE procedure ncio_write_serial_int32_1d 
      MODULE procedure ncio_write_serial_real8_1d 
      MODULE procedure ncio_write_serial_logical_1d 
      MODULE procedure ncio_write_serial_int8_2d 
      MODULE procedure ncio_write_serial_int16_2d 
      MODULE procedure ncio_write_serial_int32_2d 
      MODULE procedure ncio_write_serial_real4_2d 
      MODULE procedure ncio_write_serial_real8_2d 
      MODULE procedure ncio_write_serial_int32_3d 
      MODULE procedure ncio_write_serial_real8_3d 
      MODULE procedure ncio_write_serial_real8_4d 
      MODULE procedure ncio_write_serial_real8_5d 
   END INTERFACE ncio_write_serial

   PUBLIC :: ncio_write_time

   INTERFACE ncio_write_serial_time
      MODULE procedure ncio_write_serial_real8_1d_time 
      MODULE procedure ncio_write_serial_real8_2d_time 
      MODULE procedure ncio_write_serial_real8_3d_time 
      MODULE procedure ncio_write_serial_real8_4d_time 
   END INTERFACE ncio_write_serial_time

CONTAINS

   ! ----
   SUBROUTINE nccheck (status)
      INTEGER, INTENT(IN) :: status

      IF (status /= NF90_NOERR) THEN
         print *, trim(nf90_strerror(status))
         stop 2
      ENDIF
   END SUBROUTINE nccheck

   ! ----
   SUBROUTINE check_ncfile_exist (filename)

      USE spmd_task
      IMPLICIT NONE

      CHARACTER(len=256), INTENT(IN) :: filename
      ! Local Variables
      LOGICAL :: fexists

      inquire (file=trim(filename), exist=fexists)
      IF (.not. fexists) THEN 
         write(*,*) trim(filename), ' does not exist.'
#ifdef USEMPI
         CALL mpi_abort (p_comm_glb, p_err)
#else
         stop 2
#endif
      ENDIF

   END SUBROUTINE check_ncfile_exist

   ! ----
   SUBROUTINE ncio_create_file (filename)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)  :: filename

      ! Local Variables
      INTEGER :: ncid

      CALL nccheck( nf90_create(trim(filename), ior(NF90_CLOBBER,NF90_NETCDF4), ncid) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_create_file

   ! ----
   SUBROUTINE ncio_put_attr_str (filename, varname, attrname, attrval)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)  :: filename, varname, attrname, attrval

      ! Local Variables
      INTEGER :: ncid, varid

      CALL nccheck( nf90_open (trim(filename), NF90_WRITE, ncid) )
      CALL nccheck (nf90_inq_varid (ncid, trim(varname), varid))
      CALL nccheck (nf90_redef (ncid))
      CALL nccheck (nf90_put_att (ncid, varid, trim(attrname), trim(attrval)))
      CALL nccheck (nf90_enddef (ncid))
      CALL nccheck( nf90_close (ncid))

   END SUBROUTINE ncio_put_attr_str
   
   ! ----
   SUBROUTINE ncio_get_attr_str (filename, varname, attrname, attrval)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)  :: filename, varname, attrname
      CHARACTER(len=*), intent(out) :: attrval

      ! Local Variables
      INTEGER :: ncid, varid

      CALL nccheck( nf90_open (trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck (nf90_inq_varid (ncid, trim(varname), varid))
      CALL nccheck (nf90_get_att   (ncid, varid, trim(attrname), attrval))
      CALL nccheck( nf90_close (ncid))

   END SUBROUTINE ncio_get_attr_str

   ! ----
   SUBROUTINE ncio_get_attr_real8 (filename, varname, attrname, attrval)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)  :: filename, varname, attrname
      REAL(r8), intent(out) :: attrval

      ! Local Variables
      INTEGER :: ncid, varid

      CALL nccheck( nf90_open (trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck (nf90_inq_varid (ncid, trim(varname), varid))
      CALL nccheck (nf90_get_att (ncid, varid, trim(attrname), attrval))
      CALL nccheck (nf90_close (ncid))

   END SUBROUTINE ncio_get_attr_real8

   ! ----
   SUBROUTINE ncio_put_attr_real8 (filename, varname, attrname, attrval)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)  :: filename, varname, attrname
      REAL(r8),         intent(in)  :: attrval

      ! Local Variables
      INTEGER :: ncid, varid

      CALL nccheck( nf90_open (trim(filename), NF90_WRITE, ncid) )
      CALL nccheck (nf90_inq_varid (ncid, trim(varname), varid))
      CALL nccheck (nf90_redef (ncid))
      CALL nccheck (nf90_put_att (ncid, varid, trim(attrname), attrval))
      CALL nccheck (nf90_enddef (ncid))
      CALL nccheck( nf90_close (ncid))

   END SUBROUTINE ncio_put_attr_real8


   !---------------------------------------------------------
   SUBROUTINE ncio_inquire_varsize (filename, dataname, varsize)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, allocatable, intent(out) :: varsize(:)

      ! Local variables
      INTEGER :: ncid, varid, ndims, idm
      INTEGER, allocatable :: dimids(:)

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )

      CALL nccheck( nf90_inquire_variable(ncid, varid, ndims = ndims) )
      allocate (dimids(ndims))
      CALL nccheck( nf90_inquire_variable(ncid, varid, dimids = dimids) )

      allocate (varsize(ndims))
      DO idm = 1, ndims
         CALL nccheck( nf90_inquire_dimension(ncid, dimids(idm), len = varsize(idm)) )
      ENDDO 
      
      CALL nccheck( nf90_close(ncid) )
      deallocate (dimids)

   END SUBROUTINE ncio_inquire_varsize

   !---------------------------------------------------------
   LOGICAL function ncio_var_exist (filename, dataname)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname

      ! Local variables
      INTEGER :: ncid, varid, status

      status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
      IF (status == nf90_noerr) THEN
         status = nf90_inq_varid(ncid, trim(dataname), varid)
         ncio_var_exist = (status == nf90_noerr)
         CALL nccheck( nf90_close(ncid) )
      ELSE
         ncio_var_exist = .false.
      ENDIF

      IF (.not. ncio_var_exist) THEN
         write(*,*) 'Warning: ', trim(dataname), ' not found in ', trim(filename)
      ENDIF

   END FUNCTION ncio_var_exist

   !---------------------------------------------------------
   SUBROUTINE ncio_inquire_length (filename, dataname, length)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, intent(out) :: length

      ! Local variables
      INTEGER :: ncid, varid, ndims
      INTEGER, allocatable :: dimids(:)

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )

      CALL nccheck( nf90_inquire_variable(ncid, varid, ndims = ndims) )
      allocate (dimids(ndims))
      CALL nccheck( nf90_inquire_variable(ncid, varid, dimids = dimids) )
      CALL nccheck( nf90_inquire_dimension(ncid, dimids(ndims), len = length) )
      
      CALL nccheck( nf90_close(ncid) )
      deallocate (dimids)

   END SUBROUTINE ncio_inquire_length
   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_int32_0d (filename, dataname, rdata)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, intent(out) :: rdata 

      ! Local variables
      INTEGER :: ncid, varid

      CALL check_ncfile_exist (filename)

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_read_serial_int32_0d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_real8_0d (filename, dataname, rdata)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), intent(out) :: rdata 

      ! Local variables
      INTEGER :: ncid, varid

      CALL check_ncfile_exist (filename)

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_read_serial_real8_0d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_int8_1d (filename, dataname, rdata)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER(1), allocatable, intent(out) :: rdata (:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_int8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_int32_1d (filename, dataname, rdata)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, allocatable, intent(out) :: rdata (:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_int32_1d


   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_real8_1d (filename, dataname, rdata)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), allocatable, intent(out) :: rdata (:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_real8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_int8_2d (filename, dataname, rdata)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER(1), allocatable, intent(out) :: rdata (:,:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_int8_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_int16_2d (filename, dataname, rdata)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER(2), allocatable, intent(out) :: rdata (:,:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_int16_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_int32_2d (filename, dataname, rdata)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, allocatable, intent(out) :: rdata (:,:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_int32_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_real4_2d (filename, dataname, rdata)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(4), allocatable, intent(out) :: rdata (:,:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_real4_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_real8_2d (filename, dataname, rdata)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), allocatable, intent(out) :: rdata (:,:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_real8_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_int32_3d (filename, dataname, rdata)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, allocatable, intent(out) :: rdata (:,:,:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2), varsize(3)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_int32_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_real8_3d (filename, dataname, rdata)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), allocatable, intent(out) :: rdata (:,:,:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2), varsize(3)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_real8_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_real8_4d (filename, dataname, rdata)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), allocatable, intent(out) :: rdata (:,:,:,:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2), varsize(3), varsize(4)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_real8_4d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_serial_real8_5d (filename, dataname, rdata)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), allocatable, intent(out) :: rdata (:,:,:,:,:)

      ! Local variables
      INTEGER :: ncid, varid
      INTEGER, allocatable :: varsize(:)

      CALL check_ncfile_exist (filename)

      CALL ncio_inquire_varsize(filename, dataname, varsize)
      allocate (rdata (varsize(1), varsize(2), varsize(3), varsize(4), varsize(5)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata) )
      CALL nccheck( nf90_close(ncid) )

      deallocate (varsize)

   END SUBROUTINE ncio_read_serial_real8_5d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_bcast_serial_int32_0d (filename, dataname, rdata)

      USE netcdf
      USE spmd_task
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, intent(out) :: rdata 

      IF (p_is_master) THEN
         CALL ncio_read_serial_int32_0d (filename, dataname, rdata)
      ENDIF
         
#ifdef USEMPI
      CALL mpi_bcast (rdata, 1, MPI_INTEGER, p_root, p_comm_glb, p_err)
#endif 

   END SUBROUTINE ncio_read_bcast_serial_int32_0d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_bcast_serial_real8_0d (filename, dataname, rdata)

      USE spmd_task
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), intent(out) :: rdata 

      IF (p_is_master) THEN
         CALL ncio_read_serial_real8_0d (filename, dataname, rdata)
      ENDIF
      
#ifdef USEMPI
      CALL mpi_bcast (rdata, 1, MPI_REAL8, p_root, p_comm_glb, p_err)
#endif 

   END SUBROUTINE ncio_read_bcast_serial_real8_0d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_bcast_serial_int32_1d (filename, dataname, rdata)

      USE spmd_task
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, allocatable, intent(out) :: rdata (:)
      INTEGER :: vlen

      IF (p_is_master) THEN
         CALL ncio_read_serial_int32_1d(filename, dataname, rdata)
         vlen = size(rdata)
      ENDIF
      
#ifdef USEMPI
      CALL mpi_bcast (vlen, 1, MPI_INTEGER, p_root, p_comm_glb, p_err)
      IF (.not. p_is_master)  allocate (rdata (vlen))
      CALL mpi_bcast (rdata, vlen, MPI_INTEGER, p_root, p_comm_glb, p_err)
#endif

   END SUBROUTINE ncio_read_bcast_serial_int32_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_bcast_serial_int32_2d (filename, dataname, rdata)

      USE spmd_task
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, allocatable, intent(out) :: rdata (:,:)
      INTEGER :: vsize(2)

      IF (p_is_master) THEN
         CALL ncio_read_serial_int32_2d(filename, dataname, rdata)
         vsize = shape(rdata)
      ENDIF
      
#ifdef USEMPI
      CALL mpi_bcast (vsize, 2, MPI_INTEGER, p_root, p_comm_glb, p_err)
      IF (.not. p_is_master)  allocate (rdata (vsize(1), vsize(2)))
      CALL mpi_bcast (rdata, vsize(1)*vsize(2), MPI_INTEGER, p_root, p_comm_glb, p_err)
#endif 

   END SUBROUTINE ncio_read_bcast_serial_int32_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_bcast_serial_real8_1d (filename, dataname, rdata)

      USE netcdf
      USE spmd_task
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), allocatable, intent(out) :: rdata (:)
      INTEGER :: vlen

      IF (p_is_master) THEN
         CALL ncio_read_serial_real8_1d(filename, dataname, rdata)
         vlen = size(rdata)
      ENDIF

#ifdef USEMPI
      CALL mpi_bcast (vlen, 1, MPI_INTEGER, p_root, p_comm_glb, p_err)
      IF (.not. p_is_master)  allocate (rdata (vlen))
      CALL mpi_bcast (rdata, vlen, MPI_REAL8, p_root, p_comm_glb, p_err)
#endif 

   END SUBROUTINE ncio_read_bcast_serial_real8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_bcast_serial_real8_2d (filename, dataname, rdata)

      USE netcdf
      USE spmd_task
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), allocatable, intent(out) :: rdata (:,:)
      INTEGER :: vsize(2)

      IF (p_is_master) THEN
         CALL ncio_read_serial_real8_2d(filename, dataname, rdata)
         vsize = shape(rdata)
      ENDIF
      
#ifdef USEMPI
      CALL mpi_bcast (vsize, 2, MPI_INTEGER, p_root, p_comm_glb, p_err)
      IF (.not. p_is_master)  allocate (rdata (vsize(1),vsize(2)))
      CALL mpi_bcast (rdata, vsize(1)*vsize(2), MPI_REAL8, p_root, p_comm_glb, p_err)
#endif 

   END SUBROUTINE ncio_read_bcast_serial_real8_2d
   
   !---------------------------------------------------------
   SUBROUTINE ncio_read_bcast_serial_real8_3d (filename, dataname, rdata)

      USE netcdf
      USE spmd_task
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), allocatable, intent(out) :: rdata (:,:,:)
      INTEGER :: vsize(3)

      IF (p_is_master) THEN
         CALL ncio_read_serial_real8_3d(filename, dataname, rdata)
         vsize = shape(rdata)
      ENDIF
      
#ifdef USEMPI
      CALL mpi_bcast (vsize, 3, MPI_INTEGER, p_root, p_comm_glb, p_err)
      IF (.not. p_is_master)  allocate (rdata (vsize(1),vsize(2),vsize(3)))
      CALL mpi_bcast (rdata, vsize(1)*vsize(2)*vsize(3), MPI_REAL8, p_root, p_comm_glb, p_err)
#endif 

   END SUBROUTINE ncio_read_bcast_serial_real8_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_bcast_serial_real8_4d (filename, dataname, rdata)

      USE netcdf
      USE spmd_task
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), allocatable, intent(out) :: rdata (:,:,:,:)
      INTEGER :: vsize(4)

      IF (p_is_master) THEN
         CALL ncio_read_serial_real8_4d(filename, dataname, rdata)
         vsize = shape(rdata)
      ENDIF
      
#ifdef USEMPI
      CALL mpi_bcast (vsize, 4, MPI_INTEGER, p_root, p_comm_glb, p_err)
      IF (.not. p_is_master)  allocate (rdata (vsize(1),vsize(2),vsize(3),vsize(4)))
      CALL mpi_bcast (rdata, vsize(1)*vsize(2)*vsize(3)*vsize(4), MPI_REAL8, p_root, p_comm_glb, p_err)
#endif 

   END SUBROUTINE ncio_read_bcast_serial_real8_4d

      !---------------------------------------------------------
   SUBROUTINE ncio_read_bcast_serial_real8_5d (filename, dataname, rdata)

      USE netcdf
      USE spmd_task
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), allocatable, intent(out) :: rdata (:,:,:,:,:)
      INTEGER :: vsize(5)

      IF (p_is_master) THEN
         CALL ncio_read_serial_real8_5d(filename, dataname, rdata)
         vsize = shape(rdata)
      ENDIF
      
#ifdef USEMPI
      CALL mpi_bcast (vsize, 5, MPI_INTEGER, p_root, p_comm_glb, p_err)
      IF (.not. p_is_master)  allocate (rdata (vsize(1),vsize(2),vsize(3),vsize(4),vsize(5)))
      CALL mpi_bcast (rdata, vsize(1)*vsize(2)*vsize(3)*vsize(4)*vsize(5), MPI_REAL8, p_root, p_comm_glb, p_err)
#endif 

   END SUBROUTINE ncio_read_bcast_serial_real8_5d
   
   ! -------------------------------
   SUBROUTINE ncio_read_bcast_serial_logical_1d (filename, dataname, rdata)

      USE spmd_task
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      LOGICAL, allocatable, intent(out) :: rdata (:)
      INTEGER :: vlen
      INTEGER(1), allocatable :: rdata_byte(:)

      IF (p_is_master) THEN
         CALL ncio_read_serial_int8_1d(filename, dataname, rdata_byte)
         vlen = size(rdata_byte)
         
         allocate(rdata(vlen))
         rdata = (rdata_byte == 1)

         deallocate (rdata_byte)
      ENDIF
      
#ifdef USEMPI
      CALL mpi_bcast (vlen, 1, MPI_INTEGER, p_root, p_comm_glb, p_err)
      IF (.not. p_is_master)  allocate (rdata (vlen))
      CALL mpi_bcast (rdata, vlen, MPI_LOGICAL, p_root, p_comm_glb, p_err)
#endif

   END SUBROUTINE ncio_read_bcast_serial_logical_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_part_serial_int32_2d (filename, dataname, datastt, dataend, rdata)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, intent(in) :: datastt(2), dataend(2)
      INTEGER, allocatable, intent(out) :: rdata (:,:)

      ! Local variables
      INTEGER :: ncid, varid

      CALL check_ncfile_exist (filename)

      allocate (rdata (datastt(1):dataend(1), datastt(2):dataend(2)) )

      CALL nccheck( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
      CALL nccheck( nf90_get_var(ncid, varid, rdata, &
         (/datastt(1),datastt(2)/), (/dataend(1)-datastt(1)+1, dataend(2)-datastt(2)+1/)) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_read_part_serial_int32_2d

   ! -------------------------------
   SUBROUTINE ncio_define_dimension_int32 (filename, dimname, dimlen)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dimname
      INTEGER, intent(in) :: dimlen

      ! Local variables
      INTEGER :: ncid, dimid, status
      INTEGER :: varid

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )

      status = nf90_inq_dimid(ncid, trim(dimname), dimid)
      IF (status /= NF90_NOERR) THEN
         CALL nccheck (nf90_redef(ncid))
         IF (dimlen == 0) THEN
            CALL nccheck( nf90_def_dim(ncid, trim(dimname), NF90_UNLIMITED, dimid) )
         ELSE
            CALL nccheck( nf90_def_dim(ncid, trim(dimname), dimlen, dimid) )
         ENDIF
         if (trim(dimname) .eq. 'lon') then
            !print *, 'lon-def'
            call nccheck( nf90_def_var(ncid, 'lon', nf90_float, (/dimid/), varid) )
            call nccheck( nf90_put_att(ncid, varid, 'long_name','longitude') )
            call nccheck( nf90_put_att(ncid, varid, 'units','degrees_east') )
         elseif (trim(dimname) .eq.'lat') then
            !print *, 'lat-def'
            call nccheck( nf90_def_var(ncid, 'lat', nf90_float, (/dimid/), varid) )
            call nccheck( nf90_put_att(ncid, varid, 'long_name','latitude') )
            call nccheck( nf90_put_att(ncid, varid, 'units','degrees_north') )
         elseif (trim(dimname) .eq.'lat_cama') then
               !print *, 'lat-def'
               call nccheck( nf90_def_var(ncid, 'lat_cama', nf90_float, (/dimid/), varid) )
               call nccheck( nf90_put_att(ncid, varid, 'long_name','latitude') )
               call nccheck( nf90_put_att(ncid, varid, 'units','degrees_north') )
         elseif (trim(dimname) .eq.'lon_cama') then
            call nccheck( nf90_def_var(ncid, 'lon_cama', nf90_float, (/dimid/), varid) )
            call nccheck( nf90_put_att(ncid, varid, 'long_name','longitude') )
            call nccheck( nf90_put_att(ncid, varid, 'units','degrees_east') )                            
         endif 
         CALL nccheck (nf90_enddef(ncid))
      ENDIF

      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_define_dimension_int32

   ! -------------------------------
   SUBROUTINE ncio_define_dimension_int64 (filename, dimname, dimlen)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dimname
      INTEGER*8, intent(in) :: dimlen

      ! Local variables
      INTEGER :: ncid, dimid, status
      INTEGER :: varid

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )

      status = nf90_inq_dimid(ncid, trim(dimname), dimid)
      IF (status /= NF90_NOERR) THEN
         CALL nccheck (nf90_redef(ncid))
         IF (dimlen == 0) THEN
            CALL nccheck( nf90_def_dim(ncid, trim(dimname), NF90_UNLIMITED, dimid) )
         ELSE
            CALL nccheck( nf90_def_dim(ncid, trim(dimname), int(dimlen), dimid) )
         ENDIF
         if (trim(dimname) .eq. 'lon') then
            !print *, 'lon-def'
            call nccheck( nf90_def_var(ncid, 'lon', nf90_float, (/dimid/), varid) )
            call nccheck( nf90_put_att(ncid, varid, 'long_name','longitude') )
            call nccheck( nf90_put_att(ncid, varid, 'units','degrees_east') )
         elseif (trim(dimname) .eq.'lat') then
            !print *, 'lat-def'
            call nccheck( nf90_def_var(ncid, 'lat', nf90_float, (/dimid/), varid) )
            call nccheck( nf90_put_att(ncid, varid, 'long_name','latitude') )
            call nccheck( nf90_put_att(ncid, varid, 'units','degrees_north') )
         elseif (trim(dimname) .eq.'lat_cama') then
               !print *, 'lat-def'
               call nccheck( nf90_def_var(ncid, 'lat_cama', nf90_float, (/dimid/), varid) )
               call nccheck( nf90_put_att(ncid, varid, 'long_name','latitude') )
               call nccheck( nf90_put_att(ncid, varid, 'units','degrees_north') )
         elseif (trim(dimname) .eq.'lon_cama') then
            call nccheck( nf90_def_var(ncid, 'lon_cama', nf90_float, (/dimid/), varid) )
            call nccheck( nf90_put_att(ncid, varid, 'long_name','longitude') )
            call nccheck( nf90_put_att(ncid, varid, 'units','degrees_east') )          
                                 
         endif 
         CALL nccheck (nf90_enddef(ncid))
      ENDIF

      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_define_dimension_int64

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int32_0d (filename, dataname, wdata)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, intent(in) :: wdata 

      ! Local variables
      INTEGER :: ncid, varid, status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, varid = varid))
         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int32_0d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_0d (filename, dataname, wdata)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), intent(in) :: wdata 

      ! Local variables
      INTEGER :: ncid, varid, status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         CALL nccheck (nf90_redef(ncid))
         CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, varid = varid))
         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_0d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int8_1d (filename, dataname, wdata, dimname, compress)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER(1), intent(in) :: wdata (:)

      CHARACTER(len=*), intent(in), optional :: dimname
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid, status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. present(dimname)) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dimname), dimid))

         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_BYTE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_BYTE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int32_1d (filename, dataname, wdata, dimname, compress)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, intent(in) :: wdata (:)

      CHARACTER(len=*), intent(in), optional :: dimname
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid, status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. present(dimname)) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dimname), dimid))

         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int32_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_1d (filename, dataname, wdata, dimname, compress)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), intent(in) :: wdata (:)

      CHARACTER(len=*), intent(in), optional :: dimname
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid, status

      write(*,*) trim(dataname), trim(dimname)
      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. present(dimname)) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 
         
         CALL nccheck (nf90_inq_dimid(ncid, trim(dimname), dimid))

         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_logical_1d (filename, dataname, wdata, dimname, compress)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      LOGICAL, intent(in) :: wdata (:)

      CHARACTER(len=*), intent(in)  :: dimname
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER(1), allocatable :: wdata_byte(:)

      allocate(wdata_byte(size(wdata)))
      where(wdata) 
         wdata_byte = 1
      elsewhere
         wdata_byte = 0
      endwhere

      IF (present(compress)) THEN 
         CALL ncio_write_serial_int8_1d (filename, dataname, wdata_byte, dimname, compress)
      ELSE
         CALL ncio_write_serial_int8_1d (filename, dataname, wdata_byte, dimname)
      ENDIF

      deallocate(wdata_byte)

   END SUBROUTINE ncio_write_serial_logical_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int8_2d (filename, dataname, wdata, &
         dim1name, dim2name, compress)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER(1), intent(in) :: wdata (:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(2), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_BYTE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_BYTE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int8_2d
   
   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int16_2d (filename, dataname, wdata, &
         dim1name, dim2name, compress)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER(2), intent(in) :: wdata (:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(2), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_SHORT, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_SHORT, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int16_2d
   
   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int32_2d (filename, dataname, wdata, &
         dim1name, dim2name, compress)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, intent(in) :: wdata (:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(2), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int32_2d
   
   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real4_2d (filename, dataname, wdata, &
         dim1name, dim2name, compress)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(4), intent(in) :: wdata (:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(2), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_FLOAT, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_FLOAT, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real4_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_2d (filename, dataname, wdata, &
         dim1name, dim2name, compress)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), intent(in) :: wdata (:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(2), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_int32_3d (filename, dataname, wdata, &
         dim1name, dim2name, dim3name, compress)

      USE netcdf
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, intent(in) :: wdata (:,:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name, dim3name
      INTEGER, intent(in), optional          :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(3), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name) .and. present(dim3name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim3name), dimid(3)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_INT, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_int32_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_3d (filename, dataname, wdata, &
         dim1name, dim2name, dim3name, compress)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), intent(in) :: wdata (:,:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name, dim3name
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(3), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name) .and. present(dim3name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim3name), dimid(3)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_4d (filename, dataname, wdata, &
         dim1name, dim2name, dim3name, dim4name, compress)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), intent(in) :: wdata (:,:,:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name, dim3name, dim4name
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(4), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name) &
            .and. present(dim3name) .and. present(dim4name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim3name), dimid(3)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim4name), dimid(4)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_4d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_5d (filename, dataname, wdata, &
         dim1name, dim2name, dim3name, dim4name, dim5name, compress)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      REAL(r8), intent(in) :: wdata (:,:,:,:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name, dim3name
      CHARACTER(len=*), intent(in), optional :: dim4name, dim5name
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(5), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name) .and. present(dim3name) &
            .and. present(dim4name) .and. present(dim5name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim3name), dimid(3)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim4name), dimid(4)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim5name), dimid(5)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata) )
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_5d

   !------------------------------
   SUBROUTINE ncio_write_time (filename, dataname, time_component, itime)
      
      USE timemanager
      IMPLICIT NONE

      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      INTEGER, intent(in)  :: time_component(3)
      INTEGER, intent(out) :: itime

      ! Local variables
      INTEGER, allocatable :: time_file(:)
      INTEGER :: ncid, varid, time_id, status
      INTEGER :: timelen, minutes

      minutes = minutes_since_1900 (time_component(1), time_component(2), time_component(3))

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status == NF90_NOERR) THEN
         CALL nccheck( nf90_inq_dimid(ncid, 'time', time_id) )
         CALL nccheck( nf90_inquire_dimension(ncid, time_id, len = timelen) )

         itime = 1
         IF (timelen > 0) THEN
            allocate (time_file (timelen))
            CALL nccheck( nf90_get_var(ncid, varid, time_file) )

            DO while (itime <= timelen)
               IF (minutes == time_file(itime)) exit
               itime = itime + 1
            ENDDO

            deallocate(time_file)
         ENDIF

      ELSE
         status = nf90_inq_dimid(ncid, 'time', time_id)
         IF (status /= NF90_NOERR) THEN
            CALL nccheck( nf90_redef(ncid) )
            CALL nccheck( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_id) )
            CALL nccheck( nf90_enddef(ncid) )
         ENDIF 

         CALL nccheck( nf90_redef(ncid) )
         CALL nccheck( nf90_def_var(ncid, trim(dataname), NF90_INT, (/time_id/), varid) )

         call nccheck( nf90_put_att(ncid, varid, 'long_name', 'time') )
         call nccheck( nf90_put_att(ncid, varid, 'units', 'minutes since 1900-1-1 0:0:0') )
         CALL nccheck( nf90_enddef(ncid) )

         itime = 1
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, minutes, (/itime/)) ) 
      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_time
     

   !----------------------------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_1d_time ( &
         filename, dataname, itime, wdata, &
         dim1name, dim2name, compress)
     
      USE netcdf
      USE precision
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      INTEGER,  intent(in) :: itime
      REAL(r8), intent(in) :: wdata(:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(2), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata, &
         (/1,itime/), (/size(wdata,1),1/)) )

      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_1d_time

   !----------------------------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_2d_time ( &
         filename, dataname, itime, wdata, &
         dim1name, dim2name, dim3name, compress)
     
      USE netcdf
      USE precision
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      INTEGER,  intent(in) :: itime
      REAL(r8), intent(in) :: wdata(:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name, dim3name
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(3), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name) .and. present(dim3name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim3name), dimid(3)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata, &
         (/1,1,itime/), (/size(wdata,1),size(wdata,2),1/)) )

      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_2d_time

   !----------------------------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_3d_time ( &
         filename, dataname, itime, wdata, &
         dim1name, dim2name, dim3name, dim4name, compress)
     
      USE netcdf
      USE precision
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      INTEGER,  intent(in) :: itime
      REAL(r8), intent(in) :: wdata(:,:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name, dim3name, dim4name
      INTEGER, intent(in), optional :: compress
      ! Local variables
      INTEGER :: ncid, varid, dimid(4), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name) &
            .and. present(dim3name) .and. present(dim4name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim3name), dimid(3)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim4name), dimid(4)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata, &
         (/1,1,1,itime/), (/size(wdata,1),size(wdata,2),size(wdata,3),1/)) )

      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_3d_time
   
   !----------------------------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_4d_time ( &
         filename, dataname, itime, wdata, &
         dim1name, dim2name, dim3name, dim4name, dim5name, compress)
     
      USE netcdf
      USE precision
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      INTEGER,  intent(in) :: itime
      REAL(r8), intent(in) :: wdata(:,:,:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name, dim3name
      CHARACTER(len=*), intent(in), optional :: dim4name, dim5name
      INTEGER, intent(in), optional :: compress

      ! Local variables
      INTEGER :: ncid, varid, dimid(5), status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status /= NF90_NOERR) THEN
         IF (.not. (present(dim1name) .and. present(dim2name) &
            .and. present(dim3name) .and. present(dim4name) .and. present(dim5name))) THEN
            write(*,*) 'Warning: no dimension name for ', trim(dataname)
            RETURN
         ENDIF 

         CALL nccheck (nf90_inq_dimid(ncid, trim(dim1name), dimid(1)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim2name), dimid(2)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim3name), dimid(3)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim4name), dimid(4)))
         CALL nccheck (nf90_inq_dimid(ncid, trim(dim5name), dimid(5)))
         
         CALL nccheck (nf90_redef(ncid))
         IF (present(compress)) THEN 
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid, &
               deflate_level = compress))
         ELSE
            CALL nccheck (nf90_def_var(ncid, trim(dataname), NF90_DOUBLE, dimid, varid))
         ENDIF 

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata, &
         (/1,1,1,1,itime/), (/size(wdata,1),size(wdata,2),size(wdata,3),size(wdata,4), 1/)) )

      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_4d_time

END MODULE ncio_serial
