#include <define.h>

MODULE ncio_serial

   USE netcdf
   IMPLICIT NONE

   ! PUBLIC subroutines

   PUBLIC :: ncio_create_file

   PUBLIC :: ncio_put_attr

   PUBLIC :: ncio_var_exist
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
      MODULE procedure ncio_read_bcast_serial_logical_1d 
   END INTERFACE ncio_read_bcast_serial

   PUBLIC :: ncio_define_dimension

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
   SUBROUTINE ncio_put_attr (filename, varname, attrname, attrval)

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

   END SUBROUTINE ncio_put_attr


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
      ENDIF
      
#ifdef USEMPI
      CALL mpi_bcast (vlen, 1, MPI_INTEGER, p_root, p_comm_glb, p_err)
      IF (.not. p_is_master)  allocate (rdata (vlen))
      CALL mpi_bcast (rdata, vlen, MPI_LOGICAL, p_root, p_comm_glb, p_err)
#endif

   END SUBROUTINE ncio_read_bcast_serial_logical_1d

   ! -------------------------------
   SUBROUTINE ncio_define_dimension (filename, dimname, dimlen)

      USE netcdf
      USE precision
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dimname
      INTEGER, intent(in) :: dimlen

      ! Local variables
      INTEGER :: ncid, dimid, status

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )

      status = nf90_inq_dimid(ncid, trim(dimname), dimid)
      IF (status /= NF90_NOERR) THEN
         CALL nccheck (nf90_redef(ncid))
         IF (dimlen == 0) THEN
            CALL nccheck( nf90_def_dim(ncid, trim(dimname), NF90_UNLIMITED, dimid) )
         ELSE
            CALL nccheck( nf90_def_dim(ncid, trim(dimname), dimlen, dimid) )
         ENDIF 
         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_define_dimension

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

   END SUBROUTINE ncio_write_serial_logical_1d

   !---------------------------------------------------------
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

      IMPLICIT NONE

      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      INTEGER, intent(in)  :: time_component(3)
      INTEGER, intent(out) :: itime

      ! Local variables
      INTEGER, allocatable :: time_file(:,:)
      INTEGER :: ncid, varid, time_id, component_id, status
      INTEGER :: timelen

      CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
      status = nf90_inq_varid(ncid, trim(dataname), varid)
      IF (status == NF90_NOERR) THEN
         CALL nccheck( nf90_inq_dimid(ncid, 'time', time_id) )
         CALL nccheck( nf90_inquire_dimension(ncid, time_id, len = timelen) )

         itime = 1
         IF (timelen > 0) THEN
            allocate (time_file (3, timelen))
            CALL nccheck( nf90_get_var(ncid, varid, time_file) )

            DO while (itime <= timelen)
               IF (all(time_component == time_file(:,itime))) exit
               itime = itime + 1
            ENDDO

            deallocate(time_file)
         ENDIF

         CALL nccheck( nf90_put_var(ncid, varid, time_component, (/1,itime/), (/3,1/)) ) 

      ELSE
         status = nf90_inq_dimid(ncid, 'time', time_id)
         IF (status /= NF90_NOERR) THEN
            CALL nccheck( nf90_redef(ncid) )
            CALL nccheck( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_id) )
            CALL nccheck( nf90_enddef(ncid) )
         ENDIF 

         status = nf90_inq_dimid(ncid, 'component', component_id)
         IF (status /= NF90_NOERR) THEN
            CALL nccheck( nf90_redef(ncid) )
            CALL nccheck( nf90_def_dim(ncid, 'component', 3, component_id) )
            CALL nccheck( nf90_enddef(ncid) )
         ENDIF 
         
         CALL nccheck( nf90_redef(ncid) )
         CALL nccheck( nf90_def_var(ncid, trim(dataname), NF90_INT, &
            (/component_id, time_id/), varid) )
         CALL nccheck( nf90_enddef(ncid) )
        
         itime = 1
         CALL nccheck( nf90_put_var(ncid, varid, time_component, (/1,itime/), (/3,1/)) ) 
      ENDIF 

      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_time

   !----------------------------------------------------------------------------
   SUBROUTINE ncio_write_serial_real8_2d_time ( &
         filename, dataname, itime, wdata, &
         dim1name, dim2name, dim3name, compress, longname, units)
     
      USE netcdf
      USE precision
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      INTEGER,  intent(in) :: itime
      REAL(r8), intent(in) :: wdata(:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name, dim3name
      INTEGER, intent(in), optional :: compress
      character (len=*), intent(in),optional :: longname
      character (len=*), intent(in),optional :: units

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

         IF (present(longname)) THEN
            CALL nccheck (nf90_put_att(ncid, varid, 'long_name', trim(longname)))
         ENDIF
         
         IF (present(units)) THEN
            CALL nccheck (nf90_put_att(ncid, varid, 'units', trim(units)))
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
         dim1name, dim2name, dim3name, dim4name, compress, longname, units)
     
      USE netcdf
      USE precision
      IMPLICIT NONE
      
      CHARACTER (len=*), intent(in) :: filename
      CHARACTER (len=*), intent(in) :: dataname
      INTEGER,  intent(in) :: itime
      REAL(r8), intent(in) :: wdata(:,:,:)
      
      CHARACTER(len=*), intent(in), optional :: dim1name, dim2name, dim3name, dim4name
      INTEGER, intent(in), optional :: compress
      character (len=*), intent(in),optional :: longname
      character (len=*), intent(in),optional :: units
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

         IF (present(longname)) THEN
            CALL nccheck (nf90_put_att(ncid, varid, 'long_name', trim(longname)))
         ENDIF
         
         IF (present(units)) THEN
            CALL nccheck (nf90_put_att(ncid, varid, 'units', trim(units)))
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
         dim1name, dim2name, dim3name, dim4name, dim5name, compress, longname, units)
     
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
      character (len=*), intent(in),optional :: longname
      character (len=*), intent(in),optional :: units

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

         IF (present(longname)) THEN
            CALL nccheck (nf90_put_att(ncid, varid, 'long_name', trim(longname)))
         ENDIF
         
         IF (present(units)) THEN
            CALL nccheck (nf90_put_att(ncid, varid, 'units', trim(units)))
         ENDIF

         CALL nccheck (nf90_enddef(ncid))
      ENDIF 

      CALL nccheck( nf90_put_var(ncid, varid, wdata, &
         (/1,1,1,1,itime/), (/size(wdata,1),size(wdata,2),size(wdata,3),size(wdata,4), 1/)) )

      CALL nccheck( nf90_close(ncid) )

   END SUBROUTINE ncio_write_serial_real8_4d_time

END MODULE ncio_serial
