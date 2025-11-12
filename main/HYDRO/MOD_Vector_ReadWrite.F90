#include <define.h>

MODULE MOD_Vector_ReadWrite
!-----------------------------------------------------------------------
! DESCRIPTION:
!
!    Read/Write data in vector form.
!
! Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   PUBLIC :: vector_gather_and_write
   PUBLIC :: vector_gather_map2grid_and_write
   PUBLIC :: vector_read_and_scatter

CONTAINS

   ! -------
   SUBROUTINE vector_gather_to_master ( &
         vector, vlen, totalvlen, data_address, wdata)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_DataType
   IMPLICIT NONE

   real(r8), intent(in) :: vector (:)
   integer,  intent(in) :: vlen
   integer,  intent(in) :: totalvlen

   type(pointer_int32_1d), intent(in) :: data_address (0:)

   real(r8), allocatable,  intent(inout) :: wdata (:)

   ! Local variables
   integer :: iwork, mesg(2), isrc, ndata
   real(r8), allocatable :: rcache(:)

      IF (totalvlen <= 0) RETURN

      IF (p_is_master) THEN
         allocate (wdata (totalvlen))
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         mesg = (/p_iam_glb, vlen/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (vlen > 0) THEN
            CALL mpi_send (vector, vlen, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
      ENDIF

      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndata))

               CALL mpi_recv (rcache, ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

               wdata(data_address(p_itis_worker(isrc))%val) = rcache

               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      wdata(data_address(0)%val) = vector
#endif

   END SUBROUTINE vector_gather_to_master

   ! -------
   SUBROUTINE vector_gather_and_write ( vector, vlen, totalvlen, data_address, &
         fileout, varname, dimname, itime_in_file, longname, units)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_DataType
   USE MOD_NetCDFSerial
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(in) :: vector (:)
   integer,  intent(in) :: vlen
   integer,  intent(in) :: totalvlen

   type(pointer_int32_1d), intent(in) :: data_address (0:)

   character(len=*), intent(in) :: fileout
   character(len=*), intent(in) :: varname
   character(len=*), intent(in) :: dimname

   integer,          intent(in), optional :: itime_in_file
   character(len=*), intent(in), optional :: longname
   character(len=*), intent(in), optional :: units

   ! Local variables
   real(r8), allocatable :: wdata(:)
   logical :: write_attr


      CALL vector_gather_to_master (vector, vlen, totalvlen, data_address, wdata)

      IF (p_is_master) THEN

         IF (present(itime_in_file)) THEN
            CALL ncio_write_serial_time (fileout, varname, itime_in_file, wdata, &
               dimname, 'time', DEF_HIST_CompressLevel)
         ELSE
            CALL ncio_write_serial (fileout, varname, wdata, &
               dimname, DEF_REST_CompressLevel)
         ENDIF

         IF (present(itime_in_file)) THEN
            write_attr = itime_in_file <= 1
         ELSE
            write_attr = .true.
         ENDIF

         IF (write_attr) THEN
            CALL ncio_put_attr (fileout, varname, 'missing_value', spval)
            IF (present(longname)) CALL ncio_put_attr (fileout, varname, 'long_name', longname)
            IF (present(units   )) CALL ncio_put_attr (fileout, varname, 'units',     units   )
         ENDIF

         deallocate (wdata)

      ENDIF

   END SUBROUTINE vector_gather_and_write

   ! -------
   SUBROUTINE vector_gather_map2grid_and_write ( &
         vector,  vlen,    totalvlen, data_address, nlon, x_vec,   nlat, y_vec,   &
         fileout, varname, lon_name,  lat_name,     itime_in_file, longname, units)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_DataType
   USE MOD_NetCDFSerial
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   character(len=*), intent(in) :: fileout
   real(r8),         intent(in) :: vector (:)
   integer,          intent(in) :: vlen
   integer,          intent(in) :: totalvlen
   character(len=*), intent(in) :: varname
   character(len=*), intent(in) :: lon_name, lat_name

   type(pointer_int32_1d), intent(in) :: data_address (0:)

   integer, intent(in) :: nlon, x_vec (:)
   integer, intent(in) :: nlat, y_vec (:)

   integer,          intent(in), optional :: itime_in_file
   character(len=*), intent(in), optional :: longname
   character(len=*), intent(in), optional :: units

   ! Local variables
   integer :: i
   real(r8), allocatable :: wdata(:), wdata2d(:,:)
   logical :: write_attr

      CALL vector_gather_to_master (vector, vlen, totalvlen, data_address, wdata)

      IF (p_is_master) THEN

         allocate (wdata2d (nlon,nlat))
         wdata2d(:,:) = spval

         DO i = 1, totalvlen
            wdata2d(x_vec(i),y_vec(i)) = wdata(i)
         ENDDO

         IF (present(itime_in_file)) THEN
            CALL ncio_write_serial_time (fileout, varname, itime_in_file, wdata2d, &
               lon_name, lat_name, 'time', DEF_HIST_CompressLevel)
         ELSE
            CALL ncio_write_serial (fileout, varname, wdata2d, &
               lon_name, lat_name, DEF_REST_CompressLevel)
         ENDIF

         IF (present(itime_in_file)) THEN
            write_attr = itime_in_file == 1
         ELSE
            write_attr = .true.
         ENDIF

         IF (write_attr) THEN
            CALL ncio_put_attr (fileout, varname, 'missing_value', spval)
            IF (present(longname)) CALL ncio_put_attr (fileout, varname, 'long_name', longname)
            IF (present(units   )) CALL ncio_put_attr (fileout, varname, 'units',     units   )
         ENDIF

         deallocate (wdata  )
         deallocate (wdata2d)

      ENDIF

   END SUBROUTINE vector_gather_map2grid_and_write

   ! -----
   SUBROUTINE vector_read_and_scatter ( &
         filein, vector, vlen, varname, data_address)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_DataType
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   character(len=*),       intent(in)    :: filein
   real(r8),  allocatable, intent(inout) :: vector (:)
   integer,                intent(in)    :: vlen
   character(len=*),       intent(in)    :: varname
   type(pointer_int32_1d), intent(in)    :: data_address (0:)

   ! Local variables
   integer :: iwork, ndata
   real(r8), allocatable :: rdata(:), rcache(:)

      IF (p_is_master) THEN
         CALL ncio_read_serial (filein, varname, rdata)
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1
            IF (allocated(data_address(iwork)%val)) THEN

               ndata = size(data_address(iwork)%val)
               allocate(rcache (ndata))
               rcache = rdata(data_address(iwork)%val)

               CALL mpi_send (rcache, ndata, MPI_REAL8, &
                  p_address_worker(iwork), mpi_tag_data, p_comm_glb, p_err)

               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_worker) THEN
         IF (vlen > 0) THEN
            IF (.not. allocated(vector)) allocate(vector(vlen))
            CALL mpi_recv (vector, vlen, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      IF (.not. allocated(vector)) allocate(vector(vlen))
      vector = rdata(data_address(0)%val)
#endif

      IF (p_is_master) deallocate(rdata)

   END SUBROUTINE vector_read_and_scatter

END MODULE MOD_Vector_ReadWrite
