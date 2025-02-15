#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_IO
!-----------------------------------------------------------------------
! DESCRIPTION:
!
!    Read/Write data in lateral hydrological processes.   
!
! Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   PUBLIC :: vector_write_basin
   PUBLIC :: vector_read_basin

CONTAINS

   ! -------
   SUBROUTINE vector_write_basin ( & 
         file_basin, vector, vlen, totalvlen, varname, dimname, data_address, &
         is_hist, itime_in_file, longname, units)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_DataType
   USE MOD_NetCDFSerial
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   character(len=*), intent(in) :: file_basin
   real(r8),         intent(in) :: vector (:)
   integer,          intent(in) :: vlen
   integer,          intent(in) :: totalvlen
   character(len=*), intent(in) :: varname
   character(len=*), intent(in) :: dimname

   type(pointer_int32_1d), intent(in) :: data_address (0:)

   logical,          intent(in), optional :: is_hist
   integer,          intent(in), optional :: itime_in_file
   character(len=*), intent(in), optional :: longname
   character(len=*), intent(in), optional :: units
   
   ! Local variables
   integer :: iwork, mesg(2), isrc, ndata
   real(r8), allocatable :: rcache(:), wdata(:)

      IF (present(is_hist)) THEN
         IF (.not. is_hist) RETURN
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
         
         allocate (wdata (totalvlen))

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
      allocate (wdata (totalvlen))
      wdata(data_address(0)%val) = vector
#endif
      
      IF (p_is_master) THEN

         IF (present(itime_in_file)) THEN
            CALL ncio_write_serial_time (file_basin, varname, itime_in_file, wdata, &
               dimname, 'time', DEF_HIST_CompressLevel)
         ELSE
            CALL ncio_write_serial (file_basin, varname, wdata, &
               dimname, DEF_REST_CompressLevel)
         ENDIF

         IF (present(itime_in_file)) THEN
            IF (itime_in_file == 1) THEN
               IF (present(longname)) THEN
                  CALL ncio_put_attr (file_basin, varname, 'long_name', longname)
               ENDIF
               IF (present(units)) THEN
                  CALL ncio_put_attr (file_basin, varname, 'units', units)
               ENDIF
               CALL ncio_put_attr (file_basin, varname, 'missing_value', spval)
            ENDIF
         ENDIF

         deallocate (wdata)

      ENDIF
         
   END SUBROUTINE vector_write_basin 

   ! -----
   SUBROUTINE vector_read_basin ( &
         file_basin, vector, vlen, varname, data_address)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_DataType
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   character(len=*),       intent(in)    :: file_basin
   real(r8),  allocatable, intent(inout) :: vector (:)
   integer,                intent(in)    :: vlen
   character(len=*),       intent(in)    :: varname
   type(pointer_int32_1d), intent(in)    :: data_address (0:)

   ! Local variables
   integer :: iwork, ndata
   real(r8), allocatable :: rdata(:), rcache(:)

      IF (p_is_master) THEN
         CALL ncio_read_serial (file_basin, varname, rdata)
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

   END SUBROUTINE vector_read_basin

END MODULE MOD_Catch_IO
#endif
