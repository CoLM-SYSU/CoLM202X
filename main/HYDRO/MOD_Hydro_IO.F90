#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_Hydro_IO
   !-----------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !    Read/Write data in lateral hydrological processes.   
   !
   ! Created by Shupeng Zhang, May 2023
   !-----------------------------------------------------------------------

   PUBLIC :: vector_write_basin

   interface vector_read_basin
      MODULE procedure vector_read_basin_real8
      MODULE procedure vector_read_basin_int32
   END interface vector_read_basin

CONTAINS

   ! -------
   subroutine vector_write_basin ( & 
         file_basin, vector, vlen, totalvlen, varname, dimname, data_address, &
         is_hist, itime_in_file, longname, units)

      use MOD_Precision
      use MOD_SPMD_Task
      use MOD_Namelist
      USE MOD_DataType
      USE MOD_NetCDFSerial
      use MOD_Vars_Global, only: spval
      implicit none

      character(len=*), intent(in) :: file_basin
      real(r8),         intent(in) :: vector (:)
      INTEGER,          intent(in) :: vlen
      INTEGER,          intent(in) :: totalvlen
      character(len=*), intent(in) :: varname
      CHARACTER(len=*), intent(in) :: dimname

      TYPE(pointer_int32_1d), intent(in) :: data_address (0:)

      LOGICAL,          intent(in), optional :: is_hist
      integer,          intent(in), optional :: itime_in_file
      character(len=*), intent(in), optional :: longname
      character(len=*), intent(in), optional :: units
      
      ! Local variables
      INTEGER :: iwork, mesg(2), isrc, ndata
      REAL(r8), allocatable :: rcache(:), wdata(:)

      IF (present(is_hist)) THEN
         if (.not. is_hist) return
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      if (p_is_worker) then
         mesg = (/p_iam_glb, vlen/)
         call mpi_send (mesg, 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err) 
         IF (vlen > 0) THEN
            call mpi_send (vector, vlen, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_err) 
         ENDIF
      ENDIF
      
      IF (p_is_master) THEN
         
         allocate (wdata (totalvlen))

         DO iwork = 0, p_np_worker-1
            call mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndata))

               call mpi_recv (rcache, ndata, MPI_REAL8, isrc, &
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
            call ncio_write_serial_time (file_basin, varname, itime_in_file, wdata, &
               dimname, 'time', DEF_HIST_COMPRESS_LEVEL)
         ELSE
            call ncio_write_serial (file_basin, varname, wdata, &
               dimname, DEF_REST_COMPRESS_LEVEL)
         ENDIF

         IF (present(itime_in_file)) THEN
            IF (itime_in_file == 1) then
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
   SUBROUTINE vector_read_basin_real8 ( &
         file_basin, vector, vlen, varname, data_address)

      use MOD_Precision
      use MOD_SPMD_Task
      USE MOD_DataType
      USE MOD_NetCDFSerial
      implicit none

      character(len=*),       intent(in)    :: file_basin
      real(r8),  allocatable, intent(inout) :: vector (:)
      INTEGER,                intent(in)    :: vlen
      character(len=*),       intent(in)    :: varname
      TYPE(pointer_int32_1d), intent(in)    :: data_address (0:)

      ! Local variables
      INTEGER :: iwork, ndata
      REAL(r8), allocatable :: rdata(:), rcache(:)

      IF (p_is_master) THEN
         call ncio_read_serial (file_basin, varname, rdata)
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1
            IF (allocated(data_address(iwork)%val)) THEN
               
               ndata = size(data_address(iwork)%val)
               allocate(rcache (ndata))
               rcache = rdata(data_address(iwork)%val)

               call mpi_send (rcache, ndata, MPI_REAL8, &
                  p_address_worker(iwork), mpi_tag_data, p_comm_glb, p_err) 

               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF
      
      IF (p_is_worker) THEN
         IF (vlen > 0) THEN
            IF (.not. allocated(vector)) allocate(vector(vlen))
            call mpi_recv (vector, vlen, MPI_REAL8, p_root, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
      ENDIF
      
      CALL mpi_barrier (p_comm_glb, p_err)
#else
      IF (.not. allocated(vector)) allocate(vector(vlen))
      vector = rdata(data_address(0)%val)
#endif

      IF (p_is_master) deallocate(rdata)

   END SUBROUTINE vector_read_basin_real8

   ! -----
   SUBROUTINE vector_read_basin_int32 ( &
         file_basin, vector, vlen, varname, data_address)

      use MOD_Precision
      use MOD_SPMD_Task
      USE MOD_DataType
      USE MOD_NetCDFSerial
      implicit none

      character(len=*),       intent(in)    :: file_basin
      integer,   allocatable, intent(inout) :: vector (:)
      INTEGER,                intent(in)    :: vlen
      character(len=*),       intent(in)    :: varname
      TYPE(pointer_int32_1d), intent(in)    :: data_address (0:)

      ! Local variables
      INTEGER :: iwork, ndata
      integer, allocatable :: rdata(:), rcache(:)

      IF (p_is_master) THEN
         call ncio_read_serial (file_basin, varname, rdata)
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1
            IF (allocated(data_address(iwork)%val)) THEN
               
               ndata = size(data_address(iwork)%val)
               allocate(rcache (ndata))
               rcache = rdata(data_address(iwork)%val)

               call mpi_send (rcache, ndata, MPI_INTEGER4, &
                  p_address_worker(iwork), mpi_tag_data, p_comm_glb, p_err) 

               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF
      
      IF (p_is_worker) THEN
         IF (vlen > 0) THEN
            IF (.not. allocated(vector)) allocate(vector(vlen))
            call mpi_recv (vector, vlen, MPI_INTEGER4, p_root, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
      ENDIF
      
      CALL mpi_barrier (p_comm_glb, p_err)
#else
      IF (.not. allocated(vector)) allocate(vector(vlen))
      vector = rdata(data_address(0)%val)
#endif

      IF (p_is_master) deallocate(rdata)

   END SUBROUTINE vector_read_basin_int32

END MODULE MOD_Hydro_IO
#endif
