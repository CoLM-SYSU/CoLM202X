#include <define.h>

#ifdef USEMPI
MODULE MOD_HistWriteBack
!----------------------------------------------------------------------------
! DESCRIPTION:
!
!     Write out data to history files by a dedicated process.
!
! Author: Shupeng Zhang, 11/2023
!----------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial

   ! type of send buffer
   type :: HistSendBufferType
      integer  :: dataid
      integer  :: datatag
      integer  :: sendreqs (3)
      integer  :: sendint4 (5)
      character(len=256)    :: sendchar (9)
      real(r8), allocatable :: senddata (:)
      type(HistSendBufferType), pointer :: next
   END type HistSendBufferType

   ! Sending Variables
   type(HistSendBufferType), pointer :: HistSendBuffer
   type(HistSendBufferType), pointer :: LastSendBuffer

   ! type of times
   type :: timenodetype
      character(len=256) :: filename
      character(len=256) :: timename
      integer :: time(3)
      integer :: req (3)
      type(timenodetype), pointer :: next
   END type timenodetype

   ! time nodes
   integer :: dataid_zero = 0
   integer :: req_zero
   type(timenodetype), pointer :: timenodes, lasttime

   ! dimension information
   logical :: SDimInited = .false.
   ! 1: grid-based; 2: catchment based; 3: unstructered
   integer :: SDimType 

   ! 1: grid-based
   integer :: nGridData, nxGridSeg, nyGridSeg
   integer, allocatable :: xGridDsp(:), xGridCnt(:)
   integer, allocatable :: yGridDsp(:), yGridCnt(:)

   integer :: nlat, nlon
   real(r8), allocatable :: lat_c(:), lat_s(:), lat_n(:)
   real(r8), allocatable :: lon_c(:), lon_w(:), lon_e(:)

   ! 2: catchment based; 3: unstructured 
   ! integer :: SDimLength 
   ! integer*8, allocatable :: vindex1(:)
   ! integer,   allocatable :: vindex2(:)

   ! Memory limits
   integer*8, parameter :: MaxHistMemSize  = 8589934592_8 ! 8*1024^3 
   integer*8, parameter :: MaxHistMesgSize = 8388608_8    ! 8*1024^2 
   
   integer*8 :: TotalMemSize = 0

   integer :: itime_in_file

   ! tags
   integer, parameter :: tag_next = 1
   integer, parameter :: tag_time = 2
   integer, parameter :: tag_dims = 3

CONTAINS

   ! -----
   SUBROUTINE hist_writeback_daemon ()

   USE MOD_Namelist, only : DEF_HIST_FREQ
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   ! Local Variables
   integer :: dataid, tag
   integer :: time(3), ndims, ndim1, ndim2, dimlens(4), compress
   integer :: i, idata, isrc, ixseg, iyseg, xdsp, ydsp, xcnt, ycnt

   integer               :: recvint4 (5)
   character(len=256)    :: recvchar (9)
   real(r8), allocatable :: datathis (:)
   
   character(len=256)    :: filename, dataname, longname, units
   character(len=256)    :: dim1name, dim2name, dim3name, dim4name, dim5name
   logical               :: fexists

   real(r8), allocatable :: wdata1d(:),  wdata2d(:,:), wdata3d(:,:,:), wdata4d(:,:,:,:)


      DO WHILE (.true.)

         CALL mpi_recv (dataid, 1, MPI_INTEGER, &
            MPI_ANY_SOURCE, tag_next, p_comm_glb_plus, p_stat, p_err)

         IF (dataid < 0) THEN

            EXIT

         ELSEIF (dataid == 0) THEN
            
            CALL mpi_recv (filename, 256, MPI_CHARACTER, &
               MPI_ANY_SOURCE, tag_time, p_comm_glb_plus, p_stat, p_err)

            CALL mpi_recv (dataname, 256, MPI_CHARACTER, &
               MPI_ANY_SOURCE, tag_time, p_comm_glb_plus, p_stat, p_err)
            
            CALL mpi_recv (time, 3, MPI_INTEGER, &
               MPI_ANY_SOURCE, tag_time, p_comm_glb_plus, p_stat, p_err)

            IF (.not. SDimInited) THEN
               
               CALL mpi_recv (SDimType, 1, MPI_INTEGER, &
                  MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)
               
               IF (SDimType == 1) THEN

                  CALL mpi_recv (nGridData, 1, MPI_INTEGER, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)
                  CALL mpi_recv (nxGridSeg, 1, MPI_INTEGER, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)
                  CALL mpi_recv (nyGridSeg, 1, MPI_INTEGER, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)
         
                  allocate (xGridDsp (nxGridSeg))
                  allocate (xGridCnt (nxGridSeg))
                  allocate (yGridDsp (nyGridSeg))
                  allocate (yGridCnt (nyGridSeg))
                  
                  CALL mpi_recv (xGridDsp, nxGridSeg, MPI_INTEGER, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)
                  CALL mpi_recv (xGridCnt, nxGridSeg, MPI_INTEGER, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)
                  CALL mpi_recv (yGridDsp, nyGridSeg, MPI_INTEGER, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)
                  CALL mpi_recv (yGridCnt, nyGridSeg, MPI_INTEGER, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)

                  CALL mpi_recv (nlat, 1, MPI_INTEGER, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)

                  allocate(lat_c(nlat))
                  allocate(lat_s(nlat))
                  allocate(lat_n(nlat))

                  CALL mpi_recv (lat_c, nlat, MPI_REAL8, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)
                  CALL mpi_recv (lat_s, nlat, MPI_REAL8, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)
                  CALL mpi_recv (lat_n, nlat, MPI_REAL8, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)

                  CALL mpi_recv (nlon, 1, MPI_INTEGER, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)

                  allocate(lon_c(nlon))
                  allocate(lon_w(nlon))
                  allocate(lon_e(nlon))
               
                  CALL mpi_recv (lon_c, nlon, MPI_REAL8, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)
                  CALL mpi_recv (lon_w, nlon, MPI_REAL8, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)
                  CALL mpi_recv (lon_e, nlon, MPI_REAL8, &
                     MPI_ANY_SOURCE, tag_dims, p_comm_glb_plus, p_stat, p_err)

               ENDIF

               SDimInited = .true.

            ENDIF
           
            inquire (file=filename, exist=fexists)
            IF (.not. fexists) THEN
               
               CALL ncio_create_file (trim(filename))
               
               CALL ncio_define_dimension(filename, 'time', 0)
               
               IF (SDimType == 1) THEN
                  
                  CALL ncio_define_dimension(filename, 'lat', nlat)
                  CALL ncio_define_dimension(filename, 'lon', nlon)

                  CALL ncio_write_serial (filename, 'lat',   lat_c, 'lat')
                  CALL ncio_write_serial (filename, 'lon',   lon_c, 'lon')
                  CALL ncio_write_serial (filename, 'lat_s', lat_s, 'lat')
                  CALL ncio_write_serial (filename, 'lat_n', lat_n, 'lat')
                  CALL ncio_write_serial (filename, 'lon_w', lon_w, 'lon')
                  CALL ncio_write_serial (filename, 'lon_e', lon_e, 'lon')

                  CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
                  CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')
                  CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
                  CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')

               ENDIF

               CALL ncio_write_colm_dimension (filename)

            ENDIF
               
            CALL ncio_write_time (filename, dataname, time, itime_in_file, DEF_HIST_FREQ)

         ELSE

            !--------------------------------
            ! reveive and write history data.
            !--------------------------------

            ! (1) data header
            tag = dataid*10
            CALL mpi_recv (recvint4(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, tag, p_comm_glb_plus, p_stat, p_err)

            ndims    = recvint4(1)
            compress = recvint4(2)

            CALL mpi_recv (recvchar(1:9), 256*9, MPI_CHARACTER, &
               MPI_ANY_SOURCE, tag, p_comm_glb_plus, p_stat, p_err)

            filename = recvchar(1)
            dataname = recvchar(2) 
            dim1name = recvchar(3)
            dim2name = recvchar(4)
            dim3name = recvchar(5)
            dim4name = recvchar(6)
            dim5name = recvchar(7)
            longname = recvchar(8)
            units    = recvchar(9)

            ! (2) data 
            tag = dataid*10+1

            IF (SDimType == 1) THEN
               DO idata = 1, nGridData

                  CALL mpi_recv (recvint4(1:5), 5, MPI_INTEGER, MPI_ANY_SOURCE, &
                     tag, p_comm_glb_plus, p_stat, p_err)

                  isrc  = recvint4(1)
                  ixseg = recvint4(2)
                  iyseg = recvint4(3)
                  ndim1 = recvint4(4)
                  ndim2 = recvint4(5)

                  xdsp = xGridDsp(ixseg) 
                  ydsp = yGridDsp(iyseg) 
                  xcnt = xGridCnt(ixseg)
                  ycnt = yGridCnt(iyseg)

                  SELECTCASE (ndims)
                  CASE (2)

                     dimlens = (/nlon, nlat, 0, 0/)

                     IF (.not. allocated(wdata2d)) THEN
                        allocate (wdata2d (nlon,nlat))
                     ENDIF

                     allocate (datathis(xcnt*ycnt))
                     CALL mpi_recv (datathis, xcnt*ycnt, MPI_REAL8, &
                        isrc, tag, p_comm_glb_plus, p_stat, p_err)

                     wdata2d(xdsp+1:xdsp+xcnt, ydsp+1:ydsp+ycnt) = &
                        reshape(datathis,(/xcnt,ycnt/))

                  CASE (3)

                     dimlens = (/ndim1, nlon, nlat, 0/)

                     IF (.not. allocated(wdata3d)) THEN
                        allocate (wdata3d (ndim1,nlon,nlat))
                     ENDIF

                     allocate (datathis(ndim1*xcnt*ycnt))
                     CALL mpi_recv (datathis, ndim1*xcnt*ycnt, MPI_REAL8, &
                        isrc, tag, p_comm_glb_plus, p_stat, p_err)

                     wdata3d(:,xdsp+1:xdsp+xcnt, ydsp+1:ydsp+ycnt) = &
                        reshape(datathis,(/ndim1,xcnt,ycnt/))

                  CASE (4)

                     dimlens = (/ndim1, ndim2, nlon, nlat/)

                     IF (.not. allocated(wdata4d)) THEN
                        allocate (wdata4d (ndim1,ndim2,nlon,nlat))
                     ENDIF

                     allocate (datathis(ndim1*ndim2*xcnt*ycnt))
                     CALL mpi_recv (datathis, ndim1*ndim2*xcnt*ycnt, MPI_REAL8, &
                        isrc, tag, p_comm_glb_plus, p_stat, p_err)

                     wdata4d(:,:,xdsp+1:xdsp+xcnt, ydsp+1:ydsp+ycnt) = &
                        reshape(datathis,(/ndim1,ndim2,xcnt,ycnt/))

                  ENDSELECT

                  deallocate (datathis)

               ENDDO 

            ENDIF

            IF (ndims >= 1) CALL ncio_define_dimension (filename, dim1name, dimlens(1))
            IF (ndims >= 2) CALL ncio_define_dimension (filename, dim2name, dimlens(2))
            IF (ndims >= 3) CALL ncio_define_dimension (filename, dim3name, dimlens(3))
            IF (ndims >= 4) CALL ncio_define_dimension (filename, dim4name, dimlens(4))

            SELECTCASE (ndims)
            CASE (1)

               CALL ncio_write_serial_time (filename, dataname, itime_in_file, wdata1d, &
                  dim1name, dim2name, compress)

               deallocate(wdata1d)
            CASE (2)

               IF (.not. &
                  ((trim(dataname) == 'landarea') .or. (trim(dataname) == 'landfraction'))) THEN

                  CALL ncio_write_serial_time (filename, dataname, itime_in_file, wdata2d, &
                     dim1name, dim2name, dim3name, compress)

               ELSEIF (itime_in_file == 1) THEN

                  CALL ncio_write_serial (filename, dataname, wdata2d, dim1name, dim2name, compress)

               ENDIF

               deallocate(wdata2d)
            CASE (3)

               CALL ncio_write_serial_time (filename, dataname, itime_in_file, wdata3d, &
                  dim1name, dim2name, dim3name, dim4name, compress)

               deallocate(wdata3d)
            CASE (4)

               CALL ncio_write_serial_time (filename, dataname, itime_in_file, wdata4d, &
                  dim1name, dim2name, dim3name, dim4name, dim5name, compress)

               deallocate(wdata4d)
            ENDSELECT

            IF (itime_in_file == 1) THEN
               CALL ncio_put_attr (filename, dataname, 'long_name', longname)
               CALL ncio_put_attr (filename, dataname, 'units', units)
               CALL ncio_put_attr (filename, dataname, 'missing_value', spval)
            ENDIF

            write(*,'(3A,I0,2A)') 'HIST WriteBack: ', trim(basename(filename)), &
               ' (time ', itime_in_file, '): ', trim(dataname)

         ENDIF
      
      ENDDO

   END SUBROUTINE hist_writeback_daemon

   ! -----
   SUBROUTINE hist_writeback_latlon_time (filename, timename, time, HistConcat)

   USE MOD_Namelist
   USE MOD_Grid
   IMPLICIT NONE
   
   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: timename
   integer, intent(in)  :: time(3)
   type(grid_concat_type), intent(in) :: HistConcat

   ! Local Variables
   integer :: i

      CALL hist_writeback_append_timenodes (filename, timename, time)

      CALL mpi_isend (dataid_zero, 1, MPI_INTEGER, &
         p_address_writeback, tag_next, p_comm_glb_plus, req_zero, p_err)

      CALL mpi_isend (lasttime%filename, 256, MPI_CHARACTER, &
         p_address_writeback, tag_time, p_comm_glb_plus, lasttime%req(1), p_err)
      
      CALL mpi_isend (lasttime%timename, 256, MPI_CHARACTER, &
         p_address_writeback, tag_time, p_comm_glb_plus, lasttime%req(2), p_err)

      CALL mpi_isend (lasttime%time, 3, MPI_INTEGER, &
         p_address_writeback, tag_time, p_comm_glb_plus, lasttime%req(3), p_err)


      IF (.not. SDimInited) THEN

         SDimType = 1
         CALL mpi_send (SDimType, 1, MPI_INTEGER, p_address_writeback, tag_dims, p_comm_glb_plus, p_err)

         nGridData = HistConcat%ndatablk
         nxGridSeg = HistConcat%nxseg
         nyGridSeg = HistConcat%nyseg  

         allocate (xGridDsp (nxGridSeg))
         allocate (xGridCnt (nxGridSeg))
         allocate (yGridDsp (nyGridSeg))
         allocate (yGridCnt (nyGridSeg))

         DO i = 1, nxGridSeg
            xGridDsp(i) = HistConcat%xsegs(i)%gdsp
            xGridCnt(i) = HistConcat%xsegs(i)%cnt
         ENDDO

         DO i = 1, nyGridSeg
            yGridDsp(i) = HistConcat%ysegs(i)%gdsp
            yGridCnt(i) = HistConcat%ysegs(i)%cnt
         ENDDO

         CALL mpi_send (nGridData, 1, MPI_INTEGER, p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (nxGridSeg, 1, MPI_INTEGER, p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (nyGridSeg, 1, MPI_INTEGER, p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         
         CALL mpi_send (xGridDsp, nxGridSeg, MPI_INTEGER, p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (xGridCnt, nxGridSeg, MPI_INTEGER, p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (yGridDsp, nyGridSeg, MPI_INTEGER, p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (yGridCnt, nyGridSeg, MPI_INTEGER, p_address_writeback, tag_dims, p_comm_glb_plus, p_err)

         nlat = HistConcat%ginfo%nlat
         nlon = HistConcat%ginfo%nlon
         allocate(lat_c(nlat));  lat_c = HistConcat%ginfo%lat_c
         allocate(lat_s(nlat));  lat_s = HistConcat%ginfo%lat_s
         allocate(lat_n(nlat));  lat_n = HistConcat%ginfo%lat_n
         allocate(lon_c(nlon));  lon_c = HistConcat%ginfo%lon_c
         allocate(lon_w(nlon));  lon_w = HistConcat%ginfo%lon_w
         allocate(lon_e(nlon));  lon_e = HistConcat%ginfo%lon_e

         CALL mpi_send (nlat,     1, MPI_INTEGER, p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (lat_c, nlat, MPI_REAL8,   p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (lat_s, nlat, MPI_REAL8,   p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (lat_n, nlat, MPI_REAL8,   p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (nlon,     1, MPI_INTEGER, p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (lon_c, nlon, MPI_REAL8,   p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (lon_w, nlon, MPI_REAL8,   p_address_writeback, tag_dims, p_comm_glb_plus, p_err)
         CALL mpi_send (lon_e, nlon, MPI_REAL8,   p_address_writeback, tag_dims, p_comm_glb_plus, p_err)

         SDimInited = .true.

      ENDIF

      CALL hist_writeback_clean_timenodes ()

   END SUBROUTINE hist_writeback_latlon_time

   ! -----
   SUBROUTINE hist_writeback_append_timenodes (filename, timename, time)

   IMPLICIT NONE

   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: timename
   integer, intent(in)  :: time(3)

      IF (.not. associated(timenodes)) THEN
         allocate (timenodes)
         lasttime => timenodes
      ELSE
         allocate (lasttime%next)
         lasttime => lasttime%next
      ENDIF

      lasttime%filename = filename
      lasttime%timename = timename
      lasttime%time = time
      lasttime%next => null()

   END SUBROUTINE hist_writeback_append_timenodes 

   ! -----
   SUBROUTINE hist_writeback_clean_timenodes

   IMPLICIT NONE
   
   ! Local Variables
   logical :: senddone
   integer :: stat(MPI_STATUS_SIZE,3)
   type(timenodetype), pointer :: tempnode

      
      DO WHILE (associated(timenodes%next))

         CALL MPI_TestAll (3, timenodes%req, senddone, stat, p_err)
         
         IF (senddone) THEN
            tempnode  => timenodes
            timenodes => timenodes%next
            deallocate(tempnode)
         ELSE
            EXIT
         ENDIF
      ENDDO

   END SUBROUTINE hist_writeback_clean_timenodes

   ! -----
   SUBROUTINE hist_writeback_var_header (dataid, filename, dataname, &
         ndims,    dim1name, dim2name, dim3name, dim4name, dim5name, &
         compress, longname, units)
      
   IMPLICIT NONE

   integer,          intent(in) :: dataid
   character(len=*), intent(in) :: filename, dataname
   integer,          intent(in) :: ndims
   character(len=*), intent(in) :: dim1name, dim2name, dim3name, dim4name, dim5name
   integer,          intent(in) :: compress
   character(len=*), intent(in) :: longname, units

   ! Local Variables
   logical :: senddone
   integer :: sendstat(MPI_STATUS_SIZE,3)
   type(HistSendBufferType), pointer :: TempSendBuffer

      ! append sending buffer
      IF (.not. associated(HistSendBuffer)) THEN
         allocate (HistSendBuffer)
         LastSendBuffer => HistSendBuffer
      ELSE
         allocate (LastSendBuffer%next)
         LastSendBuffer => LastSendBuffer%next
      ENDIF

      LastSendBuffer%next => null()
      
      ! clean sending buffer and free memory
      DO WHILE (associated(HistSendBuffer%next))

         CALL MPI_Testall (3, HistSendBuffer%sendreqs, senddone, sendstat, p_err)
         
         IF (senddone) THEN
            
            TempSendBuffer => HistSendBuffer
            HistSendBuffer => HistSendBuffer%next

            IF (allocated(TempSendBuffer%senddata)) THEN
               deallocate (TempSendBuffer%senddata)
            ENDIF
            deallocate (TempSendBuffer)
         ELSE
            EXIT
         ENDIF
      ENDDO

      LastSendBuffer%dataid  = dataid
      LastSendBuffer%datatag = dataid*10
      
      LastSendBuffer%sendint4(1:2) = (/ndims, compress/)

      LastSendBuffer%sendchar(1) = filename
      LastSendBuffer%sendchar(2) = dataname
      LastSendBuffer%sendchar(3) = dim1name
      LastSendBuffer%sendchar(4) = dim2name
      LastSendBuffer%sendchar(5) = dim3name
      LastSendBuffer%sendchar(6) = dim4name
      LastSendBuffer%sendchar(7) = dim5name
      LastSendBuffer%sendchar(8) = longname
      LastSendBuffer%sendchar(9) = units
      
      CALL mpi_isend (LastSendBuffer%dataid, 1, MPI_INTEGER, &
         p_address_writeback, tag_next, p_comm_glb_plus, LastSendBuffer%sendreqs(1), p_err)
      
      CALL mpi_isend (LastSendBuffer%sendint4(1:2), 2, MPI_INTEGER, &
         p_address_writeback, LastSendBuffer%datatag, p_comm_glb_plus, LastSendBuffer%sendreqs(2), p_err)

      CALL mpi_isend (LastSendBuffer%sendchar, 256*9, MPI_CHARACTER, &
         p_address_writeback, LastSendBuffer%datatag, p_comm_glb_plus, LastSendBuffer%sendreqs(3), p_err)

   END SUBROUTINE hist_writeback_var_header

   ! -----
   SUBROUTINE hist_writeback_var ( dataid, ixseg, iyseg, &
         wdata1d, wdata2d, wdata3d, wdata4d)
      
   IMPLICIT NONE

   integer,  intent(in) :: dataid, ixseg, iyseg

   real(r8), intent(in), optional :: wdata1d(:)
   real(r8), intent(in), optional :: wdata2d(:,:)
   real(r8), intent(in), optional :: wdata3d(:,:,:)
   real(r8), intent(in), optional :: wdata4d(:,:,:,:)

   ! Local Variables
   integer :: totalsize, ndim1, ndim2
   type(HistSendBufferType), pointer :: TempSendBuffer

      ! append sending buffer
      IF (.not. associated(HistSendBuffer)) THEN
         allocate (HistSendBuffer)
         LastSendBuffer => HistSendBuffer
         TotalMemSize = 0
      ELSE
         allocate (LastSendBuffer%next)
         LastSendBuffer => LastSendBuffer%next
      ENDIF
      
      LastSendBuffer%next   => null()
      
      ! clean sending buffer and free memory
      DO WHILE ((TotalMemSize > MaxHistMemSize) .and. associated(HistSendBuffer%next))

         CALL MPI_Waitall (2, HistSendBuffer%sendreqs(1:2), p_stat, p_err)
         
         TotalMemSize = TotalMemSize - size(HistSendBuffer%senddata)

         TempSendBuffer => HistSendBuffer
         HistSendBuffer => HistSendBuffer%next
         deallocate(TempSendBuffer%senddata)
         deallocate(TempSendBuffer)

      ENDDO

      LastSendBuffer%datatag = dataid*10+1

      ndim1 = 0
      ndim2 = 0

      IF (present(wdata1d)) THEN
         
         totalsize = size(wdata1d)
         allocate(LastSendBuffer%senddata(totalsize))
         LastSendBuffer%senddata = wdata1d

      ELSEIF (present(wdata2d)) THEN
         
         totalsize = size(wdata2d)
         allocate(LastSendBuffer%senddata(totalsize))
         LastSendBuffer%senddata = reshape(wdata2d, (/totalsize/))

      ELSEIF (present(wdata3d)) THEN

         ndim1 = size(wdata3d,1)
         totalsize = size(wdata3d)
         allocate(LastSendBuffer%senddata(totalsize))
         LastSendBuffer%senddata = reshape(wdata3d, (/totalsize/))
         
      ELSEIF (present(wdata4d)) THEN

         ndim1 = size(wdata4d,1)
         ndim2 = size(wdata4d,2)
         totalsize = size(wdata4d)
         allocate(LastSendBuffer%senddata(totalsize))
         LastSendBuffer%senddata = reshape(wdata4d, (/totalsize/))
      
      ENDIF
         
      TotalMemSize = TotalMemSize + totalsize

      LastSendBuffer%sendint4(1:5) = (/p_iam_glb_plus, ixseg, iyseg, ndim1, ndim2/)

      CALL mpi_isend (LastSendBuffer%sendint4(1:5), 5, MPI_INTEGER, &
         p_address_writeback, LastSendBuffer%datatag, p_comm_glb_plus, LastSendBuffer%sendreqs(1), p_err)

      CALL mpi_isend (LastSendBuffer%senddata, totalsize, MPI_REAL8, &
         p_address_writeback, LastSendBuffer%datatag, p_comm_glb_plus, LastSendBuffer%sendreqs(2), p_err)

   END SUBROUTINE hist_writeback_var

   ! -----
   SUBROUTINE hist_writeback_exit ()

   IMPLICIT NONE
   
   ! Local Variables
   integer :: dataid, nreq
   type(timenodetype),       pointer :: tempnode
   type(HistSendBufferType), pointer :: TempSendBuffer

      lasttime => null()
      DO WHILE (associated(timenodes))

         CALL MPI_WaitAll (3, timenodes%req, p_stat, p_err)
         
         tempnode  => timenodes
         timenodes => timenodes%next
         deallocate(tempnode)
      ENDDO

      LastSendBuffer => null()
      DO WHILE (associated(HistSendBuffer))
            
         IF (allocated(HistSendBuffer%senddata)) THEN
            CALL MPI_Waitall (2, HistSendBuffer%sendreqs(1:2), p_stat, p_err)
            deallocate(HistSendBuffer%senddata)
         ELSE
            CALL MPI_Waitall (3, HistSendBuffer%sendreqs(1:3), p_stat, p_err)
         ENDIF
         
         TempSendBuffer => HistSendBuffer
         HistSendBuffer => HistSendBuffer%next
         deallocate (TempSendBuffer)
         
      ENDDO

      IF (allocated(xGridDsp)) deallocate(xGridDsp)
      IF (allocated(yGridDsp)) deallocate(yGridDsp)
      IF (allocated(xGridCnt)) deallocate(xGridCnt)
      IF (allocated(yGridCnt)) deallocate(yGridCnt)
      IF (allocated(lat_c   )) deallocate(lat_c   )
      IF (allocated(lat_s   )) deallocate(lat_s   )
      IF (allocated(lat_n   )) deallocate(lat_n   )
      IF (allocated(lon_c   )) deallocate(lon_c   )
      IF (allocated(lon_w   )) deallocate(lon_w   )
      IF (allocated(lon_e   )) deallocate(lon_e   )
      
      
      IF (.not. p_is_writeback) THEN
         CALL mpi_barrier (p_comm_glb, p_err)
      ENDIF

      IF (p_is_master) THEN
         dataid = -1
         CALL mpi_send (dataid, 1, MPI_INTEGER, p_address_writeback, tag_next, p_comm_glb_plus, p_err)
      ENDIF
      
      CALL mpi_barrier (p_comm_glb_plus, p_err)

   END SUBROUTINE hist_writeback_exit

   ! ----
   character(len=256) FUNCTION basename (fullname)
      
   IMPLICIT NONE
   character(len=*), intent(in) :: fullname

   ! Local variables
   integer :: i, n

      i = len_trim (fullname) 
      DO WHILE (i > 0)
         IF (fullname(i:i) == '/') EXIT
         i = i - 1
      ENDDO

      IF (i > 0) THEN
         basename = fullname(i+1:)
      ELSE
         basename = fullname
      ENDIF

   END FUNCTION basename

END MODULE MOD_HistWriteBack
#endif
