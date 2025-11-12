#include <define.h>

#ifdef USEMPI
MODULE MOD_HistWriteBack
!----------------------------------------------------------------------------
! !DESCRIPTION:
!
!     Write out data to history files by a dedicated process.
!
!  Author: Shupeng Zhang, 11/2023
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
      character(len=256) :: filelast
      character(len=256) :: timename
      integer :: time(3)
      integer :: req (4)
      type(timenodetype), pointer :: next
   END type timenodetype

   ! time nodes
   integer :: dataid_zero = 0
   integer :: req_zero
   type(timenodetype), pointer :: timenodes, lasttime

   ! dimension information
   logical :: SDimInited = .false.

   ! 1: grid-based
   integer :: nGridData, nxGridSeg, nyGridSeg
   integer, allocatable :: xGridDsp(:), xGridCnt(:)
   integer, allocatable :: yGridDsp(:), yGridCnt(:)

   integer :: nlat, nlon
   real(r8), allocatable :: lat_c(:), lat_s(:), lat_n(:)
   real(r8), allocatable :: lon_c(:), lon_w(:), lon_e(:)

   ! Memory limits
   integer*8, parameter :: MaxHistMemSize = 1073741824_8 ! 1024^3
   integer*8 :: TotalMemSize = 0

   integer :: itime_in_file_wb

   ! tags
   integer, parameter :: tag_next = 1
   integer, parameter :: tag_time = 2
   integer, parameter :: tag_dims = 3

CONTAINS

   SUBROUTINE hist_writeback_daemon ()

   USE MOD_Namelist, only: DEF_HIST_FREQ
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   ! Local Variables
   integer :: dataid, tag
   integer :: time(3), ndims, ndim1, ndim2, dimlens(4), compress
   integer :: i, idata, isrc, ixseg, iyseg, xdsp, ydsp, xcnt, ycnt, idim1, idim2

   integer               :: recvint4 (5)
   character(len=256)    :: recvchar (9)
   real(r8), allocatable :: datathis (:)

   character(len=256)    :: filename, filelast, dataname, longname, units
   character(len=256)    :: dim1name, dim2name, dim3name, dim4name, dim5name
   logical               :: fexists

   real(r8), allocatable :: wdata1d(:),  wdata2d(:,:), wdata3d(:,:,:), wdata4d(:,:,:,:)
   real(r8), allocatable :: tmp3d(:,:,:), tmp4d(:,:,:,:)


      DO WHILE (.true.)

         CALL mpi_recv (dataid, 1, MPI_INTEGER, &
            MPI_ANY_SOURCE, tag_next, p_comm_glb_plus, p_stat, p_err)

         IF (dataid < 0) THEN

            EXIT

         ELSEIF (dataid == 0) THEN

            CALL mpi_recv (filename, 256, MPI_CHARACTER, &
               MPI_ANY_SOURCE, tag_time, p_comm_glb_plus, p_stat, p_err)

            CALL mpi_recv (filelast, 256, MPI_CHARACTER, &
               MPI_ANY_SOURCE, tag_time, p_comm_glb_plus, p_stat, p_err)

            CALL mpi_recv (dataname, 256, MPI_CHARACTER, &
               MPI_ANY_SOURCE, tag_time, p_comm_glb_plus, p_stat, p_err)

            CALL mpi_recv (time, 3, MPI_INTEGER, &
               MPI_ANY_SOURCE, tag_time, p_comm_glb_plus, p_stat, p_err)

            IF (.not. SDimInited) THEN

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

               SDimInited = .true.

            ENDIF

            inquire (file=filename, exist=fexists)
            IF ((.not. fexists) .or. (trim(filename) /= trim(filelast))) THEN

               CALL ncio_create_file (trim(filename))

               CALL ncio_define_dimension(filename, 'time', 0)

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

               CALL ncio_write_colm_dimension (filename)

            ENDIF

            CALL ncio_write_time (filename, dataname, time, itime_in_file_wb, DEF_HIST_FREQ)

         ELSE

            !--------------------------------
            ! receive and write history data.
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
               CASE (2:3)

                  dimlens = (/nlon, nlat, 0, 0/)

                  IF (.not. allocated(wdata2d)) THEN
                     allocate (wdata2d (nlon,nlat))
                  ENDIF

                  allocate (datathis(xcnt*ycnt))
                  CALL mpi_recv (datathis, xcnt*ycnt, MPI_REAL8, &
                     isrc, tag, p_comm_glb_plus, p_stat, p_err)

                  wdata2d(xdsp+1:xdsp+xcnt, ydsp+1:ydsp+ycnt) = &
                     reshape(datathis,(/xcnt,ycnt/))

               CASE (4)

                  dimlens = (/ndim1, nlon, nlat, 0/)

                  IF (.not. allocated(wdata3d)) THEN
                     allocate (wdata3d (nlon,nlat,ndim1))
                     allocate (tmp3d   (ndim1,nlon,nlat))
                  ENDIF

                  allocate (datathis(ndim1*xcnt*ycnt))
                  CALL mpi_recv (datathis, ndim1*xcnt*ycnt, MPI_REAL8, &
                     isrc, tag, p_comm_glb_plus, p_stat, p_err)

                  tmp3d = reshape(datathis,(/ndim1,xcnt,ycnt/))
                  DO idim1 = 1, ndim1
                     wdata3d(xdsp+1:xdsp+xcnt, ydsp+1:ydsp+ycnt, idim1) = tmp3d(idim1, :, :)
                  ENDDO

               CASE (5)

                  dimlens = (/ndim1, ndim2, nlon, nlat/)

                  IF (.not. allocated(wdata4d)) THEN
                     allocate (wdata4d (nlon,nlat,ndim1,ndim2))
                     allocate (tmp4d   (ndim1,ndim2,nlon,nlat))
                  ENDIF

                  allocate (datathis(ndim1*ndim2*xcnt*ycnt))
                  CALL mpi_recv (datathis, ndim1*ndim2*xcnt*ycnt, MPI_REAL8, &
                     isrc, tag, p_comm_glb_plus, p_stat, p_err)

                  tmp4d = reshape(datathis,(/ndim1, ndim2, xcnt, ycnt/))

                  DO idim1 = 1, ndim1
                     DO idim2 = 1, ndim2
                        wdata4d(xdsp+1:xdsp+xcnt, ydsp+1:ydsp+ycnt, idim1, idim2) = &
                           tmp4d(idim1, idim2, :, :)
                     ENDDO
                  ENDDO

               ENDSELECT

               deallocate (datathis)

            ENDDO


            IF (ndims >= 4) CALL ncio_define_dimension (filename, dim3name, dimlens(1))
            IF (ndims >= 5) CALL ncio_define_dimension (filename, dim4name, dimlens(2))

            SELECTCASE (ndims)
            CASE (2) ! for variables with [lon,lat]

               CALL ncio_write_serial (filename, dataname, wdata2d, dim1name, dim2name, compress)

               deallocate(wdata2d)
            CASE (3) ! for variables with [lon,lat,time]

               CALL ncio_write_serial_time (filename, dataname, itime_in_file_wb, wdata2d, &
                  dim1name, dim2name, dim3name, compress)

               deallocate(wdata2d)
            CASE (4) ! for variables with [lon,lat,dim3,time]

               CALL ncio_write_serial_time (filename, dataname, itime_in_file_wb, wdata3d, &
                  dim1name, dim2name, dim3name, dim4name, compress)

               deallocate(tmp3d  )
               deallocate(wdata3d)
            CASE (5) ! for variables with [lon,lat,dim3,dim4,time]

               CALL ncio_write_serial_time (filename, dataname, itime_in_file_wb, wdata4d, &
                  dim1name, dim2name, dim3name, dim4name, dim5name, compress)

               deallocate(tmp4d  )
               deallocate(wdata4d)
            ENDSELECT

            IF (itime_in_file_wb <= 1) THEN
               CALL ncio_put_attr (filename, dataname, 'long_name', longname)
               CALL ncio_put_attr (filename, dataname, 'units', units)
               CALL ncio_put_attr (filename, dataname, 'missing_value', spval)
            ENDIF

            write(*,'(3A,I0,2A)') 'HIST WriteBack: ', trim(basename(filename)), &
               ' (time ', itime_in_file_wb, '): ', trim(dataname)

         ENDIF

      ENDDO

   END SUBROUTINE hist_writeback_daemon

   ! -----
   SUBROUTINE hist_writeback_latlon_time (filename, filelast, timename, time, HistConcat)

   USE MOD_Namelist
   USE MOD_Grid
   IMPLICIT NONE

   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: filelast
   character (len=*), intent(in) :: timename
   integer, intent(in)  :: time(3)
   type(grid_concat_type), intent(in) :: HistConcat

   ! Local Variables
   integer :: i
   logical :: senddone
   integer :: sendstat(MPI_STATUS_SIZE,4)
   type(timenodetype), pointer :: tempnode

      IF (.not. associated(timenodes)) THEN
         allocate (timenodes)
         lasttime => timenodes
      ELSE
         allocate (lasttime%next)
         lasttime => lasttime%next
      ENDIF

      lasttime%filename = filename
      lasttime%filelast = filelast
      lasttime%timename = timename
      lasttime%time = time
      lasttime%next => null()

      CALL mpi_isend (dataid_zero, 1, MPI_INTEGER, &
         p_address_writeback, tag_next, p_comm_glb_plus, req_zero, p_err)

      CALL mpi_isend (lasttime%filename, 256, MPI_CHARACTER, &
         p_address_writeback, tag_time, p_comm_glb_plus, lasttime%req(1), p_err)

      CALL mpi_isend (lasttime%filelast, 256, MPI_CHARACTER, &
         p_address_writeback, tag_time, p_comm_glb_plus, lasttime%req(2), p_err)

      CALL mpi_isend (lasttime%timename, 256, MPI_CHARACTER, &
         p_address_writeback, tag_time, p_comm_glb_plus, lasttime%req(3), p_err)

      CALL mpi_isend (lasttime%time, 3, MPI_INTEGER, &
         p_address_writeback, tag_time, p_comm_glb_plus, lasttime%req(4), p_err)


      IF (.not. SDimInited) THEN

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

         CALL mpi_send (nGridData, 1, MPI_INTEGER, p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (nxGridSeg, 1, MPI_INTEGER, p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (nyGridSeg, 1, MPI_INTEGER, p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)

         CALL mpi_send (xGridDsp, nxGridSeg, MPI_INTEGER, p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (xGridCnt, nxGridSeg, MPI_INTEGER, p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (yGridDsp, nyGridSeg, MPI_INTEGER, p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (yGridCnt, nyGridSeg, MPI_INTEGER, p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)

         nlat = HistConcat%ginfo%nlat
         nlon = HistConcat%ginfo%nlon
         allocate(lat_c(nlat));  lat_c = HistConcat%ginfo%lat_c
         allocate(lat_s(nlat));  lat_s = HistConcat%ginfo%lat_s
         allocate(lat_n(nlat));  lat_n = HistConcat%ginfo%lat_n
         allocate(lon_c(nlon));  lon_c = HistConcat%ginfo%lon_c
         allocate(lon_w(nlon));  lon_w = HistConcat%ginfo%lon_w
         allocate(lon_e(nlon));  lon_e = HistConcat%ginfo%lon_e

         CALL mpi_send (nlat,     1, MPI_INTEGER, p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (lat_c, nlat, MPI_REAL8,   p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (lat_s, nlat, MPI_REAL8,   p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (lat_n, nlat, MPI_REAL8,   p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (nlon,     1, MPI_INTEGER, p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (lon_c, nlon, MPI_REAL8,   p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (lon_w, nlon, MPI_REAL8,   p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)
         CALL mpi_send (lon_e, nlon, MPI_REAL8,   p_address_writeback, tag_dims, &
                        p_comm_glb_plus, p_err)

         SDimInited = .true.

      ENDIF

      DO WHILE (associated(timenodes%next))

         CALL MPI_TestAll (4, timenodes%req, senddone, sendstat(:,1:4), p_err)

         IF (senddone) THEN
            tempnode  => timenodes
            timenodes => timenodes%next
            deallocate(tempnode)
         ELSE
            EXIT
         ENDIF
      ENDDO

   END SUBROUTINE hist_writeback_latlon_time

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

         CALL MPI_Testall (3, HistSendBuffer%sendreqs, senddone, sendstat(:,1:3), p_err)

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
         p_address_writeback, LastSendBuffer%datatag, &
         p_comm_glb_plus, LastSendBuffer%sendreqs(2), p_err)

      CALL mpi_isend (LastSendBuffer%sendchar, 256*9, MPI_CHARACTER, &
         p_address_writeback, LastSendBuffer%datatag, &
         p_comm_glb_plus, LastSendBuffer%sendreqs(3), p_err)

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
   logical :: senddone
   integer :: sendstat(MPI_STATUS_SIZE,2)
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

      LastSendBuffer%next => null()

      ! clean sending buffer and free memory
      DO WHILE (associated(HistSendBuffer%next))

         IF (TotalMemSize > MaxHistMemSize) THEN
            CALL MPI_Waitall (2, HistSendBuffer%sendreqs(1:2), sendstat(:,1:2), p_err)
         ELSE
            CALL MPI_Testall (2, HistSendBuffer%sendreqs(1:2), senddone, sendstat(:,1:2), p_err)
            IF (.not. senddone) EXIT
         ENDIF

         TotalMemSize = TotalMemSize - size(HistSendBuffer%senddata)

         TempSendBuffer => HistSendBuffer
         HistSendBuffer => HistSendBuffer%next
         deallocate(TempSendBuffer%senddata)
         deallocate(TempSendBuffer)

      ENDDO

      LastSendBuffer%datatag = dataid*10+1

      ndim1 = 0
      ndim2 = 0

      IF (present(wdata2d)) THEN

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
         p_address_writeback, LastSendBuffer%datatag, &
         p_comm_glb_plus, LastSendBuffer%sendreqs(1), p_err)

      CALL mpi_isend (LastSendBuffer%senddata, totalsize, MPI_REAL8, &
         p_address_writeback, LastSendBuffer%datatag, &
         p_comm_glb_plus, LastSendBuffer%sendreqs(2), p_err)

   END SUBROUTINE hist_writeback_var

   ! -----
   SUBROUTINE hist_writeback_exit ()

   IMPLICIT NONE

   ! Local Variables
   integer :: dataid, nreq
   integer :: sendstat(MPI_STATUS_SIZE,4)
   type(timenodetype),       pointer :: tempnode
   type(HistSendBufferType), pointer :: TempSendBuffer

      lasttime => null()
      DO WHILE (associated(timenodes))

         CALL MPI_WaitAll (4, timenodes%req, sendstat(:,1:4), p_err)

         tempnode  => timenodes
         timenodes => timenodes%next
         deallocate(tempnode)
      ENDDO

      LastSendBuffer => null()
      DO WHILE (associated(HistSendBuffer))

         IF (allocated(HistSendBuffer%senddata)) THEN
            CALL MPI_Waitall (2, HistSendBuffer%sendreqs(1:2), sendstat(:,1:2), p_err)
            deallocate(HistSendBuffer%senddata)
         ELSE
            CALL MPI_Waitall (3, HistSendBuffer%sendreqs(1:3), sendstat(:,1:3), p_err)
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
         CALL mpi_send (dataid, 1, MPI_INTEGER, p_address_writeback, &
                        tag_next, p_comm_glb_plus, p_err)
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
