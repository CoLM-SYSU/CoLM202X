#include <define.h>

#ifdef USEMPI
module MOD_HistWriteBack
   !----------------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !     Write out data to history files by a dedicated process.
   !
   ! Author: Shupeng Zhang, 11/2023
   !----------------------------------------------------------------------------

   use MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial

   ! type of send buffer
   type :: HistSendBufferType
      integer :: dataid
      integer :: nreq
      integer,  allocatable :: sendreqs (:)
      integer,  allocatable :: sendstat (:,:)
      integer               :: sendint4 (8)
      character(len=256)    :: sendchar (7)
      real(r8), allocatable :: senddata (:)
      type(HistSendBufferType), pointer :: next
   END type HistSendBufferType

   ! Sending Variables
   type(HistSendBufferType), pointer :: HistSendBuffer
   type(HistSendBufferType), pointer :: ThisSendBuffer

   integer*8, parameter :: MaxHistMemSize  = 8589934592_8 ! 8*1024^3 
   integer*8, parameter :: MaxHistMesgSize = 8388608_8    ! 8*1024^2 
   
   integer*8 :: TotalMemSize = 0

contains

   ! -----
   subroutine hist_writeback_daemon ()

      IMPLICIT NONE

      ! Local Variables
      integer :: dataid, tag
      integer :: totalsize, ndata, ndims, dimlens(4), itime, compress
      integer :: i, idata, datadisp, datasize
      integer               :: recvint4(8)
      character(len=256)    :: recvchar(7), filename, dataname, dimnames(5)
      real(r8), allocatable :: recvdata(:), datathis(:)
      real(r8), allocatable :: wdata1d(:), wdata2d(:,:), wdata3d(:,:,:), wdata4d(:,:,:,:)


      DO WHILE (.true.)

         CALL mpi_recv (dataid, 1, MPI_INTEGER, &
            MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb_plus, p_stat, p_err)

         IF (dataid < 0) THEN
            EXIT
         ENDIF

         tag = dataid * 10000
         CALL mpi_recv (recvint4(1:8), 8, MPI_INTEGER, &
            MPI_ANY_SOURCE, tag, p_comm_glb_plus, p_stat, p_err)

         ndims    = recvint4(1)
         dimlens  = recvint4(2:5)
         itime    = recvint4(6)
         compress = recvint4(7)
         ndata    = recvint4(8)

         tag = tag + 1
         CALL mpi_recv (recvchar, 256*7, MPI_CHARACTER, &
            MPI_ANY_SOURCE, tag, p_comm_glb_plus, p_stat, p_err)
      
         filename = recvchar(1)
         dataname = recvchar(2) 
         dimnames = recvchar(3:7)
         
         totalsize = product(dimlens(1:ndims))

         allocate (recvdata (totalsize))

         datadisp = 0
         DO idata = 1, ndata

            tag = tag + 1
            CALL mpi_recv (datasize, 1, MPI_INTEGER, &
               MPI_ANY_SOURCE, tag, p_comm_glb_plus, p_stat, p_err)

            tag = tag + 1
            allocate (datathis (datasize))
            CALL mpi_recv (datathis, datasize, MPI_REAL8, &
               MPI_ANY_SOURCE, tag, p_comm_glb_plus, p_stat, p_err)

            recvdata(datadisp+1:datadisp+datasize) = datathis

            datadisp = datadisp + datasize
            deallocate(datathis)
         ENDDO 

         select case (ndims)
         case (1)
            allocate(wdata1d (dimlens(1)))
            wdata1d = recvdata

            call ncio_write_serial_time (filename, dataname, itime, wdata1d, &
               dimnames(1), dimnames(2), compress)

            deallocate(wdata1d)
         case (2)
            allocate(wdata2d (dimlens(1),dimlens(2)))
            wdata2d = reshape(recvdata, dimlens(1:2))

            call ncio_write_serial_time (filename, dataname, itime, wdata2d, &
               dimnames(1), dimnames(2), dimnames(3), compress)

            deallocate(wdata2d)
         case (3)
            allocate(wdata3d (dimlens(1),dimlens(2),dimlens(3)))
            wdata3d = reshape(recvdata, dimlens(1:3))

            call ncio_write_serial_time (filename, dataname, itime, wdata3d, &
               dimnames(1), dimnames(2), dimnames(3), dimnames(4), compress)

            deallocate(wdata3d)
         case (4)
            allocate(wdata4d (dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
            wdata4d = reshape(recvdata, dimlens(1:4))

            call ncio_write_serial_time (filename, dataname, itime, wdata4d, &
               dimnames(1), dimnames(2), dimnames(3), dimnames(4), dimnames(5), compress)
            
            deallocate(wdata4d)
         END select

         deallocate (recvdata)
      
      ENDDO

   END subroutine hist_writeback_daemon

   ! -----
   SUBROUTINE hist_writeback ( dataid,        &
         filename, dataname, itime,   dimnames, compress, &
         wdata1d,  wdata2d,  wdata3d, wdata4d)
      
      IMPLICIT NONE

      integer, intent(in)          :: dataid, itime, compress
      character(len=*), intent(in) :: filename, dataname, dimnames(5)
      real(r8), intent(in), optional :: wdata1d(:)
      real(r8), intent(in), optional :: wdata2d(:,:)
      real(r8), intent(in), optional :: wdata3d(:,:,:)
      real(r8), intent(in), optional :: wdata4d(:,:,:,:)

      ! Local Variables
      logical :: found
      integer :: totalsize, ndims, dimlens(4), ndata, idata, datadisp, datasize
      integer :: nreq, ireq, tag
      integer, pointer :: preqs(:)

      IF (.not. associated(HistSendBuffer)) THEN
         allocate (HistSendBuffer)
         HistSendBuffer%dataid =  -1
         HistSendBuffer%next   => null()
         TotalMemSize = 0
      ENDIF
      
      found = .false.

      ThisSendBuffer => HistSendBuffer
      DO WHILE (.true.)
         IF (ThisSendBuffer%dataid < 0) THEN
            found = .true.
         ELSE
            nreq = ThisSendBuffer%nreq
            CALL MPI_Testall (nreq, ThisSendBuffer%sendreqs(1:nreq), found, &
               ThisSendBuffer%sendstat(:,1:nreq), p_err)
         ENDIF

         IF ((.not. found) .and. associated(ThisSendBuffer%next)) THEN
            ThisSendBuffer => ThisSendBuffer%next
         ELSE
            EXIT
         ENDIF
      ENDDO

      IF (.not. found) THEN
         IF (TotalMemSize < MaxHistMemSize) THEN
            allocate (ThisSendBuffer%next)
            ThisSendBuffer => ThisSendBuffer%next
            ThisSendBuffer%next => null()
         ELSE
            ThisSendBuffer => HistSendBuffer
            nreq = ThisSendBuffer%nreq
            CALL MPI_Waitall (nreq, ThisSendBuffer%sendreqs(1:nreq), p_stat, p_err)
         ENDIF
      ENDIF
            
      ThisSendBuffer%dataid = dataid

      IF (allocated(ThisSendBuffer%sendreqs)) THEN
         deallocate(ThisSendBuffer%sendreqs)
      ENDIF

      IF (allocated(ThisSendBuffer%sendstat)) THEN
         deallocate(ThisSendBuffer%sendstat)
      ENDIF

      IF (allocated(ThisSendBuffer%senddata)) THEN
         TotalMemSize = TotalMemSize - size(ThisSendBuffer%senddata)
         deallocate(ThisSendBuffer%senddata)
      ENDIF

      IF (present(wdata1d)) THEN
         
         totalsize = size(wdata1d)
         ndims = 1
         dimlens(1) = totalsize
         
         allocate(ThisSendBuffer%senddata(totalsize))
         ThisSendBuffer%senddata = wdata1d

      ELSEIF (present(wdata2d)) THEN
         
         totalsize = size(wdata2d)
         ndims = 2
         dimlens(1:2) = shape(wdata2d)

         allocate(ThisSendBuffer%senddata(totalsize))
         ThisSendBuffer%senddata = reshape(wdata2d, (/totalsize/))

      ELSEIF (present(wdata3d)) THEN

         totalsize = size(wdata3d)
         ndims = 3
         dimlens(1:3) = shape(wdata3d)

         allocate(ThisSendBuffer%senddata(totalsize))
         ThisSendBuffer%senddata = reshape(wdata3d, (/totalsize/))

      ELSEIF (present(wdata4d)) THEN

         totalsize = size(wdata4d)
         ndims = 4 
         dimlens(1:4) = shape(wdata4d)

         allocate(ThisSendBuffer%senddata(totalsize))
         ThisSendBuffer%senddata = reshape(wdata4d, (/totalsize/))

      ENDIF
         
      TotalMemSize = TotalMemSize + totalsize

      ndata = totalsize/MaxHistMesgSize
      IF (mod(totalsize,MaxHistMesgSize) /= 0) THEN
         ndata = ndata + 1
      ENDIF

      ThisSendBuffer%sendint4 = (/ &
         ndims, dimlens(1:4), itime, compress, ndata /)

      ThisSendBuffer%sendchar(1)   = filename
      ThisSendBuffer%sendchar(2)   = dataname
      ThisSendBuffer%sendchar(3:7) = dimnames(1:5)
      
      ThisSendBuffer%nreq = ndata*2+3
      
      nreq = ThisSendBuffer%nreq
      allocate (ThisSendBuffer%sendreqs (nreq))
      allocate (ThisSendBuffer%sendstat (MPI_STATUS_SIZE,nreq))

      preqs => ThisSendBuffer%sendreqs

      ireq = 1
      CALL mpi_isend (dataid, 1, MPI_INTEGER, &
         0, mpi_tag_mesg, p_comm_glb_plus, preqs(ireq), p_err)

      ireq = ireq + 1
      tag  = dataid * 10000
      CALL mpi_isend (ThisSendBuffer%sendint4, 8, MPI_INTEGER, &
         0, tag, p_comm_glb_plus, preqs(ireq), p_err)

      tag  = tag  + 1
      ireq = ireq + 1
      CALL mpi_isend (ThisSendBuffer%sendchar, 256*7, MPI_CHARACTER, &
         0, tag, p_comm_glb_plus, preqs(ireq), p_err)

      datadisp = 0
      DO idata = 1, ndata
         datasize = min(totalsize-datadisp, MaxHistMesgSize)

         ireq = ireq + 1
         tag  = tag  + 1
         CALL mpi_isend (datasize, 1, MPI_INTEGER, &
            0, tag, p_comm_glb_plus, preqs(ireq), p_err)
         
         ireq = ireq + 1
         tag  = tag  + 1
         CALL mpi_isend (ThisSendBuffer%senddata(datadisp+1:datadisp+datasize), datasize, &
            MPI_REAL8, 0, tag, p_comm_glb_plus, preqs(ireq), p_err)

         datadisp = datadisp + datasize
      ENDDO

   END SUBROUTINE hist_writeback

   ! -----
   subroutine hist_writeback_exit ()

      IMPLICIT NONE
   
      integer :: nreq, dataid
      type(HistSendBufferType), pointer :: pp

      ThisSendBuffer => HistSendBuffer
      DO WHILE (associated(ThisSendBuffer))
            
         IF (ThisSendBuffer%dataid > 0) THEN
            nreq = ThisSendBuffer%nreq
            CALL MPI_Waitall (nreq, ThisSendBuffer%sendreqs(1:nreq), p_stat, p_err)
         ENDIF
         
         pp => ThisSendBuffer%next
         
         IF (allocated(ThisSendBuffer%sendreqs)) THEN
            deallocate(ThisSendBuffer%sendreqs)
         ENDIF
         IF (allocated(ThisSendBuffer%sendstat)) THEN
            deallocate(ThisSendBuffer%sendstat)
         ENDIF
         IF (allocated(ThisSendBuffer%senddata)) THEN
            deallocate(ThisSendBuffer%senddata)
         ENDIF

         deallocate(ThisSendBuffer)

         ThisSendBuffer => pp
      ENDDO
      
      IF (.not. p_is_writeback) THEN
         CALL mpi_barrier (p_comm_glb, p_err)
      ENDIF

      IF (p_is_master) THEN
         dataid = -1
         CALL mpi_send (dataid, 1, MPI_INTEGER, 0, mpi_tag_mesg, p_comm_glb_plus, p_err)
      ENDIF
      
      CALL mpi_barrier (p_comm_glb_plus, p_err)

   END subroutine hist_writeback_exit


end module MOD_HistWriteBack
#endif
