#include <define.h>

MODULE MOD_AggregationRequestData

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Aggregation Utilities.
!
!    On IO processes, a data daemon is running to provide data
!       at fine resolutions for worker processes.
!    On worker processes, request is sent to IO processes and
!       data is returned from IO processes.
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   IMPLICIT NONE

   PUBLIC :: aggregation_request_data

#ifdef USEMPI
   PUBLIC :: aggregation_data_daemon
   PUBLIC :: aggregation_worker_done
#endif

! ---- subroutines ----
CONTAINS

#ifdef USEMPI
   SUBROUTINE aggregation_data_daemon (grid_in, &
         data_r8_2d_in1, data_r8_2d_in2, data_r8_2d_in3, data_r8_2d_in4, &
         data_r8_2d_in5, data_r8_2d_in6,                                 &
         data_r8_3d_in1, n1_r8_3d_in1  , data_r8_3d_in2, n1_r8_3d_in2,   &
         data_i4_2d_in1, data_i4_2d_in2)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType

   IMPLICIT NONE

   type (grid_type), intent(in) :: grid_in

   ! 2D REAL data
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in1
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in2
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in3
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in4
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in5
   type (block_data_real8_2d), intent(in), optional :: data_r8_2d_in6

   ! 3D REAL data
   integer, intent(in), optional :: n1_r8_3d_in1
   type (block_data_real8_3d), intent(in), optional :: data_r8_3d_in1

   integer, intent(in), optional :: n1_r8_3d_in2
   type (block_data_real8_3d), intent(in), optional :: data_r8_3d_in2

   ! 2D INTEGER data
   type (block_data_int32_2d), intent(in), optional :: data_i4_2d_in1
   type (block_data_int32_2d), intent(in), optional :: data_i4_2d_in2

   ! Local Variables
   integer :: nreq, ireq, rmesg(2), isrc, idest
   integer :: xblk, yblk, xloc, yloc
   integer,  allocatable :: ylist(:), xlist(:)

   real(r8), allocatable :: sbuf_r8_1d(:), sbuf_r8_2d(:,:)
   integer , allocatable :: sbuf_i4_1d(:)

   logical,  allocatable :: worker_done (:)

      IF (p_is_io) THEN

         allocate (worker_done (0:p_np_worker-1))

         worker_done(:) = .false.
         DO WHILE (any(.not. worker_done))

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc = rmesg(1)
            nreq = rmesg(2)

            IF (nreq > 0) THEN

               allocate (xlist (nreq))
               allocate (ylist (nreq))

               CALL mpi_recv (xlist, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (ylist, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               idest = isrc

               allocate (sbuf_r8_1d (nreq))

               IF (present(data_r8_2d_in1)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_r8_2d_in2)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_r8_2d_in3)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in3%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_r8_2d_in4)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in4%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_r8_2d_in5)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in5%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_r8_2d_in6)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_1d(ireq) = data_r8_2d_in6%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_1d, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF


               deallocate (sbuf_r8_1d)

               IF (present(data_r8_3d_in1) .and. present(n1_r8_3d_in1)) THEN

                  allocate (sbuf_r8_2d (n1_r8_3d_in1,nreq))
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_2d(:,ireq) = data_r8_3d_in1%blk(xblk,yblk)%val(:,xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_2d, n1_r8_3d_in1*nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)

                  deallocate (sbuf_r8_2d)
               ENDIF

               IF (present(data_r8_3d_in2) .and. present(n1_r8_3d_in2)) THEN

                  allocate (sbuf_r8_2d (n1_r8_3d_in2,nreq))
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_r8_2d(:,ireq) = data_r8_3d_in2%blk(xblk,yblk)%val(:,xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_r8_2d, n1_r8_3d_in2*nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)

                  deallocate (sbuf_r8_2d)
               ENDIF

               allocate (sbuf_i4_1d (nreq))

               IF (present(data_i4_2d_in1)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_i4_1d(ireq) = data_i4_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_i4_1d, nreq, MPI_INTEGER, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_i4_2d_in2)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuf_i4_1d(ireq) = data_i4_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  CALL mpi_send (sbuf_i4_1d, nreq, MPI_INTEGER, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               deallocate (sbuf_i4_1d)

               deallocate (ylist)
               deallocate (xlist)

            ELSE
               worker_done(p_itis_worker(isrc)) = .true.
            ENDIF

         ENDDO

         deallocate (worker_done)

      ENDIF

   END SUBROUTINE aggregation_data_daemon

#endif

   !----------------------------------------------------
   SUBROUTINE aggregation_request_data (     &
         pixelset, iset, grid_in, zip, area, &
         data_r8_2d_in1, data_r8_2d_out1, &
         data_r8_2d_in2, data_r8_2d_out2, &
         data_r8_2d_in3, data_r8_2d_out3, &
         data_r8_2d_in4, data_r8_2d_out4, &
         data_r8_2d_in5, data_r8_2d_out5, &
         data_r8_2d_in6, data_r8_2d_out6, &
         data_r8_3d_in1, data_r8_3d_out1, n1_r8_3d_in1, lb1_r8_3d_in1, &
         data_r8_3d_in2, data_r8_3d_out2, n1_r8_3d_in2, lb1_r8_3d_in2, &
         data_i4_2d_in1, data_i4_2d_out1, &
         data_i4_2d_in2, data_i4_2d_out2, &
         filledvalue_i4)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_Pixelset
   USE MOD_Utils

   IMPLICIT NONE

   type (pixelset_type), intent(in) :: pixelset
   integer, intent(in)  :: iset

   type (grid_type), intent(in) :: grid_in
   logical, intent(in) :: zip

   real(r8), allocatable, intent(out), optional :: area(:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in1
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out1 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in2
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out2 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in3
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out3 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in4
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out4 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in5
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out5 (:)

   type (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in6
   real(r8), allocatable,      intent(out), optional :: data_r8_2d_out6 (:)

   integer, intent(in), optional :: n1_r8_3d_in1, lb1_r8_3d_in1
   type (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in1
   real(r8), allocatable,      intent(out), optional :: data_r8_3d_out1 (:,:)

   integer, intent(in), optional :: n1_r8_3d_in2, lb1_r8_3d_in2
   type (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in2
   real(r8), allocatable,      intent(out), optional :: data_r8_3d_out2 (:,:)

   type (block_data_int32_2d), intent(in),  optional :: data_i4_2d_in1
   integer, allocatable,       intent(out), optional :: data_i4_2d_out1 (:)

   type (block_data_int32_2d), intent(in),  optional :: data_i4_2d_in2
   integer, allocatable,       intent(out), optional :: data_i4_2d_out2 (:)

   integer, intent(in), optional :: filledvalue_i4

   ! Local Variables
   integer :: totalreq, ireq, nreq, smesg(2), isrc, idest, iproc
   integer :: ilon, ilat, xblk, yblk, xloc, yloc, iloc, nx, ny, ix, iy, ig
   integer :: ie, ipxstt, ipxend, npxl, ipxl, lb1, xgrdthis, ygrdthis
   integer,  allocatable :: ylist(:), xlist(:), ipt(:), ibuf(:), rbuf_i4_1d(:)
   integer,  allocatable :: xsorted(:), ysorted(:), xy2d(:,:)
   real(r8), allocatable :: area2d(:,:), rbuf_r8_1d(:), rbuf_r8_2d(:,:)
   logical,  allocatable :: msk(:)


      ie     = pixelset%ielm  (iset)
      ipxstt = pixelset%ipxstt(iset)
      ipxend = pixelset%ipxend(iset)
      npxl   = ipxend - ipxstt + 1

      IF (zip) THEN

         allocate (xsorted(npxl))
         allocate (ysorted(npxl))

         nx = 0; ny = 0
         DO ipxl = ipxstt, ipxend
            xgrdthis = grid_in%xgrd(mesh(ie)%ilon(ipxl))
            ygrdthis = grid_in%ygrd(mesh(ie)%ilat(ipxl))
            CALL insert_into_sorted_list1 (xgrdthis, nx, xsorted, iloc)
            CALL insert_into_sorted_list1 (ygrdthis, ny, ysorted, iloc)
         ENDDO

         allocate (xy2d (nx,ny));     xy2d(:,:) = 0

         IF (present(area)) THEN
            allocate(area2d(nx,ny));  area2d(:,:) = 0.
         ENDIF

         DO ipxl = ipxstt, ipxend
            xgrdthis = grid_in%xgrd(mesh(ie)%ilon(ipxl))
            ygrdthis = grid_in%ygrd(mesh(ie)%ilat(ipxl))

            ix = find_in_sorted_list1(xgrdthis, nx, xsorted)
            iy = find_in_sorted_list1(ygrdthis, ny, ysorted)

            xy2d(ix,iy) = xy2d(ix,iy) + 1

            IF (present(area)) THEN
               area2d(ix,iy) = area2d(ix,iy) + areaquad (&
                  pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                  pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
            ENDIF
         ENDDO

         totalreq = count(xy2d > 0)

         allocate (xlist (totalreq))
         allocate (ylist (totalreq))

         IF (present(area)) allocate(area(totalreq))

         ig = 0
         DO ix = 1, nx
            DO iy = 1, ny
               IF (xy2d(ix,iy) > 0) THEN
                  ig = ig + 1
                  xlist(ig) = xsorted(ix)
                  ylist(ig) = ysorted(iy)
                  IF (present(area)) area (ig) = area2d(ix,iy)
               ENDIF
            ENDDO
         ENDDO

         deallocate (xsorted, ysorted, xy2d)
         IF (present(area)) deallocate (area2d)

      ELSE

         allocate(xlist (npxl))
         allocate(ylist (npxl))

         IF (present(area)) allocate (area (npxl))

         totalreq = npxl
         DO ipxl = ipxstt, ipxend
            xlist(ipxl-ipxstt+1) = grid_in%xgrd(mesh(ie)%ilon(ipxl))
            ylist(ipxl-ipxstt+1) = grid_in%ygrd(mesh(ie)%ilat(ipxl))
            IF (present(area)) THEN
               area(ipxl-ipxstt+1) = areaquad (&
                  pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                  pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
            ENDIF
         ENDDO

      ENDIF

      IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1))  allocate (data_r8_2d_out1 (totalreq))
      IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2))  allocate (data_r8_2d_out2 (totalreq))
      IF (present(data_r8_2d_in3) .and. present(data_r8_2d_out3))  allocate (data_r8_2d_out3 (totalreq))
      IF (present(data_r8_2d_in4) .and. present(data_r8_2d_out4))  allocate (data_r8_2d_out4 (totalreq))
      IF (present(data_r8_2d_in5) .and. present(data_r8_2d_out5))  allocate (data_r8_2d_out5 (totalreq))
      IF (present(data_r8_2d_in6) .and. present(data_r8_2d_out6))  allocate (data_r8_2d_out6 (totalreq))

      IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
         IF (present(lb1_r8_3d_in1)) THEN
            lb1 = lb1_r8_3d_in1
         ELSE
            lb1 = 1
         ENDIF
         allocate (data_r8_3d_out1 (lb1:lb1-1+n1_r8_3d_in1,totalreq))
      ENDIF

      IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
         IF (present(lb1_r8_3d_in2)) THEN
            lb1 = lb1_r8_3d_in2
         ELSE
            lb1 = 1
         ENDIF
         allocate (data_r8_3d_out2 (lb1:lb1-1+n1_r8_3d_in2,totalreq))
      ENDIF

      IF (present(data_i4_2d_in1) .and. present(data_i4_2d_out1)) THEN
         allocate (data_i4_2d_out1 (totalreq))
         IF (present(filledvalue_i4)) THEN
            data_i4_2d_out1 = filledvalue_i4
         ENDIF
      ENDIF

      IF (present(data_i4_2d_in2) .and. present(data_i4_2d_out2)) THEN
         allocate (data_i4_2d_out2 (totalreq))
         IF (present(filledvalue_i4)) THEN
            data_i4_2d_out2 = filledvalue_i4
         ENDIF
      ENDIF

#ifdef USEMPI

      allocate (ipt (totalreq))
      allocate (msk (totalreq))

      ipt(:) = -1

      DO ireq = 1, totalreq
         xblk = grid_in%xblk(xlist(ireq))
         yblk = grid_in%yblk(ylist(ireq))
         ipt(ireq) = gblock%pio(xblk,yblk)
      ENDDO

      DO iproc = 0, p_np_io-1
         msk = (ipt == p_address_io(iproc))
         nreq = count(msk)

         IF (nreq > 0) THEN

            smesg = (/p_iam_glb, nreq/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

            allocate (ibuf (nreq))

            ibuf = pack(xlist(1:totalreq), msk)
            CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            ibuf = pack(ylist(1:totalreq), msk)
            CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            isrc = idest

            allocate (rbuf_r8_1d (nreq))

            IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out1)
            ENDIF

            IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out2)
            ENDIF

            IF (present(data_r8_2d_in3) .and. present(data_r8_2d_out3)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out3)
            ENDIF

            IF (present(data_r8_2d_in4) .and. present(data_r8_2d_out4)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out4)
            ENDIF

            IF (present(data_r8_2d_in5) .and. present(data_r8_2d_out5)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out5)
            ENDIF

            IF (present(data_r8_2d_in6) .and. present(data_r8_2d_out6)) THEN
               CALL mpi_recv (rbuf_r8_1d, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_1d, msk, data_r8_2d_out6)
            ENDIF

            deallocate (rbuf_r8_1d)

            IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
               allocate (rbuf_r8_2d (n1_r8_3d_in1,nreq))
               CALL mpi_recv (rbuf_r8_2d, n1_r8_3d_in1*nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_2d, msk, data_r8_3d_out1)
               deallocate (rbuf_r8_2d)
            ENDIF

            IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
               allocate (rbuf_r8_2d (n1_r8_3d_in2,nreq))
               CALL mpi_recv (rbuf_r8_2d, n1_r8_3d_in2*nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_r8_2d, msk, data_r8_3d_out2)
               deallocate (rbuf_r8_2d)
            ENDIF

            allocate (rbuf_i4_1d (nreq))
            IF (present(data_i4_2d_in1) .and. present(data_i4_2d_out1)) THEN
               CALL mpi_recv (rbuf_i4_1d, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_i4_1d, msk, data_i4_2d_out1)
            ENDIF

            IF (present(data_i4_2d_in2) .and. present(data_i4_2d_out2)) THEN
               CALL mpi_recv (rbuf_i4_1d, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf_i4_1d, msk, data_i4_2d_out2)
            ENDIF

            deallocate (rbuf_i4_1d)

            deallocate (ibuf)
         ENDIF
      ENDDO

      deallocate (xlist)
      deallocate (ylist)
      deallocate (ipt  )
      deallocate (msk  )

#else

      DO ireq = 1, totalreq

         xblk = grid_in%xblk(xlist(ireq))
         yblk = grid_in%yblk(ylist(ireq))
         xloc = grid_in%xloc(xlist(ireq))
         yloc = grid_in%yloc(ylist(ireq))

         IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1)) THEN
            data_r8_2d_out1(ireq) = data_r8_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2)) THEN
            data_r8_2d_out2(ireq) = data_r8_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in3) .and. present(data_r8_2d_out3)) THEN
            data_r8_2d_out3(ireq) = data_r8_2d_in3%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in4) .and. present(data_r8_2d_out4)) THEN
            data_r8_2d_out4(ireq) = data_r8_2d_in4%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in5) .and. present(data_r8_2d_out5)) THEN
            data_r8_2d_out5(ireq) = data_r8_2d_in5%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in6) .and. present(data_r8_2d_out6)) THEN
            data_r8_2d_out6(ireq) = data_r8_2d_in6%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
            data_r8_3d_out1(:,ireq) = data_r8_3d_in1%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
            data_r8_3d_out2(:,ireq) = data_r8_3d_in2%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

         IF (present(data_i4_2d_in1) .and. present(data_i4_2d_out1)) THEN
            data_i4_2d_out1(ireq) = data_i4_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_i4_2d_in2) .and. present(data_i4_2d_out2)) THEN
            data_i4_2d_out2(ireq) = data_i4_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

      ENDDO

#endif

   END SUBROUTINE aggregation_request_data

#ifdef USEMPI

   SUBROUTINE aggregation_worker_done ()

   USE MOD_SPMD_Task

   IMPLICIT NONE

   integer :: smesg(2), iproc, idest

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            smesg = (/p_iam_glb, -1/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
         ENDDO
      ENDIF

   END SUBROUTINE aggregation_worker_done

#endif


   SUBROUTINE fillnan (vec, fill, defval)

   USE MOD_Precision
   USE MOD_UserDefFun, only: isnan_ud
   IMPLICIT NONE

   real(r8), intent(inout) :: vec(:)
   logical,  intent(in)    :: fill
   real(r8), intent(in)    :: defval

   ! local variables
   integer  :: i, n
   real(r8) :: s

      n = 0
      s = 0.
      DO i = lbound(vec,1), ubound(vec,1)
         IF (.not. isnan_ud(vec(i))) THEN
            n = n + 1
            s = s + vec(i)
         ENDIF
      ENDDO

      IF ((n > 0) .and. (n < size(vec))) THEN
         s = s/n
         DO i = lbound(vec,1), ubound(vec,1)
            IF (isnan_ud(vec(i))) vec(i) = s
         ENDDO
      ENDIF

      IF ((n == 0) .and. fill) THEN
         vec(:) = defval
      ENDIF

   END SUBROUTINE fillnan

END MODULE MOD_AggregationRequestData
