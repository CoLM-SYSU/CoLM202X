#include <define.h>

MODULE mod_aggregation

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

      USE precision
      USE spmd_task
      USE mod_grid
      USE mod_data_type

      IMPLICIT NONE

      TYPE (grid_type), intent(in) :: grid_in
      
      ! 2D REAL data
      TYPE (block_data_real8_2d), intent(in), optional :: data_r8_2d_in1
      TYPE (block_data_real8_2d), intent(in), optional :: data_r8_2d_in2
      TYPE (block_data_real8_2d), intent(in), optional :: data_r8_2d_in3
      TYPE (block_data_real8_2d), intent(in), optional :: data_r8_2d_in4
      TYPE (block_data_real8_2d), intent(in), optional :: data_r8_2d_in5
      TYPE (block_data_real8_2d), intent(in), optional :: data_r8_2d_in6

      ! 3D REAL data
      INTEGER, intent(in), optional :: n1_r8_3d_in1
      TYPE (block_data_real8_3d), intent(in), optional :: data_r8_3d_in1
      
      INTEGER, intent(in), optional :: n1_r8_3d_in2
      TYPE (block_data_real8_3d), intent(in), optional :: data_r8_3d_in2
      
      ! 2D INTEGER data
      TYPE (block_data_int32_2d), intent(in), optional :: data_i4_2d_in1
      TYPE (block_data_int32_2d), intent(in), optional :: data_i4_2d_in2
      
      ! Local Variables
      INTEGER :: nreq, ireq, rmesg(2), isrc, idest
      INTEGER :: xblk, yblk, xloc, yloc
      INTEGER,  allocatable :: ylist(:), xlist(:)

      REAL(r8), allocatable :: sbuf_r8_1d(:), sbuf_r8_2d(:,:)
      INTEGER , allocatable :: sbuf_i4_1d(:)

      LOGICAL,  allocatable :: worker_done (:)

      IF (p_is_io) THEN
         
         allocate (worker_done (0:p_np_worker-1))

         worker_done(:) = .false.
         DO while (any(.not. worker_done))

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
   SUBROUTINE aggregation_request_data (  &
         pixelset, iset, grid_in, area,   &
         data_r8_2d_in1, data_r8_2d_out1, &
         data_r8_2d_in2, data_r8_2d_out2, &
         data_r8_2d_in3, data_r8_2d_out3, &
         data_r8_2d_in4, data_r8_2d_out4, &
         data_r8_2d_in5, data_r8_2d_out5, &
         data_r8_2d_in6, data_r8_2d_out6, &
         data_r8_3d_in1, data_r8_3d_out1, n1_r8_3d_in1, &
         data_r8_3d_in2, data_r8_3d_out2, n1_r8_3d_in2, &
         data_i4_2d_in1, data_i4_2d_out1, &
         data_i4_2d_in2, data_i4_2d_out2, &
         filledvalue_i4)

      USE precision
      USE spmd_task
      USE mod_block
      USE mod_pixel
      USE mod_grid
      USE mod_data_type
      USE mod_mesh
      USE mod_pixelset
      USE mod_utils

      IMPLICIT NONE

      TYPE (pixelset_type), intent(in) :: pixelset
      INTEGER, intent(in)  :: iset

      TYPE (grid_type), intent(in) :: grid_in
      
      REAL(r8), allocatable, intent(out), optional :: area(:)

      TYPE (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in1
      REAL(r8), allocatable,      intent(out), optional :: data_r8_2d_out1 (:)

      TYPE (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in2
      REAL(r8), allocatable,      intent(out), optional :: data_r8_2d_out2 (:)

      TYPE (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in3
      REAL(r8), allocatable,      intent(out), optional :: data_r8_2d_out3 (:)

      TYPE (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in4
      REAL(r8), allocatable,      intent(out), optional :: data_r8_2d_out4 (:)

      TYPE (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in5
      REAL(r8), allocatable,      intent(out), optional :: data_r8_2d_out5 (:)

      TYPE (block_data_real8_2d), intent(in),  optional :: data_r8_2d_in6
      REAL(r8), allocatable,      intent(out), optional :: data_r8_2d_out6 (:)

      INTEGER, intent(in), optional :: n1_r8_3d_in1
      TYPE (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in1
      REAL(r8), allocatable,      intent(out), optional :: data_r8_3d_out1 (:,:)

      INTEGER, intent(in), optional :: n1_r8_3d_in2
      TYPE (block_data_real8_3d), intent(in),  optional :: data_r8_3d_in2
      REAL(r8), allocatable,      intent(out), optional :: data_r8_3d_out2 (:,:)
      
      TYPE (block_data_int32_2d), intent(in),  optional :: data_i4_2d_in1
      INTEGER, allocatable,       intent(out), optional :: data_i4_2d_out1 (:)

      TYPE (block_data_int32_2d), intent(in),  optional :: data_i4_2d_in2
      INTEGER, allocatable,       intent(out), optional :: data_i4_2d_out2 (:)

      INTEGER, intent(in), optional :: filledvalue_i4

      ! Local Variables
      INTEGER :: nreq, smesg(2), isrc, idest, iproc
      INTEGER :: ilon, ilat, xblk, yblk, xloc, yloc
      INTEGER :: ipxl, ie, ipxstt, ipxend, lb1
      INTEGER,  allocatable :: ylist(:), xlist(:), ipt(:), ibuf(:)

      REAL(r8), allocatable :: rbuf_r8_1d(:), rbuf_r8_2d(:,:)
      INTEGER , allocatable :: rbuf_i4_1d(:)

      LOGICAL,  allocatable :: msk(:)


      ie     = pixelset%ielm  (iset)
      ipxstt = pixelset%ipxstt(iset)
      ipxend = pixelset%ipxend(iset)

      IF (present (area)) THEN
         allocate (area (ipxstt:ipxend))
         DO ipxl = ipxstt, ipxend
            area(ipxl) = areaquad (&
               pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
               pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
         ENDDO
      ENDIF

      IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1)) THEN
         allocate (data_r8_2d_out1 (ipxstt:ipxend))
      ENDIF

      IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2)) THEN
         allocate (data_r8_2d_out2 (ipxstt:ipxend))
      ENDIF

      IF (present(data_r8_2d_in3) .and. present(data_r8_2d_out3)) THEN
         allocate (data_r8_2d_out3 (ipxstt:ipxend))
      ENDIF

      IF (present(data_r8_2d_in4) .and. present(data_r8_2d_out4)) THEN
         allocate (data_r8_2d_out4 (ipxstt:ipxend))
      ENDIF

      IF (present(data_r8_2d_in5) .and. present(data_r8_2d_out5)) THEN
         allocate (data_r8_2d_out5 (ipxstt:ipxend))
      ENDIF

      IF (present(data_r8_2d_in6) .and. present(data_r8_2d_out6)) THEN
         allocate (data_r8_2d_out6 (ipxstt:ipxend))
      ENDIF

      IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
         lb1 = data_r8_3d_in1%lb1
         allocate (data_r8_3d_out1 (lb1:lb1-1+n1_r8_3d_in1,ipxstt:ipxend))
      ENDIF

      IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
         lb1 = data_r8_3d_in2%lb1
         allocate (data_r8_3d_out2 (lb1:lb1-1+n1_r8_3d_in2,ipxstt:ipxend))
      ENDIF

      IF (present(data_i4_2d_in1) .and. present(data_i4_2d_out1)) THEN
         allocate (data_i4_2d_out1 (ipxstt:ipxend))
         IF (present(filledvalue_i4)) THEN
            data_i4_2d_out1 = filledvalue_i4
         ENDIF
      ENDIF

      IF (present(data_i4_2d_in2) .and. present(data_i4_2d_out2)) THEN
         allocate (data_i4_2d_out2 (ipxstt:ipxend))
         IF (present(filledvalue_i4)) THEN
            data_i4_2d_out2 = filledvalue_i4
         ENDIF
      ENDIF

#ifdef USEMPI

      allocate (xlist (ipxstt:ipxend))
      allocate (ylist (ipxstt:ipxend))
      allocate (ipt   (ipxstt:ipxend))
      allocate (msk   (ipxstt:ipxend))

      ipt(:) = -1

      DO ipxl = ipxstt, ipxend
         xlist(ipxl) = grid_in%xgrd(mesh(ie)%ilon(ipxl))
         ylist(ipxl) = grid_in%ygrd(mesh(ie)%ilat(ipxl))

         IF ((xlist(ipxl) < 0) .or. (ylist(ipxl) < 0)) THEN
            write(*,*) 'Warning: request data out of range.'
         ELSE
            xblk = grid_in%xblk(xlist(ipxl))
            yblk = grid_in%yblk(ylist(ipxl))
            ipt(ipxl) = gblock%pio(xblk,yblk)
         ENDIF

      ENDDO

      DO iproc = 0, p_np_io-1
         msk = (ipt == p_address_io(iproc))
         nreq = count(msk)

         IF (nreq > 0) THEN

            smesg = (/p_iam_glb, nreq/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

            allocate (ibuf (nreq))

            ibuf = pack(xlist, msk)
            CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            ibuf = pack(ylist, msk)
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
      
      DO ipxl = ipxstt, ipxend

         ilon = grid_in%xgrd(mesh(ie)%ilon(ipxl))
         ilat = grid_in%ygrd(mesh(ie)%ilat(ipxl))

         IF ((ilon < 0) .or. (ilat < 0)) THEN
            write(*,*) 'Warning: request data out of range.'
            cycle
         ENDIF

         xblk = grid_in%xblk(ilon)
         yblk = grid_in%yblk(ilat)
         xloc = grid_in%xloc(ilon)
         yloc = grid_in%yloc(ilat)

         IF (present(data_r8_2d_in1) .and. present(data_r8_2d_out1)) THEN
            data_r8_2d_out1(ipxl) = data_r8_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in2) .and. present(data_r8_2d_out2)) THEN
            data_r8_2d_out2(ipxl) = data_r8_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in3) .and. present(data_r8_2d_out3)) THEN
            data_r8_2d_out3(ipxl) = data_r8_2d_in3%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in4) .and. present(data_r8_2d_out4)) THEN
            data_r8_2d_out4(ipxl) = data_r8_2d_in4%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in5) .and. present(data_r8_2d_out5)) THEN
            data_r8_2d_out5(ipxl) = data_r8_2d_in5%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_2d_in6) .and. present(data_r8_2d_out6)) THEN
            data_r8_2d_out6(ipxl) = data_r8_2d_in6%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in1) .and. present(data_r8_3d_out1) .and. present(n1_r8_3d_in1)) THEN
            data_r8_3d_out1(:,ipxl) = data_r8_3d_in1%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

         IF (present(data_r8_3d_in2) .and. present(data_r8_3d_out2) .and. present(n1_r8_3d_in2)) THEN
            data_r8_3d_out2(:,ipxl) = data_r8_3d_in2%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

         IF (present(data_i4_2d_in1) .and. present(data_i4_2d_out1)) THEN
            data_i4_2d_out1(ipxl) = data_i4_2d_in1%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_i4_2d_in2) .and. present(data_i4_2d_out2)) THEN
            data_i4_2d_out2(ipxl) = data_i4_2d_in2%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF


      ENDDO
      
#endif

   END SUBROUTINE aggregation_request_data

#ifdef USEMPI

   SUBROUTINE aggregation_worker_done ()

      USE spmd_task

      IMPLICIT NONE 

      INTEGER :: smesg(2), iproc, idest
   
      IF (p_is_worker) THEN 
         DO iproc = 0, p_np_io-1
            smesg = (/p_iam_glb, -1/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
         ENDDO
      ENDIF

   END SUBROUTINE aggregation_worker_done

#endif

END MODULE mod_aggregation
