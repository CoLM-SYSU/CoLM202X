#include <define.h>

MODULE mod_aggregation_lc

   IMPLICIT NONE

   PUBLIC :: aggregation_lc_request_data

#ifdef USEMPI
   PUBLIC :: aggregation_lc_data_daemon
   PUBLIC :: aggregation_lc_worker_done
#endif

! ---- subroutines ----
CONTAINS
   
#ifdef USEMPI
   SUBROUTINE aggregation_lc_data_daemon (grid_in, &
         data_in1, data_in2, data_in3, data_in4)

      USE precision
      USE spmd_task
      USE mod_grid
      USE mod_data_type

      IMPLICIT NONE

      TYPE (grid_type),           intent(in) :: grid_in
      TYPE (block_data_real8_2d), intent(in) :: data_in1
      TYPE (block_data_real8_2d), intent(in), optional :: data_in2
      TYPE (block_data_real8_2d), intent(in), optional :: data_in3
      TYPE (block_data_real8_2d), intent(in), optional :: data_in4
      
      ! Local Variables
      INTEGER :: nreq, ireq, rmesg(2), isrc, idest
      INTEGER :: xblk, yblk, xloc, yloc
      INTEGER,  allocatable :: ylist(:), xlist(:)
      REAL(r8), allocatable :: sbuff(:)
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
               
               allocate (sbuff (nreq))

               DO ireq = 1, nreq
                  xblk = grid_in%xblk(xlist(ireq))
                  yblk = grid_in%yblk(ylist(ireq))
                  xloc = grid_in%xloc(xlist(ireq))
                  yloc = grid_in%yloc(ylist(ireq))

                  sbuff(ireq) = data_in1%blk(xblk,yblk)%val(xloc,yloc)
               ENDDO

               idest = isrc
               CALL mpi_send (sbuff, nreq, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
               
               IF (present(data_in2)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuff(ireq) = data_in2%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  idest = isrc
                  CALL mpi_send (sbuff, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_in3)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuff(ireq) = data_in3%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  idest = isrc
                  CALL mpi_send (sbuff, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               IF (present(data_in4)) THEN
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     sbuff(ireq) = data_in4%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  idest = isrc
                  CALL mpi_send (sbuff, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF

               deallocate (ylist)
               deallocate (xlist)
               deallocate (sbuff)

            ELSE
               worker_done(p_itis_worker(isrc)) = .true.
            ENDIF

         ENDDO

         deallocate (worker_done)

      ENDIF 

   END SUBROUTINE aggregation_lc_data_daemon

#endif

   !----------------------------------------------------
   SUBROUTINE aggregation_lc_request_data (ipatch, grid_in, data_in1, data_out1, areall, &
         data_in2, data_out2, data_in3, data_out3, data_in4, data_out4)

      USE precision
      USE spmd_task
      USE mod_block
      USE mod_pixel
      USE mod_grid
      USE mod_data_type
      USE mod_landbasin
      USE mod_pixelset
      USE mod_landpatch
      USE mod_utils

      IMPLICIT NONE

      INTEGER, intent(in) :: ipatch
      TYPE (grid_type),           intent(in) :: grid_in
      TYPE (block_data_real8_2d), intent(in) :: data_in1
      REAL(r8), allocatable, intent(out) :: data_out1(:)
      REAL(r8), allocatable, intent(out), optional :: areall(:)

      TYPE (block_data_real8_2d), intent(in), optional :: data_in2
      REAL(r8), allocatable, intent(out), optional :: data_out2(:)

      TYPE (block_data_real8_2d), intent(in), optional :: data_in3
      REAL(r8), allocatable, intent(out), optional :: data_out3(:)

      TYPE (block_data_real8_2d), intent(in), optional :: data_in4
      REAL(r8), allocatable, intent(out), optional :: data_out4(:)

      ! Local Variables
      INTEGER :: nreq, smesg(2), isrc, idest, iproc
      INTEGER :: ilon, ilat, xblk, yblk, xloc, yloc
      INTEGER :: npxl, ipxl, iu, istt, iend
      INTEGER,  allocatable :: ylist(:), xlist(:), ipt(:), ibuf(:)
      REAL(r8), allocatable :: rbuf(:)
      LOGICAL,  allocatable :: msk(:)


      iu   = landpatch%iunt(ipatch)
      istt = landpatch%istt(ipatch)
      iend = landpatch%iend(ipatch)

      npxl = iend - istt + 1
      allocate (data_out1 (istt:iend))

      IF (present (areall)) THEN
         allocate (areall (istt:iend))
      ENDIF

      IF (present(data_in2) .and. present(data_out2)) THEN
         allocate (data_out2 (istt:iend))
      ENDIF

      IF (present(data_in3) .and. present(data_out3)) THEN
         allocate (data_out3 (istt:iend))
      ENDIF

      IF (present(data_in4) .and. present(data_out4)) THEN
         allocate (data_out4 (istt:iend))
      ENDIF

#ifdef USEMPI

      allocate (xlist (istt:iend))
      allocate (ylist (istt:iend))
      allocate (ipt   (istt:iend))
      allocate (msk   (istt:iend))

      DO ipxl = istt, iend
         xlist(ipxl) = grid_in%xgrd(landbasin(iu)%ilon(ipxl))
         ylist(ipxl) = grid_in%ygrd(landbasin(iu)%ilat(ipxl))

         xblk = grid_in%xblk(xlist(ipxl))
         yblk = grid_in%yblk(ylist(ipxl))
         ipt(ipxl) = gblock%pio(xblk,yblk)
      ENDDO

      DO iproc = 0, p_np_io-1
         msk = (ipt == p_address_io(iproc))
         nreq = count(msk)

         IF (nreq > 0) THEN

            smesg = (/p_iam_glb, nreq/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

            allocate (ibuf (nreq))
            allocate (rbuf (nreq))

            ibuf = pack(xlist, msk)
            CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            ibuf = pack(ylist, msk)
            CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            isrc = idest

            CALL mpi_recv (rbuf, nreq, MPI_REAL8, &
               isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
            CALL unpack_inplace (rbuf, msk, data_out1)

            IF (present(data_in2) .and. present(data_out2)) THEN
               CALL mpi_recv (rbuf, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf, msk, data_out2)
            ENDIF

            IF (present(data_in3) .and. present(data_out3)) THEN
               CALL mpi_recv (rbuf, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf, msk, data_out3)
            ENDIF

            IF (present(data_in4) .and. present(data_out4)) THEN
               CALL mpi_recv (rbuf, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (rbuf, msk, data_out4)
            ENDIF

            deallocate (ibuf)
            deallocate (rbuf)
         ENDIF
      ENDDO

      deallocate (xlist)
      deallocate (ylist)
      deallocate (ipt  )
      deallocate (msk  )

#else
      
      DO ipxl = istt, iend

         ilon = grid_in%xgrd(landbasin(iu)%ilon(ipxl))
         ilat = grid_in%ygrd(landbasin(iu)%ilat(ipxl))
         xblk = grid_in%xblk(ilon)
         yblk = grid_in%yblk(ilat)
         xloc = grid_in%xloc(ilon)
         yloc = grid_in%yloc(ilat)

         data_out1(ipxl) = data_in1%blk(xblk,yblk)%val(xloc,yloc)
            
         IF (present(data_in2) .and. present(data_out2)) THEN
            data_out2(ipxl) = data_in2%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_in3) .and. present(data_out3)) THEN
            data_out3(ipxl) = data_in3%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data_in4) .and. present(data_out4)) THEN
            data_out4(ipxl) = data_in4%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

      ENDDO
      
#endif

      IF (present(areall)) THEN
         DO ipxl = istt, iend
            areall(ipxl) = areaquad (&
               pixel%lat_s(landbasin(iu)%ilat(ipxl)), pixel%lat_n(landbasin(iu)%ilat(ipxl)), &
               pixel%lon_w(landbasin(iu)%ilon(ipxl)), pixel%lon_e(landbasin(iu)%ilon(ipxl)) )
         ENDDO
      ENDIF

   END SUBROUTINE aggregation_lc_request_data


#ifdef USEMPI

   SUBROUTINE aggregation_lc_worker_done ()

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

   END SUBROUTINE aggregation_lc_worker_done

#endif

END MODULE mod_aggregation_lc
