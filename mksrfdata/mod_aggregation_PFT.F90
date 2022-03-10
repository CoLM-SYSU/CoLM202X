#include <define.h>

MODULE mod_aggregation_pft

   USE GlobalVars, only : N_PFT
   IMPLICIT NONE

   PUBLIC :: aggregation_pft_request_data

#ifdef USEMPI
   PUBLIC :: aggregation_pft_data_daemon
   PUBLIC :: aggregation_pft_worker_done
#endif

! ---- subroutines ----
CONTAINS
   
#ifdef USEMPI
   SUBROUTINE aggregation_pft_data_daemon (grid, pctpft, data2, data3)

      USE precision
      USE spmd_task
      USE mod_grid
      USE mod_data_type

      IMPLICIT NONE

      TYPE (grid_type),           intent(in) :: grid
      TYPE (block_data_real8_3d), intent(in) :: pctpft

      TYPE (block_data_real8_2d), intent(in), optional :: data2
      TYPE (block_data_real8_3d), intent(in), optional :: data3
      
      ! Local Variables
      INTEGER :: nreq, ireq, rmesg(2), isrc, idest
      INTEGER :: xblk, yblk, xloc, yloc
      INTEGER,  allocatable :: ylist(:), xlist(:)
      REAL(r8), allocatable :: pbuff(:,:), sbuff2(:), sbuff3(:,:)
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
               
               allocate (pbuff (N_PFT,nreq))

               DO ireq = 1, nreq
                  xblk = grid%xblk(xlist(ireq))
                  yblk = grid%yblk(ylist(ireq))
                  xloc = grid%xloc(xlist(ireq))
                  yloc = grid%yloc(ylist(ireq))

                  pbuff(:,ireq) = pctpft%blk(xblk,yblk)%val(:,xloc,yloc)
               ENDDO

               idest = isrc
               CALL mpi_send (pbuff, N_PFT*nreq, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate(pbuff)

               IF (present(data2)) THEN
                  allocate (sbuff2(nreq))
                  DO ireq = 1, nreq
                     xblk = grid%xblk(xlist(ireq))
                     yblk = grid%yblk(ylist(ireq))
                     xloc = grid%xloc(xlist(ireq))
                     yloc = grid%yloc(ylist(ireq))

                     sbuff2(ireq) = data2%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  idest = isrc
                  CALL mpi_send (sbuff2, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)

                  deallocate(sbuff2)
               ENDIF

               IF (present(data3)) THEN
                  allocate (sbuff3(N_PFT,nreq))
                  DO ireq = 1, nreq
                     xblk = grid%xblk(xlist(ireq))
                     yblk = grid%yblk(ylist(ireq))
                     xloc = grid%xloc(xlist(ireq))
                     yloc = grid%yloc(ylist(ireq))

                     sbuff3(:,ireq) = data3%blk(xblk,yblk)%val(:,xloc,yloc)
                  ENDDO

                  idest = isrc
                  CALL mpi_send (sbuff3, N_PFT*nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)

                  deallocate(sbuff3)
               ENDIF

               deallocate (ylist)
               deallocate (xlist)
            ELSE
               worker_done(p_itis_worker(isrc)) = .true.
            ENDIF

         ENDDO

         deallocate (worker_done)

      ENDIF 

   END SUBROUTINE aggregation_pft_data_daemon

#endif

   !----------------------------------------------------
   SUBROUTINE aggregation_pft_request_data ( &
         ipatch, grid, pctpft, pctout, area, data2, dout2, data3, dout3)

      USE precision
      USE spmd_task
      USE mod_block
      USE mod_pixel
      USE mod_grid
      USE mod_data_type
      USE mod_landunit
      USE mod_pixelset
      USE mod_landpatch
      USE mod_utils

      IMPLICIT NONE

      INTEGER, intent(in) :: ipatch
      TYPE (grid_type),           intent(in) :: grid
      TYPE (block_data_real8_3d), intent(in) :: pctpft
      REAL(r8), allocatable, intent(out) :: pctout(:,:)
      REAL(r8), allocatable, intent(out), optional :: area(:)
      
      TYPE (block_data_real8_2d), intent(in), optional :: data2
      REAL(r8), allocatable, intent(out), optional :: dout2(:)

      TYPE (block_data_real8_3d), intent(in), optional :: data3
      REAL(r8), allocatable, intent(out), optional :: dout3(:,:)

      ! Local Variables
      INTEGER :: nreq, smesg(2), isrc, idest, iproc
      INTEGER :: ilon, ilat, xblk, yblk, xloc, yloc
      INTEGER :: npxl, ipxl, iu, istt, iend
      INTEGER,  allocatable :: ylist(:), xlist(:), ipt(:), ibuf(:)
      REAL(r8), allocatable :: rbuf2(:), rbuf3(:,:)
      LOGICAL,  allocatable :: msk(:)


      iu   = landpatch%iunt(ipatch)
      istt = landpatch%istt(ipatch)
      iend = landpatch%iend(ipatch)

      npxl = iend - istt + 1
      allocate (pctout (0:N_PFT-1,npxl))

      IF (present(data2)) allocate (dout2 (istt:iend))
      IF (present(data3)) allocate (dout3 (0:N_PFT-1,istt:iend))
      IF (present (area)) allocate (area (istt:iend))

#ifdef USEMPI

      allocate (xlist (istt:iend))
      allocate (ylist (istt:iend))
      allocate (ipt   (istt:iend))
      allocate (msk   (istt:iend))

      DO ipxl = istt, iend
         xlist(ipxl) = grid%xgrd(landunit(iu)%ilon(ipxl))
         ylist(ipxl) = grid%ygrd(landunit(iu)%ilat(ipxl))

         xblk = grid%xblk(xlist(ipxl))
         yblk = grid%yblk(ylist(ipxl))
         ipt(ipxl) = gblock%pio(xblk,yblk)
      ENDDO

      DO iproc = 0, p_np_io-1
         msk = (ipt == p_address_io(iproc))
         nreq = count(msk)

         IF (nreq > 0) THEN

            smesg = (/p_iam_glb, nreq/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

            allocate (ibuf       (nreq))
            allocate (rbuf3 (N_PFT,nreq))

            ibuf = pack(xlist, msk)
            CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            ibuf = pack(ylist, msk)
            CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            isrc = idest
            CALL mpi_recv (rbuf3, N_PFT*nreq, MPI_REAL8, &
               isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

            CALL unpack_inplace (rbuf3, msk, pctout)

            IF (present(data2) .and. present(dout2)) THEN
               allocate(rbuf2(nreq))
               CALL mpi_recv (rbuf2, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               CALL unpack_inplace (rbuf2, msk, dout2)
               deallocate(rbuf2)
            ENDIF

            IF (present(data3) .and. present(dout3)) THEN
               CALL mpi_recv (rbuf3, N_PFT*nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               CALL unpack_inplace (rbuf3, msk, dout3)
            ENDIF

            deallocate (ibuf )
            deallocate (rbuf3)
         ENDIF
      ENDDO

      deallocate (xlist)
      deallocate (ylist)
      deallocate (ipt  )
      deallocate (msk  )

#else
      
      DO ipxl = istt, iend

         ilon = grid%xgrd(landunit(iu)%ilon(ipxl))
         ilat = grid%ygrd(landunit(iu)%ilat(ipxl))
         xblk = grid%xblk(ilon)
         yblk = grid%yblk(ilat)
         xloc = grid%xloc(ilon)
         yloc = grid%yloc(ilat)

         pctout(:,ipxl) = pctpft%blk(xblk,yblk)%val(:,xloc,yloc)

         IF (present(data2) .and. present(dout2)) then
            dout2(ipxl) = data2%blk(xblk,yblk)%val(xloc,yloc)
         ENDIF

         IF (present(data3) .and. present(dout3)) then
            dout3(:,ipxl) = data3%blk(xblk,yblk)%val(:,xloc,yloc)
         ENDIF

      ENDDO
      
#endif

      pctout = max(pctout, 0.)

      IF (present(area)) THEN
         DO ipxl = istt, iend
            area(ipxl) = areaquad (&
               pixel%lat_s(landunit(iu)%ilat(ipxl)), pixel%lat_n(landunit(iu)%ilat(ipxl)), &
               pixel%lon_w(landunit(iu)%ilon(ipxl)), pixel%lon_e(landunit(iu)%ilon(ipxl)) )
         ENDDO
      ENDIF

   END SUBROUTINE aggregation_pft_request_data


#ifdef USEMPI

   SUBROUTINE aggregation_pft_worker_done ()

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

   END SUBROUTINE aggregation_pft_worker_done

#endif

END MODULE mod_aggregation_pft
