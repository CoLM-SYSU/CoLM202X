#include <define.h>

MODULE mod_aggregation_generic

   IMPLICIT NONE

   PUBLIC :: aggregation_gen_request_data

#ifdef USEMPI
   PUBLIC :: aggregation_gen_data_daemon
   PUBLIC :: aggregation_gen_worker_done
#endif

CONTAINS
! ---- subroutines ----
   
#ifdef USEMPI
   SUBROUTINE aggregation_gen_data_daemon (grid_in, data_r8, data_i4, data_r8_3d, ndim1)

      USE precision
      USE spmd_task
      USE mod_grid
      USE mod_data_type

      IMPLICIT NONE

      TYPE (grid_type),           intent(in) :: grid_in
      TYPE (block_data_real8_2d), intent(in), optional :: data_r8
      TYPE (block_data_int32_2d), intent(in), optional :: data_i4
      TYPE (block_data_real8_3d), intent(in), optional :: data_r8_3d
      INTEGER, optional :: ndim1
      
      ! Local Variables
      INTEGER :: nreq, ireq, rmesg(2), isrc, idest
      INTEGER :: xblk, yblk, xloc, yloc
      INTEGER,  allocatable :: ylist(:), xlist(:)
      LOGICAL,  allocatable :: worker_done (:)
      
      REAL(r8), allocatable :: buff_r8(:)
      INTEGER , allocatable :: buff_i4(:)
      REAL(r8), allocatable :: buff_r8_3d(:,:)

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
               
               IF (present(data_r8)) THEN
                  allocate (buff_r8 (nreq))
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     buff_r8(ireq) = data_r8%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  idest = isrc
                  CALL mpi_send (buff_r8, nreq, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)

                  deallocate (buff_r8)
               ENDIF

               IF (present(data_i4)) THEN
                  allocate (buff_i4 (nreq))
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     buff_i4(ireq) = data_i4%blk(xblk,yblk)%val(xloc,yloc)
                  ENDDO

                  idest = isrc
                  CALL mpi_send (buff_i4, nreq, MPI_INTEGER4, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
                  
                  deallocate (buff_i4)
               ENDIF

               IF (present(data_r8_3d) .and. present(ndim1)) THEN
                  allocate (buff_r8_3d (ndim1,nreq))
                  DO ireq = 1, nreq
                     xblk = grid_in%xblk(xlist(ireq))
                     yblk = grid_in%yblk(ylist(ireq))
                     xloc = grid_in%xloc(xlist(ireq))
                     yloc = grid_in%yloc(ylist(ireq))

                     buff_r8_3d(:,ireq) = data_r8_3d%blk(xblk,yblk)%val(:,xloc,yloc)
                  ENDDO

                  idest = isrc
                  CALL mpi_send (buff_r8_3d, nreq*ndim1, MPI_REAL8, &
                     idest, mpi_tag_data, p_comm_glb, p_err)

                  deallocate (buff_r8_3d)
               ENDIF

               deallocate (ylist)
               deallocate (xlist)

            ELSE
               worker_done(p_itis_worker(isrc)) = .true.
            ENDIF

         ENDDO

         deallocate (worker_done)

      ENDIF 

   END SUBROUTINE aggregation_gen_data_daemon

#endif

   !----------------------------------------------------
   SUBROUTINE aggregation_gen_request_data ( &
         grid_in, xlist, ylist, data_r8, out_r8, data_i4, out_i4, &
         data_r8_3d, out_r8_3d, ndim1)

      USE precision
      USE spmd_task
      USE mod_block
      USE mod_pixel
      USE mod_grid
      USE mod_data_type
      USE mod_mesh
      USE mod_pixelset
      USE mod_utils
      USE GlobalVars, only : spval

      IMPLICIT NONE

      TYPE(grid_type), intent(in) :: grid_in
      INTEGER,         intent(in) :: xlist(:), ylist(:)

      TYPE (block_data_real8_2d), intent(in),  optional :: data_r8
      REAL(r8), allocatable,      intent(out), optional :: out_r8(:)

      TYPE (block_data_int32_2d), intent(in),  optional :: data_i4
      INTEGER, allocatable,       intent(out), optional :: out_i4(:)

      TYPE (block_data_real8_3d), intent(in),  optional :: data_r8_3d
      REAL(r8), allocatable,      intent(out), optional :: out_r8_3d(:,:)
      INTEGER, optional :: ndim1

      ! Local Variables
      INTEGER :: nreq, smesg(2), isrc, idest, iproc
      INTEGER :: ilon, ilat, xblk, yblk, xloc, yloc
      INTEGER :: npxl, ipxl
      INTEGER,  allocatable :: xlist_g(:), ylist_g(:)
      INTEGER,  allocatable :: ipt(:), ibuf(:)
      LOGICAL,  allocatable :: msk(:)
      
      REAL(r8), allocatable :: buff_r8(:)
      INTEGER , allocatable :: buff_i4(:)
      REAL(r8), allocatable :: buff_r8_3d(:,:)


      npxl = size(xlist)

      IF (present(data_r8) .and. present(out_r8)) THEN
         allocate (out_r8(npxl))
         out_r8(:) = spval
      ENDIF

      IF (present(data_i4) .and. present(out_i4)) THEN
         allocate (out_i4(npxl))
         out_i4(:) = -1
      ENDIF

      IF (present(data_r8_3d) .and. present(out_r8_3d) .and. present(ndim1)) THEN
         allocate (out_r8_3d(ndim1,npxl))
         out_r8_3d(:,:) = spval
      ENDIF

#ifdef USEMPI

      allocate (xlist_g (npxl))
      allocate (ylist_g (npxl))

      allocate (ipt (npxl))
      allocate (msk (npxl))

      ipt(:) = -1

      DO ipxl = 1, npxl
         xlist_g(ipxl) = grid_in%xgrd(xlist(ipxl))
         ylist_g(ipxl) = grid_in%ygrd(ylist(ipxl))
         IF ((xlist_g(ipxl) > 0) .and. (ylist_g(ipxl) > 0)) THEN
            xblk = grid_in%xblk(xlist_g(ipxl))
            yblk = grid_in%yblk(ylist_g(ipxl))
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

            ibuf = pack(xlist_g, msk)
            CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            ibuf = pack(ylist_g, msk)
            CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

            isrc = idest

            IF (present(data_r8) .and. present(out_r8)) THEN
               allocate (buff_r8 (nreq))
               CALL mpi_recv (buff_r8, nreq, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (buff_r8, msk, out_r8)
               deallocate (buff_r8)
            ENDIF

            IF (present(data_i4) .and. present(out_i4)) THEN
               allocate (buff_i4 (nreq))
               CALL mpi_recv (buff_i4, nreq, MPI_INTEGER4, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (buff_i4, msk, out_i4)
               deallocate (buff_i4)
            ENDIF

            IF (present(data_r8_3d) .and. present(out_r8_3d) .and. present(ndim1)) THEN
               allocate (buff_r8_3d (ndim1,nreq))
               CALL mpi_recv (buff_r8_3d, nreq*ndim1, MPI_REAL8, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL unpack_inplace (buff_r8_3d, msk, out_r8_3d)
               deallocate (buff_r8_3d)
            ENDIF

            deallocate (ibuf)
         ENDIF
      ENDDO

      deallocate (xlist_g)
      deallocate (ylist_g)
      deallocate (ipt)
      deallocate (msk)

#else
      
      DO ipxl = 1, npxl

         ilon = grid_in%xgrd(xlist(ipxl))
         ilat = grid_in%ygrd(ylist(ipxl))

         IF ((ilon > 0) .and. (ilat > 0)) THEN
            xblk = grid_in%xblk(ilon)
            yblk = grid_in%yblk(ilat)
            xloc = grid_in%xloc(ilon)
            yloc = grid_in%yloc(ilat)

            IF (present(data_r8) .and. present(out_r8)) THEN
               out_r8(ipxl) = data_r8%blk(xblk,yblk)%val(xloc,yloc)
            ENDIF

            IF (present(data_i4) .and. present(out_i4)) THEN
               out_i4(ipxl) = data_i4%blk(xblk,yblk)%val(xloc,yloc)
            ENDIF

            IF (present(data_r8_3d) .and. present(out_r8_3d) .and. present(ndim1)) THEN
               out_r8_3d(:,ipxl) = data_r8_3d%blk(xblk,yblk)%val(:,xloc,yloc)
            ENDIF
         ENDIF
      ENDDO
      
#endif

   END SUBROUTINE aggregation_gen_request_data


#ifdef USEMPI

   SUBROUTINE aggregation_gen_worker_done ()

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

   END SUBROUTINE aggregation_gen_worker_done

#endif

END MODULE mod_aggregation_generic
