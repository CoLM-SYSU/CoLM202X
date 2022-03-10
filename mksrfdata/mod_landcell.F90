#include <define.h>

MODULE mod_landcell

   USE precision
   USE mod_pixelset
   USE mod_grid
   IMPLICIT NONE

   ! ---- Instance ----
   INTEGER :: numcell
   TYPE(pixelset_type) :: landcell

#ifdef CATCHMENT 
   TYPE(grid_type) :: ghband
#endif

CONTAINS
   
#ifdef GRIDBASED
   ! -------------------------------
   SUBROUTINE landcell_build 

      USE precision
      USE spmd_task
      USE mod_landunit
      IMPLICIT NONE

      ! Local Variables
      INTEGER :: icell, ncell_glb

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land celles :'
      ENDIF

      IF (p_is_worker) THEN

         numcell = numunit

         allocate (landcell%unum (numcell))
         allocate (landcell%iunt (numcell))
         allocate (landcell%istt (numcell))
         allocate (landcell%iend (numcell))
         allocate (landcell%ltyp (numcell))

         DO icell = 1, numcell
            landcell%unum(icell) = landunit(icell)%num
            landcell%iunt(icell) = icell
            landcell%istt(icell) = 1 
            landcell%iend(icell) = landunit(icell)%npxl 
         ENDDO

      ENDIF

      landcell%nset = numcell
      CALL landcell%set_vecgs 

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         CALL mpi_reduce (numcell, ncell_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', ncell_glb, ' cells on worker.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numunit, ' cells.'
#endif

   END SUBROUTINE landcell_build
#endif


#ifdef CATCHMENT 
   ! -------------------------------
   SUBROUTINE landcell_build ()

      USE precision
      USE spmd_task
      USE mod_utils
      USE mod_block
      USE mod_pixel
      USE mod_grid
      USE mod_data_type
      USE mod_landunit
      USE mod_hydro_data
      USE mod_namelist

      IMPLICIT NONE

      ! Local Variables
      TYPE (block_data_int32_2d) :: hydrdata
      INTEGER :: nreq, rmesg(2), smesg(2), isrc, idest, iproc
      INTEGER :: ilon_g, ilat_g, iloc
      INTEGER :: ig, xblk, yblk, xloc, yloc
      INTEGER :: npxl, ipxl
      INTEGER :: iu
      INTEGER, allocatable :: xlist(:), ylist(:), sbuf(:), rbuf(:)
      INTEGER, allocatable :: ltype(:), ipt(:), order(:)
      INTEGER, allocatable :: unum_tmp(:), ltyp_tmp(:), istt_tmp(:), iend_tmp(:), iunt_tmp(:)
      LOGICAL, allocatable :: msk(:)
      LOGICAL, allocatable :: worker_done(:)
      INTEGER :: ncell_glb

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land celles :'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN
         CALL allocate_block_data (ghband, hydrdata)
         CALL hydro_data_read (DEF_dir_hydrodata, 'hunit', ghband, hydrdata)
      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN 

         allocate (worker_done (0:p_np_worker-1))
         worker_done(:) = .false.

         DO while (.true.)

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc = rmesg(1)
            nreq = rmesg(2)

            IF (nreq > 0) THEN

               allocate (xlist (nreq))
               allocate (ylist (nreq))
               allocate (sbuf  (nreq))

               CALL mpi_recv (xlist, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (ylist, nreq, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO ig = 1, nreq
                  xblk = ghband%xblk(xlist(ig))
                  yblk = ghband%yblk(ylist(ig))
                  xloc = ghband%xloc(xlist(ig))
                  yloc = ghband%yloc(ylist(ig))

                  sbuf(ig) = hydrdata%blk(xblk,yblk)%val(xloc,yloc)
               ENDDO
            
               idest = isrc
               CALL mpi_send (sbuf, nreq, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (ylist)
               deallocate (xlist)
               deallocate (sbuf )

            ELSE
               worker_done(p_itis_worker(isrc)) = .true.
            ENDIF

            IF (all(worker_done)) exit
         ENDDO

         deallocate (worker_done)
      ENDIF 
#endif

      IF (p_is_worker) THEN

         allocate (unum_tmp (numunit*DEF_max_hband))
         allocate (iunt_tmp (numunit*DEF_max_hband))
         allocate (ltyp_tmp (numunit*DEF_max_hband))
         allocate (istt_tmp (numunit*DEF_max_hband))
         allocate (iend_tmp (numunit*DEF_max_hband))

         numcell = 0

         DO iu = 1, numunit
         
            npxl = landunit(iu)%npxl 
            
            allocate (ltype (1:npxl))

#ifdef USEMPI
            allocate (xlist (1:npxl))
            allocate (ylist (1:npxl))
            allocate (ipt   (1:npxl))
            allocate (msk   (1:npxl))

            DO ipxl = 1, npxl
               ilon_g = ghband%xgrd(landunit(iu)%ilon(ipxl))
               ilat_g = ghband%ygrd(landunit(iu)%ilat(ipxl))
               xlist(ipxl) = ilon_g
               ylist(ipxl) = ilat_g

               xblk = ghband%xblk(xlist(ipxl))
               yblk = ghband%yblk(ylist(ipxl))
               ipt(ipxl) = gblock%pio(xblk,yblk)
            ENDDO

            DO iproc = 0, p_np_io-1
               msk = (ipt == p_address_io(iproc))
               nreq = count(msk)

               IF (nreq > 0) THEN
               
                  smesg = (/p_iam_glb, nreq/)
                  idest = p_address_io(iproc)
                  CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

                  allocate (sbuf (nreq))
                  allocate (rbuf (nreq))

                  sbuf = pack(xlist, msk)
                  CALL mpi_send (sbuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

                  sbuf = pack(ylist, msk)
                  CALL mpi_send (sbuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

                  isrc = idest
                  CALL mpi_recv (rbuf, nreq, MPI_INTEGER, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

                  CALL unpack_inplace (rbuf, msk, ltype)

                  deallocate (sbuf)
                  deallocate (rbuf)
               ENDIF
            ENDDO
            
            deallocate (xlist)
            deallocate (ylist)
            deallocate (ipt  )
            deallocate (msk  )
#else
            DO ipxl = 1, npxl
               ilon_g = ghband%xgrd(landunit(iu)%ilon(ipxl))
               ilat_g = ghband%ygrd(landunit(iu)%ilat(ipxl))
               xblk   = ghband%xblk(ilon_g)
               yblk   = ghband%yblk(ilat_g)
               xloc   = ghband%xloc(ilon_g)
               yloc   = ghband%yloc(ilat_g)

               ltype(ipxl) = hydrdata%blk(xblk,yblk)%val(xloc,yloc)
            ENDDO
#endif
               
            allocate (order (1:npxl))
            order = (/ (ipxl, ipxl = 1, npxl) /)

            CALL quicksort (npxl, ltype, order)
               
            landunit(iu)%ilon(1:npxl) = landunit(iu)%ilon(order)
            landunit(iu)%ilat(1:npxl) = landunit(iu)%ilat(order)
            
            DO ipxl = 1, npxl
               IF (ipxl == 1) THEN
                  numcell = numcell + 1 
                  unum_tmp (numcell) = landunit(iu)%num
                  iunt_tmp (numcell) = iu
                  ltyp_tmp (numcell) = ltype(ipxl)
                  istt_tmp (numcell) = ipxl
               ELSEIF (ltype(ipxl) /= ltype(ipxl-1)) THEN
                  iend_tmp(numcell) = ipxl - 1

                  numcell = numcell + 1
                  unum_tmp (numcell) = landunit(iu)%num
                  iunt_tmp (numcell) = iu
                  ltyp_tmp (numcell) = ltype(ipxl)
                  istt_tmp (numcell) = ipxl
               ENDIF
            ENDDO
            iend_tmp(numcell) = npxl
            
            deallocate (ltype)
            deallocate (order)

         ENDDO

         allocate (landcell%unum (numcell))
         allocate (landcell%iunt (numcell))
         allocate (landcell%ltyp (numcell))
         allocate (landcell%istt (numcell))
         allocate (landcell%iend (numcell))
         
         landcell%unum = unum_tmp(1:numcell)  
         landcell%iunt = iunt_tmp(1:numcell)  
         landcell%ltyp = ltyp_tmp (1:numcell)  
         landcell%istt = istt_tmp (1:numcell)
         landcell%iend = iend_tmp (1:numcell)

         deallocate (ltyp_tmp)
         deallocate (istt_tmp)
         deallocate (iend_tmp)
         deallocate (iunt_tmp)
         deallocate (unum_tmp)

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            smesg = (/p_iam_glb, -1/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
         ENDDO
#endif


      ENDIF

      landcell%nset = numcell
      CALL landcell%set_vecgs 
         
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         CALL mpi_reduce (numcell, ncell_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', ncell_glb, ' cells on worker.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numunit, ' cells.'
#endif

   END SUBROUTINE landcell_build
#endif

END MODULE mod_landcell
