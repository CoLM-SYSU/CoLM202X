#include <define.h>

MODULE mod_hydrounit

   USE precision
   USE mod_pixelset
   USE mod_grid
   IMPLICIT NONE

   ! ---- Instance ----
   INTEGER :: numhru
   TYPE(pixelset_type) :: hydrounit

#ifdef CATCHMENT 
   TYPE(grid_type) :: ghydru
#endif

CONTAINS

#if (defined GRIDBASED || defined SinglePoint)
   ! -------------------------------
   SUBROUTINE hydrounit_build 

      USE precision
      USE spmd_task
      USE mod_landbasin
      IMPLICIT NONE

      ! Local Variables
      INTEGER :: ihru, nhru_glb

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land hydro units:'
      ENDIF

      IF (p_is_worker) THEN

         numhru = numbasin

         allocate (hydrounit%unum (numhru))
         allocate (hydrounit%iunt (numhru))
         allocate (hydrounit%istt (numhru))
         allocate (hydrounit%iend (numhru))
         allocate (hydrounit%ltyp (numhru))

         DO ihru = 1, numhru
            hydrounit%unum(ihru) = landbasin(ihru)%num
            hydrounit%iunt(ihru) = ihru
            hydrounit%istt(ihru) = 1 
            hydrounit%iend(ihru) = landbasin(ihru)%npxl 
         ENDDO

      ENDIF

      hydrounit%nset = numhru
      CALL hydrounit%set_vecgs 

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         CALL mpi_reduce (numhru, nhru_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', nhru_glb, ' hydro units on worker.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numbasin, ' hydro units.'
#endif

   END SUBROUTINE hydrounit_build
#endif


#ifdef CATCHMENT 
   ! -------------------------------
   SUBROUTINE hydrounit_build ()

      USE precision
      USE spmd_task
      USE mod_utils
      USE mod_block
      USE mod_pixel
      USE mod_grid
      USE mod_data_type
      USE mod_landbasin
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
      INTEGER :: nhru_glb

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land hydro units :'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN
         CALL allocate_block_data (ghydru, hydrdata)
         CALL hydro_data_read (DEF_dir_hydrodata, 'hunit', ghydru, hydrdata)
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
                  xblk = ghydru%xblk(xlist(ig))
                  yblk = ghydru%yblk(ylist(ig))
                  xloc = ghydru%xloc(xlist(ig))
                  yloc = ghydru%yloc(ylist(ig))

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

         allocate (unum_tmp (numbasin*DEF_max_hband))
         allocate (iunt_tmp (numbasin*DEF_max_hband))
         allocate (ltyp_tmp (numbasin*DEF_max_hband))
         allocate (istt_tmp (numbasin*DEF_max_hband))
         allocate (iend_tmp (numbasin*DEF_max_hband))

         numhru = 0

         DO iu = 1, numbasin
         
            npxl = landbasin(iu)%npxl 
            
            allocate (ltype (1:npxl))

#ifdef USEMPI
            allocate (xlist (1:npxl))
            allocate (ylist (1:npxl))
            allocate (ipt   (1:npxl))
            allocate (msk   (1:npxl))

            DO ipxl = 1, npxl
               ilon_g = ghydru%xgrd(landbasin(iu)%ilon(ipxl))
               ilat_g = ghydru%ygrd(landbasin(iu)%ilat(ipxl))
               xlist(ipxl) = ilon_g
               ylist(ipxl) = ilat_g

               xblk = ghydru%xblk(xlist(ipxl))
               yblk = ghydru%yblk(ylist(ipxl))
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
               ilon_g = ghydru%xgrd(landbasin(iu)%ilon(ipxl))
               ilat_g = ghydru%ygrd(landbasin(iu)%ilat(ipxl))
               xblk   = ghydru%xblk(ilon_g)
               yblk   = ghydru%yblk(ilat_g)
               xloc   = ghydru%xloc(ilon_g)
               yloc   = ghydru%yloc(ilat_g)

               ltype(ipxl) = hydrdata%blk(xblk,yblk)%val(xloc,yloc)
            ENDDO
#endif
               
            allocate (order (1:npxl))
            order = (/ (ipxl, ipxl = 1, npxl) /)

            CALL quicksort (npxl, ltype, order)
               
            landbasin(iu)%ilon(1:npxl) = landbasin(iu)%ilon(order)
            landbasin(iu)%ilat(1:npxl) = landbasin(iu)%ilat(order)
            
            DO ipxl = 1, npxl
               IF (ipxl == 1) THEN
                  numhru = numhru + 1 
                  unum_tmp (numhru) = landbasin(iu)%num
                  iunt_tmp (numhru) = iu
                  ltyp_tmp (numhru) = ltype(ipxl)
                  istt_tmp (numhru) = ipxl
               ELSEIF (ltype(ipxl) /= ltype(ipxl-1)) THEN
                  iend_tmp(numhru) = ipxl - 1

                  numhru = numhru + 1
                  unum_tmp (numhru) = landbasin(iu)%num
                  iunt_tmp (numhru) = iu
                  ltyp_tmp (numhru) = ltype(ipxl)
                  istt_tmp (numhru) = ipxl
               ENDIF
            ENDDO
            iend_tmp(numhru) = npxl
            
            deallocate (ltype)
            deallocate (order)

         ENDDO

         allocate (hydrounit%unum (numhru))
         allocate (hydrounit%iunt (numhru))
         allocate (hydrounit%ltyp (numhru))
         allocate (hydrounit%istt (numhru))
         allocate (hydrounit%iend (numhru))
         
         hydrounit%unum = unum_tmp(1:numhru)  
         hydrounit%iunt = iunt_tmp(1:numhru)  
         hydrounit%ltyp = ltyp_tmp (1:numhru)  
         hydrounit%istt = istt_tmp (1:numhru)
         hydrounit%iend = iend_tmp (1:numhru)

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

      hydrounit%nset = numhru
      CALL hydrounit%set_vecgs 
         
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         CALL mpi_reduce (numhru, nhru_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', nhru_glb, ' hydro units on worker.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numbasin, ' hydro units.'
#endif

   END SUBROUTINE hydrounit_build
#endif

END MODULE mod_hydrounit
