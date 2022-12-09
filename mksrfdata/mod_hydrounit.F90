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

#if (defined GRIDBASED || defined SinglePoint || defined UNSTRUCTURED)
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

         allocate (hydrounit%bindex (numhru))
         allocate (hydrounit%ibasin (numhru))
         allocate (hydrounit%ipxstt (numhru))
         allocate (hydrounit%ipxend (numhru))
         allocate (hydrounit%ltyp   (numhru))

         DO ihru = 1, numhru
            hydrounit%bindex(ihru) = landbasin(ihru)%indx
            hydrounit%ibasin(ihru) = ihru
            hydrounit%ipxstt(ihru) = 1 
            hydrounit%ipxend(ihru) = landbasin(ihru)%npxl 
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
      write(*,'(A,I12,A)') 'Total: ', numhru, ' hydro units.'
#endif

   END SUBROUTINE hydrounit_build
#endif


#ifdef CATCHMENT 
   ! -------------------------------
   SUBROUTINE hydrounit_build ()

      USE precision
      USE spmd_task
      USE ncio_serial
      USE mod_utils
      USE mod_block
      USE mod_pixel
      USE mod_grid
      USE mod_data_type
      USE mod_landbasin
      USE mod_catchment_data
      USE mod_namelist

      IMPLICIT NONE

      ! Local Variables
      CHARACTER(len=256) :: file_drainage_network
      INTEGER , allocatable :: varsize (:)
      INTEGER :: maxnumhru
      TYPE (block_data_int32_2d) :: hydrdata
      INTEGER :: nreq, rmesg(2), smesg(2), isrc, idest, iproc
      INTEGER :: ilon_g, ilat_g, iloc
      INTEGER :: ig, xblk, yblk, xloc, yloc
      INTEGER :: npxl, ipxl
      INTEGER :: iu
      INTEGER, allocatable :: xlist(:), ylist(:), sbuf(:), rbuf(:)
      INTEGER, allocatable :: ltype(:), ipt(:), order(:)
      INTEGER, allocatable :: bindex_tmp(:), ltyp_tmp(:), ipxstt_tmp(:), ipxend_tmp(:), ibasin_tmp(:)
      LOGICAL, allocatable :: msk(:)
      LOGICAL, allocatable :: worker_done(:)
      INTEGER :: nhru_glb

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land hydro units :'
         
         IF (catchment_data_in_one_file) THEN
            file_drainage_network = trim(DEF_path_catchment_data) 
         ELSE
            file_drainage_network = trim(DEF_path_catchment_data) // '/' // 'drainage_network.nc'
         ENDIF

         CALL ncio_inquire_varsize (file_drainage_network, 'hydrounit_index', varsize)
         maxnumhru = varsize(1) 
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
      CALL mpi_bcast (maxnumhru, 1, mpi_integer, p_root, p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN
         CALL allocate_block_data (ghydru, hydrdata)
      ENDIF

      CALL catchment_data_read (DEF_path_catchment_data, 'ihydrounit2d', ghydru, hydrdata, &
         catchment_data_in_one_file)

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

         allocate (bindex_tmp (numbasin*maxnumhru))
         allocate (ibasin_tmp (numbasin*maxnumhru))
         allocate (ltyp_tmp   (numbasin*maxnumhru))
         allocate (ipxstt_tmp (numbasin*maxnumhru))
         allocate (ipxend_tmp (numbasin*maxnumhru))

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
                  bindex_tmp (numhru) = landbasin(iu)%indx
                  ibasin_tmp (numhru) = iu
                  ltyp_tmp (numhru) = ltype(ipxl)
                  ipxstt_tmp (numhru) = ipxl
               ELSEIF (ltype(ipxl) /= ltype(ipxl-1)) THEN
                  ipxend_tmp(numhru) = ipxl - 1

                  numhru = numhru + 1
                  bindex_tmp (numhru) = landbasin(iu)%indx
                  ibasin_tmp (numhru) = iu
                  ltyp_tmp (numhru) = ltype(ipxl)
                  ipxstt_tmp (numhru) = ipxl
               ENDIF
            ENDDO
            ipxend_tmp(numhru) = npxl
            
            deallocate (ltype)
            deallocate (order)

         ENDDO

         allocate (hydrounit%bindex (numhru))
         allocate (hydrounit%ibasin (numhru))
         allocate (hydrounit%ltyp (numhru))
         allocate (hydrounit%ipxstt (numhru))
         allocate (hydrounit%ipxend (numhru))
         
         hydrounit%bindex = bindex_tmp(1:numhru)  
         hydrounit%ibasin = ibasin_tmp(1:numhru)  
         hydrounit%ltyp = ltyp_tmp (1:numhru)  
         hydrounit%ipxstt = ipxstt_tmp (1:numhru)
         hydrounit%ipxend = ipxend_tmp (1:numhru)

         deallocate (ltyp_tmp)
         deallocate (ipxstt_tmp)
         deallocate (ipxend_tmp)
         deallocate (ibasin_tmp)
         deallocate (bindex_tmp)

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
      write(*,'(A,I12,A)') 'Total: ', numhru, ' hydro units.'
#endif

   END SUBROUTINE hydrounit_build
#endif

END MODULE mod_hydrounit
