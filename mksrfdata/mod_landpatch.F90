#include <define.h>

MODULE mod_landpatch

   USE precision
   USE mod_grid
   USE mod_pixelset
   USE GlobalVars
   USE LC_Const
#ifdef SinglePoint
   USE mod_single_srfdata
#endif
   IMPLICIT NONE

   ! ---- Instance ----
   INTEGER :: numpatch
   TYPE(grid_type)     :: gpatch
   TYPE(pixelset_type) :: landpatch

#if (defined CROP) 
   TYPE(grid_type) :: gcrop
   REAL(r8), allocatable :: pctcrop   (:)
   INTEGER,  allocatable :: cropclass (:)
#endif

CONTAINS
   
   ! -------------------------------
   SUBROUTINE landpatch_build ()

      USE precision
      USE spmd_task
      USE mod_utils
      USE mod_block
      USE mod_pixel
      USE mod_grid
      USE mod_data_type
      USE mod_landbasin
      USE mod_hydrounit
      USE mod_namelist
      USE ncio_block
      USE mod_modis_data
#if (defined CROP) 
      USE mod_pixelsetshadow
#endif

      IMPLICIT NONE

      ! Local Variables
      CHARACTER(len=256) :: file_patch
      TYPE (block_data_int32_2d) :: patchdata
      INTEGER :: nreq, rmesg(2), smesg(2), isrc, idest, iproc
      INTEGER :: ilon_g, ilat_g, iloc
      INTEGER :: ig, xblk, yblk, xloc, yloc
      INTEGER :: npxl, ipxl
      INTEGER :: iu, ihru, ipxstt, ipxend, ipatch
      INTEGER, allocatable :: xlist(:), ylist(:), sbuf(:), rbuf(:)
      INTEGER, allocatable :: ltype(:), ipt(:), order(:)
      INTEGER, allocatable :: bindex_tmp(:), ltyp_tmp(:), ipxstt_tmp(:), ipxend_tmp(:), ibasin_tmp(:)
      LOGICAL, allocatable :: msk(:)
      LOGICAL, allocatable :: worker_done(:)
      INTEGER :: npatch_glb
#if (defined CROP) 
      TYPE(block_data_real8_3d) :: cropdata
      INTEGER :: cropfilter(1)
#endif
#ifdef USE_DOMINANT_PATCHTYPE
      INTEGER :: dominant_type
      INTEGER, allocatable :: npxl_ltype (:)
#endif

      INTEGER :: iblk, jblk

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land patches :'
      ENDIF

#if (defined SinglePoint && defined PFT_CLASSIFICATION && defined CROP)
      IF ((SITE_landtype == 12) .and. (USE_SITE_pctcrop)) THEN

         numpatch = count(SITE_pctcrop > 0.)      

         allocate (pctcrop  (numpatch))
         allocate (cropclass(numpatch))
         cropclass = pack(SITE_croptyp, SITE_pctcrop > 0.)
         pctcrop   = pack(SITE_pctcrop, SITE_pctcrop > 0.)

         pctcrop = pctcrop / sum(pctcrop)

         allocate (landpatch%bindex (numpatch))
         allocate (landpatch%ibasin (numpatch))
         allocate (landpatch%ipxstt (numpatch))
         allocate (landpatch%ipxend (numpatch))
         allocate (landpatch%ltyp   (numpatch))

         landpatch%bindex(:) = 1
         landpatch%ibasin(:) = 1
         landpatch%ipxstt(:) = 1 
         landpatch%ipxend(:) = 1
         landpatch%ltyp  (:) = 12

         landpatch%nset = numpatch
         CALL landpatch%set_vecgs 

         RETURN
      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifndef SinglePoint
      IF (p_is_io) THEN
         CALL allocate_block_data (gpatch, patchdata)
      
         file_patch = trim(DEF_dir_rawdata) // '/landtype_update.nc'
         CALL ncio_read_block (file_patch, 'landtype', gpatch, patchdata)

      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN 

         allocate (worker_done (0:p_np_worker-1))
         worker_done(:) = .false.

         DO while (.not. all(worker_done))

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
                  xblk = gpatch%xblk(xlist(ig))
                  yblk = gpatch%yblk(ylist(ig))
                  xloc = gpatch%xloc(xlist(ig))
                  yloc = gpatch%yloc(ylist(ig))

                  sbuf(ig) = patchdata%blk(xblk,yblk)%val(xloc,yloc)
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
         ENDDO

         deallocate (worker_done)
      ENDIF 
#endif
#endif

      IF (p_is_worker) THEN

         IF (numhru > 0) THEN
            allocate (bindex_tmp (numhru*N_land_classification))
            allocate (ibasin_tmp (numhru*N_land_classification))
            allocate (ltyp_tmp   (numhru*N_land_classification))
            allocate (ipxstt_tmp (numhru*N_land_classification))
            allocate (ipxend_tmp (numhru*N_land_classification))
         ENDIF

         numpatch = 0

         DO ihru = 1, numhru
         
            iu     = hydrounit%ibasin(ihru)
            ipxstt = hydrounit%ipxstt(ihru)
            ipxend = hydrounit%ipxend(ihru)
            npxl   = ipxend - ipxstt + 1 
            
            allocate (ltype (ipxstt:ipxend))

#ifndef SinglePoint
#ifdef USEMPI
            allocate (xlist (ipxstt:ipxend))
            allocate (ylist (ipxstt:ipxend))
            allocate (ipt   (ipxstt:ipxend))
            allocate (msk   (ipxstt:ipxend))

            DO ipxl = ipxstt, ipxend
               ilon_g = gpatch%xgrd(landbasin(iu)%ilon(ipxl))
               ilat_g = gpatch%ygrd(landbasin(iu)%ilat(ipxl))
               xlist(ipxl) = ilon_g
               ylist(ipxl) = ilat_g

               xblk = gpatch%xblk(xlist(ipxl))
               yblk = gpatch%yblk(ylist(ipxl))
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
            DO ipxl = ipxstt, ipxend
               ilon_g = gpatch%xgrd(landbasin(iu)%ilon(ipxl))
               ilat_g = gpatch%ygrd(landbasin(iu)%ilat(ipxl))
               xblk   = gpatch%xblk(ilon_g)
               yblk   = gpatch%yblk(ilat_g)
               xloc   = gpatch%xloc(ilon_g)
               yloc   = gpatch%yloc(ilat_g)

               ltype(ipxl) = patchdata%blk(xblk,yblk)%val(xloc,yloc)
            ENDDO
#endif
#else
            ltype(:) = SITE_landtype
#endif

#ifdef CATCHMENT
            IF (hydrounit%ltyp(ihru) == 0) THEN
#if (defined IGBP_CLASSIFICATION || defined PFT_CLASSIFICATION || defined PC_CLASSIFICATION) 
               ltype(ipxstt:ipxend) = 17
#elif (defined USGS_CLASSIFICATION)
               ltype(ipxstt:ipxend) = 16
#endif
            ENDIF
#endif

            allocate (order (ipxstt:ipxend))
            order = (/ (ipxl, ipxl = ipxstt, ipxend) /)

#ifdef PFT_CLASSIFICATION
            ! For classification of plant function types, merge all land types with soil ground 
            DO ipxl = ipxstt, ipxend
               IF (ltype(ipxl) > 0) THEN
                  IF (patchtypes(ltype(ipxl)) == 0) THEN
#if (defined CROP) 
                     !12  Croplands
                     !14  Cropland/Natural Vegetation Mosaics  ?
                     IF (ltype(ipxl) /= 12) THEN
                        ltype(ipxl) = 1
                     ENDIF
#else
                     ltype(ipxl) = 1
#endif
                  ENDIF
               ENDIF
            ENDDO
#endif
            
            CALL quicksort (npxl, ltype, order)
               
            landbasin(iu)%ilon(ipxstt:ipxend) = landbasin(iu)%ilon(order)
            landbasin(iu)%ilat(ipxstt:ipxend) = landbasin(iu)%ilat(order)

#ifdef USE_DOMINANT_PATCHTYPE
            allocate (npxl_ltype (0:maxval(ltype)))
            npxl_ltype(:) = 0
            DO ipxl = ipxstt, ipxend
               npxl_ltype(ltype(ipxl)) = npxl_ltype(ltype(ipxl)) + 1
            ENDDO 

            IF (any(ltype > 0)) THEN
               iloc = findloc(ltype > 0, .true., dim=1) + ipxstt - 1
               dominant_type = maxloc(npxl_ltype(1:), dim=1)
               ltype(iloc:ipxend) = dominant_type
            ENDIF

            deallocate(npxl_ltype)
#endif
            
            DO ipxl = ipxstt, ipxend
               IF (ipxl == ipxstt) THEN
                  numpatch = numpatch + 1 
                  ibasin_tmp(numpatch) = iu
                  bindex_tmp(numpatch) = landbasin(iu)%indx
                  ltyp_tmp  (numpatch) = ltype(ipxl)
                  ipxstt_tmp(numpatch) = ipxl
               ELSEIF (ltype(ipxl) /= ltype(ipxl-1)) THEN
                  ipxend_tmp(numpatch) = ipxl - 1

                  numpatch = numpatch + 1
                  ibasin_tmp(numpatch) = iu
                  bindex_tmp(numpatch) = landbasin(iu)%indx
                  ltyp_tmp  (numpatch) = ltype(ipxl)
                  ipxstt_tmp(numpatch) = ipxl
               ENDIF
            ENDDO
            ipxend_tmp(numpatch) = ipxend

            deallocate (ltype)
            deallocate (order)

         ENDDO
         
         IF (numpatch > 0) THEN
            allocate (landpatch%ibasin (numpatch))
            allocate (landpatch%bindex (numpatch))
            allocate (landpatch%ltyp   (numpatch))
            allocate (landpatch%ipxstt (numpatch))
            allocate (landpatch%ipxend (numpatch))

            landpatch%ibasin = ibasin_tmp(1:numpatch)  
            landpatch%bindex = bindex_tmp(1:numpatch)  
            landpatch%ltyp   = ltyp_tmp  (1:numpatch)  
            landpatch%ipxstt = ipxstt_tmp(1:numpatch)
            landpatch%ipxend = ipxend_tmp(1:numpatch)
         ENDIF

         IF (numhru > 0) THEN
            deallocate (ltyp_tmp)
            deallocate (ipxstt_tmp)
            deallocate (ipxend_tmp)
            deallocate (ibasin_tmp)
            deallocate (bindex_tmp)
         ENDIF

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            smesg = (/p_iam_glb, -1/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
         ENDDO
#endif

      ENDIF

      landpatch%nset = numpatch

      CALL landpatch%set_vecgs

#ifdef LANDONLY
      IF ((p_is_worker) .and. (numpatch > 0)) THEN
         allocate(msk(numpatch))
         msk = (landpatch%ltyp /= 0)
      ELSE
         allocate(msk(1))
         msk = .true.
      ENDIF
         
      CALL landpatch%pset_pack (msk, numpatch)

      deallocate(msk)
#endif 

#if (defined CROP) 
      IF (p_is_io) THEN
         file_patch = trim(DEF_dir_rawdata) // '/global_0.5x0.5.MOD2005_V4.5_CFT_mergetoclmpft.nc'
         CALL allocate_block_data (gcrop, cropdata, N_CFT)
         CALL ncio_read_block (file_patch, 'PCT_CFT', gcrop, N_CFT, cropdata)
      ENDIF

      cropfilter = (/ 12 /)

      CALL pixelsetshadow_build (landpatch, gcrop, cropdata, N_CFT, 1, cropfilter, &
         pctcrop, cropclass)

      numpatch = landpatch%nset
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         CALL mpi_reduce (numpatch, npatch_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npatch_glb, ' patches on worker.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numpatch, ' patches.'
#endif

   END SUBROUTINE landpatch_build

END MODULE mod_landpatch
