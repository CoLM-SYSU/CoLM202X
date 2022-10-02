#include <define.h>

MODULE mod_landpatch

   USE precision
   USE mod_grid
   USE mod_pixelset
   USE GlobalVars
   USE LC_Const
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
      INTEGER :: iu, ihru, istt, iend, ipatch
      INTEGER, allocatable :: xlist(:), ylist(:), sbuf(:), rbuf(:)
      INTEGER, allocatable :: ltype(:), ipt(:), order(:)
      INTEGER, allocatable :: unum_tmp(:), ltyp_tmp(:), istt_tmp(:), iend_tmp(:), iunt_tmp(:)
      LOGICAL, allocatable :: msk(:)
      LOGICAL, allocatable :: worker_done(:)
      INTEGER :: npatch_glb
#if (defined CROP) 
      TYPE(block_data_real8_3d) :: cropdata
      INTEGER :: cropfilter(1)
      INTEGER :: ic, jc 
#endif

      INTEGER :: iblk, jblk

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land patches :'
      ENDIF

#ifdef SinglePoint
   
#ifdef PFT_CLASSIFICATION
      IF (patchtypes(SITE_landtype) == 0) THEN
#ifdef CROP
         IF (SITE_landtype == 12) THEN
            numpatch = count(SITE_pctcrop > 0.)      
            allocate (pctcrop  (numpatch))
            allocate (cropclass(numpatch))
            jc = 0
            DO ic = 1, N_CFT
               IF (SITE_pctcrop(ic) > 0.) THEN
                  jc = jc + 1
                  cropclass(jc) = ic
                  pctcrop  (jc) = SITE_pctcrop(ic) / sum(SITE_pctcrop)
               ENDIF
            ENDDO
         ELSE
            numpatch = 1
         ENDIF
#else
         numpatch = 1
#endif
      ENDIF
#else
      numpatch = 1
#endif

      allocate (landpatch%unum (numpatch))
      allocate (landpatch%iunt (numpatch))
      allocate (landpatch%istt (numpatch))
      allocate (landpatch%iend (numpatch))
      allocate (landpatch%ltyp (numpatch))

      DO ipatch = 1, numpatch 
         landpatch%unum(ipatch) = landbasin(1)%num
         landpatch%iunt(ipatch) = 1
         landpatch%istt(ipatch) = 1 
         landpatch%iend(ipatch) = 1
#ifdef PFT_CLASSIFICATION
         IF (patchtypes(SITE_landtype) == 0) THEN
#ifdef CROP
            IF (SITE_landtype == 12) THEN
               landpatch%ltyp(ipatch) = SITE_landtype
            ELSE
               landpatch%ltyp(ipatch) = 1
            ENDIF
#else
            landpatch%ltyp(ipatch) = 1
#endif
         ELSE
            landpatch%ltyp(ipatch) = SITE_landtype
         ENDIF
#else
         landpatch%ltyp(ipatch) = SITE_landtype
#endif
      ENDDO
      
      landpatch%nset = 1
      CALL landpatch%set_vecgs 

      RETURN
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN
         CALL allocate_block_data (gpatch, patchdata)
      
#ifdef USGS_CLASSIFICATION
         file_patch = trim(DEF_dir_rawdata) // '/landtype_usgs_update.nc'
         CALL ncio_read_block (file_patch, 'landtype', gpatch, patchdata)
#endif

#ifdef IGBP_CLASSIFICATION
         file_patch = trim(DEF_dir_rawdata) // '/landtype_igbp_update.nc'
         CALL ncio_read_block (file_patch, 'landtype', gpatch, patchdata)
#endif

#ifdef PFT_CLASSIFICATION
         file_patch = trim(DEF_dir_rawdata) // '/landtype_igbp_update.nc'
         CALL ncio_read_block (file_patch, 'landtype', gpatch, patchdata)
#endif

#ifdef PC_CLASSIFICATION
         file_patch = trim(DEF_dir_rawdata) // '/landtype_igbp_update.nc'
         CALL ncio_read_block (file_patch, 'landtype', gpatch, patchdata)
#endif
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

      IF (p_is_worker) THEN

         allocate (unum_tmp (numhru*N_land_classification))
         allocate (iunt_tmp (numhru*N_land_classification))
         allocate (ltyp_tmp (numhru*N_land_classification))
         allocate (istt_tmp (numhru*N_land_classification))
         allocate (iend_tmp (numhru*N_land_classification))

         numpatch = 0

         DO ihru = 1, numhru
         
            iu   = hydrounit%iunt(ihru)
            istt = hydrounit%istt(ihru)
            iend = hydrounit%iend(ihru)
            npxl = iend - istt + 1 
            
            allocate (ltype (istt:iend))

#ifdef USEMPI
            allocate (xlist (istt:iend))
            allocate (ylist (istt:iend))
            allocate (ipt   (istt:iend))
            allocate (msk   (istt:iend))

            DO ipxl = istt, iend
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
            DO ipxl = istt, iend
               ilon_g = gpatch%xgrd(landbasin(iu)%ilon(ipxl))
               ilat_g = gpatch%ygrd(landbasin(iu)%ilat(ipxl))
               xblk   = gpatch%xblk(ilon_g)
               yblk   = gpatch%yblk(ilat_g)
               xloc   = gpatch%xloc(ilon_g)
               yloc   = gpatch%yloc(ilat_g)

               ltype(ipxl) = patchdata%blk(xblk,yblk)%val(xloc,yloc)
            ENDDO
#endif
               
            allocate (order (istt:iend))
            order = (/ (ipxl, ipxl = istt, iend) /)

#ifdef PFT_CLASSIFICATION
            ! For classification of plant function types, merge all land types with soil ground 
            DO ipxl = istt, iend
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
            
            CALL quicksort (npxl, ltype, order)
               
            landbasin(iu)%ilon(istt:iend) = landbasin(iu)%ilon(order)
            landbasin(iu)%ilat(istt:iend) = landbasin(iu)%ilat(order)
            
            DO ipxl = istt, iend
               IF (ipxl == istt) THEN
                  numpatch = numpatch + 1 
                  iunt_tmp(numpatch) = iu
                  unum_tmp(numpatch) = landbasin(iu)%num
                  ltyp_tmp(numpatch) = ltype(ipxl)
                  istt_tmp(numpatch) = ipxl
               ELSEIF (ltype(ipxl) /= ltype(ipxl-1)) THEN
                  iend_tmp(numpatch) = ipxl - 1

                  numpatch = numpatch + 1
                  iunt_tmp(numpatch) = iu
                  unum_tmp(numpatch) = landbasin(iu)%num
                  ltyp_tmp(numpatch) = ltype(ipxl)
                  istt_tmp(numpatch) = ipxl
               ENDIF
            ENDDO
            iend_tmp(numpatch) = iend
#endif

#if (defined USGS_CLASSIFICATION || defined IGBP_CLASSIFICATION || defined PC_CLASSIFICATION) 
            CALL quicksort (npxl, ltype, order)
               
            landbasin(iu)%ilon(istt:iend) = landbasin(iu)%ilon(order)
            landbasin(iu)%ilat(istt:iend) = landbasin(iu)%ilat(order)
            
            DO ipxl = istt, iend
               IF (ipxl == istt) THEN
                  numpatch = numpatch + 1 
                  iunt_tmp(numpatch) = iu
                  unum_tmp(numpatch) = landbasin(iu)%num
                  ltyp_tmp(numpatch) = ltype(ipxl)
                  istt_tmp(numpatch) = ipxl
               ELSEIF (ltype(ipxl) /= ltype(ipxl-1)) THEN
                  iend_tmp(numpatch) = ipxl - 1

                  numpatch = numpatch + 1
                  iunt_tmp(numpatch) = iu
                  unum_tmp(numpatch) = landbasin(iu)%num
                  ltyp_tmp(numpatch) = ltype(ipxl)
                  istt_tmp(numpatch) = ipxl
               ENDIF
            ENDDO
            iend_tmp(numpatch) = iend
#endif
            
            deallocate (ltype)
            deallocate (order)

         ENDDO
         
         allocate (landpatch%iunt (numpatch))
         allocate (landpatch%unum (numpatch))
         allocate (landpatch%ltyp (numpatch))
         allocate (landpatch%istt (numpatch))
         allocate (landpatch%iend (numpatch))

         landpatch%iunt = iunt_tmp(1:numpatch)  
         landpatch%unum = unum_tmp(1:numpatch)  
         landpatch%ltyp = ltyp_tmp(1:numpatch)  
         landpatch%istt = istt_tmp(1:numpatch)
         landpatch%iend = iend_tmp(1:numpatch)

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

      landpatch%nset = numpatch

      CALL landpatch%set_vecgs

#ifdef LANDONLY
      IF (p_is_worker) THEN
         allocate(msk(numpatch))
         msk = (landpatch%ltyp /= 0)
      ENDIF
         
      CALL landpatch%pset_pack (msk, numpatch)
#endif 

#if (defined CROP) 
      IF (p_is_io) THEN
!         file_patch = trim(DEF_dir_rawdata) // '/global_0.5x0.5.MOD2005_V4.5_CFT_mergetoclmpft.nc'
         file_patch = trim(DEF_dir_rawdata) // '/global_0.5x0.5.MOD2005_V4.5_CFT_lf-merged-20220930.nc'
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
