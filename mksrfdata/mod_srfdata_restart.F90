#include <define.h>  

MODULE mod_srfdata_restart

   IMPLICIT NONE

   INTEGER, parameter, PRIVATE :: rcompress = 1

   ! ----- subroutines -----
   PUBLIC :: mesh_save_to_file
   PUBLIC :: mesh_load_from_file
   
   PUBLIC :: pixelset_save_to_file
   PUBLIC :: pixelset_load_from_file
   
CONTAINS
   
   ! -----------------------
   SUBROUTINE mesh_save_to_file (lc_year, dir_landdata)

      USE spmd_task
      USE ncio_serial
      USE mod_mesh
      USE mod_block
      USE mod_utils
      IMPLICIT NONE

      INTEGER         , intent(in) :: lc_year
      CHARACTER(len=*), intent(in) :: dir_landdata

      ! Local variables
      CHARACTER(len=256) :: filename, fileblock, cyear
      INTEGER :: ie, je, nelm, elen, iblk, jblk, iworker, i
      INTEGER, allocatable :: nelm_worker(:), ndsp_worker(:)
      INTEGER, allocatable :: elmindx(:)
      INTEGER, allocatable :: npxlall(:)
      INTEGER, allocatable :: elmpixels(:,:,:)
      REAL(r8),allocatable :: lon(:), lat(:)

      ! add parameter input for time year
      write(cyear,'(i4.4)') lc_year
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif 
      IF (p_is_master) THEN
         write(*,*) 'Saving land elements ...'
         CALL system('mkdir -p ' // trim(dir_landdata) // '/mesh/' // trim(cyear))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif 

      filename = trim(dir_landdata) // '/mesh/' //trim(cyear) // '/mesh.nc'

      DO jblk = 1, gblock%nyblk
         DO iblk = 1, gblock%nxblk

#ifdef USEMPI 
            IF (p_is_worker) THEN
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
#endif
                  nelm = 0
                  elen = 0
                  DO ie = 1, numelm
                     IF ((mesh(ie)%xblk == iblk) .and. (mesh(ie)%yblk == jblk)) THEN
                        nelm = nelm + 1
                        elen = max(elen, mesh(ie)%npxl)
                     ENDIF 
                  ENDDO
                  
#ifdef USEMPI 
                  CALL mpi_allreduce (MPI_IN_PLACE, elen, 1, MPI_INTEGER, MPI_MAX, p_comm_group, p_err)
#endif

                  IF (nelm > 0) THEN

                     allocate (elmindx (nelm))
                     allocate (npxlall (nelm))
                     allocate (elmpixels (2,elen,nelm))

                     je = 0
                     DO ie = 1, numelm
                        IF ((mesh(ie)%xblk == iblk) .and. (mesh(ie)%yblk == jblk)) THEN
                           je = je + 1
                           elmindx(je) = mesh(ie)%indx
                           npxlall(je) = mesh(ie)%npxl

                           elmpixels(1,1:npxlall(je),je) = mesh(ie)%ilon
                           elmpixels(2,1:npxlall(je),je) = mesh(ie)%ilat
                        ENDIF
                     ENDDO
                  ENDIF

#ifdef USEMPI 
                  CALL mpi_gather (nelm, 1, MPI_INTEGER, &
                     MPI_INULL_P, 1, MPI_INTEGER, p_root, p_comm_group, p_err)
                  
                  CALL mpi_gatherv (elmindx, nelm, MPI_INTEGER, &
                     MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER, & ! insignificant on workers
                     p_root, p_comm_group, p_err)
         
                  CALL mpi_gatherv (npxlall, nelm, MPI_INTEGER, &
                     MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER, & ! insignificant on workers
                     p_root, p_comm_group, p_err)

                  DO ie = 1, nelm
                     CALL mpi_send (elmpixels(:,:,ie), 2*elen, MPI_INTEGER, &
                        p_root, mpi_tag_data, p_comm_group, p_err) 
                  ENDDO
               ENDIF
            ENDIF
#endif

#ifdef USEMPI 
            IF (p_is_io) THEN
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                 
                  elen = 0
                  CALL mpi_allreduce (MPI_IN_PLACE, elen, 1, MPI_INTEGER, MPI_MAX, p_comm_group, p_err)

                  allocate (nelm_worker (0:p_np_group-1))
                  nelm_worker(0) = 0                  
                  CALL mpi_gather (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     nelm_worker, 1, MPI_INTEGER, p_root, p_comm_group, p_err)
                  
                  nelm = sum(nelm_worker)

                  allocate (ndsp_worker(0:p_np_group-1))
                  ndsp_worker(0) = 0
                  DO iworker = 1, p_np_group-1
                     ndsp_worker(iworker) = ndsp_worker(iworker-1) + nelm_worker(iworker-1)
                  ENDDO

                  allocate (elmindx (nelm))
                  CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     elmindx, nelm_worker(0:), ndsp_worker(0:), MPI_INTEGER, &
                     p_root, p_comm_group, p_err)
         
                  allocate (npxlall (nelm))
                  CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     npxlall, nelm_worker(0:), ndsp_worker(0:), MPI_INTEGER, &
                     p_root, p_comm_group, p_err)
         
                  allocate (elmpixels (2, elen, nelm))
                  DO iworker = 1, p_np_group-1 
                     DO ie = ndsp_worker(iworker)+1, ndsp_worker(iworker)+nelm_worker(iworker)
                        CALL mpi_recv (elmpixels(:,:,ie), 2*elen, MPI_INTEGER, &
                           iworker, mpi_tag_data, p_comm_group, p_stat, p_err)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
#endif

            IF (p_is_io) THEN
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (nelm > 0) THEN
                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     CALL ncio_create_file (fileblock)

                     CALL ncio_define_dimension (fileblock, 'element',nelm)
                     CALL ncio_define_dimension (fileblock, 'np_max', elen)
                     CALL ncio_define_dimension (fileblock, 'ncoor',  2   )

                     CALL ncio_write_serial (fileblock, 'elmindex', elmindx, 'element')
                     CALL ncio_write_serial (fileblock, 'npxl',     npxlall, 'element')
                     CALL ncio_write_serial (fileblock, 'pixel',    elmpixels, &
                        'ncoor', 'np_max', 'element', compress = 1)
                  ENDIF
               ENDIF
            ENDIF
                     
            IF (allocated (elmindx))   deallocate(elmindx)
            IF (allocated (npxlall))   deallocate(npxlall)
            IF (allocated (elmpixels)) deallocate(elmpixels)
            
            IF (allocated (nelm_worker)) deallocate(nelm_worker)
            IF (allocated (ndsp_worker)) deallocate(ndsp_worker)

#ifdef USEMPI
            CALL mpi_barrier (p_comm_group, p_err)
#endif 
         ENDDO
      ENDDO

      IF (p_is_master) THEN

         CALL ncio_create_file (filename)
        
         CALL ncio_define_dimension (filename, 'xblk', gblock%nxblk)
         CALL ncio_define_dimension (filename, 'yblk', gblock%nyblk)
         CALL ncio_write_serial (filename, 'nelm_blk', nelm_blk, 'xblk', 'yblk')

         CALL ncio_define_dimension (filename, 'longitude', gridmesh%nlon)
         CALL ncio_define_dimension (filename, 'latitude' , gridmesh%nlat)

         allocate (lon (gridmesh%nlon))
         allocate (lat (gridmesh%nlat))

         DO i = 1, gridmesh%nlon
            lon(i) = (gridmesh%lon_w(i) + gridmesh%lon_e(i)) * 0.5
            IF (gridmesh%lon_w(i) > gridmesh%lon_e(i)) THEN
               lon(i) = lon(i) + 180.0
               CALL normalize_longitude (lon(i))
            ENDIF
         ENDDO
         CALL ncio_write_serial (filename, 'longitude', lon, 'longitude')

         DO i = 1, gridmesh%nlat
            lat(i) = (gridmesh%lat_s(i) + gridmesh%lat_n(i)) * 0.5
         ENDDO
         CALL ncio_write_serial (filename, 'latitude', lat, 'latitude')

         deallocate (lon)
         deallocate (lat)
         
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif 

      IF (p_is_master) write(*,*) 'SAVE land elements done.'

   END SUBROUTINE mesh_save_to_file

   !------------------------------------
   SUBROUTINE mesh_load_from_file (lc_year, dir_landdata)

      USE spmd_task
      USE mod_namelist
      USE mod_block
      USE ncio_serial
      USE mod_mesh
      IMPLICIT NONE

      INTEGER         , intent(in) :: lc_year
      CHARACTER(len=*), intent(in) :: dir_landdata

      ! Local variables
      CHARACTER(len=256) :: filename, fileblock, cyear
      INTEGER :: iblkme, iblk, jblk, ie, nelm, ndsp
      INTEGER, allocatable :: elmindx(:), npxl(:), pixels(:,:,:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif 
       
      IF (p_is_master) THEN
         write(*,*) 'Loading land elements ...'
      ENDIF
         
      ! add parameter input for time year
      write(cyear,'(i4.4)') lc_year
      filename = trim(dir_landdata) // '/mesh/' // trim(cyear) // '/mesh.nc'
      CALL ncio_read_bcast_serial (filename, 'nelm_blk', nelm_blk)

      IF (p_is_io) THEN

         numelm = sum(nelm_blk, mask = gblock%pio == p_iam_glb)

         print*, numelm
         IF (numelm > 0) THEN

            allocate (mesh (numelm))

            ndsp = 0
            DO iblkme = 1, gblock%nblkme 
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               nelm = nelm_blk(iblk,jblk)

               IF (nelm > 0) THEN

                  CALL get_filename_block (filename, iblk, jblk, fileblock)
                  CALL ncio_read_serial (fileblock, 'elmindex', elmindx)
                  CALL ncio_read_serial (fileblock, 'npxl',  npxl  )
                  CALL ncio_read_serial (fileblock, 'pixel', pixels)

                  DO ie = 1, nelm
                     mesh(ie+ndsp)%indx = elmindx(ie)
                     mesh(ie+ndsp)%npxl = npxl(ie)
                     mesh(ie+ndsp)%xblk = iblk
                     mesh(ie+ndsp)%yblk = jblk

                     allocate (mesh(ie+ndsp)%ilon (npxl(ie)))
                     allocate (mesh(ie+ndsp)%ilat (npxl(ie)))

                     mesh(ie+ndsp)%ilon = pixels(1,1:npxl(ie),ie)
                     mesh(ie+ndsp)%ilat = pixels(2,1:npxl(ie),ie)
                  ENDDO

                  ndsp = ndsp + nelm
               ENDIF
            ENDDO
         ENDIF
         
         IF (allocated(elmindx)) deallocate(elmindx)
         IF (allocated(npxl  ))  deallocate(npxl   )
         IF (allocated(pixels))  deallocate(pixels )
            
      ENDIF

#ifdef CLMDEBUG 
      IF (p_is_io) write(*,'(I10,A,I4)') numelm, ' elements on group ', p_iam_io
#endif

#ifdef USEMPI
      CALL scatter_mesh_from_io_to_worker
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      
      IF (p_is_master) THEN
         write(*,*) 'Loading land elements done.'
      ENDIF
      
   END SUBROUTINE mesh_load_from_file 

   !------------------------------------------------
   SUBROUTINE pixelset_save_to_file (lc_year, dir_landdata, psetname, pixelset)

      USE spmd_task
      USE mod_block
      USE ncio_vector
      USE mod_pixelset
      IMPLICIT NONE

      INTEGER         ,    intent(in) :: lc_year
      CHARACTER(len=*),    intent(in) :: dir_landdata
      CHARACTER(len=*),    intent(in) :: psetname
      TYPE(pixelset_type), intent(in) :: pixelset

      ! Local variables
      CHARACTER(len=256)   :: filename, cyear

      write(cyear,'(i4.4)') lc_year
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,*) 'Saving Pixel Sets ' // trim(psetname) // ' ...'
         CALL system('mkdir -p ' // trim(dir_landdata) // '/' // trim(psetname) // '/' // trim(cyear))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      filename = trim(dir_landdata) // '/' // trim(psetname) // '/' // trim(cyear) // '/' // trim(psetname) // '.nc'

      CALL ncio_create_file_vector (filename, pixelset)
      CALL ncio_define_dimension_vector (filename, pixelset, trim(psetname))

      CALL ncio_write_vector (filename, 'eindex', trim(psetname), pixelset, pixelset%eindex, rcompress)
      CALL ncio_write_vector (filename, 'ipxstt', trim(psetname), pixelset, pixelset%ipxstt, rcompress)
      CALL ncio_write_vector (filename, 'ipxend', trim(psetname), pixelset, pixelset%ipxend, rcompress)
      CALL ncio_write_vector (filename, 'settyp', trim(psetname), pixelset, pixelset%settyp, rcompress)
     
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) write(*,*) 'SAVE Pixel Sets ' // trim(psetname) // ' done.'

   END SUBROUTINE pixelset_save_to_file


   !---------------------------
   SUBROUTINE pixelset_load_from_file (lc_year, dir_landdata, psetname, pixelset, numset)

      USE spmd_task
      USE mod_block
      USE ncio_serial
      USE ncio_vector
      USE mod_mesh
      USE mod_pixelset
      IMPLICIT NONE

      INTEGER         ,    intent(in) :: lc_year
      CHARACTER(len=*),    intent(in) :: dir_landdata
      CHARACTER(len=*),    intent(in) :: psetname
      TYPE(pixelset_type), intent(inout) :: pixelset
      INTEGER, intent(out) :: numset

      ! Local variables
      CHARACTER(len=256) :: filename, fileblock, cyear
      INTEGER :: iset, nset, ndsp, iblkme, iblk, jblk, ie, je, nave, nres, left, iproc
      INTEGER :: nsend, nrecv
      INTEGER, allocatable :: rbuff(:), iworker(:), sbuff(:)
      LOGICAL, allocatable :: msk(:)
      LOGICAL :: fexists
     
      write(cyear,'(i4.4)') lc_year
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
         
      IF (p_is_master) THEN
         write(*,*) 'Loading Pixel Sets ' // trim(psetname) // ' ...'
      ENDIF
      
      filename = trim(dir_landdata) // '/' // trim(psetname) // '/' // trim(cyear) // '/' // trim(psetname) // '.nc'

      IF (p_is_io) THEN

         pixelset%nset = 0
                  
         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
         
            CALL get_filename_block (filename, iblk, jblk, fileblock)

            inquire (file=trim(fileblock), exist=fexists)
            IF (fexists) THEN 
               CALL ncio_inquire_length (fileblock, 'eindex', nset)
               pixelset%nset = pixelset%nset + nset
            ENDIF

         ENDDO

         IF (pixelset%nset > 0) THEN

            allocate (pixelset%eindex (pixelset%nset))

            ndsp = 0
            DO iblkme = 1, gblock%nblkme 
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               CALL get_filename_block (filename, iblk, jblk, fileblock)
               inquire (file=trim(fileblock), exist=fexists)
               IF (fexists) THEN 

                  CALL ncio_read_serial (fileblock, 'eindex', rbuff) 

                  nset = size(rbuff)
                  pixelset%eindex(ndsp+1:ndsp+nset) = rbuff

                  ndsp = ndsp + nset
               ENDIF

            ENDDO
         ENDIF
      ENDIF
         

#ifdef USEMPI
      IF (p_is_io) THEN
         IF (pixelset%nset > 0) THEN
            allocate (iworker (pixelset%nset))
            allocate (msk     (pixelset%nset))

            ie = 1
            je = 1
            iblk = mesh(ie)%xblk
            jblk = mesh(ie)%yblk
            DO iset = 1, pixelset%nset
               DO WHILE (pixelset%eindex(iset) /= mesh(ie)%indx)
                  ie = ie + 1
                  je = je + 1
                  IF ((mesh(ie)%xblk /= iblk) .or. (mesh(ie)%yblk /= jblk)) THEN
                     je = 1
                     iblk = mesh(ie)%xblk
                     jblk = mesh(ie)%yblk
                  ENDIF
               ENDDO

               nave = nelm_blk(iblk,jblk) / (p_np_group-1)
               nres = mod(nelm_blk(iblk,jblk), p_np_group-1)
               left = (nave+1) * nres
               IF (je <= left) THEN
                  iworker(iset) = (je-1) / (nave+1) + 1
               ELSE
                  iworker(iset) = (je-left-1) / nave + 1 + nres
               ENDIF
            ENDDO

            DO iproc = 1, p_np_group-1
               msk = (iworker == iproc)
               nsend = count(msk)
               CALL mpi_send (nsend, 1, MPI_INTEGER, iproc, mpi_tag_size, p_comm_group, p_err) 

               IF (nsend > 0) THEN
                  allocate (sbuff(nsend))
                  sbuff = pack(pixelset%eindex, msk)
                  CALL mpi_send (sbuff, nsend, MPI_INTEGER, iproc, mpi_tag_data, p_comm_group, p_err) 
                  deallocate (sbuff)
               ENDIF
            ENDDO
         ELSE
            DO iproc = 1, p_np_group-1
               nsend = 0
               CALL mpi_send (nsend, 1, MPI_INTEGER, iproc, mpi_tag_size, p_comm_group, p_err) 
            ENDDO
         ENDIF

      ENDIF

      IF (p_is_worker) THEN
         
         CALL mpi_recv (nrecv, 1, MPI_INTEGER, p_root, mpi_tag_size, p_comm_group, p_stat, p_err)
         
         pixelset%nset = nrecv
         IF (nrecv > 0) THEN
            allocate (pixelset%eindex (nrecv))
            CALL mpi_recv (pixelset%eindex, nrecv, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_group, p_stat, p_err)
         ENDIF
      ENDIF
#endif

      
      CALL pixelset%set_vecgs

      CALL ncio_read_vector (filename, 'ipxstt', pixelset, pixelset%ipxstt)
      CALL ncio_read_vector (filename, 'ipxend', pixelset, pixelset%ipxend)
      CALL ncio_read_vector (filename, 'settyp', pixelset, pixelset%settyp)

      IF (p_is_worker) THEN
         IF (pixelset%nset > 0) THEN

            allocate (pixelset%ielm (pixelset%nset))
            ie = 1
            DO iset = 1, pixelset%nset
               DO WHILE (pixelset%eindex(iset) /= mesh(ie)%indx)
                  ie = ie + 1
               ENDDO
               pixelset%ielm(iset) = ie
            ENDDO

         ENDIF
      ENDIF

      numset = pixelset%nset

#ifdef CLMDEBUG 
      IF (p_is_io)  write(*,*) numset, trim(psetname), ' on group', p_iam_io
#endif

   END SUBROUTINE pixelset_load_from_file

END MODULE mod_srfdata_restart
