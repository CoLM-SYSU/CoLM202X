#include <define.h>  

MODULE mod_srfdata_restart

   IMPLICIT NONE

   INTEGER, parameter, PRIVATE :: rcompress = 1

   ! ----- subroutines -----
   PUBLIC :: landbasin_save_to_file
   PUBLIC :: landbasin_load_from_file
   
   PUBLIC :: pixelset_save_to_file
   PUBLIC :: pixelset_load_from_file
   
CONTAINS
   
   ! -----------------------
   SUBROUTINE landbasin_save_to_file (dir_landdata)

      USE spmd_task
      USE ncio_serial
      USE mod_landbasin
      USE mod_block
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: dir_landdata

      ! Local variables
      CHARACTER(len=256) :: filename, fileblock
      INTEGER :: iu, ju, nbasin, ulen, iblk, jblk, iworker
      INTEGER, allocatable :: nbasin_worker(:), ndsp_worker(:)
      INTEGER, allocatable :: basinnum(:)
      INTEGER, allocatable :: npxlall(:)
      INTEGER, allocatable :: basinpixels(:,:,:)
      INTEGER, allocatable :: rbuf(:,:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif 
      IF (p_is_master) THEN
         write(*,*) 'Saving land basins ...'
         CALL system('mkdir -p ' // trim(dir_landdata) // '/landbasin')
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif 

      filename = trim(dir_landdata) // '/landbasin/landbasin.nc'

      DO jblk = 1, gblock%nyblk
         DO iblk = 1, gblock%nxblk

#ifdef USEMPI 
            IF (p_is_worker) THEN
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
#endif
                  nbasin = 0
                  ulen  = 0
                  DO iu = 1, numbasin
                     IF ((landbasin(iu)%xblk == iblk) .and. (landbasin(iu)%yblk == jblk)) THEN
                        nbasin = nbasin + 1
                        ulen  = max(ulen, landbasin(iu)%npxl)
                     ENDIF 
                  ENDDO
                  
#ifdef USEMPI 
                  CALL mpi_allreduce (MPI_IN_PLACE, ulen, 1, MPI_INTEGER, MPI_MAX, p_comm_group, p_err)
#endif

                  IF (nbasin > 0) THEN

                     allocate (basinnum (nbasin))
                     allocate (npxlall (nbasin))
                     allocate (basinpixels (2,ulen,nbasin))

                     ju = 0
                     DO iu = 1, numbasin
                        IF ((landbasin(iu)%xblk == iblk) .and. (landbasin(iu)%yblk == jblk)) THEN
                           ju = ju + 1
                           basinnum(ju) = landbasin(iu)%num
                           npxlall(ju) = landbasin(iu)%npxl

                           basinpixels(1,1:npxlall(ju),ju) = landbasin(iu)%ilon
                           basinpixels(2,1:npxlall(ju),ju) = landbasin(iu)%ilat
                        ENDIF
                     ENDDO
                  ENDIF

#ifdef USEMPI 
                  CALL mpi_gather (nbasin, 1, MPI_INTEGER, &
                     nbasin_worker, 1, MPI_INTEGER, p_root, p_comm_group, p_err)
                  
                  CALL mpi_gatherv (basinnum, nbasin, MPI_INTEGER, &
                     basinnum, nbasin_worker, ndsp_worker, MPI_INTEGER, &
                     p_root, p_comm_group, p_err)
         
                  CALL mpi_gatherv (npxlall, nbasin, MPI_INTEGER, &
                     npxlall, nbasin_worker, ndsp_worker, MPI_INTEGER, &
                     p_root, p_comm_group, p_err)

                 ! CALL mpi_gatherv (basinpixels, nbasin*2*ulen, MPI_INTEGER, &
                 !    basinpixels, nbasin_worker, ndsp_worker, MPI_INTEGER, &
                 !    p_root, p_comm_group, p_err)
                  DO iu = 1, nbasin
                     CALL mpi_send (basinpixels(:,:,iu), 2*ulen, MPI_INTEGER, &
                        p_root, mpi_tag_data, p_comm_group, p_err) 
                  ENDDO
               ENDIF
            ENDIF
#endif

#ifdef USEMPI 
            IF (p_is_io) THEN
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                 
                  ulen = 0
                  CALL mpi_allreduce (MPI_IN_PLACE, ulen, 1, MPI_INTEGER, MPI_MAX, p_comm_group, p_err)

                  allocate (nbasin_worker (0:p_np_group-1))
                  nbasin_worker(0) = 0                  
                  CALL mpi_gather (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     nbasin_worker, 1, MPI_INTEGER, p_root, p_comm_group, p_err)
                  
                  nbasin = sum(nbasin_worker)

                  allocate (ndsp_worker(0:p_np_group-1))
                  ndsp_worker(0) = 0
                  DO iworker = 1, p_np_group-1
                     ndsp_worker(iworker) = ndsp_worker(iworker-1) + nbasin_worker(iworker-1)
                  ENDDO

                  allocate (basinnum (nbasin))
                  CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     basinnum, nbasin_worker(0:), ndsp_worker(0:), MPI_INTEGER, &
                     p_root, p_comm_group, p_err)
         
                  allocate (npxlall (nbasin))
                  CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     npxlall, nbasin_worker(0:), ndsp_worker(0:), MPI_INTEGER, &
                     p_root, p_comm_group, p_err)
         
                  allocate (basinpixels (2, ulen, nbasin))
                  ! CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                  !    basinpixels, nbasin_worker*2*ulen, ndsp_worker*2*ulen, MPI_INTEGER, &
                  !    p_root, p_comm_group, p_err)
                  DO iworker = 1, p_np_group-1 
                     DO iu = ndsp_worker(iworker)+1, ndsp_worker(iworker)+nbasin_worker(iworker)
                        CALL mpi_recv (basinpixels(:,:,iu), 2*ulen, MPI_INTEGER, &
                           iworker, mpi_tag_data, p_comm_group, p_stat, p_err)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
#endif

            IF (p_is_io) THEN
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (nbasin > 0) THEN
                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     CALL ncio_create_file (fileblock)

                     CALL ncio_define_dimension (fileblock, 'landbasin', nbasin)
                     CALL ncio_define_dimension (fileblock, 'np_max',   ulen )
                     CALL ncio_define_dimension (fileblock, 'ncoor',    2    )

                     CALL ncio_write_serial (fileblock, 'basinnum', basinnum,   'landbasin')
                     CALL ncio_write_serial (fileblock, 'npxl',    npxlall,   'landbasin')
                     CALL ncio_write_serial (fileblock, 'pixel',   basinpixels, &
                        'ncoor', 'np_max', 'landbasin', compress = 1)
                  ENDIF
               ENDIF
            ENDIF
                     
            IF (allocated (basinnum))     deallocate(basinnum)
            IF (allocated (npxlall))     deallocate(npxlall)
            IF (allocated (basinpixels))  deallocate(basinpixels)
            
            IF (allocated (nbasin_worker))  deallocate(nbasin_worker)
            IF (allocated (ndsp_worker ))  deallocate(ndsp_worker )

#ifdef USEMPI
            CALL mpi_barrier (p_comm_group, p_err)
#endif 
         ENDDO
      ENDDO

      IF (p_is_master) THEN
         CALL ncio_create_file (filename)
         CALL ncio_define_dimension (filename, 'xblk', gblock%nxblk)
         CALL ncio_define_dimension (filename, 'yblk', gblock%nyblk)
         CALL ncio_write_serial (filename, 'nbasin_blk', nbasin_blk, 'xblk', 'yblk')
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif 

      IF (p_is_master) write(*,*) 'SAVE land basins done.'

   END SUBROUTINE landbasin_save_to_file

   !------------------------------------
   SUBROUTINE landbasin_load_from_file (dir_landdata)

      USE spmd_task
      USE mod_namelist
      USE mod_block
      USE ncio_serial
      USE mod_landbasin
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: dir_landdata

      ! Local variables
      CHARACTER(len=256) :: filename, fileblock
      INTEGER :: iblkme, iblk, jblk, iu, nbasin, ndsp
      INTEGER, allocatable :: basinnum(:), npxl(:), pixels(:,:,:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif 
      
      IF (p_is_master) THEN
         write(*,*) 'Loading land basins ...'
      ENDIF
         
      filename = trim(dir_landdata) // '/landbasin/landbasin.nc'
      CALL ncio_read_bcast_serial (filename, 'nbasin_blk', nbasin_blk)

      IF (p_is_io) THEN

         numbasin = sum(nbasin_blk, mask = gblock%pio == p_iam_glb)

         IF (numbasin > 0) THEN

            allocate (landbasin (numbasin))

            ndsp = 0
            DO iblkme = 1, nblkme 
               iblk = xblkme(iblkme)
               jblk = yblkme(iblkme)

               nbasin = nbasin_blk(iblk,jblk)

               IF (nbasin > 0) THEN

                  CALL get_filename_block (filename, iblk, jblk, fileblock)
                  CALL ncio_read_serial (fileblock, 'basinnum', basinnum)
                  CALL ncio_read_serial (fileblock, 'npxl',    npxl   )
                  CALL ncio_read_serial (fileblock, 'pixel',   pixels )

                  DO iu = 1, nbasin
                     landbasin(iu+ndsp)%num  = basinnum(iu)
                     landbasin(iu+ndsp)%npxl = npxl(iu)
                     landbasin(iu+ndsp)%xblk = iblk
                     landbasin(iu+ndsp)%yblk = jblk

                     allocate (landbasin(iu+ndsp)%ilon (npxl(iu)))
                     allocate (landbasin(iu+ndsp)%ilat (npxl(iu)))

                     landbasin(iu+ndsp)%ilon = pixels(1,1:npxl(iu),iu)
                     landbasin(iu+ndsp)%ilat = pixels(2,1:npxl(iu),iu)
                  ENDDO

                  ndsp = ndsp + nbasin
               ENDIF
            ENDDO
         ENDIF
         
         IF (allocated(basinnum)) deallocate(basinnum)
         IF (allocated(npxl   )) deallocate(npxl   )
         IF (allocated(pixels )) deallocate(pixels )
            
      ENDIF

#ifdef CLMDEBUG 
      IF (p_is_io)     write(*,*) numbasin, ' basins on group', p_iam_io
#endif

#ifdef USEMPI
      CALL scatter_landbasin_from_io_to_worker
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      
      IF (p_is_master) THEN
         write(*,*) 'Loading landbasin done.'
      ENDIF
      
   END SUBROUTINE landbasin_load_from_file 

   !------------------------------------------------
   SUBROUTINE pixelset_save_to_file (dir_landdata, psetname, pixelset)

      USE spmd_task
      USE mod_block
      USE ncio_vector
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*),    intent(in) :: dir_landdata
      CHARACTER(len=*),    intent(in) :: psetname
      TYPE(pixelset_type), intent(in) :: pixelset

      ! Local variables
      CHARACTER(len=256)   :: filename

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,*) 'Saving Pixel Sets ' // trim(psetname) // ' ...'
         CALL system('mkdir -p ' // trim(dir_landdata) // '/' // trim(psetname))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      filename = trim(dir_landdata) // '/' // trim(psetname) // '/' // trim(psetname) // '.nc'

      CALL ncio_create_file_vector (filename, pixelset)
      CALL ncio_define_pixelset_dimension (filename, pixelset)

      CALL ncio_write_vector (filename, 'unum', 'vector', pixelset, pixelset%unum, rcompress)
      CALL ncio_write_vector (filename, 'istt', 'vector', pixelset, pixelset%istt, rcompress)
      CALL ncio_write_vector (filename, 'iend', 'vector', pixelset, pixelset%iend, rcompress)
      CALL ncio_write_vector (filename, 'ltyp', 'vector', pixelset, pixelset%ltyp, rcompress)
      
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) write(*,*) 'SAVE Pixel Sets ' // trim(psetname) // ' done.'

   END SUBROUTINE pixelset_save_to_file


   !---------------------------
   SUBROUTINE pixelset_load_from_file (dir_landdata, psetname, pixelset, numset)

      USE spmd_task
      USE mod_block
      USE ncio_serial
      USE ncio_vector
      USE mod_landbasin
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*),    intent(in) :: dir_landdata
      CHARACTER(len=*),    intent(in) :: psetname
      TYPE(pixelset_type), intent(inout) :: pixelset
      INTEGER, intent(out) :: numset

      ! Local variables
      CHARACTER(len=256) :: filename, fileblock
      INTEGER :: iset, nset, ndsp, iblkme, iblk, jblk, iu, ju, nave, nres, left, iproc
      INTEGER :: nsend, nrecv
      INTEGER, allocatable :: rbuff(:), iworker(:), sbuff(:)
      LOGICAL, allocatable :: msk(:)
      LOGICAL :: fexists
      
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
         
      IF (p_is_master) THEN
         write(*,*) 'Loading Pixel Sets ' // trim(psetname) // ' ...'
      ENDIF
      
      filename = trim(dir_landdata) // '/' // trim(psetname) // '/' // trim(psetname) // '.nc'

      IF (p_is_io) THEN

         pixelset%nset = 0
                  
         DO iblkme = 1, nblkme 
            iblk = xblkme(iblkme)
            jblk = yblkme(iblkme)
         
            CALL get_filename_block (filename, iblk, jblk, fileblock)

            inquire (file=trim(fileblock), exist=fexists)
            IF (fexists) THEN 
               CALL ncio_inquire_length (fileblock, 'unum', nset)
               pixelset%nset = pixelset%nset + nset
            ENDIF

         ENDDO

         IF (pixelset%nset > 0) THEN

            allocate (pixelset%unum (pixelset%nset))

            ndsp = 0
            DO iblkme = 1, nblkme 
               iblk = xblkme(iblkme)
               jblk = yblkme(iblkme)

               CALL get_filename_block (filename, iblk, jblk, fileblock)
               inquire (file=trim(fileblock), exist=fexists)
               IF (fexists) THEN 

                  CALL ncio_read_serial (fileblock, 'unum', rbuff) 

                  nset = size(rbuff)
                  pixelset%unum(ndsp+1:ndsp+nset) = rbuff

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

            iu = 1
            ju = 1
            iblk = landbasin(iu)%xblk
            jblk = landbasin(iu)%yblk
            DO iset = 1, pixelset%nset
               DO WHILE (pixelset%unum(iset) /= landbasin(iu)%num)
                  iu = iu + 1
                  ju = ju + 1
                  IF ((landbasin(iu)%xblk /= iblk) .or. (landbasin(iu)%yblk /= jblk)) THEN
                     ju = 1
                     iblk = landbasin(iu)%xblk
                     jblk = landbasin(iu)%yblk
                  ENDIF
               ENDDO

               nave = nbasin_blk(iblk,jblk) / (p_np_group-1)
               nres = mod(nbasin_blk(iblk,jblk), p_np_group-1)
               left = (nave+1) * nres
               IF (ju <= left) THEN
                  iworker(iset) = (ju-1) / (nave+1) + 1
               ELSE
                  iworker(iset) = (ju-left-1) / nave + 1 + nres
               ENDIF
            ENDDO

            DO iproc = 1, p_np_group-1
               msk = (iworker == iproc)
               nsend = count(msk)
               CALL mpi_send (nsend, 1, MPI_INTEGER, iproc, mpi_tag_size, p_comm_group, p_err) 

               IF (nsend > 0) THEN
                  allocate (sbuff(nsend))
                  sbuff = pack(pixelset%unum, msk)
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
            allocate (pixelset%unum (nrecv))
            CALL mpi_recv (pixelset%unum, nrecv, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_group, p_stat, p_err)
         ENDIF
      ENDIF
#endif

      
      CALL pixelset%set_vecgs

      CALL ncio_read_vector (filename, 'istt', pixelset, pixelset%istt)
      CALL ncio_read_vector (filename, 'iend', pixelset, pixelset%iend)
      CALL ncio_read_vector (filename, 'ltyp', pixelset, pixelset%ltyp)

      IF (p_is_worker) THEN
         IF (pixelset%nset > 0) THEN

            allocate (pixelset%iunt (pixelset%nset))
            iu = 1
            DO iset = 1, pixelset%nset
               DO WHILE (pixelset%unum(iset) /= landbasin(iu)%num)
                  iu = iu + 1
               ENDDO
               pixelset%iunt(iset) = iu
            ENDDO

         ENDIF
      ENDIF

      numset = pixelset%nset

#ifdef CLMDEBUG 
      IF (p_is_io)     write(*,*) numset, trim(psetname), ' on group', p_iam_io
#endif

   END SUBROUTINE pixelset_load_from_file

END MODULE mod_srfdata_restart
