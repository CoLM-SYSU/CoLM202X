#include <define.h>  

MODULE mod_srfdata_restart

   IMPLICIT NONE

   INTEGER, parameter, PRIVATE :: rcompress = 1

   ! ----- subroutines -----
   PUBLIC :: landunit_save_to_file
   PUBLIC :: landunit_load_from_file
   
   PUBLIC :: pixelset_save_to_file
   PUBLIC :: pixelset_load_from_file
   
CONTAINS
   
   ! -----------------------
   SUBROUTINE landunit_save_to_file (dir_landdata)

      USE spmd_task
      USE ncio_serial
      USE mod_landunit
      USE mod_block
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: dir_landdata

      ! Local variables
      CHARACTER(len=256) :: filename, fileblock
      INTEGER :: iu, ju, nunit, ulen, iblk, jblk, iworker
      INTEGER, allocatable :: nunit_worker(:), ndsp_worker(:)
      INTEGER, allocatable :: unitnum(:)
      INTEGER, allocatable :: npxlall(:)
      INTEGER, allocatable :: unitpixels(:,:,:)
      INTEGER, allocatable :: rbuf(:,:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif 

      IF (p_is_master) THEN
         write(*,*) 'Saving land units ...'
      ENDIF

      filename = trim(dir_landdata) // '/landunit.nc'

      DO jblk = 1, gblock%nyblk
         DO iblk = 1, gblock%nxblk

#ifdef USEMPI 
            IF (p_is_worker) THEN
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
#endif
                  nunit = 0
                  ulen  = 0
                  DO iu = 1, numunit
                     IF ((landunit(iu)%xblk == iblk) .and. (landunit(iu)%yblk == jblk)) THEN
                        nunit = nunit + 1
                        ulen  = max(ulen, landunit(iu)%npxl)
                     ENDIF 
                  ENDDO
                  
#ifdef USEMPI 
                  CALL mpi_allreduce (MPI_IN_PLACE, ulen, 1, MPI_INTEGER, MPI_MAX, p_comm_group, p_err)
#endif

                  IF (nunit > 0) THEN

                     allocate (unitnum (nunit))
                     allocate (npxlall (nunit))
                     allocate (unitpixels (2,ulen,nunit))

                     ju = 0
                     DO iu = 1, numunit
                        IF ((landunit(iu)%xblk == iblk) .and. (landunit(iu)%yblk == jblk)) THEN
                           ju = ju + 1
                           unitnum(ju) = landunit(iu)%num
                           npxlall(ju) = landunit(iu)%npxl

                           unitpixels(1,1:npxlall(ju),ju) = landunit(iu)%ilon
                           unitpixels(2,1:npxlall(ju),ju) = landunit(iu)%ilat
                        ENDIF
                     ENDDO
                  ENDIF

#ifdef USEMPI 
                  CALL mpi_gather (nunit, 1, MPI_INTEGER, &
                     nunit_worker, 1, MPI_INTEGER, p_root, p_comm_group, p_err)
                  
                  CALL mpi_gatherv (unitnum, nunit, MPI_INTEGER, &
                     unitnum, nunit_worker, ndsp_worker, MPI_INTEGER, &
                     p_root, p_comm_group, p_err)
         
                  CALL mpi_gatherv (npxlall, nunit, MPI_INTEGER, &
                     npxlall, nunit_worker, ndsp_worker, MPI_INTEGER, &
                     p_root, p_comm_group, p_err)

                 ! CALL mpi_gatherv (unitpixels, nunit*2*ulen, MPI_INTEGER, &
                 !    unitpixels, nunit_worker, ndsp_worker, MPI_INTEGER, &
                 !    p_root, p_comm_group, p_err)
                  DO iu = 1, nunit
                     CALL mpi_send (unitpixels(:,:,iu), 2*ulen, MPI_INTEGER, &
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

                  allocate (nunit_worker (0:p_np_group-1))
                  nunit_worker(0) = 0                  
                  CALL mpi_gather (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     nunit_worker, 1, MPI_INTEGER, p_root, p_comm_group, p_err)
                  
                  nunit = sum(nunit_worker)

                  allocate (ndsp_worker(0:p_np_group-1))
                  ndsp_worker(0) = 0
                  DO iworker = 1, p_np_group-1
                     ndsp_worker(iworker) = ndsp_worker(iworker-1) + nunit_worker(iworker-1)
                  ENDDO

                  allocate (unitnum (nunit))
                  CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     unitnum, nunit_worker(0:), ndsp_worker(0:), MPI_INTEGER, &
                     p_root, p_comm_group, p_err)
         
                  allocate (npxlall (nunit))
                  CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     npxlall, nunit_worker(0:), ndsp_worker(0:), MPI_INTEGER, &
                     p_root, p_comm_group, p_err)
         
                  allocate (unitpixels (2, ulen, nunit))
                  ! CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                  !    unitpixels, nunit_worker*2*ulen, ndsp_worker*2*ulen, MPI_INTEGER, &
                  !    p_root, p_comm_group, p_err)
                  DO iworker = 1, p_np_group-1 
                     DO iu = ndsp_worker(iworker)+1, ndsp_worker(iworker)+nunit_worker(iworker)
                        CALL mpi_recv (unitpixels(:,:,iu), 2*ulen, MPI_INTEGER, &
                           iworker, mpi_tag_data, p_comm_group, p_stat, p_err)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
#endif

            IF (p_is_io) THEN
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (nunit > 0) THEN
                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     CALL ncio_create_file (fileblock)

                     CALL ncio_define_dimension (fileblock, 'landunit', nunit)
                     CALL ncio_define_dimension (fileblock, 'np_max',   ulen )
                     CALL ncio_define_dimension (fileblock, 'ncoor',    2    )

                     CALL ncio_write_serial (fileblock, 'unitnum', unitnum,   'landunit')
                     CALL ncio_write_serial (fileblock, 'npxl',    npxlall,   'landunit')
                     CALL ncio_write_serial (fileblock, 'pixel',   unitpixels, &
                        'ncoor', 'np_max', 'landunit', compress = 1)
                  ENDIF
               ENDIF
            ENDIF
                     
            IF (allocated (unitnum))     deallocate(unitnum)
            IF (allocated (npxlall))     deallocate(npxlall)
            IF (allocated (unitpixels))  deallocate(unitpixels)
            
            IF (allocated (nunit_worker))  deallocate(nunit_worker)
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
         CALL ncio_write_serial (filename, 'nunits_blk', nunits_blk, 'xblk', 'yblk')
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif 

      IF (p_is_master) write(*,*) 'SAVE land units done.'

   END SUBROUTINE landunit_save_to_file

   !------------------------------------
   SUBROUTINE landunit_load_from_file (dir_landdata)

      USE spmd_task
      USE mod_namelist
      USE mod_block
      USE ncio_serial
      USE mod_landunit
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: dir_landdata

      ! Local variables
      CHARACTER(len=256) :: filename, fileblock
      INTEGER :: iblk, jblk, iu, nunit, ndsp
      INTEGER, allocatable :: unitnum(:), npxl(:), pixels(:,:,:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif 
      
      IF (p_is_master) THEN
         write(*,*) 'Loading land units ...'
      ENDIF
         
      filename = trim(dir_landdata) // '/landunit.nc'
      CALL ncio_read_bcast_serial (filename, 'nunits_blk', nunits_blk)

      IF (p_is_io) THEN

         numunit = sum(nunits_blk, mask = gblock%pio == p_iam_glb)

         IF (numunit > 0) THEN

            allocate (landunit (numunit))

            ndsp = 0
            DO jblk = 1, gblock%nyblk
               DO iblk = 1, gblock%nxblk
                  IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                     nunit = nunits_blk(iblk,jblk)

                     IF (nunit > 0) THEN
                     
                        CALL get_filename_block (filename, iblk, jblk, fileblock)
                        CALL ncio_read_serial (fileblock, 'unitnum', unitnum)
                        CALL ncio_read_serial (fileblock, 'npxl',    npxl   )
                        CALL ncio_read_serial (fileblock, 'pixel',   pixels )

                        DO iu = 1, nunit
                           landunit(iu+ndsp)%num  = unitnum(iu)
                           landunit(iu+ndsp)%npxl = npxl(iu)
                           landunit(iu+ndsp)%xblk = iblk
                           landunit(iu+ndsp)%yblk = jblk

                           allocate (landunit(iu+ndsp)%ilon (npxl(iu)))
                           allocate (landunit(iu+ndsp)%ilat (npxl(iu)))

                           landunit(iu+ndsp)%ilon = pixels(1,1:npxl(iu),iu)
                           landunit(iu+ndsp)%ilat = pixels(2,1:npxl(iu),iu)
                        ENDDO

                        ndsp = ndsp + nunit
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         
         IF (allocated(unitnum)) deallocate(unitnum)
         IF (allocated(npxl   )) deallocate(npxl   )
         IF (allocated(pixels )) deallocate(pixels )
            
      ENDIF

#ifdef CLMDEBUG 
      IF (p_is_io)     write(*,*) numunit, ' units on group', p_iam_io
#endif

#ifdef USEMPI
      CALL scatter_landunit_from_io_to_worker
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      
      IF (p_is_master) THEN
         write(*,*) 'Loading landunit done.'
      ENDIF
      
   END SUBROUTINE landunit_load_from_file 

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
      ENDIF

      filename = trim(dir_landdata) // '/' // trim(psetname) // '.nc'

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
      USE mod_landunit
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*),    intent(in) :: dir_landdata
      CHARACTER(len=*),    intent(in) :: psetname
      TYPE(pixelset_type), intent(inout) :: pixelset
      INTEGER, intent(out) :: numset

      ! Local variables
      CHARACTER(len=256) :: filename, fileblock
      INTEGER :: iset, nset, ndsp, iblk, jblk, iu, ju, nave, nres, left, iproc
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
      
      filename = trim(dir_landdata) // '/' // trim(psetname) // '.nc'

      IF (p_is_io) THEN

         pixelset%nset = 0
                  
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
         
                  CALL get_filename_block (filename, iblk, jblk, fileblock)
                  
                  inquire (file=trim(fileblock), exist=fexists)
                  IF (fexists) THEN 
                     CALL ncio_inquire_length (fileblock, 'unum', nset)
                     pixelset%nset = pixelset%nset + nset
                  ENDIF

               ENDIF
            ENDDO
         ENDDO

         IF (pixelset%nset > 0) THEN

            allocate (pixelset%unum (pixelset%nset))

            ndsp = 0
            DO jblk = 1, gblock%nyblk
               DO iblk = 1, gblock%nxblk
                  IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     inquire (file=trim(fileblock), exist=fexists)
                     IF (fexists) THEN 

                        CALL ncio_read_serial (fileblock, 'unum', rbuff) 

                        nset = size(rbuff)
                        pixelset%unum(ndsp+1:ndsp+nset) = rbuff

                        ndsp = ndsp + nset
                     ENDIF

                  ENDIF
               ENDDO
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
            iblk = landunit(iu)%xblk
            jblk = landunit(iu)%yblk
            DO iset = 1, pixelset%nset
               DO WHILE (pixelset%unum(iset) /= landunit(iu)%num)
                  iu = iu + 1
                  ju = ju + 1
                  IF ((landunit(iu)%xblk /= iblk) .or. (landunit(iu)%yblk /= jblk)) THEN
                     ju = 1
                     iblk = landunit(iu)%xblk
                     jblk = landunit(iu)%yblk
                  ENDIF
               ENDDO

               nave = nunits_blk(iblk,jblk) / (p_np_group-1)
               nres = mod(nunits_blk(iblk,jblk), p_np_group-1)
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
               DO WHILE (pixelset%unum(iset) /= landunit(iu)%num)
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
