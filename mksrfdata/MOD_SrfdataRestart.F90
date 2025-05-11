#include <define.h>

MODULE MOD_SrfdataRestart
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    This module includes subroutines to read/write data of mesh and pixelsets.
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   IMPLICIT NONE

   ! ----- subroutines -----
   PUBLIC :: mesh_save_to_file
   PUBLIC :: mesh_load_from_file

   PUBLIC :: pixelset_save_to_file
   PUBLIC :: pixelset_load_from_file

CONTAINS

   ! -----------------------
   SUBROUTINE mesh_save_to_file (dir_landdata, lc_year)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   USE MOD_Block
   USE MOD_Utils
   IMPLICIT NONE

   character(len=*), intent(in) :: dir_landdata
   integer         , intent(in) :: lc_year

   ! Local variables
   character(len=256) :: filename, fileblock, cyear
   integer :: ie, je, nelm, totlen, tothis, iblk, jblk, iworker, i
   integer,   allocatable :: nelm_worker(:), ndsp_worker(:)
   integer*8, allocatable :: elmindx(:)
   integer,   allocatable :: npxlall(:)
   integer,   allocatable :: elmpixels(:,:)
   real(r8),  allocatable :: lon(:), lat(:)

   integer :: nsend, nrecv, ndone, ndsp

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
                  totlen = 0
                  DO ie = 1, numelm
                     IF ((mesh(ie)%xblk == iblk) .and. (mesh(ie)%yblk == jblk)) THEN
                        nelm = nelm + 1
                        totlen = totlen + mesh(ie)%npxl
                     ENDIF
                  ENDDO

                  IF (nelm > 0) THEN

                     allocate (elmindx (nelm))
                     allocate (npxlall (nelm))
                     allocate (elmpixels (2,totlen))

                     je = 0
                     ndsp = 0
                     DO ie = 1, numelm
                        IF ((mesh(ie)%xblk == iblk) .and. (mesh(ie)%yblk == jblk)) THEN
                           je = je + 1
                           elmindx(je) = mesh(ie)%indx
                           npxlall(je) = mesh(ie)%npxl

                           elmpixels(1,ndsp+1:ndsp+npxlall(je)) = mesh(ie)%ilon
                           elmpixels(2,ndsp+1:ndsp+npxlall(je)) = mesh(ie)%ilat

                           ndsp = ndsp + npxlall(je)
                        ENDIF
                     ENDDO
                  ENDIF

#ifdef USEMPI
                  CALL mpi_gather (nelm, 1, MPI_INTEGER, &
                     MPI_INULL_P, 1, MPI_INTEGER, p_root, p_comm_group, p_err)

                  CALL mpi_gatherv (elmindx, nelm, MPI_INTEGER8, &
                     MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER8, & ! insignificant on workers
                     p_root, p_comm_group, p_err)

                  CALL mpi_gatherv (npxlall, nelm, MPI_INTEGER, &
                     MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER, & ! insignificant on workers
                     p_root, p_comm_group, p_err)

                  ndone = 0
                  DO WHILE (ndone < totlen)
                     nsend = max(min(totlen-ndone, MesgMaxSize/8), 1)
                     CALL mpi_send (nsend, 1, &
                        MPI_INTEGER, p_root, mpi_tag_size, p_comm_group, p_err)
                     CALL mpi_send (elmpixels(:,ndone+1:ndone+nsend), 2*nsend, &
                        MPI_INTEGER, p_root, mpi_tag_data, p_comm_group, p_err)
                     ndone = ndone + nsend
                  ENDDO
               ENDIF
            ENDIF
#endif

#ifdef USEMPI
            IF (p_is_io) THEN
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

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
                  CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER8, &
                     elmindx, nelm_worker(0:), ndsp_worker(0:), MPI_INTEGER8, &
                     p_root, p_comm_group, p_err)

                  allocate (npxlall (nelm))
                  CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                     npxlall, nelm_worker(0:), ndsp_worker(0:), MPI_INTEGER, &
                     p_root, p_comm_group, p_err)

                  totlen = sum(npxlall)
                  allocate (elmpixels (2, totlen))

                  ndone = 0
                  DO iworker = 1, p_np_group-1

                     ndsp   = ndsp_worker(iworker)
                     tothis = ndone + sum(npxlall(ndsp+1:ndsp+nelm_worker(iworker)))

                     DO WHILE (ndone < tothis)

                        CALL mpi_recv (nrecv, 1, &
                           MPI_INTEGER, iworker, mpi_tag_size, p_comm_group, p_stat, p_err)
                        CALL mpi_recv (elmpixels(:,ndone+1:ndone+nrecv), 2*nrecv, &
                           MPI_INTEGER, iworker, mpi_tag_data, p_comm_group, p_stat, p_err)

                        ndone = ndone + nrecv
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
                     CALL ncio_define_dimension (fileblock, 'ncoor',  2   )
                     CALL ncio_define_dimension (fileblock, 'pixel',  totlen)

                     CALL ncio_write_serial (fileblock, 'elmindex',  elmindx, 'element')
                     CALL ncio_write_serial (fileblock, 'elmnpxl',   npxlall, 'element')
                     CALL ncio_write_serial (fileblock, 'elmpixels', elmpixels, &
                        'ncoor', 'pixel', compress = 1)
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
   SUBROUTINE mesh_load_from_file (dir_landdata, lc_year)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   IMPLICIT NONE

   integer         , intent(in) :: lc_year
   character(len=*), intent(in) :: dir_landdata

   ! Local variables
   character(len=256) :: filename, fileblock, cyear
   integer :: iblkme, iblk, jblk, ie, nelm, ndsp, pdsp
   integer*8, allocatable :: elmindx(:)
   integer,   allocatable :: datasize(:)
   integer,   allocatable :: npxl(:), pixels(:,:), pixels2d(:,:,:)

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

         IF (numelm > 0) THEN

            allocate (mesh (numelm))

            ndsp = 0
            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               nelm = nelm_blk(iblk,jblk)

               IF (nelm > 0) THEN

                  CALL get_filename_block (filename, iblk, jblk, fileblock)
                  CALL ncio_read_serial (fileblock, 'elmindex',  elmindx)
                  CALL ncio_read_serial (fileblock, 'elmnpxl',   npxl   )

                  CALL ncio_inquire_varsize (fileblock, 'elmpixels', datasize)
                  IF (size(datasize) == 3) THEN
                     CALL ncio_read_serial (fileblock, 'elmpixels', pixels2d)
                  ELSE
                     CALL ncio_read_serial (fileblock, 'elmpixels', pixels)
                  ENDIF

                  pdsp = 0
                  DO ie = 1, nelm
                     mesh(ie+ndsp)%indx = elmindx(ie)
                     mesh(ie+ndsp)%npxl = npxl(ie)
                     mesh(ie+ndsp)%xblk = iblk
                     mesh(ie+ndsp)%yblk = jblk

                     allocate (mesh(ie+ndsp)%ilon (npxl(ie)))
                     allocate (mesh(ie+ndsp)%ilat (npxl(ie)))

                     IF (size(datasize) == 3) THEN
                        mesh(ie+ndsp)%ilon = pixels2d(1,1:npxl(ie),ie)
                        mesh(ie+ndsp)%ilat = pixels2d(2,1:npxl(ie),ie)
                     ELSE
                        mesh(ie+ndsp)%ilon = pixels(1,pdsp+1:pdsp+npxl(ie))
                        mesh(ie+ndsp)%ilat = pixels(2,pdsp+1:pdsp+npxl(ie))
                        pdsp = pdsp + npxl(ie)
                     ENDIF
                  ENDDO

                  ndsp = ndsp + nelm
               ENDIF
            ENDDO
         ENDIF

         IF (allocated(elmindx ))  deallocate(elmindx )
         IF (allocated(npxl    ))  deallocate(npxl    )
         IF (allocated(datasize))  deallocate(datasize)
         IF (allocated(pixels  ))  deallocate(pixels  )
         IF (allocated(pixels2d))  deallocate(pixels2d)

      ENDIF

#ifdef CoLMDEBUG
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
   SUBROUTINE pixelset_save_to_file (dir_landdata, psetname, pixelset, lc_year)

   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_NetCDFVector
   USE MOD_Pixelset
   IMPLICIT NONE

   character(len=*),    intent(in) :: dir_landdata
   character(len=*),    intent(in) :: psetname
   type(pixelset_type), intent(in) :: pixelset
   integer         ,    intent(in) :: lc_year

   ! Local variables
   character(len=256)   :: filename, cyear

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

      CALL ncio_write_vector (filename, 'eindex', trim(psetname), pixelset, pixelset%eindex, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (filename, 'ipxstt', trim(psetname), pixelset, pixelset%ipxstt, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (filename, 'ipxend', trim(psetname), pixelset, pixelset%ipxend, DEF_Srfdata_CompressLevel)
      CALL ncio_write_vector (filename, 'settyp', trim(psetname), pixelset, pixelset%settyp, DEF_Srfdata_CompressLevel)

      IF (pixelset%has_shared) THEN
         CALL ncio_write_vector (filename, 'pctshared', trim(psetname), pixelset, pixelset%pctshared, DEF_Srfdata_CompressLevel)
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) write(*,*) 'SAVE Pixel Sets ' // trim(psetname) // ' done.'

   END SUBROUTINE pixelset_save_to_file


   !---------------------------
   SUBROUTINE pixelset_load_from_file (dir_landdata, psetname, pixelset, numset, lc_year)

   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_NetCDFSerial
   USE MOD_NetCDFVector
   USE MOD_Mesh
   USE MOD_Pixelset
   IMPLICIT NONE

   integer         ,    intent(in) :: lc_year
   character(len=*),    intent(in) :: dir_landdata
   character(len=*),    intent(in) :: psetname
   type(pixelset_type), intent(inout) :: pixelset
   integer, intent(out) :: numset

   ! Local variables
   character(len=256) :: filename, fileblock, blockname, cyear
   integer :: iset, nset, ndsp, iblkme, iblk, jblk, ie, je, nave, nres, left, iproc
   integer :: nsend, nrecv
   integer*8, allocatable :: rbuff(:), sbuff(:)
   integer,   allocatable :: iworker(:)
   logical,   allocatable :: msk(:)
   logical :: fexists, fexists_any

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

         fexists_any = .false.

         DO iblkme = 1, gblock%nblkme
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

#if (defined VectorInOneFileS || defined VectorInOneFileP)
            CALL get_blockname (iblk, jblk, blockname)
            CALL ncio_inquire_length_grp (filename, 'eindex', &
               trim(psetname)//'_'//trim(blockname), nset)
            pixelset%nset = pixelset%nset + nset
#else
            CALL get_filename_block (filename, iblk, jblk, fileblock)

            inquire (file=trim(fileblock), exist=fexists)
            IF (fexists) THEN
               CALL ncio_inquire_length (fileblock, 'eindex', nset)
               pixelset%nset = pixelset%nset + nset
            ENDIF

            fexists_any = fexists_any .or. fexists
#endif
         ENDDO

#if (defined VectorInOneFileS || defined VectorInOneFileP)
         fexists_any = pixelset%nset > 0
#endif

#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, fexists_any, 1, MPI_LOGICAL, MPI_LOR, p_comm_io, p_err)
#endif
         IF (.not. fexists_any) THEN
            write(*,*) 'Warning : restart file ' //trim(filename)// ' not found.'
            CALL CoLM_stop ()
         ENDIF

         IF (pixelset%nset > 0) THEN

            allocate (pixelset%eindex (pixelset%nset))

            ndsp = 0
            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

#if (defined VectorInOneFileS || defined VectorInOneFileP)
               CALL get_blockname (iblk, jblk, blockname)
               CALL ncio_inquire_length_grp (filename, 'eindex', &
                  trim(psetname)//'_'//trim(blockname), nset)

               IF (nset > 0) THEN

                  CALL ncio_read_serial_grp_int64_1d (filename, 'eindex', &
                     trim(psetname)//'_'//trim(blockname), rbuff)

                  pixelset%eindex(ndsp+1:ndsp+nset) = rbuff

                  ndsp = ndsp + nset
               ENDIF
#else
               CALL get_filename_block (filename, iblk, jblk, fileblock)
               inquire (file=trim(fileblock), exist=fexists)
               IF (fexists) THEN

                  CALL ncio_read_serial (fileblock, 'eindex', rbuff)

                  nset = size(rbuff)
                  pixelset%eindex(ndsp+1:ndsp+nset) = rbuff

                  ndsp = ndsp + nset
               ENDIF
#endif

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
                  CALL mpi_send (sbuff, nsend, MPI_INTEGER8, iproc, mpi_tag_data, p_comm_group, p_err)
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
            CALL mpi_recv (pixelset%eindex, nrecv, MPI_INTEGER8, &
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

         ELSE
            write(*,*) 'Warning: 0 ',trim(psetname), ' on worker :', p_iam_glb
         ENDIF
      ENDIF

      numset = pixelset%nset

      pixelset%has_shared = .false.
      IF (p_is_worker) THEN
         DO iset = 1, pixelset%nset-1
            IF ((pixelset%ielm(iset) == pixelset%ielm(iset+1)) &
               .and. (pixelset%ipxstt(iset) == pixelset%ipxstt(iset+1))) THEN
               pixelset%has_shared = .true.
               exit
            ENDIF
         ENDDO
      ENDIF

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, pixelset%has_shared, 1, MPI_LOGICAL, &
         MPI_LOR, p_comm_glb, p_err)
#endif

      IF (pixelset%has_shared) THEN
         CALL ncio_read_vector (filename, 'pctshared', pixelset, pixelset%pctshared)
      ENDIF

#ifdef CoLMDEBUG
      IF (p_is_io)  write(*,*) numset, trim(psetname), ' on group', p_iam_io
#endif

   END SUBROUTINE pixelset_load_from_file

END MODULE MOD_SrfdataRestart
