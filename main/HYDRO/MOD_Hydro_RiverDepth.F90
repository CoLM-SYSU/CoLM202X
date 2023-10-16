#include <define.h>

#ifdef LATERAL_FLOW

MODULE MOD_Hydro_RiverDepth

   USE MOD_Precision
   IMPLICIT NONE

   real(r8), parameter :: cH_rivdpt   = 0.1
   real(r8), parameter :: pH_rivdpt   = 0.5
   real(r8), parameter :: B0_rivdpt   = 0.0
   real(r8), parameter :: Bmin_rivdpt = 1.0

CONTAINS

   SUBROUTINE calc_riverdepth_from_runoff ()
      
      USE MOD_SPMD_Task
      USE MOD_Namelist
      USE MOD_DataType
      USE MOD_NetCDFSerial
      USE MOD_NetCDFBlock
      USE MOD_Block
      USE MOD_Mesh
      USE MOD_Grid
      USE MOD_Mapping_Grid2Pset
      USE MOD_LandElm
      USE MOD_ElmVector
      USE MOD_Hydro_RiverLakeNetwork
      USE MOD_Hydro_HillslopeNetwork
      USE MOD_Hydro_IO
      IMPLICIT NONE

      ! Local Variables
      character(len=256) :: file_rnof, file_rivdpt
      type(grid_type)    :: grid_rnof
      type(block_data_real8_2d)    :: f_rnof
      type(mapping_grid2pset_type) :: mg2p_rnof

      real(r8), allocatable :: bsnrnof(:) , bsndis(:)
      integer,  allocatable :: nups_riv(:), iups_riv(:), b_up2down(:)

      integer :: i, j, ithis, ib, jb, iblkme
      integer :: iwork, mesg(2), isrc, ndata
      real(r8), allocatable :: rcache(:)


      file_rnof = trim(DEF_dir_runtime) // '/runoff_clim.nc'

      CALL grid_rnof%define_from_file (file_rnof, 'lat', 'lon')

      call mg2p_rnof%build (grid_rnof, landelm)

      IF (p_is_io) THEN
         CALL allocate_block_data (grid_rnof, f_rnof)
         CALL ncio_read_block (file_rnof, 'ro', grid_rnof, f_rnof)

         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)
            do j = 1, grid_rnof%ycnt(jb)
               do i = 1, grid_rnof%xcnt(ib)
                  f_rnof%blk(ib,jb)%val(i,j) = max(f_rnof%blk(ib,jb)%val(i,j), 0.)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (p_is_worker) THEN
         IF (numelm > 0) allocate (bsnrnof (numelm))
      ENDIF

      call mg2p_rnof%map_aweighted (f_rnof, bsnrnof)

      IF (p_is_worker) THEN
         IF (numelm > 0) THEN
            bsnrnof = bsnrnof /24.0/3600.0 ! from m/day to m/s
            DO i = 1, numelm
               ! total runoff in basin, from m/s to m3/s
               IF (lake_id(i) <= 0) THEN
                  bsnrnof(i) = bsnrnof(i) * sum(hillslope_network(i)%area)
               ELSE
                  bsnrnof(i) = bsnrnof(i) * sum(lakes(i)%area0)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      if (p_is_worker) then
         mesg = (/p_iam_glb, numelm/)
         call mpi_send (mesg, 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err) 
         IF (numelm > 0) THEN
            call mpi_send (bsnrnof, numelm, MPI_REAL8, p_root, mpi_tag_data, p_comm_glb, p_err) 
         ENDIF
      ENDIF
      
      IF (p_is_master) THEN
         
         allocate (bsnrnof (totalnumelm))

         DO iwork = 0, p_np_worker-1
            call mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndata))

               call mpi_recv (rcache, ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
               
               bsnrnof(elm_data_address(p_itis_worker(isrc))%val) = rcache

               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF
      
      CALL mpi_barrier (p_comm_glb, p_err)
#else
      bsnrnof(elm_data_address(0)%val) = bsnrnof
#endif


      IF (p_is_master) THEN

         allocate (nups_riv (totalnumelm))
         allocate (iups_riv (totalnumelm))
         allocate (b_up2down(totalnumelm))

         allocate (bsndis   (totalnumelm))

         nups_riv(:) = 0
         DO i = 1, totalnumelm
            j = riverdown(i)
            IF (j > 0) THEN
               nups_riv(j) = nups_riv(j) + 1
            ENDIF
         ENDDO

         iups_riv(:) = 0
         ithis = 0
         DO i = 1, totalnumelm
            IF (iups_riv(i) == nups_riv(i)) THEN
               ithis = ithis + 1
               b_up2down(ithis) = i

               j = riverdown(i)
               IF (j > 0) THEN
                  iups_riv(j) = iups_riv(j) + 1
                  DO WHILE (iups_riv(j) == nups_riv(j))
                     IF (j < i) THEN
                        ithis = ithis + 1
                        b_up2down(ithis) = j
                     ENDIF
                     j = riverdown(j)
                     IF (j > 0) THEN
                        iups_riv(j) = iups_riv(j) + 1
                     ELSE
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF
            ELSE
               CYCLE
            ENDIF
         ENDDO

         bsndis(:) = 0.
         DO i = 1, totalnumelm
            j = b_up2down(i)
            bsndis(j) = bsndis(j) + bsnrnof(j)
            IF (riverdown(j) > 0) THEN
               bsndis(riverdown(j)) = bsndis(riverdown(j)) + bsndis(j)
            ENDIF
         ENDDO

         DO i = 1, totalnumelm
            riverdpth(i) = max(cH_rivdpt * (bsndis(i)**pH_rivdpt) + B0_rivdpt, Bmin_rivdpt)
         ENDDO

         file_rivdpt = trim(DEF_dir_restart) // '/' // trim(DEF_CASE_NAME) //'_riverdepth.nc'
         CALL ncio_create_file      (file_rivdpt)
         CALL ncio_define_dimension (file_rivdpt, 'basin', totalnumelm)
         call ncio_write_serial     (file_rivdpt, 'riverdepth', riverdpth, 'basin')

      ENDIF

      CALL vector_read_basin (file_rivdpt, riverdpth, numelm, 'riverdepth', elm_data_address)

      IF (allocated (bsnrnof  )) deallocate(bsnrnof  )
      IF (allocated (bsndis   )) deallocate(bsndis   )
      IF (allocated (nups_riv )) deallocate(nups_riv )
      IF (allocated (iups_riv )) deallocate(iups_riv )
      IF (allocated (b_up2down)) deallocate(b_up2down)

   END SUBROUTINE calc_riverdepth_from_runoff

END MODULE MOD_Hydro_RiverDepth
#endif
