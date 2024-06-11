#include <define.h>

#ifdef CROP
MODULE MOD_LandCrop

!------------------------------------------------------------------------------------
! DESCRIPTION:
!
!    Build crop patches.
!
! Created by Shupeng Zhang, Sep 2023
!    porting codes from Hua Yuan's OpenMP version to MPI parallel version.
!------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   IMPLICIT NONE

   ! ---- Instance ----
   type(grid_type) :: gcrop
   integer,  allocatable :: cropclass (:)
   real(r8), allocatable :: pctshrpch (:)

CONTAINS

   ! -------------------------------
   SUBROUTINE landcrop_build (lc_year)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_DataType
   USE MOD_LandElm
#ifdef CATCHMENT
   USE MOD_LandHRU
#endif
   USE MOD_LandPatch
   USE MOD_NetCDFBlock
   USE MOD_PixelsetShared
   USE MOD_5x5DataReadin
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year

   ! Local Variables
   character(len=255) :: cyear, file_patch, dir_5x5, suffix
   integer :: npatch_glb
   type(block_data_real8_2d) :: pctcrop_xy
   type(block_data_real8_3d) :: pctshared_xy
   type(block_data_real8_3d) :: cropdata
   integer :: sharedfilter(1), cropfilter(1)
   integer :: iblkme, ib, jb
   real(r8), allocatable :: pctshared  (:)
   integer , allocatable :: classshared(:)

      write(cyear,'(i4.4)') lc_year
      IF (p_is_master) THEN
         write(*,'(A)') 'Making patches (crop shared) :'
      ENDIF

#if (defined SinglePoint && defined CROP)
      IF ((SITE_landtype == CROPLAND) .and. (USE_SITE_pctcrop)) THEN

         numpatch = count(SITE_pctcrop > 0.)

         allocate (pctshrpch (numpatch))
         allocate (cropclass(numpatch))
         cropclass = pack(SITE_croptyp, SITE_pctcrop > 0.)
         pctshrpch = pack(SITE_pctcrop, SITE_pctcrop > 0.)

         pctshrpch = pctshrpch / sum(pctshrpch)

         allocate (landpatch%eindex (numpatch))
         allocate (landpatch%ipxstt (numpatch))
         allocate (landpatch%ipxend (numpatch))
         allocate (landpatch%settyp (numpatch))
         allocate (landpatch%ielm   (numpatch))

         landpatch%eindex(:) = 1
         landpatch%ielm  (:) = 1
         landpatch%ipxstt(:) = 1
         landpatch%ipxend(:) = 1
         landpatch%settyp(:) = CROPLAND
         
         landpatch%has_shared = .true.
         allocate (landpatch%pctshared(numpatch))
         landpatch%pctshared = pctshrpch

         landpatch%nset = numpatch
         CALL landpatch%set_vecgs

         RETURN
      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN

         dir_5x5 = trim(DEF_dir_rawdata) // '/plant_15s'
         suffix  = 'MOD'//trim(cyear)

         CALL allocate_block_data (gpatch, pctcrop_xy)
         CALL read_5x5_data (dir_5x5, suffix, gpatch, 'PCT_CROP', pctcrop_xy)
         
         CALL allocate_block_data (gpatch, pctshared_xy, 2)
         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)
            pctshared_xy%blk(ib,jb)%val(1,:,:) = 1. - pctcrop_xy%blk(ib,jb)%val/100.
            pctshared_xy%blk(ib,jb)%val(2,:,:) = pctcrop_xy%blk(ib,jb)%val/100.
         ENDDO
      ENDIF
      
      sharedfilter = (/ 1 /)

      CALL pixelsetshared_build (landpatch, gpatch, pctshared_xy, 2, sharedfilter, &
         pctshared, classshared)

      IF (p_is_worker) THEN
         IF (landpatch%nset > 0) THEN
            WHERE (classshared == 2) landpatch%settyp = CROPLAND
         ENDIF
      ENDIF

      IF (p_is_io) THEN
         file_patch = trim(DEF_dir_rawdata) // '/global_CFT_surface_data.nc'
         CALL allocate_block_data (gcrop, cropdata, N_CFT)
         CALL ncio_read_block (file_patch, 'PCT_CFT', gcrop, N_CFT, cropdata)
      ENDIF

      cropfilter = (/ CROPLAND /)
      
      CALL pixelsetshared_build (landpatch, gcrop, cropdata, N_CFT, cropfilter, &
         pctshrpch, cropclass, fracin = pctshared)
      
      numpatch = landpatch%nset

      landpatch%has_shared = .true.
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate(landpatch%pctshared(numpatch))
            landpatch%pctshared = pctshrpch
         ENDIF
      ENDIF

      IF (allocated(pctshared  )) deallocate(pctshared  )
      IF (allocated(classshared)) deallocate(classshared)

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numpatch, npatch_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npatch_glb, ' patches (with crop).'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numpatch, ' patches.'
#endif

      CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
#ifdef CATCHMENT
      CALL hru_patch%build (landhru, landpatch, use_frac = .true.)
#endif

      CALL write_patchfrac (DEF_dir_landdata, lc_year)
   
   END SUBROUTINE landcrop_build

END MODULE MOD_LandCrop
#endif
