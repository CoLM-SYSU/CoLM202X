#include <define.h>

#ifdef PFT_CLASSIFICATION

MODULE mod_landpft

   USE mod_namelist
   USE mod_pixelset
   USE LC_const
   USE GlobalVars
   IMPLICIT NONE

   ! ---- Instance ----
   INTEGER :: numpft
   TYPE(pixelset_type) :: landpft

   INTEGER , allocatable :: pft2patch   (:)  !patch index of a PFT
   INTEGER , allocatable :: patch_pft_s (:)  !start PFT index of a patch
   INTEGER , allocatable :: patch_pft_e (:)  !end PFT index of a patch

   ! ---- PUBLIC routines ----
   PUBLIC :: landpft_build

CONTAINS

   ! -------------------------------
   SUBROUTINE landpft_build (lc_year)

      USE precision
      USE spmd_task
      USE mod_grid
      USE mod_data_type
      USE mod_namelist
      USE mod_5x5_data
      USE mod_landpatch
      USE mod_aggregation
      USE LC_const

      IMPLICIT NONE

      INTEGER, intent(in) :: lc_year
      ! Local Variables
      CHARACTER(len=256) :: dir_5x5, suffix, cyear
      TYPE (block_data_real8_3d) :: pctpft
      REAL(r8), allocatable :: pctpft_patch(:,:), pctpft_one(:,:)
      REAL(r8), allocatable :: area_one(:)
      INTEGER  :: ipatch, ipft, npatch, npft
      REAL(r8) :: sumarea
      LOGICAL, allocatable :: patchmask (:)
      INTEGER  :: npft_glb

      ! add parameter input for time year
      write(cyear,'(i4.4)') lc_year
      IF (p_is_master) THEN
         write(*,'(A)') 'Making land plant function type tiles :'
      ENDIF

#ifdef SinglePoint
      IF (USE_SITE_pctpfts) THEN
         IF (landpatch%settyp(1) == 1) THEN
            numpft = count(SITE_pctpfts > 0.)
#ifdef CROP
         ELSEIF (landpatch%settyp(1) == CROPLAND) THEN
            numpft = numpatch
#endif
         ELSE
            numpft = 0
         ENDIF

         allocate (patch_pft_s (numpatch))
         allocate (patch_pft_e (numpatch))

         IF (numpft > 0) THEN
            allocate (landpft%eindex (numpft))
            allocate (landpft%settyp (numpft))
            allocate (landpft%ipxstt (numpft))
            allocate (landpft%ipxend (numpft))
            allocate (landpft%ielm   (numpft))

            landpft%ielm  (:) = 1
            landpft%eindex(:) = 1
            landpft%ipxstt(:) = 1
            landpft%ipxend(:) = 1

            allocate(pft2patch (numpft))

            IF (landpatch%settyp(1) == 1) THEN
               landpft%settyp = pack(SITE_pfttyp, SITE_pctpfts > 0.)

               pft2patch  (:) = 1
               patch_pft_s(:) = 1
               patch_pft_e(:) = numpft
#ifdef CROP
            ELSEIF (landpatch%settyp(1) == CROPLAND) THEN
               DO ipft = 1, numpft
                  landpft%settyp(ipft) = cropclass(ipft) + N_PFT - 1
                  pft2patch   (ipft) = ipft
                  patch_pft_s (ipft) = ipft
                  patch_pft_e (ipft) = ipft
               ENDDO
#endif
            ENDIF
         ELSE
            write(*,*) 'Warning : land type ', landpatch%settyp(1), ' for PFT_CLASSIFICATION'
            patch_pft_s(:) = -1
            patch_pft_e(:) = -1
         ENDIF

         landpft%nset = numpft
         CALL landpft%set_vecgs

         RETURN
      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      if (p_is_io) then

         call allocate_block_data (gpatch, pctpft, N_PFT_modis, lb1 = 0)
         CALL flush_block_data (pctpft, 1.0)

         dir_5x5 = trim(DEF_dir_rawdata) // '/plant_15s_clim'
         suffix  = 'MOD'//trim(cyear)
         CALL read_5x5_data_pft (dir_5x5, suffix, gpatch, 'PCT_PFT', pctpft)

#ifdef USEMPI
         CALL aggregation_data_daemon (gpatch, data_r8_3d_in1 = pctpft, n1_r8_3d_in1 = N_PFT_modis)
#endif
      end if


      if (p_is_worker) then

         IF (numpatch > 0) THEN
            allocate (pctpft_patch (0:N_PFT-1,numpatch))
            allocate (patchmask (numpatch))

            pctpft_patch(:,:) = 0
            patchmask (:) = .true.
         ENDIF

         DO ipatch = 1, numpatch
            IF (landpatch%settyp(ipatch) == 1) THEN

               CALL aggregation_request_data (landpatch, ipatch, gpatch, area = area_one, &
                  data_r8_3d_in1 = pctpft, data_r8_3d_out1 = pctpft_one, n1_r8_3d_in1 = N_PFT_modis, lb1_r8_3d_in1 = 0)

               sumarea = sum(area_one)

               DO ipft = 0, N_PFT-1
                  pctpft_patch(ipft,ipatch) = sum(pctpft_one(ipft,:) * area_one) / sumarea
               ENDDO

               IF (sum(pctpft_patch(:,ipatch)) <= 0.0) THEN
                  patchmask(ipatch) = .false.
               ENDIF
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif

         IF (numpatch > 0) THEN
            npatch = count(patchmask)
            numpft = count(pctpft_patch > 0.)
#ifdef CROP
            numpft = numpft + count(landpatch%settyp == CROPLAND)
#endif
            IF (npatch > 0) THEN
               allocate (patch_pft_s (npatch))
               allocate (patch_pft_e (npatch))
            ENDIF
         ELSE
            numpft = 0
         ENDIF

         IF (numpft > 0) THEN

            allocate (pft2patch   (numpft))

            allocate (landpft%eindex (numpft))
            allocate (landpft%settyp (numpft))
            allocate (landpft%ipxstt (numpft))
            allocate (landpft%ipxend (numpft))
            allocate (landpft%ielm   (numpft))

            npft = 0
            npatch = 0
            DO ipatch = 1, numpatch
               IF (patchmask(ipatch)) THEN
                  npatch = npatch + 1

                  IF (landpatch%settyp(ipatch) == 1) THEN
                     patch_pft_s(npatch) = npft + 1
                     patch_pft_e(npatch) = npft + count(pctpft_patch(:,ipatch) > 0)

                     DO ipft = 0, N_PFT-1
                        IF (pctpft_patch(ipft,ipatch) > 0) THEN
                           npft = npft + 1

                           landpft%ielm  (npft) = landpatch%ielm  (ipatch)
                           landpft%eindex(npft) = landpatch%eindex(ipatch)
                           landpft%ipxstt(npft) = landpatch%ipxstt(ipatch)
                           landpft%ipxend(npft) = landpatch%ipxend(ipatch)
                           landpft%settyp(npft) = ipft

                           pft2patch(npft) = npatch
                        ENDIF
                     ENDDO
#ifdef CROP
                  ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
                     npft = npft + 1
                     patch_pft_s(npatch) = npft
                     patch_pft_e(npatch) = npft

                     landpft%ielm  (npft) = landpatch%ielm  (ipatch)
                     landpft%eindex(npft) = landpatch%eindex(ipatch)
                     landpft%ipxstt(npft) = landpatch%ipxstt(ipatch)
                     landpft%ipxend(npft) = landpatch%ipxend(ipatch)
                     landpft%settyp(npft) = cropclass(ipatch) + N_PFT - 1

                     pft2patch(npft) = npatch
#endif
                  ELSE
                     patch_pft_s(npatch) = -1
                     patch_pft_e(npatch) = -1
                  ENDIF
               ENDIF
            ENDDO

         ENDIF

      ENDIF

      CALL landpatch%pset_pack(patchmask, numpatch)

      landpft%nset = numpft
      CALL landpft%set_vecgs

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numpft, npft_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npft_glb, ' plant function type tiles.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numpft, ' plant function type tiles.'
#endif

#ifdef SinglePoint
      allocate  (SITE_pfttyp(numpft))
      SITE_pfttyp(:) = landpft%settyp
#endif

      IF (allocated(pctpft_patch)) deallocate (pctpft_patch)
      IF (allocated(pctpft_one  )) deallocate (pctpft_one  )
      IF (allocated(area_one    )) deallocate (area_one    )
      IF (allocated(patchmask   )) deallocate (patchmask   )

   END SUBROUTINE landpft_build

   ! ----------------------
   SUBROUTINE map_patch_to_pft

      USE spmd_task
      USE mod_landpatch
      IMPLICIT NONE

      INTEGER :: ipatch, ipft

      IF (p_is_worker) THEN

         IF ((numpatch <= 0) .or. (numpft <= 0)) return

         IF (allocated(patch_pft_s)) deallocate(patch_pft_s)
         IF (allocated(patch_pft_e)) deallocate(patch_pft_e)
         IF (allocated(pft2patch  )) deallocate(pft2patch  )

         allocate (patch_pft_s (numpatch))
         allocate (patch_pft_e (numpatch))
         allocate (pft2patch   (numpft  ))

         ipft = 1
         DO ipatch = 1, numpatch
            IF (landpatch%settyp(ipatch) == 1) THEN

               patch_pft_s(ipatch) = ipft

               DO WHILE (ipft <= numpft)
                  IF ((landpft%eindex(ipft) == landpatch%eindex(ipatch))  &
                     .and. (landpft%ipxstt(ipft) == landpatch%ipxstt(ipatch))) THEN
                     pft2patch  (ipft  ) = ipatch
                     patch_pft_e(ipatch) = ipft
                     ipft = ipft + 1
                  ELSE
                     exit
                  ENDIF
               ENDDO
#ifdef CROP
            ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
               patch_pft_s(ipatch) = ipft
               patch_pft_e(ipatch) = ipft
               pft2patch  (ipft  ) = ipatch
               ipft = ipft + 1
#endif
            ELSE
               patch_pft_s(ipatch) = -1
               patch_pft_e(ipatch) = -1
            ENDIF

         ENDDO

      ENDIF

   END SUBROUTINE map_patch_to_pft

END MODULE mod_landpft
#endif
