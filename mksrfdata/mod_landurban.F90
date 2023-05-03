#include <define.h>

#ifdef URBAN_MODEL

MODULE mod_landurban

   USE mod_grid
   USE mod_pixelset
   USE GlobalVars, only : N_URB
   IMPLICIT NONE

   ! ---- Instance ----
   TYPE(grid_type) :: gurban

   INTEGER :: numurban
   TYPE(pixelset_type) :: landurban

   INTEGER , allocatable :: urban2patch (:)  !patch index of a urban
   INTEGER , allocatable :: patch2urban (:)  !urban index of a patch

   ! ---- PUBLIC routines ----
   PUBLIC :: landurban_build
   PUBLIC :: map_patch_to_urban

CONTAINS

   ! -------------------------------
   SUBROUTINE landurban_build ()

      USE precision
      USE spmd_task
      USE mod_grid
      USE mod_data_type
      USE mod_namelist
      USE mod_5x5_data
      USE mod_mesh
      USE mod_landpatch
      USE mod_aggregation
      USE mod_utils

      IMPLICIT NONE

      ! Local Variables
      CHARACTER(len=256) :: dir_urban
      TYPE (block_data_int32_2d) :: data_ur_class
      INTEGER, allocatable :: ibuff(:), types(:), order(:)

      INTEGER :: ipatch, jpatch, iurban
      INTEGER :: ie, ipxstt, ipxend, npxl, ipxl
      INTEGER :: nurb_glb

      INTEGER :: numpatch_
      INTEGER, allocatable :: eindex_(:)
      INTEGER, allocatable :: ipxstt_(:)
      INTEGER, allocatable :: ipxend_(:)
      INTEGER, allocatable :: settyp_(:)
      INTEGER, allocatable :: ielm_  (:)

      INTEGER :: numurban_
      INTEGER, allocatable :: urbclass (:)

      CHARACTER(len=256) :: suffix

      IF (p_is_master) THEN
         write(*,'(A)') 'Making urban type tiles :'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      if (p_is_io) then

         dir_urban = trim(DEF_dir_rawdata) // '/urb'

         CALL allocate_block_data (gurban, data_ur_class)
         CALL flush_block_data (data_ur_class, 0)

         suffix = 'URB2015'
#ifdef USE_LCZ
         CALL read_5x5_data (dir_urban, suffix, gurban, 'LCZ', data_ur_class)
#else
         CALL read_5x5_data (dir_urban, suffix, gurban, 'URBAN_DENSITY_CLASS', data_ur_class)
#endif

#ifdef USEMPI
         CALL aggregationn_data_daemon (gurban, data_i4_2d_in1 = data_ur_class)
#endif
      end if

      if (p_is_worker) then

         numpatch_ = numpatch + count(landpatch%settyp == 13) * (N_URB-1)

         allocate (eindex_(numpatch_))
         allocate (ipxstt_(numpatch_))
         allocate (ipxend_(numpatch_))
         allocate (settyp_(numpatch_))
         allocate (ielm_  (numpatch_))

         numurban_ = count(landpatch%settyp == 13) * N_URB
         IF (numurban_ > 0) THEN
            allocate (urbclass(numurban_))
         ENDIF

         jpatch = 0
         DO ipatch = 1, numpatch
            IF (landpatch%settyp(ipatch) == 13) THEN

               ie     = landpatch%ielm  (ipatch)
               ipxstt = landpatch%ipxstt(ipatch)
               ipxend = landpatch%ipxend(ipatch)

               CALL aggregation_request_data (landpatch, ipatch, gurban, &
                  data_i4_2d_in1 = data_ur_class, data_i4_2d_out1 = ibuff)

#ifndef USE_LCZ
               where (ibuff < 1 .or. ibuff > 3)
                  ibuff = 3
               END where
#endif

               npxl = ipxend - ipxstt + 1

               types(:) = ibuff
               deallocate (ibuff)
               allocate (types (ipxstt:ipxend))

               allocate (order (ipxstt:ipxend))
               order = (/ (ipxl, ipxl = ipxstt, ipxend) /)

               CALL quicksort (npxl, types, order)

               mesh(ie)%ilon(ipxstt:ipxend) = mesh(ie)%ilon(order)
               mesh(ie)%ilat(ipxstt:ipxend) = mesh(ie)%ilat(order)

               DO ipxl = ipxstt, ipxend
                  IF (ipxl /= ipxstt) THEN
                     IF (types(ipxl) /= types(ipxl-1)) THEN
                        ipxend_(jpatch) = ipxl - 1
                     ELSE
                        cycle
                     ENDIF
                  ENDIF

                  jpatch = jpatch + 1
                  eindex_(jpatch) = mesh(ie)%indx
                  settyp_(jpatch) = 13
                  ipxstt_(jpatch) = ipxl
                  ielm_  (jpatch) = ie

                  iurban = iurban + 1
                  urbclass(iurban) = types(ipxl)
               ENDDO
               ipxend_(jpatch) = ipxend

            ELSE
               jpatch = jpatch + 1
               eindex_(jpatch) = landpatch%eindex(ipatch)
               ipxstt_(jpatch) = landpatch%ipxstt(ipatch)
               ipxend_(jpatch) = landpatch%ipxend(ipatch)
               settyp_(jpatch) = landpatch%settyp(ipatch)
               ielm_  (jpatch) = landpatch%ielm  (ipatch)
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif

         numpatch = jpatch

         deallocate (landpatch%eindex)
         deallocate (landpatch%ipxstt)
         deallocate (landpatch%ipxend)
         deallocate (landpatch%settyp)
         deallocate (landpatch%ielm  )

         allocate (landpatch%eindex (numpatch))
         allocate (landpatch%ipxstt (numpatch))
         allocate (landpatch%ipxend (numpatch))
         allocate (landpatch%settyp (numpatch))
         allocate (landpatch%ielm   (numpatch))

         landpatch%eindex = eindex_(1:jpatch)
         landpatch%ipxstt = ipxstt_(1:jpatch)
         landpatch%ipxend = ipxend_(1:jpatch)
         landpatch%settyp = settyp_(1:jpatch)
         landpatch%ielm   = ielm_  (1:jpatch)

         IF (numpatch > 0) THEN
            numurban = count(landpatch%settyp == 13)
         ELSE
            numurban = 0
         ENDIF

         IF (numurban > 0) THEN
            allocate (landurban%eindex (numurban))
            allocate (landurban%settyp (numurban))
            allocate (landurban%ipxstt (numurban))
            allocate (landurban%ipxend (numurban))
            allocate (landurban%ielm   (numurban))

            landurban%eindex = pack(landpatch%eindex, landpatch%settyp == 13)
            landurban%ipxstt = pack(landpatch%ipxstt, landpatch%settyp == 13)
            landurban%ipxend = pack(landpatch%ipxend, landpatch%settyp == 13)
            landurban%ielm   = pack(landpatch%ielm  , landpatch%settyp == 13)

            landurban%settyp = urbclass(1:numurban)
         ENDIF

      ENDIF

      CALL landpatch%set_vecgs
      CALL landurban%set_vecgs

      CALL map_patch_to_urban

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numurban, nurb_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', nurb_glb, ' urban tiles.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numurban, ' urban tiles.'
#endif

      IF (allocated(ibuff)) deallocate (ibuff)
      IF (allocated(types)) deallocate (types)
      IF (allocated(order)) deallocate (order)

      IF (allocated(eindex_)) deallocate (eindex_)
      IF (allocated(ipxstt_)) deallocate (ipxstt_)
      IF (allocated(ipxend_)) deallocate (ipxend_)
      IF (allocated(settyp_)) deallocate (settyp_)
      IF (allocated(ielm_  )) deallocate (ielm_  )

      IF (allocated(urbclass)) deallocate (urbclass)

   END SUBROUTINE landurban_build

   ! ----------------------
   SUBROUTINE map_patch_to_urban

      USE spmd_task
      USE mod_landpatch
      IMPLICIT NONE

      INTEGER :: ipatch, iurban

      IF (p_is_worker) THEN

         IF ((numpatch <= 0) .or. (numurban <= 0)) return

         allocate (patch2urban (numpatch))
         allocate (urban2patch (numurban))

         iurban = 0
         DO ipatch = 1, numpatch
            IF (landpatch%settyp(ipatch) == 13) THEN
               iurban = iurban + 1
               patch2urban(ipatch) = iurban
               urban2patch(iurban) = ipatch
            ELSE
               patch2urban(ipatch) = -1
            ENDIF
         ENDDO

      ENDIF

   END SUBROUTINE map_patch_to_urban

END MODULE mod_landurban
#endif
