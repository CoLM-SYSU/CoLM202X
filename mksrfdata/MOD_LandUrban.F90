#include <define.h>

MODULE MOD_LandUrban

   !------------------------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !    Build pixelset "landurban".
   !
   ! Original authors: Hua Yuan and Wenzong Dong, 2022, OpenMP version.
   !
   ! REVISIONS:
   ! Wenzong Dong, Hua Yuan, Shupeng Zhang, 05/2023: porting codes to MPI parallel version
   !------------------------------------------------------------------------------------

   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_Vars_Global, only: N_URB, URBAN
   IMPLICIT NONE

   ! ---- Instance ----
   TYPE(grid_type) :: gurban

   INTEGER :: numurban
   TYPE(pixelset_type) :: landurban

   INTEGER , allocatable :: urban_reg   (:)  !region index of a urban
   INTEGER , allocatable :: urban2patch (:)  !patch index of a urban
   INTEGER , allocatable :: patch2urban (:)  !urban index of a patch

   ! ---- PUBLIC routines ----
   PUBLIC :: landurban_build
   PUBLIC :: map_patch_to_urban

CONTAINS

   ! -------------------------------
   SUBROUTINE landurban_build (lc_year)

      USE MOD_Precision
      USE MOD_SPMD_Task
      USE MOD_NetCDFBlock
      USE MOD_Grid
      USE MOD_DataType
      USE MOD_Namelist
      USE MOD_5x5DataReadin
      USE MOD_Mesh
      USE MOD_LandPatch
      USE MOD_LandElm
#ifdef CATCHMENT
      USE MOD_LandHRU
#endif
#if (defined CROP)
      USE MOD_PixelsetShadow
#endif
      USE MOD_AggregationRequestData
      USE MOD_Utils

      IMPLICIT NONE

      INTEGER, intent(in) :: lc_year
      ! Local Variables
      CHARACTER(len=256) :: dir_urban
      TYPE (block_data_int32_2d) :: data_urb_class ! urban type index

#if (defined CROP)
      TYPE(block_data_real8_3d) :: cropdata
      INTEGER            :: cropfilter(1)
      CHARACTER(len=256) :: file_patch
#endif
      ! local vars
      INTEGER, allocatable :: ibuff(:), types(:), order(:)

      ! index
      INTEGER :: ipatch, jpatch, iurban
      INTEGER :: ie, ipxstt, ipxend, npxl, ipxl
      INTEGER :: nurb_glb, npatch_glb

      ! local vars for landpath and landurban
      INTEGER :: numpatch_
      INTEGER, allocatable :: eindex_(:)
      INTEGER, allocatable :: ipxstt_(:)
      INTEGER, allocatable :: ipxend_(:)
      INTEGER, allocatable :: settyp_(:)
      INTEGER, allocatable :: ielm_  (:)

      INTEGER :: numurban_
      INTEGER, allocatable :: urbclass (:)

      CHARACTER(len=256) :: suffix, cyear

      IF (p_is_master) THEN
         write(*,'(A)') 'Making urban type tiles :'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ! allocate and read the grided LCZ/NCAR urban type
      if (p_is_io) then

         dir_urban = trim(DEF_dir_rawdata) // '/urban_type'

         CALL allocate_block_data (gurban, data_urb_class)
         CALL flush_block_data (data_urb_class, 0)

         !write(cyear,'(i4.4)') int(lc_year/5)*5
         suffix = 'URBTYP'
IF (DEF_URBAN_type_scheme == 1) THEN
         ! NOTE!!!
         ! region id is assigned in aggreagation_urban.F90 now
         CALL read_5x5_data (dir_urban, suffix, gurban, 'URBAN_DENSITY_CLASS', data_urb_class)
ELSE IF (DEF_URBAN_type_scheme == 2) THEN
         CALL read_5x5_data (dir_urban, suffix, gurban, 'LCZ_DOM', data_urb_class)
ENDIF

#ifdef USEMPI
         CALL aggregation_data_daemon (gurban, data_i4_2d_in1 = data_urb_class)
#endif
      end if

      if (p_is_worker) then

         IF (numpatch > 0) THEN
            ! a temporary numpatch with max urban patch
            numpatch_ = numpatch + count(landpatch%settyp == URBAN) * (N_URB-1)

            allocate (eindex_(numpatch_))
            allocate (ipxstt_(numpatch_))
            allocate (ipxend_(numpatch_))
            allocate (settyp_(numpatch_))
            allocate (ielm_  (numpatch_))

            ! max urban patch number
            numurban_ = count(landpatch%settyp == URBAN) * N_URB
            IF (numurban_ > 0) THEN
               allocate (urbclass(numurban_))
            ENDIF
         ENDIF

         jpatch = 0
         iurban = 0

         ! loop for temporary numpatch to filter duplicate urban patch
         DO ipatch = 1, numpatch
            IF (landpatch%settyp(ipatch) == URBAN) THEN

               !???
               ie     = landpatch%ielm  (ipatch)
               ipxstt = landpatch%ipxstt(ipatch)
               ipxend = landpatch%ipxend(ipatch)

               CALL aggregation_request_data (landpatch, ipatch, gurban, &
                  data_i4_2d_in1 = data_urb_class, data_i4_2d_out1 = ibuff)

IF (DEF_URBAN_type_scheme == 1) THEN
               ! Some urban patches and NCAR data are inconsistent (NCAR has no urban ID),
               ! so the these points are assigned by the 3(medium density), or can define by ueser
               where (ibuff < 1 .or. ibuff > 3)
                  ibuff = 3
               END where
ELSE IF(DEF_URBAN_type_scheme == 2) THEN
               ! Same for NCAR, fill the gap LCZ class of urban patch if LCZ data is non-urban
               where (ibuff > 10 .or. ibuff == 0)
                  ibuff = 9
               END where
ENDIF

               npxl = ipxend - ipxstt + 1

               allocate (types (ipxstt:ipxend))

               types(:) = ibuff

               deallocate (ibuff)

               allocate (order (ipxstt:ipxend))
               order = (/ (ipxl, ipxl = ipxstt, ipxend) /)

               ! change order vars, types->regid
               ! add region information, because urban type may be same,
               ! but from different region in this urban patch
               ! relative code is changed
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
                  settyp_(jpatch) = URBAN
                  ipxstt_(jpatch) = ipxl
                  ielm_  (jpatch) = ie

                  iurban = iurban + 1
                  urbclass(iurban) = types(ipxl)
               ENDDO

               ipxend_(jpatch) = ipxend

               deallocate (types)
               deallocate (order)

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

         IF (numpatch > 0) THEN
            ! update landpath with new patch number
            ! all urban type patch are included
            IF (allocated (landpatch%eindex)) deallocate (landpatch%eindex)
            IF (allocated (landpatch%ipxstt)) deallocate (landpatch%ipxstt)
            IF (allocated (landpatch%ipxend)) deallocate (landpatch%ipxend)
            IF (allocated (landpatch%settyp)) deallocate (landpatch%settyp)
            IF (allocated (landpatch%ielm  )) deallocate (landpatch%ielm  )

            allocate (landpatch%eindex (numpatch))
            allocate (landpatch%ipxstt (numpatch))
            allocate (landpatch%ipxend (numpatch))
            allocate (landpatch%settyp (numpatch))
            allocate (landpatch%ielm   (numpatch))

            ! update all information of landpatch
            landpatch%eindex = eindex_(1:jpatch)
            landpatch%ipxstt = ipxstt_(1:jpatch)
            landpatch%ipxend = ipxend_(1:jpatch)
            landpatch%settyp = settyp_(1:jpatch)
            landpatch%ielm   = ielm_  (1:jpatch)
         ENDIF

         ! update urban patch number
         IF (numpatch > 0) THEN
            numurban = count(landpatch%settyp == URBAN)
         ELSE
            numurban = 0
         ENDIF

         IF (numurban > 0) THEN
            allocate (landurban%eindex (numurban))
            allocate (landurban%settyp (numurban))
            allocate (landurban%ipxstt (numurban))
            allocate (landurban%ipxend (numurban))
            allocate (landurban%ielm   (numurban))

            ! copy urban path information from landpatch for landurban
            landurban%eindex = pack(landpatch%eindex, landpatch%settyp == URBAN)
            landurban%ipxstt = pack(landpatch%ipxstt, landpatch%settyp == URBAN)
            landurban%ipxend = pack(landpatch%ipxend, landpatch%settyp == URBAN)
            landurban%ielm   = pack(landpatch%ielm  , landpatch%settyp == URBAN)

            ! assign urban region id and type for each urban patch
            landurban%settyp = urbclass(1:numurban)
         ENDIF

         ! update land patch with urban type patch
         ! set numurban
         landurban%nset = numurban
         landpatch%nset = numpatch
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

#if (defined CROP)
      IF (p_is_io) THEN
         !file_patch = trim(DEF_dir_rawdata) // '/global_0.5x0.5.MOD2005_V4.5_CFT_mergetoclmpft.nc'
         file_patch = trim(DEF_dir_rawdata) // '/global_0.5x0.5.MOD2005_V4.5_CFT_lf-merged-20220930.nc'
         CALL allocate_block_data (gcrop, cropdata, N_CFT)
         CALL ncio_read_block (file_patch, 'PCT_CFT', gcrop, N_CFT, cropdata)
      ENDIF

      cropfilter = (/ CROPLAND /)

      CALL pixelsetshadow_build (landpatch, gcrop, cropdata, N_CFT, cropfilter, &
         pctcrop, cropclass)

      numpatch = landpatch%nset
#endif

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numpatch, npatch_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npatch_glb, ' patches.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numpatch, ' patches.'
#endif

#if (defined CROP)
      CALL elm_patch%build (landelm, landpatch, use_frac = .true., shadowfrac = pctcrop)
#else
      CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
#endif

#ifdef CATCHMENT
#if (defined CROP)
      CALL hru_patch%build (landhru, landpatch, use_frac = .true., shadowfrac = pctcrop)
#else
      CALL hru_patch%build (landhru, landpatch, use_frac = .true.)
#endif
#endif

      CALL write_patchfrac (DEF_dir_landdata, lc_year)

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

      USE MOD_SPMD_Task
      USE MOD_LandPatch
      IMPLICIT NONE

      INTEGER :: ipatch, iurban

      IF (p_is_worker) THEN

         IF ((numpatch <= 0) .or. (numurban <= 0)) return

         IF (allocated(patch2urban)) deallocate(patch2urban)
         IF (allocated(urban2patch)) deallocate(urban2patch)
         allocate (patch2urban (numpatch))
         allocate (urban2patch (numurban))

         iurban = 0
         DO ipatch = 1, numpatch
            IF (landpatch%settyp(ipatch) == URBAN) THEN
               iurban = iurban + 1
               patch2urban(ipatch) = iurban
               urban2patch(iurban) = ipatch
            ELSE
               patch2urban(ipatch) = -1
            ENDIF
         ENDDO

      ENDIF

   END SUBROUTINE map_patch_to_urban

END MODULE MOD_LandUrban
