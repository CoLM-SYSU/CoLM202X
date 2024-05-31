#include <define.h>

MODULE MOD_LandUrban
!-----------------------------------------------------------------------
!
! !DESCRIPTION:
!
!  Build pixelset "landurban".
!
!  Original authors: Hua Yuan and Wenzong Dong, 2021, OpenMP version.
!
!
! !REVISIONS:
!
!  05/2023, Wenzong Dong, Hua Yuan, Shupeng Zhang: porting codes to MPI
!           parallel version.
!
!-----------------------------------------------------------------------

   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_Vars_Global, only: N_URB, URBAN
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

   IMPLICIT NONE

   ! ---- Instance ----
   type(grid_type) :: gurban

   integer :: numurban
   type(pixelset_type) :: landurban

   integer , allocatable :: urban_reg   (:)  !region index of a urban
   integer , allocatable :: urban2patch (:)  !patch index of a urban
   integer , allocatable :: patch2urban (:)  !urban index of a patch

   ! ---- PUBLIC routines ----
   PUBLIC :: landurban_build
   PUBLIC :: map_patch_to_urban

CONTAINS

   ! -------------------------------
   SUBROUTINE landurban_build (lc_year)

   USE MOD_Precision
   USE MOD_Vars_Global
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
   USE MOD_AggregationRequestData
   USE MOD_Utils

   IMPLICIT NONE

   integer, intent(in) :: lc_year
   ! Local Variables
   character(len=256) :: dir_urban
   type (block_data_int32_2d) :: data_urb_class ! urban type index

   ! local vars
   integer, allocatable :: ibuff(:), types(:), order(:)

   ! index
   integer :: ipatch, jpatch, iurban
   integer :: ie, ipxstt, ipxend, npxl, ipxl
   integer :: nurb_glb, npatch_glb

   ! local vars for landpath and landurban
   integer :: numpatch_
   integer*8, allocatable :: eindex_(:)
   integer,   allocatable :: ipxstt_(:)
   integer,   allocatable :: ipxend_(:)
   integer,   allocatable :: settyp_(:)
   integer,   allocatable :: ielm_  (:)

   integer  :: numurban_
   integer  :: iurb, ib, imiss
   integer  :: buff_count(N_URB)
   real(r8) :: buff_p(N_URB)

   integer , allocatable :: urbclass (:)
   real(r8), allocatable :: area_one (:)

   character(len=256) :: suffix, cyear

      IF (p_is_master) THEN
         write(*,'(A)') 'Making urban type tiles :'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ! allocate and read the grided LCZ/NCAR urban type
      IF (p_is_io) THEN

         dir_urban = trim(DEF_dir_rawdata) // '/urban_type'

         CALL allocate_block_data (gurban, data_urb_class)
         CALL flush_block_data (data_urb_class, 0)

         ! read urban type data
         suffix = 'URBTYP'
IF (DEF_URBAN_type_scheme == 1) THEN
         CALL read_5x5_data (dir_urban, suffix, gurban, 'URBAN_DENSITY_CLASS', data_urb_class)
ELSE IF (DEF_URBAN_type_scheme == 2) THEN
         CALL read_5x5_data (dir_urban, suffix, gurban, 'LCZ_DOM', data_urb_class)
ENDIF

#ifdef USEMPI
         CALL aggregation_data_daemon (gurban, data_i4_2d_in1 = data_urb_class)
#endif
      ENDIF

      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN
            ! a temporary numpatch with max urban patch number
            numpatch_ = numpatch + count(landpatch%settyp == URBAN) * (N_URB-1)

            allocate (eindex_ (numpatch_ ))
            allocate (ipxstt_ (numpatch_ ))
            allocate (ipxend_ (numpatch_ ))
            allocate (settyp_ (numpatch_ ))
            allocate (ielm_   (numpatch_ ))

            ! max urban patch number (temporary)
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

               ie     = landpatch%ielm  (ipatch)
               ipxstt = landpatch%ipxstt(ipatch)
               ipxend = landpatch%ipxend(ipatch)

               CALL aggregation_request_data (landpatch, ipatch, gurban, zip = .false., area = area_one, &
                  data_i4_2d_in1 = data_urb_class, data_i4_2d_out1 = ibuff)

               ! when there is missing urban types
               !NOTE@tungwz: need duoble check below and add appropriate annotations
               ! check if there is urban pixel without URBAN ID
               imiss = count(ibuff<1 .or. ibuff>N_URB)
               IF (imiss > 0) THEN
                  ! Calculate the relative ratio of each urban types by excluding urban pixels withoht URBAN ID
                  WHERE (ibuff<1 .or. ibuff>N_URB)
                     area_one = 0
                  END WHERE

                  buff_p = 0
                  IF (sum(area_one) > 0) THEN
                     DO ib = 1, size(area_one)
                        IF (ibuff(ib)>1 .and. ibuff(ib)<N_URB) THEN
                           iurb         = ibuff(ib)
                           buff_p(iurb) = buff_p(iurb) + area_one(ib)
                        ENDIF
                     ENDDO
                     buff_p(:) = buff_p(:)/sum(area_one)
                  ENDIF

                  ! The number of URBAN ID of each type is assigned to urban pixels without URBAN ID in relative proportion
                  DO iurb = 1, N_URB-1
                     buff_count(iurb) = int(buff_p(iurb)*imiss)
                  ENDDO
                  buff_count(N_URB) = imiss - sum(buff_count(1:N_URB-1))

                  ! Some urban patches and NCAR/LCZ data are inconsistent (NCAR/LCZ has no urban ID),
                  ! so the these points are assigned
                  IF (all(buff_count==0)) THEN
                     ! If none of the urban pixels have an URBAN ID, they are assigned directly
                     IF (DEF_URBAN_type_scheme == 1) THEN
                        ibuff = 3
                     ELSEIF (DEF_URBAN_type_scheme == 2) THEN
                        ibuff = 9
                     ENDIF
                  ELSE
                     ! Otherwise, URBAN ID are assigned based on the previously calculated number
                     DO ib = 1, size(ibuff)
                        IF (ibuff(ib)<1 .or. ibuff(ib)>N_URB) THEN
                           type_loop: DO iurb = 1, N_URB
                              IF (buff_count(iurb) > 0) THEN
                                 ibuff(ib)        = iurb
                                 buff_count(iurb) = buff_count(iurb) - 1
                                 EXIT type_loop
                              ENDIF
                           ENDDO type_loop
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF

               npxl = ipxend - ipxstt + 1

               allocate (types (ipxstt:ipxend))

               types(:) = ibuff

               deallocate (ibuff)

               allocate (order (ipxstt:ipxend))
               order = (/ (ipxl, ipxl = ipxstt, ipxend) /)

               ! change order vars, types->regid ? still types below
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
                        CYCLE
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

#ifdef SinglePoint

      allocate  ( SITE_urbtyp    (numurban) )
      allocate  ( SITE_lucyid    (numurban) )

IF (.not. USE_SITE_urban_paras) THEN
      allocate  ( SITE_fveg_urb  (numurban) )
      allocate  ( SITE_htop_urb  (numurban) )
      allocate  ( SITE_flake_urb (numurban) )

      allocate  ( SITE_popden    (numurban) )
      allocate  ( SITE_froof     (numurban) )
      allocate  ( SITE_hroof     (numurban) )
      allocate  ( SITE_hwr       (numurban) )
      allocate  ( SITE_fgper     (numurban) )
      allocate  ( SITE_fgimp     (numurban) )
ENDIF

      allocate  ( SITE_em_roof   (numurban) )
      allocate  ( SITE_em_wall   (numurban) )
      allocate  ( SITE_em_gimp   (numurban) )
      allocate  ( SITE_em_gper   (numurban) )
      allocate  ( SITE_t_roommax (numurban) )
      allocate  ( SITE_t_roommin (numurban) )
      allocate  ( SITE_thickroof (numurban) )
      allocate  ( SITE_thickwall (numurban) )

      allocate  ( SITE_cv_roof   (nl_roof)  )
      allocate  ( SITE_cv_wall   (nl_wall)  )
      allocate  ( SITE_cv_gimp   (nl_soil)  )
      allocate  ( SITE_tk_roof   (nl_roof)  )
      allocate  ( SITE_tk_wall   (nl_wall)  )
      allocate  ( SITE_tk_gimp   (nl_soil)  )

      allocate  ( SITE_alb_roof  (2, 2)     )
      allocate  ( SITE_alb_wall  (2, 2)     )
      allocate  ( SITE_alb_gimp  (2, 2)     )
      allocate  ( SITE_alb_gper  (2, 2)     )

      SITE_urbtyp(:) = landurban%settyp
#endif

#ifndef CROP
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

      CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
#ifdef CATCHMENT
      CALL hru_patch%build (landhru, landpatch, use_frac = .true.)
#endif
      CALL write_patchfrac (DEF_dir_landdata, lc_year)
#endif

      IF (allocated (ibuff   )) deallocate (ibuff    )
      IF (allocated (types   )) deallocate (types    )
      IF (allocated (order   )) deallocate (order    )

      IF (allocated (eindex_ )) deallocate (eindex_  )
      IF (allocated (ipxstt_ )) deallocate (ipxstt_  )
      IF (allocated (ipxend_ )) deallocate (ipxend_  )
      IF (allocated (settyp_ )) deallocate (settyp_  )
      IF (allocated (ielm_   )) deallocate (ielm_    )

      IF (allocated (urbclass)) deallocate (urbclass )
      IF (allocated (area_one)) deallocate (area_one )

   END SUBROUTINE landurban_build

   ! ----------------------
   SUBROUTINE map_patch_to_urban

   USE MOD_SPMD_Task
   USE MOD_LandPatch
   IMPLICIT NONE

   integer :: ipatch, iurban

      IF (p_is_worker) THEN

         IF ((numpatch <= 0) .or. (numurban <= 0)) RETURN

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
