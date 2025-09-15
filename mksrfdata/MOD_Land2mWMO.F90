#include <define.h>

MODULE MOD_Land2mWMO

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Build a virtual patch "land2mWMO" for output 2 m WMO temperature.
!
!  Created by Wenzong Dong and Hua Yuan, 2025/08
!
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_LandPatch
   USE MOD_LandElm
   USE MOD_Vars_Global
   USE MOD_Const_LC
   IMPLICIT NONE

   ! ---- Instance ----
   integer, allocatable :: wmo_patch  (:)  !2m wmo patch index
   integer, allocatable :: wmo_source (:)  !source patch of a wmo patch

CONTAINS

   ! -------------------------------
   SUBROUTINE land2mwmo_build (lc_year)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Utils
   USE MOD_UserDefFun
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_Namelist
   USE MOD_NetCDFBlock

   IMPLICIT NONE

   integer, intent(in) :: lc_year
   character(len=255) :: cyear

   integer :: numpatch_
   integer :: numset
   integer :: iset
   integer :: spatch, epatch, ipatch, jpatch
   integer :: ipth, numpxl, numpth
   integer :: src_pth, pthtype, maxpxl
   integer*8, allocatable :: eindex_(:)
   integer,   allocatable :: settyp_(:), ipxstt_(:), ipxend_(:), ielm_(:)
   integer,   allocatable :: locpth(:)

   integer :: npatch_glb
   integer :: numwmo

      write(cyear,'(i4.4)') lc_year
      IF (p_is_master) THEN
         write(*,'(A)') 'Making land 2 m wmo patches:'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_worker) THEN

         numset = numelm

         numwmo     = 0
         numpatch_  = 0
         jpatch     = 0

         ! Count for 2 m WMO patches need to be set virtually
         DO iset = 1, numset

            numpth = count(landpatch%eindex==landelm%eindex(iset))
            IF (allocated(locpth)) deallocate(locpth)
            allocate(locpth(numpth))
            locpth = pack([(ipth, ipth=1, numpatch)], &
                     landpatch%eindex==landelm%eindex(iset))

            spatch = minval(locpth) ! elm_patch%substt(iset)
            epatch = maxval(locpth) ! elm_patch%subend(iset)

            maxpxl  = 0
            src_pth = -1

            DO ipatch = spatch, epatch

               pthtype = patchtypes(landpatch%settyp(ipatch))
               numpxl  = landpatch%ipxend(ipatch) - landpatch%ipxstt(ipatch) + 1

               IF (numpxl>maxpxl .and. pthtype==0) THEN
                  maxpxl  = numpxl
                  src_pth = ipatch
               ENDIF
            ENDDO

            ! a new 2m WMO patch
            IF (src_pth /= -1) THEN
               wmo_source (iset) = src_pth
               numwmo = numwmo + 1
               landelm%settyp(iset) = 1
            ELSE
               wmo_source (iset) = -1
            ENDIF

         ENDDO

         ! allocate new temporal patches memory
         IF (numpatch > 0) THEN
            ! a numpatch with WMO patch number
            numpatch_ = numpatch + numwmo

            allocate (eindex_ (numpatch_ ))
            allocate (ipxstt_ (numpatch_ ))
            allocate (ipxend_ (numpatch_ ))
            allocate (settyp_ (numpatch_ ))
            allocate (ielm_   (numpatch_ ))

         ENDIF

         numwmo = 0

         ! set for new 2 m WMO patch
         DO iset = 1, numset
            numpth = count(landpatch%eindex==landelm%eindex(iset))

            IF (allocated(locpth)) deallocate(locpth)
            allocate(locpth(numpth))

            locpth = pack([(ipth, ipth=1, numpatch)], &
                     landpatch%eindex==landelm%eindex(iset))

            spatch = minval(locpth) ! elm_patch%substt(iset)
            epatch = maxval(locpth) ! elm_patch%subend(iset)

            DO ipatch = spatch, epatch
               jpatch = jpatch + 1
               eindex_(jpatch) = landpatch%eindex(ipatch)
               settyp_(jpatch) = landpatch%settyp(ipatch)
               ipxstt_(jpatch) = landpatch%ipxstt(ipatch)
               ipxend_(jpatch) = landpatch%ipxend(ipatch)
               ielm_  (jpatch) = landpatch%ielm  (ipatch)
            ENDDO

            IF (wmo_source(iset) > 0) THEN
               jpatch = jpatch + 1
               eindex_(jpatch) = landpatch%eindex(epatch)
               settyp_(jpatch) = landpatch%settyp(wmo_source(iset))
               ipxstt_(jpatch) = -1
               ipxend_(jpatch) = -1
               ielm_  (jpatch) = landpatch%ielm  (epatch)

               ! update the newly 2m WMO source/patch index
               wmo_patch (iset) = jpatch
               wmo_source(iset) = wmo_source(iset) + numwmo
               numwmo           = numwmo + 1
            ENDIF
         ENDDO

         ! 2m WMO patch number check
         IF (jpatch .ne. numpatch_) THEN
            write(*,'(A)') 'Count land 2 m WMO patches error! See MOD_Land2mWMO.F90.'
         ENDIF

         ! set the new patch number
         numpatch = numpatch_

         ! allocate and save the new patches info
         IF (numpatch > 0) THEN
            ! update landpath with new patch number
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
            landpatch%eindex = eindex_(1:numpatch)
            landpatch%ipxstt = ipxstt_(1:numpatch)
            landpatch%ipxend = ipxend_(1:numpatch)
            landpatch%settyp = settyp_(1:numpatch)
            landpatch%ielm   = ielm_  (1:numpatch)

         ENDIF
      ENDIF

      landpatch%nset = numpatch

      CALL landpatch%set_vecgs

      IF (allocated (eindex_ )) deallocate (eindex_ )
      IF (allocated (settyp_ )) deallocate (settyp_ )
      IF (allocated (ipxstt_ )) deallocate (ipxstt_ )
      IF (allocated (ipxend_ )) deallocate (ipxend_ )
      IF (allocated (ielm_   )) deallocate (ielm_   )
      IF (allocated (locpth  )) deallocate (locpth  )

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

      CALL write_patchfrac (DEF_dir_landdata, lc_year)

   END SUBROUTINE land2mwmo_build

   SUBROUTINE land2mwmo_init
   USE MOD_Mesh
   IMPLICIT NONE

      allocate (wmo_patch  (numelm))
      allocate (wmo_source (numelm))

      wmo_patch  = -1
      wmo_source = -1

   END SUBROUTINE land2mwmo_init

   SUBROUTINE land2mwmo_final
   IMPLICIT NONE

      IF (allocated (wmo_patch )) deallocate (wmo_patch )
      IF (allocated (wmo_source)) deallocate (wmo_source)

   END SUBROUTINE land2mwmo_final

END MODULE MOD_Land2mWMO
