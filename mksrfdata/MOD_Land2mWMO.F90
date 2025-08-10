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
   integer :: spatch, epatch, ipatch
   integer :: ipth, numpxl, numpth
   integer :: src_pth, pthtype, maxpxl
   integer*8, allocatable :: eindex_(:)
   integer,   allocatable :: settyp_(:), ipxstt_(:), ipxend_(:), ielm_(:)
   integer,   allocatable :: wmopth (:), locpth(:)

   integer :: npatch_glb
   integer :: numwmo

      write(cyear,'(i4.4)') lc_year
      IF (p_is_master) THEN
         write(*,'(A)') 'Making land 2 m wmo patches :'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_worker) THEN

         numset = numelm

         allocate (wmo_source (numset))
         allocate (wmopth     (numset))

         wmo_source = -1
         wmopth     = -1
         numwmo     = 0

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
            pthtype = patchtypes(landpatch%settyp(spatch))
            IF (spatch /= epatch) THEN
               DO ipatch = spatch, epatch
                  numpxl = landpatch%ipxend(ipatch) - landpatch%ipxstt(ipatch) + 1

                  IF (numpxl>maxpxl .and. pthtype==0) THEN
                     maxpxl  = numpth
                     src_pth = ipatch
                  ENDIF
               ENDDO
            ELSE
               IF (pthtype == 0) THEN
                  src_pth = spatch
               ELSE
                  src_pth = -1
               ENDIF
            ENDIF

            IF (src_pth /= -1) numwmo = numwmo + 1
            wmo_source (iset) = src_pth
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

         ! set for new 2 m WMO patch
         DO iset = 1, numset
            numpth = count(landpatch%eindex==landelm%eindex(iset))

            IF (allocated(locpth)) deallocate(locpth)
            allocate(locpth(numpth))

            locpth = pack([(ipth, ipth=1, numpatch)], &
                     landpatch%eindex==landelm%eindex(iset))

            spatch = minval(locpth) ! elm_patch%substt(iset)
            epatch = maxval(locpth) ! elm_patch%subend(iset)

            !TODO@Wenzong: there may be problem with ipatch index?
            DO ipatch = spatch, epatch
               eindex_(ipatch) = landpatch%eindex(ipatch)
               settyp_(ipatch) = landpatch%settyp(ipatch)
               ipxstt_(ipatch) = landpatch%ipxstt(ipatch)
               ipxend_(ipatch) = landpatch%ipxend(ipatch)
               ielm_  (ipatch) = landpatch%ielm  (ipatch)
            ENDDO

            !TODO@Wenzong: there may be problem, set 2m wmo patch all the time?
            IF (wmo_source(iset) > 0) THEN
               eindex_(epatch+1) = landpatch%eindex(epatch)
               settyp_(epatch+1) = landpatch%settyp(wmo_source(iset))
               ipxstt_(epatch+1) = -1
               ipxend_(epatch+1) = -1
               ielm_  (epatch+1) = landpatch%ielm  (epatch)
               wmopth (iset)     = epatch+1
            ENDIF
         ENDDO

         ! allocate and save the new patches info
         IF (numpatch > 0) THEN
            ! update landpath with new patch number
            IF (allocated (landpatch%eindex)) deallocate (landpatch%eindex)
            IF (allocated (landpatch%ipxstt)) deallocate (landpatch%ipxstt)
            IF (allocated (landpatch%ipxend)) deallocate (landpatch%ipxend)
            IF (allocated (landpatch%settyp)) deallocate (landpatch%settyp)
            IF (allocated (landpatch%ielm  )) deallocate (landpatch%ielm  )

            allocate (landpatch%eindex (numpatch_))
            allocate (landpatch%ipxstt (numpatch_))
            allocate (landpatch%ipxend (numpatch_))
            allocate (landpatch%settyp (numpatch_))
            allocate (landpatch%ielm   (numpatch_))

            ! update all information of landpatch
            landpatch%eindex = eindex_(1:numpatch_)
            landpatch%ipxstt = ipxstt_(1:numpatch_)
            landpatch%ipxend = ipxend_(1:numpatch_)
            landpatch%settyp = settyp_(1:numpatch_)
            landpatch%ielm   = ielm_  (1:numpatch_)

            landelm%wmopth   = wmopth (1:numset   )
         ENDIF
      ENDIF

      landpatch%nset = numpatch_

      CALL landpatch%set_vecgs

#if (!defined(URBAN_MODEL) && !defined(CROP))
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
#endif
   END SUBROUTINE land2mwmo_build

END MODULE MOD_Land2mWMO
