#include <define.h>

MODULE MOD_Pixelset

!------------------------------------------------------------------------------------
! !DESCRIPTION:
!
!    Pixelset refers to a set of pixels in CoLM.
!
!    In CoLM, the global/regional area is divided into a hierarchical structure:
!    1. If GRIDBASED or UNSTRUCTURED is defined, it is
!       ELEMENT >>> PATCH
!    2. If CATCHMENT is defined, it is
!       ELEMENT >>> HRU >>> PATCH
!    If Plant FUNCTION Type classification is used, PATCH is further divided into PFT.
!    If Plant Community classification is used,     PATCH is further divided into PC.
!
!    In CoLM, the land surface is first divided into pixels, which are rasterized
!    points defined by fine-resolution data. Then ELEMENT, PATCH, HRU, PFT, PC
!    are all consists of pixels, and hence they are all pixelsets.
!
!    The highest level pixelset in CoLM is ELEMENT, all other pixelsets are subsets
!    of ELEMENTs.
!    In a pixelset, pixels are sorted to make pixels in its subsets consecutive.
!    Thus a subset can be represented by starting pixel index and ending pixel index
!    in an ELEMENT.
!
!                Example of hierarchical pixelsets
!        ************************************************ <-- pixels in an ELEMENT
!        |<------------------- ELEMENT ---------------->| <-- level 1
!        |   subset 1  |       subset 2      | subset 3 | <-- level 2
!        |s11|   s12   | s21 |   s22   | s23 |    s31   | <-- level 3
!
!    "Vector" is a collection of data when each pixelset in a given level is associated
!    with a value, representing its averaged physical, chemical or biological state.
!
!    "Vector" is usually defined on worker process, while its IO is through IO process.
!    To read,  vector is first loaded from files by IO and then scattered from IO to worker.
!    To write, vector is first gathered from worker to IO and then saved to files by IO.
!
!  Created by Shupeng Zhang, May 2023
!------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_DataType
   IMPLICIT NONE

   ! ---- data types ----
   type :: vec_gather_scatter_type

      ! for worker and io
      integer, allocatable :: vlen(:,:)

      ! for worker
      integer, allocatable :: vstt(:,:)
      integer, allocatable :: vend(:,:)

      ! for io
      integer, allocatable :: vcnt(:,:,:)
      integer, allocatable :: vdsp(:,:,:)

   CONTAINS
      final  :: vec_gather_scatter_free_mem

   END type vec_gather_scatter_type

   ! ---- data types ----
   type :: pixelset_type

      integer :: nset

      integer*8, allocatable :: eindex(:)  ! global index of element to which pixelset belongs

      integer, allocatable :: ipxstt(:)    ! start local index of pixel in the element
      integer, allocatable :: ipxend(:)    ! end   local index of pixel in the element
      integer, allocatable :: settyp(:)    ! type of pixelset

      integer, allocatable :: ielm(:)      ! local index of element to which pixelset belongs

      integer :: nblkgrp                   ! number of blocks for this process's working group
      integer, allocatable :: xblkgrp (:)  ! block index in longitude for this process's group
      integer, allocatable :: yblkgrp (:)  ! block index in latitude  for this process's group

      type(vec_gather_scatter_type) :: vecgs ! for vector gathering and scattering

      logical :: has_shared = .false.
      real(r8), allocatable :: pctshared (:)

   CONTAINS
      procedure, PUBLIC :: set_vecgs         => vec_gather_scatter_set
      procedure, PUBLIC :: get_lonlat_radian => pixelset_get_lonlat_radian
      procedure, PUBLIC :: pset_pack         => pixelset_pack
      procedure, PUBLIC :: forc_free_mem     => pixelset_forc_free_mem
      final :: pixelset_free_mem

   END type pixelset_type

   ! ---- data types ----
   type :: subset_type

      integer,  allocatable :: substt(:)
      integer,  allocatable :: subend(:)
      real(r8), allocatable :: subfrc(:)

   CONTAINS
      procedure, PUBLIC :: build => subset_build
      final :: subset_free_mem

   END type subset_type

   ! ---- data types ----
   type :: superset_type

      integer,  allocatable :: sup(:)

   CONTAINS
      procedure, PUBLIC :: build => superset_build
      final :: superset_free_mem

   END type superset_type

CONTAINS

   ! --------------------------------
   SUBROUTINE pixelset_get_lonlat_radian (this, rlon, rlat)

   USE MOD_Precision
   USE MOD_Utils
   USE MOD_Pixel
   USE MOD_Mesh

   IMPLICIT NONE
   CLASS(pixelset_type) :: this

   real(r8), intent(inout) :: rlon(:), rlat(:)

   ! Local Variables
   integer :: iset, ie, ipxstt, ipxend, npxl, ipxl
   real(r8), allocatable :: area(:)

      DO iset = 1, this%nset

         ie = this%ielm(iset)

         ipxstt = this%ipxstt (iset)
         ipxend = this%ipxend (iset)

         ! for 2m WMO patch, use all pixels
         IF (ipxstt == -1) THEN
            ipxstt = 1
            ipxend = mesh(ie)%npxl
         ENDIF

         allocate (area (ipxstt:ipxend))
         DO ipxl = ipxstt, ipxend
            area(ipxl) = areaquad (&
               pixel%lat_s(mesh(ie)%ilat(ipxl)), &
               pixel%lat_n(mesh(ie)%ilat(ipxl)), &
               pixel%lon_w(mesh(ie)%ilon(ipxl)), &
               pixel%lon_e(mesh(ie)%ilon(ipxl)) )
         ENDDO

         npxl = ipxend - ipxstt + 1
         rlat(iset) = get_pixelset_rlat ( &
            npxl, mesh(ie)%ilat(ipxstt:ipxend), area)
         rlon(iset) = get_pixelset_rlon ( &
            npxl, mesh(ie)%ilon(ipxstt:ipxend), area)

         deallocate (area)

      ENDDO

   END SUBROUTINE pixelset_get_lonlat_radian

   ! --------------------------------
   FUNCTION get_pixelset_rlat (npxl, ilat, area) result(rlat)

   USE MOD_Precision
   USE MOD_Vars_Global, only: pi
   USE MOD_Pixel
   IMPLICIT NONE

   real(r8) :: rlat

   integer,  intent(in) :: npxl
   integer,  intent(in) :: ilat(npxl)
   real(r8), intent(in) :: area(npxl)

   ! Local variables
   integer :: ipxl

      rlat = 0.0
      DO ipxl = 1, npxl
         rlat = rlat + (pixel%lat_s(ilat(ipxl)) + pixel%lat_n(ilat(ipxl))) * 0.5 * area(ipxl)
      ENDDO
      rlat = rlat / sum(area) * pi/180.0

   END FUNCTION get_pixelset_rlat

   ! --------------------------------
   FUNCTION get_pixelset_rlon (npxl, ilon, area) result(rlon)

   USE MOD_Precision
   USE MOD_Utils
   USE MOD_Vars_Global, only: pi
   USE MOD_Pixel
   IMPLICIT NONE

   real(r8) :: rlon

   integer,  intent(in) :: npxl
   integer,  intent(in) :: ilon(npxl)
   real(r8), intent(in) :: area(npxl)

   ! Local variables
   integer  :: ipxl
   real(r8) :: lon, lon0, area_done

      lon = 0.0
      area_done = 0.0
      DO ipxl = 1, npxl

         IF (pixel%lon_w(ilon(ipxl)) > pixel%lon_e(ilon(ipxl))) THEN
            lon0 = (pixel%lon_w(ilon(ipxl)) + pixel%lon_e(ilon(ipxl)) + 360.0) * 0.5
         ELSE
            lon0 = (pixel%lon_w(ilon(ipxl)) + pixel%lon_e(ilon(ipxl))) * 0.5
         ENDIF

         CALL normalize_longitude (lon0)

         IF (lon - lon0 > 180._r8) THEN
            lon = lon * area_done + (lon0 + 360._r8) * area(ipxl)
         ELSEIF (lon - lon0 < -180._r8) THEN
            lon = lon * area_done + (lon0 - 360._r8) * area(ipxl)
         ELSE
            lon = lon * area_done + lon0 * area(ipxl)
         ENDIF

         area_done = area_done + area(ipxl)
         lon = lon / area_done

         CALL normalize_longitude(lon)

      ENDDO

      rlon = lon * pi/180.0

   END FUNCTION get_pixelset_rlon

   ! --------------------------------
   SUBROUTINE pixelset_free_mem (this)

   IMPLICIT NONE
   type (pixelset_type) :: this

      IF (allocated(this%eindex)) deallocate(this%eindex)
      IF (allocated(this%ipxstt)) deallocate(this%ipxstt)
      IF (allocated(this%ipxend)) deallocate(this%ipxend)
      IF (allocated(this%settyp)) deallocate(this%settyp)

      IF (allocated(this%ielm  )) deallocate(this%ielm  )

      IF (allocated(this%xblkgrp)) deallocate(this%xblkgrp)
      IF (allocated(this%yblkgrp)) deallocate(this%yblkgrp)

      IF (allocated(this%pctshared)) deallocate(this%pctshared)

   END SUBROUTINE pixelset_free_mem

   ! --------------------------------
   SUBROUTINE pixelset_forc_free_mem (this)

   IMPLICIT NONE

   class(pixelset_type) :: this

      IF (allocated(this%eindex )) deallocate(this%eindex )
      IF (allocated(this%ipxstt )) deallocate(this%ipxstt )
      IF (allocated(this%ipxend )) deallocate(this%ipxend )
      IF (allocated(this%settyp )) deallocate(this%settyp )

      IF (allocated(this%ielm   )) deallocate(this%ielm   )

      IF (allocated(this%xblkgrp)) deallocate(this%xblkgrp)
      IF (allocated(this%yblkgrp)) deallocate(this%yblkgrp)

      IF (allocated(this%pctshared)) deallocate(this%pctshared)

   END SUBROUTINE pixelset_forc_free_mem

   ! --------------------------------
   SUBROUTINE copy_pixelset(pixel_from, pixel_to)

   IMPLICIT NONE

   type(pixelset_type), intent(in)  :: pixel_from
   type(pixelset_type), intent(out) :: pixel_to

      pixel_to%nset    = pixel_from%nset
      pixel_to%eindex  = pixel_from%eindex
      pixel_to%ipxstt  = pixel_from%ipxstt
      pixel_to%ipxend  = pixel_from%ipxend
      pixel_to%settyp  = pixel_from%settyp
      pixel_to%ielm    = pixel_from%ielm

      pixel_to%nblkgrp = pixel_from%nblkgrp
      pixel_to%xblkgrp = pixel_from%xblkgrp
      pixel_to%yblkgrp = pixel_from%yblkgrp

      IF (pixel_from%has_shared) THEN
         pixel_to%pctshared = pixel_from%pctshared
      ENDIF

   END SUBROUTINE

   ! --------------------------------
   SUBROUTINE vec_gather_scatter_set (this)

   USE MOD_Block
   USE MOD_SPMD_Task
   USE MOD_Mesh
   IMPLICIT NONE

   class(pixelset_type)  :: this

   ! Local variables
   integer :: iproc
   integer :: iset, ie, xblk, yblk, iblk, jblk, scnt, iblkgrp, iblkall
   logical, allocatable :: nonzero(:,:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (.not. allocated (this%vecgs%vlen)) THEN
         allocate (this%vecgs%vlen (gblock%nxblk, gblock%nyblk))
         this%vecgs%vlen(:,:) = 0
      ENDIF

      IF (p_is_worker) THEN

         IF (.not. allocated (this%vecgs%vstt)) THEN
            allocate (this%vecgs%vstt (gblock%nxblk, gblock%nyblk))
            allocate (this%vecgs%vend (gblock%nxblk, gblock%nyblk))
         ENDIF

         this%vecgs%vstt(:,:) = 0
         this%vecgs%vend(:,:) = -1

         ie = 1
         xblk = 0
         yblk = 0
         DO iset = 1, this%nset
            DO WHILE (this%eindex(iset) /= mesh(ie)%indx)
               ie = ie + 1
            ENDDO

            IF ((mesh(ie)%xblk /= xblk) .or. (mesh(ie)%yblk /= yblk)) THEN
               xblk = mesh(ie)%xblk
               yblk = mesh(ie)%yblk
               this%vecgs%vstt(xblk,yblk) = iset
            ENDIF

            this%vecgs%vend(xblk,yblk) = iset
         ENDDO

         this%vecgs%vlen = this%vecgs%vend - this%vecgs%vstt + 1

#ifdef USEMPI
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN

                  scnt = this%vecgs%vlen(iblk,jblk)
                  CALL mpi_gather (scnt, 1, MPI_INTEGER, &
                     MPI_INULL_P, 1, MPI_INTEGER, p_root, p_comm_group, p_err)

               ENDIF
            ENDDO
         ENDDO
#endif
      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN

         IF (.not. allocated(this%vecgs%vcnt)) THEN
            allocate (this%vecgs%vcnt (0:p_np_group-1,gblock%nxblk,gblock%nyblk))
            allocate (this%vecgs%vdsp (0:p_np_group-1,gblock%nxblk,gblock%nyblk))
         ENDIF

         this%vecgs%vcnt(:,:,:) = 0
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                  scnt = 0
                  CALL mpi_gather (scnt, 1, MPI_INTEGER, &
                     this%vecgs%vcnt(:,iblk,jblk), 1, MPI_INTEGER, &
                     p_root, p_comm_group, p_err)

                  this%vecgs%vdsp(0,iblk,jblk) = 0
                  DO iproc = 1, p_np_group-1
                     this%vecgs%vdsp(iproc,iblk,jblk) = &
                        this%vecgs%vdsp(iproc-1,iblk,jblk) + this%vecgs%vcnt(iproc-1,iblk,jblk)
                  ENDDO

                  this%vecgs%vlen(iblk,jblk) = sum(this%vecgs%vcnt(:,iblk,jblk))

               ENDIF
            ENDDO
         ENDDO
      ENDIF
#endif

      IF (p_is_io .or. p_is_worker) THEN
         allocate (nonzero (gblock%nxblk,gblock%nyblk))

         nonzero = this%vecgs%vlen > 0
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, nonzero, gblock%nxblk * gblock%nyblk, &
            MPI_LOGICAL, MPI_LOR, p_comm_group, p_err)
#endif

         this%nblkgrp = count(nonzero)
         IF (allocated(this%xblkgrp)) deallocate(this%xblkgrp)
         IF (allocated(this%yblkgrp)) deallocate(this%yblkgrp)
         allocate (this%xblkgrp (this%nblkgrp))
         allocate (this%yblkgrp (this%nblkgrp))

         iblkgrp = 0
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (nonzero(iblk,jblk)) THEN
                  iblkgrp = iblkgrp + 1
                  this%xblkgrp(iblkgrp) = iblk
                  this%yblkgrp(iblkgrp) = jblk
               ENDIF
            ENDDO
         ENDDO

         deallocate(nonzero)
      ENDIF

   END SUBROUTINE vec_gather_scatter_set

   ! --------------------------------
   SUBROUTINE pixelset_pack (this, mask, nset_packed)

   USE MOD_SPMD_Task
   IMPLICIT NONE
   class(pixelset_type) :: this
   logical, intent(in)  :: mask(:)
   integer, intent(out) :: nset_packed

   integer*8, allocatable :: eindex_(:)
   integer,   allocatable :: ipxstt_(:)
   integer,   allocatable :: ipxend_(:)
   integer,   allocatable :: settyp_(:)
   integer,   allocatable :: ielm_  (:)

   real(r8),  allocatable :: pctshared_(:)
   integer :: s, e

      IF (p_is_worker) THEN

         IF (this%nset > 0) THEN
            IF (count(mask) < this%nset) THEN

               allocate (eindex_(this%nset))
               allocate (ipxstt_(this%nset))
               allocate (ipxend_(this%nset))
               allocate (settyp_(this%nset))
               allocate (ielm_  (this%nset))

               eindex_ = this%eindex
               ipxstt_ = this%ipxstt
               ipxend_ = this%ipxend
               settyp_ = this%settyp
               ielm_   = this%ielm

               deallocate (this%eindex)
               deallocate (this%ipxstt)
               deallocate (this%ipxend)
               deallocate (this%settyp)
               deallocate (this%ielm  )

               IF (this%has_shared) THEN
                  allocate   (pctshared_(this%nset))
                  pctshared_ = this%pctshared
                  deallocate (this%pctshared)
               ENDIF

               this%nset = count(mask)

               IF (this%nset > 0) THEN

                  allocate (this%eindex(this%nset))
                  allocate (this%ipxstt(this%nset))
                  allocate (this%ipxend(this%nset))
                  allocate (this%settyp(this%nset))
                  allocate (this%ielm  (this%nset))

                  this%eindex = pack(eindex_, mask)
                  this%ipxstt = pack(ipxstt_, mask)
                  this%ipxend = pack(ipxend_, mask)
                  this%settyp = pack(settyp_, mask)
                  this%ielm   = pack(ielm_  , mask)

                  IF (this%has_shared) THEN

                     this%pctshared = pack(pctshared_, mask)

                     s = 1
                     DO WHILE (s < this%nset)
                        e = s
                        DO WHILE (e < this%nset)
                           IF ((this%ielm(e+1) == this%ielm(s)) &
                              .and. (this%ipxstt(e+1) == this%ipxstt(s))) THEN
                              e = e + 1
                           ELSE
                              EXIT
                           ENDIF
                        ENDDO

                        IF (e > s) THEN
                           this%pctshared(s:e) = this%pctshared(s:e)/sum(this%pctshared(s:e))
                        ENDIF

                        s = e + 1
                     ENDDO

                  ENDIF

               ENDIF

               deallocate (eindex_)
               deallocate (ipxstt_)
               deallocate (ipxend_)
               deallocate (settyp_)
               deallocate (ielm_  )

               IF (this%has_shared) THEN
                  deallocate (pctshared_)
               ENDIF

            ENDIF
         ENDIF

      ENDIF

      CALL this%set_vecgs

      nset_packed = this%nset

   END SUBROUTINE pixelset_pack

   ! --------------------------------
   SUBROUTINE vec_gather_scatter_free_mem (this)

   IMPLICIT NONE
   type (vec_gather_scatter_type) :: this

      IF (allocated(this%vlen))  deallocate (this%vlen)
      IF (allocated(this%vstt))  deallocate (this%vstt)
      IF (allocated(this%vend))  deallocate (this%vend)
      IF (allocated(this%vcnt))  deallocate (this%vcnt)
      IF (allocated(this%vdsp))  deallocate (this%vdsp)

   END SUBROUTINE vec_gather_scatter_free_mem

   ! --------------------------------
   SUBROUTINE subset_build (this, superset, subset, use_frac)

   USE MOD_Mesh
   USE MOD_Pixel
   USE MOD_Utils
   IMPLICIT NONE

   CLASS(subset_type) :: this

   type (pixelset_type), intent(in) :: superset
   type (pixelset_type), intent(in) :: subset
   logical, intent(in) :: use_frac

   ! Local Variables
   integer :: isuperset, isubset, ielm, ipxl, istt, iend

      IF (superset%nset <= 0) RETURN

      IF (superset%has_shared) THEN
         write(*,*) 'Warning: superset has shared area.'
      ENDIF

      IF (allocated(this%substt)) deallocate(this%substt)
      IF (allocated(this%subend)) deallocate(this%subend)

      allocate (this%substt (superset%nset))
      allocate (this%subend (superset%nset))

      this%substt =  0
      this%subend = -1

      isuperset = 1
      isubset   = 1
      DO WHILE (isubset <= subset%nset)
         IF (     (subset%eindex(isubset) == superset%eindex(isuperset)) &
            .and. (subset%ipxstt(isubset) >= superset%ipxstt(isuperset) .or. &
                   subset%ipxstt(isubset) == -1 ) &
            .and. (subset%ipxend(isubset) <= superset%ipxend(isuperset) .or. &
                   subset%ipxend(isubset) == -1 ) ) THEN

            IF (this%substt(isuperset) == 0) THEN
               this%substt(isuperset) = isubset
            ENDIF

            this%subend(isuperset) = isubset

            isubset = isubset + 1
         ELSE
            isuperset = isuperset + 1
         ENDIF
      ENDDO

      IF (use_frac) THEN

         IF (allocated(this%subfrc)) deallocate(this%subfrc)

         IF (subset%nset <= 0) RETURN

         allocate (this%subfrc (subset%nset))

         DO isubset = 1, subset%nset
            ielm = subset%ielm(isubset)
            this%subfrc(isubset) = 0
            DO ipxl = subset%ipxstt(isubset), subset%ipxend(isubset)
               IF (ipxl == -1) CYCLE
               this%subfrc(isubset) = this%subfrc(isubset) &
                  + areaquad (&
                  pixel%lat_s(mesh(ielm)%ilat(ipxl)), &
                  pixel%lat_n(mesh(ielm)%ilat(ipxl)), &
                  pixel%lon_w(mesh(ielm)%ilon(ipxl)), &
                  pixel%lon_e(mesh(ielm)%ilon(ipxl)) )
            ENDDO
            IF (subset%has_shared) THEN
               this%subfrc(isubset) = this%subfrc(isubset) * subset%pctshared(isubset)
            ENDIF
         ENDDO

         DO isuperset = 1, superset%nset
            IF (this%substt(isuperset) /= 0) THEN
               istt = this%substt(isuperset)
               iend = this%subend(isuperset)
               this%subfrc(istt:iend) = this%subfrc(istt:iend) / sum(this%subfrc(istt:iend))
            ENDIF
         ENDDO

      ENDIF

   END SUBROUTINE subset_build

   ! --------------------------------
   SUBROUTINE subset_free_mem (this)

   IMPLICIT NONE
   type (subset_type) :: this

      IF (allocated(this%substt))  deallocate (this%substt)
      IF (allocated(this%subend))  deallocate (this%subend)
      IF (allocated(this%subfrc))  deallocate (this%subfrc)

   END SUBROUTINE subset_free_mem

   ! --------------------------------
   SUBROUTINE superset_build (this, superset, subset)

   IMPLICIT NONE

   CLASS(superset_type) :: this

   type (pixelset_type), intent(in) :: superset
   type (pixelset_type), intent(in) :: subset

   ! Local Variables
   integer :: isuperset, isubset

      IF (subset%nset <= 0) RETURN

      IF (allocated(this%sup)) deallocate(this%sup)

      allocate (this%sup (subset%nset))

      isuperset = 1
      isubset   = 1
      DO WHILE (isubset <= subset%nset)
         IF (     (subset%eindex(isubset) == superset%eindex(isuperset)) &
            .and. (subset%ipxstt(isubset) >= superset%ipxstt(isuperset)) &
            .and. (subset%ipxend(isubset) <= superset%ipxend(isuperset))) THEN

            this%sup(isubset) = isuperset

            isubset = isubset + 1
         ELSE
            isuperset = isuperset + 1
         ENDIF
      ENDDO

   END SUBROUTINE superset_build

   ! --------------------------------
   SUBROUTINE superset_free_mem (this)

   IMPLICIT NONE
   type (superset_type) :: this

      IF (allocated(this%sup))  deallocate (this%sup)

   END SUBROUTINE superset_free_mem

END MODULE MOD_Pixelset
