#include <define.h>

MODULE mod_pixelset

   USE precision
   USE mod_data_type
   IMPLICIT NONE

   ! ---- data types ----
   TYPE :: vec_gather_scatter_type
      
      ! for worker and io
      INTEGER, allocatable :: vlen(:,:)

      ! for worker
      INTEGER, allocatable :: vstt(:,:)
      INTEGER, allocatable :: vend(:,:)

      ! for io
      INTEGER, allocatable :: vcnt(:,:,:)
      INTEGER, allocatable :: vdsp(:,:,:)

   CONTAINS
      final  :: vec_gather_scatter_free_mem

   END TYPE vec_gather_scatter_type

   ! ---- data types ----
   TYPE :: pixelset_type
     
      INTEGER :: nset

      INTEGER, allocatable :: unum(:)
      INTEGER, allocatable :: iunt(:)
      INTEGER, allocatable :: istt(:)
      INTEGER, allocatable :: iend(:)
      
      INTEGER, allocatable :: ltyp(:)

      LOGICAL, allocatable  :: nonzero(:,:)

      TYPE(vec_gather_scatter_type) :: vecgs

   CONTAINS 
      procedure, PUBLIC :: set_vecgs         => vec_gather_scatter_set
      procedure, PUBLIC :: get_lonlat_radian => pixelset_get_lonlat_radian
      procedure, PUBLIC :: pset_pack         => pixelset_pack
      final :: pixelset_free_mem

   END TYPE pixelset_type

CONTAINS
   
   ! --------------------------------
   SUBROUTINE pixelset_get_lonlat_radian (this, rlon, rlat)

      USE precision
      USE mod_utils
      USE mod_pixel
      USE mod_landbasin

      IMPLICIT NONE
      CLASS(pixelset_type) :: this

      REAL(r8), intent(inout) :: rlon(:), rlat(:)

      ! Local Variables
      INTEGER :: iset, iu, istt, iend, npxl, ipxl
      REAL(r8), allocatable :: area(:)

      DO iset = 1, this%nset

         iu = this%iunt(iset)

         istt = this%istt (iset)
         iend = this%iend (iset)

         allocate (area (istt:iend))
         DO ipxl = istt, iend 
            area(ipxl) = areaquad (&
               pixel%lat_s(landbasin(iu)%ilat(ipxl)), &
               pixel%lat_n(landbasin(iu)%ilat(ipxl)), &
               pixel%lon_w(landbasin(iu)%ilon(ipxl)), &
               pixel%lon_e(landbasin(iu)%ilon(ipxl)) )
         ENDDO

         npxl = iend - istt + 1
         rlat(iset) = get_pixelset_rlat ( &
            npxl, landbasin(iu)%ilat(istt:iend), area)
         rlon(iset) = get_pixelset_rlon ( &
            npxl, landbasin(iu)%ilon(istt:iend), area)

         deallocate (area)

      ENDDO

   END SUBROUTINE pixelset_get_lonlat_radian 

   ! --------------------------------
   FUNCTION get_pixelset_rlat (npxl, ilat, area) result(rlat)

      USE precision
      USE MathConstants, only : pi
      USE mod_pixel
      IMPLICIT NONE
      
      REAL(r8) :: rlat

      INTEGER,  intent(in) :: npxl
      INTEGER,  intent(in) :: ilat(npxl)
      REAL(r8), intent(in) :: area(npxl)

      ! Local variables
      INTEGER :: ipxl

      rlat = 0.0
      DO ipxl = 1, npxl
         rlat = rlat + (pixel%lat_s(ilat(ipxl)) + pixel%lat_n(ilat(ipxl))) * 0.5 * area(ipxl)
      ENDDO
      rlat = rlat / sum(area) * pi/180.0

   END FUNCTION get_pixelset_rlat

   ! --------------------------------
   FUNCTION get_pixelset_rlon (npxl, ilon, area) result(rlon)

      USE precision
      USE mod_utils
      USE MathConstants, only : pi
      USE mod_pixel
      IMPLICIT NONE

      REAL(r8) :: rlon

      INTEGER,  intent(in) :: npxl
      INTEGER,  intent(in) :: ilon(npxl)
      REAL(r8), intent(in) :: area(npxl)

      ! Local variables
      INTEGER  :: ipxl
      REAL(r8) :: lon, lon0, area_done

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
      TYPE (pixelset_type) :: this

      IF (allocated(this%unum)) deallocate(this%unum)
      IF (allocated(this%iunt)) deallocate(this%iunt)
      IF (allocated(this%istt)) deallocate(this%istt)
      IF (allocated(this%iend)) deallocate(this%iend)
      IF (allocated(this%ltyp)) deallocate(this%ltyp)
      
      IF (allocated(this%nonzero)) deallocate(this%nonzero)

   END SUBROUTINE pixelset_free_mem
   
   ! --------------------------------
   SUBROUTINE vec_gather_scatter_set (this)

      USE mod_block
      USE spmd_task
      USE mod_landbasin
      IMPLICIT NONE

      class(pixelset_type)  :: this

      ! Local variables
      INTEGER :: iproc
      INTEGER :: iset, iu, xblk, yblk, iblk, jblk, scnt
      
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
         
      IF (.not. allocated (this%vecgs%vlen)) THEN
         allocate (this%vecgs%vlen (gblock%nxblk, gblock%nyblk))
      ENDIF

      IF (p_is_worker) THEN

         IF (.not. allocated (this%vecgs%vstt)) THEN
            allocate (this%vecgs%vstt (gblock%nxblk, gblock%nyblk))
            allocate (this%vecgs%vend (gblock%nxblk, gblock%nyblk))
         ENDIF

         this%vecgs%vlen(:,:) = 0
         this%vecgs%vstt(:,:) = 0
         this%vecgs%vend(:,:) = -1

         iu = 1
         xblk = 0
         yblk = 0
         DO iset = 1, size(this%unum)
            DO WHILE (this%unum(iset) /= landbasin(iu)%num)
               iu = iu + 1
            ENDDO 

            IF ((landbasin(iu)%xblk /= xblk) .or. (landbasin(iu)%yblk /= yblk)) THEN
               xblk = landbasin(iu)%xblk
               yblk = landbasin(iu)%yblk
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
                     this%vecgs%vcnt, 1, MPI_INTEGER, p_root, p_comm_group, p_err)

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
         IF (.not. allocated(this%nonzero)) THEN
            allocate (this%nonzero (gblock%nxblk,gblock%nyblk))
         ENDIF

         this%nonzero = this%vecgs%vlen > 0
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, this%nonzero, gblock%nxblk * gblock%nyblk, &
            MPI_LOGICAL, MPI_LOR, p_comm_group, p_err)
#endif
      ENDIF

   END SUBROUTINE vec_gather_scatter_set 
   
   ! --------------------------------
   SUBROUTINE pixelset_pack (this, mask, nset_packed)

      USE spmd_task
      IMPLICIT NONE
      class(pixelset_type) :: this
      LOGICAL, intent(in)  :: mask(:)
      INTEGER, intent(out) :: nset_packed

      INTEGER, allocatable :: unum1(:)
      INTEGER, allocatable :: iunt1(:)
      INTEGER, allocatable :: istt1(:)
      INTEGER, allocatable :: iend1(:)
      INTEGER, allocatable :: ltyp1(:)

      IF (p_is_worker) THEN 
      
         IF (count(mask) < this%nset) THEN

            allocate (unum1(this%nset))
            allocate (iunt1(this%nset))
            allocate (istt1(this%nset))
            allocate (iend1(this%nset))
            allocate (ltyp1(this%nset))

            unum1 = this%unum
            iunt1 = this%iunt
            istt1 = this%istt
            iend1 = this%iend
            ltyp1 = this%ltyp

            deallocate (this%unum)
            deallocate (this%iunt)
            deallocate (this%istt)
            deallocate (this%iend)
            deallocate (this%ltyp)

            this%nset = count(mask)

            IF (this%nset > 0) THEN

               allocate (this%unum(this%nset))
               allocate (this%iunt(this%nset))
               allocate (this%istt(this%nset))
               allocate (this%iend(this%nset))
               allocate (this%ltyp(this%nset))

               this%unum = pack(unum1, mask)
               this%iunt = pack(iunt1, mask)
               this%istt = pack(istt1, mask)
               this%iend = pack(iend1, mask)
               this%ltyp = pack(ltyp1, mask)

            ENDIF

            deallocate (unum1)
            deallocate (iunt1)
            deallocate (istt1)
            deallocate (iend1)
            deallocate (ltyp1)

         ENDIF
      
      ENDIF
         
      CALL this%set_vecgs

      nset_packed = this%nset

   END SUBROUTINE pixelset_pack

   ! --------------------------------
   SUBROUTINE vec_gather_scatter_free_mem (this)
      
      USE mod_block
      USE mod_data_type
      USE spmd_task
      IMPLICIT NONE

      TYPE (vec_gather_scatter_type) :: this
      ! Local variables
      INTEGER :: i, j

      IF (allocated(this%vlen))  deallocate (this%vlen)
      IF (allocated(this%vstt))  deallocate (this%vstt)
      IF (allocated(this%vend))  deallocate (this%vend)
      IF (allocated(this%vcnt))  deallocate (this%vcnt)
      IF (allocated(this%vdsp))  deallocate (this%vdsp)
   
   END SUBROUTINE vec_gather_scatter_free_mem

END MODULE mod_pixelset
