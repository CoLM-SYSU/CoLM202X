#include <define.h>

#if (defined UNSTRUCTURED || defined CATCHMENT) 
MODULE mod_unstructured_mesh

   USE precision
   IMPLICIT NONE
   
   INTEGER, allocatable :: hru_patch_s (:)  ! start patch index of a hydrounit 
   INTEGER, allocatable :: hru_patch_e (:)  ! end   patch index of a hydrounit 
   INTEGER, allocatable :: patch2hru (:)    ! hydrounit index of a patch
   
#ifdef CATCHMENT
   INTEGER, allocatable :: cat_hru_s (:)  !start hydrounit index of a catchment
   INTEGER, allocatable :: cat_hru_e (:)  !end   hydrounit index of a catchment
#endif

   REAL(r8), allocatable :: wtpatch_hru (:)

CONTAINS
   
   ! ----------
   SUBROUTINE unstructured_mesh_init ()

      IMPLICIT NONE

#ifdef CATCHMENT
      CALL map_catchment_to_hydrounit ()
#endif
      CALL map_hydrounit_to_patch ()

   END SUBROUTINE unstructured_mesh_init 

   ! ----------
#ifdef CATCHMENT
   SUBROUTINE map_catchment_to_hydrounit

      USE spmd_task
      USE mod_landbasin
      USE mod_hydrounit
      IMPLICIT NONE

      INTEGER :: icat, ihru
      
      IF (p_is_worker) THEN

         IF (numhru <= 0) return

         allocate (cat_hru_s (numbasin))
         allocate (cat_hru_e (numbasin))

         icat = hydrounit%ibasin(1)
         cat_hru_s(icat) = 1
         DO ihru = 1, numhru
            IF (hydrounit%ibasin(ihru) /= icat) THEN
               cat_hru_e(icat) = ihru - 1
               icat = hydrounit%ibasin(ihru)
               cat_hru_s(icat) = ihru
            ENDIF
         ENDDO
         cat_hru_e(icat) = numhru

      ENDIF

   END SUBROUTINE map_catchment_to_hydrounit
#endif

   ! ----------
   SUBROUTINE map_hydrounit_to_patch

      USE spmd_task
      USE mod_utils
      USE mod_pixel
      USE mod_landbasin
      USE mod_hydrounit
      USE mod_landpatch
      IMPLICIT NONE

      INTEGER :: ihru, ipatch, ibasin, ipxl
      
      IF (p_is_worker) THEN

         allocate (hru_patch_s (numhru))
         allocate (hru_patch_e (numhru))
         IF (numpatch > 0) THEN
            allocate (patch2hru (numpatch))
         ENDIF

         ipatch = 1
         DO ihru = 1, numhru
            IF (ipatch > numpatch) THEN
               hru_patch_s(ihru) = -1
               hru_patch_e(ihru) = -1
            ELSE
               IF ((landpatch%bindex(ipatch) /= hydrounit%bindex(ihru)) &
                  .or. (landpatch%ipxstt(ipatch) < hydrounit%ipxstt(ihru)) &
                  .or. (landpatch%ipxend(ipatch) > hydrounit%ipxend(ihru))) THEN
                  hru_patch_s(ihru) = -1
                  hru_patch_e(ihru) = -1
               ELSE
                  hru_patch_s(ihru) = ipatch
                  DO WHILE (  (landpatch%bindex(ipatch) == hydrounit%bindex(ihru)) &
                        .and. (landpatch%ipxstt(ipatch) >= hydrounit%ipxstt(ihru)) &
                        .and. (landpatch%ipxend(ipatch) <= hydrounit%ipxend(ihru))) 

                     patch2hru(ipatch) = ihru

                     ipatch = ipatch + 1
                     IF (ipatch > numpatch) exit
                  ENDDO
                  hru_patch_e(ihru) = ipatch - 1
               ENDIF
            ENDIF
         ENDDO

         IF (numpatch > 0) THEN
            allocate (wtpatch_hru (numpatch))
         ENDIF
         
         DO ipatch = 1, numpatch
            
            wtpatch_hru(ipatch) = 0

            ibasin = landpatch%ibasin(ipatch)
            DO ipxl = landpatch%ipxstt(ipatch), landpatch%ipxend(ipatch) 
               wtpatch_hru(ipatch) = wtpatch_hru(ipatch) & 
                  + areaquad (&
                  pixel%lat_s(landbasin(ibasin)%ilat(ipxl)), &
                  pixel%lat_n(landbasin(ibasin)%ilat(ipxl)), &
                  pixel%lon_w(landbasin(ibasin)%ilon(ipxl)), &
                  pixel%lon_e(landbasin(ibasin)%ilon(ipxl)) )
            ENDDO
         ENDDO

         DO ihru = 1, numhru
            IF ((hru_patch_s(ihru) /= -1) .and. (hru_patch_e(ihru) /= -1)) THEN
               wtpatch_hru(hru_patch_s(ihru):hru_patch_e(ihru)) =  &
                  wtpatch_hru(hru_patch_s(ihru):hru_patch_e(ihru)) &
                  / sum(wtpatch_hru(hru_patch_s(ihru):hru_patch_e(ihru)))
            ENDIF
         ENDDO

      ENDIF

   END SUBROUTINE map_hydrounit_to_patch

   ! ----------
   SUBROUTINE unstructured_mesh_final ()

      IMPLICIT NONE

      ! Local Variables
      INTEGER :: icat

      IF (allocated(hru_patch_s)) deallocate(hru_patch_s)
      IF (allocated(hru_patch_e)) deallocate(hru_patch_e)
      IF (allocated(patch2hru  )) deallocate(patch2hru  )
      
#ifdef CATCHMENT
      IF (allocated(cat_hru_s)) deallocate(cat_hru_s)
      IF (allocated(cat_hru_e)) deallocate(cat_hru_e)
#endif

   END SUBROUTINE unstructured_mesh_final

END MODULE mod_unstructured_mesh
#endif
