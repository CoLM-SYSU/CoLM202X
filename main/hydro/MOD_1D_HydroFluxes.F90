#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_1D_HydroFluxes
   
   USE precision
   IMPLICIT NONE

   ! -- fluxes --
   REAL(r8), allocatable :: rsubs_pch (:)  ! mm/s

   REAL(r8), allocatable :: rsurf_hru (:)  ! m/s
   REAL(r8), allocatable :: rsubs_hru (:)  ! m/s

   REAL(r8), allocatable :: riverheight_ta (:) ! time step average of river height [m]
   REAL(r8), allocatable :: rivermomtem_ta (:) ! time step average of river momentum [m^2/s]
   REAL(r8), allocatable :: riverveloct_ta (:) ! time step average of river velocity [m/s]

   REAL(r8), allocatable :: dpond_hru_ta (:) ! time step average of surface water depth [m]
   REAL(r8), allocatable :: momtm_hru_ta (:) ! time step average of surface water momentum [m^2/s]
   REAL(r8), allocatable :: veloc_hru_ta (:) ! time step average of surface water veloctiy [m/s]

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_1D_HydroFluxes
   PUBLIC :: deallocate_1D_HydroFluxes

CONTAINS 

  SUBROUTINE allocate_1D_HydroFluxes

     USE spmd_task
     USE mod_mesh,      only : numelm
     USE mod_landhru,   only : numhru
     USE mod_landpatch, only : numpatch
     IMPLICIT NONE

     INTEGER :: numbasin

     numbasin = numelm

     IF (p_is_worker) THEN
        IF (numpatch > 0) THEN
           allocate (rsubs_pch (numpatch))
        ENDIF
        IF (numbasin > 0) THEN
           allocate (riverheight_ta (numbasin))
           allocate (rivermomtem_ta (numbasin))
           allocate (riverveloct_ta (numbasin))
        ENDIF
        IF (numhru > 0) THEN
           allocate (rsurf_hru (numhru))
           allocate (rsubs_hru (numhru))
           allocate (dpond_hru_ta (numhru))
           allocate (momtm_hru_ta (numhru))
           allocate (veloc_hru_ta (numhru))
        ENDIF
     ENDIF

  END SUBROUTINE allocate_1D_HydroFluxes

  SUBROUTINE deallocate_1D_HydroFluxes

     IMPLICIT NONE

     IF (allocated(rsubs_pch)) deallocate(rsubs_pch)

     IF (allocated(rsurf_hru)) deallocate(rsurf_hru)
     IF (allocated(rsubs_hru)) deallocate(rsubs_hru)
     
     IF (allocated(riverheight_ta)) deallocate(riverheight_ta)
     IF (allocated(rivermomtem_ta)) deallocate(rivermomtem_ta)
     IF (allocated(riverveloct_ta)) deallocate(riverveloct_ta)
     
     IF (allocated(dpond_hru_ta)) deallocate(dpond_hru_ta)
     IF (allocated(momtm_hru_ta)) deallocate(momtm_hru_ta)
     IF (allocated(veloc_hru_ta)) deallocate(veloc_hru_ta)

  END SUBROUTINE deallocate_1D_HydroFluxes


END MODULE MOD_1D_HydroFluxes
#endif
