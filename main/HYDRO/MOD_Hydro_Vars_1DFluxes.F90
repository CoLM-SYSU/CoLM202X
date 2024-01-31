#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Hydro_Vars_1DFluxes
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!   
!   1D fluxes in lateral hydrological processes.
!
! Created by Shupeng Zhang, May 2023
!-------------------------------------------------------------------------------------
   
   USE MOD_Precision
   IMPLICIT NONE

   ! -- fluxes --
   real(r8), allocatable :: xsubs_bsn (:)  ! subsurface lateral flow between basins                      [m/s]
   real(r8), allocatable :: xsubs_hru (:)  ! subsurface lateral flow between hydrological response units [m/s]
   real(r8), allocatable :: xsubs_pch (:)  ! subsurface lateral flow between patches inside one HRU      [m/s]

   real(r8), allocatable :: wdsrf_bsn_ta (:) ! time step average of river height   [m]
   real(r8), allocatable :: momen_riv_ta (:) ! time step average of river momentum [m^2/s]
   real(r8), allocatable :: veloc_riv_ta (:) ! time step average of river velocity [m/s]

   real(r8), allocatable :: wdsrf_hru_ta (:) ! time step average of surface water depth    [m]
   real(r8), allocatable :: momen_hru_ta (:) ! time step average of surface water momentum [m^2/s]
   real(r8), allocatable :: veloc_hru_ta (:) ! time step average of surface water veloctiy [m/s]
  
   real(r8), allocatable :: xwsur (:) ! surface water exchange [mm h2o/s]
   real(r8), allocatable :: xwsub (:) ! subsurface water exchange [mm h2o/s]
   
   real(r8), allocatable :: discharge (:) ! river discharge [m^3/s]

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_1D_HydroFluxes
   PUBLIC :: deallocate_1D_HydroFluxes

CONTAINS 

   SUBROUTINE allocate_1D_HydroFluxes

   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only : spval
   USE MOD_Mesh,      only : numelm
   USE MOD_LandHRU,   only : numhru
   USE MOD_LandPatch, only : numpatch
   IMPLICIT NONE

   integer :: numbasin

      numbasin = numelm

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (xsubs_pch (numpatch)) ; xsubs_pch (:) = spval
            allocate (xwsur     (numpatch)) ; xwsur     (:) = spval
            allocate (xwsub     (numpatch)) ; xwsub     (:) = spval
         ENDIF
         IF (numbasin > 0) THEN
            allocate (xsubs_bsn    (numbasin)) ; xsubs_bsn    (:) = spval
            allocate (wdsrf_bsn_ta (numbasin)) ; wdsrf_bsn_ta (:) = spval
            allocate (momen_riv_ta (numbasin)) ; momen_riv_ta (:) = spval
            allocate (veloc_riv_ta (numbasin)) ; veloc_riv_ta (:) = spval
            allocate (discharge    (numbasin)) ; discharge    (:) = spval
         ENDIF
         IF (numhru > 0) THEN
            allocate (xsubs_hru    (numhru)) ; xsubs_hru    (:) = spval
            allocate (wdsrf_hru_ta (numhru)) ; wdsrf_hru_ta (:) = spval
            allocate (momen_hru_ta (numhru)) ; momen_hru_ta (:) = spval
            allocate (veloc_hru_ta (numhru)) ; veloc_hru_ta (:) = spval
         ENDIF
      ENDIF

   END SUBROUTINE allocate_1D_HydroFluxes

   SUBROUTINE deallocate_1D_HydroFluxes

   IMPLICIT NONE

      IF (allocated(xsubs_pch)) deallocate(xsubs_pch)
      IF (allocated(xsubs_hru)) deallocate(xsubs_hru)
      IF (allocated(xsubs_bsn)) deallocate(xsubs_bsn)
      
      IF (allocated(wdsrf_bsn_ta)) deallocate(wdsrf_bsn_ta)
      IF (allocated(momen_riv_ta)) deallocate(momen_riv_ta)
      IF (allocated(veloc_riv_ta)) deallocate(veloc_riv_ta)
      
      IF (allocated(wdsrf_hru_ta)) deallocate(wdsrf_hru_ta)
      IF (allocated(momen_hru_ta)) deallocate(momen_hru_ta)
      IF (allocated(veloc_hru_ta)) deallocate(veloc_hru_ta)
      
      IF (allocated(xwsur)) deallocate(xwsur)
      IF (allocated(xwsub)) deallocate(xwsub)

      IF (allocated(discharge)) deallocate(discharge)

   END SUBROUTINE deallocate_1D_HydroFluxes

END MODULE MOD_Hydro_Vars_1DFluxes
#endif
