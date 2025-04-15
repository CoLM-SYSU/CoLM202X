#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_Vars_1DFluxes
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
   real(r8), allocatable :: xsubs_elm (:)  ! subsurface lateral flow between basins                      [m/s]
   real(r8), allocatable :: xsubs_hru (:)  ! subsurface lateral flow between hydrological response units [m/s]
   real(r8), allocatable :: xsubs_pch (:)  ! subsurface lateral flow between patches inside one HRU      [m/s]

   real(r8), allocatable :: wdsrf_bsn_ta (:) ! time step average of river height   [m]
   real(r8), allocatable :: momen_riv_ta (:) ! time step average of river momentum [m^2/s]
   real(r8), allocatable :: veloc_riv_ta (:) ! time step average of river velocity [m/s]
   real(r8), allocatable :: discharge_ta (:) ! river discharge [m^3/s]

   real(r8), allocatable :: wdsrf_bsnhru_ta (:) ! time step average of surface water depth    [m]
   real(r8), allocatable :: momen_bsnhru_ta (:) ! time step average of surface water momentum [m^2/s]
   real(r8), allocatable :: veloc_bsnhru_ta (:) ! time step average of surface water veloctiy [m/s]

   real(r8), allocatable :: xwsur   (:) ! surface water exchange [mm h2o/s]
   real(r8), allocatable :: xwsub   (:) ! subsurface water exchange [mm h2o/s]
   real(r8), allocatable :: fldarea (:) ! fraction of flooded area [-]

   real(r8), allocatable :: ntacc_bsn (:)

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_1D_CatchFluxes
   PUBLIC :: deallocate_1D_CatchFluxes

CONTAINS

   SUBROUTINE allocate_1D_CatchFluxes

   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   USE MOD_Mesh,        only: numelm
   USE MOD_LandHRU,     only: numhru
   USE MOD_LandPatch,   only: numpatch
   USE MOD_Catch_BasinNetwork, only: numbasin, numbsnhru
   IMPLICIT NONE

      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN
            allocate (xsubs_pch (numpatch)) ; xsubs_pch (:) = spval
            allocate (xwsur     (numpatch)) ; xwsur     (:) = spval
            allocate (xwsub     (numpatch)) ; xwsub     (:) = spval
            allocate (fldarea   (numpatch)) ; fldarea   (:) = spval
         ENDIF

         IF (numelm > 0) THEN
            allocate (xsubs_elm (numelm)) ; xsubs_elm(:) = spval
         ENDIF

         IF (numbasin > 0) THEN
            allocate (wdsrf_bsn_ta (numbasin)) ; wdsrf_bsn_ta (:) = spval
            allocate (momen_riv_ta (numbasin)) ; momen_riv_ta (:) = spval
            allocate (veloc_riv_ta (numbasin)) ; veloc_riv_ta (:) = spval
            allocate (discharge_ta (numbasin)) ; discharge_ta (:) = spval
         ENDIF

         IF (numhru > 0) THEN
            allocate (xsubs_hru (numhru)); xsubs_hru(:) = spval
         ENDIF

         IF (numbsnhru > 0) THEN
            allocate (wdsrf_bsnhru_ta (numbsnhru)) ; wdsrf_bsnhru_ta (:) = spval
            allocate (momen_bsnhru_ta (numbsnhru)) ; momen_bsnhru_ta (:) = spval
            allocate (veloc_bsnhru_ta (numbsnhru)) ; veloc_bsnhru_ta (:) = spval
         ENDIF

         IF (numbasin > 0) allocate (ntacc_bsn (numbasin))
         IF (numbasin > 0) ntacc_bsn(:) = 0.

      ENDIF

   END SUBROUTINE allocate_1D_CatchFluxes

   SUBROUTINE deallocate_1D_CatchFluxes

   IMPLICIT NONE

      IF (allocated(xsubs_elm)) deallocate(xsubs_elm)
      IF (allocated(xsubs_hru)) deallocate(xsubs_hru)
      IF (allocated(xsubs_pch)) deallocate(xsubs_pch)

      IF (allocated(wdsrf_bsn_ta)) deallocate(wdsrf_bsn_ta)
      IF (allocated(momen_riv_ta)) deallocate(momen_riv_ta)
      IF (allocated(veloc_riv_ta)) deallocate(veloc_riv_ta)
      IF (allocated(discharge_ta)) deallocate(discharge_ta)

      IF (allocated(wdsrf_bsnhru_ta)) deallocate(wdsrf_bsnhru_ta)
      IF (allocated(momen_bsnhru_ta)) deallocate(momen_bsnhru_ta)
      IF (allocated(veloc_bsnhru_ta)) deallocate(veloc_bsnhru_ta)

      IF (allocated(xwsur  )) deallocate(xwsur  )
      IF (allocated(xwsub  )) deallocate(xwsub  )
      IF (allocated(fldarea)) deallocate(fldarea)

      IF (allocated(ntacc_bsn)) deallocate(ntacc_bsn)

   END SUBROUTINE deallocate_1D_CatchFluxes

END MODULE MOD_Catch_Vars_1DFluxes
#endif
