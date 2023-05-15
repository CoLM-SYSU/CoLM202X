#include <define.h>

MODULE MOD_1D_UrbanFluxes

! -------------------------------
! Created by Hua Yuan, 12/2020
! -------------------------------

  USE precision
  IMPLICIT NONE
  SAVE

! -----------------------------------------------------------------
! Fluxes
! -----------------------------------------------------------------
  !REAL(r8), allocatable :: sabroof   (:) !solar absorption of roof [W/m2]
  !REAL(r8), allocatable :: sabwsun   (:) !solar absorption of sunlit wall [W/m2]
  !REAL(r8), allocatable :: sabwsha   (:) !solar absorption of shaded wall [W/m2]
  !REAL(r8), allocatable :: sabgimp   (:) !solar absorption of impervious [W/m2]
  !REAL(r8), allocatable :: sabgper   (:) !solar absorption of pervious [W/m2]

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_1D_UrbanFluxes
  PUBLIC :: deallocate_1D_UrbanFluxes

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_1D_UrbanFluxes
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM 1d [numurban] variables
  ! --------------------------------------------------------------------

     USE precision
     USE GlobalVars
     IMPLICIT NONE

     !allocate (sabroof (numurban))
     !allocate (sabwsun (numurban))
     !allocate (sabwsha (numurban))
     !allocate (sabgimp (numurban))
     !allocate (sabgper (numurban))

  END SUBROUTINE allocate_1D_UrbanFluxes

  SUBROUTINE deallocate_1D_UrbanFluxes
  ! --------------------------------------------------------------------
  ! deallocates memory for CLM 1d [numurban] variables
  ! --------------------------------------------------------------------

     !deallocate (sabroof )
     !deallocate (sabwsun )
     !deallocate (sabwsha )
     !deallocate (sabgimp )
     !deallocate (sabgper )

  END SUBROUTINE deallocate_1D_UrbanFluxes

END MODULE MOD_1D_UrbanFluxes
! ---------- EOP ------------
