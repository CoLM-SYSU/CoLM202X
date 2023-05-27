MODULE MathConstants

!=======================================================================
! mathematical constants
!=======================================================================

  USE MOD_Precision
  IMPLICIT NONE

  !NOTE, yuan: can be moved to GlobalVars.F90
  PUBLIC
  REAL(r8), parameter :: pi = 3.14159265358979323_r8
  REAL(r8), parameter :: deg2rad = 1.745329251994330e-2_r8

END MODULE MathConstants
