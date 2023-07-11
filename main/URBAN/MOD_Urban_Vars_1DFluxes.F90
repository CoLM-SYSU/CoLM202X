#include <define.h>

#if (defined URBAN_MODEL)
MODULE MOD_Urban_Vars_1DFluxes

! -------------------------------
! Created by Hua Yuan, 12/2020
! -------------------------------

  USE MOD_Precision
  IMPLICIT NONE
  SAVE

! -----------------------------------------------------------------
! Fluxes
! -----------------------------------------------------------------
  !REAL(r8), allocatable :: sabroof     (:) !solar absorption of roof [W/m2]
  !REAL(r8), allocatable :: sabwsun     (:) !solar absorption of sunlit wall [W/m2]
  !REAL(r8), allocatable :: sabwsha     (:) !solar absorption of shaded wall [W/m2]
  !REAL(r8), allocatable :: sabgimp     (:) !solar absorption of impervious [W/m2]
  !REAL(r8), allocatable :: sabgper     (:) !solar absorption of pervious [W/m2]

   REAL(r8), allocatable :: fsen_roof   (:) !sensible heat flux from roof [W/m2]
   REAL(r8), allocatable :: fsen_wsun   (:) !sensible heat flux from sunlit wall [W/m2]
   REAL(r8), allocatable :: fsen_wsha   (:) !sensible heat flux from shaded wall [W/m2]
   REAL(r8), allocatable :: fsen_gimp   (:) !sensible heat flux from impervious road [W/m2]
   REAL(r8), allocatable :: fsen_gper   (:) !sensible heat flux from pervious road [W/m2]
   REAL(r8), allocatable :: fsen_urbl   (:) !sensible heat flux from urban vegetation [W/m2]

   REAL(r8), allocatable :: lfevp_roof  (:) !latent heat flux from roof [W/m2]
   REAL(r8), allocatable :: lfevp_gimp  (:) !latent heat flux from impervious road [W/m2]
   REAL(r8), allocatable :: lfevp_gper  (:) !latent heat flux from pervious road [W/m2]
   REAL(r8), allocatable :: lfevp_urbl  (:) !latent heat flux from urban vegetation [W/m2]

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_1D_UrbanFluxes
  PUBLIC :: deallocate_1D_UrbanFluxes
  PUBLIC :: set_1D_UrbanFluxes

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_1D_UrbanFluxes
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM 1d [numurban] variables
  ! --------------------------------------------------------------------

     USE MOD_Precision
     USE MOD_SPMD_Task
     USE MOD_LandUrban
     IMPLICIT NONE

     IF (p_is_worker) THEN
        IF (numurban > 0) THEN
          !allocate (sabroof        (numurban))
          !allocate (sabwsun        (numurban))
          !allocate (sabwsha        (numurban))
          !allocate (sabgimp        (numurban))
          !allocate (sabgper        (numurban))
           allocate (fsen_roof      (numurban))
           allocate (fsen_wsun      (numurban))
           allocate (fsen_wsha      (numurban))
           allocate (fsen_gimp      (numurban))
           allocate (fsen_gper      (numurban))
           allocate (fsen_urbl      (numurban))

           allocate (lfevp_roof     (numurban))
           allocate (lfevp_gimp     (numurban))
           allocate (lfevp_gper     (numurban))
           allocate (lfevp_urbl     (numurban))
         ENDIF
      ENDIF

  END SUBROUTINE allocate_1D_UrbanFluxes

  SUBROUTINE deallocate_1D_UrbanFluxes
  ! --------------------------------------------------------------------
  ! deallocates memory for CLM 1d [numurban] variables
  ! --------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_LandUrban

      IF (p_is_worker) THEN
         IF (numurban > 0) THEN

           !deallocate (sabroof      )
           !deallocate (sabwsun      )
           !deallocate (sabwsha      )
           !deallocate (sabgimp      )
           !deallocate (sabgper      )
            deallocate (fsen_roof    )
            deallocate (fsen_wsun    )
            deallocate (fsen_wsha    )
            deallocate (fsen_gimp    )
            deallocate (fsen_gper    )
            deallocate (fsen_urbl    )

            deallocate (lfevp_roof   )
            deallocate (lfevp_gimp   )
            deallocate (lfevp_gper   )
            deallocate (lfevp_urbl   )

         ENDIF
      ENDIF

  END SUBROUTINE deallocate_1D_UrbanFluxes

  SUBROUTINE set_1D_UrbanFluxes(Values, Nan)
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM 1d [numurban] variables
  ! --------------------------------------------------------------------

     USE MOD_Precision
     USE MOD_SPMD_Task
     USE MOD_LandUrban
     IMPLICIT NONE
     REAL(r8),INTENT(in) :: Values
     REAL(r8),INTENT(in) :: Nan

     IF (p_is_worker) THEN
        IF (numurban > 0) THEN
          !sabroof        (:) = Values
          !sabwsun        (:) = Values
          !sabwsha        (:) = Values
          !sabgimp        (:) = Values
          !sabgper        (:) = Values
           fsen_roof      (:) = Values
           fsen_wsun      (:) = Values
           fsen_wsha      (:) = Values
           fsen_gimp      (:) = Values
           fsen_gper      (:) = Values
           fsen_urbl      (:) = Values

           lfevp_roof     (:) = Values
           lfevp_gimp     (:) = Values
           lfevp_gper     (:) = Values
           lfevp_urbl     (:) = Values
         ENDIF
      ENDIF

  END SUBROUTINE set_1D_UrbanFluxes

END MODULE MOD_Urban_Vars_1DFluxes
#endif
! ---------- EOP ------------
