#include <define.h>

MODULE GlobalVars
!-------------------------------------------------------------------------------
! PURPOSE:
!       define some global variables 
!-------------------------------------------------------------------------------
   USE precision
   IMPLICIT NONE
   SAVE
   
#ifdef USGS_CLASSIFICATION
   ! GLCC USGS number of land cover category
   INTEGER, parameter :: N_land_classification = 24 
#else 
   ! MODIS IGBP number of land cover category
   INTEGER, parameter :: N_land_classification = 17 
#endif
   
   ! number of plant functional types
   INTEGER, parameter :: N_PFT    = 16 
   
   ! vertical layer number
   INTEGER, parameter :: maxsnl   = -5
   INTEGER, parameter :: nl_soil  = 10
   INTEGER, parameter :: nl_lake  = 10
   INTEGER, parameter :: nl_roof  = 5
   INTEGER, parameter :: nl_wall  = 5
   
   INTEGER, parameter :: numurban = 1  !total number of Urban patches of grids
   INTEGER :: numpatch                 !total number of patches of grids
   INTEGER :: numpft                   !total number of PFT patches of grids
   INTEGER :: numpc                    !total number of PC patches of grids

   REAL(r8) :: z_soi (1:nl_soil)       !node depth [m]
   REAL(r8) :: z_soih(1:nl_soil)       !interface level below a zsoi level [m]
   REAL(r8) :: zi_soi(1:nl_soil)       !interface level below a zsoi level [m]
   REAL(r8) :: dz_soi(1:nl_soil)       !soil node thickness [m]
  
   REAL(r8), parameter :: spval = -1.e36_r8  !missing value
   REAL(r8), parameter :: PI    = 4*atan(1.) !pi value

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Init_GlovalVars

CONTAINS 

   SUBROUTINE Init_GlovalVars

      IMPLICIT NONE
      
      INTEGER :: nsl

      DO nsl = 1, nl_soil
         z_soi(nsl) = 0.025*(exp(0.5*(nsl-0.5))-1.)  !node depths
      ENDDO

      dz_soi(1) = 0.5*(z_soi(1)+z_soi(2))            !=zsoih(1)
      dz_soi(nl_soil) = z_soi(nl_soil)-z_soi(nl_soil-1)
      DO nsl = 2, nl_soil-1
         ! thickness between two interfaces
         dz_soi(nsl) = 0.5*(z_soi(nsl+1)-z_soi(nsl-1)) 
      ENDDO

      z_soih(nl_soil) = z_soi(nl_soil) + 0.5*dz_soi(nl_soil)
      DO nsl = 1, nl_soil-1
         z_soih(nsl) = 0.5*(z_soi(nsl)+z_soi(nsl+1)) !interface depths
      ENDDO

      zi_soi(1) = dz_soi(1)
      DO nsl = 2, nl_soil
         zi_soi(nsl) = zi_soi(nsl-1) + dz_soi(nsl)
      ENDDO

   END SUBROUTINE Init_GlovalVars

END MODULE GlobalVars
! ---------- EOP ------------
