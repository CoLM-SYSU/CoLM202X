#include <define.h>

MODULE GlobalVars
!-------------------------------------------------------------------------------
!
! !DESCRIPTION:
! Define some global variables
!
! REVISIONS:
! Hua Yuan, 08/2019: initial version partly adapted from CoLM2014
! TODO ...
!
! !USES:
   USE precision
   IMPLICIT NONE
   SAVE

#ifdef USGS_CLASSIFICATION
   ! GLCC USGS number of land cover category
   INTEGER, parameter :: N_land_classification = 24
   ! GLCC USGS land cover named index (could be added IF needed)
   INTEGER, parameter :: URBAN    = 1
   INTEGER, parameter :: WATERBODY= 16
#else
   ! MODIS IGBP number of land cover category
   INTEGER, parameter :: N_land_classification = 17
   ! MODIS IGBP land cover named index (could be added IF needed)
   INTEGER, parameter :: WETLAND  = 11
   INTEGER, parameter :: CROPLAND = 12
   INTEGER, parameter :: URBAN    = 13
   INTEGER, parameter :: GLACIERS = 15
   INTEGER, parameter :: WATERBODY= 17
#endif

   ! number of plant functional types
#ifndef CROP
   INTEGER, parameter :: N_PFT    = 16
   INTEGER, parameter :: N_CFT    = 0
#else
   INTEGER, parameter :: N_PFT    = 15
   INTEGER, parameter :: N_CFT    = 64
#endif

#ifdef USE_LCZ
   INTEGER, parameter :: N_URB    = 10
#else
   INTEGER, parameter :: N_URB    = 3
#endif

   ! vertical layer number
   INTEGER, parameter :: maxsnl   = -5
   INTEGER, parameter :: nl_soil  = 10
   INTEGER, parameter :: nl_soil_full  = 15

   INTEGER, parameter :: nl_lake  = 10
   INTEGER, parameter :: nl_roof  = 10
   INTEGER, parameter :: nl_wall  = 10
   INTEGER, parameter :: nvegwcs  = 4  ! number of vegetation water potential nodes

   ! bgc variables
   integer, parameter :: ndecomp_pools        = 7
   integer, parameter :: ndecomp_transitions  = 10
   integer, parameter :: npcropmin            = 17
   real(r8),parameter :: zmin_bedrock         = 0.4
   integer, parameter :: nbedrock             = 10
   integer, parameter :: ndecomp_pools_vr     = ndecomp_pools * nl_soil

   ! crop index
   integer, parameter :: noveg                = 0
   integer, parameter :: nbrdlf_evr_shrub     = 9
   integer, parameter :: nbrdlf_dcd_brl_shrub = 11
   integer, parameter :: nc3crop              = 15
   integer, parameter :: nc3irrig             = 16
   integer, parameter :: ntmp_corn            = 17 ! temperate_corn
   integer, parameter :: nirrig_tmp_corn      = 18 ! irrigated temperate corn
   integer, parameter :: nswheat              = 19 ! spring wheat
   integer, parameter :: nirrig_swheat        = 20 ! irrigated spring wheat
   integer, parameter :: nwwheat              = 21 ! winter wheat
   integer, parameter :: nirrig_wwheat        = 22 ! irrigated winter wheat
   integer, parameter :: ntmp_soybean         = 23 ! temperate soybean
   integer, parameter :: nirrig_tmp_soybean   = 24 ! irrigated temperate soybean
   integer, parameter :: ncotton              = 41 ! cotton
   integer, parameter :: nirrig_cotton        = 42 ! irrigated cotton
   integer, parameter :: nrice                = 61 ! rice
   integer, parameter :: nirrig_rice          = 62 ! irrigated rice
   integer, parameter :: nsugarcane           = 67 ! sugarcane
   integer, parameter :: nirrig_sugarcane     = 68 ! irrigated sugarcane
   integer, parameter :: nmiscanthus          = 71 ! miscanthus
   integer, parameter :: nirrig_miscanthus    = 72 ! irrigated miscanthus
   integer, parameter :: nswitchgrass         = 73 ! switchgrass
   integer, parameter :: nirrig_switchgrass   = 74 ! irrigated switchgrass
   integer, parameter :: ntrp_corn            = 75 ! tropical corn
   integer, parameter :: nirrig_trp_corn      = 76 ! irrigated tropical corn
   integer, parameter :: ntrp_soybean         = 77 ! tropical soybean
   integer, parameter :: nirrig_trp_soybean   = 78 ! irrigated tropical soybean

   !TODO: need moved when coupling urban model
   !INTEGER, parameter :: numurban = 1  !total number of Urban patches of grids

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

!     ndecomp_pools_vr = ndecomp_pools * nl_soil

   END SUBROUTINE Init_GlovalVars

END MODULE GlobalVars
! ---------- EOP ------------
