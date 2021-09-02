#include <define.h>

PROGRAM CLMINI
! ======================================================================
! Initialization of Land Characteristic Parameters and Initial State Variables
!
! Reference:
!     [1] Dai et al., 2003: The Common Land Model (CoLM).
!         Bull. of Amer. Meter. Soc., 84: 1013-1023
!     [2] Dai et al., 2004: A two-big-leaf model for canopy temperature,
!         photosynthesis and stomatal conductance. Journal of Climate
!     [3] Dai et al., 2014: The Terrestrial Modeling System (TMS).
!
!     Created by Yongjiu Dai Februay 2004
!     Revised by Yongjiu Dai Februay 2014
! ======================================================================

      USE precision
      USE timemanager
      USE GlobalVars
      USE LC_Const
      USE PFT_Const
      
      IMPLICIT NONE

! ----------------local variables ---------------------------------
      CHARACTER(LEN=256) :: casename ! case name
      CHARACTER(LEN=256) :: dir_model_landdata
      CHARACTER(LEN=256) :: dir_restart_hist
      CHARACTER(LEN=256) :: dir_infolist
      INTEGER :: s_year      ! starting date for run in year
      INTEGER :: s_julian    ! starting date for run in julian day
      INTEGER :: s_month     ! starting month for run 
      INTEGER :: s_day       ! starting day for run 
      INTEGER :: s_seconds   ! starting time of day for run in seconds
      INTEGER :: idate(3)    ! starting date
      LOGICAL :: greenwich   ! true: greenwich time, false: local time

      INTEGER :: lon_points  ! number of longitude points on model grid
      INTEGER :: lat_points  ! number of latitude points on model grid

      CHARACTER(LEN=256) :: finfolist      ! file name of run information

!  Required by atmospheric models's initialization (such as GRAPES/WRF/RSM/EMSs) 
      REAL(r8), allocatable :: tg_xy   (:,:)
      REAL(r8), allocatable :: albvb_xy(:,:)
      REAL(r8), allocatable :: albvd_xy(:,:)
      REAL(r8), allocatable :: albnb_xy(:,:)
      REAL(r8), allocatable :: albnd_xy(:,:)
      REAL(r8), allocatable :: trad_xy (:,:)
      REAL(r8), allocatable :: rib_xy  (:,:)
      REAL(r8), allocatable :: fm_xy   (:,:)
      REAL(r8), allocatable :: fh_xy   (:,:)
      REAL(r8), allocatable :: fq_xy   (:,:)

      namelist /clminiexp/ casename,dir_model_landdata,&
                           dir_restart_hist,dir_infolist,&
                           lon_points,lat_points,greenwich,&
                           s_year,s_month,s_day,s_seconds
! ----------------------------------------------------------------------
      read (5,clminiexp)

      CALL Init_GlovalVars
      CAll Init_LC_Const
      CAll Init_PFT_Const

      numpatch = 0
      numpft   = 0
      numpc    = 0
      
      CALL monthday2julian(s_year,s_month,s_day,s_julian)
      idate(1) = s_year; idate(2) = s_julian; idate(3) = s_seconds

      allocate ( tg_xy   (lon_points,lat_points) )
      allocate ( albvb_xy(lon_points,lat_points) )
      allocate ( albvd_xy(lon_points,lat_points) )
      allocate ( albnb_xy(lon_points,lat_points) )
      allocate ( albnd_xy(lon_points,lat_points) )
      allocate ( trad_xy (lon_points,lat_points) )
      allocate ( rib_xy  (lon_points,lat_points) )
      allocate ( fm_xy   (lon_points,lat_points) )
      allocate ( fh_xy   (lon_points,lat_points) )
      allocate ( fq_xy   (lon_points,lat_points) )

      CALL initialize (casename,dir_model_landdata,dir_restart_hist,&
                       idate,greenwich,lon_points,lat_points,&
                       tg_xy,albvb_xy,albvd_xy,albnb_xy,albnd_xy,&
                       trad_xy,rib_xy,fm_xy,fh_xy,fq_xy)

      finfolist = trim(dir_infolist)//'clmini.infolist'
      open(100,file=trim(finfolist),form='formatted')
      write(100,*) 'numpatch  = ', numpatch   !1.1
      write(100,*) 'numpft    = ', numpft     !1.2
      write(100,*) 'numpc     = ', numpc      !1.3
      IF ( greenwich ) THEN
      write(100,*) 'greenwich =       .true.' !2
      ELSE
      write(100,*) 'greenwich =      .false.' !2
      ENDIF
      write(100,*) 's_year    = ', s_year     !3
      write(100,*) 's_month   = ', s_month    !4
      write(100,*) 's_day     = ', s_day      !5
      write(100,*) 's_seconds = ', s_seconds  !6
      write(100,*) '/'
      CLOSE(100)

      deallocate ( tg_xy    )
      deallocate ( albvb_xy )
      deallocate ( albvd_xy )
      deallocate ( albnb_xy )
      deallocate ( albnd_xy )
      deallocate ( trad_xy  )
      deallocate ( rib_xy   )
      deallocate ( fm_xy    )
      deallocate ( fh_xy    )
      deallocate ( fq_xy    )

      write(6,*) 'CLM Initialization Execution Completed'

END PROGRAM CLMINI
! ----------------------------------------------------------------------
! EOP
