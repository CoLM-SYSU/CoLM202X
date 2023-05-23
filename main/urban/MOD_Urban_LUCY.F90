#include <define.h>

MODULE MOD_Urban_LUCY
  ! -----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Anthropogenic model to calculate anthropogenic heat flux for the rest
  !
  ! ORIGINAL:
  ! Wenzong Dong, May, 2022
  !
  ! -----------------------------------------------------------------------
  ! !USE
  USE precision
  USE MOD_Vars_Global
  USE MOD_Vars_PhysicalConst
  USE timemanager, only: julian2monthday, isleapyear
  IMPLICIT NONE
  SAVE
  PRIVATE :: timeweek, gmt2local
  PUBLIC  :: LUCY

CONTAINS

  ! -----------------------------------------------------------------------
  SUBROUTINE LUCY( idate       , deltim  , patchlonr, fix_holiday, &
                   week_holiday, hum_prof, wdh_prof , weh_prof   , popcell, &
                   vehicle     , Fahe    , vehc     , meta        )

  ! !DESCRIPTION:
  ! Anthropogenic heat fluxes other than building heat were calculated
  !
  ! The calling sequence is:
  ! -> gmt2local: convert GMT time to local time
  ! -> timeweek : calculate the day of the week
  !
  ! REFERENCES:
  ! 1) Grimmond, C. S. B. (1992). The suburban energy balance: Methodological considerations and results
  ! for a mid-latitude west coast city under winter and spring conditions. International Journal of Climatology,
  ! 12(5), 481–497. https://doi.org/10.1002/joc.3370120506
  ! 2) Allen, L., Lindberg, F., & Grimmond, C. S. B. (2011). Global to city scale urban anthropogenic
  ! heat flux: Model and variability. International Journal of Climatology, 31(13),
  ! 1990–2005. https://doi.org/10.1002/joc.2210
  !
  ! -----------------------------------------------------------------------

   IMPLICIT NONE

   INTEGER , intent(in) :: &
      idate(3)   ! calendar (year, julian day, seconds)

   REAL(r8), intent(in) :: &
      fix_holiday(365), &! Fixed public holidays, holiday(0) or workday(1)
      week_holiday(7)    ! week holidays

   REAL(r8), intent(in) :: &
      deltim      , &! seconds in a time step [second]
      patchlonr   , &! longitude of patch [radian]
      hum_prof(24), &! Diurnal metabolic heat profile [W/person]
      wdh_prof(24), &! Diurnal traffic flow profile of weekday
      weh_prof(24), &! Diurnal traffic flow profile of weekend
      popcell     , &! population density [person per square kilometer]
      vehicle(3)     ! vehicle numbers per thousand people

   REAL(r8) :: &
      vehc_prof(24,2), &
      carscell, &! cars numbers per thousand people
      frescell, &! freights numbers per thousand people
      mbkscell   ! motobikes numbers per thousand people

   REAL(r8), intent(out) :: &
      Fahe, &! flux from metabolic and vehicle
      vehc, &! flux from vehicle
      meta   ! flux from metabolic

   REAL(r8) ::  &
      londeg   , &! longitude of path [degree]
      car_sp   , &! distance traveled [km]
      traf_frac, &! vehicle heat profile of hour [-]
      meta_prof, &! metabolic heat profile of hour [-]
      carflx   , &! flux from car [W/m2]
      motflx   , &! flux from motorbike [W/m2]
      freflx      ! flux from freight [W/m2]


   ! local vars
   INTEGER :: &
         ldate(3), &! local time (year, julian day, seconds)
         iweek   , &! day of week
         ihour   , &! hour of day
         day     , &! day of mmonth
         month   , &! month of year
         day_inx , &! holiday index, day=1(workday), day=1(holiday)
         EC      , &! emission factor of car [J/m]
         EF      , &! emission factor of freight [J/m]
         EM         ! emission factor of motorbike [J/m]

   ! set vehicle distance traveled
   car_sp = 50

   ! emission factor Sailor and Lu (2004),
   ! all vehicle are set to same value
   EC = 3975
   EM = 3975
   EF = 3975

#ifndef USE_POINT_DATA
   ! convert GMT time to local time
   londeg = patchlonr*180/PI
   CALL gmt2local(idate, londeg, ldate)
#endif

   vehc_prof(:,1) = wdh_prof
   vehc_prof(:,2) = weh_prof

   CALL julian2monthday(ldate(1), ldate(2), month, day)
   CALL timeweek(ldate(1), month, day, iweek)

   ihour = CEILING(ldate(3)*1./3600)

   IF (day==366)  day=365
   IF (fix_holiday(day)==0 .or. week_holiday(iweek)==0) THEN
      day_inx = 1
   ELSE
      day_inx = 2
   ENDIF

   ! set traffic flow to be used of this time step
   traf_frac = vehc_prof(ihour,day_inx)
   ! set heat release per people of this time step
   meta_prof = hum_prof (ihour)

   carscell = vehicle(1)
   mbkscell = vehicle(2)
   frescell = vehicle(3)

   ! heat release of metabolism [W/m2]
   meta = popcell*meta_prof/1e6
   ! heat release of cars [W/m2]
   IF (carscell > 0) THEN
      carflx = carscell*popcell/1000
      carflx = carflx*traf_frac &
               *EC*(car_sp*1000)/1e6
      carflx = carflx/3600
   ENDIF

   ! heat release of motorbikes [W/m2]
   IF (mbkscell > 0) THEN
      motflx = mbkscell*popcell/1000
      motflx = motflx*traf_frac &
               *EM*(car_sp*1000)/1e6
      motflx   = motflx/3600
   ENDIF

   ! heat release of freight [W/m2]
   IF (frescell > 0)THEN
      freflx = frescell*popcell/1000
      freflx = freflx*traf_frac &
               *EF*(car_sp*1000)/1e6
      freflx = freflx/3600
   ENDIF

   ! total vehicle heat flux
   vehc = carflx + motflx + freflx
   ! total anthropogenic heat flux exclude building part
   Fahe = meta + vehc

  END Subroutine LUCY

  !TODO: write the below to timemanager.F90 @Wenzong
  ! -----------------------------------------------------------------------
  SUBROUTINE gmt2local(idate, long, ldate)

  ! !DESCRIPTION:
  ! A subroutine to calculate local time
  ! !PURPOSE
  ! Convert GMT time to local time in global run
  ! -----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER , intent(in ) :: idate(3)
    REAL(r8), intent(in ) :: long
    INTEGER , intent(out) :: ldate(3)

    INTEGER  :: maxday
    REAL(r8) :: tdiff

    tdiff = long/15.*3600

    ldate(3) = idate(3) + tdiff

    IF (ldate(3) < 0) THEN

      ldate(3) = 86400 + ldate(3)
      ldate(2) = idate(2) - 1

      IF (ldate(2) < 1) THEN
         ldate(1) = idate(1) - 1
         IF ( isleapyear(ldate(1)) ) THEN
            ldate(2) = 366
         ELSE
            ldate(2) = 365
         ENDIF
      ENDIF

    ELSE IF (ldate(3) > 86400) THEN

      ldate(3) = ldate(3) - 86400
      ldate(2) = idate(2) + 1

      IF ( isleapyear(ldate(1)) ) THEN
         maxday = 366
      ELSE
         maxday = 365
      ENDIF

      IF(ldate(2) > maxday) THEN
         ldate(1) = idate(1) + 1
         ldate(2) = 1
      ENDIF
    ELSE
      ldate(2) = idate(2)
      ldate(1) = idate(1)
    ENDIF
   END SUBROUTINE gmt2local

  ! -----------------------------------------------------------------------
  SUBROUTINE timeweek(year, month, day, iweek)

  ! !DESCRIPTION:
  ! A subroutine to calculate day of week
  ! !PURPOSE
  ! Calculate day of week to determine if the day is week holiday
  ! -----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, intent(in ) :: year, month
    INTEGER, intent(out) :: iweek, day

    INTEGER :: myear, mmonth
    INTEGER :: judy1, judy2, judy3
    INTEGER :: yy, mm, dd, y12, y34
    INTEGER :: A, B, C, D, i

    INTEGER, save :: DayOfMonth(12) = [31,28,31,30,31,30,31,31,30,31,30,31]

    judy1 = mod(year, 400)
    judy2 = mod(year, 100)
    judy3 = mod(year, 4  )

    IF (judy2 == 0) THEN
      IF (judy1 == 0) THEN
         DayOfMonth(2) = 29
      ELSE
         DayOfMonth(2) = 28
      ENDIF
    ELSE
      IF (judy3 == 0) THEN
         DayOfMonth(2) = 29
      ELSE
         DayOfMonth(2) = 28
      ENDIF
    ENDIF

    IF (month==1 .or. month==2) THEN
      mmonth = month + 12
      myear  = year  - 1
    ENDIF

    y12 = myear/100
    y34 = myear - y12*100

    A = int(y34/4.)
    B = int(y12/4.)
    C = y12*2
    D = int(26*(mmonth+1)/10.)

    iweek = abs(mod((y34+A+B-C+D+day-1), 7))

    DO i=1, month-1
       day = day + DayOfMonth(i)
    ENDDO

    IF (iweek == 0) THEN
       iweek = 7
    ENDIF

  END SUBROUTINE timeweek

END MODULE MOD_Urban_LUCY
