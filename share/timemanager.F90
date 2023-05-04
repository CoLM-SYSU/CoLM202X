#include <define.h>

! --------------------------------------------------------
MODULE timemanager

!
! !DESCRIPTION:
! Time manager module: to provide some basic operations for time stamp
!
! Created by Hua Yuan, 04/2014
!
! REVISIONS:
! 06/28/2017, Hua Yuan: added issame() and monthday2julian()
! TODO...
! --------------------------------------------------------


   USE precision
   IMPLICIT NONE

   TYPE :: timestamp
      INTEGER :: year, day, sec
   END TYPE timestamp

   INTERFACE ASSIGNMENT (=)
      MODULE procedure assignidate
      MODULE procedure assigntstamp
   END INTERFACE

   INTERFACE OPERATOR (+)
      MODULE procedure addsec
   END INTERFACE

   INTERFACE OPERATOR (-)
      MODULE procedure subtstamp
   END INTERFACE

   INTERFACE OPERATOR (<=)
      MODULE procedure lessequal
   END INTERFACE

   INTERFACE OPERATOR (<)
      MODULE procedure lessthan
   END INTERFACE

   INTERFACE OPERATOR (==)
      MODULE procedure isnull
      MODULE procedure besame
   END INTERFACE

   INTERFACE calendarday
      MODULE procedure calendarday_date
      MODULE procedure calendarday_stamp
   END INTERFACE

   LOGICAL, SAVE :: isgreenwich
   public get_calday

CONTAINS

   SUBROUTINE initimetype(greenwich)

      IMPLICIT NONE
      LOGICAL, intent(in) :: greenwich

      isgreenwich = greenwich

#ifndef SinglePoint
      IF (.not. isgreenwich) THEN
         write(*,*) 'Warning: Please Use Greenwich time for non-SinglePoint case.'
         isgreenwich = .true.
      ENDIF
#endif

   END SUBROUTINE initimetype

   SUBROUTINE assignidate(tstamp, idate)

      IMPLICIT NONE
      TYPE(timestamp), intent(inout) :: tstamp
      INTEGER,         intent(in)    :: idate(3)

      tstamp%year = idate(1)
      tstamp%day  = idate(2)
      tstamp%sec  = idate(3)

   END SUBROUTINE assignidate

   SUBROUTINE assigntstamp(tstamp1, tstamp2)

      IMPLICIT NONE
      TYPE(timestamp), intent(out) :: tstamp1
      TYPE(timestamp), intent(in)  :: tstamp2

      tstamp1%year = tstamp2%year
      tstamp1%day  = tstamp2%day
      tstamp1%sec  = tstamp2%sec

   END SUBROUTINE assigntstamp

   FUNCTION addsec(tstamp, sec)

      IMPLICIT NONE
      TYPE(timestamp), intent(in) :: tstamp
      INTEGER,         intent(in) :: sec
      TYPE(timestamp) :: addsec
      INTEGER         :: maxday

      addsec = tstamp
      addsec%sec = addsec%sec + sec
      IF (addsec%sec > 86400) THEN
         addsec%sec = addsec%sec - 86400
         IF( isleapyear(addsec%year) ) THEN
            maxday = 366
         ELSE
            maxday = 365
         ENDIF
         addsec%day = addsec%day + 1
         IF(addsec%day > maxday) THEN
            addsec%year = addsec%year + 1
            addsec%day = 1
         ENDIF
      ENDIF
      RETURN

   END FUNCTION addsec

   FUNCTION subtstamp(tstamp1, tstamp2)

      IMPLICIT NONE
      TYPE(timestamp), intent(in) :: tstamp1
      TYPE(timestamp), intent(in) :: tstamp2
      INTEGER :: subtstamp

      subtstamp = tstamp1%sec - tstamp2%sec
      IF (subtstamp < 0) THEN
         subtstamp = subtstamp + 86400
      ENDIF
      RETURN

   END FUNCTION subtstamp

   LOGICAL FUNCTION lessequal(tstamp1, tstamp2)

      IMPLICIT NONE
      TYPE(timestamp), intent(in) :: tstamp1
      TYPE(timestamp), intent(in) :: tstamp2

      INTEGER(kind=4) :: ts1, ts2

      ts1 = tstamp1%year*1000 + tstamp1%day
      ts2 = tstamp2%year*1000 + tstamp2%day

      lessequal = .false.

      IF (ts1 < ts2) lessequal = .true.

      IF (ts1==ts2 .AND. tstamp1%sec<=tstamp2%sec) THEN
         lessequal = .true.
      ENDIF

      RETURN

   END FUNCTION lessequal

   LOGICAL FUNCTION lessthan(tstamp1, tstamp2)

      IMPLICIT NONE
      TYPE(timestamp), intent(in) :: tstamp1
      TYPE(timestamp), intent(in) :: tstamp2

      INTEGER(kind=4) :: ts1, ts2

      ts1 = tstamp1%year*1000 + tstamp1%day
      ts2 = tstamp2%year*1000 + tstamp2%day

      lessthan = .false.

      IF (ts1 < ts2) lessthan = .true.

      IF (ts1==ts2 .AND. tstamp1%sec<tstamp2%sec) THEN
         lessthan = .true.
      ENDIF

      RETURN

   END FUNCTION lessthan

   LOGICAL FUNCTION isnull(tstamp, nullstr)

      IMPLICIT NONE
      TYPE(timestamp), intent(in) :: tstamp
      CHARACTER(4),    intent(in) :: nullstr

      IF (tstamp%year < 0 .OR. tstamp%day < 0 .OR. tstamp%sec < 0) THEN
         isnull = .true.
      ELSE
         isnull = .false.
      ENDIF
      RETURN

   END FUNCTION isnull

   LOGICAL FUNCTION besame(tstamp1, tstamp2)

      IMPLICIT NONE
      TYPE(timestamp), intent(in) :: tstamp1
      TYPE(timestamp), intent(in) :: tstamp2

      INTEGER :: idate1(3), idate2(3)

      idate1(1) = tstamp1%year
      idate1(2) = tstamp1%day
      idate1(3) = tstamp1%sec
      idate2(1) = tstamp2%year
      idate2(2) = tstamp2%day
      idate2(3) = tstamp2%sec

      CALL adj2end(idate1)
      CALL adj2end(idate2)

      IF (idate1(1)==idate2(1) .and. &
          idate1(2)==idate2(2) .and. &
          idate1(3)==idate2(3)) THEN
         besame = .true.
      ELSE
         besame = .false.
      ENDIF
      RETURN

   END FUNCTION besame

   LOGICAL FUNCTION isleapyear(year)

      IMPLICIT NONE
      INTEGER, intent(in) :: year

      IF( (mod(year,4)==0 .AND. mod(year,100)/=0) .OR. &
         mod(year,400)==0 ) THEN
         isleapyear = .true.
      ELSE
         isleapyear = .false.
      ENDIF
      RETURN
   END FUNCTION isleapyear

   SUBROUTINE julian2monthday(year, day, month, mday)

      IMPLICIT NONE
      INTEGER, intent(in)  :: year, day
      INTEGER, intent(out) :: month, mday

      INTEGER :: i, months(0:12)

      IF ( isleapyear(year) ) THEN
         months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
      ELSE
         months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      ENDIF

    ! calculate month and day values
      DO i = 1, 12
         IF (day .LE. months(i)) THEN
            month = i; exit
         ENDIF
      ENDDO
      mday = day - months(i-1)

   END SUBROUTINE julian2monthday

   SUBROUTINE monthday2julian(year, month, mday, day)

      IMPLICIT NONE
      INTEGER, intent(in)  :: year, month, mday
      INTEGER, intent(out) :: day

      INTEGER :: i, months(0:12)

      IF ( isleapyear(year) ) THEN
         months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
      ELSE
         months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      ENDIF

    ! calculate julian day
      day  = months(month-1) + mday

   END SUBROUTINE monthday2julian

   LOGICAL FUNCTION isendofhour(idate, sec)

      IMPLICIT NONE
      INTEGER, intent(in) :: idate(3)
      REAL(r8),intent(in) :: sec

      INTEGER :: hour1, hour2

      hour1 = (idate(3)-1)/3600
      hour2 = (idate(3)+int(sec)-1)/3600

      isendofhour = (hour1 /= hour2)

   END FUNCTION isendofhour

   LOGICAL FUNCTION isendofday(idate, sec)

      IMPLICIT NONE
      INTEGER, intent(in) :: idate(3)
      REAL(r8),intent(in) :: sec

      TYPE(timestamp) :: tstamp1
      TYPE(timestamp) :: tstamp2

      tstamp1 = idate
      tstamp2 = tstamp1 + int(sec)

      IF (tstamp2%day /= tstamp1%day) THEN
         isendofday = .true.
      ELSE
         isendofday = .false.
      ENDIF
      RETURN

   END FUNCTION isendofday

   LOGICAL FUNCTION isendofmonth(idate, sec)

      IMPLICIT NONE
      INTEGER, intent(in) :: idate(3)
      REAL(r8),intent(in) :: sec

      TYPE(timestamp) :: tstamp1
      TYPE(timestamp) :: tstamp2
      INTEGER :: month1, month2, day

      tstamp1 = idate
      tstamp2 = tstamp1 + int(sec)

      CALL julian2monthday(tstamp1%year, tstamp1%day, month1, day)
      CALL julian2monthday(tstamp2%year, tstamp2%day, month2, day)

      IF (month1 /= month2) THEN
         isendofmonth = .true.
      ELSE
         isendofmonth = .false.
      ENDIF
      RETURN

   END FUNCTION isendofmonth

   LOGICAL FUNCTION isendofyear(idate, sec)

      IMPLICIT NONE
      INTEGER, intent(in) :: idate(3)
      REAL(r8),intent(in) :: sec

      TYPE(timestamp) :: tstamp1
      TYPE(timestamp) :: tstamp2
      INTEGER :: month1, month2, day

      tstamp1 = idate
      tstamp2 = tstamp1 + int(sec)

      IF (tstamp1%year /= tstamp2%year) THEN
         isendofyear = .true.
      ELSE
         isendofyear = .false.
      ENDIF
      RETURN

   END FUNCTION isendofyear

   SUBROUTINE adj2begin(idate)

      IMPLICIT NONE
      INTEGER, intent(inout) :: idate(3)

      IF (idate(3) == 86400) THEN
         idate(3) = 0
         idate(2) = idate(2) + 1
         IF (isleapyear(idate(1)) .AND. idate(2)==367) THEN
            idate(1) = idate(1) + 1; idate(2) = 1
         ENDIF
         IF ( .NOT. isleapyear(idate(1)) .AND. idate(2)==366) THEN
            idate(1) = idate(1) + 1; idate(2) = 1
         ENDIF
      ENDIF

   END SUBROUTINE adj2begin

   SUBROUTINE adj2end(idate)

      IMPLICIT NONE
      INTEGER, intent(inout) :: idate(3)

      IF (idate(3) == 0) THEN
         idate(3) = 86400
         idate(2) = idate(2) - 1
         IF (idate(2) == 0) THEN
            idate(1) = idate(1) - 1
            IF ( isleapyear(idate(1)) ) THEN
               idate(2) = 366
            ELSE
               idate(2) = 365
            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE adj2end

   SUBROUTINE localtime2gmt(idate, long)

      IMPLICIT NONE
      INTEGER, intent(inout) :: idate(3)
      REAL(r8),intent(in)    :: long

      INTEGER  maxday
      REAL(r8) tdiff

      tdiff = long/15.*3600.
      idate(3) = idate(3) - int(tdiff)

      IF (idate(3) < 0) THEN

         idate(3) = 86400 + idate(3)
         idate(2) = idate(2) - 1

         IF (idate(2) < 1) THEN
            idate(1) = idate(1) - 1
            IF ( isleapyear(idate(1)) ) THEN
               idate(2) = 366
            ELSE
               idate(2) = 365
            ENDIF
         ENDIF
      ENDIF

      IF (idate(3) > 86400) THEN

         idate(3) = idate(3) - 86400
         idate(2) = idate(2) + 1

         IF ( isleapyear(idate(1)) ) THEN
            maxday = 366
         ELSE
            maxday = 365
         ENDIF

         IF(idate(2) > maxday) THEN
            idate(1) = idate(1) + 1
            idate(2) = 1
         ENDIF
      ENDIF

   END SUBROUTINE localtime2gmt

   SUBROUTINE ticktime(deltim, idate)

      IMPLICIT NONE

      REAL(r8),INTENT(in)    :: deltim
      INTEGER, INTENT(inout) :: idate(3)
      INTEGER maxday

      idate(3) = idate(3) + nint(deltim)
      IF (idate(3) > 86400) THEN

         idate(3) = idate(3) - 86400
         idate(2) = idate(2) + 1

         IF ( isleapyear(idate(1)) ) THEN
            maxday = 366
         ELSE
            maxday = 365
         ENDIF

         IF(idate(2) > maxday) THEN
            idate(1) = idate(1) + 1
            idate(2) = 1
         ENDIF
      ENDIF

   END SUBROUTINE ticktime

   REAL(r8) FUNCTION calendarday_date(date, long)

      IMPLICIT NONE
      INTEGER, intent(in) :: date(3)
      REAL(r8),optional   :: long

      INTEGER idate(3)
      REAL(r8) longitude

      idate(:) = date(:)

      IF (.NOT. present(long)) THEN
         longitude = 0._r8
      ELSE
         longitude = long
      ENDIF

      IF ( .not. isgreenwich ) THEN
         CALL localtime2gmt(idate, longitude)
      ENDIF

      calendarday_date = float(idate(2)) + float(idate(3))/86400.
      RETURN

   END FUNCTION calendarday_date

   REAL(r8) FUNCTION calendarday_stamp(stamp, long)

      IMPLICIT NONE
      TYPE(timestamp), intent(in) :: stamp
      REAL(r8),        optional   :: long

      INTEGER idate(3)
      REAL(r8) longitude

      idate(1) = stamp%year
      idate(2) = stamp%day
      idate(3) = stamp%sec

      IF (.NOT. present(long)) THEN
         longitude = 0._r8
      ELSE
         longitude = long
      ENDIF
      IF ( .not. isgreenwich ) THEN
         CALL localtime2gmt(idate, longitude)
      ENDIF

      calendarday_stamp = float(idate(2)) + float(idate(3))/86400.
      RETURN

   END FUNCTION calendarday_stamp

   INTEGER FUNCTION get_calday(mmdd,isleap)

      IMPLICIT NONE
      INTEGER, intent(in) :: mmdd
      LOGICAL, intent(in) :: isleap

      INTEGER, dimension(0:12) :: daysofmonth_noleap &
                               = (/0,31,28,31,30,31,30,31,31,30,31,30,31/)
      INTEGER, dimension(0:12) :: daysofmonth_leap &
                               = (/0,31,29,31,30,31,30,31,31,30,31,30,31/)
      INTEGER imonth, iday

      imonth = mmdd / 100
      iday   = mod(mmdd,100)
      IF(isleap)THEN
         get_calday = sum(daysofmonth_leap(0:imonth-1)) + iday
      ELSE
         get_calday = sum(daysofmonth_noleap(0:imonth-1)) + iday
      END IF
      RETURN
   END FUNCTION get_calday

   INTEGER FUNCTION minutes_since_1900 (year, julianday, second)

      IMPLICIT NONE
      INTEGER, intent(in) :: year, julianday, second

      INTEGER :: minutes_allyear = 8760
      INTEGER :: refyear(10) = (/1, 1900, 1950, 1980, 1990, 2000, 2005, 2010, 2015, 2020/)
      INTEGER :: refval (10) = (/-998776800,0,26297280,42075360,47335680,52594560,55225440,&
                                 57854880,60484320,63113760/)
      INTEGER :: iref, iyear

      iref = findloc(refyear <= year, .true., back=.true., dim=1)
      minutes_since_1900 = refval(iref)
      DO iyear = refyear(iref), year-1
         IF (isleapyear(iyear)) THEN
            minutes_since_1900 = minutes_since_1900 + 527040
         ELSE
            minutes_since_1900 = minutes_since_1900 + 525600
         ENDIF
      ENDDO

      minutes_since_1900 = minutes_since_1900 + (julianday-1) * 1440
      minutes_since_1900 = minutes_since_1900 + second/60

   END FUNCTION minutes_since_1900


END MODULE timemanager
