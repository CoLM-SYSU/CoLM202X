c     **************************************************************************
c     * This library provides subroutines related to time and file name        *
c     * handling                                                               *
c     *                                                                        *
c     **************************************************************************

c     ---------------------------------------------------------------------------
c     Concatenate a date string
c     ---------------------------------------------------------------------------

      subroutine datestring(datestr,yyyy,mm,dd,hh)

c     Declaration of subroutine parameters
      integer yyyy,mm,dd,hh
      character*20 datestr

c     Auxiliary parameters
      integer f1,f2,i0
      integer yy,ce

      i0=ichar('0')
      datestr=''

      yy=mod(yyyy,100)
      ce=int(yyyy/100)

      if ((ce.ne.19).and.(ce.ne.20)) then
         print*,'Invalid year... Stop'
         stop
      endif

      if (yy.ge.0) then
         f1=yy/10
         f2=mod(yy,10)
         if (ce.eq.19) then
            datestr=trim(datestr)//'19'//char(f1+i0)//char(f2+i0)
         else if (ce.eq.20) then
            datestr=trim(datestr)//'20'//char(f1+i0)//char(f2+i0)
         endif
      endif
      if (mm.gt.0) then
         f1=mm/10
         f2=mod(mm,10)
         datestr=trim(datestr)//char(f1+i0)//char(f2+i0)
      endif

      if (dd.gt.0) then
         f1=dd/10
         f2=mod(dd,10)
         datestr=trim(datestr)//char(f1+i0)//char(f2+i0)
      endif

      if (hh.ge.0) then
         f1=hh/10
         f2=mod(hh,10)
         datestr=trim(datestr)//'_'//char(f1+i0)//char(f2+i0)
      endif

      return
      end


c     ---------------------------------------------------------------------------
c     Calculates the new date when diff (in hours) is added to date1.
c     ---------------------------------------------------------------------------

      subroutine newdate(date1,diff,date2)
C
C     date1	int	input	array contains a date in the form
C				year,month,day,hour,step
C     diff	real	input	timestep in hours to go from date1
C     date2	int	output	array contains new date in the same form

      integer   date1(5),date2(5)
      integer   idays(12)       ! array containing the days of the monthes
      real	diff
      logical	yearchange

      data idays/31,28,31,30,31,30,31,31,30,31,30,31/

      yearchange=.false.

      if ((mod(date1(1),4).eq.0).and.(date1(2).le.2)) idays(2)=29

      date2(1)=date1(1)
      date2(2)=date1(2)
      date2(3)=date1(3)
      date2(4)=date1(4)
      date2(5)=0
      date2(4)=date1(4)+int(diff)+date1(5)

      if (date2(4).ge.24) then
        date2(3)=date2(3)+int(date2(4)/24)
        date2(4)=date2(4)-int(date2(4)/24)*24
      endif
      if (date2(4).lt.0) then
        if (mod(date2(4),24).eq.0) then
          date2(3)=date2(3)-int(abs(date2(4))/24)
          date2(4)=date2(4)+int(abs(date2(4))/24)*24
        else
          date2(3)=date2(3)-(1+int(abs(date2(4))/24))
          date2(4)=date2(4)+(1+int(abs(date2(4))/24))*24
        endif
      endif

  100 if (date2(3).gt.idays(date2(2))) then
        if ((date2(2).eq.2).and.(mod(date2(1),4).eq.0)) idays(2)=29
        date2(3)=date2(3)-idays(date2(2))
        if (idays(2).eq.29) idays(2)=28
        date2(2)=date2(2)+1
        if (date2(2).gt.12) then
*         date2(1)=date2(1)+int(date2(2)/12)
*         date2(2)=date2(2)-int(date2(2)/12)*12
          date2(1)=date2(1)+1
          date2(2)=date2(2)-12
        endif
        if (date2(2).lt.1) then
          date2(1)=date2(1)-(1+int(abs(date2(2)/12)))
          date2(2)=date2(2)+(1+int(abs(date2(2)/12)))*12
        endif
        goto 100
      endif     
  200 if (date2(3).lt.1) then
        date2(2)=date2(2)-1
        if (date2(2).gt.12) then
          date2(1)=date2(1)+int(date2(2)/12)
          date2(2)=date2(2)-int(date2(2)/12)*12
        endif
        if (date2(2).lt.1) then
          date2(1)=date2(1)-(1+int(abs(date2(2)/12)))
          date2(2)=date2(2)+(1+int(abs(date2(2)/12)))*12
        endif
        if ((date2(2).eq.2).and.(mod(date2(1),4).eq.0)) idays(2)=29
        date2(3)=date2(3)+idays(date2(2))
        if (idays(2).eq.29) idays(2)=28
        goto 200
      endif

      if (date2(2).gt.12) then
        date2(1)=date2(1)+int(date2(2)/12)
        date2(2)=date2(2)-int(date2(2)/12)*12
      endif
      if (date2(2).lt.1) then
        date2(1)=date2(1)-(1+int(abs(date2(2)/12)))
        date2(2)=date2(2)+(1+int(abs(date2(2)/12)))*12
      endif

      if (date2(1).lt.1000) then
      if (date2(1).ge.100) date2(1)=date2(1)-100
      endif

      return
      end

c     ---------------------------------------------------------------------------
c     Convert <hh.mm> -> <frac> and <frac> -> <hhmm>
c     ---------------------------------------------------------------------------

      subroutine hhmm2frac (hhmm,frac)

      real hhmm
      real frac
      
      frac=real(int(hhmm))+
     >     100.*(hhmm-real(int(hhmm)))/60.

      end


      subroutine frac2hhmm (frac,hhmm)

      real hhmm
      real frac

      real hh,mm
      
      hh =  real(int(frac))
      mm =  60. * (frac-real(int(frac)))

      hhmm = hh + 0.01 * mm

      end

c     ---------------------------------------------------------------------------
c     Time difference between two dates
c     ---------------------------------------------------------------------------

      subroutine timediff(date1,date2,diff)
 
C     New version with hour and minutes! (for hour and step [in hours]
C     use the routine oldtimediff!)
C
C     Calculates the time difference in hours (and minutes) for the two
C     dates specified by the two arrays date1 and date2.
C     They are expected to contain the following date information:
C     year      month   day     hour    minute.
C
C     date1     array specifying the first date
C     date2     array specifying the second date
C     diff      time differenc between date1 and date2 in hours
C
 
      integer   date1(5),date2(5)
      integer   idays(12)       ! array containing the days of the monthes
      real      diff
      integer   ixday,imdiff,ihdiff,iddiff,j
      integer   yy,yy1,yy2
 
      idays(1)=31
      idays(2)=28
      idays(3)=31
      idays(4)=30
      idays(5)=31
      idays(6)=30
      idays(7)=31
      idays(8)=31
      idays(9)=30
      idays(10)=31
      idays(11)=30
      idays(12)=31
 
C     Check format of year (YYYY or YY - in case of YY assume 19YY)

      if (date1(1).lt.100) date1(1)=1900+date1(1)
      if (date2(1).lt.100) date2(1)=1900+date2(1)

C     Determine if the period between date1 and date2 contains a Feb.29
 
      ixday=0   ! extra day flag
 
      yy1=min(date1(1),date2(1))
      yy2=max(date1(1),date2(1))
      if (yy1.eq.yy2) then
        if (mod(yy1,4).eq.0) then
          idays(2)=29
        endif
      else
        if (mod(yy1,4).eq.0) then
          if (((yy1.eq.date1(1)).and.(date1(2).le.2)).or.
     >        ((yy1.eq.date2(1)).and.(date2(2).le.2))) then
            ixday=ixday+1
          endif
        endif
        if (mod(yy2,4).eq.0) then
          if (((yy2.eq.date1(1)).and.(date1(2).gt.2)).or.
     >        ((yy2.eq.date2(1)).and.(date2(2).gt.2))) then
            ixday=ixday+1
          endif
        endif
        if (yy2-yy1.gt.1) then
          do yy=yy1+1,yy2-1
            if (mod(yy,4).eq.0) then
              ixday=ixday+1
            endif
          enddo
        endif
      endif
 
      ihdiff=0  ! diff. in hours between date1/date2
      iddiff=0  ! diff. in days  between date1/date2
 
      if (date1(1).gt.date2(1)) then            ! compare years
        do j=date2(1),date1(1)-1
          iddiff=iddiff+365
        enddo
        iddiff=iddiff+ixday
      else if (date1(1).lt.date2(1)) then
        do j=date1(1),date2(1)-1
          iddiff=iddiff-365
        enddo
        iddiff=iddiff-ixday
      endif
 
      if (date1(2).gt.date2(2)) then            ! compare monthes
        do j=date2(2),date1(2)-1
          iddiff=iddiff+idays(j)
        enddo
      else if (date1(2).lt.date2(2)) then
        do j=date1(2),date2(2)-1
          iddiff=iddiff-idays(j)
        enddo
      endif
 
      iddiff=iddiff+date1(3)-date2(3)
      ihdiff=iddiff*24+date1(4)-date2(4)
      imdiff=ihdiff*60+date1(5)-date2(5)
 
      ihdiff=imdiff/60
      imdiff=mod(imdiff,60)
 
      diff=real(ihdiff)+real(imdiff)/100.
 
      return
      end
