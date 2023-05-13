      program calcnewdate
C     ===================

      implicit none

      integer   date1(5),date2(5)
      integer	iargc
      real      diff
      character*(80) arg,yychar,cdat
      character*(2)  cdate(4)
      integer	flag,nc

c     check for sufficient requested arguments
      if (iargc().ne.2) then
        print*,'USAGE: newtime date (format (YY)YYMMDD_HH) timestep'
        call exit(1)
      endif

c     read and transform input
      call getarg(1,arg)
      call lenchar(arg,nc)
      if (nc.eq.9) then
        yychar=''
        read(arg(1:2),'(i2)',err=120)date1(1)
        read(arg(3:4),'(i2)',err=120)date1(2)
        read(arg(5:6),'(i2)',err=120)date1(3)
        read(arg(8:9),'(i2)',err=120)date1(4)
      else if (nc.eq.11) then
        yychar=arg(1:2)
        read(arg(3:4),'(i2)',err=120)date1(1)
        read(arg(5:6),'(i2)',err=120)date1(2)
        read(arg(7:8),'(i2)',err=120)date1(3)
        read(arg(10:11),'(i2)',err=120)date1(4)
      else
        print*,'USAGE: newtime date (format (YY)YYMMDD_HH) timestep'
        call exit(1)
      endif

      call getarg(2,arg)
      call checkchar(arg,".",flag)
      if (flag.eq.0) arg=trim(arg)//"."
      read(arg,'(f10.2)') diff

      call newdate(date1,diff,date2)

      if ((date2(1).lt.date1(1)).and.
     >    (diff.gt.0.).and.
     >    (yychar.eq.'19')) yychar='20'

c      if ((date2(1).lt.date1(1)).and.
c     >    (diff.gt.0.).and.
c     >    (yychar.eq.'20')) yychar='21'

      if (date2(1).lt.0) date2(1)=date2(1)+100

      if ((date2(1).gt.date1(1)).and.
     >    (diff.lt.0.).and.
     >    (yychar.eq.'20')) yychar='19'

c      if ((date2(1).gt.date1(1)).and.
c     >    (diff.lt.0.).and.
c     >    (yychar.eq.'21')) yychar='20'


      if (date2(1).lt.10) then
        write(cdate(1),'(a,i1)')'0',date2(1)
      else
        write(cdate(1),'(i2)')date2(1)
      endif
      if (date2(2).lt.10) then
        write(cdate(2),'(a,i1)')'0',date2(2)
      else
        write(cdate(2),'(i2)')date2(2)
      endif
      if (date2(3).lt.10) then
        write(cdate(3),'(a,i1)')'0',date2(3)
      else
        write(cdate(3),'(i2)')date2(3)
      endif
      if (date2(4).lt.10) then
        write(cdate(4),'(a,i1)')'0',date2(4)
      else
        write(cdate(4),'(i2)')date2(4)
      endif

      cdat=trim(yychar)//cdate(1)//cdate(2)//cdate(3)//'_'//cdate(4)
      write(*,'(a)')trim(cdat)

      goto 200

 120  write(*,*)"*** error: date must be in format (YY)YYMMDD_HH ***"

 200  continue

      end

      subroutine checkchar(string,char,flag)
C     ======================================

      character*(*)	string
      character*(1)	char
      integer	n,flag

      flag=0
      do n=1,len(string)
        if (string(n:n).eq.char) then
          flag=n
          return
        endif
      enddo
      end

      subroutine lenchar(string,lstr)
C     ===============================

      character*(*)     string
      integer   n,lstr

      do n=1,len(string)
        if (string(n:n).eq."") then
          lstr=n-1
          goto 100
        endif
      enddo
 100  continue
      end
