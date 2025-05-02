MODULE CMF_UTILS_MOD
!==========================================================
!* PURPOSE: Shared ulitity functions/subroutines for CaMa-Flood
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
   USE PARKIND1,                only: JPIM,   JPRB, JPRM, JPRD
   USE YOS_CMF_INPUT,           only: LOGNAM, DMIS, RMIS, NX,NY
   USE YOS_CMF_MAP,             only: NSEQMAX, NSEQALL
   IMPLICIT NONE
CONTAINS
!####################################################################
! map related subroutines & functions
!-- vecP2mapR     : convert 1D vector data -> 2D map data (real*4)
!-- vecD2mapD    : convert 1D vector data -> 2D map data (real*8)
!-- mapR2vecD     : convert 2D map data -> 1D vector data (REAL*4)
!-- mapP2vecP    : convert 2D map data -> 1D vector data (REAL*8)
!-- mapI2vecI    : convert 2D map data -> 1D vector data (Integer)
!
! time related subroutines & functions
! -- MIN2DATE  : calculate DATE of KMIN from base time (YYYY0,MM0,DD0)
! -- DATE2MIN  : convert (YYYYMMDD,HHMM) to KMIN from base time (YYYY0,MM0,DD0)
! -- SPLITDATE : splite date (YYYYMMDD) to (YYYY,MM,DD)
! -- SPLITHOUR : split hour (HHMM) to (HH,MM)
! -- IMDAYS    : function to calculate days in a monty IMDAYS(IYEAR,IMON)
!
! endian conversion
!-- CONV_END    : Convert 2D Array endian (REAL4)
!-- CONV_ENDI   : Convert 2D Array endian (Integer)
!-- ENDIAN4R    : byte swap (REAL*4)
!-- ENDIAN4I    : byte swap (Integer)
!
! file I/O
!-- INQUIRE_FID : inruire unused file FID
!-- NCERROR     : netCDF I/O wrapper
!-- CMF_CheckNaN: check the value is NaN or not
!####################################################################
   SUBROUTINE vecD2mapR(D2VEC,R2MAP)
   USE YOS_CMF_MAP,             only: I1SEQX,I1SEQY
   IMPLICIT NONE
   !* input/output
   real(KIND=JPRB),intent(in)      :: D2VEC(NSEQMAX,1)
   real(KIND=JPRM),intent(out)     :: R2MAP(NX,NY)
   !* local variable
   integer(KIND=JPIM),SAVE         ::  IX,IY,ISEQ
!$OMP THREADPRIVATE                (IX,IY)
!================================================
      R2MAP(:,:) = RMIS
!$OMP PARALLEL DO
      DO ISEQ=1,NSEQALL
         IX=I1SEQX(ISEQ)
         IY=I1SEQY(ISEQ)
         R2MAP(IX,IY) = real(D2VEC(ISEQ,1),KIND=JPRM)
      ENDDO
!$OMP END PARALLEL DO

   END SUBROUTINE vecD2mapR
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE vecD2mapD(D2VEC,D2MAP)
   USE YOS_CMF_MAP,             only: I1SEQX,I1SEQY
   IMPLICIT NONE
   !* input/output
   real(KIND=JPRB),intent(in)      :: D2VEC(NSEQMAX,1)
   real(KIND=JPRB),intent(out)     :: D2MAP(NX,NY)
   !* local variable
   integer(KIND=JPIM),SAVE         ::  IX,IY,ISEQ
!$OMP THREADPRIVATE                (IX,IY)
!================================================
      D2MAP(:,:) = DMIS
!$OMP PARALLEL DO
      DO ISEQ=1,NSEQALL
         IX=I1SEQX(ISEQ)
         IY=I1SEQY(ISEQ)
         D2MAP(IX,IY) = D2VEC(ISEQ,1)
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE vecD2mapD
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE vecP2mapP(P2VEC,P2MAP)
   USE YOS_CMF_MAP,             only: I1SEQX,I1SEQY
   IMPLICIT NONE
   !* input/output
   real(KIND=JPRD),intent(in)      :: P2VEC(NSEQMAX,1)
   real(KIND=JPRD),intent(out)     :: P2MAP(NX,NY)
   !* local variable
   integer(KIND=JPIM),SAVE         ::  IX,IY,ISEQ
!$OMP THREADPRIVATE                (IX,IY)
      !================================================
      P2MAP(:,:) = DMIS
!$OMP PARALLEL DO
      DO ISEQ=1,NSEQALL
         IX=I1SEQX(ISEQ)
         IY=I1SEQY(ISEQ)
         P2MAP(IX,IY) = P2VEC(ISEQ,1)
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE vecP2mapP
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE vecP2mapR(P2VEC,R2MAP)
   USE YOS_CMF_MAP,             only: I1SEQX,I1SEQY
   IMPLICIT NONE
   !* input/output
   real(KIND=JPRD),intent(in)      :: P2VEC(NSEQMAX,1)
   real(KIND=JPRM),intent(out)     :: R2MAP(NX,NY)
   !* local variable
   integer(KIND=JPIM),SAVE         ::  IX,IY,ISEQ
!$OMP THREADPRIVATE                (IX,IY)
      !================================================
      R2MAP(:,:) = RMIS
!$OMP PARALLEL DO
      DO ISEQ=1,NSEQALL
         IX=I1SEQX(ISEQ)
         IY=I1SEQY(ISEQ)
         R2MAP(IX,IY) = real(P2VEC(ISEQ,1),4)
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE vecP2mapR
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE mapR2vecD(R2TEMP,D2VAR)
   USE YOS_CMF_MAP,             only: I1SEQX,I1SEQY
   IMPLICIT NONE
   !* input/output
   real(KIND=JPRM),intent(in)      :: R2TEMP(NX,NY)
   real(KIND=JPRB),intent(out)     :: D2VAR(NSEQMAX,1)
   !* local variable
   integer(KIND=JPIM),SAVE         :: IX,IY, ISEQ
!$OMP THREADPRIVATE               (IX,IY)
!================================================
!$OMP PARALLEL DO
      DO ISEQ=1,NSEQALL
         IX=I1SEQX(ISEQ)
         IY=I1SEQY(ISEQ)
         D2VAR(ISEQ,1) = real(R2TEMP(IX,IY),KIND=JPRB)
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE mapR2vecD
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE mapD2vecD(D2TEMP,D2VAR)
   USE YOS_CMF_MAP,             only: I1SEQX,I1SEQY
   IMPLICIT NONE
   !* input/output
   real(KIND=JPRB),intent(in)      :: D2TEMP(NX,NY)
   real(KIND=JPRB),intent(out)     :: D2VAR(NSEQMAX,1)
   !* local variable
   integer(KIND=JPIM),SAVE         :: IX,IY, ISEQ
!$OMP THREADPRIVATE               (IX,IY)
!================================================
!$OMP PARALLEL DO
      DO ISEQ=1,NSEQALL
         IX=I1SEQX(ISEQ)
         IY=I1SEQY(ISEQ)
         D2VAR(ISEQ,1) = D2TEMP(IX,IY)
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE mapD2vecD
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE mapP2vecP(P2TEMP,P2VAR)
   USE YOS_CMF_MAP,             only: I1SEQX,I1SEQY
   IMPLICIT NONE
   !* input/output
   real(KIND=JPRD),intent(in)      :: P2TEMP(NX,NY)
   real(KIND=JPRD),intent(out)     :: P2VAR(NSEQMAX,1)
   !* local variable
   integer(KIND=JPIM),SAVE         :: IX,IY, ISEQ
!$OMP THREADPRIVATE               (IX,IY)
!================================================
!$OMP PARALLEL DO
      DO ISEQ=1,NSEQALL
         IX=I1SEQX(ISEQ)
         IY=I1SEQY(ISEQ)
         P2VAR(ISEQ,1) = P2TEMP(IX,IY)
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE mapP2vecP
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE mapP2vecD(P2TEMP,D2VAR)
   USE YOS_CMF_MAP,             only: I1SEQX,I1SEQY
   IMPLICIT NONE
   !* input/output
   real(KIND=JPRD),intent(in)      :: P2TEMP(NX,NY)
   real(KIND=JPRB),intent(out)     :: D2VAR(NSEQMAX,1)
   !* local variable
   integer(KIND=JPIM),SAVE         :: IX,IY, ISEQ
!$OMP THREADPRIVATE               (IX,IY)
!================================================
!$OMP PARALLEL DO
      DO ISEQ=1,NSEQALL
         IX=I1SEQX(ISEQ)
         IY=I1SEQY(ISEQ)
         D2VAR(ISEQ,1) = P2TEMP(IX,IY)
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE mapP2vecD
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE mapI2vecI(I2TEMP,I2VAR)
   USE YOS_CMF_MAP,             only: I1SEQX,I1SEQY
   IMPLICIT NONE
   !* input/output
   integer(KIND=JPIM),intent(in)   :: I2TEMP(NX,NY)
   integer(KIND=JPIM),intent(out)  :: I2VAR(NSEQMAX,1)
   !* local variable
   integer(KIND=JPIM),SAVE         :: IX,IY,ISEQ
!$OMP THREADPRIVATE               (IX,IY)
!================================================
!$OMP PARALLEL DO
      DO ISEQ=1,NSEQALL
         IX=I1SEQX(ISEQ)
         IY=I1SEQY(ISEQ)
         I2VAR(ISEQ,1) = I2TEMP(IX,IY)
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE mapI2vecI
!####################################################################





   !####################################################################
   ! time related subroutines & functions
   ! -- MIN2DATE  : calculate DATE of KMIN from base time (YYYY0,MM0,DD0)
   ! -- DATE2MIN  : convert (YYYYMMDD,HHMM) to KMIN from base time (YYYY0,MM0,DD0)
   ! -- SPLITDATE : splite date (YYYYMMDD) to (YYYY,MM,DD)
   ! -- SPLITHOUR : split hour (HHMM) to (HH,MM)
   ! -- IMDAYS    : function to calculate days in a monty IMDAYS(IYEAR,IMON)
   !==========================================================
   SUBROUTINE MIN2DATE(IMIN,YYYYMMDD,HHMM)
   !  Return YYYYMMDD and HHMM for IMIN
   USE YOS_CMF_TIME,            only: YYYY0, MM0, DD0
   IMPLICIT NONE
   ! local
   integer(KIND=JPIM),intent(in)   :: IMIN      !!  input minutes
   integer(KIND=JPIM),intent(out)  :: YYYYMMDD
   integer(KIND=JPIM),intent(out)  :: HHMM
   integer(KIND=JPIM)              :: YYYY,MM,DD,HH,MI,NDAYS,NDM,ID
   integer(KIND=JPIM)              :: D2MIN                   ! minutes in one day
   parameter                         (D2MIN=1440)
      !================================================
      YYYYMMDD = 0
      HHMM     = 0

      NDAYS = IMIN/D2MIN              !! days  in IMIN : 1440 = (minutes in a day)
      MI    = MOD(IMIN,D2MIN)
      HH    = INT(MI/60)              !! hours in IMIN
      MI    = MOD(MI,60)              !! mins  in IMIN

      YYYY  = YYYY0
      MM    = MM0
      DD    = DD0
      NDM   = IMDAYS(YYYY,MM)      !! number of days in a month

      ! write(LOGNAM,*)  YYYY,MM,DD
      DO ID=1,NDAYS
         DD=DD+1
         IF ( DD .gt. NDM ) THEN
            MM=MM+1
            DD=1
            IF ( MM .gt. 12 ) THEN
               MM=1
               YYYY=YYYY+1
            ENDIF
            NDM=IMDAYS(YYYY,MM)
         ENDIF
      ENDDO

      HHMM     = HH*100+MI
      YYYYMMDD = YYYY*10000+MM*100+DD
   END SUBROUTINE MIN2DATE
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   FUNCTION DATE2MIN(YYYYMMDD,HHMM)
   ! convert (YYYYMMDD,HHMM) to KMIN from base time (YYYY0,MM0,DD0)
   USE YOS_CMF_TIME,            only: YYYY0
   IMPLICIT NONE
   integer(KIND=JPIM)              :: DATE2MIN
   integer(KIND=JPIM),intent(in)   :: YYYYMMDD
   integer(KIND=JPIM),intent(in)   :: HHMM
   integer(KIND=JPIM)              :: YYYY,MM,DD,HH,MI
   integer(KIND=JPIM)              :: IY,IM
   integer(KIND=JPIM)              :: D2MIN                   ! minutes in one day
   parameter                         (D2MIN=1440)
      !================================================
      DATE2MIN = 0
      CALL SPLITDATE(YYYYMMDD,YYYY,MM,DD)
      HH = HHMM/100                          !! hour
      MI = HHMM-HH*100                       !! minute
      !============================
      IF ( YYYY .lt. YYYY0) THEN
         write(LOGNAM,*) 'DATE2MIN: YYYY .lt. YYYY0: Date Problem', YYYY,YYYY0
         STOP
      ENDIF
      IF ( MM.lt.1 .or. MM .gt. 12 ) THEN
         write(LOGNAM,*) 'DATE2MIN: MM:    Date Problem', YYYYMMDD, HHMM
         STOP
      ENDIF
      IF ( DD.lt.1 .or. DD .gt. IMDAYS(YYYY,MM)) THEN
         write(LOGNAM,*) 'DATE2MIN: DD:    Date Problem', YYYYMMDD, HHMM
         STOP
      ENDIF
      IF ( HH.lt.0 .or. HH .gt. 24) THEN
         write(LOGNAM,*) 'DATE2MIN: HH:    Date Problem', YYYYMMDD, HHMM
         STOP
      ENDIF
      IF ( MI.lt.0 .or. MI .gt. 60) THEN
         write(LOGNAM,*) 'DATE2MIN: MI:    Date Problem', YYYYMMDD, HHMM
         STOP
      ENDIF

      IY=YYYY0
      DO WHILE (IY .lt. YYYY)
         DO IM=1,12
            DATE2MIN=DATE2MIN+IMDAYS(IY,IM)*D2MIN
         ENDDO
         IY=IY+1
      ENDDO
      IM=1
      DO WHILE (IM .lt. MM )
         DATE2MIN=DATE2MIN+IMDAYS(IY,IM)*D2MIN
         IM=IM+1
      ENDDO

      DATE2MIN = DATE2MIN + (DD-1)*D2MIN
      DATE2MIN = DATE2MIN + HH*60 + MI

   END FUNCTION DATE2MIN
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE SPLITDATE(YYYYMMDD,YYYY,MM,DD)
   ! sprit YYYYMMDD to (YYYY,MM,DD)
   IMPLICIT NONE
   integer(KIND=JPIM),intent(in)   :: YYYYMMDD
   integer(KIND=JPIM),intent(out)  :: YYYY,MM,DD
      !================================================
      YYYY =  YYYYMMDD/10000
      MM   = (YYYYMMDD - YYYY*10000) / 100
      DD   =  YYYYMMDD -(YYYY*10000+MM*100)
   END SUBROUTINE SPLITDATE
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE SPLITHOUR(HHMM,HH,MI)
   ! sprit YYYYMMDD to (YYYY,MM,DD)
   IMPLICIT NONE
   integer(KIND=JPIM),intent(in)   :: HHMM
   integer(KIND=JPIM),intent(out)  :: HH,MI
      !================================================
      HH=INT(HHMM/100)
      MI=INT(HHMM-HH*100)
   END SUBROUTINE SPLITHOUR
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   FUNCTION IMDAYS(IYEAR,IMON)
   !! days in month
   USE YOS_CMF_INPUT,           only: LLEAPYR
   IMPLICIT NONE
   integer(KIND=JPIM)              :: IMDAYS
   integer(KIND=JPIM),intent(in)   :: IYEAR
   integer(KIND=JPIM),intent(in)   :: IMON
   integer(KIND=JPIM)              :: ND(12)
   DATA ND /31,28,31,30,31,30,31,31,30,31,30,31/
      !================================================
      IMDAYS=ND(IMON)
      IF ( IMON == 2 .and. LLEAPYR ) THEN
         IF ( MOD(IYEAR,400) == 0 .or. (MOD(IYEAR,100) .NE. 0 .and. MOD(IYEAR,4) .eq. 0 )) IMDAYS=29
      ENDIF
   END FUNCTION IMDAYS
!==========================================================




   !####################################################################
   ! endian conversion
   !-- CONV_END    : Convert 2D Array endian (REAL4)
   !-- CONV_ENDI   : Convert 2D Array endian (Integer)
   !-- ENDIAN4R    : byte swap (real*4)
   !-- ENDIAN4I    : byte swap (Integer)
   !####################################################################
   SUBROUTINE CONV_END(R2TEMP,NX,NY)
   !-- Convert 2D Array endian (REAL4)
   IMPLICIT NONE
   !* input/output
   integer(KIND=JPIM),intent(in)   :: NX,NY
   real(KIND=JPRM),intent(inout)   :: R2TEMP(NX,NY)
   !* local variables
   integer(KIND=JPIM)              :: IY,IX
   !================================================
      DO IY=1, NY
         DO IX=1, NY
            CALL ENDIAN4R(R2TEMP(IX,IY))
         ENDDO
      ENDDO
   END SUBROUTINE CONV_END
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE CONV_ENDI(I2TEMP,NX,NY)
   !-- Convert 2D Array endian (integer)
   IMPLICIT NONE
   !+ input/output
   integer(KIND=JPIM),intent(in)     :: NX,NY
   integer(KIND=JPIM),intent(inout)  :: I2TEMP(NX,NY)
   !* local variables
   integer(KIND=JPIM)                :: IY,IX
      !================================================
      DO IY=1, NY
         DO IX=1, NY
            CALL ENDIAN4I(I2TEMP(IX,IY))
         ENDDO
      ENDDO
   END SUBROUTINE CONV_ENDI
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE ENDIAN4R( realIn )
   !! Byte Swap 
   !
   ! Adpated from: http://www.cgd.ucar.edu/cas/software/endian.html
   !           FILE: SUBR_native_4byte_real.f90
   !     SUBPROGRAM: native_4byte_real
   !
   !         AUTHOR: David Stepaniak, NCAR/CGD/CAS
   ! DATE INITIATED: 29 April 2003 
   !  LAST MODIFIED: 29 April 2003
   IMPLICIT NONE
   !* input/output
   real(KIND=JPRM), intent(inout)  :: realIn
   !* Local variables (generic 32 bit integer spaces):
   integer                         :: i_element
   integer                         :: i_element_br
      !================================================
      ! Transfer 32 bits of realIn to generic 32 bit integer space:
      i_element_br=0
      i_element = TRANSFER( realIn, 0 )
      ! Reverse order of 4 bytes in 32 bit integer space:
      CALL MVBITS( i_element, 24, 8, i_element_br, 0  )
      CALL MVBITS( i_element, 16, 8, i_element_br, 8  )
      CALL MVBITS( i_element,  8, 8, i_element_br, 16 )
      CALL MVBITS( i_element,  0, 8, i_element_br, 24 )

      ! Transfer reversed order bytes to 32 bit real space (realOut):
      realIn = TRANSFER( i_element_br, 0.0 )
   END SUBROUTINE ENDIAN4R
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE ENDIAN4I(IntIn)
   !! Byte Swap 
   IMPLICIT NONE
   !* input/output
   integer(KIND=JPIM), intent(inout)    :: IntIn
   ! Local variables
   integer                              :: i_element
   integer                              :: i_element_br
      !================================================
      ! Transfer 32 bits of realIn to generic 32 bit integer space:
      i_element_br=0
      i_element = TRANSFER( IntIn, 0 )
      ! Reverse order of 4 bytes in 32 bit integer space:
      CALL MVBITS( i_element, 24, 8, i_element_br, 0  )
      CALL MVBITS( i_element, 16, 8, i_element_br, 8  )
      CALL MVBITS( i_element,  8, 8, i_element_br, 16 )
      CALL MVBITS( i_element,  0, 8, i_element_br, 24 )

      intIn = i_element_br
   END SUBROUTINE ENDIAN4I
   !####################################################################





   !####################################################################
   ! file I/O
   !-- INQUIRE_FID : inruire unused file FID
   !-- NCERROR     : netCDF I/O wrapper
   !####################################################################
   FUNCTION INQUIRE_FID() RESULT(FID)
   IMPLICIT NONE
   !* input/output
   integer :: FID ! FILE ID
   !* local variable
   logical :: I_OPENED ! FILE ID IS ALREADY USED or not?
   !================================================
      DO FID = 10, 999
         INQUIRE(FID,OPENED=I_OPENED)
         IF ( .not. I_OPENED ) RETURN
      ENDDO
   END FUNCTION INQUIRE_FID
   !==========================================================
   !+
   !+
   !+
   !==========================================================
#ifdef UseCDF_CMF
   SUBROUTINE NCERROR(STATUS,STRING)
   !! NETCDF error handling 
   USE NETCDF
   IMPLICIT NONE
   integer,intent(in)                     :: STATUS
   character(LEN=*),intent(in),OPTIONAL   :: STRING
   !================================================
   IF ( STATUS /= 0 ) THEN
      write(LOGNAM,*)  TRIM(NF90_STRERROR(STATUS))
      IF( PRESENT(STRING) ) write(LOGNAM,*) TRIM(STRING)
      write(LOGNAM,*) 'PROGRAM STOP ! '
      STOP 10
   ENDIF
   END SUBROUTINE NCERROR
#endif
!####################################################################

!####################################################################
FUNCTION CMF_CheckNanB(VAR,zero) RESULT(FLAG)
  implicit none
  REAL(KIND=JPRB)      :: VAR, zero
  LOGICAL              :: FLAG
  FLAG = .false.
  if(VAR*zero/=zero)then
    FLAG = .true.
  endif
END FUNCTION CMF_CheckNanB
!####################################################################

END MODULE CMF_UTILS_MOD
