MODULE MOD_IncompleteGamma

!     ALGORITHM 654, COLLECTED ALGORITHMS FROM ACM.
!     THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!     VOL. 13, NO. 3, P. 318.

CONTAINS

      SUBROUTINE GRATIO (A, X, ANS, QANS, IND)
!  ----------------------------------------------------------------------
!         EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
!                       P(A,X) AND Q(A,X)
!
!                         ----------
!
!      IT IS ASSUMED THAT A AND X ARE NONNEGATIVE, WHERE A AND X
!      ARE NOT BOTH 0.
!
!      ANS AND QANS ARE VARIABLES. GRATIO ASSIGNS ANS THE VALUE
!      P(A,X) AND QANS THE VALUE Q(A,X). IND MAY BE ANY INTEGER.
!      IF IND = 0 THEN THE USER IS REQUESTING AS MUCH ACCURACY AS
!      POSSIBLE (UP TO 14 SIGNIFICANT DIGITS). OTHERWISE, IF
!      IND = 1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT OF THE
!      6-TH SIGNIFICANT DIGIT, AND IF IND .NE. 0,1 THEN ACCURACY
!      IS REQUESTED TO WITHIN 1 UNIT OF THE 3RD SIGNIFICANT DIGIT.
!
!      ERROR RETURN ...
!         ANS IS ASSIGNED THE VALUE 2 WHEN A OR X IS NEGATIVE,
!      WHEN A*X = 0, OR WHEN P(A,X) AND Q(A,X) ARE INDETERMINANT.
!      P(A,X) AND Q(A,X) ARE COMPUTATIONALLY INDETERMINANT WHEN
!      X IS EXCEEDINGLY CLOSE TO A AND A IS EXTREMELY LARGE.
!  ----------------------------------------------------------------------
!      WRITTEN BY ALFRED H. MORRIS, JR.
!         NAVAL SURFACE WEAPONS CENTER
!         DAHLGREN, VIRGINIA
!      --------------------
      REAL J, L, ACC0(3), BIG(3), E00(3), X00(3), WK(20)
      REAL D0(13), D1(12), D2(10), D3(8), D4(6), D5(4), D6(2)
!     --------------------
      DATA ACC0(1)/5.E-15/, ACC0(2)/5.E-7/, ACC0(3)/5.E-4/
      DATA BIG(1)/20.0/, BIG(2)/14.0/, BIG(3)/10.0/
      DATA E00(1)/.25E-3/, E00(2)/.25E-1/, E00(3)/.14/
      DATA X00(1)/31.0/, X00(2)/17.0/, X00(3)/9.7/
!     --------------------
!     ALOG10 = LN(10)
!     RT2PIN = 1/SQRT(2*PI)
!     RTPI   = SQRT(PI)
!     --------------------
      DATA ALOG10/2.30258509299405/
      DATA RT2PIN/.398942280401433/
      DATA RTPI  /1.77245385090552/
      DATA THIRD /.333333333333333/
!     --------------------
      DATA D0(1) / .833333333333333E-01/, D0(2) /-.148148148148148E-01/, &
           D0(3) / .115740740740741E-02/, D0(4) / .352733686067019E-03/, &
           D0(5) /-.178755144032922E-03/, D0(6) / .391926317852244E-04/, &
           D0(7) /-.218544851067999E-05/, D0(8) /-.185406221071516E-05/, &
           D0(9) / .829671134095309E-06/, D0(10)/-.176659527368261E-06/, &
           D0(11)/ .670785354340150E-08/, D0(12)/ .102618097842403E-07/, &
           D0(13)/-.438203601845335E-08/
!     --------------------
      DATA D10   /-.185185185185185E-02/, D1(1) /-.347222222222222E-02/, &
           D1(2) / .264550264550265E-02/, D1(3) /-.990226337448560E-03/, &
           D1(4) / .205761316872428E-03/, D1(5) /-.401877572016461E-06/, &
           D1(6) /-.180985503344900E-04/, D1(7) / .764916091608111E-05/, &
           D1(8) /-.161209008945634E-05/, D1(9) / .464712780280743E-08/, &
           D1(10)/ .137863344691572E-06/, D1(11)/-.575254560351770E-07/, &
           D1(12)/ .119516285997781E-07/
!     --------------------
      DATA D20   / .413359788359788E-02/, D2(1) /-.268132716049383E-02/, &
           D2(2) / .771604938271605E-03/, D2(3) / .200938786008230E-05/, &
           D2(4) /-.107366532263652E-03/, D2(5) / .529234488291201E-04/, &
           D2(6) /-.127606351886187E-04/, D2(7) / .342357873409614E-07/, &
           D2(8) / .137219573090629E-05/, D2(9) /-.629899213838006E-06/, &
           D2(10)/ .142806142060642E-06/
!     --------------------
      DATA D30   / .649434156378601E-03/, D3(1) / .229472093621399E-03/, &
           D3(2) /-.469189494395256E-03/, D3(3) / .267720632062839E-03/, &
           D3(4) /-.756180167188398E-04/, D3(5) /-.239650511386730E-06/, &
           D3(6) / .110826541153473E-04/, D3(7) /-.567495282699160E-05/, &
           D3(8) / .142309007324359E-05/
!     --------------------
      DATA D40   /-.861888290916712E-03/, D4(1) / .784039221720067E-03/, &
           D4(2) /-.299072480303190E-03/, D4(3) /-.146384525788434E-05/, &
           D4(4) / .664149821546512E-04/, D4(5) /-.396836504717943E-04/, &
           D4(6) / .113757269706784E-04/
!     --------------------
      DATA D50   /-.336798553366358E-03/, D5(1) /-.697281375836586E-04/, &
           D5(2) / .277275324495939E-03/, D5(3) /-.199325705161888E-03/, &
           D5(4) / .679778047793721E-04/
!     --------------------
      DATA D60   / .531307936463992E-03/, D6(1) /-.592166437353694E-03/, &
           D6(2) / .270878209671804E-03/
!     --------------------
      DATA D70   / .344367606892378E-03/
!     --------------------
!     ****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
!            FLOATING POINT NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
!
!                    E = SPMPAR(1)
                    E = epsilon(1.)
!
!     --------------------
      IF (A .LT. 0.0 .OR. X .LT. 0.0) GO TO 400
      IF (A .EQ. 0.0 .AND. X .EQ. 0.0) GO TO 400
      IF (A*X .EQ. 0.0) GO TO 331
!
      IOP = IND + 1
      IF (IOP .NE. 1 .AND. IOP .NE. 2) IOP = 3
      ACC = AMAX1(ACC0(IOP),E)
      E0 = E00(IOP)
      X0 = X00(IOP)
!
!            SELECT THE APPROPRIATE ALGORITHM
!
      IF (A .GE. 1.0) GO TO 10
      IF (A .EQ. 0.5) GO TO 320
      IF (X .LT. 1.1) GO TO 110
      T1 = A*ALOG(X) - X
      U = A*EXP(T1)
      IF (U .EQ. 0.0) GO TO 310
      R = U*(1.0 + GAM1(A))
      GO TO 170
!
   10 IF (A .GE. BIG(IOP)) GO TO 20
      IF (A .GT. X .OR. X .GE. X0) GO TO 11
         TWOA = A + A
         M = INT(TWOA)
         IF (TWOA .NE. FLOAT(M)) GO TO 11
         I = M/2
         IF (A .EQ. FLOAT(I)) GO TO 140
         GO TO 150
   11 T1 = A*ALOG(X) - X
      R = EXP(T1)/GAMMA(A)
      GO TO 30
!
   20 L = X/A
      IF (L .EQ. 0.0) GO TO 300
      S = 0.5 + (0.5 - L)
      Z = RLOG(L)
      IF (Z .GE. 700.0/A) GO TO 330
      Y = A*Z
      RTA = SQRT(A)
      IF (ABS(S) .LE. E0/RTA) GO TO 250
      IF (ABS(S) .LE. 0.4) GO TO 200
!
      T = (1.0/A)**2
      T1 = (((0.75*T - 1.0)*T + 3.5)*T - 105.0)/(A*1260.0)
      T1 = T1 - Y
      R = RT2PIN*RTA*EXP(T1)
!
   30 IF (R .EQ. 0.0) GO TO 331
      IF (X .LE. AMAX1(A,ALOG10)) GO TO 50
      IF (X .LT. X0) GO TO 170
      GO TO 80
!
!                 TAYLOR SERIES FOR P/R
!
   50 APN = A + 1.0
      T = X/APN
      WK(1) = T
      DO N = 2,20
         APN = APN + 1.0
         T = T*(X/APN)
         IF (T .LE. 1.E-3) GO TO 60
         WK(N) = T
      ENDDO
      N = 20
!
   60 SUM = T
      TOL = 0.5*ACC
   61 APN = APN + 1.0
      T = T*(X/APN)
      SUM = SUM + T
      IF (T .GT. TOL) GO TO 61
!
      MAX = N - 1
      DO M = 1,MAX
         N = N - 1
         SUM = SUM + WK(N)
      ENDDO
      ANS = (R/A)*(1.0 + SUM)
      QANS = 0.5 + (0.5 - ANS)
      RETURN
!
!                 ASYMPTOTIC EXPANSION
!
   80 AMN = A - 1.0
      T = AMN/X
      WK(1) = T
      DO N = 2,20
         AMN = AMN - 1.0
         T = T*(AMN/X)
         IF (ABS(T) .LE. 1.E-3) GO TO 90
         WK(N) = T
      ENDDO
      N = 20
!
   90 SUM = T
   91 IF (ABS(T) .LE. ACC) GO TO 100
      AMN = AMN - 1.0
      T = T*(AMN/X)
      SUM = SUM + T
      GO TO 91
!
  100 MAX = N - 1
      DO M = 1,MAX
         N = N - 1
         SUM = SUM + WK(N)
      ENDDO
      QANS = (R/X)*(1.0 + SUM)
      ANS = 0.5 + (0.5 - QANS)
      RETURN
!
!             TAYLOR SERIES FOR P(A,X)/X**A
!
  110 AN = 3.0
      C = X
      SUM = X/(A + 3.0)
      TOL = 3.0*ACC/(A + 1.0)
  111    AN = AN + 1.0
         C = -C*(X/AN)
         T = C/(A + AN)
         SUM = SUM + T
         IF (ABS(T) .GT. TOL) GO TO 111
      J = A*X*((SUM/6.0 - 0.5/(A + 2.0))*X + 1.0/(A + 1.0))
!
      Z = A*ALOG(X)
      H = GAM1(A)
      G = 1.0 + H
      IF (X .LT. 0.25) GO TO 120
         IF (A .LT. X/2.59) GO TO 135
         GO TO 130
  120 IF (Z .GT. -.13394) GO TO 135
!
  130 W = EXP(Z)
      ANS = W*G*(0.5 + (0.5 - J))
      QANS = 0.5 + (0.5 - ANS)
      RETURN
!
  135 L = REXP(Z)
      W = 0.5 + (0.5 + L)
      QANS = (W*J - L)*G - H
      IF (QANS .LT. 0.0) GO TO 310
      ANS = 0.5 + (0.5 - QANS)
      RETURN
!
!             FINITE SUMS FOR Q WHEN A .GE. 1
!                 AND 2*A IS AN INTEGER
!
  140 SUM = EXP(-X)
      T = SUM
      N = 1
      C = 0.0
      GO TO 160
!
  150 RTX = SQRT(X)
      SUM = ERFC1(0,RTX)
      T = EXP(-X)/(RTPI*RTX)
      N = 0
      C = -0.5
!
  160 IF (N .EQ. I) GO TO 161
         N = N + 1
         C = C + 1.0
         T = (X*T)/C
         SUM = SUM + T
         GO TO 160
  161 QANS = SUM
      ANS = 0.5 + (0.5 - QANS)
      RETURN
!
!              CONTINUED FRACTION EXPANSION
!
  170 TOL = AMAX1(5.0*E,ACC)
      A2NM1 = 1.0
      A2N = 1.0
      B2NM1 = X
      B2N = X + (1.0 - A)
      C = 1.0
  171 A2NM1 = X*A2N + C*A2NM1
      B2NM1 = X*B2N + C*B2NM1
      AM0 = A2NM1/B2NM1
      C = C + 1.0
      CMA = C - A
      A2N = A2NM1 + CMA*A2N
      B2N = B2NM1 + CMA*B2N
      AN0 = A2N/B2N
      IF (ABS(AN0 - AM0) .GE. TOL*AN0) GO TO 171
!
      QANS = R*AN0
      ANS = 0.5 + (0.5 - QANS)
      RETURN
!
!                GENERAL TEMME EXPANSION
!
  200 IF (ABS(S) .LE. 2.0*E .AND. A*E*E .GT. 3.28E-3) GO TO 400
      C = EXP(-Y)
      W = 0.5*ERFC1(1,SQRT(Y))
      U = 1.0/A
      Z = SQRT(Z + Z)
      IF (L .LT. 1.0) Z = -Z
      IF (IOP < 2) THEN
         GO TO 210
      ELSEIF (IOP == 2) THEN
         GO TO 220
      ELSE
         GO TO 230
      ENDIF
!
  210 IF (ABS(S) .LE. 1.E-3) GO TO 260
      C0 = ((((((((((((D0(13) * Z + D0(12)) * Z + D0(11)) * Z    &
          + D0(10)) * Z + D0(9)) * Z + D0(8)) * Z + D0(7)) * Z   &
          + D0(6)) * Z + D0(5)) * Z + D0(4)) * Z + D0(3)) * Z    &
          + D0(2)) * Z + D0(1)) * Z - THIRD
      C1 = (((((((((((D1(12) * Z + D1(11)) * Z + D1(10)) * Z     &
          + D1(9)) * Z + D1(8)) * Z + D1(7)) * Z + D1(6)) * Z    &
          + D1(5)) * Z + D1(4)) * Z + D1(3)) * Z + D1(2)) * Z    &
          + D1(1)) * Z + D10
      C2 = (((((((((D2(10) * Z + D2(9)) * Z + D2(8)) * Z         &
          + D2(7)) * Z + D2(6)) * Z + D2(5)) * Z + D2(4)) * Z    &
          + D2(3)) * Z + D2(2)) * Z + D2(1)) * Z + D20
      C3 = (((((((D3(8) * Z + D3(7)) * Z + D3(6)) * Z            &
          + D3(5)) * Z + D3(4)) * Z + D3(3)) * Z + D3(2)) * Z    &
          + D3(1)) * Z + D30
      C4 = (((((D4(6) * Z + D4(5)) * Z + D4(4)) * Z + D4(3)) * Z &
          + D4(2)) * Z + D4(1)) * Z + D40
      C5 = (((D5(4) * Z + D5(3)) * Z + D5(2)) * Z + D5(1)) * Z   &
          + D50
      C6 = (D6(2) * Z + D6(1)) * Z + D60
      T  = ((((((D70*U + C6)*U + C5)*U + C4)*U + C3)*U + C2)*U   &
                      + C1)*U + C0
      GO TO 240
!
  220 C0 = (((((D0(6) * Z + D0(5)) * Z + D0(4)) * Z + D0(3)) * Z &
          + D0(2)) * Z + D0(1)) * Z - THIRD
      C1 = (((D1(4) * Z + D1(3)) * Z + D1(2)) * Z + D1(1)) * Z   &
          + D10
      C2 = D2(1) * Z + D20
      T  = (C2*U + C1)*U + C0
      GO TO 240
!
  230 T  = ((D0(3) * Z + D0(2)) * Z + D0(1)) * Z - THIRD
!
  240 IF (L .LT. 1.0) GO TO 241
      QANS = C*(W + RT2PIN*T/RTA)
      ANS = 0.5 + (0.5 - QANS)
      RETURN
  241 ANS = C*(W - RT2PIN*T/RTA)
      QANS = 0.5 + (0.5 - ANS)
      RETURN
!
!               TEMME EXPANSION FOR L = 1
!
  250 IF (A*E*E .GT. 3.28E-3) GO TO 400
      C = 0.5 + (0.5 - Y)
      W = (0.5 - SQRT(Y)*(0.5 + (0.5 - Y/3.0))/RTPI)/C
      U = 1.0/A
      Z = SQRT(Z + Z)
      IF (L .LT. 1.0) Z = -Z
      IF (IOP < 2) THEN
         GO TO 260
      ELSEIF (IOP == 2) THEN
         GO TO 270
      ELSE
         GO TO 280
      ENDIF
!
  260 C0 = ((((((D0(7) * Z + D0(6)) * Z + D0(5)) * Z + D0(4)) * Z  &
          + D0(3)) * Z + D0(2)) * Z + D0(1)) * Z - THIRD
      C1 = (((((D1(6) * Z + D1(5)) * Z + D1(4)) * Z + D1(3)) * Z   &
          + D1(2)) * Z + D1(1)) * Z + D10
      C2 = ((((D2(5) * Z + D2(4)) * Z + D2(3)) * Z + D2(2)) * Z    &
          + D2(1)) * Z + D20
      C3 = (((D3(4) * Z + D3(3)) * Z + D3(2)) * Z + D3(1)) * Z     &
          + D30
      C4 = (D4(2) * Z + D4(1)) * Z + D40
      C5 = (D5(2) * Z + D5(1)) * Z + D50
      C6 = D6(1) * Z + D60
      T  = ((((((D70*U + C6)*U + C5)*U + C4)*U + C3)*U + C2)*U     &
                      + C1)*U + C0
      GO TO 240
!
  270 C0 = (D0(2) * Z + D0(1)) * Z - THIRD
      C1 = D1(1) * Z + D10
      T  = (D20*U + C1)*U + C0
      GO TO 240
!
  280 T  = D0(1) * Z - THIRD
      GO TO 240
!
!                     SPECIAL CASES
!
  300 ANS = 0.0
      QANS = 1.0
      RETURN
!
  310 ANS = 1.0
      QANS = 0.0
      RETURN
!
  320 IF (X .GE. 0.25) GO TO 321
      ANS = ERF(SQRT(X))
      QANS = 0.5 + (0.5 - ANS)
      RETURN
  321 QANS = ERFC1(0,SQRT(X))
      ANS = 0.5 + (0.5 - QANS)
      RETURN
!
  330 IF (ABS(S) .LE. 2.0*E) GO TO 400
  331 IF (X .LE. A) GO TO 300
      GO TO 310
!
!                     ERROR RETURN
!
  400 ANS = 2.0
      RETURN

      END SUBROUTINE



      FUNCTION ERF(X)
!     ******************************************************************
!     EVALUATION OF THE REAL ERROR FUNCTION
!     ******************************************************************
      DIMENSION A(4),B(4),P(8),Q(8),R(5),S(5)
      DATA A(1)/-1.65581836870402E-4/, A(2)/3.25324098357738E-2/, &
           A(3)/1.02201136918406E-1/,  A(4)/1.12837916709552E00/
      DATA B(1)/4.64988945913179E-3/,  B(2)/7.01333417158511E-2/, &
           B(3)/4.23906732683201E-1/,  B(4)/1.00000000000000E00/
      DATA P(1)/-1.36864857382717E-7/, P(2)/5.64195517478974E-1/, &
           P(3)/7.21175825088309E00/,  P(4)/4.31622272220567E01/, &
           P(5)/1.52989285046940E02/,  P(6)/3.39320816734344E02/, &
           P(7)/4.51918953711873E02/,  P(8)/3.00459261020162E02/
      DATA Q(1)/1.00000000000000E00/,  Q(2)/1.27827273196294E01/, &
           Q(3)/7.70001529352295E01/,  Q(4)/2.77585444743988E02/, &
           Q(5)/6.38980264465631E02/,  Q(6)/9.31354094850610E02/, &
           Q(7)/7.90950925327898E02/,  Q(8)/3.00459260956983E02/
      DATA R(1)/2.10144126479064E00/,  R(2)/2.62370141675169E01/, &
           R(3)/2.13688200555087E01/,  R(4)/4.65807828718470E00/, &
           R(5)/2.82094791773523E-1/
      DATA S(1)/9.41537750555460E01/,  S(2)/1.87114811799590E02/, &
           S(3)/9.90191814623914E01/,  S(4)/1.80124575948747E01/, &
           S(5)/1.00000000000000E00/
      DATA C/5.64189583547756E-1/
!     -------------------
      AX=ABS(X)
      X2=AX*AX
      IF (AX.GE.0.5) GO TO 10
      TOP=((A(1)*X2+A(2))*X2+A(3))*X2+A(4)
      BOT=((B(1)*X2+B(2))*X2+B(3))*X2+B(4)
      ERF=X*TOP/BOT
      RETURN
!
   10 IF (AX.GT.4.0) GO TO 20
      TOP=((((((P(1)*AX+P(2))*AX+P(3))*AX+P(4))*AX+P(5))*AX &
                      +P(6))*AX+P(7))*AX+P(8)
      BOT=((((((Q(1)*AX+Q(2))*AX+Q(3))*AX+Q(4))*AX+Q(5))*AX &
                      +Q(6))*AX+Q(7))*AX+Q(8)
      ERF=1.0-EXP(-X2)*TOP/BOT
      IF (X.LT.0.0) ERF=-ERF
      RETURN
!
   20 ERF=1.0
      IF (AX.GE.5.54) GO TO 21
      T=1.0/X2
      TOP=(((R(1)*T+R(2))*T+R(3))*T+R(4))*T+R(5)
      BOT=(((S(1)*T+S(2))*T+S(3))*T+S(4))*T+S(5)
      ERF=C-TOP/(X2*BOT)
      ERF=1.0-EXP(-X2)*ERF/AX
   21 IF (X.LT.0.0) ERF=-ERF
      RETURN

      END FUNCTION



      REAL FUNCTION ERFC1(IND,X)
! ----------------------------------------------------------------------
!     EVALUATION OF THE REAL COMPLEMENTARY ERROR FUNCTION
!
!        ERFC1(IND,X) = ERFC(X)            IF IND = 0
!        ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
! ----------------------------------------------------------------------
      DIMENSION A(4),B(4),P(8),Q(8),R(5),S(5)
      DATA A(1)/-1.65581836870402E-4/, A(2)/3.25324098357738E-2/, &
           A(3)/1.02201136918406E-1/,  A(4)/1.12837916709552E00/
      DATA B(1)/4.64988945913179E-3/,  B(2)/7.01333417158511E-2/, &
           B(3)/4.23906732683201E-1/,  B(4)/1.00000000000000E00/
      DATA P(1)/-1.36864857382717E-7/, P(2)/5.64195517478974E-1/, &
           P(3)/7.21175825088309E00/,  P(4)/4.31622272220567E01/, &
           P(5)/1.52989285046940E02/,  P(6)/3.39320816734344E02/, &
           P(7)/4.51918953711873E02/,  P(8)/3.00459261020162E02/
      DATA Q(1)/1.00000000000000E00/,  Q(2)/1.27827273196294E01/, &
           Q(3)/7.70001529352295E01/,  Q(4)/2.77585444743988E02/, &
           Q(5)/6.38980264465631E02/,  Q(6)/9.31354094850610E02/, &
           Q(7)/7.90950925327898E02/,  Q(8)/3.00459260956983E02/
      DATA R(1)/2.10144126479064E00/,  R(2)/2.62370141675169E01/, &
           R(3)/2.13688200555087E01/,  R(4)/4.65807828718470E00/, &
           R(5)/2.82094791773523E-1/
      DATA S(1)/9.41537750555460E01/,  S(2)/1.87114811799590E02/, &
           S(3)/9.90191814623914E01/,  S(4)/1.80124575948747E01/, &
           S(5)/1.00000000000000E00/
      DATA C/5.64189583547756E-1/
!     -------------------
      AX=ABS(X)
      X2=AX*AX
      IF (AX.GE.0.47) GO TO 10
      TOP=((A(1)*X2+A(2))*X2+A(3))*X2+A(4)
      BOT=((B(1)*X2+B(2))*X2+B(3))*X2+B(4)
      ERFC1=1.0-X*TOP/BOT
      IF (IND.NE.0) ERFC1=EXP(X2)*ERFC1
      RETURN
!
   10 IF (AX.GT.4.0) GO TO 20
      TOP=((((((P(1)*AX+P(2))*AX+P(3))*AX+P(4))*AX+P(5))*AX &
                      +P(6))*AX+P(7))*AX+P(8)
      BOT=((((((Q(1)*AX+Q(2))*AX+Q(3))*AX+Q(4))*AX+Q(5))*AX &
                      +Q(6))*AX+Q(7))*AX+Q(8)
      ERFC1=TOP/BOT
      IF (IND.EQ.0) GO TO 11
      IF (X.LT.0.0) ERFC1=2.0*EXP(X2)-ERFC1
      RETURN
   11 ERFC1=EXP(-X2)*ERFC1
      IF (X.LT.0.0) ERFC1=2.0-ERFC1
      RETURN
!
   20 IF (X.LE.-5.33) GO TO 30
      T=1.0/X2
      TOP=(((R(1)*T+R(2))*T+R(3))*T+R(4))*T+R(5)
      BOT=(((S(1)*T+S(2))*T+S(3))*T+S(4))*T+S(5)
      ERFC1=(C-TOP/(X2*BOT))/AX
      IF (IND.EQ.0) GO TO 11
      IF (X.LT.0.0) ERFC1=2.0*EXP(X2)-ERFC1
      RETURN
!
   30 ERFC1=2.0
      IF (IND.NE.0) ERFC1=EXP(X2)*ERFC1
      RETURN

      END FUNCTION



      REAL FUNCTION REXP(X)
!     ------------------------------------------------------------------
!     COMPUTATION OF EXP(X) - 1
!     ------------------------------------------------------------------
      DATA P1/ .914041914819518E-09/, P2/ .238082361044469E-01/, &
           Q1/-.499999999085958E+00/, Q2/ .107141568980644E+00/, &
           Q3/-.119041179760821E-01/, Q4/ .595130811860248E-03/
!     ------------------
      IF (ABS(X) .GT. 0.15) GO TO 10
      REXP = X*(((P2*X + P1)*X + 1.0)/((((Q4*X + Q3)*X + Q2)*X &
                      + Q1)*X + 1.0))
      RETURN
!
   10 W = EXP(X)
      IF (X .GT. 0.0) GO TO 20
         REXP = (W - 0.5) - 0.5
         RETURN
   20 REXP = W*(0.5 + (0.5 - 1.0/W))
      RETURN
      END FUNCTION



      REAL FUNCTION RLOG(X)
!     -------------------
!     COMPUTATION OF  X - 1 - LN(X)
!     -------------------
      DATA A/.566749439387324E-01/
      DATA B/.456512608815524E-01/
!     -------------------
      DATA P0/ .333333333333333E+00/, P1/-.224696413112536E+00/, &
           P2/ .620886815375787E-02/
      DATA Q1/-.127408923933623E+01/, Q2/ .354508718369557E+00/
!     -------------------
      IF (X .LT. 0.61 .OR. X .GT. 1.57) GO TO 100
      IF (X .LT. 0.82) GO TO 10
      IF (X .GT. 1.18) GO TO 20
!
!              ARGUMENT REDUCTION
!
      U = (X - 0.5) - 0.5
      W1 = 0.0
      GO TO 30
!
   10 U = DBLE(X) - 0.7D0
      U = U/0.7
      W1 = A - U*0.3
      GO TO 30
!
   20 U = 0.75D0*DBLE(X) - 1.D0
      W1 = B + U/3.0
!
!               SERIES EXPANSION
!
   30 R = U/(U + 2.0)
      T = R*R
      W = ((P2*T + P1)*T + P0)/((Q2*T + Q1)*T + 1.0)
      RLOG = 2.0*T*(1.0/(1.0 - R) - R*W) + W1
      RETURN
!
!
  100 R = (X - 0.5) - 0.5
      RLOG = R - ALOG(X)
      RETURN

      END FUNCTION



      REAL FUNCTION GAMMA(A)
!-----------------------------------------------------------------------
!
!         EVALUATION OF THE GAMMA FUNCTION FOR REAL ARGUMENTS
!
!                           -----------
!
!     GAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT
!     BE COMPUTED.
!
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!          NAVAL SURFACE WEAPONS CENTER
!          DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
      REAL P(7), Q(7)
      DOUBLE PRECISION D, G, Z, LNX
!--------------------------
!     D = 0.5*(LN(2*PI) - 1)
!--------------------------
      DATA PI /3.1415926535898/
      DATA D /.41893853320467274178D0/
!--------------------------
      DATA P(1)/ .539637273585445E-03/,  P(2)/ .261939260042690E-02/, &
           P(3)/ .204493667594920E-01/,  P(4)/ .730981088720487E-01/, &
           P(5)/ .279648642639792E+00/,  P(6)/ .553413866010467E+00/, &
           P(7)/ 1.0/
      DATA Q(1)/-.832979206704073E-03/,  Q(2)/ .470059485860584E-02/, &
           Q(3)/ .225211131035340E-01/,  Q(4)/-.170458969313360E+00/, &
           Q(5)/-.567902761974940E-01/,  Q(6)/ .113062953091122E+01/, &
           Q(7)/ 1.0/
!--------------------------
      DATA R1/.820756370353826E-03/, R2/-.595156336428591E-03/, &
           R3/.793650663183693E-03/, R4/-.277777777770481E-02/, &
           R5/.833333333333333E-01/
!--------------------------
      GAMMA = 0.0
      X = A
      IF (ABS(A) .GE. 15.0) GO TO 60
!-----------------------------------------------------------------------
!            EVALUATION OF GAMMA(A) FOR ABS(A) .LT. 15
!-----------------------------------------------------------------------
      T = 1.0
      M = INT(A) - 1
!
!     LET T BE THE PRODUCT OF A-J WHEN A .GE. 2
!
      IF (M < 0) THEN
         GO TO 20
      ELSEIF (M == 0) THEN
         GO TO 12
      ELSE
         GO TO 10
      ENDIF

   10 DO J = 1,M
        X = X - 1.0
        T = X*T
     ENDDO
   12 X = X - 1.0
      GO TO 40
!
!     LET T BE THE PRODUCT OF A+J WHEN A .LT. 1
!
   20 T = A
      IF (A .GT. 0.0) GO TO 30
      M = - M - 1
      IF (M .EQ. 0) GO TO 22
         DO J = 1,M
            X = X + 1.0
            T = X*T
         ENDDO
   22 X = (X + 0.5) + 0.5
      T = X*T
      IF (T .EQ. 0.0) RETURN
!
   30 CONTINUE
!
!     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
!     CODE MAY BE OMITTED IF DESIRED.
!
      IF (ABS(T) .GE. 1.E-30) GO TO 40
!      IF (ABS(T)*SPMPAR(3) .LE. 1.0001) RETURN
      IF (ABS(T)*HUGE(1.) .LE. 1.0001) RETURN
      GAMMA = 1.0/T
      RETURN
!
!     COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1
!
   40 TOP = P(1)
      BOT = Q(1)
      DO I = 2,7
         TOP = P(I) + X*TOP
         BOT = Q(I) + X*BOT
      ENDDO
      GAMMA = TOP/BOT
!
!     TERMINATION
!
      IF (A .LT. 1.0) GO TO 50
      GAMMA = GAMMA*T
      RETURN
   50 GAMMA = GAMMA/T
      RETURN
!-----------------------------------------------------------------------
!            EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15
!-----------------------------------------------------------------------
   60 IF (ABS(A) .GE. 1.E3) RETURN
      IF (A .GT. 0.0) GO TO 70
      X = -A
      N = X
      T = X - N
      IF (T .GT. 0.9) T = 1.0 - T
      S = SIN(PI*T)/PI
      IF (MOD(N,2) .EQ. 0) S = -S
      IF (S .EQ. 0.0) RETURN
!
!     COMPUTE THE MODIFIED ASYMPTOTIC SUM
!
   70 T = 1.0/(X*X)
      G = ((((R1*T + R2)*T + R3)*T + R4)*T + R5)/X
!
!     ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X)
!     BUT LESS ACCURACY WILL NORMALLY BE OBTAINED.
!
      LNX = GLOG(X)
!
!     FINAL ASSEMBLY
!
      Z = X
      G = (D + G) + (Z - 0.5D0)*(LNX - 1.D0)
      W = G
      T = G - DBLE(W)
      IF (W .GT. 0.99999*EXPARG(0)) RETURN
      GAMMA = EXP(W)*(1.0 + T)
      IF (A .LT. 0.0) GAMMA = (1.0/(GAMMA*S))/X
      RETURN

      END FUNCTION



      DOUBLE PRECISION FUNCTION GLOG(X)
!     -------------------
!     EVALUATION OF LN(X) FOR X .GE. 15
!     -------------------
      REAL X
      DOUBLE PRECISION Z, W(163)
!     -------------------
      DATA C1/.286228750476730/, C2/.399999628131494/, &
           C3/.666666666752663/
!     -------------------
!     W(J) = LN(J + 14) FOR EACH J
!     -------------------
      DATA W(1) /.270805020110221007D+01/,                                   &
           W(2) /.277258872223978124D+01/, W(3) /.283321334405621608D+01/,   &
           W(4) /.289037175789616469D+01/, W(5) /.294443897916644046D+01/,   &
           W(6) /.299573227355399099D+01/, W(7) /.304452243772342300D+01/,   &
           W(8) /.309104245335831585D+01/, W(9) /.313549421592914969D+01/,   &
           W(10)/.317805383034794562D+01/, W(11)/.321887582486820075D+01/,   &
           W(12)/.325809653802148205D+01/, W(13)/.329583686600432907D+01/,   &
           W(14)/.333220451017520392D+01/, W(15)/.336729582998647403D+01/,   &
           W(16)/.340119738166215538D+01/, W(17)/.343398720448514625D+01/,   &
           W(18)/.346573590279972655D+01/, W(19)/.349650756146648024D+01/,   &
           W(20)/.352636052461616139D+01/, W(21)/.355534806148941368D+01/,   &
           W(22)/.358351893845611000D+01/, W(23)/.361091791264422444D+01/,   &
           W(24)/.363758615972638577D+01/, W(25)/.366356164612964643D+01/,   &
           W(26)/.368887945411393630D+01/, W(27)/.371357206670430780D+01/,   &
           W(28)/.373766961828336831D+01/, W(29)/.376120011569356242D+01/,   &
           W(30)/.378418963391826116D+01/
      DATA W(31)/.380666248977031976D+01/,                                   &
           W(32)/.382864139648909500D+01/, W(33)/.385014760171005859D+01/,   &
           W(34)/.387120101090789093D+01/, W(35)/.389182029811062661D+01/,   &
           W(36)/.391202300542814606D+01/, W(37)/.393182563272432577D+01/,   &
           W(38)/.395124371858142735D+01/, W(39)/.397029191355212183D+01/,   &
           W(40)/.398898404656427438D+01/, W(41)/.400733318523247092D+01/,   &
           W(42)/.402535169073514923D+01/, W(43)/.404305126783455015D+01/,   &
           W(44)/.406044301054641934D+01/, W(45)/.407753744390571945D+01/,   &
           W(46)/.409434456222210068D+01/, W(47)/.411087386417331125D+01/,   &
           W(48)/.412713438504509156D+01/, W(49)/.414313472639153269D+01/,   &
           W(50)/.415888308335967186D+01/, W(51)/.417438726989563711D+01/,   &
           W(52)/.418965474202642554D+01/, W(53)/.420469261939096606D+01/,   &
           W(54)/.421950770517610670D+01/, W(55)/.423410650459725938D+01/,   &
           W(56)/.424849524204935899D+01/, W(57)/.426267987704131542D+01/,   &
           W(58)/.427666611901605531D+01/, W(59)/.429045944114839113D+01/,   &
           W(60)/.430406509320416975D+01/
     DATA  W(61)/.431748811353631044D+01/,                                   &
           W(62)/.433073334028633108D+01/, W(63)/.434380542185368385D+01/,   &
           W(64)/.435670882668959174D+01/, W(65)/.436944785246702149D+01/,   &
           W(66)/.438202663467388161D+01/, W(67)/.439444915467243877D+01/,   &
           W(68)/.440671924726425311D+01/, W(69)/.441884060779659792D+01/,   &
           W(70)/.443081679884331362D+01/, W(71)/.444265125649031645D+01/,   &
           W(72)/.445434729625350773D+01/, W(73)/.446590811865458372D+01/,   &
           W(74)/.447733681447820647D+01/, W(75)/.448863636973213984D+01/,   &
           W(76)/.449980967033026507D+01/, W(77)/.451085950651685004D+01/,   &
           W(78)/.452178857704904031D+01/, W(79)/.453259949315325594D+01/,   &
           W(80)/.454329478227000390D+01/, W(81)/.455387689160054083D+01/,   &
           W(82)/.456434819146783624D+01/, W(83)/.457471097850338282D+01/,   &
           W(84)/.458496747867057192D+01/, W(85)/.459511985013458993D+01/,   &
           W(86)/.460517018598809137D+01/, W(87)/.461512051684125945D+01/,   &
           W(88)/.462497281328427108D+01/, W(89)/.463472898822963577D+01/,   &
           W(90)/.464439089914137266D+01/
     DATA  W(91) /.465396035015752337D+01/,                                  &
           W(92) /.466343909411206714D+01/, W(93) /.467282883446190617D+01/, &
           W(94) /.468213122712421969D+01/, W(95) /.469134788222914370D+01/, &
           W(96) /.470048036579241623D+01/, W(97) /.470953020131233414D+01/, &
           W(98) /.471849887129509454D+01/, W(99) /.472738781871234057D+01/, &
           W(100)/.473619844839449546D+01/, W(101)/.474493212836325007D+01/, &
           W(102)/.475359019110636465D+01/, W(103)/.476217393479775612D+01/, &
           W(104)/.477068462446566476D+01/, W(105)/.477912349311152939D+01/, &
           W(106)/.478749174278204599D+01/, W(107)/.479579054559674109D+01/, &
           W(108)/.480402104473325656D+01/, W(109)/.481218435537241750D+01/, &
           W(110)/.482028156560503686D+01/, W(111)/.482831373730230112D+01/, &
           W(112)/.483628190695147800D+01/, W(113)/.484418708645859127D+01/, &
           W(114)/.485203026391961717D+01/, W(115)/.485981240436167211D+01/, &
           W(116)/.486753445045558242D+01/, W(117)/.487519732320115154D+01/, &
           W(118)/.488280192258637085D+01/, W(119)/.489034912822175377D+01/, &
           W(120)/.489783979995091137D+01/
     DATA  W(121)/.490527477843842945D+01/,                                  &
           W(122)/.491265488573605201D+01/, W(123)/.491998092582812492D+01/, &
           W(124)/.492725368515720469D+01/, W(125)/.493447393313069176D+01/, &
           W(126)/.494164242260930430D+01/, W(127)/.494875989037816828D+01/, &
           W(128)/.495582705760126073D+01/, W(129)/.496284463025990728D+01/, &
           W(130)/.496981329957600062D+01/, W(131)/.497673374242057440D+01/, &
           W(132)/.498360662170833644D+01/, W(133)/.499043258677873630D+01/, &
           W(134)/.499721227376411506D+01/, W(135)/.500394630594545914D+01/, &
           W(136)/.501063529409625575D+01/, W(137)/.501727983681492433D+01/, &
           W(138)/.502388052084627639D+01/, W(139)/.503043792139243546D+01/, &
           W(140)/.503695260241362916D+01/, W(141)/.504342511691924662D+01/, &
           W(142)/.504985600724953705D+01/, W(143)/.505624580534830806D+01/, &
           W(144)/.506259503302696680D+01/, W(145)/.506890420222023153D+01/, &
           W(146)/.507517381523382692D+01/, W(147)/.508140436498446300D+01/, &
           W(148)/.508759633523238407D+01/, W(149)/.509375020080676233D+01/, &
           W(150)/.509986642782419842D+01/
     DATA  W(151)/.510594547390058061D+01/,                                  &
           W(152)/.511198778835654323D+01/, W(153)/.511799381241675511D+01/, &
           W(154)/.512396397940325892D+01/, W(155)/.512989871492307347D+01/, &
           W(156)/.513579843705026176D+01/, W(157)/.514166355650265984D+01/, &
           W(158)/.514749447681345304D+01/, W(159)/.515329159449777895D+01/, &
           W(160)/.515905529921452903D+01/, W(161)/.516478597392351405D+01/, &
           W(162)/.517048399503815178D+01/, W(163)/.517614973257382914D+01/
!
      IF (X .GE. 178.0) GO TO 10
      N = X
      T = (X - N)/(X + N)
      T2 = T*T
      Z = (((C1*T2 + C2)*T2 + C3)*T2 + 2.0)*T
      GLOG = W(N - 14) + Z
      RETURN
!
   10 GLOG = ALOG(X)
      RETURN

      END FUNCTION



      REAL FUNCTION EXPARG (IDUMMY)
!--------------------------------------------------------------------
!     COMPUTATION OF THE LARGEST ARGUMENT W FOR WHICH EXP(W)
!     MAY BE COMPUTED. (ONLY AN APPROXIMATE VALUE IS NEEDED.)
!--------------------------------------------------------------------
!      EXPARG = 0.99999*ALOG(SPMPAR(3))
      EXPARG = 0.99999*ALOG(HUGE(1.))
      RETURN
      END



      REAL FUNCTION GAM1(A)
!     ------------------------------------------------------------------
!     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5
!     ------------------------------------------------------------------
      REAL P(7), Q(5), R(9)
!     -------------------
      DATA P(1)/ .577215664901533E+00/, P(2)/-.409078193005776E+00/, &
           P(3)/-.230975380857675E+00/, P(4)/ .597275330452234E-01/, &
           P(5)/ .766968181649490E-02/, P(6)/-.514889771323592E-02/, &
           P(7)/ .589597428611429E-03/
!     -------------------
      DATA Q(1)/ .100000000000000E+01/, Q(2)/ .427569613095214E+00/, &
           Q(3)/ .158451672430138E+00/, Q(4)/ .261132021441447E-01/, &
           Q(5)/ .423244297896961E-02/
!     -------------------
      DATA R(1)/-.422784335098468E+00/, R(2)/-.771330383816272E+00/, &
           R(3)/-.244757765222226E+00/, R(4)/ .118378989872749E+00/, &
           R(5)/ .930357293360349E-03/, R(6)/-.118290993445146E-01/, &
           R(7)/ .223047661158249E-02/, R(8)/ .266505979058923E-03/, &
           R(9)/-.132674909766242E-03/
!     -------------------
      DATA S1  / .273076135303957E+00/, S2  / .559398236957378E-01/
!     -------------------
      T = A
      D = A - 0.5
      IF (D .GT. 0.0) T = D - 0.5
      IF (T < 0) THEN
         GO TO 30
      ELSEIF (T == 0) THEN
         GO TO 10
      ELSE
         GO TO 20
      ENDIF
!
   10 GAM1 = 0.0
      RETURN
!
   20 TOP = (((((P(7)*T + P(6))*T + P(5))*T + P(4))*T + P(3))*T &
                        + P(2))*T + P(1)
      BOT = (((Q(5)*T + Q(4))*T + Q(3))*T + Q(2))*T + 1.0
      W = TOP/BOT
      IF (D .GT. 0.0) GO TO 21
         GAM1 = A*W
         RETURN
   21 GAM1 = (T/A)*((W - 0.5) - 0.5)
      RETURN
!
   30 TOP = (((((((R(9)*T + R(8))*T + R(7))*T + R(6))*T + R(5))*T &
                          + R(4))*T + R(3))*T + R(2))*T + R(1)
      BOT = (S2*T + S1)*T + 1.0
      W = TOP/BOT
      IF (D .GT. 0.0) GO TO 31
         GAM1 = A*((W + 0.5) + 0.5)
         RETURN
   31 GAM1 = T*W/A
      RETURN

      END FUNCTION

END MODULE MOD_IncompleteGamma
