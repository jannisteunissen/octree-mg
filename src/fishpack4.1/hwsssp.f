C
C     file hwsssp.f
C
      SUBROUTINE HWSSSP (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,
     1                   BDPF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                  copyright (c) 1999 by UCAR                   *
C     *                                                               *
C     *       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH         *
C     *                                                               *
C     *                      all rights reserved                      *
C     *                                                               *
C     *                      FISHPACK version 4.1                     *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C
C DIMENSION OF           BDTS(N+1),    BDTF(N+1), BDPS(M+1), BDPF(M+1),
C ARGUMENTS              F(IDIMF,N+1), W(SEE ARGUMENT LIST)
C
C LATEST REVISION        NOVEMBER 1988
C
C PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION TO
C                        THE HELMHOLTZ EQUATION IN SPHERICAL
C                        COORDINATES AND ON THE SURFACE OF THE UNIT
C                        SPHERE (RADIUS OF 1).  THE EQUATION IS
C
C                          (1/SIN(THETA))(D/DTHETA)(SIN(THETA)
C                          (DU/DTHETA)) + (1/SIN(THETA)**2)(D/DPHI)
C                          (DU/DPHI)  + LAMBDA*U = F(THETA,PHI)
C
C                        WHERE THETA IS COLATITUDE AND PHI IS
C                        LONGITUDE.
C
C USAGE                  CALL HWSSSP (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,
C                                     N,NBDCND,BDPS,BDPF,ELMBDA,F,
C                                     IDIMF,PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT               TS,TF
C
C                          THE RANGE OF THETA (COLATITUDE), I.E.,
C                          TS .LE. THETA .LE. TF. TS MUST BE LESS
C                          THAN TF.  TS AND TF ARE IN RADIANS.
C                          A TS OF ZERO CORRESPONDS TO THE NORTH
C                          POLE AND A TF OF PI CORRESPONDS TO
C                          THE SOUTH POLE.
C
C                          * * * IMPORTANT * * *
C
C                          IF TF IS EQUAL TO PI THEN IT MUST BE
C                          COMPUTED USING THE STATEMENT
C                          TF = PIMACH(DUM). THIS INSURES THAT TF
C                          IN THE USER'S PROGRAM IS EQUAL TO PI IN
C                          THIS PROGRAM WHICH PERMITS SEVERAL TESTS
C                          OF THE INPUT PARAMETERS THAT OTHERWISE
C                          WOULD NOT BE POSSIBLE.
C
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (TS,TF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE M+1 GRID POINTS IN THE
C                          THETA-DIRECTION GIVEN BY
C                          THETA(I) = (I-1)DTHETA+TS FOR
C                          I = 1,2,...,M+1, WHERE
C                          DTHETA = (TF-TS)/M IS THE PANEL WIDTH.
C                          M MUST BE GREATER THAN 5
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT THETA = TS AND THETA = TF.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THETA = TF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE DERIVATIVE OF
C                               THE SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               SPECIFIED AT THETA = TS AND
C                               THETA = TF (SEE NOTES 1,2 BELOW).
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW)
C                               AND THE SOLUTION IS SPECIFIED AT
C                               THETA = TF.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE SOLUTION
C                               IS SPECIFIED AT THETA = TF.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO THETA
C                               IS SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 7  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE SOLUTION IS
C                               IS UNSPECIFIED AT THETA = TF = PI.
C                          = 8  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW) AND
C                               THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TF = PI.
C                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THETA = TF = PI.
C
C                          NOTES:
C                          IF TS = 0, DO NOT USE MBDCND = 3,4, OR 8,
C                          BUT INSTEAD USE MBDCND = 5,6, OR 9  .
C
C                          IF TF = PI, DO NOT USE MBDCND = 2,3, OR 6,
C                          BUT INSTEAD USE MBDCND = 7,8, OR 9  .
C
C                        BDTS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TS.  WHEN MBDCND = 3,4, OR 8,
C
C                          BDTS(J) = (D/DTHETA)U(TS,PHI(J)),
C                          J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTS IS
C                          A DUMMY VARIABLE.
C
C                        BDTF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
C                          THAT SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TF.  WHEN MBDCND = 2,3, OR 6,
C
C                          BDTF(J) = (D/DTHETA)U(TF,PHI(J)),
C                          J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTF IS
C                          A DUMMY VARIABLE.
C
C                        PS,PF
C                          THE RANGE OF PHI (LONGITUDE), I.E.,
C                          PS .LE. PHI .LE. PF.  PS MUST BE LESS
C                          THAN PF.  PS AND PF ARE IN RADIANS.
C                          IF PS = 0 AND PF = 2*PI, PERIODIC
C                          BOUNDARY CONDITIONS ARE USUALLY PRESCRIBED.
C
C                          * * * IMPORTANT * * *
C
C                          IF PF IS EQUAL TO 2*PI THEN IT MUST BE
C                          COMPUTED USING THE STATEMENT
C                          PF = 2.*PIMACH(DUM). THIS INSURES THAT
C                          PF IN THE USERS PROGRAM IS EQUAL TO
C                          2*PI IN THIS PROGRAM WHICH PERMITS TESTS
C                          OF THE INPUT PARAMETERS THAT OTHERWISE
C                          WOULD NOT BE POSSIBLE.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (PS,PF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS
C                          IN THE PHI-DIRECTION GIVEN BY
C                          PHI(J) = (J-1)DPHI+PS  FOR
C                          J = 1,2,...,N+1, WHERE
C                          DPHI = (PF-PS)/N IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 4
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT PHI = PS AND PHI = PF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN PHI,
C                               I.U., U(I,J) = U(I,N+J).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               PHI = PS AND PHI = PF
C                               (SEE NOTE BELOW).
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               PHI = PS (SEE NOTE BELOW)
C                               AND THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PHI = PF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PHI = PS AND PHI = PF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO PHI IS SPECIFIED
C                               AT PS AND THE SOLUTION IS SPECIFIED
C                               AT PHI = PF
C
C                          NOTE:
C                          NBDCND = 1,2, OR 4 CANNOT BE USED WITH
C                          MBDCND = 5,6,7,8, OR 9.  THE FORMER INDICATES
C                          THAT THE SOLUTION IS SPECIFIED AT A POLE, THE
C                          LATTER INDICATES THAT THE SOLUTION IS NOT
C                          SPECIFIED.  USE INSTEAD  MBDCND = 1 OR 2.
C
C                        BDPS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO PHI AT
C                          PHI = PS.  WHEN NBDCND = 3 OR 4,
C
C                            BDPS(I) = (D/DPHI)U(THETA(I),PS),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDPS IS
C                          A DUMMY VARIABLE.
C
C                        BDPF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE
C                          OF THE SOLUTION WITH RESPECT TO PHI AT
C                          PHI = PF.  WHEN NBDCND = 2 OR 3,
C
C                            BDPF(I) = (D/DPHI)U(THETA(I),PF),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDPF IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSSSP WILL
C                          ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUE OF THE RIGHT SIDE OF THE HELMHOLTZ
C                          EQUATION AND BOUNDARY VALUES (IF ANY).
C                          F MUST BE DIMENSIONED AT LEAST (M+1)*(N+1).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(THETA(I),PHI(J)).
C
C                          ON THE BOUNDARIES F IS DEFINED AS FOLLOWS:
C                          FOR J = 1,2,...,N+1 AND I = 1,2,...,M+1
C
C                          MBDCND   F(1,J)            F(M+1,J)
C                          ------   ------------      ------------
C
C                            1      U(TS,PHI(J))      U(TF,PHI(J))
C                            2      U(TS,PHI(J))      F(TF,PHI(J))
C                            3      F(TS,PHI(J))      F(TF,PHI(J))
C                            4      F(TS,PHI(J))      U(TF,PHI(J))
C                            5      F(0,PS)           U(TF,PHI(J))
C                            6      F(0,PS)           F(TF,PHI(J))
C                            7      U(TS,PHI(J))      F(PI,PS)
C                            8      F(TS,PHI(J))      F(PI,PS)
C                            9      F(0,PS)           F(PI,PS)
C
C                          NBDCND   F(I,1)            F(I,N+1)
C                          ------   --------------    --------------
C
C                            0      F(THETA(I),PS)    F(THETA(I),PS)
C                            1      U(THETA(I),PS)    U(THETA(I),PF)
C                            2      U(THETA(I),PS)    F(THETA(I),PF)
C                            3      F(THETA(I),PS)    F(THETA(I),PF)
C                            4      F(THETA(I),PS)    U(THETA(I),PF)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION U
C                          AND THE RIGHT SIDE F AT A CORNER THEN THE
C                          SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSSSP.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F. IDIMF MUST BE
C                          AT LEAST M+1  .
C
C                        W
C                          A ONE-DIMENSIONAL ARRAY THAT MUST BE
C                          PROVIDED BY THE USER FOR WORK SPACE.
C                          W MAY REQUIRE UP TO
C                          4*(N+1)+(16+INT(LOG2(N+1)))(M+1) LOCATIONS
C                          THE ACTUAL NUMBER OF LOCATIONS USED IS
C                          COMPUTED BY HWSSSP AND IS OUTPUT IN
C                          LOCATION W(1). INT( ) DENOTES THE
C                          FORTRAN INTEGER FUNCTION.
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (THETA(I),PHI(J)),  I = 1,2,...,M+1  AND
C                          J = 1,2,...,N+1  .
C
C                        PERTRB
C                          IF ONE SPECIFIES A COMBINATION OF PERIODIC,
C                          DERIVATIVE OR UNSPECIFIED BOUNDARY
C                          CONDITIONS FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HWSSSP THEN COMPUTES
C                          THIS SOLUTION, WHICH IS A LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION IS NOT UNIQUE AND IS
C                          UNNORMALIZED. THE VALUE OF PERTRB SHOULD
C                          BE SMALL COMPARED TO THE RIGHT SIDE F.
C                          OTHERWISE , A SOLUTION IS OBTAINED TO AN
C                          ESSENTIALLY DIFFERENT PROBLEM. THIS
C                          COMPARISON SHOULD ALWAYS BE MADE TO INSURE
C                          THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 8,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          = 0  NO ERROR
C                          = 1  TS.LT.0 OR TF.GT.PI
C                          = 2  TS.GE.TF
C                          = 3  MBDCND.LT.1 OR MBDCND.GT.9
C                          = 4  PS.LT.0 OR PS.GT.PI+PI
C                          = 5  PS.GE.PF
C                          = 6  N.LT.5
C                          = 7  M.LT.5
C                          = 8  NBDCND.LT.0 OR NBDCND.GT.4
C                          = 9  ELMBDA.GT.0
C                          = 10 IDIMF.LT.M+1
C                          = 11 NBDCND EQUALS 1,2 OR 4 AND MBDCND.GE.5
C                          = 12 TS.EQ.0 AND MBDCND EQUALS 3,4 OR 8
C                          = 13 TF.EQ.PI AND MBDCND EQUALS 2,3 OR 6
C                          = 14 MBDCND EQUALS 5,6 OR 9 AND TS.NE.0
C                          = 15 MBDCND.GE.7 AND TF.NE.PI
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HWSSSP, THE
C                          USER SHOULD TEST IERROR AFTER A CALL.
C
C                        W
C                          CONTAINS INTERMEDIATE VALUES THAT MUST NOT
C                          BE DESTROYED IF HWSSSP WILL BE CALLED AGAIN
C                          WITH INTL = 1. W(1) CONTAINS THE REQUIRED
C                          LENGTH OF W .
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       GENBUN, GNBNAUX, AND COMF
C FILES                  FROM FISHPACK
C
C LANGUAGE               FORTRAN
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980.
C
C PORTABILITY            FORTRAN 77
C
C ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
C                        AND THEN CALLS GENBUN TO SOLVE THE SYSTEM.
C
C TIMING                 FOR LARGE  M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO
C                          M*N*(LOG2(N)
C                        BUT ALSO DEPENDS ON INPUT PARAMETERS NBDCND
C                        AND MBDCND.
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A LOSS
C                        OF NO MORE THAN THREE SIGNIFICANT DIGITS FOR N
C                        AND M AS LARGE AS 64.  MORE DETAILS ABOUT
C                        ACCURACY CAN BE FOUND IN THE DOCUMENTATION FOR
C                        SUBROUTINE GENBUN WHICH IS THE ROUTINE THAT
C                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
C
C REFERENCES             P. N. SWARZTRAUBER, "THE DIRECT SOLUTION OF
C                        THE DISCRETE POISSON EQUATION ON THE SURFACE OF
C                        A SPHERE", S.I.A.M. J. NUMER. ANAL.,15(1974),
C                        PP 212-215.
C
C                        SWARZTRAUBER,P. AND R. SWEET, "EFFICIENT
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
C                        ELLIPTIC EQUATIONS", NCAR TN/IA-109, JULY,
C                        1975, 138 PP.
C***********************************************************************
      DIMENSION       F(IDIMF,1) ,BDTS(*)    ,BDTF(*)    ,BDPS(*)    ,
     1                BDPF(*)    ,W(*)
C
      NBR = NBDCND+1
      PI = PIMACH(DUM)
      TPI = 2.*PI
      IERROR = 0
      IF (TS.LT.0. .OR. TF.GT.PI) IERROR = 1
      IF (TS .GE. TF) IERROR = 2
      IF (MBDCND.LT.1 .OR. MBDCND.GT.9) IERROR = 3
      IF (PS.LT.0. .OR. PF.GT.TPI) IERROR = 4
      IF (PS .GE. PF) IERROR = 5
      IF (N .LT. 5) IERROR = 6
      IF (M .LT. 5) IERROR = 7
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) IERROR = 8
      IF (ELMBDA .GT. 0.) IERROR = 9
      IF (IDIMF .LT. M+1) IERROR = 10
      IF ((NBDCND.EQ.1 .OR. NBDCND.EQ.2 .OR. NBDCND.EQ.4) .AND.
     1    MBDCND.GE.5) IERROR = 11
      IF (TS.EQ.0. .AND.
     1    (MBDCND.EQ.3 .OR. MBDCND.EQ.4 .OR. MBDCND.EQ.8)) IERROR = 12
      IF (TF.EQ.PI .AND.
     1    (MBDCND.EQ.2 .OR. MBDCND.EQ.3 .OR. MBDCND.EQ.6)) IERROR = 13
      IF ((MBDCND.EQ.5 .OR. MBDCND.EQ.6 .OR. MBDCND.EQ.9) .AND.
     1    TS.NE.0.) IERROR = 14
      IF (MBDCND.GE.7 .AND. TF.NE.PI) IERROR = 15
      IF (IERROR.NE.0 .AND. IERROR.NE.9) RETURN
      CALL HWSSS1 (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,BDPF,
     1             ELMBDA,F,IDIMF,PERTRB,W,W(M+2),W(2*M+3),W(3*M+4),
     2             W(4*M+5),W(5*M+6),W(6*M+7))
      W(1) = W(6*M+7)+FLOAT(6*(M+1))
      RETURN
      END
      SUBROUTINE HWSSS1 (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,
     1                   BDPF,ELMBDA,F,IDIMF,PERTRB,AM,BM,CM,SN,SS,
     2                   SINT,D)
      DIMENSION       F(IDIMF,*) ,BDTS(*)    ,BDTF(*)    ,BDPS(*)    ,
     1                BDPF(*)    ,AM(*)      ,BM(*)      ,CM(*)      ,
     2                SS(*)      ,SN(*)      ,D(*)       ,SINT(*)
C
      PI = PIMACH(DUM)
      TPI = PI+PI
      HPI = PI/2.
      MP1 = M+1
      NP1 = N+1
      FN = N
      FM = M
      DTH = (TF-TS)/FM
      HDTH = DTH/2.
      TDT = DTH+DTH
      DPHI = (PF-PS)/FN
      TDP = DPHI+DPHI
      DPHI2 = DPHI*DPHI
      EDP2 = ELMBDA*DPHI2
      DTH2 = DTH*DTH
      CP = 4./(FN*DTH2)
      WP = FN*SIN(HDTH)/4.
      DO 102 I=1,MP1
         FIM1 = I-1
         THETA = FIM1*DTH+TS
         SINT(I) = SIN(THETA)
         IF (SINT(I)) 101,102,101
  101    T1 = 1./(DTH2*SINT(I))
         AM(I) = T1*SIN(THETA-HDTH)
         CM(I) = T1*SIN(THETA+HDTH)
         BM(I) = -AM(I)-CM(I)+ELMBDA
  102 CONTINUE
      INP = 0
      ISP = 0
C
C BOUNDARY CONDITION AT THETA=TS
C
      MBR = MBDCND+1
      GO TO (103,104,104,105,105,106,106,104,105,106),MBR
  103 ITS = 1
      GO TO 107
  104 AT = AM(2)
      ITS = 2
      GO TO 107
  105 AT = AM(1)
      ITS = 1
      CM(1) = AM(1)+CM(1)
      GO TO 107
  106 AT = AM(2)
      INP = 1
      ITS = 2
C
C BOUNDARY CONDITION THETA=TF
C
  107 GO TO (108,109,110,110,109,109,110,111,111,111),MBR
  108 ITF = M
      GO TO 112
  109 CT = CM(M)
      ITF = M
      GO TO 112
  110 CT = CM(M+1)
      AM(M+1) = AM(M+1)+CM(M+1)
      ITF = M+1
      GO TO 112
  111 ITF = M
      ISP = 1
      CT = CM(M)
C
C COMPUTE HOMOGENEOUS SOLUTION WITH SOLUTION AT POLE EQUAL TO ONE
C
  112 ITSP = ITS+1
      ITFM = ITF-1
      WTS = SINT(ITS+1)*AM(ITS+1)/CM(ITS)
      WTF = SINT(ITF-1)*CM(ITF-1)/AM(ITF)
      MUNK = ITF-ITS+1
      IF (ISP) 116,116,113
  113 D(ITS) = CM(ITS)/BM(ITS)
      DO 114 I=ITSP,M
         D(I) = CM(I)/(BM(I)-AM(I)*D(I-1))
  114 CONTINUE
      SS(M) = -D(M)
      IID = M-ITS
      DO 115 II=1,IID
         I = M-II
         SS(I) = -D(I)*SS(I+1)
  115 CONTINUE
      SS(M+1) = 1.
  116 IF (INP) 120,120,117
  117 SN(1) = 1.
      D(ITF) = AM(ITF)/BM(ITF)
      IID = ITF-2
      DO 118 II=1,IID
         I = ITF-II
         D(I) = AM(I)/(BM(I)-CM(I)*D(I+1))
  118 CONTINUE
      SN(2) = -D(2)
      DO 119 I=3,ITF
         SN(I) = -D(I)*SN(I-1)
  119 CONTINUE
C
C BOUNDARY CONDITIONS AT PHI=PS
C
  120 NBR = NBDCND+1
      WPS = 1.
      WPF = 1.
      GO TO (121,122,122,123,123),NBR
  121 JPS = 1
      GO TO 124
  122 JPS = 2
      GO TO 124
  123 JPS = 1
      WPS = .5
C
C BOUNDARY CONDITION AT PHI=PF
C
  124 GO TO (125,126,127,127,126),NBR
  125 JPF = N
      GO TO 128
  126 JPF = N
      GO TO 128
  127 WPF = .5
      JPF = N+1
  128 JPSP = JPS+1
      JPFM = JPF-1
      NUNK = JPF-JPS+1
      FJJ = JPFM-JPSP+1
C
C SCALE COEFFICIENTS FOR SUBROUTINE GENBUN
C
      DO 129 I=ITS,ITF
         CF = DPHI2*SINT(I)*SINT(I)
         AM(I) = CF*AM(I)
         BM(I) = CF*BM(I)
         CM(I) = CF*CM(I)
  129 CONTINUE
      AM(ITS) = 0.
      CM(ITF) = 0.
      ISING = 0
      GO TO (130,138,138,130,138,138,130,138,130,130),MBR
  130 GO TO (131,138,138,131,138),NBR
  131 IF (ELMBDA) 138,132,132
  132 ISING = 1
      SUM = WTS*WPS+WTS*WPF+WTF*WPS+WTF*WPF
      IF (INP) 134,134,133
  133 SUM = SUM+WP
  134 IF (ISP) 136,136,135
  135 SUM = SUM+WP
  136 SUM1 = 0.
      DO 137 I=ITSP,ITFM
         SUM1 = SUM1+SINT(I)
  137 CONTINUE
      SUM = SUM+FJJ*(SUM1+WTS+WTF)
      SUM = SUM+(WPS+WPF)*SUM1
      HNE = SUM
  138 GO TO (146,142,142,144,144,139,139,142,144,139),MBR
  139 IF (NBDCND-3) 146,140,146
  140 YHLD = F(1,JPS)-4./(FN*DPHI*DTH2)*(BDPF(2)-BDPS(2))
      DO 141 J=1,NP1
         F(1,J) = YHLD
  141 CONTINUE
      GO TO 146
  142 DO 143 J=JPS,JPF
         F(2,J) = F(2,J)-AT*F(1,J)
  143 CONTINUE
      GO TO 146
  144 DO 145 J=JPS,JPF
         F(1,J) = F(1,J)+TDT*BDTS(J)*AT
  145 CONTINUE
  146 GO TO (154,150,152,152,150,150,152,147,147,147),MBR
  147 IF (NBDCND-3) 154,148,154
  148 YHLD = F(M+1,JPS)-4./(FN*DPHI*DTH2)*(BDPF(M)-BDPS(M))
      DO 149 J=1,NP1
         F(M+1,J) = YHLD
  149 CONTINUE
      GO TO 154
  150 DO 151 J=JPS,JPF
         F(M,J) = F(M,J)-CT*F(M+1,J)
  151 CONTINUE
      GO TO 154
  152 DO 153 J=JPS,JPF
         F(M+1,J) = F(M+1,J)-TDT*BDTF(J)*CT
  153 CONTINUE
  154 GO TO (159,155,155,157,157),NBR
  155 DO 156 I=ITS,ITF
         F(I,2) = F(I,2)-F(I,1)/(DPHI2*SINT(I)*SINT(I))
  156 CONTINUE
      GO TO 159
  157 DO 158 I=ITS,ITF
         F(I,1) = F(I,1)+TDP*BDPS(I)/(DPHI2*SINT(I)*SINT(I))
  158 CONTINUE
  159 GO TO (164,160,162,162,160),NBR
  160 DO 161 I=ITS,ITF
         F(I,N) = F(I,N)-F(I,N+1)/(DPHI2*SINT(I)*SINT(I))
  161 CONTINUE
      GO TO 164
  162 DO 163 I=ITS,ITF
         F(I,N+1) = F(I,N+1)-TDP*BDPF(I)/(DPHI2*SINT(I)*SINT(I))
  163 CONTINUE
  164 CONTINUE
      PERTRB = 0.
      IF (ISING) 165,176,165
  165 SUM = WTS*WPS*F(ITS,JPS)+WTS*WPF*F(ITS,JPF)+WTF*WPS*F(ITF,JPS)+
     1      WTF*WPF*F(ITF,JPF)
      IF (INP) 167,167,166
  166 SUM = SUM+WP*F(1,JPS)
  167 IF (ISP) 169,169,168
  168 SUM = SUM+WP*F(M+1,JPS)
  169 DO 171 I=ITSP,ITFM
         SUM1 = 0.
         DO 170 J=JPSP,JPFM
            SUM1 = SUM1+F(I,J)
  170    CONTINUE
         SUM = SUM+SINT(I)*SUM1
  171 CONTINUE
      SUM1 = 0.
      SUM2 = 0.
      DO 172 J=JPSP,JPFM
         SUM1 = SUM1+F(ITS,J)
         SUM2 = SUM2+F(ITF,J)
  172 CONTINUE
      SUM = SUM+WTS*SUM1+WTF*SUM2
      SUM1 = 0.
      SUM2 = 0.
      DO 173 I=ITSP,ITFM
         SUM1 = SUM1+SINT(I)*F(I,JPS)
         SUM2 = SUM2+SINT(I)*F(I,JPF)
  173 CONTINUE
      SUM = SUM+WPS*SUM1+WPF*SUM2
      PERTRB = SUM/HNE
      DO 175 J=1,NP1
         DO 174 I=1,MP1
            F(I,J) = F(I,J)-PERTRB
  174    CONTINUE
  175 CONTINUE
C
C SCALE RIGHT SIDE FOR SUBROUTINE GENBUN
C
  176 DO 178 I=ITS,ITF
         CF = DPHI2*SINT(I)*SINT(I)
         DO 177 J=JPS,JPF
            F(I,J) = CF*F(I,J)
  177    CONTINUE
  178 CONTINUE
      CALL GENBUN (NBDCND,NUNK,1,MUNK,AM(ITS),BM(ITS),CM(ITS),IDIMF,
     1             F(ITS,JPS),IERROR,D)
      IF (ISING) 186,186,179
  179 IF (INP) 183,183,180
  180 IF (ISP) 181,181,186
  181 DO 182 J=1,NP1
         F(1,J) = 0.
  182 CONTINUE
      GO TO 209
  183 IF (ISP) 186,186,184
  184 DO 185 J=1,NP1
         F(M+1,J) = 0.
  185 CONTINUE
      GO TO 209
  186 IF (INP) 193,193,187
  187 SUM = WPS*F(ITS,JPS)+WPF*F(ITS,JPF)
      DO 188 J=JPSP,JPFM
         SUM = SUM+F(ITS,J)
  188 CONTINUE
      DFN = CP*SUM
      DNN = CP*((WPS+WPF+FJJ)*(SN(2)-1.))+ELMBDA
      DSN = CP*(WPS+WPF+FJJ)*SN(M)
      IF (ISP) 189,189,194
  189 CNP = (F(1,1)-DFN)/DNN
      DO 191 I=ITS,ITF
         HLD = CNP*SN(I)
         DO 190 J=JPS,JPF
            F(I,J) = F(I,J)+HLD
  190    CONTINUE
  191 CONTINUE
      DO 192 J=1,NP1
         F(1,J) = CNP
  192 CONTINUE
      GO TO 209
  193 IF (ISP) 209,209,194
  194 SUM = WPS*F(ITF,JPS)+WPF*F(ITF,JPF)
      DO 195 J=JPSP,JPFM
         SUM = SUM+F(ITF,J)
  195 CONTINUE
      DFS = CP*SUM
      DSS = CP*((WPS+WPF+FJJ)*(SS(M)-1.))+ELMBDA
      DNS = CP*(WPS+WPF+FJJ)*SS(2)
      IF (INP) 196,196,200
  196 CSP = (F(M+1,1)-DFS)/DSS
      DO 198 I=ITS,ITF
         HLD = CSP*SS(I)
         DO 197 J=JPS,JPF
            F(I,J) = F(I,J)+HLD
  197    CONTINUE
  198 CONTINUE
      DO 199 J=1,NP1
         F(M+1,J) = CSP
  199 CONTINUE
      GO TO 209
  200 RTN = F(1,1)-DFN
      RTS = F(M+1,1)-DFS
      IF (ISING) 202,202,201
  201 CSP = 0.
      CNP = RTN/DNN
      GO TO 205
  202 IF (ABS(DNN)-ABS(DSN)) 204,204,203
  203 DEN = DSS-DNS*DSN/DNN
      RTS = RTS-RTN*DSN/DNN
      CSP = RTS/DEN
      CNP = (RTN-CSP*DNS)/DNN
      GO TO 205
  204 DEN = DNS-DSS*DNN/DSN
      RTN = RTN-RTS*DNN/DSN
      CSP = RTN/DEN
      CNP = (RTS-DSS*CSP)/DSN
  205 DO 207 I=ITS,ITF
         HLD = CNP*SN(I)+CSP*SS(I)
         DO 206 J=JPS,JPF
            F(I,J) = F(I,J)+HLD
  206    CONTINUE
  207 CONTINUE
      DO 208 J=1,NP1
         F(1,J) = CNP
         F(M+1,J) = CSP
  208 CONTINUE
  209 IF (NBDCND) 212,210,212
  210 DO 211 I=1,MP1
         F(I,JPF+1) = F(I,JPS)
  211 CONTINUE
  212 RETURN
C
C REVISION HISTORY---
C
C SEPTEMBER 1973    VERSION 1
C APRIL     1976    VERSION 2
C JANUARY   1978    VERSION 3
C DECEMBER  1979    VERSION 3.1
C FEBRUARY  1985    DOCUMENTATION UPGRADE
C NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
C-----------------------------------------------------------------------
      END
