C
C     file hwscsp.f
C
      SUBROUTINE HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,
     1                   BDRS,BDRF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
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
C DIMENSION OF           BDTS(N+1),     BDTF(N+1), BDRS(M+1), BDRF(M+1),
C ARGUMENTS              F(IDIMF,N+1),  W(SEE ARGUMENT LIST)
C
C LATEST REVISION        NOVEMBER 1988
C
C PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION
C                        TO THE MODIFIED HELMHOLTZ EQUATION IN
C                        SPHERICAL COORDINATES ASSUMING AXISYMMETRY
C                        (NO DEPENDENCE ON LONGITUDE).  THE EQUATION
C                        IS
C
C                          (1/R**2)(D/DR)((R**2)(D/DR)U) +
C
C                          (1/(R**2)SIN(THETA))(D/DTHETA)
C
C                          (SIN(THETA)(D/DTHETA)U) +
C
C                          (LAMBDA/(RSIN(THETA))**2)U = F(THETA,R).
C
C                        THIS TWO DIMENSIONAL MODIFIED HELMHOLTZ
C                        EQUATION RESULTS FROM THE FOURIER TRANSFORM
C                        OF THE THREE DIMENSIONAL POISSON EQUATION.
C
C USAGE                  CALL HWSCSP (INTL,TS,TF,M,MBDCND,BDTS,BDTF,
C                                     RS,RF,N,NBDCND,BDRS,BDRF,ELMBDA,
C                                     F,IDIMF,PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT               INTL
C                          = 0  ON INITIAL ENTRY TO HWSCSP OR IF ANY
C                               OF THE ARGUMENTS RS, RF, N, NBDCND
C                               ARE CHANGED FROM A PREVIOUS CALL.
C                          = 1  IF RS, RF, N, NBDCND ARE ALL UNCHANGED
C                               FROM PREVIOUS CALL TO HWSCSP.
C
C                          NOTE:
C                          A CALL WITH INTL=0 TAKES APPROXIMATELY
C                          1.5 TIMES AS MUCH TIME AS A CALL WITH
C                          INTL = 1  .  ONCE A CALL WITH INTL = 0
C                          HAS BEEN MADE THEN SUBSEQUENT SOLUTIONS
C                          CORRESPONDING TO DIFFERENT F, BDTS, BDTF,
C                          BDRS, BDRF CAN BE OBTAINED FASTER WITH
C                          INTL = 1 SINCE INITIALIZATION IS NOT
C                          REPEATED.
C
C                        TS,TF
C                          THE RANGE OF THETA (COLATITUDE), I.E.,
C                          TS .LE. THETA .LE. TF. TS MUST BE LESS
C                          THAN TF.  TS AND TF ARE IN RADIANS. A TS OF
C                          ZERO CORRESPONDS TO THE NORTH POLE AND A
C                          TF OF PI CORRESPONDS TO THE SOUTH POLE.
C
C                          **** IMPORTANT ****
C
C                          IF TF IS EQUAL TO PI THEN IT MUST BE
C                          COMPUTED USING THE STATEMENT
C                          TF = PIMACH(DUM). THIS INSURES THAT TF
C                          IN THE USER'S PROGRAM IS EQUAL TO PI IN
C                          THIS PROGRAM WHICH PERMITS SEVERAL TESTS
C                          OF THE  INPUT PARAMETERS THAT OTHERWISE
C                          WOULD NOT BE POSSIBLE.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (TS,TF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE M+1 GRID POINTS
C                          IN THE THETA-DIRECTION GIVEN BY
C                          THETA(K) = (I-1)DTHETA+TS FOR
C                          I = 1,2,...,M+1, WHERE DTHETA = (TF-TS)/M
C                          IS THE PANEL WIDTH.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT THETA = TS AND  THETA = TF.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THETA = TF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS AND THETA = TF
C                               (SEE NOTES 1,2 BELOW).
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW) AND
C                               SOLUTION IS SPECIFIED AT THETA = TF.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE SOLUTION IS
C                                SPECIFIED AT THETA = TF.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO THETA
C                               IS SPECIFIED AT THETA = TF
C                               (SEE NOTE 2 BELOW).
C                          = 7  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = TS AND THE SOLUTION IS
C                                UNSPECIFIED AT THETA = TF = PI.
C                          = 8  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = TS (SEE NOTE 1 BELOW)
C                               AND THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TF = PI.
C                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
C                               THETA = TS = 0 AND THETA = TF = PI.
C
C                          NOTE 1:
C                          IF TS = 0, DO NOT USE MBDCND = 3,4, OR 8,
C                          BUT INSTEAD USE MBDCND = 5,6, OR 9  .
C
C                          NOTE 2:
C                          IF TF = PI, DO NOT USE MBDCND = 2,3, OR 6,
C                          BUT INSTEAD USE MBDCND = 7,8, OR 9  .
C
C                        BDTS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TS.  WHEN MBDCND = 3,4, OR 8,
C
C                            BDTS(J) = (D/DTHETA)U(TS,R(J)),
C                            J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTS IS
C                          A DUMMY VARIABLE.
C
C                        BDTF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO THETA AT
C                          THETA = TF.  WHEN MBDCND = 2,3, OR 6,
C
C                          BDTF(J) = (D/DTHETA)U(TF,R(J)),
C                          J = 1,2,...,N+1  .
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDTF IS
C                          A DUMMY VARIABLE.
C
C                        RS,RF
C                          THE RANGE OF R, I.E., RS .LE. R .LT. RF.
C                          RS MUST BE LESS THAN RF.  RS MUST BE
C                          NON-NEGATIVE.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (RS,RF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS IN THE
C                          R-DIRECTION GIVEN BY R(J) = (J-1)DR+RS
C                          FOR J = 1,2,...,N+1, WHERE DR = (RF-RS)/N
C                          IS THE PANEL WIDTH.
C                          N MUST BE GREATER THAN 2
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT R = RS AND R = RF.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               R = RS AND R = RF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               R = RS AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO R
C                               IS SPECIFIED AT R = RF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = RS AND R = RF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               RS AND THE SOLUTION IS SPECIFIED AT
C                               R = RF.
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = RS = 0 (SEE NOTE BELOW)  AND THE
C                               SOLUTION IS SPECIFIED AT R = RF.
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = RS = 0 (SEE NOTE BELOW) AND THE
C                               DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO R IS SPECIFIED AT R = RF.
C
C                          NOTE:
C                          NBDCND = 5 OR 6 CANNOT BE USED WITH
C                          MBDCND = 1,2,4,5, OR 7.  THE FORMER
C                          INDICATES THAT THE SOLUTION IS UNSPECIFIED
C                          AT R = 0, THE LATTER INDICATES THAT THE
C                          SOLUTION IS SPECIFIED).
C                          USE INSTEAD   NBDCND = 1 OR 2  .
C
C                        BDRS
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
C                          THE SOLUTION WITH RESPECT TO R AT R = RS.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDRS(I) = (D/DR)U(THETA(I),RS),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDRS IS
C                          A DUMMY VARIABLE.
C
C                        BDRF
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
C                          THAT SPECIFIES THE VALUES OF THE
C                          DERIVATIVE OF THE SOLUTION WITH RESPECT
C                          TO R AT R = RF.
C
C                          WHEN NBDCND = 2,3, OR 6,
C                            BDRF(I) = (D/DR)U(THETA(I),RF),
C                            I = 1,2,...,M+1  .
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDRF IS
C                          A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HWSCSP WILL
C                          ATTEMPT TO FIND A SOLUTION.  IF NBDCND = 5
C                          OR 6 OR  MBDCND = 5,6,7,8, OR 9, ELMBDA
C                          MUST BE ZERO.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
C                          LEAST (M+1)*(N+1), SPECIFYING VALUES OF THE
C                          RIGHT SIDE OF THE HELMHOLTZ EQUATION AND
C                          BOUNDARY VALUES (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I = 2,3,...,M AND J = 2,3,...,N
C                          F(I,J) = F(THETA(I),R(J)).
C
C                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
C                          FOR J=1,2,...,N+1,  I=1,2,...,M+1,
C
C                          MBDCND   F(1,J)            F(M+1,J)
C                          ------   ----------        ----------
C
C                            1      U(TS,R(J))        U(TF,R(J))
C                            2      U(TS,R(J))        F(TF,R(J))
C                            3      F(TS,R(J))        F(TF,R(J))
C                            4      F(TS,R(J))        U(TF,R(J))
C                            5      F(0,R(J))         U(TF,R(J))
C                            6      F(0,R(J))         F(TF,R(J))
C                            7      U(TS,R(J))        F(PI,R(J))
C                            8      F(TS,R(J))        F(PI,R(J))
C                            9      F(0,R(J))         F(PI,R(J))
C
C                            NBDCND   F(I,1)            F(I,N+1)
C                            ------   --------------    --------------
C
C                              1      U(THETA(I),RS)    U(THETA(I),RF)
C                              2      U(THETA(I),RS)    F(THETA(I),RF)
C                              3      F(THETA(I),RS)    F(THETA(I),RF)
C                              4      F(THETA(I),RS)    U(THETA(I),RF)
C                              5      F(TS,0)           U(THETA(I),RF)
C                              6      F(TS,0)           F(THETA(I),RF)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
C                          U AND THE RIGHT SIDE F AT A CORNER THEN
C                          THE SOLUTION MUST BE SPECIFIED.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HWSCSP.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M+1  .
C
C                        W
C                          A ONE-DIMENSIONAL ARRAY THAT MUST BE
C                          PROVIDED BY THE USER FOR WORK SPACE.
C                          ITS LENGTH CAN BE COMPUTED FROM THE
C                          FORMULA BELOW WHICH DEPENDS ON THE VALUE
C                          OF NBDCND
C
C                          IF NBDCND=2,4 OR 6 DEFINE NUNK=N
C                          IF NBDCND=1 OR 5   DEFINE NUNK=N-1
C                          IF NBDCND=3        DEFINE NUNK=N+1
C
C                          NOW SET K=INT(LOG2(NUNK))+1 AND
C                          L=2**(K+1) THEN W MUST BE DIMENSIONED
C                          AT LEAST (K-2)*L+K+5*(M+N)+MAX(2*N,6*M)+23
C
C                          **IMPORTANT**
C                          FOR PURPOSES OF CHECKING, THE REQUIRED
C                          LENGTH OF W IS COMPUTED BY HWSCSP AND
C                          STORED IN W(1) IN FLOATING POINT FORMAT.
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (THETA(I),R(J)),  I = 1,2,...,M+1,
C                                            J = 1,2,...,N+1  .
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
C                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
C                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
C                          CALCULATED AND SUBTRACTED FROM F, WHICH
C                          ENSURES THAT A SOLUTION EXISTS.  HWSCSP
C                          THEN COMPUTES THIS SOLUTION, WHICH IS A
C                          LEAST SQUARES SOLUTION TO THE ORIGINAL
C                          APPROXIMATION. THIS SOLUTION IS NOT UNIQUE
C                          AND IS UNNORMALIZED. THE VALUE OF PERTRB
C                          SHOULD BE SMALL COMPARED TO THE RIGHT SIDE
C                          F. OTHERWISE , A SOLUTION IS OBTAINED TO
C                          AN ESSENTIALLY DIFFERENT PROBLEM. THIS
C                          COMPARISON SHOULD ALWAYS BE MADE TO INSURE
C                          THAT A MEANINGFUL SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 10,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          = 1  TS.LT.0. OR TF.GT.PI
C                          = 2  TS.GE.TF
C                          = 3  M.LT.5
C                          = 4  MBDCND.LT.1 OR MBDCND.GT.9
C                          = 5  RS.LT.0
C                          = 6  RS.GE.RF
C                          = 7  N.LT.5
C                          = 8  NBDCND.LT.1 OR NBDCND.GT.6
C                          = 9  ELMBDA.GT.0
C                          = 10 IDIMF.LT.M+1
C                          = 11 ELMBDA.NE.0 AND MBDCND.GE.5
C                          = 12 ELMBDA.NE.0 AND NBDCND EQUALS 5 OR 6
C                          = 13 MBDCND EQUALS 5,6 OR 9 AND TS.NE.0
C                          = 14 MBDCND.GE.7 AND TF.NE.PI
C                          = 15 TS.EQ.0 AND MBDCND EQUALS 3,4 OR 8
C                          = 16 TF.EQ.PI AND MBDCND EQUALS 2,3 OR 6
C                          = 17 NBDCND.GE.5 AND RS.NE.0
C                          = 18 NBDCND.GE.5 AND MBDCND EQUALS 1,2,4,5 OR
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSLIBY INCORRECT CALL TO HWSCSP, THE
C                          USER SHOULD TEST IERROR AFTER A CALL.
C
C                        W
C                          CONTAINS INTERMEDIATE VALUES THAT MUST NOT
C                          BE DESTROYED IF HWSCSP WILL BE CALLED AGAIN
C                          WITH INTL = 1.  W(1) CONTAINS THE NUMBER
C                          OF LOCATIONS WHICH W MUST HAVE
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       BLKTRI, AND COMF FROM FISHPACK
C FILES
C
C LANGUAGE               FORTRAN
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980.
C
C PORTABILITY            FORTRAN 77.
C
C ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
C                        AND THEN CALLS BLKTRI TO SOLVE THE SYSTEM.
C
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, "EFFICIENT
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
C                        ELLIPTIC EQUATIONS"
C                          NCAR TN/IA-109, JULY, 1975, 138 PP.
C***********************************************************************
      DIMENSION       F(IDIMF,1) ,BDTS(*)    ,BDTF(*)    ,BDRS(*)    ,
     1                BDRF(*)    ,W(*)
C
      PI = PIMACH(DUM)
      IERROR = 0
      IF (TS.LT.0. .OR. TF.GT.PI) IERROR = 1
      IF (TS .GE. TF) IERROR = 2
      IF (M .LT. 5) IERROR = 3
      IF (MBDCND.LT.1 .OR. MBDCND.GT.9) IERROR = 4
      IF (RS .LT. 0.) IERROR = 5
      IF (RS .GE. RF) IERROR = 6
      IF (N .LT. 5) IERROR = 7
      IF (NBDCND.LT.1 .OR. NBDCND.GT.6) IERROR = 8
      IF (ELMBDA .GT. 0.) IERROR = 9
      IF (IDIMF .LT. M+1) IERROR = 10
      IF (ELMBDA.NE.0. .AND. MBDCND.GE.5) IERROR = 11
      IF (ELMBDA.NE.0. .AND. (NBDCND.EQ.5 .OR. NBDCND.EQ.6)) IERROR = 12
      IF ((MBDCND.EQ.5 .OR. MBDCND.EQ.6 .OR. MBDCND.EQ.9) .AND.
     1    TS.NE.0.) IERROR = 13
      IF (MBDCND.GE.7 .AND. TF.NE.PI) IERROR = 14
      IF (TS.EQ.0. .AND.
     1    (MBDCND.EQ.4 .OR. MBDCND.EQ.8 .OR. MBDCND.EQ.3)) IERROR = 15
      IF (TF.EQ.PI .AND.
     1    (MBDCND.EQ.2 .OR. MBDCND.EQ.3 .OR. MBDCND.EQ.6)) IERROR = 16
      IF (NBDCND.GE.5 .AND. RS.NE.0.) IERROR = 17
      IF (NBDCND.GE.5 .AND. (MBDCND.EQ.1 .OR. MBDCND.EQ.2 .OR.
     1                                    MBDCND.EQ.5 .OR. MBDCND.EQ.7))
     2    IERROR = 18
      IF (IERROR.NE.0 .AND. IERROR.NE.9) RETURN
      NCK = N
      GO TO (101,103,102,103,101,103),NBDCND
  101 NCK = NCK-1
      GO TO 103
  102 NCK = NCK+1
  103 L = 2
      K = 1
  104 L = L+L
      K = K+1
      IF (NCK-L) 105,105,104
  105 L = L+L
      NP1 = N+1
      MP1 = M+1
      I1 = (K-2)*L+K+MAX0(2*N,6*M)+13
      I2 = I1+NP1
      I3 = I2+NP1
      I4 = I3+NP1
      I5 = I4+NP1
      I6 = I5+NP1
      I7 = I6+MP1
      I8 = I7+MP1
      I9 = I8+MP1
      I10 = I9+MP1
      W(1) = FLOAT(I10+M)
      CALL HWSCS1 (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS,
     1             BDRF,ELMBDA,F,IDIMF,PERTRB,W(2),W(I1),W(I2),W(I3),
     2             W(I4),W(I5),W(I6),W(I7),W(I8),W(I9),W(I10))
      RETURN
      END
      SUBROUTINE HWSCS1 (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,
     1                   BDRS,BDRF,ELMBDA,F,IDIMF,PERTRB,W,S,AN,BN,CN,
     2                   R,AM,BM,CM,SINT,BMH)
      DIMENSION       F(IDIMF,*) ,BDRS(*)    ,BDRF(*)    ,BDTS(*)    ,
     1                BDTF(*)    ,AM(*)      ,BM(*)      ,CM(*)      ,
     2                AN(*)      ,BN(*)      ,CN(*)      ,S(*)       ,
     3                R(*)       ,SINT(*)    ,BMH(*)     ,W(*)
      PI = PIMACH(DUM)
      EPS = EPMACH(DUM)
      MP1 = M+1
      DTH = (TF-TS)/FLOAT(M)
      TDT = DTH+DTH
      HDTH = DTH/2.
      SDTS = 1./(DTH*DTH)
      DO 102 I=1,MP1
         THETA = TS+FLOAT(I-1)*DTH
         SINT(I) = SIN(THETA)
         IF (SINT(I)) 101,102,101
  101    T1 = SDTS/SINT(I)
         AM(I) = T1*SIN(THETA-HDTH)
         CM(I) = T1*SIN(THETA+HDTH)
         BM(I) = -(AM(I)+CM(I))
  102 CONTINUE
      NP1 = N+1
      DR = (RF-RS)/FLOAT(N)
      HDR = DR/2.
      TDR = DR+DR
      DR2 = DR*DR
      CZR = 6.*DTH/(DR2*(COS(TS)-COS(TF)))
      DO 103 J=1,NP1
         R(J) = RS+FLOAT(J-1)*DR
         AN(J) = (R(J)-HDR)**2/DR2
         CN(J) = (R(J)+HDR)**2/DR2
         BN(J) = -(AN(J)+CN(J))
  103 CONTINUE
      MP = 1
      NP = 1
C
C BOUNDARY CONDITION AT PHI=PS
C
      GO TO (104,104,105,105,106,106,104,105,106),MBDCND
  104 AT = AM(2)
      ITS = 2
      GO TO 107
  105 AT = AM(1)
      ITS = 1
      CM(1) = CM(1)+AM(1)
      GO TO 107
  106 ITS = 1
      BM(1) = -4.*SDTS
      CM(1) = -BM(1)
C
C BOUNDARY CONDITION AT PHI=PF
C
  107 GO TO (108,109,109,108,108,109,110,110,110),MBDCND
  108 CT = CM(M)
      ITF = M
      GO TO 111
  109 CT = CM(M+1)
      AM(M+1) = AM(M+1)+CM(M+1)
      ITF = M+1
      GO TO 111
  110 ITF = M+1
      AM(M+1) = 4.*SDTS
      BM(M+1) = -AM(M+1)
  111 WTS = SINT(ITS+1)*AM(ITS+1)/CM(ITS)
      WTF = SINT(ITF-1)*CM(ITF-1)/AM(ITF)
      ITSP = ITS+1
      ITFM = ITF-1
C
C BOUNDARY CONDITION AT R=RS
C
      ICTR = 0
      GO TO (112,112,113,113,114,114),NBDCND
  112 AR = AN(2)
      JRS = 2
      GO TO 118
  113 AR = AN(1)
      JRS = 1
      CN(1) = CN(1)+AN(1)
      GO TO 118
  114 JRS = 2
      ICTR = 1
      S(N) = AN(N)/BN(N)
      DO 115 J=3,N
         L = N-J+2
         S(L) = AN(L)/(BN(L)-CN(L)*S(L+1))
  115 CONTINUE
      S(2) = -S(2)
      DO 116 J=3,N
         S(J) = -S(J)*S(J-1)
  116 CONTINUE
      WTNM = WTS+WTF
      DO 117 I=ITSP,ITFM
         WTNM = WTNM+SINT(I)
  117 CONTINUE
      YPS = CZR*WTNM*(S(2)-1.)
C
C BOUNDARY CONDITION AT R=RF
C
  118 GO TO (119,120,120,119,119,120),NBDCND
  119 CR = CN(N)
      JRF = N
      GO TO 121
  120 CR = CN(N+1)
      AN(N+1) = AN(N+1)+CN(N+1)
      JRF = N+1
  121 WRS = AN(JRS+1)*R(JRS)**2/CN(JRS)
      WRF = CN(JRF-1)*R(JRF)**2/AN(JRF)
      WRZ = AN(JRS)/CZR
      JRSP = JRS+1
      JRFM = JRF-1
      MUNK = ITF-ITS+1
      NUNK = JRF-JRS+1
      DO 122 I=ITS,ITF
         BMH(I) = BM(I)
  122 CONTINUE
      ISING = 0
      GO TO (132,132,123,132,132,123),NBDCND
  123 GO TO (132,132,124,132,132,124,132,124,124),MBDCND
  124 IF (ELMBDA) 132,125,125
  125 ISING = 1
      SUM = WTS*WRS+WTS*WRF+WTF*WRS+WTF*WRF
      IF (ICTR) 126,127,126
  126 SUM = SUM+WRZ
  127 DO 129 J=JRSP,JRFM
         R2 = R(J)**2
         DO 128 I=ITSP,ITFM
            SUM = SUM+R2*SINT(I)
  128    CONTINUE
  129 CONTINUE
      DO 130 J=JRSP,JRFM
         SUM = SUM+(WTS+WTF)*R(J)**2
  130 CONTINUE
      DO 131 I=ITSP,ITFM
         SUM = SUM+(WRS+WRF)*SINT(I)
  131 CONTINUE
      HNE = SUM
  132 GO TO (133,133,133,133,134,134,133,133,134),MBDCND
  133 BM(ITS) = BMH(ITS)+ELMBDA/SINT(ITS)**2
  134 GO TO (135,135,135,135,135,135,136,136,136),MBDCND
  135 BM(ITF) = BMH(ITF)+ELMBDA/SINT(ITF)**2
  136 DO 137 I=ITSP,ITFM
         BM(I) = BMH(I)+ELMBDA/SINT(I)**2
  137 CONTINUE
      GO TO (138,138,140,140,142,142,138,140,142),MBDCND
  138 DO 139 J=JRS,JRF
         F(2,J) = F(2,J)-AT*F(1,J)/R(J)**2
  139 CONTINUE
      GO TO 142
  140 DO 141 J=JRS,JRF
         F(1,J) = F(1,J)+TDT*BDTS(J)*AT/R(J)**2
  141 CONTINUE
  142 GO TO (143,145,145,143,143,145,147,147,147),MBDCND
  143 DO 144 J=JRS,JRF
         F(M,J) = F(M,J)-CT*F(M+1,J)/R(J)**2
  144 CONTINUE
      GO TO 147
  145 DO 146 J=JRS,JRF
         F(M+1,J) = F(M+1,J)-TDT*BDTF(J)*CT/R(J)**2
  146 CONTINUE
  147 GO TO (151,151,153,153,148,148),NBDCND
  148 IF (MBDCND-3) 155,149,155
  149 YHLD = F(ITS,1)-CZR/TDT*(SIN(TF)*BDTF(2)-SIN(TS)*BDTS(2))
      DO 150 I=1,MP1
         F(I,1) = YHLD
  150 CONTINUE
      GO TO 155
  151 RS2 = (RS+DR)**2
      DO 152 I=ITS,ITF
         F(I,2) = F(I,2)-AR*F(I,1)/RS2
  152 CONTINUE
      GO TO 155
  153 DO 154 I=ITS,ITF
         F(I,1) = F(I,1)+TDR*BDRS(I)*AR/RS**2
  154 CONTINUE
  155 GO TO (156,158,158,156,156,158),NBDCND
  156 RF2 = (RF-DR)**2
      DO 157 I=ITS,ITF
         F(I,N) = F(I,N)-CR*F(I,N+1)/RF2
  157 CONTINUE
      GO TO 160
  158 DO 159 I=ITS,ITF
         F(I,N+1) = F(I,N+1)-TDR*BDRF(I)*CR/RF**2
  159 CONTINUE
  160 CONTINUE
      PERTRB = 0.
      IF (ISING) 161,170,161
  161 SUM = WTS*WRS*F(ITS,JRS)+WTS*WRF*F(ITS,JRF)+WTF*WRS*F(ITF,JRS)+
     1      WTF*WRF*F(ITF,JRF)
      IF (ICTR) 162,163,162
  162 SUM = SUM+WRZ*F(ITS,1)
  163 DO 165 J=JRSP,JRFM
         R2 = R(J)**2
         DO 164 I=ITSP,ITFM
            SUM = SUM+R2*SINT(I)*F(I,J)
  164    CONTINUE
  165 CONTINUE
      DO 166 J=JRSP,JRFM
         SUM = SUM+R(J)**2*(WTS*F(ITS,J)+WTF*F(ITF,J))
  166 CONTINUE
      DO 167 I=ITSP,ITFM
         SUM = SUM+SINT(I)*(WRS*F(I,JRS)+WRF*F(I,JRF))
  167 CONTINUE
      PERTRB = SUM/HNE
      DO 169 J=1,NP1
         DO 168 I=1,MP1
            F(I,J) = F(I,J)-PERTRB
  168    CONTINUE
  169 CONTINUE
  170 DO 172 J=JRS,JRF
         RSQ = R(J)**2
         DO 171 I=ITS,ITF
            F(I,J) = RSQ*F(I,J)
  171    CONTINUE
  172 CONTINUE
      IFLG = INTL
  173 CALL BLKTRI (IFLG,NP,NUNK,AN(JRS),BN(JRS),CN(JRS),MP,MUNK,
     1             AM(ITS),BM(ITS),CM(ITS),IDIMF,F(ITS,JRS),IERROR,W)
      IFLG = IFLG+1
      IF (IFLG-1) 174,173,174
  174 IF (NBDCND) 177,175,177
  175 DO 176 I=1,MP1
         F(I,JRF+1) = F(I,JRS)
  176 CONTINUE
  177 IF (MBDCND) 180,178,180
  178 DO 179 J=1,NP1
         F(ITF+1,J) = F(ITS,J)
  179 CONTINUE
  180 XP = 0.
      IF (ICTR) 181,188,181
  181 IF (ISING) 186,182,186
  182 SUM = WTS*F(ITS,2)+WTF*F(ITF,2)
      DO 183 I=ITSP,ITFM
         SUM = SUM+SINT(I)*F(I,2)
  183 CONTINUE
      YPH = CZR*SUM
      XP = (F(ITS,1)-YPH)/YPS
      DO 185 J=JRS,JRF
         XPS = XP*S(J)
         DO 184 I=ITS,ITF
            F(I,J) = F(I,J)+XPS
  184    CONTINUE
  185 CONTINUE
  186 DO 187 I=1,MP1
         F(I,1) = XP
  187 CONTINUE
  188 RETURN
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
