C
C     file hw3crt.f
C
      SUBROUTINE HW3CRT (XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,MBDCND,BDYS,
     1                   BDYF,ZS,ZF,N,NBDCND,BDZS,BDZF,ELMBDA,LDIMF,
     2                   MDIMF,F,PERTRB,IERROR,W)
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
C DIMENSION OF           BDXS(MDIMF,N+1),    BDXF(MDIMF,N+1),
C ARGUMENTS              BDYS(LDIMF,N+1),    BDYF(LDIMF,N+1),
C                        BDZS(LDIMF,M+1),    BDZF(LDIMF,M+1),
C                        F(LDIMF,MDIMF,N+1), W(SEE ARGUMENT LIST)
C
C LATEST REVISION        NOVEMBER 1988
C
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
C                        EQUATION IN CARTESIAN COORDINATES.  THIS
C                        EQUATION IS
C
C                          (D/DX)(DU/DX) + (D/DY)(DU/DY) +
C                          (D/DZ)(DU/DZ) + LAMBDA*U = F(X,Y,Z) .
C
C USAGE                  CALL HW3CRT (XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,
C                                     MBDCND,BDYS,BDYF,ZS,ZF,N,NBDCND,
C                                     BDZS,BDZF,ELMBDA,LDIMF,MDIMF,F,
C                                     PERTRB,IERROR,W)
C
C ARGUMENTS
C
C ON INPUT               XS,XF
C
C                          THE RANGE OF X, I.E. XS .LE. X .LE. XF .
C                          XS MUST BE LESS THAN XF.
C
C                        L
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (XS,XF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE L+1 GRID POINTS
C                          IN THE X-DIRECTION GIVEN BY
C                          X(I) = XS+(I-1)DX FOR I=1,2,...,L+1,
C                          WHERE DX = (XF-XS)/L IS THE PANEL WIDTH.
C                          L MUST BE AT LEAST 5.
C
C                        LBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT X = XS AND X = XF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN X,
C                               I.E. U(L+I,J,K) = U(I,J,K).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               X = XS AND X = XF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               X = XS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO X IS
C                               SPECIFIED AT X = XF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED AT
C                               X = XS AND X = XF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED AT
C                               X = XS AND THE SOLUTION IS SPECIFIED
C                               AT X=XF.
C
C                        BDXS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE DERIVATIVE OF THE SOLUTION
C                          WITH RESPECT TO X AT X = XS.
C
C                          WHEN LBDCND = 3 OR 4,
C
C                            BDXS(J,K) = (D/DX)U(XS,Y(J),Z(K)),
C                            J=1,2,...,M+1,      K=1,2,...,N+1.
C
C                          WHEN LBDCND HAS ANY OTHER VALUE, BDXS
C                          IS A DUMMY VARIABLE. BDXS MUST BE
C                          DIMENSIONED AT LEAST (M+1)*(N+1).
C
C                        BDXF
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE DERIVATIVE OF THE SOLUTION
C                          WITH RESPECT TO X AT X = XF.
C
C                          WHEN LBDCND = 2 OR 3,
C
C                            BDXF(J,K) = (D/DX)U(XF,Y(J),Z(K)),
C                            J=1,2,...,M+1,      K=1,2,...,N+1.
C
C                          WHEN LBDCND HAS ANY OTHER VALUE, BDXF IS
C                          A DUMMY VARIABLE.  BDXF MUST BE
C                          DIMENSIONED AT LEAST (M+1)*(N+1).
C
C                        YS,YF
C                          THE RANGE OF Y, I.E. YS .LE. Y .LE. YF.
C                          YS MUST BE LESS THAN YF.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (YS,YF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE M+1 GRID POINTS IN
C                          THE Y-DIRECTION GIVEN BY Y(J) = YS+(J-1)DY
C                          FOR J=1,2,...,M+1,
C                          WHERE DY = (YF-YS)/M IS THE PANEL WIDTH.
C                          M MUST BE AT LEAST 5.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Y = YS AND Y = YF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E.
C                               U(I,M+J,K) = U(I,J,K).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               Y = YS AND Y = YF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               Y = YS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO Y IS
C                               SPECIFIED AT Y = YF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = YS AND Y = YF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               AT Y = YS AND THE SOLUTION IS
C                               SPECIFIED AT Y=YF.
C
C                        BDYS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Y AT Y = YS.
C
C                          WHEN MBDCND = 3 OR 4,
C
C                            BDYS(I,K) = (D/DY)U(X(I),YS,Z(K)),
C                            I=1,2,...,L+1,      K=1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDYS
C                          IS A DUMMY VARIABLE. BDYS MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(N+1).
C
C                        BDYF
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Y AT Y = YF.
C
C                          WHEN MBDCND = 2 OR 3,
C
C                            BDYF(I,K) = (D/DY)U(X(I),YF,Z(K)),
C                            I=1,2,...,L+1,      K=1,2,...,N+1.
C
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDYF
C                          IS A DUMMY VARIABLE. BDYF MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(N+1).
C
C                        ZS,ZF
C                          THE RANGE OF Z, I.E. ZS .LE. Z .LE. ZF.
C                          ZS MUST BE LESS THAN ZF.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (ZS,ZF) IS SUBDIVIDED.
C                          HENCE, THERE WILL BE N+1 GRID POINTS
C                          IN THE Z-DIRECTION GIVEN BY
C                          Z(K) = ZS+(K-1)DZ FOR K=1,2,...,N+1,
C                          WHERE DZ = (ZF-ZS)/N IS THE PANEL WIDTH.
C                          N MUST BE AT LEAST 5.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Z = ZS AND Z = ZF.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.
C                               U(I,J,N+K) = U(I,J,K).
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               Z = ZS AND Z = ZF.
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               Z = ZS AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO Z IS
C                               SPECIFIED AT Z = ZF.
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Z IS SPECIFIED AT
C                               Z = ZS AND Z = ZF.
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Z IS SPECIFIED AT
C                               Z = ZS AND THE SOLUTION IS SPECIFIED
C                               AT Z=ZF.
C
C                        BDZS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Z AT Z = ZS.
C
C                          WHEN NBDCND = 3 OR 4,
C
C                            BDZS(I,J) = (D/DZ)U(X(I),Y(J),ZS),
C                            I=1,2,...,L+1,      J=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDZS
C                          IS A DUMMY VARIABLE. BDZS MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(M+1).
C
C                        BDZF
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE DERIVATIVE OF THE
C                          SOLUTION WITH RESPECT TO Z AT Z = ZF.
C
C                          WHEN NBDCND = 2 OR 3,
C
C                            BDZF(I,J) = (D/DZ)U(X(I),Y(J),ZF),
C                            I=1,2,...,L+1,      J=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDZF
C                          IS A DUMMY VARIABLE. BDZF MUST BE
C                          DIMENSIONED AT LEAST (L+1)*(M+1).
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION. IF LAMBDA .GT. 0, A SOLUTION
C                          MAY NOT EXIST.  HOWEVER, HW3CRT WILL
C                          ATTEMPT TO FIND A SOLUTION.
C
C                        LDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE
C                          ARRAYS F,BDYS,BDYF,BDZS,AND BDZF AS IT
C                          APPEARS IN THE PROGRAM CALLING HW3CRT.
C                          THIS PARAMETER IS USED TO SPECIFY THE
C                          VARIABLE DIMENSION OF THESE ARRAYS.
C                          LDIMF MUST BE AT LEAST L+1.
C
C                        MDIMF
C                          THE COLUMN (OR SECOND) DIMENSION OF THE
C                          ARRAY F AND THE ROW (OR FIRST) DIMENSION
C                          OF THE ARRAYS BDXS AND BDXF AS IT APPEARS
C                          IN THE PROGRAM CALLING HW3CRT.  THIS
C                          PARAMETER IS USED TO SPECIFY THE VARIABLE
C                          DIMENSION OF THESE ARRAYS.
C                          MDIMF MUST BE AT LEAST M+1.
C
C                        F
C                          A THREE-DIMENSIONAL ARRAY OF DIMENSION AT
C                          AT LEAST (L+1)*(M+1)*(N+1), SPECIFYING THE
C                          VALUES OF THE RIGHT SIDE OF THE HELMHOLZ
C                          EQUATION AND BOUNDARY VALUES (IF ANY).
C
C                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
C                          FOR I=2,3,...,L,  J=2,3,...,M,
C                          AND K=2,3,...,N
C                          F(I,J,K) = F(X(I),Y(J),Z(K)).
C
C                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
C                          FOR J=1,2,...,M+1,  K=1,2,...,N+1,
C                          AND I=1,2,...,L+1
C
C                          LBDCND      F(1,J,K)         F(L+1,J,K)
C                          ------   ---------------   ---------------
C
C                            0      F(XS,Y(J),Z(K))   F(XS,Y(J),Z(K))
C                            1      U(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
C                            2      U(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))
C                            3      F(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))
C                            4      F(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
C
C                          MBDCND      F(I,1,K)         F(I,M+1,K)
C                          ------   ---------------   ---------------
C
C                            0      F(X(I),YS,Z(K))   F(X(I),YS,Z(K))
C                            1      U(X(I),YS,Z(K))   U(X(I),YF,Z(K))
C                            2      U(X(I),YS,Z(K))   F(X(I),YF,Z(K))
C                            3      F(X(I),YS,Z(K))   F(X(I),YF,Z(K))
C                            4      F(X(I),YS,Z(K))   U(X(I),YF,Z(K))
C
C                          NBDCND      F(I,J,1)         F(I,J,N+1)
C                          ------   ---------------   ---------------
C
C                            0      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZS)
C                            1      U(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
C                            2      U(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)
C                            3      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)
C                            4      F(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
C
C                          NOTE:
C                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
C                          U AND THE RIGHT SIDE F ON A BOUNDARY,
C                          THEN THE SOLUTION MUST BE SPECIFIED.
C
C                        W
C                          A ONE-DIMENSIONAL ARRAY THAT MUST BE
C                          PROVIDED BY THE USER FOR WORK SPACE.
C                          THE LENGTH OF W MUST BE AT LEAST
C                          30 + L + M + 5*N + MAX(L,M,N) +
C                          7*(INT((L+1)/2) + INT((M+1)/2))
C
C
C
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J,K) OF THE
C                          FINITE DIFFERENCE APPROXIMATION FOR THE
C                          GRID POINT (X(I),Y(J),Z(K)) FOR
C                          I=1,2,...,L+1, J=1,2,...,M+1,
C                          AND K=1,2,...,N+1.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
C                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
C                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
C                          CALCULATED AND SUBTRACTED FROM F, WHICH
C                          ENSURES THAT A SOLUTION EXISTS.  PWSCRT
C                          THEN COMPUTES THIS SOLUTION, WHICH IS A
C                          LEAST SQUARES SOLUTION TO THE ORIGINAL
C                          APPROXIMATION.  THIS SOLUTION IS NOT
C                          UNIQUE AND IS UNNORMALIZED.  THE VALUE OF
C                          PERTRB SHOULD BE SMALL COMPARED TO THE
C                          THE RIGHT SIDE F.  OTHERWISE, A SOLUTION
C                          IS OBTAINED TO AN ESSENTIALLY DIFFERENT
C                          PROBLEM.  THIS COMPARISON SHOULD ALWAYS
C                          BE MADE TO INSURE THAT A MEANINGFUL
C                          SOLUTION HAS BEEN OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 12,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C                          =  1  XS .GE. XF
C                          =  2  L .LT. 5
C                          =  3  LBDCND .LT. 0 .OR. LBDCND .GT. 4
C                          =  4  YS .GE. YF
C                          =  5  M .LT. 5
C                          =  6  MBDCND .LT. 0 .OR. MBDCND .GT. 4
C                          =  7  ZS .GE. ZF
C                          =  8  N .LT. 5
C                          =  9  NBDCND .LT. 0 .OR. NBDCND .GT. 4
C                          = 10  LDIMF .LT. L+1
C                          = 11  MDIMF .LT. M+1
C                          = 12  LAMBDA .GT. 0
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HW3CRT, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       POIS3D, FFTPACK, AND COMF FROM FISHPACK
C FILES
C
C LANGUAGE               FORTRAN
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY 1980.
C
C PORTABILITY            FORTRAN 77
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
C                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS AND
C                        THEN CALLS POIS3D TO SOLVE THE SYSTEM.
C
C TIMING                 FOR LARGE L, M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO
C                          L*M*N*(LOG2(L)+LOG2(M)+5),
C                        BUT ALSO DEPENDS ON INPUT PARAMETERS LBDCND
C                        AND MBDCND.
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
C                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
C                        DIGITS FOR L, M AND N AS LARGE AS 32.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        ROUTINE POIS3D WHICH IS THE ROUTINE THAT
C                        ACTUALLY SOLVES THE FINITE DIFFERENCE
C                        EQUATIONS.
C
C REFERENCES             NONE
C***********************************************************************
      DIMENSION       BDXS(MDIMF,*)          ,BDXF(MDIMF,*)          ,
     1                BDYS(LDIMF,*)          ,BDYF(LDIMF,*)          ,
     2                BDZS(LDIMF,*)          ,BDZF(LDIMF,*)          ,
     3                F(LDIMF,MDIMF,*)       ,W(*)
C
C     CHECK FOR INVALID INPUT.
C
      IERROR = 0
      IF (XF .LE. XS) IERROR = 1
      IF (L .LT. 5) IERROR = 2
      IF (LBDCND.LT.0 .OR. LBDCND.GT.4) IERROR = 3
      IF (YF .LE. YS) IERROR = 4
      IF (M .LT. 5) IERROR = 5
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) IERROR = 6
      IF (ZF .LE. ZS) IERROR = 7
      IF (N .LT. 5) IERROR = 8
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) IERROR = 9
      IF (LDIMF .LT. L+1) IERROR = 10
      IF (MDIMF .LT. M+1) IERROR = 11
      IF (IERROR .NE. 0) GO TO 188
      DY = (YF-YS)/FLOAT(M)
      TWBYDY = 2./DY
      C2 = 1./(DY**2)
      MSTART = 1
      MSTOP = M
      MP1 = M+1
      MP = MBDCND+1
      GO TO (104,101,101,102,102),MP
  101 MSTART = 2
  102 GO TO (104,104,103,103,104),MP
  103 MSTOP = MP1
  104 MUNK = MSTOP-MSTART+1
      DZ = (ZF-ZS)/FLOAT(N)
      TWBYDZ = 2./DZ
      NP = NBDCND+1
      C3 = 1./(DZ**2)
      NP1 = N+1
      NSTART = 1
      NSTOP = N
      GO TO (108,105,105,106,106),NP
  105 NSTART = 2
  106 GO TO (108,108,107,107,108),NP
  107 NSTOP = NP1
  108 NUNK = NSTOP-NSTART+1
      LP1 = L+1
      DX = (XF-XS)/FLOAT(L)
      C1 = 1./(DX**2)
      TWBYDX = 2./DX
      LP = LBDCND+1
      LSTART = 1
      LSTOP = L
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      GO TO (122,109,109,112,112),LP
  109 LSTART = 2
      DO 111 J=MSTART,MSTOP
         DO 110 K=NSTART,NSTOP
            F(2,J,K) = F(2,J,K)-C1*F(1,J,K)
  110    CONTINUE
  111 CONTINUE
      GO TO 115
  112 DO 114 J=MSTART,MSTOP
         DO 113 K=NSTART,NSTOP
            F(1,J,K) = F(1,J,K)+TWBYDX*BDXS(J,K)
  113    CONTINUE
  114 CONTINUE
  115 GO TO (122,116,119,119,116),LP
  116 DO 118 J=MSTART,MSTOP
         DO 117 K=NSTART,NSTOP
            F(L,J,K) = F(L,J,K)-C1*F(LP1,J,K)
  117    CONTINUE
  118 CONTINUE
      GO TO 122
  119 LSTOP = LP1
      DO 121 J=MSTART,MSTOP
         DO 120 K=NSTART,NSTOP
            F(LP1,J,K) = F(LP1,J,K)-TWBYDX*BDXF(J,K)
  120    CONTINUE
  121 CONTINUE
  122 LUNK = LSTOP-LSTART+1
C
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
C
      GO TO (136,123,123,126,126),MP
  123 DO 125 I=LSTART,LSTOP
         DO 124 K=NSTART,NSTOP
            F(I,2,K) = F(I,2,K)-C2*F(I,1,K)
  124    CONTINUE
  125 CONTINUE
      GO TO 129
  126 DO 128 I=LSTART,LSTOP
         DO 127 K=NSTART,NSTOP
            F(I,1,K) = F(I,1,K)+TWBYDY*BDYS(I,K)
  127    CONTINUE
  128 CONTINUE
  129 GO TO (136,130,133,133,130),MP
  130 DO 132 I=LSTART,LSTOP
         DO 131 K=NSTART,NSTOP
            F(I,M,K) = F(I,M,K)-C2*F(I,MP1,K)
  131    CONTINUE
  132 CONTINUE
      GO TO 136
  133 DO 135 I=LSTART,LSTOP
         DO 134 K=NSTART,NSTOP
            F(I,MP1,K) = F(I,MP1,K)-TWBYDY*BDYF(I,K)
  134    CONTINUE
  135 CONTINUE
  136 CONTINUE
C
C     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.
C
      GO TO (150,137,137,140,140),NP
  137 DO 139 I=LSTART,LSTOP
         DO 138 J=MSTART,MSTOP
            F(I,J,2) = F(I,J,2)-C3*F(I,J,1)
  138    CONTINUE
  139 CONTINUE
      GO TO 143
  140 DO 142 I=LSTART,LSTOP
         DO 141 J=MSTART,MSTOP
            F(I,J,1) = F(I,J,1)+TWBYDZ*BDZS(I,J)
  141    CONTINUE
  142 CONTINUE
  143 GO TO (150,144,147,147,144),NP
  144 DO 146 I=LSTART,LSTOP
         DO 145 J=MSTART,MSTOP
            F(I,J,N) = F(I,J,N)-C3*F(I,J,NP1)
  145    CONTINUE
  146 CONTINUE
      GO TO 150
  147 DO 149 I=LSTART,LSTOP
         DO 148 J=MSTART,MSTOP
            F(I,J,NP1) = F(I,J,NP1)-TWBYDZ*BDZF(I,J)
  148    CONTINUE
  149 CONTINUE
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
  150 CONTINUE
      IWB = NUNK+1
      IWC = IWB+NUNK
      IWW = IWC+NUNK
      DO 151 K=1,NUNK
         I = IWC+K-1
         W(K) = C3
         W(I) = C3
         I = IWB+K-1
         W(I) = -2.*C3+ELMBDA
  151 CONTINUE
      GO TO (155,155,153,152,152),NP
  152 W(IWC) = 2.*C3
  153 GO TO (155,155,154,154,155),NP
  154 W(IWB-1) = 2.*C3
  155 CONTINUE
      PERTRB = 0.
C
C     FOR SINGULAR PROBLEMS ADJUST DATA TO INSURE A SOLUTION WILL EXIST.
C
      GO TO (156,172,172,156,172),LP
  156 GO TO (157,172,172,157,172),MP
  157 GO TO (158,172,172,158,172),NP
  158 IF (ELMBDA) 172,160,159
  159 IERROR = 12
      GO TO 172
  160 CONTINUE
      MSTPM1 = MSTOP-1
      LSTPM1 = LSTOP-1
      NSTPM1 = NSTOP-1
      XLP = (2+LP)/3
      YLP = (2+MP)/3
      ZLP = (2+NP)/3
      S1 = 0.
      DO 164 K=2,NSTPM1
         DO 162 J=2,MSTPM1
            DO 161 I=2,LSTPM1
               S1 = S1+F(I,J,K)
  161       CONTINUE
            S1 = S1+(F(1,J,K)+F(LSTOP,J,K))/XLP
  162    CONTINUE
         S2 = 0.
         DO 163 I=2,LSTPM1
            S2 = S2+F(I,1,K)+F(I,MSTOP,K)
  163    CONTINUE
         S2 = (S2+(F(1,1,K)+F(1,MSTOP,K)+F(LSTOP,1,K)+F(LSTOP,MSTOP,K))/
     1                                                          XLP)/YLP
         S1 = S1+S2
  164 CONTINUE
      S = (F(1,1,1)+F(LSTOP,1,1)+F(1,1,NSTOP)+F(LSTOP,1,NSTOP)+
     1    F(1,MSTOP,1)+F(LSTOP,MSTOP,1)+F(1,MSTOP,NSTOP)+
     2                                   F(LSTOP,MSTOP,NSTOP))/(XLP*YLP)
      DO 166 J=2,MSTPM1
         DO 165 I=2,LSTPM1
            S = S+F(I,J,1)+F(I,J,NSTOP)
  165    CONTINUE
  166 CONTINUE
      S2 = 0.
      DO 167 I=2,LSTPM1
         S2 = S2+F(I,1,1)+F(I,1,NSTOP)+F(I,MSTOP,1)+F(I,MSTOP,NSTOP)
  167 CONTINUE
      S = S2/YLP+S
      S2 = 0.
      DO 168 J=2,MSTPM1
         S2 = S2+F(1,J,1)+F(1,J,NSTOP)+F(LSTOP,J,1)+F(LSTOP,J,NSTOP)
  168 CONTINUE
      S = S2/XLP+S
      PERTRB = (S/ZLP+S1)/((FLOAT(LUNK+1)-XLP)*(FLOAT(MUNK+1)-YLP)*
     1                                              (FLOAT(NUNK+1)-ZLP))
      DO 171 I=1,LUNK
         DO 170 J=1,MUNK
            DO 169 K=1,NUNK
               F(I,J,K) = F(I,J,K)-PERTRB
  169       CONTINUE
  170    CONTINUE
  171 CONTINUE
  172 CONTINUE
      NPEROD = 0
      IF (NBDCND .EQ. 0) GO TO 173
      NPEROD = 1
      W(1) = 0.
      W(IWW-1) = 0.
  173 CONTINUE
      CALL POIS3D (LBDCND,LUNK,C1,MBDCND,MUNK,C2,NPEROD,NUNK,W,W(IWB),
     1             W(IWC),LDIMF,MDIMF,F(LSTART,MSTART,NSTART),IR,W(IWW))
C
C     FILL IN SIDES FOR PERIODIC BOUNDARY CONDITIONS.
C
      IF (LP .NE. 1) GO TO 180
      IF (MP .NE. 1) GO TO 175
      DO 174 K=NSTART,NSTOP
         F(1,MP1,K) = F(1,1,K)
  174 CONTINUE
      MSTOP = MP1
  175 IF (NP .NE. 1) GO TO 177
      DO 176 J=MSTART,MSTOP
         F(1,J,NP1) = F(1,J,1)
  176 CONTINUE
      NSTOP = NP1
  177 DO 179 J=MSTART,MSTOP
         DO 178 K=NSTART,NSTOP
            F(LP1,J,K) = F(1,J,K)
  178    CONTINUE
  179 CONTINUE
  180 CONTINUE
      IF (MP .NE. 1) GO TO 185
      IF (NP .NE. 1) GO TO 182
      DO 181 I=LSTART,LSTOP
         F(I,1,NP1) = F(I,1,1)
  181 CONTINUE
      NSTOP = NP1
  182 DO 184 I=LSTART,LSTOP
         DO 183 K=NSTART,NSTOP
            F(I,MP1,K) = F(I,1,K)
  183    CONTINUE
  184 CONTINUE
  185 CONTINUE
      IF (NP .NE. 1) GO TO 188
      DO 187 I=LSTART,LSTOP
         DO 186 J=MSTART,MSTOP
            F(I,J,NP1) = F(I,J,1)
  186    CONTINUE
  187 CONTINUE
  188 CONTINUE
      RETURN
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
