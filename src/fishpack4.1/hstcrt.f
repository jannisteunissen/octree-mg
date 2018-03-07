C
C     file hstcrt.f
C
      SUBROUTINE HSTCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
     1                   ELMBDA,F,IDIMF,PERTRB,IERROR,W)
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
C DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N),
C ARGUMENTS              W(SEE ARGUMENT LIST)
C
C LATEST REVISION        NOVEMBER 1988
C
C PURPOSE                 SOLVES THE STANDARD FIVE-POINT FINITE
C                         DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
C                         EQUATION
C                           (D/DX)(DU/DX) + (D/DY)(DU/DY) + LAMBDA*U
C                           = F(X,Y)
C                         ON A STAGGERED GRID IN CARTESIAN COORDINATES.
C
C USAGE                   CALL HSTCRT (A,B,M,MBDCND,BDA,BDB,C,D
C                                      N,NBDCND,BDC,BDD,ELMBDA,
C                                      F,IDIMF,PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT
C
C                        A,B
C                          THE RANGE OF X, I.E. A .LE. X .LE. B.
C                          A MUST BE LESS THAN B.
C
C                        M
C                          THE NUMBER OF GRID POINTS IN THE
C                          INTERVAL (A,B).  THE GRID POINTS
C                          IN THE X-DIRECTION ARE GIVEN BY
C                          X(I) = A + (I-0.5)DX FOR I=1,2,...,M
C                          WHERE DX =(B-A)/M.  M MUST BE GREATER
C                          THAN 2.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT X = A AND X = B.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN X,
C                               U(M+I,J) = U(I,J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               X = A AND X = B.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               X = A AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO X
C                               IS SPECIFIED AT X = B.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED
C                               AT X = A  AND X = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO X IS SPECIFIED
C                               AT X = A  AND THE SOLUTION IS
C                               SPECIFIED AT X = B.
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N
C                          THAT SPECIFIES THE BOUNDARY VALUES
C                          (IF ANY) OF THE SOLUTION AT X = A.
C
C                          WHEN MBDCND = 1 OR 2,
C                            BDA(J) = U(A,Y(J)) ,         J=1,2,...,N.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DX)U(A,Y(J)) ,   J=1,2,...,N.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N
C                          THAT SPECIFIES THE BOUNDARY VALUES
C                          OF THE SOLUTION AT X = B.
C
C                          WHEN MBDCND = 1 OR 4
C                            BDB(J) = U(B,Y(J)) ,        J=1,2,...,N.
C
C                          WHEN MBDCND = 2 OR 3
C                            BDB(J) = (D/DX)U(B,Y(J)) ,  J=1,2,...,N.
C
C                        C,D
C                          THE RANGE OF Y, I.E. C .LE. Y .LE. D.
C                          C MUST BE LESS THAN D.
C
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE Y-DIRECTION
C                          ARE GIVEN BY Y(J) = C + (J-0.5)DY,
C                          J=1,2,...,N, WHERE DY = (D-C)/N.
C                          N MUST BE GREATER THAN 2.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Y = C   AND Y = D.
C
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E.
C                               U(I,J) = U(I,N+J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT Y = C
C                               AND Y = D.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT Y = C
C                               AND THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = D.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = C AND Y = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO Y IS SPECIFIED AT
C                               Y = C AND THE SOLUTION IS SPECIFIED
C                               AT Y = D.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Y = C.
C
C                          WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(X(I),C) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DY)U(X(I),C),   I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Y = D.
C
C                          WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(X(I),D) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DY)U(X(I),D) ,  I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
C                          A SOLUTION MAY NOT EXIST. HOWEVER,
C                          HSTCRT WILL  ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE RIGHT SIDE OF THE
C                          HELMHOLTZ EQUATION.  FOR I=1,2,...,M
C                          AND J=1,2,...,N
C
C                            F(I,J) = F(X(I),Y(J)) .
C
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTCRT.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.
C                          IDIMF MUST BE AT LEAST M.
C
C                        W
C                          A ONE-DIMENSIONAL ARRAY THAT MUST BE
C                          PROVIDED BY THE USER FOR WORK SPACE.
C                          W MAY REQUIRE UP TO 13M + 4N +
C                          M*INT(LOG2(N)) LOCATIONS. THE ACTUAL NUMBER
C                          OF LOCATIONS USED IS COMPUTED BY HSTCRT
C                          AND IS RETURNED IN THE LOCATION W(1).
C
C ON OUTPUT              F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (X(I),Y(J)) FOR  I=1,2,...,M, J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
C                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
C                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
C                          CALCULATED AND SUBTRACTED FROM F, WHICH
C                          ENSURES THAT A SOLUTION EXISTS.  HSTCRT
C                          THEN COMPUTES THIS SOLUTION, WHICH IS A
C                          LEAST SQUARES SOLUTION TO THE ORIGINAL
C                          APPROXIMATION.  THIS SOLUTION PLUS ANY
C                          CONSTANT IS ALSO A SOLUTION; HENCE, THE
C                          SOLUTION IS NOT UNIQUE.  THE VALUE OF
C                          PERTRB SHOULD BE SMALL COMPARED TO THE
C                          RIGHT SIDE F.  OTHERWISE, A SOLUTION IS
C                          OBTAINED TO AN ESSENTIALLY DIFFERENT PROBLEM.
C                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
C                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT TO NUMBERS 0 AND  6,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C
C                          =  1  A .GE. B
C
C                          =  2  MBDCND .LT. 0 OR MBDCND .GT. 4
C
C                          =  3  C .GE. D
C
C                          =  4  N .LE. 2
C
C                         =  5  NBDCND .LT. 0 OR NBDCND .GT. 4
C
C                         =  6  LAMBDA .GT. 0
C
C                         =  7  IDIMF .LT. M
C
C                         =  8  M .LE. 2
C
C                         SINCE THIS IS THE ONLY MEANS OF INDICATING
C                         A POSSIBLY INCORRECT CALL TO HSTCRT, THE
C                         USER SHOULD TEST IERROR AFTER THE CALL.
C
C                        W
C                          W(1) CONTAINS THE REQUIRED LENGTH OF W.
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       COMF, GENBUN, GNBNAUX, AND POISTG
C FILES                  FROM FISHPACK
C
C LANGUAGE               FORTRAN
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
C                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
C                        IN JANUARY 1980.
C
C PORTABILITY            FORTRAN 77
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-DIFFERENCE
C                        EQUATIONS, INCORPORATES BOUNDARY DATA, ADJUSTS
C                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR
C                        AND CALLS EITHER POISTG OR GENBUN WHICH SOLVES
C                        THE LINEAR SYSTEM OF EQUATIONS.
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A
C                        LOSS OF NO MORE THAN FOUR SIGNIFICANT DIGITS
C                        FOR N AND M AS LARGE AS 64.  MORE DETAILED
C                        INFORMATION ABOUT ACCURACY CAN BE FOUND IN
C                        THE DOCUMENTATION FOR PACKAGE POISTG WHICH
C                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
C
C REFERENCES             U. SCHUMANN AND R. SWEET,"A DIRECT METHOD
C                        FOR THE SOLUTION OF POISSON'S EQUATION WITH
C                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
C                        ARBITRARY SIZE," J. COMP. PHYS. 20(1976),
C                        PP. 171-182.
C***********************************************************************
      DIMENSION       F(IDIMF,1) ,BDA(*)     ,BDB(*)     ,BDC(*)     ,
     1                BDD(*)     ,W(*)
C
C     CHECK FOR INVALID PARAMETERS.
C
      IERROR = 0
      IF (A .GE. B) IERROR = 1
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) IERROR = 2
      IF (C .GE. D) IERROR = 3
      IF (N .LE. 2) IERROR = 4
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) IERROR = 5
      IF (IDIMF .LT. M) IERROR = 7
      IF (M .LE. 2) IERROR = 8
      IF (IERROR .NE. 0) RETURN
      NPEROD = NBDCND
      MPEROD = 0
      IF (MBDCND .GT. 0) MPEROD = 1
      DELTAX = (B-A)/FLOAT(M)
      TWDELX = 1./DELTAX
      DELXSQ = 2./DELTAX**2
      DELTAY = (D-C)/FLOAT(N)
      TWDELY = 1./DELTAY
      DELYSQ = DELTAY**2
      TWDYSQ = 2./DELYSQ
      NP = NBDCND+1
      MP = MBDCND+1
C
C     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
C
      ID2 = M
      ID3 = ID2+M
      ID4 = ID3+M
      S = (DELTAY/DELTAX)**2
      ST2 = 2.*S
      DO 101 I=1,M
         W(I) = S
         J = ID2+I
         W(J) = -ST2+ELMBDA*DELYSQ
         J = ID3+I
         W(J) = S
  101 CONTINUE
C
C     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
C
      GO TO (111,102,102,104,104),MP
  102 DO 103 J=1,N
         F(1,J) = F(1,J)-BDA(J)*DELXSQ
  103 CONTINUE
      W(ID2+1) = W(ID2+1)-W(1)
      GO TO 106
  104 DO 105 J=1,N
         F(1,J) = F(1,J)+BDA(J)*TWDELX
  105 CONTINUE
      W(ID2+1) = W(ID2+1)+W(1)
  106 GO TO (111,107,109,109,107),MP
  107 DO 108 J=1,N
         F(M,J) = F(M,J)-BDB(J)*DELXSQ
  108 CONTINUE
      W(ID3) = W(ID3)-W(1)
      GO TO 111
  109 DO 110 J=1,N
         F(M,J) = F(M,J)-BDB(J)*TWDELX
  110 CONTINUE
      W(ID3) = W(ID3)+W(1)
  111 CONTINUE
C
C     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
C
      GO TO (121,112,112,114,114),NP
  112 DO 113 I=1,M
         F(I,1) = F(I,1)-BDC(I)*TWDYSQ
  113 CONTINUE
      GO TO 116
  114 DO 115 I=1,M
         F(I,1) = F(I,1)+BDC(I)*TWDELY
  115 CONTINUE
  116 GO TO (121,117,119,119,117),NP
  117 DO 118 I=1,M
         F(I,N) = F(I,N)-BDD(I)*TWDYSQ
  118 CONTINUE
      GO TO 121
  119 DO 120 I=1,M
         F(I,N) = F(I,N)-BDD(I)*TWDELY
  120 CONTINUE
  121 CONTINUE
      DO 123 I=1,M
         DO 122 J=1,N
            F(I,J) = F(I,J)*DELYSQ
  122    CONTINUE
  123 CONTINUE
      IF (MPEROD .EQ. 0) GO TO 124
      W(1) = 0.
      W(ID4) = 0.
  124 CONTINUE
      PERTRB = 0.
      IF (ELMBDA) 133,126,125
  125 IERROR = 6
      GO TO 133
  126 GO TO (127,133,133,127,133),MP
  127 GO TO (128,133,133,128,133),NP
C
C     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
C     WILL EXIST.
C
  128 CONTINUE
      S = 0.
      DO 130 J=1,N
         DO 129 I=1,M
            S = S+F(I,J)
  129    CONTINUE
  130 CONTINUE
      PERTRB = S/FLOAT(M*N)
      DO 132 J=1,N
         DO 131 I=1,M
            F(I,J) = F(I,J)-PERTRB
  131    CONTINUE
  132 CONTINUE
      PERTRB = PERTRB/DELYSQ
C
C     SOLVE THE EQUATION.
C
  133 CONTINUE
      IF (NPEROD .EQ. 0) GO TO 134
      CALL POISTG (NPEROD,N,MPEROD,M,W(1),W(ID2+1),W(ID3+1),IDIMF,F,
     1             IERR1,W(ID4+1))
      GO TO 135
  134 CONTINUE
      CALL GENBUN (NPEROD,N,MPEROD,M,W(1),W(ID2+1),W(ID3+1),IDIMF,F,
     1             IERR1,W(ID4+1))
  135 CONTINUE
      W(1) = W(ID4+1)+3.*FLOAT(M)
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
