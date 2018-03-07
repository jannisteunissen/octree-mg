C
C     file hstcyl.f
C
      SUBROUTINE HSTCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
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
C PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
C                        DIFFERENCE APPROXIMATION ON A STAGGERED
C                        GRID TO THE MODIFIED HELMHOLTZ EQUATION
C                        IN CYLINDRICAL COORDINATES. THIS EQUATION
C
C                          (1/R)(D/DR)(R(DU/DR)) + (D/DZ)(DU/DZ)
C
C                          + LAMBDA*(1/R**2)*U = F(R,Z)
C
C                        IS A TWO-DIMENSIONAL MODIFIED HELMHOLTZ
C                        EQUATION RESULTING FROM THE FOURIER TRANSFORM
C                        OF A THREE-DIMENSIONAL POISSON EQUATION.
C
C USAGE                  CALL HSTCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
C                                     PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT               A,B
C
C                          THE RANGE OF R, I.E. A .LE. R .LE. B.
C                          A MUST BE LESS THAN B AND A MUST BE
C                          BE NON-NEGATIVE.
C
C                        M
C                          THE NUMBER OF GRID POINTS IN THE INTERVAL
C                          (A,B).  THE GRID POINTS IN THE R-DIRECTION
C                          R-DIRECTION ARE GIVEN BY
C                          R(I) = A + (I-0.5)DR FOR I=1,2,...,M
C                          WHERE DR =(B-A)/M.
C                          M MUST BE GREATER THAN 2.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT R = A AND R = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT R = A
C                               (SEE NOTE BELOW) AND R = B.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT R = A
C                               (SEE NOTE BELOW) AND THE DERIVATIVE
C                               OF THE SOLUTION WITH RESPECT TO R IS
C                               SPECIFIED AT R = B.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND R = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE BELOW) AND THE
C                               SOLUTION IS SPECIFIED AT R = B.
C
C                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE SOLUTION IS
C                               SPECIFIED AT R = B.
C
C                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
C                               R = A = 0 AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO R IS SPECIFIED
C                               AT R = B.
C
C                          NOTE:
C                          IF A = 0, DO NOT USE MBDCND = 1,2,3, OR 4,
C                          BUT INSTEAD USE MBDCND = 5 OR 6.
C                          THE RESULTING APPROXIMATION GIVES THE ONLY
C                          MEANINGFUL BOUNDARY CONDITION,
C                          I.E. DU/DR = 0.
C                          (SEE D. GREENSPAN, 'INTRODUCTORY NUMERICAL
C                          ANALYSIS OF ELLIPTIC BOUNDARY VALUE
C                          PROBLEMS,' HARPER AND ROW, 1965, CHAPTER 5.)
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES (IF ANY)
C                          OF THE SOLUTION AT R = A.
C
C                          WHEN MBDCND = 1 OR 2,
C                            BDA(J) = U(A,Z(J)) ,       J=1,2,...,N.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DR)U(A,Z(J)) ,   J=1,2,...,N.
C
C                          WHEN MBDCND = 5 OR 6, BDA IS A DUMMY
C                          VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT R = B.
C
C                          WHEN MBDCND = 1,4,OR 5,
C                            BDB(J) = U(B,Z(J)) ,        J=1,2,...,N.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DR)U(B,Z(J)) ,   J=1,2,...,N.
C
C                        C,D
C                          THE RANGE OF Z, I.E. C .LE. Z .LE. D.
C                          C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE Z-DIRECTION
C                          ARE GIVEN BY Z(J) = C + (J-0.5)DZ,
C                          J=1,2,...,N, WHERE DZ = (D-C)/N.
C                          N MUST BE GREATER THAN 2.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT Z = C  AND Z = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.
C                               U(I,J) = U(I,N+J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT Z = C
C                               AND Z = D.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT Z = C
C                               AND THE DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO Z IS SPECIFIED AT Z = D.
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO Z IS SPECIFIED AT Z = C
C                               AND Z = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION WITH
C                               RESPECT TO Z IS SPECIFIED AT Z = C AND
C                               THE SOLUTION IS SPECIFIED AT Z = D.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Z = C.
C
C                          WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(R(I),C) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DZ)U(R(I),C),    I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT Z = D.
C
C                          WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(R(I),D) ,       I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) = (D/DZ)U(R(I),D) ,   I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE MODIFIED
C                          HELMHOLTZ EQUATION.  IF LAMBDA IS GREATER
C                          THAN 0, A SOLUTION MAY NOT EXIST.
C                          HOWEVER, HSTCYL WILL ATTEMPT TO FIND A
C                          SOLUTION.  LAMBDA MUST BE ZERO WHEN
C                          MBDCND = 5 OR 6.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
C                          THE VALUES OF THE RIGHT SIDE OF THE
C                          MODIFIED HELMHOLTZ EQUATION.
C                          FOR I=1,2,...,M   AND J=1,2,...,N
C                            F(I,J) = F(R(I),Z(J)) .
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTCYL.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
C                          BE AT LEAST M.
C
C                        W
C                          A ONE-DIMENSIONAL ARRAY THAT MUST BE
C                          PROVIDED BY THE USER FOR WORK SPACE. W MAY
C                          REQUIRE UP TO 13M + 4N + M*INT(LOG2(N))
C                          LOCATIONS.  THE ACTUAL NUMBER OF LOCATIONS
C                          USED IS COMPUTED BY HSTCYL AND IS RETURNED
C                          IN THE LOCATION W(1).
C
C ON OUTPUT
C
C                        F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (R(I),Z(J)) FOR  I=1,2,...,M, J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT, CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HSTCYL THEN COMPUTES
C                          THIS SOLUTION, WHICH IS A LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL APPROXIMATION.
C                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
C                          A SOLUTION; HENCE, THE SOLUTION IS NOT
C                          UNIQUE.  THE VALUE OF PERTRB SHOULD BE
C                          SMALL COMPARED TO THE RIGHT SIDE F.
C                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
C                          ESSENTIALLY DIFFERENT PROBLEM.
C                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
C                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
C                          OBTAINED.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS. EXCEPT TO NUMBERS 0 AND 11,
C                          A SOLUTION IS NOT ATTEMPTED.
C
C                          =  0  NO ERROR
C
C                          =  1  A .LT. 0
C
C                          =  2  A .GE. B
C
C                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6
C
C                          =  4  C .GE. D
C
C                          =  5  N .LE. 2
C
C                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4
C
C                          =  7  A = 0 AND MBDCND = 1,2,3, OR 4
C
C                          =  8  A .GT. 0 AND MBDCND .GE. 5
C
C                          =  9  M .LE. 2
C
C                          = 10  IDIMF .LT. M
C
C                          = 11  LAMBDA .GT. 0
C
C                          = 12  A=0, MBDCND .GE. 5, ELMBDA .NE. 0
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HSTCYL, THE
C                          USER SHOULD TEST IERROR AFTER THE CALL.
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
C                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR AND
C                        CALLS EITHER POISTG OR GENBUN WHICH SOLVES THE
C                        LINEAR SYSTEM OF EQUATIONS.
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               THE SOLUTION PROCESS RESULTS IN A LOSS
C                        OF NO MORE THAN FOUR SIGNIFICANT DIGITS
C                        FOR N AND M AS LARGE AS 64.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        SUBROUTINE POISTG WHICH IS THE ROUTINE THAT
C                        ACTUALLY SOLVES THE FINITE DIFFERENCE
C                        EQUATIONS.
C
C REFERENCES             U. SCHUMANN AND R. SWEET, "A DIRECT METHOD FOR
C                        THE SOLUTION OF POISSON'S EQUATION WITH NEUMANN
C                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
C                        ARBITRARY SIZE," J. COMP. PHYS. 20(1976),
C                        PP. 171-182.
C***********************************************************************
      DIMENSION       F(IDIMF,1) ,BDA(*)     ,BDB(*)     ,BDC(*)     ,
     1                BDD(*)     ,W(*)
C
      IERROR = 0
      IF (A .LT. 0.) IERROR = 1
      IF (A .GE. B) IERROR = 2
      IF (MBDCND.LE.0 .OR. MBDCND.GE.7) IERROR = 3
      IF (C .GE. D) IERROR = 4
      IF (N .LE. 2) IERROR = 5
      IF (NBDCND.LT.0 .OR. NBDCND.GE.5) IERROR = 6
      IF (A.EQ.0. .AND. MBDCND.NE.5 .AND. MBDCND.NE.6) IERROR = 7
      IF (A.GT.0. .AND. MBDCND.GE.5) IERROR = 8
      IF (IDIMF .LT. M) IERROR = 10
      IF (M .LE. 2) IERROR = 9
      IF (A.EQ.0. .AND. MBDCND.GE.5 .AND. ELMBDA.NE.0.) IERROR = 12
      IF (IERROR .NE. 0) RETURN
      DELTAR = (B-A)/FLOAT(M)
      DLRSQ = DELTAR**2
      DELTHT = (D-C)/FLOAT(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND+1
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      IWB = M
      IWC = IWB+M
      IWR = IWC+M
      DO 101 I=1,M
         J = IWR+I
         W(J) = A+(FLOAT(I)-0.5)*DELTAR
         W(I) = (A+FLOAT(I-1)*DELTAR)/(DLRSQ*W(J))
         K = IWC+I
         W(K) = (A+FLOAT(I)*DELTAR)/(DLRSQ*W(J))
         K = IWB+I
         W(K) = ELMBDA/W(J)**2-2./DLRSQ
  101 CONTINUE
C
C     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
C
      GO TO (102,102,104,104,106,106),MBDCND
  102 A1 = 2.*W(1)
      W(IWB+1) = W(IWB+1)-W(1)
      DO 103 J=1,N
         F(1,J) = F(1,J)-A1*BDA(J)
  103 CONTINUE
      GO TO 106
  104 A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1)+W(1)
      DO 105 J=1,N
         F(1,J) = F(1,J)+A1*BDA(J)
  105 CONTINUE
  106 CONTINUE
      GO TO (107,109,109,107,107,109),MBDCND
  107 W(IWC) = W(IWC)-W(IWR)
      A1 = 2.*W(IWR)
      DO 108 J=1,N
         F(M,J) = F(M,J)-A1*BDB(J)
  108 CONTINUE
      GO TO 111
  109 W(IWC) = W(IWC)+W(IWR)
      A1 = DELTAR*W(IWR)
      DO 110 J=1,N
         F(M,J) = F(M,J)-A1*BDB(J)
  110 CONTINUE
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
  111 A1 = 2./DLTHSQ
      GO TO (121,112,112,114,114),NP
  112 DO 113 I=1,M
         F(I,1) = F(I,1)-A1*BDC(I)
  113 CONTINUE
      GO TO 116
  114 A1 = 1./DELTHT
      DO 115 I=1,M
         F(I,1) = F(I,1)+A1*BDC(I)
  115 CONTINUE
  116 A1 = 2./DLTHSQ
      GO TO (121,117,119,119,117),NP
  117 DO 118 I=1,M
         F(I,N) = F(I,N)-A1*BDD(I)
  118 CONTINUE
      GO TO 121
  119 A1 = 1./DELTHT
      DO 120 I=1,M
         F(I,N) = F(I,N)-A1*BDD(I)
  120 CONTINUE
  121 CONTINUE
C
C     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
C     SOLUTION.
C
      PERTRB = 0.
      IF (ELMBDA) 130,123,122
  122 IERROR = 11
      GO TO 130
  123 GO TO (130,130,124,130,130,124),MBDCND
  124 GO TO (125,130,130,125,130),NP
  125 CONTINUE
      DO 127 I=1,M
         A1 = 0.
         DO 126 J=1,N
            A1 = A1+F(I,J)
  126    CONTINUE
         J = IWR+I
         PERTRB = PERTRB+A1*W(J)
  127 CONTINUE
      PERTRB = PERTRB/(FLOAT(M*N)*0.5*(A+B))
      DO 129 I=1,M
         DO 128 J=1,N
            F(I,J) = F(I,J)-PERTRB
  128    CONTINUE
  129 CONTINUE
  130 CONTINUE
C
C     MULTIPLY I-TH EQUATION THROUGH BY  DELTHT**2
C
      DO 132 I=1,M
         W(I) = W(I)*DLTHSQ
         J = IWC+I
         W(J) = W(J)*DLTHSQ
         J = IWB+I
         W(J) = W(J)*DLTHSQ
         DO 131 J=1,N
            F(I,J) = F(I,J)*DLTHSQ
  131    CONTINUE
  132 CONTINUE
      LP = NBDCND
      W(1) = 0.
      W(IWR) = 0.
C
C     CALL GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
C
      IF (NBDCND .EQ. 0) GO TO 133
      CALL POISTG (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
      GO TO 134
  133 CALL GENBUN (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
  134 CONTINUE
      W(1) = W(IWR+1)+3.*FLOAT(M)
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
