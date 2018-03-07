C
C     file hstplr.f
C
      SUBROUTINE HSTPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
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
C                        GRID TO THE HELMHOLTZ EQUATION IN POLAR
C                        COORDINATES.  THE EQUATION IS
C
C                           (1/R)(D/DR)(R(DU/DR)) +
C                           (1/R**2)(D/DTHETA)(DU/DTHETA) +
C                           LAMBDA*U = F(R,THETA)
C
C USAGE                  CALL HSTPLR (A,B,M,MBDCND,BDA,BDB,C,D,N,
C                                     NBDCND,BDC,BDD,ELMBDA,F,
C                                     IDIMF,PERTRB,IERROR,W)
C
C ARGUMENTS
C ON INPUT               A,B
C
C                          THE RANGE OF R, I.E. A .LE. R .LE. B.
C                          A MUST BE LESS THAN B AND A MUST BE
C                          NON-NEGATIVE.
C
C                        M
C                          THE NUMBER OF GRID POINTS IN THE INTERVAL
C                          (A,B).  THE GRID POINTS IN THE R-DIRECTION
C                          ARE GIVEN BY R(I) = A + (I-0.5)DR FOR
C                          I=1,2,...,M WHERE DR =(B-A)/M.
C                          M MUST BE GREATER THAN 2.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT R = A AND R = B.
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT R = A
C                               AND R = B.
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT R = A
C                               AND THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT R = B.
C                               (SEE NOTE 1 BELOW)
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               R = A (SEE NOTE 2 BELOW) AND R = B.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO R IS SPECIFIED AT
C                               SPECIFIED AT R = A (SEE NOTE 2 BELOW)
C                               AND THE SOLUTION IS SPECIFIED AT R = B.
C
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
C                          NOTE 1:
C                          IF A = 0, MBDCND = 2, AND NBDCND = 0 OR 3,
C                          THE SYSTEM OF EQUATIONS TO BE SOLVED IS
C                          SINGULAR.  THE UNIQUE SOLUTION IS
C                          IS DETERMINED BY EXTRAPOLATION TO THE
C                          SPECIFICATION OF U(0,THETA(1)).
C                          BUT IN THIS CASE THE RIGHT SIDE OF THE
C                          SYSTEM WILL BE PERTURBED BY THE CONSTANT
C                          PERTRB.
C
C                          NOTE 2:
C                          IF A = 0, DO NOT USE MBDCND = 3 OR 4,
C                          BUT INSTEAD USE MBDCND = 1,2,5, OR 6.
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES (IF ANY) OF
C                          THE SOLUTION AT R = A.
C
C                          WHEN MBDCND = 1 OR 2,
C                            BDA(J) = U(A,THETA(J)) ,     J=1,2,...,N.
C
C                          WHEN MBDCND = 3 OR 4,
C                            BDA(J) = (D/DR)U(A,THETA(J)) ,
C                            J=1,2,...,N.
C
C                          WHEN MBDCND = 5 OR 6, BDA IS A DUMMY
C                          VARIABLE.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT R = B.
C
C                          WHEN MBDCND = 1,4, OR 5,
C                            BDB(J) = U(B,THETA(J)) ,     J=1,2,...,N.
C
C                          WHEN MBDCND = 2,3, OR 6,
C                            BDB(J) = (D/DR)U(B,THETA(J)) ,
C                            J=1,2,...,N.
C
C                        C,D
C                          THE RANGE OF THETA, I.E. C .LE. THETA .LE. D.
C                          C MUST BE LESS THAN D.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
C                          (C,D).  THE UNKNOWNS IN THE THETA-
C                          DIRECTION ARE GIVEN BY THETA(J) = C +
C                          (J-0.5)DT,   J=1,2,...,N, WHERE
C                          DT = (D-C)/N.  N MUST BE GREATER THAN 2.
C
C                        NBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
C                          AT THETA = C  AND THETA = D.
C
C                          = 0  IF THE SOLUTION IS PERIODIC IN THETA,
C                               I.E. U(I,J) = U(I,N+J).
C
C                          = 1  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = C AND THETA = D
C                               (SEE NOTE BELOW).
C
C                          = 2  IF THE SOLUTION IS SPECIFIED AT
C                               THETA = C AND THE DERIVATIVE OF THE
C                               SOLUTION WITH RESPECT TO THETA IS
C                               SPECIFIED AT THETA = D
C                               (SEE NOTE BELOW).
C
C                          = 3  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = C AND THETA = D.
C
C                          = 4  IF THE DERIVATIVE OF THE SOLUTION
C                               WITH RESPECT TO THETA IS SPECIFIED
C                               AT THETA = C AND THE SOLUTION IS
C                               SPECIFIED AT THETA = D
C                               (SEE NOTE BELOW).
C
C                          NOTE:
C                          WHEN NBDCND = 1, 2, OR 4, DO NOT USE
C                          MBDCND = 5 OR 6 (THE FORMER INDICATES THAT
C                          THE SOLUTION IS SPECIFIED AT R =  0; THE
C                          LATTER INDICATES THE SOLUTION IS UNSPECIFIED
C                          AT R = 0).  USE INSTEAD MBDCND = 1 OR 2.
C
C                        BDC
C                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT THETA = C.
C
C                          WHEN NBDCND = 1 OR 2,
C                            BDC(I) = U(R(I),C) ,        I=1,2,...,M.
C
C                          WHEN NBDCND = 3 OR 4,
C                            BDC(I) = (D/DTHETA)U(R(I),C),
C                            I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
C                          SPECIFIES THE BOUNDARY VALUES OF THE
C                          SOLUTION AT THETA = D.
C
C                          WHEN NBDCND = 1 OR 4,
C                            BDD(I) = U(R(I),D) ,         I=1,2,...,M.
C
C                          WHEN NBDCND = 2 OR 3,
C                            BDD(I) =(D/DTHETA)U(R(I),D), I=1,2,...,M.
C
C                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
C
C                        ELMBDA
C                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
C                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
C                          A SOLUTION MAY NOT EXIST.  HOWEVER, HSTPLR
C                          WILL ATTEMPT TO FIND A SOLUTION.
C
C                        F
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT SIDE OF THE HELMHOLTZ
C                          EQUATION.
C
C                          FOR I=1,2,...,M AND J=1,2,...,N
C                            F(I,J) = F(R(I),THETA(J)) .
C
C                          F MUST BE DIMENSIONED AT LEAST M X N.
C
C                        IDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
C                          F AS IT APPEARS IN THE PROGRAM CALLING
C                          HSTPLR.  THIS PARAMETER IS USED TO SPECIFY
C                          THE VARIABLE DIMENSION OF F.
C                          IDIMF MUST BE AT LEAST M.
C
C                        W
C                          A ONE-DIMENSIONAL ARRAY THAT MUST BE
C                          PROVIDED BY THE USER FOR WORK SPACE.
C                          W MAY REQUIRE UP TO 13M + 4N +
C                          M*INT(LOG2(N)) LOCATIONS.
C                          THE ACTUAL NUMBER OF LOCATIONS USED IS
C                          COMPUTED BY HSTPLR AND IS RETURNED IN
C                          THE LOCATION W(1).
C
C
C ON OUTPUT
C
C                        F
C                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
C                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
C                          (R(I),THETA(J)) FOR I=1,2,...,M,
C                          J=1,2,...,N.
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
C                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
C                          SPECIFIED FOR A POISSON EQUATION
C                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
C                          PERTRB IS A CONSTANT CALCULATED AND
C                          SUBTRACTED FROM F, WHICH ENSURES THAT A
C                          SOLUTION EXISTS.  HSTPLR THEN COMPUTES THIS
C                          SOLUTION, WHICH IS A LEAST SQUARES SOLUTION
C                          TO THE ORIGINAL APPROXIMATION.
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
C                          =  7  A = 0 AND MBDCND = 3 OR 4
C
C                          =  8  A .GT. 0 AND MBDCND .GE. 5
C
C                          =  9  MBDCND .GE. 5 AND NBDCND .NE. 0 OR 3
C
C                          = 10  IDIMF .LT. M
C
C                          = 11  LAMBDA .GT. 0
C
C                          = 12  M .LE. 2
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING
C                          A POSSIBLY INCORRECT CALL TO HSTPLR, THE
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
C PORTABILITY            FORTRAN 77.
C
C ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-
C                        DIFFERENCE EQUATIONS, INCORPORATES BOUNDARY
C                        DATA, ADJUSTS THE RIGHT SIDE WHEN THE SYSTEM
C                        IS SINGULAR AND CALLS EITHER POISTG OR GENBUN
C                        WHICH SOLVES THE LINEAR SYSTEM OF EQUATIONS.
C
C TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
C                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
C                        DIGITS FOR N AND M AS LARGE AS 64.
C                        MORE DETAILED INFORMATION ABOUT ACCURACY
C                        CAN BE FOUND IN THE DOCUMENTATION FOR
C                        ROUTINE POISTG WHICH IS THE ROUTINE THAT
C                        ACTUALLY SOLVES THE FINITE DIFFERENCE
C                        EQUATIONS.
C
C REFERENCES             U. SCHUMANN AND R. SWEET, "A DIRECT METHOD
C                        FOR THE SOLUTION OF POISSON'S EQUATION WITH
C                        NEUMANN BOUNDARY CONDITIONS ON A STAGGERED
C                        GRID OF ARBITRARY SIZE," J. COMP. PHYS.
C                        20(1976), PP. 171-182.
C***********************************************************************
      DIMENSION       F(IDIMF,1)
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
C
      IERROR = 0
      IF (A .LT. 0.) IERROR = 1
      IF (A .GE. B) IERROR = 2
      IF (MBDCND.LE.0 .OR. MBDCND.GE.7) IERROR = 3
      IF (C .GE. D) IERROR = 4
      IF (N .LE. 2) IERROR = 5
      IF (NBDCND.LT.0 .OR. NBDCND.GE.5) IERROR = 6
      IF (A.EQ.0. .AND. (MBDCND.EQ.3 .OR. MBDCND.EQ.4)) IERROR = 7
      IF (A.GT.0. .AND. MBDCND.GE.5) IERROR = 8
      IF (MBDCND.GE.5 .AND. NBDCND.NE.0 .AND. NBDCND.NE.3) IERROR = 9
      IF (IDIMF .LT. M) IERROR = 10
      IF (M .LE. 2) IERROR = 12
      IF (IERROR .NE. 0) RETURN
      DELTAR = (B-A)/FLOAT(M)
      DLRSQ = DELTAR**2
      DELTHT = (D-C)/FLOAT(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND+1
      ISW = 1
      MB = MBDCND
      IF (A.EQ.0. .AND. MBDCND.EQ.2) MB = 6
C
C     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
C
      IWB = M
      IWC = IWB+M
      IWR = IWC+M
      DO 101 I=1,M
         J = IWR+I
         W(J) = A+(FLOAT(I)-0.5)*DELTAR
         W(I) = (A+FLOAT(I-1)*DELTAR)/DLRSQ
         K = IWC+I
         W(K) = (A+FLOAT(I)*DELTAR)/DLRSQ
         K = IWB+I
         W(K) = (ELMBDA-2./DLRSQ)*W(J)
  101 CONTINUE
      DO 103 I=1,M
         J = IWR+I
         A1 = W(J)
         DO 102 J=1,N
            F(I,J) = A1*F(I,J)
  102    CONTINUE
  103 CONTINUE
C
C     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
C
      GO TO (104,104,106,106,108,108),MB
  104 A1 = 2.*W(1)
      W(IWB+1) = W(IWB+1)-W(1)
      DO 105 J=1,N
         F(1,J) = F(1,J)-A1*BDA(J)
  105 CONTINUE
      GO TO 108
  106 A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1)+W(1)
      DO 107 J=1,N
         F(1,J) = F(1,J)+A1*BDA(J)
  107 CONTINUE
  108 GO TO (109,111,111,109,109,111),MB
  109 A1 = 2.*W(IWR)
      W(IWC) = W(IWC)-W(IWR)
      DO 110 J=1,N
         F(M,J) = F(M,J)-A1*BDB(J)
  110 CONTINUE
      GO TO 113
  111 A1 = DELTAR*W(IWR)
      W(IWC) = W(IWC)+W(IWR)
      DO 112 J=1,N
         F(M,J) = F(M,J)-A1*BDB(J)
  112 CONTINUE
C
C     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
C
  113 A1 = 2./DLTHSQ
      GO TO (123,114,114,116,116),NP
  114 DO 115 I=1,M
         J = IWR+I
         F(I,1) = F(I,1)-A1*BDC(I)/W(J)
  115 CONTINUE
      GO TO 118
  116 A1 = 1./DELTHT
      DO 117 I=1,M
         J = IWR+I
         F(I,1) = F(I,1)+A1*BDC(I)/W(J)
  117 CONTINUE
  118 A1 = 2./DLTHSQ
      GO TO (123,119,121,121,119),NP
  119 DO 120 I=1,M
         J = IWR+I
         F(I,N) = F(I,N)-A1*BDD(I)/W(J)
  120 CONTINUE
      GO TO 123
  121 A1 = 1./DELTHT
      DO 122 I=1,M
         J = IWR+I
         F(I,N) = F(I,N)-A1*BDD(I)/W(J)
  122 CONTINUE
  123 CONTINUE
C
C     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
C     SOLUTION.
C
      PERTRB = 0.
      IF (ELMBDA) 133,125,124
  124 IERROR = 11
      GO TO 133
  125 GO TO (133,133,126,133,133,126),MB
  126 GO TO (127,133,133,127,133),NP
  127 CONTINUE
      ISW = 2
      DO 129 J=1,N
         DO 128 I=1,M
            PERTRB = PERTRB+F(I,J)
  128    CONTINUE
  129 CONTINUE
      PERTRB = PERTRB/(FLOAT(M*N)*0.5*(A+B))
      DO 131 I=1,M
         J = IWR+I
         A1 = PERTRB*W(J)
         DO 130 J=1,N
            F(I,J) = F(I,J)-A1
  130    CONTINUE
  131 CONTINUE
      A2 = 0.
      DO 132 J=1,N
         A2 = A2+F(1,J)
  132 CONTINUE
      A2 = A2/W(IWR+1)
  133 CONTINUE
C
C     MULTIPLY I-TH EQUATION THROUGH BY  R(I)*DELTHT**2
C
      DO 135 I=1,M
         J = IWR+I
         A1 = DLTHSQ*W(J)
         W(I) = A1*W(I)
         J = IWC+I
         W(J) = A1*W(J)
         J = IWB+I
         W(J) = A1*W(J)
         DO 134 J=1,N
            F(I,J) = A1*F(I,J)
  134    CONTINUE
  135 CONTINUE
      LP = NBDCND
      W(1) = 0.
      W(IWR) = 0.
C
C     CALL POISTG OR GENBUN TO SOLVE THE SYSTEM OF EQUATIONS.
C
      IF (LP .EQ. 0) GO TO 136
      CALL POISTG (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
      GO TO 137
  136 CALL GENBUN (LP,N,1,M,W,W(IWB+1),W(IWC+1),IDIMF,F,IERR1,W(IWR+1))
  137 CONTINUE
      W(1) = W(IWR+1)+3.*FLOAT(M)
      IF (A.NE.0. .OR. MBDCND.NE.2 .OR. ISW.NE.2) GO TO 141
      A1 = 0.
      DO 138 J=1,N
         A1 = A1+F(1,J)
  138 CONTINUE
      A1 = (A1-DLRSQ*A2/16.)/FLOAT(N)
      IF (NBDCND .EQ. 3) A1 = A1+(BDD(1)-BDC(1))/(D-C)
      A1 = BDA(1)-A1
      DO 140 I=1,M
         DO 139 J=1,N
            F(I,J) = F(I,J)+A1
  139    CONTINUE
  140 CONTINUE
  141 CONTINUE
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
