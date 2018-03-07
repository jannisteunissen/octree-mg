C
C     file sepx4.f
C
      SUBROUTINE SEPX4 (IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     1                  NBDCND,BDC,BDD,COFX,GRHS,USOL,IDMN,W,PERTRB,
     2                  IERROR)
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
C DIMENSION OF           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
C ARGUMENTS              USOL(IDMN,N+1),     GRHS(IDMN,N+1),
C                        W (SEE ARGUMENT LIST)
C
C LATEST REVISION        NOVEMBER 1988
C
C PURPOSE                SEPX4 SOLVES FOR EITHER THE SECOND-ORDER
C                        FINITE DIFFERENCE APPROXIMATION OR A
C                        FOURTH-ORDER APPROXIMATION TO A SEPARABLE
C                        ELLIPTIC EQUATION
C
C                          AF(X)*UXX+BF(X)*UX+CF(X)*U+UYY = G(X,Y)
C
C                        ON A RECTANGLE (X GREATER THAN OR EQUAL TO
C                        A AND LESS THAN OR EQUAL TO B, Y GREATER THAN
C                        OR EQUAL TO C AND LESS THAN OR EQUAL TO D).
C                        ANY COMBINATION OF PERIODIC OR MIXED BOUNDARY
C                        CONDITIONS IS ALLOWED.  IF BOUNDARY
C                        CONDITIONS IN THE X DIRECTION ARE PERIODIC
C                        (SEE MBDCND=0 BELOW) THEN THE COEFFICIENTS
C                        MUST SATISFY
C
C                          AF(X)=C1,BF(X)=0,CF(X)=C2 FOR ALL X.
C
C                        HERE C1,C2 ARE CONSTANTS, C1.GT.0.
C
C                        THE POSSIBLE BOUNDARY CONDITIONS ARE:
C                        IN THE X-DIRECTION:
C                          (0) PERIODIC, U(X+B-A,Y)=U(X,Y) FOR
C                              ALL Y,X
C                          (1) U(A,Y), U(B,Y) ARE SPECIFIED FOR ALL Y
C                          (2) U(A,Y), DU(B,Y)/DX+BETA*U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C                          (3) DU(A,Y)/DX+ALPHA*U(A,Y),DU(B,Y)/DX+
C                              BETA*U(B,Y) ARE SPECIFIED FOR ALL Y
C                          (4) DU(A,Y)/DX+ALPHA*U(A,Y),U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C
C                        IN THE Y-DIRECTION:
C                          (0) PERIODIC, U(X,Y+D-C)=U(X,Y) FOR ALL X,Y
C                          (1) U(X,C),U(X,D) ARE SPECIFIED FOR ALL X
C                          (2) U(X,C),DU(X,D)/DY ARE SPECIFIED FOR
C                              ALL X
C                          (3) DU(X,C)/DY,DU(X,D)/DY ARE SPECIFIED FOR
C                              ALL X
C                          (4) DU(X,C)/DY,U(X,D) ARE SPECIFIED FOR
C                              ALL X
C
C USAGE                  CALL SEPX4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,
C                                   BETA,C,D,N,NBDCND,BDC,BDD,COFX,
C                                   GRHS,USOL,IDMN,W,PERTRB,IERROR)
C
C ARGUMENTS
C ON INPUT               IORDER
C                          = 2 IF A SECOND-ORDER APPROXIMATION IS
C                              SOUGHT
C                          = 4 IF A FOURTH-ORDER APPROXIMATION IS
C                              SOUGHT
C
C                        A,B
C                          THE RANGE OF THE X-INDEPENDENT VARIABLE,
C                          I.E., X IS GREATER THAN OR EQUAL TO A
C                          AND LESS THAN OR EQUAL TO B.  A MUST BE
C                          LESS THAN B.
C
C                        M
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (A,B) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE M+1 GRID POINTS IN THE X-
C                          DIRECTION GIVEN BY XI=A+(I-1)*DLX
C                          FOR I=1,2,...,M+1 WHERE DLX=(B-A)/M IS
C                          THE PANEL WIDTH.  M MUST BE LESS THAN
C                          IDMN AND GREATER THAN 5.
C
C                        MBDCND
C                          INDICATES THE TYPE OF BOUNDARY CONDITION
C                          AT X=A AND X=B
C                          = 0 IF THE SOLUTION IS PERIODIC IN X, I.E.,
C                              U(X+B-A,Y)=U(X,Y)  FOR ALL Y,X
C                          = 1 IF THE SOLUTION IS SPECIFIED AT X=A
C                              AND X=B, I.E., U(A,Y) AND U(B,Y) ARE
C                              SPECIFIED FOR ALL Y
C                          = 2 IF THE SOLUTION IS SPECIFIED AT X=A
C                              AND THE BOUNDARY CONDITION IS MIXED AT
C                              X=B, I.E., U(A,Y) AND
C                              DU(B,Y)/DX+BETA*U(B,Y) ARE SPECIFIED
C                              FOR ALL Y
C                          = 3 IF THE BOUNDARY CONDITIONS AT X=A AND
C                              X=B ARE MIXED, I.E.,
C                              DU(A,Y)/DX+ALPHA*U(A,Y) AND
C                              DU(B,Y)/DX+BETA*U(B,Y) ARE SPECIFIED
C                              FOR ALL Y
C                          = 4 IF THE BOUNDARY CONDITION AT X=A IS
C                              MIXED AND THE SOLUTION IS SPECIFIED
C                              AT X=B, I.E., DU(A,Y)/DX+ALPHA*U(A,Y)
C                              AND U(B,Y) ARE SPECIFIED FOR ALL Y
C
C                        BDA
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF
C                          DU(A,Y)/DX+ ALPHA*U(A,Y) AT X=A, WHEN
C                          MBDCND=3 OR 4.
C                          BDA(J) = DU(A,YJ)/DX+ALPHA*U(A,YJ),
C                          J=1,2,...,N+1
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
C                          A DUMMY PARAMETER.
C
C                        ALPHA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN CASE
C                          OF A MIXED BOUNDARY CONDITION AT X=A
C                          (SEE ARGUMENT BDA).  IF MBDCND IS NOT EQUAL
C                          TO EITHER 3 OR 4, THEN ALPHA IS A DUMMY
C                          PARAMETER.
C
C                        BDB
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
C                          SPECIFIES THE VALUES OF
C                          DU(B,Y)/DX+ BETA*U(B,Y) AT X=B.
C                          WHEN MBDCND=2 OR 3
C                          BDB(J) = DU(B,YJ)/DX+BETA*U(B,YJ),
C                          J=1,2,...,N+1
C                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
C                          A DUMMY PARAMETER.
C
C                        BETA
C                          THE SCALAR MULTIPLYING THE SOLUTION IN
C                          CASE OF A MIXED BOUNDARY CONDITION AT X=B
C                          (SEE ARGUMENT BDB).  IF MBDCND IS NOT EQUAL
C                          TO 2 OR 3, THEN BETA IS A DUMMY PARAMETER.
C
C                        C,D
C                          THE RANGE OF THE Y-INDEPENDENT VARIABLE,
C                          I.E., Y IS GREATER THAN OR EQUAL TO C AND
C                          LESS THAN OR EQUAL TO D.  C MUST BE LESS
C                          THAN D.
C
C                        N
C                          THE NUMBER OF PANELS INTO WHICH THE
C                          INTERVAL (C,D) IS SUBDIVIDED.  HENCE,
C                          THERE WILL BE N+1 GRID POINTS IN THE Y-
C                          DIRECTION GIVEN BY YJ=C+(J-1)*DLY FOR
C                          J=1,2,...,N+1 WHERE DLY=(D-C)/N IS THE
C                          PANEL WIDTH.  IN ADDITION, N MUST BE
C                          GREATER THAN 4.
C
C                        NBDCND
C                          INDICATES THE TYPES OF BOUNDARY CONDITIONS
C                          AT Y=C AND Y=D
C                          = 0 IF THE SOLUTION IS PERIODIC IN Y,
C                              I.E., U(X,Y+D-C)=U(X,Y) FOR ALL X,Y
C                          = 1 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND Y = D, I.E., U(X,C)  AND U(X,D)
C                              ARE SPECIFIED FOR ALL X
C                          = 2 IF THE SOLUTION IS SPECIFIED AT Y=C
C                              AND THE BOUNDARY CONDITION IS MIXED
C                              AT Y=D, I.E., DU(X,C)/DY AND U(X,D)
C                              ARE SPECIFIED FOR ALL X
C                          = 3 IF THE BOUNDARY CONDITIONS ARE MIXED
C                              AT Y=CAND Y=D I.E.,
C                              DU(X,D)/DY AND DU(X,D)/DY ARE
C                              SPECIFIED FOR ALL X
C                          = 4 IF THE BOUNDARY CONDITION IS MIXED
C                              AT Y=C AND THE SOLUTION IS SPECIFIED
C                              AT Y=D, I.E. DU(X,C)/DY+GAMA*U(X,C)
C                              AND U(X,D) ARE SPECIFIED FOR ALL X
C
C                        BDC
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIES THE VALUE DU(X,C)/DY AT Y=C.
C
C                          WHEN NBDCND=3 OR 4
C                            BDC(I) = DU(XI,C)/DY I=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
C                          A DUMMY PARAMETER.
C
C                        BDD
C                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
C                          SPECIFIED THE VALUE OF DU(X,D)/DY AT Y=D.
C
C                          WHEN NBDCND=2 OR 3
C                            BDD(I)=DU(XI,D)/DY I=1,2,...,M+1.
C
C                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
C                          A DUMMY PARAMETER.
C
C                        COFX
C                          A USER-SUPPLIED SUBPROGRAM WITH PARAMETERS
C                          X, AFUN, BFUN, CFUN WHICH RETURNS THE
C                          VALUES OF THE X-DEPENDENT COEFFICIENTS
C                          AF(X), BF(X), CF(X) IN THE ELLIPTIC
C                          EQUATION AT X.  IF BOUNDARY CONDITIONS IN
C                          THE X DIRECTION ARE PERIODIC THEN THE
C                          COEFFICIENTS MUST SATISFY AF(X)=C1,BF(X)=0,
C                          CF(X)=C2 FOR ALL X.  HERE C1.GT.0
C                          AND C2 ARE CONSTANTS.
C
C                          NOTE THAT COFX MUST BE DECLARED EXTERNAL
C                          IN THE CALLING ROUTINE.
C
C                        GRHS
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT-HAND SIDE OF THE
C                          ELLIPTIC EQUATION, I.E.,GRHS(I,J)=G(XI,YI),
C                          FOR I=2,...,M, J=2,...,N.  AT THE
C                          BOUNDARIES, GRHS IS DEFINED BY
C
C                          MBDCND   GRHS(1,J)   GRHS(M+1,J)
C                          ------   ---------   -----------
C                            0      G(A,YJ)     G(B,YJ)
C                            1         *           *
C                            2         *        G(B,YJ)  J=1,2,...,N+1
C                            3      G(A,YJ)     G(B,YJ)
C                            4      G(A,YJ)        *
C
C                          NBDCND   GRHS(I,1)   GRHS(I,N+1)
C                          ------   ---------   -----------
C                            0      G(XI,C)     G(XI,D)
C                            1         *           *
C                            2         *        G(XI,D)  I=1,2,...,M+1
C                            3      G(XI,C)     G(XI,D)
C                            4      G(XI,C)        *
C
C                          WHERE * MEANS THESE QUANTITES ARE NOT USED.
C                          GRHS SHOULD BE DIMENSIONED IDMN BY AT LEAST
C                          N+1 IN THE CALLING ROUTINE.
C
C                        USOL
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE SOLUTION ALONG THE BOUNDARIES.
C                          AT THE BOUNDARIES, USOL IS DEFINED BY
C
C                          MBDCND   USOL(1,J)   USOL(M+1,J)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(A,YJ)     U(B,YJ)
C                            2      U(A,YJ)        *     J=1,2,...,N+1
C                            3         *           *
C                            4         *        U(B,YJ)
C
C                          NBDCND   USOL(I,1)   USOL(I,N+1)
C                          ------   ---------   -----------
C                            0         *           *
C                            1      U(XI,C)     U(XI,D)
C                            2      U(XI,C)        *     I=1,2,...,M+1
C                            3         *           *
C                            4         *        U(XI,D)
C
C                          WHERE * MEANS THE QUANTITES ARE NOT USED
C                          IN THE SOLUTION.
C
C                          IF IORDER=2, THE USER MAY EQUIVALENCE GRHS
C                          AND USOL TO SAVE SPACE.  NOTE THAT IN THIS
C                          CASE THE TABLES SPECIFYING THE BOUNDARIES
C                          OF THE GRHS AND USOL ARRAYS DETERMINE THE
C                          BOUNDARIES UNIQUELY EXCEPT AT THE CORNERS.
C                          IF THE TABLES CALL FOR BOTH G(X,Y) AND
C                          U(X,Y) AT A CORNER THEN THE SOLUTION MUST
C                          BE CHOSEN.
C                          FOR EXAMPLE, IF MBDCND=2 AND NBDCND=4,
C                          THEN U(A,C), U(A,D),U(B,D) MUST BE CHOSEN
C                          AT THE CORNERS IN ADDITION TO G(B,C).
C
C                          IF IORDER=4, THEN THE TWO ARRAYS, USOL AND
C                          GRHS, MUST BE DISTINCT.
C
C                          USOL SHOULD BE DIMENSIONED IDMN BY AT LEAST
C                          N+1 IN THE CALLING ROUTINE.
C
C                        IDMN
C                          THE ROW (OR FIRST) DIMENSION OF THE ARRAYS
C                          GRHS AND USOL AS IT APPEARS IN THE PROGRAM
C                          CALLING SEPELI.  THIS PARAMETER IS USED
C                          TO SPECIFY THE VARIABLE DIMENSION OF GRHS
C                          AND USOL.  IDMN MUST BE AT LEAST 7 AND
C                          GREATER THAN OR EQUAL TO M+1.
C
C                        W
C                          A ONE-DIMENSIONAL ARRAY THAT MUST BE
C                          PROVIDED BY THE USER FOR WORK SPACE.
C                          10*N+(16+INT(LOG2(N+1)))*(M+1)+11 WILL
C                          SUFFICE AS A LENGTH FOR W.  THE ACTUAL
C                          LENGTH OF W IN THE CALLING ROUTINE
C                          MUST BE SET IN W(1) (SEE IERROR=11).
C
C
C ON OUTPUT              USOL
C                          CONTAINS THE APPROXIMATE SOLUTION TO THE
C                          ELLIPTIC EQUATION. USOL(I,J) IS THE
C                          APPROXIMATION TO U(XI,YJ) FOR I=1,2...,M+1
C                          AND J=1,2,...,N+1.  THE APPROXIMATION HAS
C                          ERROR O(DLX**2+DLY**2) IF CALLED WITH
C                          IORDER=2 AND O(DLX**4+DLY**4) IF CALLED
C                          WITH IORDER=4.
C
C                        W
C                          CONTAINS INTERMEDIATE VALUES THAT MUST NOT
C                          BE DESTROYED IF SEPELI IS CALLED AGAIN
C                          WITH INTL=1.  IN ADDITION W(1) CONTAINS
C                          THE EXACT MINIMAL LENGTH (IN FLOATING POINT)
C                          REQUIRED FOR THE WORK SPACE (SEE IERROR=11).
C
C                        PERTRB
C                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
C                          BOUNDARY CONDITIONS (I.E., ALPHA=BETA=0 IF
C                          MBDCND=3) IS SPECIFIED AND IF CF(X)=0 FOR
C                          ALL X THEN A SOLUTION TO THE DISCRETIZED
C                          MATRIX EQUATION MAY NOT EXIST
C                          (REFLECTING THE NON-UNIQUENESS OF SOLUTIONS
C                          TO THE PDE).
C                          PERTRB IS A CONSTANT CALCULATED AND
C                          SUBTRACTED FROM THE RIGHT HAND SIDE OF THE
C                          MATRIX EQUATION INSURING THE EXISTENCE OF A
C                          SOLUTION.  SEPX4 COMPUTES THIS SOLUTION
C                          WHICH IS A WEIGHTED MINIMAL LEAST SQUARES
C                          SOLUTION TO THE ORIGINAL PROBLEM.  IF
C                          SINGULARITY IS NOT DETECTED PERTRB=0.0 IS
C                          RETURNED BY SEPX4.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS OR FAILURE TO FIND A SOLUTION
C
C                          =  0 NO ERROR
C                          =  1 IF A GREATER THAN B OR C GREATER
C                               THAN D
C                          =  2 IF MBDCND LESS THAN 0 OR MBDCND
C                               GREATER THAN 4
C                          =  3 IF NBDCND LESS THAN 0 OR NBDCND
C                               GREATER THAN 4
C                          =  4 IF ATTEMPT TO FIND A SOLUTION FAILS.
C                               (THE LINEAR SYSTEM GENERATED IS NOT
C                               DIAGONALLY DOMINANT.)
C                          =  5 IF IDMN IS TOO SMALL (SEE DISCUSSION
C                               OF IDMN)
C                          =  6 IF M IS TOO SMALL OR TOO LARGE
C                               (SEE DISCUSSION OF M)
C                          =  7 IF N IS TOO SMALL (SEE DISCUSSION OF N)
C                          =  8 IF IORDER IS NOT 2 OR 4
C                          =  9 IF INTL IS NOT 0 OR 1
C                          = 10 IF AFUN IS LESS THAN OR EQUAL TO ZERO
C                               FOR SOME INTERIOR MESH POINT XI SOME
C                               INTERIOR MESH POINT (XI,YJ)
C                          = 11 IF THE WORK SPACE LENGTH INPUT IN W(1)
C                               IS LESS THAN THE EXACT MINIMAL WORK
C                               SPACE LENGTH REQUIRED OUTPUT IN W(1).
C                          = 12 IF MBDCND=0 AND AF(X)=CF(X)=CONSTANT
C                               OR BF(X)=0 FOR ALL X IS NOT TRUE.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C REQUIRED LIBRARY       COMF, GENBUN, GNBNAUX, AND SEPAUX
C FILES                  FROM FISHPACK
C
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       NONE
C FILES
C
C LANGUAGE               FORTRAN
C
C HISTORY                SEPX4 WAS DEVELOPED AT NCAR BY JOHN C.
C                        ADAMS OF THE SCIENTIFIC COMPUTING DIVISION
C                        IN OCTOBER 1978.  THE BASIS OF THIS CODE IS
C                        NCAR ROUTINE SEPELI.  BOTH PACKAGES WERE
C                        RELEASED ON NCAR'S PUBLIC LIBRARIES IN
C                        JANUARY 1980.
C
C PORTABILITY            FORTRAN 77
C
C ALGORITHM              SEPX4 AUTOMATICALLY DISCRETIZES THE SEPARABLE
C                        ELLIPTIC EQUATION WHICH IS THEN SOLVED BY A
C                        GENERALIZED CYCLIC REDUCTION ALGORITHM IN THE
C                        SUBROUTINE POIS.  THE FOURTH ORDER SOLUTION
C                        IS OBTAINED USING THE TECHNIQUE OF DEFFERRED
C                        CORRECTIONS REFERENCED BELOW.
C
C TIMING                 WHEN POSSIBLE, SEPX4 SHOULD BE USED INSTEAD
C                        OF PACKAGE SEPELI.  THE INCREASE IN SPEED
C                        IS AT LEAST A FACTOR OF THREE.
C
C REFERENCES             KELLER, H.B., NUMERICAL METHODS FOR TWO-POINT
C                        BOUNDARY-VALUE PROBLEMS, BLAISDEL (1968),
C                        WALTHAM, MASS.
C
C                        SWARZTRAUBER, P., AND R. SWEET (1975):
C                        EFFICIENT FORTRAN SUBPROGRAMS FOR THE
C                        SOLUTION OF ELLIPTIC PARTIAL DIFFERENTIAL
C                        EQUATIONS.  NCAR TECHNICAL NOTE
C                          NCAR-TN/IA-109, PP. 135-137.
C***********************************************************************
      DIMENSION       GRHS(IDMN,1)           ,USOL(IDMN,1)
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
C
C     CHECK INPUT PARAMETERS
C
      CALL C4KPRM(IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,IDMN,IERROR)
      IF (IERROR .NE. 0) RETURN
C
C     COMPUTE MINIMUM WORK SPACE AND CHECK WORK SPACE LENGTH INPUT
C
      L = N+1
      IF (NBDCND .EQ. 0) L = N
      K = M+1
      L = N+1
C     ESTIMATE LOG BASE 2 OF N
      LOG2N=INT(ALOG(FLOAT(N+1))/ALOG(2.0)+0.5)
      LENGTH=4*(N+1)+(10+LOG2N)*(M+1)
      IERROR = 11
      LINPUT = INT(W(1)+0.5)
      LOUTPT = LENGTH+6*(K+L)+1
      W(1) = FLOAT(LOUTPT)
      IF (LOUTPT .GT. LINPUT) RETURN
      IERROR = 0
C
C     SET WORK SPACE INDICES
C
      I1 = LENGTH+2
      I2 = I1+L
      I3 = I2+L
      I4 = I3+L
      I5 = I4+L
      I6 = I5+L
      I7 = I6+L
      I8 = I7+K
      I9 = I8+K
      I10 = I9+K
      I11 = I10+K
      I12 = I11+K
      I13 = 2
      CALL S4ELIP(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     1            NBDCND,BDC,BDD,COFX,W(I1),W(I2),W(I3),
     2            W(I4),W(I5),W(I6),W(I7),W(I8),W(I9),W(I10),W(I11),
     3            W(I12),GRHS,USOL,IDMN,W(I13),PERTRB,IERROR)
      RETURN
      END
      SUBROUTINE S4ELIP(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     1                  NBDCND,BDC,BDD,COFX,AN,BN,CN,DN,UN,ZN,AM,BM,
     2                  CM,DM,UM,ZM,GRHS,USOL,IDMN,W,PERTRB,IERROR)
C
C     S4ELIP SETS UP VECTORS AND ARRAYS FOR INPUT TO BLKTRI
C     AND COMPUTES A SECOND ORDER SOLUTION IN USOL.  A RETURN JUMP TO
C     SEPELI OCCURRS IF IORDER=2.  IF IORDER=4 A FOURTH ORDER
C     SOLUTION IS GENERATED IN USOL.
C
      DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     ,
     1                W(*)
      DIMENSION       GRHS(IDMN,1)           ,USOL(IDMN,1)
      DIMENSION       AN(*)      ,BN(*)      ,CN(*)      ,DN(*)      ,
     1                UN(*)      ,ZN(*)
      DIMENSION       AM(*)      ,BM(*)      ,CM(*)      ,DM(*)      ,
     1                UM(*)      ,ZM(*)
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      LOGICAL         SINGLR
      EXTERNAL        COFX
C
C     SET PARAMETERS INTERNALLY
C
      KSWX = MBDCND+1
      KSWY = NBDCND+1
      K = M+1
      L = N+1
      AIT = A
      BIT = B
      CIT = C
      DIT = D
      DLY=(DIT-CIT)/FLOAT(N)
C
C     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
C     AND NON-SPECIFIED BOUNDARIES.
C
      DO  20 I=2,M
         DO  10 J=2,N
      USOL(I,J)=DLY**2*GRHS(I,J)
   10    CONTINUE
   20 CONTINUE
      IF (KSWX.EQ.2 .OR. KSWX.EQ.3) GO TO  40
      DO  30 J=2,N
      USOL(1,J)=DLY**2*GRHS(1,J)
   30 CONTINUE
   40 CONTINUE
      IF (KSWX.EQ.2 .OR. KSWX.EQ.5) GO TO  60
      DO  50 J=2,N
      USOL(K,J)=DLY**2*GRHS(K,J)
   50 CONTINUE
   60 CONTINUE
      IF (KSWY.EQ.2 .OR. KSWY.EQ.3) GO TO  80
      DO  70 I=2,M
      USOL(I,1)=DLY**2*GRHS(I,1)
   70 CONTINUE
   80 CONTINUE
      IF (KSWY.EQ.2 .OR. KSWY.EQ.5) GO TO 100
      DO  90 I=2,M
      USOL(I,L)=DLY**2*GRHS(I,L)
   90 CONTINUE
  100 CONTINUE
      IF (KSWX.NE.2 .AND. KSWX.NE.3 .AND. KSWY.NE.2 .AND. KSWY.NE.3)
     1USOL(1,1)=DLY**2*GRHS(1,1)
      IF (KSWX.NE.2 .AND. KSWX.NE.5 .AND. KSWY.NE.2 .AND. KSWY.NE.3)
     1USOL(K,1)=DLY**2*GRHS(K,1)
      IF (KSWX.NE.2 .AND. KSWX.NE.3 .AND. KSWY.NE.2 .AND. KSWY.NE.5)
     1USOL(1,L)=DLY**2*GRHS(1,L)
      IF (KSWX.NE.2 .AND. KSWX.NE.5 .AND. KSWY.NE.2 .AND. KSWY.NE.5)
     1USOL(K,L)=DLY**2*GRHS(K,L)
      I1 = 1
C
C     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
C
      MP=1
      IF(KSWX.EQ.1) MP=0
      NP=NBDCND
C
C     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
C     IN NINT,MINT
C
      DLX = (BIT-AIT)/FLOAT(M)
      MIT = K-1
      IF (KSWX .EQ. 2) MIT = K-2
      IF (KSWX .EQ. 4) MIT = K
      DLY = (DIT-CIT)/FLOAT(N)
      NIT = L-1
      IF (KSWY .EQ. 2) NIT = L-2
      IF (KSWY .EQ. 4) NIT = L
      TDLX3 = 2.0*DLX**3
      DLX4 = DLX**4
      TDLY3 = 2.0*DLY**3
      DLY4 = DLY**4
C
C     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI
C
      IS = 1
      JS = 1
      IF (KSWX.EQ.2 .OR. KSWX.EQ.3) IS = 2
      IF (KSWY.EQ.2 .OR. KSWY.EQ.3) JS = 2
      NS = NIT+JS-1
      MS = MIT+IS-1
C
C     SET X - DIRECTION
C
      DO 110 I=1,MIT
         XI = AIT+FLOAT(IS+I-2)*DLX
         CALL COFX (XI,AI,BI,CI)
         AXI = (AI/DLX-0.5*BI)/DLX
         BXI = -2.*AI/DLX**2+CI
         CXI = (AI/DLX+0.5*BI)/DLX
      AM(I)=DLY**2*AXI
      BM(I)=DLY**2*BXI
      CM(I)=DLY**2*CXI
  110 CONTINUE
C
C     SET Y DIRECTION
C
      DYJ=1.0
      EYJ=-2.0
      FYJ=1.0
      DO 120 J=1,NIT
         AN(J) = DYJ
         BN(J) = EYJ
         CN(J) = FYJ
  120 CONTINUE
C
C     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
C
      AX1 = AM(1)
      CXM = CM(MIT)
      GO TO (170,130,150,160,140),KSWX
C
C     DIRICHLET-DIRICHLET IN X DIRECTION
C
  130 AM(1) = 0.0
      CM(MIT) = 0.0
      GO TO 170
C
C     MIXED-DIRICHLET IN X DIRECTION
C
  140 AM(1) = 0.0
      BM(1) = BM(1)+2.*ALPHA*DLX*AX1
      CM(1) = CM(1)+AX1
      CM(MIT) = 0.0
      GO TO 170
C
C     DIRICHLET-MIXED IN X DIRECTION
C
  150 AM(1) = 0.0
      AM(MIT) = AM(MIT)+CXM
      BM(MIT) = BM(MIT)-2.*BETA*DLX*CXM
      CM(MIT) = 0.0
      GO TO 170
C
C     MIXED - MIXED IN X DIRECTION
C
  160 CONTINUE
      AM(1) = 0.0
      BM(1) = BM(1)+2.*DLX*ALPHA*AX1
      CM(1) = CM(1)+AX1
      AM(MIT) = AM(MIT)+CXM
      BM(MIT) = BM(MIT)-2.*DLX*BETA*CXM
      CM(MIT) = 0.0
  170 CONTINUE
C
C     ADJUST IN Y DIRECTION UNLESS PERIODIC
C
      DY1 = AN(1)
      FYN = CN(NIT)
      GAMA=0.0
      XNU=0.0
      GO TO (220,180,200,210,190),KSWY
C
C     DIRICHLET-DIRICHLET IN Y DIRECTION
C
  180 CONTINUE
      AN(1) = 0.0
      CN(NIT) = 0.0
      GO TO 220
C
C     MIXED-DIRICHLET IN Y DIRECTION
C
  190 CONTINUE
      AN(1) = 0.0
      BN(1) = BN(1)+2.*DLY*GAMA*DY1
      CN(1) = CN(1)+DY1
      CN(NIT) = 0.0
      GO TO 220
C
C     DIRICHLET-MIXED IN Y DIRECTION
C
  200 AN(1) = 0.0
      AN(NIT) = AN(NIT)+FYN
      BN(NIT) = BN(NIT)-2.*DLY*XNU*FYN
      CN(NIT) = 0.0
      GO TO 220
C
C     MIXED - MIXED DIRECTION IN Y DIRECTION
C
  210 CONTINUE
      AN(1) = 0.0
      BN(1) = BN(1)+2.*DLY*GAMA*DY1
      CN(1) = CN(1)+DY1
      AN(NIT) = AN(NIT)+FYN
      BN(NIT) = BN(NIT)-2.0*DLY*XNU*FYN
      CN(NIT) = 0.0
  220 IF (KSWX .EQ. 1) GO TO 270
C
C     ADJUST USOL ALONG X EDGE
C
      DO 260 J=JS,NS
         IF (KSWX.NE.2 .AND. KSWX.NE.3) GO TO 230
         USOL(IS,J) = USOL(IS,J)-AX1*USOL(1,J)
         GO TO 240
  230    USOL(IS,J) = USOL(IS,J)+2.0*DLX*AX1*BDA(J)
  240    IF (KSWX.NE.2 .AND. KSWX.NE.5) GO TO 250
         USOL(MS,J) = USOL(MS,J)-CXM*USOL(K,J)
         GO TO 260
  250    USOL(MS,J) = USOL(MS,J)-2.0*DLX*CXM*BDB(J)
  260 CONTINUE
  270 IF (KSWY .EQ. 1) GO TO 320
C
C     ADJUST USOL ALONG Y EDGE
C
      DO 310 I=IS,MS
         IF (KSWY.NE.2 .AND. KSWY.NE.3) GO TO 280
         USOL(I,JS) = USOL(I,JS)-DY1*USOL(I,1)
         GO TO 290
  280    USOL(I,JS) = USOL(I,JS)+2.0*DLY*DY1*BDC(I)
  290    IF (KSWY.NE.2 .AND. KSWY.NE.5) GO TO 300
         USOL(I,NS) = USOL(I,NS)-FYN*USOL(I,L)
         GO TO 310
  300    USOL(I,NS) = USOL(I,NS)-2.0*DLY*FYN*BDD(I)
  310 CONTINUE
  320 CONTINUE
C
C     SAVE ADJUSTED EDGES IN GRHS IF IORDER=4
C
      IF (IORDER .NE. 4) GO TO 350
      DO 330 J=JS,NS
         GRHS(IS,J) = USOL(IS,J)
         GRHS(MS,J) = USOL(MS,J)
  330 CONTINUE
      DO 340 I=IS,MS
         GRHS(I,JS) = USOL(I,JS)
         GRHS(I,NS) = USOL(I,NS)
  340 CONTINUE
  350 CONTINUE
      IORD = IORDER
      PERTRB = 0.0
C
C     CHECK IF OPERATOR IS SINGULAR
C
      CALL C4KSNG(MBDCND,NBDCND,ALPHA,BETA,COFX,SINGLR)
C
C     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
C     IF SINGULAR
C
      IF (SINGLR) CALL SEPTRI (MIT,AM,BM,CM,DM,UM,ZM)
      IF (SINGLR) CALL SEPTRI (NIT,AN,BN,CN,DN,UN,ZN)
C
C     ADJUST RIGHT HAND SIDE IF NECESSARY
C
  360 CONTINUE
      IF (SINGLR) CALL SEPORT (USOL,IDMN,ZN,ZM,PERTRB)
C
C     COMPUTE SOLUTION
C
C     SAVE ADJUSTED RIGHT HAND SIDE IN GRHS
      DO 444 J=JS,NS
      DO 444 I=IS,MS
      GRHS(I,J)=USOL(I,J)
  444 CONTINUE
      CALL GENBUN(NP,NIT,MP,MIT,AM,BM,CM,IDMN,USOL(IS,JS),IEROR,W)
C     CHECK IF ERROR DETECTED IN POIS
C     THIS CAN ONLY CORRESPOND TO IERROR=12
      IF(IEROR.EQ.0) GO TO 224
C     SET ERROR FLAG IF IMPROPER COEFFICIENTS INPUT TO POIS
      IERROR=12
      RETURN
  224 CONTINUE
      IF (IERROR .NE. 0) RETURN
C
C     SET PERIODIC BOUNDARIES IF NECESSARY
C
      IF (KSWX .NE. 1) GO TO 380
      DO 370 J=1,L
         USOL(K,J) = USOL(1,J)
  370 CONTINUE
  380 IF (KSWY .NE. 1) GO TO 400
      DO 390 I=1,K
         USOL(I,L) = USOL(I,1)
  390 CONTINUE
  400 CONTINUE
C
C     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
C     NORM IF OPERATOR IS SINGULAR
C
      IF (SINGLR) CALL SEPMIN (USOL,IDMN,ZN,ZM,PRTRB)
C
C     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
C     NOT FLAGGED
C
      IF (IORD .EQ. 2) RETURN
      IORD = 2
C
C     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
C
      CALL D4FER(COFX,IDMN,USOL,GRHS)
      GO TO 360
      END
      SUBROUTINE C4KPRM(IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,IDMN,IERR
     1OR)
C
C     THIS PROGRAM CHECKS THE INPUT PARAMETERS FOR ERRORS
C
C
C
C     CHECK DEFINITION OF SOLUTION REGION
C
      IERROR = 1
      IF (A.GE.B .OR. C.GE.D) RETURN
C
C     CHECK BOUNDARY SWITCHES
C
      IERROR = 2
      IF (MBDCND.LT.0 .OR. MBDCND.GT.4) RETURN
      IERROR = 3
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) RETURN
C
C     CHECK FIRST DIMENSION IN CALLING ROUTINE
C
      IERROR = 5
      IF (IDMN .LT. 7) RETURN
C
C     CHECK M
C
      IERROR = 6
      IF (M.GT.(IDMN-1) .OR. M.LT.6) RETURN
C
C     CHECK N
C
      IERROR = 7
      IF (N .LT. 5) RETURN
C
C     CHECK IORDER
C
      IERROR = 8
      IF (IORDER.NE.2 .AND. IORDER.NE.4) RETURN
C
C     CHECK INTL
C
C
C     CHECK THAT EQUATION IS ELLIPTIC
C
      DLX = (B-A)/FLOAT(M)
      DO  30 I=2,M
         XI = A+FLOAT(I-1)*DLX
         CALL COFX (XI,AI,BI,CI)
      IF (AI.GT.0.0) GO TO 10
      IERROR=10
      RETURN
   10 CONTINUE
   30 CONTINUE
C
C     NO ERROR FOUND
C
      IERROR = 0
      RETURN
      END
      SUBROUTINE C4KSNG(MBDCND,NBDCND,ALPHA,BETA,COFX,SINGLR)
C
C     THIS SUBROUTINE CHECKS IF THE PDE   SEPELI
C     MUST SOLVE IS A SINGULAR OPERATOR
C
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      LOGICAL         SINGLR
      SINGLR = .FALSE.
C
C     CHECK IF THE BOUNDARY CONDITIONS ARE
C     ENTIRELY PERIODIC AND/OR MIXED
C
      IF ((MBDCND.NE.0 .AND. MBDCND.NE.3) .OR.
     1    (NBDCND.NE.0 .AND. NBDCND.NE.3)) RETURN
C
C     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
C
      IF (MBDCND .NE. 3) GO TO  10
      IF (ALPHA.NE.0.0 .OR. BETA.NE.0.0) RETURN
   10 CONTINUE
C
C     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
C     ARE ZERO
C
      DO  30 I=IS,MS
         XI = AIT+FLOAT(I-1)*DLX
         CALL COFX (XI,AI,BI,CI)
         IF (CI .NE. 0.0) RETURN
   30 CONTINUE
C
C     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
C
      SINGLR = .TRUE.
      RETURN
      END
      SUBROUTINE D4FER(COFX,IDMN,USOL,GRHS)
C
C     THIS SUBROUTINE FIRST APPROXIMATES THE TRUNCATION ERROR GIVEN BY
C     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY WHERE
C     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 ON THE INTERIOR AND
C     AT THE BOUNDARIES IF PERIODIC(HERE UXXX,UXXXX ARE THE THIRD
C     AND FOURTH PARTIAL DERIVATIVES OF U WITH RESPECT TO X).
C     TX IS OF THE FORM AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
C     AT X=A OR X=B IF THE BOUNDARY CONDITION THERE IS MIXED.
C     TX=0.0 ALONG SPECIFIED BOUNDARIES.  TY HAS SYMMETRIC FORM
C     IN Y WITH X,AFUN(X),BFUN(X) REPLACED BY Y,DFUN(Y),EFUN(Y).
C     THE SECOND ORDER SOLUTION IN USOL IS USED TO APPROXIMATE
C     (VIA SECOND ORDER FINITE DIFFERENCING) THE TRUNCATION ERROR
C     AND THE RESULT IS ADDED TO THE RIGHT HAND SIDE IN GRHS
C     AND THEN TRANSFERRED TO USOL TO BE USED AS A NEW RIGHT
C     HAND SIDE WHEN CALLING BLKTRI FOR A FOURTH ORDER SOLUTION.
C
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,
     1                AIT        ,BIT        ,CIT        ,DIT        ,
     2                MIT        ,NIT        ,IS         ,MS         ,
     3                JS         ,NS         ,DLX        ,DLY        ,
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4
      DIMENSION       GRHS(IDMN,1)           ,USOL(IDMN,1)
C
C
C     COMPUTE TRUNCATION ERROR APPROXIMATION OVER THE ENTIRE MESH
C
         DO  30 I=IS,MS
            XI = AIT+FLOAT(I-1)*DLX
            CALL COFX (XI,AI,BI,CI)
         DO 30 J=JS,NS
C
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
C
            CALL SEPDX (USOL,IDMN,I,J,UXXX,UXXXX)
            CALL SEPDY (USOL,IDMN,I,J,UYYY,UYYYY)
            TX = AI*UXXXX/12.0+BI*UXXX/6.0
             TY=UYYYY/12.0
C
C     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
C
            IF (KSWX.EQ.1 .OR. (I.GT.1 .AND. I.LT.K)) GO TO  10
            TX = AI/3.0*(UXXXX/4.0+UXXX/DLX)
   10       IF (KSWY.EQ.1 .OR. (J.GT.1 .AND. J.LT.L)) GO TO  20
            TY = (UYYYY/4.0+UYYY/DLY)/3.0
   20 GRHS(I,J)=GRHS(I,J)+DLY**2*(DLX**2*TX+DLY**2*TY)
   30    CONTINUE
C
C     RESET THE RIGHT HAND SIDE IN USOL
C
      DO  60 I=IS,MS
         DO  50 J=JS,NS
            USOL(I,J) = GRHS(I,J)
   50    CONTINUE
   60 CONTINUE
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
