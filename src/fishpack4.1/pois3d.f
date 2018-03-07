C
C     file pois3d.f
C
      SUBROUTINE POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF,
     1                   MDIMF,F,IERROR,W)
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
C DIMENSION OF           A(N), B(N), C(N), F(LDIMF,MDIMF,N),
C ARGUMENTS              W(SEE ARGUMENT LIST)
C
C LATEST REVISION        NOVEMBER 1988
C
C PURPOSE                SOLVES THE LINEAR SYSTEM OF EQUATIONS
C                        FOR UNKNOWN X VALUES, WHERE I=1,2,...,L,
C                        J=1,2,...,M, AND K=1,2,...,N
C
C                        C1*(X(I-1,J,K) -2.*X(I,J,K) +X(I+1,J,K)) +
C                        C2*(X(I,J-1,K) -2.*X(I,J,K) +X(I,J+1,K)) +
C                        A(K)*X(I,J,K-1) +B(K)*X(I,J,K)+ C(K)*X(I,J,K+1)
C                        = F(I,J,K)
C
C                        THE INDICES K-1 AND K+1 ARE EVALUATED MODULO N,
C                        I.E. X(I,J,0)=X(I,J,N) AND X(I,J,N+1)=X(I,J,1).
C                        THE UNKNOWNS
C                        X(0,J,K), X(L+1,J,K), X(I,0,K), AND X(I,M+1,K)
C                        ARE ASSUMED TO TAKE ON CERTAIN PRESCRIBED
C                        VALUES DESCRIBED BELOW.
C
C USAGE                  CALL POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,
C                        N,A,B,C,LDIMF,MDIMF,F,IERROR,W)
C
C ARGUMENTS
C
C ON INPUT
C                        LPEROD
C                          INDICATES THE VALUES THAT X(0,J,K) AND
C                          X(L+1,J,K) ARE ASSUMED TO HAVE.
C                          = 0  X(0,J,K)=X(L,J,K), X(L+1,J,K)=X(1,J,K)
C                          = 1  X(0,J,K) = 0,      X(L+1,J,K) = 0
C                          = 2  X(0,J,K)=0,        X(L+1,J,K)=X(L-1,J,K)
C                          = 3  X(0,J,K)=X(2,J,K), X(L+1,J,K)=X(L-1,J,K)
C                          = 4  X(0,J,K)=X(2,J,K), X(L+1,J,K) = 0.
C
C                        L
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
C                          L MUST BE AT LEAST 3.
C
C                        C1
C                          REAL CONSTANT IN THE ABOVE LINEAR SYSTEM
C                          OF EQUATIONS TO BE SOLVED.
C
C                        MPEROD
C                          INDICATES THE VALUES THAT X(I,0,K) AND
C                          X(I,M+1,K) ARE ASSUMED TO HAVE.
C                          = 0  X(I,0,K)=X(I,M,K), X(I,M+1,K)=X(I,1,K)
C                          = 1  X(I,0,K)=0,        X(I,M+1,K)=0
C                          = 2  X(I,0,K)=0,        X(I,M+1,K)=X(I,M-1,K)
C                          = 3  X(I,0,K)=X(I,2,K)  X(I,M+1,K)=X(I,M-1,K)
C                          = 4  X(I,0,K)=X(I,2,K)  X(I,M+1,K)=0
C
C                        M
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
C                          M MUST BE AT LEAST 3.
C
C                        C2
C                          REAL CONSTANT IN THE ABOVE LINEAR SYSTEM
C                          OF EQUATIONS TO BE SOLVED.
C
C                        NPEROD
C                          = 0  IF A(1) AND C(N) ARE NOT ZERO.
C                          = 1  IF A(1) = C(N) = 0.
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE K-DIRECTION.
C                          N MUST BE AT LEAST 3.
C
C                        A, B, C
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH N THAT
C                          SPECIFY THE COEFFICIENTS IN THE LINEAR
C                          EQUATIONS GIVEN ABOVE.
C
C                          IF NPEROD = 0 THE ARRAY ELEMENTS MUST NOT
C                          DEPEND UPON INDEX K, BUT MUST BE CONSTANT.
C                          SPECIFICALLY,THE SUBROUTINE CHECKS THE
C                          FOLLOWING CONDITION
C                            A(K) = C(1)
C                            C(K) = C(1)
C                            B(K) = B(1)
C                          FOR K=1,2,...,N.
C
C                        LDIMF
C                          THE ROW (OR FIRST) DIMENSION OF THE THREE-
C                          DIMENSIONAL ARRAY F AS IT APPEARS IN THE
C                          PROGRAM CALLING POIS3D.  THIS PARAMETER IS
C                          USED TO SPECIFY THE VARIABLE DIMENSION
C                          OF F.  LDIMF MUST BE AT LEAST L.
C
C                        MDIMF
C                          THE COLUMN (OR SECOND) DIMENSION OF THE THREE
C                          DIMENSIONAL ARRAY F AS IT APPEARS IN THE
C                          PROGRAM CALLING POIS3D.  THIS PARAMETER IS
C                          USED TO SPECIFY THE VARIABLE DIMENSION
C                          OF F.  MDIMF MUST BE AT LEAST M.
C
C                        F
C                          A THREE-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
C                          OF EQUATIONS GIVEN ABOVE.  F MUST BE
C                          DIMENSIONED AT LEAST L X M X N.
C
C                        W
C                          A ONE-DIMENSIONAL ARRAY THAT MUST BE PROVIDED
C                          BY THE USER FOR WORK SPACE.  THE LENGTH OF W
C                          MUST BE AT LEAST
C                            30 + L + M + 2*N + MAX(L,M,N) +
C                            7*(INT((L+1)/2) + INT((M+1)/2)).
C
C ON OUTPUT
C
C                        F
C                          CONTAINS THE SOLUTION X.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBER ZERO, A
C                          SOLUTION IS NOT ATTEMPTED.
C                          = 0  NO ERROR
C                          = 1  IF LPEROD .LT. 0 OR .GT. 4
C                          = 2  IF L .LT. 3
C                          = 3  IF MPEROD .LT. 0 OR .GT. 4
C                          = 4  IF M .LT. 3
C                          = 5  IF NPEROD .LT. 0 OR .GT. 1
C                          = 6  IF N .LT. 3
C                          = 7  IF LDIMF .LT. L
C                          = 8  IF MDIMF .LT. M
C                          = 9  IF A(K) .NE. C(1) OR C(K) .NE. C(1)
C                               OR B(I) .NE.B(1) FOR SOME K=1,2,...,N.
C                          = 10 IF NPEROD = 1 AND A(1) .NE. 0
C                               OR C(N) .NE. 0
C
C                          SINCE THIS IS THE ONLY MEANS OF INDICATING A
C                          POSSIBLY INCORRECT CALL TO POIS3D, THE USER
C                          SHOULD TEST IERROR AFTER THE CALL.
C
C SPECIAL CONDITIONS     NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       COMF AND FFTPACK FROM FISHPACK
C FILES
C
C LANGUAGE               FORTRAN
C
C HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
C                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
C                        LIBRARIES IN JANUARY, 1980.
C
C PORTABILITY            FORTRAN 77
C
C ALGORITHM              THIS SUBROUTINE SOLVES THREE-DIMENSIONAL BLOCK
C                        TRIDIAGONAL LINEAR SYSTEMS ARISING FROM FINITE
C                        DIFFERENCE APPROXIMATIONS TO THREE-DIMENSIONAL
C                        POISSON EQUATIONS USING THE FFT PACKAGE
C                        FFTPACK WRITTEN BY PAUL SWARZTRAUBER.
C
C TIMING                 FOR LARGE L, M AND N, THE OPERATION COUNT
C                        IS ROUGHLY PROPORTIONAL TO
C                          L*M*N*(LOG2(L)+LOG2(M)+5)
C                        BUT ALSO DEPENDS ON INPUT PARAMETERS LPEROD
C                        AND MPEROD.
C
C ACCURACY               TO MEASURE THE ACCURACY OF THE ALGORITHM A
C                        UNIFORM RANDOM NUMBER GENERATOR WAS USED TO
C                        CREATE A SOLUTION ARRAY X FOR THE SYSTEM GIVEN
C                        IN THE 'PURPOSE' SECTION WITH
C                          A(K) = C(K) = -0.5*B(K) = 1,  K=1,2,...,N
C                        AND, WHEN NPEROD = 1
C                          A(1) = C(N) = 0
C                          A(N) = C(1) = 2.
C
C                        THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN
C                        SYSTEM AND, USING DOUBLE PRECISION, A RIGHT
C                        SIDE Y WAS COMPUTED.  USING THIS ARRAY Y
C                        SUBROUTINE POIS3D WAS CALLED TO PRODUCE AN
C                        APPROXIMATE SOLUTION Z.  RELATIVE ERROR
C
C                        E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K
C
C                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
C                        OVER I=1,2,...,L, J=1,2,...,M AND K=1,2,...,N.
C                        VALUES OF E ARE GIVEN IN THE TABLE BELOW FOR
C                        SOME TYPICAL VALUES OF L,M AND N.
C
C                        L(=M=N)   LPEROD    MPEROD       E
C                        ------    ------    ------     ------
C
C                          16        0         0        1.E-13
C                          15        1         1        4.E-13
C                          17        3         3        2.E-13
C                          32        0         0        2.E-13
C                          31        1         1        2.E-12
C                          33        3         3        7.E-13
C
C REFERENCES              NONE
C ********************************************************************
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,
     1                F(LDIMF,MDIMF,1)       ,W(*)       ,SAVE(6)
C
      LP = LPEROD+1
      MP = MPEROD+1
      NP = NPEROD+1
C
C     CHECK FOR INVALID INPUT.
C
      IERROR = 0
      IF (LP.LT.1 .OR. LP.GT.5) IERROR = 1
      IF (L .LT. 3) IERROR = 2
      IF (MP.LT.1 .OR. MP.GT.5) IERROR = 3
      IF (M .LT. 3) IERROR = 4
      IF (NP.LT.1 .OR. NP.GT.2) IERROR = 5
      IF (N .LT. 3) IERROR = 6
      IF (LDIMF .LT. L) IERROR = 7
      IF (MDIMF .LT. M) IERROR = 8
      IF (NP .NE. 1) GO TO 103
      DO 101 K=1,N
         IF (A(K) .NE. C(1)) GO TO 102
         IF (C(K) .NE. C(1)) GO TO 102
         IF (B(K) .NE. B(1)) GO TO 102
  101 CONTINUE
      GO TO 104
  102 IERROR = 9
  103 IF (NPEROD.EQ.1 .AND. (A(1).NE.0. .OR. C(N).NE.0.)) IERROR = 10
  104 IF (IERROR .NE. 0) GO TO 122
      IWYRT = L+1
      IWT = IWYRT+M
      IWD = IWT+MAX0(L,M,N)+1
      IWBB = IWD+N
      IWX = IWBB+N
      IWY = IWX+7*((L+1)/2)+15
      GO TO (105,114),NP
C
C     REORDER UNKNOWNS WHEN NPEROD = 0.
C
  105 NH = (N+1)/2
      NHM1 = NH-1
      NODD = 1
      IF (2*NH .EQ. N) NODD = 2
      DO 111 I=1,L
         DO 110 J=1,M
            DO 106 K=1,NHM1
               NHPK = NH+K
               NHMK = NH-K
               W(K) = F(I,J,NHMK)-F(I,J,NHPK)
               W(NHPK) = F(I,J,NHMK)+F(I,J,NHPK)
  106       CONTINUE
            W(NH) = 2.*F(I,J,NH)
            GO TO (108,107),NODD
  107       W(N) = 2.*F(I,J,N)
  108       DO 109 K=1,N
               F(I,J,K) = W(K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      SAVE(1) = C(NHM1)
      SAVE(2) = A(NH)
      SAVE(3) = C(NH)
      SAVE(4) = B(NHM1)
      SAVE(5) = B(N)
      SAVE(6) = A(N)
      C(NHM1) = 0.
      A(NH) = 0.
      C(NH) = 2.*C(NH)
      GO TO (112,113),NODD
  112 B(NHM1) = B(NHM1)-A(NH-1)
      B(N) = B(N)+A(N)
      GO TO 114
  113 A(N) = C(NH)
  114 CONTINUE
      CALL POS3D1 (LP,L,MP,M,N,A,B,C,LDIMF,MDIMF,F,W,W(IWYRT),W(IWT),
     1             W(IWD),W(IWX),W(IWY),C1,C2,W(IWBB))
      GO TO (115,122),NP
  115 DO 121 I=1,L
         DO 120 J=1,M
            DO 116 K=1,NHM1
               NHMK = NH-K
               NHPK = NH+K
               W(NHMK) = .5*(F(I,J,NHPK)+F(I,J,K))
               W(NHPK) = .5*(F(I,J,NHPK)-F(I,J,K))
  116       CONTINUE
            W(NH) = .5*F(I,J,NH)
            GO TO (118,117),NODD
  117       W(N) = .5*F(I,J,N)
  118       DO 119 K=1,N
               F(I,J,K) = W(K)
  119       CONTINUE
  120    CONTINUE
  121 CONTINUE
      C(NHM1) = SAVE(1)
      A(NH) = SAVE(2)
      C(NH) = SAVE(3)
      B(NHM1) = SAVE(4)
      B(N) = SAVE(5)
      A(N) = SAVE(6)
  122 CONTINUE
      RETURN
      END
      SUBROUTINE POS3D1 (LP,L,MP,M,N,A,B,C,LDIMF,MDIMF,F,XRT,YRT,T,D,
     1                   WX,WY,C1,C2,BB)
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,
     1                F(LDIMF,MDIMF,1)       ,XRT(*)     ,YRT(*)     ,
     2                T(*)       ,D(*)       ,WX(*)      ,WY(*)      ,
     3                BB(*)
      PI = PIMACH(DUM)
      LR = L
      MR = M
      NR = N
C
C     GENERATE TRANSFORM ROOTS
C
      LRDEL = ((LP-1)*(LP-3)*(LP-5))/3
      SCALX = LR+LRDEL
      DX = PI/(2.*SCALX)
      GO TO (108,103,101,102,101),LP
  101 DI = 0.5
      SCALX = 2.*SCALX
      GO TO 104
  102 DI = 1.0
      GO TO 104
  103 DI = 0.0
  104 DO 105 I=1,LR
         XRT(I) = -4.*C1*(SIN((FLOAT(I)-DI)*DX))**2
  105 CONTINUE
      SCALX = 2.*SCALX
      GO TO (112,106,110,107,111),LP
  106 CALL SINTI (LR,WX)
      GO TO 112
  107 CALL COSTI (LR,WX)
      GO TO 112
  108 XRT(1) = 0.
      XRT(LR) = -4.*C1
      DO 109 I=3,LR,2
         XRT(I-1) = -4.*C1*(SIN(FLOAT((I-1))*DX))**2
         XRT(I) = XRT(I-1)
  109 CONTINUE
      CALL RFFTI (LR,WX)
      GO TO 112
  110 CALL SINQI (LR,WX)
      GO TO 112
  111 CALL COSQI (LR,WX)
  112 CONTINUE
      MRDEL = ((MP-1)*(MP-3)*(MP-5))/3
      SCALY = MR+MRDEL
      DY = PI/(2.*SCALY)
      GO TO (120,115,113,114,113),MP
  113 DJ = 0.5
      SCALY = 2.*SCALY
      GO TO 116
  114 DJ = 1.0
      GO TO 116
  115 DJ = 0.0
  116 DO 117 J=1,MR
         YRT(J) = -4.*C2*(SIN((FLOAT(J)-DJ)*DY))**2
  117 CONTINUE
      SCALY = 2.*SCALY
      GO TO (124,118,122,119,123),MP
  118 CALL SINTI (MR,WY)
      GO TO 124
  119 CALL COSTI (MR,WY)
      GO TO 124
  120 YRT(1) = 0.
      YRT(MR) = -4.*C2
      DO 121 J=3,MR,2
         YRT(J-1) = -4.*C2*(SIN(FLOAT((J-1))*DY))**2
         YRT(J) = YRT(J-1)
  121 CONTINUE
      CALL RFFTI (MR,WY)
      GO TO 124
  122 CALL SINQI (MR,WY)
      GO TO 124
  123 CALL COSQI (MR,WY)
  124 CONTINUE
      IFWRD = 1
      IS = 1
  125 CONTINUE
C
C     TRANSFORM X
C
      DO 141 J=1,MR
         DO 140 K=1,NR
            DO 126 I=1,LR
               T(I) = F(I,J,K)
  126       CONTINUE
            GO TO (127,130,131,134,135),LP
  127       GO TO (128,129),IFWRD
  128       CALL RFFTF (LR,T,WX)
            GO TO 138
  129       CALL RFFTB (LR,T,WX)
            GO TO 138
  130       CALL SINT (LR,T,WX)
            GO TO 138
  131       GO TO (132,133),IFWRD
  132       CALL SINQF (LR,T,WX)
            GO TO 138
  133       CALL SINQB (LR,T,WX)
            GO TO 138
  134       CALL COST (LR,T,WX)
            GO TO 138
  135       GO TO (136,137),IFWRD
  136       CALL COSQF (LR,T,WX)
            GO TO 138
  137       CALL COSQB (LR,T,WX)
  138       CONTINUE
            DO 139 I=1,LR
               F(I,J,K) = T(I)
  139       CONTINUE
  140    CONTINUE
  141 CONTINUE
      GO TO (142,164),IFWRD
C
C     TRANSFORM Y
C
  142 CONTINUE
      DO 158 I=1,LR
         DO 157 K=1,NR
            DO 143 J=1,MR
               T(J) = F(I,J,K)
  143       CONTINUE
            GO TO (144,147,148,151,152),MP
  144       GO TO (145,146),IFWRD
  145       CALL RFFTF (MR,T,WY)
            GO TO 155
  146       CALL RFFTB (MR,T,WY)
            GO TO 155
  147       CALL SINT (MR,T,WY)
            GO TO 155
  148       GO TO (149,150),IFWRD
  149       CALL SINQF (MR,T,WY)
            GO TO 155
  150       CALL SINQB (MR,T,WY)
            GO TO 155
  151       CALL COST (MR,T,WY)
            GO TO 155
  152       GO TO (153,154),IFWRD
  153       CALL COSQF (MR,T,WY)
            GO TO 155
  154       CALL COSQB (MR,T,WY)
  155       CONTINUE
            DO 156 J=1,MR
               F(I,J,K) = T(J)
  156       CONTINUE
  157    CONTINUE
  158 CONTINUE
      GO TO (159,125),IFWRD
  159 CONTINUE
C
C     SOLVE TRIDIAGONAL SYSTEMS IN Z
C
      DO 163 I=1,LR
         DO 162 J=1,MR
            DO 160 K=1,NR
               BB(K) = B(K)+XRT(I)+YRT(J)
               T(K) = F(I,J,K)
  160       CONTINUE
            CALL TRID (NR,A,BB,C,T,D)
            DO 161 K=1,NR
               F(I,J,K) = T(K)
  161       CONTINUE
  162    CONTINUE
  163 CONTINUE
      IFWRD = 2
      IS = -1
      GO TO 142
  164 CONTINUE
      DO 167 I=1,LR
         DO 166 J=1,MR
            DO 165 K=1,NR
               F(I,J,K) = F(I,J,K)/(SCALX*SCALY)
  165       CONTINUE
  166    CONTINUE
  167 CONTINUE
      RETURN
      END
      SUBROUTINE TRID (MR,A,B,C,Y,D)
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,Y(*)       ,
     1                D(*)
      M = MR
      MM1 = M-1
      Z = 1./B(1)
      D(1) = C(1)*Z
      Y(1) = Y(1)*Z
      DO 101 I=2,MM1
         Z = 1./(B(I)-A(I)*D(I-1))
         D(I) = C(I)*Z
         Y(I) = (Y(I)-A(I)*Y(I-1))*Z
  101 CONTINUE
      Z = B(M)-A(M)*D(MM1)
      IF (Z .NE. 0.) GO TO 102
      Y(M) = 0.
      GO TO 103
  102 Y(M) = (Y(M)-A(M)*Y(MM1))/Z
  103 CONTINUE
      DO 104 IP=1,MM1
         I = M-IP
         Y(I) = Y(I)-D(I)*Y(I+1)
  104 CONTINUE
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
