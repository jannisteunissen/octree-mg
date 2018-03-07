C
C     file poistg.f
C
      SUBROUTINE POISTG (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,W)
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
C DIMENSION OF           A(M),  B(M),  C(M),  Y(IDIMY,N),
C ARGUMENTS              W(SEE ARGUMENT LIST)
C
C LATEST REVISION        NOVEMBER 1988
C
C PURPOSE                SOLVES THE LINEAR SYSTEM OF EQUATIONS
C                        FOR UNKNOWN X VALUES, WHERE I=1,2,...,M
C                        AND J=1,2,...,N
C
C                        A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
C                        + X(I,J-1) - 2.*X(I,J) + X(I,J+1)
C                        = Y(I,J)
C
C                        THE INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
C                        I.E. X(0,J) = X(M,J) AND X(M+1,J) = X(1,J), AND
C                        X(I,0) MAY BE EQUAL TO X(I,1) OR -X(I,1), AND
C                        X(I,N+1) MAY BE EQUAL TO X(I,N) OR -X(I,N),
C                        DEPENDING ON AN INPUT PARAMETER.
C
C USAGE                  CALL POISTG (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,
C                                     IERROR,W)
C
C ARGUMENTS
C
C ON INPUT
C
C                        NPEROD
C                          INDICATES VALUES WHICH X(I,0) AND X(I,N+1)
C                          ARE ASSUMED TO HAVE.
C                          = 1 IF X(I,0) = -X(I,1) AND X(I,N+1) = -X(I,N
C                          = 2 IF X(I,0) = -X(I,1) AND X(I,N+1) =  X(I,N
C                          = 3 IF X(I,0) =  X(I,1) AND X(I,N+1) =  X(I,N
C                          = 4 IF X(I,0) =  X(I,1) AND X(I,N+1) = -X(I,N
C
C                        N
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
C                          N MUST BE GREATER THAN 2.
C
C                        MPEROD
C                          = 0 IF A(1) AND C(M) ARE NOT ZERO
C                          = 1 IF A(1) = C(M) = 0
C
C                        M
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
C                          M MUST BE GREATER THAN 2.
C
C                        A,B,C
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
C                          SPECIFY THE COEFFICIENTS IN THE LINEAR
C                          EQUATIONS GIVEN ABOVE.  IF MPEROD = 0 THE
C                          ARRAY ELEMENTS MUST NOT DEPEND ON INDEX I,
C                          BUT MUST BE CONSTANT.  SPECIFICALLY, THE
C                          SUBROUTINE CHECKS THE FOLLOWING CONDITION
C                            A(I) = C(1)
C                            B(I) = B(1)
C                            C(I) = C(1)
C                          FOR I = 1, 2, ..., M.
C
C                        IDIMY
C                          THE ROW (OR FIRST) DIMENSION OF THE TWO-
C                          DIMENSIONAL ARRAY Y AS IT APPEARS IN THE
C                          PROGRAM CALLING POISTG.  THIS PARAMETER IS
C                          USED TO SPECIFY THE VARIABLE DIMENSION OF Y.
C                          IDIMY MUST BE AT LEAST M.
C
C                        Y
C                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
C                          VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
C                          OF EQUATIONS GIVEN ABOVE.
C                          Y MUST BE DIMENSIONED AT LEAST M X N.
C
C                        W
C                          A ONE-DIMENSIONAL WORK ARRAY THAT MUST BE
C                          PROVIDED BY THE USER FOR WORK SPACE.  W MAY
C                          REQUIRE UP TO 9M + 4N + M(INT(LOG2(N)))
C                          LOCATIONS.  THE ACTUAL NUMBER OF LOCATIONS
C                          USED IS COMPUTED BY POISTG AND RETURNED IN
C                          LOCATION W(1).
C
C ON OUTPUT
C
C                        Y
C                          CONTAINS THE SOLUTION X.
C
C                        IERROR
C                          AN ERROR FLAG THAT INDICATES INVALID INPUT
C                          PARAMETERS.  EXCEPT FOR NUMBER ZERO, A
C                          SOLUTION IS NOT ATTEMPTED.
C                          = 0  NO ERROR
C                          = 1  IF M .LE. 2
C                          = 2  IF N .LE. 2
C                          = 3  IDIMY .LT. M
C                          = 4  IF NPEROD .LT. 1 OR NPEROD .GT. 4
C                          = 5  IF MPEROD .LT. 0 OR MPEROD .GT. 1
C                          = 6  IF MPEROD = 0 AND A(I) .NE. C(1)
C                               OR B(I) .NE. B(1) OR C(I) .NE. C(1)
C                               FOR SOME I = 1, 2, ..., M.
C                          = 7  IF MPEROD .EQ. 1 .AND.
C                               (A(1).NE.0 .OR. C(M).NE.0)
C
C                        W
C                          W(1) CONTAINS THE REQUIRED LENGTH OF W.
C
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED LIBRARY       GNBNAUX AND COMF FROM FISHPACK
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
C ALGORITHM              THIS SUBROUTINE IS AN IMPLEMENTATION OF THE
C                        ALGORITHM PRESENTED IN THE REFERENCE BELOW.
C
C TIMING                 FOR LARGE M AND N, THE EXECUTION TIME IS
C                        ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
C
C ACCURACY               TO MEASURE THE ACCURACY OF THE ALGORITHM A
C                        UNIFORM RANDOM NUMBER GENERATOR WAS USED TO
C                        CREATE A SOLUTION ARRAY X FOR THE SYSTEM GIVEN
C                        IN THE 'PURPOSE' SECTION ABOVE, WITH
C                          A(I) = C(I) = -0.5*B(I) = 1,    I=1,2,...,M
C                        AND, WHEN MPEROD = 1
C                          A(1) = C(M) = 0
C                          B(1) = B(M) =-1.
C
C                        THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN
C                        SYSTEM AND, USING DOUBLE PRECISION, A RIGHT SID
C                        Y WAS COMPUTED.  USING THIS ARRAY Y SUBROUTINE
C                        POISTG WAS CALLED TO PRODUCE AN APPROXIMATE
C                        SOLUTION Z.  THEN THE RELATIVE ERROR, DEFINED A
C                          E = MAX(ABS(Z(I,J)-X(I,J)))/MAX(ABS(X(I,J)))
C                        WHERE THE TWO MAXIMA ARE TAKEN OVER I=1,2,...,M
C                        AND J=1,2,...,N, WAS COMPUTED.  VALUES OF E ARE
C                        GIVEN IN THE TABLE BELOW FOR SOME TYPICAL VALUE
C                        OF M AND N.
C
C                        M (=N)    MPEROD    NPEROD      E
C                        ------    ------    ------    ------
C
C                          31        0-1       1-4     9.E-13
C                          31        1         1       4.E-13
C                          31        1         3       3.E-13
C                          32        0-1       1-4     3.E-12
C                          32        1         1       3.E-13
C                          32        1         3       1.E-13
C                          33        0-1       1-4     1.E-12
C                          33        1         1       4.E-13
C                          33        1         3       1.E-13
C                          63        0-1       1-4     3.E-12
C                          63        1         1       1.E-12
C                          63        1         3       2.E-13
C                          64        0-1       1-4     4.E-12
C                          64        1         1       1.E-12
C                          64        1         3       6.E-13
C                          65        0-1       1-4     2.E-13
C                          65        1         1       1.E-11
C                          65        1         3       4.E-13
C
C REFERENCES             SCHUMANN, U. AND R. SWEET,"A DIRECT METHOD
C                        FOR THE SOLUTION OF POISSON"S EQUATION WITH
C                        NEUMANN BOUNDARY CONDITIONS ON A STAGGERED
C                        GRID OF ARBITRARY SIZE," J. COMP. PHYS.
C                        20(1976), PP. 171-182.
C *********************************************************************
      DIMENSION       Y(IDIMY,1)
      DIMENSION       W(*)       ,B(*)       ,A(*)       ,C(*)
C
      IERROR = 0
      IF (M .LE. 2) IERROR = 1
      IF (N .LE. 2) IERROR = 2
      IF (IDIMY .LT. M) IERROR = 3
      IF (NPEROD.LT.1 .OR. NPEROD.GT.4) IERROR = 4
      IF (MPEROD.LT.0 .OR. MPEROD.GT.1) IERROR = 5
      IF (MPEROD .EQ. 1) GO TO 103
      DO 101 I=1,M
         IF (A(I) .NE. C(1)) GO TO 102
         IF (C(I) .NE. C(1)) GO TO 102
         IF (B(I) .NE. B(1)) GO TO 102
  101 CONTINUE
      GO TO 104
  102 IERROR = 6
      RETURN
  103 IF (A(1).NE.0. .OR. C(M).NE.0.) IERROR = 7
  104 IF (IERROR .NE. 0) RETURN
      IWBA = M+1
      IWBB = IWBA+M
      IWBC = IWBB+M
      IWB2 = IWBC+M
      IWB3 = IWB2+M
      IWW1 = IWB3+M
      IWW2 = IWW1+M
      IWW3 = IWW2+M
      IWD = IWW3+M
      IWTCOS = IWD+M
      IWP = IWTCOS+4*N
      DO 106 I=1,M
         K = IWBA+I-1
         W(K) = -A(I)
         K = IWBC+I-1
         W(K) = -C(I)
         K = IWBB+I-1
         W(K) = 2.-B(I)
         DO 105 J=1,N
            Y(I,J) = -Y(I,J)
  105    CONTINUE
  106 CONTINUE
      NP = NPEROD
      MP = MPEROD+1
      GO TO (110,107),MP
  107 CONTINUE
      GO TO (108,108,108,119),NPEROD
  108 CONTINUE
      CALL POSTG2 (NP,N,M,W(IWBA),W(IWBB),W(IWBC),IDIMY,Y,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD .EQ. 4) GO TO 120
  109 CONTINUE
      GO TO (123,129),MP
  110 CONTINUE
C
C     REORDER UNKNOWNS WHEN MP =0
C
      MH = (M+1)/2
      MHM1 = MH-1
      MODD = 1
      IF (MH*2 .EQ. M) MODD = 2
      DO 115 J=1,N
         DO 111 I=1,MHM1
            MHPI = MH+I
            MHMI = MH-I
            W(I) = Y(MHMI,J)-Y(MHPI,J)
            W(MHPI) = Y(MHMI,J)+Y(MHPI,J)
  111    CONTINUE
         W(MH) = 2.*Y(MH,J)
         GO TO (113,112),MODD
  112    W(M) = 2.*Y(M,J)
  113    CONTINUE
         DO 114 I=1,M
            Y(I,J) = W(I)
  114    CONTINUE
  115 CONTINUE
      K = IWBC+MHM1-1
      I = IWBA+MHM1
      W(K) = 0.
      W(I) = 0.
      W(K+1) = 2.*W(K+1)
      GO TO (116,117),MODD
  116 CONTINUE
      K = IWBB+MHM1-1
      W(K) = W(K)-W(I-1)
      W(IWBC-1) = W(IWBC-1)+W(IWBB-1)
      GO TO 118
  117 W(IWBB-1) = W(K+1)
  118 CONTINUE
      GO TO 107
  119 CONTINUE
C
C     REVERSE COLUMNS WHEN NPEROD = 4.
C
      IREV = 1
      NBY2 = N/2
      NP = 2
  120 DO 122 J=1,NBY2
         MSKIP = N+1-J
         DO 121 I=1,M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
  121    CONTINUE
  122 CONTINUE
      GO TO (108,109),IREV
  123 CONTINUE
      DO 128 J=1,N
         DO 124 I=1,MHM1
            MHMI = MH-I
            MHPI = MH+I
            W(MHMI) = .5*(Y(MHPI,J)+Y(I,J))
            W(MHPI) = .5*(Y(MHPI,J)-Y(I,J))
  124    CONTINUE
         W(MH) = .5*Y(MH,J)
         GO TO (126,125),MODD
  125    W(M) = .5*Y(M,J)
  126    CONTINUE
         DO 127 I=1,M
            Y(I,J) = W(I)
  127    CONTINUE
  128 CONTINUE
  129 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR W ARRAY.
C
      W(1) = IPSTOR+IWP-1
      RETURN
      END
      SUBROUTINE POSTG2 (NPEROD,N,M,A,BB,C,IDIMQ,Q,B,B2,B3,W,W2,W3,D,
     1                   TCOS,P)
C
C     SUBROUTINE TO SOLVE POISSON'S EQUATION ON A STAGGERED GRID.
C
C
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                K(4)       ,P(*)
      EQUIVALENCE     (K(1),K1)  ,(K(2),K2)  ,(K(3),K3)  ,(K(4),K4)
      NP = NPEROD
      FNUM = 0.5*FLOAT(NP/3)
      FNUM2 = 0.5*FLOAT(NP/2)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      IF (NR .LE. 3) GO TO 142
  101 CONTINUE
      JR = 2*I2R
      NROD = 1
      IF ((NR/2)*2 .EQ. NR) NROD = 0
      JSTART = 1
      JSTOP = NLAST-JR
      IF (NROD .EQ. 0) JSTOP = JSTOP-I2R
      I2RBY2 = I2R/2
      IF (JSTOP .GE. JSTART) GO TO 102
      J = JR
      GO TO 115
  102 CONTINUE
C
C     REGULAR REDUCTION.
C
      IJUMP = 1
      DO 114 J=JSTART,JSTOP,JR
         JP1 = J+I2RBY2
         JP2 = J+I2R
         JP3 = JP2+I2RBY2
         JM1 = J-I2RBY2
         JM2 = J-I2R
         JM3 = JM2-I2RBY2
         IF (J .NE. 1) GO TO 106
         CALL COSGEN (I2R,1,FNUM,0.5,TCOS)
         IF (I2R .NE. 1) GO TO 104
         DO 103 I=1,MR
            B(I) = Q(I,1)
            Q(I,1) = Q(I,2)
  103    CONTINUE
         GO TO 112
  104    DO 105 I=1,MR
            B(I) = Q(I,1)+0.5*(Q(I,JP2)-Q(I,JP1)-Q(I,JP3))
            Q(I,1) = Q(I,JP2)+Q(I,1)-Q(I,JP1)
  105    CONTINUE
         GO TO 112
  106    CONTINUE
         GO TO (107,108),IJUMP
  107    CONTINUE
         IJUMP = 2
         CALL COSGEN (I2R,1,0.5,0.0,TCOS)
  108    CONTINUE
         IF (I2R .NE. 1) GO TO 110
         DO 109 I=1,MR
            B(I) = 2.*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  109    CONTINUE
         GO TO 112
  110    DO 111 I=1,MR
            FI = Q(I,J)
            Q(I,J) = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = FI+Q(I,J)-Q(I,JM3)-Q(I,JP3)
  111    CONTINUE
  112    CONTINUE
         CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
         DO 113 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  113    CONTINUE
C
C     END OF REDUCTION FOR REGULAR UNKNOWNS.
C
  114 CONTINUE
C
C     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
C
      J = JSTOP+JR
  115 NLAST = J
      JM1 = J-I2RBY2
      JM2 = J-I2R
      JM3 = JM2-I2RBY2
      IF (NROD .EQ. 0) GO TO 125
C
C     ODD NUMBER OF UNKNOWNS
C
      IF (I2R .NE. 1) GO TO 117
      DO 116 I=1,MR
         B(I) = Q(I,J)
         Q(I,J) = Q(I,JM2)
  116 CONTINUE
      GO TO 123
  117 DO 118 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  118 CONTINUE
      IF (NRODPR .NE. 0) GO TO 120
      DO 119 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)
  119 CONTINUE
      IP = IP-MR
      GO TO 122
  120 CONTINUE
      DO 121 I=1,MR
         Q(I,J) = Q(I,J)-Q(I,JM1)+Q(I,JM2)
  121 CONTINUE
  122 IF (LR .EQ. 0) GO TO 123
      CALL COSGEN (LR,1,FNUM2,0.5,TCOS(KR+1))
  123 CONTINUE
      CALL COSGEN (KR,1,FNUM2,0.5,TCOS)
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 124 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
  124 CONTINUE
      KR = KR+I2R
      GO TO 141
  125 CONTINUE
C
C     EVEN NUMBER OF UNKNOWNS
C
      JP1 = J+I2RBY2
      JP2 = J+I2R
      IF (I2R .NE. 1) GO TO 129
      DO 126 I=1,MR
         B(I) = Q(I,J)
  126 CONTINUE
      TCOS(1) = 0.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      IP = 0
      IPSTOR = MR
      DO 127 I=1,MR
         P(I) = B(I)
         B(I) = B(I)+Q(I,N)
  127 CONTINUE
      TCOS(1) = -1.+2.*FLOAT(NP/2)
      TCOS(2) = 0.
      CALL TRIX (1,1,MR,A,BB,C,B,TCOS,D,W)
      DO 128 I=1,MR
         Q(I,J) = Q(I,JM2)+P(I)+B(I)
  128 CONTINUE
      GO TO 140
  129 CONTINUE
      DO 130 I=1,MR
         B(I) = Q(I,J)+.5*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  130 CONTINUE
      IF (NRODPR .NE. 0) GO TO 132
      DO 131 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  131 CONTINUE
      GO TO 134
  132 CONTINUE
      DO 133 I=1,MR
         B(I) = B(I)+Q(I,JP2)-Q(I,JP1)
  133 CONTINUE
  134 CONTINUE
      CALL COSGEN (I2R,1,0.5,0.0,TCOS)
      CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
      IP = IP+MR
      IPSTOR = MAX0(IPSTOR,IP+MR)
      DO 135 I=1,MR
         II = IP+I
         P(II) = B(I)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         B(I) = P(II)+Q(I,JP2)
  135 CONTINUE
      IF (LR .EQ. 0) GO TO 136
      CALL COSGEN (LR,1,FNUM2,0.5,TCOS(I2R+1))
      CALL MERGE (TCOS,0,I2R,I2R,LR,KR)
      GO TO 138
  136 DO 137 I=1,I2R
         II = KR+I
         TCOS(II) = TCOS(I)
  137 CONTINUE
  138 CALL COSGEN (KR,1,FNUM2,0.5,TCOS)
      CALL TRIX (KR,KR,MR,A,BB,C,B,TCOS,D,W)
      DO 139 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)+B(I)
  139 CONTINUE
  140 CONTINUE
      LR = KR
      KR = KR+JR
  141 CONTINUE
      NR = (NLAST-1)/JR+1
      IF (NR .LE. 3) GO TO 142
      I2R = JR
      NRODPR = NROD
      GO TO 101
  142 CONTINUE
C
C      BEGIN SOLUTION
C
      J = 1+JR
      JM1 = J-I2R
      JP1 = J+I2R
      JM2 = NLAST-I2R
      IF (NR .EQ. 2) GO TO 180
      IF (LR .NE. 0) GO TO 167
      IF (N .NE. 3) GO TO 156
C
C     CASE N = 3.
C
      GO TO (143,148,143),NP
  143 DO 144 I=1,MR
         B(I) = Q(I,2)
         B2(I) = Q(I,1)+Q(I,3)
         B3(I) = 0.
  144 CONTINUE
      GO TO (146,146,145),NP
  145 TCOS(1) = -1.
      TCOS(2) = 1.
      K1 = 1
      GO TO 147
  146 TCOS(1) = -2.
      TCOS(2) = 1.
      TCOS(3) = -1.
      K1 = 2
  147 K2 = 1
      K3 = 0
      K4 = 0
      GO TO 150
  148 DO 149 I=1,MR
         B(I) = Q(I,2)
         B2(I) = Q(I,3)
         B3(I) = Q(I,1)
  149 CONTINUE
      CALL COSGEN (3,1,0.5,0.0,TCOS)
      TCOS(4) = -1.
      TCOS(5) = 1.
      TCOS(6) = -1.
      TCOS(7) = 1.
      K1 = 3
      K2 = 2
      K3 = 1
      K4 = 1
  150 CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 151 I=1,MR
         B(I) = B(I)+B2(I)+B3(I)
  151 CONTINUE
      GO TO (153,153,152),NP
  152 TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  153 DO 154 I=1,MR
         Q(I,2) = B(I)
         B(I) = Q(I,1)+B(I)
  154 CONTINUE
      TCOS(1) = -1.+4.*FNUM
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 155 I=1,MR
         Q(I,1) = B(I)
  155 CONTINUE
      JR = 1
      I2R = 0
      GO TO 188
C
C     CASE N = 2**P+1
C
  156 CONTINUE
      DO 157 I=1,MR
         B(I) = Q(I,J)+Q(I,1)-Q(I,JM1)+Q(I,NLAST)-Q(I,JM2)
  157 CONTINUE
      GO TO (158,160,158),NP
  158 DO 159 I=1,MR
         B2(I) = Q(I,1)+Q(I,NLAST)+Q(I,J)-Q(I,JM1)-Q(I,JP1)
         B3(I) = 0.
  159 CONTINUE
      K1 = NLAST-1
      K2 = NLAST+JR-1
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(NLAST))
      TCOS(K2) = 2.*FLOAT(NP-2)
      CALL COSGEN (JR,1,0.5-FNUM,0.5,TCOS(K2+1))
      K3 = (3-NP)/2
      CALL MERGE (TCOS,K1,JR-K3,K2-K3,JR+K3,0)
      K1 = K1-1+K3
      CALL COSGEN (JR,1,FNUM,0.5,TCOS(K1+1))
      K2 = JR
      K3 = 0
      K4 = 0
      GO TO 162
  160 DO 161 I=1,MR
         FI = (Q(I,J)-Q(I,JM1)-Q(I,JP1))/2.
         B2(I) = Q(I,1)+FI
         B3(I) = Q(I,NLAST)+FI
  161 CONTINUE
      K1 = NLAST+JR-1
      K2 = K1+JR-1
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K1+1))
      CALL COSGEN (NLAST,1,0.5,0.0,TCOS(K2+1))
      CALL MERGE (TCOS,K1,JR-1,K2,NLAST,0)
      K3 = K1+NLAST-1
      K4 = K3+JR
      CALL COSGEN (JR,1,0.5,0.5,TCOS(K3+1))
      CALL COSGEN (JR,1,0.0,0.5,TCOS(K4+1))
      CALL MERGE (TCOS,K3,JR,K4,JR,K1)
      K2 = NLAST-1
      K3 = JR
      K4 = JR
  162 CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 163 I=1,MR
         B(I) = B(I)+B2(I)+B3(I)
  163 CONTINUE
      IF (NP .NE. 3) GO TO 164
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  164 DO 165 I=1,MR
         Q(I,J) = B(I)+.5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         B(I) = Q(I,J)+Q(I,1)
  165 CONTINUE
      CALL COSGEN (JR,1,FNUM,0.5,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 166 I=1,MR
         Q(I,1) = Q(I,1)-Q(I,JM1)+B(I)
  166 CONTINUE
      GO TO 188
C
C     CASE OF GENERAL N WITH NR = 3 .
C
  167 CONTINUE
      DO 168 I=1,MR
         B(I) = Q(I,1)-Q(I,JM1)+Q(I,J)
  168 CONTINUE
      IF (NROD .NE. 0) GO TO 170
      DO 169 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  169 CONTINUE
      GO TO 172
  170 DO 171 I=1,MR
         B(I) = B(I)+Q(I,NLAST)-Q(I,JM2)
  171 CONTINUE
  172 CONTINUE
      DO 173 I=1,MR
         T = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         Q(I,J) = T
         B2(I) = Q(I,NLAST)+T
         B3(I) = Q(I,1)+T
  173 CONTINUE
      K1 = KR+2*JR
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K1+1))
      K2 = K1+JR
      TCOS(K2) = 2.*FLOAT(NP-2)
      K4 = (NP-1)*(3-NP)
      K3 = K2+1-K4
      CALL COSGEN (KR+JR+K4,1,FLOAT(K4)/2.,1.-FLOAT(K4),TCOS(K3))
      K4 = 1-NP/3
      CALL MERGE (TCOS,K1,JR-K4,K2-K4,KR+JR+K4,0)
      IF (NP .EQ. 3) K1 = K1-1
      K2 = KR+JR
      K4 = K1+K2
      CALL COSGEN (KR,1,FNUM2,0.5,TCOS(K4+1))
      K3 = K4+KR
      CALL COSGEN (JR,1,FNUM,0.5,TCOS(K3+1))
      CALL MERGE (TCOS,K4,KR,K3,JR,K1)
      K4 = K3+JR
      CALL COSGEN (LR,1,FNUM2,0.5,TCOS(K4+1))
      CALL MERGE (TCOS,K3,JR,K4,LR,K1+K2)
      CALL COSGEN (KR,1,FNUM2,0.5,TCOS(K3+1))
      K3 = KR
      K4 = KR
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 174 I=1,MR
         B(I) = B(I)+B2(I)+B3(I)
  174 CONTINUE
      IF (NP .NE. 3) GO TO 175
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  175 DO 176 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+Q(I,J)
  176 CONTINUE
      CALL COSGEN (JR,1,FNUM,0.5,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      IF (JR .NE. 1) GO TO 178
      DO 177 I=1,MR
         Q(I,1) = B(I)
  177 CONTINUE
      GO TO 188
  178 CONTINUE
      DO 179 I=1,MR
         Q(I,1) = Q(I,1)-Q(I,JM1)+B(I)
  179 CONTINUE
      GO TO 188
  180 CONTINUE
C
C     CASE OF GENERAL N AND NR = 2 .
C
      DO 181 I=1,MR
         II = IP+I
         B3(I) = 0.
         B(I) = Q(I,1)+P(II)
         Q(I,1) = Q(I,1)-Q(I,JM1)
         B2(I) = Q(I,1)+Q(I,NLAST)
  181 CONTINUE
      K1 = KR+JR
      K2 = K1+JR
      CALL COSGEN (JR-1,1,0.0,1.0,TCOS(K1+1))
      GO TO (182,183,182),NP
  182 TCOS(K2) = 2.*FLOAT(NP-2)
      CALL COSGEN (KR,1,0.0,1.0,TCOS(K2+1))
      GO TO 184
  183 CALL COSGEN (KR+1,1,0.5,0.0,TCOS(K2))
  184 K4 = 1-NP/3
      CALL MERGE (TCOS,K1,JR-K4,K2-K4,KR+K4,0)
      IF (NP .EQ. 3) K1 = K1-1
      K2 = KR
      CALL COSGEN (KR,1,FNUM2,0.5,TCOS(K1+1))
      K4 = K1+KR
      CALL COSGEN (LR,1,FNUM2,0.5,TCOS(K4+1))
      K3 = LR
      K4 = 0
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 185 I=1,MR
         B(I) = B(I)+B2(I)
  185 CONTINUE
      IF (NP .NE. 3) GO TO 186
      TCOS(1) = 2.
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
  186 DO 187 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  187 CONTINUE
  188 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      J = NLAST-JR
      DO 189 I=1,MR
         B(I) = Q(I,NLAST)+Q(I,J)
  189 CONTINUE
      JM2 = NLAST-I2R
      IF (JR .NE. 1) GO TO 191
      DO 190 I=1,MR
         Q(I,NLAST) = 0.
  190 CONTINUE
      GO TO 195
  191 CONTINUE
      IF (NROD .NE. 0) GO TO 193
      DO 192 I=1,MR
         II = IP+I
         Q(I,NLAST) = P(II)
  192 CONTINUE
      IP = IP-MR
      GO TO 195
  193 DO 194 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)-Q(I,JM2)
  194 CONTINUE
  195 CONTINUE
      CALL COSGEN (KR,1,FNUM2,0.5,TCOS)
      CALL COSGEN (LR,1,FNUM2,0.5,TCOS(KR+1))
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 196 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)+B(I)
  196 CONTINUE
      NLASTP = NLAST
  197 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR .EQ. 0) GO TO 210
      JSTART = 1+JR
      KR = KR-JR
      IF (NLAST+JR .GT. N) GO TO 198
      KR = KR-JR
      NLAST = NLAST+JR
      JSTOP = NLAST-JSTEP
      GO TO 199
  198 CONTINUE
      JSTOP = NLAST-JR
  199 CONTINUE
      LR = KR-JR
      CALL COSGEN (JR,1,0.5,0.0,TCOS)
      DO 209 J=JSTART,JSTOP,JSTEP
         JM2 = J-JR
         JP2 = J+JR
         IF (J .NE. JR) GO TO 201
         DO 200 I=1,MR
            B(I) = Q(I,J)+Q(I,JP2)
  200    CONTINUE
         GO TO 203
  201    CONTINUE
         DO 202 I=1,MR
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  202    CONTINUE
  203    CONTINUE
         IF (JR .NE. 1) GO TO 205
         DO 204 I=1,MR
            Q(I,J) = 0.
  204    CONTINUE
         GO TO 207
  205    CONTINUE
         JM1 = J-I2R
         JP1 = J+I2R
         DO 206 I=1,MR
            Q(I,J) = .5*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  206    CONTINUE
  207    CONTINUE
         CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
         DO 208 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  208    CONTINUE
  209 CONTINUE
      NROD = 1
      IF (NLAST+I2R .LE. N) NROD = 0
      IF (NLASTP .NE. NLAST) GO TO 188
      GO TO 197
  210 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
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
