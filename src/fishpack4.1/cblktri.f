C
C     file cblktri.f
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
      SUBROUTINE CBLKTR (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,      
     +                   IERROR,W)
C
C                                                                       
C DIMENSION OF           AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N),
C ARGUMENTS              W(SEE ARGUMENT LIST)                           
C                                                                       
C LATEST REVISION        NOVEMBER 1988
C                                                                       
C PURPOSE                CBLKTR SOLVES A SYSTEM OF LINEAR EQUATIONS     
C                        OF THE FORM                                    
C                                                                       
C                        AN(J)*X(I,J-1) + AM(I)*X(I-1,J) +              
C                        (BN(J)+BM(I))*X(I,J) + CN(J)*X(I,J+1) +        
C                        CM(I)*X(I+1,J) = Y(I,J)                        
C                                                                       
C                        FOR I = 1,2,...,M  AND  J = 1,2,...,N.         
C                                                                       
C                        I+1 AND I-1 ARE EVALUATED MODULO M AND         
C                        J+1 AND J-1 MODULO N, I.E.,                    
C                                                                       
C                        X(I,0) = X(I,N),  X(I,N+1) = X(I,1),           
C                        X(0,J) = X(M,J),  X(M+1,J) = X(1,J).           
C                                                                       
C                        THESE EQUATIONS USUALLY RESULT FROM THE        
C                        DISCRETIZATION OF SEPARABLE ELLIPTIC           
C                        EQUATIONS.  BOUNDARY CONDITIONS MAY BE         
C                        DIRICHLET, NEUMANN, OR PERIODIC.               
C                                                                       
C                        CBLKTRI IS A COMPLEX VERSION OF PACKAGE        
C                        BLKTRI ON ULIB.                                
C                                                                       
C USAGE                  CALL CBLKTR (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,    
C                                     CM,IDIMY,Y,IERROR,W)              
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON INPUT               IFLG                                           
C                                                                       
C                          = 0  INITIALIZATION ONLY.                    
C                               CERTAIN QUANTITIES THAT DEPEND ON NP,   
C                               N, AN, BN, AND CN ARE COMPUTED AND      
C                               STORED IN THE WORK ARRAY W.             
C                                                                       
C                          = 1  THE QUANTITIES THAT WERE COMPUTED       
C                               IN THE INITIALIZATION ARE USED          
C                               TO OBTAIN THE SOLUTION X(I,J).          
C                                                                       
C                               NOTE:                                   
C                               A CALL WITH IFLG=0 TAKES                
C                               APPROXIMATELY ONE HALF THE TIME         
C                               AS A CALL WITH IFLG = 1.                
C                               HOWEVER, THE INITIALIZATION DOES        
C                               NOT HAVE TO BE REPEATED UNLESS NP,      
C                               N, AN, BN, OR CN CHANGE.                
C                                                                       
C                        NP                                             
C                          = 0  IF AN(1) AND CN(N) ARE NOT ZERO,        
C                               WHICH CORRESPONDS TO PERIODIC           
C                               BOUNARY CONDITIONS.                     
C                                                                       
C                          = 1  IF AN(1) AND CN(N) ARE ZERO.            
C                                                                       
C                        N                                              
C                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.   
C                          N MUST BE GREATER THAN 4.                    
C                          THE OPERATION COUNT IS PROPORTIONAL TO       
C                          MNLOG2(N), HENCE N SHOULD BE SELECTED        
C                          LESS THAN OR EQUAL TO M.                     
C                                                                       
C                        AN,BN,CN                                       
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH N           
C                          THAT SPECIFY THE COEFFICIENTS IN THE         
C                          LINEAR EQUATIONS GIVEN ABOVE.                
C                                                                       
C                        MP                                             
C                          = 0  IF AM(1) AND CM(M) ARE NOT ZERO,        
C                               WHICH CORRESPONDS TO PERIODIC           
C                               BOUNDARY CONDITIONS.                    
C                                                                       
C                          = 1  IF AM(1) = CM(M) = 0  .                 
C                                                                       
C                        M                                              
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.   
C                           M MUST BE GREATER THAN 4.                   
C                                                                       
C                        AM,BM,CM                                       
C                          COMPLEX ONE-DIMENSIONAL ARRAYS OF LENGTH M   
C                          THAT SPECIFY THE COEFFICIENTS IN THE LINEAR  
C                          EQUATIONS GIVEN ABOVE.                       
C                                                                       
C                        IDIMY                                          
C                          THE ROW (OR FIRST) DIMENSION OF THE          
C                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS        
C                          IN THE PROGRAM CALLING CBLKTR.               
C                          THIS PARAMETER IS USED TO SPECIFY THE        
C                          VARIABLE DIMENSION OF Y.                     
C                          IDIMY MUST BE AT LEAST M.                    
C                                                                       
C                        Y                                              
C                          A COMPLEX TWO-DIMENSIONAL ARRAY THAT         
C                          SPECIFIES THE VALUES OF THE RIGHT SIDE OF    
C                          THE LINEAR SYSTEM OF EQUATIONS GIVEN ABOVE.  
C                          Y MUST BE DIMENSIONED Y(IDIMY,N) WITH        
C                          IDIMY .GE. M.                                
C                                                                       
C                        W                                              
C                          A ONE-DIMENSIONAL ARRAY THAT MUST BE         
C                          PROVIDED BY THE USER FOR WORK SPACE.         
C                                                                       
C                          IF NP=1 DEFINE K=INT(LOG2(N))+1 AND          
C                          SET L=2**(K+1) THEN W MUST HAVE DIMENSION    
C                          (K-2)*L+K+5+MAX(2N,12M)                      
C                                                                       
C                          IF NP=0 DEFINE K=INT(LOG2(N-1))+1 AND        
C                          SET L=2**(K+1) THEN W MUST HAVE DIMENSION    
C                          (K-2)*L+K+5+2N+MAX(2N,12M)                   
C                                                                       
C                          **IMPORTANT**                                
C                          FOR PURPOSES OF CHECKING, THE REQUIRED       
C                          DIMENSION OF W IS COMPUTED BY CBLKTR AND     
C                          STORED IN W(1) IN FLOATING POINT FORMAT.     
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON OUTPUT              Y                                              
C                          CONTAINS THE SOLUTION X.                     
C                                                                       
C                        IERROR                                         
C                          AN ERROR FLAG THAT INDICATES INVALID         
C                          INPUT PARAMETERS.  EXCEPT FOR NUMBER ZER0,   
C                          A SOLUTION IS NOT ATTEMPTED.                 
C                                                                       
C                        = 0  NO ERROR.                                 
C                        = 1  M IS LESS THAN 5                          
C                        = 2  N IS LESS THAN 5                          
C                        = 3  IDIMY IS LESS THAN M.                     
C                        = 4  CBLKTR FAILED WHILE COMPUTING RESULTS     
C                             THAT DEPEND ON THE COEFFICIENT ARRAYS     
C                             AN, BN, CN.  CHECK THESE ARRAYS.          
C                        = 5  AN(J)*CN(J-1) IS LESS THAN 0 FOR SOME J.  
C                                                                       
C                             POSSIBLE REASONS FOR THIS CONDITION ARE   
C                             1. THE ARRAYS AN AND CN ARE NOT CORRECT   
C                             2. TOO LARGE A GRID SPACING WAS USED      
C                                IN THE DISCRETIZATION OF THE ELLIPTIC  
C                                EQUATION.                              
C                             3. THE LINEAR EQUATIONS RESULTED FROM A   
C                                PARTIAL DIFFERENTIAL EQUATION WHICH    
C                                WAS NOT ELLIPTIC.                      
C                                                                       
C                        W                                              
C                           CONTAINS INTERMEDIATE VALUES THAT MUST      
C                           NOT BE DESTROYED IF CBLKTR WILL BE CALLED   
C                           AGAIN WITH IFLG=1. W(1) CONTAINS THE        
C                           NUMBER OF LOCATIONS REQUIRED BY W IN        
C                           FLOATING POINT FORMAT.                      
C                                                                       
C                                                                       
C SPECIAL CONDITIONS     THE ALGORITHM MAY FAIL IF ABS(BM(I)+BN(J))     
C                        IS LESS THAN ABS(AM(I))+ABS(AN(J))+            
C                        ABS(CM(I))+ABS(CN(J))                          
C                        FOR SOME I AND J. THE ALGORITHM WILL ALSO      
C                        FAIL IF AN(J)*CN(J-1) IS LESS THAN ZERO FOR    
C                        SOME J.                                        
C                        SEE THE DESCRIPTION OF THE OUTPUT PARAMETER    
C                        IERROR.                                        
C                                                                       
C                                                                       
C I/O                    NONE                                           
C                                                                       
C PRECISION              SINGLE                                         
C                                                                       
C REQUIRED LIBRARY       COMF FROM FISHPACK
C FILES
C                                                                       
C LANGUAGE               FORTRAN                                        
C                                                                       
C HISTORY                WRITTEN BY PAUL SWARZTRAUBER AT NCAR IN        
C                        THE EARLY 1970'S.  REWRITTEN AN RELEASED       
C                        ON NCAR'S PUBLIC SOFTWARE LIBRARIES IN         
C                        JANUARY, 1980.                                 
C                                                                       
C ALGORITHM              GENERALIZED CYCLIC REDUCTION                   
C                        (SEE REFERENCE BELOW)                          
C                                                                       
C PORTABILITY            FORTRAN 77
C                        THE APPROXIMATE MACHINE ACCURACY IS COMPUTED   
C                        IN FUNCTION EPMACH                             
C                                                                       
C REFERENCES             SWARZTRAUBER,P. AND R. SWEET, 'EFFICIENT       
C                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF        
C                        ELLIPTIC EQUATIONS'                            
C                        NCAR TN/IA-109, JULY, 1975, 138 PP.            
C                                                                       
C                        SWARZTRAUBER P. N.,A DIRECT METHOD FOR         
C                        THE DISCRETE SOLUTION OF SEPARABLE             
C                        ELLIPTIC EQUATIONS, S.I.A.M.                   
C                        J. NUMER. ANAL.,11(1974) PP. 1136-1150.        
C                                                                       
C***********************************************************************
      DIMENSION       AN(*)      ,BN(*)      ,CN(*)      ,AM(*)      ,  
     1                BM(*)      ,CM(*)      ,Y(IDIMY,1) ,W(*)          
      EXTERNAL        PROC       ,PROCP      ,CPROC      ,CPROCP        
      COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        ,  
     1                NM         ,NCMPLX     ,IK                        
      COMPLEX         AM         ,BM         ,CM         ,Y             
C                                                                       
C TEST M AND N FOR THE PROPER FORM                                      
C                                                                       
      NM = N                                                            
      M2 = M+M                                                          
      IERROR = 0                                                        
      IF (M-5) 101,102,102                                              
  101 IERROR = 1                                                        
      GO TO 119                                                         
  102 IF (NM-3) 103,104,104                                             
  103 IERROR = 2                                                        
      GO TO 119                                                         
  104 IF (IDIMY-M) 105,106,106                                          
  105 IERROR = 3                                                        
      GO TO 119                                                         
  106 NH = N                                                            
      NPP = NP                                                          
      IF (NPP) 107,108,107                                              
  107 NH = NH+1                                                         
  108 IK = 2                                                            
      K = 1                                                             
  109 IK = IK+IK                                                        
      K = K+1                                                           
      IF (NH-IK) 110,110,109                                            
  110 NL = IK                                                           
      IK = IK+IK                                                        
      NL = NL-1                                                         
      IWAH = (K-2)*IK+K+6                                               
      IF (NPP) 111,112,111                                              
C                                                                       
C     DIVIDE W INTO WORKING SUB ARRAYS                                  
C                                                                       
  111 IW1 = IWAH                                                        
      IWBH = IW1+NM                                                     
      W(1) = FLOAT(IW1-1+MAX0(2*NM,12*M))                               
      GO TO 113                                                         
  112 IWBH = IWAH+NM+NM                                                 
      IW1 = IWBH                                                        
      W(1) = FLOAT(IW1-1+MAX0(2*NM,12*M))                               
      NM = NM-1                                                         
C                                                                       
C SUBROUTINE COMP B COMPUTES THE ROOTS OF THE B POLYNOMIALS             
C                                                                       
  113 IF (IERROR) 119,114,119                                           
  114 IW2 = IW1+M2                                                      
      IW3 = IW2+M2                                                      
      IWD = IW3+M2                                                      
      IWW = IWD+M2                                                      
      IWU = IWW+M2                                                      
      IF (IFLG) 116,115,116                                             
  115 CALL CCOMPB (NL,IERROR,AN,BN,CN,W(2),W(IWAH),W(IWBH))              
      GO TO 119                                                         
  116 IF (MP) 117,118,117                                               
C                                                                       
C SUBROUTINE BLKTR1 SOLVES THE LINEAR SYSTEM                            
C                                                                       
  117 CALL CBLKT1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),W(IW1),W(IW2),   
     1             W(IW3),W(IWD),W(IWW),W(IWU),PROC,CPROC)              
      GO TO 119                                                         
  118 CALL CBLKT1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),W(IW1),W(IW2),   
     1             W(IW3),W(IWD),W(IWW),W(IWU),PROCP,CPROCP)            
  119 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE CBLKT1 (N,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,B,W1,W2,W3,WD,   
     1                   WW,WU,PRDCT,CPRDCT)                            
C                                                                       
C CBLKT1 SOLVES THE LINEAR SYSTEM                                       
C                                                                       
C B  CONTAINS THE ROOTS OF ALL THE B POLYNOMIALS                        
C W1,W2,W3,WD,WW,WU  ARE ALL WORKING ARRAYS                             
C PRDCT IS EITHER PROCP OR PROC DEPENDING ON WHETHER THE BOUNDARY       
C CONDITIONS IN THE M DIRECTION ARE PERIODIC OR NOT                     
C CPRDCT IS EITHER CPROCP OR CPROC WHICH ARE CALLED IF SOME OF THE ZEROS
C OF THE B POLYNOMIALS ARE COMPLEX.                                     
C                                                                       
C                                                                       
      DIMENSION       AN(*)      ,BN(*)      ,CN(*)      ,AM(*)      ,  
     1                BM(*)      ,CM(*)      ,B(*)       ,W1(*)      ,  
     2                W2(*)      ,W3(*)      ,WD(*)      ,WW(*)      ,  
     3                WU(*)      ,Y(IDIMY,1)                            
      COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        ,  
     1                NM         ,NCMPLX     ,IK                        
      COMPLEX         AM         ,BM         ,CM         ,Y          ,  
     1                W1         ,W2         ,W3         ,WD         ,  
     2                WW         ,WU                                    
C                                                                       
C BEGIN REDUCTION PHASE                                                 
C                                                                       
      KDO = K-1                                                         
      DO 109 L=1,KDO                                                    
         IR = L-1                                                       
         I2 = 2**IR                                                     
         I1 = I2/2                                                      
         I3 = I2+I1                                                     
         I4 = I2+I2                                                     
         IRM1 = IR-1                                                    
         CALL CINDXB (I2,IR,IM2,NM2)                                     
         CALL CINDXB (I1,IRM1,IM3,NM3)                                   
         CALL CINDXB (I3,IRM1,IM1,NM1)                                   
         CALL PRDCT (NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,Y(1,I2),W3, 
     1               M,AM,BM,CM,WD,WW,WU)                               
         IF = 2**K                                                      
         DO 108 I=I4,IF,I4                                              
            IF (I-NM) 101,101,108                                       
  101       IPI1 = I+I1                                                 
            IPI2 = I+I2                                                 
            IPI3 = I+I3                                                 
            CALL CINDXC (I,IR,IDXC,NC)                                   
            IF (I-IF) 102,108,108                                       
  102       CALL CINDXA (I,IR,IDXA,NA)                                   
            CALL CINDXB (I-I1,IRM1,IM1,NM1)                              
            CALL CINDXB (IPI2,IR,IP2,NP2)                                
            CALL CINDXB (IPI1,IRM1,IP1,NP1)                              
            CALL CINDXB (IPI3,IRM1,IP3,NP3)                              
            CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W3,W1,M,AM,  
     1                  BM,CM,WD,WW,WU)                                 
            IF (IPI2-NM) 105,105,103                                    
  103       DO 104 J=1,M                                                
               W3(J) = (0.,0.)                                          
               W2(J) = (0.,0.)                                          
  104       CONTINUE                                                    
            GO TO 106                                                   
  105       CALL PRDCT (NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,         
     1                  Y(1,IPI2),W3,M,AM,BM,CM,WD,WW,WU)               
            CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W3,W2,M,AM,  
     1                  BM,CM,WD,WW,WU)                                 
  106       DO 107 J=1,M                                                
               Y(J,I) = W1(J)+W2(J)+Y(J,I)                              
  107       CONTINUE                                                    
  108    CONTINUE                                                       
  109 CONTINUE                                                          
      IF (NPP) 132,110,132                                              
C                                                                       
C     THE PERIODIC CASE IS TREATED USING THE CAPACITANCE MATRIX METHOD  
C                                                                       
  110 IF = 2**K                                                         
      I = IF/2                                                          
      I1 = I/2                                                          
      CALL CINDXB (I-I1,K-2,IM1,NM1)                                     
      CALL CINDXB (I+I1,K-2,IP1,NP1)                                     
      CALL CINDXB (I,K-1,IZ,NZ)                                          
      CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,Y(1,I),W1,M,AM,  
     1            BM,CM,WD,WW,WU)                                       
      IZR = I                                                           
      DO 111 J=1,M                                                      
         W2(J) = W1(J)                                                  
  111 CONTINUE                                                          
      DO 113 LL=2,K                                                     
         L = K-LL+1                                                     
         IR = L-1                                                       
         I2 = 2**IR                                                     
         I1 = I2/2                                                      
         I = I2                                                         
         CALL CINDXC (I,IR,IDXC,NC)                                      
         CALL CINDXB (I,IR,IZ,NZ)                                        
         CALL CINDXB (I-I1,IR-1,IM1,NM1)                                 
         CALL CINDXB (I+I1,IR-1,IP1,NP1)                                 
         CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W1,W1,M,AM,BM,  
     1               CM,WD,WW,WU)                                       
         DO 112 J=1,M                                                   
            W1(J) = Y(J,I)+W1(J)                                        
  112    CONTINUE                                                       
         CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,W1,M,AM,   
     1               BM,CM,WD,WW,WU)                                    
  113 CONTINUE                                                          
      DO 118 LL=2,K                                                     
         L = K-LL+1                                                     
         IR = L-1                                                       
         I2 = 2**IR                                                     
         I1 = I2/2                                                      
         I4 = I2+I2                                                     
         IFD = IF-I2                                                    
         DO 117 I=I2,IFD,I4                                             
            IF (I-I2-IZR) 117,114,117                                   
  114       IF (I-NM) 115,115,118                                       
  115       CALL CINDXA (I,IR,IDXA,NA)                                   
            CALL CINDXB (I,IR,IZ,NZ)                                     
            CALL CINDXB (I-I1,IR-1,IM1,NM1)                              
            CALL CINDXB (I+I1,IR-1,IP1,NP1)                              
            CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W2,W2,M,AM,  
     1                  BM,CM,WD,WW,WU)                                 
            DO 116 J=1,M                                                
               W2(J) = Y(J,I)+W2(J)                                     
  116       CONTINUE                                                    
            CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W2,W2,M,   
     1                  AM,BM,CM,WD,WW,WU)                              
            IZR = I                                                     
            IF (I-NM) 117,119,117                                       
  117    CONTINUE                                                       
  118 CONTINUE                                                          
  119 DO 120 J=1,M                                                      
         Y(J,NM+1) = Y(J,NM+1)-CN(NM+1)*W1(J)-AN(NM+1)*W2(J)            
  120 CONTINUE                                                          
      CALL CINDXB (IF/2,K-1,IM1,NM1)                                     
      CALL CINDXB (IF,K-1,IP,NP)                                         
      IF (NCMPLX) 121,122,121                                           
  121 CALL CPRDCT (NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1),         
     1             Y(1,NM+1),M,AM,BM,CM,W1,W3,WW)                       
      GO TO 123                                                         
  122 CALL PRDCT (NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1),          
     1            Y(1,NM+1),M,AM,BM,CM,WD,WW,WU)                        
  123 DO 124 J=1,M                                                      
         W1(J) = AN(1)*Y(J,NM+1)                                        
         W2(J) = CN(NM)*Y(J,NM+1)                                       
         Y(J,1) = Y(J,1)-W1(J)                                          
         Y(J,NM) = Y(J,NM)-W2(J)                                        
  124 CONTINUE                                                          
      DO 126 L=1,KDO                                                    
         IR = L-1                                                       
         I2 = 2**IR                                                     
         I4 = I2+I2                                                     
         I1 = I2/2                                                      
         I = I4                                                         
         CALL CINDXA (I,IR,IDXA,NA)                                      
         CALL CINDXB (I-I2,IR,IM2,NM2)                                   
         CALL CINDXB (I-I2-I1,IR-1,IM3,NM3)                              
         CALL CINDXB (I-I1,IR-1,IM1,NM1)                                 
         CALL PRDCT (NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,W1,W1,M,AM, 
     1               BM,CM,WD,WW,WU)                                    
         CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W1,W1,M,AM,BM,  
     1               CM,WD,WW,WU)                                       
         DO 125 J=1,M                                                   
            Y(J,I) = Y(J,I)-W1(J)                                       
  125    CONTINUE                                                       
  126 CONTINUE                                                          
C                                                                       
      IZR = NM                                                          
      DO 131 L=1,KDO                                                    
         IR = L-1                                                       
         I2 = 2**IR                                                     
         I1 = I2/2                                                      
         I3 = I2+I1                                                     
         I4 = I2+I2                                                     
         IRM1 = IR-1                                                    
         DO 130 I=I4,IF,I4                                              
            IPI1 = I+I1                                                 
            IPI2 = I+I2                                                 
            IPI3 = I+I3                                                 
            IF (IPI2-IZR) 127,128,127                                   
  127       IF (I-IZR) 130,131,130                                      
  128       CALL CINDXC (I,IR,IDXC,NC)                                   
            CALL CINDXB (IPI2,IR,IP2,NP2)                                
            CALL CINDXB (IPI1,IRM1,IP1,NP1)                              
            CALL CINDXB (IPI3,IRM1,IP3,NP3)                              
            CALL PRDCT (NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,W2,W2,M, 
     1                  AM,BM,CM,WD,WW,WU)                              
            CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W2,W2,M,AM,  
     1                  BM,CM,WD,WW,WU)                                 
            DO 129 J=1,M                                                
               Y(J,I) = Y(J,I)-W2(J)                                    
  129       CONTINUE                                                    
            IZR = I                                                     
            GO TO 131                                                   
  130    CONTINUE                                                       
  131 CONTINUE                                                          
C                                                                       
C BEGIN BACK SUBSTITUTION PHASE                                         
C                                                                       
  132 DO 144 LL=1,K                                                     
         L = K-LL+1                                                     
         IR = L-1                                                       
         IRM1 = IR-1                                                    
         I2 = 2**IR                                                     
         I1 = I2/2                                                      
         I4 = I2+I2                                                     
         IFD = IF-I2                                                    
         DO 143 I=I2,IFD,I4                                             
            IF (I-NM) 133,133,143                                       
  133       IMI1 = I-I1                                                 
            IMI2 = I-I2                                                 
            IPI1 = I+I1                                                 
            IPI2 = I+I2                                                 
            CALL CINDXA (I,IR,IDXA,NA)                                   
            CALL CINDXC (I,IR,IDXC,NC)                                   
            CALL CINDXB (I,IR,IZ,NZ)                                     
            CALL CINDXB (IMI1,IRM1,IM1,NM1)                              
            CALL CINDXB (IPI1,IRM1,IP1,NP1)                              
            IF (I-I2) 134,134,136                                       
  134       DO 135 J=1,M                                                
               W1(J) = (0.,0.)                                          
  135       CONTINUE                                                    
            GO TO 137                                                   
  136       CALL PRDCT (NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),Y(1,IMI2),   
     1                  W1,M,AM,BM,CM,WD,WW,WU)                         
  137       IF (IPI2-NM) 140,140,138                                    
  138       DO 139 J=1,M                                                
               W2(J) = (0.,0.)                                          
  139       CONTINUE                                                    
            GO TO 141                                                   
  140       CALL PRDCT (NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),Y(1,IPI2),   
     1                  W2,M,AM,BM,CM,WD,WW,WU)                         
  141       DO 142 J=1,M                                                
               W1(J) = Y(J,I)+W1(J)+W2(J)                               
  142       CONTINUE                                                    
            CALL PRDCT (NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,Y(1,I), 
     1                  M,AM,BM,CM,WD,WW,WU)                            
  143    CONTINUE                                                       
  144 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      FUNCTION CBSRH (XLL,XRR,IZ,C,A,BH,F,SGN)                           
      DIMENSION       A(*)       ,C(*)       ,BH(*)                     
      COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        ,  
     1                NM         ,NCMPLX     ,IK                        
      XL = XLL                                                          
      XR = XRR                                                          
      DX = .5*ABS(XR-XL)                                                
  101 X = .5*(XL+XR)                                                    
      IF (SGN*F(X,IZ,C,A,BH)) 103,105,102                               
  102 XR = X                                                            
      GO TO 104                                                         
  103 XL = X                                                            
  104 DX = .5*DX                                                        
      IF (DX-CNV) 105,105,101                                           
  105 CBSRH = .5*(XL+XR)                                                 
      RETURN                                                            
      END                                                               
      SUBROUTINE CCOMPB (N,IERROR,AN,BN,CN,B,AH,BH)                      
C                                                                       
C     CCOMPB COMPUTES THE ROOTS OF THE B POLYNOMIALS USING SUBROUTINE    
C     CTEVLS WHICH IS A MODIFICATION THE EISPACK PROGRAM TQLRAT.         
C     IERROR IS SET TO 4 IF EITHER CTEVLS FAILS OR IF A(J+1)*C(J) IS     
C     LESS THAN ZERO FOR SOME J.  AH,BH ARE TEMPORARY WORK ARRAYS.      
C                                                                       
      DIMENSION       AN(*)      ,BN(*)      ,CN(*)      ,B(*)       ,  
     1                AH(*)      ,BH(*)                                 
      COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        ,  
     1                NM         ,NCMPLX     ,IK                        
      EPS = EPMACH(DUM)                                                 
      BNORM = ABS(BN(1))                                                
      DO 102 J=2,NM                                                     
         BNORM = AMAX1(BNORM,ABS(BN(J)))                                
         ARG = AN(J)*CN(J-1)                                            
         IF (ARG) 119,101,101                                           
  101    B(J) = SIGN(SQRT(ARG),AN(J))                                   
  102 CONTINUE                                                          
      CNV = EPS*BNORM                                                   
      IF = 2**K                                                         
      KDO = K-1                                                         
      DO 108 L=1,KDO                                                    
         IR = L-1                                                       
         I2 = 2**IR                                                     
         I4 = I2+I2                                                     
         IPL = I4-1                                                     
         IFD = IF-I4                                                    
         DO 107 I=I4,IFD,I4                                             
            CALL CINDXB (I,L,IB,NB)                                      
            IF (NB) 108,108,103                                         
  103       JS = I-IPL                                                  
            JF = JS+NB-1                                                
            LS = 0                                                      
            DO 104 J=JS,JF                                              
               LS = LS+1                                                
               BH(LS) = BN(J)                                           
               AH(LS) = B(J)                                            
  104       CONTINUE                                                    
            CALL CTEVLS (NB,BH,AH,IERROR)                                
            IF (IERROR) 118,105,118                                     
  105       LH = IB-1                                                   
            DO 106 J=1,NB                                               
               LH = LH+1                                                
               B(LH) = -BH(J)                                           
  106       CONTINUE                                                    
  107    CONTINUE                                                       
  108 CONTINUE                                                          
      DO 109 J=1,NM                                                     
         B(J) = -BN(J)                                                  
  109 CONTINUE                                                          
      IF (NPP) 117,110,117                                              
  110 NMP = NM+1                                                        
      NB = NM+NMP                                                       
      DO 112 J=1,NB                                                     
         L1 = MOD(J-1,NMP)+1                                            
         L2 = MOD(J+NM-1,NMP)+1                                         
         ARG = AN(L1)*CN(L2)                                            
         IF (ARG) 119,111,111                                           
  111    BH(J) = SIGN(SQRT(ARG),-AN(L1))                                
         AH(J) = -BN(L1)                                                
  112 CONTINUE                                                          
      CALL CTEVLS (NB,AH,BH,IERROR)                                      
      IF (IERROR) 118,113,118                                           
  113 CALL CINDXB (IF,K-1,J2,LH)                                         
      CALL CINDXB (IF/2,K-1,J1,LH)                                       
      J2 = J2+1                                                         
      LH = J2                                                           
      N2M2 = J2+NM+NM-2                                                 
  114 D1 = ABS(B(J1)-B(J2-1))                                           
      D2 = ABS(B(J1)-B(J2))                                             
      D3 = ABS(B(J1)-B(J2+1))                                           
      IF ((D2 .LT. D1) .AND. (D2 .LT. D3)) GO TO 115                    
      B(LH) = B(J2)                                                     
      J2 = J2+1                                                         
      LH = LH+1                                                         
      IF (J2-N2M2) 114,114,116                                          
  115 J2 = J2+1                                                         
      J1 = J1+1                                                         
      IF (J2-N2M2) 114,114,116                                          
  116 B(LH) = B(N2M2+1)                                                 
      CALL CINDXB (IF,K-1,J1,J2)                                         
      J2 = J1+NMP+NMP                                                   
      CALL CPPADD (NM+1,IERROR,AN,CN,B(J1),B(J1),B(J2))                  
  117 RETURN                                                            
  118 IERROR = 4                                                        
      RETURN                                                            
  119 IERROR = 5                                                        
      RETURN                                                            
      END                                                               
      SUBROUTINE CPROC (ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,W,YY) 
C                                                                       
C PROC APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND      
C STORES THE RESULT IN Y                                                
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X              
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY      
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS       
C NA IS THE LENGTH OF THE ARRAY AA                                      
C X,Y THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y        
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX                
C M  IS THE ORDER OF THE MATRIX                                         
C D,W ARE WORK ARRAYS                                                   
C ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE              
C                                                                       
      COMPLEX         Y          ,D          ,W          ,BD         ,  
     1                CRT        ,DEN        ,Y1         ,Y2         ,  
     2                X          ,A          ,B          ,C             
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,X(*)       ,  
     1                Y(*)       ,D(*)       ,W(*)       ,BD(*)      ,  
     2                BM1(*)     ,BM2(*)     ,AA(*)      ,YY(*)         
      DO 101 J=1,M                                                      
         Y(J) = X(J)                                                    
  101 CONTINUE                                                          
      MM = M-1                                                          
      ID = ND                                                           
      M1 = NM1                                                          
      M2 = NM2                                                          
      IA = NA                                                           
  102 IFLG = 0                                                          
      IF (ID) 109,109,103                                               
  103 CRT = BD(ID)                                                      
      ID = ID-1                                                         
C                                                                       
C BEGIN SOLUTION TO SYSTEM                                              
C                                                                       
      D(M) = A(M)/(B(M)-CRT)                                            
      W(M) = Y(M)/(B(M)-CRT)                                            
      DO 104 J=2,MM                                                     
         K = M-J                                                        
         DEN = B(K+1)-CRT-C(K+1)*D(K+2)                                 
         D(K+1) = A(K+1)/DEN                                            
         W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN                            
  104 CONTINUE                                                          
      DEN = B(1)-CRT-C(1)*D(2)                                          
      IF (CABS(DEN)) 105,106,105                                        
  105 Y(1) = (Y(1)-C(1)*W(2))/DEN                                       
      GO TO 107                                                         
  106 Y(1) = (1.,0.)                                                    
  107 DO 108 J=2,M                                                      
         Y(J) = W(J)-D(J)*Y(J-1)                                        
  108 CONTINUE                                                          
  109 IF (M1) 110,110,112                                               
  110 IF (M2) 121,121,111                                               
  111 RT = BM2(M2)                                                      
      M2 = M2-1                                                         
      GO TO 117                                                         
  112 IF (M2) 113,113,114                                               
  113 RT = BM1(M1)                                                      
      M1 = M1-1                                                         
      GO TO 117                                                         
  114 IF (ABS(BM1(M1))-ABS(BM2(M2))) 116,116,115                        
  115 RT = BM1(M1)                                                      
      M1 = M1-1                                                         
      GO TO 117                                                         
  116 RT = BM2(M2)                                                      
      M2 = M2-1                                                         
  117 Y1 = (B(1)-RT)*Y(1)+C(1)*Y(2)                                     
      IF (MM-2) 120,118,118                                             
C                                                                       
C MATRIX MULTIPLICATION                                                 
C                                                                       
  118 DO 119 J=2,MM                                                     
         Y2 = A(J)*Y(J-1)+(B(J)-RT)*Y(J)+C(J)*Y(J+1)                    
         Y(J-1) = Y1                                                    
         Y1 = Y2                                                        
  119 CONTINUE                                                          
  120 Y(M) = A(M)*Y(M-1)+(B(M)-RT)*Y(M)                                 
      Y(M-1) = Y1                                                       
      IFLG = 1                                                          
      GO TO 102                                                         
  121 IF (IA) 124,124,122                                               
  122 RT = AA(IA)                                                       
      IA = IA-1                                                         
      IFLG = 1                                                          
C                                                                       
C SCALAR MULTIPLICATION                                                 
C                                                                       
      DO 123 J=1,M                                                      
         Y(J) = RT*Y(J)                                                 
  123 CONTINUE                                                          
  124 IF (IFLG) 125,125,102                                             
  125 RETURN                                                            
      END                                                               
      SUBROUTINE CPROCP (ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,U,YY)
C                                                                       
C CPROCP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND    
C STORES THE RESULT IN Y                                                
C                                                                       
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS       
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY      
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X              
C NA IS THE LENGTH OF THE ARRAY AA                                      
C X,Y THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y        
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX                
C M  IS THE ORDER OF THE MATRIX                                         
C D,U ARE WORK ARRAYS                                                   
C ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE              
C                                                                       
      COMPLEX         Y          ,D          ,U          ,V          ,  
     1                DEN        ,BH         ,YM         ,AM         ,  
     2                Y1         ,Y2         ,YH         ,BD         ,  
     3                CRT        ,X          ,A          ,B          ,C 
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,X(*)       ,  
     1                Y(*)       ,D(*)       ,U(*)       ,BD(*)      ,  
     2                BM1(*)     ,BM2(*)     ,AA(*)      ,YY(*)         
      DO 101 J=1,M                                                      
         Y(J) = X(J)                                                    
  101 CONTINUE                                                          
      MM = M-1                                                          
      MM2 = M-2                                                         
      ID = ND                                                           
      M1 = NM1                                                          
      M2 = NM2                                                          
      IA = NA                                                           
  102 IFLG = 0                                                          
      IF (ID) 111,111,103                                               
  103 CRT = BD(ID)                                                      
      ID = ID-1                                                         
      IFLG = 1                                                          
C                                                                       
C BEGIN SOLUTION TO SYSTEM                                              
C                                                                       
      BH = B(M)-CRT                                                     
      YM = Y(M)                                                         
      DEN = B(1)-CRT                                                    
      D(1) = C(1)/DEN                                                   
      U(1) = A(1)/DEN                                                   
      Y(1) = Y(1)/DEN                                                   
      V = C(M)                                                          
      IF (MM2-2) 106,104,104                                            
  104 DO 105 J=2,MM2                                                    
         DEN = B(J)-CRT-A(J)*D(J-1)                                     
         D(J) = C(J)/DEN                                                
         U(J) = -A(J)*U(J-1)/DEN                                        
         Y(J) = (Y(J)-A(J)*Y(J-1))/DEN                                  
         BH = BH-V*U(J-1)                                               
         YM = YM-V*Y(J-1)                                               
         V = -V*D(J-1)                                                  
  105 CONTINUE                                                          
  106 DEN = B(M-1)-CRT-A(M-1)*D(M-2)                                    
      D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN                               
      Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2))/DEN                               
      AM = A(M)-V*D(M-2)                                                
      BH = BH-V*U(M-2)                                                  
      YM = YM-V*Y(M-2)                                                  
      DEN = BH-AM*D(M-1)                                                
      IF (CABS(DEN)) 107,108,107                                        
  107 Y(M) = (YM-AM*Y(M-1))/DEN                                         
      GO TO 109                                                         
  108 Y(M) = (1.,0.)                                                    
  109 Y(M-1) = Y(M-1)-D(M-1)*Y(M)                                       
      DO 110 J=2,MM                                                     
         K = M-J                                                        
         Y(K) = Y(K)-D(K)*Y(K+1)-U(K)*Y(M)                              
  110 CONTINUE                                                          
  111 IF (M1) 112,112,114                                               
  112 IF (M2) 123,123,113                                               
  113 RT = BM2(M2)                                                      
      M2 = M2-1                                                         
      GO TO 119                                                         
  114 IF (M2) 115,115,116                                               
  115 RT = BM1(M1)                                                      
      M1 = M1-1                                                         
      GO TO 119                                                         
  116 IF (ABS(BM1(M1))-ABS(BM2(M2))) 118,118,117                        
  117 RT = BM1(M1)                                                      
      M1 = M1-1                                                         
      GO TO 119                                                         
  118 RT = BM2(M2)                                                      
      M2 = M2-1                                                         
C                                                                       
C MATRIX MULTIPLICATION                                                 
C                                                                       
  119 YH = Y(1)                                                         
      Y1 = (B(1)-RT)*Y(1)+C(1)*Y(2)+A(1)*Y(M)                           
      IF (MM-2) 122,120,120                                             
  120 DO 121 J=2,MM                                                     
         Y2 = A(J)*Y(J-1)+(B(J)-RT)*Y(J)+C(J)*Y(J+1)                    
         Y(J-1) = Y1                                                    
         Y1 = Y2                                                        
  121 CONTINUE                                                          
  122 Y(M) = A(M)*Y(M-1)+(B(M)-RT)*Y(M)+C(M)*YH                         
      Y(M-1) = Y1                                                       
      IFLG = 1                                                          
      GO TO 102                                                         
  123 IF (IA) 126,126,124                                               
  124 RT = AA(IA)                                                       
      IA = IA-1                                                         
      IFLG = 1                                                          
C                                                                       
C SCALAR MULTIPLICATION                                                 
C                                                                       
      DO 125 J=1,M                                                      
         Y(J) = RT*Y(J)                                                 
  125 CONTINUE                                                          
  126 IF (IFLG) 127,127,102                                             
  127 RETURN                                                            
      END                                                               
      SUBROUTINE CINDXA (I,IR,IDXA,NA)                                   
      COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        ,  
     1                NM         ,NCMPLX     ,IK                        
      NA = 2**IR                                                        
      IDXA = I-NA+1                                                     
      IF (I-NM) 102,102,101                                             
  101 NA = 0                                                            
  102 RETURN                                                            
      END                                                               
      SUBROUTINE CINDXB (I,IR,IDX,IDP)                                   
C                                                                       
C B(IDX) IS THE LOCATION OF THE FIRST ROOT OF THE B(I,IR) POLYNOMIAL    
C                                                                       
      COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        ,  
     1                NM         ,NCMPLX     ,IK                        
      IDP = 0                                                           
      IF (IR) 107,101,103                                               
  101 IF (I-NM) 102,102,107                                             
  102 IDX = I                                                           
      IDP = 1                                                           
      RETURN                                                            
  103 IZH = 2**IR                                                       
      ID = I-IZH-IZH                                                    
      IDX = ID+ID+(IR-1)*IK+IR+(IK-I)/IZH+4                             
      IPL = IZH-1                                                       
      IDP = IZH+IZH-1                                                   
      IF (I-IPL-NM) 105,105,104                                         
  104 IDP = 0                                                           
      RETURN                                                            
  105 IF (I+IPL-NM) 107,107,106                                         
  106 IDP = NM+IPL-I+1                                                  
  107 RETURN                                                            
      END                                                               
      SUBROUTINE CINDXC (I,IR,IDXC,NC)                                   
      COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        ,  
     1                NM         ,NCMPLX     ,IK                        
      NC = 2**IR                                                        
      IDXC = I                                                          
      IF (IDXC+NC-1-NM) 102,102,101                                     
  101 NC = 0                                                            
  102 RETURN                                                            
      END                                                               
      SUBROUTINE CPPADD (N,IERROR,A,C,CBP,BP,BH)                         
C                                                                       
C     CPPADD COMPUTES THE EIGENVALUES OF THE PERIODIC TRIDIAGONAL 
C     MATRIX WITH COEFFICIENTS AN,BN,CN                                        
C                                                                       
C     N IS THE ORDER OF THE BH AND BP POLYNOMIALS                           
C     ON OUTPUT BP CONTAINS THE EIGENVALUES                             
C     CBP IS THE SAME AS BP EXCEPT TYPE COMPLEX                             
C     BH IS USED TO TEMPORARILY STORE THE ROOTS OF THE B HAT POLYNOMIAL     
C       WHICH ENTERS THROUGH BP                                               
C                                                                       
      COMPLEX         CF         ,CX         ,FSG        ,HSG        ,  
     1                DD         ,F          ,FP         ,FPP        ,  
     2                CDIS       ,R1         ,R2         ,R3         ,  
     3                CBP                                               
      DIMENSION       A(*)       ,C(*)       ,BP(*)      ,BH(*)      ,  
     1                CBP(*)                                            
      COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        ,  
     1                NM         ,NCMPLX     ,IK                        
      EXTERNAL        PSGF       ,PPSPF      ,PPSGF                     
      SCNV = SQRT(CNV)                                                  
      IZ = N                                                            
      IZM = IZ-1                                                        
      IZM2 = IZ-2                                                       
      IF (BP(N)-BP(1)) 101,142,103                                      
  101 DO 102 J=1,N                                                      
         NT = N-J                                                       
         BH(J) = BP(NT+1)                                               
  102 CONTINUE                                                          
      GO TO 105                                                         
  103 DO 104 J=1,N                                                      
         BH(J) = BP(J)                                                  
  104 CONTINUE                                                          
  105 NCMPLX = 0                                                        
      MODIZ = MOD(IZ,2)                                                 
      IS = 1                                                            
      IF (MODIZ) 106,107,106                                            
  106 IF (A(1)) 110,142,107                                             
  107 XL = BH(1)                                                        
      DB = BH(3)-BH(1)                                                  
  108 XL = XL-DB                                                        
      IF (PSGF(XL,IZ,C,A,BH)) 108,108,109                               
  109 SGN = -1.                                                         
      CBP(1) = CMPLX(CBSRH(XL,BH(1),IZ,C,A,BH,PSGF,SGN),0.)              
      IS = 2                                                            
  110 IF = IZ-1                                                         
      IF (MODIZ) 111,112,111                                            
  111 IF (A(1)) 112,142,115                                             
  112 XR = BH(IZ)                                                       
      DB = BH(IZ)-BH(IZ-2)                                              
  113 XR = XR+DB                                                        
      IF (PSGF(XR,IZ,C,A,BH)) 113,114,114                               
  114 SGN = 1.                                                          
      CBP(IZ) = CMPLX(CBSRH(BH(IZ),XR,IZ,C,A,BH,PSGF,SGN),0.)            
      IF = IZ-2                                                         
  115 DO 136 IG=IS,IF,2                                                 
         XL = BH(IG)                                                    
         XR = BH(IG+1)                                                  
         SGN = -1.                                                      
         XM = CBSRH(XL,XR,IZ,C,A,BH,PPSPF,SGN)                           
         PSG = PSGF(XM,IZ,C,A,BH)                                       
         IF (ABS(PSG)-EPS) 118,118,116                                  
  116    IF (PSG*PPSGF(XM,IZ,C,A,BH)) 117,118,119                       
C                                                                       
C     CASE OF A REAL ZERO                                               
C                                                                       
  117    SGN = 1.                                                       
         CBP(IG) = CMPLX(CBSRH(BH(IG),XM,IZ,C,A,BH,PSGF,SGN),0.)         
         SGN = -1.                                                      
         CBP(IG+1) = CMPLX(CBSRH(XM,BH(IG+1),IZ,C,A,BH,PSGF,SGN),0.)     
         GO TO 136                                                      
C                                                                       
C     CASE OF A MULTIPLE ZERO                                           
C                                                                       
  118    CBP(IG) = CMPLX(XM,0.)                                         
         CBP(IG+1) = CMPLX(XM,0.)                                       
         GO TO 136                                                      
C                                                                       
C     CASE OF A COMPLEX ZERO                                            
C                                                                       
  119    IT = 0                                                         
         ICV = 0                                                        
         CX = CMPLX(XM,0.)                                              
  120    FSG = (1.,0.)                                                  
         HSG = (1.,0.)                                                  
         FP = (0.,0.)                                                   
         FPP = (0.,0.)                                                  
         DO 121 J=1,IZ                                                  
            DD = 1./(CX-BH(J))                                          
            FSG = FSG*A(J)*DD                                           
            HSG = HSG*C(J)*DD                                           
            FP = FP+DD                                                  
            FPP = FPP-DD*DD                                             
  121    CONTINUE                                                       
         IF (MODIZ) 123,122,123                                         
  122    F = (1.,0.)-FSG-HSG                                            
         GO TO 124                                                      
  123    F = (1.,0.)+FSG+HSG                                            
  124    I3 = 0                                                         
         IF (CABS(FP)) 126,126,125                                      
  125    I3 = 1                                                         
         R3 = -F/FP                                                     
  126    I2 = 0                                                         
         IF (CABS(FPP)) 132,132,127                                     
  127    I2 = 1                                                         
         CDIS = CSQRT(FP**2-2.*F*FPP)                                   
         R1 = CDIS-FP                                                   
         R2 = -FP-CDIS                                                  
         IF (CABS(R1)-CABS(R2)) 129,129,128                             
  128    R1 = R1/FPP                                                    
         GO TO 130                                                      
  129    R1 = R2/FPP                                                    
  130    R2 = 2.*F/FPP/R1                                               
         IF (CABS(R2) .LT. CABS(R1)) R1 = R2                            
         IF (I3) 133,133,131                                            
  131    IF (CABS(R3) .LT. CABS(R1)) R1 = R3                            
         GO TO 133                                                      
  132    R1 = R3                                                        
  133    CX = CX+R1                                                     
         IT = IT+1                                                      
         IF (IT .GT. 50) GO TO 142                                      
         IF (CABS(R1) .GT. SCNV) GO TO 120                              
         IF (ICV) 134,134,135                                           
  134    ICV = 1                                                        
         GO TO 120                                                      
  135    CBP(IG) = CX                                                   
         CBP(IG+1) = CONJG(CX)                                          
  136 CONTINUE                                                          
      IF (CABS(CBP(N))-CABS(CBP(1))) 137,142,139                        
  137 NHALF = N/2                                                       
      DO 138 J=1,NHALF                                                  
         NT = N-J                                                       
         CX = CBP(J)                                                    
         CBP(J) = CBP(NT+1)                                             
         CBP(NT+1) = CX                                                 
  138 CONTINUE                                                          
  139 NCMPLX = 1                                                        
      DO 140 J=2,IZ                                                     
         IF (AIMAG(CBP(J))) 143,140,143                                 
  140 CONTINUE                                                          
      NCMPLX = 0                                                        
      DO 141 J=2,IZ                                                     
         BP(J) = REAL(CBP(J))                                           
  141 CONTINUE                                                          
      GO TO 143                                                         
  142 IERROR = 4                                                        
  143 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE PROC (ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,W,U)   
C                                                                       
C PROC APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND      
C STORES THE RESULT IN Y                                                
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS       
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY      
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X              
C NA IS THE LENGTH OF THE ARRAY AA                                      
C X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y       
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX                
C M  IS THE ORDER OF THE MATRIX                                         
C D,W,U ARE WORKING ARRAYS                                              
C IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE                
C                                                                       
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,X(*)       ,  
     1                Y(*)       ,D(*)       ,W(*)       ,BD(*)      ,  
     2                BM1(*)     ,BM2(*)     ,AA(*)      ,U(*)          
      COMPLEX         X          ,Y          ,A          ,B          ,  
     1                C          ,D          ,W          ,U          ,  
     2                DEN                                               
      DO 101 J=1,M                                                      
         W(J) = X(J)                                                    
         Y(J) = W(J)                                                    
  101 CONTINUE                                                          
      MM = M-1                                                          
      ID = ND                                                           
      IBR = 0                                                           
      M1 = NM1                                                          
      M2 = NM2                                                          
      IA = NA                                                           
  102 IF (IA) 105,105,103                                               
  103 RT = AA(IA)                                                       
      IF (ND .EQ. 0) RT = -RT                                           
      IA = IA-1                                                         
C                                                                       
C SCALAR MULTIPLICATION                                                 
C                                                                       
      DO 104 J=1,M                                                      
         Y(J) = RT*W(J)                                                 
  104 CONTINUE                                                          
  105 IF (ID) 125,125,106                                               
  106 RT = BD(ID)                                                       
      ID = ID-1                                                         
      IF (ID .EQ. 0) IBR = 1                                            
C                                                                       
C BEGIN SOLUTION TO SYSTEM                                              
C                                                                       
      D(M) = A(M)/(B(M)-RT)                                             
      W(M) = Y(M)/(B(M)-RT)                                             
      DO 107 J=2,MM                                                     
         K = M-J                                                        
         DEN = B(K+1)-RT-C(K+1)*D(K+2)                                  
         D(K+1) = A(K+1)/DEN                                            
         W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN                            
  107 CONTINUE                                                          
      DEN = B(1)-RT-C(1)*D(2)                                           
      W(1) = (1.,0.)                                                    
      IF (CABS(DEN)) 108,109,108                                        
  108 W(1) = (Y(1)-C(1)*W(2))/DEN                                       
  109 DO 110 J=2,M                                                      
         W(J) = W(J)-D(J)*W(J-1)                                        
  110 CONTINUE                                                          
      IF (NA) 113,113,102                                               
  111 DO 112 J=1,M                                                      
         Y(J) = W(J)                                                    
  112 CONTINUE                                                          
      IBR = 1                                                           
      GO TO 102                                                         
  113 IF (M1) 114,114,115                                               
  114 IF (M2) 111,111,120                                               
  115 IF (M2) 117,117,116                                               
  116 IF (ABS(BM1(M1))-ABS(BM2(M2))) 120,120,117                        
  117 IF (IBR) 118,118,119                                              
  118 IF (ABS(BM1(M1)-BD(ID))-ABS(BM1(M1)-RT)) 111,119,119              
  119 RT = RT-BM1(M1)                                                   
      M1 = M1-1                                                         
      GO TO 123                                                         
  120 IF (IBR) 121,121,122                                              
  121 IF (ABS(BM2(M2)-BD(ID))-ABS(BM2(M2)-RT)) 111,122,122              
  122 RT = RT-BM2(M2)                                                   
      M2 = M2-1                                                         
  123 DO 124 J=1,M                                                      
         Y(J) = Y(J)+RT*W(J)                                            
  124 CONTINUE                                                          
      GO TO 102                                                         
  125 RETURN                                                            
      END                                                               
      SUBROUTINE PROCP (ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,U,W)  
C                                                                       
C PROCP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND     
C STORES THE RESULT IN Y        PERIODIC BOUNDARY CONDITIONS            
C                                                                       
C BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS       
C ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY      
C AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X              
C NA IS THE LENGTH OF THE ARRAY AA                                      
C X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y       
C A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX                
C M  IS THE ORDER OF THE MATRIX                                         
C D,U,W ARE WORKING ARRAYS                                              
C IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE                
C                                                                       
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,X(*)       ,  
     1                Y(*)       ,D(*)       ,U(*)       ,BD(*)      ,  
     2                BM1(*)     ,BM2(*)     ,AA(*)      ,W(*)          
      COMPLEX         X          ,Y          ,A          ,B          ,  
     1                C          ,D          ,U          ,W          ,  
     2                DEN        ,YM         ,V          ,BH         ,AM
      DO 101 J=1,M                                                      
         Y(J) = X(J)                                                    
         W(J) = Y(J)                                                    
  101 CONTINUE                                                          
      MM = M-1                                                          
      MM2 = M-2                                                         
      ID = ND                                                           
      IBR = 0                                                           
      M1 = NM1                                                          
      M2 = NM2                                                          
      IA = NA                                                           
  102 IF (IA) 105,105,103                                               
  103 RT = AA(IA)                                                       
      IF (ND .EQ. 0) RT = -RT                                           
      IA = IA-1                                                         
      DO 104 J=1,M                                                      
         Y(J) = RT*W(J)                                                 
  104 CONTINUE                                                          
  105 IF (ID) 128,128,106                                               
  106 RT = BD(ID)                                                       
      ID = ID-1                                                         
      IF (ID .EQ. 0) IBR = 1                                            
C                                                                       
C BEGIN SOLUTION TO SYSTEM                                              
C                                                                       
      BH = B(M)-RT                                                      
      YM = Y(M)                                                         
      DEN = B(1)-RT                                                     
      D(1) = C(1)/DEN                                                   
      U(1) = A(1)/DEN                                                   
      W(1) = Y(1)/DEN                                                   
      V = C(M)                                                          
      IF (MM2-2) 109,107,107                                            
  107 DO 108 J=2,MM2                                                    
         DEN = B(J)-RT-A(J)*D(J-1)                                      
         D(J) = C(J)/DEN                                                
         U(J) = -A(J)*U(J-1)/DEN                                        
         W(J) = (Y(J)-A(J)*W(J-1))/DEN                                  
         BH = BH-V*U(J-1)                                               
         YM = YM-V*W(J-1)                                               
         V = -V*D(J-1)                                                  
  108 CONTINUE                                                          
  109 DEN = B(M-1)-RT-A(M-1)*D(M-2)                                     
      D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN                               
      W(M-1) = (Y(M-1)-A(M-1)*W(M-2))/DEN                               
      AM = A(M)-V*D(M-2)                                                
      BH = BH-V*U(M-2)                                                  
      YM = YM-V*W(M-2)                                                  
      DEN = BH-AM*D(M-1)                                                
      IF (CABS(DEN)) 110,111,110                                        
  110 W(M) = (YM-AM*W(M-1))/DEN                                         
      GO TO 112                                                         
  111 W(M) = (1.,0.)                                                    
  112 W(M-1) = W(M-1)-D(M-1)*W(M)                                       
      DO 113 J=2,MM                                                     
         K = M-J                                                        
         W(K) = W(K)-D(K)*W(K+1)-U(K)*W(M)                              
  113 CONTINUE                                                          
      IF (NA) 116,116,102                                               
  114 DO 115 J=1,M                                                      
         Y(J) = W(J)                                                    
  115 CONTINUE                                                          
      IBR = 1                                                           
      GO TO 102                                                         
  116 IF (M1) 117,117,118                                               
  117 IF (M2) 114,114,123                                               
  118 IF (M2) 120,120,119                                               
  119 IF (ABS(BM1(M1))-ABS(BM2(M2))) 123,123,120                        
  120 IF (IBR) 121,121,122                                              
  121 IF (ABS(BM1(M1)-BD(ID))-ABS(BM1(M1)-RT)) 114,122,122              
  122 RT = RT-BM1(M1)                                                   
      M1 = M1-1                                                         
      GO TO 126                                                         
  123 IF (IBR) 124,124,125                                              
  124 IF (ABS(BM2(M2)-BD(ID))-ABS(BM2(M2)-RT)) 114,125,125              
  125 RT = RT-BM2(M2)                                                   
      M2 = M2-1                                                         
  126 DO 127 J=1,M                                                      
         Y(J) = Y(J)+RT*W(J)                                            
  127 CONTINUE                                                          
      GO TO 102                                                         
  128 RETURN                                                            
      END                                                               
      SUBROUTINE CTEVLS (N,D,E2,IERR)                                    
C                                                                       
      INTEGER         I          ,J          ,L          ,M          ,  
     1                N          ,II         ,L1         ,MML        ,  
     2                IERR                                              
      REAL            D(N)       ,E2(N)                                 
      REAL            B          ,C          ,F          ,G          ,  
     1                H          ,P          ,R          ,S          ,  
     2                MACHEP                                            
C                                                                       
C     REAL SQRT,ABS,SIGN                                                
C                                                                       
      COMMON /CCBLK/  NPP        ,K          ,MACHEP     ,CNV        ,  
     1                NM         ,NCMPLX     ,IK                        
C                                                                       
C     THIS SUBROUTINE IS A MODIFICATION OF THE EISPACK SUBROUTINE TQLRAT
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.                
C                                                                       
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC              
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.                     
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          
C                                                                       
C        E2 CONTAINS THE                SUBDIAGONAL ELEMENTS OF THE     
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY. 
C                                                                       
C      ON OUTPUT-                                                       
C                                                                       
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND          
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE            
C          THE SMALLEST EIGENVALUES,                                    
C                                                                       
C        E2 HAS BEEN DESTROYED,                                         
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C                     DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C                                                                       
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                                                                       
C                **********                                             
C                                                                       
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 115                                           
C                                                                       
      DO 101 I=2,N                                                      
         E2(I-1) = E2(I)*E2(I)                                          
  101 CONTINUE                                                          
C                                                                       
      F = 0.0                                                           
      B = 0.0                                                           
      E2(N) = 0.0                                                       
C                                                                       
      DO 112 L=1,N                                                      
         J = 0                                                          
         H = MACHEP*(ABS(D(L))+SQRT(E2(L)))                             
         IF (B .GT. H) GO TO 102                                        
         B = H                                                          
         C = B*B                                                        
C                                                                       
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ********** 
C                                                                       
  102    DO 103 M=L,N                                                   
            IF (E2(M) .LE. C) GO TO 104                                 
C                                                                       
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT              
C                THROUGH THE BOTTOM OF THE LOOP **********              
C                                                                       
  103    CONTINUE                                                       
C                                                                       
  104    IF (M .EQ. L) GO TO 108                                        
  105    IF (J .EQ. 30) GO TO 114                                       
         J = J+1                                                        
C                                                                       
C     ********** FORM SHIFT **********                                  
C                                                                       
         L1 = L+1                                                       
         S = SQRT(E2(L))                                                
         G = D(L)                                                       
         P = (D(L1)-G)/(2.0*S)                                          
         R = SQRT(P*P+1.0)                                              
         D(L) = S/(P+SIGN(R,P))                                         
         H = G-D(L)                                                     
C                                                                       
         DO 106 I=L1,N                                                  
            D(I) = D(I)-H                                               
  106    CONTINUE                                                       
C                                                                       
         F = F+H                                                        
C                                                                       
C     ********** RATIONAL QL TRANSFORMATION **********                  
C                                                                       
         G = D(M)                                                       
         IF (G .EQ. 0.0) G = B                                          
         H = G                                                          
         S = 0.0                                                        
         MML = M-L                                                      
C                                                                       
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********             
C                                                                       
         DO 107 II=1,MML                                                
            I = M-II                                                    
            P = G*H                                                     
            R = P+E2(I)                                                 
            E2(I+1) = S*R                                               
            S = E2(I)/R                                                 
            D(I+1) = H+S*(H+D(I))                                       
            G = D(I)-E2(I)/G                                            
            IF (G .EQ. 0.0) G = B                                       
            H = G*P/R                                                   
  107    CONTINUE                                                       
C                                                                       
         E2(L) = S*G                                                    
         D(L) = H                                                       
C                                                                       
C     ********** GUARD AGAINST UNDERFLOWED H **********                 
C                                                                       
         IF (H .EQ. 0.0) GO TO 108                                      
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 108                        
         E2(L) = H*E2(L)                                                
         IF (E2(L) .NE. 0.0) GO TO 105                                  
  108    P = D(L)+F                                                     
C                                                                       
C     ********** ORDER EIGENVALUES **********                           
C                                                                       
         IF (L .EQ. 1) GO TO 110                                        
C                                                                       
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********               
C                                                                       
         DO 109 II=2,L                                                  
            I = L+2-II                                                  
            IF (P .GE. D(I-1)) GO TO 111                                
            D(I) = D(I-1)                                               
  109    CONTINUE                                                       
C                                                                       
  110    I = 1                                                          
  111    D(I) = P                                                       
  112 CONTINUE                                                          
C                                                                       
      IF (ABS(D(N)) .GE. ABS(D(1))) GO TO 115                           
      NHALF = N/2                                                       
      DO 113 I=1,NHALF                                                  
         NTOP = N-I                                                     
         DHOLD = D(I)                                                   
         D(I) = D(NTOP+1)                                               
         D(NTOP+1) = DHOLD                                              
  113 CONTINUE                                                          
      GO TO 115                                                         
C                                                                       
C     ********** SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS **********              
C                                                                       
  114 IERR = L                                                          
  115 RETURN                                                            
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
