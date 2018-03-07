C
C     file sepaux.f
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
C PACKAGE SEPAUX         CONTAINS NO USER ENTRY POINTS.                 
C                                                                       
C LATEST REVISION        NOVEMBER 1988
C                                                                       
C PURPOSE                THIS PACKAGE CONTAINS AUXILIARY ROUTINES FOR   
C                        NCAR PUBLIC SOFTWARE PACKAGES SUCH AS SEPELI   
C                        AND SEPX4.                                     
C                                                                       
C USAGE                  SINCE THIS PACKAGE CONTAINS NO USER ENTRIES,   
C                        NO USAGE INSTRUCTIONS OR ARGUMENT DESCRIPTIONS 
C                        ARE GIVEN HERE.                                
C                                                                       
C SPECIAL CONDITIONS     NONE                                           
C                                                                       
C I/O                    NONE                                           
C                                                                       
C PRECISION              SINGLE                                         
C                                                                       
C REQUIRED LIBRARY       NONE                                           
C FILES                                                                 
C                                                                       
C LANGUAGE               FORTRAN                                        
C                                                                       
C HISTORY                DEVELOPED IN THE LATE 1970'S BY JOHN C. ADAMS  
C                        OF NCAR'S SCIENTTIFIC COMPUTING DIVISION.      
C                                                                       
C PORTABILITY            FORTRAN 77                                     
C **********************************************************************
      SUBROUTINE SEPORT (USOL,IDMN,ZN,ZM,PERTRB)                        
C                                                                       
C     THIS SUBROUTINE ORTHOGONALIZES THE ARRAY USOL WITH RESPECT TO     
C     THE CONSTANT ARRAY IN A WEIGHTED LEAST SQUARES NORM               
C                                                                       
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,  
     1                AIT        ,BIT        ,CIT        ,DIT        ,  
     2                MIT        ,NIT        ,IS         ,MS         ,  
     3                JS         ,NS         ,DLX        ,DLY        ,  
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4          
      DIMENSION       USOL(IDMN,1)           ,ZN(*)      ,ZM(*)         
      ISTR = IS                                                         
      IFNL = MS                                                         
      JSTR = JS                                                         
      JFNL = NS                                                         
C                                                                       
C     COMPUTE WEIGHTED INNER PRODUCTS                                   
C                                                                       
      UTE = 0.0                                                         
      ETE = 0.0                                                         
      DO  20 I=IS,MS                                                    
         II = I-IS+1                                                    
         DO  10 J=JS,NS                                                 
            JJ = J-JS+1                                                 
            ETE = ETE+ZM(II)*ZN(JJ)                                     
            UTE = UTE+USOL(I,J)*ZM(II)*ZN(JJ)                           
   10    CONTINUE                                                       
   20 CONTINUE                                                          
C                                                                       
C     SET PERTURBATION PARAMETER                                        
C                                                                       
      PERTRB = UTE/ETE                                                  
C                                                                       
C     SUBTRACT OFF CONSTANT PERTRB                                      
C                                                                       
      DO  40 I=ISTR,IFNL                                                
         DO  30 J=JSTR,JFNL                                             
            USOL(I,J) = USOL(I,J)-PERTRB                                
   30    CONTINUE                                                       
   40 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE SEPMIN (USOL,IDMN,ZN,ZM,PERTB)                         
C                                                                       
C     THIS SUBROUTINE ORHTOGONALIZES THE ARRAY USOL WITH RESPECT TO     
C     THE CONSTANT ARRAY IN A WEIGHTED LEAST SQUARES NORM               
C                                                                       
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,  
     1                AIT        ,BIT        ,CIT        ,DIT        ,  
     2                MIT        ,NIT        ,IS         ,MS         ,  
     3                JS         ,NS         ,DLX        ,DLY        ,  
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4          
      DIMENSION       USOL(IDMN,1)           ,ZN(*)      ,ZM(*)         
C                                                                       
C     ENTRY AT SEPMIN OCCURRS WHEN THE FINAL SOLUTION IS                
C     TO BE MINIMIZED WITH RESPECT TO THE WEIGHTED                      
C     LEAST SQUARES NORM                                                
C                                                                       
      ISTR = 1                                                          
      IFNL = K                                                          
      JSTR = 1                                                          
      JFNL = L                                                          
C                                                                       
C     COMPUTE WEIGHTED INNER PRODUCTS                                   
C                                                                       
      UTE = 0.0                                                         
      ETE = 0.0                                                         
      DO  20 I=IS,MS                                                    
         II = I-IS+1                                                    
         DO  10 J=JS,NS                                                 
            JJ = J-JS+1                                                 
            ETE = ETE+ZM(II)*ZN(JJ)                                     
            UTE = UTE+USOL(I,J)*ZM(II)*ZN(JJ)                           
   10    CONTINUE                                                       
   20 CONTINUE                                                          
C                                                                       
C     SET PERTURBATION PARAMETER                                        
C                                                                       
      PERTRB = UTE/ETE                                                  
C                                                                       
C     SUBTRACT OFF CONSTANT PERTRB                                      
C                                                                       
      DO  40 I=ISTR,IFNL                                                
         DO  30 J=JSTR,JFNL                                             
            USOL(I,J) = USOL(I,J)-PERTRB                                
   30    CONTINUE                                                       
   40 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE SEPTRI (N,A,B,C,D,U,Z)                                 
C                                                                       
C     THIS SUBROUTINE SOLVES FOR A NON-ZERO EIGENVECTOR CORRESPONDING   
C     TO THE ZERO EIGENVALUE OF THE TRANSPOSE OF THE RANK               
C     DEFICIENT ONE MATRIX WITH SUBDIAGONAL A, DIAGONAL B, AND          
C     SUPERDIAGONAL C , WITH A(1) IN THE (1,N) POSITION, WITH           
C     C(N) IN THE (N,1) POSITION, AND ALL OTHER ELEMENTS ZERO.          
C                                                                       
      DIMENSION       A(N)       ,B(N)       ,C(N)       ,D(N)       ,  
     1                U(N)       ,Z(N)                                  
      BN = B(N)                                                         
      D(1) = A(2)/B(1)                                                  
      V = A(1)                                                          
      U(1) = C(N)/B(1)                                                  
      NM2 = N-2                                                         
      DO  10 J=2,NM2                                                    
         DEN = B(J)-C(J-1)*D(J-1)                                       
         D(J) = A(J+1)/DEN                                              
         U(J) = -C(J-1)*U(J-1)/DEN                                      
         BN = BN-V*U(J-1)                                               
         V = -V*D(J-1)                                                  
   10 CONTINUE                                                          
      DEN = B(N-1)-C(N-2)*D(N-2)                                        
      D(N-1) = (A(N)-C(N-2)*U(N-2))/DEN                                 
      AN = C(N-1)-V*D(N-2)                                              
      BN = BN-V*U(N-2)                                                  
      DEN = BN-AN*D(N-1)                                                
C                                                                       
C     SET LAST COMPONENT EQUAL TO ONE                                   
C                                                                       
      Z(N) = 1.0                                                        
      Z(N-1) = -D(N-1)                                                  
      NM1 = N-1                                                         
      DO  20 J=2,NM1                                                    
         K = N-J                                                        
         Z(K) = -D(K)*Z(K+1)-U(K)*Z(N)                                  
   20 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE SEPDX (U,IDMN,I,J,UXXX,UXXXX)                          
C                                                                       
C     THIS PROGRAM COMPUTES SECOND ORDER FINITE DIFFERENCE              
C     APPROXIMATIONS TO THE THIRD AND FOURTH X                          
C     PARTIAL DERIVATIVES OF U AT THE (I,J) MESH POINT                  
C                                                                       
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,  
     1                AIT        ,BIT        ,CIT        ,DIT        ,  
     2                MIT        ,NIT        ,IS         ,MS         ,  
     3                JS         ,NS         ,DLX        ,DLY        ,  
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4          
      DIMENSION       U(IDMN,1)                                         
      IF (I.GT.2 .AND. I.LT.(K-1)) GO TO  50                            
      IF (I .EQ. 1) GO TO  10                                           
      IF (I .EQ. 2) GO TO  30                                           
      IF (I .EQ. K-1) GO TO  60                                         
      IF (I .EQ. K) GO TO  80                                           
C                                                                       
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A                  
C                                                                       
   10 IF (KSWX .EQ. 1) GO TO  20                                        
      UXXX = (-5.0*U(1,J)+18.0*U(2,J)-24.0*U(3,J)+14.0*U(4,J)-          
     1                                               3.0*U(5,J))/(TDLX3)
      UXXXX = (3.0*U(1,J)-14.0*U(2,J)+26.0*U(3,J)-24.0*U(4,J)+          
     1                                      11.0*U(5,J)-2.0*U(6,J))/DLX4
      RETURN                                                            
C                                                                       
C     PERIODIC AT X=A                                                   
C                                                                       
   20 UXXX = (-U(K-2,J)+2.0*U(K-1,J)-2.0*U(2,J)+U(3,J))/(TDLX3)         
      UXXXX = (U(K-2,J)-4.0*U(K-1,J)+6.0*U(1,J)-4.0*U(2,J)+U(3,J))/DLX4 
      RETURN                                                            
C                                                                       
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A+DLX              
C                                                                       
   30 IF (KSWX .EQ. 1) GO TO  40                                        
      UXXX = (-3.0*U(1,J)+10.0*U(2,J)-12.0*U(3,J)+6.0*U(4,J)-U(5,J))/   
     1       TDLX3                                                      
      UXXXX = (2.0*U(1,J)-9.0*U(2,J)+16.0*U(3,J)-14.0*U(4,J)+6.0*U(5,J)-
     1                                                      U(6,J))/DLX4
      RETURN                                                            
C                                                                       
C     PERIODIC AT X=A+DLX                                               
C                                                                       
   40 UXXX = (-U(K-1,J)+2.0*U(1,J)-2.0*U(3,J)+U(4,J))/(TDLX3)           
      UXXXX = (U(K-1,J)-4.0*U(1,J)+6.0*U(2,J)-4.0*U(3,J)+U(4,J))/DLX4   
      RETURN                                                            
C                                                                       
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR         
C                                                                       
   50 CONTINUE                                                          
      UXXX = (-U(I-2,J)+2.0*U(I-1,J)-2.0*U(I+1,J)+U(I+2,J))/TDLX3       
      UXXXX = (U(I-2,J)-4.0*U(I-1,J)+6.0*U(I,J)-4.0*U(I+1,J)+U(I+2,J))/ 
     1        DLX4                                                      
      RETURN                                                            
C                                                                       
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B-DLX              
C                                                                       
   60 IF (KSWX .EQ. 1) GO TO  70                                        
      UXXX = (U(K-4,J)-6.0*U(K-3,J)+12.0*U(K-2,J)-10.0*U(K-1,J)+        
     1                                                 3.0*U(K,J))/TDLX3
      UXXXX = (-U(K-5,J)+6.0*U(K-4,J)-14.0*U(K-3,J)+16.0*U(K-2,J)-      
     1                                     9.0*U(K-1,J)+2.0*U(K,J))/DLX4
      RETURN                                                            
C                                                                       
C     PERIODIC AT X=B-DLX                                               
C                                                                       
   70 UXXX = (-U(K-3,J)+2.0*U(K-2,J)-2.0*U(1,J)+U(2,J))/TDLX3           
      UXXXX = (U(K-3,J)-4.0*U(K-2,J)+6.0*U(K-1,J)-4.0*U(1,J)+U(2,J))/   
     1        DLX4                                                      
      RETURN                                                            
C                                                                       
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B                  
C                                                                       
   80 UXXX = -(3.0*U(K-4,J)-14.0*U(K-3,J)+24.0*U(K-2,J)-18.0*U(K-1,J)+  
     1                                                 5.0*U(K,J))/TDLX3
      UXXXX = (-2.0*U(K-5,J)+11.0*U(K-4,J)-24.0*U(K-3,J)+26.0*U(K-2,J)- 
     1                                    14.0*U(K-1,J)+3.0*U(K,J))/DLX4
      RETURN                                                            
      END                                                               
      SUBROUTINE SEPDY (U,IDMN,I,J,UYYY,UYYYY)                          
C                                                                       
C     THIS PROGRAM COMPUTES SECOND ORDER FINITE DIFFERENCE              
C     APPROXIMATIONS TO THE THIRD AND FOURTH Y                          
C     PARTIAL DERIVATIVES OF U AT THE (I,J) MESH POINT                  
C                                                                       
      COMMON /SPLP/   KSWX       ,KSWY       ,K          ,L          ,  
     1                AIT        ,BIT        ,CIT        ,DIT        ,  
     2                MIT        ,NIT        ,IS         ,MS         ,  
     3                JS         ,NS         ,DLX        ,DLY        ,  
     4                TDLX3      ,TDLY3      ,DLX4       ,DLY4          
      DIMENSION       U(IDMN,*)                                         
      IF (J.GT.2 .AND. J.LT.(L-1)) GO TO  50                            
      IF (J .EQ. 1) GO TO  10                                           
      IF (J .EQ. 2) GO TO  30                                           
      IF (J .EQ. L-1) GO TO  60                                         
      IF (J .EQ. L) GO TO  80                                           
C                                                                       
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C                  
C                                                                       
   10 IF (KSWY .EQ. 1) GO TO  20                                        
      UYYY = (-5.0*U(I,1)+18.0*U(I,2)-24.0*U(I,3)+14.0*U(I,4)-          
     1                                                 3.0*U(I,5))/TDLY3
      UYYYY = (3.0*U(I,1)-14.0*U(I,2)+26.0*U(I,3)-24.0*U(I,4)+          
     1                                      11.0*U(I,5)-2.0*U(I,6))/DLY4
      RETURN                                                            
C                                                                       
C     PERIODIC AT X=A                                                   
C                                                                       
   20 UYYY = (-U(I,L-2)+2.0*U(I,L-1)-2.0*U(I,2)+U(I,3))/TDLY3           
      UYYYY = (U(I,L-2)-4.0*U(I,L-1)+6.0*U(I,1)-4.0*U(I,2)+U(I,3))/DLY4 
      RETURN                                                            
C                                                                       
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C+DLY              
C                                                                       
   30 IF (KSWY .EQ. 1) GO TO  40                                        
      UYYY = (-3.0*U(I,1)+10.0*U(I,2)-12.0*U(I,3)+6.0*U(I,4)-U(I,5))/   
     1       TDLY3                                                      
      UYYYY = (2.0*U(I,1)-9.0*U(I,2)+16.0*U(I,3)-14.0*U(I,4)+6.0*U(I,5)-
     1                                                      U(I,6))/DLY4
      RETURN                                                            
C                                                                       
C     PERIODIC AT Y=C+DLY                                               
C                                                                       
   40 UYYY = (-U(I,L-1)+2.0*U(I,1)-2.0*U(I,3)+U(I,4))/TDLY3             
      UYYYY = (U(I,L-1)-4.0*U(I,1)+6.0*U(I,2)-4.0*U(I,3)+U(I,4))/DLY4   
      RETURN                                                            
C                                                                       
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR         
C                                                                       
   50 CONTINUE                                                          
      UYYY = (-U(I,J-2)+2.0*U(I,J-1)-2.0*U(I,J+1)+U(I,J+2))/TDLY3       
      UYYYY = (U(I,J-2)-4.0*U(I,J-1)+6.0*U(I,J)-4.0*U(I,J+1)+U(I,J+2))/ 
     1        DLY4                                                      
      RETURN                                                            
C                                                                       
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D-DLY              
C                                                                       
   60 IF (KSWY .EQ. 1) GO TO  70                                        
      UYYY = (U(I,L-4)-6.0*U(I,L-3)+12.0*U(I,L-2)-10.0*U(I,L-1)+        
     1                                                 3.0*U(I,L))/TDLY3
      UYYYY = (-U(I,L-5)+6.0*U(I,L-4)-14.0*U(I,L-3)+16.0*U(I,L-2)-      
     1                                     9.0*U(I,L-1)+2.0*U(I,L))/DLY4
      RETURN                                                            
C                                                                       
C     PERIODIC AT Y=D-DLY                                               
C                                                                       
   70 CONTINUE                                                          
      UYYY = (-U(I,L-3)+2.0*U(I,L-2)-2.0*U(I,1)+U(I,2))/TDLY3           
      UYYYY = (U(I,L-3)-4.0*U(I,L-2)+6.0*U(I,L-1)-4.0*U(I,1)+U(I,2))/   
     1        DLY4                                                      
      RETURN                                                            
C                                                                       
C     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D                  
C                                                                       
   80 UYYY = -(3.0*U(I,L-4)-14.0*U(I,L-3)+24.0*U(I,L-2)-18.0*U(I,L-1)+  
     1                                                 5.0*U(I,L))/TDLY3
      UYYYY = (-2.0*U(I,L-5)+11.0*U(I,L-4)-24.0*U(I,L-3)+26.0*U(I,L-2)- 
     1                                    14.0*U(I,L-1)+3.0*U(I,L))/DLY4
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
