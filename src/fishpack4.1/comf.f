C
C     file comf.f
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
C PACKAGE COMF           THE ENTRIES IN THIS PACKAGE ARE LOWLEVEL       
C                        ENTRIES, SUPPORTING FISHPACK ENTRIES BLKTRI
C                        AND CBLKTRI. THAT IS, THESE ROUTINES ARE       
C                        NOT CALLED DIRECTLY BY USERS, BUT RATHER       
C                        BY ENTRIES WITHIN BLKTRI AND CBLKTRI.          
C                        DESCRIPTION OF ENTRIES EPMACH AND PIMACH       
C                        FOLLOW BELOW.                                  
C                                                                       
C LATEST REVISION        NOVEMBER 1988
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
C ********************************************************************  
C                                                                       
C FUNCTION EPMACH (DUM)                                                 
C                                                                       
C PURPOSE                TO COMPUTE AN APPROXIMATE MACHINE ACCURACY     
C                        EPSILON ACCORDING TO THE FOLLOWING DEFINITION: 
C                        EPSILON IS THE SMALLEST NUMBER SUCH THAT       
C                        (1.+EPSILON).GT.1.)                            
C                                                                       
C USAGE                  EPS = EPMACH (DUM)                             
C                                                                       
C ARGUMENTS                                                             
C ON INPUT               DUM                                            
C                          DUMMY VALUE                                  
C                                                                       
C ARGUMENTS                                                             
C ON OUTPUT              NONE                                           
C                                                                       
C HISTORY                THE ORIGINAL VERSION, WRITTEN WHEN THE         
C                        BLKTRI PACKAGE WAS CONVERTED FROM THE          
C                        CDC 7600 TO RUN ON THE CRAY-1, CALCULATED      
C                        MACHINE ACCURACY BY SUCCESSIVE DIVISIONS       
C                        BY 10.  USE OF THIS CONSTANT CAUSED BLKTRI     
C                        TO COMPUTE SOLUTIONS ON THE CRAY-1 WITH FOUR   
C                        FEWER PLACES OF ACCURACY THAN THE VERSION      
C                        ON THE 7600.  IT WAS FOUND THAT COMPUTING      
C                        MACHINE ACCURACY BY SUCCESSIVE DIVISIONS       
C                        OF 2 PRODUCED A MACHINE ACCURACY 29% LESS      
C                        THAN THE VALUE GENERATED BY SUCCESSIVE         
C                        DIVISIONS BY 10, AND THAT USE OF THIS          
C                        MACHINE CONSTANT IN THE BLKTRI PACKAGE         
C                        RECOVERED THE ACCURACY THAT APPEARED TO        
C                        BE LOST ON CONVERSION.                         
C                                                                       
C ALGORITHM              COMPUTES MACHINE ACCURACY BY SUCCESSIVE        
C                        DIVISIONS OF TWO.                              
C                                                                       
C PORTABILITY            THIS CODE WILL EXECUTE ON MACHINES OTHER       
C                        THAN THE CRAY1, BUT THE RETURNED VALUE MAY     
C                        BE UNSATISFACTORY.  SEE HISTORY ABOVE.         
C ********************************************************************  
C                                                                       
C 11/26/2011
C This routine commented out because it also appears in fftpack
C giving a compile-time conflict when making libfishpack.a.
C                                                                       
C FUNCTION PIMACH (DUM)                                                 
C                                                                       
C PURPOSE                TO SUPPLY THE VALUE OF THE CONSTANT PI         
C                        CORRECT TO MACHINE PRECISION WHERE             
C                        PI=3.141592653589793238462643383279502884197   
C                             1693993751058209749446                    
C                                                                       
C USAGE                  PI = PIMACH (DUM)                              
C                                                                       
C ARGUMENTS                                                             
C ON INPUT               DUM                                            
C                          DUMMY VALUE                                  
C                                                                       
C ARGUMENTS                                                             
C ON OUTPUT              NONE                                           
C                                                                       
C ALGORITHM              THE VALUE OF PI IS SET TO 4.*ATAN(1.0)
C                                                                       
C PORTABILITY            THIS ENTRY IS PORTABLE, BUT USERS SHOULD       
C                        CHECK TO SEE WHETHER GREATER ACCURACY IS       
C                        REQUIRED.                                      
C                                                                       
C***********************************************************************
      FUNCTION EPMACH (DUM)                                             
      COMMON /VALUE/  V                                                 
      EPS = 1.                                                          
  101 EPS = EPS/2.                                                      
      CALL STRWRD (EPS+1.)                                               
      IF (V-1.) 102,102,101                                             
  102 EPMACH = 100.*EPS                                                 
      RETURN                                                            
      END                                                               
      SUBROUTINE STRWRD (X)                                              
      COMMON /VALUE/  V                                                 
      V = X                                                             
      RETURN                                                            
      END                                                               
C     FUNCTION PIMACH (DUM)                                             
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C                                                                       
C     PIMACH = 4.*ATAN(1.0)
C     RETURN                                                            
C     END                                                               
      FUNCTION PPSGF (X,IZ,C,A,BH)                                      
      DIMENSION       A(*)       ,C(*)       ,BH(*)                     
      SUM = 0.                                                          
      DO 101 J=1,IZ                                                     
         SUM = SUM-1./(X-BH(J))**2                                      
  101 CONTINUE                                                          
      PPSGF = SUM                                                       
      RETURN                                                            
      END                                                               
      FUNCTION PPSPF (X,IZ,C,A,BH)                                      
      DIMENSION       A(*)       ,C(*)       ,BH(*)                     
      SUM = 0.                                                          
      DO 101 J=1,IZ                                                     
         SUM = SUM+1./(X-BH(J))                                         
  101 CONTINUE                                                          
      PPSPF = SUM                                                       
      RETURN                                                            
      END                                                               
      FUNCTION PSGF (X,IZ,C,A,BH)                                       
      DIMENSION       A(*)       ,C(*)       ,BH(*)                     
      FSG = 1.                                                          
      HSG = 1.                                                          
      DO 101 J=1,IZ                                                     
         DD = 1./(X-BH(J))                                              
         FSG = FSG*A(J)*DD                                              
         HSG = HSG*C(J)*DD                                              
  101 CONTINUE                                                          
      IF (MOD(IZ,2)) 103,102,103                                        
  102 PSGF = 1.-FSG-HSG                                                 
      RETURN                                                            
  103 PSGF = 1.+FSG+HSG                                                 
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
