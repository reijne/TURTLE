      SUBROUTINE WRTTXT(IPRT)                                                   
      COMMON /KEYWRD/ KEYWRD                                                    
      COMMON /TITLES/ KOMENT,TITLE                                              
      CHARACTER KEYWRD*241, KOMENT*81, TITLE*81                                 
      DO 10 I=81,2,-1                                                           
   10 IF(KEYWRD(I:I).NE.' ')GOTO 20                                             
   20 WRITE(IPRT,'(A)')KEYWRD(:I)                                               
      IF(INDEX(KEYWRD(1:81),' +')+INDEX(KEYWRD(1:81),'&')+                      
     1INDEX(KEYWRD,'SETUP').NE.0) THEN                                          
         DO 30 I=161,82,-1                                                      
   30    IF(KEYWRD(I:I).NE.' ')GOTO 40                                          
   40    IF(KEYWRD(81:81).NE.' ')THEN                                           
            WRITE(IPRT,'(1X,A)')KEYWRD(81:I)                                    
         ELSE                                                                   
            WRITE(IPRT,'(A)')KEYWRD(81:I)                                       
         ENDIF                                                                  
      ENDIF                                                                     
      IF(INDEX(KEYWRD(81:241),' +')+INDEX(KEYWRD(81:241),'&') +                 
     1INDEX(KEYWRD,'SETUP').NE.0) THEN                                          
         DO 50 I=241,161,-1                                                     
   50    IF(KEYWRD(I:I).NE.' ')GOTO 60                                          
   60    IF(KEYWRD(161:161).NE.' ')THEN                                         
            WRITE(IPRT,'(1X,A)')KEYWRD(161:I)                                   
         ELSE                                                                   
            WRITE(IPRT,'(A)')KEYWRD(161:I)                                      
         ENDIF                                                                  
      ENDIF                                                                     
      DO 70 I=81,2,-1                                                           
   70 IF(KOMENT(I:I).NE.' ')GOTO 80                                             
   80 IF(INDEX(KOMENT,' NULL ').EQ.0)WRITE(IPRT,'(A)')KOMENT(:I)                
      DO 90 I=81,2,-1                                                           
   90 IF(TITLE(I:I).NE.' ')GOTO 100                                             
  100 IF(INDEX(TITLE,' NULL ').EQ.0)WRITE(IPRT,'(A)')TITLE(:I)                  
      END                                                                       
