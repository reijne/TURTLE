      SUBROUTINE MKNIT(N,NIT,IJ,ISYM,JSYM)                                      
      Parameter(NDI12=530)
      DIMENSION NIT(*)                                                          
      INTEGER*2 IJ(8),ISYM(8),JSYM(36)                                          
      integer*4 ipq
      Dimension IPQ(NDI12+1)                                               
cdebug      write(6,*)"isym"
cdebug      write(6,*) isym
cdebug      write(6,*)"jsym"
cdebug      write(6,*) jsym
      ipq(1)=0 
      Do i=1,ndi12
       ipq(i+1)=ipq(i)+i
      End Do
      DO 30 I=1,667                                                             
        NIT(I)=0
30    continue                                                                  
      IW=0                                                                      
      LG=0                                                                      
      DO 102 I=1,N
       NAI=IJ(I) 
       NIP=ISYM(I)
       DO 102 J=1,I
        NAJ=IJ(J) 
        NJP=ISYM(J)
         IF(NJP.GT.NIP) GO TO 103          
           NIJ=IPQ(NIP)+NJP                                        
           GO TO 104    
103      continue
         NIJ=IPQ(NJP)+NIP       
104      continue
         NIJ=JSYM(NIJ)     
         DO 102 K=1,I   
           NAK=IJ(K)        
           NKP=ISYM(K)                       
           LIM=K    
           IF(I.EQ.K) LIM=J      
      DO 102 L=1,LIM                                                            
      IW=IW+1                                                                   
      NIT(IW)=LG                                                                
      NLP=ISYM(L)                                                               
      IF(NLP.GT.NKP) GO TO 105                                                  
      NJJ=IPQ(NKP)+NLP                                                          
      GO TO 412                                                                 
105   NJJ=IPQ(NLP)+NKP                                                          
412   NJJ=JSYM(NJJ)                                                            
      IF(NIJ.NE.NJJ) GO TO 102                                                  
      IF(I.EQ.J) GO TO 106                                                      
      ISU=NAI*NAJ                                                               
      IF(I.EQ.K) GO TO 107                                                      
      JSU=NAK*IJ(L)                                                             
108   LG=LG+ISU*JSU                                                             
      GO TO 102                                                                 
106   ISU=IPQ(NAI+1)                                                            
      IF(I.EQ.K) GO TO 107                                                      
      JSU=IPQ(NAK+1)                                                            
      GO TO 108                                                                 
107   LG=LG+ISU*(ISU+1)/2                                                       
102   CONTINUE                                                                  
      IW=IW+1                                                                   
      NIT(IW)=LG                                                                
      NIT(667)=LG                                                               
      RETURN                                                                    
      END                                                                       
