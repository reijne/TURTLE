cvp                                                                             
c   aufbau des feldes igew (gewichte der vertizes im verzweigungsdiagramm)
c    unter verwendung der reversed lexical ordering nach
c     duch u. karwowski, comp. phys. rep. vol. 2 (1985) p.122
      subroutine weight(imult,nbeta,inopen)
c
      include"sga.inc"
c
c   hilfsfeld zur rekursiven berechnung der gewichte
      integer idumm(0:mxsdg+1)
c  
c initialisierung von igew
      do 10 j=1,mxsdg
        do 10 i=0,nopmax
          igew(i,j)=0
   10 continue  
      do 20 i=0,mxsdg+1
        idumm(i)=0
   20 continue  
      idumm(imult)=1
      igew(0,imult)=1
c  
c aufbau von igew
c  schleife ueber die einfach besetzten schalen 
c   ismin, ismax: minimaler und maximaler wert der multiplizitaet bei 
c                 der jeweiligen offenen schale
      ismin=imult
      ismax=imult
      do 1000 i=1,inopen
c    berechnung der oberen und unteren werte von 2s+1, die fuer i moeglich sind
        ismin=ismin-1
        if (ismin.lt.1) ismin=2
        if (i.le.nbeta) then
          ismax=ismax+1
         else
          ismax=ismax-1
        endif
c    berechnung der gewichte durch addition der werte aus i-1
        do 2010 j=ismin,ismax,2
          igew(i,j)=idumm(j-1)+idumm(j+1)
c         igew(i,j)=igew(i-1,j-1)+igew(i-1,j+1)
 2010   continue  
c    update von idumm
        do 2020 j=0,mxsdg+1
          idumm(j)=0
 2020   continue  
        do 2030 j=ismin,ismax,2
          idumm(j)=igew(i,j)
 2030   continue  
c  
 1000 continue  
c  
c  
      return                                                                    
      end                                                                       
