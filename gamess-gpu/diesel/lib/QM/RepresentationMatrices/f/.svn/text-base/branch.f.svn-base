cvp                                                                             
c   aufbau des feldes idiagr (abbildung des verzweigungsdiagramms)  
c    hier verwendung der reversed lexical ordering nach
c     duch u. karwowski, comp. phys. rep. vol. 2 (1985) p.122
      subroutine branch(neig,nbeta,inopen)
c  
      include"sga.inc"
c  
c  vorbelegung von idiagr
      do 100 j=1,inopen
        do 100 i=1,neig
          idiagr(i,j)=1
  100 continue
      do 105 j=1,mxiswh+1
        igraph(j)=0
  105 continue
c
c  letzte spinfkt in reversed lexical ordering
      do 110 i=1,nbeta
        idiagr(neig,i)=-1
        igraph(i)=i
  110 continue
      igraph(nbeta+1)=inopen+1
c  testausschrift der pfade
C        write(6,*) ' s- haben die positionen: ',
C     &      (igraph(ll),ll=1,nbeta)
c  
c  falls nur eine spineigenfkt vorhanden ist, sprung ans ende
      if (neig.eq.1) goto 9000  
c  
c  aufbau des gesamten graphen durch sukzessives  
c   vertauschen der s- nach hinten
c  
c  icntp: zaehlt alle permutationen  
c  icnts: zaehlt alle permutationen, deren graph oberhalb der achse
      icntp=1
      icnts=1
c  nicht schoen, aber selten (versuch, ausfuehrlich zu sein)
c  beginn der pseudo-schleife
  200 continue
      icntp=icntp+1
c  scannen von igraph und ausfuehren der permutation
c   imerk ist flag, ob schon permutiert wurde
      imerk=0
      do 210 i=1,nbeta
c      falls 'luecke' zwischen zwei s- operatoren, dann vertauschung
c      um eins nach rechts
        if ((igraph(i+1)-igraph(i)).gt.1) then
c        nur dann vertauschen, falls vorher noch nicht permutiert
          if (imerk.eq.0) then
            imerk=1
            igraph(i)=igraph(i)+1
c         die vorhergehenden s- operatoren werden linksbuendig geschrieben
            do 215 j=1,i-1
  215       igraph(j)=j
          endif
        endif
  210 continue
c  testausschrift der pfade
cc      write(6,*) ' s- haben die positionen: ',
cc   &      (igraph(ll),ll=1,nbeta)
c  igraph wird expandiert (ginge auch in schleife 210 mit weniger operationen)
      do 220 i=1,inopen
        ipath(i)=1
  220 continue
      do 230 i=1,nbeta
        ipath(igraph(i))=-1
  230 continue
c  test, ob graph oberhalb der achse liegt
c   jmerk: flag, jmerk=0: graph liegt oberhalb
c   jcnt<0, graph erreicht die untere haelfte
      jmerk=0
      jcnt=0
      do 240 i=inopen,1,-1
        if (jmerk.eq.0) then
          jcnt=jcnt+ipath(i)
          if (jcnt.lt.0) jmerk=1
        endif
  240 continue
c  nur falls jmerk=0, wird pfad abgespeichert
      if (jmerk.eq.0) then
        do 250 i=1,inopen
          idiagr(neig-icnts,i)=ipath(i)
  250   continue
c  testausschrift der pfade
C        write(6,*) ' s- haben die positionen: ',
C     &      (igraph(ll),ll=1,nbeta)
c    # der spineigenfunktionen wird hochgezaehlt
        icnts=icnts+1
      endif
c
c  ausstieg aus pseudo-loop, falls letzte permutation erreicht
      if (igraph(1).eq.(inopen-nbeta+1)) goto 300 
c
c  ruecksprung (gulp)
      goto 200  
c  ende der pseudo-schleife
  300 continue
c  
c  fehlerabfrage 
      if (icnts.ne.neig) then
        write(6,*) ' spineigenfunktionen falsch gezaehlt '  
        write(6,*) ' generiert: ',icnts,' vorhanden: ',neig
        stop 
      endif
c  
 9000 continue
c  
      return                                                                    
      end                                                                       
