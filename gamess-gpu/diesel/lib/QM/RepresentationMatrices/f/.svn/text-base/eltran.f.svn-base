cvp                                                                             
c   aufbau der felder iausd, diagtr, ausdtr: enthalten die positionen
c    der ausserdiagonalelemente sowie die werte der diagonal- und
c    ausserdiagonalelemente der elementaren transpositionen,
c    hier verwendung der reversed lexical ordering nach
c     duch u. karwowski, comp. phys. rep. vol. 2 (1985) p.122
      subroutine eltran(imult,neig,inopen)
c  
      include"sga.inc"
c  
c  vorbelegung von iausd, diagtr, ausdtr
      do 100 i=1,nopmax-1
        do 100 j=1,mxeigs
          iausd(j,i)=0
          diagtr(j,i)=-1.0d0
          ausdtr(j,i)=0.0d0
  100 continue
c
c  
c  falls nur eine spineigenfkt vorhanden ist, sprung ans ende
      if (neig.eq.1) goto 9000  
c  
c  schleife ueber alle spineigenfunktionen  
      do 1000 i=1,neig
c  schleife ueber die offenen schalen
c    issum: zwischensumme des spins (2*s)
        issum=imult-1-idiagr(i,1)
        do 2000 j=2,inopen
          issum=issum-idiagr(i,j)
c      felder werden nur bei (s+,s-)-kombination gesetzt
c      betrachtung der transposition (j-1,j)
          if (idiagr(i,j-1).eq.1.and.idiagr(i,j).eq.-1) then
c      berechnung der nummer l der wechselwirkenden spineigenfunktion
c        i1 und i2 sind gewichte von vertizes
c        i2 wird null gesetzt, falls vertex nicht im diagramm liegt
            if ((issum-1).gt.0) then
               i2=igew(j-2,issum-1)
             else
               i2=0
            endif
            i1=igew(j-1,issum)
c
            l=i+i1-i2
c        die transposition (j-1,j) wird mit j-1 abgespeichert
            iausd(i,j-1)=l
            iausd(l,j-1)=i
c
cx          xdiag=-1.0d0/(dble(issum+1))
            xdiag=1.0d0/(dble(issum+1))
            xndiag=-1.0d0*sqrt(1.0d0-xdiag*xdiag)
cx          xndiag=sqrt(1.0d0-xdiag*xdiag)
            diagtr(i,j-1)=-1.0d0*xdiag
            diagtr(l,j-1)=xdiag
            ausdtr(i,j-1)=xndiag
            ausdtr(l,j-1)=xndiag
          endif
c      (s-,s+)-kombination, die auf der abszisse liegen, muessen 
c      noch beruecksichtigt werden
          if (issum.eq.0.and.idiagr(i,j-1).eq.-1) then
            diagtr(i,j-1)=1.0d0
          endif
 2000   continue
c
 1000 continue
c  
 9000 continue
c  
      return                                                                    
      end                                                                       
