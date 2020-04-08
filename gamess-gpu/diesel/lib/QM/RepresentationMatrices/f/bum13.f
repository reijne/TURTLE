c
c berechnung der matrixelemente fuer dk=1, p=3 mit hilfe von sga
      subroutine bum13(smint,jjnr,jjqr,kml,kmj,inopen,ibob)
* nahezu identisch zu volker's hmsg13
* macht bislang aber nur die u-mat fuer den cb-teil
* und nicht fuer den ex-teil

cvp
cvp parameter und globale common-bloecke sind im include-file
      include "sga.inc"
      include "global"
c
c ! integralsumme ist bereits in skiny real*8 !
c   integralsumme fuer dk=1,p=3
      real*8 smint
c  
c  wzwei= wurzel zwei
      real*8 wzwei
      wzwei=dsqrt(2.0d0)
C       write(8,*) 'bum13, ibob jjnr',ibob,jjnr
c  
c   berechnung der darstellungsmatrizen mit hilfe von young-tableaux
      do 3010 kk2=1,kml
        do 3020 kk1=1,kml
          umat1(kk1,kk2)=0.0d0
 3020   continue
 3010 continue
c   alle darstellungsmatrizen werden mit wurzel von zwei multipliziert
      do 3025 kk2=1,kml
c       umat1(kk2,kk2)=-1.0d0*wzwei
        umat1(kk2,kk2)=wzwei
 3025 continue
c
c mo-nummern in abhaengigkeit von qr: a<b
* das imoq feld wird in chfeld.f (aufgerufen bei der berechnung
* des referenzraumes) aufgebaut. es enthaelt die information, welche
* positionen der offenen schalen zu beruecksichtigen sind, fuer einen
* gewissen qr, ql fall.
* bedeutung imoq feld, siehe bum01
      monrb=imoq(jjqr,1)
      monra=imoq(jjqr,2)
c
c belegung von umat1 mit u(p0), p0=(b...s)(a...s)
      itranl=inopen-1
      itranr=monra
      itrinc=-1
      call cycle1(kml,kmj,itranl,itranr,itrinc)
c
      itranl=inopen-1
      itranr=monrb
      itrinc=-1
      call cycle1(kml,kmj,itranl,itranr,itrinc)
c
c  
c  belegung von umat mit u(p0)*integralsumme
      do 3050 kk2=1,kmj
        do 3060 kk1=1,kml
          umat(kk1,kk2)=umat1(kk1,kk2)*smint
 3060   continue
 3050 continue
* bisher ist die cb + einelektronen fall abgeschlossen,

* folgender teil geht ueber in den ex-teil, 
* der position und r-fall abhaengig ist

* beginn der berechnung der darstellungsmatrizen fuer
* die exchange integrale
c   schleife ueber alle einzelnen summanden mit zweielektronenintegralen
* iwiev ist anzahl der eintraege fuer eine saf der generated
* iwiev ist anzahl der eintraege fuer eine saf der main
* bedeutung ibob, siehe bum.f
      if(ibob.eq.0) then
*alt       iwiev=kml*(inopen-2)
       iwiev=kmj*(inopen-2)
      else if (ibob.eq.1) then
*alt       iwiev=kmj*(inopen-2)
       iwiev=kml*(inopen-2)
      else 
       write(6,*) 'fehler in ibob fuer bum13'
      endif
C      write(6,*) 'hallo',iwiev
* iposi, zaehler fuer position
      iposi=0
* nposi, anzahl der summanden fuer offene schalen
      nposi=inopen-2
c    iis: zaehler fuer die permutationen p=(ai)p0
      iis=1
      do 3200 iit=1,inopen
*      write(6,*) 'hallo 3200'
c  
c permutation wird nur fuer i =/= a,b berechnet  
      if (iit.ne.monra.and.iit.ne.monrb) then
c      zaehler iis fuer die permutationen laeuft also von 2 bis inopen-1
C      write(8,*) 'hallo, iit.ne.',iit,monra,monrb
        iis=iis+1
        iposi=iposi+1
C      write(8,*) 'iposi',iposi
c  
c  initialisierung von umat2 mit u(p0)
      do 3110 kk2=1,kmj
        do 3120 kk1=1,kml
          umat2(kk1,kk2)=umat1(kk1,kk2)
 3120   continue
 3110 continue
c  
c fuer r=1  
      if (jjnr.eq.1) then
*      write(6,*) 'hallo jjnr.eq.1'
c berechnung von p=(bi)p0
      if (iit.lt.monrb) then
        itranl=monrb-1
        itranr=iit
        itrinc=-1
        call cycle2(kml,kmj,itranl,itranr,itrinc)
        itranl=iit+1
        itranr=monrb-1
        itrinc=1
        call cycle2(kml,kmj,itranl,itranr,itrinc)
      endif
      if (iit.gt.monrb) then
        itranl=monrb
        itranr=iit-1
        itrinc=1
        call cycle2(kml,kmj,itranl,itranr,itrinc)
        itranl=iit-2
        itranr=monrb
        itrinc=-1
        call cycle2(kml,kmj,itranl,itranr,itrinc)
      endif
c
*      do 3130 kk2=1,kmj
*        do 3140 kk1=1,kml
c  hier -umat, da zweielektronenintregale negativ genommen werden
*          umat(kk1,kk2)=umat(kk1,kk2)-umat2(kk1,kk2)*sac(iis-1)
* 3140   continue
* 3130 continue
* abgespeichert wird umat, damit man spaeter einfach
* aufaddiert
       if(ibob.eq.0) then
*       if(ibob.eq.1) then
C      write(6,*) 'hallo vor umspeichern ibob 0 '
C      do kk1=1,kml
C       write(6,*) (umat2(kk1,kk2),kk2=1,kmj)
C      enddo
*        do kk2=1,kmj
*         do kk1=1,kml
         do kk1=1,kml
        do kk2=1,kmj
*alt          iem=iposi+nposi*(kk1-1)+(kk2-1)*iwiev
          iem=iposi+nposi*(kk2-1)+(kk1-1)*iwiev
          emat(iem)=umat2(kk1,kk2)
C          write(8,*) 'iem 1',iem,emat(iem),iposi,nposi,iwiev
         enddo
        enddo
       else if (ibob.eq.1) then
*       else if (ibob.eq.0) then
C      write(6,*) 'hallo vor umspeichern ibob 1 '
C      do kk1=1,kml
C       write(6,*) (umat2(kk1,kk2),kk2=1,kmj)
C      enddo
*        do kk2=1,kmj
*         do kk1=1,kml
         do kk1=1,kml
        do kk2=1,kmj
*alt          iem=iposi+nposi*(kk2-1)+(kk1-1)*iwiev
          iem=iposi+nposi*(kk1-1)+(kk2-1)*iwiev
          emat(iem)=umat2(kk1,kk2)
C          write(8,*) 'iem 2',iem,emat(iem),iposi,nposi,iwiev
         enddo
        enddo
       endif
c  ende r=1
      endif
c
c  
c fuer r=2  
      if (jjnr.eq.2) then
c berechnung von p=(ai)p0
      if (iit.lt.monra) then
        itranl=monra-1
        itranr=iit
        itrinc=-1
        call cycle2(kml,kmj,itranl,itranr,itrinc)
        itranl=iit+1
        itranr=monra-1
        itrinc=1
        call cycle2(kml,kmj,itranl,itranr,itrinc)
      endif
      if (iit.gt.monra) then
        itranl=monra
        itranr=iit-1
        itrinc=1
        call cycle2(kml,kmj,itranl,itranr,itrinc)
        itranl=iit-2
        itranr=monra
        itrinc=-1
        call cycle2(kml,kmj,itranl,itranr,itrinc)
      endif
c
*      do 3230 kk2=1,kmj
*        do 3240 kk1=1,kml
c  hier -umat, da zweielektronenintregale negativ genommen werden
*          umat(kk1,kk2)=umat(kk1,kk2)-umat2(kk1,kk2)*sac(iis-1)
* 3240   continue
* 3230 continue
* abgespeichert wird umat, damit man spaeter einfach
* aufaddiert
       if(ibob.eq.0) then
*alt        do kk2=1,kmj
*alt         do kk1=1,kml
         do kk1=1,kml
        do kk2=1,kmj
*alt          iem=iposi+nposi*(kk1-1)+(kk2-1)*iwiev
          iem=iposi+nposi*(kk2-1)+(kk1-1)*iwiev
          emat(iem)=umat2(kk1,kk2)
C          write(8,*) 'iem',iem,emat(iem)
         enddo
        enddo
       else if (ibob.eq.1) then
*alt        do kk2=1,kmj
*alt         do kk1=1,kml
         do kk1=1,kml
        do kk2=1,kmj
*alt          iem=iposi+nposi*(kk2-1)+(kk1-1)*iwiev
          iem=iposi+nposi*(kk1-1)+(kk2-1)*iwiev
          emat(iem)=umat2(kk1,kk2)
C          write(8,*) 'iem',iem,emat(iem)
         enddo
        enddo
       endif
c
c  ende r=2
      endif
c
* ausdruck der matrix zum testen
*       write(6,*) 'umat2 probeausdruck',iposi
*       do kk1=1,kml
*        write(6,993) (umat2(kk1,iraus),iraus=1,kmj)
*       enddo
      endif
c
 3200 continue
c  
c umat enthaelt die ergebnismatrix
c  
993   format(10(F8.4))
      return
      end
