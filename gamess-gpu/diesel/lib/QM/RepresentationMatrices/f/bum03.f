c
c berechnung der matrixelemente fuer dk=0, p=3 mit hilfe von sga
      subroutine bum03(smint,jjnr,jjqr,jjql,kml,ibob,inopen)
* nahezu identisch zu volker's hmsg03
* macht bislang aber nur die u-mat fuer den cb-teil
* und nicht fuer den ex-teil
* ibob  ....
* fuer dk=0: ibob=0 qr>=ql umat(kk1,kk2)  kk1 mit ql (also gen)
*                                         kk2 mit qr (also main)
*            ibob=1 qr< ql umat(kk1,kk2)  kk1 mit qr (also main)
*                                         kk2 mit ql (also gen)
cvp
cvp parameter und globale common-bloecke sind im include-file
      include "sga.inc"
      include "global"
c
c ! integralsumme ist bereits in skiny real*8 !
c   integralsumme fuer dk=0,p=3
      real*8 smint
c  
c
      do 3010 kk2=1,kml
        do 3020 kk1=1,kml
          umat1(kk1,kk2)=0.0d0
 3020   continue
 3010 continue
      do 3025 kk2=1,kml
        umat1(kk2,kk2)=1.0d0
 3025 continue
c fallunterscheidung, ob qr >= ql oder qr < ql wird hier
c  wie bisher ueber ibob aufgebaut.
c  hierbei muessen die q-faelle der groesse nach sortiert werden.
c  alternativ koennten die
c  inversen darstellungsmatrizen berechnet werden (entspricht vertauschung
c  der konfigurationen), ohne die q-faelle zu sortieren.

      if (jjqr.ge.jjql) then
        ibob=0
        monra=jjqr
        monrc=jjql
       else
        ibob=1
        monra=jjql
        monrc=jjqr
      endif
c
c    p0=(a...c), hier ist wegen fallunterscheidung a > c
c    umat1 enthaelt permutation fuer einelektronenintegral
      itranl=monrc
      itranr=monra-1
      itrinc=1
      call cycle1(kml,kml,itranl,itranr,itrinc)
c
c  setzen der einelektronenterme
c  r=1
      if (jjnr.eq.1) then
c
      do 3110 kk2=1,kml
        do 3120 kk1=1,kml
          umat(kk1,kk2)=umat1(kk1,kk2)*smint
 3120   continue
 3110 continue
c
      endif
c
c  r=2
      if (jjnr.eq.2) then
c
      do 3210 kk2=1,kml
        do 3220 kk1=1,kml
          umat(kk1,kk2)=-umat1(kk1,kk2)*smint
 3220   continue
 3210 continue
c
      endif
c
* bisher ist die cb + einelektronen fall abgeschlossen,
* folgender teil geht ueber in den ex-teil, 
* der position und r-fall abhaengig ist


* iwiev ist anzahl der eintraege fuer eine saf der generated
       iwiev=kml*(inopen-1)
* iposi, zaehler fuer position
      iposi=0
* nposi, anzahl der summanden fuer offene schalen
      nposi=inopen-1
c
c   schleife ueber alle einzelnen summanden
      do 3400 iis=2,inopen
* hochzaehlen fuer zaehler fuer position
       iposi=iposi+1
c
c   initialisierung von umat2
      do 3410 kk2=1,kml
        do 3410 kk1=1,kml
          umat2(kk1,kk2)=umat1(kk1,kk2)
 3410 continue
c
c   umat2 enthaelt p1=(ai)p0
      if (iis.le.monra) then
        itranl=monra-1
        itranr=iis-1
        itrinc=-1
        call cycle2(kml,kml,itranl,itranr,itrinc)
        itranl=iis
        itranr=monra-1
        itrinc=1
        call cycle2(kml,kml,itranl,itranr,itrinc)
       else
        itranl=monra
        itranr=iis-1
        itrinc=1
        call cycle2(kml,kml,itranl,itranr,itrinc)
        itranl=iis-2
        itranr=monra
        itrinc=-1
        call cycle2(kml,kml,itranl,itranr,itrinc)
      endif
c
c  r=1
      if (jjnr.eq.1) then
c
*      do 3161 kk2=1,kml
*        do 3171 kk1=1,kml
c
c   zweielektronenterme: u(p1), hier -umat2 wg. vorzeichen von sac
*          umat(kk1,kk2)=umat(kk1,kk2)-umat2(kk1,kk2)*sac(iis-1)
* 3171   continue
* 3161 continue
*       if(ibob.eq.0) then
       if(ibob.eq.1) then
*      write(6,*) 'hallo vor umspeichern ibob 0 '
        do kk2=1,kml
         do kk1=1,kml
          iem=iposi+nposi*(kk2-1)+(kk1-1)*iwiev
          emat(iem)=umat2(kk1,kk2)
*          write(6,*) 'iem',iem
         enddo
        enddo
*       else if (ibob.eq.1) then
       else if (ibob.eq.0) then
*      write(6,*) 'hallo vor umspeichern ibob 1 '
        do kk2=1,kml
         do kk1=1,kml
          iem=iposi+nposi*(kk1-1)+(kk2-1)*iwiev
          emat(iem)=umat2(kk1,kk2)
*          write(6,*) 'iem',iem
         enddo
        enddo
       else
        write(6,*) 'fehler in ibob, siehe bum03'
       endif
c
      endif
c
c  r=2
      if (jjnr.eq.2) then
c
*      do 3261 kk2=1,kml
*        do 3271 kk1=1,kml
c
c   zweielektronenterme: u(p0+p1), hier -umatx wg. vorzeichen von sac
*          umat(kk1,kk2)=umat(kk1,kk2)
*     &                -(umat1(kk1,kk2)+umat2(kk1,kk2))*sac(iis-1)
* 3271   continue
* 3261 continue
       do kk2=1,kml
        do kk1=1,kml
         umat2(kk1,kk2)=umat1(kk1,kk2)+umat2(kk1,kk2)
        enddo
       enddo
 
*       if(ibob.eq.0) then
       if(ibob.eq.1) then
*      write(6,*) 'hallo vor umspeichern ibob 0 '
        do kk2=1,kml
         do kk1=1,kml
          iem=iposi+nposi*(kk2-1)+(kk1-1)*iwiev
          emat(iem)=umat2(kk1,kk2)
*          write(6,*) 'iem',iem
         enddo
        enddo
*       else if (ibob.eq.1) then
       else if (ibob.eq.0) then
*      write(6,*) 'hallo vor umspeichern ibob 1 '
        do kk2=1,kml
         do kk1=1,kml
          iem=iposi+nposi*(kk1-1)+(kk2-1)*iwiev
          emat(iem)=umat2(kk1,kk2)
*          write(6,*) 'iem',iem
         enddo
        enddo
       else
        write(6,*) 'fehler in ibob, siehe bum03'
       endif
c
      endif
* ausdruck der matrix zum testen
*       write(6,*) 'umat2 probeausdruck',iposi
*       if(ibob.eq.0) then
*        write(6,*) 'vorne gen, hinten main'
*       else
*        write(6,*) 'vorne main, hinten gen'
*       endif
*       do kk1=1,kml
*        write(6,993) (umat2(kk1,iraus),iraus=1,kml)
*       enddo
c
 3400 continue
c
c
993   format(10(F8.4))
      return
      end
