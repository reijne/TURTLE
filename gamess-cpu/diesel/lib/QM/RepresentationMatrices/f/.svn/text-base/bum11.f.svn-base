c
c berechnung der matrixelemente fuer dk=1, p=1 mit hilfe von sga
      subroutine bum11(cbint,exint,jjnr,jjqr,jjql,kml,kmj,inopen)
* nahezu identisch zu volkers, einzige ausnahme ist, dass die
* u mat nicht gleich mit dem echten cb bzw ex vermatscht wird
* b.e. 19.01.95
cvp
cvp parameter und globale common-bloecke sind im include-file
      include "sga.inc"
      include "global"
c
c ! integrale sind bereits in skiny real*8 !
c   integrale fuer dk=1, p=1
      real*8 cbint,exint
c  
c  wzwei= wurzel zwei
      real*8 wzwei
      wzwei=dsqrt(2.0d0)
c  
c
c   berechnung der darstellungsmatrizen mit hilfe von young-tableaux
      do 1010 kk2=1,kml
        do 1020 kk1=1,kml
          umat(kk1,kk2)=0.0d0
 1020   continue
 1010 continue
      do 1025 kk2=1,kml
        umat(kk2,kk2)=-1.0d0*wzwei
 1025 continue
c
c mo-nummern in abhaengigkeit von qr: a<b<c, ql gibt mo-nr. d
* das imo3q feld wird in chfeld.f (aufgerufen bei der berechnung
* des referenzraumes) aufgebaut. es enthaelt die information, welche
* positionen der offenen schalen zu beruecksichtigen sind, fuer einen
* gewissen qr, ql fall.

      monrd=jjql
      monrc=imo3q(jjqr,1)
      monrb=imo3q(jjqr,2)
      monra=imo3q(jjqr,3)
c
c
c    p0=(c...s)(b...s)(a...d_l), s=inopen
      if (monrd.ge.monra) then
        itranl=monrd-1
        itranr=monra
        itrinc=-1
       else
        itranl=monrd
        itranr=monra-1
        itrinc=1
      endif
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
      itranl=inopen-1
      itranr=monrb
      itrinc=-1
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
      itranl=inopen-1
      itranr=monrc
      itrinc=-1
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
      do 1030 kk2=1,kml
        do 1040 kk1=1,kml
          umat1(kk1,kk2)=umat(kk1,kk2)*exint
 1040   continue
 1030 continue
      do 1050 kk2=1,kml
        do 1060 kk1=1,kml
          umat(kk1,kk2)=umat(kk1,kk2)*cbint
 1060   continue
 1050 continue
c
c   darstellungsmatrizen fuer r=1
c    pa=(ac)p0, pb=(ab)p0
c
      if (jjnr.eq.1) then
c
c
      itranl=monrc-1
      itranr=monra
      itrinc=-1
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
      itranl=monra+1
      itranr=monrc-1
      itrinc=1
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
      itranl=monrb-1
      itranr=monra
      itrinc=-1
      call cycle1(kml,kmj,itranl,itranr,itrinc)
c
      itranl=monra+1
      itranr=monrb-1
      itrinc=1
      call cycle1(kml,kmj,itranl,itranr,itrinc)
c
cbe      do 1080 kk2=1,kmj
cbe        do 1090 kk1=1,kml
c     hier -umat1, da die integrale mit -1 multipliziert werden
cbe          umat(kk1,kk2)=umat(kk1,kk2)-umat1(kk1,kk2)
cbe 1090   continue
cbe 1080 continue
c
      endif
c
c   darstellungsmatrizen fuer r=2
c    pa=(ac)p0, pb=p0
c
      if (jjnr.eq.2) then
c
c
      itranl=monrc-1
      itranr=monra
      itrinc=-1
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
      itranl=monra+1
      itranr=monrc-1
      itrinc=1
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
cbe      do 1180 kk2=1,kmj
cbe        do 1190 kk1=1,kml
c     hier -umat1, da die integrale mit -1 multipliziert werden
cbe          umat(kk1,kk2)=umat(kk1,kk2)-umat1(kk1,kk2)
cbe 1190   continue
cbe 1180 continue
c
      endif
c
c   darstellungsmatrizen fuer r=3
c    pa=(ab)p0, pb=p0
c
      if (jjnr.eq.3) then
c
c
      itranl=monrb-1
      itranr=monra
      itrinc=-1
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
      itranl=monra+1
      itranr=monrb-1
      itrinc=1
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
cbe      do 1280 kk2=1,kmj
cbe        do 1290 kk1=1,kml
c     hier -umat1, da die integrale mit -1 multipliziert werden
cbe          umat(kk1,kk2)=umat(kk1,kk2)-umat1(kk1,kk2)
cbe 1290   continue
cbe 1280 continue
c
      endif
c
c
      return
      end
