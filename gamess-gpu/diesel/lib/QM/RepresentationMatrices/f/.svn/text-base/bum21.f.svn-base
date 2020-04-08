c
c berechnung der matrixelemente fuer dk=2 mit hilfe von sga
      subroutine bum21(cbint,exint,jjnr,jjqr,kml,kmj,inopen)
* nahezu identisch zu hmsg21 von volker
cvp
cvp parameter und globale common-bloecke sind im include-file
      include "sga.inc"
      include "global"
c
c ! integrale sind bereits in skiny real*8 !
c   integrale fuer dk=2
      real*8 cbint,exint
c  
c  
c
c   berechnung der darstellungsmatrizen mit hilfe von young-tableaux
      do 1010 kk2=1,kml
        do 1020 kk1=1,kml
          umat(kk1,kk2)=0.0d0
 1020   continue
 1010 continue
      do 1025 kk2=1,kml
        umat(kk2,kk2)=2.0d0
 1025 continue
c
c mo-nummern in abhaengigkeit von qr: a<b<c<d
      monrd=imo4q(jjqr,1)
      monrc=imo4q(jjqr,2)
      monrb=imo4q(jjqr,3)
      monra=imo4q(jjqr,4)
c
c    p0=(d...s)(c...s)(b...s)(a...s), s=inopen
      itranl=inopen-1
      itranr=monra
      itrinc=-1
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
      itranl=inopen-1
      itranr=monrd
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
c    pa=(bc)p0, pb=(bd)p0
c
      if (jjnr.eq.1) then
c
c
      itranl=monrc-1
      itranr=monrb
      itrinc=-1
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
      itranl=monrb+1
      itranr=monrc-1
      itrinc=1
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
      itranl=monrd-1
      itranr=monrb
      itrinc=-1
      call cycle1(kml,kmj,itranl,itranr,itrinc)
c
      itranl=monrb+1
      itranr=monrd-1
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
c    pa=p0, pb=(bd)p0
c
      if (jjnr.eq.2) then
c
c
      itranl=monrd-1
      itranr=monrb
      itrinc=-1
      call cycle1(kml,kmj,itranl,itranr,itrinc)
c
      itranl=monrb+1
      itranr=monrd-1
      itrinc=1
      call cycle1(kml,kmj,itranl,itranr,itrinc)
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
c    pa=p0, pb=(bc)p0
c
      if (jjnr.eq.3) then
c
c
      itranl=monrc-1
      itranr=monrb
      itrinc=-1
      call cycle1(kml,kmj,itranl,itranr,itrinc)
c
      itranl=monrb+1
      itranr=monrc-1
      itrinc=1
      call cycle1(kml,kmj,itranl,itranr,itrinc)
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
