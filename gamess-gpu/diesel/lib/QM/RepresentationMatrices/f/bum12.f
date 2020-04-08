c
c berechnung der matrixelemente fuer dk=1, p=2 mit hilfe von sga
      subroutine bum12(cbint,jjqr,kml,kmj,inopen)
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
c   berechnung der darstellungsmatrizen mit hilfe von young-tableaux
      do 2010 kk2=1,kml
        do 2020 kk1=1,kml
          umat(kk1,kk2)=0.0d0
 2020   continue
 2010 continue
c   alle darstellungsmatrizen werden mit wurzel von zwei multipliziert
      do 2025 kk2=1,kml
c       umat(kk2,kk2)=-1.0d0*wzwei
        umat(kk2,kk2)=wzwei*cbint
 2025 continue
c
c mo-nummern in abhaengigkeit von qr: a<b
* das imoq feld wird in chfeld.f (aufgerufen bei der berechnung
* des referenzraumes) aufgebaut. es enthaelt die information, welche
* positionen der offenen schalen zu beruecksichtigen sind, fuer einen
* gewissen qr, ql fall.
      monrb=imoq(jjqr,1)
      monra=imoq(jjqr,2)
c
      itranl=inopen-1
      itranr=monra
      itrinc=-1
      call cycle(kml,kmj,itranl,itranr,itrinc)
c
      itranl=inopen-1
      itranr=monrb
      itrinc=-1
      call cycle(kml,kmj,itranl,itranr,itrinc)

      return
      end
