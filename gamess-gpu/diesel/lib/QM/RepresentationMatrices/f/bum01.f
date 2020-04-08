c
c berechnung der matrixelemente fuer dk=0, p=1 mit hilfe von sga
      subroutine bum01(cbint,exint,jjnr,jjqr,jjql,kml,ibob,inopen)
* nahezu identisch zu volkers hmsg01, einzige ausnahme ist, dass die
* u mat nicht gleich mit dem echten cb bzw ex vermatscht wird
* b.e. 19.01.95
cvp
cvp parameter und globale common-bloecke sind im include-file
      include "sga.inc"
      include "global"
c
c ! integrale sind bereits in skiny real*8 !
c   integrale fuer dk=0, p=1
      real*8 cbint,exint
C      write(8,*) 'aufruf bum01'
C      write(8,*) 'cbint,exint,jjnr,jjqr,jjql,kml,ibob,inopen'
C      write(8,*) cbint,exint,jjnr,jjqr,jjql,kml,ibob,inopen
C      write(8,*)'p=1, r=',jjnr,'ql=',jjql,'qr=',jjqr
c  
c
c   berechnung der darstellungsmatrizen mit hilfe von young-tableaux
      do 1010 kk2=1,kml
        do 1020 kk1=1,kml
          umat(kk1,kk2)=0.0d0
 1020   continue
 1010 continue
      do 1025 kk2=1,kml
        umat(kk2,kk2)=1.0d0
 1025 continue
c
c fallunterscheidung, ob qr >= ql oder qr < ql wird hier bei aufbau
c  wird wie bisher ueber ibob aufgebaut. 
c  hierbei muessen die q-faelle der groesse nach sortiert werden.
c  alternativ koennten die
c  inversen darstellungsmatrizen berechnet werden (entspricht vertauschung
c  der konfigurationen), ohne die q-faelle zu sortieren.
* das imoq feld wird in chfeld.f (aufgerufen bei der berechnung
* des referenzraumes) aufgebaut. es enthaelt die information, welche
* positionen der offenen schalen zu beruecksichtigen sind, fuer einen
* gewissen qr, ql fall.
      if (jjqr.ge.jjql) then
        ibob=0
        monra=imoq(jjqr,1)
        monrb=imoq(jjqr,2)
        monrc=imoq(jjql,1)
        monrd=imoq(jjql,2)
       else
        ibob=1
        monra=imoq(jjql,1)
        monrb=imoq(jjql,2)
        monrc=imoq(jjqr,1)
        monrd=imoq(jjqr,2)
      endif
c
c   belegung von umat mit u(p0)
c    p0=(b...d)(a...c) mit a>=c
      itranl=monrc
      itranr=monra-1
      itrinc=1
      call cycle(kml,kml,itranl,itranr,itrinc)
c
      if (monrb.ge.monrd) then
        itranl=monrd
        itranr=monrb-1
        itrinc=1
       else
        itranl=monrd-1
        itranr=monrb
        itrinc=-1
      endif
      call cycle(kml,kml,itranl,itranr,itrinc)
c
c   belegung von umat1 mit u(p1)
c    p1=(a,b)pa mit a>=b
c  
      do 1060 kk2=1,kml
        do 1070 kk1=1,kml
          umat1(kk1,kk2)=umat(kk1,kk2)
 1070   continue
 1060 continue
c  
      itranl=monra-1
      itranr=monrb
      itrinc=-1
      call cycle1(kml,kml,itranl,itranr,itrinc)
c
      itranl=monrb+1
      itranr=monra-1
      itrinc=1
      call cycle1(kml,kml,itranl,itranr,itrinc)
c
c   ergebnismatrix fuer r=1 auf umat
c    pa=p0, pb=p1
c
      if (jjnr.eq.1) then
c
c   hier -umat1, da die ex-integrale mit -1 multipliziert werden
cbe      do 1180 kk2=1,kml
cbe        do 1190 kk1=1,kml
cbe          umat(kk1,kk2)=umat(kk1,kk2)*cbint-umat1(kk1,kk2)*exint
cbe 1190   continue
cbe 1180 continue
cbe beim rausschreiben von umat1 muss man hier umat1 und
cbe nicht -umat1 herausgeben, da in bum.f das exint auf 1.0
cbe und nicht wie im skiny auf -1.0d0 gesetzt wird
c  
      do kk2=1,kml
       do kk1=1,kml
        umat(kk1,kk2)=umat(kk1,kk2)*cbint
        umat1(kk1,kk2)=umat1(kk1,kk2)*exint
       enddo
      enddo
      endif
c
c   ergebnismatrix fuer r=2 auf umat
c    pa=p0+p1, pb=-p0
c
      if (jjnr.eq.2) then
c
c   
cbe      do 1280 kk2=1,kml
cbe        do 1290 kk1=1,kml
cbe          umat(kk1,kk2)=(umat(kk1,kk2)+umat1(kk1,kk2))*cbint
cbe     &                  +umat(kk1,kk2)*exint
cbe 1290   continue
cbe 1280 continue
      do kk2=1,kml
        do kk1=1,kml
          umat2(kk1,kk2)=(umat(kk1,kk2)+umat1(kk1,kk2))*cbint
        enddo
      enddo
      do kk2=1,kml
        do kk1=1,kml
          umat1(kk1,kk2)=-umat(kk1,kk2)*exint
        enddo
      enddo
      do kk2=1,kml
        do kk1=1,kml
          umat(kk1,kk2)=umat2(kk1,kk2)
        enddo
      enddo
      
c  
      endif
c
c   ergebnismatrix fuer r=3 auf umat
c    pa=p0+p1, pb=-p1
c
      if (jjnr.eq.3) then
c
c   
cbe      do 1380 kk2=1,kml
cbe        do 1390 kk1=1,kml
cbe          umat(kk1,kk2)=(umat(kk1,kk2)+umat1(kk1,kk2))*cbint
cbe     &                  +umat1(kk1,kk2)*exint
cbe 1390   continue
cbe 1380 continue
      do kk2=1,kml
        do kk1=1,kml
          umat(kk1,kk2)=(umat(kk1,kk2)+umat1(kk1,kk2))*cbint
        enddo
      enddo
      do kk2=1,kml
        do kk1=1,kml
          umat1(kk1,kk2)=-umat1(kk1,kk2)*exint
        enddo
      enddo
c  
      endif
c  
      return
      end
