c
c berechnung der u matrix fuer dk=0, p=2 mit hilfe von sga
      subroutine bum02(smint,jjqr,jjql,kml,ibob,inopen)
* nahezu identisch zu volkers hmsg02
cvp
cvp parameter und globale common-bloecke sind im include-file
      include "sga.inc"
      include "global"
c
c ! integral ist bereits in skiny real*8 !
c   integral fuer dk=0, p=2
      real*8 smint
c  
c   berechnung der darstellungsmatrizen mit hilfe von young-tableaux
      do 2010 kk2=1,kml
        do 2020 kk1=1,kml
          umat(kk1,kk2)=0.0d0
 2020   continue
 2010 continue
      do 2025 kk2=1,kml
        umat(kk2,kk2)=-1.0d0*smint
 2025 continue
c
cf fallunterscheidung, ob qr >= ql oder qr < ql wird hier bei aufbau
cf  der darstellungsmatrix getroffen
cf    if (jjqr.ge.jjql) then
cf      itranl=jjql
cf      itranr=jjqr-1
cf      itrinc=1
cf     else
cf      itranl=jjql-1
cf      itranr=jjqr
cf      itrinc=-1
cf    endif
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
c    p0=(a...c), hier ist wegen fallunterscheidung a > c
      itranl=monrc
      itranr=monra-1
      itrinc=1
      call cycle(kml,kml,itranl,itranr,itrinc)
c
c  
      return
      end
