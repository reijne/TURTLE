c
c berechnung der matrixelemente fuer dk=0, p=1 mit hilfe von sga
      subroutine bum04(kml)
* einfachster fall, enthaelt nur die diagonale
* b.e. 23.01.95
cvp
cvp parameter und globale common-bloecke sind im include-file
      include "sga.inc"
      include "global"
c
c ! integrale sind bereits in skiny real*8 !
c   integrale fuer dk=0, p=1
      real*8 cbint,exint
c  
c
      do 1025 kk2=1,kml
        umat(kk2,kk2)=1.0d0
        umat1(kk2,kk2)=1.0d0
 1025 continue
      return
      end
