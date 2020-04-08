cvp
cvp berechnung einer neuen darstellungsmatrix durch multiplikation
c    mit einem elementaren zyklus
      subroutine cycle1(kml,kmj,itranl,itranr,itrinc)
cvp
cvp parameter und globale common-bloecke sind im include-file
      include "sga.inc"
      include "global"
      real*8 diag1,diag2,ausd1,ausd2
cvp
c  darstellungsmatrizen werden von rechts anmultipliziert
c     do 2040 itran=jjqr-1,jjql,-1
c       do 2050 kk2=1,kml
c         lww=iausd(kk2,itran)
c         diag1=diagtr(kk2,itran)
c         if (lww.gt.kk2) then
c           ausd2=ausdtr(kk2,itran)
c           do 2031 kk1=1,kml
c              usp1(kk1)=umat(kk1,kk2)*diag1
c    &                      +umat(kk1,lww)*ausd2
c              usp2(kk1)=-1.0d0*umat(kk1,lww)*diag1
c    &                      +umat(kk1,kk2)*ausd2
c2031       continue
c           do 2032 kk1=1,kml
c              umat(kk1,kk2)=usp1(kk1)
c              umat(kk1,lww)=usp2(kk1)
c2032       continue
c         endif
c         if (lww.eq.0) then
c           do 2033 kk1=1,kml
c              umat(kk1,kk2)=umat(kk1,kk2)*diag1
c2033       continue
c         endif
c2050   continue
c2040 continue
c  darstellungsmatrizen werden von links anmultipliziert
c     do 2140 itran=jjql,jjqr-1
      do 2140 itran=itranl,itranr,itrinc
c   schleife ueber die verschiedenen spalten von umat1
cc      do 2150 kk2=1,kml
        do 2150 kk2=1,kmj
          do 2160 kk1=1,kml
            usp1(kk1)=umat1(kk1,kk2)
 2160     continue
          do 2170 kk1=1,kml
            diag1=diagtr(kk1,itran)
            usp1(kk1)=usp1(kk1)*diag1
 2170     continue
          do 2180 kk1=1,kml
            lww=iausd(kk1,itran)
            if (lww.gt.kk1) then
              ausd2=ausdtr(lww,itran)
              usp1(kk1)=usp1(kk1)+umat1(lww,kk2)*ausd2
              usp1(lww)=usp1(lww)+umat1(kk1,kk2)*ausd2
            endif
 2180     continue
          do 2190 kk1=1,kml
            umat1(kk1,kk2)=usp1(kk1)
 2190     continue
 2150   continue
 2140 continue
c
      return
      end
