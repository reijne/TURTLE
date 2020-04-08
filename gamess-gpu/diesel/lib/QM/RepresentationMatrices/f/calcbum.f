      subroutine calcbum(idsk,ipfa,irfa1,iqrfa,iqlfa,kml,
     *               kmj,inopen,ibob)
* module_id   : parkeu.x            author      : b.e.
* generation  : 19.01.96            update      :
* calls to    : bumxy               called by   : various programs
* description

* bum : berechnet u matrix
* steuert den aufruf routinen bumxy, die die
* darstellungsmatrizen fuer die einzelnen faelle
* berechnen
* bumxy, x steht fuer dk fall, y steht fuer r fall

* ACHTUNG ACHTUNG ACHTUNG ACHTUNG ACHTUNG ACHTUNG
* die ansteuerung geschieht so, dass die cb und ex teile
* addiert werden muessen!!!!!!!!!!!!
* u1*cb + u2*ex

* uebergabeparameter
*====================
* idsk ..... betrachteter dk fall
* ipfa ..... betrachteter p fall
* irfa ..... betrachteter R fall
* iqrfa .... betrachteter qr fall, 
*            dk=/=0 mit hoehere sk verknuepft
*            dk = 0 mit main (=festgehaltener) verknuepft
* iqlfa .... betrachteter ql fall, 
*            dk=/=0 mit niedere sk verknuepft
*            dk = 0 mit gen (=laufender) verknuepft
* kml   .... anzahl der saf
* kmj   .... anzahl der saf
*              kml > kmj
* inopen.... anzahl offener schalen der groesseren sk
* ibob  .... 
* fuer dk=1 spielt ibob nur fuer die exch der summe eine
* rolle, damit man sie richtig ins emat feld einsortiert
* ibob regelt ob gen, oder main der hoeheren angehoert
* ibob = 0: gen gehoert tieferen an, steht also hinten
* ibob = 1: gen gehoert hoeheren an, steht also vorne
* fuer dk=0: ibob=0 umat(kk1,kk2)  kk1 mit ql (also gen)
*                                  kk2 mit qr (also main)
*            ibob=1 umat(kk1,kk2)  kk1 mit qr (also main)
*                                  kk2 mit ql (also gen)
* regelung fuer die umat's:
* dk =/= 0, klare sache, die umat's sind recheckmatrix:
* zeile stehe fuer hohere sk, spalte fuer niedere
*
* dk=0, qr >= ql ==> ibob = 0: umat(kk1,kk2) kk1 mit ql (= gen)
*                                            kk2 mit qr (= main)
* dk=0, qr < ql  ==> ibob = 1: umat(kk1,kk2) kk1 mit qr (= main)
*                                            kk2 mit ql (= gen)
* 
* alle fall bezeichnungen beziehen sich auf Tabelle 3.2 der
* Diss von V.Ple\ss\


      include "sga.inc"
      include "global"


* trace option
C      if(itrace.eq.1) then
C        write(iftrac,*) 'subroutine bum called'
C      endif
*        write(jtape,*)
C        write(6,*) '***************************'
C        write(6,*) 'Informationen aus dem bum  '
C        write(6,*) 'bum:idsk,ipfa,irfa,iqrfa,iqlfa,kml,kmj,inopen'
C        write(6,*) idsk,ipfa,irfa,iqrfa,iqlfa,kml,kmj,inopen

* umdefinition von irfa
      irfa=irfa1

      do kk1=1,kml
       do kk2=1,kmj
        umat(kk1,kk2)=0
       enddo
      enddo


* dk=0 fall
      if(idsk.eq.0) then
       if(ipfa.eq.1) then
        cbint=1.0d0
        exint=1.0d0
        if(irfa.eq.6) irfa=1
        if(irfa.eq.5) irfa=2
        if(irfa.eq.4) irfa=3
        call bum01(cbint,exint,irfa,iqrfa,iqlfa,kml,ibob,inopen)
*        write(6,*) 'ibob',ibob
       else if (ipfa.eq.2) then
        cbint=1.0d0
        call bum02(cbint,iqrfa,iqlfa,kml,ibob,inopen)
*        write(6,*) 'ibob',ibob
       else if (ipfa.eq.3) then
        sm=1.0d0
          iwiev=(inopen-1)*kml*kml
          do ilauf=1,iwiev
           emat(ilauf)=0.0d0
          enddo
        call bum03(sm,irfa,iqrfa,iqlfa,kml,ibob,inopen)
* ausprinten von emat
        iwiev1=kml*kml
        ianf=0
*        write(6,*) 'emat probeausdruck'
*        do ilauf=1,iwiev1
*         iend=ianf+inopen-1
*         ianf=ianf+1
*         write(6,993) (emat(iraus),iraus=ianf,iend)
*         ianf=iend
*        enddo
* ausprinten von emat zuende
*        write(6,*) 'ibob',ibob
       else if (ipfa.eq.4) then
        call bum04(kml)
       else if (ipfa.eq.5) then
        write(6,*) 'tut noch nicht'
       else
        write(6,*) 'Fehlerhafter P-Fall fuer dk=0, siehe bum'
       endif
      endif

* dk=1 fall
      if(idsk.eq.1) then
        if(ipfa.eq.1) then
          if(irfa.ge.4) then
* codierung der r-faelle 4-6 durch minus zeichen fuer cb und ex
            cbint=-1.0d0
            exint=-1.0d0
            irfa=irfa-3
          else if(irfa.lt.4) then
            cbint= 1.0d0
            exint=1.0d0
          else
            write(6,*) 'fehler im R-fall, siehe bum'
          endif
*          write(6,*) 'cbint,exint',cbint,exint
          call bum11(cbint,exint,irfa,iqrfa,iqlfa,kml,kmj,inopen)

* p=2
        else if (ipfa.eq.2) then
          if(irfa.eq.2) then
* codierung der r-fall 2 durch minus zeichen fuer cb
            cbint=1.0d0
            call bum12(cbint,iqrfa,kml,kmj,inopen)
* codierung der r-fall 1 durch plus zeichen fuer cb
          else if(irfa.eq.1) then
            cbint=-1.0d0
            call bum12(cbint,iqrfa,kml,kmj,inopen)
          else
            write(6,*) 'fehler im R-fall, siehe bum'
          endif
* p=3
        else if (ipfa.eq.3) then
*          write(6,*) 'umat enthaelt die u fuer den cb teil'
          cbint=1.0d0
          exint=1.0d0
* fuer dk=1 spielt ibob nur fuer die exch der summe eine
* rolle, damit man sie richtig ins emat feld einsortiert
* ibob regelt ob gen, oder main der hoeheren angehoert
* ibob = 0: gen gehoert tieferen an, steht also hinten
* ibob = 1: gen gehoert hoeheren an, steht also vorne
* der versuch, es genauso zu machen wie fuer dk=0
* ist leider gescheitert, so'n schei..
* also hier ist es umgekehrt, (wenn man die main als fest)
* bezeichnet
* hier test
CCCCCCCCCCCCCCCCCCCCCC
C!!!!!!!!!!!!!!!!!!!!!
C          ibob=0
CCCCCCCCCCCCCCCCCCCCCC
          iwiev=(inopen-2)*kml*kmj
          do ilauf=1,iwiev
           emat(ilauf)=0.0d0
          enddo
          call bum13(cbint,irfa,iqrfa,kml,kmj,inopen,ibob)
          iwiev1=kml*kmj
          ianf=0
*          write(6,*) 'emat probeausdruck'
*          do ilauf=1,iwiev1
*           iend=ianf+inopen-2
*           ianf=ianf+1
*           write(6,*) ianf,iend
*           write(6,993) (emat(iraus),iraus=ianf,iend)
*           ianf=iend
*          enddo
* ausprinten von emat
        else
          write(6,*) 'fehler in den P-faellen fuer dk=1 in bum'
        endif
      endif

* dk=2 fall
      if(idsk.eq.2) then
       cbint=1.0d0
       exint=1.0d0
       call bum21(cbint,exint,irfa,iqrfa,kml,kmj,inopen)
      endif
1000  continue
C      write(6,*) 'die umat fuer cb aus bum'
C      do kk1=1,kml
C       write(6,993) (umat(kk1,kk2),kk2=1,kmj)
C        write(6,*) (umat(kk1,kk2),kk2=1,kmj)
C      enddo
*      write(6,*) 'die umat fuer ex, sodass umat*ex addiert wird'
*      do kk1=1,kml
*        write(6,993) (umat1(kk1,kk2),kk2=1,kmj)
*        write(6,*) (umat1(kk1,kk2),kk2=1,kmj)
*      enddo

993   format(10(F8.4))
      return
      end
