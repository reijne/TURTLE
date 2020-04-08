      implicit real*8(a-h,o-z)                                                  
C	implicit none

	include "sgapars.inc"
	

c
c  common-block fuer aufbau der sga-table
c   diagtr: diagonalelemente der elementaren transpositionen (j-1,j):
c           werden mit j-1 angesteuert, ferner nummer der spineigenfkt.
c           als index
c   ausdtr: ausserdiagonalelemente der elementaren transpositionen (j-1,j):
c           werden mit j-1 angesteuert, ferner nummer der spineigenfkt.
c           als index
c   iausd: iausd(i,j-1)=l --> die transposition (j-1,j) verknuepft die
c          spineigenfunktionen i und l miteinander (<i|(j-1,j)|j> =|=0)
c   idiagr enthaelt -1 fuer s- und +1 fuer s+ im verzweigungsdiagramm
c    fuer die neig spineigenfunktionen des graphen
c   igraph beschreibt einzelnen weg im graphen:
c    positionen der s- werden gespeichert
c   ipath: als +1 und -1 expandierte form von igraph 
c   igew: enthaelt die gewichte der vertizes im verzweigungsdiagramm
c         igew nicht dicht, d.h. nullen sind enthalten
      common/sgatab/ diagtr,ausdtr,iausd,idiagr,igraph,ipath
     &            ,igew
c   mehr als mxiswh s- operatoren sind nicht moeglich
      integer idiagr(mxeigs,nopmax),igraph(mxiswh+1)
     &         ,ipath(nopmax)
     &         ,igew(0:nopmax,mxsdg)
     &         ,iausd(mxeigs,nopmax-1)
      real*8 diagtr(mxeigs,nopmax-1),ausdtr(mxeigs,nopmax-1)
c
c
