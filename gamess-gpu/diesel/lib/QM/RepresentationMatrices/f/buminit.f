	subroutine buminit(multiplicity, matdim, maxop, pUMat, pUMat1,
     * pEMat
c, php5dar
     * , verbose)
	
	
	
	include "sga.inc"
	include "global"
	
	include "pointers.inc"
	
	
	dimension ibinom(0:nopmax+1,0:nopmax+1)
	integer ibinom, matdim, multiplicity
	integer verbose
	
	
	
	
C	write(6,*) 'HALLO'
	
	pUMat = loc(umat(1,1))
	pUMat1 = loc(umat1(1,1))
	pEMat = loc(emat(1))
c	php5dar = loc(hp5dar(1,1,1))
	umat(1,1)= 12765
	umat(2,1)= 9876
	umat(1,2)= 5555
	matdim = mxeigs
	maxop = nopmax*(nopmax-1)/2
	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C	Code zur Initialisierung der SGA-Table aus refmat/refdg.f
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
	
C	iouts = 44
C	rewind(iouts)
Ccvp  lesen der sga-table von ft44
C      mult=multiplicity
C*      write(6,*) ' multiplizitaet: ',mult
C      do 7010 i=1,mult
C        read(iouts) imult,maxsk,ieigs,iausd,diagtr,ausdtr
C 7010 continue

	if (verbose.gt.0) then
		write(6,*) 'creating table...'
	endif
	call sga(multiplicity, maxsk, ieigs, verbose)
	if (verbose.gt.0) then
		write(6,*) 'ready.'
	endif

C	write(6,*) multiplicity, imult
C      if (multiplicity.ne.imult) stop ' i/o von sga-table falsch '


c testausschrift
C      write(6,*) ' imult=,maxsk=,ieigs= ',imult,maxsk,ieigs
C      do 7020 i=1,ieigs(2)
C        write(6,*) ' funktion nr. ',i
C        write(6,*) ' iausd=,diagtr=,ausdtr= '
C     &            ,(iausd(i,j),j=1,3)
C     &            ,(diagtr(i,j),j=1,3)
C     &            ,(ausdtr(i,j),j=1,3)
C 7020 continue
	
	

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                Folgender Block ist aus refmat/chfeld.f
C                   (initialisieren der iqmo-Felder)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC





cvp aufbau von feldern, die konfigurationsvergleich erleichtern
cvp sollen

cvp
cvp aufbau des ibinom-feldes
cvp ibinom enthaelt die binomialkoeffizienten
cvp (i uber 0) = 1
      do 110 i=0,nopmax+1
  110 ibinom(i,0)=1
      do 115 i=1,nopmax+1
      do 115 j=1,nopmax+1
  115 ibinom(i,j)=0
cvp (0 uber 1) = 0
      do 120 i=0,nopmax+1
  120 ibinom(i,1)=i
      do 130 i=2,nopmax+1
      do 130 j=2,i
  130 ibinom(i,j)=ibinom(i-1,j)+ibinom(i-1,j-1)
  
  
cvp
cvp belegung der zuordnung q-fall <--> mo-nummern fuer dk=0,p=1
cvp iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
cvp folgt der q-fall
cvp imoq(q-fall,1): hoeheres mo
cvp imoq(q-fall,2): niedrigeres mo
      do 140 i=2,nopmax
      do 140 j=1,i-1
      iqmo(i,j)=ibinom(i,2)-i+j+1
      imoq(iqmo(i,j),1)=i
  140 imoq(iqmo(i,j),2)=j
cvp
cvp belegung der zuordnung q-fall --> mo-nummern fuer dk=1,p=1
cvp imo3q(q-fall,1): hoehstes mo
cvp imo3q(q-fall,2): mittleres mo
cvp imo3q(q-fall,3): niedrigstes mo
      do 160 i=3,nopmax
      do 160 j=2,i-1
      do 160 k=1,j-1
      ndum=ibinom(i-1,3)+ibinom(j-1,2)+k
      imo3q(ndum,1)=i
      imo3q(ndum,2)=j
  160 imo3q(ndum,3)=k
cvp
cvp belegung der zuordnung q-fall --> mo-nummern fuer dk=2
cvp imo4q(q-fall,1): hoehstes mo
cvp imo4q(q-fall,2): zweites mo
cvp imo4q(q-fall,3): drittes mo
cvp imo4q(q-fall,4): niedrigstes mo
      do 170 i=4,nopmax
      do 170 j=3,i-1
      do 170 k=2,j-1
      do 170 l=1,k-1
      ndum=ibinom(i-1,4)+ibinom(j-1,3)+ibinom(k-1,2)+l
      imo4q(ndum,1)=i
      imo4q(ndum,2)=j
      imo4q(ndum,3)=k
  170 imo4q(ndum,4)=l


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	
	
	return
	end
