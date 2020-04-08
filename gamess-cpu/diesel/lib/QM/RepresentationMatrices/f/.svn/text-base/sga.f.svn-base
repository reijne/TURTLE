	subroutine sga(multiplicity, mxskp, ieigs, verbose)
cvp                                                                             
c   programm zur berechnung einer sga-table:
c   fuer verschiedene multiplizitaeten werden fuer die hoechste
c   moegliche superkategorie die darstellungsmatrizen der elementaren
c   transpositionen bzgl. genealogisch erzeugter spineigenfunktionen 
c   berechnet. die spineigenfunktionen sind in reversed lexical ordering
c   angeordnet:
c     duch w. u. karwowski j., comp. phys. rep. vol. 2 (1985) p.122
c     rettrup s., int. journ. quant. chem. vol. 29 (1986) p.119
c
c  (vp 30.8.1994)
cvp                                                                             
c 
      include "sga.inc"
c  
c  
c mxsk: nummer der hoechsten sk fuer jedes s  
c ieigs: anzahl der  spineigenfunktionen
      integer mxsk(mxmult),ieigs(mxiswh)
	integer verbose
c  
c itape: input-file
c iout:  standard-output-file
c iouts: output fuer sga-table  
      itape=5
      iout=6
      iouts=44
c einlesen der maximalen spinquantenzahl als 2s+1 sowie der zu jedem
c  s gehoerenden maximalen sk  
C      read(itape,*) mxspin
		mxspin = 7
		mxsk(1) = 7
		mxsk(2) = 10
		mxsk(3) = 6
		mxsk(4) = 6
		mxsk(5) = 5
		mxsk(6) = 5
		mxsk(7) = 5
C		mxsk(1) = 5
C		mxsk(2) = 8
C		mxsk(3) = 4
C		mxsk(4) = 4
C		mxsk(5) = 3
C		mxsk(6) = 3
C		mxsk(7) = 3
      if (mxspin.gt.mxmult) mxspin=mxmult
	if (verbose.gt.0) then
      write(iout,*) 
      write(iout,*) ' max. multiplicity (2s+1): '
     &               ,mxspin
	endif
C      read(itape,*) (mxsk(i),i=1,mxspin)
      do 10 i=1,mxspin
        if (mxsk(i).lt.0) then
          write(iout,*) ' negative superkategorie nicht moeglich '
          write(iout,*) ' superkategorie auf eins gesetzt '
          mxsk(i)=1
        endif
        if (mxsk(i).gt.mxiswh) then
          mxsk(i)=mxiswh
        endif
   10 continue
c     write(iout,*) ' max. superkategorien fuer multiplizitaeten '
c     write(iout,*) (mxsk(i),i=1,mxspin)
c  
c schleife ueber alle spinquantenzahlen


C      do 1000 imult=1,mxspin
C      do 1000 imult=1,multiplicity

	imult = multiplicity



c  inopen: anzahl der offenen schalen, minimale anzahl ist 2s
        inopen=imult-1
c  initialisierung von ieigs
        do 1005 j=2,mxiswh
          ieigs(j)=0
 1005   continue
c  ieigs: anzahl der  spineigenfunktionen, in der untersten sk 1
        ieigs(1)=1
c  nbeta: anzahl der s- operatoren in der spinfkt.
        nbeta=0
c    rekusive berechnung der anzahl der spineigenfunktionen
        maxisk=1
        ixx=1
        do 1010 j=2,mxsk(imult)
          ixx=ixx*4*(inopen+2)*(inopen+1)
     &     / ( (inopen+3)*(inopen+3)-imult*imult )
c         write(6,*) ' imult=,j=,ixx= ',imult,j,ixx
c      inkrementierung nur falls maximale anzahl der spineigenfkt.
c       noch nicht ueberschritten
          if (ixx.le.mxeigs.and.(inopen+2).le.nopmax) then
            ieigs(j)=ixx
            inopen=inopen+2
            maxisk=j
            nbeta=nbeta+1
          endif
 1010   continue
        mxsk(imult)=maxisk
c    neig max. dimension der darstellung fuer geg. spin
c    nur fuer diese dimension wird die sga-table aufgebaut
        neig=ieigs(mxsk(imult))
c  
c  ausschrift der dimensionen der darstellung
	if (verbose.gt.0) then
        write(iout,*)
        write(iout,*) ' actual multiplicity (2s+1): ',imult
        write(iout,*) ' max. open shells: ',inopen
        write(iout,*) ' max. dimension of representation: ',neig
	endif
c  
c  falls nur eine sk bei dieser multiplizitaet dann sprung zum write  
        if (mxsk(imult).eq.1) goto 9000
c  
c    erzeugung des verzweigungsdiagramms der spineigenfunktionen
        call branch(neig,nbeta,inopen)
c  
c    aufbau des feldes igew (gewichte der vertizes im verzweigungsdiagramm)
        call weight(imult,nbeta,inopen)
c  ausschrift des verzweichungsdiagramms
C        write(iout,*) ' gewichte im verzweigungsdiagramm '
C        do 1020 j=imult+nbeta,1,-1
C          write(iout,7010) (igew(i,j),i=0,inopen)
C 1020   continue
 7010   format(20i4)
c  
c  berechnung der matrixelemente der elementaren transpositionen  
        call eltran(imult,neig,inopen)
c
c  ausschrift der informationen  
        do 1030 j=1,neig
          do 1040 k=1,inopen-1
            if (iausd(j,k).ne.0) then
C              write(iout,7020) j,k,k+1,iausd(j,k)
            endif
 1040     continue
 1030   continue
 7020   format(' spineigenfkt. ',i4,' verknuepft ueber (',i2,',',i2,
     &         ') mit spineigenfkt. ',i4)
c  
 9000   continue
c  
c  ausschrift der sga-table auf file iouts 
C       write(6,*) ' imult=,mxsk=,ieigs= ',imult,mxsk(imult),ieigs
C        write(iouts) imult,mxsk(imult),ieigs,iausd,diagtr,ausdtr
        mxskp = mxsk(imult)
c  
C 1000 continue
c  
	return
      end                                                                       
