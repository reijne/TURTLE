      subroutine calcdarhp5d(isk,kml,inopen,hp5dar)
* module_id   : parkeu.x            author      : b.e.
* generation  : 07.12.95            update      :
* calls to    :                     called by   : i2e2
* description : calculates the representation matrices for diagonalcase
cbe uebernommen von hmp5a.f von V.Pless und leicht modifiziert
cvp
cvp parameter und globale common-bloecke sind im include-file
cbe sga.inc, global und ginter enthalten
      include "sga.inc"
      include "global"
c   hilfsgroesse fuer skalarprodukt  
      real*8 xsum,ausd2
      real*8 hp5dar(*)
      ind3(i,j)=max0(i,j)*(max0(i,j)-1)/2+min0(i,j)

C      write(6,*) 'isk,kml,inopen'
c   offene schalen: (ij|ij)*u(ij) fuer i>j
        if (inopen>nopmax) then
           write (6,*) 'Too many open shells in calcdarhp5: ',inopen
           write (6,*) 'increase nopmax in sgapars.inc'
           stop "FATAL ERROR"
        end if
        if (kml>mxeigs) then
           write (6,*) 'Too many spin eigenfunctions: ',kml
           write (6,*) 'increase mxeigs in sgapars.inc'
           stop  "FATAL ERROR" 
        end if
        klauf=0
        index=inopen*(inopen-1)/2
        do 3100 i=2,inopen
          do 3110 j=1,i-1
           klauf=klauf+1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            kstar=jstar+iconf*i
            lstar=jstar+iconf*j
c     initialisierung von umat
            do 7735 ii1=1,kml
            do 7735 ii2=1,kml
              umat(ii2,ii1)=0.0d0
 7735       continue
            do 7736 ii2=1,kml
              umat(ii2,ii2)=1.0d0
 7736       continue
c
c    berechnung von u(ij)=u(j...i)u(i-1...j)
c     alternativ waere z.b. u(i-1...j)^t*u(i,i-1)*u(i-1...j)
            itranl=j
            itranr=i-2
            itrinc=1
            call cycle(kml,kml,itranl,itranr,itrinc)
            itranl=i-1
            itranr=j
            itrinc=-1
            call cycle(kml,kml,itranl,itranr,itrinc)

C	write(6,*) 'UMat:'
	do 987 iii=1, kml
C		write(6,*) (umat(iii,jjj), jjj=1, kml)
c	do 987 jjj=1, kml
           do 987 jjj=1, iii
              
           hp5dar(klauf+(ind3(iii,jjj)-1)*index)=umat(iii,jjj)
c		hp5dar(klauf,iii,jjj) = umat(iii,jjj)
987	continue

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 3110     continue
 3100   continue




 
 5555   continue
        return
        end
      subroutine calcdarhp5f(isk,kml,inopen,hp5dar)
* module_id   : parkeu.x            author      : b.e.
* generation  : 07.12.95            update      :
* calls to    :                     called by   : i2e2
* description : calculates the representation matrices for diagonalcase
cbe uebernommen von hmp5a.f von V.Pless und leicht modifiziert
cvp
cvp parameter und globale common-bloecke sind im include-file
cbe sga.inc, global und ginter enthalten
      include "sga.inc"
      include "global"
c   hilfsgroesse fuer skalarprodukt  
      real*8 xsum,ausd2
      real*4 hp5dar(*)
      ind3(i,j)=max0(i,j)*(max0(i,j)-1)/2+min0(i,j)

C      write(6,*) 'isk,kml,inopen'
c   offene schalen: (ij|ij)*u(ij) fuer i>j
        if (inopen>nopmax) then
           write (6,*) 'Too many open shells in calcdarhp5: ',inopen
           write (6,*) 'increase nopmax in sgapars.inc'
           stop "FATAL ERROR"
        end if
        if (kml>mxeigs) then
           write (6,*) 'Too many spin eigenfunctions: ',kml
           write (6,*) 'increase mxeigs in sgapars.inc'
           stop  "FATAL ERROR" 
        end if
        klauf=0
        index=inopen*(inopen-1)/2
        do 3100 i=2,inopen
          do 3110 j=1,i-1
           klauf=klauf+1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            kstar=jstar+iconf*i
            lstar=jstar+iconf*j
c     initialisierung von umat
            do 7735 ii1=1,kml
            do 7735 ii2=1,kml
              umat(ii2,ii1)=0.0d0
 7735       continue
            do 7736 ii2=1,kml
              umat(ii2,ii2)=1.0d0
 7736       continue
c
c    berechnung von u(ij)=u(j...i)u(i-1...j)
c     alternativ waere z.b. u(i-1...j)^t*u(i,i-1)*u(i-1...j)
            itranl=j
            itranr=i-2
            itrinc=1
            call cycle(kml,kml,itranl,itranr,itrinc)
            itranl=i-1
            itranr=j
            itrinc=-1
            call cycle(kml,kml,itranl,itranr,itrinc)

C	write(6,*) 'UMat:'
	do 987 iii=1, kml
C		write(6,*) (umat(iii,jjj), jjj=1, kml)
c	do 987 jjj=1, kml
           do 987 jjj=1, iii
              
           hp5dar(klauf+(ind3(iii,jjj)-1)*index)=umat(iii,jjj)
c		hp5dar(klauf,iii,jjj) = umat(iii,jjj)
987	continue

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 3110     continue
 3100   continue




 
 5555   continue
        return
        end
