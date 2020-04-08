	subroutine	fmain(refcsfs, csfs)
	integer	refcsfs, csfs
	
	parameter (n=4, dim=1000)

	real*8	x, y
	dimension	x(dim*n), y(dim*n)

	real*8	h
	dimension	h(dim*dim)
	
	integer ind



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C	y = H*x
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C	Aufbau der Vektoren

	do 100 i=1,csfs
		x(i+(1-1)*csfs) = i
		x(i+(2-1)*csfs) = 2*i
		x(i+(3-1)*csfs) = 3*i
		x(i+(4-1)*csfs) = 4*i
100	continue
	write(*,*) "HALLO"

C	y = H*x
	call hmult(n, x, y)

C	Ausgabe der Vektoren
	do 200 i=1,csfs
		write(*,*)  y(i+(1-1)*csfs),y(i+(2-1)*csfs),
     * 	y(i+(3-1)*csfs),y(i+(4-1)*csfs)
200	continue


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C	Referenzraummatrix, Indexzuordung
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


	call refmat(h)
	
	do 300 i=1,refcsfs
		write(*,*)  (h(j + (i-1)*refcsfs), j=1, refcsfs)
300	continue

	do 400 i=1,refcsfs
		call refindex(i, ind)
		write(*,*)  ind
400	continue
	
	return
	end	
	
	
	
	
