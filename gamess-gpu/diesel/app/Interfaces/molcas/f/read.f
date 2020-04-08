	program read


	implicit none
	real*8	Work
	dimension Work(100000)
	character*8 Label
	integer*4 nBas, nInts, iRc, iOpt, iComp, iSymLbl, xSymLbl
	integer*4 nsym, iph0, ij, isym, jsym, ijsym, ibas, mxjbas, jbas
	dimension nBas(8)

	iRc=-1
	iOpt=0 
	Call OpnOne(iRC,iOpt,'OneInt',10) 
	If (iRc.ne.0) CallErrMsg ('Rd1Int',001) 
	
	iRc=-1
	iOpt=0 
	iComp=0 
	iSymLbl=1 
	Label='nSym    ' 
	Call RdOne(iRc,iOpt,Label,iComp,nSym,iSymLbl) 
	If (iRc.ne.0) Call ErrMsg('Rd1Int',002) 
	write(*,*) "nSym=", nSym
	iRc=-1
	iOpt=0 
	iComp=0 
	iSymLbl=1
	Label='nBas    ' 
	Call RdOne(iRc,iOpt,Label,iComp,nBas,iSymLbl) 
	If (iRc.ne.0) Call ErrMsg('Rd1Int',003) 
	write(*,*) "nBas=", nBas
	iRc=-1
	iOpt=1 
	iComp=2 
	iSymLbl=1 
	Label='Mltpl  1'
C	Label='Velocity'
C	Label='MltplS 0'
C	Label='Kinetic   ' 
C	Label='OneHam    ' 
	Call RdOne(iRc,iOpt,Label,iComp,nInts,iSymLbl) 
	If (iRc.ne.0) Call ErrMsg('Rd1Int',004) 
	write(*,*) "nInts=", nInts

	iRc=-1
	iOpt=6 
	iSymLbl=1 
C	Label='Mltpl  1'
C	Label='Kinetic   ' 
C	Label='OneHam    ' 
	Call RdOne(iRc,iOpt,Label,iComp,Work(0),iSymLbl) 
	If (iRc.ne.0) Call ErrMsg('Rd1Int',004) 
	xSymLbl = iSymLbl


	ij = 0
	Do iSym = 1, nSym
		Do jSym = 1, iSym
			ijSym = XOR(iSym-1, jSym-1)
			If ( IAND(2**ijSym,xSymLbl).ne.0 ) then
				Do iBas=1,nBas(jSym)
					mxjBas=nBas(iSym)
					if ( iSym.eq.jSym ) mxjBas = iBas
					do jBas=1,mxjBas
						write (*,*) iSym, jSym, jBas, iBas, Work(ij), ij
						ij = ij + 1
					end do
				end do
			end if
		end do
	end do

	
	iRc=-1
	iOpt=0 
	Call ClsOne(iRC,iOpt) 
	If(iRc.ne.0) Call ErrMsg('Rd1Int',005) 
	
	
	STOP
	END
