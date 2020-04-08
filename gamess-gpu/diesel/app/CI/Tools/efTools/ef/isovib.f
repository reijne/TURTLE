      program isovib
      
c     ghost
c     this little routine takes info from
c     file hessc and hesscinfo and recalculates
c     the hessian eigenvalues with atomic masses
c     as desired.
c     has to be on file <altmass> !!
c                       ^^^^^^^^^
c     Syntax:
c
c     l1: empty or comment!!! but NEEDS TO BE THERE!
c     l2: atomic no.   mass
c     _______________________________________________



      implicit double precision (a-h,o-z)
      include 'param'
      common /mass/ ams(107),dummy(107)
      common /elemts/ elemnt(107)
      dimension  amass(107)
      dimension  anewmass(107)
      integer*4 numat,nvarn,intflag
      character*80 dummyline
      character*80 title
      character*2 elemnt
      dimension labels(numatm)
      DIMENSION HESS(maxpar,maxpar)
      DIMENSION EIGVAL(maxpar)
      DIMENSION UC(maxpar,maxpar),HESSC(maxpar*(maxpar+1)/2)
      dimension aInputmass(maxpar)
      logical existance
      logical vibdata

c     . say hello!
      print*  
      print*, '#########################################'
      print*, ' '
      print*, ' ISOVIB '
      print*, ' '
      print*, ' Designed for your convenience. '
      print*, ' ghost@silly.thch.uni-bonn.de '
      print*, '                           3.98'
      print*, '  '
      print*, '  '
      print*, '  '
      print*, '  '
      print*, ' Use this little program in combination'
      print*, ' with ef.x to recalculate an IR spectrum '
      print*, ' that has been generated with keyword '
      print*, ' HESSC in the ef.x-input file with  '
      print*, ' different masses. '
      print*, '  '
      print*, ' The altered masses have to be on a file'
      print*, ' named <altmass> !'
      print*, ' '
      print*, ' Its syntax is:'
      print*, ' 1st line:  TITLE or VOID but MUST EXIST!!!'
      print*, ' following: atomic no.    mass'
      print*, '  '
      print*, '#########################################'
      print*

c     . do we have an input?
c     . check for existance of file altmass
      inquire(file='altmass',exist=existance)
      if (.not.existance) then
         print*,'ERROR: File altmass DOES NOT EXIST'
         print*
         print*,'You need to give altered masses in file'
         print*,'altmass, which does not exist.'
         print*
         print*,'Its syntax is:'
         print*
         print*,'one comment line'
         print*,'atomic no.     mass in a.u.'
         print*
         print*
         print*,'sorry - greetings from ghost!        exiting.'
         stop
      else
         continue
      endif


c     . initialize fields
      do i=1,maxpar
         aNewmass(i) = 0.0d0
         aInputmass(i) = 0.0d0
         do j=1,maxpar
            hess(i,j) = 0.0d0
         enddo
      enddo
      


c     . read dump info from file vibdata
      inquire(file='vibdata',exist=vibdata)
      if (.not.vibdata) then
         write(*,*) 'ERROR!'
         write(*,*) 'I cannot find the dumpfile vibdata'
         write(*,*) 'Please provide it by running ef.x'
         write(*,*) 'with keyword HESSC '
         write(*,*) ' '
         write(*,*) 'I apologise for bailing out.    ghost'
         stop
      else
         open(unit=58,form='formatted',file='vibdata')
         read(58,'(A80)') dummyline
         read(58,'(A80)') dummyline
         read(58,'(3I4)') numat,nvarn,intflag
         do i=1,numat
            read(58,'(I3)') labels(i)
         enddo
         do i=1,nvarn
            read(58,'(F24.15)') (hess(i,j),j=1,nvarn)
         enddo
         close(58)
      endif
      

c     . first find out how many lines there are in the input
c       file altmass
      open (unit=47,file='altmass',form='formatted')
      rewind(47)
      read(47,'(A80)') dummyline
      nCounter=0
      do nTotLines = 1,107
         read (47,*,end=100) nDummy,aDummy
         nCounter = nCounter+1
      enddo
 100  continue
     
c     . now read alternate masses
      rewind(47)
      read(47,'(A80)') title
      if (len(dummyline).gt.1) then
         write(*,*) 
         write (*,*) '--- TITLE: --- '
         write(*,*) title
         print*
         print*
         print*
      endif
      write(*,101) nCounter
      write(*,102) numat
      write(*,104) nvarn
      write(*,106) intflag
 101  format ('No of altered isotopic masses:',I3)
 102  format ('You have',I3,' atoms')
 104  format ('Thus, the no. of DOF (nvarn) is',I4)
 106  format ('The intensity flag had been set to ',I1)
      print*
      print*
      write(*,*) 'You wished to change masses for:'
      do 120 nLine=1,numat
         read (47,*,end=130) nAtomnumber,aInputmass(nAtomnumber)
         if (aInputmass(nLine).ne.ams(nAtomnumber)) then
            write(*,*) elemnt(nAtomnumber),' old value:',ams(nAtomnumber
     $           )
            write(*,*) '   new value:',aInputmass(nAtomnumber)
            write(*,*) ' '
            ams(nAtomnumber) = aInputmass(nAtomnumber)
         endif
 120  continue
 130  continue
      close(47)
      print*
      print*
      print*
      print* , 'THE RESULTS ARE:'
      print*



c     this is crucial...: 
      natoms=numat
c     ..............................................................
c     . now copy the rest of the hessc.f routine
c       quick and dirty...
c     ..............................................................
c     . provide for atomic masses
      j=0
      do 291 i=1,natoms
         if(labels(i).lt.99)then
            j=j+1
            atmass=ams(labels(i))
            do 292 k=3*j-2,3*j
               amass(k)=atmass
 292        continue
         endif
 291  continue



c     . NO DUMP HERE !!!


      
C     . SYMMETRIZE HESSIAN and divide the elements of hessian
c       with the square root of the product of the according
c       masses
      DO 230 I=1,nvarn
         DO 220 J=1,I
            HESS(I,J)=(HESS(I,J)+HESS(J,I))*0.50D0
	    hess(i,j)=hess(i,j)/sqrt(amass(i)*amass(j))
            HESS(J,I)=HESS(I,J)
 220     CONTINUE
 230  CONTINUE
      IJ=0
      DO 290 I=1,nvarn
         DO 280 J=1,I
            IJ=IJ+1
            HESSC(IJ)=HESS(J,I)
 280     continue
 290  CONTINUE
      
      CALL HQRII(HESSC,nvarn,nvarn,EIGVAL,UC)
      
c     . now the tough part: UNITS
c     . kcal to J
      f1=4184.0d0
c     . Avogadro
      f2=6.0221367d23
c     . \AA to m
      f3=1.0d-10
      f=f1/(f3*f3*f2*100.0d0)
      negative=0
      do 294 i=1,nvarn
         if(eigval(1).lt.0) negative=1
         if(eigval(i).lt.0) eigval(i)=-eigval(i)
         eigval(i)=f*eigval(i)
         eigval(i)=1302.832d0*sqrt(eigval(i))
 294  continue
      
      
c     TAKE OUT THE INTENSITY SECTION FOR THE TIME BEING
cc     . memorize all components until we're finished distorting
cc     . then transform \mu(X) to d(mu)/d(Q)  (Q=Modes)
c      
c      if (intflag.eq.1) then
c         call transformmu(dipxl,dipyl,dipzl,dipxr,dipyr,dipzr,uc,numat
c     $        ,nvarn,ddipQx,ddipQy,ddipQz,cstep,amass)
c         call calcint(ddipQx,ddipQy,ddipQz,numat,cstep,aIntensity)
c      endif
c      


c     . NEITHER NEED WE DUMP ON HESSIANC HERE.
c
c      open(unit=4,file='hessianc')
c      rewind (4)
c      write(4,*) 'hessian dim=',nvarn
c      if(negative.eq.1) then
c         write(4,*) 'lowest eigenvalue (NEGATIVE) =',eigval(1)
c      else
c         write(4,*) 'lowest eigenvalue (POSITIVE) =', eigval(1)
c      endif
c      do i=1,nvarn
c         write(4,'(6F12.3)')(hess(i,j),j=1,i)
c      enddo
c      close (4)


      call hessout(hess,1,numat)
      
      
c     say goodbye!
      write(*,999)
 999  format(///,'It has been a pleasure working for you -  bye !',/)

      
      stop
      end
