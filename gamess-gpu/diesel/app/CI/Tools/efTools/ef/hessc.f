      SUBROUTINE HESSIANC(coord,numat)
C     this routine is the master control file for numerical
C     derivation of frequencies and intensities
C     ghost@silly.thch.uni-bonn.de
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param'
      integer dir,double
      COMMON /KEYWRD/KEYWRD
      common /mass/ ams1(107),ams2(107)
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     .     NA(NUMATM),NB(NUMATM),NC(NUMATM)
      CHARACTER*241 KEYWRD
      COMMON /ELEMTS/ ELEMNT(107)
      character*2 elemnt
c     ici
c     character*2 chardummy
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR), IDUMY, XDUM(MAXPAR)
      COMMON /SIGMA2/ GNEXT1(MAXPAR), GRAD(MAXPAR)
      DIMENSION HESS(maxpar,maxpar)
      DIMENSION EIGVAL(maxpar),labelm(maxpar)
      DIMENSION UC(maxpar,maxpar),HESSC(maxpar*(maxpar+1)/2)
      dimension coord(3,numatm),rfcrd(3*numatm)
      dimension amass(numatm)
      dimension dipxl(maxpar)
      dimension dipyl(maxpar)
      dimension dipzl(maxpar)
      dimension dipxr(maxpar)
      dimension dipyr(maxpar)
      dimension dipzr(maxpar)
      dimension ddipQx(maxpar)
      dimension ddipQy(maxpar)
      dimension ddipQz(maxpar)
      dimension aIntensity(3*numatm)
      logical hesscinfo
      logical vibdata
      logical isoav   
c     cstep will be given in AA, rfcrd in a.u.
      
      isoav=.false.

      F=1.0d0/0.529177267260d0
      write(*,*) '************************************************'
      WRITE(*,*) '  NUMER.HESSIAN DETERMINATION IN CART.COORDS.'
      write(*,*) '                              '
      write(*,*) '                               M. Gleichmann'
      write(*,*) '      Add-Ons and much more -- by cm & ghost    '
      write(*,*) ' '
      write(*,*) '************************************************'
     $     
      write(6,*) 'NUMBER OF ATOMS=',numat
      
      nvarn=3*numat
      
      cstep=0.010d0
      I=INDEX(KEYWRD,' CSTEP=')
      IF(I.NE.0) cstep=READA(KEYWRD,I)
c     cstep is in AA
      
      
      double=0
      i=index(keywrd,' DOUBLE=')
      if(i.ne.0) double=READA(KEYWRD,I)
      
      if(double.ne.0.and.double.ne.1)then
         write(*,*) 
         write(*,*)'*******WARNING********'
         write(*,*)'RECONSIDER THE VALUE OF DOUBLE, BUMMER'
         write(*,*)'cart. hessian cannot be calculated'
         write(*,*)'program stops'
         write(*,*)'*    END OF WARNING      * '
         write(*,*) 
         goto 999
      endif
      
      i=index(keywrd,'ISOAV')
      if (i.ne.0) isoav=.true.
      
      intflag=0
      i=index(keywrd,'INT')
      if (double.eq.0.and.i.ne.0) then
         write(*,*)' '
         write(*,*) 'WARNING: I am ignoring your value of DOUBLE ',
     $        'and setting it to 1.'
      endif
      if(i.ne.0) then 
         intflag = 1
         double=1
         write(*,*) ' '
         write(*,*) ' '
         write(*,*) '________________________________'
         write(*,*) '  I WILL CALCULATE INTENSITIES '
         write(*,*) ' '
         write(*,*) '  Ref.:Komornicki and Jaffe,'
         write(*,*) '      JPC 71 (1979) 2150.'
         write(*,*) ' '
         write(*,*) '                 ghost & cm 2.98 '
         write(*,*) '________________________________'
         write(*,*) ' '
      endif
      
      dir=0
      i=index(keywrd,' DIR=')
      if(i.ne.0) dir=reada(keywrd,i)
      if(dir.gt.(3*natoms)) then
         write(*,*) 
         write(*,*)'*******WARNING********'
         write(*,*)'RECONSIDER THE VALUE OF DIR, BUDDY'
         write(*,*)'cart. hessian cannot be calculated'
         write(*,*)'program stops'
         write(*,*)'*    END OF WARNING      * '
         write(*,*) 
         goto 999
      endif
      
      write(*,*) 'KEYWORD DOUBLE =',double
      write(*,*) 
      write(*,*) 'STEPSIZE CSTEP IN ANGSTROM= ',CSTEP
      write(*,*) 'STEPSIZE CSTEP IN A.U.= ',CSTEP/0.52917726d0
      write(*,*) 
      
cmmg  enough initialization data
      
      write(6,*) 'THE REFERENCE GEOMETRY FOLLOWS IN Angstrom'
      write(*,*) 
      DO 241 I=1,natoms
c     READ(11,'(8X,3F14.6,5X,A2)') (REFCRD(J,I),J=1,3),el(I)
         write(*,'(8X,3F14.6,5X,A2)') (COORD(J,I),J=1,3),
     .        elemnt(labels(I))
 241  CONTINUE
      write(*,*) 
      write(6,*) 'END OF REFERENCE GEOMETRY '
      write(*,*) 
      kl=0
c     cici
c     open(file='cemacoord',unit=71,form='formatted')
c     read(71,*) chardummy
      do 215 i=1,natoms
c     read(71,*) (coord(j,i),j=1,3)
c     write(*,*) (coord(j,i),j=1,3)
	 do 214 j=1,3
	    kl=kl+1
c     ici
C     rfcrd(kl)=coord(j,i)
	    rfcrd(kl)=F*coord(j,i)
	    labelm(kl)=labels(i)
 214     continue
 215  continue
c     close(71)
      
      
      
c     ------------------ restart case down from here.
c     in the case of a restart, file <hesscinfo> will be read
c     and the remaining elements for file <hessianc> will
c     be calculated
      
      if(dir.ne.0) then
         inquire(file='hesscinfo',exist=hesscinfo)
         if(.not.hesscinfo) then
            write(*,*)'*******WARNING********'
            write(*,*)'FILE <hesscinfo> CANNOT BE READ'
            write(*,*)'program stops'
            write(*,*)'*   END OF WARNING   *'
            goto 999
         endif
         open(unit=12,file='hesscinfo',status='old')
         rewind 12
         read(12,*) cstep
         write(*,*) 'cstep read from file <hesscinof>'
         write(*,*) 'cstep in angstrom=',cstep
         write(*,*) 'cstep in a.u.=',cstep/0.52917726d0
         write(*,*) ' '
         write(*,*)'read gradient from reference geometry'
         write(*,*) ' '
         read(12,*) i
         read(12,'(5f14.8)') (grad(j),j=1,nvarn)
         
C     k=index for hessian;kk=index for number of calculated
C     gradients apart from reference geometry
         k=0
         kk=0
         do 110 i=1,dir
            if(labelm(i).lt.99)then
               k=k+1
               kk=kk+1
               read(12,*) j
               read(12,'(5f14.8)') (gnext1(j),j=1,nvarn)
               write(*,*) ' '
               write(*,*)'read gradient from elongation geometry',i
               do 101 j=1,nvarn
                  hess(k,j)=(gnext1(j)-grad(j))/cstep
 101           continue
               if(double.eq.1) then
                  kk=kk+1
                  read(12,*) j
                  read(12,'(5f14.8)') (gnext1(j),j=1,nvarn)
                  write(*,*
     $                 )'read gradient from second elongation geometry'
     $                 ,i
                  do 102 j=1,nvarn
                     hess(k,j)=0.5*(hess(k,j)+(-gnext1(j)+grad(j))/cstep
     $                    )
 102              continue
               endif
            endif
 110     continue
         write(*,*) ' '
         close(12)
c     still restart routine. ghost.
         do 304 i=dir+1,(3*natoms)
            if(labelm(i).lt.99)then
               k=k+1
               kk=kk+1
               write(*,*)' '
               write(*,*)'******* '
               write(*,*) 'START WITH ELONGATION OF DIRECTION.',i
               write(*,*)'******* '
               write(*,*)' '
               rfcrd(i)=rfcrd(i)+cstep/0.52917726d0
               call compfgc(gnext1,rfcrd)
               call cout(nvarn,gnext1,i,kk)
               rfcrd(i)=rfcrd(i)-cstep/0.52917726d0
               do 303 j=1,nvarn
                  hess(k,j)=(gnext1(j)-grad(j))/cstep
 303           continue
               if(double.eq.1) then
                  kk=kk+1
                  write(*,*)' '
                  write(*,*)
     $                 'START WITH SECOND ELONGATION OF DIRECTION.',i
                  write(*,*)'******* '
                  rfcrd(i)=rfcrd(i)-cstep/0.52917726d0
                  call compfgc(gnext1,rfcrd)
                  call cout(nvarn,gnext1,i,kk)
                  rfcrd(i)=rfcrd(i)+cstep/0.52917726d0
                  do 301 j=1,nvarn
                     hess(k,j)=0.5*(hess(k,j)+(-gnext1(j)+grad(j))/cstep
     $                    )
 301              continue
               endif
               if(k.eq.nvarn)then
                  write(*,*)' '
                  write(*,*)'+++++++ '
                  write(*,*) 'FINISHED WITH ALL ELONGATIONS'
                  write(*,*)'+++++++ '
                  write(*,*)' '
               endif
            endif
 304     continue
      endif
c     end of restart case

c---------------------- NO RESTART FROM HERE ----
      if(dir.eq.0)then
         open (unit=12,file='hesscinfo')
         rewind 12
         write(12,*) cstep
         close(12)
         call compfgc(grad,rfcrd)
         call cout(nvarn,grad,0,0)
c     . distortions of geometry
C     . k=index for hessian;kk=index for number of calculated
C       gradients away from reference geometry
         k=0
         kk=0
         call readdip(DipRefx,DipRefy,DipRefz)
         
         do 210 i=1,3*natoms
            if(labelm(i).lt.99)then
               k=k+1 
               kk=kk+1 
               write(*,*)' '
               write(*,*)'******* '
               write(*,*) 'START WITH ELONGATION OF DIRECTION.',i
               write(*,*)'******* '
               write(*,*)' '
               rfcrd(i)=rfcrd(i)+cstep/0.52917726d0
               call compfgc(gnext1,rfcrd)
               call cout(nvarn,gnext1,i,kk)
C---  ghost first the 'left' distortion:
               if (intflag.eq.1) then
                  call readdip(dipxl(i),dipyl(i),dipzl(i))       
               endif
               rfcrd(i)=rfcrd(i)-cstep/0.52917726d0
               do 200 j=1,nvarn
                  hess(k,j)=(gnext1(j)-grad(j))/cstep
 200           continue
               if(double.eq.1) then
                  kk=kk+1 
                  write(*,*)' '
                  write(*,*)
     $                 'START WITH SECOND ELONGATION OF DIRECTION.',i
                  write(*,*)'******* '
                  write(*,*)' '
                  rfcrd(i)=rfcrd(i)-cstep/0.52917726d0
                  call compfgc(gnext1,rfcrd)
                  call cout(nvarn,gnext1,i,kk)
C---  ghost now the 'right' distortion:
                  if (intflag.eq.1) then
                     call readdip(dipxr(i),dipyr(i),dipzr(i))
                  endif
                  rfcrd(i)=rfcrd(i)+cstep/0.52917726d0
                  do 201 j=1,nvarn
                     hess(k,j)=0.5*(hess(k,j)+(-gnext1(j)+grad(j))/cstep
     $                    )
 201              continue
               endif
               if(k.eq.nvarn)then
                  write(*,*)' '
                  write(*,*)'+++++++ '
                  write(*,*) 'FINISHED WITH ALL ELONGATIONS'
                  write(*,*)'+++++++ '
                  write(*,*)' '
               endif
            endif
 210     continue
         
      endif
c     end of "no restart" case
      
      
c     . write outputheader
      print*
      print*
      print*, '      - VIBRATIONAL SPECTRUM OUPUT -'
      print*
      print*



c     . provide for atomic masses
      j=0
      do 291 i=1,natoms
         if(labels(i).lt.99)then
            j=j+1
            if(isoav) then
               atmass=ams2(labels(i))
            else
               atmass=ams1(labels(i))
            endif
            do 292 k=3*j-2,3*j
               amass(k)=atmass
 292        continue
            write(*,*) 'atom',i,' mass used ',atmass
         endif
 291  continue


c     ghost
c     . dump everything we need to re-calculate
c       the vib spectrum with isovib on file vibdata.
      inquire(file='vibdata',exist=vibdata)
      if (.not.vibdata) then
         open(unit=58,form='formatted',file='vibdata')
      else
         call system ('mv vibdata vibdata.old')
         print*, 'PLEASE NOTE! I moved file vibdata to vibdata.old'
         open(unit=58,form='formatted',file='vibdata')
      endif
      write(58,*) 'This file contains variables, vecs, and matrices'
      write(58,*) 'natoms,nvarn,intflag,labels,hess'
      write(58,'(3I4)') natoms,nvarn,intflag
      write(58,'(I3)') (labels(i),i=1,natoms)
      do i=1,nvarn
         write(58,'(F24.15)') (hess(i,j),j=1,nvarn)
      enddo
      close(58)
      
      
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
      
c     . memorize all components until we're finished distorting
c     . then transform \mu(X) to d(mu)/d(Q)  (Q=Modes)
      
      if (intflag.eq.1) then
         call transformmu(dipxl,dipyl,dipzl,dipxr,dipyr,dipzr,uc,numat
     $        ,nvarn,ddipQx,ddipQy,ddipQz,cstep,amass)
         call calcint(ddipQx,ddipQy,ddipQz,numat,cstep,aIntensity)
      endif
      
      call hessout(hess,1,numat)
      
      
 999  continue
      return
      end
