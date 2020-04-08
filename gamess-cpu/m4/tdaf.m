c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/tdaf.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   tda    =
c ******************************************************
c ******************************************************
      subroutine tamctl(nmaxa,nmaxb,nmaxc,nmaxd,nmaxe,nmaxf,nmaxg,
     * nmaxh,nmaxi,nmaxj,nmaxcj,nmaxgj,
     * indvec,jndvec,kndvec,lndvec,
     * iv,jv,kv,iconve,ivecb,iveca,
     * v,ga,sta,g,gg,gvec,aux,aa,bb,iwr)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iv,jv,kv,iconve,iveca,ivecb
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *             selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *             xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
c
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      dimension iv(nmaxgj),jv(nmaxgj),kv(nmaxgj),iconve(nmaxgj)
      dimension iveca(8,nmaxi,nmaxcj),ivecb(8,nmaxi,nmaxj)
      dimension v(nmaxe),ga(nmaxg),sta(nmaxg,nmaxc)
      dimension g(nmaxb,nmaxb),gg(*),gvec(nmaxb)
      dimension aux(nmaxc,nmaxb),aa(nmaxa,nmaxa),bb(nmaxa,nmaxa)
c
    1 format (/1x,104('*'))
    2 format (1x,'for symmetry type',i5,' the dimension of a = ',i5,
     *' : the dimension of b = ',i5//)
c
c     set maximum number of orbitals for self-energy matrix
      if (nmaxj.gt.nmaxb) nmaxj=nmaxb
      nmaxbg=nmaxb+nmaxg
      maxbg2=nmaxbg+nmaxbg
      mxbg22=maxbg2+nmaxc
      mwork=(nmaxb*nmaxb/2)+1
c
      conver=27.2116d0
      cons=conver*conver
      timfac=1.66666667d-04
      icdiag=0
      rewind nf12
      rewind nf13
      rewind nf14
      rewind nf15
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      call dimen(nmaxb,nmaxc,nmaxd,nmaxg,nmaxh,iwr)
      write (iwr,1)
      fac1= dsqrt(3.0d0)*0.5d0
c     fac2=1.5d0
      fac=1.0d0/ dsqrt(2.0d0)
      fad= dsqrt(1.5d0)
      rewind nf15
      rewind nf16
      do 40 isym=1,nsym
      do 20 i=1,norb
      iorb=indorb(i)
      if (idorb(iorb).eq.isym) go to 30
   20 continue
      go to 40
   30 write(iwr,3000)isym
3000  format(/20x,16('*')/
     *20x,'symmetry type',i3/
     *20x,16('*')/)
      write (iwr,2) isym,ndimat(isym),ndimbt(isym)
      call diaga(nmaxa,nmaxc,nmaxe,nmaxg,nmaxgj,nmaxj,iv,jv,kv,iconve,
     * ivecb,v,ga,sta,aa,aa,
     * indvec,jndvec,kndvec,lndvec,iwr)
      call diagb(nmaxa,nmaxb,nmaxc,nmaxe,nmaxg,nmaxgj,nmaxj,iv,jv,kv
     *,iconve,ivecb,v,ga,sta,aa,
     * indvec,jndvec,kndvec,lndvec,iwr)
   40 continue
      rewind nf13
      rewind nf15
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      do 100 isym=1,nsym
      do 50 i=1,norb
      iorb=indorb(i)
      if (idorb(iorb).eq.isym) go to 60
   50 continue
      go to 100
   60 write(iwr,3000)isym
      call gammaa(nmaxa,nmaxb,nmaxbg,nmaxc,nmaxf,nmaxg,nmaxh,nmaxi,
     * iv,jv,kv,iconve,iveca,  v,g,gg,gvec,aux,aa,
     * indvec,jndvec,kndvec,lndvec,iwr)
      call gammab(nmaxa,nmaxb,nmaxbg,nmaxc,nmaxf,nmaxg,nmaxi,
     * iv,jv,kv,iconve,ivecb,  v,g,gg,gvec,aux,aa,
     * indvec,jndvec,kndvec,lndvec,iwr)
  100 continue
      rewind nf12
      rewind nf13
      rewind nf15
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      if (tsele.eq.ssele) call compol(nmaxbg,nmaxc,maxbg2,
     * iv,jv,gg(1),gg(maxbg2+1))
      rewind nf13
      rewind nf15
      rewind nf17
      if (tsele.ne.ssele) go to 1000
c     dynamic core for polser
      m2=nmaxc+nmaxc
      l2=nmaxc*(nmaxc+1)/2
      l3=nmaxc*nmaxc
      i1=1
      i2=i1+l2
      i3=i2+l2
      i4=i3+l2
      i5=i4+m2
      i6=i5+nmaxc
      i7=i6+m2
      i8=i7+m2
      i9=i8+nmaxc
      i10=i9+l3
      i11=i10+maxbg2
      i12=i11+maxbg2*nmaxc
      i13=i12+(maxbg2-1)/lenwrd()+1
      i14=i13+mxbg22
      i15=i14+mxbg22
c     i16=i15+mxbg22
c
      call polser(nmaxa,nmaxb,nmaxbg,maxbg2,mxbg22,nmaxc,m2,l2,
     * nmaxf,nmaxg,nmaxi,
     * v,g,gg(mwork),iv,jv,aa,bb,
     * gg(i1),gg(i2),gg(i3),gg(i4),gg(i5),gg(i6),gg(i7),
     * gg(i8),gg(i9),gg(i10),gg(i11),gg(i12),gg(i13),gg(i14),gg(i15)
_IF(ibm,vax)
     *,iwr )
_ELSE
     *,indvec,jndvec,kndvec,lndvec,iwr )
_ENDIF
1000  rewind nf12
      rewind nf13
      rewind nf15
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      if (tsele.eq.ssele) go to 2000
c
      if(ztdaop(1).eq.'diagonal')then
c
c  diagonal approximation
c
      call diadig(nmaxa,nmaxb,nmaxc,nmaxh,g,gvec,aux,iwr)
c
      rewind nf12
      rewind nf13
      rewind nf15
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      endif
c
      call comdig(nmaxa,nmaxb,nmaxc,nmaxf,nmaxi,iveca,ivecb,
     * v,g,gg(mwork),gvec,aux,aa,bb
_IF(ibm,vax)
     *,iwr )
_ELSE
     * ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
2000  rewind nf13
      rewind nf15
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      if (selctd.ne.sseld) go to 200
      call quart(nmaxb,nmaxf,v,g,gvec,
     * indvec,jndvec,kndvec,lndvec,iwr)
  200 continue
      call ltimet(iwr)
c
c
c     use of data sets
c     a) diagonalization method
c        1) file 13 written in datgab, read in comdig, adiagr
c        2) file 14 written in input, read in correl
c        3) file 15 written in diaga, read in gammaa, datgaa
c        4) file 16 a) written in diagb, read in gammab, datgab
c                   b) used in comdig
c        5) file 17 written in datgaa, read in comdig, adiagr
c        6) file 18 used in comdig
c        7) file 19 used in comdig
c
c     b) pole search method
c        1) file 12 intermediate results for each symmetry (polser)
c        2) file 13 a) written in patgab, read in compol
c                   b) results of pole search for all symmetries
c                      used in polser
c        3) file 14 written in input, read in pcorre
c        4) file 15 a) written in diagax, read in gammaa, patgaa
c                   b) written in compol, read in polser, psadia
c        5) file 16 a) written in diagbx, read in gammab, patgab
c                   b) configuration identification for all symmetries
c                      used in polser
c        6) file 17 a) written in patgaa, read in compol
c                   b) results of pole search for one symmetry
c                       written in psa, read in polser
c        7) file 18 used in polser
c        8) file 19 used in polser
c
c
c     use of common blocks
c     tamda  : v       : integrals
c              ga      : diagonal elements of a- and b-block matrices
c                        used in selection
c              sta     : corresponding amplitudes
c              g       : a- and b-block matrices, dyson equation matrix
c                        and eigenvectors after diagonalization
c              gvec    : eigenvalues after diagonalization of g
c              aux     : eigenvectors of a given symmetry, strengths
c     kjinv  :           addressing vectors for integrals
c     work   : sigdia  : work area
c              epsil   : hartree-fock orbital energies
c              aa,bb   : constant diagrams, work area
c     york   : bmvec   : identification vector of configuration types
c                        used for output
c              fmul    : multiplicative constants for a-diagrams in case
c                        inversion fails to aid convergence in the
c                        iteration
c              range   : upper and lower range of pole search . 0 is in
c                        general lower limit for ionization calculations
c              alimit  : diagonal affinity poles below alimit are
c                        combined into centroids. affinity poles are
c                        taken to be negative.
c              blimit  : diagonal ionization poles above blimit are
c                        combined into centroids. ionization poles are
c                        taken to be positive.
c     data   : nocc    : number of occupied orbitals
c              nvirt   : number of virtual orbitals
c              nstart  = nocc+1
c              ntot    : total number of orbitals
c              nvirtx  : number of virtual orbitals which are taken into
c                        account in nondiagonally approximated poles
c              ntotx   : total number of orbitals which are taken into
c                        account in nondiagonally approximated poles
c              ninta   : number of integrals in the first part of the
c                        calculations ( contains the integrals with
c                        orbitals from ntotx+1 to ntot as well)
c              nintb   : number of integrals in the second part of the
c                        calculation
c              nsym    : number of symmetry species
c              imul    : group multiplication table
c              idorb   : symmetry identification of orbitals
c              norb    : number of orbitals for which tda calculation
c                        is to be performed
c              indorb  : their indices
c              iv,jv,kv: orbital indices of 2h1p and 2p1h configurations
c              iconve  : configuration identification vector
c              xtyp    : type of configuration, singlet or triplet
c     datb   : conver  : conversion factor from atomic units to
c                        electron volts
c              cons    = conver**2
c              fac     : constant
c              fac1    : constant
c              fad     : constant
c              word    : input : adiagr or no adiag
c              amvec   : vector of important components of eigenvectors
c              mvec    : their indices in the eigenvector
c              ndima   : dimension of a-block of given symmetry ( first
c                        part, singlet coupled configurations )
c              ndimb   : dimension of b-block of given symmetry ( first
c                        part, singlet coupled configuartions )
c              ndimat  : total dimension of a-block of given symmetry
c              ndimbt  : total dimension of b-block of given symmetry
c              ndimct  : number of affinity poles taken into account
c                        exactly in the diagonalization
c              nmaxa to nmaxj : maximum dimension as explained in tammda
c              isym    : current symmetry species
c              ndim    : dimension of matrices to be diagonalized
c              npola   : number of poles taken into account for the a-
c                        block in the case of pole search
c              npolb   : number of poles taken into account for the b-
c                        block in the case of pole search
c              npolt   : number of solutions used in polser
c              ixxx    : identification which solution is the first
c                        ionization pole of the green's function.
c                        this is determined in comdig/polser and used
c                        there and in correl/pcorre.
c              itemp   : identifies important solutions ( with large
c                        pole strength) used for output in comdig
c                        for microfiche only.
c     selcts : selcta  : input : selection or no select
c              selctb  : input from unit 20 : closed, correl or open
c              xxsel   : identification word : open
c              yysel   : identification word : correl
c              zzsel   : identification word : closed
c              sselx   : identification word : selectio
c              ssely   : identification word : degenera
c              selctc  : input : degeneracy or no degen
c              selctd  : input : quart or no quart
c              sseld   : identification word : quart
c              ssele   : identification word : psa
c              tsele   : input : psa or no psa
c     degen  : lsym    : actually used number of symmetry species
c                        (can be equal to nsym without harm).
c              nsymg   : gives the number of subparts of each symmetry
c                        i.e. the number of cases a prticular represen-
c                        tation of the non-degenerate group occurs in
c                        the full group including the degenerate
c                        representations.
c              iuse    : gives the information whether the dyson
c                        equation is to be solved (1) or not (0)
c                        for the given symmetry and subpart
c              isgnum  : gives the number of orbitals belonging
c                        to the given symmetry and subpart and with
c                        indices less or equal to ntot
c              isgind  : gives the orbital indices belonging to the
c                        given symmetry and subpart. here the orbital
c                        indices refer to the actual number in the
c                        entire set of orbitals.
c              jsgjnd  : gives the orbital indices belonging to the
c                        given symmetry and subpart. here the orbital
c                        indices refer to the number in the symmetry
c                        species i only.
c              ivecb   : identifies diagonalized poles after deletion
c                        of weak ionization poles centroids are added
c                        here.
c              iveca   : identifies diagonalized poles after deletion
c                        of weak affinity poles. centroids are added
c                        here.
c              idimb   : number of ionization poles taken into account
c                        exactly in the diagonalization method
c                        plus the number of centroids. determined in
c                        datgab used in comdig.
c              jdimb   : number of ionization poles taken into account
c                        exactly in the diagonalization method. this
c                        is determined in datgab and used in comdig
c                        only to determine the first ionization
c                        solution.
c              idima   : number of affinity poles taken into account
c                        exactly in the diagonalization method
c                        (maximally nmaxc) plus the number of centroids
c                        determined in datgaa and used in comdig.
c              ksgnum  : gives the number of orbitals belonging to the
c                        given symmetry and subpart and with indices
c                        less or equal to ntotx.
c              ndeg    : degeneracy of orbitals ( 1, 2 or 3).
c              indeg(i,j) : gives for each orbital i which orbitals j
c                        are degenerate with it.
c              idega   : parameter for the degeneracy treatment which
c                        is detemrined form the input data. it has
c                        the standard value 0 if there is no degeneracy
c                        or if there is degeneracy and only selected
c                        components are calculated in the dyson equation
c                        it has the value of 1 if degeneracy occurs and
c                        if all subparts of the occurring symmetries are
c                        going to be solved for in the dyson equation.
c     xtime  : timfac  : time conversion factor from 1/100 of a second
c                        to minutes.
c              timei   : timing data
c              icdiag  : counting parameter for diagonalizations
c              idimdi  : dimension of matrices which are diagonalized
c
c
      return
      end
      subroutine intda (nmv,iwr)
c
      implicit REAL  (a-h,o-z)
c
c      character*1 xtyp
      character*8  selct,asym,zsel
      character*8  xnod,asel
c     character*8  seld,sele,selx,sely,xsel,ysel
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
      character*8 adorb
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
     *,pop(maxorb),potn(2),ncolo(6),ncolp(maxorb,4),
     * nacta(maxorb,2),nactz(6),norb2(maxorb,3),norb22(4),
     * tempin(2),atemp,btemp,isymm,nxx,lsymm,insymm(8),iiuse(8,6),
     * jdorb(100),ncont,
     * nmaxa,nmaxb,nmaxc,nmaxd,nmaxe,nmaxf,nmaxg,nmaxh,nmaxi,nmaxj
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(100)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      dimension adorb(100)
      dimension iamul(8,8),ibmul(8,8)
c
    1 format (/1x,104('-')/)
    3 format (//' number of occupied orbitals = ',i3/
     + ' number of virtual orbitals = ',i3)
    5 format (//' orbital energies in atomic units'//)
    6 format (10(2x,f10.6))
    7 format (///' satellite bands are to be calculated for   ',i6,'
     * orbitals'//)
    8 format (' their indeces are'//)
    9 format (20(3x,i3))
   10 format (//1x,'the maximum index of a virtual orbital',
     +' for which non diagonal'/1x,
     + 'a- and b-block matrix elements are going to be',
     +' calculated is',i7/
     +1x,'the maximum number of orbitals for the nondiagonal'/1x,
     +'calculation is thus',i7//)
   11 format (//' orbital energies in ev'//)
   12 format (' number of symmetry species',5x,i3/)
   13 format (' group multiplication table '//)
   14 format (10(2x,i4))
   15 format (/1x,
     +'symmetry and degeneracy information for the orbitals'/)
   16 format (' orbital   ',i4,'   class',i4,'   degeneracy   ',i4)
   17 format (//3x,i8,' integrals are to be read in and processed')
c  18 format (7f10.6)
   22 format (1x,
     +'the number of integrals in the calculation of nondiagonal'/
     +1x,'a- and b-block matrix elements is',i9//)
   23 format (/1x,
     +'the calculation is done with the inclusion of the energy',
     +' independent diagrams'/)
   24 format (/1x,
     +'the calculation is done without the energy',
     +' independent diagrams')
c  26 format (//' the tda calculation will be done with the separate c
c    *alculation of the ground state correlation energy'//)
   27 format (//1x,
     + 'the tda calculation will be done without the separate'/
     +1x,'calculation of the ground state correlation energy'//)
c  28 format (//' the label of the data set on tape is   ',5a8)
   29 format (1x,
     +'the calculation is done with possible configuration selection',
     +' in the a- and b-blocks'/)
   30 format (1x,
     +'configuration selection is only permitted in the a-block'/)
   31 format (1x,
     +'the actual symmetry group of the molecule is assumed not'/
     +1x,'to contain degenerate representations'/)
   32 format (1x,
     +'the actual symmetry group of the molecule does contain',
     *' degenerate representations'/)
   33 format (' actual number of symmetry species used = ',i3)
   34 format (' the number of subparts of symmetry species ',i3,
     *' equals ',i3)
   35 format (' for symmetry species ',i3,' and subpart ',i3,
     *' the dyson equation will be solved')
   36 format (' for symmetry species ',i3,' and subpart ',i3,
     *' the dyson equation will not be solved')
   37 format (1x,
     +'the number of orbitals belonging to symmetry species * ',i3/
     +' and subpart ',i3,' is ',i5)
   38 format (/' the orbital indeces belonging to symmetry species '
     *,i3,' and subpart ',i3,' are  : '//3x,30i4)
   39 format (/1x,
     + 'reduced orbital indeces now referring only to the numbering'/
     + ' within a given symmetry species')
   40 format (/1x,
     +' the reduced orbital indeces belonging to symmetry species ',
     +i3,' and subpart ',i3,' are  : '//3x,30i4)
   41 format (/' quartett states will be calculated')
   42 format (/' quartett states will not be calculated')
   44 format (/1x,
     +'the constant diagrams are only calculated with the',
     +' hartree-fock greens function')
   45 format (/1x,
     +'the constant diagrams are calculated in the second step',
     *' with the hartree-fock greens function'/1x,
     +'and finally with the full greens function')
   46 format (//1x,
     +'the matrix to be inverted in the calculation of the constant'/
     +1x,'diagrams has dimension ',i5/
     +' and exceeds the dimension nmaxb = ',i5/
     +' the number of orbitals has been reduced to ',i5/
     +' in the entire calculation'/)
   47 format (//' for symmetry ',i3,' the number of orbitals ',i3/
     *' exceeds the maximum allowed one ',i3/1x,
     +'the total number of orbitals is reduced to ',i3/
     +' and the number of virtual orbitals is reduced to ',i3/)
   48 format (1x,
     + 'a degeneracy larger than 3 is found for orbital',i3,
     *' degeneracy ',i3)
   49 format (1x,
     + 'the dyson equation is to be solved for more than one'/1x,
     + 'component of a degenerate representation but not for all'/
     +1x,'components'/1x,
     + 'with the calculation of the constant diagrams or the ground'/
     * 1x,'state correlation energy the program cannot do this.'/1x,
     * 'modify the input for iuse to either calculate all components'/
     + 'or only exactly one for each degenerate representation.'/
     *' the error was detected for orbital   ',i3)
c  50 format (14(1x,a4))
   51 format (' orbital ',i4,' class',i4,' degeneracy ',i4,
     *' representation in the degenerate group ',a4)
   52 format (/1x,
     +'the solution of the dyson equation is done with the',
     +' pole search method'/)
   53 format (/' lower limit of range for pole search =  ',f12.2,'ev'/
     *' upper limit of range for pole search =  ',f12.2,'ev'/)
   54 format (/' self-energy affinity poles below ',f12.2,
     + ' ev will be combined into centroids'/
     + ' self energy ionization poles above ',f12.2,
     + ' ev will be combined into centroids'/
     + ' both quantities are used only in the case of the pole search'/)
   55 format (/1x,
     + 'the solution of the dyson equation is done with '/
     + 'the normal diagonalization method'/)
c
      data iamul  /1,2,3,4,5,6,7,8,
     *             2,1,4,3,6,5,8,7,
     *             3,4,1,2,7,8,5,6,
     *             4,3,2,1,8,7,6,5,
     *             5,6,7,8,1,2,3,4,
     *             6,5,8,7,2,1,4,3,
     *             7,8,5,6,3,4,1,2,
     *             8,7,6,5,4,3,2,1/
c
      data ibmul  /1,2,3,4,5,6,7,8,
     *             2,1,5,6,3,4,8,7,
     *             3,5,1,7,2,8,4,6,
     *             4,6,7,1,8,2,3,5,
     *             5,3,2,8,1,7,6,4,
     *             6,4,8,2,7,1,5,3,
     *             7,8,4,3,6,5,1,2,
     *             8,7,6,5,4,3,2,1/
c
      data selct  /'adiagram'/
      data asym   /'a       '/
c     data xsel   /'open    '/
c     data ysel   /'correl  '/
      data zsel   /'closed  '/
c     data selx   /'selectio'/
c     data sely   /'degenera'/
c     data seld   /'quart   '/
c     data sele   /'psa     '/
      data xnod   /'no adiag'/
      nsym=8
      range(1)=tempin(1)
      range(2)=tempin(2)
      alimit=atemp
      blimit=btemp
      asel=' '
      if(isymm.le.1)asel=asym
      nx=nxx
      lsym=lsymm
      do 40010 i=1,nsym
      nsymg(i)=insymm(i)
      do 40010 j=1,nmaxi
40010 iuse(i,j)=iiuse(i,j)
      do 40020 i=1,100
40020 adorb(i)=bmvec(i)
      if(selctc.ne.ssely)go to 40030
      do 40040 i=1,100
40040 idorb(i)=jdorb(i)
40030 idega=0
      do 100 i=1,nsym
      do 100 j=1,nmaxi
      isgnum(i,j)=0
      ksgnum(i,j)=0
      idima(i,j)=0
      idimb(i,j)=0
      jdimb(i,j)=0
      do 100 k=1,nmaxc
      isgind(i,j,k)=0
      jsgjnd(i,j,k)=0
  100 continue
      numend=-2
      if (nx.eq.0) word=xnod
      if (nx.eq.0) numend=-2
      if (nx.eq.1) numend=-1
      if (nx.eq.2) numend=2
      if (nx.eq.3) numend=6
      if (word.ne.selct) numend=-2
c     ntotk=ntot
      do 120 i=1,ntot
 120  sigdia(i)=epsil(i)*conver
      if (asel.eq.asym) go to 200
      do 180 i=1,nsym
      do 180 j=1,nsym
      imul(i,j)=ibmul(i,j)
  180 continue
      go to 240
  200 do 220 i=1,nsym
      do 220 j=1,nsym
      imul(i,j)=iamul(i,j)
  220 continue
  240 continue
  260 continue
      do 340 isym=1,nsym
      ic=0
      do 280 i=1,ntot
      if (idorb(i).ne.isym) go to 280
      ic=ic+1
  280 continue
      if (ic.le.nmaxc) go to 340
      ntot=ntot-1
      nvirt=nvirt-1
      do 300 nn=1,6
      if ( dabs(epsil(ntot)-epsil(ntot+1)).gt.0.1d-04) go to 320
      ntot=ntot-1
      nvirt=nvirt-1
  300 continue
  320 continue
      write (iwr,47) isym,ic,nmaxc,ntot,nvirt
      go to 260
  340 continue
c
c
c     range(i) sets the lower and upper limit for the pole search
c     procedure in ev.
c
c     alimit sets the lower bound on the value of the self-energy
c     affinity poles which are included. poles beyond this value
c     are combined to form centroids. affinity poles are (with the
c     exception of positive electron affinities ) of negative energy.
c
c     blimit sets the upper bound on the value of the self-energy
c     ionization poles which are included. poles beyond this value
c     are combined to form centroids. ionization poles are of
c     positive energy.
c
c     both values are to be given in ev.
c
      if (selctc.eq.ssely) go to 360
      go to 680
  360 continue
c
c     lsym is the number of symmetry species actually occurring
c
c     nsymg(i) gives for each symmetry species i the number of subparts,
c     i.e. the number of cases a particular representation of the
c     non-degenerate subgroup occurs in the full group including the
c     degenerate representations
c
c     the dyson equation will be solved separately for each subpart !
c
c
c     the array iuse(i,j) gives for each symmetry species i the
c     information if the dyson equation is going to be solved for
c     subpart j of species i:
c     if the dyson equation is not to be solved for subpart j,
c     iuse(i,j)=0, otherwise it is not zero.
c
c     the order of the subparts must be the following if the quantities
c     isgnum and isgind are generated via input of the vector adorb(i)
c     in the subroutine insym :
c     1) nondegenerate representations of the degenerate group
c     2) doubly degenerate representations of the degenerate group whose
c        components are both correlated with the same representation
c        of the nondegenerate group
c     3) doubly degenerate representations of the degenerate group whose
c        components are both correlated with different representations
c        of the nondegenerate group
c     4) triply degenerate representations
c
c     otherwise the ordering follows the order of the orbitals as
c     specified in the vector idorb ( and adorb ).
c
c
      do 420 i=1,ntot
      x=epsil(i)
      js=i-5
      if (js.lt.1) js=1
      je=i+5
      if (je.gt.ntot) je=ntot
      ic=0
      do 400 j=js,je
      if ( dabs(epsil(j)-x).gt.0.1d-04) go to 400
      ic=ic+1
      indeg(i,ic)=j
  400 continue
      ndeg(i)=ic
      if (ic.le.3) go to 420
      write (iwr,48) i,ic
      call caserr('invalid orbital degeneracy in tda')
  420 continue
      ic=0
      do 470 i=1,nsym
      do 440 j=1,ntot
      if (idorb(j).eq.i) go to 450
  440 continue
      go to 470
  450 continue
      jend=nsymg(i)
      do 460 j=1,jend
      if (iuse(i,j).eq.1) go to 460
      ic=1
  460 continue
  470 continue
      if (ic.eq.0) idega=1
      if (ic.eq.0) go to 520
      do 510 ix=1,norb
      i=indorb(ix)
      m=ndeg(i)
      if (m.eq.1) go to 510
      if (m.eq.3) go to 490
      ma=indeg(i,1)
      if (ma.eq.i) ma=indeg(i,2)
      do 480 j=1,norb
      if (ma.ne.indorb(j)) go to 480
      write (iwr,49) i
      call caserr('invalid degeneracy specification in tda')
  480 continue
      go to 510
  490 continue
      ma=indeg(i,1)
      mb=indeg(i,2)
      mc=indeg(i,3)
      ic=0
      do 500 j=1,norb
      if (indorb(j).eq.ma) ic=ic+1
      if (indorb(j).eq.mb) ic=ic+1
      if (indorb(j).eq.mc) ic=ic+1
  500 continue
      if (ic.eq.1) go to 510
      write (iwr,49) i
      call caserr('invalid degeneracy specification in tda')
  510 continue
  520 continue
c
c
      call insym(adorb,nmaxi,iwr)
c
c
c     the vector adorb contains the symmetry identification of the
c     orbitals in the degenerate group. this is used to generate the
c     quantities isgnum and isgind. it avoids complicated input of
c     these quantities.
c
c     do 540 i=1,lsym
c     jend=nsymg(i)
c     read (5,2) (isgnum(i,j),j=1,jend)
c 540 continue
c
c     the array isgnum(i,j) gives for each symmetry species i and
c     part j the number of orbitals belonging to this subpart
c
c     do 560 i=1,lsym
c     jend=nsymg(i)
c     do 560 j=1,jend
c     kend=isgnum(i,j)
c     read (5,2) (isgind(i,j,k),k=1,kend)
c 560 continue
c
c     the array isgind(i,j,k) gives for each symmetry species i and
c     subpart j the orbital indeces ( their number is isgnum(i,j))
c     belonging to this subpart. here the orbital indeces refer to
c     their actual number in the entire set of orbitals. later an this
c     will be reduced to refer to the number in the symmetry species
c     i only. this index will be stored in the array jsgjnd(i,j,k).
c
      do 660 i=1,lsym
      jend=nsymg(i)
      ic=0
      do 640 k=1,ntot
      if (idorb(k).ne.i) go to 640
      ic=ic+1
      do 620 j=1,jend
      lend=isgnum(i,j)
      do 580 l=1,lend
      if (k.eq.isgind(i,j,l)) go to 600
  580 continue
      go to 620
  600 jsgjnd(i,j,l)=ic
  620 continue
  640 continue
  660 continue
      go to 760
  680 continue
c
c     preparation of the input for the case no degenerate
c     representations are existing.
c
      lsym=nsym
      do 720 i=1,lsym
      nsymg(i)=1
      iuse(i,1)=1
      ic=0
      do 700 k=1,ntot
      if (idorb(k).ne.i) go to 700
      ic=ic+1
      isgind(i,1,ic)=k
      jsgjnd(i,1,ic)=ic
  700 continue
      isgnum(i,1)=ic
      if (ic.eq.0) iuse(i,1)=0
  720 continue
      do 740 i=1,ntot
      ndeg(i)=1
      indeg(i,1)=i
  740 continue
  760 continue
      if (ntot.le.nmaxa) go to 800
      ntot=nmaxa
c
c     check for degeneracy
c
      do 780 nn=1,6
      if ( dabs(epsil(ntot)-epsil(ntot+1)).gt.0.1d-04) go to 800
      ntot=ntot-1
  780 continue
  800 continue
c
c     check whether the number of orbitals has to be reduced in the case
c     that the constant diagrams are to be included.
c
      if (word.ne.selct) go to 1060
  820 continue
      ic=0
      do 940 isym=1,nsym
      do 840 i=1,norb
      ind=indorb(i)
      if (idorb(ind).eq.isym) go to 860
  840 continue
      go to 940
  860 continue
      jend=nsymg(isym)
      do 920 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 920
      kend=isgnum(isym,jp)
      do 900 i=1,kend
      do 880 j=i,kend
      ic=ic+1
  880 continue
  900 continue
  920 continue
  940 continue
      if (ic.le.nmaxb) go to 1040
      ntot=ntot-1
      do 960 nn=1,10
      if ( dabs(epsil(ntot)-epsil(ntot+1)).gt.0.1d-04) go to 980
      ntot=ntot-1
  960 continue
  980 continue
      do 1020 i=1,lsym
      jend=nsymg(i)
      do 1020 j=1,jend
      kend=isgnum(i,j)
      jc=0
      do 1000 k=1,kend
      if (isgind(i,j,k).le.ntot) go to 1000
      jc=jc+1
 1000 continue
      isgnum(i,j)=isgnum(i,j)-jc
 1020 continue
      write (iwr,46) ic,nmaxb,ntot
      go to 820
 1040 continue
 1060 continue
 1080 continue
      isum=0
      do 1120 i=1,ntot
      ii=idorb(i)
      do 1120 j=1,i
      jj=idorb(j)
      ijm=imul(ii,jj)
      do 1120 k=1,i
      kk=idorb(k)
      le=k
      if (i.eq.k) le=j
      do 1100 l=1,le
      ll=idorb(l)
      if (ijm.ne.imul(kk,ll)) go to 1100
      isum=isum+1
 1100 continue
 1120 continue
      if (isum.le.nmaxe) go to 1180
      ntot=ntot-1
c
c     check for degeneracy
c
      do 1140 nn=1,6
      if ( dabs(epsil(ntot)-epsil(ntot+1)).gt.0.1d-04) go to 1160
      ntot=ntot-1
 1140 continue
 1160 continue
      go to 1080
 1180 continue
      nmv=isum
      do 1220 i=1,lsym
      jend=nsymg(i)
      do 1220 j=1,jend
      kend=isgnum(i,j)
      ic=0
      do 1200 k=1,kend
      if (isgind(i,j,k).le.ntot) go to 1200
      ic=ic+1
 1200 continue
      isgnum(i,j)=isgnum(i,j)-ic
 1220 continue
      if (selctb.eq.zsel) go to 1260
      call caserr('invalid scftype for tda')
 1260 write(iwr,27)
      nvirt=ntot-nocc
      do 1445 i=1,norb
      ik=i
      if (indorb(i).gt.ntot) go to 1446
 1445 continue
      go to 1447
 1446 continue
      norb=ik-1
      if (norb.eq.0) call caserr(
     *'integral constraints preclude band orbital inclusion')
 1447 continue
      if (selcta.eq.sselx) write (iwr,29)
      if (selcta.ne.sselx) write (iwr,30)
      if (word.eq.selct) write (iwr,23)
      if (word.ne.selct) write (iwr,24)
      if (word.eq.selct.and.numend.le.2) write (iwr,44)
      if (word.eq.selct.and.numend.gt.2) write (iwr,45)
      if (ssele.ne.tsele) go to 1450
      write (iwr,52)
      write (iwr,53) range(1),range(2)
      range(1)=range(1)/conver
      range(2)=range(2)/conver
      write (iwr,54) alimit,blimit
      alimit=alimit/conver
      blimit=blimit/conver
      go to 1455
 1450 continue
      write (iwr,55)
 1455 continue
      write (iwr,3) nocc,nvirt
      write (iwr,5)
      write (iwr,6) (epsil(i),i=1,ntot)
      write (iwr,11)
      write (iwr,6) (sigdia(i),i=1,ntot)
      write (iwr,7) norb
      write (iwr,8)
      write (iwr,9) (indorb(i),i=1,norb)
      write (iwr,1)
      write (iwr,12) nsym
      write (iwr,13)
      do 1460 i=1,nsym
      write (iwr,14) (imul(j,i),j=1,nsym)
 1460 continue
      write (iwr,15)
      if (selctc.eq.ssely) go to 1485
      do 1480 i=1,ntot
      write (iwr,16) i,idorb(i),ndeg(i)
 1480 continue
      go to 1495
 1485 continue
      do 1490 i=1,ntot
      write (iwr,51) i,idorb(i),ndeg(i),adorb(i)
 1490 continue
 1495 continue
      write (iwr,17) isum
      write (iwr,1)
      if (selctc.eq.ssely) go to 1500
      write (iwr,31)
      go to 1520
 1500 write (iwr,32)
 1520 write (iwr,33) lsym
      do 1540 i=1,lsym
      write (iwr,34) i,nsymg(i)
 1540 continue
      do 1560 i=1,lsym
      jend=nsymg(i)
      do 1560 j=1,jend
      if (iuse(i,j).ne.0) write (iwr,35) i,j
      if (iuse(i,j).eq.0) write (iwr,36) i,j
 1560 continue
      do 1580 i=1,lsym
      jend=nsymg(i)
      do 1580 j=1,jend
      write (iwr,37) i,j,isgnum(i,j)
 1580 continue
      do 1600 i=1,lsym
      jend=nsymg(i)
      do 1600 j=1,jend
      kend=isgnum(i,j)
      if(kend.ne.0)
     *write (iwr,38) i,j,(isgind(i,j,k),k=1,kend)
 1600 continue
      write (iwr,39)
      do 1620 i=1,lsym
      jend=nsymg(i)
      do 1620 j=1,jend
      kend=isgnum(i,j)
      if(kend.ne.0)
     *write (iwr,40) i,j,(jsgjnd(i,j,k),k=1,kend)
 1620 continue
      ninta=isum
      ntotx=ntot
      nvirtx=nvirt
      if (isum.le.nmaxf) go to 1740
 1640 continue
      ntotx=ntotx-1
      nvirtx=nvirtx-1
      do 1660 nn=1,6
      if ( dabs(epsil(ntotx)-epsil(ntotx+1)).gt.0.1d-04) go to 1680
      ntotx=ntotx-1
      nvirtx=nvirtx-1
 1660 continue
 1680 continue
      isum=0
      do 1720 i=1,ntotx
      ii=idorb(i)
      do 1720 j=1,i
      jj=idorb(j)
      ijm=imul(ii,jj)
      do 1720 k=1,i
      kk=idorb(k)
      le=k
      if (i.eq.k) le=j
      do 1700 l=1,le
      ll=idorb(l)
      if (ijm.ne.imul(kk,ll)) go to 1700
      isum=isum+1
 1700 continue
 1720 continue
      if (isum.gt.nmaxf) go to 1640
 1740 continue
      nintb=isum
      do 1745 i=1,norb
      ik=i
      if (indorb(i).gt.ntotx) go to 1746
 1745 continue
      go to 1747
 1746 continue
      norb=ik-1
      if (norb.eq.0) call caserr(
     *'integral constraints preclude band orbital inclusion')
 1747 continue
      write (iwr,10) nvirtx,ntotx
      write (iwr,22) nintb
      if (selctd.eq.sseld) write (iwr,41)
      if (selctd.ne.sseld) write (iwr,42)
      do 1860 isym=1,nsym
      do 1760 i=1,ntot
      if (idorb(i).eq.isym) go to 1780
 1760 continue
      go to 1860
 1780 continue
      jend=nsymg(isym)
      do 1840 j=1,jend
      kend=isgnum(isym,j)
      do 1800 k=1,kend
      kkeep=k
      if (isgind(isym,j,k).gt.ntotx) go to 1820
 1800 continue
      kkeep=kkeep+1
 1820 continue
      ksgnum(isym,j)=kkeep-1
 1840 continue
 1860 continue
      return
      end
      subroutine insym (adorb,nmaxi,iwr)
c
      implicit REAL  (a-h,o-z)
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
      character*8 adorb,aka,akb,akc,bka,bkb
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      dimension adorb(100)
c
    1 format (1x,
     +'error in symmetry specification, more than two nondegenerate'/
     +1x,'representations of the degenerate group are correlated'/1x,
     *'with the same representation of the nondegenerate subgroub'/
     +1x,'this occurs for symmetry ',i3,' and orbitals ',2i5)
    2 format (' error in symmetry specification'/
     +1x,'more than two doubly degenerate representations of the'/
     +1x,'degenerate group are correlated with the same two identical'/
     +1x,'representations in the nondegenerate subgroup'/
     +1x,'this occurs for symmetry ',i3,' and orbitals ',2i5)
    3 format (' error in symmetry specification'/
     +1x,'more than three doubly degenerate representations of the'/
     +1x,'degenerate group are correlated with the same two different'/
     +1x,'representations in the nondegenerate subgroup'/
     +1x,'this occurs for symmetry ',i3,' and orbitals ',2i5)
    4 format (' error in symmetry specification for symmetry ',i3/
     *1x,'two or more components of a triply degenerate representation'/
     +1x,'are correlated with the same nondegenerate representation.'/
     *' this occurs for orbital ',i3)
    5 format (' error in symmetry specification'/
     +1x,'more than two triply degenerate representations of the'/
     +1x,'degenerate group are correlated with the same representation'/
     +1x,'in the nondegenerate subgroup'/
     +1x,'this occurs for symmetry ',i3,' and orbitals ',2i5)
    6 format (' error in symmetry specification'/
     +1x,'in the analysis it is found for symmetry ',i3/
     +1x,'that the number of subparts is ',i3,
     +' contrary to the input ',i3)
c   7 format ('   the number of orbitals belonging to symmetry species
c    * ',i3,'   and subpart   ',i3,'   is   ',i3)
c   8 format (/'   the orbital indices belonging to symmetry species   '
c    *,i3,'   and subpart   ',i3,'   are  :  '//3x,30i4)
    9 format (' for symmetry ',i3,' the number of subparts ',i3,
     *' exceeds the maximum allowed one ',i3)
c
c
      do 2000 isym=1,nsym
      jp=nsymg(isym)
c
c     investigation of the nondegenerate representations of the
c     degenerate group
c
      do 20 i=1,ntot
      if (idorb(i).eq.isym) go to 40
   20 continue
      go to 2000
   40 continue
      jcase=0
      icase=0
      icount=0
      ic=0
      jka=0
      aka=' '
      do 100 i=1,ntot
      if (idorb(i).ne.isym) go to 100
      if (ndeg(i).ne.1) go to 100
      icase=icase+1
      if (icase.eq.1) aka=adorb(i)
      do 80 j=i,ntot
      if (idorb(j).ne.isym) go to 80
      if (ndeg(j).ne.1) go to 80
      if (adorb(i).eq.adorb(j)) go to 60
      jcase=jcase+1
      if (jcase.eq.1) jka=j
      go to 80
   60 continue
      ic=ic+1
      if (ic.eq.1) icount=icount+1
      isgind(isym,icount,ic)=j
   80 continue
      if (ic.gt.0) go to 110
  100 continue
      go to 120
  110 continue
      isgnum(isym,icount)=ic
  120 continue
      if (jka.eq.0) go to 300
      icount=icount+1
      ic=0
      do 200 i=jka,ntot
      if (idorb(i).ne.isym) go to 200
      if (ndeg(i).ne.1) go to 200
      if (adorb(i).eq.aka) go to 200
      do 180 j=i,ntot
      if (idorb(j).ne.isym) go to 180
      if (ndeg(j).ne.1) go to 180
      if (adorb(j).eq.aka) go to 180
      if (adorb(i).eq.adorb(j)) go to 160
      write (iwr,1) isym,i,j
      goto 3000
  160 continue
      ic=ic+1
      isgind(isym,icount,ic)=j
  180 continue
      if (ic.gt.0) go to 220
  200 continue
      go to 300
  220 isgnum(isym,icount)=ic
  300 continue
c
c     investigation the doubly degenerate representations of the
c     degenerate group whose components are both correlated with the
c     same representation of the nondegenerate subgroup.
c
      icase=0
      jcase=0
      ic=0
      aka=' '
      akb=' '
      do 400 i=1,ntot
      if (idorb(i).ne.isym) go to 400
      if (ndeg(i).ne.2) go to 400
      ma=indeg(i,1)
      mma=indeg(i,2)
      if (idorb(ma).ne.idorb(mma)) go to 400
      if (i.eq.1) go to 320
      if ( dabs(epsil(i)-epsil(i-1)).lt.0.1d-04) go to 400
  320 continue
      icase=icase+1
      if (icase.ne.1) go to 330
      aka=adorb(ma)
      akb=adorb(mma)
  330 continue
      do 380 j=i,ntot
      if (idorb(j).ne.isym) go to 380
      if (ndeg(j).ne.2) go to 380
      mb=indeg(j,1)
      mmb=indeg(j,2)
      if (idorb(mb).ne.idorb(mmb)) go to 380
      if (idorb(ma).ne.idorb(mb)) go to 380
      if (j.eq.1) go to 340
      if ( dabs(epsil(j)-epsil(j-1)).lt.0.1d-04) go to 380
  340 continue
      if ((adorb(ma).eq.adorb(mb).and.adorb(mma).eq.adorb(mmb)).or.
     *(adorb(ma).eq.adorb(mmb).and.adorb(mma).eq.adorb(mb))) go to 350
      jcase=jcase+1
      if (jcase.eq.1) jka=j
      go to 380
  350 continue
      ic=ic+1
      if (ic.eq.1) icount=icount+1
      if (adorb(ma).eq.adorb(mb).and.adorb(mma).eq.adorb(mmb)) go to 360
      isgind(isym,icount,ic)=mmb
      isgind(isym,icount+1,ic)=mb
      go to 380
  360 isgind(isym,icount,ic)=mb
      isgind(isym,icount+1,ic)=mmb
  380 continue
      if (ic.eq.0) go to 400
      isgnum(isym,icount)=ic
      isgnum(isym,icount+1)=ic
      go to 420
  400 continue
      go to 600
  420 continue
      icount=icount+1
      if (jka.eq.0) go to 600
      icount=icount+1
      ic=0
      do 500 i=jka,ntot
      if (idorb(i).ne.isym) go to 500
      if (ndeg(i).ne.2) go to 500
      ma=indeg(i,1)
      mma=indeg(i,2)
      if (idorb(ma).ne.idorb(mma)) go to 500
      if (adorb(ma).eq.aka.or.adorb(ma).eq.akb) go to 500
      if (i.eq.1) go to 430
      if ( dabs(epsil(i)-epsil(i-1)).lt.0.1d-04) go to 500
  430 continue
      do 480 j=i,ntot
      if (idorb(j).ne.isym) go to 480
      if (ndeg(j).ne.2) go to 480
      mb=indeg(j,1)
      mmb=indeg(j,2)
      if (idorb(mb).ne.idorb(mmb)) go to 480
      if (idorb(ma).ne.idorb(mb)) go to 480
      if (j.eq.1) go to 440
      if ( dabs(epsil(j)-epsil(j-1)).lt.0.1d-04) go to 480
  440 continue
      if (adorb(mb).eq.aka.or.adorb(mb).eq.akb) go to 480
      if ((adorb(ma).eq.adorb(mb).and.adorb(mma).eq.adorb(mmb)).or.
     *(adorb(ma).eq.adorb(mmb).and.adorb(mma).eq.adorb(mb))) go to 450
      write (iwr,2) isym,i,j
      go to 3000
  450 continue
      ic=ic+1
      if (adorb(ma).eq.adorb(mb).and.adorb(mma).eq.adorb(mmb)) go to 460
      isgind(isym,icount,ic)=mmb
      isgind(isym,icount+1,ic)=mb
      go to 480
  460 isgind(isym,icount,ic)=mb
      isgind(isym,icount+1,ic)=mmb
  480 continue
      if (ic.eq.0) go to 500
      isgnum(isym,icount)=ic
      isgnum(isym,icount+1)=ic
      go to 520
  500 continue
      go to 600
  520 continue
      icount=icount+1
  600 continue
c
c     investigation of the doubly degenerate representations of the
c     degenerate group whose components are correlated with different
c     representations of the nondegenerate subgroup.
c
      icase=0
      jcase=0
      ic=0
      aka=' '
      akb=' '
      ika=0
      ikb=0
      do 700 i=1,ntot
      if (ndeg(i).ne.2) go to 700
      ma=indeg(i,1)
      mma=indeg(i,2)
      if (idorb(ma).eq.idorb(mma)) go to 700
      if (idorb(ma).eq.isym.or.idorb(mma).eq.isym) go to 610
      go to 700
  610 continue
      if (i.eq.1) go to 620
      if ( dabs(epsil(i)-epsil(i-1)).lt.0.1d-04) go to 700
  620 continue
      icase=icase+1
      if (icase.ne.1) go to 630
      aka=adorb(ma)
      akb=adorb(mma)
      ika=idorb(ma)
      ikb=idorb(mma)
  630 continue
      do 680 j=i,ntot
      if (ndeg(j).ne.2) go to 680
      mb=indeg(j,1)
      mmb=indeg(j,2)
      if (idorb(mb).eq.idorb(mmb)) go to 680
      if (idorb(mb).eq.isym.or.idorb(mmb).eq.isym) go to 640
      go to 680
  640 continue
      if (j.eq.1) go to 650
      if ( dabs(epsil(j)-epsil(j-1)).lt.0.1d-04) go to 680
  650 continue
      if ((idorb(ma).eq.idorb(mb).and.idorb(mma).eq.idorb(mmb)).or.
     *(idorb(ma).eq.idorb(mmb).and.idorb(mma).eq.idorb(mb))) go to 660
      go to 680
  660 continue
      if ((adorb(ma).eq.adorb(mb).and.adorb(mma).eq.adorb(mmb)).or.
     *(adorb(ma).eq.adorb(mmb).and.adorb(mma).eq.adorb(mb))) go to 670
      jcase=jcase+1
      if (jcase.eq.1) jka=j
      go to 680
  670 continue
      ic=ic+1
      if (ic.eq.1) icount=icount+1
      if (idorb(mb).eq.isym) go to 675
      isgind(isym,icount,ic)=mmb
      go to 680
  675 isgind(isym,icount,ic)=mb
  680 continue
      if (ic.eq.0) go to 700
      isgnum(isym,icount)=ic
      go to 710
  700 continue
      go to 1000
  710 continue
      if (jka.eq.0) go to 1000
      icount=icount+1
      ic=0
      kcase=0
      lcase=0
      bka=' '
      bkb=' '
      jkb=0
      do 800 i=jka,ntot
      if (ndeg(i).ne.2) go to 800
      ma=indeg(i,1)
      mma=indeg(i,2)
      if (idorb(ma).eq.idorb(mma)) go to 800
      if (idorb(ma).eq.isym.or.idorb(mma).eq.isym) go to 715
      go to 800
  715 continue
      if ((idorb(ma).eq.ika.and.idorb(mma).eq.ikb).or.(idorb(mma).eq.
     *ika.and.idorb(ma).eq.ikb)) go to 718
      go to 800
  718 continue
      if (i.eq.1) go to 720
      if ( dabs(epsil(i)-epsil(i-1)).lt.0.1d-04) go to 800
  720 continue
      if ((adorb(ma).eq.aka.and.adorb(mma).eq.akb).or.(adorb(mma).eq.
     *aka.and.adorb(ma).eq.akb)) go to 800
      kcase=kcase+1
      if (kcase.ne.1) go to 730
      bka=adorb(ma)
      bkb=adorb(mma)
  730 continue
      do 780 j=i,ntot
      if (ndeg(j).ne.2) go to 780
      mb=indeg(j,1)
      mmb=indeg(j,2)
      if (idorb(mb).eq.idorb(mmb)) go to 780
      if (idorb(mb).eq.isym.or.idorb(mmb).eq.isym) go to 740
      go to 780
  740 continue
      if (j.eq.1) go to 745
      if ( dabs(epsil(j)-epsil(j-1)).lt.0.1d-04) go to 780
  745 continue
      if ((idorb(ma).eq.idorb(mb).and.idorb(mma).eq.idorb(mmb)).or.
     *(idorb(ma).eq.idorb(mmb).and.idorb(mma).eq.idorb(mb))) go to 750
      go to 780
  750 continue
      if ((adorb(mb).eq.aka.and.adorb(mmb).eq.akb).or.
     *(adorb(mmb).eq.aka.and.adorb(mb).eq.akb)) go to 780
      if ((adorb(ma).eq.adorb(mb).and.adorb(mma).eq.adorb(mmb)).or.
     *(adorb(ma).eq.adorb(mmb).and.adorb(mma).eq.adorb(mb))) go to 760
      lcase=lcase+1
      if (lcase.eq.1) jkb=j
      go to 780
  760 continue
      ic=ic+1
      if (idorb(mb).eq.isym) go to 770
      isgind(isym,icount,ic)=mmb
      go to 780
  770 isgind(isym,icount,ic)=mb
  780 continue
      if (ic.eq.0) go to 800
      isgnum(isym,icount)=ic
      go to 810
  800 continue
      go to 1000
  810 continue
      if (jkb.eq.0) go to 1000
      icount=icount+1
      ic=0
      do 900 i=jkb,ntot
      if (ndeg(i).ne.2) go to 900
      ma=indeg(i,1)
      mma=indeg(i,2)
      if (idorb(ma).eq.idorb(mma)) go to 900
      if (idorb(ma).eq.isym.or.idorb(mma).eq.isym) go to 815
      go to 900
  815 continue
      if ((idorb(ma).eq.ika.and.idorb(mma).eq.ikb).or.
     *(idorb(mma).eq.ika.and.idorb(ma).eq.ikb)) go to 818
      go to 900
  818 continue
      if (i.eq.1) go to 820
      if ( dabs(epsil(i)-epsil(i-1)).lt.0.1d-04) go to 900
  820 continue
      if ((adorb(ma).eq.aka.and.adorb(mma).eq.akb).
     *or.(adorb(mma).eq.aka.and.adorb(ma).eq.akb)) go to 900
      if ((adorb(ma).eq.bka.and.adorb(mma).eq.bkb).
     *or.(adorb(mma).eq.bka.and.adorb(ma).eq.bkb)) go to 900
      do 880 j=i,ntot
      if (ndeg(j).ne.2) go to 880
      mb=indeg(j,1)
      mmb=indeg(j,2)
      if (idorb(mb).eq.idorb(mmb)) go to 880
      if (idorb(mb).eq.isym.or.idorb(mmb).eq.isym) go to 840
      go to 880
  840 continue
      if (j.eq.1) go to 845
      if ( dabs(epsil(j)-epsil(j-1)).lt.0.2d-04) go to 880
  845 continue
      if ((idorb(ma).eq.idorb(mb).and.idorb(mma).eq.idorb(mmb)).or.
     *(idorb(ma).eq.idorb(mmb).and.idorb(mma).eq.idorb(mb))) go to 850
      go to 880
  850 continue
      if ((adorb(mb).eq.aka.and.adorb(mmb).eq.akb).or.
     *(adorb(mmb).eq.aka.and.adorb(mmb).eq.akb)) go to 880
      if ((adorb(mb).eq.bka.and.adorb(mmb).eq.bkb).or.
     *(adorb(mmb).eq.bka.and.adorb(mmb).eq.bkb)) go to 880
      if ((adorb(ma).eq.adorb(mb).and.adorb(mma).eq.adorb(mmb)).or.
     *(adorb(ma).eq.adorb(mmb).and.adorb(mma).eq.adorb(mb))) go to 860
      write (iwr,3) isym,i,j
      go to 3000
  860 continue
      ic=ic+1
      if (idorb(mb).eq.isym) go to 870
      isgind(isym,icount,ic)=mmb
      go to 880
  870 isgind(isym,icount,ic)=mb
  880 continue
      if (ic.eq.0) go to 900
      isgnum(isym,icount)=ic
      go to 1000
  900 continue
 1000 continue
c
c     investigation of the triply degenerate representations of the
c     degenerate group.
c
      do 1005 i=1,ntot
      if (ndeg(i).ne.3) go to 1005
      if (i.eq.1) go to 1001
      if ( dabs(epsil(i)-epsil(i-1)).lt.0.1d-04) go to 1005
 1001 continue
      ma=indeg(i,1)
      mb=indeg(i,2)
      mc=indeg(i,3)
      if (idorb(ma).eq.idorb(mb).or.idorb(ma).eq.idorb(mc).or.
     *idorb(mb).eq.idorb(mc)) go to 1002
      go to 1005
 1002 write (iwr,4) isym,i
      go to 3000
 1005 continue
      icase=0
      jcase=0
      kcase=0
      ic=0
      ika=0
      ikb=0
      ikc=0
      aka=' '
      akb=' '
      akc=' '
      nrep=0
      ist=1
 1006 continue
      if (nrep.eq.1) ist=mc+1
      if (ist.gt.ntot) nrep=0
      if (ist.gt.ntot) go to 1210
      do 1200 i=ist,ntot
      if (ndeg(i).ne.3) go to 1200
      ma=indeg(i,1)
      mb=indeg(i,2)
      mc=indeg(i,3)
      if (idorb(ma).eq.isym.or.idorb(mb).eq.isym.or.idorb(mc).eq.isym)
     * go to 1010
      go to 1200
 1010 continue
      if (i.eq.1) go to 1020
      if ( dabs(epsil(i)-epsil(i-1)).lt.0.1d-04) go to 1200
 1020 continue
      icase=icase+1
      if (icase.ne.1) go to 1030
      aka=adorb(ma)
      akb=adorb(mb)
      akc=adorb(mc)
      ika=idorb(ma)
      ikb=idorb(mb)
      ikc=idorb(mc)
 1030 continue
      if (nrep.eq.0) go to 1034
      jc=0
      do 1032 k=1,3
      ka=indeg(i,k)
      bka=adorb(ka)
      if (bka.eq.aka) jc=jc+1
      if (bka.eq.akb) jc=jc+1
      if (bka.eq.akc) jc=jc+1
 1032 continue
      if (jc.eq.3) go to 1200
      kcase=kcase+1
      if (kcase.ne.1) go to 1034
      aka=adorb(ma)
      akb=adorb(mb)
      akc=adorb(mc)
      ika=idorb(ma)
      ikb=idorb(mb)
      ikc=idorb(mc)
 1034 continue
      do 1180 j=i,ntot
      if (ndeg(j).ne.3) go to 1180
      mma=indeg(j,1)
      mmb=indeg(j,2)
      mmc=indeg(j,3)
      if (idorb(mma).eq.isym.or.idorb(mmb).eq.isym.or.idorb(mmc).eq.isym
     *) go to 1040
      go to 1180
 1040 continue
      if (j.eq.1) go to 1050
      if ( dabs(epsil(j)-epsil(j-1)).lt.0.1d-04) go to 1180
 1050 continue
      jc=0
      do 1070 k=1,3
      ka=indeg(i,k)
      kka=idorb(ka)
      do 1060 l=1,3
      la=indeg(j,l)
      if (kka.eq.idorb(la)) jc=jc+1
 1060 continue
 1070 continue
      if (jc.lt.3) go to 1180
      if (jc.eq.3) go to 1080
      write (iwr,4) isym,j
      go to 3000
 1080 continue
      jc=0
      do 1100 k=1,3
      ka=indeg(i,k)
      bka=adorb(ka)
      do 1090 l=1,3
      la=indeg(j,l)
      if (bka.eq.adorb(la)) jc=jc+1
 1090 continue
 1100 continue
      if (jc.eq.3) go to 1110
      jcase=jcase+1
      if (jcase.eq.1) jka=j
      go to 1180
 1110 continue
      ic=ic+1
      if (ic.eq.1) icount=icount+1
      if (idorb(mma).eq.isym) go to 1130
      if (idorb(mmb).eq.isym) go to 1120
      isgind(isym,icount,ic)=mmc
      go to 1180
 1120 continue
      isgind(isym,icount,ic)=mmb
      go to 1180
 1130 continue
      isgind(isym,icount,ic)=mma
 1180 continue
      if (ic.eq.0) go to 1200
      isgnum(isym,icount)=ic
      nrep=nrep+1
      if (nrep.eq.1) go to 1182
      jka=0
      nrep=0
      go to 1900
 1182 continue
      ic=0
      if (jka.eq.0) go to 1006
      go to 1210
 1200 continue
      nrep=0
      go to 1900
 1210 continue
      if (jka.eq.0) go to 1900
      icount=icount+1
      ic=0
      do 1400 i=jka,ntot
      if (ndeg(i).ne.3) go to 1400
      ma=indeg(i,1)
      mb=indeg(i,2)
      mc=indeg(i,3)
      if (idorb(ma).eq.isym.or.idorb(mb).eq.isym.or.idorb(mc).eq.isym)
     * go to 1220
      go to 1400
 1220 continue
      jc=0
      do 1222 k=1,3
      ka=indeg(i,k)
      kka=idorb(ka)
      if (kka.eq.ika) jc=jc+1
      if (kka.eq.ikb) jc=jc+1
      if (kka.eq.ikc) jc=jc+1
 1222 continue
      if (jc.ne.3) go to 1400
      if (i.eq.1) go to 1230
      if ( dabs(epsil(i)-epsil(i-1)).lt.0.1d-04) go to 1400
 1230 continue
      jc=0
      do 1240 k=1,3
      ka=indeg(i,k)
      bka=adorb(ka)
      if (bka.eq.aka) jc=jc+1
      if (bka.eq.akb) jc=jc+1
      if (bka.eq.akc) jc=jc+1
 1240 continue
      if (jc.eq.3) go to 1400
      do 1380 j=i,ntot
      if (ndeg(j).ne.3) go to 1380
      mma=indeg(j,1)
      mmb=indeg(j,2)
      mmc=indeg(j,3)
      if (idorb(mma).eq.isym.or.idorb(mmb).eq.isym.or.idorb(mmc).eq.isym
     *) go to 1250
      go to 1380
 1250 continue
      if (j.eq.1) go to 1260
      if ( dabs(epsil(j)-epsil(j-1)).lt.0.1d-04) go to 1380
 1260 continue
      jc=0
      do 1280 k=1,3
      ka=indeg(i,k)
      kka=idorb(ka)
      do 1270 l=1,3
      la=indeg(j,l)
      if (kka.eq.idorb(la)) jc=jc+1
 1270 continue
 1280 continue
      if (jc.lt.3) go to 1380
      if (jc.eq.3) go to 1290
      write (iwr,4) isym,j
      go to 3000
 1290 continue
      jc=0
      do 1300 l=1,3
      la=indeg(j,l)
      bka=adorb(la)
      if (bka.eq.aka) jc=jc+1
      if (bka.eq.akb) jc=jc+1
      if (bka.eq.akc) jc=jc+1
 1300 continue
      if (jc.eq.3) go to 1380
      jc=0
      do 1320 k=1,3
      ka=indeg(i,k)
      bka=adorb(ka)
      do 1310 l=1,3
      la=indeg(j,l)
      if (bka.eq.adorb(la)) jc=jc+1
 1310 continue
 1320 continue
      if (jc.eq.3) go to 1330
      write (iwr,5) isym,i,j
      go to 3000
 1330 continue
      ic=ic+1
      if (idorb(mma).eq.isym) go to 1350
      if (idorb(mmb).eq.isym) go to 1340
      isgind(isym,icount,ic)=mmc
      go to 1380
 1340 isgind(isym,icount,ic)=mmb
      go to 1380
 1350 isgind(isym,icount,ic)=mma
 1380 continue
      if (ic.eq.0) go to 1400
      isgnum(isym,icount)=ic
      go to 1900
 1400 continue
 1900 continue
      if (jka.ne.0) go to 1910
      if (nrep.eq.0) go to 1910
      ic=0
      go to 1006
 1910 continue
      if (icount.eq.jp) go to 1920
      write (iwr,6) isym,icount,jp
      go to 3000
 1920 if (icount.le.nmaxi) go to 2000
      write (iwr,9) isym,icount,nmaxi
      go to 3000
 2000 continue
      return
3000  call caserr('error in tda symmetry specification')
      return
      end
      subroutine tda(q)
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/statis)
INCLUDE(common/segm)
INCLUDE(common/iofile)
      dimension q(*)
      begin=cpulft(1)
      write(iwr,1100)
 1100 format(/1x,104('-')/
     *40x,28('*')/
     *40x,'ionization spectra tda module'/
     *40x,28('*')//)
c
c     allocate core
c     allocate all available memory
c
      i10 = igmem_alloc_all(lword)
c
      if(lword.lt.1)call caserr(
     *'insufficient memory for tda module')
      call initda(q(i10),q(i10),lword)
      at=cpulft(1)
      write(iwr,1300) at
 1300 format(//5x,'end of tda module at',f8.2,' secs'//)
      call clredx
      call timana(16)
c
c     reset core allocation
c
      call gmem_free(i10)
c
      return
      end
      subroutine initda(q,iq,lword)
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxorb1=maxorb+1)
INCLUDE(common/machin)
      common /tdaa/nocc,nunocc,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *nsym,imul(10,10),idorb(100),norb,indorb(30)
INCLUDE(common/discc)
INCLUDE(common/runlab)
INCLUDE(common/iofile)
      common/junk/sigdia(mxtda4),epsil(mxtda4),
     +pop(maxorb),potn,core,ncolo(3),
     *ncore,mapcie(maxorb),mapcif(maxorb),nact,mapaie(maxorb),
     *mapaei(maxorb),iqsec,
     *iworka(maxorb),iworkb(maxorb),nactc(5),isecor,
     *norba,norbaa(maxorb),norbb,norb3(maxorb),ntoto,iwork(maxorb),noscr
     *,range(2),alimit,blimit,isymm,nx,lsym,nsymm(8),iuse(8,6),jdorb(100
     *),ncont,nmaxa,nmaxb,nmaxc,nmaxd,nmaxe,nmaxf,nmaxg,nmaxh,
     * nmaxi,nmaxj,
     *evalue(maxorb),eocc(mxorb1),nbas,newb,ncol,ieval,ipop,ispp
     *,nirr,mult(8,8),isymao(maxorb),isymmo(maxorb),nsp,mapie(maxorb)
INCLUDE(common/filel)
INCLUDE(common/restar)
      common/junkc/zopi(13),zztag(100),
     *zjob,zdate,ztime,zprog,ztype,zspace(14),ztext(10)
      common/tdab/conver,cons
      common/craypk/timfac,time1
INCLUDE(common/restri)
      dimension q(*),iq(*)
      data xstar,thresh,two /'*',1.0d-5,2.0d0/
c     data zblank,xdash,ydomo,yvmo /' ','-','domo','vmo'/
      data m0,m10,m117/0,10,117/
      data m51,m29/51,29/
      tim1=cpulft(1)
      nav = lenwrd()
c     lwor4=lword*nav
      conver=27.2116d0
      cons=conver*conver
      call secget(isect(470),1005,iblka)
      call readi(iworka,mach(14)*nav,iblka,idaf)
      call readis(norba,mach(16)*nav,idaf)
      call rdchrs(zopi,m117,idaf)
      call secget(isecor,1004,iblka)
      call rdedx(pop,mach(15),iblka,idaf)
      write(iwr,5004) yed(idaf),ibl3d,isecor
 5004 format(/5x,'dumpfile resides on ',a4,' at block ',i6//
     *          5x,'core hamiltonian to be retrieved from section ',i4/)
      write(iwr,2)nact,ncore,potn,core
    2 format(5x,'header block information :'/
     *       5x,'number of active orbitals ',i4/
     *       5x,'number of core   orbitals ',i4/
     *       5x,'nuclear repulsion energy ',e21.12/
     *       5x,'core energy              ',e21.12///)
      if(nact.gt.1) go to 10005
80000 call caserr('invalid number of active orbitals')
80010 call caserr('parameter error in tda preprocessor')
80020 call caserr('invalid number of satellite bands requested')
80030 call caserr('insufficient memory for tda module')
10005 do 82003 i=1,nact
      j=mapaie(i)
82003 mapaei(j)=i
      lfile=m6file
      do 60001 i=1,lfile
      lotape(i)=m6tape(i)
      liblk  (i)=m6blk(i)
60001 llblk  (i)=liblk(i)-m6last(i)
      ntot=ntoto
      if(ntot.ne.0)go to 5005
      ntot=nact
      do 5010 i=1,ntot
      iworka(i)=i
5010  mapie(i)=i
      go to 8460
5005  do 82004 i=1,ntot
82004 iworka(i)=mapaei(iwork(i))
      call rddir(ntot,iworka)
      write(iwr,15002)
      write(iwr,9) (iworka(i),i=1,ntot)
      if(ntot.le.1.or.ntot.gt.100) go to 80000
      if(ntot.gt.nact.or.iworka(ntot).gt.nact) go to 80000
      do 7010 i=1,nact
 7010 mapie(i)=0
      do 7011 i=1,ntot
 7011 mapie(iworka(i))=i
 8460 norb=norba
      if(norb)60006,80020,8463
8463  call rdba(nact,norb,indorb,norbaa,mapaei,iworka,mapie)
      write(iwr,7) norb
      write(iwr,9) (indorb(i),i=1,norb)
    7 format(//5x,'satellite bands are to be calculated ',
     *                'for the following',i4,'  orbitals :')
60006 call secget(isect(490),m51,iblka)
      osym=.true.
      call readi(nirr,mach(13)*nav,iblka,idaf)
      call setsto(ntot,m0,idorb)
      call secget(iqsec,3,iblka)
      call rdchr(zjob,m29,iblka,idaf)
      call reads(evalue,mach(8),idaf)
      write(iwr,89005) iqsec,ztype,ztime,zdate,zjob,ztext,nbas,ncol
89005 format(/5x,'scf mo specifications restored from section ',i4,
     *          ' of dumpfile'//
     *       5x,'header block information :'/
     *       5x,a7,' vectors created at ',
     *          a8,' on ',a8,' in the job ',a8/
     *       5x,'with the title :',5x,10a8/
     *       5x,'no. of gtos        ',i4/
     *       5x,'no. of scf vectors ',i4//)
       call symvec(q,isymao,isymmo,nbas,ncol,iblka)
      ic=0
      do 44 i=1,nact
      ida=mapie(i)
      if(ida.eq.0)go to 44
      ic=ic+1
      if(idorb(ida).ne.0)go to 80010
      idorb(ida)=isymmo(mapaie(ida))
 44   continue
      if(ic.ne.ntot)go to 80010
      nocc=0
      nunocc=0
      do 89006 i=1,ntot
      jscf=mapaie(iworka(i))
      epsil(i)=evalue(jscf)
      qio=eocc(jscf)
      if( dabs(qio-two).le.thresh) nocc=nocc+1
      if(qio.le.thresh) nunocc=nunocc+1
89006 continue
      if(nocc+nunocc.ne.ntot) go to 80010
      if(nunocc.le.0) go to 80010
      nstart=nocc+1
      if(norb.gt.0)go to 89100
      norb=nocc
      do 89099 i=1,nocc
89099 indorb(i)=i
    3 format(  5x,'case :',5x,10a8
     *       //5x,'number of active orbitals         ',i4/
     *        /5x,'number of doubly occupied orbitals',i4/
     *        /5x,'number of virtual orbitals        ',i4)
    5 format (//50x,'scf orbital energies (au)')
    6 format(/10x,8f14.7)
    9 format(/20(3x,i3))
   13 format(5x,'no symmetry is assumed for the present case'//)
  144 format(4x,120a1)
   17 format (//5x,i7,' integrals are to be read in and processed')
   24 format(//5x,'insufficient storage available, increase by',
     *      i10,' words'//)
   25 format(/5x,'main core not used:',2x,i10,' words'//)
89100 write(iwr,144) (xstar,i=1,93)
      write(iwr,3)ztitle,ntot,nocc,nunocc
      write(iwr,144) (xstar,i=1,93)
15002 format(///5x,'the following mo''s included in active list :'/)
      write(iwr,5)
      write(iwr,6) (epsil(i),i=1,ntot)
      if(.not.osym) write(iwr,13)
      temp=epsil(1)*conver
      if( dabs(range(2)).lt.1.0d-4)range(2)= dabs(temp)*1.1d0
c
      call intda(nmv,iwr)
c
c     dynamic store for tammda
c
      ixx=ntot*(ntot+1)/2
      iyy=ixx*(ntot+ntot+1)/3
      iln=ixx
      if(iln.le.3) iln=iln+1
      ju1=1
      ju2=ju1+ntot
      ju3=ju2+ntot
      ju4=ju3+ixx
      ju5=ju4+iyy
      nmaxbg=nmaxb+nmaxg
      nmaxcj=nmaxc+nmaxj
      nmaxgj=nmaxg+nmaxj
_IF(cray)
      ju6=ju5+nmaxbg
      ju7=ju6+nmaxbg
      ju8=ju7+nmaxbg
      ju9=ju8+nmaxbg
      ju10=ju9+8*nmaxi*nmaxj
      ju11=ju10+8*nmaxi*nmaxcj
_ELSE
      ju6=ju5+(nmaxbg-1)/nav+1
      ju7=ju6+(nmaxbg-1)/nav+1
      ju8=ju7+(nmaxbg-1)/nav+1
      ju9=ju8+(nmaxbg-1)/nav+1
      ju10=ju9+(8*nmaxi*nmaxj-1)/nav+1
      ju11=ju10+(8*nmaxi*nmaxcj-1)/nav+1
_ENDIF
      juaa=(ju11-1)/nav+1
      jubb=juaa+nmaxa*nmaxa
      juv =jubb+nmaxa*nmaxa
      juga=juv+nmaxe
      justa=juga+nmaxg
      last1=justa+nmaxg*nmaxc
      jug=juv+nmaxf
      jugvec=jug+nmaxb*nmaxb
      juaux=jugvec+nmaxb
      last2=juaux +nmaxb*nmaxc
c
      last=max(last1,last2)
      if(lword.ge.last)go to 96
      need=last-lword
      write(iwr,24) need
      go to 80030
   96 nfree=lword-last
      write(iwr,25) nfree
      write(iwr,1)
 1    format(/1x,'transformed integral files'/1x,26('-'))
      call filprn(m6file,m6blk,m6last,m6tape)
      write(iwr,17) nmv
      call inhelp(iq(ju1),iq(ju2),iq(ju3),iq(ju4),ntot,idorb,imul,m10)
_IF(ibm,vax)
      call linkor(iq(ju4),iq(ju3),iq(ju2),iq(ju1))
_ENDIF
      call rdvint(q(juv),iworka(ntot),idorb,mapie,imul,m10,nmv,
     * iq(ju4),iq(ju3),iq(ju2),iq(ju1),iwr)
      time1=cpulft(1)-tim1
      call tamctl(nmaxa,nmaxb,nmaxc,nmaxd,nmaxe,nmaxf,nmaxg,nmaxh,nmaxi,
     *nmaxj,nmaxcj,nmaxgj,iq(ju4),iq(ju3),iq(ju2),iq(ju1),
     * iq(ju5),iq(ju6),iq(ju7),iq(ju8),iq(ju9),iq(ju10),
     * q(juv),q(juga),q(justa),q(jug),q(jug),q(jugvec),q(juaux),
     * q(juaa),q(jubb),iwr)
      return
      end
      subroutine rdvint(v,imax,idorb,mapie,imul,ndim,nmv,
     *indvec,jndvec,kndvec,lndvec,iwr)
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/blkin/gin(510),numi
INCLUDE(common/atmblk)
_IF(ibm,vax)
      common/craypk/i205(340),j205(340),k205(340),l205(340)
_ELSE
      common/craypk/i205(1360)
_ENDIF
INCLUDE(common/filel)
      dimension v(*),mapie(*),imul(ndim,*),idorb(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
c
      ic=0
      call vclr(v,1,nmv)
      call setsto(1360,0,i205)
      itot=0
      do 8 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      call find(iunit)
      lbl=llblk(ifile)
    6 lbl=lbl+1
      call get(gin(1),mw)
      if(mw.eq.0) go to 8
      itot=itot+numi
      if(lbl.ne.0) call find(iunit)
_IF(ibm,vax)
      call upak8v(gin(num2e+1),i205)
_ELSE
      int4=1
      call unpack(gin(num2e+1),lab816,i205,numlab)
_ENDIF
      do 5 num=1,numi
_IF(ibm,vax)
      i=i205(num)
_ELSEIF(littleendian)
      i=i205(int4+1)
_ELSE
      i=i205(int4  )
_ENDIF
      if(i.le.imax) go to 10
      if(ic.ge.nmv) go to 12
      go to 7
   10 i=mapie(i)
      if(i.eq.0) go to 5
_IF(ibm,vax)
      j=mapie(j205(num))
_ELSEIF(littleendian)
      j=mapie(i205(int4))
_ELSE
      j=mapie(i205(int4+1))
_ENDIF
      if(j.eq.0) go to 5
_IF(ibm,vax)
      k=mapie(k205(num))
_ELSEIF(littleendian)
      k=mapie(i205(int4+3))
_ELSE
      k=mapie(i205(int4+2))
_ENDIF
      if(k.eq.0) go to 5
_IF(ibm,vax)
      l=mapie(l205(num))
_ELSEIF(littleendian)
      l=mapie(i205(int4+2))
_ELSE
      l=mapie(i205(int4+3))
_ENDIF
      if(l.eq.0) go to 5
      ijm=imul(idorb(i),idorb(j))
      klm=imul(idorb(k),idorb(l))
      if(ijm.ne.klm) go to 5
      ic=ic+1
      ind=indvx(i,j,k,l)
      v(ind)=gin(num)
_IF(ibm,vax)
    5 continue
_ELSE
    5 int4=int4+4
_ENDIF
    7 if(lbl.ne.0) go to 6
    8 continue
   12 write(iwr,999) itot
  999 format(//5x,i7,'  2-electron integrals scanned'//)
      if(ic.ne.nmv) write(iwr,1000) ic
 1000 format(//5x,'*** warning ***'/
     *         5x,'only ',i7,' symmetry allowed integrals found on',
     *            ' transformed integral file'/
     *         5x,'remaining integrals assumed zero'//)
      return
      end
      subroutine eideg( idiff,ideg,ipol,bmat,lmat,deriv,rvec,work,eigen,
     * zero,maxn,maxm,ndim,nd,mdim,ivec,ifile,iprint,neigen,iauf )
c
c***********************************************************************
c
c     this subroutine calculates the eigenvectors for the given
c     eigenvalue ipol (case of degenertae elements of c ). for ivec = 0
c     only the eigenvalue is written on file ifile. if ivec = 1 the
c     eigenvectors is also calculated and written on file ifile
c     (blocked with nd).
c
c     input
c     ipol              number of the pole of lmat
c     ideg              number of linearly independent vectors
c     idiff             difference between the degree of degeneracy
c                       of the corresponding diagonal element of cmat
c                       and ideg
c     eigen             degenerate diagonal element of cmat, which
c                       is simultaneously eigenvalue of the bordered
c                       matrix
c
c     output
c     eigenvalues and eigenvectors of the bordered matrix on file ifile
c
c***********************************************************************
c
      implicit REAL  (a-h,o-z)
      REAL  lmat
c
      dimension bmat(ndim,mdim),rvec(ndim,ndim),lmat(nd),work(ndim),
     *          deriv(nd)
c
    1 format(/1x,'the eigenvector of the ',i5,
     + 'th eigenvalue could not be calculated,'/
     +' as the dimension of "amat"  ndim = ',i5,
     *' is smaller than the number n = ',i4/
     +' of the degenerate elements of the matrix  "cmat" .')
c
      data d0,d1 /0.0d0,1.0d0/
c
c***********************************************************************
c
c     write the eigenvalue on file ifile without calculating the
c     eigenvector
c
c***********************************************************************
c
      if(ivec .eq. 1) go to 40
      do 20 i = 1,idiff
      neigen = neigen + 1
      write(ifile) eigen,eigen,eigen
   20 continue
      return
c
c***********************************************************************
c
c     determine whether a calculation of the eigenvector is possible
c
c***********************************************************************
c
c     if idiff + ideg   .gt.  ndim, all components are set equal to
c     1 and a warning is issued
c
   40 if(iauf .eq. 1) go to 1200
      maxall = idiff + ideg
      if(maxall .le. ndim) go to 100
      maxi = maxn + maxm
      maxj = maxn + 1
      do 80 j = maxj,maxall
      neigen = neigen + 1
      write(ifile) eigen,eigen,eigen
      iz = 0
      do 60 i = 1,maxi
      iz = iz + 1
      if(iz .eq. nd+1) iz = 1
      lmat(iz) = d1
      if(iz .eq. nd) write(ifile) lmat
   60 continue
      if(iz .ne. nd) write(ifile) lmat
      write(iprint,1) neigen,ndim,maxall
   80 continue
      return
c
c***********************************************************************
c
c     calculation of the eigenvector components by solving the system
c     of linear equations
c
c***********************************************************************
c
c     1. gaussian elimination procedure  (triangular form of the
c        coefficient matrix)
c
  100 do 110 i = 1,ndim
      do 110 j = 1,ndim
      rvec(i,j) = d0
  110 continue
      do 120 i = 1,maxn
      do 120 j = 1,maxall
      rvec(i,j) = bmat(i,ipol+j-1)
  120 continue
      if(maxn .eq. 1) go to 320
      do 300 i = 2,maxn
      do 290 j = i,maxn
  140 if( dabs(rvec(i-1,i-1)) .le. zero) go to 160
      xfak = rvec(j,i-1) / rvec(i-1,i-1)
      go to 240
  160 do 180 k = j,maxn
      kk = k
      if( dabs(rvec(k,i-1)) .gt. zero) go to 200
  180 continue
      go to 290
  200 do 220 k = 1,maxall
      xfak = rvec(i-1,k)
      rvec(i-1,k) = rvec(kk,k)
      rvec(kk,k) = xfak
  220 continue
      go to 140
  240 do 260 k = 1,maxall
      rvec(j,k) = rvec(j,k) - rvec(i-1,k) * xfak
  260 continue
  290 continue
  300 continue
c
c     2. separate the solutions which decouple completely
c
  320 ishift = 0
      do 380 i = 1,maxall
      xfak = d0
      deriv(i) = d1
      do 340 k = 1,maxn
      xfak = xfak + rvec(k,i)*rvec(k,i)
  340 continue
      if(xfak .gt. zero) go to 380
      do 360 j = 1,maxall
      work(j) = d0
  360 continue
      work(i) = d1
      deriv(i) = d0
      ishift = ishift + 1
      write(ifile) (work(j),j=1,maxall)
  380 continue
      iz = 0
      do 420 i = 1,maxall
      if(deriv(i) .eq. d0) go to 420
      iz = iz + 1
      do 400 j = 1,maxn
      rvec(j,iz) = rvec(j,i)
  400 continue
  420 continue
c
c     3. determine the linearly independent solutions
c
      id = idiff - ishift
      if(id .ne. 0) go to 460
      do 430 i = 1,ishift
      backspace ifile
  430 continue
      do 440 i = 1,ishift
      read(ifile) (rvec(k,i),k=1,maxall)
  440 continue
      do 450 i = 1,ishift
      backspace ifile
  450 continue
      if(ishift .eq. idiff) go to 900
  460 maxa = maxall - ishift
      do 520 i = 1,id
      do 470 j = 1,maxall
      work(j) = d0
  470 continue
      work(maxn+i) = d1
      do 500 jj = 1,maxn
      j = maxn - jj + 1
      kmin = j + 1
      xfak = d0
      do 480 k = kmin,maxa
      xfak = xfak + rvec(j,k) * work(k)
  480 continue
      work(j) = -xfak / rvec(j,j)
  500 continue
      write(ifile) (work(j),j=1,maxall)
  520 continue
c
c     4. schmidt orthogonalization of the soultions
c
      do 600 i = 1,id
      backspace ifile
  600 continue
      read(ifile) (rvec(k,1),k=1,maxall)
      dnorm = d0
      do 620 i = 1,maxa
      dnorm = dnorm + rvec(i,1) * rvec(i,1)
  620 continue
      dnorm =  dsqrt( dnorm )
      do 640 i = 1,maxa
      rvec(i,1) = rvec(i,1) / dnorm
  640 continue
      if(id .eq. 1) go to 780
      do 760 j = 2,id
      read(ifile) (rvec(k,j),k=1,maxall)
      imax = j - 1
      do 700 i = 1,imax
      xfak = d0
      do 680 k = 1,maxa
      xfak = xfak + rvec(k,j) * rvec(k,i)
  680 continue
      do 700 k = 1,maxa
      rvec(k,j) = rvec(k,j) - xfak * rvec(k,i)
  700 continue
      dnorm = d0
      do 720 k = 1,maxa
      dnorm = dnorm + rvec(k,j) * rvec(k,j)
  720 continue
      dnorm =  dsqrt( dnorm )
      do 740 k = 1,maxa
      rvec(k,j) = rvec(k,j) / dnorm
  740 continue
  760 continue
  780 if(id .eq. idiff) go to 860
      do 800 i = 1,idiff
      backspace ifile
  800 continue
      do 820 i = 1,ishift
      read(ifile) (rvec(k,i+id),k=1,maxall)
  820 continue
      do 840 i = 1,ishift
      backspace ifile
  840 continue
      go to 900
  860 do 880 i = 1,idiff
      backspace ifile
  880 continue
c
c     5. write the eigenvalues and eigenvectors on file ifile
c
  900 do 1000 i = 1,id
      neigen = neigen + 1
      write(ifile) eigen,eigen,eigen
      iz = 0
      imax = maxn + ipol - 1
      do 920 j = 1,imax
      iz = iz + 1
      if(iz .eq. nd+1) iz = 1
      lmat(iz) =  d0
      if(iz .eq. nd) write(ifile) lmat
  920 continue
      jz = 0
      do 940 j = 1,maxall
      iz = iz + 1
      if(iz .eq. nd+1) iz = 1
      if(deriv(j) .ne. d0) jz = jz + 1
      if(deriv(j) .eq. d0) lmat(iz) = d0
      if(deriv(j) .ne. d0) lmat(iz) = rvec(jz,i)
      if(iz .eq. nd) write(ifile) lmat
  940 continue
      if(ipol+maxall .eq. maxm) go to 980
      imin = ipol + maxall
      do 960 j = imin,maxm
      iz = iz + 1
      if(iz .eq. nd+1) iz = 1
      lmat(iz) = d0
      if(iz .eq. nd) write(ifile) lmat
  960 continue
  980 if(iz .ne. nd) write(ifile) lmat
 1000 continue
c
      if(ishift .eq. 0) return
      id = id  + 1
      do 1180 i = id,idiff
      neigen = neigen + 1
      write(ifile) eigen,eigen,eigen
      iz = 0
      imax = maxn + ipol - 1
      do 1100 j = 1,imax
      iz = iz + 1
      if(iz .eq. nd+1) iz = 1
      lmat(iz) =  d0
      if(iz .eq. nd) write(ifile) lmat
 1100 continue
      do 1120 j = 1,maxall
      iz = iz + 1
      if(iz .eq. nd+1) iz = 1
      lmat(iz) = rvec(j,i)
      if(iz .eq. nd) write(ifile) lmat
 1120 continue
      if(ipol+maxall .eq. maxm) go to 1160
      imin = ipol + maxall
      do 1140 j = imin,maxm
      iz = iz + 1
      if(iz .eq. nd+1) iz = 1
      lmat(iz) = d0
      if(iz .eq. nd) write(ifile) lmat
 1140 continue
 1160 if(iz .ne. nd) write(ifile) lmat
 1180 continue
c
      return
c
c***********************************************************************
c
c     special case that all coupling matrix elements are zero
c
c***********************************************************************
c
 1200 maxall = maxn + maxm
      write(ifile) eigen,eigen,eigen
      iz = 0
      do 1220 i = 1,maxall
      iz = iz + 1
      if(iz .eq. nd+1) iz = 1
      lmat(iz) = d0
      if(i .eq. ipol+maxn) lmat(iz) = d1
      if(iz .eq. nd) write(ifile) lmat
 1220 continue
      if(iz .ne. nd) write(ifile) lmat
c
      return
      end
      subroutine eigrsf (a,n,jobn,d,z,iz,wk,ier)
c
c-----------------------------------------------------------------------
c
c   purpose             - eigenvalues and (optionally) eigenvectors of
c                           a real symmetric matrix
c
c   usage               - call eigrsf(a,n,jobn,d,z,iz,wk,ier)
c
c   arguments    a      - input real symmetric matrix of order n,
c                           whose eigenvalues and eigenvectors
c                           are to be computed. input a is
c                           destroyed if ijob is equal to 0 or 1.
c                n      - input order of the matrix a.
c                jobn   - input option parameter.  if jobn.ge.10
c                         a is assumed to be in full storage mode
c                         (in this case, a must be dimensioned exactly
c                         n by n in the calling program).
c                         if jobn.lt.10 then a is assumed to be in
c                         symmetric storage mode.  define
c                         ijob=mod(jobn,10).  then when
c                           ijob = 0, compute eigenvalues only
c                           ijob = 1, compute eigenvalues and eigen-
c                             vectors.
c                           ijob = 2, compute eigenvalues, eigenvectors
c                             and performance index.
c                           ijob = 3, compute performance index only.
c                           if the performance index is computed, it is
c                           returned in wk(1). the routines have
c                           performed (well, satisfactorily, poorly) if
c                           wk(1) is (less than 1, between 1 and 100,
c                           greater than 100).
c                d      - output vector of length n,
c                           containing the eigenvalues of a.
c                z      - output n by n matrix containing
c                           the eigenvectors of a.
c                           the eigenvector in column j of z corres-
c                           ponds to the eigenvalue d(j).
c                           if ijob = 0, z is not used.
c                iz     - input row dimension of matrix z exactly as
c                           specified in the dimension statement in the
c                           calling program.
c                wk     - work area, the length of wk depends
c                           on the value of ijob, when
c                           ijob = 0, the length of wk is at least n.
c                           ijob = 1, the length of wk is at least n.
c                           ijob = 2, the length of wk is at least
c                             n(n+1)/2+n.
c                           ijob = 3, the length of wk is at least 1.
c                ier    - error parameter (output)
c                         terminal error
c                           ier = 128+j, indicates that eqrt2s failed
c                             to converge on eigenvalue j. eigenvalues
c                             and eigenvectors 1,...,j-1 have been
c                             computed correctly, but the eigenvalues
c                             are unordered. the performance index
c                             is set to 1000.0
c                         warning error (with fix)
c                           in the following, ijob = mod(jobn,10).
c                           ier = 66, indicates ijob is less than 0 or
c                             ijob is greater than 3. ijob set to 1.
c                           ier = 67, indicates ijob is not equal to
c                             zero, and iz is less than the order of
c                             matrix a. ijob is set to zero.
c
c   reqd. routines - ehobks,ehouss,eqrt2s,uertst,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c
c-----------------------------------------------------------------------
c
c                                  specifications for arguments
      integer            n,jobn,iz,ier
      REAL  a(1),d(1),wk(1),z(iz,1)
c                                  specifications for local variables
      integer            ijob,ir,jr,ij,ji,np1
      integer            jer,na,nd,iiz,ibeg,il,kk,lk,i,j,k,l
      REAL  anorm,asum,pi,sumz,sumr,an,s,ten,rdelp,zero,
     1                   one,thous
      REAL smach
c     data           rdelp/z3410000000000000/
      data           zero,one/0.0d0,1.0d0/,ten/10.0d0/,thous/1000.0d0/
c                                  initialize error parameters
c                                  first executable statement
      ier = 0
      rdelp = smach(1)
      jer = 0
      if (jobn.lt.10) go to 15
c                                  convert to symmetric storage mode
      k = 1
      ji = n-1
      ij = 1
      do 10 j=1,n
         do 5 i=1,j
            a(k) = a(ij)
            ij = ij+1
            k = k+1
    5    continue
         ij = ij + ji
         ji = ji - 1
   10 continue
   15 ijob = mod(jobn,10)
      if (ijob.ge.0.and.ijob.le.3) go to 20
c                                  warning error - ijob is not in the
c                                    range
      ier = 66
      ijob = 1
      go to 25
   20 if (ijob.eq.0) go to 35
   25 if (iz.ge.n) go to 30
c                                  warning error - iz is less than n
c                                    eigenvectors can not be computed,
c                                    ijob set to zero
      ier = 67
      ijob = 0
   30 if (ijob.eq.3) go to 75
   35 na = (n*(n+1))/2
      if (ijob.ne.2) go to 45
      do 40 i=1,na
         wk(i) = a(i)
   40 continue
c                                  save input a if ijob = 2
   45 nd = 1
      if (ijob.eq.2) nd = na+1
c                                  reduce a to symmetric tridiagonal
c                                    form
      call ehouss (a,n,d,wk(nd),wk(nd))
      iiz = 1
      if (ijob.eq.0) go to 60
      iiz = iz
c                                  set z to the identity matrix
      do 55 i=1,n
         do 50 j=1,n
            z(i,j) = zero
   50    continue
         z(i,i) = one
   55 continue
c                                  compute eigenvalues and eigenvectors
   60 call eqrt2s (d,wk(nd),n,z,iiz,jer)
      if (ijob.eq.0) go to 9000
      if (jer.gt.128) go to 65
c                                  back transform eigenvectors
      call ehobks (a,n,1,n,z,iz)
   65 if (ijob.le.1) go to 9000
c                                  move input matrix back to a
      do 70 i=1,na
         a(i) = wk(i)
   70 continue
      wk(1) = thous
      if (jer.ne.0) go to 9000
c                                  compute 1 - norm of a
   75 anorm = zero
      ibeg = 1
      do 85 i=1,n
         asum = zero
         il = ibeg
         kk = 1
         do 80 l=1,n
            asum = asum+ dabs(a(il))
            if (l.ge.i) kk = l
            il = il+kk
   80    continue
         anorm = dmax1(anorm,asum)
         ibeg = ibeg+i
   85 continue
      if (anorm.eq.zero) anorm = one
c                                  compute performance index
      pi = zero
      do 100 i=1,n
         ibeg = 1
         s = zero
         sumz = zero
         do 95 l=1,n
            lk = ibeg
            kk = 1
            sumz = sumz+ dabs(z(l,i))
            sumr = -d(i)*z(l,i)
            do 90 k=1,n
               sumr = sumr+a(lk)*z(k,i)
               if (k.ge.l) kk = k
               lk = lk+kk
   90       continue
            s = s+ dabs(sumr)
            ibeg = ibeg+l
   95    continue
         if (sumz.eq.zero) go to 100
         pi = dmax1(pi,s/sumz)
  100 continue
      an = n
      pi = pi/(anorm*ten*an*rdelp)
      wk(1) = pi
      if (jobn.lt.10) go to 9000
c                                  convert back to full storage mode
      np1 = n+1
      ij = (n-1)*np1 + 2
      k = (n*(np1))/2
      do 110 jr=1,n
         j = np1-jr
         do 105 ir=1,j
            ij = ij-1
            a(ij) = a(k)
            k = k-1
  105    continue
         ij = ij-jr
  110 continue
      ji = 0
      k = n-1
      do 120 i=1,n
         ij = i-n
         do 115 j=1,i
            ij = ij+n
            ji = ji+1
            a(ij) = a(ji)
  115    continue
         ji = ji + k
         k = k-1
  120 continue
 9000 continue
      if (ier.ne.0) call uertst (ier,'eigrsf')
      if (jer.eq.0) go to 9005
      ier = jer
      call uertst (ier,'eigrsf')
 9005 return
      end
      subroutine ehobks (a,n,m1,m2,z,iz)
c
c-----------------------------------------------------------------------
c
c   purpose             - back transformation to form the eigenvectors
c                           of the original symmetric matrix from the
c                           eigenvectors of the tridiagonal matrix
c
c   usage               - call ehobks (a,n,m1,m2,z,iz)
c
c   arguments    a      - the array contains the details of the house-
c                           holder reduction of the original matrix a
c                           as generated by routine ehouss.(input)
c                n      - order of the real symmetric matrix.(input)
c                m1     - m1 and m2 are two input scalars such that
c                           eigenvectors m1 to m2 of the tridiagonal
c                           matrix a have been found and normalized
c                           according to the euclidean norm.
c                m2     - see above - m1
c                z      - a two dimensional array of size n x (m2-m1+1)
c                           which contains eigenvectors m1 to m2 of
c                           tridiagonal matrix t, normalized according
c                           to euclidean norm. input z can be produced
c                           by routine eqrt2s, the resultant
c                           matrix overwrites the input z.(input/output)
c                iz     - row dimension of matrix z exactly as
c                           specified in the dimension statement in the
c                           calling program. (input)
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c-----------------------------------------------------------------------
c
c
      dimension          a(*),z(iz,*)
      REAL  a,z,h,s
c                                  first executable statement
      if (n .eq. 1) go to 30
      do 25 i=2,n
         l = i-1
         ia = (i*l)/2
         h = a(ia+i)
         if (h.eq.0.d0) go to 25
c                                  derives eigenvectors m1 to m2 of
c                                  the original matrix from eigenvectors
c                                  m1 to m2 of the symmetric
c                                  tridiagonal matrix
         do 20 j = m1,m2
            s = 0.0d0
            do 10 k = 1,l
               s = s+a(ia+k)*z(k,j)
   10       continue
            s = s/h
            do 15 k=1,l
               z(k,j) = z(k,j)-s*a(ia+k)
   15       continue
   20    continue
   25 continue
   30 return
      end
      subroutine ehouss (a,n,d,e,e2)
c
c-----------------------------------------------------------------------
c
c   purpose             - reduction of a symmetric matrix to symmetric
c                           tridiagonal form using a householder
c                           reduction
c
c   usage               - call ehouss (a,n,d,e,e2)
c
c   arguments    a      - the given n x n, real symmetric matrix a,
c                           where a is stored in symmetric storage mode.
c                           the input a is replaced by the details of
c                           the householder reduction of a.
c                n      - input order of a and the length of d, e, and
c                           e2.
c                d      - the output array of length n, giving the
c                           diagonal elements of the tridiagonal matrix.
c                e      - the output array of length n, giving the sub-
c                           diagonal in the last (n-1) elements, e(1) is
c                           set to zero.
c                e2     - output array of length n.  e2(i) = e(i)**2.
c
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c-----------------------------------------------------------------------
c
c
      dimension          a(*),d(n),e(n),e2(n)
      REAL  a,d,e,e2,zero,h,scale,one,scale1,f,g,hh
      data               zero/0.0d0/,one/1.0d0/
c                                  first executable statement
      np1 = n+1
      nn = (n*np1)/2-1
      nbeg = nn+1-n
      do 70 ii = 1,n
         i = np1-ii
         l = i-1
         h = zero
         scale = zero
         if (l .lt. 1) go to 10
c                                  scale row (algol tol then not needed)
         nk = nn
         do 5 k = 1,l
            scale = scale+ dabs(a(nk))
            nk = nk-1
    5    continue
         if (scale .ne. zero) go to 15
   10    e(i) = zero
         e2(i) = zero
         go to 65
   15    nk = nn
         scale1 = one/scale
         do 20 k = 1,l
            a(nk) = a(nk)*scale1
            h = h+a(nk)*a(nk)
            nk = nk-1
   20    continue
         e2(i) = scale*scale*h
         f = a(nn)
         g = - dsign( dsqrt(h),f)
         e(i) = scale*g
         h = h-f*g
         a(nn) = f-g
         if (l .eq. 1) go to 55
         f = zero
         jk1 = 1
         do 40 j = 1,l
            g = zero
            ik = nbeg+1
            jk = jk1
c                                  form element of a*u
            do 25 k = 1,j
               g = g+a(jk)*a(ik)
               jk = jk+1
               ik = ik+1
   25       continue
            jp1 = j+1
            if (l .lt. jp1) go to 35
            jk = jk+j-1
            do 30 k = jp1,l
               g = g+a(jk)*a(ik)
               jk = jk+k
               ik = ik+1
   30       continue
c                                  form element of p
   35       e(j) = g/h
            f = f+e(j)*a(nbeg+j)
            jk1 = jk1+j
   40    continue
         hh = f/(h+h)
c                                  form reduced a
         jk = 1
         do 50 j = 1,l
            f = a(nbeg+j)
            g = e(j)-hh*f
            e(j) = g
            do 45 k = 1,j
               a(jk) = a(jk)-f*e(k)-g*a(nbeg+k)
               jk = jk+1
   45       continue
   50    continue
   55    do 60 k = 1,l
            a(nbeg+k) = scale*a(nbeg+k)
   60    continue
   65    d(i) = a(nbeg+i)
         a(nbeg+i) = h*scale*scale
         nbeg = nbeg-i+1
         nn = nn-i
   70 continue
      return
      end
      subroutine eqrt2s (d,e,n,z,iz,ier)
c
c-----------------------------------------------------------------------
c
c   purpose             - eigenvalues and (optionally) eigenvectors of
c                           a symmetric tridiagonal matrix using the
c                           ql method.
c
c   usage               - call eqrt2s (d,e,n,z,iz,ier)
c
c   arguments    d      - on input, the vector d of length n contains
c                           the diagonal elements of the symmetric
c                           tridiagonal matrix t.
c                           on output, d contains the eigenvalues of
c                           t in ascending order.
c                e      - on input, the vector e of length n contains
c                           the sub-diagonal elements of t in position
c                           2,...,n. on output, e is destroyed.
c                n      - order of tridiagonal matrix t.(input)
c                z      - on input, z contains the identity matrix of
c                           order n.
c                           on output, z contains the eigenvectors
c                           of t. the eigenvector in column j of z
c                           corresponds to the eigenvalue d(j).
c                iz     - input row dimension of matrix z exactly as
c                           specified in the dimension statement in the
c                           calling program. if iz is less than n, the
c                           eigenvectors are not computed. in this case
c                           z is not used.
c                ier    - error parameter
c                         terminal error
c                           ier = 128+j, indicates that eqrt2s failed
c                             to converge on eigenvalue j. eigenvalues
c                             and eigenvectors 1,...,j-1 have been
c                             computed correctly, but the eigenvalues
c                             are unordered.
c
c   reqd. routines - uertst,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c   remarks      routine eqrt2s is designed to accept output
c                  vectors d and e from routine ehouss as input
c                  d and e of eqrt2s. given a symmetric tridiagonal
c                  matrix, t, vector d contains the diagonal elements
c                  of t and vector e contains the subdiagonal elements
c                  of t. see the ehouss document.
c
c-----------------------------------------------------------------------
c
c
      dimension          d(*),e(*),z(iz,*)
      REAL  d,e,z,b,c,f,g,h,p,r,s,rdelp,one,zero
      REAL smach
c     data               rdelp/z3410000000000000/
      data               zero,one/0.0d0,1.0d0/
c                                  move the last n-1 elements
c                                  of e into the first n-1 locations
c                                  first executable statement
      ier  = 0
      rdelp = smach(1)
      if (n .eq. 1) go to 9005
      do 5  i=2,n
         e(i-1) = e(i)
    5 continue
      e(n) = zero
      b = zero
      f = zero
      do  60  l=1,n
         j = 0
         h = rdelp*( dabs(d(l))+ dabs(e(l)))
         if (b.lt.h) b = h
c                                  look for small sub-diagonal element
         do 10  m=l,n
            k=m
            if ( dabs(e(k)) .le. b) go to 15
   10    continue
   15    m = k
         if (m.eq.l) go to 55
   20    if (j .eq. 30) go to 85
         j = j+1
         l1 = l+1
         g = d(l)
         p = (d(l1)-g)/(e(l)+e(l))
         r =  dsqrt(p*p+one)
         d(l) = e(l)/(p+ dsign(r,p))
         h = g-d(l)
         do 25 i = l1,n
            d(i) = d(i)-h
   25    continue
         f = f+h
c                                  ql transformation
         p = d(m)
         c = one
         s = zero
         mm1 = m-1
         mm1pl = mm1+l
         if (l.gt.mm1) go to 50
         do 45 ii=l,mm1
            i = mm1pl-ii
            g = c*e(i)
            h = c*p
            if ( dabs(p).lt. dabs(e(i))) go to 30
            c = e(i)/p
            r =  dsqrt(c*c+one)
            e(i+1) = s*p*r
            s = c/r
            c = one/r
            go to 35
   30       c = p/e(i)
            r =  dsqrt(c*c+one)
            e(i+1) = s*e(i)*r
            s = one/r
            c = c*s
   35       p = c*d(i)-s*g
            d(i+1) = h+s*(c*g+s*d(i))
            if (iz .lt. n) go to 45
c                                  form vector
            do 40 k=1,n
               h = z(k,i+1)
               z(k,i+1) = s*z(k,i)+c*h
               z(k,i) = c*z(k,i)-s*h
   40       continue
   45    continue
   50    e(l) = s*p
         d(l) = c*p
         if ( dabs(e(l)) .gt.b) go to 20
   55    d(l) = d(l) + f
   60 continue
c                                  order eigenvalues and eigenvectors
      do  80  i=1,n
         k = i
         p = d(i)
         ip1 = i+1
         if (ip1.gt.n) go to 70
         do 65  j=ip1,n
            if (d(j) .ge. p) go to 65
            k = j
            p = d(j)
   65    continue
   70    if (k.eq.i) go to 80
         d(k) = d(i)
         d(i) = p
         if (iz .lt. n) go to 80
         do 75 j = 1,n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
   75    continue
   80 continue
      go to 9005
   85 ier = 128+l
c9000 continue
      call uertst(ier,'eqrt2s')
 9005 return
      end
      subroutine vcvtfs (a,n,ia,b)
c
c-----------------------------------------------------------------------
c
c   purpose             - storage mode conversion of matrices (full
c                           to symmetric)
c
c   usage               - call vcvtfs (a,n,ia,b)
c
c   arguments    a      - input matrix of dimension n by n. a contains
c                           a symmetric matrix stored in full mode.
c                n      - order of matrix a. (input)
c                ia     - row dimension of matrix a exactly as
c                           specified in the dimension statement in the
c                           calling program. (input)
c                b      - output vector of dimension n*(n+1)/2
c                           containing matrix a in symmetric storage
c                           mode.
c
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c-----------------------------------------------------------------------
c
c
      REAL  a(ia,1),b(1)
c                                  first executable statement
      k = 1
      do 10 i = 1,n
         do 5 j = 1,i
            b(k) = a(j,i)
            k = k+1
    5    continue
   10 continue
      return
      end
      subroutine vcvtsf (a,n,b,ib)
c
c-----------------------------------------------------------------------
c
c   purpose             - storage mode conversion of matrices
c                           (symmetric to full)
c
c   usage               - call vcvtsf (a,n,b,ib)
c
c   arguments    a      - input vector of length n*(n+1)/2 containing
c                           an n by n symmetric matrix stored in
c                           symmetric storage mode.
c                n      - order of matrix a. (input)
c                b      - output matrix of dimension n by n containing
c                           matrix a in full storage mode.
c                ib     - row dimension of matrix b exactly as
c                           specified in the dimension statement in the
c                           calling program. (input)
c
c
c-----------------------------------------------------------------------
c
c
      REAL  a(*),b(ib,*)
c                                  first executable statement
      i1 = (n*(n+1))/2
      j = n
    5 jp1 = j+1
      do 10 k = 1,j
         i = jp1-k
         b(i,j) = a(i1)
         i1 = i1-1
   10 continue
      j = j-1
      if (j .ge. 1) go to 5
      if (n.lt.2) go to 25
      do 20 i = 2,n
         i1 = i-1
         do 15 j = 1,i1
            b(i,j) = b(j,i)
   15    continue
   20 continue
   25 return
      end
      subroutine eipol( pol,ipol,ideg,amat,bmat,cmat,lmat,rvec,deriv,
     1 const,eigen,diag,work,degen,maxn,maxm,ndim,nd,ndim2,mdim )
c
c***********************************************************************
c
c     the subroutine eipol calculates the eigenvalues of the matrix l
c     at the diagonal elements of c.
c
c     input
c     ipol                number of the pole of lmat
c     pol                 cmat(ipol)
c     ideg                degeneracy of the pole pol
c
c     output
c     diag(i)             eigenvalues ordered with respect to increasing
c                         magnitude. the diverging eigenvalues are put
c                         equal to -d10 * i, where i runs from 1 to ideg
c     eigen(i)            numerator of the i-th diverging eigenvalue
c     const(i)            additive constant for this eigenvalue
c     ideg                number of linearly independent vectors of the
c                         set of coupling matrix elements.
c
c***********************************************************************
c
      implicit REAL  (a-h,o-z)
      REAL  lmat
c
      dimension amat(nd),bmat(ndim,mdim),cmat(mdim),lmat(nd),work(ndim),
     *    rvec(ndim,ndim),deriv(nd),const(ndim2),eigen(ndim2),diag(ndim)
c
      data d0,d05,d1,d2,d10,d01 /0.0d0,0.5d0,1.0d0,2.0d0,1.0d10,1.0d-10/
c
c***********************************************************************
c
c     special case for maxn=1
c
c***********************************************************************
c
      if(maxn .gt. 1) go to 60
      diag(1) = -d10
      const(2) = amat(1)
      do 20 i = 1,maxm
      if( dabs(pol-cmat(i)) .lt. degen) go to 20
      const(2) = const(2)  +  bmat(1,i) * bmat(1,i) / ( pol - cmat(i) )
 20   continue
      eigen(2) = d0
      do 40 i = 1,ideg
      eigen(2) = eigen(2) + bmat(1,ipol+i-1) ** 2
 40   continue
      ideg = 1
      return
c
c***********************************************************************
c
c     calculation of the diverging eigenvalues and eigenvectors for
c     maxn.gt.1 and ideg.eq.1
c
c***********************************************************************
c
   60 if(ideg .gt. 1) go to 340
c
      const(maxn+1) = d0
      eigen(maxn+1) = d0
      do 80 i = 1,maxn
      eigen(maxn+1) = eigen(maxn+1) + bmat(i,ipol) ** 2
   80 continue
      dnorm =  dsqrt( eigen(maxn+1) )
c
c     calaculation of the eigenvectors by schmidt orthogonalization
c
      do 100 i = 1,maxn
      rvec(i,1) = bmat(i,ipol) / dnorm
  100 continue
      do 140 j = 2,maxn
      do 120 i = 1,maxn
      rvec(i,j) = d0
  120 continue
      rvec(j,j) = d1
  140 continue
      do 220 j = 2,maxn
      imax = j - 1
      do 160 i = 1,imax
      xfak = rvec(j,i)
      do 160 k = 1,maxn
      rvec(k,j) = rvec(k,j) - xfak * rvec(k,i)
  160 continue
      dnorm = d0
      do 180 k = 1,maxn
      dnorm = dnorm + rvec(k,j) ** 2
  180 continue
      dnorm =  dsqrt( dnorm )
      do 200 k = 1,maxn
      rvec(k,j) = rvec(k,j) / dnorm
  200 continue
  220 continue
c
c     repetition of the orthogonalization procedure if maxn si greater
c     than 10
c
      if(maxn .le. 10) go to 400
c
      do 320 j = 2,maxn
      imax = j - 1
      do 260 i = 1,imax
      xfak = d0
      do 240 k = 1,maxn
      xfak = xfak + rvec(k,j) * rvec(k,i)
  240 continue
      do 260 k = 1,maxn
      rvec(k,j) = rvec(k,j) - xfak * rvec(k,i)
  260 continue
      dnorm = d0
      do 280 k = 1,maxn
      dnorm = dnorm + rvec(k,j) ** 2
  280 continue
      dnorm =  dsqrt( dnorm )
      do 300 k = 1,maxn
      rvec(k,j) = rvec(k,j) / dnorm
  300 continue
  320 continue
c
      go to 400
c
c***********************************************************************
c
c     calculation of the diverging eigenvalue and eigenvectors for
c     ideg .gt. 1
c     note: the negative of the matrix is diagonlized in order that
c     the eigenvalues which are not zero are found at the beginning of
c     the array diag
c
c***********************************************************************
c
  340 nz = 0
      do 360 i = 1,maxn
      do 360 j = 1,i
      nz = nz + 1
      lmat(nz) = d0
      do 360 k = 1,ideg
      lmat(nz) = lmat(nz) - bmat(i,ipol+k-1) * bmat(j,ipol+k-1)
  360 continue
      call eigrsf(lmat,maxn,1,diag,rvec,ndim,work,jer)
      imax = ideg
      if(imax .gt. maxn) imax = maxn
      ideg = 0
      do 380 i = 1,imax
      if( dabs(diag(i)) .ge. d01) ideg = ideg + 1
      eigen(maxn+i) = -diag(i)
  380 continue
c
c***********************************************************************
c
c     calculation of the matrix l tilde
c
c***********************************************************************
c
  400 nz = 0
      do 440 i = 1,maxn
      do 440 j = 1,i
      nz = nz + 1
      deriv(nz) = d0
      do 420 k = 1,maxm
      if( dabs(pol-cmat(k)) .lt. degen) go to 420
      deriv(nz) = deriv(nz) + bmat(i,k) * bmat(j,k) / (pol - cmat(k))
  420 continue
      deriv(nz) = deriv(nz) + amat(nz)
  440 continue
c
c***********************************************************************
c
c     calculation of   rvec * l * rvec+
c
c***********************************************************************
c
      if(ideg .ge. maxn) go to 540
      nz = 0
      imin = ideg + 1
      do 520 i = imin,maxn
      do 520 j = imin,i
      nz = nz + 1
      lmat(nz) = d0
      lz = 0
      do 520 k = 1,maxn
      do 500 l = 1,k
      lz = lz + 1
      add = deriv(lz) * (rvec(k,i)*rvec(l,j)+rvec(l,i)*rvec(k,j))
      lmat(nz) = lmat(nz)  +  add
  500 continue
      lmat(nz) = lmat(nz)  -  add * d05
  520 continue
c
c***********************************************************************
c
c     calculation of the additive constant for the diverging eigenvalues
c
c***********************************************************************
c
  540 do 580 i = 1,ideg
      const(maxn+i) = d0
      nz = 0
      do 580 k = 1,maxn
      do 560 l = 1,k
      nz = nz + 1
      add = deriv(nz) * (rvec(k,i) * rvec(l,i)) * d2
      const(maxn+i) = const(maxn+i)  +  add
  560 continue
      const(maxn+i) = const(maxn+i)  -  add * d05
  580 continue
c
c***********************************************************************
c
c     diagonalization to obtain the non diverging eigenvalues
c
c***********************************************************************
c
      if(ideg .eq. maxn) go to 640
c
      jmax = maxn - ideg
      call eigrsf(lmat,jmax,0,diag,rvec,ndim,work,jer)
c
      do 600 j = 1,jmax
      work(j) = diag(j)
  600 continue
      do 620 j = 1,jmax
      diag(ideg+j) = work(j)
  620 continue
  640 do 660 i = 1,ideg
      iz = imax - i + 1
      diag(i) = -  dfloat(iz) * d10
  660 continue
c
      return
      end
      subroutine eival(omega,ieig,ipol,amat,bmat,cmat,lmat,rvec,deriv,
     1 diag,work,maxn,maxm,ndim,nd,mdim )
c
c***********************************************************************
c
c     this subroutine calculates the eigenvalues of the matrix l for
c     arbitrary values of omega.
c
c     input
c     omega                energy value for which the matrix lmat is
c                          diagonalized
c     ipol                 number of the eigenvalue whose derivative
c                          is to be calculated
c     ieig = 0             without calculation of the derivative
c     ieig = 1             with calculation of the derivative
c
c     output
c     diag(i)             i-th eigenvalue of lmat
c     work(ipol)          derivative of the eigenvalue ipol
c
c***********************************************************************
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
      REAL  lmat
c
      dimension amat(nd),bmat(ndim,mdim),cmat(mdim),lmat(nd),
     * rvec(ndim,ndim),deriv(nd),diag(ndim),work(ndim)
c
      data d0,d2 /0.0d0,2.0d0/
c
    1 format (/1x,
     +'omega value coincides with pole,calculation must be stopped')
c
c***********************************************************************
c
c     calculation of the matrix lmat and of its derivatve
c
c***********************************************************************
c
      nz = 0
      do 20 i = 1,maxn
      do 20 j = 1,i
      nz = nz + 1
      lmat(nz) = d0
      deriv(nz) = d0
      do 10 k = 1,maxm
      if ( dabs(omega-cmat(k)).eq.0.0d0) go to 100
      lmat(nz) = lmat(nz)  +  bmat(i,k) * bmat(j,k) / (omega - cmat(k))
      if ((omega-cmat(k))**2.le.1.0d-30) go to 5
      deriv(nz) = deriv(nz)  -  bmat(i,k) * bmat(j,k) /
     *                          ((omega - cmat(k)) ** 2)
      go to 10
    5 deriv(nz)=deriv(nz)-1.0d+60
   10 continue
      lmat(nz) = lmat(nz) + amat(nz)
   20 continue
c
      call eigrsf(lmat,maxn,ieig,diag,rvec,ndim,work,jer)
c
      if(ieig .eq. 0) return
c
c***********************************************************************
c
c     calculation of the derivative
c
c***********************************************************************
c
      work(ipol) = d0
      nz = 0
      do 40 j = 1,maxn
      do 30 k = 1,j
      nz = nz + 1
      add = rvec(j,ipol) * rvec(k,ipol) * deriv(nz)
      work(ipol) = work(ipol)  +  d2 * add
   30 continue
      work(ipol) = work(ipol) - add
   40 continue
      go to 200
  100 continue
      write (iwr,1)
_IF(rs6000,hp700,hp800,linux)
c
c use C exit routine
c
      call exitc(-1)
_ELSE
      call exit(-1)
_ENDIF
  200 continue
c
      return
      end
      subroutine eivec( ipol,bmat,cmat,lmat,rvec,work,eigen,deriv,maxn,
     1 maxm,ndim,nd,mdim,ivec,ifile,neigen )
c
c***********************************************************************
c
c     tis subroutine calculates the eigenvector belonging to the given
c     eigenvalue eigen with the derivative deriv. for ivec=0 only the
c     eigenvalue is written on the file ifile. for ivec=1 also the
c     eigenvector is calculated and is written on the file ifile
c     (blocked with nd).
c
c     input
c     ipol                 number of the eigenvalue of lmat
c     eigen                eigenvalue ipol of lmat
c     deriv                derivative of eigen
c     rvec                 matrix of eigenvectors of lmat
c
c     output
c     eigenvalues and eigenvectors of the bordered matrix on file ifile
c
c***********************************************************************
c
      implicit REAL  (a-h,o-z)
      REAL  lmat
c
      dimension bmat(ndim,mdim),cmat(mdim),rvec(ndim,ndim),lmat(nd),
     1 work(ndim)
c
      data d0,d1 /0.0d0,1.0d0/
c
      neigen = neigen + 1
      write(ifile) eigen,eigen,eigen
      if(ivec .eq. 0) return
c
c***********************************************************************
c
c     calculation of the eigenvector components with respect to the
c     submatrix a
c
c***********************************************************************
c
      do 10 i = 1,maxn
      work(i) = rvec(i,ipol) /  dsqrt(d1 - deriv)
      lmat(i) = work(i)
   10 continue
c
c***********************************************************************
c
c     calculation of the remaining eigenvector components
c
c***********************************************************************
c
      iz = maxn
      do 50 i = 1,maxm
      iz = iz + 1
      if(iz .eq. nd+1) iz = 1
      lmat(iz) =  d0
      do 20 j = 1,maxn
      lmat(iz) = lmat(iz) + bmat(j,i) * work(j)
   20 continue
      if ( dabs(eigen-cmat(i)).le.1.0d-12) go to 30
      lmat(iz) = lmat(iz) / (eigen - cmat(i))
      go to 40
   30 lmat(iz)=d1
   40 if(iz .eq. nd) write(ifile) lmat
   50 continue
      if(iz .ne. nd) write(ifile) lmat
c
      return
      end
      subroutine inter( ipol,pol1,pol2,eigen1,eigen2,x2,amat,bmat,cmat,
     * lmat,rvec,deriv,diag,work,crit,maxn,maxm,ndim,nd,ndim2,mdim,
     * maxpsa,iprint,iauf,ieradd )
c
c***********************************************************************
c
c     this subroutine calculates the points of intersection of the
c     eigenvalue ipol of the matrix lmat with the straight line y=omega
c     via the method of interval bisecting. this method is used in case
c     the newton-raphson procedure fails.
c
c     this subroutine is called in two different ways:
c     iauf = 1     the interval bisecting is stopped, when the values y1
c                  and y2 lie in the interval pol1,pol2. return is to
c                  the newton-raphson procedure in the calling program.
c
c     iauf = 0     the interval bisecting is continued until convergence
c                  is reached. no return is made to the newton-raphson
c                  procedure in the calling program.
c
c     input
c     ipol                 number of the eigenvalue
c     pol1                 left hand margin of the interval
c     pol2                 right hand margin of the interval
c     eigen1               eigenvalue at left hand margin of interval
c     eigen2               eigenvalue at right hand margin of interval
c
c     output
c     x2                   pole of the green's function
c
c***********************************************************************
c
      implicit REAL  (a-h,o-z)
      REAL  lmat
c
      dimension amat(nd),bmat(ndim,mdim),cmat(mdim),lmat(nd),
     *          rvec(ndim,ndim),deriv(nd),diag(ndim),work(ndim)
c
    1 format(//' ****   after ',i4,
     +' iterations it was not possible to determine the point of'/
     +' intersection in the given interval.',7x,'****'
     *//' ****   interval - limits:',47x,f15.4,',',f15.4,14x,'****'
     *//' ****   eigenvalues of the matrix lmat at the interval-limits:'
     *,11x,f15.4,',',f15.4,14x,'****'
     *//' ****   last iteration steps:',44x,f15.4,',',f15.4,14x,
     * '****'
     *//' ****   corresponding eigenvalues of lmat:',31x,f15.4,',',
     * f15.4,14x,'****')
c
    2 format(//' ****   for the ',i4, 'th eigenvalue of the'/
     +1x,'matrix lmat the value determined by iteration x = ',f10.4,
     +15x,'****'/
     +' **** coincides with a diagonal element of the matrix cmat.',
     +57x,'****'
     *//' ****   interval - limits:',47x,f15.4,',',f15.4,14x,'****'
     *//' ****   eigenvalues of the matrix lmat at the interval-limits:'
     *,11x,f15.4,',',f15.4,14x,'****')
c
    3 format (//1x,
     +'****   two poles are closer together than 1.0e-04 a.u.'/
     +1x,'no interval bisecting is therefore performed ****'
     +/' ***   this occurs for the poles :',2f12.7,43x,'****')
c
      data d0,d05 /0.0d0,0.5d0/
c
      imax = 10 * maxpsa
      iauf = iauf - 1
c
c***********************************************************************
c
c     determine the original interval
c
c***********************************************************************
c
      x1 = pol1
      x2 = pol2
      if ( dabs(x1-x2).lt.1.0d-04) go to 140
      delta = (x1 - x2) * d05
      y1 = eigen1
      y2 = eigen2
c
c***********************************************************************
c
c     interval bisecting
c
c***********************************************************************
c
      do 10 i = 1,imax
      diff1 = x1 - y1
      x2 = x1 - delta
      if(x2 .eq. pol1  .or.  x2 .eq. pol2) go to 120
c
      call eival(x2,0,ipol,amat,bmat,cmat,lmat,rvec,deriv,diag,work,
     * maxn,maxm,ndim,nd,mdim )
c
      y2 = diag(ipol)
      diff2 = x2 - y2
      if( dabs(diff2) .lt. crit) go to 20
      delta = delta * d05
c
c***********************************************************************
c
c     determine whether the values lie in the interval pol1,pol2
c
c***********************************************************************
c
      if(y1 .gt. pol1  .and.  y2 .gt. pol1  .and.  y1 .lt. pol2  .and.
     * y2 .lt. pol2  .and.  iauf .eq. 0) go to 20
c
c***********************************************************************
c
c     reset the margins of the interval
c
c***********************************************************************
c
      if(diff1*diff2 .lt. d0) go to 10
      x1 = x2
      y1 = y2
   10 continue
c
      goto 100
   20 call eival(x2,1,ipol,amat,bmat,cmat,lmat,rvec,deriv,diag,work,
     * maxn,maxm,ndim,nd,mdim )
      return
c
c***********************************************************************
c
c     error information output
c
c***********************************************************************
c
  100 ieradd = 1
      write(iprint,1) imax,pol1,pol2,eigen1,eigen2,x1,x2,y1,y2
      return
c
  120 ieradd = 1
      write(iprint,2) ipol,x2,pol1,pol2,eigen1,eigen2
      return
c
  140 ieradd = 1
      write (iprint,3) pol1,pol2
c
      return
      end
      subroutine leqt2f (a,m,n,ia,b,idgt,wkarea,ier)
c
c-----------------------------------------------------------------------
c
c   purpose             - linear equation solution - full storage
c                           mode - high accuracy solution
c
c   usage               - call leqt2f (a,m,n,ia,b,idgt,wkarea,ier)
c
c   arguments    a      - input matrix of dimension n by n containing
c                           the coefficient matrix of the equation
c                           ax = b.
c                m      - number of right-hand sides. (input)
c                n      - order of a and number of rows in b. (input)
c                ia     - row dimension of a and b exactly as specified
c                           in the dimension statement in the calling
c                           program. (input)
c                b      - input matrix of dimension n by m containing
c                           the right-hand sides of the equation ax = b.
c                         on output, the n by m matrix of solutions
c                           replaces b.
c                idgt   - input option.
c                         if idgt is greater than 0, the elements of
c                           a and b are assumed to be correct to idgt
c                           decimal digits and the routine performs
c                           an accuracy test.
c                         if idgt equals 0, the accuracy test is
c                           bypassed.
c                         on output, idgt contains the approximate
c                           number of digits in the answer which
c                           were unchanged after improvement.
c                wkarea - work area of dimension greater than or equal
c                           to n**2+3n.
c                ier    - error parameter. (output)
c                         warning error
c                           ier = 34 indicates that the accuracy test
c                             failed. the computed solution may be in
c                             error by more than can be accounted for
c                             by the uncertainty of the data. this
c                             warning can be produced only if idgt is
c                             greater than 0 on input. (see the
c                             chapter l prelude for further discussion.)
c                         terminal error
c                           ier = 129 indicates that the matrix is
c                             algorithmically singular. (see the
c                             chapter l prelude).
c                           ier = 131 indicates that the matrix is too
c                             ill-conditioned for iterative improvement
c                             to be effective.
c
c   reqd. routines - single/ludatf,luelmf,lureff,uertst,ugetio
c                  - double/ludatf,luelmf,lureff,uertst,ugetio,
c                      vxadd,vxmul,vxsto
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c-----------------------------------------------------------------------
c
c
      dimension          a(ia,*),b(ia,*),wkarea(*)
      REAL  a,b,wkarea,d1,d2,wa
c                                  first executable statement
c                                  initialize ier
      ier=0
      jer=0
      j = n*n+1
      k = j+n
      mm = k+n
      kk = 0
      mm1 = mm-1
      jj=1
      do 5 l=1,n
         do 5 i=1,n
            wkarea(jj)=a(i,l)
            jj=jj+1
    5 continue
c                                  decompose a
      call ludatf (wkarea(1),a,n,n,idgt,d1,d2,wkarea(j),wkarea(k),
     1   wa,ier)
      if (ier.gt.128) go to 25
      if (idgt .eq. 0 .or. ier .ne. 0) kk = 1
      do 15 i = 1,m
c                                  performs the elimination part of
c                                  ax = b
         call luelmf (a,b(1,i),wkarea(j),n,n,wkarea(mm))
c                                  refinement of solution to ax = b
         if (kk .ne. 0)
     *   call lureff (wkarea(1),b(1,i),a,wkarea(j),n,n,wkarea(mm),idgt,
     *               wkarea(k),wkarea(k),jer)
         do 10 ii=1,n
            b(ii,i) = wkarea(mm1+ii)
   10    continue
         if (jer.ne.0) go to 20
   15 continue
      go to 25
   20 ier = 131
   25 jj=1
      do 30 j = 1,n
         do 30 i = 1,n
            a(i,j)=wkarea(jj)
            jj=jj+1
   30 continue
      if (ier .eq. 0) go to 9005
c9000 continue
      call uertst (ier,'leqt2f')
 9005 return
      end
      subroutine ludatf (a,lu,n,ia,idgt,d1,d2,ipvt,equil,wa,ier)
c
c-----------------------------------------------------------------------
c
c   purpose             - l-u decomposition by the crout algorithm
c                           with optional accuracy test.
c
c   usage               - call ludatf (a,lu,n,ia,idgt,d1,d2,ipvt,
c                           equil,wa,ier)
c
c   arguments    a      - input matrix of dimension n by n containing
c                           the matrix to be decomposed.
c                lu     - real output matrix of dimension n by n
c                           containing the l-u decomposition of a
c                           rowwise permutation of the input matrix.
c                           for a description of the format of lu, see
c                           example.
c                n      - input scalar containing the order of the
c                           matrix a.
c                ia     - input scalar containing the row dimension of
c                           matrices a and lu exactly as specified in
c                           the calling program.
c                idgt   - input option.
c                           if idgt is greater than zero, the non-zero
c                           elements of a are assumed to be correct to
c                           idgt decimal places.  ludatf performs an
c                           accuracy test to determine if the computed
c                           decomposition is the exact decomposition
c                           of a matrix which differs from the given
c                           one by less than its uncertainty.
c                         if idgt is equal to zero, the accuracy test
c                           is bypassed.
c                d1     - output scalar containing one of the two
c                           components of the determinant. see
c                           description of parameter d2, below.
c                d2     - output scalar containing one of the
c                           two components of the determinant. the
c                           determinant may be evaluated as (d1)(2**d2).
c                ipvt   - output vector of length n containing the
c                           permutation indices. see document
c                           (algorithm).
c                equil  - output vector of length n containing
c                           reciprocals of the absolute values of
c                           the largest (in absolute value) element
c                           in each row.
c                wa     - accuracy test parameter, output only if
c                           idgt is greater than zero.
c                           see element documentation for details.
c                ier    - error parameter. (output)
c                         terminal error
c                           ier = 129 indicates that matrix a is
c                             algorithmically singular. (see the
c                             chapter l prelude).
c                         warning error
c                           ier = 34 indicates that the accuracy test
c                             failed.  the computed solution may be in
c                             error by more than can be accounted for
c                             by the uncertainty of the data.  this
c                             warning can be produced only if idgt is
c                             greater than 0 on input.  see chapter l
c                             prelude for further discussion.
c
c   reqd. routines - uertst,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c   remarks      a test for singularity is made at two levels:
c                1.  a row of the original matrix a is null.
c                2.  a column becomes null in the factorization process.
c
c-----------------------------------------------------------------------
c
c
      dimension          a(ia,*),lu(ia,*),ipvt(*),equil(*)
      REAL  a,lu,d1,d2,equil,wa,zero,one,four,sixtn,sixth,
     *                   rn,wrel,biga,big,p,sum,ai,wi,t,test,q
      data               zero,one,four,sixtn,sixth/0.d0,1.d0,4.d0,
     *                   16.0d0,.0625d0/
c                                  first executable statement
c                                  initialization
      ier = 0
      rn = n
      wrel = zero
      d1 = one
      d2 = zero
      biga = zero
      do 10 i=1,n
         big = zero
         do 5 j=1,n
            p = a(i,j)
            lu(i,j) = p
            p =  dabs(p)
            if (p .gt. big) big = p
    5    continue
         if (big .gt. biga) biga = big
         if (big .eq. zero) go to 110
         equil(i) = one/big
   10 continue
      do 105 j=1,n
         jm1 = j-1
         if (jm1 .lt. 1) go to 40
c                                  compute u(i,j), i=1,...,j-1
         do 35 i=1,jm1
            sum = lu(i,j)
            im1 = i-1
            if (idgt .eq. 0) go to 25
c                                  with accuracy test
            ai =  dabs(sum)
            wi = zero
            if (im1 .lt. 1) go to 20
            do 15 k=1,im1
               t = lu(i,k)*lu(k,j)
               sum = sum-t
               wi = wi+ dabs(t)
   15       continue
            lu(i,j) = sum
   20       wi = wi+ dabs(sum)
            if (ai .eq. zero) ai = biga
            test = wi/ai
            if (test .gt. wrel) wrel = test
            go to 35
c                                  without accuracy
   25       if (im1 .lt. 1) go to 35
            do 30 k=1,im1
               sum = sum-lu(i,k)*lu(k,j)
   30       continue
            lu(i,j) = sum
   35    continue
   40    p = zero
c                                  compute u(j,j) and l(i,j), i=j+1,...,
         do 70 i=j,n
            sum = lu(i,j)
            if (idgt .eq. 0) go to 55
c                                  with accuracy test
            ai =  dabs(sum)
            wi = zero
            if (jm1 .lt. 1) go to 50
            do 45 k=1,jm1
               t = lu(i,k)*lu(k,j)
               sum = sum-t
               wi = wi+ dabs(t)
   45       continue
            lu(i,j) = sum
   50       wi = wi+ dabs(sum)
            if (ai .eq. zero) ai = biga
            test = wi/ai
            if (test .gt. wrel) wrel = test
            go to 65
c                                  without accuracy test
   55       if (jm1 .lt. 1) go to 65
            do 60 k=1,jm1
               sum = sum-lu(i,k)*lu(k,j)
   60       continue
            lu(i,j) = sum
   65       q = equil(i)* dabs(sum)
            if (p .ge. q) go to 70
            p = q
            imax = i
   70    continue
c                                  test for algorithmic singularity
         if (rn+p .eq. rn) go to 110
         if (j .eq. imax) go to 80
c                                  interchange rows j and imax
         d1 = -d1
         do 75 k=1,n
            p = lu(imax,k)
            lu(imax,k) = lu(j,k)
            lu(j,k) = p
   75    continue
         equil(imax) = equil(j)
   80    ipvt(j) = imax
         d1 = d1*lu(j,j)
   85    if ( dabs(d1) .le. one) go to 90
         d1 = d1*sixth
         d2 = d2+four
         go to 85
   90    if ( dabs(d1) .ge. sixth) go to 95
         d1 = d1*sixtn
         d2 = d2-four
         go to 90
   95    continue
         jp1 = j+1
         if (jp1 .gt. n) go to 105
c                                  divide by pivot element u(j,j)
         p = lu(j,j)
         do 100 i=jp1,n
            lu(i,j) = lu(i,j)/p
  100    continue
  105 continue
c                                  perform accuracy test
      if (idgt .eq. 0) go to 9005
      p = 3*n+3
      wa = p*wrel
      if (wa+10.d0**(-idgt) .ne. wa) go to 9005
      ier = 34
      go to 9000
c                                  algorithmic singularity
  110 ier = 129
      d1 = zero
      d2 = zero
 9000 continue
c                                  print error
      call uertst(ier,'ludatf')
 9005 return
      end
      subroutine lureff (a,b,ul,ipvt,n,ia,x,idgt,res,dx,ier)
c
c-----------------------------------------------------------------------
c
c   purpose             - refinement of solution to linear equations -
c                           full storage mode
c
c   usage               - call lureff (a,b,ul,ipvt,n,ia,x,idgt,res,dx,
c                           ier)
c
c   arguments    a      - the coefficient matrix, ax=b, where a
c                           is n x n. (input)
c                b      - the right hand side, a vector of size n.
c                           (input)
c                ul     - a given n x n matrix, ul is the lu
c                           decomposition of a as supplied by 
c                           routine ludatf. (input)
c                ipvt   - a given vector of pivot indices of size n as
c                           supplied by routine ludatf. (input)
c                n      - order of a and ul, and also the length of
c                           b, ipvt, x, res, and dx. (input)
c                ia     - row dimension of a and ul exactly as
c                           specified in the dimension statement in the
c                           calling program. (input)
c                x      - an input vector of size n, x is an estimate
c                           to the solution of ax=b. on output, the
c                           improved result overwrites the input vector
c                           x.
c                idgt   - approximate number of digits in the answer
c                           which were unchanged after improvement.
c                           (output)
c                res    - the residual vector of size n, used as a work
c                           vector.
c                dx     - a work vector of size n.
c                ier    - error parameter. (output)
c                         terminal error
c                           ier=129 indicates iterative improvement
c                             failed. matrix is too ill conditioned.
c
c   reqd. routines - single/luelmf,uertst,ugetio
c                  - double/luelmf,uertst,ugetio,vxadd,vxmul,
c                      vxsto
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c-----------------------------------------------------------------------
c
c
      dimension          a(ia,*),ul(ia,*),b(*),x(*),res(*),dx(*),ipvt(*)
      dimension          accxt(2)
      REAL  a,accxt,b,ul,x,res,dx,zero,xnorm,dxnorm
      data               itmax/75/,zero/0.d0/
c                                  first executable statement
      ier=0
      xnorm = zero
      do 10 i=1,n
         xnorm = dmax1(xnorm, dabs(x(i)))
   10 continue
      if (xnorm .ne. zero) go to 20
      idgt = 50
      go to 9005
   20 do 45 iter=1,itmax
         do 30 i=1,n
      accxt(1) = 0.0d0
      accxt(2) = 0.0d0
            call vxadd(b(i),accxt)
            do 25 j=1,n
               call vxmul(-a(i,j),x(j),accxt)
   25       continue
            call vxsto(accxt,res(i))
   30    continue
         call luelmf (ul,res,ipvt,n,ia,dx)
         dxnorm = zero
         xnorm = zero
         do 35 i=1,n
            x(i) = x(i) + dx(i)
            dxnorm = dmax1(dxnorm, dabs(dx(i)))
            xnorm = dmax1(xnorm, dabs(x(i)))
   35    continue
         if (iter .ne. 1) go to 40
         idgt = 50
         if (dxnorm .ne. zero) idgt = -dlog10(dxnorm/xnorm)
   40    if (xnorm+dxnorm .eq. xnorm) go to 9005
   45 continue
c                                  iteration did not converge
      ier = 129
c9000 continue
      call uertst(ier,'lureff')
 9005 return
      end
      subroutine uertst (ier,name)
c
c-----------------------------------------------------------------------
c
c   purpose             - print a message reflecting an error condition
c
c   usage               - call uertst (ier,name)
c
c   arguments    ier    - error parameter. (input)
c                           ier = i+j where
c                             i = 128 implies terminal error,
c                             i =  64 implies warning with fix, and
c                             i =  32 implies warning.
c                             j = error code relevant to calling
c                                 routine.
c                name   - a six character literal string giving the
c                           name of the calling routine. (input)
c
c   precision/hardware  - single/all
c
c   reqd. routines - ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c   remarks      the error message produced by uertst is written
c                onto the standard output unit. the output unit
c                number can be determined by calling ugetio as
c                follows..   call ugetio(1,nin,nout).
c                the output unit number can be changed by calling
c                ugetio as follows..
c                                nin = 0
c                                nout = new output unit number
c                                call ugetio(3,nin,nout)
c                see the ugetio document for more details.
c
c-----------------------------------------------------------------------
c
c                                  specifications for arguments
      integer            ier
      character*6        name
c                                  specifications for local variables
      character*6        namset,nameq,ieq
      data               namset /'uerset'/
      data               nameq  /'      '/
      data               ieq    /'='/
c                                  first executable statement
      data               level/4/,ieqdf/0/
      if (ier.gt.999) go to 25
      if (ier.lt.-32) go to 55
      if (ier.le.128) go to 5
      if (level.lt.1) go to 30
c                                  print terminal message
      call ugetio(1,nin,iounit)
      if (ieqdf.eq.1) write(iounit,35) ier,nameq,ieq,name
      if (ieqdf.eq.0) write(iounit,35) ier,name
      go to 30
    5 if (ier.le.64) go to 10
      if (level.lt.2) go to 30
c                                  print warning with fix message
      call ugetio(1,nin,iounit)
      if (ieqdf.eq.1) write(iounit,40) ier,nameq,ieq,name
      if (ieqdf.eq.0) write(iounit,40) ier,name
      go to 30
   10 if (ier.le.32) go to 15
c                                  print warning message
      if (level.lt.3) go to 30
      call ugetio(1,nin,iounit)
      if (ieqdf.eq.1) write(iounit,45) ier,nameq,ieq,name
      if (ieqdf.eq.0) write(iounit,45) ier,name
      go to 30
   15 continue
c                                  check for uerset call
      if (name.ne.namset) go to 25
      levold = level
      level = ier
      ier = levold
      if (level.lt.0) level = 4
      if (level.gt.4) level = 4
      go to 30
   25 continue
      if (level.lt.4) go to 30
c                                  print non-defined message
      call ugetio(1,nin,iounit)
      if (ieqdf.eq.1) write(iounit,50) ier,nameq,ieq,name
      if (ieqdf.eq.0) write(iounit,50) ier,name
   30 ieqdf = 0
      return
   35 format(' *** terminal error',10x,'(ier = ',i3,
     1       ') from routine ',3a2,a1,3a2)
   40 format(' *** warning with fix error  (ier = ',i3,
     1       ') from routine ',3a2,a1,3a2)
   45 format('  *** warning error',11x,'(ier = ',i3,
     1       ') from routine ',3a2,a1,3a2)
   50 format(' *** undefined error',9x,'(ier = ',i5,
     1       ') from routine ',3a2,a1,3a2)
c                                  save p for p = r case
c                                    p is the page name
c                                    r is the routine name
   55 ieqdf = 1
      nameq = name
      return
      end
      subroutine ugetio(iopt,nin,nout)
c
c-----------------------------------------------------------------------
c
c   purpose             - to retrieve current values and to set new
c                           values for input and output unit
c                           identifiers.
c
c   usage               - call ugetio(iopt,nin,nout)
c
c   arguments    iopt   - option parameter. (input)
c                           if iopt=1, the current input and output
c                           unit identifier values are returned in nin
c                           and nout, respectively.
c                           if iopt=2 (3) the internal value of
c                           nin (nout) is reset for subsequent use.
c                nin    - input unit identifier.
c                           output if iopt=1, input if iopt=2.
c                nout   - output unit identifier.
c                           output if iopt=1, input if iopt=3.
c
c   precision/hardware  - single/all
c
c   reqd. routines - none required
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c   remarks      each routine that performs input and/or output
c                operations calls ugetio to obtain the current unit
c                identifier values. if ugetio is called with iopt=2 or 3
c                new unit identifier values are established. subsequent
c                input/output is performed on the new units.
c
c-----------------------------------------------------------------------
c
c                                  specifications for arguments
      integer            iopt,nin,nout
INCLUDE(common/iofile)
c                                  first executable statement
      if (iopt.eq.3) go to 10
      if (iopt.eq.2) go to 5
      if (iopt.ne.1) go to 9005
      nin = ird
      nout = iwr
      go to 9005
    5 ird = nin
      go to 9005
   10 iwr = nout
 9005 return
      end
      subroutine vxadd(a,acc)
c
c-----------------------------------------------------------------------
c
c   purpose             - extended precision add
c
c   usage               - call vxadd (a,acc)
c
c   arguments    a      - real              number to be added to the
c                           accumulator. (input)
c                acc    - accumulator. (input and output)
c                           acc is a real              vector of length
c                           2. on output, acc contains the sum of
c                           input acc and a.
c
c
c   reqd. routines - none required
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c   remarks      vxadd adds the real              number a to the
c                extended precision accumulator, acc. the subroutine
c                assumes that an extended precision number is already in
c                the accumulator. therefore, before the first call to
c                vxadd, acc(1) and acc(2) must be set to zero.
c
c-----------------------------------------------------------------------
c
c
c                                  specifications for arguments
      REAL  a,acc(2)
c                                  specifications for local variables
      REAL  x,y,z,zz
c                                  first executable statement
      x = acc(1)
      y = a
      if ( dabs(acc(1)).ge. dabs(a)) go to 1
      x = a
      y = acc(1)
c                                  compute z+zz = acc(1)+a exactly
    1 z = x+y
      zz = (x-z)+y
c                                  compute zz+acc(2) using double
c                                    precision arithmetic
      zz = zz+acc(2)
c                                  compute acc(1)+acc(2) = z+zz exactly
      acc(1) = z+zz
      acc(2) = (z-acc(1))+zz
      return
      end
      subroutine vxmul (a,b,acc)
c
c-----------------------------------------------------------------------
c
c   purpose             - extended precision multiply
c
c   usage               - call vxmul (a,b,acc)
c
c   arguments    a      - input real              number
c                b      - input real              number
c                acc    - accumulator. (input and output)
c                           acc is a real              vector of length
c                           2.  on output, acc contains the sum of
c                           input acc and a*b.
c
c   reqd. routines - vxadd
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c   remarks      vxmul adds the product a*b to the extended precision
c                accumulator, acc. the subroutine assumes that an
c                extended precision number is already in the
c                accumulator.  therefore, before the first call to
c                vxmul, acc(1) and acc(2) must be set to zero.
c
c-----------------------------------------------------------------------
c
c                                  specifications for arguments
      REAL  a,b,acc(2)
c                                  specifications for local variables
      REAL  x,ha,ta,hb,tb
      integer            ix(2),i
_IF(cray)
      logical          lx(8),li(4)
_ELSE
      logical*1          lx(8),li(4)
_ENDIF
      equivalence        (x,lx(1),ix(1)),(i,li(1))
      data               i/0/
c                                  split a = ha+ta
c                                        b = hb+tb
c                                  first executable statement
      x = a
      li(4) = lx(5)
      ix(2) = 0
      i = (i/16)*16
      lx(5) = li(4)
      ha=x
      ta=a-ha
      x = b
      li(4) = lx(5)
      ix(2) = 0
      i = (i/16)*16
      lx(5) = li(4)
      hb = x
      tb = b-hb
c                                  compute ha*hb,ha*tb,ta*hb, and ta*tb
c                                    and call vxadd to accumulate the
c                                    sum
      x = ta*tb
      call vxadd(x,acc)
      x = ha*tb
      call vxadd(x,acc)
      x = ta*hb
      call vxadd(x,acc)
      x = ha*hb
      call vxadd(x,acc)
      return
      end
      subroutine vxsto (acc,d)
c
c-----------------------------------------------------------------------
c
c   purpose             - real              store.
c
c   usage               - call vxsto(acc,d)
c
c   arguments    acc    - accumulator. (input)
c                           acc is a real              vector of length
c                           2. acc is assumed to be the result of
c                           calling vxadd or vxmul to perform extended
c                           precision operations.
c                d      - real              scalar. (output)
c                           on output, d contains a real
c                           approximation to the value of the extended
c                           precision accumulator.
c
c
c   reqd.  routines - none required
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c-----------------------------------------------------------------------
c
c                                  specifications for arguments
      REAL  acc(2),d
c                                  first executable statement
      d = acc(1)+acc(2)
      return
      end
      subroutine psa( amat,bmat,cmat,lmat,rvec,deriv,const,eigen,poleig,
     * nrpol,diag,work,range,crit,ndim,nd,ndim2,mdim,maxn,maxm,iprint,
     * ifile,maxpsa,ivec,ier,neigen )
c
c***********************************************************************
c
c     this subroutine calculates the eigenvalues ( ivec = 0 ) and the
c     eigenvectors ( ivec = 1 ) of a symmetric bordered matrix of the
c     form:
c
c                    amat(maxn,maxn)    bmat(maxn,maxm)
c                                          cmat(maxm)
c
c     in the interval range(1),range(2). (note: range(1) .lt. range(2))
c     the eigenvalues are ordered according to magnitude and are
c     written on file ifile. the eigenvectors are blocked with block-
c     length nd.
c
c***********************************************************************
c
c                                 i n p u t
c
c     amat(i)          n x n matrix a in symmetric storage mode
c     bmat(i,j)        n x m matrix of coupling matrix elements b
c     cmat(i)          m x m diagonal matrix c
c     const(i)         additive constant for the eigenvalues diverging
c                      at the diagonal elements of c
c                      i .le. maxn: left hand margin
c                      i .gt. maxn: right hand margin
c     crit             convergence criterion
c     deriv(i)         derivative of l matrix in symmetric storage mode
c     diag(i)          eigenvalues of l
c     eigen(i)         values of the eigenvalues diverging at the
c                      diagonal elements of c
c                      i .le. maxn: left hand margin
c                      i .gt. maxn: right hand margin
c     ifile            logical file number for output unit
c     ier              gives the information how many of the eigen-
c                      values could not be calculated (i.e. found)
c     iprint           printfile for the possible error messages
c                      it is advisable not to use the standard unit 6
c     ivec             ivec = 0   calculate only the eigenvalues
c                      ivec = 1   calculate eigenvalues and eigenvectors
c     lmat(i)          omega dependent n x n matrix l in symmetric
c                      storage mode
c     maxm             actual dimension of c
c                      maxm .le. mdim
c     maxn             actual dimension of a
c                      maxn .le. ndim
c     maxpsa           maximum number of iterations
c     mdim             dimension of cmat
c     nd               (ndim+1) * ndim / 2
c     ndim             dimension of amat
c     ndim2            2 * ndim
c     neigen           number of calculated eigenvalues
c     nrpol(i)         auxiliary vector for storing the diagonal
c                      elements of c
c     poleig(i)        eigenvalues of l at the diagonal elements of c
c                      cmatde
c                      i .le. maxn: left hand margin
c                      i .gt. maxn: right hand margin
c     range(i)         range(1): lower limit for the calculation of the
c                      eigenvalues
c                      range(2): upper limit for the calculation of the
c                      eigenvalues
c     rvec(i,j)        eigenvector matrix
c     work(i)          auxiliary vector for the diagonalization routine
c                      "eigrs"
c                      this is used simultaneously for the derivatives
c                      of the eigenvalues of l
c     zero             criterion for the degeneracy of the diagonal
c                      elements of c
c                      default value:  1.0e-12
c
c***********************************************************************
c
      implicit REAL  (a-h,o-z)
      REAL  lmat
c
      dimension amat(nd),bmat(ndim,mdim),cmat(mdim),lmat(nd),
     *          rvec(ndim,ndim),deriv(nd),const(ndim2),eigen(ndim2),
     *          poleig(ndim2),nrpol(mdim),diag(ndim),work(ndim),range(2)
c
    1 format(////////
     + ' ****     range(1) = ',f20.5,'  is larger than range(2)',
     +   f20.5,'       ***'//////)
c
      data d0,d05,d1,d4,d10 /0.0d0,0.5d0,1.0d0,4.0d0,1.0d10/
c
      max2n = (maxn+1) * maxn / 2
      neigen = 0
      ier = 0
      zero = 1.0d-12
      rewind ifile
c
c***********************************************************************
c
c     set range for the pole search
c
c***********************************************************************
c
      if(range(1) .gt. range(2)) go to 3000
c
c     determine all diagonal elements of c in the interval range(1)
c     to range(2)
c
      mxpol = 0
      do 20 i = 1,maxm
      if(cmat(i) .gt. range(2)) go to 20
      if(cmat(i) .lt. range(1)) go to 20
      mxpol = mxpol + 1
      nrpol(mxpol) = i
   20 continue
c
c     special treatment of the case no poles are found in the range
c
      if (mxpol.gt.0) go to 35
      do 30 i=1,maxm
      if (cmat(i).le.0.0d0) go to 30
      mxpol=1
      nrpol(mxpol)=i
      go to 35
   30 continue
   35 continue
c
c     consideration of boundary points
c
      minpol = nrpol(1) - 1
      maxpol = nrpol(mxpol) + 1
c
c***********************************************************************
c
c     determine the degeneracy
c
c***********************************************************************
c
      idegl = 0
      if(minpol .eq. 0) go to 60
      do 40 j = 1,minpol
      if( dabs(cmat(minpol) - cmat(j)) .gt. zero) go to 40
      idegl = idegl + 1
   40 continue
      minpol = nrpol(1) - idegl
   60 ideg = 0
      iadd = 1
      iz = 0
      do 120 i = 1,mxpol
      if(ideg .gt. 1) go to 100
      iz = iz + iadd
      ipol = nrpol(iz)
      jmin = iz
      ideg = 0
      iadd = 1
      jdeg = 1
      do 80 j = jmin,mxpol
      jpol = nrpol(j)
      if( dabs(cmat(ipol) - cmat(jpol)) .gt. zero) go to 80
      ideg = ideg + 1
   80 continue
      nrpol(iz) = nrpol(iz)*100 + ideg
      go to 120
  100 jdeg = jdeg + 1
      nrpol(i) = -10
      if(jdeg .ne. ideg) go to 120
      iadd = ideg
      ideg = 1
  120 continue
      idegr = 0
      if(maxpol .eq. maxm+1) go to 160
      do 140 j = maxpol,maxm
      if( dabs(cmat(maxpol) - cmat(j)) .gt. zero) go to 140
      idegr = idegr + 1
  140 continue
c
c***********************************************************************
c
c     calculation of the eigenvalues at the smallest diagonal element
c     of c in the interval range(1),range(2) or at the left hand margin
c
c***********************************************************************
c
  160 ipol = minpol
      ideg2 = idegl
      if(ipol .eq. 0) go to 300
      pol2 = cmat(ipol)
      ideg = ideg2
      call eipol( pol2,ipol,ideg2,amat,bmat,cmat,lmat,rvec,deriv,const,
     *            eigen,diag,work,zero,maxn,maxm,ndim,nd,ndim2,mdim )
      imax = ideg2
      if(imax .gt. maxn) imax = maxn
      imin = maxn - imax
      if(maxn .eq. 1) imin = 0
      if(imax .eq. maxn  .or.  maxn .eq. 1) go to 240
      do 220 i = 1,imin
      poleig(i) = diag(imax+i)
  220 continue
  240 do 260 i = 1,imax
      poleig(imin+i) =  dfloat(i) * d10
  260 continue
      do 280 i = 1,imax
      iz = maxn - i + 1
      const(iz) = const(maxn+i)
      eigen(iz) = eigen(maxn+i)
  280 continue
      go to 400
c
  300 ideg2 = 0
      do 320 i = 1,max2n
      lmat(i) = amat(i)
  320 continue
      call eigrsf(lmat,maxn,0,diag,rvec,ndim,work,jer)
      pol2 = (diag(1) + cmat(1) ) * d4
      polx = (diag(maxn) + cmat(maxm) ) * d4
      do 340 i = 1,25
      call eival(pol2,0,1,amat,bmat,cmat,lmat,rvec,deriv,diag,work,maxn,
     *           maxm,ndim,nd,mdim )
      if(diag(maxn) .gt. polx) polx = polx * d4
      if(diag(1) .gt. pol2) go to 360
      pol2 = pol2 * d4
  340 continue
  360 do 380 i = 1,maxn
      poleig(i) = diag(i)
      eigen(i) = d0
      const(i) = d0
  380 continue
c
c***********************************************************************
c
c     loop over diagonal elements of c, which are in the range
c     range(1) to range(2)
c
c***********************************************************************
c
  400 mrange = mxpol + 1
      do 1000 ieigen = 1,mrange
      if(ieigen .eq. mrange) go to 480
      if(nrpol(ieigen) .lt. 0) go to 1000
      if(ieigen+1 .eq. mrange) go to 420
      if(nrpol(ieigen+1) .lt. 0) go to 460
  420 xnull = d0
      ipol = nrpol(ieigen) / 100
      do 440 i = 1,maxn
      xnull = xnull + bmat(i,ipol) * bmat(i,ipol)
  440 continue
      iauf = 0
      if(xnull .le. zero) iauf = 1
      if(iauf .eq. 0) go to 460
      ipol = nrpol(ieigen) / 100
      call eideg( idiff,ideg2,ipol,bmat,lmat,deriv,rvec,work,cmat(ipol),
     *          zero,maxn,maxm,ndim,nd,mdim,ivec,ifile,iprint,neigen,1 )
      go to 1000
  460 pol1 = pol2
      ideg1 = ideg2
      ipol = nrpol(ieigen) / 100
      ideg2 = nrpol(ieigen) - 100 * ipol
      go to 500
  480 pol1 = pol2
      ideg1 = ideg2
      ideg2 = idegr
      ipol = maxpol
c
c     store eigenvalues at the right hand margin to the left hand margin
c
  500 if(ieigen .eq. 1) go to 600
      imax = ideg1
      if(imax .gt. maxn) imax = maxn
      imin = maxn - imax
      if(maxn .eq. 1) imin = 0
      if(imax .eq. maxn  .or.  maxn .eq. 1) go to 540
      do 520 i = 1,imin
      poleig(i) = poleig(maxn+imax+i)
  520 continue
  540 do 560 i = 1,imax
      poleig(i+imin) =  dfloat(i) * d10
  560 continue
      do 580 i = 1,imax
      iz = maxn - i + 1
      const(iz) = const(maxn+i)
      eigen(iz) = eigen(maxn+i)
  580 continue
c
c     calcuation of the eigenvalues at the right hand pole
c
  600 if(ipol .eq. maxm+1) go to 620
      ideg = ideg2
      pol2 = cmat(ipol)
      call eipol( pol2,ipol,ideg2,amat,bmat,cmat,lmat,rvec,deriv,const,
     *            eigen,diag,work,zero,maxn,maxm,ndim,nd,ndim2,mdim )
      go to 660
c
  620 pol2 = cmat(maxm) * d4
      if(cmat(1) .gt. range(1)) pol2 = polx
      ideg = ideg2
      do 640 i = 1,25
      call eival(pol2,0,1,amat,bmat,cmat,lmat,rvec,deriv,diag,work,maxn,
     *           maxm,ndim,nd,mdim )
      if(diag(maxn) .lt. pol2) go to 660
      pol2 = pol2 * d4
  640 continue
c
  660 do 680 i = 1,maxn
      poleig(maxn+i) = diag(i)
  680 continue
c
c***********************************************************************
c
c     calculate of the points of intersection in the interval pol1,pol2
c
c***********************************************************************
c
      do 900 j = 1,maxn
      if( (pol1-poleig(j)) * (pol2-poleig(maxn+j)) .gt. d0 ) go to 900
c
c     calculation of starting values
c
      if(j .le. ideg2) go to 700
      if(j .gt. maxn-ideg1) go to 720
c
c     1.) starting values for the non diverging eigenvalues
c
      x1 = pol1
      x2 = pol2
      y1 = poleig(j)
      y2 = poleig(maxn+j)
      slope = (y2-y1) / (x2-x1)
      ycut = y1  -  slope * x1
      xneu = ycut / (d1 - slope)
      go to 800
c
c     2.) starting values for the eigenvalues diverging to -infinity
c
  700 x1 = const(maxn+j) + pol2
      x2 = const(maxn+j) - pol2
      root =  dsqrt( x2**2 + d4*eigen(maxn+j) )
      omega1 = (x1 + root) * d05
      omega2 = (x1 - root) * d05
      xneu = pol1 - (pol1-pol2) * d05
      if(omega1 .lt. pol2  .and.  omega1 .gt. pol1) xneu = omega1
      if(omega2 .lt. pol2  .and.  omega2 .gt. pol1) xneu = omega2
      if(xneu .eq. omega1  .or.  xneu .eq. omega2)  go to 800
      if(ideg2 .gt. ideg1) go to 800
      x1 = const(maxn-j+1) + pol1
      x2 = const(maxn-j+1) - pol1
      root =  dsqrt( x2**2 + d4*eigen(maxn-j+1) )
      omega1 = (x1 + root) * d05
      omega2 = (x1 - root) * d05
      if(omega1 .gt. pol1  .and.  omega1 .lt. pol2) xneu = omega1
      if(omega2 .gt. pol1  .and.  omega2 .lt. pol2) xneu = omega2
      go to 800
c
c     3.) starting values for the eigenvalues diverging to +infinity
c
  720 x1 = const(j) + pol1
      x2 = const(j) - pol1
      root =  dsqrt( x2**2 + d4*eigen(j) )
      omega1 = (x1 + root) * d05
      omega2 = (x1 - root) * d05
      xneu = pol1 - (pol1-pol2) * d05
      if(omega1 .gt. pol1  .and.  omega1 .lt. pol2) xneu = omega1
      if(omega2 .gt. pol1  .and.  omega2 .lt. pol2) xneu = omega2
      if(xneu .eq. omega1  .or.  xneu .eq. omega2) go to 800
      if(ideg1 .gt. ideg2) go to 800
      x1 = const(2*maxn-j+1) + pol2
      x2 = const(2*maxn-j+1) - pol2
      root =  dsqrt( x2**2 + d4*eigen(2*maxn-j+1) )
      omega1 = (x1 + root) * d05
      omega2 = (x1 - root) * d05
      if(omega1 .lt. pol2  .and.  omega1 .gt. pol1) xneu = omega1
      if(omega2 .lt. pol2  .and.  omega2 .gt. pol1) xneu = omega2
c
c***********************************************************************
c
c     iteration
c
c***********************************************************************
c
  800 iauf = 1
      ieradd = 0
      difalt = d10
      do 860 k = 1,maxpsa
      call eival(xneu,1,j,amat,bmat,cmat,lmat,rvec,deriv,diag,work,maxn,
     *           maxm,ndim,nd,mdim )
      difneu =  dabs( xneu-diag(j) )
      if(difneu .lt. crit) go to 840
      if(difneu .gt. difalt) go to 820
      difalt = difneu
      xalt = xneu
      xneu = xalt  -  (diag(j) - xalt) / (work(j) - d1)
      if(xneu .gt. pol1  .and.  xneu .lt. pol2) go to 860
  820 call inter( j,pol1,pol2,poleig(j),poleig(maxn+j),xneu,amat,bmat,
     * cmat,lmat,rvec,deriv,diag,work,crit,maxn,maxm,ndim,nd,ndim2,mdim,
     * maxpsa,iprint,iauf,ieradd )
      if(ieradd .ne. 0) go to 880
      go to 860
  840 call eivec( j,bmat,cmat,lmat,rvec,work,diag(j),work(j),maxn,maxm,
     *            ndim,nd,mdim,ivec,ifile,neigen )
      go to 900
c
  860 continue
      ieradd=1
  880 ier = ier + ieradd
  900 continue
      idiff = ideg - ideg2
      if(idiff .le. 0) go to 1000
c     call eideg( idiff,ideg2,ipol,bmat,lmat,deriv,rvec,work,pol2,zero,
c    *            maxn,maxm,ndim,nd,mdim,ivec,ifile,iprint,neigen,0 )
      ier=ier+idiff
 1000 continue
c
      return
c
c***********************************************************************
c
c     error message printout
c
c***********************************************************************
c
 3000 ier = maxn + maxm
      write(iprint,1) range(1),range(2)
c
      return
      end
       subroutine   sdiag2 (m,n,a,d,x)
c
c
c
c      computation of all eigenvalues and eigenvectors of a real
c      symmetric matrix by the method of qr transformations.
c      if the euclidean norm of the rows varies   s t r o n g l y
c      most accurate results may be obtained by permuting rows and
c      columns to give an arrangement with increasing norms of rows.
c
c      two machine constants must be adjusted appropriately,
c      eps = minimum of all x such that 1+x is greater than 1 on the
c            computer,
c      tol = inf / eps  with inf = minimum of all positive x represen-
c            table within the computer.
c      a dimension statement e(mxtda4) may also be changed appropriately
c
c      input
c
c      (m)   not larger than 500,  corresponding value of the actual
c            dimension statement a(m,m), d(m), x(m,m),
c      (n)   not larger than (m), order of the matrix,
c      (a)   the matrix to be diagonalized, its lower triangle has to
c            be given as  ((a(i,j), j=1,i), i=1,n),
c
c      output
c
c      (d)   components d(1), ..., d(n) hold the computed eigenvalues
c            in ascending sequence. the remaining components of (d) are
c            unchanged,
c      (x)   the computed eigenvector corresponding to the j-th eigen-
c            value is stored as column (x(i,j), i=1,n). the eigenvectors
c            are normalized and orthogonal to working accuracy. the
c            remaining entries of (x) are unchanged.
c
c      array (a) is unaltered. however, the actual parameters
c      corresponding to (a) and (x)  may be identical, ''overwriting''
c      the eigenvectors on (a).
c
c      leibniz-rechenzentrum, munich 1965
c
c
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
      common/craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
      dimension   a(m,m), d(m), x(m,m)
      dimension e(mxtda4)
c
c     correct adjustment for ibm 360/91 real
c
      eps=2.5d-16
      tol=2.5d-63
c
c     counting parameter for timing in tamm-dancoff
      icdiag=icdiag+1
      if (icdiag.gt.500) icdiag=500
      tim13=cpulft(1)
c
      if(n.eq.1) go to 400
      do 10 i=1,n
      do 10 j=1,i
   10 x(i,j)=a(i,j)
c
c     householder's reduction
c     simulation of loop do 150 i=n,2,(-1)
c
      do 150 ni=2,n
      ii=n+2-ni
c     fake loop for recursive address calculation
      do 150 i=ii,ii
      l=i-2
      h=0.0d0
      g=x(i,i-1)
      if(l) 140,140,20
   20 do 30 k=1,l
   30 h=h+x(i,k)**2
      s=h+g*g
      if(s.ge.tol) go to 50
      h=0.0d0
      go to 140
   50 if(h) 140,140,60
   60 l=l+1
      f=g
      g= dsqrt(s)
      if(f) 75,75,70
   70 g=-g
   75 h=s-f*g
      x(i,i-1)=f-g
      f=0.0d0
c
      do 110 j=1,l
      x(j,i)=x(i,j)/h
      s=0.0d0
      do 80 k=1,j
   80 s=s+x(j,k)*x(i,k)
      j1=j+1
      if(j1.gt.l) go to 100
      do 90 k=j1,l
   90 s=s+x(k,j)*x(i,k)
  100 e(j)=s/h
  110 f=f+s*x(j,i)
c
      f=f/(h+h)
c
      do 120 j=1,l
  120 e(j)=e(j)-f*x(i,j)
c
      do 130 j=1,l
      f=x(i,j)
      s=e(j)
      do 130 k=1,j
  130 x(j,k)=x(j,k)-f*e(k)-x(i,k)*s
c
  140 d(i)=h
  150 e(i-1)=g
c
c     accumulation of transformation matrices
c
      d(1)=x(1,1)
      x(1,1)=1.0d0
      do 220 i=2,n
      l=i-1
      if(d(i)) 200,200,170
  170 do 190 j=1,l
      s=0.0d0
      do 180 k=1,l
  180 s=s+x(i,k)*x(k,j)
      do 190 k=1,l
  190 x(k,j)=x(k,j)-s*x(k,i)
  200 d(i)=x(i,i)
      x(i,i)=1.0d0
      do 220 j=1,l
      x(i,j)=0.0d0
  220 x(j,i)=0.0d0
c
c     diagonalization of the tridiagonal matrix
c
      b=0.0d0
      f=0.0d0
      e(n)=0.0d0
c
      do 340 l=1,n
      h=eps*( dabs(d(l))+ dabs(e(l)))
      if (h.gt.b) b=h
c
c     test for splitting
c
      do 240 j=l,n
      if ( dabs(e(j)).le.b) goto 250
  240 continue
c
c     test for convergence
c
  250 if(j.eq.l) go to 340
c
c     shift from upper 2*2 minor
c
  260 p=(d(l+1)-d(l))*0.5d0/e(l)
      r= dsqrt(p*p +1.0d0)
      if(p) 270,280,280
  270 p=p-r
      go to 290
  280 p=p+r
  290 h=d(l)-e(l)/p
      do 300 i=l,n
  300 d(i)=d(i)-h
      f=f+h
c
c     qr transformation
c
      p=d(j)
      c=1.0d0
      s=0.0d0
c
c     simulation of loop do 330 i=j-1,l,(-1)
c
      j1=j-1
      do 330 ni=l,j1
      ii=l+j1-ni
c     fake loop for recursive address calculation
      do 330 i=ii,ii
      g=c*e(i)
      h=c*p
c
c     protection against underflow of exponents
c
      if ( dabs(p).lt. dabs(e(i))) goto 310
      c=e(i)/p
      r= dsqrt(c*c+1.0d0)
      e(i+1)=s*p*r
      s=c/r
      c=1.0d0/r
      go to 320
  310 c=p/e(i)
      r= dsqrt(c*c+1.0d0)
      e(i+1)=s*e(i)*r
      s=1.0d0/r
      c=c/r
  320 p=c*d(i)-s*g
      d(i+1)=h+s*(c*g+s*d(i))
      do 330 k=1,n
      h=x(k,i+1)
      x(k,i+1)=x(k,i)*s+h*c
  330 x(k,i)=x(k,i)*c-h*s
c
      e(l)=s*p
      d(l)=c*p
      if ( dabs(e(l)).gt.b) go to 260
c
c     convergence
c
  340 d(l)=d(l)+f
c
c     ordering of eigenvalues
c
      ni=n-1
      do 380i=1,ni
      k=i
      p=d(i)
      j1=i+1
      do 360j=j1,n
      if(d(j).ge.p) goto 360
      k=j
      p=d(j)
  360 continue
      if (k.eq.i) goto 380
      d(k) =d(i)
      d(i)=p
      do 370 j=1,n
      p=x(j,i)
      x(j,i)=x(j,k)
  370 x(j,k)=p
  380 continue
      go to 410
c
c     special treatment of case n = 1
c
  400 d(1)=a(1,1)
      x(1,1)=1.0d0
  410 continue
c
      time13(icdiag)=cpulft(1)-tim13
      idimdi(icdiag)=n
      return
      end
      subroutine quart(nmaxb,nmaxf,v,g,gvec,
     *  indvec,jndvec,kndvec,lndvec,iwr)
c
      implicit REAL  (a-h,o-z)
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(100)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      dimension v(nmaxf),g(nmaxb,nmaxb),gvec(nmaxb)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
    1 format (/1x,104('-')/)
    2 format (/' energies of quartett states ( in ev ) for symmetry',
     + i3/)
    3 format (10(2x,f10.4))
    4 format (1x,
     +'the calculation of the quartett states is done without'/
     +1x,'configuration selection'/' the ',i3,
     +1x,'lowest orbitals are taken into account in the calculation')
c
      tim12=cpulft(1)
      iend=nocc-1
      if (nocc.eq.1) go to 1010
      nt=ntotx
   10 continue
      do 80 isym=1,nsym
      do 20 i=1,norb
      iorb=indorb(i)
      if (idorb(iorb).eq.isym) go to 30
   20 continue
      go to 80
   30 continue
      ic=0
      do 60 ia=nstart,nt
      iaa=idorb(ia)
      iim=imul(iaa,isym)
      do 50 i=1,iend
      ii=idorb(i)
      js=i+1
      do 40 j=js,nocc
      jj=idorb(j)
      if (iim.ne.imul(ii,jj)) go to 40
      ic=ic+1
   40 continue
   50 continue
   60 continue
      if (ic.le.nmaxb) go to 80
      nt=nt-1
c
c     check for degeneracy
c
      do 70 nn=1,6
      if ( dabs(epsil(nt)-epsil(nt+1)).gt.0.1d-03) go to 75
      nt=nt-1
   70 continue
   75 continue
      go to 10
   80 continue
      write (iwr,1)
      write (iwr,4) nt
      do 1000 isym=1,nsym
      do 120 i=1,norb
      iorb=indorb(i)
      if (idorb(iorb).eq.isym) go to 130
  120 continue
      go to 1000
  130 continue
      do 140 i=1,nmaxb
      do 140 j=1,nmaxb
      g(j,i)=0.0d0
  140 continue
      ic=0
      do 900 ia=nstart,nt
      iaa=idorb(ia)
      iim=imul(iaa,isym)
      do 890 i=1,iend
      ii=idorb(i)
      js=i+1
      do 880 j=js,nocc
      jj=idorb(j)
      if (iim.ne.imul(ii,jj)) go to 880
      ic=ic+1
      jc=0
      do 800 iap=nstart,nt
      iaap=idorb(iap)
      iimp=imul(iaap,isym)
      do 790 ip=1,iend
      iip=idorb(ip)
      jsp=ip+1
      do 780 jp=jsp,nocc
      jjp=idorb(jp)
      if (iimp.ne.imul(iip,jjp)) go to 780
      jc=jc+1
      if (jc.gt.ic) go to 780
      x=0.0d0
      if (ia.eq.iap) go to 200
      go to 300
  200 if (i.gt.ip) go to 210
      ind = indvx (j,jp,i,ip)
      jnd = indvx (j,ip,jp,i)
      go to 290
  210 if (j.le.jp) go to 240
      ind = indvx (j,jp,i,ip)
      if (jp.lt.i) go to 220
      jnd = indvx (j,ip,jp,i)
      go to 290
  220 if (i.eq.j) go to 230
      jnd = indvx (j,ip,i,jp)
      go to 290
  230 jnd = indvx (i,jp,j,ip)
      go to 290
  240 ind = indvx (jp,j,i,ip)
      jnd = indvx (jp,i,j,ip)
  290 continue
      x=x-v(ind)+v(jnd)
  300 if (i.eq.ip) go to 310
      go to 350
  310 if (j.lt.jp) go to 320
      ind = indvx (ia,iap,j,jp)
      go to 330
  320 ind = indvx (ia,iap,jp,j)
  330 x=x+v(ind)
  350 if (j.eq.jp) go to 360
      go to 400
  360 if (i.lt.ip) go to 370
      ind = indvx (ia,iap,i,ip)
      go to 380
  370 ind = indvx (ia,iap,ip,i)
  380 x=x+v(ind)
  400 if (i.eq.jp) go to 410
      go to 450
  410 ind = indvx (ia,iap,j,ip)
      x=x-v(ind)
  450 if (j.eq.ip) go to 460
      go to 500
  460 ind = indvx (ia,iap,jp,i)
      x=x-v(ind)
  500 g(ic,jc)=x
  780 continue
  790 continue
  800 continue
  880 continue
  890 continue
  900 continue
      do 910 i=1,ic
      do 910 j=i,ic
      g(i,j)=g(j,i)
  910 continue
      ic=0
      do 940 ia=nstart,nt
      iaa=idorb(ia)
      iim=imul(iaa,isym)
      epsi=epsil(ia)
      do 930 i=1,iend
      ii=idorb(i)
      epsj=epsi-epsil(i)
      js=i+1
      do 920 j=js,nocc
      jj=idorb(j)
      if (iim.ne.imul(ii,jj)) go to 920
      epsk=epsj-epsil(j)
      ic=ic+1
      g(ic,ic)=g(ic,ic)-epsk
  920 continue
  930 continue
  940 continue
      ndim=ic
      call sdiag2(nmaxb,ndim,g,gvec,g)
      write (iwr,2) isym
      do 950 i=1,ndim
      gvec(i)=gvec(i)*conver
  950 continue
      write (iwr,3) (gvec(i),i=1,ndim)
 1000 continue
 1010 continue
      time12=cpulft(1)-tim12
c
      return
      end
      subroutine psadia (numc,nmaxa,nmaxb,mxbg22,nmaxc,nmaxf,
     +                   v,g,wkarea,aa,bb,b,veca,vecx
_IF(ibm,vax)
     * ,iwr)
_ELSE
     * ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
c
      implicit REAL  (a-h,o-z)
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
_IFN(ibm,vax)
      common/linkup/iiii,jjjj,kkkk,llll
_ENDIF
      dimension v(nmaxf),g(nmaxb,nmaxb),wkarea(*),b(nmaxb,*),
     * aa(nmaxa,nmaxa),bb(nmaxa,nmaxa),veca(mxbg22),vecx(mxbg22,nmaxc)
c
_IFN(ibm,vax)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
_ENDIF
      tim14=cpulft(1)
      ic=0
      jc=0
      do 60 jsym=1,nsym
      do 20 i=1,norb
      ind=indorb(i)
      if (idorb(ind).eq.jsym) go to 30
   20 continue
      go to 60
   30 continue
      jend=nsymg(jsym)
      do 50 j=1,jend
      if (iuse(jsym,j).eq.0) go to 50
      kkend=ksgnum(jsym,j)
      do 40 ii=1,kkend
      ip=isgind(jsym,j,ii)
      do 40 jj=ii,kkend
      iq=isgind(jsym,j,jj)
      if (ip.le.nocc.and.iq.gt.nocc) ic=ic+1
      if (ip.gt.nocc.or.iq.le.nocc) jc=jc+1
   40 continue
   50 continue
   60 continue
      ickeep=ic
c     jckeep=jc
      kckeep=ic+jc
      do 120 i=1,nmaxb
      b(i,1)=0.0d0
      do 120 j=1,nmaxb
      g(j,i)=0.0d0
  120 continue
      rewind nf15
      do 900 isym=1,nsym
      do 310 i=1,norb
      ind=indorb(i)
      if (idorb(ind).eq.isym) go to 320
  310 continue
      go to 900
  320 continue
      jend=nsymg(isym)
      do 890 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 880
      kend=ksgnum(isym,jp)
      read (nf15) ix,iy,iz,jz,kz,it,(veca(i),i=1,kz),
     *          ((vecx(i,k),i=1,kz),k=1,it)
      nd=iz
      do 550 i=1,kend
      k=isgind(isym,jp,i)
      if (k.gt.nocc) go to 550
      ix=ndeg(k)
      xa=epsil(k)
      do 540 j=i,kend
      l=isgind(isym,jp,j)
      jx=ndeg(l)
      xb=epsil(l)
      x=0.0d0
      if (l.gt.nocc) go to 350
      do 340 n=1,nd
      x=x-vecx(n,j)*vecx(n,i)/((-veca(n)-xb)*(-veca(n)-xa))
  340 continue
      go to 370
  350 continue
      do 360 n=1,nd
      x=x+vecx(n,j)*vecx(n,i)/((-veca(n)-xa)*(xb-xa))
  360 continue
  370 continue
      ic=0
      jc=ickeep
      do 500 jsym=1,nsym
      do 400 ii=1,norb
      ind=indorb(ii)
      if (idorb(ind).eq.jsym) go to 410
  400 continue
      go to 500
  410 continue
      jjend=nsymg(jsym)
      do 480 jjp=1,jjend
      if (iuse(jsym,jjp).eq.0) go to 480
      kkend=ksgnum(jsym,jjp)
      do 470 ii=1,kkend
      ip=isgind(jsym,jjp,ii)
      do 460 jj=ii,kkend
      iq=isgind(jsym,jjp,jj)
      if (ip.le.nocc.and.iq.gt.nocc) go to 411
      if (ip.gt.nocc.or.iq.le.nocc) go to 412
  411 continue
      ic=ic+1
      kc=ic
      go to 413
  412 continue
      jc=jc+1
      kc=jc
  413 continue
      if (idega.eq.1) go to 440
      xc=0.0d0
      do 430 m=1,ix
      ia=indeg(k,m)
      iaa=idorb(ia)
      do 420 mm=1,jx
      ib=indeg(l,mm)
      if (iaa.ne.idorb(ib)) go to 420
      if (ia.eq.ib) go to 415
_IF(ibm,vax)
      ind1= indvor (ip,iq,ib,ia)
      ind2= indvor (ip,ib,iq,ia)
      ind3= indvor (ip,ia,iq,ib)
_ELSE
      call indvor (ip,iq,ib,ia)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ib,iq,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ia,iq,ib)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+4.0d0*v(ind1)-v(ind2)-v(ind3)
      go to 420
  415 continue
_IF(ibm,vax)
      ind1= indvor (ip,iq,ia,ia)
      ind2= indvor (ip,ia,iq,ia)
_ELSE
      call indvor (ip,iq,ia,ia)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ia,iq,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+2.0d0*v(ind1)-v(ind2)
  420 continue
  430 continue
      go to 450
  440 continue
      if (k.eq.l) go to 445
_IF(ibm,vax)
      ind1= indvor (ip,iq,l,k)
      ind2= indvor (ip,l,iq,k)
      ind3= indvor (ip,k,iq,l)
_ELSE
      call indvor (ip,iq,l,k)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,l,iq,k)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,k,iq,l)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=4.0d0*v(ind1)-v(ind2)-v(ind3)
      go to 450
  445 continue
_IF(ibm,vax)
      ind1= indvor (ip,iq,k,k)
      ind2= indvor (ip,k,iq,k)
_ELSE
      call indvor (ip,iq,k,k)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,k,iq,k)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=2.0d0*v(ind1)-v(ind2)
  450 continue
      b(kc,1)=b(kc,1)+x*xc
  460 continue
  470 continue
  480 continue
  500 continue
  540 continue
  550 continue
      nd=jz
      do 800 i=1,kend
      k=isgind(isym,jp,i)
      ix=ndeg(k)
      xa=epsil(k)
      do 790 j=i,kend
      l=isgind(isym,jp,j)
      if (l.le.nocc) go to 790
      jx=ndeg(l)
      xb=epsil(l)
      x=0.0d0
      if (k.gt.nocc) go to 610
      do 600 n=1,nd
      nn=n+iz
      x=x+vecx(nn,j)*vecx(nn,i)/((-veca(nn)-xb)*(xb-xa))
  600 continue
      go to 630
  610 continue
      do 620 n=1,nd
      nn=n+iz
      x=x+vecx(nn,j)*vecx(nn,i)/((-veca(nn)-xb)*(-veca(nn)-xa))
  620 continue
  630 continue
      ic=0
      jc=ickeep
      do 750 jsym=1,nsym
      do 650 ii=1,norb
      ind=indorb(ii)
      if (idorb(ind).eq.jsym) go to 660
  650 continue
      go to 750
  660 continue
      jjend=nsymg(jsym)
      do 740 jjp=1,jjend
      if (iuse(jsym,jjp).eq.0) go to 740
      kkend=ksgnum(jsym,jjp)
      do 730 ii=1,kkend
      ip=isgind(jsym,jjp,ii)
      do 720 jj=ii,kkend
      iq=isgind(jsym,jjp,jj)
      if (ip.le.nocc.and.iq.gt.nocc) go to 661
      if (ip.gt.nocc.or.iq.le.nocc) go to 662
  661 continue
      ic=ic+1
      kc=ic
      go to 663
  662 continue
      jc=jc+1
      kc=jc
  663 continue
      if (idega.eq.1) go to 700
      xc=0.0d0
      do 690 m=1,ix
      ia=indeg(k,m)
      iaa=idorb(ia)
      do 680 mm=1,jx
      ib=indeg(l,mm)
      if (iaa.ne.idorb(ib)) go to 680
      if (ia.eq.ib) go to 670
_IF(ibm,vax)
      ind1= indvor (ip,iq,ib,ia)
      ind2= indvor (ip,ib,iq,ia)
      ind3= indvor (ip,ia,iq,ib)
_ELSE
      call indvor (ip,iq,ib,ia)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ib,iq,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ia,iq,ib)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+4.0d0*v(ind1)-v(ind2)-v(ind3)
      go to 680
  670 continue
_IF(ibm,vax)
      ind1= indvor (ip,iq,ia,ia)
      ind2= indvor (ip,ia,iq,ia)
_ELSE
      call indvor (ip,iq,ia,ia)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ia,iq,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+2.0d0*v(ind1)-v(ind2)
  680 continue
  690 continue
      go to 710
  700 continue
      if (k.eq.l) go to 705
_IF(ibm,vax)
      ind1= indvor (ip,iq,l,k)
      ind2= indvor (ip,l,iq,k)
      ind3= indvor (ip,k,iq,l)
_ELSE
      call indvor (ip,iq,l,k)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,l,iq,k)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,k,iq,l)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=4.0d0*v(ind1)-v(ind2)-v(ind3)
      go to 710
  705 continue
_IF(ibm,vax)
      ind1= indvor (ip,iq,k,k)
      ind2= indvor (ip,k,iq,k)
_ELSE
      call indvor (ip,iq,k,k)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,k,iq,k)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=2.0d0*v(ind1)-v(ind2)
  710 continue
      b(kc,1)=b(kc,1)+x*xc
  720 continue
  730 continue
  740 continue
  750 continue
  790 continue
  800 continue
      go to 890
  880 continue
      read (nf15) dummy
  890 continue
  900 continue
      rewind nf15
      do 910 i=1,nmaxb
      do 910 j=1,nmaxb
      g(j,i)=0.0d0
  910 continue
      ic=0
      jc=ickeep
      do 300 isym=1,nsym
      do 140 i=1,norb
      ind=indorb(i)
      if (idorb(ind).eq.isym) go to 150
  140 continue
      go to 300
  150 continue
      jend=nsymg(isym)
      do 290 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 290
      kend=ksgnum(isym,jp)
      do 270 i=1,kend
      ip=isgind(isym,jp,i)
      do 260 j=i,kend
      iq=isgind(isym,jp,j)
      if (ip.le.nocc.and.iq.gt.nocc) go to 151
      if (ip.gt.nocc.or.iq.le.nocc) go to 152
  151 continue
      ic=ic+1
      kc=ic
      go to 153
  152 continue
      jc=jc+1
      kc=jc
  153 continue
      lc=0
      do 250 jsym=1,nsym
      do 170 ix=1,norb
      jnd=indorb(ix)
      if (idorb(jnd).eq.jsym) go to 180
  170 continue
      go to 250
  180 continue
      jjend=nsymg(jsym)
      do 240 jjp=1,jjend
      if (iuse(jsym,jjp).eq.0) go to 240
      kkend=ksgnum(jsym,jjp)
      do 230 ii=1,kkend
      k=isgind(jsym,jjp,ii)
      if (k.gt.nocc) go to 240
      do 220 jj=ii,kkend
      l=isgind(jsym,jjp,jj)
      if (l.le.nocc) go to 220
      lc=lc+1
      if (idega.eq.0) go to 192
_IF(ibm,vax)
      ind1= indvor (iq,ip,l,k)
      ind2= indvor (l,ip,iq,k)
      ind3= indvor (l,iq,ip,k)
_ELSE
      call indvor (iq,ip,l,k)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (l,ip,iq,k)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (l,iq,ip,k)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      g(kc,lc)=g(kc,lc)-(4.0d0*v(ind1)-v(ind2)-v(ind3))/
     *         (epsil(k)-epsil(l))
      go to 200
  192 continue
      ix=ndeg(k)
      jx=ndeg(l)
      xa=epsil(k)
      xb=epsil(l)
      xc=0.0d0
      do 196 m=1,ix
      ia=indeg(k,m)
      iaa=idorb(ia)
      do 194 mm=1,jx
      ib=indeg(l,mm)
      if (iaa.ne.idorb(ib)) go to 194
_IF(ibm,vax)
      ind1= indvor (ip,iq,ib,ia)
      ind2= indvor (ib,ip,iq,ia)
      ind3= indvor (ib,iq,ip,ia)
_ELSE
      call indvor (ip,iq,ib,ia)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,ip,iq,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,iq,ip,ia)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+4.0d0*v(ind1)-v(ind2)-v(ind3)
  194 continue
  196 continue
      g(kc,lc)=g(kc,lc)-xc/(xa-xb)
  200 continue
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
  270 continue
  290 continue
  300 continue
c
      do 301 i=1,ickeep
      g(i,i)=g(i,i)+1.0d0
  301 continue
c
      call leqt2f (g,1,ickeep,nmaxb,b,6,wkarea,ier)
      if (ier.eq.34.or.ier.eq.129) go to 1200
c
c     calculate the second part
c
      is=ickeep+1
      do 912 i=is,kckeep
      x=0.0d0
      do 911 j=1,ickeep
      x=x-g(i,j)*b(j,1)
  911 continue
      b(i,1)=b(i,1)+x
  912 continue
c
c     store constant diagrams in matrix aa(i,j)
c
      ic=0
      jc=ickeep
      do 1000 isym=1,nsym
      do 950 i=1,norb
      ind=indorb(i)
      if (idorb(ind).eq.isym) go to 960
  950 continue
      go to 1000
  960 continue
      jend=nsymg(isym)
      do 990 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 990
      kend=ksgnum(isym,jp)
      do 980 i=1,kend
      ip=isgind(isym,jp,i)
      do 970 j=i,kend
      iq=isgind(isym,jp,j)
      if (ip.le.nocc.and.iq.gt.nocc) go to 961
      if (ip.gt.nocc.or.iq.le.nocc) go to 962
  961 continue
      ic=ic+1
      kc=ic
      go to 963
  962 continue
      jc=jc+1
      kc=jc
  963 continue
      aa(ip,iq)=b(kc,1)
      aa(iq,ip)=b(kc,1)
  970 continue
  980 continue
  990 continue
 1000 continue
      do 1100 i=1,nmaxa
      do 1100 j=1,nmaxa
      bb(j,i)=0.0d0
 1100 continue
      go to 1300
 1200 continue
      numc=4
      fmul(1)=0.4d0
      fmul(2)=0.6d0
      fmul(3)=0.8d0
      do 1220 i=numc,numend
      fmul(i)=1.0d0
 1220 continue
      do 1240 i=1,nmaxa
      do 1240 j=1,nmaxa
      aa(j,i)=0.0d0
      bb(j,i)=0.0d0
 1240 continue
 1300 continue
c     test begin
c     write (iwr,1)
c   1 format ('1')
c     write (iwr,2)
c   2 format ('   a-diagrams computed with h-f greens function'//)
c     do 1050 i=1,ntotx
c     do 1010 j=1,ntotx
c     sigdia(j)=aa(j,i)*conver
c1010 continue
c     write (iwr,3) (sigdia(j),j=1,ntotx)
c1050 continue
c   3 format (/10(2x,f10.6))
c     test  end
      time14=cpulft(1)-tim14
      return
      end
      subroutine adiagr (numc,nmaxa,nmaxb,nmaxc,nmaxf,
     +                   v,g,wkarea,gvec,aux,aa,bb,b
_IF(ibm,vax)
     + ,iwr)
_ELSE
     + ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
c
      implicit REAL  (a-h,o-z)
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
_IFN(ibm,vax)
      common/linkup/iiii,jjjj,kkkk,llll
_ENDIF
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
      dimension v(nmaxf),g(nmaxb,nmaxb),wkarea(*)
      dimension gvec(nmaxb),aux(nmaxc,nmaxb),b(nmaxb,*)
      dimension aa(nmaxa,nmaxa),bb(nmaxa,nmaxa)
c
_IFN(ibm,vax)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
_ENDIF
      tim14=cpulft(1)
      ic=0
      jc=0
      do 60 jsym=1,nsym
      do 20 i=1,norb
      ind=indorb(i)
      if (idorb(ind).eq.jsym) go to 30
   20 continue
      go to 60
   30 continue
      jend=nsymg(jsym)
      do 50 j=1,jend
      if (iuse(jsym,j).eq.0) go to 50
      kkend=ksgnum(jsym,j)
      do 40 ii=1,kkend
      ip=isgind(jsym,j,ii)
      do 40 jj=ii,kkend
      iq=isgind(jsym,j,jj)
      if (ip.le.nocc.and.iq.gt.nocc) ic=ic+1
      if (ip.gt.nocc.or.iq.le.nocc) jc=jc+1
   40 continue
   50 continue
   60 continue
      ickeep=ic
c     jckeep=jc
      kckeep=ic+jc
      do 120 i=1,nmaxb
      b(i,1)=0.0d0
      do 120 j=1,nmaxb
      g(j,i)=0.0d0
  120 continue
      rewind nf17
      rewind nf13
      ic=0
      jc=ickeep
      do 300 isym=1,nsym
      do 140 i=1,norb
      ind=indorb(i)
      if (idorb(ind).eq.isym) go to 150
  140 continue
      go to 300
  150 continue
      jend=nsymg(isym)
      do 290 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 290
      kend=ksgnum(isym,jp)
      do 270 i=1,kend
      ip=isgind(isym,jp,i)
      do 260 j=i,kend
      iq=isgind(isym,jp,j)
      if (ip.le.nocc.and.iq.gt.nocc) go to 151
      if (ip.gt.nocc.or.iq.le.nocc) go to 152
  151 continue
      ic=ic+1
      kc=ic
      go to 153
  152 continue
      jc=jc+1
      kc=jc
  153 continue
      lc=0
      do 250 jsym=1,nsym
      do 170 ix=1,norb
      jnd=indorb(ix)
      if (idorb(jnd).eq.jsym) go to 180
  170 continue
      go to 250
  180 continue
      jjend=nsymg(jsym)
      do 240 jjp=1,jjend
      if (iuse(jsym,jjp).eq.0) go to 240
      kkend=ksgnum(jsym,jjp)
      do 230 ii=1,kkend
      k=isgind(jsym,jjp,ii)
      if (k.gt.nocc) go to 240
      do 220 jj=ii,kkend
      l=isgind(jsym,jjp,jj)
      if (l.le.nocc) go to 220
      lc=lc+1
      if (idega.eq.0) go to 192
_IF(ibm,vax)
      ind1= indvor (iq,ip,l,k)
      ind2= indvor (l,ip,iq,k)
      ind3= indvor (l,iq,ip,k)
_ELSE
      call indvor (iq,ip,l,k)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (l,ip,iq,k)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (l,iq,ip,k)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      g(kc,lc)=g(kc,lc)-(4.0d0*v(ind1)-v(ind2)-v(ind3))/
     *         (epsil(k)-epsil(l))
      go to 200
  192 continue
      ix=ndeg(k)
      jx=ndeg(l)
      xa=epsil(k)
      xb=epsil(l)
      xc=0.0d0
      do 196 m=1,ix
      ia=indeg(k,m)
      iaa=idorb(ia)
      do 194 mm=1,jx
      ib=indeg(l,mm)
      if (iaa.ne.idorb(ib)) go to 194
_IF(ibm,vax)
      ind1= indvor (ip,iq,ib,ia)
      ind2= indvor (ib,ip,iq,ia)
      ind3= indvor (ib,iq,ip,ia)
_ELSE
      call indvor (ip,iq,ib,ia)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,ip,iq,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,iq,ip,ia)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+4.0d0*v(ind1)-v(ind2)-v(ind3)
  194 continue
  196 continue
      g(kc,lc)=g(kc,lc)-xc/(xa-xb)
  200 continue
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
  270 continue
  290 continue
  300 continue
      do 900 isym=1,nsym
      do 310 i=1,norb
      ind=indorb(i)
      if (idorb(ind).eq.isym) go to 320
  310 continue
      go to 900
  320 continue
      jend=nsymg(isym)
      do 890 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 880
      kend=ksgnum(isym,jp)
      read (nf17) nd,(gvec(i),i=1,nd)
      read (nf17) is,nd,((aux(j,i),j=1,is),i=1,nd)
      do 550 i=1,kend
      k=isgind(isym,jp,i)
      if (k.gt.nocc) go to 550
      ix=ndeg(k)
      xa=epsil(k)
      do 540 j=i,kend
      l=isgind(isym,jp,j)
      jx=ndeg(l)
      xb=epsil(l)
      x=0.0d0
      if (l.gt.nocc) go to 350
      do 340 n=1,nd
      x=x-aux(j,n)*aux(i,n)/((gvec(n)-xb)*(gvec(n)-xa))
  340 continue
      go to 370
  350 continue
      do 360 n=1,nd
      x=x+aux(j,n)*aux(i,n)/((gvec(n)-xa)*(xb-xa))
  360 continue
  370 continue
      ic=0
      jc=ickeep
      do 500 jsym=1,nsym
      do 400 ii=1,norb
      ind=indorb(ii)
      if (idorb(ind).eq.jsym) go to 410
  400 continue
      go to 500
  410 continue
      jjend=nsymg(jsym)
      do 480 jjp=1,jjend
      if (iuse(jsym,jjp).eq.0) go to 480
      kkend=ksgnum(jsym,jjp)
      do 470 ii=1,kkend
      ip=isgind(jsym,jjp,ii)
      do 460 jj=ii,kkend
      iq=isgind(jsym,jjp,jj)
      if (ip.le.nocc.and.iq.gt.nocc) go to 411
      if (ip.gt.nocc.or.iq.le.nocc) go to 412
  411 continue
      ic=ic+1
      kc=ic
      go to 413
  412 continue
      jc=jc+1
      kc=jc
  413 continue
      if (idega.eq.1) go to 440
      xc=0.0d0
      do 430 m=1,ix
      ia=indeg(k,m)
      iaa=idorb(ia)
      do 420 mm=1,jx
      ib=indeg(l,mm)
      if (iaa.ne.idorb(ib)) go to 420
      if (ia.eq.ib) go to 415
_IF(ibm,vax)
      ind1= indvor (ip,iq,ib,ia)
      ind2= indvor (ip,ib,iq,ia)
      ind3= indvor (ip,ia,iq,ib)
_ELSE
      call indvor (ip,iq,ib,ia)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ib,iq,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ia,iq,ib)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+4.0d0*v(ind1)-v(ind2)-v(ind3)
      go to 420
  415 continue
_IF(ibm,vax)
      ind1= indvor (ip,iq,ia,ia)
      ind2= indvor (ip,ia,iq,ia)
_ELSE
      call indvor (ip,iq,ia,ia)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ia,iq,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+2.0d0*v(ind1)-v(ind2)
  420 continue
  430 continue
      go to 450
  440 continue
      if (k.eq.l) go to 445
_IF(ibm,vax)
      ind1= indvor (ip,iq,l,k)
      ind2= indvor (ip,l,iq,k)
      ind3= indvor (ip,k,iq,l)
_ELSE
      call indvor (ip,iq,l,k)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,l,iq,k)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,k,iq,l)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=4.0d0*v(ind1)-v(ind2)-v(ind3)
      go to 450
  445 continue
_IF(ibm,vax)
      ind1= indvor (ip,iq,k,k)
      ind2= indvor (ip,k,iq,k)
_ELSE
      call indvor (ip,iq,k,k)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,k,iq,k)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=2.0d0*v(ind1)-v(ind2)
  450 continue
      b(kc,1)=b(kc,1)+x*xc
  460 continue
  470 continue
  480 continue
  500 continue
  540 continue
  550 continue
      read (nf13) nd,(gvec(i),i=1,nd)
      read (nf13) is,nd,((aux(j,i),j=1,is),i=1,nd)
      do 800 i=1,kend
      k=isgind(isym,jp,i)
      ix=ndeg(k)
      xa=epsil(k)
      do 790 j=i,kend
      l=isgind(isym,jp,j)
      if (l.le.nocc) go to 790
      jx=ndeg(l)
      xb=epsil(l)
      x=0.0d0
      if (k.gt.nocc) go to 610
      do 600 n=1,nd
      x=x+aux(j,n)*aux(i,n)/((gvec(n)-xb)*(xb-xa))
  600 continue
      go to 630
  610 continue
      do 620 n=1,nd
      x=x+aux(j,n)*aux(i,n)/((gvec(n)-xb)*(gvec(n)-xa))
  620 continue
  630 continue
      ic=0
      jc=ickeep
      do 750 jsym=1,nsym
      do 650 ii=1,norb
      ind=indorb(ii)
      if (idorb(ind).eq.jsym) go to 660
  650 continue
      go to 750
  660 continue
      jjend=nsymg(jsym)
      do 740 jjp=1,jjend
      if (iuse(jsym,jjp).eq.0) go to 740
      kkend=ksgnum(jsym,jjp)
      do 730 ii=1,kkend
      ip=isgind(jsym,jjp,ii)
      do 720 jj=ii,kkend
      iq=isgind(jsym,jjp,jj)
      if (ip.le.nocc.and.iq.gt.nocc) go to 661
      if (ip.gt.nocc.or.iq.le.nocc) go to 662
  661 continue
      ic=ic+1
      kc=ic
      go to 663
  662 continue
      jc=jc+1
      kc=jc
  663 continue
      if (idega.eq.1) go to 700
      xc=0.0d0
      do 690 m=1,ix
      ia=indeg(k,m)
      iaa=idorb(ia)
      do 680 mm=1,jx
      ib=indeg(l,mm)
      if (iaa.ne.idorb(ib)) go to 680
      if (ia.eq.ib) go to 670
_IF(ibm,vax)
      ind1= indvor (ip,iq,ib,ia)
      ind2= indvor (ip,ib,iq,ia)
      ind3= indvor (ip,ia,iq,ib)
_ELSE
      call indvor (ip,iq,ib,ia)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ib,iq,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ia,iq,ib)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+4.0d0*v(ind1)-v(ind2)-v(ind3)
      go to 680
  670 continue
_IF(ibm,vax)
      ind1= indvor (ip,iq,ia,ia)
      ind2= indvor (ip,ia,iq,ia)
_ELSE
      call indvor (ip,iq,ia,ia)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,ia,iq,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+2.0d0*v(ind1)-v(ind2)
  680 continue
  690 continue
      go to 710
  700 continue
      if (k.eq.l) go to 705
_IF(ibm,vax)
      ind1= indvor (ip,iq,l,k)
      ind2= indvor (ip,l,iq,k)
      ind3= indvor (ip,k,iq,l)
_ELSE
      call indvor (ip,iq,l,k)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,l,iq,k)
      ind2=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,k,iq,l)
      ind3=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=4.0d0*v(ind1)-v(ind2)-v(ind3)
      go to 710
  705 continue
_IF(ibm,vax)
      ind2= indvor (ip,k,iq,k)
      ind1= indvor (ip,iq,k,k)
_ELSE
      call indvor (ip,iq,k,k)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,k,iq,k)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=2.0d0*v(ind1)-v(ind2)
  710 continue
      b(kc,1)=b(kc,1)+x*xc
  720 continue
  730 continue
  740 continue
  750 continue
  790 continue
  800 continue
      go to 890
  880 continue
      read (nf17) dummy
      read (nf17) dummy
      read (nf13) dummy
      read (nf13) dummy
  890 continue
  900 continue
      rewind nf17
      rewind nf13
c
      do 301 i=1,ickeep
      g(i,i)=g(i,i)+1.0d0
  301 continue
c
      call leqt2f (g,1,ickeep,nmaxb,b,6,wkarea,ier)
      if (ier.eq.34.or.ier.eq.129) go to 1200
c
c     calculate the second part
c
      is=ickeep+1
      do 912 i=is,kckeep
      x=0.0d0
      do 911 j=1,ickeep
      x=x-g(i,j)*b(j,1)
  911 continue
      b(i,1)=b(i,1)+x
  912 continue
c
c     store constant diagrams in matrix aa(i,j)
c
      ic=0
      jc=ickeep
      do 1000 isym=1,nsym
      do 950 i=1,norb
      ind=indorb(i)
      if (idorb(ind).eq.isym) go to 960
  950 continue
      go to 1000
  960 continue
      jend=nsymg(isym)
      do 990 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 990
      kend=ksgnum(isym,jp)
      do 980 i=1,kend
      ip=isgind(isym,jp,i)
      do 970 j=i,kend
      iq=isgind(isym,jp,j)
      if (ip.le.nocc.and.iq.gt.nocc) go to 961
      if (ip.gt.nocc.or.iq.le.nocc) go to 962
  961 continue
      ic=ic+1
      kc=ic
      go to 963
  962 continue
      jc=jc+1
      kc=jc
  963 continue
      aa(ip,iq)=b(kc,1)
      aa(iq,ip)=b(kc,1)
  970 continue
  980 continue
  990 continue
 1000 continue
      do 1100 i=1,nmaxa
      do 1100 j=1,nmaxa
      bb(j,i)=0.0d0
 1100 continue
      go to 1300
 1200 continue
      numc=4
      fmul(1)=0.4d0
      fmul(2)=0.6d0
      fmul(3)=0.8d0
      do 1220 i=numc,numend
      fmul(i)=1.0d0
 1220 continue
      do 1240 i=1,nmaxa
      do 1240 j=1,nmaxa
      aa(j,i)=0.0d0
      bb(j,i)=0.0d0
 1240 continue
 1300 continue
_IF()
c     test begin
c     write (iwr,1)
c   1 format ('1')
c     write (iwr,2)
c   2 format ('   a-diagrams computed with h-f greens function'//)
c     do 1050 i=1,ntotx
c     do 1010 j=1,ntotx
c     sigdia(j)=aa(j,i)*conver
c1010 continue
c     write (iwr,3) (sigdia(j),j=1,ntotx)
c1050 continue
c   3 format (/10(2x,f10.6))
c     test  end
_ENDIF
      time14=cpulft(1)-tim14
      return
      end
_EXTRACT(tdaf1,hp800)
_IF1(gbr)c $head$ tdaf_1
      subroutine comdig(nmaxa,nmaxb,nmaxc,nmaxf,nmaxi,iveca,ivecb,
     +                  v,g,wkarea,gvec,aux,aa,bb
_IF(ibm,vax)
     +,iwr )
_ELSE
     + ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iveca,ivecb
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
      character*8    selct,typ1,typ2,typ3,typ4
c     character*8    x
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
_IFN(ibm,vax)
      common/linkup/iiii,jjjj,kkkk,llll
_ENDIF
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
      dimension v(nmaxf),g(nmaxb,nmaxb),wkarea(*)
      dimension gvec(nmaxb),aux(nmaxc,nmaxb),iveca(8,nmaxi,*)
      dimension ivecb(8,nmaxi,*)
      dimension aa(nmaxa,nmaxa),bb(nmaxa,nmaxa)
_IFN(ibm,vax)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
_ENDIF
      data selct  /'adiagram'/
      data typ1   /'m'/
      data typ2   /'i'/
      data typ3   /'a'/
      data typ4   /'c'/
_IFN(ibm,vax)
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
_ENDIF
    1 format (/1x,104('-')/)
    2 format (//'   convergence obtained at iteration step   ',i3)
    3 format (' energies and pole strengths in nondiagonal',
     +1x,'approximation for symmetry ',i3,' and subpart ',i3//)
    4 format (' solution ',i3,3x,f9.3,10f9.6/(29x,10f9.6))
    5 format (' the following notation is used : m : orbital index, ',
     + 'i : 2h1p index, a : 2p1h index and c : centroid'//)
    6 format (1x,
     +'large components of eigenvectors in nondiagonal approximation')
    7 format (' eigenvector',i4,'::',10(1x,a1,i3,':',f5.2))
    8 format (/' diagonal of matrix of constant diagrams'/)
    9 format (10(2x,f10.6))
   10 format (/' convergence is not obtained'/
     +1x,'the intermediate results are printed after the iterations'//)
   11 format (///' results of first step without constant diagrams')
   12 format (/' for symmetry ',i3,' in the interval from  ',f10.5,
     *'  to  ',f10.5,' with solution of the dyson equation  ',f10.5/
     *' the sum of the squared 1h and 2h1p components of this solution',
     *' is  ',f10.5)
   13 format (//' correlation energy in first step'/
     +1x,'(without constant diagrams) = ',f10.6,
     + ' a.u. = ',f12.6,' ev'///)
   14 format (/' correlation energy including constant diagrams =',
     *f10.6,' a.u. = ',f12.6,' ev'///)
   15 format (//' sum of pole strengths of occupied orbital number',
     *i5,' in the given symmetry : ',f6.4)
   16 format (' energies of states with a pole strength larger'/
     +1x,'than 0.008 for any occupied orbital'//)
   17 format (//' results of second step where the constant'/1x,
     +1x,'diagrams are calculated with the hartree-fock',
     +' greens function')
_IF()
c  18 format (' set of important eigenvectors'//)
c  19 format (//'  vector nr',10(4x,i4,4x))
c  20 format (' eigenvalue',10(f12.5))
c  21 format (/)
c  22 format (2x,i3,2x,a1,i3,10(f12.7))
c  23 format (//'   in the tabulation of the eigenvectors below three qu
c    *antities are given for the row index identification'/
c    *'   the first one is a simple running index, the second one denote
c    *s the type of configuration which is delt with'/
c    *'   and the third one is again an index referring to the group spe
c    *cified by the second quantity'/
c    *'   m denotes a molecular orbital (occupied or virtual) and the in
c    *dex following it is the index of this orbital'/
c    *'   i denotes a 2h1p configuration, the index following it refers
c    *to the table given earlier under :'/
c    *'     eigenvalues and large components of the eigenvectors of the
c    *matrix for b (for the given symmetry) which again refers to the:'/
c    *'     index identification of the 2h1p configurations (for the giv
c    *en symmetry)')
c  24 format('   a denotes a 2p1h configuration, the index following it
c    *refers to the table given earlier under:'/
c    *'     eigenvalues and large components of the eigenvectors of the
c    *matrix for a (for the given symmetry) which again refers to :'/
c    *'     index identification of the 2p1h configurations (for the giv
c    *en symmetry)'/
c    *'   a few of the i and most of the a terms will be missing,they ha
c    *ve been among others incorporated into centroids'/
c    *'   c denotes a centroid, no index identification for these effect
c    *ive configurations can be given'/
c    *'   as they are composite quantities'//)
_ENDIF
   25 format (/' results of first step with constant diagrams'/
     +' calculated with the hartree-fock greens function')
c
      tim11=cpulft(1)
      time11=0.0d0
      iunit=nf18
      junit=nf19
      do 50 i=1,nsym
      do 50 j=1,nmaxi
      ixxx(i,j)=0
   50 continue
      numc=1
      do 60 i=1,numend
      fmul(i)=1.0d0
   60 continue
      do 70 i=1,nmaxa
      do 70 j=1,nmaxa
      aa(j,i)=0.0d0
      bb(j,i)=0.0d0
   70 continue
      if ( numend.eq.-2) nxcont=1
      if ( numend.eq.-1) nxcont=2
      if ( numend.eq.2)  nxcont=3
      if ( numend.eq.6)  nxcont=4
      if ( numend.lt.2) numend=1
      do 2000 num=1,numend
      iprint=0
      if ( num.ne.1.or.nxcont.ne.2) go to 4000
 4010 call adiagr(numc,nmaxa,nmaxb,nmaxc,nmaxf,
     * v,g,wkarea,gvec,aux,aa,bb,bb
_IF(ibm,vax)
     *,iwr )
_ELSE
     * ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
      go to 4020
 4000 if(num.eq.2)go to 4010
 4020 rewind nf13
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      do 600 isym=1,nsym
      do 100 i=1,norb
      ind=indorb(i)
      if (idorb(ind).eq.isym) go to 110
  100 continue
      go to 600
  110 continue
      jend=nsymg(isym)
      do 560 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 540
      do 120 i=1,nmaxb
      do 120 j=1,nmaxb
      g(j,i)=0.0d0
  120 continue
      isum=isgnum(isym,jp)
      ndim=idimb(isym,jp)+idima(isym,jp)+isum
      kend=isgnum(isym,jp)
      kkend=ksgnum(isym,jp)
      do 140 k=1,kend
      ind=isgind(isym,jp,k)
      do 140 l=k,kend
      jnd=isgind(isym,jp,l)
      g(k,l)=-aa(ind,jnd)
      g(l,k)=-aa(ind,jnd)
      if(k.eq.l) g(l,l)=g(l,l)-epsil(jnd)
  140 continue
      if (iprint.ne.0) go to 160
      write (iwr,8)
      do 150 i=1,ntot
      sigdia(i)=aa(i,i)*conver
  150 continue
      write (iwr,9) (sigdia(i),i=1,ntot)
      iprint=1
  160 continue
      read (nf17) nd,(gvec(i),i=1,nd)
      read (nf17) is,nd,((aux(j,i),j=1,is),i=1,nd)
      xa=-gvec(1)
      ks=idimb(isym,jp)+isum+1
      ke=ks+nd-1
      nc=0
      do 180 k=ks,ke
      nc=nc+1
      g(k,k)=-gvec(nc)
  180 continue
      do 200 l=1,kend
      nc=0
      do 200 k=ks,ke
      nc=nc+1
      xx=aux(l,nc)
      g(l,k)=xx
      g(k,l)=xx
  200 continue
      read (nf13) nd,(gvec(ir),ir=1,nd)
      read (nf13) is,nd,((aux(jr,ir),jr=1,is),ir=1,nd)
      ndd=jdimb(isym,jp)
      xb=-gvec(ndd)
      ls=isum+1
      le=nd+isum
      nc=0
      do 220 l=ls,le
      nc=nc+1
      g(l,l)=-gvec(nc)
  220 continue
      do 240 k=1,kend
      nc=0
      do 240 l=ls,le
      nc=nc+1
      xx=aux(k,nc)
      g(k,l)=xx
      g(l,k)=xx
  240 continue
      time11=time11+cpulft(1)-tim11
      call sdiag2(nmaxb,ndim,g,gvec,g)
      tim11=cpulft(1)
      isuma=0
      do 260 k=1,kend
      if (isgind(isym,jp,k).gt.nocc) go to 260
      isuma=isuma+1
  260 continue
c
c     the separation of the affinity and ionization poles of the green's
c     function is done in the following way:
c
c     if the energy is less than that of the first affinity pole of the
c     self-energy, the solution is an affinity pole of the green's
c     function. the interval between the first affinity pole and the
c     first ionization pole of the self-energy is investigated. let n
c     be the eigenvector of the green's function in this interval. if
c     the sum of the squared 1-hole and 2-hole-1-particle components of
c     this eigenvector ( c(p,n)**2 ) is less than 0.60 the solution is
c     taken as an affinity pole, otherwise as an ionization pole of
c     the green's function.
c
      do 280 nn=1,ndim
      nc=nn
      if (gvec(nn).le.xa) go to 280
      go to 300
  280 continue
  300 continue
      nsa=isum+1
      nsb=idimb(isym,jp)+isum
      xa=xa*conver
      xb=xb*conver
  320 continue
      xx=0.0d0
      if (isuma.eq.0) go to 350
      do 340 ip=1,isuma
      xx=xx+g(ip,nc)**2
  340 continue
  350 continue
      do 360 ip=nsa,nsb
      xx=xx+g(ip,nc)**2
  360 continue
      xc=gvec(nc)*conver
      write (iwr,12) isym,xa,xb,xc,xx
      if (xx.gt.0.60d0) go to 380
      nc=nc+1
      if (gvec(nc)*conver.gt.xb) go to 380
      go to 320
  380 continue
      ixxx(isym,jp)=nc
      if (word.ne.selct) go to 480
      if (num.lt.2) go to 480
      if (numend.eq.2) go to 480
      do 470 i=1,kkend
      ind=isgind(isym,jp,i)
      ix=ndeg(ind)
      do 460 j=1,kkend
      jnd=isgind(isym,jp,j)
      jx=ndeg(jnd)
      xx=0.0d0
      do 390 nn=nc,ndim
      xx=xx+g(i,nn)*g(j,nn)
  390 continue
      if (ind.eq.jnd.and.ind.le.nocc) xx=xx-1.0d0
      do 450 isymp=1,nsym
      ipp=nsymg(isymp)
      do 440 ip=1,ipp
      if (iuse(isymp,ip).eq.0) go to 440
      lend=ksgnum(isymp,ip)
      do 430 k=1,lend
      knd=isgind(isymp,ip,k)
      do 420 l=k,lend
      lnd=isgind(isymp,ip,l)
      if (idega.eq.1) go to 412
      xa=0.0d0
      do 410 m=1,ix
      ia=indeg(ind,m)
      iaa=idorb(ia)
      do 400 mm=1,jx
      ib=indeg(jnd,mm)
      if (iaa.ne.idorb(ib)) go to 400
_IF(ibm,vax)
      ind1= indvor (knd,lnd,ia,ib)
      ind2= indvor (knd,ib,lnd,ia)
_ELSE
      call indvor (knd,lnd,ia,ib)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (knd,ib,lnd,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xa=xa+v(ind1)*2.0d0-v(ind2)
  400 continue
  410 continue
      go to 418
  412 continue
_IF(ibm,vax)
      ind1= indvor (knd,lnd,ind,jnd)
      ind2= indvor (knd,jnd,lnd,ind)
_ELSE
      call indvor (knd,lnd,ind,jnd)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (knd,jnd,lnd,ind)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xa=v(ind1)*2.0d0-v(ind2)
  418 continue
      bb(knd,lnd)=bb(knd,lnd)+xa*xx
      bb(lnd,knd)=bb(knd,lnd)
  420 continue
  430 continue
  440 continue
  450 continue
  460 continue
  470 continue
  480 continue
      do 500 i=1,ndim
      gvec(i)=gvec(i)*conver
  500 continue
      write (iunit) ndim,(gvec(j),j=1,ndim)
      write (nf16) ndim,(gvec(j),j=1,ndim)
      do 520 i=1,ndim
      write (nf16) ndim,(g(j,i),j=1,ndim)
  520 continue
      go to 560
  540 continue
      read (nf17) dummy
      read (nf17) dummy
      read (nf13) dummy
      read (nf13) dummy
  560 continue
  600 continue
      rewind nf13
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      if(num.gt.1) go to 1000
      do 620 i=1,ntot
      do 620 j=i,ntot
      aa(i,j)=bb(i,j)*fmul(num)
      aa(j,i)=aa(i,j)
  620 continue
      do 630 i=1,nmaxa
      do 630 j=1,nmaxa
      bb(j,i)=0.0d0
  630 continue
      it=iunit
      iunit=junit
      junit=it
      if ( nxcont.ne.2) write (iwr,11)
      if ( nxcont.eq.2) write (iwr,25)
      ecor=0.0d0
      do 950 isym=1,nsym
      do 650 iorb=1,norb
      ind=indorb(iorb)
      if (idorb(ind).eq.isym) go to 660
  650 continue
      go to 950
  660 continue
      jend=nsymg(isym)
      do 940 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 940
      kend=isgnum(isym,jp)
      mx=isgind(isym,jp,1)
      mxx=ndeg(mx)
      factor=1.0d0
      if (selctc.eq.ssely.and.idega.eq.0) factor= dfloat(mxx)
      isuma=0
      do 680 k=1,kend
      if (isgind(isym,jp,k).gt.nocc) go to 680
      isuma=isuma+1
  680 continue
      read (nf16) nd,(gvec(j),j=1,nd)
      do 700 i=1,nd
      read (nf16) ne,(g(j,i),j=1,ne)
  700 continue
      write (iwr,1)
      write (iwr,3) isym,jp
      isum=kend
      do 720 j=1,nd
      do 720 i=1,isum
      aux(i,j)=g(i,j)* dabs(g(i,j))
  720 continue
      do 740 i=1,nd
      if (gvec(i).lt.-10.0d0) go to 740
      write (iwr,4) i,gvec(i),(aux(j,i),j=1,isum)
  740 continue
      do 760 j=1,isuma
      sum=0.0d0
      do 750 i=1,nd
      if (gvec(i).le.0.0d0) go to 750
      sum=sum+ dabs(aux(j,i))
  750 continue
      write (iwr,15) j,sum
  760 continue
      write (iwr,1)
      write (iwr,16)
      icw=0
      do 800 i=1,nd
      if (gvec(i).le.0.0d0) go to 800
      do 780 j=1,isuma
      if ( dabs(aux(j,i)).ge.0.008d0) go to 790
  780 continue
      go to 800
  790 write (iwr,4) i,gvec(i),(aux(j,i),j=1,isuma)
      icw=icw+1
      itemp(icw)=i
  800 continue
      write (iwr,1)
      write (iwr,6)
      write (iwr,5)
      do 900 i=1,nd
      ic=0
      do 880 j=1,nd
      if ( dabs(g(j,i)).lt.0.15d0) go to 880
      ic=ic+1
      if (j.gt.isum) go to 820
      l=isgind(isym,jp,j)
      bmvec(ic)=typ1
      mvec(ic)=l
      go to 860
  820 continue
      if (j.gt.(isum+idimb(isym,jp))) go to 840
      jx=j-isum
      l=ivecb(isym,jp,jx)
      if (l.eq.0) bmvec(ic)=typ4
      if (l.ne.0) bmvec(ic)=typ2
      mvec(ic)=l
      go to 860
  840 continue
      jx=j-isum-idimb(isym,jp)
      l=iveca(isym,jp,jx)
      if (l.eq.0) bmvec(ic)=typ4
      if (l.ne.0) bmvec(ic)=typ3
      mvec(ic)=l
  860 continue
      amvec(ic)=g(j,i)
      if (ic.lt.10) go to 880
      write (iwr,7) i,(bmvec(k),mvec(k),amvec(k),k=1,10)
      ic=0
  880 continue
      if (ic.eq.0) go to 900
      if (ic.eq.1.and. dabs(amvec(1)).ge.0.99d0) go to 900
      write (iwr,7) i,(bmvec(k),mvec(k),amvec(k),k=1,ic)
  900 continue
c     write eigenvectors on microfiche
      iend=0
  901 istart=iend+1
      iend=iend+10
      if (istart.gt.icw) go to 906
      if (iend.gt.icw) iend=icw
      do 905 j=1,nd
      if (j.gt.isum) go to 902
      l=isgind(isym,jp,j)
c     x=typ1
      go to 904
  902 continue
      if (j.gt.(isum+idimb(isym,jp))) go to 903
      jx=j-isum
      l=ivecb(isym,jp,jx)
c     if (l.eq.0) x=typ4
c     if (l.ne.0) x=typ2
      go to 904
  903 continue
      jx=j-isum-idimb(isym,jp)
      l=iveca(isym,jp,jx)
c     if (l.eq.0) x=typ4
c     if (l.ne.0) x=typ3
  904 continue
  905 continue
      go to 901
  906 continue
      if (selctb.eq.zzsel) go to 940
      do 910 i=1,nd
      gvec(i)=gvec(i)/conver
  910 continue
      call correl (ecor,factor,isum,nd,jp,
     * nmaxa,nmaxb,nmaxf,v,g,gvec,bb
_IF(ibm,vax)
     *,iwr )
_ELSE
     * ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
  940 continue
  950 continue
      rewind nf16
      if (selctb.eq.zzsel) go to 960
      ecore=ecor*conver
      write (iwr,13) ecor,ecore
  960 continue
      go to 2000
 1000 continue
      if (num.ne.2) go to 1300
      write (iwr,17)
      do 1200 isym=1,nsym
      do 1100 iorb=1,norb
      ind=indorb(iorb)
      if (idorb(ind).eq.isym) go to 1120
 1100 continue
      go to 1200
 1120 continue
      jend=nsymg(isym)
      do 1190 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 1190
      kend=isgnum(isym,jp)
      read (nf16) nd,(gvec(j),j=1,nd)
      do 1140 i=1,nd
      read (nf16) ne,(g(j,i),j=1,ne)
 1140 continue
      write (iwr,1)
      write (iwr,3) isym,jp
      isum=kend
      do 1160 j=1,nd
      do 1160 i=1,isum
      aux(i,j)=g(i,j)* dabs(g(i,j))
 1160 continue
      do 1180 i=1,nd
      if (gvec(i).lt.-10.0d0) go to 1180
      write (iwr,4) i,gvec(i),(aux(j,i),j=1,isum)
 1180 continue
 1190 continue
 1200 continue
      rewind nf16
 1300 continue
      do 1600 isym=1,nsym
      do 1400 iorb=1,norb
      ind=indorb(iorb)
      if (idorb(ind).eq.isym) go to 1420
 1400 continue
      go to 1600
 1420 continue
      if (num.le.numc) go to 1460
      jend=nsymg(isym)
      do 1450 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 1450
      read (iunit) nd,(aux(1,j),j=1,nd)
      read (junit) nd,(aux(2,j),j=1,nd)
      do 1440 i=1,nd
      if (aux(1,i).lt.0.0d0) go to 1440
      if (aux(1,i).gt.( dabs(epsil(1))*conver+30.0d0)) go to 1440
      if ( dabs(aux(1,i)-aux(2,i)).gt.0.1d0) go to 1460
c
c     0.1 ev is chosen as the criterion for convergence in the energies
c
 1440 continue
 1450 continue
      go to 1600
 1460 continue
      rewind nf18
      rewind nf19
      do 1480 i=1,ntot
      do 1480 j=i,ntot
      aa(i,j)=bb(i,j)*fmul(num)
      aa(j,i)=aa(i,j)
 1480 continue
      do 1500 i=1,nmaxa
      do 1500 j=1,nmaxa
      bb(j,i)=0.0d0
 1500 continue
      it=iunit
      iunit=junit
      junit=it
      go to 2000
 1600 continue
      go to 2100
 2000 continue
      if (numend.eq.1) go to 5000
      write (iwr,10)
      go to 2120
 2100 continue
      write (iwr,2) num
 2120 continue
      rewind nf13
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      ecor=0.0d0
      do 3000 isym=1,nsym
      do 2200 iorb=1,norb
      ind=indorb(iorb)
      if (idorb(ind).eq.isym) go to 2220
 2200 continue
      go to 3000
 2220 continue
      jend=nsymg(isym)
      do 2900 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 2900
      kend=isgnum(isym,jp)
      mx=isgind(isym,jp,1)
      mxx=ndeg(mx)
      factor=1.0d0
      if (selctc.eq.ssely.and.idega.eq.0) factor= dfloat(mxx)
      isuma=0
      do 2250 k=1,kend
      if (isgind(isym,jp,k).gt.nocc) go to 2250
      isuma=isuma+1
 2250 continue
      read (nf16) nd,(gvec(j),j=1,nd)
      do 2300 i=1,nd
      read (nf16) ne,(g(j,i),j=1,ne)
 2300 continue
      write (iwr,1)
      write (iwr,3) isym,jp
      isum=kend
      do 2320 j=1,nd
      do 2320 i=1,isum
      aux(i,j)=g(i,j)* dabs(g(i,j))
 2320 continue
      do 2340 i=1,nd
      if (gvec(i).lt.-10.0d0) go to 2340
      write (iwr,4) i,gvec(i),(aux(j,i),j=1,isum)
 2340 continue
      do 2380 j=1,isuma
      sum=0.0d0
      do 2360 i=1,nd
      if (gvec(i).le.0.0d0) go to 2360
      sum=sum+ dabs(aux(j,i))
 2360 continue
      write (iwr,15) j,sum
 2380 continue
      write (iwr,1)
      write (iwr,16)
      do 2500 i=1,nd
      if (gvec(i).le.0.0d0) go to 2500
      do 2420 j=1,isuma
      if ( dabs(aux(j,i)).ge.0.008d0) go to 2450
 2420 continue
      go to 2500
 2450 write (iwr,4) i,gvec(i),(aux(j,i),j=1,isuma)
 2500 continue
      write (iwr,1)
      write (iwr,6)
      write (iwr,5)
      do 2800 i=1,nd
      ic=0
      do 2720 j=1,nd
      if ( dabs(g(j,i)).lt.0.15d0) go to 2720
      ic=ic+1
      if (j.gt.isum) go to 2600
      l=isgind(isym,jp,j)
      bmvec(ic)=typ1
      mvec(ic)=l
      go to 2700
 2600 continue
      if (j.gt.(isum+idimb(isym,jp))) go to 2620
      jx=j-isum
      l=ivecb(isym,jp,jx)
      if (l.eq.0) bmvec(ic)=typ4
      if (l.ne.0) bmvec(ic)=typ2
      mvec(ic)=l
      go to 2700
 2620 continue
      jx=j-isum-idimb(isym,jp)
      l=iveca(isym,jp,jx)
      if (l.eq.0) bmvec(ic)=typ4
      if (l.ne.0) bmvec(ic)=typ3
      mvec(ic)=l
 2700 continue
      amvec(ic)=g(j,i)
      if (ic.lt.10) go to 2720
      write (iwr,7) i,(bmvec(k),mvec(k),amvec(k),k=1,10)
      ic=0
 2720 continue
      if (ic.eq.0) go to 2800
      if (ic.eq.1.and. dabs(amvec(1)).ge.0.99d0) go to 2800
      write (iwr,7) i,(bmvec(k),mvec(k),amvec(k),k=1,ic)
 2800 continue
      if (selctb.eq.zzsel) go to 2900
      do 2850 i=1,nd
      gvec(i)=gvec(i)/conver
 2850 continue
      call correl(ecor,factor,isum,nd,jp,
     * nmaxa,nmaxb,nmaxf,v,g,gvec,bb
_IF(ibm,vax)
     *,iwr )
_ELSE
     * ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
 2900 continue
 3000 continue
      rewind nf16
      if (selctb.eq.zzsel) go to 3100
      ecore=ecor*conver
      write (iwr,14) ecor,ecore
 3100 continue
 5000 continue
      time11=time11+cpulft(1)-tim11
      return
      end
      subroutine diadig(nmaxa,nmaxb,nmaxc,nmaxh,g,gvec,aux,iwr)
c
c test code for diagonal tda calculation - probably meaningless
c if a-diagrams are included
c
      implicit REAL  (a-h,o-z)
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
c     character*8    selct
c     character*8    typ1,typ2,typ3,typ4
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
      common/linkup/iiii,jjjj,kkkk,llll
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
      dimension g(nmaxb,nmaxb)
      dimension gvec(nmaxb),aux(nmaxc,nmaxb)
c     data selct  /'adiagram'/
c     data typ1   /'m'/
c     data typ2   /'i'/
c     data typ3   /'a'/
c     data typ4   /'c'/
c
    1 format ('1')
    3 format (1x,'energies and pole strengths in diagonal',
     +           ' approximation  for orbital   ',i5//)
    4 format (5(3x,i3,1x,f9.3,1x,f5.3))
    6 format (1x,'large components of eigenvectors in diagonal',
     +           ' approximation'//)
    7 format ('   eigenvector nr',i5,'::',10(1x,i3,':',f5.2))
c
      tim11=cpulft(1)
      time11=0.0d0
c     iunit=nf18
c     junit=nf19

      write(iwr,*)
      write(iwr,*)
      write(iwr,*)'in diadig -- THIS IS EXPERIMENTAL --- BEWARE'
      write(iwr,*)
      write(iwr,*)
      do 600 isym = 1,nsym
         write(iwr,*)'symmetry ',isym
         jend=nsymg(isym)
         do 560 jp=1,jend
            if (iuse(isym,jp).eq.1) then
               write(iwr,*)'subpart',jp
               ic = 0
               ndim=idimb(isym,jp)+idima(isym,jp)+1
               write(iwr,*)'ndimct, ndimat ',ndimct(isym),ndimat(isym)
               write(iwr,*)'dimensions ',idima(isym,jp), idimb(isym,jp),
     &              ndim
               write(iwr,*)'jdimb',jdimb(isym,jp)

               if(nmaxh.lt.idima(isym,jp))then
                  write(iwr,*)'extra column for approx poles'
                  ndim = ndim + 1
               endif
               do 100 i=1,norb
                  ind=indorb(i)
                  if (idorb(ind).eq.isym) then
                     ic = ic + 1
                     write(iwr,*)'orbital, ic ',ind,ic
                     rewind nf13
                     rewind nf16
                     rewind nf17
                     rewind nf18
                     rewind nf19
                     do 120 iloop=1,nmaxb
                        do 120 jloop=1,nmaxb
                           g(jloop,iloop)=0.0d0
 120                 continue

                     g(1,1) = -epsil(ind)
c dummy reads to reach correct symmetry                     
                     iskip = 0
                     do 701 isym2 = 1,nsym
                        jend2=nsymg(isym2)
                        do 702 jp2=1,jend2
                           if((isym2 .eq. isym .and. jp2. eq.jp) .or.
     &                        (iskip.eq.1)     )then
                              iskip = 1
                           else if (iuse(isym2,jp2).eq.1) then
                              read (nf17) dummy
                              read (nf17) dummy
                              read (nf13) dummy
                              read (nf13) dummy
                           endif
 702                    continue
 701                 continue

                     read (nf17) nd,(gvec(iloop),iloop=1,nd)
                     read (nf17) is,nd,((aux(jloop,iloop),
     &                    jloop=1,is),iloop=1,nd)

                     write(iwr,*)'is,nd,idima for a block',is,nd,
     &                    idima(isym,jp)

                     if(nd.gt.nmaxh)then
c                        write(iwr,*)'extra block'
c                        nex = nd - nmaxh 
c                        do jjj = 1,nex
c                           write(iwr,909)gvec(nmaxh + jjj),
c     &                          (aux(jjj,nmaxh + kkk),kkk=1,nex)
c                        enddo
                        nx = nmaxh + ic
                        g(ndim,ndim)=-gvec(nx)
                        xx = aux(ic,nx)
                        g(1,ndim) = xx
                        g(ndim,1) = xx
                     endif
c 909                 format(1x,8e10.4)

                     ks=idimb(isym,jp)+2
                     ke=ks+nd-1
                     nc=0
                     do 180 k=ks,ke
                        nc=nc+1
                        g(k,k)=-gvec(nc)
                        xx=aux(ic,nc)
                        g(1,k)=xx
                        g(k,1)=xx
 180                 continue

                     read (nf13) nd,(gvec(ir),ir=1,nd)
                     read (nf13) is,nd,((aux(jr,ir),jr=1,is),ir=1,nd)

                     write(iwr,*)'is,nd,idimb for b block',is,nd,
     &                    idimb(isym,jp)
                     ls=2
                     le=nd+1
                     nc=0
                     do 220 l=ls,le
                        nc=nc+1
                        g(l,l)=-gvec(nc)
                        xx = aux(ic,nc)
                        g(1,l) = xx
                        g(l,1) = xx
 220                 continue

                     time11=time11+cpulft(1)-tim11
                     call sdiag2(nmaxb,ndim,g,gvec,g)
                     tim11=cpulft(1)

                     write (iwr,3) i
                     do 200 j=1,ndim
                        gvec(j)=gvec(j)*conver
                        aux(1,j)=g(1,j)**2
 200                 continue
                     iend=0
 210                 istart=iend+1
                     iend=iend+5
                     if (istart.gt.ndim) go to 230
                     if (iend.gt.ndim) iend=ndim
                     if (gvec(iend).lt.-10.0e0) go to 210
                     write (iwr,4) (j,gvec(j),aux(1,j),j=istart,iend)
                     go to 210
 230                 continue
                     write (iwr,1)
                     write (iwr,6)
                     do 450 l=1,ndim
                        icount=0
                        do 400 j=1,ndim
                           if ( abs(g(j,l)).lt.0.15e0) go to 400
                           icount=icount+1
                           mvec(icount)=j
                           amvec(icount)=g(j,l)
                           if (icount.lt.10) go to 400
                           write(iwr,7) l,(mvec(k),amvec(k),k=1,10)
                           icount=0
 400                    continue
                        if (icount.eq.0) go to 450
                        if (icount.eq.1.and. 
     &                       abs(amvec(1)).ge.0.99e0) go to 450
                        write(iwr,7) l,(mvec(k),amvec(k),k=1,icount)
 450                 continue
                  endif
 100           continue
            endif
 560     continue
 600  continue
c
      return
      end
      subroutine compol(nmaxbg,nmaxc,nmax,iv,jv,veca,vecb)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iv,jv
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(3),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
      dimension iv(nmaxbg),jv(nmaxbg),veca(nmax),vecb(nmax,nmaxc)
c     this subroutine serves to combine the affinity poles and
c     amplitudes from file 17 and the ionization poles and
c     amplitudes from file 13. these two data sets are combined
c     for each symmetry species and subelement into one data set
c     which is written onto file 15.
c     only the configuration identification vectors are kept
c     separately.
c
c
      rewind nf13
      rewind nf15
      rewind nf17
      do 100 isym=1,nsym
      do 20 i=1,norb
      iorb=indorb(i)
      if (idorb(iorb).eq.isym) go to 30
   20 continue
      go to 100
   30 continue
      jend=nsymg(isym)
      do 50 j=1,jend
      read (nf17) ix,iy,iz,it,(veca(i),i=1,iz),
     *((vecb(i,k),i=1,iz),k=1,it),(iv(i),i=1,iz)
      read (nf13) ix,iy,jz,it,(veca(iz+i),i=1,jz),
     *((vecb(iz+i,k),i=1,jz),k=1,it),(jv(i),i=1,jz)
      kz=iz+jz
      write (nf15) ix,iy,iz,jz,kz,it,(veca(i),i=1,kz),
     * ((vecb(i,k),i=1,kz),k=1,it),(iv(i),i=1,iz),(jv(i),i=1,jz)
   50 continue
  100 continue
      rewind nf13
      rewind nf15
      rewind nf17
c
      return
      end
      subroutine correl (e,factor,isum,nd,jp,
     *nmaxa,nmaxb,nmaxf,v,g,gvec,bb
_IF(ibm,vax)
     *,iwr )
_ELSE
     * ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
c
      implicit REAL  (a-h,o-z)
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
_IFN(ibm,vax)
      common/linkup/iiii,jjjj,kkkk,llll
_ENDIF
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
      dimension v(nmaxf),g(nmaxb,nmaxb),gvec(nmaxb),bb(nmaxa,nmaxa)
      dimension xorb(mxtda2)
      equivalence (xorb(1),sigdia(1))
_IFN(ibm,vax)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
_ENDIF
c
    1 format (//1x,
     +'orbital contribution to the correlation energy from ip '/
     +1x,'and residues of the green function part for symmetry  ',i3,
     *' subpart ',i3/(10(2x,f10.6)))
    2 format(1x,
     +'orbital contribution to the correlation energy from green '/
     +'function times orbital energy part for symmetry ',i3,
     *' subpart ',i3/(10(2x,f10.6)))
    3 format (1x,
     +'orbital contribution to the correlation energy from the '/
     +1x,'diagonal part of the interaction times'/1x,
     +'green function part for symmetry ',i3,' subpart ',i3/
     + (10(2x,f10.6)))
    4 format(1x,
     +'total part of this (diagonal) correlation contribution ',f10.6)
    5 format(1x,
     +'total part of this correlation including the nondiagonal part '
     + ,f10.6)
    6 format(' total part of this correlation energy ',f10.6)
    7 format(1x,
     +'correlation energy contribution of the core interaction part ',
     + f10.6)
    8 format  (/1x,
     +'the total correlation energy after adding the contribution ',
     + 'from '/'symmetry type ',i3,' subpart ',i3,' is ',f10.6)
c
      ex=0.0d0
      do 20 i=1,nmaxa
      xorb(i)=0.0d0
   20 continue
      nc=ixxx(isym,jp)
      do 40 i=1,isum
      x=0.0d0
      do 30 j=nc,nd
      x=x-gvec(j)*(g(i,j)**2)
   30 continue
      x=x*factor
      xorb(i)=x
      ex=ex+x
   40 continue
      write (iwr,1) isym,jp,(xorb(i),i=1,isum)
      write (iwr,6) ex
      e=e+ex
      ex=0.0d0
      do 50 i=1,nmaxa
      do 50 j=1,nmaxa
      bb(j,i)=0.0d0
   50 continue
      do 80 k=1,isum
      i=isgind(isym,jp,k)
      do 70 l=k,isum
      j=isgind(isym,jp,l)
      x=0.0d0
      do 60 nn=nc,nd
      x=x+g(k,nn)*g(l,nn)
   60 continue
      bb(i,j)=x
      bb(j,i)=x
      if (i.eq.j.and.i.le.nocc) bb(i,i)=bb(i,i)-1.0d0
   70 continue
   80 continue
      do 85 i=1,nmaxa
      xorb(i)=0.0d0
   85 continue
      do 90 k=1,isum
      i=isgind(isym,jp,k)
      x=bb(i,i)*epsil(i)
      x=x*factor
      xorb(k)=x
      ex=ex+x
      if (i.gt.nocc) go to 90
      ex=ex-epsil(i)*factor
      xorb(k)=xorb(k)-epsil(i)*factor
   90 continue
      write (iwr,2) isym,jp,(xorb(i),i=1,isum)
      write (iwr,6) ex
      e=e+ex
      ex=0.0d0
      do 95 i=1,nmaxa
      xorb(i)=0.0d0
   95 continue
      do 120 l=1,isum
      i=isgind(isym,jp,l)
      do 110 m=l,isum
      j=isgind(isym,jp,m)
      x=0.0d0
      do 100 k=1,nocc
_IF(ibm,vax)
      ind = indvor (i,j,k,k)
      jnd = indvor (i,k,j,k)
_ELSE
      call  indvor (i,j,k,k)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (i,k,j,k)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-2.0d0*v(ind)+v(jnd)
  100 continue
      if (i.ne.j) x=x+x
      x=x*factor
      ex=ex+x*bb(i,j)
      if (i.eq.j) xorb(l)=xorb(l)+x*bb(i,i)
  110 continue
  120 continue
      write (iwr,3) isym,jp,(xorb(i),i=1,isum)
      x=0.0d0
      do 125 i=1,isum
      x=x+xorb(i)
  125 continue
      write (iwr,4) x
      write (iwr,5) ex
      e=e+ex
      ex=0.0d0
      rewind nf14
      read (nf14) nx,((g(j,i),j=1,nx),i=1,nx)
      rewind nf14
      do 130 i=1,ntot
      do 130 j=1,ntot
      ex=ex+g(j,i)*bb(j,i)
  130 continue
      ex=ex*factor
      write (iwr,7) ex
      e=e+ex
      ex=0.0d0
      write (iwr,8) isym,jp,e
      do 300 i=1,nmaxa
      do 300 j=1,nmaxa
      bb(j,i)=0.0d0
  300 continue
      return
      end
      subroutine datgaa (jkeep,jkeepx,
     * nmaxb,nmaxc,nmaxg,nmaxh,nmaxi,
     *iveca,g,gvec,aux,veca,vecb,axx,iwr)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iveca
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
      dimension iveca(8,nmaxi,*)
      dimension veca(nmaxg),vecb(nmaxg,nmaxc),vecc(mxtda4)
      dimension g(nmaxb,nmaxb),gvec(nmaxb),aux(nmaxc,nmaxb)
      dimension axx(nmaxb,6)
c
c     this subroutine analyzes the eigenvalues and eigenvectors of the
c     gammaa matrix for a with respect to the subparts of the given
c     symmetry and sorts them into different files which are written on
c     file 17. simultaneously the new dimensions are calculated as well
c     as a new index vector which connects the previous (i,ia,ib) index
c     with the symmetry subpart index.
c
    1 format (/1x,104('-')/)
    6 format (//' list of poles of the self-energy'//(10(1x,f11.5)))
    8 format (' poles of affinity block of self-energy for symmetry  '
     *,i3,' and subpart  ',i3//)
    9 format (10(1x,f11.5))
   10 format (/' residues of affinity block of the self-energy',
     +' for orbital ',i3/)
c
      tim77=cpulft(1)
      jend=nsymg(isym)
      do 50 i=1,ndim
      sigdia(i)=gvec(i)*conver
   50 continue
      write (iwr,6) (sigdia(i),i=1,ndim)
      write (iwr,1)
      crit=3.0d-04
      do 60 i=1,nmaxb
      sigdia(i)=0.0d0
      do 60 j=1,nmaxb
      g(j,i)=0.0d0
   60 continue
c
c     loop over all subparts of symmetry isym.
c     analyze only the solutions which are fully taken into account in
c     the dyson equation and not those which are used only to form
c     centroids.
c
      icont=0
      do 200 j=1,jend
      ic=0
      kend=isgnum(isym,j)
c     loop over all solutions
      do 160 i=1,ndim
      sum=0.0d0
      do 100 k=1,kend
      ind=jsgjnd(isym,j,k)
      sum=sum+ dabs(aux(ind,i))
  100 continue
      if (sum.lt.crit) go to 160
      ic=ic+1
      if (ic.gt.nmaxh) go to 180
      iveca(isym,j,ic)=i
  160 continue
      idima(isym,j)=ic
      go to 200
  180 idima(isym,j)=nmaxh
      icont=1
  200 continue
      do 210 j=1,nmaxb
      vecc(j)=0.0d0
  210 continue
c
c     calculate centroids in combination if necessary with the diagonal
c     data.
c
      if (ndim.ne.ndimat(isym)) go to 320
c
c     go to 320 if centroids from the diagonal part are to be included
c
      if (icont.eq.0) go to 700
c
c     go to 700 if there will not be any centroids
c     continue here if centroids are only included from the nondiagonal
c     calculation
c
      read (nf15) dummy
      do 300 i=1,jkeep
      p=0.0d0
      e=0.0d0
      do 240 k=1,ndim
      do 220 jp=1,jend
      ix=idima(isym,jp)
      do 220 kp=1,ix
      if (k.eq.iveca(isym,jp,kp)) go to 240
  220 continue
      x=aux(i,k)**2
      p=p+x
      e=e+gvec(k)*x
  240 continue
      e=e/p
      p= dsqrt(p)
      sigdia(i)=e
      vecc(i)=p
  300 continue
      go to 500
  320 continue
      read(nf15)jn,kn,in,ln,(veca(i),i=1,in),((vecb(i,j),i=1,in),j=1,ln)
      do 450 i=1,jkeep
      p=0.0d0
      e=0.0d0
      if (i.gt.jkeepx) go to 380
      do 360 k=1,ndim
      do 340 jp=1,jend
      ix=idima(isym,jp)
      do 340 kp=1,ix
      if (k.eq.iveca(isym,jp,kp)) go to 360
  340 continue
      x=aux(i,k)**2
      p=p+x
      e=e+gvec(k)*x
  360 continue
  380 continue
      do 400 k=1,in
      x=vecb(k,i)**2
      p=p+x
      e=e-veca(k)*x
  400 continue
      e=e/p
      p= dsqrt(p)
      sigdia(i)=e
      vecc(i)=p
  450 continue
  500 continue
      do 520 i=1,nmaxb
      do 520 j=1,nmaxb
      g(j,i)=0.0d0
  520 continue
      do 600 j=1,jend
      jx=(j-1)*nmaxc
      ic=0
      kend=isgnum(isym,j)
      do 580 i=1,ndim
      sum=0.0d0
      do 540 k=1,kend
      ind=jsgjnd(isym,j,k)
      sum=sum+ dabs(aux(ind,i))
  540 continue
      if (sum.lt.crit) go to 580
      ic=ic+1
      if (ic.gt.nmaxh) go to 600
      do 560 l=1,kend
      ind=jsgjnd(isym,j,l)
      g(jx+l,ic)=aux(ind,i)
  560 continue
      axx(ic,j)=gvec(i)
  580 continue
  600 continue
      do 660 j=1,jend
      kend=isgnum(isym,j)
      jx=(j-1)*nmaxc
      kc=0
      ic=idima(isym,j)
      do 640 i=1,jkeep
      do 620 k=1,kend
      if (jsgjnd(isym,j,k).ne.i) go to 620
      kc=kc+1
      axx(ic+kc,j)=sigdia(i)
      g(jx+kc,ic+kc)=vecc(i)
      iveca(isym,j,ic+kc)=0
  620 continue
  640 continue
      idima(isym,j)=idima(isym,j)+kc
  660 continue
      go to 820
  700 continue
      read (nf15) dummy
      do 720 i=1,nmaxb
      do 720 j=1,nmaxb
      g(j,i)=0.0d0
  720 continue
      do 800 j=1,jend
      jx=(j-1)*nmaxc
      ic=0
      kend=isgnum(isym,j)
      do 780 i=1,ndim
      sum=0.0d0
      do 740 k=1,kend
      ind=jsgjnd(isym,j,k)
      sum=sum+ dabs(aux(ind,i))
  740 continue
      if (sum.lt.crit) go to 780
      ic=ic+1
      if (ic.gt.nmaxh) go to 800
      do 760 l=1,kend
      ind=jsgjnd(isym,j,l)
      g(jx+l,ic)=aux(ind,i)
  760 continue
      axx(ic,j)=gvec(i)
  780 continue
  800 continue
  820 continue
c
c     write poles of the affinity block of the self-energy and the
c     associated amplitudes of the residues of the self-energy on
c     file 17, two data sets for each subpart of the given symmetry.
c
      do 840 j=1,jend
      ic=isgnum(isym,j)
      ix=idima(isym,j)
      jx=(j-1)*nmaxc
      write (nf17) ix,(axx(i,j),i=1,ix)
      write (nf17) ic,ix,((g(jx+i,k),i=1,ic),k=1,ix)
  840 continue
      do 940 j=1,jend
      kend=isgnum(isym,j)
      ix=idima(isym,j)
      write (iwr,8) isym,j
      jx=(j-1)*nmaxc
      do 910 i=1,ix
      gvec(i)=axx(i,j)*conver
  910 continue
      write (iwr,9) (gvec(i),i=1,ix)
      do 930 k=1,kend
      ind=isgind(isym,j,k)
      write (iwr,10) ind
      do 920 l=1,ix
      aux(1,l)=cons*(g(jx+k,l)**2)
  920 continue
      write (iwr,9) (aux(1,l),l=1,ix)
  930 continue
  940 continue
      do 990 i=1,nmaxb
      gvec(i)=0.0d0
      sigdia(i)=0.0d0
      do 970 j=1,nmaxc
      aux(j,i)=0.0d0
  970 continue
      do 980 k=1,nmaxb
      g(k,i)=0.0d0
  980 continue
  990 continue
      time7(isym)=cpulft(1)-tim77
      return
      end
      subroutine datgab(nmaxb,nmaxc,nmaxi,ivecb,
     * g,gvec,aux,axx,iwr)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 ivecb
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
c
      dimension g(nmaxb,nmaxb),gvec(nmaxb),aux(nmaxc,nmaxb)
      dimension axx(nmaxb,nmaxi),ivecb(8,nmaxi,*)
c
c     this subroutine analyzes the eigenvalues and eigenvectors of the
c     gammaa matrix for b with respect to the subparts of the given
c     symmetry and sorts them into different files which are written on
c     file 13. simultaneously the new dimensions are calculated as well
c     as a new index vector which connects the previous (i,j,ia) index
c     with the symmetry subpart index.
c
    1 format (/1x,104('-'))
    6 format (//' list of poles of the self-energy'//(10(1x,f11.5)))
    8 format (1x,
     +'poles of ionization block of self-energy for symmetry',i3,
     + ' and subpart  ',i3//)
    9 format (10(1x,f11.5))
   10 format (/1x,
     + 'residues of ionization block of the self-energy for orbital ',
     + i3/)
c
      tim10=cpulft(1)
      jend=nsymg(isym)
c
      do 50 i=1,ndim
      sigdia(i)=gvec(i)*conver
   50 continue
      write (iwr,6) (sigdia(i),i=1,ndim)
      write (iwr,1)
      crit=3.0d-04
      do 60 i=1,nmaxb
      sigdia(i)=0.0d0
      do 60 j=1,nmaxb
      g(j,i)=0.0d0
   60 continue
c
c     loop over all subparts of symmetry isym
c
      do 200 j=1,jend
      jx=(j-1)*nmaxc
      ic=0
      kend=isgnum(isym,j)
c     loop over all solutions
      do 160 i=1,ndim
      sum=0.0d0
      do 100 k=1,kend
      ind=jsgjnd(isym,j,k)
      sum=sum+ dabs(aux(ind,i))
  100 continue
      if (sum.lt.crit) go to 160
      ic=ic+1
      ivecb(isym,j,ic)=i
      do 120 l=1,kend
      ind=jsgjnd(isym,j,l)
      g(jx+l,ic)=aux(ind,i)
  120 continue
      axx(ic,j)=gvec(i)
  160 continue
      jdimb(isym,j)=ic
      idimb(isym,j)=ic
  200 continue
c
c     add centroids if there any to the data sets of the subparts j of
c     symmetry isym
c
      if (ndim.eq.ndimbt(isym)) go to 280
      read (nf16) ic,jc,(aux(1,i),i=1,ic),(aux(2,j),j=1,ic)
      do 260 j=1,jend
      kend=isgnum(isym,j)
      jx=(j-1)*nmaxc
      kc=0
      ix=idimb(isym,j)
      do 240 i=1,ic
      do 220 k=1,kend
      if (jsgjnd(isym,j,k).ne.i) go to 220
      kc=kc+1
      axx(ix+kc,j)=aux(1,i)
      g(jx+kc,ix+kc)=aux(2,i)
      ivecb(isym,j,ix+kc)=0
  220 continue
  240 continue
      idimb(isym,j)=idimb(isym,j)+kc
  260 continue
      go to 300
  280 read (nf16) dummy
  300 continue
c
c     write poles of the ionization block of the self-energy and the
c     associated amplitudes of the residues of the self-energy on
c     file 13, two data sets for each subpart of the given symmetry.
c
      do 400 j=1,jend
      ic=isgnum(isym,j)
      ix=idimb(isym,j)
      jx=(j-1)*nmaxc
      write (nf13) ix,(axx(i,j),i=1,ix)
      write (nf13) ic,ix,((g(jx+i,k),i=1,ic),k=1,ix)
  400 continue
      do 600 j=1,jend
      kend=isgnum(isym,j)
      ix=idimb(isym,j)
      write (iwr,8) isym,j
      jx=(j-1)*nmaxc
      do 510 i=1,ix
      gvec(i)=axx(i,j)*conver
  510 continue
      write (iwr,9) (gvec(i),i=1,ix)
      do 540 k=1,kend
      ind=isgind(isym,j,k)
      write (iwr,10) ind
      do 530 l=1,ix
      aux(1,l)=cons*(g(jx+k,l)**2)
  530 continue
      write (iwr,9) (aux(1,l),l=1,ix)
  540 continue
  600 continue
      do 630 i=1,nmaxb
      gvec(i)=0.0d0
      sigdia(i)=0.0d0
      do 610 j=1,nmaxc
      aux(j,i)=0.0d0
  610 continue
      do 620 k=1,nmaxb
      g(k,i)=0.0d0
  620 continue
  630 continue
      time10(isym)=cpulft(1)-tim10
      return
      end
      subroutine diagax (kdim,ldim,jkeep,
     * nmaxa,nmaxc,nmaxg,nmaxgj,nmaxj,
     * iv,jv,kv,iconve,itv,jtv,ktv,ga,sta,aa,iwr)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iv,jv,kv,iconve
_IFN1(c)      integer*2 itv,jtv,ktv
      character*1 xtyp,txtyp,ctemp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
c
      dimension txtyp(mxtda4)
      dimension ga(nmaxg),sta(nmaxg,nmaxc),aa(nmaxa,nmaxa)
      dimension iv(nmaxgj),jv(nmaxgj),kv(nmaxgj),iconve(nmaxgj)
      dimension itv(nmaxj),jtv(nmaxj),ktv(nmaxj)
c
c     npola gives the number of poles taken into account for the
c           a block
c     npolb gives the number of poles taken into account for the
c           b block
c     npolt gives the number of poles taken into account for both
c           together
c
    2 format(/' the number of diagonally approximated poles',
     +1x,'for the a-block for symmetry ',i5,' is',i8)
    3 format(/' index identification of these poles ordered ',
     +'according to energy'/)
    4 format (/6(2x,i4,' = (',i2,',',i2,',',i2,')',a1))
      dummy=0.0d0
c
c
c     data handling for the case of pole search solution to the
c     dyson equation
c
c
      if (ntotx.eq.ntot.and.ldim.le.nmaxj) go to 800
c
c     store indices of configurations of poles taken into account in a
c     nondiagonal form in one set of vectors and the indices of the
c     configurations taken into account diagonally only in another
c     one. write both on disk.
c
c     change sign of poles taken into account diagonally to agree with
c     the sign convention in the routine gammaa
c
      ic=0
      jc=0
      do 560 i=1,kdim
      if (kv(i).gt.ntotx) go to 570
      ic=ic+1
      jc=jc+1
      if (jv(i).ne.kv(i)) ic=ic+1
      if (ic.gt.(nmaxj-2)) go to 570
  560 continue
  570 continue
      kc=ic-jc
      do 580 k=1,jc
      itv(k)=iv(k)
      jtv(k)=jv(k)
      ktv(k)=kv(k)
      txtyp(k)=xtyp(k)
  580 continue
      do 590 k=1,kc
      ix=jc+k
      jx=kdim+k
      itv(ix)=iv(jx)
      jtv(ix)=jv(jx)
      ktv(ix)=kv(jx)
      txtyp(ix)=xtyp(jx)
  590 continue
      istart=jc+1
      in=kdim-jc
      jstart=kdim+kc+1
      jn=ldim-jstart+1
      kn=in+jn
      do 610 i=istart,kdim
      ix=i-jc
      iv(ix)=iv(i)
      jv(ix)=jv(i)
      kv(ix)=kv(i)
      xtyp(ix)=xtyp(i)
      ga(ix)=-ga(i)
      do 600 j=1,jkeep
      sta(ix,j)=sta(i,j)
  600 continue
  610 continue
      do 630 j=jstart,ldim
      jj=j-kdim-kc+in
      iv(jj)=iv(j)
      jv(jj)=jv(j)
      kv(jj)=kv(j)
      xtyp(jj)=xtyp(j)
      ga(jj)=-ga(j)
      do 620 l=1,jkeep
      sta(jj,l)=sta(j,l)
  620 continue
  630 continue
c
c     order diagonal poles, strengths and indices
c
      iend=kn-1
      do 660 i=1,iend
      js=i+1
      do 650 j=js,kn
      if (ga(i).le.ga(j)) go to 650
      it=iv(i)
      iv(i)=iv(j)
      iv(j)=it
      it=jv(i)
      jv(i)=jv(j)
      jv(j)=it
      it=kv(i)
      kv(i)=kv(j)
      kv(j)=it
      ctemp=xtyp(i)
      xtyp(i)=xtyp(j)
      xtyp(j)=ctemp
      temp=ga(i)
      ga(i)=ga(j)
      ga(j)=temp
      do 640 k=1,jkeep
      temp=sta(i,k)
      sta(i,k)=sta(j,k)
      sta(j,k)=temp
  640 continue
  650 continue
  660 continue
c
c     calculation of centroids if there are any in the preset limits
c
      if (ga(kn).le.-alimit) go to 740
      do 680 k=1,kn
      kkp=k
      if (ga(k).le.-alimit) go to 680
      go to 690
  680 continue
  690 continue
      if ((kn-kkp+1).le.jkeep) go to 740
      do 710 i=1,jkeep
      p=0.0d0
      e=0.0d0
      do 700 k=kkp,kn
      x=sta(k,i)**2
      p=p+x
      e=e+ga(k)*x
  700 continue
      e=e/p
      p= dsqrt(p)
      sigdia(i)=e
      aa(1,i)=p
  710 continue
      do 720 i=kkp,kn
      ga(i)=0.0d0
      do 720 j=1,jkeep
      sta(i,j)=0.0d0
  720 continue
      do 730 i=1,jkeep
      ix=kkp-1+i
      ga(ix)=sigdia(i)
      sta(ix,i)=aa(1,i)
  730 continue
      kn=kkp-1+jkeep
  740 continue
      do 750 i=1,kn
      ix=ic+i
      iconve(ix)=ix
  750 continue
      npola(isym,1)=kn
      do 760 i=1,ic
      iconve(i)=i
  760 continue
c
c     structure of data set on unit 15 for block a in the case of
c     pole search:
c
c             first record is identical to the case in diaga routine
c             second record:
c                    kn mumber of diagonal poles (including centroids)
c                    jkeep number of orbitals of symmetry isym
c                    ga,sta energies (with sign reversed) and
c                           strengths of diagonal poles
c
      write (nf15) jc,kc,ic,(itv(i),i=1,ic),(jtv(i),i=1,ic),
     *(ktv(i),i=1,ic),(txtyp(i),i=1,ic),(iconve(i),i=1,ic)
c
      write(nf15)kn,jkeep,(ga(i),i=1,kn),((sta(i,j),i=1,kn),j=1,jkeep)
c
c     list all configurations of the diagonally approximated poles
c
      write (iwr,2) isym,kn
      write (iwr,3)
      write (iwr,4) (iconve(i+ic),iv(i),jv(i),kv(i),xtyp(i),i=1,kn)
      go to 1000
  800 continue
      npola(isym,1)=0
c
c     case all poles are taken into account in the diagonalization
c
      do 810 i=1,ldim
      iconve(i)=i
  810 continue
      kc=ldim-kdim
      write (nf15) kdim,kc,ldim,(iv(i),i=1,ldim),(jv(i),i=1,ldim),
     *(kv(i),i=1,ldim),(xtyp(i),i=1,ldim),(iconve(i),i=1,ldim)
      write (nf15) (dummy,i=1,10)
c
 1000 continue
c
      return
      end
      subroutine diagb(
     * nmaxa,nmaxb,nmaxc,nmaxe,nmaxg,nmaxgj,nmaxj,
     * iv,jv,kv,iconve,ivecc,v,ga,sta,aa,
     * indvec,jndvec,kndvec,lndvec,iwr)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iv,jv,kv,iconve,ivecc
c
      character*1  styp,ttyp,ctemp
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
_IFN(ibm,vax)
      common/linkup/iiii,jjjj,kkkk,llll
_ENDIF
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
      dimension iv(nmaxgj),jv(nmaxgj),kv(nmaxgj),iconve(nmaxgj)
      dimension ivecc(*)
      dimension v(nmaxe),ga(nmaxg),sta(nmaxg,nmaxc),aa(nmaxa,nmaxa)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
c     calculation of diagonal matrix elements of b-block, reordering
c     according to magnitude, checking for maximum indeces,
c     determination which 2h1p configurations are included in the
c     nondiagonal b-block, calculation of centroids
c
      data styp  /'s'/
      data ttyp  /'t'/
c
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
      tim4=cpulft(1)
            dummy=0.0d0
c
      do 30 i=1,nmaxg
      ga(i)=0.0d0
      iv(i)=0
      jv(i)=0
      kv(i)=0
      do 20 j=1,nmaxc
      sta(i,j)=0.0d0
   20 continue
   30 continue
      ic=0
      icc=ndimb(isym)
      do 100 ia=nstart,ntot
      iaa=idorb(ia)
      iam=imul(isym,iaa)
      epsa=epsil(ia)
      do 90 i=1,nocc
      ii=idorb(i)
      epsi=epsil(i)
      do 80 j=i,nocc
      jj=idorb(j)
      if (iam.ne.imul(ii,jj)) go to 80
      epsj=epsil(j)
      ic=ic+1
      iv(ic)=ia
      jv(ic)=i
      kv(ic)=j
      xtyp(ic)=styp
      if (i.eq.j) go to 40
      icc=icc+1
      iv(icc)=ia
      jv(icc)=i
      kv(icc)=j
      xtyp(icc)=ttyp
   40 continue
      if (i.ne.j) go to 60
      ind1= indvx (i,i,i,i)
      ind2= indvx (ia,ia,i,i)
      ind3= indvx (ia,i,ia,i)
      ga(ic)=epsa-epsi-epsj-(-v(ind1)+2.0d0*v(ind2)-v(ind3))
      jc=0
      do 50 ip=1,ntot
      if (idorb(ip).ne.isym) go to 50
      jc=jc+1
_IF(ibm,vax)
      ind = indvor (ia,i,ip,i)
_ELSE
      call indvor (ia,i,ip,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xa=v(ind)
      xc=0.0d0
      del=-2.0d0*epsil(i)
      do 700 ib=nstart,ntot
      ind = indvx (ib,i,ib,i)
_IF(ibm,vax)
      jnd = indvor (ia,ib,ib,ip)
_ELSE
      call  indvor (ia,ib,ib,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc-(v(ind)*v(jnd))/(2.0d0*epsil(ib)+del)
  700 continue
      xd=0.0d0
      if (nvirt.eq.1) go to 730
      ibend=ntot-1
      do 720 ib=nstart,ibend
      kcs=ib+1
      ibb=idorb(ib)
      eps=del+epsil(ib)
      do 710 kc=kcs,ntot
      kcc=idorb(kc)
      if(ibb.ne.kcc) go to 710
      ind = indvx (kc,i,ib,i)
_IF(ibm,vax)
      jnd = indvor (kc,ia,ib,ip)
      knd = indvor (kc,ip,ia,ib)
_ELSE
      call  indvor (kc,ia,ib,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (kc,ip,ia,ib)
      knd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xd=xd-(v(ind)*(v(jnd)+v(knd)))/(epsil(kc)+eps)
  710 continue
  720 continue
  730 continue
      del=epsil(ia)-epsil(i)
      xe=0.0d0
      ibm=imul(iaa,ii)
      do 750 k=1,nocc
      kk=idorb(k)
      eps=del-epsil(k)
      do 740 ib=nstart,ntot
      ibb=idorb(ib)
      if (ibm.ne.imul(kk,ibb)) go to 740
_IF(ibm,vax)
      ind = indvor (ia,i,ib,k)
      jnd = indvor (ib,k,ip,i)
      knd = indvor (ib,ip,k,i)
      lnd = indvor (ia,k,ib,i)
_ELSE
      call indvor (ia,i,ib,k)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,k,ip,i)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,ip,k,i)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,k,ib,i)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xe=xe+((x2*(-2.0d0*x1+x4)+x3*(x1+x4)))/(eps+epsil(ib))
  740 continue
  750 continue
      sta(ic,jc)=xa+xc+xd+xe
   50 continue
      go to 80
   60 continue
      ind1= indvx (j,j,i,i)
      ind2= indvx (j,i,j,i)
      ind3= indvx (ia,ia,j,j)
      ind4= indvx (ia,j,ia,j)
      ind5= indvx (ia,ia,i,i)
      ind = indvx (ia,i,ia,i)
      ga(ic)=epsa-epsi-epsj-(-v(ind1)-v(ind2)+v(ind3)-0.5d0*(v(ind4)+
     *       v(ind))+v(ind5))
      ga(icc)=epsa-epsi-epsj-(-v(ind1)+v(ind2)+v(ind3)-1.5d0*(v(ind4)+
     *        v(ind))+v(ind5))
      jc=0
      do 70 ip=1,ntot
      if (idorb(ip).ne.isym) go to 70
      jc=jc+1
_IF(ibm,vax)
      ind1= indvor (ia,j,ip,i)
      ind2= indvor (ia,i,ip,j)
_ELSE
      call indvor (ia,j,ip,i)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,i,ip,j)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xa=v(ind1)
      xb=v(ind2)
      xc=0.0d0
      del=-epsil(i)-epsil(j)
      if(ii.ne.jj) go to 610
      do 600 ib=nstart,ntot
      ind = indvx (ib,j,ib,i)
_IF(ibm,vax)
      jnd = indvor (ib,ia,ib,ip)
_ELSE
      call  indvor (ib,ia,ib,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc-(v(ind)*v(jnd))/(2.0d0*epsil(ib)+del)
  600 continue
      xc=xc+xc
  610 continue
      xd=0.0d0
      xe=0.0d0
      if (nvirt.eq.1) go to 650
      ibend=ntot-1
      do 640 ib=nstart,ibend
      kcs=ib+1
      ibb=idorb(ib)
      eps=epsil(ib)+del
      do 630 kc=kcs,ntot
      kcc=idorb(kc)
      if (iam.ne.imul(ibb,kcc)) go to 630
      ind = indvx (kc,j,ib,i)
      jnd = indvx (kc,i,ib,j)
_IF(ibm,vax)
      knd = indvor (kc,ia,ib,ip)
      lnd = indvor (kc,ip,ia,ib)
_ELSE
      call indvor (kc,ia,ib,ip)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (kc,ip,ia,ib)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xx=1.0d0/(eps+epsil(kc))
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xd=xd-xx*((x1+x2)*(x3+x4))
      xe=xe-xx*((x1-x2)*(x3-x4))
  630 continue
  640 continue
  650 continue
      del=epsil(ia)-epsil(j)
      xf=0.0d0
      xg=0.0d0
      ibm=imul(iaa,jj)
      do 670 k=1,nocc
      kk=idorb(k)
      eps=del-epsil(k)
      do 660 ib=nstart,ntot
      ibb=idorb(ib)
      if (imul(kk,ibb).ne.ibm) go to 660
_IF(ibm,vax)
      ind = indvor (ia,j,ib,k)
      jnd = indvor (ib,k,ip,i)
      knd = indvor (ib,ip,k,i)
      lnd = indvor (ia,k,ib,j)
_ELSE
      call indvor (ia,j,ib,k)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,k,ip,i)
      jnd = indvx(iiii,jjjj,kkkk,llll)
       call indvor (ib,ip,k,i)
       knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,k,ib,j)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xx=1.0d0/(eps+epsil(ib))
      xf=xf+xx*(x2*(-2.0d0*x1+x4)+x3*(x1+x4))
      xg=xg+xx*(x2*(-2.0d0*x1+x4)+x3*(x1-x4))
  660 continue
  670 continue
      del=epsil(ia)-epsil(i)
      xh=0.0d0
      xi=0.0d0
      ibm=imul(iaa,ii)
      do 690 k=1,nocc
      kk=idorb(k)
      eps=del-epsil(k)
      do 680 ib=nstart,ntot
      ibb=idorb(ib)
      if (imul(kk,ibb).ne.ibm) go to 680
_IF(ibm,vax)
      ind = indvor (ia,i,ib,k)
      jnd = indvor (ib,k,ip,j)
      knd = indvor (ib,ip,k,j)
      lnd = indvor (ia,k,ib,i)
_ELSE
      call indvor (ia,i,ib,k)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,k,ip,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,ip,k,j)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,k,ib,i)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xx=1.0d0/(eps+epsil(ib))
      xh=xh+xx*(x2*(-2.0d0*x1+x4)+x3*(x1+x4))
      xi=xi+xx*(x2*(2.0d0*x1-x4)+x3*(-x1+x4))
  680 continue
  690 continue
      sta(ic,jc)=fac*(xa+xb+xc+xd+xf+xh)
      sta(icc,jc)=fad*(xa-xb+xe+xg+xi)
   70 continue
      jkeep=jc
   80 continue
   90 continue
  100 continue
      kdim=ndimb(isym)
      ldim=ndimbt(isym)
c
c     ordering of diagonal elements
c
      iend=kdim-1
      do 130 i=1,iend
      js=i+1
      do 120 j=js,kdim
      if (ga(i).le.ga(j)) go to 120
      temp=ga(i)
      ga(i)=ga(j)
      ga(j)=temp
      it=iv(i)
      iv(i)=iv(j)
      iv(j)=it
      it=jv(i)
      jv(i)=jv(j)
      jv(j)=it
      it=kv(i)
      kv(i)=kv(j)
      kv(j)=it
      ctemp=xtyp(i)
      xtyp(i)=xtyp(j)
      xtyp(j)=ctemp
      do 110 k=1,jkeep
      temp=sta(i,k)
      sta(i,k)=sta(j,k)
      sta(j,k)=temp
  110 continue
  120 continue
  130 continue
      ls=kdim+1
      ic=0
      do 200 k=1,kdim
      i=jv(k)
      j=kv(k)
      if (i.eq.j) go to 200
      ic=ic+1
      ia=iv(k)
      icc=kdim+ic
      do 190 l=ls,ldim
      if (ia.eq.iv(l).and.i.eq.jv(l).and.j.eq.kv(l)) go to 150
      go to 190
  150 continue
      if (icc.eq.l) go to 200
      it=iv(l)
      iv(l)=iv(icc)
      iv(icc)=it
      it=jv(l)
      jv(l)=jv(icc)
      jv(icc)=it
      it=kv(l)
      kv(l)=kv(icc)
      kv(icc)=it
      ctemp=xtyp(l)
      xtyp(l)=xtyp(icc)
      xtyp(icc)=ctemp
      temp=ga(l)
      ga(l)=ga(icc)
      ga(icc)=temp
      do 160 m=1,jkeep
      temp=sta(l,m)
      sta(l,m)=sta(icc,m)
      sta(icc,m)=temp
  160 continue
      go to 200
  190 continue
  200 continue
c
c     second reordering to separate those configurations which have
c     virtual orbital indeces larger than the number which is going
c     to be included in the nondiagonal b-block from the other ones.
c     the blocks i.le.j.and.i.lt.j.shoulf be kept separately.
c
      if (ntotx.eq.ntot) go to 400
      icc=0
      iend=kdim-1
      do 210 i=1,kdim
      if (iv(i).le.ntotx) go to 210
      icc=icc+1
  210 continue
      ic=kdim-icc
      jstart=kdim+1
      jend=ldim-1
      jcc=0
      do 220 j=jstart,ldim
      if (iv(j).le.ntotx) go to 220
      jcc=jcc+1
  220 continue
      jc=ldim-kdim-jcc
      do 250 i=1,iend
      if (iv(i).le.ntotx) go to 250
      if (i.gt.ic) go to 260
      js=i+1
      do 240 j=js,kdim
      if (iv(j).gt.ntotx) go to 240
      it=iv(i)
      iv(i)=iv(j)
      iv(j)=it
      it=jv(i)
      jv(i)=jv(j)
      jv(j)=it
      it=kv(i)
      kv(i)=kv(j)
      kv(j)=it
      ctemp=xtyp(i)
      xtyp(i)=xtyp(j)
      xtyp(j)=ctemp
      temp=ga(i)
      ga(i)=ga(j)
      ga(j)=temp
      do 230 l=1,jkeep
      temp=sta(i,l)
      sta(i,l)=sta(j,l)
      sta(j,l)=temp
  230 continue
      go to 250
  240 continue
  250 continue
  260 continue
      do 300 i=jstart,jend
      if (iv(i).le.ntotx) go to 300
      if (i.gt.(jc+kdim)) go to 310
      js=i+1
      do 290 j=js,ldim
      if (iv(j).gt.ntotx) go to 290
      it=iv(i)
      iv(i)=iv(j)
      iv(j)=it
      it=jv(i)
      jv(i)=jv(j)
      jv(j)=it
      it=kv(i)
      kv(i)=kv(j)
      kv(j)=it
      ctemp=xtyp(i)
      xtyp(i)=xtyp(j)
      xtyp(j)=ctemp
      temp=ga(i)
      ga(i)=ga(j)
      ga(j)=temp
      do 280 l=1,jkeep
      temp=sta(i,l)
      sta(i,l)=sta(j,l)
      sta(j,l)=temp
  280 continue
      go to 300
  290 continue
  300 continue
  310 continue
c
c     reorder the part with indeces larger than the number of
c     orbitals which is going to be included in the nondiagonal b-block
c     according to the magnitude of the terms in the subblock i.le.j.
c     the other blocks are kept in order.
c
      if (ic.eq.(kdim-1)) go to 395
      iend=kdim-1
      istart=ic+1
      do 340 i=istart,iend
      js=i+1
      do 330 j=js,kdim
      if (ga(i).le.ga(j)) go to 330
      temp=ga(i)
      ga(i)=ga(j)
      ga(j)=temp
      do 320 l=1,jkeep
      temp=sta(i,l)
      sta(i,l)=sta(j,l)
      sta(j,l)=temp
  320 continue
      it=iv(i)
      iv(i)=iv(j)
      iv(j)=it
      it=jv(i)
      jv(i)=jv(j)
      jv(j)=it
      it=kv(i)
      kv(i)=kv(j)
      kv(j)=it
      ctemp=xtyp(i)
      xtyp(i)=xtyp(j)
      xtyp(j)=ctemp
  330 continue
  340 continue
      ks=kdim+jc+1
      kc=0
      do 390 k=istart,kdim
      i=jv(k)
      j=kv(k)
      if (i.eq.j) go to 390
      kc=kc+1
      ia=iv(k)
      kcc=kdim+kc+jc
      do 380 l=ks,ldim
      if (ia.eq.iv(l).and.i.eq.jv(l).and.j.eq.kv(l)) go to 350
      go to 380
  350 continue
      if (kcc.eq.l) go to 390
      it=iv(l)
      iv(l)=iv(kcc)
      iv(kcc)=it
      it=jv(l)
      jv(l)=jv(kcc)
      jv(kcc)=it
      it=kv(l)
      kv(l)=kv(icc)
      kv(kcc)=it
      ctemp=xtyp(l)
      xtyp(l)=xtyp(icc)
      xtyp(icc)=ctemp
      temp=ga(l)
      ga(l)=ga(kcc)
      ga(kcc)=temp
      do 360 m=1,jkeep
      temp=sta(l,m)
      sta(l,m)=sta(kcc,m)
      sta(kcc,m)=temp
  360 continue
      go to 390
  380 continue
  390 continue
  395 continue
  400 continue
      if (tsele.eq.ssele) go to 850
c
c     computation of the maximum possible dimension for the
c     nondiagonal b-block.
c
      nrest=jkeep+ndimct(isym)
      if (ndimct(isym).ne.ndimat(isym)) nrest=nrest+jkeep
      if (ntotx.ne.ntot) go to 420
      if ((nmaxb-nrest).ge.ndimbt(isym)) go to 430
  420 nrest=nrest+jkeep
      number=nmaxb-nrest
      go to 440
  430 continue
      number=nmaxb-nrest
  440 continue
c
c     the nondiagonal b-block should have a dimension smaller or
c     equal to the quantity number
c
      if (ntotx.eq.ntot.and.ndimbt(isym).le.number) go to 820
      ic=0
      jc=0
      do 510 i=1,kdim
      if (iv(i).gt.ntotx) go to 520
      ic=ic+1
      jc=jc+1
      if (jv(i).ne.kv(i)) ic=ic+1
      if (ic.gt.(number-2)) go to 520
  510 continue
  520 continue
      kc=ic-jc
      do 530 k=1,kc
      iv(jc+k)=iv(kdim+k)
      jv(jc+k)=jv(kdim+k)
      kv(jc+k)=kv(kdim+k)
  530 continue
      istart=jc+1
      in=kdim-jc
      jstart=kdim+kc+1
      jn=ldim-jstart+1
      kn=in+jn
      do 550 i=istart,kdim
      ga(i-jc)=ga(i)
      do 540 j=1,jkeep
      sta(i-jc,j)=sta(i,j)
  540 continue
  550 continue
      do 570 j=jstart,ldim
      jj=j-kdim-kc+in
      ga(jj)=ga(j)
      do 560 l=1,jkeep
      sta(jj,l)=sta(j,l)
  560 continue
  570 continue
c
c     calculation of centroids
c
      do 590 i=1,jkeep
      p=0.0d0
      e=0.0d0
      do 580 k=1,kn
      x=sta(k,i)**2
      p=p+x
      e=e+ga(k)*x
  580 continue
      e=e/p
      p= dsqrt(p)
      sigdia(i)=-e
      sta(1,i)=p
  590 continue
      do 800 i=1,jkeep
      ga(i)=sigdia(i)
  800 continue
      do 810 i=1,ic
      iconve(i)=i
  810 continue
      write(nf16)jc,kc,ic,(iv(i),i=1,ic),(jv(i),i=1,ic),(kv(i),i=1,ic),
     *           (xtyp(i),i=1,ic),(iconve(i),i=1,ic)
c
c     structure of data set on unit 16 for block b:
c
c     1st record:
c                jc=number of configurations kept in the block i.le.j
c                kc=number of configurations kept in the block i.lt.j
c                ic=total number of configurations ic=jc+kc
c                iv(i),jv(i),kv(i) are the ia,i,j indeces of the 2h1p
c                                  configurations which are kept.
c                xtyp(i): type of spin coupling: s or t
c                iconve(i): numbering of configurations
c     2nd record:
c                jkeep=number of orbitals of symmetry isym
c                ga(i),i=1,jkeep are the centroids one for each orbital
c                                of symmetry isym
c                sta(1,j),j=1,jkeep are the strengths for each orbital
c
      write (nf16) jkeep,isym,(ga(i),i=1,jkeep),(sta(1,j),j=1,jkeep)
      go to 840
  820 continue
c
c     the complete vectors iv,jv,kv are stored since all elements will
c     be included in the nondiagonal b-block.
c     the diagonal elements will not be stored as they are not needed.
c     a dummy record is written instead.
c
      kc=ldim-kdim
      do 830 i=1,ldim
      iconve(i)=i
  830 continue
      write (nf16) kdim,kc,ldim,(iv(i),i=1,ldim),(jv(i),i=1,ldim),
     *      (kv(i),i=1,ldim),(xtyp(i),i=1,ldim),(iconve(i),i=1,ldim)
      write (nf16) (dummy,i=1,10)
  840 continue
      go to 1000
c
c
c     data handling for the case of pole search solution to the
c     dyson equation
c
c
  850 continue
c
      call diagbx (kdim,ldim,jkeep,
     * nmaxa,nmaxb,nmaxc,nmaxg,nmaxgj,nmaxj,iv,jv,kv,iconve,
     * ivecc,ivecc(nmaxj+1),ivecc(nmaxj+nmaxj+1),ga,sta,aa,iwr)
c
 1000 continue
c
      time4(isym)=cpulft(1)-tim4
      return
      end
      subroutine diagbx (kdim,ldim,jkeep,
     *             nmaxa,nmaxb,nmaxc,nmaxg,nmaxgj,nmaxj,
     *             iv,jv,kv,iconve,itv,jtv,ktv,
     *             ga,sta,aa,iwr)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iv,jv,kv,iconve,itv,jtv,ktv
      character*1 xtyp,txtyp,ctemp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
c
      dimension txtyp(mxtda4)
      dimension iv(nmaxgj),jv(nmaxgj),kv(nmaxgj),iconve(nmaxgj)
      dimension itv(nmaxj),jtv(nmaxj),ktv(nmaxj)
      dimension ga(nmaxg),sta(nmaxg,nmaxc),aa(nmaxa,nmaxa)
c
c
c     npola gives the number of poles taken into account for the
c           a block
c     npolb gives the number of poles taken into account for the
c           b block
c     npolt gives the number of poles taken into account for both
c           together
c
c
    2 format ( //1x,
     + 'the number of diagonally approximated poles for the '/
     + ' b-block for symmetry ',i5,'   is',i8)
    3 format (//1x,
     + 'index identification of these poles ordered',
     +1x,'according to energy'/)
    4 format (/6(2x,i4,' = (',i2,',',i2,',',i2,')',a1))
c
c
c
c     data handling for the case of pole search solution to the
c     dyson equation
c
c
      dummy=0.0d0
      if (ntotx.eq.ntot.and.ldim.le.nmaxb) go to 900
c
c     store indices of configurations of poles taken into account in a
c     nondiagonal form in one set of vectors and the indices of the
c     configurations taken into account diagonally only in another
c     one. write both on disk.
c
c     change sign of poles taken into account diagonally to agree with
c     the sign convention in the routine gammab
c
      ic=0
      jc=0
      do 660 i=1,kdim
      if (iv(i).gt.ntotx) go to 670
      ic=ic+1
      jc=jc+1
      if (jv(i).ne.kv(i)) ic=ic+1
      if (ic.gt.(nmaxb-2)) go to 670
  660 continue
  670 continue
      kc=ic-jc
      do 680 k=1,jc
      itv(k)=iv(k)
      jtv(k)=jv(k)
      ktv(k)=kv(k)
      txtyp(k)=xtyp(k)
  680 continue
      do 690 k=1,kc
      ix=jc+k
      jx=kdim+k
      itv(ix)=iv(jx)
      jtv(ix)=jv(jx)
      ktv(ix)=kv(jx)
      txtyp(ix)=xtyp(jx)
  690 continue
      istart=jc+1
      in=kdim-jc
      jstart=kdim+kc+1
      jn=ldim-jstart+1
      kn=in+jn
      do 710 i=istart,kdim
      ix=i-jc
      iv(ix)=iv(i)
      jv(ix)=jv(i)
      kv(ix)=kv(i)
      xtyp(ix)=xtyp(i)
      ga(ix)=-ga(i)
      do 700 j=1,jkeep
      sta(ix,j)=sta(i,j)
  700 continue
  710 continue
      do 730 j=jstart,ldim
      jj=j-kdim-kc+in
      iv(jj)=iv(j)
      jv(jj)=jv(j)
      kv(jj)=kv(j)
      xtyp(jj)=xtyp(j)
      ga(jj)=-ga(j)
      do 720 l=1,jkeep
      sta(jj,l)=sta(j,l)
  720 continue
  730 continue
c
c     order diagonal poles, strengths and indices
c
      iend=kn-1
      do 760 i=1,iend
      js=i+1
      do 750 j=js,kn
      if (ga(i).le.ga(j)) go to 750
      it=iv(i)
      iv(i)=iv(j)
      iv(j)=it
      it=jv(i)
      jv(i)=jv(j)
      jv(j)=it
      it=kv(i)
      kv(i)=kv(j)
      kv(j)=it
      ctemp=xtyp(i)
      xtyp(i)=xtyp(j)
      xtyp(j)=ctemp
      temp=ga(i)
      ga(i)=ga(j)
      ga(j)=temp
      do 740 k=1,jkeep
      temp=sta(i,k)
      sta(i,k)=sta(j,k)
      sta(j,k)=temp
  740 continue
  750 continue
  760 continue
c
c     calculation of centroids if there are any in the preset limits
c
      if (ga(1).ge.-blimit) go to 850
      do 780 k=1,kn
      kkp=k
      if (ga(k).le.-blimit) go to 780
      go to 790
  780 continue
  790 continue
      kkp=kkp-1
      if (kkp.le.jkeep) go to 850
      do 810 i=1,jkeep
      p=0.0d0
      e=0.0d0
      do 800 k=1,kkp
      x=sta(k,i)**2
      p=p+x
      e=e+ga(k)*x
  800 continue
      e=e/p
      p= dsqrt(p)
      sigdia(i)=e
      aa(1,i)=p
  810 continue
      do 820 i=1,kkp
      ga(i)=0.0d0
      do 820 j=1,jkeep
      sta(i,j)=0.0d0
  820 continue
      do 830 i=1,jkeep
      ga(i)=sigdia(i)
      sta(i,i)=aa(1,i)
  830 continue
c
c     compress vectors after centroid formation
c
      istart=kkp+1
      do 840 i=istart,kn
      ix=i-kkp+jkeep
      ga(ix)=ga(i)
      do 840 j=1,jkeep
      sta(ix,j)=sta(i,j)
  840 continue
      kn=kn-kkp+jkeep
  850 continue
      do 860 i=1,kn
      ix=ic+i
      iconve(ix)=ix
  860 continue
      do 870 i=1,ic
      iconve(i)=i
  870 continue
      npolb(isym,1)=kn
c
c     structure of data set on unit 16 for block b in the case of
c     pole search:
c
c             first record is identical to the case in diagb routine
c             second record:
c                    kn mumber of diagonal poles (including centroids)
c                    jkeep number of orbitals of symmetry isym
c                    ga,sta energies (with sign reversed) and
c                           strengths of diagonal poles
c
      write (nf16) jc,kc,ic,(itv(i),i=1,ic),(jtv(i),i=1,ic),
     *(ktv(i),i=1,ic),(txtyp(i),i=1,ic),(iconve(i),i=1,ic)
c
      write(nf16)kn,jkeep,(ga(i),i=1,kn),((sta(i,j),i=1,kn),j=1,jkeep)
c
c     list all configurations of the diagonally approximated poles
c
      write (iwr,2) isym,kn
      write (iwr,3)
      write (iwr,4) (iconve(i+ic),iv(i),jv(i),kv(i),xtyp(i),i=1,kn)
      go to 1000
  900 continue
      npolb(isym,1)=0
c
c     case all poles are taken into account in the diagonalization
c
      do 910 i=1,ldim
      iconve(i)=i
  910 continue
      kc=ldim-kdim
      write (nf16) kdim,kc,ldim,(iv(i),i=1,ldim),(jv(i),i=1,ldim),
     *(kv(i),i=1,ldim),(xtyp(i),i=1,ldim),(iconve(i),i=1,ldim)
      write (nf16) (dummy,i=1,10)
c
 1000 continue
c
      return
      end
      subroutine dimen(nmaxb,nmaxc,nmaxd,nmaxg,nmaxh,iwr)
c
      implicit REAL  (a-h,o-z)
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
    1 format (//1x,'for symmetry',i3,
     +    '  the number of diagonal elements of block a',i7,
     +    '  exceeds the maximum allowed one',i7/
     +    '  the number of virtual orbitals has been reduced to',i5,
     +    '  and the total number of orbitals to ',i5/
     +    '  for the entire calculation'//)
    2 format (//1x,'for symmetry',i3,
     +    '  the number of diagonal elements of block b',i7,
     +    '  exceeds the maximum allowed one',i7/
     +    '  the number of virtual orbitals has been reduced to',i5,
     +    '  and the total number of orbitals to ',i5/
     +    '  for the entire calculation'//)
    3 format (//1x,'for symmetry',i3,
     +    '  the number of orbitals   ',i3,
     +    '  exceeds the maximum allowed one   ',i3/
     +    '  the total number of orbitals is reduced to',i3,
     +    '  and the number of virtual orbitals is reduced to   ',i3///)
    4 format (//1x,'the number of orbitals for which calculations',
     +    '  are to be performed exceeds the maximum allowed:',
     +    2(3x,i5))
c
      tim2=cpulft(1)
 1000 continue
      do 300 isym=1,nsym
      ndima(isym)=0
      ndimb(isym)=0
      ndimat(isym)=0
      ndimbt(isym)=0
      do 5 i=1,norb
      iorb=indorb(i)
      if (idorb(iorb).eq.isym) go to 6
    5 continue
      go to 300
    6 continue
      jkeep=0
      do 7 i=1,ntot
      if (idorb(i).ne.isym) go to 7
      jkeep=jkeep+1
    7 continue
      ic=0
      do 30 i=1,nocc
      ii=idorb(i)
      ij=imul(ii,isym)
      do 20 ia=nstart,ntot
      iaa=idorb(ia)
      do 10 ib=ia,ntot
      ibb=idorb(ib)
      if (ij.ne.imul(iaa,ibb)) go to 10
      ndima(isym)=ndima(isym)+1
      if (ia.ne.ib) ic=ic+1
   10 continue
   20 continue
   30 continue
      ndimat(isym)=ndima(isym)+ic
      ic=0
      do 60 ia=nstart,ntot
      iaa=idorb(ia)
      ij=imul(iaa,isym)
      do 50 i=1,nocc
      ii=idorb(i)
      do 40 j=i,nocc
      jj=idorb(j)
      if (ij.ne.imul(ii,jj)) go to 40
      ndimb(isym)=ndimb(isym)+1
      if (i.ne.j) ic=ic+1
   40 continue
   50 continue
   60 continue
      ndimbt(isym)=ndimb(isym)+ic
      if (ndimat(isym).le.nmaxg) go to 90
      ntot=ntot-1
      nvirt=nvirt-1
      do 70 nn=1,6
      if ( dabs(epsil(ntot)-epsil(ntot+1)).gt.0.1d-04) go to 80
      ntot=ntot-1
      nvirt=nvirt-1
   70 continue
   80 continue
      if (ntot.lt.ntotx) ntotx=ntot
      if (nvirt.lt.nvirtx) nvirtx=nvirt
      write (iwr,1) isym,ndimat(isym),nmaxg,nvirt,ntot
      go to 1000
   90 continue
      ndimct(isym)=nmaxh
      if (ndimct(isym).gt.ndimat(isym)) ndimct(isym)=ndimat(isym)
c
c     ndimct is the number of affinity poles to be included completely
c
      if (selcta.eq.sselx) go to 95
      nrest=jkeep+ndimct(isym)
      if (ndimct(isym).ne.ndimat(isym)) nrest=nrest+jkeep
      if ((nmaxb-nrest).ge.ndimbt(isym)) go to 120
      nxx=nmaxb-nrest
      go to 96
   95 if (ndimbt(isym).le.nmaxg) go to 120
      nxx=nmaxg
   96 continue
      ntot=ntot-1
      nvirt=nvirt-1
      do 100 nn=1,6
      if ( dabs(epsil(ntot)-epsil(ntot+1)).gt.0.1d-04) go to 110
      ntot=ntot-1
      nvirt=nvirt-1
  100 continue
  110 continue
      write (iwr,2) isym,ndimbt(isym),nxx,nvirt,ntot
      if (ntot.lt.ntotx) ntotx=ntot
      if (nvirt.lt.nvirtx) nvirtx=nvirt
      go to 1000
  120 continue
  300 continue
      do 320 isym=1,nsym
      ic=0
      do 310 i=1,ntot
      if (idorb(i).ne.isym) go to 310
      ic=ic+1
  310 continue
      if (ic.le.nmaxc) go to 320
      ntot=ntot-1
      nvirt=nvirt-1
      do 314 nn=1,6
      if ( dabs(epsil(ntot)-epsil(ntot+1)).gt.0.1d-04) go to 316
      ntot=ntot-1
      nvirt=nvirt-1
  314 continue
  316 continue
      if (ntot.lt.ntotx) ntotx=ntot
      if (nvirt.lt.nvirtx) nvirtx=nvirt
      write (iwr,3) isym,ic,nmaxc,ntot,nvirt
      go to 1000
  320 continue
      if (norb.le.nmaxd) go to 330
      write (iwr,4) norb,nmaxd
      call caserr('invalid number of orbitals in tda')
  330 continue
      do 350 isym=1,lsym
      jend=nsymg(isym)
      do 350 j=1,jend
      kend=isgnum(isym,j)
      ic=0
      do 340 k=1,kend
      if (isgind(isym,j,k).le.ntot) go to 340
      ic=ic+1
  340 continue
      isgnum(isym,j)=isgnum(isym,j)-ic
  350 continue
c
      time2=cpulft(1)-tim2
      return
      end
      subroutine gammaa ( nmaxa,nmaxb,nmaxbg,nmaxc,
     * nmaxf,nmaxg,nmaxh,nmaxi, iv,jv,kv,iconve,iveca,
     * v,g,gg,gvec,aux,aa,
     * indvec,jndvec,kndvec,lndvec,iwr)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iv,jv,kv,iconve,iveca
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
_IFN(ibm,vax)
      common/linkup/iiii,jjjj,kkkk,llll
_ENDIF
c
      dimension iv(nmaxbg),jv(nmaxbg),kv(nmaxbg),iconve(nmaxbg)
      dimension v(nmaxf),g(nmaxb,nmaxb),gg(*),gvec(nmaxb)
      dimension aux(nmaxc,nmaxb),aa(nmaxa,nmaxa)
      dimension iveca(8,nmaxi,1)
c
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
c
    1 format (/1x,104('-'))
    2 format (1x,'the number of configurations taken into account',
     + ' in the a-block for symmetry',i5,'   is',i8)
    3 format (//1x,'index identification of these configurations'/
     +1x,'ordered according to the selection criterion'/
     +1x,'s denotes singlet spin coupling, t triplet spin coupling',
     +1x,'in the configurations'/)
    4 format (/6(2x,i4,' = (',i2,',',i2,',',i2,')',a1))
    5 format (1x,'eigenvalues and large components of',
     +' eigenvectors of matrix for a'
     +/1x,'only the lowest 20 solutions are listed'/1x,
     +'any negative energy solutions arising from linear dependencies'/
     +1x,'in the basis are deleted later'/
     +1x,'if the energies are beyond -10ev'//)
    6 format ('  nr   ',i3,'  ev  ',f10.4,'  ::',10(1x,i3,': ',f5.2))
    7 format (30x,10(1x,i3,': ',f5.2))
c
      timm=cpulft(1)
      do 20 i=1,nmaxb
      do 20 j=1,nmaxb
      g(j,i)=0.0d0
   20 continue
      read (nf15) jn,kn,in,(iv(i),i=1,in),(jv(i),i=1,in),(kv(i),i=1,in),
     *          (xtyp(i),i=1,in),(iconve(i),i=1,in)
      do 100 n=1,jn
      i=iv(n)
      ia=jv(n)
      ib=kv(n)
      do 90 m=n,jn
      ip=iv(m)
      iap=jv(m)
      ibp=kv(m)
      x=0.0d0
      if (i.ne.ip) go to 30
_IF(ibm,vax)
      ind = indvor (ib,ibp,ia,iap)
      jnd = indvor (ib,iap,ibp,ia)
_ELSE
      call indvor (ib,ibp,ia,iap)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,iap,ibp,ia)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)+v(jnd)
   30 if (ia.ne.iap) go to 40
_IF(ibm,vax)
      ind = indvor (ib,ibp,i,ip)
      jnd = indvor (ib,i,ibp,ip)
_ELSE
      call  indvor (ib,ibp,i,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,i,ibp,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)+0.5d0*v(jnd)
   40 if (ib.ne.ibp) go to 50
_IF(ibm,vax)
      ind = indvor (ia,iap,i,ip)
      jnd = indvor (ia,i,iap,ip)
_ELSE
      call indvor (ia,iap,i,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ia,i,iap,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)+0.5d0*v(jnd)
   50 if (ia.ne.ibp) go to 60
_IF(ibm,vax)
      ind = indvor (ib,iap,i,ip)
      jnd = indvor (ib,i,iap,ip)
_ELSE
      call  indvor (ib,iap,i,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,i,iap,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)+0.5d0*v(jnd)
   60 if (ib.ne.iap) go to 70
_IF(ibm,vax)
      ind = indvor (ibp,ia,i,ip)
      jnd = indvor (ibp,ip,ia,i)
_ELSE
      call indvor (ibp,ia,i,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ibp,ip,ia,i)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)+0.5d0*v(jnd)
   70 continue
      if (ia.eq.ib.and.iap.eq.ibp) go to 80
      go to 85
   80 continue
      x=x*0.5d0
      g(n,m)=x
      go to 90
   85 if (ia.eq.ib.or.iap.eq.ibp) x=x*fac
      g(n,m)=x
   90 continue
  100 continue
      ms=jn+1
      do 200 n=1,jn
      i=iv(n)
      ia=jv(n)
      ib=kv(n)
      do 190 m=ms,in
      ip=iv(m)
      iap=jv(m)
      ibp=kv(m)
      x=0.0d0
      if (ia.ne.iap) go to 110
_IF(ibm,vax)
      ind = indvor (ib,i,ibp,ip)
_ELSE
      call indvor (ib,i,ibp,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)
  110 if (ib.ne.ibp) go to 120
_IF(ibm,vax)
      ind = indvor (ia,i,iap,ip)
_ELSE
      call indvor (ia,i,iap,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)
  120 if (ia.ne.ibp) go to 130
_IF(ibm,vax)
      ind = indvor (ib,i,iap,ip)
_ELSE
      call indvor (ib,i,iap,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)
  130 if (ib.ne.iap) go to 140
_IF(ibm,vax)
      ind = indvor (ibp,ip,ia,i)
_ELSE
      call indvor (ibp,ip,ia,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)
  140 continue
      x=x*fac1
      if (ia.eq.ib) x=x*fac
      g(n,m)=x
  190 continue
  200 continue
      do 300 n=ms,in
      i=iv(n)
      ia=jv(n)
      ib=kv(n)
      do 290 m=n,in
      ip=iv(m)
      iap=jv(m)
      ibp=kv(m)
      x=0.0d0
      if (i.ne.ip) go to 210
_IF(ibm,vax)
      ind = indvor (ib,ibp,ia,iap)
      jnd = indvor (ib,iap,ibp,ia)
_ELSE
      call indvor (ib,ibp,ia,iap)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,iap,ibp,ia)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)-v(jnd)
  210 if (ia.ne.iap) go to 220
_IF(ibm,vax)
      ind = indvor (ib,ibp,i,ip)
      jnd = indvor (ib,i,ibp,ip)
_ELSE
      call indvor (ib,ibp,i,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,i,ibp,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)+1.5d0*v(jnd)
  220 if (ib.ne.ibp) go to 230
_IF(ibm,vax)
      ind = indvor (ia,iap,i,ip)
      jnd = indvor (ia,i,iap,ip)
_ELSE
      call indvor (ia,iap,i,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ia,i,iap,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)+1.5d0*v(jnd)
  230 if (ia.ne.ibp) go to 240
_IF(ibm,vax)
      ind = indvor (ib,iap,i,ip)
      jnd = indvor (ib,i,iap,ip)
_ELSE
      call indvor (ib,iap,i,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,i,iap,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)-1.5d0*v(jnd)
  240 if (ib.ne.iap) go to 250
_IF(ibm,vax)
      ind = indvor (ibp,ia,i,ip)
      jnd = indvor (ibp,ip,ia,i)
_ELSE
      call indvor (ibp,ia,i,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ibp,ip,ia,i)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)-1.5d0*v(jnd)
  250 continue
      g(n,m)=x
  290 continue
  300 continue
      do 310 j=1,in
      i=iv(j)
      ia=jv(j)
      ib=kv(j)
      g(j,j)=g(j,j)-epsil(i)+epsil(ia)+epsil(ib)
  310 continue
      ie=in-1
      do 330 i=1,ie
      js=i+1
      do 320 j=js,in
      g(j,i)=g(i,j)
  320 continue
  330 continue
      write (iwr,2) isym,in
      write (iwr,3)
      write (iwr,4) ( i,iv(i),jv(i),kv(i) ,xtyp(i),i=1,in)
      write (iwr,1)
      write (iwr,5)
      time5(isym)=cpulft(1)-timm
      ndim=in
      call sdiag2 (nmaxb,ndim,g,gvec,g)
      timm=cpulft(1)
      iend=20
      if (ndim.lt.20) iend=ndim
      do 340 i=1,iend
      jc=0
      ix=i
      x=gvec(i)*conver
      ic=0
      do 335 j=1,ndim
      if ( dabs(g(j,i)).lt.0.2d0) go to 335
      ic=ic+1
      mvec(ic)=j
      amvec(ic)=g(j,i)
      if (ic.lt.8) go to 335
      write (iwr,6) ix,x,(mvec(k),amvec(k),k=1,ic)
      jc=1
      ic=0
  335 continue
      if (ic.eq.0) go to 340
      if (jc.eq.0) write (iwr,6) ix,x,(mvec(k),amvec(k),k=1,ic)
      if (jc.eq.1) write (iwr,7) (mvec(k),amvec(k),k=1,ic)
  340 continue
      jkeepx=0
      do 345 i=1,ntotx
      if (idorb(i).ne.isym) go to 345
      jkeepx=jkeepx+1
  345 continue
      jkeep=0
      do 350 i=1,ntot
      if (idorb(i).ne.isym) go to 350
      jkeep=jkeep+1
  350 continue
c
c     calculation of the amplitudes of the residues of the self-energy
c
      do 355 i=1,nmaxb
      do 355 j=1,nmaxc
      aux(j,i)=0.0d0
  355 continue
      ic=0
      do 490 ip=1,ntotx
      if (idorb(ip).ne.isym) go to 490
      ic=ic+1
      jc=0
      do 450 m=1,jn
      i=iv(m)
      ia=jv(m)
      ib=kv(m)
      ii=idorb(i)
      iaa=idorb(ia)
      ibb=idorb(ib)
      iim=imul(ii,isym)
      if (ia.eq.ib) go to 400
      jc=jc+1
      if (ia.lt.ip) go to 360
      ind = indvx (ib,i,ia,ip)
      jnd = indvx (ib,ip,ia,i)
      go to 390
  360 if (ip.ge.ib) go to 370
      ind = indvx (ib,i,ip,ia)
      jnd = indvx (ib,ip,ia,i)
      go to 390
  370 ind = indvx (ip,ia,ib,i)
      jnd = indvx (ip,ib,ia,i)
  390 continue
      xa=v(ind)
      xb=v(jnd)
      xc=0.0d0
      del=-epsil(ia)-epsil(ib)
      if (iaa.ne.ibb) go to 810
      do 800 j=1,nocc
      ind = indvx (ib,j,ia,j)
_IF(ibm,vax)
      jnd = indvor (ip,j,i,j)
_ELSE
      call  indvor (ip,j,i,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+(v(ind)*v(jnd))/(2.0d0*epsil(j)+del)
  800 continue
      xc=xc+xc
  810 continue
      xd=0.0d0
      xe=0.0d0
      if (nocc.eq.1) go to 850
      jend=nocc-1
      do 840 j=1,jend
      ks=j+1
      jj=idorb(j)
      eps=epsil(j)+del
      do 830 k=ks,nocc
      kk=idorb(k)
      if(iim.ne.imul(jj,kk)) go to 830
      ind = indvx (ib,k,ia,j)
      jnd = indvx (ib,j,ia,k)
_IF(ibm,vax)
      knd = indvor (ip,j,i,k)
      lnd = indvor (ip,k,i,j)
_ELSE
      call indvor (ip,j,i,k)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,k,i,j)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xx=1.0d0/(eps+epsil(k))
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xd=xd+((x1+x2)*(x3+x4))*xx
      xe=xe+((x1-x2)*(x3-x4))*xx
  830 continue
  840 continue
  850 continue
      del=epsil(i)-epsil(ib)
      xf=0.0d0
      xg=0.0d0
      ibm=imul(ibb,ii)
      do 870 j=1,nocc
      jj=idorb(j)
      eps=del+epsil(j)
      do 860 kc=nstart,ntotx
      kcc=idorb(kc)
      if (imul(jj,kcc).ne.ibm) go to 860
_IF(ibm,vax)
      ind = indvor (ib,i,kc,j)
      jnd = indvor (ia,ip,kc,j)
      knd = indvor (kc,ia,ip,j)
      lnd = indvor (kc,i,ib,j)
_ELSE
      call indvor (ib,i,kc,j)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ia,ip,kc,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (kc,ia,ip,j)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (kc,i,ib,j)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xx=1.0d0/(eps-epsil(kc))
      xf=xf-xx*(x2*(-2.0d0*x1+x4)+x3*(x1+x4))
      xg=xg-xx*(x2*(-2.0d0*x1+x4)+x3*(x1-x4))
  860 continue
  870 continue
      del=epsil(i)-epsil(ia)
      xh=0.0d0
      xi=0.0d0
      ibm=imul(iaa,ii)
      do 890 j=1,nocc
      jj=idorb(j)
      eps=del+epsil(j)
      do 880 kc=nstart,ntotx
      kcc=idorb(kc)
      if (imul(jj,kcc).ne.ibm) go to 880
_IF(ibm,vax)
      ind = indvor (ia,i,kc,j)
      jnd = indvor (ib,ip,kc,j)
      knd = indvor (ib,kc,ip,j)
      lnd = indvor (ia,j,kc,i)
_ELSE
      call indvor (ia,i,kc,j)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,ip,kc,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,kc,ip,j)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,j,kc,i)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xx=1.0d0/(eps-epsil(kc))
      xh=xh-xx*(x2*(-2.0d0*x1+x4)+x3*(x1+x4))
      xi=xi-xx*(x2*(2.0d0*x1-x4)+x3*(-x1+x4))
  880 continue
  890 continue
      x=(xa+xb+xc+xd+xf+xh)*fac
      xx=(xa-xb+xe+xg+xi)*fad
      jx=jn+jc
      do 895 icount=1,ndim
      aux(ic,icount)=aux(ic,icount)+x*g(m,icount)+xx*g(jx,icount)
  895 continue
      go to 450
  400 continue
      if (ia.lt.ip) go to 420
      if (ip.lt.i) go to 410
      ind = indvx (ia,ip,ib,i)
      go to 430
  410 ind = indvx (ib,i,ia,ip)
      go to 430
  420 ind = indvx (ip,ia,ib,i)
  430 continue
      xa=v(ind)
      xc=0.0d0
      del=-2.0d0*epsil(ia)
      do 900 j=1,nocc
      ind = indvx (ia,j,ia,j)
_IF(ibm,vax)
      jnd = indvor (ip,j,i,j)
_ELSE
      call  indvor (ip,j,i,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+(v(ind)*v(jnd))/(2.0d0*epsil(j)+del)
  900 continue
      xd=0.0d0
      if (nocc.eq.1) go to 930
      jend=nocc-1
      do 920 j=1,jend
      ks=j+1
      jj=idorb(j)
      eps=epsil(j)+del
      do 910 k=ks,nocc
      kk=idorb(k)
      if (jj.ne.kk) go to 910
      ind = indvx (ia,k,ia,j)
_IF(ibm,vax)
      jnd = indvor (ip,j,i,k)
      knd = indvor (ip,k,i,j)
_ELSE
      call  indvor (ip,j,i,k)
      jnd = indvx(iiii,jjjj,kkkk,llll)
       call indvor (ip,k,i,j)
       knd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xd=xd+(v(ind)*(v(jnd)+v(knd)))/(epsil(k)+eps)
  910 continue
  920 continue
  930 continue
      del=epsil(i)-epsil(ia)
      xe=0.0d0
      ibm=imul(ii,iaa)
      do 950 j=1,nocc
      jj=idorb(j)
      eps=del+epsil(j)
      do 940 kc=nstart,ntotx
      kcc=idorb(kc)
      if (ibm.ne.imul(jj,kcc)) go to 940
_IF(ibm,vax)
      ind = indvor (ia,i,kc,j)
      jnd = indvor (ia,ip,kc,j)
      knd = indvor (ia,kc,ip,j)
      lnd = indvor (ia,j,kc,i)
_ELSE
      call indvor (ia,i,kc,j)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ia,ip,kc,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,kc,ip,j)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,j,kc,i)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xe=xe-((x2*(-2.0d0*x1+x4)+x3*(x1+x4)))/(eps-epsil(kc))
  940 continue
  950 continue
      x=xa+xc+xd+xe
      do 960 icount=1,ndim
      aux(ic,icount)=aux(ic,icount)+x*g(m,icount)
  960 continue
  450 continue
  490 continue
      time6(isym)=cpulft(1)-timm
c
      write(iwr,1)
      if (tsele.ne.ssele) call datgaa (jkeep,jkeepx,
     * nmaxb,nmaxc,nmaxg,nmaxh,nmaxi,
     * iveca,g,gvec,aux,gg,gg(nmaxg+1),aa,iwr)
      if (tsele.eq.ssele) call patgaa (jkeep,jkeepx,
     * nmaxb,nmaxbg,nmaxc,nmaxg,nmaxi,iconve,iveca
     *,g,gvec,aux,gg(1),gg(nmaxbg+1),gg(nmaxbg*(nmaxc+1)+1),aa,iwr)
c
      return
      end
      subroutine gammab(nmaxa,nmaxb,nmaxbg,nmaxc,nmaxf,nmaxg,nmaxi,
     *                  iv,jv,kv,iconve,ivecb,
     *                  v,g,gg,gvec,aux,aa,indvec,jndvec,kndvec,lndvec
     +                  ,iwr)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iv,jv,kv,iconve,ivecb
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
_IFN(ibm,vax)
      common/linkup/iiii,jjjj,kkkk,llll
_ENDIF
c
      dimension iv(nmaxbg),jv(nmaxbg),kv(nmaxbg),iconve(nmaxbg)
      dimension ivecb(8,nmaxi,*)
      dimension v(nmaxf),g(nmaxb,nmaxb),gg(1),gvec(nmaxb)
      dimension aux(nmaxc,nmaxb),aa(nmaxa,nmaxa)
c
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
    1 format (/1x,104('-'))
    2 format (1x,'the number of configurations taken into account',
     +           ' in the b-block for symmetry',i5,'   is',i8)
    3 format (//1x,'index identification of the 2h1p configurations',
     +             ' ordered according to the selection criterion'/
     +             ' s denotes singlet coupling, t triplet coupling',
     +             ' in the configurations'/)
    4 format (/6(2x,i4,' = (',i2,',',i2,',',i2,')',a1))
    5 format(/1x,'eigenvalues and large components of the ',
     +           'eigenvectors of the matrix for b'/
     +        ' only the eigenvalues between 0 and min(e(i))-30 ev',
     +        ' for valence orbital i are listed'//)
    6 format ('  nr   ',i3,'  ev  ',f10.4,'  ::',10(1x,i3,': ',f5.2))
    7 format (30x,10(1x,i3,': ',f5.2))
c
      timm=cpulft(1)
      do 20 i=1,nmaxb
      do 20 j=1,nmaxb
      g(j,i)=0.0d0
   20 continue
      read (nf16) jn,kn,in,(iv(i),i=1,in),(jv(i),i=1,in),(kv(i),i=1,in),
     *          (xtyp(i),i=1,in),(iconve(i),i=1,in)
      do 100 n=1,jn
      ia=iv(n)
      i=jv(n)
      j=kv(n)
      do 90 m=n,jn
      iap=iv(m)
      ip=jv(m)
      jp=kv(m)
      x=0.0d0
      if (ia.ne.iap) go to 30
_IF(ibm,vax)
      ind = indvor (jp,j,ip,i)
      jnd = indvor (jp,i,j,ip)
_ELSE
      call indvor (jp,j,ip,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (jp,i,j,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)-v(jnd)
   30 if (i.ne.ip) go to 40
_IF(ibm,vax)
      ind = indvor (iap,ia,jp,j)
      jnd = indvor (iap,jp,ia,j)
_ELSE
      call indvor (iap,ia,jp,j)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (iap,jp,ia,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)-0.5d0*v(jnd)
   40 if (j.ne.jp) go to 50
_IF(ibm,vax)
      ind = indvor (iap,ia,ip,i)
      jnd = indvor (iap,ip,ia,i)
_ELSE
      call indvor (iap,ia,ip,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (iap,ip,ia,i)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)-0.5d0*v(jnd)
   50 if (i.ne.jp) go to 60
_IF(ibm,vax)
      ind = indvor (iap,ia,j,ip)
      jnd = indvor (iap,ip,ia,j)
_ELSE
      call indvor (iap,ia,j,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (iap,ip,ia,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)-0.5d0*v(jnd)
   60 if (j.ne.ip) go to 70
_IF(ibm,vax)
      ind = indvor (iap,ia,jp,i)
      jnd = indvor (iap,jp,ia,i)
_ELSE
      call indvor (iap,ia,jp,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (iap,jp,ia,i)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)-0.5d0*v(jnd)
   70 continue
      if (i.eq.j.and.ip.eq.jp) go to 80
      go to 85
   80 continue
      x=x*0.5d0
      g(n,m)=x
      go to 90
   85 if (i.eq.j.or.ip.eq.jp) x=x*fac
      g(n,m)=x
   90 continue
  100 continue
      ms=jn+1
      do 200 n=1,jn
      ia=iv(n)
      i=jv(n)
      j=kv(n)
      do 190 m=ms,in
      iap=iv(m)
      ip=jv(m)
      jp=kv(m)
      x=0.0d0
      if (i.ne.ip) go to 110
_IF(ibm,vax)
      ind = indvor (iap,jp,ia,j)
_ELSE
      call indvor (iap,jp,ia,j)
      ind = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)
  110 if (j.ne.jp) go to 120
_IF(ibm,vax)
      ind = indvor (iap,ip,ia,i)
_ELSE
      call indvor (iap,ip,ia,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)
  120 if (i.ne.jp) go to 130
_IF(ibm,vax)
      ind = indvor (iap,ip,ia,j)
_ELSE
      call indvor (iap,ip,ia,j)
      ind = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)
  130 if (j.ne.ip) go to 140
_IF(ibm,vax)
      ind = indvor (iap,jp,ia,i)
_ELSE
      call indvor (iap,jp,ia,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)
  140 continue
      x=x*fac1
      if (i.eq.j) x=x*fac
      g(n,m)=x
  190 continue
  200 continue
      do 300 n=ms,in
      ia=iv(n)
      i=jv(n)
      j=kv(n)
      do 290 m=n,in
      iap=iv(m)
      ip=jv(m)
      jp=kv(m)
      x=0.0d0
      if (ia.ne.iap) go to 210
_IF(ibm,vax)
      ind = indvor (jp,j,ip,i)
      jnd = indvor (jp,i,j,ip)
_ELSE
      call indvor (jp,j,ip,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (jp,i,j,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)+v(jnd)
  210 if (i.ne.ip) go to 220
_IF(ibm,vax)
      ind = indvor (iap,ia,jp,j)
      jnd = indvor (iap,jp,ia,j)
_ELSE
      call indvor (iap,ia,jp,j)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (iap,jp,ia,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)-1.5d0*v(jnd)
  220 if (j.ne.jp) go to 230
_IF(ibm,vax)
      ind = indvor (iap,ia,ip,i)
      jnd = indvor (iap,ip,ia,i)
_ELSE
      call indvor (iap,ia,ip,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (iap,ip,ia,i)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x+v(ind)-1.5d0*v(jnd)
  230 if (i.ne.jp) go to 240
_IF(ibm,vax)
      ind = indvor (iap,ia,j,ip)
      jnd = indvor (iap,ip,ia,j)
_ELSE
      call indvor (iap,ia,j,ip)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (iap,ip,ia,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)+1.5d0*v(jnd)
  240 if (j.ne.ip) go to 250
_IF(ibm,vax)
      ind = indvor (iap,ia,jp,i)
      jnd = indvor (iap,jp,ia,i)
_ELSE
      call indvor (iap,ia,jp,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (iap,jp,ia,i)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-v(ind)+1.5d0*v(jnd)
  250 continue
      g(n,m)=x
  290 continue
  300 continue
      do 310 k=1,in
      ia=iv(k)
      i=jv(k)
      j=kv(k)
      g(k,k)=g(k,k)-epsil(ia)+epsil(i)+epsil(j)
  310 continue
      ie=in-1
      do 330 i=1,ie
      js=i+1
      do 330 j=js,in
      g(j,i)=g(i,j)
c 320 continue
  330 continue
      write (iwr,2) isym,in
      write (iwr,3)
      write (iwr,4)  (i,iv(i),jv(i),kv(i) ,xtyp(i),i=1,in)
      write (iwr,1)
      write (iwr,5)
      time8(isym)=cpulft(1)-timm
      ndim=in
      call sdiag2(nmaxb,ndim,g,gvec,g)
      timm=cpulft(1)
      xa=epsil(1)-1.0d0
      if (xa.lt.-3.0d0) xa=-3.0d0
      do 350 i=1,ndim
      if (gvec(i).lt.xa) go to 350
      jc=0
      ix=i
      x=gvec(i)*conver
      ic=0
      do 340 j=1,ndim
      if ( dabs(g(j,i)).lt.0.2d0) go to 340
      ic=ic+1
      mvec(ic)=j
      amvec(ic)=g(j,i)
      if (ic.lt.8) go to 340
      write (iwr,6) ix,x,(mvec(k),amvec(k),k=1,ic)
      jc=1
      ic=0
  340 continue
      if (ic.eq.0) go to 350
      if (jc.eq.0) write (iwr,6) ix,x,(mvec(k),amvec(k),k=1,ic)
      if (jc.eq.1) write (iwr,7) (mvec(k),amvec(k),k=1,ic)
  350 continue
      jkeep=0
      do 360 i=1,ntot
      if (idorb(i).ne.isym) go to 360
      jkeep=jkeep+1
  360 continue
c
c     calculation of amplitudes of the residues of the self-energy
c
      do 370 i=1,nmaxb
      do 370 j=1,nmaxc
      aux(j,i)=0.0d0
  370 continue
      ic=0
      do 490 ip=1,ntotx
      if (idorb(ip).ne.isym) go to 490
      ic=ic+1
      jc=0
      do 450 m=1,jn
      ia=iv(m)
      i=jv(m)
      j=kv(m)
      ii=idorb(i)
      jj=idorb(j)
      iaa=idorb(ia)
      iim=imul(iaa,isym)
      if (i.eq.j) go to 400
      jc=jc+1
_IF(ibm,vax)
      ind = indvor (ia,j,ip,i)
      jnd = indvor (ia,i,ip,j)
_ELSE
      call indvor (ia,j,ip,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ia,i,ip,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xa=v(ind)
      xb=v(jnd)
      xc=0.0d0
      del=-epsil(i)-epsil(j)
      if(ii.ne.jj) go to 810
      do 800 ib=nstart,ntotx
      ind = indvx (ib,j,ib,i)
_IF(ibm,vax)
      jnd = indvor (ib,ia,ib,ip)
_ELSE
      call  indvor (ib,ia,ib,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc-(v(ind)*v(jnd))/(2.0d0*epsil(ib)+del)
  800 continue
      xc=xc+xc
  810 continue
      xd=0.0d0
      xe=0.0d0
      if (nvirtx.eq.1) go to 850
      ibend=ntotx-1
      do 840 ib=nstart,ibend
      kcs=ib+1
      ibb=idorb(ib)
      eps=epsil(ib)+del
      do 830 kc=kcs,ntotx
      kcc=idorb(kc)
      if (iim.ne.imul(ibb,kcc)) go to 830
      ind = indvx (kc,j,ib,i)
      jnd = indvx (kc,i,ib,j)
_IF(ibm,vax)
      knd = indvor (kc,ia,ib,ip)
      lnd = indvor (kc,ip,ia,ib)
_ELSE
      call indvor (kc,ia,ib,ip)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (kc,ip,ia,ib)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xx=1.0d0/(eps+epsil(kc))
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xd=xd-xx*((x1+x2)*(x3+x4))
      xe=xe-xx*((x1-x2)*(x3-x4))
  830 continue
  840 continue
  850 continue
      del=epsil(ia)-epsil(j)
      xf=0.0d0
      xg=0.0d0
      ibm=imul(iaa,jj)
      do 870 k=1,nocc
      kk=idorb(k)
      eps=del-epsil(k)
      do 860 ib=nstart,ntotx
      ibb=idorb(ib)
      if (imul(kk,ibb).ne.ibm) go to 860
_IF(ibm,vax)
      ind = indvor (ia,j,ib,k)
      jnd = indvor (ib,k,ip,i)
      knd = indvor (ib,ip,k,i)
      lnd = indvor (ia,k,ib,j)
_ELSE
      call indvor (ia,j,ib,k)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,k,ip,i)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,ip,k,i)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,k,ib,j)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xx=1.0d0/(eps+epsil(ib))
      xf=xf+xx*(x2*(-2.0d0*x1+x4)+x3*(x1+x4))
      xg=xg+xx*(x2*(-2.0d0*x1+x4)+x3*(x1-x4))
  860 continue
  870 continue
      del=epsil(ia)-epsil(i)
      xh=0.0d0
      xi=0.0d0
      ibm=imul(iaa,ii)
      do 890 k=1,nocc
      kk=idorb(k)
      eps=del-epsil(k)
      do 880 ib=nstart,ntotx
      ibb=idorb(ib)
      if (imul(kk,ibb).ne.ibm) go to 880
_IF(ibm,vax)
      ind = indvor (ia,i,ib,k)
      jnd = indvor (ib,k,ip,j)
      knd = indvor (ib,ip,k,j)
      lnd = indvor (ia,k,ib,i)
_ELSE
      call indvor (ia,i,ib,k)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,k,ip,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,ip,k,j)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,k,ib,i)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xx=1.0d0/(eps+epsil(ib))
      xh=xh+xx*(x2*(-2.0d0*x1+x4)+x3*(x1+x4))
      xi=xi+xx*(x2*(2.0d0*x1-x4)+x3*(-x1+x4))
  880 continue
  890 continue
      x=(xa+xb+xc+xd+xf+xh)*fac
      xx=(xa-xb+xe+xg+xi)*fad
      do 895 icount=1,ndim
      aux(ic,icount)=aux(ic,icount)+x*g(m,icount)+
     *                              xx*g(jn+jc,icount)
  895 continue
      go to 450
  400 continue
      if (ia.lt.ip) go to 420
      if (ip.lt.i) go to 410
      ind = indvx (ia,i,ip,i)
      go to 430
  410 ind = indvx (ia,i,i,ip)
      go to 430
  420 ind = indvx (ip,i,ia,i)
  430 continue
      xa=v(ind)
      xc=0.0d0
      del=-2.0d0*epsil(i)
      do 900 ib=nstart,ntotx
      ind = indvx (ib,i,ib,i)
_IF(ibm,vax)
      jnd = indvor (ia,ib,ib,ip)
_ELSE
      call  indvor (ia,ib,ib,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc-(v(ind)*v(jnd))/(2.0d0*epsil(ib)+del)
  900 continue
      xd=0.0d0
      if (nvirtx.eq.1) go to 930
      ibend=ntotx-1
      do 920 ib=nstart,ibend
      kcs=ib+1
      ibb=idorb(ib)
      eps=del+epsil(ib)
      do 910 kc=kcs,ntotx
      kcc=idorb(kc)
      if(ibb.ne.kcc) go to 910
      ind = indvx (kc,i,ib,i)
_IF(ibm,vax)
      jnd = indvor (kc,ia,ib,ip)
      knd = indvor (kc,ip,ia,ib)
_ELSE
      call  indvor (kc,ia,ib,ip)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (kc,ip,ia,ib)
      knd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xd=xd-(v(ind)*(v(jnd)+v(knd)))/(epsil(kc)+eps)
  910 continue
  920 continue
  930 continue
      del=epsil(ia)-epsil(i)
      xe=0.0d0
      ibm=imul(iaa,ii)
      do 950 k=1,nocc
      kk=idorb(k)
      eps=del-epsil(k)
      do 940 ib=nstart,ntotx
      ibb=idorb(ib)
      if (ibm.ne.imul(kk,ibb)) go to 940
_IF(ibm,vax)
      ind = indvor (ia,i,ib,k)
      jnd = indvor (ib,k,ip,i)
      knd = indvor (ib,ip,k,i)
      lnd = indvor (ia,k,ib,i)
_ELSE
      call indvor (ia,i,ib,k)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,k,ip,i)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,ip,k,i)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,k,ib,i)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xe=xe+((x2*(-2.0d0*x1+x4)+x3*(x1+x4)))/(eps+epsil(ib))
  940 continue
  950 continue
      x=xa+xc+xd+xe
      do 960 icount=1,ndim
      aux(ic,icount)=aux(ic,icount)+x*g(m,icount)
  960 continue
  450 continue
  490 continue
      time9(isym)=cpulft(1)-timm
c
      if (tsele.ne.ssele) call datgab(
     *      nmaxb,nmaxc,nmaxi,ivecb,g,gvec,aux,aa,iwr)
      jkeepx=0
      do 970 i=1,ntotx
      if (idorb(i).ne.isym) go to 970
      jkeepx=jkeepx+1
  970 continue
      if (tsele.eq.ssele) call patgab(jkeep,jkeepx,
     * nmaxb,nmaxbg,nmaxc,nmaxg,nmaxi,iconve,ivecb,g,gvec,aux,
     * gg(1),gg(nmaxbg+1),gg(nmaxbg*(nmaxc+1)+1),aa,iwr)
      return
      end
      subroutine ltimet(iwr)
c
      implicit REAL  (a-h,o-z)
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
      character*8 selct
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      data         selct / 'adiagram' /
c
    1 format (//1x,104('*')/)
    2 format (2x,f8.1,' seconds in input routine')
    3 format (2x,f8.1,' seconds in dimen routine')
    4 format (2x,f8.1,' seconds in diaga routine for symmetry   ',i3)
    5 format (2x,f8.1,' seconds in diagb routine for symmetry   ',i3)
    6 format (2x,f8.1,' seconds in gammaa routine for matrix',
     +                ' construction for symmetry   ',i3)
    7 format (2x,f8.1,' seconds in gammaa routine for calculation',
     +                ' of strengths for symmetry   ',i3)
    8 format (2x,f8.1,' seconds in datgaa routine for symmetry   ',i3)
    9 format (2x,f8.1,' seconds in gammab routine for matrix',
     +                ' construction for symmetry   ',i3)
   10 format (2x,f8.1,' seconds in gammab routine for calculation',
     +                ' of strengths for symmetry   ',i3)
   11 format (2x,f8.1,' seconds in datgab routine for symmetry   ',i3)
   12 format (2x,f8.1,' seconds in comdig routine',
     +                ' without the diagonalizations')
   13 format (2x,f8.1,' seconds in calculation of quartet states')
   14 format (2x,f8.1,' seconds in diagonalization of a matrix',
     +                ' of dimension   ',i6)
   15 format (2x,f8.1,' seconds in calculation of constant diagrams',
     + ' with hartree-fock GF, this is included above ')
   16 format (//' total computation time (seconds)  =  ',f10.1)
   17 format (2x,f8.1,' seconds in obtaining solutions via ',
     +                'the pole search procedure')
c
      write (iwr,1)
      write (iwr,2) time1
      write (iwr,3) time2
      do 60 isym=1,nsym
      do 40 i=1,norb
      iorb=indorb(i)
      if (idorb(iorb).eq.isym) go to 50
   40 continue
      go to 60
   50 continue
      write (iwr,4) time3(isym),isym
      write (iwr,5) time4(isym),isym
      write (iwr,6) time5(isym),isym
      write (iwr,7) time6(isym),isym
      write (iwr,8) time7(isym),isym
      write (iwr,9) time8(isym),isym
      write (iwr,10) time9(isym),isym
      write (iwr,11) time10(isym),isym
   60 continue
      if (tsele.eq.ssele) go to 220
      write (iwr,12) time11
      if (word.ne.selct) go to 70
      write (iwr,15) time14
   70 continue
      if (selctd.ne.sseld) go to 80
      write (iwr,13) time12
   80 continue
      do 100 i=1,icdiag
      write (iwr,14) time13(i),idimdi(i)
  100 continue
      x=time1+time2+time11
      if (selctd.eq.sseld) x=x+time12
      do 140 isym=1,nsym
      do 110 i=1,norb
      iorb=indorb(i)
      if (idorb(iorb).eq.isym) go to 120
  110 continue
      go to 140
  120 continue
      x=x+time3(isym)+time4(isym)+time5(isym)+time6(isym)+time7(isym)
     *   +time8(isym)+time9(isym)+time10(isym)
  140 continue
      do 160 i=1,icdiag
      x=x+time13(i)
  160 continue
      write (iwr,16) x
      go to 300
  220 continue
      write (iwr,17) time11
      if (word.ne.selct) go to 240
      write (iwr,15) time14
  240 continue
      if (selctd.ne.sseld) go to 250
      write (iwr,13) time12
  250 continue
      do 255 i=1,icdiag
      write (iwr,14) time13(i),idimdi(i)
  255 continue
      x=time1+time2+time11
      if (selctd.eq.sseld) x=x+time12
      do 280 isym=1,nsym
      do 260 i=1,norb
      iorb=indorb(i)
      if (idorb(iorb).eq.isym) go to 270
  260 continue
      go to 280
  270 continue
      x=x+time3(isym)+time4(isym)+time5(isym)+time6(isym)+time7(isym)
     *   +time8(isym)+time9(isym)+time10(isym)
  280 continue
      do 290 i=1,icdiag
      x=x+time13(i)
  290 continue
      write (iwr,16) x
  300 continue
      return
      end
      subroutine pcorre (e,factor,isum,nd,jp,
     *        nmaxa,nmaxb,maxbg2,mxbg22,nmaxc,nmaxf,
     *        v,g,bb,veca,bmat
_IF(ibm,vax)
     * ,iwr )
_ELSE
     * ,indvec,jndvec,kndvec,lndvec,iwr )
_ENDIF
c
c
      implicit REAL  (a-h,o-z)
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
_IFN(ibm,vax)
      common/linkup/iiii,jjjj,kkkk,llll
_ENDIF
      dimension xorb(mxtda2)
      dimension v(nmaxf),g(nmaxb,nmaxb),bb(nmaxa,nmaxa)
      dimension veca(mxbg22),bmat(nmaxc,maxbg2)
      equivalence (xorb(1),sigdia(1))
_IFN(ibm,vax)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
_ENDIF
c
    1 format(//1x,'orbital contribution to the correlation energy',
     + ' from ip and residues of the green function part for symmetry '
     +  ,i3,'   subpart   ',i3/(10(2x,f10.6)))
    2 format(1x,'orbital contribution to the correlation energy',
     +  ' from green function times orbital energy part for symmetry ',
     +   i3,'   subpart   ',i3/(10(2x,f10.6)))
    3 format(1x,'orbital contribution to the correlation energy',
     +  ' from the diagonal part of the interaction times'/
     +   ' green function part for symmetry   ',i3,'   subpart   ',i3/
     +    (10(2x,f10.6)))
    4 format(1x,'total part of this ( diagonal ) correlation ',
     +    'contribution   ',f10.6)
    5 format(1x,'total part of this correlation including the ',
     +    'nondiagonal part   ',f10.6)
    6 format(1x,'total part of this correlation energy   ',f10.6)
    7 format(1x,'correlation energy contribution of the core ',
     +     'interaction part   ',f10.6)
    8 format(/1x,'the total correlation energy after adding',
     + ' the contribution from symmetry type ',i3,' subpart ',i3,
     + ' is ',f10.6)
c
      ex=0.0d0
      do 20 i=1,nmaxa
      xorb(i)=0.0d0
   20 continue
      nc=ixxx(isym,jp)
      do 40 i=1,isum
      x=0.0d0
      do 30 j=nc,nd
      x=x-veca(j)*(bmat(i,j)**2)
   30 continue
      x=x*factor
      xorb(i)=x
      ex=ex+x
   40 continue
      write (iwr,1) isym,jp,(xorb(i),i=1,isum)
      write (iwr,6) ex
      e=e+ex
      ex=0.0d0
      do 50 i=1,nmaxa
      do 50 j=1,nmaxa
      bb(j,i)=0.0d0
   50 continue
      do 80 k=1,isum
      i=isgind(isym,jp,k)
      do 70 l=k,isum
      j=isgind(isym,jp,l)
      x=0.0d0
      do 60 nn=nc,nd
      x=x+bmat(k,nn)*bmat(l,nn)
   60 continue
      bb(i,j)=x
      bb(j,i)=x
      if (i.eq.j.and.i.le.nocc) bb(i,i)=bb(i,i)-1.0d0
   70 continue
   80 continue
      do 85 i=1,nmaxa
      xorb(i)=0.0d0
   85 continue
      do 90 k=1,isum
      i=isgind(isym,jp,k)
      x=bb(i,i)*epsil(i)
      x=x*factor
      xorb(k)=x
      ex=ex+x
      if (i.gt.nocc) go to 90
      ex=ex-epsil(i)*factor
      xorb(k)=xorb(k)-epsil(i)*factor
   90 continue
      write (iwr,2) isym,jp,(xorb(i),i=1,isum)
      write (iwr,6) ex
      e=e+ex
      ex=0.0d0
      do 95 i=1,nmaxa
      xorb(i)=0.0d0
   95 continue
      do 120 l=1,isum
      i=isgind(isym,jp,l)
      do 110 m=l,isum
      j=isgind(isym,jp,m)
      x=0.0d0
      do 100 k=1,nocc
_IF(ibm,vax)
      ind = indvor (i,j,k,k)
      jnd = indvor (i,k,j,k)
_ELSE
      call indvor (i,j,k,k)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (i,k,j,k)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x=x-2.0d0*v(ind)+v(jnd)
  100 continue
      if (i.ne.j) x=x+x
      x=x*factor
      ex=ex+x*bb(i,j)
      if (i.eq.j) xorb(l)=xorb(l)+x*bb(i,i)
  110 continue
  120 continue
      write (iwr,3) isym,jp,(xorb(i),i=1,isum)
      x=0.0d0
      do 125 i=1,isum
      x=x+xorb(i)
  125 continue
      write (iwr,4) x
      write (iwr,5) ex
      e=e+ex
      ex=0.0d0
      rewind nf14
      read (nf14) nx,((g(j,i),j=1,nx),i=1,nx)
      rewind nf14
      do 130 i=1,ntot
      do 130 j=1,ntot
      ex=ex+g(j,i)*bb(j,i)
  130 continue
      ex=ex*factor
      write (iwr,7) ex
      e=e+ex
      ex=0.0d0
      write (iwr,8) isym,jp,e
      do 300 i=1,nmaxa
      do 300 j=1,nmaxa
      bb(j,i)=0.0d0
  300 continue
      return
      end
      subroutine polser(nmaxa,nmaxb,nmaxbg,maxbg2,mxbg22,
     *                  nmaxc,nmaxc2,l2,nmaxf,nmaxg,nmaxi
     *, v,g,wkarea, iv,jv, aa,bb,amat,deriv,almat,const,
     *  diag,xeigen,poleig,worka,rvec,cmat,bmat,nrpol,
     *  veca,vecb,vecc
_IF(ibm,vax)
     *,iwr )
_ELSE
     * ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iv,jv
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
      character*8 selct
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      common /cgraph/ nsolc,ndummy,en(200),ps(200),norbin(200),iccc(100)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
_IFN(ibm,vax)
      common/linkup/iiii,jjjj,kkkk,llll
_ENDIF
c
      dimension v(nmaxf),g(nmaxb,nmaxb),wkarea(*)
      dimension iv(nmaxbg),jv(nmaxbg)
      dimension aa(nmaxa,nmaxa),bb(nmaxa,nmaxa)
      dimension amat(l2),deriv(l2),almat(l2),const(nmaxc2),
     * diag(nmaxc),xeigen(nmaxc2),poleig(nmaxc2),worka(nmaxc),
     * rvec(nmaxc,nmaxc)
      dimension cmat(maxbg2),bmat(nmaxc,maxbg2),nrpol(maxbg2)
      dimension veca(mxbg22),vecb(mxbg22),vecc(mxbg22)
_IFN(ibm,vax)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
_ENDIF
      data selct  /'adiagram'/
_IFN(ibm,vax)
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
_ENDIF
c
c
c   1 format (//1x,104('*')/)
    2 format (//1x,'convergence obtained at iteration step   ',i3)
    3 format (//1x,'results of first step with constant diagrams',
     +     ' calculated with the hartree-fock greens function')
    8 format (/1x,'diagonal of matrix of constant diagrams'/)
    9 format (10(2x,f10.6))
   10 format (/1x,'convergence is not obtained, the intermediate',
     +     ' results are printed after the iterations'//)
   11 format (//1x,'results of first step without constant diagrams')
   12 format (/1x,'for symmetry  ',i3,'  in the interval from  ',f10.5,
     *' to ',f10.5,' with solution of the dyson equation  ',f10.5/
     *1x,'the sum of the squared 1h and 2h1p components of ',
     +   'this solution is  ',f10.5)
   13 format (//1x,'correlation energy in first step ',
     +   '(without constant diagrams) = ',f10.6,' a.u. = ',f12.6,
     +   ' ev'///)
   14 format (//1x,'correlation energy including constant diagrams = ',
     +     f10.6,' a.u. = ',f12.6,' ev'///)
   17 format (//1x,'results of second step where the constant',
     + ' diagrams are calculated with the hartree-fock greens function')
   25 format(/1x,
     + 'the number of solutions which could not be found in the pole',
     + ' search procedure is'/1x,
     + 'for symmetry',i5,' and subpart',i5,'   equal to',i6/1x,
     + 'this is due to an insufficient number of iterations',
     + ' or to the fact that the solutions coincide with a pole'/1x,
     + 'it is assumed that the residues of the missing solutions ',
     + 'are zero'//)
c
      tim11=cpulft(1)
      iunit=nf18
      junit=nf19
      iccc(1) = 1
      do 80 i = 2,ntot
      iccc(i) = i
      j = i-1
      if (  dabs(epsil(i)-epsil(j)) .lt. 0.1d-04 ) iccc(i) = iccc(j)
   80 continue
      do 100 i=1,nsym
      do 100 j=1,nmaxi
      ixxx(i,j)=0
  100 continue
      numc=1
      do 110 i=1,numend
      fmul(i)=1.0d0
  110 continue
      do 120 i=1,nmaxa
      do 120 j=1,nmaxa
      aa(j,i)=0.0d0
      bb(j,i)=0.0d0
  120 continue
      if ( numend.eq.-2) nxcont=1
      if ( numend.eq.-1) nxcont=2
      if ( numend.eq.2)  nxcont=3
      if ( numend.eq.6)  nxcont=4
      if ( numend.lt.2) numend=1
      do 2000 num=1,numend
      nsolc = 0
      iprint=0
      if ( num.ne.1.or.nxcont.ne.2) go to 4000
 4010 call psadia(numc,nmaxa,nmaxb,mxbg22,
     * nmaxc,nmaxf,v,g,wkarea,aa,bb,bb,veca,bmat
_IF(ibm,vax)
     * ,iwr)
_ELSE
     * ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
      go to 4020
 4000 if ( num.eq.2) go to 4010
 4020 rewind nf12
      rewind nf13
      rewind nf15
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      do 600 isym=1,nsym
      do 200 i=1,norb
      ind=indorb(i)
      if (idorb(ind).eq.isym) go to 210
  200 continue
      go to 600
  210 continue
      jend=nsymg(isym)
      do 580 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 570
      do 220 i=1,nmaxb
      do 220 j=1,nmaxb
      g(j,i)=0.0d0
  220 continue
      kend=isgnum(isym,jp)
      kkend=ksgnum(isym,jp)
      nz=0
      do 230 k=1,kend
      ind=isgind(isym,jp,k)
      do 230 l=1,k
      jnd=isgind(isym,jp,l)
      nz=nz+1
      amat(nz)=-aa(ind,jnd)
      if (k.eq.l) amat(nz)=amat(nz)-epsil(jnd)
  230 continue
      read (nf15) ix,iy,iz,jz,kz,it,(cmat(i),i=1,kz),
     *     ((bmat(k,i),i=1,kz),k=1,it),(iv(i),i=1,iz),(jv(i),i=1,jz)
      xa=cmat(iz)
      xb=cmat(iz+1)
      write (nf16) iz,jz,(iv(i),i=1,iz),(jv(i),i=1,jz)
      if (iprint.ne.0) go to 250
      write (iwr,8)
      do 240 i=1,ntot
      sigdia(i)=aa(i,i)*conver
  240 continue
      write (iwr,9) (sigdia(i),i=1,ntot)
      iprint=1
  250 continue
c
c     set variables for call to pole search program
c
      crita=1.0d-06
      ifile=nf17
cps 
c     jprint=7
      jprint=6
      ivec=1
      maxm=kz
      maxn=kend
      maxpsa=nmaxc
      mdim=maxbg2
      ndimp=nmaxc
      nd=(ndimp*(ndimp+1))/2
      ndkeep=nd
      ndim2=2*ndimp
      call psa (amat,bmat,cmat,almat,rvec,deriv,const,xeigen,poleig,
     *          nrpol,diag,worka,range,crita,ndimp,nd,ndim2,mdim,maxn,
     *          maxm,jprint,ifile,maxpsa,ivec,ier,neigen)
      rewind nf17
      if(ier.eq.0) go to 280
      write (iwr,25) isym,jp,ier
  280 continue
      npolt(isym,jp)=neigen
      npolb(isym,jp)=kz
      je=(nmaxg+nmaxb)*2
      do 290 j=1,je
      cmat(j)=0.0d0
      do 290 i=1,nmaxc
      bmat(i,j)=0.0d0
  290 continue
      isuma=0
      do 300 k=1,kend
      if (isgind(isym,jp,k).gt.nocc) go to 300
      isuma=isuma+1
  300 continue
      idim=maxn+maxm
c
c     read eigenvalues and eigenvectors from file 17, copy onto file 13.
c     file 13 will contain the results for all symmetries from the pole
c     search.
c
      do 350 n=1,neigen
      read (nf17) eigen,x,y
      cmat(n)=eigen
      ic=0
  310 continue
      read (nf17) (vecb(i+ic),i=1,nd)
      ic=ic+nd
      if (idim.le.ic) go to 320
      go to 310
  320 continue
      write (nf13) neigen,n,idim,eigen,(vecb(i),i=1,idim)
      write (nf12) neigen,n,idim,eigen,(vecb(i),i=1,idim)
      do 330 ip=1,kend
      bmat(ip,n)=vecb(ip)
  330 continue
  350 continue
      rewind nf17
      rewind nf12
c
c     the separation of the affinity and ionization poles of the green's
c     function is done in the following way:
c
c     if the energy is less than that of the first affinity pole of the
c     self-energy, the solution is an affinity pole of the green's
c     function. the interval between the first affinity pole and the
c     first ionization pole of the self-energy is investigated. let n
c     be the eigenvector of the green's function in this interval. if
c     the sum of the squared 1-hole and 2-hole-1-particle components of
c     this eigenvector ( c(p,n)**2 ) is less than 0.60 the solution is
c     taken as an affinity pole, otherwise as an ionization pole of
c     the green's function.
c
      do 370 nn=1,neigen
      nc=nn
      if (cmat(nn).le.xa) go to 360
      go to 380
  360 continue
      read (nf12) dummy
  370 continue
  380 continue
      nsa=kend+iz+1
      nsb=kend+iz+jz
      xa=xa*conver
      xb=xb*conver
  390 continue
      xx=0.0d0
      if (isuma.eq.0) go to 410
      do 400 ip=1,isuma
      xx=xx+bmat(ip,nc)**2
  400 continue
  410 continue
      read (nf12) nx,ny,mx,x,(vecb(i),i=1,mx)
      do 420 ip=nsa,nsb
      xx=xx+vecb(ip)**2
  420 continue
      xc=cmat(nc)*conver
      write (iwr,12) isym,xa,xb,xc,xx
      if (xx.gt.0.60d0) go to 430
      nc=nc+1
      if (cmat(nc)*conver.gt.xb) go to 430
      go to 390
  430 continue
      ixxx(isym,jp)=nc
      rewind nf12
      if (word.ne.selct) go to 550
      if (num.lt.2) go to 550
      if (numend.eq.2) go to 550
      do 540 i=1,kkend
      ind=isgind(isym,jp,i)
      ix=ndeg(ind)
      do 530 j=1,kkend
      jnd=isgind(isym,jp,j)
      jx=ndeg(jnd)
      xx=0.0d0
      do 440 nn=nc,neigen
      xx=xx+bmat(i,nn)*bmat(j,nn)
  440 continue
      if (ind.eq.jnd.and.ind.le.nocc) xx=xx-1.0d0
      do 520 isymp=1,nsym
      ipp=nsymg(isymp)
      do 510 ip=1,ipp
      if (iuse(isymp,ip).eq.0) go to 510
      lend=ksgnum(isymp,ip)
      do 500 k=1,lend
      knd=isgind(isymp,ip,k)
      do 490 l=k,lend
      lnd=isgind(isymp,ip,l)
      if (idega.eq.1) go to 470
      xa=0.0d0
      do 460 m=1,ix
      ia=indeg(ind,m)
      iaa=idorb(ia)
      do 450 mm=1,jx
      ib=indeg(jnd,mm)
      if (iaa.ne.idorb(ib)) go to 450
_IF(ibm,vax)
      ind1= indvor (knd,lnd,ia,ib)
      ind2= indvor (knd,ib,lnd,ia)
_ELSE
      call indvor (knd,lnd,ia,ib)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (knd,ib,lnd,ia)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xa=xa+v(ind1)*2.0d0-v(ind2)
  450 continue
  460 continue
      go to 480
  470 continue
_IF(ibm,vax)
      ind1= indvor (knd,lnd,ind,jnd)
      ind2= indvor (knd,jnd,lnd,ind)
_ELSE
      call indvor (knd,lnd,ind,jnd)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (knd,jnd,lnd,ind)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xa=v(ind1)*2.0d0-v(ind2)
  480 continue
      bb(knd,lnd)=bb(knd,lnd)+xa*xx
      bb(lnd,knd)=bb(knd,lnd)
  490 continue
  500 continue
  510 continue
  520 continue
  530 continue
  540 continue
  550 continue
      do 560 i=1,neigen
      cmat(i)=cmat(i)*conver
  560 continue
      write (iunit) neigen,(cmat(j),j=1,neigen)
      go to 580
  570 continue
      read (nf15) dummy
  580 continue
  600 continue
      rewind nf12
      rewind nf13
      rewind nf15
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      if(num.gt.1) go to 1000
      do 620 i=1,ntot
      do 620 j=i,ntot
      aa(i,j)=bb(i,j)*fmul(num)
      aa(j,i)=aa(i,j)
  620 continue
      do 630 i=1,nmaxa
      do 630 j=1,nmaxa
      bb(j,i)=0.0d0
  630 continue
      it=iunit
      iunit=junit
      junit=it
      if ( nxcont.ne.2) write (iwr,11)
      if ( nxcont.eq.2) write (iwr,3)
      ecor=0.0d0
      do 950 isym=1,nsym
      do 650 iorb=1,norb
      ind=indorb(iorb)
      if (idorb(ind).eq.isym) go to 660
  650 continue
      go to 950
  660 continue
      jend=nsymg(isym)
      do 940 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 940
      kend=isgnum(isym,jp)
      mx=isgind(isym,jp,1)
      mxx=ndeg(mx)
      factor=1.0d0
      if (selctc.eq.ssely.and.idega.eq.0) factor= dfloat(mxx)
      isuma=0
      do 680 k=1,kend
      if (isgind(isym,jp,k).gt.nocc) go to 680
      isuma=isuma+1
  680 continue
      isumb=isuma+1
      if (isumb.gt.kend) isumb=kend
      call outtda (jp,isuma,isumb,1,
     * nmaxbg,maxbg2,mxbg22,nmaxc,
     * iv,jv,veca,vecb,vecc,bmat,iwr)
      if (selctb.eq.zzsel) go to 720
      neigen=npolt(isym,jp)
      call pcorre (ecor,factor,kend,neigen,jp,
     * nmaxa,nmaxb,maxbg2,mxbg22,nmaxc,nmaxf,
     * v,g,bb,veca,bmat
_IF(ibm,vax)
     * ,iwr )
_ELSE
     * ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
  720 continue
  940 continue
  950 continue
      if (selctb.eq.zzsel) go to 960
      ecore=ecor*conver
      write (iwr,13) ecor,ecore
  960 continue
c
c     call graph
c
      rewind nf12
      rewind nf13
      rewind nf15
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      go to 2000
 1000 continue
      if (num.ne.2) go to 1300
      write (iwr,17)
      do 1200 isym=1,nsym
      do 1100 iorb=1,norb
      ind=indorb(iorb)
      if (idorb(ind).eq.isym) go to 1120
 1100 continue
      go to 1200
 1120 continue
      jend=nsymg(isym)
      do 1190 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 1190
      kend=isgnum(isym,jp)
      isuma=0
      do 1130 k=1,kend
      if (isgind(isym,jp,k).gt.nocc) go to 1130
      isuma=isuma+1
 1130 continue
      isumb=isuma+1
      if (isumb.gt.kend) isumb=kend
      call outtda (jp,isuma,isumb,2,
     * nmaxbg,maxbg2,mxbg22,nmaxc,iv,jv,
     * veca,vecb,vecc,bmat,iwr)
 1190 continue
 1200 continue
 1300 continue
      do 1600 isym=1,nsym
      do 1400 iorb=1,norb
      ind=indorb(iorb)
      if (idorb(ind).eq.isym) go to 1420
 1400 continue
      go to 1600
 1420 continue
      if (num.le.numc) go to 1460
      jend=nsymg(isym)
      do 1450 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 1450
      read (iunit) nx,(veca(i),i=1,nx)
      read (junit) nx,(vecb(i),i=1,nx)
      do 1440 i=1,nx
      if (veca(i).lt.0.0d0) go to 1440
      if (veca(i).gt.( dabs(epsil(1))*conver+30.0d0)) go to 1440
      if ( dabs(veca(i)-vecb(i)).gt.0.1d0) go to 1460
c
c     0.1 ev is chosen as the criterion for convergence in the energies
c
 1440 continue
 1450 continue
      go to 1600
 1460 continue
      rewind nf18
      rewind nf19
      do 1480 i=1,ntot
      do 1480 j=i,ntot
      aa(i,j)=bb(i,j)*fmul(num)
      aa(j,i)=aa(i,j)
 1480 continue
      do 1500 i=1,nmaxa
      do 1500 j=1,nmaxa
      bb(j,i)=0.0d0
 1500 continue
      it=iunit
      iunit=junit
      junit=it
      go to 2000
 1600 continue
      go to 2100
 2000 continue
      nsolc=0
      if (numend.eq.1) go to 5000
      write (iwr,10)
      go to 2120
 2100 continue
      nsolc=0
      write (iwr,2) num
 2120 continue
      rewind nf12
      rewind nf13
      rewind nf15
      rewind nf16
      rewind nf17
      rewind nf18
      rewind nf19
      ecor=0.0d0
      do 3000 isym=1,nsym
      do 2200 iorb=1,norb
      ind=indorb(iorb)
      if (idorb(ind).eq.isym) go to 2220
 2200 continue
      go to 3000
 2220 continue
      jend=nsymg(isym)
      do 2900 jp=1,jend
      if (iuse(isym,jp).eq.0) go to 2900
      kend=isgnum(isym,jp)
      mx=isgind(isym,jp,1)
      mxx=ndeg(mx)
      factor=1.0d0
      if (selctc.eq.ssely.and.idega.eq.0) factor= dfloat(mxx)
      isuma=0
      do 2250 k=1,kend
      if (isgind(isym,jp,k).gt.nocc) go to 2250
      isuma=isuma+1
 2250 continue
      isumb=isuma+1
      if (isumb.gt.kend) isumb=kend
      call outtda (jp,isuma,isumb,1,
     * nmaxbg,maxbg2,mxbg22,nmaxc,iv,jv,
     * veca,vecb,vecc,bmat,iwr)
      if (selctb.eq.zzsel) go to 2900
      neigen=npolt(isym,jp)
      call pcorre (ecor,factor,kend,neigen,jp,
     * nmaxa,nmaxb,maxbg2,mxbg22,nmaxc,nmaxf,
     * v,g,bb,veca,bmat
_IF(ibm,vax)
     * ,iwr )
_ELSE
     * ,indvec,jndvec,kndvec,lndvec,iwr)
_ENDIF
 2900 continue
 3000 continue
      if (selctb.eq.zzsel) go to 3100
      ecore=ecor*conver
      write (iwr,14) ecor,ecore
 3100 continue
c
c     call graph
c
 5000 continue
c
      time11=cpulft(1)-tim11
      return
      end
      subroutine diaga(
     * nmaxa,nmaxc,nmaxe,nmaxg,nmaxgj,nmaxj,
     * iv,jv,kv,iconve,ivecc,
     * v,ga,sta,aa,ax ,indvec,jndvec,kndvec,lndvec,iwr)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iv,jv,kv,iconve,ivecc
      character*1 xtyp
      character*1 styp,ttyp,ctemp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
_IFN(ibm,vax)
      common/linkup/iiii,jjjj,kkkk,llll
_ENDIF
c
      dimension iv(nmaxgj),jv(nmaxgj),kv(nmaxgj),iconve(nmaxgj)
      dimension ivecc(*)
      dimension v(nmaxe),ga(nmaxg),sta(nmaxg,nmaxc),aa(nmaxa,nmaxa)
      dimension ax(*)
      dimension indvec(*),jndvec(*),kndvec(*),lndvec(*)
c
c     calculation of diagonal matrix elements of a-block, reordering
c     according to magnitude , checking for maximum indeces,
c     determination which 2p1h configurations are included in the
c     nondiagonal a-block.
c     centroids are calculated after the diagonalization of the a-block
c
c
      data styp  /'s'/
      data ttyp  /'t'/
c
      indvx(i,j,k,l) =  lndvec(l) +
     *                  indvec(k+jndvec(j+kndvec(i)))
c
      tim3=cpulft(1)
      dummy=0.0d0
      do 40 i=1,nmaxg
      ga(i)=0.0d0
      iv(i)=0
      jv(i)=0
      kv(i)=0
      do 20 j=1,nmaxc
      sta(i,j)=0.0d0
   20 continue
   40 continue
      ic=0
      do 60 i=1,ntot
      if (idorb(i).ne.isym) go to 60
      ic=ic+1
   60 continue
      jkeep=ic
      ic=0
      icc=ndima(isym)
      do 500 i=1,nocc
      ii=idorb(i)
      iim=imul(ii,isym)
      epsi=epsil(i)
      do 480 ia=nstart,ntot
      iaa=idorb(ia)
      epsa=epsil(ia)
      do 460 ib=ia,ntot
      ibb=idorb(ib)
      if (iim.ne.imul(iaa,ibb)) go to 460
      epsb=epsil(ib)
      ic=ic+1
      iv(ic)=i
      jv(ic)=ia
      kv(ic)=ib
      xtyp(ic)=styp
      if (ia.eq.ib) go to 80
      icc=icc+1
      iv(icc)=i
      jv(icc)=ia
      kv(icc)=ib
      xtyp(icc)=ttyp
   80 continue
      if (ia.ne.ib) go to 240
      ind1= indvx (ia,ia,ia,ia)
      ind2= indvx (ia,ia,i,i)
      ind3= indvx (ia,i,ia,i)
      ga(ic)=epsi-epsa-epsb-(v(ind1)-2.0d0*v(ind2)+v(ind3))
      jc=0
      do 220 ip=1,ntot
      if (idorb(ip).ne.isym) go to 220
      jc=jc+1
_IF(ibm,vax)
      ind = indvor (ia,ip,ia,i)
_ELSE
      call indvor (ia,ip,ia,i)
      ind = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xa=v(ind)
      xc=0.0d0
      del=-2.0d0*epsil(ia)
      do 100 j=1,nocc
      ind = indvx (ia,j,ia,j)
_IF(ibm,vax)
      jnd = indvor (ip,j,i,j)
_ELSE
      call  indvor (ip,j,i,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+(v(ind)*v(jnd))/(2.0d0*epsil(j)+del)
  100 continue
      xd=0.0d0
      if (nocc.eq.1) go to 160
      jend=nocc-1
      do 140 j=1,jend
      ks=j+1
      jj=idorb(j)
      eps=epsil(j)+del
      do 120 k=ks,nocc
      kk=idorb(k)
      if (jj.ne.kk) go to 120
      ind = indvx (ia,k,ia,j)
_IF(ibm,vax)
      jnd = indvor (ip,j,i,k)
      knd = indvor (ip,k,i,j)
_ELSE
      call  indvor (ip,j,i,k)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,k,i,j)
      knd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xd=xd+(v(ind)*(v(jnd)+v(knd)))/(epsil(k)+eps)
  120 continue
  140 continue
  160 continue
      del=epsil(i)-epsil(ia)
      xe=0.0d0
      ibm=imul(ii,iaa)
      do 200 j=1,nocc
      jj=idorb(j)
      eps=del+epsil(j)
      do 180 kc=nstart,ntot
      kcc=idorb(kc)
      if (ibm.ne.imul(jj,kcc)) go to 180
_IF(ibm,vax)
      ind = indvor (ia,i,kc,j)
      jnd = indvor (ia,ip,kc,j)
      knd = indvor (ia,kc,ip,j)
      lnd = indvor (ia,j,kc,i)
_ELSE
      call indvor (ia,i,kc,j)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ia,ip,kc,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,kc,ip,j)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,j,kc,i)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xe=xe-((x2*(-2.0d0*x1+x4)+x3*(x1+x4)))/(eps-epsil(kc))
  180 continue
  200 continue
      sta(ic,jc)=xa+xc+xd+xe
  220 continue
      go to 460
  240 continue
      ind1= indvx (ib,ib,ia,ia)
      ind2= indvx (ib,ia,ib,ia)
      ind3= indvx (ib,ib,i,i)
      ind4= indvx (ib,i,ib,i)
      ind5= indvx (ia,ia,i,i)
      ind = indvx (ia,i,ia,i)
      ga(ic)=epsi-epsa-epsb-(v(ind1)+v(ind2)-v(ind3)+0.5d0*(v(ind4)+
     *v(ind))-v(ind5))
      ga(icc)=epsi-epsa-epsb-(v(ind1)-v(ind2)-v(ind3)+1.5d0*(v(ind4)+
     *v(ind))-v(ind5))
      jc=0
      do 440 ip=1,ntot
      if (idorb(ip).ne.isym) go to 440
      jc=jc+1
_IF(ibm,vax)
      ind1= indvor (ib,i,ia,ip)
      ind2= indvor (ib,ip,ia,i)
_ELSE
      call indvor (ib,i,ia,ip)
      ind1=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,ip,ia,i)
      ind2=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xa=v(ind1)
      xb=v(ind2)
      xc=0.0d0
      del=-epsil(ia)-epsil(ib)
      if (iaa.ne.ibb) go to 280
      do 260 j=1,nocc
      ind = indvx (ib,j,ia,j)
_IF(ibm,vax)
      jnd = indvor (ip,j,i,j)
_ELSE
      call  indvor (ip,j,i,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xc=xc+(v(ind)*v(jnd))/(2.0d0*epsil(j)+del)
  260 continue
      xc=xc+xc
  280 continue
      xd=0.0d0
      xe=0.0d0
      if (nocc.eq.1) go to 340
      jend=nocc-1
      do 320 j=1,jend
      ks=j+1
      jj=idorb(j)
      eps=epsil(j)+del
      do 300 k=ks,nocc
      kk=idorb(k)
      if(iim.ne.imul(jj,kk)) go to 300
      ind = indvx (ib,k,ia,j)
      jnd = indvx (ib,j,ia,k)
_IF(ibm,vax)
      knd = indvor (ip,j,i,k)
      lnd = indvor (ip,k,i,j)
_ELSE
      call indvor (ip,j,i,k)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ip,k,i,j)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      xx=1.0d0/(eps+epsil(k))
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xd=xd+((x1+x2)*(x3+x4))*xx
      xe=xe+((x1-x2)*(x3-x4))*xx
  300 continue
  320 continue
  340 continue
      del=epsil(i)-epsil(ib)
      xf=0.0d0
      xg=0.0d0
      ibm=imul(ibb,ii)
      do 380 j=1,nocc
      jj=idorb(j)
      eps=del+epsil(j)
      do 360 kc=nstart,ntot
      kcc=idorb(kc)
      if (imul(jj,kcc).ne.ibm) go to 360
_IF(ibm,vax)
      ind = indvor (ib,i,kc,j)
      jnd = indvor (ia,ip,kc,j)
      knd = indvor (kc,ia,ip,j)
      lnd = indvor (kc,i,ib,j)
_ELSE
      call indvor (ib,i,kc,j)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ia,ip,kc,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (kc,ia,ip,j)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (kc,i,ib,j)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xx=1.0d0/(eps-epsil(kc))
      xf=xf-xx*(x2*(-2.0d0*x1+x4)+x3*(x1+x4))
      xg=xg-xx*(x2*(-2.0d0*x1+x4)+x3*(x1-x4))
  360 continue
  380 continue
      del=epsil(i)-epsil(ia)
      xh=0.0d0
      xi=0.0d0
      ibm=imul(iaa,ii)
      do 420 j=1,nocc
      jj=idorb(j)
      eps=del+epsil(j)
      do 400 kc=nstart,ntot
      kcc=idorb(kc)
      if (imul(jj,kcc).ne.ibm) go to 400
_IF(ibm,vax)
      ind = indvor (ia,i,kc,j)
      jnd = indvor (ib,ip,kc,j)
      knd = indvor (ib,kc,ip,j)
      lnd = indvor (ia,j,kc,i)
_ELSE
      call indvor (ia,i,kc,j)
      ind = indvx(iiii,jjjj,kkkk,llll)
      call  indvor (ib,ip,kc,j)
      jnd = indvx(iiii,jjjj,kkkk,llll)
      call indvor (ib,kc,ip,j)
      knd=indvx(iiii,jjjj,kkkk,llll)
      call indvor (ia,j,kc,i)
      lnd=indvx(iiii,jjjj,kkkk,llll)
_ENDIF
      x1=v(ind)
      x2=v(jnd)
      x3=v(knd)
      x4=v(lnd)
      xx=1.0d0/(eps-epsil(kc))
      xh=xh-xx*(x2*(-2.0d0*x1+x4)+x3*(x1+x4))
      xi=xi-xx*(x2*(2.0d0*x1-x4)+x3*(-x1+x4))
  400 continue
  420 continue
      sta(ic,jc)=fac*(xa+xb+xc+xd+xf+xh)
      sta(icc,jc)=fad*(xa-xb+xe+xg+xi)
  440 continue
  460 continue
  480 continue
  500 continue
      kdim=ndima(isym)
      ldim=ndimat(isym)
c
c     calculation of the selection criterion 'sum over i' of the
c     absolute value of the residues ( not the amplitudes ) of the
c     self-energy (sta(ip,i)**2) divided by the poles (ga(ip)) minus
c     the orbital energies (epsil(i)) where the summation is over all
c     orbitals i with -50 ev <= epsil(i)<=0 ev.
c
      crits=-50.0d0/conver
      crite=0.0d0
      ic=-1
      do 520 i=1,ntot
      if (idorb(i).ne.isym) go to 520
      ic=ic+1
      istart=i
      if (epsil(i).gt.crits) go to 540
  520 continue
  540 continue
      do 560 i=1,ntot
      if (idorb(i).ne.isym) go to 560
      if (epsil(i).gt.crite) go to 580
      iend=i
  560 continue
  580 continue
      do 620 ip=1,ldim
      xa=ga(ip)
      x=0.0d0
      jc=ic
      do 600 i=istart,iend
      if (idorb(i).ne.isym) go to 600
      jc=jc+1
      x=x+ dabs((sta(ip,jc)**2)/(xa-epsil(i)))
  600 continue
      ax(ip)=x
  620 continue
      ic=0
      do 640 j=1,kdim
      ia=jv(j)
      ib=kv(j)
      if (ia.eq.ib) go to 640
      i=iv(j)
      ic=ic+1
      icc=kdim+ic
c     if (ax(ic).ge.ax(icc)) go to 640
c     ax(ic)=ax(icc)
c ps jan 93 changed ic to j on next two lines
      if (ax(j).ge.ax(icc)) go to 640
      ax(j)=ax(icc)
  640 continue
c
c     ordering of diagonal elements
c
      iend=kdim-1
      do 700 i=1,iend
      js=i+1
      do 680 j=js,kdim
      if (ax(i).ge.ax(j)) go to 680
      temp=ga(i)
      ga(i)=ga(j)
      ga(j)=temp
      temp=ax(i)
      ax(i)=ax(j)
      ax(j)=temp
      it=iv(i)
      iv(i)=iv(j)
      iv(j)=it
      it=jv(i)
      jv(i)=jv(j)
      jv(j)=it
      it=kv(i)
      kv(i)=kv(j)
      kv(j)=it
      ctemp=xtyp(i)
      xtyp(i)=xtyp(j)
      xtyp(j)=ctemp
      do 660 k=1,jkeep
      temp=sta(i,k)
      sta(i,k)=sta(j,k)
      sta(j,k)=temp
  660 continue
  680 continue
  700 continue
      ks=kdim+1
      ic=0
      do 780 j=1,kdim
      ia=jv(j)
      ib=kv(j)
      if (ia.eq.ib) go to 780
      ic=ic+1
      i=iv(j)
      icc=kdim+ic
      do 760 k=ks,ldim
      if (i.eq.iv(k).and.ia.eq.jv(k).and.ib.eq.kv(k)) go to 720
      go to 760
  720 continue
      if (icc.eq.k) go to 780
      it=iv(k)
      iv(k)=iv(icc)
      iv(icc)=it
      it=jv(k)
      jv(k)=jv(icc)
      jv(icc)=it
      it=kv(k)
      kv(k)=kv(icc)
      kv(icc)=it
      ctemp=xtyp(k)
      xtyp(k)=xtyp(icc)
      xtyp(icc)=ctemp
      temp=ga(k)
      ga(k)=ga(icc)
      ga(icc)=temp
      temp=ax(k)
      ax(k)=ax(icc)
      ax(icc)=temp
      do 740 l=1,jkeep
      temp=sta(k,l)
      sta(k,l)=sta(icc,l)
      sta(icc,l)=temp
  740 continue
      go to 780
  760 continue
  780 continue
c
c     second reordering to separate those configurations which have
c     virtual orbital indeces larger than the number which is going
c     to be included in the nondiagonal a-block from the other ones.
c     the blocks ia.le.ib and ia.lt.ib should be kept separately.
c
      if (ntotx.eq.ntot) go to 1200
      icc=0
      iend=kdim-1
      do 800 i=1,kdim
      if (kv(i).le.ntotx) go to 800
      icc=icc+1
  800 continue
      ic=kdim-icc
      jstart=kdim+1
      jend=ldim-1
      jcc=0
      do 820 j=jstart,ldim
      if (kv(j).le.ntotx) go to 820
      jcc=jcc+1
  820 continue
      jc=ndimat(isym)-ndima(isym)-jcc
      do 880 i=1,iend
      if (kv(i).le.ntotx) go to 880
      if (i.gt.ic) go to 900
      js=i+1
      do 860 j=js,kdim
      if (kv(j).gt.ntotx) go to 860
      it=iv(i)
      iv(i)=iv(j)
      iv(j)=it
      it=jv(i)
      jv(i)=jv(j)
      jv(j)=it
      it=kv(i)
      kv(i)=kv(j)
      kv(j)=it
      ctemp=xtyp(i)
      xtyp(i)=xtyp(j)
      xtyp(j)=ctemp
      temp=ga(i)
      ga(i)=ga(j)
      ga(j)=temp
      do 840 l=1,jkeep
      temp=sta(i,l)
      sta(i,l)=sta(j,l)
      sta(j,l)=temp
  840 continue
      go to 880
  860 continue
  880 continue
  900 continue
      do 960 i=jstart,jend
      if (kv(i).le.ntotx) go to 960
      if (i.gt.(jc+kdim)) go to 980
      js=i+1
      do 940 j=js,ldim
      if (kv(j).gt.ntotx) go to 940
      it=iv(i)
      iv(i)=iv(j)
      iv(j)=it
      it=jv(i)
      jv(i)=jv(j)
      jv(j)=it
      it=kv(i)
      kv(i)=kv(j)
      kv(j)=it
      ctemp=xtyp(i)
      xtyp(i)=xtyp(j)
      xtyp(j)=ctemp
      temp=ga(i)
      ga(i)=ga(j)
      ga(j)=temp
      do 920 l=1,jkeep
      temp=sta(i,l)
      sta(i,l)=sta(j,l)
      sta(j,l)=temp
  920 continue
      go to 960
  940 continue
  960 continue
  980 continue
c
c     reorder the part with indeces larger than the number which is
c     going to be included in the nondiagonal a-block according to the
c     magnitude of the terms in the subblock ia.le.ib. the other
c     blocks are kept in order.
c
      if (ic.eq.(kdim-1)) go to 1140
      iend=kdim-1
      istart=ic+1
      do 1040 i=istart,iend
      js=i+1
      do 1020 j=js,kdim
      if (ga(i).ge.ga(j)) go to 1020
      temp=ga(i)
      ga(i)=ga(j)
      ga(j)=temp
      do 1000 l=1,jkeep
      temp=sta(i,l)
      sta(i,l)=sta(j,l)
      sta(j,l)=temp
 1000 continue
      it=iv(i)
      iv(i)=iv(j)
      iv(j)=it
      it=jv(i)
      jv(i)=jv(j)
      jv(j)=it
      it=kv(i)
      kv(i)=kv(j)
      kv(j)=it
      ctemp=xtyp(i)
      xtyp(i)=xtyp(j)
      xtyp(j)=ctemp
 1020 continue
 1040 continue
      ks=kdim+jc+1
      kc=0
      do 1120 j=istart,kdim
      ia=jv(j)
      ib=kv(j)
      if (ia.eq.ib) go to 1120
      kc=kc+1
      i=iv(j)
      kcc=kdim+kc+jc
      do 1100 k=ks,ldim
      if (i.eq.iv(k).and.ia.eq.jv(k).and.ib.eq.kv(k)) go to 1060
      go to 1100
 1060 continue
      if (kcc.eq.k) go to 1120
      it=iv(k)
      iv(k)=iv(kcc)
      iv(kcc)=it
      it=jv(k)
      jv(k)=jv(kcc)
      jv(kcc)=it
      it=kv(k)
      kv(k)=kv(kcc)
      kv(kcc)=it
      ctemp=xtyp(k)
      xtyp(k)=xtyp(kcc)
      xtyp(kcc)=ctemp
      temp=ga(k)
      ga(k)=ga(kcc)
      ga(kcc)=temp
      do 1080 l=1,jkeep
      temp=sta(k,l)
      sta(k,l)=sta(kcc,l)
      sta(kcc,l)=temp
 1080 continue
      go to 1120
 1100 continue
 1120 continue
 1140 continue
 1200 continue
      if (tsele.eq.ssele) go to 1500
      if (ntotx.eq.ntot.and.ldim.le.nmaxj) go to 1380
      ic=0
      jc=0
      do 1220 i=1,kdim
      if (kv(i).gt.ntotx) go to 1240
      ic=ic+1
      jc=jc+1
      if (jv(i).ne.kv(i)) ic=ic+1
      if (ic.gt.(nmaxj-2)) go to 1240
 1220 continue
 1240 continue
      kc=ic-jc
      do 1260 k=1,kc
      iv(jc+k)=iv(kdim+k)
      jv(jc+k)=jv(kdim+k)
      kv(jc+k)=kv(kdim+k)
 1260 continue
      istart=jc+1
      in=kdim-jc
      jstart=kdim+kc+1
      jn=ldim-jstart+1
      kn=in+jn
      do 1300 i=istart,kdim
      ga(i-jc)=ga(i)
      do 1280 j=1,jkeep
      sta(i-jc,j)=sta(i,j)
 1280 continue
 1300 continue
      do 1340 j=jstart,ldim
      jj=j-kdim-kc+in
      ga(jj)=ga(j)
      do 1320 l=1,jkeep
      sta(jj,l)=sta(j,l)
 1320 continue
 1340 continue
      do 1360 i=1,ic
      iconve(i)=i
 1360 continue
      write(nf15)jc,kc,ic,(iv(i),i=1,ic),(jv(i),i=1,ic),(kv(i),i=1,ic),
     *           (xtyp(i),i=1,ic),(iconve(i),i=1,ic)
c
c     structure of data set on unit 15 for block a:
c
c     1st record:
c                jc=number of configurations kept in the block ia.le.ib
c                kc=number of configurations kept in the block ia.lt.ib
c                ic=total number of configurations   ic=jc+kc
c                iv(i),jv(i),kv(i) are the i,ia,ib indeces of the 2p1h
c                                  configurations which are kept.
c                xtyp(i):type of spin coupling : s or t
c                iconve(i): numbering of configurations
c
c     2nd record:
c                in=number of diagonal elements and their strengths in
c                   the block ia.le.ib, which must later be used in the
c                   calculation of centroids.
c                jn=number of diagonal elements and their strengths in
c                   the block ia.lt.ib, which must later be used in the
c                   the calculation of centroids.
c                kn=in+jn
c                list of diagonal elements of the a-block and their
c                strengths.
c
c     note: all elements which are going to be included in the
c     nondiagonal a-block have been deleted.
c
      write(nf15)in,jn,kn,jkeep,(ga(i),i=1,kn),((sta(i,j),i=1,kn),j=1,
     *jkeep)
      go to 1420
 1380 continue
      kc=ldim-kdim
      do 1400 i=1,ldim
      iconve(i)=i
 1400 continue
      write (nf15) kdim,kc,ldim,(iv(i),i=1,ldim),(jv(i),i=1,ldim),
     *      (kv(i),i=1,ldim),(xtyp(i),i=1,ldim),(iconve(i),i=1,ldim)
      write (nf15) (dummy,i=1,10)
c
c     the diagonal elements are not needed if all configurations are
c     included in the nondiagonal a-block.
c
 1420 continue
      go to 2000
c
c
c     data handling for the case of pole search solution to the
c     dyson equation
c
c
 1500 continue
c
      call diagax (kdim,ldim,jkeep,
     * nmaxa,nmaxc,nmaxg,nmaxgj,nmaxj,iv,jv,kv,iconve,
     * ivecc,  ivecc(nmaxj+1), ivecc(nmaxj+nmaxj+1),
     * ga,sta,aa,iwr)
c
 2000 continue
c
      time3(isym)=cpulft(1)-tim3
      return
      end
_IF1(i)@process opt(1)
      subroutine patgaa (jkeep,jkeepx,
     *         nmaxb,nmaxbg,nmaxc,nmaxg,nmaxi,
     *         iconve,iveca,g,gvec,aux,veca,vecb,vecc,axx,iwr)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iconve,iveca
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
      dimension iconve(nmaxbg),iveca(8,nmaxi,*)
      dimension g(nmaxb,nmaxb),gvec(nmaxb),aux(nmaxc,nmaxb)
      dimension veca(nmaxbg),vecb(nmaxbg,nmaxc),vecc(nmaxbg,nmaxc)
      dimension axx(nmaxb,nmaxi)
c
c     this subroutine analyzes the eigenvalues and eigenvectors of the
c     gammaa matrix for a with respect to the subparts of the given
c     symmetry and sorts them into different files which are written on
c     file 17. simultaneously the new dimensions are calculated as well
c     as a new index vector which connects the previous (i,ia,ib) index
c     with the symmetry subpart index.
c
c     this version is used for the pole search procedure.
c
    1 format (/1x,104('-')/)
    2 format (/)
    6 format (//1x,'list of poles of the self-energy'/
     + 1x,'any negative energy solutions arising from linear',
     +    ' dependencies in the basis are deleted later'/
     + 1x,'if they are beyond -10ev'//(10(1x,f11.5)))
    8 format (1x,
     +   'poles of affinity block of self-energy for symmetry ',i3,
     +   ' and subpart ',i3/
     +   1x,'only the lowest 100 poles are listed'//)
    9 format (10(1x,f11.5))
   10 format (/1x,'residues of affinity block of the self-energy',
     +            ' for orbital ',i3/)
c
      tim7=cpulft(1)
      jend=nsymg(isym)
c
      do 50 i=1,ndim
      sigdia(i)=gvec(i)*conver
   50 continue
      write (iwr,6) (sigdia(i),i=1,ndim)
      write (iwr,1)
      crit=3.0d-04
      do 60 i=1,nmaxb
      sigdia(i)=0.0d0
      do 60 j=1,nmaxb
      g(j,i)=0.0d0
   60 continue
c
c     loop over all subparts of symmetry isym.
c     analyze only the solutions which are fully taken into account in
c     the dyson equation and not those which are used only to form
c     centroids.
c
      do 200 j=1,jend
      ic=0
      kend=isgnum(isym,j)
c     loop over all solutions
      do 160 i=1,ndim
      sum=0.0d0
      do 100 k=1,kend
      ind=jsgjnd(isym,j,k)
      sum=sum+ dabs(aux(ind,i))
  100 continue
      if (sum.lt.crit) go to 160
      if (gvec(i).lt.-0.4d0) go to 160
      ic=ic+1
      iveca(isym,j,ic)=i
  160 continue
      idima(isym,j)=ic
  200 continue
      iend=nmaxb+nmaxg
      do 250 i=1,iend
      veca(i)=0.0d0
      do 250 j=1,nmaxc
      vecb(i,j)=0.0d0
      vecc(i,j)=0.0d0
  250 continue
c
c     get diagonal poles if there are any
c
      if (ndim.eq.ndimat(isym)) go to 600
      istart=nmaxb+1
      iend=istart+npola(isym,1)-1
      ikeep=npola(isym,1)
      read (nf15) kn,ln,(veca(i),i=istart,iend),
     *          ((vecb(i,j),i=istart,iend),j=1,ln)
c
c     sort vecb into the symmetry subelements
c
      kc=0
      do 300 j=1,jend
      kend=isgnum(isym,j)
      do 290 k=1,kend
      ind=jsgjnd(isym,j,k)
      do 280 jj=1,jkeep
      if (jj.eq.ind) go to 260
      go to 280
  260 continue
      do 270 i=istart,iend
      vecc(i,k+kc)=vecb(i,jj)
  270 continue
  280 continue
  290 continue
      kc=kc+kend
  300 continue
c
c     sort aux into the symmetry subelements, store in vecc.
c
      kc=0
      do 350 j=1,jend
      kend=isgnum(isym,j)
      do 340 k=1,kend
      ind=jsgjnd(isym,j,k)
      do 330 jj=1,jkeepx
      if (jj.eq.ind) go to 310
      go to 330
  310 continue
      do 320 i=1,ndim
      vecc(i,k+kc)=aux(jj,i)
  320 continue
  330 continue
  340 continue
      kc=kc+kend
  350 continue
      kpx=ndim+npola(isym,1)
      do 360 i=1,kpx
      iconve(i)=i
  360 continue
c
c     delete weak poles from aux part of vecc and from gvec and sort
c     energies into axx.
c     delete negative affinity poles which may arise from linear
c     dependencies in the basis if the energies are beyond -10 ev.
c
      kc=0
      do 400 j=1,jend
      ic=0
      kend=isgnum(isym,j)
      do 390 i=1,ndim
      sum=0.0d0
      do 370 k=1,kend
      sum=sum+ dabs(vecc(i,k+kc))
  370 continue
      if (sum.lt.crit) go to 390
      if (gvec(i).lt.-0.40d0) go to 390
      ic=ic+1
      do 380 k=1,kend
      vecc(ic,k+kc)=vecc(i,k+kc)
  380 continue
      axx(ic,j)=gvec(i)
  390 continue
      kc=kc+kend
      npola(isym,j)=ic
  400 continue
c
c     compress vecc
c
      kc=0
      do 430 j=1,jend
      kend=isgnum(isym,j)
      ix=npola(isym,j)+1
      iz=ix+ikeep-1
      is=nmaxb
      do 420 i=ix,iz
      is=is+1
      do 410 k=1,kend
      vecc(i,k+kc)=vecc(is,k+kc)
  410 continue
  420 continue
      kc=kc+kend
      npola(isym,j)=npola(isym,j)+ikeep
  430 continue
c
c     sort energies of axx and veca into vecb.
c
c     change sign to make affinity poles negative.
c
      do 470 j=1,jend
      ix=npola(isym,j)-ikeep
      do 450 i=1,ix
      vecb(i,j)=-axx(i,j)
  450 continue
      ix=ix+1
      iz=ix+ikeep-1
      is=nmaxb
      do 460 i=ix,iz
      is=is+1
      vecb(i,j)=-veca(is)
  460 continue
  470 continue
c
c     reorder poles
c
      kc=0
      do 570 j=1,jend
      kend=isgnum(isym,j)
      ix=npola(isym,j)
      ie=ix-1
      jx=ix-ikeep
      do 510 m=1,jx
      iconve(m)=iveca(isym,j,m)
  510 continue
      ms=jx+1
      me=ikeep+jx
      mc=ndim
      do 520 m=ms,me
      mc=mc+1
      iconve(m)=mc
  520 continue
      do 550 i=1,ie
      ls=i+1
      do 540 l=ls,ix
      if (vecb(i,j).le.vecb(l,j)) go to 540
      temp=vecb(i,j)
      vecb(i,j)=vecb(l,j)
      vecb(l,j)=temp
      do 530 k=1,kend
      temp=vecc(i,k+kc)
      vecc(i,k+kc)=vecc(l,k+kc)
      vecc(l,k+kc)=temp
  530 continue
      it=iconve(i)
      iconve(i)=iconve(l)
      iconve(l)=it
  540 continue
  550 continue
      do 560 i=1,ix
      if (vecb(i,j).le.0.4d0) go to 560
      ik=i
      go to 562
  560 continue
      go to 565
  562 continue
      ik=ik-1
      npola(isym,j)=ik
      ix=ik
  565 continue
c
c     store energies and amplitudes on file 17 one record for each
c     symmetry subspecies
c
      write (nf17) isym,j,ix,kend,(vecb(i,j),i=1,ix),
     *           ((vecc(i,k+kc),i=1,ix),k=1,kend),(iconve(i),i=1,ix)
c
c     structure of data set on file 17
c
c          isym :        symmetry species
c          j:            subspecies of this symmetry species
c          ix :          number of affinity poles for isym,j
c          kend :        number of orbitals for isym,j
c          vecb(i,j) :   energies i=1,ix
c          vecc(i,k) :   amplitudes i=1,ix,k=1,kend
c          iconve(i) :   configuration identification vector. its value
c                        is smaller or equal to the number of
c                        diagonalized poles or larger than this number
c                        independent of position because of reordering.
c                        missing values refer to deleted weak poles or
c                        to those appearing in another subspecies
c                        (only for diagonalized poles).
c
      kc=kc+kend
  570 continue
      go to 900
  600 continue
c
c     treatment of case no diagonal poles are included.
c
c     sort aux into symmetry subelements, store in vecc.
c
      read (nf15) dummy
      kc=0
      do 650 j=1,jend
      kend=isgnum(isym,j)
      do 640 k=1,kend
      ind=jsgjnd(isym,j,k)
      do 630 jj=1,jkeep
      if (jj.eq.ind) go to 610
      go to 630
  610 continue
      do 620 i=1,nmaxb
      vecc(i,k+kc)=aux(jj,i)
  620 continue
  630 continue
  640 continue
      kc=kc+kend
  650 continue
c
c     delete weak poles from vecc and from gvec and sort energies into
c     axx.
c     delete negative affinity poles which may arise from linear
c     dependencies in the basis if the energies are beyond -10ev.
c
      kc=0
      do 700 j=1,jend
      ic=0
      kend=isgnum(isym,j)
      do 680 i=1,ndim
      sum=0.0d0
      do 660 k=1,kend
      sum=sum+ dabs(vecc(i,k+kc))
  660 continue
      if (sum.lt.crit) go to 680
      if (gvec(i).lt.-0.40d0) go to 680
      ic=ic+1
      do 670 k=1,kend
      vecc(ic,k+kc)=vecc(i,k+kc)
  670 continue
      axx(ic,j)=gvec(i)
  680 continue
      kc=kc+kend
      npola(isym,j)=ic
  700 continue
      ikeep=0
c
c     sort energies of axx into vecb
c
c     change sign to make affinity poles negative
c
      do 730 j=1,jend
      ix=npola(isym,j)
      do 710 i=1,ix
      vecb(i,j)=-axx(i,j)
  710 continue
  730 continue
c
c     reorder poles and write poles and amplitudes on file 17
c
      kc=0
      do 800 j=1,jend
      kend=isgnum(isym,j)
      ix=npola(isym,j)
      do 750 m=1,ix
      iconve(m)=iveca(isym,j,m)
  750 continue
      ie=ix-1
      do 780 i=1,ie
      ls=i+1
      do 770 l=ls,ix
      if (vecb(i,j).le.vecb(l,j)) go to 770
      temp=vecb(i,j)
      vecb(i,j)=vecb(l,j)
      vecb(l,j)=temp
      do 760 k=1,kend
      temp=vecc(i,k+kc)
      vecc(i,k+kc)=vecc(l,k+kc)
      vecc(l,k+kc)=temp
  760 continue
      it=iconve(i)
      iconve(i)=iconve(l)
      iconve(l)=it
  770 continue
  780 continue
      write (nf17) isym,j,ix,kend,(vecb(i,j),i=1,ix),
     *           ((vecc(i,k+kc),i=1,ix),k=1,kend),(iconve(i),i=1,ix)
      kc=kc+kend
  800 continue
  900 continue
c
c     output for affinity poles and strengths
c
      kc=0
      do 950 j=1,jend
      if (j.ne.1) write (iwr,2)
      kend=isgnum(isym,j)
      is=1
      ie=npola(isym,j)
      if (ie.gt.100) is=ie-99
      write (iwr,8) isym,j
      do 910 i=is,ie
      gvec(i)=vecb(i,j)*conver
  910 continue
      write (iwr,9) (gvec(i),i=is,ie)
      do 930 k=1,kend
      ind=isgind(isym,j,k)
      write (iwr,10) ind
      do 920 l=is,ie
      veca(l)=cons*(vecc(l,k+kc)**2)
  920 continue
      write (iwr,9) (veca(l),l=is,ie)
  930 continue
      kc=kc+kend
  950 continue
      do 980 i=1,nmaxb
      gvec(i)=0.0d0
      sigdia(i)=0.0d0
      do 960 j=1,nmaxc
      aux(j,i)=0.0d0
  960 continue
      do 970 k=1,nmaxb
      g(k,i)=0.0d0
  970 continue
  980 continue
c
      time7(isym)=cpulft(1)-tim7
c
      return
      end
_IF1(i)@process opt(1)
      subroutine patgab (jkeep,jkeepx,
     * nmaxb,nmaxbg,nmaxc,nmaxg,nmaxi,iconve,ivecb,
     * g,gvec,aux,veca,vecb,vecc,axx,iwr)
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iconve,ivecb
      character*1 xtyp
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
c
      dimension iconve(nmaxbg),ivecb(8,nmaxi,*)
      dimension g(nmaxb,nmaxb),gvec(nmaxb),aux(nmaxc,nmaxb)
      dimension veca(nmaxbg),vecb(nmaxbg,nmaxc),vecc(nmaxbg,nmaxc)
      dimension axx(nmaxb,nmaxi)
c
c     this subroutine analyzes the eigenvalues and eigenvectors of the
c     gammaa matrix for b with respect to the subparts of the given
c     symmetry and sorts them into different files which are written on
c     file 13. simultaneously the new dimensions are calculated as well
c     as a new index vector which connects the previous (i,j,ia) index
c     with the symmetry subpart index.
c
c     this version is used for the pole search procedure.
c
    1 format (/1x,104('-')/)
    2 format (/)
    6 format (//1x,'list of poles of the self-energy'//(10(1x,f11.5)))
    8 format (1x,
     +  'poles of ionization block of self-energy for symmetry ',i3,
     +  ' and subpart ',i3/
     +   1x,'only the lowest 100 poles are listed'//)
    9 format (10(1x,f11.5))
   10 format (/1x,
     + 'residues of ionization block of the self-energy for orbital ',
     +  i3/)
c
      tim10=cpulft(1)
      jend=nsymg(isym)
c
      write (iwr,1)
      do 50 i=1,ndim
      sigdia(i)=gvec(i)*conver
   50 continue
      write (iwr,6) (sigdia(i),i=1,ndim)
      write (iwr,1)
      crit=3.0d-04
      do 60 i=1,nmaxb
      sigdia(i)=0.0d0
      do 60 j=1,nmaxb
      g(j,i)=0.0d0
   60 continue
c
c     loop over all subparts of symmetry isym
c
      do 200 j=1,jend
      ic=0
      kend=isgnum(isym,j)
c     loop over all solutions
      do 160 i=1,ndim
      sum=0.0d0
      do 100 k=1,kend
      ind=jsgjnd(isym,j,k)
      sum=sum+ dabs(aux(ind,i))
  100 continue
      if (sum.lt.crit) go to 160
      ic=ic+1
      ivecb(isym,j,ic)=i
  160 continue
      jdimb(isym,j)=ic
      idimb(isym,j)=ic
  200 continue
      iend=nmaxb+nmaxg
      do 250 i=1,iend
      veca(i)=0.0d0
      do 250 j=1,nmaxc
      vecb(i,j)=0.0d0
      vecc(i,j)=0.0d0
  250 continue
c
c     get diagonal poles if there are any
c
      if (ndim.eq.ndimbt(isym)) go to 600
      istart=nmaxb+1
      iend=istart+npolb(isym,1)-1
      ikeep=npolb(isym,1)
      read (nf16) kn,ln,(veca(i),i=istart,iend),
     *          ((vecb(i,j),i=istart,iend),j=1,ln)
c
c     sort vecb into the symmetry subelements
c
      kc=0
      do 300 j=1,jend
      kend=isgnum(isym,j)
      do 290 k=1,kend
      ind=jsgjnd(isym,j,k)
      do 280 jj=1,jkeep
      if (jj.eq.ind) go to 260
      go to 280
  260 continue
      do 270 i=istart,iend
      vecc(i,k+kc)=vecb(i,jj)
  270 continue
  280 continue
  290 continue
      kc=kc+kend
  300 continue
c
c     sort aux into the symmetry subelements, store in vecc.
c
      kc=0
      do 350 j=1,jend
      kend=isgnum(isym,j)
      do 340 k=1,kend
      ind=jsgjnd(isym,j,k)
      do 330 jj=1,jkeepx
      if (jj.eq.ind) go to 310
      go to 330
  310 continue
      do 320 i=1,ndim
      vecc(i,k+kc)=aux(jj,i)
  320 continue
  330 continue
  340 continue
      kc=kc+kend
  350 continue
      kpx=ndim+npolb(isym,1)
      do 360 i=1,kpx
      iconve(i)=i
  360 continue
c
c     delete weak poles from aux part of vecc and from gvec and sort
c     energies into axx.
c
      kc=0
      do 400 j=1,jend
      ic=0
      kend=isgnum(isym,j)
      do 390 i=1,ndim
      sum=0.0d0
      do 370 k=1,kend
      sum=sum+ dabs(vecc(i,k+kc))
  370 continue
      if (sum.lt.crit) go to 390
      ic=ic+1
      do 380 k=1,kend
      vecc(ic,k+kc)=vecc(i,k+kc)
  380 continue
      axx(ic,j)=gvec(i)
  390 continue
      kc=kc+kend
      npolb(isym,j)=ic
  400 continue
c
c     compress vecc after deletion of weak poles
c
      kc=0
      do 430 j=1,jend
      kend=isgnum(isym,j)
      ix=npolb(isym,j)+1
      iz=ix+ikeep-1
      is=nmaxb
      do 420 i=ix,iz
      is=is+1
      do 410 k=1,kend
      vecc(i,k+kc)=vecc(is,k+kc)
  410 continue
  420 continue
      kc=kc+kend
      npolb(isym,j)=npolb(isym,j)+ikeep
  430 continue
c
c     sort energies of axx and veca into vecb.
c
c     change sign to make ionization poles positive.
c
      do 470 j=1,jend
      ix=npolb(isym,j)-ikeep
      do 450 i=1,ix
      vecb(i,j)=-axx(i,j)
  450 continue
      ix=ix+1
      iz=ix+ikeep-1
      is=nmaxb
      do 460 i=ix,iz
      is=is+1
      vecb(i,j)=-veca(is)
  460 continue
  470 continue
c
c     reorder poles
c
      kc=0
      do 570 j=1,jend
      kend=isgnum(isym,j)
      ix=npolb(isym,j)
      ie=ix-1
      jx=ix-ikeep
      do 510 m=1,jx
      iconve(m)=ivecb(isym,j,m)
  510 continue
      ms=jx+1
      me=ikeep+jx
      mc=ndim
      do 520 m=ms,me
      mc=mc+1
      iconve(m)=mc
  520 continue
      do 550 i=1,ie
      ls=i+1
      do 540 l=ls,ix
      if (vecb(i,j).le.vecb(l,j)) go to 540
      temp=vecb(i,j)
      vecb(i,j)=vecb(l,j)
      vecb(l,j)=temp
      do 530 k=1,kend
      temp=vecc(i,k+kc)
      vecc(i,k+kc)=vecc(l,k+kc)
      vecc(l,k+kc)=temp
  530 continue
      it=iconve(i)
      iconve(i)=iconve(l)
      iconve(l)=it
  540 continue
  550 continue
c
c     store energies and amplitudes on file 13 one record for each
c     symmetry subspecies
c
      write (nf13) isym,j,ix,kend,(vecb(i,j),i=1,ix),
     *           ((vecc(i,k+kc),i=1,ix),k=1,kend),(iconve(i),i=1,ix)
c
c     structure of data set on file 13
c
c          isym :        symmetry species
c          j:            subspecies of this symmetry species
c          ix :          number of ionization poles for isym,j
c          kend :        number of orbitals for isym,j
c          vecb(i,j) :   energies i=1,ix
c          vecc(i,k) :   amplitudes i=1,ix,k=1,kend
c          iconve(i) :   configuration identification vector. its value
c                        is smaller or equal to the number of
c                        diagonalized poles or larger than this number
c                        independent of position because of reordering.
c                        missing values refer to deleted weak poles or
c                        to those appearing in another subspecies
c                        (only for diagonalized poles).
c
      kc=kc+kend
  570 continue
      go to 900
  600 continue
c
c     treatment of case no diagonal poles are included.
c
c     sort aux into symmetry subelements, store in vecc.
c
      read (nf16) dummy
      kc=0
      do 650 j=1,jend
      kend=isgnum(isym,j)
      do 640 k=1,kend
      ind=jsgjnd(isym,j,k)
      do 630 jj=1,jkeep
      if (jj.eq.ind) go to 610
      go to 630
  610 continue
      do 620 i=1,nmaxb
      vecc(i,k+kc)=aux(jj,i)
  620 continue
  630 continue
  640 continue
      kc=kc+kend
  650 continue
c
c     delete weak poles from vecc and from gvec and sort energies into
c     axx.
c
      kc=0
      do 700 j=1,jend
      ic=0
      kend=isgnum(isym,j)
      do 680 i=1,ndim
      sum=0.0d0
      do 660 k=1,kend
      sum=sum+ dabs(vecc(i,k+kc))
  660 continue
      if (sum.lt.crit) go to 680
      ic=ic+1
      do 670 k=1,kend
      vecc(ic,k+kc)=vecc(i,k+kc)
  670 continue
      axx(ic,j)=gvec(i)
  680 continue
      kc=kc+kend
      npolb(isym,j)=ic
  700 continue
      ikeep=0
c
c     sort energies of axx into vecb
c
c     change sign to make ionization poles positive.
c
      do 730 j=1,jend
      ix=npolb(isym,j)
      do 710 i=1,ix
      vecb(i,j)=-axx(i,j)
  710 continue
  730 continue
c
c     reorder poles and write poles and amplitudes on file 13
c
      kc=0
      do 800 j=1,jend
      kend=isgnum(isym,j)
      ix=npolb(isym,j)
      do 750 m=1,ix
      iconve(m)=ivecb(isym,j,m)
  750 continue
      ie=ix-1
      do 780 i=1,ie
      ls=i+1
      do 770 l=ls,ix
      if (vecb(i,j).le.vecb(l,j)) go to 770
      temp=vecb(i,j)
      vecb(i,j)=vecb(l,j)
      vecb(l,j)=temp
      do 760 k=1,kend
      temp=vecc(i,k+kc)
      vecc(i,k+kc)=vecc(l,k+kc)
      vecc(l,k+kc)=temp
  760 continue
      it=iconve(i)
      iconve(i)=iconve(l)
      iconve(l)=it
  770 continue
  780 continue
      write (nf13) isym,j,ix,kend,(vecb(i,j),i=1,ix),
     *           ((vecc(i,k+kc),i=1,ix),k=1,kend),(iconve(i),i=1,ix)
      kc=kc+kend
  800 continue
  900 continue
c
c     output for ionization poles and strengths
c
      kc=0
      do 950 j=1,jend
      if (j.ne.1) write (iwr,2)
      kend=isgnum(isym,j)
      ix=npolb(isym,j)
      if (ix.gt.100) ix=100
      write (iwr,8) isym,j
      do 910 i=1,ix
      gvec(i)=vecb(i,j)*conver
  910 continue
      write (iwr,9) (gvec(i),i=1,ix)
      do 930 k=1,kend
      ind=isgind(isym,j,k)
      write (iwr,10) ind
      do 920 l=1,ix
      veca(l)=cons*(vecc(l,k+kc)**2)
  920 continue
      write (iwr,9) (veca(l),l=1,ix)
  930 continue
      kc=kc+kend
  950 continue
      do 980 i=1,nmaxb
      gvec(i)=0.0d0
      sigdia(i)=0.0d0
      do 960 j=1,nmaxc
      aux(j,i)=0.0d0
  960 continue
      do 970 k=1,nmaxb
      g(k,i)=0.0d0
  970 continue
  980 continue
c
      time10(isym)=cpulft(1)-tim10
      return
      end
_IFN(ibm,vax)
      subroutine indvor(i,j,k,l)
      implicit REAL  (a-h,o-z)
      common/linkup/ia,ib,ic,id
      ia=max(i,j)
      ib=min(i,j)
      ic=max(k,l)
      id=min(k,l)
      if (ia.gt.ic) return
      if (ia.lt.ic) then
      it=ia
      ia=ic
      ic=it
      it=ib
      ib=id
      id=it
      else
      if (ib.lt.id) then
      it=ib
      ib=id
      id=it
      endif
      endif
      return
      end
_ENDIF
      subroutine ver_tdaf(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/tdaf.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
_ENDEXTRACT
_EXTRACT(outtda,mips4)
      subroutine outtda (jp,isuma,isumb,ncont,
     * nmaxbg,maxbg2,mxbg22,nmaxc,
     * iv,jv,veca,vecb,vecc,bmat,iwr )
c
      implicit REAL  (a-h,o-z)
_IFN1(c)      integer*2 iv,jv
      character*1 xtyp
      character*1 typ1,typ2,typ3
c     character*1 typ4
      character*8    selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,bmvec,ztdaop
INCLUDE(common/sizes)
      common /junk/  sigdia(mxtda4),epsil(mxtda4)
      common /york/  fmul(15),range(2),alimit,blimit
      common /tdaa/  nocc,nvirt,nstart,ntot,nvirtx,ntotx,ninta,nintb,
     *               nsym,imul(10,10),idorb(100),norb,indorb(30)
      common /tdab/ conver,cons,fac,fad,fac1,amvec(10),mvec(10),
     *               ndima(8),ndimb(8),ndimat(8),ndimbt(8),ndimct(8),
     *               isym,ndim,numend,npola(8,6),npolb(8,6),
     *               npolt(8,6),ixxx(8,6),itemp(mxtda4)
      common /junkc/selcta,selctb,xxsel,yysel,zzsel,sselx,ssely,selctc,
     *               selctd,sseld,ssele,tsele,word,ztdaop(4),bmvec(10),
     *               xtyp(mxtda1)
      common /miscop/lsym,nsymg(8),iuse(8,6),isgnum(8,6),
     *               isgind(8,6,mxtda3),jsgjnd(8,6,mxtda3),
     *               idimb(8,6),jdimb(8,6),idima(8,6),ksgnum(8,6),
     *               ndeg(mxtda2),indeg(mxtda2,3),idega
      common /craypk/ timfac,time1,time2,time3(8),time4(8),time5(8),
     *               time6(8),time7(8),time8(8),time9(8),time10(8),
     *               time11,time12,time13(500),time14,icdiag,idimdi(500)
      common /cgraph/ nsolc,ndummy,en(200),ps(200),norbin(200),iccc(100)
      common/ftape/nf12,nf13,nf14,nf15,nf16,nf17,nf18,nf19
c
      dimension iv(nmaxbg),jv(nmaxbg)
      dimension veca(mxbg22),vecb(mxbg22),vecc(mxbg22)
      dimension bmat(nmaxc,maxbg2)
c
    1 format (/1x,104('-')/)
    3 format (1x,'energies and pole strengths in nondiagonal',
     +           ' approximation for symmetry ',i3,
     +           ' and subpart   ',i3//)
    4 format (1x,'solution   ',i3,3x,f9.3,10(1x,f9.6)/(29x,10(1x,f9.6)))
    5 format (1x,'following notation is used : m : orbital index', 
     +           ' i : 2h1p index, a : 2p1h index and c : centroid'//)
    6 format (1x,'large components of eigenvectors in ',
     +           'nondiagonal approximation')
    7 format (1x,'eigenvector',i4,'::',10(1x,a1,i3,': ',f5.2))
   15 format (/1x,'sum of pole strengths of occupied orbital number',
     *i5,'   in the given symmetry :   ',f6.4)
   16 format (1x,'energies of states with a pole strength larger',
     +           ' than 0.005 for any occupied orbital are only listed'
     +            //)
c
      data typ1   /'m'/
      data typ2   /'i'/
      data typ3   /'a'/
c     data typ4   /'c'/
c
c     this subroutine prints the output for the pole search procedure
c
      write(6,1)
      write (iwr,3) isym,jp
      write (iwr,16)
      rewind nf12
      neigen=npolt(isym,jp)
      kend=isgnum(isym,jp)
      keigen=neigen
      do 100 i=1,nmaxc
      sigdia(i)=0.0d0
  100 continue
      do 500 k=1,neigen
      read (nf13) meigen,n,idim,eigen,(vecb(i),i=1,idim)
      veca(k)=eigen
      eigen=eigen*conver
      do 120 i=1,kend
      bmat(i,k)=vecb(i)
  120 continue
      if (eigen.gt.-20.0d0) go to 140
      keigen=keigen-1
      go to 500
  140 continue
      write (nf12) neigen,k,idim,eigen,(vecb(i),i=1,idim)
      do 200 i=1,kend
      vecc(i)=vecb(i)* dabs(vecb(i))
  200 continue
      if ( eigen .le. 5.0d0 .or. eigen .ge. 70.0d0 ) goto 210
      pi = 0.0d0
      vcmax = 0.0d0
      jmax  = 0
      ie = isuma
      if ( isuma .lt. 1 ) ie = 1
      do 201 i=1,ie
      dvci =  dabs(vecc(i))
      pi = pi + dvci
      if ( dvci .le. vcmax ) goto 201
      vcmax = dvci
      jmax  = i
  201 continue
      if ( pi .lt. 0.01d0 ) goto 210
      if ( pi .gt. 0.03d0 .and. vcmax .lt. 0.8d0*pi ) jmax = 0
      nsolc = nsolc + 1
      norbin(nsolc) = 0
      en    (nsolc) = eigen
      ps    (nsolc) = pi
      if ( jmax .ne. 0 ) goto 203
      norbin(nsolc) = 0
      goto 210
  203 ix = isgind(isym,jp,jmax)
      norbin(nsolc) = iccc(ix)
  210 if (eigen.lt.21.0d0) go to 240
      do 220 i=1,isumb
      if ( dabs(vecc(i)).ge.0.005d0) go to 240
  220 continue
      go to 440
  240 continue
      isumx = isumb
      if ( eigen .lt. 4.0d0 ) isumx = kend
      write (iwr,4) k,eigen,(vecc(j),j=1,isumx)
  440 continue
      if (eigen.le.0.0d0) go to 480
      do 460 i=1,isuma
      sigdia(i)=sigdia(i)+ dabs(vecc(i))
  460 continue
  480 continue
  500 continue
      ndif=neigen-keigen
      rewind nf12
      do 520 i=1,isuma
      write (iwr,15) i,sigdia(i)
  520 continue
      if (ncont.eq.2) go to 1100
      write (iwr,1)
      read (nf16) iz,jz,(iv(i),i=1,iz),(jv(i),i=1,jz)
      write (iwr,6)
      write (iwr,5)
      do 1000 k=1,keigen
      kp=k+ndif
      read (nf12) neigen,n,idim,eigen,(vecb(i),i=1,idim)
      do 540 i=1,kend
      vecc(i)=vecb(i)**2
  540 continue
      if ( eigen .lt. 4.0d0 ) go to 580
      do 560 i=1,isumb
      if (vecc(i).gt.0.005d0) go to 580
  560 continue
      go to 1000
  580 continue
      jc=0
      do 800 j=1,idim
      if ( dabs(vecb(j)).lt.0.15d0) go to 800
      jc=jc+1
      if (j.gt.kend) go to 640
      bmvec(jc)=typ1
      mvec(jc)=isgind(isym,jp,j)
      go to 680
  640 continue
      if (j.gt.(iz+kend)) go to 660
      bmvec(jc)=typ3
      mvec(jc)=iv(j-kend)
      go to 680
  660 continue
      bmvec(jc)=typ2
      mvec(jc)=jv(j-iz-kend)
  680 continue
      amvec(jc)=vecb(j)
      if (jc.lt.10) go to 800
      write (iwr,7) kp,(bmvec(l),mvec(l),amvec(l),l=1,jc)
      jc=0
  800 continue
      if (jc.eq.0) go to 820
      if (jc.eq.1.and. dabs(amvec(1)).gt.0.99d0) go to 820
      write (iwr,7) kp,(bmvec(l),mvec(l),amvec(l),l=1,jc)
  820 continue
 1000 continue
 1100 continue
      rewind nf12
c
      return
      end
_ENDEXTRACT
