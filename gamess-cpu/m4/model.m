      subroutine model(rx,ix,index)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
INCLUDE(common/sizes)
      parameter(mxatm3 = maxat3 * 3)
_IFN1(cuf)      integer*2 ix
c
c --- this is the interface driver routine to perform the
c     molecular mechanical calculation to get the energy
c     and gradient contribution due the surrounding atoms
c     to the atoms treated by ab initio method. it generates
c     the z matrix for the ab initio atoms and orient the
c     the entire mechanics coordinates along the ab initio coordinates.
c
c     the nature of the task is controlled by the arguement index.
c
c         index   ...  0   initialisation of the run.
c                          the model molecular topology and the
c                          coordinates are read and other variab
c                          les are initialised and it is written
c                          into the direct access files
c
c                 ...  1   the z matrix is generated for the ab initio
c                          atoms from the mechanics coordinates
c
c                 ...  2   the energy contribution due to the
c                          surroundings is calculated
c
c                 ...  3   the gradient contribution due to the
c                          surroundings is calculated
c
c                 ...  4   the ab initio atoms are minimised by the
c                          molecular mechanics method
c
c                 ...  5   the surroundings are minimised by the
c                          molecular mechanics method
c
c                 ...  6   the second derivatives for the ab initio
c                          atoms are calculated by molecular
c                          mechanics method and transformed to
c                          the internal coordinates
c
c                 ...  7   the second derivatives for the ab initio
c                          atoms are calculated and the cartesian
c                          derivatives are output
c
c                  ... 8   the surrounding atoms are shaken by
c                          dynamics
c
c                  ... 9   Tidy up !!!
c
c     ----- the memory is dynamically allocated and to understand
c           the location of various arrays consult the routine
c           lcmem0 and lcmem1 -----
c
INCLUDE(common/phycon)
INCLUDE(common/modj)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
      common/czmat/ianz(maxnz),iz(maxnz,4),bl(maxnz),alpha(maxnz),
     & beta(maxnz),lbl(maxnz),lalpha(maxnz),lbeta(maxnz),nz,nvar
INCLUDE(common/cntl1)
INCLUDE(common/finish)
INCLUDE(common/restar)
INCLUDE(common/funct)
INCLUDE(common/g80nb)
      common/junk/de(3*maxat),cg80(3*maxat)
INCLUDE(common/infob)
INCLUDE(common/moda)
INCLUDE(common/modf)
INCLUDE(common/modmin)
INCLUDE(common/modmd)
INCLUDE(common/modfov)
INCLUDE(common/mod2d)
      common/memla/i02,i04,i06,i08,i10,i12,i14,i16,i18,i20,
     +             i22,i24,i26,i28,i30,i32,i34,i36,i38,i40,
     +             i42,i44,i46,i48,i50,i52,i54,i56,i58,i60,
     +             i62,i64,i66,i68,i70,i72,i74,i76,i78,i80
      common/memlb/l05,l10,l15,l20,l25,l30,l35,l40,l45,l50,
     +             l55,l60,l65,l70,l75,l80,l85,l90,l95,l100,
     +             k02,k04,k06,k08,k10,k12,k14,k16,k18,k20,
     +             k22,k24,k26,k28,k30,k32,k34,k36,k38,k40
INCLUDE(common/runhed)
INCLUDE(common/nbterm)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      save nr,nr3,len173
      dimension rx(*),ix(*),ene(30)
      dimension zdnam(3)
      character *7 fnm
      character *5 snm
      data fnm/"model.m"/
      data snm/"model"/
      data atocal/627.5251387752d0/, chgfac/18.2223d+00/
      data m51/51/
      data zdnam/'de/dx','de/dy','de/dz'/
c
      istep = index + 1
      nav = lenwrd()
      go to (1100,1100,1300,1300,1400,1400,1500,1500,1600,1700) , istep
c ...
c     ----- istep = 1  initialisation -----
c
c     ----- nword is the number of i*2 words per double word. if
c           the i*2 arrays are abandoned then nword must be suitably
c           defined for sucessful memory management. memb is the
c           maximum memory limit -----
c ...
_IFN1(cuf) 1100 nword = 4
_IF1(cuf) 1100 nword = 1
c ...
c     ----- read the necessary data to initiate the run -----
c ...
      write(iwr,9699)
      if ( kform .gt. 0 ) then
            write(iwr,9700)
      else
            write(iwr,9702)
      end if
      call rdparm(rx,ix)
      dxm  = accin(3)
      dx0  = accin(3) / 40.0d0
      drms = accin(2) /  6.0d0
      dele = cccnv
c ...
c     ----- print data characterizing the md-run -----
c ...
      nrp  = natom
      nrpt = npm * nrp
      nr   = nrpt + nsm * nram
      write(iwr,9600) (ititl(i),i=1,20)
      if ( idiel .eq. 0 ) then
           write(iwr,9620)
      else
           write(iwr,9622) dielc
      end if
      if ( ntmin .ne. 0 ) then
           write(iwr,9630)imcyc
      else
           write(iwr,9632)jm
      end if
      if ( ntmin .eq. 3 ) then
           write(iwr,9640)
           if ( init .eq. 0 ) then
                write(iwr,9642)tempi
           else
                if ( ntx .eq. 0 ) write(iwr,9644)
                if ( ntx .eq. 1 ) write(iwr,9646)
           end if
           if ( ntb ) 1110, 1120, 1130
 1110      write(iwr,9650)
           go to 1135
 1120      write(iwr,9652)
           go to 1135
 1130      if ( bbeta .eq. 90.0d0 ) then
                write(iwr,9654)
           else
                write(iwr,9656)
           endif
 1135      continue
           write(iwr,9660)( box(i),i=1,3)
           if ( opress ) write(iwr,9662)
           if ( otemp  ) then
                write(iwr,9664) tempo
                write(iwr,9666) tautp
           end if
           write(iwr,9668)tims
      end if
      write(iwr,9670)
      go to (1141,1142,1143,1144,1145,1146,1147,1148),ntf
 1141 write(iwr,9671)
      go to 1149
 1142 write(iwr,9672)
      go to 1149
 1143 write(iwr,9673)
      go to 1149
 1144 write(iwr,9674)
      go to 1149
 1145 write(iwr,9675)
      go to 1149
 1146 write(iwr,9676)
      go to 1149
 1147 write(iwr,9677)
      go to 1149
 1148 write(iwr,9678)
 1149 continue
      if ( ntnb .eq. 0 ) then
           write(iwr,9680)
      else
           write(iwr,9682)nsnb
      end if
      write(iwr,9684)cut,scnb,scee
      write(iwr,9686)dx0,dxm,dele,drms
      if ( ntxo .eq. 0 ) then
           write(iwr,9688)
      else if (ntxo .eq. 1) then
           write(iwr,9690)
      else
           write(iwr,9689)
      end if
      if ( ntmin .eq. 3 ) then
           if ( ntx .eq. 0 ) then
                write(iwr,9692)
           else
                write(iwr,9694)
           end if
      else
           if ( ntx .eq. 0 ) then
                write(iwr,9696)
           else
                write(iwr,9698)
           end if
      end if
      if ( osubs ) write(iwr,9701)
      if ( oaminc ) write(iwr,9704)
      cut = cut * cut
c ...
c     ----- evaluate some constants -----
c ...
      nr  = (npm * nrp) + (nsm * nram)
      nr3 = 3 * nr
_IF1(ivu)      call opfrt(ntapes(1),0,ntx,1)
c ...
c     read in coordinates form fortran channel 21.
c ...
      call getcmn(nr,rx(l30),ntx)
_IF1(v)      close(unit=ntapes(1),dispose='save')
_IF1(iu)      close(unit=ntapes(1))
c ...
c     is the input z-matrix to be changed?
c ...
      if ( osubs ) call genzmt(rx(l30))
      natg3 = 3 * nat
c ...
c    change new coordinates from a.u. into angstroms and
c    transfer to local array cg80.
c ...
_IF1(civu)      call scaler(natg3,toang(1),cg80,c)
_IFN1(civu)      call vsmul(c,1,toang(1),cg80,1,natg3)
c ...
c    align the mechanics and ab initio coordinate frames.
c    aligned coordinates are returned in rx(l30).
c ...
      call alignc(natom,mapmod,cg80,rx(l30),rx(l35),rx(l40),rx(l45))
c ...
c    ----- write the mechanics data into the dumpfile
c          store the common modfov -----
c    ----- non rotated coordinates. -----
c ...
c     ----- replace the charges if   ocharg   is on -----
c     ----- before storing mechanics data
c ...
      write(iwr,9730)
      write(iwr,9740) ( i , mapmod(i) , i = 1 , nat )
      if ( ocharg ) then
           chgold = 0.0d+00
           chgnew = 0.0d+00
           ii     = l15 - 1
           do 8010 i = 1 , nat
                  tchg(i)      = tchg(i) * chgfac
                  mapat        = iabs(mapmod(i))
                  chgold       = chgold  + rx(ii+mapat)
                  chgnew       = chgnew  + tchg(i)
                  rx(ii+mapat) = tchg(i)
8010       continue
           chgold = chgold / chgfac
           chgnew = chgnew / chgfac
           write(iwr,9750) chgold, chgnew
      end if
c
      if ( irest .le. 0 ) call amread(rx,ix,0,1,1,1,nprint)
      call resfov(1)
c ...
c      symmetry routines.
c ...
      call symm(iwr,rx)
c ...
c   the internal coordinates created from the z-matrix-specified 
c   connectivity and  the  mechanics cartesian  coordinates have 
c   now been c   generated ( genzmt ).
c   the  coordinate frame of  the "molecular  mechanics" system  
c   has now  been c   transformed  so  as  to  be  aligned
c   with  the coordinate  frame  of  the "ab initio"  system (alignc).
c   Now, it  is  necessary to  determine  the  symmetry of
c   the  "ab initio" system and  make the  appropriate  
c   transformations ( symm ).  After this  transformation  of  the
c   "ab initio"  coordinate frame it is again necessary  to realign
c   the "molecular mechanics" and "ab initio" coordinate frames. so...
c ...
      if ( nosymm .eq. 0 ) then
           write(iwr,1210)
           call rotc(rx(l30),rx(l30),nr,1)
           write(iwr,9760)
           do 8000 i = 1 , nat
              write(iwr,9765)( c(j,i) , j = 1 , 3 )
8000       continue
_IF1(civu)      call scaler(natg3,toang(1),cg80,c)
_IFN1(civu)      call vsmul(c,1,toang(1),cg80,1,natg3)
      end if
c ...
c     ----- generate the nonbonded list for the ab initio
c           atoms to be used
c           in one electron integrals -----
c ...
      call genbel(natom,natbel,mapmod,ix(i60),0)
      call stuf80(natom,mapmod,ix(i60),rx(l30),cg80,0)
_IF1(ivu)      call opfrt(ntapes(2),-1,ntxo,0)
      call minrit(nres,rx(l30),rx(l10),rx(l05),ix(i02))
_IF1(v)      close(unit=ntapes(2),dispose='save')
_IF1(iu)      close(unit=ntapes(2))
c ...
c     ----- if ntmin.gt.0 branch for minimisation -----
c ...
      if ( ntmin .le. 0) go to 100
      if ( ntmin .eq. 1) istep = 5
      if ( ntmin .eq. 2) istep = 6
      if ( ntmin .eq. 3) istep = 8
      go to (1400,1400,1600) , ntmin
c ...
c     ----- write the mechanics data into the dumpfile and the
c           mechanics coordinates into the restart file
c           if no minimisation -----
c ...
  100 natmod = natom
      call nbg80(natom,nres,ix(i08),ix(i10),ix(i02),rx(l35),rx(l30),cut)
c
c     ----- write the nbact array which is in rx(l35)
c           onto the dumpfile -----
c     ----- write mechanics coordinates ... rx(l30) to dumpfile. -----
c ...
      call secput(isect(472),m51,lensec(nr3),ibl172)
      skal = 1.0d0 / toang(1)
      call dscal(nr3,skal,rx(l30),1)
      call wrt3(rx(l30),nr3,ibl172,idaf)
      len173 = (natom-1) / nav + 1
      call secput(isect(473),m51,lensec(len173),ibl173)
      call wrt3(rx(l35),len173,ibl173,idaf)
c
c     ----- write the charges for 1e integrals -----
c
      skal = 1.0d0 / chgfac
      call dscal(natom,skal,rx(l15),1)
      call secput(isect(474),m51,lensec(natom),ibl174)
      call wrt3(rx(l15),natom,ibl174,idaf)
      if ( nprint .ne. -5 ) then
           write(iwr,9722)ibl172,lensec(nr3)
           write(iwr,9723)ibl173,lensec(len173)
           write(iwr,9724)ibl174,lensec(natom)
      end if
      return
c ...
c     ----- istep = 3 and 4   energy of the qm atoms
c           due to surroundings -----
c ...
c     ----- define the ab initio atoms as belly
c           without the junction atoms --
c ...
 1300 if ( istep .eq. 3 ) write(iwr,9518)
      if ( istep .eq. 4 ) write(iwr,9528)
      call amread(rx,ix,1,1,1,1,nprint)
      call resfov(0)
      call genbel(natom,natbel,mapmod,ix(i60),0)
      npro = natbel
c ...
c     ----- zero the ab initio charges if the surrounding is already
c           included into the 1e hamiltonian -----
c ...
      if ( .not. oaminc ) go to 1340
      ii = l15 - 1
      do 1320 i = 1 , nat
             iat = mapmod(i)
             if ( iat .lt. 0 ) go to 1320
             rx(ii+iat) = 0.0d+00
 1320 continue
c ...
c     ----- inject the ab initio atoms into mechanics -----
c ...
 1340 if ( nosymm .eq. 0 ) then
           write(iwr,1210)
           call rotc(rx(l30),rx(l30),nr,1)
      end if
      natg3 = 3 * nat
_IF1(civu)      call scaler(natg3,toang(1),cg80,c)
_IFN1(civu)      call vsmul(c,1,toang(1),cg80,1,natg3)
      call stuf80(natom,mapmod,ix(i60),rx(l30),cg80,0)
c ...
c     ----- keep only the junction bonds etc -----
c ...
      ogauss = .true.
      call varadj(ix,ogauss)
c ...
c     ----- now evaluate the energy and gradients -----
c ...
      ntnb = 1
      ndrv = 1
      call bforce(rx,ix,ene,ndrv,nbnum,ogauss)
c ...
c     ----- if only gradient evaluation skip energy -----
c ...
      if ( istep .eq. 4 ) go to 1380
c ...
c           convert form kcals. to atomic units.
c ...
      emod = ene(1) / atocal
      enuc = 0.0d+00
c ...
c     ----- if the surrounding is included in the 1e hamiltonian
c           evaluate the nuclear repulsion energy. -----
c ...
      if ( .not. oaminc ) go to 1360
c
c       call scaler(natg3,toang(1),cg80,c)
c
      call secget(isect(473),m51,ibl173)
      call rdedx(rx(l40),len173,ibl173,idaf)
_IF1()      call zvnuc(natom,mapmod,enuc,cg80,rx(l30),imass,rx(l15),rx(l40))
      call zvnuc(natom,enuc,cg80,rx(l30),imass,rx(l15),rx(l40))
      enuc = enuc * toang(1) / chgfac
 1360 continue
      enrgy = enrgy + emod + enuc
      write(iwr,9108) emod,enuc
      go to 1399
c      return
c ...
c     ----- if istep.eq.4 transfer the gradients -----
c
c     ----- store the coordinates at each gradient evaluation -----
c ...
 1380  continue
_IF1(ivu)      call opfrt(ntapes(2),-1,ntxo,0)
      call minrit(nres,rx(l30),rx(l10),rx(l05),ix(i02))
_IF1(v)      close(unit=ntapes(2),dispose='save')
_IF1(iu)      close(unit=ntapes(2))
c ...
c     ----- transfer the gradient    rx(l35) --> cg80    -----
c ...
      call gath80(natom,mapmod,ix(i60),rx(l35),cg80,0)
      dff  = dsqrt( ddot(natg3,cg80,1,cg80,1) / dfloat(natg3))
      gmax = grdmax(natg3,cg80)
      call prmode(nrep,1,dff,gmax,ene)
c ...
c     convert to atomic units ( length ) by mulitplying by toang(1)
c     convert to atomic units ( energy ) by dividing by atocal
c ...
      fact = toang(1) / atocal
c ...
c     ----- update the gradient -----
c ...
      call dscal(natg3,fact,cg80,1)
c
      do 8020 i = 1 , nat
             i3 = 3 * i - 3
             if ( mapmod(i) .le. 0 ) go to 240
             j31 = i3 + 1
             j32 = i3 + 2
             j33 = i3 + 3
             de(j31) = de(j31) -  cg80(j31) 
             de(j32) = de(j32) -  cg80(j32)
             de(j33) = de(j33) -  cg80(j33) 
 240         continue
8020   continue
c
      write(iwr,1391)
      mmax = 0
      mst = 0
      nst = 0
 1390 mmin = mmax + 1
      mmax = mmax + 7
      if (mmax .gt. nat )  mmax = nat
      write(iwr,1322)
      write(iwr,1324) ( i , i = mmin , mmax )
      write(iwr,1322)
      mst = nst +1
      nst = mmax * 3
      do 8030 i = 1 , 3
         j = i - 1
         write(iwr,1325) zdnam(i),( cg80(loop+j),loop=mst,nst,3 )
8030  continue
      if (mmax .lt. nat ) go to 1390
c ...
c
c     ----- write mechanics coordinates ... rx(l30) to dumpfile
c ...
 1399 call secput(isect(472),m51,lensec(nr3),ibl172)
      skal = 1.0d0 / toang(1)
      call dscal(nr3,skal,rx(l30),1)
      call wrt3(rx(l30),nr3,ibl172,idaf)
      return
c ...
c     ----- istep = 5 and 6   surrounding or ab initio atoms are
c           minimised -----
c ...
1400   continue
_IF1(ivu)      call opfrt(ntapes(2),-1,ntxo,0)
      call mdbox(ntb,bbeta,box)
      if ( istep .eq. 5 ) igath = -1
      if ( istep .eq. 6 ) igath =  1
      if ( istep .eq. 5 ) write(iwr,9538)
      if ( istep .eq. 6 ) write(iwr,9548)
c ...
c     ----- gather the atoms depending on istep -----
c ...
      call amread(rx,ix,1,1,1,1,nprint)
      call resfov(0)
      call genbel(natom,natbel,mapmod,ix(i60),igath)
      npro = natbel
      ogauss = .false.
      call varadj(ix,ogauss)
      write(iwr,9218) npro
c
c     ----- gather the belly atoms in rx(l45) for the minimiser -----
c
      call gath80(natom,mapmod,ix(i60),rx(l30),rx(l45),igath)
      call bfast(rx,ix,
     +           rx(l30),rx(l45),rx(l35),rx(l40),rx(l65),rx(l10),rx(l05)
     +          ,ix(i02),ix(i60),oerstp,oconv,mapmod,igath)
      call minrit(nres,rx(l30),rx(l10),rx(l05),ix(i02))
      call amread(rx,ix,0,1,0,0,nprint)
      if(oconv) write(iwr,9128)
      if(oerstp) write(iwr,9138)
      return
c ...
c     ----- calculation of the hessian from model
c ...
 1500 write(iwr,9558)
      call amread(rx,ix,1,1,1,1,nprint)
      call resfov(0)
c ...
c     ----- calculate the cartesian hessian including the junction
c           atoms -----
c ...
      call cal2d(rx,ix,mapmod,nz)
      nvarx = 3 * nat
      nz3   = 3 * nz
      okcal  = .true.
c ...
c     ----- the first and second cartesian derivatives are now stored
c           in rx(l05) and rx(nz3+l05) respectively. if istep .eq. 8
c           output it and skip the internal coordinate transformation
c ...
      if ( istep .eq. 7 ) go to 520
c      call d2xout(nvarx,rx(1),rx(nz3+1),okcal)
      call d2xout(nvarx,rx(l05),rx(nz3+l05),okcal)
      go to 540
  520 continue
      call ztranf(rx,nvarx,okcal,kprint)
  540 continue
      return
c ...
c     ----- molecular dynamics calculation on the surroundings -----
c ..
 1600 write(iwr,9568)
      call amread(rx,ix,1,1,1,1,nprint)
      call resfov(0)
      call genbel(natom,natbel,mapmod,ix(i60),-1)
      npro = natbel
      ogauss = .false.
      call varadj(ix,ogauss)
      write(iwr,9218) npro
c ...
c     ----- enter into the dynamics -----
c ...
      call runmd(rx,ix,natom,npro,nres,rx(l10),
     * ix(i02),rx(l05),rx(l30),rx(l35)
     &           ,rx(l40),rx(l20),ix(i60),ogauss)
c ...
c     ----- all done return -----
c ...
      return
c ...
c     ------- tidy up memory -----
c ...
 1700 call lcmem2
      return
 1210 format(/4x,'transformation of mechanics coordinates after ',
     + 'orientation of ab initio system coordinates',/)
 1391 format(/4x,19('='),/4x,'mechanics gradients',/4x,19('='))
 1322 format(/)
 1324 format(5x,'atom',5x,7(i2,15x))
 1325 format(7x,a5,2x,7(e15.8,2x))
 9108 format(/4x,'energy due to the surroundings =',f16.8,
     +       /4x,'nuclear repulsion energy       =',f16.8)
 9138 format(/4x,' minimisation aborted ... error return')
 9128 format(/4x,'successful completion of molecular mechanics run.')
 9218 format(/4x,'number of surrounding atoms =',i5)
_IF1() 9508 format(/4x,28('*'),/4x,'model generation of z matrix',/4x,
_IF1()     +            28('*'),//)
 9518 format(/4x,28('*'),/4x,'mechanics energy calculation',
     +       /4x,28('*'),//)
 9528 format(/4x,18('*'),/4x,'mechanics gradient',/4x,18('*'),//)
 9538 format(/4x,37('*'),/4x,
     *     'mechanics minimisation of surrounding',
     +       /4x,37('*'),//)
 9548 format(/4x,40('*'),/4x,
     +     'mechanics minimisation of ab intio atoms',
     +       /4x,40('*'),//)
 9558 format(/4x,29('*'),/4x,'mechanics hessian calculation',
     +       /4x,29('*'),//)
 9568 format(/4x,32('*'),/4x,'dynamical shaking of surrounding',
     +       /4x,32('*'),//)
 9600 format(/4x,
     +'topology file created in model job with the title :- ',//4x,
     +       20a4)
 9620 format(/4x,'a distance-dependent dielectric funtion is to used',
     +' in calculating the electrostatic energy.')
 9622 format(/4x,'a constant dielectric funtion is to used in',
     +' calculating the electrostatic energy :- ',f10.5)
 9630 format(/4x,
     +'maximum number of cycles of mechanics optimisation = ',
     + i5)
 9632 format(/4x,
     +'maximum number of cycles of quantum optimisation = ',
     + i5)
 9640 format(/4x,'control data for molecular dynamics run ')
 9642 format(/4x,
     +'initial velocities calculated assuming a maxwellian',
     +           ' distribution at ',f10.5, ' kelvin.')
 9644 format(/4x,'initial velocities are to be read through fortran',
     +           ' channel 21 [formatted file]')
 9646 format(/4x,'initial velocities are to be read through fortran',
     +           ' channel 21 [unformatted file]')
 9650 format(/4x,'periodicity is to be applied in the dynamics run.',/
     +           ' the box is to be an octahedral prism.')
 9652 format(/4x,
     +' no periodicity is to be applied during the dynamics run.')
 9654 format(/4x,'periodicity is to be applied in the dynamics run.',/
     +           ' the box is to be a rectangluar prism.')
 9656 format(/4x,'periodicity is to be applied in the dynamics run.',/
     +           ' the box is to be a monclinic prism.')
 9660 format(/4x,'the lengths of the sides of the box are to be :-',/
     +           f10.5,5x,f10.5,5x,f10.5)
 9662 format(/4x,
     +'constant pressure is to be applied during the dynamics run.')
 9664 format(/4x,
     +'a constant temperature of ',f6.2,' kelvin is to be maintained
     + during the dynamics run.')
 9666 format(/4x,'coupling coefficient for scaling velocities to',
     +           ' maintain constant temperature = ',f10.5)
 9668 format(/4x,'time step for dynamics run = ',f10.5,' picoseconds.')
 9670 format(/4x,'method of force evaluation :-')
 9671 format('+',36x,'complete interaction is calculated.')
 9672 format('+',36x,
     + 'bond interactions involving h-atoms are to be omitted.')
 9673 format('+',36x,'all the bond interactions are to be omitted.')
 9674 format('+',36x,
     +'angles involving h-atoms and all bond interactions are to be',
     +' omitted')
 9675 format('+',36x,
     +'all bond and angle interactions are to be omitted.')
 9676 format('+',36x,'dihedrals involving h-atoms and all bonds and ',
     +           'all angles are to be omitted.')
 9677 format('+',36x,'all bond, angle and dihedral interactions are to',
     +           ' be omitted.')
 9678 format('+',36x,'all bond, angle, dihedral and non-bonded ',
     +           'interactions are to be omitted.')
 9680 format(/4x,'no pair list will be generated.')
 9682 format(/4x,'a pair list will be generated every ',i4,' cycles.')
 9684 format(/4x,
     +'cutoff distance for non-bonded pairs            = ',f10.5,
     +' angs.',/,4x,
     +'scale factor for 1-4 vdw interactions           = ',f10.5,/,4x,
     +'scale factor for 1-4 electrostatic interactions = ',f10.5)
 9686 format(/4x,
     +'initial step length in mechanics optimsation ',
     +'                        = ',f10.5,/,4x,
     +'maximum step length allowed in mechanics optimisation ',
     +'               = ',f10.5,/,4x,
     +'convergence criterium for energy in mechanics optimisation ',
     +'          = ',f10.5,/,4x,
     +'convergence criterium for norm of gradient in mechanics ',
     +'optimisation = ',f10.5,/)
 9688 format(/4x,'final / restart coordinates are to be written in a ',
     +'formatted form.')
 9689 format(/4x,'final / restart coordinates are to be written in a ',
     +'pdb formatted form.')
 9690 format(/4x,'final / restart coordinates are to be written in an ',
     +'unformatted form.')
 9692 format(/4x,
     +'initial velocities will be read form a formatted file.')
 9694 format(/4x,
     +'initial velocities will be read form an unformatted file.')
 9696 format(/4x,
     +'initial coordinates will be read form a formatted file.')
 9698 format(/4x,
     +'initial coordinates will be read form an unformatted file.')
 9701 format(/4x,'z-matrix variables will be regenerated from',
     +' atomic coordinates.',/)
 9699 format(/5x,25('*'),/5x,'initiate mechanics module',
     +       /5x,25('*'),/)
 9700 format(/4x,'reading formatted ',
     +'m o l e c u l a r   t o p o l o g y   f i l e',/)
 9702 format(/4x,'reading unformatted  ( binary )',
     +'m o l e c u l a r   t o p o l o g y   f i l e',/)
 9704 format(/4x,'partial charges will be included in one-electron',
     +' hamiltonian.')
 9722 format(/1x,
     +'writing dumpfile section 472 to block ',i5/1x,
     +'length of section 472               = ',i5/)
 9723 format(/1x,
     +'writing dumpfile section 473 to block ',i5/1x,
     +'length of section 473               = ',i5/)
 9724 format(/1x,
     +'writing dumpfile section 474 to block ',i5/1x,
     +'length of section 474               = ',i5/)
 9730 format(/5x,'ab initio -> mechanics mapping of atoms',/
     +5x,39('='),/)
 9740 format(/4x, 8(1x,'(',i5, '->', i5,')',1x,',') )
 9750 format(/4x,'the charges of the ab initio atoms have',
     +'been replaced.',/
     +    ,/4x,'net charge (old) =',f10.5,' net charge (new) =',f10.5,/)
 9760 format(/4x,'coordinates of ab initio system after',
     +' orientation.',//,11x,'x',15x,'y',15x,'z',/)
 9765 format(' ',e15.8,5x,e15.8,5x,e15.8)
      end
      subroutine rdparm(rx,ix)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ix
INCLUDE(common/sizes)
c
INCLUDE(common/modfov)
      common/memla/i02,i04,i06,i08,i10,i12,i14,i16,i18,i20,
     +             i22,i24,i26,i28,i30,i32,i34,i36,i38,i40,
     +             i42,i44,i46,i48,i50,i52,i54,i56,i58,i60,
     +             i62,i64,i66,i68,i70,i72,i74,i76,i78,i80
      common/memlb/l05,l10,l15,l20,l25,l30,l35,l40,l45,l50,
     +             l55,l60,l65,l70,l75,l80,l85,l90,l95,l100,
     +             k02,k04,k06,k08,k10,k12,k14,k16,k18,k20,
     +             k22,k24,k26,k28,k30,k32,k34,k36,k38,k40
INCLUDE(common/nbterm)
INCLUDE(common/runhed)
INCLUDE(common/iofile)
INCLUDE(common/modj)
c
      dimension ix(*),rx(*)
c
c     ----- read the molecular topology -----
c
      if(kform.gt.0) go to 300
c
c     ----- binary format -----
c
_IF1(ivu)      call opfrt(ntapes(3),0,1,1)
      read(ntapes(3)) (ititl(i),i=1,20)
      read(ntapes(3)) natom,ntypes,nbonh,mbona,ntheth,mtheta,
     +                nphih,mphia,
     +                nhparm,nparm,nnb,nres,nbona,ntheta,nphia,numbnd,
     +                numang,numphi,natyp,nphb,idum,idum,idum,idum,
     +                idum,idum,idum,idum
c
c     ----- now all the data is available for memory allocation.
c           do the partition and load the arrays. both belly
c           and constrained minimisation share the same arrays -----
c
      if(nrcp.gt.0 .and. npro.gt.0) go to 900
      call lcmem0
      call lcmem1
c
c     nat3 = 3*natom
      ntype = ntypes*ntypes
      nttyp = ntypes*(ntypes+1)/2
c
c     ----- read the symbols and the charges and the masses -----
c
      call getsng(ntapes(3),natom,rx(l10),kform)
      call getdbl(ntapes(3),natom,rx(l15),rx(l15))
      call getdbl(ntapes(3),natom,rx(l20),rx(l20))
c
      read(ntapes(3)) (ix(i04+i-1),i = 1,natom)
      read(ntapes(3)) (ix(i+i08-1),i = 1,natom)
      read(ntapes(3)) (ix(i+i06-1),i = 1,ntype)
      call getsng(ntapes(3),nres,rx(l05),kform)
      read(ntapes(3)) (ix(i+i02-1),i=1,nres)
      ix(i02+nres) = natom+1
c
c     ----- read the parameters -----
c
      call getdbl(ntapes(3),numbnd,rx(k02),rx(k02))
      call getdbl(ntapes(3),numbnd,rx(k04),rx(k04))
      call getdbl(ntapes(3),numang,rx(k06),rx(k06))
      call getdbl(ntapes(3),numang,rx(k08),rx(k08))
      call getdbl(ntapes(3),numphi,rx(k10),rx(k10))
      call getdbl(ntapes(3),numphi,rx(k12),rx(k12))
      call getdbl(ntapes(3),numphi,rx(k14),rx(k14))
      call getdbl(ntapes(3),natyp,rx(k16),rx(k16))
      call getdbl(ntapes(3),nttyp,rx(k24),rx(k24))
      call getdbl(ntapes(3),nttyp,rx(k26),rx(k26))
c
c     ----- read the information regarding bonding -----
c
      read(ntapes(3)) (ix(i+i12-1),ix(i+i14-1),ix(i+i16-1),i = 1,nbonh)
      read(ntapes(3)) (ix(i+i18-1),ix(i+i20-1),ix(i+i22-1),i = 1,nbona)
      read(ntapes(3)) (ix(i+i24-1),ix(i+i26-1),ix(i+i28-1),ix(i+i30-1),
     +          i = 1,ntheth)
      read(ntapes(3)) (ix(i+i32-1),ix(i+i34-1),ix(i+i36-1),ix(i+i38-1),
     +          i = 1,ntheta)
      read(ntapes(3)) (ix(i+i40-1),ix(i+i42-1),ix(i+i44-1),ix(i+i46-1),
     +          ix(i+i48-1),i = 1,nphih)
      read(ntapes(3)) (ix(i+i50-1),ix(i+i52-1),ix(i+i54-1),ix(i+i56-1),
     +          ix(i+i58-1),i = 1,nphia)
      read(ntapes(3)) (ix(i+i10-1),i=1,nnb)
c
c     ----- read the h-bond parameters -----
c
      call getdbl(ntapes(3),nphb,rx(k18),rx(k18))
      call getdbl(ntapes(3),nphb,rx(k20),rx(k20))
      call getdbl(ntapes(3),nphb,rx(k22),rx(k22))
c
c     ----- read isymbl,itree,join and irotat arrays -----
c
      call getsng(ntapes(3),natom,rx(l55),kform)
      call getsng(ntapes(3),natom,rx(l60),kform)
      read(ntapes(3)) (ix(i+i62-1),i=1,natom)
      read(ntapes(3)) (ix(i+i64-1),i=1,natom)
      go to 320
c
c     ----- formatted input -----
c
  300 continue
_IF1(ivu)      call opfrt(ntapes(3),0,0,1)
      read(ntapes(3),9008) (ititl(i),i=1,20)
      read(ntapes(3),9018) natom,ntypes,nbonh,mbona,
     +              ntheth,mtheta,nphih,mphia,
     +              nhparm,nparm,nnb,nres,nbona,ntheta,nphia,numbnd,
     +              numang,numphi,natyp,nphb,idum,idum,idum,idum,
     +              idum,idum,idum,idum
c
c     ----- now all the data is available for memory allocation.
c           do the partition and load the arrays. both belly
c           and constrained minimisation share the same arrays -----
c
      if(nrcp.gt.0 .and. npro.gt.0) go to 900
      call lcmem0
      call lcmem1
c
c     nat3 = 3*natom
      ntype = ntypes*ntypes
      nttyp = ntypes*(ntypes+1)/2
c
c     ----- read the symbols and the charges and the masses -----
c
      call getsng(ntapes(3),natom,rx(l10),kform)
      read(ntapes(3),9028) (rx(l15+i-1),i=1,natom)
      read(ntapes(3),9028) (rx(l20+i-1),i=1,natom)
c
      read(ntapes(3),9018) (ix(i04+i-1),i = 1,natom)
      read(ntapes(3),9018) (ix(i+i08-1),i = 1,natom)
      read(ntapes(3),9018) (ix(i+i06-1),i = 1,ntype)
      call getsng(ntapes(3),nres,rx(l05),kform)
      read(ntapes(3),9018) (ix(i+i02-1),i=1,nres)
      ix(i02+nres) = natom+1
c
c     ----- read the parameters -----
c
      read(ntapes(3),9028) (rx(k02+i-1),i=1,numbnd)
      read(ntapes(3),9028) (rx(k04+i-1),i=1,numbnd)
      read(ntapes(3),9028) (rx(k06+i-1),i=1,numang)
      read(ntapes(3),9028) (rx(k08+i-1),i=1,numang)
      read(ntapes(3),9028) (rx(k10+i-1),i=1,numphi)
      read(ntapes(3),9028) (rx(k12+i-1),i=1,numphi)
      read(ntapes(3),9028) (rx(k14+i-1),i=1,numphi)
      read(ntapes(3),9028) (rx(k16+i-1),i=1,natyp)
      read(ntapes(3),9028) (rx(k24+i-1),i=1,nttyp)
      read(ntapes(3),9028) (rx(k26+i-1),i=1,nttyp)
c
c     ----- read the information regarding bonding -----
c
      read(ntapes(3),9018)
     +     (ix(i+i12-1),ix(i+i14-1),ix(i+i16-1),i = 1,nbonh)
      read(ntapes(3),9018)
     +     (ix(i+i18-1),ix(i+i20-1),ix(i+i22-1),i = 1,nbona)
      read(ntapes(3),9018)
     +     (ix(i+i24-1),ix(i+i26-1),ix(i+i28-1),ix(i+i30-1),
     +               i = 1,ntheth)
      read(ntapes(3),9018)
     +     (ix(i+i32-1),ix(i+i34-1),ix(i+i36-1),ix(i+i38-1),
     +               i = 1,ntheta)
      read(ntapes(3),9018)
     +     (ix(i+i40-1),ix(i+i42-1),ix(i+i44-1),ix(i+i46-1),
     +               ix(i+i48-1),i = 1,nphih)
      read(ntapes(3),9018)
     +     (ix(i+i50-1),ix(i+i52-1),ix(i+i54-1),ix(i+i56-1),
     +               ix(i+i58-1),i = 1,nphia)
      read(ntapes(3),9018)(ix(i+i10-1),i=1,nnb)
c ...
c     ----- read the h-bond parameters -----
c ...
      read(ntapes(3),9028) (rx(k18+i-1),i=1,nphb)
      read(ntapes(3),9028) (rx(k20+i-1),i=1,nphb)
      read(ntapes(3),9028) (rx(k22+i-1),i=1,nphb)
c ...
c     ----- read isymbl,itree,join and irotat arrays -----
c ...
      call getsng(ntapes(3),natom,rx(l55),kform)
      call getsng(ntapes(3),natom,rx(l60),kform)
      read(ntapes(3),9018) (ix(i+i62-1),i=1,natom)
      read(ntapes(3),9018) (ix(i+i64-1),i=1,natom)
c ...
c     ----- calculate the inverse of the masses -----
c ...
  320 continue
      ii = l20 - 1
      do 120 i = 1 , natom
             rxiii = rx(i+ii)
             rxiii = 1.0d0 / rxiii
             rx(i+ii) = rxiii
  120 continue
c ...
c     ----- scale the charges if dielc .gt. 1.0e0 -----
c ...
      if ( dielc .le. 1.0d0 ) go to 160
      rdumd = 1.0d0 / dsqrt(dielc)
      call dscal(natom,rdumd,rx(l15),1)
  160 continue
_IF1(v)      close(unit=ntapes(3),dispose='save')
_IF1(iu)      close(unit=ntapes(3))
      return
  900 continue
      write(iwr,9108)
 9008 format(20a4)
 9018 format(12i6)
 9028 format(5e16.8)
 9108 format(/5x,'warning ... both belly and constraint flags are '
     +         , 'on ... not allowed')
      call caserr('invalid flags in model')
      return
      end
      subroutine rotc(a,b,n,iway)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
INCLUDE(common/sizes)
INCLUDE(common/phycon)
INCLUDE(common/molsym)
      dimension a(*),b(*)
      trxb = trx * toang(1)
      tryb = try * toang(1)
      trzb = trz * toang(1)
      i3  = 0
      if ( iway .ge. 0 ) then
c ...
c forward transformation.
c ...
           do 8000 i = 1 , n
                  a1      = a(i3+1)
                  a2      = a(i3+2)
                  a3      = a(i3+3)
                  a1      = a1 + trxb
                  a2      = a2 + tryb
                  a3      = a3 + trzb
                  b(i3+1) = tr(1,1)*a1 + tr(1,2)*a2 + tr(1,3)*a3
                  b(i3+2) = tr(2,1)*a1 + tr(2,2)*a2 + tr(2,3)*a3
                  b(i3+3) = tr(3,1)*a1 + tr(3,2)*a2 + tr(3,3)*a3
                  i3      = i3 + 3
8000  continue
      else
c ...
c backward transformation.
c ...
           do 8010 i = 1 , n
                  a1      = a(i3+1)
                  a2      = a(i3+2)
                  a3      = a(i3+3)
                  b(i3+1) = tr(1,1)*a1 + tr(2,1)*a2 + tr(3,1)*a3
                  b(i3+2) = tr(1,2)*a1 + tr(2,2)*a2 + tr(3,2)*a3
                  b(i3+3) = tr(1,3)*a1 + tr(2,3)*a2 + tr(3,3)*a3
                  a1      = a1 - trxb
                  a2      = a2 - tryb
                  a3      = a3 - trzb
                  i3      = i3 + 3
8010       continue
      end if
      return
      end
c
      subroutine resfov(index)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
c
c     ----- routine to restore the common /modfov/ for each
c           entry into model -----
c
INCLUDE(common/modfov)
      dimension ihc(2)
      equivalence (ihc(1),natom)
      data len/50/
      lenn=len/lenwrd()
      if ( index .le. 0 ) then
c
c     ----- restore the common -----
c
_IF1(cuf)      call fmove(jmodf,ihc,lenn)
_IFN1(cuf)      call icopy(lenn,jmodf,1,ihc,1)
      else
c
c     ----- store the common -----
c
_IF1(cuf)      call fmove(ihc,jmodf,lenn)
_IFN1(cuf)      call icopy(lenn,ihc,1,jmodf,1)
      endif
c
      return
      end
      subroutine amread(rx,ix,index,ireca,irecb,irecc,nprint)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
_IFN1(cuf)      integer*2 ix
c ...
c     ----- routine to load the mechanics arrays into memory -----
c ...
INCLUDE(common/iofile)
      common/memla/i02,i04,i06,i08,i10,i12,i14,i16,i18,i20,
     +             i22,i24,i26,i28,i30,i32,i34,i36,i38,i40,
     +             i42,i44,i46,i48,i50,i52,i54,i56,i58,i60,
     +             i62,i64,i66,i68,i70,i72,i74,i76,i78,i80
      common/memlb/l05,l10,l15,l20,l25,l30,l35,l40,l45,l50,
     +             l55,l60,l65,l70,l75,l80,l85,l90,l95,l100,
     +             k02,k04,k06,k08,k10,k12,k14,k16,k18,k20,
     +             k22,k24,k26,k28,k30,k32,k34,k36,k38,k40
INCLUDE(common/mod2d)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      dimension rx(*),ix(*)
      data  m51 /51/
c
      len175   = l35 - l05
      len176   = l65 - l50
      ln177a   = k28 - k02
      ln177b   = (i80-i02+nword-1)/nword
      oreca    = ireca .gt. 0
      orecb    = irecb .gt. 0
      orecc    = irecc .gt. 0
      nav = lenwrd()
c ...
c     ----- load or unload depending on index -----
c ...
      if ( index .gt. 0 ) then
c ...
c     ----- load the arrays from the dumpfile.
c ...
       if ( oreca ) then
           call secget(isect(475),m51,ibl175)
           call rdedx(rx(l05),len175,ibl175,idaf)
       end if
       if ( orecb ) then
           call secget(isect(476),m51,ibl176)
           call rdedx(rx(l50),len176,ibl176,idaf)
       end if
       if ( orecc ) then
           call secget(isect(477),m51,ibl177)
           call rdedx(rx(k02),ln177a,ibl177,idaf)
           call readis(ix(i02),ln177b*nav,idaf)
       end if
c ...
c     ----- unload the arrays into the direct access files -----
c ...
      else
       if ( oreca ) then
           call secput(isect(475),m51,lensec(len175),ibl175)
           if ( nprint .ne. -5 ) write(iwr,1175) ibl175,lensec(len175)
           call wrt3 (rx(l05),len175,ibl175,idaf)
       end if
       if ( orecb ) then
       call secput(isect(476),m51,lensec(len176),ibl176)
       if ( nprint .ne. -5 ) write(iwr,1176) ibl176,lensec(len176)
       call wrt3 (rx(l50),len176,ibl176,idaf)
       end if
       if ( orecc ) then
       len177 = lensec(ln177a) + lensec(ln177b)
       call secput(isect(477),m51,len177,ibl177)
       if( nprint .ne. -5 ) write(iwr,1177) ibl177,len177
            call wrt3 (rx(k02),ln177a,ibl177,idaf)
            call wrt3s(ix(i02),ln177b,idaf)
       end if
      endif
      return
 1175 format(/1x,
     +'writing dumpfile section 475 to block ',i5/1x,
     +'length of section 475               = ',i5/)
 1176 format(/1x,
     +'writing dumpfile section 476 to block ',i5/1x,
     +'length of section 476               = ',i5/)
 1177 format(/1x,
     +'writing dumpfile section 477 to block ',i5/1x,
     +'length of section 477               = ',i5/)
      end
_IF(ibm,vax,cyber205)
      subroutine opfrt(nf,istat,iform,ir)
      character*12 form,stat
c
c     ----- routine to open different types of files for sequential
c           accessing. the file name is standard default or properly
c           assigned ones -----
c
c           nf           io unit number
c
c           istat        status key
c                   = 0  status is 'old'
c                   = 1  status is 'new'
c                   =-1  status is 'unknown'
c
c           iform        formatted or unformatted access
c                   = 0  form is 'formatted'
c                   = 1  form is 'unformatted'
c                   = 2  form is 'formatted pdb'
c
c           ir           read and write access
c                   = 0  both read and write access
c                   = 1  readonly access
c
      dimension form(2),stat(3)
      data form/'formatted   ','unformatted '/
      data stat/'old         ','new         ','unknown     '/
c
c     ----- setup the status expression -----
c
      if(istat) 140,100,120
  120 is = 2
      go to 160
  100 is = 1
      go to 160
  140 is = 3
  160 continue
c
c     ----- setup the format expression -----
c
      if(iform.eq.1) then
        ifo = 2
      else
        ifo = 1
      endif
c
c     ----- now open the file -----
c
      if(ir.gt.0) go to 180
c
_IF1(v)      open(unit = nf,status = stat(is),form = form(ifo),shared)
_IFN1(v)      open(unit = nf,status = stat(is),form = form(ifo))
      go to 200
  180 continue
_IF1(xv)      open(unit = nf,status = stat(is),form = form(ifo),readonly)
_IFN1(xv)      open(unit = nf,status = stat(is),form = form(ifo))
  200 continue
      return
      end
_ENDIF
      subroutine getdbl(nf,n,rx,rxd)
c rx changed to *8, to replace implicit (ps)
      REAL  rxd, rx
      dimension rx(*),rxd(*)
      read(nf) (rx(i),i=1,n)
      ni = n + 1
      do 100 i = 1 , n
             idum = ni - i
             rxd(idum) = rx(idum)
  100 continue
      return
      end
      subroutine getsng(nf,n,ix,kform)
      dimension ix(*)
      if(kform.le.0) read(nf) (ix(i),i=1,n)
      if(kform.gt.0) read(nf,9008) (ix(i),i=1,n)
      return
 9008 format(20a4)
      end
      subroutine lcmem0
      implicit REAL  (a-h,o-z)
c
c     ----- routine to determine the memory for various
c           arrays according to their needs -----
c           note that lcmem0 and lcmem1 perform the
c           same role as the original locmem
c
INCLUDE(common/gmempara)
INCLUDE(common/modmin)
INCLUDE(common/modfov)
      common/memla/i02,i04,i06,i08,i10,i12,i14,i16,i18,i20,
     +             i22,i24,i26,i28,i30,i32,i34,i36,i38,i40,
     +             i42,i44,i46,i48,i50,i52,i54,i56,i58,i60,
     +             i62,i64,i66,i68,i70,i72,i74,i76,i78,i80
      common/memlb/l05,l10,l15,l20,l25,l30,l35,l40,l45,l50,
     +             l55,l60,l65,l70,l75,l80,l85,l90,l95,l100,
     +             k02,k04,k06,k08,k10,k12,k14,k16,k18,k20,
     +             k22,k24,k26,k28,k30,k32,k34,k36,k38,k40
INCLUDE(common/mod2d)
INCLUDE(common/iofile)
c
      dimension la(2),lb(2)
      equivalence (la(1),i02),(lb(1),l05)
      character *7 fnm
      character *6 snm
      data fnm/"model.m"/
      data snm/"lcmem0"/
c
c     ----- length of arrays in dynamic common -----
c           these are refered by the eventual pointers l* set in lcmem1-----
c
c           lbres   ...  l05
c           igraph  ...  l10
c           cg      ...  l15
c           amass   ...  l20
c           rxchrg   ...  l25
c           c       ...  l30
c           f       ...  l35
c           v       ...  l40
c           rxr      ...  l45
c           wt      ...  l50
c           isymbl  ...  l55
c           itree   ...  l60
c           work    ...  l65
c
c     ----- arrays refered by the pointers k* -----
c
c           rk      ...  k02
c           req     ...  k04
c           tk      ...  k06
c           teq     ...  k08
c           pk      ...  k10
c           pn      ...  k12
c           phase   ...  k14
c           solty   ...  k16
c           asol    ...  k18
c           bsol    ...  k20
c           hbcut   ...  k22
c           cn1     ...  k24
c           cn2     ...  k26
c
c     ----- arrays refered by pointers i* -----
c           for belly minimisation igres == irotat
c
c           ipres   ...  i02      ita     ...  i32
c           iac     ...  i04      jta     ...  i34
c           no      ...  i06      kta     ...  i36
c           iblo    ...  i08      icta    ...  i38
c           inb     ...  i10      iph     ...  i40
c           ibh     ...  i12      jph     ...  i42
c           jbh     ...  i14      kph     ...  i44
c           icbh    ...  i16      lph     ...  i46
c           iba     ...  i18      icph    ...  i48
c           jba     ...  i20      ipa     ...  i50
c           icba    ...  i22      jpa     ...  i52
c           ith     ...  i24      kpa     ...  i54
c           jth     ...  i26      lpa     ...  i56
c           kth     ...  i28      icpa    ...  i58
c           icth    ...  i30      igroup  ...  i60
c           join    ...  i62
c           irotat  ...  i64
c           iar1    ...  i78
c           iar2    ...  i80
c
c     ----- now make the partition -----
c
      l0 = nres+5
      l1 = natom+5
      l2 = 5
      if(npro.gt.0 .or. nrcp.gt.0) l2 = l1
      l4 = 18*l1
      if(ntrun.le.1) l4 = 9*l1
      m0 = ntypes*ntypes+5
      m1 = nnb+5
      m2 = nbonh+5
      m3 = nbona+5
      m4 = ntheth+5
      m5 = ntheta+5
      m6 = nphih+5
      m7 = nphia+5
c
      k1 = numbnd+6
      k2 = numang+6
      k3 = numphi+6
      k4 = natyp+6
      k5 = nphb+10
      k6 = ntypes*(ntypes+1)/2+10
c
      l05 = 0
      l10 = l05+l0
      l15 = l10+l1
      l20 = l15+l1
      l25 = l20+l1
      l30 = l25+l1
      l35 = l30+3*l1
      l40 = l35+3*l1
      l45 = l40+3*l1
      l50 = l45+3*l1
      l55 = l50+l1
      l60 = l55+l1
      l65 = l60+l1
      l70 = l65+l4
c
      k02 = l70
      k04 = k02+k1
      k06 = k04+k1
      k08 = k06+k2
      k10 = k08+k2
      k12 = k10+k3
      k14 = k12+k3
      k16 = k14+k3
      k18 = k16+k4
      k20 = k18+k5
      k22 = k20+k5
      k24 = k22+k5
      k26 = k24+k6
      k28 = k26+k6
c
c     ----- nword = number of integer*2 words per double word -----
c
      i02 = (k28-1)*nword+1
      i04 = i02+l0
      i06 = i04+l1
      i08 = i06+m0
      i10 = i08+l1
      i12 = i10+m1
      i14 = i12+m2
      i16 = i14+m2
      i18 = i16+m2
      i20 = i18+m3
      i22 = i20+m3
      i24 = i22+m3
      i26 = i24+m4
      i28 = i26+m4
      i30 = i28+m4
      i32 = i30+m4
      i34 = i32+m5
      i36 = i34+m5
      i38 = i36+m5
      i40 = i38+m5
      i42 = i40+m6
      i44 = i42+m6
      i46 = i44+m6
      i48 = i46+m6
      i50 = i48+m6
      i52 = i50+m7
      i54 = i52+m7
      i56 = i54+m7
      i58 = i56+m7
      i60 = i58+m7
      i62 = i60+l2
      i64 = i62+l1
      i78 = i64+l1
      i80 = i78+l1
c
c     ----- check for the availability of memory -----
c
c  .. mfg 4/4/90
c  .. the previous core allocation was in error
c     and would not work on a 64-bit integer machine
c     here is first stab at correcting things
c     i80 now assigned arbitrary 20000 words
      lasti   = i80+20000
      maxnb =   20000
      length  = (k28-l05) + (lasti-i02+nword-1)/nword
c     now allocate total memory required; the actual offsets
c     will be determined in lcmem1
      l05 = igmem_alloc_inf(length,fnm,snm,"l05",IGMEM_DEBUG)
      end
      subroutine lcmem1
      implicit REAL  (a-h,o-z)
c
c     ----- routine to partition the memory for various
c           arrays according to their needs -----
c           first pointer has been set in lcmem0 (l05)
c
INCLUDE(common/modmin)
INCLUDE(common/modfov)
      common/memla/i02,i04,i06,i08,i10,i12,i14,i16,i18,i20,
     +             i22,i24,i26,i28,i30,i32,i34,i36,i38,i40,
     +             i42,i44,i46,i48,i50,i52,i54,i56,i58,i60,
     +             i62,i64,i66,i68,i70,i72,i74,i76,i78,i80
      common/memlb/l05,l10,l15,l20,l25,l30,l35,l40,l45,l50,
     +             l55,l60,l65,l70,l75,l80,l85,l90,l95,l100,
     +             k02,k04,k06,k08,k10,k12,k14,k16,k18,k20,
     +             k22,k24,k26,k28,k30,k32,k34,k36,k38,k40
INCLUDE(common/mod2d)
INCLUDE(common/iofile)
c
      dimension la(2),lb(2)
      equivalence (la(1),i02),(lb(1),l05)
c
c     ----- location of arrays in dynamic common -----
c           these are refered by the pointers l* -----
c
c           lbres   ...  l05
c           igraph  ...  l10
c           cg      ...  l15
c           amass   ...  l20
c           rxchrg   ...  l25
c           c       ...  l30
c           f       ...  l35
c           v       ...  l40
c           rxr      ...  l45
c           wt      ...  l50
c           isymbl  ...  l55
c           itree   ...  l60
c           work    ...  l65
c
c     ----- arrays refered by the pointers k* -----
c
c           rk      ...  k02
c           req     ...  k04
c           tk      ...  k06
c           teq     ...  k08
c           pk      ...  k10
c           pn      ...  k12
c           phase   ...  k14
c           solty   ...  k16
c           asol    ...  k18
c           bsol    ...  k20
c           hbcut   ...  k22
c           cn1     ...  k24
c           cn2     ...  k26
c
c     ----- arrays refered by pointers i* -----
c           for belly minimisation igres == irotat
c
c           ipres   ...  i02      ita     ...  i32
c           iac     ...  i04      jta     ...  i34
c           no      ...  i06      kta     ...  i36
c           iblo    ...  i08      icta    ...  i38
c           inb     ...  i10      iph     ...  i40
c           ibh     ...  i12      jph     ...  i42
c           jbh     ...  i14      kph     ...  i44
c           icbh    ...  i16      lph     ...  i46
c           iba     ...  i18      icph    ...  i48
c           jba     ...  i20      ipa     ...  i50
c           icba    ...  i22      jpa     ...  i52
c           ith     ...  i24      kpa     ...  i54
c           jth     ...  i26      lpa     ...  i56
c           kth     ...  i28      icpa    ...  i58
c           icth    ...  i30      igroup  ...  i60
c           join    ...  i62
c           irotat  ...  i64
c           iar1    ...  i78
c           iar2    ...  i80
c
c     ----- now make the partition -----
c
      l0 = nres+5
      l1 = natom+5
      l2 = 5
      if(npro.gt.0 .or. nrcp.gt.0) l2 = l1
      l4 = 18*l1
      if(ntrun.le.1) l4 = 9*l1
      m0 = ntypes*ntypes+5
      m1 = nnb+5
      m2 = nbonh+5
      m3 = nbona+5
      m4 = ntheth+5
      m5 = ntheta+5
      m6 = nphih+5
      m7 = nphia+5
c
      k1 = numbnd+6
      k2 = numang+6
      k3 = numphi+6
      k4 = natyp+6
      k5 = nphb+10
      k6 = ntypes*(ntypes+1)/2+10
c
c     l05 set in lcmem0
c
      l10 = l05+l0
      l15 = l10+l1
      l20 = l15+l1
      l25 = l20+l1
      l30 = l25+l1
      l35 = l30+3*l1
      l40 = l35+3*l1
      l45 = l40+3*l1
      l50 = l45+3*l1
      l55 = l50+l1
      l60 = l55+l1
      l65 = l60+l1
      l70 = l65+l4
c
      k02 = l70
      k04 = k02+k1
      k06 = k04+k1
      k08 = k06+k2
      k10 = k08+k2
      k12 = k10+k3
      k14 = k12+k3
      k16 = k14+k3
      k18 = k16+k4
      k20 = k18+k5
      k22 = k20+k5
      k24 = k22+k5
      k26 = k24+k6
      k28 = k26+k6
c
c     ----- nword = number of integer*2 words per double word -----
c
      i02 = (k28-1)*nword+1
      i04 = i02+l0
      i06 = i04+l1
      i08 = i06+m0
      i10 = i08+l1
      i12 = i10+m1
      i14 = i12+m2
      i16 = i14+m2
      i18 = i16+m2
      i20 = i18+m3
      i22 = i20+m3
      i24 = i22+m3
      i26 = i24+m4
      i28 = i26+m4
      i30 = i28+m4
      i32 = i30+m4
      i34 = i32+m5
      i36 = i34+m5
      i38 = i36+m5
      i40 = i38+m5
      i42 = i40+m6
      i44 = i42+m6
      i46 = i44+m6
      i48 = i46+m6
      i50 = i48+m6
      i52 = i50+m7
      i54 = i52+m7
      i56 = i54+m7
      i58 = i56+m7
      i60 = i58+m7
      i62 = i60+l2
      i64 = i62+l1
      i78 = i64+l1
      i80 = i78+l1
c
c     ----- check for the availability of memory -----
c
c     lasti   = i80+20000
c     maxnb =   20000
c     length  = (k28-l05) + (lasti-i02+nword-1)/nword
c
      return
      end
      subroutine lcmem2
      implicit REAL  (a-h,o-z)
c
c     ----- routine to determine the memory for various
c           arrays according to their needs -----
c           note that lcmem0 and lcmem1 perform the
c           same role as the original locmem
c
      common/memlb/l05,l10,l15,l20,l25,l30,l35,l40,l45,l50,
     +             l55,l60,l65,l70,l75,l80,l85,l90,l95,l100,
     +             k02,k04,k06,k08,k10,k12,k14,k16,k18,k20,
     +             k22,k24,k26,k28,k30,k32,k34,k36,k38,k40
c
      character*7 fnm
      character*6 snm
      data fnm,snm/"model.m","lcmem2"/
      call gmem_free_inf(l05,fnm,snm,"l05")
      end
      subroutine getcmn(nr,rx,ntx)
      implicit REAL  (a-h,o-z)
c
c     ----- routine to read the coordinates for minimisation -----
c
INCLUDE(common/runhed)
INCLUDE(common/iofile)
      dimension rx(*)
      nr3 = 3 * nr
c                              formatted input
      if(ntx.gt.0) go to 100
      write(iwr,3000)
      read(ntapes(1),1000) ititl1
      read(ntapes(1),6000) natom
      if(natom.ne.nr) go to 300
      read(ntapes(1),5000) (rx(i),i = 1,nr3)
      go to 200
c                              unformatted input
  100 write(iwr,4000)
      read(ntapes(1)) ititl1
      read(ntapes(1)) natom
      if(natom.ne.nr) go to 300
      read(ntapes(1)) (rx(i3),i3 = 1,nr3)
c
  200 continue
      write(iwr,2000) ititl1
      return
  300 continue
      write(iwr,7000)
      call caserr('f-mismatch')
 1000 format(20a4)
 2000 format(/4x,
     +'atomic coordinates created in model job with title:-',//4x,20a4)
 3000 format(/4x,'reading formatted ',
     +'  a t o m i c   c o o r d i n a t e s  file',/)
 4000 format(/4x,'reading unformatted (binary) ',
     +'  a t o m i c   c o o r d i n a t e s  file',/)
 5000 format(6f12.7)
 6000 format(i5)
 7000 format(t2,
     +'f-mismatch, number of coords does not match the',
     +       ' parm value')
      return
      end
      subroutine minrit(nres,rx,igraph,lbres,ipres)
      implicit REAL  (a-h,o-z)
_IF(parallel)
      logical opg_root
_ENDIF
_IFN1(cuf)      integer*2 ipres
      character *4 sol, anms
c
c     ----- routine to write final coordinates after
c           exiting from minmisation -----
c
INCLUDE(common/moda)
INCLUDE(common/modf)
INCLUDE(common/modmin)
INCLUDE(common/runhed)
INCLUDE(common/iofile)
      dimension rx(*),igraph(*),lbres(*),ipres(*),anms(5)
      data sol/'sol '/
c
      nr   = npm * nrp + nsm * nram
      nr3  = 3 * nr
      nrpt = npm * nrp
      do loop = 1,5
       anms(loop) = '    '
      enddo
_IF(parallel)
      if ( ntxo .gt. 1 .and. opg_root()) then
_ELSE
      if ( ntxo .gt. 1 ) then
_ENDIF
c
      write(ntapes(2),40) ititl
      write(ntapes(2),220) nr
      i0 = 0
      j0 = 0
      if ( nrpt .le. 0 ) go to 271
      do 270 nn = 1 , npm
         do 268 j = 1 , nres
            j0 = j0 + 1
            j1 = ipres(j)
            j2 = ipres(j+1) - 1
            do 100 k = j1 , j2
                   jj = i0 + k
                   jc = 3 * ( jj - 1 )
                 if(ntxo.lt.0) 
     +                 write(ntapes(2),292) jj,igraph(k),lbres(j),j0,
     +                                         (rx(jc+m),m = 1 , 3 )
                 if(ntxo.eq.2) 
     +                 write(ntapes(2),240) jj,igraph(k),lbres(j),j0,
     +                                         (rx(jc+m), m = 1 , 3 )
  100       continue
  268    continue
         i0 = i0 + nrp
  270 continue
  271 if ( nsm .le. 0 ) go to 281
      do 280 nn = 1 , nsm
         j0 = j0 + 1
         do 278 j = 1 , nram
            jj = i0 + j
            jc = 3 * (jj-1)
          if(ntxo.lt.0)
     +          write(ntapes(2),292)jj,anms(j),sol,j0,(rx(jc+m),m=1,3)
          if(ntxo.eq.2)
     +          write(ntapes(2),240)jj,anms(j),sol,j0,(rx(jc+m),m=1,3)
  278    continue
         i0 = i0 + nram
  280 continue
  281 if ( ntp .le. 0 ) go to 300
      write(ntapes(2),293) ( box(m) , m = 1 , 3 )
c
_IF(parallel)
      else if (ntxo.eq.1 .and. opg_root() ) then
_ELSE
      else if (ntxo.eq.1 ) then
_ENDIF
c
      write(ntapes(2)) ititl
c #### value of tt not set .....
c      set to zero to placate compiler
      tt=0.0d0
      write(ntapes(2)) nr, tt
      write(ntapes(2)) ( rx(i3) , i3 = 1 , nr3 )
      if ( ntp .gt. 0 ) write(ntapes(2)) ( box(m) , m = 1 , 3 )
c ...
c     ----- dump other information for the pdb generation -----
c ...
      write(ntapes(2)) nres, ( ipres(i) , i = 1 , nres + 1 )
      write(ntapes(2)) ( lbres(i)  , i = 1 , nres )
      write(ntapes(2)) ( igraph(i) , i = 1 , nr )
c
_IF(parallel)
      else if(opg_root()) then
_ELSE
      else
_ENDIF
c
      write(ntapes(2),40) ititl
c #### value of tt not set .....
c      set to zero to placate compiler
      tt=0.0d0
      write(ntapes(2),220) nr, tt
      write(ntapes(2),296) ( rx(loop) , loop = 1 , nr3 )
      if ( ntp .gt. 0 ) write(ntapes(2),296) ( box(m) , m = 1 , 3 )
c
      endif
  300 continue
_IF(parallel)
      if(opg_root()) rewind ntapes(2)
_ELSE
      rewind ntapes(2)
_ENDIF
      return
   40 format(20a4)
  220 format(i5,f10.5)
  240 format('atom',2x,i5,1x,a4,1x,a3,2x,i4,4x,3f8.3)
  292 format(i5,1x,2a4,1x,i5,3f8.3,3f8.4)
  293 format(3f10.5)
  296 format(6f12.7)
      end
      subroutine varadj(ix,ogauss)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ix
      logical ogauss
c
c     ----- routine to do the necessary accomodations for protein
c           belly minimisations -----
c
INCLUDE(common/modfov)
      common/memla/i02,i04,i06,i08,i10,i12,i14,i16,i18,i20,
     +             i22,i24,i26,i28,i30,i32,i34,i36,i38,i40,
     +             i42,i44,i46,i48,i50,i52,i54,i56,i58,i60,
     +             i62,i64,i66,i68,i70,i72,i74,i76,i78,i80
      common/memlb/l05,l10,l15,l20,l25,l30,l35,l40,l45,l50,
     +             l55,l60,l65,l70,l75,l80,l85,l90,l95,l100,
     +             k02,k04,k06,k08,k10,k12,k14,k16,k18,k20,
     +             k22,k24,k26,k28,k30,k32,k34,k36,k38,k40
c
      dimension ix(*)
c
c     ----- delete bonds which are in the belly alone or within
c           the ab initio part -----
c
      mbonh = nbonh
      call setbon(nbonh,mbonh,ix(i12),ix(i14),ix(i16),ix(i60),ogauss)
      call setbon(nbona,mbona,ix(i18),ix(i20),ix(i22),ix(i60),ogauss)
c
c     ----- delete the bonds which are in the belly alone -----
c
      mtheth = ntheth
      call setang(ntheth,mtheth,ix(i24),ix(i26),ix(i28),ix(i30),
     +            ix(i60),ogauss)
      call setang(ntheta,mtheta,ix(i32),ix(i34),ix(i36),ix(i38),
     +            ix(i60),ogauss)
c
c     ----- delete the dihedrals -----
c
      mphih = nphih
      call setdih(nphih,mphih,ix(i40),ix(i42),ix(i44),ix(i46),
     +            ix(i48),ix(i60),ogauss)
      call setdih(nphia,mphia,ix(i50),ix(i52),ix(i54),ix(i56),
     +            ix(i58),ix(i60),ogauss)
c
c     ----- find the total number of active atoms and residues -----
c
      call setatm(natom,nres,npro,ix(i02),ix(i60),ix(i64))
c
c     ----- all done return -----
c
      return
      end
      subroutine setbon(nb,mb,ib,jb,icb,igrp,ogauss)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ib,jb,icb,igrp
      logical ogauss,osame
      dimension ib(*),jb(*),icb(*),igrp(*)
      nba = 0
      nca = 0
      do 100 i = 1 , nb
             iat      = ib(i) / 3 + 1
             jat      = jb(i) / 3 + 1
             if ( (igrp(iat).le.0) .and. (igrp(jat).le.0) ) go to 100
             osame =  (igrp(iat).eq.igrp(jat)) .and. ogauss 
             if ( osame ) go to 100
             if ( i .gt. mb ) nca = nca + 1
             nba      = nba + 1
             ib(nba)  = ib(i)
             jb(nba)  = jb(i)
             icb(nba) = icb(i)
  100 continue
      nb = nba
      mb = nba - nca
      return
      end
      subroutine setang(nt,mt,it,jt,kt,ict,igrp,ogauss)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 it,jt,kt,ict,igrp
      logical ogauss,osame
      dimension it(*),jt(*),kt(*),ict(*),igrp(*)
c
      nta = 0
      mta = 0
      do 100 i = 1,nt
      iat = it(i)/3+1
      jat = jt(i)/3+1
      kat = kt(i)/3+1
      if(igrp(iat).le.0 .and. igrp(jat).le.0 .and. igrp(kat).le.0)
     +   go to 100
      osame = igrp(iat).eq.igrp(jat) .and. igrp(jat).eq.igrp(kat) .and.
     +       ogauss
      if(osame) go to 100
      if(i.gt.mt) mta = mta+1
      nta = nta+1
      it(nta) = it(i)
      jt(nta) = jt(i)
      kt(nta) = kt(i)
      ict(nta) = ict(i)
  100 continue
      nt = nta
      mt = nta-mta
      return
      end
      subroutine setdih(np,mp,ip,jp,kp,lp,icp,igrp,ogauss)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ip,jp,kp,lp,icp,igrp
      logical osame,ogauss
      dimension ip(*),jp(*),kp(*),lp(*),icp(*),igrp(*)
c
      npa = 0
      mpa = 0
      do 100 i = 1,np
      iat = ip(i)/3+1
      jat = jp(i)/3+1
      i4kp = kp(i)
      i4lp = lp(i)
      kat = iabs(i4kp)/3+1
      lat = iabs(i4lp)/3+1
      if(igrp(iat).le.0 .and. igrp(jat).le.0 .and. igrp(kat).le.0
     +   .and. igrp(lat).le.0) go to 100
      osame = igrp(iat).eq.igrp(jat) .and. igrp(jat).eq.igrp(kat)
     +       .and. igrp(kat).eq.igrp(lat) .and. ogauss
      if(osame) go to 100
      if(i.gt.mp) mpa = mpa+1
      npa = npa+1
      ip(npa) = ip(i)
      jp(npa) = jp(i)
      kp(npa) = kp(i)
      lp(npa) = lp(i)
      icp(npa) = icp(i)
  100 continue
      np = npa
      mp = npa-mpa
      return
      end
      subroutine setatm(nat,nres,natb,ipres,igrp,igres)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ipres,igrp,igres
      logical active
      dimension ipres(*),igrp(*),igres(*)
c
      idum = 0
      do 100 i = 1,nat
      if(igrp(i).le.0) go to 100
      idum = idum+1
  100 continue
      natb = idum
c
      do 120 i = 1,nres
      i1 = ipres(i)
      i2 = ipres(i+1)-1
      igres(i) = 0
      active = .false.
      do 140 j = i1,i2
      if(igrp(j).le.0) go to 140
      active = .true.
  140 continue
      if(active) igres(i) = 1
  120 continue
      return
      end
      subroutine genbel(natoms,natbel,mapmod,igrp,iflag)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 igrp
INCLUDE(common/sizes)
INCLUDE(common/infoa)
c ...
c     ----- routine to mark the ab initio atoms as bellys. the quantum
c           dummy atoms attached at the qm/mm junctions which are
c           flagged as negative in mapmod is ignored for energy
c           calculation depending on iflag -----
c ...
      dimension mapmod(*),igrp(*)
_IF1(cu)      call setsto(natoms,0,igrp)
_IF1(f)      call vfill(0,igrp,1,natoms)
_IFN1(cuf)      do 100 loop=1,natoms
_IFN1(cuf)  100 igrp(loop)=0
      if ( iflag ) 140,140,120
c ...
c             ----- iflag .gt. 0  all ab initio atoms
c                   are marked as bellys -----
c ...
  120 continue
      do 8000 i = 1 , nat
         iat       = iabs( mapmod(i) )
         igrp(iat) = 1
8000  continue
      natb = nat
      go to 200
c
c     ----- iflag .eq. 0 junction atoms are left -----
c
  140 natb = 0
      do 8010 i = 1 , nat
         iat       = mapmod(i)
         if ( iat .le. 0 ) go to 240
         natb      = natb + 1
         igrp(iat) = 1
 240     continue
8010  continue
      if ( iflag .lt. 0 ) go to 160
      go to 200
c
c     ----- iflag .lt. 0 the ab initio atoms compliments are taken -----
c
  160 do  8020 i = 1 , natoms
          idum = igrp(i)
          if(idum) 280,280,300
  280     igrp(i) = 1
          go to 260
  300     igrp(i) = 0
  260     continue
8020  continue
      natb = natoms - natb
  200 natbel = natb
      return
      end
      subroutine gath80(natoms,mapmod,igrp,rxp,rx,iflag)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
_IFN1(cuf)      integer*2 igrp
INCLUDE(common/sizes)
INCLUDE(common/infoa)
c
c     ----- routine to gather the belly atoms . the quantum
c           dummy atoms attached at the qm/mm junctions which are
c           flagged as negative in mapmod is ignored for energy
c           calculation depending on iflag -----
c
c           iflag  = 1   all the qm atoms are considered
c                  = 0   dummy junction atoms are ignored
c                  =-1   the compliments of iflag = 0 is considered
c
      dimension mapmod(*),igrp(*),rxp(*),rx(*)
c
      if(iflag.gt.0)  then
c
c     ----- iflag .gt. 0  all ab initio atoms are marked as bellys -----
c
       i3 = 0
       do 220 i = 1 , nat
             iat = iabs(mapmod(i))
             j3 = 3 * iat - 3
             rx(i3+1) = rxp(j3+1)
             rx(i3+2) = rxp(j3+2)
             rx(i3+3) = rxp(j3+3)
             i3 = i3+3
  220  continue
c *****************
c      write(6,290) (mapmod(loop),loop=1,nat)
c290   format(1x,' mapmod ',10i5)
c      i3=0
c      do 270 loop=1,nat
c      write(6,280)loop,rx(i3+1),rx(i3+2),rx(i3+3)
c280   format(1x,i4,3(1x,f15.7))
c270   i3=i3+3
c *****************
c
c     ----- iflag .eq. 0 junction atoms are left -----
c
      else if (iflag.eq.0) then
       do 240 i = 1 , nat
        iat = mapmod(i)
           if ( iat .gt. 0 ) then
             i3 = 3 * i - 3
             j3 = 3 * iat - 3
             rx(i3+1) = rxp(j3+1)
             rx(i3+2) = rxp(j3+2)
             rx(i3+3) = rxp(j3+3)
           endif
  240  continue
c
c     ----- iflag .lt. 0 the ab initio atoms compliments are taken -----
c
      else
       i3 = 0
       do 260 i = 1 , natoms
          if ( igrp(i) .gt. 0 ) then
             j3 = 3 * i - 3
             rx(i3+1) = rxp(j3+1)
             rx(i3+2) = rxp(j3+2)
             rx(i3+3) = rxp(j3+3)
             i3 = i3+3
        endif
  260  continue
c
      endif
      return
      end
      subroutine stuf80(natoms,mapmod,igrp,rxp,rx,iflag)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 igrp
INCLUDE(common/sizes)
INCLUDE(common/infoa)
c ...
c     ----- routine to stuff the belly atoms . the quantum
c           dummy atoms attached at the qm/mm junctions which are
c           flagged as negative in mapmod is ignored for energy
c           calculation depending on iflag -----
c
c           iflag  =  1   all the qm atoms are considered
c                  =  0   dummy junction atoms are ignored
c                  = -1   the compliments of iflag = 0 is considered
c ...
      dimension mapmod(*),igrp(*),rxp(*),rx(*)
      if ( iflag ) 160,140,120
c ...
c     ----- iflag .gt. 0  all ab initio atoms are marked as bellys -----
c ...
  120 i3 = 0
      do 220 i = 1 , nat
          iat = iabs(mapmod(i))
          j3  = 3 * iat - 3
          rxp(j3+1) = rx(i3+1)
          rxp(j3+2) = rx(i3+2)
          rxp(j3+3) = rx(i3+3)
          i3 = i3 + 3
  220 continue
      go to 200
c ...
c     ----- iflag .eq. 0 junction atoms are left -----
c ...
  140 do 240 i = 1 , nat
          iat = mapmod(i)
          if ( iat .le. 0 )  go to 240
          i3 = 3 * i - 3
          j3 = 3 * iat - 3
          rxp(j3+1) = rx(i3+1)
          rxp(j3+2) = rx(i3+2)
          rxp(j3+3) = rx(i3+3)
  240 continue
      go to 200
c ...
c     ----- iflag .lt. 0 the ab initio atoms compliments are taken -----
c ...
  160 i3 = 0
      do 260 i = 1 , natoms
          if ( igrp(i) .le. 0 ) go to 260
          j3 = 3 * i - 3
          rxp(j3+1) = rx(i3+1)
          rxp(j3+2) = rx(i3+2)
          rxp(j3+3) = rx(i3+3)
          i3 = i3 + 3
  260 continue
  200 continue
      return
      end
      subroutine cal2d(rx,ix,mapmod,nz)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
_IFN1(cuf)      integer*2 ix
      logical ogauss
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/infoa)
c
c     ----- routine to get the second derivatives of the belly
c           by numerical differentiation. at return the first and
c           second derivatives are stored starting at rx(1) -----
c
INCLUDE(common/moda)
INCLUDE(common/modf)
INCLUDE(common/modfov)
      common/memla/i02,i04,i06,i08,i10,i12,i14,i16,i18,i20,
     +             i22,i24,i26,i28,i30,i32,i34,i36,i38,i40,
     +             i42,i44,i46,i48,i50,i52,i54,i56,i58,i60,
     +             i62,i64,i66,i68,i70,i72,i74,i76,i78,i80
      common/memlb/l05,l10,l15,l20,l25,l30,l35,l40,l45,l50,
     +             l55,l60,l65,l70,l75,l80,l85,l90,l95,l100,
     +             k02,k04,k06,k08,k10,k12,k14,k16,k18,k20,
     +             k22,k24,k26,k28,k30,k32,k34,k36,k38,k40
INCLUDE(common/mod2d)
c
      dimension ene(30),ix(*),rx(*)
_IF1()      dimension ener(30)
      dimension mapmod(*)
      character *7 fnm
      character *5 snm
      data fnm/"model.m"/
      data snm/"cald2"/
c
      data half,done,rel/0.50d0,1.0d0,1.0d-06/
c
      n    = 3 * nat
      nat3 = 3 * natom
      nz3  = 3 * nz
c
c     ----- partition the work array -----
c
      ih10 = l65
      ih20 = ih10 + nat3
c     ih30 = ih20 + nat3
c
c     ----- gather the ab initio atoms in rx(l45) -----
c
      call genbel(natom,natbel,mapmod,ix(i60),1)
      npro   = natbel
      ogauss = .false.
      call varadj(ix,ogauss)
      call gath80(natom,mapmod,ix(i60),rx(l30),rx(l45),1)
c
c     ----- calculate the initial gradients -----
c
      ntnb = 1
      ndrv = 1
      k    = 0
      call bforce(rx,ix,ene,ndrv,npair,ogauss)
      call gath80(natom,mapmod,ix(i60),rx(l35),rx(ih10),1)
c
c     ----- allocate the memory for hessian and extract the forces ----
c
      length    = n*n
      ihloc = igmem_alloc_inf(length,fnm,snm,"ihloc",IGMEM_DEBUG)
c
c     ----- now calculate the hessian -----
c
      ntnb  = 0
      l45m1 = l45 - 1
      do 100 j = 1 , n
             l45m1j     = l45m1 + j
             hh         = rel
             rxhold     = rx(l45m1j)
             rx(l45m1j) = rx(l45m1j) + hh
c
c     ----- evaluate the gradients -----
c
             call stuf80(natom,mapmod,ix(i60),rx(l30),rx(l45),1)
             call bforce(rx,ix,ene,ndrv,npair,ogauss)
             call gath80(natom,mapmod,ix(i60),rx(l35),rx(ih20),1)
c
             rx(l45m1j) = rxhold
             ihloc1     = ihloc - 1
             ih101      = ih10 - 1
             ih201      = ih20 - 1
             if ( isw .ne. 1 ) then
c
c     ----- central differences -----
c
             rx(l45m1j) = rxhold - hh
c
             call stuf80(natom,mapmod,ix(i60),rx(l30),rx(l45),1)
             call bforce(rx,ix,ene,ndrv,npair,ogauss)
             call gath80(natom,mapmod,ix(i60),rx(l35),rx(ih10),1)
c
             rx(l45m1j) = rxhold
             rhh        = half / hh
             do 120 i = 1 , n
                    k = k + 1
                    rx(ihloc1+k) = ( rx(ih101+i) - rx(ih201+i) ) * rhh
  120        continue
c
             else
c
c     ----- forward differences -----
c
             rhh        = done / hh
             do 160 i = 1 , n
                    k = k + 1
                    rx(ihloc1+k) = ( rx(ih101+i) - rx(ih201+i) ) * rhh
c            write(6,1999)j,i,ih101,ih201, rx(ih101+i),rx(ih201+i)
c1999        format(1x,'j,i,ih101,ih201,rx101,rx201 ',4i5,2f12.5)
  160        continue
             endif
  100 continue
c
c     ----- if central difference calculate the original gradient -----
c
      if ( isw .ne. 1 ) then
       call stuf80(natom,mapmod,ix(i60),rx(l30),rx(l45),1)
       call bforce(rx,ix,ene,ndrv,npair,ogauss)
       call gath80(natom,mapmod,ix(i60),rx(l35),rx(ih10),1)
      endif
c
c     ----- transfer the first and second derivatives at rx(1) -----
c
      ii    = ih10 - 1
      l05m1 = l05 - 1
      do 220 i = 1 , n
             rx(l05m1+i) = -rx(ii+i)
  220 continue
      nn     = n * n
      ii     = ihloc - 1
      nz3l05 = nz3 + l05m1
      do 240 i = 1 , nn
             rx(nz3l05+i) = rx(ii+i)
  240 continue
c
c     ----- all done return -----
c
      call gmem_free_inf(ihloc,fnm,snm,"ihloc")
c
      return
      end
      subroutine nbg80(natom,nres,iblo,inb,ipres,nbact,rx,cut)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
_IFN1(cuf)      integer*2 iblo,inb,ipres
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/molsym)
c ....
c     ----- routine to generate the excluded atoms
c           list for the ab initio
c           atoms to be used in one electron integeral generation -----
c ....
INCLUDE(common/g80nb)
      dimension iblo(*),inb(*),nbact(*),ipres(*),rx(*)
c ....
c     ----- set the non bonded active array ----
c ....
      call setsto(natom,0,nbact)
c ...
c     ----- calculate the cmass of the ab initio atoms
c           and flag the ab initio atoms -----
c ...
c
       rxc = -trx
       ryc = -try
       rzc = -trz
       do 8000 i = 1 , nat
          nbact(iabs(mapmod(i))) = -1
8000   continue
_IF1()c
_IF1()c      rxc = 0.0e0
_IF1()c      ryc = 0.0e0
_IF1()c      rzc = 0.0e0
_IF1()c      do 110 i = 1 , nat
_IF1()c          iat        = iabs(mapmod(i))
_IF1()c          nbact(iat) = -1
_IF1()c          iat3       = 3 * iat - 3
_IF1()c          rxc        = rxc + rx(iat3+1)
_IF1()c          ryc        = ryc + rx(iat3+2)
_IF1()c          rzc        = rzc + rx(iat3+3)
_IF1()c  110 continue
_IF1()c      rxc = rxc / float(nat)
_IF1()c      ryc = ryc / float(nat)
_IF1()c      rzc = rzc / float(nat)
c ...
c     ----- flag those atoms whose distances are .gt. cut
c           from the center of mass of ab initio atoms -----
c ...
      do 300 i = 1 , nres
             i1 = ipres(i)
             i2 = ipres(i+1) - 1
             j3 = 3 * i1 - 3
             do 240 j  = i1 , i2
              r  = (rx(j3+1)-rxc)**2+(rx(j3+2)-ryc)**2+(rx(j3+3)-rzc)**2
              j3 = j3 + 3
              if ( r .le. cut ) go to 300
  240        continue
c
c     ----- the residue has to be skipped -----
c
             do 280 j = i1 , i2
                    nbact(j) = -1
  280        continue
  300 continue
c
c     ----- form the excluded atoms for ab initio -----
c
      nnx  = 0
      nxa = 0
      do 200 k = 1 , nat
             iat = iabs(mapmod(k))
             ipnonb(k) = nnx + 1
             idum = 0
             do 120 i = 1 , iat-1
                    ix = iblo(i)
                    do 140 j = 1 , ix
                           if ( inb(idum+j) .ne. iat ) go to 140
                           if ( nbact(i) .lt. 0 ) go to 140
                           nnx  = nnx  + 1
                           nxa = nxa + 1
                           inonb(nnx) = i
                           nbact(i)  = 1
  140               continue
                    idum = idum + ix
  120        continue
             ix = iblo(iat)
             do 160 j = 1 , ix
                    jat = inb(idum+j)
                    if( jat .le. 0 ) go to 160
                    if( nbact(jat) .lt. 0)  go to 160
                    nnx  = nnx  + 1
                    nxa = nxa + 1
                    inonb(nnx)  = jat
                    nbact(jat) = 1
  160        continue
             nnx = nnx + 1
             inonb(nnx) = 0
  200 continue
      ipnonb(nat+1) = nnx + 1
      nnbg           = nxa
      write(iwr,9008) nnbg
c      write(iwr,9018) (nbact(i),i=1,natom)
c      write(iwr,9028)
c      write(iwr,9018) (inonb(i),i=1,nnx)
      return
 9008 format(/5x,'number of atoms excluded from the ab initio set =',
     &i5,/)
c9018 format(10i5)
c9028 format(//)
      end
      subroutine zvnuc(natmod,enuc,cc,cmod,ian,czmod,nbact)
_IF1()      subroutine zvnuc(natmod,mapmod,enuc,cc,cmod,ian,czmod,nbact)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
INCLUDE(common/sizes)
      logical oskip
INCLUDE(common/infoa)
c
c     ----- routine to evaluate the energy due to the nuclear
c           and surrounding interactions -----
c
_IF1()      dimension mapmod(*),cc(*),cmod(*),ian(*),czmod(*),nbact(*)
      dimension cc(*),cmod(*),ian(*),czmod(*),nbact(*)
      dum  = 0.0d0
      iat3 = -3
      do 100 i = 1 , nat
             iat3 = iat3 + 3
             iat  = i
_IF1()             if(mapmod(i).lt.0) go to 100
             angz = dfloat(ian(i))
             gx   = cc(iat3+1)
             gy   = cc(iat3+2)
             gz   = cc(iat3+3)
             jc3  = -3
             do 120 j = 1 , natmod
                    jc3  = jc3 + 3
                    jat  = j
                    call skip80(iat,iat,jat,nbact,oskip)
                    if(oskip) go to 120
                    anvz = czmod(j)
                    cx   = gx - cmod(jc3+1)
                    cy   = gy - cmod(jc3+2)
                    cz   = gz - cmod(jc3+3)
                    rij  = dsqrt(cx*cx + cy*cy + cz*cz)
                    dum  = dum + angz * anvz / rij
  120        continue
  100 continue
      enuc = dum
      return
      end
      subroutine ztranf(rx,nvarx,okcal,idump)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
c     ----- routine to transform the cartesian second derivative to
c           second derivative. the first and second derivatives are
c           stored starting from rx(l05) on entering this routine -----
c
      common/czmat/ianz(maxnz),iz(maxnz,4),bl(maxnz),alpha(maxnz),
     & beta(maxnz),lbl(maxnz),lalpha(maxnz),lbeta(maxnz),nz,nvar
INCLUDE(common/prints)
      common/memlb/l05
INCLUDE(common/phycon)
      dimension rx(*)
      character *7 fnm
      character *6 snm
      data fnm/"model.m"/
      data snm/"ztranf"/
      data atocal/627.5251387752d0/
c ...
c     ----- evalutate some constants ----
c ...
      nparm = 3 * nz - 6
      l1    = 3 * nz
c     l2    = l1 * ( l1 + 1 ) / 2
      l3    = l1 * l1
c
c     ----- grow memory for c ...tranff -----
c
      i10 = l05
      i20 = i10 + l1
c
c     determine memory requiremnts, using final pointer tags
c
      i30   = 0
      i40   = i30 + l1
      i50   = i40 + l1
      ii10  = i50 + l3
      ii20  = ii10 + 12 * nparm
      ii30  = ii20 +  4 * nparm
      ii40  = ii30 + nparm * nparm
      ii50  = ii40 +  5 * nz
      ii60  = ii50 +  3 * nz
      ii70  = ii60 +  3 * nz
_IF(unicos)
c     allow for minv usage on cray
      i70sp = max((nparm+nparm),3*nz)
_ELSE
      i70sp = 3 * nz
_ENDIF
      ii80  = ii70 +  i70sp
      last  = ii80 +  3 * nz
c
      i30 = igmem_alloc_inf(last,fnm,snm,"i30",IGMEM_DEBUG)
c
c     now determine pointers
c
      i40   = i30 + l1
      i50   = i40 + l1
      ii10  = i50 + l3
      ii20  = ii10 + 12 * nparm
      ii30  = ii20 +  4 * nparm
      ii40  = ii30 + nparm * nparm
      ii50  = ii40 +  5 * nz
      ii60  = ii50 +  3 * nz
      ii70  = ii60 +  3 * nz
_IF(unicos)
c     allow for minv usage on cray
      i70sp = max((nparm+nparm),3*nz)
_ELSE
      i70sp = 3 * nz
_ENDIF
      ii80  = ii70 +  i70sp
      last  = ii80 +  3 * nz
c ...
c     ----- linearise the cartesian 2nd derivatives which is in
c           x(i20) -----
c ...
c     call prsq(rx(i20),nvarx,nvarx,nvarx)
      call linear(rx(i20),rx(i20),nvarx,nvarx)
      if( okcal ) then
          fca = toang(1) / atocal
          fcb = fca  * toang(1)
          nn  = nvarx * ( nvarx + 1 ) / 2
      call dscal(nvarx,fca,rx(i10),1)
      call dscal(nn,fcb,rx(i20),1)
      end if
c ...
c     ----- now call the transformation routine ----
c ...
      if(oprint(44)) idump = 2
      call tranff(nvarx,rx(i10),rx(i20),
     + rx(i30),rx(i40),rx(i50),
     + rx(ii10),rx(ii20),rx(ii30),rx(ii40),rx(ii50),
     + rx(ii60),rx(ii70),rx(ii80), nparm,idump)
c ...
c     ----- write out the second derivative matrix -----
c ...
      call squar(rx(i20),rx(i50),nparm,nparm,0)
      call prsq(rx(i50),nparm,nparm,nparm)
c
c     now transform second-derivatives over internal
c     coordinates to second-derivatives over actual z-matriz
c     variables.
c
      ii20=ii10+nparm
      ii30=ii20+nvar
      call zputf(nparm,rx(i20),rx(i50),rx(ii10),rx(ii20),
     *           rx(ii30),idump)
      call squar(rx(i50),rx(i20),nvar,nvar,0)
      if(idump.gt.0) call prsq(rx(i20),nvar,nvar,nvar)
c ...
c     ----- return the fast memory -----
c ...
      call gmem_free_inf(i30,fnm,snm,"i30")
      idump=0
      return
      end
      subroutine zputf(nparm,ffx,frcnst,ftmp1,ftmp2,fftmp,iprint)
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/czmat/ianz(maxnz),iz(maxnz,4),bl(maxnz),alpha(maxnz),
     & beta(maxnz),lbl(maxnz),lalpha(maxnz),lbeta(maxnz),nz,nvar
INCLUDE(common/iofile)
c
c     arguments:
c
c     ffx    ... input matrix of second derivatives (lower
c                triangle.)
c     frcnst ... output matrix of second derivatives over
c                actual variables.
c     ftmp1  ... scratch vector of length (nparm).
c     ftmp2  ... scratch vector of length (nvar).
c     fftmp  ... scratch vector of length (nparm*nvar).
c
      dimension ffx(*),frcnst(*)
      dimension ftmp1(*),ftmp2(*),fftmp(*)
c     data zero/0.d0/
 2001 format(32h from putff, contents of frcnst )
c
c     statement function for linear indexing.
      lind(i,j)=(i*(i-1))/2+j
c
c     clear the output array.
      nvartt=(nvar*(nvar+1))/2
      call vclr(frcnst,1,nvartt)
c
      do 40 i=1,nparm
      do 10 j=1,i
      ij=lind(i,j)
   10 ftmp1(j)=ffx(ij)
      do 20 j=i,nparm
      ij=lind(j,i)
   20 ftmp1(j)=ffx(ij)
      call putf(nz,lbl,lalpha,lbeta,nvar,ftmp1,ftmp2,iprint-1)
      do 30 j=1,nvar
      ij=i+nparm*(j-1)
   30 fftmp(ij)=ftmp2(j)
   40 continue
c
      do 80 i=1,nvar
      ii=nparm*(i-1)
      call putf(nz,lbl,lalpha,lbeta,nvar,fftmp(ii+1),ftmp2,
     $          iprint-1)
      do 60 j=1,i
      ij=lind(i,j)
   60 frcnst(ij)=frcnst(ij)+ftmp2(j)
   80 continue
c
c     possibly print frcnst.
      if(iprint.le.0)go to 90
      write(iwr,2001)
      call ltoutd(nvar,frcnst,1)
   90 continue
c
      return
      end
      subroutine tranff(nvarx,fx,ffx,ftmp1,ftmp2,fftmp,
     +  rxxi10,rxxi20,rxxi30,rxxi40,rxxi50,rxxi60,rxxi70,rxxi80,
     +           nparm,idump)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
INCLUDE(common/sizes)
c
c**********************************************************************
c     routine to transform a matrix of second derivatives over
c     cartesian coordinates to second derivatives over internal
c     coordintes.
c     this routine is to be used in conjunction with tranf.
c
c     arguments:
c
c     nparm  ... number of z-matrix degrees of freedom.
c     fx     ... input vector containing first derivatives over
c                internal coordinates.
c     ffx    ... input array contaning second derivatives over
c                cartesian coordinates.
c                at end, contains second derivatives over
c                internal coordinates.
c     ftmp1  ... scratch vector of length (3*nz).
c     ftmp2  ... scratch vector of length (nparm).
c     fftmp  ... scratch vector of length (3*nz*nparm).
c     ia     ... the matrix indexing array.
c
c     this routine needs the following space in core:
c       1.  ib     ... integer b-matrix, length 4*nparm.
c       2.  b      ... b-matrix                 12*nparm.
c       3.  g      ... g-matrix                 nparm**2.
c       4.  xm     ... scr for formbg           5*nz
c       5.  cc     ... scr for formbg           3*natoms (use nz).
c       6.  cz     ... scr for formbg           3*nz.
c       7.  ll     ... scr for formbg           nparm  (use flt pt).
c       8.  mm     ... scr for formbg           nparm  (use flt pt).
c**********************************************************************
c
      common/czmat/ianz(maxnz),iz(maxnz,4),bl(maxnz),alpha(maxnz),
     & beta(maxnz),lbl(maxnz),lalpha(maxnz),lbeta(maxnz),nz,nnvar
INCLUDE(common/phycon)
      dimension fx(*),ffx(*),ftmp1(*),ftmp2(*),fftmp(*)
      dimension rxxi10(*),rxxi20(*),rxxi30(*),rxxi40(*)
      dimension rxxi50(*),rxxi60(*),rxxi70(*),rxxi80(*)
c     data dzero/0.d0/ 
      data pt5/0.5d0/, del1/0.00001d0/
c
c     ----- initialization -----
c
      nz3   = 3 * nz
      nnprm = ( nparm * ( nparm + 1 ) ) / 2
c
c     ----- obtain initial ib,b,g matrix -----
c
      call formbg(maxnz,nz,ianz,iz,bl,alpha,beta,nparm,
     + rxxi10,rxxi20,rxxi30,rxxi40,rxxi50,rxxi60,
     + rxxi70,rxxi80,idump)
c
c     ----- transform the input matrix -----
c
c
c     ----- transform first suffix -----
c
c
      do 40 i = 1 , nvarx
c
c     ----- copy out one row into temp array ftmp1 -----
c
            ii = i*(i-1)/2
            do 10 j = 1,i
                  ij = ii+j
                  ftmp1(j) = ffx(ij)
   10       continue
            do 20 j = i,nvarx
                  ij = j*(j-1)/2+i
                  ftmp1(j) = ffx(ij)
   20       continue
c
c     ----- transform current row -----
c
       call tranf(nparm,nz,ianz,
     *    ftmp1,ftmp2,rxxi20,rxxi10,rxxi30,rxxi70)
c
c     ----- pack result into scratch array -----
c
            do 30 j = 1,nparm
                  ij = i+nz3*(j-1)
                  fftmp(ij) = ftmp2(j)
   30       continue
   40 continue
c
c     ----- transform second suffix -----
c
      call vclr(ffx,1,nnprm)
      do 80 i = 1,nparm
            ii = nz3*(i-1)
            iii=i*(i-1)/2
c
c     ----- transform current column -----
c
       call tranf(nparm,nz,ianz,
     *    fftmp(ii+1),ftmp2,rxxi20,rxxi10,rxxi30,rxxi70)
c
c     ----- pack this back into the original array (now stored symmetric
c
            do 60 j = 1,i
                  ij = iii+j
                  ffx(ij) = ffx(ij)+ftmp2(j)
   60       continue
            do 70 j = i,nparm
                  ij = j*(j-1)/2+i
                  ffx(ij) = ffx(ij)+ftmp2(j)
   70       continue
   80 continue
      call prtri(ffx,nparm)
c
c     -----  determine effect of differentiating the coordinate
c            tranformation . this is done by numerically
c            differentiating the b,ib,g matrices -----
c
c     ----- loop over the number of z-matrix degrees of freedom -----
c
      do 180 i = 1,nparm
             ii = i*(i-1)/2
c
c     ----- determine what to increment -----
c
             if(i.gt.(nz-1)) go to 110
             ipz = i+1
             jpz = 1
             del = del1
             go to 130
  110        if(i.gt.(2*nz-3)) go to 120
             ipz = i-nz+3
             jpz = 2
             del = del1
             go to 130
  120        continue
             ipz = i-2*nz+6
             jpz = 3
             del = del1
  130        continue
c
c     ----- determine step up -----
c
      delc = del
      call zmmod(bl,alpha,beta,ipz,jpz,delc,0)
      call formbg(maxnz,nz,ianz,iz,bl,alpha,beta,nparm,
     + rxxi10,rxxi20,rxxi30,rxxi40,rxxi50,rxxi60,
     + rxxi70,rxxi80,idump)
      call tranf(nparm,nz,ianz,
     *    fx,ftmp1,rxxi20,rxxi10,rxxi30,rxxi70)
c
c     ----- determine step back -----
c
      delc = -(del+del)
       call zmmod(bl,alpha,beta,ipz,jpz,delc,0)
      call formbg(maxnz,nz,ianz,iz,bl,alpha,beta,nparm,
     + rxxi10,rxxi20,rxxi30,rxxi40,rxxi50,rxxi60,
     + rxxi70,rxxi80,idump)
       call tranf(nparm,nz,ianz,
     *    fx,ftmp2,rxxi20,rxxi10,rxxi30,rxxi70)
c
c     ----- restore modified z-matrix element to its original state ---
c
            call zmmod(bl,alpha,beta,ipz,jpz,del,0)
c
c     ----- compute derivatives -----
c
      do 140 j = 1 , nparm
             ftmp1(j) = ( ftmp1(j) - ftmp2(j) ) / ( del1 + del1 )
  140 continue
* ********
c      write(6,9997)i
c9997  format(1x,'i= ',i3,1x,'ftmp1')
c      write(6,9998)(ftmp1(loop),loop=1,nparm)
c9998  format(1x,5d15.6)
* ********
c
c     ----- add contributions into ffx -----
c
      do 150 j = 1 , i
             ij = ii + j
             ffx(ij) = ffx(ij) - ftmp1(j)
  150 continue
      do 160 j = i , nparm
             ij = j*(j-1)/2 + i
             ffx(ij) = ffx(ij) - ftmp1(j)
  160 continue
  180 continue
      call dscal(nnprm,pt5,ffx,1)
      return
      end
      subroutine zmmod(bl,alpha,beta,ipz,jpz,del,idump)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
INCLUDE(common/sizes)
      common/czmat/ianz(maxnz),iz(maxnz,4),bll(3,maxnz),
     * lbl(maxnz),lalpha(maxnz),lbeta(maxnz),nz,nvar
INCLUDE(common/phycon)
INCLUDE(common/iofile)
      dimension bl(*),alpha(*),beta(*)
c
      go to(10,20,30),jpz
   10 bl(ipz) = bl(ipz) + del
      go to 40
   20 alpha(ipz) = alpha(ipz) + del
      go to 40
   30 beta(ipz) = beta(ipz) + del
   40 if(idump.gt.0) then
       call sprintxz(maxnz,nz,ianz,iz,bl,alpha,beta,toang(1),iwr)
      endif
      return
      end
      subroutine hesadj(hx,ia,ifroz,ndim,nvar)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
c
c     ----- routine to remove the frozen variables from the total
c           hessian matrix -----
c
      dimension hx(*),ia(*),ifroz(*)
c
c     ----- set up the matrix pointers -----
c
      do 100 i = 1,ndim
             ia(i) = i * ( i - 1 ) / 2
  100 continue
c
c     ----- now remove the frozen variables from the hessian -----
c
      k    = 0
      idum = 0
      do 120 i = 1 , ndim
             if ( ifroz(i) .ne. 0 ) go to 120
             idum = idum + 1
             do 140 j = 1 , i
                    if ( ifroz(j) .ne. 0 ) go to 140
                    ij    = ia(i) + j
                    k     = k + 1
                    hx(k) = hx(ij)
  140        continue
  120 continue
      nvar = idum
      return
      end
      subroutine d2xout(n,fx,ffx,okcal)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
c
c     ----- routine to output the cartesian second derivatives -----
c
INCLUDE(common/iofile)
INCLUDE(common/restar)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
INCLUDE(common/phycon)
      dimension fx(*),ffx(*)
      data m51/51/
      data atocal/627.5251387752d0/
c
c     ----- convert to atomic units if necessary -----
c
      call linear(ffx,ffx,n,n)
      if( okcal ) then
          fca = toang(1) / atocal
      call dscal(n ,fca,fx(1),1)
          fcb = fca  * toang(1)
          nn  = n * ( n + 1 ) / 2
      call dscal(nn,fcb,ffx(1),1)
      end if
      call squar(ffx,ffx,n,n,0)
      call prsq(ffx,n,n,n)
      call pusql(ffx,n,n,n)
      len1 = n * n
      len178 =lensec(len1)
      call secput(isect(478),m51,len178,ibl178)
      if ( nprint .ne. -5 ) write(iwr,1000) ibl178,len178
      call wrt3(ffx(1),len1,ibl178,idaf)
      return
 1000 format(/1x,
     +'writing dumpfile section 478 to block ',i5/1x,
     +'length of section 478               = ',i5/)
      end
      subroutine alignc(natom,mapmod,rxg80,rxmod,work,wt,rxtm)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
c ...
c  rxg80 coordinates of ab initio atoms ( from zmatrix )
c  rxmod      "      "  mechanics atoms ( including 'ab initio' atoms )
c  work       "       "   "     "   ( only 'ab initio' atoms )
c ...
c
c     ----- routine to align the mechanics coordinates along the
c           ab initio coordinates -----
c
      dimension rxg80(*),rxmod(*),mapmod(*),work(*),wt(*),rxtm(*)
      dimension rot(9),igdum(4)
c ...
c     ----- gather the ab initio coordinates  -----
c ...
      call gath80(natom,mapmod,igdum,rxmod,work,1)
c ...
c     ----- move the origin of mechanics and ab inito atoms
c           to the first atoms -----
c ...
      cx = work(1)
      cy = work(2)
      cz = work(3)
      gx = rxg80(1)
      gy = rxg80(2)
      gz = rxg80(3)
c
      i3 = 0
      do 120 i = 1 , nat
             work(i3+1)  = work(i3+1)  - cx
             work(i3+2)  = work(i3+2)  - cy
             work(i3+3)  = work(i3+3)  - cz
             rxg80(i3+1) = rxg80(i3+1) - gx
             rxg80(i3+2) = rxg80(i3+2) - gy
             rxg80(i3+3) = rxg80(i3+3) - gz
             i3          = i3 + 3
  120 continue
c
c     ----- now make a least square fit of the two sets -----
c
      do 140 i = 1 , nat
             wt(i) = 0.0d0
             if ( mapmod(i) .le. 0 ) go to 140
             wt(i) = 1.0d0
  140 continue
c
      call iterst(nat,rxg80,work,rxtm,wt,rot,esq)
c
c     ----- now transform the entire mechanics coordinate -----
c
      i3 = 0
      do 160 i = 1,natom
      rxmod(i3+1) = rxmod(i3+1)-cx
      rxmod(i3+2) = rxmod(i3+2)-cy
      rxmod(i3+3) = rxmod(i3+3)-cz
      ax = rot(1)*rxmod(i3+1)+rot(4)*rxmod(i3+2)+rot(7)*rxmod(i3+3)
      ay = rot(2)*rxmod(i3+1)+rot(5)*rxmod(i3+2)+rot(8)*rxmod(i3+3)
      az = rot(3)*rxmod(i3+1)+rot(6)*rxmod(i3+2)+rot(9)*rxmod(i3+3)
      rxmod(i3+1) = ax+gx
      rxmod(i3+2) = ay+gy
      rxmod(i3+3) = az+gz
      i3 = i3+3
  160 continue
c
c     ----- write out the coordinates -----
c
      write(iwr,9008)
      write(iwr,9018) nat
      write(iwr,9028)
      i3 = 0
      do 220 i = 1,nat
      rxg80(i3+1) = rxg80(i3+1)+gx
      rxg80(i3+2) = rxg80(i3+2)+gy
      rxg80(i3+3) = rxg80(i3+3)+gz
      kl          = mapmod(i)
      k           = iabs(kl)
      k3          = 3 * k - 3
      write(iwr,9038)i,(rxg80(i3+j),j=1,3),kl,(rxmod(k3+j),j=1,3)
      i3 = i3+3
  220 continue
      write(iwr,9048) esq
c
c     ----- all done return -----
c
      return
 9008 format(/5x,'initial alignment of mechanics and ab initio'
     * ,'coordinates.')
 9018 format(/5x,'number of ab initio atoms =',i5)
 9028 format(/26x,'ab initio',50x,'mechanics',/)
 9038 format(i5,2x,e15.8,2x,e15.8,2x,e15.8,5x,
     +       i5,2x,e15.8,2x,e15.8,2x,e15.8)
 9048 format(/5x,'square of the deviations =',e15.8)
      end
      subroutine iterst(n,x0,x1,xt,wt,rotf,esq)
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
c
c     ----- routine to rotate coordinates in x1 to give least
c          squares fit to coords in x0 -----
c
c          xt  ...  temporary storage
c          x1  ...  rotated coordinates
c
      dimension x1(3,*),x0(3,*),xt(3,*),wt(*),rotf(3,*)
      dimension aa(3,3),t(3),rot(3,3),rott(3,3)
c
c     ----- evaluate some constants -----
c
      jflag = 0
      tol = 1.0d-12
c     nm = n
      nb = n
      an = dfloat(n)
c     ainv = 1.0d0
c
c     ---- coordinate inversion for handedness comparison. ----
c
   20 continue
      if(jflag .eq. 0) go to 24
c
      do 22 j = 1, 3
      do 21 i = 1, n
      xt(j,i) = x1(j,i)
   21 continue
   22 continue
c
      do 23 i = 1, n
      x1(1,i) = -xt(1,i)
   23 continue
      do 200 i = 1,3
      do 200 j = 1,3
      rott(j,i) = rot(j,i)
  200 continue
c
c     ---- calculation of correlation matrix. ----
c
   24 continue
      do 30 i = 1, 3
      do 30 j = 1, 3
      aa(i,j) = 0.0d0
   30 continue
c
      do 32 k = 1, n
      do 31 i = 1, 3
      xx = wt(k)*x1(i,k)
      do 31 j = 1, 3
      aa(i,j) = aa(i,j) + xx*x0(j,k)
   31 continue
   32 continue
c
      do 41 i = 1, 3
      do 40 j = 1, 3
      rot(i,j) = 0.0d0
   40 continue
      rot(i,i) = 1.0d0
   41 continue
c
c     ---- from here to 70, iterative rotation scheme. ----
c
      ict = 0
      go to 51
   50 continue
      ix = ix + 1
      if(ix .lt. 4) go to 52
      if(iflag .eq. 0) go to 70
c
   51 continue
      iflag = 0
      ix = 1
c
   52 continue
      ict = ict + 1
      if(ict .gt. 200) go to 90
      iy = ix + 1
      if(iy .eq. 4) iy = 1
      iz = 6 - ix - iy
      sig = aa(iz,iy) - aa(iy,iz)
      gam = aa(iy,iy) + aa(iz,iz)
      sg = dsqrt(sig*sig + gam*gam)
      if(sg .eq. 0.0d0) go to 50
      sg = 1.0d0 / sg
      if(dabs(sig) .le. tol*dabs(gam)) go to 50
c
      do 60 k = 1, 3
      bb = gam*aa(iy,k) + sig*aa(iz,k)
      cc = gam*aa(iz,k) - sig*aa(iy,k)
      aa(iy,k) = bb*sg
      aa(iz,k) = cc*sg
      bb = gam*rot(iy,k) + sig*rot(iz,k)
      cc = gam*rot(iz,k) - sig*rot(iy,k)
      rot(iy,k) = bb*sg
      rot(iz,k) = cc*sg
   60 continue
c
      iflag = 1
      go to 50
   70 continue
c
c     ----- rotate the coordinates -----
c
      do 83 i = 1,nb
      do 81 j = 1,3
      aaa = 0.0d0
      do 80 k = 1,3
      aaa = aaa+rot(j,k)*x1(k,i)
   80 continue
      t(j) = aaa
   81 continue
      do 82 j = 1,3
      x1(j,i) = t(j)
   82 continue
   83 continue
c
   90 continue
c
c     ---- rms computation -----
c
      esq = 0.0d0
      do 100 i = 1, n
      do 100 j = 1, 3
      esq = esq + wt(i)*(x0(j,i) - x1(j,i)) ** 2
  100 continue
c
      esq = dsqrt(esq/an)
c
c     ---- handedness comparison. ----
c
      if(jflag .eq. 1) go to 110
      jflag = 1
      ert = esq
      go to 20
c
  110 continue
      if(esq .gt. ert) go to 400
c
c     ----- form the total rotation matrix -----
c
      rott(1,1) = -rott(1,1)
      rott(1,2) = -rott(1,2)
      rott(1,3) = -rott(1,3)
      do 220 i = 1,3
      do 220 j = 1,3
      dum = 0.0d0
      do 240 k = 1,3
      dum = dum+rot(j,k)*rott(k,i)
  240 continue
      rotf(j,i) = dum
  220 continue
      go to 440
c
  400 continue
      do 120 i = 1, n
      do 120 j = 1, 3
      x1(j,i) = xt(j,i)
  120 continue
      do 420 i = 1,3
      do 420 j = 1,3
      rotf(j,i) = rott(j,i)
  420 continue
      esq = ert
      write(iwr,9008)
  440 continue
      return
 9008 format(/5x,'handedness is not changed in iterst')
      end
_EXTRACT(bfast,hp800)
      subroutine bfast(rx,ix,
     +                 xt,x,fw,fg,w,igraph,lbres,ipres,igrp,erstop,
     +                 convgd,mapmod,igath)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ipres,igrp
_IFN1(cuf)      integer*2 ix
      logical erstop,convgd,skip,newstr,steep,gauss
INCLUDE(common/sizes)
c
c     ----- new conjugate gradient minimiser ------
c
INCLUDE(common/iofile)
INCLUDE(common/statis)
INCLUDE(common/moda)
INCLUDE(common/modf)
INCLUDE(common/modmin)
INCLUDE(common/modfov)
c
      dimension rx(*),ix(*)
      dimension x(*),fg(*),w(*),igraph(*),lbres(*),ipres(*)
      dimension xt(*),igrp(*),fw(*)
      dimension ene(30),mapmod(*)
c
      data maxlin,mxfcon,mstcyc/ 10, 4 , 10 /
      data dxstm,crits,stepm/1.0d-08,1.0d-09,1.0d-12/
c ...
c                   timing statistics.
c ...
      call cpuwal(begin,ebegin)
c
c     ----- evaluate some constants -----
c
      skale  = 2.0d0
      nr     = npro
      n      = 3 * nr
      ier    = 0
      ntnb   = 1
      fnq    = dsqrt(dfloat(n))
c     imes   = 0
      rms    = 0.0d0
      erstop = .false.
      convgd = .false.
      steep  = .false.
      skip   = .false.
      gauss  = .false.
      newstr = .false.
      nstcyc = 0
      fold   = 0.0d+00
      dxst   = dx0
c
c     ----- partition the working array -----
c
      irsdx  = n
      irsdg  = irsdx  + n
      iginit = irsdg  + n
      ixopt  = iginit + n
      igopt  = ixopt  + n
c
c     ----- set some parameters to begin the calculation -----
c
      iterc  = 0
      ncalls = 0
      iterfm = iterc
c
c     ----- let the initial search direction be minus the gradient
c           vector. iterrs gives the iteration number of the most
c           recent restart , but is set to zero when steepest descent
c           direction is used -----
c
c    5 continue
c ...
c                         keep check on time.
c ...
    5 call timrem(tlefts)
c
      ncalls = ncalls + 1
      if ( ncalls .eq. ncalls / nsnb * nsnb ) ntnb = 1
      if ( ntnb .eq. 1 .and. ncalls .gt. 1 ) steep = .true.
c
c     ----- calculate the force and energy -----
c
      call stuf80(natom,mapmod,igrp,xt,x,igath)
      call bforce(rx,ix,ene,1,ndum,gauss)
      call gath80(natom,mapmod,igrp,fw,fg,igath)
c
      f      = ene(1)
      fcon   = ene(nrep+1)
      ntnb   = 0
c
c     ----- print the intermediate results -----
c
      if(.not.((ncalls.eq.ncalls/ntpr*ntpr).or.ncalls.eq.1)) go to 200
      fdmax  = grdmax(n,fg)
_IF1(u)      if ( ncalls .eq. 1 ) rms = dsqrt( vecsum(fg,fg,n) ) / fnq
_IFN1(u)      if ( ncalls .eq. 1 ) rms = dnrm2( n,fg,1)  / fnq
      call prmode(nrep,ncalls,rms,fdmax,ene)
      call minrit(nres,xt,igraph,lbres,ipres)
      if(skip) go to 980
c
c     ----- do some steepest step whenever line search fails or
c           non-bonded update -----
c
  200 continue
      if ( .not. steep ) go to 400
      nstcyc = nstcyc + 1
      if(nstcyc.gt.mstcyc) go to 380
c
      sum    = ddot(n,fg,1,fg,1)
      rms    = dsqrt(sum)/fnq
      if(dxst.le.crits) dxst = dxstm
      dxst   = dxst/2.0d+00
      if(f.lt.fold) dxst = dxst*2.4d+00
      dxsth  = dxst/sum
      if(nstcyc.gt.1 .and. f.gt.fmin) go to 300
      fmin   = f
      fcmin  = fcon
      nfopt  = ncalls
      do 320 i = 1 , n
             w(ixopt+i) = x(i)
             w(igopt+i) = -fg(i)
  320 continue
  300 continue
c
c     ----- check for convergence -----
c
      if(rms.le.drms) go to 920
      if(ncalls.ne.imcyc) go to 340
      ier    = 131
      go to 900
  340 continue
      fold   = f
      do 360 i = 1 , n
             x(i)   = x(i) + dxsth * fg(i)
  360 continue
c ...
c                   check timing before starting a new cycle.
c ...
      call timrem(tlefti)
      if ( tlefti .gt. skale*(tlefts-tlefti)/ncalls ) go to 5
c ...
c                   tough. ran out of time.
c ...
      write(iwr,9300)
      go to 900
c
  380 continue
      steep  = .false.
      newstr = .true.
      nstcyc = 0
c
c     ----- start of conjugate gradient steps -----
c
  400 continue
      do 420 i = 1 , n
             fg(i)  = -fg(i)
  420 continue
c
      if ( .not. newstr .and. ncalls .ge. 2 ) go to 20
   10 do 15 i = 1 , n
            w(i)   = -fg(i)
   15 continue
      iterrs = 0
      if(newstr) iterc = 0
      if(iterc.gt.0) go to 80
c
   20 gnew  = 0.0d0
      sum   = 0.0d0
_IF1()      gnew = ambdot(n,w,fg)
_IF1()      sum  = ambdot(n,fg,fg)
      gnew  = ddot(n,w,1,fg,1)
      sum   = ddot(n,fg,1,fg,1)
      rms   = dsqrt(sum) / fnq
      if ( newstr .or. ncalls .eq. 1 ) go to 35
      fch   = f - fmin
c
c     ----- store the values of x, f and g, if they are the best that
c           have been calculated so far. test for convergence ----
c
      if ( fch ) 35, 30, 50
   30 if ( gnew / gmin .lt. -1.0d0 ) go to 45
   35 fmin  = f
      fcmin = fcon
      gsqrd = sum
      nfopt = ncalls
      do 40 i = 1 , n
            w(ixopt+i) = x(i)
            w(igopt+i) = fg(i)
   40 continue
   45 if ( rms .le. drms ) go to 920
c
c     ----- test if the value of maxcyc allows another call of funct ---
c
   50 if(ncalls.ne.imcyc) go to 55
      ier   = 131
      go to 900
   55 if(.not.newstr .and. ncalls.gt.1) go to 100
c
c     ----- set dfpr to dx0 . stmin is usually the step-length
c           of the most recent line search that gives the least value
c           of f -----
c
      dfpr = dxst*gsqrd
      stmin = dxst
      newstr = .false.
c
c     ----- begin the iteration -----
c
   80 iterc = iterc+1
c
      finit = f
      ginit = 0.0d0
      do 85 i = 1 , n
            w(iginit+i) = fg(i)
   85 continue
      ginit = ddot(n,w,1,fg,1)
      if ( ginit .ge. 0.0d0 ) go to 165
      gmin = ginit
      sbound = -1.0d0
      nfbeg = ncalls
      iretry = -1
c
      stepch = dmin1(stmin,dabs(dfpr/ginit))
      stmin = stepm
c
   90 step = stmin+stepch
      dxst = step
      work = 0.0d0
      do 95 i=1,n
      x(i) = w(ixopt+i)+stepch*w(i)
   95 work = dmax1(work,dabs(x(i)-w(ixopt+i)))
      if(work.gt.0.0d0) then
c ...
c                   check timing before starting a new cycle.
c ...
           call timrem(tlefti)
           if ( tlefti .gt. skale*(tlefts-tlefti)/ncalls ) go to 5
c ...
c                   tough. ran out of time.
c ...
           write(iwr,9300)
           go to 900
      end if
c
c     ----- terminate the line search if stepch is effectively zero ---
c
      if(ncalls.gt.nfbeg+1) go to 115
      if(dabs(gmin/ginit)-0.2d0) 170,170,115
c
  100 work = (fch+fch)/stepch-gnew-gmin
      ddspln = (gnew-gmin)/stepch
      if(ncalls.gt.nfopt) sbound = step
      if(ncalls.gt.nfopt) go to 105
      if(gmin*gnew.le.0.0d0) sbound = stmin
      stmin = step
      gmin = gnew
      stepch = -stepch
  105 if(fch.ne.0.0d0) ddspln = ddspln+(work+work)/stepch
c
c     ----- test for convergence of the line search, but force at least
c           two steps to be taken in order not to lose quadratic
c           termination -----
c
      if(gmin.eq.0.0d0) go to 170
      if(ncalls.le.nfbeg+1) go to 120
      if(dabs(gmin/ginit).le.0.2d0) go to 170
c
c     ----- apply the test that depends on the parameter maxlin -----
c
  110 if(ncalls.lt.nfopt+maxlin) go to 120
c
c     ----- possible non bonded update. make a restart -----
c
  115 write(iwr,9088)
      steep = .true.
      go to 170
c
  120 stepch = 0.5d0*(sbound-stmin)
      if(sbound.lt.-0.5d0) stepch = 9.0d0*stmin
      gspln = gmin+stepch*ddspln
      if(gmin*gspln.lt.0.0d0) stepch = stepch*gmin/(gmin-gspln)
      go to 90
c
c     ----- calculate the value of bbeta in the new direction -----
c
  125 sum = 0.0d0
_IF1()      sum = ambdot(n,fg,w(iginit+1))
      sum = ddot(n,fg,1,w(iginit+1),1)
      bbeta = (gsqrd-sum)/(gmin-ginit)
c
c     ----- test that the new search direction can be made downhill.
c           if not then try to improve the accuracy of the line
c           search -----
c
      if(dabs(bbeta*gmin).le.0.2d0*gsqrd) go to 135
      iretry = iretry+1
      if(iretry.le.0) go to 110
c
  135 if(f.lt.finit) iterfm = iterc
      if(iterc.lt.iterfm+mxfcon) go to 140
      write(iwr,9088)
      steep = .true.
      go to 170
  140 dfpr = stmin*ginit
c
c     ----- branch if a restart procedure is required due to the
c           iteration number or due to the scalar product of
c           consecutive gradients -----
c
      if(iretry.gt.0) go to 10
      if(iterrs.eq.0) go to 155
      if(iterc-iterrs.ge.n) go to 155
      if(dabs(sum).ge.0.2d0*gsqrd) go to 155
c
c     ----- calculate gama in the new search direction. gamden is
c           set by the restart procedure -----
c
      gama = 0.0d0
      sum = 0.0d0
_IF1()      gama = ambdot(n,fg,w(irsdg+1))
_IF1()      sum  = ambdot(n,fg,w(irsdx+1))
      gama = ddot(n,fg,1,w(irsdg+1),1)
      sum  = ddot(n,fg,1,w(irsdx+1),1)
      gama = gama/gamden
c
c     ----- restart if the new search direction is not sufficiently
c           downhill ----
c
      if(dabs(bbeta*gmin+gama*sum).ge.0.2d0*gsqrd) go to 155
c
c     ----- calculate the new search direction -----
c
      do 150 i=1,n
  150 w(i) = -fg(i)+bbeta*w(i)+gama*w(irsdx+i)
      go to 80
c
c     ----- apply the restart procedure -----
c
  155 gamden = gmin-ginit
      do 160 i=1,n
      w(irsdx+i) = w(i)
      w(irsdg+i) = fg(i)-w(iginit+i)
  160 w(i) = -fg(i)+bbeta*w(i)
      iterrs = iterc
      go to 80
c
c     ----- set ier to indicate that the search direction is uphill ---
c
  165 continue
      steep = .true.
      write(iwr,9088)
c
c     ----- ensure that f, x and g are optimal -----
c
  170 if(ncalls.eq.nfopt) go to 180
      f = fmin
      fcon = fcmin
      do 175 i=1,n
      x(i) = w(ixopt+i)
  175 fg(i) = w(igopt+i)
  180 continue
      if ( steep ) then
c ...
c                   check timing before starting a new cycle.
c ...
           call timrem(tlefti)
           if ( tlefti .gt. skale*(tlefts-tlefti)/ncalls ) go to 5
c ...
c                   tough. ran out of time.
c ...
           write(iwr,9300)
           go to 900
      end if
      if ( ier .eq. 0 ) go to 125
  900 continue
      if(ier.eq.129) write(iwr,9018)
      if(ier.eq.130) write(iwr,9028)
      if(ier.eq.131) write(iwr,9038)
      if(ier.eq.132) write(iwr,9048)
      erstop = .true.
      go to 940
c
  920 continue
      convgd = .true.
c
c     ----- write the final results -----
c
  940 continue
      write(iwr,9268)
      fdmax = grdmax(n,fg)
      call prmode(nrep,ncalls,rms,fdmax,ene)
c ...
c                   timing statistics.
c ...
      call timana(20)
      return
  980 continue
c ...
c                   timing statistics.
c ...
      write(iwr,9068)
      call timana(20)
      return
 9018 format('  line search abandoned ... problem with g')
 9028 format('  search direction is uphill ')
 9038 format('  maximum number of f evaluation exceeded')
 9048 format('  the value of f could not be reduced')
 9068 format(/4x,'warning ... time limit exceeded ... to be restarted')
 9088 format(/4x,' .... restarted due to linmin failure ...')
 9268 format(//20x,'final results',/)
 9300 format(//10x,44('*')//10x,'*** warning ***'/
     +10x,'mechanics optimisation has not converged yet'/
     +10x,'this job must be restarted'/10x,44('*')//)
      end
_ENDEXTRACT
      subroutine prmode(nrep,nstep,dff,fdmax,ene)
      implicit REAL  (a-h,o-z)
      dimension ene(*)
INCLUDE(common/iofile)
c
c     ----- partition energy terms. -----
c
      ebond = ene(6) + ene(7)
      eangle = ene(8) + ene(9)
      edihed = ene(10) + ene(13)
      enb = ene(2)
      eel = ene(3) + ene(5)
      ehbond = ene(4)
      enb14 = ene(11) + ene(14)
      eel14 = ene(12) + ene(15)
      econst = ene(16) + ene(17) + ene(18) + ene(19)
c
      write(iwr,9018)
      write(iwr,9028) nstep,ene(1),ene(nrep+1),dff,fdmax
      write(iwr,9038) ebond,eangle,edihed
      write(iwr,9048) enb,eel,ehbond
      write(iwr,9058) enb14,eel14,econst
c
 9018 format (//,'   nstep',5x,'energy',10x,'econstr',11x,'rms',12x,
     +        'gmax')
 9028 format(1x,i5,4(5x,1pe11.4),/)
 9038 format (1x,'bond    = ',1pe11.4,5x,'angle   = ',1pe11.4,5x,
     +        'dihed      = ',1pe11.4)
 9048 format (1x,'nonbond = ',1pe11.4,5x,'eel     = ',1pe11.4,5x,
     +        'hbond      = ',1pe11.4)
 9058 format (1x,'1-4 nb  = ',1pe11.4,5x,'1-4 eel = ',1pe11.4,5x,
     +        'constraint = ',1pe11.4,/)
      return
      end
_IF1()      REAL  function ambdot(n,dx,dy)
_IF1()      implicit REAL  (a-h,o-z)
_IF1()c
_IF1()c     ----- evaluates dot product of dx and dy -----
_IF1()c
_IF1()      dimension dx(2),dy(2)
_IF1()c
_IF1()      dum = 0.0d+00
_IF1()      if(n.le.0) go to 60
_IF1()c
_IF1()c     ----- clean-up loop so remaining vector length is a
_IF1()c           multiple of 20 -----
_IF1()c
_IF1()      m = n - (n/20)*20
_IF1()      if(m .eq. 0) go to 40
_IF1()      do 30 i = 1,m
_IF1()      dum = dum + dx(i)*dy(i)
_IF1()   30 continue
_IF1()      if(n .lt. 20) go to 60
_IF1()   40 mp1 = m + 1
_IF1()      do 50 i = mp1,n,20
_IF1()      dum = dum+dx(i)*dy(i)+dx(i+1)*dy(i+1)+dx(i+2)*dy(i+2)+
_IF1()     +          dx(i+3)*dy(i+3)+dx(i+4)*dy(i+4)+dx(i+5)*dy(i+5)+
_IF1()     +          dx(i+6)*dy(i+6)+dx(i+7)*dy(i+7)+dx(i+8)*dy(i+8)+
_IF1()     +          dx(i+9)*dy(i+9)+dx(i+10)*dy(i+10)+dx(i+11)*dy(i+11)+
_IF1()     +          dx(i+12)*dy(i+12)+dx(i+13)*dy(i+13)+dx(i+14)*dy(i+14)+
_IF1()     +          dx(i+15)*dy(i+15)+dx(i+16)*dy(i+16)+dx(i+17)*dy(i+17)+
_IF1()     +          dx(i+18)*dy(i+18)+dx(i+19)*dy(i+19)
_IF1()   50 continue
_IF1()   60 continue
_IF1()      ambdot = dum
_IF1()      return
_IF1()      end
      function grdmax(n,g)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
      dimension g(*)
      dum = 0.0d+00
      do 1 i = 1,n
             gi = dabs( g(i) )
             if ( gi .gt. dum ) dum = gi
 1    continue
      grdmax = dum
      return
      end
      subroutine bforce(x,ix,ene,ndrv,nbnum,gauss)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ix
      logical dield,gauss
c
c     ----- this routine will calculate the force (negative of first
c           derivative) and the energy -----
c
c           ndrv      flag for evaluating energy related stuff
c                = 0  only energy will be calculated
c                = 1  energy and force calculated
c
INCLUDE(common/moda)
INCLUDE(common/modf)
INCLUDE(common/modfov)
      common/memla/i02,i04,i06,i08,i10,i12,i14,i16,i18,i20,
     +             i22,i24,i26,i28,i30,i32,i34,i36,i38,i40,
     +             i42,i44,i46,i48,i50,i52,i54,i56,i58,i60,
     +             i62,i64,i66,i68,i70,i72,i74,i76,i78,i80
      common/memlb/l05,l10,l15,l20,l25,l30,l35,l40,l45,l50,
     +             l55,l60,l65,l70,l75,l80,l85,l90,l95,l100,
     +             k02,k04,k06,k08,k10,k12,k14,k16,k18,k20,
     +             k22,k24,k26,k28,k30,k32,k34,k36,k38,k40
INCLUDE(common/iofile)
INCLUDE(common/segm)
INCLUDE(common/mod2d)
INCLUDE(common/statis)
INCLUDE(common/nbterm)
      dimension ene(*),ix(*),x(*)
c ...
c                   timing analysis calls.
c ...
      call cpuwal(begin,ebegin)
c
c     ----- set the box related variables -----
c
      call mdbox(ntb,bbeta,box)
c
      nr = npm * nrp + nsm * nram
      if ( nr .ne. natom ) then
       write(iwr,428) nr,natom
       call caserr('atom mismatch in bforce')
      endif
      nr3   = 3 * nr
c     imes  = 0
      nrep  = 15
      nrepc = 5
      dield = idiel .le. 0
c
c     ----- make pairlist, when required -----
c
      if ( ntnb .ne. 0 ) then
c
c     ----- generate the pairlist depending on ntn -----
c
       if ( ntn .le. 1) then
        write(iwr,448)
        call bnbpar(natom,npair,ix(i06),ix(i04),ix(i08),ix(i10),
     +  ix(i78),ix(i80),nhb,x(l30),ix(i60),cut,ntypes,gauss)
       else if ( ntn .eq. 2 ) then
        write(iwr,458)
        call bresnr(nres,npair,ix(i02),x(l30),ix(i64),ix(i78),
     +            ix(i80),cut)
       else
        write(iwr,478)
        call bresna(nres,npair,nhb,ix(i02),ix(i04),ix(i06),
     +  ix(i08),ix(i10),ix(i78),ix(i80),x(l30),ix(i64),cut,
     +          ntypes,ix(i60),gauss)
c
       endif
       write(iwr,438) npair,nhb
       nbnum  = npair
c
c
        if(npair.gt.maxnb) then
         write(iwr,404) npair,maxnb
         call caserr('insufficient memory in bforce')
        endif
c
c     ----- zero out the energies and forces -----
c
      end if
      do 18 i = 1 , nrep+nrepc
   18 ene(i) = 0.0d+00
      do 20 i = 1 , nr3
   20 x(l35+i-1) = 0.0d0
      dum = 0.0d+00
c
c     ----- calculate the nb contributions including h-bonds -----
c
      if ( ntf .ge. 8 ) go to 30
      if ( ntn - 2 ) 34, 38, 34
   34 continue
      call nonbon(natom,ix(i78),ix(i80),ix(i04),ix(i06),x(l30),
     +      x(l35),x(k24),x(k26),x(k18),x(k20),x(k22),x(l15),x(l25),
     +      ene(2),ene(4),ene(3),dield,ntypes,ndrv)
c     write(6,1000)
c1000  format(1x,'*** nonbon')
c     write(6,1010)(x(l35+loop-1),loop=1,nr3)
c1010  format(1x,3(1x,f13.6))
      go to 30
   38 continue
      call resnbn(natom,nres,ix(i02),ix(i78),ix(i80),ix(i04),ix(i06),
     +     ix(i08),ix(i10),x(l30),x(l35),x(k24),x(k26),x(k18),x(k20),
     +     x(k22),x(l15),x(l25),ene(2),ene(4),ene(3),dield,
     +     ntypes,ndrv)
c      write(6,1020)
c1020  format(1x,'*** resnbn')
c      write(6,1010)(x(l35+loop-1),loop=1,nr3)
   30 continue
c
c     ----- calculate other contributions -----
c
      if ( ntf .gt. 8 ) go to 41
      go to (41, 42, 43, 44, 45, 46, 150, 165) ,ntf
   41 if ( nbonh .eq. 0 ) go to 42
      call bond(nbonh,ix(i12),ix(i14),ix(i16),x(k02),x(k04),ntb,
     +          x(l30),x(l35),ene(6),nbonh,dum,ndrv)
c      write(6,1030)
c1030  format(1x,'*** bond-nbonh')
c      write(6,1010)(x(l35+loop-1),loop=1,nr3)
      if ( ntf .eq. 12 ) go to 150
   42 if(nbona.eq.0) go to 43
      call bond(nbona,ix(i18),ix(i20),ix(i22),x(k02),x(k04),ntb,
     +          x(l30),x(l35),ene(7),mbona,ene(17),ndrv)
c      write(6,1040)
c1040  format(1x,'*** bond-nbona')
c      write(6,1010)(x(l35+loop-1),loop=1,nr3)
      if(ntf.eq.13) go to 150
   43 if(ntheth.eq.0) go to 44
      call angl(ntheth,ix(i24),ix(i26),ix(i28),ix(i30),x(k06),x(k08),
     +          ntb,x(l30),x(l35),ene(8),ntheth,dum,ndrv)
c      write(6,1050)
c1050  format(1x,'*** angl-ntheth')
c      write(6,1010)(x(l35+loop-1),loop=1,nr3)
      if(ntf.eq.14) go to 150
   44 if(ntheta.eq.0) go to 45
      call angl(ntheta,ix(i32),ix(i34),ix(i36),ix(i38),x(k06),x(k08),
     +          ntb,x(l30),x(l35),ene(9),mtheta,ene(18),ndrv)
c      write(6,1060)
c1060  format(1x,'*** angl-ntheta')
c      write(6,1010)(x(l35+loop-1),loop=1,nr3)
      if(ntf.eq.15) go to 150
   45 if(nphih.eq.0) go to 46
      call ephi(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48),x(k10),
     +     x(k12),x(k14),x(k24),x(k26),ntypes,x(l15),ix(i06),
     +     ix(i04),ntb,x(l30),x(l35),ene(10),ene(11),ene(12),scnb,scee,
     +     nphih,dum,dield,ndrv)
c      write(6,1070)
c1070  format(1x,'*** ephi-nphih')
c      write(6,1010)(x(l35+loop-1),loop=1,nr3)
      if(ntf.eq.16) go to 150
   46 if(nphia.eq.0) go to 150
      call ephi(nphia,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58),x(k10),
     +     x(k12),x(k14),x(k24),x(k26),ntypes,x(l15),ix(i06),
     +     ix(i04),ntb,x(l30),x(l35),ene(13),ene(14),ene(15),scnb,scee,
     +     mphia,ene(19),dield,ndrv)
c      write(6,1080)
c1080  format(1x,'*** ephi-nphia')
c      write(6,1010)(x(l35+loop-1),loop=1,nr3)
  150 continue
c
c     ----- calculate the position constraint energy -----
c
c     if(natc.gt.0)
c    +   call xconst(natc,ene(20),ix(i60),x,x(l35),x(l45),x(l50))
c
c     ----- calculate total energy -----
c
      do 152 m = 2,nrep
  152 ene(1) = ene(1)+ene(m)
      do 154 m = 2,nrepc
  154 ene(nrep+1) = ene(nrep+1)+ene(nrep+m)
      ene(1) = ene(1)+ene(20)
c
  165 continue
c
      call timana(19)
c
      return
  404 format (/4x,'too many non-bonded pairs ... stopped in force',
     +        /4x,2i10)
_IF1()  412 format(/4x,'too many dihedrals ... stopped in force',/4x,8i6)
  428 format(/4x,'number of atoms do not match in force',/4x,2i5)
  438 format(/4x,'number of non-bonded pairs     =',i10,
     +       /4x,'number of hbonded pairs        =',i10)
  448 format(/4x,'non-bonded pairs are generated in atom base')
  458 format(/4x,'non-bonded pairs are generated in residue base',
     +       /4x,'and stored as residue pairs')
  478 format(/4x,'non-bonded pairs are generated in residue base',
     +       /4x,'and are stored as atom pairs')
      end
      subroutine mdbox(ntb,beta,box)
      implicit REAL  (a-h,o-z)
      logical new
c
c     ----- routine to set the box related variables -----
c
      common/setbox/boxx(3),betax,boxh(3),boxhm,boxhm2,boxoh,boxoq,
     +              cosb,cosb2,ntm,ntbb
      dimension box(*)
      data new/.false./
c
      betax   = beta
      boxx(1) = box(1)
      boxx(2) = box(2)
      boxx(3) = box(3)
      ntbb    = ntb
      if(new) go to 3
      new     = .true.
      done    = 1.d0
      pye     = 4.d0* datan(done)
      conv    = 1.8d2/pye
      ntm     = 0
      if ( ntb ) 4, 8, 2
    2 betar   = beta/conv
      cosb    = dcos(betar)
      cosb2   = cosb+cosb
      if ( dabs(cosb) .ge. 1.d-4 ) ntm = 1
      go to 4
    3 if ( iabs(ntb) .ne. 2 ) go to 8
    4 boxhm   = box(1) * 0.5d0
      do 5 m = 1,3
           boxh(m) = box(m)*0.5d0
    5 if(boxh(m).lt.boxhm) boxhm = boxh(m)
      boxhm2  = boxhm**2
      if(ntb.gt.0) go to 8
      boxoh   = boxh(1)
      boxoq   = box(1)*0.75d0
      boxhm2  = boxhm2*0.75d0
    8 continue
      return
      end
      subroutine percon(rij2,xij)
      implicit REAL  (a-h,o-z)
      common/setbox/box(3),beta,boxh(3),boxhm,boxhm2,boxoh,boxoq,
     +              cosb,cosb2,ntm,ntb
      dimension xij(*)
c
c     ----- periodic boundary condition through minimum image
c           convention -----
c
      if(rij2.lt.boxhm2) go to 280
      if(ntm.ne.0) call traco (1,0,xij,beta,1)
      do 180 m = 1,3
      if(xij(m).lt.boxh(m)) go to 200
      xij(m) = xij(m)-box(m)
      go to 180
  200 if(xij(m).ge.-boxh(m)) go to 180
      xij(m) = xij(m)+box(m)
  180 continue
      if(ntm.ne.0) call traco (1,0,xij,beta,-1)
      if(ntb.gt.0) go to 240
      if(dabs(xij(1))+ dabs(xij(2))+ dabs(xij(3)).lt.boxoq) go to 240
      do 220 m = 1,3
  220 xij(m) = xij(m)- dsign(boxoh,xij(m))
  240 rij2 = 0.d0
      do 260 m = 1,3
  260 rij2 = rij2+xij(m)**2
c
  280 continue
      return
      end
      subroutine traco(nr,npa,x,beta,it)
      implicit REAL  (a-h,o-z)
      logical new
c
c     ----- routine to make transformation from cartesian to oblique
c           coordiantes and vice versa -----
c
      dimension x(*)
      data new/.false./
c
c     ----- initialize -----
c
      if(new) go to 8
      new   = .true.
      done  = 1.d0
      pye   = 4.d0* datan(done)
      conv  = 1.8d2/pye
      betar = beta/conv
      cosb  = dcos(betar)
      sinb  = dsin(betar)
    8 continue
c
c     ----- perform the transformation -----
c
      if(it) 41,51,31
c
   31 i3 = 3*npa+1
      do 40 j = 1,nr
      xh = x(i3+2)/sinb
      x(i3+2) = xh
      x(i3) = x(i3)-cosb*xh
      i3 = i3+3
   40 continue
      go to 51
   41 i3 = 3*npa+1
      do 50 j = 1,nr
      x(i3) = x(i3)+cosb*x(i3+2)
      x(i3+2) = x(i3+2)*sinb
      i3 = i3+3
   50 continue
   51 return
      end
      subroutine bond(nb,ib,jb,icb,rk,req,ntb,x,f,eb,mb,ecn,ndrv)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ib,jb,icb
c
c     ----- routine to get bond energy and forces for the potential
c           of cb*(b-b0)**2
c
      dimension ib(*),jb(*),icb(*),rk(*),req(*),x(*),f(*)
      dimension xij(3)
c
      eb  = 0.0d+00
      ecn = 0.0d+00
c
      do 30 n = 1 , nb
            i3   = ib(n)
            j3   = jb(n)
            ipc  = icb(n)
            rij2 = 0.0d0
            do 32 m = 1 , 3
                  xij(m) = x(i3+m) - x(j3+m)
                  rij2   = rij2 + xij(m)**2
   32       continue
            if ( ntb .eq. 0 ) go to 41
c
c     ----- call the periodic condition -----
c
            call percon(rij2,xij)
   41       continue
            rij  = dsqrt(rij2)
            db   = rij - req(ipc)
            df   = rk(ipc) * db
            ebh  = df * db
            eb   = eb + ebh
            if ( n .gt. mb ) ecn = ecn + ebh
c
            if ( ndrv .le. 0 ) go to 30
c
c     ----- update the force array -----
c
_IF1()   55       continue
            df   = 2.0d0 * df / rij
            do 90 m = 1 , 3
                  xh      = xij(m) * df
                  f(i3+m) = f(i3+m) - xh
                  f(j3+m) = f(j3+m) + xh
   90       continue
   30 continue
      return
      end
      subroutine angl(nba,it,jt,kt,ict,tk,teq,ntb,x,f,eba,mba,ecn,ndrv)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 it,jt,kt,ict
c
c     ----- routine to get the bond energies and forces for the
c           potential of the type ct*(t-t0)**2
c
      dimension it(*),jt(*),kt(*),ict(*),tk(*),teq(*),x(*),f(*)
      dimension xxij(3),xxkj(3)
c
      equivalence (xxij(1),xij),(xxij(2),yij),(xxij(3),zij)
      equivalence (xxkj(1),xkj),(xxkj(2),ykj),(xxkj(3),zkj)
c
      eba = 0.0d+00
      ecn = 0.0d+00
c
c      write(6,*)' subrotuine angl '
c
      do 300 n = 1,nba
      i3 = it(n)
      j3 = jt(n)
      k3 = kt(n)
      ic = ict(n)
c
c     ----- calculation of the bond angle -----
c
      rij2 = 0.d0
      rkj2 = 0.d0
      do 32 m = 1,3
      xj = x(j3+m)
c
c      write(6,*) 'xj ,x(i3+m) ,x(k3+m) ,i3 ,k3 ,m'
c      write(6,10)xj,x(i3+m),x(k3+m),i3,k3,m
c   10 format(' ',f15.10,5x,f15.10,5x,f15.10,5x,i4,5x,i4,5x,i4)
c
      xxij(m) = x(i3+m)-xj
      xxkj(m) = x(k3+m)-xj
      rij2 = rij2+xxij(m)**2
   32 rkj2 = rkj2+xxkj(m)**2
      if(ntb.eq.0) go to 51
c
c     ----- set the periodic boundary condition -----
c
      call percon(rij2,xxij)
      call percon(rkj2,xxkj)
   51 continue
      rij = dsqrt(rij2)
      rkj = dsqrt(rkj2)
      rrik = rij*rkj
      cst = (xij*xkj+yij*ykj+zij*zkj)/rrik
      if(cst.gt. 1.0d0) cst = 1.0d0
      if(cst.lt.-1.0d0) cst = -1.0d0
      at = dacos(cst)
c
c     ----- calculation of the energy -----
c
      da = at-teq(ic)
      df = tk(ic)*da
      ebah = df*da
      eba = eba+ebah
      if(n.gt.mba) ecn = ecn+ebah
c
c     ----- calculation of the force -----
c
      if(ndrv.le.0) go to 300
_IF1()   85 continue
      snt = dsin(at)
      if( dabs(snt).lt.1.d-07) snt = 1.d-07
      st = -(df+df)/snt
      sth = st*cst
      cik = st/rrik
      cii = sth/rij2
      ckk = sth/rkj2
      dt1 = cik*xkj-cii*xij
      dt2 = cik*ykj-cii*yij
      dt3 = cik*zkj-cii*zij
      dt7 = cik*xij-ckk*xkj
      dt8 = cik*yij-ckk*ykj
      dt9 = cik*zij-ckk*zkj
      dt4 = -dt1-dt7
      dt5 = -dt2-dt8
      dt6 = -dt3-dt9
c
      f(i3+1) = f(i3+1)-dt1
      f(i3+2) = f(i3+2)-dt2
      f(i3+3) = f(i3+3)-dt3
      f(j3+1) = f(j3+1)-dt4
      f(j3+2) = f(j3+2)-dt5
      f(j3+3) = f(j3+3)-dt6
      f(k3+1) = f(k3+1)-dt7
      f(k3+2) = f(k3+2)-dt8
      f(k3+3) = f(k3+3)-dt9
c
c     ----- finish loop over angles -----
c
  300 continue
      return
      end
      subroutine ephi(nphi,ip,jp,kp,lp,icp,pk,pn,phase,cn1,cn2,ntypes,
     +                cg,no,iac,ntb,x,f,ep,enbp,eelp,snb,see,mphi,
     +                ecn,dield,ndrv)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ip,jp,kp,lp,icp,no,iac
      logical dield
c
INCLUDE(common/iofile)
      dimension ip(*),jp(*),kp(*),lp(*),icp(*),pk(*),pn(*),phase(*)
      dimension cn1(*),cn2(*),cg(*),no(*),iac(*),x(*),f(*)
      dimension dc(6),t(6)
      dimension xxij(3),xxkj(3),xxkl(3),xxil(3)
c
      equivalence (xxij(1),xij),(xxij(2),yij),(xxij(3),zij)
      equivalence (xxkj(1),xkj),(xxkj(2),ykj),(xxkj(3),zkj)
      equivalence (xxkl(1),xkl),(xxkl(2),ykl),(xxkl(3),zkl)
c
      data pi/3.1415926d0/
c
      ep = 0.0d+00
      ecn = 0.0d+00
      enbp = 0.0d+00
      eelp = 0.0d+00
      scnb = snb
      scee = see
c
      do 200 iphi = 1,nphi
      i3 = ip(iphi)
      j3 = jp(iphi)
      k3 = kp(iphi)
      l3 = lp(iphi)
      ic = icp(iphi)
      kdiv = 1
c
c     ----- l3.lt.0 for improper torsion angles
c           k3.lt.0 for torsions with end atoms that should not
c           have nonbonded terms calculated between them -----
c
      if(l3.lt.0) kdiv = 0
      if(k3.lt.0) kdiv = 0
      k3 = iabs(k3)
      l3 = iabs(l3)
c
c     ----- calculation of the dihedral -----
c
      do 332 m = 1,3
      xj = x(j3+m)
      xk = x(k3+m)
      xxij(m) = x(i3+m)-xj
      xxkj(m) = xk-xj
      xxkl(m) = xk-x(l3+m)
  332 continue
      if(ntb.eq.0) go to 361
c
c     ----- apply periodic boundary condition -----
c
      rij2 = xxij(1)**2+xxij(2)**2+xxij(3)**2
      rkj2 = xxkj(1)**2+xxkj(2)**2+xxkj(3)**2
      rkl2 = xxkl(1)**2+xxkl(2)**2+xxkl(3)**2
      call percon(rij2,xxij)
      call percon(rkj2,xxkj)
      call percon(rkl2,xxkl)
  361 continue
c
      dx = yij*zkj-zij*ykj
      dy = zij*xkj-xij*zkj
      dz = xij*ykj-yij*xkj
      gx = zkj*ykl-ykj*zkl
      gy = xkj*zkl-zkj*xkl
      gz = ykj*xkl-xkj*ykl
c
      bi = dx*dx+dy*dy+dz*dz
      bk = gx*gx+gy*gy+gz*gz
      ct = dx*gx+dy*gy+dz*gz
c
c     ----- branch if linear dihedral -----
c
      if(bk.lt.0.01d0.or.bi.lt.0.01d0) go to 121
c
      bi = dsqrt(bi)
      bk = dsqrt(bk)
      z1 = 1.0d0/bi
      z2 = 1.0d0/bk
      ct = ct*z1*z2
      if(ct.gt.  1.0d0 ) ct =  1.0d0
      if(ct.lt.(-1.0d0)) ct = -1.0d0
      ap = dacos(ct)
c
      s = xkj*(dz*gy-dy*gz)+ykj*(dx*gz-dz*gx)+zkj*(dy*gx-dx*gy)
c
      if(s.lt.0.0d0) ap = -ap
      ap = pi-ap
      cphi = dcos(ap)
      sphi = dsin(ap)
c
c     ----- ct value above is actually -cosphi; here we change
c           its sign -----
c
      ct = -ct
c
c     ----- calculate the energy and derivatives -----
c
  120 continue
c
c     ----- get df = first der. of potential w/respect to cosphi -----
c
c           the torsional potential is assumed to have the form:
c            e = pk(ic) * (1.0+phase*cos(pn(ic)*phi)
c            where phase = 1.0 or -1.0, and pn = 1,2,3,4, or 6
c
c     ----- energy terms for dihedrals are expressed in terms of
c           cosphi in order to eliminate problem for planar angles ----
c
      iper = idint(dabs(pn(ic))+0.0001d0)
      if(iper.eq.idint(dabs(pn(ic))-0.0001d0)) go to 400
      if(iper.gt.6 .or. iper.lt.1) go to 400
      go to (30,35,40,45,50,55),iper
c
c     ----- here for iper = 1 -----
c
   30 e = ct
      df = 1.0d0
      go to 60
c
c     ----- here for iper = 2 -----
c
   35 e = 2.0d0*ct*ct-1.0d0
      df = 4.0d0*ct
      go to 60
c
c     ----- here for iper = 3 -----
c
   40 ct2 = ct*ct
      e = ct*(4.0d0*ct2-3.0d0)
      df  = 12.0d0*ct2-3.0d0
      go to 60
c
c     ----- here for iper = 4 -----
c
   45 ct2 = ct*ct
      e = 1.0d0+ct2*8.0d0*(ct2-1.0d0)
      df = 16.0d0*ct*(ct2+ct2-1.0d0)
      go to 60
   50 go to 400
c
c     ----- here for iper = 6 -----
c
   55 ct2 = ct*ct
      e = ct2*(ct2*(ct2*32.0d0-48.0d0)+18.0d0)-1.0d0
      df = ct*(ct2*(ct2*192.0d0-192.0d0)+36.0d0)
c
   60 arg = pk(ic)
      if(phase(ic).eq.0.0d0) go to 65
      arg = -arg
c
c     ----- if phase angle is other than 0 or pi then assume this angle
c           has the old energy form e = pk*( 1.0+cos(pn*phi-phase) -----
c
      if(dabs(3.14159d0-phase(ic)).lt.0.01d0) go to 65
      arg = -arg
      apn = dabs(pn(ic))
      argum = apn*ap-phase(ic)
      e = dcos(argum)
      dfp = -apn*dsin(argum)
      df = -dfp/sphi
      write(iwr,402) iphi,i3,j3,k3,l3
c
   65 e = pk(ic)+arg*e
      df = df*arg
c
      ep = ep+e
      if(iphi.gt.mphi) ecn = ecn+e
c
c     ----- skip the derivative sections if only enrergy needed -----
c
      if(ndrv.le.0) go to 301
c
c     ----- now do torsional first derivatives -----
c
c     ----- first, set up t array -----
c
      t(1) = dx
      t(2) = dy
      t(3) = dz
      t(4) = -gx
      t(5) = -gy
      t(6) = -gz
c
c     ----- now, set up array dc = first der. of cosphi w/respect
c           to the cartesian differences t -----
c
      z11 = z1*z1
      z12 = z1*z2
      z22 = z2*z2
      do 201 i = 1,3
      dc(i) = t(i+3)*z12-cphi*t(i)*z11
  201 dc(i+3) = t(i)*z12-cphi*t(i+3)*z22
c
c     ----- update the first derivative array -----
c
      dr1 = df*(dc(3)*ykj - dc(2)*zkj)
      dr2 = df*(dc(1)*zkj - dc(3)*xkj)
      dr3 = df*(dc(2)*xkj - dc(1)*ykj)
      dr4 = df*(dc(6)*ykj - dc(5)*zkj)
      dr5 = df*(dc(4)*zkj - dc(6)*xkj)
      dr6 = df*(dc(5)*xkj - dc(4)*ykj)
      drx = df*(-dc(2)*zij + dc(3)*yij + dc(5)*zkl - dc(6)*ykl)
      dry = df*( dc(1)*zij - dc(3)*xij - dc(4)*zkl + dc(6)*xkl)
      drz = df*(-dc(1)*yij + dc(2)*xij + dc(4)*ykl - dc(5)*xkl)
      f(i3+1) = f(i3+1) - dr1
      f(i3+2) = f(i3+2) - dr2
      f(i3+3) = f(i3+3) - dr3
      f(j3+1) = f(j3+1) - drx + dr1
      f(j3+2) = f(j3+2) - dry + dr2
      f(j3+3) = f(j3+3) - drz + dr3
      f(k3+1) = f(k3+1) + drx + dr4
      f(k3+2) = f(k3+2) + dry + dr5
      f(k3+3) = f(k3+3) + drz + dr6
      f(l3+1) = f(l3+1) - dr4
      f(l3+2) = f(l3+2) - dr5
      f(l3+3) = f(l3+3) - dr6
c
  301 if(pn(ic).gt.0) go to 121
      ic = ic+1
      go to 120
  121 continue
c
c     ----- compute nonbonded interactions -----
c
      if(kdiv.eq.0) go to 200
c
      ril2 = 0.0d0
      do 510 m = 1,3
      xxil(m) = x(i3+m)-x(l3+m)
      ril2 = ril2+xxil(m)**2
  510 continue
      if(ntb.eq.0) go to 520
c
c     ----- apply periodic boundary condition -----
c
      call percon(ril2,xxil)
  520 continue
c
      xd = xxil(1)
      yd = xxil(2)
      zd = xxil(3)
      ii = (i3+3)/3
      jj = (l3+3)/3
      ia1 = iac(ii)
      ia2 = iac(jj)
      index = ntypes*(ia1-1)+ia2
      ic = no(index)
c
c     ----- calculate the 14-eel energy -----
c
      r2 = 1.0d0/ril2
      if(dield) go to 540
c
c     ----- constant dielectric -----
c
      r1 = dsqrt(r2)
      g = cg(ii)*cg(jj)*r1
      eelp = eelp+g
      df2 = -g/scee
      go to 560
c
c     ----- distance depenent dielectric -----
c
  540 continue
      g = cg(ii)*cg(jj)*r2
      eelp = eelp+g
      df2 = -(g+g)/scee
  560 continue
c
c     ----- do calculate vdw for 1-4 h-bonds -----
c
      if(ic.gt.0) go to 530
      ibig = max(ia1,ia2)
      isml = min(ia1,ia2)
      ic = ibig*(ibig-1)/2+isml
  530 continue
      r6 = r2*r2*r2
      r12 = r6*r6
      f1 = cn1(ic)*r12
      f2 = cn2(ic)*r6
      e = f1-f2
      enbp = enbp+e
c
c     ----- scale the force appropriately -----
c
      df1 = (-12.0d0*f1+6.0d0*f2)/scnb
c
c     ----- skip if derivatives not needed -----
c
      if(ndrv.le.0) go to 200
c
c     ----- seperate scale factor for nb and ee -----
c
      df = (df1+df2)*r2
      xa = xd*df
      ya = yd*df
      za = zd*df
      f(i3+1) = f(i3+1)-xa
      f(i3+2) = f(i3+2)-ya
      f(i3+3) = f(i3+3)-za
      f(l3+1) = f(l3+1)+xa
      f(l3+2) = f(l3+2)+ya
      f(l3+3) = f(l3+3)+za
c
  200 continue
c
c     ----- scale of the 1-4 vdw and eel by 2 is hardwired -----
c
      enbp = enbp/scnb
      eelp = eelp/scee
      return
  400 write(iwr,402) iphi,i3,j3,k3,l3,phase(ic)
  402 format (' dihedral angle error: ',5i4,f10.4)
_IF1()  403 format(' dihedral ',4i5,'  has energy = ',f10.5)
      call caserr('dihedral angle error')
      return
      end
      subroutine nonbon(natom,iar1,iar2,iac,ico,x,f,cn1,
     +                  cn2,asol,bsol,hbcut,cg,xchrg,enb,ehb,eel,
     +                  dield,ntypes,ndrv)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 iar1,iar2,iac,ico
      logical dield
c
      common/setbox/box(3),beta,boxh(3),boxhm,boxhm2,boxoh,boxoq,
     +              cosb,cosb2,ntm,ntb
c
      dimension iar1(*),iar2(*),iac(*),ico(*),x(*),f(*),cg(*)
      dimension xchrg(*),cn1(*),cn2(*),asol(*),bsol(*),hbcut(*)
      dimension xij(3)
_IF1()      dimension xijp(3)
c
      enb = 0.0d+00
      eel = 0.0d+00
      ehb = 0.0d+00
      lpair = 0
      lim = natom-1
c
c     ----- transform cartesian to oblique if necessary -----
c
      if(ntm.ne.0) call traco(natom,0,x,beta,1)
c
      do 540 i = 1,lim
      npr = iar1(i)
      if(npr.eq.0) go to 579
      i3 = 3*i-3
c     cgi = cg(i)
      iaci = iac(i)
      jaci = ntypes*(iaci-1)
      do 541 jj = 1,npr
      lpair = lpair+1
      j = iar2(lpair)
      j3 = 3*j-3
      do 550 m = 1,3
      xij(m) = x(i3+m)-x(j3+m)
  550 continue
      s = xij(1)**2 + xij(2)**2 + xij(3)**2
      if(ntb.eq.0) go to 569
      do 562 m = 1,3
      if(xij(m) .lt. boxh(m)) go to 563
      xij(m) = xij(m)-box(m)
      go to 562
  563 if(xij(m) .ge. -boxh(m)) go to 562
      xij(m) = xij(m)+box(m)
  562 continue
_IF1()  565 continue
      s = xij(1)**2+xij(2)**2+xij(3)**2
      if(ntb) 567,569,568
  567 s = s+dmin1(0.0d0,boxoq-dabs(xij(1))-dabs(xij(2))-dabs(xij(3)))
     +          *box(1)
      go to 569
  568 if(ntm.ne.0) s = s+cosb2*xij(1)*xij(3)
  569 continue
      iacj = iac(j)
      index = jaci+iacj
      ic = ico(index)
      if(ntm.ne.0) call traco(1,0,xij,beta,-1)
c
c     ----- claculate the electrostaic energy -----
c
      r2 = 1.0d0/s
      if(dield) go to 320
c
c     ----- constant dielectric -----
c
      r1 = dsqrt(r2)
      g = cg(i)*cg(j)*r1
      eel = eel+g
      df2 = -g
      go to 340
c
c     ----- distance dependent dielctric -----
c
  320 continue
      g = cg(i)*cg(j)*r2
      eel = eel+g
      df2 = -(g+g)
  340 continue
c
c     ----- vdw energy ... branch to h-bond if ic .lt. 0 -----
c
      if(ic.le.0) go to 200
      r6 = r2*r2*r2
      r12 = r6*r6
      f1 = cn1(ic)*r12
      f2 = cn2(ic)*r6
      e = f1-f2
      enb = enb + e
      df1 = (-12.0d0*f1 + 6.0d0*f2)
      df = (df1+df2)*r2
      go to 260
c
c     ----- h-bond pairs 10-12 potential -----
c
  200 continue
      ic = iabs(ic)
      df1 = 0.0d0
      if(s.gt.hbcut(ic)) go to 220
c      r10 = r2**5
      r10 = r2*r2*r2*r2*r2
      f1 = asol(ic)*r10*r2
      f2 = bsol(ic)*r10
      ehb = ehb+f1-f2
      df1 = (-12.0d0*f1+10.0d0*f2)
  220 continue
      df = (df1+df2)*r2
c
c     ----- update the force array -----
c
  260 if(ndrv.le.0) go to 541
c
      xa = xij(1)*df
      ya = xij(2)*df
      za = xij(3)*df
      f(i3+1) = f(i3+1)-xa
      f(i3+2) = f(i3+2)-ya
      f(i3+3) = f(i3+3)-za
      f(j3+1) = f(j3+1)+xa
      f(j3+2) = f(j3+2)+ya
      f(j3+3) = f(j3+3)+za
c
  541 continue
  579 continue
  540 continue
c
c     ----- transform the oblique coordinates to cartesian if needed
c
      if(ntm.ne.0) call traco(natom,0,x,beta,-1)
      return
      end
      subroutine resdis(ia,ib,ja,jb,x,val)
      implicit REAL  (a-h,o-z)
      dimension x(*)
c
c     ----- calculate the smallest distance between the
c           residues -----
c
      dum = 1.0d+10
      do 100 i = ia , ib
             i3 = 3 * ( i - 1 )
             xi = x(i3+1)
             yi = x(i3+2)
             zi = x(i3+3)
             do 120 j = ja , jb
                    j3 = 3 * ( j - 1 )
                    xj = x(j3+1) - xi
                    yj = x(j3+2) - yi
                    zj = x(j3+3) - zi
                    r  = xj * xj + yj * yj + zj * zj
                    if ( r .lt. dum ) dum = r
  120        continue
  100 continue
      val = dum
      return
      end
      subroutine resnbn(natom,nres,ipres,iara,iarb,iac,ico,iblo,inb,
     +                  x,f,cn1,cn2,asol,bsol,hbcut,cg,xchrg,enb,ehb,
     +                  eel,dield,ntypes,ndrv)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ipres,iara,iarb,iac,ico,iblo,inb
      logical dield
c
      common/setbox/box(3),beta,boxh(3),boxhm,boxhm2,boxoh,boxoq,
     +              cosb,cosb2,ntm,ntb
c
      dimension iara(*),iarb(*),iac(*),ico(*),iblo(*),inb(*)
      dimension ipres(*),x(*),f(*),cg(*),xchrg(*),cn1(*),cn2(*)
      dimension asol(*),bsol(*),hbcut(*)
      dimension ifj(100),xij(3)
c
      enb = 0.0d+00
      eel = 0.0d+00
      ehb = 0.0d+00
      ji = 0
      jloc = 0
c
c     ----- transform cartesian to oblique if necessary -----
c
      if(ntm.ne.0) call traco(natom,0,x,beta,1)
c
c     ----- calculate nb energies and force contribution -----
c
      do 800 ii = 1, nres
      ip1 = ipres(ii)
      ip2 = ipres(ii+1)-1
      lpr = iara(ii)
      if(lpr.eq.0) go to 800
c
      do 840 i = ip1, ip2
      nx = iblo(i)
c
      do 600 j = 1, nx
  600 ifj(j) = inb(ji+j)
c
      ifj(nx+1) = 0
      ji = ji+nx
      lj = 1
      ityi = iac(i)
c
c     ----- if iac = 0 then exit the loop and pick up next i -----
c
      if(ityi .eq. 0) go to 780
      iaci = ntypes*(ityi-1)
      i3 = 3*(i-1)
c
c     ----- inner do loop to pick up atom j -----
c
      do 760 jk = 1,lpr
      jj = iarb(jloc+jk)
      jp1 = ipres(jj)
      jp2 = ipres(jj+1)-1
      if(ii .eq. jj) jp1 = i+1
c
      do 740 j = jp1, jp2
c
      if(j .ne. ifj(lj)) go to 620
      lj = lj+1
      go to 740
c
c     ----- now calculate the nb and h-bond contributions -----
c
  620 continue
      iacj = iac(j)
      if(iacj .eq. 0) go to 740
      j3 = 3*(j-1)
      index = iaci+iacj
      ic = ico(index)
c
      do 300 m = 1,3
      xij(m) = x(i3+m)-x(j3+m)
  300 continue
      s = xij(1)**2 + xij(2)**2 + xij(3)**2
      if(ntb.eq.0) go to 320
      do 340 m = 1,3
      if(xij(m) .lt. boxh(m)) go to 360
      xij(m) = xij(m)-box(m)
      go to 340
  360 if(xij(m) .ge. -boxh(m)) go to 340
      xij(m) = xij(m)+box(m)
  340 continue
_IF1()  380 continue
      s = xij(1)**2+xij(2)**2+xij(3)**2
      if(ntb) 400,320,420
  400 s = s+dmin1(0.0d0,boxoq-dabs(xij(1))-dabs(xij(2))-dabs(xij(3)))
     +          *box(1)
      go to 320
  420 if(ntm.ne.0) s = s+cosb2*xij(1)*xij(3)
  320 continue
      if(ntm.ne.0) call traco(1,0,xij,beta,-1)
c
c     ----- calculate the electrostatic energy -----
c
      r2 = 1.0d0/s
      if(dield) go to 520
c
c     ----- constant dielectric -----
c
      r1 = dsqrt(r2)
      g = cg(i)*cg(j)*r1
      eel = eel+g
      df2 = -g
      go to 540
c
c     ----- distance dependent dielectric -----
c
  520 continue
      g = cg(i)*cg(j)*r2
      eel = eel+g
      df2 = -(g+g)
  540 continue
c
c     ----- branch to h-bond if ic .lt. 0 -----
c
      if(ic.le.0) go to 200
      r6 = r2*r2*r2
      r12 = r6*r6
      f1 = cn1(ic)*r12
      f2 = cn2(ic)*r6
      e = f1-f2
      enb = enb + e
      df1 = (-12.0d0*f1 + 6.0d0*f2)
      df = (df1+df2)*r2
      go to 260
c
c     ----- h-bond pairs 10-12 potential -----
c
  200 continue
      ic = iabs(ic)
      df1 = 0.0d0
      if(s.gt.hbcut(ic)) go to 220
c      r10 = r2**5
      r10 = r2*r2*r2*r2*r2
      f1 = asol(ic)*r10*r2
      f2 = bsol(ic)*r10
      ehb = ehb+f1-f2
      df1 = (-12.0d0*f1+10.0d0*f2)
  220 continue
      df = (df1+df2)*r2
  260 continue
c
c     ----- calculate the force if needed -----
c
      if(ndrv.le.0) go to 740
c
      xa = xij(1)*df
      ya = xij(2)*df
      za = xij(3)*df
      f(i3+1) = f(i3+1)-xa
      f(i3+2) = f(i3+2)-ya
      f(i3+3) = f(i3+3)-za
      f(j3+1) = f(j3+1)+xa
      f(j3+2) = f(j3+2)+ya
      f(j3+3) = f(j3+3)+za
c
  740 continue
  760 continue
  780 continue
  840 continue
      jloc = jloc+lpr
  800 continue
c
c     ----- transform the oblique coordinates to cartesian if needed
c
      if(ntm.ne.0) call traco(natom,0,x,beta,-1)
      return
      end
      subroutine bnbpar(natom,npair,ico,iac,iblo,inb,iar1,iar2,nhb,
     +                  x,igrp,cut,ntypes,gauss)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ico,iac,iblo,inb,iar1,iar2,igrp
      logical frozi,frozj,gauss,same
_IFN1(cuf)      integer*2 ifj
c
      common/setbox/box(3),beta,boxh(3),boxhm,boxhm2,boxoh,boxoq,
     +              cosb,cosb2,ntm,ntb
c
      dimension ico(*),iac(*),iblo(*),inb(*),iar1(*),iar2(*),x(*)
      dimension ifj(100),xij(3),igrp(*)
c
c     ----- transform cartesian to oblique if necessary -----
c
      if(ntm.ne.0) call traco(natom,0,x,beta,1)
c
      lhb = 0
      lim = natom-1
      ji = 0
      lpair = 0
c
      do 100 i = 1,lim
      frozi = igrp(i).le.0
      i3 = 3*i-3
      ity1 = iac(i)
      nx = iblo(i)
      do 10 k = 1,nx
   10 ifj(k) = inb(ji+k)
      ifj(nx+1) = 0
      ji = ji+nx
      i1 = i+1
      lj = 1
      lpr = 0
      if(ity1.eq.0) go to 140
      iaci = ntypes*(ity1-1)
c
      do 40 j = i1,natom
      if(j.ne.ifj(lj)) go to 20
      lj = lj+1
      go to 40
   20 ity2 = iac(j)
      if(ity2.eq.0) go to 40
      frozj = igrp(j).le.0
c
      if(frozi .and. frozj) go to 40
c
c     ----- if g80 atoms skip it -----
c
      same = igrp(i).eq.igrp(j) .and. gauss
      if(same) go to 40
c
      j3 = 3*j-3
      xij(1) = x(i3+1)-x(j3+1)
      xij(2) = x(i3+2)-x(j3+2)
      xij(3) = x(i3+3)-x(j3+3)
c
      s = xij(1)**2 + xij(2)**2 + xij(3)**2
      if(ntb.eq.0) go to 49
      do 42 mm = 1,3
      if(xij(mm) .lt. boxh(mm)) go to 43
      xij(mm) = xij(mm)-box(mm)
      go to 42
   43 if(xij(mm) .ge. -boxh(mm)) go to 42
      xij(mm) = xij(mm)+box(mm)
   42 continue
_IF1()   45 continue
      s = xij(1)**2+xij(2)**2+xij(3)**2
      if(ntb) 47,48,49
   47 s = s+dmin1(0.0d0,boxoq-dabs(xij(1))-dabs(xij(2))-dabs(xij(3)))
     +          *box(1)
      go to 49
   48 if(ntm.ne.0) s = s+cosb2*xij(1)*xij(3)
   49 continue
c
      if(s.gt.cut) go to 40
      index = ity2+iaci
      index = ico(index)
      if(index.lt.0) lhb = lhb+1
      lpr = lpr+1
      lpair = lpair+1
      iar2(lpair) = j
   40 continue
  140 continue
      iar1(i) = lpr
  100 continue
      nhb = lhb
      npair = lpair
c
c     ----- transform the oblique to cartesian -----
c
      if(ntm.ne.0) call traco(natom,0,x,beta,-1)
      return
      end
      subroutine bresnr(nres,npair,ipres,x,igres,iara,iarb,cut)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ipres,igres,iara,iarb
      logical frozi,frozj
c
c     ----- routine to generate the non-bonded residue pairs
c           and store them sequentially -----
c
      dimension ipres(*),igres(*),iara(*),iarb(*),x(*)
c
      lpair = 0
      do 100 i = 1,nres
             iara(i) = 0
             i1      = ipres(i)
             i2      = ipres(i+1)-1
             frozi   = igres(i).le.0
             lpr     = 0
             do 120 j = i,nres
                    frozj       = igres(j).le.0
                    if(frozi .and. frozj) go to 120
                    j1          = ipres(j)
                    j2          = ipres(j+1)-1
                    call resdis(i1,i2,j1,j2,x,val)
                    if(val.gt.cut) go to 120
                    lpr         = lpr+1
                    lpair       = lpair+1
                    iarb(lpair) = j
  120        continue
             iara(i) = lpr
  100 continue
      npair = lpair
      return
      end
      subroutine bresna(nres,npair,nhb,ipres,iac,ico,iblo,inb,
     +                  iar1,iar2,x,igres,cut,ntypes,igrp,gauss)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 iar1,iar2,iac,ico,iblo,inb,ipres,igres,igrp,ifj
_IFN1(cuf)      logical*1 skip
_IF1(cuf)      logical skip
      logical frozi,frozj,same,gauss
c
      dimension iar1(*),iar2(*),iac(*),ico(*)
      dimension iblo(*),inb(*),ipres(*),x(*),igres(*),igrp(*)
      dimension ifj(30),skip(2500)
c
c     ----- calculate nb pairs -----
c
      lhb = 0
      lpair = 0
      ji = 0
c
      do 300 ii = 1,nres
      ip1 = ipres(ii)
      ip2 = ipres(ii+1)-1
      frozi = igres(ii).le.0
c
c     ----- generate the nb residues to be included for the
c           current residue ii -----
c
      do 280 jj = ii,nres
      skip(jj) = .false.
      frozj = igres(jj).le.0
      if(frozi .and. frozj) skip(jj) = .true.
      if(skip(jj)) go to 280
c
      jp1 = ipres(jj)
      jp2 = ipres(jj+1)-1
      call resdis(ip1,ip2,jp1,jp2,x,val)
      if(val.gt.cut) skip(jj) = .true.
  280 continue
c
      do 200 i = ip1,ip2
      lpr = 0
      nx = iblo(i)
c
      do 220 j = 1,nx
      ifj(j) = inb(ji+j)
  220 continue
c
      ifj(nx+1) = 0
      ji = ji+nx
      lj = 1
      ityi = iac(i)
c
c     ----- if iac = 0 then exit the loop and pick up next i -----
c
      if(ityi.eq.0) go to 180
      iaci = ntypes*(ityi-1)
c
c     ----- inner do loop to pick up atom j -----
c
      do 210 jj = ii,nres
      jp1 = ipres(jj)
      jp2 = ipres(jj+1)-1
c
      if(ii.eq.jj) jp1 = i+1
c
      do 160 j = jp1,jp2
c
      if(j.ne.ifj(lj)) go to 140
      lj = lj+1
      go to 160
c
  140 continue
c
c     ----- skip if the residue is not active -----
c
      if(skip(jj)) go to 160
c
c     ----- if g80 pairs skip it -----
c
      same = igrp(i).eq.igrp(j) .and. gauss
      if(same) go to 160
c
      iacj = iac(j)
      if(iacj.eq.0) go to 160
      index = iaci+iacj
      ic = ico(index)
c
c     ----- check whether the pair is a nb or hb -----
c
      if(ic.lt.0) lhb = lhb+1
      lpr = lpr+1
      lpair = lpair+1
      iar2(lpair) = j
  160 continue
  210 continue
  180 continue
      iar1(i) = lpr
  200 continue
  300 continue
      nhb = lhb
      npair = lpair
      return
      end
      subroutine runmd(rx,ix,natom,natbel,nres,
     +                 igraph,ipres,lbres,x,f,v,
     +                 winv,igrp,gauss)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ipres,igrp
_IFN1(cuf)      integer*2 ix
      logical gauss
c
INCLUDE(common/iofile)
INCLUDE(common/moda)
INCLUDE(common/modf)
INCLUDE(common/modmd)
INCLUDE(common/statis)
c
      dimension rx(*),ix(*)
      dimension igraph(*),ipres(*),lbres(*),x(*),winv(*),f(*),v(*)
      dimension igrp(*)
      dimension enert(70),enert2(70),ener(70),fac(4)
      call cpuwal(begin,ebegin)
c
c     ----- initialize some variables -----
c
      skale = 2.0d0
      nrpt  = natom
c     nrpt3 = 3 * nrpt
      nr    = nrpt
      nr3   = 3 * nr
      ntnb  = 1
c
c     ----- if the run is fresh get the velocity by maxwellian
c           distribution -----
c
      if(init.gt.0) go to 300
      call setvel(nr,v,winv,tempi)
      go to 320
  300 continue
_IF1(ivu)      call opfrt(ntapes(1),-1,ntx,1)
      call getmdc(nr,init,ntx,v,v)
_IF1(v)      close(unit=ntapes(4),dispose='save')
_IF1(iu)      close(unit=ntapes(4))
c      call closc(21,0)
  320 continue
_IF1(ivu)      call opfrt(ntapes(2),-1,ntxo,0)
c
c     ----- cleanup the velocity for frozen atoms -----
c
      call bellyf(nr,igrp,v)
c
      ndfp = 3*natbel
      ndf = ndfp
c
c     onet = 1.d0/3.d0
c
      tmas = 0.d0
      do 4 n = 1,nr
    4 tmas = tmas+1.0d0/winv(n)
c
      boltz2 = 8.31441d-3/2.d0
      boltz2 = boltz2/4.184d0
      fac(1) = boltz2*ndf
      fac(2) = boltz2*ndfp
      fac(3) = 1.0d-05
      if(ndfp.eq.0) fac(2) = 1.d-6
      ekinp0 = fac(2)*tempo
      ekin0 = ekinp0
      dttp = tims/tautp
      dt5 = tims/2.d0
c
      nren = 53
      nrek = 4
      nrep = 15
      nrepc = 5
      nstep = 0
c
      do 17 i = 1,nrek
   17 ener(i) = 0.d0
      ener(5) = 1.d0
      ener(6) = 1.d0
      do 18 m = 1,3
   18 ener(m+6) = box(m)
      do 19 i = 10,nren
   19 ener(i) = 0.d0
      do 20 i = 1,nren
      enert(i) = 0.d0
   20 enert2(i) = 0.d0
      nfor = 0
c
c     ----- make a first dynamics step -----
c
      if(init.eq.1) go to 51
c
      nfor = 1
      goto 100
c
   42 nfor = 0
      ntnb = 0
      i3 = 0
      ener(3) = 0.d0
      do 44 j = 1,nrp
      win = winv(j)
      winf = win*dt5
      do 43 m = 1,3
      i3 = i3+1
      ener(3) = ener(3)+v(i3)**2/win
   43 v(i3) = v(i3)-f(i3)*winf
   44 continue
      ener(3) = ener(3)/2.d0
      ener(2) = ener(3)+ener(4)
      ener(1) = ener(2)+ener(23)
      ekgp = ener(3)
      if(ekgp.lt.1.d-4) ekgp = ekinp0
c
   51 eold3 = 0.d0
      i3 = 0
      do 54 j = 1,nrp
      win = winv(j)
      do 53 m = 1,3
      i3 = i3+1
   53 eold3 = eold3+v(i3)**2/win
   54 continue
      eold3 = eold3/2.d0
      if(init.ne.1) goto 63
      ekgp = eold3
      goto 100
c
c     ----- print the initial energies and temperatures -----
c
   63 continue
_IF1()      dff = ambdot(nr3,f,f)
      dff = ddot(nr3,f,1,f,1)
      dff = dsqrt(dff/dfloat(ndfp))
      gmax = grdmax(nr3,f)
      call prmode(nrep,nstep,dff,gmax,ener(23))
      do 64 m = 2,nrek
   64 ener(m) = ener(m)/fac(m-1)
      write(iwr,9158) ener(2)
      do 65 m = 2,nrek
   65 ener(m) = ener(m)*fac(m-1)
      write(iwr,9108) nstep,t,ener(1),ener(2)
      if ( nstlim .eq. 0 ) then
c ...
c timing statistics.
c ...
            call timana(21)
            return
      end if
c
c     ----- dynamics step -----
c
  100 continue
c ...
c                   keeping watch on the time.
c ...
      call timrem(tlefts)
c
c     ----- calculate the force and cleanup the frozen part -----
c
      call bforce(rx,ix,ener(23),1,nbnum,gauss)
      call bellyf(nr,igrp,f)
c
      if(nfor.eq.1) goto 42
c
_IF1()      scalp =  dsqrt(1.d0+dttp*(ekin0/(ekgp+ekgs)-1.d0))
      scalp =  dsqrt(1.d0+dttp*(ekin0/ekgp-1.d0))
      i3 = 0
      do 133 j = 1,nrp
      win = winv(j)*tims
      do 132 m = 1,3
      i3 = i3+1
  132 v(i3) = (v(i3)+f(i3)*win)*scalp
  133 continue
      enew3 = 0.d0
      i3 = 0
      do 143 j = 1,nrp
      win = winv(j)
      do 142 m = 1,3
      i3 = i3+1
  142 enew3 = enew3+v(i3)**2/win
  143 continue
      enew3 = enew3/2.d0
      do 151 i3 = 1,nr3
  151 x(i3) = x(i3)+v(i3)*tims
      ener(3) = (eold3+enew3)/2.d0
      ener(2) = ener(3)+ener(4)
      ener(1) = ener(2)+ener(23)
      eold3 = enew3
      ekgp = enew3
c
      nstep = nstep+1
      t = t+tims
      do 210 m = 1,nren
      enert(m) = enert(m)+ener(m)
  210 enert2(m) = enert2(m)+ener(m)**2
c
      ntnb = 0
      if(nstep.eq.nstep/nsnb*nsnb) ntnb = 1
c
c     ----- output per ... steps -----
c
      if(nstep.ne.nstep/mdpr*mdpr) go to 350
_IF1()      dff = ambdot(nr3,f,f)
      dff = ddot(nr3,f,1,f,1)
      dff = dsqrt(dff/dfloat(ndfp))
      gmax = grdmax(nr3,f)
      call prmode(nrep,nstep,dff,gmax,ener(23))
      write(iwr,9108) nstep,t,ener(1),ener(2)
      call mdwrit(nr,nres,ipres,igraph,lbres,x,v,ntxo,t)
c
  350 if ( nstep .lt. nstlim ) then
           call timrem(tlefti)
           if ( tlefti .gt. skale*(tlefts-tlefti)/nstep ) go to 100
           write(iwr,9300)
      end if
c
c     ----- print averages -----
c
      tspan = nstep
      do 510 m = 1,nren
      enert(m) = enert(m)/tspan
      enert2(m) = enert2(m)/tspan-enert(m)**2
      if(enert2(m).lt.0.d0) enert2(m) = 0.d0
  510 enert2(m) =  dsqrt(enert2(m))
      write(iwr,9118)
      write(iwr,9108) nstep,t,enert(1),enert(2)
      write(iwr,9128)
      write(iwr,9108) nstep,t,enert2(1),enert2(2)
      do 520 m = 2,nrek
      enert(m) = enert(m)/fac(m-1)
  520 enert2(m) = enert2(m)/fac(m-1)
      write(iwr,9138) enert(2)
      write(iwr,9148) enert2(2)
c     temp = enert(2)
c ...
c timing statistics.
c ...
      call timana(21)
      return
 9148 format(5x,'rms fluctuation temperature =',f10.3, 'k')
 9138 format(5x,'average temperature =',f10.3,' k')
 9158 format(10x,'initial temeprature =',f10.3,' k')
 9108 format(5x,'nstep =',i5,2x,'time =',f10.3,2x,'te =',f11.2,2x,
     +       'ke =',f11.2)
 9118 format(10x,' e n e r g y  a v e r a g e s ')
 9128 format(10x,' f l u c t u a t i o n s ')
 9300 format(//10x,41('*')//10x,'***warning ***'/
     +10x,'molecular dynamics has not terminated yet',
     +10x,'this job must be restarted'/10x,41('*')//)
      end
      subroutine bellyf(nat,igrp,f)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 igrp
      dimension igrp(*),f(*)
      do 100 i = 1,nat
      if(igrp(i).gt.0) go to 100
      i3 = 3*i-3
      f(i3+1) = 0.0d0
      f(i3+2) = 0.0d0
      f(i3+3) = 0.0d0
  100 continue
      return
      end
_IF(f77)
      subroutine setvel(natom,v,winv,tempi)
      implicit REAL  (a-h,o-z)
_IF(hpux11)
      external rand
_ENDIF
c
      dimension v(*),winv(*)
c
c     ----- take the velocities from a maxwellian,
c           when required -----
c
      ig = 17737717
      iseed = ig
_IF(alliant,titan,apollo,sun)
      dumf = drand(iseed)
_ELSEIF(vms)
      dumf = decc$rand(iseed)
_ELSEIF(ksr)
      dumf = random(iseed)
_ELSEIF(rs6000)
      dumf = srand(iseed)
_ELSEIF(cray)
c cray ranf needs no seeed
_ELSE
c
c for machines without intrinsic rand(i)
c enable function source in util3.m
c
      dumf = rand(iseed)
_ENDIF
      boltz = 8.31441d-3/4.184d0
      boltz = boltz*tempi
      i = 0
      do 163 j = 1,natom
      sd =  dsqrt(boltz*winv(j))
      do 162 m = 1,3
_IF(cyber205,fps,ibm,convex,vax)
      call gauss(0.d0,sd,vg,iseed)
_ELSE
      call gauss(0.d0,sd,vg)
_ENDIF
      i = i+1
  162 v(i) = vg
  163 continue
      return
      end
c
_IF(cyber205,fps,ibm,convex,vax)
      subroutine gauss(am,sd,v,iseed)
_ELSE
      subroutine gauss(am,sd,v)
_ENDIF
      implicit REAL  (a-h,o-z)
_IF1()      REAL  ran
_IF(hpux11)
      external rand
_ENDIF
c
c     ----- routine to get normally distributed random number with
c           a given mean and standard deviation (ibm routine) -----
c
      a = 0.d0
      do 50 i = 1,12
_IF(cyber205,fps,ibm,convex,vax)
       y = ran(iseed)
_ELSEIF(cray)
       y = ranf()
_ELSEIF(alliant,titan,apollo,sun)
       y = drand(0)
_ELSEIF(vms)
       y = decc$rand(0)
_ELSEIF(ksr)
       y = random()
_ELSEIF(rs6000)
       y = rand()
_ELSE
c hp uses this to enable EXTNAME
       y = rand(0)
_ENDIF
   50 a = a+y
      v = (a-6.d0)*sd+am
      return
      end
_ELSE
!
!  version built on the F90 intrinsics
!  PS, July 07 
!
      subroutine setvel(natom,v,winv,tempi)
      implicit none

      dimension v(*),winv(*)
c
c     ----- take the velocities from a maxwellian,
c           when required -----
c
      integer i, j, m, is, natom
      real(kind=8) v, sd, vg, winv, tempi, boltz

! if in future the seed needs to be set explicitly
! this code should help
!      integer, dimension(:), allocatable :: iseed
!      call random_seed(size=is)
!      allocate(iseed(is))
!      iseed(1) = 17737717
!      call random_seed(put=iseed)

      call random_seed()

      boltz = 8.31441_8 - 3.0_8/4.184_8
      boltz = boltz*tempi
      i = 0
      do j = 1,natom
         sd =  dsqrt(boltz*winv(j))
         do m = 1,3
            call gauss(0.0_8,sd,vg)
            i = i+1
            v(i) = vg
         enddo
      enddo
      return
      end
c
      subroutine gauss(am,sd,v)
      implicit none
      real (kind=8) y,a,v,am,sd
      integer i
c
c     ----- routine to get normally distributed random number with
c           a given mean and standard deviation (ibm routine) -----
c
      a = 0.0_8
      do i = 1,12
         call random_number(y)
         a = a+y
      enddo
      v = (a-6.0_8)*sd+am
      return
      end
_ENDIF
      subroutine getmdc(natom,init,ntx,x,v)
      implicit REAL  (a-h,o-z)
INCLUDE(common/runhed)
INCLUDE(common/iofile)
      dimension x(*),v(*)
c
      if(ntx.gt.0) go to 160
      read(ntapes(1),9008) ititl
      write(iwr,9008) ititl
      read(ntapes(1),9018) matom
      if(matom.ne.natom) go to 200
      i3 = 0
      do 100 i = 1,matom
      read(ntapes(1),9028) (x(i3+j),j=1,3),(v(i3+j),j=1,3)
      i3 = i3+3
  100 continue
      go to 180
c
  160 continue
      read(ntapes(1)) ititl
      write(iwr,9008) ititl
      read(ntapes(1)) matom
      if(matom.ne.natom) go to 200
      read(ntapes(1)) (x(i3),i3=1,3*matom)
      if(init.ne.0) read(ntapes(1)) (v(i3),i3=1,3*matom)
  180 continue
      return
  200 continue
      write(iwr,9108)
      call caserr('mismatch in number of atoms')
      return
 9008 format(20a4)
 9018 format(10i5)
 9028 format(20x,3f8.3,3f8.4)
 9108 format(5x,'***** number of atoms do not match in getmdc *****')
      end
      subroutine mdwrit(natom,nres,ipres,igraph,lbres,x,v,ntxo,t)
      implicit REAL  (a-h,o-z)
_IFN1(cuf)      integer*2 ipres
_IF(parallel)
      logical opg_root
_ENDIF
INCLUDE(common/runhed)
INCLUDE(common/iofile)
      dimension ipres(*),igraph(*),lbres(*),x(*),v(*)
c
_IF(parallel)
      if(opg_root()) then
_ENDIF
      if(ntxo.gt.0) go to 200
c
c     ----- formatted output -----
c
      write(ntapes(2),9008) ititl
      write(ntapes(2),9018) natom,t
      do 100 i = 1,nres
      i1 = ipres(i)
      i2 = ipres(i+1)-1
      do 120 j = i1,i2
      j3 = 3*j-3
      write(ntapes(2),9028) j,igraph(j),lbres(i),i,
     *(x(j3+k),k=1,3),(v(j3+k),  k = 1,3)
  120 continue
  100 continue
      go to 220
c
c     ----- binary output -----
c
  200 continue
      write(ntapes(2)) ititl
      write(ntapes(2)) natom,t
      write(ntapes(2)) (x(i3),i3=1,3*natom)
      write(ntapes(2)) (v(i3),i3=1,3*natom)
  220 continue
c      call iwind(ntapes(2))
      rewind ntapes(2)
_IF(parallel)
      endif
_ENDIF
      return
 9008 format(20a4)
 9018 format(i5,f10.3)
 9028 format(i5,1x,a4,1x,a4,i5,3f8.3,3f8.4)
      end
      subroutine genzmt(coord)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
INCLUDE(common/sizes)
INCLUDE(common/runlab)
      common/czmat/ianz(maxnz),iz(maxnz,4),bl(maxnz),alpha(maxnz),
     * beta(maxnz),lbl(maxnz),lalpha(maxnz),lbeta(maxnz),nz,nvar
INCLUDE(common/csubst)
INCLUDE(common/infoa)
INCLUDE(common/infob)
INCLUDE(common/phycon)
INCLUDE(common/iofile)
INCLUDE(common/modj)
INCLUDE(common/g80nb)
      common/junk/scr(1)
c
c     ----- routine to regenerate the z-matrix
c           for the ab initio atoms -----
c
      dimension coord(*)
      data dzero/0.0d0/
c     data chgfac/18.2223d+00/
      data zblank /'        '/
      pi    = 4.0d0*datan(1.0d0)
c
c     ----- now regenerate the z matrix -----
c
c ...
c       check for presence of dummy atoms.
c ...
      do 8000 i = 1 , nz
         if ( ianz(i) .lt. 0 ) then
              write(iwr,1001)
              call caserr('invalid geometry specification in model')
         end if
8000  continue
             bl(1)    = dzero
             alpha(1) = dzero
             beta(1)  = dzero
      do 300 i = 2 , nz
c
      if(i.le.2) then
             iat      = iabs(mapmod(i))
             jat      = iabs(mapmod(iz(i,1)))
             call dismod(iat,jat,coord,dis)
             if(mapmod(i).le.0) dis = bl(i)*toang(1)
             bl(i)    = dis / toang(1)
c
      else if(i.eq.3) then
c
             iat      = iabs(mapmod(i))
             jat      = iabs(mapmod(iz(i,1)))
             kat      = iabs(mapmod(iz(i,2)))
             call dismod(iat,jat,coord,dis)
             call angmod(iat,jat,kat,coord,ang,0)
             if(mapmod(i).le.0) dis = bl(i)*toang(1)
             bl(i)    = dis/toang(1)
             alpha(i) = ang
c
      else
             iat      = iabs(mapmod(i))
             jat      = iabs(mapmod(iz(i,1)))
             kat      = iabs(mapmod(iz(i,2)))
             lat      = iabs(mapmod(iz(i,3)))
             call dismod(iat,jat,coord,dis)
             call angmod(iat,jat,kat,coord,ang,0)
             call dihmod(iat,jat,kat,lat,coord,dih,0)
             if(mapmod(i).le.0) dis = bl(i)*toang(1)
             bl(i)    = dis/toang(1)
             alpha(i) = ang
             beta(i)  = dih
      endif
  300 continue
c
c     ----- insert three dummy atoms to make sure that the all the
c           degrees of freedom of g80 atoms can be varied -----
c
      do 400 i = 1 , nz
             idum         = nz + 4 - i
             jdum         = nz + 1 - i
             iz(idum,1)   = iz(jdum,1) + 3
             iz(idum,2)   = iz(jdum,2) + 3
             iz(idum,3)   = iz(jdum,3) + 3
             iz(idum,4)   = iz(jdum,4)
             bl(idum)     = bl(jdum)
             alpha(idum)  = alpha(jdum)
             beta(idum)   = beta(jdum)
             ianz(idum)   = ianz(jdum)
c ...
c update rest of z-matrix data.
c ...
             lbl(idum)    = lbl(jdum)
             lalpha(idum) = lalpha(jdum)
             lbeta(idum)  = lbeta(jdum)
c
  400 continue
c
      if(ifree.eq.0) nvar = nvar+6
      if(ifree.eq.1) nvar = nvar+3
c
_IF1()      do i = 1 , nvar
_IF1()             idum         = nvar + 4 - i
_IF1()             jdum         = idum - 3
_IF1()             intvec(idum) = intvec(jdum)
_IF1()             fpvec(idum)  = fpvec(jdum)
_IF1()      end do
      do 8010 i = 1 , nat
             idum          = nat + 4 - i
             jdum          = idum - 3
             zaname(idum)  = zaname(jdum)
             czan(idum)    = czan(jdum)
c-weights No longer use this mass array
c             amas(idum)   = amas(jdum)
             imass(idum)     = imass(jdum)
             nuct(idum)    = nuct(jdum)
             ipseud(idum) = ipseud(jdum)
             symz(idum)    = symz(jdum)
8010  continue
      do 420 i = 1 , 3
             iz(i,1)   = 0
             iz(i,2)   = 0
             iz(i,3)   = 0
             iz(i,4)   = 0
             bl(i)     = dzero
             alpha(i)  = dzero
             beta(i)   = dzero
             ianz(i)   = -1
c ...
c update rest of z-matrix data.
c ...
             lbl(i)    = 0
             lalpha(i) = 0
             lbeta(i)  = 0
_IF1()             intvec(i) = 0
_IF1()             fpvec(i)  = dzero
             zaname(i) = zblank
c
  420 continue
c ...
c update rest of z-matrix data.
c ...
c
      iz(2,1)  = 1
      bl(2)    = 1.0d+00/toang(1)
      iz(3,1)  = 2
      iz(3,2)  = 1
      bl(3)    = bl(2)
      alpha(3) = pi/2.0d+00
      iz(4,1)  = 3
      iz(4,2)  = 2
      iz(4,3)  = 1
      iz(4,4)  = 0
      bl(4)    = bl(2)
      alpha(4) = alpha(3)
      beta(4)  = pi
      iz(5,2)  = 3
      iz(5,3)  = 2
      iz(5,4)  = 0
      alpha(5) = alpha(3)
      beta(5)  = beta(4)
      iz(6,3)  = 3
      iz(6,4)  = 0
      beta(6)  = beta(4)
c
c ... should all dummy variables be free?
c
      if(ifree.eq.0) then
      lbl(4) = nvar -5
      lalpha(4) = nvar - 4
      lbeta(4) = nvar - 3
      lalpha(5) = nvar - 2
      lbeta(5) = nvar - 1
      lbeta(6) = nvar
      endif
      if(ifree.eq.1) then
      lbl(4) = 0
      lalpha(4) = 0
      lbeta(4) = 0
      lalpha(5) = nvar - 2
      lbeta(5) = nvar - 1
      lbeta(6) = nvar
      endif
c
      call vclr(values,1,nvar)
      nz = nz + 3
      do 8020 i = 1 , nz
         lbli   = lbl(i)
         lalphi = lalpha(i)
         lbetai = lbeta(i)
         if( lbli .ne. 0 ) then
             values(lbli)   = bl(i)
         end if
         if( lalphi .ne. 0 ) then
             values(lalphi) = alpha(i)
         end if
         if( lbetai .ne. 0 ) then
             values(lbetai) = beta(i)
         end if
8020  continue
      do 8030 i = 1 , nvar
         valui     = values(i)
         cmin10(i) = valui
         cmin20(i) = valui
8030  continue
c      nz       = nz + 3
c
c     ----- print out the new z-matrix -----
c
       write(iwr,1002)
c ...
c      bond lengths are in atomic units at this stage.
c      bond angles and dihedral angles are in radians.
c      sprintxz does not change the units.
c ...
      call sprintxz(maxnz,nz,ianz,iz,bl,alpha,beta,toang(1),iwr)
      n1    = 1
      n2    = n1 + nz
      n3    = n2 + nz
      n4    = n3 + nz
      n5    = n4 + nz
      n6    = n5 + nz
      otest = .true.
      call stocxz(maxnz,nz,ianz,iz,bl,alpha,beta,otest,natom,imass,c,
     1          scr(n6),scr(n1),scr(n2),scr(n3),scr(n4),scr(n5),
     2          iwr,oerro)
      if(oerro) call caserr(
     * 'fatal error detected in z-matrix routines')
      j = 0
      do 30 i = 1 , nz
            if ( ianz(i) )30, 31, 31
 31         j         = j + 1
            zaname(j) = zaname(i)
            czan(j)   = czan(i)
 30   continue
      j = natom
      nzp1 = nz + 1
      if(nzp1.gt.nat) goto 22
      do 20 i = nzp1 , nat
            j         = j + 1
            imass(j)    = imass(i)
            zaname(j) = zaname(i)
            czan(j)   = czan(i)
            c(1,j)    = c(1,i)
            c(2,j)    = c(2,i)
            c(3,j)    = c(3,i)
 20   continue
 22   continue
      nat=j
c store contents of czan in symz and ctranz
      call dcopy(nat,czan,1,symz,1)
      call dcopy(nat,czan,1,czanr,1)
      write(iwr,101)
 101  format(/
     *30x,'coordinates (a.u.) - prior to orientation'/
     *1x,72('-')/
     *17x,'atom',16x,'x',14x,'y',14x,'z'/)
      do 102 i=1,nat
 102  write(iwr,103)i,zaname(i),(c(j,i),j=1,3)
 103  format(11x,i4,2x,a8,2x,3f15.6)
      write(iwr,104)
 104  format(/1x,72('-')/ )
      return
c
 1001 format(/4x,'dummy atoms are not permitted',
     +' in the ab initio / mechanics interface')
 1002 format(/34x,21('*'),/34x,'regenerated z-matrix.',/34x,21('*'))
c
      end
      subroutine dismod(i,j,rx,r)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
INCLUDE(common/sizes)
      dimension rx(*)
      i3  = 3 * i - 3
      j3  = 3 * j - 3
      rxd = rx(i3+1) - rx(j3+1)
      ryd = rx(i3+2) - rx(j3+2)
      rzd = rx(i3+3) - rx(j3+3)
      r   = dsqrt(rxd*rxd + ryd*ryd + rzd*rzd)
      return
      end
      subroutine angmod(i,j,k,rx,val,index)
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x), character*4 (y)
INCLUDE(common/sizes)
      dimension rx(*)
c      dimension rxa(3),rxb(3),d(3)
c
c     ----- computes angle between atoms i j k -----
c           index=0 angle returned in radians
c                =1 angle returned in degrees
c
      data convv/57.2957795131d0/
c
      i3     = 3 * i - 3
      j3     = 3 * j - 3
      k3     = 3 * k - 3
      j31    = j3 + 1
      j32    = j3 + 2
      j33    = j3 + 3
      rxa1 = rx(i3+1) - rx(j31)
      rxa2 = rx(i3+2) - rx(j32)
      rxa3 = rx(i3+3) - rx(j33)
      rxb1 = rx(k3+1) - rx(j31)
      rxb2 = rx(k3+2) - rx(j32)
      rxb3 = rx(k3+3) - rx(j33)
      dot    = rxa1*rxb1 + rxa2*rxb2 + rxa3*rxb3
      duma   = rxa1*rxa1 + rxa2*rxa2 + rxa3*rxa3
      dumb   = rxb1*rxb1 + rxb2*rxb2 + rxb3*rxb3
      dot    = dot / dsqrt(duma*dumb)
      if ( dot .gt.  1.0d0 ) dot =  1.0d0
      if ( dot .lt. -1.0d0 ) dot = -1.0d0
      val    = dacos(dot)
      if ( index .eq. 0 ) go to 100
      val    = val * convv
  100 continue
      return
      end
      subroutine dihmod(i,j,k,l,x,val,index)
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
c
c     ----- routine to calculate the torsion -----
c
      dimension x(*)
_IF1()      dimension y(2),z(2)
c
      data convv/57.2957795131d0/
c
      val = 0.0d0
      i3 = 3*i-3
      j3 = 3*j-3
      k3 = 3*k-3
      l3 = 3*l-3
      xij = x(i3+1)-x(j3+1)
      yij = x(i3+2)-x(j3+2)
      zij = x(i3+3)-x(j3+3)
      xkj = x(k3+1)-x(j3+1)
      ykj = x(k3+2)-x(j3+2)
      zkj = x(k3+3)-x(j3+3)
      xkl = x(k3+1)-x(l3+1)
      ykl = x(k3+2)-x(l3+2)
      zkl = x(k3+3)-x(l3+3)
c
      dx = yij*zkj-zij*ykj
      dy = zij*xkj-xij*zkj
      dz = xij*ykj-yij*xkj
      gx = zkj*ykl-ykj*zkl
      gy = xkj*zkl-zkj*xkl
      gz = ykj*xkl-xkj*ykl
c
      bi = dx*dx+dy*dy+dz*dz
      bk = gx*gx+gy*gy+gz*gz
      ct = dx*gx+dy*gy+dz*gz
c
c     ----- branch if linear dihedral -----
c
      if(bk.lt.0.01d0.or.bi.lt.0.01d0) go to 120
c
      ct = ct/dsqrt(bi*bk)
      if(ct.gt.  1.0d0 ) ct = 1.0d0
      if(ct.lt.(-1.0d0)) ct = -1.0d0
c
c      ----- ct = -cosphi  change its sign -----
c
      ct = -ct
      ct = dacos(ct)*convv
      sig = xij*gx+yij*gy+zij*gz
      if(sig.gt.0.d0) ct  = 3.60d+02-ct
      val = ct
      if(index.le.0) val = val/convv
      return
  120 continue
      write(iwr,9008) i,j,k,l
      return
 9008 format(5x,'warning ... linear dihedral ',4i5,/)
      end
      subroutine ver_model(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/model.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
