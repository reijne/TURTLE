c 
c  $Author: jmht $
c  $Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
c  $Locker:  $
c  $Revision: 6317 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/cphf.m,v $
c  $State: Exp $
c
c the following are omited for the parallel 2nd deriv version
c   chfndr:
c   chfdrv:
c   pfockc:
c   pdens:
c   chfcla:
c   symrhs:
c   ovlcl:
c  
      subroutine chfeqv(q,eps,au,u,work,b,cc,uu,uau,maxc,skipp,
     & npx,irmax,oconv)
      implicit real*8  (a-h,o-z)
c
c     new version of chfeq - vector algorithm -
c
      logical skipp
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
      common/blkin/g(510),nword
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      common/small/alpha(50,50),aa(50,50),wks1(50),
     + wks2(50),iblu(50),iblut(50),iblau(50),iatms(100)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      character*10 charwall
      dimension q(*)
      dimension eps(mn),au(mn,npx),u(mn,npx),b(maxc,npx)
      dimension skipp(100),cc(maxc,npx),uu(maxc,npx)
      dimension work(irmax,mn),uau(maxc,maxc,npx)
c
c     OCONV indicates whether a linear system has converged. 
c     Its dimension is taken to be equal to that of skipp.
c
c     logical oconv(100)
      logical oconv(*)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      data smal/1.0d-13/
      data tich/1.0d-24/
c
c      iterative solution of simultaneous equations to give
c      coupled hartree-fock first order wavefunction
c
      if (maxc.gt.50) maxc = 50
      uconv = 10.0d0**(-iconvv)
      lenblk = lensec(mn)
c
c      rhs and solution stored on scratchfile
c      rhs beginning at block iblkb
c      sol. beginning at block iblku
c
      iblkb = iblks + npstar*lenblk
      iblast = iblks + lenblk*np*2
      iblku = iblkb + lenblk*np
      iblk4 = jblk(1)
      idev4 = nofile(1)
      iaa = 50
      ifail = 0
      write(iwr,6000) cpulft(1) ,charwall()
      if (oprn(12)) write (iwr,6010)
      if (oprn(13)) write (iwr,6020)
c      initialise arrays
c
      nprang = npfin - npstar
      call vclr(uau,1,maxc*maxc*nprang)
      call vclr(u,1,mn*nprang)
      call vclr(b,1,maxc*nprang)
      call vclr(au,1,mn*nprang)
      call vclr(uu,1,maxc*nprang)
      npi = 0
      do 20 ko = npstar + 1 , npfin
c
c      get rhs for this perturbation
c
         oconv(ko) = .false.
         if (skipp(ko)) then
            if(odebug(2)) write (iwr,6030) ko
         else
            npi = npi + 1
            call rdedx(au(1,npi),mn,iblkb,ifils)
            if (oprn(12)) write (iwr,6040) ko , (au(i,npi),i=1,mn)
            iatms(npi) = (ko-1)/3+1
         end if
         iblkb = iblkb + lenblk
 20   continue
c
c      get zeroth order estimate
c
      if (npi.eq.0) go to 290
      do 40 ko = 1 , npi
         do 30 i = 1 , mn
            v = -au(i,ko)*eps(i)
            if (dabs(v).le.tich) v = 0.0d0
            b(1,ko) = b(1,ko) + v*v
            u(i,ko) = v
 30      continue
         if (b(1,ko).le.0.0d0) oconv(ko) = .true.
 40   continue
      iblu(1) = iblast
      ibadd = lensec(mn*npi)
      iblast = iblast + ibadd
      call wrt3(u,mn*npi,iblu(1),ifils)
c
c      start of iterative solution of chf equations
c      50 iterations are allowed ---  usually less than 10
c      are necessary
c
      if(oprn(6)) then
       write(iwr,6100)
      endif
      do 270 no = 1 , maxc
         nox = no
c
c      read in the combinations of the 2-electron integrals
c      corresponding to the hessian ( 'a-matrix' )
c
         call rdedx(u,mn*npi,iblu(no),ifils)
         call vclr(au,1,mn*npi)
         call search(iblk4,idev4)
c
c      a-matrix on file idev4=nofile(1)=secondary mainfile= ed4 (default
c      starting block jblk(1) = 1 (default)
c
         ifi = 1
         ila = min(irmax,mn)
         n0 = 1
c
c     get the first block and unpack labels
c
         call find(idev4)
         if (lcpf .or. cicv .or. cicx) then
 50         call vclr(work,1,mn*irmax)
            call get2a(work,ifi,ila,mn,idev4)
            call mxmau2(work,u,au,mn,ifi,ila,ila-ifi+1,npi,irmax)
            if (ifi.le.mn) go to 50
         else
 60         call vclr(work,1,mn*irmax)
            call get1a(work,ifi,ila,mn,n0,idev4)
            call mxmau1(work,u,au,mn,ifi,ila,ila-ifi+1,npi,irmax)
            if (ifi.le.mn) go to 60
         end if
c
c     add DFT contributions
c
         call au_dft(q,u,au,npi,iatms)
c
c     scale au by difference of eigenvalues
c
         do 80 nn = 1 , npi
            do 70 i = 1 , mn
               au(i,nn) = au(i,nn)*eps(i)
 70         continue
 80      continue
c
c      store au for this iteration
c
         iblau(no) = iblast
         iblast = iblast + ibadd
         call wrt3(au,mn*npi,iblau(no),ifils)
c
c      uu = dot product of u vectors
c
         do 90 nn = 1 , npi
            uu(no,nn) = ddot(mn,u(1,nn),1,u(1,nn),1)
 90      continue
c
c      uau is dot product of u with au
c      last column of uau
c      read au for each iteration
c
         do 110 noo = 1 , no
            call rdedx(au,mn*npi,iblau(noo),ifils)
            do 100 nn = 1 , npi
               uau(no,noo,nn) = ddot(mn,u(1,nn),1,au(1,nn),1)
 100        continue
 110     continue
c
c      last row of uau
c      read u for each iteration
c
         do 130 noo = 1 , no - 1
            call rdedx(u,mn*npi,iblu(noo),ifils)
            do 120 nn = 1 , npi
               uau(noo,no,nn) = ddot(mn,u(1,nn),1,au(1,nn),1)
 120        continue
 130     continue
c
c      nag routine to solve a small set of simultaneous equations
c
         nnn = no
         do 170 nn = 1 , npi
            if (.not.oconv(nn)) then
               do 150 nuu = 1 , no
                  do 140 noo = 1 , no
                     alpha(noo,nuu) = uau(noo,nuu,nn)
 140              continue
 150           continue
               do 160 noo = 1 , no
                  alpha(noo,noo) = alpha(noo,noo) + uu(noo,nn)
 160           continue
               ifail = 0
               call f04atf(alpha,iaa,b(1,nn),nnn,cc(1,nn),aa,iaa,
     +                     wks1,wks2,ifail)
            endif
 170     continue
c
c      form new solution vectors, overwriting au
c
         call vclr(au,1,mn*npi)
         do 200 j = 1 , no
            call rdedx(u,mn*npi,iblu(j),ifils)
            do 190 nn = 1 , npi
               ccjn = cc(j,nn)
               if (dabs(ccjn).gt.tich.and..not.oconv(nn)) then
                  do 180 i = 1 , mn
                     au(i,nn) = au(i,nn) + ccjn*u(i,nn)
 180              continue
               end if
 190        continue
 200     continue
c
c      carry any converged solutions forward
c      note that u is reused for the convergence check
c
         if (no.ne.1) then
            call rdedx(u,mn*npi,iblut(no-1),ifils)
            do nn=1,npi
               if (oconv(nn)) then
                  do i=1,mn
                     au(i,nn) = u(i,nn)
                  enddo
               endif
            enddo
         endif
c
c      write total solution onto file
c
         iblut(no) = iblast
         iblast = iblast + ibadd
         call wrt3(au,mn*npi,iblut(no),ifils)
c
c      check for convergence
c
         if (no.ne.1) then
            gmax = 0.0d0
            do 210 nn = 1 , npi
c
               if (.not.oconv(nn)) then
                  call vsub(au(1,nn),1,u(1,nn),1,u(1,nn),1,mn)
                  gnorm = ddot(mn,u(1,nn),1,u(1,nn),1)/dble(mn)
                  gnorm = dsqrt(gnorm)
                  gmax = dmax1(gmax,gnorm)
                  oconv(nn) = gnorm.le.uconv
               endif
 210        continue
            if (oprn(6)) write (iwr,6050) no , gmax
            if (gmax.le.uconv) then
               write (iwr,6090)
               go to 280
            end if
         end if
c
c     form new expansion vectors
c
         call rdedx(au,mn*npi,iblau(no),ifils)
         do 240 mo = 1 , no
            call rdedx(u,mn*npi,iblu(mo),ifils)
            do 230 nn = 1 , npi
               if (.not.oconv(nn)) then
                  fac = uau(mo,no,nn)/uu(mo,nn)
                  do 220 i = 1 , mn
                     au(i,nn) = au(i,nn) - fac*u(i,nn)
 220              continue
               endif
 230        continue
 240     continue
         gmax = 0.0d0
         do 260 nn = 1 , npi
            if (.not.oconv(nn)) then
               do 250 i = 1 , mn
                  if (dabs(au(i,nn)).le.tich) au(i,nn) = 0.0d0
 250           continue
               gnorm = ddot(mn,au(1,nn),1,au(1,nn),1)/dble(mn)
               gnorm = dsqrt(gnorm)
               gmax = dmax1(gmax,gnorm)
            endif
 260     continue
         if (oprn(6)) write (iwr,6060) gmax
         if (gmax.le.smal) then
           write (iwr,*) ' converged - new expansion vector negligible '
           go to 280
         end if
c
c
         iblu(no+1) = iblast
         iblast = iblast + ibadd
         call wrt3(au,mn*npi,iblu(no+1),ifils)
c
c
 270  continue
      write (iwr,*) ' no full convergence after 50 iterations '
      write (iwr,*) ' this will require changes to the program ! '
c
 280  call timit(3)
      write (iwr,6070) nox , cpulft(1) ,charwall()
      call rdedx(u,mn*npi,iblut(nox),ifils)
 290  npi = 0
      do 300 ko = npstar + 1 , npfin
         if (skipp(ko)) then
            call vclr(au,1,mn)
            call wrt3(au,mn,iblku,ifils)
         else
            npi = npi + 1
            call wrt3(u(1,npi),mn,iblku,ifils)
            if (oprn(13)) write (iwr,6080) npi , (u(i,npi),i=1,mn)
         end if
         iblku = iblku + lenblk
 300  continue
      return
 6000 format(/1x,
     +'commence iterative solution of chf equations at ',f8.2,
     +' seconds',a10,' wall')
 6010 format (//1x,'print right-hand-side to chf equations')
 6020 format (//1x,'print solutions to chf equations')
 6030 format (1x,'perturbation',i5,' omitted')
 6040 format (//1x,'perturbation  ',i4//(5x,5f16.8))
 6050 format (i10,5x,f15.10)
 6060 format (30x,f20.15)
 6070 format (/1x,
     + 'chf converged at iteration',i4/1x,
     + 'chf complete at ',f8.2,' seconds',a10,' wall')
 6080 format (//1x,'solution  ',i4//(5x,5f16.8))
 6090 format(/1x,'chf converged - wavefunctions stationary')
 6100 format(/
     +  6x,'iteration',9x,'tester',2x,'expansion vector norm'/
     +  6x,47('=')/)
      end
      function memreq_chfeqv(q,skipp,npx)
      implicit real*8  (a-h,o-z)
c
c     new version of chfeq - vector algorithm -
c
      logical skipp
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
      common/blkin/g(510),nword
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      common/small/alpha(50,50),aa(50,50),wks1(50),
     + wks2(50),iblu(50),iblut(50),iblau(50),iatms(100)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      dimension q(*)
c     dimension eps(mn),au(mn,npx),u(mn,npx),b(maxc,npx)
      dimension skipp(100)
c     dimension cc(maxc,npx),uu(maxc,npx)
c     dimension work(irmax,mn),uau(maxc,maxc,npx)
c
c     OCONV indicates whether a linear system has converged. 
c     Its dimension is taken to be equal to that of skipp.
c
c     logical oconv(100)
c     logical oconv(*)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      data smal/1.0d-13/
      data tich/1.0d-24/

      memreq_chfeqv = 0
      npi = 0
      do 20 ko = npstar + 1 , npfin
c
c      get rhs for this perturbation
c
         if (skipp(ko)) then
         else
            npi = npi + 1
            iatms(npi) = (ko-1)/3+1
         end if
 20   continue
c
c     add DFT contributions
c
         memreq_chfeqv = memreq_au_dft(q,q,q,npi,iatms)
      return
 6000 format(/1x,
     +'commence iterative solution of chf equations at ',f8.2,
     +' seconds',a10,' wall')
 6010 format (//1x,'print right-hand-side to chf equations')
 6020 format (//1x,'print solutions to chf equations')
 6030 format (1x,'perturbation',i5,' omitted')
 6040 format (//1x,'perturbation  ',i4//(5x,5f16.8))
 6050 format (i10,5x,f15.10)
 6060 format (30x,f20.15)
 6070 format (/1x,
     + 'chf converged at iteration',i4/1x,
     + 'chf complete at ',f8.2,' seconds',a10,' wall')
 6080 format (//1x,'solution  ',i4//(5x,5f16.8))
 6090 format(/1x,'chf converged - wavefunctions stationary')
 6100 format(/
     +  6x,'iteration',9x,'tester',2x,'expansion vector norm'/
     +  6x,47('=')/)
      end
      subroutine au_dft(q,u,au,npi,iatms)
      implicit none
c
c     Add the DFT contributions to the matrix vector products
c
c     Parameters:
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      integer m8
      parameter (m8=8)
c
c     Input:
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
c     COMMON/DRIVE_DFT/ contains options that specify how GAMESS-UK
c     should drive the CCP1 DFT module.
c
      integer KS_AO,     KS_MO,     KS_AOMO
      integer KS_GRD_AO, KS_GRD_MO, KS_GRD_AOMO
      integer KS_RHS_AO, KS_RHS_MO, KS_RHS_AOMO
      integer KS_LHS_AO, KS_LHS_MO, KS_LHS_AOMO
      integer KS_DT_AO,  KS_DT_MO,  KS_DT_AOMO
      integer KS_HES_AO, KS_HES_MO, KS_HES_AOMO
      integer KS_DX_AO,  KS_DX_MO,  KS_DX_AOMO
      parameter (KS_AO    =01,KS_MO    =02,KS_AOMO    =03)
      parameter (KS_GRD_AO=04,KS_GRD_MO=05,KS_GRD_AOMO=06)
      parameter (KS_RHS_AO=07,KS_RHS_MO=08,KS_RHS_AOMO=09)
      parameter (KS_LHS_AO=10,KS_LHS_MO=11,KS_LHS_AOMO=12)
      parameter (KS_DT_AO =13,KS_DT_MO =14,KS_DT_AOMO =15)
      parameter (KS_HES_AO=16,KS_HES_MO=17,KS_HES_AOMO=18)
      parameter (KS_DX_AO =19,KS_DX_MO =20,KS_DX_AOMO =21)
c
c     if ks_bas.eq.KS_AO           call CD_energy_ao
c     if ks_bas.eq.KS_MO           call CD_energy_mo
c     if ks_bas.eq.KS_AOMO         call CD_energy
c
c     if ks_grd_bas.eq.KS_GRD_AO   call CD_force_ao
c     if ks_grd_bas.eq.KS_GRD_MO   call CD_force_mo
c     if ks_grd_bas.eq.KS_GRD_AOMO call CD_force
c
c     if ks_rhs_bas.eq.KS_RHS_AO   call CD_chf_rhs_ao
c     if ks_rhs_bas.eq.KS_RHS_MO   call CD_chf_rhs_mo
c     if ks_rhs_bas.eq.KS_RHS_AOMO call CD_chf_rhs
c
c     if ks_lhs_bas.eq.KS_LHS_AO   call CD_chf_lhs_ao
c     if ks_lhs_bas.eq.KS_LHS_MO   call CD_chf_lhs_mo
c     if ks_lhs_bas.eq.KS_LHS_AOMO call CD_chf_lhs
c
c     if ks_dt_bas.eq.KS_DT_AO     call CD_dksm_ao
c     if ks_dt_bas.eq.KS_DT_MO     call CD_dksm_mo
c     if ks_dt_bas.eq.KS_DT_AOMO   call CD_dksm
c
c     if ks_hes_bas.eq.KS_HES_AO   call CD_hess_ao
c     if ks_hes_bas.eq.KS_HES_MO   call CD_hess_mo
c     if ks_hes_bas.eq.KS_HES_AOMO call CD_hess
c
c     if ks_dx_bas.eq.KS_DX_AO     call CD_dksm_exp_ao
c     if ks_dx_bas.eq.KS_DX_MO     call CD_dksm_exp_mo
c     if ks_dx_bas.eq.KS_DX_AOMO   call CD_dksm_exp
c
      integer ks_bas,     ks_grd_bas, ks_hes_bas
      integer ks_rhs_bas, ks_lhs_bas, ks_dt_bas
      integer ks_dx_bas
      logical ogeompert
c
      common/drive_dft/ks_bas, ks_grd_bas, 
     &       ks_rhs_bas, ks_lhs_bas, ks_dt_bas,
     &       ks_hes_bas, ks_dx_bas,
     &       ogeompert
      integer npi,iatms(npi)
      real*8 u(mn,npi)
      real*8 au(mn,npi)
c
c     Workspace:
c
      real*8 q(*)
c
c     Functions:
c
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
      integer igmem_alloc_inf
      integer igmem_null
c
c     Local:
c
      integer iavc,ierror,iblok,inull
c
c     Code:
c
      if (.not.CD_active()) return
c
c     get the MO coefficients
c
      iavc = igmem_alloc_inf(num*num,'cphf.m','au_dft','alpha-vectors',
     &                       IGMEM_NORMAL)
      call secget(isect(8),m8,iblok)
      iblok = iblok + mvadd
      call rdedx(q(iavc),num*ncoorb,iblok,ifild)
c
      inull = igmem_null()
      ierror = CD_chf_lhs_mo(q,q,npi,iatms,ogeompert,ncoorb,nocca,0,
     &                       q(iavc),q(inull),u,q(inull),au,q(inull),
     &                       .false.,iwr)
c
      call gmem_free_inf(iavc,'cphf.m','au_dft','alpha-vectors')
      end
      integer function memreq_au_dft(q,u,au,npi,iatms)
      implicit none
c
c     Add the DFT contributions to the matrix vector products
c
c     Parameters:
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      integer m8
      parameter (m8=8)
c
c     Input:
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
c     COMMON/DRIVE_DFT/ contains options that specify how GAMESS-UK
c     should drive the CCP1 DFT module.
c
      integer KS_AO,     KS_MO,     KS_AOMO
      integer KS_GRD_AO, KS_GRD_MO, KS_GRD_AOMO
      integer KS_RHS_AO, KS_RHS_MO, KS_RHS_AOMO
      integer KS_LHS_AO, KS_LHS_MO, KS_LHS_AOMO
      integer KS_DT_AO,  KS_DT_MO,  KS_DT_AOMO
      integer KS_HES_AO, KS_HES_MO, KS_HES_AOMO
      integer KS_DX_AO,  KS_DX_MO,  KS_DX_AOMO
      parameter (KS_AO    =01,KS_MO    =02,KS_AOMO    =03)
      parameter (KS_GRD_AO=04,KS_GRD_MO=05,KS_GRD_AOMO=06)
      parameter (KS_RHS_AO=07,KS_RHS_MO=08,KS_RHS_AOMO=09)
      parameter (KS_LHS_AO=10,KS_LHS_MO=11,KS_LHS_AOMO=12)
      parameter (KS_DT_AO =13,KS_DT_MO =14,KS_DT_AOMO =15)
      parameter (KS_HES_AO=16,KS_HES_MO=17,KS_HES_AOMO=18)
      parameter (KS_DX_AO =19,KS_DX_MO =20,KS_DX_AOMO =21)
c
c     if ks_bas.eq.KS_AO           call CD_energy_ao
c     if ks_bas.eq.KS_MO           call CD_energy_mo
c     if ks_bas.eq.KS_AOMO         call CD_energy
c
c     if ks_grd_bas.eq.KS_GRD_AO   call CD_force_ao
c     if ks_grd_bas.eq.KS_GRD_MO   call CD_force_mo
c     if ks_grd_bas.eq.KS_GRD_AOMO call CD_force
c
c     if ks_rhs_bas.eq.KS_RHS_AO   call CD_chf_rhs_ao
c     if ks_rhs_bas.eq.KS_RHS_MO   call CD_chf_rhs_mo
c     if ks_rhs_bas.eq.KS_RHS_AOMO call CD_chf_rhs
c
c     if ks_lhs_bas.eq.KS_LHS_AO   call CD_chf_lhs_ao
c     if ks_lhs_bas.eq.KS_LHS_MO   call CD_chf_lhs_mo
c     if ks_lhs_bas.eq.KS_LHS_AOMO call CD_chf_lhs
c
c     if ks_dt_bas.eq.KS_DT_AO     call CD_dksm_ao
c     if ks_dt_bas.eq.KS_DT_MO     call CD_dksm_mo
c     if ks_dt_bas.eq.KS_DT_AOMO   call CD_dksm
c
c     if ks_hes_bas.eq.KS_HES_AO   call CD_hess_ao
c     if ks_hes_bas.eq.KS_HES_MO   call CD_hess_mo
c     if ks_hes_bas.eq.KS_HES_AOMO call CD_hess
c
c     if ks_dx_bas.eq.KS_DX_AO     call CD_dksm_exp_ao
c     if ks_dx_bas.eq.KS_DX_MO     call CD_dksm_exp_mo
c     if ks_dx_bas.eq.KS_DX_AOMO   call CD_dksm_exp
c
      integer ks_bas,     ks_grd_bas, ks_hes_bas
      integer ks_rhs_bas, ks_lhs_bas, ks_dt_bas
      integer ks_dx_bas
      logical ogeompert
c
      common/drive_dft/ks_bas, ks_grd_bas, 
     &       ks_rhs_bas, ks_lhs_bas, ks_dt_bas,
     &       ks_hes_bas, ks_dx_bas,
     &       ogeompert
      integer npi,iatms(npi)
      real*8 u(mn,npi)
      real*8 au(mn,npi)
c
c     Workspace:
c
      real*8 q(*)
c
c     Functions:
c
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
      integer igmem_alloc_inf, igmem_overhead
      integer igmem_null
c
c     Local:
c
      integer iavc,ierror,iblok,inull
c
c     Code:
c
      memreq_au_dft = 0
      if (.not.CD_active()) return
      inull = igmem_null()
      iavc = inull
      memreq_au_dft = num*num + igmem_overhead() +
     &  CD_memreq_chf_lhs_mo(q,q,npi,iatms,ogeompert,ncoorb,nocca,0,
     &                       q(iavc),q(inull),u,q(inull),au,q(inull),
     &                       .false.,iwr)
      end
      subroutine symrhs(qq,ibstar,skipp,mapnr,iso,nshels)
c
c    symmetrise r.h.s. ( d2h and subgroups only )
c
      implicit real*8  (a-h,o-z)
      logical skipp
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension skipp(3,nat),mapnr(*)
      dimension qq(*),iso(nshels,*)
c
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
      common/mpshl/ns(maxorb)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c IN_CLUDE(common/cphf)
c
      common/bufb/ptr(3,144),ict(maxat,8)
      common/symmos/imos(8,maxorb)
      character *8 grhf
      data one/1.0d0/
      data grhf/'grhf'/
c
c
c
      call rdedx(ptr(1,1),nw196(1),ibl196(1),ifild)
      nav = lenwrd()
      call readi(iso,nw196(5)*nav,ibl196(5),ifild)
      do 40 ii = 1 , nshell
         ic = katom(ii)
         do 30 it = 1 , nt
            id = iso(ii,it)
            ict(ic,it) = katom(id)
 30      continue
 40   continue
c
      do 80 n = 1 , nat
         do 50 i = 1 , 3
            skipp(i,n) = .false.
 50      continue
         do 70 nop = 1 , nt
            if (ict(n,nop).gt.n) then
               do 60 i = 1 , 3
                  skipp(i,n) = .true.
 60            continue
               go to 80
            end if
 70      continue
 80   continue
      do 110 n = 1 , nat
         do 90 nop = 1 , nt
            if (ict(n,nop).ne.n) go to 110
 90      continue
         do 100 i = 1 , 3
            skipp(i,n) = .true.
 100     continue
c        nuniq = n
         go to 120
 110  continue
 120  continue
cjiri-cphf
c - ignore transl invariance for first
c     if(itransinvar.eq.0) then
c     skipp(1,1)=.false.
c     skipp(2,1)=.false.
c     skipp(3,1)=.false.
c     endif
cjiriend
      ntpls1 = noccb + 1
      nplus1 = nocca + 1
      nsoc = noccb - nocca
      nvirta = nsa4 - noccb
      ibll = lensec(mn)
      iblku = ibstar
      an = one/dble(nt)
      ioff = mn + 1
      nat3 = nat*3
c     read in u vectors
      do 130 n = 1 , nat3
         call rdedx(qq(ioff),mn,iblku,ifils)
         iblku = iblku + ibll
         ioff = ioff + mn
 130  continue
c     loop over vectors
      iblku = ibstar
      do 370 n = 1 , nat
         do 360 nc = 1 , 3
            ioff = ((n-1)*3+nc)*mn
c     copy vector for atom n, component nc into work area
            do 140 i = 1 , mn
               qq(i) = qq(ioff+i)
 140        continue
c     work along the elements of this vector
c loop over double-single and double-virtual
            if (scftyp.eq.grhf) then
               ij = 0
               do 180 i = 1 , nsa4
                  do 170 ia = 1 , i
                     if (ns(i).ne.ns(ia)) then
                        ij = ij + 1
c     loop over symmetry operations
c     except identity
                        do 160 iop = 2 , nt
                           isign = imos(iop,i)*imos(iop,ia)
                           sign = dble(isign)
                           niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                           ioff = (niop-1)*3*mn
                           npnc = (iop-1)*3 + nc
                           do 150 k = 1 , 3
                              ioff = ioff + mn
                              qq(ij) = ptr(k,npnc)*sign*qq(ioff+ij)
     +                                 + qq(ij)
 150                       continue
 160                    continue
                     end if
 170              continue
 180           continue
               if (.not.(.not.lcpf .and. .not.cicv .and. (.not.cicx)))
     +             then
                  do 220 i = 1 , nsa4
                     do 210 ia = 1 , i
                        iia = i*(i-1)/2 + ia
                        if (mapnr(iia).lt.0) then
                           ij = mnnr - mapnr(iia)
                           do 200 iop = 2 , nt
                              isign = imos(iop,i)*imos(iop,ia)
                              sign = dble(isign)
                              niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                              ioff = (niop-1)*3*mn
                              npnc = (iop-1)*3 + nc
                              do 190 k = 1 , 3
                                 ioff = ioff + mn
                                 qq(ij) = ptr(k,npnc)*sign*qq(ioff+ij)
     +                              + qq(ij)
 190                          continue
 200                       continue
                        end if
 210                 continue
 220              continue
               end if
            else
               if (nocca.ne.0) then
                  do 260 i = 1 , nocca
                     do 250 ia = nplus1 , nsa4
                        ij = (ia-nocca-1)*nocca + i
c     loop over symmetry operations
c     except identity
                        do 240 iop = 2 , nt
                           isign = imos(iop,i)*imos(iop,ia)
                           sign = dble(isign)
                           niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                           ioff = (niop-1)*3*mn
                           npnc = (iop-1)*3 + nc
                           do 230 k = 1 , 3
                              ioff = ioff + mn
                              qq(ij) = ptr(k,npnc)*sign*qq(ioff+ij)
     +                                 + qq(ij)
 230                       continue
 240                    continue
 250                 continue
 260              continue
                  if (.not.(.not.lcpf .and. .not.cicv .and. (.not.cicx))
     +                ) then
                     do 300 i = 1 , nsa4
                        do 290 ia = 1 , i
                           iia = i*(i-1)/2 + ia
                           if (mapnr(iia).lt.0) then
                              ij = mnnr - mapnr(iia)
                              do 280 iop = 2 , nt
                                 isign = imos(iop,i)*imos(iop,ia)
                                 sign = dble(isign)
                                 niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                                 ioff = (niop-1)*3*mn
                                 npnc = (iop-1)*3 + nc
                                 do 270 k = 1 , 3
                                    ioff = ioff + mn
                                    qq(ij) = ptr(k,npnc)
     +                                 *sign*qq(ioff+ij) + qq(ij)
 270                             continue
 280                          continue
                           end if
 290                    continue
 300                 continue
                  end if
               end if
               if (noccb.ne.nocca) then
c open shell only - loop over single-virtual
                  do 340 i = nplus1 , noccb
                     do 330 ia = ntpls1 , nsa4
                        ij = nvirta*nocca + (ia-nsoc-1)*nsoc + i - nocca
c     loop over symmetry operations
c     except identity
                        do 320 iop = 2 , nt
                           isign = imos(iop,i)*imos(iop,ia)
                           sign = dble(isign)
                           niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                           ioff = (niop-1)*3*mn
                           npnc = (iop-1)*3 + nc
                           do 310 k = 1 , 3
                              ioff = ioff + mn
                              qq(ij) = ptr(k,npnc)*sign*qq(ioff+ij)
     +                                 + qq(ij)
 310                       continue
 320                    continue
 330                 continue
 340              continue
               end if
            end if
            do 350 i = 1 , mn
               qq(i) = an*qq(i)
 350        continue
            call wrt3(qq(1),mn,iblku,ifils)
            iblku = iblku + ibll
 360     continue
 370  continue
      return
      end
      subroutine chfcla(qq,iqq)
c
c    assemble coupled hartree fock contribution to closed shell
c    scf second derivatives
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      logical out
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      common/maxlen/maxq
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/bufb/ioffn(maxat*3),icol(maxorb)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      dimension qq(*),iqq(*)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
c     contibution from first derivative of density matrix
c     with first derivative of integrals
c
      nat3 = nat*3
      nlen = nat3*nat3
      iof = lenrel(nw196(5))
      iofs = nw196(5) + lenint(nat*nt)
      i10 = iofs + 1
      i11 = i10 + nlen
      maxw = maxq - nlen - nx - iofs
      if (maxw.lt.(nx*3)) call caserr(' insufficient core')
      do 20 i = 1 , nat3
         icol(i) = (i-1)*nat3 + iofs
 20   continue
      call vclr(qq(i10),1,nlen)
      ltri = ikyp(ncoorb)
      lenblk = lensec(ltri)
      out = odebug(6)
c     density matrix derivatives section 15
      iposd = iochf(15)
      npass = 1
      maxnuc = 0
      nadd = nat
      ntot = nx*nat3
 30   if (ntot.le.maxw) then
         do 50 ipass = 1 , npass
            minnuc = maxnuc + 1
            maxnuc = maxnuc + nadd
            if (maxnuc.gt.nat) maxnuc = nat
            nuc = maxnuc - minnuc + 1
            nuc3 = nuc*3
            k = (minnuc-1)*3 + 1
            ioffn(k) = nlen + nx + i10
            do 40 i = 1 , nuc3
               ioffn(k+1) = ioffn(k) + nx
               k = k + 1
 40         continue
 50      continue
         maxnuc = 0
         do 90 ipass = 1 , npass
c
c     derivatives of fock matrix (m.o. basis, no wavefunction term)
c     at section 13
            iposf = iochf(13)
            minnuc = maxnuc + 1
            maxnuc = maxnuc + nadd
            if (maxnuc.gt.nat) maxnuc = nat
            nuc = maxnuc - minnuc + 1
            nuc3 = nuc*3
            do 60 k = 1 , nuc3
               ioff = ioffn(k)
c
c     read perturbed density matrices
c
               call rdedx(qq(ioff),ltri,iposd,ifockf)
               iposd = iposd + lenblk
 60         continue
            do 80 nn = 1 , nat3
c
c     read derivative fock matrix
c
               call rdedx(qq(i11),ltri,iposf,ifockf)
               iposf = iposf + lenblk
               k = (minnuc-1)*3
               do 70 kk = 1 , nuc3
                  ioff = ioffn(k+kk)
c
c     contribution to second derivative
c
                  dum = tracep(qq(i11),qq(ioff),ncoorb)
                  qq(k+kk+icol(nn)) = dum
 70            continue
 80         continue
 90      continue
         call rdedx(qq(1),nw196(5),ibl196(5),ifild)
         if (out) then
            call dr2sym(qq(i10),qq(i11),iqq(1),iqq(iof+1),nat,nat3,
     +      nshell)
            write (iwr,6010)
            call prnder(qq(i10),nat3,iwr)
         end if
         call secget(isect(60),60,isec46)
         call rdedx(qq(i11),nlen,isec46,ifild)
         call vadd(qq(i11),1,qq(i10),1,qq(i11),1,nlen)
         call dr2sym(qq(i11),qq(i10),iqq(1),iqq(iof+1),nat,nat3,
     +               nshell)
         call wrt3(qq(i11),nlen,isec46,ifild)
         if (out) then
            write (iwr,6020)
            call prnder(qq(i11),nat3,iwr)
         end if
c        return
      else
         npass = npass + 1
         nadd = nat/npass + 1
         ntot = nadd*3*nx
         go to 30
      end if
 6010 format (//' coupled hartree-fock contribution')
 6020 format (//' total so far')
      end
      subroutine chfopa(qq,iqq)
c
c    assemble coupled hartree fock contribution to open shell
c    ( high-spin) scf second derivatives
c
      implicit real*8  (a-h,o-z)
      logical out
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      common/maxlen/maxq
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/bufb/ioffn(maxat*3),icol(maxorb)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      dimension qq(*),iqq(*)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
c     contibution from first derivative of density matrix
c     with first derivative of integrals
c
      nat3 = nat*3
      nlen = nat3*nat3
      iof = lenrel(nw196(5))
      iofs = nw196(5) + lenint(nat*nt)
      i10 = iofs + 1
      i11 = i10 + nlen
      i12 = i11 + nx
      maxw = maxq - nlen - nx - nx - iofs
      if (maxw.lt.(nx*6)) call caserr(' insufficient core')
      do 20 i = 1 , nat3
         icol(i) = (i-1)*nat3 + i10 - 1
 20   continue
      call vclr(qq(i10),1,nlen)
      ltri = ikyp(ncoorb)
      lenblk = lensec(ltri)
      out = odebug(6)
c     density matrix derivatives section 15
      iposd = iochf(15)
      npass = 1
      maxnuc = 0
      nadd = nat
      ntott = nx*nat3*2
 30   if (ntott.le.maxw) then
         do 50 ipass = 1 , npass
            minnuc = maxnuc + 1
            maxnuc = maxnuc + nadd
            if (maxnuc.gt.nat) maxnuc = nat
            nuc = maxnuc - minnuc + 1
            nuc3 = nuc*3
            k = (minnuc-1)*3 + 1
            ioffn(k) = nlen + nx + nx + i10
            do 40 i = 1 , nuc3
               ioffn(k+1) = ioffn(k) + nx
               k = k + 1
 40         continue
 50      continue
         maxnuc = 0
         do 110 ipass = 1 , npass
c
c     derivatives of fock matrix (m.o. basis, no wavefunction term)
c     at section 13. followed by derivatives of k matrix,(ka).
            iposf = iochf(13)
            iposk = iposf + 3*nat*lenblk
c  derivatives of overlap matrices at section 14
            iposs = iochf(14)
c  derivatives of density matrices at section 15
            iposd = iochf(15)
            minnuc = maxnuc + 1
            maxnuc = maxnuc + nadd
            if (maxnuc.gt.nat) maxnuc = nat
            nuc = maxnuc - minnuc + 1
            nuc3 = nuc*3
            nword3 = nuc3*nx
            do 60 k = 1 , nuc3
               ioff = ioffn(k)
c     read perturbed density matrix
               call rdedx(qq(ioff),ltri,iposd,ifockf)
c     read derivative overlap matrix
               call rdedx(qq(ioff+nword3),ltri,iposs,ifockf)
               iposd = iposd + lenblk
               iposs = iposs + lenblk
 60         continue
            do 100 nn = 1 , nat3
c     read derivative fock matrix
               call rdedx(qq(i11),ltri,iposf,ifockf)
c     read derivative k matrix
               call rdedx(qq(i12),ltri,iposk,ifockf)
               do 80 i = 1 , ncoorb
                  foci = 0.0d0
                  if (i.le.nb) foci = 1.0d0
                  if (i.le.na) foci = 2.0d0
                  do 70 j = 1 , i
                     focj = 0.0d0
                     if (j.le.nb) focj = 1.0d0
                     if (j.le.na) focj = 2.0d0
                     ijlen = i*(i-1)/2 + j - 1
                     qq(i11+ijlen) = qq(i11+ijlen) - (2.0d0-foci-focj)
     +                               *qq(i12+ijlen)
                     qq(i12+ijlen) = foci*focj*qq(i12+ijlen)
 70               continue
 80            continue
               iposf = iposf + lenblk
               iposk = iposk + lenblk
               k = (minnuc-1)*3
               do 90 kk = 1 , nuc3
                  ioff = ioffn(k+kk)
                  dum1 = tracep(qq(i11),qq(ioff),ncoorb)
                  dum2 = tracep(qq(i12),qq(ioff+nword3),nb)
                  qq(k+kk+icol(nn)) = dum1 + dum2
 90            continue
 100        continue
 110     continue
         call rdedx(qq(1),nw196(5),ibl196(5),ifild)
         if (out) then
            call dr2sym(qq(i10),qq(i11),iqq(1),iqq(iof+1),nat,nat3,
     +                  nshell)
            write (iwr,6010)
            call prnder(qq(i10),nat3,iwr)
         end if
         call secget(isect(60),60,isec46)
         call rdedx(qq(i11),nlen,isec46,ifild)
         call vadd(qq(i11),1,qq(i10),1,qq(i11),1,nlen)
         call dr2sym(qq(i11),qq(i10),iqq(1),iqq(iof+1),nat,nat3,
     +               nshell)
         call wrt3(qq(i11),nlen,isec46,ifild)
         if (out) then
            write (iwr,6020)
            call prnder(qq(i11),nat3,iwr)
         end if
         return
      else
         npass = npass + 1
         nadd = nat/npass + 1
         ntott = nadd*6*nx
         go to 30
      end if
 6010 format (//' coupled hartree-fock contribution')
 6020 format (//' total so far')
      end
      subroutine chfndr(q,iq)
c
c    driving routine for nuclear displacement chf routines
c    -----------------------------------------------------
c
      implicit real*8  (a-h,o-z)
      logical lstop,skipp,acore
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/lsort/skipp(3*maxat)
      common/maxlen/maxq
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
      logical lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis
      common/scfblk/en,etot,ehf,sh1(2),sh2(2),gap1(2),gap2(2),
     1              d12,d13,d23,canna,cannb,cannc,fx,fy,fz,
     2              lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis,
     3              ncyc,ischm,lock,maxit,nconv,npunch,lokcyc
      common/mpshl/ns(maxorb)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)
c
c     COMMON/DRIVE_DFT/ contains options that specify how GAMESS-UK
c     should drive the CCP1 DFT module.
c
      integer KS_AO,     KS_MO,     KS_AOMO
      integer KS_GRD_AO, KS_GRD_MO, KS_GRD_AOMO
      integer KS_RHS_AO, KS_RHS_MO, KS_RHS_AOMO
      integer KS_LHS_AO, KS_LHS_MO, KS_LHS_AOMO
      integer KS_DT_AO,  KS_DT_MO,  KS_DT_AOMO
      integer KS_HES_AO, KS_HES_MO, KS_HES_AOMO
      integer KS_DX_AO,  KS_DX_MO,  KS_DX_AOMO
      parameter (KS_AO    =01,KS_MO    =02,KS_AOMO    =03)
      parameter (KS_GRD_AO=04,KS_GRD_MO=05,KS_GRD_AOMO=06)
      parameter (KS_RHS_AO=07,KS_RHS_MO=08,KS_RHS_AOMO=09)
      parameter (KS_LHS_AO=10,KS_LHS_MO=11,KS_LHS_AOMO=12)
      parameter (KS_DT_AO =13,KS_DT_MO =14,KS_DT_AOMO =15)
      parameter (KS_HES_AO=16,KS_HES_MO=17,KS_HES_AOMO=18)
      parameter (KS_DX_AO =19,KS_DX_MO =20,KS_DX_AOMO =21)
c
c     if ks_bas.eq.KS_AO           call CD_energy_ao
c     if ks_bas.eq.KS_MO           call CD_energy_mo
c     if ks_bas.eq.KS_AOMO         call CD_energy
c
c     if ks_grd_bas.eq.KS_GRD_AO   call CD_force_ao
c     if ks_grd_bas.eq.KS_GRD_MO   call CD_force_mo
c     if ks_grd_bas.eq.KS_GRD_AOMO call CD_force
c
c     if ks_rhs_bas.eq.KS_RHS_AO   call CD_chf_rhs_ao
c     if ks_rhs_bas.eq.KS_RHS_MO   call CD_chf_rhs_mo
c     if ks_rhs_bas.eq.KS_RHS_AOMO call CD_chf_rhs
c
c     if ks_lhs_bas.eq.KS_LHS_AO   call CD_chf_lhs_ao
c     if ks_lhs_bas.eq.KS_LHS_MO   call CD_chf_lhs_mo
c     if ks_lhs_bas.eq.KS_LHS_AOMO call CD_chf_lhs
c
c     if ks_dt_bas.eq.KS_DT_AO     call CD_dksm_ao
c     if ks_dt_bas.eq.KS_DT_MO     call CD_dksm_mo
c     if ks_dt_bas.eq.KS_DT_AOMO   call CD_dksm
c
c     if ks_hes_bas.eq.KS_HES_AO   call CD_hess_ao
c     if ks_hes_bas.eq.KS_HES_MO   call CD_hess_mo
c     if ks_hes_bas.eq.KS_HES_AOMO call CD_hess
c
c     if ks_dx_bas.eq.KS_DX_AO     call CD_dksm_exp_ao
c     if ks_dx_bas.eq.KS_DX_MO     call CD_dksm_exp_mo
c     if ks_dx_bas.eq.KS_DX_AOMO   call CD_dksm_exp
c
      integer ks_bas,     ks_grd_bas, ks_hes_bas
      integer ks_rhs_bas, ks_lhs_bas, ks_dt_bas
      integer ks_dx_bas
      logical ogeompert
c
      common/drive_dft/ks_bas, ks_grd_bas, 
     &       ks_rhs_bas, ks_lhs_bas, ks_dt_bas,
     &       ks_hes_bas, ks_dx_bas,
     &       ogeompert
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
      logical ogeompert_save
      dimension iq(*), q(*)
      character *8 grhf,oscf,closed
      data grhf/'grhf'/
      data oscf/'oscf'/
      data closed/'closed'/
c
      ierror = CD_set_2e()
      nav = lenwrd()
      lenx = lensec(nx)*nat*15
      call cpuwal(begin,ebegin)
      call wrt3z(iblks,ifils,lenx)
      write (iwr,6010)
      if (scftyp.eq.oscf) then
         if((dabs(canna-2.0d0).gt.1.0d-6).or.
     +      (dabs(cannb).gt.1.0d-6).or.
     +      (dabs(cannc+2.0d0).gt.1.0d-6)) then
            write (iwr,6020)
            call caserr('stop')
         end if
      end if
      np = nat*3
      if (scftyp.eq.grhf) then
c
         mtype = 0
         call secget(isect(53),mtype,iblok)
         call readi(nact,lds(isect(53))*nav,iblok,ifild)
c
         call derlag(q,erga,ergb,fjk)
         call bfnshl(ns,nsa4)
         imap = igmem_alloc(lenint(nx))
         iimap = lenrel(imap-1)+1
         call ijmapr(nsa4,ns,iq(iimap),mn,mnx)
         ieps = igmem_alloc(mn)
c
         lenblk = lensec(mn)
         ibeta = iblks + lenblk*nat*6
         ibzeta = ibeta + lensec(nx)
c
         i01 = igmem_alloc(nx)
         i1  = igmem_alloc(nx*njk1)
         i2  = igmem_alloc_all(maxa)
         call lgrhfm(q(i01),ibeta,erga,ergb,fjk)
         call lgrhf(q(i1),ibzeta,erga,ergb,fjk)
         call chfgrs(q(ieps),iq(iimap),q(i01),q(i1),q(i2),
     +               ibeta,ibzeta,maxa)
         call gmem_free(i2)
         call gmem_free(i1)
         call gmem_free(i01)
c
         i01 = igmem_alloc(mn*nat)
         i1  = igmem_alloc(mn*nat)
         nxt3 = nx*(5+njk1) + mn*3
         nxt1 = mn*nat+3*max(num*num,nx*nat)+nx*(njk1+1)
         maxa = igmem_max_memory()
         if (nxt1.le.maxa) then
            i2  = igmem_alloc(mn*nat)
            i3  = igmem_alloc(max(num*num,nx*nat))
            i4  = igmem_alloc(max(num*num,nx*nat))
            i5  = igmem_alloc(max(num*num,nx*nat))
            i6  = igmem_alloc(nx)
            i7  = igmem_alloc(nx*njk1)
            i8  = i3
            i9  = i4
            i10 = i5
            acore = .true.
         else
            call caserr("What a mess! Out store rhsgvb not coded???")
            if (nxt3.gt.maxq) then
               write (iwr,6030)
               call caserr('stop')
            end if
            i6 = i2 + mn
            i7 = i6 + nx
            i3 = i7 + nx*njk1
            i4 = i3 + num*num
            i5 = i4 + num*num
            i8 = i5 + num*num
            i9 = i8 + nx
            i10 = i9 + nx
            acore = .false.
         end if
         call rhsgvb(iq(iimap),q(i01),q(i1),q(i2),q(i3),q(i4),
     +        q(i5),q(i6),q(i7),q(i8),q(i9),q(i10),acore,
     +        erga,ergb,ibeta,ibzeta)
         if (acore) then
            call gmem_free(i7)
            call gmem_free(i6)
            call gmem_free(i5)
            call gmem_free(i4)
            call gmem_free(i3)
            call gmem_free(i2)
         else
         endif
         call gmem_free(i1)
         call gmem_free(i01)
c
         i01 = igmem_alloc(mn*(3*nat+1))
         iso = igmem_alloc(nw196(5))
         iiso = lenrel(iso-1)+1
         call symrhs(q(i01),iblks,skipp,iq(iimap),iq(iiso),nshell)
         call gmem_free(iso)
         call gmem_free(i01)
         lenblk = lensec(mn)
         iblku = iblks + np*lenblk
         lstop = .false.
c
c        solve linear equations
c
         ogeompert_save = ogeompert
         ogeompert = .true.
         call chfdrv(q(ieps),lstop,skipp)
         ogeompert = ogeompert_save
         call gmem_free(ieps)
         call gmem_free(imap)
      else
         mn = noccb*nvirta + (noccb-nocca)*nocca
         call grhfbl(scftyp)
         call bfnshl(ns,nsa4)
c
c      sort out a-matrix
c
         i01 = igmem_alloc_all(maxa)
         if (scftyp.eq.closed) call chfcls(q(i01),maxa)
         if (scftyp.eq.oscf) call chfops(q(i01),maxa)
         call gmem_free(i01)
c
         mnmx = mn
         ieps = igmem_alloc(mnmx)
         ibx  = igmem_alloc(mnmx)
         iby  = igmem_alloc(mnmx)
         ibz  = igmem_alloc(mnmx)
         ieval= igmem_alloc(num)
c
c     r.h.s of equations
c
         if (scftyp.eq.closed) then
            call start_time_period(TP_2D_CHFRHS)
            isx  = igmem_alloc(nx)
            isy  = igmem_alloc(nx)
            isz  = igmem_alloc(nx)
            call rhscl(q(ieps),q(ibx),q(iby),q(ibz),q(ieval),
     +                 q(isx),q(isy),q(isz))
            call gmem_free(isz)
            call gmem_free(isy)
            call gmem_free(isx)
            call rhscl_dft(q,iq)
            call end_time_period(TP_2D_CHFRHS)
         endif
         if (scftyp.eq.oscf) then
            ibase = igmem_alloc_all(maxa)
            call rhsrhf(q(ieps),q(ibx),q(iby),q(ibz),q(ieval),q(ibase))
            call gmem_free(ibase)
         endif
         call gmem_free(ieval)
         call gmem_free(ibz)
         call gmem_free(iby)
         call gmem_free(ibx)
         iso  = igmem_alloc(nw196(5))
         iiso = lenrel(iso-1)+1
         imap = igmem_alloc(lenint(nx))
         iimap= lenrel(imap-1)+1
         mnmx = mn
         iwrk = igmem_alloc(mnmx*(3*nat+1))
         call start_time_period(TP_2D_SYMMRHS)
         call symrhs(q(iwrk),iblks,skipp,iq(iimap),iq(iiso),nshell)
         call end_time_period(TP_2D_SYMMRHS)
         call gmem_free(iwrk)
         call gmem_free(imap)
         call gmem_free(iso)
c
         lenblk = lensec(mn)
         iblku = iblks + np*lenblk
         lstop = .false.
c
c        solve linear equations
c
         ogeompert_save = ogeompert
         ogeompert = .true.
         call start_time_period(TP_2D_CHFDRV)
         call chfdrv(q(ieps),lstop,skipp)
         call end_time_period(TP_2D_CHFDRV)
         ogeompert = ogeompert_save
         call gmem_free(ieps)
      end if
c
      if (lstop) then
         call revise
         write (iwr,6040)
         call timana(22)
         call clenms('stop')
      else
         iso  = igmem_alloc(nw196(5))
         iiso = lenrel(iso-1)+1
         iwrk = igmem_alloc(3*mn*(nat+1))
         call start_time_period(TP_2D_SYMMU)
         call symmu(q(iwrk),iblku,skipp,iq(iiso),nshell)
         call end_time_period(TP_2D_SYMMU)
         call gmem_free(iwrk)
         call gmem_free(iso)
c
c     perturbed density matrices
c
         iwrk = igmem_alloc(2*nx+mn)
         call start_time_period(TP_2D_PDENS)
         call pdens(q(iwrk),lstop)
         call end_time_period(TP_2D_PDENS)
         call gmem_free(iwrk)
         call start_time_period(TP_2D_PFOCK)
         if (scftyp.eq.closed) then
            call pfockc(q)
            call pdksmc(q,iq)
         endif
         if (scftyp.eq.oscf) call pfocko(q)
         call end_time_period(TP_2D_PFOCK)
         call revise
         call clredx
         call delfil(nofile(1))
         call timana(22)
      end if
      ierror = CD_reset_2e()
 6010 format (/' solve chf equations for nuclear motions')
 6020 format (/1x,'oscf chf equations only work with',
     +        ' canonicalisation 2.0, 0.0, -2.0'/1x,
     +        '------     see manual for details --------')
 6030 format (//' insufficient store for gcphf ( gderci ) ')
 6040 format (//1x,'insufficient time to finish chf equations'//1x,
     +        'restart job'//)
      end
      subroutine pdksmc(q,iq)
      implicit none
c
c     Add the DFT wavefunction contribution to the perturbed Fock
c     matrix
c
c     Parameters:
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer m8
      parameter (m8=8)
c
c     Input:
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      integer mp2gcycl
      logical opg_grad, opg_alloc, opg_alloc_sch, opg_alloc_dmat
      logical opg_grad_com
      logical opg_pertbns_sep
      common/mp2grd/mp2gcycl, opg_grad, opg_alloc, opg_alloc_sch, 
     +              opg_alloc_dmat, opg_grad_com,
     +              opg_pertbns_sep
c
c
c     Workspace:
c
      integer iq(*)
      real*8 q(*)
c
c     Local:
c
      integer ibt, ibd, ltri, iatms,iiatms, iblok, inull
      integer nat3, i, j, k, l, ierror, n, m, lennew
      integer ida, ifa, icc, ida_t
      integer ij, ija
      integer ntot, nnow
      character *6 fnm
      character *6 snm
c
c     Functions:
c
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      integer lensec
      integer igmem_null
      integer igmem_alloc_inf, lenrel
      data fnm/'cphf.m'/
      data snm/'pdksmc'/
c
c     Code:
c
      if (.not.CD_active()) return
      inull = igmem_null()
c
      if (opg_pertbns_sep) then
c
c        Do everything in batches of 1 coordinate at a time
c        if this code works we should always use it just with 
c        different batch sizes...
c
         ltri = ikyp(ncoorb)
         nat3 = 3*nat
         ntot = 1
         ida = igmem_alloc_inf(ltri*ntot,fnm,snm,
     &                         'pert-density',IGMEM_NORMAL)
         ifa = igmem_alloc_inf(ltri*ntot,fnm,snm,
     &                         'pert-fock',IGMEM_NORMAL)
         icc = igmem_alloc_inf(num*num,fnm,snm,
     &                         'alpha-vectors',IGMEM_NORMAL)
         lennew = lensec(ltri)
         iatms = igmem_alloc_inf(nat3,fnm,snm,
     &                           'perturbations',IGMEM_DEBUG)
         iiatms = lenrel(iatms-1)+1
c
         call secget(isect(8),m8,iblok)
         iblok = iblok + mvadd
         call rdedx(q(icc),num*ncoorb,iblok,ifild)
c
         nnow  = 0
         n     = 0
         m     = 0
         ibd   = iochf(15)
         ibt   = iochf(16)
         ida_t = ida
         do i = 1, nat3
c
c           collect upto ntot perturbed density matrices
c
            n    = n + 1
            nnow = nnow + 1
            iq(iiatms+nnow-1) = (i-1)/3+1
            call rdedx(q(ida_t),ltri,ibd,ifockf)
c
            call dscal(ltri,0.5d0,q(ida_t),1)
            do j = 1, nocca
               q(ida_t-1+j*(j+1)/2) = 0.5d0*q(ida_t-1+j*(j+1)/2)
            enddo
c
            ida_t = ida_t + ltri
            ibd   = ibd   + lennew
c
            if (nnow.eq.ntot.or.i.eq.nat3) then
c
c              got enough perturbed density matrices so do the work
c
               call vclr(q(ifa),1,ltri*nnow)
               ierror = CD_chf_dksm_mo(iq,q,nnow,q(iatms),ncoorb,
     &                  nocca,0,q(icc),q(inull),q(ida),
     &                  q(inull),q(ifa),q(inull),.false.,iwr)
               do j = 1, nnow
                  m = m + 1
                  call rdedx(q(ida),ltri,ibt,ifockf)
                  ija = 0
                  do k = 1, nsa4
                     do l = 1, k
                        ij  = iky(mapie(k)) + mapie(l) - 1
                        q(ida+ij) = q(ida+ij) + q(ifa+(j-1)*ltri+ija)
                        ija = ija + 1
                     enddo
                  enddo
                  call wrt3(q(ida),ltri,ibt,ifockf)
                  ibt = ibt + lennew
               enddo
c
c              reset counters for the next batch
c
               nnow  = 0
               ida_t = ida
            endif
         enddo
         call gmem_free_inf(iatms,fnm,snm,'perturbations')
         call gmem_free_inf(icc,fnm,snm,'alpha-vectors')
         call gmem_free_inf(ifa,fnm,snm,'pert-fock')
         call gmem_free_inf(ida,fnm,snm,'pert-density')
         call revise
         call clredx
         return
      endif
c
      ltri = ikyp(ncoorb)
      nat3 = 3*nat
      ida = igmem_alloc_inf(ltri*nat3,fnm,snm,'pert-density',
     &                      IGMEM_NORMAL)
      ifa = igmem_alloc_inf(ltri*nat3,fnm,snm,'pert-fock',
     &                      IGMEM_NORMAL)
      icc = igmem_alloc_inf(num*num,fnm,snm,'alpha-vectors',
     &                      IGMEM_NORMAL)
      lennew = lensec(ltri)
      iatms = igmem_alloc_inf(nat3,fnm,snm,
     &                        'perturbations',IGMEM_DEBUG)
      iiatms = lenrel(iatms-1)+1
c
      call secget(isect(8),m8,iblok)
      iblok = iblok + mvadd
      call rdedx(q(icc),num*ncoorb,iblok,ifild)
c
      call vclr(q(ifa),1,ltri*nat3)
c
      ibd   = iochf(15)
      ida_t = ida
      n = 0
      do i = 1, nat
         do j = 1, 3
            n = n + 1
            iq(iiatms+n-1) = i
            call rdedx(q(ida_t),ltri,ibd,ifockf)
            ida_t = ida_t + ltri
            ibd   = ibd   + lennew
         enddo
      enddo
      call dscal(ltri*nat3,0.5d0,q(ida),1)
      do i = 1, nat3
         do j = 1, nocca
            q(ida-1+(i-1)*ltri+j*(j+1)/2)
     &      = 0.5d0*q(ida-1+(i-1)*ltri+j*(j+1)/2)
         enddo
      enddo
c
      ierror = CD_chf_dksm_mo(iq,q,nat3,q(iatms),ncoorb,nocca,0,
     &                        q(icc),q(inull),q(ida),q(inull),
     &                        q(ifa),q(inull),.false.,iwr)
c
      ibt = iochf(16)
      do k = 1, nat3
         call rdedx(q(icc),ltri,ibt,ifockf)
         ija = 0
         do i = 1, nsa4
            do j = 1, i
               ij  = iky(mapie(i)) + mapie(j) - 1
               q(icc+ij) = q(icc+ij) + q(ifa+(k-1)*ltri+ija)
               ija = ija + 1
            enddo
         enddo
         call wrt3(q(icc),ltri,ibt,ifockf)
         ibt = ibt + lennew
      enddo
c
      call gmem_free_inf(iatms,fnm,snm,'perturbations')
      call gmem_free_inf(icc,fnm,snm,'alpha-vectors')
      call gmem_free_inf(ifa,fnm,snm,'pert-fock')
      call gmem_free_inf(ida,fnm,snm,'pert-density')
      call revise
      call clredx
      end
      subroutine rhscl_dft(q,iq)
      implicit none
c
c     Adds the DFT contributions onto the right-hand-sides.
c
c     Parameters:
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer m8
      parameter(m8=8)
c
c     Input:
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c     COMMON/DRIVE_DFT/ contains options that specify how GAMESS-UK
c     should drive the CCP1 DFT module.
c
      integer KS_AO,     KS_MO,     KS_AOMO
      integer KS_GRD_AO, KS_GRD_MO, KS_GRD_AOMO
      integer KS_RHS_AO, KS_RHS_MO, KS_RHS_AOMO
      integer KS_LHS_AO, KS_LHS_MO, KS_LHS_AOMO
      integer KS_DT_AO,  KS_DT_MO,  KS_DT_AOMO
      integer KS_HES_AO, KS_HES_MO, KS_HES_AOMO
      integer KS_DX_AO,  KS_DX_MO,  KS_DX_AOMO
      parameter (KS_AO    =01,KS_MO    =02,KS_AOMO    =03)
      parameter (KS_GRD_AO=04,KS_GRD_MO=05,KS_GRD_AOMO=06)
      parameter (KS_RHS_AO=07,KS_RHS_MO=08,KS_RHS_AOMO=09)
      parameter (KS_LHS_AO=10,KS_LHS_MO=11,KS_LHS_AOMO=12)
      parameter (KS_DT_AO =13,KS_DT_MO =14,KS_DT_AOMO =15)
      parameter (KS_HES_AO=16,KS_HES_MO=17,KS_HES_AOMO=18)
      parameter (KS_DX_AO =19,KS_DX_MO =20,KS_DX_AOMO =21)
c
c     if ks_bas.eq.KS_AO           call CD_energy_ao
c     if ks_bas.eq.KS_MO           call CD_energy_mo
c     if ks_bas.eq.KS_AOMO         call CD_energy
c
c     if ks_grd_bas.eq.KS_GRD_AO   call CD_force_ao
c     if ks_grd_bas.eq.KS_GRD_MO   call CD_force_mo
c     if ks_grd_bas.eq.KS_GRD_AOMO call CD_force
c
c     if ks_rhs_bas.eq.KS_RHS_AO   call CD_chf_rhs_ao
c     if ks_rhs_bas.eq.KS_RHS_MO   call CD_chf_rhs_mo
c     if ks_rhs_bas.eq.KS_RHS_AOMO call CD_chf_rhs
c
c     if ks_lhs_bas.eq.KS_LHS_AO   call CD_chf_lhs_ao
c     if ks_lhs_bas.eq.KS_LHS_MO   call CD_chf_lhs_mo
c     if ks_lhs_bas.eq.KS_LHS_AOMO call CD_chf_lhs
c
c     if ks_dt_bas.eq.KS_DT_AO     call CD_dksm_ao
c     if ks_dt_bas.eq.KS_DT_MO     call CD_dksm_mo
c     if ks_dt_bas.eq.KS_DT_AOMO   call CD_dksm
c
c     if ks_hes_bas.eq.KS_HES_AO   call CD_hess_ao
c     if ks_hes_bas.eq.KS_HES_MO   call CD_hess_mo
c     if ks_hes_bas.eq.KS_HES_AOMO call CD_hess
c
c     if ks_dx_bas.eq.KS_DX_AO     call CD_dksm_exp_ao
c     if ks_dx_bas.eq.KS_DX_MO     call CD_dksm_exp_mo
c     if ks_dx_bas.eq.KS_DX_AOMO   call CD_dksm_exp
c
      integer ks_bas,     ks_grd_bas, ks_hes_bas
      integer ks_rhs_bas, ks_lhs_bas, ks_dt_bas
      integer ks_dx_bas
      logical ogeompert
c
      common/drive_dft/ks_bas, ks_grd_bas, 
     &       ks_rhs_bas, ks_lhs_bas, ks_dt_bas,
     &       ks_hes_bas, ks_dx_bas,
     &       ogeompert
c
c     Workspace:
c
      real*8 q(*)
      integer iq(*)
c
c     Functions:
c
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      integer igmem_alloc_inf, igmem_null
      integer lensec, lenrel
c
c     Local:
c
      integer nat3,lennew,newblk,mnblk,is,ib,ic,iblok,iatms,iiatms
      integer ibs,ics,icb,i,j,ierror,iblkb
      integer mu, nu, ni, nd, nk, num2
      integer id, iscr, itmp, inull
      integer m7
      parameter(m7=7)
      character*6 fnm
      character*9 snm
      data fnm/'cphf.m'/
      data snm/'rhscl_dft'/
      if (.not.CD_active()) return

      nat3   = nat*3
      lennew = iky(ncoorb)+ncoorb
      newblk = lensec(lennew)
      mnblk  = lensec(mn)
      inull  = igmem_null()
      if (ks_rhs_bas.eq.KS_RHS_MO) then
         is = igmem_alloc_inf(lennew*nat3,fnm,snm,
     &                        'pert-overlap',IGMEM_NORMAL)
         ib = igmem_alloc_inf(mn*nat3,fnm,snm,
     &                        'right-hand-sides',IGMEM_NORMAL)
         ic = igmem_alloc_inf(num*num,fnm,snm,
     &                        'alpha-vectors',IGMEM_NORMAL)
c
c        Get the MO-coefficients
c
         call secget(isect(8),m8,iblok)
         iblok = iblok + mvadd
         call rdedx(q(ic),num*ncoorb,iblok,ifild)
c
c        Get the derivative overlap matrices
c
         ibs = iochf(14)
         ics = is
         do i = 1, nat3
            call rdedx(q(ics),lennew,ibs,ifockf)
            ics = ics + lennew
            ibs = ibs + newblk
         enddo
c
c        Calculate the RHS contributions
c
         iatms = igmem_alloc_inf(nat3,fnm,snm,
     &                           'perturbations',IGMEM_DEBUG)
         iiatms = lenrel(iatms-1)+1
         do i = 1, nat
            iq(iiatms+3*(i-1)+0) = i
            iq(iiatms+3*(i-1)+1) = i
            iq(iiatms+3*(i-1)+2) = i
         enddo
         call vclr(q(ib),1,mn*nat3)
         ierror = CD_chf_rhs_mo(iq,q,nat3,iq(iiatms),ncoorb,nocca,0,
     &                          q(ic),q(inull),q(is),q(inull),
     &                          q(ib),q(inull),.false.,iwr)
         call gmem_free_inf(iatms,fnm,snm,'perturbations')
c
c        Add the DFT RHS contributions onto the Hartree-Fock parts
c
         iblkb = iblks
         icb = ib
         do i = 1, nat3
            call rdedx(q(ic),mn,iblkb,ifils)
            do j = 0, mn-1
               q(ic+j)=q(ic+j)+q(icb+j)
            enddo
            call wrt3(q(ic),mn,iblkb,ifils)
            iblkb = iblkb + mnblk
            icb   = icb   + mn
         enddo
         call clredx
         call gmem_free_inf(ic,fnm,snm,'alpha-vectors')
         call gmem_free_inf(ib,fnm,snm,'right-hand-sides')
         call gmem_free_inf(is,fnm,snm,'pert-overlap')
      else if (ks_rhs_bas.eq.KS_RHS_AO) then
         num2 = num*num
         ib = igmem_alloc_inf(num2*nat3,fnm,snm,'right-hand-side-ao',
     &                        IGMEM_NORMAL)
         is = igmem_alloc_inf(num2*nat3,fnm,snm,'pert-overlap-ao',
     &                        IGMEM_NORMAL)
         ic = igmem_alloc_inf(num2,fnm,snm,'vectors',
     &                        IGMEM_NORMAL)
c
c        get the mo-coefficients
c
         call secget(isect(8),m8,iblok)
         iblok = iblok + mvadd
         call rdedx(q(ic),num*ncoorb,iblok,ifild)
c
c        load the derivative overlap matrices and transform them
c        to AO basis.
c
         itmp = igmem_alloc_inf(lennew,fnm,snm,'temp',IGMEM_DEBUG)
         iscr = igmem_alloc_inf(num,fnm,snm,'scratch',IGMEM_DEBUG)
         call vclr(q(is),1,num*num*nat3)
         ibs = iochf(14)
         do i = 1, nat3
            call rdedx(q(itmp),lennew,ibs,ifockf)
            ibs = ibs + newblk
            do ni = 1, nocc
               call vclr(q(iscr),1,num)
               do nu = 1, num
                  do nk = 1, nocc
                     q(iscr-1+nu)=q(iscr-1+nu)
     &                           -0.5d0*q(ic-1+nu+(nk-1)*num)
     &                           *q(itmp-1+iky(max(ni,nk))+min(ni,nk))
                  enddo
               enddo
               do nu = 1, num
                  do mu = 1, num
                     q(is-1+mu+(nu-1)*num+(i-1)*num2) 
     &               = q(is-1+mu+(nu-1)*num+(i-1)*num2)
     &               + q(ic-1+nu+(ni-1)*num)*q(iscr-1+mu)
     &               + q(ic-1+mu+(ni-1)*num)*q(iscr-1+nu)
                  enddo
               enddo
            enddo
         enddo
         call gmem_free_inf(iscr,fnm,snm,'scratch')
         call gmem_free_inf(itmp,fnm,snm,'temp')
         call gmem_free_inf(ic,fnm,snm,'vectors')
c
c        load the density matrix
c
         id = igmem_alloc_inf(lennew,fnm,snm,'density',IGMEM_NORMAL)
         call secget(isect(7),m7,iblok)
         call rdedx(q(id),nx,iblok,ifild)
c
c        call DFT module
c
         iatms = igmem_alloc_inf(nat3,fnm,snm,'perturbations',
     &                           IGMEM_DEBUG)
         iiatms = lenrel(iatms-1)+1
         do i = 1, nat
            iq(iiatms+3*(i-1)+0) = i
            iq(iiatms+3*(i-1)+1) = i
            iq(iiatms+3*(i-1)+2) = i
         enddo
         call vclr(q(ib),1,num2*nat3)
         ierror = CD_chf_rhs_ao(iq,q,nat3,iq(iiatms),q(id),q(inull),
     &            q(is),q(inull),q(ib),q(inull),.false.,iwr)
         call gmem_free_inf(iatms,fnm,snm,'perturbations')
         call gmem_free_inf(id,fnm,snm,'density')
         call gmem_free_inf(is,fnm,snm,'pert-overlap-ao')
c
c        load the vectors again
c
         ic = igmem_alloc_inf(num*ncoorb,fnm,snm,'vectors',IGMEM_NORMAL)
         call secget(isect(8),m8,iblok)
         iblok = iblok + mvadd
         call rdedx(q(ic),num*ncoorb,iblok,ifild)
c
c        transform and store DFT results
c
         itmp = igmem_alloc_inf(mn,fnm,snm,'temp',IGMEM_DEBUG)
         iscr = igmem_alloc_inf(num,fnm,snm,'scratch',IGMEM_DEBUG)
         iblkb = iblks
         do i = 1, nat3
            call rdedx(q(itmp),mn,iblkb,ifils)
            do ni = 1, nocc
               call vclr(q(iscr),1,num)
               do mu = 1, num
                  do nu = 1, num
                     q(iscr-1+mu) = q(iscr-1+mu)
     &               + q(ic-1+nu+(ni-1)*num)*
     &                 q(ib-1+mu+(nu-1)*num+(i-1)*num2)
                  enddo
               enddo
               do nd = 1, nvirt
                  do mu = 1, num
                     q(itmp-1+ni+(nd-1)*nocc) = q(itmp-1+ni+(nd-1)*nocc)
     &               + q(ic-1+mu+(nocc+nd-1)*num)*q(iscr-1+mu)
                  enddo
               enddo
            enddo
            call wrt3(q(itmp),mn,iblkb,ifils)
            iblkb = iblkb + mnblk
         enddo
         call clredx
         call gmem_free_inf(iscr,fnm,snm,'scratch')
         call gmem_free_inf(itmp,fnm,snm,'temp')
         call gmem_free_inf(ic,fnm,snm,'vectors')
         call gmem_free_inf(ib,fnm,snm,'right-hand-side-ao')

      else if (ks_rhs_bas.eq.KS_RHS_AOMO) then
         call caserr('rhscl_dft: KS_RHS_AOMO implemented yet!')
      else
         call caserr('rhscl_dft: invalid option')
      endif
      end
      subroutine efdens(u,d)
c
c   constructs first order density matrices from electric-field
c   perturbed orbitals
c  ------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension u(*),d(*)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      common/mpshl/ns(maxorb)
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      character *8 pnames
      common/tdhfx/pnames(50)
c
c
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      character *8 grhf
      data grhf/'grhf'/
c
      length = lensec(mn)
      iblku = iblks + np*length
      nsoc = noccb - nocca
      ndp1 = nocca + 1
      nsp1 = noccb + 1
      if (odebug(3)) write (iwr,6010)
      nw = iky(ncoorb+1)
      iblll = lensec(nw)
      icomp = 31
      do 110 npert = 1 , 9
         if (.not.(opskip(npert))) then
            call rdedx(u,mn,iblku,ifils)
            iblku = iblku + length
            do 20 j = 1 , nw
               d(j) = 0.0d0
 20         continue
            if (scftyp.eq.grhf) then
               nr = 0
               do 40 i = 1 , nsa4
                  do 30 j = 1 , i
                     if (ns(i).ne.ns(j)) then
                        nr = nr + 1
                        ij = iky(mapie(i)) + mapie(j)
                        d(ij) = u(nr)*(fjk(ns(j))-fjk(ns(i)))
                     end if
 30               continue
 40            continue
            else
               if (nsoc.ne.0 .and. nocca.ne.0) then
                  do 60 i = ndp1 , noccb
                     do 50 j = 1 , nocca
                        ij = iky(mapie(i)) + mapie(j)
                        mt = (i-nocca-1)*nocca + j
                        d(ij) = u(mt)
 50                  continue
 60               continue
               end if
               if (nocca.ne.0) then
                  do 80 i = nsp1 , nsa4
                     do 70 j = 1 , nocca
                        ij = iky(mapie(i)) + mapie(j)
                        mt = (i-nocca-1)*nocca + j
                        d(ij) = u(mt) + u(mt)
 70                  continue
 80               continue
               end if
               if (nsoc.ne.0) then
                  do 100 i = nsp1 , nsa4
                     do 90 j = ndp1 , noccb
                        ij = iky(mapie(i)) + mapie(j)
                        mt = nocca*nvirta + (i-nsoc-1)*nsoc + j - nocca
                        d(ij) = u(mt)
 90                  continue
 100              continue
               end if
            end if
            call secput(isect(icomp),icomp,iblll,iblok)
            lds(isect(icomp)) = nw
            if (odebug(3)) call prtris(d,ncoorb,iwr)
            call wrt3(d,nw,iblok,ifild)
         end if
         icomp = icomp + 1
 110  continue
      call revind
      call revise
      call clredx
      return
 6010 format (//1x,'total perturbed density matrices in mo basis')
      end
      subroutine pdens(ss,lstop)
c
c    perturbed density matrices (nuclear motions) in m.o. basis
c
      implicit real*8  (a-h,o-z)
      logical lstop
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
      common/mpshl/ns(maxorb)
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
      dimension ss(*)
      character *8 grhf,oscf
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      logical open
      data grhf/'grhf'/
      data oscf/'oscf'/
      data zero/0.0d0/
      open = scftyp.eq.oscf
c
c
      if (odebug(3)) write (iwr,6010)
      iblll = lensec(mn)
      iblku = iblks + iblll*np
      if (lstop) then
c
c     dump because running out of time
c
         iochf(15) = iochf(1)
         iposd = iochf(1)
         do 20 n = 1 , npfin
            call rdedx(ss,mn,iblku,ifils)
            call wrt3(ss,mn,iposd,ifockf)
            iblku = iblku + iblll
            iposd = iposd + iblll
 20      continue
         npstar = npfin
         irest = 1
         return
      else
c
c
         nsoc = noccb - nocca
         ndpls1 = nocca + 1
         ntpls1 = noccb + 1
         nw = iky(ncoorb+1)
c
c
         iposs = iochf(14)
         iblls = lensec(nw)
c
c     output to section 15 of fockfile.
c     overlap matrix derivatives (m.o. basis) on section 14
c     perturbed wavefunctions on scratchfile
c
         iochf(15) = iochf(1)
         iposd = iochf(15)
         iposdo = iposd + 3*nat*iblls
         call wrt3z(iposd,ifockf,3*nat*iblls)
         if (open) call wrt3z(iposd,ifockf,3*nat*iblls)
         nu = nx + nx
         nopen = nu + mn
         do 170 i = 1 , nat
            do 160 j = 1 , 3
               call rdedx(ss,nw,iposs,ifockf)
               call rdedx(ss(nu+1),mn,iblku,ifils)
               iposs = iposs + iblls
               iblku = iblku + iblll
               do 30 k = 1 , nw
                  if (open) ss(nopen+k) = zero
                  ss(nx+k) = zero
 30            continue
               if (scftyp.eq.grhf) then
                  nr = 0
                  kl = 0
                  do 50 k = 1 , nsa4
                     do 40 l = 1 , k
                        kl = iky(mapie(k)) + mapie(l)
                        ss(nx+kl) = -0.5d0*ss(kl)*
     +                              (fjk(ns(k))+fjk(ns(l)))
                        if (ns(k).ne.ns(l)) then
                           nr = nr + 1
c                          ssss = ss(nx+kl)
                           ss(nx+kl) = ss(nx+kl) + 0.5d0*ss(nu+nr)
     +                                  *(fjk(ns(l))-fjk(ns(k)))
                        end if
 40                  continue
 50               continue
               else
                  if (nocca.ne.0) then
                     do 70 k = 1 , nocca
                        do 60 l = 1 , k
                           kl = iky(mapie(k)) + mapie(l)
                           ss(nx+kl) = -2.0d0*ss(kl)
 60                     continue
 70                  continue
                  end if
                  if (nsoc.ne.0 .and. nocca.ne.0) then
                     do 90 k = ndpls1 , noccb
                        do 80 l = 1 , nocca
                           kl = iky(mapie(k)) + mapie(l)
                           mt = (k-nocca-1)*nocca + l
                           if (open) ss(nopen+kl) = -ss(kl) - ss(nu+mt)
                           ss(nx+kl) = -ss(kl) + ss(nu+mt)
 80                     continue
 90                  continue
                  end if
                  if (nsoc.ne.0) then
                     do 110 k = ndpls1 , noccb
                        do 100 l = ndpls1 , k
                           kl = iky(mapie(k)) + mapie(l)
                           if (open) ss(nopen+kl) = -ss(kl)
                           ss(nx+kl) = -ss(kl)
 100                    continue
 110                 continue
                  end if
                  if (nocca.ne.0) then
                     do 130 k = ntpls1 , num
                        do 120 l = 1 , nocca
                           kl = iky(mapie(k)) + mapie(l)
                           mt = (k-nocca-1)*nocca + l
                           ss(nx+kl) = 2.0d0*ss(nu+mt)
 120                    continue
 130                 continue
                  end if
                  if (nsoc.ne.0) then
                     do 150 k = ntpls1 , num
                        do 140 l = ndpls1 , noccb
                           kl = iky(mapie(k)) + mapie(l)
                           mt = nocca*nvirta + (k-nsoc-1)*nsoc + l -
     +                          nocca
                           if (open) ss(nopen+kl) = ss(nu+mt)
                           ss(nx+kl) = ss(nu+mt)
 140                    continue
 150                 continue
                  end if
               end if
               if (odebug(3)) call prtris(ss(nx+1),ncoorb,iwr)
               call wrt3(ss(nx+1),nw,iposd,ifockf)
               if (open) call wrt3(ss(nopen+1),nw,iposdo,ifockf)
               iposd = iposd + iblls
               if (open) iposdo = iposdo + iblls
 160        continue
 170     continue
         iochf(1) = max(iposd,iposdo)
         call clredx
         irest = 0
c
c
         return
      end if
 6010 format (//1x,'total perturbed density matrices in mo basis')
      end
      subroutine chfdrv(eps,lstop,skipp)
      implicit real*8  (a-h,o-z)
c
c     driving routine for c.h.f.
c     calls chfeqs or chfeq/chfeqv
c
      logical lstop,skipp
      dimension skipp(*)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
      logical scalar
      common/rpaoptions/scalar
c
      character *8 open
      character *6 fnm
      character *6 snm
      dimension eps(mn)
      data open/'open'/
      data fnm/'cphf.m'/
      data snm/'chfdrv'/
c
c
      npstar = 0
      npfin = np
      maxmem = igmem_max_memory()
      memovh = igmem_overhead()
      reqmem = 4*mn+2*mn*mn
      if (mn.le.250 .and. reqmem.lt.maxmem. and. scalar .and. 
     &    .not.CD_active()) then
c
c      use in-core solution
c
         i1 = igmem_alloc_inf(mn,fnm,snm,'b',IGMEM_DEBUG)
         i2 = igmem_alloc_inf(mn,fnm,snm,'cc',IGMEM_DEBUG)
         i3 = igmem_alloc_inf(mn,fnm,snm,'wks1',IGMEM_DEBUG)
         i4 = igmem_alloc_inf(mn,fnm,snm,'wks2',IGMEM_DEBUG)
         i5 = igmem_alloc_inf(mn*mn,fnm,snm,'alpha',IGMEM_DEBUG)
         i6 = igmem_alloc_inf(mn*mn,fnm,snm,'aa',IGMEM_DEBUG)
c
         call chfeqs(eps,qq(ivoff+i1),qq(ivoff+i2),qq(ivoff+i3),
     +               qq(ivoff+i4),qq(ivoff+i5),qq(ivoff+i6),mn,skipp)
c
         call gmem_free_inf(i6,fnm,snm,'aa')
         call gmem_free_inf(i5,fnm,snm,'alpha')
         call gmem_free_inf(i4,fnm,snm,'wks2')
         call gmem_free_inf(i3,fnm,snm,'wks1')
         call gmem_free_inf(i2,fnm,snm,'cc')
         call gmem_free_inf(i1,fnm,snm,'b')
c
       else
c
c       use iterative method
c
         maxc = 50
         if ((mp2 .or. mp3) .and. scftyp.eq.open) then
c
c         use older chfeq
c
            i1 = igmem_alloc_inf(mn,fnm,snm,'b',IGMEM_DEBUG)
            i2 = igmem_alloc_inf(mn,fnm,snm,'utotal',IGMEM_DEBUG)
            i3 = igmem_alloc_inf(mn,fnm,snm,'uold',IGMEM_DEBUG)
            i4 = igmem_alloc_inf(mn*maxc,fnm,snm,'au',IGMEM_NORMAL)
            i5 = igmem_alloc_inf(mn*maxc,fnm,snm,'u',IGMEM_NORMAL)
            write (iwr,*) 'using older chfeq'
            call chfeq(qq(ivoff+1),eps,qq(ivoff+i1),qq(ivoff+i2),
     +                 qq(ivoff+i3),qq(ivoff+i4),qq(ivoff+i5),
     +                 maxc,lstop,skipp)
            call gmem_free_inf(i5,fnm,snm,'u')
            call gmem_free_inf(i4,fnm,snm,'au')
            call gmem_free_inf(i3,fnm,snm,'uold')
            call gmem_free_inf(i2,fnm,snm,'utotal')
            call gmem_free_inf(i1,fnm,snm,'b')
            if (lstop) call clenms('run out of time')
         else
 20         npx = npfin - npstar
            i1 = 0
            i2 = i1 + mn*npx        + memovh
            i3 = i2 + mn*npx        + memovh
            i4 = i3 + maxc*npx      + memovh
            i5 = i4 + maxc*npx      + memovh
            i6 = i5 + maxc*npx      + memovh
            i7 = i6 + maxc*maxc*npx + memovh
            i8 = i7 + npx           + memovh
            ileft = maxmem - i8 - memreq_chfeqv(qq(ivoff+1),skipp,npx)
     &            - memovh
c
c           irmax is the number of columns of the CPHF matrix that can
c           be loaded into core at a time. The whole matrix has 
c           dimension mn*mn.
c
            irmax = ileft/mn
            if (irmax.ge.mn) then
c
c              We can load the whole matrix in core so that is fine.
c
               irmax = mn
            else if (irmax.lt.10) then
c
c              We want to be able to load at least 10 columns at a 
c              time. So keep reducing the number of systems of linear
c              equations until we have enough memory left to do this.
c              If we cannot get 10 columns in core then abort the
c              calculation.
c
               npfin = npfin - 1
               if (npfin.gt.npstar) go to 20
               nreq1 = i8 + mn*10 + memovh
               nreq2 = i8 + mn*32 + memovh
               write (iwr,6010) maxmem , nreq1 , nreq2
               call caserr('insufficient core')
               go to 99999
            else if (irmax.lt.32) then
c
c              We ideally want to be able to load 32 columns at a time.
c              So keep reducing the number of systems of linear 
c              equations until we have enough memory to do this or 
c              until there is 1 system left.
c
               if ((npfin-1).gt.npstar) then
                  npfin = npfin - 1
                  go to 20
               endif
            end if
            i1 = igmem_alloc_inf(mn*npx,fnm,snm,'au',IGMEM_DEBUG)
            i2 = igmem_alloc_inf(mn*npx,fnm,snm,'u',IGMEM_DEBUG)
            i3 = igmem_alloc_inf(maxc*npx,fnm,snm,'b',IGMEM_DEBUG)
            i4 = igmem_alloc_inf(maxc*npx,fnm,snm,'cc',IGMEM_DEBUG)
            i5 = igmem_alloc_inf(maxc*npx,fnm,snm,'uu',IGMEM_DEBUG)
            i6 = igmem_alloc_inf(maxc*maxc*npx,fnm,snm,'uau',
     +                           IGMEM_DEBUG)
            i7 = igmem_alloc_inf(irmax*mn,fnm,snm,'work',IGMEM_DEBUG)
            i8 = igmem_alloc_inf(npx,fnm,snm,'oconv',IGMEM_DEBUG)
            call chfeqv(qq(ivoff+1),eps,qq(ivoff+i1),qq(ivoff+i2),
     +                  qq(ivoff+i7),qq(ivoff+i3),qq(ivoff+i4),
     +                  qq(ivoff+i5),qq(ivoff+i6),
     +                  maxc,skipp,npx,irmax,qq(ivoff+i8))
            call gmem_free_inf(i8,fnm,snm,'oconv')
            call gmem_free_inf(i7,fnm,snm,'work')
            call gmem_free_inf(i6,fnm,snm,'uau')
            call gmem_free_inf(i5,fnm,snm,'uu')
            call gmem_free_inf(i4,fnm,snm,'cc')
            call gmem_free_inf(i3,fnm,snm,'b')
            call gmem_free_inf(i2,fnm,snm,'u')
            call gmem_free_inf(i1,fnm,snm,'au')
            if (npfin.lt.np) then
               npstar = npfin
               npfin = np
               go to 20
            end if
         end if
      end if
99999 return
 6010 format (//1x,'insufficient store for chf equations'//1x,
     +        'store available ',i8/1x,'required - at least ',i8,
     +        ' and preferably ',i8)
      end
      subroutine derlag(q,xerg,yerg,f)
c
c      derivative lagrangian   (integral derivatives only)
c      general scf case
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
      common/mpshl/inshel(maxorb)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      dimension f(11),xerg(11,11),yerg(11,11),q(2)
      call bfnshl(inshel,nsa4)
      nat3 = nat*3
      nfok = nat3 + nat3*njk*2
      i1 = igmem_alloc(num*num)
      i2 = igmem_alloc(nfok*nx)
      call search(iochf(13),ifockf)
      ij = i2
      do 20 n = 1 , nfok
         call reads(q(ij),nx,ifockf)
         call actmot(q(ij),nsa4,mapie,iky)
         ij = ij + nx
 20   continue
      iochf(17) = iochf(1)
      iochf(1) = iochf(1) + lensec(nsa4*nsa4)*nat3
      call search(iochf(17),ifockf)
      lenbig = nx*nat3*2
      do 80 n = 1 , nat3
         call vclr(q(i1),1,num*num)
         ij = i2 + (n-1)*nx
         do 40 i = 1 , nsa4
            do 30 j = 1 , i
               q(i1-1+(j-1)*nsa4+i) = q(ij)*f(inshel(i))*0.5d0
               q(i1-1+(i-1)*nsa4+j) = q(ij)*f(inshel(j))*0.5d0
               ij = ij + 1
 30         continue
 40      continue
         ij = i2 + nx*nat3 + (n-1)*nx
         do 70 i = 1 , nsa4
            do 60 j = 1 , i
               ijj = ij
               ijk = ij + nx*nat3
               do 50 k = 1 , njk
                  q(i1-1+(j-1)*nsa4+i) 
     +            = q(i1-1+(j-1)*nsa4+i) + xerg(inshel(i),k)
     +                           *q(ijj) + yerg(inshel(i),k)*q(ijk)
                  if (i.ne.j) then
                     q(i1-1+(i-1)*nsa4+j) 
     +               = q(i1-1+(i-1)*nsa4+j) + xerg(inshel(j),k)
     +                              *q(ijj) + yerg(inshel(j),k)*q(ijk)
                  end if
                  ijj = ijj + lenbig
                  ijk = ijk + lenbig
 50            continue
               ij = ij + 1
 60         continue
 70      continue
         if (odebug(16)) call prsqm(q(i1),nsa4,nsa4,nsa4,iwr)
         call wrt3s(q(i1),nsa4*nsa4,ifockf)
 80   continue
      call gmem_free(i2)
      call gmem_free(i1)
c
      return
      end
      subroutine effock(q)
c
c     perturbed fock operators for electric field perturbations
c
      implicit real*8  (a-h,o-z)
      logical mpir,mpol,mpir0,mpir1
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/blkin/g(510),nword
      common/craypk/labs(1360)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
      dimension q(*)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/maxlen/maxq
c
      character *8 pnames
      common/tdhfx/pnames(50)
c
c
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
      character *8 dipdmp,polmp2,infrar
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      data        half /       0.5d0/
      data dipdmp /'dipder'/,infrar/'infrared'/
      data polmp2/'polariza'/
      ltri = ikyp(ncoorb)
      length = lensec(ltri)
      icomp = 31
      ifout = 70
c
      np = 0
      maxp = 9
      if (mp2) maxp = 3
      do 20 i = 1 , maxp
         if (ione(i+3).ne.0) np = np + 1
 20   continue
      mpir0 = runtyp.eq.dipdmp .and. mp2
      mpir1 = runtyp.eq.infrar .and. mp2
      mpir = mpir0 .or.mpir1
      mpol = runtyp.eq.polmp2 .and. mp2
      if (mpir .or. mpol) then
         npert = 3
c
c    >>>>> eventually npert = np
c
         nat3 = 3*nat
         nij = nocca*(nocca+1)/2
         nab = nvirta*(nvirta+1)/2
c        ntri = ncoorb*(ncoorb+1)/2
         mpblk(1) = 1
         mpblk(2) = mpblk(1) + lensec(nat3*nat3)
         mpblk(3) = mpblk(2) + lensec(nij*nab)
         mpblk(4) = mpblk(3) + lensec(nij*nab)
         mpblk(5) = mpblk(4) + lensec(ncoorb)
         mpblk(6) = mpblk(5) + lensec(ncoorb*ncoorb*npert)
         mpblk(7) = mpblk(6) + lensec(ncoorb*ncoorb*npert)
         call wrt3z(1,1,mpblk(7))
      end if
c
      call setsto(1360,0,labs)
c
      do 100 n = 1 , 9
         if (.not.(opskip(n))) then
            if (odebug(4)) write (iwr,6010)
c
c     get position of perturbed density matrix on dumpfile
c
            jtype = 0
            call secget(isect(icomp),jtype,iblok)
            call rdedx(q(ltri+1),lds(isect(icomp)),iblok,ifild)
            ija = 0
            do 40 i = 1 , nsa4
               do 30 j = 1 , i
                  ija = ija + 1
                  ij = iky(mapie(i)) + mapie(j)
                  q(ija) = half*q(ij+ltri)
 30            continue
 40         continue
            call vclr(q(ltri+1),1,ltri)
            do 50 i = 1 , ncoorb
               ii = iky(i+1)
               q(ii) = q(ii)*0.5d0
 50         continue
c     total density matrix in q(1)
c
c     scan 2-electron integrals
c
            do 70 i = 1 , mmfile
               iunit = nufile(i)
               call search(kblk(i),iunit)
               call find(iunit)
 60            call get(g,m)
               if (m.ne.0) then
                  if (o255i) then
                     call sgmata(q(ltri+1),q(1))
                  else
                     call sgmata_255(q(ltri+1),q(1))
                  endif
                  call find(iunit)
                  go to 60
               end if
 70         continue
c     get one-electron integrals , and form complete
c     perturbed fock matrix
c
            jtype = 0
            call secget(ipsec(n),jtype,iblok)
            call rdedx(q(1),lds(ipsec(n)),iblok,ifild)
            ij = ltri
            do 90 i = 1 , nsa4
               do 80 j = 1 , i
                  ij = ij + 1
                  ij1 = iky(mapie(i)) + mapie(j)
                  q(ij1) = q(ij1) + q(ij)
 80            continue
 90         continue
c
c     store on the fockfile (ed1)
c
            call secput(isect(ifout),ifout,length,iblok)
            lds(isect(ifout)) = ltri
            call wrt3(q,ltri,iblok,ifild)
            if (odebug(4)) call prtris(q,ncoorb,iwr)
         end if
         icomp = icomp + 1
         ifout = ifout + 1
 100  continue
      call revind
      if (mpir .or. mpol) then
c
c    >>>>>>>>>>>> note npert = 3 temporarily, but
c    >>>>>>>>>>>> but should be pert = np to get
c    >>>>>>>>>>>> quadrupole perturbations as well
c
         npert = 3
         i1 = ltri + 1
         i2 = i1 + ncoorb
         i3 = i2 + ncoorb*ncoorb*npert
         i4 = i3 + ncoorb*ncoorb*npert
         if (i4.gt.maxq) call caserr('insufficient core for makeuf')
         call umatef(q(1),q(i1),q(i2),q(i3),ltri,ncoorb,npert)
      end if
      return
 6010 format (//1x,'perturbed total fock matrices')
      end
      subroutine pfockc(q)
c
c     adds derivative wavefunction term to derivative integral
c     term to make complete derivative of fock operator
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/blkin/g(510),nword
      common/craypk/labs(1)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      common/maxlen/maxq
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      dimension q(*)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      nat3 = nat*3
      ltri = ikyp(ncoorb)
      ibm1 = kblk(1)
      ibm2 = nblk(1)
      idevm = nufile(1)
      length = lensec(ltri)
      ibt = iochf(1)
c
c     perturbed density matrices at section 15
c     derivatives of integrals at section 13
c     complete fock matrix output at section 16 (m.o. basis)
c
      ibh = iochf(13)
      ibd = iochf(15)
      iochf(16) = ibt
c
c      modifications for multipass version
c
      mfok = igmem_max_memory()/(ltri+ltri)
c     mfok is maximum no. of fock matrices per pass , therefore
      npass = (nat3-1)/mfok + 1
      nadd = min(nat3,mfok)
      mi = 1
      ma = nadd
      mzero = igmem_alloc(nadd*ltri)
      mhalf = igmem_alloc(nadd*ltri)
c     mtwo = mhalf + mhalf
c     mthalf = mtwo + mhalf
c     write(iwr,*)' nadd,mi,ma,mhalf',nadd,mi,ma,mhalf
      call setsto(1360,0,labs)
      do 150 ipass = 1 , npass
         i = mhalf 
c    read in batch of density matrices
         do 20 n = mi , ma
            call rdedx(q(i),ltri,ibd,ifockf)
            i = i + ltri
            ibd = ibd + length
 20      continue
c
c     map to active only
c
         ija = 0
         do 50 i = 1 , nsa4
            do 40 j = 1 , i
               ija = ija + 1
               ij = iky(mapie(i)) + mapie(j)
               k = mzero-1
               l = mhalf-1
               do 30 n = mi , ma
                  q(ija+k) = 0.5d0*q(ij+l)
                  k = k + ltri
                  l = l + ltri
 30            continue
 40         continue
 50      continue
c
c     multiply diagonal of density matrix by 0.5
c
         do 70 i = 1 , nsa4
            ii = mzero+iky(i+1)-1
            k = 0
            do 60 n = mi , ma
               q(ii+k) = 0.5d0*q(ii+k)
               k = k + ltri
 60         continue
 70      continue
         call vclr(q(mhalf),1,ltri*nadd)
c     construct fock matrix in q(mhalf)
c
         call search(ibm1,idevm)
         call find(idevm)
         do 80 ib = ibm1 , ibm2
            call get(g,nw)
            if (nword.eq.0) go to 90
            if (nw.eq.0) go to 90
            call find(idevm)
            call sgmatm(q(mhalf),q(mzero),nadd,ltri)
 80      continue
c
c     get term involving integral derivatives
c     and form complete derivative fock matrix
c
 90      i = mzero
         do 100 n = mi , ma
            call rdedx(q(i),ltri,ibh,ifockf)
            i = i + ltri
            ibh = ibh + length
 100     continue
         ija = 0
         do 130 i = 1 , nsa4
            do 120 j = 1 , i
               ija = ija + 1
               ij = iky(mapie(i)) + mapie(j)
               k = mzero-1
               l = mhalf-1
               do 110 n = mi , ma
                  q(ij+k) = q(ij+k) + q(ija+l)
                  k = k + ltri
                  l = l + ltri
 110           continue
 120        continue
 130     continue
         i = mzero
         do 140 n = mi , ma
            call wrt3(q(i),ltri,ibt,ifockf)
            if (odebug(4)) then
               write (iwr,6010) n
               call prtris(q(i),ncoorb,iwr)
            end if
            i = i + ltri
            ibt = ibt + length
 140     continue
         mi = mi + nadd
         ma = ma + nadd
         ma = min(ma,nat3)
 150  continue
      iochf(1) = ibt
      call revise
      call clredx
      if(odebug(30)) then
       write (iwr,6020) iochf(13),iochf(14),iochf(15),iochf(16)
      endif
      call gmem_free(mhalf)
      call gmem_free(mzero)
      if (.not.mp2) return
c
      i2 = igmem_alloc(nw196(5))
      i1 = igmem_alloc((nat3+1)*ltri)
      call symfck(q(i1),q(i2),nshell)
      call gmem_free(i1)
      call gmem_free(i2)
c
      i0 = igmem_alloc(ltri)
      i1 = igmem_alloc(ltri)
      i2 = igmem_alloc(ncoorb)
      i3 = igmem_alloc(ncoorb*ncoorb*nat3)
      call umat(q(i0),q(i1),q(i2),q(i3),ltri,ncoorb,nat3)
      call gmem_free(i3)
      call gmem_free(i2)
      call gmem_free(i1)
      call gmem_free(i0)
      call revise
      return
 6010 format (//5x,'perturbed fock matrix in m.o. basis -- perturbation'
     +        ,i6/)
 6020  format(/1x,'hamfile summary'/
     +         1x,'section 13 at block ',i5/
     +         1x,'section 14 at block ',i5/
     +         1x,'section 15 at block ',i5/
     +         1x,'section 16 at block ',i5/)
      end
      subroutine pfocko(q)
c
c     perturbed density matrices for high-spin open-shell
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      common/blkin/g(510),nword
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
      dimension q(*)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      data zero , half /0.0d0, 0.5d0/
      nat3 = nat*3
      ltri = ikyp(ncoorb)
      ltri0 = igmem_alloc(ltri)
      ltri1 = igmem_alloc(ltri)
      ltri2 = igmem_alloc(ltri)
      ltri3 = igmem_alloc(ltri)
      ibm1 = kblk(1)
      ibm2 = nblk(1)
      idevm = nufile(1)
      length = lensec(ltri)
      ibt = iochf(1)
c
c     perturbed density matrices at section 15
c     derivatives of integrals at section 13
c     complete fock matrix output at section 16 (m.o. basis)
c
      ibh = iochf(13)
      ibd = iochf(15)
      ibdo = ibd + 3*nat*length
      ibk = ibh + 3*nat*length
      iochf(16) = ibt
      do 100 n = 1 , nat3
         call rdedx(q(ltri1),ltri,ibd,ifockf)
c  assume open shell perturbed density matrices are stored
c  after total perturbed density matrices on same section of
c  same file
         call rdedx(q(ltri3),ltri,ibdo,ifockf)
c
c      map active elements only into q(1),q(ltri2+1)
c
         ija = 0
         do 30 i = 1 , nsa4
            do 20 j = 1 , i
               ij = iky(mapie(i)) + mapie(j) - 1
               q(ltri0+ija) = half*q(ij+ltri)
               q(ltri2+ija) = half*q(ij+ltri3)
               ija = ija + 1
 20         continue
 30      continue
         do 40 i = 1 , ncoorb
            ii = iky(i+1)-1
            q(ltri0+ii) = q(ltri0+ii)*0.5d0
            q(ltri2+ii) = q(ltri2+ii)*0.5d0
 40      continue
         do 50 i = 0 , ltri-1
            q(i+ltri1) = zero
            q(i+ltri3) = zero
 50      continue
c  calculate 1/2g:d and put it in q(l+1) - q(2l)
c  calculate 1/2j:do and put it in q(3l+1) - q(4l)
c********************************
c     ibm version of proc2 produces f and k as results.
c
c      the cray version requires the diagonal elements of
c     the density matrices *0.5 and produces f and +k as
c     results.
c
c**********************************
         call search(ibm1,idevm)
         call find(idevm)
         do 60 ib = ibm1 , ibm2
            call get(g,nw)
            if (nword.eq.0) go to 70
            if (nw.eq.0) go to 70
            call find(idevm)
            call proc2f(q(ltri1),q(ltri0),q(ltri3),q(ltri2),
     +                 1.0d0,.true.,.true.)
 60      continue
c  read fa into q(1) to q(l)
c  read 1/2ka into q(2l+1) to q(3l)
 70      call rdedx(q(ltri0),ltri,ibh,ifockf)
         call rdedx(q(ltri2),ltri,ibk,ifockf)
         ij = 0
         do 90 i = 1 , nsa4
            foci = 2.0d0
            if (i.gt.nocca) foci = 1.0d0
            if (i.gt.noccb) foci = 0.0d0
            do 80 j = 1 , i
               focj = 2.0d0
               if (j.gt.nocca) focj = 1.0d0
               if (j.gt.noccb) focj = 0.0d0
               ijt = iky(mapie(i)) + mapie(j) - 1
               ij1 = ijt + ltri0
               ij0 = ij  + ltri1
               ij2 = ijt + ltri2
               ij3 = ij  + ltri2
               q(ij1) = (foci+focj)*(q(ij1)+q(ij0))
     +                  *0.5d0 - (foci*(2.0d0-foci)+focj*(2.0d0-focj))
     +                  *(q(ij3)+q(ij2))*0.5d0
               ij = ij + 1
c
 80         continue
 90      continue
         call wrt3(q(ltri0),ltri,ibt,ifockf)
         ibh = ibh + length
         ibd = ibd + length
         ibk = ibk + length
         ibdo = ibdo + length
         ibt = ibt + length
         if (odebug(4)) write (iwr,6010) n
         if (odebug(4)) call prtris(q(ltri0),ncoorb,iwr)
 100  continue
      call gmem_free(ltri3)
      call gmem_free(ltri2)
      call gmem_free(ltri1)
      call gmem_free(ltri0)
      iochf(1) = ibt
      call revise
      call clredx
      return
 6010 format (//5x,'perturbed fock matrix in m.o. basis -- perturbation'
     +        ,i6/)
      end
      subroutine get1a(a,ifi,ila,mn,n0,ifil)
      implicit real*8  (a-h,o-z)
c
c     used in chfeqv
c
      common/blkin/g(510),nword
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/craypk/labs(680),lab1(340),lab2(340)
      dimension a(ila-ifi+1,mn)
c
      ifi1 = ifi - 1
 20   if (n0.eq.1) then
         call get(g,nw)
c
         if (nw.eq.0) then
c...         EOF disable => next reads (see below)
            n0 = -1
            nword = 0
         end if 
c
         if (nword.eq.0) go to 60
         call unpack(g(num2ep+1),lab1632,labs,numlabp)
         n2 = 1
         do n = 1 , nword
            lab1(n) = labs(n2)
            lab2(n) = labs(n2+1)
            n2 = n2 + 2
         enddo
      end if
      if (n0.lt.0) go to 60
      nnext = n0
      if (lab1(nword).le.ila) then
         do 40 n = nnext , nword
            a(lab1(n)-ifi1,lab2(n)) = g(n)
 40      continue
      else
         do 50 n = nnext , nword
            n0 = n
            if (lab1(n).gt.ila) go to 60
            a(lab1(n)-ifi1,lab2(n)) = g(n)
 50      continue
      end if
      call find(ifil)
      n0 = 1
      go to 20
 60   do 80 i = ifi + 1 , ila
         do 70 j = ifi , i - 1
            a(j-ifi1,i) = a(i-ifi1,j)
 70      continue
 80   continue
      return
      end
      subroutine get2a(a,ifi,ila,mn,ifil)
      implicit real*8  (a-h,o-z)
c
c     used in chfeqv
c
      common/blkin/g(510),nword
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/craypk/labs(1360)
      dimension a(mn,ila-ifi+1)
c
      iblk1 = 1
      call search(iblk1,ifil)
      call find(ifil)
 20   call get(g,nw)
      if (nw.eq.0 .or. nword.eq.0) return
      call unpack(g(num2ep+1),lab1632,labs,numlabp)
      do 30 n = 1 , nword
         lab1 = labs(n+n-1)
         lab2 = labs(n+n)
         if (lab2.ge.ifi .and. lab2.le.ila) a(lab1,lab2-ifi+1) = g(n)
 30   continue
      call find(ifil)
      go to 20
      end
      subroutine lgrhf(a,ibzeta,alphax,betax,fx)
c
c     generalised lagrangian for grhf
c
      implicit real*8  (a-h,o-z)
      dimension a(*),alphax(11,11),betax(11,11),fx(11)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      common/craypk/labs(1360)
      common/maxlen/maxq
      common/blkin/g(510),nword
      common/mpshl/ns(maxorb)
      logical exist
c     data m21/21/
c
      call bfnshl(ns,nsa4)
c
c     read m.o. t+v matrix off dumpfile
c
      call secloc(isect(21),exist,isec21)
      if (exist) then
         call rdedx(a(1),lds(isect(21)),isec21,ifild)
      else
         call caserr('transformed t+v matrix required')
      end if
c
c     multiply by 1-electron occupation numbers
c
      call actmot(a,nsa4,mapie,iky)
      ij = 0
      do 40 i = 1 , nsa4
         do 30 j = 1 , i
            ij = ij + 1
            fij = 0.5d0*a(ij)
            nn = ij
            do 20 n = 1 , njk1
               a(nn) = fij*fx(n)
               nn = nn + nx
 20         continue
 30      continue
 40   continue
c
c     clear label buffer, labs (i205)
      call setsto(1360,0,labs)
c
c     loop over the transformed two-electron integrals
c     which are input from ed6 ( default)
c
      do 230 ifile = 1 , mmfile
         mblkk = kblk(ifile)
         idevm = nufile(ifile)
c        lblkm = nblk(ifile)
c
         call search(mblkk,idevm)
         call find(idevm)
c
c     read block of integrals into /blkin/
c
 50      call get(g(1),nw)
         if (nword.gt.0) then
            if (nw.gt.0) then
               call find(idevm)
c
c     loop over integrals in a block
c
               call unpack(g(num2e+1),lab816,labs,numlab)
               do 220 int = 1 , nword
c
c     unpack the labels
c
                  kk2 = (int+int) + (int+int)
                  i = labs(kk2-2)
                  j = labs(kk2-3)
                  k = labs(kk2  )
                  l = labs(kk2-1)
                  gg = g(int)
                  it = ijkltp(i,j,k,l)
                  go to (60,60,60,80,120,160,140,180,200,100,100,220,
     +                   220,220) , it
c     <ii/ii>,<ii/il>,<il/ll>
 60               il = iky(i) + l
                  do 70 n = 1 , njk1
                     a(il) = a(il) + gg*(alphax(n,ns(j))+betax(n,ns(j)))
                     il = il + nx
 70               continue
                  go to 220
c     <ii/kk>
 80               kk = iky(k) + k
                  do 90 n = 1 , njk1
                     a(kk) = a(kk) + gg*alphax(n,ns(i))
                     kk = kk + nx
 90               continue
c     <ij/kk>
 100              ij = iky(i) + j
                  do 110 n = 1 , njk1
                     a(ij) = a(ij) + gg*alphax(n,ns(k))
                     ij = ij + nx
 110              continue
                  go to 220
c      <ij/ij>
 120              ii = iky(i) + i
                  do 130 n = 1 , njk1
                     a(ii) = a(ii) + gg*betax(n,ns(j))
                     ii = ii + nx
 130              continue
c     <ij/il>
 140              jl = iky(j) + l
                  do 150 n = 1 , njk1
                     a(jl) = a(jl) + gg*betax(n,ns(i))
                     jl = jl + nx
 150              continue
                  go to 220
c    <ii/kl>
 160              kl = iky(k) + l
                  do 170 n = 1 , njk1
                     a(kl) = a(kl) + gg*alphax(n,ns(i))
                     kl = kl + nx
 170              continue
                  go to 220
c     <ij/jl>
 180              il = iky(i) + l
                  do 190 n = 1 , njk1
                     a(il) = a(il) + gg*betax(n,ns(j))
                     il = il + nx
 190              continue
                  go to 220
c     <ij/kj>
 200              ik = iky(i) + k
                  do 210 n = 1 , njk1
                     a(ik) = a(ik) + gg*betax(n,ns(j))
                     ik = ik + nx
 210              continue
 220           continue
               go to 50
            end if
         end if
 230  continue
      ij = 1
      do 240 n = 1 , njk1
         if (odebug(16)) call prtris(a(ij),nsa4,iwr)
         ij = ij + nx
 240  continue
      call wrt3(a,nx*njk1,ibzeta,ifils)
      return
      end
      subroutine lgrhfm(a,ibeta,alphax,betax,fx)
c
c      lagrangian in mo basis --- from integral list
c      for grhf
c
      implicit real*8  (a-h,o-z)
      dimension a(*),alphax(11,11),betax(11,11),fx(11)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      common/maxlen/maxq
      common/mpshl/ns(maxorb)
      common/craypk/labs(1360)
      common/blkin/g(510),nword
      logical exist
c     data m21/21/
c
      call bfnshl(ns,nsa4)
c
c     read m.o. t+v matrix off dumpfile
c
      call secloc(isect(21),exist,isec21)
      if (exist) then
        call rdedx(a(1),lds(isect(21)),isec21,ifild)
      else
        call caserr('transformed t+v matrix required')
      end if
c
c     multiply by 1-electron occupation numbers
c
      call actmot(a,nsa4,mapie,iky)
      ij = 0
      do 30 i = 1 , nsa4
         fi = 0.5d0*fx(ns(i))
         do 20 j = 1 , i
            ij = ij + 1
            a(ij) = a(ij)*fi
 20      continue
 30   continue
c
c     clear label buffer, labs (i205)
      call setsto(1360,0,labs)
c
c     loop over the transformed two-electron integrals
c     which are input from ed6 ( default)
c
      do 140 ifile = 1 , mmfile
         mblkk = kblk(ifile)
         idevm = nufile(ifile)
c
         call search(mblkk,idevm)
         call find(idevm)
c
c     read block of integrals into /blkin/
c
 40      call get(g(1),nw)
         if (nword.gt.0) then
            if (nw.gt.0) then
               call find(idevm)
c
c     loop over integrals in a block
c
               call unpack(g(num2e+1),lab816,labs,numlab)
               do 130 int = 1 , nword
c
c     unpack the labels
c
                  kk2 = (int+int) + (int+int)
                  i = labs(kk2-2)
                  j = labs(kk2-3)
                  k = labs(kk2  )
                  l = labs(kk2-1)
                  gg = g(int)
                  it = ijkltp(i,j,k,l)
                  go to (50,50,50,60,80,100,90,110,120,70,70,130,130,
     +                   130) , it
 50               il = iky(i) + l
                  a(il) = a(il)
     +                    + gg*(alphax(ns(i),ns(j))+betax(ns(i),ns(j)))
                  go to 130
 60               kk = iky(k) + k
                  a(kk) = a(kk) + gg*alphax(ns(k),ns(i))
 70               ij = iky(i) + j
                  a(ij) = a(ij) + gg*alphax(ns(i),ns(k))
                  go to 130
 80               ii = iky(i) + i
                  a(ii) = a(ii) + gg*betax(ns(i),ns(j))
 90               jl = iky(j) + l
                  a(jl) = a(jl) + gg*betax(ns(j),ns(i))
                  go to 130
 100              kl = iky(k) + l
                  a(kl) = a(kl) + gg*alphax(ns(k),ns(i))
                  go to 130
 110              il = iky(i) + l
                  a(il) = a(il) + gg*betax(ns(i),ns(j))
                  go to 130
 120              ik = iky(i) + k
                  a(ik) = a(ik) + gg*betax(ns(i),ns(j))
 130           continue
               go to 40
            end if
         end if
 140  continue
      if (odebug(16)) call prtris(a,nsa4,iwr)
      call wrt3(a,nx,ibeta,ifils)
      return
      end
      subroutine umat(f,s,e,u,ltri,norb,nat3)
c
c     produces full u matrix (including  redundant terms)
c     for nuclear perturbations
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension f(ltri),s(ltri),e(norb),u(norb,norb,nat3)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      ind(i,j) = iky(max(i,j)) + min(i,j)
c
      ifile1 = 1
c
c
      call secget(isect(9),9,isec9)
      call rdedx(e,ncoorb,isec9,ifild)
c     perturbed fock matrix at section 16 , m.o. basis
c     derivatives of overlap at section 14
c
      ibf = iochf(16)
      ibs = iochf(14)
      length = lensec(ltri)
      call vclr(u(1,1,nat3-2),1,3*norb*norb)
      do 40 n = 1 , nat3 - 3
         call rdedx(f,ltri,ibf,ifockf)
         call rdedx(s,ltri,ibs,ifockf)
         ibf = ibf + length
         ibs = ibs + length
         do 30 i = 1 , ncoorb
            do 20 j = 1 , ncoorb
               diff = e(j) - e(i)
c
               if (dabs(diff).gt.1.0d-3) then
                  u(i,j,n) = (f(ind(i,j))-s(ind(i,j))*e(j))/diff
               else
                  u(i,j,n) = -s(ind(i,j))*0.5d0
               end if
c
               if (dabs(u(i,j,n)).le.1.0d-15) u(i,j,n) = 0.0d0
c
 20         continue
 30      continue
 40   continue
c
      call wrt3(u,norb*norb*nat3,mpblk(6),ifile1)
c
      ibf = iochf(16)
      do 70 n = 1 , nat3 - 3
         call rdedx(f,ltri,ibf,ifockf)
         ibf = ibf + length
c
         ipiq = 0
         do 60 ipp = 1 , norb
            do 50 iq = 1 , ipp
               ipiq = ipiq + 1
               diff = e(ipp) - e(iq)
               if (dabs(diff).gt.1.0d-3) then
                  u(ipp,iq,n) = 0.0d0
               else
                  u(ipp,iq,n) = f(ipiq) + u(ipp,iq,n)*e(ipp)
     +                          + u(iq,ipp,n)*e(iq)
               end if
               u(iq,ipp,n) = u(ipp,iq,n)
 50         continue
 60      continue
 70   continue
      call wrt3(e,norb,mpblk(4),ifile1)
      call wrt3(u,norb*norb*nat3,mpblk(5),ifile1)
      return
      end
      subroutine umatef(f,e,eder,u,ltri,norb,npert)
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension f(ltri),e(norb),eder(norb,norb,npert),
     1   u(norb,norb,npert)
c
c     makes all the u-vector elements, including redundant pairs
c     for electric field perturbations, and perturbed fock
c     operator elements
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
      character*10 charwall
c
      ind(i,j) = iky(max(i,j)) + min(i,j)
c
c
      call secget(isect(9),9,isec9)
      call rdedx(e,lds(isect(9)),isec9,ifild)
c
c      perturbed fock matrices from section 70 thru 72 of dumpfile
c
      ifout = 70
      do 80 n = 1 , npert
         if (ione(n+3).ne.0) then
            call secget(isect(ifout),ifout,ibf)
            call rdedx(f,ltri,ibf,ifild)
c
c     perturbed eigenvalues
c
            call vclr(eder(1,1,n),1,ncoorb*ncoorb)
            do 30 i = 1 , nocca
               do 20 j = 1 , i
                  ij = iky(i) + j
                  eder(i,j,n) = f(ij)
                  eder(j,i,n) = eder(i,j,n)
 20            continue
 30         continue
            do 50 i = nocca + 1 , ncoorb
               do 40 j = nocca + 1 , i
                  ij = iky(i) + j
                  eder(i,j,n) = f(ij)
                  eder(j,i,n) = eder(i,j,n)
 40            continue
 50         continue
            call vclr(u(1,1,n),1,ncoorb*ncoorb)
            do 70 ia = nocca + 1 , ncoorb
               do 60 i = 1 , nocca
                  diff = e(ia) - e(i)
                  u(ia,i,n) = -f(ind(ia,i))/diff
                  u(i,ia,n) = -u(ia,i,n)
 60            continue
 70         continue
         end if
         ifout = ifout + 1
 80   continue
      ifile1 = 1
      call wrt3(e,ncoorb,mpblk(4),ifile1)
      call wrt3(eder,ncoorb*ncoorb*npert,mpblk(5),ifile1)
      call wrt3(u,ncoorb*ncoorb*npert,mpblk(6),ifile1)
      write (iwr,6010) cpulft(1) ,charwall()
 6010 format(/1x,'derivative eigenvalue generation complete at ',
     +  f8.2,' seconds',a10,' wall')
      return
      end
      subroutine nrmap(mnr,nocca,noccb,nsa,iky,mapie)
      implicit real*8  (a-h,o-z)
      dimension mnr(*),mapie(*),iky(*)
c
c     nsa is total number of active m.o.'s
c     nocca is number of doubly occupied
c     noccb is total occupied
c
c     maps from lower triangle to
c     non-redundant pairs
c
      nsoc = noccb - nocca
      nvirta = nsa - noccb
      ntpls1 = noccb + 1
      ndpls1 = nocca + 1
      if (nocca.ne.0) then
         do 30 j = 1 , nocca
            do 20 i = ndpls1 , nsa
               it = (i-ndpls1)*nocca + j
               mnr(it) = iky(mapie(i)) + mapie(j)
 20         continue
 30      continue
      end if
      if (noccb.ne.nocca) then
         do 50 j = ndpls1 , noccb
            do 40 i = ntpls1 , nsa
               it = nvirta*nocca + (i-nsoc-1)*nsoc + j - nocca
               mnr(it) = iky(mapie(i)) + mapie(j)
 40         continue
 50      continue
      end if
      return
      end
      subroutine mxmau1(a,u,au,mn,ifi,ila,irange,np,irmax)
      implicit real*8  (a-h,o-z)
c
c     used in chfeqv
c
      dimension a(irange,mn),u(mn,np),au(mn,np)
      ir = irange
      jr = ila
      if (ir.le.0 .or. jr.le.0 .or. np.le.0) return
      call mxmb(a,1,irange,u,1,mn,au(ifi,1),1,mn,ir,jr,np)
      jr = ir
      ir = ifi - 1
      if (ir.gt.0) call mxmb(a,irange,1,u(ifi,1),1,mn,au,1,mn,ir,jr,np)
      ifi = ila + 1
      ila = ila + irmax
      if (ila.gt.mn) ila = mn
      return
      end
      subroutine mxmau2(a,u,au,mn,ifi,ila,ir,np,irmax)
      implicit real*8  (a-h,o-z)
c
c     used in chfeqv
c
      dimension a(mn,ir),u(mn,np),au(mn,np)
      if (ir.le.0 .or. np.le.0) return
      call mxmb(a,1,mn,u(ifi,1),1,mn,au,1,mn,mn,ir,np)
      ifi = ila + 1
      ila = ila + irmax
      if (ila.gt.mn) ila = mn
      return
      end
      subroutine symfck(qq,iso,nshels)
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension qq(*),iso(nshels,*)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/mpshl/ns(maxorb)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/bufb/ptr(3,144),ict(maxat,8)
      common/symmos/imos(8,maxorb)
      data one/1.0d0/
c
      ind(i,j) = iky(max(i,j)) + min(i,j)
c
c     write(iwr,*) 'entered mpfsym..number of operations ',nt
      if (nt.eq.1) return
c
      nav = lenwrd()
c
      call rdedx(ptr(1,1),nw196(1),ibl196(1),ifild)
      call readi(iso,nw196(5)*nav,ibl196(5),ifild)
      do 40 ii = 1 , nshell
         ic = katom(ii)
         do 30 it = 1 , nt
            id = iso(ii,it)
            ict(ic,it) = katom(id)
 30      continue
 40   continue
c
      ntri = ncoorb*(ncoorb+1)/2
c     n2 = ncoorb*ncoorb
      iblen = lensec(ntri)
      iblok = iochf(16)
      an = one/dble(nt)
      ioff = ntri + 1
      nat3 = nat*3
c     n3n = nat3
c     read in m.o. perturbed fock matrix
      do 50 n = 1 , nat3
         call rdedx(qq(ioff),ntri,iblok,ifockf)
         iblok = iblok + iblen
         ioff = ioff + ntri
 50   continue
c     loop over matrices
      iblok = iochf(16)
      do 130 n = 1 , nat
         do 120 nc = 1 , 3
            ioff = ((n-1)*3+nc)*ntri
c     copy matrix for atom n component nc into work area
c     this is equivalent to identity
            do 60 i = 1 , ntri
               qq(i) = qq(ioff+i)
 60         continue
c     work along the elements of this matrix
            do 100 iip = 1 , ncoorb
               do 90 iiq = 1 , iip
                  ipq = ind(iip,iiq)
c     loop over symmetry operations
c     except identity
                  do 80 iop = 2 , nt
                     isign = imos(iop,iip)*imos(iop,iiq)
                     sign = dble(isign)
                     niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                     ioff = (niop-1)*3*ntri
                     npnc = (iop-1)*3 + nc
                     do 70 k = 1 , 3
                        ioff = ioff + ntri
                        qq(ipq) = ptr(k,npnc)*sign*qq(ioff+ipq)
     +                            + qq(ipq)
 70                  continue
 80               continue
 90            continue
 100        continue
            do 110 i = 1 , ntri
               qq(i) = an*qq(i)
 110        continue
            call wrt3(qq(1),ntri,iblok,ifockf)
            iblok = iblok + iblen
 120     continue
 130  continue
      return
      end
      subroutine wamat(a,mnmn,ifil,nst,nfin,odebug,iww)
c
c     writes out a-matrix ( hessian ) as constructed by chfcls
c     ( and equivalent open-shell routines)
c     -------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension a(*)
      logical odebug
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      common/craypk/labs(1360)
      common/blkin/g(510),nint
      data m0/0/
c
      call setsto(680,0,labs)
c
      nint = m0
      nstp1 = nst + 1
      itri = m0
      do 30 i = nstp1 , nfin
         do 20 j = 1 , i
            itri = itri + 1
            val = a(itri)
            if (val.ne.0.0d0) then
               nint = nint + 1
               g(nint) = val
               labs(nint+nint-1) = i
               labs(nint+nint) = j
               if (nint.ge.num2e) then
                  call pack(g(num2ep+1),lab1632,labs,numlabp)
                  call put(g,m511,ifil)
                  call setsto(680,0,labs)
                  nint = m0
               end if
               if (.not.(.not.lcpf .and. .not.cicv .and. (.not.cicx)))
     +             then
                  if (i.ne.j) then
                     nint = nint + 1
                     g(nint) = val
                     labs(nint+nint-1) = j
                     labs(nint+nint) = i
                     if (nint.ge.num2e) then
                        call pack(g(num2ep+1),lab1632,labs,numlabp)
                        call put(g,m511,ifil)
                        call setsto(680,0,labs)
                        nint = m0
                     end if
                  end if
               end if
            end if
 20      continue
 30   continue
      if (nint.ne.m0) then
         call pack(g(num2ep+1),lab1632,labs,numlabp)
         call put(g,m511,ifil)
         call setsto(680,0,labs)
         nint = m0
      end if
      if (.not.(lcpf .or. cicv .or. cicx)) then
         if (nfin.eq.mnmn) call put(g,m0,ifil)
      end if
      if (odebug) write (iww,6010)
      if (odebug) call prtris(a,mnmn,iww)
      return
 6010 format (///1x,'a-matrix')
      end
      subroutine polasm(b,u,prop,ndim,npdim)
c
c    assemble and print out polarisabilities
c
      implicit real*8  (a-h,o-z)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      character *8 pnames
      common/tdhfx/pnames(50)
c
c
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      dimension prop(npdim,npdim),orb(100,6)
      dimension b(ndim,np),u(ndim,np)
      character *8 closed
      data closed /'closed'/
c
c
      call vclr(prop,1,npdim*npdim)
c
c     read in pertubation vectors in b
c
      call rdedv(b,ndim,np,iblks,ifils)
c
c     read the solution vectors in u
c
      call rdedvs(u,ndim,np,ifils)
c
      do 30 i = 1 , np
         do 20 j = 1 , np
            zz = -4.0d0*ddot(mn,u(1,i),1,b(1,j),1)
            prop(i,j) = zz
 20      continue
 30   continue
c
c
      if(oprn(40)) call prpol0(prop,npdim)
      if (scftyp.ne.closed) return
      if (nocca.gt.100) return
      call vclr(orb,1,600)
      nnn = 0
      do 70 k = 1 , 3
         do 60 l = 1 , k
            nnn = nnn + 1
            mmm = 0
            do 50 j = nocca + 1 , nsa4
               do 40 i = 1 , nocca
                  mmm = mmm + 1
                  orb(i,nnn) = orb(i,nnn) - 4.0d0*u(mmm,k)*b(mmm,l)
 40            continue
 50         continue
 60      continue
 70   continue
      if(oprn(26) .and. oprn(40)) then
       write (iwr,6010)
       do 80 i = 1 , nocca
          write (iwr,6020) i , (orb(i,j),j=1,6)
 80    continue
      endif
      return
 6010 format (/
     + 10x,'***********************************************'/
     + 10x,'* m.o. contributions to dipole polarizability *'/
     + 10x,'***********************************************'//
     + 1x,'m.o.',8x,'xx',16x,'xy',16x,'yy',16x,'xz',16x,'yz',16x,
     +        'zz'/)
 6020 format (1x,i4,6f18.8)
      end
      subroutine poldrv(q,iq)
c
c      driving routine for polarisability calculations
c      -----------------------------------------------
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      logical lmag,ldynam,lbig,lstop,skipp
      dimension skipp(100)
c
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c
      character *8 pnames
      common/tdhfx/pnames(50)
c
c
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common /maxlen/ maxq
      common/mpshl/ns(maxorb)
      logical lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis
      common/scfblk/en,etot,ehf,sh1(2),sh2(2),gap1(2),gap2(2),
     1              d12,d13,d23,canna,cannb,cannc,fx,fy,fz,
     2              lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis,
     3              ncyc,ischm,lock,maxit,nconv,npunch,lokcyc
      common /small / abcis(30),weight(30)
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
c     COMMON/DRIVE_DFT/ contains options that specify how GAMESS-UK
c     should drive the CCP1 DFT module.
c
      integer KS_AO,     KS_MO,     KS_AOMO
      integer KS_GRD_AO, KS_GRD_MO, KS_GRD_AOMO
      integer KS_RHS_AO, KS_RHS_MO, KS_RHS_AOMO
      integer KS_LHS_AO, KS_LHS_MO, KS_LHS_AOMO
      integer KS_DT_AO,  KS_DT_MO,  KS_DT_AOMO
      integer KS_HES_AO, KS_HES_MO, KS_HES_AOMO
      integer KS_DX_AO,  KS_DX_MO,  KS_DX_AOMO
      parameter (KS_AO    =01,KS_MO    =02,KS_AOMO    =03)
      parameter (KS_GRD_AO=04,KS_GRD_MO=05,KS_GRD_AOMO=06)
      parameter (KS_RHS_AO=07,KS_RHS_MO=08,KS_RHS_AOMO=09)
      parameter (KS_LHS_AO=10,KS_LHS_MO=11,KS_LHS_AOMO=12)
      parameter (KS_DT_AO =13,KS_DT_MO =14,KS_DT_AOMO =15)
      parameter (KS_HES_AO=16,KS_HES_MO=17,KS_HES_AOMO=18)
      parameter (KS_DX_AO =19,KS_DX_MO =20,KS_DX_AOMO =21)
c
c     if ks_bas.eq.KS_AO           call CD_energy_ao
c     if ks_bas.eq.KS_MO           call CD_energy_mo
c     if ks_bas.eq.KS_AOMO         call CD_energy
c
c     if ks_grd_bas.eq.KS_GRD_AO   call CD_force_ao
c     if ks_grd_bas.eq.KS_GRD_MO   call CD_force_mo
c     if ks_grd_bas.eq.KS_GRD_AOMO call CD_force
c
c     if ks_rhs_bas.eq.KS_RHS_AO   call CD_chf_rhs_ao
c     if ks_rhs_bas.eq.KS_RHS_MO   call CD_chf_rhs_mo
c     if ks_rhs_bas.eq.KS_RHS_AOMO call CD_chf_rhs
c
c     if ks_lhs_bas.eq.KS_LHS_AO   call CD_chf_lhs_ao
c     if ks_lhs_bas.eq.KS_LHS_MO   call CD_chf_lhs_mo
c     if ks_lhs_bas.eq.KS_LHS_AOMO call CD_chf_lhs
c
c     if ks_dt_bas.eq.KS_DT_AO     call CD_dksm_ao
c     if ks_dt_bas.eq.KS_DT_MO     call CD_dksm_mo
c     if ks_dt_bas.eq.KS_DT_AOMO   call CD_dksm
c
c     if ks_hes_bas.eq.KS_HES_AO   call CD_hess_ao
c     if ks_hes_bas.eq.KS_HES_MO   call CD_hess_mo
c     if ks_hes_bas.eq.KS_HES_AOMO call CD_hess
c
c     if ks_dx_bas.eq.KS_DX_AO     call CD_dksm_exp_ao
c     if ks_dx_bas.eq.KS_DX_MO     call CD_dksm_exp_mo
c     if ks_dx_bas.eq.KS_DX_AOMO   call CD_dksm_exp
c
      integer ks_bas,     ks_grd_bas, ks_hes_bas
      integer ks_rhs_bas, ks_lhs_bas, ks_dt_bas
      integer ks_dx_bas
      logical ogeompert
c
      common/drive_dft/ks_bas, ks_grd_bas, 
     &       ks_rhs_bas, ks_lhs_bas, ks_dt_bas,
     &       ks_hes_bas, ks_dx_bas,
     &       ogeompert
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
      logical ogeompert_save
      dimension iq(*)
      dimension q(*)
      character *8 closed,oscf,grhf,open
      data closed,oscf,grhf,open/'closed','oscf','grhf','open'/
c
      ierror = CD_set_2e()
      lenx = lensec(nx)*nat*15
      call cpuwal(begin,ebegin)
      call wrt3z(iblks,ifils,lenx)
      if (scftyp.eq.oscf) then
         if((dabs(canna-2.0d0).gt.1.0d-6).or.
     +      (dabs(cannb).gt.1.0d-6).or.
     +      (dabs(cannc+2.0d0).gt.1.0d-6)) then
            write (iwr,6010)
            call caserr('stop')
         end if
      end if
c...  structure of section
c     a) /tdhf/ 1 block
      lds(isect(52)) = 31 + lenint(192)
      ldsect(isect(52)) = 50
c     b) trans energies at least 1 block
c     c) trans moments at least 1 block
c     d) properties nfreq+1 blocks
c
c   ======= next bit mainly for setting up dynamic props.
c...  restore if necessary
      ifreq = 0
      if (irestp.ne.0) then
         call secget(isect(52),52,iblkj)
         call rdchr(pnames,ldsect(isect(52)),iblkj,ifild)
         call reads(freq,lds(isect(52)),ifild)
         write (iwr,6020) isect(52) , irestp , ifreq
      end if
c
      lbig = npole.lt.0 .or. npole.eq.999
      if (irestp.eq.4) then
         ieps = igmem_alloc(mn)
         call tdchf(q(ieps),lbig)
         call gmem_free(ieps)
      else
         if (npole.eq.999) npole = 0
         npole = iabs(npole)
         ldynam = nfreq.ne.0 .or. npole.ne.0
         lmag = ldynam
         npdim = max(np,9)
         len = 2 + lensec(npole) + lensec(npole*np) + (iabs(nfreq)+1)
     +         *lensec(npdim*npdim)
         call secput(isect(52),52,len,iblkj)
         call revind
         call wrtc(pnames,ldsect(isect(52)),iblkj,ifild)
         call wrt3s(freq,lds(isect(52)),ifild)
c
c
c  ======================================
         imap  = igmem_alloc(lenint(nx))
         iimap = lenrel(imap-1)+1
         if (scftyp.eq.open) then
            write(*,*)'UHF polarizabilities are not supported'
            call caserr('No UHF polarizabilities')
         endif
         if (scftyp.ne.grhf) then
            if (noccb.lt.nocca) then
c              Original CADPAC code assumes noccb >= nocca
               call caserr('Code broken try GRHF instead')
            endif
            mn = noccb*nvirta + (noccb-nocca)*nocca
            call grhfbl(scftyp)
            call bfnshl(ns,nsa4)
            call nrmapo(iq(iimap),nocca,noccb,nsa4,iky,mapie)
         else
            call bfnshl(ns,nsa4)
            call ijmapr(nsa4,ns,iq(iimap),mn,mnx)
         end if
         ieps = igmem_alloc(mn)
c
c     sort out required integrals ( form a-matrix )
c
         i01 = igmem_alloc_all(maxa)
         if (scftyp.eq.closed) call chfcls(q(i01),maxa)
         if (scftyp.eq.oscf) call chfops(q(i01),maxa)
         call gmem_free(i01)
         if (scftyp.eq.grhf) then
            lenblk = lensec(mn)
            ibeta = iblks + lenblk*nat*6
            ibzeta = ibeta + lensec(nx)
c           i2 = i1 + nx*njk1
            i01 = igmem_alloc(nx)
            i1  = igmem_alloc(nx*njk1)
            i2  = igmem_alloc_all(maxa)
            call lgrhfm(q(i01),ibeta,erga,ergb,fjk)
            call lgrhf(q(i1),ibzeta,erga,ergb,fjk)
            call chfgrs(q(ieps),iq(iimap),q(i01),q(i1),q(i2),
     +                  ibeta,ibzeta,maxa)
            call gmem_free(i2)
            call gmem_free(i1)
            call gmem_free(i01)
         end if
         i01 = igmem_alloc_all(maxa)
         if (lmag) call chficl(q(i01),maxa)
         call gmem_free(i01)
c
c
c   right-hand-side of equations
c
         i01 = igmem_alloc(mn)
         i1  = igmem_alloc(nx)
         call rhsemp(q(ieps),q(i01),q(i1),iq(iimap))
         call gmem_free(i1)
         call gmem_free(i01)
c    solve chf equations
c
         lstop = .false.
         do 30 n = 1 , np
            skipp(n) = .false.
 30      continue
         ogeompert_save = ogeompert
         ogeompert = .false.
         call chfdrv(q(ieps),lstop,skipp)
         ogeompert = ogeompert_save
         if (lstop) then
            write (iwr,*)
     +         'insufficient time - restart polarizability from'
     +         , ' beginning'
            call timana(22)
            call clenms('stop')
         end if
         if (ifreq.gt.0) then
            call tdchf(q(ieps),lbig)
         else
c
c     construct perturbed density matrices
c
            if (ldens .or. ldiag) then
               i01 = igmem_alloc(mn)
               i1  = igmem_alloc(nx)
               call efdens(q(i01),q(i1))
               call gmem_free(i1)
               call gmem_free(i01)
            endif
c
c     assemble polarizability tensors
c
            npdim = max(np,9)
            i1 = i01 + mn*np
            i2 = i1 + mn*np
c           i3 = i2 + npdim*npdim
            i01 = igmem_alloc(mn*np)
            i1  = igmem_alloc(mn*np)
            i2  = igmem_alloc(npdim*npdim)
            call polasm(q(i01),q(i1),q(i2),mn,npdim)
            leng = lensec(mn)*np
            call secput(isect(65),65,leng,iblok)
            lds(isect(65)) = mn*np
            call search(iblok,ifild)
            do 40 i = 1 , np
               call wrt3s(q(i1+(i-1)*mn),mn,ifild)
 40         continue
            call revind
            call gmem_free(i2)
            call gmem_free(i1)
            call gmem_free(i01)
c
c     construct complete perturbed fock matrix if requested
c
            if (ldiag) then
               i01 = igmem_alloc_all(maxa)
               call effock(q(i01))
               call gmem_free(i01)
            endif
c           call timit(1)
            if (npole.ne.0 .or. nfreq.ne.0) then
               if (scftyp.ne.closed) then
                  write (iwr,6030)
               else
                  call tdchf(q(ieps),lbig)
               end if
            end if
         end if
         call gmem_free(ieps)
         call gmem_free(imap)
      end if
      irestp = 0
      nfreq = 0
      npole = 0
      ierror = CD_reset_2e()
      call revise
      call secget(isect(52),52,iblkj)
      call wrtc(pnames,ldsect(isect(52)),iblkj,ifild)
      call wrt3s(freq,lds(isect(52)),ifild)
      call delfil(nofile(1))
      call timana(22)
      return
 6010 format (/1x,'oscf chf equations only work with',
     +        ' canonicalisation 2.0, 0.0, -2.0'/1x,
     +        '------     see manual for details --------')
 6020 format (' restart information restored from section',i4,
     +        ' of dumpfile'/' irestp =',i3/' ifreq =',i3/)
 6030 format (//1x,'dynamic properties only available for',
     +        ' closed-shell systems')
      end
      subroutine prpol1(prop,npdim,mini,maxi,nstart,fac,iw)
      implicit real*8  (a-h,o-z)
c
c     another printing routine for polarisabilities
c
      character *8 pbuff
c
      character *8 pnames
      common/tdhfx/pnames(50)
c
c
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
      dimension prop(npdim,npdim)
      dimension pbuff(6),pout(6)
      nnn = 0
      do 20 i = mini , maxi
         if (.not.(opskip(i))) then
            nnn = nnn + 1
            pbuff(nnn) = pnames(i)
         end if
 20   continue
      if (nnn.eq.0) return
      write (iw,6010) (pbuff(i),i=1,nnn)
      nj = 0
      do 40 j = 1 , 50
         if (.not.(opskip(j))) then
            nj = nj + 1
            nnn = 0
            do 30 i = mini , maxi
               if (.not.(opskip(i))) then
                  nnn = nnn + 1
                  pout(nnn) = prop(nnn+nstart,nj)/fac
               end if
 30         continue
            write (iw,6020) pnames(j) , (pout(i),i=1,nnn)
         end if
 40   continue
      nstart = nstart + nnn
      return
 6010 format (//3x,'perturbed',25x,'perturbation'/3x,'operator',4x,
     +        6(4x,a8,5x))
 6020 format (/2x,a8,2x,6f17.8)
      end
      subroutine prpol0(prop,npdim)
c
      implicit real*8  (a-h,o-z)
c
c     printing routine for polarisabilities
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      character *1 tag
c
      character *8 pnames
      common/tdhfx/pnames(50)
c
c
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      dimension prop(npdim,npdim)
      dimension alpha(6,6),tag(3)
c
      data zero,three,con1/0.0d0,3.0d0,0.16487762d0/
      data tag/'x','y','z'/
c
      call secget(isect(52),52,iblkj)
      iblkj = iblkj + 2 + lensec(npole) + lensec(npole*np)
     +        + ifreq*lensec(npdim*npdim)
      call wrt3(prop,npdim*npdim,iblkj,ifild)
      if (ifreq.ne.0 .and. freq(ifreq).ne.0.0d0) then
         write (iwr,6010) freq(ifreq)
         call prfreq(2,freq(ifreq))
         if (opunch(5)) write (ipu,6020) freq(ifreq)
      else
         if (oprn(26) .and. oprn(40)) write (iwr,6030)
         if (opunch(5)) write (ipu,6040)
      end if
      if (opunch(5)) write (ipu,6050) ((prop(i,j),i=1,np),j=1,np)
      if (ogen) then
         nstart = 0
         fac = 1.0d0
         mini = 1
         maxi = 6
 20      call prpol1(prop,npdim,mini,maxi,nstart,fac,iwr)
         mini = mini + 6
         if (mini.gt.np) return
         maxi = maxi + 6
         if (maxi.gt.np) maxi = np
         go to 20
      else
c
c     construct dipole-polarizability
c
         do 40 i = 1 , 3
            do 30 j = 1 , 3
               alpha(i,j) = zero
 30         continue
 40      continue
         npd = 0
         np1 = 0
         do 60 i = 1 , 3
            if (.not.(opskip(i))) then
               np1 = np1 + 1
               npd = npd + 1
               np2 = 0
               do 50 j = 1 , i
                  if (.not.(opskip(j))) then
                     np2 = np2 + 1
                     alpha(i,j) = prop(np1,np2)
                     alpha(j,i) = prop(np2,np1)
                  end if
 50            continue
            end if
 60      continue
         if (np1.ne.0) then
            if (oprn(26) .and. oprn(40)) write (iwr,6070) (tag(i),i=1,3)
            do 70 i = 1 , 3
               if (oprn(26) .and. oprn(40)) 
     +         write (iwr,6060) tag(i) , (alpha(i,j),j=1,3)
 70         continue
c
c     convert polarizability to s.i. units
c
            do 90 i = 1 , 3
               do 80 j = 1 , 3
                  alpha(i,j) = con1*alpha(i,j)
 80            continue
 90         continue
            if (oprn(26) .and. oprn(40)) write (iwr,6080) (tag(i),i=1,3)
            do 100 i = 1 , 3
               if (oprn(26) .and. oprn(40)) 
     +         write (iwr,6060) tag(i) , (alpha(i,j),j=1,3)
 100        continue
         end if
c
c
c    construct the dipole-quadrupole polarizability (a tensor)
c
         do 120 i = 1 , 6
            do 110 j = 1 , 6
               alpha(i,j) = zero
 110        continue
 120     continue
         np1 = 0
         do 140 i = 1 , 3
            if (.not.(opskip(i))) then
               np1 = np1 + 1
               np2 = npd
               do 130 j = 1 , 6
                  if (.not.(opskip(j+3))) then
                     np2 = np2 + 1
                     alpha(i,j) = prop(np1,np2)
                  end if
 130           continue
            end if
 140     continue
         if (np2.eq.npd) return
         if (oprn(26) .and. oprn(40) )
     +       write (iwr,6090)
         do 150 i = 1 , 3
            if (oprn(26) .and. oprn(40) )
     +       write (iwr,6100) (alpha(i,j),j=1,6)
 150     continue
c
c     quadrupole-quadrupole polarizability (c tensor)
c
         do 170 i = 1 , 6
            do 160 j = 1 , 6
               alpha(i,j) = zero
 160        continue
 170     continue
         np1 = npd
         do 190 i = 1 , 6
            if (.not.(opskip(i+3))) then
               np1 = np1 + 1
               np2 = npd
               do 180 j = 1 , 6
                  if (.not.(opskip(j+3))) then
                     np2 = np2 + 1
                     alpha(i,j) = prop(np1,np2)/three
                  end if
 180           continue
            end if
 190     continue
         if (oprn(26) .and. oprn(40) )
     +       write (iwr,6110)
         do 200 i = 1 , 6
            if (oprn(26) .and. oprn(40) )
     +       write (iwr,6100) (alpha(i,j),j=1,6)
 200     continue
         return
      end if
 6010 format (//10x,'****************************************'/10x,
     +        'second order properties at omega squared =',f20.10/10x,
     +        '****************************************'/)
 6020 format (1x,'dynamic polarizabilty at omega squared=',f16.6)
 6030 format (//
     +   30x,'**********************************'/
     +   30x,'* static second order properties *'/
     +   30x,'**********************************'/)
 6040 format (1x,'static polarizability')
 6050 format (1x,3e20.12)
 6060 format (/10x,a1,3f15.7)
 6070 format (//10x,'========================='/
     +          10x,'= polarizability tensor ='/
     +          10x,'========================='//
     +          10x,'in atomic units (bohr**3)'//
     +              18x,a1,14x,a1,14x,a1)
 6080 format (//10x,
     +        'in s.i. units (10**-40 farad meter**2)'//
     +     18x,a1,14x,a1,14x,a1)
 6090 format (//10x,'==============================================='/
     +          10x,'= dipole-quadrupole polarizability (a tensor) ='/
     +          10x,'==============================================='//
     +          10x,'in atomic units (bohr**4)'//)
 6100 format (10x,6f15.7)
 6110 format (//
     +  10x,'==================================================='/
     +  10x,'= quadrupole-quadrupole polarizability (c tensor) ='/
     +  10x,'==================================================='//
     +  10x,'in atomic units (bohr**5)'//)
      end
      subroutine chfopo(a,maxa,mn,ifil,iblk,lblk,ifils,iblks,
     &   nocca,nsoc,nvirta)
      implicit real*8  (a-h,o-z)
c
c     writes out results from sorto
c
      dimension a(mn,*)
      common/craypk/labs(1360)
      common/blkin/g(510),nint
      common/maxlen/maxq
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      data m0/0/
      npass = ((mn*mn)/maxa) + 1
      ncol = maxa/mn
      if (ncol.gt.mn) ncol = mn
      nst = 1
      nfin = ncol
      call search(iblk,ifil)
      lblk = 0
      do 80 np = 1 , npass
         call vclr(a,1,mn*ncol)
         call search(iblks,ifils)
 20      call find(ifils)
         call get(g,nw)
         if (nw.eq.0) then
            nint = 0
            nsd1 = nocca*nsoc
            nsv = nsd1 + nocca*nvirta
            nsd1 = nsd1 + 1
c           mini = 1
            do 50 j = nst , nfin
               do 40 i = 1 , j
                  val = a(i,j-nst+1)
                  if (val.ne.0.0d0) then
                     if (i.lt.nsd1 .or. i.gt.nsv) val = val*0.5d0
                     nint = nint + 1
                     g(nint) = val
                     labs(nint+nint-1) = j
                     labs(nint+nint) = i
                     if (nint.ge.num2e) then
                     call pack(g(num2ep+1),lab1632,labs,numlabp)
                     call put(g,m511,ifil)
                     do 30 iiii = 1 , 680
                        labs(iiii) = 0
 30                  continue
                     nint = 0
                     lblk = lblk + 1
                     end if
                  end if
 40            continue
 50         continue
            if (nint.ne.0) then
               call pack(g(num2ep+1),lab1632,labs,numlabp)
               call put(g,m511,ifil)
               do 60 iiii = 1 , 680
                  labs(iiii) = 0
 60            continue
               lblk = lblk + 1
               nint = 0
            end if
            nst = nst + ncol
            nfin = nfin + ncol
            if (nfin.gt.mn) nfin = mn
         else
            call unpack(g(num2ep+1),lab1632,labs,numlabp)
            do 70 n = 1 , nint
               i = labs(n+n-1)
               j = labs(n+n)
               if (j.ge.nst .and. j.le.nfin) then
                  a(i,j-nst+1) = a(i,j-nst+1) + g(n)
               end if
 70         continue
            go to 20
         end if
 80   continue
      nint = 0
      call pack(g(num2ep+1),16,labs,numlabp)
      call put(g,m0,ifil)
      lblk = lblk + 1
      return
      end
      subroutine rhscl(eps,bx,by,bz,eval,sx,sy,sz)
      implicit real*8  (a-h,o-z)
c
c     r h s of chf equations ( nuclear displacements )
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/maxlen/maxq
      common/craypk/labs(1360)
      common/small/labout(1360)
      common/blkin/gin(510),nword
      common/out/gout(510),nw
      dimension eps(*),bx(*),by(*),bz(*),eval(*),sx(*),sy(*),sz(*)
      data smal/1.0d-20/
      data zero,one,half,four/0.0d0,1.0d0,0.5d0,4.0d0/
c
      hf_wght = CD_HF_exchange_weight()
c
c     sorts out scf eigenvalues and perturbation matrix elements
c
      do 20 iiii = 1 , 1360
         labout(iiii) = 0
 20   continue
      nplus1 = nocca + 1
      call secget(isect(9),9,isec9)
c
c     read in eigenvalues
c
      call rdedx(eval,lds(isect(9)),isec9,ifild)
      do 40 ii = 1 , nocca
         i1 = mapie(ii)
         do 30 jj = nplus1 , nsa4
            j1 = mapie(jj)
            it = (jj-nocca-1)*nocca + ii
            eps(it) = one/(eval(j1)-eval(i1))
 30      continue
 40   continue
c
c     sort out integrals involving one virtual m.o.
c
      i = 0
      j = 0
      k = 0
      l = 0
      ifili = nufile(1)
      ifilo = ifils
      ib1 = kblk(1)
      ib2 = nblk(1)
      nw = 0
      nat3 = nat*3
      iblll = lensec(mn)
      iblkb = iblks
      ib3 = iblkb + iblll*nat3
      call search(ib1,ifili)
      call wrt3z(iblkb,ifilo,ib3)
      call search(ib3,ifilo)
      do 70 ibl = ib1 , ib2
         call find(ifili)
         call get(gin,nn)
c
c     complete list of two-electron integrals coming
c     in from transformed mainfile.
c     those of form <aj/kl> going out
c     on scratchfile (ed7)
c
         if (nn.eq.0) go to 80
        call unpack(gin(num2e+1),lab816,labs,numlab)
         do 60 n = 1 , nword
            n4 = (n+n) + (n+n)

            i = labs(n4-2)
            j = labs(n4-3)
            k = labs(n4  )
            l = labs(n4-1)
            if (i.gt.nocca) then
               if (j.le.nocca .and. k.le.nocca .and. l.le.nocca) then
                  nw = nw + 1
                  gout(nw) = gin(n)
                  n4 = (nw+nw) + (nw+nw)
                  labout(n4-3) = i
                  labout(n4-2) = j
                  labout(n4-1) = k
                  labout(n4) = l
                  if (nw.ge.num2e) then
c
c
c     writing out a block
c
                     call pack(gout(num2e+1),lab816,labout,numlab)
                     call put(gout,m511,ifilo)
                     do 50 iiii = 1 , 1360
                        labout(iiii) = 0
 50                  continue
                     ib3 = ib3 + 1
                     nw = 0
                  end if
               end if
            end if
 60      continue
 70   continue
 80   if (nw.ne.0) then
         call pack(gout(num2e+1),lab816,labout,numlab)
         call put(gout,m511,ifilo)
         do 90 iiii = 1 , 1360
            labout(iiii) = 0
 90      continue
         ib3 = ib3 + 1
      end if
      ib4 = ib3 - 1
c
c
c     sorting completed
c
c
      ib3 = iblkb + iblll*nat3
c      derivatives of overlap matrix section 14 of fockfile
c
      ibs = iochf(14)
      lennew = iky(ncoorb+1)
      newblk = lensec(lennew)
      np = nat3
      i = 0
      j = 0
      k = 0
      l = 0
      do 130 n = 1 , nat
         do 100 jj = 1 , mn
            bx(jj) = zero
            by(jj) = zero
            bz(jj) = zero
 100     continue
         call rdedx(sx,lennew,ibs,ifockf)
         call reads(sy,lennew,ifockf)
         call reads(sz,lennew,ifockf)
         ibs = ibs + newblk*3
c
c     have just read overlap derivatives sx,sy,sz for
c     one atom.  now scan list of integrals <aj/kl>
c
         call search(ib3,ifilo)
         do 120 ibl = ib3 , ib4
            call find(ifilo)
            call get(gin,nn)
            call unpack(gin(num2e+1),lab816,labs,numlab)
            do 110 int = 1 , nword
               n4 = (int+int) + (int+int)
               i = labs(n4-3)
               j = labs(n4-2)
               k = labs(n4-1)
               l = labs(n4)
               gg = -gin(int)
               if (k.eq.l) gg = gg*half
               gg4 = gg*four
               ii = (i-nocca-1)*nocca
               nij = ii + j
               nik = ii + k
               nil = ii + l
               j1 = mapie(j)
               k1 = mapie(k)
               l1 = mapie(l)
               nkl = iky(k1) + l1
               njk = iky(j1) + k1
               njl = iky(j1) + l1
               if (k1.gt.j1) njk = iky(k1) + j1
               if (l1.gt.j1) njl = iky(l1) + j1
               bx(nij) = bx(nij) + gg4*sx(nkl)
               by(nij) = by(nij) + gg4*sy(nkl)
               bz(nij) = bz(nij) + gg4*sz(nkl)
               gg = gg*hf_wght
               bx(nik) = bx(nik) - gg*sx(njl)
               by(nik) = by(nik) - gg*sy(njl)
               bz(nik) = bz(nik) - gg*sz(njl)
               bx(nil) = bx(nil) - gg*sx(njk)
               by(nil) = by(nil) - gg*sy(njk)
               bz(nil) = bz(nil) - gg*sz(njk)
 110        continue
 120     continue
c
c     have just formed that part of r.h.s. involving
c     product of s' with <ij/kl>
c     store this on stratchfile
c
         call wrt3(bx,mn,iblkb,ifils)
         call wrt3s(by,mn,ifils)
         call wrt3s(bz,mn,ifils)
         iblkb = iblkb + iblll*3
 130  continue
c
c
c
      ibh = iochf(13)
      ibs = iochf(14)
      iblkb = iblks
      do 180 l = 1 , nat3
c
c     read perturbed fock matrix elements (integral contribution
c     only) from fockfile at section iochf(13)
c
         call rdedx(sx,lennew,ibh,ifockf)
c
c     and the overlap derivatives again
c
         call rdedx(sy,lennew,ibs,ifockf)
         call rdedx(bx,mn,iblkb,ifils)
         ibh = ibh + newblk
         ibs = ibs + newblk
         do 150 j = 1 , nocca
            j1 = mapie(j)
            do 140 i = nplus1 , nsa4
               i1 = mapie(i)
               it = iky(i1) + j1
               mt = (i-nocca-1)*nocca + j
               bx(mt) = bx(mt) + sx(it) - eval(j1)*sy(it)
 140        continue
 150     continue
         do 170 j = 1 , nocca
            do 160 i = nplus1 , nsa4
               mt = (i-nocca-1)*nocca + j
               if (dabs(bx(mt)).lt.smal) bx(mt) = zero
 160        continue
 170     continue
c
c
c     write complete r.h.s. to scratchfile
c
         call wrt3(bx,mn,iblkb,ifils)
         iblkb = iblkb + iblll
 180  continue
      call clredx
      return
      end
      subroutine rhsgvb(mapnr,bx,by,bz,epsx,epsy,epsz,eta,zeta,sx,
     1                        sy,sz,acore,alpha,beta,iblok,iblok2)
c
c      right hand side of general scf chf equations
c
      implicit real*8  (a-h,o-z)
      dimension epsx(*),bx(*),by(*),bz(*),epsy(*),sx(*),sy(*),sz(*)
      dimension epsz(*),eta(*),alpha(11,*),beta(11,*),
     1          zeta(nx,njk1)
      dimension mapnr(*)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
      logical acore
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/maxlen/maxq
      common/blkin/gin(510),nword
      common/craypk/labs(1360)
      integer xat,x
      common/mpshl/ns(maxorb)
c
      data half/0.5d0/
      data two / 2.0d0/
c
      ioff(i) = i*(i-1)/2
c
c     sorts out scf eigenvalues and perturbation matrix elements
c
c
      do 20 iiii = 1 , 1360
         labs(iiii) = 0
 20   continue
      iblkb = iblks
      call vclr(bx,1,mn*nat)
      call vclr(by,1,mn*nat)
      call vclr(bz,1,mn*nat)
c     mtype = 0
      ibh = iochf(17)
      ibs = iochf(14)
      iblll = lensec(mn)
      numsq = ncoorb*ncoorb
      newblk = lensec(numsq)
c
      if (.not.acore) call caserr(' increase core to allow for b ')
c
      do 50 xat = 1 , nat
         call rdedx(sx,numsq,ibh,ifockf)
         call reads(sy,numsq,ifockf)
         call reads(sz,numsq,ifockf)
c
         ibh = ibh + newblk*3
c
         ij = 0
         do 40 i = 1 , num
            do 30 j = 1 , i
               ij = ij + 1
               if (mapnr(ij).gt.0) then
                  mp = mapnr(ij) + (xat-1)*mn
                  bx(mp) = -two*(sx(i+(j-1)*num)-sx(j+(i-1)*num))
                  by(mp) = -two*(sy(i+(j-1)*num)-sy(j+(i-1)*num))
                  bz(mp) = -two*(sz(i+(j-1)*num)-sz(j+(i-1)*num))
               end if
 30         continue
 40      continue
 50   continue
c
      newblk = lensec(nx)
      do 60 xat = 1 , nat
         mx = 1 + (xat-1)*nx
         call rdedx(sx(mx),nx,ibs,ifockf)
         call reads(sy(mx),nx,ifockf)
         call reads(sz(mx),nx,ifockf)
         ibs = ibs + newblk*3
 60   continue
c
c   for historical reasons the scratchfile is set up in the order
c
c   b : mnnr*nat*3
c
c   xir - xri : xir is the ci lagrangian over independent pairs : mnnr
c
c   zir : solution of atz = xir-xri : mnnr
c
c   eta : nx
c
c   zeta : nx*njk1
c
c     lenblk = lensec(mn)
      call rdedx(eta,nx,iblok,ifils)
c
      call rdedx(zeta,nx*njk1,iblok2,ifils)
c
      do 100 xat = 1 , nat
         ijx = 0
         do 90 i = 1 , ncoorb
            do 80 j = 1 , i
               ijx = ijx + 1
               if (mapnr(ijx).gt.0) then
                  mp = mapnr(ijx) + (xat-1)*mn
                  do 70 k = 1 , ncoorb
                     jkx = ioff(j) + k
                     ikx = ioff(i) + k
                     if (k.gt.j) jkx = ioff(k) + j
                     if (k.gt.i) ikx = ioff(k) + i
                     bx(mp) = bx(mp) + sx(jkx+(xat-1)*nx)
     +                        *(eta(ikx)-zeta(ikx,ns(j)))
     +                        - sx(ikx+(xat-1)*nx)
     +                        *(eta(jkx)-zeta(jkx,ns(i)))
                     by(mp) = by(mp) + sy(jkx+(xat-1)*nx)
     +                        *(eta(ikx)-zeta(ikx,ns(j)))
     +                        - sy(ikx+(xat-1)*nx)
     +                        *(eta(jkx)-zeta(jkx,ns(i)))
                     bz(mp) = bz(mp) + sz(jkx+(xat-1)*nx)
     +                        *(eta(ikx)-zeta(ikx,ns(j)))
     +                        - sz(ikx+(xat-1)*nx)
     +                        *(eta(jkx)-zeta(jkx,ns(i)))
 70               continue
               end if
 80         continue
 90      continue
 100  continue
      i = 0
      j = 0
      k = 0
      l = 0
      ifili = nufile(1)
      ib1 = kblk(1)
      ib2 = nblk(1)
c     nat3 = nat*3
      call search(ib1,ifili)
      do 200 ibl = ib1 , ib2
         call find(ifili)
         call get(gin,nw)
c
c     complete list of two-electron integrals coming
c     in from transformed mainfile.
c
         if (nw.ne.0) then
        call unpack(gin(num2e+1),lab816,labs,numlab)
            do 190 n = 1 , nword
               i = labs(4*n-2)
               j = labs(4*n-3)
               k = labs(4*n  )
               l = labs(4*n-1)
               ij = ioff(i) + j
               kl = ioff(k) + l
               ik = ioff(i) + k
               jl = ioff(j) + l
               if (l.gt.j) jl = ioff(l) + j
               il = ioff(i) + l
               jk = ioff(j) + k
               if (k.gt.j) jk = ioff(k) + j
c
c
c
c   type 6 : i = j = k = l
c
               if (ijkltp(i,j,k,l).ne.1) then
c
c
c   l virtual implies 3 or more virtual orbitals
c
                  if (ns(l).le.njk) then
c
c   type 2 : i = j = k > l
c
                     if (ijkltp(i,j,k,l).eq.2) then
                        if (mapnr(il).gt.0) then
                           mp = mapnr(il)
                           do 110 x = 1 , nat
                              bx(mp+(x-1)*mn) = bx(mp+(x-1)*mn)
     +                           + sx((x-1)*nx+ij)
     +                           *two*(alpha(ns(i),ns(i))
     +                           -alpha(ns(i),ns(l))+beta(ns(i),ns(i))
     +                           -beta(ns(i),ns(l)))*gin(n)
                              by(mp+(x-1)*mn) = by(mp+(x-1)*mn)
     +                           + sy((x-1)*nx+ij)
     +                           *two*(alpha(ns(i),ns(i))
     +                           -alpha(ns(i),ns(l))+beta(ns(i),ns(i))
     +                           -beta(ns(i),ns(l)))*gin(n)
                              bz(mp+(x-1)*mn) = bz(mp+(x-1)*mn)
     +                           + sz((x-1)*nx+ij)
     +                           *two*(alpha(ns(i),ns(i))
     +                           -alpha(ns(i),ns(l))+beta(ns(i),ns(i))
     +                           -beta(ns(i),ns(l)))*gin(n)
 110                       continue
                        end if
                        go to 190
                     end if
c
c
c   type 3 : i > j = k = l
c
                     if (ijkltp(i,j,k,l).eq.3) then
                        if (mapnr(il).gt.0) then
                           mp = mapnr(il)
                           do 120 x = 1 , nat
                              bx(mp+(x-1)*mn) = bx(mp+(x-1)*mn)
     +                           + sx((x-1)*nx+jk)
     +                           *two*(alpha(ns(i),ns(l))
     +                           -alpha(ns(l),ns(l))+beta(ns(i),ns(l))
     +                           -beta(ns(l),ns(l)))*gin(n)
                              by(mp+(x-1)*mn) = by(mp+(x-1)*mn)
     +                           + sy((x-1)*nx+jk)
     +                           *two*(alpha(ns(i),ns(l))
     +                           -alpha(ns(l),ns(l))+beta(ns(i),ns(l))
     +                           -beta(ns(l),ns(l)))*gin(n)
                              bz(mp+(x-1)*mn) = bz(mp+(x-1)*mn)
     +                           + sz((x-1)*nx+jk)
     +                           *two*(alpha(ns(i),ns(l))
     +                           -alpha(ns(l),ns(l))+beta(ns(i),ns(l))
     +                           -beta(ns(l),ns(l)))*gin(n)
 120                       continue
                        end if
                        go to 190
                     end if
c
c
c   type 4 : i = j > k = l
c
                     if (ijkltp(i,j,k,l).ge.4) then
                        gg = gin(n)
                        if (i.eq.j) gg = gg*half
                        if (k.eq.l) gg = gg*half
                        if (i.eq.k .and. j.eq.l) gg = gg*half
c
c   type 5 : i = k > j = l
c   type 6 : i = j > k > l
c   type 7 : i = k > j > l
c   type 8 : i > j = k > l
c   type 10 : i > j > k = l
c   type 11 : i > k = l > j
c   type 12 : i > j > k > l
c   type 13 : i > k > j > l
c   type 14 : i > k > l > j
c
                        if (mapnr(ij).gt.0) then
                           mp = mapnr(ij)
                           do 130 x = 1 , nat
                              bx(mp+(x-1)*mn) = bx(mp+(x-1)*mn)
     +                           + sx((x-1)*nx+kl)
     +                           *two*aijklx(i,j,k,l,ns,alpha)*gg
                              by(mp+(x-1)*mn) = by(mp+(x-1)*mn)
     +                           + sy((x-1)*nx+kl)
     +                           *two*aijklx(i,j,k,l,ns,alpha)*gg
                              bz(mp+(x-1)*mn) = bz(mp+(x-1)*mn)
     +                           + sz((x-1)*nx+kl)
     +                           *two*aijklx(i,j,k,l,ns,alpha)*gg
 130                       continue
                        end if
                        if (mapnr(kl).gt.0) then
                           mp = mapnr(kl)
                           do 140 x = 1 , nat
                              bx(mp+(x-1)*mn) = bx(mp+(x-1)*mn)
     +                           + sx((x-1)*nx+ij)
     +                           *two*aijklx(k,l,i,j,ns,alpha)*gg
                              by(mp+(x-1)*mn) = by(mp+(x-1)*mn)
     +                           + sy((x-1)*nx+ij)
     +                           *two*aijklx(k,l,i,j,ns,alpha)*gg
                              bz(mp+(x-1)*mn) = bz(mp+(x-1)*mn)
     +                           + sz((x-1)*nx+ij)
     +                           *two*aijklx(k,l,i,j,ns,alpha)*gg
 140                       continue
                        end if
                        if (mapnr(il).gt.0) then
                           mp = mapnr(il)
                           do 150 x = 1 , nat
                              bx(mp+(x-1)*mn) = bx(mp+(x-1)*mn)
     +                           + sx((x-1)*nx+jk)
     +                           *aijklx(i,l,j,k,ns,beta)*gg
                              by(mp+(x-1)*mn) = by(mp+(x-1)*mn)
     +                           + sy((x-1)*nx+jk)
     +                           *aijklx(i,l,j,k,ns,beta)*gg
                              bz(mp+(x-1)*mn) = bz(mp+(x-1)*mn)
     +                           + sz((x-1)*nx+jk)
     +                           *aijklx(i,l,j,k,ns,beta)*gg
 150                       continue
                        end if
                        if (mapnr(ik).gt.0) then
                           mp = mapnr(ik)
                           do 160 x = 1 , nat
                              bx(mp+(x-1)*mn) = bx(mp+(x-1)*mn)
     +                           + sx((x-1)*nx+jl)
     +                           *aijklx(i,k,j,l,ns,beta)*gg
                              by(mp+(x-1)*mn) = by(mp+(x-1)*mn)
     +                           + sy((x-1)*nx+jl)
     +                           *aijklx(i,k,j,l,ns,beta)*gg
                              bz(mp+(x-1)*mn) = bz(mp+(x-1)*mn)
     +                           + sz((x-1)*nx+jl)
     +                           *aijklx(i,k,j,l,ns,beta)*gg
 160                       continue
                        end if
                        if (mapnr(jl).gt.0) then
                           mp = mapnr(jl)
                           jsw = j
                           lsw = l
                           if (j.lt.l) then
                              jsw = l
                              lsw = j
                           end if
                           do 170 x = 1 , nat
                              bx(mp+(x-1)*mn) = bx(mp+(x-1)*mn)
     +                           + sx((x-1)*nx+ik)
     +                           *aijklx(jsw,lsw,i,k,ns,beta)*gg
                              by(mp+(x-1)*mn) = by(mp+(x-1)*mn)
     +                           + sy((x-1)*nx+ik)
     +                           *aijklx(jsw,lsw,i,k,ns,beta)*gg
                              bz(mp+(x-1)*mn) = bz(mp+(x-1)*mn)
     +                           + sz((x-1)*nx+ik)
     +                           *aijklx(jsw,lsw,i,k,ns,beta)*gg
 170                       continue
                        end if
                        if (mapnr(jk).gt.0) then
                           mp = mapnr(jk)
                           jsw = j
                           ksw = k
                           if (j.lt.k) then
                              jsw = k
                              ksw = j
                           end if
                           do 180 x = 1 , nat
                              bx(mp+(x-1)*mn) = bx(mp+(x-1)*mn)
     +                           + sx((x-1)*nx+il)
     +                           *aijklx(jsw,ksw,i,l,ns,beta)*gg
                              by(mp+(x-1)*mn) = by(mp+(x-1)*mn)
     +                           + sy((x-1)*nx+il)
     +                           *aijklx(jsw,ksw,i,l,ns,beta)*gg
                              bz(mp+(x-1)*mn) = bz(mp+(x-1)*mn)
     +                           + sz((x-1)*nx+il)
     +                           *aijklx(jsw,ksw,i,l,ns,beta)*gg
 180                       continue
                        end if
                     end if
                  end if
               end if
c     end if
c
c
c
 190        continue
         end if
c
 200  continue
c
      do 210 x = 1 , nat
         call wrt3(bx(1+(x-1)*mn),mn,iblkb,ifils)
         call wrt3s(by(1+(x-1)*mn),mn,ifils)
         call wrt3s(bz(1+(x-1)*mn),mn,ifils)
         iblkb = iblkb + iblll*3
 210  continue
c
      return
      end
      subroutine rhsrhf(eps,b,bb,u,eval,qq)
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c     r h s of chf equations ( nuclear displacements)
c
c     logical out
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/maxlen/maxq
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
      common/craypk/labs(1360)
      common/small/labout(1360)
      common/blkin/gin(510),nword
      common/out/gout(510),nw
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      dimension eps(*),b(*),bb(*),u(*),eval(*),qq(*)
      data smal/1.0d-20/
      data zero,one,two,half,four/0.0d0,1.0d0,2.0d0,0.5d0,4.0d0/
c
      ind(i,j) = iky(max(i,j)) + min(i,j)
c
c     sorts out scf eigenvalues and perturbation matrix elements
c
      do 20 iiii = 1 , 1360
         labout(iiii) = 0
 20   continue
c     out = oprn(12)
      nsoc = noccb - nocca
      nvirta = nsa4 - noccb
      nsd = nocca*nsoc
      mn = nocca*nsoc + noccb*nvirta
      ndpls1 = nocca + 1
      ntpls1 = noccb + 1
      call secget(isect(9),9,isec)
      call rdedx(eval,lds(isect(9)),isec,ifild)
      if (nocca.ne.0) then
         do 40 jj = 1 , nocca
            do 30 ii = ndpls1 , nsa4
               it = (ii-ndpls1)*nocca + jj
               eps(it) = one/(eval(mapie(ii))-eval(mapie(jj)))
 30         continue
 40      continue
      end if
      if (noccb.ne.nocca) then
         do 60 jj = ndpls1 , noccb
            do 50 ii = ntpls1 , nsa4
               it = nvirta*nocca + (ii-nsoc-1)*nsoc + jj - nocca
               eps(it) = one/(eval(mapie(ii))-eval(mapie(jj)))
 50         continue
 60      continue
      end if
c
c     sort out integrals contributing to b
c
      i = 0
      j = 0
      k = 0
      l = 0
      ifili = nufile(1)
      ifilo = ifils
      iblkb = iblks
      iblll = lensec(mn)
      nat3 = nat*3
      ib1 = kblk(1)
      ib2 = nblk(1)
      nw = 0
      ib3 = iblks + iblll*nat3
      call search(ib1,ifili)
      call wrt3z(iblkb,ifilo,ib3)
      call search(ib3,ifilo)
      do 90 ibl = ib1 , ib2
         call find(ifili)
         call get(gin,nn)
         if (nn.eq.0) go to 100
         if (nword.eq.0) go to 100
         call unpack(gin(num2e+1),lab816,labs,numlab)
         do 80 n = 1 , nword
            n4 = (n+n) + (n+n)
            i = labs(n4-2)
            j = labs(n4-3)
            k = labs(n4  )
            l = labs(n4-1)

            if (j.le.noccb .and. k.le.noccb .and. i.gt.nocca) then
               if (j.le.nocca .or. k.gt.nocca .or. i.gt.noccb) then
                  if (j.le.nocca .or. i.gt.noccb .or. l.le.nocca) then
                     nw = nw + 1
                     gout(nw) = gin(n)
                     n4 = (nw+nw) + (nw+nw)
                     labout(n4-3) = i
                     labout(n4-2) = j
                     labout(n4-1) = k
                     labout(n4) = l
                     if (nw.ge.num2e) then
                        call pack(gout(num2e+1),lab816,labout,numlab)
                        call put(gout,m511,ifilo)
                        do 70 iiii = 1 , 1360
                           labout(iiii) = 0
 70                     continue
                        ib3 = ib3 + 1
                        nw = 0
                     end if
                  end if
               end if
            end if
 80      continue
 90   continue
 100  if (nw.ne.0) then
         call pack(gout(num2e+1),lab816,labout,numlab)
         call put(gout,m511,ifilo)
         do 110 iiii = 1 , 1360
            labout(iiii) = 0
 110     continue
         ib3 = ib3 + 1
      end if
      ib4 = ib3 - 1
      ib3 = iblks + iblll*nat3
c
c
c     derivatives of fock matrix elements at section 13
c     derivatives of overlap at section 14
      ibh = iochf(13)
      ibs = iochf(14)
      lennew = iky(ncoorb+1)
      newblk = lensec(lennew)
      np = nat3
      i = 0
      j = 0
      k = 0
      l = 0
      do 150 n = 1 , nat
         do 120 jj = 1 , mn
            b(jj) = zero
            bb(jj) = zero
            u(jj) = zero
 120     continue
         call rdedx(qq,lennew,ibs,ifockf)
         ibs = ibs + newblk
         call rdedx(qq(lennew+1),lennew,ibs,ifockf)
         ibs = ibs + newblk
         call rdedx(qq(lennew+lennew+1),lennew,ibs,ifockf)
         ibs = ibs + newblk
         call actmot(qq(1),nsa4,mapie,iky)
         call actmot(qq(lennew+1),nsa4,mapie,iky)
         call actmot(qq(lennew+lennew+1),nsa4,mapie,iky)
         call search(ib3,ifilo)
         do 140 ibl = ib3 , ib4
            call find(ifilo)
            call get(gin,nn)
            call unpack(gin(num2e+1),lab816,labs,numlab)
            do 130 int = 1 , nword
               n4 = (int+int) + (int+int)
               i = labs(n4-3)
               j = labs(n4-2)
               k = labs(n4-1)
               l = labs(n4)
               gg = -gin(int)
               if (j.le.nocca) then
                  if (k.le.nocca) then
                     if (k.eq.l) gg = gg*half
                     ggd = gg*four
                     ggo = -gg
                     ii = (i-nocca-1)*nocca
                     nij = ii + j
                     nik = ii + k
                     nil = ii + l
                     nkl = iky(k) + l
                     njk = iky(j) + k
                     njl = iky(j) + l
                     if (k.gt.j) njk = iky(k) + j
                     if (l.gt.j) njl = iky(l) + j
                     b(nij) = b(nij) + ggd*qq(nkl)
                     bb(nij) = bb(nij) + ggd*qq(lennew+nkl)
                     u(nij) = u(nij) + ggd*qq(lennew+lennew+nkl)
                     b(nik) = b(nik) + ggo*qq(njl)
                     bb(nik) = bb(nik) + ggo*qq(lennew+njl)
                     u(nik) = u(nik) + ggo*qq(lennew+lennew+njl)
                     b(nil) = b(nil) + ggo*qq(njk)
                     bb(nil) = bb(nil) + ggo*qq(lennew+njk)
                     u(nil) = u(nil) + ggo*qq(lennew+lennew+njk)
                  else if (l.le.nocca) then
                     if (i.le.noccb) then
c  (sd/sd)
                        ggd = two*gg
                        ii = (i-nocca-1)*nocca
                        nij = ii + j
                        nkl = iky(k) + l
                        b(nij) = b(nij) + ggd*qq(nkl)
                        bb(nij) = bb(nij) + ggd*qq(lennew+nkl)
                        u(nij) = u(nij) + ggd*qq(lennew+lennew+nkl)
                        if (i.ne.k .or. j.ne.l) then
                           kk = (k-nocca-1)*nocca
                           nkl = kk + l
                           nij = ind(i,j)
                           b(nkl) = b(nkl) + ggd*qq(nij)
                           bb(nkl) = bb(nkl) + ggd*qq(lennew+nij)
                           u(nkl) = u(nkl) + ggd*qq(lennew+lennew+nij)
                        end if
                     else
c  (vd/sd)
                        ggo = -gg
                        ggd = two*gg
                        ggos = -half*gg
                        iis = (i-nocca-1)*nocca
                        ii = nvirta*nocca + (i-nsoc-1)*nsoc - nocca
                        nik = ii + k
                        nij = iis + j
                        nil = iis + l
                        njl = iky(j) + l
                        if (l.gt.j) njl = iky(l) + j
                        b(nik) = b(nik) + ggo*qq(njl)
                        bb(nik) = bb(nik) + ggo*qq(lennew+njl)
                        u(nik) = u(nik) + ggo*qq(lennew+lennew+njl)
                        b(nij) = b(nij) + ggd*qq(ind(k,l))
                        bb(nij) = bb(nij) + ggd*qq(lennew+ind(k,l))
                        u(nij) = u(nij) + ggd*qq(lennew+lennew+ind(k,l))
                        b(nil) = b(nil) + ggos*qq(ind(j,k))
                        bb(nil) = bb(nil) + ggos*qq(lennew+ind(k,j))
                        u(nil) = u(nil)
     +                           + ggos*qq(lennew+lennew+ind(k,j))
                     end if
                  else if (i.le.noccb) then
c  (sd/ss)
                     if (k.eq.l) gg = gg*half
                     ggd = gg*two
                     ii = (i-nocca-1)*nocca
                     nij = ii + j
                     nkl = iky(k) + l
                     b(nij) = b(nij) + ggd*qq(nkl)
                     bb(nij) = bb(nij) + ggd*qq(lennew+nkl)
                     u(nij) = u(nij) + ggd*qq(lennew+lennew+nkl)
                  else
c  (vd/ss)
                     ggos = -gg
                     if (k.eq.l) gg = gg*half
                     ggd = gg*two
                     ii = (i-nocca-1)*nocca
                     iis = nvirta*nocca + (i-nsoc-1)*nsoc - nocca
                     nik = iis + k
                     nil = iis + l
                     nij = ii + j
                     nkl = iky(k) + l
                     b(nij) = b(nij) + ggd*qq(nkl)
                     bb(nij) = bb(nij) + ggd*qq(lennew+nkl)
                     u(nij) = u(nij) + ggd*qq(lennew+lennew+nkl)
                     b(nik) = b(nik) + ggos*qq(ind(j,l))
                     bb(nik) = bb(nik) + ggos*qq(lennew+ind(j,l))
                     u(nik) = u(nik) + ggos*qq(lennew+lennew+ind(j,l))
                     if (k.ne.l) then
                        b(nil) = b(nil) + ggos*qq(ind(j,k))
                        bb(nil) = bb(nil) + ggos*qq(lennew+ind(j,k))
                        u(nil) = u(nil)
     +                           + ggos*qq(lennew+lennew+ind(j,k))
                     end if
                  end if
               else if (i.le.noccb) then
                  if (i.eq.j) gg = gg*half
                  ggd = gg*two
                  kk = (k-nocca-1)*nocca
                  nkl = kk + l
                  nij = iky(i) + j
                  b(nkl) = b(nkl) + ggd*qq(nij)
                  bb(nkl) = bb(nkl) + ggd*qq(lennew+nij)
                  u(nkl) = u(nkl) + ggd*qq(lennew+lennew+nij)
               else if (k.le.nocca) then
c  (vs/dd)
                  ggos = -half*gg
                  if (k.eq.l) gg = gg*half
                  ggd = gg*four
                  ii = nvirta*nocca + (i-nsoc-1)*nsoc - nocca
                  iis = (i-nocca-1)*nocca
                  nik = iis + k
                  nil = iis + l
                  nij = ii + j
                  nkl = iky(k) + l
                  b(nij) = b(nij) + ggd*qq(nkl)
                  bb(nij) = bb(nij) + ggd*qq(lennew+nkl)
                  u(nij) = u(nij) + ggd*qq(lennew+lennew+nkl)
                  b(nik) = b(nik) + ggos*qq(ind(j,l))
                  bb(nik) = bb(nik) + ggos*qq(lennew+ind(j,l))
                  u(nik) = u(nik) + ggos*qq(lennew+lennew+ind(j,l))
                  if (k.ne.l) then
                     b(nil) = b(nil) + ggos*qq(ind(k,j))
                     bb(nil) = bb(nil) + ggos*qq(lennew+ind(k,j))
                     u(nil) = u(nil) + ggos*qq(lennew+lennew+ind(k,j))
                  end if
               else if (l.le.nocca) then
c (vs/sd)
                  ggo = -half*gg
                  ggd = two*gg
                  ggos = -gg
                  iis = nvirta*nocca + (i-nsoc-1)*nsoc - nocca
                  ii = (i-nocca-1)*nocca
                  nil = ii + l
                  nij = iis + j
                  nik = iis + k
                  njk = iky(j) + k
                  if (k.gt.j) njk = iky(k) + j
                  b(nil) = b(nil) + ggo*qq(njk)
                  bb(nil) = bb(nil) + ggo*qq(lennew+njk)
                  u(nil) = u(nil) + ggo*qq(lennew+lennew+njk)
                  b(nij) = b(nij) + ggd*qq(ind(k,l))
                  bb(nij) = bb(nij) + ggd*qq(lennew+ind(k,l))
                  u(nij) = u(nij) + ggd*qq(lennew+lennew+ind(k,l))
                  b(nik) = b(nik) + ggos*qq(ind(j,l))
                  bb(nik) = bb(nik) + ggos*qq(lennew+ind(j,l))
                  u(nik) = u(nik) + ggos*qq(lennew+lennew+ind(j,l))
               else
                  if (k.eq.l) gg = gg*half
                  ggd = gg*two
                  ggo = -gg
                  ii = nvirta*nocca + (i-nsoc-1)*nsoc - nocca
                  nij = ii + j
                  nik = ii + k
                  nil = ii + l
                  nkl = iky(k) + l
                  njk = iky(j) + k
                  njl = iky(j) + l
                  if (k.gt.j) njk = iky(k) + j
                  if (l.gt.j) njl = iky(l) + j
                  b(nij) = b(nij) + ggd*qq(nkl)
                  bb(nij) = bb(nij) + ggd*qq(lennew+nkl)
                  u(nij) = u(nij) + ggd*qq(lennew+lennew+nkl)
                  b(nik) = b(nik) + ggo*qq(njl)
                  bb(nik) = bb(nik) + ggo*qq(lennew+njl)
                  u(nik) = u(nik) + ggo*qq(lennew+lennew+njl)
                  b(nil) = b(nil) + ggo*qq(njk)
                  bb(nil) = bb(nil) + ggo*qq(lennew+njk)
                  u(nil) = u(nil) + ggo*qq(lennew+lennew+njk)
               end if
 130        continue
 140     continue
         call wrt3(b,mn,iblkb,ifils)
         iblkb = iblkb + iblll
         call wrt3(bb,mn,iblkb,ifils)
         iblkb = iblkb + iblll
         call wrt3(u,mn,iblkb,ifils)
         iblkb = iblkb + iblll
 150  continue
c
c
      ibs = iochf(14)
      iblkb = iblks
      ioffk = lennew + lennew + 1
      ioffk1 = ioffk + lennew
      ioffk0 = ioffk - 1
      ioffka = ioffk1 - 1
      call secget(isect(43),43,iblok)
      call rdedx(qq(ioffk),lds(isect(43)),iblok,ifild)
      lenblk = lensec(nx)
      iblkk = iochf(13) + nat3*lenblk
      do 260 l = 1 , nat3
         call rdedx(qq,lennew,ibh,ifockf)
         call rdedx(qq(lennew+1),lennew,ibs,ifockf)
         call rdedx(b,mn,iblkb,ifils)
         call rdedx(qq(ioffk1),lennew,iblkk,ifockf)
         call actmot(qq(1),nsa4,mapie,iky)
         call actmot(qq(lennew+1),nsa4,mapie,iky)
         call actmot(qq(ioffk1),nsa4,mapie,iky)
         iblkk = iblkk + lenblk
         ibh = ibh + newblk
         ibs = ibs + newblk
         if (nocca.ne.0) then
            do 190 j = 1 , nocca
               do 180 i = ndpls1 , nsa4
                  xxk = 1.0d0
                  if (i.gt.noccb) xxk = 2.0d0
                  it = iky(i) + j
                  mt = (i-nocca-1)*nocca + j
                  b(mt) = b(mt) + qq(it) - eval(j)*qq(lennew+it)
                  if (i.le.noccb) b(mt) = b(mt) + qq(ioffka+ind(i,j))
                  do 160 k = 1 , nocca
                     ik = iky(i) + k
                     b(mt) = b(mt) + xxk*qq(ind(j,k)+ioffk0)
     +                       *qq(lennew+ik)
 160              continue
                  if (nocca.ne.noccb) then
                     if (i.gt.noccb) then
                        do 170 k = ndpls1 , noccb
                           ik = iky(i) + k
                           b(mt) = b(mt) + qq(ind(k,j)+ioffk0)
     +                             *qq(lennew+ik)
 170                    continue
                     end if
                  end if
 180           continue
 190        continue
         end if
         if (nocca.ne.noccb) then
            do 230 j = ndpls1 , noccb
               do 220 i = ntpls1 , nsa4
                  it = iky(i) + j
                  mt = nvirta*nocca + (i-nsoc-1)*nsoc + j - nocca
                  b(mt) = b(mt) - qq(ind(i,j)+ioffka) + qq(it) - eval(j)
     +                    *qq(lennew+it)
                  do 200 k = 1 , nocca
                     ik = iky(i) + k
                     jk = iky(j) + k
                     b(mt) = b(mt) + 2.0d0*qq(ind(j,k)+ioffk0)
     +                       *qq(lennew+ik) + qq(ind(k,i)+ioffk0)
     +                       *qq(lennew+jk)
 200              continue
                  do 210 k = ndpls1 , noccb
                     ik = iky(i) + k
                     b(mt) = b(mt) + qq(ind(j,k)+ioffk0)*qq(lennew+ik)
 210              continue
 220           continue
 230        continue
         end if
c
c
c
         nsv = nsd + nocca*nvirta
         nsd1 = nsd + 1
         do 240 i = nsd1 , nsv
            b(i) = b(i) + b(i)
 240     continue
         do 250 i = 1 , mn
            b(i) = b(i)*0.5d0
            if (dabs(b(i)).lt.smal) b(i) = zero
 250     continue
         call wrt3(b,mn,iblkb,ifils)
         iblkb = iblkb + iblll
 260  continue
      do 270 i = nsd1 , nsv
         eps(i) = eps(i)*0.5d0
 270  continue
      do 280 i = 1 , mn
         eps(i) = eps(i) + eps(i)
 280  continue
      return
      end
      subroutine rhsemp(eps,b,x,mnr)
c
c   constructs right-hand-side of chf for electric and
c   and magnetic perturbations (called RHSIN in CADPAC)
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension eps(*),mnr(*),b(*),x(*)
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/mpshl/ns(maxorb)
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
      character *8 pnames
      common/tdhfx/pnames(50)
c
c
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      character*10 charwall
      character *8 grhf
      data grhf/'grhf'/
      data maxper/50/
c
      data zero,one/0.0d0,1.0d0/
c
c     sorts out scf eigenvalues and perturbation matrix elements
c
c     read in orbital energies
c
      iblkb = iblks
      if (scftyp.ne.grhf) then
         call secget(isect(9),9,iblok)
         call rdedx(x,lds(isect(9)),iblok,ifild)
c
         nsoc = noccb - nocca
         ntpls1 = noccb + 1
         ndpls1 = nocca + 1
         if (nocca.ne.0) then
            do 30 j = 1 , nocca
               do 20 i = ndpls1 , nsa4
                  it = (i-ndpls1)*nocca + j
                  eps(it) = one/(x(mapie(i))-x(mapie(j)))
 20            continue
 30         continue
         end if
         if (noccb.ne.nocca) then
            do 50 j = ndpls1 , noccb
               do 40 i = ntpls1 , nsa4
                  it = nvirta*nocca + (i-nsoc-1)*nsoc + j - nocca
                  eps(it) = one/(x(mapie(i))-x(mapie(j)))
 40            continue
 50         continue
            nsd = nocca*nsoc
            nsv = nsd + nocca*nvirta
            nsd1 = nsd + 1
            do 60 i = nsd1 , nsv
               eps(i) = eps(i)*0.5d0
 60         continue
            do 70 i = 1 , mn
               eps(i) = eps(i) + eps(i)
 70         continue
         end if
      end if
      iblll = lensec(mn)
      lennew = ikyp(nsa4)
c     nword = ikyp(ncoorb)
      np = 0
      do 110 l = 1 , maxper
         if (.not.(opskip(l))) then
            do 80 i = 1 , mn
               b(i) = zero
 80         continue
            if(odebug(30)) write (iwr,6010) pnames(l) , ipsec(l)
            np = np + 1
            jtype = 0
            call secget(ipsec(l),jtype,isec)
            call rdedx(x,lds(ipsec(l)),isec,ifild)
            call ijconr(x,b,lennew,mnr)
            ij = 0
            do 100 i = 1 , nsa4
               do 90 j = 1 , i
                  ij = ij + 1
                  if (ns(i).ne.ns(j)) then
                     nr = mnr(ij)
                     b(nr) = b(nr)*(fjk(ns(j))-fjk(ns(i)))*0.5d0
                  end if
 90            continue
 100        continue
            call wrt3(b,mn,iblkb,ifils)
            iblkb = iblkb + iblll
         end if
 110  continue
      write (iwr,6020) cpulft(1) ,charwall()
      return
 6010 format (1x,'perturbation',4x,a4,' from section',i6)
 6020 format (/1x,'construction of r.h.s. of equations',' complete at',
     +        1x,f8.2,' seconds',a10,' wall')
      end
      subroutine ovlcl(qq,iqq)
c
c    assemble overlap contribution to closed shell scf second
c    derivatives
c   term involving derivative of lagrangian
c   closed shell case
c
      implicit real*8  (a-h,o-z)
      logical out
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      common/bufb/e(maxorb),e1(maxorb)
      dimension qq(*),iqq(*)
      common/maxlen/maxq
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      out = odebug(6)
      call secget(isect(9),9,isec9)
      call rdedx(e,lds(isect(9)),isec9,ifild)
      nat3 = nat*3
      nlen = nat3*nat3
      iof = lenrel(nw196(5))
      iofs = nw196(5) + lenint(nat*nt)
      i10 = iofs + 1
      i11 = i10 + nlen
      ioff = i11 + nx
      noc1 = na + 1
c
c     perturbed fock matrix at section 16 , m.o. basis
c     derivatives of overlap at section 14
c
      ibf = iochf(16)
      ibs = iochf(14)
      ltri = ikyp(ncoorb)
      length = lensec(ltri)
      if (odebug(6)) then
         write (iwr,6010)
         do 30 n = 1 , nat3
            call rdedx(qq(i11),ltri,ibf,ifockf)
            call rdedx(qq(ioff),ltri,ibs,ifockf)
            ibf = ibf + length
            ibs = ibs + length
            do 20 i = 1 , ncoorb
               ii = ikyp(i) - 1
c
c     perturbed eigenvalues ( used as check only)
c
               e1(i) = qq(i11+ii) - e(i)*qq(ioff+ii)
 20         continue
            write (iwr,6020) (e1(i),i=1,ncoorb)
 30      continue
      end if
      ibf = iochf(16)
c
c     perturbed density matrix at section 15
c
      ibd = iochf(15)
      do 90 n = 1 , nat3
         n10 = i10 - 1 + n
c
c     complete derivative of fock matrix (only elements
c     with two occupied orbitals are needed)
c
         call rdedx(qq(i11),ltri,ibf,ifockf)
c
c     perturbed density matrix
c
         call rdedx(qq(ioff),ltri,ibd,ifockf)
c
c     derivative lagrangian elements for two occupied orbitals
c
         do 50 i = 1 , na
            ii = iky(i) - 1
            do 40 j = 1 , i
               ij = ii + j
               qq(ioff+ij) = qq(ioff+ij)*(e(i)+e(j)) + qq(i11+ij)
     +                       + qq(i11+ij)
 40         continue
 50      continue
c
c     derivative lagrangian elements for one occupied and one virtual m.
c
         do 70 i = noc1 , ncoorb
            ii = iky(i) - 1
            do 60 j = 1 , na
               ij = ii + j
               qq(ioff+ij) = qq(ioff+ij)*e(j)
 60         continue
 70      continue
c
c
         ibf = ibf + length
         ibd = ibd + length
         ibs = iochf(14)
         do 80 m = 1 , nat3
c
c     derivative overlap matrix
c
            call rdedx(qq(i11),ltri,ibs,ifockf)
            ibs = ibs + length
c
c     take product with overlap derivatives
c
            qq(n10+(m-1)*nat3) = -tracep(qq(i11),qq(ioff),ncoorb)
 80      continue
 90   continue
      call rdedx(qq(1),nw196(5),ibl196(5),ifild)
      if (out) then
         call dr2sym(qq(i10),qq(i11),iqq(1),iqq(iof+1),nat,nat3,
     +               nshell)
         write (iwr,6030)
         call prnder(qq(i10),nat3,iwr)
      end if
      call secget(isect(60),60,isec46)
      call rdedx(qq(i11),nlen,isec46,ifild)
      call vadd(qq(i11),1,qq(i10),1,qq(i11),1,nlen)
      call wrt3(qq(i11),nlen,isec46,ifild)
      return
 6010 format (//5x,'perturbed eigenvalues'//)
 6020 format (//(5x,6f16.8))
 6030 format (//' contribution from derivative of lagrangian')
      end
      subroutine ovlop(qq,iqq)
c
c    assemble overlap contribution to high-spin open-shell
c    scf second derivatives
c
c   term involving derivative of lagrangian
c   open shell case
c
      implicit real*8  (a-h,o-z)
      logical out
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      common/bufb/e(maxorb),e1(maxorb)
      dimension qq(*),iqq(*)
      common/maxlen/maxq
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      out = odebug(6)
      call secget(isect(9),9,isec9)
      call rdedx(e,lds(isect(9)),isec9,ifild)
      nat3 = nat*3
      nlen = nat3*nat3
      iof = lenrel(nw196(5))
      iofs = nw196(5) + lenint(nat*nt)
      i10 = iofs + 1
      i11 = i10 + nlen
      ioff = i11 + nx
      ioff2 = ioff + nx
      ioff3 = ioff2 + nx
      ioff4 = ioff3 + nx
      ioff5 = ioff4 + nx
      ioff6 = ioff5 + nx
      ioff7 = ioff6 + nx
      if (maxq.lt.(ioff7+nx+510+iofs)) call caserr(' insufficient core')
      ntpls1 = nb + 1
c     ndpls1 = na + 1
c     noc1 = nocca + 1
c
c     perturbed fock matrix at section 16 , m.o. basis
c     derivatives of overlap at section 14
c
      ibf = iochf(16)
      ibs = iochf(14)
      ltri = ikyp(ncoorb)
      length = lensec(ltri)
      if (odebug(6)) then
         write (iwr,6010)
         do 30 n = 1 , nat3
            call rdedx(qq(i11),ltri,ibf,ifockf)
            call rdedx(qq(ioff),ltri,ibs,ifockf)
            ibf = ibf + length
            ibs = ibs + length
            do 20 i = 1 , ncoorb
               ii = ikyp(i) - 1
c
c     perturbed eigenvalues ( used as check only)
c
               e1(i) = qq(i11+ii) - e(i)*qq(ioff+ii)
 20         continue
            write (iwr,6020) (e1(i),i=1,ncoorb)
 30      continue
      end if
      ibf = iochf(16)
      ibss = iochf(14)
      ibd = iochf(15)
c  part of perturbed fock matrix as calculated by fd2 at section 16
c  overlap derivatives at section 14
c  perturbed density matrices at section 15
      call secget(isect(43),43,iblok)
      call rdedx(qq(ioff3),ltri,iblok,ifild)
c  read in 1/2k
      do 90 n = 1 , nat3
         n10 = i10 - 1 + n
         call rdedx(qq(nlen),ltri,ibf,ifockf)
         call rdedx(qq(ioff),ltri,ibd,ifockf)
         call rdedx(qq(ioff2),ltri,ibss,ifockf)
c  form dk and kd
         call mxmtri(qq(ioff4),qq(ioff),qq(ioff3),ltri,ncoorb)
         call mxmtri(qq(ioff5),qq(ioff3),qq(ioff),ltri,ncoorb)
c  form kfs and sfk
         call mxmftr(qq(ioff6),qq(ioff3),qq(ioff2),ltri,ncoorb,na,
     +  nb)
         call mxmftr(qq(ioff7),qq(ioff2),qq(ioff3),ltri,ncoorb,na,
     +  nb)
c  two occupied orbitals
         do 50 i = 1 , nb
            foci = 2.0d0
            if (i.gt.na) foci = 1.0d0
            ii = iky(i) - 1
            do 40 j = 1 , i
               focj = 2.0d0
               if (j.gt.na) focj = 1.0d0
               ij = ii + j
               qq(ioff+ij) = qq(ioff+ij)*(e(i)+e(j))
     +                       /2.0d0 - qq(ioff2+ij)*(foci*e(i)+focj*e(j))
     +                       /2.0d0 + (foci-focj)
     +                       *(qq(ioff4+ij)-qq(ioff5+ij))
     +                       /2.0d0 + (foci+focj)
     +                       *(qq(ioff6+ij)+qq(ioff7+ij))
     +                       /2.0d0 + qq(i11+ij)
 40         continue
 50      continue
c  one occupied and one virtual m.o.
         do 70 i = ntpls1 , ncoorb
            ii = iky(i) - 1
            do 60 j = 1 , nb
               focj = 2.0d0
               if (j.gt.na) focj = 1.0d0
               ij = ii + j
               qq(ioff+ij) = qq(ioff+ij)*e(j) - focj*qq(ioff4+ij)
 60         continue
 70      continue
         ibf = ibf + length
         ibd = ibd + length
         ibss = ibss + length
         ibsb = iochf(14)
         do 80 m = 1 , nat3
            call rdedx(qq(i11),ltri,ibsb,ifockf)
            ibsb = ibsb + length
c
c     take product with overlap derivatives
c
            qq(n10+(m-1)*nat3) = -tracep(qq(i11),qq(ioff),ncoorb)
 80      continue
 90   continue
      if (out) then
         call rdedx(qq(1),nw196(5),ibl196(5),ifild)
         call dr2sym(qq(i10),qq(i11),iqq(1),iqq(iof+1),nat,nat3,
     +               nshell)
         write (iwr,6030)
         call prnder(qq(i10),nat3,iwr)
      end if
      call secget(isect(60),60,isec46)
      call rdedx(qq(i11),nlen,isec46,ifild)
      call vadd(qq(i11),1,qq(i10),1,qq(i11),1,nlen)
      call wrt3(qq(i11),nlen,isec46,ifild)
      return
 6010 format (//5x,'perturbed eigenvalues'//)
 6020 format (//(5x,6f16.8))
 6030 format (//' contribution from derivative of lagrangian')
      end
      subroutine sgmatm(fock,p,nfok,ltri)
c
c     closed shell fock operator construction
c     vectorised version to produce
c     many fock operators simulataneously
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension p(*),fock(*)
      common/blkin/gg(510),mword
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
      logical ohf_exch
      common/craypk/labs(1360)
c
      iword = 1
      call unpack(gg(num2e+1),lab816,labs,numlab)
c
      hf_wght  = CD_HF_exchange_weight()
      ohf_exch = (.not.CD_active()).or.CD_HF_exchange()
      do 40 iw = 1 , mword

         i = labs(iword+1)
         j = labs(iword  )
         k = labs(iword+3)
         l = labs(iword+2)
         gik = gg(iw)
         g2 = gik + gik
         g4 = g2 + g2
         ikyi = iky(i)
         ikyj = iky(j)
         ikyk = iky(k)
         ik = ikyi + k
         il = ikyi + l
         ij = ikyi + j
         jk = ikyj + k
         jl = ikyj + l
         kl = ikyk + l
         ioff = 0
         do 20 n = 1 , nfok
            aij = g4*p(kl+ioff) + fock(ij+ioff)
            fock(kl+ioff) = g4*p(ij+ioff) + fock(kl+ioff)
            fock(ij+ioff) = aij
            ioff = ioff + ltri
 20      continue
c... exchange
         if (ohf_exch) then
            g2  = hf_wght*g2
            gik = hf_wght*gik
         else
            iword = iword+4
            goto 40
         endif
         gil = gik
         if (i.eq.k .or. j.eq.l) gik = g2
         if (j.eq.k) gil = g2
         if (j.lt.k) then
            jk = ikyk + j
            if (j.lt.l) then
               jl = iky(l) + j
            end if
         end if
         ioff = 0
         do 30 n = 1 , nfok
            ajk = fock(jk+ioff) - gil*p(il+ioff)
            ail = fock(il+ioff) - gil*p(jk+ioff)
            aik = fock(ik+ioff) - gik*p(jl+ioff)
            fock(jl+ioff) = fock(jl+ioff) - gik*p(ik+ioff)
            fock(jk+ioff) = ajk
            fock(il+ioff) = ail
            fock(ik+ioff) = aik
            ioff = ioff + ltri
 30      continue
         iword = iword + 4
 40   continue
      return
      end
      subroutine chfeqs(eps,b,cc,wks1,wks2,alpha,aa,ndim,skipp)
      implicit real*8  (a-h,o-z)
c
c     simultaneous equations - small case
c
      logical skipp
      dimension skipp(100)
c
      dimension eps(ndim),b(ndim),cc(ndim),wks1(ndim),wks2(ndim),
     &    alpha(ndim,ndim),aa(ndim,ndim)
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
      common/craypk/labs(1360)
      common/blkin/g(510),nword
c
      data zero,one/0.0d0,1.0d0/
c
c
      lenblk = lensec(mn)
c
c     rhs of equations input from scratchfile
c
      iblkb = iblks
      iblku = iblkb + lenblk*np
      idev4 = nofile(1)
      iblk4 = jblk(1)
      if (oprn(12)) write (iwr,6010)
      if (oprn(13)) write (iwr,6020)
      ifail = 0
c
      do 30 i = 1 , mn
         do 20 j = 1 , mn
            alpha(j,i) = zero
 20      continue
 30   continue
      do 40 i = 1 , mn
         alpha(i,i) = one/eps(i)
 40   continue
c
c     read thru 2-electron integrals ( a-matrix )
c     file idev4 = nofile(1) = ed4 (default)
c     iblk4 =starting block
c
      call search(iblk4,idev4)
      call find(idev4)
 50   call get(g(1),nw)
      if (nw.gt.0) then
         if (nword.gt.0) then
c
c     use a block of integrals
c
            call find(idev4)
            call unpack(g(num2ep+1),lab1632,labs,numlabp)
            do 60 i = 1 , nword
               lab1 = labs(i+i-1)
               lab2 = labs(i+i)
               gg = g(i)
               alpha(lab1,lab2) = alpha(lab1,lab2) + gg
               if (.not.(lcpf .or. cicv .or. cicx .or. lab1.eq.lab2))
     +             then
                  alpha(lab2,lab1) = alpha(lab2,lab1) + gg
               end if
 60         continue
            go to 50
         end if
      end if
c
c     loop over the perturbations
c
      if (odebug(2)) write (iwr,6030)
      if (odebug(2)) call prsqm(alpha,mn,mn,mn,iwr)
      do 80 i = 1 , mn
         do 70 j = 1 , mn
            alpha(i,j) = alpha(i,j)*eps(i)
 70      continue
 80   continue
c
      do 110 ko = 1 , np
c     get rhs
c
         call rdedx(b,mn,iblkb,ifils)
c      write(iwr,978)(b(i),i=1,mn)
c978   format(' b in chfeqs = ',5f15.10)
         if (oprn(12)) write (iwr,6040) ko , (b(i),i=1,mn)
         iblkb = iblkb + lenblk
         do 90 i = 1 , mn
            cc(i) = 0.0d0
 90      continue
         if (odebug(2).and.skipp(ko)) write (iwr,6050) ko
         if (.not.(skipp(ko))) then
            do 100 i = 1 , mn
               b(i) = -b(i)*eps(i)
 100        continue
c
c    use nag routine to solve
c
            call f04atf(alpha,mn,b,mn,cc,aa,mn,wks1,wks2,ifail)
         end if
c
c    write solution onto scratchfile
c
         call wrt3(cc,mn,iblku,ifils)
         iblku = iblku + lenblk
         if (oprn(12) .or. oprn(13)) then
            write (iwr,6060) ko , (cc(i),i=1,mn)
         end if
 110  continue
      return
 6010 format (//1x,'print right-hand-side of chf equations')
 6020 format (//1x,'print solution to chf equations')
 6030 format (//1x,'a-matrix routine chfeqs')
 6040 format (//1x,'perturbation  ',i4//(5x,5f16.8))
 6050 format (/1x,'perturbation',i5,' omitted')
 6060 format (//1x,'solution  ',i4//(5x,5f16.8))
      end
      subroutine mxmtri(a,b,c,ltri,n)
c
c  produces lower triangle of the product matrix a=bc
c  where b and c are symmetric matrices input as lower triangles
c  n is the dimension of the square matrices,
c  ltri is number of elements in lower triangle.
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension a(ltri),b(ltri),c(ltri)
      do 20 i = 1 , ltri
         a(i) = 0.0d0
 20   continue
      do 50 i = 1 , n
         ii = iky(i)
         do 40 k = 1 , n
            ik = ii + k
            if (k.gt.i) ik = iky(k) + i
            bik = b(ik)
            if (bik.ne.0.0d0) then
               kk = iky(k)
               do 30 j = 1 , i
                  jk = kk + j
                  if (j.gt.k) jk = iky(j) + k
                  a(ii+j) = a(ii+j) + bik*c(jk)
 30            continue
            end if
 40      continue
 50   continue
      return
      end
      subroutine mxmftr(a,b,c,ltri,n,ndoc,ntot)
c
c  produces lower triangle of the product matrix a=bfc
c  where b and c are symmetric matrices input as lower triangles
c  and f is the occupation number of the kth orbital (f=0,1,2)
c  n is the dimension of the square matrices,
c  ltri is number of elements in lower triangle.
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension a(ltri),b(ltri),c(ltri)
      do 20 i = 1 , ltri
         a(i) = 0.0d0
 20   continue
      do 50 i = 1 , n
         ii = iky(i)
         do 40 k = 1 , ntot
            fock = 2.0d0
            if (k.gt.ndoc) fock = 1.0d0
            ik = ii + k
            if (k.gt.i) ik = iky(k) + i
            bfik = b(ik)*fock
            if (bfik.ne.0.0d0) then
               kk = iky(k)
               do 30 j = 1 , i
                  jk = kk + j
                  if (j.gt.k) jk = iky(j) + k
                  a(ii+j) = a(ii+j) + bfik*c(jk)
 30            continue
            end if
 40      continue
 50   continue
      return
      end
      subroutine chfcls(a,maxa)
c
c   sorts out a-matrix ( hessian ) for closed shell chf
c   ------------------------------------------------------
c
      implicit real*8  (a-h,o-z)
      dimension a(*)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
      common/craypk/labs(1360)
      common/blkin/g(510),nword
      character*10 charwall
c
      hf_wght = CD_HF_exchange_weight()
      call search(jblk(1),nofile(1))
c
c     work out number of passes
c     a matrix is mn*(mn+1)/2 ; core available is maxa
c
      call setsto(1360,0,labs)
      if (.not.lcpf .and. .not.cicv .and. (.not.cicx)) mnnr = mn
      nst = 0
      nmin = 0
      nfin = mnnr
      do 20 nn = 1 , mnnr
         last = nn*(nn+1)/2
         if (last.gt.maxa) then
            nfin = nn - 1
            go to 30
         end if
 20   continue
 30   call vclr(a,1,maxa)
c
c     loop over the transformed two-electron integrals
c     which are input from ed6 ( default)
c
      do 60 ifile = 1 , mmfile
         mblkk = kblk(ifile)
         idevm = nufile(ifile)
c        lblkm = nblk(ifile)
c
c
         call search(mblkk,idevm)
 40      call find(idevm)
c
c     read block of integrals into /blkin/
c
         call get(g(1),nw)
         if (nword.gt.0) then
            if (nw.gt.0) then
c
c     loop over integrals in a block
c
               call unpack(g(num2e+1),lab816,labs,numlab)
               do 50 kk = 1 , nword
c
c     unpack the labels
c
                  kk2 = kk + kk + kk + kk
                  i = labs(kk2-2)
                  j = labs(kk2-3)
                  k = labs(kk2  )
                  l = labs(kk2-1)
                  gg = g(kk)
                  if (i.gt.nocca) then
                     if (l.le.nocca) then
                        if (j.le.nocca) then
                           if (k.gt.nocca) then
c
c     type (xo/xo)
c
                              naa = (i-nocca-1)*nocca + j
                              nbb = (k-nocca-1)*nocca + l
                              if (naa.lt.nbb) then
                                 nswop = naa
                                 naa = nbb
                                 nbb = nswop
                              end if
                              if (naa.gt.nst .and. naa.le.nfin) then
                                 ntri = naa*(naa-1)/2 + nbb - nmin
                                 a(ntri) = a(ntri) + 4.0d0*gg
                              end if
                              naa = (i-nocca-1)*nocca + l
                              nbb = (k-nocca-1)*nocca + j
                              if (naa.lt.nbb) then
                                 nswop = naa
                                 naa = nbb
                                 nbb = nswop
                              end if
                              if (naa.gt.nst .and. naa.le.nfin) then
                                 ntri = naa*(naa-1)/2 + nbb - nmin
                                 a(ntri) = a(ntri) - gg*hf_wght
                              end if
                           end if
                        else if (k.le.nocca) then
c
c     type (xx/oo)
c
                           naa = (i-nocca-1)*nocca + k
                           nbb = (j-nocca-1)*nocca + l
                           if (naa.lt.nbb) then
                              nswop = naa
                              naa = nbb
                              nbb = nswop
                           end if
c
c     only that part of triangle between a(nst+1,1) and a(nfin,nfin)
c     is constructed in this pass
c
                           if (naa.gt.nst .and. naa.le.nfin) then
                              ntri = naa*(naa-1)/2 + nbb - nmin
                              a(ntri) = a(ntri) - gg*hf_wght
                           end if
                           if (i.ne.j .and. k.ne.l) then
                              naa = (i-nocca-1)*nocca + l
                              nbb = (j-nocca-1)*nocca + k
                              if (naa.lt.nbb) then
                                 nswop = naa
                                 naa = nbb
                                 nbb = nswop
                              end if
                              if (naa.gt.nst .and. naa.le.nfin) then
                                 ntri = naa*(naa-1)/2 + nbb - nmin
                                 a(ntri) = a(ntri) - gg*hf_wght
                              end if
                           end if
                        end if
                     end if
                  end if
c
 50            continue
               go to 40
            end if
         end if
 60   continue
      call wamat(a,mnnr,nofile(1),nst,nfin,odebug(2),iwr)
      if (nfin.eq.mnnr) then
c
         nword = 0
         mblk(1) = iposun(nofile(1)) - 1
         if (nprint.ne.-5) then
            dum = cpulft(1)
            write (iwr,6010) dum ,charwall()
         end if
         return
      else
         nst = nfin
         nfin = mnnr
         nmin = nst*(nst+1)/2
         nstp1 = nst + 1
         do 70 nn = nstp1 , mnnr
            last = nn*(nn+1)/2 - nmin
            if (last.gt.maxa) then
               nfin = nn - 1
               go to 30
            end if
 70      continue
      end if
      go to 30
 6010 format (/1x,'construction of a-matrix complete at ',f8.2,
     +        ' seconds',a10,' wall'/)
      end
      subroutine chfgrs(eps,mapnr,eta,zeta,a,iblok,iblok2,maxa)
c
c     a sorting routine for hessian for general scf chf equations
c
      implicit real*8  (a-h,o-z)
      dimension eta(*),zeta(nx,njk1),mapnr(*),a(*),eps(*)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/maxlen/maxq
      common/blkin/g(510),nword
      common/craypk/labs(1360)
      common/mpshl/ns(maxorb)
      character*10 charwall
      data prec,fattor/1.0d-9,2.0d0/
      data two, one /2.0d0, 1.0d0 /
      data half / 0.5d0 /
      ioff(i) = i*(i-1)/2
c
c     lennew = iky(nsa4+1)
c
      if (.not.lcpf .and. .not.cicv .and. (.not.cicx)) mnnr = mn
      do 20 iiii = 1 , 1360
         labs(iiii) = 0
 20   continue
c     lenblk = lensec(mn)
      call rdedx(eta,nx,iblok,ifils)
      call rdedx(zeta,nx*njk1,iblok2,ifils)
c
c      integrals are read from ed6 (nufile(1))
c      result is output on ed4 (nofile(1))
c      only integrals of type (xo/xo) and (xx/oo) are
c      required though more can be present on ed6.
c
c      input integrals in canonical order
c      lower triangle of a-matrix is output
c
c      several passes through integral file allowed
c
      call search(jblk(1),nofile(1))
c
c     work out number of passes
c     a matrix is mn*(mn+1)/2 ; core available is maxa
c
      npass = (mnnr*(mnnr+1)/2)/maxa + 1
      if(odebug(39)) then
       write(iwr,6020) npass
       write(iwr,6040) maxa, mnnr
      endif
      nst = 0
      nmina = 1
      nmin = 0
      nfin = mnnr
      nfinx = mnnr*(mnnr+1)/2
      do 30 nn = 1 , mnnr
         last = nn*(nn+1)/2
         if (last.gt.maxa) then
            nfin = nn - 1
            nfinx = nfin*(nfin+1)/2
            write(iwr,*) ' sortj : number of passes gt 1 '
            go to 40
         end if
 30   continue
 40   call vclr(a,1,maxa)
c
c     loop over the transformed two-electron integrals
c     which are input from ed6 ( default)
c
      do 70 ifile = 1 , mmfile
         mblkk = kblk(ifile)
         idevm = nufile(ifile)
c        lblkm = nblk(ifile)
c
         call search(mblkk,idevm)
 50      call find(idevm)
c
c     read block of integrals into /blkin/
c
         call get(g(1),nw)
         if (nword.gt.0) then
            if (nw.gt.0) then
c
c     loop over integrals in a block
c
               call unpack(g(num2e+1),lab816,labs,numlab)
               do 60 kk = 1 , nword
c
c     unpack the labels
c
                  kk2 = kk + kk + kk + kk

                  i = labs(kk2-2)
                  j = labs(kk2-3)
                  k = labs(kk2  )
                  l = labs(kk2-1)
                  gg = -g(kk)
                  ij = ioff(i) + j
                  kl = ioff(k) + l
                  ik = ioff(i) + k
                  jl = ioff(j) + l
                  if (l.gt.j) jl = ioff(l) + j
                  il = ioff(i) + l
                  jk = ioff(j) + k
                  if (k.gt.j) jk = ioff(k) + j
c
                  if (ijkltp(i,j,k,l).lt.4) go to 60
                  if (i.eq.j) gg = gg*half
                  if (k.eq.l) gg = gg*half
                  if (i.eq.k .and. j.eq.l) gg = gg*half
c
                  if (mapnr(ij).gt.0 .and. mapnr(kl).gt.0) then
                     ntri = ioff(mapnr(ij)) + mapnr(kl)
                     if (ntri.gt.nmin .and. ntri.le.nfinx) then
                        ntri = ntri - nmin
                        a(ntri) = a(ntri) + two*aijkl(i,j,k,l,ns,erga)
     +                            *gg
                     end if
                  end if
c
                  if (mapnr(ik).gt.0 .and. mapnr(jl).gt.0) then
                     ntri = ioff(mapnr(ik)) + mapnr(jl)
                     if (ntri.gt.nmin .and. ntri.le.nfinx) then
                        ntri = ntri - nmin
                        jsw = j
                        lsw = l
                        if (j.lt.l) then
                           jsw = l
                           lsw = j
                        end if
                        a(ntri) = a(ntri) + aijkl(i,k,jsw,lsw,ns,ergb)
     +                            *gg
                     end if
                  end if
c
                  if (mapnr(il).gt.0 .and. mapnr(jk).gt.0) then
                     ntri = ioff(mapnr(il)) + mapnr(jk)
                     if (mapnr(il).lt.mapnr(jk)) ntri = ioff(mapnr(jk))
     +                   + mapnr(il)
                     if (ntri.gt.nmin .and. ntri.le.nfinx) then
                        ntri = ntri - nmin
                        jsw = j
                        ksw = k
                        if (j.lt.k) then
                           jsw = k
                           ksw = j
                        end if
                        a(ntri) = a(ntri) + aijkl(i,l,jsw,ksw,ns,ergb)
     +                            *gg
                     end if
                  end if
c
 60            continue
               go to 50
            end if
         end if
 70   continue
c
      do 80 ixj = nmina , nfin
         ixx = ixj*(ixj+1)/2 - nmin
         a(ixx) = a(ixx)*two
 80   continue
      do 120 ix = 1 , nsa4
         do 110 jx = 1 , ix
            ijx = ioff(ix) + jx
            if (mapnr(ijx).ne.0) then
               do 100 kx = 1 , ix
                  lmax = kx
                  if (ix.eq.kx) lmax = jx
                  do 90 lx = 1 , lmax
                     klx = ioff(kx) + lx
                     if (mapnr(klx).ne.0) then
                        ntri = ioff(mapnr(ijx)) + mapnr(klx)
                        if (ntri.gt.nmin .and. ntri.le.nfinx) then
                           ntri = ntri - nmin
                           ikx = ioff(ix) + kx
                           ilx = ioff(ix) + lx
                           jlx = ioff(jx) + lx
                           jkx = ioff(jx) + kx
                           if (lx.gt.jx) jlx = ioff(lx) + jx
                           if (kx.gt.jx) jkx = ioff(kx) + jx
                           if (lx.eq.jx) a(ntri) = a(ntri)
     +                         - (eta(ikx)-zeta(ikx,ns(jx)))
                           if (jx.eq.kx) a(ntri) = a(ntri)
     +                         + (eta(ilx)-zeta(ilx,ns(jx)))
                           if (kx.eq.ix) a(ntri) = a(ntri)
     +                         - (eta(jlx)-zeta(jlx,ns(ix)))
                           if (lx.eq.ix) a(ntri) = a(ntri)
     +                         + (eta(jkx)-zeta(jkx,ns(ix)))
                        end if
                     end if
 90               continue
 100           continue
            end if
 110     continue
 120  continue
c
c
c    this sets up the array eps(mn) to hold the diagonal values of
c    a which are used to divide b in the solution of (1-a)u=b
c
      do 140 ixj = nmina , nfin
         ixx = ixj*(ixj+1)/2 - nmin
c     new partition of the a matrix into diagonal-off
c     diagonal components
         aixx = a(ixx)
         aaixx = dabs(aixx)
         ixi0 = ixx - ixj + 1
         ixi1 = ixx - 1
c    scan last row of a matrix and find largest absolute value
         aaizz = 0.d0
         do 130 ixi = ixi0 , ixi1
            aiyy = a(ixi)
            aaiyy = dabs(aiyy)
            if (aaiyy.gt.aaizz) then
               aaizz = aaiyy
            end if
 130     continue
         if (aaixx.lt.prec) then
            diago = 1.0d0
            eps(ixj) = one/diago
            a(ixx) = a(ixx) - diago
         else if (aaizz.lt.prec) then
            eps(ixj) = one/aixx
            a(ixx) = 0.d0
         else if (aaixx/aaizz.gt.fattor) then
            eps(ixj) = one/aixx
            a(ixx) = 0.d0
         else
            diago = aixx*aaizz*fattor
            eps(ixj) = one/diago
            a(ixx) = a(ixx) - diago
         end if
 140  continue
      call wamat(a,mnnr,nofile(1),nst,nfin,odebug(2),iwr)
      if (nfin.eq.mnnr) then
c
         nword = 0
         mblk(1) = iposun(nofile(1)) - 1
         write (iwr,6010) cpulft(1) ,charwall()
         return
      else
         nst = nfin
         nmina = nfin + 1
         nfin = mnnr
         nfinx = mnnr*(mnnr+1)/2
         nmin = nst*(nst+1)/2
         nstp1 = nst + 1
         do 150 nn = nstp1 , mnnr
            last = nn*(nn+1)/2 - nmin
            if (last.gt.maxa) then
               nfin = nn - 1
               nfinx = nfin*(nfin+1)/2
               go to 160
            end if
 150     continue
 160     write(iwr,6030)cpulft(1) ,charwall()
         go to 40
      end if
 6010 format (/1x,'construction of a-matrix complete at',
     +  f8.2,' seconds',a10,' wall'/)
 6020 format (/1x,'rough estimate of number of passes ',i3)
 6030 format (/1x,'commence next pass at ',f8.2,' seconds',a10,' wall')
 6040 format (1x,'maxa,mnnr = ' , 2i10)
      end
      subroutine chfops(akm,maxa)
      implicit real*8  (a-h,o-z)
      dimension akm(*)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c     hessian sorting routine for open shell chf
c
      common/craypk/labs(1360)
      common/small/labout(1360)
      common/blkin/g(510),nword
      common/out/a(510),nn
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
      common/maxlen/maxq
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      data m0,m1/0,1/
      data m101,m102,m103,m104,m105,m106,m107,m108,m109,m110,m111,m112/
     $21,22,23,24,25,26,27,28,29,30,31,32/
      data m131,m132,m133,m134,m135,m136,m137,m138,m139,m140,m141,m142/
     $41,42,43,44,45,46,47,48,49,50,51,52/
      data m151,m152,m153,m154,m155,m156,m157,m158,m159,m160,m161,m162/
     $61,62,63,64,65,66,67,68,69,70,71,72/
      data m113/80/
      data fd,fs,fv/2.0d0,1.0d0,0.0d0/
c
      ind(i,j) = iky(max(i,j)) + min(i,j)
c
c     read the k-matrix
c
      call secget(isect(43),43,iblok)
      call rdedx(akm(1),lds(isect(43)),iblok,ifild)
c
c      map to active orbitals only
c
      call actmot(akm,nsa4,mapie,iky)
c
c      routine sorto sorts the two-electron integrals
c      the integrals are multiplied
c      by the required weighting factor and the indices ii
c      and jj indicate the array elements to which each
c      integral contributes  ---- i.e. sets up a small
c      formula tape.
c      open shell case. ii is i,j. jj is k,l. i.e. ordering of
c      indices reversed from that in p. saxe's paper.
c      entire matrix is output, not lower triangle.
c      sorter is for real perturbations
c
      do 20 iiii = 1 , 1360
         labout(iiii) = 0
 20   continue
      call setsto(1360,0,labs)
      fsdvd = -1.0d0
      fvdvd = -1.0d0
      fvsvd = -1.0d0
      fsdsd = -(fs+fd-1.0d0)/2.0d0
      fvdsd = -(fv+fd-1.0d0)/2.0d0
      fvssd = -(fv+fs-1.0d0)/2.0d0
      fsdvs = -(3.0d0-fs-fd)/2.0d0
      fvdvs = -(3.0d0-fv-fd)/2.0d0
      fvsvs = -(3.0d0-fv-fs)/2.0d0
c
c     integrals are read in from the transformed mainfile (ed6)
c     and are written out to secondary mainfile (ed4 = idev4)
c
      idev4 = ifils
      iblk4 = iblks
      call search(iblk4,idev4)
      nsoc = noccb - nocca
      nvirta = nsa4 - noccb
      ipos = m0
c
c     read through integral file
c
      do 470 ifile = 1 , mmfile
         mblkk = kblk(ifile)
         idevm = nufile(ifile)
         lblkm = nblk(ifile)
         call search(mblkk,idevm)
         do 460 ib = mblkk , lblkm
            call find(idevm)
            call get(g(1),nw)
            if (nw.le.m0) go to 470
            call unpack(g(num2e+1),lab816,labs,numlab)
            do 450 kk = 1 , nword
c
c     unpack integral labels
c
               kk2 = kk + kk + kk + kk
               i = labs(kk2-2)
               j = labs(kk2-3)
               k = labs(kk2  )
               l = labs(kk2-1)
               gg = g(kk)
               if (k.gt.noccb) then
                  if (l.gt.noccb) go to 450
                  if (j.gt.noccb) go to 450
                  if (l.le.nocca) then
                     if (j.gt.nocca) go to 240
c
c     type (vd/vd)
c
                     naa = (i-nocca-m1)*nocca + j
                     nbb = (k-nocca-m1)*nocca + l
                     ipos = ipos + m1
                     label = m103
                     xx = 4.0d0
                     if (i.eq.k .or. j.eq.l) xx = 4.0d0 + fvdvd
                     a(ipos) = xx*gg
                     if (i.eq.k) a(ipos) = a(ipos) + 2.0d0*akm(ind(j,l))
                     if (j.eq.l) a(ipos) = a(ipos) + 2.0d0*akm(ind(i,k))
                     labout(ipos+ipos-1) = naa
                     labout(ipos+ipos) = nbb
                     if (ipos.ge.num2e) go to 400
                     go to 90
                  else
                     if (j.le.nocca) go to 240
c
c     type (vs/vs)
c
                     naa = nvirta*nocca + (i-nsoc-m1)*nsoc + j - nocca
                     nbb = nvirta*nocca + (k-nsoc-m1)*nsoc + l - nocca
                     ipos = ipos + m1
                     label = m109
                     xx = 2.0d0
                     if (i.eq.k .or. j.eq.l) xx = 2.0d0 + fvsvs
                     a(ipos) = xx*gg
                     if (i.eq.k) a(ipos) = a(ipos) + akm(ind(j,l))
                     if (j.eq.l) a(ipos) = a(ipos) + akm(ind(i,k))
                     labout(ipos+ipos-1) = naa
                     labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
                     go to 280
                  end if
               else if (k.le.nocca) then
                  if (j.le.nocca) go to 450
                  if (i.le.noccb) then
c
c     type (ss/dd)
c
                     naa = (i-nocca-m1)*nocca + k
                     nbb = (j-nocca-m1)*nocca + l
                     ipos = ipos + m1
                     label = m104
                     a(ipos) = fsdsd*gg
                     labout(ipos+ipos-1) = naa
                     labout(ipos+ipos) = nbb
                     if (ipos.ge.num2e) go to 400
                     go to 120
                  else if (j.le.noccb) then
c
c     type (vs/dd)
c
                     naa = (i-nocca-m1)*nocca + k
                     nbb = (j-nocca-m1)*nocca + l
                     ipos = ipos + m1
                     label = m106
                     a(ipos) = fvdsd*gg
                     labout(ipos+ipos-1) = naa
                     labout(ipos+ipos) = nbb
                     if (ipos.ge.num2e) go to 400
                     go to 180
                  else
c
c     type (vv/dd)
c
                     naa = (i-nocca-m1)*nocca + k
                     nbb = (j-nocca-m1)*nocca + l
                     ipos = ipos + m1
                     label = m110
                     a(ipos) = fvdvd*gg
                     labout(ipos+ipos-1) = naa
                     labout(ipos+ipos) = nbb
                     if (ipos.ge.num2e) go to 400
                     go to 310
                  end if
               else if (l.le.nocca) then
                  if (i.le.noccb) then
                     if (j.gt.nocca) go to 450
c
c     type (sd/sd)
c
                     naa = (i-nocca-m1)*nocca + j
                     nbb = (k-nocca-m1)*nocca + l
                     ipos = ipos + m1
                     label = m101
                     xx = 2.0d0
                     if (i.eq.k .or. j.eq.l) xx = 2.0d0 + fsdsd
                     a(ipos) = xx*gg
                     if (i.eq.k) a(ipos) = a(ipos) + akm(ind(j,l))
                     if (j.eq.l) a(ipos) = a(ipos) + akm(ind(i,k))
                     labout(ipos+ipos-1) = naa
                     labout(ipos+ipos) = nbb
                     if (ipos.ge.num2e) go to 400
                  else if (j.le.noccb) then
                     if (j.le.nocca) then
c
c     type (vd/sd)
c
                        naa = (i-nocca-m1)*nocca + j
                        nbb = (k-nocca-m1)*nocca + l
                        ipos = ipos + m1
                        label = m102
                        xx = 2.0d0
                        if (j.eq.l) xx = 2.0d0 + fvdsd
                        a(ipos) = xx*gg
                        if (j.eq.l) a(ipos) = a(ipos) + akm(ind(i,k))
                        labout(ipos+ipos-1) = naa
                        labout(ipos+ipos) = nbb
                        if (ipos.ge.num2e) go to 400
                        go to 60
                     else
c
c     type (vs/sd)
c
                        naa = nvirta*nocca + (i-nsoc-m1)*nsoc + 
     +                        j - nocca
                        nbb = (k-nocca-m1)*nocca + l
                        ipos = ipos + m1
                        label = m107
                        xx = 2.0d0
                        if (j.eq.k) xx = 2.0d0 + fvssd
                        a(ipos) = xx*gg
                        if (j.eq.k) a(ipos) = a(ipos) + akm(ind(i,l))
                        labout(ipos+ipos-1) = naa
                        labout(ipos+ipos) = nbb
                        if (ipos.ge.num2e) go to 400
                        go to 210
                     end if
                  else
c
c     type (vv/sd)
c
                     naa = nocca*nvirta + (i-nsoc-m1)*nsoc + k - nocca
                     nbb = (j-nocca-m1)*nocca + l
                     ipos = ipos + m1
                     label = m111
                     a(ipos) = fvsvd*gg
                     labout(ipos+ipos-1) = naa
                     labout(ipos+ipos) = nbb
                     if (ipos.ge.num2e) go to 400
                     go to 340
                  end if
               else
                  if (i.le.noccb) go to 450
                  if (j.le.noccb) then
                     if (j.gt.nocca) go to 450
c
c     type (vd/ss)
c
                     naa = nvirta*nocca + (i-nsoc-m1)*nsoc + k - nocca
                     nbb = (l-nocca-m1)*nocca + j
                     ipos = ipos + m1
                     label = m105
                     a(ipos) = fvssd*gg
                     labout(ipos+ipos-1) = naa
                     labout(ipos+ipos) = nbb
                     if (ipos.ge.num2e) go to 400
                     go to 150
                  else
c
c     type (vv/ss)
c
                     naa = nocca*nvirta + (i-nsoc-m1)*nsoc + k - nocca
                     nbb = nvirta*nocca + (j-nsoc-m1)*nsoc + l - nocca
                     ipos = ipos + m1
                     label = m112
                     a(ipos) = fvsvs*gg
                     labout(ipos+ipos-1) = naa
                     labout(ipos+ipos) = nbb
                     if (ipos.ge.num2e) go to 400
                     go to 370
                  end if
               end if
 30            if (naa.eq.nbb) go to 450
               ipos = ipos + m1
               label = m131
               a(ipos) = xx*gg
               if (i.eq.k) a(ipos) = a(ipos) + akm(ind(j,l))
               if (j.eq.l) a(ipos) = a(ipos) + akm(ind(i,k))
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 40            if (i.eq.k .or. j.eq.l) go to 450
               naa = (i-nocca-m1)*nocca + l
               nbb = (k-nocca-m1)*nocca + j
               ipos = ipos + m1
               label = m151
               a(ipos) = fsdsd*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 50            ipos = ipos + m1
               label = m113
               a(ipos) = fsdsd*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
               go to 450
 60            ipos = ipos + m1
               label = m132
               xx = 4.0d0
               if (j.eq.l) xx = 4.0d0 + fsdvd
               a(ipos) = xx*gg
               if (j.eq.l) a(ipos) = a(ipos) + 2.0d0*akm(ind(k,i))
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 70            if (j.eq.l) go to 450
               naa = (i-nocca-m1)*nocca + l
               nbb = (k-nocca-m1)*nocca + j
               ipos = ipos + m1
               label = m152
               a(ipos) = fvdsd*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 80            ipos = ipos + m1
               label = m113
               a(ipos) = fsdvd*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
               go to 450
 90            if (naa.eq.nbb) go to 450
               ipos = ipos + 1
               label = m133
               a(ipos) = xx*gg
               if (i.eq.k) a(ipos) = a(ipos) + 2.0d0*akm(ind(j,l))
               if (j.eq.l) a(ipos) = a(ipos) + 2.0d0*akm(ind(i,k))
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 100           if (i.eq.k .or. j.eq.l) go to 450
               naa = (i-nocca-m1)*nocca + l
               nbb = (k-nocca-m1)*nocca + j
               ipos = ipos + m1
               label = m153
               a(ipos) = fvdvd*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 110           ipos = ipos + m1
               label = m113
               a(ipos) = fvdvd*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
               go to 450
 120           if (naa.eq.nbb) go to 450
               ipos = ipos + m1
               label = m134
               a(ipos) = fsdsd*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 130           if (i.eq.j .or. k.eq.l) go to 450
               naa = (i-nocca-m1)*nocca + l
               nbb = (j-nocca-m1)*nocca + k
               ipos = ipos + m1
               label = m154
               a(ipos) = fsdsd*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 140           ipos = ipos + 1
               label = m113
               a(ipos) = fsdsd*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
               go to 450
 150           ipos = ipos + 1
               label = m135
               a(ipos) = fsdvs*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 160           if (k.eq.l) go to 450
               naa = nvirta*nocca + (i-nsoc-m1)*nsoc + l - nocca
               nbb = (k-nocca-m1)*nocca + j
               ipos = ipos + m1
               label = m155
               a(ipos) = fvssd*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 170           ipos = ipos + m1
               label = m113
               a(ipos) = fsdvs*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
               go to 450
 180           ipos = ipos + m1
               label = m136
               a(ipos) = fsdvd*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 190           if (k.eq.l) go to 450
               naa = (i-nocca-m1)*nocca + l
               nbb = (j-nocca-m1)*nocca + k
               ipos = ipos + m1
               label = m156
               a(ipos) = fvdsd*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 200           ipos = ipos + m1
               label = m113
               a(ipos) = fsdvd*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
               go to 450
 210           ipos = ipos + 1
               label = m137
               if (j.eq.k) xx = 2.0d0 + fsdvs
               a(ipos) = xx*gg
               if (j.eq.k) a(ipos) = a(ipos) + akm(ind(l,i))
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 220           if (j.eq.k) go to 450
               naa = nvirta*nocca + (i-nsoc-m1)*nsoc + k - nocca
               nbb = (j-nocca-m1)*nocca + l
               ipos = ipos + m1
               label = m157
               a(ipos) = fvssd*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 230           ipos = ipos + m1
               label = m113
               a(ipos) = fsdvs*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
               go to 450
c
c     type (vs/vd) and (vd/vs)
c
 240           naa = nvirta*nocca + (i-nsoc-m1)*nsoc + j - nocca
               nbb = (k-nocca-m1)*nocca + l
               if (j.lt.l) naa = nvirta*nocca + (k-nsoc-m1)*nsoc + l -
     +                          nocca
               if (j.lt.l) nbb = (i-nocca-m1)*nocca + j
               ipos = ipos + m1
               label = m108
               xx = 4.0d0
               if (i.eq.k) xx = 4.0d0 + fvsvd
               a(ipos) = xx*gg
               if (i.eq.k) a(ipos) = a(ipos) + 2.0d0*akm(ind(j,l))
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 250           ipos = ipos + m1
               label = m138
               xx = 2.0d0
               if (i.eq.k) xx = 2.0d0 + fvdvs
               a(ipos) = xx*gg
               if (i.eq.k) a(ipos) = a(ipos) + akm(ind(l,j))
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 260           if (i.eq.k) go to 450
               naa = nvirta*nocca + (k-nsoc-m1)*nsoc + j - nocca
               nbb = (i-nocca-m1)*nocca + l
               if (j.lt.l) naa = nvirta*nocca + (i-nsoc-m1)*nsoc + l -
     +                          nocca
               if (j.lt.l) nbb = (k-nocca-m1)*nocca + j
               ipos = ipos + m1
               label = m158
               a(ipos) = fvsvd*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 270           ipos = ipos + m1
               label = m113
               a(ipos) = gg*fvdvs
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
               go to 450
 280           if (naa.eq.nbb) go to 450
               label = m139
               ipos = ipos + m1
               a(ipos) = xx*gg
               if (i.eq.k) a(ipos) = a(ipos) + akm(ind(j,l))
               if (j.eq.l) a(ipos) = a(ipos) + akm(ind(i,k))
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 290           if (i.eq.k .or. j.eq.l) go to 450
               naa = nvirta*nocca + (i-nsoc-m1)*nsoc + l - nocca
               nbb = nvirta*nocca + (k-nsoc-m1)*nsoc + j - nocca
               ipos = ipos + m1
               label = m159
               a(ipos) = fvsvs*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 300           ipos = ipos + m1
               label = m113
               a(ipos) = fvsvs*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
               go to 450
 310           if (naa.eq.nbb) go to 450
               ipos = ipos + m1
               label = m140
               a(ipos) = fvdvd*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 320           if (i.eq.j .or. k.eq.l) go to 450
               naa = (i-nocca-m1)*nocca + l
               nbb = (j-nocca-m1)*nocca + k
               ipos = ipos + m1
               label = m160
               a(ipos) = fvdvd*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 330           ipos = ipos + m1
               label = m113
               a(ipos) = fvdvd*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
               go to 450
 340           ipos = ipos + m1
               label = m141
               a(ipos) = fvdvs*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 350           if (i.eq.j) go to 450
               naa = nvirta*nocca + (j-nsoc-m1)*nsoc + k - nocca
               nbb = (i-nocca-m1)*nocca + l
               ipos = ipos + m1
               label = m161
               a(ipos) = fvsvd*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 360           ipos = ipos + m1
               label = m113
               a(ipos) = fvdvs*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
               go to 450
 370           if (naa.eq.nbb) go to 450
               ipos = ipos + m1
               label = m142
               a(ipos) = fvsvs*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.ge.num2e) go to 400
 380           if (i.eq.j .or. k.eq.l) go to 450
               naa = nvirta*nocca + (i-nsoc-m1)*nsoc + l - nocca
               nbb = nvirta*nocca + (j-nsoc-m1)*nsoc + k - nocca
               ipos = ipos + m1
               label = m162
               a(ipos) = fvsvs*gg
               labout(ipos+ipos-1) = naa
               labout(ipos+ipos) = nbb
               if (ipos.ge.num2e) go to 400
 390           ipos = ipos + m1
               label = m113
               a(ipos) = fvsvs*gg
               labout(ipos+ipos-1) = nbb
               labout(ipos+ipos) = naa
               if (ipos.lt.num2e) go to 450
c
c
 400           nn = num2e
               call pack(a(num2ep+1),lab1632,labout,numlabp)
               call put(a,m511,idev4)
               do 410 iiii = 1 , 680
                  labout(iiii) = 0
 410           continue
               ipos = m0
               iblk4 = iblk4 + m1
               label1 = label/20
               label2 = label - 20*label1
               go to (420,430,440,450) , label1
 420           go to (30,60,90,120,150,180,210,250,280,310,340,370) ,
     +                label2
 430           go to (40,70,100,130,160,190,220,260,290,320,350,380) ,
     +                label2
 440           go to (50,80,110,140,170,200,230,270,300,330,360,390) ,
     +                label2
 450        continue
 460     continue
 470  continue
c
c     write out block
c
      nn = ipos
      if (nn.gt.0) then
               call pack(a(num2ep+1),lab1632,labout,numlabp)
         call put(a,m511,idev4)
         do 480 iiii = 1 , 680
            labout(iiii) = 0
 480     continue
      end if
      nn = m0
       call pack(a(num2ep+1),lab1632,labout,numlabp)
      call put(a,m0,idev4)
      call chfopo(akm,maxa,mn,nofile(1),jblk(1),mblk(1),ifils,iblks,
     +  nocca,nsoc,nvirta)
      call clredx
      return
      end
      subroutine symmu(qq,ibstar,skipp,iso,nshels)
c
c    symmetrise u-vectors ( chf solutions )
c    d2h and subgroups only
c
      implicit real*8  (a-h,o-z)
      logical skipp
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension skipp(3,nat)
      dimension qq(*),iso(nshels,*)
c
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      common/mpshl/ns(maxorb)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/bufb/ptr(3,144),ict(maxat,8)
      common/symmos/imos(8,maxorb)
      character *8 grhf
      data grhf/'grhf'/
      data one/1.0d0/
c
      nav = lenwrd()
      call readi(iso,nw196(5)*nav,ibl196(5),ifild)
      call rdedx(ptr(1,1),nw196(1),ibl196(1),ifild)
c
      do 40 ii = 1 , nshell
         ic = katom(ii)
         do 30 it = 1 , nt
            id = iso(ii,it)
            ict(ic,it) = katom(id)
 30      continue
 40   continue
c
      nuniq = 0
      do 60 n = 1 , nat
         do 50 nop = 1 , nt
            if (ict(n,nop).ne.n) go to 60
 50      continue
         nuniq = n
         go to 70
 60   continue
 70   ntpls1 = noccb + 1
      nplus1 = nocca + 1
      nsoc = noccb - nocca
      nvirta = nsa4 - noccb
      ibll = lensec(mn)
      iblku = ibstar
      ioff = mn*3 + 1
      nat3 = nat*3
c     read in u vectors
      do 80 n = 1 , nat3
         call rdedx(qq(ioff),mn,iblku,ifils)
         iblku = iblku + ibll
         ioff = ioff + mn
 80   continue
      mn3 = mn*3
c
c
c     loop over vectors
      do 270 n = 1 , nat
         if (.not.(skipp(1,n))) then
            ioff = n*mn3
c     copy vectors for atom n into work area
c
            do 90 i = 1 , mn3
               qq(i) = qq(ioff+i)
 90         continue
c
c     zero all elements related to components of atom n
c     by symmetry
c
            nsame = 0
            do 110 iop = 1 , nt
               niop = ict(n,iop)
               if (niop.eq.n) nsame = nsame + 1
               ioff = niop*mn3
               do 100 i = 1 , mn3
                  qq(ioff+i) = 0.0d0
 100           continue
 110        continue
            nsame = max(nsame,1)
            an = one/dble(nsame)
c
c     work along the elements of this vector
c loop over double-single and double-virtual
c
            if (scftyp.eq.grhf) then
               ij = 0
               do 160 i = 1 , nsa4
                  do 150 ia = 1 , i
                     if (ns(i).ne.ns(ia)) then
                        ij = ij + 1
c     loop over symmetry operations
                        do 140 iop = 1 , nt
                           isign = imos(iop,i)*imos(iop,ia)
                           sign = dble(isign)*an
                           niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                           do 130 nc = 1 , 3
                              ioff = (nc-1)*mn
                              npnc = (iop-1)*3 + nc
                              do 120 k = 1 , 3
                                 iof2 = (niop*3+k-1)*mn
                                 qq(iof2+ij) = qq(iof2+ij)
     +                              + sign*ptr(k,npnc)*qq(ioff+ij)
 120                          continue
 130                       continue
 140                    continue
                     end if
 150              continue
 160           continue
            else
               if (nocca.ne.0) then
                  do 210 i = 1 , nocca
                     do 200 ia = nplus1 , nsa4
                        ij = (ia-nocca-1)*nocca + i
c     loop over symmetry operations
                        do 190 iop = 1 , nt
                           isign = imos(iop,i)*imos(iop,ia)
                           sign = dble(isign)*an
                           niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                           do 180 nc = 1 , 3
                              ioff = (nc-1)*mn
                              npnc = (iop-1)*3 + nc
                              do 170 k = 1 , 3
                                 iof2 = (niop*3+k-1)*mn
                                 qq(iof2+ij) = qq(iof2+ij)
     +                              + sign*ptr(k,npnc)*qq(ioff+ij)
 170                          continue
 180                       continue
 190                    continue
 200                 continue
 210              continue
               end if
               if (noccb.ne.nocca) then
c open shell only - loop over single-virtual
                  do 260 i = nplus1 , noccb
                     do 250 ia = ntpls1 , nsa4
                        ij = nvirta*nocca + (ia-nsoc-1)*nsoc + i - nocca
c     loop over symmetry operations
                        do 240 iop = 1 , nt
                           isign = imos(iop,i)*imos(iop,ia)
                           sign = dble(isign)*an
                           niop = ict(n,iop)
c     niop is the atom equivalent to n under operation
                           do 230 nc = 1 , 3
                              ioff = (nc-1)*mn
                              npnc = (iop-1)*3 + nc
                              do 220 k = 1 , 3
                                 iof2 = (niop*3+k-1)*mn
                                 qq(iof2+ij) = qq(iof2+ij)
     +                              + sign*ptr(k,npnc)*qq(ioff+ij)
 220                          continue
 230                       continue
 240                    continue
 250                 continue
 260              continue
               end if
            end if
         end if
 270  continue
c
c
c     translational invariance
c
      iof2 = nuniq*mn3
      if (nuniq.ne.0) then
         do 280 i = 1 , mn3
            qq(iof2+i) = 0.0d0
 280     continue
         do 300 n = 1 , nat
            if (n.ne.nuniq) then
               ioff = n*mn3
               do 290 i = 1 , mn3
                  qq(iof2+i) = qq(iof2+i) - qq(ioff+i)
 290           continue
            end if
 300     continue
      end if
c
c     write it all out again
c
      iblku = ibstar
      ioff = mn3 + 1
      lenuu = ibll*nat3
      call secput(isect(66),66,lenuu,iblko)
      lds(isect(66)) = mn*nat3
      call revind
      do 310 n = 1 , nat3
         call wrt3(qq(ioff),mn,iblku,ifils)
         call wrt3(qq(ioff),mn,iblko,ifild)
         iblku = iblku + ibll
         iblko = iblko + ibll
         ioff = ioff + mn
 310  continue
      return
      end
      subroutine bfnshl(inshel,nsa)
c
c      determine what shell basis function belongs to
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension inshel(*)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
c      loop over shells
c
      do 30 i = 1 , njk1
         ns = nbshel(i)
         il = ilfshl(i)
         do 20 nb = 1 , ns
            n = iactiv(il+nb)
c
c     basis function n is in shell i
c
            inshel(n) = i
 20      continue
 30   continue
      do 40 i = 1 , nsa
         inshel(i) = inshel(mapie(i))
 40   continue
      return
      end
      subroutine wgrhf(q,iq,alpha,beta)
c
c     'w matrix'   required in assembly of second derivatives
c      for general open-shell case
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      logical vir,out,occ
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      common/mpshl/ns(maxorb)
      dimension q(*),iq(*),alpha(11,*),beta(11,*)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      common/craypk/labs(1360)
      common/blkin/g(510),nword
      common/maxlen/maxq
c
      vir(i) = ns(i).eq.njk1
      occ(i) = ns(i).ne.njk1
      out = odebug(16)
      call bfnshl(ns,nsa4)
      nat3 = nat*3
      nlen = nat3*nat3
      nsq = nsa4*nsa4
      netaa = nat3*nsq
      nzeta = njk1*nx
      novlp = nat3*nx
      iof = lenrel(nw196(5))
      iofs = nw196(5) + lenint(nat*nt)
c
      i10 = iofs + 1
      ioa1 = i10 + nlen
      iob1 = ioa1 + netaa
      ioc1 = iob1 + novlp
      iod1 = ioc1 + nx
      ioe1 = iod1 + nzeta
      nreq = ioe1 + nat3*mn
      if (nreq.gt.maxq) call caserr('insufficient core')
      call vclr(q(i10),1,nlen)
      ioa = ioa1
      call search(iochf(17),ifockf)
      do 20 n = 1 , nat3
         call reads(q(ioa),nsq,ifockf)
         ioa = ioa + nsq
 20   continue
      iob = iob1
      call search(iochf(14),ifockf)
      do 30 n = 1 , nat3
         call reads(q(iob),nx,ifockf)
         call actmot(q(iob),nsa4,mapie,iky)
         iob = iob + nx
 30   continue
      ioa = ioa1 - 1
      do 60 n = 1 , nat3
         do 50 i = 1 , nsa4
            do 40 j = 1 , i
               sum = q(ioa+(i-1)*nsa4+j)
               q(ioa+(i-1)*nsa4+j) = 2.0d0*q(ioa+(j-1)*nsa4+i)
               q(ioa+(j-1)*nsa4+i) = sum + sum
 40         continue
 50      continue
         ioa = ioa + nsq
 60   continue
      ioa = ioa1 - 1
      do 100 na = 1 , nat3
         iob = iob1 - 1
         do 90 nb = 1 , nat3
            sum = 0.0d0
            do 80 i = 1 , nsa4
               if (.not.(vir(i))) then
                  do 70 j = 1 , i
                     if (.not.(vir(j))) then
                        ij = iky(i) + j
                        sij = q(iob+ij)
                        if (i.eq.j) sij = sij*0.5d0
                        sum = sum +
     +                        sij*(q(ioa+(i-1)*nsa4+j)+
     +                             q(ioa+(j-1)*nsa4+i) )
                     end if
 70               continue
               end if
 80         continue
            q((na-1)*nat3+nb+i10-1) = q((na-1)*nat3+nb+i10-1) - sum
            iob = iob + nx
 90      continue
         ioa = ioa + nsq
 100  continue
      lenb = lensec(mn)
      ibeta = iblks + lenb*nat*6
      ibzeta = ibeta + lensec(nx)
      ioa = ioa1 - 1
      iob = iob1 - 1
      ioc = ioc1 - 1
      iod = iod1 - 1
      call lgrhfm(q(ioc1),ibeta,alpha,beta,fjk)
      call lgrhf(q(iod1),ibzeta,alpha,beta,fjk)
      do 140 n = 1 , nat3
         do 130 k = 1 , nsa4
            if (.not.(vir(k))) then
               do 120 i = 1 , nsa4
                  do 110 j = 1 , i
                     ik = min(i,k) + iky(max(i,k))
                     jk = min(j,k) + iky(max(j,k))
                     if (occ(j)) q(ioa+(j-1)*nsa4+i) = 
     +                     q(ioa+(j-1)*nsa4+i)
     +                   - 2.0d0*q(ioc+jk)*q(iob+ik)
     +                   - (q(ioc+ik)+q(iod+ik+(ns(j)-1)*nx))*q(iob+jk)
                     if (i.ne.j) then
                        if (occ(i)) q(ioa+(i-1)*nsa4+j)
     +                      = q(ioa+(i-1)*nsa4+j) - 2.0d0*q(ioc+ik)
     +                      *q(iob+jk)
     +                      - (q(ioc+jk)+q(iod+jk+(ns(i)-1)*nx))
     +                      *q(iob+ik)
                     end if
 110              continue
 120           continue
            end if
 130     continue
         ioa = ioa + nsq
         iob = iob + nx
 140  continue
      do 230 ifile = 1 , mmfile
         mblkk = kblk(ifile)
         idevm = nufile(ifile)
c        lblkm = nblk(ifile)
         call search(mblkk,idevm)
         call find(idevm)
 150     call get(g(1),nw)
         if (nword.gt.0) then
            if (nw.gt.0) then
               call find(idevm)
               call unpack(g(num2e+1),lab816,labs,numlab)
               do 220 int = 1 , nword
                  kk2 = (int+int) + (int+int)
                  i = labs(kk2-2)
                  j = labs(kk2-3)
                  k = labs(kk2  )
                  l = labs(kk2-1)
                  gg = -g(int)
                  if (i.eq.j) gg = gg*0.5d0
                  if (k.eq.l) gg = gg*0.5d0
                  if (i.eq.k .and. j.eq.l) gg = gg*0.5d0
                  ioa = ioa1 - 1
                  iob = iob1 - 1
                  if (occ(k) .and. occ(l)) then
                     kl = iky(k) + l
                     do 160 n = 1 , nat3
                        ijw = ioa + (j-1)*nsa4 + i
                        jiw = ioa + (i-1)*nsa4 + j
                        if (occ(j)) q(ijw) = q(ijw) + gg*q(iob+kl)
     +                      *(alpha(ns(j),ns(k))+alpha(ns(j),ns(l)))
     +                      *2.0d0
                        if (occ(i)) q(jiw) = q(jiw) + gg*q(iob+kl)
     +                      *(alpha(ns(i),ns(k))+alpha(ns(i),ns(l)))
     +                      *2.0d0
                        ioa = ioa + nsq
                        iob = iob + nx
 160                 continue
                     ioa = ioa1 - 1
                     iob = iob1 - 1
                  end if
                  if (occ(i) .and. occ(j)) then
                     ij = iky(i) + j
                     do 170 n = 1 , nat3
                        klw = ioa + (l-1)*nsa4 + k
                        lkw = ioa + (k-1)*nsa4 + l
                        if (occ(l)) q(klw) = q(klw) + gg*q(iob+ij)
     +                      *(alpha(ns(l),ns(i))+alpha(ns(l),ns(j)))
     +                      *2.0d0
                        if (occ(k)) q(lkw) = q(lkw) + gg*q(iob+ij)
     +                      *(alpha(ns(k),ns(i))+alpha(ns(k),ns(j)))
     +                      *2.0d0
                        ioa = ioa + nsq
                        iob = iob + nx
 170                 continue
                     ioa = ioa1 - 1
                     iob = iob1 - 1
                  end if
                  if (occ(i) .and. occ(k)) then
                     ik = iky(i) + k
                     do 180 n = 1 , nat3
                        jlw = ioa + (l-1)*nsa4 + j
                        ljw = ioa + (j-1)*nsa4 + l
                        if (occ(l)) q(jlw) = q(jlw) + gg*q(iob+ik)
     +                      *(beta(ns(l),ns(i))+beta(ns(l),ns(k)))
                        if (occ(j)) q(ljw) = q(ljw) + gg*q(iob+ik)
     +                      *(beta(ns(j),ns(i))+beta(ns(j),ns(k)))
                        ioa = ioa + nsq
                        iob = iob + nx
 180                 continue
                     ioa = ioa1 - 1
                     iob = iob1 - 1
                  end if
                  if (occ(j) .and. occ(l)) then
                     jl = min(j,l) + iky(max(j,l))
                     do 190 n = 1 , nat3
                        ikw = ioa + (k-1)*nsa4 + i
                        kiw = ioa + (i-1)*nsa4 + k
                        if (occ(k)) q(ikw) = q(ikw) + gg*q(iob+jl)
     +                      *(beta(ns(k),ns(j))+beta(ns(k),ns(l)))
                        if (occ(i)) q(kiw) = q(kiw) + gg*q(iob+jl)
     +                      *(beta(ns(i),ns(j))+beta(ns(i),ns(l)))
                        ioa = ioa + nsq
                        iob = iob + nx
 190                 continue
                     ioa = ioa1 - 1
                     iob = iob1 - 1
                  end if
                  if (occ(i) .and. occ(l)) then
                     il = iky(i) + l
                     do 200 n = 1 , nat3
                        jkw = ioa + (k-1)*nsa4 + j
                        kjw = ioa + (j-1)*nsa4 + k
                        if (occ(k)) q(jkw) = q(jkw) + gg*q(iob+il)
     +                      *(beta(ns(k),ns(i))+beta(ns(k),ns(l)))
                        if (occ(j)) q(kjw) = q(kjw) + gg*q(iob+il)
     +                      *(beta(ns(j),ns(i))+beta(ns(j),ns(l)))
                        ioa = ioa + nsq
                        iob = iob + nx
 200                 continue
                     ioa = ioa1 - 1
                     iob = iob1 - 1
                  end if
                  if (occ(j) .and. occ(k)) then
                     jk = min(j,k) + iky(max(j,k))
                     do 210 n = 1 , nat3
                        ilw = ioa + (l-1)*nsa4 + i
                        liw = ioa + (i-1)*nsa4 + l
                        if (occ(l)) q(ilw) = q(ilw) + gg*q(iob+jk)
     +                      *(beta(ns(l),ns(j))+beta(ns(l),ns(k)))
                        if (occ(i)) q(liw) = q(liw) + gg*q(iob+jk)
     +                      *(beta(ns(i),ns(j))+beta(ns(i),ns(k)))
                        ioa = ioa + nsq
                        iob = iob + nx
 210                 continue
                     ioa = ioa1 - 1
                     iob = iob1 - 1
                  end if
 220           continue
               go to 150
            end if
         end if
 230  continue
      ioe = ioe1 - 1
      ioa = ioa1 - 1
      call secget(isect(66),66,isec66)
      call search(isec66,ifild)
      do 240 n = 1 , nat3
         call reads(q(ioe+1),mn,ifild)
         ioe = ioe + mn
 240  continue
      ioe = ioe1 - 1
      do 280 na = 1 , nat3
         iob = iob1 - 1
         ioe = ioe1 - 1
         do 270 nb = 1 , nat3
            sum = 0.0d0
            ij = 0
            nr = 0
            do 260 i = 1 , nsa4
               do 250 j = 1 , i
                  ij = ij + 1
                  ijw = ioa + (j-1)*nsa4 + i
                  jiw = ioa + (i-1)*nsa4 + j
                  ss = q(iob+ij)
                  if (i.eq.j) ss = ss*0.5d0
                  sum = sum - (q(ijw)+q(jiw))*ss
                  if (ns(i).ne.ns(j)) then
                     nr = nr + 1
                     sum = sum + (q(ijw)-q(jiw))*q(ioe+nr)
                  end if
 250           continue
 260        continue
            q((nb-1)*nat3+na+i10-1) = q((nb-1)*nat3+na+i10-1) + sum
            iob = iob + nx
            ioe = ioe + mn
 270     continue
         ioa = ioa + nsq
 280  continue
      call rdedx(q(1),nw196(5),ibl196(5),ifild)
      if (out) then
         write (iwr,6010)
         call dr2sym(q(i10),q(ioa1),iq(1),iq(iof+1),nat,nat3,
     +               nshell)
         call prnder(q(i10),nat3,iwr)
      end if
      call secget(isect(60),60,isec46)
      call rdedx(q(ioa1),nlen,isec46,ifild)
      do 290 i = 1 , nlen
         q(ioa1-1+i) = q(ioa1-1+i) + q(i10-1+i)
 290  continue
      call dr2sym(q(ioa1),q(i10),iq(1),iq(iof+1),nat,nat3,
     +            nshell)
      call wrt3(q(ioa1),nlen,isec46,ifild)
      if (out) then
         write (iwr,6020)
         call prnder(q(ioa1),nat3,iwr)
      end if
      return
 6010 format (//' coupled hartree-fock contribution')
 6020 format (//' total so far')
      end
      function aijkl(i,j,k,l,ns,alpha)
      implicit real*8  (a-h,o-z)
      dimension ns(*),alpha(11,11)
      aijkl = alpha(ns(j),ns(k)) + alpha(ns(i),ns(l))
     +        - alpha(ns(i),ns(k)) - alpha(ns(j),ns(l))
      return
      end
      function aijklx(i,j,k,l,ns,alpha)
      implicit real*8  (a-h,o-z)
      dimension ns(*),alpha(11,11)
      aijklx = -alpha(ns(j),ns(k)) + alpha(ns(i),ns(l))
     +         + alpha(ns(i),ns(k)) - alpha(ns(j),ns(l))
      return
      end
      subroutine ijmapr(num,inshel,nrmap,nr,nrx)
c
      implicit real*8  (a-h,o-z)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
      parameter (maxbfn=255)
      common/jspec/nspj(maxbfn)
c
      dimension nrmap(*),inshel(*)
c
      ioff(i) = i*(i-1)/2
c
      nr = 0
      ij = 0
      do 30 i = 1 , num
         do 20 j = 1 , i
            ij = ij + 1
            nrmap(ij) = 0
            if (inshel(i).ne.inshel(j)) then
               nr = nr + 1
               nrmap(ij) = nr
            end if
c
c     result : i and j in same shell - pair (i,j) is
c              redundant - nrmap(ij)=0
c
c            : i and j in different shells - pair(i,j)
c              non-redundant - nrmap(ij)=nr -
c              nr = position of non-redundant pair
c
 20      continue
 30   continue
c
      nrx = 0
      do 50 i = ncore + 1 , nocc
         do 40 j = 1 , ncore
            ij = ioff(i) + j
            if (nrmap(ij).eq.0) then
               nrx = nrx - 1
               nrmap(ij) = nrx
            end if
 40      continue
 50   continue
c
      do 70 i = nupact + 1 , ncoorb
         do 60 j = nocc + 1 , nupact
            ij = ioff(i) + j
            if (nrmap(ij).le.0) then
               nrx = nrx - 1
               nrmap(ij) = nrx
            end if
 60      continue
 70   continue
      do 90 i = ncore + 1 , nupact
         do 80 j = ncore + 1 , i - 1
            if (nspj(i).ne.nspj(j)) then
               ij = ioff(i) + j
               if (nrmap(ij).eq.0) then
                  nrx = nrx - 1
                  nrmap(ij) = nrx
               end if
            end if
 80      continue
 90   continue
c
      nrx = -nrx
c
c     nrx = no of additional independent pairs
c
      return
      end
      function ijkltp(i,j,k,l)
c=================================================================
c     classify integrals by coincidences between labels
c
      implicit real*8  (a-h,o-z)
      if (j.lt.k) then
         if (j.lt.l) then
            if (k.ne.l) then
               ijkltp = 14
               return
            else
               ijkltp = 11
               return
            end if
         else if (j.eq.l) then
            if (i.ne.k) then
               ijkltp = 9
               return
            else
               ijkltp = 5
               return
            end if
         else if (i.ne.k) then
            ijkltp = 13
            return
         else
            ijkltp = 7
            return
         end if
      else if (j.eq.k) then
         if (i.ne.j) then
            if (k.ne.l) then
               ijkltp = 8
               return
            else
               ijkltp = 3
               return
            end if
         else if (k.ne.l) then
            ijkltp = 2
            return
         else
            ijkltp = 1
            return
         end if
      else if (i.ne.j) then
         if (k.ne.l) then
            ijkltp = 12
            return
         else
            ijkltp = 10
            return
         end if
      else if (k.ne.l) then
         ijkltp = 6
         return
      else
         ijkltp = 4
         return
      end if
      end
      subroutine actmot(a,nsa,mapie,iky)
      implicit real*8  (a-h,o-z)
c     condenses triangular array so that it contains only
c     elements over active m.o.s
c
      dimension a(*),mapie(*),iky(*)
      ija = 0
      do 30 i = 1 , nsa
         do 20 j = 1 , i
            ija = ija + 1
            ij = iky(mapie(i)) + mapie(j)
            if (ija.gt.ij) go to 40
            a(ija) = a(ij)
 20      continue
 30   continue
      return
 40   call caserr('error in actmot')
      end
      subroutine ver_cphf(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/cphf.m,v $
     +     "/
      data revision /"$Revision: 6317 $"/
      data date /"$Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
