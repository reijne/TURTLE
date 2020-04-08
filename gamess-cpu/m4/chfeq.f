      subroutine chfeq(q,eps,b,utotal,uold,au,u,maxc,lstop,skipp)
      implicit real*8  (a-h,o-z)
      character *8 fkd
      dimension q(*)
      dimension eps(*),b(*),utotal(*),uold(*),au(mn,*),u(mn,*)
c
c     simultaneous equations for chf - large case
c
      logical lstop,skipp
      dimension skipp(100)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
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
      common/craypk/labs(1360)
      common/blkin/g(510),nword
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      common/small/alpha(50,50),aa(50,50),bee(50),cc(50),wks1(50),
     + wks2(50)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      character*10 charwall
c
      data fkd/'fockder'/
c
c     iterative solution of simultaneous equations to give
c     coupled hartree-fock first order wavefunction
c
      if (maxc.gt.50) maxc = 50
      uconv = 10.0d0**(-iconvv)
      lenblk = lensec(mn)
c
c     rhs and solution stored on scratchfile
c     rhs beginning at block iblkb
c     sol. beginning at block iblku
c
      iblkb = iblks + npstar*lenblk
      iblku = iblkb + lenblk*np
      iblk4 = jblk(1)
      idev4 = nofile(1)
      call timit(3)
      tim0 = tim
      tim1 = tim
      iaa = 50
      ifail = 0
      write(iwr,6000) cpulft(1) ,charwall()
c
c     loop over all perturbations taken one at a time
c
      if (oprn(12)) write (iwr,6010)
      if (oprn(13)) write (iwr,6020)
      npst = npstar + 1
      do 150 ko = npst , np
c
c     get rhs for this perturbation
c
         call rdedx(b,mn,iblkb,ifils)
         if (oprn(12)) write (iwr,6030) ko , (b(i),i=1,mn)
         iblkb = iblkb + lenblk
c
c
c     initialise arrays
c
         call vclr(utotal,1,mn)
         if (odebug(2).and.skipp(ko)) write (iwr,6060) ko
         if (.not.(skipp(ko))) then
            call vclr(alpha,1,iaa*iaa)
            call vclr(u,1,mn*maxc)
            call vclr(bee,1,iaa)
            call vclr(au,1,mn*maxc)
c
c     get zeroth order estimate
c
            do 20 i = 1 , mn
               v = -b(i)*eps(i)
               bee(1) = bee(1) + v*v
               uold(i) = v
               u(i,1) = v
 20         continue
c
c     start of iterative solution of chf equations
c     50 iterations are allowed ---  usually less than 10
c     are neccessary
c
            if(oprn(6))write(iwr,6100)ko
            do 110 no = 1 , maxc
c
c     read in the combinations of the 2-electron integrals
c     corresponding to the hessian ( 'a-matrix' )
c
               call search(iblk4,idev4)
c
c     a-matrix on file idev4=nofile(1)=secondary mainfile= ed4 (default)
c     starting block jblk(1) = 1 (default)
c
               call find(idev4)
 30            call get(g(1),nw)
c     have got one block
c
               if (nw.gt.0) then
                  if (nword.gt.0) then
                     call find(idev4)
                     call unpack(g(num2ep+1),lab1632,labs,numlabp)
                     if (lcpf .or. cicv .or. cicx) then
                        do 40 i = 1 , nword
                           lab1 = labs(i+i-1)
                           lab2 = labs(i+i)
                           au(lab1,no) = au(lab1,no) + g(i)*u(lab2,no)
 40                     continue
                     else
                        do 50 i = 1 , nword
                           lab1 = labs(i+i-1)
                           lab2 = labs(i+i)
                           tt = au(lab1,no) + g(i)*u(lab2,no)
                           au(lab2,no) = au(lab2,no) + g(i)*u(lab1,no)
                           au(lab1,no) = tt
 50                     continue
c
c     go back for another block
c
                     end if
                     go to 30
                  end if
               end if
c
c
               do 60 i = 1 , mn
                  au(i,no) = au(i,no)*eps(i)
 60            continue
               alpha(no,no) = 
     *      ddot(mn,u(1,no),1,u(1,no),1)+ddot(mn,u(1,no),1,au(1,no),1)
               if (no.gt.1) then
                  no1 = no - 1
                  do 70 noo = 1 , no1
                     alpha(noo,no) = ddot(mn,u(1,noo),1,au(1,no),1)
                     alpha(no,noo) = ddot(mn,u(1,no),1,au(1,noo),1)
 70               continue
                  do 80 noo = 1 , no
                     bee(noo) = ddot(mn,u(1,noo),1,u(1,1),1)
 80               continue
                  nnn = no
c
c     nag routine to solve a small set of simultaneous equations
c
                  call f04atf(alpha,iaa,bee,nnn,cc,aa,iaa,wks1,wks2,
     +                        ifail)
               else
                  cc(1) = bee(1)/alpha(1,1)
               end if
               call mxmb(u,1,mn,cc,1,no,utotal,1,mn,mn,no,1)
c
c     check for convergence
c
               call vsub(utotal,1,uold,1,uold,1,mn)
               sumsq = ddot(mn,uold,1,uold,1)/dble(mn)
               sumsq = dsqrt(sumsq)
               if (oprn(6)) write (iwr,6040) no , sumsq
               if (sumsq.le.uconv) go to 130
               call vclr(uold,1,mn)
c
c    form new estimate of solution
c
               do 100 mo = 1 , no
                  alp = 
     +        ddot(mn,u(1,mo),1,au(1,no),1) /
     +        ddot(mn,u(1,mo),1,u(1,mo),1)
                  do 90 i = 1 , mn
                     uold(i) = uold(i) + alp*u(i,mo)
 90               continue
 100           continue
c
               call vsub(au(1,no),1,uold,1,u(1,no+1),1,mn)
               if (no.ge.maxc) go to 120
               call dcopy(mn,utotal,1,uold,1)
               call vclr(utotal,1,mn)
 110        continue
 120        write (iwr,6050) ko
         end if
         go to 140
 130     call timit(3)
         write (iwr,6070) no , ko , cpulft(1) ,charwall()
 140     call wrt3(utotal,mn,iblku,ifils)
         iblku = iblku + lenblk
         if (oprn(12) .or. oprn(13)) then
            write (iwr,6080) ko , (utotal(i),i=1,mn)
         end if
         dtim = tim - tim1
         tim1 = tim
         if ((timlim-tim).le.(dtim+dtim)) then
            if (fkder.eq.fkd .and. ko.ne.np) then
               lstop = .true.
               npfin = ko
               ti = tim0
               write (iwr,6090)
               return
            end if
         end if
 150  continue
      ti = tim0
      return
 6000 format(/1x,
     +'commence iterative solutions of chf equations at ',
     + f8.2,' seconds',a10,' wall')
 6010 format (//1x,'print right-hand-side to chf equations')
 6020 format (//1x,'print solutions to chf equations')
 6030 format (//1x,'perturbation  ',i4//(5x,5f16.8))
 6040 format (10x,i5,f15.10)
 6050 format (/10x,'convergence not achieved  ,  component',i4)
 6060 format (1x,'perturbation',i5,' omitted')
 6070 format (/1x,
     + 'chf converged at iteration',i4/1x,
     + 'chf for perturbation',i4,' complete at ',f8.2,' seconds'
     + ,a10,' wall')
 6080 format (//1x,'solution  ',i4//(5x,5f16.8))
 6090 format (//1x,'running out of time!'//)
 6100 format(/6x,'perturbation ',i4/
     +        6x,'iteration',9x,'tester'/
     +        6x,24('=')/)
      end
