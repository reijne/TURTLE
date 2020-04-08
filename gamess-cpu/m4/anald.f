c 
c  $Author: mrdj $
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/anald.m,v $
c  $State: Exp $
c  
      subroutine vcdden(u,c,bd,t,num,nocca,ncoorb,ndim,ncomp)
c
c     construction of a perturbed density which multiplies
c     derivative of overlap in vcd
c
      implicit real*8  (a-h,o-z)
      dimension u(ndim,ncomp),c(num,ncoorb),bd(num,num,ncomp)
      dimension t(num)
      data z/0.0d0/
      call vclr(bd,1,num*num*ncomp)
      do 70 n = 1 , ncomp
         do 60 iq = 1 , num
            call vclr(t,1,nocca)
            do 30 j = nocca + 1 , ncoorb
               ciqj = c(iq,j)
               if (ciqj.ne.z) then
                  jj = (j-nocca-1)*nocca
                  ciqj = ciqj + ciqj
                  do 20 i = 1 , nocca
                     t(i) = t(i) + ciqj*u(jj+i,n)
 20               continue
               end if
 30         continue
            do 50 i = 1 , nocca
               ti = t(i)
               if (ti.ne.z) then
                  do 40 ipp = 1 , num
                     bd(ipp,iq,n) = bd(ipp,iq,n) + c(ipp,i)*ti
 40               continue
               end if
 50         continue
 60      continue
 70   continue
      return
      end
      subroutine tdchf2 (eps,diag,v1,v2,v3,zm,aa,alpha,rhs,bb,ndim,mv)
c
      implicit real*8  (a-h,o-z)
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
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
      dimension eps(ndim),diag(ndim),v1(ndim,np),v2(ndim,np),v3(ndim,np)
     1         ,zm(mv,mv),aa(mv,mv),alpha(mv,np),rhs(mv,np),bb(mv,np)
      data accur /1.0d-6/
c
c     npv = np*mv
      npdim = np*ndim
      mv2 = mv**2
      lenblk = lensec(ndim)
      iblkr = iblks + lenblk*np
      iblkb = iblkr + lenblk*np
      idev4 = nofile(1)
      iblk4 = jblk(1)
      iconv = 0
      iter = 0
      nrstrt = 0
      nlast = 0
c     nstart = 1
      nvec = np
      nend = np
      call vclr(zm,1,mv2)
c
c...  load rhs & trial vectors
      call rdedv(v1,ndim,np,iblks,ifils)
      call search(iblkhi,idev4)
      call tdchhv(v1,v2,eps,ndim,np,idev4)
      do 30 ko = 1 , np
         do 20 i = 1 , ndim
            v2(i,ko) = v2(i,ko)/(diag(i)-freq(ifreq))
 20      continue
 30   continue
c..   we use last solution as trial vectors
      call rdedvs(v1,ndim,np,ifils)
c...  dump modified rhs
      call vwrt3(v2,ndim,np,iblkr,ifils)
c
c..   start of iteration
 40   iter = iter + 1
c
c...   orthogonalise
      icount = 0
 50   icount = icount + 1
      call tdchfo(v1,v3,bb,ndim,nvec,nlast,ifils,iblkb,1)
      zzz = 1.0d0
      do 80 ko = 1 , nvec
         if (ko.ne.1) then
            kom = ko - 1
            do 60 lo = 1 , kom
               bb(lo,1) = -ddot(ndim,v1(1,lo),1,v1(1,ko),1)
 60         continue
            call mxmb(v1,1,ndim,bb,1,0,v1(1,ko),1,0,ndim,kom,1)
         end if
         zz = 1.0d0/dsqrt(ddot(ndim,v1(1,ko),1,v1(1,ko),1))
         zzz = dmax1(zzz,zz)
         do 70 i = 1 , ndim
            v1(i,ko) = v1(i,ko)*zz
 70      continue
 80   continue
      if (zzz.gt.1.0d4 .and. icount.le.1) go to 50
c
c
c...  sigma vector
      call search(iblk4,idev4)
      call tdchhv(v1,v2,eps,ndim,nvec,idev4)
      call tdchhv(v2,v3,eps,ndim,nvec,idev4)
      do 100 ko = 1 , nvec
         do 90 i = 1 , ndim
            v2(i,ko) = (v3(i,ko)-freq(ifreq)*v1(i,ko))
     +                 /(diag(i)-freq(ifreq))
 90      continue
 100  continue
c...  write vectors
      do 110 ko = 1 , nvec
         call wrt3s(v1(1,ko),ndim,ifils)
         call wrt3s(v2(1,ko),ndim,ifils)
 110  continue
c
c...  small matrix
      call search(iblkb,ifils)
      do 140 ko = 1 , nend
         call reads(v3,ndim,ifils)
         do 120 lo = 1 , nvec
            zm(ko,lo+nlast) = ddot(ndim,v2(1,lo),1,v3,1)
 120     continue
         call reads(v3,ndim,ifils)
         do 130 lo = 1 , nvec
            zm(lo+nlast,ko) = ddot(ndim,v3,1,v1(1,lo),1)
 130     continue
 140  continue
c
c...  elements of right hand side
      call rdedv(v2,ndim,np,iblkr,ifils)
      do 160 io = 1 , np
         do 150 ko = 1 , nvec
            rhs(ko+nlast,io) = -ddot(ndim,v2(1,io),1,v1(1,ko),1)
 150     continue
 160  continue
c
c...  solve
      call f04aef(zm,mv,rhs,mv,nend,np,alpha,mv,v3,aa,mv,bb,mv,0)
c
c...  construct residual
      do 170 i = 1 , nend
         call reads(v3,ndim,ifils)
         call reads(v3,ndim,ifils)
         call mxmb(v3,1,0,alpha(i,1),1,mv,v2,1,ndim,ndim,1,np)
 170  continue
c      call filirt (v2,v3,bb,ndim,nvec,npole,ifild,iblkj,0)
c
c...  convergence test & construction of new vectors
      newvec = 0
      do 180 ivec = 1 , np
         bb(ivec,1) = dsqrt(ddot(ndim,v2(1,ivec),1,v2(1,ivec),1))
         if (bb(ivec,1).gt.accur) then
c
c...  new vector in series required
            newvec = newvec + 1
            call dcopy(ndim,v2(1,ivec),1,v1(1,newvec),1)
         end if
 180  continue
      write (iwr,6010) iter , (bb(ivec,1),ivec=1,np)
      if (newvec.eq.0) then
c
c...  final processing
         write (iwr,6020) iter , nrstrt , nend
         iconv = 1
      else if (newvec+nend.le.mv) then
c
c...  return for more
c        nprev = nvec
         nvec = newvec
         nlast = nend
c        nstart = nlast + 1
         nend = nlast + nvec
         go to 40
      end if
c
c...  entry here if not enough room for more expansion vectors
      call vclr(v1,1,npdim)
      call search(iblkb,ifils)
      do 190 ivec = 1 , nend
         call reads(v3,ndim,ifils)
         call mxmb(v3,1,ndim,alpha(ivec,1),1,mv,v1,1,ndim,ndim,1,np)
         call reads(v3,ndim,ifils)
 190  continue
      nlast = 0
      nend = np
      nvec = np
c     nstart = 1
      nrstrt = nrstrt + 1
      call vclr(zm,1,mv2)
      if (iconv.ne.1) go to 40
      call vwrt3(v1,ndim,np,iblkr,ifils)
      return
 6010 format (' iteration',i3,5x,'convergence:',11e8.1/30x,11e8.1)
 6020 format (1x,'convergence reached in',i3,' iterations (',i1,
     +        ' restarts,',i3,' expansion vectors)')
      end
      subroutine tdchf1 (eps,diag,v1,v2,v3,zm,zm1,vr,vi,rr,ri,intger
     1          ,ndim,mv)
c
      implicit real*8  (a-h,o-z)
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
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
c
      dimension eps(ndim),diag(ndim)
     1         ,v1(ndim,npole),v2(ndim,npole),v3(ndim,npole)
     2         ,zm(mv,mv),zm1(mv,mv),vr(mv,mv),vi(mv,mv)
     3         ,rr(mv),ri(mv),intger(mv)
      data accur /1.0d-6/
      data delta /0.1d0/
c
      npdim = npole*ndim
c     mv2 = mv**2
      lenblk = lensec(ndim)
      iblkb = iblks + lenblk*np*2
      idev4 = nofile(1)
      iblk4 = jblk(1)
      nrstrt = 0
      iter = 0
      nvec = npole
      nend = nvec
      nlast = 0
      nstart = 1
c
c...  select trial vectors
c
      call vclr(v2,1,ndim)
      call vclr(v1,1,npdim)
      do 30 j = 1 , npole
         zz = 1.0d25
         do 20 i = 1 , ndim
            if (diag(i).lt.zz .and. v2(i,1).le.0.05d0) then
               kk = i
               zz = diag(i)
            end if
 20      continue
         v2(kk,1) = dble(j) + 0.1d0
         v1(kk,j) = 1.0d0
 30   continue
      call search(iblkb,ifils)
c
c
c...  start of iteration
 40   iter = iter + 1
c..   construct sigma for these trial vectors
      call search(iblk4,idev4)
      call tdchhv(v1,v3,eps,ndim,nvec,idev4)
      call tdchhv(v3,v2,eps,ndim,nvec,idev4)
c
c
c
c..  write out vectors onto scratchfile .. assume positioned ok
      do 50 i = 1 , nvec
         call wrt3s(v1(1,i),ndim,ifils)
         call wrt3s(v2(1,i),ndim,ifils)
 50   continue
c
c
c...  update interaction matrix m
c..   read expansion vectors in
      call search(iblkb,ifils)
      do 80 i = 1 , nend
         call reads(v3,ndim,ifils)
         do 60 j = nstart , nend
            zm(i,j) = ddot(ndim,v3,1,v2(1,j-nlast),1)
 60      continue
c..   read sigma vector in
         call reads(v3,ndim,ifils)
         do 70 j = nstart , nend
            zm(j,i) = ddot(ndim,v1(1,j-nlast),1,v3,1)
 70      continue
 80   continue
c
c..   initialise hermiticity fudge factor
      epsiln = -delta
 90   epsiln = epsiln + delta
      do 110 i = 1 , nend
         do 100 j = 1 , nend
            zm1(i,j) = (zm(i,j)+epsiln*zm(j,i))/(epsiln+1.0d0)
 100     continue
 110  continue
c
c..   solve eigenproblem
c
      call f02agf(zm1,mv,nend,rr,ri,vi,mv,vr,mv,intger,0)
c
c..   complex roots ??
      do 120 i = 1 , nend
         if (dabs(ri(i)).gt.1.0d-10) go to 90
 120  continue
c
c..  sort on real parts of eigenvalues
      ifail = 0
      call m01ajf(rr,vr,intger,ri,nend,mv,ifail)
      do 130 i = 1 , nend
         call dcopy(nend,vi(1,intger(i)),1,vr(1,i),1)
 130  continue
c
      write (iwr,6010) iter , epsiln , (rr(i),i=1,npole)
c
c
c...   append correction vectors
      if (nend+nvec.gt.mv) then
c
c ====================================================================
c     have run out of expansion vectors
         igoto = 1
      else
         call search(iblkb,ifils)
         call vclr(v1,1,npdim)
         do 160 l = 1 , nend
            call reads(v3,ndim,ifils)
            call reads(v2,ndim,ifils)
            do 150 j = 1 , npole
               do 140 i = 1 , ndim
                  v1(i,j) = v1(i,j) + (v2(i,1)-rr(j)*v3(i,1))*vr(l,j)
 140           continue
 150        continue
 160     continue
         do 180 j = 1 , npole
            do 170 i = 1 , ndim
               v1(i,j) = v1(i,j)/(rr(j)-diag(i))
 170        continue
 180     continue
c
c
c..   orthonormalise
         nvec = npole
         icount = 0
 190     icount = icount + 1
         call tdchfo(v1,v3,vi,ndim,nvec,nend,ifils,iblkb,1)
         nvv = nvec
         nvec = 0
         do 220 i = 1 , nvv
            if (nvec.ne.0) then
               do 200 j = 1 , nvec
                  zz = -ddot(ndim,v1(1,i),1,v1(1,j),1)
                  call mxmb(v1(1,j),1,ndim,zz,1,1,v1(1,i),1,ndim,ndim,1,
     +                      1)
 200           continue
            end if
            zz = dsqrt(ddot(ndim,v1(1,i),1,v1(1,i),1))
            if (zz.gt.accur) then
c..   new expansion vector therefore required
               nvec = nvec + 1
               do 210 j = 1 , ndim
                  v1(j,nvec) = v1(j,i)/zz
 210           continue
            end if
 220     continue
         if (nvec.eq.0) then
c
c ====================================================================
c     convergence reached
            write (iwr,6020) iter , nrstrt , nend
            igoto = 0
         else
c...  re-do orthogonality because of round off
            if (icount.eq.1) go to 190
            nlast = nend
            nstart = nlast + 1
            nend = nlast + nvec
            go to 40
         end if
      end if
c
c...  assemble solution
      call vclr(v1,1,npdim)
      call search(iblkb,ifils)
      do 230 i = 1 , nend
         call reads(v2,ndim,ifils)
         call reads(v3,ndim,ifils)
         call mxmb(v2,1,ndim,vr(i,1),1,mv,v1,1,ndim,ndim,1,npole)
 230  continue
      if (igoto.eq.1) then
         nrstrt = nrstrt + 1
         nend = npole
         nlast = 0
         nstart = 1
         nvec = npole
         call search(iblkb,ifils)
         write (iwr,6030)
         go to 40
      else
c
c...  normalise
         call search(iblk4,idev4)
         call tdchhv(v1,v2,eps,ndim,npole,idev4)
         do 250 i = 1 , npole
            yy = 1.0d0/dsqrt(ddot(ndim,v2(1,i),1,v2(1,i),1))
            zz = 1.0d0/dsqrt(ddot(ndim,v1(1,i),1,v2(1,i),1))
            do 240 j = 1 , ndim
               v2(j,i) = v2(j,i)*yy
               v1(j,i) = v1(j,i)*zz
 240        continue
 250     continue
c
c...  relies on mv.ge.np
         call tdchfa(v1,rr,v2,ri,v3,ndim,npole)
         return
      end if
 6010 format (' iteration',i4,10x,'epsilon =',f4.1,10x,'eigenvalues:',
     +        10(/5f20.10))
 6020 format (1x,'convergence reached in',i3,' iterations (',i1,
     +        ' restarts,',i3,' expansion vectors)'/)
 6030 format (1x,'maximum number of expansion vectors reached; restart')
      end
      subroutine amtrm0 (b,prop,nprop)
c
      implicit real*8  (a-h,o-z)
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
      common/rtdata/whi(45),wlow(45),rhi(45),rlow(45),rfac(45),
     1amps(9),ipoint(9),madd(20)
      common/small/a(30),w(30)
      character *8 pn
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
      dimension prop(np,np)
      dimension b(nc6,nprop)
      dimension pn(50),iang(50)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      data pi/3.141592653589793d0/
      nb(i) = (i-1)/511 + 1
c
      maxper = 50
      nn = 0
      do 20 i = 1 , maxper
         if (.not.(opskip(i))) then
            nn = nn + 1
            pn(nn) = pnames(i)
            iang(nn) = ipang(i)
         end if
 20   continue
      write (iwr,6010)
      write (iwr,6020) (pn(i),i=1,nn)
c
c
      call secget(isect(52),52,iblkj)
      iblkj = iblkj + 2 + nb(npole) + nb(npole*np) + nb(np*np)
      call search(iblkj,ifild)
      do 50 i = 1 , nc6
         call reads(prop,np*np,ifild)
         l = 0
         do 40 j = 1 , np
            do 30 k = 1 , j
               l = l + 1
               b(i,l) = prop(j,k)
 30         continue
 40      continue
 50   continue
c
      go to (60,80) , ic6
c
c...  gauss-legendre using omega = w0 (1+t)/(1-t)
 60   ng2 = nc6/2
      do 70 i1 = 1 , ng2
         i2 = nc6 + 1 - i1
         ww = wlow(ipoint(ng2)+i1)*w0/pi
         rr = dsqrt(rlow(ipoint(ng2)+i1))
         w(i1) = ww/(1.0d0+rr)**2
         w(i2) = ww/(1.0d0-rr)**2
 70   continue
      go to 100
c
c...  mid-point using omega = w0 tan t
 80   do 90 i = 1 , nc6
         w(i) = w0/((dcos(dble(2*i-1)*pi/(dble(4*nc6))))
     +          **2*dble(4*nc6))
 90   continue
c
 100  write (iwr,6030)
      do 160 iordc = nc6min , nc6max
         if (oc6(iordc)) then
            iord = iordc - 2
            ij = 0
            write (iwr,6040) iord
            do 150 i = 1 , np
               do 140 j = 1 , i
                  ij = ij + 1
                  kl = 0
                  do 130 k = 1 , i
                     lmax = k
                     if (i.eq.k) lmax = j
                     do 120 l = 1 , lmax
                        kl = kl + 1
                        iiii = iang(i) + iang(j) + iang(k) + iang(l)
                        if (iiii.eq.iord) then
                           zz = 0.0d0
                           do 110 i1 = 1 , nc6
                              zz = zz + w(i1)*b(i1,ij)*b(i1,kl)
 110                       continue
                           if (dabs(zz).gt.1.0d-8) write (iwr,6050) 
     +                      pn(i) , pn(j) , pn(k) , pn(l) , zz
                        end if
 120                 continue
 130              continue
 140           continue
 150        continue
         end if
 160  continue
      if (opunch(5)) write (ipu,6060) (w(i),i=1,nc6)
      return
 6010 format (//1x,'perturbations'/
     +          1x,'*************'/)
 6020 format (1x,10a8)
 6030 format (1x,'dispersion energy integrals calculated by numerical',
     +        ' quadrature'/)
 6040 format (//1x,'****************************************'/
     +        /1x,'* terms for total angular momentum',i4,' *'/
     +        /1x,'****************************************')
 6050 format (4x,a8,1x,a8,4x,a8,1x,a8,f18.8)
 6060 format (1x,'integration weights for dispersion'/(1x,3e20.12))
      end
      subroutine chfidr(q,iq)
c
c     chf driver for imaginary perturbations....particularly
c     magnetizabilities
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
      logical lstop,skipp
      dimension skipp(100)
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
      common /maxlen/ maxq
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
      logical ogeompert_save
      common/mpshl/ns(maxorb)
      dimension iq(*),q(*)
c     character *8 grhf
      dimension dia(6,6)
c     data grhf/'grhf'/
c
c
      lenx = lensec(nx)
      call wrt3z(iblks,ifils,lenx*100)
      write (iwr,6010) gx , gy , gz
      i1 = nx + 1
      i2 = nx + i1
      i3 = nx + i2
      i4 = nx + i3
      i5 = i4 + nx
      i6 = i5 + nx
      i7 = i6 + nx
      i8 = i7 + nx
      if (i8.gt.maxq) then
         write (iwr,6020) i8 , maxq
         call caserr('insufficient core')
      endif
c
c     second moment integrals to give diamagnetic term
c
      i0 = igmem_alloc(nx)
      i1 = igmem_alloc(nx)
      i2 = igmem_alloc(nx)
      i3 = igmem_alloc(nx)
      i4 = igmem_alloc(nx)
      i5 = igmem_alloc(nx)
      i6 = igmem_alloc(nx)
      i7 = igmem_alloc(nx)
c
      call qmints(q(i0),q(i1),q(i2),q(i3),q(i4),q(i5))
      call onepdm(q(i6),q(i7))
      dia(1,1) = -tracep(q(i0),q(i6),num)
      dia(2,2) = -tracep(q(i1),q(i6),num)
      dia(3,3) = -tracep(q(i2),q(i6),num)
      dia(1,2) = 0.25d0*tracep(q(i3),q(i6),num)
      dia(1,3) = 0.25d0*tracep(q(i4),q(i6),num)
      dia(2,3) = 0.25d0*tracep(q(i5),q(i6),num)
      xx = 0.25d0*(dia(2,2)+dia(3,3))
      yy = 0.25d0*(dia(1,1)+dia(3,3))
      zz = 0.25d0*(dia(1,1)+dia(2,2))
      dia(1,1) = xx
      dia(2,2) = yy
      dia(3,3) = zz
      dia(2,1) = dia(1,2)
      dia(3,1) = dia(1,3)
      dia(3,2) = dia(2,3)
c
      call gmem_free(i7)
      call gmem_free(i6)
      call gmem_free(i5)
      call gmem_free(i4)
c
c      angular momentum integrals
c
      i5 = i4 + num*num
      i4 = igmem_alloc(num*num)
      i5 = igmem_alloc(num)
c
      call amints(q(i0),q(i1),q(i2),q(i3),q(i4),q(i5),.true.)
      call lmints(q(i0),q(i1),q(i2),q(i3),q(i4),q(i5),.true.)
c
      call gmem_free(i5)
      call gmem_free(i4)
      call gmem_free(i3)
      call gmem_free(i2)
      call gmem_free(i1)
      call gmem_free(i0)
c
c
      mn = noccb*nvirta + (noccb-nocca)*nocca
      call grhfbl(scftyp)
      call bfnshl(ns,nsa4)
      i0 = lenrel(mn) + 1
      imap  = igmem_alloc(lenint(nx))
      iimap = lenrel(imap-1)+1
      call nrmapo(iq(iimap),nocca,noccb,nsa4,iky,mapie)
      i01 = mn + lenint(nx) + 1
      i1 = i01 + nx
      nreq = mn*(np+np+1) + np*np + lenint(nx)
      if (nreq.gt.maxq) call caserr('insufficient core')
      call search(jblk(1),nofile(1))
c
c     magnetic hessian ( a-matrix for imaginary perturbations)
c
      i01 = igmem_alloc_all(maxa)
      call chficl(q(i01),maxa)
      call gmem_free(i01)
c
c     right-hand-side of equations
c
      ieps = igmem_alloc(mn)
      i01  = igmem_alloc(mn)
      i1   = igmem_alloc(nx)
      call rhsemp(q(ieps),q(i01),q(i1),iq(iimap))
      call gmem_free(i1)
      call gmem_free(i01)
c     solve chf equations
c
      do 20 n = 1 , np
         skipp(n) = .false.
 20   continue
      lstop = .false.
      ogeompert_save = ogeompert
      ogeompert = .false.
      call chfdrv(q(ieps),lstop,skipp)
      ogeompert = ogeompert_save
      if (lstop) then
         write (iwr,*) ' insufficient time.'
         call clenms('restart job')
      end if
c
      call gmem_free(ieps)
      call gmem_free(imap)
c
      i2 = i01 + mn*np
      i3 = i2 + nx
      i4 = i3 + nx
      i5 = i4 + num
      i5m = 566 + lenint(302)
      if (i5.lt.i5m) i5 = i5m
      write (iwr,6030)
c
c     magnetizability
c
      i01 = igmem_alloc(mn*np)
      i2  = igmem_alloc(mn*np)
      call magasm(q(i01),q(i2),dia,mn)
      lenuu = np*lensec(mn)
      call secput(isect(67),67,lenuu,iblko)
      do 30 ki = 1 , np
         call wrt3(q(i2+(ki-1)*mn),mn,iblko,ifild)
         iblko = iblko + lensec(mn)
 30   continue
      call revind
      call gmem_free(i2)
c
      i2  = igmem_alloc(nx)
      i3  = igmem_alloc(nx)
      imap  = igmem_alloc(lenint(nx))
      iimap = lenrel(imap-1)+1
      call nrmapo(iq(iimap),nocca,noccb,nsa4,iky,mapie)
c
c     low frequency limit of g-tensor
c
      if (lgten) write (iwr,6040)
      if (lgten) call gtensl(q(i01),q(i2),q(i3),iq(iimap),mn)
      call gmem_free(imap)
      call gmem_free(i3)
      call gmem_free(i2)
c
c     vibrational circular dichroism parameters
c
      i2  = igmem_alloc(num*num)
      i3  = igmem_alloc(mn)
      i4  = igmem_alloc(num)
      i5  = igmem_alloc(num*num*3)
      if (lvcd) then
         write (iwr,6050)
         call vcdpar(q(i01),q(i2),q(i3),q(i4),q(i5),mn)
         write (iwr,6060)
         call vrepdd(q(i01),q(i2),q(i3),q(i4),q(i5),mn)
      end if
      call gmem_free(i5)
      call gmem_free(i4)
      call gmem_free(i3)
      call gmem_free(i2)
      call gmem_free(i01)
c
c
      do 40 i = 1 , 9
         opskip(i) = .true.
 40   continue
      call timit(3)
      return
 6010 format (//1x,'gauge origin',3f12.6)
 6020 format (//1x,'store required',i8,'  available',i8)
 6030 format (/
     +  20x,'*********************'/
     +  20x,'** magnetisability **'/
     +  20x,'*********************'/)
 6040 format (/
     +  20x,'**********************'/
     +  20x,'** optical activity **'/
     +  20x,'**********************'/)
 6050 format (/
     +  20x,'*********************'/
     +  20x,'** vcd parameters  **'/
     +  20x,'*********************'/)
 6060 format (/
     +10x,'*************************************************'/
     +10x,'* velocity representation of dipole derivatives *'/
     +10x,'*************************************************'/)
      end
      subroutine abzez (n,a,b,z,eig,wkspce,ifail)
c
c..   solve ab z = e z
c..   a,b symmetric pos def
c..   z may overwrite a
c..   z normalised thro b
      implicit real*8  (a-h,o-z)
      dimension a(n,n),b(n,n),z(n,n),wkspce(n,2),eig(n)
      tol = x02adf(dum)
      eps = x02aaf(dum)
c
      call f01bdf(n,a,n,b,n,wkspce(1,1),ifail)
c
      call f01ajf(n,tol,a,n,eig,wkspce(1,2),z,n)
c
      call f02amf(n,eps,eig,wkspce(1,2),z,n,ifail)
c
      call f01aff(n,1,n,b,n,wkspce(1,1),z,n)
c
      return
      end
      subroutine tdchfd (eps,diag,ndim)
c
      implicit real*8  (a-h,o-z)
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
      common/blkin/g(510),nword
      dimension eps(ndim),diag(ndim)
c
      call search(jblk(1),nofile(1))
      do 20 i = 1 , ndim
         diag(i) = 1.0d0/eps(i)
 20   continue
      icount = 0
 30   icount = icount + 1
      call find(nofile(1))
 40   call get(g,nw)
      if (nw.eq.0) then
         if (icount.eq.1) go to 30
         do 50 i = 1 , ndim
            diag(i) = diag(i)/eps(i)
 50      continue
         return
      else
         call find(nofile(1))
         call unpack(g(num2ep+1),lab1632,labs,numlabp)
         do 60 i = 1 , nword
            lab = labs(i+i-1)
            if (lab.eq.labs(i+i)) then
               diag(lab) = diag(lab) + g(i)
            end if
 60      continue
         go to 40
      end if
      end
      subroutine tdchf(eps,lbig)
c
c     time dependent coupled hartree-fock routines
c
      implicit real*8  (a-h,o-z)
      logical  lbig
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
      common /small / abcis(30),weight(30)
      dimension eps(mn)
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
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
      logical qskip,lmix
      character *8 qnames
      common/qercmx/qnames(50)
      common/qercom/frqq(31),nfrqq,npoleb,iq6(36),
     & iqsec(50),qskip(50),iqang(50),npb,ijunq(4),lmix,ifileb,iblkb
c
      call secget(isect(52),52,iblkj)
      irestp = max(irestp,1)
c
      maxq  = igmem_max_memory()
      npdim = max(np,9)
      if (irestp.eq.4) go to 70
      if (irestp.lt.2) then
         npd = 0
         do 20 i = 1 , 3
            if (.not.opskip(i)) npd = npd + 1
 20      continue
         i1 = 1
         i2 = i1 + mn*npd
         nreq = i2 + mn*npd
         if (nreq.le.maxq .and. npd.ne.0) then
            i1 = igmem_alloc(mn*npd)
            i2 = igmem_alloc(mn*npd)
            call trksum(eps,qq(ivoff+i1),qq(ivoff+i2),npd,mn)
            call gmem_free(i2)
            call gmem_free(i1)
            call timit(3)
            irestp = 2
         end if
      end if
c
c
c..   attempt to do rpa in core
      i1 = 1
      i2 = i1 + mn*np
      i3 = i2 + mn*2
      i4 = i3 + mn
      i5 = i4 + mn**2
      nreq = i5 + mn**2
      idiag = -1
      write (iwr,6010) nreq , maxq
c      lbig=.true.
      if (nreq.gt.maxq .or. lbig) then
c
c...  large matrix routines
         idiag = igmem_alloc(mn)
         call tdchfd(eps,qq(ivoff+idiag),mn)
c
         nfreq = iabs(nfreq)
         if (npole.eq.0 .or. irestp.ge.3) then
            irestp = max(irestp,3)
         else
            npole = npole + 1
 30         npole = npole - 1
            i2 = i1 + mn
            i3 = i2 + mn*npole
            i4 = i3 + mn*npole
            i5 = i4 + mn*npole
            mv = 20 + 5*npole
            mvmin = max(np,10+2*npole)
            if (i5.gt.maxq) go to 30
 40         mv = mv - 1
            i6 = i5 + mv**2
            i7 = i6 + mv**2
            i8 = i7 + mv**2
            i9 = i8 + mv**2
            i10 = i9 + mv
            i11 = i10 + mv
            nreq = i11 + mv
            if (nreq.gt.maxq) go to 40
            if (mv.lt.mvmin) go to 30
            mv = min(mv,mn)
            i2  = igmem_alloc(mn*npole)
            i3  = igmem_alloc(mn*npole)
            i4  = igmem_alloc(mn*npole)
            i5  = igmem_alloc(mv*mv)
            i6  = igmem_alloc(mv*mv)
            i7  = igmem_alloc(mv*mv)
            i8  = igmem_alloc(mv*mv)
            i9  = igmem_alloc(mv)
            i10 = igmem_alloc(mv)
            i11 = igmem_alloc(mv)
            call tdchf1(eps,qq(ivoff+idiag),qq(ivoff+i2),qq(ivoff+i3),
     +           qq(ivoff+i4),qq(ivoff+i5),qq(ivoff+i6),qq(ivoff+i7),
     +           qq(ivoff+i8),qq(ivoff+i9),qq(ivoff+i10),qq(ivoff+i11),
     +           mn,mv)
            call timit(3)
            call gmem_free(i11)
            call gmem_free(i10)
            call gmem_free(i9)
            call gmem_free(i8)
            call gmem_free(i7)
            call gmem_free(i6)
            call gmem_free(i5)
            call gmem_free(i4)
            call gmem_free(i3)
            call gmem_free(i2)
c
            irestp = max(irestp,3)
         end if
      else
         npole = max(npole,1)
         i1 = igmem_alloc(mn*np)
         i2 = igmem_alloc(2*mn)
         i3 = igmem_alloc(mn)
         i4 = igmem_alloc(mn*mn)
         i5 = igmem_alloc(mn*mn)
         call tdchf0(eps,qq(ivoff+i1),qq(ivoff+i2),qq(ivoff+i3),
     +        qq(ivoff+i4),qq(ivoff+i5),mn)
         call gmem_free(i5)
         call gmem_free(i4)
         call gmem_free(i3)
         call gmem_free(i2)
         call gmem_free(i1)
         call timit(3)
         go to 70
      end if
c
c
 50   if (ifreq.lt.nfreq .and. irestp.lt.4) then
         ifreq = ifreq + 1
         nreqmx = 0
         mvmax = min(mn,15+5*np)
         mv = mvmax
         mvmin = min(mn-1,10+2*np)
 60      if (mv.le.mvmin) then
            if (idiag.gt.0) call gmem_free(idiag)
            write (iwr,6020) maxq , nreq , nreqmx
            return
         else
            mv = mv - 1
            i2 = 1
            i3 = i2 + mn*np
            i4 = i3 + mn*np
            i5 = i4 + mn*np
            i6 = i5 + mv**2
            i7 = i6 + mv**2
            i8 = i7 + mv*np
            i9 = i8 + mv*np
            nreq = i9 + mv*np
            nreqmx = max(nreqmx,nreq)
            if (nreq.gt.maxq) go to 60
            i2 = igmem_alloc(mn*np)
            i3 = igmem_alloc(mn*np)
            i4 = igmem_alloc(mn*np)
            i5 = igmem_alloc(mv*mv)
            i6 = igmem_alloc(mv*mv)
            i7 = igmem_alloc(mv*np)
            i8 = igmem_alloc(mv*np)
            i9 = igmem_alloc(mv*np)
            call tdchf2(eps,qq(ivoff+idiag),qq(ivoff+i2),qq(ivoff+i3),
     +           qq(ivoff+i4),qq(ivoff+i5),qq(ivoff+i6),
     +           qq(ivoff+i7),qq(ivoff+i8),qq(ivoff+i9),mn,mv)
            call gmem_free(i9)
            call gmem_free(i8)
            call gmem_free(i7)
            call gmem_free(i6)
            call gmem_free(i5)
            call gmem_free(i4)
            i4 = igmem_alloc(npdim*npdim)
            call polasm(qq(ivoff+i2),qq(ivoff+i3),qq(ivoff+i4),mn,npdim)
            call gmem_free(i4)
            call gmem_free(i3)
            call gmem_free(i2)
            call wrtc(pnames,ldsect(isect(52)),iblkj,ifild)
            call wrt3s(freq,lds(isect(52)),ifild)
            call timit(3)
            if (ifreq.lt.nfreq) then
               if ((tx*1.5d0).lt.(timlim-tim)) go to 50
               call revise
               write (iwr,6030)
               call clenms('restart job')
            end if
         end if
      end if
 70   irestp = 4
      if (idiag.gt.0) call gmem_free(idiag)
      nprop = (npa*(npa+1))/2
      if (ic6.eq.0) return
      i0 = igmem_alloc(nprop*nc6) 
      if (.not.lmix) then
         i1 = igmem_alloc(np*np)
         call amtrm0(qq(ivoff+i0),qq(ivoff+i1),nprop)
         call gmem_free(i1)
      else
         npropb = npb*(npb+1)/2
         i1 = igmem_alloc(npropb*nc6)
         i2 = igmem_alloc(npa*npa)
         i3 = igmem_alloc(npb*npb)
         call amtrm1(eps(1),eps(i1),eps(i2),eps(i3),nprop,npropb)
         call gmem_free(i3)
         call gmem_free(i2)
         call gmem_free(i1)
      end if
      call gmem_free(i0)
      return
 6010 format (1x,
     + 'memory for solving tdhf problem in core =',i8,6x,
     +        '(',i8,' available)'/)
 6020 format (/
     +  1x,'*******************************************************'/
     +  1x,'insufficient storage for iterative solution of chf eqns'/
     +  1x,'*******************************************************'/
     +  1x,'store available',i9/1x,'required - at least',i9,
     +        ' and preferably',i9//)
 6030 format (/10x,'***** insufficient time to continue *****')
      end
      subroutine eflint(dd,helf,natg,wrtout)
c
c    electric fields at nuclei
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
      logical iandj,wrtout,norm,double,allout
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/rtwt/xx,u(12),w(12),nroots
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
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
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
      common/junk/desp(3,maxat),
     +  xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
     + ,cx,cy,cz,ijn(225),ijx(225),ijy(225),ijz(225),
     + xin(250),yin(250),zin(250),dfac(225),dij(225)
      dimension helf(225,natg,3)
      dimension dd(*)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/small/de(3,maxat),dn(3,maxat)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      common/blkin/gout(510),ninth
c
      character *2 dnam
      dimension dnam(3)
      data ndim/5/
      data dnam/'ex','ey','ez'/
      data pi212 /1.1283791670955d0/
      data zero,one,two /0.0d0,1.0d0,2.0d0/
c
c
      write (iwr,6010) (dnam(i),i=1,3)
      i1 = 1 + nx
      call onepdm(dd(1),dd(i1))
c
      allout = odebug(19)
      if (allout) write (iwr,6020)
      if (allout) write (iwr,6030)
      tol = 2.30258d0*itol
      do 30 i = 1 , 3
         do 20 n = 1 , nat
            de(i,n) = zero
 20      continue
 30   continue
      iposf1 = iblks
      ninth = 1
      if (wrtout) call search(iposf1,ifils)
      norm = normf.ne.1 .or. normp.ne.1
c     ----- ishell
      do 250 ii = 1 , nshell
         i = katom(ii)
         xi = c(1,i)
         yi = c(2,i)
         zi = c(3,i)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 240 jj = 1 , ii
            j = katom(jj)
            xj = c(1,j)
            yj = c(2,j)
            zj = c(3,j)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
            nroots = (lit+ljt+1)/2
            nd2 = ndim*ndim
            nd2r = nd2*nroots
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
            iandj = ii.eq.jj
c     ----- prepare indices for pairs of (i,j) functions
            ij = 0
            max = maxj
            do 50 i = mini , maxi
               if (iandj) max = i
               do 40 j = minj , max
                  ij = ij + 1
                  ijn(ij) = iky(loci+i) + locj + j
                  dfac(ij) = two
                  if (iandj .and. i.eq.j) dfac(ij) = one
 40            continue
 50         continue
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,iandj,ndim,
     +           1,1)
            do 80 i = 1 , ij
               do 70 j = 1 , nat
                  do 60 k = 1 , 3
                     helf(i,j,k) = zero
 60               continue
 70            continue
 80         continue
c     ----- i primitive
            jgmax = j2
            do 170 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               if (iandj) jgmax = ig
               do 160 jg = j1 , jgmax
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     double = iandj .and. ig.ne.jg
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                               csj,cpj,cdj,cfj,cgj,
     +                           mini,maxi,minj,maxj,iandj,double,norm)
c     ..... hellmann-feynman term .....
                     dum = pi212/aa
                     dum = dum + dum
                     do 90 i = 1 , ij
                        dij(i) = dij(i)*dum
 90                  continue
                     aax = aa*ax
                     aay = aa*ay
                     aaz = aa*az
                     do 150 ic = 1 , nat
                        cx = c(1,ic)
                        cy = c(2,ic)
                        cz = c(3,ic)
                        xx = aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
                        if (nroots.le.3) call rt123
                        if (nroots.eq.4) call roots4
                        if (nroots.eq.5) call roots5
                        mm = 0
                        do 120 k = 1 , nroots
                           uu = aa*u(k)
                           ww = w(k)
                           ww = ww*uu
                           tt = aa + uu
                           t = dsqrt(tt)
                           tinv = one/tt
                           x0 = (aax+uu*cx)*tinv
                           y0 = (aay+uu*cy)*tinv
                           z0 = (aaz+uu*cz)*tinv
                           in = -ndim + mm
                           do 110 i = 1 , lit
                              in = in + ndim
                              ni = i
                              do 100 j = 1 , ljt
                                 jn = in + j
                                 nj = j
                                 call stvin2()
                                 xin(jn) = xint
                                 yin(jn) = yint
                                 zin(jn) = zint*ww
                                 call dvint
                                 xin(jn+nd2r) = xint
                                 yin(jn+nd2r) = yint
                                 zin(jn+nd2r) = zint*ww
 100                          continue
 110                       continue
                           mm = mm + nd2
 120                    continue
                        do 140 i = 1 , ij
                           nnx = ijx(i)
                           nny = ijy(i)
                           nnz = ijz(i)
                           dumx = zero
                           dumy = zero
                           dumz = zero
                           mm = 0
                           do 130 k = 1 , nroots
                              dumx = dumx + xin(nnx+mm+nd2r)*yin(nny+mm)
     +                               *zin(nnz+mm)
                              dumy = dumy + xin(nnx+mm)*yin(nny+mm+nd2r)
     +                               *zin(nnz+mm)
                              dumz = dumz + xin(nnx+mm)*yin(nny+mm)
     +                               *zin(nnz+mm+nd2r)
                              mm = mm + nd2
 130                       continue
                           dum = dij(i)
                           dumx = dumx*dum
                           dumy = dumy*dum
                           dumz = dumz*dum
                           helf(i,ic,1) = helf(i,ic,1) + dumx
                           helf(i,ic,2) = helf(i,ic,2) + dumy
                           helf(i,ic,3) = helf(i,ic,3) + dumz
                           dum = dd(ijn(i))*dfac(i)
                           de(1,ic) = de(1,ic) + dum*dumx
                           de(2,ic) = de(2,ic) + dum*dumy
                           de(3,ic) = de(3,ic) + dum*dumz
 140                    continue
 150                 continue
                  end if
 160           continue
 170        continue
            if (allout) then
               iij = 0
               do 200 i = mini , maxi
                  max = maxj
                  if (iandj) max = i
                  ipp = loci + i
                  do 190 j = minj , max
                     iij = iij + 1
                     jp = locj + j
                     do 180 k = 1 , nat
                        write (iwr,6040) ipp , jp , k ,
     +                                  (helf(iij,k,l),l=1,3) ,
     +                                  dd(ijn(iij))
 180                 continue
 190              continue
 200           continue
            end if
            if (wrtout) then
               do 230 i = 1 , ij
                  do 220 j = 1 , nat
                     do 210 k = 1 , 3
                        gout(ninth) = helf(i,j,k)
                        ninth = ninth + 1
 210                 continue
                     if (ninth.ge.m511) then
                        ninth = ninth - 1
                        call put(gout,m511,ifils)
                        ninth = 1
                        iposf1 = iposf1 + 1
                     end if
 220              continue
 230           continue
            end if
c     ----- end of 'primitive' loops -----
 240     continue
 250  continue
c     ----- end of 'shell' loops -----
      if (wrtout) then
         ninth = ninth - 1
         call put(gout,m511,ifils)
         iposf1 = iposf1 + 1
      end if
      call eflnuc(de,helf,nat,czan,c)
      do 260 i = 1 , nat
         write (iwr,6050) zaname(i) , (de(j,i),j=1,3)
 260  continue
      return
 6010 format (/
     +        20x,'****************************'/
     +        20x,'* electric field at nuclei *'/
     +        20x,'****************************'//
     +   14x,'atom',12x,a2,12x,a2,12x,a2///)
 6020 format (/10x,'electric field integrals')
 6030 format (//5x,'i',4x,'j',4x,'k',15x,'ex',18x,'ey',18x,'ez',18x,
     +        'dij')
 6040 format (1x,3i5,5x,4f20.8)
 6050 format (10x,a8,3(f14.7))
      end
      subroutine eflmo(q)
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
      logical iandj
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
      common/blkin/gout(510),ninth
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      parameter (maxn3=maxat*3)
      dimension ioff(maxn3)
      dimension q(*)
c
      nat3 = nat*3
c     length = lensec(nx)
      nbsq = num*ncoorb
      ntot = nx*nat3
c     nreq = ntot + nx + nbsq
c
      do 20 i = 1 , nat3
         ioff(i) = (i-1)*nx
 20   continue
      call vclr(q,1,ntot)
      nint = 1
      call search(iblks,ifils)
      call find(ifils)
      call get(gout,nw)
      do 70 ii = 1 , nshell
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
         do 60 jj = 1 , ii
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
            iandj = ii.eq.jj
            max = maxj
            nn = 0
            do 50 i = mini , maxi
               if (iandj) max = i
               do 40 j = minj , max
                  nn = nn + 1
                  ij = iky(loci+i) + locj + j
                  do 30 nnn = 1 , nat3
                     ioffn = ioff(nnn) + ij
                     q(ioffn) = q(ioffn) + gout(nint)
                     nint = nint + 1
                     if (nint.ge.m511) then
                        call find(ifils)
                        call get(gout,nw)
                        nint = 1
                     end if
 30               continue
 40            continue
 50         continue
 60      continue
 70   continue
c
c
c     ao electric field integrals in q(1) to q(ntot)
c
c     get vectors
c
      iv = ntot + 1
      id = iv + nbsq
      m = 0
      call secget(isect(8),m,iblok)
      call rdedx(q(iv),nbsq,iblok+mvadd,ifild)
      call vclr(q(id),1,nx)
      do 80 i = 1 , nat3
         ioffn = ioff(i) + 1
         call dcopy(nx,q(ioffn),1,q(id),1)
         call qhq1(q(ioffn),q(iv),ilifq,ncoorb,q(id),iky,num)
 80   continue
c
c      transformed integrals in q(1) to q(ntot)
c
      return
      end
      subroutine tdchfo (v,w,x,ndim,nv,nw,idev,ibl,nskip)
      implicit real*8  (a-h,o-z)
      dimension  v(ndim,nv),w(ndim),x(nv)
      call search(ibl,idev)
      if (nv.le.0 .or. nw.le.0 .or. ndim.le.0) return
      do 40 iw = 1 , nw
         call reads(w,ndim,idev)
         do 20 iv = 1 , nv
            x(iv) = -ddot(ndim,v(1,iv),1,w,1)
 20      continue
         call mxmb(w,1,0,x,0,1,v,1,ndim,ndim,1,nv)
         if (nskip.gt.0) then
            do 30 i = 1 , nskip
               call reads(w,ndim,idev)
 30         continue
         end if
 40   continue
      return
      end
      subroutine gtensl(u,da,db,mapnr,ndim)
c
c     low frequency limit of g-tensor (optical activity)
c
c     mixed length-velocity polarisability
      implicit real*8  (a-h,o-z)
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
      logical exist
c
      dimension gg(6,6),mapnr(*)
      dimension u(ndim,np),da(*),db(*)
      character*3 b,e,p
      dimension b(3),e(3),p(3)
      data p/'px','py','pz'/
      data b/'bx','by','bz'/
      data e/'ex','ey','ez'/
c
c
      call rdedv(u,ndim,np,iblks,ifils)
c
c     read the solution vectors in u
c
      call rdedvs(u,ndim,np,ifils)
c
      do 30 i = 1 , 6
         do 20 j = 1 , 6
            gg(i,j) = 0.0d0
 20      continue
 30   continue
      lennew = ncoorb*(ncoorb+1)/2
      do 50 j = 1 , 3
         call secloc(30+j,exist,iblock)
         if (exist) then
            call rdedx(da,lennew,iblock,ifild)
            call ijconr(da,db,lennew,mapnr)
c
c     this is g(e,b), and g(e,a)
c
            do 40 k = 1 , 6
               gg(j,k) = -2.0d0*ddot(mn,db,1,u(1,k),1)
 40         continue
         end if
 50   continue
      write (iwr,6020)
      write (iwr,6010) (b(i),i=1,3)
      do 60 i = 1 , 3
         write (iwr,6030) e(i) , (gg(i,j),j=1,3)
 60   continue
      do 80 i = 1 , 3
         do 70 j = 1 , 3
            gg(i,j) = 0.08724947d0*gg(i,j)
 70      continue
 80   continue
      write (iwr,6040)
      write (iwr,6010) (b(i),i=1,3)
      do 90 i = 1 , 3
         write (iwr,6030) e(i) , (gg(i,j),j=1,3)
 90   continue
      write (iwr,6060)
      write (iwr,6010) (p(i),i=1,3)
      do 100 i = 1 , 3
         write (iwr,6030) e(i) , (gg(i,j),j=4,6)
 100  continue
      return
 6010 format (20x,3(a3,13x))
 6020 format (/
     +    20x,'***********************************'/
     +    20x,'* low frequency limit of g-tensor *'/
     +    20x,'*    (atomic units --- bohr**4)   *'/
     +    20x,'***********************************'/)
 6030 format (/5x,a3,4x,3f16.8)
 6040 format (/
     +    20x,'******************************************'/
     +    20x,'* in s.i. units - 10**-50 farad meter**3 *'/
     +    20x,'******************************************'/)
 6060 format (/
     +20x,'****************************************************'/
     +20x,'* polarisability in length-velocity representation *'/
     +20x,'****************************************************'/)
      end
      subroutine tdchhm (a,eps,ndim,idev4)
      implicit real*8  (a-h,o-z)
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
      common/blkin/g(510),nword
      dimension a(ndim,ndim),eps(ndim)
c
      ndim2 = ndim**2
      call vclr(a,1,ndim2)
      do 20 j = 1 , ndim
         a(j,j) = 1.0d0/eps(j)
 20   continue
      call find(idev4)
 30   call get(g,nw)
      if (nw.eq.0) then
         do 50 i = 2 , ndim
            im = i - 1
            do 40 j = 1 , im
               a(i,j) = a(i,j) + a(j,i)
               a(j,i) = a(i,j)
 40         continue
 50      continue
         return
      else
         call find(idev4)
         call unpack(g(num2ep+1),lab1632,labs,numlabp)
         do 60 i = 1 , nword
            lab1 = labs(i+i-1)
            lab2 = labs(i+i)
            a(lab1,lab2) = a(lab1,lab2) + g(i)
 60      continue
         go to 30
      end if
      end
      subroutine tdchhv (v1,v2,eps,ndim,nvec,idev4)
      implicit real*8  (a-h,o-z)
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
       common/blkin/g(510),nword
       dimension v1(ndim,nvec),v2(ndim,nvec),eps(ndim)
c
      do 30 ivec = 1 , nvec
         do 20 j = 1 , ndim
            v2(j,ivec) = v1(j,ivec)/eps(j)
 20      continue
 30   continue
      call find(idev4)
 40   call get(g,nw)
      if (nw.eq.0) return
      call find(idev4)
         call unpack(g(num2ep+1),lab1632,labs,numlabp)
      do 60 i = 1 , nword
         lab1 = labs(i+i-1)
         lab2 = labs(i+i)
         do 50 ivec = 1 , nvec
            t = v2(lab1,ivec) + v1(lab2,ivec)*g(i)
            v2(lab2,ivec) = v2(lab2,ivec) + v1(lab1,ivec)*g(i)
            v2(lab1,ivec) = t
 50      continue
 60   continue
      go to 40
      end
      subroutine hfdder(q)
c
c   dipole derivatives as derivative of hellman-feynman force
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
      dimension q(*)
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
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
c     character *8 closed
      character *8 comp
      dimension comp(3)
      common/small/ dipd(3,3,maxat),dipn(3,3,maxat)
c     data closed/'closed'/
      data comp/'d/dx','d/dy','d/dz'/
c     nuclear term
      do 50 n = 1 , nat
         do 30 i = 1 , 3
            do 20 j = 1 , 3
               dipn(i,j,n) = 0.0d0
               dipd(i,j,n) = 0.0d0
 20         continue
 30      continue
         do 40 i = 1 , 3
            dipn(i,i,n) = czan(n)
 40      continue
 50   continue
c
      nat3 = nat*3
      ntot = nx*nat3
      id1 = ntot + 1
      id2 = id1 + nx
      id3 = id2 + nx
      id4 = id3 + nx
      nreq = id4 + 100*nat*3
      if (nreq.gt.maxq) then
         write (iwr,6010) nreq , maxq
         call caserr('insufficient core')
      end if
c
c     calculate electric field integrals
c
      call eflint(q,q(id4),nat,.true.)
c
c     transform to mo basis
c
      call eflmo(q)
      m = 0
      call secget(isect(31),m,iblok)
      call search(iblok,ifild)
      ltri = iky(ncoorb+1)
c
c      read in perturbed density matrices ( dipole perturbation)
c
      call reads(q(id1),ltri,ifild)
      call reads(q(id2),ltri,ifild)
      call reads(q(id3),ltri,ifild)
c
c
c
      ioff = 1
      do 70 n = 1 , nat
         do 60 ncomp = 1 , 3
            dipd(1,ncomp,n) = dipn(1,ncomp,n)
     +                        + tracep(q(ioff),q(id1),ncoorb)*czan(n)
            dipd(2,ncomp,n) = dipn(2,ncomp,n)
     +                        + tracep(q(ioff),q(id2),ncoorb)*czan(n)
            dipd(3,ncomp,n) = dipn(3,ncomp,n)
     +                        + tracep(q(ioff),q(id3),ncoorb)*czan(n)
            ioff = ioff + nx
 60      continue
 70   continue
c
c
      write (iwr,6020)
      do 90 n = 1 , nat
         write (iwr,6030)
         do 80 nc = 1 , 3
         write (iwr,6040) zaname(n),comp(nc),(dipd(nn,nc,n),nn=1,3)
 80      continue
 90   continue
c
      call dmdrot(c,dipd,nat)
      call timit(3)
      return
 6010 format (//10x,'core required',i10,' core available',i10)
 6020 format (/
     +   20x,'****************************************'/
     +   20x,'*  hellmann-feynman dipole derivatives *'/
     +   20x,'****************************************'//
     +   30x,'x',15x,'y',15x,'z'/)
 6030 format (/)
 6040 format (5x,a8,2x,a8,3f16.8)
      end
      subroutine hfhpol(q)
c
c     driving routine for scf hyperpolarisability calculation
c
      implicit real*8  (a-h,o-z)
      logical exist
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
      common/maxlen/maxq
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      logical okillv
      integer isecvv, itypvv
      common /vectrn/ isecvv,itypvv,okillv
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
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
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
      logical ciopt, ciforc, mp2, hfgr, bfgs, ump2, lmeth2
      logical ump3, rmp3, ordmo, mp2w, loptor, ladp, lcpf
      logical lopti, lmcscf, lforce, lci, lcart, lmcdat
      logical lfdtrn, unit7, lcontr, lvcd, lgten
      logical ldenom, ignore, ldens, lset, ladapt, lsym, latmol
      logical berny, llibry, limpt, fpres, oss, ldiag, lskip
      logical opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common /restrl/ ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,
     +rmp3,ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
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
c
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      common/scfopt/maxcyc,mconv,nconv
c
      common/junke/ibl5,ibl54,ibl56,maxt,ires,ipass,nteff,
     1     npass1,npass2,lentry,nbuck,mloww,mhi,ntri,iacc,iontrn
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
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
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
      common/crnams/iangc(40),iangs(24)
      character *8 pnamc,ac6,pnams
      character *4 bfnam,atnam
      common/crnamx/ac6(6),pnamc(40),pnams(24),bfnam(20),atnam(37)
      dimension q(*)
      dimension beta(3,3,3)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
      character *8 closed
      data closed/'closed'/
c
      if (lci .or. lmcscf .or. mp3 .or. mp2) call caserr(
     +'hyperpolarisability unavailable for correlated wavefunction')
      if (scftyp.ne.closed) call caserr(
     +    'hyperpolarisability for closed shells only')
c
      ldiag = .true.
      fkder=' '
c
c     evaluate integrals
c
      if (opass2) then
       if(nopk.eq.0.or.nopk.ne.nopkr.or.iofsym.eq.1.or.
     +    iofrst.ne.iofsym) then
        write (iwr,6030)
        opass2 = .false.
       endif
      endif
      nopk = 1
      iofsym = 0
      isecvv = isect(8)
      itypvv = 8
      nconv = max(nconv,7)
      if (mprest.lt.1) then
         call integ(q)
         mprest = 1
         call revise
      end if
      if (mprest.lt.2) then
         call scfrun(q)
         mprest = 2
         call revise
      end if
      if (opass6) then
            if (iscftp.lt.4)
     +      call caserr(' less restricted transformation required')
      end if
c
c  restrict the transformation
c
      iscftp = 4
c
      if (ldiag) iscftp = max(iscftp,4)
      npass1 = max(npas41,1)
      npass2 = max(npas42,1)
      iconvv = max(iconvv,9)
      lword4 = maxq
c
c   do the 4-index transformation and set
c
      if (mprest.ge.3) oprn(4) = .false.
c
      ibase = igmem_alloc_all(mword)
c     mmaxq = maxq
      maxq = mword
      lword4 = maxq
c
      call indxsy(q(ibase))
      oprn(4) = .false.
      if (mprest.lt.3) then
         call cpuwal(begin,ebegin)
         call indx2t(q(ibase))
         call indx4t(q(ibase))
         mprest = 3
         call revise
         call timana(11)
      end if
      if (iscftp.lt.4)
     + call caserr(' less restricted transformation required')
      ipol = 1
      np = 0
      maxp = 9
      if (.not.ogen) then
         do 30 i = 1 , maxp
            if (ione(i+3).ne.0) then
               np = np + 1
               opskip(i) = .false.
               ipsec(i) = isect(i+21)
            end if
 30      continue
         npa = np
      end if
      write (iwr,6020)
      call dipmom(q(ibase))
      call gmem_free(ibase)
      call poldrv(q,q)
      mprest = 0
      irest = 0
      do 40 i = 1 , 9
         opskip(i) = .true.
 40   continue
      call revise
c
      do 50 i = 70 , 72
         isee = i
         call secloc(isect(isee),exist,iblok)
         if (.not.exist) then
            write (iwr,6010)
            call caserr('derivative fock operators not present')
         end if
 50   continue
      ltri = iky(ncoorb+1)
      i0 = igmem_alloc(ltri)
      i1 = igmem_alloc(ltri)
      i2 = igmem_alloc(ltri)
      i3 = igmem_alloc(ltri)
      i4 = igmem_alloc(ncoorb)
      i5 = igmem_alloc(ncoorb*ncoorb)
      i6 = igmem_alloc(ncoorb*ncoorb)
      call diphyp(beta,q(i0),q(i1),q(i2),q(i3),q(i4),q(i5),q(i6),ltri,
     +  ncoorb,3)
      call gmem_free(i6)
      call gmem_free(i5)
      call gmem_free(i4)
      call gmem_free(i3)
      call gmem_free(i2)
      call gmem_free(i1)
      call gmem_free(i0)
      ldiag = .false.
 6020 format (//
     +1x,'************************************'/
     +1x,'scf dipole moment and polarisability'/
     +1x,'************************************'/)
6010  format(/
     +1x,'********************************************'/
     +1x,'hyperpolarisability calculation not possible'/
     +1x,'********************************************'//
     +1x,'the derivative fock operators are not present'/)
 6030 format(/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/
     + 1x,'= integrals must NOT be in supermatrix form ='/
     + 1x,'=        requested bypass is ignored        ='/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/)
      return
      end
      subroutine diphyp(beta,fa,fb,fc,x,e,ua,ub,ltri,norb,npert)
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
      dimension x(ltri),fa(ltri),fb(ltri),fc(ltri),e(norb),
     1   beta(npert,npert,npert),ua(norb,norb),ub(norb,norb)
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
      character *3 cx
      dimension cx(3)
      data cx/'x','y','z'/
c
      ind(i,j) = iky(max(i,j)) + min(i,j)
c
      call secget(isect(9),9,isec9)
      call rdedx(e,lds(isect(9)),isec9,ifild)
      call vclr(beta,1,npert*npert*npert)
      call vclr(ua,1,ncoorb*ncoorb)
c
c      perturbed fock matrices from section 2 thru 4 of fockfile
c
      do 140 na = 1 , 3
         iseca1 = na + 69
         call secget(isect(iseca1),iseca1,iblf)
         call rdedx(fa,ltri,iblf,ifild)
         do 30 i = nocc + 1 , ncoorb
            do 20 j = 1 , nocc
               diff = e(j) - e(i)
               if (dabs(diff).gt.1.0d-10) then
                  ua(i,j) = fa(ind(i,j))/diff
                  ua(j,i) = -ua(i,j)
               end if
 20         continue
 30      continue
         call vclr(ub,1,ncoorb*ncoorb)
         do 130 nb = 1 , 3
            isecb1 = nb + 69
            call secget(isect(isecb1),isecb1,iblf)
            call rdedx(fb,ltri,iblf,ifild)
            do 50 i = nocc + 1 , ncoorb
               do 40 j = 1 , nocc
                  diff = e(j) - e(i)
                  if (dabs(diff).gt.1.0d-10) then
                     ub(i,j) = fb(ind(i,j))/diff
                     ub(j,i) = -ub(i,j)
                  end if
 40            continue
 50         continue
            call vclr(x,1,ltri)
            do 80 i = 1 , nocc
               do 70 j = 1 , i
                  ij = iky(i) + j
                  do 60 k = nocc + 1 , ncoorb
                     x(ij) = x(ij) - ua(i,k)*ub(j,k) - ua(j,k)*ub(i,k)
 60               continue
 70            continue
 80         continue
            do 110 i = nocc + 1 , ncoorb
               do 100 j = nocc + 1 , i
                  ij = iky(i) + j
                  do 90 k = 1 , nocc
                     x(ij) = x(ij) + ua(i,k)*ub(j,k) + ua(j,k)*ub(i,k)
 90               continue
 100           continue
 110        continue
            do 120 nc = 1 , 3
               isecc1 = nc + 69
               call secget(isect(isecc1),isecc1,iblf)
               call rdedx(fc,ltri,iblf,ifild)
               abc = -tracep(fc,x,ncoorb)
               abc = abc + abc
               beta(na,nb,nc) = beta(na,nb,nc) + abc
               beta(nc,na,nb) = beta(nc,na,nb) + abc
               beta(nb,nc,na) = beta(nb,nc,na) + abc
 120        continue
 130     continue
 140  continue
      write (iwr,6010)
      do 170 i = 1 , 3
         do 160 j = 1 , 3
            do 150 k = 1 , 3
               write (iwr,6020) cx(i) , cx(j) , cx(k) , beta(i,j,k)
 150        continue
 160     continue
 170  continue
      call timit(3)
      return
 6010 format (/
     +  20x,'********************************'/
     +  20x,'*  dipole hyperpolarizability  *'/
     +  20x,'*         (atomic units)       *'/
     +  20x,'********************************'/)
 6020 format (20x,3a3,f18.8)
      end
      subroutine vcmrul(r,gg,nat,gauge)
c
c     rotational sum rules in vcd parameters
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension r(3,maxat),gg(3,maxat,6),rot(9),gauge(3)
      character*2 head
      dimension head(9)
      data head/'xx','xy','xz','yx','yy','yz','zx','zy','zz'/
      n = 0
      do 40 ibet = 1 , 3
         do 30 ialp = 1 , 3
            n = n + 1
            rot(n) = 0.0d0
            do 20 i = 1 , nat
               rot(n) = rot(n)
     +                  + 4.0d0*(gg(ialp,i,ibet)+gg(ialp,i,ibet+3))
 20         continue
 30      continue
 40   continue
      write (iwr,6010)
      do 100 i=1,9
 100  write (iwr,6020) head(i),rot(i)
      n = 0
      do 90 ibet = 1 , 3
         do 80 idel = 1 , 3
            n = n + 1
            rot(n) = 0.0d0
            do 70 ialp = 1 , 3
               do 60 igam = 1 , 3
                  skew = etijk(ialp,igam,idel)
                  if (skew.ne.0.0d0) then
                     do 50 i = 1 , nat
                        rot(n) = rot(n) + skew*(r(igam,i)-gauge(igam))
     +                           *gg(ialp,i,ibet)
 50                  continue
                  end if
 60            continue
 70         continue
 80      continue
 90   continue
      write (iwr,6030)
      do 110 i=1,9
 110  write (iwr,6020) head(i),rot(i)
      return
 6010 format (/
     +   10x,'*******************************'/
     +   10x,'* vcd translational sum rules *'/
     +   10x,'*          (dipole)           *'/
     +   10x,'*******************************'/)
 6020 format (10x,a3,3x,f11.6)
 6030 format (/
     +   10x,'****************************'/
     +   10x,'* vcd rotational sum rules *'/
     +   10x,'*     (magnetisability)    *'/
     +   10x,'****************************'/)
      end
      subroutine amtrm1(ba,bb,propa,propb,nprop,npropb)
c
      implicit real*8  (a-h,o-z)
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
      common/rtdata/whi(45),wlow(45),rhi(45),rlow(45),rfac(45),
     1amps(9),ipoint(9),madd(20)
      common /small / a(30),w(30)
      character *8 pn
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
      logical qskip,lmix
      character *8 qnames,qn
      common/qercmx/qnames(50)
      common/qercom/frqq(31),nfrqq,npoleb,iq6(36),
     & iqsec(50),qskip(50),iqang(50),npb,ijunq(4),lmix,ifileb,iblkb
c
      dimension propa(npa,npa),propb(npb,npb)
      dimension ba(nc6,nprop),bb(nc6,npropb)
      dimension pn(50),iang(50)
      dimension qn(50),jang(50)
      data pi/3.141592653589793d0/
c
c
      nb(i) = (i-1)/511 + 1
c
      maxper = 50
      nn = 0
      do 20 i = 1 , maxper
         if (.not.(opskip(i))) then
            nn = nn + 1
            pn(nn) = pnames(i)
            iang(nn) = ipang(i)
         end if
 20   continue
      write (iwr,6010)
      write (iwr,6020) (pn(i),i=1,nn)
      nn = 0
      do 30 i = 1 , maxper
         if (.not.(qskip(i))) then
            nn = nn + 1
            qn(nn) = qnames(i)
            jang(nn) = iqang(i)
         end if
 30   continue
      write (iwr,6030)
      write (iwr,6020) (qn(i),i=1,nn)
c
c
      call secget(isect(52),52,iblkj)
      iblkj = iblkj + 2 + nb(npole) + nb(npole*npa) + nb(npa*npa)
      call search(iblkj,ifild)
      do 60 i = 1 , nc6
         call reads(propa,npa*npa,ifild)
         l = 0
         do 50 j = 1 , npa
            do 40 k = 1 , j
               l = l + 1
               ba(i,l) = propa(j,k)
 40         continue
 50      continue
 60   continue
c
      iblkj = iblkb + 2 + nb(npoleb) + nb(npoleb*npb) + nb(npb*npb)
      call search(iblkj,ifilb)
      do 90 i = 1 , nc6
         call reads(propb,npb*npb,ifilb)
         l = 0
         do 80 j = 1 , npb
            do 70 k = 1 , j
               l = l + 1
               bb(i,l) = propb(j,k)
 70         continue
 80      continue
 90   continue
      go to (100,120) , ic6
c
c...  gauss-legendre using omega = w0 (1+t)/(1-t)
 100  ng2 = nc6/2
      do 110 i1 = 1 , ng2
         i2 = nc6 + 1 - i1
         ww = wlow(ipoint(ng2)+i1)*w0/pi
         rr = dsqrt(rlow(ipoint(ng2)+i1))
         w(i1) = ww/(1.0d0+rr)**2
         w(i2) = ww/(1.0d0-rr)**2
 110  continue
      go to 140
c
c...  mid-point using omega = w0 tan t
 120  do 130 i = 1 , nc6
         w(i) = w0/((dcos(dble(2*i-1)*pi/(dble(4*nc6))))
     +          **2*dble(4*nc6))
 130  continue
c
 140  write (iwr,6040)
      do 200 iorc = nc6min , nc6max
         if (oc6(iorc)) then
            iord = iorc - 2
            ij = 0
            write (iwr,6050) iord
            do 190 i = 1 , npa
               do 180 j = 1 , i
                  ij = ij + 1
                  kl = 0
                  do 170 k = 1 , npb
                     do 160 l = 1 , k
                        kl = kl + 1
                        iiii = iang(i) + iang(j) + jang(k) + jang(l)
                        if (iiii.eq.iord) then
                           zz = 0.0d0
                           do 150 i1 = 1 , nc6
                              zz = zz + w(i1)*ba(i1,ij)*bb(i1,kl)
 150                       continue
                           if (dabs(zz).gt.1.0d-8) write (iwr,6060) 
     +                      pn(i) , pn(j) , qn(k) , qn(l) , zz
                        end if
 160                 continue
 170              continue
 180           continue
 190        continue
         end if
 200  continue
      return
 6010 format (/1x,'**** perturbations on molecule 1')
 6020 format (1x,10a8)
 6030 format (/1x,'**** perturbations on molecule 2')
 6040 format (1x,'dispersion energy integrals calculated by numerical',
     +        ' quadrature'/)
 6050 format (//
     +    1x,'************************************'/
     +    1x,'terms for total angular momentum',i4/
     +    1x,'************************************'//)
 6060 format (4x,a8,1x,a8,4x,a8,1x,a8,f18.8)
      end
      subroutine eflnuc(de,drg,nat,zan,c)
c
c    nuclear components of fields at nuclei
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
      dimension  zan(maxat),c(3,maxat),de(3,maxat)
      dimension drg(nat,*)
      data zero,one/0.0d0,1.0d0/
      drg(1,1) = zero
      do 40 k = 2 , nat
         drg(k,k) = zero
         k1 = k - 1
         do 30 l = 1 , k1
            rkl = zero
            do 20 i = 1 , 3
               rkl = rkl + (c(i,k)-c(i,l))**2
 20         continue
            drg(k,l) = one/rkl
            drg(l,k) = dsqrt(rkl)
 30      continue
 40   continue
      do 90 kk = 1 , 3
         do 60 k = 2 , nat
            km1 = k - 1
            do 50 l = 1 , km1
               zal = zan(l)
               pkl = (c(kk,k)-c(kk,l))/drg(l,k)
               de(kk,k) = de(kk,k) + pkl*drg(k,l)*zal
 50         continue
 60      continue
         nat1 = nat - 1
         do 80 k = 1 , nat1
            kp1 = k + 1
            do 70 l = kp1 , nat
               zal = zan(l)
               pkl = (c(kk,k)-c(kk,l))/drg(k,l)
               de(kk,k) = de(kk,k) + pkl*drg(l,k)*zal
 70         continue
 80      continue
 90   continue
      return
      end
      subroutine pdrasm(pold,tempol,dd)
c
c     polarizability derivative routines - assembly section
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
c     logical out
      logical norm
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,
     1t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
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
      common/specal/ndenin,iflden,iblden
      common/scrtch/
     +    ptr(3,144),ddtr(6,288),ftr(10,480),gtr(15,720),
     +    dij(100),
     +    xin(25),yin(25),zin(25),xd(25),yd(25),zd(25),
     +    xdip(25),ydip(25),zdip(25),
     +    xdipd(25),ydipd(25),zdipd(25),
     +    ijx(100),ijy(100),ijz(100)
      dimension dd(*)
      dimension pold(3,3,*),tempol(3,3,*)
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
      logical exist
      dimension dtemp(3,3,3)
      character *8 closed
      data closed/'closed'/
c
      data ndim/5/
      data zer,one /0.0d0,1.0d0/
c
c     calculate derivatives of the dipole moment
c
      if (scftyp.ne.closed) call caserr(
     +   'polarisability derivatives only for closed-shell scf')
      do 20 i = 70 , 72
         isee = i
         call secloc(isect(isee),exist,iblok)
         if (.not.exist) then
            write (iwr,6010)
6010  format(/
     +1x,'*************************************************'/
     +1x,'the derivative fock operators were not calculated'/
     +1x,'in preceding polarisability calculation'/
     +1x,'*************************************************'/)
            call caserr('repeat chf with keepfock option set')
         end if
 20   continue
      if (iscftp.lt.4) then
         write (iwr,6020)
         call caserr(' rerun transformation and chf')
      end if
      origx = gx
      origy = gy
      origz = gz
      call search(iblden,iflden)
      ioff = 1
      do 30 i = 1 , ndenin
         call reads(dd(ioff),nx,iflden)
         ioff = ioff + nx
 30   continue
      ioffd = ioff
      iofft = ioffd + nx
      ioffv = iofft + num
      m = 0
      call secget(isect(8),m,iblok)
      call rdedx(dd(ioffv),num*ncoorb,iblok+mvadd,ifild)
      ioff = 1
      do 40 i = 1 , ndenin
         call dcopy(nx,dd(ioff),1,dd(ioffd),1)
         call demoao(dd(ioffd),dd(ioff),dd(ioffv),dd(iofft),num,ncoorb,
     +  num)
         ioff = ioff + nx
 40   continue
      tol = 2.30258d0*itol
c     out = odebug(21)
      norm = normf.ne.1 .or. normp.ne.1
      nat3 = nat*3
c     nuclear term
      do 70 n = 1 , nat3
         do 60 i = 1 , 3
            do 50 j = 1 , 3
               tempol(i,j,n) = zer
               pold(i,j,n) = zer
 50         continue
 60      continue
 70   continue
c     ----- ishell
      do 180 ii = 1 , nshell
         iat = katom(ii)
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 170 jj = 1 , nshell
            jat = katom(jj)
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
c           nroots = (lit+ljt+1)/2
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c     ----- prepare indices for pairs of (i,j) functions
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,.false.,
     +           ndim,1,1)
            do 80 icp = 1 , ndenin
               dtemp(icp,1,1) = 0.0d0
               dtemp(icp,2,2) = 0.0d0
               dtemp(icp,3,3) = 0.0d0
               dtemp(icp,1,2) = 0.0d0
               dtemp(icp,1,3) = 0.0d0
               dtemp(icp,2,1) = 0.0d0
               dtemp(icp,2,3) = 0.0d0
               dtemp(icp,3,1) = 0.0d0
               dtemp(icp,3,2) = 0.0d0
 80         continue
c     ----- i primitive
            do 150 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               do 140 jg = j1 , j2
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                               csj,cpj,cdj,cfj,cgj,
     +                           mini,maxi,minj,maxj,.false.,.false.,
     +                           norm)
c     ----- overlap
                     t = dsqrt(aa)
                     tinv = one/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     lit1 = lit + 1
                     in = -ndim
                     do 100 i = 1 , lit1
                        in = in + ndim
                        ni = i
                        do 90 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call dmsint()
                           xin(jn) = xint0*tinv
                           yin(jn) = yint0*tinv
                           zin(jn) = zint0*tinv
                           xdip(jn) = xintx*tinv
                           ydip(jn) = yinty*tinv
                           zdip(jn) = zintz*tinv
 90                     continue
 100                 continue
                     call oneld(xin,yin,zin,xd,yd,zd,ai,lit,ljt,1,ndim)
                     call oneld(xdip,ydip,zdip,xdipd,ydipd,zdipd,ai,lit,
     +                          ljt,1,ndim)
c     ----- calculate derivatives of dipole matrix -----
                     n = 0
                     do 130 i = mini , maxi
                        in = loci + i
                        do 120 j = minj , maxj
                           n = n + 1
                           jn = locj + j
                           nn = iky(in) + jn
                           if (jn.gt.in) nn = iky(jn) + in
                           ioff = 0
                           do 110 icp = 1 , ndenin
                              dum = dd(nn+ioff)*dij(n)
                              dum = dum + dum
                              nnx = ijx(n)
                              nny = ijy(n)
                              nnz = ijz(n)
                              dtemp(icp,1,1) = dtemp(icp,1,1)
     +                           + dum*xdipd(nnx)*yin(nny)*zin(nnz)
                              dtemp(icp,2,2) = dtemp(icp,2,2)
     +                           + dum*xin(nnx)*ydipd(nny)*zin(nnz)
                              dtemp(icp,3,3) = dtemp(icp,3,3)
     +                           + dum*xin(nnx)*yin(nny)*zdipd(nnz)
                              dtemp(icp,1,2) = dtemp(icp,1,2)
     +                           + dum*xdip(nnx)*yd(nny)*zin(nnz)
                              dtemp(icp,2,1) = dtemp(icp,2,1)
     +                           + dum*xd(nnx)*ydip(nny)*zin(nnz)
                              dtemp(icp,1,3) = dtemp(icp,1,3)
     +                           + dum*xdip(nnx)*yin(nny)*zd(nnz)
                              dtemp(icp,3,1) = dtemp(icp,3,1)
     +                           + dum*xd(nnx)*yin(nny)*zdip(nnz)
                              dtemp(icp,2,3) = dtemp(icp,2,3)
     +                           + dum*xin(nnx)*ydip(nny)*zd(nnz)
                              dtemp(icp,3,2) = dtemp(icp,3,2)
     +                           + dum*xin(nnx)*yd(nny)*zdip(nnz)
                              ioff = ioff + nx
 110                       continue
 120                    continue
 130                 continue
                  end if
 140           continue
 150        continue
            do 160 icp = 1 , ndenin
               tempol(icp,1,iat*3-2) = tempol(icp,1,iat*3-2)
     +                                 - dtemp(icp,1,1)
               tempol(icp,2,iat*3-1) = tempol(icp,2,iat*3-1)
     +                                 - dtemp(icp,2,2)
               tempol(icp,3,iat*3) = tempol(icp,3,iat*3)
     +                               - dtemp(icp,3,3)
               tempol(icp,1,iat*3-1) = tempol(icp,1,iat*3-1)
     +                                 - dtemp(icp,1,2)
               tempol(icp,1,iat*3) = tempol(icp,1,iat*3)
     +                               - dtemp(icp,1,3)
               tempol(icp,2,iat*3-2) = tempol(icp,2,iat*3-2)
     +                                 - dtemp(icp,2,1)
               tempol(icp,3,iat*3-2) = tempol(icp,3,iat*3-2)
     +                                 - dtemp(icp,3,1)
               tempol(icp,2,iat*3) = tempol(icp,2,iat*3)
     +                               - dtemp(icp,2,3)
               tempol(icp,3,iat*3-1) = tempol(icp,3,iat*3-1)
     +                                 - dtemp(icp,3,2)
 160        continue
 170     continue
 180  continue
c
      do 210 naa = 1 , 3
         do 200 nbb = 1 , 3
            do 190 ncc = 1 , nat3
               pold(naa,nbb,ncc) = tempol(naa,nbb,ncc)
     +                             + tempol(nbb,naa,ncc)
 190        continue
 200     continue
 210  continue
      call pdrres(dd,tempol,pold)
      return
6020  format(/
     +1x,'***************************************',
     +1x,'too great a restriction used in 4-index',
     +1x,'for polarizability derivatives use'/,
     +1x,'restrict 99 or restrict 4 only'/,
     +1x,'***************************************')
      end
      subroutine pdrsym(iso,pold,qudd,qd,ict,ndim,nshels)
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
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/scrtch/ptr(3,144),dtr(6,288)
      dimension map(6)
      dimension qudd(6,3,ndim),qd(6,3,ndim),ict(ndim,48)
      dimension pold(3,3,*),iso(nshels,*)
      data map/1,4,2,5,6,3/
      if (nt.eq.1) return
      np = 0
      do 50 i = 1 , 3
         do 40 j = 1 , i
            np = np + 1
            n3 = 0
            do 30 k = 1 , nat
               do 20 l = 1 , 3
                  n3 = n3 + 1
                  qudd(map(np),l,k) = pold(i,j,n3)
 20            continue
 30         continue
 40      continue
 50   continue
      zero = 0.0d0
c     one = 1.0d0
      do 80 nd = 1 , 6
         do 70 ncc = 1 , 3
            do 60 n = 1 , nat
               qd(nd,ncc,n) = zero
 60         continue
 70      continue
 80   continue
c
c     ----- set transformation table: atoms versus symmetry operations.
c
      do 110 ii = 1 , nshell
         ic = katom(ii)
         do 100 it = 1 , nt
            id = iso(ii,it)
            ict(ic,it) = katom(id)
 100     continue
 110  continue
c
c
      do 170 n = 1 , nat
         do 160 nd = 1 , 6
            do 150 ncc = 1 , 3
               do 140 nop = 1 , nt
                  icnu = ict(n,nop)
                  nnn = 6*(nop-1)
                  nn = 3*(nop-1)
                  do 130 ndd = 1 , 6
                     do 120 ncd = 1 , 3
                        qd(nd,ncc,n) = qd(nd,ncc,n) + qudd(ndd,ncd,icnu)
     +                                 *dtr(ndd,nnn+nd)*ptr(ncd,nn+ncc)
 120                 continue
 130              continue
 140           continue
 150        continue
 160     continue
 170  continue
      dum = dble(nt)
      do 200 n = 1 , nat
         do 190 i = 1 , 6
            do 180 j = 1 , 3
               qudd(i,j,n) = qd(i,j,n)/dum
 180        continue
 190     continue
 200  continue
      np = 0
      do 240 i = 1 , 3
         do 230 j = 1 , i
            np = np + 1
            n3 = 0
            do 220 k = 1 , nat
               do 210 l = 1 , 3
                  n3 = n3 + 1
                  pold(i,j,n3) = qudd(map(np),l,k)
                  pold(j,i,n3) = pold(i,j,n3)
 210           continue
 220        continue
 230     continue
 240  continue
      return
      end
      subroutine pder00(iso,beta,gamm,fa,fb,fc,sc,x,w,e,ua,ub,
     +  qudd,qd,ict,nat,
     +  ltri,norb,npert,nat3,nshels)
c
      implicit real*8  (a-h,o-z)
      dimension x(ltri),fa(ltri),fb(ltri),fc(ltri),e(norb),
     1   beta(npert,npert,nat3),ua(norb,norb),ub(norb,norb)
      dimension gamm(npert,npert,nat3)
      dimension qudd(6,3,nat),qd(6,3,nat),ict(nat,48)
      dimension sc(ltri),w(ltri),iso(nshels,*)
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
      ind(i,j) = iky(max(i,j)) + min(i,j)
c
      lenb = lensec(ltri)
c
      call secget(isect(9),9,isec9)
      call rdedx(e,lds(isect(9)),isec9,ifild)
      call vclr(beta,1,npert*npert*nat3)
      call vclr(ua,1,ncoorb*ncoorb)
c
c      perturbed fock matrices from section 70 thru 72 of dumpfile
c
      do 140 naa = 1 , npert
         iseca1 = naa + 69
         call secget(isect(iseca1),iseca1,iblf)
         call rdedx(fa,ltri,iblf,ifild)
         do 30 i = nocc + 1 , ncoorb
            do 20 j = 1 , nocc
               diff = e(j) - e(i)
               if (dabs(diff).gt.1.0d-10) then
                  ua(i,j) = fa(ind(i,j))/diff
                  ua(j,i) = -ua(i,j)
               end if
 20         continue
 30      continue
         call vclr(ub,1,ncoorb*ncoorb)
         do 130 nbb = 1 , npert
            isecb1 = nbb + 69
            call secget(isect(isecb1),isecb1,iblf)
            call rdedx(fb,ltri,iblf,ifild)
            do 50 i = nocc + 1 , ncoorb
               do 40 j = 1 , nocc
                  diff = e(j) - e(i)
                  if (dabs(diff).gt.1.0d-10) then
                     ub(i,j) = fb(ind(i,j))/diff
                     ub(j,i) = -ub(i,j)
                  end if
 40            continue
 50         continue
            call vclr(w,1,ltri)
            call vclr(x,1,ltri)
            do 80 i = 1 , nocc
               do 70 j = 1 , i
                  ij = iky(i) + j
                  do 60 k = nocc + 1 , ncoorb
                     tem = ua(i,k)*ub(j,k) + ua(j,k)*ub(i,k)
                     x(ij) = x(ij) - tem
                     w(ij) = w(ij) - 0.5d0*(e(i)+e(j))*tem
 60               continue
 70            continue
 80         continue
            do 110 i = nocc + 1 , ncoorb
               do 100 j = 1 , i
                  ij = iky(i) + j
                  do 90 k = 1 , nocc
                     tem = ua(i,k)*ub(j,k) + ua(j,k)*ub(i,k)
                     x(ij) = x(ij) + tem
                     w(ij) = w(ij) + e(k)*tem
 90               continue
 100           continue
 110        continue
            ibf = iochf(16)
            ibs = iochf(14)
            do 120 ncc = 1 , nat3
               call rdedx(fc,ltri,ibf,ifockf)
               call rdedx(sc,ltri,ibs,ifockf)
               abc = -tracep(fc,x,ncoorb) + tracep(sc,w,ncoorb)
               abc = abc + abc
               beta(naa,nbb,ncc) = beta(naa,nbb,ncc) + abc
               ibf = ibf + lenb
               ibs = ibs + lenb
 120        continue
 130     continue
 140  continue
      call pdrsym(iso,beta,qudd,qd,ict,nat,nshels)
      do 170 i = 1 , npert
         do 160 j = 1 , npert
            do 150 k = 1 , nat3
               gamm(i,j,k) = gamm(i,j,k) + beta(i,j,k)
 150        continue
 160     continue
 170  continue
      return
      end
      subroutine pder01(iso,beta,gamm,fa,fb,fc,sb,x,w,e,ua,ub,
     +     qudd,qd,ict,nat,
     +     ltri,norb,npert,nat3,nshels)
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
      dimension x(ltri),fa(ltri),fb(ltri),fc(ltri),e(norb),
     1   beta(npert,npert,nat3),ua(norb,norb),ub(norb,norb)
      dimension gamm(npert,npert,nat3)
      dimension sb(ltri),w(ltri),iso(nshels,*)
      dimension qudd(6,3,nat),qd(6,3,nat),ict(48,nat)
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
      lenb = lensec(ltri)
c
      call secget(isect(9),9,isec9)
      call rdedx(e,lds(isect(9)),isec9,ifild)
      call vclr(beta,1,npert*npert*nat3)
      call vclr(ua,1,ncoorb*ncoorb)
c
c      perturbed fock matrices from section 2 thru 4 of fockfile
c
      do 160 naa = 1 , npert
         iseca1 = naa + 69
         call secget(isect(iseca1),iseca1,iblf)
         call rdedx(fa,ltri,iblf,ifild)
         do 30 i = nocc + 1 , ncoorb
            do 20 j = 1 , nocc
               diff = e(j) - e(i)
               if (dabs(diff).gt.1.0d-10) then
                  ua(i,j) = fa(ind(i,j))/diff
                  ua(j,i) = -ua(i,j)
               end if
 20         continue
 30      continue
         call vclr(ub,1,ncoorb*ncoorb)
         ibf = iochf(16)
         ibs = iochf(14)
         do 150 nbb = 1 , nat3
            call rdedx(fb,ltri,ibf,ifockf)
            call rdedx(sb,ltri,ibs,ifockf)
            do 50 i = 1 , ncoorb
               do 40 j = 1 , i
                  ub(i,j) = -0.5d0*sb(ind(i,j))
                  ub(j,i) = ub(i,j)
 40            continue
 50         continue
            do 70 i = nocc + 1 , ncoorb
               do 60 j = 1 , nocc
                  diff = e(j) - e(i)
                  ij = ind(i,j)
                  if (dabs(diff).gt.1.0d-10) then
                     ub(i,j) = (fb(ij)-e(j)*sb(ij))/diff
                     ub(j,i) = -(ub(i,j)+sb(ij))
                  end if
 60            continue
 70         continue
            call vclr(w,1,ltri)
            call vclr(x,1,ltri)
            do 100 i = 1 , nocc
               do 90 j = 1 , i
                  ij = iky(i) + j
                  do 80 k = 1 , ncoorb
                     tem = ua(i,k)*ub(j,k) + ua(j,k)*ub(i,k)
                     x(ij) = x(ij) - tem
 80               continue
 90            continue
 100        continue
            do 130 i = 1 , ncoorb
               do 120 j = 1 , i
                  ij = iky(i) + j
                  do 110 k = 1 , nocc
                     tem = ua(i,k)*ub(j,k) + ua(j,k)*ub(i,k)
                     x(ij) = x(ij) + tem
 110              continue
 120           continue
 130        continue
            do 140 ncc = 1 , npert
               iseca1 = ncc + 69
               call secget(isect(iseca1),iseca1,iblf)
               call rdedx(fc,ltri,iblf,ifild)
               abc = -tracep(fc,x,ncoorb)
               abc = abc + abc
               beta(naa,ncc,nbb) = beta(naa,ncc,nbb) + abc
               beta(ncc,naa,nbb) = beta(ncc,naa,nbb) + abc
 140        continue
            ibf = ibf + lenb
            ibs = ibs + lenb
 150     continue
 160  continue
      call pdrsym(iso,beta,qudd,qd,ict,nat,nshels)
      do 190 i = 1 , npert
         do 180 j = 1 , npert
            do 170 k = 1 , nat3
               gamm(i,j,k) = gamm(i,j,k) + beta(i,j,k)
 170        continue
 180     continue
 190  continue
      return
      end
      subroutine pder02(iso,beta,gamm,da,fb,fc,vec,work,
     +    qudd,qd,ict,nat,
     +    ltri,norb,nbasis,npert,nat3,nshels)
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
      character*2 cx,cxi,cxj
      dimension fc(ltri),iso(nshels,*)
      dimension work(nbasis),vec(nbasis,norb)
      dimension da(ltri),fb(ltri),
     1   beta(npert,npert,nat3)
      dimension gamm(npert,npert,nat3)
      dimension qudd(6,3,nat),qd(6,3,nat),ict(48,nat)
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
       common/specal/ndenin,iflden,iblden
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
       dimension cx(9)
       dimension cxi(3),cxj(3),buff(3),kbuff(3)
       data cx/'x','y','z','xx','yy','zz','xy','xz','yz'/
c
      ind(i,j) = iky(max(i,j)) + min(i,j)
c
      m = 0
      call secget(isect(8),m,iblok)
      call rdedx(vec,norb*nbasis,iblok+mvadd,ifild)
c     lenb = lensec(ltri)
c
c
      call vclr(beta,1,npert*npert*nat3)
c
c      perturbed fock matrices from section 2 thru 4 of fockfile
      call search(iblden,ifild)
c
      do 60 naa = 1 , npert
         call reads(da,ltri,ifild)
         ibf = iochf(18)
         call search(ibf,ifockf)
         do 50 nbb = 1 , ndenin
            do 40 ncc = 1 , nat3
               call reads(fc,ltri,ifockf)
               call qhq2(fb,fc,vec,work,iky,ncoorb,nbasis,nbasis)
               abc = 0.0d0
               do 30 i = nocc + 1 , ncoorb
                  do 20 j = 1 , nocc
                     ij = ind(i,j)
                     abc = abc - da(ij)*fb(ij)
 20               continue
 30            continue
               beta(naa,nbb,ncc) = beta(naa,nbb,ncc) + abc
               beta(nbb,naa,ncc) = beta(nbb,naa,ncc) + abc
 40         continue
 50      continue
 60   continue
      call pdrsym(iso,beta,qudd,qd,ict,nat,nshels)
      write (iwr,6010)
      if (opunch(4)) write (ipu,6020)
      if (opunch(4)) write (ipu,6030) title
      l = 0
      do 90 i = 1 , npert
         do 80 j = 1 , ndenin
            do 70 k = 1 , nat3
               gamm(i,j,k) = gamm(i,j,k) + beta(i,j,k)
               if (opunch(4)) write (ipu,6040) cx(i) , cx(j) , k ,
     +                               gamm(i,j,k)
               l = l + 1
               cxi(l) = cx(i)
               cxj(l) = cx(j)
               kbuff(l) = k
               buff(l) = gamm(i,j,k)
               if(l.lt.3) go to 70
               write (iwr,6050) (cxi(l) , cxj(l) , kbuff(l) , 
     +         buff(l),l =1,3)
               l = 0
 70         continue
 80      continue
 90   continue
        if(l.ne.0)
     +  write (iwr,6050) (cxi(m) , cxj(m) , kbuff(m) , 
     +         buff(m),m =1,l)
      leng = lensec(npert*npert*nat3)
      lds(isect(55)) = npert*npert*nat3
      ityp = 55
      call secput(isect(55),ityp,leng,iblok)
      call wrt3(gamm,npert*npert*nat3,iblok,ifild)
      call revind
      return
 6010 format (/
     + 20x,'****************************************'/
     + 20x,'* cartesian polarizability derivatives *'/
     + 20x,'*           (atomic units)             *'/
     + 20x,'****************************************'/)
 6020 format (1x,'cartesian polarizability derivatives')
 6030 format (1x,10a8)
 6040 format (1x,2(a2,6x),i8,d20.12)
 6050 format (3(10x,2a3,i4,f15.7))
      end
      subroutine pdrres(q,beta,gamm)
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
      common/maxlen/maxq
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/scrtch/ptr(3,144),dtr(6,288)
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
      character*10 charwall
      parameter (maxn3=maxat*3)
      dimension beta(3,3,maxn3),gamm(3,3,maxn3)
      dimension q(*)
c
      ltri = iky(ncoorb+1)
      i0 = 1 + nw196(5)
      i1 = ltri + i0
      i2 = i1 + ltri
      i3 = i2 + ltri
      i4 = i3 + ltri
      i5 = i4 + ltri
      i6 = i5 + ltri
      i7 = i6 + ncoorb
      i8 = i7 + ncoorb*ncoorb
      i9 = i8 + ncoorb*ncoorb
      i10 = i9 + 18*nat
      i11 = i10 + 18*nat
      last = i11 + 48*nat
      if (last.gt.maxq) call caserr('insufficient core')
      nat3 = nat*3
c
c     ------------------------------------------------------
c     read in transformation matrices for s,p,d,f functions.
c     ------------------------------------------------------
c
      call rdedx(ptr,nw196(1),ibl196(1),ifild)
      call rdedx(dtr,nw196(2),ibl196(2),ifild)
      call rdedx(q,nw196(5),ibl196(5),ifild)
c
      call pder00(q,beta,gamm,q(i0),q(i1),q(i2),q(i3),q(i4),q(i5),q(i6)
     +  ,q(i7),q(i8),q(i9),q(i10),q(i11),nat,ltri,ncoorb,3,nat3,nshell)
      call pder01(q,beta,gamm,q(i0),q(i1),q(i2),q(i3),q(i4),q(i5),q(i6)
     +  ,q(i7),q(i8),q(i9),q(i10),q(i11),nat,ltri,ncoorb,3,nat3,nshell)
      call pder02(q,beta,gamm,q(i0),q(i1),q(i2),q(i3),q(i5),q(i9),
     +  q(i10),q(i11),nat,ltri,ncoorb,ncoorb,3,nat3,nshell)
      write(iwr,6010) cpulft(1),charwall()
 6010  format(/1x,'assembly of polarisability derivatives complete at',
     + 1x,f8.2,' seconds',a10,' wall')
      call timit(3)
      return
      end
      subroutine tdchfa (vec,eig,tm,wkspce,prop,ndim,nsol)
      implicit real*8  (a-h,o-z)
      logical full
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
      character *8 pn
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
      dimension prop(np,np),sum(6)
      dimension vec(ndim,nsol),eig(nsol),tm(np,nsol),wkspce(ndim)
      dimension pn(50),iang(50)
c
      npdim = max(np,9)
      maxper = 50
      nn = 0
      do 20 i = 1 , maxper
         if (.not.(opskip(i))) then
            nn = nn + 1
            pn(nn) = pnames(i)
            iang(nn) = ipang(i)
         end if
 20   continue
c
c
      full = nsol.eq.ndim
      npolep = npole*np
      nsolp = nsol*np
c     lenblk = lensec(ndim)
c
c...  compute transition moments
      do 40 i = 1 , nsol
         eig(i) = dsqrt(eig(i))
         zz = dsqrt(eig(i)*2.0d0)
         do 30 j = 1 , ndim
            vec(j,i) = vec(j,i)*zz
 30      continue
 40   continue
      call search(iblks,ifils)
      call vclr(tm,1,nsolp)
      do 50 ko = 1 , np
         call reads(wkspce,ndim,ifils)
c     punch*,wkspce
         call mxmb(vec,ndim,1,wkspce,1,ndim,tm(ko,1),np,1,nsol,ndim,1)
 50   continue
      write (iwr,6010)
      do 60 i = 1 , npole
         write (iwr,6020) i
         call prfreq(2,eig(i))
         write (iwr,6030) (tm(ko,i),ko=1,np)
 60   continue
      npdim = max(np,9)
      call secget(isect(52),52,iblkj)
      iblkj = iblkj + 2
      call wrt3(eig,npole,iblkj,ifild)
      call wrt3s(tm,npolep,ifild)
c
      if (full) then
c...  polarisability
         if (nfreq.ne.0) then
            write (iwr,6040)
            do 100 ifreq = 1 , nfreq
               do 90 jo = 1 , np
                  do 80 ko = 1 , jo
                     zz = 0.0d0
                     do 70 l = 1 , ndim
                        zz = zz + tm(jo,l)*tm(ko,l)*eig(l)
     +                       /(eig(l)**2-freq(ifreq))
 70                  continue
                     prop(jo,ko) = zz*2.0d0
                     prop(ko,jo) = zz*2.0d0
 80               continue
 90            continue
               call prpol0(prop,npdim)
 100        continue
         end if
c
c...  sums
         write (iwr,6050)
         do 130 ko = 1 , np
            call vclr(sum,1,6)
            do 120 j = 1 , ndim
               omeg = eig(j)**2
               f = 2.0d0*eig(j)*tm(ko,j)**2
               do 110 i = 1 , 6
                  sum(i) = sum(i) + f
                  f = f/omeg
 110           continue
 120        continue
            write (iwr,6060) pn(ko) , sum
 130     continue
c
c...  dispersion integrals
         write (iwr,6070)
         ndisp = ndim
         do 200 iordc = nc6min , nc6max
            if (oc6(iordc)) then
               iord = iordc - 2
               write (iwr,6080) iord
               do 190 io = 1 , np
                  do 180 jo = 1 , io
                     do 170 ko = 1 , io
                        lomax = ko
                        if (io.eq.ko) lomax = jo
                        do 160 lo = 1 , lomax
                           iiii = iang(io) + iang(jo) + iang(ko)
     +                            + iang(lo)
                           if (iiii.eq.iord) then
                              zz = 0.0d0
                              do 150 i = 1 , ndisp
                                 do 140 j = 1 , ndisp
                                    zz = zz + tm(io,i)*tm(jo,i)*tm(ko,j)
     +                                 *tm(lo,j)/(eig(i)+eig(j))
 140                             continue
 150                          continue
                              if (dabs(zz).gt.1.0d-8) write (iwr,6090)
     +                            pn(io) , pn(jo) , pn(ko) , pn(lo) , zz
                           end if
 160                    continue
 170                 continue
 180              continue
 190           continue
            end if
 200     continue
      end if
      return
 6010 format (/
     +  20x,'********************************'/
     +  20x,'*** analysis of rpa solution ***'/
     +  20x,'********************************'//
     + 1x,'transition energies and matrix elements:'//
     + 1x,'========================================')
 6020 format (1x,'transition',i4,5x,'transition energy:')
 6030 format (' transition matrix elements:',3f20.10,10(/28x,3f20.10))
 6040 format (1x,'polarizabilities from rpa solution'//)
 6050 format (/
     + 30x,'========================'/
     + 30x,'= energy weighted sums ='/
     + 30x,'========================'//
     + 15x,'s(0)',16x,'s(-2)',15x,
     +        's(-4)',15x,'s(-6)',15x,'s(-8)',15x,'s(-10)'/)
 6060 format (1x,a8,1x,6f20.7)
 6070 format (/
     + 10x,'================================================='/
     + 10x,'= dispersion energy integrals from rpa solution ='/
     + 10x,'================================================='/)
 6080 format (//
     + 20x,'****************************************'/
     + 20x,'* terms for total angular momentum',i4,' *'/
     + 20x,'****************************************'//)
 6090 format (4x,a8,1x,a8,4x,a8,1x,a8,f18.8)
      end
      subroutine tdchf0 (eps,b,wkspce,eig,hr,hi,ndim)
c
      implicit real*8  (a-h,o-z)
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
c
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
      dimension eps(ndim),b(ndim,np),wkspce(ndim,2),eig(ndim)
     1,    hr(ndim,ndim),hi(ndim,ndim)
c
      call search(jblk(1),nofile(1))
      call tdchhm(hr,eps,ndim,nofile(1))
      call tdchhm(hi,eps,ndim,nofile(1))
c
      ifail = 0
      call abzez(ndim,hi,hr,hi,eig,wkspce,ifail)
c
c     this may fail if np.gt.ndim
      call tdchfa(hi,eig,b,wkspce,hr,ndim,ndim)
      return
      end
      subroutine chficl(a,maxa)
      implicit real*8  (a-h,o-z)
      dimension a(*)
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/craypk/labs(1360)
      common/blkin/g(510),nword
      character*10 charwall
c
c      routine chficl constructs the 'a-matrix' for
c      a closed shell system (imaginary perturbations)
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
      iblkhi = iposun(nofile(1))
c
      nst = 0
      nmin = 0
      nfin = mn
      do 20 nn = 1 , mn
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
                  kk2 = (kk+kk) + (kk+kk)
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
                              na = (i-nocca-1)*nocca + l
                              nb = (k-nocca-1)*nocca + j
                              if (na.lt.nb) then
                                 nswop = na
                                 na = nb
                                 nb = nswop
                              end if
                              if (na.gt.nst .and. na.le.nfin) then
                                 ntri = na*(na-1)/2 + nb - nmin
                                 a(ntri) = a(ntri) + gg
                              end if
                           end if
                        else if (k.le.nocca) then
c
c     type (xx/oo)
c
                           na = (i-nocca-1)*nocca + k
                           nb = (j-nocca-1)*nocca + l
                           if (na.lt.nb) then
                              nswop = na
                              na = nb
                              nb = nswop
                           end if
c
c     only that part of triangle between a(nst+1,1) and a(nfin,nfin)
c     is constructed in this pass
c
                           if (na.gt.nst .and. na.le.nfin) then
                              ntri = na*(na-1)/2 + nb - nmin
                              a(ntri) = a(ntri) - gg
                           end if
                           if (i.ne.j .and. k.ne.l) then
                              na = (i-nocca-1)*nocca + l
                              nb = (j-nocca-1)*nocca + k
                              if (na.lt.nb) then
                                 nswop = na
                                 na = nb
                                 nb = nswop
                              end if
                              if (na.gt.nst .and. na.le.nfin) then
                                 ntri = na*(na-1)/2 + nb - nmin
                                 a(ntri) = a(ntri) - gg
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
      call wamat(a,mn,nofile(1),nst,nfin,odebug(2),iwr)
      if (nfin.eq.mn) then
c
         write (iwr,6010) cpulft(1),charwall()
         return
      else
         nst = nfin
         nfin = mn
         nmin = nst*(nst+1)/2
         nstp1 = nst + 1
         do 70 nn = nstp1 , mn
            last = nn*(nn+1)/2 - nmin
            if (last.gt.maxa) then
               nfin = nn - 1
               go to 30
            end if
 70      continue
      end if
      go to 30
 6010 format (/1x,'construction of magnetic hessian complete at',
     + f8.2,' seconds',a10,' wall')
      end
      subroutine trksum (eps,v1,v2,npd,ndim)
c
      implicit real*8  (a-h,o-z)
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
      dimension v1(ndim,npd),v2(ndim,npd),eps(ndim)
c
      write (iwr,6010)
c
      call rdedv(v1,ndim,npd,iblks,ifils)
      call search(iblkhi,nofile(1))
      call tdchhv(v1,v2,eps,ndim,npd,nofile(1))
      ko = 0
      do 20 i = 1 , 3
         if (.not.(opskip(i))) then
            ko = ko + 1
            sum = 4.0d0*ddot(ndim,v1(1,ko),1,v2(1,ko),1)
            percnt = dble(nocca*2)
            percnt = 100.0d0*(sum-percnt)/percnt
            write (iwr,6020) pnames(i) , sum , percnt
         end if
 20   continue
      return
 6010 format (/
     + 20x,'*********************************'/
     + 20x,'* thomas-reich-kuhn dipole sums *'/
     + 20x,'*********************************'/)
 6020 format (1x,'sum rule for ',a8,'  =',f15.7,5x,'(',f8.4,'%)')
      end
      subroutine magasm(b,u,dia,ndim)
c
c     assembles magnetisability
c
      implicit real*8  (a-h,o-z)
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
      dimension prop(6,6),dia(6,6)
      dimension b(ndim,np),u(ndim,np)
c
c     read in pertubation vectors in b
c
      call rdedv(b,ndim,np,iblks,ifils)
c
c     read the solution vectors in u
c
      call rdedvs(u,ndim,np,ifils)
c
      write (iwr,6010)
      do 20 i = 1 , 3
         write (iwr,6020) (dia(i,j),j=1,3)
 20   continue
      do 40 i = 1 , 3
         do 30 j = 1 , 3
            prop(i,j) = -4.0d0*ddot(mn,u(1,i),1,b(1,j),1)
 30      continue
 40   continue
      write (iwr,6030)
      do 50 i = 1 , 3
         write (iwr,6020) (prop(i,j),j=1,3)
 50   continue
      do 70 i = 1 , 3
         do 60 j = 1 , 3
            prop(i,j) = prop(i,j) + dia(i,j)
 60      continue
 70   continue
      write (iwr,6040)
      do 80 i = 1 , 3
         write (iwr,6020) (prop(i,j),j=1,3)
 80   continue
      do 100 i = 1 , 3
         do 90 j = 1 , 3
            prop(i,j) = 7.8910322d0*prop(i,j)
 90      continue
 100  continue
      write (iwr,6050)
      do 110 i = 1 , 3
         write (iwr,6020) (prop(i,j),j=1,3)
 110  continue
      return
 6010 format (/
     + 20x,'*************************************'/
     + 20x,'* diamagnetic susceptibility tensor *'/
     + 20x,'*           (atomic units)          *'/
     + 20x,'*************************************'/)
 6020 format (/10x,3f18.7)
 6030 format (/
     + 20x,'**************************************'/
     + 20x,'* paramagnetic susceptibility tensor *'/
     + 20x,'*            (atomic units)          *'/
     + 20x,'**************************************'/)
 6040 format (/
     + 20x,'************************'/
     + 20x,'* total susceptibility *'/
     + 20x,'*    (atomic units)    *'/
     + 20x,'************************'/)
 6050 format (/
     + 20x,'*******************************************'/
     + 20x,'* in s.i. units - 10**-29 joule tesla**-2 *'/
     + 20x,'*******************************************'/)
      end
      subroutine vcdovl(vdcd,bd,nat3,ncomp)
c
c     overlap contribution to vcd
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
      logical out,norm
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
      common/junk/xyz(3,maxat),
     +  xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
     + ,cx,cy,cz,idum(225),ijx(225),ijy(225),ijz(225),
     + xin(25),yin(25),zin(25),xd(25),yd(25),zd(25),dij(225),
     + sx(225),sy(225),sz(225)
      dimension sss(225,3)
      dimension vdcd(nat3,ncomp),bd(num,num,ncomp)
      equivalence (sx(1),sss(1,1))
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      data ndim/5/
      data zero,one /0.0d0,1.0d0/
c
c     ----- calculate derivatives if the overlap matrix -----
c
      tol = 2.30258d0*itol
      out = odebug(1)
      norm = normf.ne.1 .or. normp.ne.1
c     ----- ishell
      do 130 ii = 1 , nshell
         iat = katom(ii)
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         lit1 = lit + 1
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 120 jj = 1 , nshell
            jat = katom(jj)
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c     ----- prepare indices for pairs of (i,j) functions
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,.false.,
     +           ndim,1,1)
            do 20 i = 1 , ij
               sx(i) = zero
               sy(i) = zero
               sz(i) = zero
 20         continue
c     ----- i primitive
            do 70 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               do 60 jg = j1 , j2
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                               csj,cpj,cdj,cfj,cgj,
     +                           mini,maxi,minj,maxj,.false.,.false.,
     +                           norm)
c     ----- overlap
                     t = dsqrt(aa)
                     tinv = one/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     in = -ndim
                     do 40 i = 1 , lit1
                        in = in + ndim
                        ni = i
                        do 30 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call stvin2()
                           xin(jn) = xint*tinv
                           yin(jn) = yint*tinv
                           zin(jn) = zint*tinv
 30                     continue
 40                  continue
c
                     call oneld(xin,yin,zin,xd,yd,zd,ai,lit,ljt,1,ndim)
c
                     do 50 i = 1 , ij
                        mx = ijx(i)
                        my = ijy(i)
                        mz = ijz(i)
                        sx(i) = sx(i) + dij(i)*xd(mx)*yin(my)*zin(mz)
                        sy(i) = sy(i) + dij(i)*xin(mx)*yd(my)*zin(mz)
                        sz(i) = sz(i) + dij(i)*xin(mx)*yin(my)*zd(mz)
 50                  continue
                  end if
c     ----- end of primitive loops -----
 60            continue
 70         continue
c     ----- calculate derivatives of overlap matrix -----
            iatl = (iat-1)*3
            n = 0
            do 110 i = mini , maxi
               in = loci + i
               do 100 j = minj , maxj
                  n = n + 1
                  jn = locj + j
                  do 90 k = 1 , 3
                     do 80 l = 1 , ncomp
                        vdcd(iatl+k,l) = vdcd(iatl+k,l) + sss(n,k)
     +                                  *bd(in,jn,l)
                        if (out) write (iwr,6010) ii , jj , in , jn , 
     +                               k , l , sss(n,k) , bd(in,jn,l) ,
     +                               vdcd(iatl+k,l)
 80                  continue
 90               continue
 100           continue
 110        continue
 120     continue
 130  continue
      return
 6010 format (1x,6i10,3f20.10)
      end
      subroutine vcdpar(u,da,db,temp,denb,ndim)
c                vcdpar(u,da,db,temp,denb,mapnr,ndim)
c
c     v.c.d. parameters
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
      parameter (maxn3=maxat*3)
      logical exist
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
      common/junk/xyz(3,maxat),xint(13),nij(2),cxyz(3),ijx(4,225),
     + xin(6,25),dij(4,225),
     + gg(maxn3,6)
c
      dimension gga(maxn3,3)
      equivalence (gga(1,1),gg(1,4))
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      dimension u(ndim,*),da(*),db(*),denb(num,num,3)
      dimension temp(*)
c     dimension mapnr(*)
c
c
      nat3 = nat*3
c
c     read the solution vectors in u
c
      call secloc(isect(67),exist,iblok)
      if (.not.exist) call caserr(
     +  'perturbed wavefunction (magnetic) not present')
      call search(iblok,ifild)
      do 20 i = 1 , 3
         call reads(u(1,i),mn,ifild)
 20   continue
      if (odebug(1)) write (iwr,6010) (u(i,1),i=1,mn)
      if (odebug(1)) write (iwr,6010) (u(i,2),i=1,mn)
      if (odebug(1)) write (iwr,6010) (u(i,3),i=1,mn)
c
      do 40 i = 1 , 3
         do 30 j = 1 , nat3
            gga(j,i) = 0.0d0
            gg(j,i) = 0.0d0
 30      continue
 40   continue
c
c     lennew = ncoorb*(ncoorb+1)/2
c
c     the orbital part is same as g-tensor
c
      call secloc(isect(66),exist,iblk66)
      if (.not.exist) call caserr(
     +      'perturbed wavefunction (nuc. displac.) not present')
      call search(iblk66,ifild)
      do 50 j = 1 , nat3
         call reads(db,mn,ifild)
         if (odebug(1)) write (iwr,6010) (db(i),i=1,mn)
         gg(j,1) = 2.0d0*ddot(mn,db,1,u(1,1),1)
         gg(j,2) = 2.0d0*ddot(mn,db,1,u(1,2),1)
         gg(j,3) = 2.0d0*ddot(mn,db,1,u(1,3),1)
 50   continue
c
c---------------------------------------------------------------
c
c
c     overlap part ---- requires derivative of overlap and
c     a perturbed denisty matrix
c
c     get mo's
c
      mtyp = 0
      call secget(isect(8),mtyp,iblok)
      call rdedx(da,num*ncoorb,iblok+mvadd,ifild)
c
c     construct density
c
      call vcdden(u,da,denb,temp,num,nocca,ncoorb,mn,3)
c
c     overlap derivative
c
      call vcdovl(gga,denb,maxn3,3)
c------------------------------------------------------------------
      do 70 j = 1 , nat3
         do 60 i = 1 , 3
            gg(j,i) = gg(j,i) + gga(j,i)
 60      continue
 70   continue
      write (iwr,6020)
      do 90 jj = 1 , nat
         write (iwr,6040)
         do 80 i = 1 , 3
            j = (jj-1)*3 + i
            write (iwr,6030) (gg(j,k),k=1,3)
 80      continue
 90   continue
c----------------convert units---------------------
      do 110 j = 1 , nat3
         do 100 i = 1 , 3
            gga(j,i) = gg(j,i)/1.2439d5
 100     continue
 110  continue
      write (iwr,6050)
      do 130 jj = 1 , nat
         write (iwr,6040)
         do 120 i = 1 , 3
            j = (jj-1)*3 + i
            write (iwr,6060) (gga(j,k),k=1,3)
 120     continue
 130  continue
c--------------------------nuclear term-------------------
c
      call vcdnuc(gga,maxn3,3,gx)
      write (iwr,6070)
      do 150 jj = 1 , nat
         write (iwr,6040)
         do 140 i = 1 , 3
            j = (jj-1)*3 + i
            write (iwr,6030) (gga(j,k),k=1,3)
 140     continue
 150  continue
c-----------------sum rules------------------------------
      call vcmrul(c,gg,nat,gx)
      m56 = 56
      m2 = 2
      m540 = 6*maxn3
      call secput(isect(56),m56,m2,iblok)
      call wrt3(gg,m540,iblok,ifild)
      lds(isect(56)) = m540
      call revind
c
c-------------------punch option------------------------
c
      if (opunch(1)) then
         write (ipu,6080)
         write (ipu,6090) (title(i),i=1,10)
         write (ipu,6100)
         do 160 i = 1 , nat
            write (ipu,6130) c(1,i) , c(2,i) , c(3,i)
 160     continue
         write (ipu,6110)
         write (ipu,6130) gx , gy , gz
         write (ipu,6120)
         do 170 j = 1 , nat3
            write (ipu,6130) (gg(j,k),k=1,3)
 170     continue
         write (ipu,6140)
         do 180 j = 1 , nat3
            write (ipu,6130) (gga(j,k),k=1,3)
 180     continue
      end if
c
      return
 6010 format (///(1x,5f16.8))
 6020 format (/
     + 20x,'***********************************'/
     + 20x,'* total wavefunction overlap term *'/
     + 20x,'*         (atomic units)          *'/
     + 20x,'***********************************'//
     +        16x,'bx',16x,'by',16x,'bz'/)
 6030 format (10x,3f16.8)
 6040 format (/)
 6050 format (/
     + 20x,'***********************************'/
     + 20x,'* total wavefunction overlap term *'/
     + 20x,'*     units 1/(angstrom*tesla)    *'/
     + 20x,'***********************************'//
     +        16x,'bx',16x,'by',16x,'bz'/)
 6060 format (10x,3e16.8)
 6070 format (/
     + 20x,'******************'/
     + 20x,'*  nuclear term  *'/
     + 20x,'* (atomic units) *'/
     + 20x,'******************'//
     + 16x,'bx',16x,'by',16x,'bz'/)
 6080 format (1x,'vcd data')
 6090 format (1x,10a8)
 6100 format (1x,'current geometry')
 6110 format (1x,'current gauge origin')
 6120 format (/1x,'vcd electronic contribution, atomic units')
 6130 format (1x,3e20.12)
 6140 format (/1x,'vcd nuclear contribution, atomic units')
      end
      subroutine vcdnuc(gg,ndim,ncomp,gauge)
c
c      nuclear term
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
      dimension gg(ndim,ncomp),gauge(3)
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
      do 30 i = 1 , ndim
         do 20 j = 1 , ncomp
            gg(i,j) = 0.0d0
 20      continue
 30   continue
c
      do 70 i = 1 , ncomp
         do 60 j = 1 , 3
            do 50 k = 1 , 3
               skew = etijk(i,j,k)
               if (skew.ne.0.0d0) then
                  do 40 n = 1 , nat
                     gg((n-1)*3+k,i) = gg((n-1)*3+k,i)
     +                                 + skew*(c(j,n)-gauge(j))*czan(n)
     +                                 *0.25d0
 40               continue
               end if
 50         continue
 60      continue
 70   continue
      return
      end
      subroutine vrepdd(u,da,db,temp,denb,ndim)
c                vrepdd(u,da,db,temp,denb,mapnr,ndim)
c
c     velocity representation dipole derivatives
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
      parameter (maxn3=maxat*3)
      logical exist
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
      common/junk/xyz(3,maxat),xint(13),nij(2),cxyz(3),ijx(4,225),
     + xin(6,25),dij(4,225),
     + gg(maxn3,6)
      dimension gga(maxn3,3)
      equivalence (gga(1,1),gg(1,4))
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      dimension u(ndim,*),da(*),db(*),denb(num,num,3)
c     dimension mapnr(*)
      dimension temp(*)
c
c
      nat3 = nat*3
      call secloc(isect(67),exist,iblok)
      if (.not.exist) call caserr(
     +  'first order wavefunction (magnetic) not present')
      call search(iblok,ifild)
c
c     read the solution vectors in u
c
      do 20 i = 1 , 3
         call reads(u(1,i),mn,ifild)
 20   continue
      do 30 i = 1 , 3
         call reads(u(1,i),mn,ifild)
 30   continue
c
c     ( the odd overwriting in above lines is
c       deliberate - it ensures the correct 3 perturbations
c       are in u )
c
      if (odebug(1)) write (iwr,6010) (u(i,1),i=1,mn)
      if (odebug(1)) write (iwr,6010) (u(i,2),i=1,mn)
      if (odebug(1)) write (iwr,6010) (u(i,3),i=1,mn)
c
      do 50 i = 1 , 3
         do 40 j = 1 , nat3
            gga(j,i) = 0.0d0
            gg(j,i) = 0.0d0
 40      continue
 50   continue
c     lennew = ncoorb*(ncoorb+1)/2
c
c     the orbital part is same as g-tensor
c     and the vcd tensor
      call secloc(isect(66),exist,iblk66)
      if (.not.exist) call caserr(
     +            'first order wavefunction (nuc. displac.) not present'
     +            )
      call search(iblk66,ifild)
      do 60 j = 1 , nat3
c     read the perturbed density matrix for nuclear
c     displacements
         call reads(db,mn,ifild)
         if (odebug(1)) write (iwr,6010) (db(i),i=1,mn)
         gg(j,1) = 2.0d0*ddot(mn,db,1,u(1,1),1)
         gg(j,2) = 2.0d0*ddot(mn,db,1,u(1,2),1)
         gg(j,3) = 2.0d0*ddot(mn,db,1,u(1,3),1)
 60   continue
c
c
c     overlap part ---- requires derivative of overlap and
c     a perturbed density matrix . this part is exactly
c     the same as vcd term ( only difference is the u matrix )
c
c     get mo's
c
      mtyp = 0
      call secget(isect(8),mtyp,iblok)
      call rdedx(da,num*ncoorb,iblok+mvadd,ifild)
c
c     construct density
c
      call vcdden(u,da,denb,temp,num,nocca,ncoorb,mn,3)
c
c     overlap derivative
c
      call vcdovl(gga,denb,maxn3,3)
c
c------------------------electronic total---------------
c
      do 80 j = 1 , nat3
         do 70 i = 1 , 3
            gg(j,i) = -2.0d0*(gg(j,i)+gga(j,i))
 70      continue
 80   continue
      write (iwr,6040)
      do 100 jj = 1 , nat
         write (iwr,6030)
         do 90 i = 1 , 3
            j = (jj-1)*3 + i
            write (iwr,6020) (gg(j,k),k=1,3)
 90      continue
 100  continue
c
c ------------------add nuclear term-----------------------
c
      write (iwr,6050)
      do 120 jj = 1 , nat
         write (iwr,6030)
         do 110 i = 1 , 3
            j = (jj-1)*3 + i
            gg(j,i) = gg(j,i) + czan(jj)
            write (iwr,6020) (gg(j,k),k=1,3)
 110     continue
 120  continue
      call vcdrul(c,gg,nat,maxn3)
      if (opunch(1)) then
         write (ipu,6060)
         do 130 i = 1 , nat3
            write (ipu,6070) gg(i,1) , gg(i,2) , gg(i,3)
 130     continue
      end if
c
      return
 6010 format (/(1x,5f16.8))
c
 6020 format (10x,3f16.8)
 6030 format (/)
 6040 format (/
     + 20x,'*******************'/
     + 20x,'* electronic term *'/
     + 20x,'*  (atomic units) *'/
     + 20x,'*******************'//
     + 16x,'px',16x,'py',16x,'pz'/)
 6050 format (/
     + 20x,'************************'/
     + 20x,'* total (atomic units) *'/
     + 20x,'************************'//
     + 16x,'px',16x,'py',16x,'pz'/)
 6060 format (1x,'velocity formalism dipole derivatives')
 6070 format (1x,3e20.12)
      end
      subroutine vcdrul(r,dipd,nat,maxn3)
c
c --------------dipole rotational sum rules----------
c  ( like diprot but with allowance for changed storage )
c
      implicit real*8  (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension r(3,nat),dipd(maxn3,3),rot(9)
      character*2 head
      dimension head(9)
      data head/'xx','xy','xz','yx','yy','yz','zx','zy','zz'/
      n = 0
      do 40 ialp = 1 , 3
         do 30 ibet = 1 , 3
            n = n + 1
            rot(n) = 0.0d0
            do 20 i = 1 , nat
               rot(n) = rot(n) + dipd((i-1)*3+ibet,ialp)
 20         continue
 30      continue
 40   continue
      write (iwr,6010)
      do 100 i=1,9
 100  write (iwr,6030) head(i),rot(i)
      n = 0
      do 90 ibet = 1 , 3
         do 80 idel = 1 , 3
            n = n + 1
            rot(n) = 0.0d0
            do 70 ialp = 1 , 3
               do 60 igam = 1 , 3
                  skew = etijk(ialp,igam,idel)
                  if (skew.ne.0.0d0) then
                     do 50 i = 1 , nat
                        rot(n) = rot(n) + skew*r(igam,i)
     +                           *dipd((i-1)*3+ibet,ialp)
 50                  continue
                  end if
 60            continue
 70         continue
 80      continue
 90   continue
      write (iwr,6020)
      do 110 i =1,9
 110  write (iwr,6030) head(i),rot(i)
      return
 6010 format (/
     + 20x,'**********************************'/
     + 20x,'* dipole translational sum rules *'/
     + 20x,'**********************************'/)
 6020 format (/
     + 20x,'*******************************'/
     + 20x,'* dipole rotational sum rules *'/
     + 20x,'*******************************'/)
 6030 format (10x,a3,3x,f11.6)
      end
      subroutine prfreq (icode,zz)
c..  icode=1 --> frequency in a.u.
c..  icode=2 --> frequency squared in a.u.
      implicit real*8  (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character *4  i,blank,char
      data i,blank/'i ','  '/
      z = dabs(zz) + 1.0d-25
      if (icode.eq.2) z = dsqrt(z)
      ev = z*27.2107d0
      cm = z*219474.6d0
      ang = 1.0d8/cm
      char = blank
      if (icode.eq.2 .and. zz.lt.0.0d0) char = i
      write (iwr,6010) z , char , ev , char , cm , char , ang
      return
 6010 format (1x,f20.7,a4,'hartrees,',f20.8,a4,'ev,',f12.2,a4,
     +        'wavenumbers,',f12.2,' angstroms')
      end
      subroutine qhq2(h,f,c,t,iky,ncol,nrow,ndim)
      implicit real*8  (a-h,o-z)
c     two index transformation
      dimension h(*),f(*),c(*),t(*),iky(*)
      data zero/0.0d0/
      nj = 0
      ij = 0
      do 40 j = 1 , ncol
         do 20 k = 1 , nrow
            ckj = c(k+nj)
            kk = iky(k) + 1
            klen = k
            t(k) = ddot(klen,f(kk),1,c(nj+1),1)
            kless1 = k - 1
            if (kless1.gt.0 .and. ckj.ne.zero) then
                call daxpy(kless1,ckj,f(kk),1,t,1)
            endif
 20      continue
         ni = 0
         do 30 i = 1 , j
            ij = ij + 1
            h(ij) = ddot(nrow,t(1),1,c(ni+1),1)
            ni = ni + ndim
 30      continue
         nj = nj + ndim
 40   continue
      return
      end
      subroutine vwrt3 (v,n,nv,ibl,idev)
      implicit real*8  (a-h,o-z)
      dimension v(n,nv)
      data ilen/511/
      call search(ibl,idev)
      do 30 iv = 1 , nv
         i = 1
         k = n
 20      nw = min(k,ilen)
         call put(v(i,iv),nw,idev)
         i = i + nw
         k = k - nw
         if (k.gt.0) go to 20
 30   continue
      return
      end
      subroutine ver_anald(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/anald.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
