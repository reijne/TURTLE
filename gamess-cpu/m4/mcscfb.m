c 
c  $Author: jvl $
c  $Date: 2009-05-06 17:18:58 +0200 (Wed, 06 May 2009) $
c  $Locker:  $
c  $Revision: 5957 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mcscfb.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   mcorbopt     =
c ******************************************************
c ******************************************************
      subroutine sigma (q,vec,sigm,iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/jobopt)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/mcaddr)
      dimension q(*),vec(*),sigm(*)
c
c..   control routine to construct hessian on vector
c
c...  isigma = 1    full hessian
c              2    orbital hessian
c              3    ci hamiltonian
c
c...  iaugmx = 1    augmented matrix
c
      call accnt('sigma',1)
c
      goto (10,20,30),isigma
c
c
10    continue
      call full1 (q(1),vec,sigm,q(izint),q(igam),q(izint1),
     1           q(igam1),q(icivec),iwrite)
      n = nvar
      goto 40
c
20    continue
      call orb1 (q(1),vec,sigm,q(izint),q(igam),q(izint1),q(igam1),
     +           iwrite )
      n = nrot
      goto 40
c
30    call ci1 (q(1),vec,sigm,q(idiag),q(izint))
      n = nci
      goto 40
c
40    if (iaugmx.eq.1 .and. isigma.ne.1) sigm(n+1) = 0.0d0
      call accnt(' ',1)
      return
      end
      subroutine full0 (q,grad,diag,zint,gam,civec)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
INCLUDE(common/mcaddr)
      dimension q(*),grad(nvar),diag(nvar),zint(ne),gam(ne),civec(nci)
      call vclr(grad,1,nvar)
c
      call orb0 (q(1),grad,diag,zint,gam,civec)
c
      call fci0 (q(1),civec,zint,grad(nrot+1),diag(nrot+1))
      return
      end
      subroutine full1 (q,vec,sigma,zint,gam,zint1,gam1,civec,
     +                  iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
INCLUDE(common/mcaddr)
      dimension vec(nvar),sigma(nvar),zint(ne),gam(ne)
     1        ,zint1(ne),gam1(ne)
     2    ,civec(nci)
      dimension q(*),izpos(8)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c      call outvec (vec,nvar,'vector given to full1')
      signh =  dfloat(isignh)
      call vclr(sigma,1,nvar)
c
      if (ncoremc.gt.0) then
      call denst1 (q,civec,vec(nrot+1),gam1)
      zz = -ddot(nci,civec,1,vec(nrot+1),1)
      call daxpy(nact**2
     +  ,zz,gam(ic1d+1),1,gam1(ic1d+1),1)
      end if
c
      call orb1 (q(1),vec,sigma,zint,gam,zint1,gam1,iwrite)
c
      call fci1 (q(1),vec(nrot+1),sigma(nrot+1),zint,zint1,gam,gam1
     1 ,civec)
c
      call accnt('full1',2)
c...  contribution of orbital part of sigma from ci part of vec
c...  damp coupling
c      if (slast.gt.1.0e-1) call sscal (ne,0.98,gam1,1)
      ibase = icorr(0)
      izz = icorr(nbasis*nact)
      iqq = icorr(nact**2)
      call vclr(q(izz),1,iqq-izz)
      isingl = ne + 1
      do 20 i=1,nact
      do 20 j=1,i
      do 10 k=1,nact
      do 10 l=1,nact
      q(iqq+(k-1)*nact+l-1) = gam1(ind(ilifa(l)+k,ilifa(i)+j))
     1         + signh * gam1(ind(ilifa(k)+l,ilifa(i)+j))
      if (i.ne.j) q(iqq+(k-1)*nact+l-1) = q(iqq+(k-1)*nact+l-1)
     1                 + gam1(ind(ilifa(l)+k,ilifa(j)+i))
     2         + signh * gam1(ind(ilifa(k)+l,ilifa(j)+i))
10    continue
c
      call mxmb (zint(isingl),1,nbasis, q(iqq),1,nact, q(izz),1,nbasis,
     1  nbasis,nact,nact)
c
20    isingl = isingl + nbasis*nact
c
      do 30 k=1,nact
      do 30 l=1,nact
30    q(iqq+(k-1)*nact+l-1) = gam1(ic1d+ilifa(l)+k)
     1         + signh * gam1(ic1d+ilifa(k)+l)
      call mxmb (zint(isingl),1,nbasis, q(iqq),1,nact, q(izz),1,nbasis,
     1 nbasis,nact,nact)
c
c... slot in by symmetry for assembly routine
      call corlsr (iqq)
      do 50 isyma=1,nirrr
      loop=nsymm(isyma)*nprm(isyma)
      izpos(isyma) = icorr(loop)
      if(loop.gt.0)
     *call vclr(q(izpos(isyma)),1,loop)
      ioffs = -1+izpos(isyma)
      do 40 ib=ncor(isyma)+1,nprm(isyma)
      ibb = iorbsm(ib-1+istart(isyma)) - nst
      do 40 ia=1,nsymm(isyma)
      iaa = iorbsm(ia-1+istart(isyma))
c60    q(ioffs+(ib-1)*nsymm(isyma)+ia) = - q(izz-1+iaa+ibb*nbasis)
40    q(ioffs+(ib-1)*nsymm(isyma)+ia) =  q(izz-1+iaa+ibb*nbasis)
50    continue
      call orbasm (q(1),izpos,sigma,isignh)
      call corlsr (ibase)
c
c...  projectors
      zz = -ddot(nci,civec,1,vec(nrot+1),1)
      call daxpy(nci
     + ,zz,q(igrad+nrot),1,sigma(nrot+1),1)
      call orth (nci,sigma(nrot+1),civec)
c
c      call outvec (sigma,nvar,'vector returned by full1')
      call accnt(' ',2)
      return
      end
      function hamilt (gam,zint)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
      dimension gam(ne),zint(ne)
      val = ddot(ne,gam,1,zint,1)
      ijij=0
      do 1 ij=1,nact**2
      ijij = ijij + ij
1     val = val + (-0.5d0)*gam(ijij)*zint(ijij)
      hamilt = val
      return
      end
      subroutine hesini(q,iwrite)
      implicit REAL  (a-h,o-z)
      logical btest
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/mcaddr)
      dimension q(*)
c
      call accnt('hesini',1)
      n = nrot+nci
      if (isigma.eq.2) n = nrot
      m=n
      isigm1=isigma
c...  this bit of code to test core available & give realistic error msg
      junk = 4*ne + nci + 2*n + 1
      if (isigma.eq.1) junk = junk + nbasis*nact*((nact*(nact+1))/2+1)
      if (iwrnr.eq.1) then
      isigma=2
      m=nrot
      junk = 3*ne + 2*nstate*nci + 6*nrot + 6*maxbas**2 + 5
      if(igwgt.gt.0) junk = junk + nrot
      do 10 isym=1,nirrr
      junk = junk + 2*maxbas**2
      junk = junk + 3*nsymm(isym)**2
      junk = junk + 2*nsymm(isym)*nprm(isym)
10    junk = junk + 2*nact**2
      end if
      if(iter.eq.0.and.mcprin) write(iwrite,15) junk
15    format(/' estimated minimum core required:',i8/)
      junk = icorr(junk)
      call corlsr (junk)
c
c     shift = 0
      igrad = icorr(m+1)
      idiag = icorr(m+1)
      q(igrad+m) = 0.0d0
      q(idiag+m) = 0.0d0
c...  above nonsense for augmented matrix approach
      izint = icorr(ne)
      call vclr(q(izint),1,ne)
c..   for single excitation integrals
      if (isigma.eq.1) junk  = icorr(nbasis*nact*((nact*(nact+1))/2+1))
      igam  = icorr(ne)
      izint1 = icorr(ne)
      if(iwrnr.eq.0) igam1  = icorr(ne)
      icivec = icorr(nci*nstate)
      goto (20,30), isigma
c
20    call full0 (q(1),q(igrad),q(idiag),q(izint),q(igam),q(icivec))
      goto 40
c
30    call orb0 (q(1),q(igrad),q(idiag),q(izint),q(igam),q(icivec))
      goto 40
c
40    energy = core + hamilt (q(igam),q(izint))
      gradnt = dnrm2(m,q(igrad),1)
      if (btest(iprint,7)) call outvec (q(idiag),n,'diagonal elements of
     1 hessian')
      if (btest(iprint,8)) call outvec (q(igrad),m,'gradient')
c
      if (btest(iprint,7)) then
      isave = iaugmx
      iaugmx = 0
      ivec = icorr(m)
      ihes = icorr(m**2)
      do 50 i=1,m
      call vclr(q(ivec),1,m)
      q(ivec-1+i)=1.0d0
50    call sigma(q(1),q(ivec),q(ihes+(i-1)*m),iwrite)
      call outsqr (q(ihes),m,m,m,'hessian matrix')
      call corlsr (ivec)
      iaugmx = isave
      end if
      if(iwrnr.eq.1) call corlsr(icivec)
      isigma=isigm1
      call accnt('other',1)
      return
      end
      subroutine orb0(q,grad,diag,zint,gam,civec)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common /block / nosymm,nosym(8),iosym(8,maxorb)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
INCLUDE(common/mcaddr)
      common /intbuf/ intpos,intfil,intmod
INCLUDE(common/mcff)
      common/lsort/v(maxorb),ijpos(8),ikpos(8),itpos(8),igpos(8)
     1               ,ifpos(8),iprev(8),jprev(8)
      logical cigr,cicv,mpgr,mcgr,cicx,umpgr
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,
     +  cicv,mnnr,
     +  mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,isecnd,isecsy,
     +  irlagr,iadfrc,nfc,intlgr,umpgr,ispaer(26),nd2mo,
     +  nncore,ncact,nvr,ifilh,iblkh,iblk1,
     +  iblk22,ntot,nupact,ijr3,rspare(70)
      dimension q(*),grad(nrot),diag(nrot),zint(ne),civec(*),gam(ne)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      call accnt('orb0',2)
      ibase = icorr(0)
      call cget(civec,nstate)
      call densav (q(1),civec,gam)
c
      do 10 isym=1,nirrr
      ijpos(isym) = icorr(maxbas**2)
      ikpos(isym) = icorr(maxbas**2)
      ifpos(isym) = icorr(nsymm(isym)*nprm(isym))
      itpos(isym) = icorr(nsymm(isym)*nprm(isym))
      igpos(isym) = icorr(nsymm(isym)**2)
10    continue
      ilpos = icorr(maxbas**2)
      iltpos = icorr(maxbas**2)
      ippos = icorr(nact**2)
      iqpos = icorr(nact**2)
      call vclr(q(ijpos(1)),1,iqpos+nact**2-ijpos(1))
c
      signh =  dfloat(isignh)
c
c...  loop over integrals
      intpos=0
      intfil=num6
ctest iblf  =iblk6
      intmod=0
      call search (iblk6,num6)
      do 110 i=1,nprim
      isymi = itype(i)
      do 20 is=1,nirrr
20    jprev(is) = 0
      do 110 j=1,i
      isymj = itype(j)
      isymij = mults(isymi,isymj)
      nsymij = nsymm(isymi) * nsymm(isymj)
      call loadjk (q(1),ijpos,ikpos,i,j,zint)
      if (i.le.ncoremc .and. i.ne.j) goto 100
      if (j.le.ncoremc) then
c...  l matrix required
      call trnsps(q(ikpos(isymj)),q(ilpos),nsymm(isymj),nsymm(isymi))
      if (isignh.eq.1) then
       call dscal(nsymij,(-1.0d0),q(ilpos),1)
          endif
      if(isignh.eq.1) then
      call daxpy(nsymij
     +  ,4.0d0,q(ikpos(isymi)),1,q(ilpos),1)
          endif
      call vsub(q(ilpos),1,q(ijpos(isymi)),1,q(ilpos),1,nsymij)
      if(i.ne.j)
     *call trnsps(q(ilpos),q(iltpos),nsymm(isymi),nsymm(isymj))
      if (i.le.ncoremc) then
c...  core/core ii
c...  t(p,i) = 2 l(p,p)
      nsi=nsymm(isymi)         
      ix1=itpos(isymi)+jprev(isymi)
      call daxpy(nsi,2.0d0,q(ilpos),nsi+1,q(ix1),1)
c
      else
c...  active/core t=i
c...  t(v,j) = -signh *2 *u(v,j) = -signh *2 *gam(v,t)*l(j,t)
      if (isymi.eq.isymj) then
      ill = iltpos+jprev(isymi)+ncor(isymi)
      ivv=itpos(isymi)+jprev(isymi)+ncor(isymi)
      do 30 iv=1,nact
      if (itypea(iv).ne.isymi) goto 30
      q(ivv) = q(ivv) + (-2.0d0)*signh*gam(ic1d+(i-nst)*nact+iv)
     1    * q(ill)
      ill = ill + 1
      ivv = ivv+1
30    continue
      end if
      end if
c
      else
c...  active i,j
      onedel=1.0d0
      if (i.eq.j) onedel=0.0d0
      twodel = 1.0d0+onedel
c
      if (isymij.eq.1) then
c...   contrib to core+active fock matrix (g)
      fj = twodel*gam(ic1d+(i-nst)*nact+j-ncoremc)
      fk = fj*(-0.5d0)
      do 40 isymp=1,nirrr
      nsym2=nsymm(isymp)*nsymm(isymp)
      if(nsym2.eq.0)go to 40
      call daxpy(nsym2
     +  ,fj,q(ijpos(isymp)),1,q(igpos(isymp)),1)
      call daxpy(nsym2
     +  ,fk,q(ikpos(isymp)),1,q(igpos(isymp)),1)
40    continue
      end if
c
c...  p(u,t) = gam(u,t,i,j) + onedel*gam(u,t,j,i)
c...  q(u,t) =(gam(t,i,j,u) + sign  *gam(t,i,u,j) ) * twodel
         ivx=ilifa(i-ncoremc)+j-ncoremc
         do 90 isymt=1,nirrr
         isymu=mults(isymt,isymij)
         itu=-1
         do 60 it=1,nact
         if (isymt.ne.itypea(it)) goto 60
         do 50 iu=1,nact
         if (isymu.ne.itypea(iu)) goto 50
         itu=itu+1
         q(ippos+itu) = gam(ind(ivx,ilifa(it)+iu))
     1                 +onedel * gam(ind(ivx,ilifa(iu)+it))
         q(iqpos+itu) = (gam(ind(ilifa(it)+i-ncoremc,
     +                   ilifa(j-ncoremc)+iu)) + signh * 
     +                   gam(ind(ilifa(it)+i-ncoremc,
     +                   ilifa(iu)+j-ncoremc)) )*twodel
50       continue
60       continue
c
c
c...  f(p,t) = j(p,u)*gam(u,t)
         call mxmb (q(ijpos(isymt)+ncor(isymu)*nsymm(isymt))
     1 ,1,nsymm(isymt), q(ippos),1,nactt(isymu),
     2 q(ifpos(isymt)+ncor(isymt)*nsymm(isymt)),1,nsymm(isymt),
     3 nsymm(isymt),nactt(isymu),nactt(isymt))
c
      if (isymt.eq.isymu) then
c...  t(t,u) = -sign*(p(t,u)*j(t,u) + q(t,u)*k(t,u))
      itad = itpos(isymt) + ncor(isymt)*(nsymm(isymt)+1) - 1
      ipad = ippos-1
      iqad = iqpos-1
      ijad = itad - itpos(isymt) + ijpos(isymt)
      ikad = itad - itpos(isymt) + ikpos(isymt)
      do 80 it=1,nactt(isymt)
      do 70 iu=1,nactt(isymt)
70    q(itad+iu) = q(itad+iu) - signh * (q(ipad+iu)*q(ijad+iu)
     1                                 + q(iqad+iu)*q(ikad+iu) )
      ipad = ipad + nactt(isymt)
      iqad = iqad + nactt(isymt)
      ijad = ijad + nsymm(isymt)
      ikad = ikad + nsymm(isymt)
80    itad = itad + nsymm(isymt)
c
c...  t(p,t) = j(p,p)*p(t,t) + k(p,p)*q(t,t)
      call mxmb (q(ijpos(isymt)),nsymm(isymt)+1,0,
     1 q(ippos),0,nactt(isymt)+1,
     2 q(itpos(isymt)+ncor(isymt)*nsymm(isymt)),1,nsymm(isymt),
     3 nsymm(isymt),1,nactt(isymt) )
      call mxmb (q(ikpos(isymt)),nsymm(isymt)+1,0,
     1 q(iqpos),0,nactt(isymt)+1,
     2 q(itpos(isymt)+ncor(isymt)*nsymm(isymt)),1,nsymm(isymt),
     3 nsymm(isymt),1,nactt(isymt) )
      end if
90    continue
      end if
100   jprev(isymj) = jprev(isymj)+nsymm(isymj)
110   continue
c
c...  process core hamiltonian
      call loadjk (q(1),ijpos,ikpos,0,0,zint)
c
c...  contribution to g matrix
      do 120 isymp=1,nirrr
      nsym2=nsymm(isymp)*nsymm(isymp)
      ig=igpos(isymp)
      if(nsym2.gt.0)
     *call vadd(q(ig),1,q(ijpos(isymp)),1,q(ig),1,nsym2)
120   continue
c
c...  form 1 pdm by symmetry
      do 170 isymt=1,nirrr
      itu = ippos-1
      do 140 it=1,nact
      if (itypea(it).ne.isymt) goto 140
      do 130 iu=1,nact
      if (itypea(iu).ne.isymt) goto 130
      itu = itu + 1
      q(itu) = gam(ic1d+ilifa(it)+iu)
130   continue
140   continue
c
c...  f(p,t) = j(p,u) * gam(u,t)
      call mxmb(q(ijpos(isymt)+ncor(isymt)*nsymm(isymt)),1,nsymm(isymt),
     1 q(ippos),1,nactt(isymt),
     2 q(ifpos(isymt)+ncor(isymt)*nsymm(isymt)),1,nsymm(isymt),
     3 nsymm(isymt),nactt(isymt),nactt(isymt) )
c
c...  t(p,t) = j(p,p)*gam(t,t)
      call mxmb (q(ijpos(isymt)),nsymm(isymt)+1,0,
     1 q(ippos),0,nactt(isymt)+1,
     2 q(itpos(isymt)+ncor(isymt)*nsymm(isymt)),1,nsymm(isymt),
     3 nsymm(isymt),1,nactt(isymt) )
c
c...  t(t,u) = -sign * j(t,u)*gam(t,u)
      itad = itpos(isymt) + ncor(isymt)*(nsymm(isymt)+1) - 1
      ipad = ippos-1
      ijad = itad - itpos(isymt) + ijpos(isymt)
      do 160 it=1,nactt(isymt)
      do 150 iu=1,nactt(isymt)
150   q(itad+iu) = q(itad+iu) - signh * q(ipad+iu)*q(ijad+iu)
      ipad = ipad + nactt(isymt)
      ijad = ijad + nsymm(isymt)
160   itad = itad + nsymm(isymt)
170   continue
c
c
c...  processing of g matrix
      do 190 isymp=1,nirrr
      if(nsymm(isymp).le.0)go to 1000
c...  double up/symmetrise g(p,q) = g(p,q) + g(q,p)
      igp = igpos(isymp)
      ijp = ijpos(isymp)
      call trnsps(q(igp),q(ijp),nsymm(isymp),nsymm(isymp))
      call vadd(q(ijp),1,q(igp),1,q(igp),1,nsymm(isymp)**2)
c...  f(p,i) = g(p,i)
1000  if (ncor(isymp).gt.0) then
      ncns=ncor(isymp)*nsymm(isymp)       
      call dcopy(ncns,q(igpos(isymp)),1,q(ifpos(isymp)),1)
c
c...  t(p,i) = t(p,i)  + g(p,p)
      itt = itpos(isymp)
      do 180 i=1,ncor(isymp)
      igg = igpos(isymp)
      do 180 ip=1,nsymm(isymp)
      q(itt) = q(itt) + q(igg)
      igg = igg + nsymm(isymp)+1
180   itt = itt + 1
      end if
190   continue
c
c...  skip h0
      junk=icorr(ltrimo)
      call intin(q(junk),ltrimo)
      call corlsr(junk)
c...  write out g matrix onto end of mo integral interface
      call search (iposun(num6)-1,num6)
      intmod = -1
      do 200 isymp=1,nirrr
      igg = igpos(isymp)
      do 200 ip=1,nsymm(isymp)
      call intou (q(igg),ip)
200   igg = igg + nsymm(isymp)
c...  write out f matrix
      irlagr = iposun(num6)
      intlgr = intpos
      do 210 isymp=1,nirrr
210   call intou (q(ifpos(isymp)),nsymm(isymp)*nprm(isymp))
      call intend
c
c
c...  t(p,ti) = t(p,ti) - f(ti,ti)
c...  care of double counting?
      do 240 isymp=1,nirrr
      itt = itpos(isymp) - 1
      iff = ifpos(isymp)
      do 230 it=1,nprm(isymp)
      do 220 ip=1,nsymm(isymp)
220   q(itt+ip) = q(itt+ip) - q(iff)
      itt = itt + nsymm(isymp)
230   iff = iff + nsymm(isymp) + 1
240   continue
c
c...  load diagonal elements
      call vclr(diag,1,nrot)
      call orbasm (q(1),itpos,diag,-1)
_IF(cray,ibm,vax)
      do 25 i=1,nrot
   25 diag(i) = -diag(i)
_ELSE
      call vneg(diag,1,diag,1,nrot)
_ENDIF
c
c...  construct gradient
      call vclr(grad,1,nrot)
      if (isignh.eq.1) then
      call orbasm (q(1),ifpos,grad,1)
      gradnt = dnrm2(nrot,grad,1)
      end if
c
      call corlsr (ibase)
      call accnt(' ',2)
      return
      end
      subroutine orb1 (q,vec,sigma,zint,gam,zint1,gam1,iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
      common /intbuf/ intpos,intfil,intmod
INCLUDE(common/mcff)
      common/lsort/v(maxorb),ijpos(8),ikpos(8),izpos(8),ixpos(8),w(31)
     1 ,iprev(8),jprev(8),iupos(8),impos(8),inpos(8)
      dimension q(*),vec(nrot),sigma(nrot),zint(ne),gam(ne),zint1(ne)
     1              ,gam1(ne)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      ibase = icorr(0)
      call accnt('orb1',2)
      enext = potnuc+efreez
      do 10 isym=1,nirrr
      ijpos(isym) = icorr(maxbas**2)
      ikpos(isym) = icorr(maxbas**2)
      izpos(isym) = icorr(nsymm(isym)*nprm(isym))
      ixpos(isym) = icorr(nsymm(isym)**2)
      iupos(isym) = ixpos(isym)
      if (iwrnr.eq.1) iupos(isym) = icorr(nsymm(isym)**2)
      if (isigma.eq.1) impos(isym) = icorr(nact**2)
      if (isigma.eq.1) inpos(isym) = icorr(nact**2)
10    continue
      ilpos = icorr(maxbas**2)
      iltpos = icorr(maxbas**2)
      iwork1 = icorr(maxbas**2)
      iwork2 = icorr(maxbas**2)
      ippos = icorr(nact**2)
      iqpos = icorr(nact**2)
      call vclr(q(ibase),1,iqpos+nact**2-ibase)
c
c...  load the x or t=exp(x)-1 matrix
      call xasm(q(1),ixpos,vec)
      if (iwrnr.eq.1) then
c...  load u for non-linear
      do 20 isym=1,nirrr
      loop=nsymm(isym)
      if(loop.le.0)go to 20
      if (lto(10)) call outsqr(q(ixpos(isym)),loop,loop
     1     ,loop,'t matrix')
      ix1=iupos(isym)
      call vclr(q(ix1),1,loop**2)
      call vfill(1.0d0,q(ix1),loop+1,loop)
      call vadd(q(ix1),1,q(ixpos(isym)),1,q(ix1),1,loop**2)
 20   continue
      end if
c
c...  +1 for real hessian, -1 for imaginary
      signh =  dfloat(isignh)
c
c...  damp orb-ci coupling
      damp = 1.0d0
c      if (slast.gt.1.e-1 .and. iwrnr.ne.1 .and. isignh.eq.1) damp =0.98
c
c... loop over integrals
      call vclr(zint1,1,ne)
      intpos = 0
      intfil = num6
ctest iblf   = iblk6
      do 30 is=1,nirrr
30    iprev(is)=0
      do 210 i=1,nprim
      isymi = itype(i)
      niprev = iprev(isymi)*nsymm(isymi)
      do 40 is=1,nirrr
40    jprev(is)=0
      do 200 j=1,i
      if (lto(10)) write(iwrite,1000)i,j
 1000  format(' integrals for ',2i5)
      isymj = itype(j)
      njprev = jprev(isymj)*nsymm(isymj)
      isymij = mults(isymi,isymj)
      nsymij = nsymm(isymi)*nsymm(isymj)
      call loadjk (q(1),ijpos,ikpos,i,j,zint)
c
      if (j.gt.ncoremc) goto 80
c...  j is core, so form l matrix
      call 
     *trnsps(q(ikpos(isymj)),q(ilpos),nsymm(isymj),nsymm(isymi))
         if (isignh.eq.1) then
      call dscal(nsymij,-1.0d0,q(ilpos),1)
          endif
      if(isignh.eq.1) then
      call daxpy(nsymij
     +  ,4.0d0,q(ikpos(isymi)),1,q(ilpos),1)
          endif
      call vsub(q(ilpos),1,q(ijpos(isymi)),1,q(ilpos),1,nsymij)
      if(i.ne.j)
     *call trnsps(q(ilpos),q(iltpos),nsymm(isymi),nsymm(isymj) )
c
      if (lto(10)) call outsqr (q(ilpos),nsymm(isymi),nsymm(isymi)
     1,nsymm(isymj),'l matrix')
      if (i.gt.ncoremc) goto 60
c...  ij = core-core: z(p,i) = 2*l(p,q)*x(q,j)
c... if i.ne.j       z(p,j) = 2*l(q,p)*x(p,i)
      call fmove (q(ixpos(isymj)+njprev),q(iwork1),nsymm(isymj))
      call dscal(nsymm(isymj),2.0d0,q(iwork1),1)
      enext = enext + q(izpos(isymi)+niprev+iprev(isymi))
c      print*,enext
      call mxmb (q(ilpos),1,nsymm(isymi), q(iwork1),1,0,
     1          q(izpos(isymi)+niprev),1,0, nsymm(isymi),nsymm(isymj),1)
      enext = enext - q(izpos(isymi)+niprev+iprev(isymi))
c      print*,enext
      fac = -0.5d0
      if (i.ne.j) fac = -1.0d0
      enext = enext +fac*q(ilpos+iprev(isymi)+jprev(isymj)*nsymm(isymi))
c      print*,enext
      if (i.eq.j) goto 50
      call fmove (q(ixpos(isymi)+niprev),q(iwork1),nsymm(isymi))
      call dscal(nsymm(isymi),2.0d0,q(iwork1),1)
      enext = enext + q(izpos(isymj)+njprev+jprev(isymj))
c      print*,enext
      call mxmb (q(iltpos),1,nsymm(isymj), q(iwork1),1,0,
     1         q(izpos(isymj)+njprev),1,0, nsymm(isymj),nsymm(isymi),1)
      enext = enext - q(izpos(isymj)+njprev+jprev(isymj))
c      print*,enext
50    continue
      goto 200
c
c...  i.gt.ncoremc.ge.j
c...  z(p,t) =l(p,q) * x(q,j)? * gamma(t,i)?
60    call mxmaa(q(ilpos),1,nsymm(isymi),
     1 q(ixpos(isymj)+jprev(isymj)*nsymm(isymj)),1,nsymm(isymj),
     2 q(iwork1),1,0, nsymm(isymi),nsymm(isymj),1)
c
cpaul test occurrence of virtuals of the right symmetry
      if (nsecc(isymi).eq.0) goto 200
      ii=1
      do 70 k=1,nact
      if (itypea(k).ne.isymi) goto 70
c...  for ci part of full hessian
      val = 0.5d0*q(iwork1+ii+ncor(isymi)-1)
      ix1 = ixpos(isymi)+(ii+ncor(isymi)-1)*nsymm(isymi)
      if (iwrnr.eq.1) val = -0.25d0 * q(ilpos +
     1       jprev(isymj)*nsymm(isymi)+ii-1+ncor(isymi) ) +
     2  ddot(nsymm(isymi),q(iwork1),1,q(ix1),1)
      zint1 (ic1d+ilifa(i-ncoremc)+k) = 
     +        zint1(ic1d+ilifa(i-ncoremc)+k)+val
      zint1 (ic1d+ilifa(k)+i-ncoremc) = 
     +        zint1(ic1d+ilifa(k)+i-ncoremc)+val
      w(ii) = gam1(ic1d+ilifa(i-ncoremc) + k) * damp
      v(ii) = gam(ic1d+ilifa(k)+i-ncoremc)
      ii=ii+1
70    continue
c
      enext=enext+dsum(nprm(isymi),q(izpos(isymi)),nsymm(isymi)+1)
c      print*,enext
      call mxmb (q(iwork1),1,0, v,0,1,
     1 q(izpos(isymi)+ncor(isymi)*nsymm(isymi)),1,nsymm(isymi),
     2 nsymm(isymi),1,nactt(isymi))
      ix1 = ilpos+jprev(isymj)*nsymm(isymi)+ncor(isymi)
      enext=enext-dsum(nprm(isymi),q(izpos(isymi)),nsymm(isymi)+1)
     1 - 0.5d0*ddot(nactt(isymi),q(ix1),1,v,1)
c      print*,enext
c
c..  contrib to orb-ci part
      if (isigma.eq.1 .and. iwrnr.ne.1) then
c..  z(p,j) = l(v,p) * gam1 (i,v)
      call mxmb (q(iltpos+ncor(isymi)*nsymm(isymj)),1,nsymm(isymj),
     1 w,1,0, q(izpos(isymj)+njprev),1,0, nsymm(isymj),nactt(isymj),1)
      end if
c...  z(p,j) = l(q,p)x(q,t)*gam(t,i)?
      call mxmaa(q(ixpos(isymi)+ncor(isymi)*nsymm(isymi)),1,nsymm(isymi)
     1,v,1,0, q(iwork1),1,0, nsymm(isymi),nactt(isymi),1)
      enext = enext + q(izpos(isymj)+njprev+jprev(isymj))
c      print*,enext
      call mxmb (q(iltpos),1,nsymm(isymj), q(iwork1),1,0,
     1 q(izpos(isymj)+jprev(isymj)*nsymm(isymj)),1,0,
     2 nsymm(isymj),nsymm(isymi),1)
      enext = enext - q(izpos(isymj)+njprev+jprev(isymj))
c      print*,enext
      goto 200
c
c...  i,j active
80    ivx = ilifa(i-ncoremc)+j-ncoremc
      ixv = ilifa(j-ncoremc) + i-ncoremc
      onedel=1.0d0
      if (i.eq.j) onedel=0.0d0
      if (isigma.eq.1) then
      do 90 isymt=1,nirrr
      if (.not.lto(10)) goto 91
      call outsqr(q(ijpos(isymt)),nsymm(isymt),nsymm(isymt),nsymm(isymt)
     1,'coulomb')
      call outsqr(q(ikpos(isymt)),nsymm(isymt),nsymm(isymt),nsymm(isymt)
     1,'exchange')
 91   continue
      call vclr(q(impos(isymt)),1,nact**2)
90    continue
      end if
      do 140 isymt=1,nirrr
      isymu = mults(isymt,isymij)
c
c...  set up density matrices for this symmetry
c...  p(t,u) = gam(u,t,i,j) + onedel*gam(u,t,j,i)
c...  q(t,u) = gam(u,i,j,t) + sign*gam(u,i,t,j)
      itu=-1
      do 110 iu=1,nact
      if (isymu.ne.itypea(iu)) goto 110
      do 100 it=1,nact
      if (isymt.ne.itypea(it)) goto 100
      itu = itu + 1
      q(ippos+itu) = gam(ind(ivx,ilifa(it)+iu))
     1    + onedel * gam(ind(ivx,ilifa(iu)+it))
      q(iqpos+itu) = gam(ind(ilifa(iu)+i-ncoremc,
     +               ilifa(j-ncoremc)+it))
     1    + signh  * gam(ind(ilifa(iu)+i-ncoremc,
     +               ilifa(it)+j-ncoremc))
100   continue
110   continue
c
c....   c o u l o m b   i n t e g r a l s
c...  z(p,u) =j(p,q)*u(q,t)? * p(t,u)
      call mxmaa(q(ijpos(isymu)),1,nsymm(isymu),
     1 q(iupos(isymt)+ncor(isymt)*nsymm(isymt)),1,nsymm(isymt),
     2 q(iwork1),1,nsymm(isymu), nsymm(isymu),nsymm(isymt),nactt(isymt))
      call mxmb (q(iwork1),1,nsymm(isymu),q(ippos),1,nactt(isymt),
     1 q(izpos(isymu)+ncor(isymu)*nsymm(isymu)),1,nsymm(isymu),
     2 nsymm(isymu),nactt(isymt),nactt(isymu) )
c
      if (isigma .eq. 1 ) then
c...  load these transformed integrals for the ci part
      if (iwrnr.eq.1) then
      call mxmb (q(iupos(isymu)+ncor(isymu)*nsymm(isymu)),nsymm(isymu),1
     1,q(iwork1),1,nsymm(isymu), q(impos(isymu)),1,nactt(isymu),
     2 nactt(isymu),nsymm(isymu),nactt(isymt) )
c      if (lto(10)) print*,'m for symmetry ',isymu
c      if (lto(10)) calloutsqr(q(impos(isymu)),nactt(isymu),nactt(isymu)
c     <,nactt(isymt),'m')
      else
      itu = iwork1-1+ncor(isymu)
      do 130 it=1,nactt(isymt)
      do 120 iu=1,nactt(isymu)
      itu = itu + 1
      q(impos(isymu)+(iu-1)+(it-1)*nactt(isymu)) =
     1q(impos(isymu)+(iu-1)+(it-1)*nactt(isymu)) + q(itu)
120   q(impos(isymt)+(it-1)+(iu-1)*nactt(isymt)) =
     1q(impos(isymt)+(it-1)+(iu-1)*nactt(isymt)) + signh*q(itu)
130   itu = itu + nsecc(isymu)+ncor(isymu)
      end if
      end if
c
c....  e x c h a n g e   i n t e g r a l s
c...  z(p,u) =k(p,q)*x(q,t)? * q(t,u)
      enext = enext + dsum(nprm(isymu),q(izpos(isymu)),nsymm(isymu)+1)
c      print*,enext
      call mxmaa(q(ikpos(isymu)),1,nsymm(isymu),
     1 q(ixpos(isymt)+ncor(isymt)*nsymm(isymt)),1,nsymm(isymt),
     2 q(iwork1),1,nsymm(isymu), nsymm(isymu),nsymm(isymt),nactt(isymt))
      call mxmb (q(iwork1),1,nsymm(isymu),q(iqpos),1,nactt(isymt),
     1 q(izpos(isymu)+ncor(isymu)*nsymm(isymu)),1,nsymm(isymu),
     2 nsymm(isymu),nactt(isymt),nactt(isymu) )
      enext = enext - dsum(nprm(isymu),q(izpos(isymu)),nsymm(isymu)+1)
c      print*,enext
c
c...  store t(dagger).k.t = n in the case of non-linear with coupling
      if (iwrnr.eq.1 .and. isigma.eq.1)  then
      call mxmaa(q(ixpos(isymu)+ncor(isymu)*nsymm(isymu)),nsymm(isymu),1
     1,q(iwork1),1,nsymm(isymu), q(inpos(isymu)),1,nactt(isymu),
     2 nactt(isymu),nsymm(isymu),nactt(isymt) )
      if (i.eq.j) then
      ix1 = inpos(isymu)
      nnn = nactt(isymu)*nactt(isymt)
       call dscal(nnn,0.5d0,q(ix1),1)
         endif
c      if (lto(10)) print*,'n for symmetry ',isymu
c      if (lto(10)) calloutsqr(q(inpos(isymu)),nactt(isymu),nactt(isymu)
c     <,nactt(isymt),'n')
      end if
c
      if (i.ne.j) then
c...  z(p,t) =k(q,p)*x(q,u)? * q(t,u)
      enext = enext + dsum(nprm(isymt),q(izpos(isymt)),nsymm(isymt)+1)
c      print*,enext
      call mxmaa(q(ikpos(isymu)),nsymm(isymu),1,
     1 q(ixpos(isymu)+ncor(isymu)*nsymm(isymu)),1,nsymm(isymu),
     2 q(iwork1),1,nsymm(isymt), nsymm(isymt),nsymm(isymu),nactt(isymu))
      call mxmb (q(iwork1),1,nsymm(isymt),q(iqpos),nactt(isymt),1,
     1 q(izpos(isymt)+ncor(isymt)*nsymm(isymt)),1,nsymm(isymt),
     2 nsymm(isymt),nactt(isymu),nactt(isymt) )
      enext = enext - dsum(nprm(isymt),q(izpos(isymt)),nsymm(isymt)+1)
c      print*,enext
      end if
140   continue
c..   put the transformed integrals into zint1 as appropriate
      do 190 isymt=1,nirrr
      isymu = mults(isymt,isymij)
      if (isigma.eq.1) then
      itu = impos(isymt) - 1
      do 160 iu=1,nact
      if (itypea(iu).ne.isymu) goto 160
      do 150 it=1,nact
      if (itypea(it).ne.isymt) goto 150
      itu = itu + 1
      int1 = ind (ivx,ilifa(it)+iu)
      int2 = ind (ixv,ilifa(it)+iu)
      zint1(int1) = zint1(int1) + q(itu)
      if (int1.ne.int2) zint1(int2) = zint1(int2) + q(itu)
150   continue
160   continue
      end if
c
      if (iwrnr.eq.1) then
      if (isigma.eq.1 ) then
      itu = inpos(isymt) - 1
      do 180 iu=1,nact
      if (itypea(iu).ne.isymu) goto 180
      do 170 it=1,nact
      if (itypea(it).ne.isymt) goto 170
      itu = itu + 1
      int1 = ind (ilifa(it)+i-ncoremc,ilifa(iu)+j-ncoremc)
      zint1(int1) = zint1(int1) + q(itu)
      int1 = ind (ilifa(i-ncoremc)+it,ilifa(iu)+j-ncoremc)
      zint1(int1) = zint1(int1) + q(itu)
      int1 = ind (ilifa(it)+i-ncoremc,ilifa(j-ncoremc)+iu)
      zint1(int1) = zint1(int1) + q(itu)
      int1 = ind (ilifa(i-ncoremc)+it,ilifa(j-ncoremc)+iu)
      zint1(int1) = zint1(int1) + q(itu)
170   continue
180   continue
      end if
c..   contributions to core
      if (isymt.eq.isymu .and. ncor(isymt).gt.0) then
      nsym2=nsymm(isymt)**2
      ikp=ikpos(isymt)
      ijp=ijpos(isymt)
      call trnsps(q(ikp),q(iwork2),nsymm(isymt),nsymm(isymt))
      call vadd(q(ikp),1,q(iwork2),1,q(ikp),1,nsym2)
      call daxpy(nsym2,-4.0d0,q(ijp),1,q(ikp),1)
      call mxmaa(q(ikpos(isymt)),1,nsymm(isymt),
     1 q(iupos(isymt)),1,nsymm(isymt), q(iwork2),1,nsymm(isymt),
     2 nsymm(isymt),nsymm(isymt),ncor(isymt) )
      val = - gam(ic1d+ilifa(i-ncoremc)+j-ncoremc)
      if(i.eq.j) val=val*0.5d0
      nsnc = nsymm(isymt)*ncor(isymt)
      call daxpy(nsnc
     +  ,val,q(iwork2),1,q(izpos(isymt)),1)
      if (isigma.eq.1) then
       val = -ddot(nsnc,q(iwork2),1,q(iupos(isymt)),1)*0.5d0
      zint1(ic1d+ilifa(i-ncoremc)+j-ncoremc) =
     1    zint1(ic1d+ilifa(i-ncoremc)+j-ncoremc) + val
      if (i.ne.j) zint1(ic1d+ilifa(j-ncoremc)+i-ncoremc) =
     1    zint1(ic1d+ilifa(j-ncoremc)+i-ncoremc) + val
      end if
      end if
      end if
190   continue
c
      goto 200
200   jprev(isymj) = jprev(isymj) + 1
210   iprev(isymi) = iprev(isymi) + 1
c
c...  process core hamiltonian
      call loadjk (q(1),ijpos,ikpos,0,0,zint)
      do 260 isymt=1,nirrr
      if (iwrnr.eq.1 .and. ncor(isymt).gt.0) then
      call mxmaa(q(ijpos(isymt)),1,nsymm(isymt),
     1 q(iupos(isymt)),1,nsymm(isymt), q(iwork1),1,nsymm(isymt),
     2 nsymm(isymt),nsymm(isymt),ncor(isymt) )
      nsnc = nsymm(isymt)*ncor(isymt)
      call daxpy(nsnc
     +  ,2.0d0,q(iwork1),1,q(izpos(isymt)),1)
      end if
      itu=ippos-1
      do 230 it=1,nact
      if (itypea(it).ne.isymt) goto 230
      do 220 iu=1,nact
      if (itypea(iu).ne.isymt) goto 220
      itu=itu+1
      q(itu) = gam(ic1d + ilifa(it) + iu)
220   continue
230   continue
c
c...  z(p,u) =j(p,q)*u(q,t)? * gam(t,u)
      call mxmaa(q(ijpos(isymt)),1,nsymm(isymt),
     1 q(iupos(isymt)+ncor(isymt)*nsymm(isymt)),1,nsymm(isymt),
     2 q(iwork1),1,nsymm(isymt), nsymm(isymt),nsymm(isymt),nactt(isymt))
      call mxmb (q(iwork1),1,nsymm(isymt),q(ippos),1,nactt(isymt),
     1 q(izpos(isymt)+ncor(isymt)*nsymm(isymt)),1,nsymm(isymt),
     2 nsymm(isymt),nactt(isymt),nactt(isymt) )
c
      if (isigma.eq.1) then
c...  load transformed integrals for ci part
      itu = iwork1-1+ncor(isymt)
      itt=iwork1-nsymm(isymt)
      do 250 it=1,nact
      if (itypea(it).ne.isymt) goto 250
      itt = itt + nsymm(isymt)
      iuu = iupos(isymt) + nsymm(isymt)*(ncor(isymt)-1)
      do 240 iu=1,nact
      if (itypea(iu).ne.isymt) goto 240
      iuu = iuu + nsymm(isymt)
      itu = itu + 1
      val = q(itu)
      if (iwrnr.eq.1) val=ddot(nsymm(isymt),q(itt),1,q(iuu),1)
      intad = ic1d + ilifa(it)+iu
      zint1(intad) = zint1(intad) + val
      if (iwrnr.eq.1) goto 240
      intad = ic1d + ilifa(iu)+it
      zint1(intad) = zint1(intad) + signh*val
240   continue
      itu = itu + nsecc(isymt)+ncor(isymt)
250   continue
      end if
c
260   continue
c
      if (iwrnr.ne.1) then
c...  process g (core + active fock) matrix
      call loadjk (q(1),ijpos,ikpos,-1,-1,zint)
      do 270 isymt=1,nirrr
270   call mxmb (q(ijpos(isymt)),1,nsymm(isymt),
     1 q(iupos(isymt)),1,nsymm(isymt),
     2 q(izpos(isymt)),1,nsymm(isymt),
     3 nsymm(isymt),nsymm(isymt),ncor(isymt))
c
c...  process f (fock) matrix
      do 290 isymt=1,nirrr
      l = nsymm(isymt)*nprm(isymt)
      call intin (q(iwork1),l)
c
c...  w(p,q) = -0.5 (f(p,q)+f(q,p))
      call vclr(q(iwork2),1,nsymm(isymt)**2)
      call fmove (q(iwork1),q(iwork2),l)
      call dscal(l,(-0.5d0),q(iwork2),1)
      do 280 it=1,nprm(isymt)
      do 280 ip=1,nsymm(isymt)
280   q(iwork2-1+it+(ip-1)*nsymm(isymt)) =
     1q(iwork2-1+it+(ip-1)*nsymm(isymt)) + (-0.5d0) *
     2q(iwork1-1+ip+(it-1)*nsymm(isymt))
c
c...  z(p,t) = w(p,q) * x(q,t)
      call mxmb (q(iwork2),1,nsymm(isymt),
     1 q(ixpos(isymt)),1,nsymm(isymt),
     2 q(izpos(isymt)),1,nsymm(isymt),
     3 nsymm(isymt),nsymm(isymt),nprm(isymt))
c...  z(a,t) = w(a,u) * w(u,t)
      call mxmb (q(ixpos(isymt)+nprm(isymt)),1,nsymm(isymt),
     1 q(iwork2),1,nsymm(isymt),
     2 q(izpos(isymt)+nprm(isymt)),1,nsymm(isymt),
     3 nsecc(isymt),nprm(isymt),nprm(isymt) )
c
290   continue
      end if
c
      if (iwrnr.eq.1) then
c      print*,'energy before ub term ',enext
c...  form u(dagger)*b
      do 300 isym=1,nirrr
c     call outsqr(q(izpos(isym)),nsymm(isym),nsymm(isym),nprm(isym),'z')
c     calloutsqr(q(iupos(isym)),nsymm(isym),nsymm(isym),nsymm(isym),'u')
      call mxmaa(q(iupos(isym)),nsymm(isym),1,
     1 q(izpos(isym)),1,nsymm(isym),
     2 q(iwork2),1,nsymm(isym), nsymm(isym),nsymm(isym),nprm(isym) )
      call fmove (q(iwork2),q(izpos(isym)),nsymm(isym)*nprm(isym))
c...  contrib to second order energy
      enext = enext + dsum(nprm(isym),q(iwork2),nsymm(isym)+1)
 300  continue
c      print*,'energy after ub term  ',enext
      enext = enext - ddot(ic1d,gam,1,zint,1)
c      print*,enext
      ijij = 0
      do 310 ij=1,nact**2
      ijij = ijij + ij
310   enext = enext + 0.5d0*gam(ijij)*zint(ijij)
c      print*,enext
      end if
c
      call vclr(sigma,1,nrot)
      call orbasm (q(1),izpos,sigma,isignh)
c
c...  double on diagonal of transformed 2elec
      if (isigma.eq.1) then
      if (damp.ne.1.0d0) then
       call dscal(ne,damp,zint1,1)
           endif
      ijij = 0
      do 320 ij=1,nact**2
      ijij = ijij + ij
320   zint1(ijij) = zint1(ijij)*2.0d0
      if (iwrnr.eq.1) then
c...  must subtract (ij/kl) since otherwise got twice.
      do 330 ijkl=1,ic1d
330   zint1(ijkl) = zint1(ijkl) - zint(ijkl)
c     ee = hamilt(gam,zint1) + potnuc
c      print*,'gamma.zint1 = ',ee
      end if
      end if
      call accnt(' ',2)
      call corlsr (ibase)
      return
      end
      subroutine orb2 (q,vec,zint,zint1,loadz)
      implicit REAL  (a-h,o-z)
      parameter (maxinc=5000)
INCLUDE(common/sizes)
INCLUDE(common/multic)
      common/coruse/ ncruse
INCLUDE(common/count)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
      common /intbu2/ iblff2,intpo2,intfi2,intmo2,intnw2,intb2,g2(511)
      common /intbuf/ intpos,intfil,intmod
INCLUDE(common/mcff)
      common /cpos / iupos(8),iapos(8)
      common/lsort/v(maxorb),ijpos(8),ikpos(8),izpos(8),ixpos(8),w(31)
     1 ,iprev(8),jprev(8),idum(8),impos(8),inpos(8)
     2 ,iypos(8)
      dimension q(*),vec(nrot),zint(ne),zint1(ne)
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
c     timeio=0
      ibase = icorr(0)
c      call second(cporb2)
      call accnt('orb2',2)
      enext = potnuc+efreez
      do 10 isym=1,nirrr
      ijpos(isym) = icorr(maxbas**2)
      ikpos(isym) = icorr(maxbas**2)
      izpos(isym) = icorr(nsymm(isym)*nprm(isym))
      ixpos(isym) = icorr(nsymm(isym)**2)
      iypos(isym) = ixpos(isym)
      if (ideltr.eq.0) iypos(isym)=icorr(nsymm(isym)**2)
      if (isigma.eq.1) impos(isym) = icorr(nact**2)
      if (isigma.eq.1) inpos(isym) = icorr(nact**2)
10    continue
      ilpos = icorr(maxbas**2)
      iltpos = icorr(maxbas**2)
      iwork1 = icorr(2*maxbas**2)
      iwork0 = iwork1
      call vclr(q(ibase),1,iwork1+maxbas**2-ibase)
c
c...  load the x or t=exp(x)-1 matrix
      call xasm(q(1),ixpos,vec)
c...  load u for non-linear
      do 20 isym=1,nirrr
      n=nsymm(isym)
      if(n.eq.0) go to 20
      iu=iupos(isym)
      ix=ixpos(isym)
      iy=iypos(isym)
      if(ideltr.eq.0) then
c...  form u(r)
      call vclr(q(iy),1,n**2)
      call vfill(1.0d0,q(iy),n+1,n)
      call vadd(q(iy),1,q(ix),1,q(iy),1,n**2)
      call mxmaa(q(iu),1,n,q(iy),1,n,q(iwork1),1,n,n,n,n)
       call dcopy(n**2,q(iwork1),1,q(iy),1)
      call vclr(q(ix),1,n**2)
      call vfill(-1.0d0,q(ix),n+1,n)
      call vadd(q(iy),1,q(ix),1,q(ix),1,n**2)
c...  replace current u
      call fmove(q(iy),q(iu),n**2)
      else
c...  form u(r)*dr
      call mxmaa(q(iu),1,n,q(ix),1,n,q(iwork1),1,n,n,n,n)
      call fmove(q(iwork1),q(iy),n**2)
      end if
c...  now we have y=u*dr (ideltr.ne.0)  or u(r) (ideltr.eq.0)
c...              x=u*dr (ideltr.ne.0)  or t(r) (ideltr.eq.0)
c...  x overlays y for ideltr=1
20    continue
c
      left=min(0,icorrm()-maxbas**2)
      increm=min(maxinc,left)
      junk = icorr(-increm)
      left = left-increm
      maxcor = iwork1+increm
c
c... loop over integrals
      if(isigma.eq.1) call vclr(zint1,1,ne)
      intpos = 0
      intfil = num6
ctest iblf   = iblk6
      intpo2 = 0
      intfi2 = num4
ctest iblff2 = iblk4
      do 30 is=1,nirrr
30    iprev(is)=0
      do 210 i=1,nprim
      isymi = itype(i)
      niprev = iprev(isymi)*nsymm(isymi)
      nsymi=nsymm(isymi)
      ncori=ncor(isymi)
      iii=i-ncoremc
      if(iii.gt.0) iil=ilifa(iii)
      ixi=ixpos(isymi)
      izi=izpos(isymi)
      iji=ijpos(isymi)
      iki=ikpos(isymi)
      do 40 is=1,nirrr
40    jprev(is)=0
      do 200 j=1,i
      isymj = itype(j)
      nsymj=nsymm(isymj)
      njprev = jprev(isymj)*nsymj
      isymij = mults(isymi,isymj)
      nsymij = nsymi*nsymj
      jjj=j-ncoremc
      if(jjj.gt.0) jjl=ilifa(jjj)
      ixj=ixpos(isymj)
      izj=izpos(isymj)
      ikj=ikpos(isymj)
c      call second(cpu)
      if(loadz.eq.1) then
      call loadjk(q(1),ijpos,ikpos,i,j,zint)
      else
      call loadop(q(1),ijpos,ikpos,i,j)
      end if
c      call second(cpu1)
c      timeio=timeio+cpu1-cpu
c
      if (j.gt.ncoremc) goto 80
c...  j is core, so form l matrix
      call trnsps(q(ikj),q(ilpos),nsymj,nsymi)
      if (isignh.eq.1) then
       call dscal(nsymij,-1.0d0,q(ilpos),1)
           endif
      call daxpy(nsymij,4.0d0,q(iki),1,q(ilpos),1)
      call vsub(q(ilpos),1,q(iji),1,q(ilpos),1,nsymij)
      if(i.ne.j)
     * call trnsps(q(ilpos),q(iltpos),nsymi,nsymj)
c
      if (i.gt.ncoremc) goto 60
c...  ij = core-core: z(p,j) = 2*l(p,q)*x(q,j)
c...  if i.ne.j       z(q,j) = 2*l(q,p)*x(p,j)
      call fmove (q(ixj+njprev),q(iwork1),nsymj)
      call dscal(nsymj,2.0d0,q(iwork1),1)
      enext = enext + q(izi+niprev+iprev(isymi))
      call mxmb (q(ilpos),1,nsymi, q(iwork1),1,0,
     1           q(izi+niprev),1,0, nsymi,nsymj,1)
      enext = enext - q(izi+niprev+iprev(isymi))
      fac = -0.5d0
      if (i.ne.j) fac = -1.0d0
      enext = enext +fac*q(ilpos+iprev(isymi)+jprev(isymj)*nsymi)
      if (i.eq.j) goto 50
      call fmove (q(ixi+niprev),q(iwork1),nsymi)
      call dscal(nsymi,2.0d0,q(iwork1),1)
      enext = enext + q(izj+njprev+jprev(isymj))
      call mxmb (q(iltpos),1,nsymj, q(iwork1),1,0,
     1           q(izj+njprev),1,0, nsymj,nsymi,1)
      enext = enext - q(izj+njprev+jprev(isymj))
50    continue
      goto 200
c
c...  i.gt.ncoremc.ge.j
c...  z(p,j) =l(p,q) * x(q,j)?
60    call mxmaa(q(ilpos),1,nsymi,
     1 q(ixj+njprev),1,nsymj,q(iwork1),1,0,nsymi,nsymj,1)
c
      if(isigma.eq.1) then
      iic=ncori
      do 70 k=1,nact
      if (itypea(k).ne.isymi) goto 70
      val = -0.25d0 * q(ilpos + jprev(isymj)*nsymi + iic) +
     1  ddot(nsymi,q(iwork1),1,q(ixi+iic*nsymi),1)
      zint1 (ic1d+iil+k) = zint1(ic1d+iil+k) + val
      zint1 (ic1d+ilifa(k)+iii) = zint1(ic1d+ilifa(k)+iii) + val
      iic=iic+1
70    continue
      end if
c
      iwork1=iwork1+nsymi
      call fmove(q(ilpos+jprev(isymj)*nsymi+ncori),q(iwork1),
     1   nactt(isymi))
      iwork1=iwork1+nactt(isymi)
c
      if(iwork1.gt.maxcor) then
      if(left.gt.0) then
      increm = min(maxinc,left)
      left = left - increm
      junk = icorr(-increm)
      maxcor = maxcor + increm
      else
      call intou2(q(iwork0),iwork1-iwork0)
      ncruse=max(ncruse,iwork1)
      iwork1=iwork0
      end if
      end if
c
c...  z(p,j) = l(q,p) * x(q,t)
      call mxmaa(q(iltpos),1,nsymj,q(ixi+ncori*nsymi),1,nsymi,
     1   q(iwork1),1,nsymj,nsymj,nsymi,nactt(isymi))
      iwork1=iwork1+nsymj*nactt(isymi)
      if(iwork1.gt.maxcor) then
      if(left.gt.0) then
      increm = min(maxinc,left)
      left = left - increm
      junk = icorr(-increm)
      maxcor = maxcor + increm
      else
      call intou2(q(iwork0),iwork1-iwork0)
      ncruse=max(ncruse,iwork1)
      iwork1=iwork0
      end if
      end if
      goto 200
c
c...  i,j active
80    ivx = iil+jjj
      ixv = jjl+iii
      if (isigma.eq.1) then
      do 90 isymt=1,nirrr
      call vclr(q(impos(isymt)),1,nact**2)
90    continue
      end if
c
      do 140 isymt=1,nirrr
      isymu = mults(isymt,isymij)
      ncu=nactt(isymu)
      if(ncu.eq.0) goto 140
      nct=nactt(isymt)
      if(nct.eq.0) goto 140
      nnu=nsymm(isymu)
      nnt=nsymm(isymt)
      iju=ijpos(isymu)
      iku=ikpos(isymu)
      ixu=ixpos(isymu)
c     ixt=ixpos(isymt)
      iyu=iypos(isymu)
      iyt=iypos(isymt)
      icu=ncor(isymu)
      ict=ncor(isymt)
c
c...   c o u l o m b   i n t e g r a l s
c
c...  z(p,t) =j(p,q)*u(q,t)?
      call mxmaa(q(iju),1,nnu,q(iyt+ict*nnt),1,nnt,
     1 q(iwork1),1,nnu,nnu,nnt,nct)
c
      if (isigma.eq.1) then
c...  load these transformed integrals for the ci part
      call mxmb (q(iyu+icu*nnu),nnu,1
     1,q(iwork1),1,nnu, q(impos(isymu)),1,ncu,ncu,nnu,nct)
      end if
c
      iwork1=iwork1+nnu*nct
      if(iwork1.gt.maxcor) then
      if(left.gt.0) then
      increm = min(maxinc,left)
      left = left - increm
      junk = icorr(-increm)
      maxcor = maxcor + increm
      else
      call intou2(q(iwork0),iwork1-iwork0)
      ncruse=max(ncruse,iwork1)
      iwork1=iwork0
      end if
      end if
c
c...  e x c h a n g e   i n t e g r a l s
c
c...  z(p,t) =k(p,q)*x(q,t)?
      call mxmaa(q(iku),1,nnu,q(ixpos(isymt)+ict*nnt),1,nnt,
     1 q(iwork1),1,nnu, nnu,nnt,nct)
c
c...  store t(dagger).k.t = n in the case of non-linear with coupling
      if (isigma.eq.1)  then
      call mxmaa(q(ixu+icu*nnu),nnu,1,
     1  q(iwork1),1,nnu, q(inpos(isymu)),1,ncu,ncu,nnu,nct)
      if (i.eq.j) then
      iiin = inpos(isymu)
      call dscal(ncu*nct,0.5d0,q(iiin),1)
          endif
      end if
c
      iwork1=iwork1+nnu*nct
      if(iwork1.gt.maxcor) then
      if(left.gt.0) then
      increm = min(maxinc,left)
      left = left - increm
      junk = icorr(-increm)
      maxcor = maxcor + increm
      else
      call intou2(q(iwork0),iwork1-iwork0)
      ncruse=max(ncruse,iwork1)
      iwork1=iwork0
      end if
      end if
c
      if (i.ne.j) then
c...  z(p,u) =k(q,p)*x(q,u)?
      call mxmaa(q(iku),nnu,1,q(ixu+icu*nnu),1,nnu,
     1 q(iwork1),1,nnt, nnt,nnu,ncu)
      iwork1=iwork1+nnt*ncu
      if(iwork1.gt.maxcor) then
      if(left.gt.0) then
      increm = min(maxinc,left)
      left = left - increm
      junk = icorr(-increm)
      maxcor = maxcor + increm
      else
      call intou2(q(iwork0),iwork1-iwork0)
      ncruse=max(ncruse,iwork1)
      iwork1=iwork0
      end if
      end if
      end if
c
140   continue
c
c..   put the transformed integrals into zint1 as appropriate
      do 190 isymt=1,nirrr
      isymu = mults(isymt,isymij)
      if(nactt(isymt).eq.0) goto 185
      if(nactt(isymu).eq.0) goto 185
      if (isigma.eq.1) then
      itu = impos(isymt) - 1
      do 160 iu=1,nact
      if (itypea(iu).ne.isymu) goto 160
      do 150 it=1,nact
      if (itypea(it).ne.isymt) goto 150
      itu = itu + 1
      int1 = ind (ivx,ilifa(it)+iu)
      int2 = ind (ixv,ilifa(it)+iu)
      zint1(int1) = zint1(int1) + q(itu)
      if (int1.ne.int2) zint1(int2) = zint1(int2) + q(itu)
150   continue
160   continue
c
      itu = inpos(isymt) - 1
      do 180 iu=1,nact
      if (itypea(iu).ne.isymu) goto 180
      iljj=ilifa(iu)+jjj
      jjiu=jjl+iu
      do 170 it=1,nact
      if (itypea(it).ne.isymt) goto 170
      itu = itu + 1
      int1 = ind (ilifa(it)+iii,iljj)
      zint1(int1) = zint1(int1) + q(itu)
      int1 = ind (iil+it,iljj)
      zint1(int1) = zint1(int1) + q(itu)
      int1 = ind (ilifa(it)+iii,jjiu)
      zint1(int1) = zint1(int1) + q(itu)
      int1 = ind (iil+it,jjiu)
      zint1(int1) = zint1(int1) + q(itu)
170   continue
180   continue
      end if
c..   contributions to core
185   if (isymt.eq.isymu .and. ncor(isymt).gt.0) then
      nct=ncor(isymt)
      nnt=nsymm(isymt)
      ijt=ijpos(isymt)
      ikt=ikpos(isymt)
      iyt=iypos(isymt)
      call trnsps(q(ikt),q(iwork1),nnt,nnt)
      call vadd(q(ikt),1,q(iwork1),1,q(ikt),1,nnt**2)
      call daxpy(nnt**2,-4.0d0,q(ijt),1,q(ikt),1)
      call mxmaa(q(ikt),1,nnt,
     1 q(iyt),1,nnt,q(iwork1),1,nnt,nnt,nnt,nct)
c
      if (isigma.eq.1) then
      val = -0.5d0*ddot(nct*nnt,q(iwork1),1,q(iyt),1)
      zint1(ic1d+iil+jjj) = zint1(ic1d+iil+jjj) + val
      if (i.ne.j) zint1(ic1d+jjl+iii) = zint1(ic1d+jjl+iii) + val
      end if
c
      iwork1=iwork1+nnt*nct
      if(iwork1.gt.maxcor) then
      if(left.gt.0) then
      increm = min(maxinc,left)
      left = left - increm
      junk = icorr(-increm)
      maxcor = maxcor + increm
      else
      call intou2(q(iwork0),iwork1-iwork0)
      ncruse=max(ncruse,iwork1)
      iwork1=iwork0
      end if
      end if
c
      end if
190   continue
c
200   jprev(isymj) = jprev(isymj) + 1
210   iprev(isymi) = iprev(isymi) + 1
c
c...  process core hamiltonian
      if(loadz.eq.1) then
      call loadjk (q(1),ijpos,ikpos,0,0,zint)
      else
      call loadop (q(1),ijpos,ikpos,0,0)
      end if
      do 260 isymt=1,nirrr
      nnt=nsymm(isymt)
      nct=ncor(isymt)
      iyt=iypos(isymt)
      if (nct.gt.0) then
      call mxmaa(q(ijpos(isymt)),1,nnt,
     1 q(iyt),1,nnt, q(iwork1),1,nnt,nnt,nnt,nct)
      iwork1=iwork1+nnt*nct
      if(iwork1.gt.maxcor) then
      if(left.gt.0) then
      increm = min(maxinc,left)
      left = left - increm
      junk = icorr(-increm)
      maxcor = maxcor + increm
      else
      call intou2(q(iwork0),iwork1-iwork0)
      ncruse=max(ncruse,iwork1)
      iwork1=iwork0
      end if
      end if
c
      end if
c
      if(nactt(isymt).eq.0) goto 260
c...  z(p,t) =j(p,q)*u(q,t)?
      call mxmaa(q(ijpos(isymt)),1,nnt,
     1 q(iyt+nct*nnt),1,nnt,q(iwork1),1,nnt,nnt,nnt,nactt(isymt))
c
      if (isigma.eq.1) then
c...  load transformed integrals for ci part
      itu=iwork1+nct-1
      itt=iwork1-nnt
      do 250 it=1,nact
      if (itypea(it).ne.isymt) goto 250
      itt = itt + nnt
      iuu = iyt + nnt*(nct-1)
      do 240 iu=1,nact
      if (itypea(iu).ne.isymt) goto 240
      iuu = iuu + nnt
      itu = itu + 1
      val=ddot(nnt,q(itt),1,q(iuu),1)
      intad = ic1d + ilifa(it)+iu
      zint1(intad) = zint1(intad) + val
240   continue
      itu = itu + nsecc(isymt)+nct
250   continue
      end if
c
      iwork1=iwork1+nnt*nactt(isymt)
      if(iwork1.gt.maxcor) then
      if(left.gt.0) then
      increm = min(maxinc,left)
      left = left - increm
      junk = icorr(-increm)
      maxcor = maxcor + increm
      else
      call intou2(q(iwork0),iwork1-iwork0)
      ncruse=max(ncruse,iwork1)
      iwork1=iwork0
      end if
      end if
c
260   continue
c..   save core part of z-matrix
      do 280 isym=1,nirrr
      if(nprm(isym).eq.0) goto 280
      call dcopy(nsymm(isym)*nprm(isym),q(izpos(isym)),1,q(iwork1),1)
      iwork1=iwork1+nsymm(isym)*nprm(isym)
      if(iwork1.gt.maxcor) then
      if(left.gt.0) then
      increm = min(maxinc,left)
      left = left - increm
      junk = icorr(-increm)
      maxcor = maxcor + increm
      else
      call intou2(q(iwork0),iwork1-iwork0)
      ncruse=max(ncruse,iwork1)
      iwork1=iwork0
      end if
      end if
280   continue
c
      if(iwork1.gt.iwork0) call intou2(q(iwork0),iwork1-iwork0)
      ncruse=max(ncruse,iwork1)
      call inten2
c
c...  double on diagonal of transformed 2elec
      if (isigma.eq.1) then
      ijij = 0
      do 320 ij=1,nact**2
      ijij = ijij + ij
320   zint1(ijij) = zint1(ijij)*2.0d0
c...  must subtract (ij/kl) since otherwise got twice.
      do 330 ijkl=1,ic1d
330   zint1(ijkl) = zint1(ijkl) - zint(ijkl)
      end if
      call accnt(' ',2)
      call corlsr (ibase)
_IF1()c      call second(cpu)
_IF1()c      cporb2=cpu-cporb2
_IF1()c      oi=timeio/cporb2*100
_IF1()c      write(iwrite,141) cporb2,timeio,oi
_IF1()c141   format(' cpu times orb2=',f10.6,'  i/o=',f10.6,' = ',f10.2,' %')
      return
      end
      subroutine orb3 (q,vec,sigma,zint,gam)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
      common /intbu2/ iblf,intpos,intfil,intmod,intnw,intb2,g(511)
      common /cpos / iupos(8),iapos(8)
      common/lsort/v(maxorb),ijpos(8),ikpos(8),izpos(8),ixpos(8),w(31)
     1 ,iprev(8),jprev(8),idum(8),impos(8),inpos(8)
     2 ,iypos(8)
      dimension q(*),vec(nrot),sigma(nrot),gam(ne),zint(ne)
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
      ibase = icorr(0)
      call accnt('orb3',2)
c     timeio=0
c     call second(cpu)
      do 10 isym=1,nirrr
      ixpos(isym) = icorr(nsymm(isym)**2)
10    izpos(isym) = icorr(nsymm(isym)*nprm(isym))
      iwork1 = icorr(maxbas**2)
      iwork2 = icorr(maxbas**2)
      ippos = icorr(nact**2)
      iqpos = icorr(nact**2)
      call vclr(q(ibase),1,iqpos+nact**2-ibase)
c
c
c... loop over integrals
      intpos = 0
ctest iblf   = iblk4
      intmod = 1
      do 30 is=1,nirrr
30    iprev(is)=0
      do 210 i=1,nprim
      isymi = itype(i)
      if(i.le.ncoremc) goto 210
c     niprev = iprev(isymi)*nsymm(isymi)
      nsymi=nsymm(isymi)
      ncori=ncor(isymi)
      iii=i-ncoremc
      if(iii.gt.0) iil=ilifa(iii)
      izi=izpos(isymi)
      do 40 is=1,nirrr
40    jprev(is)=0
      do 200 j=1,i
      isymj = itype(j)
      nsymj = nsymm(isymj)
      njprev = jprev(isymj)*nsymj
      isymij = mults(isymi,isymj)
      jjj=j-ncoremc
      if(jjj.gt.0) jjl=ilifa(jjj)
      izj=izpos(isymj)
c
      if(j.gt.ncoremc) goto 80
c
c...  core-valence
      ii=1
      do 70 k=1,nact
      if(itypea(k).ne.isymi) goto 70
      v(ii) = gam(ic1d+iii+ilifa(k))
      ii=ii+1
70    continue
c
      enext=enext+dsum(nprm(isymi),q(izi),nsymi+1)
      call intin2(q(iwork1),nsymi)
      call mxmb (q(iwork1),1,0,v,0,1,q(izi+ncori*nsymi),1,nsymi,
     2 nsymi,1,nactt(isymi))
      call intin2(q(iwork1),nactt(isymi))
      enext=enext-dsum(nprm(isymi),q(izi),nsymi+1)
     1  -0.5d0*ddot(nactt(isymi),q(iwork1),1,v,1)
c
      enext = enext + q(izj+njprev+jprev(isymj))
      call intin2(q(iwork1),nsymj*nactt(isymi))
      call mxmb(q(iwork1),1,nsymj,
     1   v,1,0,q(izj+njprev),1,0,nsymj,nactt(isymi),1)
      enext = enext - q(izj+njprev+jprev(isymj))
      goto 200
c
c...  i,j active
80    ivx = iil+jjj
c     ixv = jjl+iii
      do 140 isymt=1,nirrr
      isymu = mults(isymt,isymij)
      ncu=nactt(isymu)
      if(ncu.eq.0) goto 140
      nct=nactt(isymt)
      if(nct.eq.0) goto 140
      nnu=nsymm(isymu)
      nnt=nsymm(isymt)
      icu=ncor(isymu)
      ict=ncor(isymt)
      izu=izpos(isymu)
      izt=izpos(isymt)
c
c...  set up density matrices for this symmetry
c...  p(t,u) = gam(u,t,i,j) + gam(u,t,j,i)
c...  q(t,u) = gam(u,i,j,t) + gam(u,i,t,j)
      itu=-1
      do 110 iu=1,nact
      if (isymu.ne.itypea(iu)) goto 110
      iliu=ilifa(iu)
      ilii=iliu+iii
      if(i.ne.j) then
      do 100 it=1,nact
      if (isymt.ne.itypea(it)) goto 100
      itu = itu + 1
      q(ippos+itu) = gam(ind(ivx,ilifa(it)+iu))
     1             + gam(ind(ivx,iliu+it))
      q(iqpos+itu) = gam(ind(ilii,jjl+it))
     1             + gam(ind(ilii,ilifa(it)+jjj))
100   continue
      else
      do 101 it=1,nact
      if (isymt.ne.itypea(it)) goto 101
      itu = itu + 1
      q(ippos+itu) = gam(ind(ivx,ilifa(it)+iu))
      q(iqpos+itu) = gam(ind(ilii,jjl+it))
     1             + gam(ind(ilii,ilifa(it)+jjj))
101   continue
      end if
110   continue
c
c...  c o u l o m b   i n t e g r a l s
c
      call intin2(q(iwork1),nnu*nct)
      call mxmb (q(iwork1),1,nnu,q(ippos),1,nct,
     1 q(izu+icu*nnu),1,nnu,nnu,nct,ncu)
c
c...  e x c h a n g e   i n t e g r a l s
c
      enext = enext + dsum(nprm(isymu),q(izu),nnu+1)
      call intin2(q(iwork1),nnu*nct)
      call mxmb (q(iwork1),1,nnu,q(iqpos),1,nct,
     1 q(izu+icu*nnu),1,nnu,nnu,nct,ncu)
      enext = enext - dsum(nprm(isymu),q(izu),nnu+1)
c
      if (i.ne.j) then
      enext = enext + dsum(nprm(isymt),q(izt),nnt+1)
      call intin2(q(iwork1),nnt*ncu)
      call mxmb (q(iwork1),1,nnt,q(iqpos),nct,1,
     1 q(izt+ict*nnt),1,nnt,nnt,ncu,nct)
      enext = enext - dsum(nprm(isymt),q(izt),nnt+1)
      end if
c
140   continue
c
c..   contributions to core
      if(ncoremc.eq.0) goto 200
      do 190 isymt=1,nirrr
      if(ncor(isymt).eq.0) goto 190
      if(mults(isymt,isymij).ne.isymt) goto 190
      nct=ncor(isymt)
      nnt=nsymm(isymt)
      call intin2(q(iwork1),nnt*nct)
      val = - gam(ic1d+iil+jjj)
      if (i.eq.j) val=0.5d0*val
      call daxpy(nnt*nct
     +  ,val,q(iwork1),1,q(izpos(isymt)),1)
190   continue
c
200   jprev(isymj) = jprev(isymj) + 1
210   iprev(isymi) = iprev(isymi) + 1
c
c...  process core hamiltonian
      do 260 isymt=1,nirrr
      nnt=nsymm(isymt)
      nct=ncor(isymt)
      if (nct.gt.0) then
      call intin2(q(iwork1),nnt*nct)
      call daxpy(nnt*nct
     +  ,2.0d0,q(iwork1),1,q(izpos(isymt)),1)
      end if
c
      if(nactt(isymt).eq.0) goto 260
      nac=nactt(isymt)
      itu=ippos-1
      do 230 it=1,nact
      if (itypea(it).ne.isymt) goto 230
      do 220 iu=1,nact
      if (itypea(iu).ne.isymt) goto 220
      itu=itu+1
      q(itu) = gam(ic1d + ilifa(it) + iu)
220   continue
230   continue
c
      call intin2(q(iwork1),nnt*nac)
      call mxmb (q(iwork1),1,nnt,q(ippos),1,nac,
     1 q(izpos(isymt)+nct*nnt),1,nnt,nnt,nac,nac)
c
260   continue
c
c...  form u(dagger)*b
c
      if(ideltr.ne.0) call xasm(q(1),ixpos,vec)
      do 300 isym=1,nirrr
      m=nprm(isym)
      if(m.eq.0) goto 300
      n=nsymm(isym)
      ia=iapos(isym)
      iu=iupos(isym)
      iz=izpos(isym)
      call intin2(q(iwork1),n*m)
      call vadd(q(iz),1,q(iwork1),1,q(iz),1,n*m)
      call mxmaa(q(iu),n,1,q(iz),1,n,q(iwork1),1,n,n,n,m)
      call dcopy(n*m,q(iwork1),1,q(iz),1)
      if(ideltr.eq.0) then
      call dcopy(n*m,q(iwork1),1,q(ia),1)
      else
      ix=ixpos(isym)
      call vadd(q(iz),1,q(ia),1,q(iz),1,n*m)
c...  form -0.5*[a+a(degger)]
      call expda(q(ia),q(iwork2),n,m)
c...  form -0.5*[a+a(degger)]*dr
      call mxmaa(q(iwork2),1,n,q(ix),1,n,q(iwork1),1,n,n,n,n)
c...  antisymmetrize external part
      call reduca(q(iwork1),q(iz),n,m)
      end if
c...  contrib to second order energy
      enext = enext + dsum(m,q(iz),n+1)
300   continue
      enext = enext - ddot(ic1d,gam,1,zint,1)
c
      ijij = 0
      do 310 ij=1,nact**2
      ijij = ijij + ij
310   enext = enext + 0.5d0*gam(ijij)*zint(ijij)
c
      call vclr(sigma,1,nrot)
      call orbasm (q(1),izpos,sigma,isignh)
c
_IF1()c      call second(cpu1)
_IF1()c      cporb3=cpu1-cpu
_IF1()c      oi=timeio/cporb3*100
_IF1()c      write(iwrite,141) cporb3,timeio,oi
_IF1()c141   format(' cpu times orb3=',f10.6,'  i/o=',f10.6,' = ',f10.2,' %')
      call accnt(' ',2)
      call corlsr (ibase)
      return
      end
      subroutine loadi (q,qj)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mcff)
INCLUDE(common/multic)
INCLUDE(common/syminf)
      common /intbuf/ intpos,intfil,intmod
INCLUDE(common/mcscra)
      dimension q(*),qj(*)
      intpos=0
      intfil=num6
ctest iblf  =iblk6
      ibuff = icorr(max(maxbas**2,ltrimo))
      do 50 i=1,nprim
      do 50 j=1,i
      ij=jadr(i,j)
      isymij = mults(itype(i),itype(j))
      do 30 isyma=1,nirrr
      isymb = mults(isyma,isymij)
      if (isymb-isyma) 10,20,30
10    call intin (q(ibuff),nsymm(isyma)*nsymm(isymb))
      ib = ibuff
      do 11 ii=1,nprm(isymb)
      call fmove (q(ib),qj(ij),nprm(isyma))
      ib = ib + nsymm(isyma)
11    ij = ij + nprm(isyma)
      goto 30
20    call intin (q(ibuff),(nsymm(isyma)*(nsymm(isyma)+1))/2)
      call fmove (q(ibuff),qj(ij),(nprm(isyma)*(nprm(isyma)+1))/2)
      ij = ij + (nprm(isyma)*(nprm(isyma)+1))/2
30    continue
c.....skip exchange
      do 40 isyma=1,nirrr
      isymb = mults(isyma,isymij)
40    call intin (q(ibuff),nsymm(isyma)*nsymm(isymb))
50    continue
c....skip core hamiltonian
      call intin (q(ibuff),ltrimo)
c....bare h0
      call intin (q(ibuff),ltrimo)
      ij=ifc
      ib=ibuff
      do 70 isym=1,nirrr
      l1 = (nprm(isym)*(nprm(isym)+1))/2
      l2 = (nsymm(isym)*(nsymm(isym)+1))/2
      call fmove (q(ib),qj(ij),l1)
      ib=ib+l2
70    ij=ij+l1
      call corlsr (ibuff)
      return
      end
      subroutine loadjk (q,ijpos,ikpos,i,j,zint)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
      dimension q(*),ijpos(8),ikpos(8)
      dimension zint(ne)
      ind(k,l)=max(k,l)*(max(k,l)-1)/2+min(k,l)
      isingl = ne + (ind(i-ncoremc,j-ncoremc)-1)*nbasis*nact
      if (i.eq.0) isingl = ne + ind(nact,nact)*nact*nbasis
      isymij = 1
      if (i.gt.0) isymij = mults(itype(i),itype(j))
      ibuff = icorr(-maxbas**2)
c
c...  skip h0
      if(i.lt.0) then
      if(ltrimo.gt.maxbas**2) junk=icorr(-ltrimo+maxbas**2)
      call intin(q(ibuff),ltrimo)
      end if
c...  j matrix
      do 30 isyma=1,nirrr
      isymb = mults(isyma,isymij)
      if (isymb-isyma) 10,20,30
10    call intin (q(ijpos(isyma)),nsymm(isyma)*nsymm(isymb))
      ija=ijpos(isyma)
      ijb=ijpos(isymb)
      call trnsps(q(ija),q(ijb),nsymm(isyma),nsymm(isymb))
      goto 30
20    call intin (q(ibuff),(nsymm(isyma)*(nsymm(isyma)+1))/2)
      call square (q(ijpos(isyma)),q(ibuff),nsymm(isyma),nsymm(isyma))
30    continue
c
      if (i.gt.0) then
c...  exchange
      do 40 isyma=1,nirrr
      isymb = mults(isyma,isymij)
      call intin (q(ibuff),nsymm(isyma)*nsymm(isymb))
      ikb=ikpos(isymb)
      call trnsps(q(ibuff),q(ikb),nsymm(isyma),nsymm(isymb))
 40   continue
      end if
      if (j.gt.ncoremc) then
c...  load coulomb integrals j(t,u) -> zint(t,u,ij)
      ij = ilifa(i-ncoremc)+j-ncoremc
      ji = ilifa(j-ncoremc)+i-ncoremc
      do 60 isyma=1,nirrr
      isymb = mults (isyma,isymij)
      iad = ijpos(isyma)+ncor(isymb)*nsymm(isyma)-1
      do 50 ib=ncor(isymb)+1,nprm(isymb)
      ibb = iorbsm(ib-1+istart(isymb)) - ncoremc
      do 55 ia=ncor(isyma)+1,nprm(isyma)
      iaa = iorbsm(ia-1+istart(isyma)) - ncoremc
      zint(ind(ij,ilifa(ibb)+iaa)) = q(iad+ia)
55    zint(ind(ji,ilifa(ibb)+iaa)) = q(iad+ia)
      if (isigma.eq.1.and.iwrnr.eq.0) then
      do 65 ia = 1,nsymm(isyma)
      iaa = iorbsm(ia-1+istart(isyma))
65    zint (isingl + iaa + (ibb-1)*nbasis) = q(iad+ia)
      end if
50    iad = iad + nsymm(isyma)
60    continue
      else if (i.eq.0) then
      do 61 isyma=1,nirrr
      iad = ijpos(isyma)+ncor(isyma)*nsymm(isyma)-1
      do 51 ib=ncor(isyma)+1,nprm(isyma)
      ibb = iorbsm(ib-1+istart(isyma)) - ncoremc
      do 56 ia=ncor(isyma)+1,nprm(isyma)
      iaa = iorbsm(ia-1+istart(isyma)) - ncoremc
56    zint(ilifa(ibb)+iaa+ic1d) = q(iad+ia)
      if (isigma.eq.1.and.iwrnr.eq.0) then
      do 66 ia=1,nsymm(isyma)
      iaa = iorbsm(ia-1+istart(isyma))
66    zint(isingl + iaa + (ibb-1)*nbasis) = q(iad+ia)
      end if
51    iad = iad + nsymm(isyma)
61    continue
      end if
      call corlsr (ibuff)
      return
      end
      subroutine loadop (q,ijpos,ikpos,i,j)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
      dimension q(*),ijpos(8),ikpos(8)
      isymij = 1
      if (i.gt.0) isymij = mults(itype(i),itype(j))
      ibuff = icorr(maxbas**2)
c
c...  j matrix
      do 30 isyma=1,nirrr
      isymb = mults(isyma,isymij)
      if (isymb-isyma) 10,20,30
10    call intin (q(ijpos(isyma)),nsymm(isyma)*nsymm(isymb))
      ija=ijpos(isyma)
      ijb=ijpos(isymb)
      call trnsps(q(ija),q(ijb),nsymm(isyma),nsymm(isymb))
      goto 30
20    call intin (q(ibuff),(nsymm(isyma)*(nsymm(isyma)+1))/2)
      call square (q(ijpos(isyma)),q(ibuff),nsymm(isyma),nsymm(isyma))
30    continue
c
      if (i.gt.0) then
c...  exchange
      do 40 isyma=1,nirrr
      isymb = mults(isyma,isymij)
      call intin (q(ibuff),nsymm(isyma)*nsymm(isymb))
      ikb=ikpos(isymb)
      call trnsps(q(ibuff),q(ikb),nsymm(isyma),nsymm(isymb))
 40   continue
      end if
      call corlsr (ibuff)
      return
      end
      subroutine mkzint (qj,zint)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/mcscra)
      dimension qj(*),zint(ne)
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      call vclr(zint,1,ne)
      do 51 i=1,nact
      do 50 j=1,i
      iad=jadr(i+ncoremc,j+ncoremc)-1
      ij = ilifa(i)+j
      ji = ilifa(j)+i
      isymi = itypea(i)
      isymj = itypea(j)
      isymij=mults(isymi,isymj)
      do 40 isymk=1,nirrr
      isyml = mults(isymk,isymij)
      if (isymk.lt.isyml) goto 40
      maxl=nprm(isymk)
      do 31 k=1,nprm(isyml)
      kk = iorbs(k,isyml)
      if (isymk.eq.isyml) maxl=k
      if(kk.eq.0) goto 31
      k1=ilifa(kk)
      do 30 l=1,maxl
      ll = iorbs(l,isymk)
      if(ll.eq.0) goto 30
      kl = k1+ll
      lk = ilifa(ll)+kk
      zint (ind(ij,kl)) = qj(iad+l)
      zint (ind(ij,lk)) = qj(iad+l)
      zint (ind(ji,kl)) = qj(iad+l)
      zint (ind(ji,lk)) = qj(iad+l)
30    continue
31    iad=iad+maxl
40    continue
50    continue
51    continue
      iad=ifc-1
      do 80 isymk=1,nirrr
      do 70 k=1,nprm(isymk)
      kk = iorbs(k,isymk)
      if(kk.eq.0) goto 70
      do 60 l=1,k
      ll = iorbs(l,isymk)
      if(ll.eq.0) goto 60
      zint (ic1d+ilifa(kk)+ll) = qj(iad+l)
      zint (ic1d+ilifa(ll)+kk) = qj(iad+l)
60    continue
70    iad=iad+k
80    continue
      return
      end
      subroutine orbasm (q,ifpos,vec,isign)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
      common /block / nosymm,nosym(8),iosym(8,maxorb)
      dimension ifpos(8)
      dimension q(*),vec(nrot)
      signh =  dfloat(isign)
      ivec=0
      do 40 isym=1,nirrr
c
      ired = 0
      itt = ncor(isym)-1
      do 20 it=1,nact
      if (itypea(it).eq.isym) itt = itt + 1
      ii = -1
      do 10 i=1,it-1+ncoremc
      if(itype(i).eq.isym) ii=ii+1
      if (itype(i).ne.itypea(it)) goto 10
      if (i.gt.ncoremc) then
      ired = ired + 1
      if (irottu(ired).eq.0) goto 10
      end if
      if (itype(i).ne.isym) goto 10
      ivec = ivec + 1
      vec(ivec) = vec(ivec) + signh * q(ifpos(isym)+itt+ii*nsymm(isym))
     1                              - q(ifpos(isym)+ii+itt*nsymm(isym))
      if(iosym(isym,itt+1).ne.iosym(isym,ii+1)) vec(ivec)=0
10    continue
20    continue
      do 30 ia=nprm(isym)+1,nsymm(isym)
      do 30 iqq=1,nprm(isym)
      ivec = ivec + 1
      vec(ivec) =vec(ivec)+signh*q(ifpos(isym)+ia-1+(iqq-1)*nsymm(isym))
30    if(iosym(isym,ia).ne.iosym(isym,iqq)) vec(ivec)=0
40    continue
      return
      end
      subroutine orblab (num,i1,j1,isym)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
      ivec=0
      do 40 isym=1,nirrr
      ired = 0
      j1 = ncor(isym)
      do 20 it=1,nact
      if (itypea(it).eq.isym) j1 = j1 + 1
      i1 = 0
      do 10 i=1,it-1+ncoremc
      if(itype(i).eq.isym) i1=i1+1
      if (itype(i).ne.itypea(it)) goto 10
      if (i.gt.ncoremc) then
      ired = ired + 1
      if (irottu(ired).eq.0) goto 10
      end if
      if (itype(i).ne.isym) goto 10
      ivec = ivec + 1
      if (ivec.eq.num) return
10    continue
20    continue
      do 30 j1=nprm(isym)+1,nsymm(isym)
      do 30 i1=1,nprm(isym)
      ivec = ivec + 1
      if (ivec.eq.num) return
30    continue
40    continue
      call caserr('illegal orbital rotation analysed')
      end
      subroutine xasm (q,ixpos,vec)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
      common /block / nosymm,nosym(8),iosym(8,maxorb)
      integer ixpos(8)
      dimension q(*),vec(nrot)
      sign= -  dfloat(isignh)
      ivec=0
      ntexp=max(2,ntexp)
      do 160 isym=1,nirrr
      ns=nsymm(isym)
      ni=nprm(isym)
      ns1=ns+1
      ns2=ns**2
      ni1=ni+1
      ni2=ni**2
      ired=0
      itext=ixpos(isym)
      itint=icorr(ni2)
      call vclr(q(itint),1,ni2)
      call vclr(q(itext),1,ns2)
      itt=ncor(isym)-1
      do 20 it=1,nact
      if(itypea(it).eq.isym)itt=itt+1
      ii=-1
      do 10 i=1,it-1+ncoremc
      if(itype(i).eq.isym) ii=ii+1
      if(itype(i).ne.itypea(it)) goto 10
      if(i.gt.ncoremc) then
      ired=ired+1
      if(irottu(ired).eq.0) goto 10
      end if
      if (itype(i).ne.isym) goto 10
      ivec=ivec+1
      if(iosym(isym,itt+1).ne.iosym(isym,ii+1)) vec(ivec)=0
      if(iuprod.eq.0) q(itext+itt+ii*ns)=vec(ivec)
      q(itint+itt+ii*ni)=vec(ivec)
10    continue
20    continue
c
      it1=itext-1-ns
      do 30 ia=ni1,ns
      do 30 iqq=1,ni
      ivec=ivec+1
      if(iosym(isym,ia).ne.iosym(isym,iqq)) vec(ivec)=0
30    q(it1+ia+ns*iqq)=vec(ivec)
      if(iuprod.ne.0) then
      it2=itint-1-ni
      do 40 i=1,ni
      do 40 j=1,i
40    q(it2+j+ni*i) = sign * q(it2+i+ni*j)
      end if
      do 50 i=1,ns
      do 50 j=1,i
50    q(it1+j+ns*i) = sign * q(it1+i+ns*j)
c
      if (iwrnr.eq.1.and.ideltr.eq.0) then
c..   exponentiate to form exp(x)-1
      zz = dasum(ns2,q(itext),1)
      if(iuprod.ne.0) zz=zz+dasum(ni2,q(itint),1)
      fac = zz
      nterm = 2
60    nterm = nterm + 1
      fac = fac*zz/ dfloat(nterm)
      if (fac.gt.1.0d-8) goto 60
      len=ns2
      if(iuprod.ne.0) len=max(len,ni2)
      iu = icorr(len)
      iv = icorr(len)
      call vclr(q(iu),1,ns2)
      nterm=min(nterm,ntexp)
c.....form t(ext)
      do 70 i=nterm,2,-1
      call vclr(q(iv),1,ns2)
      call daxpy(ns2,(1.0d0/ dfloat(i)),q(iu),1,q(iv),1)
      call vclr(q(iu),1,ns2)
      call vfill(1.0d0,q(iu),ns1,ns)
70    call mxmb (q(itext),1,ns, q(iv),1,ns, q(iu),1,ns, ns,ns,ns)
      call fmove (q(itext),q(iv),ns2)
      call mxmaa(q(iv),1,ns,q(iu),1,ns,q(itext),1,ns,ns,ns,ns)
      if(iuprod.ne.0) then
c.....form u(ext)
      ij=itext-ns1
      do 80 i=1,ns
80    q(ij+i*ns1)=q(ij+i*ns1)+1.0d0
      call vclr(q(iu),1,ni2)
c.....form u(int)
      do 90 i=nterm,2,-1
      call vclr(q(iv),1,ni2)
      call daxpy(ni2,(1.0d0/ dfloat(i)),q(iu),1,q(iv),1)
      call vclr(q(iu),1,ni2)
      call vfill(1.0d0,q(iu),ni1,ni)
90    call mxmb (q(itint),1,ni, q(iv),1,ni, q(iu),1,ni, ni,ni,ni)
      call fmove (q(itint),q(iv),ni2)
      call mxmaa(q(iv),1,ni,q(iu),1,ni,q(itint),1,ni,ni,ni,ni)
      ij=itint-ni1
      do 100 i=1,ni
100   q(ij+i*ni1)=q(ij+i*ni1)+1.0d0
      if(iuprod.gt.1) goto 130
c.....form u(int)*u(ext)
      call mxmaa(q(itint),1,ni,q(itext),1,ns,q(iu),1,ni,ni,ni,ns)
      jii=iu-1
      ji=itext-1
      do 120 j=1,ns
      do 110 i=1,ni
110   q(ji+i)=q(jii+i)
      ji=ji+ns
120   jii=jii+ni
      goto 140
c.....form u(ext)*u(int)
130   call mxmaa(q(itext),1,ns,q(itint),1,ni,q(iu),1,ns,ns,ni,ni)
      call dcopy(ns*ni,q(iu),1,q(itext),1)
140   ij=itext-ns1
c.....form t=u-1
      do 150 i=1,ns
150   q(ij+i*ns1)=q(ij+i*ns1)-1.0d0
      end if
      end if
      call corlsr(itint)
160   continue
      return
      end
      subroutine uasm (q,ixpos,vec)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
      common /block / nosymm,nosym(8),iosym(8,maxorb)
      integer ixpos(8)
      dimension q(*),vec(nrot)
      ivec=0
      do 160 isym=1,nirrr
      ns=nsymm(isym)
      ni=nprm(isym)
      ni1=ni+1
c     ns1=ns+1
      ired=0
      iu=ixpos(isym)
      itt=ncor(isym)-1
      do 20 it=1,nact
      if(itypea(it).eq.isym) itt=itt+1
      ii=-1
      do 10 i=1,it-1+ncoremc
      if(itype(i).eq.isym) ii=ii+1
      if(itype(i).ne.itypea(it)) goto 10
      if(i.gt.ncoremc) then
      ired=ired+1
      if(irottu(ired).eq.0) goto 10
      end if
      if (itype(i).ne.isym) goto 10
      ivec=ivec+1
      vec(ivec)=0.5d0*(q(iu+itt+ii*ns)-q(iu+ii+itt*ns))
      if(iosym(isym,itt+1).ne.iosym(isym,ii+1)) vec(ivec)=0
10    continue
20    continue
c
      do 30 ia=ni1,ns
      do 30 iqq=1,ni
      ivec=ivec+1
      vec(ivec)=0.5d0*(q(iu+(ia-1)+ns*(iqq-1))
     1     -q(iu+(iqq-1)+ns*(ia-1)))
30    if(iosym(isym,ia).ne.iosym(isym,iqq)) vec(ivec)=0
160   continue
c
      return
      end
      subroutine intu (q,ixpos,vec,anorm,iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common /block / nosymm,nosym(8),iosym(8,maxorb)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
      dimension q(*),ixpos(8),ijpos(8)
      dimension vec(nrot)
      sign= -  dfloat(isignh)
      ivec=0
      ntexp=max(2,ntexp)
      anorm=0
      do 60 isym=1,nirrr
      ired=0
      call vclr(q(ixpos(isym)),1,nprm(isym)**2)
      itt=ncor(isym)-1
      do 20 it=1,nact
      if(itypea(it).eq.isym)itt=itt+1
      ii=-1
      do 10 i=1,it-1+ncoremc
      if(itype(i).eq.isym) ii=ii+1
      if(itype(i).ne.itypea(it)) goto 10
      if(i.gt.ncoremc) then
      ired=ired+1
      if(irottu(ired).eq.0) goto 10
      end if
      if (itype(i).ne.isym) goto 10
      ivec=ivec+1
      if(iosym(isym,itt+1).ne.iosym(isym,ii+1)) vec(ivec)=0
      q(ixpos(isym)+itt+ii*nprm(isym))=vec(ivec)
      anorm=anorm+vec(ivec)**2
      vec(ivec)=-vec(ivec)
10    continue
20    continue
c
      call vclr(vec(ivec+1)
     +  ,1,(nsymm(isym)-nprm(isym))*nprm(isym))
      ivec=ivec+(nsymm(isym)-nprm(isym))*nprm(isym)
      do 30 i=1,nprm(isym)
      do 30 j=1,i
30    q(ixpos(isym)+(j-1)+nprm(isym)*(i-1)) =
     1 sign * q(ixpos(isym)+(i-1)+nprm(isym)*(j-1))
c
c..   exponentiate to form exp(x)
      ns = nprm(isym)
      ns2 = ns**2
      zz = dasum(ns2,q(ixpos(isym)),1)
      fac = zz
      nterm = 2
40    nterm = nterm + 1
      fac = fac*zz/ dfloat(nterm)
      if (fac.gt.1.0d-8) goto 40
      iu = icorr(ns2)
      iv = icorr(ns2)
      call vclr(q(iu),1,ns2)
      nterm=min(nterm,ntexp)
      do 50 i=nterm,2,-1
      call vclr(q(iv),1,ns2)
      call daxpy(ns2,(1.0d0/ dfloat(i)),q(iu),1,q(iv),1)
      call vclr(q(iu),1,ns2)
      call vfill(1.0d0,q(iu),ns+1,ns)
50    call mxmb (q(ixpos(isym)),1,ns, q(iv),1,ns, q(iu),1,ns, ns,ns,ns)
      call fmove (q(ixpos(isym)),q(iv),ns2)
      call vclr(q(ixpos(isym)),1,ns2)
      call vfill(1.0d0,q(ixpos(isym)),nprm(isym)+1,nprm(isym))
      call mxmb (q(iv),1,ns,q(iu),1,ns,q(ixpos(isym)),1,ns,ns,ns,ns)
c      if (lto(4)) then
c      call outsqr(q(ixpos(isym)),nprm(isym),nprm(isym),nprm(isym),'u')
c      end if
      call corlsr(iu)
60    continue
      anorm=dsqrt(anorm)
      if (mcprin) write(iwrite,70)anorm
70    format(/' length of internal x-vector:',f12.8/)
c...  transform orbitals
      do 80 isym=1,nirrr
80    ijpos(isym) = icorr(nsymao(isym)**2)
      ia=icorr(maxbas**2)
      call qget (q(1),ijpos)
      do 100 isym=1,nirrr
      ja=ijpos(isym)+ifreez(isym)*nsymao(isym)
      call fmove (q(ja),q(ia),nsymao(isym)*nsymm(isym))
100   call mxmaa(q(ia),1,nsymao(isym),
     1 q(ixpos(isym)),1,nprm(isym),
     2 q(ja),1,nsymao(isym), nsymao(isym),nprm(isym),nprm(isym) )
      call qput (q(1),ijpos,iwrite)
      call corlsr (ijpos(1))
c
      return
      end
      subroutine updui (q,u,x)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/multic)
INCLUDE(common/mcscra)
      dimension q(*),u(*),x(*)
      call accnt('updui',2)
      ix=icorr(nprim**2)
      iu=icorr(nprim**2)
      iv=icorr(nprim**2)
      ioff=0
      ix1=ix-1
      do 50 isym=1,nirrr
      ns=nprm(isym)
      ns2=ns**2
      call vclr(q(ix),1,ns2)
      i1=ns
      do 20 i=2,ns
      ij=i1
      ji=i
      do 10 j=1,i-1
      ij=ij+1
      ijj=igrd(ij+ioff)
      if(ijj.eq.0) goto 10
      q(ix1+ij)=x(ijj)
      q(ix1+ji)=-x(ijj)
10    ji=ji+ns
20    i1=i1+ns
c
c      call outsqr (q(ix),nprm(isym),nprm(isym),nprm(isym),'x matrix')
c      call outsqr (u(ioff+1),nprm(isym),nprm(isym),nprm(isym)
c     >,'original u matrix')
c
      zz=dasum(ns2,q(ix),1)
      fac=zz
      nterm=3
30    nterm=nterm+1
      fac=fac*zz/ dfloat(nterm)
      if (fac.gt.1.0d-10) goto 30
      call vclr(q(iu),1,ns2)
      do 40 i=nterm,1,-1
      call vclr(q(iv),1,ns2)
      call daxpy(ns2,(1.0d0/ dfloat(i)),q(iu),1,q(iv),1)
      call vclr(q(iu),1,ns2)
      call vfill(1.0d0,q(iu),ns+1,ns)
40    call mxmb (q(ix),1,ns,q(iv),1,ns,q(iu),1,ns,ns,ns,ns)
      call fmove (u(ioff+1),q(iv),ns2)
      call mxmaa(q(iv),1,ns,q(iu),1,ns,u(ioff+1),1,ns,ns,ns,ns)
c      call outsqr (u(ioff+1),nprm(isym),nprm(isym),nprm(isym)
c     >,'resultant u matrix')
50    ioff=ioff+ns2
      call corlsr (ix)
      call accnt(' ',2)
      return
      end
      subroutine orbtra(q,iupos,iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
      dimension q(*),iupos(8),iqpos(8)
      ibase = icorr(0)
      do 10 i=1,nirrr
10    iqpos(i) = icorr(nsymao(i)**2)
      ia=icorr(maxbas**2)
c
      call qget (q(1),iqpos)
c
      do 40 isym=1,nirrr
      iu=iupos(isym)
      ns=nsymao(isym)
      ja=iqpos(isym)+ifreez(isym)*ns
      if (lto(4)) call outsqr (q(iqpos(isym)),ns,ns,ns,'old mos')
      if (lto(4)) call outsqr (q(iu),ns,ns,ns,'transformation matrix')
      call fmove(q(ja),q(ia),nsymm(isym)*ns)
      call mxmaa(q(ia),1,ns, q(iu),1,nsymm(isym),
     1 q(ja),1,ns, ns,nsymm(isym),nsymm(isym))
40    continue
c
c..   orthogonalise the new mos and write to dump
      call qput (q(1),iqpos,iwrite)
      call corlsr (ibase)
      return
      end
      subroutine ofset(iext)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
      common /block / nosymm,nosym(8),iosym(8,maxorb)
INCLUDE(common/mcscra)
c
c.....offsets and dimensions for internal optimization
c
      nmax=0
      numax=0
      do 10 isy=1,nirrr
      nt(isy)=nsymm(isy)
      if(iext.eq.0) nt(isy)=nprm(isy)
10    locc(isy)=0
      ii=0
      do 15 i=1,nprim
      isy=itype(i)
      locc(isy)=locc(isy)+1
      iorb(i)=locc(isy)
      itri(i)=ii
15    ii=ii+i
c      print*,'iorbs'
      do 17 isy=1,nirrr
      jj=0
      do 16 j=1,nprim
      if(itype(j).ne.isy) goto 16
      jj=jj+1
      iorbs(jj,isy)=j-ncoremc
      if(j.le.ncoremc) iorbs(jj,isy)=0
c      print 31,isy,jj,iorbs(jj,isy)
16    continue
17    continue
      nlu=0
      do 30 is=1,nirrr
      ioffu(is)=nlu
      nlu=nlu+locc(is)**2
      nmax=max(nmax,nt(is))
      numax=max(numax,locc(is)**2)
      nld=0
      nlq=0
      do 30 js=1,nirrr
      ks=mults(is,js)
      ioffq(is,js)=nlq
      nlq=nlq+nt(js)*nt(ks)
      nlqq(is)=nlq
      if(ks-js) 20,25,30
20    ioffd(is,js)=nld
      nld=nld+nt(js)*nt(ks)
      goto 30
25    ioffd(is,js)=nld
      nld=nld+nt(js)*(nt(js)+1)/2
30    nldd(is)=nld
      ntdg=nldd(1)
      ntqg=nlqq(1)
      jaa=1
      kaa=1
c
c.....operator adresses
c
c      print*,'operator adresses'
      if(ncoremc.ne.0) then
c.....core-core and valence-core
      do 50 i=1,nprim
      je=min(ncoremc,i)
      do 50 j=1,je
      isy=mults(itype(i),itype(j))
      jadr(i,j)=jaa
      jadr(j,i)=jaa
      kadr(i,j)=kaa
c      print 31,i,j,jaa,kaa
c31   format(1x,2i3,2i6)
      jaa=jaa+nldd(isy)
50    kaa=kaa+nlqq(isy)
      end if
c.....valence-valence ordered to symmetry
      do 45 isy=1,nirrr
      jadrs(isy)=jaa
      kadrs(isy)=kaa
      do 40 i=ncoremc+1,nprim
      is=itype(i)
      js=mults(is,isy)
      do 35 j=ncoremc+1,i
      if(itype(j).ne.js) goto 35
      jadr(i,j)=jaa
      jadr(j,i)=jaa
      kadr(i,j)=kaa
c      print 31,i,j,jaa,kaa
      jaa=jaa+nldd(isy)
      kaa=kaa+nlqq(isy)
35    continue
40    continue
      if(isy.ne.1) goto 45
      ifc=jaa
      jaa=jaa+ntdg
45    continue
      lenj=jaa-1
      lenk=kaa-1
c      print*,'lenj,lenk',lenj,lenk
c      print*,'jadrs',(jadrs(i),i=1,nirrr)
c      print*,'kadrs',(kadrs(i),i=1,nirrr)
      nop=nprim*(nprim+1)/2
c
      ivec=0
      do 60 i=1,nlu
60    igrd(i)=0
c      print*,'igrd'
      do 80 isym=1,nirrr
      ired=0
      ns=nprm(isym)
      do 70 j=ncoremc+1,nprim
      jj=iorb(j)
      ij=ioffu(isym)+(jj-1)*ns
      ji=ioffu(isym)+jj-ns
      do 65 i=1,j-1
      if(itype(i).ne.itype(j)) goto 65
      if(i.gt.ncoremc) then
      ired=ired+1
      if(irottu(ired).eq.0) goto 65
      end if
      if (itype(i).ne.isym) goto 65
      ii=iorb(i)
      if(iosym(isym,ii).ne.iosym(isym,jj)) goto 65
      ivec=ivec+1
      igrd(ij+ii)=ivec
      igrd(ji+ii*ns)=ivec
65    continue
70    continue
      ij=ioffu(isym)
      do 80 i=1,ns
c      print 75,(igrd(ij+j),j=1,ns)
c75   format(1x,19i4)
80    ij=ij+ns
      nrotti=ivec
c
cpaul for single determinant scf calculations
c
      if (nci.eq.1) nrotti = 0
c
      return
      end
      subroutine symden(gam,n)
      implicit REAL  (a-h,o-z)
      dimension gam(*)
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
      ijkl=0
      ij=0
      do 100 i=1,n
      ji=i
      do 90 j=1,n
      ij=ij+1
      kl=0
      do 80 k=1,i
      do 60 l=1,j
      ijkl=ijkl+1
      jikl=ind(ji,kl+l)
      if(jikl.le.ijkl) goto 60
      gam(ijkl)=0.5d0*(gam(ijkl)+gam(jikl))
      gam(jikl)=gam(ijkl)
60    continue
      if(k.eq.i) goto 80
      lk=k+j*n
      do 70 l=j+1,n
      ijkl=ijkl+1
      ijlk=ind(ij,lk)
      if(ijlk.le.ijkl) goto 70
      gam(ijkl)=0.5d0*(gam(ijkl)+gam(ijlk))
      gam(ijlk)=gam(ijkl)
70    lk=lk+n
80    kl=kl+n
      gam(ind(ji,ji))=gam(ijkl)
90    ji=ji+n
100   continue
      return
      end
      subroutine traop(u,aijkl,bijkl,q,nwl,nused,iext,iwrite)
      implicit REAL  (a-h,o-z)
      logical read,start
INCLUDE(common/sizes)
INCLUDE(common/mcff)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
INCLUDE(common/mcscra)
      dimension u(*),q(*),aijkl(*),bijkl(*)
c     character*8 string(2)
c     data string/'j-oper','k-oper'/
      data half/0.5d0/
c
      if (iext.eq.1) call accnt('traop',2)
      if (iext.ne.1) call accnt('traopi',2)
c
c.....offsets
c
      if(iext.ne.0) call ofset(iext)
c
c.....check core available
c
      nused=0
      n2max=2*nmax**2
      nwopt=ntdg+ntqg+max(n2max,nmax+numax)+numax*ntqg
      nwr=nwopt-nwl
      if(nwr.gt.0) write(iwrite,50) nwl,nwopt,nwr
50    format(/
     + 1x,'traop: warning - not all operators will fit in core'/
     + 1x,'space provided:',i8,' : optimal space - ',i8/
     + 1x,'to minimize IO increase work area by ',i9, 'words.')
c
c.....transform ho
c
c      print*,'initial core energy:',core
      core=potnuc+efreez
      naa=1
      if(iext.ne.0) then
      naa=ntdg+1
      call rdfrt(q,ntdg,1,-1)
      call transi(u,q,1,1,q(naa))
      if(ncoremc.gt.0) core=core+ctracem(q)
      call wtfrt(q,ntdg,1,-1)
      else
_IF1()c      call priop(aijkl(ifc),1,1,'h0 init ',0,iwrite)
      call fmove(aijkl(ifc),bijkl(ifc),ntdg)
      call transi(u,bijkl(ifc),1,1,q(naa))
_IF1()c      call priop(bijkl(ifc),1,1,'h0 trans',0,iwrite)
      end if
c
c.....loop over operator symmetries
c
      do 210 isa=1,nirrr
      do 210 jsa=1,isa
      if(locc(isa).eq.0.or.locc(jsa).eq.0) goto 210
      isyop=mults(isa,jsa)
c
c.....transformation of coulomb operators first
c
      nlg=nldd(isyop)
      nlr=nlg
      ifil=1
      np=1
c
70    iva=naa+nlg
      nbf=iva+numax
      j00=max(nbf+nmax,iva+n2max)
      iva=iva-1
      ia=1
      ja=1
      start=.true.
      read=.true.
80    iop=0
c
c.....loop over transformed operators
c
      do 190 k=1,nprim
      do 190 l=1,k
      iop=iop+1
      ks=max(itype(k),itype(l))
      ls=min(itype(k),itype(l))
      if(ks.ne.isa.or.ls.ne.jsa) goto 190
      if(start) call vclr(q(naa),1,nlg)
      if(.not.start) call rdfrt(q(naa),nlg,ifil,iadw+iop)
      ka=ioffu(itype(k))+(iorb(k)-1)*locc(itype(k))
      la=ioffu(itype(l))+(iorb(l)-1)*locc(itype(l))
      iop1=ia*(ia-1)/2+ja-1
      jaa=j00
      nvec=0
      j1=ja
c
c.....loop over original operators
c
      do 130 i=ia,nprim
      do 120 j=j1,i
      iop1=iop1+1
      is=max(itype(i),itype(j))
      js=min(itype(i),itype(j))
      if(is.ne.isa.or.js.ne.jsa) goto 120
      if(.not.read) goto 90
      if(jaa+nlg.gt.nwl+1) goto 230
      if(iext.eq.0) call fmove(aijkl(jadr(i,j)),q(jaa),nlg)
      if(iext.eq.1) call rdfrt(q(jaa),nlg,ifil,iadr+iop1)
_IF1()c      call priop(q(jaa),isyop,np,string(ifil),10*i+j,iwrite)
      if(np.ne.0) goto 90
      if(itype(i).ne.isa) call transp(q(jaa),q(jaa),isyop,q(nbf))
      if(isa.ne.jsa.or.i.eq.j) goto 90
      if(jaa+2*nlg.gt.nwl+1) goto 230
      call transp(q(jaa),q(jaa+nlg),1,q(nbf))
90    ii=iorb(i)
      jj=iorb(j)
      fak=0.0d0
      if(itype(i).eq.itype(k)) fak=u(ka+ii)*u(la+jj)
      if(np.ne.0.or.isa.ne.jsa) goto 100
      jaa=jaa+nlg
      nused=max(nused,jaa-1)
      nvec=nvec+1
      q(iva+nvec)=fak
      fak=0.0d0
      if(i.eq.j) goto 110
100   if(itype(j).eq.itype(k)) fak=fak+u(ka+jj)*u(la+ii)
      if(i.eq.j) fak=fak*half
      jaa=jaa+nlg
      nused=max(nused,jaa-1)
      nvec=nvec+1
      q(iva+nvec)=fak
110   if(jaa+nlr.gt.nwl) goto 140
120   continue
      j1=1
130   continue
c
      je=nprim
      goto 150
140   ie=i
      je=j
150   call mxmb(q(j00),1,nlg,q(iva+1),1,0,q(naa),1,0,nlg,nvec,1)
      if(je.ne.nprim) goto 160
c.....final transformation u(dagger)*op*u
      call transi(u,q(naa),isyop,np,q(iva+1))
      if(np.eq.0.and.itype(k).ne.isa) call transp(q(naa),q(naa),isyop,
     1  q(nbf))
160   read=.false.
      if(iext.eq.0) call fmove(q(naa),bijkl(jadr(k,l)),nlg)
      if(iext.eq.1) call wtfrt(q(naa),nlg,ifil,iadw+iop)
_IF1()c      call priop(q(naa),isyop,np,string(ifil),10*k+l,iwrite)
      if(je.ne.nprim.or.k.ne.l) goto 190
c.....contributions to core fock operator
      if(k.gt.ncoremc.or.iext.eq.0) goto 190
      fak=2.0d0
      if(np.eq.0) then
      call redsub(q(naa),q)
      else
      call daxpy(ntdg,2.0d0,q(naa),1,q,1)
      end if
190   continue
      if(je.eq.nprim) goto 200
      if(iext.eq.0) goto 230
c.....next batch of original operators
      start=.false.
      read=.true.
      maxv=locc(isa)*locc(jsa)
      if(isa.eq.jsa.and.ifil.eq.1) maxv=locc(isa)*(locc(isa)+1)/2
      nwopt=j00-1+maxv*nlg
      ia=ie
      ja=je+1
      if(ja.le.ie) goto 80
      ja=1
      ia=ia+1
      goto 80
c
c.....repeat for exchange operators (iext=1)
c
200   if(ifil.eq.2.or.iext.eq.0) goto 210
      ifil=2
      np=0
      nlg=nlqq(isyop)
      nlr=nlg
      if(isyop.eq.1) nlr=2*nlg
      goto 70
210   continue
      if(iext.ne.0) goto 220
c
c.....replace old coulomb operators and generate exchange operators
c
      call dcopy(lenj,bijkl,1,aijkl,1)
      call sortjk(aijkl,bijkl)
      if(ncoremc.gt.0) then
c
c.....generate fc
c
      core=core+ctracem(aijkl(ifc))
      do 215 i=1,ncoremc
      call daxpy(ntdg
     + ,2.0d0,aijkl(jadr(i,i)),1,aijkl(ifc),1)
215   call redsub(bijkl(kadr(i,i)),aijkl(ifc))
      core=core+ctracem(aijkl(ifc))
_IF1()c      call priop(aijkl(ifc),1,1,'fc trans',0,iwrite)
c      print*,'new core energy',core
      end if
      call accnt(' ',2)
      return
220   continue
      if(ncoremc.gt.0) core=core+ctracem(q)
      call wtfrt(q,ntdg,1,0)
      call accnt(' ',2)
      return
230   write(iwrite,240)
240   format(/1x,'not enough core in traop ')
      call caserr('insufficient memory in traop')
      return
      end
      subroutine sortjk(opj,opk)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/multic)
INCLUDE(common/mcscra)
      dimension opj(*),opk(*)
c
c.....routine generators internal exchange operators from coul operators
c
c.....loop over exchange operators k(ij)
c
      do 190 i=1,nprim
      is=itype(i)
c     ii=iorb(i)
      do 180 j=1,i
      js=itype(j)
      isyk=mults(is,js)
      nj=nprm(js)
      jj=iorb(j)
      j2=2*jj
c
c.....loop over symmetry blocks of k(ij)
c
      do 150 ks=1,nirrr
      ls=mults(isyk,ks)
      isyj=mults(is,ks)
      nl=nprm(ls)
      nk=nprm(ks)
      ka=kadr(i,j)+ioffq(isyk,ks)-nk
      if(ls-js) 50,80,120
50    j1=ioffd(isyj,js)+jj-nj-1
      do 70 k=1,nprim
      if(itype(k).ne.ks) goto 70
      ja=j1+jadr(k,i)
      do 60 l=1,nl
60    opk(ka+l*nk)=opj(ja+l*nj)
      ka=ka+1
70    continue
      goto 150
80    j1=ioffd(1,js)+itri(jj)-1
      do 110 k=1,nprim
      if(itype(k).ne.ks) goto 110
      ja=j1+jadr(k,i)
      do 90 l=1,jj
90    opk(ka+l*nk)=opj(ja+l)
      ja=ja+j2
      do 100 l=jj+1,nl
      opk(ka+l*nk)=opj(ja)
100   ja=ja+l
      ka=ka+1
110   continue
      goto 150
120   j1=ioffd(isyj,ls)+(jj-1)*nl-1
      do 140 k=1,nprim
      if(itype(k).ne.ks) goto 140
      ja=j1+jadr(k,i)
      do 130 l=1,nl
130   opk(ka+l*nk)=opj(ja+l)
      ka=ka+1
140   continue
150   continue
_IF1()c      call priop(opj(jadr(i,j)),isyk,1,'j-oper  ',10*i+j,iwrite)
_IF1()c      call priop(opk(kadr(i,j)),isyk,0,'k-oper  ',10*i+j,iwrite)
180   continue
190   continue
      return
      end
      subroutine transi(u,op,isy,np,q)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/multic)
INCLUDE(common/mcscra)
      dimension u(*),op(*),q(*)
c
c.....transforms internal part of operators
c
c.....isy=symmetry of operator
c.....np=0 square operator (k-op)
c.....np=1 symmetric operator (j-op, lower triangle stored)
c
c.....q: work array. for np=0 length of q = largest symm. block of u
c.....               for np=1 length of q = 2 * largest sym. block of u.
c
c
      ja=1
      do 100 isa=1,nirrr
      isb=mults(isy,isa)
      if(isb.gt.isa.and.np.ne.0) goto 100
      na=nt(isa)
      nb=nt(isb)
      if(na.eq.0.or.nb.eq.0) goto 100
      la=locc(isa)
      lb=locc(isb)
      ia=ioffu(isa)+1
      ib=ioffu(isb)+1
      if(isa.eq.isb.and.np.eq.1) goto 50
      if(lb.eq.0) goto 10
      call mxmaa(op(ja),1,na,u(ib),1,lb,q,1,na,na,lb,lb)
10    if(nb.gt.lb) call fmove(op(ja+lb*na),q(lb*na+1),na*(nb-lb))
      if(la.gt.0)callmxmaa(q,na,1,u(ia),1,la,op(ja),na,1,nb,la,la)
      iaa=0
      jaa=ja-1
      if(la.ge.na) goto 40
      do 30 j=1,nb
      do 20 i=la+1,na
20    op(jaa+i)=q(iaa+i)
      jaa=jaa+na
30    iaa=iaa+na
40    ja=ja+na*nb
      goto 100
50    if(la.eq.0) goto 90
      nbf=na**2+1
      call square(q(nbf),op(ja),na,na)
      call mxmaa(q(nbf),1,na,u(ib),1,lb,q,1,na,na,lb,lb)
      if(nb.gt.lb) call fmove(q(nbf+lb*na),q(lb*na+1),na*(nb-lb))
      call mxmaa(q,na,1,u(ia),1,la,q(nbf),na,1,nb,la,la)
      if(la.ge.na) goto 80
      iaa=0
      jaa=nbf-1
      do 70 j=1,nb
      do 60 i=la+1,na
60    q(jaa+i)=q(iaa+i)
      iaa=iaa+na
70    jaa=jaa+na
80    call reduce(q(nbf),op(ja),na)
90    ja=ja+na*(na+1)/2
100   continue
      return
      end
      subroutine amat(fop,i,j,ieqj,ai,aj,ni,nj)
      implicit REAL  (a-h,o-z)
      logical ieqj
      dimension fop(ni,nj),ai(ni,ni),aj(nj,nj)
c
c.....routine calculates a(k,i)=<k/f(ij)/j>
c
_IF(cray,ibm,vax)
      do 10 k=1,ni
   10 ai(k,i)=ai(k,i)+fop(k,j)
_ELSE
      call vadd(ai(1,i),1,fop(1,j),1,ai(1,i),1,ni)
_ENDIF
      if(ieqj) return
      do 20 k=1,nj
20    aj(k,j)=aj(k,j)+fop(i,k)
      return
      end
      subroutine goper(q,opj,opk,fop,gop,dijkl,dikjl,diljk,gam,grad,
     * hes)
      implicit REAL  (a-h,o-z)
      logical ieqj,isy1
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/mcscra)
      dimension opj(*),opk(*),gop(*),fop(*),dijkl(*),dikjl(*),diljk(*)
      dimension q(*),iapos(8),gam(*),grad(*),hes(*)
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
c.....routine calculates f and g operators
c
      call accnt('goper',2)
      energy=core
      ifc1=ifc-1
      do 5 isy=1,nirrr
5     iapos(isy)=icorr(nprm(isy)**2)
      call vclr(hes,1,nrotti**2)
      call vclr(q(iapos(1)),1,nlu)
      if(ncoremc.gt.0) then
c
c...  form gc
c
      igc=icorr(ntqg)
      call dcopy(ntdg,opj(ifc),1,fop,1)
      do 15 i=1,nact
      ii=i+ncoremc
      ij=ilifa(i)+ic1d
      is=itypea(i)
      do 10 j=1,i
      if(itypea(j).ne.is) goto 10
      if(gam(ij+j).eq.0) goto 10
      fac=gam(ij+j)
      if(i.ne.j) fac=2.0d0*fac
      jj=j+ncoremc
      call daxpy(ntdg,fac,opj(jadr(ii,jj)),1,fop,1)
      call redadd(opk(kadr(ii,jj)),fop,-0.5d0*gam(ij+j),i-j)
10    continue
15    continue
c
c.....core-valence part of a-matrix
c
      do 25 is=1,nirrr
      if(ncor(is).eq.0) goto 25
      nc=ncor(is)
      ni=nprm(is)
      ig=ioffq(1,is)+igc
      id=ioffd(1,is)+1
      call square(q(ig),fop(id),ni,ni)
      call fmove(q(ig),q(iapos(is)),ni*nc)
25    continue
      end if
c
      do 180 i=1,nprim
      do 180 j=1,i
      is=itype(i)
      js=itype(j)
      ii=iorb(i)
      jj=iorb(j)
      ieqj=i.eq.j
      ni=nprm(is)
      nj=nprm(js)
      isyg=mults(is,js)
      isy1=isyg.eq.1
      leng=ni*nj
      lenf=leng
      ijj=(jj-1)*ni+ii
      if(isy1) then
      ijj=ind(ii,jj)
      ifcij=ifc1+ioffd(1,is)+ijj
      lenf=ni*(ni+1)/2
      end if
      fac=2.0d0
      if(ieqj) fac=1.0d0
      if(j.gt.ncoremc) goto 50
      if(i.gt.ncoremc) goto 30
c
c.....core-core
c
      call loper(opj(jadr(i,j)),opk(kadr(i,j)),gop,isyg,is,js,ni,nj)
      if(ieqj) call vadd(gop,1,q(igc+ioffq(1,is)),1,gop,1,ni*ni)
      goto 170
c
c.....core-valence
c
30    il=ilifa(i-ncoremc)+ic1d
      call vclr(gop,1,ni*nj)
      do 40 ll=1,nact
      if(itypea(ll).ne.is) goto 40
      l=ncoremc+ll
      if(gam(il+ll).eq.0) goto 40
      call loper(opj(jadr(l,j)),opk(kadr(l,j)),fop,isyg,is,js,ni,nj)
      call daxpy(ni*nj,0.5d0*gam(il+ll),fop,1,gop,1)
40    continue
      goto 170
c
c.....valence-valence
c
50    iii=i-ncoremc
      jjj=j-ncoremc
      il=ilifa(iii)
      jl=ilifa(jjj)
      ij=il+jjj
      iv=0
c
c.....set up density matrix blocks required
c
      do 70 k=1,nact
      ks=itypea(k)
      ls=mults(isyg,ks)
      if(nactt(ls).eq.0) goto 70
      kl=ilifa(k)
      jk=jl+k
      ik=il+k
      do 60 l=1,k
      if(itypea(l).ne.ls) goto 60
      iv=iv+1
      dijkl(iv)=gam(ind(ij,kl+l))
      dikjl(iv)=gam(ind(ik,jl+l))
      diljk(iv)=gam(ind(il+l,jk))
60    continue
      if(.not.isy1) goto 70
c.....corrections for k=l
      dijkl(iv)=0.5d0*dijkl(iv)
      diljk(iv)=0.0d0
70    continue
      iv1=iv
      if(isy1) then
c.....contribution of core fock operator
      iv1=iv+1
      dijkl(iv1)=0.5d0*gam(ic1d+ij)
      energy=energy+fac*dijkl(iv1)*opj(ifcij)
      end if
c
c.....operators f(ij) * 0.5
c
      iblock=max(is,js)
      ja=jadrs(isyg)+ioffd(isyg,iblock)
      call mxmaa(opj(ja),1,nldd(isyg),dijkl,1,0,fop,1,0,lenf,iv1,1)
c
c.....contributions of f(ij) to a-matrix
c
      if(is-js) 100,80,90
80    call square(gop,fop,ni,ni)
      call amat(gop,ii,jj,ieqj,q(iapos(is)),q(iapos(js)),ni,nj)
      call vclr(fop,1,leng)
      goto 110
90    call amat(fop,ii,jj,ieqj,q(iapos(is)),q(iapos(js)),ni,nj)
      call fmove(fop,gop,leng)
      call vclr(fop,1,leng)
      goto 110
100   call amat(fop,jj,ii,ieqj,q(iapos(js)),q(iapos(is)),nj,ni)
      call vclr(gop,1,leng)
c
c.....contributions of exchange operators k(kl) to g(ij) * 0.5
c
110   ka=kadrs(isyg)+ioffq(isyg,is)
      call mxmb(opk(ka),1,nlqq(isyg),dikjl,1,0,gop,1,0,leng,iv,1)
c
c.....contributions of exchange operators k(kl) to g(ji) * 0.5
c
      ka=kadrs(isyg)+ioffq(isyg,js)
      call mxmb(opk(ka),1,nlqq(isyg),diljk,1,0,fop,1,0,leng,iv,1)
      call tradd(gop,fop,ni,nj)
c
c.....contribution of present g to hessian
c
170   call hesi(gop,hes,ii,jj,is,js,ni,nj,ieqj)
c
180   continue
c
c.....contributions of a matrix to hessian and gradient
c
      do 200 is=1,nirrr
      ii=iapos(is)+(nprm(is)+1)*ncor(is)
      energy=energy+dsum(nactt(is),q(ii),nprm(is)+1)
200   continue
      call gradi(q(iapos(1)),grad,hes)
      call accnt(' ',2)
      call corlsr(iapos(1))
      return
      end
      subroutine loper(opj,opk,opl,isy,is,js,ni,nj)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mcscra)
      dimension opj(*),opk(*),opl(*)
c
c.....forms operators l(ij)= 4*k(ij) - k(ji) - j(ij)
c
      if(is-js) 10,20,30
10    ja=ioffd(isy,js)+1
      call trnsps(opj(ja),opl,nj,ni)
      goto 40
20    ja=ioffd(isy,is)+1
      call square(opl,opj(ja),ni,ni)
      goto 40
30    ja=ioffd(isy,is)+1
      call fmove(opj(ja),opl,ni*nj)
40    ijq=ioffq(isy,is)
      jiq=ioffq(isy,js)-nj+1
      ij=0
      do 60 j=1,nj
      do 50 i=1,ni
50    opl(ij+i)=-opl(ij+i)+4.0d0*opk(ijq+i)-opk(jiq+i*nj)
      ij=ij+ni
      ijq=ijq+ni
60    jiq=jiq+1
      return
      end
      subroutine gradi(a,grad,h)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/mcscra)
      dimension grad(*),a(*),h(nrotti,nrotti)
c
c.....calculate gradient and symmetrize a
c
_IF1()c      call priop(a,1,0,'amat',0,iwrite)
      do 200 is=1,nirrr
      ni=locc(is)
      j1=ioffu(is)
      i1=j1+ni
      do 120 i=2,ni
      ji=j1+i
      ij=i1
      do 110 j=1,i-1
      ij=ij+1
      ijj=igrd(ij)
      if(ijj.ne.0) grad(ijj)=a(ij)-a(ji)
      a(ij)=0.5d0*(a(ij)+a(ji))
      a(ji)=a(ij)
110   ji=ji+ni
120   i1=i1+ni
200   continue
c
c.....contributions of a+a(dagger) to hessian matrix
c
      do 100 is=1,nirrr
      ni=locc(is)
      i1=ioffu(is)+ni
      do 90 ii=2,ni
      j1=ioffu(is)
      do 80 jj=1,ii-1
      ijj=igrd(i1+jj)
      if(ijj.eq.0) goto 80
      ik=i1
      jk=j1
      do 70 kk=1,ni
      ik=ik+1
      jk=jk+1
      ikk=igrd(ik)
      if(ikk.eq.0) goto 30
      if(ii-kk) 10,30,20
10    h(ikk,ijj)=h(ikk,ijj)+a(jk)
      goto 30
20    h(ikk,ijj)=h(ikk,ijj)-a(jk)
30    jkk=igrd(jk)
      if(jkk.eq.0) goto 70
      if(jj-kk) 50,70,40
40    h(jkk,ijj)=h(jkk,ijj)+a(ik)
      goto 70
50    h(jkk,ijj)=h(jkk,ijj)-a(ik)
70    continue
80    j1=j1+ni
90    i1=i1+ni
100   continue
      return
      end
      subroutine hesi(gij,h,i,j,is,js,ni,nj,ieqj)
      implicit REAL  (a-h,o-z)
      logical ieqj,neg
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/mcscra)
      dimension gij(ni,nj),h(nrotti,nrotti)
c
c.....contributions of operator g(ij) to hessian matrix
c
      ik=ioffu(is)+i
      j1=ioffu(js)+j
      do 90 k=1,ni
      ikk=igrd(ik)
      if(ikk.eq.0) goto 90
      jl=j1
      neg=k.gt.i
      do 70 l=1,nj
      if(l.eq.j) goto 60
      jll=igrd(jl)
      if(jll.eq.0) goto 70
      if(neg) goto 50
      h(ikk,jll)=h(ikk,jll)+gij(k,l)
      if(.not.ieqj) h(jll,ikk)=h(jll,ikk)+gij(k,l)
      goto 70
50    h(ikk,jll)=h(ikk,jll)-gij(k,l)
      if(.not.ieqj) h(jll,ikk)=h(jll,ikk)-gij(k,l)
      goto 70
60    neg=.not.neg
70    jl=jl+nj
90    ik=ik+ni
      return
      end
      subroutine diismc(q,r,g,wgt,b,n,nwgt,ndim,ndel,bfak,ndmax,iblk1)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (lseg=511)
INCLUDE(common/jobopt)
      common /mcopt / var,varc,thzr,one
INCLUDE(common/multic)
      dimension q(*),r(*),g(*),wgt(*),b(*)
c
c      call accnt('diis',2)
      ibuf=icorr(lseg)
      n1=n+1
      do 5 i=1,nwgt
  5   g(i)=g(i)/wgt(i)
      one=1.0d0
      singlr=1.0d-5
      ifil=numscr
      nblk=lensec(n)
      ic=ndim*(ndim+1)/2
      ndim1=ndim
      if(ndim.eq.0) then
      ndel=0
      bfak=1.0d0
      goto 20
      end if
      iblk=iblk1+2*ndel*nblk
      do 15 i=1,ndim
      ia=1
      b(ic+i)=0
      do 10 ib=1,nblk
      l=min(n1-ia,lseg)
      call rdedx(q(ibuf),l,iblk+ib,ifil)
      b(ic+i)=b(ic+i)+bfak*ddot(l,g(ia),1,q(ibuf),1)
10    ia=ia+lseg
15    iblk=iblk+2*nblk
      call rdedx(b,ic,iblk+1,ifil)
20    ndim=ndim+1
      b(ic+ndim)=bfak*ddot(n,g,1,g,1)
      fak=one/b(ic+ndim)
      call dscal(ic+ndim,fak,b,1)
      bfak=bfak*fak
      iblk=iblk1+2*(ndim1+ndel)*nblk+1
      call search(iblk,ifil)
      call wrt3s(g,n,ifil)
      call wrt3s(r,n,ifil)
      if(ndim.eq.1) goto 110
      iv=ic+ndim
      ix=iv+ndim
      iv2=ix+ndim**2+1
c     lmax=iv2+ndim-1
c     m = 1
      if (ndim.gt.ndmax) goto 60
       call square(b(ix+1),b,ndim,ndim)
c      call outtri(b,ndim,'b-matrix')
      ifail = 0
      call f02abf(b(ix+1),ndim,ndim,b(iv+1),b(ix+1),ndim,b(iv2),ifail)
c      call outvec(b(iv+1),ndim,'eigenvalues of b-matrix')
      do 40 i=1,ndim
40    if(b(iv+i).gt.singlr) goto 50
50    if(i.eq.1) goto 110
c      m=i-1
c60    ndel=ndel+m
c      ijo=m*(m+3)/2
c      ijn=0
c      ndim=ndim-m
c      if(ndim.gt.1) goto 80
c      print 70, (b(iv+i),i=1,ndim+m)
c70    format(/'warning: b-matrix in diis badly conditioned.'//
c     1  ' eigenvalues:',  (t15,9e12.4))
c...  abandon diis & set ndim to zero
60    ndim = 0
      ndel = 0
      bfak = one
      goto 160
c80    do 100 i=1,ndim
c      do 90 j=1,i
c90    b(ijn+j)=b(ijo+j)
c      ijn=ijn+i
c100   ijo=ijo+i+m
c      if (ndim.gt.1) goto 30
110   call wrt3s(b,ndim*(ndim+1)/2,ifil)
      if(ndim.eq.1) goto 160
      ij=ix
      do 130 i=1,ndim
      ci=0
      do 120 j=1,ndim
120   ci=ci+b(ij+j)
      b(iv+i)=ci/b(iv+i)
130   ij=ij+ndim
      call mxmaa(b(ix+1),1,ndim,b(iv+1),1,0,b,1,0,ndim,ndim,1)
      cc=one/dsum(ndim,b,1)
      call dscal(ndim,cc,b,1)
c      call outvec(b,ndim,'scaling coefficients')
      call dscal(n,b(ndim),r,1)
      call dscal(n,b(ndim),g,1)
      iblk=iblk1+2*ndel*nblk
      ndim1 = ndim-1
      do 150 i=1,ndim1
      ia=1
      do 140 ib=1,nblk
      l=min(n1-ia,lseg)
      iblk=iblk+1
      call rdedx(q(ibuf),l,iblk,ifil)
      call daxpy(l,b(i),q(ibuf),1,g(ia),1)
140   ia=ia+lseg
      ia=1
      do 150 ib=1,nblk
      l=min(n1-ia,lseg)
      iblk=iblk+1
      call rdedx(q(ibuf),l,iblk,ifil)
      call daxpy(l,b(i),q(ibuf),1,r(ia),1)
150   ia=ia+lseg
160   do 170 i=1,nwgt
170   g(i)=g(i)*wgt(i)
      call corlsr(ibuf)
c      call accnt(' ',2)
      return
      end
c ******************************************************
c ******************************************************
c             =   mcwvfn     =
c ******************************************************
c ******************************************************
      subroutine wvfn(q,iq,iwrite,ipunch)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      logical btest
INCLUDE(common/multic)
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
      common/coruse/ ncruse
      common /mccore/ intrel,lword,ltop,lmax,lmin,mreal
      common/ cpos / iupos(8),iapos(8)
      integer ixpos(8)
      character*72 outbuf
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
INCLUDE(common/mcff)
INCLUDE(common/mcaddr)
      common /restar/ nprint(6),iresga
      dimension q(*),iq(*)
c
      ibase=icorr(0)
      ncruse=0
      iresga=0
      write (outbuf,10)
10    format('mcscf convergence:')
      call wto (outbuf,72)
      write (outbuf,20)
20    format('it.',10x,'energy',19x,'energy change',5x,'grad',6x,
     1'step')
      call wto (outbuf,72)
      cpulas = cpulft(1)
c
      isignh = 1
      iretrn=1
30    cpu = cpulft(1)
      iterp = iter+1
      if (iterp.le.maxcyc) goto 50
      write(iwrite,40)
 40   format(/1x,51('*')/
     *        ' ** mcscf ****  maximum number of iterations reached'/
     *1x,51('*')/)
      iretrn=0
c      iretrn=2
      iresga=3
      write(iwrite,185)lmin,ncruse,lmax
      return
c
50    if ( (cpu-cpulas)*safty(1)+safty(2) .lt. cpulft(0) ) goto 70
      write(iwrite,60)
 60   format(/1x,48('*')/
     *' ** mcscf ****  insufficient cpu time to continue'/1x,48('*')/)
      iretrn=0
c      iretrn=1
      if (mcprin) write(iwrite,185)lmin,ncruse,lmax
      iresga=3
      return
c
70    cpulas = cpu
      if (mcprin) then
       dumt = seccpu()
       write(iwrite,80) iterp,dumt
      endif
 80   format(/1x,79('*')/
     *1x,'** mcscf ****  start of iteration',i3,' at time',f9.2/)
c
c======= i n t e g r a l    t r a n s f o r m a t i o n ================
c
      iblf=1
      ifinit=0
      call mctran(q(1),iq(1))
c ====== g r a d i e n t  &  h e s s i a n   i n i t i a l i s a t i o n
c
      isigma = 1
      if (btest(itinfo(iterp),2)) isigma = 2
      n = nrot+nstate*nci
      if (isigma.eq.2) n = nrot
      iaugmx = 0
c
c...  optimise on internal space
      maxaug=min(maxaug,nrot-1)
      if(augvar.le.0) maxaug=0
      iblkp=2*maxaug*lensec(nrot+1)
     1     +2*maxdis*lensec(nstate*nci+nrot)
     >     +lensec(maxdis*(maxdis+1)/2)+iblk8
      idgot=0
      if (btest(itinfo(iterp),7).or.btest(itinfo(iterp),1))
     >     call optint(q(1),iq(1),idgot,iwrite,ipunch)
      if (iter.eq.0) idgot=0
c
      iwrnr=0
      if (btest(itinfo(iterp),3)) iwrnr=1
      call hesini(q(1),iwrite)
c
      write (outbuf,90) energy,gradnt
      call wto(outbuf,72)
      write(iwrite,90)energy,gradnt
90    format(' energy =',f20.10,10x,'gradient =',e11.2)
      if (gradnt.lt.conv .or.  dabs(energy-elast).lt.econv) iretrn=0
      if(energy-elast.gt.thrdiv) then
c
c....  no convergence
c
      write(iwrite,91)
 91   format(/1x,22('*')/
     *' *** no convergence ***'/
     *1x,22('*')/)
      iretrn=0
      iresga=3
      return
      end if
c
c =========== o p t i m i s a t i o n ==================================
c....  if convergence, do method for iteration 40
      if (iretrn .eq. 0) itinfo (iterp) = itinfo(40)
c
c...  decide on method
      method = 5
      do 100 i=3,6
      if (btest(itinfo(iterp),i)) method = i
100   continue
      isol=icorr(n+1)
      goto (110,140,150,160), method-2
c
c...  werner-meyer non-linear
110   continue
      call corlsr(isol)
      do 115 isym=1,nirrr
      iupos(isym)=icorr(nsymm(isym)**2)
      call vclr(q(iupos(isym)),1,nsymm(isym)**2)
      call vfill(1.0d0,q(iupos(isym)),nsymm(isym)+1,nsymm(isym))
 115  continue
      nittra=0
c     anormo=1.d5
      isol = icorr(n+1)
      icivec=isol+nrot
120   call optwrn (q(1),iq(1),q(isol),n,nittra,idgot,iwrite)
      stepl=dnrm2(nrot,q(isol),1)
      if(iterp.lt.nitrep) goto 135
      if(nittra.ge.maxrep) goto 135
      nittra=nittra+1
      call corlsr(icivec)
      do 130 isym=1,nirrr
130   ixpos(isym)=icorr(nprm(isym)**2)
      call intu(q(1),ixpos,q(isol),anorm,iwrite)
      if(anorm.lt.1.d-5) goto 135
      ixx=iadr
      iadr=iadw
      iadw=ixx
      if(ifinit.eq.0) call finit
      if(iadw.eq.0) call freset
      if(iadw.ne.0) call offset
      l=icorrm()
      iad=icorr(-l)
      call traop(q(ixpos(1)),q(iad),q(iad),q(iad),l,luse,1,iwrite)
      ncruse=max(ncruse,iad+luse-1)
      call corlsr(icivec)
      if(isigma.ne.2) icivec=icorr(nci*nstate)
      goto 120
135   call orbtra(q(1),iupos,iwrite)
      goto 175
c
c...  augmented hessian
140   call optaug (q(1),q(isol),n,iwrite)
      goto 170
c
c..   straight nr
150   call optnr (q(1),q(isol),n,iwrite)
      goto 170
c
c...  no optimisation
160   continue
      call vclr(q(isol),1,n)
      enext = energy
      goto 170
c
c ======================================================================
c
170   continue
      if (btest(iprint,11)) call outvec (q(isol),n,'computed step')
      stepl = dnrm2(n,q(isol),1)
c
      call mcupda (q(1),q(isol),iwrite)
c
175   if (mcprin) write(iwrite,180)stepl
180   format(/' step length =',e11.2)
      call mcchek(q(1),iwrite)
      if (stepl.lt.sconv .or.  dabs(enext-energy).lt.econv) iretrn = 0
      call corlsr(ibase)
      if(iretrn.eq.0.and.mcprin) write(iwrite,185)lmin,ncruse,lmax
185   format(/1x,'core required: ',i8,' (min)'/
     +        1x,'               ',i8,' (maximum used)'/
     +        1x,'             ',i10,' (maximum allocated)')
c
      call mcprnt(q(1),iprint,iwrite,ipunch)
c
c...  iteration completed
      if (mcprin) then
       write(iwrite,190)
190    format(/' convergence progress'
     1        /' --------------------'
     2       //' iteration',17x,'gradient',23x,'energy',20x
     3       ,'change in energy')
       if (iter.ge.2) write(iwrite,200)iter-1,glast2,elast2
       if (iter.ge.2) write(iwrite,200)iter  ,glast ,elast ,elast-elast2
       if (iter.eq.1) write(iwrite,200)iter  ,glast ,elast
       if (iter.ge.1) write(iwrite,200)iterp ,gradnt,energy,energy-elast
     1  ,'variational'
       if (iter.eq.0) write(iwrite,200)iterp ,gradnt,energy
       write(iwrite,200)iterp+1,0.0d0,enext ,enext-energy,'(projected)'
200    format(i6,4x,2f30.12,:,f30.12,:,a15)
      endif
      write (outbuf,210) iterp,energy,energy-elast,gradnt,stepl
210   format(i2,2f25.12,2e10.2)
      call wto (outbuf,72)
c
      iter   = iterp
      glast2 = glast
      glast  = gradnt
      elast2 = elast
      elast  = energy
      slast  = stepl
      call mcdump
      if (iretrn.ne.0) call freset
      if (iretrn.ne.0) goto 30
c...  convergence reached
      ifwvfn = 1
      write(iwrite,220)
 220  format(/1x,34('*')/
     *' ** mcscf ****  convergence reached'/1x,34('*')/)
      iresga=0
      return
      end
      subroutine aughes (q,sol,grad,vec,ddr,itera,iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (radtol=0.001d0, cnvfc2=0.01d0,mv=30,acc=0.05d0)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
INCLUDE(common/mcaddr)
      common /scra  / zm(mv,mv),zn(mv+1,mv+1),p(mv),beta(mv+1,mv+1),
     1              gamma(mv+1),scrach(mv+1),r(mv+1),eig(mv+1),wks(mv+1)
      dimension q(*),sol(*),grad(*),vec(*)
c
      radius = 0.5d0
c     nc=1
      naug = nrot + 1
      n = nrot
      nblk=lensec(nrot+1)
      iblkr=iblk8
      iblkg=nblk+iblk8
      nblk=nblk*2
      ifil = numscr
      iblkif=iblk8
      nalph = 0
      alpha = 1.0d0
      nvec = 0
      eig(1) = 0.0d0
      beta0 = 1.0d0
      dl = 0.0d0
      itera = 0
      nconv = 0
      test = 1.0d0
c
      sol1=sol(naug)
      grad(naug)=0.0d0
      sol(naug)=0.0d0
      call dcopy(naug,grad,1,q(igrad),1)
      goto 120
c
c...  start of davidson iteration
40    itera = itera + 1
      nalph = 0
      dalpha = 1.0d0
c     dll = 1.0d30
c...  set up davidson matrix for a given value of alpha
50    do 60 i=1,nvec
      do 60 j=1,nvec
60    zn(i,j) = zm(i,j) / alpha
      do 70 i=1,nvec
      zn(i,nvec+1) = p(i)
70    zn(nvec+1,i) = p(i)
      zn(nvec+1,nvec+1) = 0.0d0
      ifail=0
      call f02abf (zn,mv+1,nvec+1,eig,beta,mv+1,scrach,ifail)
c      call outsqr (beta,mv+1,nvec+1,nvec+1,'beta')
c...  find the norm of the step vector
      beta0 = beta(nvec+1,1)
      dl =  dabs( 1.0d0/(alpha*beta0)) *
     1     dnrm2(nvec,beta,1)
c      print*,'alpha,dl ',alpha,dl
      if ((alpha.eq.1.0d0.and.dl.le.radius)
     1 .or. ( dabs((dl-radius)/radius).le.radtol)
     2 .or. nalph.gt.20) goto 120
c..   this shift not ok, so we must update it
c     dll = dl
      nalph = nalph + 1
      do 110 i=1,nvec
      do 100 j=1,nvec
100   zn(j,i) = zm(j,i) - (1.0d0/beta0)*beta(j,1)*p(i)
      zn(i,i) = zn(i,i) - eig(1)*alpha
110   r(i) = 2.0d0*eig(1)*beta(i,1)
c      call outsqr (zn,mv+1,nvec,nvec,'n')
      ifail = 1
      call f04arf (zn,mv+1,r,nvec,gamma,wks,ifail)
      if (ifail.eq.1) call caserr('a not pos def in updating alpha')
      if (ifail.eq.2)write(iwrite,1000)
1000  format(               ' ill conditioning in updating alpha')
c      call outvec (gamma,nvec,'gamma')
      dalpha = ((1-beta0**2)/ddot(nvec,gamma,1,beta,1))
     1              * (1.0d0- dl/radius)
      if (alpha+dalpha.lt.1.0d0) dalpha = 0.5d0-0.5d0*alpha
      alpha = alpha + dalpha
      goto 50
c...  we come here when we have a shift that is ok
120   if(dl.gt.0.0d0) test= dabs(beta(nvec,1)/(dl*alpha*beta0))
      if (test.lt.cnvfc2.or.itera.ge.maxaug)  nconv=1
c...  must now compute residual
      call search (iblkif,ifil)
      call vclr(sol,1,naug)
      call vclr(grad,1,naug)
      do 150 i=1,nvec
      call reads (vec,naug,ifil)
      call daxpy(naug,beta(i,1),vec,1,sol,1)
      call reads (vec,naug,ifil)
      call daxpy(naug,beta(i,1)/alpha,vec,1,grad,1)
150   continue
c...  contribution from the 1st exp vector
      call daxpy(naug,beta0,q(igrad),1,grad,1)
      call daxpy(naug,-eig(1),sol,1,grad,1)
      grad(naug) = ddot(nvec,beta,1,p,1) - eig(1)*beta0
      res=dnrm2(naug,grad,1)
      if(itera.eq.0) thresh=acc*res
      if(res.lt.thresh) nconv=1
c      call outvec (grad,naug,'residual')
c...  now premultiply by approx inverse hessian
      if(nconv.ne.0.and.iter.ne.0) goto 170
      fac=beta0*alpha
      do 160 i=1,nrot
      rr=sol(i)
      sol(i)=0
      if(q(idiag-1+i).gt.1.d19) goto 160
      if(nconv.ne.0.and.iter.eq.0) then
      sol(i)=rr
      if( dabs(rr).lt.1.d-5) goto 160
      buff = (grad(i)-q(igrad-1+i)*fac)/rr
      if(buff.gt.0) q(idiag-1+i) = buff
      goto 160
      end if
      buff = q(idiag-1+i)/alpha - eig(1)
      sol(i) = grad(i) / buff
160   continue
c      call outvec (sol,naug,'premultiplied residual')
170   if(itera.ge.1.and.ipri.ne.0) then
       dumt = seccpu()
       write(iwrite,180)itera,dumt,res,dl,eig(1),alpha,nalph,test
      endif
180   format(1x,i6,f9.2,2f14.8,f44.8,f8.2,' (',i1,')',f10.4)
      if (nconv.eq.1) goto 240
      if (nvec.ge.mv) call caserr('too many iterations in davidson')
c...  construct new expansion vector, orthogonal to all previous
      sol(naug) = 0.0d0
215   iblk=iblkr
      zz=1.0d0/dnrm2(nrot,sol,1)
      call dscal(nrot,zz,sol,1)
      do 220 ivec=1,nvec
      call rdedx(vec,naug,iblk,ifil)
      zz = -ddot(nrot,sol,1,vec,1)
      call daxpy(naug,zz,vec,1,sol,1)
220   iblk=iblk+nblk
      zz = 1.0d0/dnrm2(nrot,sol,1)
      call dscal(nrot,zz,sol,1)
      if(zz.gt.1.d3) then
      write(iwrite,225)zz
225   format(' large normalization factor in aughes:', e12.2)
      goto 215
      end if
      iblk=iblkg
      do 230 i=1,nvec
      call rdedx(grad,naug,iblk,ifil)
      iblk=iblk+nblk
      zm(i,nvec+1) = ddot(naug,grad,1,sol,1)
230   zm(nvec+1,i) = zm(i,nvec+1)
      nvec = nvec + 1
      call wrt3s (sol,naug,ifil)
      isigma=0
      call orb2(q(1),sol,q(izint),q(izint1),0)
      call orb3(q(1),sol,grad,q(izint),q(igam))
      call vsub(grad,1,q(igrad),1,grad,1,nrot)
      grad(naug)=0.0d0
      call wrt3s (grad,naug,ifil)
      zm(nvec,nvec) = ddot(naug,grad,1,sol,1)
      p(nvec) = ddot(naug,q(igrad),1,sol,1)
      goto 40
c
c...  convergence reached
c
240   fac=1.0d0/(beta0*alpha)
      call dscal(nrot,fac,sol,1)
      ddr = dnrm2(n,sol,1)
      sol(naug)=sol1
c
      return
      end
      subroutine augi (q,hes,grad,x,n)
      implicit REAL  (a-h,o-z)
      parameter (rad=0.5d0,tol=0.05d0)
      dimension q(*),hes(n,n),grad(n),x(n)
c      call outsqr (hes,n,n,n,'hessian')
c      call outvec (grad,n,'gradient')
      n1=n+1
      imat1 = icorr(n1**2)
      imat2 = icorr(n1**2)
      ivec1 = icorr(n1)
      ivec2 = icorr(n1)
      ivec3 = icorr(n1)
      ivec4 = icorr(n1)
      alpha = 1.0d0
      call vclr(q(imat1+n1*n),1,n1)
      do 20 it=1,50
      ij1=imat1-1
      do 50 i=1,n
      do 40 j=1,n
40    q(ij1+j)=hes(j,i)
      q(ij1+n1)=0.0d0
50    ij1=ij1+n1
      call dscal(n1**2,(1.0d0/alpha),q(imat1),1)
      call dcopy(n,grad,1,q(imat1+n),n1)
      call dcopy(n,grad,1,q(imat1+n*n1),1)
      q(imat1+n1**2-1) = 0.0d0
c      call outsqr (q(imat1),n1,n1,n1,'augmented hessian')
c
      ifail=0
      call f02abf(q(imat1),n1,n1,q(ivec1),q(imat2),n1,q(ivec2),ifail)
      e = q(ivec1)
      call fmove (q(imat2),x,n)
      fac = 1.0d0/ (alpha*q(imat2+n))
      call dscal(n,fac,x,1)
      dl = dnrm2(n,x,1)
      if ( (it.eq.1.and.dl.lt.rad)
     1 .or.( dabs(dl-rad).lt.tol) ) goto 30
      call fmove(hes,q(imat1),n*n)
      do 10 i=1,n
      q(imat1+(i-1)*(n+1)) = q(imat1+(i-1)*(n+1)) - e*alpha
      do 10 j=1,n
10    q(imat1+(i-1)+(j-1)*n) = q(imat1+(i-1)+(j-1)*n)
     1 - alpha*x(i)*grad(j)
      call fmove (x,q(ivec1),n)
      call dscal(n,(2.0d0*e),q(ivec1),1)
      ifail=0
      call f04atf (q(imat1),n,q(ivec1),n,q(ivec2),q(imat2),n
     1,q(ivec3),q(ivec4),ifail)
      zl = ddot(n,x,1,q(ivec2),1)
      dalpha = (dl**2/zl)*(1.0d0-dl/rad)
      if (alpha+dalpha.lt.1.0d0) dalpha = 0.5d0-0.5d0*alpha
20    alpha = alpha+dalpha
30    call corlsr (imat1)
c      call outvec (x,n,'augi solution')
      return
      end
      subroutine filort (c,s,n,nvec,num,iblk)
      implicit REAL  (a-h,o-z)
      dimension c(n),s(n)
      call search (iblk,num)
      do 10 ivec=1,nvec
      call reads (s,n,num)
      zz = - ddot(n,s,1,c,1)
      call daxpy(n,zz,s,1,c,1)
10    continue
      cc = 1.0d0/dnrm2(n,c,1)
      call dscal(n,cc,c,1)
      return
      end
      subroutine optaug (q,sol,n,iwrite)
c     radtol - fractional tolerance for step length in getting shift
c     cnvfc1 - convergence threshold for desired eigenvector, measured
c              as fractional contribution of latest expansion vector
c              to am solution
c     cnvfc2 - convergence threshold for lower roots
c     mv     - max no. of expansion vectors
c     maxrt - max. no. of states
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (radtol=0.001d0, cnvfc1=0.005d0,cnvfc2=0.01d0,mv=20,
     * maxrt=5)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
INCLUDE(common/mcaddr)
      logical lconv(maxrt)
      common /scra  / zm(mv,mv),zn(mv+1,mv+1),p(mv),beta(mv+1,mv+1),
     1              gamma(mv+1),scrach(mv+1),r(mv+1),eig(mv+1),wks(mv+1)
     2              ,tester(maxrt),accur(maxrt),ovlap(mv+1,mv+1)
      dimension q(*),sol(n)
      ibase = icorr(0)
c
c...  update step radius by evaluating trust ratio
      if (radius .gt. 0.0d0) then
      ratio = (energy-elast) / (enext-elast)
      if ( dabs(ratio-1.0d0) .gt. trust1) radius = radius*tfac1
      if ( dabs(ratio-1.0d0) .lt. trust2) radius = radius*tfac2
      else
      ratio = 0.0d0
      radius = -radius
      end if
      radius = dmin1(radius,0.75d0)
      if (mcprin) write(iwrite,10)ratio,radius
10    format(' entering augmented matrix optimisation',10x
     1,'trust ratio  =',f15.8,5x,' maximum step =',e12.4)
c
      call vclr(ovlap,1,(mv+1)**2)
      call vfill(1.0d0,ovlap,mv+2,mv+1)
c...  disable sparsity code for small cases to avoid nasties
      spars = sparse
      if (n.lt.50) spars = 0.0d0
      ispars = 0
      isparo = 0
c     ibordr = igrad
      iaugmx = 1
      nc = nci
      if (isigma.ne.1) nc = 1
      naug = nrot + nc
      n = nrot
      if (isigma.eq.1) n = nvar
      numv = num4
      nums = numscr
      iblkv=iblk4
      iblks=iblk8
      isig = icorr(naug)
      nroot = 1
      if (isigma.eq.1) nroot = nstate
      itrack = nroot
      nalph = 0
      alpha = 1.0d0
      nvec = 0
      eig(nroot) = 0.0d0
      beta0 = 1.0d0
      dl = 0.0d0
      itera = 0
c
c...  required accuracy of solutions
      do 20 i=1,nroot
      accur(i) = cnvfc2
20    tester(i) = 1.0d0
      accur(nroot) = cnvfc1
c
      if (mcprin) write(iwrite,30)
30    format(' ',5x,'time',3x,'convergence',4x,'next',5x,'damp'
     1,9x,'norm',5x,'eigenvalue')
      goto 120
c
c...    start of davidson iteration
40    itera = itera + 1
      nalph = 0
      dalpha = 1.0d0
c     dll = 1.0d30
c...   set up davidson matrix for a given value of alpha
50    do 60 i=1,nvec
      do 60 j=1,nvec
60    zn(i,j) = zm(i,j) / alpha
      do 70 i=1,nvec
      zn(i,nvec+1) = p(i)
70    zn(nvec+1,i) = p(i)
      zn(nvec+1,nvec+1) = 0.0d0
c      call outsqr (zn,mv+1,nvec+1,nvec+1,'n')
c      call outsqr (ovlap,mv+1,nvec+1,nvec+1,'ovlap')
      ifail=0
      if (spars.le.0.0d0) then
      call f02abf (zn,mv+1,nvec+1,eig,beta,mv+1,scrach,ifail)
      else
      call f02aef (zn,mv+1,ovlap,mv+1,nvec+1,eig,beta,mv+1,scrach,wks
     * ,ifail)
      end if
      do 80 i=1,nroot
80    tester(i) =  dabs(beta(nvec,i))
c      call outsqr (beta,mv+1,nvec+1,nvec+1,'beta')
c...  find the norm of the step vector
      beta0 = beta(nvec+1,nroot)
      do 90 i=1,nvec
      gamma(i) = ddot(nvec,ovlap(1,i),1,beta(1,nroot),1)
 90   continue
      dl =  dabs( 1.0d0/(alpha*beta0)) *
     1   dsqrt(ddot(nvec,gamma,1,beta(1,nroot),1) )
c      print*,'alpha,dl ',alpha,dl
      if ( (alpha.eq.1.0d0.and.dl.le.radius)
     1 .or. ( dabs((dl-radius)/radius).le.radtol)
     > .or. nalph.gt.20
     3 .or. nvec.lt.nroot)
     4 goto 120
c..   this shift not ok, so we must update it
c     dll = dl
      nalph = nalph + 1
      do 110 i=1,nvec
      do 100 j=1,nvec
100   zn(j,i) = zm(j,i) - (1.0d0/beta0)*beta(j,nroot)*p(i)
      zn(i,i) = zn(i,i) - eig(nroot)*alpha
110   r(i) = 2.0d0*eig(nroot)*beta(i,nroot)
c      call outsqr (zn,mv+1,nvec,nvec,'n')
      ifail = 1
      call f04arf (zn,mv+1,r,nvec,gamma,wks,ifail)
      if (ifail.eq.1) call caserr('a not pos def in updating alpha')
      if (ifail.eq.2)write(iwrite,1000)
1000  format(               ' ill conditioning in updating alpha')
c      call outvec (gamma,nvec,'gamma')
      dalpha = ((1-beta0**2)/ddot(nvec,gamma,1,beta(1,nroot),1))
     1              * (1.0d0- dl/radius)
      if (alpha+dalpha.lt.1.0d0) dalpha = 0.5d0-0.5d0*alpha
      alpha = alpha + dalpha
      goto 50
c...  we come here when we have a shift that is ok
c...  must decide which vector to iterate next - itrack
120   if(dl.gt.0.0d0)
     +     tester(nroot)=dabs(beta(nvec,nroot)/(dl*alpha*beta0))
      test = 0.0d0
      nconv = 0
      do 130 i=1,nroot
      test = dmax1(test,tester(i))
      lconv(i) = tester(i).le.accur(i)
      if (lconv(i)) nconv = nconv + 1
130   continue
      if (nconv.eq.nroot) goto 170
140   itrack = mod(itrack,nroot)+1
      if (lconv(itrack)) goto 140
c...  must now compute residual
      call search (iblkv,numv)
      call search (iblks,nums)
      call vclr(sol,1,n)
      do 150 i=1,nvec
      call reads (q(isig),naug,numv)
      fac = - eig(itrack)*beta(i,itrack)
      call daxpy(naug,fac,q(isig),1,sol,1)
      call reads (q(isig),naug,nums)
      fac = beta(i,itrack) / alpha
      call daxpy(naug,fac,q(isig),1,sol,1)
150   continue
c...  contribution from the 1st exp vector i.e. ci vector
      call daxpy(n,beta0,q(igrad),1,sol,1)
      fac = ddot(nvec,beta(1,itrack),1,p,1) - eig(itrack)*beta0
      if (isigma.eq.1) then
      call daxpy(nci,fac,q(icivec),1,sol(nrot+1),1)
      else
      sol(naug) = fac
      end if
c      call outvec (sol,naug,'residual')
c...  now premultiply by approx inverse hessian
      do 160 i=1,naug
      buff = q(idiag-1+i)/alpha - eig(itrack)
      if ( dabs(buff).lt.1.0d-16) buff=1.0d30
160   sol(i) = sol(i) / buff
c      call outvec (sol,naug,'premultiplied residual')
170   if (nconv.eq.nroot) itrack = 0
      if (mcprin) then
       dumt = seccpu()
       write(iwrite,180)itera,dumt,test,itrack,alpha,nalph,dl,
     + (eig(i),i=1,nroot)
      endif
180   format(i3,f9.2,f14.8,i7,f10.3,'(',i1,')',f10.5,5f10.4)
      if (nconv.eq.nroot) goto 240
      if (nvec.ge.mv) call caserr('too many iterations in davidson')
c...  construct new expansion vector, orthogonal to all previous
      if (isigma.eq.1) then
      cc = - ddot(nci,sol(nrot+1),1,q(icivec),1)
      call daxpy(nci,cc,q(icivec),1,sol(nrot+1),1)
      else
      sol(nrot+1) = 0.0d0
      end if
      call filort (sol,q(isig),naug,nvec,numv,iblkv)
c =================================
c..   this next block of code to make the vector sparse
      if (spars.gt.0.0d0.and. nvec.ne.0) then
      do 190 i=1,naug
      if ( dabs(sol(i)).ge.spars) goto 190
      sol(i) = 0.0d0
      if (i.le.nrot) isparo = isparo + 1
      ispars = ispars + 1
190   continue
c...  vector must be orthogonal to ci
      if (isigma.eq.1) then
      bc = 0.0d0
      cc = 0.0d0
      do 200 i=1,nci
      if (sol(i+nrot).ne.0.0d0) cc = cc + q(icivec-1+i)**2
200   bc = bc + sol(i+nrot)*q(icivec-1+i)
      if (cc.gt.0.0d0) bc = - bc/cc
      do 210 i=1,nci
      if (sol(i+nrot).ne.0.0d0) sol(i+nrot) = sol(i+nrot)
     1                      + q(icivec-1+i) * bc
210   continue
      end if
      cc = dnrm2(n,sol,1)
      if (cc.le.0) call caserr('singularity in sparsity algorithm')
      cc = 1.0d0/cc
      call dscal(n,cc,sol,1)
c...  extend the overlap matrix
      call search (iblkv,numv)
      do 220 i=1,nvec
      call reads (q(isig),naug,numv)
      ovlap(i,nvec+1) = ddot(naug,sol,1,q(isig),1)
220   ovlap(nvec+1,i) = ovlap(i,nvec+1)
      end if
c =================================
      call search (iblks,nums)
      do 230 i=1,nvec
      call reads (q(isig),naug,nums)
      zm(i,nvec+1) = ddot(naug,q(isig),1,sol,1)
230   zm(nvec+1,i) = zm(i,nvec+1)
      nvec = nvec + 1
      call wrt3s (sol,naug,numv)
      call sigma (q(1),sol,q(isig),iwrite)
      call wrt3s (q(isig),naug,nums)
      zm(nvec,nvec) = ddot(naug,q(isig),1,sol,1)
      p(nvec) = ddot(naug,q(igrad),1,sol,1)
      goto 40
c
c...  convergence reached -- pick it all up
240   continue
      call search (iblkv,numv)
      call vclr(sol,1,n)
      do 250 i=1,nvec
      call reads (q(isig),naug,numv)
      fac = beta(i,nroot)/(beta0*alpha)
      call daxpy(n,fac,q(isig),1,sol,1)
250   continue
      dl = dnrm2(n,sol,1)
      if (dl.gt.radius) then
       call dscal(n,(radius/dl),sol,1)
           endif
c...   search for other low-lying eigenvalues of matrix
      eigmax = dmax1( dabs(eig(nroot)),0.01d0)
      write(iwrite,260)
260   format(' other low lying eigenvalues, and approximate accuracy')
      i = nroot+1
270    write(iwrite,280)i,eig(i), dabs(beta(nvec,i))
280   format(i3,f20.7,e15.3)
      i = i+1
      if (i.le.nvec.and.eig(i).lt.eigmax) goto 270
      sparc =  dfloat(ispars-isparo)  /  dfloat(nc*(nvec-1))
      sparo =  dfloat(isparo)  /  dfloat(nrot*(nvec-1))
c
      if (mcprin) write(iwrite,290)sparo,sparc
290   format(' augmented matrix optimisation finished'
     1,10x,'average sparsity = ',2pf5.1,'% (orbital rotations)'
     2,f10.1,'% (ci rotations)')
      call corlsr (ibase)
c
      enext = energy + eig(nroot) * (1.0d0/alpha + dl**2*alpha )
      return
      end
      subroutine optint(q,iq,idgot,iwrite,ipunch)
      implicit REAL  (a-h,o-z)
      logical orbopt,canoni,cdone,prnt
      logical btest
INCLUDE(common/sizes)
INCLUDE(common/jobopt)
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
INCLUDE(common/multic)
INCLUDE(common/syminf)
      common/coruse/ ncruse
      common/mceig/eigmc(5)
INCLUDE(common/mcscra)
      dimension q(*),iq(*),ijpos(8),iupos(8),eig(5)
c
      prnt=mcprin
      idgot=0
      mxitc=maxitc+2
      cdone=.false.
      if (btest(itinfo(iter+1),1).or.iguess.le.0) then
      mxitc=maxitc+3
      orbopt = .false.
      maxdav = 20
      canoni=iguess.eq.0.and.nrottu.eq.0.and.icang.ne.0
      if(canoni) maxdav=2
      else
      if (nrotit+nrottu.eq.0) return
      orbopt = .true.
      maxdav = 1
      canoni=.false.
      end if
c
      call accnt('internal',1)
      if (mcprin) then
       dumt = seccpu()
       write(iwrite,10)dumt
      endif
10    format(/1x,'entering internal optimisation at time',f9.2)
c
c...  disc offset on file 8 for ci2 workspace - after diis
c...  load transformation matrix
      ibase = icorr(0)
      do 20 isym=1,nirrr
      junk = nprm(isym)**2
      iupos(isym) = icorr(junk)
      call vclr(q(iupos(isym)),1,junk)
      call vfill(1.0d0,q(iupos(isym)),nprm(isym)+1,nprm(isym))
20    continue
25    call ofset(0)
c
      icivec=icorr(nstate*nci)
      igam=icorr(ne)
      izint=igam
c
c...  load initial internal integrals
      ijipos=icorr(lenj)
      call loadi (q(1),q(ijipos))
c
c...  initialise for rdfrt/wtfrt
      if (orbopt.or.cdone) call finit
c
c...  load ci vector
      if (iguess.ne.0) then
      call cget(q(icivec),nstate)
      end if
c
c...  begin optimisation
      ndim=0
      ddr=0.0d0
      ddc=0.0d0
      varc=1.0d0
      itoo=0
      if (orbopt.and.iguess.gt.0) call densav (q(1),q(icivec),q(igam))
      if (orbopt.and.iguess.gt.0) call symden(q(igam),nact)
c
      if(prnt) write(iwrite,40)
      if(prnt) prnt=.false.
40    format(// 12x,'time',2x,'orb. gradient  orb. change',
     1   '   ci gradient    ci change         energy      diis'/
     2 1x,94('-')/)
      do 160 itc=1,mxitc
      call ofset(0)
      ijjpos=icorr(lenj)
      ikkpos=icorr(max(lenj,lenk))
c
      if(.not.orbopt) goto 125
      do 120 ito=1,maxito
      l = icorrm()
      iad = icorr(-l)
c...  transform internal integrals
      call fmove(q(ijipos),q(ijjpos),lenj)
      call traop (q(iupos(1)),q(ijjpos),q(ikkpos),q(iad),l,luse,0,
     +            iwrite)
      ncruse=max(ncruse,iad+luse-1)
      call corlsr (iad)
      ihes=icorr(nrotti**2)
      igrad=icorr(nrotti)
      ig=icorr(nprim**2)
      if=icorr(ntqg)
      ix=icorr(nrotti)
      id1=icorr(nact**2+1)
      id2=icorr(nact**2)
      id3=icorr(nact**2)
c...  calculate gradient and hessian matrix
      call goper(
     *  q(1),    q(ijjpos),q(ikkpos),q(if),q(ig),q(id1),q(id2),q(id3),
     1           q(igam),q(igrad),q(ihes))
c...  augmented hessian procedure
      call augi (q(1),q(ihes),q(igrad),q(ix),nrotti)
      varr = 2.0d0*dnrm2(nrotti,q(igrad),1)
      if(ito.eq.1) then
      if (mcprin) then
       dumt = seccpu()
       write(iwrite,60)itc,itoo,dumt,varr,ddr,varc,
     +                 ddc,energy,ndim
      endif
60    format(1x,2i3,f9.2,4f14.8,1f16.8,i6)
      if (varr+varc.le.1.0d-5) goto 170
      if(itc.eq.mxitc) goto 170
      ddr=0
      ddc=0
      end if
      ddr = ddr + dnrm2(nrotti,q(ix),1)
      itoo=ito
      call updui (q(1),q(iupos(1)),q(ix))
      call corlsr (ihes)
      if (ddr.lt.1.0d-5) goto 130
120   continue
      goto 130
125   l = icorrm()
      iad = icorr(-l)
c...  dummy transform internal integrals
      call fmove(q(ijipos),q(ijjpos),lenj)
      call traop (q(iupos(1)),q(ijjpos),q(ikkpos),q(iad),l,luse,0,
     +            iwrite)
      ncruse=max(ncruse,iad+luse-1)
      call corlsr (iad)
130   if (isigma.ne.1.and.orbopt) goto 170
      varc=0
      if(itc.gt.mxitc-2) goto 160
      call mkzint(q(ijjpos),q(izint))
      call corlsr (ijjpos)
      isg = icorr(nci*nstate)
      call ci2(q(1),iq(1),q(icivec),q(isg),q(izint),eig,idgot,
     +         maxdav,ciacc,iwrite)
      idgot=1
      varc = dnrm2(nci*nstate,q(isg),1)
      if (idsci.eq.0.or.varc.lt.1.d-7.or.(.not.orbopt)) goto 150
      if (varc.gt.disvar) then
      ndim=0
      goto 150
      end if
      call cphase (q(icivec),q(isg))
      call anilp (q(icivec))
      nd1=ndim+1
      junk=nd1*(nd1+1)/2+2*nd1+nd1**2
      ids=icorr(junk)
      call diismc(q(1),q(icivec),q(isg),q(isg),q(ids),nci*nstate,0,
     1   ndim,ndel,bfak,maxdis,iblk8-1)
      call creatp (q(icivec))
      call orthci(q(icivec),nci,nstate,iwrite)
150   call dci (q(icivec),q(isg),eig,iwrite)
      ddc=dnrm2(nci*nstate,q(isg),1)
      call corlsr (ijjpos)
c
c...  new density matrix
c
      if(canoni) goto 170
      if (btest(itinfo(iter+1),7).and..not.orbopt) then
      if (nrotit+nrottu.eq.0) goto 170
      if(.not.cdone) call finit
      orbopt = .true.
      maxdav = 1
      canoni=.false.
      end if
      if(.not.orbopt) goto 170
      call densav (q(1),q(icivec),q(igam))
      call symden(q(igam),nact)
c
160   continue
c
c...  dump ci vector
170   call cput(q(icivec),nstate)
      call corlsr (icivec)
      if(canoni) then
      call canorb(q(1),iq(1),q(iupos(1)),1,1,0,0,0,iwrite,
     +            ipunch)
c      call ofset(0)
_IF1()c      call priop(q(iupos(1)),1,0,'u-canon',0,iwrite)
      if (mcprin) write(iwrite,185)
185   format(/' internal orbitals canonicalized'/)
      cdone=.true.
      canoni=.false.
      maxdav=20
      iguess=0
      idgot=0
      goto 25
      end if
c
      if (orbopt.or.cdone) then
c...  transform orbitals
      do 190 isym=1,nirrr
190   ijpos(isym) = icorr(nsymao(isym)**2)
      ka = icorr(maxbas**2)
      call qget (q(1),ijpos)
      do 210 isym=1,nirrr
      ja=ijpos(isym)+ifreez(isym)*nsymao(isym)
      call fmove (q(ja),q(ka),nsymao(isym)*nsymm(isym))
210   call mxmaa(q(ka),1,nsymao(isym),
     1 q(iupos(isym)),1,nprm(isym),
     2 q(ja),1,nsymao(isym), nsymao(isym),nprm(isym),nprm(isym) )
      call qput (q(1),ijpos,iwrite)
      call corlsr (icivec)
c
c...  transform operators
      l = icorrm()
      iad = icorr(-l)
      call traop (q(iupos(1)),q(iad),q(iad),q(iad),l,luse,1,iwrite)
      end if
c
      if (mcprin) then
      if(nstate.gt.1.or.(.not.orbopt))write(iwrite,219)(eig(i)+core,i=1,
     * nstate)
      endif
      do i=1,nstate
         eigmc(i) = eig(i)+core
      end do
219   format(/' final state energies:',5f22.12)
      ncruse=max(ncruse,iad+luse-1)
      call corlsr (ibase)
      if (mcprin) then
       dumt = seccpu()
       write(iwrite,220)dumt
      endif
220   format(/' end of internal optimisation at time  ',f9.2)
      call accnt(' ',1)
      return
      end
      subroutine optnr (q,sol,n,iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/mcaddr)
      dimension q(*),sol(n)
      write(iwrite,1000)
1000  format ( ' entering newton-raphson optimization' )
c
      accur = 0.0d0
      do 1011 i=1,n
1011  accur = dmin1(accur, dabs(q(idiag-1+i)))
      accur = 0.000001d0*gradnt / dmax1(accur,1.0d-2)
      accur = dmin1( dabs(accur),1.0d-4)
      accur = dmax1(accur,1.0d-12)
c
      north = 0
      itrial = 0
      if (isigma.eq.1) then
      north = 1
      call vclr(sol,1,nrot)
      call cget(sol(nrot+1),1)
      call wrt3 (sol,nvar,iblk4,num4)
      end if
      call mcpopl(q(1),q(igrad),sol,q(idiag),n,num4,iblk4,accur,north,
     *itrial,1,iwrite)
c
c..   nr eqns h.s = -g  so put in - sign into soln
      call dscal(n,-1.0d0,sol,1)
      enext = energy +  ddot(n,q(igrad),1,sol,1)
      return
      end
      subroutine optwrn (q,iq,sol,n,nittra,idgot,iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common /block / nosymm,nosym(8),iosym(8,maxorb)
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
INCLUDE(common/count)
      logical skip,copt,dispol,augmnt,cnotcv
      common/ cpos / iupos(8),iapos(8)
      common/lsort/ v(maxorb),ijpos(8),ikpos(8),izpos(8),ixpos(8),w(31)
     1 ,iprev(8),jprev(8),idum(8),impos(8),inpos(8)
      common /mcopt / var,varc,thzr,one,dec
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
INCLUDE(common/mcaddr)
      dimension q(*),iq(*),sol(n),eig(5)
      call accnt('optwrn',1)
      if(disvar.le.0) maxdis=0
      if(maxaug.eq.0) idsci=min(idsci,1)
      if(nittra.eq.0.and.mcprin) then
       dumt = seccpu()
       write(iwrite,10)dumt
      endif
10    format(/' entering werner-knowles optimization at time:',f8.2)
      if(nittra.eq.0.and.iter.eq.0.and.mcprin) 
     +    write(iwrite,11)idsci,disvar,
     +    maxdis,augvar,maxaug
11    format(/' diis scheme',i2,':     start threshold:',f8.2,
     1 '  max. steps:',i4/
     2 ' ah scheme:         start threshold:',f8.2,'  max. steps:',i4)
      pi=3.141592654d0
c
c...  required accuracy of solution
c
      accur=dmin1(0.1d0*gradnt**2,gradnt*0.01d0)
      if(iter.ge.1) accur=dmin1(accur,1.d-5)
      accur=dmax1(accur,varmin)
c
      cnotcv=.true.
      copt=isigma.eq.1
      loadz=0
      if(copt) loadz=1
      isigm1=isigma
      if (copt.and.nittra.eq.0.and.iter.eq.0.and.mcprin) 
     + write(iwrite,30)copvar,icimax, icimx1,icimx2,cishft
30    format(' ci optimization:   start threshold:',f8.2,
     1   '  max. steps:',3i4//
     2       ' shift parameter for q-space optimization = ',f6.2)
      if(nittra.eq.0.and.iter.eq.0.and.mcprin) 
     + write(iwrite,31)drdamp,drmax
31    format(' damping parameter for delta r =            ',f6.2/
     1       ' maximum step for delta r =                 ',f6.2/)
c
      iaugmx = 0
      iwrnr  = 1
      dispol=.false.
c
c...  core allocation
c
      n = nrot
      if(copt) n=n+nci*nstate
      ibase  = icorr(0)
      ird = icorr(nrot)
      iro = icorr(nrot+1)
      isg = icorr(n+1)
      iwg = isg
      if(igwgt.ne.0) iwg = icorr(nrot)
      if(disvar.le.0) maxdis=0
      do 1 isym=1,nirrr
   1  iapos(isym)=icorr(nsymm(isym)*nprm(isym))
      ig=igrad-1
      id=idiag-1
c
c...  get trial solution
c
      iblk1=2*maxaug*lensec(nrot+1)+iblk8-1
      negd2e=0
      rmax1 = pi / (2.0d0*gfak1)
      if(nittra.gt.0) goto 80
      do 70 i=1,nrot
      call orblab (i,i1,j1,isym)
      if(iosym(isym,i1).ne.iosym(isym,j1)) then
      q(ig+i)=0
      q(id+i)=1.d20
      sol(i)=0
      if(ipri.ne.0) write(iwrite,51)i1,j1,isym
51    format(' start: orbital rotation=',i3,',',i3,' (sym ',i3,
     1   ')  eliminated.')
      goto 70
      end if
      if( dabs(q(ig+i)).lt.1.d-9) q(ig+i)=0
      e1=q(ig+i)
      e2=q(id+i)
      if( dabs(e2).lt.1.d-7) goto 60
      sol(i)=-e1/e2
      if( dabs(sol(i)).le.rmax1.and.e2.gt.0) goto 70
      if(e2.lt.0.0d0) negd2e=negd2e+1
      alpha =  datan(gfak1*e1/e2)
      if (e2.lt.0.0d0) alpha =  dsign(pi- dabs(alpha),e1)
      sol(i) = -alpha/gfak1
      if( dabs(sol(i)).gt.drmax) sol(i)= dsign(drmax,sol(i))
      q(id+i) =  dabs(e2/ dcos(alpha))
      write(iwrite,50)i1,j1,isym,e1,e2,q(id+i),sol(i)
50    format(' start: orbital rotation=',i3,',',i3,' (sym ',i1,')',
     > '   e1=',e10.3,'   e2=',e10.3,'   e2new=',e10.3,'   r=',e10.3)
      goto 70
60    sol(i)=0
      write(iwrite,61)
61    format(1x,'***warning, small diagonal hessian element')
      q(id+i)= dabs(e1)*10.0d0
      if(e1.eq.0) q(id+i)=1.0d0
      write(iwrite,50)i1,j1,isym,e1,e2,q(id+i),sol(i)
70    continue
c
80    if(ipri.eq.0) goto 90
       call outvec (q(igrad),nrot,'orbital gradient')
       call outvec (q(idiag),nrot,'diagonal elements of hessian')
       call outvec (sol,nrot,'first estimate of orbital rotations')
90     continue
      call vclr(q(iro),1,nrot)
      call vclr(q(ird),1,nrot)
      if (copt) then
c
c...  read ci vectors
c
      call cget(sol(nrot+1),nstate)
      if(ipri.ne.0)
     1 call outvec (sol(nrot+1),nci,'ci-coefficients')
      end if
      var=dnrm2(nrot,q(igrad),1)
      ddr=dnrm2(nrot,sol,1)
c
c.....start values
c
      accurr=varmax
      if(iter.ge.1.and.negd2e.eq.0) accurr=0.0d0
      if(nittra.eq.0.and.mcprin) write(iwrite,20)accur,accurr
20    format(/' convergence thresholds:',f12.8,
     1      ' (gradient)',f12.8,' (step)')
      accurc=0.1d0*accur
      if(accurr.eq.0.0d0) maxci=1
      if(accurr.eq.0.0d0) icstep=1
      icmax=icimax
      if(iter.gt.0)   icmax=icimx1
      if(nittra.ne.0) icmax=icimx2
      if(isigma.ne.1) icmax=0
      varr=var
      varc=0
      ddc=0
      dd=ddr
      ddl=ddr
      ndiso=0
      ndim=0
      idlin=0
      ideltr=0
      iclin=2
      if (iter.eq.0 .and. nittra.eq.0) iclin=icstrt
      icopt=0
      eold=energy
      depi=1.0d0
      itaug=0
      skip=.true.
      augmnt=.false.
c
      if (negd2e.gt.0) write(iwrite,120) negd2e
120   format(/i3,' negative diagonal elements in hessian matrix')
c
c...  weight factors for diis
c
      if(igwgt.le.0) then
      iwgt=0
      else
      iwgt=nrot
      call fmove(q(idiag),q(iwg),nrot)
      end if
c
c.....perform microiterations
c
      if(mcprin) write(iwrite,130)
130   format(/  12x,'time',2x,'orb. gradient  orb. change',
     1   '   ci gradient    ci change         energy      diis',
     2   '         de            de(ci)'/1x,127('-')/)
      do 220 linit=1,itmaxr
      if(dd.le.disvar) idlin=max(linit,idstrt,idlin)
      dispol=linit.eq.idlin
      augmnt=augmnt.or.dd.le.augvar
      idtyp=min(1,idsci)
      if(augmnt) idtyp=min(2,idsci)
      copt=linit.ge.iclin.and.dd.lt.copvar.and.icopt.lt.icmax.and.cnotcv
      isigma=0
      if(copt) isigma=1
      n=nrot
      if(copt) then
      n=nrot+nci*nstate
      iclin=linit+icstep
      end if
c
c.....one index integral transformation
c
      call vsub(sol,1,q(ird),1,q(isg),1,nrot)
      call dcopy(nrot,sol,1,q(ird),1)
      call orb2(q(1),q(isg),q(izint),q(izint1),loadz)
      loadz=0
      ddc=0
      dec=0
      if(copt) then
c
c.....ci gradient
c
      icopt=icopt+1
      if (iter.eq.0.and.nittra.eq.0.and.icopt.eq.3) idgot=0
      call ci2(q(1),iq(1),sol(nrot+1),q(isg+nrot),q(izint1),eig,
     +         idgot,maxci,accurc,iwrite)
      idgot=1
      varc=dnrm2(nci*nstate,q(isg+nrot),1)
c
      if(dispol.and.idtyp.eq.2) then
      if(varc.lt.accurc .or. maxci.gt.1) goto 135
      if(varc.gt.disvar) then
      ndim=0
      goto 135
      end if
      idlin=linit+idstep
c.....extrapolate q-space
      call cphase(sol(nrot+1),q(isg+nrot))
      call anilp(sol(nrot+1))
      if(ndiso.ge.0.or.ndim.eq.maxdis) ndim=0
      ndiso=-1
      nd1=ndim+1
      ids=icorr(nd1*(nd1+1)/2+2*nd1+nd1**2)
      call diismc(q(1),sol(nrot+1),q(isg+nrot),q(iwg),q(ids),nci*nstate
     1    ,iwgt,ndim,ndel,bfak,maxdis,iblk1)
      call corlsr(ids)
c.....add p-space and orthonormalize
      call creatp(sol(nrot+1))
      call orthci(sol(nrot+1),nci,nstate,iwrite)
      end if
c.....delta c
135   call dci(sol(nrot+1),q(isg+nrot),eig,iwrite)
      ddc=dnrm2(nci*nstate,q(isg+nrot),1)
      cnotcv= dabs(dec).gt.econv*1.0d-3
      call densav(q(1),sol(nrot+1),q(igam))
      end if
c
c.....orbital gradient
c
      call orb3(q(1),q(isg),q(isg),q(izint),q(igam))
c
c...  test for convergence
c
      varr=dnrm2(nrot,q(isg),1)
      var=varr+varc
      depo=depi
      dep=enext-eold
      depi= dabs(dep)
      if(linit.ge.3) then
      if(icmax.gt.0.and.icopt.eq.0) goto 180
      dep2= dabs((energy-enext)**3)
      dep2=dmax1(dep2,econv*1.d-3)
      dep2=dmin1(dep2,1.d-6)
      if(depi+depo.lt.dep2) goto 230
      if(var.lt.accur .or. dmax1(dd,ddl).lt.accurr) goto 230
      end if
c
c...  augmented hessian method update of orbital rotations
c
180   if(augmnt) then
      ideltr=1
      call vclr(sol,1,nrot)
      call vclr(q(ird),1,nrot)
      ee = enext
      call aughes(q(1),sol,q(isg),q(iro),ddr,itaug,iwrite)
      enext = ee
      ideltr=0
      goto 190
      end if
c
c...  diis extrapolation
c
      if(linit.le.ipri) call outvec(q(isg),n,'gradient before diis')
c.....q(isg)=gradient at this stage
      if(.not.dispol.or.idtyp.gt.1) goto 150
      idlin=linit+idstep
c.....extrapolate only on q-space
      ndis=nrot
      if(copt.and.idsci.ge.1) ndis=n
      if((ndim.eq.maxdis).or.(ndis.ne.ndiso)) ndim=0
      ndiso=ndis
      if(ndis.gt.nrot) call cphase(sol(nrot+1),q(isg+nrot))
      if(ndis.gt.nrot) call anilp(sol(nrot+1))
      nd1=ndim+1
      ids=icorr(nd1*(nd1+1)/2+2*nd1+nd1**2)
      call diismc(q(1),sol,q(isg),q(iwg),q(ids),ndis,iwgt,ndim,
     1  ndel,bfak,maxdis,iblk1)
      call corlsr(ids)
      if(ndis.gt.nrot) then
c.....add p-space and orthonormalize
      call creatp(sol(nrot+1))
      call orthci(sol(nrot+1),nci,nstate,iwrite)
      end if
      if(linit.gt.ipri) goto 150
      call outvec(q(isg),n,'gradient after diis')
      call outvec(sol,n,'solution after diis')
c
c.....delta r
c
150   call rvec(sol,q(iro),q(isg),q(igrad),q(idiag),nrot,skip)
      skip=.false.
      ddr=dnrm2(nrot,q(isg),1)
      if(linit.gt.ipri) goto 190
      call outvec(q(idiag),nrot,'diagonal elements of hessian')
      call outvec(q(isg),n,'change of solution')
      call outvec(sol,n,'new solution')
 190  if(mcprin) then
       dumt = seccpu()
       write(iwrite,200)linit,itaug,dumt,varr,ddr,varc,
     +                  ddc,enext,ndim,dep, dec
      endif
200   format(1x,2i3,f9.2,4f14.8,1f16.8,i6,1x,2f16.10)
      ddl=dd
      dd=ddr+ddc
      if(icopt.ge.icmax.or..not.cnotcv) varc=0
220   eold=enext
c
c...  approximate r=0.5*[u-u(dagger)]
c
230   call uasm(q(1),iupos,sol)
      isigma=isigm1
c
      if (icopt.gt.0) then
c...  replace dumpfile copy of ci vectors; thus return no change in sol
      call cput(sol(nrot+1),nstate)
      call vclr(sol(nrot+1),1,nci*nstate)
c...  print energies
      eshift = enext - ddot(nstate,weight,1,eig,1)
      if (nstate.gt.1) write(iwrite,250)(eig(i)+eshift,i=1,nstate)
250   format(/' state energies:',5f17.8)
      end if
c
      if(ipri.eq.0) goto 260
      call outvec(sol,n,'solution')
c
260   call corlsr(ibase)
      iwrnr = 0
      ideltr=0
      if (lto(5)) then
      write(iwrite,1000)
1000  format(' returning from optwrn; solution: ')
      do 245 i=1,nrot
      call orblab (i,i1,j1,isym)
245     write(iwrite,246)i1,j1,isym,sol(i)
246   format(i4,',',i3,' (sym ',i1,')',f20.12)
      end if
      call accnt(' ',1)
      return
      end
      subroutine mcpopl (q,r,s,diag,n,num,iblnum,accur,north,itrial,
     + iprint,iwrite)
c
c      pople's linear equation solver - version for a single rhs
c      access to matrix through subroutine sigma
c      (approximate) diagonal elements supplied through diag
c      on entry, rhs in r; preserved on exit
c      on exit, solution in s
c      atmol work file given by num; if north > 0 on entry, the
c      routine will search for the solution in the space orthogonal
c      to the first north vectors on file num.
c      if itrial=1, trial solution has been supplied in s
c
c      external routines required:
c      icorr etc. - dynamic allocation of storage
c      f04atf     - nag routine to solve linear equations
c      mxmb       - matrix multiplier
c      sdot       - scaler product
c      sscal      - scale a vector by a scalar
c      zero       - zero out a vector
c      outsqr,outvec - routines to print matrix and vector
c
c      search,read,reads,wrt3,wrt3s - atmol i/o
c                                               pjk, 10-2-84
      implicit REAL  (a-h,o-z)
      parameter (maxdim=20)
      common /scra  / zm(maxdim,maxdim),p(maxdim),alpha(maxdim)
      dimension q(*),r(n),s(n),diag(n)
c
      if (itrial.ne.1) call vclr(s,1,n)
      call vclr(zm,1,maxdim*maxdim)
      if (north.ge.n) return
      nrstrt = 0
      ivec = 0
      iter = 0
      iblock =  lensec(n)     * north + iblnum
      iv = icorr(n)
      iconv = 2
c
      if (itrial.ne.1) then
_IF(cray,ibm,vax)
      do 10 i=1,n
   10 s(i) = r(i)/diag(i)
_ELSE
      call vdiv(r,1,diag,1,s,1,n)
_ENDIF
      end if
c
c...  start of iteration
20    iter = iter + 1
c
c...  create expansion vector from last sigma
      call fmove (s,q(iv),n)
30    call search (iblnum,num)
      ibuff = icorr(n)
      do 40 i=1,north+ivec
      call reads (q(ibuff),n,num)
      zz = - ddot(n,q(ibuff),1,q(iv),1)
40    call mxmb (q(ibuff),1,0, zz,0,0, q(iv),1,0, n,1,1)
      call corlsr (ibuff)
      zz = 1.0d0/dnrm2(n,q(iv),1)
      call dscal(n,zz,q(iv),1)
c...  check on complete span of space
      if (zz.gt.1.0d10) goto 100
c...  check on near annihilation & repeat orthogonality if so
      if (zz.gt.10.0d0) goto 30
      ivec = ivec + 1
      if (ivec.gt.1) zm(ivec,ivec-1) = ddot(n,q(iv),1,s,1)
      call wrt3s (q(iv),n,num)
c...  get rhs
      p(ivec) = 0.0d0
      do 50 j=1,n
50    p(ivec) = p(ivec) + q(iv-1+j)*r(j)/diag(j)
      call search (iblnum,num)
      do 60 i=1,north
      call reads (s,n,num)
      zz = - ddot(n,s,1,r,1)
      do 60 j=1,n
60    p(ivec) = p(ivec) + q(iv-1+j)*zz*s(j)/diag(j)
c
      call sigma (q(1),q(iv),s,iwrite)
_IF(cray,ibm,vax)
      do 70 i=1,n
   70 s(i) = s(i) / diag(i)
_ELSE
      call vdiv(s,1,diag,1,s,1,n)
_ENDIF
      if (iprint.gt.3) then
      call outvec (q(iv),n,'expansion vector')
      call outvec(s,n,'sigma')
      end if
c
c...  interaction with all previous expansion vectors
      call search (iblock,num)
      do 80 jvec=1,ivec
      call reads (q(iv),n,num)
      zm(jvec,ivec) = ddot(n,q(iv),1,s,1)
80    continue
c
c...  solve in the subspace
      ibuff = icorr (ivec*(ivec+2))
      ifail=0
      call f04atf (zm,maxdim,p,ivec,alpha,q(ibuff+ivec*2),ivec,
     1  q(ibuff),q(ibuff+ivec),ifail)
      call corlsr (ibuff)
      if(iprint.gt.2) then
      call outsqr (zm,maxdim,ivec,ivec,'interaction matrix')
      call outvec (p,ivec,'right hand side')
      call outvec (alpha,ivec,'solution')
      end if
      actest =  dabs(alpha(ivec))
      if (iprint.ge.1) write(iwrite,90)iter,actest
90    format(' iteration',i3,5x,'convergence:',e10.2)
      if (actest.lt.accur) goto 100
      iconv = 1
      if (ivec.ge.maxdim) goto 100
      iconv = 2
      go to 20
c
c...  assemble solution
100   call search (iblock,num)
      call vclr(s,1,n)
      do 110 jvec=1,ivec
      call reads (q(iv),n,num)
110   call mxmb (q(iv),1,0, alpha(jvec),0,0, s,1,0, n,1,1)
      go to (120,130),iconv
c
c...  convergence not reached - we try again
120   if (iprint.gt.0) write(iwrite,1000)
 1000 format(' max. no. of vectors reached; restarting')
      nrstrt = nrstrt + 1
      ivec = 0
      go to 20
c
c...  convergence reached - tidy up and quit
130   if (iprint.gt.0) write(iwrite,140)iter,nrstrt,ivec
140   format(' convergence reached in',i3,' iterations ('
     1,i1,' restarts,',i3,' expansion vectors)')
      if (iprint.gt.1) call outvec (s,n,'solution')
      call corlsr (iv)
      end
      subroutine rvec(r,rold,g,gold,d,n,skip)
      implicit REAL  (a-h,o-z)
      logical skip
INCLUDE(common/sizes)
      dimension r(*),rold(*),g(*),gold(*),d(*)
INCLUDE(common/jobopt)
      common /mcopt / var,varc,thzr,one
INCLUDE(common/multic)
c
      gfaki=one/gfak2
c     gmax=0
      do 50 i=1,n
      e2=d(i)
      if(e2.gt.1.d19) then
      r(i)=0
      g(i)=0
      goto 50
      end if
      if(skip) goto 10
      dro=r(i)-rold(i)
      if( dabs(dro).gt.thzr) e2=(g(i)-gold(i))/dro
      if(e2.lt.0.0d0) e2=d(i)
      e2=dmax1(e2,d(i)*gfaki)
      e2=dmin1(e2,d(i)*gfak2)
      d(i)=0.8d0*d(i)+0.2d0*e2
10    dr=-drdamp*g(i)/e2
      rold(i)=r(i)
      gold(i)=g(i)
      drmx=drmax
      if(skip) goto 40
      if(irdamp.eq.0) goto 30
      if( dabs(r(i)).lt.thzr) goto 30
      if( dsign(one,dr).eq. dsign(one,r(i))) goto 20
      drmx=dmin1(drmax, dabs(r(i))*1.7d0)
      goto 30
20    drmx=dmin1(drmax, dabs(r(i)))
30    if( dabs(dro).lt.thzr) goto 40
      if( dsign(one,dr).eq. dsign(one,dro)) goto 40
      drmx= dabs(dro)*gfak3
40    if( dabs(dr).gt.drmx) dr= dsign(drmx,dr)
      g(i)=dr
50    continue
      call vadd(r,1,g,1,r,1,n)
      return
      end
      subroutine mcchek(q,iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
      logical change
      dimension q(*),spin(5),dum(1)
      if (iexc.ge.0) return
      change = .false.
      iv1 = icorr(nci*nstate)
      iv2 = icorr(nci*nstate)
      call cget(q(iv1),nstate)
      na=0
      nb=0
      do 20 i=1,nact
      ii = ibfcod(i)
      if (ii.eq.5.or.ii.eq.7) na=na+1
      if (ii.eq.5.or.ii.eq.8) nb=nb+1
20    continue
      mulmin = na-nb+3
      mulmax = na+nb+1
      ss =  dfloat(na-nb)*0.5d0
      ss = ss * (ss+1.0d0)
      do 60 mul=mulmin,mulmax,2
      s =  dfloat(mul-1)*0.5d0
      call detci (q(1),q(1),dum,dum, q(iv1),dum,q(iv2),11,nstate)
c      call outvec (q(iv1),nci*nstate,'c')
c      call outvec (q(iv2),nci*nstate,'s**2 c')
      dev=0.0d0
      do 30 istate=1,nstate
      iiiv1 = iv1+(istate-1)*nci
      iiiv2 = iv2+(istate-1)*nci
      spin(istate) = ddot(nci,q(iiiv1),1,q(iiiv2),1)
30    dev = dmax1(dev, dabs(ss-spin(istate)))
      if (dev.le.1.0d-9) goto 70
      write(iwrite,40)(spin(istate),istate=1,nstate)
40    format(' *** warning *** <s**2> =',5f13.9)
      write(iwrite,50)s
50    format(' annihilate s = ',f4.1,' component')
      call daxpy(nci*nstate
     + ,(-s*(s+1.0d0)),q(iv1),1,q(iv2),1)
      call dcopy(nci*nstate,q(iv2),1,q(iv1),1)
c      call outvec (q(iv1),nci*nstate,'new c')
      call orthci (q(iv1),nci,nstate,iwrite)
c      call outvec (q(iv1),nci*nstate,'orthogonalised new c')
      change = .true.
60    continue
70    if (.not.change) goto 100
      write(iwrite,80)(spin(istate),istate=1,nstate)
80    format(' final values of <s**2>  ',5f13.9)
      call cput(q(iv1),nstate)
100   call corlsr (iv1)
      return
      end
      subroutine mctran(q,iq)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mcff)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
INCLUDE(common/atmblk)
      common /mctrns/ master,isym12,isym1,itran1,itran2,ntran
     >               ,ires,maxt,ibl5,ibl54,ibl56,nbuck,nteff
     >               ,mloww,mhi,ntri,nwb,ibuck
     >               ,mark(200),nwbuck(200),ioffpr(8,8)
      common /stak / btri,mlow,nstack,iblock
      common /bufb  / nwbnwb,lnklnk,gsort(5118)
      common/junk/klsort(2,3412)
      common/craypk/ijout(680)
      common/sortpk/klout(680)
      common /mccore/ intrel
      dimension q(*),iq(*)
      call accnt('4 index ',1)
      ifcpos=icorr(ltri)
      ifinit=0
      iblf=1
      call setsto(680,0,ijout)
      call setsto(680,0,klout)
      call tran1(q(1),q(ifcpos))
c
c...  stage 2 sort k matrix and calc c
      if(o255i) then
       ibl5=((iblsrt-1)*2)/3
      else
       ibl5=((iblsrt-1)*2)/4
      endif
      jword=icorrm()-lenrec*nirrr*2
     1               -lensq*nirrr-lensq
      maxt=jword/lentri
      nword=icorrm()/(1+2/intrel)-intrel
      ires=min(nword/ibl5,200)
      if (ires.lt.1 .or. maxt.lt.1) call caserr('not enough store')
      ntran=nprim*nbasao
      npass=0
10    npass=npass+1
      nteff=(ntran-1)/npass+1
      if ((nteff-1)/ires.ge.maxt) goto 10
      npass=min(npass,ntran)
      if (npass.gt.200) call caserr('too many passes needed')
      master=0
      do 20 ipass=1,npass
      iad1=icorr(-nword)
      iad2=icori(-(nword+nword))
      call trasor (q(iad1),iq(iad2),num6,iblk6)
      call corlsr(iad1)
20    call tran2 (q(1),ipass,npass)
c...  stage 3 sort d matrix and then transform j,k operators to mo
      ntran=(nprim*(nprim+1))/2
      nteff=ntran
      master=0
      jword=icorrm()-lensq*nirrr*2-lensq-nblkq-lentra(1)
      maxt=jword/(lentri*3)
      nword=icorrm()/(1+2/intrel)-intrel
      ires=min(nword/ibl5,200)
      if ((ntran-1)/ires.ge.maxt) call caserr('not enough store')
      iad1=icorr(-nword)
      iad2=icori(-(nword+nword))
      call trasor (q(iad1),iq(iad2),num4,iblk4)
      call corlsr (iad1)
      call tran3(q(1),q(ifcpos))
      call corlsr(ifcpos)
      call accnt(' ',1)
      itrsfm=1
      return
      end
      subroutine tran1(q,fc)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
      common /mctrns/ master,isym12,isym1,itran1,itran2,ntran
     >               ,ires,maxt,ibl5,ibl54,ibl56,nbuck,nteff
     >               ,mloww,mhi,ntri,nwb,ibuck
     >               ,mark(200),nwbuck(200),ioffpr(8,8)
      common /stak/ btri,mlow,nstack,iblock
      common /intbuf/ intpos,intfil,intmod
INCLUDE(common/mcff)
      common /lsort/ gkout(510),mwordk,ispk
     1              ,gjout(510),mwordj,icode
      common/craypk/ijout(680)
      common/sortpk/klout(680)
INCLUDE(common/atmblk)
      dimension iqpos(8),ihpos(8),iipos(8),ijpos(8),impos(8)
      dimension q(*),fc(*)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      call accnt('tran1   ',2)
c..   allocate core
      do 10 isym=1,nirrr
      junk=nsymao(isym)**2
      iqpos(isym)=icorr(junk)
      ihpos(isym)=icorr(lenrec)
      iipos(isym)=icorr(lenrec)
      ijpos(isym)=icorr(lenrec)
10    impos(isym)=icorr(lenprm)
      ibuff=icorr(lentri)
      ig=icorr(lensq)
      call qget (q(1),iqpos)
      do 15 i=1,nirrr
      loop=nsymm(i)*nsymao(i)
      ia=iqpos(i)
      ja=iqpos(i)+ifreez(i)*nsymao(i)
 15   if(ja.gt.ia. and .loop.gt.0)
     * call fmove(q(ja),q(ia),loop)
c
      intpos=0
      intfil=num2
      iblf=iblk2
      mwordk=0
_IFN1(iv)      int4k=-1 
      mwordj=0
_IFN1(iv)      int4j=-1 
      icode=1
      call intin(fc,ltri)
      call search (iblk6,num6)
      call search (iblk4,num4)
      labirr=0
      do 180 isympq = 1,nirrr
      labpq=0
      labir = labirr
      do 170 isymp=1,nirrr
      isymq = mults (isymp,isympq)
      if (isymq.gt.isymp) goto 170
      nsymp = nsymao(isymp)
      nsymq = nsymao(isymq)
      iqmax = nsymq
      do 160 ip=1,nsymp
      ipm = ip -1
      if (isympq.eq.1) iqmax=ip
      do 150 iqq=1,iqmax
c     labpq=labpq+65536
      labpq=labpq+1
      labir=labirr
c...  load integrals (../pq) from disc
      call intin (q(ibuff),lentra(isympq))
      ibuff1 = ibuff
      do 20 isymr=1,nirrr
      call vclr(q(ihpos(isymr)),1,lenrec)
      call vclr(q(iipos(isymr)),1,lenrec)
      call vclr(q(ijpos(isymr)),1,lenrec)
20    continue
c ====================================================================
      if (isympq-1)700,700,200
c
 700  if (isymp.eq.1) go to 400
c....  isympq=1 code
      do 30 isymr = 1,isymp-1
      if(nprm(isymr)*nsymao(isymr).le.0)go to 30
      call trankh (q(ihpos(isymr)),q(ibuff1),q(iqpos(isymr))
     1            ,nprm(isymr),nsymao(isymr),1,nsymao(isymr),q(ig) )
      call tranki (q(iipos(isymr)),q(ibuff1),q(iqpos(isymr))
     1            ,nprm(isymr),nsymao(isymr),1,nsymao(isymr) )
 30   ibuff1 = ibuff1 + (nsymao(isymr)*(nsymao(isymr)+1))/2
c
c....  symr = symp
 400  if(nprm(isymp)*ipm.le.0)go to 500
      call trankh (q(ihpos(isymp)),q(ibuff1),q(iqpos(isymp))
     1            ,nprm(isymp),nsymp,1,ipm, q(ig) )
      call tranki (q(iipos(isymp)),q(ibuff1),q(iqpos(isymp))
     1            ,nprm(isymp),nsymp,1,ipm )
c
c... r=p special code
 500  call mxmb (q(iqpos(isymp)),nsymp,1, q(ibuff1+(ip*ipm)/2),1,nsymp,
     1           q(ihpos(isymp)+ipm),nsymp,1, nprm(isymp),iqq,1)
      call mxmb (q(ibuff1+(ip*ipm)/2),1,nsymp,
     1 q(iqpos(isymp)+ipm),1,nsymp,
     2           q(iipos(isymp)),1,nsymp, iqq,1,nprm(isymp))
      call mxmb (q(iqpos(isymp)+iqq-1),nsymp,1,
     1           q(ibuff1+(ip*ipm)/2+iqq-1),1,nsymp,
     2           q(ijpos(isymp)+ipm),nsymp,1,  nprm(isymp),ip-iqq+1,1)
c
c.... r=p+1,n code
      if(nprm(isymp)*(nsymp-ip).le.0)go to 600
      call trankh (q(ijpos(isymp)),q(ibuff1),q(iqpos(isymp))
     1            ,nprm(isymp),nsymp,ip+1,nsymp ,q(ig) )
600   ibuff1 = ibuff1 + (nsymp*(nsymp+1))/2
c...    symr > symp code
      do 40 isymr=isymp+1,nirrr
      if(nprm(isymr)*nsymao(isymr).le.0)go to 40
      call trankh (q(ijpos(isymr)),q(ibuff1),q(iqpos(isymr))
     1            ,nprm(isymr),nsymao(isymr),1,nsymao(isymr) ,q(ig) )
40    ibuff1 = ibuff1 + (nsymao(isymr)*(nsymao(isymr)+1))/2
c
c
c ==================
      go to 300
c...  isympq .ne. 1 code
200   do 70 isymr=1,isymp
      isyms = mults (isymr,isympq)
      if (isyms.gt.isymr) goto 70
      nsymr = nsymao(isymr)
      nsyms = nsymao(isyms)
      nprimr = nprm(isymr)
      nprims = nprm(isyms)
      ilastr = nsymr
      if (isymr.eq.isymp) ilastr=ip-1
      call mxmb (q(ibuff1),nsyms,1, q(iqpos(isyms)),1,nsyms,
     1           q(ihpos(isymr)),1,nsymr, ilastr,nsyms,nprims)
      call mxmb (q(ibuff1),1,nsyms, q(iqpos(isymr)),1,nsymr,
     1           q(iipos(isyms)),1,nsyms, nsyms,ilastr,nprimr)
      if (isymr.lt.isymp) goto 60
c...  r=p special code
      call mxmb(q(ibuff1+ipm*nsyms),1,nsyms,q(iqpos(isymp)+ipm),1,nsymp,
     1           q(iipos(isyms)),1,nsyms, iqq,1,nprimr)
      do 50 i=1,nprims
      ioffs = ipm+(i-1)*nsymp
      q(ihpos(isymp)+ioffs) = q(ihpos(isymp)+ioffs) +
     1 ddot(iqq,q(ibuff1+ipm*nsyms),1,q(iqpos(isyms)+(i-1)*nsyms),1)
      ix1=ibuff1+ipm*nsyms+iqq-1
50    q(ijpos(isymp)+ioffs) = q(ijpos(isymp)+ioffs) +
     1 ddot(nsyms-iqq+1,q(ix1),1,q(iqpos(isyms)+(i-1)*nsyms+iqq-1),1)
60    if (isymr.ne.isymp) ibuff1 = ibuff1 + nsyms*nsymr
70    continue
c
c... symr > symp code
      do 80 isymr=isymp,nirrr
      isyms=mults(isymr,isympq)
      if (isyms.gt.isymr) goto 80
      ifrstr=0
      if (isymr.eq.isymp) ifrstr=ip
      call mxmb (q(ibuff1+ifrstr*nsymao(isyms)),nsymao(isyms),1,
     1           q(iqpos(isyms)),1,nsymao(isyms),
     2           q(ijpos(isymr)+ifrstr),1,nsymao(isymr),
     3           nsymao(isymr)-ifrstr,nsymao(isyms),nprm(isyms) )
      ibuff1 = ibuff1 + nsymao(isyms)*nsymao(isymr)
80    continue
c ====================================================================
c
c...  we now have h,i,j matrices for this p,q
c...  assemble and write k
300   do 110 isymr=1,isymp
      irmax = nsymao(isymr)
      if (isymr.eq.isymp) irmax = ip
      isyms = mults(isymr,isympq)
      nprms = nprm(isyms)
      ihh = ihpos(isymr) - 1
      iii = iipos(isymr) - 1
      do 100 i=1,nprms
      do 90 ir=1,irmax
      mwordk=mwordk+1
      gkout(mwordk)=q(ihh+ir)+q(iii+ir)
c     kout(mwordk) = labir+ir+labpq
_IFN1(iv)      int4k=int4k+2
_IFN1(iv)      klout(int4k  )=labpq
_IFN1(iv)      klout(int4k+1)=labir+ir
c      print*,'mwordk,gkout,kout,labir,ir,labpq'
c     >,mwordk,gkout(mwordk),kout(mwordk),labir,ir,labpq
      if (gkout(mwordk).eq.0.0d0) then
       mwordk=mwordk-1
       int4k=int4k-2
      endif
      if (mwordk.eq.num2e) then
        call trout (gkout,klout,mwordk,0,num6)
        int4k=-1
      endif
      mwordk=mod(mwordk,num2e)
90    continue
      labir = labir + nsymao(isymr)
      ihh = ihh + nsymao(isymr)
100   iii = iii + nsymao(isymr)
110   continue
c
c...  now compute m matrix (coulomb in ao basis) l will be at i
      do 120 isymr=1,nirrr
      isymi=mults(isymr,isympq)
      loop = nprm(isymi)*nsymao(isymr)
      if(loop.le.0) go to 120
      call vadd(q(ihpos(isymr)),1,q(ijpos(isymr)),1,q(iipos(1)),1,loop)
      call mxmaa(q(iqpos(isymr)),nsymao(isymr),1,
     1  q(iipos(1)),1,nsymao(isymr),
     2 q(impos(isymr)),1,nprm(isymr),
     3 nprm(isymr),nsymao(isymr),nprm(isymi))
120   continue
c
      do 140 isymi=1,nirrr
      nprimi=nprm(isymi)
      isymj=mults(isymi,isympq)
      if (isymj.gt.isymi) goto 140
      nprimj=nprm(isymj)
      jmax=nprimj
      if (nprimi.eq.0 .or. nprimj.eq.0) goto 140
      do 130 i=1,nprimi
      ii=iorbsm(istart(isymi)-1+i)
      if (isympq.eq.1) jmax=i
      do 130 j=1,jmax
      jj=iorbsm(istart(isymj)-1+j)
      mwordj=mwordj+1
      gjout(mwordj)=q(impos(isymj)+j-1+(i-1)*nprimj)
     1             +q(impos(isymi)+i-1+(j-1)*nprimi)
c     jout(mwordj)=ind(ii,jj)+labpq
_IFN1(iv)      int4j=int4j+2
_IFN1(iv)      ijout(int4j  )=labpq
_IFN1(iv)      ijout(int4j+1)=ind(ii,jj)
      if (gjout(mwordj).eq.0.0d0) then
       mwordj=mwordj-1
       int4j=int4j-2
      endif
      if (mwordj.eq.num2e) then
       call trout (gjout,ijout,mwordj,icode,num4)
       int4j = -1
      endif
130   mwordj=mod(mwordj,num2e)
140   continue
c
150   continue
160   continue
170   continue
180   labirr=labir
      call corlsr (iqpos(1))
      if (mwordk.ge.1) call trout (gkout,klout,mwordk,0,num6)
      call put (gkout,0,num6)
      if (mwordj.ge.1) call trout (gjout,ijout,mwordj,icode,num4)
      ntran=nprim*nbasao
      call accnt(' ',2)
      return
      end
      subroutine tran2(q,ipass,npass)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
INCLUDE(common/atmblk)
      common /mctrns/ master,isym12,isym1,itran1,itran2,ntran
     >               ,ires,maxt,ibl5,ibl54,ibl56,nbuck,nteff
     >               ,mloww,mhi,ntri,nwb,ibuck
     >               ,mark(200),nwbuck(200),ioffpr(8,8)
      common /stak/ btri,mlow,nstack,iblock
      common /bufb/ nwbnwb,lnklnk,gsort(5118)
      common /junk/ klsort(2,3412)
      common /lsort/ gout(510),mword,icode
      common /sortpk/ ijout(680)
      dimension q(*),iqpos(8),iapos(8)
c
      call accnt('tran2',2)
c...  addressing of p-r result
      do 10 isympr=1,nirrr
      iad=lentri
      do 10 isymp=1,nirrr
      ioffpr(isympr,isymp)=iad
10    iad=iad+nsymao(isymp)*nsymao(mults(isymp,isympr))
      do 20 isym=1,nirrr
      iapos(isym)=icorr(lenrec)
20    iqpos(isym)=icorr(nsymao(isym)**2)
      call qget (q(1),iqpos)
      do 25 i=1,nirrr
      ia=iqpos(i)+ifreez(i)*nsymao(i)
      ja=iqpos(i)
      loop = nsymm(i)*nsymao(i)
      if(ia.gt.ja.and.loop.gt.0)
     *call fmove(q(ia),q(ja),loop)
 25   continue
      ig=icorr(lensq)
      nsize=ntri*lentri
      ibuff=icorr(-nsize)
      ibuffm=ibuff-1
      isymir=isym12
      isymr=isym1
      i = itran1
      ir = itran2
      if (ipass.gt.1) goto 30
      isymir = nirrr
      do 1 isymi=1,nirrr
      do 1 isymrr=1,nirrr
      if (nsymao(isymrr).le.0 .or. nprm(isymi).le.0)  goto 1
      isymir = min(isymir,mults(isymi,isymrr))
1     continue
      do 2 isymr=1,nirrr
      if (nprm(mults(isymir,isymr)).gt.0) goto 3
2     continue
3     continue
      i=1
      ir=0
30    isymi=mults(isymr,isymir)
      mword=0
_IFN1(iv)      int4 = -1
      icode=0
c
      call stopbk
_IF(parallel)
c **** MPP
      call closbf(0)
      call setbfa(-1)
c **** MPP
_ENDIF
      do 190 ibuck=1,nbuck
      mhigh=min(mhi,mloww+ntri)
      mtri=mhigh-mloww
      mloww=mloww+1
      nsize=mtri*lentri
      call vclr(q(ibuff),1,nsize)
      mk=mark(ibuck)
40    if (mk.eq.9999999) goto 60
      call mcsrti (mk,nk)
      do 50 iword=1,nk
50    q(ibuffm+(klsort(1,iword)-mloww)*lentri+klsort(2,iword)) 
     +          = gsort(iword)
      goto 40
60    mapp=ibuff
      do 180 itri=1,mtri
      master=master+1
70    ir=ir+1
      if (ir.le.nsymao(isymr)) goto 100
80    ir=0
      i=i+1
      if (i.le.nprm(isymi)) goto 70
90    i=0
      isymr=isymr+1
      isymi=mults(isymr,isymir)
      if (isymr.le.nirrr) goto 80
      isymr=0
      isymir=isymir+1
      if (isymir.le.nirrr) goto 90
100   nsymr=nsymao(isymr)
c     nprimi=nprm(isymi)
      ii=iorbsm(istart(isymi)-1+i)
c
      map = mapp
      do 110 isymp=1,nirrr
      if (isymp.lt.isymr.and.isymir.eq.1) map = map +
     1 (nsymao(isymp)*(nsymao(isymp)+1))/2
      if (isymp.lt.isymr.and.mults(isymp,isymir).lt.isymp)
     1       map = map + nsymao(isymp)*nsymao(mults(isymir,isymp))
      call vclr(q(iapos(isymp)),1,lenrec)
110   continue
      do 120 isymp=isymr,nirrr
      isymq=mults(isymp,isymir)
      if (isymq.gt.isymp) goto 120
      nsymp=nsymao(isymp)
      nsymq = nsymao(isymq)
      nprimj=nprm(isymq)
      if (isymir.eq.1) then
c.... isymir.eq.1 code
      irmin = 1
      if (isymp.eq.isymr) irmin = ir
      call trankh (q(iapos(isymp)),q(map),q(iqpos(isymp))
     1            ,nprm(isymp),nsymao(isymp),irmin,nsymp ,q(ig) )
      call tranki (q(iapos(isymp)),q(map),q(iqpos(isymp))
     1            ,nprm(isymp),nsymao(isymp),irmin,nsymp )
      map = map + (nsymp*(nsymp+1))/2
c
      else
c....  isymir.ne.1 code
      irmax=nsymp
      if (isymp.eq.isymr) irmax=nsymp+1-ir
      irmin=0
      if (isymp.eq.isymr) irmin=ir-1
      call mxmb (q(map+irmin*nsymq),nsymq,1,
     1          q(iqpos(isymq)),1,nsymq, q(iapos(isymp)+irmin),1,nsymp,
     2           irmax,nsymq,nprm(isymq) )
      call mxmb (q(map+irmin*nsymq),1,nsymq,
     1           q(iqpos(isymp)+irmin),1,nsymp,
     2           q(iapos(isymq)),1,nsymq, nsymq,irmax,nprm(isymp))
      map=map+nsymp*nsymq
c
c
      end if
120   continue
c..   now assemble and output c as contributions to d
      do 170 isymp=1,nirrr
      isymj=mults(isymp,isymir)
      isymij = mults(isymi,isymj)
      nprimj=nprm(isymj)
      nsymp=nsymao(isymp)
      do 160 j=1,nprimj
      ia1 = (j-1)*nsymp-1 + iapos(isymp)
      jj=iorbsm(istart(isymj)-1+j)
      if (ii.gt.jj) goto 140
c..   i.le.j  output as jipr
      do 130 ip=1,nsymp
      if (q(ia1+ip).eq.0.0d0) goto 130
      mword=mword+1
      gout(mword)=q(ia1+ip)
c     iout(mword) = (jj*(jj-1))/2+ii
c    1          +(ioffpr(isymij,isymp)+(ir-1)*nsymp+ip)*65536
      int4 = int4 + 2
      ijout(int4  ) = ioffpr(isymij,isymp)+(ir-1)*nsymp+ip
      ijout(int4+1) = (jj*(jj-1))/2+ii
      if (mword.eq.num2e) then
       call trout(gout,ijout,mword,icode,num4)
       int4 = -1
      endif
130   mword=mod(mword,num2e)
140   if (jj.gt.ii) goto 160
c...  j.le.i  output as ijrp
      do 150 ip=1,nsymp
      if (q(ia1+ip).eq.0.0d0) goto 150
      mword=mword+1
      gout(mword)=q(ia1+ip)
c     iout(mword) = (ii*(ii-1))/2+jj
c    1          +(ioffpr(isymij,isymr)+(ip-1)*nsymr+ir)*65536
      int4 = int4 + 2
      ijout(int4  ) = ioffpr(isymij,isymr)+(ip-1)*nsymr+ir
      ijout(int4+1) = (ii*(ii-1))/2+jj
      if (mword.eq.num2e) then
       call trout(gout,ijout,mword,icode,num4)
       int4 = -1
      endif
150   mword=mod(mword,num2e)
160   continue
170   continue
180   mapp = mapp + lentri
190   mloww=mhigh
      if (mword.ge.1) call trout(gout,ijout,mword,icode,num4)
      if (ipass.eq.npass) call put (gout,0,num4)
      isym12=isymir
      isym1=isymr
      itran1=i
      itran2=ir
      call corlsr(iapos(1))
      call accnt(' ',2)
      return
      end
      subroutine tran3(q,fc)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/prnprn)
      common /mctrns/ master,isym12,isym1,itran1,itran2,ntran
     >               ,ires,maxt,ibl5,ibl54,ibl56,nbuck,nteff
     >               ,mloww,mhi,ntri,nwb,ibuck
     >               ,mark(200),nwbuck(200),ioffpr(8,8)
      common /stak/ btri,mlow,nstack,iblock
      common /bufb/ nwbnwb,lnklnk,gsort(5118)
      common /junk/klsort(2,3412)
INCLUDE(common/jobopt)
      common /intbuf/ intpos,intfil,intmod
INCLUDE(common/mcff)
      dimension q(*),fc(*)
      dimension iqpos(8),ijpos(8),ikpos(8),i1pos(8)
      call accnt('tran3',2)
      ibase = icorr(0)
      intpos=-1
      intfil=num6
ctest iblf  =iblk6
c
c...  alloc storage for j,k transformation
      do 10 isym=1,nirrr
      junk=nsymao(isym)**2
      iqpos(isym)=icorr(junk)
      i1pos(isym)=icorr(junk)
      if(junk.gt.0)call vclr(q(i1pos(isym)),1,junk)
      ijpos(isym)=icorr(lensq)
10    ikpos(isym)=icorr(lensq)
      irpos=icorr(lensq)
      ig=icorr(lensq)
      call qget(q(1),iqpos)
      do 15 i=1,nirrr
      ia=iqpos(i)+ifreez(i)*nsymao(i)
      ja=iqpos(i)
      loop=nsymm(i)*nsymao(i)
      if(ia.gt.ja.and.loop.gt.0)
     * call fmove(q(ia),q(ja),loop)
 15   continue
c
_IF(parallel)
c **** MPP
      call closbf(0)
      call setbfa(-1)
c **** MPP
_ENDIF
      ibuck=0
      map=0
      itri=9999999
      mhi=(nprim*(nprim+1))/2
      mloww=-ntri+1
      lentr3=lentri*3
      nsize=lentr3*ntri
      ibuff=icorr(-nsize)
      do 180 i=1,nprim
      do 180 j=1,i
      itri=itri+1
      map=map+lentr3
      if (itri.le.ntri) goto 50
c..   read in next bucket
      mloww=mloww+ntri
      ibuck=ibuck+1
      call vclr(q(ibuff),1,nsize)
      mk=mark(ibuck)
20    if (mk.eq.9999999) goto 40
      call mcsrti(mk,nk)
      do 30 iword=1,nk
      map=ibuff-1+(klsort(1,iword)-mloww)*lentr3+klsort(2,iword)
      q(map)=q(map)+gsort(iword)
30    continue
      goto 20
40    map=ibuff
      itri=1
50    continue
      isymi=itype(i)
      isymj=itype(j)
      isymij=mults(isymi,isymj)
      mapp=map
      do 90 isymp=1,nirrr
      isymr=mults(isymp,isymij)
      nsymr=nsymm(isymr)
      nsymp=nsymm(isymp)
      nar=nsymao(isymr)
      nap=nsymao(isymp)
c...  transform exchange operator kij
      call mxmaa(q(map+ioffpr(isymij,isymp)),1,nap,
     1           q(iqpos(isymr)),1,nar,
     2           q(irpos),1,nap ,nap,nar,nsymr)
      call mxmaa(q(irpos),nap,1,
     1           q(iqpos(isymp)),1,nap,
     2           q(ikpos(isymr)),1,nsymr, nsymr,nap,nsymp)
c
c...  transform coulomb operator jij
      if (isymr-isymp) 60,70,90
60    call mxmaa(q(mapp),nar,1,
     1           q(iqpos(isymr)),1,nar,
     2           q(irpos),1,nap,  nap,nar,nsymr)
      mapp=mapp+nap*nar
      goto 80
70    continue
      call vclr(q(irpos),1,lensq)
      call trankh (q(irpos),q(mapp),q(iqpos(isymp))
     1            ,nsymm(isymp),nap,1,nap ,q(ig) )
      mapp=mapp+(nap*(nap+1))/2
80    call mxmaa(q(irpos),nap,1,
     1           q(iqpos(isymp)),1,nap,
     2           q(ijpos(isymp)),nsymp,1,  nsymr,nap,nsymp)
90    continue
      if (isymij.eq.1) then
c...  symmetrise j operator
      do 110 isyma=1,nirrr
      nsyma=nsymm(isyma)
      iab=ijpos(isyma)-nsyma-1
      do 100 ia=1,nsyma
      do 100 ib=1,ia
      q(iab+ib*nsyma+ia)=q(iab+ib*nsyma+ia)+q(iab+ia*nsyma+ib)
100   q(iab+ia*nsyma+ib)=q(iab+ib*nsyma+ia)
110   continue
      end if
c
      if (j.le.ncoremc .and. i.eq.j) then
c..   addition to core hamiltonian 2j(ii)-k(ii)
      do 120 isyma=1,nirrr
      loop=nsymm(isyma)*nsymm(isyma)
      if(loop.eq.0) go to 120
      call daxpy(loop
     + ,2.0d0,q(ijpos(isyma)),1,q(i1pos(isyma)),1)
      call daxpy(loop
     + ,-1.0d0,q(ikpos(isyma)),1,q(i1pos(isyma)),1)
120   continue
c
      end if
c...   put j&k matrices
      if(odebug(34)) call prisq(q(1),ijpos,isymij,1,'j-oper',10*i+j)
      if(odebug(34)) call prisq(q(1),ikpos,isymij,0,'k-oper',10*i+j)
      do 160 isyma=1,nirrr
      isymb=mults(isyma,isymij)
      if (isyma-isymb) 160,130,150
c...  j matrix, triangle
130   do 140 ia=1,nsymm(isyma)
140   call intou (q(ijpos(isyma)+(ia-1)*nsymm(isyma)),ia)
      goto 160
c...  j matrix, rectangle
150   call intou (q(ijpos(isyma)),nsymm(isyma)*nsymm(isymb))
160   continue
c
c...  k matrix
      do 170 isyma=1,nirrr
      isymb=mults(isyma,isymij)
170   call intou (q(ikpos(isyma)),nsymm(isyma)*nsymm(isymb) )
180   continue
c =========================================================
c
c...  do frozen core hamiltonian
      core = potnuc+efreez
      iad=1
      do 220 isym=1,nirrr
      n = nsymm(isym)
      na = nsymao(isym)
      if(n.le.0. or. na.le.0) go to 2202
      call vclr(q(irpos),1,na*n)
      call trankh (q(irpos),fc(iad),q(iqpos(isym)),n,na,
     1    1,na,q(ig))
      call mxmaa(q(irpos),na,1, q(iqpos(isym)),1,na,
     1q(ijpos(isym)),1,n, n,na,n)
      iab=ijpos(isym)-n-1
      do 200 ia=1,n
      do 200 ib=1,ia
      q(iab+ia*n+ib)=q(iab+ia*n+ib)+q(iab+ib*n+ia)
200   q(iab+ib*n+ia)=q(iab+ia*n+ib)
       if(odebug(34)) call prisq(q(1),ijpos,1,1,'h0',0)
       if(odebug(34)) call prisq(q(1),i1pos,1,1,'2jc-kc',0)
c
c...  now do core energy
2202  if(ncor(isym).le.0) go to 2201
      core = core + 2.0d0*dsum(ncor(isym),q(ijpos(isym)),n+1)
     1            +       dsum(ncor(isym),q(i1pos(isym)),n+1)
c...  so now it is ok to add in two electron
2201  if(n.le.0) go to 220
      call vadd(q(ijpos(isym)),1,q(i1pos(isym)),1,q(i1pos(isym)),1,n*n)
      do 210 ia=1,n
210   call intou (q(i1pos(isym)+(ia-1)*n),ia)
220   iad=iad+na*(na+1)/2
c.....dump one electron hamiltonian in mo basis
      if(odebug(34))  call prisq(q(1),i1pos,1,1,'fc',0)
      ij=0
      do 250 isym=1,nirrr
      ijj=ijpos(isym)-1
      do 240 i=1,nsymm(isym)
      do 230 j=1,i
230   fc(ij+j)=q(ijj+j)
      ij=ij+i
240   ijj=ijj+nsymm(isym)
250   continue
      call intou(fc,ltrimo)
c      print*,'core energy=',core
c
      call intend
      call corlsr(ibase)
      call accnt(' ',2)
      return
      end
_EXTRACT(trasor,mips4)
      subroutine trasor (g,nijkl,iunit,iblku)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/prnprn)
      common /mctrns/ master,isym12,isym1,itran1,itran2,ntran
     >               ,ires,maxt,ibl5,ibl54,ibl56,nbuck,nteff
     >               ,mloww,mhi,ntri,nwb,ibuck
     >               ,mark(200),nwbuck(200),ioffpr(8,8)
      common /stak / btri,mlow,nstack,iblock
INCLUDE(common/jobopt)
INCLUDE(common/multic)
INCLUDE(common/atmblk)
c
      common /lsort/ gin(510),mword,mwor,ijin(680)
      dimension g(ibl5,ires),nijkl(ibl5+ibl5,ires)
      call accnt('transort',2)
c...   determine base and limit triangles for this pass
      mloww=master
      mhi = min(mloww+nteff,ntran)
      mtri=mhi-mloww
c...   determine minimum no. of bucks.
      nbuck=ires
10    ntri=(mtri-1)/nbuck+1
      if(ntri.gt.maxt)goto 20
      nbuck=nbuck-1
      if(nbuck)10,20,10
20    nbuck=nbuck+1
      ntri=(mtri-1)/nbuck+1
c...   ntri=max. no. of triangles controlled by 1 bucket
c...   nbuck=number of buckets
      do 30 ibuck=1,nbuck
      mark(ibuck)=9999999
30    nwbuck(ibuck)=0
_IF(parallel)
c **** MPP
      call closbf(0)
      call setbfa(-1)
c **** MPP
_ENDIF
      iblock=0
      mlow=mloww+1
c...   start loop over secondary mainfile blocks
      call search(iblku,iunit)
      call find(iunit)
40    call get(gin(1),m)
      if (m.eq.0) goto 60
      call find(iunit)
c...   process input block
      call upcktr (gin(num2ep+1),ijin)
      iword2 = 1
      do 50 iword=1,mword
c **** I dont think the liitlendian stuff is now required
      ktri = ijin(iword2)
      itri = ijin(iword2+1)
      if (odebug(34)) then
      write(6,9998) iword, gin(iword), itri, ktri
 9998 format('trasor: iword, gin, itri, ktri = ', i4, f15.8, 2x, 2i8)
      endif
      iword2 = iword2 + 2
      if(itri.gt.mhi.or.itri.lt.mlow)goto 50
      ibuck=(itri-mlow)/ntri+1
      nwb=nwbuck(ibuck)+1
      g(nwb,ibuck)=gin(iword)
      nijkl(nwb+nwb-1,ibuck)=itri
      nijkl(nwb+nwb  ,ibuck)=ktri
      if(nwb.eq.ibl5) call mcsrto (g(1,ibuck),nijkl(1,ibuck))
      nwbuck(ibuck)=nwb
50    continue
      goto 40
60    continue
c...   secondary mainfile now swept
c...   clear up output
      do 70 ibuck=1,nbuck
      nwb=nwbuck(ibuck)
      if (nwb.ne.0) call mcsrto (g(1,ibuck),nijkl(1,ibuck))
70    continue
      call stopbk
      call accnt(' ',2)
      return
      end
_ENDEXTRACT
      subroutine mcupda (q,sol,iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
      dimension q(*),sol(*),ixpos(8),iqpos(8)
      ibase = icorr(0)
      do 10 i=1,nirrr
      iqpos(i) = icorr(nsymao(i)**2)
10    ixpos(i) = icorr(nsymao(i)**2)
      iu = icorr(maxbas**2)
      iv = icorr(maxbas**2)
c     iw = icorr(maxbas**2)
      call qget (q(1),iqpos)
      call xasm (q(1),ixpos,sol)
      do 40 isym=1,nirrr
      ns = nsymm(isym)
      ns2 = ns**2
      if(ns.le.0) go to 400
c...  u = exp(x)
c...  estimate no. of terms wanted in exponential
      zz = dasum(ns2,q(ixpos(isym)),1)
      fac = zz
      nterm = 2
20    nterm = nterm + 1
      fac = fac*zz/ dfloat(nterm)
      if (fac.gt.1.0d-4) goto 20
c...  so now expand exponential (from the inside out)
      call vclr(q(iu),1,ns2)
c      print 30,nterm
c30    format(' expand orbital rotation with',i5,'terms')
      do 30 i=nterm,1,-1
      call vclr(q(iv),1,ns2)
      call daxpy(ns2,(1.0d0/ dfloat(i)),q(iu),1,q(iv),1)
      call vclr(q(iu),1,ns2)
      call vfill(1.0d0,q(iu),ns+1,ns)
30    call mxmb (q(ixpos(isym)),1,ns, q(iv),1,ns, q(iu),1,ns,ns,ns,ns)
c
      if (lto(4)) call outsqr (q(iqpos(isym)),ns,ns,ns,'old mos')
      if (lto(4)) call outsqr (q(iu),ns,ns,ns,'transformation matrix')
400   loop =         ifreez(isym)*nsymao(isym)
      ia = iqpos(isym) + loop
      ja = ixpos(isym) + loop
      if(loop.gt.0)
     *call fmove(q(iqpos(isym)),q(ixpos(isym)),loop)
      if(ns.ne.0.and.nsymao(isym).ne.0)
     *call mxmaa(q(ia),1,nsymao(isym), q(iu),1,ns,
     1 q(ja),1,nsymao(isym),nsymao(isym),ns,ns)
40    continue
c
c..   orthogonalise the new mos and write to dump
      call qput (q(1),ixpos,iwrite)
      call corlsr (ibase)
c
c...   deal with ci update if there is one.
      if (isigma.ne.1) return
      ibuff = icorr(nci)
      call cget (q(ibuff),1)
      zz = -ddot(nci,sol(nrot+1),1,q(ibuff),1)
      call daxpy(nci,zz,q(ibuff),1,sol(nrot+1),1)
      zl = dnrm2(nci,sol(nrot+1),1)
      cosl =  dcos(zl)
      if ( dabs(zl).gt.1.0d-14) zl =  dsin(zl)/zl
      call dscal(nci,cosl,q(ibuff),1)
      call daxpy(nci,zl,sol(nrot+1),1,q(ibuff),1)
      zz = 1.0d0/dnrm2(nci,q(ibuff),1)
      call dscal(nci,zz,q(ibuff),1)
      call cput (q(ibuff),1)
      call corlsr (ibuff)
      return
      end
c ******************************************************
c ******************************************************
c             =   mcanal     =
c ******************************************************
c ******************************************************
      subroutine mcanal(q,iq,iwrite,ipunch)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mcscra)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
INCLUDE(common/prints)
      dimension q(*),iq(*)
c
      if (nstate.gt.1) call mceprn(q(1),iwrite)
      if(nstate.eq.1) write(iwrite,10)enext
10    format(/' total energy',t22,f22.12, ' a.u.')
c
      if (n1elec.gt.0) call mcprop(q(1),iwrite)
      if (iguess.le.0) goto 20
c
      iuu=icorr(nprim**2)
      call canorb(q(1),iq(1),q(iuu),icani,icant,icana,isecn,1,iwrite,
     +            ipunch)
      call corlsr(iuu)
c
20    call mcprnt(q(1),lprint,iwrite,ipunch)
c
c      call mcpun(q(1))
c
      return
      end
      subroutine mcprorb(q)
      implicit REAL (a-h,o-z)
      dimension q(*)
c
INCLUDE(common/sizes)
      parameter (mxorb1 = maxorb+1)
INCLUDE(common/iofile)
INCLUDE(common/restri)
INCLUDE(common/prints)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
      common/blkorbs/eig(maxorb),pop(mxorb1),nba,new,ncoll,
     * jeig,jocc,ipad
      character*8 zcomm
      common/junkc/zcomm(29)
INCLUDE(common/machin)
INCLUDE(common/atmol3)
INCLUDE(common/segm)
INCLUDE(common/mapper)
INCLUDE(common/multic)
c      common/multic/radius(41),irad(25),ncoremc,nact,jrad(5),
c     +              nfreez,krad(13),itype(maxorb),lrad(75),isecn
INCLUDE(common/harmon)
INCLUDE(common/natorb)
c
      data m3,m29,m10,m16/3,29,10,16/
c
      l0 = newbas0
      l1 = num
      l3 = num*num
c     decide on high or low format for vector print
      if (oprint(20)) then
       lprnt = l1
      else
       lprnt = nfreez+ncoremc+nact
       lprnt = min(lprnt+5,l0)
      endif
c
      nwf = mach(8)
      nwctr = mach(9)
c
c     ----- set pointers for partitioning of core -----
c
      length = max(10,l3)
      i10 = icorr(length)
c
      mout = mouta
      if (moutb.gt.0) mout = moutb
      if (isecn.gt.0) mout = isecn
      write(iwr,6061) mout
6061  format(/,' vectors restored from section',i4)
      call secget(mout,m3,iblvec)
      call rdchr(zcomm(1),m29,iblvec,idaf)
      call reads(eig,nwf,idaf)
      iblvec = iblvec + lensec(nwctr) + lensec(nwf) + lensec(m29)
      call rdedx(q(i10),l3,iblvec,idaf)
      call tdown(q(i10),ilifq,q(i10),ilifq,l0)
      write (iwr,6060)
 6060 format (/10x,42('=')/10x,'mcscf ''natural'' orbitals',
     1        ' in standard basis',/10x,42('='))
      call prev(q(i10),pop,lprnt,l1,l1)
      call corlsr(i10)
c
      return
      end
      subroutine blkop(q,isym,i1,m,orb,n1,opr,n2,u,n3,iwrite)
c.....orbitals (i1+i,i=1,m) in symmetry isym are reordered according to
c.....local symmetry.
c.....orb, opr only symmetry block isym given
c.....(i1 does not incude frozen cores, but frozen core orbitals are
c.....assumed to be in orb)
c.....n1=first dimension in orb
c.....(if n1=0 orbitals and iosym not reordered)
c.....n2=first dimension in opr
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      logical ord
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
      common /block / nosymm,nosym(8),iosym(8,maxorb)
      dimension q(*),iof(20),lblk(20)
      dimension opr(*),orb(*),u(n3,n3)
c
      itemp=icorr(n2**2)
      call fmove(opr,q(itemp),n2**2)
      ls=nosym(isym)
      if(ls.le.1) return
      ord=n1.gt.0
      do 10 is=1,ls
10    lblk(is)=0
      do 20 i=1,m
      is=iosym(isym,i1+i)
20    lblk(is)=lblk(is)+1
      do 30 is=1,ls
30    if(lblk(is).eq.m) return
      if(n3.ne.0) then
      do 35 i=1,m
      do 35 j=1,m
35    u(i1+i,i1+j)=0.0d0
      end if
c      if(ord) call outsqr(orb,n1,n1,n1,'initial orbitals')
c      call outsqr(opr,n2,n2,n2,'initial operator')
c      print 40,isym,(iosym(isym,i),i=1,n2)
c40    format(/' initial orbital symmetries in siym',i2,' :',(t40,40i2))
      iof(1)=0
      do 60 is=2,ls
60    iof(is)=iof(is-1)+lblk(is-1)
c.....reorder orbitals and columns of operator
      j1=icorr(n1*m)
      j2=icorr(n2**2)
      k1=i1+ifreez(isym)-1
      k2=i1-1
      do 70 i=1,m
      is=iosym(isym,i1+i)
      if(ord) call fmove(orb((k1+i)*n1+1),q(j1+iof(is)*n1),n1)
      call fmove(opr((k2+i)*n2+1),q(j2+iof(is)*n2),n2)
      iof(is)=iof(is)+1
70    if(n3.ne.0) u(i1+i,i1+iof(is))=1.0d0
c      if(n3.gt.0) call outsqr(u,n3,n3,n3,'transformation matrix')
      if(ord) call fmove(q(j1),orb((k1+1)*n1+1),n1*m)
      call fmove(q(j2),opr(i1*n2+1),m*n2)
      call fmove(opr,q(j2),n2*n2)
      call corlsr(j1)
c.....reorder rows of operator
      do 80 is=1,ls
80    iof(is)=iof(is)-lblk(is)+i1-n2
      i2=j2+i1-n2-1
      do 90 i=1,m
      is=iosym(isym,i1+i)
      iof(is)=iof(is)+1
      i2=i2+1
      do 90 j=1,n2
90    opr(iof(is)+j*n2)=q(i2+j*n2)
      call corlsr(j2)
c      if(ord) call outsqr(orb,n1,n1,n1,'final orbitals')
c      call outsqr(opr,n2,n2,n2,'final operator')
c.....check blocking
      nfault=0
      i=i1
      do 120 is=1,nosym(isym)
      do 120 ii=1,lblk(is)
      j1=i*n2+i1
      i=i+1
      if(ord) iosym(isym,i)=is
      j=i1
      do 110 js=1,nosym(isym)
      if(is.eq.js) then
      j=j+lblk(js)
      goto 110
      end if
      do 100 jj=1,lblk(js)
      j=j+1
      if( dabs(opr(j1+jj)).gt.1.d-6) nfault=nfault+1
c      if( abs(opr(j1+jj)).gt.1.e-6) print 101,i,j,is,js,ii,jj,n2
c     1   ,j1+jj,opr(j1+jj)
c101  format(' symm fault: i,j,is,js,ii,jj,n2,ij,val:',8i4,e12.4)
100   opr(j1+jj)=0
110   j1=j1+lblk(js)
120   continue
c      if(ord) print 130,isym,(iosym(isym,i),i=1,n2)
c130  format(/' final orbital symmetries in sym',i2,' : ',(t40,40i2))
      if(nfault.ne.0) then
      write(iwrite,140) nfault,isym
140   format(/' *** warning (blkop):',i5,
     1   ' non zero elements in sym',i3)
c      call outsqr(q(itemp),n2,n2,n2,'initial operator')
c      call outsqr(opr,n2,n2,n2,'final operator')
      end if
      call corlsr(itemp)
      return
      end
      subroutine canorb(q,iq,u,itycor,ityact,ityext,istor,iprnt,
     +                  iwrite,ipunch)
      implicit REAL  (a-h,o-z)
      logical fock
      character*7 typ
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
INCLUDE(common/mcscra)
      common /block / nosymm,nosym(8),iosym(8,maxorb)
      common /intbuf/ intpos,intfil,intmod
INCLUDE(common/mcff)
      dimension q(*),iq(*),iqpos(8),ijpos(8),ikpos(8),u(*)
c
      if(itycor+ityact+ityext.eq.0) return
      ibase=icorr(0)
      fock=itycor.eq.1.or.ityact.eq.1.or.ityext.eq.1
      if(fock.and.itrsfm.eq.0) call mctran(q(1),iq(1))
      if(iguess.eq.0) call caserr('no ci vector provided')
      call ofset(1)
      call vclr(u,1,nlu)
      do 40 is=1,nirrr
      if(nprm(is).gt.0) then
      call vfill(1.0d0,u(ioffu(is)+1),nprm(is)+1,nprm(is))
      endif
 40   continue
c
c.....first order density matrix
c
      iu=icorr(nlu)
      id1=icorr(nact**2)
      call vclr(q(id1),1,nact**2)
      icivec=icorr(nci*nstate)
      call cget(q(icivec),nstate)
      igam=icorr(ne)
      ivec=icivec-nci
c     nbl=lensec(nci)
      do 10 istate=1,nstate
      call denst1(q(1),q(ivec+istate*nci),q(ivec+istate*nci),q(igam))
      call daxpy(nact**2
     +  ,weight(istate),q(igam+ic1d),1,q(id1),1)
10    continue
      call corlsr(igam)
      if(icinat.le.0) call corlsr(icivec)
c
c.....fock operator if required
c
      ieig=icorr(nbasao)
      call vclr(q(ieig),1,nbasao)
      igc=icorr(ntqg)
      call vclr(q(igc),1,ntqg)
      if(fock) then
      do 20 i=1,nirrr
20    ijpos(i)=icorr(nsymm(i)**2)
      do 30 i=1,nirrr
30    ikpos(i)=icorr(nsymm(i)**2)
      intpos = 0
      intfil = num6
      iblf   = iblk6
      do 50 i=1,nprim
      is=itype(i)
      do 50 j=1,i
      call loadop(q(1),ijpos,ikpos,i,j)
      if(itype(j).ne.is) goto 50
      if(i.le.ncoremc.or.j.le.ncoremc) goto 50
      ij=id1+ilifa(i-ncoremc)+j-ncoremc-1
      if(q(ij).eq.0) goto 50
      fac=q(ij)
      if(i.ne.j) fac=fac+fac
      call daxpy(ntqg,fac,q(ijpos(1)),1,q(igc),1)
      call daxpy(ntqg,-0.5d0*q(ij),q(ikpos(1)),1,q(igc),1)
      if(i.eq.j) goto 50
      call transp(q(ikpos(1)),q(ikpos(1)),1,q(ijpos(1)))
      call daxpy(ntqg,-0.5d0*q(ij),q(ikpos(1)),1,q(igc),1)
50    continue
      call loadop(q(1),ijpos,ikpos,0,0)
      call vadd(q(igc),1,q(ijpos(1)),1,q(igc),1,ntqg)
_IF1()c      call priop(q(igc),1,0,'gc-matrix',0,iwrite)
      end if
      ien=ieig-1
      do 56 isym=1,nirrr
      ien=ien+ifreez(isym)
      n=nsymm(isym)
      n1=n+1
      igd=igc+ioffq(1,isym)-n1
      do 55 i=1,n
55    q(ien+i)=q(igd+i*n1)
56    ien=ien+n
      call corlsr(ijpos(1))
c
c.....pseudo canonical core orbitals
c
      do 60 i=1,nirrr
60    iqpos(i)=icorr(nsymao(i)**2)
      iocc=icorr(nbasao)
      iscr=icorr(maxbas)
      ibuf=icorr(maxbas**2)
      call qget(q(1),iqpos)
      call vclr(q(iu),1,nlu)
      iuu=iu
      do 80 i=1,nirrr
      if(nprm(i).le.0) go to 80
      call vfill(1.0d0,q(iuu),nprm(i)+1,nprm(i))
       iuu=iuu+nprm(i)**2
 80   continue
      if(itycor.ne.0) then
      ien=ieig
      do 100 isym=1,nirrr
      ig=ioffq(1,isym)+igc
      n=nsymm(isym)
      m=ncor(isym)
      ien=ien+ifreez(isym)
      if(m.le.1) goto 100
      iuu=ioffu(isym)+1
      n1=nsymao(isym)
      call blkop(q(1),isym,0,m,q(iqpos(isym)),n1,q(ig),n,u(iuu),m,
     +           iwrite)
      iv=iu+ioffu(isym)
      ifail=0
      call f02abf(q(ig),n,m,q(ien),q(iv),nprm(isym),q(iscr),ifail)
      n=nsymao(isym)
      ii=iqpos(isym)+ifreez(isym)*n
      call fmove(q(ii),q(ibuf),n*m)
      i=nprm(isym)
      call mxmaa(q(ibuf),1,n,q(iv),1,i,q(ii),1,n,n,m,m)
      call mxmaa(u(iuu),1,i,q(iv),1,i,q(ibuf),1,i,i,i,i)
      call fmove(q(ibuf),u(iuu),i**2)
100   ien=ien+nsymm(isym)
      end if
c
c.....pseudo canonical external orbitals
c
      if(ityext.eq.1) then
      ivec=icorr(maxbas**2)
      ien=ieig
      do 200 isym=1,nirrr
      n=nsymm(isym)
      m=n-nprm(isym)
      ien=ien+ifreez(isym)+nprm(isym)
      if(m.le.1) goto 200
      n1=nsymao(isym)
      ig=ioffq(1,isym)+igc
      call blkop(q(1),isym,nprm(isym),m,q(iqpos(isym)),n1,q(ig),n,
     * q(ivec),0,iwrite)
      ig=igc+ioffq(1,isym)+nprm(isym)*(n+1)
      ifail=0
      call f02abf(q(ig),n,m,q(ien),q(ivec),m,q(iscr),ifail)
      n=nsymao(isym)
      ii=iqpos(isym)+(ifreez(isym)+nprm(isym))*n
      call fmove(q(ii),q(ibuf),n*m)
      call mxmaa(q(ibuf),1,n,q(ivec),1,m,q(ii),1,n,n,m,m)
200   ien=ien+m
      call corlsr(ivec)
      end if
c
c.....natural or pseudo canonical active orbitals
c
      idc=icorr(nlu)
      id=id1-1
      call vclr(q(idc),1,nlu)
      do 250 i=1,nact
      is=itypea(i)
      ii=idc+ioffu(is)+(iorb(i+ncoremc)-1)*nprm(is)-1
      do 240 j=1,nact
      if(itypea(j).ne.is) goto 240
      q(ii+iorb(j+ncoremc))=-q(id+ilifa(i)+j)
240   continue
250   continue
      ioc=iocc
      ien=ieig
      do 300 isym=1,nirrr
      n=nsymm(isym)
      m=nactt(isym)
      ien=ien+ifreez(isym)+ncor(isym)
      if(m.eq.0) goto 300
      i=nprm(isym)
      n1=nsymao(isym)
      iuu=ioffu(isym)+1
      ig=ioffq(1,isym)+igc
      id=ioffu(isym)+idc
      if(ityact.ne.0) then
      if(fock) call blkop(q(1),isym,ncor(isym),m,q(iqpos(isym)),0,q(ig),
     1n, u(iuu),0,iwrite)
      call blkop(q(1),isym,ncor(isym),m,q(iqpos(isym)),n1,q(id),i,
     * u(iuu),i,iwrite)
      end if
      ia=ioffq(1,isym)+ncor(isym)*(n+1)
      ib=ioffu(isym)+ncor(isym)*(i+1)
      iv=iu+ib
      i1=i+1
      i2=n+1
      k1=ioc-1
      k2=ien-1
      ii=idc+ib-i1
      jj=igc+ia-i2
c.....diagonal elements of density and fock operator
      do 285 j=1,m
      q(k2+j)=q(jj+j*i2)
285   q(k1+j)=q(ii+j*i1)
      if(ityact.eq.1) then
c.....diagonalize fock operator
      ifail=0
      call f02abf(q(igc+ia),n,m,q(ien),q(iv),i,q(iscr),ifail)
      call mxmaa(q(idc+ib),1,i,q(iv),1,i,q(ibuf),1,m,m,m,m)
      j1=ioc-1
      end if
      if(ityact.eq.2) then
c.....diagonalize density matrix
      ifail=0
      call f02abf(q(idc+ib),i,m,q(ioc),q(iv),i,q(iscr),ifail)
      call mxmaa(q(igc+ia),1,n,q(iv),1,i,q(ibuf),1,m,m,m,m)
      j1=ien-1
      end if
      ioc=ioc+m
      if(ityact.eq.0) goto 300
c.....transform diagonal elements of operator not used
      do 290 j=1,m
      q(j1+j)=ddot(m,q(ibuf+(j-1)*m),1,q(iv+(j-1)*i),1)
290   continue
c.....transform orbitals and form total transformation matrix
      n=nsymao(isym)
      ii=iqpos(isym)+(ifreez(isym)+ncor(isym))*n
      call fmove(q(ii),q(ibuf),n*m)
      call mxmaa(q(ibuf),1,n,q(iv),1,i,q(ii),1,n,n,m,m)
      call mxmaa(u(iuu),1,i,q(iu+ioffu(isym)),1,i,q(ibuf),1,i,i,i,i)
      call fmove(q(ibuf),u(iuu),i*i)
300   ien=ien+nsymm(isym)-ncor(isym)
      if(iprnt.eq.0) goto 1000
c
c.....print transformed orbitals
c
      if (mcprin) write(iwrite,310)
310   format(//' natural or pseudo canonical orbitals'/
     1         ' ------------------------------------'//
     2 ' nr sym   typ          occ       energy      coefficients')
      ioc=iocc
      ien=ieig-1
      char=0.0d0
      do 330 isy=1,nirrr
      isym=isy
      n=nsymao(isym)
      m=min(n,ifreez(isym)+nprm(isym)+2)
      ij=iqpos(isym)-1
      do 320 i=1,m
      occ=2.0d0
      typ='frozen '
      if(i.gt.ifreez(isym)) typ='core   '
      if(i.gt.ifreez(isym)+ncor(isym)) typ='active '
      if(i.gt.ifreez(isym)+nprm(isym)) typ='virtual'
      if(typ.eq.'virtual') occ=0.0d0
      if(typ.eq.'active') occ=-q(ioc)
      char=char+occ
      if(typ.eq.'active') ioc=ioc+1
      qmfg=q(ien+i)
      if (mcprin) write(iwrite,340)
     *i,isym,typ,occ,qmfg,   (q(ij+j),j=1,n)
340   format(/1x,i2,i3,3x,a7,2f12.5,(t41,5f15.8))
320   ij=ij+n
      ien=ien+nsymao(isym)
330   continue
      if (mcprin) write(iwrite,350)char
350   format(/' total charge:',f17.12)
c
c.....dump transformed orbitals
c
      if(istor.gt.0) then
      isave=iblkq
      iblkq=iblkn
      call hedout(q(ieig),q(iocc))
      call qput(q(1),iqpos,iwrite)
      if(iblkn.gt.0.and.mcprin) write(iwrite,385)isecn
385   format(/' natural orbital dump at section',
     1   t46,i3)
      iblkq=isave
      end if
c
c.....ci if required
c
      if(icinat.gt.0) then
      if (mcprin) write(iwrite,390)
390   format(/' ci diagonalisation for new orbitals'/
     1        ' -----------------------------------'/)
      call ofset(0)
      izint=icorr(ne)
      ijjpos=icorr(lenj)
      ikkpos=icorr(lenk)
      call loadi(q(1),q(ijjpos))
      l=icorrm()
      iad=icorr(-l)
      call traop(q(iu),q(ijjpos),q(ikkpos),q(iad),l,luse,0,iwrite)
      call mkzint(q(ijjpos),q(izint))
      call corlsr(ijjpos)
      isg=icorr(nci*nstate)
      isave=iguess
      iguess=0
      call ci2(q(1),iq(1),q(icivec),q(isg),q(izint),q(ieig),0,10,1.0d-5,
     +         iwrite)
      iguess=isave
      call corlsr(izint)
      if(isecnc.ne.0) then
      isav1=isecd
      isav2=iblkdmc
      isav3=iblkc
      isav4=isec
      isav5=iblkq
      if(istor.gt.0) then
      isec=isecn
      iblkq=iblkn
      end if
      isecd=isecnc
      call cidmpi
c.....new dump section refers to natural orbitals and new ci vector
      call mcdump
      call cput(q(icivec),nstate)
      if (mcprin) write(iwrite,387)isecnc
387   format(/' job info and ci vector dump at section',t46,i3)
      isecd=isav1
      iblkdmc=isav2
      iblkc=isav3
      isec=isav4
      iblkq=isav5
      end if
      call prici(q(1),iq(1),icivec,iwrite,ipunch)
      if (mcprin) write(iwrite,360)(q(ieig+i-1)+core,i=1,nstate)
360   format(/' final state energies:  ',(t22,5f18.10))
      if (mcprin) write(iwrite,361)
361   format()
      end if
1000  call corlsr(ibase)
      return
      end
      subroutine prici(q,iq,ivec,iwrite,ipunch)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/detcic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
INCLUDE(common/prnprn)
      common /couple/ surd(511)
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1              ,jorb(511),korb(511),lorb(511)
      dimension q(*),iq(*),icga(31),icgb(31)
      character*256 name
      data cipunch /1.0d-10/
      nocc=nact*nci
      if (iexc.lt.0) nocc=2*nact+1
      iocc = icori(nocc)
      if (iexc.ge.0) then
       do i=1,nocc
        iq(i+iocc-1) = 0
       enddo
       call rdedx(surd,511,iblkft,numft)
       do i=1,nblkd1
       call find(numft)
       call get (g,nw)
       call unpkft(nw)
        do iw=1,nw
         io=iorb(iw)
         iii=ii(iw)
         iq(iocc-1+(io-1)*nci+iii) = surd(ival(iw))+.0001d0
        enddo
       enddo
c       do 1111 i=1,nci
c1111   write (6,*) i, (iq(iocc-1+(j-1)*nci+i),j=1,nact)
c
      endif
c
      idone = icori(nci)-1
      do 350 istate=1,nstate
      write(iwrite,300)istate
300   format(/' ci vector for state',i2/
     1       ' ---------------------'//
     2   '     nr          coefficient     occupancy'/)
      do i=1,nci
       iq(idone+i) = 0
      enddo
      mmax = min(nci,20)
      do 330 iii=1,mmax
      z = -1.0d0
      do 320 j=1,nci
      if ( dabs(q(ivec+(istate-1)*nci-1+j)).lt.z .or.
     1      iq(idone+j).ne.0) goto 320
      z =  dabs(q(ivec+(istate-1)*nci-1+j))
      i = j
320   continue
      iq (idone+i) = 1
      mfg = ivec + (istate-1)*nci -1 +i
      if (iexc.lt.0) then
       call detocc (i,icga,icgb)
       nact2=2*nact
       do j=1,nact2
        iq(iocc+j)=0
       enddo
       do j=1,na
        ia=2*icga(j)-1
        iq(iocc+ia)=1
       enddo
       do j=1,nb
        ib=2*icgb(j)
        iq(iocc+ib)=1
       enddo
       write(iwrite,335)
     1  i,q(mfg),(iq(iocc+j),j=1,nact2)
335    format(i7,f25.14,6x,31(2i1,1x))
      else
       write(iwrite,340)
     1  i,q(mfg),(iq(iocc-1+(io-1)*nci+i)
     2           ,io=1,nact)
      end if
330   continue
c     should we punch the coefficients?
      if(opunch(17)) then
       write(name,'(A13,I3.3)')'civecs.ascii_',istate
       open(ipunch,file=name,form='formatted',
     +      status='unknown')
       rewind ipunch
       do k=1,nci
        mfg = ivec + (istate-1)*nci -1 + k
        if(dabs(q(mfg)).ge.cipunch) then
         if (iexc.lt.0) then
          call detocc (k,icga,icgb)
          nact2=2*nact
          do j=1,nact2
           iq(iocc+j)=0
          enddo
          do j=1,na
           ia=2*icga(j)-1
           iq(iocc+ia)=1
          enddo
          do j=1,nb
           ib=2*icgb(j)
           iq(iocc+ib)=1
          enddo
          write(ipunch,335)k,q(mfg),
     +     (iq(iocc+j),j=1,nact2)
         else
          write(ipunch,340)k,q(mfg),
     +     (iq(iocc-1+(io-1)*nci+i),io=1,nact)
         end if
        endif
       enddo
       close(ipunch,status='keep')
      endif
340   format(i7,f25.14,5x,63(i2,1x))
350   continue
      call corlsi (iocc)
      return
      end
      subroutine mceprn(q,iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/syminf)
INCLUDE(common/multic)
INCLUDE(common/jobopt)
      dimension q(*),energ(5)
      ibase = icorr(0)
      ivec = icorr(nci*nstate)
      isigma = icorr(nci)
      idiag = icorr(nci)
      izint = icorr(ne)
c     lbl = (nci-1)/511+1
      call ci0 (q(1),q(izint),q(idiag))
      if(iguess.eq.0) call caserr('no ci vector provided')
      call cget(q(ivec),nstate)
      do 10 istate=1,nstate
      call vclr(q(isigma),1,nci)
      call ci1 (q(1),q(ivec),q(isigma),q(idiag),q(izint))
      energ(istate) = ddot(nci,q(ivec),1,q(isigma),1)+core
      ivec=ivec+nci
10    continue
      if(mcprin) write(iwrite,20)(energ(istate),istate=1,nstate)
20    format(/1x,'final state energies: ',5f20.12)
      call corlsr(ibase)
      return
      end
      subroutine mcprop(q,iwrite)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (minpr=1,maxpr=6)
      logical many
      character*22 string(2)
      character*20 chars(minpr:maxpr)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
      common/mceig/eigmc(5)
      dimension ppmat(3,5,5)
INCLUDE(common/multic)
      dimension q(*),dt(5),pmat(5,5),dmat(5,maxpr)
      data string/' ','(state averaged value)'/
      data debye/2.54158d0/
      data chars /'overlap'
     1           ,'kinetic energy'
     2           ,'one electron energy'
     3           ,'x dipole'
     4           ,'y dipole'
     5           ,'z dipole'
     6/
      call accnt('mcprop',1)
      call vclr(ppmat,1,3*5*5)
      ibase = icorr(0)
c     many = ssum(nstate-1,weight,1).gt.0.0e0
      many = nstate.gt.1
c     npri = nprim+nfreez
c     nco = ncoremc+nfreez
      npr2 = 0
      do 10 i=1,nirrr
10    npr2 = npr2 + (nactt(i)*(nactt(i)+1))/2
      iden = icorr(npr2 * (nstate*(nstate+1))/2 )
c
c...  form (transition) density matrices
      igam = icorr(ne)
      ivec = icorr(nci*nstate)
      if(iguess.eq.0) call caserr('no ci vector provided')
      call cget (q(ivec),nstate)
c      call outvec (q(ivec),nci*nstate,'result of cget')
      ii = iden
      do 50 istat1 = 1,nstate
      ivec1 = (istat1-1)*nci+ivec
      do 50 istat2 = 1,istat1
      ivec2 = (istat2-1)*nci+ivec
      call denst1 (q(1),q(ivec1),q(ivec2),q(igam))
c...  symmetry pack density matrix
      do 40 isym=1,nirrr
      do 30 i=1,nact
      if (itypea(i).ne.isym) goto 30
      do 20 j=1,i
      if (itypea(j).ne.isym) goto 20
      q(ii) = q(igam+ic1d+(i-1)+(j-1)*nact)
      if (i.ne.j) q(ii) = q(ii) + q(igam+ic1d+(j-1)+(i-1)*nact)
      ii = ii + 1
20    continue
30    continue
40    continue
c      call outsqr (q(igam+ic1d),nact,nact,nact,'raw density')
c      call outvec (q(ii-npr2),npr2,'symmetrised density')
50    continue
      call corlsr (igam)
c
c...  so now we can do properties
      iprop = icorr(npr2)
      call vclr(dt,1,nstate)
      call vclr(dmat,1,maxpr*5)
c     dtav = 0.0d0
      vt = -1.0d0
      do 160 i=1,n1elec
      ipr = i1elec(i)
      if (ipr.gt.maxpr .or. ipr.lt.minpr) then
      write(iwrite,70)ipr
70    format (/1x,'*** requested property number',i8,' is out of range')
      goto 160
      end if
      call pget (q(1),q(iprop),ipr,znuc,zcor)
c      call outvec (q(iprop),npr2,'property integrals')
      ii = iden
      sign=1.0d0
      if(ipr.ge.4.and.ipr.le.6) sign=-1.0d0
      ityp=0
      nonz=0
      do 90 istate=1,nstate
      if(weight(istate).ne.0) ityp=ityp+1
      do 80 jstate=1,istate
      pmat(jstate,istate) = sign*ddot(npr2,q(ii),1,q(iprop),1)
      pmat(istate,jstate) = pmat(jstate,istate)
      if(istate.eq.jstate) goto 80
      if(pmat(istate,jstate).ne.0) nonz=nonz+1
80    ii = ii + npr2
      pmat(istate,istate) = pmat(istate,istate) + znuc + zcor*sign
      if(pmat(istate,istate).ne.0) nonz=nonz+1
90    dmat(istate,ipr) = pmat(istate,istate)
c
c...  dipole
c
      if (ipr.ge.4.and.ipr.le.6) then
      do 100 istate=1,nstate
100   dt(istate)=dt(istate)+dmat(istate,ipr)**2
      end if
c..   kinetic energy
      pav=ddot(nstate,dmat(1,ipr),1,weight,1)
      if (ipr.eq.2. and .pav.ne.0.0d0) vt = (pav-enext)/pav
      if (many.and.nonz.ne.0) then
      write(iwrite,120)chars(ipr)
120   format(/' transition matrix for the property  ',a20/)
      do 130 istate=1,nstate
130   write(iwrite,140)(pmat(istate,j),j=1,nstate)
c
      if (ipr.ge.4.and.ipr.le.6) then
         if (nstate.gt.5) call caserr('extend mcprop dimensions')
c...   pmat**2*(delta-e)*2/3
c...   pmat**2*(delta-e)*2/3
         do istate=1,nstate
            do j=1,istate-1
               ppmat(ipr-3,istate,j) = (pmat(istate,j)**2)*
     1                          abs(eigmc(istate)-eigmc(j))/1.5d0
            end do
         end do
      end if
c
140   format(1x,5f16.6)
      end if
160   continue
c
      write(iwrite,122)
122   format(/,1x,88(1h-),/,1x,
     1 ' Excitation energies and  oscillator strengths',/,13x,
     2 'Delta E   Delta E(eV)',6x,'Fx',12x,'Fy',12x,'Fz',12x,'F')
141   format(2i4,f12.6,f10.2,4f14.6)
      do istate=1,nstate
         do j=1,istate - 1
            write(iwrite,141) j,istate,eigmc(istate)-eigmc(j),
     1            (eigmc(istate)-eigmc(j))*27.21165d0,
     2            ppmat(1,istate,j),ppmat(2,istate,j),ppmat(3,istate,j),
     3            ppmat(1,istate,j)+ppmat(2,istate,j)+ppmat(3,istate,j)
         end do
      end do
      write(iwrite,'(1x,88(1h-))')
c
      do 170 istate=1,nstate
      dt(istate) =  dsqrt(dt(istate))
      istat=iroot1+istate-1
      if (mcprin) write(iwrite,60)istat
60    format(/' one electron properties for state',i2/
     1        ' -----------------------------------'/)
      do 165 i=1,n1elec
      ipr=i1elec(i)
      if(ipr.lt.4.or.ipr.gt.6) then
      if (mcprin) write(iwrite,115)chars(ipr),dmat(istate,ipr)
      else
      if(dabs(dabs(dmat(istate,ipr))-dt(istate)).lt.1.0d-8)dt(istate)=0
      if (mcprin) write(iwrite,115)
     *chars(ipr),dmat(istate,ipr),dmat(istate,ipr)*debye
115   format(1x,a20,f16.6,' a.u.',:,f16.6,' debye',5x,a22)
      end if
165   continue
      if(dt(istate).ne.0.and.mcprin) write(iwrite,115)
     +       'total dipole moment ',
     +        dt(istate),dt(istate)*debye
170   continue
      if (vt.gt.0.0d0.and.mcprin) write(iwrite,190)vt,string(ityp)
190   format(/' virial ratio',f24.6,5x,a22)
      call corlsr (ibase)
      call accnt(' ',1)
      return
      end
      subroutine pget (q,p,iprop,znuc,zcor)
c...  p over active only, symmetry packed
c     zcor = 2*tr[over cor+fzc](p)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
      dimension q(*),p(*)
      dimension ix(8),ivec(8)
      ipp = 1
      do 10 isym=1,nirrr
10    ix(isym) = icorr((nsymao(isym)*(nsymao(isym)+1))/2)
      zcor=0
      call get1 (q(1),q(ix(1)),iprop,znuc)
c
      do 20 isym=1,nirrr
c      call outtri (q(ix(isym)),nsymao(isym),'sym block of ao prop')
20    ivec(isym) = icorr(nsymao(isym)**2)
      call qget (q(1),ivec)
c
      do 50 isym=1,nirrr
      n =nsymao(isym)
      if (n.eq.0) go to 50
c     nn = n**2
      nnt = (n*(n+1))/2
      m = nactt(isym)
      mc = ncor(isym) + ifreez(isym)
      do i=1,n
         q((i*(i+1))/2+ix(isym)-1) = q((i*(i+1))/2+ix(isym)-1)*0.5d0
      end do
c
c...  core contribution
c
      if (mc.ne.0) then
         id = icorr(nnt)
         ii = ivec(isym)
         idd = id
         do 30 i=1,n
            call mxmaa(q(ivec(isym)),1,n, q(ii),n,0, q(idd),1,0, i,mc,1)
            idd = idd+i
30       ii = ii+1
         zcor = zcor + 4.0d0*ddot(nnt,q(ix(isym)),1,q(id),1)
         call corlsr (id)
      end if
c
      if (m.eq.0) go to 50
c
c.... active transformation
c
      iy = icorr(n**2)
      ir = icorr(n*m)
      ip = icorr(m**2)
      call vclr(q(ir),1,n*m)
      call trankh (q(ir),q(ix(isym)),q(ivec(isym)+mc*n),m,n,1,n,q(iy))
c      call outsqr (q(ir),n,n,m,'r matrix')
      call mxmaa(q(ir),n,1, q(ivec(isym)+mc*n),1,n, q(ip),1,m, m,n,m)
c      call outsqr (q(ip),m,m,m,'p after mxmaa')
      do 40 i=1,m
      do 40 j=1,i
      p(ipp) = q(ip+(i-1)+(j-1)*m) + q(ip+(j-1)+(i-1)*m)
40    ipp = ipp+1
      call corlsr (iy)
50    continue
      call corlsr (ix(1))
      return
      end
      subroutine mcprnt(q,iprin,iwrite,ipunch)
c
c     printing routine. called by mcanal (lprint) & wvfn (iprint)
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      logical btest
INCLUDE(common/multic)
INCLUDE(common/detcic)
INCLUDE(common/jobopt)
INCLUDE(common/syminf)
      common /couple/ surd(511)
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1              ,jorb(511),korb(511),lorb(511)
      common /intbuf/  intpos,intfil,intmod
INCLUDE(common/mcff)
INCLUDE(common/prnprn)
      character*24 buff
      character*9 word
      dimension q(*),iqpos(8),id(8),occno(31),scra(31),orbnat(31)
      data ifzcpr /0/
c
      ibase = icorr(0)
c
      if (btest(iprin,12).or.btest(iprin,10).or.
     +    opunch(17)) then
      ivec = icorr(nci*nstate)
      if(iguess.eq.0) call caserr('no ci vector provided')
      call cget(q(ivec),nstate)
      end if
      if (btest(iprin,9).or.btest(iprin,12).or.btest(iprin,13)) then
      do 20 isym=1,nirrr
20    iqpos(isym) = icorr(nsymao(isym)**2)
      call qget (q(1),iqpos)
      end if
c
c... =integral=
      if (btest(iprin,3)) then
      intpos=0
      intfil=num6
      iblf  =iblk6
      intmod=0
      ibuff = icorr(maxbas**2)
      write(iwrite,1000)
1000  format(' two electron integrals'
     1     /' ----------------------')
      do 80 i=1,nprim
      do 80 j=1,i
      isymij = mults(itype(i),itype(j))
      write(iwrite,1001)i,j
1001  format  (' active orbital pair',2i3)
      do 60 isyma=1,nirrr
      write (buff,30) ' coulomb',isyma
30    format(a8,' for symmetry',i2)
      isymb = mults(isyma,isymij)
      if (nsymm(isyma)*nsymm(isymb).eq.0) goto 60
      if (isymb-isyma) 40,50,60
40    call intin (q(ibuff),nsymm(isyma)*nsymm(isymb))
      call outsqr (q(ibuff),nsymm(isyma),nsymm(isyma),nsymm(isymb),buff)
      goto 60
50    call intin (q(ibuff),(nsymm(isyma)*(nsymm(isyma)+1))/2)
      call outtri (q(ibuff),nsymm(isyma),buff)
60    continue
      do 70 isyma=1,nirrr
      isymb = mults(isyma,isymij)
      if (nsymm(isyma)*nsymm(isymb).eq.0) goto 70
      call intin (q(ibuff),nsymm(isyma)*nsymm(isymb))
      write (buff,30) 'exchange',isyma
      call outsqr (q(ibuff),nsymm(isyma),nsymm(isyma),nsymm(isymb),buff)
70    continue
80    continue
      write(iwrite,1002)
 1002 format(' one electron integrals'
     1     /' ----------------------')
      icount = 0
90    do 110 isyma=1,nirrr
      if (nsymm(isyma).eq.0) goto 110
      call intin (q(ibuff),(nsymm(isyma)*(nsymm(isyma)+1))/2)
      write (buff,100) isyma
100   format('integrals for symmetry',i2)
      call outtri (q(ibuff),nsymm(isyma),buff)
110   continue
      if (icount.eq.1) goto 120
      icount = icount + 1
      write(iwrite,1003)
 1003 format ( ' g matrix'/' --------')
      goto 90
 120  write(iwrite,1004)
 1004 format(' lagrangian matrix'/' -----------------')
      print *,' not available anymore due to anal mcgrad switch '
c      do 130 isyma=1,nirrr
c      call intin (q(ibuff),nsymm(isyma)*nprm(isyma))
c      write (buff,100) isyma
c130   call outsqr (q(ibuff),nsymm(isyma),nsymm(isyma),nprm(isyma),buff)
c      call corlsr (ibuff)
      end if
c
c...  =orbitals=
      if (btest(iprin,9)) then
      if(mcprin)write(iwrite,140)
140   format(/1x,'molecular orbital coefficients'
     1      /' ------------------------------')
      do 150 isym=1,nirrr
150   id(isym) = 0
      mmax = nfreez+nprim
      if (btest(iprin,13)) mmax=nbasao
      do 170 i=1,mmax
      if (i.gt.nfreez) isym=itype(i-nfreez)
      if (i.le.nfreez) isym=ifzsym(i)
      kk=id(isym)+1
      id(isym)=kk
      if (ifzcpr.ne.0.and.i.le.nfreez) goto 170
      word='(active)'
      if (i.le.ncoremc+nfreez) word='(core)'
      if (i.le.nfreez) word='(frozen)'
      if (i.gt.nprim+nfreez) word = '(virtual)'
      if(mcprin)write(iwrite,160)
     *i,word,isym,(q(iqpos(isym)+(kk-1)*nsymao(isym)-1+j),
     1j=1,nsymao(isym))
160   format(/i4,1x,a9,1x,'sym',i2,2x,
     1'coeffs',70(t29,9f11.7/) )
170   continue
      ifzcpr=1
      end if
c
c... =natorb=
      if (btest(iprin,12)) then
      igam = icorr(ne)
      do 260 istat1=1,nstate
      do 260 istat2=1,istat1
      call denst1(q(1),q(ivec+(istat1-1)*nci),q(ivec+(istat2-1)*nci),
     * q(igam))
      if (istat1.eq.istat2) then
      if(mcprin)write(iwrite,180)istat1
180   format(/1x,'one-particle density matrix for state',i2/
     1      ' ---------------------------------------')
      if(mcprin)call outsqr (q(igam+ic1d),nact,nact,nact,' ')
c...  natural orbital analysis
      if(mcprin)write(iwrite,190)
190   format(/1x,'natural orbitals'/' ----------------')
      ifail=0
      call f02abf (q(igam+ic1d),nact,nact,occno,q(igam),nact,scra
     * ,ifail)
      do 240 i=1,nact
      occ = occno(nact+1-i)
      call fmove (q(igam+(nact-i)*nact),orbnat,nact)
_IF(vax)
      zmax = 0.0d0
      do 1 j=1,nact
      if ( dabs(orbnat(j)).le.zmax) goto 1
      mmax = j
      zmax =  dabs(orbnat(j))
   1  continue
_ELSE
      mmax=idamax(nact,orbnat,1)
      zmax=dabs(orbnat(mmax))
_ENDIF
      isym = itypea(mmax)
      if (zmax.lt.0.0d0) then
       call dscal(nact,-1.0d0,orbnat,1)
         endif
      kk=ncor(isym)+ifreez(isym)
      call vclr(g,1,nsymao(isym))
      do 210 j=1,nact
      if (itypea(j).ne.isym) goto 210
      ix1 = iqpos(isym)+kk*nsymao(isym)
      call daxpy(nsymao(isym),orbnat(j),q(ix1),1,g,1)
      kk = kk + 1
210   continue
      if (mcprin) then
       mfg = i+nfreez+ncoremc
       write(iwrite,220)mfg ,occ,isym,(orbnat(j),j=1,nact)
       write(iwrite,230)(g(j),j=1,nsymao(isym))
220    format(/i4,f10.7,1x,'sym',i2,2x,
     1 'e-vec',8(t29,9f11.7/) )
230    format(/22x,'coeffs',70(t29,9f11.7/) )
      endif
240   continue
      else
      write(iwrite,250)istat1,istat2
250   format(/1x,'one-particle transition density matrix for states',2i2
     1 /    ' -----------------------------------------------------')
      call outsqr (q(igam+ic1d),nact,nact,nact,' ')
      end if
260   continue
      call corlsr (igam)
      end if
c
c... =civector=
      if (btest(iprin,10).or.opunch(17)) call prici(q(1),q(1),
     +                                   ivec,iwrite,ipunch)
c
c... =formulae=
      if (btest(iprin,2) .and. iexc.ge.0) then
      call rdedx(surd,511,iblkft,numft)
      write(iwrite,360)
360   format(/1x,'the formula tape'/' ----------------'/
     1/1x,'1-particle formulae:'
     2 /9x,'walks',10x,'value',10x,'i',3x,'j')
      call find (numft)
370   call get (g,nw)
      if (nw.eq.0) goto 400
      call find (numft)
      call unpkft (nw)
      do 380 i=1,nw
      iva=ival(i)
380   write(iwrite,390) ii(i),jj(i),surd(iva),iorb(i),jorb(i)
390   format(2i8,f15.7,5x,4i4)
      goto 370
400   write(iwrite,410)
410   format(/1x,'2-particle formulae:'
     1 /9x,'walks',10x,'value',10x,'i',3x,'j',3x,'k',3x,'l')
      call find(numft)
420   call get (g,nw)
      if (nw.eq.0) return
      call find (numft)
      call unpkft (nw)
      do 430 i=1,nw
      iva=ival(i)
 430  write(iwrite,390)
     * ii(i),jj(i),surd(iva),iorb(i),jorb(i),korb(i),lorb(i)
      goto 420
      end if
c
      call corlsr (ibase)
c
      return
      end
      subroutine hedout(qeig,qocc)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (mxorb1=maxorb+1)
      character *7 typ
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jinfo)
INCLUDE(common/machin)
      logical iftran
      character *8 param,gitle
      common /junkc / param(19),gitle(10)
INCLUDE(common/runlab)
      common/blkorbs/eigout(maxorb),pop(mxorb1),
     1                    nbas,newb,ncoll,ivalue,ioccc,ipad
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     2        ctran(mxorb3),iftran,iftri
      dimension qeig(*),qocc(*),mapie(maxorb)
      data m29/29/
      param(1)=zanam
      param(2)=zdate
      param(3)=ztime
      param(4)=' gamess'
      param(5)='mcscf'
      do 70 i=1,10
 70   gitle(i)=ztitle(i)
      nbas=nbasao
      newb=nbas
      ncoll=nbas
      ivalue=1
      ioccc=1
      iad=1
      do 10 isym=1,nirrr
      do 10 i=1,nbasao
      if(i-nfreez)11,11,12
 12   is=itype(i-nfreez)
      go to 13
 11   is=ifzsym(i)
 13   if(is.ne.isym)go to 10
      mapie(iad)=i
      iad=iad+1
 10   continue
      iad=0
      ien=0
      ioc=1
      do 20 isym=1,nirrr
      n=nsymao(isym)
      do 30 i=1,n
      occ=2.0d0
      typ='frozen '
      if(i.gt.ifreez(isym))typ='core   '
      if(i.gt.ifreez(isym)+ncor(isym))typ='active '
      if(i.gt.ifreez(isym)+nprm(isym))typ='virtual'
      if(typ.eq.'virtual')occ=0.0d0
      if(typ.ne.'active')go to 40
      occ=-qocc(ioc)
      ioc=ioc+1
 40   val=qeig(ien+i)
      iad=iad+1
      eigout(mapie(iad))=val
 30   pop(mapie(iad))=occ
 20   ien=ien+n
      pop(mxorb1)=enext
      iblko=iblkq-(lensec(m29)+lensec(mach(8))+lensec(mach(9)))
      call wrtc(param,m29,iblko,num3)
      call wrt3s(eigout,mach(8),num3)
      nav = lenwrd()
      call wrt3is(ilifc,mach(9)*nav,num3)
c
      return
      end
      subroutine mcgrad(q,iq)
      implicit REAL  (a-h,o-z)
      logical otran,exist
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/harmon)
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),ctran(mxorb3
     *),otran
INCLUDE(common/multic)
INCLUDE(common/syminf)
      logical cigr,cicv,mpgr,mcgr,cicx,umpgr
      common/cigrad/cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,mnnr,
     1  mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,isecnd,isecsy,
     2  irlagr,iadfrc,nfc,intlgr,umpgr,ispaer(26),nd2mo,
     3  newcor,ncact,nvr,ifilh,iblkh,iblk1,
     4  iblk22,ntot,nupact,ijr3,rspare(70)
      common/intbuf/intpos,intfil,intmod,intnw,gax(511)
      common/lsort/g(511),ii(511),jj(511),ival(511),iorb(511),
     1 jorb(511),korb(511),lorb(511)
      common/couple/eig(511)
INCLUDE(common/atmblk)
      common/scfblk/en,etot
      common/restri/nfilez(63),lda(508),isect(508),ldx(508)
      dimension iqpos(8),q(*),iq(*)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      etot=enext
      nav = lenwrd()
_IF1()c     write(iwr,119)enext
_IF1()c119   format(/' total energy',t22,f22.12)
c
c     isecdd=101
c     isecll=102
c     iscigr=103
c     isecmo=107
c     isecnd=105
c     isecsy=106
c     ifil2d=5
c     iblk2d=1
      mtyp=0
      lncigr=70+lenint(60)
c
c dimensions
c
      nb2=nbasao*nbasao
      nb22=(nb2+nbasao)/2
c
      nfc=ncoremc+nfreez
      iadfrc=nfc*nbasao
c     nfca=nfc+nact
      newcor=nfc
      ncact=nact
c
c density matrix is first built over active mo's and then over frozen
c and core,and transformed to the ao basis
c
      ng1ns=nact*nact
c
c lagrangian matrix is built over primary orbitals (core+active) (row)
c and (core+active+virtual) (columns)
c
      nlag=(nprim+nfreez)*nbasao
c
c the dimension is actually larger to store an intermediate matrix in
c the transformation to the ao basis
c
c addresses
c
      lenciv=nci*nstate
      ibase=icorr(0)
      ivec=icorr(lenciv)
      indnwa=icori(nact)
      icmo=icorr(nb2)
      isymnw=icori(nbasao)
      indnw=icori(nbasao*nirrr)
c
c ivec address of the ci vector
c icmo address of the mo coefficients
c isymnw address of symmetry index vector of the reordered mo's
c indnw address of permutation vector
c (isym-1)*nbasao+i defines the position in the reordered mo matrix
c of the ith orbital of symmetry isym.
c
      do 1 isym=1,nirrr
1     iqpos(isym)=icorr(nsymao(isym)*nsymao(isym))
c
c initialise the mo matrix to zero
c
      itmp = iqpos(nirrr)+nsymao(nirrr)*nsymao(nirrr)
      call vclr(q(1),1,itmp)
c
c read in the symmetry blocked mo matrix
c
      call qget(q(1),iqpos)
c
c reorder the mo matrix in the order frozen,core,active,virtual
c and initialise the indexes for frozen,core and active orbitals
c
      if=0
      ic=nfreez
      ia=ncoremc+ic
      iv=nact+ia
c
c loop over symmetries
c
      do 2 isym=1,nirrr
c
c count type of each orbital within the symmetry
c
      nf=ifreez(isym)
      nc=ncor(isym)+nf
      na=nactt(isym)+nc
      isymnd=(isym-1)*nbasao - 1
c
c the number of orbitals of a particular symmetry
c
      norb=nsymao(isym)
      if(norb.ne.0) then
      do 3 i=1,norb
c
c the orbital type in that symmetry
c
      if(i.le.nf)goto 4
      if(i.le.nc)goto 5
      if(i.le.na)goto 6
      iiorb=iv
      iv=iv+1
      goto 7
4     iiorb=if
      if=if+1
      goto 7
5     iiorb=ic
      ic=ic+1
      goto 7
6     iiorb=ia
      ia=ia+1
c
c define the appropriate component of the permutation vector
c and the new symmetry vector
c
7     iq(indnw+isymnd+i)=iiorb+1
      iq(isymnw+iiorb)=isym
c     iperm=indnw+isymnd+i
c     isyvec=isymnw+iiorb
c
c define the addresses of the orbital in the symmetry blocked and in
c the new cmo matrix
c
      iold=iqpos(isym)+(i-1)*norb
      inew=icmo+iiorb*nbasao+istart(isym)-1
      call dcopy(norb,q(iold),1,q(inew),1)
3     continue
c
      endif
2     continue
c
c release cpu
c
      call corlsr(iqpos(1))
c
c     new indexes to permute density matrix elements
c
      newind=0
      do 14 isym=1,nirrr
      do 46 iact=1,nact
      if(itypea(iact).eq.isym)then
      newind=newind+1
      iq(indnwa-1+iact)=newind
      end if
   46 continue
   14 continue
c
c address of the ao symmetry adaption matrix
c
      icsym=icorr(nb2)
c
c generate symmetry adaption matrix
c
c     mtyp=0
c     call secget(isect(51),mtyp,iblok)
c     iblok=iblok+mvadd
c     call rdedx(q(icsym),nb2,iblok,num3)
c
      if(otran)call caserr(
     *'the basis set is not symmetry adapted')
      call vclr(q(icsym),1,nb2)
      do 100 i=1,nbasao
      iii=ilifq(i)+icsym-1
 100  q(iii+i)=1.0d0
      call tdown(q(icsym),ilifq,q(icsym),ilifq,nbasao)
c
      icpr=icorr(nb2)
c
c icsym address of the symmetry adaption matrix
c icpr address of product
c
      call vclr(q(icpr),1,nb2)
      call mxmb(q(icsym),1,nbasao,q(icmo),1,nbasao,q(icpr),1,nbasao,
     1 nbasao,nbasao,nbasao)
      call dcopy(nb2,q(icpr),1,q(icmo),1)
      call corlsr(icsym)
c
c the mo expression in the ao basis in icmo.
c address of the one electron density matrix ig1ns.we require more
c cpu to perform the transformation to the ao basis
c
      ng2ns=ng1ns*(ng1ns+1)/2
      ig2ns=icorr(ng2ns)
      ig1ns=icorr(ng1ns)
c
c read ci vector (fundamental state)
c
      call cget(q(ivec),nstate)
      call vclr(q(ig2ns),1,ng2ns+ng1ns)
      call densav (q(1),q(ivec),q(ig2ns))
      call symden (q(ig2ns),nact)
c
c now form the full one electron density matrix over active orbitals
c
      ng1s=ng1ns
      ig1s=icorr(ng1s)
      call vclr(q(ig1s),1,ng1s)
c
c now fill density matrix with the active components and symmetrise.
c
      do 15 i=1,nact
      is=iq(indnwa-1+i)
      do 15 j=1,i
      js=iq(indnwa-1+j)
      ijs=(js-1)*nact+is-1
      jis=(is-1)*nact+js-1
      qav=q(ig1ns+(i-1)*nact+j-1)
      q(ig1s+ijs)=qav
15    q(ig1s+jis)=qav
      call dcopy(ng1s,q(ig1s),1,q(ig1ns),1)
      call corlsr(ig1s)
c
c transform to the ao basis and get cpu for the product
c
      ninter=nbasao*nact
      interm=icorr(ninter)
      call vclr(q(interm),1,ninter)
      call mxmb(q(icmo+nfc*nbasao),1,nbasao,q(ig1ns),1,nact,q(interm)
     1 ,1,nbasao,nbasao,nact,nact)
      do js=0, ninter-1
        q(ig1ns+js) = q(interm+js)
      enddo
      call corlsr(ig1ns+ninter)
c
c get cpu to perform the final matrix multiplication
c
      ifinal=icorr(nb2)
      call vclr(q(ifinal),1,nb2)
      call mxmb(q(ig1ns),1,nbasao,q(icmo+nfc*nbasao),nbasao,1,q(ifinal)
     1 ,1,nbasao,nbasao,nact,nbasao)
      do js=0, nb2-1
        q(ig1ns+js) = q(ifinal+js)
      enddo
      call corlsr(ig1ns+nb2)
c
c make a triangular matrix
c
      ilast=icorr(nb22)
      ij=0
      do 16 i=1,nbasao
      do 16 j=1,i
      q(ilast+ij)=q(ig1ns+(j-1)*nbasao+i-1)
   16 ij=ij+1
c
c write on tape:the mo matrix (ao basis) ,permutation and symmetry
c vectors and the density matrix
c
c transfer to file the ci gradient where lncigr is the length
c of the cigrad common block
c
      mtyp=0
      length=lensec(lncigr)
      call secput(iscigr,mtyp,length,iblok)
      call revind
      call wrt3(cigr,lncigr,iblok,num3)
c
      mtyp=0
      length=lensec(nb2)+mvadd
      call secput(isecmo,mtyp,length,iblok)
      call revind
c     call outsqr(q(icmo),nbasao,nbasao,nbasao,'mo coefficients')
      call wrt3(q(icmo),nb2,iblok+mvadd,num3)
      mtyp=0
      lenth1=nbasao*nirrr
      length=lensec(lenth1)
      call secput(isecnd,mtyp,length,iblok)
      call revind
      call wrt3i(iq(indnw),lenth1*nav,iblok,num3)
      mtyp=0
      length=lensec(nbasao)
      call secput(isecsy,mtyp,length,iblok)
      call revind
      call wrt3i(iq(isymnw),nbasao*nav,iblok,num3)
      mtyp=0
      length=lensec(nb22)
      call secput(isecdd,mtyp,length,iblok)
      call revind
c     call outtri(q(ilast),nbasao,'transformed one pdm')
      call wrt3(q(ilast),nb22,iblok,num3)
c
c release cpu not needed
c
      call corlsr(ig1ns)
c
c now permute the full two electron density
c
      ng2s = (nact*(nact+1))/2
      ng2s = (ng2s*(ng2s+1))/2
      ig2s=icorr(ng2s)
      do 29 i=1,nact
      is=iq(indnwa-1+i)
      do 29 j=1,i
      js=iq(indnwa-1+j)
      ijs=ind(is,js)
      ij=(j-1)*nact+i
      do 29 k=1,i
      ks=iq(indnwa-1+k)
c
      lm=k
c
      if(i.eq.k)lm=j
      do 29 l=1,lm
      ls=iq(indnwa-1+l)
      kls=ind(ks,ls)
      ijkls=ind(ijs,kls)
      kl=(l-1)*nact+k
      ijkl=ind(ij,kl)-1
29    q(ig2s+ijkls-1)=q(ig2ns+ijkl)*0.5d0
c     if(nprint.ne.-4)goto 410
c     write(6,400)(q(ig2s+i-1),i=1,ng2s)
c400   format(' 2-pdm =',6f15.8)
c      print*,'writing two particle density matrix to block',
c     <iblk2d,' on file ',ifil2d
c      call outtri (q(ig2s),(nact*(nact+1))/2,'tpdm')
c     print *,' in mcgrad  iblk2d  ifil2d ng2s   ',iblk2d,ifil2d,ng2s
      call wrt3(q(ig2s),ng2s,iblk2d,ifil2d)
      call corlsr (ig2ns)
c
c find lagrangian on integral file
c
      call search(irlagr,num6)
      call find(num6)
      call get(gax,intnw)
      if(intnw.eq.0)call caserr('end of file reached')
      intnw=intnw+1
      intpos=intlgr
c
c find maximum dimension of symmetry block
c
      ndmax=0
      do 17 i=1,nirrr
17    ndmax=max(ndmax,nsymm(i)*nprm(i))
      ilag=icorr(nlag)
c     ilgblk=icorr(ndmax)
      ilgblk = icorr(nbasao**2)
      intmod=0
c
c initialise full lagrangian
c
      call vclr(q(ilag),1,nlag)
      do 18 isym=1,nirrr
      norb=nsymm(isym)
      norbpr=nprm(isym)
      nblok=norb*norbpr
      if (oharm) nblok=nsym0(isym)*norbpr
      call intin(q(ilgblk),nblok)
      if (oharm)  call exphvs(q(ilgblk),isym)
c**** debug
c      write(6,2998)
c      write(6,3000)(q(ilgblk+i-1),i=1,nblok)
c2998  format(' q(ilgblk) = ')
c3000  format(1x,8f9.5)
c**** debug
c
c lagrangian block in q(ilgblk)
c find addresses in full lagrangian matrix and change sign
c
      ishf=ifreez(isym)
      iadpr=(isym-1)*nbasao+indnw
      if(norbpr.eq.0)goto 18
      do 19 j=1,norbpr
      jtr=iq(iadpr+j+ishf-1)-nfreez
      ilagbb=(j-1)*norb
      ilagb=(jtr-1)*nbasis
      do 19 i=1,norb
      itr=iq(iadpr+i+ishf-1)-nfreez
      indold=ilgblk+ilagbb+i-1
      indnew=ilag+ilagb+itr-1
19    q(indnew)=-q(indold)
18    continue
      call corlsr(ilgblk)
c**** debug
c      write(6,2999)
c      write(6,3000)(q(ilag+i-1),i=1,nlag)
c2999  format(' q(ilag) = ')
c**** debug
c
c     consider frozen orbitals
c
      if(nfreez.gt.0) then
       mtyp = 0
       call secloc(isect(9),exist,iblok)
       if (exist) then
         call rdedx(eig,lda(isect(9)),iblok,num3)
       else
         call caserr('eigenvalues not found')
       end if
       minter=nb2
       jnterm=icorr(minter)
       call vclr(q(jnterm),1,minter)
       do 1000 j=1,nfreez
       ij = jnterm+(j-1)*nbasao+j-1
       q(ij) = -2.0d0*eig(j)
 1000  continue
       mi = 0
       mj = nfreez*(nbasao+1)
       do 120 i=1,nbasis
       call dcopy(nbasis,q(ilag+mi),1,q(jnterm+mj),1)
       mi=mi+nbasis
  120  mj=mj+nbasao
c**** debug
c      mj=(nact+ncoremc+nfreez)*nbasao
c      write(6,2995)
c      write(6,3000)(q(jnterm+i-1),i=1,mj)
c2995  format(' q(jnterm) = ')
c**** debug
c
c perform the first part of the transformation to the ao basis
c
       ninter=nb2
       interm=icorr(ninter)
       call vclr(q(interm),1,ninter)
       call mxmb(q(icmo),1,nbasao,q(jnterm),1,nbasao,
     1  q(interm),1,nbasao,nbasao,nbasao,nprim+nfreez)
       do js=0, ninter-1
         q(ilag+js) = q(interm+js)
       enddo
       call vclr(q(interm),1,nb2)
       call mxmb(q(ilag),1,nbasao,q(icmo),nbasao,1
     1  ,q(interm),1,nbasao,nbasao,nprim+nfreez,nbasao)
      else
c
c perform the first part of the transformation to the ao basis
c
      ninter=nb2
      interm=icorr(ninter)
      call vclr(q(interm),1,ninter)
      call mxmb(q(icmo+nfreez*nbasao),1,nbasao,q(ilag),1,nbasis,
     1 q(interm),1,nbasao,nbasao,nbasis,nprim)
      do js=0, ninter-1
        q(ilag+js) = q(interm+js)
      enddo
      call vclr(q(interm),1,nb2)
      call mxmb(q(ilag),1,nbasao,q(icmo+nfreez*nbasao),nbasao,1
     1 ,q(interm),1,nbasao,nbasao,nprim,nbasao)
c
      endif
c
c make lagrangian triangular over ao basis
c
      ij=0
      do 20 j=1,nbasao
      do 20 i=1,j
      q(ilag+ij)=q(interm+(j-1)*nbasao+i-1)
   20 ij=ij+1
c
c write on dumpfile
c
      mtyp=0
      length=lensec(nb2)
      call secput(isecll,mtyp,length,iblok)
      call revind
c****debug
c      call outtri(q(ilag),nbasao,'transformed lagrangian')
c****debug
      call wrt3(q(ilag),nb22,iblok,num3)
      call corlsr(icmo)
c
      call clredx
      call corlsr(ibase)
      return
      end
_IF()
      subroutine priop(a,isy,np,string,nr,iwrite)
      implicit REAL (a-h,o-z)
      dimension a(*)
      character*(*) string
INCLUDE(common/sizes)
INCLUDE(common/mcff)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/jobopt)
INCLUDE(common/mcscra)
c
      write(iwrite,10) string,nr
10    format(/1x,a8,i4)
      ioff=0
      do 70 isa=1,nirrr
      isb=mults(isy,isa)
      na=nt(isa)
      nb=nt(isb)
      if(na.eq.0.or.nb.eq.0) goto 70
      if(np.ne.0.and.isb.gt.isa) goto 70
      write(iwrite,20) isy,isa
20    format(' block',i2,'.',i1)
      if(np.ne.0.and.isa.eq.isb) goto 50
      ij=ioff+1-na
      ioff=ioff+na*nb
      do 40 i=1,na
      write(iwrite,30) (a(ij+j*na),j=1,nb)
30    format(1x,10f13.8)
40    ij=ij+1
      goto 70
50    ij=ioff
      do 60 i=1,na
      write(iwrite,30) (a(ij+j),j=1,i)
60    ij=ij+i
      ioff=ij
70    continue
      return
      end
_ENDIF
      subroutine reduca(a,b,n,m)
      implicit REAL  (a-h,o-z)
      dimension a(n,n),b(n,m)
      do 20 j=1,m
      do 10 i=1,m
  10  b(i,j)=b(i,j)+a(i,j)
      do 20 i=m+1,n
  20  b(i,j)=b(i,j)+a(i,j)-a(j,i)
      return
      end
      subroutine expda(a,b,n,m)
      implicit REAL  (a-h,o-z)
      dimension a(n,m),b(n,n)
      call vclr(b,1,n**2)
      do 20 j=1,m
      do 10 i=1,m
  10  b(i,j)=-0.5d0*(a(i,j)+a(j,i))
      do 20 i=m+1,n
      b(i,j)=-0.5d0*a(i,j)
  20  b(j,i)=b(i,j)
      return
      end
      subroutine trankh (z,g,q,nocc,n,minr,maxr,gg)
c
c     kernel for "h" & "j" matrix and others in 4-index
c     g is in triangular form
c     the logic of the code is
c
c     i=1,nocc
c       s = 1,maxr
c         r = max(minr,s),maxr
c           z(r,i) = g(rs) * q(s,i) + z(r,i)
c     since r>=s, copy g into work array gg for linear addressing
c     on a scalar machine, probably best to eliminate copying
c     and use inner product method.
c
      implicit REAL  (a-h,o-z)
      dimension z(*),g(*),q(*),gg(*)
      if (nocc*(maxr-minr+1).le.0) return
      do 1 ir=minr,maxr
      do 1 is=1,ir
1     gg(ir+(is-1)*maxr) = g((ir*(ir-1))/2+is)
c
      do 2 i=1,nocc
      do 2 is=1,maxr
      if (q(is+(i-1)*n).eq.0.0d0) goto 2
      mr = max(minr,is)
      call daxpy(maxr-mr+1,q(is+(i-1)*n),gg(mr+(is-1)*maxr),1
     * ,z(mr+(i-1)*n),1)
2     continue
      return
      end
      subroutine tranki (z,g,q,nocc,n,minr,maxr)
c
c     kernel for "i" matrix and others in 4-index
c     g is in triangular form
c     the logic of the code is
c
c     i=1,nocc
c       r=minr,maxr
c         s=1,r
c           z(s,i) = g(rs) * q(r,i) + z(s,i)
c
      implicit REAL  (a-h,o-z)
      dimension z(*),g(*),q(*)
      if (nocc*(maxr-minr+1).le.0) return
      do 2 i=1,nocc
      do 2 ir=minr,maxr
      if (q(ir+(i-1)*n).eq.0.0d0) goto 2
      irs = (ir*(ir-1))/2+1
      call daxpy(ir
     +  ,q(ir+(i-1)*n),g(irs),1,z(1+(i-1)*n),1)
2     continue
      return
      end
      subroutine rdfrt (q,len,jork,jknum)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c...  read a j or k operator from 1st half of file
INCLUDE(common/mcff)
INCLUDE(common/multic)
INCLUDE(common/syminf)
      common /lsort  / g(511)
      dimension q(*)
c
      iadd=0
      if(jknum) 25,20,5
5     nop=nprim*(nprim+1)/2
      iop=jknum
      if(iop.gt.nop) then
      iop=iop-nop
      iadd=iword2
      end if
      if (jork.eq.2) goto 10
c...  must be coulomb operator
      call readg (q,len,jad(iop)+iadd,num6)
      return
c...  exchange operator
10    call readg (q,len,kad(iop)+iadd,num6)
      return
c...  core fock matrix fc if code 0
20    if(iadr.ne.0) iadd=iword2
      call readg (q,len,lad+iadd,num6)
      return
c.....bare hamiltonian if code<0
25    if(iadr.ne.0) iadd=iword2
      call readg(q,len,lad+iadd+ltrimo,num6)
      return
      end
      subroutine wtfrt (q,len,jork,jknum)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c...  write a j or k operator from 1st half of file
INCLUDE(common/mcff)
INCLUDE(common/multic)
INCLUDE(common/syminf)
      common /lsort  / g(511)
      dimension q(*)
c
      iadd=0
      if(jknum) 45,40,26
26    iop=jknum
      nop=nprim*(nprim+1)/2
      if(iop.gt.nop) then
      iop=iop-nop
      iadd=iword2
      end if
      if (jork.eq.2) goto 30
c...  must be coulomb operator
      call wrtg (q,len,jad(iop)+iadd,num6)
      return
c...  exchange operator
30    call wrtg (q,len,kad(iop)+iadd,num6)
      return
c...  fc matrix if code 0
40    if(iadw.ne.0) iadd=iword2
      call wrtg (q,len,lad+iadd,num6)
      return
c...  bare hamiltonian for code<0
45    if(iadw.ne.0) iadd=iword2
      call wrtg(q,len,lad+iadd+ltrimo,num6)
      return
      end
      subroutine finit
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c...  read a j or k operator from 1st half of file
INCLUDE(common/iofile)
INCLUDE(common/mcff)
INCLUDE(common/multic)
INCLUDE(common/syminf)
      common /lsort  / g(511)
      ifinit=1
      iadr=0
      iadw=nprim*(nprim+1)/2
c...initialise for direct access on mo integrals
      call search (iblk6,num6)
      if (nprim.gt.mcprim) then
        write(iwr,'(1x,"MCPRIM = ",i4)')mcprim
        call caserr(
     +  'no. of active orbitals in excess of allowed maximum (MCPRIM)')
      endif
      do 70 is=1,nirrr
      lj(is) = 0
      lk(is) = 0
      do 70 is1=1,nirrr
      is2 = mults(is,is1)
      if (is1-is2) 70,60,50
50    lj(is) = lj(is) + nsymm(is1)*nsymm(is2)
      goto 70
60    lj(is) = lj(is) + (nsymm(is1)*(nsymm(is1)+1))/2
70    lk(is) = lk(is) + nsymm(is1)*nsymm(is2)
c      print*,'** lj,lk ',lj,lk
c
      iword2 = 0
      ij = 0
      do 80 i=1,nprim
      isymi = itype(i)
      do 80 j=1,i
      isymj = itype(j)
      ij = ij+1
      jad(ij) = iword2
      iword2 = iword2 + lj(mults(isymi,isymj))
      kad(ij) = iword2
      iword2 = iword2 + lk(mults(isymi,isymj))
c      print*,'80 loop; i,j,jad(ij),kad(ij) ',i,j,jad(ij),kad(ij)
80     continue
      lad = iword2
      iword2 = iword2 + 3*ltrimo
      do 90 is=1,nirrr
90    iword2 = iword2 + nprm(is)*nsymm(is)
c
c.... this marks the end of the list
      iblf1 = iblk6
      iblf2 = lensec(iword2)+iblk6
c.... assume that we will next use second list
      nbl = iblf2-iblf1
      iword2 = nbl*511
c...  we must write dummy records of the correct length
      call search (iblf2,num6)
      do 100 i=1,nbl
100   call put (g,511,num6)
c
      call offset
c
      return
      end
      subroutine offset
      implicit REAL  (a-h,o-z)
c...  read a j or k operator from 1st half of file
INCLUDE(common/sizes)
INCLUDE(common/mcff)
c
      iblf = iblf2
c...   disconnect atmol buffers
_IFN1(iv)       call kilbuf
      return
      end
      subroutine freset
      implicit REAL  (a-h,o-z)
c...  read a j or k operator from 1st half of file
INCLUDE(common/sizes)
INCLUDE(common/mcff)
c
      iblf = iblf1
_IFN1(iv)       call kilbuf
      return
      end
      subroutine redadd(q,a,fac,np)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/mcscra)
      dimension q(*),a(*)
c
      ijd=0
      ijq=0
      do 60 is=1,nirrr
      n=nt(is)
      jiq=ijq-n+1
      do 50 i=1,n
      if(np.eq.0) goto 20
      do 10 j=1,i
10    a(ijd+j)=a(ijd+j)+fac*(q(ijq+j)+q(jiq+j*n))
      jiq=jiq+1
      goto 40
20    do 30 j=1,i
30    a(ijd+j)=a(ijd+j)+fac*q(ijq+j)
40    ijd=ijd+i
50    ijq=ijq+n
60    continue
      return
      end
      subroutine redsub(q,a)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/mcscra)
      dimension a(*),q(*)
      ijd=0
      ijq=0
      do 30 is=1,nirrr
      n=nt(is)
      do 20 i=1,n
      do 10 j=1,i
10    a(ijd+j)=a(ijd+j)-q(ijq+j)
      ijd=ijd+i
20    ijq=ijq+n
30    continue
      return
      end
      subroutine reduce(q,a,n)
      implicit REAL  (a-h,o-z)
      dimension a(*),q(*)
      if(n.eq.0) return
      ijd=0
      ijq=0
      do 20 i=1,n
      do 10 j=1,i
10    a(ijd+j)=q(ijq+j)
      ijd=ijd+i
20    ijq=ijq+n
      return
      end
      subroutine tradd(a,b,n,m)
      implicit REAL  (a-h,o-z)
      dimension a(n,m),b(m,n)
      do 20 i=1,n
      do 20 j=1,m
20    a(i,j)=a(i,j)+b(j,i)
      return
      end
      subroutine transp(a,b,isy,vec)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/multic)
INCLUDE(common/syminf)
INCLUDE(common/mcscra)
      dimension a(*),b(*),vec(*)
c
c.....transposes matrix a. result in array b. a and b can be the same
c
      do 50 isa=1,nirrr
      isb=mults(isy,isa)
      if(isb.gt.isa) goto 50
      na=nt(isa)
      nb=nt(isb)
      if(na.eq.0.or.nb.eq.0) goto 50
      ij=ioffq(isy,isa)
      ji=ioffq(isy,isb)+1-nb
      je=na
      do 40 i=1,nb
      if(isa.eq.isb) je=i
      do 10 j=1,je
10    vec(j)=a(ij+j)
      do 20 j=1,je
20    b(ij+j)=a(ji+j*nb)
      do 30 j=1,je
30    b(ji+j*nb)=vec(j)
      ij=ij+na
40    ji=ji+1
50    continue
      return
      end
      subroutine wrtg (a,ll,iword,iunit)
      implicit REAL  (a-h,o-z)
c...   write l words from a to offset iword on file iunit
      dimension buf(512),ibuf(512),a(*)
      equivalence (buf(1),ibuf(1))
      l=ll
      ia=1
      ibl = iword/511
      iwordb = iword-ibl*511
      ibl = ibl+1
      if (iwordb.eq.0) goto 500
c..   first word not on block boundary
      call getda (iunit,buf,1,ibl)
      len = min(511-iwordb,l)
      call fmove (a,buf(iwordb+1),len)
      call putda (iunit,buf,1,ibl)
      ibl=ibl+1
      l=l-len
      ia=ia+len
500   if (l.lt.511) goto 600
      call fmove (a(ia),buf,511)
      call putda (iunit,buf,1,ibl)
      ibl=ibl+1
      l=l-511
      ia=ia+511
      goto 500
600   if (l.eq.0) return
      call getda (iunit,buf,1,ibl)
      call fmove (a(ia),buf,l)
      call putda (iunit,buf,1,ibl)
      return
      end
      subroutine putda (iunit,buf,lbl,ibl)
      implicit REAL  (a-h,o-z)
      dimension buf(511)
      call search (ibl,iunit)
      do 20 i=1,lbl
20    call put (buf,511,iunit)
      return
      end
      subroutine getda (iunit,buf,lbl,ibl)
      implicit REAL  (a-h,o-z)
      dimension buf(511)
      call search (ibl,iunit)
      do 10 i=1,lbl
      call find (iunit)
10    call get (buf,nw)
      return
      end
      subroutine readg (a,ll,iword,iunit)
      implicit REAL  (a-h,o-z)
c...   read l words into a from offset iword on file iunit
      dimension buf(512),ibuf(512),a(*)
      equivalence (buf(1),ibuf(1))
      l=ll
      ia=1
      ibl = iword/511
      iwordb = iword-ibl*511
      ibl = ibl+1
      if (iwordb.eq.0) goto 100
c..   first word not on block boundary
      call getda (iunit,buf,1,ibl)
      ibl=ibl+1
      len = min(511-iwordb,l)
      call fmove (buf(iwordb+1),a,len)
      l=l-len
      ia=ia+len
100   if (l.lt.511) goto 200
      call getda (iunit,buf,1,ibl)
      ibl=ibl+1
      l=l-511
      call fmove (buf,a(ia),511)
      ia=ia+511
      goto 100
200   if (l.eq.0) return
      call getda (iunit,buf,1,ibl)
      call fmove (buf,a(ia),l)
      return
      end
      subroutine intou2 (q,len)
      implicit REAL  (a-h,o-z)
INCLUDE(common/count)
      common /intbu2/ iblf,intpos,intfil,intmod,intnw,intb2,g(511)
      dimension q(*)
      iou2=iou2+1
      lenou2=lenou2+len
      if(len.le.0) return
      iq=1
      left=len
      if (intpos.gt.0) goto 10
      intmod=-1
      call search (iblf,intfil)
      intpos=1
10    lnow=min(left,512-intpos)
      call fmove (q(iq),g(intpos),lnow)
      iq=iq+lnow
      intpos=intpos+lnow
      left=left-lnow
      if (intpos.lt.512) return
      call put (g,511,intfil)
      intpos=1
      if (left.gt.0) goto 10
      return
      end
      subroutine inten2
      implicit REAL  (a-h,o-z)
INCLUDE(common/count)
      common /intbu2/ iblf,intpos,intfil,intmod,intnw,intb2,g(511)
      if (intmod.ge.0) goto 20
      if (intpos.gt.1) call put (g,511,intfil)
      call put (g,0,intfil)
      intmod=0
      intpos=1
      return
20    if (intmod.eq.1) call get(g,intnw)
      return
      end
      subroutine intin2 (q,len)
      implicit REAL  (a-h,o-z)
INCLUDE(common/count)
      common /intbu2/ iblf,intpos,intfil,intmod,intnw,intb2,g(511)
      dimension q(*)
cjk
cjk      write(6,57) intpos,intfil,intmod,intnw
cjk57    format('enter intin2',4i8)
      inn2=inn2+1
      lenin2=lenin2+len
      if(len.le.0) return
      iq=1
      left=len
      if (intpos.gt.0) goto 30
      intpos=1
      call search (iblf,intfil)
      call find (intfil)
      call get (g,intnw)
      if (intnw.eq.0) call caserr('intin2: end of file reached')
      intnw=intnw+1
      if (intmod.eq.1) call find (intfil)
30    lnow=min(left,intnw-intpos)
cjk
cjk      write(6,*) iq,intpos,left,lnow,intnw,len
      call fmove (g(intpos),q(iq),lnow)
      iq = iq + lnow
      intpos=intpos+lnow
      left=left-lnow
      if (intpos.lt.intnw) goto 50
      if (intmod.eq.0) call find(intfil)
      call get(g,intnw)
cjk40    if (intnw.eq.0) call caserr('intin2: end of file reached')
      intnw=intnw+1
      if (intmod.eq.1) call find(intfil)
      intpos=1
      if (left.gt.0) go to 30
50    continue
      return
      end
      subroutine ver_mcscfb(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mcscfb.m,v $
     +     "/
      data revision /"$Revision: 5957 $"/
      data date /"$Date: 2009-05-06 17:18:58 +0200 (Wed, 06 May 2009) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
