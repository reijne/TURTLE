c
c...   stash => stash_nbo  (occurs also in secd_parallel)
c      jvl,2003
c*****************************************************************************
c
c
c
c                  n  b  o     p  r  o  g  r  a  m
c
c                   (system independent routines)
c
c
c
c              last program modification:  october 22, 1991
c
c
c        !!! cray compilation requires 64 bit (-i64) integers !!!
c        (see, in particular, sr jobopt, sr nbopen, and sr debyte)
c
c  GAMESS-UK version based on code from CCL FTP site, August 1998
c  Only implemented for RHF and only on some platforms
c
c*****************************************************************************
c
c  main subroutine:
c
c      subroutine nbo(core,nboopt,memory)
c
c  job initialization routines: (called by sr nbo)
c
c      subroutine nboset(nboopt)
c      subroutine NBOjobopt(nboopt)
c      subroutine nbodim(memory)
c
c  nao/nbo/nlmo formation routines: (called by sr nbo)
c
c      subroutine naodrv(dm,t,a)
c      subroutine naosim(dm,t,a)
c      subroutine dmnao(dm,t,a)
c      subroutine dmsim(dm,t,a)
c      subroutine nbodrv(dm,t,a,memory)
c
c  routines called by the nao drivers:
c
c      subroutine simtrm(a,s,v,ndim,n,iwmulp,iwcubf)
c      subroutine mulana(bs,vmayer,bmayer,iwmulp,iwcubf)
c      subroutine dfgorb(renorm,dm,t,itran,iwcubf,itopt,lfnpr)
c      subroutine nao(t,s,occ,blk,sblk,eval,c,evect,eval2,listao,nblock)
c      subroutine naoanl(dm,spnao,bindex,bindt,bmo,ovpop,f,enao)
c      subroutine frmtmo(t,tmo,c,scr,index,iflg)
c
c  routines called by sr nao:
c
c      subroutine loadav(listao,nl,m,s,ndim,a,b,mxaolm)
c      subroutine atdiag(n,a,b,eval,c)
c      subroutine setbas(lstocc,lstemt,nocc,nemt,iat,l,nl,nf,ndim)
c      subroutine newwts(s,t,wt)
c      subroutine worth(s,t,blk,list,ndim,nbas,n,occ,eval,bigblk)
c      subroutine NBOshmdt(t,s,ndim,nbas,nocc,lstocc,nemt,lstemt,sblk)
c      subroutine newryd(t,s,tpnao,dmblk,sblk,evect,occ,eval,eval2,
c     +                       list,irpnao)
c      subroutine rydiag(t,s,tpnao,dmblk,sblk,occ,eval,evect,eval2,
c     +                    iorb,nc,nm,nstart,nrydc,larc,list,irpnao)
c      subroutine rydsel(lstemt,nemt,nsel1,list1,nsel2,list2,wt)
c      subroutine rediag(dm,t,tpnao,eval,blk,c,irank,irpnao)
c      subroutine redblk(t,tpnao,il,dm,blk,eval,c,nf,iorb,nc,irank,irpnao)
c
c  routines called by the nbo/nlmo drivers:
c
c      subroutine nathyb(dm,t,guide,bndocc,pol,q,v,blk,c,eval,borb,
c     +                                       p,ta,hyb,va,vb,topo)
c      subroutine chsdrv(dm,t,guide,bndocc,pol,q,v,blk,c,eval,borb,
c     +                                       p,ta,hyb,va,vb,topo)
c      subroutine choose(dm,t,guide,bndocc,pol,q,v,blk,c,eval,borb,
c     +                                  p,ta,hyb,va,vb,topo,iflg)
c      subroutine srtnbo(t,bndocc)
c      subroutine xcited(dm,t,hyb,thyb,s,occ,scr,iscr)
c      subroutine anlyze(t,bndocc,hyb,hycoef,thyb)
c      subroutine htype(hyb,ltyp,mxao,nh,coef,pct,nl,isgn)
c      subroutine frmhyb(hyb,thyb,coef,hycoef,kl,ku,nhyb)
c      subroutine hybdir(bndocc,atcoor,thyb,tbnd,scr)
c      subroutine hybcmp(xyz,pcent,ihyb,jctr,hyb)
c      subroutine fndmol(iatoms)
c      subroutine nbocla(bndocc,accthr)
c      subroutine fnboan(bndocc,f,molnbo)
c      subroutine nbosum(f,bndocc,list,lista,scr)
c      subroutine getdel(ibo,occ,thr1,thr2,nl,list,del,deloc,iflg)
c      subroutine dlcstr(ibo,il,nl,list,ml,istr)
c      subroutine nlmo(n,a,eval,evec,tsym,reson,nocc,ialarm)
c      subroutine lmoanl(t,s,reson,occ,ts,border,owbord,atlmo,siab,nocc,nab)
c      subroutine dipanl(dm,t,c,tnbo,dx,dy,dz,scr,index)
c      subroutine dipele(dxyz,c,t,scr,eta,nocc,index)
c      subroutine NBOdipnuc(dx,dy,dz,atcoor,eta,nocc)
c
c  routines called by sr nathyb, sr choose:
c
c      subroutine core(dm,t,borb,pol,q,hyb,bndocc,ibd,detail,lfnpr)
c      function iwprj(nctr)
c      subroutine deplet(dm,t,q,pol,borb,bndocc,nbd)
c      subroutine NBOload(dm,iat1,iat2,iat3,blk,nb)
c      subroutine prjexp(borb,iat1,iat2,iat3,q,p,pk,hyb,va,vb,hybexp)
c      subroutine stash_nbo(borb,ibd,iat1,iat2,iat3,pol,q,hyb)
c      subroutine orthyb(q,s,ta,eval,c,ialarm,iflg)
c      subroutine frmprj(p,ia,q,nk,pk,vk,pi)
c      subroutine augmnt(p,blk,c,eval,dm,ta,borb,v,larc,ia,nocc,norb)
c      subroutine repol(dm,q,pol,blk,eval,c,nbd)
c      subroutine formt(t,q,pol)
c      subroutine cycles(iter,thresh,guide,bndocc,topo,icont)
c
c  routines called by sr nlmo:
c
c      subroutine symuni(tsym,a,cos,sin,ovlp,blk,eval,nrot,
c     +           niuniq,njuniq,ilist,jlist,noff,ioff,joff,ndim)
c      subroutine NBOsymort(s,t,blk,ndim,n,eval)
c
c  nbo energetic analysis routines:
c
c      subroutine nboean(a,memory,nboopt,idone)
c      subroutine nbodel(a,memory,idone)
c      subroutine delete(f,trf,ndim,idel,len,itype,ndel,ntrunc,done,
c     +                  ispin)
c      subroutine newdm(dm,u,eig,ndim,idel,len,ndel,itype,nmoocc,ispin)
c      subroutine rnkeig(rank,eig,n,ndim,arcrnk)
c      subroutine simltr(n,ndim,f,u,r,s,kntrol)
c
c  nbo direct access file (daf) routines:
c
c      subroutine nbfile(new,error)
c      subroutine nbopen(new,error)
c      subroutine nbwrit(ix,nx,idar)
c      subroutine nbread(ix,nx,idar)
c      subroutine nbclos(seq)
c      subroutine nbinqr(idar)
c
c      subroutine fetitl(title)
c      subroutine fee0(edel,etot)
c      subroutine sve0(edel)
c      subroutine fecoor(atcoor)
c      subroutine fesraw(s)
c      subroutine fedraw(dm,scr)
c      subroutine fefao(f,iwfock)
c      subroutine feaomo(t,it)
c      subroutine fedxyz(dxyz,i)
c      subroutine svnbo(t,occ,iscr)
c      subroutine fenbo(t,occ,iscr,nelec)
c      subroutine fetnbo(t)
c      subroutine svpnao(t)
c      subroutine fepnao(t)
c      subroutine svsnao(s)
c      subroutine fesnao(s)
c      subroutine svtnab(t)
c      subroutine fetnab(t)
c      subroutine svtlmo(t)
c      subroutine fetlmo(t)
c      subroutine svtnho(t)
c      subroutine fetnho(t)
c      subroutine svppao(dm)
c      subroutine feppao(dm)
c      subroutine svtnao(t)
c      subroutine fetnao(t)
c      subroutine svnlmo(t)
c      subroutine fenlmo(t)
c      subroutine svdnao(dm)
c      subroutine fednao(dm)
c      subroutine svfnbo(f)
c      subroutine fefnbo(f)
c      subroutine svnewd(dm)
c      subroutine fenewd(dm)
c      subroutine feinfo(icore,iswean)
c      subroutine febas(nshell,nexp,iscr)
c
c  free format input routines:
c
c      subroutine strtin(lfnin)
c      subroutine rdcrd
c      subroutine ifld(int,error)
c      subroutine rfld(real,error)
c      subroutine hfld(keywd,leng,endd)
c      subroutine fndfld
c      function equal(ia,ib,l)
c
c  other system-independent i/o routines:
c
c      subroutine geninp(newdaf)
c      subroutine nboinp(nboopt,idone)
c      subroutine corinp(iess,icor)
c      subroutine chsinp(iess,ichs)
c      subroutine delinp(nboopt,idone)
c
c      subroutine rdcore(jcore)
c      subroutine wrppna(t,occ,iflg)
c      subroutine rdppna(t,occ,iflg)
c      subroutine wrtnao(t,iflg)
c      subroutine rdtnao(dm,t,scr,iflg)
c      subroutine wrtnab(t,iflg)
c      subroutine rdtnab(t,dm,bndocc,scr,iflg)
c      subroutine wrtnbo(t,bndocc,iflg)
c      subroutine wrnlmo(t,dm,iflg)
c      subroutine wrbas(scr,iscr,lfn)
c      subroutine wrarc(scr,iscr,lfn)
c
c      subroutine aout(a,mr,nr,nc,title,index,iflg)
c      subroutine aprint(a,mr,nr,nc,title,index,mcol)
c      subroutine awrite(a,mr,nr,nc,title,lfn)
c      subroutine aread(a,mr,nr,nc,job,lfn,error)
c      subroutine altout(a,mr,mc,nr,nc)
c      subroutine keypar(string,len,iflg,lfn,read,error)
c      function ioinqr(iflg)
c      subroutine lblao
c      subroutine lblnao
c      subroutine lblnbo
c      subroutine lblnho(inho,inbo,ictr,nctr)
c
c  general utility routines:
c
c      subroutine angles(x,y,z,theta,phi)
c      function bdfind(iat,jat)
c      subroutine chem(nat,natoms,lista,nl,istr)
c      subroutine consol(aut,alt,ndim,n)
c      subroutine convin(ij,len,ik,error)
c      subroutine convrt(n,nc1,nc2)
c      subroutine copy(a,b,ndim,nr,nc)
c      subroutine cortbl(iat,icore,iecp)
c      subroutine debyte(i,ibyte)
c      subroutine halt(word)
c      subroutine idigit(kint,ik,nd,maxd)
c      function ihtyp(ibo,jbo)
c      subroutine NBOjacobi(n,a,eivu,eivr,ndim,nvdim,icontr)
c      subroutine limtrn(t,m,a,b,ndim,nbas,ncdim,nc,iopt)
c      subroutine matmlt(a,b,v,ndim,n)
c      subroutine matml2(a,b,v,ndim,n)
c      function nameat(iz)
c      subroutine normlz(a,s,m,n)
c      subroutine NBOorder(rank,list,n,ndim,arcrnk)
c      subroutine NBOpack(t,ndim,nbas,l2)
c      subroutine rank(eig,n,ndim,arcrnk)
c      subroutine simtrn(a,t,v,ndim,n)
c      subroutine simtrs(a,s,v,ndim,n)
c      subroutine NBOtransp(a,ndim,n)
c      subroutine NBOunpack(t,ndim,nbas,l2)
c      subroutine valtbl(iat,ival)
c      function veclen(x,n,ndim)
c
c      subroutine lineq(a,x,b,scr,n,m,ndim,mdim,zertol,eps,maxit,lfnpr,
c     +                 ierr)
c      subroutine factor(a,w,d,ipivot,n,ndim,zertol,iflag)
c      subroutine fndsol(a,x,b,w,r,e,ipivot,n,ndim,eps,maxit,lfnpr,ierr)
c      subroutine subst(x,w,b,ipivot,n,ndim)
c
c*****************************************************************************
      subroutine nbo(core,memory,nboopt)
c*****************************************************************************
c
c  input:
c     core       core memory to be dynamically allocated for storage needs.
c     memory     the number of real*8 words available in core.
c     nboopt(10) list of nbo options as summarized below:
c
c     nboopt(1)  = -2       do nothing
c                = -1       natural population analysis (npa) only
c                =  0       perform npa/nbo/nlmo analyses
c                =  1       perform npa/nbo/nlmo analyses, don't read keywords
c                =  2       perform one fock matrix deletion, forming new dm
c                =  3       evaluate and print the energy change from deletion
c
c     nboopt(2)  =  0       scf density
c                =  1       mp first order density
c                =  3       mp2 density
c                =  4       mp3 density
c                =  5       mp4 density
c                =  6       ci one-particle density
c                =  7       ci density
c                =  8       qci/cc density
c                =  9       density correct to second order
c
c     nboopt(3)  =  1       transform dipole moment matrices to nbo/nlmo bases
c
c     nboopt(4)  =  1       allow strongly resonant lewis structures
c                           (force the resonance keyword)
c
c     nboopt(5)  =  1       spin-annihilated uhf (auhf) wavefunction
c
c     nboopt(6-9)           unused
c
c     nboopt(10) =  0       general version of the nbo program (gennbo)
c                =  1       ampac version
c                =  6       gamess version
c                =  7       hondo version
c                =  8x      gaussian 8x version
c------------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical newdaf,error,seq
c
c  nbo common blocks:
c
cold      parameter(maxatm = 99,maxbas = 500)
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbao/lctr(maxbas),lang(maxbas)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension core(memory),nboopt(10)
c
c  if nboopt(1).eq.-2, then no nbo analysis was requested:
c
      if(nboopt(1).eq.-2) return
c
c  set default options:
c
      call nboset(nboopt)
c
c  if this is the general version of the program, read the $gennbo keylist:
c
       if(nboopt(10).eq.0) then
        call geninp(newdaf)
      else
        newdaf = .true.
      end if
c
c  search the input file for the $nbo keylist:
c
      call nboinp(nboopt,idone)
      if(idone.eq.1) return
c
c  read in job options from the $nbo keylist:
c
      call NBOjobopt(nboopt)
c
c  check filename and open sequential files:
c
      call nbfile(newdaf,error)
      if(error) return
c
c  open the nbo direct access file:
c
      call nbopen(newdaf,error)
      if(error) then
        write(lfnpr,900)
        return
      end if
c
c  fetch atoms, basis, and wave function information:
c
      call feaoin(core,core,nboopt)
      if(complx) return
c
c  write the job title to the output file:
c
      call fetitl(core)
      write(lfnpr,910) (core(i),i=1,8)
c
c  set up dimensioning information and determine if enough space is available:
c
      call nbodim(memory)
c
c  set up basic storage:
c
c  core(ndm) :  ndim by ndim matrix to store density matrix
c  core(nt)  :  ndim by ndim matrix to hold overlap or transformation matrices
c  core(nscr):  scratch storage, dynamically allocated according needs
c
      n2   = ndim*ndim
      ndm  = 1
      nt   = ndm + n2
      nscr = nt  + n2
      mem  = memory - nscr + 1
c
c  read in input overlap and density matrices, ao basis:
c
      alpha = .false.
      beta  = .false.
      ispin = 0
      call fedraw(core(ndm),core(nscr))
c
c  simulate the natural population analysis if the input basis is orthogonal:
c
      if(ortho) then
        call naosim(core(ndm),core(nt),core(nscr))
c
c  load the overlap matrix into core(nt) and perform the natural population
c  analysis:
c
      else
        call fesraw(core(nt))
        call naodrv(core(ndm),core(nt),core(nscr))
      end if
c
c  note: core(ndm) now contains the total density matrix in the nao basis
c        and core(nt) contains the ao to nao transformation
c
c  perform closed shell nbo analysis:
c
      if(.not.open) then
        call nbodrv(core(ndm),core(nt),core(nscr),mem)
      else
c
c  perform open shell nbo analysis:
c
c  first, analyze alpha density matrix:
c
        alpha = .true.
        beta  = .false.
        ispin = 2
        if(ortho) then
          call dmsim(core(ndm),core(nt),core(nscr))
        else
          call dmnao(core(ndm),core(nt),core(nscr))
        end if
        call nbodrv(core(ndm),core(nt),core(nscr),mem)
c
c  now, analyze beta density matrix:
c
        alpha = .false.
        beta  = .true.
        ispin = -2
        if(ortho) then
          call dmsim(core(ndm),core(nt),core(nscr))
        else
          call dmnao(core(ndm),core(nt),core(nscr))
        end if
        call nbodrv(core(ndm),core(nt),core(nscr),mem)
      end if
c
c  close the nbo direct access file and other external files:
c
      seq = .true.
      call nbclos(seq)
c
      return
c
  900 format(/1x,'nbo direct access file could not be opened.  nbo ',
     + 'program aborted.')
  910 format(/1x,'job title: ',8a8)
      end
c*****************************************************************************
c
c  job initialization routines: (called by sr nbo)
c
c      subroutine nboset(nboopt)
c      subroutine NBOjobopt(nboopt)
c      subroutine nbodim(memory)
c
c*****************************************************************************
      subroutine nboset(nboopt)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension nboopt(10)
c
      parameter (maxatm = 750, maxbas = 4096)
      parameter(maxfil = 40)
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nblbl/nlew,nval,lbl(10,maxbas,4)
      common/nbname/filenm,nfile,ifile(maxfil)
      character*80 filenm
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
      data tenth,half/0.1d0,0.5d0/
c
c  set default job options:  (modifications to these defaults
c  should not be made here but later in this subroutine)
c
c  use the bond-order matrix, not the occupation matrix (expectation
c  values of the density operator)
c
      iwdm   =  1
      iw3c   =  0
      iwapol =  0
      iwhybs =  0
      iwpnao =  0
      iwtnao =  0
      iwtnab =  0
      iwtnbo =  0
c
c  use the fock matrix, if there is one:
c
      call sectst(isect(42),iwfock)
c
c  set to the desired print level + 10:
c
      iprint = 12
      ipseud =  0
      iwdetl =  0
      iwmulp =  0
      ichoos =  0
      kopt   =  0
      jcore  =  0
      iwcubf = 0
      open   = .false.
      ortho  = .false.
      uhf    = .false.
      auhf   = .false.
      rohf   = .false.
      ci     = .false.
      mcscf  = .false.
      complx = .false.
      do 10 i = 1,60
        jprint(i) = 0
   10 continue
c
      lfnao  =  31
      lfnpna =  32
      lfnnao =  33
      lfnpnh =  34
      lfnnho =  35
      lfnpnb =  36
      lfnnbo =  37
      lfnpnl =  38
      lfnnlm =  39
      lfnmo  =  40
      lfndm  =  41
      lfnnab =  42
      lfnppa =  43
      lfnarc =  47
c
c  set positive in routine NBOjobopt if chosen by the user:
c
      lfndaf = -48
      lfndef =  49
c
c  setting nval negative indicates that this variable has not
c  been determined yet:
c
      nval   = -1
c
c  initialize the character string used to create filenames:
c
      filenm(1:4) = 'file'
      do 50 i = 5,80
        filenm(i:i) = char(32)
   50 continue
c
c  that some thresholds are .lt.0 indicates that these variables have not
c  been set by the user:
c
      thrset =  -1.9d0
      prjset =  -0.2d0
      accthr =  -tenth
      crtset =   1.999d0
      e2thr  =  -half
      athr   =  -1.000d0
      pthr   = -25.000d0
      ethr   =  -0.100d0
      dthr   =  -0.020d0
      dlthr  =  -1.000d0
      chsthr =  -0.100d0
c
c  set job options according to nboopt:
c
c  skip the computation of the nbos?
c
      if(nboopt(1).eq.-1) jprint(1) = 1
c
c  turn off $choose and $core keylists if $nbo keylist is not to
c  be read:
c
      if(nboopt(1).eq.1) ichoos = -1
      if(nboopt(1).eq.1) jcore  = -1
c
c  force dipole analysis?
c
      if(nboopt(3).ne.0) then
        jprint(46) = 1
      end if
c
c  force resonance keyword?
c
      if(nboopt(4).ne.0) jprint(14) = 1
c
c  program version:
c
      jprint(2) = nboopt(10)
c
      return
      end
c******************************************************************************
      subroutine NBOjobopt(nboopt)
c******************************************************************************
      implicit real*8 (a-h,o-z)
      logical error,end,equal,nextwd,read,equal2
      dimension nboopt(10),inttmp(80)
c
      parameter(keylen = 9)
      parameter(maxfil = 40)
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nbname/filenm,nfile,ifile(maxfil)
      character*80 filenm
c
      dimension keywd(keylen),k3cbnd(6),kepert(6),klfnpr(5),kthrsh(6),
     * kdetl(6),kmula(5),kmulor(6),kprjth(6),knbnlm(7),
     * kaopre(6),knlmo(4),knaomo(5),knbomo(5),knobnd(6),knpa(3),
     * kskipb(6),krpnao(5),kbndid(6),knlmmo(6),kreson(5),kppnao(7),
     * kaonao(5),knanbo(6),kaonbo(5),kaonlm(6),kfnbo(4),kfnlmo(5),
     * kdmnbo(5),kdmnlm(6),kprint(5),knanlm(7),kspnao(5),kspnho(5),
     * kspnbo(5),kaoinf(6),kaopnb(6),kaomo(4),knaonh(6),knhnlm(7),
     * kaonho(5),kfnho(4),kaopnh(6),kfnao(4),knhonb(6),kspnlm(6),
     * knrt(3),kdmnho(5),kdmnao(5),kplot(4),kaopnl(7),kdiao(4),
     * kbend(4),knhomo(5),ksao(3),kfao(3),kdmao(4),kboao(4),kdinlm(6),
     * knbosm(6),knbo(3),kdipol(6),kdinao(5),kdinho(5),kdinbo(5),
     * knbdaf(6),karchv(7),kfile(4),kpolar(6),knrtop(6),knrtrf(6),
     * kchsth(6),knrtdt(6),knrtth(6)
c
      dimension kalt(4),kbfgs(4),kpowel(6),ksap(3)
c
      data k3cbnd/1h3,1hc,1hb,1ho,1hn,1hd/,klfnpr/1hl,1hf,1hn,1hp,1hr/,
     * kthrsh/1ht,1hh,1hr,1he,1hs,1hh/,kepert/1he,1h2,1hp,1he,1hr,1ht/,
     * kplot/1hp,1hl,1ho,1ht/,kdetl/1hd,1he,1ht,1ha,1hi,1hl/,
     * kmula/1hm,1hu,1hl,1ha,1ht/,kmulor/1hm,1hu,1hl,1ho,1hr,1hb/,
     * kprjth/1hp,1hr,1hj,1ht,1hh,1hr/,kaopre/1ha,1ho,1hp,1hn,1ha,1ho/,
     * knlmo/1hn,1hl,1hm,1ho/,knpa/1hn,1hp,1ha/,knbo/1hn,1hb,1ho/,
     * knaomo/1hn,1ha,1ho,1hm,1ho/,knbomo/1hn,1hb,1ho,1hm,1ho/,
     * knobnd/1hn,1ho,1hb,1ho,1hn,1hd/,kskipb/1hs,1hk,1hi,1hp,1hb,1ho/,
     * krpnao/1hr,1hp,1hn,1ha,1ho/,kbndid/1hb,1hn,1hd,1hi,1hd,1hx/,
     * knlmmo/1hn,1hl,1hm,1ho,1hm,1ho/,kreson/1hr,1he,1hs,1ho,1hn/,
     * kppnao/1hp,1ha,1ho,1hp,1hn,1ha,1ho/,kaonao/1ha,1ho,1hn,1ha,1ho/,
     * knanbo/1hn,1ha,1ho,1hn,1hb,1ho/,kaonbo/1ha,1ho,1hn,1hb,1ho/
c
      data kaonlm/1ha,1ho,1hn,1hl,1hm,1ho/,kfnbo/1hf,1hn,1hb,1ho/,
     * kfnlmo/1hf,1hn,1hl,1hm,1ho/,kprint/1hp,1hr,1hi,1hn,1ht/,
     * kdmnbo/1hd,1hm,1hn,1hb,1ho/,kdmnlm/1hd,1hm,1hn,1hl,1hm,1ho/,
     * knanlm/1hn,1ha,1ho,1hn,1hl,1hm,1ho/,kaomo/1ha,1ho,1hm,1ho/,
     * kspnao/1hs,1hp,1hn,1ha,1ho/,kspnho/1hs,1hp,1hn,1hh,1ho/,
     * kspnbo/1hs,1hp,1hn,1hb,1ho/,kfnao/1hf,1hn,1ha,1ho/,
     * kaoinf/1ha,1ho,1hi,1hn,1hf,1ho/,kaopnb/1ha,1ho,1hp,1hn,1hb,1ho/,
     * kaonho/1ha,1ho,1hn,1hh,1ho/,kfnho/1hf,1hn,1hh,1ho/,
     * kaopnh/1ha,1ho,1hp,1hn,1hh,1ho/,knrt/1hn,1hr,1ht/,
     * knbnlm/1hn,1hb,1ho,1hn,1hl,1hm,1ho/,kdiao/1hd,1hi,1ha,1ho/,
     * kdmnho/1hd,1hm,1hn,1hh,1ho/,kdmnao/1hd,1hm,1hn,1ha,1ho/,
     * kbend/1hb,1he,1hn,1hd/,knbosm/1hn,1hb,1ho,1hs,1hu,1hm/,
     * knhomo/1hn,1hh,1ho,1hm,1ho/,ksao/1hs,1ha,1ho/,kfao/1hf,1ha,1ho/
c
      data kdmao/1hd,1hm,1ha,1ho/,kboao/1hb,1ho,1ha,1ho/,
     * kdipol/1hd,1hi,1hp,1ho,1hl,1he/,knaonh/1hn,1ha,1ho,1hn,1hh,1ho/,
     * knhnlm/1hn,1hh,1ho,1hn,1hl,1hm,1ho/,kdinao/1hd,1hi,1hn,1ha,1ho/,
     * knhonb/1hn,1hh,1ho,1hn,1hb,1ho/,kspnlm/1hs,1hp,1hn,1hl,1hm,1ho/,
     * kaopnl/1ha,1ho,1hp,1hn,1hl,1hm,1ho/,kdinho/1hd,1hi,1hn,1hh,1ho/,
     * kdinbo/1hd,1hi,1hn,1hb,1ho/,kdinlm/1hd,1hi,1hn,1hl,1hm,1ho/,
     * knbdaf/1hn,1hb,1ho,1hd,1ha,1hf/,
     * karchv/1ha,1hr,1hc,1hh,1hi,1hv,1he/,kfile/1hf,1hi,1hl,1he/,
     * kpolar/1ha,1hp,1ho,1hl,1ha,1hr/,knrtop/1hn,1hr,1ht,1ho,1hp,1ht/,
     * knrtrf/1hn,1hr,1ht,1hr,1he,1hf/,kchsth/1hc,1hh,1hs,1ht,1hh,1hr/,
     * knrtdt/1hn,1hr,1ht,1hd,1ht,1hl/,
     * knrtth/1hn,1hr,1ht,1ht,1hh,1hr/
c
      data kalt/1h$,1he,1hn,1hd/,kbfgs/1hb,1hf,1hg,1hs/,
     * kpowel/1hp,1ho,1hw,1he,1hl,1hl/,ksap/1hs,1ha,1hp/
c
      data zero,one/0.0d0,1.0d0/
      data ifull,ival,ilew/4hfull,3hval,3hlew/
      data iprnt,iwrit,iread/4hprnt,4hwrit,4hread/
      data ia,ib,ip/1ha,1hb,1hp/
c
c  read in job options, in a keyword directed manner:
c
      numopt = 0
      lennm  = 0
      if(nboopt(1).eq.1) goto 4500
c
c  begin loop to identify keyword "keywd":
c
      nextwd = .true.
  100 leng = keylen
      if(nextwd) call hfld(keywd,leng,end)
      nextwd = .true.
      if((leng.eq.0).or.end) go to 4500
      if(equal2(keywd,kalt,4)) go to 4500
      numopt = numopt + 1
c
c  keyword: 3cbond -- search for three-center bonds
c   (default is to search only for one- and two-center nbos)
      if(.not.equal(keywd,k3cbnd,6)) go to 500
        iw3c = 1
        go to 100
c  keyword: lfnpr -- specify output lfn
  500 if(.not.equal(keywd,klfnpr,5)) go to 510
      call ifld(lfnpr,error)
        if(error) call halt('lfnpr')
        go to 100
c  keyword: thresh -- specify fixed occupancy threshold for nbo search
  510 if(.not.equal(keywd,kthrsh,6)) go to 540
      call rfld(thrset,error)
        if(error) call halt('thresh')
        go to 100
c  keyword: detail -- print details of nbo search procedure
  540 if(.not.equal(keywd,kdetl,6))  go to 550
        iwdetl = 1
        go to 100
c  keyword: mulat -- print mulliken populations by atom
  550 if(.not.equal(keywd,kmula,5))  go to 560
        iwmulp = 1
        go to 100
c  keyword: mulorb -- print mulliken populations by orbital and atom
  560 if(.not.equal(keywd,kmulor,6)) go to 580
        iwmulp = 2
        go to 100
c  keyword: prjthr -- user sets value of projection threshold for nbo search
c           for rejecting linearly dependent hybrids
  580 if(.not.equal(keywd,kprjth,6)) go to 610
      call rfld(prjset,error)
        if(error) call halt('prjthr')
        go to 100
c  keyword: fnbo -- print nbo fock matrix
  610 if(.not.equal(keywd,kfnbo,4)) go to 620
        jprint(37) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(37),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(37).eq.ival) jprint(37) = ifull
        end if
        go to 100
c  keyword: aopnao -- output raw ao to pnao transformation
  620 if(.not.equal(keywd,kaopre,6)) go to 640
        jprint(44) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(44),lfnpna,read,error)
          if(error) nextwd = .false.
          if(jprint(44).eq.ival) jprint(44) = ifull
          if(jprint(44).eq.ilew) jprint(44) = ifull
        end if
        go to 100
c  keyword: nlmomo -- compute and print nlmo to mo transf.
  640 if(.not.equal(keywd,knlmmo,6)) go to 650
        jprint(13) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(13),lfndef,read,error)
          if(error) nextwd = .false.
        end if
        go to 100
c  keyword: nlmo -- compute and print nlmos
  650 if(.not.equal(keywd,knlmo,4))  go to 660
        if(leng.ne.4) go to 660
        jprint(8) = 1
        go to 100
c  keyword: naomo -- compute and print nao to mo transf.
  660 if(.not.equal(keywd,knaomo,5)) go to 670
        jprint(9) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(9),lfndef,read,error)
          if(error) nextwd = .false.
        end if
        go to 100
c  keyword: nbomo -- compute and print nbo to mo transf.
  670 if(.not.equal(keywd,knbomo,5)) go to 680
        jprint(45) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(45),lfndef,read,error)
          if(error) nextwd = .false.
        end if
        go to 100
c  keyword: nobond -- compute only one-center nbos
  680 if(.not.equal(keywd,knobnd,6)) go to 690
        jprint(10) = 1
        go to 100
c  keyword: skipbo -- skip nbo procedure
  690 if(.not.equal(keywd,kskipb,6)) go to 700
        jprint(1) = 1
        go to 100
c  keyword: rpnao -- compute revised pure ao to pnao transf.
  700 if(.not.equal(keywd,krpnao,5)) go to 710
        jprint(11) = 1
        go to 100
c  keyword: bndidx -- print bond indices
  710 if(.not.equal(keywd,kbndid,6)) go to 730
        jprint(12) = 1
        go to 100
c  keyword: resonance -- allow strongly "non-lewis" nbo occupancies
c   (overrides automatic shutdown of nbo procedure in strongly
c    delocalized cases)
  730 if(.not.equal(keywd,kreson,5)) go to 740
        jprint(14) = 1
        go to 100
c  keyword: paopnao -- i/o with pao to pnao transformation
  740 if(.not.equal(keywd,kppnao,7)) go to 750
        iwpnao = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .true.
          call keypar(keywd,leng,iwpnao,lfnppa,read,error)
          if(error) nextwd = .false.
          if(iwpnao.eq.ival) iwpnao = ifull
          if(iwpnao.eq.ilew) iwpnao = ifull
        end if
        go to 100
c  keyword: aonao -- i/o with ao to nao transformation
  750 if(.not.equal(keywd,kaonao,5)) go to 760
        iwtnao = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .true.
          call keypar(keywd,leng,iwtnao,lfnnao,read,error)
          if(error) nextwd = .false.
          if(iwtnao.eq.ival) iwtnao = ifull
          if(iwtnao.eq.ilew) iwtnao = ifull
        end if
        go to 100
c  keyword: naonbo -- i/o with nao to nbo transformation
  760 if(.not.equal(keywd,knanbo,6)) go to 770
        iwtnab = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .true.
          call keypar(keywd,leng,iwtnab,lfnnab,read,error)
          if(error) nextwd = .false.
          if(iwtnab.eq.ival) iwtnab = ifull
        end if
        go to 100
c  keyword: aonbo -- output ao to nbo transf. information
  770 if(.not.equal(keywd,kaonbo,5)) go to 780
        iwtnbo = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,iwtnbo,lfnnbo,read,error)
          if(error) nextwd = .false.
          if(iwtnbo.eq.ival) iwtnbo = ifull
        end if
        go to 100
c  keyword: fnlmo -- print nlmo fock matrix
  780 if(.not.equal(keywd,kfnlmo,5)) go to 790
        jprint(15) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(15),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(15).eq.ival) jprint(15) = ifull
        end if
        go to 100
c  keyword: dmnbo -- print nbo density matrix
  790 if(.not.equal(keywd,kdmnbo,5)) go to 800
        jprint(16) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(16),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(16).eq.ival) jprint(16) = ifull
        end if
        go to 100
c  keyword: dmnlmo -- print nlmo density matrix
  800 if(.not.equal(keywd,kdmnlm,6)) go to 810
        jprint(17) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(17),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(17).eq.ival) jprint(17) = ifull
        end if
        go to 100
c  keyword: aonlmo -- compute and output ao to nlmo transf.
  810 if(.not.equal(keywd,kaonlm,6)) go to 820
        jprint(23) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(23),lfnnlm,read,error)
          if(error) nextwd = .false.
          if(jprint(23).eq.ival) jprint(23) = ifull
        end if
        go to 100
c  keyword: print -- read in print option level "iprint"
  820 if(.not.equal(keywd,kprint,5)) go to 830
        call ifld(iprint,error)
        if(error) call halt('print')
        go to 100
c  keyword: naonlmo -- print nao to nlmo transformation matrix
  830 if(.not.equal(keywd,knanlm,7)) go to 840
        jprint(18) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(18),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(18).eq.ival) jprint(18) = ifull
        end if
        go to 100
c  keyword: spnao -- print s-pnao overlap matrix
  840 if(.not.equal(keywd,kspnao,5)) go to 850
        jprint(19) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(19),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(19).eq.ival) jprint(19) = ifull
          if(jprint(19).eq.ilew) jprint(19) = ifull
        end if
        go to 100
c  keyword: spnho -- print s-pnho overlap matrix
  850 if(.not.equal(keywd,kspnho,5)) go to 860
        jprint(20) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(20),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(20).eq.ival) jprint(20) = ifull
          if(jprint(20).eq.ilew) jprint(20) = ifull
        end if
        go to 100
c  keyword: nhonlmo -- output the nho to nlmo transformation
  860 if(.not.equal(keywd,knhnlm,7)) go to 870
        jprint(24) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(24),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(24).eq.ival) jprint(24) = ifull
        end if
        go to 100
c  keyword: spnbo -- print s-pnbo overlap matrix
  870 if(.not.equal(keywd,kspnbo,5)) go to 880
        jprint(21) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(21),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(21).eq.ival) jprint(21) = ifull
        end if
        go to 100
c  keyword: aoinfo -- write basis set info
  880 if(.not.equal(keywd,kaoinf,6)) go to 910
        jprint(22) = lfnao
        call ifld(itemp,error)
        if(.not.error) jprint(22) = abs(itemp)
        go to 100
c  keyword: aopnbo -- write ao to pnbo transformation
  910 if(.not.equal(keywd,kaopnb,6)) go to 920
        jprint(25) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(25),lfnpnb,read,error)
          if(error) nextwd = .false.
          if(jprint(25).eq.ival) jprint(25) = ifull
        end if
        go to 100
c  keyword: aomo -- write ao to mo transformation
  920 if(.not.equal(keywd,kaomo,4)) go to 930
        jprint(26) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(26),lfnmo,read,error)
          if(error) nextwd = .false.
        end if
        go to 100
c  keyword: dmao -- write ao density matrix
  930 if(.not.equal(keywd,kdmao,4)) go to 940
        jprint(27) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(27),lfndm,read,error)
          if(error) nextwd = .false.
          if(jprint(27).eq.ival) jprint(27) = ifull
          if(jprint(27).eq.ilew) jprint(27) = ifull
        end if
        go to 100
c  keyword: aonho -- write ao to nho transformation
  940 if(.not.equal(keywd,kaonho,5)) go to 950
        jprint(28) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(28),lfnnho,read,error)
          if(error) nextwd = .false.
          if(jprint(28).eq.ival) jprint(28) = ifull
          if(jprint(28).eq.ilew) jprint(28) = ifull
        end if
        go to 100
c  keyword: fnho -- print nho fock matrix
  950 if(.not.equal(keywd,kfnho,4)) go to 960
        jprint(29) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(29),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(29).eq.ival) jprint(29) = ifull
          if(jprint(29).eq.ilew) jprint(29) = ifull
        end if
        go to 100
c  keyword: aopnho -- write ao to pnho transformation
  960 if(.not.equal(keywd,kaopnh,6)) go to 970
        jprint(30) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(30),lfnpnh,read,error)
          if(error) nextwd = .false.
          if(jprint(30).eq.ival) jprint(30) = ifull
          if(jprint(30).eq.ilew) jprint(30) = ifull
        end if
        go to 100
c  keyword: fnao -- print nao fock matrix
  970 if(.not.equal(keywd,kfnao,4)) go to 990
        jprint(31) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(31),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(31).eq.ival) jprint(31) = ifull
          if(jprint(31).eq.ilew) jprint(31) = ifull
        end if
        go to 100
c  keyword: naonho -- output the nao to nho transformation
  990 if(.not.equal(keywd,knaonh,6)) go to 1010
        jprint(33) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(33),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(33).eq.ival) jprint(33) = ifull
          if(jprint(33).eq.ilew) jprint(33) = ifull
        end if
        go to 100
c  keyword: dmnho -- print nho density matrix
 1010 if(.not.equal(keywd,kdmnho,5)) go to 1020
        jprint(34) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(34),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(34).eq.ival) jprint(34) = ifull
          if(jprint(34).eq.ilew) jprint(34) = ifull
        end if
        go to 100
c  keyword: dmnao -- print nao density matrix
 1020 if(.not.equal(keywd,kdmnao,5)) go to 1040
        jprint(35) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(35),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(35).eq.ival) jprint(35) = ifull
          if(jprint(35).eq.ilew) jprint(35) = ifull
        end if
        go to 100
c  keyword: bend -- print nho directionality and bond bending info
 1040 if(.not.equal(keywd,kbend,4)) go to 1050
        jprint(36) = 1
        call rfld(temp,error)
        if(error) go to 100
        athr = abs(temp)
        call rfld(temp,error)
        if(error) go to 100
        pthr = abs(temp)
        if(pthr.lt.one) pthr = one
        call rfld(temp,error)
        if(error) go to 100
        ethr = abs(temp)
        go to 100
c  keyword: nhomo -- compute and print nho to mo transf.
 1050 if(.not.equal(keywd,knhomo,5)) go to 1060
        jprint(38) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(38),lfndef,read,error)
          if(error) nextwd = .false.
        end if
        go to 100
c  keyword: sao -- print ao overlap matrix
 1060 if(.not.equal(keywd,ksao,3)) go to 1070
        jprint(39) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(39),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(39).eq.ival) jprint(39) = ifull
          if(jprint(39).eq.ilew) jprint(39) = ifull
        end if
        go to 100
c  keyword: fao -- print ao fock matrix
 1070 if(.not.equal(keywd,kfao,3)) go to 1080
        jprint(40) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(40),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(40).eq.ival) jprint(40) = ifull
          if(jprint(40).eq.ilew) jprint(40) = ifull
        end if
        go to 100
c  keyword: nhonbo -- output nho to nbo transformation
 1080 if(.not.equal(keywd,knhonb,6)) go to 1090
        jprint(41) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(41),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(41).eq.ival) jprint(41) = ifull
        end if
        go to 100
c  keyword: boao -- print ao bond-order matrix
 1090 if(.not.equal(keywd,kboao,4)) go to 1100
        jprint(42) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(42),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(42).eq.ival) jprint(42) = ifull
          if(jprint(42).eq.ilew) jprint(42) = ifull
        end if
        go to 100
c  keyword: e2pert -- 2nd-order perturbative analysis of the nbo fock matrix
 1100 if(.not.equal(keywd,kepert,6)) go to 1110
        jprint(3) = 1
        call rfld(temp,error)
        if(error) go to 100
        e2thr = abs(temp)
        go to 100
c  keyword: plot -- write ao basis, density, and transforms for plotting
 1110 if(.not.equal(keywd,kplot,4)) go to 1120
        jprint(43) = 1
        go to 100
c  keyword: npa -- print the natural population analysis
 1120 if(.not.equal(keywd,knpa,3)) go to 1130
        jprint(4) = 1
        go to 100
c  keyword: nbosum -- print the nbo summary
 1130 if(.not.equal(keywd,knbosm,6)) go to 1140
        jprint(6) = 1
        go to 100
c  keyword: nbo -- print the nbo analysis
 1140 if(.not.equal(keywd,knbo,3)) go to 1150
        if(leng.ne.3) go to 1150
        jprint(5) = 1
        go to 100
c  keyword: dipole -- print nbo/nlmo dipole analysis:
 1150 if(.not.equal(keywd,kdipol,6)) go to 1160
        jprint(46) = 1
        call rfld(temp,error)
        if(error) go to 100
        dthr = abs(temp)
        go to 100
c  keyword: nbonlmo -- print nbo to nlmo transformation matrix
 1160 if(.not.equal(keywd,knbnlm,7)) go to 1170
        jprint(47) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(47),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(47).eq.ival) jprint(47) = ifull
        end if
        go to 100
c  keyword: spnlmo -- output the pnlmo overlap matrix
 1170 if(.not.equal(keywd,kspnlm,6)) go to 1180
        jprint(48) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(48),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(48).eq.ival) jprint(48) = ifull
        end if
        go to 100
c  keyword: aopnlmo -- output the ao-pnlmo transformation matrix
 1180 if(.not.equal(keywd,kaopnl,7)) go to 1190
        jprint(49) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(49),lfnpnl,read,error)
          if(error) nextwd = .false.
          if(jprint(49).eq.ival) jprint(49) = ifull
        end if
        go to 100
c  keyword: diao -- output the ao dipole integrals
 1190 if(.not.equal(keywd,kdiao,4)) go to 1200
        jprint(50) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(50),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(50).eq.ival) jprint(50) = ifull
          if(jprint(50).eq.ilew) jprint(50) = ifull
        end if
        go to 100
c  keyword: dinao -- output the nao dipole integrals
 1200 if(.not.equal(keywd,kdinao,5)) go to 1210
        jprint(51) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(51),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(51).eq.ival) jprint(51) = ifull
          if(jprint(51).eq.ilew) jprint(51) = ifull
        end if
        go to 100
c  keyword: dinho -- output the nho dipole integrals
 1210 if(.not.equal(keywd,kdinho,5)) go to 1220
        jprint(52) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(52),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(52).eq.ival) jprint(52) = ifull
          if(jprint(52).eq.ilew) jprint(52) = ifull
        end if
        go to 100
c  keyword: dinbo -- output the nbo dipole integrals
 1220 if(.not.equal(keywd,kdinbo,5)) go to 1230
        jprint(53) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(53),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(53).eq.ival) jprint(53) = ifull
        end if
        go to 100
c  keyword: dinlmo -- output the nlmo dipole integrals
 1230 if(.not.equal(keywd,kdinlm,6)) go to 1240
        jprint(54) = ifull
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          read = .false.
          call keypar(keywd,leng,jprint(54),lfndef,read,error)
          if(error) nextwd = .false.
          if(jprint(54).eq.ival) jprint(54) = ifull
        end if
        go to 100
c  keyword: nbodaf -- choose an alternate daf
 1240 if(.not.equal(keywd,knbdaf,6)) go to 1250
        lfndaf = abs(lfndaf)
        call ifld(itemp,error)
        if(.not.error) lfndaf = abs(itemp)
        go to 100
c  keyword: archive -- write the archive file
 1250 if(.not.equal(keywd,karchv,7)) go to 1260
        jprint(7) = lfnarc
        call ifld(itemp,error)
        if(.not.error) jprint(7) = abs(itemp)
        go to 100
c  keyword: file -- select alternate filename
 1260 if(.not.equal(keywd,kfile,4)) go to 1270
        leng = 80
        call hfld(inttmp,leng,end)
        if(.not.end) lennm = leng
        go to 100
c  keyword: apolar -- enforce apolar bonds:
 1270 if(.not.equal(keywd,kpolar,6)) go to 1290
        iwapol = 1
        go to 100
c  keyword: nrtopt -- optimize nrt weights:
 1290 if(.not.equal(keywd,knrtop,6)) go to 1300
        if(jprint(14).eq.0) jprint(14) = 1
        if(jprint(32).eq.0) jprint(32) = 1
        jprint(55) = ib
        leng = keylen
        call hfld(keywd,leng,end)
        if(.not.end) then
          if(equal(keywd,kbfgs,4)) then
            jprint(55) = ib
          else if(equal(keywd,kpowel,6)) then
            jprint(55) = ip
          else if(equal(keywd,ksap,3)) then
            jprint(55) = 1234567
            if(leng.gt.3) call convin(keywd(4),leng-3,jprint(55),error)
            if(error) call halt('nrtopt')
            jprint(55) = abs(jprint(55))
          else if(equal(keywd,ksap,2)) then
            jprint(55) = -1234567
            if(leng.gt.2) call convin(keywd(3),leng-2,jprint(55),error)
            if(error) call halt('nrtopt')
            jprint(55) = -abs(jprint(55))
          else
            nextwd = .false.
          end if
        end if
	go to 100
c  keyword: nrtref -- number of reference structures in nrt analysis
 1300 if(.not.equal(keywd,knrtrf,6)) go to 1310
        if(jprint(14).eq.0) jprint(14) = 1
        if(jprint(32).eq.0) jprint(32) = 1
        call ifld(itemp,error)
        if(error) go to 100
        jprint(56) = max(1,abs(itemp))
        go to 100
c  keyword: chsthr -- set the occupancy threshold in choose
 1310 if(.not.equal(keywd,kchsth,6)) go to 1320
        chsthr = abs(chsthr)
        call rfld(temp,error)
        if(error) go to 100
        chsthr = abs(temp)
        go to 100
c  keyword: nrtdtl -- detail nrt analysis
 1320 if(.not.equal(keywd,knrtdt,6)) go to 1340
        if(jprint(14).eq.0) jprint(14) = 1
        if(jprint(32).eq.0) jprint(32) = 1
        jprint(57) = 1
	call ifld(itemp,error)
        if(error) go to 100
        jprint(57) = max(1,abs(itemp))
        go to 100
c  keyword: nrtthr -- set threshold for delocalization list
 1340 if(.not.equal(keywd,knrtth,6)) go to 1360
        if(jprint(14).eq.0) jprint(14) = 1
        if(jprint(32).eq.0) jprint(32) = 1
        dlthr = abs(dlthr)
	call rfld(temp,error)
	if(error) go to 100
	dlthr = abs(temp)
        go to 100
c  keyword: nrt -- perform natural resonance theory analysis:
c  (note that we should check this keyword after we check the other
c  nrt keywords, like nrtopt.  otherwise, keyword conflicts can
c  occur.)
 1360 if(.not.equal(keywd,knrt,3)) go to 1370
        jprint(14) = 1
        jprint(32) = 1
        call ifld(itemp,error)
        if(.not.error) jprint(32) = itemp
        go to 100
 1370 go to 4800
c------------------------------------------------------------------------------
 4500 continue
c
c  if option file was selected, extract the filename from hollerith
c  array inttmp:
c
      if(lennm.ne.0) then
        idiv = ib - ia
        do 4510 i = 1,lennm
          filenm(i:i) = char(mod((inttmp(i)-ia)/idiv,256) + 65)
 4510   continue
        do 4520 i = lennm+1,80
          filenm(i:i) = char(32)
 4520   continue
      end if
c------------------------------------------------------------------------------
c
c  if the print level is set to zero and no other options were entered,
c  completely shut off program printing:
c
      if(numopt.eq.1.and.iprint.eq.0) iprint = -1
c
c  check for job options that are currently incompatable:
c
      if((iwdm.eq.0).and.(iwmulp.ne.0)) go to 4900
c
c  check for job options that are strictly incompatible:
c
      if(ortho) then
        iwtnao     = 0
        jprint(9)  = 0
        jprint(11) = 0
        jprint(18) = 0
        jprint(19) = 0
        jprint(20) = 0
        jprint(21) = 0
        jprint(25) = 0
        jprint(30) = 0
        jprint(31) = 0
        jprint(33) = 0
        jprint(35) = 0
        jprint(39) = 0
        jprint(44) = 0
        jprint(48) = 0
        jprint(49) = 0
        jprint(51) = 0
      end if
c------------------------------------------------------------------------------
c
c  start printing nbo output:
c
      if(iprint.ge.0) then
        write(lfnpr,6000)
        if(numopt.gt.0) write(lfnpr,6010)
c------------------------------------------------------------------------------
 6000 format(/1x,79('*')/,13x,
     * 'N A T U R A L   A T O M I C   O R B I T A L   A N D'/,
     * 10x,'N A T U R A L   B O N D   O R B I T A L   ',
     *    'A N A L Y S I S',/1x,79('*'))
 6010 format(1x)
c------------------------------------------------------------------------------
c
c  job control keywords:
c
        if(jprint(4).ne.0) write(lfnpr,6020)
        if(jprint(5).ne.0) write(lfnpr,6030)
        if(jprint(6).ne.0) write(lfnpr,6040)
        if(jprint(14).ne.0) write(lfnpr,6050)
        if(jprint(10).ne.0) write(lfnpr,6060)
        if(iw3c.ne.0) write(lfnpr,6070)
        if(jprint(1).ne.0) write(lfnpr,6080)
        if(jprint(8).ne.0) write(lfnpr,6090)
        if(jprint(32).ne.0) write(lfnpr,6100)
        if(jprint(55).eq.ib) then
          write(lfnpr,6110)
        else if(jprint(55).eq.ip) then
          write(lfnpr,6111)
        else if(jprint(55).lt.0) then
          write(lfnpr,6112)
        else if(jprint(55).gt.0) then
          write(lfnpr,6113)
        end if
        if(jprint(56).ne.0) write(lfnpr,6120) jprint(56)
        if(dlthr.ge.zero) write(lfnpr,6160) dlthr
        if(jprint(57).ne.0) write(lfnpr,6170) jprint(57)
c------------------------------------------------------------------------------
 6020 format(1x,'      /npa    / : print natural population analysis')
 6030 format(1x,'      /nbo    / : print natural bond orbital analysis')
 6040 format(1x,'      /nbosum / : print summary of the nbo analysis')
 6050 format(1x,'      /reson  / : allow strongly delocalized nbo ',
     *  'set')
 6060 format(1x,'      /nobond / : no two-center nbo search')
 6070 format(1x,'      /3cbond / : search for 3-center bonds')
 6080 format(1x,'      /skipbo / : skip nbo transformation step')
 6090 format(1x,'      /nlmo   / : form natural localized molecular',
     *  ' orbitals')
 6100 format(1x,'      /nrt    / : perform natural resonance theory ',
     * 'analysis')
 6110 format(1x,'      /nrtopt / : optimize resonance weights with ',
     * 'bfgs method')
 6111 format(1x,'      /nrtopt / : optimize resonance weights with ',
     * 'powell method')
 6112 format(1x,'      /nrtopt / : optimize resonance weights with ',
     * 'anneal method')
 6113 format(1x,'      /nrtopt / : optimize resonance weights with ',
     * 'anneal method + penalty')
 6120 format(1x,'      /nrtref / : number of reference structures set',
     * ' to',i3)
 6160 format(1x,'      /nrtthr / : set to ',f5.2)
 6170 format(1x,'      /nrtdtl / : set to ',i2)
c------------------------------------------------------------------------------
c
c  job threshold keywords:
c
        if(jprint(36).ne.0) write(lfnpr,6500)
        if(athr.ge.zero.or.pthr.ge.zero.or.ethr.ge.zero)
     +             write(lfnpr,6510) abs(athr),abs(pthr),abs(ethr)
        if(jprint(3).ne.0) write(lfnpr,6520)
        if(e2thr.gt.zero) write(lfnpr,6530) e2thr
        if(jprint(46).ne.0) write(lfnpr,6540)
        if(dthr.ge.zero) write(lfnpr,6550) abs(dthr)
        if(thrset.gt.zero) write(lfnpr,6560) thrset
        if(prjset.gt.zero) write(lfnpr,6570) prjset
        if(chsthr.gt.zero) write(lfnpr,6580) chsthr
c------------------------------------------------------------------------------
 6500 format(1x,'      /bend   / : print nho directionality table')
 6510 format(1x,'                  print thresholds set to (',f4.1,
     *   ',',f5.1,',',f5.2,')')
 6520 format(1x,'      /e2pert / : analyze nbo fock matrix')
 6530 format(1x,'                  print threshold set to ',f5.2)
 6540 format(1x,'      /dipole / : print nbo/nlmo dipole moment ',
     *   'analysis')
 6550 format(1x,'                  print threshold set to ',f5.2)
 6560 format(1x,'      /thresh / : set to ',f5.2)
 6570 format(1x,'      /prjthr / : set to ',f5.2)
 6580 format(1x,'      /chsthr / : set to ',f5.2)
c------------------------------------------------------------------------------
c
c  matrix output keywords:
c
        if(jprint(44).eq.ifull) then
          write(lfnpr,7000)
        else if(ioinqr(jprint(44)).eq.iprnt) then
          write(lfnpr,7002) jprint(44)
        else if(ioinqr(jprint(44)).eq.iwrit) then
          write(lfnpr,7004) abs(jprint(44))
        end if
        if(iwtnao.eq.ifull) then
          write(lfnpr,7010)
        else if(ioinqr(iwtnao).eq.iprnt) then
          write(lfnpr,7012) iwtnao
        else if(ioinqr(iwtnao).eq.iwrit) then
          write(lfnpr,7014) abs(iwtnao)
        else if(ioinqr(iwtnao).eq.iread) then
          write(lfnpr,7016) abs(iwtnao/1000)
        end if
        if(jprint(30).eq.ifull) then
          write(lfnpr,7020)
        else if(ioinqr(jprint(30)).eq.iprnt) then
          write(lfnpr,7022) jprint(30)
        else if(ioinqr(jprint(30)).eq.iwrit) then
          write(lfnpr,7024) abs(jprint(30))
        end if
        if(jprint(28).eq.ifull) then
          write(lfnpr,7030)
        else if(ioinqr(jprint(28)).eq.iprnt) then
          write(lfnpr,7032) jprint(28)
        else if(ioinqr(jprint(28)).eq.iwrit) then
          write(lfnpr,7034) abs(jprint(28))
        end if
        if(jprint(25).eq.ifull) then
          write(lfnpr,7040)
        else if(jprint(25).eq.ilew) then
          write(lfnpr,7042)
        else if(ioinqr(jprint(25)).eq.iprnt) then
          write(lfnpr,7044) jprint(25)
        else if(ioinqr(jprint(25)).eq.iwrit) then
          write(lfnpr,7046) abs(jprint(25))
        end if
        if(iwtnbo.eq.ifull) then
          write(lfnpr,7050)
        else if(iwtnbo.eq.ilew) then
          write(lfnpr,7052)
        else if(ioinqr(iwtnbo).eq.iprnt) then
          write(lfnpr,7054) iwtnbo
        else if(ioinqr(iwtnbo).eq.iwrit) then
          write(lfnpr,7056) abs(iwtnbo)
        end if
        if(jprint(49).eq.ifull) then
          write(lfnpr,7060)
        else if(jprint(49).eq.ilew) then
          write(lfnpr,7062)
        else if(ioinqr(jprint(49)).eq.iprnt) then
          write(lfnpr,7064) jprint(49)
        else if(ioinqr(jprint(49)).eq.iwrit) then
          write(lfnpr,7066) abs(jprint(49))
        end if
        if(jprint(23).eq.ifull) then
          write(lfnpr,7070)
        else if(jprint(23).eq.ilew) then
          write(lfnpr,7072)
        else if(ioinqr(jprint(23)).eq.iprnt) then
          write(lfnpr,7074) jprint(23)
        else if(ioinqr(jprint(23)).eq.iwrit) then
          write(lfnpr,7076) abs(jprint(23))
        end if
        if(jprint(26).eq.ifull) then
          write(lfnpr,7080)
        else if(jprint(26).eq.ival) then
          write(lfnpr,7082)
        else if(jprint(26).eq.ilew) then
          write(lfnpr,7084)
        else if(ioinqr(jprint(26)).eq.iprnt) then
          write(lfnpr,7086) jprint(26)
        else if(ioinqr(jprint(26)).eq.iwrit) then
          write(lfnpr,7088) abs(jprint(26))
        end if
        if(iwpnao.eq.ifull) then
          write(lfnpr,7090)
        else if(ioinqr(iwpnao).eq.iprnt) then
          write(lfnpr,7092) iwpnao
        else if(ioinqr(iwpnao).eq.iwrit) then
          write(lfnpr,7094) abs(iwpnao)
        else if(ioinqr(iwpnao).eq.iread) then
          write(lfnpr,7096) abs(iwpnao/1000)
        end if
        if(jprint(33).eq.ifull) then
          write(lfnpr,7100)
        else if(ioinqr(jprint(33)).eq.iprnt) then
          write(lfnpr,7102) jprint(33)
        else if(ioinqr(jprint(33)).eq.iwrit) then
          write(lfnpr,7104) abs(jprint(33))
        end if
        if(iwtnab.eq.ifull) then
          write(lfnpr,7110)
        else if(iwtnab.eq.ilew) then
          write(lfnpr,7112)
        else if(ioinqr(iwtnab).eq.iprnt) then
          write(lfnpr,7114) iwtnab
        else if(ioinqr(iwtnab).eq.iwrit) then
          write(lfnpr,7116) abs(iwtnab)
        else if(ioinqr(iwtnab).eq.iread) then
          write(lfnpr,7118) abs(iwtnab/1000)
        end if
        if(jprint(18).eq.ifull) then
          write(lfnpr,7120)
        else if(jprint(18).eq.ilew) then
          write(lfnpr,7122)
        else if(ioinqr(jprint(18)).eq.iprnt) then
          write(lfnpr,7124) jprint(18)
        else if(ioinqr(jprint(18)).eq.iwrit) then
          write(lfnpr,7126) abs(jprint(18))
        end if
        if(jprint(9).eq.ifull) then
          write(lfnpr,7130)
        else if(jprint(9).eq.ival) then
          write(lfnpr,7132)
        else if(jprint(9).eq.ilew) then
          write(lfnpr,7134)
        else if(ioinqr(jprint(9)).eq.iprnt) then
          write(lfnpr,7136) jprint(9)
        else if(ioinqr(jprint(9)).eq.iwrit) then
          write(lfnpr,7138) abs(jprint(9))
        end if
        if(jprint(41).eq.ifull) then
          write(lfnpr,7140)
        else if(jprint(41).eq.ilew) then
          write(lfnpr,7142)
        else if(ioinqr(jprint(41)).eq.iprnt) then
          write(lfnpr,7144) jprint(41)
        else if(ioinqr(jprint(41)).eq.iwrit) then
          write(lfnpr,7146) abs(jprint(41))
        end if
        if(jprint(24).eq.ifull) then
          write(lfnpr,7150)
        else if(jprint(24).eq.ilew) then
          write(lfnpr,7152)
        else if(ioinqr(jprint(24)).eq.iprnt) then
          write(lfnpr,7154) jprint(24)
        else if(ioinqr(jprint(24)).eq.iwrit) then
          write(lfnpr,7156) abs(jprint(24))
        end if
        if(jprint(38).eq.ifull) then
          write(lfnpr,7160)
        else if(jprint(38).eq.ival) then
          write(lfnpr,7162)
        else if(jprint(38).eq.ilew) then
          write(lfnpr,7164)
        else if(ioinqr(jprint(38)).eq.iprnt) then
          write(lfnpr,7166) jprint(38)
        else if(ioinqr(jprint(38)).eq.iwrit) then
          write(lfnpr,7168) abs(jprint(38))
        end if
        if(jprint(47).eq.ifull) then
          write(lfnpr,7170)
        else if(jprint(47).eq.ilew) then
          write(lfnpr,7172)
        else if(ioinqr(jprint(47)).eq.iprnt) then
          write(lfnpr,7174) jprint(47)
        else if(ioinqr(jprint(47)).eq.iwrit) then
          write(lfnpr,7176) abs(jprint(47))
        end if
        if(jprint(45).eq.ifull) then
          write(lfnpr,7180)
        else if(jprint(45).eq.ival) then
          write(lfnpr,7182)
        else if(jprint(45).eq.ilew) then
          write(lfnpr,7184)
        else if(ioinqr(jprint(45)).eq.iprnt) then
          write(lfnpr,7186) jprint(45)
        else if(ioinqr(jprint(45)).eq.iwrit) then
          write(lfnpr,7188) abs(jprint(45))
        end if
        if(jprint(13).eq.ifull) then
          write(lfnpr,7190)
        else if(jprint(13).eq.ival) then
          write(lfnpr,7192)
        else if(jprint(13).eq.ilew) then
          write(lfnpr,7194)
        else if(ioinqr(jprint(13)).eq.iprnt) then
          write(lfnpr,7196) jprint(13)
        else if(ioinqr(jprint(13)).eq.iwrit) then
          write(lfnpr,7198) abs(jprint(13))
        end if
        if(jprint(42).eq.ifull) then
          write(lfnpr,7200)
        else if(ioinqr(jprint(42)).eq.iprnt) then
          write(lfnpr,7202) jprint(42)
        else if(ioinqr(jprint(42)).eq.iwrit) then
          write(lfnpr,7204) abs(jprint(42))
        end if
        if(jprint(27).eq.ifull) then
          write(lfnpr,7210)
        else if(ioinqr(jprint(27)).eq.iprnt) then
          write(lfnpr,7212) jprint(27)
        else if(ioinqr(jprint(27)).eq.iwrit) then
          write(lfnpr,7214) abs(jprint(27))
        end if
        if(jprint(35).eq.ifull) then
          write(lfnpr,7220)
        else if(ioinqr(jprint(35)).eq.iprnt) then
          write(lfnpr,7222) jprint(35)
        else if(ioinqr(jprint(35)).eq.iwrit) then
          write(lfnpr,7224) abs(jprint(35))
        end if
        if(jprint(34).eq.ifull) then
          write(lfnpr,7230)
        else if(ioinqr(jprint(34)).eq.iprnt) then
          write(lfnpr,7232) jprint(34)
        else if(ioinqr(jprint(34)).eq.iwrit) then
          write(lfnpr,7234) abs(jprint(34))
        end if
        if(jprint(16).eq.ifull) then
          write(lfnpr,7240)
        else if(jprint(16).eq.ilew) then
          write(lfnpr,7242)
        else if(ioinqr(jprint(16)).eq.iprnt) then
          write(lfnpr,7244) jprint(16)
        else if(ioinqr(jprint(16)).eq.iwrit) then
          write(lfnpr,7246) abs(jprint(16))
        end if
        if(jprint(17).eq.ifull) then
          write(lfnpr,7250)
        else if(jprint(17).eq.ilew) then
          write(lfnpr,7252)
        else if(ioinqr(jprint(17)).eq.iprnt) then
          write(lfnpr,7254) jprint(17)
        else if(ioinqr(jprint(17)).eq.iwrit) then
          write(lfnpr,7256) abs(jprint(17))
        end if
        if(jprint(40).eq.ifull) then
          write(lfnpr,7260)
        else if(ioinqr(jprint(40)).eq.iprnt) then
          write(lfnpr,7262) jprint(40)
        else if(ioinqr(jprint(40)).eq.iwrit) then
          write(lfnpr,7264) abs(jprint(40))
        end if
        if(jprint(31).eq.ifull) then
          write(lfnpr,7270)
        else if(ioinqr(jprint(31)).eq.iprnt) then
          write(lfnpr,7272) jprint(31)
        else if(ioinqr(jprint(31)).eq.iwrit) then
          write(lfnpr,7274) abs(jprint(31))
        end if
        if(jprint(29).eq.ifull) then
          write(lfnpr,7280)
        else if(ioinqr(jprint(29)).eq.iprnt) then
          write(lfnpr,7282) jprint(29)
        else if(ioinqr(jprint(29)).eq.iwrit) then
          write(lfnpr,7284) abs(jprint(29))
        end if
        if(jprint(37).eq.ifull) then
          write(lfnpr,7290)
        else if(jprint(37).eq.ilew) then
          write(lfnpr,7292)
        else if(ioinqr(jprint(37)).eq.iprnt) then
          write(lfnpr,7294) jprint(37)
        else if(ioinqr(jprint(37)).eq.iwrit) then
          write(lfnpr,7296) abs(jprint(37))
        end if
        if(jprint(15).eq.ifull) then
          write(lfnpr,7300)
        else if(jprint(15).eq.ilew) then
          write(lfnpr,7302)
        else if(ioinqr(jprint(15)).eq.iprnt) then
          write(lfnpr,7304) jprint(15)
        else if(ioinqr(jprint(15)).eq.iwrit) then
          write(lfnpr,7306) abs(jprint(15))
        end if
        if(jprint(50).eq.ifull) then
          write(lfnpr,7310)
        else if(ioinqr(jprint(50)).eq.iprnt) then
          write(lfnpr,7312) jprint(50)
        else if(ioinqr(jprint(50)).eq.iwrit) then
          write(lfnpr,7314) abs(jprint(50))
        end if
        if(jprint(51).eq.ifull) then
          write(lfnpr,7320)
        else if(ioinqr(jprint(51)).eq.iprnt) then
          write(lfnpr,7322) jprint(51)
        else if(ioinqr(jprint(51)).eq.iwrit) then
          write(lfnpr,7324) abs(jprint(51))
        end if
        if(jprint(52).eq.ifull) then
          write(lfnpr,7330)
        else if(ioinqr(jprint(52)).eq.iprnt) then
          write(lfnpr,7332) jprint(52)
        else if(ioinqr(jprint(52)).eq.iwrit) then
          write(lfnpr,7334) abs(jprint(52))
        end if
        if(jprint(53).eq.ifull) then
          write(lfnpr,7340)
        else if(jprint(53).eq.ilew) then
          write(lfnpr,7342)
        else if(ioinqr(jprint(53)).eq.iprnt) then
          write(lfnpr,7344) jprint(53)
        else if(ioinqr(jprint(53)).eq.iwrit) then
          write(lfnpr,7346) abs(jprint(53))
        end if
        if(jprint(54).eq.ifull) then
          write(lfnpr,7350)
        else if(jprint(54).eq.ilew) then
          write(lfnpr,7352)
        else if(ioinqr(jprint(54)).eq.iprnt) then
          write(lfnpr,7354) jprint(54)
        else if(ioinqr(jprint(54)).eq.iwrit) then
          write(lfnpr,7356) abs(jprint(54))
        end if
        if(jprint(39).eq.ifull) then
          write(lfnpr,7360)
        else if(ioinqr(jprint(39)).eq.iprnt) then
          write(lfnpr,7362) jprint(39)
        else if(ioinqr(jprint(39)).eq.iwrit) then
          write(lfnpr,7364) abs(jprint(39))
        end if
        if(jprint(19).eq.ifull) then
          write(lfnpr,7370)
        else if(ioinqr(jprint(19)).eq.iprnt) then
          write(lfnpr,7372) jprint(19)
        else if(ioinqr(jprint(19)).eq.iwrit) then
          write(lfnpr,7374) abs(jprint(19))
        end if
        if(jprint(20).eq.ifull) then
          write(lfnpr,7380)
        else if(ioinqr(jprint(20)).eq.iprnt) then
          write(lfnpr,7382) jprint(20)
        else if(ioinqr(jprint(20)).eq.iwrit) then
          write(lfnpr,7384) abs(jprint(20))
        end if
        if(jprint(21).eq.ifull) then
          write(lfnpr,7390)
        else if(jprint(21).eq.ilew) then
          write(lfnpr,7392)
        else if(ioinqr(jprint(21)).eq.iprnt) then
          write(lfnpr,7394) jprint(21)
        else if(ioinqr(jprint(21)).eq.iwrit) then
          write(lfnpr,7396) abs(jprint(21))
        end if
        if(jprint(48).eq.ifull) then
          write(lfnpr,7400)
        else if(jprint(48).eq.ilew) then
          write(lfnpr,7402)
        else if(ioinqr(jprint(48)).eq.iprnt) then
          write(lfnpr,7404) jprint(48)
        else if(ioinqr(jprint(48)).eq.iwrit) then
          write(lfnpr,7406) abs(jprint(48))
        end if
c------------------------------------------------------------------------------
 7000 format(1x,'      /aopnao / : print the ao to pnao transformation')
 7002 format(1x,'      /aopnao / : print ',i3,' columns of the ao to ',
     *   'pnao transformation')
 7004 format(1x,'      /aopnao / : write the ao to pnao transformation',
     *   ' to lfn',i3)
 7010 format(1x,'      /aonao  / : print the ao to nao transformation')
 7012 format(1x,'      /aonao  / : print ',i3,' columns of the ao ',
     *   'to nao transformation')
 7014 format(1x,'      /aonao  / : write the ao to nao transformation ',
     *   'to lfn',i3)
 7016 format(1x,'      /aonao  / : read ao to nao transformation from ',
     *          'lfn',i3)
 7020 format(1x,'      /aopnho / : print the ao to pnho ',
     *   'transformation')
 7022 format(1x,'      /aopnho / : print ',i3,' columns of the ao to ',
     *   'pnho transformation')
 7024 format(1x,'      /aopnho / : write the ao to pnho transformation',
     *   ' to lfn',i3)
 7030 format(1x,'      /aonho  / : print the ao to nho transformation')
 7032 format(1x,'      /aonho  / : print ',i3,' columns of the ao to ',
     *   'nho transformation')
 7034 format(1x,'      /aonho  / : write the ao to nho transformation ',
     *   'to lfn',i3)
 7040 format(1x,'      /aopnbo / : print the ao to pnbo ',
     *   'transformation')
 7042 format(1x,'      /aopnbo / : print the occupied pnbos in the ao ',
     *   'basis')
 7044 format(1x,'      /aopnbo / : print ',i3,' columns of the ao to ',
     *   'pnbo transformation')
 7046 format(1x,'      /aopnbo / : write the ao to pnbo transformation',
     *   ' to lfn',i3)
 7050 format(1x,'      /aonbo  / : print the ao to nbo transformation')
 7052 format(1x,'      /aonbo  / : print the occupied nbos in the ao ',
     *   'basis')
 7054 format(1x,'      /aonbo  / : print ',i3,' columns of the ao ',
     *   'to nbo transformation')
 7056 format(1x,'      /aonbo  / : write the ao to nbo transformation ',
     *   'to lfn',i3)
 7060 format(1x,'      /aopnlmo/ : print the ao to pnlmo ',
     *   'transformation')
 7062 format(1x,'      /aopnlmo/ : print the occupied pnlmos in the ao',
     *   ' basis')
 7064 format(1x,'      /aopnlmo/ : print ',i3,' columns of the ao to ',
     *   'pnlmo transformation')
 7066 format(1x,'      /aopnlmo/ : write the ao to pnlmo transformatio',
     *   'n to lfn',i3)
 7070 format(1x,'      /aonlmo / : print the ao to nlmo ',
     *   'transformation')
 7072 format(1x,'      /aonlmo / : print the occupied nlmos in the ao ',
     *   'basis')
 7074 format(1x,'      /aonlmo / : print ',i3,' columns of the ao to ',
     *   'nlmo transformation')
 7076 format(1x,'      /aonlmo / : write the ao to nlmo transformation',
     *   ' to lfn',i3)
 7080 format(1x,'      /aomo   / : print all mos in the ao basis')
 7082 format(1x,'      /aomo   / : print core and valence mos in ',
     *   'the ao basis')
 7084 format(1x,'      /aomo   / : print the occupied mos in the ao ',
     *   'basis')
 7086 format(1x,'      /aomo   / : print ',i3,' lowest energy mos ',
     *   'in the ao basis')
 7088 format(1x,'      /aomo   / : write the ao to mo transformation ',
     *   'to lfn',i3)
 7090 format(1x,'      /paopnao/ : print the pao to pnao ',
     *   'transformation')
 7092 format(1x,'      /paopnao/ : print ',i3,' columns of the pao ',
     *   'to pnao transformation')
 7094 format(1x,'      /paopnao/ : write the pao to pnao ',
     *   'transformation to lfn',i3)
 7096 format(1x,'      /paopnao/ : read pao to pnao transformation ',
     *          'from lfn',i3)
 7100 format(1x,'      /naonho / : print the nao to nho transformation')
 7102 format(1x,'      /naonho / : print ',i3,' columns of the nao ',
     *   'to nho transformation')
 7104 format(1x,'      /naonho / : write the nao to nho transformation',
     *   ' to lfn',i3)
 7110 format(1x,'      /naonbo / : print the nao to nbo transformation')
 7112 format(1x,'      /naonbo / : print the occupied nbos in the nao ',
     *   'basis')
 7114 format(1x,'      /naonbo / : print ',i3,' columns of the nao ',
     *   'to nbo transformation')
 7116 format(1x,'      /naonbo / : write the nao to nbo transformation',
     *   ' to lfn',i3)
 7118 format(1x,'      /naonbo / : read nao to nbo transformation from',
     *          ' lfn',i3)
 7120 format(1x,'      /naonlmo/ : print the nao to nlmo ',
     *   'transformation')
 7122 format(1x,'      /naonlmo/ : print the occupied nlmos in the nao',
     *   ' basis')
 7124 format(1x,'      /naonlmo/ : print ',i3,' columns of the nao ',
     *   'to nlmo transformation')
 7126 format(1x,'      /naonlmo/ : write the nao to nlmo ',
     *   'transformation to lfn',i3)
 7130 format(1x,'      /naomo  / : print all mos in the nao basis')
 7132 format(1x,'      /naomo  / : print core and valence mos in ',
     *   'the nao basis')
 7134 format(1x,'      /naomo  / : print the occupied mos in the nao ',
     *   'basis')
 7136 format(1x,'      /naomo  / : print ',i3,' lowest energy mos ',
     *   'in the nao basis')
 7138 format(1x,'      /naomo  / : write the nao to mo transformation ',
     *   'to lfn',i3)
 7140 format(1x,'      /nhonbo / : print the nho to nbo transformation')
 7142 format(1x,'      /nhonbo / : print the occupied nbos in the nho ',
     *   'basis')
 7144 format(1x,'      /nhonbo / : print ',i3,' columns of the nho ',
     *   'to nbo transformation')
 7146 format(1x,'      /nhonbo / : write the nho to nbo transformation',
     *   ' to lfn',i3)
 7150 format(1x,'      /nhonlmo/ : print the nho to nlmo ',
     *   'transformation')
 7152 format(1x,'      /nhonlmo/ : print the occupied nlmos in the nho',
     *   ' basis')
 7154 format(1x,'      /nhonlmo/ : print ',i3,' columns of the nho ',
     *   'to nlmo transformation')
 7156 format(1x,'      /nhonlmo/ : write the nho to nlmo ',
     *   'transformation to lfn',i3)
 7160 format(1x,'      /nhomo  / : print all mos in the nho basis')
 7162 format(1x,'      /nhomo  / : print core and valence mos in ',
     *   'the nho basis')
 7164 format(1x,'      /nhomo  / : print the occupied mos in the nho ',
     *   'basis')
 7166 format(1x,'      /nhomo  / : print ',i3,' lowest energy mos ',
     *   'in the nho basis')
 7168 format(1x,'      /nhomo  / : write the nho to mo transformation ',
     *   'to lfn',i3)
 7170 format(1x,'      /nbonlmo/ : print the nbo to nlmo ',
     *   'transformation')
 7172 format(1x,'      /nbonlmo/ : print the occupied nlmos in the nbo',
     *   ' basis')
 7174 format(1x,'      /nbonlmo/ : print ',i3,' columns of the nbo ',
     *   'to nlmo transformation')
 7176 format(1x,'      /nbonlmo/ : write the nbo to nlmo ',
     *   'transformation to lfn',i3)
 7180 format(1x,'      /nbomo  / : print all mos in the nbo basis')
 7182 format(1x,'      /nbomo  / : print core and valence mos in ',
     *   'the nbo basis')
 7184 format(1x,'      /nbomo  / : print the occupied mos in the nbo ',
     *   'basis')
 7186 format(1x,'      /nbomo  / : print ',i3,' lowest energy mos ',
     *   'in the nbo basis')
 7188 format(1x,'      /nbomo  / : write the nbo to mo transformation ',
     *   'to lfn',i3)
 7190 format(1x,'      /nlmomo / : print all mos in the nlmo basis')
 7192 format(1x,'      /nlmomo / : print core and valence mos in ',
     *   'the nlmo basis')
 7194 format(1x,'      /nlmomo / : print the occupied mos in the nlmo ',
     *   'basis')
 7196 format(1x,'      /nlmomo / : print ',i3,' lowest energy mos ',
     *   'in the nlmo basis')
 7198 format(1x,'      /nlmomo / : write the nlmo to mo transformation',
     *   ' to lfn',i3)
 7200 format(1x,'      /boao   / : print the ao bond-order matrix')
 7202 format(1x,'      /boao   / : print ',i3,' columns of the ao ',
     *   'bond-order matrix')
 7204 format(1x,'      /boao   / : write the ao bond-order matrix to ',
     *   'lfn',i3)
 7210 format(1x,'      /dmao   / : print the ao density matrix')
 7212 format(1x,'      /dmao   / : print ',i3,' columns of the ao ',
     *   'density matrix')
 7214 format(1x,'      /dmao   / : write the ao density matrix to ',
     *   'lfn',i3)
 7220 format(1x,'      /dmnao  / : print the nao density matrix')
 7222 format(1x,'      /dmnao  / : print ',i3,' columns of the nao ',
     *   'density matrix')
 7224 format(1x,'      /dmnao  / : write the nao density matrix to ',
     *   'lfn',i3)
 7230 format(1x,'      /dmnho  / : print the nho density matrix')
 7232 format(1x,'      /dmnho  / : print ',i3,' columns of the nho ',
     *   'density matrix')
 7234 format(1x,'      /dmnho  / : write the nho density matrix to ',
     *   'lfn',i3)
 7240 format(1x,'      /dmnbo  / : print the nbo density matrix')
 7242 format(1x,'      /dmnbo  / : print the density matrix elements ',
     *   'of the occupied nbos')
 7244 format(1x,'      /dmnbo  / : print ',i3,' columns of the nbo ',
     *   'density matrix')
 7246 format(1x,'      /dmnbo  / : write the nbo density matrix to ',
     *   'lfn',i3)
 7250 format(1x,'      /dmnlmo / : print the nlmo density matrix')
 7252 format(1x,'      /dmnlmo / : print the density matrix elements ',
     *   'of the occupied nlmos')
 7254 format(1x,'      /dmnlmo / : print ',i3,' columns of the nlmo ',
     *   'density matrix')
 7256 format(1x,'      /dmnlmo / : write the nlmo density matrix to ',
     *   'lfn',i3)
 7260 format(1x,'      /fao    / : print the ao fock matrix')
 7262 format(1x,'      /fao    / : print ',i3,' columns of the ao ',
     *   'fock matrix')
 7264 format(1x,'      /fao    / : write the ao fock matrix to ',
     *   'lfn',i3)
 7270 format(1x,'      /fnao   / : print the nao fock matrix')
 7272 format(1x,'      /fnao   / : print ',i3,' columns of the nao ',
     *   'fock matrix')
 7274 format(1x,'      /fnao   / : write the nao fock matrix to ',
     *   'lfn',i3)
 7280 format(1x,'      /fnho   / : print the nho fock matrix')
 7282 format(1x,'      /fnho   / : print ',i3,' columns of the nho ',
     *   'fock matrix')
 7284 format(1x,'      /fnho   / : write the nho fock matrix to ',
     *   'lfn',i3)
 7290 format(1x,'      /fnbo   / : print the nbo fock matrix')
 7292 format(1x,'      /fnbo   / : print the fock matrix elements of ',
     *   'the occupied nbos')
 7294 format(1x,'      /fnbo   / : print ',i3,' columns of the nbo ',
     *   'fock matrix')
 7296 format(1x,'      /fnbo   / : write the nbo fock matrix to ',
     *   'lfn',i3)
 7300 format(1x,'      /fnlmo  / : print the nlmo fock matrix')
 7302 format(1x,'      /fnlmo  / : print the fock matrix elements of ',
     *   'the occupied nlmos')
 7304 format(1x,'      /fnlmo  / : print ',i3,' columns of the nlmo ',
     *   'fock matrix')
 7306 format(1x,'      /fnlmo  / : write the nlmo fock matrix to ',
     *   'lfn',i3)
 7310 format(1x,'      /diao   / : print the ao dipole integrals')
 7312 format(1x,'      /diao   / : print ',i3,' columns of the ao ',
     *   'dipole matrices')
 7314 format(1x,'      /diao   / : write the ao dipole integrals',
     *   ' to lfn',i3)
 7320 format(1x,'      /dinao  / : print the nao dipole integrals')
 7322 format(1x,'      /dinao  / : print ',i3,' columns of the nao ',
     *   'dipole matrices')
 7324 format(1x,'      /dinao  / : write the nao dipole integrals',
     *   ' to lfn',i3)
 7330 format(1x,'      /dinho  / : print the nho dipole integrals')
 7332 format(1x,'      /dinho  / : print ',i3,' columns of the nho ',
     *   'dipole matrices')
 7334 format(1x,'      /dinho  / : write the nho dipole integrals',
     *   ' to lfn',i3)
 7340 format(1x,'      /dinbo  / : print the nbo dipole integrals')
 7342 format(1x,'      /dinbo  / : print the dipole integrals of ',
     *   'occupied nbos')
 7344 format(1x,'      /dinbo  / : print ',i3,' columns of the nbo ',
     *   'dipole matrices')
 7346 format(1x,'      /dinbo  / : write the nbo dipole integrals',
     *   ' to lfn',i3)
 7350 format(1x,'      /dinlmo / : print the nlmo dipole integrals')
 7352 format(1x,'      /dinlmo / : print the dipole integrals of ',
     *   'occupied nlmos')
 7354 format(1x,'      /dinlmo / : print ',i3,' columns of the nlmo ',
     *   'dipole matrices')
 7356 format(1x,'      /dinlmo / : write the nlmo dipole integrals',
     *   ' to lfn',i3)
 7360 format(1x,'      /sao    / : print the ao overlap matrix')
 7362 format(1x,'      /sao    / : print ',i3,' columns of the ao ',
     *   'overlap matrix')
 7364 format(1x,'      /sao    / : write the ao overlap matrix to ',
     *   'lfn',i3)
 7370 format(1x,'      /spnao  / : print the pnao overlap matrix')
 7372 format(1x,'      /spnao  / : print ',i3,' columns of the pnao ',
     *   'overlap matrix')
 7374 format(1x,'      /spnao  / : write the pnao overlap matrix to ',
     *   'lfn',i3)
 7380 format(1x,'      /spnho  / : print the pnho overlap matrix')
 7382 format(1x,'      /spnho  / : print ',i3,' columns of the pnho ',
     *   'overlap matrix')
 7384 format(1x,'      /spnho  / : write the pnho overlap matrix to ',
     *   'lfn',i3)
 7390 format(1x,'      /spnbo  / : print the pnbo overlap matrix')
 7392 format(1x,'      /spnbo  / : print the overlap matrix elements ',
     *   'of the occupied pnbos')
 7394 format(1x,'      /spnbo  / : print ',i3,' columns of the pnbo ',
     *   'overlap matrix')
 7396 format(1x,'      /spnbo  / : write the pnbo overlap matrix to ',
     *   'lfn',i3)
 7400 format(1x,'      /spnlmo / : print the pnlmo overlap matrix')
 7402 format(1x,'      /spnlmo / : print the overlap matrix elements ',
     *   'of the occupied pnlmos')
 7404 format(1x,'      /spnlmo / : print ',i3,' columns of the pnlmo ',
     *   'overlap matrix')
 7406 format(1x,'      /spnlmo / : write the pnlmo overlap matrix to ',
     *   'lfn',i3)
c------------------------------------------------------------------------------
c
c  other output control keywords:
c
        if(lfnpr.ne.6) write(lfnpr,8000) lfnpr
        if(jprint(43).ne.0) write(lfnpr,8010)
        if(iwdetl.ne.0) write(lfnpr,8020)
        if(jprint(7).ne.0) write(lfnpr,8030) jprint(7)
        if(jprint(12).ne.0) write(lfnpr,8040)
        if(lfndaf.ge.0) write(lfnpr,8050) lfndaf
        if(jprint(22).ne.0) write(lfnpr,8060) jprint(22)
        if(iwmulp.eq.1) write(lfnpr,8070)
        if(iwmulp.eq.2) write(lfnpr,8080)
        if(iwapol.ne.0) write(lfnpr,8090)
        if(jprint(11).ne.0) write(lfnpr,8100)
        if(lennm.ne.0) write(lfnpr,8110) filenm(1:52)
c
        if(iprint.lt.10) then
          write(lfnpr,8500) iprint
        else
          iprint = iprint - 10
        end if
c------------------------------------------------------------------------------
 8000 format(1x,'      /lfnpr  / : set to',i3)
 8010 format(1x,'      /plot   / : write information for the orbital',
     *   ' plotter')
 8020 format(1x,'      /detail / : write out details of nbo search')
 8030 format(1x,'      /archive/ : write the archive file to lfn',i3)
 8040 format(1x,'      /bndidx / : print bond indices based on ',
     *  'the nao density matrix')
 8050 format(1x,'      /nbodaf / : nbo direct access file written on',
     *   ' lfn',i3)
 8060 format(1x,'      /aoinfo / : write ao information to lfn',i3)
 8070 format(1x,'      /mulat  / : print mulliken populations',
     *                ' by atom')
 8080 format(1x,'      /mulorb / : print mulliken populations',
     *                ' by orbital and atom')
 8090 format(1x,'      /apolar / : enforce apolar nbos')
 8100 format(1x,'      /rpnao  / : revise tpnao with tryd and tred')
 8110 format(1x,'      /file   / : set to ',a52)
 8500 format(1x,'      /print  / : print level set to',i3)
c------------------------------------------------------------------------------
      end if
c
c  set print level options:
c
      if(iprint.gt.0) then
        jprint(4)  =  1
        jprint(5)  =  1
      end if
c
      if(iprint.gt.1) then
        jprint(3)  =  1
        jprint(6)  =  1
        jprint(36) =  1
      end if
c
      if(iprint.gt.2) then
        jprint(8)  =  1
        jprint(12) =  1
        jprint(46) =  1
      end if
c
      if(iprint.gt.3) then
        if(jprint(7).eq.0)  jprint(7)  = lfnarc
        if(jprint(9).eq.0)  jprint(9)  = ifull
        if(jprint(13).eq.0) jprint(13) = ifull
                            jprint(14) = 1
        if(jprint(15).eq.0) jprint(15) = ifull
        if(jprint(16).eq.0) jprint(16) = ifull
        if(jprint(17).eq.0) jprint(17) = ifull
        if(jprint(18).eq.0) jprint(18) = ifull
        if(jprint(19).eq.0) jprint(19) = ifull
        if(jprint(20).eq.0) jprint(20) = ifull
        if(jprint(21).eq.0) jprint(21) = ifull
        if(jprint(24).eq.0) jprint(24) = ifull
        if(jprint(29).eq.0) jprint(29) = ifull
        if(jprint(31).eq.0) jprint(31) = ifull
        if(jprint(32).eq.0) jprint(32) = 1
        if(jprint(33).eq.0) jprint(33) = ifull
        if(jprint(34).eq.0) jprint(34) = ifull
        if(jprint(35).eq.0) jprint(35) = ifull
        if(jprint(37).eq.0) jprint(37) = ifull
        if(jprint(38).eq.0) jprint(38) = ifull
        if(jprint(39).eq.0) jprint(39) = ifull
        if(jprint(40).eq.0) jprint(40) = ifull
        if(jprint(41).eq.0) jprint(41) = ifull
        if(jprint(42).eq.0) jprint(42) = ifull
                            jprint(43) = 1
        if(jprint(45).eq.0) jprint(45) = ifull
        if(jprint(47).eq.0) jprint(47) = ifull
        if(jprint(48).eq.0) jprint(48) = ifull
        if(jprint(50).eq.0) jprint(50) = ifull
        if(jprint(51).eq.0) jprint(51) = ifull
        if(jprint(52).eq.0) jprint(52) = ifull
        if(jprint(53).eq.0) jprint(53) = ifull
        if(jprint(54).eq.0) jprint(54) = ifull
        if(jprint(55).eq.0) jprint(55) = 1
        if(iwtnab.eq.0)     iwtnab     = ifull
                            iwdetl     = 1
        if(iwdm.ne.0)       iwmulp     = 2
      end if
c
c  turn on the nlmo analysis if required:
c
      if(jprint(13).ne.0) jprint(8) = 1
      if(jprint(15).ne.0) jprint(8) = 1
      if(jprint(17).ne.0) jprint(8) = 1
      if(jprint(18).ne.0) jprint(8) = 1
      if(jprint(23).ne.0) jprint(8) = 1
      if(jprint(46).ne.0) jprint(8) = 1
      if(jprint(47).ne.0) jprint(8) = 1
      if(jprint(48).ne.0) jprint(8) = 1
      if(jprint(49).ne.0) jprint(8) = 1
      if(jprint(54).ne.0) jprint(8) = 1
c
c  take care of the plot option:
c
      if(jprint(43).ne.0) then
                            jprint(8)  =  1
        if(jprint(22).eq.0) jprint(22) =  lfnao
        if(iwtnao.eq.0)     iwtnao     = -lfnnao
        if(jprint(28).eq.0) jprint(28) = -lfnnho
        if(iwtnbo.eq.0)     iwtnbo     = -lfnnbo
        if(jprint(23).eq.0) jprint(23) = -lfnnlm
        if(jprint(26).eq.0) jprint(26) = -lfnmo
        if(jprint(27).eq.0) jprint(27) = -lfndm
        if(jprint(44).eq.0) jprint(44) = -lfnpna
        if(jprint(30).eq.0) jprint(30) = -lfnpnh
        if(jprint(25).eq.0) jprint(25) = -lfnpnb
        if(jprint(49).eq.0) jprint(49) = -lfnpnl
      end if
c
c  print hybrids if the nbo output is requested:
c
      iwhybs = jprint(5)
      return
c
c  abort program: unrecognizable keyword encountered
c
 4800 write(lfnpr,9800) (keywd(i),i=1,6)
        stop
c
c  incompatible job options have been requested:
c
 4900 continue
        write(lfnpr,9900)
        stop
c
 9800 format(1x,'error: unrecognizable keyword  >>',6a1,'<<',/,1x,
     *          'program must halt.')
 9900 format(1x,'the nbo program must stop because the options /mulat/',
     + ' and /mulorb/',/1x,'currently require the ao bond order matrix',
     + ', rather than the ao density',/1x,'matrix.  the program could ',
     + 'be modified to permit this.')
      end
c*****************************************************************************
      subroutine nbodim(memory)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      dimension nspdfg(5,2)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbbas/label(maxbas,6),lval(maxbas),imval(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbao/lctr(maxbas),lang(maxbas)
c
      data iread/4hread/
c
c  nbodim:  set up dimensioning information, lists in common/nbatom/,
c           and compare storage needs with amount of storage available
c
c  find:
c     mxaolm, the maximum number of atomic orbitals of the same symmetry
c              on a single atom
c     mxao, the maximum number of atomic orbitals per atom
c     mxbo, the maximum number of atomic orbitals per two-center or
c              three-center bond
c
      do 300 i = 1,nbas
        lm = lang(i)
        lval(i) = lm/100
        im = lm - lval(i)*100
        if(im.gt.50) im = im - 50
        imval(i) = im
  300 continue
c
      mxao   = 0
      mxao2  = 0
      mxao3  = 0
      mxaolm = 0
      llu = 0
      do 500 i = 1,natoms
        n = 0
        do 400 il = 1,5
          do 400 ityp = 1,2
  400       nspdfg(il,ityp) = 0
        do 410 j = 1,nbas
          if(lctr(j).ne.i) go to 410
          lm = lang(j)
          l = lm/100
          im = lm - l*100
c
c  if im.ne.1 (that is, if this is not the first component of the
c  ang. mom. l functions on the atom), don't count it in nspdfg:
c
          if(im.ne.1) go to 410
c
c  ityp=1 for cartesian function, =2 for true spherical harmonic:
c
          ityp = 1
          if(im.gt.50) ityp = 2
          il = l + 1
          nspdfg(il,ityp) = nspdfg(il,ityp)+1
  410     if(lctr(j).eq.i) n = n + 1
c
c  number of s orbitals= no. s orbs input + no. cartesian d and g orbs:
c
        nspdfg(1,1) = nspdfg(1,1) + nspdfg(1,2) + nspdfg(3,1) +
     +                nspdfg(5,1)
c
c  number of p orbitals= no. p orbs input + no. cartesian f orbs:
c
        nspdfg(2,1) = nspdfg(2,1) + nspdfg(2,2) + nspdfg(4,1)
c
c  number of d orbitals= no. d orbs input + no. cartesian g orbs:
c
        nspdfg(3,1) = nspdfg(3,1) + nspdfg(3,2) + nspdfg(5,1)
c
c  number of f orbitals:
c
        nspdfg(4,1) = nspdfg(4,1) + nspdfg(4,2)
c
c  number of g orbitals:
c
        nspdfg(5,1) = nspdfg(5,1) + nspdfg(5,2)
c
        do 430 il = 1,5
          if(nspdfg(il,1).le.mxaolm) go to 430
          mxaolm = nspdfg(il,1)
  430   continue
c
        norbs(i) = n
        ll(i) = llu + 1
        lu(i) = ll(i) + n - 1
        llu = lu(i)
        if(n.le.mxao) go to 460
        mxao3 = mxao2
        mxao2 = mxao
        mxao = n
        go to 500
  460   if(n.le.mxao2) go to 480
        mxao3 = mxao2
        mxao2 = n
        go to 500
  480   if(n.le.mxao3) go to 500
        mxao3 = n
  500 continue
      mxbo = mxao + mxao2
      if(iw3c.eq.1) mxbo = mxbo + mxao3
c
c  compute storage requirements and compare with available core space:
c
c  storage for density matrix (dm) and transformations (t):
c
      need0 = 2*ndim*ndim
c
c  compute storage for natural population analysis:
c
      need1 = 0
      io = ioinqr(iwtnao)
      if(io.ne.iread.and..not.ortho) then
        need  = ndim + ndim + ndim*ndim + mxaolm*mxaolm + ndim
     +        + mxaolm*mxaolm + mxaolm*mxaolm + ndim + 9*mxaolm
        need1 = max(need1,need)
      end if
c
      need  = natoms*natoms + natoms + natoms*natoms + natoms*natoms +
     +        ndim*ndim + ndim
      need1 = max(need1,need)
c
      need  = natoms*natoms + ndim*ndim + ndim
      need1 = max(need1,need)
c
      if(jprint(9).ne.0) then
        need  = natoms*natoms + ndim*ndim + ndim*ndim + ndim*(ndim+5)
        need1 = max(need1,need)
      end if
c
      need1 = need1 + need0
c
c  compute storage for natural bond orbital analysis:
c
      need2 = 0
      if(jprint(1).eq.0) then
        if(ioinqr(iwtnab).ne.iread) then
          need  = natoms*natoms + ndim + 3*ndim + mxao*ndim + ndim
     +          + mxbo*mxbo + mxbo*mxbo + mxbo + mxbo + mxao*mxao
     +          + mxao*mxao + mxao + mxao + mxao + natoms*natoms
        else
          need  = natoms*natoms + ndim + 3*ndim
        end if
        need2 = max(need2,need)
c
        if(.not.ortho) then
          need  = natoms*natoms + 4*ndim*ndim + mxao + 3*ndim
          need2 = max(need2,need)
        end if
c
        need  = natoms*natoms + ndim + mxao + ndim*ndim + ndim*ndim
     +        + ndim + ndim
        need2 = max(need2,need)
c
        need  = natoms*natoms + ndim + ndim + ndim + ndim*ndim
        need2 = max(need2,need)
c
        if(jprint(36).ne.0) then
          need  = natoms*natoms + ndim + 3*natoms + ndim*ndim
     +          + ndim*ndim + ndim
          need2 = max(need2,need)
        end if
c
        need  = natoms*natoms + ndim + ndim*ndim + ndim*ndim
     +        + ndim*(ndim+5)
        need2 = max(need2,need)
c
        if(jprint(6).ne.0) then
          need  = natoms*natoms + ndim + ndim*ndim + ndim + natoms
     +          + ndim
          need2 = max(need2,need)
        end if
c
c  compute storage for natural localized molecular orbital analysis:
c
        need3 = 0
        if(jprint(8).ne.0) then
          need  = natoms*natoms + ndim + ndim + ndim*ndim + ndim*ndim
          need3 = max(need3,need)
c
          need  = ndim + ndim + ndim + natoms*natoms + 2*natoms*natoms
     +          + ndim*natoms + ndim*natoms*(natoms-1)/2 + ndim*ndim
          need3 = max(need3,need)
c
          need  = natoms*natoms + ndim*ndim + ndim*ndim + ndim*(ndim+5)
          need3 = max(need3,need)
c
          if(jprint(46).ne.0) then
            need  = ndim*ndim + ndim*ndim + ndim*ndim + ndim*ndim
     +            + ndim*ndim + ndim*ndim + ndim + natoms*natoms
            need3 = max(need3,need)
          end if
        end if
      end if
c
c  print scratch storage requirements:
c
      if(iprint.ge.0) then
        if(jprint(1).eq.0) then
          if(jprint(8).ne.0) then
            write(lfnpr,3300) need1,need2,need3,memory
          else
            need3 = 0
            write(lfnpr,3200) need1,need2,memory
          end if
        else
          need2 = 0
          need3 = 0
          write(lfnpr,3100) need1,memory
        end if
      end if
      if(need1.gt.memory.or.need2.gt.memory.or.need3.gt.memory) goto 990
      return
c
  990 write(lfnpr,9900)
      stop
c
 3100 format(/1x,'storage needed:',i9,' in npa (',i9,' available)')
 3200 format(/1x,'storage needed:',i9,' in npa,',i9,' in nbo (',i9,
     + ' available)')
 3300 format(/1x,'storage needed:',i9,' in npa,',i9,' in nbo,',i9,
     + ' in nlmo (',i9,' available)')
 9900 format(/1x,'*** not enough core storage is available ***'/)
      end
c**************************************************************************
c
c  nao/nbo/nlmo formation routines: (called by sr nbo)
c
c      subroutine naodrv(dm,t,a)
c      subroutine naosim(dm,t,a)
c      subroutine dmnao(dm,t,a)
c      subroutine dmsim(dm,t,a)
c      subroutine nbodrv(dm,t,a,memory)
c
c**************************************************************************
      subroutine naodrv(dm,t,a)
c**************************************************************************
      implicit real*8 (a-h,o-z)
c
c  driver subroutine to calculate natural atomic orbitals (naos)
c  given 1-particle density matrix in an arbitrary atom-centered
c  atomic orbital basis set.
c
c        t = overlap matrix for the primitive ao basis
c             (on return, this is the ao to nao transformation matrix)
c       dm = density matrix in the primitive ao basis
c               (or bond-order matrix, if iwdm = 1)
c
c   the spin nature of dm is indicated by:
c    ispin =  0: spinless  (closed shell)
c    ispin = +2: alpha spin
c    ispin = -2: spin
c   (ispin is the reciprocal of the s(z) quantum no.)
c
c        a = scratch storage from the main program.  the location of a(1)
c               is in the common block /scm/ in the main program,
c               after the storage for the matrices 's','dm'
c             ('a' is the vector which is partitioned
c                   according to the storage needs of each program run)
c     atom, basis, option, nbinfo: common blocks with data transfered from
c        from the input programs.
c
c-----------------------------------------------------------------------------
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nbao/lctr(maxbas),lang(maxbas)
c
      dimension t(ndim,ndim),dm(ndim,ndim),a(1)
      character*80 title
c
      data one/1.0d0/
      data iprnt,iwrit,iread/4hprnt,4hwrit,4hread/
c
c  form labels for the raw ao basis set:
c
      call lblao
c
c  copy the ao centers and labels from /nbao/ to /nbbas/:
c
      do 5 i = 1,nbas
        lbl(i) = lctr(i)
        lorbc(i) = lang(i)
    5 continue
c
c  write out the ao basis set information:
c
      if(jprint(22).gt.0) then
        call wrbas(a,a,jprint(22))
      end if
c
c  write out the archive file:
c
      if(jprint(7).ne.0) then
        call wrarc(a,a,jprint(7))
      end if
c
c  output the ao overlap matrix:
c
      io = ioinqr(jprint(39))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        title = 'ao overlap matrix:'
        call aout(t,ndim,-nbas,nbas,title,1,jprint(39))
      end if
c
c  output the ao-mo transformation matrix:
c
      io = ioinqr(jprint(26))
      if(.not.open.and.(io.eq.iprnt.or.io.eq.iwrit)) then
        call feaomo(a,it)
        if(it.ne.0) then
          title = 'mos in the ao basis:'
          call aout(a,ndim,nbas,nbas,title,1,jprint(26))
        end if
      end if
c
c  output the ao fock matrix:
c
      io = ioinqr(jprint(40))
      if(.not.open.and.(io.eq.iprnt.or.io.eq.iwrit)) then
        call fefao(a,iwfock)
        if(iwfock.ne.0) then
          title = 'ao fock matrix:'
          call aout(a,ndim,-nbas,nbas,title,1,jprint(40))
        end if
      end if
c
c  output the ao bond-order matrix:
c
      io = ioinqr(jprint(42))
      if(iwdm.eq.1.and.(io.eq.iprnt.or.io.eq.iwrit)) then
        title = 'spinless ao bond-order matrix:'
        call aout(dm,ndim,-nbas,nbas,title,1,jprint(42))
      end if
c
c  convert the bond-order matrix to the density matrix:
c
      if(iwdm.ne.0) call simtrm(dm,t,a,ndim,nbas,iwmulp,iwcubf)
c
c  output the ao density matrix:
c
      io = ioinqr(jprint(27))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        title = 'spinless ao density matrix:'
        call aout(dm,ndim,-nbas,nbas,title,1,jprint(27))
      end if
c
c  output the ao dipole matrices:
c
      io = ioinqr(jprint(50))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        ix = 1
        call fedxyz(a,ix)
        if(ix.ne.0) then
          title = 'ao x dipole integrals:'
          call aout(a,ndim,-nbas,nbas,title,1,jprint(50))
        end if
        ix = 2
        call fedxyz(a,ix)
        if(ix.ne.0) then
          title = 'ao y dipole integrals:'
          call aout(a,ndim,-nbas,nbas,title,1,jprint(50))
        end if
        ix = 3
        call fedxyz(a,ix)
        if(ix.ne.0) then
          title = 'ao z dipole integrals:'
          call aout(a,ndim,-nbas,nbas,title,1,jprint(50))
        end if
      end if
c
c  allocate scratch common for nao routines:
c
c  a(i1) = v(ndim)  (also used for guide(natoms,natoms))
c  a(i2) = renorm(ndim)
c  a(i3) = blk(ndim,ndim)
c  a(i4) = sblk(mxaolm,mxaolm)
c  a(i5) = eval(ndim)
c  a(i6) = c(mxaolm,mxaolm)
c  a(i7) = evect(mxaolm,mxaolm)
c  a(i8) = eval2(ndim)
c  leave this last in the list since these are integers:
c  a(i9) = listao(mxaolm,9)
c
      nblock = mxaolm*mxaolm
      i1   = 1
      i2   = i1 + ndim
      i3   = i2 + ndim
      i4   = i3 + ndim*ndim
      i5   = i4 + nblock
      i6   = i5 + ndim
      i7   = i6 + nblock
      i8   = i7 + nblock
      i9   = i8 + ndim
c     iend = i9 + 9*mxaolm
c
c  read in t-nao, nao labels, the pnao overlap matrix, and compute the
c  nao density matrix: (note that t contains the pnao overlap matrix
c  after rdtnao is called)
c
      if(ioinqr(iwtnao).eq.iread) then
        call rdtnao(dm,t,a(i1),iwtnao)
        go to 580
      end if
c
c  transform all sets of cartesian d,f,g orbitals, and relabel all orbitals:
c
      call dfgorb(a(i2),dm,t,ictran,iwcubf,0,lfnpr)
c
c  store pure ao density matrix in scratch storage:
c
      call svppao(dm)
c
c  consolidate density matrix and overlap matrix in dm:
c
      call consol(dm,t,ndim,nbas)
c
c  find natural atomic orbital basis set transformation t from dm:
c  (upon return, dm contains the full nao density matrix)
c
      call nao(t,dm,a(i1),a(i3),a(i4),a(i5),a(i6),a(i7),a(i8),a(i9),
     *         nblock)
c
c  if d orbitals were transformed, transform the nao transformation t
c  so that t is the transform from the original ao's to the nao's:
c
      if(ictran.ne.0) call dfgorb(a(i2),dm,t,idtran,iwcubf,1,lfnpr)
c
c  save tnao for later use:
c
      call svtnao(t)
c
c  if d orbitals were transformed, transform the pnao transformation
c  so that it is the transform from the original ao's to the pnao's:
c
      call fepnao(a(i3))
c
c  for case that rpnaos are written to disk, set occupancy weights to -1
c  as a signal that they should be recomputed:
c
      do 260 i = 0,nbas-1
  260   a(i4+i) = -one
c
      if(ictran.ne.0) call dfgorb(a(i2),dm,a(i3),idtran,iwcubf,1,lfnpr)
c
c  compute non-orthogonal nao overlap matrix, spnao:
c
      call fesraw(t)
      call simtrs(t,a(i3),a(i4),ndim,nbas)
      call svsnao(t)
c
c  write t-nao, nao labels, and the pnao overlap matrix:
c
      if(ioinqr(iwtnao).eq.iwrit) call wrtnao(t,iwtnao)
c
c  dm is now the density matrix in the nao basis
c  t is the non-orthogonal pnao overlap matrix  (!!!)
c
  580 continue
      i1   = 1
      i2   = i1 + natoms*natoms
      i3   = i2 + natoms
      i4   = i3 + natoms*natoms
      i5   = i4 + natoms*natoms
      i6   = i5 + ndim*ndim
c     iend = i6 + ndim
      call naoanl(dm,t,a(i1),a(i2),a(i3),a(i4),a(i5),a(i6))
c
c  do not destroy the matrix at a(i1).  this holds the wiberg bond
c  index which needs to be passed to the nbo routines.
c
c  save the nao density matrix:
c
      call svdnao(dm)
c
c  form the nao labels:
c
      call lblnao
c
c  reorganize the scratch vector:
c
      i1   = 1
      i2   = i1 + natoms*natoms
c     iend = i2 + ndim*ndim
c
c  output the ao-pnao transformation matrix:
c
      io = ioinqr(jprint(44))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fepnao(t)
        title = 'pnaos in the ao basis:'
        call aout(t,ndim,nbas,nbas,title,1,jprint(44))
      end if
c
c  output the pnao overlap matrix:
c
      io = ioinqr(jprint(19))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fesnao(a(i2))
        title = 'pnao overlap matrix:'
        call aout(a(i2),ndim,-nbas,nbas,title,2,jprint(19))
      end if
c
c  fetch the ao-nao transformation from the nbo daf:
c
      call fetnao(t)
c
c  print the ao-nao transformation matrix:
c
      if(ioinqr(iwtnao).eq.iprnt) then
        title = 'naos in the ao basis:'
        call aout(t,ndim,nbas,nbas,title,1,iwtnao)
      end if
c
c  output the nao dipole matrices:
c
      io = ioinqr(jprint(51))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        ix = 1
        call fedxyz(a(i2),ix)
        if(ix.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nao x dipole integrals:'
          call aout(a(i2),ndim,-nbas,nbas,title,2,jprint(51))
        end if
        ix = 2
        call fedxyz(a(i2),ix)
        if(ix.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nao y dipole integrals:'
          call aout(a(i2),ndim,-nbas,nbas,title,2,jprint(51))
        end if
        ix = 3
        call fedxyz(a(i2),ix)
        if(ix.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nao z dipole integrals:'
          call aout(a(i2),ndim,-nbas,nbas,title,2,jprint(51))
        end if
      end if
c
c  if this is an open shell wavefunction, don't do anything more:
c
      if(open) return
c
c  output the nao-mo transformation matrix:
c
      io = ioinqr(jprint(9))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        i1   = 1
        i2   = i1 + natoms*natoms
        i3   = i2 + ndim*ndim
        i4   = i3 + ndim*ndim
c       iend = i4 + ndim*(ndim+5)
        call frmtmo(t,a(i2),a(i3),a(i4),2,jprint(9))
      end if
c
c  reorganize the scratch vector:
c
      i1   = 1
      i2   = i1 + natoms*natoms
      i3   = i2 + ndim*ndim
c     iend = i3 + ndim
c
c  output the nao fock matrix:
c
      io = ioinqr(jprint(31))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fefao(a(i2),iwfock)
        if(iwfock.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nao fock matrix:'
          call aout(a(i2),ndim,-nbas,nbas,title,2,jprint(31))
        end if
      end if
c
c  output the nao density matrix:
c
      io = ioinqr(jprint(35))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        title = 'nao density matrix:'
        call aout(dm,ndim,-nbas,nbas,title,2,jprint(35))
      end if
      return
      end
c*****************************************************************************
      subroutine naosim(dm,t,a)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbao/lctr(maxbas),lang(maxbas)
c
      dimension dm(ndim,ndim),t(ndim,ndim),a(1)
      character*80 title
c
      data zero,one/0.0d0,1.0d0/
      data iprnt,iwrit/4hprnt,4hwrit/
c
c  this routine simulates the action of the nao subprogram:
c
c  form labels for the raw ao basis set:
c
      call lblao
c
c  copy the ao centers and labels from /nbao/ to /nbbas/:
c
      do 5 i = 1,nbas
        lbl(i) = lctr(i)
        lorbc(i) = lang(i)
    5 continue
c
c  write out the ao basis set information:
c
      if(jprint(22).gt.0) then
        call wrbas(a,a,jprint(22))
      end if
c
c  write out the archive file:
c
      if(jprint(7).ne.0) then
        call wrarc(a,a,jprint(7))
      end if
c
c  output the ao-mo transformation matrix:
c
      io = ioinqr(jprint(26))
      if(.not.open.and.(io.eq.iprnt.or.io.eq.iwrit)) then
        call feaomo(a,it)
        if(it.ne.0) then
          title = 'mos in the ao basis:'
          call aout(a,ndim,nbas,nbas,title,1,jprint(26))
        end if
      end if
c
c  output the ao fock matrix:
c
      io = ioinqr(jprint(40))
      if(.not.open.and.(io.eq.iprnt.or.io.eq.iwrit)) then
        call fefao(a,iwfock)
        if(iwfock.ne.0) then
          title = 'ao fock matrix:'
          call aout(a,ndim,-nbas,nbas,title,1,jprint(40))
        end if
      end if
c
c  output the ao density matrix:
c
      io = ioinqr(jprint(27))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        title = 'spinless ao density matrix:'
        call aout(dm,ndim,-nbas,nbas,title,1,jprint(27))
      end if
c
c  output the ao dipole matrices:
c
      io = ioinqr(jprint(50))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        ix = 1
        call fedxyz(a,ix)
        if(ix.ne.0) then
          title = 'ao x dipole integrals:'
          call aout(a,ndim,-nbas,nbas,title,1,jprint(50))
        end if
        ix = 2
        call fedxyz(a,ix)
        if(ix.ne.0) then
          title = 'ao y dipole integrals:'
          call aout(a,ndim,-nbas,nbas,title,1,jprint(50))
        end if
        ix = 3
        call fedxyz(a,ix)
        if(ix.ne.0) then
          title = 'ao z dipole integrals:'
          call aout(a,ndim,-nbas,nbas,title,1,jprint(50))
        end if
      end if
c
c  initialize the ao to nao transformation matrix (unit matrix):
c
      do 20 j = 1,nbas
        do 10 i = 1,nbas
          t(i,j) = zero
   10   continue
        t(j,j) = one
   20 continue
c
c  save tnao for later use:
c
      call svtnao(t)
c
c  fill atomic orbital information lists:
c
      do 30 i = 1,nbas
        naoctr(i) = lctr(i)
        naol(i)   = lang(i)
        lstocc(i) = 1
   30 continue
c
c  perform the natural population analysis: (note that routine naoanl
c  expects to find the overlap matrix in t, which is the unit matrix
c  for orthogonal basis sets. upon return from naoanl, t is the ao to
c  nao transformation, which is still a unit matrix):
c
      i1   = 1
      i2   = i1 + natoms*natoms
      i3   = i2 + natoms
      i4   = i3 + natoms*natoms
      i5   = i4 + natoms*natoms
      i6   = i5 + ndim*ndim
c     iend = i6 + ndim
      call naoanl(dm,t,a(i1),a(i2),a(i3),a(i4),a(i5),a(i6))
c
c  do not destroy the matrix at a(i1).  this holds the wiberg bond
c  index which needs to be passed to the nbo routines.
c
c  save the nao density matrix:
c
      call svdnao(dm)
c
c  form the nao labels:
c
      call lblnao
c
c  if this is an open shell wavefunction, don't do anything more:
c
      if(open) return
c
c  output the nao-mo transformation matrix:
c
      io = ioinqr(jprint(9))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        i1   = 1
        i2   = i1 + natoms*natoms
        i3   = i2 + ndim*ndim
        i4   = i3 + ndim*ndim
c       iend = i4 + ndim*(ndim+5)
        call frmtmo(t,a(i2),a(i3),a(i4),2,jprint(9))
      end if
c
c  reorganize the scratch vector:
c
      i1   = 1
      i2   = i1 + natoms*natoms
      i3   = i2 + ndim*ndim
c     iend = i3 + ndim
c
c  output the nao fock matrix:
c
      io = ioinqr(jprint(31))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fefao(a(i2),iwfock)
        if(iwfock.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nao fock matrix:'
          call aout(a(i2),ndim,-nbas,nbas,title,2,jprint(31))
        end if
      end if
c
c  output the nao density matrix:
c
      io = ioinqr(jprint(35))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        title = 'nao density matrix:'
        call aout(dm,ndim,-nbas,nbas,title,2,jprint(35))
      end if
      return
      end
c**************************************************************************
      subroutine dmnao(dm,t,a)
c**************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nbao/lctr(maxbas),lang(maxbas)
      common/nbnao/naoc(maxbas),naoa(maxbas),ltyp(maxbas),iprin(maxbas)
c
      dimension dm(ndim,ndim),t(ndim,ndim),a(1)
      character*80 title
c
      data iprnt,iwrit/4hprnt,4hwrit/
c
c  place alpha or beta occupation matrix in dm and transform from the ao
c  to nao basis:
c
      if(alpha) then
        if(jprint(4).ne.0) write(lfnpr,2100)
      else
        do 70 i = 1,nbas
          naoctr(i) = naoc(i)
          naol(i) = naoa(i)
          lbl(i) = lctr(i)
          lorbc(i) = lang(i)
   70   continue
        call fetnao(t)
        if(jprint(4).ne.0) write(lfnpr,2200)
      end if
c
c  output the ao-mo transformation matrix:
c
      io = ioinqr(jprint(26))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call feaomo(a,it)
        if(it.ne.0) then
          title = 'mos in the ao basis:'
          call aout(a,ndim,nbas,nbas,title,1,jprint(26))
        end if
      end if
c
c  output the ao fock matrix:
c
      io = ioinqr(jprint(40))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fefao(a,iwfock)
        if(iwfock.ne.0) then
          title = 'ao fock matrix:'
          call aout(a,ndim,-nbas,nbas,title,1,jprint(40))
        end if
      end if
c
c  fetch alpha or beta dm (according to whether alpha or beta is true):
c
      call fedraw(dm,a)
c
c  output the ao bond-order matrix:
c
      io = ioinqr(jprint(42))
      if(iwdm.ne.0.and.(io.eq.iprnt.or.io.eq.iwrit)) then
        title = 'ao bond-order matrix:'
        call aout(dm,ndim,-nbas,nbas,title,1,jprint(42))
      end if
c
c  convert the bond-order matrix to the density matrix:
c
      if(iwdm.ne.0) then
        i1   = 1
        i2   = i1 + ndim*ndim
c       iend = i2 + ndim*ndim
        call fesraw(a(i1))
        call simtrm(dm,a(i1),a(i2),ndim,nbas,iwmulp,iwcubf)
      end if
c
c  output the ao density matrix:
c
      io = ioinqr(jprint(27))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        title = 'ao density matrix:'
        call aout(dm,ndim,-nbas,nbas,title,1,jprint(27))
      end if
c
c  transform dm to the nao basis:
c
      call simtrs(dm,t,a,ndim,nbas)
c
c  save the nao density matrix in scratch storage:
c
      call svdnao(dm)
c
c  print the natural population analysis for this spin case:
c
      i1   = 1
      i2   = i1 + natoms*natoms
      i3   = i2 + natoms
      i4   = i3 + natoms*natoms
      i5   = i4 + natoms*natoms
      i6   = i5 + ndim*ndim
c     iend = i6 + ndim
      call fesnao(t)
      call naoanl(dm,t,a(i1),a(i2),a(i3),a(i4),a(i5),a(i6))
c
c  note: do not destroy the wiberg bond index which is stored in the first
c  natoms*natoms elements of the scratch vector a.  this is matrix is
c  required for the nbo analysis:
c
c  note that t is now t-ao-nao:
c
c  form the nao labels:
c
      call lblnao
c
c  output the nao-mo transformation matrix:
c
      io = ioinqr(jprint(9))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        i1   = 1
        i2   = i1 + natoms*natoms
        i3   = i2 + ndim*ndim
        i4   = i3 + ndim*ndim
c       iend = i4 + ndim*(ndim+5)
        call frmtmo(t,a(i2),a(i3),a(i4),2,jprint(9))
      end if
c
c  reorganize the scratch vector:
c
      i1   = 1
      i2   = i1 + natoms*natoms
      i3   = i2 + ndim*ndim
c     iend = i3 + ndim
c
c  output the nao fock matrix:
c
      io = ioinqr(jprint(31))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fefao(a(i2),iwfock)
        if(iwfock.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nao fock matrix:'
          call aout(a(i2),ndim,-nbas,nbas,title,2,jprint(31))
        end if
      end if
c
c  output the nao density matrix:
c
      io = ioinqr(jprint(35))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        title = 'nao density matrix:'
        call aout(dm,ndim,-nbas,nbas,title,2,jprint(35))
      end if
      return
c
 2100 format(//1x,
     * '***************************************************',/1x,
     * '*******         alpha spin orbitals         *******',/1x,
     * '***************************************************')
 2200 format(//1x,
     * '***************************************************',/1x,
     * '*******         beta  spin orbitals         *******',/1x,
     * '***************************************************')
      end
c**************************************************************************
      subroutine dmsim(dm,t,a)
c**************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nbao/lctr(maxbas),lang(maxbas)
      common/nbnao/naoc(maxbas),naoa(maxbas),ltyp(maxbas),iprin(maxbas)
c
      dimension dm(ndim,ndim),t(ndim,ndim),a(1)
      character*80 title
c
      data iprnt,iwrit/4hprnt,4hwrit/
c
c  simulate the alpha/beta nao subprogram:
c
      if(alpha) then
        if(jprint(4).ne.0) write(lfnpr,2100)
      else
        do 70 i = 1,nbas
          naoctr(i) = naoc(i)
          naol(i) = naoa(i)
          lbl(i) = lctr(i)
          lorbc(i) = lang(i)
   70   continue
        call fetnao(t)
        if(jprint(4).ne.0) write(lfnpr,2200)
      end if
c
c  output the ao-mo transformation matrix:
c
      io = ioinqr(jprint(26))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call feaomo(a,it)
        if(it.ne.0) then
          title = 'mos in the ao basis:'
          call aout(a,ndim,nbas,nbas,title,1,jprint(26))
        end if
      end if
c
c  output the ao fock matrix:
c
      io = ioinqr(jprint(40))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fefao(a,iwfock)
        if(iwfock.ne.0) then
          title = 'ao fock matrix:'
          call aout(a,ndim,-nbas,nbas,title,1,jprint(40))
        end if
      end if
c
c  fetch alpha or beta dm (according to whether alpha or beta is true):
c
      call fedraw(dm,a)
c
c  output the ao density matrix:
c
      io = ioinqr(jprint(27))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        title = 'ao density matrix:'
        call aout(dm,ndim,-nbas,nbas,title,1,jprint(27))
      end if
c
c  save the nao density matrix in scratch storage:
c
      call svdnao(dm)
c
c  print the natural population analysis for this spin case:
c
      i1   = 1
      i2   = i1 + natoms*natoms
      i3   = i2 + natoms
      i4   = i3 + natoms*natoms
      i5   = i4 + natoms*natoms
      i6   = i5 + ndim*ndim
c     iend = i6 + ndim
      call naoanl(dm,t,a(i1),a(i2),a(i3),a(i4),a(i5),a(i6))
c
c  note: do not destroy the wiberg bond index which is stored in the first
c  natoms*natoms elements of the scratch vector a.  this is matrix is
c  required for the nbo analysis:
c
c  form the nao labels:
c
      call lblnao
c
c  output the nao-mo transformation matrix:
c
      io = ioinqr(jprint(9))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        i1   = 1
        i2   = i1 + natoms*natoms
        i3   = i2 + ndim*ndim
        i4   = i3 + ndim*ndim
c       iend = i4 + ndim*(ndim+5)
        call frmtmo(t,a(i2),a(i3),a(i4),2,jprint(9))
      end if
c
c  reorganize the scratch vector:
c
      i1   = 1
      i2   = i1 + natoms*natoms
      i3   = i2 + ndim*ndim
c     iend = i3 + ndim
c
c  output the nao fock matrix:
c
      io = ioinqr(jprint(31))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fefao(a(i2),iwfock)
        if(iwfock.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nao fock matrix:'
          call aout(a(i2),ndim,-nbas,nbas,title,2,jprint(31))
        end if
      end if
c
c  print the nao density matrix:
c
      io = ioinqr(jprint(35))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        title = 'nao density matrix:'
        call aout(dm,ndim,-nbas,nbas,title,2,jprint(35))
      end if
      return
c
 2100 format(//1x,
     * '***************************************************',/1x,
     * '*******         alpha spin orbitals         *******',/1x,
     * '***************************************************')
 2200 format(//1x,
     * '***************************************************',/1x,
     * '*******         beta  spin orbitals         *******',/1x,
     * '***************************************************')
      end
c**************************************************************************
      subroutine nbodrv(dm,t,a,memory)
c**************************************************************************
c
c  driver subroutine to calculate natural hybrid orbitals (nhos) and
c  natural bond orbitals (nbos) from the density matrix in the nao basis
c
c        t = scratch storage
c       dm = nao density matrix
c            the spin nature of dm is indicated by:
c            ispin =  0: spinless  (closed shell)
c            ispin = +2: alpha spin
c            ispin = -2: spin
c            (ispin is the reciprocal of the s(z) quantum no.)
c
c        a = scratch storage from the main program.  the location of a(1)
c               is in the common block /scm/ in the main program,
c               after the storage for the matrices 's','dm'
c             ('a' is the vector which is partitioned
c                   according to the storage needs of each program run)
c     atom, basis, option, nbinfo: common blocks with data transfered from
c        from the input programs.
c
c-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*80 title
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbmol/nmolec,molat(maxatm),molec(maxatm,maxatm),
     +              nmola,molata(maxatm),moleca(maxatm,maxatm)
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
c
      dimension t(ndim,ndim),dm(ndim,ndim),a(1)
c
      data iprnt,iwrit,iread/4hprnt,4hwrit,4hread/
      data zero/0.0d0/
c
c  skip nbo transformation if requested:
c
      if(jprint(1).gt.0) then
        write(lfnpr,2000)
        return
      end if
c
c  organize scratch storage vector a.  warning:  this is redefined
c  several times after the nbos are formed
c
c  a(i0)  = guide(natoms,natoms)
c  a(i1)  = bndocc(ndim)
c  a(i2)  = pol(ndim,3)
c  a(i3)  = q(mxao,ndim)
c  a(i4)  = v(ndim)
c  a(i5)  = blk(mxbo,mxbo)
c  a(i6)  = c(mxbo,mxbo)
c  a(i7)  = eval(mxbo)
c  a(i8)  = borb(mxbo)
c  a(i9)  = p(mxao,mxao)
c  a(i10) = pk(mxao,mxao)
c  a(i11) = hyb(mxao)
c  a(i12) = va(mxao)
c  a(i13) = vb(mxao)
c
      i0   = 1
      i1   = i0  + natoms*natoms
      i2   = i1  + ndim
      i3   = i2  + 3*ndim
      i4   = i3  + mxao*ndim
      i5   = i4  + ndim
      i6   = i5  + mxbo*mxbo
      i7   = i6  + mxbo*mxbo
      i8   = i7  + mxbo
      i9   = i8  + mxbo
      i10  = i9  + mxao*mxao
      i11  = i10 + mxao*mxao
      i12  = i11 + mxao
      i13  = i12 + mxao
      i14  = i13 + mxao
c     iend = i14 + natoms*natoms
c
      if(jprint(5).ne.0.and.ispin.eq.0) write(lfnpr,1400)
      if(jprint(5).ne.0.and.ispin.eq.2) write(lfnpr,1410)
      if(jprint(5).ne.0.and.ispin.eq.-2) write(lfnpr,1420)
c
c  read in t-nab, label, ibxm, transform dm, and find bndocc if iwtnab=iread:
c
      if(ioinqr(iwtnab).eq.iread) then
        call rdtnab(t,dm,a(i1),a(i2),iwtnab)
      else
c
c  search input file for $core input:
c
        if(.not.beta) then
          call corinp(jprint(2),jcore)
          call rdcore(jcore)
        end if
c
c  search input file for $choose input:
c
        if(.not.beta) then
          call chsinp(jprint(2),ichoos)
          if(open.and.ichoos.eq.1.and.jprint(32).ne.0) then
            write(lfnpr,1390)
            ichoos = 0
          end if
        end if
c
c  calculate natural hybrid orbitals and bond orbitals:
c
        if(ichoos.ne.1) call nathyb(dm,t,a(i0),a(i1),a(i2),a(i3),a(i4),
     +                            a(i5),a(i6),a(i7),a(i8),a(i9),a(i10),
     +                            a(i11),a(i12),a(i13),a(i14))
        if(ichoos.eq.1) call chsdrv(dm,t,a(i0),a(i1),a(i2),a(i3),a(i4),
     +                            a(i5),a(i6),a(i7),a(i8),a(i9),a(i10),
     +                            a(i11),a(i12),a(i13),a(i14))
c
c  if nbo search was abandoned, don't try to do anything further:
c
        if(jprint(1).lt.0) return
c
c  sort the nbos by atom:
c
        call srtnbo(t,a(i1))
c
c  form the nbo density matrix:
c
        call simtrs(dm,t,a(i2),ndim,nbas)
c
c  check nho overlaps to see if bond orbitals should be relabelled:
c
        if(.not.ortho) then
          i0   = 1
          i1   = i0 + natoms*natoms
          i2   = i1 + ndim
          i3   = i2 + mxao
          i4   = i3 + ndim*ndim
          i5   = i4 + ndim*ndim
          i6   = i5 + ndim
c         iend = i6 + ndim
          call xcited(dm,t,a(i2),a(i3),a(i4),a(i5),a(i6),a(i6))
        end if
      end if
c
c  t  now contains the nao-nbo transformation matrix
c  dm now contains the nbo density matrix
c  a(i0)  contains the wiberg bond index matrix      ! don't destroy this
c  a(i1)  contains the nbo occupancies               ! don't destroy this
c  a(i2)  is scratch space
c
c  save the nao-nbo transformation on the nbo daf:
c
      call svtnab(t)
c
c  form the nbo labels:
c
      call lblnbo
c
c  write out the analysis of bond orbitals:
c
      i0   = 1
      i1   = i0 + natoms*natoms
      i2   = i1 + ndim
      i3   = i2 + ndim
      i4   = i3 + ndim
c     iend = i4 + ndim*ndim
      call anlyze(t,a(i1),a(i2),a(i3),a(i4))
c
c  write out hybrid directionality and bond bending info:
c
      if(jprint(36).ne.0) then
        i0   = 1
        i1   = i0 + natoms*natoms
        i2   = i1 + ndim
        i3   = i2 + 3*natoms
        i4   = i3 + ndim*ndim
        i5   = i4 + ndim*ndim
c       iend = i5 + ndim
        call hybdir(a(i1),a(i2),a(i3),a(i4),a(i5))
      end if
c
c  find molecular units:
c
      call fndmol(a(i2))
c
c  classify all the nbos according to donor/acceptor type:
c
      call nbocla(a(i1),accthr)
c
c  output transformation matrices for the pnho and nho basis sets,
c  and the nho density matrix, nho fock matrix, and nho dipole matrices:
c
c  the section of the code makes use of t and dm.  these matrices
c  will be restored later:  [note: do not destroy info already stored
c  in a(i0) and a(i1)]
c
c  reorganize the scratch vector:
c
      i0   = 1
      i1   = i0 + natoms*natoms
      i2   = i1 + ndim
      i3   = i2 + ndim*ndim
      i4   = i3 + ndim*ndim
c     iend = i4 + ndim*(ndim+5)
c
c  output the ao-pnho transformation and the pnho overlap matrix:
c
      io = ioinqr(jprint(20))
      jo = ioinqr(jprint(30))
      if((io.eq.iprnt.or.io.eq.iwrit).or.
     +   (jo.eq.iprnt.or.jo.eq.iwrit)) then
        call fepnao(t)
        call fetnho(a(i2))
        call matmlt(t,a(i2),a(i3),ndim,nbas)
        call fesraw(a(i2))
        call normlz(t,a(i2),ndim,nbas)
        if(jo.eq.iprnt.or.jo.eq.iwrit) then
          title = 'pnhos in the ao basis:'
          call aout(t,ndim,nbas,nbas,title,1,jprint(30))
        end if
        if(io.eq.iprnt.or.io.eq.iwrit) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'pnho overlap matrix:'
          call aout(a(i2),ndim,-nbas,nbas,title,3,jprint(20))
        end if
      endif
c
c  form the ao-nho transformation matrix:
c
      call fetnao(t)
      call fetnho(a(i2))
      call matmlt(t,a(i2),a(i3),ndim,nbas)
c
c  output the ao-nho transformation matrix:
c
      io = ioinqr(jprint(28))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        title = 'nhos in the ao basis:'
        call aout(t,ndim,nbas,nbas,title,1,jprint(28))
      end if
c
c  output the nao-nho transformation matrix:
c
      io = ioinqr(jprint(33))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fetnho(a(i2))
        title = 'nhos in the nao basis:'
        call aout(a(i2),ndim,nbas,nbas,title,2,jprint(33))
      end if
c
c  output the nho-mo transformation matrix:
c
      io = ioinqr(jprint(38))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call frmtmo(t,a(i2),a(i3),a(i4),3,jprint(38))
      end if
c
c  output the nho density matrix:
c
      io = ioinqr(jprint(34))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fedraw(dm,a(i2))
        if(iwdm.eq.1) then
          call fesraw(a(i2))
          call simtrs(dm,a(i2),a(i3),ndim,nbas)
        end if
        call simtrs(dm,t,a(i2),ndim,nbas)
        title = 'nho density matrix:'
        call aout(dm,ndim,-nbas,nbas,title,3,jprint(34))
      end if
c
c  output the nho fock matrix:
c
      io = ioinqr(jprint(29))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fefao(a(i2),iwfock)
        if(iwfock.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nho fock matrix:'
          call aout(a(i2),ndim,-nbas,nbas,title,3,jprint(29))
        end if
      end if
c
c  output the nho dipole matrices:
c
      io = ioinqr(jprint(52))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        ix = 1
        call fedxyz(a(i2),ix)
        if(ix.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nho x dipole integrals:'
          call aout(a(i2),ndim,-nbas,nbas,title,3,jprint(52))
        end if
        ix = 2
        call fedxyz(a(i2),ix)
        if(ix.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nho y dipole integrals:'
          call aout(a(i2),ndim,-nbas,nbas,title,3,jprint(52))
        end if
        ix = 3
        call fedxyz(a(i2),ix)
        if(ix.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nho z dipole integrals:'
          call aout(a(i2),ndim,-nbas,nbas,title,3,jprint(52))
        end if
      end if
c
c  output transformation matrices for the pnbo and nbo basis sets,
c  and the nbo density matrix, nbo fock matrix, and nbo dipole matrices:
c
c  [note: do not destroy info already stored in a(i0) and a(i1)]
c
c  reorganize the scratch vector:
c
      i0   = 1
      i1   = i0 + natoms*natoms
      i2   = i1 + ndim
      i3   = i2 + ndim*ndim
      i4   = i3 + ndim*ndim
c     iend = i4 + ndim*(ndim+5)
c
c  output the ao-pnbo transformation and the pnbo overlap matrix:
c
      io = ioinqr(jprint(21))
      jo = ioinqr(jprint(25))
      if((io.eq.iprnt.or.io.eq.iwrit).or.
     +   (jo.eq.iprnt.or.jo.eq.iwrit)) then
        call fepnao(t)
        call fetnab(a(i2))
        call matmlt(t,a(i2),a(i3),ndim,nbas)
        call fesraw(a(i2))
        call normlz(t,a(i2),ndim,nbas)
        if(jo.eq.iprnt.or.jo.eq.iwrit) then
          title = 'pnbos in the ao basis:'
          call aout(t,ndim,nbas,nbas,title,1,jprint(25))
        end if
        if(io.eq.iprnt.or.io.eq.iwrit) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'pnbo overlap matrix:'
          call aout(a(i2),ndim,-nbas,nbas,title,4,jprint(21))
        end if
      end if
c
c  form the ao-nbo transformation matrix:
c
      call fetnao(t)
      call fetnab(a(i2))
      call matmlt(t,a(i2),a(i3),ndim,nbas)
c
c  save the ao-nbo transformation, nbo occs, and nbo labels on nbo daf:
c
      call svnbo(t,a(i1),a(i2))
c
c  write the ao-nbo transformation with nbo labels and occupancies:
c
      if(ioinqr(iwtnbo).eq.iwrit) call wrtnbo(t,a(i1),iwtnbo)
c
c  print the ao-nbo transformation matrix:
c
      if(ioinqr(iwtnbo).eq.iprnt) then
        title = 'nbos in the ao basis:'
        call aout(t,ndim,nbas,nbas,title,1,iwtnbo)
      end if
c
c  write the nao-nbo transformation matrix:
c
      if(ioinqr(iwtnab).eq.iwrit) then
        call fetnab(a(i2))
        call wrtnab(a(i2),iwtnab)
      end if
c
c  print the nao-nbo transformation to the output file:
c
      if(ioinqr(iwtnab).eq.iprnt) then
        call fetnab(a(i2))
        title = 'nbos in the nao basis:'
        call aout(a(i2),ndim,nbas,nbas,title,2,iwtnab)
      end if
c
c  output the nho-nbo transformation matrix:
c
      io = ioinqr(jprint(41))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call fetnho(a(i2))
        call NBOtransp(a(i2),ndim,nbas)
        call fetnab(a(i3))
        call matmlt(a(i2),a(i3),a(i4),ndim,nbas)
        title = 'nbos in the nho basis:'
        call aout(a(i2),ndim,nbas,nbas,title,3,jprint(41))
      end if
c
c  output the nbo-mo transformation matrix:
c
      io = ioinqr(jprint(45))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        call frmtmo(t,a(i2),a(i3),a(i4),4,jprint(45))
      end if
c
c  form the nbo density matrix:
c
      call fedraw(dm,a(i2))
      if(iwdm.eq.1.and..not.ortho) then
        call fesraw(a(i2))
        call simtrs(dm,a(i2),a(i3),ndim,nbas)
      end if
      call simtrs(dm,t,a(i2),ndim,nbas)
c
c  output the nbo density matrix:
c
      io = ioinqr(jprint(16))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        title = 'nbo density matrix:'
        call aout(dm,ndim,-nbas,nbas,title,4,jprint(16))
      end if
c
c  output the nbo fock matrix:
c
      call fefao(a(i2),iwfock)
      if(iwfock.ne.0) then
        call simtrs(a(i2),t,a(i3),ndim,nbas)
        call svfnbo(a(i2))
        io = ioinqr(jprint(37))
        if(io.eq.iprnt.or.io.eq.iwrit) then
          title = 'nbo fock matrix:'
          call aout(a(i2),ndim,-nbas,nbas,title,4,jprint(37))
        end if
      end if
c
c  output the nbo dipole matrices:
c
      io = ioinqr(jprint(53))
      if(io.eq.iprnt.or.io.eq.iwrit) then
        ix = 1
        call fedxyz(a(i2),ix)
        if(ix.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nbo x dipole integrals:'
          call aout(a(i2),ndim,-nbas,nbas,title,4,jprint(53))
        end if
        ix = 2
        call fedxyz(a(i2),ix)
        if(ix.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nbo y dipole integrals:'
          call aout(a(i2),ndim,-nbas,nbas,title,4,jprint(53))
        end if
        ix = 3
        call fedxyz(a(i2),ix)
        if(ix.ne.0) then
          call simtrs(a(i2),t,a(i3),ndim,nbas)
          title = 'nbo z dipole integrals:'
          call aout(a(i2),ndim,-nbas,nbas,title,4,jprint(53))
        end if
      end if
c
c  perform perturbative analysis of the nbo fock matrix:
c
      if(jprint(3).eq.1.and.iwfock.ne.0) call fnboan(a(i1),a(i2),a(i3))
c
c  print the nbo summary:
c
      if(jprint(6).eq.1) then
        i0   = 1
        i1   = i0 + natoms*natoms
        i2   = i1 + ndim
        i3   = i2 + ndim*ndim
        i4   = i3 + ndim
        i5   = i4 + natoms
c       iend = i5 + ndim
        call nbosum(a(i2),a(i1),a(i3),a(i4),a(i5))
      end if
c
c  continue with the construction of the nlmos:
c
      if(jprint(8).ne.0) then
c
c  store in a(i0) the vector reson(ndim), the squares of the diagonal
c  elements of the nbo to nlmo transformation matrix.  ialarm sounds
c  the alarm that the nlmo step is to be skipped:
c
c   dm   : nbo density         ! transformed to nlmo basis on return
c   a(i0): reson(ndim)         ! percentages of parent nbo
c   a(i1): lmoocc(ndim)        ! nlmo occupancies
c   a(i2): tnlmo(ndim,ndim)    ! nbo-nlmo transform
c   a(i3): tsym                ! scratch
c
c  (do not destroy the wiberg bond index!)
c
        i0   = 1 + natoms*natoms
        i1   = i0 + ndim
        i2   = i1 + ndim
        i3   = i2 + ndim*ndim
c       iend = i3 + ndim*ndim
        call nlmo(nbas,dm,a(i1),a(i2),a(i3),a(i0),nocc,ialarm)
        if(ialarm.ne.0) return
c
c  save the nbo to nlmo transformation matrix on the nbo daf:
c
        call svtlmo(a(i2))
c
c  form the nao to nlmo transformation in t:
c
        call fetnab(t)
        call matmlt(t,a(i2),a(i3),ndim,nbas)
c
c  set up storage for lmoanl:
c
c   a(i0): reson(ndim)
c   a(i1): lmoocc(ndim)
c   a(i2): ts(ndim)
c   a(i3): border(natoms,natoms)
c   a(i4): owbord(natoms,natoms)
c   a(i5): atlmo(nocc,natoms)
c   a(i6): siab(nocc,nab)
c
c  (do not destroy the wiberg bond index!)
c
        nab = natoms*(natoms-1)/2
        if(natoms.eq.1) nab = 1
        i0   = 1 + natoms*natoms
        i1   = i0 + ndim
        i2   = i1 + ndim
        i3   = i2 + ndim
        i4   = i3 + natoms*natoms
        i5   = i4 + natoms*natoms
        i6   = i5 + nocc*natoms
        i7   = i6 + nocc*nab
c       iend = i7 + ndim*ndim
        call copy(dm,a(i7),ndim,nbas,nbas)
        call lmoanl(t,a(i7),a(i0),a(i1),a(i2),a(i3),a(i4),a(i5),
     +              a(i6),nocc,nab)
c
c  output transformation matrices for the pnlmo and nlmo basis sets,
c  and the nlmo density matrix, nlmo fock matrix, and nlmo dipole matrices:
c
c  reorganize the scratch vector:
c
c  (do not destroy the wiberg bond index!)
c
        i0   = 1 + natoms*natoms
        i1   = i0 + ndim*ndim
        i2   = i1 + ndim*ndim
c       iend = i2 + ndim*(ndim+5)
c
c  output the ao-pnlmo transformation and the pnlmo overlap matrix:
c
        io = ioinqr(jprint(48))
        jo = ioinqr(jprint(49))
        if((io.eq.iprnt.or.io.eq.iwrit).or.
     +     (jo.eq.iprnt.or.jo.eq.iwrit)) then
          call fepnao(t)
          call fetnab(a(i0))
          call matmlt(t,a(i0),a(i1),ndim,nbas)
          call fetlmo(a(i0))
          call matmlt(t,a(i0),a(i1),ndim,nbas)
          call fesraw(a(i0))
          call normlz(t,a(i0),ndim,nbas)
          if(jo.eq.iprnt.or.jo.eq.iwrit) then
            title = 'pnlmos in the ao basis:'
            call aout(t,ndim,nbas,nbas,title,1,jprint(49))
          end if
          if(io.eq.iprnt.or.io.eq.iwrit) then
            call simtrs(a(i0),t,a(i1),ndim,nbas)
            title = 'pnlmo overlap matrix:'
            call aout(a(i0),ndim,-nbas,nbas,title,5,jprint(48))
          end if
        end if
c
c  form the ao-nlmo transformation matrix:
c
        call fetnao(t)
        call fetnab(a(i0))
        call matmlt(t,a(i0),a(i1),ndim,nbas)
        call fetlmo(a(i0))
        call matmlt(t,a(i0),a(i1),ndim,nbas)
c
c  save the ao-nlmo transformation on nbo daf:
c
        call svnlmo(t)
c
c  write out the ao-nlmo transformation matrix:
c
        io = ioinqr(jprint(23))
        if(io.eq.iwrit) call wrnlmo(t,dm,jprint(23))
c
c  print the ao-nlmo transformation matrix:
c
        if(io.eq.iprnt) then
          title = 'nlmos in the ao basis:'
          call aout(t,ndim,nbas,nbas,title,1,jprint(23))
        end if
c
c  output the nao-nlmo transformation matrix:
c
        io = ioinqr(jprint(18))
        if(io.eq.iprnt.or.io.eq.iwrit) then
          call fetnab(a(i0))
          call fetlmo(a(i1))
          call matmlt(a(i0),a(i1),a(i2),ndim,nbas)
          title = 'nlmos in the nao basis:'
          call aout(t,ndim,nbas,nbas,title,2,jprint(18))
        end if
c
c  output the nho-nlmo transformation matrix:
c
        io = ioinqr(jprint(24))
        if(io.eq.iprnt.or.io.eq.iwrit) then
          call fetnho(a(i0))
          call NBOtransp(a(i0),ndim,nbas)
          call fetnab(a(i1))
          call matmlt(a(i0),a(i1),a(i2),ndim,nbas)
          call fetlmo(a(i1))
          call matmlt(a(i0),a(i1),a(i2),ndim,nbas)
          title = 'nlmos in the nho basis:'
          call aout(a(i0),ndim,nbas,nbas,title,3,jprint(24))
        end if
c
c  output the nbo-nlmo transformation matrix:
c
        io = ioinqr(jprint(47))
        if(io.eq.iprnt.or.io.eq.iwrit) then
          call fetlmo(a(i0))
          title = 'nlmos in the nbo basis:'
          call aout(a(i0),ndim,nbas,nbas,title,4,jprint(47))
        end if
c
c  output the nlmo-mo transformation matrix:
c
        io = ioinqr(jprint(13))
        if(io.eq.iprnt.or.io.eq.iwrit) then
          call frmtmo(t,a(i0),a(i1),a(i2),5,jprint(13))
        end if
c
c  output the nlmo density matrix:
c
        io = ioinqr(jprint(17))
        if(io.eq.iprnt.or.io.eq.iwrit) then
          title = 'nlmo density matrix:'
          call aout(dm,ndim,-nbas,nbas,title,5,jprint(17))
        end if
c
c  output the nlmo fock matrix:
c
        io = ioinqr(jprint(15))
        if(io.eq.iprnt.or.io.eq.iwrit) then
          call fefao(a(i0),iwfock)
          if(iwfock.ne.0) then
            call simtrs(a(i0),t,a(i1),ndim,nbas)
            title = 'nlmo fock matrix:'
            call aout(a(i0),ndim,-nbas,nbas,title,5,jprint(15))
          end if
        end if
c
c  output the nlmo dipole matrices:
c
        io = ioinqr(jprint(54))
        if(io.eq.iprnt.or.io.eq.iwrit) then
          ix = 1
          call fedxyz(a(i0),ix)
          if(ix.ne.0) then
            call simtrs(a(i0),t,a(i1),ndim,nbas)
            title = 'nlmo x dipole integrals:'
            call aout(a(i0),ndim,-nbas,nbas,title,5,jprint(54))
          end if
          ix = 2
          call fedxyz(a(i0),ix)
          if(ix.ne.0) then
            call simtrs(a(i0),t,a(i1),ndim,nbas)
            title = 'nlmo y dipole integrals:'
            call aout(a(i0),ndim,-nbas,nbas,title,5,jprint(54))
          end if
          ix = 3
          call fedxyz(a(i0),ix)
          if(ix.ne.0) then
            call simtrs(a(i0),t,a(i1),ndim,nbas)
            title = 'nlmo z dipole integrals:'
            call aout(a(i0),ndim,-nbas,nbas,title,5,jprint(54))
          end if
        end if
c
c  perform the nbo/nlmo dipole moment analysis:
c
c  dm   :  nlmo density matrix
c  t    :  ao-nlmo transformation matrix
c  a(i1):  c(ndim,ndim)
c  a(i2):  tnbo(ndim,ndim)
c  a(i3):  dx(ndim,ndim)
c  a(i4):  dy(ndim,ndim)
c  a(i5):  dz(ndim,ndim)
c  a(i6):  scr(ndim,ndim)
c  a(i7):  index(ndim)
c
c  (do not destroy the wiberg bond index!)
c
        if(jprint(46).ne.0) then
          i1   = 1 + natoms*natoms
          i2   = i1 + ndim*ndim
          i3   = i2 + ndim*ndim
          i4   = i3 + ndim*ndim
          i5   = i4 + ndim*ndim
          i6   = i5 + ndim*ndim
          i7   = i6 + ndim*ndim
c         iend = i7 + ndim
          call dipanl(dm,t,a(i1),a(i2),a(i3),a(i4),a(i5),a(i6),a(i7))
        end if
      end if
c
c  perform natural resonance theory analysis:
c
      if(jprint(32).ne.0) then
c
c  carefully determine the maximum number of resonance structures
c  (maxres) that the scratch vector can accomodate.  assume that
c  there will be roughly 6(=nel)  elements required per atom to store
c  the topo matrices for each resonance structure: (1 for number of
c  bonds, 1 for number of lone pairs, and 4 bonded atoms -- see
c  sr topstr)
c
        nel = 6
        tot = zero
        do 80 ibas = 1,nbas
          tot = tot + dm(ibas,ibas)
   80   continue
        nelec = nint(tot)
        nlow = natoms*(natoms-1)/2
        maxref = max(jprint(56),1)
c
c  carefully determine the maximum number of resonance structures (maxres)
c  which the scratch vector can accomodate.  assume ndim is larger than
c  maxres (this is not usually the case):
c
        ic = ndim*ndim + 4*ndim + mxao*ndim + ndim + mxbo*mxbo +
     +       mxbo*mxbo + mxbo + mxbo + mxao*mxao + mxao*mxao +
     +       mxao + mxao + mxao + natoms*natoms + ndim*maxref +
     +       ndim*ndim + maxref + maxref + ndim*maxref + ndim +
     +       ndim*ndim + ndim*ndim + ndim*ndim + natoms*natoms +
     +       maxref - memory
        ib = ndim*maxref + 6*maxref + nlow*maxref + 9 + natoms*nel
        ia = 0
        maxres = int(-ic / ib)
c
c  check this assumption:
c
        if(maxres.gt.ndim) then
          ic = ic - ndim*ndim - ndim*ndim
          ia = 2
          det = sqrt(real(ib * ib - 4 * ia * ic))
          maxres = int((-real(ib) + det) / real(2 * ia))
        end if
        if(maxres.gt.ndim*ndim) then
          ic = ic - ndim*ndim
          ib = ib + 1
          ia = 2
          det = sqrt(real(ib * ib - 4 * ia * ic))
          maxres = int((-real(ib) + det) / real(2 * ia))
        end if
        len = nel * natoms * maxres
c
c  partition the scratch vector:
c
        i0  = 1
        i1  = i0 + natoms*natoms
        i2  = i1 + maxres*maxref
        i3  = i2 + maxres*maxref
        i4  = i3 + maxref
        mem = memory - i4 + 1
c       call nrtdrv(dm,t,a(i0),a(i1),a(i2),a(i3),a(i4),maxres,maxref,
c    +              nlow,len,nelec,mem)
      end if
      return
c
 1390 format(/1x,'WARNING:  the $choose keylist is incompatible with ',
     + 'the nrt analysis for open',/1x,'          shell nbo analyses.',
     + '  program execution will continue, ignoring the',/1x,'       ',
     + '   $choose keylist.')
 1400 format(//1x,'NATURAL BOND ORBITAL ANALYSIS:')
 1410 format(//1x,'NATURAL BOND ORBITAL ANALYSIS,',
     * ' Alpha spin orbitals:')
 1420 format(//1x,'NATURAL BOND ORBITAL ANALYSIS,',
     * ' Beta spin orbitals:')
 2000 format(//1x,'nbo analysis skipped by request.')
      end
c*****************************************************************************
c
c  routines called by the nao drivers:
c
c      subroutine simtrm(a,s,v,ndim,n,iwmulp,iwcubf)
c      subroutine mulana(bs,vmayer,bmayer,iwmulp,iwcubf)
c      subroutine dfgorb(renorm,dm,t,itran,iwcubf,itopt,lfnpr)
c      subroutine nao(t,s,occ,blk,sblk,eval,c,evect,eval2,listao,nblock)
c      subroutine naoanl(dm,spnao,bindex,bindt,bmo,ovpop,f,enao)
c      subroutine frmtmo(t,tmo,c,scr,index,iflg)
c
c*****************************************************************************
      subroutine simtrm(a,s,v,ndim,n,iwmulp,iwcubf)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  similarity transform a ==> s(transpose)*a*s, using scratch vector v.
c   write the diagonal elements of a*s by calling subroutine mulana if
c          iwmulp.ne.0
c     (these are the mulliken populations if s= overlap matrix
c                                       and a= bond-order matrix)
c
      dimension a(ndim,ndim),s(ndim,ndim),v(1)
      call matmlt(a,s,v,ndim,n)
      i1=ndim+1
      if(iwmulp.ne.0) call mulana(a,v(1),v(i1),iwmulp,iwcubf)
      call matml2(s,a,v,ndim,n)
      return
      end
c*****************************************************************************
      subroutine mulana(bs,vmayer,bmayer,iwmulp,iwcubf)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c
c  perform mayer-mulliken bond order analysis
c
c  print out diagonal elements of bs=b*s, where
c      b= bond-order matrix,   s= overlap matrix,   both in original ao basis
c   this constitutes a mulliken population analysis.
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      dimension bs(ndim,ndim),vmayer(natoms),bmayer(natoms,natoms),
     *          iang(5),angl(60),lang(60),cubicf(7)
      character*80 title
      data iang/1hs,1hp,1hd,1hf,1hg/
      data lang/ 51,151,152,153,251,252,253,254,255,
     *          351,352,353,354,355,356,357,
     *          451,452,453,454,455,456,457,458,459,
     *            1,101,102,103,201,202,203,204,205,206,
     *          301,302,303,304,305,306,307,308,309,310,
     *          401,402,403,404,405,406,407,408,409,410,
     *          411,412,413,414,415/
      data angl/4h    ,4hx   ,4hy   ,4hz   ,4hxy  ,4hxz  ,4hyz  ,
     *   4hx2y2,4hz2  ,4h(0) ,4h(c1),4h(s1),4h(c2),4h(s2),4h(c3),
     *   4h(s3),4h(0) ,4h(c1),4h(s1),4h(c2),4h(s2),4h(c3),4h(s3),
     *   4h(c4),4h(s4),
     *          4h    ,4hx   ,4hy   ,4hz   ,4hxx  ,4hxy  ,4hxz  ,
     *   4hyy  ,4hyz  ,4hzz  ,4hxxx ,4hxxy ,4hxxz ,4hxyy ,4hxyz ,
     *   4hxzz ,4hyyy ,4hyyz ,4hyzz ,4hzzz ,4hxxxx,4hxxxy,4hxxxz,
     *   4hxxyy,4hxxyz,4hxxzz,4hxyyy,4hxyyz,4hxyzz,4hxzzz,4hyyyy,
     *   4hyyyz,4hyyzz,4hyzzz,4hzzzz/
      data cubicf/4h(d1),4h(d2),4h(d3),4h(b) ,4h(e1),4h(e2),4h(e3)/
      data zero/0.0d0/
      if(iwcubf.eq.0) go to 20
c  if the f functions are a cubic set, insert the proper labels:
        do 10 i=1,7
          ii=i+9
   10     angl(ii)=cubicf(i)
   20 continue
      if(iwmulp.eq.1) write(lfnpr,1000)
      if(iwmulp.eq.2) write(lfnpr,1100)
      if(iwmulp.eq.2) write(lfnpr,1200)
      sumt=zero
      do 100 i=1,natoms
        vmayer(i)=zero
        do 100 j=1,natoms
  100     bmayer(i,j)=zero
      do 300 iat=1,natoms
        iz=iatno(iat)
        nam=nameat(iz)
        sumat=zero
        do 200 i=1,nbas
          if(lbl(i).ne.iat) go to 200
          lm=lorbc(i)
          l=lm/100
          il=iang(l+1)
          do 130 ilm=1,60
            if(lm.eq.lang(ilm)) go to 140
 130        continue
c
          stop
 140      continue
          occ=bs(i,i)
          sumat=sumat+occ
        if(iwmulp.eq.2) write(lfnpr,1300) i,nam,iat,il,angl(ilm),occ
        do 180 j=1,nbas
          jat=lbl(j)
          if(jat.eq.iat) go to 180
          bmayer(iat,jat)=bmayer(iat,jat)+bs(i,j)*bs(j,i)
  180     continue
  200   continue
        if(iwmulp.eq.1) write(lfnpr,1800) nam,iat,sumat
        if(iwmulp.eq.2) write(lfnpr,1900) nam,iat,sumat
  300   sumt=sumt+sumat
      if(iwmulp.ne.0) write(lfnpr,1600) sumt
      title = 'mayer-mulliken atom-atom bond order matrix:'
      call aout(bmayer,natoms,natoms,natoms,title,0,natoms)
      do 310 i=1,natoms
        do 310 j=1,natoms
  310     vmayer(i)=vmayer(i)+bmayer(i,j)
      title = 'mayer-mulliken valencies by atom:'
      call aout(vmayer,natoms,natoms,1,title,0,1)
      return
 1000 format(//1x,'total gross mulliken populations by atom:',
     * //4x,'atom #',7x,'total')
 1100 format(//1x,'input atomic orbitals, gross mulliken populations:',
     +//1x,' ao',2x,'atom #',2x,'lang',2x,'mulliken population',
     +4x,'atom #',7x,'total')
 1200 format(1x,79('-'))
 1300 format(1x,i3,3x,a2,i3,2x,a1,a4,f13.7)
 1600 format(/1x,'total number of electrons: ',f11.6)
 1800 format(5x,a2,i3,f15.7)
 1900 format(44x,a2,i3,f15.7)
      end
c*****************************************************************************
      subroutine dfgorb(renorm,dm,t,itran,iwcubf,itopt,lfnpr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/list(6,maxbas),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      dimension t(ndim,ndim),dm(ndim,ndim),a(6,6),b(6),m(6),
     *  renorm(ndim),
     *  lf(3,3),lfcub(3,3),lft(3,3),lfcubt(3,3),lg(3,3),lgt(3,3)
      data lf    /301,304,306,302,307,309,303,308,310/
      data lfcub /306,304,301,309,302,307,303,308,310/
      data lft   /151,356,352,152,357,353,153,354,351/
      data lfcubt/151,355,351,152,356,352,153,357,353/
      data lg    /402,407,409,403,408,410,405,412,414/
      data lgt   /251,455,459,252,452,456,253,453,457/
      data zero,one,two,three,four,six,eight
     *    /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,6.0d0,8.0d0/
c**********************************************************************
c
c     subroutine to transform the overlap and density matrices if d, f, g
c  orbitals are present, or transform a transformation matrix so that it
c  starts from the raw ao instead of the pure ao basis
c  this transformation will not work if dm is the bond-order matrix.
c
c         list(6,maxbas): the list of functions to be transformed
c               list(1,i),list(2,i),list(3,i) are corresponding sets of
c               d,f, or g functions.  it is assumed that, for example,
c               the third dx2 function found in the angular momenta list "lorb"
c               corresponds to the third dy2 and the third dz2 functions in
c               the list of basis functions!
c         itran=idtran+iftran+igtran
c         idtran: the number of sets of cartesian d orbitals found
c         iftran: the number of sets of cartesian f orbitals found
c         igtran: the number of sets of cartesian g orbitals found
c         a     : the transformation matrix
c
c         itopt : if zero, transform dm and s (in t) from raw ao to pure
c                                                                   ao basis
c                 if one,  pre-multiply t by the ao to pure ao transf.
c                        --- this converts a transf. that starts from pure aos
c                            to a transf. that starts from the raw aos
c
c         renorm: renormalization vector for cartesian to pure transform.
c                 (produced if itopt=0, used as input if itopt=1)
c
c**********************************************************************
      do 10 i=1,nbas
  10  lorb(i)=0
      idtran=0
      n1=0
      n2=0
      n3=0
      n4=0
      n5=0
      n6=0
c...construct list:
      do 70 ibas=1,nbas
c   dx2:
        if(lorbc(ibas).ne.201) go to 20
          n1=n1+1
          list(1,n1)=ibas
          go to 70
c   dy2:
   20   if(lorbc(ibas).ne.204) go to 30
          n2=n2+1
          list(2,n2)=ibas
          go to 70
c   dz2:
   30   if(lorbc(ibas).ne.206) go to 40
          n3=n3+1
          list(3,n3)=ibas
          go to 70
c   label dxy:
   40   if(lorbc(ibas).ne.202) go to 50
          n4=n4+1
          lorb(ibas)=251
          go to 70
c   label dxz:
   50   if(lorbc(ibas).ne.203) go to 60
          n5 =n5+1
          lorb(ibas)=252
          go to 70
c   label dyz:
   60   if(lorbc(ibas).ne.205) go to 70
          n6=n6+1
          lorb(ibas)=253
   70   continue
      if(n1.ne.n2.or.n1.ne.n3) go to 1950
      if(n1.ne.n4.or.n1.ne.n5.or.n1.ne.n6) go to 1950
      idtran=n1
      if(idtran.eq.0) go to 160
c set up transform. coeff:
c  s=r2=x2+y2+z2:
      a(1,1)= one
      a(2,1)= one
      a(3,1)= one
c  d(x2-y2):
      a(1,2)= one
      a(2,2)=-one
      a(3,2)= zero
c  d(3z2-r2)=-x2-y2+2z2:
      a(1,3)=-one
      a(2,3)=-one
      a(3,3)= two
      if(itopt.eq.0) go to 110
      do 90 j=1,3
        renor=renorm(list(j,1))
        do 90 i=1,3
   90     a(i,j)=a(i,j)*renor
      call NBOtransp(a,6,3)
  110 continue
c...loop over d sets in dlist:
      do 150 id=1,idtran
        m(1)=list(1,id)
        m(2)=list(2,id)
        m(3)=list(3,id)
c...transform s and dm:
        if(itopt.ne.0) call limtrn(t,m,a,b,ndim,nbas,6,3,-1)
        if(itopt.ne.0) go to 150
          call limtrn(t,m,a,b,ndim,nbas,6,3,0)
          call limtrn(dm,m,a,b,ndim,nbas,6,3,0)
c...set the orbital labels for the 3 orbitals transformed:
        lorb(m(1))=51
        lorb(m(2))=254
        lorb(m(3))=255
  150   continue
c**********************************************************************
  160 continue
c  f orbitals
      iftran=0
      do 400 ifblk=1,3
        n1=0
        n2=0
        n3=0
        if(iwcubf.ne.0) go to 190
          lf1=lf(1,ifblk)
          lf2=lf(2,ifblk)
          lf3=lf(3,ifblk)
          go to 200
  190   continue
          lf1=lfcub(1,ifblk)
          lf2=lfcub(2,ifblk)
          lf3=lfcub(3,ifblk)
  200   continue
c...construct the list:
        do 260 ibas=1,nbas
          if(lorbc(ibas).ne.lf1) go to 220
            n1=n1+1
            list(1,n1)=ibas
            go to 260
  220     if(lorbc(ibas).ne.lf2) go to 230
            n2=n2+1
            list(2,n2)=ibas
            go to 260
  230     if(lorbc(ibas).ne.lf3) go to 260
            n3=n3+1
            list(3,n3)=ibas
            go to 260
  260     continue
        if(n1.ne.n2.or.n1.ne.n3) go to 1960
        if(ifblk.eq.1) iftran=n1
        if((ifblk.ne.1).and.(iftran.ne.n1)) go to 1960
        if(iftran.eq.0) go to 500
        if(iwcubf.eq.0) go to 270
c set up transform. coeff, cubic f orbitals
c  px=x*r2, py=y*r2, pz=z*z2
          a(1,1)= one
          a(2,1)= one
          a(3,1)= one
c  fx(z2-y2), fy(z2-x2), fz(x2-y2)
          a(1,2)= one
          a(2,2)=-one
          a(3,2)= zero
c  fx(5z2-3r2), fy(5y2-3r2), fz(5z2-3r2)
          a(1,3)=-three
          a(2,3)=-three
          a(3,3)= two
          go to 310
  270   if(ifblk.gt.1) go to 280
c set up transform. coeff, for first f block
c  px=x*r2
          a(1,1)= one
          a(2,1)= one
          a(3,1)= one
c  fx(x2-3y2)
          a(1,2)= one
          a(2,2)=-three
          a(3,2)= zero
c  fx(5z2-r2)
          a(1,3)=-one
          a(2,3)=-one
          a(3,3)= four
          go to 310
  280   if(ifblk.eq.3) go to 290
c set up transform. coeff, for second f block
c  py=y*r2
          a(1,1)= one
          a(2,1)= one
          a(3,1)= one
c  fy(3x2-y2)
          a(1,2)= three
          a(2,2)=-one
          a(3,2)= zero
c  fy(5z2-r2)
          a(1,3)=-one
          a(2,3)=-one
          a(3,3)= four
          go to 310
  290   continue
c set up transform. coeff, for third f block
c  pz z*r2
          a(1,1)= one
          a(2,1)= one
          a(3,1)= one
c  fz(x2-y2)
          a(1,2)= one
          a(2,2)=-one
          a(3,2)= zero
c  fz(5z2-3r2)
          a(1,3)=-three
          a(2,3)=-three
          a(3,3)= two
  310   continue
      if(itopt.eq.0) go to 330
      do 320 j=1,3
        renor=renorm(list(j,1))
        do 320 i=1,3
  320     a(i,j)=a(i,j)*renor
      call NBOtransp(a,6,3)
  330 continue
c...loop over f sets in list:
        do 390 it=1,iftran
          m(1)=list(1,it)
          m(2)=list(2,it)
          m(3)=list(3,it)
c...transform s and dm, or t (if itopt.ne.0)
        if(itopt.ne.0) call limtrn(t,m,a,b,ndim,nbas,6,3,-1)
        if(itopt.ne.0) go to 340
          call limtrn(t,m,a,b,ndim,nbas,6,3,0)
          call limtrn(dm,m,a,b,ndim,nbas,6,3,0)
c...fix the orbital labels for the 3 orbitals transformed:
  340   continue
          if(iwcubf.ne.0) go to 350
            lorb(m(1))=lft(1,ifblk)
            lorb(m(2))=lft(2,ifblk)
            lorb(m(3))=lft(3,ifblk)
            go to 390
  350     continue
            lorb(m(1))=lfcubt(1,ifblk)
            lorb(m(2))=lfcubt(2,ifblk)
            lorb(m(3))=lfcubt(3,ifblk)
  390     continue
  400   continue
c   search for fxyz and relabel:
      lf1=305
      lf1t=355
      if(iwcubf.ne.0) lf1t=354
      n1=0
      do 420 ibas=1,nbas
        if(lorbc(ibas).ne.lf1) go to 420
          n1=n1+1
          lorb(ibas)=lf1t
  420     continue
      if(iftran.ne.n1) go to 1960
  500 continue
c  g orbitals
      igtran=0
      do 800 igblk=1,3
        n1=0
        n2=0
        n3=0
          lg1=lg(1,igblk)
          lg2=lg(2,igblk)
          lg3=lg(3,igblk)
c...construct the list:
        do 560 ibas=1,nbas
          lang=lorbc(ibas)
          if(lang.ne.lg1) go to 520
            n1=n1+1
            list(1,n1)=ibas
            go to 560
  520     if(lang.ne.lg2) go to 530
            n2=n2+1
            list(2,n2)=ibas
            go to 560
  530     if(lang.ne.lg3) go to 560
            n3=n3+1
            list(3,n3)=ibas
            go to 560
  560     continue
        if(n1.ne.n2.or.n1.ne.n3) go to 1970
        if(igblk.eq.1) igtran=n1
        if((igblk.ne.1).and.(igtran.ne.n1)) go to 1970
        if(igtran.eq.0) go to 1000
          if(igblk.gt.1) go to 580
c set up transform. coeff, for first g block
c  dxy=xy*r2
            a(1,1)= one
            a(2,1)= one
            a(3,1)= one
c  g(2s)
            a(1,2)= one
            a(2,2)=-one
            a(3,2)= six
c  g(4s)
            a(1,3)= one
            a(2,3)=-one
            a(3,3)= zero
            go to 610
  580     if(igblk.eq.3) go to 590
c set up transform. coeff, for second g block
c  dxz=xz*r2
            a(1,1)= one
            a(2,1)= one
            a(3,1)= one
c  g(1c)
            a(1,2)=-three
            a(2,2)=-three
            a(3,2)= four
c  g(3c)
            a(1,3)= one
            a(2,3)=-three
            a(3,3)= zero
            go to 610
  590     continue
c set up transform. coeff, for third g block
c  dyz=yz*r2
            a(1,1)= one
            a(2,1)= one
            a(3,1)= one
c  g(1s)
            a(1,2)=-three
            a(2,2)=-three
            a(3,2)= four
c  g(3s)
            a(1,3)= three
            a(2,3)=-one
            a(3,3)= zero
  610   continue
      if(itopt.eq.0) go to 630
      do 620 j=1,3
        renor=renorm(list(j,1))
        do 620 i=1,3
  620     a(i,j)=a(i,j)*renor
      call NBOtransp(a,6,3)
  630 continue
c...loop over g sets in list:
        do 690 it=1,igtran
          m(1)=list(1,it)
          m(2)=list(2,it)
          m(3)=list(3,it)
c...transform s and dm, or t (if itopt.ne.0)
          if(itopt.ne.0) call limtrn(t,m,a,b,ndim,nbas,6,3,-1)
          if(itopt.ne.0) go to 660
            call limtrn(t,m,a,b,ndim,nbas,6,3,0)
            call limtrn(dm,m,a,b,ndim,nbas,6,3,0)
c...fix the orbital labels for the 3 orbitals transformed:
  660   continue
          lorb(m(1))=lgt(1,igblk)
          lorb(m(2))=lgt(2,igblk)
          lorb(m(3))=lgt(3,igblk)
  690     continue
  800   continue
c  g orbitals --- fourth (6x6) block
        n1=0
        n2=0
        n3=0
        n4=0
        n5=0
        n6=0
c...construct the list:
        do 870 ibas=1,nbas
          lang=lorbc(ibas)
          if(lang.ne.401) go to 820
            n1=n1+1
            list(1,n1)=ibas
            go to 870
  820     if(lang.ne.411) go to 830
            n2=n2+1
            list(2,n2)=ibas
            go to 870
  830     if(lang.ne.415) go to 840
            n3=n3+1
            list(3,n3)=ibas
            go to 870
  840     if(lang.ne.404) go to 850
            n4=n4+1
            list(1,n4)=ibas
            go to 870
  850     if(lang.ne.406) go to 860
            n5=n5+1
            list(2,n5)=ibas
            go to 870
  860     if(lang.ne.413) go to 870
            n6=n6+1
            list(3,n6)=ibas
            go to 870
  870     continue
        if(igtran.ne.n1.or.n1.ne.n2.or.n1.ne.n3) go to 1970
        if(n1.ne.n4.or.n1.ne.n5.or.n1.ne.n6) go to 1970
c set up transform. coeff, for fourth g block
c  s=(r2)2
            a(1,1)= one
            a(2,1)= one
            a(3,1)= one
            a(4,1)= two
            a(5,1)= two
            a(6,1)= two
c  d(3z2-r2)
            a(1,2)=-one
            a(2,2)=-one
            a(3,2)= two
            a(4,2)=-two
            a(5,2)= one
            a(6,2)= one
c  d(x2-y2)
            a(1,3)= one
            a(2,3)=-one
            a(3,3)= zero
            a(4,3)= zero
            a(5,3)= one
            a(6,3)=-one
c  g(0)
            a(1,4)= three
            a(2,4)= three
            a(3,4)= eight
            a(4,4)= six
            a(5,4)=-six*four
            a(6,4)=-six*four
c  g(2c)
            a(1,5)=-one
            a(2,5)=-one
            a(3,5)= zero
            a(4,5)= six
            a(5,5)=-six
            a(6,5)= zero
c  g(4c)
            a(1,6)= one
            a(2,6)= one
            a(3,6)= zero
            a(4,6)=-six
            a(5,6)= zero
            a(6,6)= zero
      if(itopt.eq.0) go to 930
      do 920 j=1,6
        renor=renorm(list(j,1))
        do 920 i=1,6
  920     a(i,j)=a(i,j)*renor
      call NBOtransp(a,6,6)
  930 continue
        if(itopt.ne.0) call NBOtransp(a,6,6)
c...loop over g sets in list:
        do 960 it=1,igtran
          m(1)=list(1,it)
          m(2)=list(2,it)
          m(3)=list(3,it)
          m(4)=list(4,it)
          m(5)=list(5,it)
          m(6)=list(6,it)
c...transform s and dm:
c...transform s and dm, or t (if itopt.ne.0)
          if(itopt.ne.0) call limtrn(t,m,a,b,ndim,nbas,6,6,-1)
          if(itopt.ne.0) go to 950
            call limtrn(t,m,a,b,ndim,nbas,6,6,0)
            call limtrn(dm,m,a,b,ndim,nbas,6,6,0)
c...change the orbital labels for the 3 orbitals transformed:
  950     continue
          lorb(m(1))=51
          lorb(m(2))=254
          lorb(m(3))=255
          lorb(m(4))=451
          lorb(m(5))=454
          lorb(m(6))=458
  960     continue
c  renormalization, itopt=0 :
 1000 continue
      itran=idtran+iftran+igtran
      if(itopt.ne.0) return
      if(itran.eq.0) go to 1200
      do 1020 i=1,nbas
        x=t(i,i)
 1020   renorm(i)=one/sqrt(x)
      do 1040 i=1,nbas
        do 1040 j=1,nbas
          rij=renorm(i)*renorm(j)
          t(i,j)=t(i,j)*rij
 1040    dm(i,j)=dm(i,j)*rij
c  relabelling of non-transformed orbitals:
 1200 continue
      do 1230 i=1,nbas
        if(lorb(i).ne.0) go to 1230
        lang=lorbc(i)
        lorb(i)=lang
        l=lang/100
        idif=lang-l*100
        if(idif.gt.50) go to 1230
          lorb(i)=lorb(i)+50
 1230   continue
      return
c  error messages:
 1950 write(lfnpr,1951)
 1951 format(' unequal numbers of d function components were',
     +' found in the input.',/,' these cannot be properly transformed-',
     +'-perhaps they were improperly labelled.')
      stop
 1960 write(lfnpr,1961)
 1961 format(' unequal numbers of f function components were',
     +' found in the input.',/,' these cannot be properly transformed-',
     +'-perhaps they were improperly labelled.')
      stop
 1970 write(lfnpr,1971)
 1971 format(' unequal numbers of g function components were',
     +' found in the input.',/,' these cannot be properly transformed-',
     +'-perhaps they were improperly labelled.')
      stop
      end
c*****************************************************************************
      subroutine nao(t,s,occ,blk,sblk,eval,c,evect,eval2,listao,nblock)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c********************************************************************
c
c  main subroutine 'nao' for natural atomic orbital basis set.
c
c
c  input required:
c        s = overlap matrix elements in lower triangle (below diagonal)
c          = density matrix elements in upper triangle (including diag.)
c               (input ao's must(!) be normalized.  on return, s is the
c                full density matrix.  overlap matrix elements are lost.)
c      lbl = list of atomic centers; lbl(i) = n if orbital i is on center n
c     lorb = list of angular momentum type for each orbital;
c            lorb(i) = n if orbital i is of 'type' n.
c            n = ( 51,151,152,153)     = (s,px,py,pz)
c              = (251,252,253,254,255) = (dxy,dxz,dyz,d(x2-y2),d(3z2-r2))
c              = (351-357 for the 7 types of f orbitals)
c              = (451-459 for the 9 types of g orbitals)
c
c  output:
c        t = transformation matrix from input ao's to nao's (rows are
c            labelled by primitive ao's, columns by nao's)
c   naoctr = list of atomic centers for nao's; naoctr(i) = n if nao # i
c            is on center #n.
c     naol = list of angular momentum type for each nao, same format as "lorb"
c
c  before return:
c   lstocc = list of natural minimal basis ('occupied') orbitals;
c            lstocc(i)=n (i=1,...,nocc) means that nao #n belongs
c            to the nmb set.
c   lstemt = list of natural rydberg basis ('empty') orbitals;
c            lstemt(i)=n (i=1,...,nemt) means that nao #n belongs
c            to the nrb set.
c
c  after return:
c   lstocc(i) = 1 only if nao #i belongs to the nmb set.
c
c********************************************************************
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension t(ndim,ndim),s(ndim,ndim),occ(ndim),blk(ndim,ndim),
     +          sblk(mxaolm,mxaolm),eval(nbas),eval2(nbas),
     +          listao(mxaolm,9),c(nblock),evect(mxaolm,mxaolm)
      character*80 title
      data zero,one/0.0d0,1.0d0/
      data iprnt,iwrit,iread/4hprnt,4hwrit,4hread/
c
c  skip t-nao formation if ioinqr(iwpnao).eq.iread:
c
      if(ioinqr(iwpnao).eq.iread) go to 200
c
c  zero transformation matrix t:
c
      do 10 j = 1,nbas
        lstocc(j) = 0
        lstemt(j) = 0
        do 10 i = 1,nbas
   10     t(i,j) = zero
c
c  nf counts the accumulated orbitals:
c
      nf = 0
c
c  nocc counts the accumulated 'occupied' orbitals:
c  nemt counts the accumulated 'empty' orbitals:
c
      nocc = 0
      nemt = 0
c
c  begin main nao loop over atomic centers:
c
      do 140 icntr = 1,natoms
c
c  loop over angular momentum blocks (s,p,d,f,g).  nl counts the number
c  of orbitals in each "m" component of the "l" block:
c
        do 130 il = 1,5
          if(nf.gt.nbas) go to 130
          l = il - 1
          m = 2*l + 1
c
c  scan orbital labels to gather 'listao' of orbitals belonging to
c  proper atom and angular momentum symmetry:
c
          do 20 im = 1,m
            lang = 100*l + im + 50
            nl = 0
            do 20 i = 1,nbas
              if((lbl(i).ne.icntr).or.(lorb(i).ne.lang)) go to 20
              nl = nl + 1
              listao(nl,im) = i
   20       continue
          if(nl.eq.0) go to 140
c
c  load this list of orbitals into blk and sblk (density matrix and
c  overlap elements, resp.), and average the density matrix elements
c  over the m components of l for the atom:
c
          call loadav(listao,nl,m,s,ndim,blk,sblk,mxaolm)
c
c  solve the generalized eigenvalue problem:
c
          call atdiag(nl,blk,sblk,eval,c)
c
c  order the eigenvectors by occupancy eigenvalue:
c
          call rank(eval,nl,nl,larc)
c
c  loop over the 2*l+1 components to store t-nao data:
c
          do 120 im = 1,m
c
c  partition orbitals into 'occupied' and 'empty' sets:
c
            call setbas(lstocc,lstemt,nocc,nemt,icntr,l,nl,nf,ndim)
c
c  store the ordered eigenvectors in t:
c
            do 120 j = 1,nl
              jr = larc(j)
              nf = nf + 1
              occ(nf) = eval(j)
              do 110 i = 1,nl
                iao = listao(i,im)
                ijr = i + nl*(jr-1)
                t(iao,nf) = c(ijr)
  110         continue
c
c  make up nao orbital labels:
c
              naoctr(nf) = icntr
              naol(nf) = l*100 + im + 50
  120       continue
  130     continue
  140   continue
  200 continue
c
c  read in pre-orthogonal t-nao data:
c
      if(ioinqr(iwpnao).ne.iread) go to 300
****        originally this had 2 arguments ??
****        call rdppna(t,occ)
        call rdppna(t,occ,iwpnao)
c
c  recompute and symmetry-average weights, reorganize lstocc if the input
c  pnaos are rpnaos:
c
        if(occ(1).lt.zero) call newwts(s,t,occ)
        nocc = 0
        nemt = 0
        lang = 0
        ilbl = 1
        nlang = 0
        do 280 i = 1,nbas
          if(lstocc(i).gt.0) nocc = nocc + 1
          if((naoctr(i).ne.ilbl).or.(naol(i).ne.lang)) go to 240
            nlang = nlang + 1
            go to 250
  240     if(nlang.gt.mxaolm) mxaolm = nlang
            nlang = 1
            ilbl = naoctr(i)
            lang = naol(i)
  250     continue
          do 260 j = 1,nbas
  260       if(lstocc(j).eq.i) go to 280
          nemt = nemt + 1
          lstemt(nemt) = i
  280   continue
  300 continue
c
c  write preorthogonal t-nao data to lfnppa:
c
      if(ioinqr(iwpnao).eq.iwrit) call wrppna(t,occ,iwpnao)
c
c  save t-pnao for later use in computing the non-orthogonal overlaps
c  between naos or nbos:
c
      call svpnao(t)
      if(ioinqr(iwpnao).eq.iprnt) then
        title = 'pnaos in the pao basis:'
        call aout(t,ndim,nbas,nbas,title,-1,iwpnao)
      end if
c
c  final orthogonalization:
c
      do 450 i = 1,nbas
        do 440 j = 1,i
  440     s(j,i) = s(i,j)
  450   s(i,i) = one
      call worth(s,t,blk,lstocc,ndim,nbas,nocc,occ,eval,blk)
      if(nemt.eq.0) go to 700
      call NBOshmdt(t,s,ndim,nbas,nocc,lstocc,nemt,lstemt,blk)
c
c  put p-pao in upper triangle of s (and diagonal):
c
      call feppao(blk)
      do 460 j = 1,nbas
        do 460 i = 1,j
  460     s(i,j) = blk(i,j)
      call newryd(t,s,blk,c,sblk,evect,occ,eval,eval2,listao,
     *                                             jprint(11))
c
c  select the significant rydbergs, put in "larc".
c  put the list of the rest of the rydbergs into "listao",
c  and set the weightings of these low occupancy orbitals to one.
c  then, do a weighted orthog. among the significant rydbergs,
c  schmidt orthog. the low occ. ryds to these, and finally
c  do a lowdin orthog. among the low occ. ryds.:
c
      call rydsel(lstemt,nemt,nsel1,larc,nsel2,listao,occ)
      if(nsel1.eq.0) go to 690
      call worth(s,t,blk,larc,ndim,nbas,nsel1,occ,eval,blk)
      if(nsel2.eq.0) go to 700
  690 continue
      if(nsel1.ne.0)
     *   call NBOshmdt(t,s,ndim,nbas,nsel1,larc,nsel2,listao,blk)
      call worth(s,t,blk,listao,ndim,nbas,nsel2,occ,eval,blk)
  700 continue
      call feppao(s)
      call simtrs(s,t,occ,ndim,nbas)
      call rediag(s,t,blk,occ,sblk,c,listao,jprint(11))
c
c  return occupied list 'lstocc' of 1's or 0's:
c
      do 820 i = 1,nbas
  820   lstocc(i) = 1
      do 840 i = 1,nemt
  840   lstocc(lstemt(i)) = 0
      return
      end
c*****************************************************************************
      subroutine naoanl(dm,spnao,bindex,bindt,bmo,ovpop,f,enao)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      character*80 title
      logical first,core,allzer
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
c
c    perform the natural population analysis
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbnao/naoc(maxbas),naoa(maxbas),ltyp(maxbas),iprin(maxbas)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nblbl/nlew,nval,lbls(10,maxbas,4)
c
      dimension dm(ndim,ndim),spnao(ndim,ndim),bindex(natoms,natoms),
     *     bindt(natoms),ovpop(natoms,natoms),f(ndim,ndim),enao(ndim),
     *     jprin(maxbas),icore(4),ival(4),nwarn(maxatm),labec(20,2),
     *     occec(20),bmo(natoms,natoms)
      dimension iang(5),angl(25),lang(25),cubicf(7)
      character*8 znameat
c
      data iryd/3hryd/
      data iang/1hs,1hp,1hd,1hf,1hg/
      data lang/ 51,151,152,153,251,252,253,254,255,
     *          351,352,353,354,355,356,357,
     *          451,452,453,454,455,456,457,458,459/
      data angl/4h    ,4hx   ,4hy   ,4hz   ,4hxy  ,4hxz  ,4hyz  ,
     *   4hx2y2,4hz2  ,4h(0) ,4h(c1),4h(s1),4h(c2),4h(s2),4h(c3),
     *   4h(s3),4h(0) ,4h(c1),4h(s1),4h(c2),4h(s2),4h(c3),4h(s3),
     *   4h(c4),4h(s4)/
      data cubicf/4h(d1),4h(d2),4h(d3),4h(b) ,4h(e1),4h(e2),4h(e3)/
      data zero,tenth,two/0.0d0,0.1d0,2.0d0/
c
c  test, test2, allow, and allow2 are numbers used in determining if the
c  density matrix trace is close to being an integer.  test2 (allow2) must
c  be slightly greater than twice test (allow):
c
      data test,test2/1.0d-5,2.1d-5/
      data allow,allow2/1.0d-3,2.1d-3/
      data ichcor,ichval,ichryd/3hcor,3hval,3hryd/
c
c  if the f functions are a cubic set, insert the proper labels:
c
      if(iwcubf.eq.0) goto 20
        do 10 i = 1,7
          ii = i+9
   10     angl(ii) = cubicf(i)
   20 continue
c
c  update the nao atom-atom valency matrix:
c
      do 30 j = 1,natoms
        do 30 i = 1,natoms
          ovpop(i,j) = zero
          bmo(i,j) = zero
   30     bindex(i,j) = zero
      do 50 i = 1,nbas
        iat = naoctr(i)
        do 40 j = 1,nbas
          jat = naoctr(j)
          if(jat.ne.iat) then
            sij = spnao(i,j)
            dmij = dm(i,j)
            dmij2 = dmij*dmij
            dmsij = dmij*sij
            bindex(jat,iat) = bindex(jat,iat) + dmij2
            bmo(jat,iat) = bmo(jat,iat) + dmij
            ovpop(jat,iat) = ovpop(jat,iat) + dmsij
          end if
   40   continue
   50 continue
c
c  determine the nao orbital energies if a fock matrix exists.  use
c  spnao to store tnao:
c
      call fetnao(spnao)
      ifock = iwfock
      if(open.and..not.(alpha.or.beta)) ifock = 0
      if(ifock.eq.1) then
        call fefao(f,iwfock)
        if(iwfock.ne.0) then
          do 80 i = 1,nbas
            enrg = zero
            do 70 j = 1,nbas
              do 60 k = 1,nbas
                enrg = enrg + spnao(j,i)*f(j,k)*spnao(k,i)
   60         continue
   70       continue
            enao(i) = enrg
   80     continue
        end if
      end if
c
c  label nao's as either 'cor', 'val', or 'ryd':
c
      do 200 i = 1,nbas
        ltyp(i) = iryd
  200 continue
      iecp = 0
      do 300 nctr = 1,natoms
        call cortbl(nctr,icore,iecp)
        call valtbl(nctr,ival)
c
c  loop over s,p,d,f orbitals:
c
        do 290 l = 0,3
          ityp = iang(l+1)
          lnum = 2*l + 1
          if(icore(l+1).le.0) goto 240
c
c  label core orbitals:
c
          do 230 m = 1,icore(l+1)
            do 220 la = 1,lnum
              morb = 0
              occ = -1.0d0
              do 210 n = 1,nbas
                lm = naol(n)
                norb = lm/100
                il = iang(norb+1)
                na = mod(naol(n),50)
                if(naoctr(n).eq.nctr.and.il.eq.ityp.and.
     +            dm(n,n).gt.occ.and.ltyp(n).eq.iryd.and.
     +                                         la.eq.na) then
                      morb = n
                      occ = dm(n,n)
                end if
  210         continue
              if(morb.eq.0) then
                write(lfnpr,2500) ityp,nameat(iatno(nctr)),nctr,
     +                            (icore(i),i=1,4),m,la
                stop
              end if
              ltyp(morb) = ichcor
  220       continue
  230     continue
  240     continue
          if(ival(l+1).le.0) goto 280
c
c  label valence orbitals:
c
          do 270 m = 1,ival(l+1)
            do 260 la = 1,lnum
              morb = 0
              occ = -1.0d0
              do 250 n = 1,nbas
                lm = naol(n)
                norb = lm/100
                il = iang(norb+1)
                na = mod(naol(n),50)
                if(naoctr(n).eq.nctr.and.il.eq.ityp.and.
     +            dm(n,n).gt.occ.and.ltyp(n).eq.iryd.and.
     +                                         la.eq.na) then
                      morb = n
                      occ = dm(n,n)
                end if
  250         continue
              if(morb.eq.0) then
                write(lfnpr,2600) ityp,nameat(iatno(nctr)),nctr,
     +                            (ival(i),i=1,4),m,la
                stop
              end if
              ltyp(morb) = ichval
  260       continue
  270     continue
  280     continue
  290   continue
  300 continue
c
c  assign principal quantum numbers using the nao occupancies:
c
      do 390 i = 1,nbas
        iprin(i) = 0
  390 continue
      do 450 nctr = 1,natoms
        iecp = 1
        call cortbl(nctr,ival,iecp)
        iecp = 0
        call cortbl(nctr,icore,iecp)
        do 440 l = 0,4
          ityp = iang(l+1)
          mmax = 2*l + 1
          do 430 m = 1,mmax
            if(l.eq.4) then
              n = 3
            else
              n = ival(l+1) - icore(l+1) + l
            end if
  400       continue
              morb = 0
              occ = -1.0d0
              do 410 j = 1,nbas
                lm = naol(j)
                norb = lm/100
                il = iang(norb+1)
                na = mod(naol(j),50)
                  if(naoctr(j).eq.nctr.and.il.eq.ityp.and.
     +              dm(j,j).gt.occ.and.iprin(j).eq.0.and.
     +                                           m.eq.na) then
                        morb = j
                        occ = dm(j,j)
                  end if
  410           continue
              if(morb.eq.0) goto 420
              n = n + 1
              iprin(morb) = n
            goto 400
  420       continue
  430     continue
  440   continue
  450 continue
c
c  assign principal quantum numbers using the nao fock matrix elements:
c
      if(ifock.eq.0) goto 580
      do 490 i = 1,nbas
        jprin(i) = 0
  490 continue
      do 550 nctr = 1,natoms
        iecp = 1
        call cortbl(nctr,ival,iecp)
        iecp = 0
        call cortbl(nctr,icore,iecp)
        do 540 l = 0,4
          ityp = iang(l+1)
          mmax = 2*l + 1
          do 530 m = 1,mmax
            if(l.eq.4) then
              n = 3
            else
              n = ival(l+1) - icore(l+1) + l
            end if
  500       continue
              morb = 0
              enrg = 1.0d6
              do 510 j = 1,nbas
                lm = naol(j)
                norb = lm/100
                il = iang(norb+1)
                na = mod(naol(j),50)
                  if(naoctr(j).eq.nctr.and.il.eq.ityp.and.
     +              enao(j).lt.enrg.and.jprin(j).eq.0.and.
     +                                           m.eq.na) then
                        morb = j
                        enrg = enao(j)
                  end if
  510           continue
              if(morb.eq.0) goto 520
              n = n + 1
              jprin(morb) = n
            goto 500
  520       continue
  530     continue
  540   continue
  550 continue
  580 continue
c
c  count the total number of electrons:
c
      tot = zero
      do 600 inao = 1,nbas
        tot = tot + dm(inao,inao)
  600 continue
      nel = tot + tenth
c
c  store nel for use by the output routines:
c
      nlew = nel
c
c  check to see if the total number of electrons found is an integer:
c
      if(tot.ge.zero) then
        sumtt = tot + test
        sumti = dint(sumtt)
        sumtf = sumtt - sumti
        if(sumtf.gt.test2) then
          sumtt = tot + allow
          sumti = dint(sumtt)
          sumtf = sumtt - sumti
          if(sumtf.gt.allow2) then
            write(lfnpr,955)
            jprint(4) = -1
          else
            write(lfnpr,956)
          end if
        end if
      else
        write(lfnpr,955)
        jprint(4) = -1
      end if
c
c  write out natural population analysis:
c
      if(jprint(4).ne.0) then
        if(ifock.eq.1) then
          write(lfnpr,900)
        else
          write(lfnpr,910)
        end if
        jctr = 1
        do 700 i = 1,nbas
          nctr = naoctr(i)
          if(nctr.ne.jctr) then
            write(lfnpr,*)
            jctr = nctr
          end if
          iat = iatno(nctr)
          nam = nameat(iat)
          lm = naol(i)
          l = lm/100
          il = iang(l+1)
          do 680 ilm = 1,25
            if(lm.eq.lang(ilm)) goto 690
  680     continue
  690     continue
          occ = dm(i,i)
          if(occ.lt.zero) occ = zero
          if(ifock.eq.1) then
            write(lfnpr,920) i,nam,nctr,il,angl(ilm),ltyp(i),
     +                        jprin(i),il,occ,enao(i)
          else
            write(lfnpr,920) i,nam,nctr,il,angl(ilm),ltyp(i),
     +                        iprin(i),il,occ
          end if
  700   continue
c
c  add note about effective core potentials if used:
c
        iecp = 0
        do 710 i = 1,natoms
          iecp = iecp + iatno(i) - iznuc(i)
  710   continue
        if(ipseud.ne.0) then
          if(alpha.or.beta) iecp = iecp/2
          write(lfnpr,930) iecp
        end if
c
c  write out warnings for low occupancy core orbitals:
c
        crthrs = crtset
        if(alpha.or.beta) crthrs = crthrs - 1.0d0
        do 715 n = 1,natoms
          nwarn(n) = 0
  715   continue
        do 720 i = 1,nbas
          ictr = naoctr(i)
          if(ltyp(i).eq.ichcor.and.dm(i,i).lt.crthrs)
     +       nwarn(ictr) = nwarn(ictr) + 1
  720   continue
        first = .true.
        do 725 n = 1,natoms
          nam = nameat(iatno(n))
          if(nwarn(n).eq.1) then
            if(first) then
              write(lfnpr,931) crthrs,nam,n
              first = .false.
            else
              write(lfnpr,932) crthrs,nam,n
            end if
          else if(nwarn(n).gt.1) then
            if(first) then
              write(lfnpr,933) nwarn(n),crthrs,nam,n
              first = .false.
            else
              write(lfnpr,934) nwarn(n),crthrs,nam,n
            end if
          end if
  725   continue
c
c  write out warnings for population inversions:
c
        if(ifock.eq.1) then
          do 730 n = 1,natoms
            nwarn(n) = 0
  730     continue
          do 735 i = 1,nbas
            ictr = naoctr(i)
            if(iprin(i).ne.jprin(i)) nwarn(ictr) = 1
            iprin(i) = jprin(i)
  735     continue
          first = .true.
          do 738 n = 1,natoms
            nam = nameat(iatno(n))
            if(nwarn(n).gt.0) then
              if(first) then
                write(lfnpr,936) nam,n
                first = .false.
              else
                write(lfnpr,937) nam,n
              end if
            end if
  738     continue
        end if
c
c  summarize the natural population analysis:
c
        write(lfnpr,939)
        sumac = zero
        sumav = zero
        sumar = zero
        nomac = 0
        do 750 i = 1,natoms
          sumc = zero
          sumv = zero
          sumr = zero
c         nam = nameat(iatno(i))
          do 740 j = 1,nbas
            if(naoctr(j).eq.i) then
              occ = dm(j,j)
              if(occ.lt.zero) occ = zero
              if(ltyp(j).eq.ichcor) sumc = sumc + occ
              if(ltyp(j).eq.ichval) sumv = sumv + occ
              if(ltyp(j).eq.ichryd) sumr = sumr + occ
              if(ltyp(j).eq.ichcor) nomac = nomac + 2
            end if
  740     continue
          tot = sumc + sumv + sumr
          if(alpha.or.beta) then
            chg = iznuc(i)/2.0 - tot
          else
            chg = iznuc(i) - tot
          end if
          ecp = dble(iatno(i) - iznuc(i))
          if(alpha.or.beta) ecp = ecp/two
          write(lfnpr,940) znameat(i),i,chg,sumc+ecp,sumv,sumr,tot+ecp
          sumac = sumac + sumc
          sumav = sumav + sumv
          sumar = sumar + sumr
  750   continue
        tot = sumac + sumav + sumar
        chg = -1.0 * tot
        if(alpha.or.beta) then
          nomac = nomac/2
          do 760 i = 1,natoms
            chg = chg + iznuc(i)/2.0d0
  760     continue
        else
          do 770 i = 1,natoms
            chg = chg + iznuc(i)
  770     continue
        end if
        write(lfnpr,950) chg,sumac+dble(iecp),sumav,sumar,
     +                   tot+dble(iecp)
c
c  write out nmb and nrb populations and percentage occupancies:
c
        write(lfnpr,960)
        noma = nel
        nomav = noma - nomac
        suma = sumac + sumav
        if(ipseud.ne.0) then
          ecp = iecp
          suma = suma + ecp
          noma = noma + iecp
          write(lfnpr,970) ecp
        end if
        if(nomac.ne.0) then
          pcent = sumac/dble(nomac) * 100.0d0
          write(lfnpr,980) sumac,pcent,nomac
        else if(sumac.ne.zero) then
          pcent = zero
          write(lfnpr,980) sumac,pcent,nomac
        end if
        if(nomav.ne.0) then
          pcent = sumav/nomav * 100.0d0
          write(lfnpr,990) sumav,pcent,nomav
        else if(sumav.ne.zero) then
          pcent = zero
          write(lfnpr,990) sumav,pcent,nomav
        end if
        if(noma.ne.0) then
          pcent = suma/noma * 100.0d0
        else
          pcent = zero
        end if
        write(lfnpr,1000) suma,pcent,noma
        if(noma.ne.0) then
          pcent = sumar/noma * 100.0d0
          write(lfnpr,1010) sumar,pcent,noma
        else if(sumar.ne.zero) then
          pcent = 0.0d0
          write(lfnpr,1010) sumar,pcent,noma
        end if
c
c  write out natural electron configuration:
c
        write(lfnpr,1040)
        do 899 nctr = 1,natoms
          ict = 0
          iecp = 1
          call cortbl(nctr,icore,iecp)
          do 870 npl = 1,8
            do 860 n = 1,npl
              l = npl - n
              if(l.ge.0.and.l.lt.n) then
                if(n.gt.icore(l+1)+l) then
                  ict = ict + 1
                  labec(ict,1) = n
                  labec(ict,2) = iang(l+1)
                  occec(ict) = zero
                end if
              end if
  860       continue
  870     continue
          do 890 i = 1,nbas
            ictr = naoctr(i)
            if(ictr.eq.nctr.and.ltyp(i).ne.ichcor) then
              norb = naol(i)/100
              il = iang(norb+1)
              do 880 j = 1,ict
                if(iprin(i).eq.labec(j,1).and.
     +                   il.eq.labec(j,2)) then
                  occec(j) = occec(j) + dm(i,i)
                  goto 890
                end if
  880         continue
            end if
  890     continue
          if(labec(1,1).ne.1) then
            core = .true.
          else
            core = .false.
          end if
          thold = 5.0d-3
          jmax = ict
c
c  remove low occupancy subshells:
c
          do 893 jct = 1,ict
  891       continue
            if(occec(jct).lt.thold) then
              allzer = .true.
              do 892 kct = jct,ict-1
                labec(kct,1) = labec(kct+1,1)
                labec(kct,2) = labec(kct+1,2)
                occec(kct)   = occec(kct+1)
                if(occec(kct).ge.thold) allzer = .false.
  892         continue
              occec(ict) = zero
              if(allzer) then
                jmax = jct - 1
                goto 895
              end if
              goto 891
            end if
  893     continue
  895     continue
          nam = nameat(iatno(nctr))
          if(jmax.eq.0) then
            if(.not.core) then
              write(lfnpr,1050) nam,nctr
            else
              write(lfnpr,1060) nam,nctr
            end if
          else
            if(.not.core) then
              write(lfnpr,1050) nam,nctr,((labec(k,j),j=1,2),occec(k),
     +                          k=1,jmax)
            else
              write(lfnpr,1060) nam,nctr,((labec(k,j),j=1,2),occec(k),
     +                          k=1,jmax)
            end if
          end if
  899   continue
      end if
      if(jprint(4).lt.0) stop
c
c  write out wiberg bond index matrix if requested:
c
      if(jprint(12).ne.0) then
        title = 'wiberg bond index matrix in the nao basis:'
        call aout(bindex,natoms,natoms,natoms,title,0,natoms)
        do 3010 iat = 1,natoms
          bindt(iat) = zero
          do 3000 jat = 1,natoms
            if(iat.eq.jat) goto 3000
            bindt(iat) = bindt(iat) + bindex(jat,iat)
 3000     continue
 3010   continue
        title = 'wiberg bond index, totals by atom:'
        call aout(bindt,natoms,natoms,1,title,0,1)
c
c  write out overlap-weighted bond populations:
c
        title = 'atom-atom overlap-weighted nao bond order:'
        call aout(ovpop,natoms,natoms,natoms,title,0,natoms)
        do 3030 iat = 1,natoms
          bindt(iat) = zero
          do 3020 jat = 1,natoms
            if(iat.eq.jat) goto 3020
            bindt(iat) = bindt(iat) + ovpop(jat,iat)
 3020     continue
 3030   continue
        title(1:43)  = 'atom-atom overlap-weighted nao bond order, '
        title(44:58) = 'totals by atom:'
        call aout(bindt,natoms,natoms,1,title,0,1)
c
c  write out mo bond orders:
c
        title = 'mo bond order:'
        call aout(bmo,natoms,natoms,natoms,title,0,natoms)
        do 3050 iat = 1,natoms
          bindt(iat) = zero
          do 3040 jat = 1,natoms
            if(iat.eq.jat) goto 3040
            bindt(iat) = bindt(iat) + bmo(jat,iat)
 3040     continue
 3050   continue
        title  = 'mo atomic valencies:'
        call aout(bindt,natoms,natoms,1,title,0,1)
      end if
c
c  save nao info in common/nbnao/:
c
      do 888 i = 1,nbas
        naoc(i) = naoctr(i)
        naoa(i) = naol(i)
  888 continue
      return
c
  900 format(//,1x,
     +'NATURAL POPULATIONS:  Natural atomic orbital occupancies ',/,1x,
     +'                                                         ',/,1x,
     +'     NAO  atom    #  lang   Type(AO)    Occupancy      Energy'
     +,/,1x,
     +'--------------------------------------------------------------')
  910 format(//,3x,
     +'NATURAL POPULATIONS:  Natural atomic orbital occupancies ',/,1x,
     +'                                                         ',/,3x,
     +'   NAO  atom    #  lang   Type(AO)    Occupancy           ',/,1x,
     +'-------------------------------------------------              ')
  920 format(1x,i8,1x,a4,i6,2x,a1,a4,2x,a3,'(',i2,a1,')',4x,
     + f8.5,4x,f10.5)
  930 format(/,1x,
     +'[',i3,' electrons found in the effective core potential]')
  931 format(/,1x,
     +'WARNING:  1 low occupancy (<',f6.4,'e) core orbital  found ',
     +'on ',a2,i2)
  932 format(1x,
     +'          1 low occupancy (<',f6.4,'e) core orbital  found ',
     +'on ',a2,i2)
  933 format(/,1x,
     +'WARNING:',i3,' low occupancy (<',f6.4,'e) core orbitals found',
     +' on ',a2,i2)
  934 format(1x,
     +'        ',i3,' low occupancy (<',f6.4,'e) core orbitals found',
     +' on ',a2,i2)
  936 format(/,1x,
     +'WARNING:  population inversion found on atom ',a2,i2)
  937 format(1x,
     +'          population inversion found on atom ',a2,i2)
  939 format(//,1x,
     +'Summary of Natural Population Analysis:                  ',/,1x,
     +'                                                         ',/,1x,
     +'                                      Natural Population ',/,1x,
     +'                   Natural   ',45('-'),/,1x,2x,'Atom        #',
     +4x,'Charge',9x,'Core',6x,'Valence',4x,'Rydberg',6x,'Total',/,1x,
     +74('-'))
  940 format(1x,2x,a8,i5,2x,f9.5,4x,f9.5,3x,f9.5,2x,f9.5,3x,f9.5)
  950 format(1x,74('='),/,1x,'  * Total *',6x,f9.5,2x,f11.5,1x,f11.5,1x,
     + f10.5,1x,f11.5)
  955 format(/1x,
     +'number of electrons is not an integer!  please check your ',
     +'data.',/)
  956 format(/1x,
     +'WARNING: number of electrons is not within 1.0d-5 of an',
     +' integer.'/)
  960 format(/,1x,
     +'                                Natural Population      ',/,1x,
     +'--------------------------------------------------------')
  970 format(1x,'  Effective Core          ',f10.5)
  980 format(1x,'  Core                    ',f10.5,' (',f8.4,
     +'% of ',i3,')')
  990 format(1x,'  Valence                 ',f10.5,' (',f8.4,
     +'% of ',i3,')')
 1000 format(1x,'  Natural minimal basis   ',f10.5,' (',f8.4,
     +'% of ',i3,')')
 1010 format(1x,'  Natural rydberg basis   ',f10.5,' (',f8.4,
     +'% of ',i3,')',/,1x,
     +'--------------------------------------------------------')
 1040 format(/1x,
     +'   Atom    #          Natural Electron Configuration',/,1x,
     + 76('-'))
 1050 format(1x,2x,a4,i6,6x,6x,(13(i1,a1,'(',f5.2,')')))
 1060 format(1x,2x,a4,i6,6x,'[core]',(13(i1,a1,'(',f5.2,')')))
 2500 format(/1x,'subroutine naoanl could not find a ',a1,'-type ',
     + 'core orbital on atom ',a2,i2,'.',/,1x,'icore :',4i3,
     + '     m :',i3,'     la :',i3)
 2600 format(/1x,'subroutine naoanl could not find a ',a1,'-type ',
     + 'valence orbital on atom ',a2,i2,'.',/,1x,'ival :',4i3,
     + '     m :',i3,'     la :',i3)
      end
c*****************************************************************************
      subroutine frmtmo(t,tmo,c,scr,index,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      character*80 title
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbnao/naoc(maxbas),naoa(maxbas),ltyp(maxbas),
     +       iprin(maxbas)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      dimension t(ndim,ndim),tmo(ndim,ndim),c(ndim,ndim),
     +          scr(ndim*(ndim+5))
      character *4 basis(4)
c
      data basis/' nao',' nho',' nbo','nlmo'/
      data zero/0.0d0/
c
c  input:
c     t     --  transformation from ao basis to currect basis
c     index --  current basis = 2,3,4,5 (nao,nho,nbo,nlmo)
c     iflg  --  number of columns of tmo to print
c               or external lfn to write to
c
c  fetch the ao to mo transformation matrix:
c
      call feaomo(c,it)
      if(it.eq.0) return
c
c  find the mo transformation matrix:
c
      zertol = 1.0d-8
      eps    = 1.0d-8
      maxit  = 10
      lfn0   = 0
      call lineq(t,tmo,c,scr,nbas,nbas,ndim,ndim,zertol,eps,maxit,
     +           lfn0,ierr)
      if(ierr.ne.0) then
        write(lfnpr,910) basis(index-1)
        if(ierr.eq.1) write(lfnpr,920) basis(index-1)
        stop
      end if
c
c  make sure the largest coefficient in each column is positive:
c
      do 30 kcol = 1,nbas
        tmax = zero
        do 10 jrow = 1,nbas
          if(abs(tmo(jrow,kcol)).gt.abs(tmax)) tmax = tmo(jrow,kcol)
   10   continue
        if(tmax.lt.zero) then
          do 20 jrow = 1,nbas
            tmo(jrow,kcol) = -tmo(jrow,kcol)
   20     continue
        end if
   30 continue
c
c  write or print the mo transformation matrix:
c
      if(index.eq.2) title = 'mos in the nao basis:'
      if(index.eq.3) title = 'mos in the nho basis:'
      if(index.eq.4) title = 'mos in the nbo basis:'
      if(index.eq.5) title = 'mos in the nlmo basis:'
      call aout(tmo,ndim,nbas,nbas,title,index,iflg)
      return
c
  910 format(/1x,'error calculating the ',a4,' to mo transformation')
  920 format(1x,'the ao to ',a4,' transformation is not invertible')
      end
c****************************************************************************
c
c  routines called by sr nao:
c
c      subroutine loadav(listao,nl,m,s,ndim,a,b,mxaolm)
c      subroutine atdiag(n,a,b,eval,c)
c      subroutine setbas(lstocc,lstemt,nocc,nemt,iat,l,nl,nf,ndim)
c      subroutine newwts(s,t,wt)
c      subroutine worth(s,t,blk,list,ndim,nbas,n,occ,eval,bigblk)
c      subroutine NBOshmdt(t,s,ndim,nbas,nocc,lstocc,nemt,lstemt,sblk)
c      subroutine newryd(t,s,tpnao,dmblk,sblk,evect,occ,eval,eval2,
c     +                       list,irpnao)
c      subroutine rydiag(t,s,tpnao,dmblk,sblk,occ,eval,evect,eval2,
c     +                    iorb,nc,nm,nstart,nrydc,larc,list,irpnao)
c      subroutine rydsel(lstemt,nemt,nsel1,list1,nsel2,list2,wt)
c      subroutine rediag(dm,t,tpnao,eval,blk,c,irank,irpnao)
c      subroutine redblk(t,tpnao,il,dm,blk,eval,c,nf,iorb,nc,irank,irpnao)
c
c*****************************************************************************
      subroutine loadav(listao,nl,m,s,ndim,a,b,mxaolm)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension s(ndim,ndim),listao(mxaolm,9),a(nl,nl),b(nl,nl)
      data one,zero/1.0d0,0.0d0/
c
c  average the ao density matrix elements over the m=2*l+1 components
c  of l for a particular atom.
c  load density matrix elements (upper triangle of s, incl. diagonal)
c  into a, overlap matrix elements (lower triangle of s) into b, for
c  orbitals of 'list'
c
      do 30 j=1,nl
        do 20 i=1,j
c  find average dm element over the values of im:
          sum=zero
          do 10 im=1,m
            iao=listao(i,im)
            jao=listao(j,im)
   10       sum=sum+s(iao,jao)
          ave=sum/m
c  density matrix elements into a:
          a(i,j)=ave
          a(j,i)=ave
c  overlap matrix elements into b:
          b(i,j)=s(jao,iao)
   20     b(j,i)=b(i,j)
   30   b(j,j)=one
      return
      end
c*****************************************************************************
      subroutine atdiag(n,a,b,eval,c)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  solve generalized eigenvalue problem (a-eval*b)*c = 0.
c
c  use jacobi to diagonalize b**(-1/2)*a*b**(-1/2); a and b are destroyed.
c
      dimension a(n,n),b(n,n),eval(n),c(n,n)
      data zero,one/0.0d0,1.0d0/
c  first form b**(-1/2) and store it in b:
      call NBOjacobi(n,b,eval,c,n,n,0)
      do 10 i=1,n
   10   eval(i)=one/sqrt(eval(i))
      do 30 i=1,n
        do 30 j=1,i
          temp=zero
          do 20 k=1,n
   20       temp=temp+eval(k)*c(i,k)*c(j,k)
          b(i,j)=temp
   30     b(j,i)=temp
c  now similarity transform a with b:
      call simtrs(a,b,eval,n,n)
c  diagonalize a:
      call NBOjacobi(n,a,eval,c,n,n,1)
c  multiply b*c to get eigenvectors for original problem, store in a:
      do 50 i=1,n
        do 50 j=1,n
          temp=zero
          do 40 k=1,n
   40       temp=temp+b(i,k)*c(k,j)
   50     a(i,j)=temp
c  move final eigenvectors to c:
      call copy(a,c,n,n,n)
      return
      end
c*****************************************************************************
      subroutine setbas(lstocc,lstemt,nocc,nemt,iat,l,nl,nf,ndim)
c*****************************************************************************
c
c  select the set of natural minimal basis (nmb) orbitals for a particular
c  atom and angular symmetry type:  (up to atomic number 105)
c
c------------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension lstocc(ndim),lstemt(ndim)
      dimension icore(4),ival(4)
c
c  if g orbitals or orbitals of even higher angular symmetry are selected,
c  there are none in the nmb:
c
      if(l.ge.4) goto 100
c
c  find core and valence orbitals for this atom:
c
      iecp = 0
      call cortbl(iat,icore,iecp)
      call valtbl(iat,ival)
c
c  determine the number of shells with angular symmetry l in the nmb.
c  if there are a negative number of core orbitals, ignore them:
c
      nshell = max(icore(l+1),0) + ival(l+1)
      if(nshell.eq.0) goto 100
c
c  select sets of occupied and empty nao's:
c
      do 10 j = 1,nshell
        nocc = nocc + 1
        lstocc(nocc) = nf + j
   10 continue
      left = nl - nshell
      if(left.eq.0) return
      do 20 j = 1,left
        nemt = nemt + 1
        lstemt(nemt) = nf + nshell + j
   20 continue
      return
c
c  no nmb l-type orbitals found for this atom:
c
  100 continue
      do 110 j = 1,nl
        nemt = nemt + 1
        lstemt(nemt) = nf + j
  110 continue
      return
      end
c*****************************************************************************
      subroutine newwts(s,t,wt)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      dimension t(ndim,ndim),s(ndim,ndim),wt(ndim)
      character*80 title
c
      data zero/0.0d0/
c
c  recompute occupancy weights
      nocc=0
      do 20 i=1,nbas
        sum=zero
        do 10 j=1,nbas
          do 10 k=1,nbas
            sjk=s(j,k)
            if(j.gt.k) sjk=s(k,j)
   10       sum=sum+t(j,i)*sjk*t(k,i)
        wt(i)=sum
c  reformat list lstocc:
        if(lstocc(i).eq.0) go to 20
        nocc=nocc+1
        lstocc(nocc)=i
   20   continue
      nstart=nocc+1
      do 40 i=nstart,ndim
   40   lstocc(i)=0
c  symmetry-average weights:
      nl=1
      iorb=0
  100 iorb=iorb+nl
        if(iorb.gt.nbas) go to 600
        nl=1
        ilbl=naoctr(iorb)
        il=naol(iorb)/100
        nm=il*2+1
        imax=nbas-iorb
        do 130 iadd=1,imax
          jorb=iorb+iadd
          jorbl=naol(jorb)/100
          if(naoctr(jorb).ne.ilbl.or.jorbl.ne.il) go to 140
  130     nl=nl+1
  140   nc=nl/nm
        do 500 i=1,nc
          sum=zero
          do 300 m=1,nm
            inao=iorb+(i-1)+(m-1)*nc
  300       sum=sum+wt(inao)
          av=sum/nm
          do 400 m=1,nm
            inao=iorb+(i-1)+(m-1)*nc
  400       wt(inao)=av
  500     continue
      go to 100
c
  600 continue
      title = 'new symmetry-averaged occupancy weights:'
      call aout(wt,nbas,nbas,1,title,-1,1)
      return
c
      end
c*****************************************************************************
      subroutine worth(s,t,blk,list,ndim,nbas,n,occ,eval,bigblk)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c******************************************************************
c
c   worth: occupancy weighted orthogonalization subroutine
c
c   s:           full overlap matrix (pure ao basis)
c                 (note: upper triangle used for scratch, but restored again)
c   t:           pure ao to pre-nao transformation
c   list:        list of orbitals to be orthogonalized
c   n:           number of orbitals in list
c   occ:         list of symmetry averaged occupancy weightings
c
c   note:    blk and bigblk share the same storage but are
c               dimensioned differently.
c
c******************************************************************
      dimension s(ndim,ndim),t(ndim,ndim),blk(n,n)
      dimension occ(ndim),list(ndim),eval(ndim),bigblk(ndim,ndim)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      data zero,one/0.0d0,1.0d0/
      data ntime/0/
c
c  important constants:
c           wtthr           all weighting factors smaller than wtthr are set
c                            to the value of wtthr.
c           diagth          threshold for matrix diagonalization used in
c                            subroutine jacobi.  in jacobi, this constant
c                            is called "doneth".
c           danger          criterion for deciding that the job should be
c                            aborted due to numerical problems caused by near
c                            linear dependencies in the basis set.  all
c                            eigenvalues of the weighted overlap matrix must
c                            be greater than diagth*danger.
c
      data wtthr,diagth,danger/1.0d-3,1.0d-12,1.0d2/
c
      ntime=ntime+1
c  multiply the weight by a constant so that the maximum weight is one,
c   and set any resulting weight that is less than wtthr to the value of wtthr:
      wtmax=zero
      do 10 i=1,n
        ip=list(i)
        if(occ(ip).gt.wtmax) wtmax=occ(ip)
   10   continue
      do 20 i=1,n
        ip=list(i)
        eval(ip)=occ(ip)/wtmax
        if(eval(ip).lt.wtthr) eval(ip)=wtthr
   20   continue
c  form the weighted pre-nao vectors:
      do 30 j=1,n
        jp=list(j)
        do 30 i=1,nbas
   30     t(i,jp)=t(i,jp)*eval(jp)
c  form the weighted overlap matrix of the vectors in the upper triangle of s:
      do 70 i=1,n
        ip=list(i)
        do 70 j=1,nbas
          sij=zero
          do 40 k=1,nbas
            tki=t(k,ip)
            if(tki.eq.zero) go to 40
            sij=sij+tki*s(k,j)
   40       continue
   70     bigblk(j,i)=sij
      do 100 i=1,n
        do 100 j=1,i
          jp=list(j)
          sij=zero
          do 90 k=1,nbas
            tkj=t(k,jp)
            if(tkj.eq.zero) go to 90
            sij=sij+bigblk(k,i)*tkj
   90       continue
  100     s(j,i)=sij
c  diagonalize s-tilde (the weighted overlap matrix):
      call NBOjacobi(n,s,eval,blk,ndim,n,0)
c
c  form the inverse sqrt of the overlap matrix of these weighted vectors:
      smlest=one
      toosml=diagth*danger
      do 150 i=1,n
        eigenv=eval(i)
        if(eigenv.lt.toosml) go to 900
        eval(i)=one/sqrt(eigenv)
        if(eigenv.lt.smlest) smlest=eigenv
  150  continue
      do 170 i=1,n
        do 170 j=1,i
          sij=zero
          do 160 k=1,n
  160       sij=sij+eval(k)*blk(i,k)*blk(j,k)
  170     s(j,i)=sij
c
c  the upper triangle of s (including the diagonal)
c   now contains the -0.5 power of the weighted overlap matrix,
c   and is the weighted orthog. transform that we want.
c   now, form the total transformation:
      do 300 i=1,nbas
        do 260 j=1,n
          eval(j)=zero
          do 220 k=1,j
            kp=list(k)
            tik=t(i,kp)
            if(tik.eq.zero) go to 220
            eval(j)=eval(j)+tik*s(k,j)
  220       continue
          jp1=j+1
          do 240 k=jp1,n
            kp=list(k)
            tik=t(i,kp)
            if(tik.eq.zero) go to 240
            eval(j)=eval(j)+tik*s(j,k)
  240       continue
  260     continue
        do 280 j=1,n
          jp=list(j)
  280     t(i,jp)=eval(j)
  300   continue
c  restore full overlap matrix s:
      do 400 i=1,nbas
        im1=i-1
        do 380 j=1,im1
  380     s(j,i)=s(i,j)
  400   s(i,i)=one
      return
c
  900 write(lfnpr,1000) eigenv,toosml
      stop
c
 1000 format(//1x,'an eigenvalue of the weighted pre-nao overlap',
     +' matrix of ',1pe13.4,' has been',/1x,'found, which is lower than'
     +,' the allowed threshold of ',1pe13.4,'.  this is',/,1x,'probably',
     +' caused by either an error in the data given to the analysis',
     +' program',/,1x,'or by numerical problems caused by near linear',
     +' dependencies among the basis',/,1x,'functions.')
      end
c*****************************************************************************
      subroutine NBOshmdt(t,s,ndim,nbas,nocc,lstocc,nemt,lstemt,sblk)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  schmidt orthogonalization of column vectors in t
c     schmidt orthogonalize each empty orbital (specified in 'lstemt')
c        to the orthonormal occupied (lstocc) orbitals;
c
      dimension t(ndim,ndim),s(ndim,ndim),lstocc(ndim),lstemt(ndim),
     *              sblk(ndim,ndim)
      data zero/0.0d0/
      do 30 i=1,nbas
        do 30 j=1,nocc
          jp=lstocc(j)
          sji=zero
          do 10 k=1,nbas
   10       sji=sji+t(k,jp)*s(k,i)
   30     sblk(i,j)=sji
c   schmidt orthogonalize each unoccupied /ui> to each /vj>:
c...loop over unoccupied /ui>'s,
      do 120 i=1,nemt
        ip=lstemt(i)
c...loop over occupied /vj>'s,
        do 60 j=1,nocc
          jp=lstocc(j)
c...calculate sji = <vj/ui>,
          sji=zero
          do 40 k=1,nbas
   40       sji=sji+sblk(k,j)*t(k,ip)
c...and replace each /ui> = /ui> - sji*/vj>.
          do 50 k=1,nbas
   50       t(k,ip)=t(k,ip)-sji*t(k,jp)
   60     continue
  120   continue
      return
      end
c*****************************************************************************
      subroutine newryd(t,s,tpnao,dmblk,sblk,evect,occ,eval,eval2,
     *                       list,irpnao)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      dimension t(ndim,ndim),s(ndim,ndim),tpnao(ndim,ndim),occ(ndim),
     +          dmblk(mxaolm,mxaolm),sblk(mxaolm,mxaolm),eval(nbas),
     +          evect(mxaolm,mxaolm),eval2(nbas),list(mxaolm)
      data one/1.0d0/
c
c  compute new rydberg naos after the schmidt orthogonalization to
c  the minimal nao set has been done:
c
c  if requested (irpnao=jprint(11)=1), update pnao transformation with tryd:
c
      if(irpnao.eq.1) call fepnao(tpnao)
c
      nl=1
      iorb=0
  100 iorb=iorb+nl
        if(iorb.gt.nbas) go to 300
        nl=1
        ilbl=naoctr(iorb)
        il=naol(iorb)/100
        nm=il*2+1
        imax=nbas-iorb
        do 130 iadd=1,imax
          jorb=iorb+iadd
          jorbl=naol(jorb)/100
          if(naoctr(jorb).ne.ilbl.or.jorbl.ne.il) go to 140
  130     nl=nl+1
  140   nc=nl/nm
        nskip=0
        imax=iorb-1+nc
        do 150 i=1,nbas
          inao=lstocc(i)
          if(inao.lt.iorb.or.inao.gt.imax) go to 150
          nskip=nskip+1
  150     continue
        if(nskip.eq.nc) go to 100
        nstart=nskip+1
        nrydc=nc-nskip
        call rydiag(t,s,tpnao,dmblk,sblk,occ,eval,evect,eval2,
     *              iorb,nc,nm,nstart,nrydc,larc,list,irpnao)
c  end of loop starting at 100
        go to 100
  300 continue
c  restore s:
      do 350 i=1,nbas
        im1=i-1
        do 340 j=1,im1
  340     s(j,i)=s(i,j)
  350   s(i,i)=one
c
c  save updated t-pnao transformation:
c
      if(irpnao.eq.1) call svpnao(tpnao)
      return
      end
c*****************************************************************************
      subroutine rydiag(t,s,tpnao,dmblk,sblk,occ,eval,evect,eval2,
     *                    iorb,nc,nm,nstart,nrydc,larc,list,irpnao)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension t(ndim,ndim),s(ndim,ndim),tpnao(ndim,ndim),occ(nbas),
     *        dmblk(nrydc,nrydc),sblk(nrydc,nrydc),eval(nbas),
     *        evect(nrydc,nrydc),larc(nrydc),list(nrydc),eval2(nbas)
      data zero/0.0d0/
c
c  diagonalize one rydberg block, update t-nao (in t) and, if irpnao.eq.1,
c  update tpnao:
c
      ii=0
      do 20 i=1,nrydc
        do 20 j=1,nrydc
          dmblk(i,j)=zero
          sblk(i,j)=zero
   20     continue
      do 500 i=nstart,nc
        ii=ii+1
        do 300 m=1,nm
          inao=iorb+(i-1)+(m-1)*nc
          do 140 k=1,nbas
            dmsum=zero
            ssum=zero
            km1=k-1
            do 100 l=1,km1
              tli=t(l,inao)
              dmsum=dmsum+tli*s(l,k)
  100         ssum=ssum+tli*s(k,l)
            tki=t(k,inao)
            dmsum=dmsum+tki*s(k,k)
            ssum=ssum+tki
            kp1=k+1
            do 120 l=kp1,nbas
              tli=t(l,inao)
              dmsum=dmsum+tli*s(k,l)
  120         ssum=ssum+tli*s(l,k)
            eval(k)=dmsum
            eval2(k)=ssum
  140       continue
          jj=0
          do 240 j=nstart,i
            jj=jj+1
            jnao=iorb+(j-1)+(m-1)*nc
            dmsum=zero
            ssum=zero
            do 200 k=1,nbas
              tkj=t(k,jnao)
              dmsum=dmsum+eval(k)*tkj
  200         ssum=ssum+eval2(k)*tkj
            dmblk(ii,jj)=dmblk(ii,jj)+dmsum
            sblk(ii,jj)=sblk(ii,jj)+ssum
  240       continue
  300     continue
          do 350 jj=1,ii
            dmblk(ii,jj)=dmblk(ii,jj)/nm
            dmblk(jj,ii)=dmblk(ii,jj)
            sblk(ii,jj)=sblk(ii,jj)/nm
  350       sblk(jj,ii)=sblk(ii,jj)
  500     continue
      call atdiag(nrydc,dmblk,sblk,eval,evect)
      call rank(eval,nrydc,nrydc,larc)
      do 600 j=1,nrydc
        jc=larc(j)
        do 600 i=1,nrydc
  600     sblk(i,j)=evect(i,jc)
      do 700 m=1,nm
        jj=0
        do 680 j=nstart,nc
          jj=jj+1
          jnao=iorb+(j-1)+(m-1)*nc
          occ(jnao)=eval(jj)
          list(jj)=jnao
  680     continue
c  use limtrn to update t:
        call limtrn(t,list,sblk,dmblk,ndim,nbas,nrydc,nrydc,1)
  700 continue
c
      if(irpnao.eq.0) return
c
c  update tpnao, but do this in such a way that the intra-atomic blocks
c  of the overlap matrix in the revised pnao matrix remain diagonal
c  and that the pnaos remain normalized.   in order to accomplish this,
c  we must lowdin-orthogonalize the rydberg transformation in "sblk":
c
      call NBOsymort(evect,sblk,dmblk,nrydc,nrydc,eval)
      do 800 m=1,nm
        jj=0
        do 780 j=nstart,nc
          jj=jj+1
          list(jj)=iorb+(j-1)+(m-1)*nc
  780     continue
        call limtrn(tpnao,list,sblk,dmblk,ndim,nbas,nrydc,nrydc,1)
  800   continue
      return
      end
c*****************************************************************************
      subroutine rydsel(lstemt,nemt,nsel1,list1,nsel2,list2,wt)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension wt(ndim),list1(nbas),list2(nbas),lstemt(nbas)
      data one,wtthr/1.0d0,1.0d-4/
c
c   divide the rydberg orbitals into 2 groups:
c      list1:     rydbergs of significant occupancy ( .gt.wtthr )
c
c      list2:     rydbergs of very low occupancy ( .lt.wtthr )
c
c      wtthr is set to 0.0001
c
c    set the weightings of the rydbergs in list2 to one so that the weighted
c      orthogonalization that will later be done among these orbitals will
c      be in fact a lowdin orthog.
c
      nsel1=0
      nsel2=0
      do 100 i=1,nemt
        iryd=lstemt(i)
        if(wt(iryd).lt.wtthr) go to 50
          nsel1=nsel1+1
          list1(nsel1)=iryd
          go to 100
   50   continue
          nsel2=nsel2+1
          list2(nsel2)=iryd
          wt(iryd)=one
  100   continue
      return
      end
c*****************************************************************************
      subroutine rediag(dm,t,tpnao,eval,blk,c,irank,irpnao)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/ldeg(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      dimension dm(ndim,ndim),t(ndim,ndim),tpnao(ndim,ndim),
     +  c(mxaolm,mxaolm),eval(ndim),blk(mxaolm,mxaolm),irank(nbas)
c
c  rediagonalize the symmetry averaged dm subblocks for each angular
c  symmetry on each atom:
c
c  read in old t-pnao into tpnao so that it can be updated (if irpnao.eq.1):
c
      if(irpnao.eq.1) call fepnao(tpnao)
      nf = 0
      iorb = 0
      nl = 1
   10 iorb = iorb + nl
        if(iorb.gt.nbas) go to 100
        nl = 1
        ilbl = naoctr(iorb)
        il = naol(iorb)/100
        nm = il*2 + 1
        imax = nbas - iorb
        do 30 iadd = 1,imax
          jorb = iorb + iadd
          jorbl = naol(jorb)/100
          if((naoctr(jorb).ne.ilbl).or.(jorbl.ne.il)) go to 40
   30     nl = nl + 1
   40   nc = nl/nm
        if(nc.eq.1) go to 80
        call redblk(t,tpnao,il,dm,blk,eval,c,nf,iorb,nc,irank,
     *                 irpnao)
        go to 10
   80   do 90 m = 1,nm
          nf = nf + 1
   90     continue
        go to 10
  100 continue
      if(irpnao.eq.0) return
c
c  save new t-pnao from tpnao:
c
      call svpnao(tpnao)
      return
      end
c*****************************************************************************
      subroutine redblk(t,tpnao,il,dm,blk,eval,c,nf,iorb,nc,irank,
     *                           irpnao)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/ldeg(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      dimension dm(ndim,ndim),blk(nc,nc),c(nc,nc),eval(ndim),
     +  t(ndim,ndim),tpnao(ndim,ndim),irank(nbas)
      data zero/0.0d0/
c
c  find the rediagonalization transformation for the dm subblock for
c  the angular momentum "il" on an atom, put in t2:
c
      nm = il*2 + 1
      do 30 j = 1,nc
        do 30 i = 1,j
          sum = zero
          do 10 m = 1,nm
            inao = iorb + i-1 + (m-1)*nc
            jnao = iorb + j-1 + (m-1)*nc
   10       sum = sum + dm(inao,jnao)
          ave = sum/nm
          blk(i,j) = ave
   30     blk(j,i) = ave
      call NBOjacobi(nc,blk,eval,c,nc,nc,1)
      call rank(eval,nc,nc,larc)
      do 80 j = 1,nc
        jc = larc(j)
        do 80 i = 1,nc
   80     blk(i,j) = c(i,jc)
      do 110 m = 1,nm
        do 100 j = 1,nc
          nf = nf + 1
  100     irank(j) = nf
        call limtrn(t,irank,blk,c,ndim,nbas,nc,nc,1)
        call limtrn(dm,irank,blk,c,ndim,nbas,nc,nc,0)
        if(irpnao.eq.1) call limtrn(tpnao,irank,blk,c,ndim,nbas,nc,nc,1)
  110 continue
c
      return
      end
c****************************************************************************
c
c  routines called by the nbo/nlmo drivers:
c
c      subroutine nathyb(dm,t,guide,bndocc,pol,q,v,blk,c,eval,borb,
c     +                                       p,ta,hyb,va,vb,topo)
c      subroutine chsdrv(dm,t,guide,bndocc,pol,q,v,blk,c,eval,borb,
c     +                                       p,ta,hyb,va,vb,topo)
c      subroutine choose(dm,t,guide,bndocc,pol,q,v,blk,c,eval,borb,
c     +                                  p,ta,hyb,va,vb,topo,iflg)
c      subroutine srtnbo(t,bndocc)
c      subroutine xcited(dm,t,hyb,thyb,s,occ,scr,iscr)
c      subroutine anlyze(t,bndocc,hyb,hycoef,thyb)
c      subroutine htype(hyb,ltyp,mxao,nh,coef,pct,nl,isgn)
c      subroutine frmhyb(hyb,thyb,coef,hycoef,kl,ku,nhyb)
c      subroutine hybdir(bndocc,atcoor,thyb,tbnd,scr)
c      subroutine hybcmp(xyz,pcent,ihyb,jctr,hyb)
c      subroutine fndmol(iatoms)
c      subroutine nbocla(bndocc,accthr)
c      subroutine fnboan(bndocc,f,molnbo)
c      subroutine nbosum(f,bndocc,list,lista,scr)
c      subroutine getdel(ibo,occ,thr1,thr2,nl,list,del,deloc,iflg)
c      subroutine dlcstr(ibo,il,nl,list,ml,istr)
c      subroutine nlmo(n,a,eval,evec,tsym,reson,nocc,ialarm)
c      subroutine lmoanl(t,s,reson,occ,ts,border,owbord,atlmo,siab,nocc,nab)
c      subroutine dipanl(dm,t,c,tnbo,dx,dy,dz,scr,index)
c      subroutine dipele(dxyz,c,t,scr,eta,nocc,index)
c      subroutine NBOdipnuc(dx,dy,dz,atcoor,eta,nocc)
c
c****************************************************************************
      subroutine nathyb(dm,t,guide,bndocc,pol,q,v,blk,c,eval,borb,
     *                                       p,ta,hyb,va,vb,topo)
c****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  construct orthogonal matrix t for transformation from ao's to
c  natural hybrid bond orbitals using input density matrix dm.
c
c  required input includes:
c        dm = density matrix in orthonormal atomic orbital basis;
c                  real(1,ndim;1,ndim)
c      nbas = no. of orbitals = actual dimension of dm,s,t,naol,dmt
c    natoms = no. of atoms (not including ghosts) in the molecule
c     iatno = list of atomic numbers
c    naoctr = orbital label list.  naoctr(i)=iat if nao # i is on atom iat
c                integer(1,ndim).  naos of given atom grouped together.
c      iw3c = 1 if program is to search for 3-center bonds,
c           = 0 otherwise
c     guide = wiberg atom-atom bond index matrix, used as guide for nbo search
c
c  output:
c         t = bond orbital transformation matrix (ndim,ndim).
c                rows are labelled by naos, columns by nbos.
c     label = list of bond orbital labels
c      ibxm = permutation list of bond orbital labels (very important!)
c
      logical detail,nobond,first
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      integer ul
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       ul(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbnao/naoctr(maxbas),naol(maxbas),ltyp(maxbas),
     +       iprin(maxbas)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nbtopo/iorder(maxatm),jorder(maxatm),ntopo(maxatm,maxatm),
     +            n3ctr,i3ctr(10,3)
      dimension dm(ndim,ndim),t(ndim,ndim),v(ndim),borb(mxbo),
     * pol(ndim,3),bndocc(ndim),name(3),hybexp(3),
     * q(mxao,ndim),blk(mxbo,mxbo),eval(mxbo),c(mxbo,mxbo),
     * p(mxao,mxao),ta(mxao,mxao),hyb(mxao),va(mxao),vb(mxao),
     * guide(natoms,natoms),topo(natoms*natoms)
      data gthrsh/1.5d-1/
      data istar,iblnk/1h*,1h /
      data name/2hlp,2hbd,2h3c/
      data lry,lcr/2hry,2hcr/
      data zero,zerop,tenth,one,two,four
     * /0.d0,1.d-6,0.1d0,1.0d0,2.0d0,4.0d0/
      data twop/2.0001d0/
      data pt8,pt99/0.8d0,0.99d0/
c
c  prjinc, the amount to increase prjthr by if problems with linear
c  dependency between the hybrids arise.
c
      data prjinc/0.05d0/
c
      nopval(i) = norbs(i) - ino(i)
c
      detail = .false.
      if(iwdetl.ne.0) detail = .true.
      nobond = .false.
      if(jprint(10).ne.0) nobond = .true.
c
c  initial iteration loop:  if no satisfactory lewis structure (all
c  antibond occupancies < 0.1) for thresh = 1.90, thresh is decremented
c  up to 4 times by 0.1 in search of a better structure.  if the dm is
c  not spinless, thresh is set to 0.90 and is decremented as above.
c
      prjthr = abs(prjset)
      thresh = abs(thrset)
      if(ispin.ne.0) thresh = thresh - one
      if(nobond) thresh = one
      if(nobond.and.(ispin.ne.0)) thresh = one/two
      if(ispin.ne.0) gthrsh = gthrsh/four
c
c  determine the atom ordering for the initial search for bonds:
c
      if(natoms.eq.1) then
        iorder(1) = 1
        goto 45
      end if
c
c  find the two atoms which have the largest bond index:
c
      gmax = zero
      do 10 j = 2,natoms
        do 5 i = 1,j-1
          if(guide(i,j).gt.gmax) then
            gmax = guide(i,j)
            iat  = i
          end if
    5   continue
   10 continue
      iorder(1) = iat
c
c  add atoms to iorder according to these connectivities:
c
      icnt = 1
      inxt = icnt
      jcnt = icnt
   15 iptr = inxt
        i1st = 1
        do 20 i = 1,natoms
          topo(i) = guide(i,iorder(iptr))
   20   continue
        call rank(topo,natoms,natoms,jorder)
        jptr = 1
   25   if(topo(jptr).gt.pt8) then
          iflg = 1
          do 30 i = 1,icnt
            if(iorder(i).eq.jorder(jptr)) iflg = 0
   30     continue
          if(iflg.eq.1) then
            icnt = icnt + 1
            iorder(icnt) = jorder(jptr)
            if(i1st.eq.1) then
              i1st = 0
              inxt = icnt
            end if
          end if
        else
          goto 35
        end if
        jptr = jptr + 1
        goto 25
c
   35   continue
        if(i1st.eq.1) then
          jcnt = jcnt + 1
          inxt = jcnt
          if(inxt.gt.natoms) goto 45
          if(inxt.gt.icnt) then
            kptr = 0
   40       kptr = kptr + 1
            kflg = 1
            do 41 i = 1,icnt
              if(iorder(i).eq.kptr) kflg = 0
   41       continue
            if(kflg.eq.0) goto 40
            icnt = icnt + 1
            iorder(icnt) = kptr
          end if
        end if
      goto 15
c
   45 continue
      iter   = 0
      ialarm = 0
   50 if(ialarm.eq.0) iter = iter + 1
c
c  store density matrix in upper triangle of t:
c
      do 60 j = 1,nbas
        do 60 i = 1,j
   60     t(i,j) = dm(i,j)
c
c  zero arrays q, pol, iathy, ino, and label:
c
      do 100 i = 1,nbas
        do 70 k = 1,2
   70     label(i,k) = iblnk
        do 80 k = 3,6
   80     label(i,k) = 0
        do 90 k = 1,3
          pol(i,k) = zero
   90     iathy(i,k) = 0
        do 100 k = 1,mxao
  100     q(k,i) = zero
      do 110 i = 1,natoms
  110   ino(i) = 0
c
c  remove core orbitals from the density matrix:
c
      ibd = 0
      call core(dm,t,borb,pol,q,hyb,bndocc,ibd,detail,lfnpr)
c
c  main nho loops
c  --------------
c  doubly occupied (iocc=1) or singly occupied (iocc=2) nho's
c  if ispin.ne.0, search is only for singly occupied nbos (iocc=1):
c
      occmx = thresh
c
c  main nho loops over singles, doubles, and triples of atoms:
c
      na1 = natoms + 1
      do 310 ia1 = 1,na1
        ia = ia1 - 1
        if((ia.gt.0).and.(nopval(iorder(ia)).le.0)) go to 310
        do 300 ib1 = 1,na1
          ib = ib1 - 1
          if((ib.gt.0).and.(nopval(iorder(ib)).le.0)) go to 300
          do 290 ic1 = 2,na1
            ic = ic1 - 1
            if((ic.gt.0).and.(nopval(iorder(ic)).le.0)) go to 290
            if(ia.ne.0) go to 130
            if(ib.ne.0) go to 120
c
c  lone pairs:
c
            nctr = 1
            iat1 = iorder(ic)
            iat2 = 0
            iat3 = 0
            go to 140
c
c  bond pairs:
c
  120       continue
            if(nobond) go to 290
            nctr = 2
            iat1 = iorder(ib)
            iat2 = iorder(ic)
            iat3 = 0
            if(iat1.ge.iat2) go to 290
            if(guide(iat1,iat2).lt.gthrsh) go to 290
            go to 140
c
c  3-center bonds:
c
  130       continue
            if(iw3c.ne.1) go to 320
            nctr = 3
            iat1 = iorder(ia)
            iat2 = iorder(ib)
            iat3 = iorder(ic)
            if(iat1.ge.iat2) go to 300
            if(iat2.ge.iat3) go to 290
            if(guide(iat1,iat2).gt.gthrsh) go to 140
            if(guide(iat1,iat3).gt.gthrsh) go to 140
            if(guide(iat2,iat3).gt.gthrsh) go to 140
            go to 290
  140       continue
c
c  deplete dm of one(two) center orbitals if search for two(three)
c  center orbitals is beginning:
c
            if(iwprj(nctr).ne.0)
     *            call deplet(dm,t,q,pol,borb,bndocc,ibd)
c
c  load proper atomic blocks of dm into blk:
c
            call NBOload(dm,iat1,iat2,iat3,blk,nb)
c
c  diagonalize blk:
c
            call NBOjacobi(nb,blk,eval,c,mxbo,mxbo,1)
c
c  rank eigenvectors by occupancy eigenvalue:
c
            call rank(eval,nb,mxbo,larc)
            if(detail) write(lfnpr,1400) iat1,iat2,iat3
            if(detail) write(lfnpr,1403) thresh
            if(detail) write(lfnpr,1405) (eval(irnk),irnk=1,nb)
            iaccep = 0
            do 250 irnk = 1,nb
              ir = larc(irnk)
              occ = eval(irnk)
              do 200 i = 1,nb
  200           borb(i) = c(i,ir)
              if(detail) write(lfnpr,1410) irnk,occ
              if(detail) write(lfnpr,1420) (borb(i),i=1,nb)
c
c  throw out orbital if occupancy is less than the threshhold "occmx":
c
              if(occ.lt.occmx) go to 280
c
c  check to see that bond orbital "borb" doesn't contain previously used
c  hybrids:
c
              if(nctr.eq.1) go to 240
              call prjexp(borb,iat1,iat2,iat3,q,p,ta,hyb,va,vb,hybexp)
              if(.not.detail) go to 220
              do 210 ihyb = 1,nctr
  210           write(lfnpr,1500) ihyb,hybexp(ihyb)
  220         continue
              do 230 ihyb = 1,nctr
  230           if(hybexp(ihyb).lt.prjthr) go to 250
  240         continue
              ibd = ibd + 1
              iaccep = iaccep + 1
c
c  decompose "borb" into its constituent atomic hybrids and store in q:
c
              call stash_nbo(borb,ibd,iat1,iat2,iat3,pol,q,hyb)
c
c  construct bond orbital labels:
c
              label(ibd,1) = name(nctr)
              label(ibd,2) = iblnk
              label(ibd,3) = iaccep
              label(ibd,4) = iat1
              label(ibd,5) = iat2
              label(ibd,6) = iat3
              bndocc(ibd) = occ
              if(detail) write(lfnpr,1600) ibd,(label(ibd,i),i=1,3)
  250         continue
  280       continue
  290       continue
  300     continue
  310   continue
  320 continue
c
c  symmetric orthogonalization of principal hybrids:
c
      call orthyb(q,blk,ta,eval,c,ialarm,0)
c
c   ialarm sounds the alarm that there is linear dependency between some
c   of the hybrids. the remedy is to increase prjthr and repeat the nbo
c   search. ialarm is equal to the number of the violating atom.
c
      if(ialarm.ne.0) then
        oldprj = prjthr
        prjthr = oldprj + prjinc
        if(jprint(5).ne.0) write(lfnpr,1800) oldprj,prjthr
        if(prjthr.ge.pt99) then
          write(lfnpr,1810) ialarm
          jprint(1) = -1
          return
        end if
        goto 700
      end if
c
c  augment open-valence atoms with non-arbitrary hybrids orthogonal to
c  those found previously:
c
      do 580 ia = 1,natoms
        if(nopval(ia).le.0) go to 580
c
c  iula: upper limit of naos on atom. find nmb, the number of natural
c  minimal basis functions on the atom:
c
        lla = ll(ia)
        iula = ul(ia)
        nmb = 0
        do 470 i = lla,iula
          if(lstocc(i).eq.1) nmb = nmb + 1
  470   continue
c
c  find the number of bond, core, and lone pair hybrids on the atom, iocc:
c  also find iocclp, number of lone pair orbitals already found on ia, for
c  use in labelling the extra lone pairs below:
c
        iocc = 0
        iocclp = 0
        do 480 ib = 1,ibd
          if((label(ib,4).ne.ia).and.(label(ib,5).ne.ia).and.
     *            (label(ib,6).ne.ia)) go to 480
          iocc = iocc + 1
          if(label(ib,1).eq.name(1)) iocclp = iocclp + 1
  480   continue
c
c  nexlp: number of extra (low occupancy) lp orbitals on atom iat. (this
c  is the number of low occupancy orbitals with valence shell character)
c  set nexlp to zero if (nmb-iocc) is less than zero in order that the
c  orbitals are not miscounted!!
c
        nexlp = nmb - iocc
        if(nexlp.lt.0) nexlp = 0
        nocc = ino(ia)
        call frmprj(p,ia,q,nocc,ta,va,vb)
        norb = norbs(ia)
        naugm = norb - nocc
        call NBOaugmnt(p,blk,c,eval,dm,ta,borb,v,larc,ia,nocc,norb)
c
c  stash_nbo and label extra lone pairs that augmnt put in blk: (these ar
c  taken to be the highest occupied orbitals, which augmnt places first)
c
        do 510 iaugm = 1,nexlp
          do 500 j = 1,norb
  500       borb(j) = blk(j,iaugm)
          ibd = ibd + 1
          call stash_nbo(borb,ibd,ia,0,0,pol,q,hyb)
          label(ibd,1) = name(1)
          label(ibd,2) = iblnk
          label(ibd,3) = iaugm + iocclp
          label(ibd,4) = ia
          label(ibd,5) = 0
          label(ibd,6) = 0
  510   continue
c
c  stash_nbo and label the rydberg orbitals that augmnt put in blk:
c
        iryd = 0
        nstart = nexlp + 1
        do 540 iaugm = nstart,naugm
          do 530 j = 1,norb
  530       borb(j) = blk(j,iaugm)
          ibd = ibd + 1
          iryd = iryd + 1
          call stash_nbo(borb,ibd,ia,0,0,pol,q,hyb)
          label(ibd,1) = lry
          label(ibd,2) = istar
          label(ibd,3) = iryd
          label(ibd,4) = ia
          label(ibd,5) = 0
          label(ibd,6) = 0
  540     continue
  580   continue
c
c  include antibond labels:
c
      ibo = ibd
      do 660 i = 1,ibo
c
c  exit loop if label(i,1) is 'lp', 'ry', or 'cr':
c
        if(label(i,1).eq.name(1)) go to 660
        if(label(i,1).eq.lry) go to 660
        if(label(i,1).eq.lcr) go to 660
         nab = 1
         if(label(i,1).eq.name(3)) nab = 2
         do 650 iab = 1,nab
           ibd = ibd + 1
           do 640 j = 1,6
  640        label(ibd,j) = label(i,j)
           label(ibd,2) = istar
  650      continue
  660   continue
c
c  replace density matrix dm from t:
c
  700 continue
      do 740 j=1,nbas
        do 740 i=1,j
          dm(i,j)=t(i,j)
          dm(j,i)=dm(i,j)
          t(j,i)=zero
  740     t(i,j)=zero
c
c  remember the alarm!
c
      if(ialarm.ne.0) go to 50
c
c  miscounted bond orbitals...exit for open shell:
c
      if(ibd.ne.nbas) then
        write(lfnpr,1200) thresh,ibd,nbas
        write(lfnpr,1210) (i,(label(i,j),j=1,6),i=1,ibd)
        stop
      end if
c
c  find new polarization parameters for orthonormal hybrids:
c
      call repol(dm,q,pol,blk,eval,c,ibd)
c
c  form final t-nab (nao to nbo transformation) from orthonormal
c  hybrids:
c
      call formt(t,q,pol)
c
c  find occupancies, find total number of electrons and occupied orbitals:
c
      totele = zero
      do 800 i = 1,nbas
        occi = zero
        do 790 j = 1,nbas
          do 790 k = 1,nbas
  790       occi = occi + t(j,i) * dm(j,k) * t(k,i)
        if(abs(occi).lt.zerop) occi = zero
        if(occi.gt.twop) go to 960
        zeropm = -zerop
        if(occi.lt.zeropm) go to 960
        bndocc(i) = occi
        v(i) = occi
        totele = totele + bndocc(i)
  800 continue
      nel = totele + tenth
      if(abs(totele-nel).gt.5e-4) go to 970
      totele = nel
      nocc = nel
      if(ispin.eq.0) nocc = nocc/2 + mod(nocc,2)
c
c  make sure all but the nocc highest occupied nbos are starred:
c
      call rank(v,nbas,ndim,larc)
      do 804 i = 1,nocc
        ir = larc(i)
        label(ibxm(ir),2) = iblnk
  804 continue
      do 805 i = nocc+1,nbas
        ir = larc(i)
        label(ibxm(ir),2) = istar
  805 continue
c
c  determine whether this is a good resonance structure:
c
      call cycles(iter,thresh,guide,bndocc,topo,icont)
      if(icont.eq.0) then
        jprint(1) = -1
        return
      end if
      if(icont.eq.-1) go to 50
      if(icont.eq.1) go to 50
c
c  before final return, write out info about core orbitals which
c  were isolated in subroutine core:
c
      crthrs = crtset
      if(ispin.ne.0) crthrs = crthrs - one
      first = .true.
      do 952 iat = 1,natoms
        ilow = 0
        do 951 i = 1,nbas
          if(label(ibxm(i),1).eq.lcr.and.label(ibxm(i),4).eq.iat
     +       .and.bndocc(i).lt.crthrs) ilow = ilow + 1
  951   continue
        if(ilow.ne.0) then
          if(first) then
            first = .false.
            nam = nameat(iatno(iat))
            if(ilow.ne.1) then
              if(jprint(5).eq.1) write(lfnpr,3010) ilow,crthrs,nam,iat
            else
              if(jprint(5).eq.1) write(lfnpr,3011) ilow,crthrs,nam,iat
            end if
          else
            nam = nameat(iatno(iat))
            if(ilow.ne.1) then
              if(jprint(5).eq.1) write(lfnpr,3020) ilow,crthrs,nam,iat
            else
              if(jprint(5).eq.1) write(lfnpr,3021) ilow,crthrs,nam,iat
            end if
          end if
        end if
  952 continue
      return
c
c  problems with a bond orbital occupancy:
c
  960 write(lfnpr,1300) occi
      jprint(1) = -1
      return
c
c  total number of electrons is not an integer:
c
  970 write(lfnpr,1310) totele
      jprint(1) = -1
      return
c
 1200 format(/,1x,'for an occupancy threshold of ',f4.2,' the search',
     + ' for nbos found',/,1x,i3,' orbitals orbitals rather than ',i4)
 1210 format(3x,'label ',i3,':',a3,a1,i2,3i3)
 1300 format(/,1x,'a bond orbital with an occupancy of ',f8.5,
     + ' electrons was found!',/,1x,'please check you input data.')
 1310 format(/,1x,'the total number of electron is not an integer:',
     + f10.5,/,1x,'please check your input data.')
 1400 format(/,1x,'search of dm block between the following atoms:',
     +          3i4)
 1403 format(6x,'select orbitals with eigenvalue > ',f9.6)
 1405 format(6x,8f9.6)
 1410 format(6x,'eigenvector (',i2,') has occupancy ',f9.6,':')
 1420 format(11x,8f7.4)
 1500 format(11x,'hybrid ',i1,' in eigenvector has a projection ',
     +    'expectation of ',f6.3)
 1600 format(11x,'*** nbo accepted: number',i3,'.   label:',a2,a1,
     + '(',i2,')')
 1800 format(/4x,'prjthr will be raised from ',f6.3,' to',f6.3,
     + ' and the nbo search repeated.',/)
 1810 format(//,1x,'linearly independent hybrids for atom',i3,
     +' cannot be found.',/,1x,'the nbo program must abort.')
 3010 format(/,1x,
     +'WARNING:',i3,' low occupancy (<',f6.4,'e) core orbitals ',
     +'found on ',a2,i2)
 3011 format(/,1x,
     +'WARNING:',i3,' low occupancy (<',f6.4,'e) core orbital  ',
     +'found on ',a2,i2)
 3020 format(1x,
     +'        ',i3,' low occupancy (<',f6.4,'e) core orbitals ',
     +'found on ',a2,i2)
 3021 format(1x,
     +'        ',i3,' low occupancy (<',f6.4,'e) core orbital  ',
     +'found on ',a2,i2)
      end
c*****************************************************************************
      subroutine chsdrv(dm,t,guide,bndocc,pol,q,v,blk,c,eval,borb,
     *                                       p,ta,hyb,va,vb,topo)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      logical end,error,equal,equal2
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbnao/naoctr(maxbas),naol(maxbas),ltyp(maxbas),
     +       iprin(maxbas)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbtopo/iorder(maxatm),jorder(maxatm),ntopo(maxatm,maxatm),
     +            n3ctr,i3ctr(10,3)
      dimension dm(ndim,ndim),t(ndim,ndim),guide(natoms,natoms),
     * bndocc(ndim),pol(ndim,3),q(mxao,ndim),v(ndim),blk(mxbo,mxbo),
     * c(mxbo,mxbo),eval(mxbo),borb(mxbo),p(mxao,mxao),ta(mxao,mxao),
     * hyb(mxao),va(mxao),vb(mxao),topo(natoms,natoms)
      dimension keywd(6),klone(4),kbond(4),k3cbon(6),kalpha(5),
     * kbeta(4),ival(4),kalt(4)
      data klone/1hl,1ho,1hn,1he/,
     *     kbond/1hb,1ho,1hn,1hd/,
     *     k3cbon/1h3,1hc,1hb,1ho,1hn,1hd/,
     *     kalpha/1ha,1hl,1hp,1hh,1ha/,
     *     kbeta/1hb,1he,1ht,1ha/,
     *     ks/1hs/,kd/1hd/,kt/1ht/,kq/1hq/,
     *     kalt/1h$,1he,1hn,1hd/
c
c  search for alpha or beta character string in case of alpha or
c  beta spin density matrices:
c
      if(ispin.eq.2) then
   20   leng = 5
          call hfld(keywd,leng,end)
          if(end.and.leng.eq.0) goto 810
          if(.not.equal(keywd,kalpha,5)) goto 20
        continue
      else if(ispin.eq.-2) then
   30   leng = 5
          call hfld(keywd,leng,end)
          if(end.and.leng.eq.0) goto 820
          if(.not.equal(keywd,kbeta,4)) goto 30
        continue
      end if
c
c  fill diagonal elements of the topo matrix with nominal numbers of
c  lone pairs to be found on each atom:
c
      do 50 iat = 1,natoms
        nlp = 0
        call valtbl(iat,ival)
        do 40 l = 0,3
          nlp = nlp + ival(l+1)*(2*l + 1)
   40   continue
        ntopo(iat,iat) = 100 + nlp
   50 continue
c
c  read in chosen lone pairs, bonds, and 3-center bonds:
c
      nctr = 0
      n3ctr = 0
   60 continue
        leng = 6
        call hfld(keywd,leng,end)
        if(end.or.equal2(keywd,kalt,4)) goto 300
        nctro = nctr
        nctr = 0
        if(equal(keywd,klone,4))  nctr = 1
        if(equal(keywd,kbond,4))  nctr = 2
        if(equal(keywd,k3cbon,6)) nctr = 3
        if(nctr.eq.0) go to 1010
        if(nctr.lt.nctro) go to 1020
        goto (100,150,200), nctr
c
c  read in lone pairs:
c
  100 continue
        call ifld(iat,error)
        if(error) then
          leng = 6
          call hfld(keywd,leng,end)
          go to 60
        end if
        call ifld(num,error)
        if(error) goto 830
        ntopo(iat,iat) = num
      goto 100
c
c  read in bonds:
c
  150 continue
        leng = 1
        call hfld(keywd,leng,end)
        if(end) goto 60
        num = 0
        if(equal(keywd,ks,1)) num = 1
        if(equal(keywd,kd,1)) num = 2
        if(equal(keywd,kt,1)) num = 3
        if(equal(keywd,kq,1)) num = 4
        if(num.eq.0) goto 840
        call ifld(iat1,error)
        if(error) goto 840
        call ifld(iat2,error)
        if(error) goto 840
        iat = max(iat1,iat2)
        jat = min(iat1,iat2)
        ntopo(iat,jat) = num
        ntopo(jat,iat) = num
      goto 150
c
c  read in 3-center bonds:
c
  200 continue
        if(iw3c.ne.1) iw3c = 1
        leng = 1
        call hfld(keywd,leng,end)
        if(end) goto 60
        num = 0
        if(equal(keywd,ks,1)) num = 1
        if(equal(keywd,kd,1)) num = 2
        if(equal(keywd,kt,1)) num = 3
        if(equal(keywd,kq,1)) num = 4
        if(num.eq.0) goto 860
        call ifld(iat1,error)
        if(error) goto 860
        call ifld(iat2,error)
        if(error) goto 860
        call ifld(iat3,error)
        if(error) goto 860
        n3ctr = n3ctr + 1
        if(n3ctr.gt.10) goto 870
        i3ctr(n3ctr,1) = iat1
        i3ctr(n3ctr,2) = iat2
        i3ctr(n3ctr,3) = iat3
      goto 200
c
c  modify nominal sets of lone pairs by number of bonds and 3-center
c  bonds.
c
  300 continue
      do 330 iat = 1,natoms
        nlp = ntopo(iat,iat)
        if(nlp.lt.100) goto 330
        nlp = mod(nlp,100)
        nbd = 0
        do 310 jat = 1,natoms
          if(iat.ne.jat.and.ntopo(jat,iat).ne.0) then
            nbd = nbd + ntopo(jat,iat)
          end if
  310   continue
        do 320 kat = 1,3
          do 315 jat = 1,n3ctr
            if(i3ctr(jat,kat).eq.iat) nbd = nbd + 1
  315     continue
  320   continue
        nlp = nlp - nbd
        if(nlp.lt.0) nlp = 0
        ntopo(iat,iat) = nlp
  330 continue
c
c  use choose to find bond orbitals using ntopo and i3ctr:
c
      iflg = 0
      call choose(dm,t,guide,bndocc,pol,q,v,blk,c,eval,borb,p,ta,hyb,
     +            va,vb,topo,iflg)
      return
c
  810 write(lfnpr,1180)
      jprint(1) = -1
      return
c
  820 write(lfnpr,1190)
      jprint(1) = -1
      return
  830 write(lfnpr,1130)
      jprint(1) = -1
      return
c
  840 write(lfnpr,1140)
      jprint(1) = -1
      return
c
  860 write(lfnpr,1160)
      jprint(1) = -1
      return
c
  870 write(lfnpr,1170)
      jprint(1) = -1
      return
c
 1010 write(lfnpr,1110) (keywd(i),i=1,6)
      jprint(1) = -1
      return
c
 1020 write(lfnpr,1120)
      jprint(1) = -1
      return
c
 1110 format(/1x,'error in input of bond orbitals:',/,1x,
     * 'keyword for orbital type is not lone, bond, or 3cbond (read `',
     * 6a1,''')')
 1120 format(/1x,'error in input of bond orbitals:',/,1x,
     * 'orbital types should be in the order: lone, bond, 3cbond')
 1130 format(/1x,'error in input of bond orbitals:',/,1x,
     * 'unrecognizable characters in input of lone orbitals')
 1140 format(/1x,'error in input of bond orbitals:',/,1x,
     * 'unrecognizable characters in input of two center orbitals')
 1160 format(/1x,'error in input of bond orbitals:',/,1x,
     * 'unrecognizable characters in input of three center orbitals')
 1170 format(/1x,'too many three center bonds:',
     * '  increase parameter max3c')
 1180 format(/1x,'end of file encountered before the word alpha was ',
     * 'found')
 1190 format(/1x,'end of file encountered before the word beta was ',
     * 'found')
      end
c*****************************************************************************
c
c     subroutines called by nathyb and chsdrv for forming nbos
c
c
c*****************************************************************************
      subroutine choose(dm,t,guide,bndocc,pol,q,v,blk,c,eval,borb,
     *                                    p,ta,hyb,va,vb,topo,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  construct orthogonal matrix t for transformation from ao's to
c  natural hybrid bond orbitals using input density matrix dm
c  with the chosen bonding pattern read from lfnin
c
      logical detail,first,print,left
      integer ul
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       ul(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbnao/naoctr(maxbas),naol(maxbas),ltyp(maxbas),
     +       iprin(maxbas)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nbtopo/iorder(maxatm),jorder(maxatm),ntopo(maxatm,maxatm),
     +            n3ctr,i3ctr(10,3)
c
      dimension dm(ndim,ndim),t(ndim,ndim),guide(natoms,natoms),
     * bndocc(ndim),pol(ndim,3),q(mxao,ndim),v(ndim),blk(mxbo,mxbo),
     * c(mxbo,mxbo),eval(mxbo),borb(mxbo),p(mxao,mxao),ta(mxao,mxao),
     * hyb(mxao),va(mxao),vb(mxao),topo(natoms,natoms)
      dimension name(3),hybexp(3),ktopo(maxatm,maxatm),kflg(10)
      dimension scr(maxatm*(maxatm-1)/2),ipt(maxatm*(maxatm-1)/2)
c
      data istar,iblnk/1h*,1h /
      data name/2hlp,2hbd,2h3c/
      data lry,lcr/2hry,2hcr/
      data zero,zerop,tenth,pt99,one,two,twop
     + /0.0d0,1.0d-6,0.1d0,0.99d0,1.0d0,2.0d0,2.0001d0/
c
c  iflg is a print flag on entering choose.  if set to 0(1), choose
c  will(not) print some output to lfnpr.  on exit, if iflg is set to
c  -1, there was an error in finding the requested structure:
c
c  prjinc, the amount to increase prjthr by if problems with linear
c  dependency between the hybrids arise.
c
      data prjinc/0.05d0/
c
      nopval(i) = norbs(i) - ino(i)
c
      print = .false.
      if(iflg.eq.0) print = .true.
      if(jprint(5).eq.0) print = .false.
      detail = .false.
      if(iwdetl.ne.0) detail = .true.
      prjthr = abs(prjset)
      iter = 0
c
c  initialize ktopo and kflg arrays:  (kflg is set to 1 if the 3-center
c  bond has not been fund yet.)
c
      do 10 i = 1,natoms
        do 5 j = 1,i
          ktopo(i,j) = ntopo(i,j)
          ktopo(j,i) = ntopo(j,i)
    5   continue
   10 continue
      do 15 i = 1,n3ctr
        kflg(i) = 1
   15 continue
c
c  determine the atom ordering for the search for bond orbitals:
c
      if(natoms.eq.1) then
        iorder(1) = 1
      else
        ii = 0
        do 20 jat = 2,natoms
          do 19 iat = 1,jat-1
            ii = ii + 1
            scr(ii) = ktopo(iat,jat) - guide(iat,jat)
   19     continue
   20   continue
        nn = natoms * (natoms - 1) / 2
        call rank(scr,nn,nn,ipt)
c
c  begin search for bond orbitals where the formal bond order is much
c  greater than the corresponding wiberg bond index:
c
        ipos = 0
        jpos = 0
   21   continue
          jpos = jpos + 1
          if(jpos.gt.nn) stop 'problems with atom permutation list'
          iat  = ipt(jpos)
          jat  = 2
   22     continue
            if(jat.gt.iat) goto 23
            iat = iat - jat + 1
            jat = jat + 1
          goto 22
   23     continue
c
c  add iat and jat to the atom permutation list iorder:
c
          mflg = 0
          do 24 i = 1,ipos
            if(iorder(i).eq.iat) mflg = 1
   24     continue
          if(mflg.eq.0) then
            ipos = ipos + 1
            iorder(ipos) = iat
          end if
          mflg = 0
          do 25 i = 1,ipos
            if(iorder(i).eq.jat) mflg = 1
   25     continue
          if(mflg.eq.0) then
            ipos = ipos + 1
            iorder(ipos) = jat
          end if
        if(ipos.lt.natoms) goto 21
      end if
c
c  return to here if it should prove necessary to raise prjthr:
c
   35 continue
      iter = iter + 1
      occthr = abs(thrset)
      if(ispin.ne.0) occthr = occthr - one
      occthr = occthr + tenth
c
c  store density matrix in upper triangle of t:
c
      do 50 j = 1,nbas
        do 40 i = 1,j
          t(i,j) = dm(i,j)
   40   continue
   50 continue
c
c  zero arrays q,pol,iathy,ino,label:
c
      do 100 i = 1,nbas
        do 60 k = 1,2
          label(i,k) = iblnk
   60   continue
        do 70 k = 3,6
          label(i,k) = 0
   70   continue
        do 80 k = 1,3
          pol(i,k) = zero
          iathy(i,k) = 0
   80   continue
        do 90 k = 1,mxao
          q(k,i) = zero
   90   continue
  100 continue
      do 110 i = 1,natoms
        ino(i) = 0
  110 continue
c
c  remove core orbitals from the density matrix:
c
      ibd = 0
      call core(dm,t,borb,pol,q,hyb,bndocc,ibd,detail,lfnpr)
c
c  return here if there are still more lone pairs or bonds to be found.
c  lower the occupancy threshold for acceptance by a tenth:
c
  115 continue
      occthr = occthr - tenth
      left = .false.
c
c   ********      start directed nbo search      *********
c
c  loop over numbers of centers, removing lone pairs and 2- and 3-center
c  bonds from the density matrix according to ktopo and i3ctr:
c
      nctr = 0
  120 nctr = nctr + 1
c
c  deplete dm of one(two) center orbitals if search for two(three)
c  center orbitals is beginning:
c
        if(nctr.ne.1) call deplet(dm,t,q,pol,borb,bndocc,ibd)
c
        icntr = 0
c
c  return here for 3-c bonds and lone pairs:
c
  130   icntr = icntr + 1
          if(nctr.eq.1) then
            if(icntr.gt.natoms) goto 120
            num = ktopo(iorder(icntr),iorder(icntr))
            if(num.le.0) goto 130
            iat1 = iorder(icntr)
            iat2 = 0
            iat3 = 0
            goto 200
          else if(nctr.eq.2) then
            if(icntr.gt.natoms) goto 120
            jcntr = icntr
c
c  return here for 2-c bonds:
c
  150       jcntr = jcntr + 1
              if(jcntr.gt.natoms) goto 130
              num = ktopo(iorder(jcntr),iorder(icntr))
              if(num.eq.0) goto 150
              iat1 = min(iorder(icntr),iorder(jcntr))
              iat2 = max(iorder(icntr),iorder(jcntr))
              iat3 = 0
              goto 200
          else if(nctr.eq.3) then
            if(icntr.gt.n3ctr) goto 120
            if(kflg(icntr).eq.0) goto 130
            num = 1
            iat1 = min(i3ctr(icntr,1),i3ctr(icntr,2),i3ctr(icntr,3))
            iat3 = max(i3ctr(icntr,1),i3ctr(icntr,2),i3ctr(icntr,3))
            iat2 = i3ctr(icntr,1)
            if(iat2.eq.iat1.or.iat2.eq.iat3) iat2 = i3ctr(icntr,2)
            if(iat2.eq.iat1.or.iat2.eq.iat3) iat2 = i3ctr(icntr,3)
            goto 200
          else
            goto 300
          end if
c
c  load proper atomic blocks of dm into blk, and diagonalize blk:
c
  200 continue
      call NBOload(dm,iat1,iat2,iat3,blk,nb)
      call NBOjacobi(nb,blk,eval,c,mxbo,mxbo,1)
c
c  rank eigenvectors by occupancy eigenvalue:
c
      call rank(eval,nb,mxbo,larc)
      if(detail) write(lfnpr,1400) iat1,iat2,iat3
      if(detail) write(lfnpr,1402) num,occthr
      if(detail) write(lfnpr,1405) (eval(irnk),irnk=1,nb)
c
c  loop over eigenvalues selecting the num highest occupied:
c
      iaccep = 0
      do 250 irnk = 1,nb
        ir = larc(irnk)
        occ = eval(irnk)
        do 210 i = 1,nb
  210     borb(i) = c(i,ir)
        if(detail) write(lfnpr,1410) irnk,occ
        if(detail) write(lfnpr,1420) (borb(i),i=1,nb)
c
c  if this is a low occupancy orbital, skip the rest of these and can come
c  back to them later:
c
        if(occ.lt.occthr) then
          if(nctr.eq.1) then
            ktopo(iat1,iat1) = num - iaccep
            if(detail) write(lfnpr,1610) ktopo(iat1,iat1)
          else if(nctr.eq.2) then
            ktopo(iat1,iat2) = num - iaccep
            ktopo(iat2,iat1) = ktopo(iat1,iat2)
            if(detail) write(lfnpr,1610) ktopo(iat1,iat2)
          else
            ione = 1
            if(detail) write(lfnpr,1610) ione
          end if
          if(left) then
            if(occmax.lt.occ) occmax = occ
          else
            left = .true.
            occmax = occ
          end if
          goto 280
        end if
c
c  check to see if bond orbital "borb" contains previously used hybrids:
c
c        if(nctr.ne.1) then
          call prjexp(borb,iat1,iat2,iat3,q,p,ta,hyb,va,vb,hybexp)
          if(detail) then
            do 220 ihyb = 1,nctr
              write(lfnpr,1500) ihyb,hybexp(ihyb)
  220       continue
          end if
          do 230 ihyb = 1,nctr
            if(hybexp(ihyb).lt.prjthr) goto 250
  230     continue
c        end if
        ibd = ibd + 1
        iaccep = iaccep + 1
c
c  decompose "borb" into its constituent atomic hybrids and store in q:
c
        call stash_nbo(borb,ibd,iat1,iat2,iat3,pol,q,hyb)
c
c  construct bond orbital labels:
c
        if(nctr.eq.1) then
          ishift = ntopo(iat1,iat1) - ktopo(iat1,iat1)
        else if(nctr.eq.2) then
          ishift = ntopo(iat1,iat2) - ktopo(iat1,iat2)
        else
          ishift = 0
        end if
        label(ibd,1) = name(nctr)
        label(ibd,2) = iblnk
        label(ibd,3) = iaccep + ishift
        label(ibd,4) = iat1
        label(ibd,5) = iat2
        label(ibd,6) = iat3
        bndocc(ibd) = occ
        if(detail) write(lfnpr,1600) ibd,(label(ibd,i),i=1,3)
        if(iaccep.eq.num) then
          if(nctr.eq.1) then
            ktopo(iat1,iat1) = 0
          else if(nctr.eq.2) then
            ktopo(iat1,iat2) = 0
            ktopo(iat2,iat1) = 0
          else
            kflg(icntr) = 0
          end if
          goto 280
        end if
  250 continue
        if(iaccep.ne.num.and.nctr.eq.2.and.print)
     *            write(lfnpr,2000) prjthr,iaccep,num,iat1,iat2
        if(iaccep.ne.num.and.nctr.eq.3.and.print)
     *            write(lfnpr,2100) prjthr,iaccep,num,iat1,iat2,iat3
        iflg = -1
  280 if(nctr.eq.1.or.nctr.eq.3) then
        goto 130
      else
  290   jcntr=jcntr+1
        if(jcntr.gt.natoms) goto 130
        num=ktopo(iorder(jcntr),iorder(icntr))
        if(num.eq.0) goto 290
        iat1=iorder(icntr)
        iat2=iorder(jcntr)
        iat3=0
        goto 200
      end if
c
c   ******** end of loop for directed nbo search *********
c
  300 continue
c
c  if some orbitals were left behind, go back and fetch them:
c
      if(left) then
        occthr = occmax
        goto 115
      end if
c
c  symmetrically orthogonalize principal hybrids:
c
      call orthyb(q,blk,ta,eval,c,ialarm,iflg)
c
c  ialarm sounds the alarm that there is linear dependency between some of the
c  hybrids.  ialarm is equal to the number of the violating atom.  replenish
c  dm from t and repeat the nbo search:
c
      if(ialarm.ne.0) then
        oldprj = prjthr
        prjthr = oldprj + prjinc
        if(print) write(lfnpr,1800) oldprj,prjthr
        if(prjthr.ge.pt99) then
          if(print) write(lfnpr,1810) ialarm
          iflg = -1
          jprint(1) = -1
          return
        end if
        goto 700
      end if
c
c  augment open-valence atoms with non-arbitrary hybrids orthogonal to those
c  found previously:
c
      do 580 ia = 1,natoms
        if(nopval(ia).le.0) goto 580
c
c  find nmb, the number of natural minimal basis functions on this atom:
c
        lla = ll(ia)
        iula = ul(ia)
        nmb = 0
        do 470 i = lla,iula
          if(lstocc(i).eq.1) nmb = nmb + 1
  470   continue
c
c  find the number of bond, core, and lone pair hybrids on this atom, iocc.
c  also find iocclp, the number of lone pair orbitals already found
c  on atom ia for use in labelling the extra lone pairs below:
c
        iocc = 0
        iocclp = 0
        do 480 ib = 1,ibd
          if((label(ib,4).ne.ia).and.(label(ib,5).ne.ia).and.
     *            (label(ib,6).ne.ia)) goto 480
          iocc = iocc + 1
          if(label(ib,1).eq.name(1)) then
            iocclp = iocclp + 1
          end if
  480   continue
c
c  nexlp, the number of extra (low occupancy) lp orbitals on atom iat.
c  (this is the number of low occupancy orbitals with valence shell character)
c  set nexlp to zero if (nmb-iocc) is less than zero!!
c
        nexlp = nmb - iocc
        if(nexlp.lt.0) nexlp = 0
        nocc = ino(ia)
        call frmprj(p,ia,q,nocc,ta,va,vb)
        norb = norbs(ia)
        naugm = norb - nocc
        call NBOaugmnt(p,blk,c,eval,dm,ta,borb,v,larc,ia,nocc,norb)
c
c  stash_nbo and label extra lone pairs that augmnt put in blk:
c  (these are taken to be the highest occupied orbitals, which
c  augmnt places first)
c
        do 510 iaugm = 1,nexlp
          do 500 j = 1,norb
  500       borb(j) = blk(j,iaugm)
          ibd = ibd + 1
          call stash_nbo(borb,ibd,ia,0,0,pol,q,hyb)
          label(ibd,1) = name(1)
          label(ibd,2) = istar
          label(ibd,3) = iaugm+iocclp
          label(ibd,4) = ia
          label(ibd,5) = 0
          label(ibd,6) = 0
  510   continue
c
c  stash_nbo and label the rydberg orbitals that augmnt put in blk:
c
        iryd = 0
        nstart = nexlp + 1
        do 540 iaugm = nstart,naugm
          do 530 j = 1,norb
  530       borb(j) = blk(j,iaugm)
          ibd = ibd + 1
          iryd = iryd + 1
          call stash_nbo(borb,ibd,ia,0,0,pol,q,hyb)
          label(ibd,1) = lry
          label(ibd,2) = istar
          label(ibd,3) = iryd
          label(ibd,4) = ia
          label(ibd,5) = 0
          label(ibd,6) = 0
  540   continue
  580 continue
c
c  include antibond labels:
c
      ibo = ibd
      do 660 i = 1,ibo
c
c  exit loop if label(i,1) is 'lp', 'ry', or 'cr':
c
        if(label(i,1).eq.name(1)) goto 660
        if(label(i,1).eq.lry)     goto 660
        if(label(i,1).eq.lcr)     goto 660
          nab = 1
          if(label(i,1).eq.name(3)) nab = 2
          do 650 iab = 1,nab
            ibd = ibd + 1
            do 640 j = 1,6
  640         label(ibd,j) = label(i,j)
            label(ibd,2) = istar
  650     continue
  660 continue
      if(ibd.eq.nbas) goto 670
        write(lfnpr,2200)
        stop
  670 continue
c
c  replace density matrix dm from t:
c
  700 continue
      do 750 j = 1,nbas
        do 740 i = 1,j
          dm(i,j) = t(i,j)
          dm(j,i) = dm(i,j)
          t(j,i) = zero
          t(i,j) = zero
  740   continue
  750 continue
c
c  if the alarm sounded, repeat directed nbo search:
c
      if(ialarm.ne.0) goto 35
c
c  find new polarization parameters for orthonormal hybrids:
c
      call repol(dm,q,pol,blk,eval,c,ibd)
c
c  form final t-nab (nao to nbo transformation) from orthonormal hybrids:
c
      call formt(t,q,pol)
c
c  find occupancies, find total number of electrons and occupied orbitals:
c
      totele = zero
      do 800 i = 1,nbas
        occi = zero
        do 790 j = 1,nbas
          do 790 k = 1,nbas
  790       occi = occi + t(j,i) * dm(j,k) * t(k,i)
        if(abs(occi).lt.zerop) occi = zero
        if(occi.gt.twop) go to 960
        zeropm = -zerop
        if(occi.lt.zeropm) go to 960
        bndocc(i) = occi
        v(i) = occi
        totele = totele + bndocc(i)
  800 continue
      nel = totele + tenth
      if(abs(totele-nel).gt.5e-4) go to 965
      totele = nel
      nocc = nel
      if(ispin.eq.0) nocc = nocc/2 + mod(nocc,2)
c
c  if the number of unstarred orbitals is not equal to the number of occupied
c  mos, then simply rank the orbitals by occupancy, and `unstarr' the nocc
c  highest occupied:  (this can be dangerous!  however, many of the subsequent
c  routines assume the only nocc orbitals are starred, and therefore, this
c  mismatch must be corrected.)
c
      nostr = 0
      do 801 i = 1,nbas
        if(label(ibxm(i),2).ne.istar) nostr = nostr + 1
  801 continue
      if(nostr.ne.nocc) then
        call rank(v,nbas,ndim,larc)
        do 804 i = 1,nocc
          ir = larc(i)
          label(ibxm(ir),2) = iblnk
  804   continue
        do 805 i = nocc+1,nbas
          ir = larc(i)
          label(ibxm(ir),2) = istar
  805   continue
      end if
c
c  determine whether this is a good resonance structure:
c
      call cycles(iter,abs(thrset),guide,bndocc,topo,icont)
c
c  write out info about core orbitals which were isolated in subroutine
c  core:
c
      if(.not.print) goto 953
      crthrs = crtset
      if(ispin.ne.0) crthrs = crthrs - one
      first = .true.
      do 952 iat = 1,natoms
        ilow = 0
        do 951 i = 1,nbas
          if(label(ibxm(i),1).eq.lcr.and.label(ibxm(i),4).eq.iat
     +       .and.bndocc(i).lt.crthrs) ilow = ilow + 1
  951   continue
        if(ilow.ne.0) then
          if(first) then
            first = .false.
            nam = nameat(iatno(iat))
            if(ilow.ne.1) then
              write(lfnpr,3010) ilow,crthrs,nam,iat
            else
              write(lfnpr,3011) ilow,crthrs,nam,iat
            end if
          else
            nam = nameat(iatno(iat))
            if(ilow.ne.1) then
              write(lfnpr,3020) ilow,crthrs,nam,iat
            else
              write(lfnpr,3021) ilow,crthrs,nam,iat
            end if
          end if
        end if
  952 continue
  953 continue
      return
c
c  bad orbital occupancy:
c
  960 if(print) write(lfnpr,1300) occi
      iflg = -1
      jprint(1) = -1
      return
c
c  total number of electrons is not an integer:
c
  965 write(lfnpr,1310) totele
      iflg = -1
      jprint(1) = -1
      return
c
 1300 format(/,1x,'a bond orbital with an occupancy of ',f8.5,
     + ' electrons was found!',/,1x,'please check you input data.')
 1310 format(/,1x,'the total number of electron is not an integer:',
     + f10.5,/,1x,'please check your input data.')
 1400 format(/,1x,'search of dm block between the following atoms:',
     +          3i4)
 1402 format(6x,'select ',i2,' orbital(s) with eigenvalue > ',f9.6)
 1405 format(6x,8f9.6)
 1410 format(6x,'eigenvector (',i2,') has occupancy ',f9.6,':')
 1420 format(11x,8f7.4)
 1500 format(11x,'hybrid ',i1,' in eigenvector has a projection ',
     +    'expectation of ',f6.3)
 1600 format(11x,'*** nbo accepted: number',i3,'.   label:',a2,a1,
     + '(',i2,')')
 1610 format(1x,'still need to find',i2,' more orbital(s)')
 1800 format(/4x,'prjthr will be raised from ',f6.3,' to',f6.3,
     + ' and the nbo search repeated.',/)
 1810 format(//,1x,'linearly independent hybrids for atom',i3,
     +' cannot be found.',/,1x,'the nbo program must abort.')
 2000 format(/,1x,'at a projection threshold of',f6.3,', only ',i1,
     + ' of the ',i1,' requested bonds',/,1x,'between atoms ',i2,
     + ' and ',i2,' can be constructed.  the nbo analysis will',/,
     + 1x,'continue, augmenting the nbo set with extra lone pairs ',
     + 'on the atoms',/,1x,'as necessary.')
 2100 format(/,1x,'at a projection threshold of',f6.3,', only ',i1,
     + ' of the ',i1,' requested bonds',/,1x,'between atoms ',i2,', ',
     + i2,', and ',i2,' can be constructed.  the nbo analysis',/,1x,
     + 'will continue, augmenting the nbo set with extra lone pairs ',
     + 'on the',/,1x,'atoms as necessary.')
 2200 format(/,1x,'miscounted orbitals, program must abort')
 3010 format(/,1x,
     +'WARNING:',i3,' low occupancy (<',f6.4,'e) core orbitals ',
     +'found on ',a2,i2)
 3011 format(/,1x,
     +'WARNING:',i3,' low occupancy (<',f6.4,'e) core orbital  ',
     +'found on ',a2,i2)
 3020 format(1x,
     +'        ',i3,' low occupancy (<',f6.4,'e) core orbitals ',
     +'found on ',a2,i2)
 3021 format(1x,
     +'        ',i3,' low occupancy (<',f6.4,'e) core orbital  ',
     +'found on ',a2,i2)
      end
c*****************************************************************************
      subroutine srtnbo(t,bndocc)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical permut
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
c
      dimension t(ndim,ndim),bndocc(ndim)
      dimension name(3)
c
      data lbd,l3c,name,lstar/2hbd,2h3c,2hcr,2hlp,2hry,1h*/
c
c  reorder the nbos according to bond type and constituent atomic centers:
c
c  fix atom ordering in the nbo labels:
c
      do 100 i = 1,nbas
        nctr = 0
        do 10 j = 4,6
          if(label(i,j).ne.0) then
            nctr = nctr + 1
            larc(nctr) = label(i,j)
          end if
   10   continue
        do 30 j = 1,nctr-1
          do 20 k = 1,nctr-j
            if(larc(k).gt.larc(k+1)) then
              itemp     = larc(k)
              larc(k)   = larc(k+1)
              larc(k+1) = itemp
            end if
   20     continue
   30   continue
        do 40 j = 1,nctr
          label(i,j+3) = larc(j)
   40   continue
        do 50 j = nctr+1,3
          label(i,j+3) = 0
   50   continue
  100 continue
c
c  place the 2- and 3-center bonds first in the list of nbos: (no bonds if
c  the nobond keyword was specified)
c
      icnt = 0
      if(jprint(10).eq.0) then
        do 200 i = 1,natoms-1
          do 190 j = i+1,natoms
            if(i.ne.j) then
              k = -1
  110         k = k + 1
              do 180 l = icnt+1,nbas
                lbl1 = label(ibxm(l),1)
                lbl2 = label(ibxm(l),2)
                lbl3 = label(ibxm(l),3)
                lbl4 = label(ibxm(l),4)
                lbl5 = label(ibxm(l),5)
                lbl6 = label(ibxm(l),6)
                if((lbl1.eq.lbd.or.lbl1.eq.l3c).and.lbl2.ne.lstar) then
                  if(lbl4.eq.i.and.lbl5.eq.j.and.lbl6.eq.k) then
                    icnt = icnt + 1
                    label(ibxm(l),1)    = label(ibxm(icnt),1)
                    label(ibxm(l),2)    = label(ibxm(icnt),2)
                    label(ibxm(l),3)    = label(ibxm(icnt),3)
                    label(ibxm(l),4)    = label(ibxm(icnt),4)
                    label(ibxm(l),5)    = label(ibxm(icnt),5)
                    label(ibxm(l),6)    = label(ibxm(icnt),6)
                    label(ibxm(icnt),1) = lbl1
                    label(ibxm(icnt),2) = lbl2
                    label(ibxm(icnt),3) = lbl3
                    label(ibxm(icnt),4) = lbl4
                    label(ibxm(icnt),5) = lbl5
                    label(ibxm(icnt),6) = lbl6
                    temp         = bndocc(l)
                    bndocc(l)    = bndocc(icnt)
                    bndocc(icnt) = temp
                    do 170 m = 1,nbas
                      temp      = t(m,l)
                      t(m,l)    = t(m,icnt)
                      t(m,icnt) = temp
  170               continue
                  end if
                end if
  180         continue
              if(iw3c.ne.0.and.k.eq.0) k = j
              if(k.gt.0.and.k.lt.natoms) goto 110
            end if
  190     continue
  200   continue
      end if
c
c  next add any core, lone pair, and rydberg orbitals:
c
      do 300 ii = 1,3
        do 290 i = 1,natoms
          do 280 j = icnt+1,nbas
            lbl1 = label(ibxm(j),1)
            lbl4 = label(ibxm(j),4)
            if(lbl1.eq.name(ii).and.lbl4.eq.i) then
              icnt = icnt + 1
              do 260 k = 1,6
                itemp               = label(ibxm(j),k)
                label(ibxm(j),k)    = label(ibxm(icnt),k)
                label(ibxm(icnt),k) = itemp
  260         continue
              temp         = bndocc(j)
              bndocc(j)    = bndocc(icnt)
              bndocc(icnt) = temp
              do 270 k = 1,nbas
                temp      = t(k,j)
                t(k,j)    = t(k,icnt)
                t(k,icnt) = temp
  270         continue
            end if
  280     continue
  290   continue
  300 continue
c
c  add in any antibonds:
c
      if(jprint(10).eq.0) then
        do 400 i = 1,natoms-1
          do 390 j = i+1,natoms
            if(i.ne.j) then
              k = -1
              if(iw3c.ne.0) k = j
  310         k = k + 1
              do 380 l = icnt+1,nbas
                lbl1 = label(ibxm(l),1)
                lbl2 = label(ibxm(l),2)
                lbl3 = label(ibxm(l),3)
                lbl4 = label(ibxm(l),4)
                lbl5 = label(ibxm(l),5)
                lbl6 = label(ibxm(l),6)
                if((lbl1.eq.lbd.or.lbl1.eq.l3c).and.lbl2.eq.lstar) then
                  if(lbl4.eq.i.and.lbl5.eq.j.and.lbl6.eq.k) then
                    icnt = icnt + 1
                    label(ibxm(l),1)    = label(ibxm(icnt),1)
                    label(ibxm(l),2)    = label(ibxm(icnt),2)
                    label(ibxm(l),3)    = label(ibxm(icnt),3)
                    label(ibxm(l),4)    = label(ibxm(icnt),4)
                    label(ibxm(l),5)    = label(ibxm(icnt),5)
                    label(ibxm(l),6)    = label(ibxm(icnt),6)
                    label(ibxm(icnt),1) = lbl1
                    label(ibxm(icnt),2) = lbl2
                    label(ibxm(icnt),3) = lbl3
                    label(ibxm(icnt),4) = lbl4
                    label(ibxm(icnt),5) = lbl5
                    label(ibxm(icnt),6) = lbl6
                    temp         = bndocc(l)
                    bndocc(l)    = bndocc(icnt)
                    bndocc(icnt) = temp
                    do 370 m = 1,nbas
                      temp      = t(m,l)
                      t(m,l)    = t(m,icnt)
                      t(m,icnt) = temp
  370               continue
                  end if
                end if
  380         continue
              if(k.gt.0.and.k.lt.natoms) goto 310
            end if
  390     continue
  400   continue
      end if
c
c  lastly, make sure orbitals are ordered by serial number:
c
  410 permut = .false.
      do 500 i = 1,nbas-1
        if(label(ibxm(i),1).eq.label(ibxm(i+1),1)) then
          if(label(ibxm(i),2).eq.label(ibxm(i+1),2)) then
            if(label(ibxm(i),4).eq.label(ibxm(i+1),4)) then
              if(label(ibxm(i),5).eq.label(ibxm(i+1),5)) then
                if(label(ibxm(i),6).eq.label(ibxm(i+1),6)) then
                  if(label(ibxm(i),3).gt.label(ibxm(i+1),3)) then
                    permut = .true.
                    lbl3 = label(ibxm(i),3)
                    label(ibxm(i),3) = label(ibxm(i+1),3)
                    label(ibxm(i+1),3) = lbl3
                    temp = bndocc(i)
                    bndocc(i) = bndocc(i+1)
                    bndocc(i+1) = temp
                    do 490 j = 1,nbas
                      temp = t(j,i)
                      t(j,i) = t(j,i+1)
                      t(j,i+1) = temp
  490               continue
                  end if
                end if
              end if
            end if
          end if
        end if
  500 continue
      if(permut) goto 410
      return
      end
c*****************************************************************************
      subroutine xcited(dm,t,hyb,thyb,s,occ,scr,iscr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical first
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbnao/naoc(maxbas),naoa(maxbas),ltyp1(maxbas),
     +       iprin(maxbas)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),ltyp(maxbas),iathy(maxbas,3)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      dimension dm(ndim,ndim),t(ndim,ndim),hyb(mxao),thyb(ndim,ndim),
     +          s(ndim,ndim),occ(ndim),scr(ndim),iscr(ndim)
      dimension pct(5),iat(3)
      data llp,lbd,l3c,lcr,lry/2hlp,2hbd,2h3c,2hcr,2hry/
      data zero,tenth,one,thresh/0.0d0,0.1d0,1.0d0,1.0d-4/
      data lstar,lblnk/1h*,1h /
c
c  form a temporary nao to nho transformation matrix.  check hybrid
c  overlap to make sure the nbo's were properly labelled as lewis
c  and non-lewis orbitals:
c
c  count number of hybrids as they are written out:
c
      nhyb = 0
c
c  main loop over bond orbitals:
c
      do 200 nbond = 1,nbas
        ib = ibxm(nbond)
        lbl = label(ib,1)
        if(lbl.eq.llp.or.lbl.eq.lcr.or.lbl.eq.lry) nctr = 1
        if(lbl.eq.lbd) nctr = 2
        if(lbl.eq.l3c) nctr = 3
c
c  loop over atomic centers of bond orbital nbond:
c
        do 190 ictr = 1,nctr
          i = label(ib,ictr+3)
          kl = ll(i)
          ku = lu(i)
          do 120 k = 1,mxao
            ltyp(k) = 0
  120       hyb(k) = zero
c
c  choose sign for polarization coefficients:
c
          isgn = 1
          if(label(ib,2).ne.lstar) go to 130
          if(ictr.lt.2) go to 130
          if(ictr.eq.3) ipar3c = -ipar3c
          if(ictr.eq.3.and.ipar3c.gt.0) go to 130
          isgn = -isgn
  130     continue
c
c  extract hybrid (hyb) from transformation matrix t; ltyp(i) is the
c  orbital angular momentum quantum no. of a.o. # i:
c
          kh = 0
          do 140 k = kl,ku
            kh = kh + 1
            hyb(kh) = t(k,nbond)
  140       ltyp(kh) = naoa(k)/100
          call htype(hyb,ltyp,mxao,kh,coef,pct,nl,isgn)
          if(abs(coef).lt.thresh) go to 190
c
c  check to see if this orbital has been found before:
c
          do 160 ihyb = 1,nhyb
            temp = zero
            ih = 0
            do 150 k = kl,ku
              ih = ih + 1
              temp = temp + hyb(ih)*thyb(k,ihyb)
  150       continue
            if(abs(abs(temp)-one).lt.thresh) go to 190
            if(abs(temp).gt.thresh) then
              write(lfnpr,900) nhyb+1,nbond,ictr,temp,ihyb
              stop
            end if
  160     continue
c
c  add this hybrid to the temporary thyb:
c
          nhyb = nhyb + 1
          if(nhyb.gt.nbas) stop 'too many hybrids'
          do 170 k = 1,nbas
            thyb(k,nhyb) = zero
  170     continue
          ih = 0
          do 180 k = kl,ku
            ih = ih + 1
            thyb(k,nhyb) = hyb(ih)
  180     continue
  190   continue
  200 continue
      if(nhyb.lt.nbas) stop 'missing hybrids'
c
c  thyb now contains the temporary nao to nho transformation matrix.
c  form the non-orthogonal pnho overlap and nho to nbo transformation matrices:
c
      call fesnao(s)
      call simtrs(s,thyb,scr,ndim,nbas)
c
      call NBOtransp(thyb,ndim,nbas)
      call matmlt(thyb,t,scr,ndim,nbas)
c
c  check to see that the bonds and antibonds have the correct hybrid
c  overlap.  fix the labels if there is a problem:
c
      first = .true.
      do 300 nbond = 1,nbas
        ib = ibxm(nbond)
        lbl1 = label(ib,1)
        if(lbl1.eq.llp.or.lbl1.eq.lcr.or.lbl1.eq.lry) ictr = 1
        if(lbl1.eq.lbd) ictr = 2
        if(lbl1.eq.l3c) ictr = 3
        nctr = 0
        do 210 ihyb = 1,nhyb
          if(abs(thyb(ihyb,nbond)).gt.thresh) then
            nctr = nctr + 1
            if(nctr.gt.3) then
              write(lfnpr,910) nbond
              stop
            end if
            iat(nctr) = ihyb
          end if
  210   continue
        if(nctr.gt.ictr) then
          write(lfnpr,920) ictr,nbond,nctr
          stop
        end if
        if(nctr.gt.1) then
          isgn = 1
          do 230 jctr = 1,nctr-1
            do 220 kctr = jctr+1,nctr
              jhyb = iat(jctr)
              khyb = iat(kctr)
              temp = s(jhyb,khyb)*thyb(jhyb,nbond)*thyb(khyb,nbond)
              if(temp.lt.zero) isgn = -1
  220       continue
  230     continue
          lbl2 = label(ib,2)
          if(lbl2.eq.lblnk.and.isgn.eq.-1) then
            if(first.and.jprint(5).ne.0) write(lfnpr,930)
            first = .false.
            label(ib,2) = lstar
            if(jprint(5).ne.0) write(lfnpr,940) nbond,lbl1,lstar
          else if(lbl2.eq.lstar.and.isgn.eq.1) then
            if(first.and.jprint(5).ne.0) write(lfnpr,930)
            first = .false.
            label(ib,2) = lblnk
            if(jprint(5).ne.0) write(lfnpr,940) nbond,lbl1,lblnk
          end if
        end if
  300 continue
c
c  determine the number of occupied orbitals:
c
      tot = zero
      do 310 i = 1,nbas
        tot = tot + dm(i,i)
  310 continue
      nocc = tot + tenth
      if(ispin.eq.0) nocc = nocc/2 + mod(nocc,2)
c
c  count the number of unstarred orbitals:
c
      icnt = 0
      do 320 i = 1,nbas
        if(label(ibxm(i),2).ne.lstar) icnt = icnt + 1
  320 continue
c
c  if the number of unstarred orbitals is not equal to the number of
c  occupied orbitals, fix the orbital labels:
c
      if(icnt.ne.nocc) then
        do 330 i = 1,nbas
          occ(i) = dm(i,i)
  330   continue
        call rank(occ,nbas,ndim,iscr)
c
c  if there are more unstarred orbitals than occupied, add stars to the
c  least occupied lone pairs:
c
        if(icnt.gt.nocc) then
          idiff = icnt - nocc
          do 350 i = 1,idiff
            ip = 0
            do 340 j = 1,nbas
              jp = ibxm(iscr(j))
              if(label(jp,1).eq.llp.and.label(jp,2).ne.lstar) ip = j
  340       continue
            if(ip.eq.0) then
              write(lfnpr,950) icnt,nocc
              stop
            end if
            label(ibxm(iscr(ip)),2) = lstar
            if(jprint(5).ne.0) write(lfnpr,940) iscr(ip),
     +                         label(ibxm(iscr(ip)),1),lstar
  350     continue
c
c  remove stars from the highest occupied lone pairs/rydbergs if there are
c  too few starred orbitals:
c
        else
          idiff = nocc - icnt
          do 370 i = 1,idiff
            ip = 0
            do 360 j = nbas,1,-1
              jp = ibxm(iscr(j))
              if((label(jp,1).eq.llp.or.label(jp,1).eq.lry)
     +                         .and.label(jp,2).eq.lstar) ip = j
  360       continue
            if(ip.eq.0) then
              write(lfnpr,950) icnt,nocc
              stop
            end if
            label(ibxm(iscr(ip)),2) = lblnk
            if(jprint(5).ne.0) write(lfnpr,940) iscr(ip),
     +                         label(ibxm(iscr(ip)),1),lblnk
  370     continue
        end if
      end if
      return
c
  900 format(/1x,'hybrid ',i3,' (nbo ',i3,', center ',i2,') has a ',
     + 'non-negligible overlap of ',f8.5,/,1x,'with hybrid ',i3,'.')
  910 format(/1x,'nbo ',i3,' has hybrid contributions from more than ',
     + '3 atomic centers.')
  920 format(/1x,'error: the ',i1,'-center nbo ',i3,' has ',
     + 'contributions from ',i2,' atomic centers.')
  930 format(/1x,'          --- apparent excited state configuration ',
     + '---',/1x,'the following "inverted" nbo labels reflect the ',
     + 'actual hybrid overlap:')
  940 format(1x,'                nbo ',i3,' has been relabelled ',a2,a1)
  950 format(/1x,'unable to label the nbos properly: ',i3,' unstarred ',
     + 'orbitals',/1x,'                                   ',i3,
     + ' occupied orbitals')
      end
c*****************************************************************************
      subroutine anlyze(t,bndocc,hyb,hycoef,thyb)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      integer ul
c
c  print out details of bond-orbital transformation from matrix t.
c
c  required input:
c         t = transformation matrix from s.r. nathyb; real (1,ndim;1,ndim)
c      ndim = declared dimensionality of array t
c      nbas = no. of orbitals = actual dimension of t, naol
c      naol = integer list of orbital angular momentum type
c                naol(i)/100 = l = q.n. of atomic orbital i
c     iatno = list of atomic numbers; iatno(i) is the atomic number
c                of atom i as an integer
c    natoms = no. of atoms (not including ghosts) in the molecule
c    iwhybs = 1 if hybrid a.o. coefficients are to be printed,
c             0 otherwise
c     lfnpr = logical file number for printout.
c    naoctr = list of atomic centers of oao or nao basis orbitals
c     label = list of bond orbital labels
c      ibxm = permutation list of bond orbitals
c    bndocc = list of bond orbital occupancies
c     ispin = 0 for spinless nbos
c           = 2 for alpha spin nbos
c           =-2 for beta  spin nbos
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       ul(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbnao/naoc(maxbas),naoa(maxbas),ltyp1(maxbas),
     +       iprin(maxbas)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),ltyp(maxbas),iathy(maxbas,3)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      dimension t(ndim,ndim),hyb(mxao),bndocc(ndim),thyb(ndim,ndim),
     * pct(5),pow(5),lname(5),isp(3),nam(3),ich(3,2),hycoef(ndim)
      data llp,lbd,l3c,lcr,lry/2hlp,2hbd,2h3c,2hcr,2hry/
      data lname/1hs,1hp,1hd,1hf,1hg/
      data zero,thresh,t99,t99p/0.0d0,1.d-2,99.99d0,99.995d0/
      data tenth,hundrd,tthoth/0.1d0,100.0d0,0.0001d0/
      data lhyp,lblnk,lstar,l2blnk/1h-,1h ,1h*,2h  /
c
c  count the number of electrons:
c
      totele = zero
      do 20 i = 1,nbas
        totele = totele + bndocc(i)
   20   continue
      totele = totele + tenth
      nel = totele
      totele = nel
c
c  count the number of core orbitals and the occupancies of the core,
c  valence lewis, valence non-lewis, and extra-valence rydberg orbitals.
c  (also count the number of electrons in the ecp, if employed)
c
      mcr = 0
      occcr = zero
      occvl = zero
      occvnl = zero
      do 50 i = 1,nbas
        if(label(ibxm(i),2).ne.lstar) then
          if(label(ibxm(i),1).eq.lcr) then
            mcr = mcr + 1
            occcr = occcr + bndocc(i)
          else
            occvl = occvl + bndocc(i)
          end if
        else
          if(label(ibxm(i),1).ne.lry) then
            occvnl = occvnl + bndocc(i)
          end if
        end if
   50 continue
      occevr = totele - occcr - occvl - occvnl
      if(ispin.eq.0) then
        mcr = 2 * mcr
      end if
      mvl = nel - mcr
      mecp = 0
      if(ipseud.ne.0) then
        do 60 i = 1,natoms
          mecp = mecp + iatno(i) - iznuc(i)
   60   continue
        if(ispin.ne.0) mecp = mecp/2
      end if
      mlew = mcr + mvl + mecp
      occlew = occcr + occvl + mecp
      occnon = occvnl + occevr
c
c  write summary of nbo occupancies:
c
      if(jprint(5).eq.1.and.nel.ne.0) then
        write(lfnpr,2000)
        if(ipseud.ne.0) write(lfnpr,2010) dble(mecp)
        if(mcr.ne.0) then
          pcent = occcr/mcr * hundrd
          write(lfnpr,2020) occcr,pcent,mcr
        end if
        if(mvl.ne.0) then
          pcent = occvl/mvl * hundrd
          write(lfnpr,2030) occvl,pcent,mvl
        end if
        write(lfnpr,2040)
        pcent = occlew/mlew * hundrd
        write(lfnpr,2050) occlew,pcent,mlew
        write(lfnpr,2060)
        pcent = occvnl/mlew * hundrd
        write(lfnpr,2070) occvnl,pcent,mlew
        pcent = occevr/mlew * hundrd
        write(lfnpr,2080) occevr,pcent,mlew
        write(lfnpr,2040)
        pcent = occnon/mlew * hundrd
        write(lfnpr,2090) occnon,pcent,mlew
        write(lfnpr,2100)
      end if
c
c  write out nbos:
c
      if(jprint(5).eq.1) then
        write(lfnpr,1000)
        write(lfnpr,1100) (lhyp,j=1,79)
      end if
c
c  main loop over bond orbitals:
c
      nhyb = 0
      mhyb = 0
      ipar3c = 1
      do 180 nbond = 1,nbas
        ib = ibxm(nbond)
        lbl = label(ib,1)
        if(lbl.eq.llp.or.lbl.eq.lcr.or.lbl.eq.lry) nctr = 1
        if(lbl.eq.lbd) nctr = 2
        if(lbl.eq.l3c) nctr = 3
        do 110 i = 1,3
          ia = label(ib,i+3)
          call convrt(ia,ich(i,1),ich(i,2))
          nam(i) = l2blnk
          if(ia.gt.0) nam(i) = nameat(iatno(ia))
          isp(i) = lhyp
          if(i.ge.nctr) isp(i) = lblnk
  110     continue
c
c  loop over atomic centers of bond orbital nbond:
c
        do 170 ictr = 1,nctr
          i = label(ib,ictr+3)
          nel = nameat(iatno(i))
          kl = ll(i)
          ku = ul(i)
          do 120 k = 1,mxao
            ltyp(k) = 0
  120       hyb(k) = zero
c
c  choose sign for polarization coefficients:
c
          isgn = 1
          if(label(ib,2).ne.lstar) go to 130
          if(ictr.lt.2) go to 130
          if(ictr.eq.3) ipar3c = -ipar3c
          if(ictr.eq.3.and.ipar3c.gt.0) go to 130
          isgn = -isgn
  130     continue
c
c  extract hybrid (hyb) from transformation matrix t; ltyp(i) is the
c  orbital angular momentum quantum no. of a.o. # i:
c
          kh = 0
          do 140 k = kl,ku
            kh = kh + 1
            hyb(kh) = t(k,nbond)
  140       ltyp(kh) = naoa(k)/100
          call htype(hyb,ltyp,mxao,kh,coef,pct,nl,isgn)
c
c  find leading non-zero contribution to determine pow(l) for each l:
c
          lstd = 0
          do 160 l = 1,nl
            if(lstd.gt.0) go to 150
            pow(l) = zero
            std = pct(l)
            if(std.lt.thresh) go to 160
            lstd = l
  150       pow(l) = pct(l)/std
            if(pow(l).gt.t99p) pow(l) = t99
  160       continue
c
c  write out nho for center ictr:
c
          coefsq = coef * coef * hundrd
          nl1 = nl
          if(nl1.gt.3) nl1 = 3
          if(ictr.eq.1.and.nctr.eq.1.and.jprint(5).eq.1)
     +      write(lfnpr,1210) nbond,bndocc(nbond),
     +        (label(ib,k),k=1,3),nam(1),ich(1,1),ich(1,2),
     +        pct(1),(lname(l),pow(l),pct(l),l=2,nl1)
          if(ictr.eq.1.and.nctr.gt.1.and.jprint(5).eq.1)
     +      write(lfnpr,1220) nbond,bndocc(nbond),
     +        (label(ib,k),k=1,3),
     +        (nam(k),ich(k,1),ich(k,2),isp(k),k=1,3)
          if(nctr.ne.1.and.jprint(5).eq.1) write(lfnpr,1300) coefsq,
     +        coef,nel,i,pct(1),(lname(l),pow(l),pct(l),l=2,nl1)
          if(nl.gt.3.and.jprint(5).eq.1) write(lfnpr,1310)
     +        (lname(l),pow(l),pct(l),l=4,nl)
          if(iwhybs.ne.0.and.bndocc(nbond).gt.tthoth.and.jprint(5).eq.1)
     +        write(lfnpr,1500) (hyb(k),k=1,kh)
          call frmhyb(hyb,thyb,coef,hycoef,kl,ku,nhyb)
c
c  if this is a new hybrid, form its label:
c
          if(mhyb.ne.nhyb) then
            mhyb = nhyb
            call lblnho(nhyb,nbond,ictr,nctr)
          end if
  170   continue
  180 continue
      return
c
 1000 format(//,1x,'    (Occupancy)   Bond Orbital/ Coefficients/ ',
     + 'Hybrids')
 1100 format(1x,80a1)
 1210 format(1x,i3,'. (',f7.5,') ',a2,a1,'(',i2,')',a2,2a1,12x,
     + ' s(',f6.2,'%)',2(a1,f5.2,'(',f6.2,'%)'))
 1220 format(1x,i3,'. (',f7.5,') ',a2,a1,'(',i2,')',3(a2,3a1))
 1300 format(16x,'(',f6.2,'%)',2x,
     + f7.4,'*',a2,i2,' s(',f6.2,'%)',2(a1,f5.2,'(',f6.2,'%)'))
 1310 format(50x,2(a1,f5.2,'(',f6.2,'%)'))
 1500 format(39x,5f8.4)
 2000 format(/,1x,56('-'))
 2010 format(1x,'  Effective core          ',f9.5)
 2020 format(1x,'  Core                    ',f9.5,' (',f7.3,'% of ',
     +  i3,')')
 2030 format(1x,'  Valence Lewis           ',f9.5,' (',f7.3,'% of ',
     +  i3,')')
 2040 format(2x,18('='),7x,28('='))
 2050 format(1x,'  Total Lewis             ',f9.5,' (',f7.3,'% of ',
     +  i3,')')
 2060 format(2x,53('-'))
 2070 format(1x,'  Valence non-Lewis       ',f9.5,' (',f7.3,'% of ',
     +  i3,')')
 2080 format(1x,'  Rydberg non-Lewis       ',f9.5,' (',f7.3,'% of ',
     +  i3,')')
 2090 format(1x,'  Total non-Lewis         ',f9.5,' (',f7.3,'% of ',
     +  i3,')')
 2100 format(1x,56('-'))
      end
c*****************************************************************************
      subroutine htype(hyb,ltyp,mxao,nh,coef,pct,nl,isgn)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension hyb(mxao),ltyp(mxao),pct(5)
c
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
c  analyze input hybrid 'hyb' for polarization coefficient 'coef'
c  and percentages of each angular momentum component.
c
      data zero,thresh,hundrd/0.0d0,1.d-4,100.0d0/
c
      nl = 0
c
c  zero percentages and polarization coefficient:
c
      do 10 l1 = 1,5
   10   pct(l1) = zero
      coef = zero
c
c  loop over atomic contributions to hybrid, computing percentages
c  and polarization coefficient:
c
      do 20 i = 1,nh
        l1 = ltyp(i) + 1
        if(l1.gt.5) go to 800
        pct(l1) = pct(l1) + hyb(i)**2
   20   coef = coef + hyb(i)**2
      if(abs(coef).lt.thresh) return
c
c  calculate percentage contribution for each angular symmetry:
c
      do 30 l1 = 1,5
   30   pct(l1) = pct(l1)/coef*hundrd
      coef = sqrt(coef)
c
c  switch the sign of the coefficient if isgn is negative:
c
      if(isgn.lt.0) coef = -coef
c
c  normalize the hybrid:
c
      do 50 i = 1,nh
   50   hyb(i) = hyb(i)/coef
c
c  find the maximum number of angular momentum types (nl):
c
      do 60 i = 1,nh
        if(abs(hyb(i)).lt.thresh) go to 60
         if(ltyp(i).le.nl) go to 60
          nl = ltyp(i)
   60   continue
      nl = nl + 1
      return
c
  800 continue
      write(lfnpr,900) l1-1
      stop
c
  900 format(/1x,'ao with unknown angular symmetry, l = ',i3)
      end
c*****************************************************************************
      subroutine frmhyb(hyb,thyb,coef,hycoef,kl,ku,nhyb)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      dimension hyb(1),thyb(ndim,ndim),hycoef(ndim)
c
      data zero,one,thresh/0.0d0,1.0d0,1.0d-4/
c
c  form full nao to nho tranformation in thyb, adding one hybrid with
c  each call.  put polarization coef in hycoef for each hybrid.
c
c  make sure this hybrid isn't already in the list:
c
      if(abs(coef).lt.thresh) return
      do 20 ihyb = 1,nhyb
        temp = zero
        ih = 0
        do 10 k = kl,ku
          ih = ih + 1
          temp = temp + hyb(ih)*thyb(k,ihyb)
   10   continue
        if(abs(abs(temp)-one).lt.thresh) return
        if(abs(temp).gt.thresh) then
          write(lfnpr,900) nhyb+1,temp,ihyb
          stop
        end if
   20 continue
c
c  add this hybrid to the list:
c
      nhyb = nhyb + 1
      if(nhyb.gt.nbas) stop 'too many hybrids'
      do 50 i = 1,nbas
        thyb(i,nhyb) = zero
   50 continue
      ih = 0
      do 70 i = kl,ku
        ih = ih + 1
        thyb(i,nhyb) = hyb(ih)
   70 continue
      hycoef(nhyb) = coef
      if(nhyb.ne.nbas) return
      call svtnho(thyb)
      return
c
  900 format(/1x,'hybrid ',i3,' has a ',
     + 'non-negligible overlap of ',f8.5,' with hybrid ',i3,'.')
      end
c*****************************************************************************
      subroutine hybdir(bndocc,atcoor,thyb,tbnd,scr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),ltyp(maxbas),iathy(maxbas,3)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nbnao/naoc(maxbas),naoa(maxbas),ltyp1(maxbas),
     +       iprin(maxbas)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      dimension bndocc(ndim),atcoor(natoms*3),thyb(ndim,ndim),
     +          tbnd(ndim,ndim),scr(ndim)
      dimension istr(8),phyb(3),xyz(3,2),khyb(3),azi(2),pol(2),dev(2)
      dimension iskip(2)
c
      data lcr,llp,lry,lbd,l3c/2hcr,2hlp,2hry,2hbd,2h3c/
      data lhyp/1h-/
      data zero,one,thresh,cutoff/0.0d0,1.0d0,1.0d-4,1.0d-8/
c
c  compute hybrid directionality and bond bending for selected nbo's:
c
c  thresholds:   athr  --   angular deviation threshold
c                pthr  --   percentage p-character threshold
c                ethr  --   occupancy threshold
c
      conv = 180.0d0/(4.0d0*atan(one))
      write(lfnpr,900) abs(athr),abs(pthr),abs(ethr)
c
c  get atomic centers, nao to nho trans., and nao to nbo trans.:
c
      call fecoor(atcoor)
      call fetnho(thyb)
      call fetnab(tbnd)
      call NBOtransp(tbnd,ndim,nbas)
      call matmlt(tbnd,thyb,scr,ndim,nbas)
c
c  loop over nbos:
c
      icnt = 0
      do 100 ibas = 1,nbas
        ib = ibxm(ibas)
        lbl1 = label(ib,1)
        lbl2 = label(ib,2)
        lbl3 = label(ib,3)
        if(lbl1.eq.llp.or.lbl1.eq.lry) nctr = 1
        if(lbl1.eq.lbd) nctr = 2
c
c  skip 3-center orbitals, core orbitals, low occupancy orbitals:
c
        if(lbl1.eq.l3c) go to 100
        if(lbl1.eq.lcr) go to 100
        if(bndocc(ibas).lt.abs(ethr)) go to 100
c
c  find the hybrids which contribute to this nbo:
c
        ictr = 0
        do 10 ihyb = 1,nbas
          if(abs(tbnd(ibas,ihyb)).gt.thresh) then
            ictr = ictr + 1
            khyb(ictr) = ihyb
          end if
   10   continue
        if(ictr.ne.nctr) then
          write(lfnpr,910) nctr,ibas,ictr
          stop
        end if
c
c  make sure the hybrids are on the proper nuclear centers and compute
c  the percentage p-character in the hybrid:
c
        do 30 ictr = 1,nctr
          ihyb = khyb(ictr)
          jctr = label(ib,ictr+3)
          call hybcmp(xyz(1,ictr),phyb(ictr),ihyb,jctr,thyb(1,ihyb))
   30   continue
c
c  if these hybrids have low p-character, skip them:
c
        iskip(1) = 0
        iskip(2) = 0
        if(nctr.eq.1.and.phyb(1).lt.abs(pthr)) go to 100
        if(nctr.eq.2) then
          if(phyb(1).lt.abs(pthr)) iskip(1) = 1
          if(phyb(2).lt.abs(pthr)) iskip(2) = 1
          if(iskip(1).eq.1.and.iskip(2).eq.1) go to 100
        end if
c
c  compute the polar and azimuthal angles of each hybrid:
c
        do 70 ictr = 1,nctr
          if(iskip(ictr).eq.1) go to 70
          call angles(xyz(1,ictr),xyz(2,ictr),xyz(3,ictr),pol(ictr),
     +                azi(ictr))
   70   continue
c
c  compute the deviation from the line of nuclear centers for 2-center
c  orbitals:
c
        if(nctr.eq.2) then
          ictr = label(ib,4)
          jctr = label(ib,5)
          x = atcoor(jctr*3-2) - atcoor(ictr*3-2)
          y = atcoor(jctr*3-1) - atcoor(ictr*3-1)
          z = atcoor(jctr*3)   - atcoor(ictr*3)
          if(abs(x).lt.cutoff) x = zero
          if(abs(y).lt.cutoff) y = zero
          if(abs(z).lt.cutoff) z = zero
          r = sqrt(x*x + y*y + z*z)
          x = x / r
          y = y / r
          z = z / r
          call angles(x,y,z,theta,phi)
          proj = xyz(1,1)*x + xyz(2,1)*y + xyz(3,1)*z
          if(abs(proj-one).lt.cutoff) then
            dev(1) = zero
          else if(abs(proj+one).lt.cutoff) then
            dev(1) = 180.0d0
          else if(proj.lt.one.and.proj.gt.-one) then
            dev(1) = acos(proj) * conv
            dev(1) = abs(dev(1))
          else
            stop 'arccosine out of bounds in sr hybdir'
          end if
          proj = xyz(1,2)*x + xyz(2,2)*y + xyz(3,2)*z
          if(abs(proj-one).lt.cutoff) then
            dev(2) = 180.0d0
          else if(abs(proj+one).lt.cutoff) then
            dev(2) = zero
          else if(proj.lt.one.and.proj.gt.-one) then
            dev(2) = acos(proj) * conv
            dev(2) = abs(abs(dev(2)) - 180.0d0)
          else
            stop 'arccosine out of bounds in sr hybdir'
          end if
          if(dev(1).lt.abs(athr)) iskip(1) = 1
          if(dev(2).lt.abs(athr)) iskip(2) = 1
          if(iskip(1).eq.1.and.iskip(2).eq.1) go to 100
        end if
c
c  write out directionality info:
c
        icnt = icnt + 1
        istr(1) = lbl1
        istr(2) = lbl2
        istr(3) = lbl3
        istr(4) = nameat(iatno(label(ib,4)))
        istr(5) = label(ib,4)
        if(nctr.eq.2) then
          istr(6) = lhyp
          istr(7) = nameat(iatno(label(ib,5)))
          istr(8) = label(ib,5)
          if(iskip(1).eq.1) then
            write(lfnpr,940) ibas,(istr(i),i=1,8),theta,phi,pol(2),
     +                       azi(2),dev(2)
          else if(iskip(2).eq.1) then
            write(lfnpr,950) ibas,(istr(i),i=1,8),theta,phi,pol(1),
     +                       azi(1),dev(1)
          else
            write(lfnpr,960) ibas,(istr(i),i=1,8),theta,phi,pol(1),
     +                       azi(1),dev(1),pol(2),azi(2),dev(2)
          end if
        else
          write(lfnpr,970) ibas,(istr(i),i=1,5),pol(1),azi(1)
        end if
  100 continue
      if(icnt.eq.0) write(lfnpr,980)
      return
c
  900 format(//1x,'NHO Directionality and "Bond Bending" (deviations ',
     + 'from line of nuclear centers)',//1x,'        [Thresholds for ',
     + 'printing:  angular deviation  > ',f4.1,' degree]',/1x,
     + '                                   hybrid p-character > ',f4.1,
     + '%',/1x,'                                   orbital occupancy  ',
     + '>  ',f4.2,'e',//1x,'                      Line of Centers     ',
     + '   Hybrid 1              Hybrid 2',/1x,'                      ',
     + '---------------  -------------------   ------------------',/1x,
     + '          NBO           Theta   Phi    Theta   Phi    Dev    ',
     + 'Theta   Phi    Dev',/1x,'=====================================',
     + '==========================================')
  910 format(/1x,'error: the ',i1,'-center nbo ',i3,' has ',
     + 'contributions from ',i2,' atomic centers.')
  940 format(1x,i3,'. ',a2,a1,'(',i2,')',a2,i2,a1,a2,i2,3x,f5.1,2x,f5.1,
     + '     --     --    --     ',f5.1,2x,f5.1,1x,f5.1)
  950 format(1x,i3,'. ',a2,a1,'(',i2,')',a2,i2,a1,a2,i2,3x,f5.1,2x,f5.1,
     + 3x,f5.1,2x,f5.1,1x,f5.1,'      --     --    --')
  960 format(1x,i3,'. ',a2,a1,'(',i2,')',a2,i2,a1,a2,i2,3x,f5.1,2x,f5.1,
     + 3x,f5.1,2x,f5.1,1x,f5.1,4x,f5.1,2x,f5.1,1x,f5.1)
  970 format(1x,i3,'. ',a2,a1,'(',i2,')',a2,i2,'          --     --',4x,
     + f5.1,2x,f5.1,'   --       --     --    --')
  980 format(1x,'   none exceeding thresholds')
      end
c*****************************************************************************
      subroutine hybcmp(xyz,pcent,ihyb,jctr,hyb)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension xyz(3),hyb(1)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbnao/naoc(maxbas),naoa(maxbas),ltyp(maxbas),
     +       iprin(maxbas)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      data zero,thresh,cutoff/0.0d0,1.0d-4,1.0d-8/
c
c  add the px,py,pz components of this hybrid vectorially and determine
c  its percentage p-character:
c
      xyz(1) = zero
      xyz(2) = zero
      xyz(3) = zero
      pcent  = zero
      hnorm  = zero
c
c  make sure this hybrid is situated on the correct atom, jctr:
c
      jmax  = 1
      tmax  = abs(hyb(1))
      do 10 inao = 2,nbas
        if(abs(hyb(inao)).gt.tmax) then
          jmax = inao
          tmax = abs(hyb(inao))
        end if
   10 continue
      if(naoc(jmax).ne.jctr) then
        write(lfnpr,920) ihyb,jctr,naoc(jmax)
        stop
      end if
c
c  find the sign of the largest s-component of this hybrid:
c
      jmax  = 0
      tmax  = zero
      do 20 inao = 1,nbas
        l = naoa(inao)/100
        if(l.eq.0.and.abs(hyb(inao)).gt.tmax) then
          jmax = inao
          tmax = abs(hyb(inao))
        end if
   20 continue
c
c  if the sign of the largest s-component is negative, change the
c  phase of this hybrid:
c
      if(jmax.ne.0.and.hyb(jmax).lt.-thresh) then
        do 30 inao = 1,nbas
          hyb(inao) = -hyb(inao)
   30   continue
      endif
c
c  sum the px,py,pz components of this hybrid, determine the percent
c  p-character:
c
      do 40 inao = 1,nbas
        if(naoc(inao).eq.jctr) then
          l = naoa(inao)/100
          if(l.eq.1) then
            pcent = pcent + hyb(inao)*hyb(inao)
            m = mod(naoa(inao),50)
            xyz(m) = xyz(m) + hyb(inao)
          end if
          hnorm = hnorm + hyb(inao)*hyb(inao)
        end if
   40 continue
      if(hnorm.lt.thresh) then
        write(lfnpr,930) jctr,ihyb
        stop
      end if
      pcent = pcent/hnorm * 100.0d0
c
c  normalize the px,py,pz vector:
c
      hnorm = zero
      do 50 ix = 1,3
        if(abs(xyz(ix)).lt.cutoff) xyz(ix) = zero
        hnorm = hnorm + xyz(ix)*xyz(ix)
   50 continue
      hnorm = sqrt(hnorm)
      if(abs(hnorm).lt.cutoff) then
        pcent = zero
      else
        do 60 ix = 1,3
          xyz(ix) = xyz(ix)/hnorm
   60   continue
      end if
      return
c
  920 format(/1x,'expected to find hybrid ',i3,' on nuclear center ',
     + i2,' rather than center ',i2,'.')
  930 format(/1x,'the atomic orbitals on nuclear center ',i2,' do not ',
     + 'contribute to hybrid ',i3,'.')
      end
c*****************************************************************************
      subroutine fndmol(iatoms)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbmol/nmolec,molat(maxatm),molec(maxatm,maxatm),
     +              nmola,molata(maxatm),moleca(maxatm,maxatm)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension iatoms(natoms)
      logical bdfind
c
c  find molecular units :  modified algorithm replacing original which
c  had problems with determining molecular units for odd numberings of
c  atoms.   (e. glendening  3/12/88)
c
      nmolec = 0
      do 20 i = 1,natoms
        molat(i) = 0
        do 10 j = 1,natoms
          molec(i,j) = 0
   10   continue
   20 continue
      do 30 i = 1,natoms
        iatoms(i) = i
   30 continue
      latoms = natoms
   40 continue
        nmolec = nmolec+1
        molat(nmolec) = 1
        molec(nmolec,1) = iatoms(1)
        latoms = latoms-1
        if(latoms.ne.0) then
          do 50 i = 1,latoms
            iatoms(i) = iatoms(i+1)
   50     continue
          iat = 1
   60     continue
            i = 1
   70       continue
              if(bdfind(molec(nmolec,iat),iatoms(i))) then
                molat(nmolec) = molat(nmolec)+1
                molec(nmolec,molat(nmolec)) = iatoms(i)
                latoms = latoms-1
                if(i.le.latoms) then
                  do 80 j = i,latoms
                    iatoms(j) = iatoms(j+1)
   80             continue
                end if
              else
                i = i+1
              end if
            if(i.le.latoms) goto 70
            iat = iat+1
          if(iat.le.molat(nmolec).and.latoms.ne.0) goto 60
        end if
      if(latoms.gt.0) goto 40
c
c  sort atoms in molecular units:
c
      do 110 i = 1,nmolec
        do 100 j = 1,molat(i)-1
          do 90 k = 1,molat(i)-j
            if(molec(i,k).gt.molec(i,k+1)) then
              itemp = molec(i,k)
              molec(i,k) = molec(i,k+1)
              molec(i,k+1) = itemp
            end if
   90     continue
  100   continue
  110 continue
c
c  alpha spin: save bonding info in nmola,molata,moleca:
c
      if(ispin.eq.2) then
        nmola = nmolec
        do 610 imol = 1,nmolec
          molata(imol) = molat(imol)
          imolat = molat(imol)
          do 600 iatmol = 1,imolat
            moleca(imol,iatmol) = molec(imol,iatmol)
  600     continue
  610   continue
c
c  beta spin: make sure that beta molecular units are the same as alpha:
c
      else if(ispin.eq.-2) then
        if(nmola.ne.nmolec) go to 800
        do 730 imol = 1,nmolec
          imolat = molat(imol)
          if(imolat.ne.molata(imol)) go to 800
          do 720 iatmol = 1,imolat
            if(moleca(imol,iatmol).ne.molec(imol,iatmol)) go to 800
  720     continue
  730   continue
      end if
      return
c
  800 write(lfnpr,1800)
      nmola = -nmola
      return
c
 1800 format(/1x,'the molecular units found in the alpha and beta ',
     + 'manifolds are inequivalent.',/1x,'for labelling purposes, ',
     + 'the molecular units of the beta system will be used.')
      end
c*****************************************************************************
      subroutine nbocla(bndocc,accthr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      parameter (maxatm = 750, maxbas = 4096)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),mollst(maxbas),iathy(maxbas,3)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbmol/nmolec,molat(maxatm),molec(maxatm,maxatm),
     +              nmola,molata(maxatm),moleca(maxatm,maxatm)
      dimension bndocc(nbas)
      data lbd,l3c,lstar/2hbd,2h3c,1h*/
      data thresh,one,zero,two/1.50d0,1.0d0,0.0d0,2.0d0/
      data donthr/1.0d-1/
c
c  classify nbos according to donor/acceptor type:
c
      if(accthr.le.zero) then
        accthr = thresh
        if(ispin.ne.0) accthr = accthr - one
      end if
      if(ispin.ne.0) donthr = donthr / two
c
c  make up list mollst of which "molecule" each atom is in:
c
      do 80 iat = 1,natoms
        do 60 imol = 1,nmolec
          imolat = molat(imol)
          do 50 iatmol = 1,imolat
            if(molec(imol,iatmol).eq.iat) go to 70
   50     continue
   60   continue
        stop 'routine nbocla'
   70   mollst(iat) = imol
   80   continue
c
c  make up lists of nbo orbitals:
c    nbouni(ibas) = molecular unit
c    nbotyp(ibas) = number of centers (+10 if a low occupancy lone pair)
c                                     (+20 if an antibond/rydberg)
      do 200 ibas = 1,nbas
        ib = ibxm(ibas)
        iat = label(ib,4)
        imol = mollst(iat)
        nbouni(ibas) = imol
        lab = label(ib,1)
        nctr = 1
        if(lab.eq.lbd) nctr = 2
        if(lab.eq.l3c) nctr = 3
        nbotyp(ibas) = nctr
        if(label(ib,2).eq.lstar) go to 180
        if(bndocc(ibas).gt.accthr) go to 200
c
c  low occupancy valence orbital
c
          nbotyp(ibas) = nctr + 10
          go to 200
c
c  antibond/rydberg
c
  180   nbotyp(ibas) = nctr + 20
c
c  high occupancy ry* or bd* orbital
c
        if(bndocc(ibas).gt.donthr) nbotyp(ibas) = nctr + 10
  200 continue
      return
      end
c*****************************************************************************
      subroutine fnboan(bndocc,f,molnbo)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),mollst(maxbas),iathy(maxbas,3)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbmol/nmolec,molat(maxatm),molec(maxatm,maxatm),
     +              nmola,molata(maxatm),moleca(maxatm,maxatm)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      dimension bndocc(nbas),f(ndim,ndim),molnbo(2,nbas,nmolec)
      dimension inam(3),jnam(3),ich(3,2),jch(3,2),isp(3),jsp(3)
c
      data lbd/2hbd/,l3c/2h3c/,lblnk1/1h /,lblnk2/2h  /,lhyp/1h-/
      data hundth/0.01d0/
      data aukcal/627.51d0/,evkcal/23.060d0/
      data zero,one,two,ten/0.0d0,1.0d0,2.0d0,1.0d1/
c
c  perform 2nd order analysis of the fock matrix:
c
c  ethr1 is the threshold for printing the intramolecular perturbational
c  energies (0.5 kcal/mol for closed shell, 0.25 kcal/mol for open shell)
c  similarly, ethr2 is the intermolecular threshold, (0.05 kcal/mol).
c
      ethr1 = abs(e2thr)
      if(ispin.ne.0.and.e2thr.lt.zero) ethr1 = ethr1/two
      ethr2 = abs(e2thr)/ten
      if(ispin.ne.0.and.e2thr.lt.zero) ethr2 = ethr2/two
c
c  fetch the nbo fock matrix:
c
      ntri = ndim * (ndim+1)/2
      call fefnbo(f)
      call NBOunpack(f,ndim,nbas,ntri)
c
c  analyze fock matrix:
c
c  make up list molnbo(1,ibas,imol) of core/lp/bond nbos in molec. unit imol
c               molnbo(2,ibas,imol) of rydberg/antibond nbos in molec. imol
c
      do 200 imol = 1,nmolec
        nocc = 0
        nstar = 0
        do 110 ibas = 1,nbas
          do 100 i = 1,2
            molnbo(i,ibas,imol) = 0
  100     continue
  110   continue
        do 150 ibas = 1,nbas
          if(imol.ne.nbouni(ibas)) go to 150
          if(nbotyp(ibas).gt.20) go to 130
            nocc = nocc + 1
            molnbo(1,nocc,imol) = ibas
            if(nbotyp(ibas).lt.10) go to 150
  130     continue
            nstar = nstar + 1
            molnbo(2,nstar,imol) = ibas
  150   continue
  200 continue
c
c  determine the conversion from input energy units to kcal:
c
      if(munit.eq.0) then
        conv = aukcal
      else if(munit.eq.1) then
        conv = evkcal
      else
        conv = one
      end if
c
c  loop over pairs of units:
c
      write(lfnpr,2700) ethr1
      if(nmolec.gt.1) write(lfnpr,2710) ethr2
      if(munit.eq.0) then
        write(lfnpr,2720)
      else if(munit.eq.1) then
        write(lfnpr,2730)
      else
        write(lfnpr,2740)
      end if
      do 400 imol = 1,nmolec
        do 400 jmol = 1,nmolec
          if(imol.eq.jmol) write(lfnpr,2300) imol
          if(imol.ne.jmol) write(lfnpr,2400) imol,jmol
          ethrsh = ethr1
          if(imol.ne.jmol) ethrsh = ethr2
          nele = 0
          do 305 iocc = 1,nbas
            ibas = molnbo(1,iocc,imol)
            if(ibas.eq.0) go to 305
            ib = ibxm(ibas)
            lbl = label(ib,1)
            nctr = 1
            if(lbl.eq.lbd) nctr = 2
            if(lbl.eq.l3c) nctr = 3
            do 250 i = 1,3
              ia = label(ib,i+3)
              call convrt(ia,ich(i,1),ich(i,2))
              inam(i) = lblnk2
              if(ia.gt.0) inam(i) = nameat(iatno(ia))
              isp(i) = lhyp
              if(i.ge.nctr) isp(i) = lblnk1
  250       continue
            do 300 jstar = 1,nbas
              jbas = molnbo(2,jstar,jmol)
              if(jbas.eq.0) go to 300
              if(ibas.eq.jbas) go to 300
              de = f(jbas,jbas) - f(ibas,ibas)
              if(de.lt.hundth) go to 300
              absfij = abs(f(ibas,jbas))
              epert = (absfij**2)/de
c
c  compute occupancy factor to multiply by:
c
              totocc = bndocc(ibas)+bndocc(jbas)
              fulloc = two
              if(ispin.ne.0) fulloc = one
              occfac = totocc
              if(totocc.gt.fulloc) occfac = two * fulloc - totocc
c
c  multiply epert by sum of occupancies of nbos ibas and jbas:
c
              epert = epert * occfac
              ekcal = epert * conv
              if(ekcal.lt.ethrsh) go to 300
              nele = nele + 1
              jb = ibxm(jbas)
              lbl = label(jb,1)
              nctr = 1
              if(lbl.eq.lbd) nctr = 2
              if(lbl.eq.l3c) nctr = 3
              do 260 j = 1,3
                ja = label(jb,j+3)
                call convrt(ja,jch(j,1),jch(j,2))
                jnam(j) = lblnk2
                if(ja.gt.0) jnam(j) = nameat(iatno(ja))
                jsp(j) = lhyp
                if(j.ge.nctr) jsp(j) = lblnk1
  260         continue
              write(lfnpr,2800) ibas,(label(ib,k),k=1,3),
     *           (inam(k),ich(k,1),ich(k,2),isp(k),k=1,2),
     *            inam(3),ich(3,1),ich(3,2),
     *                           jbas,(label(jb,k),k=1,3),
     *           (jnam(k),jch(k,1),jch(k,2),jsp(k),k=1,2),
     *            jnam(3),jch(3,1),jch(3,2),
     *                          ekcal,de,absfij
  300   continue
  305   continue
        if(nele.eq.0) write(lfnpr,2500)
  400 continue
      return
c
 2300 format(/1x,'within unit ',i2)
 2400 format(/1x,'from unit ',i2,' to unit ',i2)
 2500 format(1x,'      none above threshold')
 2700 format(//,1x,'Second Order Perturbation Theory Analysis ',
     *             'of Fock Matrix in NBO Basis'//,1x,
     *          '    Threshold for printing:  ',f5.2,' kcal/mol')
 2710 format(1x,'   (intermolecular threshold:',f5.2,' kcal/mol)')
 2720 format(56x,'  E(2)  E(j)-E(i) F(i,j)'/
     * 6x,'Donor NBO (i)',14x,'Acceptor NBO (j)',7x,
     *            'kcal/mol   a.u.    a.u. ',/1x,79('='))
 2730 format(56x,'  E(2)  E(j)-E(i) F(i,j)'/
     * 6x,'Donor NBO (i)',14x,'Acceptor NBO (j)',7x,
     *            'kcal/mol   e.v.    e.v. ',/1x,79('='))
 2740 format(56x,'  E(2)  E(j)-E(i) F(i,j)'/
     * 6x,'Donor NBO (i)',14x,'Acceptor NBO (j)',7x,
     *            'kcal/mol   kcal    kcal ',/1x,79('='))
 2800 format(1x,i3,'. ',a2,a1,'(',i2,')',a2,3a1,a2,3a1,a2,2a1,
     *   '/',i3,'. ',a2,a1,'(',i2,')',a2,3a1,a2,3a1,a2,2a1,
     *       f8.2,f8.2,f9.3)
      end
c*****************************************************************************
      subroutine nbosum(f,bndocc,list,lista,scr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical first
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbmol/nmolec,molat(maxatm),molec(maxatm,maxatm),
     +              nmola,molata(maxatm),moleca(maxatm,maxatm)
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      dimension f(ndim,ndim),bndocc(ndim),list(ndim),lista(natoms,2),
     +          scr(1)
      dimension istr(80),ilab(9)
c
      data zero,eps,two,ten,hundrd/0.0d0,5.0d-6,2.0d0,1.0d1,1.0d2/
      data tenth/1.0d-1/
      data lstar,lry/1h*,2hry/
c
c  set flag to zero -- determine strong delocalizations from perturbative
c  analysis of the nbo fock matrix:
c
      iflg = 0
c
c  threshold for printing delocalizations:
c
      thr1 = abs(e2thr)
      if(ispin.ne.0) thr1 = thr1/two
      thr2 = thr1 / ten
c
c  get fock matrix if there is one:
c
      if(iwfock.ne.0) then
        ntri = ndim * (ndim+1)/2
        call fefnbo(f)
        call NBOunpack(f,ndim,nbas,ntri)
      end if
c
c  print summary heading, then loop over molecules:
c
      if(iwfock.ne.0) then
        write(lfnpr,900)
      else
        write(lfnpr,910)
      end if
      do 200 imol = 1,nmolec
c
c  determine the molecular formula, the nuclear charge, and the number of
c  ecp electrons of this molecular unit:
c
        nat  = 0
        mecp = 0
        charge = zero
        do 20 iat = 1,molat(imol)
          kat = iatno(molec(imol,iat))
          mecp = mecp + dble(kat - iznuc(molec(imol,iat)))
          charge = charge + dble(kat)
          do 10 jat = 1,nat
            if(lista(jat,1).eq.kat) then
              lista(jat,2) = lista(jat,2) + 1
              go to 20
            end if
   10     continue
          nat = nat + 1
          lista(nat,1) = kat
          lista(nat,2) = 1
   20   continue
        if(ispin.ne.0) mecp = mecp/2
        if(ispin.ne.0) charge = charge/two
        call chem(nat,natoms,lista,nl,istr)
        write(lfnpr,920) imol,(istr(i),i=1,nl)
c
c  loop over nbo's on this molecular unit:
c
        occlew = dble(mecp)
        occnon = zero
        occryd = zero
        do 190 ibas = 1,nbas
          if(nbouni(ibas).eq.imol) then
            ib = ibxm(ibas)
            ilab(1) = label(ib,1)
            ilab(2) = label(ib,2)
            ilab(3) = label(ib,3)
            iptr    = 3
            nctr    = mod(nbotyp(ibas),10)
            do 30 ictr = 1,nctr
              iptr         = iptr + 2
              ilab(iptr)   = label(ib,ictr+3)
              ilab(iptr-1) = nameat(iatno(ilab(iptr)))
   30       continue
            occ  = bndocc(ibas)
            if(ilab(1).eq.lry) then
              occryd = occryd + occ
            else if(ilab(2).eq.lstar) then
              occnon = occnon + occ
            else
              occlew = occlew + occ
            end if
c
c  if there is a fock matrix, find the orbital energy and principal
c  delocalizations:
c
            if(iwfock.ne.0) then
              enrg  = f(ibas,ibas)
              call getdel(ibas,occ,thr1,thr2,nl,list,scr,f,iflg)
              first = .true.
              il    = 0
   40         call dlcstr(ibas,il,nl,list,ml,istr)
              if(first) then
                if(nctr.eq.1) then
                  write(lfnpr,930) ibas,(ilab(i),i=1,iptr),occ,enrg,
     +                             (istr(j),j=1,ml)
                else if(nctr.eq.2) then
                  write(lfnpr,940) ibas,(ilab(i),i=1,iptr),occ,enrg,
     +                             (istr(j),j=1,ml)
                else
                  write(lfnpr,950) ibas,(ilab(i),i=1,iptr),occ,enrg,
     +                             (istr(j),j=1,ml)
                end if
                first = .false.
              else
                  write(lfnpr,960) (istr(j),j=1,ml)
              end if
              if(il.lt.nl) go to 40
c
c  otherwise only write out orbital labels and occupancy:
c
            else
              if(nctr.eq.1) then
                write(lfnpr,930) ibas,(ilab(i),i=1,iptr),occ
              else if(nctr.eq.2) then
                write(lfnpr,940) ibas,(ilab(i),i=1,iptr),occ
              else
                write(lfnpr,950) ibas,(ilab(i),i=1,iptr),occ
              end if
            end if
          end if
  190   continue
        write(lfnpr,970)
        total = occlew + occnon + occryd
c
c  make sure the total number of electrons is an integer if there is only
c  one molecular unit:
c
        if(nmolec.eq.1) then
          total  = total + tenth
          nel    = total
          total  = nel
          occryd = total - occlew - occnon
        end if
c
c  write a summary of the electron population on this molecular unit:
c
        if(abs(total-dble(nint(total))).lt.1.0d-5)
     +                    total = dble(nint(total))
        charge = charge - total
        if(total.gt.eps) then
          plew = occlew/total*hundrd
          pnon = occnon/total*hundrd
          pryd = occryd/total*hundrd
        else
          plew = zero
          pnon = zero
          pryd = zero
        end if
        write(lfnpr,980) occlew,plew
        write(lfnpr,990) occnon,pnon
        write(lfnpr,1000) occryd,pryd
        write(lfnpr,970)
        write(lfnpr,1010) imol,total,hundrd
        write(lfnpr,1020) imol,charge
        if(imol.lt.nmolec) write(lfnpr,*)
  200 continue
      return
c
  900 format(//1x,'Natural Bond Orbitals (Summary):',//53x,'Principal ',
     + 'Delocalizations',/1x,'          NBO              Occupancy  ',
     + '  Energy      (geminal,vicinal,remote)',/1x,79('='))
  910 format(//1x,'Natural Bond orbitals (Summary):',//1x,'          ',
     + 'NBO              Occupancy  ',/1x,40('-'))
  920 format(1x,'Molecular Unit ',i2,'  ',60a1)
  930 format(1x,i3,'. ',a2,a1,'(',i2,')',a2,i2,10x,f9.5,f12.5,4x,28a1)
  940 format(1x,i3,'. ',a2,a1,'(',i2,')',a2,i2,'-',a2,i2,5x,f9.5,f12.5,
     + 4x,28a1)
  950 format(1x,i3,'. ',a2,a1,'(',i2,')',a2,i2,'-',a2,i2,'-',a2,i2,f9.5,
     + f12.5,4x,28a1)
  960 format(52x,28a1)
  970 format(1x,'      -------------------------------')
  980 format(1x,'             Total Lewis',f11.5,'  (',f8.4,'%)')
  990 format(1x,'       Valence non-Lewis',f11.5,'  (',f8.4,'%)')
 1000 format(1x,'       Rydberg non-Lewis',f11.5,'  (',f8.4,'%)')
 1010 format(1x,'           Total unit ',i2,f11.5,'  (',f8.4,'%)')
 1020 format(1x,'          Charge unit ',i2,f11.5)
      end
c*****************************************************************************
      subroutine getdel(ibo,occ,thr1,thr2,nl,list,del,deloc,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      dimension list(ndim),del(ndim),deloc(ndim,ndim)
c
      data zero,one,cutoff,tenth/0.0d0,1.0d0,1.0d-4,0.1d0/
      data aukcal,evkcal/627.51,23.060/
c
c determine the conversion factor to kcal:
c
      if(munit.eq.0) then
        conv = aukcal
      else if(munit.eq.1) then
        conv = evkcal
      else
        conv = one
      end if
c
c determine the strength of each delocalization:
c
      do 10 jbo = 1,nbas
        list(jbo) = 0
        del(jbo) = zero
   10 continue
c
      nl = 0
      if(occ.lt.tenth) return
      do 20 jbo = 1,nbas
        if(ibo.ne.jbo) then
          if(nbotyp(jbo).ge.10) then
            del(jbo) = deloc(ibo,jbo)*deloc(ibo,jbo)
            if(iflg.eq.0) then
              div = abs(deloc(ibo,ibo)-deloc(jbo,jbo))
              if(div.ne.zero) then
                del(jbo) = occ * del(jbo)/div * conv
              else
                del(jbo) = zero
              end if
            end if
          end if
          if(del(jbo).gt.thr2.and.nbouni(ibo).ne.nbouni(jbo)) then
            nl = nl + 1
            list(nl) = jbo
          else if(del(jbo).gt.thr1) then
            nl = nl + 1
            list(nl) = jbo
          end if
        end if
   20 continue
c
c  sort delocalizations:
c
      do 100 i = 1,nl
        do 90 j = 1,nl-1
          kbo = list(j)
          lbo = list(j+1)
          if(del(lbo)-del(kbo).gt.cutoff) then
            list(j) = lbo
            list(j+1) = kbo
          end if
   90   continue
  100 continue
      return
      end
c*****************************************************************************
      subroutine dlcstr(ibo,il,nl,list,ml,istr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxchr = 28, maxd = 4)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
      dimension list(ndim),istr(80)
      integer ik(maxd)
c
      data icomma,ileft,iright/1h,,1h(,1h)/
c
c  build a character string (for the nbo summary table) which contains
c  the delocalization information for nbo # ibo:
c
      ml = 0
   10 il = il + 1
      if(il.gt.nl) go to 30
      call idigit(list(il),ik,nd,maxd)
      if(ml+nd+4.gt.maxchr) go to 30
      if(ml.ne.0) then
        ml = ml + 1
        istr(ml) = icomma
      end if
      do 20 i = 1,nd
        ml = ml + 1
        istr(ml) = ik(i)
   20 continue
      ml = ml + 1
      istr(ml) = ileft
      ml = ml + 1
      istr(ml) = ihtyp(ibo,list(il))
      ml = ml + 1
      istr(ml) = iright
      go to 10
c
   30 il = il - 1
      return
      end
c*****************************************************************************
      subroutine nlmo(n,a,eval,evec,tsym,reson,nocc,ialarm)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  form natural localized molecular orbitals from density matrix a.
c
c        n: actual dimension of a,evec
c     ndim: declared dimension of a,evec
c     tsym: scratch
c    reson: squares of diagonal elements of nbo to nlmo transf, times 100%
c   ialarm: alarm that the orbital occupancies are out of order and that
c           the lmo step should be avoided
c
c  these values are set:
c
c     differ = 1.0d-5
c
c     done   = 1.0d-10 (this is the parameter for convergence of the off-
c                       diagonal matrix elements.)
c
c     eps    = 1.0d-11 (this parameter has to do with the machine precision
c                       and should be set to a value between "done" and the
c                       machine precision.)
c
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical zeroj
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      dimension a(ndim,ndim),evec(ndim,1),eval(1),tsym(1),reson(ndim)
      dimension rot(2,2)
      dimension ilist(100),jlist(100),ioff(100),joff(100),iuniq(100),
     +   juniq(100)
c
c  important parameters:
c
      data differ,done,eps/1.0d-5,1.0d-10,1.0d-11/
c
c  noffmx is set to the dimension of vectors ilist,jlist,ioff,joff,iuniq,juniq:
c
      data degthr,noffmx/1.0d-3,100/
      data zero,one,ten,hundrd/0.0d0,1.0d0,10.0d0,100.0d0/
c
      write(lfnpr,8390)
      thr1 = one - degthr
      thr2 = one - degthr*5
      ntime = 0
c
c  if there is only one basis function, solve this trivial case and return:
c
      if(n.gt.1) go to 10
        evec(1,1) = one
        eval(1) = a(1,1)
        return
   10 continue
c
      do 30 j = 1,n
        do 20 i = 1,n
   20     evec(i,j) = zero
   30   evec(j,j) = one
c
c  count the number of electrons and occupied orbitals:
c
      totele = zero
      do 50 i = 1,n
   50   totele = totele + a(i,i)
      totele = totele + differ
      nocc = totele
      if(ispin.eq.0) nocc = nocc/2 + mod(nocc,2)
      nvirst = nocc + 1
c
c  check if occupancies are in order:
c
      ialarm = 0
      virmax = zero
      do 60 j = nvirst,n
        if(a(j,j).lt.virmax) go to 60
        virmax = a(j,j)
   60 continue
      occmin = hundrd
      do 70 i = 1,nocc
        if(a(i,i).gt.occmin) go to 70
        occmin = a(i,i)
   70 continue
      x = occmin - virmax
c
c  21 oct 1987.  the following feature of the program has been
c    turned off because sometimes it is not possible to diagonalize
c    the nbo density matrix when one of the `a' nbos is degenerate
c    in occupancy with one or more `b' nbos:
c
c  the "abs(x).lt.differ" part of the next line is included so that
c   nlmos can be computed when a number of orbitals are nearly
c   degenerate in occupancy, as for instance in cli6, where six
c   lithium lone pairs are degenerate but only one of them can
c   be placed in the "occupied" set of nlmos.
c     if(x.gt.zero.or.abs(x).lt.differ) go to 100
c
c  the above statement is replaced by:
c
      if(x.gt.differ) go to 100
c
c  occupancies out of order:
c
      ialarm = 1
      if(abs(x).gt.differ) go to 80
        write(lfnpr,8010)
        go to 90
   80   write(lfnpr,8000)
   90 continue
      return
c
c   start loop:
c
  100 continue
      ntime = ntime + 1
c
c  first, find element a(iocc,jemt) of largest magnitude, offtop:
c
      offtop = zero
      do 200 jemt = nvirst,n
        do 200 iocc = 1,nocc
          absaij = abs(a(iocc,jemt))
          if(absaij.lt.offtop) go to 200
          offtop = absaij
          aii = a(iocc,iocc)
          ajj = a(jemt,jemt)
  200 continue
c
c  return if convergence has been achieved:
c
      if(offtop.lt.done) go to 900
c
c  find all elements degenerate with largest one, offtop:
c  (check corresponding diagonal elements also)
c  noff: number of degenerate elements
c  ioff(k),joff(k): kth degenerate element
c
      offtst = offtop * thr1
      aiil = aii*thr2
      ajjl = ajj*thr2
      aiiu = aii/thr2
      ajju = ajj/thr2
      zeroj = .false.
      if(ajj.lt.differ) zeroj = .true.
      noff = 0
      do 250 jemt = nvirst,n
        do 250 iocc = 1,nocc
          absaij = abs(a(iocc,jemt))
          if(absaij.lt.offtst) go to 250
          aiii = a(iocc,iocc)
          ajjj = a(jemt,jemt)
          if((aiii.lt.aiil).or.(aiii.gt.aiiu)) go to 250
c
c  skip test of diag. elem. if small (.lt.differ):
c
          if(zeroj) go to 240
          if((ajjj.lt.ajjl).or.(ajjj.gt.ajju)) go to 250
  240     noff = noff + 1
          ioff(noff) = iocc
          joff(noff) = jemt
  250     continue
      if(noff.lt.noffmx) go to 260
        write(lfnpr,2500) noff,noffmx
 2500   format(//1x,'noff = ',i5,' is greater than noffmx =',i5,
     *            /5x,'  must abort nlmo procedure')
        ialarm = 1
        return
  260 continue
c
      s = ajj - aii
      abss = abs(s)
c
c  if the rotation is very close to 45 degrees, set sin and cos to 1/(root 2)
c
      test=eps*offtop
      if (abss.gt.test) go to 330
      s=.707106781d0
      c=s
      go to 340
c
c  calculation of sin and cos for rotation that is not very close to 45 degrees
  330 t=offtop/s
      s=0.25d0/ sqrt(0.25d0+t*t)
c
c    jacobi rotation angle:   cos=c ,  sin=s
      c= sqrt(0.5d0+s)
      s=2.d0*t*s/c
  340 continue
c  print statements for nlmo procedure details:
c      write(lfnpr,9903) offtop,s,c,noff
c 9903 format(' ******   offtop,s,c,noff:',3f14.9,i3)
c      write(lfnpr,9901) (ioff(i),i=1,noff)
c 9901 format(' ioff:',20i3)
c      write(lfnpr,9902) (joff(i),i=1,noff)
c 9902 format(' joff:',20i3)
c
c     simple 2 by 2 rotation, no degeneracy problems:
      if(noff.gt.1) go to 400
        iocc=ioff(1)
        jemt=joff(1)
        if(a(iocc,jemt).lt.zero) s=-s
        rot(1,1)=c
        rot(2,2)=c
        rot(1,2)=s
        rot(2,1)=-s
        ioff(2)=joff(1)
        call limtrn(a,ioff,rot,eval,ndim,n,2,2,0)
c
c     rotation completed
        do 380 i=1,n
          t=evec(i,iocc)
          evec(i,iocc)=c*t-evec(i,jemt)*s
  380     evec(i,jemt)=s*t+evec(i,jemt)*c
        go to 800
c
  400 continue
c
c  noff.gt.1:
c   compute "averaged" unitary transformation so that symmetry is preserved
c
c    construct unique lists of orbitals involved:
c
c      iuniq(l): l-th unique occupied orb.
c      niuniq:   no. of unique occ. orbs
c      ilist(l): location in the unique list (iuniq) of the i value of the
c                            l-th offdiag. element
c      juniq, njuniq, and jlist are for the empty orbitals.
c
        iuniq(1)=ioff(1)
        ilist(1)=1
        niuniq=1
        do 500 moff=2,noff
          i=ioff(moff)
          iimax=moff-1
          do 490 ii=1,iimax
            if(ioff(ii).ne.i) go to 490
            ilist(moff)=ilist(ii)
            go to 500
  490       continue
          niuniq=niuniq+1
          ilist(moff)=niuniq
          iuniq(niuniq)=i
  500     continue
c
        juniq(1)=joff(1)
        jlist(1)=niuniq+1
        njuniq=1
        do 540 moff=2,noff
          j=joff(moff)
          jjmax=moff-1
          do 530 jj=1,jjmax
            if(joff(jj).ne.j) go to 530
            jlist(moff)=jlist(jj)
            go to 540
  530       continue
          njuniq=njuniq+1
          jlist(moff)=njuniq+niuniq
          juniq(njuniq)=j
  540     continue
        nrot=niuniq+njuniq
        nrot2=nrot*nrot
        n1=nrot2+1
        n2=nrot2+n1
c  construct tsym:
        call symuni(tsym,a,c,s,tsym(n1),tsym(n2),eval,nrot,
     *              niuniq,njuniq,
     *              ilist,jlist,noff,ioff,joff,ndim)
c
c   make iuniq into a complete list of the unique orbitals, and transform
c    the nbo to nlmo transf. (evec) and the dm (a) by tsym:
        ii=niuniq
        do 700 i=1,njuniq
          ii=ii+1
  700     iuniq(ii)=juniq(i)
        call limtrn(evec,iuniq,tsym,eval,ndim,n,nrot,nrot,1)
        call limtrn(a,iuniq,tsym,eval,ndim,n,nrot,nrot,0)
c  see how much the elements were reduced:
c        do 750 moff=1,noff
c          i=ioff(moff)
c          j=joff(moff)
c          write(lfnpr,9920) i,j,(a(i,j))
c 9920     format(' i,j,aij:',2i3,f14.9)
c  750     continue
c
  800   continue
c      totele=zero
c      do 810 j=1,n
c        totele=totele+a(j,j)
c  810   continue
c      tot=nel
c      fract=totele-tot
c      write(lfnpr,7000) noff,totele,fract
      go to 100
c
c  finished: place occupancies in eval and count up electrons:
c
  900 continue
      totele = zero
      do 910 j = 1,n
        eval(j) = a(j,j)
        totele = totele + eval(j)
        x = evec(j,j)
        reson(j) = x * x * hundrd
  910 continue
      totp = totele + differ
      nel = totp
      tot = nel
      fract = abs(totele-tot)
      if(fract.gt.differ) go to 990
c
c  find the largest off-diagonal density matrix element:
c
      amax = zero
      do 960 j = 2,n
        jm1 = j - 1
        do 950 i = 1,jm1
          if(abs(a(i,j)).lt.amax) go to 950
          amax = abs(a(i,j))
  950   continue
  960 continue
      write(lfnpr,9500) amax
c
c  if this is a correlated wavefunction, return to the calling routine:
c
      if(ci.or.mcscf.or.auhf) return
c
c  for scf wavefunctions, make sure this matrix element is small:
c
      if(amax.lt.hundrd*hundrd*done) return
      write(lfnpr,9550)
      ialarm = 1
      return
c
c  non-integer number of electrons:
c
  990 write(lfnpr,9900) differ,totele
      write(lfnpr,9600)
      write(lfnpr,9610) (eval(i),i=1,nbas)
      ialarm = 1
      return
c
 8000 format(/1x,'highest occupied nbos are not at the beginning',
     +    ' of the nbo list;',/,1x,'the nlmo program is not ',
     +    'currently set up to handle this.')
 8010 format(/1x,'degeneracy between orbitals in the (a) and (b)',
     *     ' sets detected;',
     *    /1x,'nlmo program cannot always handle this situation.')
 8390 format(//1x,'natural localized molecular orbital (nlmo) ',
     *     'analysis:')
 9500 format(/1x,'maximum off-diagonal element of dm in nlmo basis:',
     *         e13.5)
 9550 format(/1x,'something went wrong in the nlmo procedure; density',
     * ' matrix of scf',/1x,'wave function has not been diagonalized')
 9600 format(/1x,'occupancies of nlmos:')
 9610 format(/1x,8f10.5)
 9900 format(/1x,'number of electrons (trace of dm, nlmo basis) is not',
     * ' within ',f10.5/' of an integer:',f10.5,' - - program abort')
      end
c*****************************************************************************
      subroutine lmoanl(t,s,reson,occ,ts,border,owbord,atlmo,
     *                  siab,nocc,nab)
c*****************************************************************************
c revision 1.2  88/03/03  11:29:56  reed
c to reduce amount of output, deleted some blank lines, commented out print
c of atom totals for bond orders, and the atomic contrib. to the nlmo is
c only printed if it is greater than 0.01%.
c
      implicit real*8 (a-h,o-z)
      integer ul
      logical closed
c
c  print out details of nao to nlmo transformation in matrix t.
c
c  required input:
c      ndim = declared dimensionality of array t
c      nbas = no. of orbitals = actual dimension of t, naol
c      naol = integer list of orbital ang. momentum type
c                naol(i)/100 = l = q.n. of atomic orbital i
c     iatno = list of atomic numbers; iatno(i) is the nuclear charge
c                of atom i as an integer
c    natoms = no. of atoms (not including ghosts) in the molecule
c    iwhybs = 1 if hybrid a.o. coefficients are to be printed,
c             0 otherwise.
c     lfnpr = logical file number for printout.
c    naoctr = list of atomic centers of oao or nao basis orbitals
c     label = list of bond orbital labels
c      ibxm = permutation list of bond orbitals
c    bndocc = list of bond orbital occupancies
c     ispin = 0 for closed shell
c           = 2 for alpha spin
c           =-2 for beta  spin
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbnao/naoctr(maxbas),naol(maxbas),ltyp1(maxbas),
     +       iprin(maxbas)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       ul(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),ltyp(maxbas),iathy(maxbas,3)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      dimension t(ndim,ndim),s(ndim,ndim),occ(ndim),reson(ndim),
     * ts(ndim),siab(nocc,nab),atlmo(nocc,natoms),
     * border(natoms,natoms),owbord(natoms,natoms),
     * pct(5),pow(5),lname(5),isp(3),nam(3),ich(3,2)
      character*80 title
      data llp,lbd,l3c,lcr,lry/2hlp,2hbd,2h3c,2hcr,2hry/
      data lname/1hs,1hp,1hd,1hf,1hg/
      data zero,hundth,t99,t99p/0.0d0,1.d-2,99.99d0,99.995d0/
      data two,tenth,hundrd,thr/2.0d0,0.1d0,100.0d0,1.0d-6/
      data lhyp,lblnk,l2blnk/1h-,1h ,2h  /
      data bothr/2.0d-3/
c
      closed=.true.
      if(ispin.ne.0) closed=.false.
      if(ispin.eq.0) write(lfnpr,8400)
      if(ispin.eq.2) write(lfnpr,8410)
      if(ispin.eq.-2) write(lfnpr,8420)
      write(lfnpr,8000)
      write(lfnpr,8100) (lhyp,j=1,79)
c  loop over occupied nlmos:
      do 900 nlmo=1,nbas
        if(occ(nlmo).lt.tenth) go to 900
        ib=ibxm(nlmo)
        lbl=label(ib,1)
        if(lbl.eq.llp.or.lbl.eq.lcr.or.lbl.eq.lry) nctr=1
        if(lbl.eq.lbd) nctr=2
        if(lbl.eq.l3c) nctr=3
        do 110 i=1,3
          ia=label(ib,i+3)
          call convrt(ia,ich(i,1),ich(i,2))
          nam(i)=l2blnk
          if(ia.gt.0) nam(i)=nameat(iatno(ia))
          isp(i)=lhyp
          if(i.ge.nctr) isp(i)=lblnk
  110     continue
c  loop over atomic centers of bond orbital nbond
        do 170 ictr=1,nctr
          isp(ictr)=lhyp
          if(ictr.eq.nctr) isp(ictr)=lblnk
          i=label(ib,ictr+3)
          nel=nameat(iatno(i))
  170     continue
          write(lfnpr,8220) nlmo,occ(nlmo),reson(nlmo),(label(ib,k),
     +                k=1,3),(nam(k),ich(k,1),ich(k,2),isp(k),k=1,3)
          if(occ(nlmo).lt.tenth.and.lbl.eq.lry) go to 900
c  loop over atoms:  (j counts over naos)
        do 700 iat=1,natoms
          nl=0
          do 200 l=1,5
  200       pct(l)=zero
          jlow=ll(iat)
          jhigh=ul(iat)
          do 300 j=jlow,jhigh
            l=naol(j)/100+1
            coef=t(j,nlmo)
            pct(l)=pct(l)+coef*coef
  300       continue
c  print out contribution from atom iat (and save in atlmo):
          nl=l
          pol=zero
          do 340 l=1,5
  340       pol=pol+pct(l)
          if(nlmo.le.nocc) atlmo(nlmo,iat)=pol
          pctpol=pol*hundrd
c  print only contributions greater than 0.01%
          if(pctpol.lt.hundth) go to 700
          do 350 l=1,5
  350       pct(l)=hundrd*pct(l)/pol
c  find leading non-zero contribution to determine pow(l) for each l
          lstd=0
          do 460 l=1,nl
            if(lstd.gt.0) go to 450
             pow(l)=zero
             std=pct(l)
             if(std.lt.hundth) go to 460
              lstd=l
  450       pow(l)=pct(l)/std
             if(pow(l).gt.t99p) pow(l)=t99
  460     continue
          nl1=nl
          nel=nameat(iatno(iat))
          if(nl1.gt.3) nl1=3
          write(lfnpr,8300)
     *        pctpol,nel,iat,pct(1),(lname(l),pow(l),pct(l),l=2,nl1)
          if(nl.gt.3) write(lfnpr,8310)
     *        (lname(l),pow(l),pct(l),l=4,nl)
  700     continue
  900   continue
c
c  now, compute hybrid overlaps siab:
c
      if(ortho) goto 2200
      call fesnao(s)
      do 1500 nlmo=1,nocc
        iab=0
        natm1=natoms-1
        do 1400 iat=1,natm1
          ialow=ll(iat)
          iahigh=ul(iat)
          do 1100 l=1,nbas
            if(l.ge.ialow.and.l.le.iahigh) go to 1100
            ts(l)=zero
            do 1050 k=ialow,iahigh
 1050         ts(l)=ts(l)+t(k,nlmo)*s(k,l)
 1100       continue
c          if(iat.gt.2) go to 1130
c          call altout(ts,1,ndim,1,ndim)
c 1130     continue
          jat0=iat+1
          do 1300 jat=jat0,natoms
            iab=iab+1
            ovp=zero
            jalow=ll(jat)
            jahigh=ul(jat)
            do 1200 l=jalow,jahigh
 1200         ovp=ovp+ts(l)*t(l,nlmo)
            anorm=sqrt(atlmo(nlmo,iat)*atlmo(nlmo,jat))
            if(anorm.lt.thr) go to 1250
            siab(nlmo,iab)=ovp/anorm
c            if(iat.gt.2) go to 1300
c            write(lfnpr,9996) jat,iab,jalow,jahigh,ovp,anorm,
c     *                    siab(nlmo,iab)
c 9996       format(1x,'jat,iab,jalow,jahigh,ovp,anorm,siab:',
c     *              /5x,4i3,3f11.6)
            go to 1300
 1250       siab(nlmo,iab)=zero
c            if(iat.gt.2) go to 1300
c            write(lfnpr,9996) jat,iab,jalow,jahigh,ovp,anorm,
c     *                    siab(nlmo,iab)
 1300       continue
 1400     continue
 1500   continue
c  now we are ready to compute bond orders!
      if(jprint(12).ne.0) then
        iab=0
        natm1=natoms-1
        write(lfnpr,9000)
        do 2000 iat=1,natm1
          jat0=iat+1
          do 1900 jat=jat0,natoms
            iab=iab+1
            sum=zero
            owsum=zero
            do 1800 nlmo=1,nocc
              alama2=atlmo(nlmo,iat)
              alamb2=atlmo(nlmo,jat)
              ovp=siab(nlmo,iab)
              bo=alama2
              if(alamb2.lt.alama2) bo=alamb2
c              write(lfnpr,8999) alama2,alamb2,bo
c 8999         format(1x,'alama2,alamb2,bo:',3f14.7)
              if(closed) bo=bo*two
              owbo=bo*ovp
              if(ovp.lt.zero) bo=-bo
              if(abs(bo).gt.bothr)
     *          write(lfnpr,9100) iat,jat,nlmo,bo,ovp
              sum=sum+bo
              owsum=owsum+owbo
 1800         continue
c            write(lfnpr,9110) sum,owsum
            border(iat,jat)=sum
            border(jat,iat)=sum
            owbord(iat,jat)=owsum
            owbord(jat,iat)=owsum
 1900       continue
 2000     continue
c  zero diagonal elements!
        do 2020 iat=1,natoms
          border(iat,iat)=zero
 2020     owbord(iat,iat)=zero
c  compute totals by atom and print results:
        do 2100 iat=1,natoms
          sum=zero
          do 2050 jat=1,natoms
            sum=sum+border(iat,jat)
 2050     continue
          ts(iat)=sum
 2100   continue
        title = 'atom-atom net linear nlmo/npa bond orders:'
        call aout(border,natoms,natoms,natoms,title,0,natoms)
        title = 'linear nlmo/npa bond orders, totals by atom:'
        call aout(ts,natoms,natoms,1,title,0,1)
      end if
 2200 continue
      return
c
 8000 format(1x,'nlmo/occupancy/percent from parent nbo/ atomic ',
     + 'hybrid contributions')
 8100 format(1x,80a1)
 8220 format(1x,i3,'. (',f7.5,') ',f8.4,'%  ',a2,a1,'(',i2,')',
     + 3(a2,3a1))
 8300 format(26x,f7.3,'% ',a2,i2,' s(',f6.2,'%)',2(a1,f5.2,'(',
     +  f6.2,'%)'))
 8310 format(50x,2(a1,f5.2,'(',f6.2,'%)'))
 8400 format(/1x,'hybridization/polarization analysis of nlmos ',
     *  'in nao basis:')
 8410 format(/1x,'hybridization/polarization analysis of nlmos ',
     *  'in nao basis, alpha spin:')
 8420 format(/1x,'hybridization/polarization analysis of nlmos ',
     *  'in nao basis, beta spin:')
 9000 format(/1x,'individual lmo bond orders greater than 0.002',
     *   ' in magnitude,'/1x,
     * 'with the overlap between the hybrids in the nlmo given:',//1x,
     *   'atom i / atom j / nlmo / bond order / hybrid overlap /')
 9100 format(1x,i4,i8,2x,i6,f14.7,f16.7)
      end
c*****************************************************************************
      subroutine dipanl(dm,t,c,tnbo,dx,dy,dz,scr,index)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical test
c
c  dm       --  nlmo density matrix (input)
c  t        --  ao to nlmo transformation matrix (input)
c  c        --  nbo to nlmo transformation matrix (retrieved from nbodaf)
c  tnbo     --  ao to nbo transformation (retrieved from nbodaf)
c  dx,dy,dz --  ao dipole matrices (retrieved from nbodaf)
c  scr      --  ndim*ndim word scratch vector
c  index    --  temporary indexing array
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbdxyz/xdip,ydip,zdip,charge(maxatm)
      common/nbatom/iatno(maxatm),ino(maxatm),norb(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      dimension dm(ndim,ndim),t(ndim,ndim),c(ndim,ndim),tnbo(ndim,ndim),
     +         dx(ndim,ndim),dy(ndim,ndim),dz(ndim,ndim),scr(ndim*ndim),
     +         index(ndim)
      dimension istr(14),couple(3)
c
      data tenten,small,zero,tenth,one,two/1.0d-10,1.0d-5,0.0d0,0.1d0,
     +                                     1.0d0,2.0d0/
      data toesu/4.803242e-10/
      data ihyph,iblnk/1h-,1h /
c
      debye = toesu / tenten
c
c  copy the nuclear charges into charge:
c
      if(alpha.or..not.open) then
        do 10 i = 1,natoms
          charge(i) = iznuc(i)
   10   continue
      end if
c
c  determine the number of occupied orbitals and make sure that the
c  occupied nlmos are at the beginning of the list:
c
      tot = zero
      do 20 i = 1,nbas
        tot = tot + dm(i,i)
        scr(i) = dm(i,i)
   20 continue
      nel = tot + tenth
      tot = nel
      nocc = nel
      if(.not.open) nocc = nocc/2 + mod(nocc,2)
c
      call rank(scr,nbas,ndim,index)
      do 30 i = 1,nocc
        if(index(i).gt.nocc) then
          write(lfnpr,1000)
          return
        end if
   30 continue
c
c  determine the occupancy factor:
c
      eta = two
      if(open) eta = one
c
c  compute the electronic contributions to the nbo bond dipole moments:
c
      call fetlmo(c)
      call fetnbo(tnbo)
      ii = 1
      call dipele(dx,c,tnbo,scr,eta,nocc,ii)
      if(ii.eq.0) return
      ii = 2
      call dipele(dy,c,tnbo,scr,eta,nocc,ii)
      if(ii.eq.0) return
      ii = 3
      call dipele(dz,c,tnbo,scr,eta,nocc,ii)
      if(ii.eq.0) return
c
c  add the nuclear contributions to these bond dipole moments:
c
      call NBOdipnuc(dx,dy,dz,scr,eta,nocc)
c
c  convert to debye:
c
      do 50 i = 1,nocc
        do 40 j = 1,nbas
          dx(j,i) = dx(j,i) * debye
          dy(j,i) = dy(j,i) * debye
          dz(j,i) = dz(j,i) * debye
   40   continue
   50 continue
c
c  print dipole analysis:
c
      xnbo  = zero
      ynbo  = zero
      znbo  = zero
      xnlmo = zero
      ynlmo = zero
      znlmo = zero
      do 100 i = 1,nocc
        if(i.eq.1) then
          if(alpha) write(lfnpr,1010)
          if(beta)  write(lfnpr,1020)
          if(.not.open) write(lfnpr,1030)
          write(lfnpr,1040) abs(dthr)
        else
          write(lfnpr,1050)
        end if
c
c  build the label for this nbo/nlmo:
c
        ib = ibxm(i)
        istr(1) = label(ib,1)
        istr(2) = label(ib,2)
        istr(3) = label(ib,3)
        do 70 j = 1,3
          j4 = 4 * j
          if(label(ib,j+3).eq.0) then
            do 60 k = j4-1,j4+2
              istr(k) = iblnk
   60       continue
          else
            if(j.ne.1) istr(j4-1) = ihyph
            istr(j4)   = nameat(iatno(label(ib,j+3)))
            call convrt(label(ib,j+3),istr(j4+1),istr(j4+2))
          end if
   70   continue
c
c  compute the nlmo bond dipole (the nbo bond dipoles are on the diagonal
c  of dx,dy,dz):
c
        x = zero
        y = zero
        z = zero
        do 80 j = 1,nbas
          x = x + dx(j,i)
          y = y + dy(j,i)
          z = z + dz(j,i)
   80   continue
c
        xnbo  = xnbo  + dx(i,i)
        ynbo  = ynbo  + dy(i,i)
        znbo  = znbo  + dz(i,i)
        xnlmo = xnlmo + x
        ynlmo = ynlmo + y
        znlmo = znlmo + z
c
c  compute the net dipole for these orbitals:
c
        tot = sqrt(dx(i,i)*dx(i,i) + dy(i,i)*dy(i,i) + dz(i,i)*dz(i,i))
        totnlm = sqrt(x*x + y*y + z*z)
c
        write(lfnpr,1060) i,(istr(j),j=1,14),x,y,z,totnlm,
     +                    dx(i,i),dy(i,i),dz(i,i),tot
c
c  print delocalization terms which are stronger than abs(dthr):
c
        icnt = 0
        do 90 j = 1,nbas
          if(j.ne.i) then
            tot = sqrt(dx(j,i)*dx(j,i) + dy(j,i)*dy(j,i)
     +                                  + dz(j,i)*dz(j,i))
            if(tot.gt.abs(dthr)) then
              icnt = icnt + 1
              index(icnt) = j
              scr(icnt) = tot
            end if
          end if
   90   continue
c
        do 95 j = 1,icnt
          do 94 k = 1,icnt-j
            if(scr(k+1)-scr(k).gt.small) then
              itemp      = index(k)
              index(k)   = index(k+1)
              index(k+1) = itemp
              temp       = scr(k)
              scr(k)     = scr(k+1)
              scr(k+1)   = temp
            end if
   94     continue
   95   continue
c
        do 96 jj = 1,icnt
          j = index(jj)
          write(lfnpr,1070) j,dx(j,i),dy(j,i),dz(j,i),scr(jj)
   96   continue
  100 continue
c
c  compute and print the correction for residual nuclear charges:
c
      if(.not.alpha) then
        call fecoor(scr)
        x = zero
        y = zero
        z = zero
        test = .false.
        do 110 i = 1,natoms
          if(abs(charge(i)).gt.small) test = .true.
          x = x + scr(3*i-2) * charge(i) * debye
          y = y + scr(3*i-1) * charge(i) * debye
          z = z + scr(3*i)   * charge(i) * debye
  110   continue
        if(test) then
          tot = sqrt(x*x + y*y + z*z)
          write(lfnpr,1080) x,y,z,tot,x,y,z,tot
          xnbo  = xnbo  + x
          ynbo  = ynbo  + y
          znbo  = znbo  + z
          xnlmo = xnlmo + x
          ynlmo = ynlmo + y
          znlmo = znlmo + z
        end if
      end if
c
c  print net dipole moments:
c
      tot = sqrt(xnbo*xnbo + ynbo*ynbo + znbo*znbo)
      totnlm = sqrt(xnlmo*xnlmo + ynlmo*ynlmo + znlmo*znlmo)
      write(lfnpr,1090) xnlmo,ynlmo,znlmo,totnlm,xnbo,ynbo,znbo,tot
c
c  compute and print the total delocalization correction:
c
      x = xnlmo - xnbo
      y = ynlmo - ynbo
      z = znlmo - znbo
      tot = sqrt(x*x + y*y + z*z)
      write(lfnpr,1100) x,y,z,tot
c
c  compute and print the nlmo coupling correction:
c
      test = .false.
      do 130 i = 1,nbas
        if(i.gt.nocc.and.abs(dm(i,i)).gt.small) test = .true.
        do 120 j = i+1,nbas
          if(abs(dm(j,i)).gt.small) test = .true.
  120   continue
  130 continue
      if(test) then
        tot = zero
        do 160 k = 1,3
          ii = k
          call fedxyz(dx,ii)
          call simtrs(dx,t,scr,ndim,nbas)
          couple(k) = zero
          do 150 i = 1,nbas
            if(i.le.nocc) then
              couple(k) = couple(k) + (eta - dm(i,i)) * dx(i,i)
            else
              couple(k) = couple(k) - dm(i,i) * dx(i,i)
            end if
            do 140 j = i+1,nbas
              couple(k) = couple(k) - two * dm(j,i) * dx(j,i)
  140       continue
  150     continue
          couple(k) = couple(k) * debye
          tot = tot + couple(k) * couple(k)
  160   continue
        tot = sqrt(tot)
        write(lfnpr,1110) xnlmo,ynlmo,znlmo,totnlm,xnlmo,ynlmo,znlmo,
     +                    totnlm,(couple(k),k=1,3),tot
        xnlmo = xnlmo + couple(1)
        ynlmo = ynlmo + couple(2)
        znlmo = znlmo + couple(3)
        totnlm = sqrt(xnlmo*xnlmo + ynlmo*ynlmo + znlmo*znlmo)
        if(alpha) write(lfnpr,1120) xnlmo,ynlmo,znlmo,totnlm
        if(beta)  write(lfnpr,1130) xnlmo,ynlmo,znlmo,totnlm
        if(.not.open) write(lfnpr,1140) xnlmo,ynlmo,znlmo,totnlm
      else
        if(alpha) write(lfnpr,1120) xnlmo,ynlmo,znlmo,totnlm,
     +                              xnlmo,ynlmo,znlmo,totnlm
        if(beta)  write(lfnpr,1130) xnlmo,ynlmo,znlmo,totnlm,
     +                              xnlmo,ynlmo,znlmo,totnlm
        if(.not.open) write(lfnpr,1140) xnlmo,ynlmo,znlmo,totnlm,
     +                                  xnlmo,ynlmo,znlmo,totnlm
      end if
c
c  save the alpha spin dipoles:
c
      if(alpha) then
        xdip = xnlmo
        ydip = ynlmo
        zdip = znlmo
      end if
c
c  print out the total dipole moment for open shell species:
c
      if(beta) then
        xnlmo  = xnlmo + xdip
        ynlmo  = ynlmo + ydip
        znlmo  = znlmo + zdip
        totnlm = sqrt(xnlmo*xnlmo + ynlmo*ynlmo + znlmo*znlmo)
        write(lfnpr,1140) xnlmo,ynlmo,znlmo,totnlm
      end if
      return
c
 1000 format(/1x,'the highest occupied nbos are not at the beginning ',
     + 'of the list.',/1x,'the dipole moment analysis is currently not',
     + ' set up to handle this.')
 1010 format(//1x,'dipole moment analysis, alpha spin:')
 1020 format(//1x,'dipole moment analysis, beta spin:')
 1030 format(//1x,'dipole moment analysis:')
 1040 format(/1x,'[print threshold: net dipole >',f5.2,' debye]',//1x,
     + '                                nlmo bond dipole            ',
     + 'nbo bond dipole',/1x,'                            ----------',
     + '---------------  ------------------------',/1x,'         ',
     + 'orbital              x     y     z   total      x     y     ',
     + 'z   total',/1x,79('='))
 1050 format(1x)
 1060 format(1x,i3,'. ',a2,a1,'(',i2,')',a2,3a1,a2,3a1,a2,2a1,1x,4f6.2,
     + 3x,4f6.2)
 1070 format(1x,44x,'deloc ',i3,':',4f6.2)
 1080 format(/1x,'  residual nuclear charge  ',4f6.2,'   ',4f6.2)
 1090 format(1x,'                           -----------------------',
     + '-----------------------------',/1x,'        net dipole moment',
     + '  ',4f6.2,'   ',4f6.2)
 1100 format(1x,'delocalization correction  ',24x,'   ',4f6.2,/1x,
     + '                           -----------------------------',
     + '-----------------------')
 1110 format(1x,'        net dipole moment  ',4f6.2,'   ',4f6.2,/1x,
     + ' nlmo coupling correction  ',4f6.2,/1x,'                  ',
     + '         -------------------------')
 1120 format(1x,'        alpha spin dipole  ',4f6.2,'   ',4f6.2)
 1130 format(1x,'         beta spin dipole  ',4f6.2,'   ',4f6.2)
 1140 format(1x,'      total dipole moment  ',4f6.2,'   ',4f6.2)
      end
c*****************************************************************************
      subroutine dipele(dxyz,c,t,scr,eta,nocc,index)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension dxyz(ndim,ndim),c(ndim,ndim),t(ndim,ndim),scr(ndim,ndim)
c
c
c  compute the electronic contribution for the x (index=1), y (=2),
c  and z (=3) components of the dipole:
c
c  get the ao dipole matrix and transform to the nbo basis:
c
      call fedxyz(dxyz,index)
      if(index.eq.0) return
      call simtrs(dxyz,t,scr,ndim,nbas)
c
c  compute the electronic contribution for doubly occupied, filled nbos:
c
      do 30 i = 1,nocc
        scr(i,i) = -eta * dxyz(i,i)
   30 continue
c
c  compute delocalization contributions for each filled nbo:
c
      do 60 i = 1,nocc
        do 50 j = 1,nbas
          if(j.ne.i) then
            scr(j,i) = c(j,i) * dxyz(i,i) - c(i,i) * dxyz(j,i)
            do 40 k = 1,nbas
                scr(j,i) = scr(j,i) - c(k,i) * dxyz(k,j)
   40       continue
            scr(j,i) = eta * c(j,i) * scr(j,i)
          end if
   50   continue
   60 continue
      call copy(scr,dxyz,ndim,nbas,nbas)
      return
      end
c*****************************************************************************
      subroutine NBOdipnuc(dx,dy,dz,atcoor,eta,nocc)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbdxyz/xdip,ydip,zdip,charge(maxatm)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      dimension dx(ndim,ndim),dy(ndim,ndim),dz(ndim,ndim),
     +          atcoor(3,natoms)
c
      data zero/0.0d0/
c
c  fetch the atomic coordinates:
c
      call fecoor(atcoor)
c
c  calculate the nuclear contributions to the dipole moment:
c
      do 20 i = 1,nocc
        nctr = mod(nbotyp(i),10)
        x = zero
        y = zero
        z = zero
        do 10 j = 1,nctr
          iat = label(ibxm(i),j+3)
          x   = x + atcoor(1,iat)
          y   = y + atcoor(2,iat)
          z   = z + atcoor(3,iat)
          charge(iat) = charge(iat) - eta/nctr
   10   continue
        x = eta * x / nctr
        y = eta * y / nctr
        z = eta * z / nctr
        dx(i,i) = dx(i,i) + x
        dy(i,i) = dy(i,i) + y
        dz(i,i) = dz(i,i) + z
   20 continue
      return
      end
c*****************************************************************************
c
c  routines called by sr nathyb, sr choose:
c
c      subroutine core(dm,t,borb,pol,q,hyb,bndocc,ibd,detail,lfnpr)
c      function iwprj(nctr)
c      subroutine deplet(dm,t,q,pol,borb,bndocc,nbd)
c      subroutine load(dm,iat1,iat2,iat3,blk,nb)
c      subroutine prjexp(borb,iat1,iat2,iat3,q,p,pk,hyb,va,vb,hybexp)
c      subroutine stash_nbo(borb,ibd,iat1,iat2,iat3,pol,q,hyb)
c      subroutine orthyb(q,s,ta,eval,c,ialarm,iflg)
c      subroutine frmprj(p,ia,q,nk,pk,vk,pi)
c      subroutine NBOaugmnt(p,blk,c,eval,dm,ta,borb,v,larc,ia,nocc,norb)
c      subroutine repol(dm,q,pol,blk,eval,c,nbd)
c      subroutine formt(t,q,pol)
c      subroutine cycles(iter,thresh,guide,bndocc,topo,icont)
c
c*****************************************************************************
      subroutine core(dm,t,borb,pol,q,hyb,bndocc,ibd,detail,lfnpr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  label core, valence, and rydberg nao's and deplete dm of the density
c  of the core orbitals
c
      logical detail,first
      parameter (maxatm = 750, maxbas = 4096)
      common/nbnao/naoctr(maxbas),naol(maxbas),ltyp(maxbas),
     +       iprin(maxbas)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ill(maxatm),
     +       iul(maxatm),iznuc(maxatm),iatcr(maxatm)
      dimension dm(ndim,ndim),t(ndim,ndim),borb(mxbo),pol(ndim,3),
     *  q(mxao,ndim),hyb(mxao),bndocc(ndim),icore(4),ival(4),iang(5)
      data zero,one/0.0d0,1.0d0/
      data iblk,icor/2h  ,2hcr/
      data ichcor,ichval,ichryd/3hcor,3hval,3hryd/
      data iang/1hs,1hp,1hd,1hf,1hg/
c
c  label nao's on each center:
c
      do 10 i = 1,nbas
        ltyp(i) = iryd
   10 continue
      iecp = 0
      do 110 nctr = 1,natoms
        call cortbl(nctr,icore,iecp)
        call valtbl(nctr,ival)
c
c  loop over s,p,d,f orbitals:
c
        do 100 l = 0,3
          ityp = iang(l+1)
          lnum = 2*l + 1
          if(icore(l+1).le.0) goto 50
c
c  label core orbitals:
c
          do 40 m = 1,icore(l+1)
            do 30 la = 1,lnum
              morb = 0
              occ = -1.0d0
              do 20 n = 1,nbas
                lm = naol(n)
                norb = lm/100
                il = iang(norb+1)
                na = mod(naol(n),50)
                if(naoctr(n).eq.nctr.and.il.eq.ityp.and.
     +            dm(n,n).gt.occ.and.ltyp(n).eq.iryd.and.
     +                                         la.eq.na) then
                      morb = n
                      occ = dm(n,n)
                end if
   20         continue
              if(morb.eq.0) then
                write(lfnpr,2500) ityp,nameat(iatno(nctr)),nctr,
     +                            (icore(i),i=1,4),m,la
                stop
              end if
              ltyp(morb) = ichcor
   30       continue
   40     continue
   50     continue
          if(ival(l+1).le.0) goto 90
c
c  label valence orbitals:
c
          do 80 m = 1,ival(l+1)
            do 70 la = 1,lnum
              morb = 0
              occ = -1.0d0
              do 60 n = 1,nbas
                lm = naol(n)
                norb = lm/100
                il = iang(norb+1)
                na = mod(naol(n),50)
                if(naoctr(n).eq.nctr.and.il.eq.ityp.and.
     +            dm(n,n).gt.occ.and.ltyp(n).eq.iryd.and.
     +                                         la.eq.na) then
                      morb = n
                      occ = dm(n,n)
                end if
   60         continue
              if(morb.eq.0) then
                write(lfnpr,2600) ityp,nameat(iatno(nctr)),nctr,
     +                            (ival(i),i=1,4),m,la
                stop
              end if
              ltyp(morb) = ichval
   70       continue
   80     continue
   90     continue
  100   continue
  110 continue
c
c  isolate core orbitals on all atoms, removing their density from the
c  density matrix:
c
      do 300 iat = 1,natoms
        nb = iul(iat) - ill(iat) + 1
        iac = 0
        first = .true.
        do 290 n = ill(iat),iul(iat)
          if(ltyp(n).eq.ichcor) then
            if(detail.and.first) then
              first = .false.
              write(lfnpr,1000) iat
            end if
            iac = iac + 1
            ibd = ibd + 1
            do 280 i = 1,nb
              borb(i) = zero
  280       continue
            borb(n-ill(iat)+1) = one
            call stash_nbo(borb,ibd,iat,0,0,pol,q,hyb)
            label(ibd,1) = icor
            label(ibd,2) = iblk
            label(ibd,3) = iac
            label(ibd,4) = iat
            bndocc(ibd)  = dm(n,n)
            if(detail) write(lfnpr,1010) iac,bndocc(ibd)
            if(detail) write(lfnpr,1020) (borb(i),i=1,nb)
            if(detail) write(lfnpr,1030) ibd,(label(ibd,i),i=1,3)
          end if
  290   continue
  300 continue
c
c  deplete the density matrix of cr orbitals:
c
      call deplet(dm,t,q,pol,borb,bndocc,ibd)
      return
c
 1000 format(/,1x,'search of dm block for core orbitals on atom:',i4)
 1010 format(6x,'eigenvector (',i2,') has occupancy ',f9.6,':')
 1020 format(11x,8f7.4)
 1030 format(11x,'*** nbo accepted: number',i3,'.   label:',a2,a1,
     + '(',i2,')')
 2500 format(/1x,'subroutine core could not find a ',a1,'-type ',
     + 'core orbital on atom ',a2,i2,'.',/,1x,'icore :',4i3,
     + '     m :',i3,'     la :',i3)
 2600 format(/1x,'subroutine core could not find a ',a1,'-type ',
     + 'valence orbital on atom ',a2,i2,'.',/,1x,'ival :',4i3,
     + '     m :',i3,'     la :',i3)
      end
c*****************************************************************************
      function iwprj(nctr)
c*****************************************************************************
      data nctr0/0/
c
c  return 0 (no projection wanted) if nctr is unchanged, 1 otherwise.
c
      iwprj=0
      if(nctr.eq.nctr0) return
       iwprj=1
       nctr0=nctr
       return
      end
c*****************************************************************************
      subroutine deplet(dm,t,q,pol,borb,bndocc,nbd)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  deplete density matrix dm of contribution from b.o.'borb':
c     dm ==> dm - occ*borb*borb(transpose).
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ill(maxatm),
     +       iul(maxatm),iznuc(maxatm),iatcr(maxatm)
      dimension dm(ndim,ndim),t(ndim,ndim),q(mxao,ndim),pol(ndim,3),
     *  borb(mxbo),bndocc(ndim)
      dimension iat(3)
c  restore dm from t
      do 10 j=1,nbas
        do 10 i=1,j
          dm(i,j)=t(i,j)
   10     dm(j,i)=dm(i,j)
c  main loop over nbd available bond orbitals:
      do 90 ibd=1,nbd
        occ=bndocc(ibd)
c  find atoms for b.o. #ibd
        nctr=0
        do 20 j=1,3
          iat(j)=label(ibd,j+3)
          if(iat(j).le.0) go to 30
          nctr=nctr+1
   20     continue
c  reconstruct borb for b.o. #ibd
   30   nelm=0
        do 40 ictr=1,nctr
          ia=iat(ictr)
          ihyb=iathy(ibd,ictr)+ill(ia)-1
          p=pol(ibd,ictr)
          nh=norbs(ia)
          do 40 ih=1,nh
            nelm=nelm+1
   40       borb(nelm)=p*q(ih,ihyb)
c  subtract occ*borb*borb(t) from dm
        nrow=0
        do 80 ictr=1,nctr
          ia=iat(ictr)
          iu=iul(ia)
          il=ill(ia)
          do 70 irow=il,iu
            nrow=nrow+1
            ncol=0
            do 60 jctr=1,nctr
              ja=iat(jctr)
              ju=iul(ja)
              jl=ill(ja)
              do 50 icol=jl,ju
                ncol=ncol+1
   50           dm(irow,icol)=dm(irow,icol)-occ*borb(nrow)*borb(ncol)
   60         continue
   70       continue
   80     continue
   90   continue
      return
      end
c*****************************************************************************
      subroutine NBOload(dm,iat1,iat2,iat3,blk,nb)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  zero the matrix 'blk' and load in atomic blocks of density
c  matrix 'dm' for the atoms listed in 'iat'
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ill(maxatm),
     +       iul(maxatm),iznuc(maxatm),iatcr(maxatm)
      dimension blk(mxbo,mxbo),dm(ndim,ndim),iat(3)
      data zero/0.0d0/
      iat(1)=iat1
      iat(2)=iat2
      iat(3)=iat3
c  zero 'blk'
      do 10 i=1,mxbo
        do 10 j=1,mxbo
   10     blk(i,j)=zero
      nrow=0
      ncol=0
      do 50 i=1,3
        ia=iat(i)
        if(ia.eq.0) go to 50
        iu=iul(ia)
        il=ill(ia)
        do 40 irow=il,iu
          nrow=nrow+1
          ncol=0
          do 30 j=1,3
            ja=iat(j)
            if(ja.eq.0) go to 30
            ju=iul(ja)
            jl=ill(ja)
            do 20 icol=jl,ju
              ncol=ncol+1
              blk(nrow,ncol)=dm(irow,icol)
   20         continue
   30       continue
   40     continue
   50   continue
      nb=nrow
      return
      end
c*****************************************************************************
      subroutine prjexp(borb,iat1,iat2,iat3,q,p,pk,hyb,va,vb,hybexp)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  determine how much of borb is composed of previously used hybrids.
c
c  return hybexp(i) = expectation value of hybrid "i" in borb over the
c                     projection operator p for the atom of the hybrid.
c
c  if no hybrid on atom i contributes to borb, hybexp(i) = zero.
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ill(maxatm),
     +       iul(maxatm),iznuc(maxatm),iatcr(maxatm)
c
      dimension iat(3),hyb(mxao),borb(mxbo),q(mxao,ndim),p(mxao,mxao),
     *  pk(mxao,mxao),va(mxao),vb(mxao),hybexp(3)
c
      data zero,one,eps/0.0d0,1.0d0,1.0d-5/
c
c  loop over atomic hybrids:
c
      iat(1) = iat1
      iat(2) = iat2
      iat(3) = iat3
      kmax   = 0
      do 50 i = 1,3
        hybexp(i) = zero
        ia = iat(i)
        if(ia.eq.0) go to 50
c
c  extract the ith atomic hybrid from borb:
c
        nu = iul(ia)
        nl = ill(ia)
        kmin = kmax + 1
        kmax = kmax + nu - nl + 1
        mj = 0
        do 10 k = kmin,kmax
          mj = mj + 1
          hyb(mj) = borb(k)
   10   continue
c
c  do hybrids from the ith atom contribute to borb?
c
        s = zero
        do 20 j = 1,mj
          s = s + hyb(j)**2
   20   continue
        if(s.lt.eps) go to 50
c
c  determine the projection expectation for this hybrid:
c
        nh = ino(ia)
        if(nh.eq.0) then
          hybexp(i) = one
        else
          call frmprj(p,ia,q,nh,pk,va,vb)
          pav = zero
          do 40 j = 1,mj
            do 30 k = 1,mj
              pav = pav + hyb(k) * p(k,j) * hyb(j)
   30       continue
   40     continue
          hybexp(i) = abs(pav) / s
        end if
   50 continue
      return
      end
c*****************************************************************************
      subroutine stash_nbo(borb,ibd,iat1,iat2,iat3,pol,q,hyb)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  decompose bond orbital 'borb' and store constituent hybrids in q
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ill(maxatm),
     +       iul(maxatm),iznuc(maxatm),iatcr(maxatm)
c
      dimension pol(ndim,3),q(mxao,ndim),borb(mxbo),iat(3),hyb(mxao)
c
      data zero/0.0d0/
c
c  loop over centers:
c
      iat(1) = iat1
      iat(2) = iat2
      iat(3) = iat3
      kmax   = 0
      do 40 i = 1,3
        ia = iat(i)
        if(ia.eq.0) go to 40
        nu = iul(ia)
        nl = ill(ia)
c
c  extract hybrid from bond orbital for atom ia:
c
        kmin = kmax + 1
        kmax = kmax + nu - nl + 1
        mj = 0
        do 10 k = kmin,kmax
          mj = mj + 1
          hyb(mj) = borb(k)
   10   continue
c
c  extract polarization coefficient, store in 'pol':
c
        psq = zero
        do 20 j = 1,mj
          psq = psq + hyb(j)**2
   20   continue
        p = sqrt(psq)
        pol(ibd,i) = p
c
c  one more hybrid for atom ia:
c
        ino(ia) = ino(ia) + 1
        ncol = ill(ia) + ino(ia) - 1
c
c  place normalized hybrid in appropriate block of q:
c
        nh = nu - nl + 1
        do 30 nrow = 1,nh
          if(p.eq.zero) then
            q(nrow,ncol) = zero
          else
            q(nrow,ncol) = hyb(nrow)/p
          end if
   30   continue
        iathy(ibd,i) = ino(ia)
   40 continue
      return
      end
c*****************************************************************************
      subroutine orthyb(q,s,ta,eval,c,ialarm,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  symmetric orthogonalization of available hybrids in q:
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       ilu(maxatm),iznuc(maxatm),iatcr(maxatm)
c
      dimension q(mxao,ndim),s(mxbo,mxbo),ta(mxao,mxao),
     *                    eval(mxbo),c(mxbo,mxbo)
c
      data zero,one/0.0d0,1.0d0/
      data toosml/1.0d-4/
c
c  toosml: "too small" -- threshold for an s matrix eigenvalue that is too
c   small and will cause numerical problems and is indicative of near-linear
c   dependency in the hybrids:
c
      ialarm = 0
      do 100 ia = 1,natoms
        il = ll(ia)
        nh = ino(ia)
        if(nh.gt.mxao) go to 800
        if(nh.le.1) go to 100
c
c  load ia-block of q into ta:
c
        do 10 j = 1,nh
          do 5 i = 1,mxao
            ta(i,j) = q(i,il+j-1)
    5     continue
   10   continue
c
c  form overlap matrix s = ta(transp)*ta:
c
        do 30 j = 1,nh
          do 25 i = j,nh
            temp = zero
            do 20 k = 1,mxao
              temp = temp + ta(k,i) * ta(k,j)
   20       continue
            s(i,j) = temp
            s(j,i) = temp
   25     continue
   30   continue
c
c  diagonalize overlap matrix:
c
        call NBOjacobi(nh,s,eval,c,mxbo,mxbo,0)
c
c  form inverse square root of s, store in s: (avoid numerical problems
c  of linear dependence ("too small" eigenvalues) by prescreening the
c  eigenvalues)
c
        do 40 i = 1,nh
          if(eval(i).lt.toosml) go to 810
          eval(i) = one / sqrt(eval(i))
   40   continue
        do 60 j = 1,nh
          do 55 i = j,nh
            temp = zero
            do 50 k = 1,nh
              temp = temp + eval(k) * c(i,k) * c(j,k)
   50       continue
            s(i,j) = temp
            s(j,i) = temp
   55     continue
   60   continue
c
c  form new tap=ta*s**(-1/2), store in c:
c
        do 80 j = 1,nh
          do 75 i = 1,mxao
            temp = zero
            do 70 k = 1,nh
              temp = temp + ta(i,k) * s(k,j)
   70       continue
            c(i,j) = temp
   75     continue
   80   continue
c
c  replace orthogonalized ta in array q:
c
        do 90 j = 1,nh
          do 85 i = 1,mxao
            q(i,il+j-1) = c(i,j)
   85     continue
   90   continue
  100 continue
c
c  symmetric orthogonalization complete:
c
      return
c
c  sound the alarm that too many hybrids were found on this atom:
c
  800 continue
      ialarm = ia
      if(iflg.eq.0) write(lfnpr,900) mxao,ia,nh
      return
c
c  sound the alarm that there are too many hybrids or that there is
c  linear dependency in the hybrids!!
c
  810 continue
      ialarm = ia
      if(iflg.eq.0) write(lfnpr,910) ia,eval(i),toosml
      return
c
  900 format(/4x,'only expected to find',i3,' hybrids on atom',i3,
     + ', but found',i3,'.')
  910 format(/4x,'the hybrids on atom',i3,' are linearly dependent.',
     + '  an eigenvalue (',f10.6,')',/4x,'of the hybrid overlap ',
     + 'matrix is too small (<',f7.5,').')
      end
c*****************************************************************************
      subroutine frmprj(p,ia,q,nk,pk,vk,pi)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  form projection matrix p to annihilate components of nk occupied
c  hybrids for atom ia.
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ill(maxatm),
     +       iul(maxatm),iznuc(maxatm),iatcr(maxatm)
c
      dimension p(mxao,mxao),vk(mxao),pi(mxao),q(mxao,ndim),
     *          pk(mxao,mxao)
c
      data zero,one/0.0d0,1.0d0/
c
c  initialize p = unit matrix:
c
      nb = norbs(ia)
      do 10 j = 1,nb
        do 5 i = 1,j
          p(i,j) = zero
          p(j,i) = zero
          if(i.eq.j) p(i,j) = one
    5   continue
   10 continue
c
c  form projection matrix p = p1*p2*...*pk*...*pnk to annihilate
c  components of the nk occupied hybrids vk:  pk = i - vk*vk(t).
c  loop over occupied hybrids vk, k = 1,...,nk:
c
      if(nk.le.0) return
c
c  extract occupied hybrid vk from array q:
c
      do 90 k = 1,nk
        icol = ill(ia) + k - 1
        do 30 i = 1,nb
          vk(i) = q(i,icol)
   30   continue
c
c  form projection matrix pk:
c
        do 40 j = 1,nb
          do 35 i = 1,j
            pk(i,j) = -vk(i) * vk(j)
            pk(j,i) = pk(i,j)
            if(i.eq.j) pk(i,j) = pk(i,j) + one
   35     continue
   40   continue
c
c  accumulate total projector p(k) = p(k-1)*pk:
c
        do 80 i = 1,nb
          do 60 j = 1,nb
            pi(j) = zero
            do 50 l = 1,nb
              pi(j) = pi(j) + p(i,l) * pk(l,j)
   50       continue
   60     continue
          do 70 j = 1,nb
            p(i,j) = pi(j)
   70     continue
   80   continue
   90 continue
      return
      end
c*****************************************************************************
      subroutine NBOaugmnt(p,blk,c,eval,dm,ta,borb,v,larc,ia,nocc,norb)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
      dimension p(mxao,mxao),ta(mxao,mxao),dm(ndim,ndim),c(mxbo,mxbo),
     + eval(mxbo),borb(mxbo),v(mxbo),blk(mxbo,mxbo),larc(nbas)
c
      data zero,eps,pt99,one/0.0d0,1.0d-5,0.99d0,1.0d0/
c
c  first, form set of "optimally diagonal" unit vectors to span rydberg space:
c
      naug = norb - nocc
      do 10 i = 1,norb
        larc(i) = 0
   10 continue
c
c  select projected nao unit vector from projector in p:
c
      do 300 iproj = 1,naug
        imax = 0
        prjmax = zero
        do 80 iao = 1,norb
          if(larc(iao).ne.0) go to 80
          proj = abs(p(iao,iao))
          if(proj.gt.pt99) go to 100
          if(proj.lt.prjmax) go to 80
          prjmax = proj
          imax = iao
   80   continue
        iao = imax
        proj = prjmax
  100   continue
c
c  put vector in borb, normalize, and save in c:
c
        sb = zero
        do 120 j = 1,norb
          b = p(iao,j)
          sb = sb + b * b
          borb(j) = b
  120   continue
        larc(iao) = iproj
        rnorm = one / sqrt(sb)
        do 130 j = 1,norb
          borb(j) = borb(j) * rnorm
  130   continue
        do 140 j = 1,norb
          c(j,iproj) = borb(j)
  140   continue
        if(iproj.eq.naug) go to 300
c
c  add borb to the projector in p:
c
        do 150 j = 1,norb
          do 145 i = 1,j
            ta(i,j) = -borb(i) * borb(j)
            ta(j,i) = ta(i,j)
            if(i.eq.j) ta(i,i) = ta(i,i) + one
  145     continue
  150   continue
        do 200 i = 1,norb
          do 180 j = 1,norb
            v(j) = zero
            do 170 l = 1,norb
              v(j) = v(j) + p(i,l) * ta(l,j)
  170       continue
  180     continue
          do 190 j = 1,norb
            p(i,j) = v(j)
  190     continue
  200   continue
  300 continue
c
c  put projected vectors in ta, ordered according to the nao parent:
c
      iaug = 0
      do 350 iao = 1,norb
        if(larc(iao).eq.0) go to 350
        iaug = iaug + 1
        itcol = larc(iao)
        do 330 j = 1,norb
          ta(j,iaug) = c(j,itcol)
  330   continue
  350 continue
c
c  load dm block for atom ia in blk:
c
      call NBOload(dm,ia,0,0,blk,norb)
c
c  form block of dm in rydberg basis in upper corner of blk:
c
      do 500 ib = 1,norb
        do 450 j = 1,naug
          sum = zero
          do 440 k = 1,norb
            sum = sum + blk(ib,k) * ta(k,j)
  440     continue
          v(j) = sum
  450   continue
        do 480 j = 1,naug
          blk(ib,j) = v(j)
  480   continue
  500 continue
      do 550 j = 1,naug
        do 520 i = 1,j
          sum = zero
          do 510 k = 1,norb
            sum = sum + ta(k,i) * blk(k,j)
  510     continue
          v(i) = sum
  520   continue
        do 530 i = 1,naug
          blk(i,j) = v(i)
  530   continue
  550 continue
      do 560 j = 1,naug
        jj = j - 1
        do 555 i = 1,jj
          blk(j,i) = blk(i,j)
  555   continue
  560 continue
c
c  diagonalize dm:
c
      call NBOjacobi(naug,blk,eval,c,mxbo,mxbo,1)
c
c  order eigenvectors by occupancy (within eps), form final rydberg vectors:
c
      do 570 i = 1,naug
        larc(i) = i
  570 continue
      naug1 = naug - 1
      do 620 i = 1,naug1
        i1 = i + 1
        do 610 j = i1,naug
          diff = eval(j) - eval(i)
          if(diff.lt.eps) go to 610
          temp = eval(i)
          eval(i) = eval(j)
          eval(j) = temp
          itemp = larc(i)
          larc(i) = larc(j)
          larc(j) = itemp
  610   continue
  620 continue
      do 700 j = 1,naug
        lj = larc(j)
        do 680 i = 1,norb
          sum = zero
          do 670 k = 1,naug
            sum = sum + ta(i,k) * c(k,lj)
  670     continue
          blk(i,j) = sum
  680   continue
  700 continue
      return
      end
c*****************************************************************************
      subroutine repol(dm,q,pol,blk,eval,c,nbd)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical print,first
c
c  diagonalize density matrix in basis of orthonormal hybrids for
c  each bond orbital to find new polarization coefficients.
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ill(maxatm),
     +       iul(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension dm(ndim,ndim),q(mxao,ndim),pol(ndim,3),
     *        blk(mxbo,mxbo),eval(mxbo),c(mxbo,mxbo)
c
      data zero,pt1,one,two/0.0d0,0.1d0,1.0d0,2.0d0/
      data lstar/1h*/
c
c  first, count number of bonds and 3c bonds:
c
      nbond = 0
      n3cb  = 0
      do 20 ib = 1,nbas
        if(label(ib,2).eq.lstar) go to 20
        if(label(ib,5).eq.0) go to 20
        nbond = nbond + 1
        if(label(ib,6).eq.0) go to 20
        n3cb = n3cb + 1
   20 continue
c
c  iab+1 is the number of the first antibond in the nbo list:
c
      iab = nbas - nbond - n3cb
c
      print = jprint(5).eq.1
      first = .true.
      apcoef = one / sqrt(two)
      do 200 ib = 1,nbd
        if(label(ib,2).eq.lstar) go to 200
        nctr = 1
        if(label(ib,5).gt.0) nctr = 2
        if(label(ib,6).gt.0) nctr = 3
        if(nctr.eq.1) go to 200
        if(iwapol.eq.0.or.nctr.eq.3) then
          do 120 i = 1,nctr
            ia  = label(ib,i+3)
            nhi = norbs(ia)
            do 115 j = 1,i
              ja  = label(ib,j+3)
              nhj = norbs(ja)
              dij = zero
              do 110 ir = 1,nhi
                irp = ill(ia)+ir-1
                cri = q(ir,ill(ia)+iathy(ib,i)-1)
                do 105 js = 1,nhj
                  jsp = ill(ja) + js - 1
                  csj = q(js,ill(ja)+iathy(ib,j)-1)
                  dij = dij+cri*csj*dm(irp,jsp)
  105           continue
  110         continue
              blk(i,j) = dij
              blk(j,i) = dij
  115       continue
  120     continue
c
c  diagonalize 'blk' and extract new polarization coefficients
c
          call NBOjacobi(nctr,blk,eval,c,mxbo,mxbo,0)
          call rank(eval,nctr,mxbo,larc)
c
c  make sure repolarization is not too drastic (take a look at the bond
c  orbital only):
c
          s = zero
          do 125 i = 1,nctr
            s = s + pol(ib,i) * c(i,larc(1))
  125     continue
          if(s.lt.pt1.and.nctr.eq.2) then
            if(first.and.print) write(lfnpr,*)
            first = .false.
            if(print) write(lfnpr,900) ib,s
            iab = iab + 1
            pol(iab,1) =  pol(ib,2)
            pol(iab,2) = -pol(ib,1)
          else
c
c  store the new polarization coefficients in pol:
c
            do 130 i = 1,nctr
              pol(ib,i) = c(i,larc(1))
  130       continue
            iab = iab + 1
            do 150 i = 1,nctr
              pol(iab,i) = c(i,larc(2))
  150       continue
            if(nctr.ne.3) go to 200
            iab = iab + 1
            do 160 i = 1,nctr
              pol(iab,i) = c(i,larc(3))
  160       continue
          end if
c
c  constrain bonds to be apolar, if requested (not set up to work with
c  3-center bonds):
c
        else
          pol(ib,1) = apcoef
          pol(ib,2) = apcoef
          iab = iab + 1
          pol(iab,1) = apcoef
          pol(iab,2) = -apcoef
        end if
  200 continue
      return
c
  900 format(1x,'WARNING: significant repolarization of nbo ',i3,' (s=',
     + f7.4,'); repol disabled.')
      end
c*****************************************************************************
      subroutine formt(t,q,pol)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      integer ul
c
c  construction of final transformation  matrix t from orthonormal
c  hybrids; rows of t labelled by naos, columns by nbos.
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),naoc(maxbas),naol(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),ibx(maxbas),iathy(maxbas,3)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       ul(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbnao/naoctr(maxbas),naoa(maxbas),ltyp(maxbas),
     +       iprin(maxbas)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      dimension t(ndim,ndim),q(mxao,ndim),pol(ndim,3)
c
      data lcr,llp,lbd,lstar,lry/2hcr,2hlp,2hbd,1h*,2hry/
      data zero/0.0d0/
c
c  reorder occupied nbos to put lone and core pairs last:
c
      ncr = 0
      nlp = 0
      nbds = 0
      do 10 nscan = 1,nbas
        if(label(nscan,2).eq.lstar) go to 10
        nbds = nbds + 1
        if(label(nscan,1).eq.llp) nlp = nlp + 1
        if(label(nscan,1).eq.lcr) ncr = ncr + 1
   10 continue
      icr = 0
      ilp = 0
      ibo = 0
      iab = 0
      do 40 ibd = 1,nbas
        if(label(ibd,2).eq.lstar) go to 30
        if(label(ibd,1).eq.lcr) go to 15
        if(label(ibd,1).eq.llp) go to 20
c
c  pair bonds:
c
        ibo = ibo + 1
        ibx(ibd) = ibo
        go to 40
c
c  core pairs:
c
   15   icr = icr + 1
        ibx(ibd) = icr + nbds - ncr - nlp
        go to 40
c
c  lone pairs and core pairs:
c
   20   ilp = ilp + 1
        ibx(ibd) = ilp + nbds - nlp
        go to 40
c
c  antibonds:
c
   30   iab = iab + 1
        ibx(ibd) = nbds + iab
   40 continue
c
c  zero transformation array:
c
      do 60 i = 1,nbas
        do 50 j = 1,nbas
          t(i,j) = zero
   50   continue
   60 continue
c
c  deposit final bond orbitals in matrix t:
c
      nbo = 0
      do 130 ibd = 1,nbas
        kbd = ibd
        if(label(ibd,2).ne.lstar) go to 100
        if(label(ibd,1).eq.lry) go to 100
        if(label(ibd,1).eq.llp) go to 100
c
c  antibond orbitals: search occupied orb. list to get proper hybrids.
c  search occupied bond orbs. for match with antibond atoms:
c
        do 90 k = 1,nbo
          do 70 i = 4,6
            if(label(k,i).ne.label(ibd,i)) go to 90
            if((label(k,3).le.0).and.(label(k,1).eq.lbd)) go to 90
   70     continue
c
c  negative irnk = label(k,3) means bond orbital was already used:
c
c  found match; set label(k,3)<0:
c
          kbd = k
          label(kbd,3) = -label(kbd,3)
          go to 100
   90   continue
c
c  couldn't find match...exit:
c
        write(lfnpr,9000) ibd,(label(ibd,jj),jj=1,6)
        stop
c
c  deposit bond orbitals in t matrix:
c
  100   continue
        do 120 i = 1,3
          ia = label(ibd,i+3)
          if(ia.eq.0) go to 120
          jl = ll(ia)
          ju = ul(ia)
          irow = 0
          icol = jl + iathy(kbd,i) - 1
          do 110 j = jl,ju
            irow = irow + 1
            jb = ibx(ibd)
  110       t(j,jb) = pol(ibd,i) * q(irow,icol)
  120     continue
        if(ibd.eq.kbd) nbo = ibd
  130   continue
c
c  restore label(i,3) > 0:
c
      do 140 i = 1,nbas
        if(label(i,3).lt.0) label(i,3) = -label(i,3)
  140   continue
c
c  set array ibxm: ibxm(ib) is the current location of b.o. # ib:
c
      do 150 ib = 1,nbas
        i = ibx(ib)
  150   ibxm(i) = ib
c
c  set phase of 1-center orbitals such that the largest s-type nao contribution
c  is positive:
c
      do 200 ib = 1,nbas
        nctr = 1
        do 160 il = 5,6
          if(label(ibxm(ib),il).ne.0) nctr = nctr + 1
  160   continue
        if(nctr.eq.1) then
          jmax = 0
          tmax = -1.0d0
          do 170 in = 1,nbas
            if(naoa(in).lt.100) then
              if(abs(t(in,ib)).gt.tmax) then
                jmax = in
                tmax = abs(t(in,ib))
              end if
            end if
  170     continue
          if(jmax.ne.0) then
            if(t(jmax,ib).lt.-1.0d-4) then
              do 180 in = 1,nbas
                t(in,ib) = -t(in,ib)
  180         continue
            end if
          end if
        end if
  200 continue
      return
c
 9000 format(/,1x,'can''t find bond/antibond match for nbo ',
     + i3,2x,a2,a1,'(',i2,')',3i4)
      end
c*****************************************************************************
      subroutine cycles(iter,thresh,guide,bndocc,topo,icont)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbthr/thrset,prjset,accthr,crtset,e2thr,athr,pthr,ethr,
     +             dthr,dlthr,chsthr
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbtopo/iorder(maxatm),jorder(maxatm),ntopo(maxatm,maxatm),
     +            n3ctr,i3ctr(10,3)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension guide(natoms,natoms),bndocc(ndim),topo(natoms,natoms)
c
      save jter,devmin,rhomin,best,rho,jbadl
c
      data lcr,lbd,l3c,llp,lstar/2hcr,2hbd,2h3c,2hlp,1h*/
      data small,zero,tenth,one,onept5,three,hundrd
     +              /1.0d-4,0.0d0,0.1d0,1.0d0,1.5d0,3.0d0,1.0d2/
      data devthr/0.1d0/
      data jtermx/9/
c
c  subroutine cycles controls the search for an acceptable resonance
c  structure:
c
c  arguments:
c        iter   : iteration counter incremented by the calling routine
c        thresh : occupancy threshold used in search for nbos
c        guide  : wiberg bond index
c        bndocc : array containing the nbo occupancies
c        topo   : bond index matrix to be compared with the wiberg indices
c        icont  : control flag (see below)
c
c  iter, guide, and bndocc are unaltered by this routine
c  thresh is modified by this routine, if the resonance keyword is selected
c  the topo matrix is constructed by this routine
c
c  control flag : (set by this routine)
c    icont =  2 : an acceptable lewis structure has been found, continue
c          =  1 : an acceptable lewis structure has been found, recompute the
c                 nbos for this structure
c          =  0 : bogus lewis structure, terminate search for nbos
c          = -1 : occupancy threshold and/or atom ordering have been
c                 changed.  repeat the search for nbos.
c
c  set atom permuting counter and minimum deviation in guide-topo:
c
      if(iter.eq.1) then
        jter  =  0
        icont = -1
      end if
      jter = jter + 1
      if(jter.eq.1) devmin = hundrd
c
c  the minimum occupancy threshold is 1.5e (0.5e for open shell):
c
      thrmin = onept5
      if(ispin.ne.0) thrmin = thrmin - one
c
c  determine the number of low occupancy orbitals in the lewis structure:
c
      ibadl  = 0
      ibadnl = 0
      sumlew = zero
      totele = zero
      do 10 i = 1,nbas
        totele = totele + bndocc(i)
        if(label(ibxm(i),2).ne.lstar) then
          sumlew = sumlew + bndocc(i)
          if(bndocc(i).lt.thresh) ibadl = ibadl + 1
        else
          if(bndocc(i).gt.abs(accthr)) ibadnl = ibadnl + 1
        end if
   10 continue
      nel    = totele + tenth
      totele = nel
      sum    = totele - sumlew
c
c  count the ecp electrons in the lewis structure:
c
      if(ipseud.ne.0) then
        mecp = 0
        do 20 iat = 1,natoms
          mecp = mecp + iatno(iat) - iznuc(iat)
   20   continue
        if(ispin.ne.0) mecp = mecp/2
        sumlew = sumlew + dble(mecp)
      end if
c
c  keep track of the best lewis structure found so far:
c
      if(jter.eq.1) rhomin = hundrd
      if(iter.eq.1.or.sum.lt.rho) then
        best  = thresh
        rho   = sum
        jbadl = ibadl
        do 25 i = 1,natoms
          jorder(i) = iorder(i)
   25   continue
      end if
c
c  count the number of core, lone pair, and bonding orbitals in this
c  resonance structure:
c
      mcr = 0
      mbd = 0
      m3c = 0
      mlp = 0
      do 30 i = 1,nbas
        if(label(i,1).eq.lcr.and.label(i,2).ne.lstar) mcr = mcr + 1
        if(label(i,1).eq.lbd.and.label(i,2).ne.lstar) mbd = mbd + 1
        if(label(i,1).eq.l3c.and.label(i,2).ne.lstar) m3c = m3c + 1
        if(label(i,1).eq.llp.and.label(i,2).ne.lstar) mlp = mlp + 1
   30 continue
c
c  build the topo matrix from lone pairs and 2- and 3-center bonds:
c
      do 50 i = 1,natoms
        do 40 j = 1,natoms
          topo(i,j) = zero
   40   continue
   50 continue
c
      do 60 i = 1,nbas
        ib   = ibxm(i)
        if(label(ib,1).ne.lcr.and.label(ib,2).ne.lstar) then
          iat1 = label(ib,4)
          nctr = 1
          iat2 = label(ib,5)
          if(iat2.ne.0) nctr = 2
          iat3 = label(ib,6)
          if(iat3.ne.0) nctr = 3
          if(nctr.eq.1) then
            topo(iat1,iat1) = topo(iat1,iat1) + one
          else if(nctr.eq.2) then
            topo(iat1,iat2) = topo(iat1,iat2) + one
            topo(iat2,iat1) = topo(iat2,iat1) + one
          else
            topo(iat1,iat2) = topo(iat1,iat2) + one/three
            topo(iat2,iat1) = topo(iat2,iat1) + one/three
            topo(iat1,iat3) = topo(iat1,iat3) + one/three
            topo(iat3,iat1) = topo(iat3,iat1) + one/three
            topo(iat2,iat3) = topo(iat2,iat3) + one/three
            topo(iat3,iat2) = topo(iat3,iat2) + one/three
          end if
        end if
   60 continue
c
c  determine the largest off-diagonal element of guide-topo:
c
      dev = zero
      do 80 j = 2,natoms
        do 70 i = 1,j-1
          if(guide(i,j)-topo(i,j).gt.dev) then
            dev = guide(i,j) - topo(i,j)
            iat = i
            jat = j
          end if
   70   continue
   80 continue
c
c  write info about this resonance structure:
c
      if(jprint(5).eq.1) then
        if(iter.eq.1) write(lfnpr,1000)
        write(lfnpr,1010) iter,jter,abs(thresh),sumlew,sum,mcr,mbd,
     +                    m3c,mlp,ibadl,ibadnl,dev
      end if
c
c  decide if this structure is acceptable:
c
c   *  accept the structure if choose was employed.
c   *  accept the structure if there is only one atom.
c   *  accept the structure if there are no low occupancy lewis orbitals
c      and dev is less than devthr.
c   *  accept the structure if the nobond option was selected.
c
c  good resonance structure:
c
      if(ibadl.eq.0.and.dev.lt.devthr) then
        if(jprint(5).eq.1) write(lfnpr,1020)
        if(jprint(5).eq.1) write(lfnpr,1030)
        icont = 2
        return
c
c  only one atom:
c
      else if(natoms.eq.1) then
        if(jprint(5).eq.1) write(lfnpr,1020)
        if(jprint(5).eq.1) write(lfnpr,1035)
        icont = 2
        return
c
c  directed nbo search:
c
      else if(ichoos.eq.1) then
        if(jprint(5).eq.1) write(lfnpr,1020)
        if(jprint(5).eq.1) write(lfnpr,1040)
        icont = 2
        return
c
c  nobond option selected:
c
      else if(jprint(10).ne.0) then
        if(jprint(5).eq.1) write(lfnpr,1020)
        if(jprint(5).eq.1) write(lfnpr,1050)
        icont = 2
        return
      end if
c
c  structure accepted due to the specification of the resonance keyword
c  or the occupancy threshold.  otherwise, accept the structure only if
c  there are no high occupancy lewis orbitals:
c
      if(icont.eq.1) then
        if(thrset.ge.zero) then
          if(jprint(5).eq.1) write(lfnpr,1020)
          if(jprint(5).eq.1) write(lfnpr,1060)
          icont = 2
        else if(jprint(14).ne.0) then
          if(jprint(5).eq.1) write(lfnpr,1020)
          if(jprint(5).eq.1) write(lfnpr,1070)
          icont = 2
        else if(ibadl.ne.0) then
          if(jprint(5).eq.1) write(lfnpr,1020)
          if(jprint(5).eq.1) write(lfnpr,1030)
          icont = 2
        end if
        return
      end if
c
c  if dev.eq.devmin.and.sum.eq.rhomin or too many atoms permutations,
c  stop atom permutations:
c
      if((abs(dev-devmin).lt.small.and.abs(sum-rhomin).lt.small).or.
     +                                 jter.ge.jtermx) then
c
c  if the occupancy threshold was set by the user, accept the best
c  structure:
c
        if(thrset.ge.zero) then
          if(abs(sum-rho).lt.small) then
            if(jprint(5).eq.1) write(lfnpr,1020)
            if(jprint(5).eq.1) write(lfnpr,1060)
            icont = 2
          else
            do 90 i = 1,natoms
              iorder(i) = jorder(i)
   90       continue
            jter  = 0
            icont = 1
          end if
c
c  if the resonance keyword was specified, pick the best resonance structure
c  for this occupancy threshold, and possibly decrement the threshold and
c  continue the search:
c
        else if(jprint(14).ne.0) then
          thresh = thresh - tenth
          if(thrmin-thresh.gt.small) then
            thresh = thresh + tenth
            if(abs(thresh-best).lt.small.and.abs(sum-rho).lt.small) then
              if(jprint(5).eq.1) write(lfnpr,1020)
              if(jprint(5).eq.1) write(lfnpr,1070)
              icont = 2
            else
              do 100 i = 1,natoms
                iorder(i) = jorder(i)
  100         continue
              thresh = best
              jter  = 0
              icont = 1
            end if
          else
            do 110 i = 1,natoms
              iorder(i) = jorder(i)
  110       continue
            jter  =  0
            icont = -1
          end if
c
c  otherwise, accept the best structure, but only if it had no lewis
c  orbitals with occupancy less than the occupancy threshold:
c
        else
          if(abs(sum-rho).lt.small.and.ibadl.eq.0) then
            if(jprint(5).eq.1) write(lfnpr,1020)
            if(jprint(5).eq.1) write(lfnpr,1030)
            icont = 2
          else if(jbadl.eq.0) then
            do 115 i = 1,natoms
              iorder(i) = jorder(i)
  115       continue
            jter  = 0
            icont = 1
          else
            if(jprint(5).eq.1) write(lfnpr,1020)
            if(jprint(5).eq.1) write(lfnpr,1080)
            icont = 0
          end if
        end if
        return
c
c  loop through atom ordering to find alternative resonance structures:
c
      else
        if(dev.lt.devmin) devmin = dev
        if(sum.lt.rhomin) rhomin = sum
        if(iat.eq.iorder(1).and.jat.eq.iorder(2)) then
          dev1 = zero
          do 130 j = 2,natoms
            do 120 i = 1,j-1
              if(guide(i,j)-topo(i,j).gt.dev1) then
                if((i.ne.iorder(1).and.j.ne.iorder(2)).and.
     +             (j.ne.iorder(1).and.i.ne.iorder(2))) then
                  dev1 = guide(i,j) - topo(i,j)
                  iat  = i
                  jat  = j
                end if
              end if
  120       continue
  130     continue
        end if
c
        jflg = 0
        do 140 i = natoms,2,-1
          if(iorder(i).eq.jat) jflg = 1
          if(jflg.eq.1) iorder(i) = iorder(i-1)
  140   continue
        iorder(1) = jat
        iflg = 0
        do 150 i = natoms,2,-1
          if(iorder(i).eq.iat) iflg = 1
          if(iflg.eq.1) iorder(i) = iorder(i-1)
  150   continue
        iorder(1) = iat
        icont = -1
      end if
      return
c
 1000 format(/1x,'                      Occupancies       Lewis ',
     + 'Structure    Low   High',/1x,'          occ.    --------',
     + '-----------  -----------------   occ   occ',/1x,' Cycle ',
     + '  Thresh.   Lewis   Non-Lewis     CR  BD  3C  LP    (L) ',
     + '  (NL)   Dev',/1x,77('='))
 1010 format(1x,i3,'(',i1,')',3x,f5.2,f12.5,f10.5,3x,4i4,2x,i4,3x,i4,
     + 3x,f5.2)
 1020 format(1x,77('-'))
 1030 format(/1x,'Structure accepted: No low occupancy Lewis orbitals')
 1035 format(/1x,'Structure accepted: only a single atom')
 1040 format(/1x,'Structure accepted: nbos selected via the $choose ',
     + 'keylist')
 1050 format(/1x,'Structure accepted: search for bonds prevented ',
     + 'by nobond keyword')
 1060 format(/1x,'Structure accepted: occupancy threshold (thresh) ',
     + 'set by user')
 1070 format(/1x,'Structure accepted: resonance keyword permits ',
     + 'strongly delocalized structure')
 1080 format(/1x,'Only strongly delocalized resonance structures can',
     + ' be found.',/1x,'The default procedure is to abort the nbo ',
     + 'search.  Include',/1x,'the resonance keyword in the $nbo ',
     + 'keylist to override this test.')
      end
c*****************************************************************************
c
c  routines called by sr nlmo:
c
c      subroutine symuni(tsym,a,cos,sin,ovlp,blk,eval,nrot,
c     +           niuniq,njuniq,ilist,jlist,noff,ioff,joff,ndim)
c      subroutine NBOsymort(s,t,blk,ndim,n,eval)
c
c*****************************************************************************
      subroutine symuni(tsym,a,cos,sin,ovlp,blk,eval,nrot,
     *           niuniq,njuniq,ilist,jlist,noff,ioff,joff,ndim)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      dimension tsym(nrot,nrot),a(ndim,ndim),blk(nrot,nrot),
     *  ovlp(nrot,nrot),eval(nrot)
      dimension ioff(noff),joff(noff),ilist(noff),jlist(noff)
      data zero,one/0.0d0,1.0d0/
      data eps/1.0d-6/
      do 40 i=1,nrot
        do 30 j=1,nrot
   30     tsym(i,j)=zero
   40     tsym(i,i)=one
      do 60 moff=1,noff
        iocc=ilist(moff)
        jemt=jlist(moff)
        do 60 i=1,nrot
          t=tsym(i,iocc)
          u=tsym(i,jemt)
          tsym(i,iocc)=cos*t-sin*u
   60     tsym(i,jemt)=sin*t+cos*u
c
c   average groups of the elements of the transformation matrix tsym
c    so that the symmetry inherent in the density matrix a is preserved,
c    making sure that the resulting "averaged" transformation is unitary
c
c
      jst=niuniq+1
      nrot=jst-1+njuniq
c
c   ave. diag. elem of occ orbs
      if(niuniq.eq.1) go to 140
      tot=zero
      do 100 i=1,niuniq
  100   tot=tot+tsym(i,i)
      ave=tot/niuniq
      do 110 i=1,niuniq
  110   tsym(i,i)=ave
c
c   ave. diag. elem of empty orbs
  140 if(njuniq.eq.1) go to 180
      tot=zero
      do 150 j=jst,nrot
  150   tot=tot+tsym(j,j)
      ave=tot/njuniq
      do 160 j=jst,nrot
  160   tsym(j,j)=ave
c
c  zero offdiag elem betw occ orbs:
  180 if(niuniq.eq.1) go to 240
      do 220 i=2,niuniq
        do 220 j=1,i
          if(i.eq.j) go to 220
          tsym(i,j)=zero
          tsym(j,i)=zero
  220     continue
c
c  zero offdiag elem betw empty orbs:
  240 if(njuniq.eq.1) go to 280
      jst2=jst+1
      do 270 i=jst2,nrot
        do 270 j=jst,i
          if(i.eq.j) go to 270
          tsym(i,j)=zero
          tsym(j,i)=zero
  270     continue
c
c  ave. offdiag elem betw occ and empty orbs (pivoted elements only):
  280 continue
      tot=zero
      do 310 moff=1,noff
        ii=ilist(moff)
        jj=jlist(moff)
  310   tot=tot+abs(tsym(ii,jj))+abs(tsym(jj,ii))
      noff2=noff*2
      ave=tot/noff2
      do 330 moff=1,noff
        ii=ilist(moff)
        jj=jlist(moff)
        tsym(ii,jj)=-ave
  330   tsym(jj,ii)= ave
c
c  now zero the non-pivoted elements:
      do 450 i=1,niuniq
        do 440 j=jst,nrot
          do 420 moff=1,noff
            if(i.eq.ilist(moff).and.j.eq.jlist(moff)) go to 440
  420       continue
          tsym(i,j)= zero
          tsym(j,i)= zero
  440     continue
  450   continue
c
c  renormalize vectors:
      do 700 j=1,nrot
        tot=zero
        do 650 i=1,nrot
  650     tot=tot+tsym(i,j)*tsym(i,j)
        rnorm=sqrt(tot)
        if(rnorm.gt.eps) go to 680
          write(lfnpr,2880) nrot,tot,eps,rnorm
 2880     format('nrot,tot,eps,rnorm:',i3,3f14.9)
          call altout(tsym,nrot,nrot,nrot,nrot)
          stop
  680   continue
        do 690 i=1,nrot
  690     tsym(i,j)=tsym(i,j)/rnorm
  700   continue
c
c  now, make sure the signs are correct:
      do 800 moff=1,noff
        i=ioff(moff)
        j=joff(moff)
        if(a(i,j).gt.zero) go to 800
          ii=ilist(moff)
          jj=jlist(moff)
          tsym(ii,jj)=-tsym(ii,jj)
          tsym(jj,ii)=-tsym(jj,ii)
  800   continue
c
c  finally, the crucial step of symmetrically orthogonalizing the vectors
c   so that the transformation is unitary:
      call NBOsymort(ovlp,tsym,blk,nrot,nrot,eval)
      return
c
      end
c*****************************************************************************
      subroutine NBOsymort(s,t,blk,ndim,n,eval)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c******************************************************************
c
c   NBOsymort: symmetric orthogonalization subroutine
c
c   s:           full overlap matrix               (destroyed!)
c   t:           vectors to be orthoged.
c   n:           number of vectors
c
c   note:    blk and bigblk share the same storage but are
c               dimensioned differently.
c            the same applies for s and sblk.
c
c******************************************************************
      dimension s(n,n),t(ndim,ndim),blk(n,n),eval(n)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      data zero,one/0.0d0,1.0d0/
c
c  important constants:
c           diagth          threshold for matrix diagonalization used in
c                            subroutine jacobi.  in jacobi, this constant
c                            is called "doneth".
c           danger          criterion for deciding that the job should be
c                            aborted due to numerical problems caused by near
c                            linear dependencies in the basis set.  all
c                            eigenvalues of the weighted overlap matrix must
c                            be greater than diagth*danger.
c
      data diagth,danger/1.0d-12,1.0d3/
c
c  form the inverse sqrt of the overlap matrix of the vectors:
      do 70 i=1,n
        do 70 j=1,n
          sij=zero
          do 40 k=1,n
   40       sij=sij+t(k,i)*t(k,j)
   70     s(i,j)=sij
      call NBOjacobi(n,s,eval,blk,n,n,0)
      smlest=one
      toosml=diagth*danger
      do 150 i=1,n
        eigenv=eval(i)
        if(eigenv.lt.toosml) go to 900
        eval(i)=one/sqrt(eigenv)
        if(eigenv.lt.smlest) smlest=eigenv
  150  continue
      do 170 i=1,n
        do 170 j=1,i
          sij=zero
          do 160 k=1,n
  160       sij=sij+eval(k)*blk(i,k)*blk(j,k)
          s(i,j)=sij
  170     s(j,i)=sij
c
c  s now contains the -0.5 power of the overlap matrix,
c   and is the orthog. transform that we want.
c   now, form the total transformation:
      do 210 i=1,n
        do 200 j=1,n
          eval(j)=zero
          do 200 k=1,n
  200       eval(j)=eval(j)+t(i,k)*s(k,j)
      do 210 j=1,n
  210   t(i,j)=eval(j)
      return
c
  900 write(lfnpr,910) eigenv,toosml
  910 format(/1x,'an eigenvalue of the overlap matrix of the ',
     *   'symmetrized jacobi transf. ',
     *   'matrix of ',e13.5,' has been found.'/1x,
     *   'this is lower than the allowed threshold of ',e13.5)
      stop
      end
c*****************************************************************************
c
c  nbo energetic analysis routines:
c
c      subroutine nboean(a,memory,nboopt,idone)
c      subroutine nbodel(a,memory,idone)
c      subroutine delete(f,trf,ndim,idel,len,itype,ndel,ntrunc,done,
c     +                  ispin)
c      subroutine newdm(dm,u,eig,ndim,idel,len,ndel,itype,nmoocc,ispin)
c      subroutine rnkeig(rank,eig,n,ndim,arcrnk)
c      subroutine simltr(n,ndim,f,u,r,s,kntrol)
c
c*****************************************************************************
      subroutine nboean(a,memory,nboopt,idone)
c*****************************************************************************
c
c     nboean: controller subroutine to do nbo energetic analysis
c               by fock matrix deletion method
c
c       a(memory) is scratch storage
c
c       nboopt(1) = 2       read in next deletion and form new dm
c                 = 3       compute energy change for this deletion
c
c       set idone to 1 if no deletions are found:
c
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical error,new,seq
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension a(memory),nboopt(10)
c
      data thrneg/-1.0d-3/
      data one,aukcal,evkcal/1.0d0,627.51,23.061/
c
c  open the old nbo daf:
c
      new = .false.
      call nbopen(new,error)
      if(error) then
        idone = 1
        return
      end if
      call feinfo(a,iswean)
c
c  if nboopt(1) = 3,  compute the energy of deletion:
c
      if(nboopt(1).eq.3) then
        call fee0(edel,etot)
        echang = edel - etot
        if(munit.eq.0) then
          conv = aukcal
        else if(munit.eq.1) then
          conv = evkcal
        else
          conv = one
        end if
        ekcal = echang * conv
        if(ekcal.lt.thrneg) write(lfnpr,2130)
        if(munit.eq.0) then
          write(lfnpr,2100) edel,etot,echang,ekcal
        else if(munit.eq.1) then
          write(lfnpr,2110) edel,etot,echang,ekcal
        else
          write(lfnpr,2120) edel,etot,echang,ekcal
        end if
        idone = 0
        seq = .false.
        call nbclos(seq)
        return
      end if
c
c  perform the nbo energetic analysis:
c
c  if iswean is set to 1, search for the $del keylist:
c
      if(iswean.eq.1) then
        call delinp(nboopt,idone)
        if(idone.eq.1) goto 900
      else if(nboopt(10).gt.80) then
        call strtin(lfnin)
      end if
c
c  rohf, mcscf, ci, and auhf wave functions are not acceptable:
c
      if(rohf.or.mcscf.or.ci.or.auhf) then
        idone = 1
        goto 900
      end if
c
      ispin = 0
      if(uhf) ispin = 2
      alpha = .false.
      beta  = .false.
      if(uhf) alpha = .true.
      call nbodel(a,memory,idone)
      if(idone.eq.1) goto 900
c
      if(uhf) then
        ispin = -2
        alpha = .false.
        beta  = .true.
        call nbodel(a,memory,idone)
      end if
c
      write(lfnpr,3000)
      seq = .false.
      call nbclos(seq)
      return
c
  900 continue
      seq = .false.
      call nbclos(seq)
      return
c
 2100 format(1x,78('-'),/,3x,
     +'energy of deletion : ',f20.9,/,3x,
     +'  total scf energy : ',f20.9,/,3x,
     +'                       -------------------',/,3x,
     +'     energy change : ',f17.6,' a.u.,   ',f13.3,' kcal/mol'/
     +1x,78('-'))
 2110 format(1x,78('-'),/,3x,
     +'energy of deletion : ',f20.9,/,3x,
     +'  total scf energy : ',f20.9,/,3x,
     +'                       -------------------',/,3x,
     +'     energy change : ',f17.6,' e.v.,   ',f13.3,' kcal/mol'/
     +1x,78('-'))
 2120 format(1x,78('-'),/,3x,
     +'energy of deletion : ',f13.3,/,3x,
     +'  total scf energy : ',f13.3,/,3x,
     +'                       -------------------',/,3x,
     +'     energy change : ',f13.3,' kcal/mol,   ',f13.3,' kcal/mol'/
     +1x,78('-'))
 2130 format(/,6x,
     +'***** WARNING *****  the variational principle has been',/,5x,
     +'  violated and the above deletion energy is invalid!!',//,5x,
     +'probable cause:  a deletion was attempted that did not ',/,5x,
     +'have as high symmetry as was employed in the integral',/,5x,
     +'and scf computation.  remedy:  redo computation without',/,5x,
     +'symmetry if this non-symmetry-conserving deletion is still',/,5x,
     +'desired.')
 3000 format(/1x,
     +'next step:  evaluate the energy of the new density matrix',/,1x,
     +'            that has been constructed from the deleted nbo',/,1x,
     +'            fock matrix by doing one scf cycle.'/)
      end
c*****************************************************************************
      subroutine nbodel(a,memory,idone)
c*****************************************************************************
c
c     nbodel: subroutine to delete bond orbital fock matrix elements for
c              a particular spin case:
c                ispin = 0     closed shell
c                        2     alpha spin
c                       -2     beta  spin
c
c     idone is set equal to 1 if there are no more deletions,
c                           0 otherwise.
c
c     a(memory) is scratch storage
c
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical done
      dimension a(memory),ich(3,2),inam(3),isp(3)
c
c  nbo common blocks:
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       iatno(maxbas),ibxm(maxbas),iscr1(2*maxbas),iscr2(2*maxbas)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      data lbd/2hbd/,l3c/2h3c/,lblnk2/2h  /,lblnk1/1h /,lhyp/1h-/
c
c   fnbo  :  nbo fock matrix (triangular)
c   trf   :  truncated fock matrix (square)
c   eigvr :  eigenvectors of ftrunc
c   dmnew :  new ao dm (from truncation) -- triangular
c   occ   :  occupation vector of bond orbitals
c   occnew:  occupation vector of bond orbitals, after deletion
c   tnbo  :  ao to nbo transformation matrix
c   scr   :  scratch vector
c
c  set up storage space:
c
c   a(n1):  occ
c   a(n2):  occnew
c   a(n3):  tnbo
c   a(n4):  fnbo, eigvr
c   a(n5):  scr, trf, dmnew
c   a(n6):  scr
c   a(n7):  idel
c
      nsq  = ndim*ndim
      n1   = 1
      n2   = n1 + ndim
      n3   = n2 + ndim
      n4   = n3 + nsq
      n5   = n4 + nsq
      n6   = n5 + nsq
      n7   = n6 + ndim
      nend = n7 + nsq
      if(nend.gt.memory) go to 950
      call fenbo(a(n3),a(n1),a(n5),nelec)
      call fefnbo(a(n4))
c
c  delete requested fock matrix elements, forming truncated fock matrix
c             in trf
c
c   idel  :  list of deleted orbitals, elements, or blocks
c   itype :  type of deletion: 1 for orbitals
c                              2 for individual matrix elements
c			       3 for zeroing intersection between two sets
c			                                 of orbitals
c                              4 for entire matrix blocks
c   ndel  :  number of orbitals, elements or blocks to be deleted
c
      call delete(a(n4),a(n5),ndim,a(n7),nsq,itype,ndel,ntrunc,done,
     +            ispin)
c
c  if no more deletions, exit program
c
      if(done) go to 900
c  diagonalize truncated fock matrix in trf
c
      call NBOjacobi(ntrunc,a(n5),a(n2),a(n4),ndim,ndim,0)
c
c  construct new density matrix in dm from eigenvectors of trf,
c   in nbo basis:
c   a(n2):  eigenvalues of trf        (entering)
c   a(n2):  new nbo orbital occupancies  (exiting)
c
      nmoocc=nelec
      if(ispin.eq.0) nmoocc=nelec/2
      call newdm(a(n5),a(n4),a(n2),ndim,a(n7),nsq,ndel,itype,nmoocc,
     +           ispin)
c
c  take transpose of t so that it can transform the density matrix
c    from the nbo basis to the unsymmetrized ao basis:
c
      call NBOtransp(a(n3),ndim,ndim)
      call simltr(ndim,ndim,a(n5),a(n3),a(n4),a(n6),1)
      call svnewd(a(n5))
c
      write(lfnpr,2200)
      write(lfnpr,2700)
      do 500 ibas=1,ndim
            ib=ibxm(ibas)
            lbl=label(ib,1)
            nctr=1
            if(lbl.eq.lbd) nctr=2
            if(lbl.eq.l3c) nctr=3
            do 350 i=1,3
              iat=label(ib,i+3)
              call convrt(iat,ich(i,1),ich(i,2))
              inam(i)=lblnk2
              if(iat.gt.0) inam(i)=nameat(iatno(iat))
              isp(i)=lhyp
              if(i.ge.nctr) isp(i)=lblnk1
  350         continue
        i=n1-1+ibas
        ii=n2-1+ibas
        occchg=a(ii)-a(i)
        write(lfnpr,2800) ibas,(label(ib,k),k=1,3),
     *         (inam(k),ich(k,1),ich(k,2),isp(k),k=1,3),
     *         a(i),a(ii),occchg
  500   continue
      idone=0
      return
c
  900 continue
      idone=1
      return
c
  950 continue
      write(lfnpr,9500) nend,memory
      idone=1
      return
c
 2200 format(/1x,'occupations of bond orbitals:')
 2700 format(/7x,'orbital',19x,'no deletions   this deletion   change',
     + /,1x,78('-'))
 2800 format(1x,i3,'. ',a2,a1,'(',i2,')',3(a2,3a1),
     *       9x,f7.5,8x,f7.5,3x,f8.5)
 9500 format(/1x,'insufficient memory in subroutine nbodel:',
     *      /5x,'memory needed: ',i10,'   memory available: ',i10,
     *      /1x,'deletions halted!')
      end
c*****************************************************************************
      subroutine delete(f,trf,ndim,idel,len,itype,ndel,ntrunc,done,
     +                  ispin)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical error,done,equal,equal2
      logical donor,accptr,list1,list2
      dimension keywd(6),f(1),trf(ndim,ndim),idel(len)
      dimension lorb(3),lele(3),lblo(3),ldel(3),lzero(4),lsame(4),
     *          lend(3),ldestr(6),ldeloc(5),lnostr(6),latom(4),
     *          lnogem(5),lnovic(5),lalt(4)
      parameter (maxatm = 750, maxbas = 4096)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       iatno(maxbas),ibxm(maxbas),iscr1(2*maxbas),iscr2(2*maxbas)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      data zero/0.0d0/,istar/1h*/
      data ldel/1hd,1he,1hl/,lzero/1hz,1he,1hr,1ho/,lend/1he,1hn,1hd/
      data lalpha,lbeta/1ha,1hb/,lsame/1hs,1ha,1hm,1he/
      data lorb,lele,lblo/1ho,1hr,1hb,1he,1hl,1he,1hb,1hl,1ho/
      data ldestr/1hd,1he,1hs,1ht,1ha,1hr/
      data lnostr/1hn,1ho,1hs,1ht,1ha,1hr/
      data ldeloc/1hd,1he,1hl,1ho,1hc/,latom/1ha,1ht,1ho,1hm/
      data lnovic/1hn,1ho,1hv,1hi,1hc/,lnogem/1hn,1ho,1hg,1he,1hm/
      data lalt/1h$,1he,1hn,1hd/
      data lg,lv/1hg,1hv/
c
c   this subroutine is called at the start of each deletion and reads
c    in from lfnin the instructions for this deletion
c
c   ntrunc= dimension of fock matrix after deletions:
      ntrunc=ndim
      write(lfnpr,8700)
c  count up number of molecular units, nchemu:
      nchemu=0
      do 1 i=1,ndim
        nunit=nbouni(i)
        if(nunit.gt.nchemu) nchemu=nunit
    1   continue
      if(ispin.eq.0) go to 10
c  if open shell, look for first letter of "alpha" or "beta" keyword:
        leng=3
        call hfld(keywd,leng,done)
        if(equal(keywd,lend,3)) done=.true.
        if(equal2(keywd,lalt,3)) done=.true.
        if(done) return
        if((ispin.eq.2).and.(keywd(1).ne.lalpha)) go to 9300
        if((ispin.eq.-2).and.(keywd(1).ne.lbeta)) go to 9400
        if(ispin.eq.2) write(lfnpr,8100)
        if(ispin.eq.-2) write(lfnpr,8200)
c  search for first 3 letters of "delete", "zero", "same", "destar",
c    "nostar", "nogem", "novic", or an end mark '**':
   10 continue
      leng=3
      call hfld(keywd,leng,done)
      if(equal(keywd,lend,3)) done=.true.
      if(equal2(keywd,lalt,3)) done=.true.
      if(done) return
c  if beta deletions are the same as the alpha deletions already read in,
c    skip to 100:
      if((ispin.eq.-2).and.equal(keywd,lsame,3)) go to 100
      if(equal(keywd,lzero,3)) go to 600
      if(equal(keywd,lnovic,3)) go to 3000
      if(equal(keywd,lnogem,3)) go to 3010
      if(equal(keywd,ldestr,3)) go to 5000
      if(equal(keywd,lnostr,3)) go to 5500
      if(.not.equal(keywd,ldel,3)) go to 9000
c  read in number of items to delete, ndel:
      call ifld(ndel,error)
      if(error) go to 9100
c  read in type of deletion and determine if it is orbital, element, or block:
c   (itype stores the deletion type)
      call hfld(keywd,leng,done)
      if(leng.lt.3) go to 9200
      if(.not.equal(keywd,lorb,3)) go to 20
      itype=1
      go to 80
  20  if(.not.equal(keywd,lele,3)) go to 30
      itype=2
      go to 80
  30  if(.not.equal(keywd,lblo,3)) go to 9200
      itype=4
  80  continue
c  nread=number of numbers that must be read
      nread=ndel*itype
c  read in orbitals,elements, or blocks:
      do 90 i=1,nread
        call ifld(idel(i),error)
        if(error) go to 9500
   90   continue
c
  100 continue
      if(itype.ne.1) go to 200
c   delete ndel orbitals, adjusting ntrunc accordingly:
        ntrunc=ndim-ndel
c   order the orbital numbers:
        call NBOorder(iscr1,idel,ndel,ndim,iscr2)
        write(lfnpr,8610) (idel(i),i=1,ndel)
c   fill trf with truncated fock matrix, deleting requested orbitals:
        iff=0
        iout=1
        ii=0
        do 140 i=1,ndim
          if(iout.gt.ndel) go to 110
          if(i.ne.idel(iout)) go to 110
            iff=iff+i
            iout=iout+1
            go to 140
  110     continue
            ii=ii+1
            jout=1
            jj=0
            do 130 j=1,i
              if(jout.gt.ndel) go to 120
              if(j.ne.idel(jout)) go to 120
                iff=iff+1
                jout=jout+1
                go to 130
  120         continue
                jj=jj+1
                iff=iff+1
                trf(ii,jj)=f(iff)
                trf(jj,ii)=f(iff)
  130         continue
  140     continue
        return
  200 continue
c  element or block deletions: start by filling trf with full nbo fock matrix:
      ii=0
      do 210 i=1,ndim
        do 210 j=1,i
          ii=ii+1
          trf(i,j)=f(ii)
          trf(j,i)=f(ii)
  210     continue
      if(itype.ne.2) go to 300
c  zero requested matrix elements:
        ndel2=ndel*2
        write(lfnpr,8620) (idel(i),i=1,ndel2)
        do 240 i=1,ndel
          i2=2*i
          id=idel(i2-1)
          jd=idel(i2)
          trf(id,jd)=zero
          trf(jd,id)=zero
  240     continue
        return
  300 continue
      if(itype.ne.4) stop
c  zero requested matrix blocks:
        do 400 id=1,ndel
          idst=(id-1)*4
          j1=idel(idst+1)
          j2=idel(idst+2)
          i1=idel(idst+3)
          i2=idel(idst+4)
          if(j1.le.j2) go to 320
            idel(idst+2)=j1
            idel(idst+1)=j2
            j1=idel(idst+1)
            j2=idel(idst+2)
  320     if(i1.le.i2) go to 330
            idel(idst+4)=i1
            idel(idst+3)=i2
            i1=idel(idst+3)
            i2=idel(idst+4)
  330     do 380 i=i1,i2
            do 380 j=j1,j2
c  skip diagonal elements:
              if(i.eq.j) go to 380
              trf(i,j)=zero
              trf(j,i)=zero
  380       continue
  400     continue
        ndel4=ndel*4
        write(lfnpr,8640) (idel(i),i=1,ndel4)
      return
c  delete intersection in fock matrix between pairs of sets of orbitals:
  600 itype=3
c  start by filling trf with full nbo fock matrix:
      ii=0
      do 610 i=1,ndim
        do 610 j=1,i
          ii=ii+1
          trf(i,j)=f(ii)
          trf(j,i)=f(ii)
  610     continue
c  read in number of pairs of sets of orbitals, ndel:
      call ifld(ndel,error)
      if(error) go to 9500
      leng=5
c  check the next word to see if it is "delocalization" instead of "block":
c  (if so, the block will be specified by molecular units instead of by blocks)
      call hfld(keywd,leng,done)
      if(equal(keywd,ldeloc,5)) go to 1000
c  check the word to see if it is "atom" instead of "block":
c   (if so, the block will be specified by orbitals on groups of atoms)
      if(equal(keywd,latom,4)) go to 1200
      nstart=0
      do 800 k=1,ndel
c  read in the number of orbitals in each set of the pair, nset1 and nset2:
c    (skip the 'by' between nset1 and nset2)
        call ifld(nset1,error)
        if(error) go to 9500
        call hfld(keywd,leng,done)
        call ifld(nset2,error)
        if(error) go to 9500
        nstart=nstart+2
        idel(nstart-1)=nset1
        idel(nstart)=nset2
c  read in the orbitals of both sets
        ntot=nset1+nset2
        do 620 i=1,ntot
          call ifld(idel(nstart+i),error)
          if(error) go to 9500
  620     continue
c  now, zero all intersecting elements between the two sets:
        nstrt2=nstart+nset1
        do 700 i=1,nset1
          id=idel(nstart+i)
          do 700 j=1,nset2
            jd=idel(nstrt2+j)
            if(id.eq.jd) go to 700
            trf(id,jd)=zero
            trf(jd,id)=zero
  700       continue
        nstart=nstart+ntot
  800   continue
      go to 4000
c
c  zeroing of delocalization within or between molecular units.
c
c   use the nbo molecular unit (nbouni) and nbo type (nbotyp) lists.
 1000 continue
      nstart=0
      do 1100 k=1,ndel
c  skip the next word ("from"):
        call hfld(keywd,leng,done)
c  read in the number of the first molecular unit, iunit1:
        call ifld(iunit1,error)
        if(error) go to 9500
c  skip the "to" and read in iunit2:
        call hfld(keywd,leng,done)
        call ifld(iunit2,error)
        if(error) go to 9500
        write(lfnpr,8300) iunit1,iunit2
        nstart=nstart+2
c  find all of the nonstar (core/"lone pair"/bond) nbos on unit iunit1:
        nset1=0
        do 1020 ibas=1,ndim
          if(nbouni(ibas).ne.iunit1) go to 1020
          if(nbotyp(ibas).gt.20) go to 1020
          nset1=nset1+1
          idel(nstart+nset1)=ibas
 1020     continue
        idel(nstart-1)=nset1
c  find all of the star (rydberg/antibond) nbos on unit iunit2:
        nset2=0
        nstrt2=nstart+nset1
        do 1040 ibas=1,ndim
          if(nbouni(ibas).ne.iunit2) go to 1040
          if(nbotyp(ibas).lt.10) go to 1040
          nset2=nset2+1
          idel(nstrt2+nset2)=ibas
 1040     continue
        idel(nstart)=nset2
        ntot=nset1+nset2
c  now, zero all intersecting elements between the two sets:
        do 1060 i=1,nset1
          id=idel(nstart+i)
          do 1060 j=1,nset2
            jd=idel(nstrt2+j)
            if(id.eq.jd) go to 1060
            trf(id,jd)=zero
            trf(jd,id)=zero
 1060       continue
        nstart=nstart+ntot
 1100   continue
      go to 4000
c
c   zeroing of delocalization between groups of atoms
c
c   use the nbo type (nbotyp) and nbo label (label) lists.
 1200 continue
      mstart=0
      nstart=0
c  skip the 'blocks' before nset1:
      call hfld(keywd,leng,done)
      do 1400 k=1,ndel
c  read in the number of atoms in each set of the pair, nset1 and nset2:
c    (skip the 'by' between nset1 and nset2)
        call ifld(mset1,error)
        if(error) go to 9500
        call hfld(keywd,leng,done)
        call ifld(mset2,error)
        if(error) go to 9500
        mstart=mstart+2
        iscr1(mstart-1)=mset1
        iscr1(mstart)=mset2
c  read in the atoms of both sets:
        mtot=mset1+mset2
        do 1220 i=1,mtot
          call ifld(iscr1(mstart+i),error)
          if(error) go to 9500
 1220     continue
        mstrt2=mstart+mset1
        write(lfnpr,8350)
        write(lfnpr,8631) (iscr1(mstart+i),i=1,mset1)
        write(lfnpr,8360)
        write(lfnpr,8631) (iscr1(mstrt2+i),i=1,mset2)
        write(lfnpr,8370)
c  construct the list of the two sets of orbitals from the atom lists,
c    placing the orbital list in idel in the standard manner for itype=3:
        nstart=nstart+2
        nset1=0
        nset2=0
        do 1300 jbas=1,ndim
          donor=.false.
          accptr=.false.
          if(nbotyp(jbas).lt.20) donor=.true.
          if(nbotyp(jbas).ge.10) accptr=.true.
          list1=.false.
          list2=.false.
c    remember to consult ibxm before getting info from label!
          jb=ibxm(jbas)
          do 1240 j=4,6
            jat=label(jb,j)
            if(jat.eq.0) go to 1240
            do 1230 i=1,mset1
              iat=iscr1(mstart+i)
              if(iat.ne.jat) go to 1230
              go to 1240
 1230         continue
            go to 1250
 1240       continue
          list1=.true.
 1250     continue
          do 1270 j=4,6
            jat=label(jb,j)
            if(jat.eq.0) go to 1270
            do 1260 i=1,mset2
              iat=iscr1(mstrt2+i)
              if(iat.ne.jat) go to 1260
              go to 1270
 1260         continue
            go to 1280
 1270       continue
          list2=.true.
 1280     continue
          if(list1.and.list2) go to 1300
          if(.not.list1.and..not.list2) go to 1300
          if(list1.and..not.donor) go to 1300
          if(list2.and..not.accptr) go to 1300
          if(list2) go to 1290
c   list1.and.donor=.true. case:
            nset1=nset1+1
            idel(nstart+nset1)=jbas
            go to 1300
c   list2.and.accptr=.true. case:
 1290     continue
            nset2=nset2+1
            iscr2(nset2)=jbas
 1300   continue
c
        idel(nstart-1)=nset1
        idel(nstart)=nset2
        ntot=nset1+nset2
c  place orbital set 2 in idel:
        nstrt2=nstart+nset1
        do 1320 i=1,nset2
 1320     idel(nstrt2+i)=iscr2(i)
c  now, zero all intersecting elements between the two sets of orbitals:
        do 1340 i=1,nset1
          id=idel(nstart+i)
          do 1340 j=1,nset2
            jd=idel(nstrt2+j)
            trf(id,jd)=zero
 1340       trf(jd,id)=zero
        mstart=mstart+ntot
        nstart=nstart+ntot
 1400   continue
      go to 4000
c
c  delete all vicinal or geminal delocalizations:
c
 3000 ivic=1
      write(lfnpr,8550)
      goto 3020
 3010 ivic=0
      write(lfnpr,8560)
 3020 continue
      itype=3
c
c  start by filling trf with full nbo fock matrix:
c
      ii=0
      do 3025 i=1,ndim
        do 3025 j=1,i
          ii=ii+1
          trf(i,j)=f(ii)
          trf(j,i)=f(ii)
 3025 continue
c
c  find the total number of blocks of the fock matrix to delete:
c
      ndel=0
      nstart=0
      do 3070 ibas=1,ndim
        ib=ibxm(ibas)
        if(label(ib,2).ne.istar) then
          nacc=0
          do 3060 jbas=1,ndim
            jb=ibxm(jbas)
            if(label(jb,2).eq.istar) then
              itmp = ihtyp(ibas,jbas)
c
c  vicinal delocalization:
c
              if(ivic.eq.1.and.itmp.eq.lv) then
                nacc=nacc+1
                idel(nstart+nacc+3)=jbas
c
c  geminal delocalization:
c
              else if(ivic.eq.0.and.itmp.eq.lg) then
                nacc=nacc+1
                idel(nstart+nacc+3)=jbas
              end if
            end if
 3060     continue
          if(nacc.gt.0) then
            ndel=ndel+1
            idel(nstart+1)=1
            idel(nstart+2)=nacc
            idel(nstart+3)=ibas
            do 3065 jb=1,nacc
              jbas=idel(nstart+jb+3)
              if(jbas.ne.ibas) then
                trf(ibas,jbas)=zero
                trf(jbas,ibas)=zero
              end if
 3065       continue
            nstart=nstart+nacc+3
            if(nstart.gt.len) stop 'increase dimension of array idel'
          end if
        end if
 3070 continue
      goto 4000
c
c  write out information from deletion, for itype=3:
 4000 continue
      indx=0
      do 4050 k=1,ndel
        nset1=idel(indx+1)
        nset2=idel(indx+2)
        indx=indx+2
        nl=indx+1
        nu=indx+nset1
        write(lfnpr,8630)
        write(lfnpr,8631) (idel(i),i=nl,nu)
        write(lfnpr,8632)
        nl=indx+nset1+1
        nu=indx+nset1+nset2
        write(lfnpr,8631) (idel(i),i=nl,nu)
        indx=nu
 4050   continue
      return
c  delete all the "star" nbos on one or more molecules:
c   (set itype=1 for orbital deletions)
 5000 continue
      itype=1
c  read in the number of molecular units to "destar":
      call ifld(nunits,error)
      if(error) go to 9500
c  skip the keyword "units":
      leng=3
      call hfld(keywd,leng,done)
c  read in the numbers of the units to destar, finding the star orbitals
c   from the lists nbouni and nbotyp:
      ndel=0
      do 5100 i=1,nunits
        call ifld(iunit,error)
        if(error) go to 9500
        write(lfnpr,8400) iunit
        do 5050 ibas=1,ndim
          if(nbouni(ibas).ne.iunit) go to 5050
          if(label(ibas,2).ne.istar) go to 5050
          ndel=ndel+1
          idel(ndel)=ibas
 5050     continue
 5100   continue
c  go and do the deletions of the ndel orbitals that are now in idel:
      go to 100
c
c  delete all star nbos:
 5500 continue
      itype=1
      ndel=0
      write(lfnpr,8500)
      do 5600 ibas=1,ndim
        if(label(ibas,2).ne.istar) go to 5600
        ndel=ndel+1
        idel(ndel)=ibas
 5600   continue
c  go and do the deletions of the ndel orbitals that are now in idel:
      go to 100
c
 8100 format(1x,' ----------- alpha spin nbo deletions ----------- '/)
 8200 format(1x,' ----------- beta  spin nbo deletions ----------- '/)
 8300 format(1x,'zero delocalization from unit ',i2,' to unit ',i2)
 8350 format(1x,'zero delocalization from nbos localized on atoms:')
 8360 format(1x,'to nbos localized on atoms:')
 8370 format(1x,'    (nbos in common to the two groups of atoms ',
     *  'left out)')
 8400 format(1x,'destar unit ',i2,': delete all rydberg/antibond',
     * ' nbos from this unit')
 8500 format(1x,'nostar: delete all rydberg/antibond nbos')
 8550 format(1x,'novic: delete all vicinal delocalizations')
 8560 format(1x,'nogem: delete all geminal delocalizations')
 8610 format(1x,'deletion of the following orbitals ',
     * 'from the nbo fock matrix:',(/1x,20i4))
 8620 format(1x,'deletion of the following nbo fock matrix ',
     * 'elements:',/(7(2x,'(',i3,',',i3,')')))
 8630 format(1x,'deletion of the nbo fock matrix elements ',
     * 'between orbitals:')
 8631 format(1x,20i4)
 8632 format(1x,'and orbitals:')
 8640 format(1x,'deletion of the following nbo fock matrix ',
     * 'blocks:',/(2(2x,'(',i3,'-',i3,'/',i3,'-',i3,')')))
 8700 format(/)
c
c  error messages:
 9000 write(lfnpr,9010) (keywd(i),i=1,3)
 9010 format(1x,'first character string does not have the',
     * ' first three letters of delete or zero:',/1x,3a1)
      stop
 9100 write(lfnpr,9110)
 9110 format(1x,'non-integer was input for number of items to delete.')
      stop
 9200 write(lfnpr,9210) (keywd(i),i=1,3)
 9210 format(1x,'no match with first three letters of the keywords ',
     * 'for deletion type'/' (orbital,element,block) found:',
     * 3a1)
      stop
 9300 write(lfnpr,9310)
 9310 format(1x,'keyword alpha (or a) not found to start alpha nbo',
     *          ' deletion input.')
      stop
 9400 write(lfnpr,9410)
 9410 format(1x,'keyword beta (or b) not found to start beta nbo',
     *          ' deletion input.')
 9500 write(lfnpr,9510)
 9510 format(' there is an error in the input of deletions.')
      stop
      end
c*****************************************************************************
      subroutine newdm(dm,u,eig,ndim,idel,len,ndel,itype,nmoocc,ispin)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      parameter (maxatm = 750, maxbas = 4096)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       iatno(maxbas),ibxm(maxbas),nrank(2*maxbas),locc(2*maxbas)
      dimension dm(1),u(ndim,ndim),eig(ndim),idel(len)
      data zero/0.0d0/,one/1.0d0/,two/2.0d0/
c  onetwo: one if open shell (ispin.ne.0), two if closed shell (doubly occ mos)
      onetwo=two
      if(ispin.ne.0) onetwo=one
c  ntrunc: dimension of truncated fock matrix
      ntrunc=ndim
      if(itype.eq.1) ntrunc=ndim-ndel
c  rank the eigenvalues 'eig' from the truncated fock matrix from lowest
c   to highest in 'nrank':
      call rnkeig(nrank,eig,ntrunc,ndim,locc)
c  put in 'locc' the locations of the 'nmoocc' lowest eigenvalues:
c   (these correspond to the doubly occupied mos)
      nocc=0
      do 20 i=1,ntrunc
        if(nrank(i).gt.nmoocc) go to 20
          nocc=nocc+1
          locc(nocc)=i
   20   continue
c  ndelor: number of deleted orbitals
      ndelor=ndim-ntrunc
c
c  construct the new nbo density matrix:
c
c  loop over rows:
      ii=0
      ij=0
      iout=1
      do 105 i=1,ndim
        if(iout.gt.ndelor) go to 40
        if(i.ne.idel(iout)) go to 40
c  zero rows of the new nbo density matrix that were zeroed
c    in the truncation, also zeroing the orbital occpancy, eig(i):
          iout=iout+1
          eig(i)=zero
          do 30 j=1,i
            ij=ij+1
   30       dm(ij)=zero
          go to 105
   40   continue
        ii=ii+1
c  loop over columns:
        jout=1
        jj=0
        do 100 j=1,i
          if(jout.gt.ndelor) go to 50
          if(j.ne.idel(jout)) go to 50
c  zero columns of the new nbo density matrix that were zeroed
c    in the truncation of the nbo fock matrix:
            jout=jout+1
            ij=ij+1
            dm(ij)=zero
            go to 100
   50     continue
c  find dm(ij) from the eigenvectors of the truncated nbo fock matrix in 'u',
c  summing over the occupied mos, and multiplying by two for double occupancy:
          jj=jj+1
          sum=zero
          do 80 k=1,nmoocc
   80       sum=sum+u(ii,locc(k))*u(jj,locc(k))
          ij=ij+1
          dm(ij)=sum*onetwo
          if(i.eq.j) eig(i)=sum*onetwo
  100   continue
  105   continue
      return
      end
c*****************************************************************************
      subroutine rnkeig(rank,eig,n,ndim,arcrnk)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  rank eigenvalues in 'eig', lowest values first, in 'rank':
c
      integer rank,arcrnk
      dimension rank(ndim),eig(ndim),arcrnk(ndim)
      do 10 i=1,n
   10   arcrnk(i)=i
      do 40 i=1,n
        if(i.eq.n) go to 30
         i1=i+1
         do 20 j=i1,n
         if(eig(j).ge.eig(i)) go to 20
           temp=eig(i)
           eig(i)=eig(j)
           eig(j)=temp
           itemp=arcrnk(i)
           arcrnk(i)=arcrnk(j)
           arcrnk(j)=itemp
   20     continue
   30    rank(arcrnk(i))=i
   40   continue
      return
      end
c*****************************************************************************
      subroutine simltr(n,ndim,f,u,r,s,kntrol)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension f(1),u(ndim,1),s(1),r(1)
c  take u(transpose)*f*u:
c     f    matrix to be transformed (packed upper triangular)
c     u    is the transformation matrix
c     r    is the matrix in which the result will be returned
c     s    is a scratch matrix of dimension n
c     kntrol....=0  result returned only in  r
c               =1  result copied into  f
c
      in=0
      do 50 i=1,n
        jn=0
        do 20 j=1,n
          sum=0.0d0
          kn=0
          do 10 k=1,n
            jk=jn+k
            if(j.lt.k) jk=kn+j
            sum=sum+f(jk)*u(k,i)
   10       kn=kn+k
          s(j)=sum
   20     jn=jn+j
        do 40 j=1,i
          sum=0.0d0
          do 30 k=1,n
   30       sum=sum+s(k)*u(k,j)
          ij=in+j
   40     r(ij)=sum
   50   in=in+i
      if(kntrol.eq.0) return
      nt=n*(n+1)/2
      do 60 i=1,nt
   60   f(i)=r(i)
      return
      end
c*****************************************************************************
c
c  nbo direct access file (daf) routines:
c
c      subroutine nbfile(new,error)
c      subroutine nbopen(new,error)
c      subroutine nbwrit(ix,nx,idar)
c      subroutine nbread(ix,nx,idar)
c      subroutine nbclos(seq)
c      subroutine nbinqr(idar)
c
c      subroutine fetitl(title)
c      subroutine fee0(edel,etot)
c      subroutine sve0(edel)
c      subroutine fecoor(atcoor)
c      subroutine fesraw(s)
c      subroutine fedraw(dm,scr)
c      subroutine fefao(f,iwfock)
c      subroutine feaomo(t,it)
c      subroutine fedxyz(dxyz,i)
c      subroutine svnbo(t,occ,iscr)
c      subroutine fenbo(t,occ,iscr,nelec)
c      subroutine fetnbo(t)
c      subroutine svpnao(t)
c      subroutine fepnao(t)
c      subroutine svsnao(s)
c      subroutine fesnao(s)
c      subroutine svtnab(t)
c      subroutine fetnab(t)
c      subroutine svtlmo(t)
c      subroutine fetlmo(t)
c      subroutine svtnho(t)
c      subroutine fetnho(t)
c      subroutine svppao(dm)
c      subroutine feppao(dm)
c      subroutine svtnao(t)
c      subroutine fetnao(t)
c      subroutine svnlmo(t)
c      subroutine fenlmo(t)
c      subroutine svdnao(dm)
c      subroutine fednao(dm)
c      subroutine svfnbo(f)
c      subroutine fefnbo(f)
c      subroutine svnewd(dm)
c      subroutine fenewd(dm)
c      subroutine feinfo(icore,iswean)
c      subroutine febas(nshell,nexp,iscr)
c
c*****************************************************************************
      subroutine nbfile(new,error)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical new,error,need,there
      character*80 temp
c
      parameter (maxfil = 40)
c
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbname/filenm,nfile,ifile(maxfil)
      character*80 filenm
c
      data iwrit,iread/4hwrit,4hread/
c
c  create a list ifile of external lfns.  first find the files that
c  will be written:
c
      error = .false.
      nfile = 0
      do 10 i = 1,999
        need = .false.
        if(iwpnao.eq.-i)     need = .true.
        if(iwtnao.eq.-i)     need = .true.
        if(iwtnab.eq.-i)     need = .true.
        if(iwtnbo.eq.-i)     need = .true.
        if(jprint(7).eq. i)  need = .true.
        if(jprint(9).eq.-i)  need = .true.
        if(jprint(13).eq.-i) need = .true.
        if(jprint(15).eq.-i) need = .true.
        if(jprint(16).eq.-i) need = .true.
        if(jprint(17).eq.-i) need = .true.
        if(jprint(18).eq.-i) need = .true.
        if(jprint(19).eq.-i) need = .true.
        if(jprint(20).eq.-i) need = .true.
        if(jprint(21).eq.-i) need = .true.
        if(jprint(22).eq. i) need = .true.
        if(jprint(23).eq.-i) need = .true.
        if(jprint(24).eq.-i) need = .true.
        if(jprint(25).eq.-i) need = .true.
        if(jprint(26).eq.-i) need = .true.
        if(jprint(27).eq.-i) need = .true.
        if(jprint(28).eq.-i) need = .true.
        if(jprint(29).eq.-i) need = .true.
        if(jprint(30).eq.-i) need = .true.
        if(jprint(31).eq.-i) need = .true.
        if(jprint(33).eq.-i) need = .true.
        if(jprint(34).eq.-i) need = .true.
        if(jprint(35).eq.-i) need = .true.
        if(jprint(37).eq.-i) need = .true.
        if(jprint(38).eq.-i) need = .true.
        if(jprint(39).eq.-i) need = .true.
        if(jprint(40).eq.-i) need = .true.
        if(jprint(41).eq.-i) need = .true.
        if(jprint(42).eq.-i) need = .true.
        if(jprint(44).eq.-i) need = .true.
        if(jprint(45).eq.-i) need = .true.
        if(jprint(47).eq.-i) need = .true.
        if(jprint(48).eq.-i) need = .true.
        if(jprint(49).eq.-i) need = .true.
        if(jprint(50).eq.-i) need = .true.
        if(jprint(51).eq.-i) need = .true.
        if(jprint(52).eq.-i) need = .true.
        if(jprint(53).eq.-i) need = .true.
        if(jprint(54).eq.-i) need = .true.
        if(need) then
          nfile = nfile + 1
          if(nfile.gt.maxfil) then
            write(lfnpr,890) maxfil
            error = .true.
            return
          end if
          ifile(nfile) = i
        end if
   10 continue
c
c  add files that may be read:
c
      mfile = nfile
      if(ioinqr(iwpnao).eq.iread) then
        mfile = mfile + 1
        if(mfile.gt.maxfil) then
          write(lfnpr,890) maxfil
          error = .true.
          return
        end if
        ifile(mfile) = iwpnao/1000
      end if
      if(ioinqr(iwtnao).eq.iread) then
        mfile = mfile + 1
        if(mfile.gt.maxfil) then
          write(lfnpr,890) maxfil
          error = .true.
          return
        end if
        ifile(mfile) = iwtnao/1000
      end if
      if(ioinqr(iwtnab).eq.iread) then
        mfile = mfile + 1
        if(mfile.gt.maxfil) then
          write(lfnpr,890) maxfil
          error = .true.
          return
        end if
        ifile(mfile) = iwtnab/1000
      end if
c
c  make sure that no files are both written and read:
c
      do 30 i = nfile+1,mfile
        do 20 j = 1,nfile
          if(abs(ifile(i)).eq.ifile(j)) then
            write(lfnpr,900) ifile(j)
            error = .true.
            return
          end if
   20   continue
   30 continue
      nfile = mfile
c
c  also check that the nbo daf has its own lfn:
c
      do 40 i = 1,nfile
        if(abs(ifile(i)).eq.abs(lfndaf)) then
          write(lfnpr,900) ifile(i)
          error = .true.
          return
        end if
   40 continue
c
c  select an alternate filename if this one is not acceptable:
c
      temp = filenm
      do 50 i = 1,80
        if(temp(i:i).eq.char(32)) then
          length = i - 1
          go to 60
        end if
   50 continue
      length = 76
   60 continue
      io = ioinqr(iwpnao)
      jo = ioinqr(iwtnao)
      ko = ioinqr(iwtnab)
      if(new.and.io.ne.iread.and.jo.ne.iread.and.ko.ne.iread) then
        do 100 i = 0,999
          len = length
          if(i.ne.0) then
            ii = i
   65       len = len + 1
            temp(len:len) = char(mod(ii,10) + 48)
            ii = ii / 10
            if(ii.ne.0) goto 65
            if(len.eq.length+2) then
              temp(len+1:len+1) = temp(len:len)
              temp(len:len) = temp(len-1:len-1)
              temp(len-1:len-1) = temp(len+1:len+1)
            else if(len.eq.length+3) then
              temp(len+1:len+1) = temp(len:len)
              temp(len:len) = temp(len-2:len-2)
              temp(len-2:len-2) = temp(len+1:len+1)
            end if
          end if
          temp(len+1:len+1) = '.'
c
c  first check the daf:
c
          k = abs(lfndaf)
          if(abs(lfndaf).lt.100) k = k * 10
          temp(len+2:len+2) = char(k/100 + 48)
          temp(len+3:len+3) = char(mod(k/10,10) + 48)
          if(abs(lfndaf).lt.100) then
            temp(len+4:len+4) = char(32)
          else
            temp(len+4:len+4) = char(mod(k,10) + 48)
          end if
          inquire(file=temp,exist=there)
          if(there) go to 100
c
c  now check the rest:
c
          do 70 j = 1,nfile
            k = abs(ifile(j))
            if(abs(ifile(j)).lt.100) k = k * 10
            temp(len+2:len+2) = char(k/100 + 48)
            temp(len+3:len+3) = char(mod(k/10,10) + 48)
            if(abs(ifile(j)).lt.100) then
              temp(len+4:len+4) = char(32)
            else
              temp(len+4:len+4) = char(mod(k,10) + 48)
            end if
            inquire(file=temp, exist=there)
            if(there) go to 100
   70     continue
          go to 200
  100   continue
        write(lfnpr,910)
        error = .true.
        return
c
c  this is a good one!!  if the filename has changed, write a warning:
c
  200   continue
        if(filenm(1:len).ne.temp(1:len)) then
          filenm(1:len) = temp(1:len)
          do 210 i = len+1,80
            filenm(i:i) = char(32)
  210     continue
          write(lfnpr,920) filenm(1:52)
        end if
        length = len
      end if
c
c  open external files:
c
      temp = filenm
      temp(length+1:length+1) = '.'
      do 300 i = 1,nfile
        k = abs(ifile(i))
        if(abs(ifile(i)).lt.100) k = k * 10
        temp(length+2:length+2) = char(k/100 + 48)
        temp(length+3:length+3) = char(mod(k/10,10) + 48)
        if(abs(ifile(i)).lt.100) then
          temp(length+4:length+4) = char(32)
        else
          temp(length+4:length+4) = char(mod(k,10) + 48)
        end if
        if(ifile(i).gt.0) then
          open(unit=ifile(i), file=temp, status='new')
        else
          open(unit=abs(ifile(i)), file=temp, status='old')
        end if
  300 continue
      return
c
  890 format(/1x,'i/o is limited to ',i2,' files.  program abort.')
  900 format(/1x,'illegal request for input and output with lfn',i3)
  910 format(/1x,'the search for an acceptable filename has failed.')
  920 format(/1x,'filename:  changed to ',a52)
      end
c*****************************************************************************
      subroutine nbopen(new,error)
c*****************************************************************************
c
c  the following records of the nbo direct access file (daf) are used:
c
c          1  ---   nbodaf common block
c          2  ---   job title
c          3  ---   natoms,ndim,nbas,munit,wavefunction flags,iswean
c          4  ---   iatno,iznuc,lctr,lang
c          5  ---   ao basis set information
c          8  ---   deletion energy, total energy
c          9  ---   atomic coordinates
c         10  ---   ao overlap matrix
c         11  ---   pnao overlap matrix
c         20  ---   ao density matrix (alpha)
c         21  ---   ao density matrix (beta)
c         22  ---   pure ao density matrix
c         23  ---   nao density matrix (alpha)
c         24  ---   nao density matrix (beta)
c         25  ---   ao density matrix with nbo deletions (alpha)
c         26  ---   ao density matrix with nbo deletions (beta)
c         27  ---   nbo occupancies (alpha)
c         28  ---   nbo occupancies (beta)
c         30  ---   ao fock matrix (alpha)
c         31  ---   ao fock matrix (beta)
c         32  ---   nao fock matrix (alpha)
c         33  ---   nao fock matrix (beta)
c         34  ---   nbo fock matrix (alpha)
c         35  ---   nbo fock matrix (beta)
c         40  ---   ao to mo transformation matrix (alpha)
c         41  ---   ao to mo transformation matrix (beta)
c         42  ---   ao to pnao transformation matrix
c         43  ---   ao to nao transformation matrix
c         44  ---   ao to nbo transformation matrix  (alpha)
c         45  ---   ao to nbo transformation matrix  (beta)
c         46  ---   ao to nlmo transformation matrix
c         47  ---   nao to nho transformation matrix
c         48  ---   nao to nbo transformation matrix
c         49  ---   nbo to nlmo transformation matrix
c         50  ---   x dipole integrals
c         51  ---   y dipole integrals
c         52  ---   z dipole integrals
c         60  ---   nbo labels (alpha)
c         61  ---   nbo labels (beta)
c-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical new,error
      character*80 temp
c
c  note that isingl is no longer a parameter (6/7/90):
c
      parameter (length = 256)
      parameter (nbdar = 100)
      parameter (maxfil = 40)
c
      common/nbodaf/inbo,nav,ionbo(nbdar)
      common/nbonav/ixdnbo(length),nbnav,isingl
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbname/filenm,nfile,ifile(maxfil)
      character*80 filenm
c
      dimension ix(nbdar+2),ixsnbo(length/2)
c
      equivalence (ixsnbo(1),ixdnbo(1))
      equivalence (ix(1),inbo)
c
      save isw,lenrec
c
      data iblnk/1h /
      data isw/0/
c
c     inbo   :  fortran file number
c     ionbo  :  indexing array mapping the logical records of the
c               nbo daf onto the physical records of the disk file
c     nav    :  number of physical records currently on the daf
c     nbdar  :  maximum number of logical records on the daf
c
      inbo = abs(lfndaf)
c
c  are we working on a 32 (isingl=2) or 64 (isingl=1) bit machine?
c
      if(isw.eq.0) then
        do 10 i = 1,4
          iblnk = iblnk / 256
   10   continue
        if(iblnk.eq.0) then
          isingl = 2
        else
          isingl = 1
        end if
c
c  determine an appropriate record length for the nbo daf:
c
        lrec   = length / 4
        lenrec = 0
        do 30 i = 1,6
          lrec = lrec * 2
          open(unit=inbo, file='nb$temp.dat', status='new',
     +         access='direct', recl=lrec, form='unformatted',
     +         err=40)
          write(inbo,rec=1,err=20) ixdnbo
c
c  if i.eq.1 at this point, err did not work properly in the preceding
c  statement (this appears to be the case for the xl fortran compiler
c  running on an ibm risc station/6000):
c
          if(i.eq.1) lrec = length * 8 / isingl
          if(isingl.eq.1) lenrec = lrec / 2
          if(isingl.eq.2) lenrec = lrec
   20     close(unit=inbo, status='delete')
          if(lenrec.ne.0) go to 50
   30   continue
c
c  problems...
c
   40   continue
        write(lfnpr,900)
        error = .true.
        return
c
   50   continue
        isw = 1
      end if
c
c  open the nbo direct access file (daf) -- typically assigned to lfn48:
c
      temp = filenm
      do 60 i = 1,80
        if(temp(i:i).eq.char(32)) then
          len = i - 1
          go to 70
        end if
   60 continue
      len = 76
   70 continue
      k = inbo
      if(inbo.lt.100) k = k * 10
      temp(len+1:len+1) = '.'
      temp(len+2:len+2) = char(k/100 + 48)
      temp(len+3:len+3) = char(mod(k/10,10) + 48)
      if(inbo.lt.100) then
        temp(len+4:len+4) = char(32)
      else
        temp(len+4:len+4) = char(mod(k,10) + 48)
      end if
c
c  if this is a new nbo daf, write common/nbodaf/ on the first record:
c
      if(new) then
        open(unit=inbo, file=temp, status='new', access='direct',
     +       recl=lenrec, form='unformatted', err=110)
        nav   = 1
        nbnav = 1
        do 80 i = 1,nbdar
          ionbo(i) = 0
   80   continue
        nf = 1
        nx = (nbdar + 2) / isingl
        call nbwrit(ix,nx,nf)
c
c  otherwise, open the old file and read in common/nbodaf/ from the
c  first record:
c
      else
        open(unit=inbo, file=temp, status='old', access='direct',
     +       recl=lenrec, form='unformatted', err=110)
        nbnav = 1
        maxix = length * isingl/2
        ldar  = nbdar + 2
        max = 0
   90   min = max + 1
        max = max + maxix
        if(max.gt.ldar) max = ldar
        if(isingl.eq.1) read(inbo,rec=nbnav) ixsnbo
        if(isingl.eq.2) read(inbo,rec=nbnav) ixdnbo
        do 100 i = min,max
          ix(i) = ixdnbo(i-min+1)
  100   continue
        nbnav = nbnav + 1
        if(max.lt.ldar) go to 90
        inbo = abs(lfndaf)
      end if
      error = .false.
      return
c
c  error encountered while opening this file:
c
  110 error = .true.
      return
c
  900 format(/1x,'routine nbopen could not determine an appropriate ',
     + 'record length.')
      end
c*****************************************************************************
      subroutine nbwrit(ix,nx,idar)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (length = 256)
      parameter (nbdar = 100)
c
      common/nbodaf/inbo,nav,ionbo(nbdar)
      common/nbonav/ixdnbo(length),nbnav,isingl
c
      dimension ix(1),ixsnbo(length/2)
c
      equivalence (ixsnbo(1),ixdnbo(1))
c
      maxix = length * isingl / 2
      ldar  = nx * isingl
      if(ionbo(idar).ne.0) go to 100
c
c  if this is the first write to the nbo daf:
c
      ionbo(idar) = nav
      nbnav = nav
c
      max = 0
   10 min = max + 1
      max = max + maxix
      if(max.gt.ldar) max = ldar
      do 20 i = min,max
   20 ixdnbo(i-min+1) = ix(i)
      if(isingl.eq.1) write(inbo,rec=nbnav) ixsnbo
      if(isingl.eq.2) write(inbo,rec=nbnav) ixdnbo
      nbnav = nbnav + 1
      if(max.lt.ldar) go to 10
      nav = nbnav
      return
c
c  or if this is a rewrite:
c
  100 continue
      nbnav = ionbo(idar)
      max = 0
  110 min = max + 1
      max = max + maxix
      if(max.gt.ldar) max = ldar
      do 120 i = min,max
  120 ixdnbo(i-min+1) = ix(i)
      if(isingl.eq.1) write(inbo,rec=nbnav) ixsnbo
      if(isingl.eq.2) write(inbo,rec=nbnav) ixdnbo
      nbnav = nbnav + 1
      if(max.lt.ldar) go to 110
      return
      end
c*****************************************************************************
      subroutine nbread(ix,nx,idar)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (length = 256)
      parameter (nbdar = 100)
c
      common/nbodaf/inbo,nav,ionbo(nbdar)
      common/nbonav/ixdnbo(length),nbnav,isingl
c
      dimension ix(1),ixsnbo(length/2)
c
      equivalence (ixsnbo(1),ixdnbo(1))
c
      nbnav = ionbo(idar)
      maxix = length * isingl / 2
      ldar  = nx * isingl
c
      max = 0
   10 min = max + 1
      max = max + maxix
      if(max.gt.ldar) max = ldar
      if(isingl.eq.1) read(inbo,rec=nbnav) ixsnbo
      if(isingl.eq.2) read(inbo,rec=nbnav) ixdnbo
      do 20 i = min,max
   20 ix(i) = ixdnbo(i-min+1)
      nbnav = nbnav + 1
      if(max.lt.ldar) go to 10
      return
      end
c*****************************************************************************
      subroutine nbclos(seq)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical seq
c
      parameter (length = 256)
      parameter (nbdar = 100)
      parameter (maxfil = 40)
c
      common/nbodaf/inbo,nav,ionbo(nbdar)
      common/nbonav/ixdnbo(length),nbnav,isingl
      common/nbname/filenm,nfile,ifile(maxfil)
      character*80 filenm
c
      dimension ix(nbdar+2)
      equivalence (ix(1),inbo)
c
c  first close the nbo direct access file, remembering to write
c  common/nbodaf/ to the first logical record:
c
      nf = 1
      nx = (nbdar + 2) / isingl
      call nbwrit(ix,nx,nf)
      close(unit=inbo, status='keep')
c
c  then close the remainder of the files used by the nbo program:
c
      do 10 i = 1,nfile
        close(unit=abs(ifile(i)), status='keep')
   10 continue
      return
      end
c*****************************************************************************
      subroutine nbinqr(idar)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (nbdar = 100)
      common/nbodaf/inbo,nav,ionbo(nbdar)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      if(idar.lt.1.or.idar.gt.nbdar) then
        write(lfnpr,900) idar,nbdar
        stop
      end if
c
      if(ionbo(idar).eq.0) idar = 0
      return
c
  900 format(/1x,'nbo daf record out of range: idar = ',i4,
     + '  nbdar = ',i4)
      end
c*****************************************************************************
      subroutine fetitl(title)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension title(10)
c
c  fetitl:  fetches the job title from the nbodaf:
c
      nfile = 2
      call nbread(title,10,nfile)
      return
      end
c*****************************************************************************
      subroutine fee0(edel,etot)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension x(2)
c
c  fee0:  fetches the deletion and total scf energy
c
      nfile = 8
      call nbread(x,2,nfile)
      edel = x(1)
      etot = x(2)
      return
      end
c*****************************************************************************
      subroutine sve0(edel)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension x(2)
c
c  sve0:  saves the deletion energy
c
      nfile = 8
      call nbread(x,2,nfile)
      x(1) = edel
      call nbwrit(x,2,nfile)
      return
      end
c*****************************************************************************
      subroutine fecoor(atcoor)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension atcoor(3*natoms)
c
c  fecoor:  fetch the atomic cartesian coordinates in angstroms.
c
      nfile = 9
      call nbread(atcoor,3*natoms,nfile)
      return
      end
c*****************************************************************************
      subroutine febas(nshell,nexp,iscr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension iscr(1)
c
c  febas:  fetches the basis set info
c
      nfile = 5
      call nbinqr(nfile)
      if(nfile.gt.0) then
        call nbread(iscr,2,nfile)
        ii = 0
        ii = ii + 1
        nshell = iscr(ii)
        ii = ii + 1
        nexp   = iscr(ii)
        len    = 2 + 3*nshell + 5*nexp
        call nbread(iscr,len,nfile)
      else
        nshell = 0
      end if
      return
      end
c*****************************************************************************
      subroutine fesraw(s)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
      dimension s(ndim,ndim)
c
c  fesraw:  fetches the overlap matrix (raw ao. basis)
c           into s(ndim,ndim) a full square matrix.
c
      nfile = 10
      l2 = ndim*(ndim+1)/2
      call nbread(s,l2,nfile)
      call NBOunpack(s,ndim,nbas,l2)
      return
      end
c*****************************************************************************
      subroutine fedraw(dm,scr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension dm(1),scr(1)
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
      data nfilea,nfileb/20,21/
c
c  fedraw:  fetches the density matrix (raw a.o. basis) in dm(ndim,ndim)
c           if alpha =.true.  fetch alpha matrix
c           if beta  =.true   fetch beta matrix.
c           if open .and. .not.(alpha .or. beta) =.true  fetch the total d.m.
c
      l2 = ndim*(ndim+1)/2
      nfile = nfilea
      if(beta) nfile = nfileb
      call nbread(dm,l2,nfile)
c
      if(.not.open) goto 300
      if(alpha.or.beta) goto 300
      call nbread(scr,l2,nfileb)
c
c  form the total density matrix:
c
      do 100 i = 1,l2
        dm(i) = dm(i) + scr(i)
  100 continue
c
  300 call NBOunpack(dm,ndim,nbas,l2)
      return
      end
c*****************************************************************************
      subroutine fefao(f,iwfock)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension f(1)
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
      data nfilea,nfileb/30,31/
c
c  fefao:  fetches the ao fock matrix
c          if alpha .eq. .true.  we want the alpha fock matrix
c          if beta .eq. .true.  we want the beta fock matrix.
c          if the requested matrix does not exist then iwfock = 0
c
      l2 = ndim*(ndim+1)/2
      nfile = nfilea
      if(beta) nfile = nfileb
      call nbinqr(nfile)
      if(nfile.gt.0) then
        call nbread(f,l2,nfile)
        call NBOunpack(f,ndim,nbas,l2)
      else
        iwfock = 0
      end if
      return
      end
c*****************************************************************************
      subroutine feaomo(t,it)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension t(1)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
c
      data nfilea,nfileb/40,41/
c
c feaomo:  fetch the ao to mo transformation matrix:
c          (it = 1, ao to mo transform is on nbo daf)
c          (it = 0, ao to mo transform is not on nbo daf)
c
      nfile = nfilea
      if (beta) nfile = nfileb
      call nbinqr(nfile)
      if(nfile.gt.0) then
        it = 1
        l3 = ndim*ndim
        call nbread(t,l3,nfile)
      else
        it = 0
      end if
      return
      end
c*****************************************************************************
      subroutine fedxyz(dxyz,i)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension dxyz(1)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
      data nfilex,nfiley,nfilez/50,51,52/
c
c  fedxyz:    fetch the ao dipole moment matrices (in angstroms)
c      i=1:  x       i=2:    y           i=3:   z
c
      if(i.eq.1) nfile = nfilex
      if(i.eq.2) nfile = nfiley
      if(i.eq.3) nfile = nfilez
c
      call nbinqr(nfile)
      if(nfile.gt.0) then
        l2 = ndim*(ndim+1)/2
        call nbread(dxyz,l2,nfile)
        call NBOunpack(dxyz,ndim,nbas,l2)
      else
        i = 0
      end if
      return
      end
c*****************************************************************************
      subroutine svnbo(t,occ,iscr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbatom/iatno(maxatm),ino(maxatm),norb(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
c
      dimension t(ndim,ndim),occ(ndim),iscr(1)
c
c  svnbo:  saves nbo information (transformation, occupancies, labels, etc.)
c          if alpha .eq. .true.  save the alpha information
c          if beta .eq. .true.  save the beta information.
c
c  save the ao to nbo transformation matrix:
c
      l1 = ndim
      l3 = ndim*ndim
      l4 = 10*ndim
      nfile = 44
      if (beta) nfile = 45
      call nbwrit(t,l3,nfile)
c
c  save nbo orbital occupancies:
c
      nfile = 27
      if (beta) nfile = 28
      call nbwrit(occ,l1,nfile)
c
c  save the lists of nbo information for later use in the deletions.
c  pack the information into iscr(10*ndim):
c
      ii = 0
      do 40 k = 1,6
        do 30 i = 1,nbas
          ii = ii + 1
          iscr(ii) = label(i,k)
   30   continue
   40 continue
      do 50 i = 1,nbas
        ii = ii + 1
        iscr(ii) = ibxm(i)
   50 continue
      do 60 i = 1,natoms
        ii = ii + 1
        iscr(ii) = iatno(i)
   60 continue
      do 70 i = 1,nbas
        ii = ii + 1
        iscr(ii) = nbouni(i)
   70 continue
      do 80 i = 1,nbas
        ii = ii + 1
        iscr(ii) = nbotyp(i)
   80 continue
c
      nfile = 60
      if (beta) nfile = 61
      call nbwrit(iscr,l4,nfile)
c
      return
      end
c*****************************************************************************
      subroutine fenbo(t,occ,iscr,nelec)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       iatno(maxbas),ibxm(maxbas),iscr1(2*maxbas),iscr2(2*maxbas)
c
      dimension t(ndim,ndim),occ(ndim),iscr(1)
c
      data zero,tenth /0.0d0,1.0d-1/
c
c  fenbo:  fetches nbo information (transformation, occupancies, labels, etc.)
c          if alpha .eq. .true.  fetch the alpha information
c          if beta .eq. .true.  fetch the beta information.
c
c  fetch the ao to nbo transformation matrix:
c
      l1 = ndim
      l3 = ndim*ndim
      l4 = ndim*10
      nfile = 44
      if (beta) nfile = 45
      call nbread(t,l3,nfile)
c
c  fetch nbo orbital occupancies:
c
      nfile = 27
      if (beta) nfile = 28
      call nbread(occ,l1,nfile)
c
c  count up the total number of electrons as an integer nelec:
c
      ele = zero
      do 10 i = 1,nbas
        ele = ele + occ(i)
   10 continue
      ele = ele + tenth
      nelec = ele
c
c  fetch the various lists of nbo information for use in the deletions.
c  unpack the information into label(maxbas,6),ibxm(maxbas),iatno(maxbas),
c  nbouni(maxbas) and nbotyp(maxbas) from iscr(10*ndim):
c
      nfile = 60
      if (beta) nfile = 61
      call nbread(iscr,l4,nfile)
c
      ii = 0
      do 40 k = 1,6
        do 30 i = 1,nbas
          ii = ii + 1
          label(i,k) = iscr(ii)
   30   continue
   40 continue
      do 50 i = 1,nbas
        ii = ii + 1
        ibxm(i) = iscr(ii)
   50 continue
      do 60 i = 1,natoms
        ii = ii + 1
        iatno(i) = iscr(ii)
   60 continue
      do 70 i = 1,nbas
        ii = ii + 1
        nbouni(i) = iscr(ii)
   70 continue
      do 80 i = 1,nbas
        ii = ii + 1
        nbotyp(i) = iscr(ii)
   80 continue
c
      return
      end
c*****************************************************************************
      subroutine fetnbo(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension t(1)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
c  fetnbo: fetch the ao to nbo transformation matrix
c
      l3 = ndim*ndim
      nfile = 44
      if (beta) nfile = 45
      call nbread(t,l3,nfile)
      return
      end
c*****************************************************************************
      subroutine svpnao(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension t(1)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
c  svpnao:  saves the ao to pnao transformation matrix.
c
      nfile = 42
      l3 = ndim*ndim
      call nbwrit(t,l3,nfile)
      return
      end
c*****************************************************************************
      subroutine fepnao(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension t(1)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
c  fepnao:  fetches the ao to pnao transformation matrix.
c
      nfile = 42
      l3 = ndim*ndim
      call nbread(t,l3,nfile)
      return
      end
c*****************************************************************************
      subroutine svsnao(s)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension s(ndim,ndim)
c
c   svsnao:  save the overlap matrix in the pnao or rpnao basis set.
c
      nfile = 11
      l2 = ndim*(ndim+1)/2
      call NBOpack(s,ndim,nbas,l2)
      call nbwrit(s,l2,nfile)
      call NBOunpack(s,ndim,nbas,l2)
      return
      end
c*****************************************************************************
      subroutine fesnao(s)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
      dimension s(ndim,ndim)
c
c   fesnao:   fetch the overlap matrix in the pnao or rpnao basis set.
c
      nfile = 11
      l2 = ndim*(ndim+1)/2
      call nbread(s,l2,nfile)
      call NBOunpack(s,ndim,nbas,l2)
      return
      end
c*****************************************************************************
      subroutine svtnab(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension t(ndim,ndim)
c
c  svtnab:  save the nao to nbo transformation matrix.
c
      nfile = 48
      l3 = ndim*ndim
      call nbwrit(t,l3,nfile)
      return
      end
c*****************************************************************************
      subroutine fetnab(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension t(ndim,ndim)
c
c  fetnab:  fetch the nao to nbo transformation matrix
c
      nfile = 48
      l3 = ndim*ndim
      call nbread(t,l3,nfile)
      return
      end
c*****************************************************************************
      subroutine svtlmo(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension t(ndim,ndim)
c
c  svtlmo:  save the nbo to nlmo transformation matrix.
c
      nfile = 49
      l3 = ndim*ndim
      call nbwrit(t,l3,nfile)
      return
      end
c*****************************************************************************
      subroutine fetlmo(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension t(ndim,ndim)
c
c  fetlmo:  fetch the nbo to nlmo transformation matrix
c
      nfile = 49
      l3 = ndim*ndim
      call nbread(t,l3,nfile)
      return
      end
c*****************************************************************************
      subroutine svtnho(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension t(ndim,ndim)
c
c  svtnho:   temporarily save the nao to nho transformation
c
      nfile = 47
      l3 = ndim*ndim
      call nbwrit(t,l3,nfile)
      return
      end
c*****************************************************************************
      subroutine fetnho(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension t(ndim,ndim)
c
c  fetnho:   fetch the nao to nho transformation
c
      nfile = 47
      l3 = ndim*ndim
      call nbread(t,l3,nfile)
      return
      end
c*****************************************************************************
      subroutine svppao(dm)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension dm(ndim,ndim)
c
c  svppao:  temporarily saves the pure ao (pao) density matrix.
c           (this is not the raw ao basis, but the basis after the
c           transformation from cartesian to pure d,f,g functions).
c
      nfile = 22
      l2 = ndim*(ndim+1)/2
      call NBOpack(dm,ndim,nbas,l2)
      call nbwrit(dm,l2,nfile)
      call NBOunpack(dm,ndim,nbas,l2)
      return
      end
c*****************************************************************************
      subroutine feppao(dm)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension dm(ndim,ndim)
c
c  feppao:  fetches the pure ao (pao) density matrix.
c           (this is not the raw ao basis, but the basis after the
c           transformation from cartesian to pure d,f,g functions).
c
      nfile = 22
      l2 = ndim*(ndim+1)/2
      call nbread(dm,l2,nfile)
      call NBOunpack(dm,ndim,nbas,l2)
      return
      end
c*****************************************************************************
      subroutine svtnao(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension t(ndim,ndim)
c
c  svtnao:  save the ao to nao transformation matrix.
c
      if(.not.ortho) then
        nfile = 43
        l3 = ndim*ndim
        call nbwrit(t,l3,nfile)
      end if
      return
      end
c*****************************************************************************
      subroutine fetnao(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension t(ndim,ndim)
c
      data zero,one/0.0d0,1.0d0/
c
c  fetnao:  fetches the ao to nao transformation matrix.
c
      if(ortho) then
        do 20 j = 1,ndim
          do 10 i = 1,ndim
            t(i,j) = zero
   10     continue
          t(j,j) = one
   20   continue
      else
        nfile = 43
        l3 = ndim*ndim
        call nbread(t,l3,nfile)
      end if
      return
      end
c*****************************************************************************
      subroutine svnlmo(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension t(ndim,ndim)
c
c  svnlmo:  save the ao to nlmo transformation matrix
c
      nfile = 46
      l3 = ndim*ndim
      call nbwrit(t,l3,nfile)
      return
      end
c*****************************************************************************
      subroutine fenlmo(t)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension t(ndim,ndim)
c
c  fenlmo:  fetch the ao to nlmo transformation matrix
c
      nfile = 46
      l3 = ndim*ndim
      call nbread(t,l3,nfile)
      return
      end
c*****************************************************************************
      subroutine svdnao(dm)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension dm(ndim,ndim)
c
c  svdnao:  save the nao density matrix
c
      if(.not.ortho) then
        nfile = 23
        if(beta) nfile = 24
        l2 = ndim*(ndim+1)/2
        call NBOpack(dm,ndim,nbas,l2)
        call nbwrit(dm,l2,nfile)
        call NBOunpack(dm,ndim,nbas,l2)
      end if
      return
      end
c*****************************************************************************
      subroutine fednao(dm)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      dimension dm(ndim,ndim)
c
c  fednao:  fetches the nao density matrix (ao dm for orthogonal basis sets)
c
      if(ortho) then
        call fedraw(dm,dm)
      else
        nfile = 23
        if(beta) nfile = 24
        l2 = ndim*(ndim+1)/2
        call nbread(dm,l2,nfile)
        call NBOunpack(dm,ndim,nbas,l2)
      end if
      return
      end
c*****************************************************************************
      subroutine svfnbo(f)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
c
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      dimension f(ndim,ndim)
c
c  svfnbo:  saves the nbo fock matrix
c
      nfile = 34
      if (beta) nfile = 35
      l2 = ndim*(ndim+1)/2
      call NBOpack(f,ndim,nbas,l2)
      call nbwrit(f,l2,nfile)
      call NBOunpack(f,ndim,nbas,l2)
      return
      end
c*****************************************************************************
      subroutine fefnbo(f)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      dimension f(1)
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
c  fefnbo:  fetches the nbo fock matrix, leaving it in triangular form!!
c           if alpha.eq.true.  we want the alpha fock matrix
c           if beta.eq.true.   we want the beta fock matrix.
c
      nfile = 34
      if (beta) nfile = 35
      l2 = ndim*(ndim+1)/2
      call nbread(f,l2,nfile)
      return
      end
c*****************************************************************************
      subroutine svnewd(dm)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      dimension dm(1)
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
c svnewd:  save the new density matrix (raw ao basis) from nbo deletion
c
      nfile = 25
      if (beta) nfile = 26
      l2 = ndim*(ndim+1)/2
      call nbwrit(dm,l2,nfile)
      return
      end
c*****************************************************************************
      subroutine fenewd(dm)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      dimension dm(1)
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
c fenewd:  fetch the new density matrix (raw ao basis)
c
      nfile = 25
      if (beta) nfile = 26
      l2 = ndim*(ndim+1)/2
      call nbread(dm,l2,nfile)
      return
      end
c*****************************************************************************
      subroutine feinfo(icore,iswean)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      dimension icore(12)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbbas/label(maxbas,6),lval(maxbas),imval(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbao/lctr(maxbas),lang(maxbas)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
c  restore wavefunction information from the nbo daf:
c
c  restore natoms, ndim, nbas, munit, wavefunction flags, iswean:
c
      nfile = 3
      call nbread(icore,12,nfile)
      natoms = icore(1)
      ndim   = icore(2)
      nbas   = icore(3)
      munit  = icore(4)
      rohf   = .false.
      if(icore(5).eq.1)  rohf  = .true.
      uhf    = .false.
      if(icore(6).eq.1)  uhf   = .true.
      ci     = .false.
      if(icore(7).eq.1)  ci    = .true.
      open   = .false.
      if(icore(8).eq.1)  open  = .true.
      mcscf  = .false.
      if(icore(9).eq.1)  mcscf = .true.
      auhf   = .false.
      if(icore(10).eq.1) auhf  = .true.
      ortho  = .false.
      if(icore(11).eq.1) ortho = .true.
      iswean = icore(12)
c
c  if iswean is 1, set icore(12) to 0 and write to nbo daf.  note, iswean is
c  set to 1 by the feaoin driver routine.  this tells the energetic analysis
c  routines to search for the $del keylist.  iswean is reset to 0 here so
c  that multiple deletions can be read from a single $del keylist:
c
      if(iswean.eq.1) then
        icore(12) = 0
        call nbwrit(icore,12,nfile)
      end if
      return
      end
c*****************************************************************************
c
c  free format input routines:
c
c      subroutine strtin(lfnin)
c      subroutine rdcrd
c      subroutine ifld(int,error)
c      subroutine rfld(real,error)
c      subroutine hfld(keywd,leng,endd)
c      subroutine fndfld
c      function equal(ia,ib,l)
c
c*****************************************************************************
c
c  user  instructions:
c
c     1. the character string "end" is the field terminating mark:
c
c     2. commas and equal signs are treated as equivalent to blanks.
c          commas, equal signs, and blanks delimit input items.
c
c     3. all characters to the right of an exclamation mark ! are treated as
c          comments, and the next card is read in when these are encountered.
c
c     4. upper and lower case characters can be read by these routines.
c          however, lower case characters are converted to upper case
c          when encountered.
c
c     5. to read in data for the first time from lfn "lfnin" (perhaps
c          after using these subroutines to read in data from another lfn),
c          or to continue reading in data from lfnin after encountering
c          a field terminating mark, call strtin(lfnin)  (start input)
c
c     6. to fetch the next non-blank string of characters from lfn lfnin,
c           call hfld(keywd,length,end),
c            where keywd   is a vector of dimension "length"  or longer,
c                  length  is the maximum number of characters to fetch,
c                  end     must be a declared logical variable.
c           upon return,
c            end=.true. if a field terminating mark was found to be the next
c                 non-blank character string.  otherwise, end=.false.
c            end=.true. and length=0 means the end-of-file was found.
c            length is changed to the actual number of characters in string
c                 if this is less than the value of length set by the calling
c                 program.
c            keywd(1) through keywd(length) contain the character string,
c                 one character per element of keywd.
c
c     7. to fetch the integer value of the next character string,
c           call ifld(int,error),
c            where int     is the variable to be read,
c                  error   must be a declared logical variable.
c            upon return,
c             if error=.false., an integer was found and placed in "int".
c             if error=.true. and int.gt.0, a field terminating mark was
c                 found as the next character string.
c             if error=.true. and int.lt.0, the next character string found
c                 was neither an integer nor a field terminating mark.
c
c     8. to fetch the real value of the next character string,
c           (an exponent is allowed, with or without an "e" or "f".
c             if no letter is present to signify the exponent field,
c             a + or - sign must start the exponent.  if no mantissa is
c             present, the exponent field must start with a letter, and
c             the mantissa is set to one.)
c           call rfld(real,error),
c            where real    is the variable to be read,
c                  error   must be a declared logical variable.
c            upon return,
c             if error=.false., a real number was found and placed in "real".
c             if error=.true. and real.gt.1, a field terminating mark was
c                 found as the next character string.
c             if error=.true. and real.lt.-1, the next character string found
c                 was neither a real number nor a field terminating mark.
c
c     9. to compare the corresponding first l elements of each of two vectors
c          ia(l) and ib(l) to see if the vectors are equivalent,
c           use the function equal(ia,ib,l).
c           equal must be declared logical in the calling program,
c            and the function value (.true. or .false.) will tell if the
c            vectors ia and ib are equal up to element l.
c        note: this function is useful for determining if a character string
c          read by a call to hfld matches a certain keyword which is stored
c          in a vector, one character per element.
c
c
c*****************************************************************************
      subroutine strtin(lfnin)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbcrd1/icd(80),look(80),length,ipt,lfn,nexp
      common/nbcrd2/point,end,next,exp
      logical point,end,next,exp
c
c  initialize input from lfn lfnin:
c
      lfn  = lfnin
      end  = .false.
      next = .true.
      call rdcrd
c
      return
      end
c*****************************************************************************
      subroutine rdcrd
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  subroutine name changed from rdcard, due to conflict with gamess:
c
      common/nbcrd1/icd(80),look(80),length,ipt,lfn,nexp
      common/nbcrd2/point,end,next,exp
      logical point,end,next,exp
c
      data ia,ichara,icharz/1ha,1hA,1hZ/
      data iblnk,iq,ii/1h ,1h`,1hi/
c
c  read in the next card at lfn:
c
      read(lfn,1000,end=800,err=800) icd
c
c  change all lower case characters to upper case:
c
      call lcase2(icd)
c     do 10 i = 1,80
c       if(icd(i).ge.ichara.and.icd(i).le.icharz) then
c         icd(i) = icd(i) + ichara - ia
c       end if
c  10 continue
c
c  treat tabs as spaces:
c
      itab = iblnk + ii - iq
      do 20 i = 1,80
	if(icd(i).eq.itab) icd(i) = iblnk
   20 continue
c
c  reset column pointer, ipt:
c
      ipt = 1
      return
c
c  end of file encountered
c
  800 continue
      end = .true.
      return
c
 1000 format(80a1)
      end
c*****************************************************************************
      subroutine ifld(int,error)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical error
c
      common/nbcrd1/icd(80),look(80),length,ipt,lfn,nexp
      common/nbcrd2/point,end,next,exp
      logical point,end,next,exp
c
      data zero,one,small/0.0d0,1.0d0,1.0d-3/
c
c  search lfn for the next string of non-blank characters, see if they
c  form an integer (if not, error=.true.) and, if so, place its numerical
c  value in "int":
c
      int = 0
      call rfld(real,error)
c
c  if decimal point or an exponent.lt.0, error = .true.:
c
      if(exp) go to 100
      if(point) go to 100
      if(nexp.lt.0) go to 100
      if(length.eq.0) go to 100
      sign = one
      if(real.lt.zero) sign = -one
      real = real + small * sign
      int = real
      return
c
  100 error = .true.
      next = .false.
      return
      end
c*****************************************************************************
      subroutine rfld(real,error)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical error,expsgn,mantis
c
      common/nbcrd1/icd(80),look(80),length,ipt,lfn,nexp
      common/nbcrd2/point,end,next,exp
      logical point,end,next,exp
c
      dimension nchar(15)
c
      data nchar/1h0,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1h.,1h+,1h-,
     +  1hd,1he/
      data zero,one,ten/0.0d0,1.0d0,10.0d0/
c
c  search lfn for the next string of non-blank characters, see if they form
c  a real number (exponent is optional) (if not, error=.true.) and, if so,
c  place its numerical value in "real":
c
      real   = zero
      sign   = one
      ndec   = 0
      isexp  = 1
      nexp   = 0
      expsgn = .false.
      exp    = .false.
      point  = .false.
      error  = .false.
      mantis = .false.
      end    = .false.
c
c  find the next string of non-blank characters, "look", of length "length":
c
      if(next) call fndfld
      if(end) go to 300
      if(length.eq.0) go to 300
c
c  find the numerical value of the characters in "look":
c
      do 200 j = 1,length
        lk = look(j)
        do 20 i = 1,15
          if(lk.eq.nchar(i)) go to 40
   20   continue
        go to 300
   40   k = i - 11
        if(k) 60,80,100
c
c  this character is a number:
c
   60     continue
          if(exp) go to 70
c
c  add digit to mantissa:
c
          mantis = .true.
          real = real * ten + dble(i - 1)
c
c  if we are to the right of a decimal point, increment the decimal counter:
c
          if(point) ndec = ndec + 1
          go to 200
c
c  add digit to exponent:
c
   70     nexp = nexp * 10 + (i - 1)
          go to 200
c
c  decimal point:
c
   80     if(point) go to 300
          point = .true.
          go to 200
c
c  exponent (+,-,d,e):
c
  100     continue
          go to (110,130,150,150), k
c
c  plus sign: if not first character, count as part of exponent:
c
  110       if(j.eq.1) go to 200
              if(expsgn) go to 200
              expsgn = .true.
              exp = .true.
              go to 200
c
c  minus sign: if not first character, count as part of exponent:
c
  130       if(j.ne.1) go to 140
              sign = -one
              go to 200
  140         isexp = -1
              if(expsgn) go to 200
              expsgn = .true.
              exp = .true.
              go to 200
c
c  d or e: start of exponent:
c
  150       if(exp) go to 300
            exp = .true.
  200  continue
c
c  set final value of real (if no mantissa, but exponent present,
c  set mantissa to one):
c
      if(exp.and..not.mantis) real = one
      real = real * sign * (ten**(-ndec+isexp*nexp))
      next = .true.
      return
c
c  no real number found, or field terminating mark:
c
  300 continue
      error = .true.
      real  = -ten
      if(end) real = ten
      return
      end
c*****************************************************************************
      subroutine hfld(keywd,leng,endd)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical endd,equal
c
      common/nbcrd1/icd(80),look(80),length,ipt,lfn,nexp
      common/nbcrd2/point,end,next,exp
      logical point,end,next,exp
c
      dimension keywd(leng),kend(3)
c
      data nbla/1h /
      data kend/1he,1hn,1hd/
c
c  search lfn and find next non-blank string of characters and place
c  in the vector "keywd".  leng, from the calling program, is maximum
c  length of string to put in the vector keywd.  if "length" is less
c  than "leng", leng is set to length upon return:
c
      if(next) call fndfld
      endd  = end
      leng1 = leng
      leng  = min(length,leng)
c
c  place leng characters into keywd:
c
      do 10 i = 1,leng
        keywd(i) = look(i)
   10 continue
c
c  fill the rest of keywd with blanks:
c
      do 20 i = leng+1,leng1
        keywd(i) = nbla
   20 continue
      next = .true.
c
c  check for end of input:
c
      if(equal(look,kend,3)) endd = .true.
      return
      end
c*****************************************************************************
      subroutine fndfld
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbcrd1/icd(80),look(80),length,ipt,lfn,nexp
      common/nbcrd2/point,end,next,exp
      logical point,end,next,exp
c
      data nbla/1h /,ncom/1h,/,nexc/1h!/,neq/1h=/
c
c  find next non-blank string of characters in lfn.  read in another line
c  of lfn until non-blank characters are found and place them in "look",
c  of length "length":
c
      if(end) go to 35
      if(ipt.ge.80) call rdcrd
      if(end) go to 35
c
c  look for start of field.  skip to next card if "!" is encountered
c  (comment field):
c
   10 continue
      do 20 ncol = ipt,80
        icard = icd(ncol)
        if(icard.eq.nexc) go to 30
        if(icard.ne.nbla.and.icard.ne.ncom.and.icard.ne.neq) go to 40
   20 continue
c
c  nothing additional found on this card, continue with the next card:
c
   30 call rdcrd
      if(.not.end) go to 10
c
c  end of file found:
c
   35 length = 0
      return
c
c  look for the end of this field, counting characters as we go and
c  storing these character in look:
c
   40 m = 0
      do 80 mcol = ncol,80
        ichar = icd(mcol)
        if(ichar.eq.nbla.or.ichar.eq.ncom.or.ichar.eq.neq) go to 100
        m = m + 1
        look(m) = ichar
   80 continue
c
c  set length to the length of the new string in look and reset ipt to
c  the next space after this string:
c
  100 length = m
      ipt = mcol
      next = .false.
      return
      end
c*****************************************************************************
      function equal(ia,ib,l)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical equal
c
      dimension ia(l),ib(l)
c
c  test if the first l elements of vectors ia and ib are equal:
c
      equal = .false.
      do 10 i = 1,l
        if(ia(i).ne.ib(i)) go to 20
   10 continue
      equal = .true.
   20 return
      end
c rh
c test if ia (keyword from input) and ib (all keys) are equal
c note that this function returns true if 1: ia(1:l) .eq. ib(1:l),
c 2: ia(2:l) .eq. ib(2:l) or 3: ia(1:l-1).eq.ib(2:l)
c thus equal2 returns true if ia=%key and ib=$key ...
c
      function equal2(ia,ib,l)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical equal2
c
      dimension ia(l),ib(l)
c
c  test if the first l elements of vectors ia and ib are equal:
c
      equal2 = .false.
      do 10 i = 1,l
        if(ia(i).ne.ib(i)) go to 20
   10 continue
      equal2 = .true.
      goto 30
   20 continue
c
      do i=2,l
        if (ia(i).ne.ib(i)) go to 25
      enddo
      equal2 = .true.
      goto 30
   25 continue
      do i=1,l-1
        if (ia(i).ne.ib(i+1)) go to 30
      enddo
      equal2 = .true.
   30 return
      end
c*****************************************************************************
c
c  other system-independent i/o routines:
c
c      subroutine geninp(newdaf)
c      subroutine nboinp(nboopt,idone)
c      subroutine corinp(iess,icor)
c      subroutine chsinp(iess,ichs)
c      subroutine delinp(nboopt,idone)
c
c      subroutine rdcore(jcore)
c      subroutine wrppna(t,occ,iflg)
c      subroutine rdppna(t,occ,iflg)
c      subroutine wrtnao(t,iflg)
c      subroutine rdtnao(dm,t,scr,iflg)
c      subroutine wrtnab(t,iflg)
c      subroutine rdtnab(t,dm,bndocc,scr,iflg)
c      subroutine wrtnbo(t,bndocc,iflg)
c      subroutine wrnlmo(t,dm,iflg)
c      subroutine wrbas(scr,iscr,lfn)
c      subroutine wrarc(scr,iscr,lfn)
c
c      subroutine aout(a,mr,nr,nc,title,index,iflg)
c      subroutine aprint(a,mr,nr,nc,title,index,mcol)
c      subroutine awrite(a,mr,nr,nc,title,lfn)
c      subroutine aread(a,mr,nr,nc,job,lfn,error)
c      subroutine altout(a,mr,mc,nr,nc)
c      subroutine keypar(string,len,iflg,lfn,read,error)
c      function ioinqr(iflg)
c      subroutine lblao
c      subroutine lblnao
c      subroutine lblnbo
c      subroutine lblnho(inho,inbo,ictr,nctr)
c
c*****************************************************************************
      subroutine geninp(newdaf)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical newdaf,end,error,equal
c
      dimension keywd(6),kgen(4),kend(4),kreuse(5),knbas(4),knatom(6),
     +      kupper(5),kopen(4),kortho(5),kbohr(4),kbodm(4),kev(2),
     +      kcubf(6)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbgen/reuse,upper,bohr,denop
      logical reuse,upper,bohr,denop
c
      data kgen/1h$,1hg,1he,1hn/,kend/1h$,1he,1hn,1hd/,
     + kreuse/1hr,1he,1hu,1hs,1he/,knbas/1hn,1hb,1ha,1hs/,
     + knatom/1hn,1ha,1ht,1ho,1hm,1hs/,kupper/1hu,1hp,1hp,1he,1hr/,
     + kopen/1ho,1hp,1he,1hn/,kortho/1ho,1hr,1ht,1hh,1ho/,
     + kbohr/1hb,1ho,1hh,1hr/,kbodm/1hb,1ho,1hd,1hm/,
     + kev/1he,1hv/kcubf/1hc,1hu,1hb,1hi,1hc,1hf/
c
c  initialize variables:
c
      nbas   = 0
      natoms = 0
      munit  = 0
      reuse  = .false.
      upper  = .false.
      bohr   = .false.
      denop  = .true.
c
c  search lfnin for $gen:
c
      rewind(lfnin)
   10 call strtin(lfnin)
      len = 6
      call hfld(keywd,len,end)
      if(len.eq.0.and.end) stop 'no $gen keylist in the input file'
      if(.not.equal(keywd,kgen,4)) goto 10
c
c  $gen has been found, now read keywords:
c
   20 len = 6
      call hfld(keywd,len,end)
      if(equal(keywd,kend,4)) goto 700
c
c  keyword reuse -- reuse data already stored on the nbo daf:
c
      if(equal(keywd,kreuse,5)) then
        reuse = .true.
        goto 20
      end if
c
c  keyword nbas -- specify the number of basis functions:
c
      if(equal(keywd,knbas,4)) then
        call ifld(nbas,error)
        if(error) stop 'error reading in number of basis functions nbas'
        goto 20
      end if
c
c  keyword natoms -- specify the number of atoms:
c
      if(equal(keywd,knatom,4)) then
        call ifld(natoms,error)
        if(error) stop 'error reading in number of atoms natoms'
        goto 20
      end if
c
c  keyword upper -- read only upper triangular portions of matrices:
c
      if(equal(keywd,kupper,5)) then
        upper = .true.
        goto 20
      end if
c
c  keyword open -- open shell species (alpha and beta matrices read):
c
      if(equal(keywd,kopen,4)) then
        open = .true.
        goto 20
      end if
c
c  keyword ortho -- orthogonal basis set (skip nao analysis):
c
      if(equal(keywd,kortho,5)) then
        ortho = .true.
        goto 20
      end if
c
c  keyword bohr -- atomic coordinates, dipole integrals in bohr:
c
      if(equal(keywd,kbohr,4)) then
        bohr = .true.
        goto 20
      end if
c
c  keyword bodm -- input bond order matrix:
c
      if(equal(keywd,kbodm,4)) then
        denop = .false.
        goto 20
      end if
c
c  keyword ev -- expectation values of the fock operator are in ev:
c
      if(equal(keywd,kev,2)) then
        munit = 1
        goto 20
      end if
c
c  keyword cubicf -- use set of cubic f functions:
c
      if(equal(keywd,kcubf,6)) then
        iwcubf = 1
        goto 20
      end if
c
c  unknown keyword -- halt program:
c
      write(lfnpr,900) keywd
      stop
c
c  end of $gen input encountered, make sure gennbo has all info needed:
c
  700 continue
      if(reuse) then
        newdaf = .false.
        return
      else
        newdaf = .true.
      endif
c
      ndim = nbas
      if(nbas.le.0) stop 'nbas must be specified in $gen keylist'
      if(nbas.gt.maxbas) stop 'increase parameter maxbas'
      if(natoms.le.0) stop 'natoms must be specified in $gen keylist'
      if(natoms.gt.maxatm) stop 'increase parameter maxatm'
      return
c
  900 format(1x,'unrecognized keyword >',6a1,'<')
      end
c*****************************************************************************
      subroutine nboinp(nboopt,idone)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical end,equal,equal2
      dimension nboopt(10)
      dimension keywd(6),knbo(4),krun(4)
c
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      data knbo/1h$,1hn,1hb,1ho/
      data krun/1hr,1hu,1hn,1ht/
c
c  if nboopt(1) = 1, don't search for keywords, just continue with
c  default options:
c
      if(nboopt(1).eq.1) then
        idone = 0
        return
      end if
c
c  if this is the gamess, hondo, or general version of the nbo program,
c  rewind the input file before searching for $nbo:
c
      irep = 1
      if(nboopt(10).eq.0) irep = 0
      if(nboopt(10).eq.6) irep = 0
      if(nboopt(10).eq.7) irep = 0
      if(irep.eq.0) rewind(lfnin)
c
c search input file for runtype ...
c
      icount=0
   3  call strtin(lfnin)
      len=4
      call hfld(keywd,len,end)
      if (equal(keywd,krun,4)) goto 10
      if(len.eq.0.and.end) goto 60
      goto 3
c
c  search input file for $nbo: (it is the second occurrence..)
c
   10 call strtin(lfnin)
      len = 6
      call hfld(keywd,len,end)
      if(equal2(keywd,knbo,4)) then
         icount=icount+1
         if (icount.eq.2) goto 50
      endif
      if(len.eq.0.and.end) goto 60
      goto 10
c
c  $nbo found -- continue with the nbo analysis:
c
   50 continue
      idone = 0
      return
c
c  end of file encountered -- stop nbo analysis, except for the general
c  version of the program (set nboopt(1) so keywords are not read):
c
   60 continue
      if(irep.eq.1) then
        rewind(lfnin)
        irep = irep + 1
        goto 10
      else if(nboopt(10).eq.0) then
        nboopt(1) = 1
        idone = 0
      else
        idone = 1
      end if
      return
      end
c*****************************************************************************
      subroutine corinp(iess,icor)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical end,equal2
c
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension keywd(6),kcor(4),kchs(4),kdel(4),knbo(4),knrt(4)
c
      data kcor/1h$,1hc,1ho,1hr/,kchs/1h$,1hc,1hh,1ho/,
     +     kdel/1h$,1hd,1he,1hl/,knbo/1h$,1hn,1hb,1ho/,
     +     knrt/1h$,1hn,1hr,1ht/
c
c  if icor is set to -1, do not read in the $core keylist:
c
      if(icor.eq.-1) return
c
c  if this is the gamess, hondo, or general version of the nbo program,
c  rewind the input file before searching for $core:
c
      irep = 1
      if(iess.eq.0) irep = 0
      if(iess.eq.6) irep = 0
      if(iess.eq.7) irep = 0
      if(irep.eq.0) rewind(lfnin)
c
c  search input file for $core:
c
   10 call strtin(lfnin)
      len = 6
      call hfld(keywd,len,end)
      if(equal2(keywd,kcor,4)) goto 50
      if(equal2(keywd,knbo,4)) goto 60
      if(equal2(keywd,kchs,4)) goto 60
      if(equal2(keywd,kdel,4)) goto 60
      if(equal2(keywd,knrt,4)) goto 60
      if(len.eq.0.and.end) goto 70
      goto 10
c
c  $core found:
c
   50 continue
      icor = 1
      return
c
c  $nbo, $choose, $del -- discontinue the search for $core (gaussian, ampac)
c        or $nrt          continue searching for $core (gennbo, gamess, hondo)
c
   60 continue
      if(irep.eq.0) goto 10
      backspace(lfnin)
      icor = 0
      return
c
c  end of file encountered:
c
   70 continue
      icor = 0
      return
      end
c*****************************************************************************
      subroutine chsinp(iess,ichs)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical end,equal2
c
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension keywd(6),kchs(4),kdel(4),knbo(4),knrt(4)
c
      data kchs/1h$,1hc,1hh,1ho/,kdel/1h$,1hd,1he,1hl/,
     +     knbo/1h$,1hn,1hb,1ho/,knrt/1h$,1hn,1hr,1ht/
c
c  if ichs is set to -1, do not search for the $choose keylist:
c
      if(ichs.eq.-1) return
c
c  if this is the gamess, hondo, or general version of the nbo program,
c  rewind the input file before searching for $choose:
c
      irep = 1
      if(iess.eq.0) irep = 0
      if(iess.eq.6) irep = 0
      if(iess.eq.7) irep = 0
      if(irep.eq.0) rewind(lfnin)
c
c  search input file for $choose:
c
   10 call strtin(lfnin)
      len = 6
      call hfld(keywd,len,end)
      if(equal2(keywd,kchs,4)) goto 50
      if(equal2(keywd,knbo,4)) goto 60
      if(equal2(keywd,kdel,4)) goto 60
      if(equal2(keywd,knrt,4)) goto 60
      if(len.eq.0.and.end) goto 70
      goto 10
c
c  $choose found:
c
   50 continue
      ichs = 1
      return
c
c  $nbo, $del found -- discontinue the search for $choose (gaussian, ampac)
c      or $nrt         continue searching for $choose (gennbo, gamess, hondo)
c
   60 continue
      if(irep.eq.0) goto 10
      backspace(lfnin)
      ichs = 0
      return
c
c  end of file encountered:
c
   70 continue
      ichs = 0
      return
      end
c*****************************************************************************
      subroutine delinp(nboopt,idone)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical end,equal2
      dimension nboopt(10)
      dimension keywd(6),kdel(4),knbo(4)
c
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      data kdel/1h$,1hd,1he,1hl/,knbo/1h$,1hn,1hb,1ho/
c
c  if this is the gamess, hondo, or general version of the nbo program,
c  rewind the input file before searching for $del:
c
      irep = 1
      if(nboopt(10).eq.0) irep = 0
      if(nboopt(10).eq.6) irep = 0
      if(nboopt(10).eq.7) irep = 0
      if(irep.eq.0) rewind(lfnin)
c
c  search input file for $del:
c
   10 call strtin(lfnin)
      len = 6
      call hfld(keywd,len,end)
      if(equal2(keywd,kdel,4)) goto 50
      if(equal2(keywd,knbo,4)) goto 60
      if(len.eq.0.and.end) goto 70
      goto 10
c
c  $del found -- continue with the nbo energetic analysis:
c
   50 continue
      idone = 0
      return
c
c  $nbo found -- discontinue the search for $del (gaussian, ampac)
c                continue searching for $del (gennbo, gamess, hondo)
c
   60 continue
      if(irep.eq.0) goto 10
      backspace(lfnin)
      idone = 1
      return
c
c  end of file encountered -- stop nbo energetic analysis
c
   70 continue
      if(irep.eq.1) then
        rewind(lfnin)
        irep = irep + 1
        goto 10
      else
        idone = 1
      end if
      return
      end
c*****************************************************************************
      subroutine rdcore(jcore)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical error
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
c  initialize the atomic core array:
c
      do 10 i = 1,natoms
        iatcr(i) = -1
   10 continue
c
c  read in modifications to the nominal core table:
c
      if(jcore.eq.1) then
        write(lfnpr,900)
   30   call ifld(ii,error)
        if(error) goto 40
        if(ii.lt.1.or.ii.gt.natoms) goto 810
        call ifld(jj,error)
        if(error) goto 820
        if(jj.lt.0) goto 830
        iatcr(ii) = jj
        goto 30
      end if
   40 continue
      return
c
  810 write(lfnpr,910) ii
      stop
c
  820 write(lfnpr,920) ii
      stop
c
  830 write(lfnpr,930) jj,ii
      stop
c
  900 format(/1x,'modified core list read from the $core keylist')
  910 format(/1x,'atom ',i4,' not found on this molecule')
  920 format(/1x,'no core orbitals selected for atom ',i4)
  930 format(/1x,i4,' core orbitals on atom ',i4,' does not make sense')
      end
c*****************************************************************************
      subroutine wrppna(t,occ,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
c
      dimension t(ndim,ndim),occ(ndim)
      character*80 title
c
c  write the pnao information to the external file abs(iflg):
c
c  note: this is the pure-ao to pnao transformation, not the raw ao
c        to pnao transform.
c
      title = 'pnaos in the pao basis:'
      call aout(t,ndim,nbas,nbas,title,-1,iflg)
c
c  write the nao orbital labels to the external file:
c
      lfn = abs(iflg)
      write(lfn,900) (naoctr(j),j=1,nbas)
      write(lfn,900) (naol(j),j=1,nbas)
      write(lfn,900) (lstocc(j),j=1,nbas)
c
c  write the pnao orbital occupancies:
c
      write(lfn,910) (occ(j),j=1,nbas)
      return
c
  900 format(1x,20i4)
  910 format(1x,5f15.9)
      end
c*****************************************************************************
      subroutine rdppna(t,occ,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension t(ndim,ndim),occ(ndim)
      dimension job(20)
      logical error
c
c  read the pnao information from the external file abs(iflg/1000)
c
c  note: this is the pure-ao to pnao transformation, not the raw ao
c        to pnao transform.
c
      lfn = abs(iflg/1000)
      write(lfnpr,900)
c
      if(ispin.ge.0) rewind(lfn)
      call aread(t,ndim,nbas,nbas,job,lfn,error)
      if(error) goto 800
      if(ispin.ge.0) write(lfnpr,910) job
      if(ispin.lt.0) write(lfnpr,920)
c
c  read in orbital labels from lfn:
c
      read(lfn,1000,end=810) (naoctr(j),j=1,nbas)
      read(lfn,1000,end=810) (naol(j),j=1,nbas)
      read(lfn,1000,end=810) (lstocc(j),j=1,nbas)
c
c  read orbital occupancies:
c
      read(lfn,1010,end=820) (occ(j),j=1,nbas)
      return
c
  800 write(lfnpr,950) lfn
      stop
c
  810 write(lfnpr,960) lfn
      stop
c
  820 write(lfnpr,970) lfn
      stop
c
  900 format(/1x,'pnao basis set from a previous calculation used:')
  910 format(1x,20a4)
  920 format(/1x,'see alpha nbo output for title of the transformation')
  950 format(/1x,'error reading pao to pnao transformation from lfn',i3)
  960 format(/1x,'error reading pnao orbital labels from lfn',i3)
  970 format(/1x,'error reading pnao orbital occupancies from lfn',i3)
 1000 format(1x,20i4)
 1010 format(1x,5f15.9)
      end
c*****************************************************************************
      subroutine wrtnao(t,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
      dimension t(ndim,ndim)
      character*80 title
c
c  note: t is the pnao overlap matrix on return to the calling routine
c
c  fetch the ao to nao transformation from the nbo daf, and write
c  it to the external file abs(iflg):
c
      call fetnao(t)
      title = 'naos in the ao basis:'
      call aout(t,ndim,nbas,nbas,title,1,iflg)
c
c  write the nao orbital labels to the external file:
c
      lfn = abs(iflg)
      write(lfn,900) (naoctr(j),j=1,nbas)
      write(lfn,900) (naol(j),j=1,nbas)
      write(lfn,900) (lstocc(j),j=1,nbas)
c
c  fetch the pnao overlap matrix from the nbo daf, and store only the
c  upper triangular portion on the external file:
c
      call fesnao(t)
      title = 'pnao overlap matrix:'
      call aout(t,ndim,-nbas,nbas,title,2,iflg)
      return
c
  900 format(1x,20i4)
      end
c*****************************************************************************
      subroutine rdtnao(dm,t,scr,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),lstemt(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension t(ndim,ndim),dm(ndim,ndim),scr(ndim)
      dimension job(20)
      logical error
c
c  note: t is the pnao overlap matrix on return to the calling routine
c        dm is the nao density matrix on return
c
c  read in ao to nao transformation from the external file abs(iflg/1000),
c  and store it on the nbo daf:
c
      lfn = abs(iflg/1000)
      write(lfnpr,900)
c
      rewind(lfn)
      call aread(t,ndim,nbas,nbas,job,lfn,error)
      if(error) goto 800
      write(lfnpr,910) job
      call svtnao(t)
c
c  transform the ao density matrix, presently in dm, to the nao basis:
c
      call simtrs(dm,t,scr,ndim,nbas)
c
c  read in orbital labels from lfn:
c
      read(lfn,1000,end=810) (naoctr(j),j=1,nbas)
      read(lfn,1000,end=810) (naol(j),j=1,nbas)
      read(lfn,1000,end=810) (lstocc(j),j=1,nbas)
c
c  read the pnao overlap from lfn, and save this matrix on the nbo daf:
c
      call aread(t,ndim,-nbas,nbas,job,lfn,error)
      if(error) goto 820
      call svsnao(t)
      return
c
  800 write(lfnpr,950) lfn
      stop
c
  810 write(lfnpr,960) lfn
      stop
c
  820 write(lfnpr,970) lfn
      stop
c
  900 format(/1x,'nao basis set from a previous calculation used:')
  910 format(1x,20a4)
  950 format(/1x,'error reading ao to nao transformation from lfn',i3)
  960 format(/1x,'error reading nao orbital labels from lfn',i3)
  970 format(/1x,'error reading pnao overlap matrix from lfn',i3)
 1000 format(1x,20i4)
      end
c*****************************************************************************
      subroutine wrtnab(t,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
c
      dimension t(ndim,ndim)
      character*80 title
c
c  write the nao to nbo transformation and nbo info to external file
c  abs(iflg):
c
      title = 'nbos in the nao basis:'
      call aout(t,ndim,nbas,nbas,title,2,iflg)
c
c  write the nbo labels:
c
      lfn = abs(iflg)
      do 10 i = 1,nbas
        write(lfn,900) (label(i,j),j=1,6),ibxm(i)
   10 continue
      return
c
  900 format(1x,a2,a1,4i3,3x,i3)
      end
c*****************************************************************************
      subroutine rdtnab(t,dm,bndocc,scr,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension t(ndim,ndim),dm(ndim,ndim),bndocc(ndim),scr(ndim)
      dimension job(20)
      logical error
c
c  read the nao to nbo transformation matrix from the external file
c  abs(iflg/1000).  also read the nbo labels, the nbo occupancies,
c  and transform the input nao density matrix to the nbo basis:
c
      lfn = abs(iflg/1000)
      write(lfnpr,900)
c
      if(ispin.ge.0) rewind(lfn)
      call aread(t,ndim,nbas,nbas,job,lfn,error)
      if(error) goto 800
      if(ispin.ge.0) write(lfnpr,910) job
      if(ispin.lt.0) write(lfnpr,920)
c
c  read the nbo labels:
c
      do 10 i = 1,nbas
        read(lfn,1000,end=810) (label(i,j),j=1,6),ibxm(i)
   10 continue
c
c  transform the nao density matrix, dm, to the nbo basis, and store the
c  nbo occupancies in bndocc:
c
      call simtrs(dm,t,scr,ndim,nbas)
      do 20 i = 1,nbas
        bndocc(i) = dm(i,i)
   20 continue
      return
c
  800 write(lfnpr,950) lfn
      stop
c
  810 write(lfnpr,960) lfn
      stop
c
  900 format(/1x,'nao to nbo transformation from a previous ',
     + 'calculation will be used:')
  910 format(1x,20a4)
  920 format(/1x,'see alpha nbo output for title of the transformation')
  950 format(/1x,'error reading nao to nbo transformation from lfn',i3)
  960 format(/1x,'error reading nbo orbital labels from lfn',i3)
 1000 format(1x,a2,a1,4i3,3x,i3)
      end
c*****************************************************************************
      subroutine wrtnbo(t,bndocc,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension t(ndim,ndim),bndocc(1)
      character*80 title
c
c  write the ao to nbo transformation matrix and nbo info to the external
c  file abs(iflg):
c
      title = 'nbos in the ao basis:'
      call aout(t,ndim,nbas,nbas,title,1,iflg)
c
c  write out the nbo occupancies:
c
      lfn = abs(iflg)
      write(lfn,900) (bndocc(j),j=1,nbas)
c
c  write out nbouni, nbotyp, label, ibxm, and iatno:
c
      write(lfn,910) (nbouni(j),j=1,nbas)
      write(lfn,910) (nbotyp(j),j=1,nbas)
      write(lfn,920) (label(j,1),j=1,nbas)
      write(lfn,920) (label(j,2),j=1,nbas)
      write(lfn,910) (label(j,3),j=1,nbas)
      write(lfn,910) (label(j,4),j=1,nbas)
      write(lfn,910) (label(j,5),j=1,nbas)
      write(lfn,910) (label(j,6),j=1,nbas)
      write(lfn,910) (ibxm(j),j=1,nbas)
      write(lfn,910) (iatno(j),j=1,natoms)
      return
c
  900 format(1x,5f15.9)
  910 format(1x,20i3)
  920 format(1x,20a3)
      end
c*****************************************************************************
      subroutine wrnlmo(t,dm,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),lbl(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension t(ndim,ndim),dm(ndim,ndim)
      character*80 title
c
c  write the ao to nlmo transformation matrix and nlmo info to the external
c  file abs(iflg):
c
      title = 'nlmos in the ao basis:'
      call aout(t,ndim,nbas,nbas,title,1,iflg)
c
c  write out the nlmo occupancies:
c
      lfn = abs(iflg)
      write(lfn,900) (dm(j,j),j=1,nbas)
c
c  write out nbouni, nbotyp, label, ibxm, and iatno:
c
      write(lfn,910) (nbouni(j),j=1,nbas)
      write(lfn,910) (nbotyp(j),j=1,nbas)
      write(lfn,920) (label(j,1),j=1,nbas)
      write(lfn,920) (label(j,2),j=1,nbas)
      write(lfn,910) (label(j,3),j=1,nbas)
      write(lfn,910) (label(j,4),j=1,nbas)
      write(lfn,910) (label(j,5),j=1,nbas)
      write(lfn,910) (label(j,6),j=1,nbas)
      write(lfn,910) (ibxm(j),j=1,nbas)
      write(lfn,910) (iatno(j),j=1,natoms)
      return
c
  900 format(1x,5f15.9)
  910 format(1x,20i3)
  920 format(1x,20a3)
      end
c*****************************************************************************
      subroutine wrbas(scr,iscr,lfn)
c*****************************************************************************
c
c  save the ao basis set information on an external file:
c
c-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbao/lctr(maxbas),lang(maxbas)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension scr(1),iscr(1)
c
c  fetch the number of shells nshell, the number of exponents nexp,
c  the ncomp, nprim, and nptr arrays, and the orbital exponents and
c  coefficients from the nbo daf:
c
      call febas(nshell,nexp,iscr)
c
c  if nshell is zero, then no basis set info has been stored in the
c  daf:
c
      if(nshell.eq.0) then
        write(lfnpr,900)
        return
      end if
c
c  partition the scratch arrays:  (note that scr and iscr occupy the same
c  space in memory)
c
c  iscr: (integer)
c
c   nshell  nexp   ncomp   nprim   nptr
c  +------+------+-------+-------+-------+-----------------------------------
c                i1      i2      i3
c
c  scr: (real)
c                                                                     atcoor
c                                           exp   cs   cp   cd   cf   title
c  ---------------------------------------+-----+----+----+----+----+--------
c                                         i4    i5   i6   i7   i8   i9
c
c  iscr(i1) : ncomp(1..nshell)
c  iscr(i2) : nprim(1..nshell)
c  iscr(i3) : nptr(1..nshell)
c  scr(i4)  : exp(1..nexp)
c  scr(i5)  : cs(1..nexp)
c  scr(i6)  : cp(1..nexp)
c  scr(i7)  : cd(1..nexp)
c  scr(i8)  : cf(1..nexp)
c  scr(i9)  : title(10) or atcoor(3*natoms)
c
      i1   = 3
      i2   = i1 + nshell
      i3   = i2 + nshell
      i4   = i3 + nshell
      i5   = i4 + nexp
      i6   = i5 + nexp
      i7   = i6 + nexp
      i8   = i7 + nexp
      i9   = i8 + nexp
c     iend = i9 + max(3*natoms,10)
c
c  fetch job title and write it to the aoinfo external file:
c
      call fetitl(scr(i9))
c
c  begin writing to the aoinfo external file:
c
      write(lfn,910) (scr(i9+i),i=0,9)
      write(lfn,920) natoms,nshell,nexp
c
c  fetch the atomic coordinates:
c
      call fecoor(scr(i9))
c
c  write atomic numbers and coordinates to external file:
c
      j = 0
      do 10 i = 1,natoms
        write(lfn,930) iatno(i),(scr(i9+j+k),k=0,2)
        j = j + 3
   10 continue
      write(lfn,940)
c
c  write out information about each shell in the basis set:
c
c     nctr(i)  --  atomic center of the ith shell
c
c     ncomp(i) --  number of components in the ith shell
c
c     nptr(i)  --  pointer for the ith shell into the primitive parameters
c                  of exp, cs, cp, cd, and cf
c
c     nprim(i) --  number of primitive functions in the ith shell
c
c     label(1..ncomp(i)) -- symmetry labels for the orbitals of this shell
c
      j1 = 1
      j2 = i1
      j3 = i3
      j4 = i2
      do 20 i = 1,nshell
        ncomp = iscr(j2)
        nprim = iscr(j3)
        nptr  = iscr(j4)
        write(lfn,950) lctr(j1),ncomp,nprim,nptr
        write(lfn,950) ((lang(j1+j)),j=0,ncomp-1)
        j1 = j1 + ncomp
        j2 = j2 + 1
        j3 = j3 + 1
        j4 = j4 + 1
   20 continue
      write(lfn,940)
c
c  write out the primitive parameters:
c
      write(lfn,960) (scr(i4+i),i=0,nexp-1)
      write(lfn,970)
      write(lfn,960) (scr(i5+i),i=0,nexp-1)
      write(lfn,970)
      write(lfn,960) (scr(i6+i),i=0,nexp-1)
      write(lfn,970)
      write(lfn,960) (scr(i7+i),i=0,nexp-1)
      write(lfn,970)
      write(lfn,960) (scr(i8+i),i=0,nexp-1)
      return
c
  900 format(/1x,'no basis set information is stored on the nbo direct',
     + ' access file.',/1x,'thus, no aoinfo file can be written.')
  910 format(1x,9a8,a7,/1x,'basis set information needed for plotting ',
     + 'orbitals',/1x,75('-'))
  920 format(1x,3i6,/1x,75('-'))
  930 format(1x,i4,3(2x,f12.9))
  940 format(1x,75('-'))
  950 format(1x,10i6)
  960 format(2x,4e18.9)
  970 format(1x)
      end
c*****************************************************************************
      subroutine wrarc(scr,iscr,lfn)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      parameter (maxd = 4)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbao/lctr(maxbas),lang(maxbas)
c
      dimension scr(1),iscr(1),ik(maxd)
      dimension kgen(7),knat(6),kbas(4),kopen(4),kortho(5),kupper(5),
     + kbodm(4),kev(2),kcubf(6),kend(4),kcal(4)
c
      data kgen/1h$,1hg,1he,1hn,1hn,1hb,1ho/,kbas/1hn,1hb,1ha,1hs/,
     +     knat/1hn,1ha,1ht,1ho,1hm,1hs/,kopen/1ho,1hp,1he,1hn/,
     +     kortho/1ho,1hr,1ht,1hh,1ho/,kupper/1hu,1hp,1hp,1he,1hr/,
     +     kbodm/1hb,1ho,1hd,1hm/,kev/1he,1hv/,kend/1h$,1he,1hn,1hd/,
     +     kcubf/1hc,1hu,1hb,1hi,1hc,1hf/,kcal/1hk,1hc,1ha,1hl/
      data kblnk,keq/1h ,1h=/
      data ablnks,acentr,alabel/8h        ,8hcenter =,8h label =/
      data anshll,anexp ,ancomp/8hnshell =,8h  nexp =,8h ncomp =/
      data anprim,anptr ,aexp  /8h nprim =,8h  nptr =,8h   exp =/
      data acs,acp,acd,acf/8h    cs =,8h    cp =,8h    cd =,8h    cf =/
      data zero/0.0d0/
c
c  write the archive file to lfn:
c
c  this routine has been written assuming nbas = ndim.  skip if
c  this condition is not satisfied:
c
      if(nbas.ne.ndim) then
        write(lfnpr,890)
        return
      end if
c
c  form the $gennbo keylist in iscr:
c
      nc = 0
      do 10 i = 1,7
        nc = nc + 1
        iscr(nc) = kgen(i)
   10 continue
      nc = nc + 1
      iscr(nc) = kblnk
      nc = nc + 1
      iscr(nc) = kblnk
c
c  add the number of atoms and basis functions:
c
      do 20 i = 1,6
        nc = nc + 1
        iscr(nc) = knat(i)
   20 continue
      nc = nc + 1
      iscr(nc) = keq
      call idigit(natoms,ik,nd,maxd)
      do 30 i = 1,nd
        nc = nc + 1
        iscr(nc) = ik(i)
   30 continue
      nc = nc + 1
      iscr(nc) = kblnk
      nc = nc + 1
      iscr(nc) = kblnk
c
      do 40 i = 1,4
        nc = nc + 1
        iscr(nc) = kbas(i)
   40 continue
      nc = nc + 1
      iscr(nc) = keq
      call idigit(nbas,ik,nd,maxd)
      do 50 i = 1,nd
        nc = nc + 1
        iscr(nc) = ik(i)
   50 continue
      nc = nc + 1
      iscr(nc) = kblnk
      nc = nc + 1
      iscr(nc) = kblnk
c
c  if open shell, add the open keyword:
c
      if(open) then
        do 60 i = 1,4
          nc = nc + 1
          iscr(nc) = kopen(i)
   60   continue
        nc = nc + 1
        iscr(nc) = kblnk
        nc = nc + 1
        iscr(nc) = kblnk
      end if
c
c  if the ao basis is orthogonal, add the ortho keyword:
c
      if(ortho) then
        do 70 i = 1,5
          nc = nc + 1
          iscr(nc) = kortho(i)
   70   continue
        nc = nc + 1
        iscr(nc) = kblnk
        nc = nc + 1
        iscr(nc) = kblnk
      end if
c
c  only upper triangular portions of symmetric matrices will be given:
c
      do 80 i = 1,5
        nc = nc + 1
        iscr(nc) = kupper(i)
   80 continue
      nc = nc + 1
      iscr(nc) = kblnk
      nc = nc + 1
      iscr(nc) = kblnk
c
c  enter the bond-order matrix, bodm, if possible:
c
      if(iwdm.eq.1) then
        do 90 i = 1,4
          nc = nc + 1
          iscr(nc) = kbodm(i)
   90   continue
        nc = nc + 1
        iscr(nc) = kblnk
        nc = nc + 1
        iscr(nc) = kblnk
      end if
c
c  add ev if the energy units are in electron volts:
c
      if(munit.eq.1) then
        nc = nc + 1
        iscr(nc) = kev(1)
        nc = nc + 1
        iscr(nc) = kev(2)
        nc = nc + 1
        iscr(nc) = kblnk
        nc = nc + 1
        iscr(nc) = kblnk
      end if
c
c  add kcal if the energy units are in kcal/mol:
c
      if(munit.eq.1) then
        nc = nc + 1
        iscr(nc) = kcal(1)
        nc = nc + 1
        iscr(nc) = kcal(2)
        nc = nc + 1
        iscr(nc) = kcal(3)
        nc = nc + 1
        iscr(nc) = kcal(4)
        nc = nc + 1
        iscr(nc) = kblnk
        nc = nc + 1
        iscr(nc) = kblnk
      end if
c
c  add cubicf if these types of orbitals are being used:
c
      if(iwcubf.ne.0) then
        do 100 i = 1,6
          nc = nc + 1
          iscr(nc) = kcubf(i)
  100   continue
        nc = nc + 1
        iscr(nc) = kblnk
        nc = nc + 1
        iscr(nc) = kblnk
      end if
c
c  add $end:
c
      do 110 i = 1,4
        nc = nc + 1
        iscr(nc) = kend(i)
  110 continue
c
c  write the $gennbo keylist to the archive file:
c
      write(lfn,900) (iscr(i),i=1,nc)
c
c  write the $nbo keylist to the archive file:
c
      write(lfn,910)
c
c  write the $coord data list to the archive file:
c
      write(lfn,920)
      call fetitl(scr)
      write(lfn,930) (scr(i),i=1,10)
      call fecoor(scr)
      j = 1
      do 120 i = 1,natoms
        write(lfn,940) iatno(i),iznuc(i),scr(j),scr(j+1),scr(j+2)
        j = j + 3
  120 continue
      write(lfn,950)
c
c  write the $basis datalist to the archive file (info from /nbao/):
c
      write(lfn,960)
      nint = 17
      str = acentr
      do 130 i = 1,(nbas-1)/nint+1
        nl  = (i - 1) * nint + 1
        nu  = min(nl+nint-1,nbas)
        write(lfn,970) str,(lctr(j),j=nl,nu)
        str = ablnks
  130 continue
      str = alabel
      do 140 i = 1,(nbas-1)/nint+1
        nl  = (i - 1) * nint + 1
        nu  = min(nl+nint-1,nbas)
        write(lfn,970) str,(lang(j),j=nl,nu)
        str = ablnks
  140 continue
      write(lfn,950)
c
c  write the $contract datalist to the archive file:
c
c  fetch the basis set info from the nbo daf:
c
      call febas(nshell,nexp,iscr)
c
c  partition the scratch vector:
c
c  iscr(i1) : ncomp(1..nshell)
c  iscr(i2) : nprim(1..nshell)
c  iscr(i3) : nptr(1..nshell)
c  scr(i4)  : exp(1..nexp)
c  scr(i5)  : cs(1..nexp)
c  scr(i6)  : cp(1..nexp)
c  scr(i7)  : cd(1..nexp)
c  scr(i8)  : cf(1..nexp)
c
      i1   = 3
      i2   = i1 + nshell
      i3   = i2 + nshell
      i4   = i3 + nshell
      i5   = i4 + nexp
      i6   = i5 + nexp
      i7   = i6 + nexp
      i8   = i7 + nexp
c     iend = i8 + nexp
c
c  if nshell is zero, then no basis set info was ever stored on
c  the daf:
c
      if(nshell.gt.0) then
c
c  write out numbers of shells and orbital exponents:
c
        write(lfn,980)
        write(lfn,970) anshll,nshell
        write(lfn,970) anexp,nexp
c
c  write out the number of components in each shell:
c
        nint = 17
        str = ancomp
        do 150 i = 1,(nshell-1)/nint+1
          nl  = (i - 1) * nint + 1
          nu  = min(nl+nint-1,nshell)
          write(lfn,970) str,(iscr(j),j=i1+nl-1,i1+nu-1)
          str = ablnks
  150   continue
c
c  write out the number of primitives in each shell:
c
        str = anprim
        do 160 i = 1,(nshell-1)/nint+1
          nl  = (i - 1) * nint + 1
          nu  = min(nl+nint-1,nshell)
          write(lfn,970) str,(iscr(j),j=i2+nl-1,i2+nu-1)
          str = ablnks
  160   continue
c
c  write out pointer array which maps orbital exponents and coefficients
c  onto each shell:
c
        str = anptr
        do 170 i = 1,(nshell-1)/nint+1
          nl  = (i - 1) * nint + 1
          nu  = min(nl+nint-1,nshell)
          write(lfn,970) str,(iscr(j),j=i3+nl-1,i3+nu-1)
          str = ablnks
  170   continue
c
c  write out orbital exponents:
c
        nreal = 4
        str   = aexp
        do 180 i = 1,(nexp-1)/nreal+1
          nl  = (i - 1) * nreal + 1
          nu  = min(nl+nreal-1,nexp)
          write(lfn,990) str,(scr(j),j=i4+nl-1,i4+nu-1)
          str = ablnks
  180   continue
c
c  write out the orbital coefficients for each angular symmetry type
c  unless there are no basis functions of that type:
c
        do 210 i = 1,4
          if(i.eq.1) then
            str = acs
            ii  = i5
          else if(i.eq.2) then
            str = acp
            ii  = i6
          else if(i.eq.3) then
            str = acd
            ii  = i7
          else if(i.eq.4) then
            str = acf
            ii  = i8
          end if
          iflg = 0
          do 190 j = ii,ii+nexp-1
            if(scr(j).ne.zero) iflg = 1
  190     continue
          if(iflg.eq.1) then
            do 200 j = 1,(nexp-1)/nreal+1
              nl  = (j - 1) * nreal + 1
              nu  = min(nl+nreal-1,nexp)
              write(lfn,990) str,(scr(k),k=ii+nl-1,ii+nu-1)
              str = ablnks
  200       continue
          end if
  210   continue
        write(lfn,950)
      end if
c
c  write the $overlap datalist unless the ao basis is orthogonal:
c
      l2 = ndim * (ndim + 1) / 2
      if(.not.ortho) then
        write(lfn,1000)
        call fesraw(scr)
        l2 = ndim * (ndim + 1) / 2
        call NBOpack(scr,ndim,nbas,l2)
        write(lfn,1010) (scr(i),i=1,l2)
        write(lfn,950)
      end if
c
c  write the $density datalist:
c
      write(lfn,1020)
      if(open) then
        alpha = .true.
        beta  = .false.
        call fedraw(scr,scr)
        call NBOpack(scr,ndim,nbas,l2)
        write(lfn,1010) (scr(i),i=1,l2)
        alpha = .false.
        beta  = .true.
        call fedraw(scr,scr)
        call NBOpack(scr,ndim,nbas,l2)
        write(lfn,1010) (scr(i),i=1,l2)
      else
        alpha = .false.
        beta  = .false.
        call fedraw(scr,scr)
        call NBOpack(scr,ndim,nbas,l2)
        write(lfn,1010) (scr(i),i=1,l2)
      end if
      write(lfn,950)
c
c  write the $fock datalist:
c
      if(open) then
        alpha = .true.
        beta  = .false.
        iwfock = 1
        call fefao(scr,iwfock)
        if(iwfock.ne.0) then
          write(lfn,1030)
          call NBOpack(scr,ndim,nbas,l2)
          write(lfn,1010) (scr(i),i=1,l2)
          alpha = .false.
          beta  = .true.
          call fefao(scr,iwfock)
          call NBOpack(scr,ndim,nbas,l2)
          write(lfn,1010) (scr(i),i=1,l2)
          write(lfn,950)
        end if
      else
        alpha = .false.
        beta  = .false.
        iwfock = 1
        call fefao(scr,iwfock)
        if(iwfock.ne.0) then
          write(lfn,1030)
          call NBOpack(scr,ndim,nbas,l2)
          write(lfn,1010) (scr(i),i=1,l2)
          write(lfn,950)
        end if
      end if
c
c  write the $lcaomo datalist:
c
      if(open) then
        alpha = .true.
        beta  = .false.
        call feaomo(scr,iaomo)
        if(iaomo.eq.1) then
          write(lfn,1040)
          write(lfn,1010) (scr(i),i=1,ndim*ndim)
          alpha = .false.
          beta  = .true.
          call feaomo(scr,iaomo)
          write(lfn,1010) (scr(i),i=1,ndim*ndim)
          write(lfn,950)
        end if
      else
        alpha = .false.
        beta  = .false.
        call feaomo(scr,iaomo)
        if(iaomo.eq.1) then
          write(lfn,1040)
          write(lfn,1010) (scr(i),i=1,ndim*ndim)
          write(lfn,950)
        end if
      end if
c
c  write the $dipole datalist:
c
      idip = 1
      call fedxyz(scr,idip)
      if(idip.ne.0) then
        write(lfn,1050)
        call NBOpack(scr,ndim,nbas,l2)
        write(lfn,1010) (scr(i),i=1,l2)
        idip = 2
        call fedxyz(scr,idip)
        call NBOpack(scr,ndim,nbas,l2)
        write(lfn,1010) (scr(i),i=1,l2)
        idip = 3
        call fedxyz(scr,idip)
        call NBOpack(scr,ndim,nbas,l2)
        write(lfn,1010) (scr(i),i=1,l2)
        write(lfn,950)
      end if
c
c  reset logicals alpha and beta:
c
      alpha = ispin.eq.2
      beta  = ispin.eq.-2
      return
c
  890 format(/1x,'the routine which writes the archive file assumes ',
     + 'nbas = ndim.  since',/1x,'this condition is not satisfied, ',
     + 'the archive file will not be written.')
  900 format(1x,78a1)
  910 format(1x,'$nbo  $end')
  920 format(1x,'$coord')
  930 format(1x,9a8,a6)
  940 format(1x,2i5,3f15.6)
  950 format(1x,'$end')
  960 format(1x,'$basis')
  970 format(1x,1x,a8,1x,17(i3,1x))
  980 format(1x,'$contract')
  990 format(1x,1x,a8,1x,4(e15.7,1x))
 1000 format(1x,'$overlap')
 1010 format(1x,1x,5e15.7)
 1020 format(1x,'$density')
 1030 format(1x,'$fock')
 1040 format(1x,'$lcaomo')
 1050 format(1x,'$dipole')
      end
c*****************************************************************************
      subroutine aout(a,mr,nr,nc,title,index,iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension a(mr,1)
      character*80 title
      dimension ishell(4)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nblbl/nlew,nval,lbl(10,maxbas,4)
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
c
      data kfull,kval,klew/4hfull,3hval,3hlew/
c
c  either write a to an external file, or print it in the output file:
c
c  input:  a     -- matrix to be printed or written out
c
c          mr    -- row dimension of matrix a in calling routine
c
c          nr    -- abs(nr) is the actual number of rows to be output
c                   [if nr is negative, iflg is negative (write), and
c                    abs(nr).eq.nc (square matrix), only the upper
c                    triangular portion is written out]
c
c          nc    -- actual number of columns in matrix a
c                   [used to determine if a is square, and as an upper
c                    limit on iflg]
c
c          title -- character*80 variable containing a matrix title
c
c          index -- index selecting appropriate output labels
c                   0 : atom labels
c                   1 : ao   labels
c                   2 : nao  labels
c                   3 : nho  labels
c                   4 : nbo  labels
c                   5 : nlmo labels
c
c          iflg  -- print/write flag
c                   negative : write to lfn abs(iflg)
c                   positive : print iflg columns of a
c                   'full'   : print the full matrix
c                   'val'    : print n columns of a, where n is the
c                              number of core + valence orbitals and
c                              is determined by this routine
c                   'lew'    : print n columns of a, where n is the
c                              number of occupied orbitals and is
c                              determined by this routine
c
      jflg = iflg
      if(jflg.eq.0) return
c
c  if jflg is full, then output the total number of columns:
c
      if(jflg.eq.kfull) jflg = abs(nc)
c
c  if jflg = val, output only the valence orbitals, determined from the
c  core and valence tables:
c
      if(jflg.eq.kval) then
        if(nval.lt.0) then
          iecp = 0
          jflg = 0
          do 30 iat = 1,natoms
            call cortbl(iat,ishell,iecp)
            do 10 i = 1,4
              mult = 2 * (i-1) + 1
              jflg = jflg + ishell(i)*mult
   10       continue
            call valtbl(iat,ishell)
            do 20 i = 1,4
              mult = 2 * (i-1) + 1
              jflg = jflg + ishell(i)*mult
   20       continue
   30     continue
        else
          jflg = nval
        end if
      end if
c
c  if jflg is lew, only output the occupied orbitals:
c
      if(jflg.eq.klew) jflg = nlew
c
c  if jflg is positive, print the matrix a in the output file:
c
      if(jflg.gt.0) call aprint(a,mr,nr,nc,title,index,jflg)
c
c  if jflg is negative but greater than -1000, write matrix a to the external
c  file abs(jflg):
c
      if(jflg.lt.0.and.jflg.gt.-1000) call awrite(a,mr,nr,nc,title,jflg)
c
      return
      end
c*****************************************************************************
      subroutine aprint(a,mr,nr,nc,title,index,mcol)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension a(mr,1)
      character*80 title
      dimension basis(5)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nblbl/nlew,nval,lbl(10,maxbas,4)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      data basis/4h ao ,4h nao,4h nho,4h nbo,4hnlmo/
      data atom,dashes/4hatom,8h--------/
      data tenth/0.1d0/
c
c  determine the number of columns of matrix a to print in the output file:
c
      ncol = mcol
      if(ncol.gt.abs(nc)) ncol = abs(nc)
c
      nn = abs(nr)
      ilabel = index
      if(ilabel.eq.5) ilabel = 4
c
      tmax = abs(a(1,1))
      do 20 j = 1,ncol
        do 10 i = 1,nn
          if(abs(a(i,j)).gt.tmax) tmax = abs(a(i,j))
   10   continue
   20 continue
      if(tmax.lt.tenth) then
        nd = 1
      else
        nd = int(log10(tmax)) + 1
      end if
c
c  print the matrix title:
c
      write(lfnpr,1000) title(1:78)
c
c  print the matrix a: (basis function labels)
c
      if(ilabel.ge.1.and.ilabel.le.4) then
        maxcol = min(10-nd,8)
        if(maxcol.lt.6) then
          call altout(a,mr,ncol,nn,ncol)
        else
          ncl = 1
          ncu = maxcol
          nloops = (ncol - 1) / maxcol + 1
          do 60 l = 1,nloops
            if(ncu.gt.ncol) ncu = ncol
            if(maxcol.eq.8) then
              write(lfnpr,900) basis(index),(j,j=ncl,ncu)
              write(lfnpr,910) (dashes,j=ncl,ncu)
              do 30 i = 1,nn
                write(lfnpr,920) i,(lbl(j,i,ilabel),j=1,10),
     +                           (a(i,k),k=ncl,ncu)
   30         continue
            else if(maxcol.eq.7) then
              write(lfnpr,901) basis(index),(j,j=ncl,ncu)
              write(lfnpr,911) (dashes,j=ncl,ncu)
              do 40 i = 1,nn
                write(lfnpr,921) i,(lbl(j,i,ilabel),j=1,10),
     +                           (a(i,k),k=ncl,ncu)
   40         continue
            else
              write(lfnpr,902) basis(index),(j,j=ncl,ncu)
              write(lfnpr,912) (dashes,dashes,j=ncl,ncu)
              do 50 i = 1,nn
                write(lfnpr,922) i,(lbl(j,i,ilabel),j=1,10),
     +                           (a(i,k),k=ncl,ncu)
   50         continue
            end if
            ncl = ncu + 1
            ncu = ncu + maxcol
   60     continue
        end if
c
c  print the matrix a: (atom labels)
c
      else if(ilabel.eq.0) then
        maxcol = min(10-nd,9)
        if(maxcol.lt.7) then
          call altout(a,mr,ncol,n,ncol)
        else
          ncl = 1
          ncu = maxcol
          nloops = (ncol - 1) / maxcol + 1
          do 160 l = 1,nloops
            if(ncu.gt.ncol) ncu = ncol
            if(maxcol.eq.9) then
              write(lfnpr,1900) atom,(j,j=ncl,ncu)
              write(lfnpr,1910) (dashes,j=ncl,ncu)
              do 130 i = 1,nn
                write(lfnpr,1920) i,nameat(iatno(i)),
     +                            (a(i,k),k=ncl,ncu)
  130         continue
            else if(maxcol.eq.8) then
              write(lfnpr,1901) atom,(j,j=ncl,ncu)
              write(lfnpr,1911) (dashes,j=ncl,ncu)
              do 140 i = 1,nn
                write(lfnpr,1921) i,nameat(iatno(i)),
     +                            (a(i,k),k=ncl,ncu)
  140         continue
            else
              write(lfnpr,1902) atom,(j,j=ncl,ncu)
              write(lfnpr,1912) (dashes,j=ncl,ncu)
              do 150 i = 1,nn
                write(lfnpr,1922) i,nameat(iatno(i)),
     +                            (a(i,k),k=ncl,ncu)
  150         continue
            end if
            ncl = ncu + 1
            ncu = ncu + maxcol
  160     continue
        end if
c
c  print the matrix a: (no labels)
c
      else
        call altout(a,mr,ncol,nn,ncol)
      end if
      return
c
  900 format(/9x,a4,3x,8(3x,i3,2x))
  901 format(/9x,a4,3x,7(4x,i3,2x))
  902 format(/9x,a4,3x,6(4x,i3,3x))
  910 format(6x,'----------',8(1x,a7))
  911 format(6x,'----------',7(1x,a8))
  912 format(6x,'----------',6(1x,a8,a1))
  920 format(1x,i3,'. ',10a1,8f8.4)
  921 format(1x,i3,'. ',10a1,7f9.4)
  922 format(1x,i3,'. ',10a1,6f10.4)
 1000 format(//1x,a78)
 1900 format(/5x,a4,9(2x,i3,3x))
 1901 format(/5x,a4,8(3x,i3,3x))
 1902 format(/5x,a4,7(3x,i3,4x))
 1910 format(5x,'----',1x,9(a6,2x))
 1911 format(5x,'----',1x,8(a7,2x))
 1912 format(5x,'----',1x,7(a8,2x))
 1920 format(1x,i3,'. ',a2,9f8.4)
 1921 format(1x,i3,'. ',a2,8f9.4)
 1922 format(1x,i3,'. ',a2,7f10.4)
      end
c*****************************************************************************
      subroutine awrite(a,mr,nr,nc,title,lfn)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension a(mr,1)
      character*80 title
      dimension xjob(10)
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
c  write the matrix a to the external file abs(lfn).  include job title,
c  matrix title, and specify the spin in needed:
c
      lfnout = abs(lfn)
      if(lfnout.eq.lfnpr) write(lfnout,890)
      if(alpha.or..not.open.or.lfnout.eq.lfnpr) then
        call fetitl(xjob)
        write(lfnout,900) xjob
        write(lfnout,910) title(1:79)
      end if
      if(alpha) write(lfnout,920)
      if(beta)  write(lfnout,930)
c
c  if this is a square matrix and nr is negative, only write the upper
c  triangular portion.  otherwise, write out the full matrix:
c
      if(abs(nr).eq.abs(nc).and.nr.lt.0) then
        write(lfnout,1000) ((a(i,j),i=1,j),j=1,abs(nr))
      else
        do 10 j = 1,abs(nc)
          write(lfnout,1000) (a(i,j),i=1,abs(nr))
   10   continue
      end if
      return
c
  890 format(/1x)
  900 format(1x,9a8,a7)
  910 format(1x,a79,/1x,79('-'))
  920 format(1x,'alpha spin')
  930 format(1x,'beta  spin')
 1000 format(1x,5f15.9)
      end
c*****************************************************************************
      subroutine aread(a,mr,nr,nc,job,lfn,error)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension a(mr,1),job(20)
      dimension itemp(20)
      logical error
c
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
c
      data idash,ialfa,ibeta/4h----,4halph,4hbeta/
c
c  read the matrix a to the external file lfn:
c
c  input:  mr    -- row dimension of matrix a in calling routine
c
c          nr    -- abs(nr) is the actual number of rows to be read
c                   [if nr is negative and abs(nr).eq.nc (square matrix),
c                    only the upper triangular portion is stored in the
c                    input file.  this routine will read the upper triangular
c                    portion and unpack it.]
c
c          nc    -- actual number of columns in matrix a
c                   [used to determine if a is square]
c
c          lfn   -- input file
c
c  output: job   -- integer array containing the job title
c                   [closed shell or alpha spin only]
c
c          error -- set to .true. if the end-of-file was encountered while
c                   reading
c
      if(alpha.or..not.open) read(lfn,1000,end=800) job
      if(.not.open) istr = idash
      if(alpha)     istr = ialfa
      if(beta)      istr = ibeta
c
   10 read(lfn,1000,end=800) itemp
      if(itemp(1).ne.istr) goto 10
c
c  if this is a square matrix and nr is negative, only read the upper
c  triangular portion.  otherwise, read the full matrix:
c
      if(abs(nr).eq.abs(nc).and.nr.lt.0) then
        read(lfn,900,end=800) ((a(i,j),i=1,j),j=1,abs(nr))
        do 30 j = 1,abs(nr)-1
          do 20 i = j+1,abs(nr)
            a(i,j) = a(j,i)
   20     continue
   30   continue
      else
        do 40 j = 1,abs(nc)
          read(lfn,900,end=800) (a(i,j),i=1,abs(nr))
   40   continue
      end if
      error = .false.
      return
c
  800 error = .true.
      return
c
  900 format(1x,5f15.9)
 1000 format(1x,20a4)
      end
c*****************************************************************************
      subroutine altout(a,mr,mc,nr,nc)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      dimension a(mr,mc)
c
c  for 80 column output:
c  list elements of array a (matrix or vector).
c     mr,mc declared row and column dimensionality,
c     nr,nc actual row and column dimensionality,
c
      ncl=1
      ncu=6
      nloops=nc/6+1
      do 20 l=1,nloops
        if(ncu.gt.nc) ncu=nc
        write(lfnpr,1100) (j,j=ncl,ncu)
        do 10 i=1,nr
   10     write(lfnpr,1200) i,(a(i,j),j=ncl,ncu)
        if(ncu.ge.nc) return
        ncl=ncu+1
   20   ncu=ncu+6
      return
 1100 format(/11x,10(i3,9x))
 1200 format(1x,i3,10f12.5)
      end
c*****************************************************************************
      subroutine keypar(string,len,iflg,lfn,read,error)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      integer string(len)
      logical read,error
c
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      data iw,ir,ip,ic,iv,il/1hw,1hr,1hp,1hc,1hv,1hl/
      data kfull,kval,klew/4hfull,3hval,3hlew/
c
c  interpret the keyword parameter string, storing the result in iflg.
c  (the default iflg should be passed to this routine through iflg)
c
c  the following strings are acceptable:
c
c    string = wnnn     means write to the external file nnn (iflg = -nnn)
c                      (if nnn is omitted, iflg = -lfn)
c
c    string = rnnn     means read from the external file nnn (iflg = -nnn*1000)
c                      (if nnn is omitted, iflg = -lfn)
c                      (read must be true to allow reading)
c
c    string = pnnnc    means print nnn columns to the output file (iflg = nnn)
c                      (if nnn is omitted, print full matrix, iflg = 'full')
c                      (the c is optional, it means columns)
c
c    string = pval     means print val columns to output file (iflg = 'val')
c                      (val is the number of core + valence orbitals)
c                      (only the v is necessary)
c
c
c    string = plew     means print lew columns to output file (iflg = 'lew'
c                      (lew is the number of occupied orbitals)
c                      (only the l is necessary)
c
c    string = other    iflg is left untouched
c
      error = .false.
c
c  process string = w..:
c
      if(string(1).eq.iw) then
        if(len.eq.1) then
          iflg = -lfn
          return
        end if
        if(len.gt.1) then
          call convin(string(2),len-1,iflg,error)
          if(error) return
          if(iflg.gt.1000) then
            write(lfnpr,900)
            write(lfnpr,910) iflg
            stop
          end if
          iflg = -iflg
        end if
c
c  process string = r..:
c
      else if(string(1).eq.ir) then
        if(.not.read) then
          error = .true.
          return
        end if
        if(len.eq.1) then
          iflg = -lfn * 1000
          return
        end if
        if(len.gt.1) then
          call convin(string(2),len-1,iflg,error)
          if(error) return
          if(iflg.gt.1000) then
            write(lfnpr,900)
            write(lfnpr,920) iflg
            stop
          end if
          iflg = -iflg * 1000
        end if
c
c  process string = p..:
c
      else if(string(1).eq.ip) then
        if(string(2).eq.iv) then
          iflg = kval
          return
        end if
        if(string(2).eq.il) then
          iflg = klew
          return
        end if
        if(len.eq.1) then
          iflg = kfull
          return
        end if
        if(len.gt.1) then
          if(string(len).ne.ic) then
            call convin(string(2),len-1,iflg,error)
          else
            call convin(string(2),len-2,iflg,error)
          end if
        end if
      else
        error = .true.
      end if
      return
c
  900 format(/1x,'the nbo program will only communicate with external ',
     + 'files 0 thru 999.')
  910 format(1x,'you''re attempting to write to file ',i6,'.')
  920 format(1x,'you''re attempting to read from file ',i6,'.')
      end
c*****************************************************************************
      function ioinqr(iflg)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      data kfull,kval,klew/4hfull,3hval,3hlew/
      data kblnk,kprnt,kwrit,kread/4h    ,4hprnt,4hwrit,4hread/
c
c  interpret iflg, determining whether the corresponding matrix should be
c  printed, written out, or read:
c
      if(iflg.eq.kfull) then
        ioinqr = kprnt
      else if(iflg.eq.kval) then
        ioinqr = kprnt
      else if(iflg.eq.klew) then
        ioinqr = kprnt
      else if(iflg.gt.0) then
        ioinqr = kprnt
      else if(iflg.lt.0.and.iflg.gt.-1000) then
        ioinqr = kwrit
      else if(iflg.lt.0) then
        ioinqr = kread
      else
        ioinqr = kblnk
      end if
      return
      end
c*****************************************************************************
      subroutine lblao
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      parameter(maxd = 2)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbao/lctr(maxbas),lang(maxbas)
      common/nblbl/nlew,nval,iaolbl(10,maxbas),naolbl(10,maxbas),
     +                       nholbl(10,maxbas),nbolbl(10,maxbas)
c
      dimension istr(maxd),iang(5),ixyz(3),ibyte(4),num(10)
c
      data iblnk/1h /
      data iang/1hs,1hp,1hd,1hf,1hg/
      data ixyz/1hx,1hy,1hz/
      data ileft,iright/1h(,1h)/
      data num/1h0,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9/
c
      do 20 iao = 1,nbas
        do 10 i = 1,10
          iaolbl(i,iao) = iblnk
   10   continue
        lbl = nameat(iatno(lctr(iao)))
        call debyte(lbl,ibyte)
        iaolbl(1,iao) = ibyte(1)
        iaolbl(2,iao) = ibyte(2)
        call idigit(lctr(iao),istr,nd,maxd)
        if(nd.eq.1) then
          iaolbl(4,iao) = istr(1)
        else
          iaolbl(3,iao) = istr(1)
          iaolbl(4,iao) = istr(2)
        end if
        iaolbl(6,iao) = ileft
        l = lang(iao)/100
        iaolbl(7,iao) = iang(l+1)
        if(l.eq.0) then
          iaolbl(8,iao) = iright
        else if(l.eq.1) then
          m = mod(lang(iao),10)
          iaolbl(8,iao) = ixyz(m)
          iaolbl(9,iao) = iright
        else if(l.eq.2.or.l.eq.3) then
          iaolbl(8,iao) = num(mod(lang(iao),10)+1)
          iaolbl(9,iao) = iright
        end if
   20 continue
      return
      end
c*****************************************************************************
      subroutine lblnao
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      parameter(maxd = 2)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbnao/naoctr(maxbas),naol(maxbas),ltyp(maxbas),
     +       iprin(maxbas)
      common/nblbl/nlew,nval,iaolbl(10,maxbas),naolbl(10,maxbas),
     +                       nholbl(10,maxbas),nbolbl(10,maxbas)
c
      dimension istr(maxd),iang(5),ixyz(3),ibyte(4),num(10)
c
      data iblnk/1h /
      data iang/1hs,1hp,1hd,1hf,1hg/
      data ixyz/1hx,1hy,1hz/
      data ileft,iright/1h(,1h)/
      data num/1h0,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9/
c
      do 20 inao = 1,nbas
        do 10 i = 1,10
          naolbl(i,inao) = iblnk
   10   continue
        lbl = nameat(iatno(naoctr(inao)))
        call debyte(lbl,ibyte)
        naolbl(1,inao) = ibyte(1)
        naolbl(2,inao) = ibyte(2)
        call idigit(naoctr(inao),istr,nd,maxd)
        if(nd.eq.1) then
          naolbl(4,inao) = istr(1)
        else
          naolbl(3,inao) = istr(1)
          naolbl(4,inao) = istr(2)
        end if
        naolbl(5,inao) = ileft
        call idigit(iprin(inao),istr,nd,maxd)
        if(nd.eq.1) then
          naolbl(7,inao) = istr(1)
        else
          naolbl(6,inao) = istr(1)
          naolbl(7,inao) = istr(2)
        end if
        l = naol(inao)/100
        naolbl(8,inao) = iang(l+1)
        if(l.eq.1) then
          m = mod(naol(inao),10)
          naolbl(9,inao) = ixyz(m)
        else if(l.eq.2.or.l.eq.3) then
          naolbl(9,inao) = num(mod(naol(inao),10)+1)
        end if
        naolbl(10,inao) = iright
   20 continue
      return
      end
c*****************************************************************************
      subroutine lblnbo
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      parameter(maxd = 2)
      integer istr(maxd),ibyte(4)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),lbl1(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nblbl/nlew,nval,iaolbl(10,maxbas),naolbl(10,maxbas),
     +       nholbl(10,maxbas),nbolbl(10,maxbas)
c
      data iblnk,ic,il,ip,ir,iy,istar,ihyp/1h ,1hc,1hl,1hp,1hr,1hy,1h*,
     + 1h-/
      data icr,ilp/2hcr,2hlp/
      data ileft,iright/1h(,1h)/
c
      do 20 inbo = 1,nbas
        do 10 i = 1,10
          nbolbl(i,inbo) = iblnk
   10   continue
        ib = ibxm(inbo)
        nctr = 1
        if(label(ib,5).ne.0) nctr = 2
        if(label(ib,6).ne.0) nctr = 3
c
c  one-center labels:
c
        if(nctr.eq.1) then
          lbl = nameat(iatno(label(ib,4)))
          call debyte(lbl,ibyte)
          nbolbl(1,inbo) = ibyte(1)
          nbolbl(2,inbo) = ibyte(2)
          call idigit(label(ib,4),istr,nd,maxd)
          if(nd.eq.1) then
            nbolbl(4,inbo) = istr(1)
          else
            nbolbl(3,inbo) = istr(1)
            nbolbl(4,inbo) = istr(2)
          end if
          nbolbl(5,inbo) = ileft
          if(label(ib,1).eq.icr) then
            nbolbl(6,inbo) = ic
            nbolbl(7,inbo) = ir
            nbolbl(8,inbo) = iright
          else if(label(ib,1).eq.ilp) then
            nbolbl(6,inbo) = il
            nbolbl(7,inbo) = ip
            if(label(ib,2).eq.istar) then
              nbolbl(8,inbo) = istar
              nbolbl(9,inbo) = iright
            else
              nbolbl(8,inbo) = iright
            end if
          else
            nbolbl(6,inbo) = ir
            nbolbl(7,inbo) = iy
            nbolbl(8,inbo) = istar
            nbolbl(9,inbo) = iright
          end if
c
c  two-center labels:
c
        else if(nctr.eq.2) then
          lbl = nameat(iatno(label(ib,4)))
          call debyte(lbl,ibyte)
          nbolbl(1,inbo) = ibyte(1)
          nbolbl(2,inbo) = ibyte(2)
          call idigit(label(ib,4),istr,nd,maxd)
          if(nd.eq.1) then
            nbolbl(4,inbo) = istr(1)
          else
            nbolbl(3,inbo) = istr(1)
            nbolbl(4,inbo) = istr(2)
          end if
          nbolbl(5,inbo) = ihyp
          lbl = nameat(iatno(label(ib,5)))
          call debyte(lbl,ibyte)
          nbolbl(6,inbo) = ibyte(1)
          nbolbl(7,inbo) = ibyte(2)
          call idigit(label(ib,5),istr,nd,maxd)
          if(nd.eq.1) then
            nbolbl(9,inbo) = istr(1)
          else
            nbolbl(8,inbo) = istr(1)
            nbolbl(9,inbo) = istr(2)
          end if
          nbolbl(10,inbo) = label(ib,2)
c
c  three-center labels:
c
        else
          call idigit(label(ib,4),istr,nd,maxd)
          if(nd.eq.1) then
            nbolbl(2,inbo) = istr(1)
          else
            nbolbl(1,inbo) = istr(1)
            nbolbl(2,inbo) = istr(2)
          end if
          nbolbl(3,inbo) = ihyp
          call idigit(label(ib,5),istr,nd,maxd)
          if(nd.eq.1) then
            nbolbl(5,inbo) = istr(1)
          else
            nbolbl(4,inbo) = istr(1)
            nbolbl(5,inbo) = istr(2)
          end if
          nbolbl(6,inbo) = ihyp
          call idigit(label(ib,6),istr,nd,maxd)
          if(nd.eq.1) then
            nbolbl(8,inbo) = istr(1)
          else
            nbolbl(7,inbo) = istr(1)
            nbolbl(8,inbo) = istr(2)
          end if
          nbolbl(9,inbo) = label(ib,2)
        end if
   20 continue
      return
      end
c*****************************************************************************
      subroutine lblnho(inho,inbo,ictr,nctr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      parameter(maxd = 2)
      integer istr(maxd),ibyte(4)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),lbl1(maxbas),
     +       lorbc(maxbas),lorb(maxbas)
      common/nblbl/nlew,nval,iaolbl(10,maxbas),naolbl(10,maxbas),
     +                       nholbl(10,maxbas),nbolbl(10,maxbas)
c
      data iblnk,ic,il,ip,ir,iy,i3,istar,ihyp/1h ,1hc,1hl,1hp,1hr,1hy,
     +                                        1h3,1h*,1h-/
      data icr,ilp/2hcr,2hlp/
      data ileft,iright/1h(,1h)/
c
      do 10 i = 1,10
        nholbl(i,inho) = iblnk
   10 continue
      ib = ibxm(inbo)
c
c  one-center labels:
c
      if(nctr.eq.1) then
        lbl = nameat(iatno(label(ib,4)))
        call debyte(lbl,ibyte)
        nholbl(1,inho) = ibyte(1)
        nholbl(2,inho) = ibyte(2)
        call idigit(label(ib,4),istr,nd,maxd)
        if(nd.eq.1) then
          nholbl(4,inho) = istr(1)
        else
          nholbl(3,inho) = istr(1)
          nholbl(4,inho) = istr(2)
        end if
        nholbl(5,inho) = ileft
        if(label(ib,1).eq.icr) then
          nholbl(6,inho) = ic
          nholbl(7,inho) = ir
          nholbl(8,inho) = iright
        else if(label(ib,1).eq.ilp) then
          nholbl(6,inho) = il
          nholbl(7,inho) = ip
          if(label(ib,2).eq.istar) then
            nholbl(8,inho) = istar
            nholbl(9,inho) = iright
          else
            nholbl(8,inho) = iright
          end if
        else
          nholbl(6,inho) = ir
          nholbl(7,inho) = iy
          nholbl(8,inho) = istar
          nholbl(9,inho) = iright
        end if
c
c  two-center and three-center labels:
c
      else
        lbl = nameat(iatno(label(ib,3+ictr)))
        call debyte(lbl,ibyte)
        nholbl(1,inho) = ibyte(1)
        nholbl(2,inho) = ibyte(2)
        call idigit(label(ib,3+ictr),istr,nd,maxd)
        if(nd.eq.1) then
          nholbl(4,inho) = istr(1)
        else
          nholbl(3,inho) = istr(1)
          nholbl(4,inho) = istr(2)
        end if
        nholbl(5,inho) = ileft
        if(nctr.eq.2) then
          lbl = nameat(iatno(label(ib,6-ictr)))
          call debyte(lbl,ibyte)
          nholbl(6,inho) = ibyte(1)
          nholbl(7,inho) = ibyte(2)
          call idigit(label(ib,6-ictr),istr,nd,maxd)
          if(nd.eq.1) then
            nholbl(9,inho) = istr(1)
          else
            nholbl(8,inho) = istr(1)
            nholbl(9,inho) = istr(2)
          end if
          nholbl(10,inho) = iright
        else
          nholbl(6,inho) = i3
          nholbl(7,inho) = ihyp
          nholbl(8,inho) = ic
          nholbl(9,inho) = iright
        end if
      end if
      return
      end
c*****************************************************************************
c
c  general utility routines:
c
c      subroutine angles(x,y,z,theta,phi)
c      function bdfind(iat,jat)
c      subroutine chem(nat,natoms,lista,nl,istr)
c      subroutine consol(aut,alt,ndim,n)
c      subroutine convin(ij,len,ik,error)
c      subroutine convrt(n,nc1,nc2)
c      subroutine copy(a,b,ndim,nr,nc)
c      subroutine cortbl(iat,icore,iecp)
c      subroutine debyte(i,ibyte)
c      subroutine halt(word)
c      subroutine idigit(kint,ik,nd,maxd)
c      function ihtyp(ibo,jbo)
c      subroutine NBOjacobi(n,a,eivu,eivr,ndim,nvdim,icontr)
c      subroutine limtrn(t,m,a,b,ndim,nbas,ncdim,nc,iopt)
c      subroutine matmlt(a,b,v,ndim,n)
c      subroutine matml2(a,b,v,ndim,n)
c      function nameat(iz)
c      subroutine normlz(a,s,m,n)
c      subroutine NBOorder(rank,list,n,ndim,arcrnk)
c      subroutine NBOpack(t,ndim,nbas,l2)
c      subroutine rank(eig,n,ndim,arcrnk)
c      subroutine simtrn(a,t,v,ndim,n)
c      subroutine simtrs(a,s,v,ndim,n)
c      subroutine NBOtransp(a,ndim,n)
c      subroutine NBOunpack(t,ndim,nbas,l2)
c      subroutine valtbl(iat,ival)
c      function veclen(x,n,ndim)
c
c      subroutine lineq(a,x,b,scr,n,m,ndim,mdim,zertol,eps,maxit,lfnpr,
c     +                 ierr)
c      subroutine factor(a,w,d,ipivot,n,ndim,zertol,iflag)
c      subroutine fndsol(a,x,b,w,r,e,ipivot,n,ndim,eps,maxit,lfnpr,ierr)
c      subroutine subst(x,w,b,ipivot,n,ndim)
c
c*****************************************************************************
      subroutine angles(x,y,z,theta,phi)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      data zero,cutoff,one/0.0d0,1.0d-8,1.0d0/
c
      conv = 180.0d0/(4.0d0*atan(one))
      if(x.eq.zero.and.y.eq.zero) then
        if(z.ge.zero) then
          theta = zero
        else
          theta = 180.0d0
        end if
        phi = zero
      else
        if(abs(z-one).lt.cutoff) then
          theta = zero
        else if(abs(z+one).lt.cutoff) then
          theta = 180.0d0
        else if(z.lt.one.and.z.gt.-one) then
          theta = acos(z) * conv
          if(theta.gt.180.0d0) theta = 360.0d0 - theta
        else
          stop 'arccosine out of bounds in sr angles'
        end if
        phi   = atan2(y,x) * conv
        if(phi.lt.zero) phi = phi + 360.0d0
        if(abs(phi-360.0d0).lt.0.05d0) phi = zero
        if(abs(theta).lt.0.05d0) phi = zero
        if(abs(theta-180.0d0).lt.0.05d0) phi = zero
      end if
      return
      end
c*****************************************************************************
      function bdfind(iat,jat)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical bdfind,ifound,jfound
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbbas/label(maxbas,6),naoctr(maxbas),naol(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
c
      data lstar/1h*/
c
c  set bdfind=.true. if there is at least one bond between atoms iat and jat
c
      do 100 ibas = 1,nbas
        ib = ibxm(ibas)
        if(label(ib,2).eq.lstar) go to 100
        if(label(ib,3).ne.1) go to 100
        ifound = .false.
        jfound = .false.
        do 50 k = 4,6
          if(label(ib,k).eq.iat) ifound = .true.
          if(label(ib,k).eq.jat) jfound = .true.
   50   continue
        if(ifound.and.jfound) go to 200
  100 continue
      bdfind = .false.
      return
  200 bdfind = .true.
      return
      end
c*****************************************************************************
      subroutine chem(nat,natoms,lista,nl,istr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension lista(natoms,2),istr(80)
c
      parameter (maxd = 4)
      dimension inum(maxd),ibyte(4)
      data ic,ih,iblnk,ileft,iright/1hc,1hh,1h ,1h(,1h)/
c
c  build the chemical formula from the list of atoms in lista:
c
c  get chemical symbols:
c
      do 10 iat = 1,nat
        lista(iat,1) = nameat(lista(iat,1))
   10 continue
c
c  alphabetize these symbols:
c
      do 30 iat = 1,nat-1
        do 20 jat = 1,nat-iat
          if(lista(jat,1).gt.lista(jat+1,1)) then
            itemp = lista(jat,1)
            lista(jat,1) = lista(jat+1,1)
            lista(jat+1,1) = itemp
            itemp = lista(jat,2)
            lista(jat,2) = lista(jat+1,2)
            lista(jat+1,2) = itemp
          end if
   20   continue
   30 continue
c
c  build chemical formula in istr:
c
c  first carbon...
c
      nl = 1
      istr(nl) = ileft
      do 50 iat = 1,nat
        call debyte(lista(iat,1),ibyte)
        if(ibyte(1).eq.iblnk.and.ibyte(2).eq.ic) then
          nl = nl + 1
          istr(nl) = ic
          if(lista(iat,2).ne.1) then
            call idigit(lista(iat,2),inum,nd,maxd)
            do 40 il = 1,nd
              nl = nl + 1
              istr(nl) = inum(il)
   40       continue
          end if
          lista(iat,2) = 0
        end if
   50 continue
c
c  then hydrogen...
c
      do 70 iat = 1,nat
        call debyte(lista(iat,1),ibyte)
        if(ibyte(1).eq.iblnk.and.ibyte(2).eq.ih) then
          nl = nl + 1
          istr(nl) = ih
          if(lista(iat,2).ne.1) then
            call idigit(lista(iat,2),inum,nd,maxd)
            do 60 il = 1,nd
              nl = nl + 1
              istr(nl) = inum(il)
   60       continue
          end if
          lista(iat,2) = 0
        end if
   70 continue
c
c  and now the rest...
c
      do 90 iat = 1,nat
        if(lista(iat,2).ne.0) then
          call debyte(lista(iat,1),ibyte)
          if(ibyte(1).ne.iblnk) then
            nl = nl + 1
            istr(nl) = ibyte(1)
          end if
          if(ibyte(2).ne.iblnk) then
            nl = nl + 1
            istr(nl) = ibyte(2)
          end if
          if(lista(iat,2).ne.1) then
            call idigit(lista(iat,2),inum,nd,maxd)
            do 80 il = 1,nd
              nl = nl + 1
              istr(nl) = inum(il)
   80       continue
          end if
          lista(iat,2) = 0
        end if
   90 continue
      nl = nl + 1
      istr(nl) = iright
      return
      end
c*****************************************************************************
      subroutine consol(aut,alt,ndim,n)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  consolidate aut, alt to a single matrix, with aut as upper triangle
c  (including diagonal) and alt as lower triangle.  store result in aut.
c
      dimension aut(ndim,ndim),alt(ndim,ndim)
      nm1=n-1
      do 10 j=1,nm1
        jp1=j+1
        do 10 i=jp1,n
   10     aut(i,j)=alt(i,j)
      return
      end
c*****************************************************************************
      subroutine convin(ij,len,ik,error)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension ij(1)
      dimension int(10)
      logical error
c
      data int/1h0,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9/
c
c  convert the array ij(len) into an integer ik:
c
      error = .false.
      if(len.le.0) then
        error = .true.
        return
      end if
c
c  make sure all elements of ij are integers:
c
      il   = 0
      mult = 1
      do 30 i = len,1,-1
        do 10 j = 1,10
          jj = j - 1
          if(ij(i).eq.int(j)) goto 20
   10   continue
        error = .true.
        return
c
   20   il = il + jj * mult
        mult = mult * 10
   30 continue
      ik = il
      return
      end
c*****************************************************************************
      subroutine convrt(n,nc1,nc2)
c*****************************************************************************
c
c  convert 2-digit integer 'n' to two literal characters 'nc1','nc2'.
c
      dimension int(10)
      data isp,int/1h ,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1h0/
c
      nc1=isp
      nc2=isp
      if(n.le.0) return
       if(n.ge.10) go to 10
        nc2=int(n)
        return
   10 n1=n/10
       if(n1.gt.9) stop 'routine convrt'
       nc1=int(n1)
       n2=n-n1*10
       if(n2.eq.0) n2=10
       nc2=int(n2)
       return
      end
c*****************************************************************************
      subroutine copy(a,b,ndim,nr,nc)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension a(ndim,1),b(ndim,1)
c
c  copy a to b:
c
      do 20 j = 1,nc
        do 10 i = 1,nr
          b(i,j) = a(i,j)
   10   continue
   20 continue
      return
      end
c*****************************************************************************
      subroutine cortbl(iat,icore,iecp)
c*****************************************************************************
c
c   core table:
c
c     determine the number of subshells of core orbitals of each angular
c     symmetry for atom number iat.  icore is an integer array lmax+1
c     long which returns the number of subshells to the calling subroutine:
c     the number of s subshells in icore(1), the number of p subshells
c     in icore(2), etc...
c
c     if the core option has been used, the core orbitals stored in the array
c     iatcr are used rather than the core orbitals of the nominal core table.
c
c     if iecp = 0 return the number of subshells, excluding subshells of
c                 an effective core potential.
c     if iecp = 1 return the number of subshells, including subshells of
c                 an effective core potential.
c
c     note: it is possible for a negative number of core orbitals be found
c     if effective core potentials are employed.  this happens when the
c     number of core electrons in the effective core potential is either
c     greater than the nominal number of core electrons or is greater than the
c     number of core electrons requested when using the core option.
c
c------------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (lmax = 3)
      integer core(57),icore(4),itemp(4),iord(16),jord(20),kord(20)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
c
      data iord/1,1,3,1,3,5,1,3,5,1,3,7,5,1,3,7/
      data jord/1,1,3,1,3,1,5,3,1,5,3,1,7,5,3,1,7,5,3,1/
      data kord/1,2,1,3,2,4,1,3,5,2,4,6,1,3,5,7,2,4,6,8/
      data core/2,0,8,1,1,8,2,2,1,12,2,3,2,6,3,3,2,1,12,3,4,3,1,6,3,4,3,
     +  2,16,3,5,4,2,10,4,5,4,2,1,6,4,5,4,3,1,16,4,6,5,3,1,10,4,6,5,3,2/
c
c  initialize arrays.  if there is no nuclear charge at this center,
c  return to calling routine:
c
      do 10 l = 0,lmax
        icore(l+1) = 0
        itemp(l+1) = 0
   10 continue
      if(iatno(iat).le.0) return
c
c  if the core option has not been used for this atom, use the nominal
c  set of core orbitals:
c
      if(jcore.ne.1.or.iatcr(iat).lt.0) then
        jat = iatno(iat)
        ii  = 0
   20   ii  = ii + 1
          jat = jat - core(ii)
          ii = ii + 1
          if(jat.le.0) then
            do 30 l = 1,core(ii)
              icore(l) = core(ii+l)
   30       continue
          else
            ii = ii + core(ii)
          end if
        if(jat.gt.0) goto 20
      else
c
c  if the core option has been used, determine the number of core
c  orbitals on this atom:
c
        ii = iatcr(iat)
        if(ii.gt.0) then
          ict = 0
   40     ict = ict + 1
          l = iord(ict)/2
          icore(l+1) = icore(l+1) + 1
          ii = ii - iord(ict)
          if(ii.gt.0) goto 40
        end if
      end if
c
c  if effective core potentials were used and iecp = 0, remove
c  the core orbitals of the ecp:
c
      if(ipseud.ne.0.and.iecp.eq.0) then
        ii = iatno(iat)
        ict = 0
   50   ict = ict + 1
          ii = ii - 2 * jord(ict)
        if(ii.gt.0) goto 50
        ii = iznuc(iat) - ii
        if(ii.le.0) stop 'zero or negative iznuc entry?'
        ict = ict + 1
   60   ict = ict - 1
          if(ict.le.0) stop 'error in sr cortbl'
          ii = ii - 2 * jord(ict)
          if(ii.ge.0) then
            l = jord(ict)/2
            if(icore(l+1).ge.kord(ict)) itemp(l+1) = itemp(l+1) + 1
          else
            ii = ii + 2 * jord(ict)
          end if
        if(ii.ne.0) goto 60
        do 70 l = 0,lmax
          icore(l+1) = itemp(l+1)
   70   continue
      end if
      return
      end
c*****************************************************************************
      subroutine debyte(i,ibyte)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension ibyte(4),kb(4)
c
      save kb,kpad,ksw
c
      data ksw/0/
      data ktmp/4habcd/
c
c  extract four hollerith characters from i, store in ibtye:
c
c  if this is the first time that this routine is called, determine
c  in which bytes of an integer word the hollerith characters reside:
c
      if(ksw.eq.0) then
        ksw   = 1
        do 10 k = 1,4
          kb(k) = 0
   10   continue
        kbyte = 0
   20   kbyte = kbyte + 1
          if(kbyte.gt.8) stop 'routine debyte is limited to integer*8'
          ktest = mod(ktmp,256)
c
c PS - changed  to map onto a-c ASCII
c
          if(ktest.eq.97) kb(1) = kbyte
          if(ktest.eq.98) kb(2) = kbyte
          if(ktest.eq.99) kb(3) = kbyte
          if(ktest.eq.100) kb(4) = kbyte
          ktmp = ktmp/256
        if(ktmp.ne.0) goto 20
        do 30 k = 1,4
          if(kb(k).eq.0) stop 'error in routine debyte'
   30   continue
c
c  determine the bit padding:
c
        kpad = 0
        kmlt = 1
        do 40 k = 1,kbyte
          if(k.ne.kb(1)) kpad = kpad + 32 * kmlt
          if(k.ne.kbyte) kmlt = kmlt * 256
   40   continue
c
        do 60 k = 1,4
          kmax  = kb(k) - 1
          kb(k) = 1
          do 50 l = 1,kmax
            kb(k) = kb(k) * 256
   50     continue
   60   continue
      end if
c
c  extract four hollerith characters from i:
c
      do 100 k = 1,4
        ibyte(k) = mod(i/kb(k),256)*kb(1) + kpad
  100 continue
      return
      end
c*****************************************************************************
      subroutine halt(word)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      data blank/1h /
c
      if(word.eq.blank) return
      write(lfnpr,1000) word
      stop
c
 1000 format(' non-integer encountered when trying to read variable ',
     + '/',a6,'/')
      end
c*****************************************************************************
      subroutine idigit(kint,ik,nd,maxd)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension ik(maxd),int(10)
c
      data iblnk,int/1h ,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1h0/
c
c  converts the integer kint into the first nd elements of hollerith array
c  ik(maxd):
c
      jint = kint
      nd   = maxd
      do 10 id = maxd,1,-1
        ii = mod(jint,10)
        if(ii.eq.0) ii = 10
        ik(id) = int(ii)
        if(ii.ne.10) nd = id
        jint = jint/10
   10 continue
      nd = maxd - nd + 1
c
c  shift integer rep in ik so that the number occupies the first nd
c  elements:
c
      do 20 id = 1,nd
        ik(id) = ik(id+maxd-nd)
   20 continue
      do 30 id = nd+1,maxd
        ik(id) = iblnk
   30 continue
      return
      end
c*****************************************************************************
      function ihtyp(ibo,jbo)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      logical bdfind
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbbas/label(maxbas,6),nbouni(maxbas),nbotyp(maxbas),
     +       lstocc(maxbas),ibxm(maxbas),larc(maxbas),iathy(maxbas,3)
c
      data iv,ig,ir/1hv,1hg,1hr/
c
c  determine whether the ibo->jbo delocalization is vicinal (ihtyp='v'),
c  geminal (ihtyp='g'), or remote (ihtyp='r'):
c
      ihtyp = ir
      if(nbouni(ibo).eq.nbouni(jbo)) then
        ictr = mod(nbotyp(ibo),10)
        ib   = ibxm(ibo)
        jctr = mod(nbotyp(jbo),10)
        jb   = ibxm(jbo)
        do 20 i = 1,ictr
          iat = label(ib,i+3)
          do 10 j = 1,jctr
            jat = label(jb,j+3)
            if(iat.eq.jat) then
              ihtyp = ig
              return
            else if(bdfind(iat,jat)) then
              ihtyp = iv
            end if
   10     continue
   20   continue
      end if
      return
      end
c*****************************************************************************
      subroutine NBOjacobi(n,a,eivu,eivr,ndim,nvdim,icontr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  diagonalize real symmetric matrix a by jacobi rotations:
c        n: actual dimension of a,eivr
c     ndim: declared dimension of a,eivr
c   icontr: control option
c
c
c       ********  modified version, march 1986  *************
c
c
c     icontr = 0: reduce all off-diagonal elements to "done" or smaller
c                                       -- this sets fulmix=.true.
c
c     icontr = 1: do the same as for icontr=0 except do not mix orbitals that
c       are degenerate to within "differ" if the offdiagonal element connecting
c       them is less than "differ".
c                                       -- this sets fulmix=.false.
c
c  for the purposes of the nao and nbo programs, these values are set:
c
c     differ = 1.0d-5
c          threshold for considering two vectors nondegenerate if
c                         icontr=1
c     done   = 1.0d-13
c          this is the parameter for convergence of the off-diagonal
c             matrix elements.  (absolute)  --- reduced from 1.0d-10
c             on 8/31/88.  a more converged fock matrix was required
c             for the nbo deletions with symmetry to work properly
c             (edg) ---
c
c     eps    = 0.5d-13
c          this parameter has to do with the machine precision and should
c             be set to a value between "done" and the machine precision.
c             --- reduced from 1.0d-11.  8/31/88 (edg) ---
c
      logical fulmix
      dimension a(ndim,1),eivr(nvdim,1),eivu(1)
c  important parameters:
      data differ,done,eps,pt99/1.0d-5,1.0d-13,0.5d-13,0.99d0/
      data zero,one,five/0.0d0,1.0d0,5.0d0/
c
      fulmix=.true.
      if(icontr.eq.1) fulmix=.false.
      if(n.gt.1) go to 10
       eivr(1,1)=one
       eivu(1)=a(1,1)
       return
   10 continue
      do 30 j=1,n
        do 20 i=1,n
   20     eivr(i,j)=zero
   30   eivr(j,j)=one
c
c        find the absolutely largest element of a
c
c  first check the off-diagonal elements:
      atop=zero
      do 50 j=2,n
        jm1=j-1
        do 50 i=1,jm1
          if(atop.gt.abs(a(i,j))) go to 50
          atop= abs(a(i,j))
   50     continue
      offtop=atop
c  now check the diagonal elements:
      do 60 j=1,n
          if(atop.gt.abs(a(j,j))) go to 60
          atop= abs(a(j,j))
   60     continue
c  if matrix is already effectively diagonal,
c              put diagonal elements in eivu and return
      if(atop.lt.done) go to 260
      if(offtop.lt.done) go to 260
c
c        calculate the stopping criterion -- dstop
c
      avgf= dble(n*(n-1)/2)
      d=0.0d0
      do 80 jj=2,n
        do 80 ii=2,jj
          s=a(ii-1,jj)/atop
   80     d=s*s+d
      dstop=(1.d-7)*d
c
c        calculate the threshold, thrsh
c
      thrsh= sqrt(d/avgf)*atop
c  to make thrsh different than any matrix element of a, multiply by 0.99
      thrsh=thrsh*pt99
      if(thrsh.lt.done) thrsh=done
c
c        start a sweep
c
   90 iflag=0
      do 250 jcol=2,n
        jcol1=jcol-1
        do 250 irow=1,jcol1
          aij=a(irow,jcol)
c
c        compare the off-diagonal element with thrsh
c
          absaij=abs(aij)
          if (absaij.lt.thrsh) go to 250
          aii=a(irow,irow)
          ajj=a(jcol,jcol)
          s=ajj-aii
          abss=abs(s)
c  don't rotate the vectors irow and jcol if irow and jcol would still
c     be degenerate within "differ":
          if(fulmix) go to 100
            if((abss.lt.differ).and.(absaij.lt.differ)) go to 250
  100     continue
c
c        check to see if the chosen rotation is less than the rounding error
c        if so , then do not rotate.
c
          test=eps*abss
          if (absaij.lt.test) go to 250
          iflag=1
c
c        if the rotation is very close to 45 degrees, set sin and cos
c        to 1/(root 2).
c
          test=eps*absaij
          if (abss.gt.test) go to 130
          s=.707106781d0
          c=s
          go to 140
c
c        calculation of sin and cos for rotation that is not very close
c        to 45 degrees
c
  130     t=aij/s
          s=0.25d0/ sqrt(0.25d0+t*t)
c
c        cos=c ,  sin=s
c
          c= sqrt(0.5d0+s)
          s=2.d0*t*s/c
c
c        calculation of the new elements of matrix a
c
  140     do 150 i=1,irow
            t=a(i,irow)
            u=a(i,jcol)
            a(i,irow)=c*t-s*u
  150       a(i,jcol)=s*t+c*u
          i2=irow+2
          if (i2.gt.jcol) go to 180
          do 170 i=i2,jcol
            t=a(i-1,jcol)
            u=a(irow,i-1)
            a(i-1,jcol)=s*u+c*t
  170       a(irow,i-1)=c*u-s*t
  180     a(jcol,jcol)=s*aij+c*ajj
          a(irow,irow)=c*a(irow,irow)-s*(c*aij-s*ajj)
          do 190 j=jcol,n
            t=a(irow,j)
            u=a(jcol,j)
            a(irow,j)=c*t-s*u
  190       a(jcol,j)=s*t+c*u
c
c        rotation completed
c
          do 210 i=1,n
            t=eivr(i,irow)
            eivr(i,irow)=c*t-eivr(i,jcol)*s
  210       eivr(i,jcol)=s*t+eivr(i,jcol)*c
c
c        calculate the new norm d and compare with dstop
c
          s=aij/atop
          d=d-s*s
          if (d.gt.dstop) go to 240
c
c        recalculate dstop and thrsh to discard rounding errors
c
          d=zero
          do 230 jj=2,n
            do 230 ii=2,jj
              s=a(ii-1,jj)/atop
  230         d=s*s+d
          dstop=(1.d-7)*d
  240     continue
          oldthr=thrsh
          thrsh= sqrt(d/avgf)*atop*pt99
          if(thrsh.lt.done) thrsh=done*pt99
          if(thrsh.gt.oldthr) thrsh=oldthr
  250     continue
      if(thrsh.lt.done) go to 260
      if(iflag.eq.1) go to 90
      thrsh=thrsh/five
      go to 90
c
c        place eigenvalues in eivu
c
  260 continue
      do 270 j=1,n
        eivu(j)=a(j,j)
  270   continue
      return
      end
c*****************************************************************************
      subroutine limtrn(t,m,a,b,ndim,nbas,ncdim,nc,iopt)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension t(ndim,ndim),m(ncdim),a(ncdim,ncdim),b(ncdim)
c...do a limited transformation of t, including only the "nc" rows and
c    columns specified in the vector "m":
c
c   iopt= 1 :  take t=t*a
c   iopt= 0 :  take t=a(transpose)*t*a
c   iopt=-1 :  take t=a(transpose)*t
c
c
      if(iopt.eq.1) go to 100
c   first, take t=a(transpose)*t, where t=s,dm
        do 30 j=1,nbas
          do 10 k=1,nc
   10     b(k)=t(m(k),j)
          do 30 i=1,nc
            sum=0.0d0
            do 20 k=1,nc
   20         sum=sum+a(k,i)*b(k)
   30       t(m(i),j)=sum
      if(iopt.eq.-1) return
c   now, take t=t*a
  100 continue
        do 160 i=1,nbas
          do 140 k=1,nc
  140       b(k)=t(i,m(k))
          do 160 j=1,nc
            sum=0.0d0
            do 150 k=1,nc
  150         sum=sum+b(k)*a(k,j)
  160       t(i,m(j))=sum
      return
      end
c*****************************************************************************
      subroutine matmlt(a,b,v,ndim,n)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension a(1),b(1),v(ndim)
      data zero/0.0d0/
c
c  multiply a*b (using scratch vector v), store result in a:
c
      ndif=ndim-n
      do 30 i=1,n
        kj=0
        ikk=i-ndim
        do 20 j=1,n
          ik=ikk
          temp=zero
          do 10 k=1,n
            ik=ik+ndim
            kj=kj+1
   10       temp=temp+a(ik)*b(kj)
          kj=kj+ndif
   20   v(j)=temp
        ij=i-ndim
        do 30 j=1,n
          ij=ij+ndim
   30     a(ij)=v(j)
      return
      end
c*****************************************************************************
      subroutine matml2(a,b,v,ndim,n)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension a(1),b(1),v(ndim)
      data zero/0.0d0/
c                       b=a(transpose)*b
c  multiply a(transpose)*b (using scratch vector v), store result in b:
c    assume a*b is a symmetric matrix, so almost half the work is saved.
c    this can be the second step in a similarity transformation of b by a.
c
      ij=0
      ijj=-ndim
      kjj=-ndim
      do 50 j=1,n
        kii=-ndim
        kjj=kjj+ndim
        do 20 i=1,j
          kii=kii+ndim
          ki=kii
          kj=kjj
          temp=zero
          do 10 k=1,n
            ki=ki+1
            kj=kj+1
   10       temp=temp+a(ki)*b(kj)
   20   v(i)=temp
        ijj=ijj+ndim
        ij=ijj
        ji=j-ndim
        jm1=j-1
        do 30 i=1,jm1
          ij=ij+1
          ji=ji+ndim
          vv=v(i)
          b(ij)=vv
   30     b(ji)=vv
        ij=ij+1
   50   b(ij)=v(j)
      return
      end
c*****************************************************************************
      function nameat(iz)
c*****************************************************************************
c
c  return atomic symbol for nuclear charge iz (.le. 103):
c

      dimension name(103)
      data ighost/2hgh/iblank/2h  /
      data name/2h h,2hhe,2hli,2hbe,2h b,2h c,2h n,2h o,2h f,2hne,
     + 2hna,2hmg,2hal,2hsi,2h p,2h s,2hcl,2har,2h k,2hca,2hsc,2hti,
     + 2h v,2hcr,2hmn,2hfe,2hco,2hni,2hcu,2hzn,2hga,2hge,2has,
     + 2hse,2hbr,2hkr,2hrb,2hsr,2h y,2hzr,2hnb,2hmo,2htc,2hru,
     + 2hrh,2hpd,2hag,2hcd,2hin,2hsn,2hsb,2hte,2h i,2hxe,2hcs,
     + 2hba,2hla,2hce,2hpr,2hnd,2hpm,2hsm,2heu,2hgd,2htb,2hdy,
     + 2hho,2her,2htm,2hyb,2hlu,2hhf,2hta,2h w,2hre,2hos,2hir,
     + 2hpt,2hau,2hhg,2htl,2hpb,2hbi,2hpo,2hat,2hrn,2hfr,2hra,
     + 2hac,2hth,2hpa,2h u,2hnp,2hpu,2ham,2hcm,2hbk,2hcf,2hes,
     + 2hfm,2hmd,2hno,2hlr/
c
      if(iz.lt.0.or.iz.gt.103) nameat = iblank
      if(iz.gt.0) nameat = name(iz)
      if(iz.eq.0) nameat = ighost
      return
      end
c*****************************************************************************
      subroutine normlz(a,s,m,n)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension a(m,m),s(m,m)
c
      data zero,one /0.0d0,1.0d0/
c
c normalize columns of a
c
      do 40 i = 1,n
        temp = zero
        do 20 j = 1,n
          do 10 k = 1,n
            temp = temp + a(j,i)*a(k,i)*s(j,k)
   10     continue
   20   continue
        factor = one/sqrt(temp)
        do 30 j = 1,n
          a(j,i) = factor * a(j,i)
   30   continue
   40 continue
      return
      end
c*****************************************************************************
      subroutine NBOorder(rank,list,n,ndim,arcrnk)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  rank positive elements of integer 'list', lowest values first.
c
      integer rank,arcrnk,temp
      dimension rank(ndim),list(ndim),arcrnk(ndim)
      do 10 i=1,n
   10   arcrnk(i)=i
      do 40 i=1,n
        if(i.eq.n)go to 30
        i1=i+1
        do 20 j=i1,n
          if(list(j).ge.list(i))go to 20
          temp=list(i)
          list(i)=list(j)
          list(j)=temp
          temp=arcrnk(i)
          arcrnk(i)=arcrnk(j)
          arcrnk(j)=temp
   20     continue
   30   rank(arcrnk(i))=i
        if(list(i).le.0) go to 50
   40   continue
      return
   50 do 60 k=i,n
        rank(arcrnk(k))=0
   60   continue
      return
      end
c*****************************************************************************
      subroutine NBOpack(t,ndim,nbas,l2)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension t(1)
c
      data zero/0.0d0/
c
c  pack:  packs a symmetric matrix t into an upper triangular matrix.
c         t should be dimensioned (ndim,ndim) in the calling routine:
c
      if(nbas.gt.ndim) stop 'nbas is greater than ndim'
      ii = 0
      do 200 j = 1,nbas
        jptr = (j-1) * ndim
        do 100 i = 1,j
          iptr = jptr + i
          ii = ii + 1
          t(ii) = t(iptr)
  100   continue
  200 continue
      if(ii.ne.l2) stop 'error in routine NBOpack'
c
      do 300 i = ii+1,ndim*ndim
        t(i) = zero
  300 continue
      return
      end
c*****************************************************************************
      subroutine rank(eig,n,ndim,arcrnk)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  order numbers in 'eig', highest values first,
c    and construct 'arcrnk':
c     arcrnk(i) is the old location of the i-th highest value in eig
c      note: upon return, eig(i) is the i-th highest value in eig
c      important: numbers in eig are not switched unless they differ
c       by more than "differ":  5.0d-8
c
      integer arcrnk
      dimension arcrnk(ndim),eig(ndim)
      data differ/5.0d-8/
      do 10 i=1,n
   10   arcrnk(i)=i
      do 40 i=1,n
        if(i.eq.n)go to 40
        i1=i+1
        do 20 j=i1,n
          if((eig(j)-eig(i)).lt.differ) go to 20
          temp=eig(i)
          eig(i)=eig(j)
          eig(j)=temp
          itemp=arcrnk(i)
          arcrnk(i)=arcrnk(j)
          arcrnk(j)=itemp
   20     continue
   40   continue
      return
      end
c*****************************************************************************
      subroutine simtrn(a,t,v,ndim,n)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  similarity transform a ==> t(transpose)*a*t, using scratch vector v.
c
      dimension a(ndim,ndim),t(ndim,ndim),v(ndim)
      call matmlt(a,t,v,ndim,n)
      call NBOtransp(a,ndim,n)
      call matmlt(a,t,v,ndim,n)
      call NBOtransp(a,ndim,n)
      return
      end
c*****************************************************************************
      subroutine simtrs(a,s,v,ndim,n)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
c
c  similarity transform a ==> s(transpose)*a*s, using scratch vector v.
c    fast version --- assumes result is a symmetric matrix
c
      dimension a(ndim,ndim),s(ndim,ndim),v(ndim)
      call matmlt(a,s,v,ndim,n)
      call matml2(s,a,v,ndim,n)
      return
      end
c*****************************************************************************
      subroutine NBOtransp(a,ndim,n)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension a(ndim,ndim)
c
c  transpose matrix a, store result in a.
c
      do 10 i=1,n
        do 10 j=1,i
          temp=a(i,j)
          a(i,j)=a(j,i)
   10     a(j,i)=temp
      return
      end
c*****************************************************************************
      subroutine NBOunpack(t,ndim,nbas,l2)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension t(1)
c
c  unpack:  unpacks an upper triangular matrix (vector l2 long) into a
c           symmetric matrix t(nbas,nbas).  note: t should be dimensioned
c           (ndim,ndim) in the calling routine.
c
c  first spread out the l2 numbers into the upper part of the whole array.
c
      j = 0
      k = 1
      iptr = (ndim + 1)*(nbas - k) + 1
      do 200 i = l2,1,-1
        t(iptr-j) = t(i)
        if(j.lt.nbas-k) then
          j = j + 1
        else
          j = 0
          k = k + 1
          iptr = (ndim + 1)*(nbas - k) + 1
        end if
  200 continue
c
c  now fill in the holes in the output array.
c
      do 400 j = 1,nbas-1
        icol = (j-1)*ndim
        do 300 i = j+1,nbas
          iptr = icol + i
          jptr = (i-1)*ndim + j
          t(iptr) = t(jptr)
  300   continue
  400 continue
c
      return
      end
c*****************************************************************************
      subroutine valtbl(iat,ival)
c*****************************************************************************
c
c   valence table:
c
c     determine the number of sets of valence orbitals of each angular
c     symmetry for atom number iat.  ival is an integer array lmax+1
c     long which returns the number of sets to the calling subroutine:
c     the number of s subshells in ival(1), the number of p subshells
c     in ival(2), etc...
c
c------------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (lmax = 3)
      dimension ival(4),icore(4),iord(20)
c
      parameter (maxatm = 750, maxbas = 4096)
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseud,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)
c
      data iord/1,1,3,1,3,1,5,3,1,5,3,1,7,5,3,1,7,5,3,1/
c
      do 10 l = 0,lmax
        ival(l+1) = 0
   10 continue
c
c  count the number of filled or partially filled subshells:
c
      ii = iatno(iat)
      if(ii.gt.0) then
        ict = 0
   20   ict = ict + 1
        l = iord(ict)/2
        ival(l+1) = ival(l+1) + 1
        ii = ii - 2*iord(ict)
        if(ii.gt.0) goto 20
      end if
c
c  remove the core subshells.  note: if there are more core orbitals
c  in the effective core potential than in the nominal core table or
c  from the core option, remove these extra core orbitals from the
c  set of valence orbitals:
c
      iecp = 1
      call cortbl(iat,icore,iecp)
      do 50 l = 0,lmax
        ival(l+1) = ival(l+1) - icore(l+1)
   50 continue
      iecp = 0
      call cortbl(iat,icore,iecp)
      do 60 l = 0,lmax
        if(icore(l+1).lt.0) then
          ival(l+1) = ival(l+1) + icore(l+1)
        end if
   60 continue
      return
      end
c*****************************************************************************
      function veclen(x,n,ndim)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension x(ndim)
      data zero/0.0d0/
c
      sum = zero
      do 10 i = 1,n
        sum = sum + x(i)*x(i)
   10 continue
      veclen = sqrt(sum)
      return
      end
c*****************************************************************************
      subroutine lineq(a,x,b,scr,n,m,ndim,mdim,zertol,eps,maxit,lfnpr,
     +                 ierr)
c*****************************************************************************
c
c  solve the system of linear equations  a * x  =  b  for matrix x
c                                        ~   ~     ~             ~
c  input
c -------
c  * coefficient matrix a of dimension (n,n) with actual
c    dimension (ndim,ndim).
c  * matrix b of dimension (n,m) with actual dimension
c    (ndim,mdim)
c  * working space scr dimensioned (ndim,ndim+5).
c  * zero tolerance zertol.
c  * threshold on euclidean norm (vector length) of the
c    error vector relative to the norm of a column of x.
c  * maximum number of iterations maxit allowed during
c    iterative improvement.
c  * logical file number lfnpr for printing during iterative
c    improvement.  set to zero to no printing is desired.
c
c  output
c --------
c  * solution x of dimension (n,m) with actual dimension
c    (ndim,mdim).
c  * euclidean norm of the final error vector, eps.
c  * number of iterations taken during interative improvement,
c    maxit.
c  * error flag :    ierr = -1   iterative improvement did not
c                                converge
c                    ierr =  0   no errors encountered
c                    ierr =  1   a matrix is not invertible
c
c------------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(ndim,ndim),x(ndim,mdim),b(ndim,mdim),
     +          scr(ndim*(ndim+5))
      data zero/0.0/
c
      if(n.lt.1) stop 'dimension n is not positive'
c
c  partition scratch space:
c
      i1 = 1
      i2 = i1 + ndim*ndim
      i3 = i2 + ndim
      i4 = i3 + ndim
      i5 = i4 + ndim
      i6 = i5 + ndim
c
c  perform gauss elimination with scaled partial pivoting:
c
      call factor(a,scr(i1),scr(i2),scr(i6),n,ndim,zertol,iflag)
      if(iflag.eq.0) then
        ierr = 1
        return
      else
        ierr = 0
      end if
c
c  loop over columns of x and b:
c
      epsmax = zero
      itsmax = 0
      do 30 kcol = 1,m
        do 10 jrow = 1,n
          scr(i4+jrow-1) = x(jrow,kcol)
          scr(i5+jrow-1) = b(jrow,kcol)
   10   continue
        its = maxit
        del = eps
c
c  use back-substitution and iterative improvement to determine
c  the solution x:
c
        call fndsol(a,scr(i4),scr(i5),scr(i1),scr(i2),scr(i3),scr(i6),
     +              n,ndim,del,its,lfnpr,ierr)
        if(ierr.ne.0) return
c
c  copy solution into x:
c
        do 20 jrow = 1,n
          x(jrow,kcol) = scr(i4+jrow-1)
   20   continue
        if(del.gt.epsmax) epsmax = del
        if(its.gt.itsmax) itsmax = its
   30 continue
c
      eps = epsmax
      maxit = itsmax
      return
      end
c*****************************************************************************
      subroutine factor(a,w,d,ipivot,n,ndim,zertol,iflag)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension a(ndim,ndim),w(ndim,ndim),d(ndim),ipivot(ndim)
      data zero,one/0.0d0,1.0d0/
c
c  initial iflag.  if iflag is 1, then an even number of interchanges
c  has been carried out.  if equal to -1, then an odd number of inter-
c  changes have taken place.  if iflag is set to zero on return to the
c  calling routine, then the matrix is not invertible:
c
      iflag = 1
c
c  copy coefficient matrix a to w:
c
      call copy(a,w,ndim,n,n)
c
c  initialize d and ipivot:
c
      do 20 i = 1,n
        ipivot(i) = i
        rowmax = zero
        do 10 j = 1,n
          if(abs(w(i,j)).gt.rowmax) rowmax = abs(w(i,j))
   10   continue
        if(rowmax.le.zertol) then
          iflag = 0
          rowmax = one
        end if
        d(i) = rowmax
   20 continue
      if(n.eq.1) return
c
c  loop over rows, factorizing matrix w:
c
      do 100 k = 1,n-1
c
c  determine the pivot row istar:
c
        colmax = abs(w(k,k))/d(k)
        istar = k
        do 30 i = k+1,n
          temp = abs(w(i,k))/d(k)
          if(temp.gt.colmax) then
            colmax = temp
            istar = i
          end if
   30   continue
        if(colmax.eq.zero) then
          iflag = 0
        else
          if(istar.gt.k) then
            iflag = -iflag
            itemp = ipivot(istar)
            ipivot(istar) = ipivot(k)
            ipivot(k) = itemp
            temp = d(istar)
            d(istar) = d(k)
            d(k) = temp
            do 40 j = 1,n
              temp = w(istar,j)
              w(istar,j) = w(k,j)
              w(k,j) = temp
   40       continue
          end if
c
c  eliminate x(k) from rows k+1,...,n:
c
          do 60 i = k+1,n
            w(i,k) = w(i,k)/w(k,k)
            do 50 j = k+1,n
              w(i,j) = w(i,j) - w(i,k)*w(k,j)
   50       continue
   60     continue
        end if
  100 continue
      if(abs(w(n,n)).le.zertol) iflag = 0
      return
      end
c*****************************************************************************
      subroutine fndsol(a,x,b,w,r,e,ipivot,n,ndim,eps,maxit,lfnpr,ierr)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension a(ndim,ndim),x(ndim),b(ndim),w(ndim,ndim),r(ndim),
     +          e(ndim),ipivot(ndim)
      data zero/0.0d0/
c
c  find initial guess for x by back substitution:
c
      call copy(b,e,ndim,n,1)
      call subst(x,w,e,ipivot,n,ndim)
      if(maxit.eq.0) return
c
c  iterate until the vector length of the error vector relative to
c  x is less than eps:
c
      rellen = zero
      iter = 0
   10 if(rellen.gt.eps) then
        iter = iter + 1
        do 30 i = 1,n
          r(i) = b(i)
          do 20 j = 1,n
            r(i) = r(i) - a(i,j)*x(j)
   20     continue
   30   continue
        call subst(e,w,r,ipivot,n,ndim)
        elen = veclen(e,n,ndim)
        xlen = veclen(x,n,ndim)
        rellen = elen/xlen
        do 40 i = 1,n
          x(i) = x(i) + e(i)
   40   continue
c
c  print out iterative improvement info:
c
        if(lfnpr.gt.0) then
          write(lfnpr,900) iter,rellen
        end if
c
c  if too many iterations have taken place, halt furthur iterations:
c
        if(iter.eq.maxit) then
          if(rellen.gt.eps) ierr = -1
          if(lfnpr.gt.0) then
            if(ierr.lt.0) then
              write(lfnpr,910)
            else
              write(lfnpr,920)
            end if
          end if
          eps = rellen
          return
        end if
c
c  error vector is converged:
c
      else
        if(lfnpr.gt.0) write(lfnpr,920)
        eps = rellen
        maxit = iter
        return
      end if
      goto 10
c
  900 format(1x,'iter = ',i3,'    relative length = ',f10.7)
  910 format(1x,'no convergence within the specified number of ',
     + 'iterations')
  920 format(1x,'the error vector is converged')
      end
c*****************************************************************************
      subroutine subst(x,w,b,ipivot,n,ndim)
c*****************************************************************************
      implicit real*8 (a-h,o-z)
      dimension x(ndim),w(ndim,ndim),b(ndim),ipivot(ndim)
      data zero/0.0d0/
c
      if(n.eq.1) then
        x(1) = b(1)/w(1,1)
        return
      end if
c
c  use multipliers stored in w and back substitution to find x:
c
      ip = ipivot(1)
      x(1) = b(ip)
      do 20 i = 2,n
        sum = zero
        do 10 j = 1,i-1
          sum = w(i,j)*x(j) + sum
   10   continue
        ip = ipivot(i)
        x(i) = b(ip) - sum
   20 continue
      x(n) = x(n)/w(n,n)
      do 40 i = n-1,1,-1
        sum = zero
        do 30 j = i+1,n
          sum = w(i,j)*x(j) + sum
   30   continue
        x(i) = (x(i) - sum)/w(i,i)
   40 continue
      return
      end
c*****************************************************************************
c
c                 e n d    o f    n b o    p r o g r a m
c
c*****************************************************************************


c***********************************************************************
c
c
c                          g  m  s  n  b  o
c
c
c                    gamess version of nbo program
c
c
c  driver routines:
c
c      subroutine runnbo(q)
c      subroutine feaoin(core,icore,nboopt)
c      subroutine delscf(a,ia)
c
c***********************************************************************
      subroutine runnbo(q)
c***********************************************************************
      implicit real*8 (a-h,o-z)
      dimension nboopt(10)
      real*8 q(*)
c
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +           lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +           lfndaf,lfndef
c
c  gamess common block:
c
cUK      common /iofile/ ir,iw,ip,is,ipk,idaf,nav,ioda(99)

c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
       integer igmem_alloc_all

cUK      common /fmcom/ core(1)
c
      lfnin = ird
      lfnpr = iwr
c
c  determine the amount of available memory for the nbo analysis.
c
cc      call valfm(icur)
cc      call gotfm(memory)

      icur=0
c
c  set nbo options.
c
      nboopt(1)  =  0
      nboopt(2)  =  0
      nboopt(3)  =  0
      nboopt(4)  =  0
      nboopt(5)  =  0
      nboopt(6)  =  0
      nboopt(7)  =  0
      nboopt(8)  =  0
      nboopt(9)  =  0
      nboopt(10) =  6
c
c  perform the npa/nbo/nlmo analyses.
c
      ibase = igmem_alloc_all(memory)
      call nbo(q(ibase),memory,nboopt)
      call gmem_free(ibase)
c
c  perform the energetic analysis.
c
   10 nboopt(1) = 2
      ibase = igmem_alloc_all(memory)
      call nboean(q(ibase),memory,nboopt,idone)
      call gmem_free(ibase)
      if(idone.ne.0) goto 20
      call delscf(q)
      nboopt(1) = 3
      ibase = igmem_alloc_all(memory)
      call nboean(q(ibase),memory,nboopt,idone)
      call gmem_free(ibase)
      goto 10
c
   20 call endnbo
      return
      end
c***********************************************************************
      subroutine feaoin(core,icore,nboopt)
c***********************************************************************
      implicit none
      real*8 core(*)
      integer icore(*),nboopt(10)
c
c ----------------------------------------------------------------------
c
c   this routine fetchs basis set information from the gamess common
c  blocks and stores it in the nbo common blocks and direct access file
c  (daf) for use by the nbo analysis.
c
c ----------------------------------------------------------------------
c
c  routine feaoin accesses the following records of the dictionary file:
c
c          2  ---   total energy
c         12  ---   ao overlap matrix
c         14  ---   ao fock matrix (alpha)
c         15  ---   ao to mo transformation matrix (alpha)
c         16  ---   ao density matrix (bond order matrix) (alpha)
c         18  ---   ao fock matrix (beta)
c         19  ---   ao to mo transformation matrix (beta)
c         20  ---   ao density matrix (bond order matrix) (beta)
c         23  ---   x dipole integrals
c         24  ---   y dipole integrals
c         25  ---   z dipole integrals
c
c ----------------------------------------------------------------------
c
c  nbo common blocks
c
c
      integer maxatm, maxbas
      parameter (maxatm = 750, maxbas = 4096)

      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho

      integer ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit

       integer iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseudX,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint
      common/nbopt/iwdm,iw3c,iwapol,iwhybs,iwpnao,iwtnao,iwtnab,
     + iwtnbo,iwfock,iwcubf,ipseudX,kopt,iprint,iwdetl,iwmulp,ichoos,
     + jcore,jprint(60)

      integer lctr, lang
      common/nbao/lctr(maxbas),lang(maxbas)

      integer iatno,ino,norbs,ll,lu,iznuc,iatcr
      common/nbatom/iatno(maxatm),ino(maxatm),norbs(maxatm),ll(maxatm),
     +       lu(maxatm),iznuc(maxatm),iatcr(maxatm)

      integer lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +        lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +        lfndaf,lfndef
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +           lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +           lfndaf,lfndef
c
c
c  gamess common blocks
c
cUK      parameter (mxgtot=5000, mxsh=1000, mxatm=50)
cUK      common /runlab/ title(10),a(mxatm),b(mxatm),bflab(2047)
cUK      common /xyzprp/ x(3),pad(35)
cUK      common /iofile/ ir,iw,ip,is,ipk,idaf,nav,ioda(99)
cUK      common /infoa / nat,ich,mul,num,nx,ne,na,nb,zan(mxatm),c(3,mxatm)
cUK      common /nshel / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),
cUK     *                cf(mxgtot),cg(mxgtot),
cUK     *                kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),
cUK     *                kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
cUK      common /scfopt/ scftyp,blktyp,maxit,mconv,nconv,npunch
cUK      common /ecp2  / clp(400),zlp(400),nlp(400),kfirst(mxatm,6),
cUK     *                klast(mxatm,6),lmax(mxatm),lpskip(mxatm),
cUK     *                izcore(mxatm)

      real*8 atitle(10)
      equivalence (ztitle(1),atitle(1))

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
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
c
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
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
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
c
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
      real*8 toang, ams
      integer ifau, nosymm, iseczz
      common/phycon/toang(30),ams(54),ifau,nosymm,iseczz

c
c  local arrays
c
      integer labels(20)
      character*8 wfns(6)
      logical o1e(6)
      real*8 array(10)
      logical wstate(6,6)
      character *8 zcomm(19), zbuff(10)
      real*8 eiga(maxorb), frocca(maxorb),eigb(maxorb), froccb(maxorb)
      real*8 zero, dum
      integer i, j, k, ii, nfile, iatom, iat, min, max, loop
      integer istate, iblkk
      integer l2, l3
      integer itmp, itvec, ieiga, ieigb, idiffa, idiffb, ncola, ncolb
      integer iblkha, iblkhb, iblkv, newbas
      integer ix1, ix2
      integer nexp
      integer lensec
c
c functions from gamess-uk
c
      integer isubst
      external isubst
c
c  obtain the following information:
c
c    rohf        =.true. if rhf open shell wavefunction
c                =.false. otherwise
c
c    uhf         =.true. if uhf wavefunction
c                =.false. otherwise
c
c    auhf        =.true. if spin-annihilated uhf wavefunction
c                =.false. otherwise
c
c    ci          =.true. if ci wavefunction
c                =.false. otherwise
c
c    open        =.true. if open shell wavefunction
c                =.false. otherwise
c
c    complx      =.true. if complex wavefunction
c                =.false. otherwise
c                (note: the program is not capable of handling this.)
c
c    natoms      number of atomic centers
c
c    ndim        dimension of matrices (overlap and density)
c
c    nbas        number of basis functions (.le.ndim)
c
c    ipseud      set to one if pseudopotentials are used.
c                (NB renamed ispeudX here due to clash with /infoa/)
c
c    iwcubf      this pertains only basis sets with f functions.
c
c                if cartesian f functions are input, set iwcubf to:
c                    0,  if these are to be transformed to the
c                        standard set of pure f functions
c                    1,  if these are to be transformed to the
c                        cubic set of pure f functions
c
c                if pure f functions are input, set to iwcubf to:
c                    0,  if these are standard f functions
c                    1,  if these are cubic f functions
c
c    iatno(i),i=1,natoms
c                list of atomic numbers
c
c    lctr(i),i=1,nbas
c                list of atomic centers of the basis functions
c                (lctr(3)=2 if basis function 3 is on atom 2)
c
c    lang(i),i=1,nbas
c                list of angular symmetry information for the basis
c                functions
c
c labels array contains nbo labels for the atomic orbitals
c
      data labels /
c
c          s
c          ---
     +     1,
c
c          px    py    pz
c          ---   ---   ---
     +     101,  102,  103,
c
c          dxx   dyy   dzz   dxy   dxz   dyz
c          ---   ---   ---   ---   ---   ---
     +     201,  204,  206,  202,  203,  205,
c
c          fxxx  fyyy  fzzz  fxxy  fxxz  fxyy  fxyz  fxzz  fyyz  fyzz
c          ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
     +     301,  307,  310,  302,  303,  304,  305,  306,  308,  309 /
c
c
c  wstate array contains the values which should be set in the nbo
c  common block /nbflag/ depending on wavefunction.
c
      data wstate /
c                     logical variable in common nbflag
c                rohf      uhf      ci      open     mcscf    auhf
c              -------   -------  ------   ------   ------   ------
c wavefunction
c        rhf
     +         .false., .false., .false., .false., .false., .false.,
c        uhf
     +         .false., .true. , .false., .true. , .false., .false.,
c        rohf
     +         .true. , .false., .false., .true. , .false., .false.,
c        gvb
     +         .true.,  .false., .false., .true. , .false., .false.,
c        mcscf
     +         .false., .false., .false., .false., .true. , .false.,
c        ci
     +         .false., .false., .true. , .false., .false., .false./
c
c
c  wavefunction types:
c
      data wfns  /'rhf','uhf','rohf','gvb','mcscf','ci'/
c
      data zero/0.0d0/
c
c  store job title on nbodaf:
c
      do 5 i = 1,10
        core(i) = atitle(i)
    5 continue
      nfile = 2
      call nbwrit(core,10,nfile)
c
c  get the number of atoms from nat and store the atomic numbers in
c  iatno and nuclear charges in iznuc.  (note: atomic numbers and
c  nuclear charges may not be equivalent if effective core potentials
c  (ecp) are used.)
c
      natoms = nat
      do 10 i = 1,nat
        iatno(i) = isubst(zaname(i))
        iznuc(i) = nint(czan(i))
        if(iznuc(i) .ne. iatno(i)) ipseudX = 1
   10 continue
c
c  katom array contains which atom the shell is on, kmin and kmax
c  determine the components in the shell by pointing to a range in the
c  labels array:
c
      ii = 0
      do 30 i = 1,nshell
        iatom = katom(i)
        min   = kmin(i)
        max   = kmax(i)
        do 20 j = min,max
          ii = ii + 1
          lctr(ii) = iatom
          lang(ii) = labels(j)
   20   continue
   30 continue
c
      nbas  = ii
      ndim  = nbas
c
c  inititialize various nbo options depending upon the wavefunction
c  type and basis set type.
c
c  first, turn off the complex orbitals, indicate that the pure set
c  of f functions is desired when transforming from the cartesian set.
c
      complx = .false.
      iwcubf = 0
      ortho  = .false.
c
c  next set up the wavefunction flags.
c
      do 50 i = 1,6
        istate = i
        if (zscftp.eq.wfns(i)) goto 60
   50 continue
      stop 'unknown wfntyp'
c
   60 rohf  = wstate(1,istate)
      uhf   = wstate(2,istate)
      ci    = wstate(3,istate)
      open  = wstate(4,istate)
      mcscf = wstate(5,istate)
      auhf  = wstate(6,istate)
c
c  no fock matrices for rohf, mcscf, or ci wavefunctions:
c
      if (rohf.or.mcscf.or.ci) iwfock = 0
c
c  expectation values of the fock operator are in atomic units:
c
      munit = 0
c
c  store natoms, ndim, nbas, munit, wavefunction flags, iswean:
c
      icore(1)  = natoms
      icore(2)  = ndim
      icore(3)  = nbas
      icore(4)  = munit
      icore(5)  = 0
      if(rohf)  icore(5)  = 1
      icore(6)  = 0
      if(uhf)   icore(6)  = 1
      icore(7)  = 0
      if(ci)    icore(7)  = 1
      icore(8)  = 0
      if(open)  icore(8)  = 1
      icore(9)  = 0
      if(mcscf) icore(9)  = 1
      icore(10) = 0
      if(auhf)  icore(10)  = 1
      icore(11) = 0
      if(ortho) icore(11) = 1
      icore(12) = 1
      nfile = 3
      call nbwrit(icore,12,nfile)
c
c  store iatno, iznuc, lctr, and lang on nbo daf:
c
      ii = 0
      do 70 i = 1,natoms
        ii = ii + 1
        icore(ii) = iatno(i)
   70 continue
      do 80 i = 1,natoms
        ii = ii + 1
        icore(ii) = iznuc(i)
   80 continue
      do 90 i = 1,nbas
        ii = ii + 1
        icore(ii) = lctr(i)
   90 continue
      do 95 i = 1,nbas
        ii = ii + 1
        icore(ii) = lang(i)
   95 continue
      nfile = 4
      call nbwrit(icore,2*natoms+2*nbas,nfile)
c
c  fetch the total energy from the dictionary file and store it on the
c  nbo daf:
c
cUK      nfile = 2
cUK      call daread(idaf,ioda,core,3,nfile,nav)
cUK      core(1) = core(3)
cUK      core(2) = core(3)

      call secget(isect(494),16,iblkk)
      call rdedx(array,10,iblkk,numdu)

      core(1) = array(3)
      core(2) = array(3)

      nfile = 8
      call nbwrit(core,2,nfile)
c
c  store the atomic coordinates on the nbo daf: (note that these
c  coordinates are used in the calculation of dipole moments. gamess
c  requires the cartesian origin to be at the center of mass!!)
c


cUK  !!! there is no translation here
cUK  (seems not to affect default analysis)

      i = 0
      do 110 iat = 1,natoms
        do 100 k = 1,3
          i = i + 1
cUK          core(i) = (c(k,iat) - x(k)) * toang(1)
          core(i) = c(k,iat)
  100   continue
  110 continue
      nfile = 9
      call nbwrit(core,3*natoms,nfile)
c
c  store the overlap matrix on the nbodaf:
c
cc      nfile = 12
      l2 = ndim*(ndim+1)/2
cc      call daread(idaf,ioda,core,l2,nfile,nav)

      do loop = 1,6
         o1e(loop) = .false.
      enddo
      o1e(1) = .true.

      call getmat(core,dum,dum,dum, dum, dum,
     +     array,ndim,o1e,ionsec)

      nfile = 10
      call nbwrit(core,l2,nfile)
c
c  store the density matrices on the nbodaf:
c
      l2 = ndim*(ndim+1)/2
      call secget(isect(497),19,ibl3pa)
      call rdedx(core,l2,ibl3pa,idaf)
      nfile = 20
      call nbwrit(core,l2,nfile)
c
      if(open) then
        ibl3pb = ibl3pa + lensec(l2) + lensec(ndim)
        call rdedx(core,l2,ibl3pb,idaf)
        nfile = 21
        call nbwrit(core,l2,nfile)
      end if
c
c  store the fock matrices on the nbodaf:
c
c  - needs more work
c
      if(iwfock.ne.0) then

         l2 = ndim*(ndim+1)/2
         call secget(isect(42),42,iblkk)
         call rdedx(core,l2,iblkk,numdu)

         nfile = 30
         call nbwrit(core,l2,nfile)
         if(open) then

            call reads(core,l2,numdu)
            nfile = 31
            call nbwrit(core,l2,nfile)
         end if
      end if
c
c  store the ao to mo transformation matrices on the nbodaf:
c
      l3 = ndim*ndim
      itmp = mouta
      itvec = 3
      call secget(itmp,itvec,iblkv)
      call getqp(zcomm,zbuff,eiga,frocca,nbas,newbas,ncola,ieiga,
     +     idiffa,maxorb,iblkv)
c
      call rdedx(core,l3,iblkv,idaf)
      if (.not.otran) call tdown(core,ilifq,core,ilifq,ncola)

      nfile = 40
      call nbwrit(core,l3,nfile)
      if(open) then
         itmp = moutb
         call secget(itmp,itvec,iblkv)
         call getqp(zcomm,zbuff,eigb,froccb,nbas,newbas,ncolb,
     +        ieigb,idiffb,maxorb,iblkv)
         call rdedx(core,l3,iblkv,idaf)
         if (.not.otran) call tdown(core,ilifq,core,ilifq,ncolb)
         nfile = 41
         call nbwrit(core,l3,nfile)
      end if
c
c  store the x,y,z dipole integrals on the nbodaf:
c

      do loop = 1,3
         o1e(loop) = .false.
      enddo
      do loop = 4,6
         o1e(loop) = .true.
      enddo

      l2 = ndim*(ndim+1)/2
      ix1 = l2 +1
      ix2 = 2*l2 +1
      call getmat(core, core, core, core, core(ix1), core(ix2),
     +     array,num,o1e,ionsec)

      do 120 i = 1,l2
         core(i) = core(i) * toang(1)
         core(ix1-1+i) = core(ix1-1+i) * toang(1)
         core(ix2-1+i) = core(ix2-1+i) * toang(1)
 120  continue

      nfile = 50
      call nbwrit(core,l2,nfile)
      nfile = 51
      call nbwrit(core(ix1),l2,nfile)
      nfile = 52
      call nbwrit(core(ix2),l2,nfile)
c
c  store the ao basis set info on the nbo daf:  (note that two integers
c  and three integer arrays are stored first.  also remember that icore
c  and core occupy the same memory.)
c
      nexp = 0
      do 150 i = 1,mxprim
         if(ex(i).eq.zero) goto 150
         nexp = i
  150 continue
      do 160 i = 1,2+3*nshell+5*nexp
         core(i) = zero
  160 continue
      icore(1) = nshell
      icore(2) = nexp
c
c  ncomp(i) -- the number of components in the ith shell:
c
      ii = 2
      do 170 i = 1,nshell
        ii = ii + 1
        icore(ii) = kmax(i) - kmin(i) + 1
  170 continue
c
c  nprim(i) -- the number of gaussian primitives in the ith shell:
c
      do 180 i = 1,nshell
        ii = ii + 1
        icore(ii) = kng(i)
  180 continue
c
c  nptr(i) -- pointer for the ith shell into the gaussian parameters,
c             exp, cs, cp, etc.:
c
      do 190 i = 1,nshell
        ii = ii + 1
        icore(ii) = kstart(i)
  190 continue
c
c  exp(i) -- orbital exponents indexed by nptr:
c
      do 200 i = 1,nexp
        ii = ii + 1
        core(ii) = ex(i)
  200 continue
c
c  cs,cp,cd,cf -- orbital coefficients:
c
      do 210 i = 1,nexp
        ii = ii + 1
        core(ii) = cs(i)
  210 continue
      do 220 i = 1,nexp
        ii = ii + 1
        core(ii) = cp(i)
  220 continue
      do 230 i = 1,nexp
        ii = ii + 1
        core(ii) = cd(i)
  230 continue
      do 240 i = 1,nexp
        ii = ii + 1
cUK - add f functions
        core(ii) = cf(i)
  240 continue
      nfile = 5
      call nbwrit(core,ii,nfile)
c
      return
      end
c***********************************************************************
      subroutine delscf(a)
c***********************************************************************
      implicit none
      logical new,error,seq
c
c nbo common blocks:
c
      integer ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit
      common/nbinfo/ispin,natoms,ndim,nbas,mxbo,mxao,mxaolm,munit

      logical rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho
      common/nbflag/rohf,uhf,ci,open,complx,alpha,beta,mcscf,auhf,ortho

      integer lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +        lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +        lfndaf,lfndef
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +           lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +           lfndaf,lfndef

c
c  gamess common blocks:
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
c for h0 at ibl7f only
c
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
c
c for direct SCF
       integer nshtri, nshblk, len171, kkk, ilen
       integer lensec
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
       integer igmem_alloc
c
c arguments
c
      real*8 a(*)
c
c
c local variables
c
      integer ipt1, ipt2,ipt3,ipt4,ipt5,ipt6
      integer i, ntri, nsq
      integer irdmat, iprefa, ifock, idmat
      real*8 en, ehf, ehfa, ehfb, ehf1, ehf2, edel
      real*8 two
c
c functions
c
      real*8 enucf, tracep
      external enucf, tracep
c
c data
c
      data two/2.0d0/
c
c-----------------------------------------------------------------------
c
c  set pointers:
c
      ntri = ndim*(ndim+1)/2
      nsq  = ndim*ndim
c
c  a(ipt1)   ---   density matrix (alpha)
c  a(ipt2)   ---   density matrix (beta)
c  a(ipt3)   ---   fock matrix (alpha)
c  a(ipt4)   ---   fock matrix (beta)
c  a(ipt5)   ---   core hamiltonian matrix
c  a(ipt6)   ---   integral buffer, scratch
c
      ilen = 1 + 6 * ntri
      if (odscf) then
       nshtri=nshell*(nshell+1)/2
       nshblk = lensec(nshtri)
       ilen = ilen + 2 * nshtri
      endif
      ipt1 = igmem_alloc(ilen)
c
      ipt2 = ipt1 + ntri
      ipt3 = ipt2 + ntri
      ipt4 = ipt3 + ntri
      ipt5 = ipt4 + ntri
      ipt6 = ipt5 + ntri
      if (odscf) then
       irdmat = ipt6 + ntri
       iprefa = irdmat + nshtri
      endif
c
c  rewind integral file (not for gamess-uk)
c
c  open the nbo direct access file
c
      new = .false.
      call nbopen(new,error)
      if(error) then
        write(lfnpr,900)
        stop
      end if
c
c  calculate nuclear repulsion energy:
c
      en = enucf(nat,czan,c)

      if(uhf) then
c
c  uhf wavefunction: fetch the nbo deletion density matrix and construct
c      and symmetrize the skeleton fock matrix:
c
c       call caserr('no uhf deletions')
        alpha = .true.
        beta = .false.
        call fenewd(a(ipt1))
c
        alpha = .false.
        beta = .true.
        call fenewd(a(ipt2))
c
c       direct or conventional uhf ?
        if (odscf) then
c       need to restore the shell screening data for jandk
          irest = 0
          m171t=171
          call secget(isect(471),m171t,ibl171)
          call rdedx(a(iprefa),nshtri,ibl171,idaf)
          call reads(a(irdmat),nshtri,idaf)
          dlnmxd=-9999999.0d0
          do  kkk=1,nshtri
           if(a(irdmat+kkk-1).gt.dlnmxd) then
            dlnmxd=a(irdmat+kkk-1)
           endif
          enddo
c
          ifock = ipt3
          idmat = ipt1
          call dhstaru(a,a(ipt1),a(ipt2),a(ipt3),a(ipt4),
     +                   a(iprefa),a(irdmat),irest)

        else
          call hstaru(a(ipt1),a(ipt2),a(ipt3),a(ipt4),nopk)
        endif
        call symh(a(ipt3),a(ipt6),iky,0,0)
        call symh(a(ipt4),a(ipt6),iky,0,0)
c
c  read in core hamiltonian matrix and calculate the hf energy:
c
        call rdedx(a(ipt5),ntri,ibl7f,num8)
        call vadd(a(ipt3),1,a(ipt5),1,a(ipt3),1,ntri)
        call vadd(a(ipt4),1,a(ipt5),1,a(ipt4),1,ntri)
        ehfa = tracep(a(ipt1),a(ipt5),nbas) +
     +         tracep(a(ipt1),a(ipt3),nbas)
        ehfb = tracep(a(ipt2),a(ipt5),nbas) +
     +         tracep(a(ipt2),a(ipt4),nbas)
        ehf = (ehfa + ehfb)/two
        edel = ehf + en
c
c  rhf wavefunction: fetch the nbo deletion density matrix and construct
c      and symmetrize the skeleton fock matrix:
c
      else
        call fenewd(a(ipt1))
c       direct or conventional scf ?
c       need to restore the shell screening data for jandk
        if (odscf) then
          irest = 0
          m171t=171
          call secget(isect(471),m171t,ibl171)
          call rdedx(a(iprefa),nshtri,ibl171,idaf)
          call reads(a(irdmat),nshtri,idaf)
          dlnmxd=-9999999.0d0
          do  kkk=1,nshtri
           if(a(irdmat+kkk-1).gt.dlnmxd) then
            dlnmxd=a(irdmat+kkk-1)
           endif
          enddo
c
          ifock = ipt3
          idmat = ipt1
          call dhstar(a,a(ifock),a(idmat),a(iprefa),a(irdmat),irest)
        else
          call hstar(a(ipt1),a(ipt3),a(ipt6),nopk)
        endif
c
        call symh(a(ipt3),a(ipt6),iky,0,0)
c
c  read in core hamiltonian matrix and calculate the hf energy:
c
        call rdedx(a(ipt5),ntri,ibl7f,num8)
        call vadd(a(ipt3),1,a(ipt5),1,a(ipt3),1,ntri)

        ehf1 = tracep(a(ipt1),a(ipt5),ndim)
        ehf2 = tracep(a(ipt1),a(ipt3),ndim)

        ehf = (ehf1 + ehf2)/two
        edel = ehf + en
      end if
c
c  save the deletion energy on the nbo direct access file and close the
c  file:
c
      call sve0(edel)
      seq = .false.
      call nbclos(seq)
c
      call gmem_free(ipt1)
c
      return
c
  900 format(/1x,'error opening the nbo direct access file in ',
     + 'subroutine delscf.')
      end
c***********************************************************************
c
c           e n d    o f    g m s n b o    r o u t i n e s
c
c***********************************************************************
      subroutine ver_nbo(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/nbo.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end

      subroutine endnbo
      implicit real*8 (a-h,o-z)
      logical end,equal
      dimension keywd(4),kenter(4)
c
      common/nbio/lfnin,lfnpr,lfnao,lfnpna,lfnnao,lfnpnh,lfnnho,lfnpnb,
     +            lfnnbo,lfnpnl,lfnnlm,lfnmo,lfndm,lfnnab,lfnppa,lfnarc,
     +            lfndaf,lfndef
c
      data kenter/1he,1hn,1ht,1he/
c
c search input file for enter (if it is there) ...
c
   3  call strtin(lfnin)
      len=4
      call hfld(keywd,len,end)
      if (equal(keywd,kenter,4)) goto 10
      if(len.eq.0.and.end) goto 10
      goto 3
  10  return
      end
      subroutine lcase2(string)
      integer*1 string(160)
c
c...  convert to lower case (use ascii table)
c
      do 125 i=1,160,4
         mark=string(i)
         if (mark.ge.65.and.mark.le.90) string(i) = mark+32
125   continue
c
      return
      end
      character*8 function znameat(i)
c
c...   return the full name of the atom from runlab for nbo
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
      znameat = zaname(i)
      return
      end
