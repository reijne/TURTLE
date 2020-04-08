      subroutine gethes(h,kk10p,q,iwork,
     +     ic5e,ic6e,skv,brill,iky,kk10ps,work,qq,gam1,gam2)
c
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
c
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      common/craypk/i205(1360)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/actlen/lena(8),icfcor(100),ilifp(maxorb)
     1                          ,lenext(8)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      integer n20cas
      integer ifsimu, isimul, iam, noci, ipople
      parameter (n20cas = 100)
      common /simul/ ifsimu,isimul(n20cas),iam(n20cas),noci(n20cas),
     +               ipople
c
c
      integer mode, nrever, ifhes, isort, iaugm, isp, ibuffn
      common /ctrl/ mode,nrever,ifhes,isort(n20cas),iaugm,
     +              isp,ibuffn
c
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
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
      common/disc/isel,iselr,iselw,irep,ichek,ipos(maxlfn)
      common/qpar/ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,
     +            nst,n1e,isym,mults,iduqpar(5)
      common/blkin /g(510),nword
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
c... ic2e, ic3e and ic4e have room for 70 orbitals. Previously 15
c... orbitals could be accommodated for casscf. Now more are needed
c... for vb.
c
      integer ic1e, ic2e, ic3e, ic4e,maxcas
      parameter (maxcas=150)
      common /qice / ic1e(maxorb,2),ic2e(maxcas*(maxcas+1)/2+1),
     #               ic3e(maxcas*(maxcas+1)/2),ic4e(maxcas*(maxcas+1)/2)
c
      common/block /length,nsymm(8),nblock(8,maxorb)
c
      integer nirr, mult, isymao, isymmo, irrr, iss
      integer mstart, nfin, nmc
      logical symm_diag
      common /gjs/ nirr,mult(8,8),isymao(maxorb),
     + isymmo(maxorb),irrr,iss,mstart(8),nfin(8),nmc,
     + symm_diag
c
c
      character*10 charwall
      dimension work(*),iky(*),q(*),ij(2),qq(*),gam1(*),gam2(*)
      dimension brill(*),h(*),ic5e(*),ic6e(*),skv(*),iwork(ma,*)
c
c
      data half,donem,two,fourm/.5d0,-1.0d0,2.0d0,-4.0d0/
      data eight,quartm/8.0d0,-0.25d0/
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      t=cpulft(1)
      if(mode.ne.2.and.nprint.ne.-5)write(iwr,10)t ,charwall()
      if(mode.eq.2.and.nprint.ne.-5)write(iwr,20)t ,charwall()
c     nblkq=ma*ma
      kk10=kk10p-1
      kk10s=kk10ps-1
      lentri=ma*(ma+1)/2
      len2=lentri+lentri
      call rdfock(q,work,work(len2+1),skv,gam1)
      m2=2
      if(oprint(8))call fout(skv,gam1,m2,iwr)
      lentri=kk10s*(kk10s+1)/2
      do 30 in1=1,kk10p
30    iky(in1)=in1*(in1-1)/2
      iky(kk10p+1)=kk10p*(kk10p+1)/2
      call vclr(h,1,lentri)
c
      do 40 kk=1,kk10s
      ic5e(kk)=ic5e(kk+1)
      ic6e(kk)=ic6e(kk+1)
40    iwork(ic5e(kk),ic6e(kk))=kk
      call vclr(brill,1,kk10)
c
      nprim2=nprim*nprim
      call vclr(qq,1,nprim2)
      lfile=m6file
      do 670 i=1,lfile
      lotape(i)=m6tape(i)
      liblk(i)  =m6blk(i)
 670  llblk(i)  =liblk(i)-m6last(i)
c
      do 680 ii=1,lfile
      idev=lotape(ii)
      call search(liblk(ii),idev)
      call find(idev)
      jj=llblk(ii)
50    jj=jj+1
      call get(g(1),nw)
      if(nw)680,680,60
60    if(jj.ne.0)call find(idev)
      call unpack(g(num2e+1),lab816,i205,numlab)
      int4=1
      do 660 kk=1,nword
      j=i205(int4  )
      i=i205(int4+1)
      l=i205(int4+2)
      k=i205(int4+3)
      int4=int4+4
      ggg=g(kk)
c
      if(k.gt.nprim)goto70
      if(j.gt.nprim)goto160
      if(i.gt.nprim)goto260
      if(l.gt.ncore)goto440
      if(k.le.ncore)goto470
      if(j.gt.ncore)goto450
      goto540
70    if(j.gt.ncore.and.l.gt.ncore)goto 90
      if(j.gt.ncore)goto 120
      if(l.gt.ncore)goto 130
c
c.....(ai|bj) ints.
c
      if(isymmo(j).ne.isymmo(i).or.isymmo(k).ne.isymmo(l))goto80
      in=ind(iwork(j,i),iwork(l,k))
      h(in)=h(in)+eight*ggg
80    if(isymmo(j).ne.isymmo(k).or.isymmo(i).ne.isymmo(l))goto660
      in=ind(iwork(l,i),iwork(j,k))
      h(in)=h(in)-(ggg+ggg)
      goto660
c
c.....(at|bu) ints.
90    isymvx=mult(isymmo(j),isymmo(l))
      do 110 it=nst,nprim
      oitl=it.eq.l
      oitj=it.eq.j
      oitii=isymmo(it).ne.isymmo(i)
      oitik=isymmo(it).ne.isymmo(k)
      do 110 iu=nst,it
      oiui=isymmo(iu).ne.isymmo(i)
      oiuk=isymmo(iu).ne.isymmo(k)
      if(mult(isymmo(it),isymmo(iu)).ne.isymvx)goto110
      if(oitii) go to 100
      in1=icfcor(max(it,j))+min(it,j)
      in2=icfcor(max(iu,l))+min(iu,l)
      dumx=gam2(ic2e(max(in1,in2))+ic3e(min(in1,in2)))*ggg
      dumy=dumx
      if(oitj   )dumx=dumx+dumx
      if(iu.ne.l)dumx=dumx*half
      if(in1.eq.in2)dumx=dumx+dumx
      if(oiuk) go to 100
      in3=iky(iwork(it,i))+iwork(iu,k)
      h(in3)=h(in3)+dumx
c
100   if(it.eq.iu.or.oitik.or.i.eq.k)go to 110
      in1=icfcor(max(it,l))+min(it,l)
      in2=icfcor(max(iu,j))+min(iu,j)
      dumx=gam2(ic2e(max(in1,in2))+ic3e(min(in1,in2)))*ggg
      dumy=dumx
      if(oitl   )dumx=dumx+dumx
      if(iu.ne.j)dumx=dumx*half
      if(in1.eq.in2)dumx=dumx+dumx
      if(oiui) go to 110
      in3=iky(iwork(it,k))+iwork(iu,i)
      h(in3)=h(in3)+dumx
110   continue
c
      if(i.ne.k.or.j.le.l)goto660
c
c......recycle (av|ax) as (ax|av)
c
      in1=j
      j=l
      l=in1
      goto90
c
c
c.....(at|bi) ints.
120   k1=i
      l1=j
      i=k
      j=l
      k=k1
      l=l1
c
c.....(ai|bt)ints.
c
130   oij=isymmo(i).ne.isymmo(j)
      ojk=isymmo(j).ne.isymmo(k)
      do 150 it=nst,nprim
      if(isymmo(it).ne.isymmo(l))goto150
      dumx=gam1(ic1e(max(it,l),1)+ic1e(min(it,l),2))*(ggg+ggg)
      if(it.eq.l)dumx=dumx+dumx
      if(isymmo(it).ne.isymmo(k).or.oij)goto140
      in=iky(iwork(it,k))+iwork(j,i)
      h(in)=h(in)+dumx
140   if(isymmo(it).ne.isymmo(i).or.ojk)goto150
      in=iky(iwork(it,i))+iwork(j,k)
      h(in)=h(in)+dumx*quartm
150   continue
      goto660
c
160   if(l.gt.ncore)goto180
      if(k.gt.ncore)goto210
c
c.....(ab|ij) integrals
c
      ojk=isymmo(j).ne.isymmo(k)
      ojl=isymmo(j).ne.isymmo(l)
      if(isymmo(i).ne.isymmo(k).or.ojl)goto170
      in1=iky(iwork(k,i))+iwork(l,j)
      h(in1)=h(in1)-(ggg+ggg)
170   if(isymmo(i).ne.isymmo(l).or.k.eq.l.or.i.eq.j.or.ojk)
     +   go to 660
c
c
c.....(ab|tu) integrals
      in2=iky(iwork(k,j))+iwork(l,i)
      if(in1.ne.in2)h(in2)=h(in2)-(ggg+ggg)
      goto660
180   isymvx=mult(isymmo(k),isymmo(l))
      ivx=iky(k-ncore)+l-ncore
      do 200 it=nst,nprim
      do 200 iu=nst,it
      isymtu=mult(isymmo(it),isymmo(iu))
      if(isymtu.ne.isymvx)goto200
      oiui=isymmo(iu).ne.isymmo(i)
      oiuj=isymmo(iu).ne.isymmo(j)
      itu=iky(it-ncore)+iu-ncore
      in=ic2e(max(itu,ivx))+ic3e(min(itu,ivx))
      dumx=gam2(in)*ggg
      dumy=dumx
      if(itu.ne.ivx)dumx=dumx*half
      if(it.eq.iu)dumx=dumx+dumx
      if(isymmo(it).ne.isymmo(i).or.oiuj)goto190
      in=ind(iwork(iu,j),iwork(it,i))
      h(in)=h(in)+dumx
190   if(isymmo(it).ne.isymmo(j).or.i.eq.j.or.it.eq.iu.or.oiui)
     +   go to 200
      in=ind(iwork(it,j),iwork(iu,i))
      h(in)=h(in)+dumx
200   continue
      goto660
c
c
c.....(ab|ti) integrals
210   continue
      if(isymmo(l).ne.isymmo(i))goto230
      in1=iwork(l,i)
      do 220 it=nst,nprim
      if(isymmo(it).ne.isymmo(j))goto220
      in2=iwork(it,j)
      in=iky(in2)+in1
      dumx=ggg*gam1(ic1e(max(k,it),1)+ic1e(min(k,it),2))
      if(it.ne.k)dumx=dumx*half
      h(in)=h(in)-dumx
220   continue
230   if(i.eq.j.or.isymmo(l).ne.isymmo(j))goto250
      in1=iwork(l,j)
      do 240 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))goto240
      in2=iwork(it,i)
      in=iky(in2)+in1
      dumx=ggg*gam1(ic1e(max(k,it),1)+ic1e(min(k,it),2))
      if(it.ne.k)dumx=dumx*half
      h(in)=h(in)-dumx
240   continue
250   continue
      goto660
c
260   if(l.gt.ncore)goto 370
      if(j.le.ncore)goto310
      if(k.le.ncore)goto340
c
c.....(at|ui) integrals
c
      isymvx=mult(isymmo(j),isymmo(k))
      do 280 it=nst,nprim
      oitk = it .ne. k
      if(isymmo(it).ne.isymmo(l))goto280
      do 270 iu=nst,nprim
      if(isymmo(iu).ne.isymmo(i))goto270
      in1=icfcor(max(it,k))+min(it,k)
c
c......(ai|tj) integrals
      in2=icfcor(max(iu,j))+min(iu,j)
      dumx=ggg*gam2(ic2e(max(in1,in2))+ic3e(min(in1,in2)))
      if(oitk   )dumx=dumx*half
      if(iu.eq.j)dumx=dumx+dumx
      if(in1.eq.in2)dumx=dumx+dumx
      in1=iky(iwork(iu,i))+iwork(l,it)
      h(in1)=h(in1)-dumx
270   continue
280   continue
c
      olj=isymmo(j).ne.isymmo(l)
      olk=isymmo(k).ne.isymmo(l)
      do 300 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))goto300
      in1=iky(iwork(it,i))
      if(isymmo(it).ne.isymmo(k).or.olj)goto290
      in2=ic1e(max(it,k),1)+ic1e(min(it,k),2)
      in3=iwork(l,j)+in1
      dumx=gam1(in2)*ggg
      if(it.ne.k)dumx=dumx*half
      h(in3)=h(in3)-dumx
c
290   if(isymmo(it).ne.isymmo(j).or.olk)goto300
      dumx=gam1(ic1e(max(it,j),1)+ic1e(min(it,j),2))*(ggg+ggg)
      in3=iwork(l,k)+in1
      if(it.eq.j)dumx=dumx+dumx
      h(in3)=h(in3)+dumx
300   continue
      goto660
310   if(k.le.ncore)goto660
      oij=isymmo(i).ne.isymmo(j)
      oil=isymmo(i).ne.isymmo(l)
      do 330 it=nst,nprim
      if(isymmo(it).ne.isymmo(k))goto330
      in1=ic1e(max(it,k),1)+ic1e(min(it,k),2)
      dumx=gam1(in1)*half
      if(it.eq.k)dumx=(gam1(in1)-two)
      dumx=dumx*ggg
      if(isymmo(it).ne.isymmo(j).or.oil)goto320
      in3=ind(iwork(l,i),iwork(j,it))
      h(in3)=h(in3)+dumx
320   if(isymmo(it).ne.isymmo(l).or.oij)goto330
      in3=ind(iwork(l,it),iwork(j,i))
      h(in3)=h(in3)+dumx*fourm
330   continue
      goto660
c
c
c......(at|ij) integrals
 340  oik=isymmo(i).ne.isymmo(k)
      oil=isymmo(i).ne.isymmo(l)
      do 360 it=nst,nprim
      if(isymmo(it).ne.isymmo(j))goto360
      in1=ic1e(max(it,j),1)+ic1e(min(it,j),2)
      dumx=gam1(in1)*ggg*half
      if(it.eq.j)dumx=(gam1(in1)-two)*ggg
      if(isymmo(it).ne.isymmo(k).or.oil)goto350
      in2=ind(iwork(l,i),iwork(k,it))
      h(in2)=h(in2)+dumx
350   if(isymmo(it).ne.isymmo(l).or.k.eq.l.or.oik)goto360
      in2=iky(iwork(k,i))+iwork(l,it)
      h(in2)=h(in2)+dumx
360   continue
      goto660
c
370   if(j.gt.ncore)goto 420
c
c.....(ai|tu) integrals
c
      ivx=icfcor(k)+l
      isymvx=ic4e(ivx)
      do 390 it=nst,nprim
      oiti=isymmo(it).ne.isymmo(i)
      oitj=isymmo(it).ne.isymmo(j)
      do 390 iu=nst,it
      itu=icfcor(it)+iu
      if(ic4e(itu).ne.isymvx)goto390
      oiui=isymmo(iu).ne.isymmo(i)
      oiuj=isymmo(iu).ne.isymmo(j)
      dumx=ggg*gam2(ic2e(max(itu,ivx))+ic3e(min(itu,ivx)))
      if(itu.ne.ivx)dumx=dumx*half
      if(oitj.or.oiui) go to 380
      in1=iky(iwork(iu,i))+iwork(j,it)
      h(in1)=h(in1)-dumx
380   if(oiti.or.oiuj) go to 390
      in1=iky(iwork(it,i))+iwork(j,iu)
      h(in1)=h(in1)-dumx
390   continue
c
      ojk=isymmo(j).ne.isymmo(k)
      ojl=isymmo(j).ne.isymmo(l)
      do 410 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))goto410
      in3=iky(iwork(it,i))
      if(isymmo(it).ne.isymmo(l).or.ojk)goto400
      dumx=gam1(ic1e(max(it,l),1)+ic1e(min(it,l),2))*ggg
      if(it.ne.l)dumx=dumx*half
      in2=in3+iwork(j,k)
      h(in2)=h(in2)-dumx
400   if(isymmo(it).ne.isymmo(k).or.k.eq.l.or.ojl)goto410
      dumx=gam1(ic1e(max(it,k),1)+ic1e(min(it,k),2))*ggg
      if(it.ne.k)dumx=dumx*half
      in2=in3+iwork(j,l)
      h(in2)=h(in2)-dumx
410   continue
c
      goto660
c
c.....(au/vx) integrals
c
420   n2=icfcor(k)+l
c
      do 430 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))go to 430
      ll=iwork(it,i)
      n1=icfcor(max(it,j))+min(it,j)
      temp=ggg*gam2(ic2e(max(n1,n2))+ic3e(min(n1,n2)))
      if(j.ne.it)temp=temp*half
      if(n1.ne.n2)temp=temp*half
      brill(ll)=brill(ll)+temp
430    continue
      goto660
c
c-----------------------------------------------------
c.....(vi/tu)---> (tu/vi)
c
440   if(j.gt.ncore)goto610
       k1=k
      l1=l
      k=i
      l=j
      i=k1
      j=l1
c
c
c.....(tu/ij) integrals
450   n1=icfcor(i)+j
      do 460 it=nst,nprim
      if(isymmo(it).ne.isymmo(l))go to 460
      ll=iwork(l,it)
      n2=icfcor(max(it,k))+min(it,k)
      temp=ggg*gam2(ic2e(max(n1,n2))+ic3e(min(n1,n2)))
      if(k.ne.it)temp=temp*half
      if(n1.ne.n2)temp=temp*half
      brill(ll)=brill(ll)-temp
460   continue
      goto660
470   if(j.le.ncore)goto660
      ivx=icfcor(i)+j
      isymvx=ic4e(ivx)
      okl=k.eq.l
      do 490 it=nst,nprim
      oitk= isymmo(it).ne.isymmo(k)
      oitl= isymmo(it).ne.isymmo(l)
      do 490 iu=nst,it
      in1=icfcor(it)+iu
      if(ic4e(in1).ne.isymvx)goto490
      oiul=isymmo(l).ne.isymmo(iu)
      oiuk=isymmo(k).ne.isymmo(iu)
      dumx=gam2(ic2e(max(ivx,in1))+ic3e(min(ivx,in1)))*ggg
      if(it.ne.iu)dumx=dumx*half
      if(ivx.eq.in1)dumx=dumx+dumx
      if(oitk.or.oiul) go to 480
      in=iky(iwork(k,it))+iwork(l,iu)
      h(in)=h(in)+dumx
480   if(oiuk.or.it.eq.iu.or.okl.or.oitl)goto490
      in=iky(iwork(k,iu))+iwork(l,it)
      h(in)=h(in)+dumx
 490  continue
      oil=isymmo(i).ne.isymmo(l)
      oik=isymmo(i).ne.isymmo(k)
      ojl=isymmo(j).ne.isymmo(l)
      ojk=isymmo(j).ne.isymmo(k)
      do 530 it=nst,nprim
      if(isymmo(it).ne.isymmo(i).or.i.eq.j)goto510
      in1=ic1e(max(it,i),1)+ic1e(min(it,i),2)
      dumx=gam1(in1)*half
      if(it.eq.i)dumx=gam1(in1)+donem
      dumx=dumx*ggg
      if(isymmo(it).ne.isymmo(k).or.ojl.or.(k.eq.l.and.j.gt.it))
     +   go to 500
      in=ind(iwork(k,it),iwork(l,j))
      h(in)=h(in)+dumx
500   if(isymmo(it).ne.isymmo(l).or.ojk.or.(k.eq.l.and.it.gt.j))
     +   go to 510
      in=iky(iwork(k,j))+iwork(l,it)
      h(in)=h(in)+dumx
510   if(isymmo(it).ne.isymmo(j))goto530
      in1=ic1e(max(it,j),1)+ic1e(min(it,j),2)
      dumx=gam1(in1)*half
      if(it.eq.j)dumx=gam1(in1)+donem
      dumx=dumx*ggg
      if(isymmo(it).ne.isymmo(k).or.oil.or.(k.eq.l.and.i.gt.it))
     +   go to 520
      in=ind(iwork(k,it),iwork(l,i))
      h(in)=h(in)+dumx
520   if(isymmo(it).ne.isymmo(l).or.oik.or.(k.eq.l.and.it.gt.i))
     +   go to 530
      in=iky(iwork(k,i))+iwork(l,it)
      h(in)=h(in)+dumx
530   continue
      goto660
c
c.....(ti|uj) integrals
c
540   ivx=icfcor(i)+k
      isymvx=ic4e(ivx)
      do 560 it=nst,nprim
      oiti=it.ne.i
      oitk=it.ne.k
      oitj=isymmo(it).ne.isymmo(j)
      oitl=isymmo(it).ne.isymmo(l)
      do 560 iu=nst,it
      if(mult(isymmo(it),isymmo(iu)).ne.isymvx)goto560
      oiuj=isymmo(iu).ne.isymmo(j)
      oiul=isymmo(iu).ne.isymmo(l)
      in1=icfcor(max(iu,i))+min(iu,i)
      in2=icfcor(max(it,k))+min(it,k)
      dumx=gam2(ic2e(max(in1,in2))+ic3e(min(in1,in2)))*ggg
      if(iu.eq.i)dumx=dumx+dumx
      if(oitk   )dumx=dumx*half
      if(in1.eq.in2)dumx=dumx+dumx
      in1=icfcor(max(iu,k))+min(iu,k)
      in2=icfcor(max(it,i))+min(it,i)
      dumy=ggg*gam2(ic2e(max(in1,in2))+ic3e(min(in1,in2)))
      if(oiti   )dumy=dumy*half
      if(iu.eq.k)dumy=dumy+dumy
      if(in1.eq.in2)dumy=dumy+dumy
      if(oitj.or.oiul         )goto550
      in3=ind(iwork(j,it),iwork(l,iu))
      h(in3)=h(in3)+dumy
550   if((j.eq.l.and.i.eq.k).or.(it.eq.iu.and.j.ne.l))goto560
      if(oitl.or.oiuj        )goto560
      in3=ind(iwork(l,it),iwork(j,iu))
      h(in3)=h(in3)+dumx
560   continue
c
      okj=isymmo(j).ne.isymmo(k)
      okl=isymmo(k).ne.isymmo(l)
      oil=isymmo(i).ne.isymmo(l)
      oij=isymmo(i).ne.isymmo(j)
      do 600 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))goto580
      in1=ic1e(max(it,i),1)+ic1e(min(it,i),2)
      dumx=gam1(in1)*half
      if(it.eq.i)dumx=gam1(in1)+donem
      dumx=dumx*ggg
      if(j.eq.l.and.it.eq.k)dumx=dumx+dumx
      if(isymmo(it).ne.isymmo(l).or.okj)goto570
      in3=ind(iwork(l,it),iwork(j,k))
      h(in3)=h(in3)+dumx
570   if(isymmo(it).ne.isymmo(j).or.okl)goto580
      in3=ind(iwork(j,it),iwork(l,k))
      h(in3)=h(in3)+dumx*fourm
580   if(isymmo(it).ne.isymmo(k).or.(i.eq.k.and.j.eq.l))goto600
      in1=ic1e(max(it,k),1)+ic1e(min(it,k),2)
      dumx=gam1(in1)*half
      if(it.eq.k)dumx=gam1(in1)+donem
      dumx=dumx*ggg
      if(j.eq.l.and.it.eq.i)dumx=dumx+dumx
      if(isymmo(it).ne.isymmo(j).or.oil)goto590
      in3=ind(iwork(j,it),iwork(l,i))
      h(in3)=h(in3)+dumx
590   if(isymmo(it).ne.isymmo(l).or.oij)goto600
      in3=ind(iwork(l,it),iwork(j,i))
      h(in3)=h(in3)+dumx*fourm
600   continue
      goto660
610   iuv=icfcor(i)+j
      ixy=icfcor(k)+l
      do 650 it=nst,nprim
      if(isymmo(it).ne.isymmo(i))goto620
      in1=ilifp(it)+i
      in2=icfcor(max(it,j))+min(it,j)
      dumx=gam2(ic2e(max(ixy,in2))+ic3e(min(ixy,in2)))*ggg
      if(it.ne.j)dumx=dumx*half
      if(in2.ne.ixy)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
620   if(isymmo(it).ne.isymmo(j).or.i.eq.j)goto630
      in1=ilifp(it)+j
      in2=icfcor(max(it,i))+min(it,i)
      dumx=gam2(ic2e(max(in2,ixy))+ic3e(min(in2,ixy)))*ggg
      if(it.ne.i)dumx=dumx*half
      if(in2.ne.ixy)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
630   if(iuv.eq.ixy)goto650
      if(isymmo(it).ne.isymmo(k))goto640
      in1=ilifp(it)+k
      in2=icfcor(max(it,l))+min(it,l)
      dumx=ggg*gam2(ic2e(max(in2,iuv))+ic3e(min(in2,iuv)))
      if(it.ne.l)dumx=dumx*half
      if(in2.ne.iuv)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
640   if(isymmo(it).ne.isymmo(l).or.k.eq.l)goto650
      in1=ilifp(it)+l
      in2=icfcor(max(it,k))+min(it,k)
      dumx=gam2(ic2e(max(in2,iuv))+ic3e(min(in2,iuv)))*ggg
      if(it.ne.k)dumx=dumx*half
      if(in2.ne.iuv)dumx=dumx*half
      qq(in1)=qq(in1)+dumx
650   continue
c
660   continue
      if(jj)50,680,50
680   continue
      call dscal(kk10,two,brill,1)
      return
10    format(' commence first derivative evaluation at ',f8.2,
     * ' seconds',a10,' wall')
20    format(' commence first and second derivative',
     1' evaluation at ',f8.2,' seconds',a10,' wall')
      end
