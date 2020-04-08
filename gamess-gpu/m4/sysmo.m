      subroutine sysmowrite(q)
      implicit REAL (a-h,o-z)
c
      dimension q(*)
INCLUDE(common/restar)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/sysmo)

      write(iwr,100)
      write(iwr,105)
      write(iwr,110)thresh1,thresh2,maxits,mgm1,mgm2
      irenorm=igmem_alloc(num)
      call write_fort3(q,q(irenorm))
      call write_fort4(q,q(irenorm))
      call write_fort11(q,q(irenorm))
      call write_fort20(q)
      call write_fort23_28(q,q(irenorm))
      call gmem_free(irenorm)
  100 format(//40x,20('*')/40x,'SYSMO file interface'
     +       /40x,20('*'))
  105 format(/5x,'SYSMO is written by P. Lazzeretti and R. Zanasi',
     + /5x,'University of Modena',//,1x,104('+'))
  110 format(/5x,'SYSMO settings:',
     + /5x,'Convergence criterion                  : ',E10.4,
     + /5x,'Discard a-matrix elements smaller than : ',E10.4,
     + /5x,'Maximum number of iterations           : ',I10,
     + /5x,'Mgm1 setting                           : ',I10,
     + /5x,'Mgm2 setting                           : ',I10)
      return
      end
      subroutine write_fort3(q,drenorm)
      implicit REAL (a-h,o-z)
c
INCLUDE(common/restar)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
INCLUDE(common/iofile)
       common/scfblk/enuc,etotal
INCLUDE(common/dump3)
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      REAL czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne_infoa, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne_infoa,na,nb,czan(maxat),
     +            c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
INCLUDE(common/periodic)
INCLUDE(common/syminf)
INCLUDE(common/harmon)
INCLUDE(common/tran)
INCLUDE(common/gjs)
INCLUDE(common/machin)
INCLUDE(common/restri)
c     common/craypk/mmmm(65),isymao(maxorb),isymmo(maxorb)
      character *4 nogr,zelement,zirrep
      dimension eta(mxprim,5),ns(8),buf(9)
      dimension ideg(3), iddeg(6),inv(8),imp(8),ndri(8)
      dimension dtras(3,8),ddtras(6,8),tcar(3,3,8),m(maxorb,8)
      character*4 its(8)
      dimension q(*),drenorm(*)
      data m51/51/
c
c eerste record
c
      nogr='game'
      ngaus=0
      do i=1,nshell
         mini=kmin(i)
         maxi=kmax(i)
         ngaus=ngaus+(maxi-mini+1)*kng(i)
      enddo
      write(3)nogr,nt,nat,num,ngaus,enuc
c
c tweede record
c
      do i=1,nat
         zelement=zelem(int(czan(i)))
         write(3)zelement,c(1,i),c(2,i),c(3,i),czan(i)
      enddo
c
c derde record
c
      ngaus2=0
      do i=1,nshell
         mini=kmin(i)
         maxi=kmax(i)
         itype=maxi-mini+1
         if (itype.eq.1) then
            ngaus2=ngaus2+1
            do kk=kstart(i),kstart(i)+kng(i)-1
             ll=kk-kstart(i)
            eta(ngaus2+ll,1)=c(1,katom(i))
            eta(ngaus2+ll,2)=c(2,katom(i))
            eta(ngaus2+ll,3)=c(3,katom(i))
            eta(ngaus2+ll,4)=ex(kk)
            eta(ngaus2+ll,5)=cs(kk)
            enddo
            ngaus2=ngaus2-1+kng(i)
         elseif(itype.eq.3) then
            ngaus2=ngaus2+1
            do kk=kstart(i),kstart(i)+kng(i)-1
              ll=kk-kstart(i)
              do jj=1,3
               eta(ngaus2+ll+(jj-1)*kng(i),1)=c(1,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),2)=c(2,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),3)=c(3,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),4)=ex(kk)
               eta(ngaus2+ll+(jj-1)*kng(i),5)=cp(kk)
              enddo
            enddo
            ngaus2=ngaus2-1+3*kng(i)
         elseif(itype.eq.4) then
c... sp shell
            ngaus2=ngaus2+1
            do kk=kstart(i),kstart(i)+kng(i)-1
              ll=kk-kstart(i)
            eta(ngaus2+ll,1)=c(1,katom(i))
            eta(ngaus2+ll,2)=c(2,katom(i))
            eta(ngaus2+ll,3)=c(3,katom(i))
            eta(ngaus2+ll,4)=ex(kk)
            eta(ngaus2+ll,5)=cs(kk)
            enddo
            ngaus2=ngaus2-1+kng(i)
            ngaus2=ngaus2+1
            do kk=kstart(i),kstart(i)+kng(i)-1
              ll=kk-kstart(i)
              do jj=1,3
               eta(ngaus2+ll+(jj-1)*kng(i),1)=c(1,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),2)=c(2,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),3)=c(3,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),4)=ex(kk)
               eta(ngaus2+ll+(jj-1)*kng(i),5)=cp(kk)
              enddo
            enddo
            ngaus2=ngaus2-1+3*kng(i)
         elseif(itype.eq.6) then
            ngaus2=ngaus2+1
            do kk=kstart(i),kstart(i)+kng(i)-1
              ll=kk-kstart(i)
              do jj=1,6
               eta(ngaus2+ll+(jj-1)*kng(i),1)=c(1,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),2)=c(2,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),3)=c(3,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),4)=ex(kk)
               eta(ngaus2+ll+(jj-1)*kng(i),5)=cd(kk)
              enddo
            enddo
            ngaus2=ngaus2-1+6*kng(i)
         elseif(itype.eq.10) then
            ngaus2=ngaus2+1
            do kk=kstart(i),kstart(i)+kng(i)-1
              ll=kk-kstart(i)
              do jj=1,10
               eta(ngaus2+ll+(jj-1)*kng(i),1)=c(1,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),2)=c(2,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),3)=c(3,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),4)=ex(kk)
               eta(ngaus2+ll+(jj-1)*kng(i),5)=cf(kk)
              enddo
            enddo
            ngaus2=ngaus2-1+10*kng(i)
         elseif(itype.eq.15) then
            ngaus2=ngaus2+1
            do kk=kstart(i),kstart(i)+kng(i)-1
              ll=kk-kstart(i)
              do jj=1,15
               eta(ngaus2+ll+(jj-1)*kng(i),1)=c(1,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),2)=c(2,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),3)=c(3,katom(i))
               eta(ngaus2+ll+(jj-1)*kng(i),4)=ex(kk)
               eta(ngaus2+ll+(jj-1)*kng(i),5)=cg(kk)
              enddo
            enddo
            ngaus2=ngaus2-1+15*kng(i)
         else
           call caserr('only s,p,d,f,g, and sp')
         endif
      enddo
      if (ngaus2.ne.ngaus)call caserr('i have lost/gained primitives')
      do i=1,ngaus
         write(3)(eta(i,j),j=1,5)
      enddo
c
c vierde record
c
      nfirst=1
      icount=0
      do i=1,nshell
         mini=kmin(i)
         maxi=kmax(i)
         itype=maxi-mini+1
         if (itype.eq.1) then
c s-function
            ityp=1
         elseif(itype.eq.3) then
c p-function
            ityp=2
         elseif(itype.eq.4) then
c sp-function
            ityp=1
         elseif (itype.eq.6) then
c d-function
            ityp=5
         elseif(itype.eq.10) then
c f-function
            ityp=11
         elseif(itype.eq.15) then
c g-function
            ityp=21
         else
            call caserr('unknown shell')
         endif
         do j=mini,maxi
           nlast=nfirst+kng(i)-1
           write(3)katom(i),ityp,nfirst,nlast,kng(i)
           nfirst=nfirst+kng(i)
           icount=icount+1
           drenorm(icount)=1.0d0
           if (ityp.ge.8.and.ityp.le.10) then
              drenorm(icount)=dsqrt(3.0d0)
           elseif (ityp.ge.14.and.ityp.le.19) then
              drenorm(icount)=dsqrt(5.0d0)
           elseif (ityp.eq.20) then
              drenorm(icount)=dsqrt(15.0d0)
           endif
           ityp=ityp+1
         enddo
      enddo
c
c vijfde record
c
      do i=1,3
        ideg(i)=i
      enddo
      call vclr(buf,1,9)
      buf(1)=1.0
      buf(5)=1
      buf(9)=1
      call GROUP1(BUF,NT,3,nogr,ideg)
      goto 115
      do i=1,6
        iddeg(i)=i
      enddo
      mis=0
      ckl=0.0d0
      skl=0.0d0
      write(3)(ideg(i),i=1,3),mis,ckl,skl,(iddeg(j),j=1,6)
c
c zesde record
c
      do i=1,nt
         inv(i)=i
      enddo
      write(3)(inv(i),i=1,nt)
c
c zevende record
c deze begrijp ik nog niet

      do i=1,3
         write(3)(dtras(i,j),j=1,nt),(ddtras(2*i-1,j),j=1,nt),
     *           (ddtras(2*i,j),j=1,nt),((tcar(i,k,j),k=1,3),j=1,nt)
      enddo
c
c achtste record
c weet ik ook nog niet
      imp(1)=1
      write(3)(imp(i),i=1,nt)
c
c negende record
c ook uitzoeken
 115  continue
      do i=1,num
         m(i,1)=i
         write(3)(m(i,j),j=1,nt)
      enddo
c
c tiende record
c
      ncombt=0
      ncombg=0
      nofun=0
      nocoe=0
      write(3)ncombt,ncombg,nofun,nocoe
c
c elfde recorde
c   
      do i=1,nt
         write(its(i),'(A3,I1)')'SYM',i
      enddo
      write(3)nt,(its(i),i=1,nt)
c
c twaalfde record
c
      do i=1,nt
         ndri(i)=1
      enddo
      write(3)(ndri(i),i=1,nt)
c
c dertiende record
c 
      do i=1,nt
         write(3)(mults(i,j),j=1,nt)
      enddo
c
c viertiende recorde
c symmetry adaption
      i10=igmem_alloc(num*num)
      call vclr(q(i10),1,num*num)
      do i=1,num
         il=ilifc(i)
         do j=1,ntran(i)
c           index=(itran(il+j)-1)*num+i
c           q(i10+index-1)=ctran(il+j) 
c           index=(itran(il+j)-1)*num+i
c           q(i10+index-1)=ctran(il+j) 
            index=(i-1)*num+itran(il+j)
            q(i10+index-1)=ctran(il+j) 
         enddo
      enddo
      do i=1,num
         index=i10+(i-1)*num-1
         write(3)(q(index+j),j=1,num)
      enddo
      call gmem_free(i10)
c
c vijftiende record
c
      nav=lenwrd()
      call secget(isect(490),m51,iblk51)
      call readi(nirr,mach(13)*nav,iblk51,idaf)
      do i=1,8
         ns(i)=0
      enddo
      do i=1,num
         ns(isymao(i))=ns(isymao(i))+1
      enddo
      write(3)nt,(ns(i),i=1,nt)
c
c zestiende record met scf informatie
c
c     nmo=newbas0
      nmo=num
      norocc=na
      nb=nt
c     do i=1,8
c        ns(i)=0
c     enddo
c     do i=1,na
c        ns(isymmo(i))=ns(isymmo(i))+1
c     enddo
      write(3)nmo,norocc,nb,(ns(i),i=1,nb)
      close(3)
      return
      end
      subroutine write_fort4(q,dren)
      implicit REAL (a-h,o-z)
      dimension q(*)
INCLUDE(common/sizes)
INCLUDE(common/dump3)
INCLUDE(common/infoa)
c
INCLUDE(common/iofile)
_IF(ibm,vax)
      common/craypk/i205(340),j205(340),k205(340),l205(340)
_ELSE 
      common/craypk/i205(1360)
_ENDIF
INCLUDE(common/filel)
INCLUDE(common/atmblk)
INCLUDE(common/zorac)
      common/blkin/g(510),nnn
c
      logical ostf(6)
      dimension potnuc(10)
      dimension labels(10000),valint(5000)
      dimension dren(*)
      ind3(i,j)=max0(i,j)*(max0(i,j)-1)/2+min0(i,j)
      l1=num*(num+1)/2
      i10=igmem_alloc(l1)
      i20=igmem_alloc(l1)
      i30=igmem_alloc(l1)
      i40=igmem_alloc(l1)
      ostf(1)=.true.
      ostf(2)=.true.
      ostf(3)=.true.
      ostf(4)=.false.
      ostf(5)=.false.
      ostf(6)=.false.
c
c i10 - s
c i20 - t
c i30 - t+v
c
      call getmat(q(i10),q(i20),q(i30),q(i40),q(i40),
     *q(i40),potnuc,num,ostf,ionsec)
c     do i=1,l1
c        q(i30+i-1)=q(i30+i-1)-q(i20+i-1)
c     enddo
c now i30 = v
c
c add zora corrections to t-integrals
c     if (ozora) call zora(q,q,q(i20),'read')
      do i=1,num
         do j=1,i
            q(i10-1+ind3(i,j))=q(i10-1+ind3(i,j))/(dren(i)*dren(j))
         enddo
      enddo
      nint=0
      irec=0
      nrec=l1/5000
      if (mod(l1,5000).ne.0)nrec=nrec+1
      iflst=0
      do i=1,num
         do j=1,i
            nint=nint+1
            index=ind3(i,j)
            valint(nint)=q(i10+index-1)
            labels(nint)=0
            labels(nint)=ishft(i,10)
            labels(nint)=ishft(ior(labels(nint),j),10)
            if (nint.eq.5000) then
               irec=irec+1
               if (nrec.eq.irec) iflst=1
               write(4)nint,iflst,(labels(kk),kk=1,nint),
     *                 (valint(kk),kk=1,nint)
               nint=0
            endif
         enddo
      enddo
      if (nint.ne.0) then
         irec=irec+1
         if (nrec.eq.irec) iflst=1
         write(4)nint,iflst,(labels(kk),kk=1,nint),
     *           (valint(kk),kk=1,nint)
         nint=0
      endif
      if (iflst.ne.1)call caserr('no last record')
      goto 111
c
      nint=0
      irec=0
      nrec=l1/5000
      if (mod(l1,5000).ne.0)nrec=nrec+1
      iflst=0
      do i=1,num
         do j=1,i
            nint=nint+1
            index=ind3(i,j)
            valint(nint)=q(i20+index-1)
            labels(nint)=0
            labels(nint)=ishft(i,10)
            labels(nint)=ishft(ior(labels(nint),j),10)
            if (nint.eq.5000) then
               irec=irec+1
               if (nrec.eq.irec) iflst=1
               write(4)nint,iflst,(labels(kk),kk=1,nint),
     *                 (valint(kk),kk=1,nint)
               nint=0
            endif
         enddo
      enddo
      if (nint.ne.0) then
         irec=irec+1
         if (nrec.eq.irec) iflst=1
         write(4)nint,iflst,(labels(kk),kk=1,nint),
     *           (valint(kk),kk=1,nint)
         nint=0
      endif
      if (iflst.ne.1)call caserr('no last record')
      nint=0
      irec=0
      nrec=l1/5000
      if (mod(l1,5000).ne.0)nrec=nrec+1
      iflst=0
      do i=1,num
         do j=1,i
            nint=nint+1
            index=ind3(i,j)
            valint(nint)=q(i30+index-1)
            labels(nint)=0
            labels(nint)=ishft(i,10)
            labels(nint)=ishft(ior(labels(nint),j),10)
            if (nint.eq.5000) then
               irec=irec+1
               if (nrec.eq.irec) iflst=1
               write(4)nint,iflst,(labels(kk),kk=1,nint),
     *                 (valint(kk),kk=1,nint)
               nint=0
            endif
         enddo
      enddo
      if (nint.ne.0) then
         irec=irec+1
         if (nrec.eq.irec) iflst=1
         write(4)nint,iflst,(labels(kk),kk=1,nint),
     *           (valint(kk),kk=1,nint)
         nint=0
      endif
      if (iflst.ne.1)call caserr('no last record')
c
c     call prtri(q(i10),num)
c     call prtri(q(i20),num)
c     call prtri(q(i30),num)
 111  call gmem_free(i40)
      call gmem_free(i30)
      call gmem_free(i20)
      call gmem_free(i10)
      close(4)
      return
c two-electron integrals
      call setsto(1360,0,i205)
c
      nint=0
      iflst=0
      do 20 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      call find(iunit)
      lbl=llblk(ifile)
 10   lbl=lbl+1
      call get(g(1),mw)
      if(mw.eq.0) go to 20
      if(lbl.ne.0) call find(iunit)
_IF(ibm,vax)
      call upak8v(g(num2e+1),i205)
_ELSE
      int4=1
      call unpack(g(num2e+1),lab816,i205,numlab)
_ENDIF
      do ik1 = 1, nnn
_IF(ibm,vax)
      labi = i205(ik1)
      labj = j205(ik1)
      labk = k205(ik1)
      labl = l205(ik1)
_ELSEIF(littleendian)
      labj = i205(int4)
      labi = i205(int4+1)
      labl = i205(int4+2)
      labk = i205(int4+3)
_ELSE
      labi = i205(int4)
      labj = i205(int4+1)
      labk = i205(int4+2)
      labl = i205(int4+3)
_ENDIF
_IFN(ibm,vax)
      int4 = int4 + 4
_ENDIF
         nint=nint+1
         valint(nint)=g(ik1)
         mu=mu_calc(labi,labj,labk,labl)
         labels(nint*2)=0
         labels(nint*2-1)=0
         labels(nint*2)=ior(ishft(ior(ishft(labl,10),mu),10),1)
         labels(nint*2-1)=ior(ishft(ior(ishft(labi,10),labj),10),labk)
         if (nint.eq.5000) then
            write(4)nint,iflst,(labels(kk),kk=1,nint*2),
     *             (valint(kk),kk=1,nint)
            nint=0
         endif
      enddo
      goto 10        
20    continue        
c
      iflst=1
      write(4)nint,iflst,(labels(kk),kk=1,nint*2),
     *       (valint(kk),kk=1,nint)
      close(4)
      NORE4=0                                                           05/07/93
      NORE14=0                                                          01/04/94
      WRITE(13,'(I10,2X,I10)')NORE4,NORE14    
      close(13)
      return
      end
      subroutine write_fort11(q,drenorm)
      implicit REAL (a-h,o-z)
      dimension q(*),drenorm(*)
INCLUDE(common/sizes)
INCLUDE(common/common)
INCLUDE(common/atmol3)
INCLUDE(common/infoa)
INCLUDE(common/harmon)
INCLUDE(common/dump3)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
      data m1/1/
      character*8 st,sg1,sg2
      data sg1,sg2/'OCCUPATO','VIRTUALE'/
      ind3(i,j)=max0(i,j)*(max0(i,j)-1)/2+min0(i,j)
      l1=num
      l2=num*(num+1)/2
      l3=num*num
      i10=igmem_alloc(l3)
      i20=igmem_alloc(l1)
      i30=igmem_alloc(l1)
      call getq(q(i10),q(i20),q(i30),num,num,m1,m1,m1,mouta,scftyp)
      call tdown(q(i10),ilifq,q(i10),ilifq,num)
      call renorm_sys(q(i10),drenorm,num)
      do i=1,num
         if (q(i30+i-1).gt.1.999) then
           st=sg1
         else
           st=sg2
c          q(i20+i-1)=1.0d0+q(i20+i-1)
         endif
         write(11)q(i20+i-1),st,(q(i10+(i-1)*num+j-1),j=1,num)
      enddo
      len2=lensec(l1)
      m9=9
      call secput(isect(9),m9,len2,iblkm9)
      call wrt3(q(i20),l1,iblkm9,idaf)
      call gmem_free(i30)
      call gmem_free(i20)
      call gmem_free(i10)
      i10=igmem_alloc(l2)
      call rdedx(q(i10),l2,ibl3pa,idaf)
      call dscal(l2,0.5d0,q(i10),1)
      do i=1,num
         write(11)(q(i10+ind3(i,j)-1),j=1,i)
      enddo
      close(11)
      call gmem_free(i10)
      return
      end
      subroutine write_fort20(q)
      implicit REAL (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
      dimension q(*)
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
      ind3(i,j)=max0(i,j)*(max0(i,j)-1)/2+min0(i,j)
      nbasis=num
      nx=nbasis*(nbasis+1)/2
      gx=0.0
      gy=0.0
      gz=0.0
c angular momentum
      l2=nbasis*(nbasis+1)/2
      ix2=igmem_alloc(nbasis*(nbasis+1)/2)
      iy2=igmem_alloc(nbasis*(nbasis+1)/2)
      iz2=igmem_alloc(nbasis*(nbasis+1)/2)
      i100 = igmem_alloc(nx)
      i101 = igmem_alloc(nbasis*nbasis)
      i102 = igmem_alloc(nbasis)
      call amints(q(ix2),q(iy2),q(iz2),q(i100),q(i101),q(i102),
     &            .true.)
c
c dit begrijp ik niet....
c
c     call dscal(l2,2.0d0,q(ix2),1)
c     call dscal(l2,2.0d0,q(iy2),1)
c     call dscal(l2,2.0d0,q(iz2),1)
c
c     ncomp=3
c     write(20)'OP19ORIG',ncomp
c     do i=1,num
c        write(20)(q(ix2+ind3(i,j)-1),j=1,i)
c     enddo
c     do i=1,num
c        write(20)(q(iy2+ind3(i,j)-1),j=1,i)
c     enddo
c     do i=1,num
c        write(20)(q(iz2+ind3(i,j)-1),j=1,i)
c     enddo
c linear momentum
      call lmints(q(ix2),q(iy2),q(iz2),q(i100),q(i101),q(i102),
     &            .true.)
c     
c     call dscal(l2,-1.0d0,q(ix2),1)
c     call dscal(l2,-1.0d0,q(iy2),1)
c     call dscal(l2,-1.0d0,q(iz2),1)
c    
c     write(20)'OP21ORIG',ncomp
c     do i=1,num
c        write(20)(q(ix2+ind3(i,j)-1),j=1,i)
c     enddo
c     do i=1,num
c        write(20)(q(iy2+ind3(i,j)-1),j=1,i)
c     enddo
c     do i=1,num
c        write(20)(q(iz2+ind3(i,j)-1),j=1,i)
c     enddo
      call gmem_free(i102)
      call gmem_free(i101)
      call gmem_free(i100)
      call gmem_free(iz2)
      call gmem_free(iy2)
      call gmem_free(ix2)
      ncomp=0
      write(20)'**ULTIMO',ncomp
      close(20)
c
      return
      end
      function mu_calc(i,j,k,l)
      implicit REAL (a-h,o-z)
C...                                                           MUIJKL
C...                     _______ J>K ------ K>L ----- I>J>K>L     1
C...                    |           \______ K=L ----- I>J>K=L     2
C...              ____ I>K ----- J=K ------ K>L ----- I>J=K>L     3
C...             |      |           \______ K=L ----- I>J=K=L     4
C...      _____ I>J     |               /-- J>L ----- I>K>J>L     5
C...     |         \    |           /-K>L - J=L ----- I>K>J=L     6
C...     |          \    ------- J<K    \__ J<L ----- I>K>L>J     7
C...     |           \              \______ K=L ----- I>K=L>J     8
C...     I             I=K ---------------- J>L ----- I=K>J>L     9
C...     |                \________________ J=L ----- I=K>J=L    10
C...     |         /-- I>K ---------------- K>L ----- I=J>K>L    11
C...      ----- I=J       \________________ K=L ----- I=J>K=L    12
C...               \__ I=K ---------------- K>L ----- I=J=K>L    13
C...                      \________________ K=L ----- I=J=K=L    14
C...
C...
      mu=0
      if (i.eq.j.and.i.gt.k.and.k.gt.l) then
         mu=11
      elseif(i.eq.j.and.i.gt.k.and.k.eq.l) then
         mu=12
      elseif(i.eq.j.and.i.eq.k.and.k.gt.l) then
         mu=13
      elseif(i.eq.j.and.i.eq.k.and.k.eq.l) then
         mu=14
      elseif(i.gt.j.and.i.eq.k.and.j.gt.l) then
         mu=9
      elseif(i.gt.j.and.i.eq.k.and.j.eq.l) then
         mu=10
      elseif(i.gt.j.and.i.gt.k.and.j.gt.k.and.k.gt.l) then
         mu=1
      elseif(i.gt.j.and.i.gt.k.and.j.gt.k.and.k.eq.l) then
         mu=2
      elseif(i.gt.j.and.i.gt.k.and.j.eq.k.and.k.gt.l) then
         mu=3
      elseif(i.gt.j.and.i.gt.k.and.j.eq.k.and.k.eq.l) then
         mu=4
      elseif(i.gt.j.and.i.gt.k.and.j.lt.k.and.k.eq.l) then
         mu=8
      elseif(i.gt.j.and.i.gt.k.and.j.lt.k.and.k.gt.l.and.j.gt.l) then
         mu=5
      elseif(i.gt.j.and.i.gt.k.and.j.lt.k.and.k.gt.l.and.j.eq.l) then
         mu=6
      elseif(i.gt.j.and.i.gt.k.and.j.lt.k.and.k.gt.l.and.j.lt.l) then
         mu=7
      else
         call caserr('stupid mu....')
      endif
      if (mu.eq.0) call caserr('stupido mu equals 0')
      mu_calc=mu
      return
      end

      SUBROUTINE MOLTRI(A,B,C,N)
C...
C...  QUESTA SUBROUTINE ESEGUE IL PRODOTTO DI MATRICI
C...             C = A X B
C...  TUTTE QUADRATE E DI DIMENSIONI N.
C...
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C...
      DIMENSION A(N,N),B(N,N),C(N,N)
      DATA ZERO/0.0D+0/
C...
      DO 1 I=1,N
      DO 1 J=1,N
      SOM=ZERO
      DO 2 K=1,N
2     SOM=SOM+A(I,K)*B(K,J)
1     C(I,J)=SOM
      RETURN
      END

      SUBROUTINE GEND(A,B,IFINE)
C...
C...  QUESTA SUBROUTINE CONTROLLA SE I RAPPRESENTATIVI DI DUE
C...  OPERAZIONI DI SIMMETRIA, CONTENUTI NELLE MATRICI A E B,
C...  SONO UGUALI (IFINE=1) O NO (IFINE=0).
C...
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C...
      DIMENSION A(9),B(9)
      DATA ACRY/0.000001D+0/
C...
      IFINE =0
      DO 1 I=1,9
      IF(ABS(A(I)-B(I)).GT.ACRY)RETURN
 1    CONTINUE
      IFINE =1
      RETURN
      END

      SUBROUTINE GROUP1(BUF,NTRN,NOTAP3,nogr,ideg)
C...
C...  QUESTA ROUTINE,IN BASE ALLE RELAZIONI DI SIMMETRIA
C...  FRA LE COMPONENTI DI QUANTITA' TENSORIALI, DETERMINA
C...  LE MATRICI DI TRASFORMAZIONE VETTORIALI, UNA PER OGNI OPERAZIONE
C...  DI SIMMETRIA, DELLA DENSITA' ELETTRONICA E DELLO
C...  HAMILTONIANO DI FOCK PERTURBATI AL PRIMO ORDINE.
C...  LA ROUTINE DETERMINA ANCHE PER QUALI ROTAZIONI IMPROPRIE
C...  E' NECESSARIO CAMBIARE DI SEGNO NEL CASO DI PERTURBAZIONI
C...  POLARI E L'INVERSA DI OGNI OPERAZIONE DI SIMETRIA.
C...
C...  LA SCELTA INIZIALE ASSUME LA GEOMETRIA MOLECOLARE DISPOSTA IN
C...  MODO TALE CHE L'ASSE CARTESIANO Z SIA SEMPRE COINCIDENTE CON
C...  L'ASSE PRINCIPALE DI SIMMETRIA DI GRUPPI ASSIALI E CHE TALE
C...  ASSE SIA SEMPRE IL PRIMO GENERATORE.
C...  POSSONO VERIFICARSI I SEGUENTI CASI DESCRITTI DAL VETTORE IDEG:
C...   IDEG = 1 2 3  NON ESISTE NESSUNA DEGENERAZIONE;
C...   IDEG = 1 1 3  L'ASSE PRINCIPALE DI SIMMETRIA PROVOCA L'INTERSCAMBIO
C...                 (MISC=0) O LA COMBINAZIONE (MISC=1) DI X E Y CON
C...                 CONSEGUENTE DEGENERAZIONE DELLE RELATIVE COMPONENTI
C...                 VETTORIALI E TENSORIALI;
C...   IDEG = 1 1 1  IN QUESTO CASO IL GRUPPO DI SIMMETRIA NON E' ASSIALE
C...                 MA DI TIPO CUBICO, IL PRIMO GENERATORE E' SCELTO IN
C...                 MODO TALE CHE TRASFORMI X IN Y E Y IN Z CON
C...                 LE CONSEGUENTI DEGENERAZIONI.
C...
c...  The equivalent formalism for the quadrupole perturbation requires the
c...  iddeg array. There are just five possible cases, which are described
c...  in the routine rmtran of mo700.
c...
c...  the spherical tensor basis used is
c...  1    Theta-zz                  (1/2)    (3zz-rr)
c...  2    (1/sqrt3)Theta(xx-yy)     (sqrt3/2)(xx -yy)
c...  3    (2/sqrt3)Theta-xy         (sqrt3)  (xy)
c...  4    (2/sqrt3)Theta-yz         (sqrt3)  (yz)
c...  5    (2/sqrt3)Theta-xz         (sqrt3)  (xz)
c...  6    spherical combination     (1/sqrt3)(xx+yy+zz)
c...
c...  iddeg = 1 2 3 4 5 6   No equivalence
c...  iddeg = 1 2 3 4 4 6   xz,yz only (axial groups without mixing of x,y)
c...  iddeg = 1 2 2 4 4 6   xx-yy,xy and xz,yz (axial groups with mixing)
c...  iddeg = 1 1 3 3 3 6   zz,xx-yy and xy,yz,xz (cubic groups)
c...  iddeg = 1 1 1 1 1 6   icosahedral groups (mo700 CANNOT deal with them)
c...
c...  ddtras is an array which holds the transformation of the unique
c...  components of the quadrupole basis, in the same way that dtras
c...  stores the information for the dipole.
c...
c...  viene anche memorizzato su file il rappresentativo sulla base
c...  cartesiana di ogni operatore di simmetria
c...
C...  POICHE' TALI INFORMAZIONI SONO NECESSARIE AI PROGRAMMI CHF
C...  ESSE SONO SCRITTE SUL FILE DI OUTPUT NOTAP3.
C...
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C...
c     COMMON/DEG/IDEG(3)
      dimension ideg(3)
      DIMENSION BUF(18),IV(120),DTRAS(3,120),C(3,3)
cmo695
      dimension tcar(3,3,120)
cmo700 ...
      dimension iddeg(6),ddtras(6,120),t1(3,3),t2(6,6,120)
      character*4 nogr,ih1,ih2
      data ih1,ih2/'IH  ','I   '/
c...
      DATA ACRY/0.000001D+0/
      DATA ZERO,half,uno,three/0.0D+0,0.5d+0,1.0d+0,3.0d+0/
      root3=sqrt(three)
      roo32=half*root3
C...
C...  DETERMINA PER OGNI OPERAZIONE DI SIMMETRIA LA SUA INVERSA.
C...  L'INDICE DELL'INVERSA E' MESSO NEL VETTORE IV.
C...
      DO 4 I=1,NTRN
      IND=(I-1)*9+1
      DO 5 J=1,NTRN
      IND1=(J-1)*9+1
      CALL MOLTRI(BUF(IND),BUF(IND1),C,3)
      CALL GEND(C,BUF(1),IFINE)
      IF(IFINE.NE.0)GOTO 6
 5    CONTINUE
      GOTO 24
 6    IV(I)=J
 4    CONTINUE
C...
C...  IN FUNZIONE DELLE RELAZIONI TRA LE COMPONENTI (VEDI IDEG)
C...  STABILISCE SE VI SONO MISCELAMENTI (MISC=1) PIUTOSTO CHE
C...  SEMPLICI INTERSCAMBI (MISC=0) FRA ESSE, NEL CASO CHE MISC=1
C...  PRELEVA I COEFFICIENTI DELLA COMBINAZIONE COSENO E SENO
C...  DALLA TRASFORMAZIONE INVERSA DEL PRIMO GENERATORE.
C...
      MISC=0
      COSENO=ZERO
      SENO=ZERO
      IF(IDEG(1)+IDEG(2)+IDEG(3)-5)1,2,1
2     IF(ABS(BUF(11)).GT.ACRY.AND.ABS(BUF(14)).GT.ACRY)MISC=1
      J=IV(2)
      IND=(J-1)*9
      COSENO=BUF(IND+1)
      SENO=BUF(IND+4)
C...
cmo700 ...
c...  initialise iddeg for computation of the quadrupole representation
c...
1     continue
      IDDEG(1)=1
      IDDEG(2)=2
      IDDEG(3)=3
      IDDEG(4)=4
      IDDEG(5)=5
      IDDEG(6)=6
c... caso cubico.
      if(ideg(3).eq.1)then
      if(nogr.eq.ih1.or.nogr.eq.ih2)then
      IDDEG(1)=1
      IDDEG(2)=1
      IDDEG(3)=1
      IDDEG(4)=1
      IDDEG(5)=1
      IDDEG(6)=6
      else
      IDDEG(1)=1
      IDDEG(2)=1
      IDDEG(3)=3
      IDDEG(4)=3
      IDDEG(5)=3
      IDDEG(6)=6
      endif
      endif
c...  caso gruppo assiale.
      if(ideg(2).eq.1.and.ideg(3).eq.3)then
      if(misc.eq.0)then
c...  no mixing
      IDDEG(1)=1
      IDDEG(2)=2
      IDDEG(3)=3
      IDDEG(4)=4
      IDDEG(5)=4
      IDDEG(6)=6
      else
c...  mixing
      IDDEG(1)=1
      IDDEG(2)=2
      IDDEG(3)=2
      IDDEG(4)=4
      IDDEG(5)=4
      IDDEG(6)=6
      endif
      endif
C...  SCRIVE SU NOTAP3 IL TIPO DI DEGENERAZIONI ASSIEME AI COEFFICIENTI
C...  DI EVENTUALI MISCELAMENTI.
C...  SCRIVE SU NOTAP3 LE OPERAZIONI INVERSE.
C...
      WRITE(NOTAP3)(IDEG(I),I=1,3),MISC,COSENO,SENO,
     +(iddeg(i),i=1,6)
      WRITE(NOTAP3)(IV(I),I=1,NTRN)
C...
C...  CALCOLA LE MATRICI DI TRASFORMAZIONE VETTORIALI DI QUANTITA'
C...  TIPO HAMILTONIANO DI FOCK O DENSITA' ELETTRONICA PERTURBATO
C...  AL PRIMO ORDINE. METTE IL RISULTATO IN DTRAS.
C...
      DO 400 I=1,3
      IF(IDEG(I).NE.I)GOTO 400
      DO 399 J=I,3
      IF(IDEG(J).NE.I)GOTO 399
      DO 388 K=1,NTRN
      IND=(K-1)*9
 388  DTRAS(J,K)=BUF(IND+(J-1)*3+I)
 399  CONTINUE
 400  CONTINUE
c...
c...  T1 is the dipole representation matrix, and T2 is the matrix for
c...  the 6 cartesian products in normalised spherical tensor form
c...
      do 3771 k=1,ntrn
      ind=(k-1)*9
      do 3661 i=1,3
      do 3661 j=1,3
      t1(i,j)=buf(ind+(j-1)*3+i)
3661  continue
      t2(1,1,k)=half*(three*t1(3,3)**2-uno)
      t2(1,2,k)=roo32*(t1(3,1)**2-t1(3,2)**2)
      t2(1,3,k)=root3*t1(3,1)*t1(3,2)
      t2(1,4,k)=root3*t1(3,2)*t1(3,3)
      t2(1,5,k)=root3*t1(3,1)*t1(3,3)
      t2(1,6,k)=zero
      t2(2,1,k)=roo32*(t1(1,3)**2-t1(2,3)**2)
      t2(2,2,k)=half*(t1(1,1)**2+t1(2,2)**2-t1(1,2)**2-t1(2,1)**2)
      t2(2,3,k)=t1(1,1)*t1(1,2)-t1(2,1)*t1(2,2)
      t2(2,4,k)=t1(1,2)*t1(1,3)-t1(2,2)*t1(2,3)
      t2(2,5,k)=t1(1,1)*t1(1,3)-t1(2,1)*t1(2,3)
      t2(2,6,k)=zero
      t2(3,1,k)=root3*t1(1,3)*t1(2,3)
      t2(3,2,k)=t1(1,1)*t1(2,1)-t1(2,1)*t1(2,2)
      t2(3,3,k)=t1(1,1)*t1(2,2)+t1(1,2)*t1(2,1)
      t2(3,4,k)=t1(1,2)*t1(2,3)+t1(1,3)*t1(2,2)
      t2(3,5,k)=t1(1,1)*t1(2,3)+t1(1,3)*t1(2,1)
      t2(3,6,k)=zero
      t2(4,1,k)=root3*t1(2,1)*t1(3,1)
      t2(4,2,k)=t1(2,1)*t1(3,1)-t1(2,2)*t1(3,2)
      t2(4,3,k)=t1(2,1)*t1(3,2)+t1(2,2)*t1(3,1)
      t2(4,4,k)=t1(2,2)*t1(3,3)+t1(2,3)*t1(3,2)
      t2(4,5,k)=t1(2,1)*t1(3,3)+t1(2,3)*t1(3,1)
      t2(4,6,k)=zero
      t2(5,1,k)=root3*t1(1,1)*t1(3,1)
      t2(5,2,k)=t1(1,1)*t1(3,1)-t1(1,2)*t1(3,2)
      t2(5,3,k)=t1(1,1)*t1(3,2)+t1(1,2)*t1(3,1)
      t2(5,4,k)=t1(1,2)*t1(3,3)+t1(2,3)*t1(3,2)
      t2(5,5,k)=t1(1,1)*t1(3,3)+t1(1,3)*t1(3,1)
      t2(5,6,k)=zero
      t2(6,1,k)=zero
      t2(6,2,k)=zero
      t2(6,3,k)=zero
      t2(6,4,k)=zero
      t2(6,5,k)=zero
      t2(6,6,k)=uno
3771  continue
      do 4001 i=1,6
      if(iddeg(i).ne.i)goto 4001
      do 3991 j=i,6
      if(iddeg(j).ne.i)goto 3991
      do 3881 k=1,ntrn
 3881 ddtras(j,k)=t2(i,j,k)
 3991 continue
 4001 continue
c...
c...  trasferisce in tcar i rappresentativi delle
c...  operazioni di simmetria sulla base cartesiana
c...
      do 150 k=1,ntrn
      ind=(k-1)*9
      do 150 i=1,3
      do 150 j=1,3
      tcar(i,j,k)=buf(ind+(j-1)*3+i)
150   continue
c...
C...  SCRIVE DTRAS, ddtras e tcar SU NOTAP3
C...
      DO 3 I=1,3
3     WRITE(NOTAP3)(DTRAS(I,J),J=1,NTRN),
     +(ddtras(2*i-1,j),j=1,ntrn),(ddtras(2*i,j),j=1,ntrn),
     +((tcar(i,k,j),k=1,3),j=1,ntrn)
C...
C...  UTILIZZA IL VETTORE IV INTRODUCENDOVI 1 PER ROTAZIONI NORMALI
C...  E -1 PER ROTAZIONI IMPROPRIE.
C...
      DO 600 I=1,NTRN
      IND=(I-1)*9+1
      CALL DET_SYSMO(BUF(IND),D)
      IF(D)601,603,602
 601  IV(I)=-1
      GOTO 600
 602  IV(I)=1
 600  CONTINUE
C...
C...  SCRIVE IL VETTORE IV SU NOTAP3.
C...
      WRITE(NOTAP3)(IV(I),I=1,NTRN)
C...
      RETURN
C...
C...  USCITE PER ERRORE
C...
 24   WRITE(6,42)I
      CALL caserr('sysmo error')
 603  WRITE(6,306)I
      CALL caserr('sysmo error')
42    FORMAT(/1X,'TRASFORMAZIONE INVERSA NON TROVATA PER',I5)
306   FORMAT(/1X,'TROVATO UN DETERMINANTE NULLO PER',I5)
      END
 
      SUBROUTINE DET_SYSMO(A,D)
C...      
C...  CALCOLA IL DETERMINANTE D DI UNA MATRICE A 3X3.
C...      
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C...                 
      DIMENSION A(3,3)
      D=A(1,1)*A(2,2)*A(3,3) 
      D=D+A(1,2)*A(2,3)*A(3,1)
      D=D+A(1,3)*A(2,1)*A(3,2)
      D=D-A(3,1)*A(2,2)*A(1,3) 
      D=D-A(3,2)*A(2,3)*A(1,1) 
      D=D-A(3,3)*A(2,1)*A(1,2)
      RETURN         
      END            

      subroutine write_fort23_28(q,drenorm)
      implicit REAL (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/common)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/vectrn)
INCLUDE(common/iofile)
INCLUDE(common/machin)
INCLUDE(common/files)
INCLUDE(common/filel)
INCLUDE(common/atmblk)
INCLUDE(common/atmol3)
INCLUDE(common/mapper)
INCLUDE(common/harmon)
INCLUDE(../dft/common/ccpdft.hf77)
      common/blkin/g(510),nnn
_IF(ibm,vax)
      common/craypk/i205(340),j205(340),k205(340),l205(340)
_ELSE 
      common/craypk/i205(1360)
_ENDIF
      common/maxlen/maxq
      common/junke/ibl5,ibl54,ibl56,maxt,ires,ipass,nteff,
     1     npass1,npass2,lentri,nbuck,mloww,mhi,ntri,iacc,iontrn
      common/data1/vlist(100,3),cenmas(100),
     * newpro,nocp,mcen,mspace,itemp(100),jtemp(100)
     * ,nacta,nactb,louta(maxorb),loutb(maxorb)
     * ,norbt,norbta(maxorb),norbtp,
     *  norbc,norbca(maxorb),norbcp,ihamcp,ihamc
     * ,norb2,nseco(maxorb),norb3,nthird(maxorb),norbg,norbga(maxorb),
     *  norba,norbaa(maxorb),
     *  omulla,omullo,omullg,omullp,omullu,mullac(maxorb),mulla,oactrn
      common/junkc/ylab(26),ztype(35),zlist(100),zcentg(100),
     * ztagc(100),zhead(496),ztita(10),zatagg(maxat),zmul(maxorb)
INCLUDE(common/file12)
      common/lsort/iijjnk(4*maxorb),
     *pop(maxorb),potn,coren,ncolo,nbasq,newb,ncore,
     *mapcie(maxorb),ilifc(maxorb),nval,mapaie(maxorb),ilifa(maxorb),
     *iqsec,jdsec
      lind(i,j)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
c
      character*8 ULTIMO
      dimension ztda(113)
      equivalence (ztda(1),ztita(1))
      logical oexch,oiter
      data m3,m9/3,9/
      data m1,m2/1,2/
c
      dimension q(*),drenorm(*)
c
      oexch=.false.
      oiter=.false.
c
      dfacex=1.0d0
      if (CD_active()) then
         dfacex = 0.0d0
      else
         oexch=.true.
      endif
      if (CD_has_HF_exchange()) then
         write(iwr,'(1x,"Using ",f7.2,"% Hartree-Fock exchange")')
     &   100.0d0*CD_has_HF_exchange_weight()
         dfacex=CD_has_HF_exchange_weight()
         oexch=.true.
      endif
c
      if (oexch) call do_4index(q)
c
c magnetic hessian
c
c
c try to see if it all fits in core, if not go to low memory version
c
      l1=na*(newbas0-na)
      l2=l1*l1
      nx=num*(num+1)/2
c
      if (oexch.and.l1.gt.100) oiter=.true.
      if (oexch) then
         need=l2+2*num*num+(num*num+1)/2+2*num+4*l1+l1*l1
         if (igmem_max_memory().lt.need.or.oiter) then
            write(iwr,1101) 
1101  format (/'  Solving linear equations iteratively '/)
            call do_cphfsysmo(q,q,dfacex,drenorm)
            return
         endif
      endif
c
      if (oexch) then
         iamat=igmem_alloc(l2)
         call vclr(q(iamat),1,l2)
      else
         iamat=igmem_alloc(l1)
         call vclr(q(iamat),1,l1)
      endif
c
      ieps=igmem_alloc(num)
      call secget(isect(9),m9,isec9)
      call rdedx(q(ieps),num,isec9,ifild)
      ind1=0
      do iq = 1,na
         do ip = na+1,newbas0
             ind1=ind1+1
             if (ind1.gt.l1)call caserr('strange0')
             if (oexch) then
                ind4=(ind1-1)*l1+ind1
             else
                ind4=ind1
             endif
             q(iamat+ind4-1)=2.0d0*(q(ieps+ip-1)-q(ieps+iq-1))
         enddo
      enddo 
      call gmem_free(ieps)
c
c
      if (.not.oexch) goto 110
c
      nvirtual=newbas0-na
      call setsto(1360,0,i205)
      do 20 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      call find(iunit)
      lbl=llblk(ifile)
 10   lbl=lbl+1
      call get(g(1),mw)
      if(mw.eq.0) go to 20
      if(lbl.ne.0) call find(iunit)
_IF(ibm,vax)
      call upak8v(g(num2e+1),i205)
_ELSE
      int4=1
      call unpack(g(num2e+1),lab816,i205,numlab)
_ENDIF
      do ik1 = 1, nnn
_IF(ibm,vax)
      i = i205(ik1)
      j = j205(ik1)
      k = k205(ik1)
      l = l205(ik1)
_ELSEIF(littleendian)
      j = i205(int4)
      i = i205(int4+1)
      l = i205(int4+2)
      k = i205(int4+3)
_ELSE
      i = i205(int4)
      j = i205(int4+1)
      k = i205(int4+2)
      l = i205(int4+3)
_ENDIF
_IFN(ibm,vax)
      int4 = int4 + 4
_ENDIF
c        
         ind1=lind(i,j)
         ind2=lind(k,l)
         if (i.gt.na.and.j.gt.na) then
c (vv|oo)
          indp=i-na
          indq=k
          indr=j-na
          inds=l
          ii1=(indq-1)*nvirtual+indp
          ii2=(inds-1)*nvirtual+indr
        q(iamat-1+l1*(ii1-1)+ii2)=q(iamat-1+l1*(ii1-1)+ii2)
     1     -2.0d0*dfacex*g(ik1)
          if (ii1.ne.ii2) then
        q(iamat-1+l1*(ii2-1)+ii1)=q(iamat-1+l1*(ii2-1)+ii1)
     1     -2.0d0*dfacex*g(ik1)
         endif
         if (i.ne.j.and.k.ne.l) then
          ii1=(inds-1)*nvirtual+indp
          ii2=(indq-1)*nvirtual+indr
        q(iamat-1+l1*(ii1-1)+ii2)=q(iamat-1+l1*(ii1-1)+ii2)
     1     -2.0d0*dfacex*g(ik1)
        q(iamat-1+l1*(ii2-1)+ii1)=q(iamat-1+l1*(ii2-1)+ii1)
     1     -2.0d0*dfacex*g(ik1)
         endif 
         elseif (i.gt.na.and.j.le.na) then
c (vo|vo)
c
          indp=i-na
          indq=l
          indr=k-na
          inds=j
          ii1=(indq-1)*nvirtual+indp
          ii2=(inds-1)*nvirtual+indr
        q(iamat-1+l1*(ii1-1)+ii2)=q(iamat-1+l1*(ii1-1)+ii2)
     1     +2.0d0*dfacex*g(ik1)
          if (ii1.ne.ii2) then
        q(iamat-1+l1*(ii2-1)+ii1)=q(iamat-1+l1*(ii2-1)+ii1)
     1     +2.0d0*dfacex*g(ik1)
         endif

         else
            print *,'strange integral',i,j,k,l
             call caserr('strange integral')
         endif
      enddo
      goto 10
20    continue
c
c     print *,'amat'
c     call prsq(q(iamat),l1,l1,l1)
110   gx=0.0
      gy=0.0
      gz=0.0
c
      ilenna=84
      i0 = igmem_alloc(nx)
      icc1=igmem_alloc(num*num)
      icc0=igmem_alloc(num*num)
      ityp=3
      ix20=igmem_alloc(num)
      ix30=igmem_alloc(num)
      call getq(q(icc0),q(ix20),q(ix30),num,num,m1,m1,m1,mouta,scftyp)
      call gmem_free(ix30)
      call gmem_free(ix20)
      call tdown(q(icc0),ilifq,q(icc0),ilifq,num)
      call renorm_sys(q(icc0),drenorm,num)
c     write(6,'(A)')'Vectors from dumpfile'
c     call prsq(q(icc0),num,num,num)
c
      write(24)num,na,newbas0
      write(24)(q(icc0+i-1),i=1,num*num)
      if (oexch) then
         iaa=igmem_alloc(l1*l1)
         iwk1=igmem_alloc(l1)
         iwk2=igmem_alloc(l1)
      endif
      ixmat=igmem_alloc(l1)
      ibmat=igmem_alloc(l1)
c
      idump=igmem_alloc(num*(num+1)/2)
      call vclr(q(idump),1,num*(num+1)/2)
c
      do iop=1,2
      if (iop.eq.1) then
         ntape=28
         qfac=-1.0d0
      else
         ntape=23
         qfac=2.0d0
      endif
      do ijk=1,3
         itype=ilenna
         call secget(isect(ilenna),itype,isec)
         call rdedx(q(i0),nx,isec,idaf)
         ilenna=ilenna+1
         call vclr(q(ibmat),1,l1)
         ili=0
         do j=1,na
          do i=na+1,newbas0
              ikl=lind(i,j)
              ili=ili+1
              q(ibmat-1+ili)=-2.0d0*q(i0-1+ikl)
          enddo
         enddo
         call vclr(q(ixmat),1,l1)
         call vclr(q(icc1),1,na*num)
         if (oexch) then
         call f04atf(q(iamat),l1,q(ibmat),l1,q(ixmat),q(iaa),l1,
     1    q(iwk1),q(iwk2),ifail)
         if (ifail.ne.0) call caserr('ifail.ne.0')
         else
            do i=1,l1
               q(ixmat-1+i)=q(ibmat-1+i)/q(iamat-1+i)
            enddo
         endif
         ccc=qfac*ddot(l1,q(ixmat),1,q(ibmat),1)
         write(24)(q(ixmat+ii-1),ii=1,l1)
         print *,'component',iop,ijk,ccc
         indx=0
         do i=1,na
            do j=na+1,newbas0
               indx=indx+1
               do k=1,num
                  q(icc1-1+num*(i-1)+k)=q(icc1-1+num*(i-1)+k)+
     1                   qfac*q(ixmat-1+indx)*q(icc0-1+num*(j-1)+k)
               enddo
            enddo
         enddo
         call write_fort_ntape(ntape,q(idump),q(icc1),num,na)
      enddo
      ULTIMO='**ULTIMO'
      write(ntape)ULTIMO
      close(ntape)
      enddo
      call gmem_free(idump)
      call gmem_free(ibmat)
      call gmem_free(ixmat)
      if (oexch) then
         call gmem_free(iwk2)
         call gmem_free(iwk1)
         call gmem_free(iaa)
      endif
      call gmem_free(icc0)
      call gmem_free(icc1)
      call gmem_free(i0)
      call gmem_free(iamat)
      return
      end

      subroutine write_fort_ntape(ntape,f1,c1,num,nocc)
      implicit REAL (a-h,o-z)
      dimension f1(*),c1(num,nocc)
      character*8 ANCORA,ULTIMO
      lind(i,j)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
      ANCORA='**ANCORA'
      ULTIMO='**ULTIMO'
      write(ntape)ANCORA
      do i=1,num
         write(ntape)(f1(lind(i,j)),j=1,i)
      enddo
      do i=1,num
         write(ntape)(f1(lind(i,j)),j=1,i)
      enddo
      do i=1,nocc
         write(ntape)(c1(j,i),j=1,num)
      enddo
      return
      end

      subroutine do_4index(q)
      implicit REAL (a-h,o-z)
      dimension q(*)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/common)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/vectrn)
INCLUDE(common/iofile)
INCLUDE(common/machin)
INCLUDE(common/files)
INCLUDE(common/filel)
INCLUDE(common/atmblk)
INCLUDE(common/atmol3)
INCLUDE(common/mapper)
INCLUDE(common/harmon)
INCLUDE(../dft/common/ccpdft.hf77)
      common/blkin/g(510),nnn
_IF(ibm,vax)
      common/craypk/i205(340),j205(340),k205(340),l205(340)
_ELSE
      common/craypk/i205(1360)
_ENDIF
      common/maxlen/maxq
      common/junke/ibl5,ibl54,ibl56,maxt,ires,ipass,nteff,
     1     npass1,npass2,lentri,nbuck,mloww,mhi,ntri,iacc,iontrn
      common/data1/vlist(100,3),cenmas(100),
     * newpro,nocp,mcen,mspace,itemp(100),jtemp(100)
     * ,nacta,nactb,louta(maxorb),loutb(maxorb)
     * ,norbt,norbta(maxorb),norbtp,
     *  norbc,norbca(maxorb),norbcp,ihamcp,ihamc
     * ,norb2,nseco(maxorb),norb3,nthird(maxorb),norbg,norbga(maxorb),
     *  norba,norbaa(maxorb),
     *  omulla,omullo,omullg,omullp,omullu,mullac(maxorb),mulla,oactrn
      common/junkc/ylab(26),ztype(35),zlist(100),zcentg(100),
     * ztagc(100),zhead(496),ztita(10),zatagg(maxat),zmul(maxorb)
INCLUDE(common/file12)
      common/lsort/iijjnk(4*maxorb),
     *pop(maxorb),potn,coren,ncolo,nbasq,newb,ncore,
     *mapcie(maxorb),ilifc(maxorb),nval,mapaie(maxorb),ilifa(maxorb),
     *iqsec,jdsec
      lind(i,j)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
c
      character*8 ULTIMO
      dimension ztda(113)
      equivalence (ztda(1),ztita(1))
      data m3,m9/3,9/
      data m1,m2/1,2/
c
c
c
      ioldiscftp=iscftp
      call filchk(n4file,n4blk,n4last,n4tape,m2)
      call filchk(n6file,n6blk,n6last,n6tape,m3)
      acc1 = 10.0d0**(-iacc4)
      write (iwr,6720) acc1
c
      call setbfa
      norbt=newbas0
      norbc=0
      nb = na
      do i = 1 , norbt
         norbta(i) = i
      enddo
      nw1 = 117
      nav=lenwrd()
      len = lensec(mach(14)) + lensec(mach(16)) + lensec(nw1)
      call secput(isect(470),1005,len,isecbl)
      call wrt3i(norbt,mach(14)*nav,isecbl,idaf)
      call wrt3is(norb2,mach(16)*nav,idaf)
      call wrtcs(ztda,nw1,idaf)
      iscftp=2
      isecvv = isect(8)
      itypvv = 8
      npass1 = max(npas41,1)
      npass2 = max(npas42,1)
      iconvv = max(iconvv,9)
      lword4 = maxq
c
c   do the 4-index transformation and set
c
      ibase = igmem_alloc_all(mword)
      maxq = mword
      lword4 = maxq
c
      call indxsy(q(ibase))
c     call indx2t(q(ibase))
      call indx4t(q(ibase))
      call revise
      call gmem_free(ibase)
c
      iscftp=ioldiscftp
c
      lfile=m6file
      do i=1,lfile
       lotape(i) = m6tape(i)
       liblk(i)  = m6blk(i)
       llblk(i)  = liblk(i) - m6last(i)
      enddo
      return
 6720 format (/1x,'4-index accuracy factor',e12.1)
      end

      subroutine do_cphfsysmo(q,iq,dfacex,drenorm)
      implicit REAL (a-h,o-z)
      dimension q(*),drenorm(*),iq(*)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/common)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/vectrn)
INCLUDE(common/iofile)
INCLUDE(common/machin)
INCLUDE(common/files)
INCLUDE(common/filel)
INCLUDE(common/atmblk)
INCLUDE(common/atmol3)
INCLUDE(common/scra7)
INCLUDE(common/sysmo)
INCLUDE(common/mapper)
INCLUDE(common/harmon)
INCLUDE(../dft/common/ccpdft.hf77)
      common/blkin/g(510),nnn
      common/bloksort/ibb,numw
_IF(ibm,vax)
      common/craypk/i205(340),j205(340),k205(340),l205(340)
_ELSE
      common/craypk/i205(1360)
_ENDIF
      common/maxlen/maxq
      common/junke/ibl5,ibl54,ibl56,maxt,ires,ipass,nteff,
     1     npass1,npass2,lentri,nbuck,mloww,mhi,ntri,iacc,iontrn
      common/data1/vlist(100,3),cenmas(100),
     * newpro,nocp,mcen,mspace,itemp(100),jtemp(100)
     * ,nacta,nactb,louta(maxorb),loutb(maxorb)
     * ,norbt,norbta(maxorb),norbtp,
     *  norbc,norbca(maxorb),norbcp,ihamcp,ihamc
     * ,norb2,nseco(maxorb),norb3,nthird(maxorb),norbg,norbga(maxorb),
     *  norba,norbaa(maxorb),
     *  omulla,omullo,omullg,omullp,omullu,mullac(maxorb),mulla,oactrn
      common/junkc/ylab(26),ztype(35),zlist(100),zcentg(100),
     * ztagc(100),zhead(496),ztita(10),zatagg(maxat),zmul(maxorb)
INCLUDE(common/file12)
      common/lsort/iijjnk(4*maxorb),
     *pop(maxorb),potn,coren,ncolo,nbasq,newb,ncore,
     *mapcie(maxorb),ilifc(maxorb),nval,mapaie(maxorb),ilifa(maxorb),
     *iqsec,jdsec
      lind(i,j)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
c
      character*8 ULTIMO
      character*3 name
      dimension ztda(113)
      equivalence (ztda(1),ztita(1))
      data m3,m9/3,9/
      data m1,m2/1,2/
      data n511/511/
c
c
c
      l1=na*(newbas0-na)
      nx=num*(num+1)/2
c
      idiag=igmem_alloc(l1)
      ieps=igmem_alloc(num)
      call secget(isect(9),m9,isec9)
      call rdedx(q(ieps),lds(isect(9)),isec9,ifild)
      ind1=0
      do jq = 1,na
         do jp = na+1,newbas0
             ind1=ind1+1
             if (ind1.gt.l1)call caserr('strange0')
             q(idiag+ind1-1)=2.0d0*(q(ieps+jp-1)-q(ieps+jq-1))
         enddo
      enddo
      call gmem_free(ieps)
c
      nvirtual=newbas0-na
c
_IF(i8)
      ibucketinfo=igmem_alloc(2*l1)
      ibuck=ibucketinfo
      call vclr(q(ibucketinfo),1,2*l1)
_ELSE
      ibucketinfo=igmem_alloc(l1)
      call vclr(q(ibucketinfo),1,l1)
      ibuck=2*(ibucketinfo-1)+1
_ENDIF
      ibuck2=ibuck+l1
c
      nbuckets = l1
      minbucketsize = 512
      memavail = igmem_max_memory()
      nsize = min(memavail/nbuckets,max(l1,minbucketsize))
      npass = 1
      if (nsize .lt. minbucketsize) then
         npass = minbucketsize/nsize
         if (mod(minbucketsize,nsize).ne.0) npass = npass + 1
         nbuckets = l1/npass
         if (mod(l1,npass).ne.0) nbuckets = nbuckets + 1
         nsize = min(memavail/nbuckets,l1)
      endif
_IF(i8)
      maxinb=nsize/2-1
_ELSE
      maxinb=2*nsize/3 - 1
_ENDIF
      write(iwr,1102)npass
1102  format (/'  number of sort passes ',i4/)

c
      ibase=igmem_alloc(memavail)
_IF(i8)
      ibase_lab=ibase
_ELSE
      ibase_lab=2*(ibase-1)+1
_ENDIF
c
      numw=nsize
      name="sor"
      call opencc(name,3,2,ier)
      if (ier.ne.0)call caserr("during open sor file")
      minbuck=1
      do ipass=1,npass
         maxbuck=min(minbuck+nbuckets-1,l1)
         call setsto(1360,0,i205)
         do 20 ifile=1,lfile
         iunit=lotape(ifile)
         call search(liblk(ifile),iunit)
         call find(iunit)
         lbl=llblk(ifile)
 10      lbl=lbl+1
         call get(g(1),mw)
         if(mw.eq.0) go to 20
         if(lbl.ne.0) call find(iunit)
_IF(ibm,vax)
         call upak8v(g(num2e+1),i205)
_ELSE
         int4=1
         call unpack(g(num2e+1),lab816,i205,numlab)
_ENDIF
         do ik1 = 1, nnn
_IF(ibm,vax)
         i = i205(ik1)
         j = j205(ik1)
         k = k205(ik1)
         l = l205(ik1)
_ELSEIF(littleendian)
         j = i205(int4)
         i = i205(int4+1)
         l = i205(int4+2)
         k = i205(int4+3)
_ELSE
         i = i205(int4)
         j = i205(int4+1)
         k = i205(int4+2)
         l = i205(int4+3)
_ENDIF
_IFN(ibm,vax)
         int4 = int4 + 4
_ENDIF
c        
         ind1=lind(i,j)
         ind2=lind(k,l)
         if (i.gt.na.and.j.gt.na) then
c (vv|oo)
          indp=i-na
          indq=k
          indr=j-na
          inds=l
          ii1=(indq-1)*nvirtual+indp
          ii2=(inds-1)*nvirtual+indr
          if (ii1.ge.minbuck.and.ii1.le.maxbuck) then
             index_buck=ii1-minbuck
             iq(ibuck-1+ii1)=iq(ibuck-1+ii1)+1
             q(ibase-1+index_buck*nsize+iq(ibuck-1+ii1))=
     1         -2.0d0*dfacex*g(ik1)
_IF(i8)
         iad_lab=ibase_lab-1+index_buck*nsize+maxinb+iq(ibuck-1+ii1)
_ELSE
         iad_lab=ibase_lab-1+2*index_buck*nsize+2*maxinb+iq(ibuck-1+ii1)
_ENDIF
             iq(iad_lab)=ii2
             if (iq(ibuck-1+ii1).ge.maxinb) then
                call putsort(q(ibase+index_buck*nsize),q(ibase+
     1             index_buck*nsize),iq(ibuck-1+ii1),iq(ibuck2-1+ii1))
                iq(ibuck-1+ii1)=0
             endif
          endif
          if (ii1.ne.ii2) then
          ii3=ii1
          ii1=ii2
          ii2=ii3
          if (ii1.ge.minbuck.and.ii1.le.maxbuck) then
             index_buck=ii1-minbuck
             iq(ibuck-1+ii1)=iq(ibuck-1+ii1)+1
             q(ibase-1+index_buck*nsize+iq(ibuck-1+ii1))=
     1         -2.0d0*dfacex*g(ik1)
_IF(i8)
         iad_lab=ibase_lab-1+index_buck*nsize+maxinb+iq(ibuck-1+ii1)
_ELSE    
         iad_lab=ibase_lab-1+2*index_buck*nsize+2*maxinb+iq(ibuck-1+ii1)
_ENDIF       
             iq(iad_lab)=ii2
             if (iq(ibuck-1+ii1).ge.maxinb) then
                call putsort(q(ibase+index_buck*nsize),q(ibase+
     1             index_buck*nsize),iq(ibuck-1+ii1),iq(ibuck2-1+ii1))
                iq(ibuck-1+ii1)=0
             endif
          endif
         endif
         if (i.ne.j.and.k.ne.l) then
          ii1=(inds-1)*nvirtual+indp
          ii2=(indq-1)*nvirtual+indr
          if (ii1.ge.minbuck.and.ii1.le.maxbuck) then
             index_buck=ii1-minbuck
             iq(ibuck-1+ii1)=iq(ibuck-1+ii1)+1
             q(ibase-1+index_buck*nsize+iq(ibuck-1+ii1))=
     1         -2.0d0*dfacex*g(ik1)
_IF(i8)
         iad_lab=ibase_lab-1+index_buck*nsize+maxinb+iq(ibuck-1+ii1)
_ELSE    
         iad_lab=ibase_lab-1+2*index_buck*nsize+2*maxinb+iq(ibuck-1+ii1)
_ENDIF       
             iq(iad_lab)=ii2
             if (iq(ibuck-1+ii1).ge.maxinb) then
                call putsort(q(ibase+index_buck*nsize),q(ibase+
     1             index_buck*nsize),iq(ibuck-1+ii1),iq(ibuck2-1+ii1))
                iq(ibuck-1+ii1)=0
             endif
          endif
          ii3=ii1
          ii1=ii2
          ii2=ii3
          if (ii1.ge.minbuck.and.ii1.le.maxbuck) then
             index_buck=ii1-minbuck
             iq(ibuck-1+ii1)=iq(ibuck-1+ii1)+1
             q(ibase-1+index_buck*nsize+iq(ibuck-1+ii1))=
     1         -2.0d0*dfacex*g(ik1)
_IF(i8)
         iad_lab=ibase_lab-1+index_buck*nsize+maxinb+iq(ibuck-1+ii1)
_ELSE    
         iad_lab=ibase_lab-1+2*index_buck*nsize+2*maxinb+iq(ibuck-1+ii1)
_ENDIF       
             iq(iad_lab)=ii2
             if (iq(ibuck-1+ii1).ge.maxinb) then
                call putsort(q(ibase+index_buck*nsize),q(ibase+
     1             index_buck*nsize),iq(ibuck-1+ii1),iq(ibuck2-1+ii1))
                iq(ibuck-1+ii1)=0
             endif
          endif
         endif 
         elseif (i.gt.na.and.j.le.na) then
c (vo|vo)
c
          indp=i-na
          indq=l
          indr=k-na
          inds=j
          ii1=(indq-1)*nvirtual+indp
          ii2=(inds-1)*nvirtual+indr
          if (ii1.ge.minbuck.and.ii1.le.maxbuck) then
             index_buck=ii1-minbuck
             iq(ibuck-1+ii1)=iq(ibuck-1+ii1)+1
             q(ibase-1+index_buck*nsize+iq(ibuck-1+ii1))=
     1         2.0d0*dfacex*g(ik1)
_IF(i8)
         iad_lab=ibase_lab-1+index_buck*nsize+maxinb+iq(ibuck-1+ii1)
_ELSE    
         iad_lab=ibase_lab-1+2*index_buck*nsize+2*maxinb+iq(ibuck-1+ii1)
_ENDIF       
             iq(iad_lab)=ii2
             if (iq(ibuck-1+ii1).ge.maxinb) then
                call putsort(q(ibase+index_buck*nsize),q(ibase+
     1             index_buck*nsize),iq(ibuck-1+ii1),iq(ibuck2-1+ii1))
                iq(ibuck-1+ii1)=0
             endif
          endif
          if (ii1.ne.ii2) then
          ii3=ii1
          ii1=ii2
          ii2=ii3
          if (ii1.ge.minbuck.and.ii1.le.maxbuck) then
             index_buck=ii1-minbuck
             iq(ibuck-1+ii1)=iq(ibuck-1+ii1)+1
             q(ibase-1+index_buck*nsize+iq(ibuck-1+ii1))=
     1         2.0d0*dfacex*g(ik1)
_IF(i8)
         iad_lab=ibase_lab-1+index_buck*nsize+maxinb+iq(ibuck-1+ii1)
_ELSE    
         iad_lab=ibase_lab-1+2*index_buck*nsize+2*maxinb+iq(ibuck-1+ii1)
_ENDIF       
             iq(iad_lab)=ii2
             if (iq(ibuck-1+ii1).ge.maxinb) then
                call putsort(q(ibase+index_buck*nsize),q(ibase+
     1             index_buck*nsize),iq(ibuck-1+ii1),iq(ibuck2-1+ii1))
                iq(ibuck-1+ii1)=0
             endif
          endif
          endif
         else
            print *,'strange integral',i,j,k,l
             call caserr('strange integral')
         endif
         enddo
         goto 10
20       continue
      do izkl=minbuck,maxbuck
         index_buck=izkl-minbuck
         call putsort(q(ibase+index_buck*nsize),q(ibase+
     1             index_buck*nsize),iq(ibuck-1+izkl),iq(ibuck2-1+izkl))
                iq(ibuck-1+izkl)=0
      enddo
      minbuck=minbuck+nbuckets
      enddo
      call gmem_free(ibase)
c
c  start backchaining and write final rows of a-matrix on ed7
c
      ibinrow=ibl7la
      ibl7_amat=ibinrow+lensec(l1)
      iblkk7=ibl7_amat
c
      ibase=igmem_alloc(numw)
      irow=igmem_alloc(l1)
_IF(i8)
      ibase_lab=ibase
_ELSE
      ibase_lab=2*(ibase-1)+1
_ENDIF
      inrow=igmem_alloc(l1)
      nxblok=0
      do ibuck=1,l1
         call vclr(q(irow),1,l1)
         iorig=iq(ibuck2-1+ibuck)
 40      if (iorig.eq.0) goto 45
         call getsort(q(ibase),q(ibase),numi,iorig)
         do ints=1,numi
_IF(i8)
            index=iq(ibase_lab+maxinb+ints-1)
_ELSE
            index=iq(ibase_lab+2*maxinb+ints-1)
_ENDIF
            q(irow-1+index)=q(irow-1+index)+q(ibase-1+ints)
         enddo
         if (iorig.ne.0) goto 40
 45      q(irow-1+ibuck)=q(irow-1+ibuck)+q(idiag-1+ibuck)
         call dumprow(q(irow),l1,iblkk7,num8,q(inrow),ibuck,nbloks)
         nxblok=nxblok+nbloks
      enddo
      ibl7la=ibl7la + lensec(l1) + nxblok*lensec(n511)
      call wrt3(q(inrow),l1,ibinrow,num8)
      call gmem_free(inrow)
      call gmem_free(irow)
      call gmem_free(ibase)
      call gmem_free(ibucketinfo)
c     call gmem_free(idiag)
      call closecc(2)
      call delcc("sor",3,ier)
      if (ier.ne.0) call caserr("during delcc of sor file")
c
      gx=0.0
      gy=0.0
      gz=0.0
c
      ilenna=84
      i0 = igmem_alloc(nx)
      icc1=igmem_alloc(num*num)
      icc0=igmem_alloc(num*num)
      ityp=3
      ix20=igmem_alloc(num)
      ix30=igmem_alloc(num)
      call getq(q(icc0),q(ix20),q(ix30),num,num,m1,m1,m1,mouta,scftyp)
      call gmem_free(ix30)
      call gmem_free(ix20)
      call tdown(q(icc0),ilifq,q(icc0),ilifq,num)
      call renorm_sys(q(icc0),drenorm,num)
c     write(6,'(A)')'Vectors from dumpfile'
c     call prsq(q(icc0),num,num,num)
c
      ixmat=igmem_alloc(l1)
      ibmat=igmem_alloc(l1)
      iamat=igmem_alloc(l1)
      call rdedx(q(iamat),l1,ibinrow,num8)
c
      idump=igmem_alloc(num*(num+1)/2)
      call vclr(q(idump),1,num*(num+1)/2)
      ndim1=mgm2
      mgmres=mgm1
      eps=thresh1
      maxit=maxits
      iwork_size=l1*(2*ndim1+mgmres+2)
      iwork=igmem_alloc(iwork_size)
c
      do iop=1,2
      if (iop.eq.1) then
         ntape=28
         qfac=-1.0d0
      else
         ntape=23
         qfac=2.0d0
      endif
      do ijk=1,3
         itype=ilenna
         call secget(isect(ilenna),itype,isec)
         call rdedx(q(i0),nx,isec,idaf)
         ilenna=ilenna+1
         call vclr(q(ibmat),1,l1)
         ili=0
         do j=1,na
          do i=na+1,newbas0
              ikl=lind(i,j)
              ili=ili+1
              q(ibmat-1+ili)=-2.0d0*q(i0-1+ikl)
          enddo
         enddo
         call vclr(q(ixmat),1,l1)
         call vclr(q(icc1),1,na*num)
c
         do i=1,l1
            q(ixmat-1+i)=q(ibmat-1+i)/q(idiag-1+i)
         enddo
         maxit=maxits
         call gmresr(oktest,l1,ndim1,mgmres,q(ibmat),q(ixmat),q(iwork),
     1       eps,'abs',maxit,resid,iflag,q(iamat),ibl7_amat)
         if (iflag.eq.0) write(iwr,'(A,E10.4,A,I5)')'Convergence reached
     1, residue = ',resid, ' # iterations = ',maxit
         if (iflag.ne.0) then
         write(iwr,'(A,E10.4)')'No convergence reached, residue = ',
     1      resid
            call caserr('No convergence')
         endif
c
         ccc=qfac*ddot(l1,q(ixmat),1,q(ibmat),1)
         write(24)(q(ixmat+ii-1),ii=1,l1)
         print *,'component',iop,ijk,ccc
         indx=0
         do i=1,na
            do j=na+1,newbas0
               indx=indx+1
               do k=1,num
                  q(icc1-1+num*(i-1)+k)=q(icc1-1+num*(i-1)+k)+
     1                   qfac*q(ixmat-1+indx)*q(icc0-1+num*(j-1)+k)
               enddo
            enddo
         enddo
         call write_fort_ntape(ntape,q(idump),q(icc1),num,na)
      enddo
      ULTIMO='**ULTIMO'
      write(ntape)ULTIMO
      close(ntape)
      enddo
      call vclr(g,1,n511)
      call gmem_free(iwork)
      call gmem_free(idump)
      call gmem_free(iamat)
      call gmem_free(ibmat)
      call gmem_free(ixmat)
      call gmem_free(icc0)
      call gmem_free(icc1)
      call gmem_free(i0)
      call gmem_free(idiag)
      return
      end

      subroutine dumprow(row,nn,iblok,num8,ibbb,irow,nbloks)
      implicit REAL (a-h,o-z)
      dimension row(nn),ibbb(nn)
      common/blkin/g(340),ig(340),nnn
INCLUDE(common/sysmo)
      integer*4 ig
      data m511/511/
c
      index=0
      ibbb(irow)=0
      nnn=0
      do i=1,nn
         if (dabs(row(i)).gt.thresh2) then
            index=index+1
            g(index)=row(i)
            ig(index)=i
            if (index.ge.340) then
               nnn=index
               ibbb(irow)=ibbb(irow)+1
               call wrt3(g,m511,iblok,num8)
               iblok=iblok+lensec(m511)
               index=0
            endif
         endif
      enddo
      if (index.ne.0) then
         nnn=index
         ibbb(irow)=ibbb(irow)+1
         call wrt3(g,m511,iblok,num8)
         iblok=iblok+lensec(m511)
         index=0
      endif
      nbloks=ibbb(irow)
      return
      end

      subroutine renorm_sys(v,drenorm,num)
      implicit REAL (a-h,o-z)
      dimension v(num,num),drenorm(num)
      do i=1,num
         do j=1,num
            v(i,j)=v(i,j)*drenorm(i)
         enddo
      enddo
      return
      end

      subroutine matvec_sys(xmat,ymat,l1,inrow,iblk7)
      implicit REAL (a-h,o-z)
      dimension xmat(l1),ymat(l1),inrow(l1)
INCLUDE(common/iofile)
INCLUDE(common/scra7)
      common/blkin/g(340),ig(340),nnn
      integer*4 ig
      data m511/511/
c
      call vclr(ymat,1,l1)
      ibl=iblk7
c
      do i=1,l1
         do nb=1,inrow(i)
            call rdedx(g,m511,ibl,num8)
            ibl=ibl+lensec(m511)
            do j=1,nnn
               ymat(i)=ymat(i)+xmat(ig(j))*g(j)
            enddo
         enddo
      enddo
      return
      end

      subroutine gmresr(oktest,n,j,mgmres,b,x,work,eps,stc,
     &                  maxits,resid,iflag,amat,ibl7_amat)
C*********************************************************
C GMRESR algorithm to solve linear system Ax = b
C
C  author:
C  m.botchev, utrecht university, december 1996
C  report bugs to botchev@cwi.nl or botchev@excite.com
C
C Copyright (c) 1996 by M.A.Botchev
C Permission to copy all or part of this work is granted,
C provided that the copies are not made or distributed
C for resale, and that the copyright notice and this
C notice are retained.
C
C Details to the algorithm may be found in
C  H.A. van der Vorst, C. Vuik, "GMRESR: a Family of Nested GMRES
C  Methods", Num. Lin. Alg. Appl., vol. 1(4), 369--386 (1994)
C
C parameter list:
C oktest is TRUE if intermediate residuals should be printed
C n      == INTEGER size of problem
C j      == INTEGER truncation parameter (work with j last vectors)
C mgmres == INTEGER dimension of the envoked GMRES
C           if mgmres.eq.0 then we get GCR algorithm, the simplest
C           version of GMRESR 
C b      == DOUBLE PRECISION righthand side vector
C x      == DOUBLE PRECISION initial guess on input,
C           (approximate) solution on output
C work   == DOUBLE PRECISION work space 1 of size n x (2*j+mgmres+2)
C eps    == DOUBLE PRECISION tolerance of stopping criterion. 
C           process is stopped as soon as 
C           (1) residual norm has been dumped by factor eps, 
C           i.e.  ||res|| / ||res0|| <= eps   OR
C           (2) maximum number of iterations maxit has been performed
C stc    == CHARACTER*3
C           Determine stopping criterion (||.|| denotes the 2-norm):
C           stc='rel'    -- relative stopping crit.: ||res|| < eps*||res0||
C           stc='abs'    -- absolute stopping crit.: ||res|| < eps
C maxits == INTEGER max. no. outer_iterative_steps/truncation_length on input
C           on output it is the actual number of total iterative steps   
C resid  == DOUBLE PRECISION residual measure (depends on stopping criterion)
C           achieved on output 
C iflag  == INTEGER on output 0 - solution found within tolerance
C                             1 - no convergence within maxits
C ----------------------------------------------------------
C subroutines used
C   matvec   == matrix-vector product y <- A*x
C  blas subroutines:
C   dscal
C   daxpy
C  blas functions:
C   ddot
C   dnrm2 
C**********************************************************


C     list of variables: arrays in alphabetical order,
C     then other variables in alphabetical order

      logical oktest
      character*3 stc

      integer i,iflag,its,nits,itsinn,j,k,maxits,mgmres,n

      REAL b(n),x(n),work(n,0 : 2*j+mgmres+2 -1),
     &     alpha,alphai,cknorm,ckres,ddot,dnrm2,eps,epsinn,
     &     res0,resnor,resid,amat(n)

C distribute work space work(n,*) among some virtual variables;
C namely, we think of columns of work as being occupied by 
C c(n,0:j-1), u(n,0:j-1), resid(n), workgmr(n,mgmres+1)
C therefore we define "shifts"

      integer c, u, workres, workgmre,ibl7_amat
C ----------------------------------------------------------------------------

      if((stc.NE.'rel').and.(stc.NE.'abs'))then
         PRINT *,'Error in VACGMRESR:'
         PRINT *,'PARAMETER STC=',stc,' SHOULD BE rel OR abs.'
         STOP
      endif

C     c occupies j columns 0...j-1:
      c = 0 
C     u occupies j columns j...2*j-1:
      u = j
C     resid occupies 1 column No. 2*j:
      workres = 2*j    
C     workgmre occupies mgmres+1 columns 2*j+1...2*j+mgmres+1:
      workgmre = 2*j+1 
C     so that we can access, say, to the (k)th column of the virtual
C     array c(n,0:j-1) as work(1,c+k),
C     virtual residual vector resid(n) is work(1,workres) and so on ...

C ***Furthermore, we build sequences c_k and u_k, k = 0,...,m-1
C but we store only the last j vectors in the following way:
C Assume j=3, then
C --------------------------------------------------------------
C  k    |  number of column of work | columns of work which are vectors
C       |  in which we store c_k    |  we actually store and use
C  0    |           0               |   c_0             u_0            ...
C  1    |           1               |   c_0  c_1        u_0  u_1       ...
C  2    |           2               |   c_0  c_1  c_2   u_0  u_1  u_2  ...
C  3    |           0               |   c_3  c_1  c_2   u_3  u_1  u_2  ...
C  4    |           1               |   c_3  c_4  c_2   u_3  u_4  u_2  ... 
C  5    |           2               |   c_3  c_4  c_5   u_3  u_4  u_5  ...
C  6    |           0               |   c_6  c_4  c_5   u_6  u_4  u_5  ...
C ...   |           ...             |      ...               ...
C This mapping is performed by function mod(k,j)

C     Reset iteration counter
      nits= 0
      its = 0
      
C     Calculate (initial) residual norm
      call matvec_sys(x,work(1,workres),n,amat,ibl7_amat)
      alpha = -1
      call daxpy(n,alpha,b,1,work(1,workres),1)
      call dscal(n,alpha,work(1,workres),1)

C     Calculate residual norm and quit if it is zero
      res0 = dnrm2(n,work(1,workres),1)
      resnor = res0
      resid = 0

      if ( res0 .eq. 0.0d0 ) then
         iflag = 0
         maxits = 0
         return  
      end if

      if (stc.eq.'abs') then
         resid=resnor
      else
         resid=resnor/res0
      endif

      if ( resid .le. eps ) then
         iflag = 0
         maxits = 0
         return
      end if 
 
C     Main iterative loop ============================
      k = -1
      do while (.true.)

         if(oktest)write(*,199)its,resid
 199     format('   its =', i4, ' resid =', d20.6)

C        Loop to increment dimension of projection subspace
         k=k+1

C        Number of step (not restart) to be done
         its = its + 1
C        write(*,'(A,i3)') '+++++++++++++++++ its ',its 

C        - - - - - - - - - - - - - - - - - - - - - - - - -
C        This part should deliver 
C        u(1,k) <-- invA * resid
C        where u(1,k) is the k-th column of array u(1:n,0:m) and
C        invA is some reasonable approximation to the inverse of A
C
C        If mgmres=0 then no inner iterations are performed 
C        to get invA, so that invA is just the identity matrix. 
C        In this case algorithm GMRESR is nothing but GCR
C
C        Otherwise for inner iterations we perform ONE restart of GMRES
C        ATTN: for this implementation it is crucial to perform only
C        ONE restart of GMRES
         if (mgmres.eq.0) then
C           u(1,k) := resid  
            call dcopy(n,work(1,workres),1,work(1,u+mod(k,j)),1)
            call matvec_sys(work(1,u+mod(k,j)),work(1,c+mod(k,j)),n,
     1                                               amat,ibl7_amat)
            nits=nits+1
         else
C           Solve linear system A*u(1,k)=resid by GMRES
C           The stopping criterion for inner iterations is 
C           always absolute but it is automatically adjusted
C           not to be stricter than the stopping criterion for the 
C           outer iterations.  For example, if stop.criterion for
C           the outer iterations is relative than absolute stop.
C           criterion for the inner iterations is (eps*res0)
C           Accuracy for inner iteration:

            if(stc.eq.'abs')then
               epsinn = eps
            else
               epsinn = eps*res0
            endif

C           After envoking gmres0 epsinn and itsinn contain actual achieved
C           accuracy and number of performed iterations respectively

            itsinn=mgmres

            call gmres0(oktest,n,mgmres,
     &           work(1,workres),work(1,u+mod(k,j)),
     &           work(1,c+mod(k,j)),work(1,workgmre),
     &           epsinn,itsinn,amat,ibl7_amat)

            nits=nits+itsinn
         end if           
C - - - - - - - - - - - - - - - - - - - - - - - - 
      
C        Inner loop to orthogonalize 
C        c(1,k) with respect to c(1,k-j),...,c(1,k-1)
C        and to update correspondingly 
C        u(1,k) with respect to u(1,k-j),...,u(1,k-1)
C        parameter j is used only here
         do i = max0(0,k-j),k-1
            alphai = ddot(n,work(1,c+mod(i,j)),1,work(1,c+mod(k,j)),1)
            call daxpy(n,-alphai,work(1,c+mod(i,j)),1,
     &           work(1,c+mod(k,j)),1)
            call daxpy(n,-alphai,work(1,u+mod(i,j)),1,
     &           work(1,u+mod(k,j)),1)
         end do

C        Normalize c(1,k) and "normalize" u(1,k)
         cknorm = dnrm2(n,work(1,c+mod(k,j)),1)
         cknorm = 1 / cknorm
         call dscal(n,cknorm,work(1,c+mod(k,j)),1)
         call dscal(n,cknorm,work(1,u+mod(k,j)),1)

C        Update current solution and residual
         ckres = ddot(n,work(1,c+mod(k,j)),1,work(1,workres),1)
         call daxpy(n, ckres,work(1,u+mod(k,j)),1,x,          1)
         call daxpy(n,-ckres,work(1,c+mod(k,j)),1,work(1,workres),1)

C        call show(n,10,x,'GMRESR       ')  

C        Calculate residual norm, check convergence
         resnor = dnrm2(n,work(1,workres),1)

         if (stc.eq.'abs') then
            resid=resnor
         else
            resid=resnor/res0
         endif

         if ( resid .le. eps ) then
            iflag = 0
            maxits = nits
            return
         end if
         if (its .ge. maxits*j) then
            iflag = 1
            maxits = nits
            return
         end if

C        print 11, '            ||res|| = ',resnor 
C 11     format(A,d)
C 13     format(i4,A,d)
      
      end do
C End of inifinite iterative loop =================
C End of GMRESR subroutine      
      end 

C=============================================================================
       subroutine gmres0(oktest,n,im,rhs,uu,cc,work0,eps,maxits,amat,
     1 ibl7_amat)

C This is the modified GMRES routine gmres0 adapted for GMRESR by 
C Mike Botchev, Utrecht University, Dec. 1996
C For detail on how to make GMRES (for GMRESR) cheaper see 
C the above-mentioned paper on GMRESR 
c*************************************************************
C This code was initially written by Youcef Saad
C then revised by Henk A. van der Vorst  
C and Mike Botchev (oct. 1996)
C ************************************************************ 
c gmres algorithm . simple version .  (may 23, 1985)
c parameter list:
c oktest == TRUE for printing intermediate results
c n      == size of problem
c im     == size of krylov subspace:  should not exceed 50 in this
c          version (can be reset in code. looking at comment below)
c rhs    == right hand side
c uu     == initial guess for vector u (see above-mentioned paper on GMRESR)
c           on input, approximate solution on output
c cc     == initial guess for vector c (see above-mentioned paper on GMRESR)
c           on input, approximate solution on output
c work0  == work space of size n x (im+1)
c eps    == tolerance for stopping criterion. process is stopped
c           as soon as ( ||.|| is the euclidean norm):
c           || current residual || <= eps  
c maxits == maximum number of iterations allowed
c           on OUTPUT: actual number of iterations performed
c ----------------------------------------------------------------
c subroutines 
c matvec      == matrix vector multiplication y <- A*x
c
c BLAS:
c dcopy       == y <-- x routine
c ddot        == dot product function
c dnrm2       == euclidean norm function
c daxpy       == y <-- y+ax routine
c dscal       == x <-- ax routine
c dtsrv       == to solve linear system with a triangular matrix
c*************************************************************
c-------------------------------------------------------------
c arnoldi size should not exceed 10 in this version..
c to reset modify maxdim. BUT:             ----------------
c maxdim was set to 10 because of two reasons:
c (1) it is assumed in this implementation that it is cheaper to
c make maxdim vector updates than to make 1 matrix-vector
c multiplication;
c (2) for large maxdim we may lose the orthogonality property
c on which this cheap implementation is based.
c Please keep it in mind changing maxdim
c-------------------------------------------------------------
      integer maxdim,maxd1,md1max
      parameter (maxdim=10, maxd1=maxdim+1, md1max=maxdim*maxd1)

      logical oktest
      integer jjj,jj1,ibl7_amat
      integer i,i1,im,its,j,k,k1,maxits,n
      REAL cc,coeff,coef1,dabs,ddot,dnrm2,dsqrt,eps,epsmac,
     &                 gam,rhs(n),ro,uu(n),work0(n,im+1),t,amat(*)

      REAL hh(maxd1,maxdim),hh1(maxd1,maxdim),c(maxdim),
     &                 s(maxdim),rs(maxd1),rs1(maxd1)

      data (( hh(jj1,jjj), jj1=1,maxd1), jjj=1,maxdim) / md1max*0.0 / ,
     &      epsmac / 1.d-16 / 
C-----------------------------------------------------------------------------

      if (im .gt. maxdim) then
         im = maxdim
         write (*,'(A,i2)') 'GMRES0: dimension has been reduced to ',im
         write (*,'(A)') ' => reset MAXDIM if you want it to be more'
         write (*,'(A)') ' BUT read comments near MAXDIM before'
      end if

      its = 0

C     ----------------------------
C     Outer loop starts here.. 
C     BUT here (for GMRESR) only 1 outer loop is allowed
C     Compute initial residual vector 
C     ----------------------------
 10   continue
C        do not calculate initial residual first restart because 
C        initial guess is always zero. 
C        make initial guess zero:
         coeff = 0.0
         call dscal(n,coeff,uu,1)
C        make initial residual right-hand side:
         call dcopy(n,rhs,1,work0,1)

	 ro = dnrm2 (n, work0, 1)
	 if ((ro .eq. 0.0d0).or.(ro .le. eps)) then
            call matvec_sys(uu, cc, n,amat,ibl7_amat)
            eps = ro
            maxits = its 
            return
         end if

         coeff = 1 / ro
         call dscal(n,coeff,work0,1)

	 if (oktest) write(*, 199) its, ro

c        initialize 1-st term  of rhs of hessenberg system..
	 rs(1) = ro
	 i = 0

 4       continue
            i=i+1
            its = its + 1
            i1 = i + 1
            call matvec_sys(work0(1,i), work0(1,i1), n,amat,ibl7_amat)
c           -----------------------------------------
c           modified gram - schmidt...
c           -----------------------------------------
            do j=1, i
               t = ddot(n, work0(1,j),1,work0(1,i1),1)
               hh(j,i) = t
               call daxpy(n, -t, work0(1,j), 1, work0(1,i1), 1)
            end do
            t = dnrm2(n, work0(1,i1), 1)
            hh(i1,i) = t
            if (t .ne. 0.0d0)then
               t = 1 / t
               call dscal(n, t, work0(1,i1), 1)
C              save new column of hh in hh1 to reproduce vector cc later on
               call dcopy(maxd1,hh(1,i),1,hh1(1,i),1)
            endif
c           done with modified gram schmidt and arnoldi step..

c           now  update factorization of hh
            if (i .ne. 1) then
c              perform previous transformations  on i-th column of h
               do k=2,i
                  k1 = k-1
                  t = hh(k1,i)
                  hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
                  hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
               end do
            endif
            gam = dsqrt(hh(i,i)**2 + hh(i1,i)**2)
            if (gam .eq. 0.0d0) gam = epsmac

c           determine next plane rotation
            c(i) = hh(i,i)/gam
            s(i) = hh(i1,i)/gam
            rs(i1) = -s(i)*rs(i)
            rs(i) =  c(i)*rs(i)

c           determine residual norm and test for convergence-
            hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
            ro = dabs(rs(i1))
            if (oktest) write(*, 199) its, ro
	 if ((i .lt. im) .and. (ro .gt. eps))  goto 4
c
c        now compute solution. first solve upper triangular system.
c
C        rs := hh(1:i,1:i) ^-1 * rs

         call dtrsv('U','N','N',i,hh,maxd1,rs,1)   
c        done with back substitution..

c        now form linear combination to get vector uu
	 do j=1, i
            t = rs(j)
	    call daxpy(n, t, work0(1,j), 1, uu,1)
         end do
C        DO NOT restart outer loop EVEN when necessary (that is crucial
C        for this implementation of GMRESR):  NEVER goto 10 !  
C     if (ro .gt. eps .and. its .lt. maxits) goto 10

C     Finally, reproduce vector cc as cc = A*uu = work0*hh1*rs:
C     rs := hh1(1:i1,1:i) * rs
      coeff = 1
      coef1 = 0
      call dgemv('N',i1,i,coeff,hh1,maxd1,rs,1,coef1,rs1,1)

C     now multiply Krylov basis vectors work0 by rs:
C     cc := work0*rs
      call dscal(n,coef1,cc,1)
      do j=1, i1
         t = rs1(j)
	 call daxpy(n, t, work0(1,j), 1, cc,1)
      end do        

 199  format('itsinn =', i4, ' res. norm =', d20.6)

      maxits=its
      eps=ro 
      return
c------------------------------- end of gmres0 ----------------------
      end
