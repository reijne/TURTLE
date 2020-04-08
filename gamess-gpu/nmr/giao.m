c
c     NMR code from HONDO 2005
c     taken with kind permission of M. Dupuis
c
      subroutine giaox(ibase)
      implicit REAL (a-h,o-z)
      parameter   (mxioda=1000)
      parameter   (mxatom=MAXAT)
      parameter   (maxorb=MAXORB)
      character*8 errmsg
      character*8 wfntyp
      character*8 scftyp
      character*8 memtyp
      character*8 hndtyp
      character*8 dsktyp
      character*8 filtyp
      character*8 pcktyp
      character*8 direct
      character*8 semid 
      character*8 scf   
      character*8 rhf   
      logical     times
      logical     some
      logical     out
      logical     dbug    
      logical     disk
INCLUDE(../m4/common/dump3)
      common/maxlen/maxcor
      common/iofile/ir,iw,ipu(3),idaf
      common/hdafile/idafh,navh,ioda(2,mxioda)
      common/output/nprint
      common/prop3c/wfntyp
      common/scfopt1/scftyp
      common/qmtopt/tolqmt,nqmt
      common/intopt/nopk,nok,nosqur
      common/hndopt/hndtyp,dsktyp,filtyp,pcktyp
      common/intfil/nintmx
      common/pcklab/labsiz
INCLUDE(../m4/common/mapper)
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     1 no(10),nco,nseto,npair
      common/scfmem/ifck,idum,iden,jdum(3),ibuf,kdum(2),need
      common/dircas/i2case
      common/dirmem/ldum(2),jden,jfck,mdum(2),ndim
      common/dirgia/ifck0,ifckx,ifcky,ifckz,idens,idim
      common/cvggia/cnva,cnvb,damp
      common/pargia/vnew,vold,maxit
      common/filgia/ift1,ift2,ift3,ift4
INCLUDE(../m4/common/vcore)
chvd
c     it seems that we need to muck about with iadvec
c     and maxvec at this point to set stuff up for
c     subroutine tfsq.
c
      common/memvec/iadvec,maxvec
chvd
      dimension    b20(6)
      dimension    cm(3)
      dimension    errmsg(3)
      data errmsg  /'program ','stop in ','- giaox-'/
      data memtyp  /'in core '/
      data direct  /'direct  '/
      data semid   /'semi-dir'/
      data scf     /'scf     '/
      data rhf     /'rhf     '/
      data zero    /0.0d+00/
      data one     /1.0d+00/
      data two     /2.0d+00/
      data cvel    /137.0359895d+00/
      data ppm     /26.62566914d+00/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      some=.true.
      some=some.or.out
c
      times=.true. 
      times=times.and.some
c
      disk  = dsktyp.ne.memtyp
c
      if(some) then
         write(iw,9999)
      endif
      call fiflsh(iw)
c
      if(wfntyp.ne.scf) then 
         if(scftyp.ne.rhf) then
            if(nseto.ne.0) then
               write(iw,9998)
               call hnderr(3,errmsg)
            endif
         endif
      endif
c
      l1= num
      l2=(num*(num+1))/2
      l3= num*num
      l0= nqmt
      m1= nqmt
      m2=(nqmt*(nqmt+1))/2
      m3= nqmt* nqmt
      noc=na
c
      ift1=11
      ift2=12
      ift3=13
      ift4=14
      call rewfil(ift1)
      call rewfil(ift2)
      call rewfil(ift3)
      call rewfil(ift4)
c
c     ----- input -----
c
      call giainp
      call fiflsh(iw)
c
c     ----- set pointers for partitioning of core -----
c
      i00=ibase
      i01=i00+l3
      i02=i01+l1
      i03=i02+nat*3
      i04=i03+nat*9
      i05=i04+nat*9
      i06=i05+nat*9
      i07=i06+nat*9
c
c     ----- more memory -----
c
      i10=i07
      i20=i10+l2
      i30=i20+l2
      i40=i30+l2
      i50=i40+l2
      i60=i50+l3
      i70=i60+l3
      i80=i70+l3
      i90=i80+l3
      i100=i90+l3
      i110=i100+l1
      i120=i110+l1
      i130=i110+nintmx*(1+labsiz)
      lasti=i130-1
c
      j10=i10
      j11=j10+l3
      j12=j11+l3
      j20=j12+l3
      j21=j20+l3
      j22=j21+l3
      j30=j22+l3
      j40=j30+l3
      j50=j40+l2
      j60=j50+max(l1,nintmx*(1+labsiz))
      lastj=j60-1
c
      last =max(lasti,lastj)
chvd
c     muck with iadvec and maxvec
c
      iadvec = last
      maxvec = maxcor
chvd
      next =last+1
      need=l3+l1+l3+3*noc+last
      if(need.gt.maxcor) then
         write(iw,9994) need,maxcor
         call hnderr(3,errmsg)
      endif
c
c     ----- scf orbitals -----
c
      call rdedx(Q(i00),l3,ibl3qa,idaf)
      call rdedx(Q(i01),l1,ibl3ea,idaf)
      call tdown(Q(i00),ilifq,Q(i00),ilifq,l0)
      if(out) then
         write(iw,*) 'c(0)'
         call prev(Q(i00),Q(i01),l0,l1,l1)
      endif
c
      if(out) then
c
c        ----- the -dpl- cphf equations as a check -----
c
         write(iw,*) '-dpl- cphf equations as a check'
         call fiflsh(iw)
c
         do i=1,3
            cm(i)=zero
         enddo
         call dipint(Q(i10),Q(i20),Q(i30),cm,l1)
         call dawrit(idafh,ioda,Q(i10),l2,33,navh)
         call dawrit(idafh,ioda,Q(i20),l2,34,navh)
         call dawrit(idafh,ioda,Q(i30),l2,35,navh)
c
         mfld=1
         nfld=1
         do ifld=mfld,nfld
            if(out) then
               write(iw,9996)
               write(iw,*) 'field direction',ifld
            endif
c
c           ----- get perturbation -----
c
            if(ifld.eq.1) then
               idar=33
            elseif(ifld.eq.2) then
               idar=34
            elseif(ifld.eq.3) then
               idar=35
            endif
            call daread(idafh,ioda,Q(i10),l2,idar)
c
c           ----- initial u(1) -----  
c
            call giau_dpl(Q(i01),Q(i00),
     +                    Q(i10),Q(i60),
     +                    Q(i70),num,na,l0,l1)
c
            dmp=damp
            iter=-1
  100       iter=iter+1
            if(dbug) then
               write(iw,*) 'iteration',iter,maxit
            endif
c
c           ----- d(1) -----
c
            call giad_dpl(Q(i00),Q(i60),
     +                    Q(i70),Q(i20),
     +                    num,na,l0,l1)
c
c            ----- check convergence -----
c
            call giacvg_dpl(Q(i20),Q(i10),
     +                      Q(i60),Q(i70),
     +                      num,na,l0,l1,iend,iter,some)
            if((iend.ne.1).or.(iend.eq.1.and.dmp.ne.zero)) then      
c
c              ----- damping -----
c
               if(iter.eq.0) then
                  call dawrit(idafh,ioda,Q(i20),l2,31,navh)
               else
                  if(iend.eq.1) then
                     dmp=zero
                  endif
                  if(dmp.ne.zero) then
                     call daread(idafh,ioda,Q(i30),l2,31)
                     do i=1,l2
                        Q(i-1+i20)=dmp*Q(i-1+i30)
     +                            +(one-dmp)*Q(i-1+i20)
                     enddo
                     call dawrit(idafh,ioda,Q(i20),l2,31,navh)
                  else
                     call dawrit(idafh,ioda,Q(i20),l2,31,navh)
                  endif
               endif
c
c              ----- f(1) -----
c
               if(dbug) then
                  write(iw,*) 'd(1)'
                  call prtr(Q(i20),l1)
               endif
               call lcpsnc
               if(hndtyp.eq.direct.or.hndtyp.eq.semid) then
                  ifck=i30
                  iden=i20
                  ibuf=i40
chvd              call hdir(0)
               else
                  if(disk) then
chvd                 call hstar(Q(i20),Q(i30),
chvd +                          Q(i120),
chvd +                          Q(i120+nintmx),
chvd 1                          nintmx,iky,nopk)
                  else
chvd                 call xstar(Q(i20),Q(i30),
chvd +                          Q(i120),
chvd +                          Q(i120+nintmx),
chvd 1                          nintmx,iky,nopk)
                  endif
               endif
               call lcpadd(Q(i30),Q(i40),l2)
c
c              ----- add h(1) to f(1) -----
c
               do i=1,l2
                  Q(i-1+i30)=Q(i-1+i30)
     +                      -(-Q(i-1+i10))
               enddo
               if(dbug) then
                  write(iw,*) 'f(1)'
                  call prtr(Q(i30),l1)
               endif
c
c              ----- new u(1) -----
c
               call giau_dpl(Q(i01),Q(i00),
     +                       Q(i30),Q(i60),
     +                       Q(i70),num,na,l0,l1)
c
               if(dbug) then
                  write(iw,*) 'end of iteration',iter,maxit
               endif
               if(iter.gt.maxit) then      
                  write(iw,9997)   
                  call hnderr(3,errmsg)
               endif
               go to 100
            else 
c
c              ----- converged ; save u(1) and d(1) -----
c
               write(iw,*) '-dpl- cphf converged'
            endif
c
         enddo
         if(out) then
c
c           ----- for debugging get -f00- -----
c
            call daread(idafh,ioda,Q(i20),l2,16)
            if(dbug) then
               write(iw,*) 'd(0)'
               call prtr(Q(i20),l1)
            endif
            call lcpsnc
            if(hndtyp.eq.direct.or.hndtyp.eq.semid) then
               ifck=i30
               iden=i20
               ibuf=i40
chvd           call hdir(0)
            else
               if(disk) then
chvd              call hstar(Q(i20),Q(i30),
chvd +                       Q(i120),Q(i120+nintmx),
chvd 1                       nintmx,iky,nopk)
               else
chvd              call xstar(Q(i20),Q(i30),
chvd +                       Q(i120),Q(i120+nintmx),
chvd 1                       nintmx,iky,nopk)
               endif
            endif
            call lcpadd(Q(i30),Q(i40),l2)
            if(dbug) then
               write(iw,*) 'f(0)'
               call prtr(Q(i30),l1)
            endif
c
         endif
c
      endif ! out
      call fiflsh(iw)
c
c     ----- the -giao- cphf equations now -----
c
      if(out) then
         write(iw,*) '-giao- cphf equations now'
      endif
      call fiflsh(iw)
      do i=1,3*nat
         Q(i+i02-1)=zero
      enddo
      do i=1,9*nat
         Q(i+i03-1)=zero
         Q(i+i04-1)=zero
         Q(i+i05-1)=zero
      enddo
c
c     ----- calculate -giao- perturbation operators -----
c           -s10-     on     -ift1-
c           -h10-     on     -ift2-
c           -h01-     on     -ift3-
c
      if(out) then
         write(iw,*) '-h20- diamagnetic susceptibility'
         call fiflsh(iw)
      endif
      call cpuwal(tim,wtim)
       tim0= tim
      wtim0=wtim
      call giah20(Q(j10),Q(i00),
     +            b20,Q(j40),l1,noc,l1)   
      call cpuwal(tim,wtim)
       tim1= tim
      wtim1=wtim
       tim20= tim1- tim0
      wtim20=wtim1-wtim0
      write(iw,9993)
      write(iw,9992) b20
      if(times) then
         write(iw,9991) tim20,wtim20
      endif
      call fiflsh(iw)
c
      if(out) then
         write(iw,*) '-h01-'
         call fiflsh(iw)
      endif
      call cpuwal(tim,wtim)
       tim0= tim
      wtim0=wtim
      call giah01(Q(j10),Q(i00),
     +            Q(j30),Q(j40),
     +            l1,noc,l0,l1)
      call cpuwal(tim,wtim)
       tim1= tim
      wtim1=wtim
       tim01= tim1- tim0
      wtim01=wtim1-wtim0
      if(times) then
         write(iw,9990) tim01,wtim01
      endif
      call fiflsh(iw)
c
      if(out) then
         write(iw,*) '-s10-'
         call fiflsh(iw)
      endif
      call cpuwal(tim,wtim)
       tim0= tim
      wtim0=wtim
      call gias10(Q(j10),Q(i00),
     +            Q(j30),Q(j40),
     +            l1,noc,l0,l1)   
      call cpuwal(tim,wtim)
       tim1= tim
      wtim1=wtim
       tim10= tim1- tim0
      wtim10=wtim1-wtim0
      if(times) then
         write(iw,9989) tim10,wtim10
      endif
      call fiflsh(iw)
c
      if(out) then
         write(iw,*) '-h10-'
         call fiflsh(iw)
      endif
      ifck0=j30
      ifckx=j20
      ifcky=j21
      ifckz=j22
      idens=j40
      idim =l1
      call cpuwal(tim,wtim)
       tim0= tim
      wtim0=wtim
      call giah10(Q(j10),Q(i00),
     +            Q(j30),Q(j40),
     +            Q(j50),Q(j30),
     +            l1,noc,l0,l1,next)   
      call cpuwal(tim,wtim)
       tim1= tim
      wtim1=wtim
       tim10= tim1- tim0
      wtim10=wtim1-wtim0
      if(times) then
         write(iw,9988) tim10,wtim10
      endif
      call fiflsh(iw)
c
c     ----- ready for the -cphf- iterations -----
c
      call cpuwal(tim,wtim)
       tim0= tim
      wtim0=wtim
c
      mfld=1
      nfld=3
      do ifld=mfld,nfld
         if(some) then
            write(iw,9996)
            write(iw,*) 'field direction',ifld
            call fiflsh(iw)
         endif
c
c     ----- initial -u(1)- from -h(1)- and -s(1)- -----
c           -u(1)-  .....  record -32- of -idafh-
c           -d(1)-  .....  record -31- of -idafh-
c
         do i=1,l3
            Q(i+j10-1)=zero
         enddo
         call giafil(ifld)
         call cpuwal(tim,wtim)
          tim1= tim
         call giau1(Q(i01),Q(j30),
     +              Q(j11),Q(j10),
     +              Q(j20),num,na,l0,l1)
         call cpuwal(tim,wtim)
         if(times) then
            write(iw,*) 'time -u1   - = ',(tim-tim1)
         endif
c
         dmp=damp
         iter=-1
  200    iter=iter+1
         if(dbug) then
            write(iw,*) 'iteration',iter,maxit
            call fiflsh(iw)
         endif
c
c     ----- d(1) -----
c
         call cpuwal(tim,wtim)
          tim1= tim
         call giav1(Q(i00),Q(j30),
     +              Q(j10),num,na,l0,l1)
         call giad1(Q(i00),Q(j10),
     +              Q(j20),num,na,l1)
         call cpuwal(tim,wtim)
         if(times) then
            write(iw,*) 'time -v1+d1- = ',(tim-tim1)
         endif
c
c      ----- check convergence -----
c
         call cpuwal(tim,wtim)
          tim1= tim
         call giacvg(Q(j20),Q(j21),
     +               Q(j30),Q(j21),
     +               num,na,l0,l1,iend,iter,some)
         call cpuwal(tim,wtim)
         if(times) then
            write(iw,*) 'time - cvg - = ',(tim-tim1)
         endif
         if((iend.ne.1).or.(iend.eq.1.and.dmp.ne.zero)) then
c
c     ----- damping -----
c
            if(iter.eq.0) then
               call dawrit(idafh,ioda,Q(j20),l3,31,navh)
            else
               if(iend.eq.1) then
                  dmp=zero
               endif
               if(dmp.ne.zero) then
                  call daread(idafh,ioda,Q(j21),l3,31)
                  do i=1,l3
                     Q(i-1+j20)=dmp*Q(i-1+j21)
     +                         +(one-dmp)*Q(i-1+j20)
                  enddo
                  call dawrit(idafh,ioda,Q(j20),l3,31,navh)
               else
                  call dawrit(idafh,ioda,Q(j20),l3,31,navh)
               endif
            endif
c
c     ----- k(1) -----
c
            if(dbug) then
               write(iw,*) 'd(1)'
               call prsq(Q(j20),l1,l1,l1)
            endif
c
            call lcpsnc
            if(hndtyp.eq.direct.or.hndtyp.eq.semid) then
               jfck=j10
               jden=j20
               ndim=l1
               i2case=6
               write(iw,*) 'not direct yet. stop'
               call hnderr(3,errmsg)
            else
               if(disk) then
                  call cpuwal(tim,wtim)
                   tim1= tim
                  call giak10(Q(j20),Q(j11),
     +                        Q(j50),Q(j50+nintmx),
     +                        nintmx,l1,nopk)
                  call cpuwal(tim,wtim)
                  if(times) then
                     write(iw,*) 'time - k10 - = ',(tim-tim1)
                  endif
               else
                  write(iw,*) 'not memory yet. stop'
                  call hnderr(3,errmsg)
               endif
            endif
            call lcpadd(Q(j11),Q(j10),l3)
c
            if(dbug) then
               write(iw,*) 'k10(1) in ao basis'
               call prsq(Q(j11),l1,l1,l1)
            endif
c
            call cpuwal(tim,wtim)
             tim1= tim
            call tfsq(Q(j10),Q(j11),
     +                Q(i00),Q(j40),
     +                l0,l1,l1)
            call cpuwal(tim,wtim)
            if(times) then
               write(iw,*) 'time -tfsq - = ',(tim-tim1)
            endif
            if(dbug) then
               write(iw,*) 'k10(1) in canonical mo basis'
               call prsq(Q(j10),l1,l1,l1)
            endif
c
c     ----- new u(1) -----
c
            call giafil(ifld)
            call cpuwal(tim,wtim)
             tim1= tim
            call giau1(Q(i01),Q(j30),
     +                 Q(j11),Q(j10),
     +                 Q(j20),
     +                 num,na,l0,l1)
            call cpuwal(tim,wtim)
            if(times) then
               write(iw,*) 'time -u1   - = ',(tim-tim1)
            endif
c
            if(dbug) then
               write(iw,*) 'end of iteration',iter,maxit
            endif
            if(iter.gt.maxit) then
               write(iw,9997)
               call hnderr(3,errmsg)
            endif
            go to 200
         else
c
c     ----- converged ; save u(1) and d(1) -----
c           -d(1)-    to    -31-
c           -u(1)-    to    -32-
c
            if(some) then
               write(iw,*) '-giao- cphf converged'
            endif
            call dawrit(idafh,ioda,Q(j20),l3,31,navh)
            call dawrit(idafh,ioda,Q(j30),l3,32,navh)
            if(out) then
               write(iw,*) 'converged -u1-'
               call prsq(Q(j30),l0,l0,l1)
               write(iw,*) 'converged -d1-'
               call prsq(Q(j20),l1,l1,l1)
            endif
         endif
c
c     ----- -h01- contribution to magnetic shielding -----
c           -d(1)-    from  -31-
c           -d(0)-    from  -16-
c
         call rdedx(Q(j40),l2,ibl3pa,idaf)
         call daread(idafh,ioda,Q(j30),l3,31)
         call cpuwal(tim,wtim)
          tim1= tim
         call giae11(Q(i02),Q(i03),
     +               Q(j10),Q(j30),
     +               Q(j40),l1,ifld)
         call cpuwal(tim,wtim)
         if(times) then
            write(iw,*) 'time - e11 - = ',(tim-tim1)
         endif
c
      enddo
      call cpuwal(tim,wtim)
       tim1= tim
      wtim1=wtim
       tim10= tim1- tim0
      wtim10=wtim1-wtim0
      if(times) then
         write(iw,9987) tim10,wtim10
      endif
c
c     ----- global sum on -e11- -----
c
      call lcpadd(Q(i03),Q(i06),9*nat)
c
c     ----- -h11- contributions to magnetic shielding -----
c
      call rdedx(Q(j40),l2,ibl3pa,idaf)
      call giah11(Q(i04),Q(i05),Q(j40))
c
      units =ppm/two    
      do i=1,9*nat
         Q(i+i03-1)= Q(i+i03-1)*units
         Q(i+i04-1)= Q(i+i04-1)*units
         Q(i+i05-1)= Q(i+i05-1)*units
      enddo
      call giaprt(Q(i03),nat,1)
      if(out) then
         call giaprt(Q(i04),nat,2)
         call giaprt(Q(i05),nat,3)
      endif
c
      do i=1,9*nat
         Q(i+i04-1)=(Q(i+i04-1)+Q(i+i05-1))
      enddo
      call giaprt(Q(i04),nat,4)
c
      do i=1,9*nat
         Q(i+i03-1)=(Q(i+i03-1)+Q(i+i04-1))
      enddo
      call giaprt(Q(i03),nat,5)
c
c     ----- all done -----
c
      write(iw,9996)
      call fiflsh(iw)
c
      return
 9999 format(
     1 10x,57('-'),/,
     2 10x,'-giao- magnetic susceptibility/chemical shielding tensors',/,
     3 10x,57('-'))
 9998 format(' -giao- implemented for closed shell systems only. stop')
 9997 format(' too many iterations in the -giao- cphf . stop. ')
 9996 format(/)
 9995 format(i7,3f15.8)
 9994 format(' not enough memory in -giaox- . stop.',/,
     1       ' need, maxcor = ',2i10)
 9993 format(/,' diamagnetic susceptibity tensor -b20- (a.u.) = ',/,
     1         ' ---------------------------------------------- ',/)
 9992 format(' b20(xx) = ',f20.10,/,
     1       ' b20(yy) = ',f20.10,/,
     2       ' b20(zz) = ',f20.10,/,
     3       ' b20(xy) = ',f20.10,/,
     4       ' b20(xz) = ',f20.10,/,
     5       ' b20(yz) = ',f20.10,/)
 9991 format(' -giah20- , cpu time = ',f10.2,' wall time = ',f10.2)
 9990 format(' -giah01- , cpu time = ',f10.2,' wall time = ',f10.2)
 9989 format(' -gias10- , cpu time = ',f10.2,' wall time = ',f10.2)
 9988 format(' -giah10- , cpu time = ',f10.2,' wall time = ',f10.2)
 9987 format(' -giak10- , cpu time = ',f10.2,' wall time = ',f10.2)
      end
      subroutine giaprt(cs,nat,iflg)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxbfn =MAXORB)
      parameter (maxorb=MAXORB)
      parameter (two   =2.0d+00)
      parameter (three =3.0d+00)
      common/iofile/ir,iw
      common/mollab/anam(mxatom),bnam(mxatom),bflab(mxbfn)
INCLUDE(../m4/common/mapper)
      dimension cs(3,3,*)
      dimension a(6),axs(3,3),eig(3)
c
      if(iflg.eq.1) then
c
c     ----- paramagnetic -----
c
         write(iw,9995) 
      elseif(iflg.eq.2) then
c
c     ----- diamagnetic -----
c
         write(iw,9996) 
      elseif(iflg.eq.3) then
c
c     ----- paramagnetic -----
c
         write(iw,9997) 
      elseif(iflg.eq.4) then
c
c     -----              -----
c
         write(iw,9998) 
      elseif(iflg.eq.5) then
c
c     ----- paramagnetic + diamagnetic -----
c
         write(iw,9999) 
      endif
c
c     ----- loop over atoms -----
c
      do iat=1,nat
         write(iw,9989) iat,anam(iat) 
c
         write(iw,9994)
         write(iw,9993) cs(1,1,iat),cs(1,2,iat),cs(1,3,iat)
         write(iw,9992) cs(2,1,iat),cs(2,2,iat),cs(2,3,iat)
         write(iw,9991) cs(3,1,iat),cs(3,2,iat),cs(3,3,iat)
c
         if(iflg.ne.5) then
            averag=(cs(1,1,iat)+cs(2,2,iat)+cs(3,3,iat))/three
            write(iw,9988) averag
         endif
c
         if(iflg.eq.5) then
            a(1)= cs(1,1,iat)
            a(2)=(cs(2,1,iat)+cs(1,2,iat))/two
            a(3)= cs(2,2,iat)
            a(4)=(cs(3,1,iat)+cs(1,3,iat))/two
            a(5)=(cs(3,2,iat)+cs(2,3,iat))/two
            a(6)= cs(3,3,iat)
            call diaaxs(a,axs,eig,iky,3,3,3)
            averag=(eig(1) + eig(2) + eig(3))/three
            aniso = eig(1) -(eig(2) + eig(3))/two
            write(iw,9987) averag,aniso
            write(iw,9986) (    j ,j=1,3)
            write(iw,9985) (eig(j),j=1,3)
            write(iw,9984)
            do i=1,3
               write(iw,9983) i,(axs(i,j),j=1,3)
            enddo
         endif
c
      enddo
c
      return
 9999 format(//,' shielding tensor (ppm) = ',
     1          ' ( paramagnetic + diamagnetic ) ',
     2        /,' ------------------------ ')
 9998 format(//,' shielding tensor (ppm) = ',
     1          ' (              p00 * h11 ) ',
     2        /,' ------------------------ ')
 9997 format(//,' shielding tensor (ppm) = ',
     1          ' ( paramagnetic p00 * h11 ) ',
     2        /,' ------------------------ ')
 9996 format(//,' shielding tensor (ppm) = ',
     1          ' (  diamagnetic p00 * h11 ) ',
     2        /,' ------------------------ ')
 9995 format(//,' shielding tensor (ppm) = ',
     1          ' ( paramagnetic p10 * h01 ) ',
     2        /,' ------------------------ ')
 9994 format(/,7x,7x,'bx',3x,7x,'by',3x,7x,'bz')
 9993 format('      x',3f12.4)
 9992 format('      y',3f12.4)
 9991 format('      z',3f12.4)
 9990 format(/)
 9989 format(/,1x,i3,2x,a8,/,1x,'-------------')               
 9988 format(12x,' isotropic = ',f12.4)
 9987 format(12x,' isotropic = ',f12.4,/,
     1       12x,'anisotropy = ',f12.4)
 9986 format(/,7x,3(7x,i1,4x))
 9985 format(' eig = ',3f12.4)
 9984 format(' vec = ') 
 9983 format(6x,i1,3f12.4)
      end
      subroutine giainp
      implicit REAL (a-h,o-z)
c
c     ----- this routine reads the -$giao- namelist information -----
c           and initializes variables controlling the cphf run.
c
      common/iofile/ir,iw
      common/cvggia/cnva,cnvb,dampd
      common/pargia/vnew,vold,itmax
      data cnv1  /1.0d-03/
      data cnv2  /1.0d-03/
      data damp  /0.2d+00/
      data maxit /50/
c
      namelist /giao/ maxit,cnv1,cnv2,damp
c
cDEBUG
      write(*,*)'*** giainp: needs fixing'
      goto 10
cDEBUG
      call rewfil(ir)
      read(ir,giao,end=10,err=10)
   10 continue
c
      itmax =maxit
      cnva  =cnv1
      cnvb  =cnv2
      dampd =damp
c
      write(iw,9999) itmax,cnva,cnvb,dampd
c
      return
 9999 format(/,' ----- $giao -----',/,
     1       /,' maxit  = ',i10  ,/,' cnva   = ',f10.5,
     2       /,' cnvb   = ',f10.5,/,' damp   = ',f10.5)
      end
      subroutine giacvg_dpl(d,h,u,u0,num,noc,norb,ndim,
     1                      iend,iter,some)
      implicit REAL (a-h,o-z)
c
c     ----- this routine checks for convergence of the cphf -----
c           calculation.  convergence is determined by two
c           criteria: the convergence of the tensor element
c           for the xx, yy, or zz direction and the
c           convergence of the individual elements of the u
c           matrix.
c
      logical some
      parameter (zero=0.0d+00)
      parameter (mxioda=1000)   
      common/cvggia/cnva,cnvb
      common/pargia/vnew,vold
      common/iofile/ir,iw
      common/hdafile/idafh,navh,ioda(2,mxioda)
      dimension d(*),h(*),u(ndim,*),u0(ndim,*)
c
      if(iter.gt.0) then    
c
c     ----- iter > 1 -----
c
         vold  = vnew
         vnew  =-prpdot(d,h,num)
         diff  = vold - vnew
         udelm = zero
         uval  = zero
         ival  = 0
         jval  = 0
c
         num3=num*num
         call daread(idafh,ioda,u0,num3,32)
c
         icon=0
         do j=noc+1,norb
            do i=1,noc
               udel = u(i,j) - u0(i,j)
               if(abs(udel).ge.cnvb) icon = icon + 1
               if(abs(udel).gt.abs(udelm)) then
                  udelm = udel
                  uval = u(i,j)
                  ival = i
                  jval = j
               endif
            enddo
         enddo
         num3=num*num
         call dawrit(idafh,ioda,u,num3,32,navh)
c
         if(abs(diff).lt.cnva.and.icon.eq.0) then 
            iend = 1
         else
            iend = 0
         endif
         if(some) then
            write(iw,9999) iter,vnew,diff,udelm,uval,ival,jval,
     1                     icon,iend
         endif
         return
c
      else
c
c     ----- iter = 1 -----
c
         vold  = zero
         vnew  =-prpdot(d,h,num)
         diff  = vold - vnew
         udelm = zero
         uval  = zero
         ival  = 0
         jval  = 0
c
         num3=num*num
         call dawrit(idafh,ioda,u,num3,32,navh)
c
         icon = noc * (norb - noc)
         iend=0
         if(some) then
            write(iw,9999) iter,vnew,diff,udelm,uval,ival,jval,
     1                     icon,iend
         endif
      endif
      return
 9999 format(1x,i3,f20.9,3f16.9,2i4,i10,i3)
      end
      subroutine giad_dpl(c,u,x,d,num,noc,norb,ndim)
      implicit REAL (a-h,o-z)
c
c     ----- construct the derivative density matrix: -----
c                  d' = 2*((c'c+) + (cc'+))
c
      parameter (zero=0.0d+00)
      parameter (one =1.0d+00)
      parameter (two =2.0d+00)
      dimension c(ndim,*),u(ndim,*),x(ndim,*),d(*)
c
      num2=(num*(num+1))/2
      do i=1,num2
         d(i)=zero
      enddo
      do j=1,norb
         do i=1,num
            dum=zero
            do k=1,norb
               dum=dum+c(i,k)*u(k,j)
            enddo
            x(i,j)=dum
         enddo
      enddo
      do i=1,noc
         ic=1
         do j=1,num
            do k=1,j
               d(ic) = d(ic)+c(j,i)*x(k,i)+x(j,i)*c(k,i)
               ic=ic+1
            enddo
         enddo
      enddo
      do i=1,num2
         d(i)=d(i)*two
      enddo
      return
      end
      subroutine giau_dpl(e,c,f,u,fc,num,noc,norb,ndim)
      implicit REAL (a-h,o-z)
c
c     ----- construct the transformation matrix -u- used -----
c           to get the first order wavefunction by
c           projecting out the 'ij' element of the
c           derivative fock matrix and weighting it
c           by the energy denominator.
c
      parameter (zero=0.0d+00)
      logical   dbug
      common/iofile/ir,iw
      dimension e(*),c(ndim,*),f(*),u(ndim,*),fc(ndim,*)
c
      dbug=.false.
c
      do j=1,norb
         do is=1,num
            fc(is,j)=zero
         enddo
      enddo
      do j=noc+1,norb
         ist = 0
         do is=1,num
            csj = c(is,j)
            yy = zero
            do it=1,is
               ist = ist+1
               fc(it,j) = fc(it,j) + f(ist)*csj
               xx = f(ist)*c(it,j)
               yy = yy + xx
            enddo
            fc(is,j) = fc(is,j) + yy - xx
         enddo
      enddo
c
c     ----- f(1) .... if debug -----
c
      if(dbug) then
         do j=1,norb
            do i=1,norb
               u(i,j)=zero
            enddo
         enddo
         do i=1,noc
            do j=noc+1,norb
               yy = zero
               do is=1,num
                  yy = yy + c(is,i)*fc(is,j)
               enddo
               u(i,j) =  yy
               u(j,i) =  yy
            enddo
         enddo
         if(dbug) then
            write(iw,*) 'f(1) in -giau-'
            call prsq(u,norb,norb,ndim)
         endif
      endif
c
c     ----- u(1) -----
c
      do j=1,norb
         do i=1,norb
            u(i,j)=zero
         enddo
      enddo
      do i=1,noc
         do j=noc+1,norb
            yy = zero
            do is=1,num
               yy = yy + c(is,i)*fc(is,j)
            enddo
            diff = e(i) - e(j)
            yy = yy / diff
            u(i,j) = -yy
            u(j,i) =  yy
         enddo
      enddo
      if(dbug) then
         write(iw,*) 'u(1) in -giau-'
         call prsq(u,norb,norb,ndim)
      endif
c
      return
      end
      subroutine giau1(e,u,s,f,h,num,noc,norb,ndim)
      implicit REAL (a-h,o-z)
c
c     ----- construct the transformation matrix -u- used -----
c           to get the first order wavefunction by
c           projecting out the 'ij' element of the
c           derivative fock matrix and weighting it
c           by the energy denominator.
c
      parameter (zero=0.0d+00)
      logical   dbug
      logical   out
      common/iofile/ir,iw
      common/filgia/ift1,ift2,ift3,ift4
      dimension e(*)
      dimension f(ndim,*),h(ndim,*)
      dimension s(ndim,*),u(ndim,*)
      data two /2.0d+00/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
c     ----- u(1) -----
c
      do j=1,norb
         do k=1,norb
            u(k,j)=zero
         enddo
      enddo
c
c     ----- read overlap derivatives -----
c
      call giard(ift1,s,ndim*ndim)
c
c     -----    overlap      term for the ( occ,occ) block -----
c
      do j=1,noc
         do k=1,noc
            u(k,j)= - s(k,j) / two 
         enddo
      enddo
      if(dbug) then
         write(iw,*) 'u(1) ( occ,occ) block in -giau1-'
         call prsq(u,norb,norb,ndim)
         call fiflsh(iw)
      endif
      if(out) then
         write(iw,*) '-eig- in -giau1-'
         call prsq(e,1,norb,norb)
      endif
c
c     ----- h(1) -----
c
      call giard(ift2,h,ndim*ndim)
c
c     ----- u contains the overlap derivative term -----
c
      do j=1,noc
         do k=noc+1,norb
            dum    = (   e(j)   - e(k)   )
            dum1   = (   h(k,j) - f(k,j) ) / dum
            dum2   = ( - s(k,j) * e(j)   ) / dum
            u(k,j) =     dum1   + dum2
            if(dbug) then
               write(iw,9999) k,j,dum
               write(iw,9998) h(k,j),f(k,j),dum1
               write(iw,9997) s(k,j),e(j)  ,dum2,u(k,j)
            endif
         enddo
      enddo
c
      if(out) then
         write(iw,*) 'u(1) ( occ,occ) block in -giau1-'
         write(iw,*) 'u(1) (virt,occ) block in -giau1-'
         call prsq(u,norb,norb,ndim)
         call fiflsh(iw)
      endif
c
      return
 9999 format(' in -giau1- k,j = ',2i5,' e(j)-e(k) = ',f15.8)
 9998 format(' h(k,j),f(k,j),dum1 = ',3f15.8)
 9997 format(' s(k,j),e(j)  ,dum2 = ',3f15.8,/,
     1       '               dum  = ',30x,f15.8)
      end
      subroutine giav1(v0,u,v1,num,noc,norb,ndim)
      implicit REAL (a-h,o-z)
      parameter (zero=0.0d+00)
      logical   out
      common/iofile/ir,iw
      dimension v0(ndim,*),u(ndim,*),v1(ndim,*)
c
      out=.false.
c
      do j=1,noc
         do i=1,num
            dum=zero
            do k=1,norb
               dum=dum+v0(i,k)*u(k,j)
            enddo
            v1(i,j)=dum
         enddo
      enddo
c
      if(out) then
         write(iw,*) '-v(1)- in -giav1-'
         call prsq(v1,norb,num,ndim)
      endif
      return
      end
      subroutine giad1(v0,v1,d1,num,noc,ndim)
      implicit REAL (a-h,o-z)
c
c     ----- construct the derivative density matrix -----
c                  d' = 2*((v'v+) - (vv'+))
c
      parameter (zero=0.0d+00)
      parameter (two =2.0d+00)
      logical   dbug
      common/iofile/ir,iw
      dimension v0(ndim,*),v1(ndim,*),d1(ndim,*)
c
      dbug=.false.
c
      do j=1,num
         do i=1,num
            d1(i,j)=zero
         enddo
      enddo
      do ioc=1,noc
         do j=1,num
            do i=1,num
               d1(i,j) = d1(i,j)
     1                 + (v0(i,ioc)*v1(j,ioc)-v1(i,ioc)*v0(j,ioc))
            enddo
         enddo
      enddo
      do j=1,num
         do i=1,num
            d1(i,j)=d1(i,j)*two
         enddo
      enddo
c
      if(dbug) then
         write(iw,*) 'd(1) in -giad1-'
         call prsq(d1,num,num,ndim)
      endif
      return
      end
      subroutine giawt(ift,h,ndim)
      implicit REAL (a-h,o-z)
      dimension h(ndim)
      write(ift) h
      return
      end
      subroutine giard(ift,h,ndim)
      implicit REAL (a-h,o-z)
      dimension h(ndim)
      read(ift) h
      return
      end
      subroutine giafil(ifld)
      implicit REAL (a-h,o-z)
c
c     ----- position the -gia- files -----
c
      logical dbug
      common/iofile/ir,iw
      common/filgia/ift1,ift2,ift3,ift4
c
      dbug=.false.
c
      rewind(ift1)
      go to (30,20,10),ifld
   10 read(ift1)
   20 read(ift1)
   30 continue
      if(dbug) then
         write(iw,*) '-ift1- positioned to -ifld- = ',ifld
         call fiflsh(iw)
      endif
c
      rewind(ift2)
      go to (130,120,110),ifld
  110 read(ift2)
  120 read(ift2)
  130 continue
      if(dbug) then
         write(iw,*) '-ift2- positioned to -ifld- = ',ifld
         call fiflsh(iw)
      endif
c
      return
      end
      subroutine giak00(gauge,f10,f00,w,t,num,noc,norb,ndim)
      implicit REAL (a-h,o-z)
      parameter (maxorb=MAXORB)
      logical dbug
      logical out
      common/iofile/ir,iw
      common/filgia/ift1,ift2,ift3,ift4
INCLUDE(../m4/common/mapper)
      dimension f10(ndim,ndim,*)
      dimension f00(ndim,ndim,*)
      dimension   w(ndim,ndim)
      dimension   t(*)
      dimension gauge(3,noc)
      data zero    /0.0d+00/
      data one     /1.0d+00/
      data two     /2.0d+00/
      data cvel    /137.0359895d+00/
c
      dbug=.false.
      out =.true. 
      out =out.or.dbug
c
c     ----- gauge dependent terms first -----
c                  kx , ky , kz 
c
      if(out) then
         write(iw,*) '-k00- gauge dependent terms'      
      endif
      call rewfil(ift2)
      do iop=1,3
         if(out) then
            write(iw,*) 'read     rx, ry, rz ... = ',iop   
         endif
         norb2=(norb*(norb+1))/2
         call giard(ift2,t,norb2)
         do iorb=1,norb
            do jorb=1,iorb
               ij=iky(iorb)+jorb
               f10(iorb,jorb,iop)=t(ij)
               f10(jorb,iorb,iop)=t(ij)
            enddo
         enddo
      enddo
c
      iexch=1
      if(out) then
         write(iw,*) 'read k0             ... = ',iexch
      endif
      call rewfil(ift3)
      norb2=(norb*(norb+1))/2
      call giard(ift3,t,norb2)
      do iorb=1,norb
         do jorb=1,iorb
            ij=iky(iorb)+jorb
              w(iorb,jorb)=t(ij)
              w(jorb,iorb)=t(ij)
         enddo
      enddo
c
      call mxm(w(1,1),ndim,f10(1,1,1),ndim,f00(1,1,1),ndim,
     1        norb,norb,norb)
      call mxm(w(1,1),ndim,f10(1,1,2),ndim,f00(1,1,2),ndim,
     1        norb,norb,norb)
      call mxm(w(1,1),ndim,f10(1,1,3),ndim,f00(1,1,3),ndim,
     1        norb,norb,norb)
c
c     ----- we are ready ... -----
c
      call rewfil(ift1)
      call rewfil(ift4)
      do ioc=1,noc
         gx=gauge(1,ioc)
         gy=gauge(2,ioc)
         gz=gauge(3,ioc)
         if(out) then
            write(iw,9999) ioc,gx,gy,gz
         endif
         call giard(ift1,f10(1,1,1),ndim*ndim)
         call giard(ift1,f10(1,1,2),ndim*ndim)
         call giard(ift1,f10(1,1,3),ndim*ndim)
         if(dbug) then
            write(iw,*) 'b(x) perturbation before -k00-'
            call prsq(f10(1,1,1),norb,norb,ndim)
            write(iw,*) 'b(y) perturbation before -k00-'
            call prsq(f10(1,1,2),norb,norb,ndim)
            write(iw,*) 'b(z) perturbation before -k00-'
            call prsq(f10(1,1,3),norb,norb,ndim)
         endif
         do ifld=1,3
            if(ifld.eq.1) then
               c1=gz
               c2=gy
               iop1=2
               iop2=3
            elseif(ifld.eq.2) then
               c1=gx
               c2=gz
               iop1=3
               iop2=1
            elseif(ifld.eq.3) then
               c1=gy
               c2=gx
               iop1=1
               iop2=3
            endif 
            velfac=one/(two*cvel)
            do jorb=1,norb
               do iorb=1,norb
                  f10(iorb,jorb,ifld)=f10(iorb,jorb,ifld)
     1      -velfac*( c1*(f00(iorb,jorb,iop1)-f00(jorb,iorb,iop1))
     2               -c2*(f00(iorb,jorb,iop2)-f00(jorb,iorb,iop2)) )
               enddo
            enddo
         enddo
         call giawt(ift4,f10(1,1,1),ndim*ndim)
         call giawt(ift4,f10(1,1,2),ndim*ndim)
         call giawt(ift4,f10(1,1,3),ndim*ndim)
         if(dbug) then
            write(iw,*) 'b(x) perturbation with -k00- gauge dependent'
            call prsq(f10(1,1,1),norb,norb,ndim)
            write(iw,*) 'b(y) perturbation with -k00- gauge dependent'
            call prsq(f10(1,1,2),norb,norb,ndim)
            write(iw,*) 'b(z) perturbation with -k00- gauge dependent'
            call prsq(f10(1,1,3),norb,norb,ndim)
         endif
      enddo
      call rewfil(ift1)
      call rewfil(ift4)
c
c     ----- gauge independent terms now -----
c
      if(out) then
         write(iw,*) '-k00- gauge independent terms'      
      endif
      call rewfil(ift3)
      read(ift3)
      do iop=2,4
         if(out) then
            write(iw,*) 'read     kx, ky, kz ... = ',iop  
         endif
         norb2=(norb*(norb+1))/2
         call giard(ift3,t,norb2)
         do iorb=1,norb
            do jorb=1,iorb
               ij=iky(iorb)+jorb
               f00(iorb,jorb,iop)=t(ij)
               f00(jorb,iorb,iop)=t(ij)
            enddo
         enddo
      enddo
c
c     -----  x ,  y ,  z -----
c
      call rewfil(ift2)
      do iop=1,3
         if(out) then
            write(iw,*) 'read     rx, ry, rz ... = ',iop   
         endif
         norb2=(norb*(norb+1))/2
         call giard(ift2,t,norb2)
         do iorb=1,norb
            do jorb=1,iorb
               ij=iky(iorb)+jorb
               f10(iorb,jorb,iop)=t(ij)
               f10(jorb,iorb,iop)=t(ij)
            enddo
         enddo
      enddo
c
c     ----- we are ready ... -----
c
      call rewfil(ift1)
      call rewfil(ift4)
      do ioc=1,noc
         gx=gauge(1,ioc)
         gy=gauge(2,ioc)
         gz=gauge(3,ioc)
         if(out) then
            write(iw,9999) ioc,gx,gy,gz
         endif
         do ifld=1,3
            if(ifld.eq.1) then
               iop1=3
               iop2=2
            elseif(ifld.eq.2) then
               iop1=1
               iop2=3
            elseif(ifld.eq.3) then
               iop1=2
               iop2=1
            endif
            call giard(ift4,w(1,1),ndim*ndim)
            do jorb=1,norb
               do iorb=1,norb
                  dum=zero
                  do lorb=1,norb
                     dum=dum+f00(iorb,lorb,iop1)*f10(lorb,jorb,iop2)
     1                      -f00(jorb,lorb,iop1)*f10(lorb,iorb,iop2)
     2                      -f00(iorb,lorb,iop2)*f10(lorb,jorb,iop1)
     3                      +f00(jorb,lorb,iop2)*f10(lorb,iorb,iop1)
                  enddo
                  velfac=one/(two*cvel)
                  w(iorb,jorb)=w(iorb,jorb)+velfac*dum
               enddo
            enddo
            call giawt(ift1,w(1,1),ndim*ndim)
            if(dbug) then
               if(ifld.eq.1) then
            write(iw,*) 'b(x) perturbation with -k00- gauge independ.'
               elseif(ifld.eq.2) then
            write(iw,*) 'b(y) perturbation with -k00- gauge independ.'
               elseif(ifld.eq.3) then
            write(iw,*) 'b(z) perturbation with -k00- gauge independ.'
               endif
               call prsq(w(1,1),norb,norb,ndim)
            endif
         enddo
      enddo
      call rewfil(ift1)
      call rewfil(ift4)
c
      return
 9999 format(' in -giak00- , gauge for -ioc- = ',i5,
     1       ' gx,gy,gz = ',3f15.9)
      end
      subroutine gias10(s10,c,s00,t,num,noc,norb,ndim)
      implicit REAL (a-h,o-z)
      logical dbug
      logical out
      common/iofile/ir,iw
      common/filgia/ift1,ift2,ift3,ift4
      dimension s10(ndim*ndim,*)
      dimension c(ndim,*)
      dimension s00(*)
      dimension t(*)
      data zero /0.0d+00/
      data one  /1.0d+00/
      data two  /2.0d+00/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      if(out) then
         write(iw,*) 'start -gias10-'
         call fiflsh(iw)
      endif
      do iop=1,3
         do ij=1,ndim*ndim
            s10(ij,iop)=zero
         enddo
      enddo
c
c     ----- overlap derivatives -s10- -----
c
      call s10int(s10(1,1),s10(1,2),s10(1,3),s00,ndim)
c
c     ----- transform to canonical basis -----
c
      do iop=1,3
         if(dbug) then
            write(iw,*) '-s10- operators in -ao- basis'
            call prsq(s10(1,iop),num,num,ndim)
         endif
         call tfsq(s10(1,iop),s10(1,iop),c,t,norb,num,ndim)
         if(out) then
            write(iw,*) '-s10- operators in -mo- basis'
            call prsq(s10(1,iop),num,num,ndim)
         endif
      enddo
c
c     ----- write out -----
c
      call rewfil(ift1)
      call giawt(ift1,s10(1,1),ndim*ndim)
      call giawt(ift1,s10(1,2),ndim*ndim)
      call giawt(ift1,s10(1,3),ndim*ndim)
      call rewfil(ift1)
c
      if(out) then
         write(iw,*) 'end -gias10-'
         call fiflsh(iw)
      endif
c
      return
      end
      subroutine s10int(s10x,s10y,s10z,s00,ndim)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      logical norm
      logical dbug
      logical out
      common/iofile/ir,iw
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(mxatom),
     1                                          c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/baspar/normf,normp,itol
      common/xyzstv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension s10x(ndim,*),s10y(ndim,*),s10z(ndim,*)
      dimension  s00(ndim,*)
      dimension sx(225),sy(225),sz(225)
      dimension  s(225)
      dimension dij(225)
      dimension  xs(6,5), ys(6,5), zs(6,5)
      dimension xxs(6,5),yys(6,5),zzs(6,5)
      dimension ijx(35),ijy(35),ijz(35)
      data zero  /0.0d+00/
      data one   /1.0e+00/
      data rln10 /2.30258d+00/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      data ijx    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      data ijy    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      data ijz    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      tol=rln10*itol
      norm=normf.ne.1.or.normp.ne.1
c
c     ----- ishell -----
c
      do ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
      litmod=lit+1
c
c     ----- jshell -----
c
      do jj=1,nshell
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
c
      qijx=yi*zj-zi*yj
      qijy=zi*xj-xi*zj
      qijz=xi*yj-yi*xj
      tijx=xi-xj
      tijy=yi-yj
      tijz=zi-zj
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
c     ----- prepare indices for pairs of (i,j) functions -----
c
      ij=0
      do i=mini,maxi
         do j=minj,maxj
            ij=ij+1
            s (ij)=zero
            sx(ij)=zero
            sy(ij)=zero
            sz(ij)=zero
         enddo    
      enddo    
c
c     ----- i primitive -----
c
      do ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
      cgi=cg(ig)
c
c     ----- j primitive -----
c
      do jg=j1,j2
      aj=ex(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.le.tol) then       
      fac= exp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      cgj=cg(jg)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      ij=0
      do 360 i=mini,maxi
      go to (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),i
  110 dum1=csi*fac
      go to 220
  120 dum1=cpi*fac
      go to 220
  130 dum1=cdi*fac
      go to 220
  140 if(norm) dum1=dum1*sqrt3
      go to 220
  150 dum1=cfi*fac
      go to 220
  160 if(norm) dum1=dum1*sqrt5
      go to 220
  170 if(norm) dum1=dum1*sqrt3
      go to 220
  180 dum1=cgi*fac
      go to 220
  190 if(norm) dum1=dum1*sqrt7
      go to 220
  200 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 220
  210 if(norm) dum1=dum1*sqrt3
  220 continue
c
      do 360 j=minj,maxj
      go to (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),j
  230 dum2=dum1*csj
      go to 350
  250 dum2=dum1*cpj
      go to 350
  260 dum2=dum1*cdj
      go to 350
  270 if(norm) dum2=dum2*sqrt3
      go to 350
  280 dum2=dum1*cfj
      go to 350
  290 if(norm) dum2=dum2*sqrt5
      go to 350
  300 if(norm) dum2=dum2*sqrt3
      go to 350
  310 dum2=dum1*cgj
      go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      ij=ij+1
  360 dij(ij)=dum2
c
c     ----- x , y , z integrals -----
c
      t = sqrt(aa1)
      x0=ax
      y0=ay
      z0=az
      do j=1,ljt
         nj=j
         do i=1,litmod
            ni=i
            call sxyz
             xs(i,j)=xint*t
             ys(i,j)=yint*t 
             zs(i,j)=zint*t 
         enddo
      enddo
c
      do j=1,ljt
         do i=1,lit
            xxs(i,j)=xs(i+1,j)
            yys(i,j)=ys(i+1,j)
            zzs(i,j)=zs(i+1,j)
         enddo
      enddo
c
      ij=0
      do i=mini,maxi
         ix=ijx(i)
         iy=ijy(i)
         iz=ijz(i)
         do j=minj,maxj
            jx=ijx(j)
            jy=ijy(j)
            jz=ijz(j)
            ij=ij+1
            dum = xs(ix,jx)* ys(iy,jy)* zs(iz,jz)
            dumx=xxs(ix,jx)* ys(iy,jy)* zs(iz,jz)
            dumy= xs(ix,jx)*yys(iy,jy)* zs(iz,jz)
            dumz= xs(ix,jx)* ys(iy,jy)*zzs(iz,jz)
         s (ij)=s (ij)+dij(ij)*(     dum                    )
         sx(ij)=sx(ij)+dij(ij)*(qijx*dum+tijy*dumz-tijz*dumy)
         sy(ij)=sy(ij)+dij(ij)*(qijy*dum+tijz*dumx-tijx*dumz)
         sz(ij)=sz(ij)+dij(ij)*(qijz*dum+tijx*dumy-tijy*dumx)
         enddo
      enddo
c
            endif
         enddo
      enddo
c
      ij=0
      do i=mini,maxi
         do j=minj,maxj
            ij=ij+1
            iloc=loci+i
            jloc=locj+j
            s00 (iloc,jloc)=s (ij)
            s10x(iloc,jloc)=sx(ij)
            s10y(iloc,jloc)=sy(ij)
            s10z(iloc,jloc)=sz(ij)
         enddo
      enddo
c
         enddo
      enddo
c
      if(out) then
         write(iw,*) '-s00 - in ao basis'
         call prsq(s00 ,num,num,ndim)
         write(iw,*) '-s10x- in ao basis'
         call prsq(s10x,num,num,ndim)
         write(iw,*) '-s10y- in ao basis'
         call prsq(s10y,num,num,ndim)
         write(iw,*) '-s10z- in ao basis'
         call prsq(s10z,num,num,ndim)
      endif
c
      return
      end
      subroutine giah10(h10,c00,f00,d00,t,h00,
     1                  num,noc,norb,ndim,next)
      implicit REAL (a-h,o-z)
      logical times
      logical dbug
      logical more
      logical out
      common/iofile/ir,iw
      common/filgia/ift1,ift2,ift3,ift4
      dimension h10(ndim*ndim,*)
      dimension c00(ndim,*)
      dimension f00(*)
      dimension d00(*)
      dimension t(*)
      dimension h00(*)
      data zero /0.0d+00/
      data one  /1.0d+00/
      data two  /2.0d+00/
c
      dbug=.false.
      more=.false.
      more=more.or.dbug
      out =.false.
      out =out.or.more
c
      times=.false.
      times=times.or.out
c
      if(out) then
         write(iw,*) 'start -giah10-'
         call fiflsh(iw)
      endif
      do iop=1,6
         do ij=1,ndim*ndim
            h10(ij,iop)=zero
         enddo
      enddo
c
c     ----- l operator -----
c
      call gial10(h10(1,1),h10(1,2),h10(1,3),h00,ndim)
      call fiflsh(iw)
c
c     ----- t + v operators -----
c
      call gitv10(h10(1,4),h10(1,5),h10(1,6),h00,ndim)
c
c     ----- sum up into 1-e operators -----
c
      do iop=1,3
         do ij=1,ndim*ndim
            h10(ij,iop)=h10(ij,iop)+h10(ij,iop+3)
         enddo
      enddo
      if(out) then
         do iop=1,3
            write(iw,*) '1-e       operators in -ao- basis'
            if(more) then
               call prsq(h10(1,iop),num,num,ndim)
            endif
            call fiflsh(iw)
         enddo
      endif
c
c     ----- write out -----
c
      call rewfil(ift2)
      call giawt(ift2,h10(1,1),ndim*ndim)
      call giawt(ift2,h10(1,2),ndim*ndim)
      call giawt(ift2,h10(1,3),ndim*ndim)
      call rewfil(ift2)
c
c     ----- now 2-e operators -----
c
      call giag10(h10(1,4),f00,d00,ndim,next)
c
c     ----- global sum for parallel processing -----
c
      do iop=4,6
         call lcpadd(h10(1,iop),h00,ndim*ndim)
      enddo
c
      if(out) then
         do iop=4,6
            write(iw,*) '      2-e operators in -ao- basis'
            if(more) then
               call prsq(h10(1,iop),num,num,ndim)
            endif
            call fiflsh(iw)
         enddo
      endif
c
c     ----- transform to canonical basis -----
c
      do iop=1,3
         do ij=1,ndim*ndim
            h10(ij,iop)=(h10(ij,iop)+h10(ij,iop+3))
         enddo
         if(out) then
            write(iw,*) '1-e + 2-e operators in -ao- basis'
            if(more) then
               call prsq(h10(1,iop),num,num,ndim)
            endif 
            call fiflsh(iw)
         endif
         call tfsq(h10(1,iop),h10(1,iop),c00,t,norb,num,ndim)
         if(out) then
            write(iw,*) '1-e + 2-e operators in -mo- basis'
            if(more) then
               call prsq(h10(1,iop),num,num,ndim)
            endif 
            call fiflsh(iw)
         endif
      enddo
c
c     ----- write out -----
c
      call rewfil(ift2)
      call giawt(ift2,h10(1,1),ndim*ndim)
      call giawt(ift2,h10(1,2),ndim*ndim)
      call giawt(ift2,h10(1,3),ndim*ndim)
      call rewfil(ift2)
      call fiflsh(iw)
c
      if(out) then
         write(iw,*) 'end -giah10-'
         call fiflsh(iw)
      endif
c
      return
      end
      subroutine gial10(h10x,h10y,h10z,s00,ndim)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      logical more
      logical out
      common/iofile/ir,iw
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      dimension h10x(ndim,*),h10y(ndim,*),h10z(ndim,*)
      dimension s00 (ndim,*)
      dimension xl(225),yl(225),zl(225)
      dimension ss(225)
      data zero  /0.0d+00/
      data one   /1.0e+00/
c
      more=.false.
      out =.false.
      out =out.or.more
c
      if(out) then
         write(iw,*) 'start -gial10-'
         call fiflsh(iw)
      endif
c
c     ----- ishell -----
c
      do ii=1,nshell
         mini=kmin(ii)
         maxi=kmax(ii)
         loci=kloc(ii)-mini
c
c     ----- jshell -----
c
         do jj=1,nshell
            minj=kmin(jj)
            maxj=kmax(jj)
            locj=kloc(jj)-minj
c
c     ----- l integrals -----
c
            call l10int(ii,jj,xl,yl,zl,ss)
c
c     ----- l operators -----
c
            ij=0
            do i=mini,maxi
               do j=minj,maxj 
                  iloc=loci+i
                  jloc=locj+j
                  ij=ij+1
                  s00 (iloc,jloc)= ss(ij)
                  h10x(iloc,jloc)=-xl(ij)
                  h10y(iloc,jloc)=-yl(ij)
                  h10z(iloc,jloc)=-zl(ij)
               enddo
            enddo
c
         enddo
      enddo
c
      if(out) then        
         write(iw,*) '- s00    - operator in ao basis'
         if(more) then
            call prsq(s00 ,num,num,ndim)
         endif
         write(iw,*) '- l10 -x - operator in ao basis'
         if(more) then
            call prsq(h10x,num,num,ndim)
         endif
         call fiflsh(iw)
         write(iw,*) '- l10 -y - operator in ao basis'
         if(more) then
            call prsq(h10y,num,num,ndim)
         endif
         call fiflsh(iw)
         write(iw,*) '- l10 -z - operator in ao basis'
         if(more) then
            call prsq(h10z,num,num,ndim)
         endif
         call fiflsh(iw)
      endif
      if(out) then
         write(iw,*) 'end   -gial10-'
         call fiflsh(iw)
      endif
c
      return
      end
      subroutine l10int(ii,jj,xl,yl,zl,ss)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      logical dbug
      logical out
      logical norm
      common/iofile/ir,iw
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/baspar/normf,normp,itol
      common/xyzstv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      dimension xl(*),yl(*),zl(*),ss(*)
      dimension cij(225)
      dimension  xs(5,6), ys(5,6), zs(5,6)
      dimension xsx(5,6),ysy(5,6),zsz(5,6)
      dimension xsd(5,6),ysd(5,6),zsd(5,6)
      dimension ijx(35),ijy(35),ijz(35)
      data zero  /0.0d+00/
      data one   /1.0e+00/
      data rln10 /2.30258d+00/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      data ijx    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      data ijy    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      data ijz    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      if(out) then
         write(iw,*) 'in    -h10int-',ii,jj
         call fiflsh(iw)
      endif
c
      tol=rln10*itol
      norm=normf.ne.1.or.normp.ne.1
c
c     ----- ishell -----
c
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell -----
c
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
c
      ljtmod=ljt+1
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
      ij=0
      do i=mini,maxi
         do j=minj,maxj
            ij=ij+1
            ss(ij)=zero
            xl(ij)=zero
            yl(ij)=zero
            zl(ij)=zero
         enddo
      enddo
c
c     ----- i primitive -----
c
      do ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
      cgi=cg(ig)
c
c     ----- j primitive -----
c
      do jg=j1,j2
      aj=ex(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.le.tol) then        
      fac= exp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      cgj=cg(jg)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      ij=0
      do 360 i=mini,maxi
      go to (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),i
  110 dum1=csi*fac
      go to 220
  120 dum1=cpi*fac
      go to 220
  130 dum1=cdi*fac
      go to 220
  140 if(norm) dum1=dum1*sqrt3
      go to 220
  150 dum1=cfi*fac
      go to 220
  160 if(norm) dum1=dum1*sqrt5
      go to 220
  170 if(norm) dum1=dum1*sqrt3
      go to 220
  180 dum1=cgi*fac
      go to 220
  190 if(norm) dum1=dum1*sqrt7
      go to 220
  200 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 220
  210 if(norm) dum1=dum1*sqrt3
  220 continue
c
      do 360 j=minj,maxj
      go to (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),j
  230 dum2=dum1*csj
      go to 350
  250 dum2=dum1*cpj
      go to 350
  260 dum2=dum1*cdj
      go to 350
  270 if(norm) dum2=dum2*sqrt3
      go to 350
  280 dum2=dum1*cfj
      go to 350
  290 if(norm) dum2=dum2*sqrt5
      go to 350
  300 if(norm) dum2=dum2*sqrt3
      go to 350
  310 dum2=dum1*cgj
      go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      ij=ij+1
  360 cij(ij)=dum2
c
c     ----- integral components -----
c
      t = sqrt(aa1)
      x0=ax
      y0=ay
      z0=az
      do j=1,ljtmod
         nj=j
         do i=1,lit
            ni=i
            call sxyz
             xs(i,j)=xint*t
             ys(i,j)=yint*t
             zs(i,j)=zint*t
         enddo
      enddo
c
      do j=1,ljt
         do i=1,lit
            xsx(i,j)=xs(i,j+1)
            ysy(i,j)=ys(i,j+1)
            zsz(i,j)=zs(i,j+1)
         enddo
      enddo
c
      call l10xyz(xsd,ysd,zsd,xs,ys,zs,lit,ljt,aj)
c
c     ----- l integrals -----
c
      ij=0
      do i=mini,maxi
         ix=ijx(i)
         iy=ijy(i)
         iz=ijz(i)
         do j=minj,maxj
            jx=ijx(j)
            jy=ijy(j)
            jz=ijz(j)
            ij=ij+1
            dum =  xs(ix,jx)* ys(iy,jy)* zs(iz,jz)
            dumx=  xs(ix,jx)*ysy(iy,jy)*zsd(iz,jz)
     1           - xs(ix,jx)*ysd(iy,jy)*zsz(iz,jz)
            dumy= xsd(ix,jx)* ys(iy,jy)*zsz(iz,jz)
     1           -xsx(ix,jx)* ys(iy,jy)*zsd(iz,jz)
            dumz= xsx(ix,jx)*ysd(iy,jy)* zs(iz,jz)
     1           -xsd(ix,jx)*ysy(iy,jy)* zs(iz,jz)
            ss(ij)=ss(ij)+cij(ij)*dum 
            xl(ij)=xl(ij)+cij(ij)*dumx
            yl(ij)=yl(ij)+cij(ij)*dumy
            zl(ij)=zl(ij)+cij(ij)*dumz
         enddo
      enddo
c
            endif
         enddo
      enddo
c
      if(out) then
         write(iw,*) 'out   -h10int-'
         call fiflsh(iw)
      endif
      return
      end
      subroutine l10xyz(dx,dy,dz,x,y,z,ni,nj,aj)
      implicit REAL (a-h,o-z)
      dimension  x(5,*), y(5,*), z(5,*)
      dimension dx(5,*),dy(5,*),dz(5,*)
c
      do i=1,ni
         dx(i,1)= (-(aj+aj)*x(i,2))
         dy(i,1)= (-(aj+aj)*y(i,2))
         dz(i,1)= (-(aj+aj)*z(i,2))
      enddo
c
      if(nj.eq.1) return
c
      do j=2,nj
         do i=1,ni
            dx(i,j)= (dble(j-1)*x(i,j-1)-(aj+aj)*x(i,j+1))
            dy(i,j)= (dble(j-1)*y(i,j-1)-(aj+aj)*y(i,j+1))
            dz(i,j)= (dble(j-1)*z(i,j-1)-(aj+aj)*z(i,j+1))
         enddo
      enddo
c
      return
      end
      subroutine gitv10(h10x,h10y,h10z,h00,ndim)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      character*8 errmsg
      logical norm
      logical dbug
      logical out
      common/output/nprint
      common/iofile/ir,iw
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(mxatom),
     1                                          c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/baspar/normf,normp,itol
      common/xyzstv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/root/xx,u(12),w(12),nroots
      dimension h10x(ndim,*),h10y(ndim,*),h10z(ndim,*)
      dimension h00 (ndim,*)
      dimension cij(225)
      dimension g(225)
      dimension v(225)
      dimension  gx(225),gy(225),gz(225)
      dimension  vx(225),vy(225),vz(225)
      dimension ijx(35),ijy(35),ijz(35)
      dimension  xs(6,7)  , ys(6,7)  , zs(6,7)
      dimension  xg(6,5)  , yg(6,5)  , zg(6,5)
      dimension  xv(6,5,5), yv(6,5,5), zv(6,5,5)
      dimension xxs(6,7)  ,yys(6,7)  ,zzs(6,7)
      dimension xxg(6,5)  ,yyg(6,5)  ,zzg(6,5)
      dimension xxv(6,5,5),yyv(6,5,5),zzv(6,5,5)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','-gitv10-'/
      data maxrys /5/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data pi212  /1.1283791670955d+00/
      data sqrt3  /1.73205080756888d+00/
      data sqrt5  /2.23606797749979d+00/
      data sqrt7  /2.64575131106459d+00/
      data ijx    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      data ijy    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      data ijz    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      if(out) then
         write(iw,*) 'start -gitv10-'
         call fiflsh(iw)
      endif
c
      tol=rln10*itol
      norm=normf.ne.1.or.normp.ne.1
c
c     ----- ishell -----
c
      do ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
      litmod=lit+1
c
c     ----- jshell -----
c
      do jj=1,nshell
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
c
      ljtmod=ljt+2
c
      qijx=yi*zj-zi*yj
      qijy=zi*xj-xi*zj
      qijz=xi*yj-yi*xj
      tijx=xi-xj
      tijy=yi-yj
      tijz=zi-zj
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
      nroots=(1+lit+ljt-2)/2+1
      if(nroots.gt.maxrys) then
         write(iw,9997) maxrys,lit,ljt,nroots
         call hnderr(3,errmsg)
      endif
c
      ij=0
      do i=mini,maxi
         do j=minj,maxj
            ij=ij+1
            g (ij)=zero
            gx(ij)=zero
            gy(ij)=zero
            gz(ij)=zero
            v (ij)=zero
            vx(ij)=zero
            vy(ij)=zero
            vz(ij)=zero
         enddo
      enddo
c
c     ----- i primitive -----
c
      do ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
      cgi=cg(ig)
c
c     ----- j primitive -----
c
      do jg=j1,j2
      aj=ex(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.le.tol) then       
      fac= exp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      cgj=cg(jg)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      ij=0
      do 360 i=mini,maxi
      go to (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),i
  110 dum1=csi*fac
      go to 220
  120 dum1=cpi*fac
      go to 220
  130 dum1=cdi*fac
      go to 220
  140 if(norm) dum1=dum1*sqrt3
      go to 220
  150 dum1=cfi*fac
      go to 220
  160 if(norm) dum1=dum1*sqrt5
      go to 220
  170 if(norm) dum1=dum1*sqrt3
      go to 220
  180 dum1=cgi*fac
      go to 220
  190 if(norm) dum1=dum1*sqrt7
      go to 220
  200 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 220
  210 if(norm) dum1=dum1*sqrt3
  220 continue
c
      do 360 j=minj,maxj
      go to (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),j
  230 dum2=dum1*csj
      go to 350
  250 dum2=dum1*cpj
      go to 350
  260 dum2=dum1*cdj
      go to 350
  270 if(norm) dum2=dum2*sqrt3
      go to 350
  280 dum2=dum1*cfj
      go to 350
  290 if(norm) dum2=dum2*sqrt5
      go to 350
  300 if(norm) dum2=dum2*sqrt3
      go to 350
  310 dum2=dum1*cgj
      go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      ij=ij+1
  360 cij(ij)=dum2
c
c     ----- kinetic energy -----
c
      t = sqrt(aa1)
      x0=ax
      y0=ay
      z0=az
      do j=1,ljtmod
         nj=j
         do i=1,litmod
            ni=i
            call sxyz
            xs(i,j)=xint*t
            ys(i,j)=yint*t
            zs(i,j)=zint*t
         enddo
      enddo
c
      do j=1,ljtmod
         do i=1,lit
            xxs(i,j)=xs(i+1,j)
            yys(i,j)=ys(i+1,j)
            zzs(i,j)=zs(i+1,j)
         enddo
      enddo
c
      call t10xyz( xg, yg, zg, xs, ys, zs,lit,ljt,aj)
      call t10xyz(xxg,yyg,zzg,xxs,yys,zzs,lit,ljt,aj)
c
      ij=0
      do i=mini,maxi
         ix=ijx(i)
         iy=ijy(i)
         iz=ijz(i)
         do j=minj,maxj
            jx=ijx(j)
            jy=ijy(j)
            jz=ijz(j)
            ij=ij+1
            dum = xg(ix,jx)* ys(iy,jy)* zs(iz,jz)
     1          + xs(ix,jx)* yg(iy,jy)* zs(iz,jz)
     2          + xs(ix,jx)* ys(iy,jy)* zg(iz,jz)
            dumx=xxg(ix,jx)* ys(iy,jy)* zs(iz,jz)
     1          +xxs(ix,jx)* yg(iy,jy)* zs(iz,jz)
     2          +xxs(ix,jx)* ys(iy,jy)* zg(iz,jz)
            dumy= xg(ix,jx)*yys(iy,jy)* zs(iz,jz)
     1          + xs(ix,jx)*yyg(iy,jy)* zs(iz,jz)
     2          + xs(ix,jx)*yys(iy,jy)* zg(iz,jz)
            dumz= xg(ix,jx)* ys(iy,jy)*zzs(iz,jz)
     1          + xs(ix,jx)* yg(iy,jy)*zzs(iz,jz)
     2          + xs(ix,jx)* ys(iy,jy)*zzg(iz,jz)
      g (ij)=g (ij)+cij(ij)*(     dum                    )
      gx(ij)=gx(ij)+cij(ij)*(qijx*dum+tijy*dumz-tijz*dumy)
      gy(ij)=gy(ij)+cij(ij)*(qijy*dum+tijz*dumx-tijx*dumz)
      gz(ij)=gz(ij)+cij(ij)*(qijz*dum+tijx*dumy-tijy*dumx)
         enddo
      enddo
c
c     ----- nuclear attraction -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do ic=1,nat
      znuc=-zan(ic)
      cx=c(1,ic)
      cy=c(2,ic)
      cz=c(3,ic)
      xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots.le.3) call rt123
      if(nroots.eq.4) call root4
      if(nroots.eq.5) call root5
      do iroot=1,nroots
         uu=u(iroot)*aa
         ww=w(iroot)*znuc
         tt=one/(aa+uu)
         t = sqrt(tt)
         x0=(aax+uu*cx)*tt
         y0=(aay+uu*cy)*tt
         z0=(aaz+uu*cz)*tt
         do j=1,ljt
            nj=j
            do i=1,litmod
               ni=i
               call sxyz
               xv(i,j,iroot)=xint
               yv(i,j,iroot)=yint
               zv(i,j,iroot)=zint*ww
            enddo
         enddo
      enddo
c
      do iroot=1,nroots
         do j=1,ljt
            do i=1,lit
               xxv(i,j,iroot)=xv(i+1,j,iroot)
               yyv(i,j,iroot)=yv(i+1,j,iroot)
               zzv(i,j,iroot)=zv(i+1,j,iroot)
            enddo
         enddo
      enddo
c
      ij=0
      do i=mini,maxi
         ix=ijx(i)
         iy=ijy(i)
         iz=ijz(i)
         do j=minj,maxj
            jx=ijx(j)
            jy=ijy(j)
            jz=ijz(j)
            dum =zero
            dumx=zero
            dumy=zero
            dumz=zero
            do iroot=1,nroots
      dum =dum + xv(ix,jx,iroot)* yv(iy,jy,iroot)* zv(iz,jz,iroot)
      dumx=dumx+xxv(ix,jx,iroot)* yv(iy,jy,iroot)* zv(iz,jz,iroot)
      dumy=dumy+ xv(ix,jx,iroot)*yyv(iy,jy,iroot)* zv(iz,jz,iroot)
      dumz=dumz+ xv(ix,jx,iroot)* yv(iy,jy,iroot)*zzv(iz,jz,iroot)
            enddo
            dum =dum *(aa1*pi212)
            dumx=dumx*(aa1*pi212)
            dumy=dumy*(aa1*pi212)
            dumz=dumz*(aa1*pi212)
            ij=ij+1
      v (ij)=v (ij)+cij(ij)*(     dum                    )
      vx(ij)=vx(ij)+cij(ij)*(qijx*dum+tijy*dumz-tijz*dumy)
      vy(ij)=vy(ij)+cij(ij)*(qijy*dum+tijz*dumx-tijx*dumz)
      vz(ij)=vz(ij)+cij(ij)*(qijz*dum+tijx*dumy-tijy*dumx)
         enddo
      enddo
c
c     ----- end loop over atoms -----
c
      enddo
c
c     ----- end loop over primitives -----
c
            endif
         enddo
      enddo
c
      ij=0
      do i=mini,maxi
         do j=minj,maxj
            ij=ij+1
            iloc=loci+i
            jloc=locj+j
            h00 (iloc,jloc)=g (ij)+v (ij)
            h10x(iloc,jloc)=gx(ij)+vx(ij)
            h10y(iloc,jloc)=gy(ij)+vy(ij)
            h10z(iloc,jloc)=gz(ij)+vz(ij)
         enddo
      enddo
c
c     ----- end loop over shells -----
c
         enddo
      enddo
c
      if(out) then
         write(iw,*) '- h00    - in ao basis'
         if(dbug) then
            call prsq(h00 ,num,num,ndim)
         endif
         write(iw,*) '- tv10-x - in ao basis'
         if(dbug) then
            call prsq(h10x,num,num,ndim)
         endif
         write(iw,*) '- tv10-y - in ao basis'
         if(dbug) then
            call prsq(h10y,num,num,ndim)
         endif
         write(iw,*) '- tv10-z - in ao basis'
         if(dbug) then
            call prsq(h10z,num,num,ndim)
         endif
         call fiflsh(iw)
         write(iw,*) 'end   -gitv10-'
         call fiflsh(iw)
      endif
      return
 9997 format(' in -gitv10- the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i2,/,' lit,ljt,nroots= ',3i3)
      end
      subroutine t10xyz(xt,yt,zt,xs,ys,zs,ni,nj,aj)
      implicit REAL (a-h,o-z)
      dimension xt(6,*),yt(6,*),zt(6,*)
      dimension xs(6,*),ys(6,*),zs(6,*)
c
      do i=1,ni
         xt(i,1)=(xs(i,1  )            -xs(i,3  )*(aj+aj))*aj
         yt(i,1)=(ys(i,1  )            -ys(i,3  )*(aj+aj))*aj
         zt(i,1)=(zs(i,1  )            -zs(i,3  )*(aj+aj))*aj
      enddo     
c
      if(nj.eq.1) return
c
      do i=1,ni
         xt(i,2)=(xs(i,2  )*dble(2+2-1)-xs(i,4  )*(aj+aj))*aj
         yt(i,2)=(ys(i,2  )*dble(2+2-1)-ys(i,4  )*(aj+aj))*aj
         zt(i,2)=(zs(i,2  )*dble(2+2-1)-zs(i,4  )*(aj+aj))*aj
      enddo
c
      if(nj.eq.2) return
c
      do j=3,nj
         do i=1,ni
         xt(i,j)=(xs(i,j  )*dble(j+j-1)-xs(i,j+2)*(aj+aj))*aj
     1           -xs(i,j-2)*dble(((j-1)*(j-2))/2)
         yt(i,j)=(ys(i,j  )*dble(j+j-1)-ys(i,j+2)*(aj+aj))*aj
     1           -ys(i,j-2)*dble(((j-1)*(j-2))/2)
         zt(i,j)=(zs(i,j  )*dble(j+j-1)-zs(i,j+2)*(aj+aj))*aj
     1           -zs(i,j-2)*dble(((j-1)*(j-2))/2)
         enddo
      enddo
c
      return
      end
      subroutine giah20(h20,c,b20,t,num,noc,ndim)
      implicit REAL (a-h,o-z)
      parameter (maxorb=MAXORB)
      logical dbug
      logical out
      common/iofile/ir,iw
INCLUDE(../m4/common/mapper)
      dimension h20(ndim*ndim,*)
      dimension c(ndim,*)
      dimension t(*)
      dimension cm(3)     
      dimension rr(6)     
      dimension b20(6)
      data cvel  /137.0359895d+00/
      data zero  /0.0d+00/
      data one   /1.0d+00/
      data two   /2.0d+00/
      data eight /8.0d+00/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      velfac=one/(eight*cvel*cvel)
c
c     ----- r**2 operators : xx yy zz xy xz yz -----
c         for simplicity, recalculate for each mo
c
      do iop=1,6
         rr(iop)=zero 
      enddo
c
      do ioc=1,noc
c
c     ----- r**2 integrals : xx yy zz xy xz yz -----
c
         cm(1)=zero
         cm(2)=zero
         cm(3)=zero
         call qdpint(h20(1,1),h20(1,2),h20(1,3),
     1               h20(1,4),h20(1,5),h20(1,6),cm,ndim)
         if(out.and.ioc.eq.1) then
            do iop=1,6
               write(iw,*) 'r**2 operator, iop,ioc = ',iop,ioc
               call prtr(h20(1,iop),num)
            enddo
         endif
c
c     ----- susceptibility contribution -----
c
         do iop=1,6
            do i=1,num
               dum=zero
               do j=1,num
                  ij=iky(max(i,j))+min(i,j)
                  dum=dum+h20(ij,iop)*c(j,ioc)
               enddo
               t(i)=dum
            enddo
            dum=zero
            do i=1,num
               dum=dum+c(i,ioc)*t(i)
            enddo
            rr(iop)=rr(iop)+dum 
         enddo
      enddo
c
      do iop=1,6
         rr(iop)=rr(iop)*two
      enddo
      if(out) then
         write(iw,9998)
         write(iw,9997) (rr(iop),iop=1,6)
      endif
c
c     ----- diamagnetic susceptibility tensor -----
c
      b20(1)=rr(2)+rr(3)
      b20(2)=rr(1)+rr(3)
      b20(3)=rr(1)+rr(2)
      b20(4)=-two*rr(4)
      b20(5)=-two*rr(5)
      b20(6)=-two*rr(6)
      if(out) then
         write(iw,9996)
         write(iw,9997) (b20(iop),iop=1,6)
      endif
c
      return
 9998 format(' in -giah20- , -r2- tensor = ')
 9997 format('    xx = ',f20.10,/,
     1       '    yy = ',f20.10,/,
     2       '    zz = ',f20.10,/,
     3       '    xy = ',f20.10,/,
     4       '    xz = ',f20.10,/,
     5       '    yz = ',f20.10,/)
 9996 format(' in -giah20- , -b20- tensor = ')
      end
      subroutine giae11(e01,e11,h01,d10,d00,ndim,ifld)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (zero=0.0d+00)
      parameter (two =2.0d+00)
      logical   out
      logical   dbug
      common/iofile/ir,iw
      common/filgia/ift1,ift2,ift3,ift4
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(mxatom),
     1                                          c(3,mxatom)
      common/lcapid/nap,iap
      dimension e01(3,*),e11(3,3,*)            
      dimension h01(ndim,*)
      dimension d10(ndim,*)
      dimension d00(*)
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      do kat=1,nat
         do kxyz=1,3
            e11(kxyz,ifld,kat)=zero
         enddo
      enddo
c
      call rewfil(ift3)
      do kat=1,nat
         if(mod(kat-1,nap)+1.eq.iap) then
            if(out) then
               write(iw,*) 'in -giae11, kat,iap = ',kat,iap
            endif
c
         do kxyz=1,3
            call giard(ift3,h01,ndim*ndim)
            dum=zero
            do j=1,num
               do i=1,num
                  dum=dum+d10(i,j)*h01(i,j)
               enddo
            enddo
            dum=dum*two
            e11(kxyz,ifld,kat)=dum
         enddo
         if(out) then
            write(iw,9999) ifld 
            if(dbug) then
               call prsq(e11(1,1,kat),3,3,3)
            endif
         endif
c
         endif
      enddo
c
      return
 9999 format(' in -giae11- -e11- for -ifld- = ',i5)
      end
      subroutine giah01(h01,c00,uuu,ttt,nbf,noc,norb,ndim)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      character*8 errmsg
      logical norm
      logical some
      logical out
      logical dbug
      common/lcapid/nap,iap
      common/filgia/ift1,ift2,ift3,ift4
      common/iofile/ir,iw
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(mxatom),
     1                                          c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/baspar/normf,normp,itol
      common/root/yy,u5(12),w5(12),nroots
      common/hfk/xx,u9(9),w9(9),mroots
      common/xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
     1                              ,cx,cy,cz
      dimension h01(ndim,ndim,*)
      dimension c00(*)
      dimension ttt(*)
      dimension uuu(*)
      dimension cij(225)
      dimension h01x(225),h01y(225),h01z(225)
      dimension ijx(35),ijy(35),ijz(35)
      dimension   xv(5,6,6),  yv(5,6,6),  zv(5,6,6)
      dimension  dxv(5,6,6), dyv(5,6,6), dzv(5,6,6)
      dimension  xvx(5,6,6), yvy(5,6,6), zvz(5,6,6)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','-giah01-'/
      data maxrys /6/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data pi212  /1.1283791670955d+00/
      data sqrt3  /1.73205080756888d+00/
      data sqrt5  /2.23606797749979d+00/
      data sqrt7  /2.64575131106459d+00/
      data ijx    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      data ijy    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      data ijz    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      some=.false.
      some=some.or.out
c
      tol=rln10*itol
      norm=normf.ne.1.or.normp.ne.1
c
      npwr=1+1
c
c     ----- kat -----
c
      call rewfil(ift3)
      do kat=1,nat
         if(mod(kat-1,nap)+1.eq.iap) then
            if(out) then
               write(iw,*) 'in -giah01-, kat,iap = ',
     1                                   kat,iap
            endif
c
c     ----- ishell -----
c
      do ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell -----
c
      do jj=1,nshell
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
c
      ljtmod=ljt+1
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt+npwr-2)/2+1
      mroots=nroots
      if(nroots.gt.maxrys) then
         write(iw,9999) maxrys,lit,ljt,nroots
         call hnderr(3,errmsg)
      endif
c
c     ----- prepare indices for pairs of (i,j) functions -----
c
      ij=0
      do i=mini,maxi
         do j=minj,maxj
            ij=ij+1
            h01x(ij)=zero
            h01y(ij)=zero
            h01z(ij)=zero
         enddo
      enddo
c
c     ----- i primitive -----
c
      do ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
      cgi=cg(ig)
c
c     ----- j primitive -----
c
      do jg=j1,j2
      aj=ex(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.le.tol) then       
      fac= exp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      cgj=cg(jg)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      ij=0
      do 360 i=mini,maxi
      go to (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),i
  110 dum1=csi*fac
      go to 220
  120 dum1=cpi*fac
      go to 220
  130 dum1=cdi*fac
      go to 220
  140 if(norm) dum1=dum1*sqrt3
      go to 220
  150 dum1=cfi*fac
      go to 220
  160 if(norm) dum1=dum1*sqrt5
      go to 220
  170 if(norm) dum1=dum1*sqrt3
      go to 220
  180 dum1=cgi*fac
      go to 220
  190 if(norm) dum1=dum1*sqrt7
      go to 220
  200 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 220
  210 if(norm) dum1=dum1*sqrt3
  220 continue
c
      do 360 j=minj,maxj 
      go to (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),j
  230 dum2=dum1*csj
      go to 350
  250 dum2=dum1*cpj
      go to 350
  260 dum2=dum1*cdj
      go to 350
  270 if(norm) dum2=dum2*sqrt3
      go to 350
  280 dum2=dum1*cfj
      go to 350
  290 if(norm) dum2=dum2*sqrt5
      go to 350
  300 if(norm) dum2=dum2*sqrt3
      go to 350
  310 dum2=dum1*cgj
      go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      ij=ij+1
      cij(ij)=dum2
  360 continue
c
c     ----- -h01- term -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
c
         znuc=one
         cx=c(1,kat)
         cy=c(2,kat)
         cz=c(3,kat)
         xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
         yy=xx
         if(nroots.le.3) call rt123
         if(nroots.eq.4) call root4
         if(nroots.eq.5) call root5
         if(nroots.eq.6) then
            call droot
         else
            do iroot=1,nroots
               u9(iroot)=u5(iroot)
               w9(iroot)=w5(iroot)
            enddo
         endif
         do iroot=1,nroots
            uu=u9(iroot)*aa
            ww=w9(iroot)*znuc
            ww=ww*(uu+uu)
            tt=one/(aa+uu)
            t = sqrt(tt)
            x0=(aax+uu*cx)*tt
            y0=(aay+uu*cy)*tt
            z0=(aaz+uu*cz)*tt
            do j=1,ljtmod
               nj=j
               do i=1,lit
                  ni=i
                  call dsxyz
                   xv(i,j,iroot)=xint
                   yv(i,j,iroot)=yint
                   zv(i,j,iroot)=zint*ww
                  call dvxyz
                  dxv(i,j,iroot)=xint
                  dyv(i,j,iroot)=yint
                  dzv(i,j,iroot)=zint*ww
               enddo
            enddo
         enddo
c
         do iroot=1,nroots
            call h01xyz(
     1            xvx(1,1,iroot),yvy(1,1,iroot),zvz(1,1,iroot),
     2             xv(1,1,iroot), yv(1,1,iroot), zv(1,1,iroot),
     3           lit,ljt,aj)
         enddo
c
         ij=0
         do i=mini,maxi
            ix=ijx(i)
            iy=ijy(i)
            iz=ijz(i)
            do j=minj,maxj
               jx=ijx(j)
               jy=ijy(j)
               jz=ijz(j)
               dumx=zero
               dumy=zero
               dumz=zero
               do iroot=1,nroots
                  dumx=dumx+
     1   xv(ix,jx,iroot)*dyv(iy,jy,iroot)*zvz(iz,jz,iroot)
     2 - xv(ix,jx,iroot)*yvy(iy,jy,iroot)*dzv(iz,jz,iroot)
                  dumy=dumy+
     1  xvx(ix,jx,iroot)* yv(iy,jy,iroot)*dzv(iz,jz,iroot)
     2 -dxv(ix,jx,iroot)* yv(iy,jy,iroot)*zvz(iz,jz,iroot)
                  dumz=dumz+
     1  dxv(ix,jx,iroot)*yvy(iy,jy,iroot)* zv(iz,jz,iroot)
     2 -xvx(ix,jx,iroot)*dyv(iy,jy,iroot)* zv(iz,jz,iroot)
               enddo     
               ij=ij+1
               dum=pi212*aa1*cij(ij)
               dumx =dumx*dum
               dumy =dumy*dum
               dumz =dumz*dum
               h01x(ij)=h01x(ij)+dumx
               h01y(ij)=h01y(ij)+dumy
               h01z(ij)=h01z(ij)+dumz
            enddo     
         enddo     
c
            endif
         enddo      
      enddo      
c
      ij=0
      do i=mini,maxi
         do j=minj,maxj
            iloc=loci+i
            jloc=locj+j
            ij=ij+1
            h01(iloc,jloc,1)= h01x(ij)
            h01(iloc,jloc,2)= h01y(ij)
            h01(iloc,jloc,3)= h01z(ij)
         enddo
      enddo
      if(dbug) then
         ij=0
         do i=mini,maxi
            do j=minj,maxj
               ij=ij+1
               write(iw,9998) ii,jj,i,j,h01x(ij),h01y(ij),
     1                                           h01z(ij)
            enddo
         enddo
      endif
c
         enddo      
      enddo      
c
c     ----- print -----
c
      if(out) then
         write(iw,9997) kat
         if(dbug) then
            call prsq(h01(1,1,1),num,num,ndim)
         endif
         write(iw,9996) kat
         if(dbug) then
            call prsq(h01(1,1,2),num,num,ndim)
         endif
         write(iw,9995) kat
         if(dbug) then
            call prsq(h01(1,1,3),num,num,ndim)
         endif
      endif
c
c     ----- write to -ift3- -----
c
      call giawt(ift3,h01(1,1,1),ndim*ndim)
      call giawt(ift3,h01(1,1,2),ndim*ndim)
      call giawt(ift3,h01(1,1,3),ndim*ndim)
c
      if(out) then
         call tfsq(uuu,h01(1,1,1),c00,ttt,norb,num,ndim)
         write(iw,9994) kat
         if(dbug) then
            call prsq(uuu,norb,norb,ndim)
         endif
         call tfsq(uuu,h01(1,1,2),c00,ttt,norb,num,ndim)
         write(iw,9993) kat
         if(dbug) then
            call prsq(uuu,norb,norb,ndim)
         endif
         call tfsq(uuu,h01(1,1,3),c00,ttt,norb,num,ndim)
         write(iw,9992) kat
         if(dbug) then
            call prsq(uuu,norb,norb,ndim)
         endif
      endif
c
         endif
      enddo
c
      call rewfil(ift3)
c
      return
 9999 format(' in -giah01- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3,/,' lit,ljt,nroots= ',3i3)
 9998 format(' ii,jj,i,j,h01x,y,z ... = ',4i4,3f15.10)
 9997 format(' -h01 x- in -ao- for atom = ',i5)
 9996 format(' -h01 y- in _ao- for atom = ',i5)
 9995 format(' -h01 z- in -ao- for atom = ',i5)
 9994 format(' -h01 x- in -mo- for atom = ',i5)
 9993 format(' -h01 y- in _mo- for atom = ',i5)
 9992 format(' -h01 z- in -mo- for atom = ',i5)
      end
      subroutine h01xyz(dx,dy,dz,x,y,z,ni,nj,aj)
      implicit REAL (a-h,o-z)
      dimension  x(5,*), y(5,*), z(5,*)
      dimension dx(5,*),dy(5,*),dz(5,*)
c
      do i=1,ni
         dx(i,1)= (-(aj+aj)*x(i,2))
         dy(i,1)= (-(aj+aj)*y(i,2))
         dz(i,1)= (-(aj+aj)*z(i,2))
      enddo
c
      if(nj.eq.1) return
c
      do j=2,nj
         do i=1,ni
            dx(i,j)= (dble(j-1)*x(i,j-1)-(aj+aj)*x(i,j+1))
            dy(i,j)= (dble(j-1)*y(i,j-1)-(aj+aj)*y(i,j+1))
            dz(i,j)= (dble(j-1)*z(i,j-1)-(aj+aj)*z(i,j+1))
         enddo
      enddo
c
      return
      end
      subroutine giah11(e11d,e11p,d0)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxshel=MXSHEL)
      parameter (mxprim=MXPRIM)
      parameter (maxorb=MAXORB)
      character*8 errmsg
      logical norm
      logical some
      logical out
      logical more
      logical dbug
      common/iofile/ir,iw
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(mxatom),
     1                                          c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/baspar/normf,normp,itol
      common/root/yy,u5(12),w5(12),nroots
      common/hfk/xx,u9(9),w9(9),mroots
      common/xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
     1                              ,cx,cy,cz
INCLUDE(../m4/common/mapper)
      dimension e11d(3,3,*),e11p(3,3,*)            
      dimension d0(*)
      dimension cij(225)
      dimension dij(225)
      dimension v00(4,225)
      dimension h01(3,225)
      dimension h11d(9,225)
      dimension h11p(9,225)
      dimension h11t(9,225)
      dimension ijx(35),ijy(35),ijz(35)
      dimension   vx0 (6,6,6),  vy0 (6,6,6),  vz0 (6,6,6)
      dimension  dvx0 (6,6,6), dvy0 (6,6,6), dvz0 (6,6,6)
      dimension   vx0d(6,6,6),  vy0d(6,6,6),  vz0d(6,6,6)
      dimension   vx1 (6,6,6),  vy1 (6,6,6),  vz1 (6,6,6)
      dimension  dvx1 (6,6,6), dvy1 (6,6,6), dvz1 (6,6,6)
      dimension   vx1x(6,6,6),  vy1y(6,6,6),  vz1z(6,6,6)
      dimension  dvx1x(6,6,6), dvy1y(6,6,6), dvz1z(6,6,6)
      dimension   vx1d(6,6,6),  vy1d(6,6,6),  vz1d(6,6,6)
      dimension  xvx1 (6,6,6), yvy1 (6,6,6), zvz1 (6,6,6)
      dimension xdvx1 (6,6,6),ydvy1 (6,6,6),zdvz1 (6,6,6)
      dimension  xvx1d(6,6,6), yvy1d(6,6,6), zvz1d(6,6,6)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','-giah11-'/
      data maxrys /6/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data pi212  /1.1283791670955d+00/
      data sqrt3  /1.73205080756888d+00/
      data sqrt5  /2.23606797749979d+00/
      data sqrt7  /2.64575131106459d+00/
      data ijx    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      data ijy    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      data ijz    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      some=.false.
      some=some.or.out
c
      tol=rln10*itol
      norm=normf.ne.1.or.normp.ne.1
c
      do kat=1,nat
         do j=1,3
            do i=1,3
               e11d(i,j,kat)=zero
               e11p(i,j,kat)=zero
            enddo
         enddo
      enddo
c
c     ----- ishell -----
c
      do ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
      litmod=lit+1
c
c     ----- jshell -----
c
      do jj=1,nshell
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
c
      ljtmod=ljt+1
c
      qijx=yi*zj-zi*yj
      qijy=zi*xj-xi*zj
      qijz=xi*yj-yi*xj
      tijx=xi   -   xj
      tijy=yi   -   yj
      tijz=zi   -   zj
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
      npwr=1+1+1
      nroots=(lit+ljt+npwr-2)/2+1
      mroots=nroots
      if(nroots.gt.maxrys) then
         write(iw,9999) maxrys,lit,ljt,nroots
         call hnderr(3,errmsg)
      endif
c
      if(out) then
         ij=0
         do i=mini,maxi
            ix=ijx(i)
            iy=ijy(i)
            iz=ijz(i)
            do j=minj,maxj
               jx=ijx(j)
               jy=ijy(j)
               jz=ijz(j)
               ij=ij+1
               do ix =1,4
                  v00(ix ,ij)=zero
               enddo
               do ix =1,3
                  h01(ix ,ij)=zero
               enddo
               do ixy=1,9
                  h11d(ixy,ij)=zero
                  h11p(ixy,ij)=zero
                  h11t(ixy,ij)=zero
               enddo
            enddo
         enddo
         write(iw,*) 'in -giah11- ii,jj,xi,xj,qijx,tijx ...'
         write(iw,*)  ii , jj
         write(iw,*) xi,yi,zi
         write(iw,*) xj,yj,zj
         write(iw,*) qijx,qijy,qijz
         write(iw,*) tijx,tijy,tijz
      endif
c
c     ----- i primitive -----
c
      do ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
      cgi=cg(ig)
c
c     ----- j primitive -----
c
      do jg=j1,j2
      aj=ex(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.le.tol) then       
      fac= exp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      cgj=cg(jg)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      ij=0
      do 360 i=mini,maxi
      go to (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),i
  110 dum1=csi*fac
      go to 220
  120 dum1=cpi*fac
      go to 220
  130 dum1=cdi*fac
      go to 220
  140 if(norm) dum1=dum1*sqrt3
      go to 220
  150 dum1=cfi*fac
      go to 220
  160 if(norm) dum1=dum1*sqrt5
      go to 220
  170 if(norm) dum1=dum1*sqrt3
      go to 220
  180 dum1=cgi*fac
      go to 220
  190 if(norm) dum1=dum1*sqrt7
      go to 220
  200 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 220
  210 if(norm) dum1=dum1*sqrt3
  220 continue
c
      do 360 j=minj,maxj 
      go to (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),j
  230 dum2=dum1*csj
      go to 350
  250 dum2=dum1*cpj
      go to 350
  260 dum2=dum1*cdj
      go to 350
  270 if(norm) dum2=dum2*sqrt3
      go to 350
  280 dum2=dum1*cfj
      go to 350
  290 if(norm) dum2=dum2*sqrt5
      go to 350
  300 if(norm) dum2=dum2*sqrt3
      go to 350
  310 dum2=dum1*cgj
      go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      iloc=max(loci+i,locj+j)
      jloc=min(loci+i,locj+j)
      nn=iky(iloc)+jloc
      den0=d0(nn)*two
      ij=ij+1
      cij(ij)=dum2
      dij(ij)=den0
  360 continue          
c
c     ----- -h11- dia + paramagnetic terms -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do kat=1,nat
c
         more=out.and.kat.eq.1
c
         znuc=one
         cx=c(1,kat)
         cy=c(2,kat)
         cz=c(3,kat)
         xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
         yy=xx
         if(nroots.le.3) call rt123
         if(nroots.eq.4) call root4
         if(nroots.eq.5) call root5
         if(nroots.eq.6) then
            call droot
         else
            do iroot=1,nroots
               u9(iroot)=u5(iroot)
               w9(iroot)=w5(iroot)
            enddo
         endif
         do iroot=1,nroots
            uu=u9(iroot)*aa
            ww=w9(iroot)*znuc
            vv=ww
            ww=ww*(uu+uu)
            tt=one/(aa+uu)
            t = sqrt(tt)
            x0=(aax+uu*cx)*tt
            y0=(aay+uu*cy)*tt
            z0=(aaz+uu*cz)*tt
            do j=1,ljtmod
               nj=j
               do i=1,litmod
                  ni=i
                  call dsxyz
c
c     ----- for 1/r -----
c
               if(more) then
                   vx0(i,j,iroot)=xint
                   vy0(i,j,iroot)=yint
                   vz0(i,j,iroot)=zint*vv
               endif
c
c     ----- for x/r**3 -----
c
                   vx1(i,j,iroot)=xint
                   vy1(i,j,iroot)=yint
                   vz1(i,j,iroot)=zint*ww
                  call dvxyz
                  dvx1(i,j,iroot)=xint
                  dvy1(i,j,iroot)=yint
                  dvz1(i,j,iroot)=zint*ww
               enddo
            enddo
         enddo
c
         do iroot=1,nroots
            do j=1,ljt
               do i=1,lit
                    vx1x(i,j,iroot)= vx1(i,j+1,iroot)
                    vy1y(i,j,iroot)= vy1(i,j+1,iroot)
                    vz1z(i,j,iroot)= vz1(i,j+1,iroot)
                   dvx1x(i,j,iroot)=dvx1(i,j+1,iroot)
                   dvy1y(i,j,iroot)=dvy1(i,j+1,iroot)
                   dvz1z(i,j,iroot)=dvz1(i,j+1,iroot)
               enddo
            enddo
            do j=1,ljtmod
               do i=1,lit
                   xvx1 (i,j,iroot)= vx1(i+1,j,iroot)
                   yvy1 (i,j,iroot)= vy1(i+1,j,iroot)
                   zvz1 (i,j,iroot)= vz1(i+1,j,iroot)
                  xdvx1 (i,j,iroot)=dvx1(i+1,j,iroot)
                  ydvy1 (i,j,iroot)=dvy1(i+1,j,iroot)
                  zdvz1 (i,j,iroot)=dvz1(i+1,j,iroot)
               enddo
            enddo
         enddo
         do iroot=1,nroots
c
c     ----- for 1/r derivatives wrt. centers -i- and -j- -----
c
         if(more) then
            call v00xyz(
     1           dvx0(1,1,iroot),dvy0(1,1,iroot),dvz0(1,1,iroot),
     2           vx0d(1,1,iroot),vy0d(1,1,iroot),vz0d(1,1,iroot),
     3            vx0(1,1,iroot), vy0(1,1,iroot), vz0(1,1,iroot),
     4           lit,ljt,ai,aj)
         endif
c
c     ----- d/dx ... operators -----
c
            call h11xyz(
     1            vx1d(1,1,iroot), vy1d(1,1,iroot), vz1d(1,1,iroot),
     2            vx1 (1,1,iroot), vy1 (1,1,iroot), vz1 (1,1,iroot),
     3           lit,ljt,aj)
            call h11xyz(
     1           xvx1d(1,1,iroot),yvy1d(1,1,iroot),zvz1d(1,1,iroot),
     2           xvx1 (1,1,iroot),yvy1 (1,1,iroot),zvz1 (1,1,iroot),
     3           lit,ljt,aj)
         enddo
c
         ij=0
         do i=mini,maxi
            ix=ijx(i)
            iy=ijy(i)
            iz=ijz(i)
            do j=minj,maxj
               jx=ijx(j)
               jy=ijy(j)
               jz=ijz(j)
               dumv  =zero
               transx=zero
               transy=zero
               transz=zero
               dumvx1=zero
               dumvy1=zero
               dumvz1=zero
               dumvx =zero
               dumvy =zero
               dumvz =zero
               dumxyz=zero
               dumxx =zero
               dumxy =zero
               dumxz =zero
               dumyx =zero
               dumyy =zero
               dumyz =zero
               dumzx =zero
               dumzy =zero
               dumzz =zero
               dum1  =zero
               dum2  =zero
               dum3  =zero
               dumx1 =zero
               dumx2 =zero
               dumx3 =zero
               dumy1 =zero
               dumy2 =zero
               dumy3 =zero
               dumz1 =zero
               dumz2 =zero
               dumz3 =zero
               do iroot=1,nroots
c
      if(more) then
c
c     ----- 1/r -----
c
      dumv  =dumv  +
     1    vx0(ix,jx,iroot)* vy0(iy,jy,iroot)* vz0(iz,jz,iroot)
c
c     ----- translation invariance for nuclear derivatives of 1/r -----
c
      transx=transx+
     1   dvx0 (ix,jx,iroot)* vy0 (iy,jy,iroot)* vz0 (iz,jz,iroot)
     2 + dvx1 (ix,jx,iroot)* vy1 (iy,jy,iroot)* vz1 (iz,jz,iroot)
     3 +  vx0d(ix,jx,iroot)* vy0 (iy,jy,iroot)* vz0 (iz,jz,iroot)
      transy=transy+
     1    vx0 (ix,jx,iroot)*dvy0 (iy,jy,iroot)* vz0 (iz,jz,iroot)
     2 +  vx1 (ix,jx,iroot)*dvy1 (iy,jy,iroot)* vz1 (iz,jz,iroot)
     3 +  vx0 (ix,jx,iroot)* vy0d(iy,jy,iroot)* vz0 (iz,jz,iroot)
      transz=transz+
     1    vx0 (ix,jx,iroot)* vy0 (iy,jy,iroot)*dvz0 (iz,jz,iroot)
     2 +  vx1 (ix,jx,iroot)* vy1 (iy,jy,iroot)*dvz1 (iz,jz,iroot)
     3 +  vx0 (ix,jx,iroot)* vy0 (iy,jy,iroot)* vz0d(iz,jz,iroot)
c
c     ----- x/r**3 -----
c
      dumvx1=dumvx1+
     1   dvx1(ix,jx,iroot)* vy1(iy,jy,iroot)* vz1(iz,jz,iroot)
      dumvy1=dumvy1+
     1    vx1(ix,jx,iroot)*dvy1(iy,jy,iroot)* vz1(iz,jz,iroot)
      dumvz1=dumvz1+
     1    vx1(ix,jx,iroot)* vy1(iy,jy,iroot)*dvz1(iz,jz,iroot)
c
      endif
c
c     ----- for h(1,1) diamagnetic -----
c
      dumxx =dumxx +
     1  dvx1x(ix,jx,iroot)* vy1 (iy,jy,iroot)* vz1 (iz,jz,iroot)
      dumxy =dumxy +
     1   vx1x(ix,jx,iroot)*dvy1 (iy,jy,iroot)* vz1 (iz,jz,iroot)
      dumxz =dumxz +
     1   vx1x(ix,jx,iroot)* vy1 (iy,jy,iroot)*dvz1 (iz,jz,iroot)
      dumyx =dumyx +
     1  dvx1 (ix,jx,iroot)* vy1y(iy,jy,iroot)* vz1 (iz,jz,iroot)
      dumyy =dumyy +
     1   vx1 (ix,jx,iroot)*dvy1y(iy,jy,iroot)* vz1 (iz,jz,iroot)
      dumyz =dumyz +
     1   vx1 (ix,jx,iroot)* vy1y(iy,jy,iroot)*dvz1 (iz,jz,iroot)
      dumzx =dumzx +
     1  dvx1 (ix,jx,iroot)* vy1 (iy,jy,iroot)* vz1z(iz,jz,iroot)
      dumzy =dumzy +
     1   vx1 (ix,jx,iroot)*dvy1 (iy,jy,iroot)* vz1z(iz,jz,iroot)
      dumzz =dumzz +
     1   vx1 (ix,jx,iroot)* vy1 (iy,jy,iroot)*dvz1z(iz,jz,iroot)
c
c     ----- for h(1,1) paramagnetic -----
c
      dum1=dum1+
     1    vx1 (ix,jx,iroot)* dvy1 (iy,jy,iroot)*  vz1d(iz,jz,iroot)
     2 -  vx1 (ix,jx,iroot)*  vy1d(iy,jy,iroot)* dvz1 (iz,jz,iroot)
      dum2=dum2+
     1    vx1d(ix,jx,iroot)*  vy1 (iy,jy,iroot)* dvz1 (iz,jz,iroot)
     2 - dvx1 (ix,jx,iroot)*  vy1 (iy,jy,iroot)*  vz1d(iz,jz,iroot)
      dum3=dum3+
     1   dvx1 (ix,jx,iroot)*  vy1d(iy,jy,iroot)*  vz1 (iz,jz,iroot)
     2 -  vx1d(ix,jx,iroot)* dvy1 (iy,jy,iroot)*  vz1 (iz,jz,iroot)
c
      dumx1=dumx1+
     1   xvx1 (ix,jx,iroot)* dvy1 (iy,jy,iroot)*  vz1d(iz,jz,iroot)
     2 - xvx1 (ix,jx,iroot)*  vy1d(iy,jy,iroot)* dvz1 (iz,jz,iroot)
      dumx2=dumx2+
     1   xvx1d(ix,jx,iroot)*  vy1 (iy,jy,iroot)* dvz1 (iz,jz,iroot)
     2 -xdvx1 (ix,jx,iroot)*  vy1 (iy,jy,iroot)*  vz1d(iz,jz,iroot)
      dumx3=dumx3+
     1  xdvx1 (ix,jx,iroot)*  vy1d(iy,jy,iroot)*  vz1 (iz,jz,iroot)
     2 - xvx1d(ix,jx,iroot)* dvy1 (iy,jy,iroot)*  vz1 (iz,jz,iroot)
      dumy1=dumy1+
     1    vx1 (ix,jx,iroot)*ydvy1 (iy,jy,iroot)*  vz1d(iz,jz,iroot)
     2 -  vx1 (ix,jx,iroot)* yvy1d(iy,jy,iroot)* dvz1 (iz,jz,iroot)
      dumy2=dumy2+
     1    vx1d(ix,jx,iroot)* yvy1 (iy,jy,iroot)* dvz1 (iz,jz,iroot)
     2 - dvx1 (ix,jx,iroot)* yvy1 (iy,jy,iroot)*  vz1d(iz,jz,iroot)
      dumy3=dumy3+
     1   dvx1 (ix,jx,iroot)* yvy1d(iy,jy,iroot)*  vz1 (iz,jz,iroot)
     2 -  vx1d(ix,jx,iroot)*ydvy1 (iy,jy,iroot)*  vz1 (iz,jz,iroot)
      dumz1=dumz1+
     1    vx1 (ix,jx,iroot)* dvy1 (iy,jy,iroot)* zvz1d(iz,jz,iroot)
     2 -  vx1 (ix,jx,iroot)*  vy1d(iy,jy,iroot)*zdvz1 (iz,jz,iroot)
      dumz2=dumz2+
     1    vx1d(ix,jx,iroot)*  vy1 (iy,jy,iroot)*zdvz1 (iz,jz,iroot)
     2 - dvx1 (ix,jx,iroot)*  vy1 (iy,jy,iroot)* zvz1d(iz,jz,iroot)
      dumz3=dumz3+
     1   dvx1 (ix,jx,iroot)*  vy1d(iy,jy,iroot)* zvz1 (iz,jz,iroot)
     2 -  vx1d(ix,jx,iroot)* dvy1 (iy,jy,iroot)* zvz1 (iz,jz,iroot)
c
               enddo     
               ij=ij+1
               dum=pi212*aa1*cij(ij)
               dumv  =dumv  *dum
               transx=transx*dum
               transy=transy*dum
               transz=transz*dum
               dumvx =dumvx *dum
               dumvy =dumvy *dum
               dumvz =dumvz *dum
               dumxx =dumxx *dum
               dumxy =dumxy *dum
               dumxz =dumxz *dum
               dumyx =dumyx *dum
               dumyy =dumyy *dum
               dumyz =dumyz *dum
               dumzx =dumzx *dum
               dumzy =dumzy *dum
               dumzz =dumzz *dum
               dum1  =dum1  *dum
               dum2  =dum2  *dum
               dum3  =dum3  *dum
               dumx1 =dumx1 *dum
               dumx2 =dumx2 *dum
               dumx3 =dumx3 *dum
               dumy1 =dumy1 *dum
               dumy2 =dumy2 *dum
               dumy3 =dumy3 *dum
               dumz1 =dumz1 *dum
               dumz2 =dumz2 *dum
               dumz3 =dumz3 *dum
c
            if(more) then
c
c     ----- 1/r and translation invariance -----
c
               v00(1,ij)=v00(1,ij)+dumv
               v00(2,ij)=v00(2,ij)+dumvx
               v00(3,ij)=v00(3,ij)+dumvy
               v00(4,ij)=v00(4,ij)+dumvz
c
            endif
c
c     ----- h(1,1) diamagnetic -----
c
               dumxyz= dumxx +dumyy +dumzz
               tmpxxd=(dumxyz-dumxx)
               tmpxyd=(      -dumxy)
               tmpxzd=(      -dumxz)
               tmpyxd=(      -dumyx)
               tmpyyd=(dumxyz-dumyy)
               tmpyzd=(      -dumyz)
               tmpzxd=(      -dumzx)
               tmpzyd=(      -dumzy)
               tmpzzd=(dumxyz-dumzz)
            if(more) then
               h11d(1,ij)=h11d(1,ij)+tmpxxd
               h11d(2,ij)=h11d(2,ij)+tmpxyd
               h11d(3,ij)=h11d(3,ij)+tmpxzd
               h11d(4,ij)=h11d(4,ij)+tmpyxd
               h11d(5,ij)=h11d(5,ij)+tmpyyd
               h11d(6,ij)=h11d(6,ij)+tmpyzd
               h11d(7,ij)=h11d(7,ij)+tmpzxd
               h11d(8,ij)=h11d(8,ij)+tmpzyd
               h11d(9,ij)=h11d(9,ij)+tmpzzd
            endif
               e11d(1,1,kat)=e11d(1,1,kat)+tmpxxd*dij(ij)
               e11d(1,2,kat)=e11d(1,2,kat)+tmpxyd*dij(ij)
               e11d(1,3,kat)=e11d(1,3,kat)+tmpxzd*dij(ij)
               e11d(2,1,kat)=e11d(2,1,kat)+tmpyxd*dij(ij)
               e11d(2,2,kat)=e11d(2,2,kat)+tmpyyd*dij(ij)
               e11d(2,3,kat)=e11d(2,3,kat)+tmpyzd*dij(ij)
               e11d(3,1,kat)=e11d(3,1,kat)+tmpzxd*dij(ij)
               e11d(3,2,kat)=e11d(3,2,kat)+tmpzyd*dij(ij)
               e11d(3,3,kat)=e11d(3,3,kat)+tmpzzd*dij(ij)
c
c     ----- h(1,1) paramagnetic -----
c
            if(more) then
               h01(1,ij)=h01(1,ij)+dum1
               h01(2,ij)=h01(2,ij)+dum2
               h01(3,ij)=h01(3,ij)+dum3
            endif
c  
c     ----- a bit unclear ... yx ... instead of ... xy ... -----
c
               tmpxxp=(qijx*dum1+tijy*dumz1-tijz*dumy1)
               tmpyxp=(qijx*dum2+tijy*dumz2-tijz*dumy2)
               tmpzxp=(qijx*dum3+tijy*dumz3-tijz*dumy3)
               tmpxyp=(qijy*dum1+tijz*dumx1-tijx*dumz1)
               tmpyyp=(qijy*dum2+tijz*dumx2-tijx*dumz2)
               tmpzyp=(qijy*dum3+tijz*dumx3-tijx*dumz3)
               tmpxzp=(qijz*dum1+tijx*dumy1-tijy*dumx1)
               tmpyzp=(qijz*dum2+tijx*dumy2-tijy*dumx2)
               tmpzzp=(qijz*dum3+tijx*dumy3-tijy*dumx3)
c   
            if(more) then
               h11p(1,ij)=h11p(1,ij)+tmpxxp
               h11p(2,ij)=h11p(2,ij)+tmpxyp
               h11p(3,ij)=h11p(3,ij)+tmpxzp
               h11p(4,ij)=h11p(4,ij)+tmpyxp
               h11p(5,ij)=h11p(5,ij)+tmpyyp
               h11p(6,ij)=h11p(6,ij)+tmpyzp
               h11p(7,ij)=h11p(7,ij)+tmpzxp
               h11p(8,ij)=h11p(8,ij)+tmpzyp
               h11p(9,ij)=h11p(9,ij)+tmpzzp
            endif
               e11p(1,1,kat)=e11p(1,1,kat)+tmpxxp*dij(ij)
               e11p(1,2,kat)=e11p(1,2,kat)+tmpxyp*dij(ij)
               e11p(1,3,kat)=e11p(1,3,kat)+tmpxzp*dij(ij)
               e11p(2,1,kat)=e11p(2,1,kat)+tmpyxp*dij(ij)
               e11p(2,2,kat)=e11p(2,2,kat)+tmpyyp*dij(ij)
               e11p(2,3,kat)=e11p(2,3,kat)+tmpyzp*dij(ij)
               e11p(3,1,kat)=e11p(3,1,kat)+tmpzxp*dij(ij)
               e11p(3,2,kat)=e11p(3,2,kat)+tmpzyp*dij(ij)
               e11p(3,3,kat)=e11p(3,3,kat)+tmpzzp*dij(ij)
c
c     ----- h(1,1) dia + para magnetic -----
c
               tmpxx =tmpxxd+tmpxxp
               tmpxy =tmpxyd+tmpxyp
               tmpxz =tmpxzd+tmpxzp
               tmpyx =tmpyxd+tmpyxp
               tmpyy =tmpyyd+tmpyyp
               tmpyz =tmpyzd+tmpyzp
               tmpzx =tmpzxd+tmpzxp
               tmpzy =tmpzyd+tmpzyp
               tmpzz =tmpzzd+tmpzzp
            if(more) then
               h11t(1,ij)=h11t(1,ij)+tmpxx 
               h11t(2,ij)=h11t(2,ij)+tmpxy 
               h11t(3,ij)=h11t(3,ij)+tmpxz 
               h11t(4,ij)=h11t(4,ij)+tmpyx 
               h11t(5,ij)=h11t(5,ij)+tmpyy 
               h11t(6,ij)=h11t(6,ij)+tmpyz 
               h11t(7,ij)=h11t(7,ij)+tmpzx 
               h11t(8,ij)=h11t(8,ij)+tmpzy 
               h11t(9,ij)=h11t(9,ij)+tmpzz 
c
               if(abs(transx).ge.1.0d-08.or.abs(transy).ge.1.0d-08
     1                                  .or.abs(transz).ge.1.0d-08) then
                  write(iw,9993)
                  call hnderr(3,errmsg)
               endif
            endif
c
            enddo     
         enddo     
c
      enddo      
c
            endif
         enddo      
      enddo      
c
      if(out) then
         write(iw,*) ' in -giah11- ii,jj  = ',
     1                             ii,jj
         do ix=1,4
            ij=0
            do i=mini,maxi
               do j=minj,maxj
                  ij=ij+1
                  write(iw,9997) ix,i,j,v00(ix,ij),cij(ij) 
               enddo
            enddo
         enddo
         do ix=1,3
            ij=0
            do i=mini,maxi
               do j=minj,maxj
                  ij=ij+1
                  write(iw,9992) ix,i,j,h01(ix,ij)
               enddo
            enddo
         enddo
         write(iw,9996)
         do ixy=1,9
            ij=0
            do i=mini,maxi
               do j=minj,maxj
                  ij=ij+1
                  write(iw,9995) ixy,i,j,h11d(ixy,ij),h11p(ixy,ij),
     1                                   h11t(ixy,ij),dij(ij)
               enddo
            enddo
         enddo
c 
         if(out) then
            do kat=1,nat
               write(iw,*) ' in -giah11- ii,jj,kat = ',ii,jj,kat
               call prsq(e11d(1,1,kat),3,3,3)
               call prsq(e11p(1,1,kat),3,3,3)
            enddo
         endif
      endif
c
         enddo      
      enddo      
c
c     ----- print -----
c
      if(some) then
         write(iw,9998) 
         do kat=1,nat
            call prsq(e11d(1,1,kat),3,3,3)
         enddo
         write(iw,9994) 
         do kat=1,nat
            call prsq(e11p(1,1,kat),3,3,3)
         enddo
      endif
c
      return
 9999 format(' in -gih11p- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3,/,' lit,ljt,nroots= ',3i3)
 9998 format(/,' e11d = p00 * h11d = ',
     1       /,' ------------------- ')
 9997 format(' -v00(ix , )- = ',3i5,2f15.10)
 9996 format(/)
 9995 format(' -h11(ixy, )- = ',3i5,4f15.10)
 9994 format(/,' e11p = p00 * h11p = ',
     1       /,' ------------------- ')
 9993 format(' something wrong with translational',
     1       ' invariance in -giah11- . stop. ')
 9992 format(' -h01(ix , )- = ',3i5,f15.10)
      end
      subroutine h11xyz(dx,dy,dz,x,y,z,ni,nj,aj)
      implicit REAL (a-h,o-z)
      dimension  x(6,*), y(6,*), z(6,*)
      dimension dx(6,*),dy(6,*),dz(6,*)
c
      do i=1,ni
         dx(i,1)= (-(aj+aj)*x(i,2))
         dy(i,1)= (-(aj+aj)*y(i,2))
         dz(i,1)= (-(aj+aj)*z(i,2))
      enddo
c
      if(nj.eq.1) return
c
      do j=2,nj
         do i=1,ni
            dx(i,j)= (dble(j-1)*x(i,j-1)-(aj+aj)*x(i,j+1))
            dy(i,j)= (dble(j-1)*y(i,j-1)-(aj+aj)*y(i,j+1))
            dz(i,j)= (dble(j-1)*z(i,j-1)-(aj+aj)*z(i,j+1))
         enddo
      enddo
c
      return
      end
      subroutine v00xyz(dx,dy,dz,xd,yd,zd,x,y,z,ni,nj,ai,aj)
      implicit REAL (a-h,o-z)
      dimension  x(6,*), y(6,*), z(6,*)
      dimension dx(6,*),dy(6,*),dz(6,*)
      dimension xd(6,*),yd(6,*),zd(6,*)
c
c     ----- derivatives with respect to xj ... -----
c
      do i=1,ni
         xd(i,1)=-(-(aj+aj)*x(i,2))
         yd(i,1)=-(-(aj+aj)*y(i,2))
         zd(i,1)=-(-(aj+aj)*z(i,2))
      enddo
      if(nj.gt.1) then  
         do j=2,nj
            do i=1,ni
               xd(i,j)=-(dble(j-1)*x(i,j-1)-(aj+aj)*x(i,j+1))
               yd(i,j)=-(dble(j-1)*y(i,j-1)-(aj+aj)*y(i,j+1))
               zd(i,j)=-(dble(j-1)*z(i,j-1)-(aj+aj)*z(i,j+1))
            enddo
         enddo
      endif
c
c     ----- derivatives with respect to xi ... -----
c
      do j=1,nj
         dx(1,j)=-(-(ai+ai)*x(2,j))
         dy(1,j)=-(-(ai+ai)*y(2,j))
         dz(1,j)=-(-(ai+ai)*z(2,j))
      enddo
      if(ni.gt.1) then  
         do j=1,nj
            do i=2,ni
               dx(i,j)=-(dble(i-1)*x(i-1,j)-(ai+ai)*x(i+1,j))
               dy(i,j)=-(dble(i-1)*y(i-1,j)-(ai+ai)*y(i+1,j))
               dz(i,j)=-(dble(i-1)*z(i-1,j)-(ai+ai)*z(i+1,j))
            enddo
         enddo
      endif
c
      return
      end
      subroutine giakd(kflg)
      implicit REAL (a-h,o-z)
      return
      end
      subroutine giakx(d,f,xx,ix,nintmx,ia,nopk)
      implicit REAL (a-h,o-z)
      return
      end
      subroutine giak(d,f,xx,ix,nintmx,iky,nopk)
      implicit REAL (a-h,o-z)
c
c     ----- -giak- forms a skeleton matrix -----
c                   f=( h* + h )/2
c
c              f(i,j)=(h**(i,j) + h**(j,i))/2
c
c     indices in labels are in standard order_
c      i.ge.j , k.ge.l , (ij).ge.(kl)
c
c     all contributions are made into lower half of
c     skeleton matrix.
c     only off-diagonal elements need be divided by two,
c     to obtain the correct f matrix.
c
      parameter (mxatom=MAXAT)
      character*8 errmsg
      character*8 pck2,pck4
      character*8 dsk4,dsk6
      integer shiftr
      logical index2
      logical dbllab
      common/iofile/ir,iw
      common/sqfile/ijk,ipk
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
      common/hndopt/hndtyp,dsktyp,filtyp,pcktyp
      common/pcklab/labsiz
      dimension d(1),f(1),iky(1)
      dimension xx(nintmx),ix(nintmx)
      dimension errmsg(3)
      data errmsg    /'program ','stop in ','- giak -'/
      data zero,pt5  /0.0d+00,0.5d+00/
      data pck2,pck4 /'-2-label','-4-label'/
      data dsk4,dsk6 /'dsk-p -o','dsk-pk-o'/
c
      shiftr(iarg,ibit)= ishft(iarg,-ibit)
      land(iarg1,iarg2)= iand(iarg1,iarg2)
      mask1(iarg)=2**iarg-1
      mask2(iarg)=2**iarg-1
      iword1(iarg)=land(shiftr(iarg,24),mask1(8))
      iword2(iarg)=land(shiftr(iarg,16),mask1(8))
      iword3(iarg)=land(shiftr(iarg, 8),mask1(8))
      iword4(iarg)=land(shiftr(iarg, 0),mask1(8))
      jword1(iarg)=land(shiftr(iarg,16),mask2(16))
      jword2(iarg)=land(shiftr(iarg, 0),mask2(16))
c
      dbllab=labsiz.eq.2
      num2=(num*(num+1))/2
c
      do 50 m=1,num2
   50    f(m)=zero
c
      if(nopk.ne.1) go to 1000
c
c     ----- integrals are not in supermatrix form (nopk=.true.) -----
c
  100 call pread(ijk,xx,ix,nxx,nintmx)
      if(nxx.eq.0) go to 350
      nint=iabs(nxx)
      if(nint.gt.nintmx) call hnderr(3,errmsg)
      if(dbllab) then
         do 200 m=1,nint
            label1=ix((m-1)*labsiz+1)
            label2=ix((m-1)*labsiz+2)
            i=jword1(label1)
            j=jword2(label1)
            k=jword1(label2)
            l=jword2(label2)
            nij=iky(i)+j
            nkl=iky(k)+l
            nik=iky(i)+k
            nil=iky(i)+l
            njk=iky(max(j,k))+min(j,k)
            njl=iky(max(j,l))+min(j,l)
            val=xx(m)
            val4=(val+val)+(val+val)
            f(nik)=f(nik)+val*d(njl)
            f(nil)=f(nil)+val*d(njk)
            f(njk)=f(njk)+val*d(nil)
            f(njl)=f(njl)+val*d(nik)
  200       continue
      else
         do 300 m=1,nint
            label=ix(m)
            i=iword1(label)
            j=iword2(label)
            k=iword3(label)
            l=iword4(label)
            nij=iky(i)+j
            nkl=iky(k)+l
            nik=iky(i)+k
            nil=iky(i)+l
            njk=iky(max(j,k))+min(j,k)
            njl=iky(max(j,l))+min(j,l)
            val=xx(m)
            val4=(val+val)+(val+val)
            f(nik)=f(nik)+val*d(njl)
            f(nil)=f(nil)+val*d(njk)
            f(njk)=f(njk)+val*d(nil)
            f(njl)=f(njl)+val*d(nik)
  300       continue
      endif
      if(nxx.gt.0) go to 100
  350 continue
      if(num.eq.1) go to 500
      do 400 m=2,num
         do 400 n=1,m-1
            nij=iky(m)+n
  400       f(nij)=f(nij)*pt5
  500 continue
      call rewfil(ijk)
      return
c
c     ----- integrals are in supermatrix form (nopk=.false.) -----
c
 1000 continue
      write(iw,9999)
      call hnderr(3,errmsg)
      return
 9999 format(' integrals may not be in supermatrix form',
     1       ' for the -giao- -k+0- operator. stop')
      end
      subroutine gias(gauge,r00,s,t,num,noc,norb,ndim)
      implicit REAL (a-h,o-z)
      parameter (maxorb=MAXORB)
      logical dbug
      logical out
      common/iofile/ir,iw
      common/filgia/ift1,ift2,ift3,ift4
INCLUDE(../m4/common/mapper)
      dimension r00(ndim,ndim,*)
      dimension   s(ndim,ndim)
      dimension   t(*)
      dimension gauge(3,noc)
      data zero    /0.0d+00/
      data one     /1.0d+00/
      data two     /2.0d+00/
      data cvel    /137.0359895d+00/
c
      dbug=.false.
      out =.true. 
      out =out.or.dbug
c
c     ----- orthonormalization derivative term -----
c
      if(out) then
         write(iw,*) 'orthonormalization derivative term'      
      endif
      call rewfil(ift2)
      do iop=1,3
         if(out) then
            write(iw,*) 'read     rx, ry, rz ... = ',iop   
         endif
         norb2=(norb*(norb+1))/2
         call giard(ift2,t,norb2)
         do iorb=1,norb
            do jorb=1,iorb
               ij=iky(iorb)+jorb
               r00(iorb,jorb,iop)=t(ij)
               r00(jorb,iorb,iop)=t(ij)
            enddo
         enddo
         if(dbug) then
            call prsq(r00(1,1,iop),norb,norb,ndim)
         endif
      enddo
c
c     ----- we are ready ... -----
c
      call rewfil(ift2)
      do ifld=1,3
         if(ifld.eq.1) then
            iop1=2
            iop2=3
            if(out) then
               write(iw,*) 'b(x) perturbation'
            endif
         elseif(ifld.eq.2) then
            iop1=3
            iop2=1
            if(out) then
               write(iw,*) 'b(y) perturbation'
            endif
         elseif(ifld.eq.3) then
            iop1=1
            iop2=3
            if(out) then
               write(iw,*) 'b(z) perturbation'
            endif
         endif 
         do korb=1,norb
            do jorb=1,norb
               s(jorb,korb)=zero
            enddo
         enddo
         velfac=one/(two*cvel)
         do koc=1,norb
            if(koc.le.noc) then
               gxk=gauge(1,koc)
               gyk=gauge(2,koc)
               gzk=gauge(3,koc)
               if(dbug) then
                  write(iw,9999) koc,gxk,gyk,gzk
               endif
            else
               gxk=zero
               gyk=zero
               gzk=zero
            endif
            do joc=1,norb
               if(joc.le.noc) then
                  gxj=gauge(1,joc)
                  gyj=gauge(2,joc)
                  gzj=gauge(3,joc)
                  if(dbug) then
                     write(iw,9999) joc,gxj,gyj,gzj
                  endif
               else
                  gxk=zero
                  gyk=zero
                  gzk=zero
               endif
               gx=gxj-gxk
               gy=gyj-gyk
               gz=gzj-gzk
               if(ifld.eq.1) then
                  c1=gz
                  c2=gy
               elseif(ifld.eq.2) then
                  c1=gx
                  c2=gz
               elseif(ifld.eq.3) then
                  c1=gy
                  c2=gx
               endif
               s(joc,koc)=( c1*r00(joc,koc,iop1)
     1                     -c2*r00(joc,koc,iop2) )*(velfac/two)
            enddo
         enddo
         call giawt(ift2,s,ndim*ndim)
         if(dbug) then
            call prsq(s,norb,norb,ndim)
         endif
      enddo
      call rewfil(ift2)
c
      return
 9999 format(' in -gias- , gauge for -occ- = ',i5,
     1       ' gx,gy,gz = ',3f15.9)
      end
      subroutine giacvg(d,h,u,u0,num,noc,norb,ndim,
     1                  iend,iter,some)
      implicit REAL (a-h,o-z)
c
c     ----- this routine checks for convergence of the cphf -----
c           calculation.  convergence is determined by two
c           criteria: the convergence of the tensor element
c           for the xx, yy, or zz direction and the
c           convergence of the individual elements of the u
c           matrix.
c
      logical dbug
      logical  out
      logical some
      parameter (zero=0.0d+00)
      parameter (mxioda=1000)   
      common/cvggia/cnva,cnvb
      common/pargia/vnew,vold
      common/iofile/ir,iw
      common/hdafile/idafh,navh,ioda(2,mxioda)
      dimension d(ndim,*),h(ndim,*),u(ndim,*),u0(ndim,*)
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      if(dbug) then
         write(iw,*) '-d(1) in -giacvg-'
         call prsq(d,num,num,ndim)
      endif
      if(out) then
         write(iw,*) '-u(1) in -giacvg-'
         call prsq(u,norb,norb,ndim)
      endif
c
      if(iter.gt.0) then    
c
c     ----- iter > 1 -----
c
         vold  = vnew
         vnew  = vold           
         diff  = vold - vnew
         udelm = zero
         uval  = zero
         ival  = 0
         jval  = 0
c
         num3=num*num
         call daread(idafh,ioda,u0,num3,32)
         if(dbug) then
            write(iw,*) '-previous u(1) in -giacvg-'
            call prsq(u0,norb,norb,ndim)
         endif
c
         icon=0
         do j=1,noc
            do i=noc+1,norb
               udel = u(i,j) - u0(i,j)
               if(abs(udel).ge.cnvb) icon = icon + 1
               if(abs(udel).gt.abs(udelm)) then
                  udelm = udel
                  uval = u(i,j)
                  ival = i
                  jval = j
               endif
            enddo
         enddo
c
         num3=num*num
         call dawrit(idafh,ioda,u,num3,32,navh)
c
         if(abs(diff).lt.cnva.and.icon.eq.0) then 
            iend = 1
         else
            iend = 0
         endif
         if(some) then
            write(iw,9999) iter,vnew,diff,udelm,uval,ival,jval,
     1                     icon,iend
            call fiflsh(iw)
         endif
         return
c
      else
c
c     ----- iter = 1 -----
c
         vold  = zero
         vnew  = vold            
         diff  = vold - vnew
         udelm = zero
         uval  = zero
         ival  = 0
         jval  = 0
c
         num3=num*num
         call dawrit(idafh,ioda,u,num3,32,navh)
c
         icon = noc * (norb - noc)
         iend=0
         if(some) then
            write(iw,9999) iter,vnew,diff,udelm,uval,ival,jval,
     1                     icon,iend
            call fiflsh(iw)
         endif
      endif
      return
 9999 format(1x,i3,f20.9,3f16.9,2i4,i10,i3)
      end
      subroutine giak10(d,f,xx,ix,nintmx,ndim,nopk)
      implicit REAL (a-h,o-z)
c
c     ----- this routine forms (num*num) fock matrix. it is
c           assumed that both  d and f matrices are unsymmetric.----
c
c     ----- at the present time  the nopk=.false. option is not
c           included.                                           -----
c     ----- for the 2e integrals indices in labels are in
c           standard order , that is,
c               i.ge.j , k.ge.l , (ij).ge.(kl)
c
      parameter   (mxatom=MAXAT)
      parameter   (zero=0.0d+00)
      character*8 errmsg
      character*8 hndtyp
      character*8 dsktyp
      character*8 filtyp
      character*8 pcktyp
      character*8 disk
      logical     dbug
      logical     out
INCLUDE(../m4/common/filel)
      common/iofile/ir,iw,ip,main
      common/blkin/gin(510),nint
      common/sqfile/ijk,ipk
      common/hndopt/hndtyp,dsktyp,filtyp,pcktyp
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(mxatom),
     1                                         c(3,mxatom)
      dimension d(ndim,*),f(ndim,*),xx(*),ix(*)
      dimension errmsg(3)
      data   disk /'disk io '/
      data errmsg /'program ','stop in ','-giak10-'/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      if(out) then
         write(iw,*) '-d(1)- in -giak10-'
         call prsq(d,num,num,ndim)
      endif
c
      do m2=1,num
         do m1=1,num
            f(m1,m2)=zero
         enddo
      enddo
c
      if(nopk.ne.1) then
c
c     ----- integrals are in supermatrix form (nopk=.false.) -----
c     ----- not implimented in fpl systems at this point, so ----
c
         write (iw,9999)
         call hnderr(3,errmsg)
      endif
      if(dsktyp.eq.disk) then
c
c     ----- integrals are not in supermatrix form (nopk=.true.) -----
c
         do 1000 ifile=1,lfile
            main=lotape(ifile)
            call search(liblk(ifile),main)
            call find(main)
            jblock=llblk(ifile)
 1004       jblock=jblock+1
            call get(gin,nw)
            if (nw.eq.0) goto 1005
            if (jblock.ne.0) call find(main)
            call k10gia(f,d,ndim)
            if (jblock) 1004,1000,1004
 1005       llblk(ifile)=liblk(ifile)-iposun(main)+1
 1000    continue
c
      else
c
c     ----- in-core method -----
c
         call k10gix(f,d,ndim)
      endif
      return
 9999 format(' for -gia- calculations, the integrals may not be',
     1       ' in supermatrix format. stop.',/,
     2       ' set -npkfil- = 1 in -$intgrl- .')
 9998 format(' in-core method not yet implemented')
      end
      subroutine k10gia(f,d,ndim)
      implicit REAL (a-h,o-z)
c
c     CPHF matrix-vector multiplication adapted to GAMESS-UK 
c     in analogy to sgmata [util3.m].
c
      parameter (pt5=0.5d+00)
      logical*1 integ,iii,jjj,kkk,lll
      dimension iii(4),jjj(4),kkk(4),lll(4)
      equivalence (iii(1),i),(jjj(1),j)
      equivalence (kkk(1),k),(lll(1),l)
      common/blkin/gg(340),integ(1360),mword
      dimension f(ndim,*),d(ndim,*)
      data i,j,k,l/4*0/
c
      iword = 1
      do iw=1,mword
         jjj(1) = integ(iword)
         iii(1) = integ(iword+1)
         lll(1) = integ(iword+2)
         kkk(1) = integ(iword+3)
         val = gg(iw)
         valh = val*pt5
         if (i.eq.j) valh=valh*pt5
         if (k.eq.l) valh=valh*pt5
         if (i.eq.k.and.j.eq.l) valh=valh*pt5
         dum    = d(j,l)*valh
         f(i,k) = f(i,k)-dum
         f(k,i) = f(k,i)+dum
         dum    = d(j,k)*valh
         f(i,l) = f(i,l)-dum
         f(l,i) = f(l,i)+dum
         dum    = d(i,l)*valh
         f(j,k) = f(j,k)-dum
         f(k,j) = f(k,j)+dum
         dum    = d(i,k)*valh
         f(j,l) = f(j,l)-dum
         f(l,j) = f(l,j)+dum
         iword = iword+4
      enddo       
      return
      end
      subroutine k10gix(f,d,num)
      implicit REAL (a-h,o-z)
      parameter (pt5=0.5d+00)
      integer shiftr
      common/lcm   /y(1)
      dimension f(num,*),d(num,*)
      dimension iy(1)
      equivalence (iy(1),y(1))
c
      shiftr(iarg,ibit)= ishft(iarg,-ibit)
      land(iarg1,iarg2)= iand(iarg1,iarg2)
      mask1(iarg)=2**iarg-1
      jword1(iarg)=land(shiftr(iarg,16),mask1(16))
      jword2(iarg)=land(shiftr(iarg, 0),mask1(16))
c
      num2=(num *(num +1))/2
      num4=(num2*(num2+1))/2
      ny=0
      ij=0
      do i=1,num
      do j=1,i
         ij=ij+1
         mkl=iy(ij)
         if(mkl.ne.0) then      
         do my=1,mkl
            val = y(ny+my+num2+num4)
            kl  =iy(ny+my+num2     )
            k=jword1(kl)
            l=jword2(kl)
            valh = val*pt5
            f(i,k) = f(i,k)-d(j,l)*valh
            f(i,l) = f(i,l)-d(j,k)*valh
            f(j,k) = f(j,k)-d(i,l)*valh
            f(j,l) = f(j,l)-d(i,k)*valh
            f(k,i) = f(k,i)-d(l,j)*valh
            f(k,j) = f(k,j)-d(l,i)*valh
            f(l,i) = f(l,i)-d(k,j)*valh
            f(l,j) = f(l,j)-d(k,i)*valh
         enddo
         endif
         ny = ny+ij
      enddo
      enddo
      return
      end
      subroutine giag10(h10,f00,d00,ndim,next)
      implicit REAL (a-h,o-z)
c
c     ----- -giao- 2-e integrals -----
c
      parameter (mxioda=1000)
      logical   times
      logical   dbug
      logical   out
INCLUDE(../m4/common/dump3)
      common/iofile/ir,iw,ipu(3),idaf
      dimension h10(ndim*ndim,*)
      dimension f00(ndim*ndim)
      dimension d00(*)
      data zero /0.0d+00/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      times=.false.
      times=times.or.out
c
      if(out) then
         write(iw,*) 'in -giag10-'
         call fiflsh(iw)
      endif
c
      do iop=1,3
         do ij=1,ndim*ndim
            h10(ij,iop)=zero
         enddo
      enddo
      do ij=1,ndim*ndim
         f00(ij)=zero
      enddo
c
c     ----- density matrix -----
c
      ndim2=(ndim*(ndim+1))/2
      call rdedx(d00,ndim2,ibl3pa,idaf)
c
c     ----- -giao- integrals -----
c
      call cpuwal(tim,wtim)
       tim0= tim
      wtim0=wtim
c
      call gia2ei(next)
c
      call cpuwal(tim,wtim)
       tim1= tim
      wtim1=wtim
       tim2e= tim1- tim0
      wtim2e=wtim1-wtim0
      if(times) then
         write(iw,9999) tim2e,wtim2e
      endif
c
      if(out) then
         write(iw,*) 'd00'
         if(dbug) then
            call prtr(d00,ndim)
         endif
         write(iw,*) 'g00'
         if(dbug) then
            call prsq(f00,ndim,ndim,ndim)
         endif
         write(iw,*) 'g10-x'
         if(dbug) then
            call prsq(h10(1,1),ndim,ndim,ndim)
         endif
         write(iw,*) 'g10-y'
         if(dbug) then
            call prsq(h10(1,2),ndim,ndim,ndim)
         endif
         write(iw,*) 'g10-z'
         if(dbug) then
         call prsq(h10(1,3),ndim,ndim,ndim)
         endif
      endif
      call fiflsh(iw)
c
      return
 9999 format(' -giag10- , cpu time = ',f10.2,' wall time = ',f10.2)
      end
