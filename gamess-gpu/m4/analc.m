c *****************************************************
c *****************************************************
c             = analc  =
c ******************************************************
c ******************************************************
**==dma.f
      subroutine dma (q)
c
c
c                      distributed multipole analysis
c
c
c  the charge distribution is analysed by treating the overlap density
c  of each pair of primitives as a multipole expansion about the centre
c  of the overlap distribution. this expansion (which terminates) is
c  then moved to the nearest of a list of expansion sites, giving a
c  non-terminating but rapidly convergent expansion. by default, the
c  expansion sites comprise the nuclei, but sites can be deleted or
c  additional ones added. an additional site at the centre of the
c  molecule is recommended for example if there is no atom there. it is
c  possible to specify the maximum rank of multipole which can be
c  generated at any site; multipoles of higher rank which would have
c  been moved to such limited sites are moved to other other sites
c  instead.
c
c  normally the scf density matrix is analysed, but it is possible to
c  specify other density matrices, generated for example by
c  perturbations arising from correlation or external fields, by
c  casscf or mcscf calculation, or by direct-ci calculations.
c
c
_IF(f90test)
      use junk_dma, onoscf=>notscf
_ENDIF
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/discc)
INCLUDE(common/iofile)
INCLUDE(common/restar)
INCLUDE(common/restrj)
INCLUDE(common/infoa)
INCLUDE(common/scra7)
INCLUDE(common/timez)
INCLUDE(common/atmol3)
INCLUDE(common/unocas)
INCLUDE(common/segm)
INCLUDE(common/mapper)
      common/bufb/
     *  eiga(maxorb),frocca(maxorb),eigb(maxorb),froccb(maxorb)
      common/junkc/zcom1(29),zcomm(19),zbuff(10)
     *,zatom(maxat)
INCLUDE(common/tran)
INCLUDE(common/limy)
_IFN(f90test)
      common /junk/cordox,cordoy,cordoz, nc, lmax, tol,
     &           tshift, mindc, maxdc, oshb, oline ,
     &           kw,ida,onoscf,oldd,ocorr,omos,onucle,ojunk
     & ,ccx(3,maxat), limit(maxat),radius(maxat)
_ENDIF
INCLUDE(common/runlab)
      dimension q(*),ytype(2)
      data ytype/'-a-','-b-'/
      data itvec/3/
      data m990/990/
c
_IFN(f90test)
      nav = lenwrd()
      nwp2 = 2*maxorb + 82 + 60/nav
      nwdma = 5 + 4*maxat + (13+maxat+nav)/nav
      i = 1
      if (iplot.gt.0) i = i + iplot*(lensec(nwp2)+lensec(140))
      call rdchr(zatom,maxat,i,num8)
      call reads(cordox,nwdma,num8)
      if (nc.le.0 .or. nc.gt.maxat)
     +    call caserr('parameter error in dma pre-processor')
_ENDIF
      odens = .false.
      if (ida.gt.0 .and. onoscf) odens = .true.
      newbas = num
      lenb = num*(num+1)/2
      l3 = num*num
      l4 = nc*121
      i10 = 1
      i20 = i10 + lenb
      i30 = i20
      if (zscftp.ne.'uhf' .and. zscftp.ne.'gvb' .and. zscftp.ne.'grhf')
     +    then
         need = i20 + max(l3,l4)
      else
         i30 = i20 + lenb
         need = i20 + max(lenb+l3,l4)
      end if
      write (iwr,6010)
      if (lwordp.lt.need) then
         write (iwr,6020) lwordp , need
         call caserr('insufficient memory for dma analysis')
      end if
      if (.not.odens) then
         if (mina.le.0) then
            call caserr('invalid section specified for input vectors')
         end if
         if(iuno.ne.0) then
          call dmano(q(i20),q(i10),iuno,iunopp,' nos')
          old = otran
          go to 20
         endif
         call secget(mina,itvec,iblkv)
         call getqp(zcomm,zbuff,eiga,frocca,nbas,newbas,ncola,ieiga,
     +              idiffa,maxorb,iblkv)
         old = otran
         write (iwr,6030) ytype(1) , mina , ibl3d , yed(ied3) ,
     +                   (zcomm(7-i),i=1,6) , zbuff
         call rdedx(q(i30),l3,iblkv,idaf)
         if (.not.otran) call tdown(q(i30),ilifq,q(i30),ilifq,ncola)
         if (zscftp.eq.'gvb' .or. zscftp.eq.'grhf') then
            call reden2(q(i30),q(i10),q(i20),frocca,newbas,lenb)
            call vadd(q(i10),1,q(i20),1,q(i10),1,lenb)
         else
            call dmtx(q(i10),q(i30),frocca,iky,num,newbas,newbas)
         end if
         if (zscftp.eq.'uhf') then
            if (minb.le.0) then
             call caserr('invalid section specified for input vectors')
            else
               call secget(minb,itvec,iblkv)
               call getqp(zcomm,zbuff,eigb,froccb,nbas,newbas,ncolb,
     +                    ieigb,idiffb,maxorb,iblkv)
               write (iwr,6030) ytype(2) , minb , ibl3d , yed(ied3) ,
     +                         (zcomm(7-i),i=1,6) , zbuff
               call rdedx(q(i30),l3,iblkv,idaf)
               if (.not.otran) call tdown(q(i30),ilifq,q(i30),ilifq,
     +             ncolb)
               call dmtx(q(i20),q(i30),froccb,iky,num,newbas,newbas)
               call vadd(q(i10),1,q(i20),1,q(i10),1,lenb)
            end if
         end if
      else
         write (iwr,6040) ida
         call secget(ida,m990,iblkv)
         call rdedx(q(i10),lenb,iblkv,idaf)
      end if
20    otran = .true.
      ibl7la= iposun(num8)
      call dma2(q(i20),q(i10),iwr)
      otran = old

_IF(f90test)
      call dealloc_junk_dma
_ENDIF

      return
 6010 format (/40x,37('=')/40x,
     +        'distributed multipole analysis module'/40x,37('=')/)
 6020 format (/10x,'core store analysis'/10x,19('-')//10x,
     +        'main core available = ',i7,' words'/10x,
     +        'main core required  = ',i7,' words')
 6030 format (//1x,a3,'vectors restored from section',i4,
     +        ' of dumpfile starting at block',i6,' of ',
     +        a4//' header block information : '/
     +        ' vectors created under account ',a8/1x,a7,
     +        'vectors created by ',a8,' program at ',a8,' on ',a8,
     +        ' in the job ',a8/' with the title : ',10a8)
 6040 format (/1x,'1-particle density matrix restored from section ',i3)
      end
**==dma2.f
      subroutine dma2(q,densty,iwr)
c
_IF(f90test)
      use junk_dma
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
      character *8 zzzzzz
INCLUDE(common/sizes)
INCLUDE(common/timez)
INCLUDE(common/infoa)
_IFN(f90test)
      character*8 ,name
      common /junkc/ zzzzzz(58),name(maxat)
      logical lshb, linear
      logical lmos,notscf ,ladd
      logical   nuclei, lcorr
      common /junk/ ox, oy, oz, nc, lmax, tol,
     &           tshift, mindc, maxdc, lshb, linear,
     &           kw,ida,notscf,ladd,lcorr,lmos,nuclei,ljunk
     & ,x(3,maxat), limit(maxat),radius(maxat), d(20,20),
     &           qx(121), binom(20,20), rt(20)
     &          ,gx(112), gy(112), gz(112)
_ENDIF
INCLUDE(common/mapper)
INCLUDE(common/runlab)
      character*10 charwall
c
      dimension q(121,*),densty(*)
c      move origin to (ox,oy,oz)
c
c
      write (iwr,6010) ox , oy , oz
      do 20 i = 1 , nat
         c(1,i) = c(1,i) - ox
         x(1,i) = x(1,i) - ox
         c(2,i) = c(2,i) - oy
         x(2,i) = x(2,i) - oy
         c(3,i) = c(3,i) - oz
         x(3,i) = x(3,i) - oz
 20   continue
c
      if (linear) then
         do 30 i = 2 , nat
            linear = linear .and.
     +        ((c(1,i)-c(1,1))**2+(c(2,i)-c(2,1))**2.lt.1.0d-6)
 30      continue
c
         if (.not.linear) write (iwr,6020)
      end if
c
c
      binom(1,1) = 1.0d0
      rt(1) = 1.0d0
      do 50 k = 2 , 20
         rt(k) = dsqrt(dfloat(k))
         binom(k,1) = k
         binom(k-1,k) = 0.0d0
         do 40 m = 2 , k
            binom(k,m) = binom(k-1,m-1) + binom(k-1,m)
 40      continue
 50   continue
      if (.not.(linear)) then
         do 70 k = 2 , 20
            do 60 m = 1 , k
               binom(k,m) = dsqrt(binom(k,m))
 60         continue
 70      continue
c
         lmax = min(lmax,10)
      end if
c
      l = 0
      do 80 i = 1 , nc
         limit(i) = min(lmax,limit(i))
         l = max(l,limit(i))
 80   continue
      lmax = l
c
c
      do 90 n = 1 , nc
         call vclr(q(1,n),1,121)
 90   continue
c
c
_IF(f90test)
      if (linear) call dmaql0(q,densty,kw)
      if (.not.linear) call dmaqlm(q,densty,kw)
_ELSE
      if (linear) call dmaql0(q,densty,kw,nuclei)
      if (.not.linear) call dmaqlm(q,densty,kw,nuclei)
_ENDIF
c
      call vclr(qx,1,121)
c
c
c punch site information
c
       call blkdms(nc, name, x, limit, radius)
c
      do 100 i = 1 , nc
c
c  print results for site i
         write (iwr,6030) i , name(i) , (x(k,i),k=1,3) , limit(i) ,
     +                   radius(i)
         call dmapq(q(1,i),limit(i),linear,iwr)
c
c punch results for site i
c
         call blkdma(q(1,i),limit(i),linear,i)
c
c  shift total for site i to origin
c
         if (.not.linear) call dmasq(q(1,i),0,limit(i),qx,lmax,x(1,i),
     +                               x(2,i),x(3,i),iwr)
         if (linear) call dmasz(q(1,i),0,limit(i),qx,lmax,x(3,i))
c
 100  continue
c
c  print total multipoles
      write (iwr,6040)
      call dmapq(qx,lmax,linear,iwr)
c
c
c     reset origin
c
      do 110 i = 1 , nat
         c(1,i) = c(1,i) + ox
         c(2,i) = c(2,i) + oy
         c(3,i) = c(3,i) + oz
 110  continue
c
c
      cpu = cpulft(1)
      write (iwr,6050) cpu,charwall()
      return
 6010 format (//1x,'the origin of the coordinate frame is moved',' to ',
     +        3f14.8/)
 6020 format (' atoms are not arranged parallel to the z axis'/
     +        ' the general routine will be used')
 6030 format (/1x,'site',i7,' tag : ',a8/1x,22('=')//'   x =',f12.6,
     +        '  y =',f12.6,'  z =',f12.6/'   maximum rank =',i3,
     +        '   relative radius =',f7.3/)
 6040 format (/5x,35('=')/5x,'total multipoles referred to origin'/5x,
     +        35('=')/)
 6050 format (/' end of distributed multipole analysis at ',f9.2,
     +        ' seconds',a10,' wall'/)
      end
**==dmain.f
      subroutine dmain
c
c                      distributed multipole analysis
c                            = input routine =
c
_IF(f90test)
      use junk_dma
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
      character *8 b
      character *4 abuf,a,word,ylabel,char,all,ytrunc,b4
      character *4 ash,bsh,yrad
INCLUDE(common/sizes)
      parameter (max2 = maxorb+maxorb)
INCLUDE(common/timez)
INCLUDE(common/runlab)
INCLUDE(common/restar)
INCLUDE(common/work)
INCLUDE(common/machin)
INCLUDE(common/direc)
c
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/infob)
c
INCLUDE(common/scra7)
_IFN(f90test)
      logical lshb, linear
      logical lmos,notscf ,ladd
      logical   nuclei, lcorr
      common/junk /vlist(400),newpro(204),lmo(max2),
     * ox,oy,oz,
     &              nc, lmax, tol,
     &           tshift, mindc, maxdc, lshb, linear,
     &           kw,ida,notscf,ladd,lcorr,lmos,nuclei,ljunk
     & ,x(3,maxat), limit(maxat),radius(maxat)
      character *8 ztype,name
      common/junkc/ylabel(26),ztype(831),name(maxat)
_ENDIF
c
c
      dimension word(16)
      data  ash /'a'/, bsh /'b'/
      data word / 'add', 'dele', 'gaug',
     &     'limi', 'line', 'gene', 'atom', 'repo',
     &     'corr', 'nonu', 'shif', 'note',
     &     'scf',  'dens', 'mos', 'radi'/
      data limdir,char,all/16,'char','all'/
      data rln10 / 2.30258d0/
c
c  a global option slat or vdw may be specified to change all atomic radii
c  to slater radii, which does give better charges for e.g.HI
c  (may be again overwritten by the radius directive)
c
c
c  directive syntax:
c
c   multipoles
c   \options\
c   start
c
c
c  options may occur in any order and are:
c
c  linear  use a faster version, which is applicable when all the
c          atoms lie in a line parallel to the z axis and only the
c          z components of the multipoles are required.  in this
c          case the maximum rank is 20. this option is revoked if the
c          molecule is found not to be linear. if the molecule is subjec
c          to an external field, or is not in a singlet sigma state, the
c          may be other non-vanishing multipole moments which will not b
c          calculated; however the ql0 will still be correct.
c
c  limit \name\ lmax
c          limit the rank of multipoles on sites with the name given
c          to lmax at most.  (contributions with higher ranks are
c          moved to other sites.) if no name is given the limit applies
c          all sites. default (and maximum) is 20 for the linear
c          version, 10 otherwise.
c
c  atoms   (default) move all contributions to nearest atom.
c
c  delete name
c          delete all sites with the name given.  delete all deletes
c          all sites. delete charge       deletes nuclear charges on
c          atoms in molecule.
c
c  radius name radius
c          specify a relative radius for all sites with the name given.
c          the actual distances from an overlap centre to the sites are
c          scaled by dividing by the relative radii of the sites, and
c          the contributions are moved to the site which is closest, in
c          terms of scaled distances, to the overlap centre. the default
c          is that all sites have relative radius 1.0.
c
c  add name x y z \lmax \radius\\
c          add a new site at (x,y,z) with the name specified.  the
c          multipole rank is limited to lmax if a value is specified,
c          and a relative radius can be specified also.
c
c  report  print the multipole contributions of each pair of primitives
c          as the calculation proceeds.
c
c  corrections
c          use the perturbation correction to the density matrix
c          instead of the scf density matrix. scf can be specified to
c          use the scf density matrix, but this is the default.
c
c  density isec
c          the density matrix is to be read from section isec of the
c          dumpfile.
c
c  mos isecv
c          the density matrix is in the m.o. basis , and the appropriate
c          m.o.'s are in section isecv. default is section nominated
c          on vectors directive
c
c  nonuclear
c          the nuclear contribution is not to be evaluated.
c
c  shift tshift \all\ \a\ \b\
c          distribute multipoles from an overlap contribution around
c          several dma sites, using a gaussian weighting function.
c          'tshift' is a cutoff parameter; maximum 1, minimum 1e-6.
c          1 (default) means distribution is to the nearest sites
c          only. all means that overlaps coincident with dma sites
c          are also distributed. a and b are different methods for
c          doing the distribution.
c
      nav = lenwrd()
c
      nwdma = 5 + 4*maxat + (13+maxat+nav)/nav
      lendma = lensec(maxat) + lensec(nwdma)
c

_IF(f90test)
      call alloc_junk_dma(maxat)
_ENDIF

c   
c...   is there global radii set ??
c
      call inpa4(yrad)
c
      lmax = 20
      ox = 0.0d0
      oy = 0.0d0
      oz = 0.0d0
      kw = 0
      ida = 0
      linear = .false.
      ladd = .false.
      lcorr = .false.
      notscf = .false.
      lmos = .true.
      nuclei = .true.
      tshift = 1.0d0
      lshb = .false.
      tol = rln10*itol
 20   mindc = 1
      maxdc = maxat
      do 40 i = 1 , nat
         do 30 j = 1 , 3
            x(j,i) = c(j,i)
 30      continue
         limit(i) = lmax
         if (yrad.eq.'slat'.or.yrad.eq.'vdw') then
            ii = czanr(i) + 0.01
_IF(ccpdft)
            radius(i) = srad(ii)
_ELSE
            call caserr('for slater radii enable dft module')
_ENDIF
         else
            radius(i) = 1.0d0
         end if
         name(i) = zaname(i)
 40   continue
      nc = nat
c
 50   call input
      call inpa4(a)
      k = locatc(word,limdir,a)
      if (k.ne.0) then
         go to (250,150,60,200,70,80,20,90,100,110,180,50,120,130,140,
     +          230) , k
      else
         k = locatc(ydd(101),limit2,a)
         if (k.ne.0) then
c
c  start
c
            jrec = jrec - 1
c
_IFN(f90test)
            call wrtc(name,maxat,ibl7la,num8)
            call wrt3s(ox,nwdma,num8)
            ibl7la= ibl7la+ lendma
_ENDIF
            return
         else
            call caserr(
     +     'unrecognised dma directive or faulty directive ordering')
         end if
      end if
c    gauge
 60   call inpf(ox)
      call inpf(oy)
      call inpf(oz)
      go to 50
c    linear
 70   linear = .true.
      go to 50
c      general
 80   linear = .false.
      go to 50
c     report
 90   kw = iwr
      go to 50
c     correcti
 100  notscf = .true.
      nuclei = .false.
      lcorr = .true.
      call inpa4(a)
      if (a.eq.word(2)) ladd = .true.
      if (ladd) nuclei = .true.
      go to 50
c     nonuclear
 110  nuclei = .false.
      go to 50
c    scf
 120  lcorr = .false.
      ladd = .false.
      nuclei = .true.
      lmos = .false.
      notscf = .false.
      go to 50
c     density
 130  notscf = .true.
      call inpi(ida)
      go to 50
c     mos
 140  lmos = .true.
      call inpi(isecv)
      go to 50
c
c
c  delete
 150  call inpa(b)
      b4 = ytrunc(b)
      if (b4.eq.all) then
c
c  delete all
         nc = 0
c
      else if (b4.ne.char) then
c
c  delete name
         k = 0
         do 170 i = 1 , nc
            if (b.ne.name(i)) then
               k = k + 1
               if (k.ne.i) then
                  name(k) = name(i)
                  do 160 j = 1 , 3
                     x(j,k) = x(j,i)
 160              continue
                  limit(k) = limit(i)
               end if
            end if
 170     continue
         nc = k
      else
c
c  delete charge
         nuclei = .false.
      end if
      go to 50
c
c  shift parameter
 180  call inpf(tshift)
      if (tshift.gt.1.0d0) tshift = 1.0d0
      if (tshift.lt.9.9d-7) tshift = 1.0d-6
 190  if (jump.le.2) go to 50
      call inpa4(abuf)
      if (abuf.eq.all) tshift = -tshift
      if (abuf.eq.bsh) lshb = .true.
      if (abuf.eq.ash) lshb = .false.
      go to 190
c
c  limit
 200  if (jump.gt.2) then
c
         call inpa(b)
         call inpi(l)
         do 210 i = 1 , nc
            if (name(i).eq.b) limit(i) = l
 210     continue
      else
         call inpi(lmax)
         lmax = min(lmax,20)
         do 220 i = 1 , nc
            limit(i) = lmax
 220     continue
      end if
      go to 50
c
c  radius
 230  call inpa(b)
      call inpf(r)
      do 240 i = 1 , nc
         if (name(i).eq.b) radius(i) = r
 240  continue
      go to 50
c
c  add
 250  call inpa(b)
      nc = nc + 1
      if (nc.gt.maxat) write (iwr,6010) maxat
      if (nc.gt.maxat) call caserr('too many sites for dma analysis')
      name(nc) = b
      do 260 i = 1 , 3
         call inpf(x(i,nc))
 260  continue
      limit(nc) = lmax
      radius(nc) = 1.0d0
      if (jump.gt.5) call inpi(limit(nc))
      if (jump.gt.6) call inpf(radius(nc))
      go to 50
 6010 format (' too many sites -- maximum is',i4)
      end
**==dmamq.f
      subroutine dmamq (qc,xp,yp,zp,iw)
c
c----------------------------------------------------------------- dmamq
c
c
c  move the set of multipoles in /junk/ to the nearest multipole site
c  (but if two or more sites are almost equidistant, move a fraction
c  to each).
c
c
_IF(f90test)
      use junk_dma, q=>qx
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
      integer t1,t2
INCLUDE(common/sizes)
_IFN(f90test)
      logical lshb, linear
      common /junk/ oxyz(3), nc, lmax, tol,
     &           tshift, mindc, maxdc, lshb, linear, kwwww(8),
     &  x(3,maxat), limit(maxat),radius(maxat), d(20,20),
     &           q(121), binom(20,20), rt(20)
     &          ,gx(112), gy(112), gz(112)
      character * 8 zzzzz,name
      common/junkc/zzzzz(58),name(maxat)
_ENDIF
INCLUDE(common/mapper)
      common/bufb/rr(maxat)
c
c
      dimension qc(121,*),m(6)
      data eps /1.0d-6/
c
c
      j = 1
      do 20 i = 1 , nc
         rr(i) = ((xp-x(1,i))**2+(yp-x(2,i))**2+(zp-x(3,i))**2)
     +           /radius(i)**2
         if (limit(i).gt.limit(j)) j = i
 20   continue
c
      lp1sq = (lmax+1)**2
      low = 0
c
 30   k = j
      do 40 i = 1 , nc
         if (rr(i).lt.rr(k) .and. limit(i).ge.low) k = i
 40   continue
c
      t1 = low**2 + 1
      t2 = (limit(k)+1)**2
c
      if (rr(k).gt.1.0d-6) then
c
         n = 1
         m(1) = k
         do 50 i = 1 , nc
            if (rr(i).le.rr(k)+eps .and. i.ne.k .and. limit(i)
     +          .eq.limit(k) .and. limit(i).ge.low) then
               n = n + 1
               m(n) = i
            end if
 50      continue
         if (n.ne.1) then
c
            an = 1.0d0/dfloat(n)
            do 60 k = t1 , t2
               q(k) = an*q(k)
 60         continue
         end if
c
         do 70 i = 1 , n
            k = m(i)
            call dmasq(q,low,limit(k),qc(1,k),lmax,xp-x(1,k),yp-x(2,k),
     +                 zp-x(3,k),iw)
 70      continue
c
         do 80 i = t1 , t2
            q(i) = 0.0d0
 80      continue
         if (limit(k).ge.lmax) return
c
         t1 = t2 + 1
         do 100 i = 1 , n
            k = m(i)
            call dmasq(qc(1,k),limit(k)+1,lmax,q,lmax,x(1,k)-xp,x(2,k)
     +                 -yp,x(3,k)-zp,iw)
            do 90 l = t1 , lp1sq
               qc(l,k) = 0.0d0
 90         continue
 100     continue
      else
         do 110 i = t1 , t2
            qc(i,k) = qc(i,k) + q(i)
            q(i) = 0.0d0
 110     continue
         if (limit(k).ge.lmax) return
      end if
c
      low = limit(k) + 1
      go to 30
c
      end
**==dmamqs.f
      subroutine dmamqs (q,xp,yp,zp,aa,iw)
c
c-----------------------------------------------------------------dmamqs
c
c
c  move the set of multipoles in qx to the multipole sites
c  using a gaussian distribution of the multipoles.
c
_IF(f90test)
      use junk_dma
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
_IFN(f90test)
      logical lshb, linear
      common /junk/ oxyz(3), nc, lmax, tol,
     &           tshift, mindc, maxdc, lshb, linear, kwwww(8),
     &  x(3,maxat), limit(maxat),radius(maxat), d(20,20),
     &           qx(121), binom(20,20), rt(20)
     &          ,gx(112), gy(112), gz(112)
      common/junkc/zzzzz(58),name(maxat)
      character *8 zzzzz,name
_ENDIF
INCLUDE(common/mapper)
c
      common/bufb/dam(maxat),
     *y(121),xt(121),fac(maxat),ex(maxat),rr(maxat),mk(maxat)
c
      dimension q(121,*)
c
      data zero /0.0d0/, eps /1.0d-8/
c
c
c
c   if p coincides with a dma site, the multipoles in qx are
c   added to that site. otherwise the rank k multipole in qx
c   is distributed among the dma sites that will accept rank
c   k multipoles. the fraction of the multipole transferred to
c   a site at distance r from p is proportional to:
c
c          exp (-aa * r**2)
c
c   where aa is the sum of the exponents of the two primitive
c   gaussians whose overlap gives the multipoles in qx. but if
c   the fraction to be transferred is less than tshift, none of
c   the multipole is transferred to that site.
c
c
      l2 = (lmax+1)**2
c
c   if only one site is available, transfer all multipoles to it.
      if (nc.ne.1) then
c
c   initialise and find distances of sites from p.
         k = 0
         tf = dabs(tshift)
         tf = dlog(tf)
         li = 0
         do 20 i = 1 , nc
            rr(i) = ((zp-x(3,i))**2+(yp-x(2,i))**2+(xp-x(1,i))**2)
     +              /radius(i)**2
            if (rr(i).lt.eps) li = i
 20      continue
         if ((li.ne.0) .and. (tshift.ge.zero)) then
c
c   dma site at position: transfer all possible multipoles to it,
c   unless all was given in shift command (tshift negative).
            lm = (limit(li)+1)**2
            do 30 j = 1 , lm
               q(j,li) = q(j,li) + qx(j)
               qx(j) = zero
 30         continue
c   return if all multipoles transferred.
            if (limit(li).eq.lmax) return
            k = limit(li) + 1
         end if
c
c   find nearest site imin that can take multipoles rank k
 40      rrm = 1.0d8
         imin = 0
         do 50 i = 1 , nc
            if ((rr(i).lt.rrm) .and. (limit(i).ge.k)) then
               imin = i
               rrm = rr(i)
            end if
 50      continue
         if (imin.eq.0) return
c
c   find all sites m which i) can take multipoles rank k
c   ii) have a gaussian factor larger than tshift with respect
c   to the nearest site imin. also find the lowest limit
c   lmin on the sites m.
         nm = 0
         lmin = limit(imin)
         do 60 i = 1 , nc
            if (limit(i).ge.k) then
               ext = (rrm-rr(i)+eps)*aa
               if (ext.ge.tf) then
                  nm = nm + 1
                  ex(nm) = ext - eps*aa
                  mk(nm) = i
                  lmin = min(lmin,limit(i))
               end if
            end if
 60      continue
c
c   if only one site m is available, transfer all multipoles
c   up to limit(m)=lmin to this site.
         if (nm.gt.1) then
c
c   more than one site m. if shift a reset lmin to k. (one rank
c   only transferred at a time). run over all sites m, forming
c   gaussian factors.
            if (.not.lshb) lmin = k
            tfac = zero
            do 70 m = 1 , nm
               fac(m) = dexp(ex(m))
               tfac = tfac + fac(m)
 70         continue
c   run over m; multiply qx multipoles ranks k-lmin by factor,
c   transfer site m, using a temporary array. if lmin = lmax,
c   go onto next site: otherwise transfer back all multipoles
c   higher than lmin to qx
            l1 = k*k + 1
            lm = (lmin+1)**2
            do 110 m = 1 , nm
               do 80 j = 1 , l2
                  y(j) = zero
                  xt(j) = zero
 80            continue
               do 90 j = l1 , lm
                  xt(j) = qx(j)*fac(m)/tfac
 90            continue
               i = mk(m)
               call dmasq(xt,k,lmin,y,lmax,xp-x(1,i),yp-x(2,i),
     +                    zp-x(3,i),iw)
               do 100 j = l1 , lm
                  q(j,i) = q(j,i) + y(j)
 100           continue
               if (lmin.ne.lmax) then
                  call dmasq(y,lmin+1,lmax,qx,lmax,x(1,i)-xp,x(2,i)-yp,
     +                       x(3,i)-zp,iw)
               end if
 110        continue
c
c  set transferred multipoles to zero. if all multipoles
c  transferred, return; otherwise reset k to lmin+1 and
c  go to 50.
            do 120 j = l1 , lm
               qx(j) = zero
 120        continue
            if (lmin.eq.lmax) return
            k = lmin + 1
         else
            call dmasq(qx,k,lmin,q(1,imin),lmax,xp-x(1,imin),
     +                 yp-x(2,imin),zp-x(3,imin),iw)
            l1 = k*k + 1
            lm = (lmin+1)**2
            do 130 j = l1 , lm
               qx(j) = zero
 130        continue
c   if all multipoles transferred, return. otherwise shift back
c   higher multipoles, reset k and go to 50.
            if (lmin.eq.lmax) return
            call dmasq(q(1,imin),lmin+1,lmax,qx,lmax,x(1,imin)-xp,
     +                 x(2,imin)-yp,x(3,imin)-zp,iw)
            l1 = (lmin+1)**2 + 1
            do 140 j = l1 , l2
               q(j,imin) = zero
 140        continue
            k = lmin + 1
         end if
         go to 40
      else
         call dmasq(qx,0,lmax,q(1,1),lmax,xp-x(1,1),yp-x(2,1),
     +              zp-x(3,1),iw)
         call vclr(qx,1,l2)
         return
      end if
c
c
      end
**==dmamz.f
      subroutine dmamz (q,p, iw)
c
c----------------------------------------------------------------- dmamz
c
c
c  move the multipole contributions in qx to the nearest site
c
c
_IF(f90test)
      use junk_dma
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
      character *1 space
INCLUDE(common/sizes)
_IFN(f90test)
      logical lshb, linear
      common /junk/ oxyz(3), nc, lmax, tol,
     &           tshift, mindc, maxdc, lshb, linear, kwwww(8),
     &  x(3,maxat), limit(maxat),radius(maxat), d(20,20),
     &           qx(121), binom(20,20), rt(20)
     &          ,gx(112), gy(112), gz(112)
      common/junkc/zzzzz(58),name(maxat)
      character *8  zzzzz,name
_ENDIF
INCLUDE(common/mapper)
      common/bufb/r(maxat)
c
      dimension q(121,*)
      dimension m(2)
      data eps /1.0d-8/, space /' '/
c
c
c  j is one of the sites with the largest value of limit.
      j = 1
      do 20 i = 1 , nc
         r(i) = dabs(x(3,i)-p)/radius(i)
         if (limit(i).gt.limit(j)) j = i
 20   continue
c
c  low is the lowest rank of multipole to be moved at this stage
      low = 0
 30   k = j
      do 40 i = 1 , nc
         if (r(i).lt.r(k) .and. limit(i).ge.low) k = i
 40   continue
c
      n = 1
      m(1) = k
      do 50 i = 1 , nc
         if (r(i).le.r(k)+eps .and. i.ne.k .and. limit(i).ge.low .and.
     +       limit(i).eq.limit(k)) then
            n = 2
            m(2) = i
         end if
 50   continue
c
c  multipoles of ranks low to limit(k) are to be moved at this stage
c
      if (iw.gt.0) write (6,6010) p , low , limit(k) ,
     +                            (space,m(i),x(3,m(i)),i=1,n)
c
      l1 = low + 1
      l2 = limit(k) + 1
      if (n.ne.1) then
         do 60 i = l1 , l2
            qx(i) = 0.5d0*qx(i)
 60      continue
      end if
c
      do 70 i = 1 , n
         k = m(i)
         call dmasz(qx,low,limit(k),q(1,k),lmax,p-x(3,k))
 70   continue
c
      do 80 i = l1 , l2
         qx(i) = 0.0d0
 80   continue
c
c  at this point the new sites carry multipoles of all ranks up to
c  lmax.  shift any of ranks higher than limit(k) back to the original
c  origin.  note that following this process qx carries multipoles of
c  all ranks from limit(k)+1 to lmax.
c
      if (limit(k).eq.lmax) return
c
      do 100 i = 1 , n
         k = m(i)
         call dmasz(q(1,k),l2,lmax,qx,lmax,x(3,k)-p)
         do 90 l = l2 , lmax
            q(l+1,k) = 0.0d0
 90      continue
 100  continue
      low = l2
      go to 30
 6010 format (' from',f7.3,': ranks',i3,' to',i3,' to be moved to site',
     +        a1,i1,' at',f7.3,a1,'and site',i2,' at',f7.3)
c
      end
**==dmamzs.f
      subroutine dmamzs (q,p, aa)
c
c-----------------------------------------------------------------dmamzs
c
c
c  move the set of multipoles in qx to the multipole sites
c  using a gaussian distribution of the multipoles.
c
c
_IF(f90test)
      use junk_dma
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
_IFN(f90test)
      logical lshb, linear
      common /junk/ oxyz(3), nc, lmax, tol,
     &           tshift, mindc, maxdc, lshb, linear, kwwww(8),
     &  x(3,maxat), limit(maxat),radius(maxat), d(20,20),
     &           qx(121), binom(20,20), rt(20)
     &          ,gx(112),gy(112),gz(112)
      character *8 zzzzz,name
      common /junkc/ zzzzz(58),name(maxat)
_ENDIF
INCLUDE(common/mapper)
c
      common/bufb/dam(maxat),
     *y(21),xt(21),fac(maxat),ex(maxat),rr(maxat),mk(maxat)
c
      dimension q(121,*)
c
      data zero /0.0d0/, eps /1.0d-8/
c
c
c
c   if p coincides with a dma site, the multipoles in qx are
c   added to that site. otherwise the rank k multipole in qx
c   is distributed among the dma sites that will accept rank
c   k multipoles. the fraction of the multipole transferred to
c   a site at distance r from p is proportional to:
c
c          exp (-aa * r**2)
c
c   where aa is the sum of the exponents of the two primitive
c   gaussians whose overlap gives the multipoles in qx. but if
c   the fraction to be transferred is less than tshift, none of
c   the multipole is transferred to that site.
c
c
      l2 = lmax + 1
c
c   if only one site is available, transfer all multipoles to it.
      if (nc.ne.1) then
c
c   initialise and find distances of sites from p.
         k = 0
         tf = dabs(tshift)
         tf = dlog(tf)
         li = 0
         do 20 i = 1 , nc
            rr(i) = ((p-x(3,i))/radius(i))**2
            if (rr(i).lt.eps) li = i
 20      continue
         if ((li.ne.0) .and. (tshift.ge.zero)) then
c
c   dma site at position: transfer all possible multipoles to it,
c   unless all was given in shift command (tshift negative).
            lm = limit(li) + 1
            do 30 j = 1 , lm
               q(j,li) = q(j,li) + qx(j)
               qx(j) = zero
 30         continue
c   return if all multipoles transferred.
            if (limit(li).eq.lmax) return
            k = limit(li) + 1
         end if
c
c   find nearest site imin that can take multipoles rank k
 40      rrm = 1.0d8
         imin = 0
         do 50 i = 1 , nc
            if ((rr(i).lt.rrm) .and. (limit(i).ge.k)) then
               imin = i
               rrm = rr(i)
            end if
 50      continue
         if (imin.eq.0) return
c
c   find all sites m which i) can take multipoles rank k
c   ii) have a gaussian factor larger than tshift with respect
c   to the nearest site imin. also find the lowest limit
c   lmin on the sites m.
         nm = 0
         lmin = limit(imin)
         do 60 i = 1 , nc
            if (limit(i).ge.k) then
               ext = (rrm-rr(i)+eps)*aa
               if (ext.ge.tf) then
                  nm = nm + 1
                  ex(nm) = ext - eps*aa
                  mk(nm) = i
                  lmin = min(lmin,limit(i))
               end if
            end if
 60      continue
c
c   if only one site m is available, transfer all multipoles
c   up to limit(m)=lmin to this site.
         if (nm.gt.1) then
c
c   more than one site m. if shift a reset lmin to k. (one rank
c   only transferred at a time). run over all sites m, forming
c   gaussian factors.
            if (.not.lshb) lmin = k
            tfac = zero
            do 70 m = 1 , nm
               fac(m) = dexp(ex(m))
               tfac = tfac + fac(m)
 70         continue
c   run over m; multiply qx multipoles ranks k-lmin by factor,
c   transfer site m, using a temporary array. if lmin = lmax,
c   go onto next site: otherwise transfer back all multipoles
c   higher than lmin to qx
            do 110 m = 1 , nm
               do 80 j = 1 , l2
                  y(j) = zero
                  xt(j) = zero
 80            continue
               do 90 j = k , lmin
                  xt(j+1) = qx(j+1)*fac(m)/tfac
 90            continue
               i = mk(m)
               call dmasz(xt,k,lmin,y,lmax,p-x(3,i))
               do 100 j = k , lmin
                  q(j+1,i) = q(j+1,i) + y(j+1)
 100           continue
               if (lmin.ne.lmax) then
                  call dmasz(y,lmin+1,lmax,qx,lmax,x(3,i)-p)
               end if
 110        continue
c
c  set transferred multipoles to zero. if all multipoles
c  transferred, return; otherwise reset k to lmin+1 and
c  go to 50.
            do 120 j = k , lmin
               qx(j+1) = zero
 120        continue
            if (lmin.eq.lmax) return
            k = lmin + 1
         else
            call dmasz(qx,k,lmin,q(1,imin),lmax,p-x(3,imin))
            do 130 j = k , lmin
               qx(j+1) = zero
 130        continue
c   if all multipoles transferred, return. otherwise shift back
c   higher multipoles, reset k and go to 50.
            if (lmin.eq.lmax) return
            call dmasz(q(1,imin),lmin+1,lmax,qx,lmax,x(3,imin)-p)
            l1 = lmin + 1
            do 140 j = l1 , lmax
               q(j+1,imin) = zero
 140        continue
            k = lmin + 1
         end if
         go to 40
      else
         call dmasz(qx,0,lmax,q(1,1),lmax,p-x(3,1))
         call vclr(qx,1,l2)
         return
      end if
c
      end
**==dmano.f
      subroutine dmano(cr,dens,iuno,iunopp,nos)
      implicit REAL  (a-h,o-z)
      dimension cr(*),dens(*)
      character*(*) nos
      character *4 nam
      character *8 com,tit
INCLUDE(common/sizes)
      parameter (mxorb1=maxorb+1)
c...
c... input space natural orbitals and generate 1-pdm
c...
INCLUDE(common/tran)
      common/blkorbs/ deig(maxorb),docc(mxorb1),nbasis,newbas,
     + ncol,jeig,jocc,ipad
INCLUDE(common/discc)
INCLUDE(common/iofile)
INCLUDE(common/machin)
      common/junkc/com(19),tit(10)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
INCLUDE(common/mapper)
      data m29/29/
      if(iuno.eq.0)
     + call caserr('invalid section for natural orbitals')
      write(iwr,98760)nos,iuno,ibl3d,yed(idaf)
98760 format(/1x,a6,'restored from section',i4,
     *' of dumpfile starting at block',i7,' of ',a4)
      call secget(iuno,3,ibluno)
      call rdchr(com,m29,ibluno,idaf)
      call reads(deig,mach(8),idaf)
      nbsq = ncol * nbasis
      nav = lenwrd()
      call readis(ilifc,mach(9)*nav,idaf)
      call reads(cr(1),nbsq,idaf)
      call ibasgn(maxorb,0,nbasis,ilifq)
c...
c... form 1-pdm
c...
      call tdown(cr(1),ilifq,cr(1),ilifq,ncol)
      if(iunopp.eq.2) then
       write(iwr,98761)nos
98761  format(//1x,104('-')/
     * //45x,a6,'(a.o. basis)'/46x,17('-'))
       call prev(cr(1),docc,ncol,newbas,newbas)
      endif
      lenbas=nbasis*(nbasis+1)/2
      call vclr(dens,1,lenbas)
      do 702 loop=1,ncol
      mix=ilifq(loop)
      kb=0
      frodo=docc(loop)
      do 702 nss=1,nbasis
      fdum = frodo*cr(mix+nss)
      call daxpy(nss,fdum,cr(mix+1),1,dens(kb+1),1)
702   kb=kb+nss
      return
      end
**==dmapq.f
      subroutine dmapq(q,lm, linear, iw)
c
c-----------------------------------------------------------------dmapq
c
c
      implicit REAL  (a-h,o-z),integer  (i-n)
      character *2 ql,qm,lq
      integer p
      logical big, linear
      dimension q(121)
      dimension ql(15), qm(31), p(31)
c
      data ql/'q1','q2','q3','q4','q5','q6','q7','q8','q9',
     &    'qa','qb','qc','qd','qe','qf'/
      data qm /'0 ','1c','1s','2c','2s','3c','3s','4c','4s','5c','5s',
     &    '6c','6s','7c','7s','8c','8s','9c','9s','ac','as',
     &    'bc','bs','cc','cs','dc','ds','ec','es','fc','fs'/
      data lq /'q'/
c
c
      if (linear) then
c
c
         write (iw,6010) q(1) , (lq,k,q(k+1),k=1,lm)
         return
      else
c
         write (iw,6020) q(1)
         k = 1
         do 30 l = 1 , lm
            ll1 = l + l + 1
            n = 0
            qsq = 0.0d0
            do 20 i = 1 , ll1
               qsq = qsq + q(k+i)**2
               if (dabs(q(k+i)).ge.5d-7) then
                  n = n + 1
                  p(n) = i
               end if
 20         continue
c
            qs = dsqrt(qsq)
            big = (qs.ge.1d3)
            if (n.gt.0 .and. big) write (iw,6030) ql(l) , qs ,
     +          (ql(l),qm(p(i)),q(p(i)+k),i=1,n)
            if (n.eq.0 .and. .not.big) write (iw,6040) ql(l) , qs
            if (n.gt.0 .and. .not.big) write (iw,6040) ql(l) , qs ,
     +          (ql(l),qm(p(i)),q(p(i)+k),i=1,n)
            k = k + ll1
 30      continue
         return
      end if
 6010 format ('0',4x,'q0 =',f14.8,4(7x,a1,i1,' =',f14.8)/'0',
     +        5(4x,a1,i1,' =',f14.8,3x)
     +        /10('0',5(3x,a1,i2,' =',1pd14.6,3x)/))
 6020 format (' ',19x,'q00  =',f11.6)
 6030 format (' |',a2,'| =',1pe11.3,6(2x,2a2,' =',1pe11.3)
     +        /(18x,6(2x,2a2,' =',1pe11.3)))
 6040 format (' |',a2,'| =',f11.6,6(2x,2a2,' =',f11.6)
     +        /(18x,6(2x,2a2,' =',f11.6)))
c
      end
**==dmaql0.f
_IF(f90test)
      subroutine dmaql0(qc,densty,iw)
_ELSE
      subroutine dmaql0(qc,densty,iw, nuclei)
_ENDIF
c
c-----------------------------------------------------------------dmaql0
c
c
c  calculate multipole moments, and shift them to the nearest site. in
c  this routine, appropriate for linear molecules, only the moments qlm
c  with m=0 are calculated and stored, and they are held in the order
c  q0, q1, q2, ...
c
c  if iw>0, print details of the progress of the calculation.
c  if nuclei, include the nuclear charges in the calculation.
c
c
c  binom(k,m) contains the binomial coefficient (k).
c                                               (m)
c  rt(k) contains sqrt(k).
c
c
_IF(f90test)
      use junk_dma
_ENDIF

      implicit REAL  (a-h,o-z),integer  (i-n)
      logical ieqj, iieqjj, ls
_IFN(f90test)
      logical nuclei
_ENDIF
c
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/hermit)
INCLUDE(common/wermit)
c
_IFN(f90test)
      logical lshb, linear
       common /junk/ oxyz(3), nc, lmax, tol,
     &           tshift, mindc, maxdc, lshb, linear,
     &           kwwwww(8),
     &  x(3,maxat), limit(maxat),radius(maxat), d(20,20),
     &           q0, q(120),
     &           binom(20,20), rt(20)
     &          ,gx(112), gy(112), gz(112)
      character *8 zzzzz,name
      common/junkc/zzzzz(58),name(maxat)
_ENDIF
INCLUDE(common/mapper)
c
INCLUDE(common/runlab)
INCLUDE(common/infoa)
c
_IF(f90test)
      real(kind=8), pointer :: q0
      real(kind=8), pointer :: q(:)
_ENDIF

      dimension qc(121,*),densty(*)
      dimension ix(20),iy(20),iz(20)
c  addresses of gauss-hermite points and weights.
      dimension mink(9), maxk(9)
      data mink /1,2,4,7,11,16,22,29,37/
      data maxk /1,3,6,10,15,21,28,36,45/
c
c  ix(m) is the power of x in the mth basis function in the list s, px,
c  py, ...
      data ix /0,  1,0,0,  2,0,0,1,1,0,  3,0,0,2,2,1,0,1,0,1/
      data iy /0,  0,1,0,  0,2,0,1,0,1,  0,3,0,1,0,2,2,0,1,1/
      data iz /0,  0,0,1,  0,0,2,0,1,1,  0,0,3,0,1,0,1,2,2,1/
c
      data l1/4/, ll1/7/
c
c
_IF(f90test)
      q0 => qx(1)
      q  => qx(2:)
_ENDIF

      ls = (tshift.lt.0.999999d0)
c
      if (iw.gt.0) write (iw,6010)
c
c  loop over pairs of atoms
      do 420 i = 1 , nat
         zi = c(3,i)
c
         q0 = 0.0d0
         do 20 m = 1 , lmax
            q(m) = 0.0d0
 20      continue
         if (.not.(.not.nuclei .or. i.lt.mindc .or. i.gt.maxdc)) then
            q0 = czan(i)
            if (.not.ls) call dmamz(qc,zi,iw)
            if (ls) call dmamzs(qc,zi,1.0d6)
         end if
c
         do 410 j = 1 , i
c
            ieqj = i.eq.j
c
            zji = zi - c(3,j)
            rr = zji**2
c
c  loop over shells for atom i
c
            do 30 ii = 1 , nshell
               if (katom(ii).eq.i) go to 40
c
c  no basis functions on atom i
 30         continue
            go to 420
c
 40         ii1 = ii
            do 50 ii = ii1 , nshell
               if (katom(ii).ne.i) go to 60
 50         continue
            ii = nshell + 1
 60         ii2 = ii - 1
c
            do 400 ii = ii1 , ii2
               i1 = kstart(ii)
               i2 = i1 + kng(ii) - 1
               lit = ktype(ii)
               mini = kmin(ii)
               maxi = kmax(ii)
               loci = kloc(ii) - mini
c
c  loop over shells for atom j
c
               do 70 jj = 1 , nshell
                  if (katom(jj).eq.j) go to 80
c
c  no basis functions on atom j
 70            continue
               go to 410
c
 80            jj1 = jj
               do 90 jj = jj1 , nshell
                  if (katom(jj).ne.j) go to 100
 90            continue
               jj = nshell + 1
 100           jj2 = jj - 1
c
               if (ieqj) jj2 = ii
               do 390 jj = jj1 , jj2
                  j1 = kstart(jj)
                  j2 = j1 + kng(jj) - 1
                  ljt = ktype(jj)
                  minj = kmin(jj)
                  maxj = kmax(jj)
                  locj = kloc(jj) - minj
                  iieqjj = ii.eq.jj
c
c  set up temporary density matrix for this pair of shells
c
                  call vclr(d,1,400)
                  do 130 ib = mini , maxi
                     m = iky(loci+ib)
                     if (iieqjj) then
                        do 110 jb = minj , ib
                           d(ib,jb) = densty(m+locj+jb)
                           d(jb,ib) = densty(m+locj+jb)
 110                    continue
                     else
                        do 120 jb = minj , maxj
                           d(ib,jb) = densty(m+locj+jb)
 120                    continue
                     end if
 130              continue
c
c  insert factors of sqrt(3) for xy, xz and yz functions, if present,
c  and factors of sqrt(5) and sqrt(15) for f functions.
                  if (maxi.ge.8) then
                     do 150 ib = 8 , 10
                        do 140 jb = minj , maxj
                           d(ib,jb) = rt(3)*d(ib,jb)
 140                    continue
 150                 continue
                     if (maxi.gt.10) then
                        do 170 jb = minj , maxj
                           do 160 ib = 14 , 19
                              d(ib,jb) = rt(5)*d(ib,jb)
 160                       continue
                           d(20,jb) = rt(15)*d(20,jb)
 170                    continue
                     end if
                  end if
c
                  if (maxj.ge.8) then
                     do 190 ib = mini , maxi
                        do 180 jb = 8 , 10
                           d(ib,jb) = rt(3)*d(ib,jb)
 180                    continue
 190                 continue
                     if (maxj.gt.10) then
                        do 210 ib = mini , maxi
                           do 200 jb = 14 , 19
                              d(ib,jb) = rt(5)*d(ib,jb)
 200                       continue
                           d(ib,20) = rt(15)*d(ib,20)
 210                    continue
                     end if
                  end if
c
c
c  i primitive
c
                  jgmax = j2
                  do 380 ig = i1 , i2
                     ai = ex(ig)
                     arri = ai*rr
c
c
c  j primitive
c
                     if (iieqjj) jgmax = ig
                     do 370 jg = j1 , jgmax
                        aj = ex(jg)
                        aa = ai + aj
                        dum = aj*arri/aa
                        if (dum.le.tol) then
                           fac = dexp(-dum)
c  introduce factor of 2 if (a) shells are different (ii.ne.jj) because
c  loops over atoms and shells count each pair only once, and (b) if
c  ii.eq.jj but ig.ne.jg, because loops over primitives count each pair
c  only once.  in fact ig.ne.jg covers both cases.  however, different
c  atoms may use the same shells and primitives if there is symmetry, an
c  there is always a factor of 2 if the atoms are different.
                           if (ig.ne.jg .or. .not.ieqj) fac = 2.0d0*fac
c
c  vectors a and b are the positions of atoms i and j relative to the
c  centre p of the overlap distribution.
                           p = aj/aa
                           za = p*zji
                           zb = za - zji
                           zp = zi - za
c
c
                           t = dsqrt(1.0d0/aa)
c
c  nq-1 is the maximum rank of multipole to which these functions
c  contribute. the integrals involve polynomials up to order 2(nq-1),
c  for which nq integration points are required.
                           nq = lit + ljt - 1
                           k1 = mink(nq)
                           k2 = maxk(nq)
                           nq2 = 2*nq
c
c
c  clear integral arrays
                           call vclr(gx,1,112)
                           call vclr(gy,1,112)
                           call vclr(gz,1,112)
                           q0 = 0.0d0
                           call vclr(q,1,lmax)
c
c  the following loop runs through the integration points, accumulating
c  in gz(ll1*l1*(ia-1)+l1*(ib-1)+iq) the quantity
c
c       sum(k) gk * (zk-za)**(ia-1) * (zk-zb)**(ib-1) * zk**(iq-1)
c
c  where gk=w(k)/sqrt(aa) and zk=h(k)/sqrt(aa), and l1=l+1 and
c  ll1=2l+1, where l is the maximum angular momentum to be handled,
c  i.e. 3 at present.
c
                           do 260 k = k1 , k2
                              s = h(k)*t
                              g = w(k)*t
                              zas = s - za
                              zbs = s - zb
c
                              ma = 0
                              paz = g
c
                              do 240 ia = 1 , lit
c
                                 mb = ma
                                 pz = paz
c
                                 do 230 ib = 1 , ljt
c
                                    pq = pz
c
                                    do 220 iq = 1 , nq
                                       gz(mb+iq) = gz(mb+iq) + pq
                                       pq = pq*s
 220                                continue
c
                                    mb = mb + ll1
                                    pz = pz*zbs
 230                             continue
c
                                 ma = ma + l1*ll1
                                 paz = paz*zas
 240                          continue
c
c  in gx it is only necessary to accumulate the sums of even powers of
c  xk:
c
c           sum(k) gk * xk**(iq-1),  iq odd,
c
c  since xa=xb=0. the same values serve for gy.
c
                              ps = g
                              do 250 iq = 1 , nq2 , 2
                                 gx(iq) = gx(iq) + ps
                                 ps = ps*(s**2)
 250                          continue
c
 260                       continue
c
                           do 270 iq = 1 , nq2 , 2
                              gy(iq) = gx(iq)
 270                       continue
c
c
c  now these basic integrals are used to construct the multipole moments
c  for the overlap density corresponding to each pair of basis functions
c  in the pair of shells.
c
                           ci = cs(ig)
                           do 360 ia = mini , maxi
                              if (ia.gt.1) ci = cp(ig)
                              if (ia.gt.4) ci = cd(ig)
                              if (ia.gt.10) ci = cf(ig)
c
                              cj = cs(jg)
                              do 350 jb = minj , maxj
                                 if (jb.gt.1) cj = cp(jg)
                                 if (jb.gt.4) cj = cd(jg)
                                 if (jb.gt.10) cj = cf(jg)
c
                                 f = -fac*ci*cj*d(ia,jb)
c
c
c  the integral of (x**i)*(y**j)*(z**k) over the current pair of atomic
c  primitive orbitals is f*gx(mx+i)*gy(my+j)*gz(mz+k), and is zero
c  unless both mx and my are odd (indexing even powers of x and y).
c
                                 mx = ix(ia) + ix(jb) + 1
                                 my = iy(ia) + iy(jb) + 1
                                 if (mod(mx,2).eq.0 .or. mod(my,2).eq.0)
     +                               go to 350
                                 mz = (iz(ia)*l1+iz(jb))*ll1 + 1
c
c
                                 go to (340,330,320,310,300,290,280) ,
     +                                  nq
c
c
 280                             q(6) = q(6)
     +                                  + 0.0625d0*f*(16.0d0*gx(mx)*gy
     +                                  (my)*gz(mz+6)
     +                                  -120.0d0*(gx(mx+2)*gy(my)+gx(mx)
     +                                  *gy(my+2))*gz(mz+4)
     +                                  +90.0d0*(gx(mx+4)*gy(my)
     +                                  +2.0d0*gx(mx+2)*gy(my+2)+gx(mx)
     +                                  *gy(my+4))*gz(mz+2)
     +                                  -(5.0d0*gx(mx+6)*gy(my)
     +                                  +15.0d0*gx(mx+4)*gy(my+2)
     +                                  +15.0d0*gx(mx+2)*gy(my+4)
     +                                  +5.0d0*gx(mx)*gy(my+6))*gz(mz))
c
c
 290                             q(5) = q(5)
     +                                  + 0.125d0*f*(8.0d0*gz(mz+5)*gx
     +                                  (mx)*gy(my)-40.0d0*gz(mz+3)
     +                                  *(gx(mx+2)*gy(my)+gx(mx)
     +                                  *gy(my+2))+15.0d0*gz(mz+1)
     +                                  *(gx(mx+4)*gy(my)+2.0d0*gx(mx+2)
     +                                  *gy(my+2)+gx(mx)*gy(my+4)))
c
c
 300                             q(4) = q(4)
     +                                  + 0.125d0*f*(8.0d0*gx(mx)*gy(my)
     +                                  *gz(mz+4)
     +                                  -24.0d0*(gx(mx+2)*gy(my)+gx(mx)
     +                                  *gy(my+2))*gz(mz+2)
     +                                  +3.0d0*(gx(mx+4)*gy(my)+gx(mx)
     +                                  *gy(my+4))*gz(mz)
     +                                  +6.0d0*(gx(mx+2)*gy(my+2)*gz(mz)
     +                                  ))
c
c
 310                             q(3) = q(3)
     +                                  + f*(gx(mx)*gy(my)*gz(mz+3)-
     +                                  1.5d0*(gx(mx+2)*gy(my)+gx(mx)
     +                                  *gy(my+2))*gz(mz+1))
c
c
 320                             q(2) = q(2)
     +                                  + f*(gx(mx)*gy(my)*gz(mz+2)-
     +                                  0.5d0*(gx(mx+2)*gy(my)+gx(mx)
     +                                  *gy(my+2))*gz(mz))
c
c
 330                             q(1) = q(1) + f*gx(mx)*gy(my)*gz(mz+1)
c
c
 340                             q0 = q0 + f*gx(mx)*gy(my)*gz(mz)
c
c  end of loop over basis functions
 350                          continue
 360                       continue
c
                           if (iw.gt.0) write (iw,6020) i , j , ii ,
     +                         jj , ig , jg , zp , q0 , (q(iq),iq=1,6)
c
c  move multipoles to expansion centre nearest to overlap centre p.
                           if (.not.ls) call dmamz(qc,zp,iw)
                           if (ls) call dmamzs(qc,zp,aa)
                        end if
c
c  end of loop over primitives
 370                 continue
 380              continue
c
                  if (iw.gt.0) write (iw,6030)
c
c  end of loop over shells
 390           continue
 400        continue
c
c  end of loop over atoms
 410     continue
 420  continue
c
      return
 6010 format ('     atoms   shells primitives  position',
     +        '   multipole contributions ...'/)
 6020 format (1x,3(i5,i4),f11.4,3x,7f12.8)
 6030 format (1x)
      end
**==dmaqlm.f
_IF(f90test)
      subroutine dmaqlm(qc,densty,iw)
_ELSE
      subroutine dmaqlm(qc,densty,iw,nuclei)
_ENDIF
c
c-----------------------------------------------------------------dmaqlm
c
c
c  calculate multipole moments, and shift them to the nearest site.
c  normal routine, for use with non-linear molecules.
c
c
c  if iw>0, print details of the progress of the calculation.
c  if nuclei, include the nuclear charges in the calculation.
c
c
c  the multipoles qlmc and qlms are sqrt(2) times the real and
c  imaginary parts of the complex multipole \ql,-m\*=(-1)**m qlm,
c  except for m=0.
c
c
c  qlm=<rlm*p>, where p is the charge density operator and the rlm
c  are
c           r00  = 1
c
c           r10  = z
c           r11c = x
c           r11s = y
c
c           r20  = (1/2)*(3*z**2-r**2)
c           r21c = sqrt(3)*x*z
c           r21s = sqrt(3)*y*z
c           r22c = sqrt(3/4)*(x**2-y**2)
c           r22s = sqrt(3)*x*y
c
c           etc.
c
c  the multipoles are stored in the order q00, q10, q11c, q11s, q20, ...
c
c
c  the program calculates multipoles up to rank 6 for each pair of
c  primitives. this is enough to cope with s, p, d and f functions.
c  multipoles up to rank lmax are retained when shifting to a new
c  origin. the maximum value allowed for lmax is currently 10.
c
c  the basis functions are:
c
c               1   2   3   4   5   6   7   8   9  10
c               s   px  py  pz dxx dyy dzz dxy dxz dyz
c
c               11  12  13  14  15  16  17  18  19  20
c               xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz
c
c
c  note that integrals involving dxy, dxz and dyz must be multiplied by
c  sqrt(3), because of the difference in normalization from dxx etc.
c  similarly those involving the f functions xxy,xxz,xyy,yyz,xzz and yzz
c  must be multiplied by sqrt(5) and those involving xyz by sqrt(15).
c
c  binom(k,m) contains the square root of the binomial coefficient (k).
c                                                                  (m)
c  rt(k) contains sqrt(k) .
c
c
_IF(f90test)
      use junk_dma, qt=>qx
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
      logical ieqj, iieqjj, ls
_IFN(f90test)
      logical nuclei
_ENDIF
c
INCLUDE(common/sizes)
c
INCLUDE(common/nshel)
c
c  integration points and weights
INCLUDE(common/hermit)
INCLUDE(common/wermit)
c
_IFN(f90test)
      logical lshb, linear
      common /junk/ oxyz(3), nc, lmax, tol,
     &           tshift, mindc, maxdc, lshb, linear,
     &           kwwww(8),
     &  x(3,maxat), limit(maxat),radius(maxat), d(20,20),
     &           q00, q10,q11c,q11s, q20,q21c,q21s,q22c,q22s,
     &           q30,q31c,q31s,q32c,q32s,q33c,q33s,
     &           q40,q41c,q41s,q42c,q42s,q43c,q43s,q44c,q44s,
     &           q50,q51c,q51s,q52c,q52s,q53c,q53s,q54c,q54s,q55c,q55s,
     &           q60,q61c,q61s,q62c,q62s,q63c,q63s,q64c,q64s,q65c,q65s,
     &                       q66c,q66s,
     &           dummy(72), binom(20,20), rt(20)
     &          ,gx(112), gy(112), gz(112)
      character *8 zzzzz,name
      common/junkc/zzzzz(58),name(maxat)
_ENDIF
INCLUDE(common/mapper)
INCLUDE(common/runlab)
INCLUDE(common/infoa)
      dimension qc(121,*),densty(*)
      dimension mink(9), maxk(9)
      dimension ix(20),iy(20),iz(20)
_IFN(f90test)
      dimension qt(121)
      equivalence (q00,qt(1))
_ENDIF
c
c  addresses of gauss-hermite points and weights.
      data mink /1,2,4,7,11,16,22,29,37/
      data maxk /1,3,6,10,15,21,28,36,45/
c
c  ix(m) is the power of x in the mth basis function in the list s, px,
c  py, ...
      data ix /0,  1,0,0,  2,0,0,1,1,0,  3,0,0,2,2,1,0,1,0,1/
      data iy /0,  0,1,0,  0,2,0,1,0,1,  0,3,0,1,0,2,2,0,1,1/
      data iz /0,  0,0,1,  0,0,2,0,1,1,  0,0,3,0,1,0,1,2,2,1/
c
      data l1/4/, ll1/7/
c
c
      ls = (tshift.lt.0.999999d0)
      lp1sq = (lmax+1)**2
c
      if (iw.gt.0) write (iw,6010)
c
c  loop over pairs of atoms
      do 390 i = 1 , nat
         xi = c(1,i)
         yi = c(2,i)
         zi = c(3,i)
         call vclr(qt,1,lp1sq)

         if (.not.(.not.nuclei .or. i.lt.mindc .or. i.gt.maxdc)) then
            q00 = czan(i)
            if (.not.ls) call dmamq(qc,xi,yi,zi,iw)
            if (ls) call dmamqs(qc,xi,yi,zi,1.0d6,iw)
         end if
c
         do 380 j = 1 , i
c
            ieqj = i.eq.j
c
            xji = xi - c(1,j)
            yji = yi - c(2,j)
            zji = zi - c(3,j)
            rr = xji**2 + yji**2 + zji**2
c
c  loop over shells for atom i
c
            do 20 ii = 1 , nshell
               if (katom(ii).eq.i) go to 30
c
c  no basis functions on atom i
 20         continue
            go to 390
c
 30         ii1 = ii
            do 40 ii = ii1 , nshell
               if (katom(ii).ne.i) go to 50
 40         continue
            ii = nshell + 1
 50         ii2 = ii - 1
c
            do 370 ii = ii1 , ii2
               i1 = kstart(ii)
               i2 = i1 + kng(ii) - 1
               lit = ktype(ii)
               mini = kmin(ii)
               maxi = kmax(ii)
               loci = kloc(ii) - mini
c
c  loop over shells for atom j
c
               do 60 jj = 1 , nshell
                  if (katom(jj).eq.j) go to 70
c
c  no basis functions on atom j
 60            continue
               go to 380
c
 70            jj1 = jj
               do 80 jj = jj1 , nshell
                  if (katom(jj).ne.j) go to 90
 80            continue
               jj = nshell + 1
 90            jj2 = jj - 1
c
               if (ieqj) jj2 = ii
               do 360 jj = jj1 , jj2
                  j1 = kstart(jj)
                  j2 = j1 + kng(jj) - 1
                  ljt = ktype(jj)
                  minj = kmin(jj)
                  maxj = kmax(jj)
                  locj = kloc(jj) - minj
                  iieqjj = ii.eq.jj
c
c  set up temporary density matrix for this pair of shells
c
                  call vclr(d,1,400)
                  do 120 ib = mini , maxi
                     m = iky(loci+ib)
                     if (iieqjj) then
                        do 100 jb = minj , ib
                           d(ib,jb) = densty(m+locj+jb)
                           d(jb,ib) = densty(m+locj+jb)
 100                    continue
                     else
                        do 110 jb = minj , maxj
                           d(ib,jb) = densty(m+locj+jb)
 110                    continue
                     end if
 120              continue
c
c  insert factors of sqrt(3) for xy, xz and yz functions, if present,
c  and factors of sqrt(5) and sqrt(15) for f functions.
                  if (maxi.ge.8) then
                     do 140 ib = 8 , 10
                        do 130 jb = minj , maxj
                           d(ib,jb) = rt(3)*d(ib,jb)
 130                    continue
 140                 continue
                     if (maxi.gt.10) then
                        do 160 jb = minj , maxj
                           do 150 ib = 14 , 19
                              d(ib,jb) = rt(5)*d(ib,jb)
 150                       continue
                           d(20,jb) = rt(15)*d(20,jb)
 160                    continue
                     end if
                  end if
c
                  if (maxj.ge.8) then
                     do 180 ib = mini , maxi
                        do 170 jb = 8 , 10
                           d(ib,jb) = rt(3)*d(ib,jb)
 170                    continue
 180                 continue
                     if (maxj.gt.10) then
                        do 200 ib = mini , maxi
                           do 190 jb = 14 , 19
                              d(ib,jb) = rt(5)*d(ib,jb)
 190                       continue
                           d(ib,20) = rt(15)*d(ib,20)
 200                    continue
                     end if
                  end if
c
c
c  i primitive
c
                  jgmax = j2
                  do 350 ig = i1 , i2
                     ai = ex(ig)
                     arri = ai*rr
c
c
c  j primitive
c
                     if (iieqjj) jgmax = ig
                     do 340 jg = j1 , jgmax
                        aj = ex(jg)
                        aa = ai + aj
                        dum = aj*arri/aa
                        if (dum.le.tol) then
                           fac = dexp(-dum)
c  introduce factor of 2 if (a) shells are different (ii.ne.jj) because
c  loops over atoms and shells count each pair only once, and (b) if
c  ii.eq.jj but ig.ne.jg, because loops over primitives count each pair
c  only once.  in fact ig.ne.jg covers both cases.  however, different
c  atoms may use the same shells and primitives if there is symmetry, an
c  there is always a factor of 2 if the atoms are different.
                           if (ig.ne.jg .or. .not.ieqj) fac = 2.0d0*fac
c
c  vectors a and b are the positions of atoms i and j relative to the
c  centre p of the overlap distribution.
                           p = aj/aa
                           xa = p*xji
                           ya = p*yji
                           za = p*zji
                           xb = xa - xji
                           yb = ya - yji
                           zb = za - zji
                           xp = xi - xa
                           yp = yi - ya
                           zp = zi - za
c
c
                           t = 1.0d0/dsqrt(aa)
c
                           if (iw.gt.0) write (iw,6020) i , j , ii ,
     +                         jj , ig , jg , xp , yp , zp
c
c  nq-1 is the maximum rank of multipole to which these functions
c  contribute. the integrals involve polynomials up to order 2(nq-1),
c  for which nq integration points are required.
                           nq = lit + ljt - 1
                           k1 = mink(nq)
                           k2 = maxk(nq)
c
c
c  clear integral arrays
                           call vclr(gx,1,112)
                           call vclr(gy,1,112)
                           call vclr(gz,1,112)
                           call vclr(qt,1,lp1sq)
c
c  the following loop runs through the integration points, accumulating
c  in gx(ll1*l1*(ia-1)+l1*(ib-1)+iq) the quantity
c
c       sum(k) gk * (xk-xa)**(ia-1) * (xk-xb)**(ib-1) * xk**(iq-1)
c
c  where gk=w(k)/sqrt(aa) and xk=h(k)/sqrt(aa), and l1=l+1 and
c  ll1=2l+1, where l is the maximum angular momentum to be handled,
c  i.e. 3 at present.
c  similar expressions for y and z integrals are formed in gy and gz.
c
                           do 240 k = k1 , k2
                              s = h(k)*t
                              g = w(k)*t
                              xas = s - xa
                              yas = s - ya
                              zas = s - za
                              xbs = s - xb
                              ybs = s - yb
                              zbs = s - zb
c
                              ma = 0
                              pax = g
                              pay = g
                              paz = g
c
                              do 230 ia = 1 , lit
c
                                 mb = ma
                                 px = pax
                                 py = pay
                                 pz = paz
c
                                 do 220 ib = 1 , ljt
c
                                    pq = 1.0d0
c
                                    do 210 iq = 1 , nq
                                       gx(mb+iq) = gx(mb+iq) + px*pq
                                       gy(mb+iq) = gy(mb+iq) + py*pq
                                       gz(mb+iq) = gz(mb+iq) + pz*pq
                                       pq = pq*s
 210                                continue
c
                                    mb = mb + ll1
                                    px = px*xbs
                                    py = py*ybs
                                    pz = pz*zbs
 220                             continue
c
                                 ma = ma + l1*ll1
                                 pax = pax*xas
                                 pay = pay*yas
                                 paz = paz*zas
 230                          continue
c
 240                       continue
c
c
c  now these basic integrals are used to construct the multipole moments
c  for the overlap density corresponding to each pair of basis functions
c  in the pair of shells.
c
                           ci = cs(ig)
                           do 330 ia = mini , maxi
                              if (ia.gt.1) ci = cp(ig)
                              if (ia.gt.4) ci = cd(ig)
                              if (ia.gt.10) ci = cf(ig)
c
                              cj = cs(jg)
                              do 320 jb = minj , maxj
                                 if (jb.gt.1) cj = cp(jg)
                                 if (jb.gt.4) cj = cd(jg)
                                 if (jb.gt.10) cj = cf(jg)
c
                                 f = -fac*ci*cj*d(ia,jb)
c
                                 mx = (ix(ia)*l1+ix(jb))*ll1 + 1
                                 my = (iy(ia)*l1+iy(jb))*ll1 + 1
                                 mz = (iz(ia)*l1+iz(jb))*ll1 + 1
c
c  now the integral of (x**i)*(y**j)*(z**k) over the current pair of
c  atomic primitive orbitals is f*gx(mx+i)*gy(my+j)*gz(mz+k).
c
                                 go to (310,300,290,280,270,260,250) ,
     +                                  nq
c
c
c  q6 terms
 250                             q60 = q60 +
     +                                 0.0625d0*f*(16.0d0*gx(mx)*gy(my)
     +                                 *gz(mz+6)
     +                                 -120.0d0*(gx(mx+2)*gy(my)+gx(mx)
     +                                 *gy(my+2))*gz(mz+4)
     +                                 +90.0d0*(gx(mx+4)*gy(my)
     +                                 +2.0d0*gx(mx+2)*gy(my+2)+gx(mx)
     +                                 *gy(my+4))*gz(mz+2)
     +                                 -(5.0d0*gx(mx+6)*gy(my)
     +                                 +15.0d0*gx(mx+4)*gy(my+2)
     +                                 +15.0d0*gx(mx+2)*gy(my+4)
     +                                 +5.0d0*gx(mx)*gy(my+6))*gz(mz))
c
                                 q61c = q61c + 0.125d0*f*rt(3)*rt(7)
     +                                  *(8.0d0*gz(mz+5)*gx(mx+1)*gy(my)
     +                                  -20.0d0*gz(mz+3)
     +                                  *(gx(mx+3)*gy(my)+gx(mx+1)
     +                                  *gy(my+2))+5.0d0*gz(mz+1)
     +                                  *(gx(mx+5)*gy(my)+2.0d0*gx(mx+3)
     +                                  *gy(my+2)+gx(mx+1)*gy(my+4)))
                                 q61s = q61s + 0.125d0*f*rt(3)*rt(7)
     +                                  *(8.0d0*gz(mz+5)*gx(mx)*gy(my+1)
     +                                  -20.0d0*gz(mz+3)
     +                                  *(gx(mx+2)*gy(my+1)+gx(mx)
     +                                  *gy(my+3))+5.0d0*gz(mz+1)
     +                                  *(gx(mx+4)*gy(my+1)
     +                                  +2.0d0*gx(mx+2)*gy(my+3)+gx(mx)
     +                                  *gy(my+5)))
c
                                 q62c = q62c + 0.03125d0*f*rt(14)*rt(15)
     +                                  *(16.0d0*gz(mz+4)
     +                                  *(gx(mx+2)*gy(my)-gx(mx)
     +                                  *gy(my+2))-16.0d0*gz(mz+2)
     +                                  *(gx(mx+4)*gy(my)-gx(mx)
     +                                  *gy(my+4))+gz(mz)
     +                                  *(gx(mx+6)*gy(my)+gx(mx+4)
     +                                  *gy(my+2)-gx(mx+2)*gy(my+4)
     +                                  -gx(mx)*gy(my+6)))
                                 q62s = q62s + 0.0625d0*f*rt(14)*rt(15)
     +                                  *(16.0d0*gz(mz+4)*gx(mx+1)
     +                                  *gy(my+1)-16.0d0*gz(mz+2)
     +                                  *(gx(mx+3)*gy(my+1)+gx(mx+1)
     +                                  *gy(my+3))+gz(mz)
     +                                  *(gx(mx+5)*gy(my+1)
     +                                  +2.0d0*gx(mx+3)*gy(my+3)
     +                                  +gx(mx+1)*gy(my+5)))
c
                                 q63c = q63c + 0.0625d0*f*rt(14)*rt(15)
     +                                  *(8.0d0*gz(mz+3)
     +                                  *(gx(mx+3)*gy(my)-3.0d0*gx(mx+1)
     +                                  *gy(my+2))-3.0d0*gz(mz+1)
     +                                  *(gx(mx+5)*gy(my)-2.0d0*gx(mx+3)
     +                                  *gy(my+2)-3.0d0*gx(mx+1)
     +                                  *gy(my+4)))
                                 q63s = q63s + 0.0625d0*f*rt(14)*rt(15)
     +                                  *(8.0d0*gz(mz+3)
     +                                  *(3.0d0*gx(mx+2)*gy(my+1)-gx(mx)
     +                                  *gy(my+3))-3.0d0*gz(mz+1)
     +                                  *(3.0d0*gx(mx+4)*gy(my+1)
     +                                  +2.0d0*gx(mx+2)*gy(my+3)-gx(mx)
     +                                  *gy(my+5)))
c
                                 q64c = q64c + 0.1875d0*f*rt(7)
     +                                  *(10.0d0*gz(mz+2)
     +                                  *(gx(mx+4)*gy(my)-6.0d0*gx(mx+2)
     +                                  *gy(my+2)+gx(mx)*gy(my+4))
     +                                  -gz(mz)
     +                                  *(gx(mx+6)*gy(my)+gx(mx)*gy
     +                                  (my+6))+5.0d0*gz(mz)
     +                                  *(gx(mx+4)*gy(my+2)+gx(mx+2)
     +                                  *gy(my+4)))
                                 q64s = q64s + 0.75d0*f*rt(7)
     +                                  *(10.0d0*gz(mz+2)
     +                                  *(gx(mx+3)*gy(my+1)-gx(mx+1)
     +                                  *gy(my+3))-gz(mz)
     +                                  *(gx(mx+5)*gy(my+1)-gx(mx+1)
     +                                  *gy(my+5)))
c
                                 q65c = q65c + 0.1875d0*rt(11)*rt(14)
     +                                  *f*gz(mz+1)
     +                                  *(gx(mx+5)*gy(my)-10.0d0*gx
     +                                  (mx+3)*gy(my+2)+5.0d0*gx(mx+1)
     +                                  *gy(my+4))
                                 q65s = q65s + 0.1875d0*rt(11)*rt(14)
     +                                  *f*gz(mz+1)
     +                                  *(5.0d0*gx(mx+4)*gy(my+1)
     +                                  -10.0d0*gx(mx+2)*gy(my+3)+gx(mx)
     +                                  *gy(my+5))
c
                                 q66c = q66c + 0.03125d0*rt(6)*rt(7)
     +                                  *rt(11)*f*gz(mz)
     +                                  *(gx(mx+6)*gy(my)
     +                                  -15.0d0*gx(mx+4)*gy(my+2)
     +                                  +15.0d0*gx(mx+2)*gy(my+4)-gx(mx)
     +                                  *gy(my+6))
                                 q66s = q66s + 0.0625d0*rt(6)*rt(7)
     +                                  *rt(11)*f*gz(mz)
     +                                  *(3.0d0*gx(mx+5)*gy(my+1)
     +                                  -10.0d0*gx(mx+3)*gy(my+3)
     +                                  +3.0d0*gx(mx+1)*gy(my+5))
c
c  q5 terms
 260                             q50 = q50 +
     +                                 0.125d0*f*(8.0d0*gz(mz+5)*gx(mx)
     +                                 *gy(my)-40.0d0*gz(mz+3)
     +                                 *(gx(mx+2)*gy(my)+gx(mx)*gy(my+2)
     +                                 )+15.0d0*gz(mz+1)
     +                                 *(gx(mx+4)*gy(my)+2.0d0*gx(mx+2)
     +                                 *gy(my+2)+gx(mx)*gy(my+4)))
c
                                 q51c = q51c + 0.125d0*rt(15)
     +                                  *f*(8.0d0*gz(mz+4)*gx(mx+1)
     +                                  *gy(my)-12.0d0*gz(mz+2)
     +                                  *(gx(mx+3)*gy(my)+gx(mx+1)
     +                                  *gy(my+2))+gz(mz)
     +                                  *(gx(mx+5)*gy(my)+2.0d0*gx(mx+3)
     +                                  *gy(my+2)+gx(mx+1)*gy(my+4)))
                                 q51s = q51s + 0.125d0*rt(15)
     +                                  *f*(8.0d0*gz(mz+4)*gx(mx)
     +                                  *gy(my+1)-12.0d0*gz(mz+2)
     +                                  *(gx(mx+2)*gy(my+1)+gx(mx)
     +                                  *gy(my+3))+gz(mz)
     +                                  *(gx(mx+4)*gy(my+1)
     +                                  +2.0d0*gx(mx+2)*gy(my+3)+gx(mx)
     +                                  *gy(my+5)))
c
                                 q52c = q52c + 0.25d0*rt(7)*rt(15)
     +                                  *f*(2.0d0*gz(mz+3)
     +                                  *(gx(mx+2)*gy(my)-gx(mx)
     +                                  *gy(my+2))-gz(mz+1)
     +                                  *(gx(mx+4)*gy(my)-gx(mx)
     +                                  *gy(my+4)))
                                 q52s = q52s + 0.5d0*rt(7)*rt(15)
     +                                  *f*(2.0d0*gz(mz+3)*gx(mx+1)
     +                                  *gy(my+1)-gz(mz+1)
     +                                  *(gx(mx+3)*gy(my+1)+gx(mx+1)
     +                                  *gy(my+3)))
c
                                 q53c = q53c + 0.0625d0*rt(7)*rt(10)
     +                                  *f*(8.0d0*gz(mz+2)
     +                                  *(gx(mx+3)*gy(my)-3.0d0*gx(mx+1)
     +                                  *gy(my+2))-gz(mz)
     +                                  *(gx(mx+5)*gy(my)-2.0d0*gx(mx+3)
     +                                  *gy(my+2)-3.0d0*gx(mx+1)
     +                                  *gy(my+4)))
                                 q53s = q53s + 0.0625d0*rt(7)*rt(10)
     +                                  *f*(8.0d0*gz(mz+2)
     +                                  *(3.0d0*gx(mx+2)*gy(my+1)-gx(mx)
     +                                  *gy(my+3))-gz(mz)
     +                                  *(3.0d0*gx(mx+4)*gy(my+1)
     +                                  +2.0d0*gx(mx+2)*gy(my+3)-gx(mx)
     +                                  *gy(my+5)))
c
                                 q54c = q54c + 0.375d0*rt(5)*rt(7)
     +                                  *f*gz(mz+1)
     +                                  *(gx(mx+4)*gy(my)-6.0d0*gx(mx+2)
     +                                  *gy(my+2)+gx(mx)*gy(my+4))
                                 q54s = q54s + 1.5d0*rt(5)*rt(7)
     +                                  *f*gz(mz+1)
     +                                  *(gx(mx+3)*gy(my+1)-gx(mx+1)
     +                                  *gy(my+3))
c
                                 q55c = q55c + 0.1875d0*rt(14)*f*gz(mz)
     +                                  *(gx(mx+5)*gy(my)
     +                                  -10.0d0*gx(mx+3)*gy(my+2)
     +                                  +5.0d0*gx(mx+1)*gy(my+4))
                                 q55s = q55s + 0.1875d0*rt(14)*f*gz(mz)
     +                                  *(5.0d0*gx(mx+4)*gy(my+1)
     +                                  -10.0d0*gx(mx+2)*gy(my+3)+gx(mx)
     +                                  *gy(my+5))
c
c  hexadecapole terms
 270                             q40 = q40 +
     +                                 0.125d0*f*(8.0d0*gx(mx)*gy(my)
     +                                 *gz(mz+4)
     +                                 -24.0d0*(gx(mx+2)*gy(my)+gx(mx)
     +                                 *gy(my+2))*gz(mz+2)
     +                                 +3.0d0*(gx(mx+4)*gy(my)+gx(mx)
     +                                 *gy(my+4))*gz(mz)
     +                                 +6.0d0*(gx(mx+2)*gy(my+2)*gz(mz))
     +                                 )
c
                                 q41c = q41c + rt(10)
     +                                  *f*(gx(mx+1)*gy(my)*gz(mz+3)
     +                                  -0.75d0*(gx(mx+3)*gy(my)
     +                                  +gx(mx+1)*gy(my+2))*gz(mz+1))
c
                                 q41s = q41s + rt(10)
     +                                  *f*(gx(mx)*gy(my+1)*gz(mz+3)
     +                                  -0.75d0*(gx(mx+2)*gy(my+1)
     +                                  +gx(mx)*gy(my+3))*gz(mz+1))
c
                                 q42c = q42c + 0.25d0*rt(5)
     +                                  *f*(6.0d0*(gx(mx+2)*gy(my)
     +                                  -gx(mx)*gy(my+2))*gz(mz+2)
     +                                  -(gx(mx+4)*gy(my)-gx(mx)
     +                                  *gy(my+4))*gz(mz))
c
                                 q42s = q42s + 0.5d0*rt(5)
     +                                  *f*(6.0d0*gx(mx+1)*gy(my+1)
     +                                  *gz(mz+2)
     +                                  -(gx(mx+3)*gy(my+1)+gx(mx+1)
     +                                  *gy(my+3))*gz(mz))
c
                                 q43c = q43c + 0.25d0*rt(10)*rt(7)
     +                                  *f*(gx(mx+3)*gy(my)
     +                                  -3.0d0*gx(mx+1)*gy(my+2))
     +                                  *gz(mz+1)
c
                                 q43s = q43s + 0.25d0*rt(10)*rt(7)
     +                                  *f*(3.0d0*gx(mx+2)*gy(my+1)
     +                                  -gx(mx)*gy(my+3))*gz(mz+1)
c
                                 q44c = q44c + 0.125d0*rt(5)*rt(7)
     +                                  *f*(gx(mx+4)*gy(my)
     +                                  -6.0d0*gx(mx+2)*gy(my+2)+gx(mx)
     +                                  *gy(my+4))*gz(mz)
c
                                 q44s = q44s + 0.5d0*rt(5)*rt(7)
     +                                  *f*(gx(mx+3)*gy(my+1)-gx(mx+1)
     +                                  *gy(my+3))*gz(mz)
c
c  octopole terms
c
 280                             q30 = q30 +
     +                                 f*(gx(mx)*gy(my)*gz(mz+3)-1.5d0*
     +                                 (gx(mx+2)*gy(my)+gx(mx)*gy(my+2))
     +                                 *gz(mz+1))
c
                                 q31c = q31c + rt(6)
     +                                  *f*(gx(mx+1)*gy(my)*gz(mz+2)
     +                                  -0.25d0*(gx(mx+3)*gy(my)
     +                                  +gx(mx+1)*gy(my+2))*gz(mz))
c
                                 q31s = q31s + rt(6)
     +                                  *f*(gx(mx)*gy(my+1)*gz(mz+2)
     +                                  -0.25d0*(gx(mx+2)*gy(my+1)
     +                                  +gx(mx)*gy(my+3))*gz(mz))
c
                                 q32c = q32c + 0.5d0*rt(15)
     +                                  *f*(gx(mx+2)*gy(my)-gx(mx)
     +                                  *gy(my+2))*gz(mz+1)
c
                                 q32s = q32s + rt(15)*f*gx(mx+1)
     +                                  *gy(my+1)*gz(mz+1)
c
                                 q33c = q33c + 0.25d0*rt(10)
     +                                  *f*(gx(mx+3)*gy(my)*gz(mz)
     +                                  -3.0d0*gx(mx+1)*gy(my+2)*gz(mz))
                                 q33s = q33s + 0.25d0*rt(10)
     +                                  *f*(3.0d0*gx(mx+2)*gy(my+1)
     +                                  *gz(mz)-gx(mx)*gy(my+3)*gz(mz))
c
c  quadrupole terms
c
 290                             q20 = q20 +
     +                                 f*(gx(mx)*gy(my)*gz(mz+2)-0.5d0*
     +                                 (gx(mx+2)*gy(my)+gx(mx)*gy(my+2))
     +                                 *gz(mz))
c
                                 q21c = q21c + f*rt(3)*gx(mx+1)*gy(my)
     +                                  *gz(mz+1)
c
                                 q21s = q21s + f*rt(3)*gx(mx)*gy(my+1)
     +                                  *gz(mz+1)
c
                                 q22c = q22c + f*0.5d0*rt(3)
     +                                  *(gx(mx+2)*gy(my)-gx(mx)
     +                                  *gy(my+2))*gz(mz)
c
                                 q22s = q22s + f*rt(3)*gx(mx+1)*gy(my+1)
     +                                  *gz(mz)
c
c  dipole terms
c
 300                             q10 = q10 + f*gx(mx)*gy(my)*gz(mz+1)
c
                                 q11c = q11c + f*gx(mx+1)*gy(my)*gz(mz)
c
                                 q11s = q11s + f*gx(mx)*gy(my+1)*gz(mz)
c
c  monopole term
c
 310                             q00 = q00 + f*gx(mx)*gy(my)*gz(mz)
c
c  end of loop over basis functions
 320                          continue
 330                       continue
c
                           k = nq**2
                           if (iw.gt.0) write (iw,6030) (qt(l),l=1,k)
c
c  move multipoles to expansion centre nearest to overlap centre p.
                           if (.not.ls) call dmamq(qc,xp,yp,zp,iw)
                           if (ls) call dmamqs(qc,xp,yp,zp,aa,iw)
                        end if
c
c  end of loop over primitives
 340                 continue
 350              continue
c
c  end of loop over shells
 360           continue
 370        continue
c
c  end of loop over atoms
 380     continue
 390  continue
c
      return
 6010 format ('     atoms   shells primitives            position'/
     +        '     multipole contributions ...'/)
 6020 format (1x,3(i5,i4),3x,3f10.5)
 6030 format (1x,f10.6/1x,3f10.6/1x,5f10.6/1x,7f10.6/1x,9f10.6/1x,
     +        11f10.6/1x,13f10.6)
      end
**==dmasdh.f
      subroutine dmasdh (x,y,z, j, r,max,iwr)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  computes regular solid harmonics r**k ckq(theta,phi) for ranks k up
c  to j, if j >= 0;
c  or irregular solid harmonics r**(-k-1) ckq(theta,phi) for ranks k up
c  to |j|, if j < 0.
c
c
_IF(f90test)
      use junk_dma, xx=>x
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
_IFN(f90test)
      logical lshb, linear
      common /junk/ oxyz(3), nc, lmax, tol,
     &           tshift, mindc, maxdc, lshb, linear, kwwww(8),
     &  xx(3,maxat), limit(maxat),radius(maxat), d(20,20),
     &           qx(121), binom(20,20), rt(20)
      character *8 zzzzz,name
      common/junkc/zzzzz(58),name(maxat)
_ENDIF
INCLUDE(common/mapper)
c
      dimension r(max)
c
c  locations in r are used as follows:
c
c        1    2    3    4    5    6    7    8    9   10   11  ...
c  kq = 00   10   11c  11s  20   21c  21s  22c  22s  30   31c ...
c
c  r(k,0) is real and is left in location k**2 + 1.
c  r(k,mc) and r(k,ms) are sqrt(2) times the real and imaginary parts
c  respectively of the complex solid harmonic r(k,-m)* = (-1)**m r(k,m),
c  and are left in locations k**2 + 2m and k**2 + 2m + 1 respectively.
c
c

      l = iabs(j)
      if ((l+1)**2.gt.max) then
c
       write (iwr,6010) l
       call caserr('insufficient memory for harmonics in dma analysis')
      else
         rr = x**2 + y**2 + z**2
         if (j.ge.0) then
c
c  regular
            r(1) = 1.0d0
            r(2) = z
            r(3) = x
            r(4) = y
            rfz = z
            rfx = x
            rfy = y
         else
c
c  irregular
            rr = 1.0d0/rr
            rfx = x*rr
            rfy = y*rr
            rfz = z*rr
            r(1) = dsqrt(rr)
            r(2) = rfz*r(1)
            r(3) = rfx*r(1)
            r(4) = rfy*r(1)
         end if
c
c  remaining values are found using recursion formulae, relating
c  the new set n to the current set k and the previous set p.
c
         k = 1
      end if
 20   n = k + 1
      ln = n*n + 1
      lk = k*k + 1
      lp = (k-1)**2 + 1
      a2kp1 = k + k + 1
c
c  obtain r(k+1,0) from r(k,0)*r(1,0) and r(k-1,0)
c
      r(ln) = (a2kp1*r(lk)*rfz-k*rr*r(lp))/(k+1)
c
      m = 1
      ln = ln + 1
      lk = lk + 1
      lp = lp + 1
      if (k.ne.1) then
c
c  obtain r(k+1,m) from r(k,m)*r(1,0) and r(k-1,m)
c
 30      r(ln) = (a2kp1*r(lk)*rfz-rt(k+m)*rt(k-m)*rr*r(lp))
     +           /(rt(n+m)*rt(n-m))
         r(ln+1) = (a2kp1*r(lk+1)*rfz-rt(k+m)*rt(k-m)*rr*r(lp+1))
     +             /(rt(n+m)*rt(n-m))
         m = m + 1
         ln = ln + 2
         lk = lk + 2
         lp = lp + 2
         if (m.lt.k) go to 30
      end if
c
c  obtain r(k+1,k) from r(k,k)*r(1,0)
c
      r(ln) = rt(n+k)*r(lk)*rfz
      r(ln+1) = rt(n+k)*r(lk+1)*rfz
      ln = ln + 2
c
c  obtain r(k+1,k+1) from r(k,k)*r(1,1)
c
      s = rt(n+k)/rt(n+n)
      r(ln) = s*(rfx*r(lk)-rfy*r(lk+1))
      r(ln+1) = s*(rfx*r(lk+1)+rfy*r(lk))
c
      k = k + 1
      if (k.lt.l) go to 20
      return
 6010 format (' insufficient array space for harmonics up to rank',i3)
      end
**==dmasq.f
      subroutine dmasq (q1,l1,m1, q2,m2, x,y,z, iw)
c
c-----------------------------------------------------------------dmasq
c
c
c  shift the multipoles q1 relative to the point (x,y,z) to the point
c  (0,0,0) and add them to the multipole expansion q2
c
c  multipoles of ranks l1 through m1 are to be tranferred from q1
c  and added to q2, keeping ranks up to m2 in the transferred expansion.
c
c  q2(l,m) = sum(k,q) sqrt | (l+m) (l-m) | * q1(k,q) * r(l-k,m-q)
c                          | (k+q) (k-q) |
c
c  where the quantities in the square root are binomial coefficients
c  and the r() are regular solid harmonics of (x,y,z).
c
_IF(f90test)
      use junk_dma, xx=>x
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
      COMPLEX      rc,qc,qz
      integer qmin, qmax, q, t1, t2
INCLUDE(common/sizes)
      dimension q1(121), q2(121), r(121)
      dimension rc(121), qc(121), qz(121)
c
_IFN(f90test)
      character *8 zzzzz,name
      common/junkc/zzzzz(58),name(maxat)
      logical lshb, linear
      common /junk/ oxyz(3), nc, lmax, tol,
     &           tshift, mindc, maxdc, lshb, linear, kwwww(8),
     &  xx(3,maxat), limit(maxat),radius(maxat), d(20,20),
     &           qx(121), binom(20,20), rt(20)
_ENDIF
INCLUDE(common/mapper)
c
c
c  real multipoles and harmonics are in the order
c  q00, q10, q11c, q11s, q20, ...
c
c  complex multipoles and harmonics are in the order
c  q00, q1-1, q10, q11, q2-2, ...
c
      data rthalf /0.7071067811865475244d0/, eps /0.0d0/
c
c
      if (l1.gt.m1 .or. l1.gt.m2) return
c
c  estimate largest significant transferred multipole.  the magnitude of
c  the contribution of the q(k,q) to q2(l,m) is of order
c  |q(k)|*r**(l-k), multiplied by binomial coefficient factors which we
c  estimate as 2**(l-k).  if the total of such estimates for rank l is
c  greater than sqrt(eps), transfers up to rank l are calculated
c  explicitly.
c

c  if eps is zero, this procedure is bypassed.
      n2 = m2
      k1 = max(1,l1)
      if (eps.ne.0.0d0 .and. m2.ne.0) then
c
         r2 = 4.0d0*(x**2+y**2+z**2)
         n2 = 0
         a = 0.0d0
         if (l1.eq.0) a = q1(1)**2
         do 30 k = k1 , m2
            a = a*r2
            if (k.le.m1) then
               t1 = k**2 + 1
               t2 = (k+1)**2
               do 20 i = t1 , t2
                  a = a + q1(i)**2
 20            continue
            end if
            if (a.gt.eps) n2 = k
 30      continue
      end if
c
c  evaluate solid harmonics in real form

      call dmasdh(x,y,z,n2,r,121,iw)
c
c  construct complex solid harmonics rc
      rc(1) = cmplx(r(1),0.0d0)
      do 50 k = 1 , n2
         kb = k**2 + k + 1
         km = k**2 + 1
         rc(kb) = cmplx(r(km),0.0d0)
         km = km + 1
         s = rthalf
         do 40 m = 1 , k
            s = -s
            rc(kb-m) = cmplx(rthalf*r(km),-rthalf*r(km+1))
            rc(kb+m) = cmplx(s*r(km),s*r(km+1))
            km = km + 2
 40      continue
 50   continue
c
c  construct complex multipoles qc corresponding to original
c  real multipoles q1
      if (l1.eq.0) qc(1) = cmplx(q1(1),0.0d0)
      if (m1.ne.0) then
         do 70 k = k1 , m1
            kb = k**2 + k + 1
            km = k**2 + 1
            qc(kb) = cmplx(q1(km),0.0d0)
            km = km + 1
            s = rthalf
            do 60 m = 1 , k
               s = -s
               qc(kb-m) = cmplx(rthalf*q1(km),-rthalf*q1(km+1))
               qc(kb+m) = cmplx(s*q1(km),s*q1(km+1))
               km = km + 2
 60         continue
 70      continue
      end if
c
c  construct shifted complex multipoles qz (only for non-negative m)
      if (l1.eq.0) qz(1) = qc(1)
      do 110 l = k1 , n2
         kmax = min(l,m1)
         lb = l**2 + l + 1
         lm = lb
         m = 0
c
 80      qz(lm) = 0.0d0
         if (l1.eq.0) qz(lm) = qc(1)*rc(lm)
         if (k1.le.kmax) then
            do 100 k = k1 , kmax
               qmin = max(-k,k-l+m)
               qmax = min(k,l-k+m)
               kb = k**2 + k + 1
               jb = (l-k)**2 + (l-k) + 1
c  special cases -- binom(n,0) is not tabulated
               if (qmin.eq.-k) qz(lm) = qz(lm) + binom(l-m,k-qmin)
     +                                  *qc(kb+qmin)*rc(jb+m-qmin)
               if (qmin.eq.-k) qmin = -k + 1
               if (qmax.eq.k) qz(lm) = qz(lm) + binom(l+m,k+qmax)
     +                                 *qc(kb+qmax)*rc(jb+m-qmax)
               if (qmax.eq.k) qmax = k - 1
               if (qmin.le.qmax) then
                  do 90 q = qmin , qmax
                     qz(lm) = qz(lm) + binom(l+m,k+q)*binom(l-m,k-q)
     +                        *qc(kb+q)*rc(jb+m-q)
 90               continue
               end if
 100        continue
         end if
c
         m = m + 1
         lm = lm + 1
         if (m.le.l) go to 80
 110  continue
c
c  construct real multipoles and add to q2
      if (l1.eq.0) q2(1) = q2(1) + dreal(qz(1))
      do 130 k = k1 , n2
         kb = k**2 + k + 1
         km = k**2 + 1
         q2(km) = q2(km) + dreal(qz(kb))
         s = 1.0d0/rthalf
         km = km + 1
         do 120 m = 1 , k
            s = -s
            q2(km) = q2(km) + s*dreal(qz(kb+m))
            q2(km+1) = q2(km+1) + s*dimag(qz(kb+m))
            km = km + 2
 120     continue
 130  continue
c
      return
      end
**==dmasz.f
      subroutine dmasz(q1,l1,m1, q2,m2, z)
c
c-----------------------------------------------------------------dmasz
c
c
c  shift those multipoles in q1 with l1 <= n <= m1 to the origin to
c  which the multipoles q2 are referred, and add the shifted values
c  to q2.
c  values are required in q2 only up to n=m2.
c
c  z gives the position of q1 relative to an origin at q2.
c
c  qn(0) = sum(s) \ (n!/s!(n-s)!) * qs(z) * z**(n-s) \
c
c  note that multipole qn is in array entry n+1.
c
c
_IF(f90test)
      use junk_dma
_ENDIF
      implicit REAL  (a-h,o-z),integer  (i-n)
      integer s, s1, s2, sp1
INCLUDE(common/sizes)
c
_IFN(f90test)
      character *8 zzzzz,name
      common/junkc/zzzzz(58),name(maxat)
      logical lshb, linear
      common /junk/ oxyz(3), nc, lmax, tol,
     &           tshift, mindc, maxdc, lshb, linear, kwwww(8),
     &  x(3,maxat), limit(maxat),radius(maxat), d(20,20),
     &           qx(121), binom(20,20), rt(20)
_ENDIF
INCLUDE(common/mapper)
c
      dimension q1(121), q2(121), zs(22)
c
      if (l1.gt.m1) return
      zs(1) = 1.0d0
      do 20 i = 1 , m2
         zs(i+1) = z*zs(i)
 20   continue
c
c  charge on q1
      if (l1.le.0) then
         q2(1) = q2(1) + q1(1)
         do 30 n = 1 , m2
            q2(n+1) = q2(n+1) + q1(1)*zs(n+1)
 30      continue
      end if
c
      if (m1.le.0) return
      s1 = max(1,l1)
      s2 = min(m1,m2)
      do 50 s = s1 , s2
         q2(s+1) = q2(s+1) + q1(s+1)
         if (s.eq.m2) return
         sp1 = s + 1
         do 40 n = sp1 , m2
            q2(n+1) = q2(n+1) + binom(n,s)*q1(s+1)*zs(n-s+1)
 40      continue
 50   continue
c
      return
      end
**==reden2.f
      subroutine reden2(trans,da,db,nconf,l1,l2)
c
c     ----- calculate the reduced density matrix -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/scfwfn)
      dimension trans(l1,*),da(l2),db(l2)
      dimension nconf(*)
      data dzero,two/0.0d0,-2.0d0/
c
c     ----- set up some control data needed to generate the data -----
c
      if (nco.ne.0) call setsto(nco,1,nconf)
      nopen = 0
      iconf = ncores + 1
      iorb = nco
      if (nseto.ne.0) then
         do 20 i = 1 , nseto
            nopen = nopen + no(i)
            nop = no(i)
            call setsto(nop,iconf,nconf(iorb+1))
            iorb = iorb + nop
            iconf = iconf + 1
 20      continue
      end if
      if (npair.ne.0) then
         npair2 = npair + npair
         do 30 i = 1 , npair2
            nconf(i+iorb) = iconf
            iconf = iconf + 1
 30      continue
      end if
      norb = nco + nopen + npair + npair
c
      call vclr(da,1,l2)
      call vclr(db,1,l2)
c
c
c     ----- get alpha part first -----
c
c     the open part of alpha will be off by 0.5
c
      iadd = 1
      do 50 i = 1 , l1
         do 40 k = 1 , norb
          dum = f(nconf(k))*trans(i,k)
          if (dum.ne.dzero) then
           call daxpy(i,dum,trans(1,k),1,da(iadd),1)
          end if
 40     continue
c
c     ----- now get beta part - first get only the open part -----
c
         iadd = iadd + i
 50   continue
      ilo = nco + 1
      ihi = nco + nopen
      if (nopen.gt.0) then
         iadd = 1
         do 70 i = 1 , l1
           do 60 k = ilo , ihi
            dum = f(nconf(k))*trans(i,k)
            if (dum.ne.dzero) then
             call daxpy(i,dum,trans(1,k),1,db(iadd),1)
            end if
 60        continue
         iadd = iadd + i
 70      continue
      end if
c
c     ----- now add this into alpha to get the correct alpha
c           density -----
c
      call vadd(da,1,db,1,da,1,l2)
c
c     ----- now subtract out two times what's in db to give the
c           correct beta -----
c
_IF1(civu)      call gtriad(l2,two,db,da,db)
_IFN1(civu)      call vsma(db,1,two,da,1,db,1,l2)
      return
      end
c
c potential fitted charges
c
      subroutine pdc1(core)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/potfit/dipx,dipy,dipz,scalex,charge,
     &     isec(mxcalc),nsec,
     &     ocut,odipol,ocharg,opschg,osym
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/runlab)
INCLUDE(common/coval)
c
c grid definition parameters
c
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
c
c data calculation parameters
c
      common/dfcalc/cdata(5,mxcalc),ictype(mxcalc),icsect(mxcalc),
     &            icgrid(mxcalc),icgsec(mxcalc),ndata(mxcalc),
     &            icstat(mxcalc),icdata(5,mxcalc),ncalc
      common/dfclc2/cscal(mxcalc),iblkp(mxcalc)
c
c largest number of equivalent atoms
c
      parameter(maxsym=8)
c
      character*10 charwall
      dimension core(*), dum(2)
      parameter (au2d = 2.541580609681d0)
      
      write(iwr,1000)
      if(ocharg)write(iwr,1010)charge
      if(odipol) then
         if(.not.ocharg) then
            charge = 0.0d0
            ocharg=.true.
            write(iwr,1020)
         endif
         write(iwr,1030)dipx,dipy,dipz
      else
         write(iwr,*)'no dipole constraint'
      endif
      if(.not.ocharg)write(iwr,*)'no total charge constraint'
      if(ocut)write(iwr,1040)scalex
      if(opschg)write(iwr,1045)
c
c   atom symmetry
c
      i0 = 1
      nav = lenwrd()
      isy = i0
      ieq = isy + maxsym*nat/nav
      ifree = ieq + (nat+1)/nav
c
      if(osym)then
         write(iwr,*)'use symmetry equivalences between atoms'
         call pdcsym(maxsym,core(isy),core(ieq),natsy,iwr)
      else
c set up unit mapping array for no symm case
         call pdcsy0(maxsym,core(isy),core(ieq),natsy)
      endif
c
c --- allocate memory for potential arrays
c
c  clear grid storage arrays 
c
      ngrid=0
      ncalc=0
c
c  running total of grid points
c
      ng=0
c
c --- read potentials from dumpfile
c
      ip=ifree
      do 10 i=1,nsec
         ier=0
         call rcasec(isec(i),icalc,idaf,ibl3d,core(ifree),.true.,ier)
         write(iwr,1050)ndata(icalc),isec(i)
         if(dabs(cscal(icalc)-1.0d0).gt.1.0d-20)then
            write(iwr,1055)cscal(icalc)
            nnn = ndata(icalc)
            ccc = 1.0d0/cscal(icalc)
            call dscal(nnn,ccc,core(ifree),1)
         endif
         ng = ng + ndata(icalc)
         ifree=ifree+ndata(icalc)
 10   continue
c
c --- read point coordinates from dumpfile
c
      ig=ifree
      ng1=0
      do 20 i=1,nsec
         ier=0
         call rgrsec(igsect(icgrid(i)),igrid,idaf,ibl3d,
     &        core(ifree),.true.,ier)
         ng1 = ng1+nptot(igrid)
         if(ngdata(igrid).ne.0)then
            write(iwr,1060)nptot(igrid),igsect(icgrid(i))
         else
c
c  compute coordinates
c            
            call setgrd(igrid,iwr)
            do 15 ii = 1,nptot(igrid)
               call getpt(ii,gx,gy,gz,dum,1,idum)
               core(ifree + 3*(ii-1) + 0) = gx
               core(ifree + 3*(ii-1) + 1) = gy
               core(ifree + 3*(ii-1) + 2) = gz
 15         continue
            write(iwr,1061)nptot(igrid),igsect(icgrid(i))
         endif
         ifree=ifree+nptot(igrid)*3
 20   continue
      if(ng.ne.ng1)call caserr('pdc1 - grid point counting error')
c
c --- allocate memory for remaining arrays
c
      isize=natsy
      if(odipol)then
         isize=isize+4
      else if(ocharg)then
         isize=isize+1
      endif
      ifit=ifree
      irad=ifit+isize
      ia = irad+nat
      ib = ia+(isize)*(isize)
      iaa =ib+isize
      iwk1=iaa+(isize)*(isize)
      iwk2=iwk1+isize
      irr=iwk2+isize
      ichs=irr + ng*nat
      ifree = ichs+natsy
c
c  set up squared covalent radii
c
      if(ocut)then
         owarn = .false.
         do 30 i = 1,nat
            izi = nint(czan(i))
            if(izi.gt.0)then
               core(irad + i - 1) = (scalex*cov(izi))**2
            else
               owarn = .true.
               core(irad + i - 1) = 0.0d0
            endif
 30      continue
         if(owarn)write(iwr,1063)
      endif
c
c  fit charges
c
      call pdc(nat,c,core(ip),core(ig),ng,core(ifit),
     &   core(ia),core(ib),core(iaa),core(iwk1),core(iwk2),
     &   isize,core(irad),core(irr),sdev,ctot,dipx1,dipy1,dipz1,
c    &   core(ifree),
     &   natsy,core(ieq),maxsym,core(isy),core(ichs),iwr)
c
      write(iwr,1065)sdev
      dipxyz = sqrt (dipx1*dipx1 + dipy1*dipy1 + dipz1*dipz1)
      write(iwr, 1070) ctot
      write(iwr, 1080) dipx1,dipy1,dipz1,dipxyz
      dipx1 = dipx1 * au2d
      dipy1 = dipy1 * au2d
      dipz1 = dipz1 * au2d
      dipxyz = dipxyz * au2d
      write(iwr,1090) dipx1,dipy1,dipz1,dipxyz
c
c write out the atom labels and their respective charges
c
      write(iwr,1100)
      do 40 i = 1,nat
         write(iwr,1110)zaname(i),(c(j,i),j=1,3),core(ifit+i-1)
 40   continue
c
c punchfile output
c
      call blkpdc(core(ifit))
c
      write(iwr,1120)cpulft(1),charwall()
      return
 1000 format(1x,104('=')//46x,'potential derived charges module'//
     &     1x,104('=')//10x,'Control Parameters'/10x,18('-'),/)
 1010 format(1x,'Total charge will be constrained to ',f10.5)
 1020 format(1x,'Due to the dipole constraint the total charge will',
     &     ' be constrained to 0.0')
 1030 format(1x,'The dipole moment will be constrained to ',
     &     '(',f10.5,',',f10.5,',',f10.5,' ) atomic units')
 1040 format(1x,'Points close to the nuclei will be excluded,',
     &     ' scale factor for the covalent radii = ',f8.3)
 1045 format(1x,'Make correction for pseudopotential cores')
 1050 format(1x,i5,' data values read from dumpfile section ',i3)
 1055 format(1x,'scale factor ',f20.6,' removed')
 1060 format(1x,i5,' sets of grid coordinates read from ',
     &     'dumpfile section ',i3)
 1061 format(1x,i5,' sets of grid coordinates generated from ',
     &     'definition on dumpfile section ',i3)
 1063 format(//1x,'WARNING - points close to dummy centres and point ',
     &     'charges will not be removed')
 1065 format(//,1x,'Fitting Results',/1x,15('-'),/,
     &  ' The Standard Deviation is ', f12.8)
 1070 format(1x,'The Total Charge is     ', f12.6)
 1080 format(/1x,'DIPOLE MOMENTS (atomic units)'/
     -        '         x ',f14.7/
     -        '         y ',f14.7/
     -        '         z ',f14.7/
     -        '     total ',f14.7)
 1090 format(/,1x,'DIPOLE MOMENTS (debye)'/
     -        '         x ',f14.7/
     -        '         y ',f14.7/
     -        '         z ',f14.7/
     -        '     total ',f14.7)
 1100 format(//21x,'Potential Derived Charges',/,21x,25('='),//,
     &     7x,'symbol',13x,'x           y           z',9x,' charge')
 1110 format(10x,a8,3f12.5,3x,f10.4)
 1120 format(//1x,104('-')//1x,'end of potential derived charges at ',
     &     f8.2,' seconds',a10,' wall')
      end
      subroutine pdcsym(maxsym,isym,neq,natsy,iwr)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
      dimension isym(maxsym,*), neq(*)
c
c GAMESS symmetry
c
INCLUDE(common/symtry)
INCLUDE(common/transf)
      dimension tempi2(3),cnew(3)
      equivalence (cnew(1),pnew)
      data toler / 1.0d-6 /
c
      call ptgrp(0)
      natsy = 0
      do 50 i = 1,nat
         call lframe(c(1,i),c(2,i),c(3,i),
     &        tempi2(1),tempi2(2),tempi2(3))
c        oeq = .false.
         do 45 j = 1, natsy
            call lframe(c(1,isym(1,j)),c(2,isym(1,j)),c(3,isym(1,j)),
     &           psmal, qsmal, rsmal)
            do 43 it = 1 , nt
               nn = 9*(it - 1)
               call trans1(nn)
               tester = dist(tempi2,cnew)
               if(tester.lt.toler)then
                  neq(j) = neq(j) + 1
                  isym(neq(j),j) = i
c                 oeq = .true.
                  goto 47
               endif
 43         continue
 45      continue
         natsy=natsy+1
         neq(natsy)=1
         isym(1,natsy)=i
         do 46 jjj=2,maxsym
            isym(jjj,natsy)=0
 46      continue
 47      continue
 50   continue
      write(iwr,*)
      write(iwr,*)' unique atom   symmetry related atoms'
      do 51 iii = 1,natsy
         write(iwr,1031)(isym(jjj,iii),jjj=1,neq(iii))
 1031    format(1x,i4,10x,7i4)
 51   continue
      return
      end
      subroutine pdcsy0(maxsym,isym,neq,natsy)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension isym(maxsym,*), neq(*)
INCLUDE(common/infoa)

      do 53 iii = 1,nat
         neq(iii)=1
         isym(1,iii)=iii
         do 52 jjj=2,maxsym
            isym(jjj,iii)=0
 52      continue
 53   continue
      natsy=nat
      end
      subroutine pdc(nat, coord, gmep, grid, ng, chg, 
     &     a,b,aa,wks1,wks2,
     &     ia,radat2,r,sdev,ctot,dipx1,dipy1,dipz1,
c    &     core,
     &     natsy,neq,nsmax,isym,chgsy,iw)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      common/potfit/dipx,dipy,dipz,scalex,charge,
     &     isec(mxcalc),nsec,
     &     ocut,odipol,ocharg,opschg,osym
      dimension coord(3,*),gmep(*),grid(3,*),chg(*)
      dimension a(ia,*),b(*),aa(ia,*),wks1(*),wks2(*)
      dimension radat2(*), r(ng,*)
c     dimension core(*)
      dimension chgsy(*), neq(*),isym(nsmax,natsy)
c
c If scaling has been asked for call rempts
c to remove lines close to atoms
c
c     if(ocut)call rempts(coord,radat2,grid,gmep,nat,ng,scalex,iw)
      if(ocut)call rempts(coord,radat2,grid,gmep,nat,ng,iw)
c
c  check point count
c
      if (ng.le.ia)
     &   call caserr('too few potential points to perform fit')
c
c If pseudo potential charges were read in modify the 
c electrostatic potential
c
      if(opschg)call pschgs(coord,grid,gmep,pschg,nat,ng)
c 
c calculate the distance from each grid point to the current atom
c
      do 10 loop = 1, nat
         do 20 i = 1, ng
            r(i, loop) = dist(coord(1,loop),grid(1,i))
 20      continue
 10   continue
c     
c For each point find the closest atom
c keep a sum of the number of such points for each atom
c and a running total of the distance squared of such a point
c
c      call info(coord, r, mdec, ng, grid, label, nat,
c     -            cstart,cend,clist,contour,ict)
c
c Construct the linear equation Ab=x
c
      n = natsy
      do 218 i = 1, n
c Construct A
         do 201 j = 1,n
            a(j,i) = 0.0d0
 201     continue
         do 208 ieq = 1,neq(i)
            do 207 j = 1, i
               do 206 jeq = 1, neq(j)
            do 205 ip = 1 , ng
                     a(j,i) = a(j,i) +
     &                    1.0d0/(r(ip,isym(ieq,i))*r(ip,isym(jeq,j)))
 205        continue
 206           continue
 207        continue
 208     continue
c Construct b
         do 211 j = 1,n
            a(i,j) = a(j,i)
 211     continue
         b(i) = 0.0d0
         do 214 ieq=1,neq(i)
         do 212 ip = 1 , ng
               b(i) = b(i) + gmep(ip)/r(ip,isym(ieq,i))
 212     continue
 214     continue
 218  continue
c
      if(ocharg) then
         do 220 i = 1, natsy
            a(i,n+1) = 0.0d0
            a(n+1,i) = 0.0d0
            do 221 ieq = 1,neq(i)
               a(i,n+1) = a(i,n+1) + 1.0d0
               a(n+1,i) = a(n+1,i) + 1.0d0
 221        continue
 220     continue
         a(n+1,n+1) = 0.0d0
         b(n+1) = charge
         n = n + 1
      endif
c
      if(odipol) then
         do 230 i = 1,natsy
            a(i,n+1)=0.0d0
            a(i,n+2)=0.0d0
            a(i,n+3)=0.0d0
            a(n+1,i)=0.0d0
            a(n+2,i)=0.0d0
            a(n+3,i)=0.0d0
            do 229 ieq=1,neq(i)
               a(i,n+1)=a(i,n+1) + coord(1,isym(ieq,i))
               a(i,n+2)=a(i,n+2) + coord(2,isym(ieq,i))
               a(i,n+3)=a(i,n+3) + coord(3,isym(ieq,i))
               a(n+1,i)=a(n+1,i) + coord(1,isym(ieq,i))
               a(n+2,i)=a(n+2,i) + coord(2,isym(ieq,i))
               a(n+3,i)=a(n+3,i) + coord(3,isym(ieq,i))
 229        continue 
 230     continue
        a(n  ,n+1)=0.0d0
        a(n  ,n+2)=0.0d0
        a(n  ,n+3)=0.0d0
        a(n+1,n  )=0.0d0
        a(n+1,n+1)=0.0d0
        a(n+1,n+2)=0.0d0
        a(n+1,n+3)=0.0d0
        a(n+2,n  )=0.0d0
        a(n+2,n+1)=0.0d0
        a(n+2,n+2)=0.0d0
        a(n+2,n+3)=0.0d0
        a(n+3,n  )=0.0d0
        a(n+3,n+1)=0.0d0
        a(n+3,n+2)=0.0d0
        a(n+3,n+3)=0.0d0
        b(n+1) = dipx
        b(n+2) = dipy
        b(n+3) = dipz
        n = n + 3
      endif
_IF(usenag)
c
c Call NAG routine to solve linear equations
c
      ifail = 0
      tol = 5.0d-8
      call f04atf(a,ia,b,n,chg,aa,ia,wks1,wks2,ifail)
c      call f04jdf(n,n,a,ia,b,tol,sigma,irank,aa,ia*ia,ifail)
c
c must check ifail value before proceeding
c
      if (ifail .ne. 0) 
     &     call caserr('nag failure in f04jdf, called from pdc')
_ELSE
c
c invert a by diagonalisation
c aa holds eigenvectors, wks1 evectors
c
      ifail=0
      call eigen(a,ia,ia,wks1,aa,wks2,ifail)
      if(ifail.ne.0)call caserr('eigen failure in pdc1')
      do 240 i = 1,ia
         if(dabs(wks1(i)).gt.1.0d-10)then
            wks1(i) = 1.0d0/wks1(i)
         else
            write(iw,*)'small evalue'
         endif
 240  continue
c
c transform the evectors into aa
c
      do 120 i = 1 , n
         do 100 j = i , n
            rx = 0.0d0
            do 90 k = 1 , n
               rx = rx + aa(i,k)*aa(j,k)*wks1(k)
 90         continue
            wks2(j) = rx
 100     continue
         do 110 j = i , n
            aa(i,j) = wks2(j)
 110     continue
 120  continue
c
      do 140 i = 1 , n
         do 130 j = 1 , i
            aa(i,j) = aa(j,i)
 130     continue
 140  continue
c
c  check inversion
c
      sum = 0.0d0
c     oerr=.false.
      do 251 ii = 1,ia
         sum = 0.0d0
         do 250 jj=1,ia
            sum = sum + a(ii,jj)*aa(jj,ii)
 250     continue
         if(dabs(sum - 1.0d0).gt.1.0d-8)then
            write(iw,*)'non unit diagonal',ii,sum
c           oerr=.true.
         endif
 251  continue
c      if(oerr)call caserr('inversion failure in pdc')
c
c  copy back to a
c
      call dcopy(ia*ia,aa,1,a,1)
c
c compute charges
c
      do 253 ii = 1,ia
         chgsy(ii) = 0.0d0
         do 252 jj = 1,ia
            chgsy(ii) = chgsy(ii) + a(jj,ii)*b(jj)
 252     continue
 253  continue
      ifail=0
_ENDIF
c
c sort out symmetry equivalent charges      
c
      do 255 ii = 1,natsy
         do 254 ieq = 1,neq(ii)
            chg(isym(ieq,ii)) = chgsy(ii)
 254     continue
 255  continue
c
c calcuate fitting error
c
      fsumq = 0.0d0
      do 270 k = 1,ng
         gnew = 0.0d0
         do 265 i = 1,nat
            gnew = gnew + chg(i)/r(k,i)
 265     continue
         fsumq = fsumq + (gmep(k) - gnew)**2
 270  continue
      sdev = sqrt(fsumq/(ng-1))

      ctot = 0.0d0
      dipx1 = 0.0d0
      dipy1 = 0.0d0
      dipz1 = 0.0d0
      do 260 i = 1, nat
         ctot = ctot + chg(i)
         dipx1 = dipx1 + coord(1,i)*chg(i)
         dipy1 = dipy1 + coord(2,i)*chg(i)
         dipz1 = dipz1 + coord(3,i)*chg(i)
 260  continue
c
      return
      end
c
      subroutine rempts (coord,radat2,grid,gmep,non,ng,
c    +                   scalex,iw)
     +                   iw)
c
c routine which takes a set of points and removes points
c which are within a certain distance from atoms.  
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension coord(3,*), radat2(*), grid(3,*), gmep(*)
      integer cnsidr
c
c statement function for calculating square distances.
c
      discal(ax1,ay1,az1,ax2,ay2,az2) = (ax1-ax2)*(ax1-ax2) + 
     +                            (ay1-ay2)*(ay1-ay2) +
     +                            (az1-az2)*(az1-az2)
c
c consider all the points from the current point which have not yet
c been considered.
c
      next = 0
      ndel = 0
      do 10 cnsidr = 1,ng
         ax1 = grid(1,cnsidr)
         ay1 = grid(2,cnsidr)
         az1 = grid(3,cnsidr)
c
c Find any atom closer to the point than radat
c
         do 20 j = 1, non
            dist = discal(ax1, ay1, az1, 
     &           coord(1,j), coord(2,j), coord(3,j))
            if (dist .lt. radat2(j)) then 
               ndel = ndel + 1
               goto 10
            endif
 20      continue
c             
c move the point in cnsidr to the available place held in 'next + 1'
c
         next = next + 1 
         grid(1,next) = grid(1,cnsidr)
         grid(2,next) = grid(2,cnsidr)
         grid(3,next) = grid(3,cnsidr)
         gmep(next)  = gmep(cnsidr)
 10   continue
c
c update the number of points now held in coordinate arrays
c
c     ic = next
      if(ndel .gt. 0) then
         write(iw,*)
         write(iw,*)'Points too close to atoms deleted ',ndel
      endif
      ng = ng - ndel
      return
      end
c
      subroutine pschgs(coord,grid,gmep,pschg,non,ic)
c
c routine to remove the pseudo potential charges from the grid
c 
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)

      dimension coord(3,*), grid(3,*), gmep(*), pschg(*)
c
c Loop arond the grid points
c
      do 900 j = 1, ic
          do 60 i = 1, non
c
c Calculate the pdc potential at this point due to atom i
c
             rr = dist(coord(1,i),grid(1,i))
c
c Remove that component of the potential
c
             gmep(j) = gmep(j) - pschg(i) / rr
  60      continue
  900 continue
      return
      end
c
c routines to initialise potential fitting code
c
c 0 - inquiry
c 1 - initialise
c
      logical function opotf(iflag)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/potft1/opotf1
      common/potfit/dipx,dipy,dipz,scalex,charge,
     &     isec(mxcalc),nsec,
     &     ocut,odipol,ocharg,opschg,osym
c     data opotft1 /.false./
      if(iflag.eq.1)then
         opotf1=.true.
         odipol=.false.
         ocharg=.false.
         opschg=.false.
         ocut=.false.
         osym=.false.
         nsec=0
      else if(iflag.eq.2)then
         opotf1=.false.
      endif
      opotf=opotf1
      return
      end
      subroutine potfsc(isec1)
c
c register request for a dumpfile section
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/potfit/dipx,dipy,dipz,scalex,charge,
     &     isec(mxcalc),nsec,
     &     ocut,odipol,ocharg,opschg,osym
      if(nsec.eq.mxcalc)
     &     call caserr('too many sections specified for potfit')
      nsec = nsec+1
      isec(nsec)=isec1
      return
      end
      subroutine potfdp(dx,dy,dz)
c
c set dipole for fitting
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/potfit/dipx,dipy,dipz,scalex,charge,
     &     isec(mxcalc),nsec,
     &     ocut,odipol,ocharg,opschg,osym
      odipol=.true.
      dipx=dx
      dipy=dy
      dipz=dz
      return
      end
      subroutine potfsy
c
c switch on symmetry for atomic charges fitting
c
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/potfit/dipx,dipy,dipz,scalex,charge,
     &     isec(mxcalc),nsec,
     &     ocut,odipol,ocharg,opschg,osym
      osym=.true.
      return
      end
      subroutine potfch(charg1)
c
c set total charge restraint
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/potfit/dipx,dipy,dipz,scalex,charge,
     &     isec(mxcalc),nsec,
     &     ocut,odipol,ocharg,opschg,osym
      ocharg=.true.
      charge=charg1
      return
      end
      subroutine potfrc(scalx1)
c
c set atomic distance cutoff
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/potfit/dipx,dipy,dipz,scalex,charge,
     &     isec(mxcalc),nsec,
     &     ocut,odipol,ocharg,opschg,osym
      ocut=.true.
      scalex=scalx1
      return
      end
      subroutine ver_analc(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/analc.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
