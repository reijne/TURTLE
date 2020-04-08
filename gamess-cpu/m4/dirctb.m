c 
c  $Author: mrdj $
c  $Date: 2010-09-17 16:11:09 +0200 (Fri, 17 Sep 2010) $
c  $Locker:  $
c  $Revision: 6195 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/dirctb.m,v $
c  $State: Exp $
c  
      function b22(hh12,hh22)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
INCLUDE(common/dir_erg)
      b=h11-hh22
      c= dsqrt(b*b+4.0d0*hh12*hh12)
      if(moddia)111,222,222
c-------   l o c k    m o d e   ------------------------
222   if(b)111,22,22
22    e=(h11+hh22+c)*0.5d0
      goto 112
c-------   m i n i m i z e    e n e r g y   ---------------
111    e=(h11+hh22-c)*0.5d0
112    b=e-hh22
       absa= dabs(hh12)
       absb= dabs(b)
       if(absb.ge.absa)goto  1
c... use first line of secular equation
        if(absa.ge.(1.0d-10))goto 2
        if(hh12)4,3,3
4       b22=(h11-e)*1.0d+10
        return
3       b22=(e-h11)*1.0d+10
        return
2       b22=(e-h11)/hh12
        return
c... use second line of secular equation
1       if(absb.lt.(1.0d-10))goto 998
        b22=hh12/b
        return
998     b22=0.0d0
        return
        end
      subroutine diabas(iconf,iwr)
c...
c... to generate main and disk store relative addresses for c-z
c... vectors in iconf(1,...) and (iconf(2,...) respectiveley
c... and to compute starting disk block addresses for c-z vectors
c... in index(26-29) and start disc block address for diagonal of
c... hamiltonian in index(199)
c...
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf,pack2
_ELSEIF(convex,i8drct)
      integer *8 iconf,pack2
_ELSEIF(bits8)
      REAL iconf
_ELSE
      REAL iconf,pack2
_ENDIF
      character *4 pin,out,diagmo
      dimension iconf(5,*),diagmo(3)
INCLUDE(common/sizes)
INCLUDE(common/helpr)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88,khb111,khb222,khb333
     *,khbiii,khbkkk
INCLUDE(common/cdcryz)
INCLUDE(common/cdcryi)
INCLUDE(common/ccntl)
INCLUDE(common/disktl)
INCLUDE(common/symci)
INCLUDE(common/auxh)
INCLUDE(common/cntrl)
      common/symcon/root2,root2i
INCLUDE(common/timanb)
INCLUDE(common/outctl)
INCLUDE(common/diactl)
INCLUDE(common/corctl)
INCLUDE(common/ccepa)
INCLUDE(common/dir_erg)
INCLUDE(common/mr_mp2)
INCLUDE(common/mr_cepctl)
INCLUDE(common/modeci)
c
_IF(i8)
      pack2 (i1,i2)  = ior(ishft(i1,32),i2)
_ENDIF
c
      data pin,out/'in','out'/
      data diagmo/'emin','lock','vmin'/
      kk = 0
      ll = 0
      kkk = 0
      lll = 0
      kbl = 0
c
      timbeg=cpulft(1)
      call walltime(wtimbeg)
      if (maxcyc.gt.0) 
     *   write(iwr,1000)timbeg,diagmo(moddia+2),maxcyc,lcall,thresh,
     * eshift
      if (maxcyc.gt.0) call prjajd('direc',iconf,iwr)
1000  format(//1x,80('=')//
     *' **** diagonalizer called at',f9.2,' secs'//
     *' mode = ',a4/
     *' maximum number of cycles =',i4/
     *' maximum davidson subspace = ',i4/
     *' threshold = ',e15.6/
     *' level shift parameter = ',e15.6/)
2000  format(1x,23('-')/1x,a4,'store mode selected'/1x,23('-'))
      if(malt) write(iwr,2323)
 2323 format(/' level shift alternation requested')
      if(moddia.eq.1)lvar=.true.
       root2= dsqrt(2.0d0)
       root2i=1.0d0/root2
c...      vacuum states
      if(nval)2,1,2
2     do 22 loop=1,nnst
_IF(cray,t3e,convex,i8,i8drct)
      iconf(1,loop)=kk
      kk=kk+iconf(3,loop)
_ELSEIF(bits8)
      call pad8(kk,iconf(1,loop))
      kk=kk+ipad(iconf(3,loop))
_ELSE
      iconf(1,loop)=pad(kk)
      kk=kk+ipad(iconf(3,loop))
_ENDIF
22    continue
c...     doublet states
 1    kbl=(nval+510)/511
_IF(cray,t3e,convex,i8,i8drct)
       iconf(2,1)=kbl
_ELSEIF(bits8)
       call pad8(kbl,iconf(2,1))
_ELSE
       iconf(2,1)=pad(kbl)
_ENDIF
        if(ndoub)4,3,4
c       note starting disk address for doublet c-z vector in iconf(2,1)
4      kbl=(ndoub+510)/511+kbl
      do 44 loop=nstrt1,nlast1
_IF(cray,t3e,convex,i8,i8drct)
      iconf(1,loop)=ll
      ll=norbe(iconf(4,loop))*iconf(3,loop)+ll
_ELSEIF(bits8)
      call pad8(ll,iconf(1,loop))
      ll=norbe(ipad(iconf(4,loop)))*ipad(iconf(3,loop))+ll
_ELSE
      iconf(1,loop)=pad(ll)
      ll=norbe(ipad(iconf(4,loop)))*ipad(iconf(3,loop))+ll
_ENDIF
44    continue
c...       n-2 states
3     if(nmin2)6,5,6
6      do 66 loop=nstrt2,nlast2
       call upack2(iconf(3,loop),nn,mm)
_IF(cray,convex,t3e,i8,i8drct)
       isi=iconf(4,loop)
_ELSE
       isi=ipad(iconf(4,loop))
_ENDIF
       mm=mm*norbsi(isi)
       nn=nn*norbtr(isi)
_IF(bits8)
       call pack2(lll,kkk,iconf(1,loop))
_ELSEIF(i8drct)
       iconf(1,loop) = pack2(lll,kkk)
_ELSE
       iconf(1,loop) = pack2(lll,kkk)
_ENDIF
       kkk=kkk+mm
       lll=lll+nn
       lbl=(mm+510)/511+kbl
_IF(bits8)
       call pack2(lbl,kbl,iconf(2,loop))
_ELSEIF(i8drct)
       iconf(2,loop) = pack2(lbl,kbl)
_ELSE
       iconf(2,loop)=pack2(lbl,kbl)
_ENDIF
66     kbl=(nn+510)/511+lbl
c... now set up disk starting block numbers for
c... c1 , z1 , c2 , z2 , diagonal of hamiltonian
c
c
5      index(199)=nxblk
ccepa  contracted vectors space on disc  (batches of 3) c.c/c.z/z.c
       nxblk = nxblk + kbl
       if (icepa.ge.0) then
          ncpcc = nmin2*2 + nmin1 + nnst - nref + 1
          ncpbl = (ncpcc*3-1)/511 + 1
          nxblk = nxblk + kbl
          kcpbl = nxblk
          nxblk = nxblk + (lcall*(lcall+1)/2)*ncpbl
       else
          ncpcc = 0
       end if
ccepa
_IF(ga)
c
c rwah create global array ed19 file (numcz) with search

       call create_ed19(kbl,iwr)
       call create_index(iconf,iconf)
       call ibasgn(2,1,kbl,index(26))
       call ibasgn(2,nxblk+kbl,kbl,index(28))
       ispbig=kbl
       ispsma=(ntotal+510)/511
       call ibasgn(118,index(29)+kbl,ispsma,index(30))
_ELSE
       call ibasgn(4,nxblk,kbl,index(26))
       nxblk = index(29)+kbl
       ispbig=kbl
       ispsma=(ntotal+510)/511
       call ibasgn(118,nxblk,ispsma,index(30))
_ENDIF
c... set up base points for case where c and z in store
      khbz=khbill
      khbc=khbz+ntotal
      if (icepa.ge.0) then
         icpsto=.true.
         kcepc=khbc
         kcepz=khbz
         kcpcc=khbc+ntotal
         kdavcz=kcpcc+5*ncpcc
         kcpcd=kdavcz+ntotal
         khbb=kcpcd+mxdvpa
      else  
         khbb=khbc+ntotal
      endif
       kumpe=khbb
       jack14=khbb+lbig
      khb1=khbb+m2exch
      khb44=khb1+mxsspa
      khb55=khb44+mxttpa
      khb66=khb55+mxstpa
      khb88=khb66+mxstpa
      khb2=khbb+m2fock
      khb3=khb2+m2coul
      khb4=khb3+mxsspa
      khb5=khb4+mxttpa
      khb7=khb5+mxstpa
      khb8=khb7+nubsq
      lack7=nubsq+mubsq+mubsq
      jack7=khb8+lack7
      jack9=khb88+lack7
       khb111=khbb+m2em1e
       khb222=khb111+mx12pa
       khb333=khb222+mubsq
      lack5=nvmax*max12
      jack5=khb333+lack5
       khbiii=khbb+m2psup
      khbkkk=khbiii+mx12pa
      kkk=max(maxpa(1),maxpa(2))
      lack12=max(nvmax*kkk,nubtr*maxpa(3))
      jack12=khbkkk+lack12
      khbzd=khbz+nval
      khbzs=khbzd+ndoub
      khbzt=khbzs+nsingl
      khbcd=khbc+nval
      khbcs=khbcd+ndoub
      khbct=khbcs+nsingl
      jack1=mxvvpa+khbb
      jack2=mxddpa+khbb
      lack4=m2em1e+mxdvpa
      jack4=lack4+khbb
      mack=nvmax*maxpa(3)+mxddpa
      lack6=max(mack,mubsi)+m2exab+mubsq
      jack6=lack6+khbb
      lack8=m2exch+mack
      jack8=khbb+lack8
      jack10=khb1+mxv2pa
      jack11=khb2+mxv2pa
      jack13=khbb+n127sq
      mxsstt=mxsspa+mxttpa
      jack3=khbb+mxsstt
      if(nval)31,30,31
30    jack4=0
      jack10=0
      jack11=0
31    if(ndoub)33,32,33
32    jack2=0
      jack4=0
      jack5=0
      jack6=0
      jack8=0
      jack12=0
33    if(nmin2)35,34,35
34    jack3=0
      jack5=0
      jack7=0
      jack9=0
      jack10=0
      jack11=0
      jack12=0
      jack13=0
35    continue
c...
c...   contracted  vectors core-partitioning / instore
c...
      mmmmm = max(jack1,jack2,jack3,jack4,jack5,jack6,jack7,jack8,
     *jack9,jack10,jack11,jack12,jack13,jack14)
      if(mmmmm .gt.ngot.and. in_store) then
       write(iwr,2003) ngot, mmmmm
2003  format(/1x,' ***** insufficient memory for in-store mode'//
     +        1x,' available = ', i10, ' words'/
     +        1x,' required  = ', i10, ' words'//)
       call caserr('insufficient memory for in-store mode')
      endif
_IF(ga)
c
c rwah parallel ga forces out store mode ....
c
      out_store=.true.
      in_store=.false.
c
_ENDIF
      if(mmmmm .gt.ngot.or. out_store) goto 999
      write(iwr,2000)pin
      return
c... out store params
999   if (out_store) then
       write(iwr,2002)
2002   format(/1x,'**** out store mode requested **** ')
      endif
      instor=.false.
ccepa  / c205n
c...  contracted vector core-partitioning for out-store mode
c...  icpsto and davidson may still be possible
c
      kcepz = khbill
      kcepc = kcepz + ntotal
c...   5 ... diag/diag**2/c.c/c.z/z.c
      kcpcc = kcepc + ntotal
c     nuse = kcpcc + ntotal
      khbb = kcpcc + 5*ncpcc
      nuse=khbb
      if (odiag_1) go to 976
      if (lcall.eq.2) goto 978
      if (nuse.le.ngot) go to 1001
c...  sorry not enough core / no icpsto or davidson
c...  if n*n davidson is possible is further decided in davids
976   icpsto = .false.
      kcepz=khbill
      kcepc=khbill+kbig
      kcpcc=kcepz+ntotal
      khbb=kcpcc + 5*ncpcc
      if (2*kbig.gt.ntotal) goto 978
      nuse=khbb
      if (nuse.le.ngot) goto 1001
978   if (icepa.ge.0) then
         lcall = 2
         kcepz=khbill
         kcepc=kcepz+kbig
         kcpcc=kcepc+kbig
         nuse = kcpcc + 5*ncpcc
         icpsto = .false.
c...     no ispsto in direcg (dircta)
      endif
c     nuse = 0
c
ccepa   / c205n
1001  kumpc=khbill+kbig
      kumpe=kumpc+kbig
      jack14=kumpe+lbig
      khgot2=khbill+(ngot-khbill)/2
      lump1=khbzd
      lump2=lump1+nval
      jack1=lump2+mxvvpa
      lump3=khbill+ndoub
      lump4=lump3+ndoub
      jack2=lump4+mxddpa
      lump5=khbill+lack4
      lump6=lump5+nval
      lump7=lump6+ndoub
      if (icepa.ge.0) lump7=lump6+max(mxdvpa,ndoub)
      lump8=lump7+nval
      jack4=lump8+ndoub
      lump9=khbill+lack6
      lump10=lump9+ndoub
      jack6=lump10+ndoub
      lump11=lack8+khbill
      lump12=lump11+ndoub
      jack8=lump12+ndoub
      lump13=khbill+jbig
      lump14=khbill+m2exch
      lump99=lump14-nlast1
      lump15=lump14+jbig
      lump16=lump15+jbig
      lump17=lump16+nval
      lump18=lump17+nval
      jack10=lump18+mxv2pa
      lump19=khbill+m2fock
      lump20=lump19+mbig
      lump21=lump20+mbig
      lump22=lump21+nval
      lump23=lump22+nval
      jack11=lump23+mxv2pa
      lump24=khbill+n127sq
      lump25=lump24+(ngot-lump24)/2
      jack13=lump24+nbig+nbig
      lump26=khbill-nlast1
      lump27=khbill+nmin2
      lump28=ngot-lump27
      lump29=lump27+mxsstt
      mgot0=(ngot-lump29)/4
      lump30=lump29+mgot0
      lump31=lump30+mgot0
      lump32=lump31+mgot0
      lump33=lump19+m2coul
      khb3=lump33+nmin2
      lump33=lump33-nlast1
      khb4=khb3+mxsspa
      khb5=khb4+mxttpa
      khb7=khb5+mxstpa
      khb8=khb7+nubsq
      lump38=khb8+lack7
      mgot1=(ngot-lump38)/4
      lump39=lump38+mgot1
      lump40=lump39+mgot1
      lump41=lump40+mgot1
      khb1=lump14+nmin2
      khb44=khb1+mxsspa
      khb55=khb44+mxttpa
      khb66=khb55+mxstpa
      khb88=khb66+mxstpa
      lump47=khb88+lack7
      mgot2=(ngot-lump47)/4
      lump48=lump47+mgot2
      lump49=lump48+mgot2
      lump50=lump49+mgot2
      lump51=khbill+m2em1e
      lump52=lump51+ndoub
      lump53=lump52+ndoub
      lump54=lump53+jbig
      khb111=lump54+jbig
      khb222=khb111+mx12pa
      khb333=khb222+mubsq
      jack5=khb333+lack5
_IF(ga)
      lump55=khbill
_ELSE
      lump55=khbill+m2psup
_ENDIF
      lump56=lump55+jbig
      lump57=lump56+jbig
      lump58=lump57+ndoub
      khbiii=lump58+ndoub
      khbkkk=khbiii+mx12pa
_IF(ga)
c     write(6,*)'khbkkk,lack12 = ', khbkkk,lack12
      lump59=khbkkk+lack12
c     write(6,*)'ngott,lump59 =', ngott,lump59
      jack12=max(ngot,lump59)
_ELSE
      jack12=khbkkk+lack12
_ENDIF
      if(nval)131,130,131
130   jack4=0
      jack10=0
      jack11=0
131   if(ndoub)133,132,133
132   jack2=0
      jack4=0
      jack5=0
      jack6=0
      jack8=0
      jack12=0
133   if(nmin2)135,134,135
134    jack5=0
       jack10=0
       jack11=0
       jack12=0
      jack13=0
135   nuse=max(jack1,jack2,jack4,jack5,jack6,jack8,jack10,jack11,
     *jack12,jack13,jack14)
      write(iwr,2000)out
      call corfai('genz')
      return
      end
      subroutine genc1(cr,iconf,iwr)
c... to generate trial vector from input data or default info
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer   iconf
_ELSEIF(convex,i8drct)
      integer *8   iconf
_ELSE
      REAL  iconf
_ENDIF
      dimension cr(*),iconf(5,*)
INCLUDE(common/trial)
INCLUDE(common/symci)
INCLUDE(common/ccntl)
INCLUDE(common/corctl)
INCLUDE(common/auxh)
      common/stoctl/khbill,khbz,khbc
INCLUDE(common/disktl)
INCLUDE(common/cntrl)
INCLUDE(common/diactl)
INCLUDE(common/mr_mp2)
INCLUDE(common/ccepa)
      common/remc_ed19/nini,ipoint
      nlow = 0
      if(instor)goto 888
c--------- out of store version ----------------------------------------
      if(mrstor)222,111,222
_IF(ga)
111   call search(index(26),numcz)
_ELSE
111   call search(index(26),numci)
_ENDIF
      call vclr(cr(khbill+1),1,lbig)
c
      if (iselvc.ne.0.or.krefmp.gt.0) then 
         if (krefmp.gt.0) then
c...     start vector was already determined for mrmp2 ==> use that
            call dcopy(nrfcsf,cr(krefmp+1),1,cr(khbill+1),1)
_IF(ga)
            call wrt3(cr(khbill+1),nval,index(26),numcz)
_ELSE
            call wrt3(cr(khbill+1),nval,index(26),numci)
_ENDIF
            call vclr(cr(khbill+1),1,nrfcsf)
         else
c...     determine start-vector by diag. h-matrix
            call setc1d(cr(khbill+1),iconf,enul,ngot-khbill,iwr)
_IF(ga)
            call wrt3(cr(khbill+1),nval,index(26),numcz)
_ELSE
            call wrt3(cr(khbill+1),nval,index(26),numci)
_ENDIF
            call vclr(cr(khbill+1),1,lbig)
         end if
c
         nlow = nval
         modev = .false.
         nvec = 0
         call setc1(cr,moded,nlow,ndoub)
      else
         call setc1(cr,modev,nlow,nval)
         call setc1(cr,moded,nlow,ndoub)
      endif
      if(nmin2)66,999,66
66    nlox=nvds
      do 100 loop=nstrt2,nlast2
_IF(cray,convex,t3e,i8,i8drct)
        isi=iconf(4,loop)
_ELSE
        isi=ipad(iconf(4,loop))
_ENDIF
         call upack2(iconf(3,loop),lamt,lams)
         call setc1(cr,modest,nlow,lams*norbsi(isi))
100      call setc1(cr,modest,nlox,lamt*norbtr(isi))
      goto 999
c---------- in store version -------------------------------------------
888   if(mrstor)444,333,444
444   call getcz(iconf,cr(khbc+1),index(26))
222   modev=nnst.eq.0
      moded=nmin1.eq.0
      modest=nmin2.eq.0
      goto 999
333   continue
      call vclr(cr(khbc+1),1,ntotal)
c
      if (iselvc.ne.0.or.krefmp.gt.0) then 
         if (krefmp.gt.0) then
c...     start vector was already determined for mrmp2 ==> use that
            call dcopy(nrfcsf,cr(krefmp+1),1,cr(khbc+1),1)
         else
c
c...   start-vector from selected part of vac. h
c...   core partitioning : see out store 
c...   note khbz is before khbc so shift c vector afterwards
c
            call setc1d(cr(khbill+1),iconf,enul,ngot-khbill,iwr)
            call dcopy(nval,cr(khbill+1),1,cr(khbc+1),1)
c...   clear the rest of c (and z to be sure)
            call vclr(cr(khbc+nval+1),1,ntotal-nval)
            call vclr(cr(khbz+1),1,ntotal)
       end if
c
        modev = .false.
        go to 998
c
      endif
c
      do 1410 loop=1,nvec
      m=ivec(loop)
      if(m.lt.1.or.m.gt.ntotal)call caserr(
     *'invalid ICSF parameter in TRIAL directive')
      cr(khbc+m)=tvec(loop)
      if(m.gt.nval)goto 1411
      modev=.false.
      goto 1410
1411  if(m.gt.nvd)goto 1412
      moded=.false.
      goto 1410
1412  modest=.false.
1410  continue
998    call dumpcz(iconf,cr(khbc+1),index(26))
999   return
      end
      subroutine setc1(cr,mode,nlow,nn)
      implicit REAL  (a-h,o-z),integer  (i-n)
      logical mode
INCLUDE(common/sizes)
      dimension cr(*)
      common/moco/jvec(64)
      common/stoctl/khbill
INCLUDE(common/trial)
INCLUDE(common/comrjb)
INCLUDE(common/disktl)
      common/remc_ed19/nini,ipoint
c
      if(nn)888,999,888
888   nhi=nlow+nn
      kh=khbill-nlow
c...  read guess from numci (only reference configs)
      if (nlow.ne.0.or.ncv.le.0) go to 555
      n = iposun(numci)
      call rdedx(cr(kh+1),ncv,iblkcv,numci)
      call search(n,numci)
      mode = .false.
      go to 444
c
555   n=0
      do 1 loop=1,nvec
      m=ivec(loop)
      if(m.le.nlow.or.m.gt.nhi)goto 1
      mode=.false.
      m=kh+m
      cr(m)=tvec(loop)
      n=n+1
      jvec(n)=m
1     continue
_IF(ga)
444   continue
      njnj=iposun(numcz)
      call wrt3s(cr(khbill+1),nn,numcz)
_ELSE
444   call wrt3s(cr(khbill+1),nn,numci)
_ENDIF
      if (nlow.ne.0.or.ncv.le.0) go to 333
      call vclr(cr(kh+1),1,ncv)
      n = 0
333   if(n)666,777,666
666   do 2 loop=1,n
2     cr(jvec(loop))=0.0d0
777   nlow=nhi
999   return
      end
      subroutine setc1d(v,iconf,enul,ngot,iwr)
      implicit REAL (a-h,o-z)
c
c...  subroutine to determine ci-start vector by diagonalizing a
c...  selected valence-space h-matrix
c...  needs 2*nval space (possibly)
c...  this routine has been simplified (back), 
c...  while a davidson has been added
c...  for larger cases ; If you want changes please contact 
c...  Joop van Lenthe (joop@chem.uu.nl)
c
      dimension v(*),iconf(5,*)
c...   v is at least nval*2 + maxpat**2
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
c
INCLUDE(common/sizes)
INCLUDE(common/trial)
INCLUDE(common/disktl)
INCLUDE(common/cntrl)
INCLUDE(common/ccntl)
INCLUDE(common/helpr)
INCLUDE(common/auxh)
INCLUDE(common/dir_erg)
      parameter (mselvc=100)
      common/junk/ igenz(4088),eig(mselvc),f(mselvc*(mselvc+1)/2),
     1             q(mselvc,mselvc),ilifq(mselvc)
c
c...  iselvc = -1 => all nref
c...  iselvc = -2 => all nval
c...   else nselvc
c...  iselvc = 1 => first  nselvc
c...  iselvc = 2 => lowest neslvc
c     
      if (ngot-nval*2-maxpat**2.lt.1) then
         write(iwr,605) nval*2+maxpat**2-ngot
         call caserr('Not enough memory for setc1d')
      endif
      n_in_set = nselvc
      n_in_set = min(n_in_set,nval)
      if (iselvc.eq.-1) then
c
c...  count # csf's 
c
         nrfcsf = 0
         do 10 i=1,nref
            call upack2(iconf(3,i),jj,ii)
            nrfcsf = nrfcsf + ii
10       continue
         nselvc = max(nrfcsf,1)
         n_in_set = max(nrfcsf,1)
      else if (iselvc.eq.-2) then
         nselvc = nval
         n_in_set = nval
      end if
c
      nselvc = min(nselvc,mselvc,nval)
c
c...  iselvc : 1  =>  select first nselv states
c...           2  =>  select lowest nselv states
c
      if (iselvc.eq.1.or.nselvc.eq.n_in_set) then
         call ibasgn(nselvc,1,1,v(nval+1))
      else
c...  read diagonal elements
         call rdedx(v,nval,index(199),numci)
c...  sort diagonal in ascending order / take perm (in isel) along
         call sora(v,v(nval+1),nval,mselvc+1)
c...  be sure to take all degenerate states but no more then <mselvc>
         nadd = 1
210      if (dabs(v(nselvc+1)-v(nselvc)).le.1.0d-5) then
            if (nselvc.eq.mselvc) nadd = -1
            nselvc = nselvc + nadd
            if (nselvc.lt.nval) go to 210
         end if
      end if
c
c...  we have selected the states to be included / get h-matrix
c
      call rhval(f,v(nval+1),v(nval*2+1),nselvc,iconf)
_IF(parallel)
      nnlk=nselvc*(nselvc+1)/2
      call pg_dgop(2233,f,nnlk,'+')
_ENDIF
c
c...  diagonalise h-matrix
c
      if (iprvc.gt.10) then
         write(iwr,601)
601      format(/'   //trial//  selected h-matrix ')
         call prtri(f,nselvc)
      end if
c
      call ibasgn(nselvc,0,mselvc,ilifq)
      call jacobi(f,iky,nselvc,q,ilifq,nselvc,eig,2,2,1.0d-12)
c
c...  jacobi call parameter one before last determines eigenvalue order
c...  0 nop / 1 unsorted / 2 increasing / 3 decreasing / 4 lock
c
      if (iprvc.gt.10) then
         write(iwr,602)
602      format(/'    //trial//  eigen-vectors  ')
         call prvc3(q,mselvc,nselvc,eig,1,iwr)
      end if
c
c...  now get the required state
c
      call vclr(v,1,nval)
      call dsctr(nselvc,q(1,istvc),v(nval+1),v)  
c
      enul = eig(istvc)

      if (nselvc.ne.n_in_set) then
        nselvc = n_in_set
        kz = nval+1
        kb = kz+nval
        kcz = kb + maxpa(4)*maxpa(4)
        maxdav = (ngot-kcz)/(2*nselvc)
        if (maxdav.lt.5) then
           write(iwr,605) 5*2*nselvc-ngot+kcz
           call caserr('Not enough memory for vacuum davidson')
        endif
        call davidv(v,v(kz),v(kb),nselvc,enul,100,thresh,
     1              eshift,v(kcz),maxdav,iconf,iwr,.false.)
      endif
c
      if (iprvc.gt.5)
     * write(iwr,600) nselvc,istvc,eig(istvc)
600   format(/' start-vector from ',i2,'-space h-matrix',/
     *        '   state ',i2,'    energy ',d22.14,/)
605   format(' Enlarge memory by ',i10,' words')
c
      return
      end
      subroutine prvc3 (c,n,m,eps,ich,iwr)
      implicit REAL (a-h,o-z)
      character *5 sub
c
c   routine for printing of column-vectors and eigenvalues.
c
c            c  : matrix of coefficients
c            n  : dimension of array c(n,m)
c            m  : x columns to be printed # rows
c            eps: vector of eigenvalues
c            ich: 1 => occ or eig / 0 => column-numbers
c
c
      dimension c(n,m), eps(m)
      dimension sub(2)
      data sub /'-----','-----'/
      ncc = 10
      nbl = (m-1)/ncc
      nlast = m-ncc*nbl
      nbl = nbl+1
      nc1 = 0
      do 60 ib=1,nbl
      if ( ib .eq. nbl ) ncc=nlast
      nc0 = nc1+1
      nc1 = nc1+ncc
      if (ich.ne.0) write(iwr,10) ( eps(ic), ic=nc0,nc1 )
      if (ich.eq.0) write(iwr,11) ( ic, ic=nc0,nc1 )
   10 format(////8x,10f12.6)
   11 format(////8x,10(4x,i4,4x))
      write(iwr,20) ( sub, i=1,ncc )
   20 format(8x,10(2x,2a5))
      write(iwr,30)
   30 format(1h )
      do 50 ia=1,n
      write(iwr,40) ia, ( c(ia,ic), ic=nc0,nc1 )
   40 format(3x,i3,2x,10f12.6)
   50 continue
   60 continue
      write(iwr,30)
c
      return
      end
      subroutine sora(v,ip,nv,ns)
c
c...  sort ns elements of the nv elemensts of v in ascending order
c...  the permutation is returned in ip
c
      implicit REAL (a-h,o-z)
      dimension v(nv),ip(nv)
c
      nn = min(ns,nv-1)
      call ibasgn(nv,1,1,ip)
c
      do 10 j=1,nn
c...     (205 minimu starts at 0 / others at 1 )
      i=idmin(nv+1-j,v(j),1) + j-1
         temp = v(j)
         itemp = ip(j)
         v(j) = v(i)
         ip(j) = ip(i)
         v(i) = temp
         ip(i) = itemp
10    continue
c
      return
      end
      subroutine sori(iv,ip,nv,ns)
c
c...  sort ns elements of the nv elemensts of iv in ascending order
c...  the permutation is returned in ip
c
      implicit REAL (a-h,o-z)
      dimension iv(nv),ip(nv)
c
      nn = min(ns,nv-1)
      call ibasgn(nv,1,1,ip)
c
      do 10 j=1,nn
         i = j
         do 5 k=j,nv
5        if (iv(k).lt.iv(i)) i = k
         iemp = iv(j)
         itemp = ip(j)
         iv(j) = iv(i)
         ip(j) = ip(i)
         iv(i) = iemp
         ip(i) = itemp
10    continue
c
      return
      end
_EXTRACT(rhval,hp800)
      subroutine rhval(h,ip,cr,nsel,iconf)
c
c...  read the valence h-matrix for nsel selected (ip) csf's
c
      implicit REAL (a-h,o-z)
      dimension h(*),cr(*),ip(nsel),iconf(5,*)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
c
INCLUDE(common/disktl)
c       (direct)     common /spew/ model,ic,ir
      common/spew/rmode(3)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
c        (note  ic >= ir)
c
c...  sort ip in ascending order
c
      call sori(ip,h,nsel,nsel)
c
      call vclr(h,1,nsel*(nsel+1)/2)
      call brtsb(index(1))
c
999   call getsb(rmode(1),1)
      if (model.ne.2) return
      call getsb(rmode(2),2)
_IF(cray,convex,t3e,i8,i8drct)
      nc=iconf(3,ic)
      nr=iconf(3,ir)
_ELSE
      nc=ipad(iconf(3,ic))
      nr=ipad(iconf(3,ir))
_ENDIF
      call getsb(cr(1),nc*nr)
_IF(cray,convex,t3e,i8,i8drct)
      icbas=iconf(1,ic)
      irbas=iconf(1,ir)
_ELSE
      icbas=ipad(iconf(1,ic))
      irbas=ipad(iconf(1,ir))
_ENDIF
c     ... since ic only goes up ...
      if (ip(nsel).lt.icbas) return
c
      do 100 icp = 1,nsel
         if (ip(icp).gt.icbas+nc) go to 999
         if (ip(icp).le.icbas) go to 100
c...  loop over de ir
         do 50 irp = 1,icp
            if (ip(irp).gt.irbas+nr) go to 100
            if (ip(irp).le.irbas) go to 50
c...  found an element
            ih = icp*(icp-1)/2 + irp
            icr = (ip(icp)-icbas-1)*nr + ip(irp) - irbas
            h(ih) = cr(icr)
c
50       continue
100   continue
c
      go to 999
c
      end
_ENDEXTRACT
      subroutine vvijkl(cr,iconf)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
INCLUDE(common/prnprn)
      dimension cr(*),iconf(5,*)
      common/stoctl/khbill,khbz,khbc,khbb
      common/spew/rmode(3)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
c...
c... (ij/kl)    vacuum/vacuum    interactions
c...
999   call getsb(rmode(2),2)
_IF(cray,convex,t3e,i8,i8drct)
      nc=iconf(3,ic)
      nr=iconf(3,ir)
_ELSE
      nc=ipad(iconf(3,ic))
      nr=ipad(iconf(3,ir))
_ENDIF
      call getsb(cr(khbb+1),nc*nr)
_IF(cray,convex,t3e,i8,i8drct)
      icbas=iconf(1,ic)
      irbas=iconf(1,ir)
_ELSE
      icbas=ipad(iconf(1,ic))
      irbas=ipad(iconf(1,ir))
_ENDIF
      call mxmb(cr(khbb+1),1,nr,cr(icbas+khbc+1),1,1,cr(irbas+khbz+1),1,
     *1,nr,nc,1)
      if(ic.ne.ir)call mxmb(cr(khbb+1),nr,1,cr(irbas+khbc+1),1,1,
     *cr(icbas+khbz+1),1,1,nc,nr,1)
      if(odebug(40)) then
       call dbgciv('vvijkl',cr(irbas+khbz+1),nr,1,1,1)
       if(ic.ne.ir)call dbgciv('vvijkl',cr(icbas+khbz+1),nc,1,1,1)
      endif
      call getsb(rmode(1),1)
      if(model.eq.2)goto 999
      return
      end
      subroutine dumpcz(iconf,cz,ibl)
c... to dump c or z vector to ed7 when in store mode ok
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension cz(*),iconf(5,*)
INCLUDE(common/cntrl)
INCLUDE(common/ccntl)
INCLUDE(common/symci)
INCLUDE(common/disktl)
      call search(ibl,numci)
      call wrt3s(cz,nval,numci)
      call wrt3s(cz(nval+1),ndoub,numci)
      if(nmin2)1,999,1
1      kk=nvd
      ll=nvds
      do 100 loop=nstrt2,nlast2
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,loop)
_ELSE
      isi=ipad(iconf(4,loop))
_ENDIF
      call upack2(iconf(3,loop),lamt,lams)
      lams=lams*norbsi(isi)
      lamt=lamt*norbtr(isi)
      call wrt3s(cz(kk+1),lams,numci)
      kk=kk+lams
      call wrt3s(cz(ll+1),lamt,numci)
      ll=ll+lamt
100   continue
999   return
      end
      subroutine getcz(iconf,cz,ibl)
c... to restore c or z vecor when instor mode ok
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension cz(*),iconf(5,*)
INCLUDE(common/cntrl)
INCLUDE(common/ccntl)
INCLUDE(common/disktl)
INCLUDE(common/symci)
      call search(ibl,numci)
      call reads(cz,nval,numci)
      call reads(cz(nval+1),ndoub,numci)
      if(nmin2)1,999,1
1      kk=nvd
       ll =nvds
       do 100 loop=nstrt2,nlast2
_IF(cray,convex,t3e,i8,i8drct)
       isi=iconf(4,loop)
_ELSE
       isi=ipad(iconf(4,loop))
_ENDIF
       call upack2(iconf(3,loop),lamt,lams)
       lams=lams*norbsi(isi)
       lamt=lamt*norbtr(isi)
       call reads(cz(kk+1),lams,numci)
       kk=kk+lams
       call reads(cz(ll+1),lamt,numci)
       ll=ll+lamt
100    continue
999   return
      end
      subroutine putcz(iconf,cz,ibl)
      implicit REAL  (a-h,o-z), integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension cz(*),iconf(5,*)
INCLUDE(common/cntrl)
INCLUDE(common/ccntl)
INCLUDE(common/symci)
INCLUDE(common/disktl)
      call search(ibl,numci)
      call wrt3s(cz,nval,numci)
      call wrt3s(cz(nval+1),ndoub,numci)
      if(nmin2)1,999,1
1     kk=nvd
      do 100 loop=nstrt2,nlast2
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,loop)
_ELSE
      isi=ipad(iconf(4,loop))
_ENDIF
      call upack2(iconf(3,loop),lamt,lams)
      lams=norbsi(isi)*lams
      lamt=norbtr(isi)*lamt
      call wrt3s(cz(kk+1),lams,numci)
      kk=kk+lams
      call wrt3s(cz(kk+1),lamt,numci)
100   kk=kk+lamt
999   return
      end
      subroutine eblock(jconf,iconf,iclo,mgot)
      implicit REAL  (a-h,o-z),integer  (i-n)
c
c...   in case o errors in outstore mode please contact
c...      Joop van Lenthe
c...   common lsort communicates with termsb
c
      integer ibang
_IF(cray,t3e,i8)
      integer iconf,jconf
_ELSEIF(convex,i8drct)
      integer *8 iconf,jconf
_ELSE
      REAL iconf,jconf
_ENDIF
      dimension jank(7)
      dimension jconf(*),iconf(5,*)
      dimension iclo(*)
INCLUDE(common/diactl)
       common/disksa/bof(511),no,iblo
INCLUDE(common/symci)
INCLUDE(common/ccntl)
INCLUDE(common/outctl)
      common/stoctl/khbill
      common/spew/rmode(5)
INCLUDE(common/disktl)
      common/lsort/ibang(255),nwr(256),icnow,ko,kblo,iw1
      common/junk /iclow(255),ichi(255),ibias(255),ibiat(255),mblock
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir),
     * (rmode(4),kkkbas),(rmode(5),itrans)
      data jank/-1,0,0,0,0,0,1/
c...
c...  establish blocking of n-2  c  and  z  vectors
c...
      mblock=0
      mwrdst=0
      ibs=0
      ibt=0
      do 1 ic=nstrt2,nlast2
      call upack2(iconf(3,ic),ntt,nss)
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,ic)
_ELSE
      isi=ipad(iconf(4,ic))
_ENDIF
      nss=norbsi(isi)*nss
      ntt=norbtr(isi)*ntt
2     mr=mwrdst+nss+ntt
      if(mr.le.mgot)goto 3
      if(mwrdst)984,989,984
984   ibiat(mblock)=mwrdss
      ichi(mblock)=ic-1
      mwrdst=0
      goto 2
3     if(mwrdst)5,4,5
4     mblock=mblock+1
      if(mblock.gt.255)goto 989
      iclow(mblock)=ic
      ibias(mblock)=ibs
      mwrdss=ibt
_IF(cray,t3e,convex,i8,i8drct)
    5 jconf(ic+lump26)=mblock
_ELSEIF(bits8)
    5 call pad8(mblock,jconf(ic+lump26))
_ELSE
    5 jconf(ic+lump26)=pad(mblock)
_ENDIF
      mwrdst=mr
      ibt=ibt-ntt
      ibs=ibs-nss
1     mwrdss=mwrdss+nss
      ibiat(mblock)=mwrdss
      ichi(mblock)=nlast2
      if(mblock.le.2)goto 12345
c...
c...      sort  the  coupling  coefficients
c...
      mblo=iblo
      mo=no
      nodel=model
      lodel=jank(nodel-3)
      iwrd=3+lodel
      iw1=iwrd+1
      icnow=0
99    lo=no
      lblo=iblo
      call getsb(rmode(2),iwrd)
      call upack2(iconf(3,ic),lamtc,lamsc)
      call upack2(iconf(3,ir),lamtr,lamsr)
      if(lodel)20,21,22
c...     (ij/kl)
20    nsx=lamsc*lamsr+lamtc*lamtr
      goto 23
c...     (ij/ab) + (ia/ib) + fock(ab)
21    if(kkkbas)84,83,84
83    nsx=0
      goto 23
84    if(ic.ne.ir)goto 20
      nsx=lamsc*(lamsc+lamtc)+lamtc*lamtc
      goto 23
c...     (ia/jb)
22    nsx=(lamsc+lamtc)*(lamsr+lamtr)
23    continue
_IF(cray,convex,t3e,i8,i8drct)
      mc=jconf(ic+lump26)
      mr=mc-jconf(ir+lump26)
_ELSE
      mc=ipad(jconf(ic+lump26))
      mr=mc-ipad(jconf(ir+lump26))
_ENDIF
      if(mc.eq.icnow)goto 29
      if(icnow)28,27,28
c...
c...   clear up buckets generated from previous set of columns
c...
   28    call termsb(jconf)
c...
c...   initiate new set of columns
c...
27    icnow=mc
      ko=lo
      kblo=lblo
      mx=lump28/mc
      call setsto(mc,0,nwr)
      call ibasgn(mc,lump27,mx,ibang)
29    nww=nwr(mr+1)
      nw=nww+nsx+iw1
      if(nw.gt.mx)go to 989
      ib=nww+ibang(mr+1)
      nwr(mr+1)=nw
      call dcopy(iw1,rmode,1,jconf(ib+1),1)
      call getsb(jconf(ib+iw1+1),nsx)
      call getsb(rmode(1),1)
      if(model.eq.nodel)goto 99
c...
c... clear up buckets generated from final set of columns
c...
      iw1=0
      call termsb(jconf)
      call possb(mo,mblo)
12345 msort=.true.
_IF(cray,convex)
      call fmove(iclow,iclo,1021)
_ELSE
      call icopy(1021,iclow,1,iclo,1)
_ENDIF
      return
989   call caserr('insufficient memory available in eblock')
      return
      end
      subroutine setpag(jconf,icl,mumper)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer jconf
_ELSEIF(convex,i8drct)
      integer *8 jconf
_ELSE
      REAL jconf
_ENDIF
c...
c... to initiate pageing process for n-2 c and z vectors
c...
      dimension icl(1021),mumper(*)
      dimension jconf(*)
      common/junk /icll(255),ichh(255),ibiass(255),ibiatt(255),mblok,
     *mump(4)
      common/pager/iparit,khbjco,nwolf,iclogo,icval,izval,
     *mxx,khbzsx,khbztx,khbcsx,khbctx,mumpzx,mumpcx,
     *myy,khbzsy,khbzty,khbcsy,khbcty,mumpzy,mumpcy
      iclogo=0
      mxx=0
      myy=0
_IF(cray,convex)
      call fmove(mumper,mump,4)
      call fmove(icl,icll,1021)
_ELSE
      call icopy(4,mumper,1,mump,1)
      call icopy(1021,icl,1,icll,1)
_ENDIF
      do 1 m=1,mblok
       il=icll(m)
_IF(cray,t3e,i8)
      call setsto(ichh(m)-il+1,m,jconf(khbjco+il))
_ELSE
      call setstz(ichh(m)-il+1,m,jconf(khbjco+il))
_ENDIF
 1    continue
      return
      end
      subroutine pagexy(cr,iconf,jconf)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf,jconf
_ELSEIF(convex,i8drct)
      integer *8  iconf,jconf
_ELSE
      REAL iconf,jconf
_ENDIF
c...
c... control routine for n-2 c and z vector pageing
c...
      dimension jconf(*),cr(*),iconf(5,*)
      common/spew/rmode(5)
      common/junk /icll(255),ichh(255),ibiass(255),ibiatt(255),mblok,
     *mump(4)
      common/pager/iparit,khbjco,nwolf,iclogo,icval,izval,
     *mxx,khbzsx,khbztx,khbcsx,khbctx,mumpzx,mumpcx,
     *myy,khbzsy,khbzty,khbcsy,khbcty,mumpzy,mumpcy
INCLUDE(common/symci)
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt,
     *khcsc,khzsc,khctc,khztc,isc,khcsr,khzsr,
     *khctr,khztr,isr,norsir,nortrr
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir),
     * (rmode(4),kkkbas),(rmode(5),itrans)
      call getsb(rmode(2),nwolf)
      if(ic.eq.iclogo)goto 333
_IF(cray,convex,t3e,i8,i8drct)
      mcc=jconf(ic+khbjco)
_ELSE
      mcc=ipad(jconf(ic+khbjco))
_ENDIF
      iclogo=ic
      if(mcc.eq.mxx)goto 222
      call putxy(cr,iconf,khbzsy,khbzty,myy)
      myy=mxx
      khbzsy=khbzsx
      khbzty=khbztx
      khbcsy=khbcsx
      khbcty=khbctx
      mumpzx=mump(iparit+2)
      mumpcx=mump(iparit+3)
      iparit=-iparit
      mumpzy=mump(iparit+2)
      mumpcy=mump(iparit+3)
      mxx=mcc
      ibs=ibiass(mcc)
      ibt=ibiatt(mcc)
      khbzsx=mumpzx+ibs
      khbztx=mumpzx+ibt
      khbcsx=mumpcx+ibs
      khbctx=mumpcx+ibt
c...   retrieve   c(x)
      call getxy(cr,iconf,cr(mumpcx+1),khbcsx,khbctx,mxx,icval,1)
c...   retrieve   z(x)
      call getxy(cr,iconf,cr(mumpzx+1),khbzsx,khbztx,mxx,izval,0)
222   call upack2(iconf(3,ic),lamtc,lamsc)
      call upack2(iconf(1,ic),ke,je)
_IF(cray,convex,t3e,i8,i8drct)
      isc=iconf(4,ic)
_ELSE
      isc=ipad(iconf(4,ic))
_ENDIF
      ns=norbsi(isc)
      nt=norbtr(isc)
      khzsc=khbzsx+je
      khcsc=khbcsx+je
      khztc=khbztx+ke
      khctc=khbctx+ke
333   call upack2(iconf(3,ir),lamtr,lamsr)
      call upack2(iconf(1,ir),ke,je)
_IF(cray,convex,t3e,i8,i8drct)
      isr=iconf(4,ir)
_ELSE
      isr=ipad(iconf(4,ir))
_ENDIF
      norsir=norbsi(isr)
      nortrr=norbtr(isr)
      lamss=lamsc*lamsr
      lamtt=lamtc*lamtr
_IF(cray,convex,t3e,i8,i8drct)
      mcc=jconf(ir+khbjco)
_ELSE
      mcc=ipad(jconf(ir+khbjco))
_ENDIF
      if(mcc.ne.mxx)goto 111
      khzsr=khbzsx+je
      khcsr=khbcsx+je
      khztr=khbztx+ke
      khctr=khbctx+ke
      return
111   if(mcc.eq.myy)goto 444
      call putxy(cr,iconf,khbzsy,khbzty,myy)
      myy=mcc
      ibs=ibiass(mcc)
      ibt=ibiatt(mcc)
      khbzsy=mumpzy+ibs
      khbzty=mumpzy+ibt
      khbcsy=mumpcy+ibs
      khbcty=mumpcy+ibt
c...   retrieve   c(y)
      call getxy(cr,iconf,cr(mumpcy+1),khbcsy,khbcty,myy,icval,1)
c...   retrieve   z(y)
      call getxy(cr,iconf,cr(mumpzy+1),khbzsy,khbzty,myy,izval,0)
444   khzsr=khbzsy+je
      khcsr=khbcsy+je
      khztr=khbzty+ke
      khctr=khbcty+ke
      return
      end
      subroutine getxy(cr,iconf,c,khs,kht,mxy,iblok,jsi)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension c(*),cr(*),iconf(5,*)
      common/junk /icll(255),ichh(255)
INCLUDE(common/symci)
      common/symcon/root2
INCLUDE(common/cdcryz)
INCLUDE(common/disktl)
      il=icll(mxy)
      ih=ichh(mxy)
      nsp=0
_IF(cray,t3e,convex,i8,i8drct)
      je=iand(iconf(2,il),mms)
_ELSEIF(bits8)
      call dand8(iconf(2,il),mms,tmp8)
      je=ipad(tmp8)
_ELSE
      je=ipad(dand(iconf(2,il),mms))
_ENDIF
_IF(ga)
      if (jsi.eq.1) then
         call search(iblok+je,numcz)
      endif
_ELSE
      call search(iblok+je,numci)
_ENDIF
      do 2 loop=il,ih
      call upack2(iconf(3,loop),ntt,nss)
      call upack2(iconf(1,loop),ke,je)
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,loop)
_ELSE
      isi=ipad(iconf(4,loop))
_ENDIF
      if(isi.eq.jsi)nsp=nsp+nss
_IF(ga)
      if (jsi.eq.1) then
         call reads(cr(khs+je+1),nss*norbsi(isi),numcz)
         call reads(cr(kht+ke+1),ntt*norbtr(isi),numcz)
      elseif(jsi.eq.0) then
         call vclr(cr(khs+je+1),1,nss*norbsi(isi))
         call vclr(cr(kht+ke+1),1,ntt*norbtr(isi))
      else
         call caserr('jsi not 0 or 1 in getxy')
      endif
_ELSE
      call reads(cr(khs+je+1),nss*norbsi(isi),numci)
      call reads(cr(kht+ke+1),ntt*norbtr(isi),numci)
_ENDIF
2     continue
      if(jsi)888,999,888
888   call renorm(c,root2,nsp)
999   return
      end
      subroutine putxy(cr,iconf,khs,kht,mxy)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension cr(*),iconf(5,*)
INCLUDE(common/disktl)
      common/pager/iparit,khbjco,nwolf,iclogo,icval,izval,
     *mxx,khbzsx,khbztx,khbcsx,khbctx,mumpzx,mumpcx,
     *myy,khbzsy,khbzty,khbcsy,khbcty,mumpzy,mumpcy
      common/junk /icll(255),ichh(255)
INCLUDE(common/cdcryz)
INCLUDE(common/symci)
      common/symcon/root2,root2i
c
      if(mxy)888,999,888
888   il=icll(mxy)
      ih=ichh(mxy)
_IF(cray,t3e,convex,i8,i8drct)
      je=iand(iconf(2,il),mms)
_ELSEIF(bits8)
      call dand8(iconf(2,il),mms,tmp8)
      je = ipad(tmp8)
_ELSE
      je=ipad(dand(iconf(2,il),mms))
_ENDIF
_IF(ga)
      call search(izval+je,numcz)
_ELSE
      call search(izval+je,numci)
_ENDIF
      do 2 loop=il,ih
      call upack2(iconf(3,loop),ntt,nss)
      call upack2(iconf(1,loop),ke,je)
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,loop)
_ELSE
      isi=ipad(iconf(4,loop))
_ENDIF
_IF(ga)
      if (isi.eq.1)call renorm(cr(khs+je+1),root2i,nss)
      call wrt3s(cr(khs+je+1),nss*norbsi(isi),numcz)
      call wrt3s(cr(kht+ke+1),ntt*norbtr(isi),numcz)
_ELSE
      call wrt3s(cr(khs+je+1),nss*norbsi(isi),numci)
      call wrt3s(cr(kht+ke+1),ntt*norbtr(isi),numci)
_ENDIF
2     continue
999   return
      end
_IF(parallel)
      subroutine genz(cr,iconf,iblkcz)
c MPP
c MPP       parallel for ipsc
c MPP       default numci here is private
c MPP       c/z vectors are on private file
c MPP
_ELSE
      subroutine genz(cr,iconf,iblkcz,h11)
_ENDIF(parallel)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
      integer *8 shiftr
      external shiftr
_ELSE
      REAL iconf
_ENDIF
c... to compute z=h*c where c is a trial ci vector
c... and to obtain self energy of c
      logical lrenor
      logical nodest,modeds,noded,nodev,logic,loop1,msortx
      dimension iblkcz(2)
      dimension cr(*),iconf(5,*)
INCLUDE(common/sizes)
      common/erg/aaabbb(nd200+22),nnnnnn(4),icyc
      common/junk /icll(255),ichh(255),ibiass(255),ibiatt(255),mblokk,
     *mump(4),icl0(1021),icl1(1021),icl2(1021)
      common/pager/iparit,khbjco,nwolf,iclogo,icval,izval,
     *mxx,khbzsx,khbztx,khbcsx,khbctx,mumpzx,mumpcx,
     *myy,khbzsy,khbzty,khbcsy,khbcty,mumpzy,mumpcy
      common/spcoef/nwwwww(9),iwwwww(9)
INCLUDE(common/timanb)
       common /lsort/ jbass3(8)
INCLUDE(common/cntrl)
INCLUDE(common/diactl)
INCLUDE(common/ccntl)
INCLUDE(common/disktl)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88
_IF(ga)
     *,khb111,khb222,khb333
     *,khbiii,khbkkk,khbddd
_ENDIF
      common/spew/rmode(5)
INCLUDE(common/cdcryz)
INCLUDE(common/cdcryi)
INCLUDE(common/table)
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt,
     *khcsc,khzsc,khctc,khztc,isc,khcsr,khzsr,khctr,khztr,isr,
     *norsir,nortrr,lamstc
      common/symcon/root2,root2i
INCLUDE(common/symci)
INCLUDE(common/auxh)
INCLUDE(common/helpr)
INCLUDE(common/outctl)
INCLUDE(common/disc)
_IF(ga)
INCLUDE(common/corctl)
INCLUDE(common/orb)
_ENDIF
_IF(parallel)
c MPP
      logical oipsci,odo,odont,oroot
      parameter (maxsp=10)
c MPP    maxsp is max # c-vectors to be treated at a time in ssabcd
c MPP
_ENDIF(parallel)
_IF(ga)
      common/rem_iabc/myext(62,nd200),myibass3(8,nd200),myniky(nd200,8),
     *mynorbe(nd200,8),iabase(62),opcivect
      logical ofailiabc,opcivect
_ENDIF
_IF(linux)
      external fget
_ENDIF
      equivalence (nspnsp,nspips(1)),(nornor,norbsi(1))
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      equivalence (rmode(4),kkkbas),(rmode(5),icfock)
       loijka=.false.
       loijab=.false.
       loiajb=.false.
       modeds=moded.and.modest
       tmlast=cpulft(1)
_IF(parallel)
c MPP
c MPP          switch getsb etc numci useage to private
c MPP          all common file accesses will be explicit
c MPP
        call pried7
c MPP
_ENDIF(parallel)
      call brtsb(index(1))
      call getsb(rmode(1),1)
      icyc=icyc+1
      if(instor)goto 887
c-----------------out of store version----------------------------------
      kh=0
      icval=iblkcz(1)
      izval=iblkcz(2)
_IF(cray,convex,t3e,i8,i8drct)
      icdoub=icval+iconf(2,1)
      izdoub=izval+iconf(2,1)
_ELSE
      icdoub=icval+ipad(iconf(2,1))
      izdoub=izval+ipad(iconf(2,1))
_ENDIF
      nodev=.not.modev
      noded=.not.moded
      nodest=.not.modest
       msortx=msort
92121 nodel=model
      goto (99999,90002,90003,90004,90005,90006,90007,90008,90009,
     *90010,90011,90012,90013),model
c...
c... vacuum / vacuum (ij/kl)
c...
90002 continue
      call vclr(cr(khbill+1),1,nval)
      if(modev)goto 92002
_IF(ga)
      call rdedx(cr(lump1+1),nval,icval,numcz)
_ELSE
      call rdedx(cr(lump1+1),nval,icval,numci)
_ENDIF
      khbc=lump1
      khbz=khbill
      khbb=lump2
      call vvijkl(cr,iconf)
_IF(ga)
91002 call wrt3(cr(khbz+1),nval,izval,numcz)
_ELSE
91002 call wrt3(cr(khbz+1),nval,izval,numci)
_ENDIF
90888 call timbrk
      goto 92121
_IF(ga)
92002 call wrt3(cr(khbill+1),nval,izval,numcz)
_ELSE
92002 call wrt3(cr(khbill+1),nval,izval,numci)
_ENDIF
93131    call possb(nwwwww(nodel-1),iwwwww(nodel-1))
      call getsb(rmode(1),1)
      goto 90888
c...
c... doublet / doublet (ij/kl)
c...
90003 continue
      call vclr(cr(khbill+1),1,ndoub)
       if(moded)goto 92003
_IF(ga)
      call rdedx(cr(lump3+1),ndoub,icdoub,numcz)
_ELSE
      call rdedx(cr(lump3+1),ndoub,icdoub,numci)
_ENDIF
      khbcd=lump3
      khbzd=khbill
      khbb=lump4
      call ddijkl(cr,iconf)
_IF(ga)
91003 call wrt3(cr(khbzd+1),ndoub,izdoub,numcz)
_ELSE
91003 call wrt3(cr(khbzd+1),ndoub,izdoub,numci)
_ENDIF
      goto 90888
_IF(ga)
92003 call wrt3(cr(khbill+1),ndoub,izdoub,numcz)
_ELSE
92003 call wrt3(cr(khbill+1),ndoub,izdoub,numci)
_ENDIF
      goto 93131
c...
c...  n-2 / n-2  (ij/kl)
c...
90004 continue
      call vclr(cr(khbill+1),1,nbig)
_IF(cray,t3e,convex,i8,i8drct)
      jc=iand(iconf(2,nstrt2),mms)
_ELSEIF(bits8)
      call dand8(iconf(2,nstrt2),mms,tmp8)
      jc=ipad(tmp8)
_ELSE
      jc=ipad(dand(iconf(2,nstrt2),mms))
_ENDIF
_IF(ga)
      call search(izval+jc,numcz)
_ELSE
      call search(izval+jc,numci)
_ENDIF
      do 91004 loop=nstrt2,nlast2
      call upack2(iconf(3,loop),lamt,lams)
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,loop)
_ELSE
      isi=ipad(iconf(4,loop))
_ENDIF
_IF(ga)
      call wrt3s(cr(khbill+1),lams*norbsi(isi),numcz)
91004 call wrt3s(cr(khbill+1),lamt*norbtr(isi),numcz)
_ELSE
      call wrt3s(cr(khbill+1),lams*norbsi(isi),numci)
91004 call wrt3s(cr(khbill+1),lamt*norbtr(isi),numci)
_ENDIF
      if(modest)goto 93131
         if(.not.msortx)call eblock(cr,iconf,icl0,mgot0)
      khbb=lump27
      khbjco=lump26
      nwolf=2
      call setpag(cr,icl0,lump29)
94004 call pagexy(cr,iconf,cr)
      call getsb(cr(khbb+1),lamss+lamtt)
      call ssijkl(cr,cr(khzsc+1),cr(khztc+1),cr(khcsc+1),cr(khctc+1),
     *cr(khzsr+1),cr(khztr+1),cr(khcsr+1),cr(khctr+1))
      call getsb(rmode(1),1)
      if(model-4)95004,94004,95004
95004 call putxy(cr,iconf,khbzsx,khbztx,mxx)
      call putxy(cr,iconf,khbzsy,khbzty,myy)
      goto 90888
c...
c... doublet / vacuum (ij/ka) and fock(ia)
c...
_IF(ga)
90005 call vclr(cr(lump5+1),1,nvd)
_ELSE
90005 call rdedx(cr(lump5+1),nvd,izval,numci)
_ENDIF
      khbz=lump5
      khbzd=lump6
      khbc=lump7
      khbcd=lump8
      khbb=khbill
_IF(ga)
      if(nodev)call rdedx(cr(khbc+1),nval,icval,numcz)
      if(noded)call rdedx(cr(khbcd+1),ndoub,icdoub,numcz)
_ELSE
      if(nodev)call rdedx(cr(khbc+1),nval,icval,numci)
      if(noded)call rdedx(cr(khbcd+1),ndoub,icdoub,numci)
_ENDIF
      call dvijka(cr,iconf)
_IF(ga)
      if(noded)call wrt3(cr(khbz+1),nval,izval,numcz)
_ELSE
      if(noded)call wrt3(cr(khbz+1),nval,izval,numci)
_ENDIF
      if(nodev)goto 91003
      goto 90888
c...
c... n-2 / doublet (ij/ka) and fock(ia)
c...
90006 if(modeds)goto 93131
      khbzd=lump51
      khbcd=lump52
      khzsc=lump53
      khcsc=lump54
      call getsb(rmode(2),3)
_IF(ga)
      if(nodest)call vclr(cr(khbzd+1),1,ndoub)
      if(noded)call rdedx(cr(khbcd+1),ndoub,icdoub,numcz)
_ELSE
      if(nodest)call rdedx(cr(khbzd+1),ndoub,izdoub,numci)
      if(noded)call rdedx(cr(khbcd+1),ndoub,icdoub,numci)
_ENDIF
      if(loijka)goto 61006
      call rdedx(cr(khbill+1),m2e,index(2),numci)
      call brtsc(index(3))
      call getsc(rmode(5),1)
61006 call upack2(iconf(2,ic),kc,jc)
      call upack2(iconf(3,ic),lamtc,lamsc)
_IF(cray,convex,t3e,i8,i8drct)
      isc=iconf(4,ic)
_ELSE
      isc=ipad(iconf(4,ic))
_ENDIF
      ns=norbsi(isc)
      nt=norbtr(isc)
      nss=ns*lamsc
      ntt=nt*lamtc
      mus=nss+ntt
      khztc=khzsc+nss
      khctc=khcsc+nss
      if(modest)goto 69006
_IF(ga)
c rwah in twee rdedx gemaakt ...
c     call rdedx(cr(khcsc+1),mus,icval+jc,numcz)
      call rdedx(cr(khcsc+1),nss,icval+jc,numcz)
      call rdedx(cr(khctc+1),ntt,icval+kc,numcz)
_ELSE
      call rdedx(cr(khcsc+1),mus,icval+jc,numci)
_ENDIF
      if(isc.eq.1)call renorm(cr(khcsc+1),root2,lamsc)
69006 if(moded)goto 73006
      jc=jc+izval
_IF(ga)
      call vclr(cr(khzsc+1),1,mus)
_ELSE
      call rdedx(cr(khzsc+1),mus,jc,numci)
_ENDIF
73006 call sdijka(cr,iconf)
      if(moded)goto 64006
      if(nss)63006,62006,63006
_IF(ga)
63006 if(isc.eq.1)call renorm(cr(khzsc+1),root2i,lamsc)
      call wrt3(cr(khzsc+1),nss,jc,numcz)
_ELSE
63006 call wrt3(cr(khzsc+1),nss,jc,numci)
_ENDIF
      if(ntt)62006,64006,62006
_IF(ga)
62006 call wrt3(cr(khztc+1),ntt,izval+kc,numcz)
_ELSE
62006 call wrt3(cr(khztc+1),ntt,izval+kc,numci)
_ENDIF
64006 if(model-6)68006,61006,68006
68006 if(nodest)goto 91003
      goto 90888
c...
c... doublet / doublet (ij/ab) + (ia/ib) + fock(ab)
c...
90007 if(moded)goto 93131
_IF(ga)
      call rdedx(cr(lump10+1),ndoub,icdoub,numcz)
      call vclr(cr(lump9+1),1,ndoub)
_ELSE
      call rdedx(cr(lump10+1),ndoub,icdoub,numci)
      call rdedx(cr(lump9+1),ndoub,izdoub,numci)
_ENDIF
      khbzd=lump9
      khbcd=lump10
      khbb=khbill
      call ddijab(cr,iconf)
      goto 91003
c...
c...  n-2 / n-2  (ij/ab) + (ia/ib) + fock(ab)
c...
90008 if(modest)goto 93131
      if(msortx)goto 96008
      call eblock(cr,iconf,icl1,mgot1)
      goto 95008
96008 if(loijab)goto 86008
95008 call rdedx(cr(khbill+1),m2fock,index(4),numci)
      call reads(cr(lump19+1),m2coul,numci)
86008 call brtsc(index(6))
      nwolf=3
      khbjco=lump33
      call setpag(cr,icl1,lump38)
81008 call pagexy(cr,iconf,cr)
      if(kkkbas)84008,83008,84008
c... fock operator
83008 call ssfock(cr)
      goto 82008
84008 kkkbas=kkkbas+kh
      call getsb(cr(khb3+1),lamss)
      call getsb(cr(khb4+1),lamtt)
      if(ic.ne.ir)goto 85008
c... (ia/ib)
      call getsb(cr(khb5+1),lamsc*lamtc)
      call ssiaib(cr)
      goto 82008
c... (ij/ab)
85008 call ssijab(cr)
82008 call getsb(rmode(1),1)
      if(model-8)95004,81008,95004
c...
c...  doublet / doublet (ia/jb)
c...
90009 if(moded)goto 91009
_IF(ga)
      call rdedx(cr(lump12+1),ndoub,icdoub,numcz)
      call vclr(cr(lump11+1),1,ndoub)
_ELSE
      call rdedx(cr(lump12+1),ndoub,icdoub,numci)
      call rdedx(cr(lump11+1),ndoub,izdoub,numci)
_ENDIF
91009 khbzd=lump11
      khbcd=lump12
      khbb=khbill
      call ddiajb(cr,iconf)
      if(noded)goto 91003
      goto 90888
c...
c... n-2  /  n-2  (ia/jb)
c...
90010 if(modest)goto 93131
      if(msortx)goto 71010
      call eblock(cr,iconf,icl2,mgot2)
      goto 61010
71010 if(loiajb)goto 51010
61010 call rdedx(cr(khbill+1),m2exch,index(8),numci)
51010 nwolf=4
      khbjco=lump99
      call setpag(cr,icl2,lump47)
81010 call pagexy(cr,iconf,cr)
      call getsb(cr(khb1+1),lamss)
      call getsb(cr(khb44+1),lamtt)
      call getsb(cr(khb55+1),lamsc*lamtr)
      call getsb(cr(khb66+1),lamsr*lamtc)
      lamstc=lamsc+lamtc
      call ssiajb(cr)
      call getsb(rmode(1),1)
      if(model-10)95004,81010,95004
c...
c...  n-2 / vacuum (ia/jb)
c...
90011   call getsb(rmode(2),3)
        call rdedx(cr(khbill+1),m2exch,index(5),numci)
        khbz=lump16
        khbc=lump17
_IF(ga)
      if(nodest)call vclr(cr(khbz+1),1,nval)
      if(nodev)call rdedx(cr(khbc+1),nval,icval,numcz)
_ELSE
      if(nodest)call rdedx(cr(khbz+1),nval,izval,numci)
      if(nodev)call rdedx(cr(khbc+1),nval,icval,numci)
_ENDIF
91011 continue
_IF(cray,convex,t3e,i8,i8drct)
      itop=iconf(4,ic)
_ELSE
      itop=ipad(iconf(4,ic))
_ENDIF
      ns=norbsi(itop)
      nt=norbtr(itop)
      call upack2(iconf(3,ic),lamtc,lamsc)
      call upack2(iconf(2,ic),kc,jc)
      nss=ns*lamsc
      ntt=nt*lamtc
      mus=nss+ntt
      khbztc=lump14+nss
      if(modest)goto 94011
_IF(ga)
c rwah twee rdedx-en van gemaakt ...
c     call rdedx(cr(lump15+1),mus,icval+jc,numcz)
      call rdedx(cr(lump15+1),nss,icval+jc,numcz)
      call rdedx(cr(lump15+nss+1),ntt,icval+kc,numcz)
_ELSE
      call rdedx(cr(lump15+1),mus,icval+jc,numci)
_ENDIF
      if(itop.eq.1)call renorm(cr(lump15+1),root2i,lamsc)
94011 jc=jc+izval
_IF(ga)
      if(nodev)call vclr(cr(lump14+1),1,mus)
_ELSE
      if(nodev)call rdedx(cr(lump14+1),mus,jc,numci)
_ENDIF
      call sviajb(cr,iconf,cr(lump14+1),cr(khbztc+1),cr(lump15+1),
     *cr(lump15+nss+1),cr(lump18+1))
      if(modev)goto 95011
      if(nss)92011,93011,92011
_IF(ga)
92011 if(itop.eq.1)call renorm(cr(lump14+1),root2i,lamsc)
      call wrt3(cr(lump14+1),nss,jc,numcz)
_ELSE
92011 call wrt3(cr(lump14+1),nss,jc,numci)
_ENDIF
      if(ntt)93011,95011,93011
_IF(ga)
93011 call wrt3(cr(khbztc+1),ntt,izval+kc,numcz)
_ELSE
93011 call wrt3(cr(khbztc+1),ntt,izval+kc,numci)
_ENDIF
95011   if(model.eq.11)goto 91011
96011   if(nodest)goto 91002
        goto 90888
c...
c...  n-2  /  vacuum  (ia/ib)
c...
90012 call getsb(rmode(2),3)
      call rdedx(cr(khbill+1),m2fock,index(4),numci)
      khbz=lump21
      khbc=lump22
_IF(ga)
      if(nodest)call vclr(cr(khbz+1),1,nval)
      if(nodev)call rdedx(cr(khbc+1),nval,icval,numcz)
_ELSE
      if(nodest)call rdedx(cr(khbz+1),nval,izval,numci)
      if(nodev)call rdedx(cr(khbc+1),nval,icval,numci)
_ENDIF
91012 continue
_IF(cray,t3e,convex,i8,i8drct)
      itop=iconf(4,ic)
      lamsc=iand(iconf(3,ic),mms)
_ELSEIF(bits8)
      itop=ipad(iconf(4,ic))
      call dand8(iconf(3,ic),mms,tmp8)
      lamsc=ipad(tmp8)
_ELSE
      itop=ipad(iconf(4,ic))
      lamsc=ipad(dand(iconf(3,ic),mms))
_ENDIF
      nss=norbsi(itop)*lamsc
_IF(cray,t3e,convex,i8,i8drct)
      jc=iand(iconf(2,ic),mms)
_ELSEIF(bits8)
      call dand8(iconf(2,ic),mms,tmp8)
      jc=ipad(tmp8)
_ELSE
      jc=ipad(dand(iconf(2,ic),mms))
_ENDIF
      if(modest)goto 93012
_IF(ga)
      call rdedx(cr(lump20+1),nss,icval+jc,numcz)
_ELSE
      call rdedx(cr(lump20+1),nss,icval+jc,numci)
_ENDIF
      if(itop.eq.1)call renorm(cr(lump20+1),root2i,lamsc)
93012 jc=jc+izval
_IF(ga)
      if(nodev)call vclr(cr(lump19+1),1,nss)
_ELSE
      if(nodev)call rdedx(cr(lump19+1),nss,jc,numci)
_ENDIF
      call sviaib(cr,iconf,cr(lump19+1),cr(lump20+1),cr(lump23+1))
_IF(ga)
      if (nodev.and.itop.eq.1)call renorm(cr(lump19+1),root2i,lamsc)
      if(nodev)call wrt3(cr(lump19+1),nss,jc,numcz)
_ELSE
      if(nodev)call wrt3(cr(lump19+1),nss,jc,numci)
_ENDIF
      if(model.eq.12)goto 91012
      goto 96011
c...
c...   n-2  /  doublet   (ia/bc)
c...
90013 if(modeds) goto 80014
      khzsc=lump55
      khcsc=lump56
      khbzd=lump57
      khbcd=lump58
      icfock=0
_IF(ga)
      if (opcivect)call getmyminmax(myminb,mymaxb)
      if(nodest)call vclr(cr(khbzd+1),1,ndoub)
      if(noded)call rdedx(cr(khbcd+1),ndoub,icdoub,numcz)
_ELSE
      if(nodest)call rdedx(cr(khbzd+1),ndoub,izdoub,numci)
      if(noded)call rdedx(cr(khbcd+1),ndoub,icdoub,numci)
_ENDIF
_IF(parallel)
c MPP       index(9) on common file
c         call comed7
_ENDIF(parallel)
      call brtsc(index(9))
_IF(parallel)
c         call pried7
c get all iabc integrals in core
      memuseint=0
      ofailiabc=.false.
      do ii=1,nint
          isi=iro(ii)
          iabase(ii)=memuseint+lump59
          memuseint=memuseint+myniky(ii,isi)
          if (memuseint+lump59.gt.ngot) ofailiabc=.true.
          if (.not.ofailiabc) then
             call getsc(cr(iabase(ii)+1),myniky(ii,isi))
          endif
      enddo
      if (ofailiabc) then
         write(6,'(A,I3,I15)')'Core needed for iabc: ',
     1                         ipg_nodeid(),memuseint
         write(6,'(A,I3,I15)')'Available: ',
     1                         ipg_nodeid(),ngot-lump59
         call caserr('Not enough core to hold all iabc integrals')
      endif
_ENDIF(parallel)
      call getsb(rmode(2),3)
_IF(ga)
      goto 53013
_ENDIF
51013 itop=kkkbas-icfock
      if(itop.eq.0) go to 53013
_IF(parallel)
c MPP       
c         call comed7
_ENDIF(parallel)
      do 54013 loop=1,itop
      isi=iro(loop+icfock)
_IF(ga)
54013 call getsc(cr(lump59+1),myniky(kkkbas,isi))
_ELSE
54013 call getsc(cr(khbill+1),niky(isi))
_ENDIF
_IF(parallel)
c         call pried7
c MPP       
_ENDIF(parallel)
_IF(ga)
      do 1 loop=1,nirr
    1 jbass3(loop)=myibass3(loop,kkkbas)+khbill
_ELSE
      do 1 loop=1,nirr
    1 jbass3(loop)=ibass3(loop,isi)+khbill
_ENDIF
      icfock=kkkbas
53013 call upack2(iconf(3,ic),lamtc,lamsc)
      call upack2(iconf(2,ic),kc,jc)
_IF(ga)
      if (opcivect) then
         odo=.false.
         if (myminb.le.jc+1.and.jc+1.le.mymaxb)odo=.true.
         if(.not.odo)goto 57000
      else
         odo=.true.
      endif
_ENDIF
_IF(cray,convex,t3e,i8,i8drct)
      itop=iconf(4,ic)
_ELSE
      itop=ipad(iconf(4,ic))
_ENDIF
      ns=norbsi(itop)
      nt=norbtr(itop)
      nss=ns*lamsc
      ntt=nt*lamtc
      khztc=khzsc+nss
      khctc=khcsc+nss
      mus=nss+ntt
      if(modest)goto 59013
_IF(ga)
c rwah twee rdedx-en ...
c     call rdedx(cr(khcsc+1),mus,icval+jc,numcz)
      call rdedx(cr(khcsc+1),nss,icval+jc,numcz)
      call rdedx(cr(khcsc+nss+1),ntt,icval+kc,numcz)
_ELSE
      call rdedx(cr(khcsc+1),mus,icval+jc,numci)
_ENDIF
      if(itop.eq.1)call renorm(cr(khcsc+1),root2i,lamsc)
59013 if(moded)goto 85013
      jc=jc+izval
_IF(ga)
      call vclr(cr(khzsc+1),1,mus)
_ELSE
      call rdedx(cr(khzsc+1),mus,jc,numci)
_ENDIF
85013 call sdiabc(cr,iconf)
         if(moded)goto 58013
_IF(ga)
      if(itop.eq.1)call renorm(cr(khzsc+1),root2i,lamsc)
      call wrt3(cr(khzsc+1),nss,jc,numcz)
      call wrt3(cr(khztc+1),ntt,kc+izval,numcz)
_ELSE
      call wrt3(cr(khzsc+1),nss,jc,numci)
      call wrt3(cr(khztc+1),ntt,kc+izval,numci)
_ENDIF
_IF(ga)
c odo is false
57000 if (.not.odo) then
         lamst=lamsc+lamtc
         icl=ic
_IF(cray,i8)
57001    lamr=iconf(3,ir)
_ELSE
57001    lamr=ipad(iconf(3,ir))
_ENDIF
         call skipsb(lamr*lamst)
         call getsb(rmode(1),1)
         if (model.ne.13) goto 58013
         call getsb(rmode(2),3)
         if (ic.ne.icl)goto 58013
         goto 57001         
      endif
58013 if(model-13)68006,53013,68006
_ELSE
58013 if(model-13)68006,51013,68006
_ENDIF
c...
c...   n-2   /   n-2 (ab/cd)
c...
_IF(parallel)
c MPP
c MPP    ab/cd integrals are on common file
c MPP    vectors are on private file
c MPP    changes try to disturb code as little as possible
c MPP    all have to read sym-1 z-vector to renormalise
c MPP    number of c-vectors is limited to spread effort
c MPP    default numci is common file
c MPP
_ENDIF(parallel)
99999 if(modest) goto 80014
_IF(parallel)
c MPP
c       iflop = iipsci()
c       call comed7
        call pried7
c MPP
_ENDIF(parallel)
      is=nstrt2
      it=is
      do 99014 loop=1,nirr
c...     singlets
      loop1=loop.eq.1
      nss=norbsi(loop)
      logic=.true.
      indexx=index(loop+9)
91014 ibasz=lump24
      ibasc=lump25
      nsp=0
_IF(parallel)
c MPP
c      odont = oipsci()
c      odo = .not.odont
c MPP
_ENDIF(parallel)
      do 98014 noop=is,nlast2
      isx=noop
_IF(cray,t3e,convex,i8,i8drct)
      if(iconf(4,noop).ne.loop)goto 97014
      lams=iand(iconf(3,noop),mms)
_ELSEIF(bits8)
      if(ipad(iconf(4,noop)).ne.loop)goto 97014
      call dand8(iconf(3,noop),mms,tmp8)
      lams = ipad(tmp8)
_ELSE
      if(ipad(iconf(4,noop)).ne.loop)goto 97014
      lams=ipad(dand(iconf(3,noop),mms))
_ENDIF
      ntt=lams*nss
      if(ntt.eq.0) go to 98014
      ib=ibasz+ntt
_IF(parallel)
c MPP
      if(ib.ge.lump25.or.(nsp+lams.gt.maxsp.and.nsp.ne.0))goto 95014
c MPP
_ELSE
      if(ib.ge.lump25)goto 95014
_ENDIF
      nsp=nsp+lams
_IF(cray,t3e,convex,i8,i8drct)
      jc=iand(iconf(2,noop),mms)
_ELSEIF(bits8)
      call dand8(iconf(2,noop),mms,tmp8)
      jc=ipad(tmp8)
_ELSE
      jc=ipad(dand(iconf(2,noop),mms))
_ENDIF
_IF(parallel)
c MPP
_IF(ga)
c     if (odo.or.loop1) call vclr(cr(ibasz+1),1,ntt)
c     if (odo) call rdedx(cr(ibasc+1),ntt,icval+jc,numcz)
      call vclr(cr(ibasz+1),1,ntt)
      call rdedx(cr(ibasc+1),ntt,icval+jc,numcz)
_ELSE
      if (odo.or.loop1) call rdedx(cr(ibasz+1),ntt,izval+jc,numcip)
      if (odo) call rdedx(cr(ibasc+1),ntt,icval+jc,numcip)
_ENDIF
c MPP
_ELSE
      call rdedx(cr(ibasz+1),ntt,izval+jc,numci)
      call rdedx(cr(ibasc+1),ntt,icval+jc,numci)
_ENDIF
      ibasz=ib
      ibasc=ibasc+ntt
98014 continue
      isx=nlst21
97014 logic=.false.
95014 if(nsp.eq.0) go to 93014
_IF(parallel)
c MPP
c     if(loop1.and.odo)call renorm(cr(lump25+1),root2i,nsp)
c     if (odo)
c    #call ssabcd(cr(lump24+1),cr(khbill+1),cr(lump25+1),nsp,nss,indexx)
      if(loop1)call renorm(cr(lump25+1),root2i,nsp)
      call ssabcd(cr(lump24+1),cr(khbill+1),cr(lump25+1),nsp,nss,indexx)
c MPP
_ELSE
      if(loop1)call renorm(cr(lump25+1),root2i,nsp)
      call ssabcd(cr(lump24+1),cr(khbill+1),cr(lump25+1),nsp,nss,indexx)
_ENDIF
      if(loop1)call renorm(cr(lump24+1),root2i,nsp)
      mus=isx-1
      ibasz=lump24
      do 92014 noop=is,mus
_IF(cray,t3e,convex,i8,i8drct)
      lams=iand(iconf(3,noop),mms)
_ELSEIF(bits8)
      call dand8(iconf(3,noop),mms,tmp8)
      lams = ipad(tmp8)
_ELSE
      lams=ipad(dand(iconf(3,noop),mms))
_ENDIF
      if(lams.eq.0)go to 92014
_IF(cray,t3e,convex,i8,i8drct)
      jc=iand(iconf(2,noop),mms)
_ELSEIF(bits8)
      call dand8(iconf(2,noop),mms,tmp8)
      jc=ipad(tmp8)
_ELSE
      jc=ipad(dand(iconf(2,noop),mms))
_ENDIF
      ntt=nss*lams
_IF(parallel)
c MPP
_IF(ga)
c     if (loop1.or.odo) call wrt3(cr(ibasz+1),ntt,izval+jc,numcz)
c     if (odo) call wrt3(cr(ibasz+1),ntt,izval+jc,numcz)
      call wrt3(cr(ibasz+1),ntt,izval+jc,numcz)
_ELSE
      if (loop1.or.odo) call wrt3(cr(ibasz+1),ntt,izval+jc,numcip)
_ENDIF
c MPP
_ELSE
      call wrt3(cr(ibasz+1),ntt,izval+jc,numci)
_ENDIF
      ibasz=ibasz+ntt
92014 continue
93014 is=isx
      if(logic)goto 91014
c... triplets
      nss=norbtr(loop)
      logic=.true.
      indexx=index(loop+17)
81014 ibasz=lump24
      ibasc=lump25
      nsp=0
_IF(parallel)
c MPP
c      odont = oipsci()
c      odo = .not.odont
c MPP
_ENDIF(parallel)
      do 88014 noop=it,nlast2
      isx=noop
_IF(cray,t3e,convex,i8,i8drct)
      if(iconf(4,noop).ne.loop)goto 87014
_IF(cray,t3e,i8drct)
      lams=shiftr(iconf(3,noop),n32)
_ELSE
      lams=ishft(iconf(3,noop),-n32)
_ENDIF
_ELSE
      if(ipad(iconf(4,noop)).ne.loop)goto 87014
_IF(bits8)
      call shiftr8(iconf(3,noop),n32,tmp8)
      lams=ipad(tmp8)
_ELSE
      lams=ipad(shiftr(iconf(3,noop),n32))
_ENDIF
_ENDIF
      ntt=lams*nss
      if(ntt)86014,88014,86014
86014 ib=ibasz+ntt
_IF(parallel)
c MPP
      if(ib.ge.lump25.or.(nsp+lams.gt.maxsp.and.nsp.ne.0))goto 85014
c MPP
_ELSE
      if(ib.ge.lump25)goto 85014
_ENDIF
      nsp=nsp+lams
_IF(cray,t3e,convex,i8,i8drct)
_IF(cray,t3e,i8drct)
      jc=shiftr(iconf(2,noop),n32)
_ELSE
      jc=ishft(iconf(2,noop),-n32)
_ENDIF
_ELSEIF(bits8)
      call shiftr8(iconf(2,noop),n32,tmp8)
      jc=ipad(tmp8)
_ELSE
      jc=ipad(shiftr(iconf(2,noop),n32))
_ENDIF
_IF(parallel)
c MPP
_IF(ga)
c     if (odo) call vclr(cr(ibasz+1),1,ntt)
c     if (odo) call rdedx(cr(ibasc+1),ntt,icval+jc,numcz)
      call vclr(cr(ibasz+1),1,ntt)
      call rdedx(cr(ibasc+1),ntt,icval+jc,numcz)
_ELSE
      if (odo) call rdedx(cr(ibasz+1),ntt,izval+jc,numcip)
      if (odo) call rdedx(cr(ibasc+1),ntt,icval+jc,numcip)
_ENDIF
c MPP
_ELSE
      call rdedx(cr(ibasz+1),ntt,izval+jc,numci)
      call rdedx(cr(ibasc+1),ntt,icval+jc,numci)
_ENDIF
      ibasz=ib
      ibasc=ibasc+ntt
88014 continue
      isx=nlst21
87014 logic=.false.
85014 if(nsp)84014,83014,84014
_IF(parallel)
c MPP
c84014 if (odo)
c     #call ssabcd(cr(lump24+1),cr(khbill+1),cr(lump25+1),nsp,nss,indexx)
84014 call ssabcd(cr(lump24+1),cr(khbill+1),cr(lump25+1),nsp,nss,indexx)
c MPP
_ELSE
84014 call ssabcd(cr(lump24+1),cr(khbill+1),cr(lump25+1),nsp,nss,indexx)
_ENDIF
      mus=isx-1
      ibasz=lump24
      do 82014 noop=it,mus
_IF(cray,t3e,convex,i8,i8drct)
_IF(cray,t3e,i8drct)
      lams=shiftr(iconf(3,noop),n32)
_ELSE
      lams=ishft(iconf(3,noop),-n32)
_ENDIF
_ELSEIF(bits8)
      call shiftr8(iconf(3,noop),n32,tmp8)
      lams = ipad(tmp8)
_ELSE
      lams=ipad(shiftr(iconf(3,noop),n32))
_ENDIF
      if(lams.eq.0)go to 82014
_IF(cray,t3e,convex,i8,i8drct)
_IF(cray,t3e,i8drct)
       jc=shiftr(iconf(2,noop),n32)
_ELSE
       jc=ishft(iconf(2,noop),-n32)
_ENDIF
_ELSEIF(bits8)
      call shiftr8(iconf(2,noop),n32,tmp8)
      jc = ipad(tmp8)
_ELSE
      jc=ipad(shiftr(iconf(2,noop),n32))
_ENDIF
      ntt=nss*lams
_IF(parallel)
c MPP
_IF(ga)
      call wrt3(cr(ibasz+1),ntt,izval+jc,numcz)
c     if (odo) call wrt3(cr(ibasz+1),ntt,izval+jc,numcz)
_ELSE
      if (odo) call wrt3(cr(ibasz+1),ntt,izval+jc,numcip)
_ENDIF
c MPP
_ELSE
      call wrt3(cr(ibasz+1),ntt,izval+jc,numci)
_ENDIF
      ibasz=ibasz+ntt
82014 continue
83014 it=isx
      if(logic)goto 81014
99014 continue
      nodel=14
      call timbrk
_IF(parallel)
c MPP  h11-computation has been moved out of genz
c MPP     reset numci to private
        call pried7
c MPP
c MPP     add z-vectors and compress
c MPP     also compress c-vector
c MPP     then calculate h11 (paralllel)
79014     continue
c MPP
      return
c MPP
_ELSE
c... compute h11
      ntt=ntotal
79014 h11=0.0d0
99114 if(ntt)99111,777,99111
99111 call search(icval,numci)
      mus=khbill
      i29x=icval
99112 call fget(cr(mus+1),m,numci)
      icval=icval+1
      ntt=ntt-m
      mus=mus+m
      if((mus+511).le.khgot2.and.ntt.gt.0)goto 99112
      nss=mus-khbill
      call rdedx(cr(mus+1),nss,izval,numci)
      h11=h11+ddot(nss,cr(khbill+1),1,cr(mus+1),1)
      izval=izval+icval-i29x
      goto 99114
_ENDIF
80014 if(nspnsp)30014,31014,30014
c... set singlet z vector normalization to standard form
c
_IF(ga)
30014 continue
_ELSE
30014 do 32014 loop=nstrt2,nlastx
_ENDIF
_IFN(ga)

_IF(cray,t3e,convex,i8,i8drct)
_IF(cray,t3e)
      mus=and(iconf(3,loop),mms)
_ELSE
      mus=iand(iconf(3,loop),mms)
_ENDIF
_ELSEIF(bits8)
      call dand8(iconf(3,loop),mms,tmp8)
      mus = ipad(tmp8)
_ELSE
      mus=ipad(dand(iconf(3,loop),mms))
_ENDIF
      if(mus)33014,32014,33014
33014 nss=mus*nornor
_IF(cray,t3e,convex,i8,i8drct)
      jc=iand(iconf(2,loop),mms)
_ELSEIF(bits8)
      call dand8(iconf(2,loop),mms,tmp8)
      jc = ipad(tmp8)
_ELSE
      jc=ipad(dand(iconf(2,loop),mms))
_ENDIF
      jc=jc+izval
      call rdedx(cr(khbill+1),nss,jc,numci)
      call renorm(cr(khbill+1),root2i,mus)
      call wrt3(cr(khbill+1),nss,jc,numci)
32014 continue
_ENDIF
31014 ntt=nval
      if(noded)ntt=nvd
      goto 79014
c--------------------in store version-----------------------------------
887   kh=khbb-khbill
      lrenor=.true.
      call vclr(cr(khbz+1),1,ntotal)
2121   nodel=model
      goto (9999,2000,3,4,5,6,7,8,9,10,11,12,13),model
c...
c...      vacuum / vacuum (ij/kl) interactions
c...
2000   if(modev)goto  4747
      call vvijkl(cr,iconf)
888   call timbrk
      goto 2121
4747  call possb(nwwwww(nodel-1),iwwwww(nodel-1))
      call getsb(rmode(1),1)
      goto 888
c...
c...       doublet / doublet (ij/kl) interactions
c...
3     if(moded)goto 4747
       call ddijkl(cr,iconf)
      goto 888
c...
c...       n-2 / n-2 (ij/kl) interactions
c...
4      if(modest)goto 4747
        call renorm(cr(khbcs+1),root2,nspnsp)
44    call getsb(rmode(2),2)
      call upack2(iconf(3,ic),lamtc,lamsc)
      call upack2(iconf(3,ir),lamtr,lamsr)
      lamss=lamsc*lamsr
      lamtt=lamtc*lamtr
      call getsb(cr(khbb+1),lamss+lamtt)
_IF(cray,convex,t3e,i8,i8drct)
      itop=iconf(4,ic)
_ELSE
      itop=ipad(iconf(4,ic))
_ENDIF
      ns=norbsi(itop)
      nt=norbtr(itop)
      call upack2(iconf(1,ic),kc,jc)
      call upack2(iconf(1,ir),kr,jr)
      call ssijkl(cr,cr(khbzs+jc+1),cr(khbzt+kc+1),cr(khbcs+jc+1),
     *cr(khbct+kc+1),cr(khbzs+jr+1),cr(khbzt+kr+1),cr(khbcs+jr+1),
     *cr(khbct+kr+1))
      call getsb(rmode(1),1)
      if(model-4)888,44,888
c...
c...      doublet / vacuum (ij/ka) and fock(ia) interactions
c...
5      call dvijka(cr,iconf)
       goto 888
c...
c...     n-2 / doublet (ij/ka) and fock(ia) interactions
c...
6      if(modeds)goto 4747
        call getsb(rmode(2),3)
        if(loijka)goto 61
        call rdedx(cr(khbb+1),m2e,index(2),numci)
        call brtsc(index(3))
        call getsc(rmode(5),1)
61      call upack2(iconf(1,ic),kc,jc)
        khzsc=khbzs+jc
        khcsc=khbcs+jc
        khztc=khbzt+kc
        khctc=khbct+kc
_IF(cray,convex,t3e,i8,i8drct)
        isc=iconf(4,ic)
_ELSE
        isc=ipad(iconf(4,ic))
_ENDIF
       ns=norbsi(isc)
        nt=norbtr(isc)
       call upack2(iconf(3,ic),lamtc,lamsc)
        call sdijka(cr,iconf)
        if(model-6)888,61,888
c...
c...     doublet / doublet (ij/ab)+(ia/ib)+fock(ab)
c...
7      if(moded)goto 4747
        call ddijab(cr,iconf)
        goto 888
c...
c...    n-2 / n-2 (ij/ab)+(ia/ib)+fock(ab)
c...
8      if(modest)goto 4747
       icl=0
c... read in (ia/ib) + (ij/ab) integrals
        if(loijab)goto 86
       call rdedx(cr(khbb+1),m2fock,index(4),numci)
       call reads(cr(khb2+1),m2coul,numci)
c... initiate input of fock operators
86     call brtsc(index(6))
81     call getsb(rmode(2),3)
       if(ic.eq.icl)goto 82
       icl=ic
      call upack2(iconf(1,ic),kc,jc)
      khzsc=khbzs+jc
      khcsc=khbcs+jc
      khztc=khbzt+kc
      khctc=khbct+kc
      call upack2(iconf(3,ic),lamtc,lamsc)
_IF(cray,convex,t3e,i8,i8drct)
      isc=iconf(4,ic)
_ELSE
      isc=ipad(iconf(4,ic))
_ENDIF
      ns=norbsi(isc)
      nt=norbtr(isc)
82    if(kkkbas)84,83,84
c... fock operator contribution
83    call ssfock(cr)
      goto 8080
84     kkkbas=kkkbas+kh
       if(icl.ne.ir)goto 85
c... (ia/ib) contribution
       call getsb(cr(khb3+1),lamsc*lamsc)
       call getsb(cr(khb4+1),lamtc*lamtc)
       call getsb(cr(khb5+1),lamsc*lamtc)
c... khbb (ia/ib) integrals
c... khb2 (ij/ab) integrals
c... khb3 b(singlet-singlet)
c... khb4 b(triplet-triplet)
c... khb5 b(singlet-triplet)
c... khb7 operator in squared form
c... khb8 scratch space
      call ssiaib(cr)
      goto 8080
c... (ij/ab) contribution
85    call upack2(iconf(3,ir),lamtr,lamsr)
       lamss=lamsc*lamsr
      lamtt=lamtc*lamtr
      call getsb(cr(khb3+1),lamss)
      call getsb(cr(khb4+1),lamtt)
_IF(cray,convex,t3e,i8,i8drct)
      isr=iconf(4,ir)
_ELSE
      isr=ipad(iconf(4,ir))
_ENDIF
      call upack2(iconf(1,ir),kc,jc)
      khzsr=khbzs+jc
      khcsr=khbcs+jc
      khztr=khbzt+kc
      khctr=khbct+kc
      norsir=norbsi(isr)
      nortrr=norbtr(isr)
      call ssijab(cr)
8080  call getsb(rmode(1),1)
       if(model-8)888,81,888
c...
c...     doublet / doublet (ia/jb) interactions
c...
9       call ddiajb(cr,iconf)
        goto 888
c...
c...     n-2 / n-2 (ia/jb) interactions
c...
10     if(modest)goto 4747
         icl=0
c... read in (ia/jb) integrals
       if(.not.loiajb)call rdedx(cr(khbb+1),m2exch,index(8),numci)
101   call getsb(rmode(2),4)
      if(ic.eq.icl)goto 102
      icl=ic
      call upack2(iconf(1,ic),kc,jc)
      khzsc=khbzs+jc
      khcsc=khbcs+jc
       khztc=khbzt+kc
       khctc=khbct+kc
       call upack2(iconf(3,ic),lamtc,lamsc)
      lamstc=lamsc+lamtc
_IF(cray,convex,t3e,i8,i8drct)
       isc=iconf(4,ic)
_ELSE
       isc=ipad(iconf(4,ic))
_ENDIF
      ns=norbsi(isc)
      nt=norbtr(isc)
102   call upack2(iconf(3,ir),lamtr,lamsr)
      call getsb(cr(khb1+1),lamsr*lamsc)
      call getsb(cr(khb44+1),lamtr*lamtc)
      call getsb(cr(khb55+1),lamsc*lamtr)
      call getsb(cr(khb66+1),lamsr*lamtc)
c... khb1  =  singlet(r)  -  singlet(c)
c... khb44 =  triplet(r)  -  triplet(c)
c... khb55 =  triplet(r)  -  singlet(c)
c... khb66 =  singlet(r)  -  triplet(c)
c... khb88 =  scratch space
_IF(cray,convex,t3e,i8,i8drct)
      isr=iconf(4,ir)
_ELSE
      isr=ipad(iconf(4,ir))
_ENDIF
      call upack2(iconf(1,ir),kc,jc)
      khzsr=khbzs+jc
      khcsr=khbcs+jc
      khztr=khbzt+kc
      khctr=khbct+kc
      norsir=norbsi(isr)
      nortrr=norbtr(isr)
      call ssiajb(cr)
      call getsb(rmode(1),1)
      if(model-10)888,101,888
c...
c...      n-2 / vacuum  (ia/jb) interactions
c...
11    call getsb(rmode(2),3)
      call rdedx(cr(khbb+1),m2exch,index(5),numci)
      if(lrenor)call renorm(cr(khbcs+1),0.5d0,nspnsp)
      lrenor=.false.
 111   continue
_IF(cray,convex,t3e,i8,i8drct)
      itop=iconf(4,ic)
_ELSE
      itop=ipad(iconf(4,ic))
_ENDIF
      ns=norbsi(itop)
      nt=norbtr(itop)
      call upack2(iconf(3,ic),lamtc,lamsc)
      call upack2(iconf(1,ic),kc,jc)
      call sviajb(cr,iconf,cr(khbzs+jc+1),cr(khbzt+kc+1),cr(khbcs+jc+1),
     *cr(khbct+kc+1),cr(khb1+1))
      if(model-11)888,111,888
c...
c...      singlet / vacuum (ia/ib) interactions
c...
12    call getsb(rmode(2),3)
      if(lrenor)call renorm(cr(khbcs+1),0.5d0,nspnsp)
       lrenor=.false.
      call rdedx(cr(khbb+1),m2fock,index(4),numci)
 112  continue
_IF(cray,t3e,convex,i8,i8drct)
      lamsc=iand(iconf(3,ic),mms)
      jc=iand(iconf(1,ic),mms)
_ELSEIF(bits8)
      call dand8(iconf(3,ic),mms,tmp8)
      lamsc = ipad(tmp8)
      call dand8(iconf(1,ic),mms,tmp8)
      jc = ipad(tmp8)
_ELSE
      lamsc=ipad(dand(iconf(3,ic),mms))
      jc=ipad(dand(iconf(1,ic),mms))
_ENDIF
      call sviaib(cr,iconf,cr(khbzs+jc+1),cr(khbcs+jc+1),cr(khb2+1))
      if(model-12)888,112,888
c...
c...      n-2 / doublet (ia/bc) interactons
c...
13    if(modeds)goto 9996
      icfock=0
_IF(parallel)
c MPP      iabc ints from common file
         call comed7
_ENDIF(parallel)
       call brtsc(index(9))
_IF(parallel)
         call pried7
c MPP
_ENDIF(parallel)
       if(lrenor)call renorm(cr(khbcs+1),0.5d0,nspnsp)
       lrenor=.false.
       call getsb(rmode(2),3)
131    itop=kkkbas-icfock
        if(itop.eq.0) go to 132
_IF(parallel)
c MPP      
         call comed7
_ENDIF(parallel)
         do 134 loop=1,itop
         isi=iro(loop+icfock)
134      call getsc(cr(khbb+1),niky(isi))
_IF(parallel)
         call pried7
c MPP
_ENDIF(parallel)
        do 2 loop=1,nirr
    2   jbass3(loop)=ibass3(loop,isi)+khbb
        icfock=kkkbas
132   continue
_IF(cray,convex,t3e,i8,i8drct)
      itop=iconf(4,ic)
_ELSE
      itop=ipad(iconf(4,ic))
_ENDIF
      ns=norbsi(itop)
      nt=norbtr(itop)
      call upack2(iconf(3,ic),lamtc,lamsc)
      call upack2(iconf(1,ic),kc,jc)
      khzsc=khbzs+jc
      khcsc=khbcs+jc
      khztc=khbzt+kc
      khctc=khbct+kc
      call sdiabc(cr,iconf)
       if(model.eq.13)goto 131
      call timbrk
9999  if(modest)goto 9996
      
       if(lrenor)call renorm(cr(khbcs+1),0.5d0,nspnsp)
c...
c...      n-2 / n-2  (ab/cd) interactions
c...
_IF(parallel)
c MPP
c MPP        in store for now only divided by symmetry and singlet/triplet
c MPP         if very unevenly divided split (cf out store mode)
c MPP        default numci is common
        call comed7
        iflop = iipsci()
_ENDIF(parallel)
      itop=1
      jtop=1
      do 47474 loop=1,nirr
c... singlets
_IF(parallel)
c MPP
      if (.not.oipsci())
     *call ssabcd(cr(khbzs+itop),cr(khbb+1),cr(khbcs+itop),
_ELSE
      call ssabcd(cr(khbzs+itop),cr(khbb+1),cr(khbcs+itop),
_ENDIF
     *nspips(loop),norbsi(loop),index(loop+9))
c... triplets
_IF(parallel)
      if (.not.oipsci())
     *call ssabcd(cr(khbzt+jtop),cr(khbb+1),cr(khbct+jtop),
_ELSE
      call ssabcd(cr(khbzt+jtop),cr(khbb+1),cr(khbct+jtop),
_ENDIF
     *nspipt(loop),norbtr(loop),index(loop+17))
      itop=itop+nsitot(loop)
47474   jtop=jtop+ntrtot(loop)
_IF(parallel)
c MPP      reset numci to private
        call pried7
c MPP
_ENDIF(parallel)
c... reset singlet c vector to standard normalization
        call renorm(cr(khbcs+1),root2,nspnsp)
       nodel=14
      call timbrk
c... compute h11
_IF(parallel)
c .. this requires checking (not smart to do it before gop)
_ENDIF(parallel)
9996   continue
       call renorm(cr(khbzs+1),root2i,nspnsp)
_IFN(parallel)
       h11=ddot(ntotal,cr(khbc+1),1,cr(khbz+1),1)
_ENDIF(not,parallel)
777    return
       end
      subroutine ddijkl(cr,iconf)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
INCLUDE(common/prnprn)
c...
c... (ij/kl) doublet/doublet interactions
c...
      dimension cr(*),iconf(5,*)
      common/spew/rmode(3)
INCLUDE(common/symci)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
999   call getsb(rmode(2),2)
_IF(cray,t3e,convex,i8,i8drct)
      nc=iconf(3,ic)
      nr=iconf(3,ir)
_ELSE
      nc=ipad(iconf(3,ic))
      nr=ipad(iconf(3,ir))
_ENDIF
       call getsb(cr(khbb+1),nc*nr)
_IF(cray,convex,t3e,i8,i8drct)
      icbas=iconf(1,ic)
      irbas=iconf(1,ir)
      ne=norbe(iconf(4,ic))
_ELSE
      icbas=ipad(iconf(1,ic))
      irbas=ipad(iconf(1,ir))
      ne=norbe(ipad(iconf(4,ic)))
_ENDIF
      call mxmb(cr(irbas+khbcd+1),1,ne,cr(khbb+1),1,nr,cr(icbas+khbzd+1)
     *,1,ne,ne,nr,nc)
      if(ic.ne.ir)call mxmb(cr(icbas+khbcd+1),1,ne,cr(khbb+1),nr,1,
     *cr(irbas+khbzd+1),1,ne,ne,nc,nr)
      if(odebug(40)) then
       call dbgciv('ddijkl',cr(icbas+khbzd+1),ne,nc,1,ne)
       if(ic.ne.ir)call dbgciv('ddijkl',cr(irbas+khbzd+1),ne,nr,1,ne)
      endif
      call getsb(rmode(1),1)
      if(model.eq.3)goto 999
      return
      end
      subroutine ssijkl(cr,zsc,ztc,csc,ctc,zsr,ztr,csr,ctr)
      implicit REAL  (a-h,o-z),integer (i-n)
INCLUDE(common/prnprn)
c...
c...  (ij/kl)    n-2/n-2     interactions
c...
      dimension zsc(*),ztc(*),csc(*),ctc(*),zsr(*),ztr(*),csr(*),ctr(*)
      dimension cr(*)
      common/spew/rmode(3)
      common/stoctl/khbill,khbz,khbc,khbb
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
c...      singlet / singlet
      if(lamss)2,1,2
2     call mxmb(csr,1,ns,cr(khbb+1),1,lamsr,zsc,1,ns,ns,lamsr,lamsc)
      if(ic.ne.ir)call mxmb(csc,1,ns,cr(khbb+1),lamsr,1,zsr,1,ns,
     *ns,lamsc,lamsr)
      if(odebug(40)) then
       call dbgciv('ssijkl',zsc,ns,lamsc,1,ns)
       if(ic.ne.ir)call dbgciv('ssijkl',zsr,ns,lamsr,1,ns)
      endif
c...      triplet / triplet
      if(lamtt)1,999,1
1      mu=lamss+khbb
      call mxmb(ctr,1,nt,cr(mu+1),1,lamtr,ztc,1,nt,nt,lamtr,lamtc)
      if(ic.ne.ir)call mxmb(ctc,1,nt,cr(mu+1),lamtr,1,ztr,1,nt,
     *nt,lamtc,lamtr)
      if(odebug(40)) then
       call dbgciv('ssijkl',ztc,nt,lamtc,1,nt)
       if(ic.ne.ir)call dbgciv('ssijkl',ztr,nt,lamtr,1,nt)
      endif
999   return
      end
      subroutine dvijka(cr,iconf)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
c...  doublet/vacuum (ij/ka) and fock(ia) contributions
      dimension cr(*),iconf(5,*)
      common/lsort/smode,scr(201)
      common/spew/rmode(5)
INCLUDE(common/diactl)
INCLUDE(common/symci)
INCLUDE(common/helpr)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh
INCLUDE(common/disktl)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir),
     *            (rmode(4),ninti),(rmode(5),icfock)
      equivalence (smode,integr)
       icl=0
      khbx=khbb+m2em1e
c... khbb    (ij/ka) integrals
c... khbx     coupling coefficients
c... read (ij/ka) integrals
      call rdedx(cr(khbb+1),m2e,index(2),numci)
c... initiate input of fock(ia) operators
      call brtsc(index(3))
      call getsc(rmode(5),1)
      loijka=.true.
55    call getsb(rmode(2),3)
      if(ic.eq.icl)goto 500
      icl=ic
_IF(cray,convex,t3e,i8,i8drct)
      itop=iconf(1,ic)
_ELSE
      itop=ipad(iconf(1,ic))
_ENDIF
      khzdub=khbzd+itop
      khcdub=khbcd+itop
_IF(cray,convex,t3e,i8,i8drct)
      ne=norbe(iconf(4,ic))
      lamd=iconf(3,ic)
_ELSE
      ne=norbe(ipad(iconf(4,ic)))
      lamd=ipad(iconf(3,ic))
_ENDIF
501    if(ic-icfock)500,1,2
c... skip redundant fock operator
_IF(cray,convex,t3e,i8,i8drct)
    2  call skipsc(norbe(iconf(4,icfock))+2)
_ELSE
    2  call skipsc(norbe(ipad(iconf(4,icfock)))+2)
_ENDIF
      goto 600
c... get fock operator
1     call getsc(smode,2)
      call getsc(cr(integr+kh+1),ne)
600    call getsc(rmode(5),1)
      goto 501
500   continue
_IF(cray,convex,t3e,i8,i8drct)
      lamv=iconf(3,ir)
      itop=iconf(1,ir)
_ELSE
      lamv=ipad(iconf(3,ir))
      itop=ipad(iconf(1,ir))
_ENDIF
      khi=ninti+kh
      call getsb(cr(khbx+1),lamv*lamd)
      if(modev)goto 502
c... construct l
      call mxmd(cr(khbx+1),lamv,1,cr(khbc+itop+1),1,1,scr(1),1,1,
     *lamd,lamv,1)
c... update z doublet
      call mxmb(cr(khi+1),1,1,scr(1),1,1,cr(khzdub+1),1,ne,
     *ne,1,lamd)
      if(odebug(40)) then
       call dbgciv('dvijka',cr(khzdub+1),ne,lamd,1,ne)
      endif
502   if(moded)goto 503
c... construct m
      call mxmd(cr(khcdub+1),ne,1,cr(khi+1),1,1,scr(1),1,1,
     *lamd,ne,1)
c... update z vacuum
      call mxmb(cr(khbx+1),1,lamv,scr(1),1,1,cr(khbz+itop+1),1,1,
     *lamv,lamd,1)
      if(odebug(40)) then
       call dbgciv('dvijka',cr(khbz+itop+1),lamv,1,1,1)
      endif
503   call getsb(rmode(1),1)
      if(model.eq.5)goto 55
      return
      end
      subroutine sdijka(cr,iconf)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      logical iscne1,isceq1
c...  n-2 / doublet (ij/ka) and fock(ia) contributions
      dimension cr(*),iconf(5,*)
INCLUDE(common/table)
      common/lsort/smode(2),scr(nd200),neii(8),
     *mcolcz(8),mrowcz(8),ibassc(8),ibassz(8),ibattc(8),ibattz(8)
      common/spew/rmode(5)
INCLUDE(common/diactl)
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt,
     *khcsc,khzsc,khctc,khztc,isc
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88,khb111,khb222,khb333
INCLUDE(common/symci)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir),
     * (rmode(4),ninti),(rmode(5),icfock)
c... khb111  coupling coefficients
      equivalence (smode(1),integr),(smode(2),nwww)
c... khb222  squared c or z matrix sub block
c... khb333  m or l
      icl=ic
      isceq1=isc.eq.1
      iscne1=.not.isceq1
      do 4444 isr=1,nirr
      isi=mult(isr,isc)
      nei=norbe(isi)
      neii(isr)=nei
      if(isi.ge.isr)goto 4455
      iss=isi
      mcolcz(isr)=nei
      mrowcz(isr)=1
      goto 4466
4455  iss=isr
      mcolcz(isr)=1
      mrowcz(isr)=norbe(isr)
4466  ibassc(isr)=ibassi(iss,isc)+khcsc
      ibassz(isr)=ibassi(iss,isc)+khzsc
      ibattc(isr)=ibastr(iss,isc)+khctc
4444  ibattz(isr)=ibastr(iss,isc)+khztc
6501   if(icfock-ic)62,63,6500
c... skip redundant fock operator
62    call getsc(smode(1),2)
       call skipsc(nwww)
       goto 6600
c... get fock operator
63     call getsc(smode(1),2)
       call getsc(cr(integr+kh+1),nwww)
6600   call getsc(rmode(5),1)
       goto 6501
6500   continue
_IF(cray,convex,t3e,i8,i8drct)
       lamr=iconf(3,ir)
       isr=iconf(4,ir)
       itop=iconf(1,ir)
_ELSE
       lamr=ipad(iconf(3,ir))
       isr=ipad(iconf(4,ir))
       itop=ipad(iconf(1,ir))
_ENDIF
      ner=norbe(isr)
       lamt=lamtc
      if(isceq1.and.ner.eq.1)lamt=0
      lamlam=lamsc+lamt
       call getsb(cr(khb111+1),lamlam*lamr)
      nei=neii(isr)
      intbas=ninti+kh
      mcol=mcolcz(isr)
      mrow=mrowcz(isr)
       if(moded)goto 1112
       if(lamr.ne.1)goto 2000
c... use partial matrix element driven scheme
       irl=ir
       call mxmd(cr(intbas+1),1,1,cr(khb111+1),1,1,
     *cr(khb333+1),1,nei,nei,1,lamlam)
3001  call getsb(rmode(1),1)
      if(model.ne.6)goto 3000
      call getsb(rmode(2),3)
      if(ic.ne.icl.or.ir.ne.irl)goto 3000
      call getsb(cr(khb111+1),lamlam)
      call mxmb(cr(ninti+kh+1),1,1,cr(khb111+1),1,1,
     *cr(khb333+1),1,nei,nei,1,lamlam)
      goto 3001
c... update z(n-2)
3000   icbas=ibassz(isr)
       icbat=ibattz(isr)
      imbas=khb333
      khbdub=khbcd+itop
      if(iscne1)goto 300
c... z(n-2) totally symmetric
      if(lamsc)303,304,303
303   do 3003 lam=1,lamsc
      call mxmd(cr(khbdub+1),1,1,cr(imbas+1),1,1,
     *cr(khb222+1),1,ner,ner,1,nei)
      imbas=imbas+nei
      call symm1a(cr(icbas+1),cr(khb222+1),ner)
3003  icbas=icbas+ns
      if(lamt)304,3111,304
304   do 3004 lam=1,lamt
      call mxmd(cr(khbdub+1),1,1,cr(imbas+1),1,1,
     *cr(khb222+1),1,ner,ner,1,nei)
      imbas=imbas+nei
      call anti1a(cr(icbat+1),cr(khb222+1),ner)
3004  icbat=icbat+nt
      goto 3111
c... z(n-2) not totally symmetric
300    if(lamsc)305,306,305
305     do 3005 lam=1,lamsc
        if(ner.ge.nei)goto 3555
       call mxmb(cr(imbas+1),1,1,cr(khbdub+1),1,1,
     *cr(icbas+1),mrow,mcol,nei,1,ner)
       goto 3444
3555  call mxmb(cr(khbdub+1),1,1,cr(imbas+1),1,1,
     *cr(icbas+1),mcol,mrow,ner,1,nei)
3444  imbas=imbas+nei
3005   icbas=icbas+ns
      if(lamt)306,3111,306
306   do 3006 lam=1,lamt
      if(ner.ge.nei)goto 3777
      call mxmb(cr(imbas+1),1,1,cr(khbdub+1),1,1,
     *cr(icbat+1),mrow,mcol,nei,1,ner)
      goto 3666
3777  call mxmb(cr(khbdub+1),1,1,cr(imbas+1),1,1,
     *cr(icbat+1),mcol,mrow,ner,1,nei)
3666  imbas=imbas+nei
3006  icbat=icbat+nt
3111  if(modest)goto 3222
c... update z(doublet)
       icbas=ibassc(isr)
       icbat=ibattc(isr)
       imbas=khb333
       khbdub=khbzd+itop
      if(iscne1)goto 400
c... c(n-2) totally symmetric
      if(lamsc)403,404,403
403   do 4003 lam=1,lamsc
      call square(cr(khb222+1),cr(icbas+1),ner,ner)
      icbas=icbas+ns
      call mxmb(cr(khb222+1),1,ner,cr(imbas+1),1,1,
     *cr(khbdub+1),1,1,ner,nei,1)
4003  imbas=imbas+nei
       if(lamt)404,3222,404
404   do 4004 lam=1,lamt
      call sqtrip(cr(khb222+1),cr(icbat+1),ner)
      icbat=icbat+nt
      call mxmb(cr(khb222+1),1,ner,cr(imbas+1),1,1,
     *cr(khbdub+1),1,1,ner,nei,1)
4004  imbas=imbas+nei
      goto 3222
c... c(n-2) not totally symmetric
400   if(lamsc)405,406,405
405   do 4005 lam=1,lamsc
      call mxmb(cr(icbas+1),mcol,mrow,cr(imbas+1),1,1,
     *cr(khbdub+1),1,1,ner,nei,1)
      icbas=icbas+ns
4005  imbas=imbas+nei
      if(lamt)406,3222,406
406   do 4006 lam=1,lamt
      call mxmb(cr(icbat+1),mcol,mrow,cr(imbas+1),1,1,
     *cr(khbdub+1),1,1,ner,nei,1)
      icbat=icbat+nt
4006  imbas=imbas+nei
3222   if(model-6)999,2001,999
c... update z(n-2)
2000  call mxmd(cr(khbcd+itop+1),1,ner,cr(khb111+1),1,lamr,
     *cr(khb333+1),1,ner,ner,lamr,lamlam)
      icbas=ibassz(isr)
      icbat=ibattz(isr)
      imbas=khb333
      if(iscne1)goto 900
c... z(n-2) totally symmetric
      if(lamsc)903,904,903
903   do 9003 lam=1,lamsc
      call mxmd(cr(imbas+1),1,1,cr(intbas+1),1,1,
     *cr(khb222+1),1,ner,ner,1,nei)
      imbas=imbas+ner
      call symm1a(cr(icbas+1),cr(khb222+1),ner)
9003  icbas=icbas+ns
      if(lamt)904,1111,904
904   do 9004 lam=1,lamt
      call mxmd(cr(imbas+1),1,1,cr(intbas+1),1,1,
     *cr(khb222+1),1,ner,ner,1,nei)
      imbas=imbas+ner
      call anti1a(cr(icbat+1),cr(khb222+1),ner)
9004  icbat=icbat+nt
      goto 1111
c... z(n-2) not totally symmetric
900   if(lamsc)905,906,905
905   do 9005 lam=1,lamsc
      if(ner.ge.nei)goto 9555
      call mxmb(cr(intbas+1),1,1,cr(imbas+1),1,1,cr(icbas+1),mrow,mcol,
     *nei,1,ner)
      if(odebug(40)) then
       call dbgciv('sdijka',cr(icbas+1),nei,ner,mrow,mcol)
      endif
      goto 9444
9555  call mxmb(cr(imbas+1),1,1,cr(intbas+1),1,1,cr(icbas+1),mcol,mrow,
     *ner,1,nei)
      if(odebug(40)) then
       call dbgciv('sdijka',cr(icbas+1),ner,nei,mcol,mrow)
      endif
9444  imbas=imbas+ner
9005  icbas=icbas+ns
      if(lamt)906,1111,906
906   do 9006 lam=1,lamt
      if(ner.ge.nei)goto 9777
      call mxmb(cr(intbas+1),1,1,cr(imbas+1),1,1,cr(icbat+1),mrow,mcol,
     *nei,1,ner)
      if(odebug(40)) then
       call dbgciv('sdijka',cr(icbat+1),nei,ner,mrow,mcol)
      endif
      goto 9666
9777  call mxmb(cr(imbas+1),1,1,cr(intbas+1),1,1,cr(icbat+1),mcol,mrow,
     *ner,1,nei)
      if(odebug(40)) then
       call dbgciv('sdijka',cr(icbat+1),ner,nei,mcol,mrow)
      endif
9666  imbas=imbas+ner
9006  icbat=icbat+nt
1111   if(modest)goto 2222
c... update z(doublet)
1112  icbas=ibassc(isr)
      icbat=ibattc(isr)
      imbas=khb333
      if(iscne1)goto 700
c... c(n-2) totally symmetric
      if(lamsc)703,704,703
703   do 7003 lam=1,lamsc
      call square(cr(khb222+1),cr(icbas+1),ner,ner)
      icbas=icbas+ns
      call mxmd(cr(khb222+1),1,ner,cr(intbas+1),1,1,
     *cr(imbas+1),1,1,ner,nei,1)
7003  imbas=imbas+ner
      if(lamt)704,750,704
704   do 7004 lam=1,lamt
      call sqtrip(cr(khb222+1),cr(icbat+1),ner)
      icbat=icbat+nt
      call mxmd(cr(khb222+1),1,ner,cr(intbas+1),1,1,
     *cr(imbas+1),1,1,ner,nei,1)
7004  imbas=imbas+ner
      goto 750
c...  c(n-2) not totally symmetric
700   if(lamsc)705,706,705
705   do 7005 lam=1,lamsc
      call mxmd(cr(icbas+1),mcol,mrow,cr(intbas+1),1,1,
     *cr(imbas+1),1,1,ner,nei,1)
      icbas=icbas+ns
7005  imbas=imbas+ner
      if(lamt)706,750,706
706   do 7006 lam=1,lamt
      call mxmd(cr(icbat+1),mcol,mrow,cr(intbas+1),1,1,
     *cr(imbas+1),1,1,ner,nei,1)
      icbat=icbat+nt
7006  imbas=imbas+ner
750   call mxmb(cr(khb333+1),1,ner,cr(khb111+1),lamr,1,
     *cr(khbzd+itop+1),1,ner,ner,lamlam,lamr)
      if(odebug(40)) then
       call dbgciv('sdijka',cr(khbzd+itop+1),ner,lamr,1,ner)
      endif
2222   call getsb(rmode(1),1)
       if(model.ne.6)goto 999
       call getsb(rmode(2),3)
2001   if(ic.eq.icl)goto 6500
999    return
       end
      subroutine sdiabc(cr,iconf)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      logical ltgtld
c...  n-2 / doublet (ia/bc) contributions
      dimension cr(*),iconf(5,*)
      common/spew/rmode(5)
       common/lsort/jbass3(8)
INCLUDE(common/diactl)
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt,
     *khcsc,khzsc,khctc,khztc
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88,khb111,khb222,khb333
     *,khbiii,khbkkk,khbddd
INCLUDE(common/sizes)
INCLUDE(common/symci)
INCLUDE(common/helpr)
_IF(ga)
      common/rem_iabc/myext(62,nd200),myibass3(8,nd200),myniky(nd200,8),
     *mynorbe(nd200,8),iabase(62),opcivect
      logical opcivect
_ENDIF
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
     *           ,(rmode(4),ninti),(rmode(5),icfock)
c... khbiii s-d  coupling coefficients
c... khbjjj t-d coupling coefficients
c... khbkkk   l or m matrices
_IF(ga)
      nproc = ipg_nnodes()
      if (opcivect)nproc=1
_ENDIF
      nsnt=ns+nt
      lamst=lamsc+lamtc
      icl=ic
333    continue
_IF(cray,convex,t3e,i8,i8drct)
      lamr=iconf(3,ir)
      itop=iconf(1,ir)
      isr=iconf(4,ir)
_ELSE
      lamr=ipad(iconf(3,ir))
      itop=ipad(iconf(1,ir))
      isr=ipad(iconf(4,ir))
_ENDIF
_IF(ga)
      ner=mynorbe(ninti,isr)
      ner2=norbe(isr)
      istq = ipg_stdub(ner,ner2,isr,ninti)
_ELSE
      ner=norbe(isr)
_ENDIF
c
_IF(ga)
      ibasp=myibass3(isr,ninti)+iabase(ninti)
_ELSE
      ibasp=jbass3(isr)
_ENDIF
      ibasq=ibasp+ns
      khzdub=khbzd+itop
      khcdub=khbcd+itop
      khbjjj=khbiii+lamr*lamsc
       call getsb(cr(khbiii+1),lamst*lamr)
      ltgtld=lamtc.gt.lamr
       if(modest)goto 1112
_IF(ga)
c... update z(doublet)
       if(lamsc)22,21,22
22    call mxmd(cr(ibasp+1),nsnt,1,cr(khcsc+1),1,ns,
     *cr(khbkkk+1),1,ner,ner,ns,lamsc)
      call mxmb(cr(khbkkk+1),1,ner,cr(khbiii+1),lamr,1,
     *cr(khzdub+istq),nproc,ner2,ner,lamsc,lamr)
21    if(ltgtld)goto 31
      if(lamtc)32,1111,32
32    call mxmd(cr(ibasq+1),nsnt,1,cr(khctc+1),1,nt,
     *cr(khbkkk+1),1,ner,ner,nt,lamtc)
      call mxmb(cr(khbkkk+1),1,ner,cr(khbjjj+1),lamr,1,
     *cr(khzdub+istq),nproc,ner2,ner,lamtc,lamr)
      goto 1111
31    call mxmd(cr(khctc+1),1,nt,cr(khbjjj+1),lamr,1,
     *cr(khbkkk+1),1,nt,nt,lamtc,lamr)
      call mxmb(cr(ibasq+1),nsnt,1,cr(khbkkk+1),1,nt,
     *cr(khzdub+istq),nproc,ner2,ner,nt,lamr)
1111   if(moded)goto 999
c... update z(singlet)
1112   if(lamsc)42,41,42
42    continue
      call mxmd(cr(khcdub+istq),nproc,ner2,cr(khbiii+1),1,lamr,
     *cr(khbkkk+1),1,ner,ner,lamr,lamsc)
      call mxmb(cr(ibasp+1),1,nsnt,cr(khbkkk+1),1,ner,
     *cr(khzsc+1),1,ns,ns,ner,lamsc)
c... update z(triplet)
41    if(ltgtld)goto 50
      if(lamtc)52,999,52
52    call mxmd(cr(khcdub+istq),nproc,ner2,cr(khbjjj+1),1,lamr,
     *cr(khbkkk+1),1,ner,ner,lamr,lamtc)
      call mxmb(cr(ibasq+1),1,nsnt,cr(khbkkk+1),1,ner,
     *cr(khztc+1),1,nt,nt,ner,lamtc)
      goto 999
50    call mxmd(cr(ibasq+1),1,nsnt,cr(khcdub+istq),nproc,ner2,
     *cr(khbkkk+1),1,nt,nt,ner,lamr)
       call mxmb(cr(khbkkk+1),1,nt,cr(khbjjj+1),1,lamr,
     *cr(khztc+1),1,nt,nt,lamr,lamtc)
999   call getsb(rmode(1),1)
_ELSE
c... update z(doublet)
       if(lamsc)22,21,22
22    call mxmd(cr(ibasp+1),nsnt,1,cr(khcsc+1),1,ns,
     *cr(khbkkk+1),1,ner,ner,ns,lamsc)
      call mxmb(cr(khbkkk+1),1,ner,cr(khbiii+1),lamr,1,
     *cr(khzdub+1),1,ner,ner,lamsc,lamr)
21    if(ltgtld)goto 31
      if(lamtc)32,1111,32
32    call mxmd(cr(ibasq+1),nsnt,1,cr(khctc+1),1,nt,
     *cr(khbkkk+1),1,ner,ner,nt,lamtc)
      call mxmb(cr(khbkkk+1),1,ner,cr(khbjjj+1),lamr,1,
     *cr(khzdub+1),1,ner,ner,lamtc,lamr)
      goto 1111
31     call mxmd(cr(khctc+1),1,nt,cr(khbjjj+1),lamr,1,
     *cr(khbkkk+1),1,nt,nt,lamtc,lamr)
      call mxmb(cr(ibasq+1),nsnt,1,cr(khbkkk+1),1,nt,
     *cr(khzdub+1),1,ner,ner,nt,lamr)
1111   if(moded)goto 999
c... update z(singlet)
1112   if(lamsc)42,41,42
42    continue
      call mxmd(cr(khcdub+1),1,ner,cr(khbiii+1),1,lamr,
     *cr(khbkkk+1),1,ner,ner,lamr,lamsc)
      call mxmb(cr(ibasp+1),1,nsnt,cr(khbkkk+1),1,ner,
     *cr(khzsc+1),1,ns,ns,ner,lamsc)
c... update z(triplet)
41    if(ltgtld)goto 50
      if(lamtc)52,999,52
52    call mxmd(cr(khcdub+1),1,ner,cr(khbjjj+1),1,lamr,
     *cr(khbkkk+1),1,ner,ner,lamr,lamtc)
      call mxmb(cr(ibasq+1),1,nsnt,cr(khbkkk+1),1,ner,
     *cr(khztc+1),1,nt,nt,ner,lamtc)
      goto 999
50    call mxmd(cr(ibasq+1),1,nsnt,cr(khcdub+1),1,ner,
     *cr(khbkkk+1),1,nt,nt,ner,lamr)
       call mxmb(cr(khbkkk+1),1,nt,cr(khbjjj+1),1,lamr,
     *cr(khztc+1),1,nt,nt,lamr,lamtc)
999   call getsb(rmode(1),1)
_ENDIF
      if(model.ne.13)goto 222
      call getsb(rmode(2),3)
_IF(ga)
      if(ic.eq.icl)goto 333
_ELSE
      if(ic.eq.icl.and.ninti.eq.icfock)goto 333
_ENDIF
222   return
      end
      subroutine ddijab(cr,iconf)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
c... doublet / doublet (ij/ab) + (ia/ib) + fock(ab) interactions
c... all these operators are hermitian
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      dimension cr(*),iconf(5,*)
      common/expanc/joker(8)
INCLUDE(common/symci)
INCLUDE(common/auxh)
INCLUDE(common/diactl)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh
INCLUDE(common/helpr)
INCLUDE(common/disktl)
      common/spew/rmode(4)
      equivalence(rmode(1),model),(rmode(2),ic),(rmode(3),ir),
     *           (rmode(4),kkkbas)
      icl=0
c... khbb   (ia/ib) + (ij/ab) integrals
c... khbw   squared operator
c... khbx   triangle of fock operator or coupling coefficients
c... khby   l or m matrices
       khbw=khbb+m2exab
       khbx=khbw+mubsq
       khby=khbx+mxddpa
c... read (ia/ib) + (ij/ab) integrals
      call rdedx(cr(khbb+1),m2fock,index(4),numci)
      loijab=.true.
      call reads(cr(khbb+m2fock+1),m2coul,numci)
c... initiate input of fock(ab) operators
      call brtsc(index(7))
55    call getsb(rmode(2),3)
      if(icl.eq.ic)goto 500
      icl=ic
_IF(cray,convex,t3e,i8,i8drct)
      itop=iconf(1,ic)
_ELSE
      itop=ipad(iconf(1,ic))
_ENDIF
      khzduc=khbzd+itop
      khcduc=khbcd+itop
_IF(cray,convex,t3e,i8,i8drct)
      lamc=iconf(3,ic)
      isc=iconf(4,ic)
_ELSE
      lamc=ipad(iconf(3,ic))
      isc=ipad(iconf(4,ic))
_ENDIF
      nec=norbe(isc)
      do 8 loop=1,isc
    8  joker(loop)=madsi(loop,isc)+kh
       joke=joker(isc)
500   if(kkkbas)900,899,900
c... fock operator contribution
899   call getsc(cr(khbx+1),iky(nec+1))
      call square(cr(khbw+1),cr(khbx+1),nec,nec)
      call mxmb(cr(khbw+1),1,nec,cr(khcduc+1),1,nec,cr(khzduc+1),1,nec,
     *nec,nec,lamc)
      if(odebug(40)) then
       call dbgciv('ddijab',cr(khzduc+1),nec,lamc,1,nec)
      endif
      goto 503
900   continue
_IF(cray,convex,t3e,i8,i8drct)
      lamr=iconf(3,ir)
_ELSE
      lamr=ipad(iconf(3,ir))
_ENDIF
       call getsb(cr(khbx+1),lamr*lamc)
_IF(cray,convex,t3e,i8,i8drct)
      isr=iconf(4,ir)
_ELSE
      isr=ipad(iconf(4,ir))
_ENDIF
       ner=norbe(isr)
_IF(cray,convex,t3e,i8,i8drct)
      itop=iconf(1,ir)
_ELSE
      itop=ipad(iconf(1,ir))
_ENDIF
      if(isc.ne.isr)goto 665
c... operator is totally symmetric
      call square(cr(khbw+1),cr(kkkbas+joke+1),nec,nec)
      kbas=khbw
      mmm2=1
      nnn2=nec
      goto 999
c... operator not totally symmetric
665   kbas=kkkbas+joker(isr)
       mmm2=ner
      nnn2=1
999   if(lamr-lamc)1999,1,2999
1     if(nec.le.ner)goto 1999
c...  m(ic)
2999   call mxmd(cr(kbas+1),1,ner,cr(khcduc+1),1,nec,
     *cr(khby+1),1,ner,ner,nec,lamc)
c... update z(ir)
      call mxmb(cr(khby+1),1,ner,cr(khbx+1),lamr,1,
     *cr(khbzd+itop+1),1,ner,ner,lamc,lamr)
      if(odebug(40)) then
       call dbgciv('ddijab',cr(khbzd+itop+1),ner,lamr,1,ner)
      endif
      if(ic.eq.ir)goto 503
      call mxmd(cr(khbcd+itop+1),1,ner,cr(khbx+1),1,lamr,
     *cr(khby+1),1,ner,ner,lamr,lamc)
c... update z(ic)
      call mxmb(cr(kbas+1),mmm2,nnn2,cr(khby+1),1,ner,
     *cr(khzduc+1),1,nec,nec,ner,lamc)
      if(odebug(40)) then
       call dbgciv('ddijab',cr(khzduc+1),nec,lamc,1,nec)
      endif
      goto 503
1999  call mxmd(cr(khcduc+1),1,nec,cr(khbx+1),lamr,1,cr(khby+1),
     *1,nec,nec,lamc,lamr)
c... update z(ir)
      call mxmb(cr(kbas+1),1,ner,cr(khby+1),1,nec,cr(khbzd+itop+1),
     *1,ner,ner,nec,lamr)
      if(odebug(40)) then
       call dbgciv('ddijab',cr(khbzd+itop+1),ner,lamr,1,ner)
      endif
      if(ic.eq.ir)goto 503
      call mxmd(cr(kbas+1),mmm2,nnn2,cr(khbcd+itop+1),1,ner,
     *cr(khby+1),1,nec,nec,ner,lamr)
c...   update  z(ic)
      call mxmb(cr(khby+1),1,nec,cr(khbx+1),1,lamr,
     *cr(khzduc+1),1,nec,nec,lamr,lamc)
      if(odebug(40)) then
       call dbgciv('ddijab',cr(khzduc+1),nec,lamc,1,nec)
      endif
503   call getsb(rmode(1),1)
      if(model.eq.7)goto 55
      return
      end
      subroutine ddiajb(cr,iconf)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
c...  doublet /doublet  (ia/jb)  interactions
c... these exchange operators not hermitian
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      dimension cr(*),iconf(5,*)
      common/expanc/joker(8),koker(8)
INCLUDE(common/symci)
INCLUDE(common/auxh)
INCLUDE(common/diactl)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh
INCLUDE(common/helpr)
INCLUDE(common/disktl)
      common/spew/rmode(5)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir),
     *            (rmode(4),kkkbas),(rmode(5),itrans)
      icl=0
c... khbb (ia/jb) integrals
c... khbx   coupling coefficients
c... khby   l or m matrices
      khbx=khbb+m2exch
      if(moded)goto 55
       khby=khbx+mxddpa
c... read (ia/jb) integrals
      call rdedx(cr(khbb+1),m2exch,index(8),numci)
      loiajb=.true.
55    call getsb(rmode(2),4)
      if(icl.eq.ic)goto 500
      icl=ic
_IF(cray,convex,t3e,i8,i8drct)
      itop=iconf(1,ic)
_ELSE
      itop=ipad(iconf(1,ic))
_ENDIF
      khzduc=khbzd+itop
      khcduc=khbcd+itop
_IF(cray,convex,t3e,i8,i8drct)
      lamc=iconf(3,ic)
      isc=iconf(4,ic)
_ELSE
      lamc=ipad(iconf(3,ic))
      isc=ipad(iconf(4,ic))
_ENDIF
      nec=norbe(isc)
      do 8 loop=1,isc
      koker(loop)=nadsq(loop,isc)+kh
    8 joker(loop)=madsq(loop,isc)+kh
500   continue
_IF(cray,convex,t3e,i8,i8drct)
      lamr=iconf(3,ir)
_ELSE
      lamr=ipad(iconf(3,ir))
_ENDIF
       call getsb(cr(khbx+1),lamr*lamc)
       if(moded)goto 503
_IF(cray,convex,t3e,i8,i8drct)
      isr=iconf(4,ir)
      itop=iconf(1,ir)
_ELSE
      isr=ipad(iconf(4,ir))
      itop=ipad(iconf(1,ir))
_ENDIF
       ner=norbe(isr)
       if(itrans)665,664,664
c... matrix transposed
664    kbas=kkkbas+joker(isr)
       mmm1=ner
       nnn1=1
       goto 999
c... matrix not transposed
665    kbas=kkkbas+koker(isr)
      nnn1=nec
      mmm1=1
999   if(lamr-lamc)1999,1,2999
1     if(nec.le.ner)goto 1999
c...  n(ic)
2999   call mxmd(cr(kbas+1),nnn1,mmm1,cr(khcduc+1),1,nec,
     *cr(khby+1),1,ner,ner,nec,lamc)
c... update z(ir)
      call mxmb(cr(khby+1),1,ner,cr(khbx+1),lamr,1,
     *cr(khbzd+itop+1),1,ner,ner,lamc,lamr)
      if(odebug(40)) then
       call dbgciv('ddiajb',cr(khbzd+itop+1),ner,lamr,1,ner)
      endif
      call mxmd(cr(khbcd+itop+1),1,ner,cr(khbx+1),1,lamr,
     *cr(khby+1),1,ner,ner,lamr,lamc)
c... update z(ic)
      call mxmb(cr(kbas+1),mmm1,nnn1,cr(khby+1),1,ner,
     *cr(khzduc+1),1,nec,nec,ner,lamc)
      if(odebug(40)) then
       call dbgciv('ddiajb',cr(khzduc+1),nec,lamc,1,nec)
      endif
      goto 503
1999  call mxmd(cr(khcduc+1),1,nec,cr(khbx+1),lamr,1,cr(khby+1),
     *1,nec,nec,lamc,lamr)
c... update z(ir)
      call mxmb(cr(kbas+1),nnn1,mmm1,cr(khby+1),1,nec,
     *cr(khbzd+itop+1),1,ner,ner,nec,lamr)
      if(odebug(40)) then
       call dbgciv('ddiajb',cr(khbzd+itop+1),ner,lamr,1,ner)
      endif
      call mxmd(cr(kbas+1),mmm1,nnn1,cr(khbcd+itop+1),1,ner,
     *cr(khby+1),1,nec,nec,ner,lamr)
c...   update  z(ic)
      call mxmb(cr(khby+1),1,nec,cr(khbx+1),1,lamr,
     *cr(khzduc+1),1,nec,nec,lamr,lamc)
      if(odebug(40)) then
       call dbgciv('ddiajb',cr(khzduc+1),nec,lamc,1,nec)
      endif
503   call getsb(rmode(1),1)
      if(model.eq.9)goto 55
      return
      end
      subroutine sviaib(cr,iconf,zs,cs,b)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
INCLUDE(common/prnprn)
      dimension b(*),zs(*),cs(*)
      dimension cr(*),iconf(5,*)
c...   singlet / vacuum (ia/ib) interactions
INCLUDE(common/symci)
      common/moco/ns,nt,lams,lamt
      common/lsort/scr(1)
INCLUDE(common/diactl)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh
      common/spew/rmode(4)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir),
     *            (rmode(4),ninti)
      icl=ic
_IF(cray,convex,t3e,i8,i8drct)
    5 lamv=iconf(3,ir)
      itop=iconf(1,ir)
_ELSE
    5 lamv=ipad(iconf(3,ir))
      itop=ipad(iconf(1,ir))
_ENDIF
      call getsb(b(1),lamv*lams)
      ipbas=ninti+kh
      if(modest)goto 70
c... k(ic)
      call mxmd(cs,norbsi(1),1,cr(ipbas+1),1,1,scr,1,1,
     +          lams,norbsi(1),1)
c... update z vacuum
      call mxmb(b,1,lamv,scr,1,1,cr(khbz+itop+1),1,1,
     *lamv,lams,1)
      if(odebug(40)) then
       call dbgciv('sviaib',cr(khbz+itop+1),lamv,1,1,1)
      endif
70    if(modev)goto 80
c... l(ic,ir)
      call mxmd(b,lamv,1,cr(khbc+itop+1),1,1,scr,1,1,
     *lams,lamv,1)
c... update z singlet
      call mxmb(cr(ipbas+1),1,1,scr,1,1,zs,1,norbsi(1),
     +          norbsi(1),1,lams)
      if(odebug(40)) then
       call dbgciv('sviaib',zs,norbsi(1),lams,1,norbsi(1))
      endif
80    call getsb(rmode(1),1)
      if(model.ne.12)goto 999
      call getsb(rmode(2),3)
      if(ic.eq.icl)goto 5
999   return
      end
      subroutine sviajb(cr,iconf,zs,zt,cs,ct,b)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
INCLUDE(common/prnprn)
      dimension b(*),zs(*),zt(*),cs(*),ct(*)
      dimension cr(*),iconf(5,*)
c...   n-2 / vacuum (ia/jb) interactions
      common/moco/ns,nt,lams,lamt
      common/lsort/scr(1)
INCLUDE(common/diactl)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh
      common/spew/rmode(4)
      equivalence (rmode(4),ninti)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      icl=ic
      lamst=lams+lamt
_IF(cray,convex,t3e,i8,i8drct)
    5 lamv=iconf(3,ir)
      itop=iconf(1,ir)
_ELSE
    5 lamv=ipad(iconf(3,ir))
      itop=ipad(iconf(1,ir))
_ENDIF
      call getsb(b(1),lamv*lamst)
      ipbas=ninti+kh
      iqbas=ipbas+ns
      if(modest)goto 70
      if(lams)98,96,98
c...  k(ic)
98    call mxmd(cs,ns,1,cr(ipbas+1),1,1,scr,1,1,lams,ns,1)
      if(lamt)96,71,96
96    call mxmd(ct,nt,1,cr(iqbas+1),1,1,scr(lams+1),1,1,lamt,nt,1)
c... update z vacuum
71    call mxmb(b,1,lamv,scr,1,1,cr(khbz+itop+1),1,1,
     *lamv,lamst,1)
      if(odebug(40)) then
       call dbgciv('sviajb',cr(khbz+itop+1),lamv,1,1,1)
      endif
70    if(modev)goto 80
c...  l(ic,ir)
      call mxmd(b,lamv,1,cr(khbc+itop+1),1,1,scr,1,1,
     *lamst,lamv,1)
      if(lams)97,99,97
c... update z singlet
97    call mxmb(cr(ipbas+1),1,1,scr,1,1,zs,1,ns,ns,1,lams)
      if(odebug(40)) then
       call dbgciv('sviajb',zs,ns,lams,1,ns)
      endif
      if(lamt)99,80,99
c... update z triplet
99    call mxmb(cr(iqbas+1),1,1,scr(lams+1),1,1,zt,1,nt,nt,1,lamt)
      if(odebug(40)) then
       call dbgciv('sviajb',zt,nt,lamt,1,nt)
      endif
80    call getsb(rmode(1),1)
      if(model.ne.11)goto 999
      call getsb(rmode(2),3)
      if(ic.eq.icl)goto 5
999   return
      end
      subroutine renorm(c,fac,nspips)
      implicit REAL  (a-h,o-z),integer (i-n)
      dimension c(*)
INCLUDE(common/table)
INCLUDE(common/symci)
      if(nspips)1,999,1
1      mu=0
      do 2 irr=1,nirr
      ne=norbe(irr)
      if(ne)3,2,3
3     do 4 ie=1,ne
      mu=mu+ie
      call dscal(nspips,fac,c(mu),norbsi(1))
 4    continue
2     continue
999   return
      end
      subroutine mulss(oper,c,lamc,norsic,isc,
     *z,lamz,norz,isz,b,mcolb,mrowb,r)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      dimension b(*),oper(*),c(*),z(*),r(*)
INCLUDE(common/table)
INCLUDE(common/helpr)
      common/expanc/ibasso(8),mcolo(8),mrowo(8)
INCLUDE(common/symci)
c
c... z(singlet)  -  c(singlet)  (ij/ab)  fock(ab) interactions
      ibbb=1
      izzz=1
      if(lamz.le.lamc)goto 444
      do 1100 lam=1,lamc
      do 1200 iooz=1,nirr
      niz=norbe(iooz)
      if(niz)1201,1200,1201
1201  jooz=mult(iooz,isz)
      njz=norbe(jooz)
      if(njz)1202,1200,1202
1202  jooo=mult(jooz,isc)
      njo=norbe(jooo)
      if(njo)1203,1200,1203
1203  lbas=niz*njz
      jooioo=jooz-iooz
      if(jooioo)1301,1302,1302
c...   z  vector  transposed
1301   nnnz=1
       mmmz=njz
       jump=ibassi(jooz,isz)
      goto 1888
c... z vector not transposed
1302   mmmz=1
       nnnz=niz
       jump=ibassi(iooz,isz)
1888   if(jooz-jooo)1211,1210,1212
c... c vector totally symmetric
1210   call square(r(lbas+1),c(izzz+ibassi(jooo,1)),njo,njo)
       if(niz.ge.njz)goto 1400
       call mxmd(r(lbas+1),1,njo,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),r,nnnz,mmmz,njz,njo,niz)
      goto 1401
1400  call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(lbas+1),
     *1,njo,r,mmmz,nnnz,niz,njo,njz)
      goto 1401
c... c vector transposed
1211  nnnc=1
      mmmc=njz
      jum=ibassi(jooz,isc)+izzz
      goto 1213
c... c vector not transposed
1212  mmmc=1
      nnnc=njo
      jum=ibassi(jooo,isc)+izzz
1213  if(niz.ge.njz)goto 1214
      call mxmd(c(jum),nnnc,mmmc,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),r,nnnz,mmmz,njz,njo,niz)
      goto 1401
1214  call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),c(jum),
     *mmmc,nnnc,r,mmmz,nnnz,niz,njo,njz)
1401  kbas=0
      if(jooioo)1403,1402,1403
1402  kbas=lbas
      lbas=iky(niz+1)
      call symm1(r(kbas+1),r,niz)
1403  call mxmb(r(kbas+1),1,0,b(ibbb),1,mrowb,z(jump+1),1,norz,
     *lbas,1,lamz)
1200  continue
      izzz=izzz+norsic
1100  ibbb=ibbb+mcolb
      goto 999
444   do 100 lam=1,lamz
      if(mrowb)112,111,112
c... fock operator contribution
 111  continue
      call dcopy(norsic,c(izzz),1,r,1)
      goto 113
112   call mxmd(c,1,norsic,b(ibbb),mcolb,mrowb,r,1,1,
     *norsic,lamc,1)
113   do 200 iooz=1,nirr
      niz=norbe(iooz)
      if(niz)201,200,201
201   jooz=mult(iooz,isz)
      njz=norbe(jooz)
      if(njz)202,200,202
202   jooo=mult(jooz,isc)
c...   z(iooz,jooz) = oper(iooz,jooo) * c(jooo,jooz)
      njo=norbe(jooo)
      if(njo)203,200,203
203   lbas=norsic
      if(jooz-jooo)211,210,212
c... c vector totally symmetric
210   call square(r(lbas+1),r(ibassi(jooo,1)+1),njo,njo)
      kbas=lbas
      lbas=lbas+njo*njo
      goto 213
c... c transposed
211   kbas=ibassi(jooz,isc)
      mmmc=njz
      nnnc=1
      goto 214
c... c not transposed
212   kbas=ibassi(jooo,isc)
213   nnnc=njo
      mmmc=1
214   if(jooz-iooz)301,300,302
c...  z vector totally symmetric
300   call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(kbas+1),
     *mmmc,nnnc,r(lbas+1),1,niz,niz,njo,njz)
      if(odebug(40)) then
       call dbgciv('mulss ',r(lbas+1),niz,njz,1,niz)
      endif
      call symm1a(z(izzz+ibassi(iooz,1)),r(lbas+1),niz)
      goto 200
c... z vector transposed
301   nnnz=1
      mmmz=njz
      jump=ibassi(jooz,isz)+izzz
      goto 888
c... z vector not transposed
302   mmmz=1
      nnnz=niz
      jump=ibassi(iooz,isz)+izzz
888   if(niz.ge.njz)goto 1000
c... evaluate z (transpose)
      call mxmb(r(kbas+1),nnnc,mmmc,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),z(jump),nnnz,mmmz,njz,njo,niz)
      if(odebug(40)) then
       call dbgciv('mulss ',z(jump),njz,niz,nnnz,mmmz)
      endif
      goto 200
c... evaluate z
1000  call mxmb(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(kbas+1),
     *mmmc,nnnc,z(jump),mmmz,nnnz,niz,njo,njz)
      if(odebug(40)) then
       call dbgciv('mulss ',z(jump),niz,njz,mmmz,nnnz)
      endif
200   continue
      ibbb=ibbb+mrowb
100   izzz=izzz+norz
999   return
      end
      subroutine multt(oper,c,lamc,norsic,isc,
     *z,lamz,norz,isz,b,mcolb,mrowb,r)
      implicit REAL  (a-h,o-z),integer (i-n)
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      logical lmode,iszisc
      dimension b(*),oper(*),c(*),z(*),r(*)
INCLUDE(common/table)
INCLUDE(common/helpr)
      common/expanc/ibasso(8),mcolo(8),mrowo(8)
INCLUDE(common/symci)
c
c... z(triplet)  -  c(triplet)  (ij/ab)  fock(ab) interactions
      ibbb=1
      izzz=1
      iszisc=isz.eq.1.or.isc.eq.1
      if(lamz.le.lamc)goto 444
      do 1100 lam=1,lamc
      do 1200 iooz=1,nirr
      niz=norbe(iooz)
      if(niz)1201,1200,1201
1201  jooz=mult(iooz,isz)
      njz=norbe(jooz)
      if(njz)1202,1200,1202
1202  if(njz.eq.1.and.iszisc)goto 1200
      jooo=mult(jooz,isc)
      njo=norbe(jooo)
      if(njo)1203,1200,1203
1203  lbas=niz*njz
      jooioo=jooz-iooz
      if(jooioo)1301,1302,1302
c...   z  vector  transposed
1301   nnnz=1
       mmmz=njz
       jump=ibastr(jooz,isz)
       lmode=.true.
      goto 1888
c... z vector not transposed
1302   mmmz=1
       nnnz=niz
       jump=ibastr(iooz,isz)
       lmode=.false.
1888   if(jooz-jooo)1211,1210,1212
c... c vector totally symmetric
1210   call sqtrip(r(lbas+1),c(izzz+ibastr(jooo,1)),njo)
       if(niz.ge.njz)goto 1400
       call mxmd(r(lbas+1),1,njo,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),r,nnnz,mmmz,njz,njo,niz)
      lmode=.not.lmode
      goto 1401
1400  call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(lbas+1),
     *1,njo,r,mmmz,nnnz,niz,njo,njz)
      goto 1401
c... c vector transposed
1211  nnnc=1
      mmmc=njz
      jum=ibastr(jooz,isc)+izzz
      lmode=.not.lmode
      goto 1213
c... c vector not transposed
1212  mmmc=1
      nnnc=njo
      jum=ibastr(jooo,isc)+izzz
1213  if(niz.ge.njz)goto 1214
      call mxmd(c(jum),nnnc,mmmc,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),r,nnnz,mmmz,njz,njo,niz)
      goto 1401
1214  call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),c(jum),
     *mmmc,nnnc,r,mmmz,nnnz,niz,njo,njz)
1401  kbas=0
      if(jooioo)1403,1402,1403
1402  kbas=lbas
      lbas=iky(niz)
      call anti1(r(kbas+1),r,niz)
1403  if(lmode)goto 14400
      call mxmb(r(kbas+1),1,0,b(ibbb),1,mrowb,z(jump+1),1,norz,
     *lbas,1,lamz)
      if(odebug(40)) then
       call dbgciv('multt ',z(jump+1),lbas,lamz,1,norz)
      endif
      goto 1200
14400 call mxmbn(r(kbas+1),1,0,b(ibbb),1,mrowb,z(jump+1),1,norz,
     *lbas,1,lamz)
      if(odebug(40)) then
       call dbgciv('multt ',z(jump+1),lbas,lamz,1,norz)
      endif
1200  continue
      izzz=izzz+norsic
1100  ibbb=ibbb+mcolb
      goto 999
444   do 100 lam=1,lamz
      if(mrowb)112,111,112
c... fock operator contribution
 111  continue
      call dcopy(norsic,c(izzz),1,r,1)
      goto 113
112   call mxmd(c,1,norsic,b(ibbb),mcolb,mrowb,r,1,1,
     *norsic,lamc,1)
113   do 200 iooz=1,nirr
      niz=norbe(iooz)
      if(niz)201,200,201
201   jooz=mult(iooz,isz)
      njz=norbe(jooz)
      if(njz)202,200,202
202   if(njz.eq.1.and.iszisc)goto 200
      jooo=mult(jooz,isc)
c...   z(iooz,jooz) = oper(iooz,jooo) * c(jooo,jooz)
      njo=norbe(jooo)
      if(njo)203,200,203
203   lbas=norsic
      lmode=.true.
      if(jooz-jooo)211,210,212
c... c vector totally symmetric
210   call sqtrip(r(lbas+1),r(ibastr(jooo,1)+1),njo)
      kbas=lbas
      lbas=lbas+njo*njo
      goto 213
c... c transposed
211   kbas=ibastr(jooz,isc)
      mmmc=njz
      nnnc=1
      lmode=.false.
      goto 214
c... c not transposed
212   kbas=ibastr(jooo,isc)
213   nnnc=njo
      mmmc=1
214   if(jooz-iooz)301,300,302
c...  z vector totally symmetric
300   call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(kbas+1),
     *mmmc,nnnc,r(lbas+1),1,niz,niz,njo,njz)
      if(lmode)goto 1300
      call anti1s(z(izzz+ibastr(iooz,1)),r(lbas+1),niz)
      goto 200
1300  call anti1a(z(izzz+ibastr(iooz,1)),r(lbas+1),niz)
      goto 200
c... z vector transposed
301   nnnz=1
      mmmz=njz
      jump=ibastr(jooz,isz)+izzz
      lmode=.not.lmode
      goto 888
c... z vector not transposed
302   mmmz=1
      nnnz=niz
      jump=ibastr(iooz,isz)+izzz
888   if(lmode)goto 777
      if(niz.ge.njz)goto 2000
c... evaluate z transpose
      call mxmbn(r(kbas+1),nnnc,mmmc,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),z(jump),nnnz,mmmz,njz,njo,niz)
      if(odebug(40)) then
       call dbgciv('multt ',z(jump),njz,niz,nnnz,mmmz)
      endif
      goto 200
c... evaluate z
2000  call mxmbn(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(kbas+1),
     *mmmc,nnnc,z(jump),mmmz,nnnz,niz,njo,njz)
      if(odebug(40)) then
       call dbgciv('multt ',z(jump),niz,njz,mmmz,nnnz)
      endif
      goto 200
777   if(niz.ge.njz)goto 1000
c... evaluate z (transpose)
      call mxmb(r(kbas+1),nnnc,mmmc,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),z(jump),nnnz,mmmz,njz,njo,niz)
      if(odebug(40)) then
       call dbgciv('multt ',z(jump),njz,niz,nnnz,mmmz)
      endif
      goto 200
c... evaluate z
1000  call mxmb(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(kbas+1),
     *mmmc,nnnc,z(jump),mmmz,nnnz,niz,njo,njz)
      if(odebug(40)) then
       call dbgciv('multt ',z(jump),niz,njz,mmmz,nnnz)
      endif
200   continue
      ibbb=ibbb+mrowb
100   izzz=izzz+norz
999   return
      end
      subroutine nulst(oper,c,lamc,norsic,isc,
     *z,lamz,norz,isz,t,lamt,nort,b,mcolb,mrowb,bb,mcolbb,mrowbb,r)
      implicit REAL  (a-h,o-z),integer  (i-n)
      logical lmode
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      dimension bb(*),t(*),b(*),oper(*),c(*),z(*),r(*)
INCLUDE(common/table)
INCLUDE(common/helpr)
      common/expanc/ibasso(8),mcolo(8),mrowo(8)
INCLUDE(common/symci)
c
c... z(singlet + triplet) - c(singlet)  (ia/jb) interactions
      if(lamc)4000,999,4000
4000  ibbbb=1
      ibbb=1
      izzz=1
      do 1100 lam=1,lamc
      do 1200 iooz=1,nirr
      niz=norbe(iooz)
      if(niz)1201,1200,1201
1201  jooz=mult(iooz,isz)
      njz=norbe(jooz)
      if(njz)1202,1200,1202
1202  jooo=mult(jooz,isc)
      njo=norbe(jooo)
      if(njo)1203,1200,1203
1203  lbas=niz*njz
      jooioo=jooz-iooz
      if(jooioo)1301,1302,1302
c...   z  vector  transposed
1301   nnnz=1
       mmmz=njz
       jump=ibassi(jooz,isz)
       jumpt=jump
       lmode=.true.
      goto 1888
c... z vector not transposed
1302   mmmz=1
       nnnz=niz
       jump=ibassi(iooz,isz)
       jumpt=ibastr(iooz,isz)
       lmode=.false.
1888   if(jooz-jooo)1211,1210,1212
c... c vector totally symmetric
1210   call square(r(lbas+1),c(izzz+ibassi(jooo,1)),njo,njo)
       if(niz.ge.njz)goto 1400
       call mxmd(r(lbas+1),1,njo,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),r,nnnz,mmmz,njz,njo,niz)
      goto 1401
1400  call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(lbas+1),
     *1,njo,r,mmmz,nnnz,niz,njo,njz)
      goto 1401
c... c vector transposed
1211  nnnc=1
      mmmc=njz
      jum=ibassi(jooz,isc)+izzz
      goto 1213
c... c vector not transposed
1212  mmmc=1
      nnnc=njo
      jum=ibassi(jooo,isc)+izzz
1213  if(niz.ge.njz)goto 1214
      call mxmd(c(jum),nnnc,mmmc,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),r,nnnz,mmmz,njz,njo,niz)
      goto 1401
1214  call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),c(jum),
     *mmmc,nnnc,r,mmmz,nnnz,niz,njo,njz)
1401  kbas=0
      if(jooioo)1403,1402,1403
1402  kbas=lbas
      lbas=iky(niz+1)
      call symm1(r(kbas+1),r,niz)
1403  call mxmb(r(kbas+1),1,0,b(ibbb),1,mrowb,z(jump+1),1,norz,
     *lbas,1,lamz)
      if(odebug(40)) then
       call dbgciv('nulst ',z(jump+1),lbas,lamz,1,norz)
      endif
      if(jooioo)1503,1502,1503
1502  if(niz.eq.1)goto 1200
      lbas=iky(niz)
      call anti1(r(kbas+1),r,niz)
1503  if(lmode)goto 15500
      call mxmb(r(kbas+1),1,0,bb(ibbbb),1,mrowbb,t(jumpt+1),1,nort,
     *lbas,1,lamt)
      if(odebug(40)) then
       call dbgciv('nulst ',t(jumpt+1),lbas,lamt,1,nort)
      endif
      goto 1200
15500 call mxmbn(r(kbas+1),1,0,bb(ibbbb),1,mrowbb,t(jumpt+1),
     *1,nort,lbas,1,lamt)
      if(odebug(40)) then
       call dbgciv('nulst ',t(jumpt+1),lbas,lamt,1,nort)
      endif
1200  continue
      izzz=izzz+norsic
      ibbbb=ibbbb+mcolbb
1100  ibbb=ibbb+mcolb
999   return
      end
      subroutine nults(oper,c,lamc,norsic,isc,
     *s,lams,nors,isz,z,lamz,norz,bb,mcolbb,mrowbb,b,mcolb,mrowb,r)
      implicit REAL  (a-h,o-z),integer  (i-n)
      logical lmode,mmode,isceq1
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      dimension bb(*),s(*),b(*),oper(*),c(*),z(*),r(*)
INCLUDE(common/table)
INCLUDE(common/helpr)
      common/expanc/ibasso(8),mcolo(8),mrowo(8)
INCLUDE(common/symci)
c
c... z(singlet + triplet)  -   c(triplet) (ia/jb) interactions
      if(lamc)4000,999,4000
4000  isceq1=isc.eq.1
      ibbbb=1
      ibbb=1
      izzz=1
      do 1100 lam=1,lamc
      do 1200 iooz=1,nirr
      niz=norbe(iooz)
      if(niz)1201,1200,1201
1201  jooz=mult(iooz,isz)
      njz=norbe(jooz)
      if(njz)1202,1200,1202
1202  if(njz.eq.1.and.isceq1)goto 1200
      jooo=mult(jooz,isc)
      njo=norbe(jooo)
      if(njo)1203,1200,1203
1203  lbas=niz*njz
      jooioo=jooz-iooz
      if(jooioo)1301,1302,1302
c...   z  vector  transposed
1301   nnnz=1
       mmmz=njz
       jump=ibastr(jooz,isz)
       jumps=jump
       lmode=.true.
      goto 1888
c... z vector not transposed
1302   mmmz=1
       nnnz=niz
       jump=ibastr(iooz,isz)
       jumps=ibassi(iooz,isz)
       lmode=.false.
1888   if(jooz-jooo)1211,1210,1212
c... c vector totally symmetric
1210   if(njo.eq.1)goto 1200
       call sqtrip(r(lbas+1),c(izzz+ibastr(jooo,1)),njo)
       if(niz.ge.njz)goto 1400
       call mxmd(r(lbas+1),1,njo,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),r,nnnz,mmmz,njz,njo,niz)
      lmode=.not.lmode
      mmode=.true.
      goto 1401
1400  call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(lbas+1),
     *1,njo,r,mmmz,nnnz,niz,njo,njz)
       mmode=.false.
      goto 1401
c... c vector transposed
1211  nnnc=1
      mmmc=njz
      jum=ibastr(jooz,isc)+izzz
      lmode=.not.lmode
      mmode=.true.
      goto 1213
c... c vector not transposed
1212  mmmc=1
      nnnc=njo
      jum=ibastr(jooo,isc)+izzz
       mmode=.false.
1213  if(niz.ge.njz)goto 1214
      call mxmd(c(jum),nnnc,mmmc,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),r,nnnz,mmmz,njz,njo,niz)
      goto 1401
1214  call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),c(jum),
     *mmmc,nnnc,r,mmmz,nnnz,niz,njo,njz)
1401  kbas=0
      if(jooioo)1403,1402,1403
1402  kbas=lbas
      if(niz.eq.1)goto 1440
      lbas=iky(niz)
      call anti1(r(kbas+1),r,niz)
1403  if(lmode)goto 14400
      call mxmb(r(kbas+1),1,0,b(ibbb),1,mrowb,z(jump+1),
     *1,norz,lbas,1,lamz)
      if(odebug(40)) then
       call dbgciv('nults ',z(jump+1),lbas,lamz,1,norz)
      endif
      goto 1404
14400 call mxmbn(r(kbas+1),1,0,b(ibbb),1,mrowb,z(jump+1),
     *1,norz,lbas,1,lamz)
      if(odebug(40)) then
       call dbgciv('nults ',z(jump+1),lbas,lamz,1,norz)
      endif
1404   if(jooioo)1503,1440,1503
1440   lbas=iky(niz+1)
       call symm1(r(kbas+1),r,niz)
1503  if(mmode)goto 15500
      call mxmb(r(kbas+1),1,0,bb(ibbbb),1,mrowbb,s(jumps+1),1,nors,
     *lbas,1,lams)
      if(odebug(40)) then
       call dbgciv('nults ',s(jumps+1),lbas,lams,1,nors)
      endif
      goto 1200
15500 call mxmbn(r(kbas+1),1,0,bb(ibbbb),1,mrowbb,s(jumps+1),1,nors,
     *lbas,1,lams)
      if(odebug(40)) then
       call dbgciv('nults ',s(jumps+1),lbas,lams,1,nors)
      endif
1200  continue
      izzz=izzz+norsic
      ibbbb=ibbbb+mcolbb
1100  ibbb=ibbb+mcolb
999   return
      end
      subroutine mulst(oper,c,lamc,norsic,isc,t,lamt,nortrc
     *,z,lamz,norz,isz,b,mcolb,mrowb,bb,mcolbb,mrowbb,r)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/prnprn)
      dimension oper(*),b(*),t(*),bb(*),c(*),z(*),r(*)
INCLUDE(common/table)
      common/expanc/ibasso(8),mcolo(8),mrowo(8)
INCLUDE(common/symci)
c
c... z(singlet)  -  c(singlet+triplet)  (ia/jb) interactions
      if(lamz)3008,9999,3008
3008  if(lamt)3009,3002,3009
3002  call mulss(oper,c,lamc,norsic,isc,z,lamz,norz,isz,b,mcolb,
     *mrowb,r)
      goto 9999
3009  ibbb=1
      izzz=1
      ibbbb=1
      mu=nortrc+norsic
      do 100 lam=1,lamz
      if (lamc.eq.0) then
         call vclr(r,1,norsic)
      else
         call mxmd(c,1,norsic,b(ibbb),mcolb,mrowb,r,1,1,
     +             norsic,lamc,1)
      end if
      call mxmd(t,1,nortrc,bb(ibbbb),mcolbb,mrowbb,r(norsic+1),1,1,
     +          nortrc,lamt,1)
      do 200 iooz=1,nirr
      niz=norbe(iooz)
      if(niz)201,200,201
201   jooz=mult(iooz,isz)
      njz=norbe(jooz)
      if(njz)202,200,202
202   jooo=mult(jooz,isc)
c...   z(iooz,jooz) = oper(iooz,jooo) * c(jooo,jooz)
      njo=norbe(jooo)
      if(njo)203,200,203
203   njonjz=njo*njz
      lbas=mu+njonjz
      if(jooz-jooo)211,210,212
c... c vector totally symmetric
210   if(njo.ne.1)goto 209
      r(mu+1)=r(ibassi(jooo,1)+1)
      goto 213
209   call sqsitr(r(mu+1),r(ibassi(jooo,1)+1),
     *r(ibastr(jooo,1)+norsic+1),njo)
      goto 213
c... c transposed
211   kbas=ibassi(jooz,isc)
      mmmc=njz
      nnnc=1
      call vsub(r(kbas+1),1,r(kbas+norsic+1),1,r(mu+1),1,njonjz)
      goto 214
c... c not transposed
212   kbas=ibassi(jooo,isc)
      call vadd(r(kbas+1),1,r(kbas+norsic+1),1,r(mu+1),1,njonjz)
213   nnnc=njo
      mmmc=1
214   if(jooz-iooz)301,300,302
c...  z vector totally symmetric
300   call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(mu+1),
     *mmmc,nnnc,r(lbas+1),1,niz,niz,njo,njz)
      if(odebug(40)) then
       call dbgciv('mulst ',r(lbas+1),niz,njz,1,niz)
      endif
      call symm1a(z(izzz+ibassi(iooz,1)),r(lbas+1),niz)
      goto 200
c... z vector transposed
301   nnnz=1
      mmmz=njz
      jump=ibassi(jooz,isz)+izzz
      goto 888
c... z vector not transposed
302   mmmz=1
      nnnz=niz
      jump=ibassi(iooz,isz)+izzz
888   if(niz.ge.njz)goto 1000
c... evaluate z (transpose)
      call mxmb(r(mu+1),nnnc,mmmc,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),z(jump),nnnz,mmmz,njz,njo,niz)
      if(odebug(40)) then
       call dbgciv('mulst ',z(jump),njz,niz,nnnz,mmmz)
      endif
      goto 200
c... evaluate z
1000  call mxmb(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(mu+1),
     *mmmc,nnnc,z(jump),mmmz,nnnz,niz,njo,njz)
      if(odebug(40)) then
       call dbgciv('mulst ',z(jump),niz,njz,mmmz,nnnz)
      endif
200   continue
      ibbbb=ibbbb+mrowbb
      ibbb=ibbb+mrowb
100   izzz=izzz+norz
9999  return
      end
      subroutine mults(oper,c,lamc,norsic,isc,t,lamt,nortrc
     *,z,lamz,norz,isz,b,mcolb,mrowb,bb,mcolbb,mrowbb,r)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/prnprn)
      dimension oper(*),b(*),t(*),bb(*),c(*),z(*),r(*)
INCLUDE(common/table)
      common/expanc/ibasso(8),mcolo(8),mrowo(8)
INCLUDE(common/symci)
c
c... z(triplet)  -  c(singlet+triplet)  (ia/jb) interactions
      if(lamz)3008,9999,3008
3008   if(lamc)3001,3000,3001
3000   call multt(oper,t,lamt,nortrc,isc,z,lamz,norz,isz,bb,
     *mcolbb,mrowbb,r)
       goto 9999
3001  ibbb=1
      izzz=1
      ibbbb=1
      mu=nortrc+norsic
      do 100 lam=1,lamz
      call mxmd(c,1,norsic,b(ibbb),mcolb,mrowb,r,1,1,
     *norsic,lamc,1)
      if(lamt)9988,9977,9988
9988  call mxmd(t,1,nortrc,bb(ibbbb),mcolbb,mrowbb,r(norsic+1),1,1,
     *nortrc,lamt,1)
9977  do 200 iooz=1,nirr
      niz=norbe(iooz)
      if(niz)201,200,201
201   jooz=mult(iooz,isz)
      njz=norbe(jooz)
      if(njz)202,200,202
202   jooo=mult(jooz,isc)
c...   z(iooz,jooz) = oper(iooz,jooo) * c(jooo,jooz)
      njo=norbe(jooo)
      if(njo)203,200,203
203    njonjz=njo*njz
       lbas=mu+njonjz
      if(jooz-jooo)211,210,212
c... c vector totally symmetric
210   if(njo.ne.1)goto 209
       r(mu+1)=r(ibassi(jooo,1)+1)
      goto 213
209   call sqsitr(r(mu+1),r(ibassi(jooo,1)+1),
     *r(ibastr(jooo,1)+norsic+1),njo)
      goto 213
c... c transposed
211   kbas=ibassi(jooz,isc)
      mmmc=njz
      nnnc=1
      call vsub(r(kbas+1),1,r(kbas+norsic+1),1,r(mu+1),1,njonjz)
      goto 214
c... c not transposed
212   kbas=ibassi(jooo,isc)
      call vadd(r(kbas+1),1,r(kbas+norsic+1),1,r(mu+1),1,njonjz)
213   nnnc=njo
      mmmc=1
214   if(jooz-iooz)301,300,302
c...  z vector totally symmetric
300   if(niz.eq.1)goto 200
      call mxmd(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(mu+1),
     *mmmc,nnnc,r(lbas+1),1,niz,niz,njo,njz)
      if(odebug(40)) then
       call dbgciv('mults ',r(lbas+1),niz,njz,1,niz)
      endif
      call anti1a(z(izzz+ibastr(iooz,1)),r(lbas+1),niz)
      goto 200
c... z vector transposed
301   if(niz.ge.njz)goto 2000
c... evaluate z (transpose)
      call mxmbn(r(mu+1),nnnc,mmmc,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),z(izzz+ibassi(jooz,isz)),1,njz,njz,njo,niz)
      if(odebug(40)) then
       call dbgciv('mults ',z(izzz+ibassi(jooz,isz)),njz,niz,1,njz)
      endif
      goto 200
c... evaluate z
2000  call mxmbn(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(mu+1),
     *mmmc,nnnc,z(izzz+ibassi(jooz,isz)),njz,1,niz,njo,njz)
      if(odebug(40)) then
       call dbgciv('mults ',z(izzz+ibassi(jooz,isz)),niz,njz,njz,1)
      endif
      goto 200
c... z vector not transposed
302   if(niz.ge.njz)goto 1000
c... evaluate z (transpose)
      call mxmb(r(mu+1),nnnc,mmmc,oper(ibasso(iooz)+1),mrowo(iooz),
     *mcolo(iooz),z(izzz+ibassi(iooz,isz)),niz,1,njz,njo,niz)
      if(odebug(40)) then
       call dbgciv('mults ',z(izzz+ibassi(iooz,isz)),njz,niz,niz,1)
      endif
      goto 200
c... evaluate z
1000  call mxmb(oper(ibasso(iooz)+1),mcolo(iooz),mrowo(iooz),r(mu+1),
     *mmmc,nnnc,z(izzz+ibassi(iooz,isz)),1,niz,niz,njo,njz)
      if(odebug(40)) then
       call dbgciv('mults ',z(izzz+ibassi(iooz,isz)),niz,njz,1,niz)
      endif
200   continue
      ibbbb=ibbbb+mrowbb
      ibbb=ibbb+mrowb
100   izzz=izzz+norz
9999  return
      end
      subroutine ssfock(cr)
c... n-2 / n-2 fock(ab) interactions
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension cr(*)
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt,
     *khcsc,khzsc,khctc,khztc,isc,khcsr,khzsr,khctr,khztr,isr,
     *norsir,nortrr
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
INCLUDE(common/sizes)
INCLUDE(common/symci)
      common/expanc/ibasso(8),mcolo(8),mrowo(8)
INCLUDE(common/helpr)
INCLUDE(common/table)
      do 1 loop=1,nirr
      nec=norbe(loop)
      if(nec)2,1,2
2     call getsc(cr(khb8+1),iky(nec+1))
      ibasso(loop)=ibassq(loop,1)
      call square(cr(khb7+ibasso(loop)+1),cr(khb8+1),nec,nec)
      mcolo(loop)=1
      mrowo(loop)=nec
1     continue
      if(lamsc)891,890,891
891   call mulss(cr(khb7+1),cr(khcsc+1),lamsc,ns,isc,
     *cr(khzsc+1),lamsc,ns,isc,cr,0,0,cr(khb8+1))
      if(lamtc)890,999,890
890   call multt(cr(khb7+1),cr(khctc+1),lamtc,nt,isc,
     *cr(khztc+1),lamtc,nt,isc,cr,0,0,cr(khb8+1))
999   return
      end
      subroutine ssiaib(cr)
      implicit REAL  (a-h,o-z),integer  (i-n)
c... n-2 / n-2 (ia/ib) interactions
      dimension cr(*)
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt,
     *khcsc,khzsc,khctc,khztc,isc,khcsr,khzsr,khctr,khztr,isr,
     *norsir,nortrr
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
      common/spew/rmode(4)
INCLUDE(common/symci)
INCLUDE(common/table)
      common/expanc/ibasso(8),mcolo(8),mrowo(8)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      equivalence (rmode(4),kkkbas)
      do 1 loop=1,nirr
      nec=norbe(loop)
      if(nec)2,1,2
2     ibasso(loop)=ibassq(loop,1)
      call square(cr(khb7+ibasso(loop)+1),cr(kkkbas+ibassi(loop,1)+1)
     *,nec,nec)
      mcolo(loop)=1
      mrowo(loop)=nec
1     continue
c... z(c) singlet
      call mulst(cr(khb7+1),cr(khcsc+1),lamsc,ns,isc,
     *cr(khctc+1),lamtc,nt,cr(khzsc+1),lamsc,ns,isc,cr(khb3+1),
     *1,lamsc,cr(khb5+1),1,lamtc,cr(khb8+1))
c... z(c) triplet
      call mults(cr(khb7+1),cr(khcsc+1),lamsc,ns,isc,cr(khctc+1),
     *lamtc,nt,cr(khztc+1),lamtc,nt,isc,cr(khb5+1),lamtc,1,
     *cr(khb4+1),1,lamtc,cr(khb8+1))
      return
      end
      subroutine ssijab(cr)
      implicit REAL  (a-h,o-z),integer  (i-n)
c... n-2 / n-2 (ij/ab) interactions
      dimension cr(*)
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt,
     *khcsc,khzsc,khctc,khztc,isc,khcsr,khzsr,khctr,khztr,isr,
     *norsir,nortrr
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
      common/spew/rmode(4)
INCLUDE(common/symci)
INCLUDE(common/table)
      common/expanc/ibasso(8),mcolo(8),mrowo(8)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir),
     *(rmode(4),kkkbas)
       iso=mult(isc,isr)
       if(iso.ne.1)goto 3
      do 1 loop=1,nirr
      nec=norbe(loop)
      if(nec)2,1,2
2     ibasso(loop)=ibassq(loop,1)
      call square(cr(khb7+ibasso(loop)+1),cr(kkkbas+ibassi(loop,1)+1)
     *,nec,nec)
      mcolo(loop)=1
      mrowo(loop)=nec
1     continue
      kbas=khb7
       goto 777
3      do 5 loop=1,nirr
      moop=mult(loop,iso)
      if(loop.ge.moop)goto 5
      ibasso(loop)=ibassi(loop,iso)
      ibasso(moop)=ibasso(loop)
      mcolo(loop)=1
      mrowo(moop)=1
      mrowo(loop)=norbe(loop)
      mcolo(moop)=norbe(loop)
5     continue
      kbas=kkkbas
777   if(lamss)891,890,891
c... z(c) singlet
891   call mulss(cr(kbas+1),cr(khcsr+1),lamsr,norsir,isr,
     *cr(khzsc+1),lamsc,ns,isc,cr(khb3+1),1,lamsr,cr(khb8+1))
c... z(r) singlet
      call mulss(cr(kbas+1),cr(khcsc+1),lamsc,ns,isc,cr(khzsr+1),
     *lamsr,norsir,isr,cr(khb3+1),lamsr,1,cr(khb8+1))
      if(lamtt)890,999,890
c... z(c) triplet
890   call multt(cr(kbas+1),cr(khctr+1),lamtr,nortrr,isr,
     *cr(khztc+1),lamtc,nt,isc,cr(khb4+1),1,lamtr,cr(khb8+1))
c... z(r) triplet
      call multt(cr(kbas+1),cr(khctc+1),lamtc,nt,isc,cr(khztr+1),
     *lamtr,nortrr,isr,cr(khb4+1),lamtr,1,cr(khb8+1))
999   return
      end
      subroutine ssiajb(cr)
      implicit REAL  (a-h,o-z),integer (i-n)
c... n-2 / n-2 (ia/jb) interactions
      dimension cr(*)
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt,
     *khcsc,khzsc,khctc,khztc,isc,khcsr,khzsr,khctr,khztr,isr,
     *norsir,nortrr,lamstc
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88
      common/spew/rmode(5)
INCLUDE(common/symci)
INCLUDE(common/table)
      common/expanc/ibasso(8),mcolo(8),mrowo(8)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      equivalence (rmode(4),kkkbas),(rmode(5),itrans)
       iso=mult(isc,isr)
        kkkbas=kkkbas+kh
       lamrc=lamsr+lamtr-lamstc
      if(itrans)665,664,664
c... operator transposed
664   do 2 loop=1,nirr
      moop=mult(loop,iso)
       ibasso(loop)=ibassq(moop,iso)
       mrowo(loop)=1
2      mcolo(loop)=norbe(moop)
       goto 666
c... operator not transposed
665    do 1 loop=1,nirr
       ibasso(loop)=ibassq(loop,iso)
      mcolo(loop)=1
1      mrowo(loop)=norbe(loop)
666    if(lamrc)890,891,891
c... z(c) - singlet r
890   call nulst(cr(kkkbas+1),cr(khcsr+1),lamsr,norsir,isr,
     *cr(khzsc+1),lamsc,ns,isc,cr(khztc+1),lamtc,nt,
     *cr(khb1+1),1,lamsr,cr(khb66+1),1,lamsr,cr(khb88+1))
c... z(c)  -  triplet(r)
      call nults(cr(kkkbas+1),cr(khctr+1),lamtr,nortrr,isr,
     *cr(khzsc+1),lamsc,ns,isc,cr(khztc+1),lamtc,nt,
     *cr(khb55+1),1,lamtr,cr(khb44+1),1,lamtr,cr(khb88+1))
      goto 892
c... z(c) singlet
891    call mulst(cr(kkkbas+1),cr(khcsr+1),lamsr,norsir,isr,
     *cr(khctr+1),lamtr,nortrr,cr(khzsc+1),lamsc,ns,isc,
     *cr(khb1+1),1,lamsr,cr(khb55+1),1,lamtr,cr(khb88+1))
c... z(c) triplet
       call mults(cr(kkkbas+1),cr(khcsr+1),lamsr,norsir,isr,
     *cr(khctr+1),lamtr,nortrr,cr(khztc+1),lamtc,nt,isc,
     *cr(khb66+1),1,lamsr,cr(khb44+1),1,lamtr,cr(khb88+1))
892   if(itrans)765,764,764
765   do 22 loop=1,nirr
      moop=mult(loop,iso)
      ibasso(loop)=ibassq(moop,iso)
      mrowo(loop)=1
22      mcolo(loop)=norbe(moop)
       goto 888
764    do 11 loop=1,nirr
        ibasso(loop)=ibassq(loop,iso)
        mcolo(loop)=1
11      mrowo(loop)=norbe(loop)
888   if(lamrc)291,291,290
c... z(r)  - singlet c
290   call nulst(cr(kkkbas+1),cr(khcsc+1),lamsc,ns,isc,
     *cr(khzsr+1),lamsr,norsir,isr,cr(khztr+1),lamtr,nortrr,
     *cr(khb1+1),lamsr,1,cr(khb55+1),lamtr,1,cr(khb88+1))
c...  z(r)  -  triplet c
      call nults(cr(kkkbas+1),cr(khctc+1),lamtc,nt,isc,
     *cr(khzsr+1),lamsr,norsir,isr,cr(khztr+1),lamtr,nortrr,
     *cr(khb66+1),lamsr,1,cr(khb44+1),lamtr,1,cr(khb88+1))
      goto 999
c...   z(r) singlet
291   call mulst(cr(kkkbas+1),cr(khcsc+1),lamsc,ns,isc,
     *cr(khctc+1),lamtc,nt,cr(khzsr+1),lamsr,norsir,isr,
     *cr(khb1+1),lamsr,1,cr(khb66+1),lamsr,1,cr(khb88+1))
c... z(r) triplet
      call mults(cr(kkkbas+1),cr(khcsc+1),lamsc,ns,isc,
     *cr(khctc+1),lamtc,nt,cr(khztr+1),lamtr,nortrr,isr,
     *cr(khb55+1),lamtr,1,cr(khb44+1),lamtr,1,cr(khb88+1))
999    return
       end
      subroutine ssabcd(z,g,c,nspipc,norsic,iblk)
      implicit REAL  (a-h,o-z), integer (i-n)
      dimension z(*),g(*),c(*)
INCLUDE(common/prnprn)
INCLUDE(common/symci)
_IF(parallel)
        logical odo,oipsci
_ENDIF
c
c...  n-2 / n-2 (ab/cd) interactions
       if(nspipc)888,999,888
888    j=1
      js=1
      jl=n127
       idim=n127
      jdim=n127
      nwin=n127sq
      call brtsb(iblk)
_IF(parallel)
       iflop=iipsci()
_ENDIF
444    itop=norsic-jl
_IF(parallel)
        odo=.not.oipsci()
_ENDIF
       if(itop)666,555,555
c... border blocks
666    jdim=n127+itop
       nwin=jdim*n127
555    is=1
_IF(parallel)
       if (.not.odo) goto 2
_ENDIF
       do 1 i=1,j
      call getsb(g(1),nwin)
      if(i.ne.j)goto 333
c... diagonal block
      idim=jdim
      goto 222
c... use g in transposed form
333   continue
      call mxmb(g,n127,1,c(is),1,norsic,z(js),1,norsic,
     *jdim,n127,nspipc)
      if(odebug(40)) then
       call dbgciv('ssabcd',z(js),jdim,nspipc,1,norsic)
      endif
c... use g straight
222   continue
      call mxmb(g,1,n127,c(js),1,norsic,z(is),1,norsic,
     *idim,jdim,nspipc)
      if(odebug(40)) then
       call dbgciv('ssabcd',z(is),idim,nspipc,1,norsic)
      endif
1      is=is+n127
2      j=j+1
       js=js+n127
      jl=jl+n127
       if(js.le.norsic)goto 444
999   return
      end
      subroutine timbrk
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/timanb)
      top=cpulft(1)
      timsec(nodel-1)=top-tmlast+timsec(nodel-1)
      tmlast=top
      call walltime(top)
      wtimsec(nodel-1)=top-wtmlast+timsec(nodel-1)
      wtmlast=top
      return
      end
      subroutine timprt(iwr)
      implicit REAL  (a-h,o-z),integer (i-n)
INCLUDE(common/timanb)
INCLUDE(common/cic)
INCLUDE(common/prints)
c
      if(oprint(28)) return
      frodo = 100.0d0
      bilbo=cpulft(1)
      top = bilbo-timbeg
      if (dabs(top).le.1.0d-3) return
      bilbo=frodo/top
      write(iwr,100)
100   format(//' diagonalizer timing analysis'/1x,28('=')//
     *' integral interaction cpu percent      wall(sec)'/
     *' -----------------------------------------------')
      do 1 loop=1,13
      top=timsec(loop)*bilbo
      frodo=frodo-top
1     write(iwr,200)text11(loop),text22(loop),top,wtimsec(loop)
200   format(1x,a8,1x,a8,f15.3,f15.3)
      write(iwr,200)text11(14),text22(14),frodo
      return
      end
      function expect(f,v,n)
      implicit REAL  (a-h,o-z), integer (i-n)
      dimension f(*),v(*)
      common/three/p(1)
      m=2
      p(1)=v(1)*f(1)
      do 1 i=2,n
      p(i)=ddot(i,v,1,f(m),1)
      call daxpy(i-1,v(i),f(m),1,p,1)
1     m=m+i
      expect=ddot(n,p,1,v,1)
      return
      end
      subroutine diagpq(d,e,ibl,norsic)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension d(*),e(*)
      common/moco/n128
INCLUDE(common/symci)
c
INCLUDE(common/disktl) 
_IF(parallel)
      logical odo,oipsci 
_ENDIF
      if(norsic)888,999,888
888   call brtsc(ibl)
_IF(parallel)
       iflop=iipsci()
_ENDIF
      j=0
      js=1
      itop=norsic
      jdim=n127
      nwin=n127sq
444   itop=itop-n127
_IF(parallel)
      odo=.not.oipsci()
_ENDIF
      if(itop)666,555,555
c... border blocks
666   jdim=n127+itop
      nwin=jdim*n127
_IF(parallel)
555   if(.not.odo)goto 2
_ELSE
555   continue
_ENDIF
      call skipsc(nwin*j)
      call getsc(e,nwin)
      call dcopy(jdim,e,n128,d(js),1)
2     js=js+jdim
      j=j+1
      if(itop)999,999,444
999   continue
      return
      end
      subroutine diagon(cr,iconf,iwr)
_IF(parallel)
c MPP
c MPP    this routine is always a pain
c MPP    since it is assumed that the diagonal is always there
c MPP    (usually not unreasonable) but wrong on distributed systems
c MPP    results are gatherd in davips (and compacted/divided)
c MPP
_ENDIF(parallel)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension cr(*),iconf(5,*)
c...
c... to build diagonal of hamiltonian
c... this routine used only if  =eigen= directive not presented
c...
c... used for both instore and outstore modes
c...
       logical luck,load,oijkl
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
INCLUDE(common/disktl)
      common/spew/rmode(5)
      common/stoctl/khbill
INCLUDE(common/cntrl)
INCLUDE(common/table)
      common/moco/n128
      common/lsort/irebas(8),multt(8),multtt(8,8),scr(1000)
INCLUDE(common/orb)
INCLUDE(common/auxh)
      common/spcoef/nwwwww(9),iwwwww(9)
INCLUDE(common/symci)
INCLUDE(common/helpr)
INCLUDE(common/ccntl)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      equivalence  (rmode(4),kkkbas),(rmode(5),icfock)
      equivalence (nornor,norbsi(1))
      equivalence (ibldub,index(198))
      luck = .false.
      load = .false.
      jakass = 0
_IF(parallel)
c MPP 
c MPP     clear index(199) / usefull for n-2
c MPP     divide work on abcd diagonal
c MPP
       call pried7
       call clrcz(iconf,index(199),iwr,numci)
       iflop = iipsci()
_ENDIF(parallel)
c
_IF(cray,convex,t3e,i8,i8drct)
      ibldub=index(199)+iconf(2,1)
_ELSE
      ibldub=index(199)+ipad(iconf(2,1))
_ENDIF
      call brtsb(index(1))
111   call getsb(rmode(1),1)
888   if(odebug(40)) then
       write(6,*)'diagon: model = ', model,ipg_nodeid()
      endif
      if (model.gt.9)goto 999
c888   goto (999,2000,3000,4000,5000,5000,6000,7000,999),model
      goto (999,2000,3000,4000,5000,5000,6000,7000,999),model
c...
c... vacuum / vacuum (ij/kl)
c...
2000   lump=khbill+nval
_IF(parallel)
c MPP
       call vclr(cr(khbill+1),1,nval)
c MPP
_ENDIF(parallel)
42     call getsb(rmode(2),2)
_IF(cray,convex,t3e,i8,i8drct)
       lams=iconf(3,ic)
_ELSE
       lams=ipad(iconf(3,ic))
_ENDIF
      if(ic.eq.ir)goto 12
_IF(cray,convex,t3e,i8,i8drct)
      klaus=lams*iconf(3,ir)
      call skipsb(klaus)
_ELSE
      call skipsb(lams*ipad(iconf(3,ir)))
_ENDIF
      goto 22
12    call getsb(cr(lump+1),lams*lams)
_IF(cray,convex,t3e,i8,i8drct)
      call dcopy(lams,cr(lump+1),lams+1,cr(iconf(1,ic)+khbill+1),1)
_ELSE
      call dcopy(lams,cr(lump+1),lams+1
     +   ,cr(ipad(iconf(1,ic))+khbill+1),1)
_ENDIF
22    call getsb(rmode(1),1)
      if(model-2)52,42,52
52    call wrt3(cr(khbill+1),nval,index(199),numci)
      if(odebug(40)) then
       call dbgvec('vv-ijkl ',cr(khbill+1),nval)
      endif
      goto 888
c...
c... doublet / doublet   (ij/kl)
c...
3000   lump=khbill+ndoub
_IF(parallel)
c MPP
      call vclr(cr(khbill+1),1,ndoub)
c MPP
_ENDIF(parallel)
13     call getsb(rmode(2),2)
_IF(cray,convex,t3e,i8,i8drct)
       lams=iconf(3,ic)
_ELSE
       lams=ipad(iconf(3,ic))
_ENDIF
       if(ic.eq.ir)goto 23
_IF(cray,convex,t3e,i8,i8drct)
       klaus=lams*iconf(3,ir)
       call skipsb(klaus)
_ELSE
       call skipsb(lams*ipad(iconf(3,ir)))
_ENDIF
       goto 33
23     call getsb(cr(lump+1),lams*lams)
_IF(cray,convex,t3e,i8,i8drct)
       icbas=khbill+iconf(1,ic)
       ne=norbe(iconf(4,ic))
_ELSE
       icbas=khbill+ipad(iconf(1,ic))
       ne=norbe(ipad(iconf(4,ic)))
_ENDIF
       lams1=lams+1
       idbas=lump
       do 43 loop=1,lams
_IF(cray,t3e)
       call setsto(ne,cr(idbas+1),cr(icbas+1))
_ELSEIF(ibm,vax)
       call setstr(ne,cr(idbas+1),cr(icbas+1))
_ELSE
       call vfill(cr(idbas+1),cr(icbas+1),1,ne)
_ENDIF
      idbas=idbas+lams1
43    icbas=icbas+ne
33    call getsb(rmode(1),1)
       if(model-3)163,13,163
163   lumq=khbill
63    call wrt3(cr(lumq+1),ndoub,ibldub,numci)
      if(odebug(40)) then
       call dbgvec('dd-ijkl ',cr(lumq+1),ndoub)
      endif
       goto 888
c...
c...  n-2 / n-2   (ij/kl)  and   (ab/cd)
c...
4000  islast=0
      lump=khbill+jbig
      lumq=lump+nubsi
      lumr=lumq+nubtr
      n128=n127+1
_IF(ga)
      ia=nstrt2
_ENDIF
14     call getsb(rmode(2),2)
      call upack2(iconf(3,ic),lamt,lams)
      if(ic.eq.ir)goto 24
      call upack2(iconf(3,ir),joop,loop)
      call skipsb(lams*loop+lamt*joop)
      goto 84
24    continue
_IF(ga)
      if (ia.eq.ic) then
         oijkl=.true.
      else
         oijkl=.false.
      endif
      call remabcd(ia,iconf,cr,lump,lumq,lumr,oijkl,islast)
      ia=ia+1
      if (.not.oijkl) goto 24
_ENDIF
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,ic)
_ELSE
      isi=ipad(iconf(4,ic))
_ENDIF
      call upack2(iconf(2,ic),kc,jc)
_IFN(ga)
      if(isi.eq.islast)goto 34
_ENDIF
      islast=isi
      norsi=norbsi(isi)
      nortr=norbtr(isi)
_IFN(ga)
c... retrieve diagonal of (ab/cd)  p  supermatrix
      call diagpq(cr(lump+1),cr(lumr+1),index(isi+9),norsi)
      if(isi.ne.1)goto 204
c... normalize diagonal of totally symmetric block
      do 304 moop=1,nirr
      ne=norbe(moop)
      if(ne)404,304,404
404   idbas=lump+ibassi(moop,1)
      do 1 loop=1,ne
      joop=iky(loop+1)+idbas
    1 cr(joop)=cr(joop)*0.5d0
304   continue
c... retrieve diagonal of  (ab/cd)  q  supermatrix
204   continue
      call diagpq(cr(lumq+1),cr(lumr+1),index(isi+17),nortr)
_ENDIF
34    if(lams)54,44,54
54    lams1=lams+1
      idbas=lumr
       icbas=khbill
      call getsb(cr(lumr+1),lams*lams)
      do 64 loop=1,lams
      frodo=cr(idbas+1)
_IF(cray,ibm,vax)
_IF1(c)cdir$ ivdep
      do 6 joop=1,norsi
    6 cr(icbas+joop)=cr(lump+joop)+frodo
_ELSE
      call vsadd(cr(lump+1),1,frodo,cr(icbas+1),1,norsi)
_ENDIF
      idbas=idbas+lams1
64    icbas=icbas+norsi
      call wrt3(cr(khbill+1),lams*norsi,index(199)+jc,numci)
      if(odebug(40)) then
       call dbgvec('n-2sijkl',cr(khbill+1),lams*norsi)
      endif
      if(lamt)44,84,44
44    lams1=lamt+1
      idbas=lumr
      icbas=khbill
      call getsb(cr(lumr+1),lamt*lamt)
      do 104 loop=1,lamt
      frodo=cr(idbas+1)
_IF(cray,ibm,vax)
_IF1(c)cdir$ ivdep
      do 7 joop=1,nortr
    7 cr(icbas+joop)=cr(lumq+joop)+frodo
_ELSE
      call vsadd(cr(lumq+1),1,frodo,cr(icbas+1),1,nortr)
_ENDIF
      idbas=idbas+lams1
104   icbas=icbas+nortr
      call wrt3(cr(khbill+1),lamt*nortr,index(199)+kc,numci)
      if(odebug(40)) then
       call dbgvec('n-2tijkl',cr(khbill+1),lamt*nortr)
      endif
84    call getsb(rmode(1),1)
_IF(ga)
      if(model-4)85,14,85
85    if (ia.gt.nlast2)goto 888
      oijkl=.false.
      do ib=ia,nlast2
         call remabcd(ib,iconf,cr,lump,lumq,lumr,oijkl,islast)
      enddo
      goto 888
_ELSE
      if(model-4)888,14,888
_ENDIF
5000  call possb(nwwwww(5),iwwwww(5))
       goto 111
c...
c...  doublet / doublet (ia/ia) and fock(aa)
c...
6000  load=.true.
       goto 88
17     lumq=lump+max(mubsi,maxpa(3)*maxpa(3))
       call rdedx(cr(lumq+1),ndoub,ibldub,numci)
       lumr=lumq+ndoub
       call brtsc(index(7))
27     call getsb(rmode(2),3)
_IF(cray,convex,t3e,i8,i8drct)
       lams=iconf(3,ic)
_ELSE
       lams=ipad(iconf(3,ic))
_ENDIF
       if(ic.eq.ir)goto 37
_IF(cray,convex,t3e,i8,i8drct)
       klaus=lams*iconf(3,ir)
       call skipsb(klaus)
_ELSE
       call skipsb(lams*ipad(iconf(3,ir)))
_ENDIF
       goto 97
37    continue
_IF(cray,convex,t3e,i8,i8drct)
       isi=iconf(4,ic)
_ELSE
       isi=ipad(iconf(4,ic))
_ENDIF
       ne=norbe(isi)
_IF(cray,convex,t3e,i8,i8drct)
       idbas=lumq+iconf(1,ic)
_ELSE
       idbas=lumq+ipad(iconf(1,ic))
_ENDIF
       if(kkkbas)57,47,57
c...  fock  operator contribution
47     call getsc(cr(lump+1),iky(ne+1))
_IF(cray)
       call gather(ne,scr,cr(lump+1),iky(2))
_ELSE
       call dgthr(ne,cr(lump+1),scr,iky(2))
_ENDIF
       do 77 loop=1,lams
      call vadd(cr(idbas+1),1,scr,1,cr(idbas+1),1,ne)
77     idbas=idbas+ne
       goto 97
c... (ia/ia) contribution
57     call getsb(cr(lump+1),lams*lams)
       call mxmb(cr(kkkbas+irebas(isi)+1),1,1,cr(lump+1),1,lams+1,
     * cr(idbas+1),1,ne,ne,1,lams)
97     call getsb(rmode(1),1)
       if(model-7)663,27,663
663   if(odebug(40)) then
       call dbgvec('dd-iaia ',cr(lumq+1),ndoub)
      endif
      go to 63
c...
c...   n-2 / n-2 (ia/ia) and fock(aa)
c...
7000  if(load)goto 18
88    call rdedx(cr(khbill+1),m2fock,index(4),numci)
c... compress exchange operators so that only diagonal present
      do 5008 loop=1,nirr
      irebas(loop)=jakass
      jak=jakass+khbill
      ne=norbe(loop)
       if(ne)128,5008,128
128   moop=khbill+ibassi(loop,1)
      do 608 joop=1,ne
      icbas=jak+joop
      idbas=moop+iky(joop+1)
_IF(cray,ibm,vax)
_IF1(c)cdir$ ivdep
       do 8 ns=1,nint
      cr(icbas)=cr(idbas)
      icbas=icbas+nornor
    8 idbas=idbas+nornor
_ELSE
      call vmov(cr(idbas),nornor,cr(icbas),nornor,nint)
_ENDIF
 608   continue
       jakass=jakass+ne
5008    continue
      lump=khbill+m2fock
      if(load)goto 17
18     call brtsc(index(6))
       moop=max(maxpa(1),maxpa(2))
       lumr=lump+max(mubsi,jbig,moop*moop)
      islast=0
      do 1008 isi=1,nirr
      do 1008 loop=1,nirr
      isf=mult(loop,isi)
      nf=norbe(isf)
      if(nf)1008,2008,1008
2008  isf=0
1008  multtt(loop,isi)=isf
28    call getsb(rmode(2),3)
      if(ic.eq.islast)goto 38
      if(islast)78,68,78
78    call upack2(iconf(2,islast),kc,jc)
      jc=index(199)+jc
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,islast)
_ELSE
      isi=ipad(iconf(4,islast))
_ENDIF
      nss=lams*norbsi(isi)
      ntt=lamt*norbtr(isi)
      call rdedx(cr(lump+1),nss+ntt,jc,numci)
      icbas=lump
      idbas=lumr
      if(isi.ne.1)goto 708
c... totally symmetric external pair
      if(nss)728,718,728
c... singlets
728   do 748 loop=1,lams
      do 748 moop=1,nirr
      ne=norbe(moop)
      if(ne)758,748,758
758   do 448 joop=1,ne
      frodo=cr(idbas+joop)
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 2 noop=1,joop
    2 cr(icbas+noop)=cr(icbas+noop)+cr(idbas+noop)+frodo
448   icbas=icbas+joop
      idbas=idbas+ne
748   continue
      if(ntt)718,738,718
c... update triplets
718   do 778 loop=1,lamt
      do 778 moop=1,nirr
      ne=norbe(moop)
      if(ne.lt.2)goto 778
      do 678 joop=2,ne
      frodo=cr(idbas+joop)
      ioop=joop-1
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 3 noop=1,ioop
    3 cr(icbas+noop)=cr(icbas+noop)+cr(idbas+noop)+frodo
678   icbas=icbas+ioop
778   idbas=idbas+ne
      goto 738
c... external pair not totally symmetric
708   do 668 loop=1,nirr
668   multt(loop)=multtt(loop,isi)
      do 788 loop=1,lamspt
      do 508 moop=1,nirr
      ne=norbe(moop)
      if(ne)518,508,518
518   isf=multt(moop)
      if(isf.lt.moop)goto 508
      nf=norbe(isf)
      iee=irebas(moop)+idbas
      iff=irebas(isf)+idbas
      do 548 joop=1,nf
      frodo=cr(iff+joop)
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
       do 4 noop=1,ne
    4  cr(icbas+noop)=cr(icbas+noop)+cr(iee+noop)+frodo
548   icbas=icbas+ne
508   continue
788   idbas=idbas+next
738   if(nss)118,108,118
118    call wrt3(cr(lump+1),nss,jc,numci)
      if(odebug(40)) then
       call dbgvec('n-2siaia',cr(lump+1),nss)
      endif
      if(ntt)108,68,108
108   call wrt3(cr(lump+nss+1),ntt,index(199)+kc,numci)
      if(odebug(40)) then
       call dbgvec('n-2tiaia',cr(lump+nss+1),ntt)
      endif
68     if(luck)goto 999
      call upack2(iconf(3,ic),lamt,lams)
      lamss=lams*lams
      lamtt=lamt*lamt
      lamst=lams*lamt
      lamspt=lams+lamt
      lamst1=lamspt-1
      lams1=lams+1
      lamt1=lamt+1
      islast=ic
38    if(ic.eq.ir)goto 138
      call upack2(iconf(3,ir),joop,loop)
      call skipsb(lams*loop+lamt*joop)
      goto 148
138   if(kkkbas)168,158,168
c... fock operator contribution
158   idbas=lumr
      do 178 moop=1,nirr
      ne=norbe(moop)
      if(ne)188,178,188
188   call getsc(cr(lump+1),iky(ne+1))
_IF(cray)
      call gather(ne,cr(idbas+1),cr(lump+1),iky(2))
_ELSE
      call dgthr(ne,cr(lump+1),cr(idbas+1),iky(2))
_ENDIF
178   idbas=idbas+ne
      if(lamst1)218,148,218
218   do 228 loop=1,lamst1
      call dcopy(next,cr(lumr+1),1,cr(idbas+1),1)
228   idbas=idbas+next
      goto 148
c... (ia/ia) contribution
168   if(lams)248,238,248
248   call getsb(cr(lump+1),lamss)
      call dcopy(lams,cr(lump+1),lams1,scr,1)
      if(lamt)238,278,238
238   call getsb(cr(lump+1),lamtt)
      call dcopy(lamt,cr(lump+1),lamt1,scr(lams1),1)
278   call mxmb(cr(kkkbas+1),1,1,scr,1,1,cr(lumr+1),1,next,
     *next,1,lamspt)
      call skipsb(lamst)
148   call getsb(rmode(1),1)
      if(model-8)288,28,288
288   luck=.true.
      goto 78
999   return
      end
      subroutine histo(cr,icr,n)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer icr
_ELSEIF(convex,i8drct)
      integer *8 icr
_ELSE
      REAL icr
_ENDIF
      dimension cr(*),icr(*)
      common /moco/ alow,xhist,diamin
      common /junk / ihist(1)
       do 1 loop=1,n
1     icr(loop)=(alow-dlog( dabs(cr(loop))+diamin))*xhist
      do 2 loop=1,n
      j=icr(loop)
2     ihist(j+1)=ihist(j+1)+1
       return
      end
      subroutine vecpri(cr,iconf,kblcf,iwr)
c...
c... ci vector printer
c... works for both instore and outstore modes
c... print occupation-patterns next to CI-coefficients
c... we do assume that cr and iconf are identical !!!!
c...
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      character * 8  vacuum,double,single,triple
      dimension cr(*),iconf(5,*)
      common/junk /ihist(13002)
c... for cdc cyber 174 use     ihist(6502)
INCLUDE(common/sizes)
INCLUDE(common/cntrl)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,
     *khbcd,khbcs,khbct
      common /moco/ alow,xhist,diamin,nbase(8)
INCLUDE(common/disktl)
INCLUDE(common/ccntl)
INCLUDE(common/orb)
INCLUDE(common/table)
INCLUDE(common/symci)
      logical malt
      common/erg/potnuc,totnuc,energ,eigen(nd200),h11,z11,
     *h22,s22,h12,s12,tester,eshift,thresh,aaa,bbb,
     *scream,gprin,z12,z22,vacpri,dubpri,stpri
     *,qester,maxcyc,moddia,mprin,malt
INCLUDE(common/mapp)
c
      logical new,modeis
      common/moco/iorb(64)
c
      data vacuum,double,single,triple/'  vacuum',' doublet',
     *' singlet',' triplet'/
c
      mzero = 0
      sos = 0.0d0
      sosm = 0.0d0
      nrefs=0
c
cjvl   read save configs back in iconf(1,nlast2+1)
cjvl   so khbill has to move (temporarily)
c
      khsave = khbill
      khbill = khbill + nlast2*5
      call rdedx(iconf(1,nlast2+1),nlast2*5,kblcf,numci)
      modeis = modei
cjvl   check if modei involves nint-next exchange if so disable
        if (.not.modei) then
           do 51 i=1,nint
              if (mapie(i).gt.nint) then
                 write(iwr,52)
52    format(/' #@#@ a reorder exchanging internal and external'
     *       ,' orbitals has taken place #@#@',/
     *       ,' #@#@ orbital numbering refers to reordered set #@#@')
                 modei = .true.
                 go to 53
              end if
51          continue
53          continue
        end if
c...
c... evaluate histogram
c...
       top=cpulft(1)
      i26=index(26)
      muck=nint
      do 2000 loop=1,nirr
      nbase(loop)=muck
2000   muck=muck+norbe(loop)
      nhist=13002
c... for cdc cyber 174  use     nhist=6502
      diamin=1.0d-5
      alow=dlog(1.02d0)
      xhist=nhist
      xhist=(xhist-2.0d0)/(alow-dlog(diamin))
_IF(cray)
      call szero(ihist,nhist)
_ELSE
      call izero(nhist,ihist,1)
_ENDIF
_IF(ga)
      call search(i26,numcz)
_ELSE
      call search(i26,numci)
_ENDIF
c... vacuum states
      if(nval)9000,1,9000
_IF(ga)
9000  call reads(cr(khbill+1),nval,numcz)
_ELSE
9000  call reads(cr(khbill+1),nval,numci)
_ENDIF
       call histo(cr(khbill+1),cr(khbill+nval+1),nval)
c... doublet states
1     if(ndoub)4,3,4
_IF(ga)
4      call reads(cr(khbill+1),ndoub,numcz)
_ELSE
4      call reads(cr(khbill+1),ndoub,numci)
_ENDIF
       call histo(cr(khbill+1),cr(khbill+ndoub+1),ndoub)
c...  singlet and triplet states
3      if(nmin2)6,5,6
6     do 700 loop=nstrt2,nlast2
      call upack2(iconf(3,loop),lamt,lams)
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,loop)
_ELSE
      isi=ipad(iconf(4,loop))
_ENDIF
_IF(ga)
      call reads(cr(khbill+1),lams*norbsi(isi),numcz)
      call reads(cr(khbill+1+lams*norbsi(isi)),lamt*norbtr(isi),numcz)
_ELSE
      lams=lams*norbsi(isi)+lamt*norbtr(isi)
      call reads(cr(khbill+1),lams,numci)
_ENDIF
      call histo(cr(khbill+1),cr(khbill+lams+1),lams)
700   continue
c...
c... histogram complete
c... now find print threshold with its aid
c...
5     muck=0
      nhist=nhist-1
      do 500 loop=1,nhist
      bilbo=loop
      muck=muck+ihist(loop)
      if(muck.ge.mprin)goto 501
500   continue
501   bilbo=dexp(-bilbo/xhist)*1.02d0-diamin
       if(bilbo.lt.gprin)bilbo=gprin
c...
c... now for actual printing
c...
      write(iwr,100)top,bilbo
 100  format(//1x,80('=')//
     *' **** ci vector printer called at',f9.2,' secs'
     *//' list of ci coefficients greater than',e14.4///
     *2x,47('*')/
     *'  coefficient internal external external internal'/
     *18x,             'spin     spin      mos     conf  -- occ.'/
     *2x,67('*')/)
      if (.not.modei) write(iwr,54)
54    format('     SPIN-COUPLING refers to REORDERED orbitals',/,
     *2x,67('*')/)
101   format(1x,f12.8,i9,1x,a8,i5,i4,i9,:,'  -- ',63i1)
_IF(ga)
      call search(i26,numcz)
_ELSE
      call search(i26,numci)
_ENDIF
c... vacuum states
      if(nnst)12,11,12
_IF(ga)
12    call reads(cr(khbill+1),nval,numcz)
_ELSE
12    call reads(cr(khbill+1),nval,numci)
_ENDIF
      hobbit=dmin1(bilbo,vacpri)
       koop=khbill
      do 900 loop=1,nnst
_IF(cray,convex,t3e,i8,i8drct)
      lams=iconf(3,loop)
_ELSE
      lams=ipad(iconf(3,loop))
_ENDIF
        new = .true.
      do 900 noop=1,lams
       koop=koop+1
      frodo=cr(koop)
      if( dabs(frodo).lt.hobbit)goto 910
        if (new) then
           call upack(iconf(1,nlast2+loop),nint,iorb)
           if (modei) then
              write(iwr,101)frodo,noop,vacuum,mzero,mzero,loop
     *                      ,(iorb(i),i=1,nint)
           else
              write(iwr,101)frodo,noop,vacuum,mzero,mzero,loop
     *                      ,(iorb(mapei(i)),i=1,nint)
           end if
           new = .false.
        else
           write(iwr,101)frodo,noop,vacuum,mzero,mzero,loop
        end if
910   continue
      if (loop .le. nref) then
      nrefs = nrefs+1
      sosm = sosm +frodo*frodo
      endif
900   sos=frodo*frodo+sos
c... doublet states
11    if(nmin1)14,13,14
_IF(ga)
14    call reads(cr(khbill+1),ndoub,numcz)
_ELSE
14    call reads(cr(khbill+1),ndoub,numci)
_ENDIF
      hobbit=dmin1(bilbo,dubpri)
       koop=khbill
      do 901 loop=nstrt1,nlast1
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,loop)
_ELSE
      isi=ipad(iconf(4,loop))
_ENDIF
      nj=norbe(isi)
      mj=nbase(isi)
_IF(cray,convex,t3e,i8,i8drct)
      lams=iconf(3,loop)
_ELSE
      lams=ipad(iconf(3,loop))
_ENDIF
      new = .true.
      do 901 noop=1,lams
      do 901 moop=1,nj
      koop=koop+1
      frodo=cr(koop)
      if( dabs(frodo).lt.hobbit)goto 901
      kc=moop+mj
        if (new) then
           call upack(iconf(1,nlast2+loop),nint,iorb)
           if (modei) then
              write(iwr,101)frodo,noop,double,mzero,kc,loop
     *                    ,(iorb(i),i=1,nint)
           else
              write(iwr,101)frodo,noop,double,mzero,mapie(kc),loop
     *                    ,(iorb(mapei(i)),i=1,nint)
           end if
           new = .false.
        else
           if (modei) then
              write(iwr,101)frodo,noop,double,mzero,kc,loop
           else
              write(iwr,101)frodo,noop,double,mzero,mapie(kc),loop
           end if
        end if
901   continue
c... singlet  triplet states
13    if(nmin2)16,999,16
c16    hobbit=dmin1(bilbo,stpri)
16    hobbit=bilbo
      do 1001 loop=nstrt2,nlast2
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,loop)
_ELSE
      isi=ipad(iconf(4,loop))
_ENDIF
      call upack2(iconf(3,loop),lamt,lams)
      new = .true.
      lamss=lams*norbsi(isi)
      lamtt=lamt*norbtr(isi)
      if(lamss)1003,1002,1003
_IF(ga)
1003  call reads(cr(khbill+1),lamss,numcz)
_ELSE
1003  call reads(cr(khbill+1),lamss,numci)
_ENDIF
       koop=khbill
      do 1031 noop=1,lams
      do 1031 jsi=1,nirr
      nj=norbe(jsi)
      if(nj)1008,1031,1008
1008  ksi=mult(jsi,isi)
      mj=nbase(jsi)
      if(ksi-jsi)1031,1007,1006
c... c vector totally symmetric
1007  do 1009 moop=1,nj
      jc=moop+mj
      do 1009 joop=1,moop
      koop=koop+1
      frodo=cr(koop)
      if( dabs(frodo).lt.hobbit)goto 1009
      kc=joop+mj
        if (new) then
           call upack(iconf(1,nlast2+loop),nint,iorb)
           if (modei) then
              write(iwr,101)frodo,noop,single,kc,jc,loop
     *                    ,(iorb(i),i=1,nint)
           else
              write(iwr,101)frodo,noop,single,mapie(kc),mapie(jc),loop
     *                    ,(iorb(mapei(i)),i=1,nint)
           end if
           new = .false.
        else
           if (modei) then
              write(iwr,101)frodo,noop,single,kc,jc,loop
           else
              write(iwr,101)frodo,noop,single,mapie(kc),mapie(jc),loop
           end if
        end if
1009  continue
      goto 1031
c... c vector not totally symmetric
1006   nk=norbe(ksi)
       if(nk)1011,1031,1011
1011  mk=nbase(ksi)
      do 1010 joop=1,nk
      kc=joop+mk
       do 1010 moop=1,nj
      koop=koop+1
      frodo=cr(koop)
       if( dabs(frodo).lt.hobbit)goto 1010
      jc=moop+mj
        if (new) then
           call upack(iconf(1,nlast2+loop),nint,iorb)
           if (modei) then
              write(iwr,101)frodo,noop,single,jc,kc,loop
     *                    ,(iorb(i),i=1,nint)
           else
              write(iwr,101)frodo,noop,single,mapie(jc),mapie(kc),loop
     *                    ,(iorb(mapei(i)),i=1,nint)
           end if
           new = .false.
        else
           if (modei) then
              write(iwr,101)frodo,noop,single,jc,kc,loop
           else
              write(iwr,101)frodo,noop,single,mapie(jc),mapie(kc),loop
           end if
        end if
1010  continue
1031  continue
      if(lamtt)1002,1001,1002
_IF(ga)
1002   call reads(cr(khbill+1),lamtt,numcz)
_ELSE
1002   call reads(cr(khbill+1),lamtt,numci)
_ENDIF
       koop=khbill
      do 1041 noop=1,lamt
      do 1041 jsi=1,nirr
      nj=norbe(jsi)
      ksi=mult(jsi,isi)
      mj=nbase(jsi)
      if(ksi-jsi)1041,1047,1046
c... c vector totally symmetric
1047   if(nj.le.1)goto  1041
      do 1049 moop=2,nj
      jc=moop+mj
      muck=moop-1
      do 1049 joop=1,muck
      koop=koop+1
      frodo=cr(koop)
      if( dabs(frodo).lt.hobbit)goto 1049
      kc=joop+mj
        if (new) then
           call upack(iconf(1,nlast2+loop),nint,iorb)
           if (modei) then
              write(iwr,101)frodo,noop,triple,kc,jc,loop
     *                    ,(iorb(i),i=1,nint)
           else
              write(iwr,101)frodo,noop,triple,mapie(kc),mapie(jc),loop
     *                    ,(iorb(mapei(i)),i=1,nint)
           end if
           new = .false.
        else
           if (modei) then
              write(iwr,101)frodo,noop,triple,kc,jc,loop
           else
              write(iwr,101)frodo,noop,triple,mapie(kc),mapie(jc),loop
           end if
        end if
1049  continue
      goto 1041
c... c vector not totally symmetric
1046  nk=norbe(ksi)
      if(nj*nk)1051,1041,1051
1051  mk=nbase(ksi)
      do 150 joop=1,nk
      kc=joop+mk
      do 150 moop=1,nj
      koop=koop+1
      frodo=cr(koop)
      if( dabs(frodo).lt.hobbit)goto 150
      jc=moop+mj
        if (new) then
           call upack(iconf(1,nlast2+loop),nint,iorb)
           if (modei) then
              write(iwr,101)frodo,noop,triple,jc,kc,loop
     *                    ,(iorb(i),i=1,nint)
           else
              write(iwr,101)frodo,noop,triple,mapie(jc),mapie(kc),loop
     *                    ,(iorb(mapei(i)),i=1,nint)
           end if
           new = .false.
        else
           if (modei) then
              write(iwr,101)frodo,noop,triple,jc,kc,loop
           else
              write(iwr,101)frodo,noop,triple,mapie(jc),mapie(kc),loop
           end if
        end if
150   continue
1041  continue
1001  continue
999   write(iwr,109)sos,nrefs,sosm
109   format(//' sum of squares of vacuum csf ci coefficients=',e21.12/
     +  ' sum of squares of ',i5,'-main ci coefficients=',e21.12/)
c
cjvl   extended ci-vector printer
c
      khbill = khsave
      modei = modeis
      call flushn(iwr)
c
      return
      end
      subroutine nosas(cr,iconf,oprint)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      logical oprint
      character *8 commm,title,direct,com,tit
      character *8 date,time,accno,ajnam,anos,bnos
INCLUDE(common/sizes)
_IF(drf)
INCLUDE(../drf/comdrf/sizesrf)
INCLUDE(../drf/comdrf/drfpar)
INCLUDE(../drf/comdrf/dafil)
INCLUDE(../drf/comdrf/darw)
_ENDIF
INCLUDE(common/drfopt)
      parameter (mxorb1=maxorb+1)
c...
c... space and spin natural orbital generator
c... use for both in and out store modes
c...
      logical lelim,o1e
      dimension cr(*),iconf(5,*),o1e(6)
      common/erg/eigen(nd200+3),h11
      common/moco/ibasso(8),ibavvv(8),lexbas(8)
      common/lsort/ilifx(maxorb),ilifs(maxorb),ilifq(maxorb),
     + lelim(maxorb),comman(maxorb),pottt(2),ncolo,nbas,newb,ncore,
     + mapcie(maxorb),ilifc(maxorb),nsact,mapaie(maxorb),
     + ilifa(maxorb),iqsec
INCLUDE(common/natorb)
INCLUDE(common/file12)
INCLUDE(common/dump3)
      common/scra /ilifd(maxorb),ntrad(maxorb),itran(mxorb3),
     * ctran(mxorb3),itrad(2),deig(maxorb),
     *docc(mxorb1),nbasis,newbas,ncolq,jeig,jocc,jtemp,
     *potnuc,dx,dy,dz,sone(508)
INCLUDE(common/timez)
INCLUDE(common/mapp)
INCLUDE(common/helpr)
INCLUDE(common/symci)
INCLUDE(common/auxh)
INCLUDE(common/orb)
      common/symcon/root2,root2i
INCLUDE(common/cdcryz)
INCLUDE(common/corctl)
      common/spew/rmode(4)
INCLUDE(common/discc)
INCLUDE(common/iofile)
INCLUDE(common/machin)
      common/runlab/commm(19),title(1)
      common/jinfo/date,time,accno,ajnam
      common/junkc/com(19),tit(10)
INCLUDE(common/disktl)
INCLUDE(common/cntrl)
      common/stoctl/khbill
INCLUDE(common/ccntl)
INCLUDE(common/table)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      equivalence (norbq,norbsq(1)),(bilbo,ir)
      equivalence (maa,norbsi(1))
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      equivalence (rmode(4),ninti)
      data direct,anos,bnos/'direct','nos','snos'/
      data m29/29/
      kexit = 0
      irl = 0
      ich = 0
      icl = 0
c
      top=cpulft(1)
      write(iwr,100)top
100   format(//1x,80('=')//
     *' **** natural orbital generator called at' ,f9.2,' secs')
      i26=index(26)
      mtrorb=iky(nint+1)
      nato=mtrorb+m1ia+maa
      natong=-nato
      khhill=khbill+nato
      khia=khhill+mtrorb
      khaa=m1ia+khia
      khvd=khaa+maa
      khd=khvd+nval
      khe=khd+ndoub
      khf=khe+maxpat*maxpat*2
      jc=maxpa(1)*nubsi
      kc=maxpa(2)*nubtr
      mix=jc+kc
      khg=khf+mix
      joke=max(jc,kc)
      khh=khg+mix
      khgg=max12*nubsq+khf
      khi=khg+kc
      khj=maxpa(2)*nubsi+khi
      khk=khj+mubsq
      khl=khk+mubsq
      khm=khl+mubsq
      jc=max(maxpa(4),maxpa(3)*nvmax)
      kc=max(mix,max12*nvmax)
      nuse=max(khm+norbq,khh+joke,khgg+kc,khf+jc)
      call corfai('natorb')
      call vclr(cr(khbill+1),1,nato+nato)
      joke=khia+1
       do 8001 loop=1,nint
       iap(loop)=joke
       joke=norbe(iro(loop))+joke
8001   cr(khhill+iky(loop+1))=2.0d0
      if(nvd)8888,7777,8888
_IF(ga)
8888  call rdedx(cr(khvd+1),nvd,i26,numcz)
_ELSE
8888  call rdedx(cr(khvd+1),nvd,i26,numci)
_ENDIF
7777  call brtsb(index(197))
6666  call getsb(rmode(1),4)
      goto (9999,8000,9000,4,5,6,6,8,9,10,10,12,12,6,15,16),model
c...
c... spin vacuum/vacuum (ii) + (ij)
c...
10    khgill=khbill
      goto 2222
c...
c... vacuum/vacuum   (ij)
c...
8000   khgill=khhill
2222   continue
_IF(cray,convex,t3e,i8,i8drct)
      lamc=iconf(3,ic)
      lamr=iconf(3,ir)
_ELSE
      lamc=ipad(iconf(3,ic))
      lamr=ipad(iconf(3,ir))
_ENDIF
      call getsb(cr(khe+1),lamc*lamr)
_IF(cray,convex,t3e,i8,i8drct)
      call mxmd(cr(khe+1),1,lamr,cr(khvd+iconf(1,ic)+1),1,1,
     *cr(khf+1),1,1,lamr,lamc,1)
      cr(khgill+ninti)=
     * ddot(lamr,cr(khf+1),1,cr(khvd+iconf(1,ir)+1),1)
_ELSE
      call mxmd(cr(khe+1),1,lamr,cr(khvd+ipad(iconf(1,ic))+1),1,1,
     *cr(khf+1),1,1,lamr,lamc,1)
      cr(khgill+ninti)=
     * ddot(lamr,cr(khf+1),1,cr(khvd+ipad(iconf(1,ir))+1),1)
_ENDIF
     * +cr(khgill+ninti)
      goto 6666
c...
c... vacuum/vacuum   (ii)
c...
9000  if(ic.eq.icl)goto 33
_IF(cray,convex,t3e,i8,i8drct)
      jc=khvd+iconf(1,ic)
      loop8=iconf(3,ic)
_ELSE
      jc=khvd+ipad(iconf(1,ic))
      loop8=ipad(iconf(3,ic))
_ENDIF
      frodo=ddot(loop8,cr(jc+1),1,cr(jc+1),1)
333   icl=ic
33    cr(khhill+ninti)=frodo*bilbo+cr(khhill+ninti)
      goto 6666
c...
c... spin doublet/doublet (ii) + (ij)
c...
12    khgill=khbill
      goto 40
c...
c... doublet/doublet   (ij)
c...
4     khgill=khhill
40    continue
_IF(cray,t3e,convex,i8,i8drct)
      lamc=iconf(3,ic)
      lamr=iconf(3,ir)
_ELSE
      lamc=ipad(iconf(3,ic))
      lamr=ipad(iconf(3,ir))
_ENDIF
      call getsb(cr(khe+1),lamc*lamr)
_IF(cray,t3e,convex,i8,i8drct)
      ne=norbe(iconf(4,ir))
      call mxmd(cr(khd+iconf(1,ir)+1),1,ne,cr(khe+1),1,lamr,
     *cr(khf+1),1,ne,ne,lamr,lamc)
      ix1 = khd+iconf(1,ic)+1
_ELSE
      ne=norbe(ipad(iconf(4,ir)))
      call mxmd(cr(khd+ipad(iconf(1,ir))+1),1,ne,cr(khe+1),1,lamr,
     *cr(khf+1),1,ne,ne,lamr,lamc)
      ix1 = khd+ipad(iconf(1,ic))+1
_ENDIF
      cr(khgill+ninti)=
     *ddot(ne*lamc,cr(khf+1),1,cr(ix1),1)
     * +cr(khgill+ninti)
      goto 6666
c...
c... doublet/doublet   (ii)
c...
5     if(ic.eq.icl)goto 33
_IF(cray,t3e,convex,i8,i8drct)
      jc=khd+iconf(1,ic)
      loop8=iconf(3,ic)
      ix1=iconf(4,ic)
_ELSE
      jc=khd+ipad(iconf(1,ic))
      loop8=ipad(iconf(3,ic))
      ix1=ipad(iconf(4,ic))
_ENDIF
      frodo=ddot(loop8*norbe(ix1),cr(jc+1),1,cr(jc+1),1)
      goto 333
c...
c... n-2/n-2   (ij)  and  (ii)
c...
6     if(ic.eq.icl)goto 66
      icl=ic
      call upack2(iconf(3,icl),lamct,lamc)
_IF(cray,t3e,convex,i8,i8drct)
      jsym=iconf(4,icl)
_ELSE
      jsym=ipad(iconf(4,icl))
_ENDIF
      norbs=norbsi(jsym)
      norbt=norbtr(jsym)
      nss=norbs*lamc
      ntt=norbt*lamct
      khff=khf+nss
_IF(ga)
      call upack2(iconf(2,icl),jc2,kc)
_ELSE
_IF(cray,t3e,convex,i8,i8drct)
_IF(cray,t3e)
      kc=and(iconf(2,icl),mms)
_ELSE
      kc=iand(iconf(2,icl),mms)
_ENDIF
_ELSEIF(bits8)
      call dand8(iconf(2,icl),mms,tmp8)
      kc=ipad(tmp8)
_ELSE
      kc=ipad(dand(iconf(2,icl),mms))
_ENDIF
_ENDIF
_IF(ga)
      call rdedx(cr(khf+1),nss,i26+kc,numcz)
      call rdedx(cr(khf+1+nss),ntt,i26+jc2,numcz)
_ELSE
      call rdedx(cr(khf+1),nss+ntt,i26+kc,numci)
_ENDIF
      frodo1=0.0d0
      frodo2=0.0d0
      if(nss.eq.0)go to 6300
      frodo1=ddot(nss,cr(khf+1),1,cr(khf+1),1)
      if(ntt.eq.0)go to 6302
6300  continue
      frodo2=ddot(ntt,cr(khff+1),1,cr(khff+1),1)
6302  frodo=frodo1+frodo2
c... check for diagonal case
66    if(model-7)644,33,61
644   loop=khhill+ninti
63    if(ir.eq.irl)goto 62
c... n-2/n-2   (ij)
      call upack2(iconf(3,ir),lamrt,lamr)
_IF(ga)
      call upack2(iconf(2,ir),kc2,jc)
_ELSE
_IF(cray,t3e,convex,i8,i8drct)
_IF(cray,t3e)
      jc=and(iconf(2,ir),mms)
_ELSE
      jc=iand(iconf(2,ir),mms)
_ENDIF
_ELSEIF(bits8)
      call dand8(iconf(2,ir),mms,tmp8)
      jc=ipad(tmp8)
_ELSE
      jc=ipad(dand(iconf(2,ir),mms))
_ENDIF
_ENDIF
      lss=norbs*lamr
      khgggg=khg+lss
_IF(ga)
      call rdedx(cr(khg+1),lss,i26+jc,numcz)
      call rdedx(cr(khg+1+lss),norbt*lamrt,i26+kc2,numcz)
_ELSE
      call rdedx(cr(khg+1),lss+norbt*lamrt,i26+jc,numci)
_ENDIF
      irl=ir
c... singlet/singlet
62    mix=lamc*lamr
      mixt=lamct*lamrt
      if(mix)68,69,68
68    call getsb(cr(khe+1),mix)
      call mxmd(cr(khg+1),1,norbs,cr(khe+1),1,lamr,
     *cr(khh+1),1,norbs,norbs,lamr,lamc)
      cr(loop)=ddot(nss,cr(khh+1),1,cr(khf+1),1)+cr(loop)
c... triplet/triplet
      if(mixt)69,6666,69
69    call getsb(cr(khe+1),mixt)
      call mxmd(cr(khgggg+1),1,norbt,cr(khe+1),1,lamrt,
     *cr(khh+1),1,norbt,norbt,lamrt,lamct)
 655  continue
       cr(loop)=ddot(ntt,cr(khh+1),1,cr(khff+1),1)+cr(loop)
      goto 6666
c... spin natorb ii or ij case
61    loop=khbill+ninti
      if(ir.ne.ic)goto 63
c... spin natorb (ii) case
      if(lamc)64,65,64
64    call getsb(cr(khe+1),lamc*lamc)
      if(lamc.ne.1)goto 640
      cr(loop)=frodo1*cr(khe+1)+cr(loop)
      goto 650
640   call mxmd(cr(khf+1),1,norbs,cr(khe+1),1,lamc,cr(khh+1),1,norbs,
     *norbs,lamc,lamc)
      cr(loop)=ddot(nss,cr(khh+1),1,cr(khf+1),1)+cr(loop)
650   if(lamct)65,6666,65
65    call getsb(cr(khe+1),lamct*lamct)
      if(lamct.ne.1)goto 660
      cr(loop)=frodo2*cr(khe+1)+cr(loop)
      goto 6666
660   call mxmd(cr(khff+1),1,norbt,cr(khe+1),1,lamct,cr(khh+1),1,norbt,
     *norbt,lamct,lamct)
      goto 655
c...
c... spin doublet/vacuum (ia)
c...
15    mix=natong
      goto 80
c...
c... doublet/vacuum   (ia)
c...
8      mix=0
80     continue
_IF(cray,t3e,convex,i8,i8drct)
      lamc=iconf(3,ic)
      lamr=iconf(3,ir)
_ELSE
      lamc=ipad(iconf(3,ic))
      lamr=ipad(iconf(3,ir))
_ENDIF
      call getsb(cr(khe+1),lamc*lamr)
_IF(cray,t3e,convex,i8,i8drct)
      call mxmd(cr(khe+1),lamr,1,cr(khvd+iconf(1,ir)+1),1,1,
_ELSE
      call mxmd(cr(khe+1),lamr,1,cr(khvd+ipad(iconf(1,ir))+1),1,1,
_ENDIF
     *cr(khf+1),1,1,lamc,lamr,1)
      ne=norbe(iro(ninti))
_IF(cray,t3e,convex,i8,i8drct)
      call mxmb(cr(khd+iconf(1,ic)+1),1,ne,cr(khf+1),1,1,
_ELSE
      call mxmb(cr(khd+ipad(iconf(1,ic))+1),1,ne,cr(khf+1),1,1,
_ENDIF
     *cr(iap(ninti)+mix),1,1,ne,lamc,1)
      goto 6666
c...
c... spin n-2/doublet (ia)
c...
16    mixt=natong
      goto 966
c...
c... n-2/doublet   (ia)
c...
9     mixt=0
966   if(ic.eq.ich)goto 99
      ich=ic
      call upack2(iconf(3,ich),lamct,lamc)
      lamcc=lamc+lamct
_IF(cray,t3e,convex,i8,i8drct)
      jsym=iconf(4,ich)
_ELSE
      jsym=ipad(iconf(4,ich))
_ENDIF
      norbs=norbsi(jsym)
      norbt=norbtr(jsym)
      mix=lamc*norbs+lamct*norbt
_IF(ga)
      call upack2(iconf(2,ich),kc,jc)
_ELSE
_IF(cray,t3e,convex,i8,i8drct)
_IF(cray,t3e)
      jc=and(iconf(2,ich),mms)
_ELSE
      jc=iand(iconf(2,ich),mms)
_ENDIF
_ELSEIF(bits8)
      call dand8(iconf(2,ich),mms,tmp8)
      jc=ipad(tmp8)
_ELSE
      jc=ipad(dand(iconf(2,ich),mms))
_ENDIF
_ENDIF
      if(jsym.ne.1)goto 98
c... n-2 vectors totally symmetric
_IF(ga)
      mixs=lamc*norbs
      mixt=lamct*norbt
      call rdedx(cr(khgg+1),mixs,i26+jc,numcz)
      call rdedx(cr(khgg+1+mixs),mixt,i26+kc,numcz)
_ELSE
      call rdedx(cr(khgg+1),mix,i26+jc,numci)
_ENDIF
      jb=khgg
      kb=khf
      if(lamc)96,97,96
96    do 95 loop=1,lamc
      do 94 isym=1,nirr
      ne=norbe(isym)
      if(ne)93,94,93
93    joke=kb+ibassq(isym,1)
      call square(cr(joke+1),cr(jb+ibassi(isym,1)+1),ne,ne)
      mix=ne+1
      call dscal( ne,root2,cr(joke+1),mix)
94    continue
      jb=jb+norbs
95    kb=kb+norbq
      if(lamct)97,92,97
97    do 91 loop=1,lamct
      do 90 isym=1,nirr
      ne=norbe(isym)
      if(ne.gt.1)
     *call sqtrip(cr(kb+ibassq(isym,1)+1),
     +            cr(jb+ibastr(isym,1)+1),ne)
90    continue
      jb=jb+norbt
91    kb=kb+norbq
92    norbs=norbq
      do 900 loop=1,nirr
900   ibasso(loop)=khf+ibassq(loop,1)
      goto 99
c... n-2 vectors not totally symmetric
_IF(ga)
98    mixs=lamc*norbs
      mixt=lamct*norbt
      call rdedx(cr(khf+1),mixs,i26+jc,numcz)
      call rdedx(cr(khf+1+mixs),mixt,i26+kc,numcz)
_ELSE
98    call rdedx(cr(khf+1),mix,i26+jc,numci)
_ENDIF
      do 1 loop=1,nirr
    1 ibasso(loop)=khf+ibassi(loop,jsym)
c... produce density matrix contributions
99    continue
_IF(cray,t3e,convex,i8,i8drct)
      lamr=iconf(3,ir)
      ksym=iconf(4,ir)
_ELSE
      lamr=ipad(iconf(3,ir))
      ksym=ipad(iconf(4,ir))
_ENDIF
      nf=norbe(ksym)
      lsym=iro(ninti)
      ne=norbe(lsym)
      lamccc=lamcc
      if((jsym+ne).eq.2)lamccc=lamc
      call getsb(cr(khe+1),lamr*lamccc)
_IF(cray,t3e,convex,i8,i8drct)
      call mxmd(cr(khd+iconf(1,ir)+1),1,nf,cr(khe+1),1,lamr,
_ELSE
      call mxmd(cr(khd+ipad(iconf(1,ir))+1),1,nf,cr(khe+1),1,lamr,
_ENDIF
     *cr(khgg+1),1,nf,nf,lamr,lamccc)
       if(lsym.le.ksym)goto 910
      mcol=nf
      mrow=1
      kb=ibasso(ksym)
      goto 909
910   mcol=1
      mrow=ne
      kb=ibasso(lsym)
909   jb=khgg
      jc=iap(ninti)+mixt
      do 905 loop=1,lamccc
      call mxmb(cr(kb+1),mcol,mrow,cr(jb+1),1,nf,cr(jc),1,1,ne,nf,1)
      kb=kb+norbs
905   jb=jb+nf
      goto 6666
9999  joke=nint
      do 1000 loop=1,nirr
      lexbas(loop)=joke
      joke=joke+norbe(loop)
      ibavvv(loop)=ibassq(loop,1)+khm
1000  ibasso(loop)=ibassi(loop,1)+khaa
      if(nmin1)1002,1001,1002
c...
c... doublet/doublet   (ab)
c...
1002  jb=khd
      do 1003 icl=nstrt1,nlast1
_IF(cray,t3e,convex,i8,i8drct)
      lamr=iconf(3,icl)
      lsym=iconf(4,icl)
_ELSE
      lamr=ipad(iconf(3,icl))
      lsym=ipad(iconf(4,icl))
_ENDIF
      ne=norbe(lsym)
      kb=ibasso(lsym)
       if(mspin1)goto 1134
c... spin
      mix=khf
      jc=jb
      kc=kb+natong
      call getsb(cr(khe+1),lamr*lamr)
      call mxmd(cr(jc+1),1,ne,cr(khe+1),1,lamr,cr(mix+1),1,ne,
     *ne,lamr,lamr)
      do 1133 loop=1,lamr
      call mxmtt(cr(mix+1),cr(jc+1),1,1,ne,1,cr(kc+1))
      jc=jc+ne
1133  mix=mix+ne
c... space
1134  do 1003 loop=1,lamr
      call mxmtt(cr(jb+1),cr(jb+1),1,1,ne,1,cr(kb+1))
1003  jb=jb+ne
1001  if(nmin2)1004,1111,1004
c...
c... n-2/n-2   (ab)
c...
c... zeroize the full square form of (ab) spin density matrix
1004  continue
      call vclr(cr(khm+1),1,norbq)
      do 1010 icl=nstrt2,nlast2
      call upack2(iconf(3,icl),lamct,lamc)
      lamr=lamct+lamc
_IF(cray,t3e,convex,i8,i8drct)
      lsym=iconf(4,icl)
_ELSE
      lsym=ipad(iconf(4,icl))
_ENDIF
      norbs=norbsi(lsym)
      norbt=norbtr(lsym)
      lss=lamc*norbs
      ltt=lamct*norbt
      khgggg=khf+lss
_IF(ga)
      call upack2(iconf(2,icl),kc,jc)
_ELSE
_IF(cray,t3e,convex,i8,i8drct)
_IF(cray,t3e)
      jc=and(iconf(2,icl),mms)
_ELSE
      jc=iand(iconf(2,icl),mms)
_ENDIF
_ELSEIF(bits8)
      call dand8(iconf(2,icl),mms,tmp8)
      jc=ipad(tmp8)
_ELSE
      jc=ipad(dand(iconf(2,icl),mms))
_ENDIF
_ENDIF
_IF(ga)
      call rdedx(cr(khf+1),lss,i26+jc,numcz)
      call rdedx(cr(khf+1+lss),ltt,i26+kc,numcz)
_ELSE
      call rdedx(cr(khf+1),ltt+lss,i26+jc,numci)
_ENDIF
      if(lsym.ne.1)goto 1011
c************ n-2 vectors totally symmetric
      if(lamc)1031,1041,1031
c... singlet contribution to space density matrix
1031  do 1033 isym=1,nirr
      ne=norbe(isym)
      if(ne)1034,1033,1034
1034  jb=khf+ibassi(isym,1)
      kc=ibasso(isym)
      joke=ne+1
      do 1032 loop=1,lamc
      call square(cr(khg+1),cr(jb+1),ne,ne)
_IF1()      kb=khg
_IF1()      do 1035 joop=1,ne
_IF1()      cr(kb+1)=cr(kb+1)*root2
_IF1() 1035 kb=kb+joke
c
      call dscal(ne,root2,cr(khg+1),joke)
      call mxmtt(cr(khg+1),cr(khg+1),1,ne,ne,ne,cr(kc+1))
1032  jb=jb+norbs
1033  continue
      if(lamct)1041,1010,1041
c... triplet contribution to space density matrix
1041  do 1043 isym=1,nirr
      ne=norbe(isym)
      if(ne.le.1)goto 1043
      jb=khgggg+ibastr(isym,1)
      kc=ibasso(isym)
      do 1042 loop=1,lamct
      call sqtrip(cr(khg+1),cr(jb+1),ne)
      call mxmtt(cr(khg+1),cr(khg+1),1,ne,ne,ne,cr(kc+1))
1042  jb=jb+norbt
1043  continue
      if(mspin1)goto 1010
c... spin density matrix
      call getsb(cr(khe+1),lamct*lamct)
      call mxmd(cr(khgggg+1),1,norbt,cr(khe+1),1,lamct,cr(khg+1),1,norbt
     *,norbt,lamct,lamct)
      call vclr(cr(khi+1),1,lamct*norbs)
      if(lamc)1401,1402,1401
1401  call getsb(cr(khe+1),lamc*lamct)
      call mxmb(cr(khf+1),1,norbs,cr(khe+1),1,lamc,cr(khi+1),1,norbs,
     *norbs,lamc,lamct)
1402  do 1483 isym=1,nirr
      ne=norbe(isym)
      if(ne.le.1)goto 1483
      nenf=ne*ne
      mix=ibavvv(isym)
      kb=khi+ibassi(isym,1)
      joke=ibastr(isym,1)
      jb=khg+joke
      kc=khgggg+joke
       joke=ne+1
      do 1480 loop=1,lamct
      call square(cr(khk+1),cr(kb+1),ne,ne)
      call dscal(ne,root2,cr(khk+1),joke)
      call sqtrip(cr(khj+1),cr(jb+1),ne)
      call vsub(cr(khk+1),1,cr(khj+1),1,cr(khl+1),1,nenf)
      call sqtrip(cr(khj+1),cr(kc+1),ne)
      call mxmb(cr(khl+1),1,ne,cr(khj+1),1,ne,cr(mix+1),1,ne,ne,ne,ne)
      kb=kb+norbs
      jb=jb+norbt
1480  kc=kc+norbt
1483  continue
      goto 1010
c************ n-2 vectors not totally symmetric
c... space density matrix
1011  do 1083 isym=1,nirr
       jsym=mult(isym,lsym)
      if(jsym.lt.isym)goto 1083
      ne=norbe(isym)
      nf=norbe(jsym)
      if(ne*nf)1084,1083,1084
1084  jc=khf+ibassi(isym,lsym)
      jb=ibasso(isym)
      kb=ibasso(jsym)
      do 1080 loop=1,lamr
      call mxmtt(cr(jc+1),cr(jc+1),1,ne,ne,nf,cr(jb+1))
      call mxmtt(cr(jc+1),cr(jc+1),ne,1,nf,ne,cr(kb+1))
1080  jc=jc+norbs
1083  continue
      if(mspin1)goto 1010
      if(lamct)1300,1010,1300
c... spin density matrix
1300  call getsb(cr(khe+1),lamct*lamct)
      call mxmd(cr(khgggg+1),1,norbs,cr(khe+1),1,lamct,cr(khg+1),
     *1,norbs,norbs,lamct,lamct)
      call vclr(cr(khi+1),1,ltt)
      if(lamc)1301,1302,1301
1301  call getsb(cr(khe+1),lamc*lamct)
      call mxmb(cr(khf+1),1,norbs,cr(khe+1),1,lamc,cr(khi+1),1,norbs,
     *norbs,lamc,lamct)
1302  do 1383 isym=1,nirr
      jsym=mult(isym,lsym)
      if(jsym.lt.isym)goto 1383
      ne=norbe(isym)
      nf=norbe(jsym)
      nenf=ne*nf
      if(nenf)1384,1383,1384
1384  mix=ibavvv(isym)
      mixt=ibavvv(jsym)
      joke=ibassi(isym,lsym)
      jb=khg+joke
      kb=khi+joke
      kc=khgggg+joke
      do 1380 loop=1,lamct
      call vsub(cr(jb +1),1,cr(kb+1),1,cr(khl+1),1,nenf)
      call mxmb(cr(khl+1),1,ne,cr(kc+1),ne,1,cr(mix+1),1,ne,ne,nf,ne)
      call vadd(cr(jb+1),1,cr(kb+1),1,cr(khl+1),1,nenf)
      call mxmb(cr(khl+1),ne,1,cr(kc+1),1,ne,cr(mixt+1),1,nf,nf,ne,nf)
      jb=jb+norbs
      kb=kb+norbs
1380  kc=kc+norbs
1383  continue
1010  continue
      do 1377 isym=1,nirr
      ne=norbe(isym)
      if(ne)1378,1377,1378
1378  call symm1a(cr(ibasso(isym)+natong+1),cr(ibavvv(isym)+1),ne)
1377  continue
c...
c... get density matrix into triangle form
c...
1111  nbsqm=ninnex*ninnex
      khf=khhill+nbsqm
      khfmtr=khf+mtrorb
      ntrorb=iky(ninnex+1)
      nuse=khf+ntrorb
      ntmt=ntrorb-mtrorb
      call corfai('natorb')
 4000 continue
      call dcopy(mtrorb,cr(khhill+kexit+1),1,cr(khf+1),1)
      if(next)2001,2000,2001
2001  continue
      call vclr(cr(khfmtr+1),1,ntmt)
c... translate the (ia) elements
      jb=khia+kexit
      do 2002 loop=1,nint
      isym=iro(loop)
      ne=norbe(isym)
      call dsctr(ne
     +  ,cr(jb+1),iky(lexbas(isym)+1),cr(khf+loop+1))
2002  jb=jb+ne
c... translate the (ab) elements
      kc=khaa+kexit
      do 2005 isym=1,nirr
      ne=norbe(isym)
      if(ne)2006,2005,2006
2006  jc=lexbas(isym)
      joke=jc+khf
      do 2007 loop=1,ne
      call dcopy(loop,cr(kc+1),1,cr(joke+iky(jc+loop)+1),1)
2007  kc=kc+loop
2005  continue
2000  if(kexit)4001,2900,4001
2900  call secget(isec3,1004,jblkk)
      do 30000 i=1,10
30000 tit(i)=title(i)
      call rdedx(comman,mach(15),jblkk,idaf)
      call vclr(docc,1,ncolo)
      nbasis=nbas
      ncolq=ncolo
      newbas=newb
      jeig=-1
      jocc=1
      com(1)=ajnam
      com(2)=date
      com(3)=time
      com(4)=direct
      com(5)=anos
      com(6)=accno
_IF(cray)
      call setsto(maxorb,.false.,lelim)
_ELSE
      call setstl(maxorb,.false.,lelim)
_ENDIF
      call ibasgn(maxorb,khhill,ninnex,ilifs)
      call ibasgn(maxorb,0,nbas,ilifq)
c... diagonalize density matrix (active mo basis)
      call jacobi(cr(khf+1),iky,ninnex,cr,ilifs,ninnex,docc(ncore+1),
     *2,3,1.0d-10)
c...
c...  now get eigen vectors from dumpfile
c...
      len1=lensec(mach(8))
      len2=lensec(mach(9))
c     nbsq=ncolo*nbas
      nbsq=nbas*nbas
      lenv=lensec(nbsq)
      nblkk=1+len1+len2+lenv
      khg=khf+nbsq
      khgg=khg+nbsq
      lenbas=nbas*(nbas+1)/2
      nnnn=ninnex*nbas
      khgggg=khgg+nnnn
      khstvd=khgggg+lenbas
      ilifvv=ilifq(ncore+1)+khf
      nuse=khstvd+lenbas
      call corfai('natorb')
      call secget(iqsec,3,jblkk)
      nav = lenwrd()
      call readi(ilifd,mach(9)*nav,jblkk+len1+1,idaf)
      call reads(cr(khg+1),nbsq,idaf)
      call vclr(cr(ilifvv+1),1,nnnn)
      if(ncore)2700,2777,2700
2700  do 2701 loop=1,ncore
      call dcopy(newb
     +  ,cr(khg+ilifc(loop)+1),1,cr(khf+ilifq(loop)+1),1)
      lelim(mapcie(loop))=.true.
2701  docc(loop)=2.0d0
2777  ncol=ninnex+ncore
      write(iwr,98764)anos
98764 format(//45x,a6,'occupation numbers'/45x,24('-'))
      write(iwr,98765)(docc(loop),loop=1,ncol)
98765 format(/10x,8f14.7)
      do 3000 loop=1,ninnex
      jc=mapaie(mapie(loop))
       lelim(jc)=.true.
3000  ilifx(loop)=ilifq(jc)
       jc=ncol
      do 3010 loop=1,ncolo
      if(lelim(loop))goto 3010
      jc=jc+1
      call dcopy(newb
     +  ,cr(khg+ilifq(loop)+1),1,cr(khf+ilifq(jc)+1),1)
3010  continue
      do 3020 loop=1,ninnex
      call dcopy(newb
     +  ,cr(khg+ilifx(loop)+1),1,cr(khgg+ilifq(loop)+1),1)
3020  continue
      call mxmb(cr(khgg+1),1,nbasis,cr(khhill+1),1,ninnex,
     *cr(ilifvv+1),1,nbasis,newb,ninnex,ninnex)
      if(ispacg)3002,3002,3003
3003  write(iwr,98760)anos,ispacg,ibl3d,yed(idaf)
98760 format(///1x,a6,'stored in section',i4,
     *' of dumpfile starting at block',i7,' of ',a4)
      call secput(ispacg,3,nblkk,iblkk)
      call wrtc(com,m29,iblkk,idaf)
c
c     store total energy
c
      docc(mxorb1) = h11
c
      call wrt3s(deig,mach(8),idaf)
      call wrt3is(ilifd,mach(9)*nav,idaf)
      call wrt3s(cr(khf+1),nbsq,idaf)
c...
c... evaluate dipole,1-electron and kinetic energies
c...
3002  call tdown(cr(khf+1),ilifq,cr(khf+1),ilifq,ncol)
      if(.not.latorb. and. oprint) then
       write(iwr,98761)anos
98761  format(//1x,80('-')/
     * //45x,a6,'(a.o. basis)'/45x,18('-'))
c      call writem(cr(khf+1),ilifq,newb,ncol)
       call prev(cr(khf+1),docc,ncol,newb,newb)
      endif
3001  lenbas=nbas*(nbas+1)/2
      call vclr(cr(khgggg+1),1,lenbas)
      do 702 loop=1,ncol
      mix=ilifq(loop)+khf
      kb=khgggg
      frodo=docc(loop)*2.0d0
      do 702 nss=1,nbas
      call daxpy(nss
     +  ,frodo*cr(mix+nss),cr(mix+1),1,cr(kb+1),1)
702   kb=kb+nss
       kb=khgggg
_IF(drf)
      if (field .ne. ' ') then
        do loop = 1, lenbas
          cr(khgggg+loop) = cr(khgggg+loop)*0.5d00
        enddo
        call dawrit(idafh,ioda,cr(khgggg+1),lenbas,16,nav)
        do loop = 1, lenbas
          cr(khgggg+loop) = cr(khgggg+loop)*2.0d00
        enddo
       endif
_ENDIF
       do 2 loop=1,nbas
       cr(kb+1)=cr(kb+1)*0.5d0
    2  kb=kb+loop+1
c
      if(kexit)781,788,781

788   continue
c
c     call getmat to establish block pointers and restore potnuc etc
c
_IF(cray)
      call setsto(6,.false.,o1e)
_ELSE
      call setstl(6,.false.,o1e)
_ENDIF
      call getmat(cr(khstvd+1),cr(khstvd+1),cr(khstvd+1),
     +            cr(khstvd+1),cr(khstvd+1),cr(khstvd+1),
     +            potnuc,nbas,o1e,isect(492))
c
c     loop over t, t+v, x, y, z
c
      call rdedx(cr(khstvd+1),lenbas,ibl3t,idaf)
      top = ddot(lenbas,cr(khstvd+1),1,cr(khgggg+1),1)
      call rdedx(cr(khstvd+1),lenbas,ibl3f,idaf)
      frodo = ddot(lenbas,cr(khstvd+1),1,cr(khgggg+1),1)
      call rdedx(cr(khstvd+1),lenbas,ibl3x,idaf)
      ex = ddot(lenbas,cr(khstvd+1),1,cr(khgggg+1),1)
      call rdedx(cr(khstvd+1),lenbas,ibl3y,idaf)
      ey = ddot(lenbas,cr(khstvd+1),1,cr(khgggg+1),1)
      call rdedx(cr(khstvd+1),lenbas,ibl3z,idaf)
      ez = ddot(lenbas,cr(khstvd+1),1,cr(khgggg+1),1)
      tx=dx+ex
      ty=dy+ey
      tz=dz+ez
      tsq= dsqrt(tx*tx+ty*ty+tz*tz)
      write(iwr,98769)h11,potnuc,frodo,top,dx,ex,tx,dy,ey,ty,dz,ez,tz
     * ,tsq
98769 format(//
     + 1x,'            total energy',f17.10/
     + 1x,'nuclear repulsion energy',f17.10/
     + 1x,'       1-electron energy',f17.10/
     + 1x,'          kinetic energy',f17.10//
     +17x,'dipole moments'/17x,14('=')//
     + 11x,'nuclear',6x,'electronic',11x,'total'/
     + 11x,7('-'),6x,10('-'),11x,5('-')/
     *' x',3f16.7/' y',3f16.7/' z',3f16.7//
     *' total dipole moment (a.u.)=',f16.7)
      if(mspin1)go to 4101
c...
c... generate spin natural orbitals
c...
      kexit=natong
      go to 4000
4001  continue
      call vclr(docc,1,ncolo)
      call jacobi(cr(khf+1),iky,ninnex,cr,ilifs,ninnex,docc,2,3,1.0d-10)
      write(iwr,98764)bnos
      write(iwr,98765)(docc(loop),loop=1,ninnex)
      call vclr(cr(khf+1),1,nnnn)
      call mxmb(cr(khgg+1),1,nbasis,cr(khhill+1),1,ninnex,
     *cr(khf+1),1,nbasis,newb,ninnex,ninnex)
      if(ncore)4203,4202,4203
4203  do 4210 loop=1,ncore
      jb=mapcie(loop)
      ix1 = khg+ilifq(jb)+1
      ix2 = khf+ilifq(ninnex+loop)+1
      call dcopy(newb,cr(ix1),1,cr(ix2),1)
 4210 continue
4202  if(isping)4201,4211,4201
 4201 com(5)=bnos
      write(iwr,98760)bnos,isping,ibl3d,yed(idaf)
      call secput(isping,3,nblkk,iblkk)
      call wrtc(com,m29,iblkk,idaf)
      call wrt3s(deig,mach(8),idaf)
      call wrt3is(ilifd,mach(9)*nav,idaf)
      call wrt3s(cr(khf+1),nbsq,idaf)
4211  if(.not.latorb. and. oprint) then
       write(iwr,98761)bnos
       call tdown(cr(khf+1),ilifq,cr(khf+1),ilifq,ninnex)
c      call writem(cr(khf+1),ilifq,newb,ninnex)
       call prev(cr(khf+1),docc,ninnex,newb,newb)
      endif
4101  if(ispind)goto 4200
      goto 3001
781   continue
_IF()
c    isect(497) should contain the fermi contact integrals
c    originally from atmol integw. At present de-activate
c    this code .. will be computed anyway if PROPERTY ATOMS
c    is included
781   call secget(isect(497),2000,iblkk)
      call rdedx(cr(khf+1),lenbas,iblkk,idaf)
      bilbo=ddot(lenbas,cr(khf+1),1,cr(khgggg+1),1)
      write(iwr,787)bilbo
787   format(//' spin density=',e21.11)
_ENDIF
4200  return
      end
_EXTRACT(sym12,hp800)
      subroutine sym12
c
c...  symmetry-checking of 1- and 2-electron integrals
c
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(convex,i8drct)
      integer *8 iior,ijkln,ior4
      integer *8 itemp8
_ELSE
      integer  iior,ijkln,ior4
_ENDIF
      logical ijk256
      integer popcnt
_IF(extleadz)
      external leadz
_ENDIF
_IF(extpopcnt)
      external popcnt
_ENDIF
INCLUDE(common/sizes)
      common/blkin/value(510),mmword
INCLUDE(common/machin)
INCLUDE(common/helpr)
INCLUDE(common/orb)
INCLUDE(common/discc)
INCLUDE(common/iofile)
INCLUDE(common/cdcryi)
INCLUDE(common/cic)
      logical lsymd,lsadap
      common/symchk/crtsym,excit,nrep,lsymd,lsadap
INCLUDE(common/file12)
INCLUDE(common/mapp)
_IF(ibm,vax)
      common/craypk/i205(340),j205(340),k205(340),l205(340)
_ELSE
      common/craypk/i205(1360)
_ENDIF
INCLUDE(common/table)
      common/expanc/mapiro(nd200)
INCLUDE(common/prints)
INCLUDE(common/atmblk)
c..
       parameter (nijkln=303,n32map=32)
       dimension ijkln(nijkln),maplab(n32map)
c..
c..      next array's to speed popcnt's / compare's up (for non cray)
c..      they are only used where it really matters
c..      and dont work if nrep still > 8
c..
      dimension ipopcn(256),ijk256(256)
c..
_IF(linux)
      external fget
_ENDIF
      data ipopcn/0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
     1           1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
     2           1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
     3           2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
     4           1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
     5           2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
     6           2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
     7           3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
     8           1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
     9           2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
     a           2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
     b           3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
     c           2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
     d           3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
     e           3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
     f           4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8/
      nin=0
      nun = 0
      top=cpulft(1)
c..
_IF(cray)
      call setsto(256,.false.,ijk256)
_ELSE
      call setstl(256,.false.,ijk256)
_ENDIF
c..      if next do-loop causes problems (overflow) lower n32map
c..      alternatively n32map may be as high as 47 on crays etc
      ii = 1
      do 10 i=1,n32map
         maplab(i) = ii
         ii = ii + ii
10    continue
c...
c... test 1-elec integrals
c...
      call setsto(ninnex,0,iro)
      call secini(ibl3d,idaf)
      call secget(isec3,1004,jblkk)
      if(.not.oprint(28)) write(iwr,1100)top,isec3,ibl3d,yed(idaf)
1100  format(//,
     *' **** symmetry check called at',f13.3,' secs'
     *//' transformed 1-electron integrals restored from section',i4,
     *' of dumpfile starting at block',i8,' of ',a4)
      jblkk=jblkk+lensec(mach(15))
      call setsto(1360,0,i205)
      call search(jblkk,idaf)
      call fget(value,kword,idaf)
      joop = iky(ninnex+1)
      do 101 iii=1,nd200
      i=mapei(iii)
      do 101 jjj=1,iii
      j=mapei(jjj)
      nin=nin+1
      if(nin.le.kword)goto 102
      call fget(value,kword,idaf)
      nin=1
102   if (i.gt.ninnex.or.j.gt.ninnex) go to 101
      joop=joop-1
      if( dabs(value(nin)).lt.crtsym)goto 1044
      iroi=iro(i)
      iroj=iro(j)
      if (iroj.eq.0) then
         nrep=nrep+1
         iroj=nrep
      end if
      if (iroi.eq.0) then
         nrep=nrep+1
         iroi=nrep
      end if
      iro(i)=iroi
      iro(j)=iroj
      if (iroi.ne.iroj)
     1 call comsym(iroi,iroj,ipopcn,ijk256,ijkln,maplab,nun)
1044  if(joop.eq.0) go to 10102
101   continue
10102 if(nrep.eq.1)goto 776
c..
      if (nrep.gt.n32map) then
          write(iwr,10103) nrep,n32map
10103   format(/' *** # reps (',i3,' ) > ',i3,' / arbitrary reduction ')
10104     irf = nrep
          call comsym(irf,irf-1,ipopcn,ijk256,ijkln,maplab,nun)
c..       (nrep in /symchk/ is reduced in comsym)
          if (nrep.gt.n32map) go to 10104
       end if
c..
c
c...  get 2-electron symmetry information
c
_IF(cray)
      call gather(ninnex,mapiro,maplab,iro)
_ELSE
      call igthr(ninnex,maplab,mapiro,iro)
_ENDIF
      do 20113 ioopn=1,nsect2
      ibl=iblk2(ioopn)
      lbl=ibl-lblk2(ioopn)
      if(lbl)66889,20113,66889
66889 iunit=notap2(ioopn)
      call search(ibl,iunit)
66    call fget(value,joop,iunit)
      if(joop)20110,20113,20110
20110 if(modei)goto 2077
      call upak8w(value(num2e+1),i205,mapei)
      goto 2088
2077  continue
_IF(ibm,vax)
      call upak8v(value(num2e+1),i205)
_ELSE
      call unpack(value(num2e+1),lab816,i205,numlab)
_ENDIF
2088  continue
_IFN1(iv)      int4=1
      do 33 iwword=1,mmword
      if ( dabs(value(iwword)).lt.crtsym) go to 33
_IF(ibm,vax)
      i=i205(iwword)
      j=j205(iwword)
      k=k205(iwword)
      l=l205(iwword)
_ELSE
      i=i205(int4  )
      j=i205(int4+1)
      k=i205(int4+2)
      l=i205(int4+3)
_ENDIF
      if(i.gt.ninnex.or.j.gt.ninnex.or.k.gt.ninnex.or.l.gt.ninnex)
     1 goto 33
      i=mapiro(i)
      j=mapiro(j)
      k=mapiro(k)
      l=mapiro(l)
_IF(cray,t3e,convex,t3d,t3e)
      ior4=or(i,or(j,or(k,l)))
_ELSEIF(i8,i8drct)
      ior4=IOR64(i,IOR64(j,IOR64(k,l)))
_ELSE
      ior4=IOR32(i,IOR32(j,IOR32(k,l)))
_ENDIF
      if (nrep.le.8) then
         nior=ipopcn(ior4+1)
      else
_IF(cray,t3e,convex,i8,i8drct)
         nior = popcnt(ior4)
_ELSEIF(bits8)
         call pad8(ior4,tmp8)
         nior = popcnt(tmp8)
_ELSE
         nior = popcnt(pad(ior4))
_ENDIF
      end if
      goto (33,200,300,400),nior
c...   nior=2
200   continue
_IF(cray,t3e,convex,t3d,t3e)
      iior=xor(i,xor(j,xor(k,l)))
_ELSEIF(i8,i8drct)
      iior=IXOR64(i,IXOR64(j,IXOR64(k,l)))
_ELSE
      iior=IXOR32(i,IXOR32(j,IXOR32(k,l)))
_ENDIF
      if(iior.ne.ior4     )go to 33
c... (ii/ij) case
_IF(cray,t3e,convex,i8,i8drct)
201   ir1 = 64 - leadz(ior4)
_IF(cray,t3e)
      ir2 = 64 - leadz(xor(ior4,maplab(ir1)))
_ELSE
      itemp8 = maplab(ir1)
      ir2 = 64 - leadz(ieor(ior4,itemp8))
_ENDIF
_ELSEIF(bits8)
201   call pad8(ior4,tmp8)
      ir1 = 64-leadz(tmp8)
      itmp = IXOR32(ior4,maplab(ir1))
      call pad8(itmp,tmp8)
      ir2 = 64 - leadz(tmp8)
_ELSE
201   ir1 = 64-leadz(pad(ior4))
      ir2 = 64 - leadz(pad(IXOR32(ior4,maplab(ir1))))
_ENDIF
      call comsym(ir1,ir2,ipopcn,ijk256,ijkln,maplab,nun)
_IF(cray)
      call gather(ninnex,mapiro,maplab,iro)
_ELSE
      call igthr(ninnex,maplab,mapiro,iro)
_ENDIF
c..
c..      check if more reps have become identical
c..
      do 260 ii=1,nun
_IF(cray,t3e,convex,i8,i8drct)
         if (popcnt(ijkln(ii)).eq.2) then
_ELSEIF(bits8)
         call pad8(ijkln(ii),tmp8)
         if (popcnt(tmp8).eq.2) then
_ELSE
         if (popcnt(pad(ijkln(ii))).eq.2) then
_ENDIF
         ior4 = ijkln(ii)
         go to 201
         end if
260   continue
      goto 33
c... nior=3  (ii/jk)
300   continue
_IF(cray,t3e,convex,t3d,t3e)
      ior4=xor(i,xor(j,xor(k,l)))
_ELSEIF(i8,i8drct)
      ior4=IXOR64(i,IXOR64(j,IXOR64(k,l)))
_ELSE
      ior4=IXOR32(i,IXOR32(j,IXOR32(k,l)))
_ENDIF
      goto 201
c... nior=4    new (ij/kl) ??? - first check for probable equivalence
400   continue
      if (nrep.gt.8) then
c..      more then 8 reps still / do it in the slow way
         do 401 ii = 1,nun
401      if (ijkln(ii).eq.ior4) go to 33
      else
         if(ijk256(ior4+1))go to 33
         ijk256(ior4+1)=.true.
      end if
      nun=nun+1
      if (nun.gt.nijkln) call caserr(' symdiv nun overflow see jvl ')
      ijkln(nun)=ior4
_IF(ibm,vax)
   33 continue
_ELSE
   33 int4=int4+4
_ENDIF
      lbl=lbl+1
      if(lbl)66,20113,66
20113 continue
776   nirr=nrep
      if(nrep.le.2)goto 77777
      if(nrep.ge.4)goto 777
      nirr=4
      goto 77777
c...
c... set - up multiplication table
c...
777   continue
c..
c..      first check if nrep must be reduced further
c..
      if (nrep.gt.8) then
         write(iwr,800) nrep
800      format(' ** # symdet reps too high (',i4,' ) / reduce')
         ir1 = nrep
         ir2 = nrep - 1
810      call comsym(ir1,ir2,ipopcn,ijk256,ijkln,maplab,nun)
c..        no need to update mapiro / check ijkln 's
         do 820 ii=1,nun
_IF(cray,t3e,convex,i8,i8drct)
            if (popcnt(ijkln(ii)).eq.2) then
               ior4 = ijkln(ii)
               ir1 = 64 - leadz(ior4)
_IF(cray,t3e)
               ir2 = 64 - leadz(xor(ior4,maplab(ir1)))
_ELSE
               itemp8 = maplab(ir1)
               ir2 = 64 - leadz(ieor(ior4,itemp8))
_ENDIF
_ELSEIF(bits8)
            call pad8(ijkln(ii),tmp8)
            if(popcnt(tmp8).eq.2) then
               ior4 = ijkln(ii)
               call pad8(ior4,tmp8)
               ir1 = 64 - leadz(tmp8)
               itmp  = IXOR32(ior4,maplab(ir1))
               call pad8(itmp,tmp8)
               ir2 = 64 - leadz(tmp8)
_ELSE
            if (popcnt(pad(ijkln(ii))).eq.2) then
               ior4 = ijkln(ii)
               ir1 = 64 - leadz(pad(ior4))
               ir2 = 64 - leadz(pad(IXOR32(ior4,maplab(ir1))))
_ENDIF
               go to 810
            end if
820      continue
c...    see if ok now
         go to 776
      end if
c..
      call setsto(54,0,mult(3,2))
      do 444 i=3,8
      mult(i,i)=1
444   mult(1,i)=i
447    call fillin (ijkln,maplab,nun)
c...   look for unassigned elements of multiplication table
      if (nrep.gt.4) nrep=8
      nirr=nrep
      do 448 i=2,nirr
      k=i-1
      do 448 j=1,k
       if (mult(j,i).eq.0) then
c...    assign next available representation
c...    the rare cases that we have d2h with 4 irreps are ignored
         newrep = 1
445      newrep = newrep + 1
         if (newrep.gt.nrep) call caserr(' symmetry assign error')
         do 446 ii=1,nrep
c...       see if this one was used
          if (newrep.eq.mult(ii,i).or.newrep.eq.mult(i,ii).or.
     1        newrep.eq.mult(ii,j).or.newrep.eq.mult(j,ii)) go to 445
446      continue
         mult(j,i)=newrep
         mult(i,j)=newrep
         go to 447
       end if
448   continue
c
c...  assign symmetry labels
c
77777 sname=anumt(nirr)
      do 77778 i=1,nirr
77778 rname(i)=anumt(i)
      if(.not.oprint(28)) then
       write(iwr,557)crtsym
557    format(/' integral threshold for symdet=',e15.7//
     * 7x,'multiplication table'/
     * 7x,'====================')
       do 558 j=1,nirr
558    write(iwr,559)(mult(i,j),i=1,nirr)
559    format(5x,8i3)
      endif
      return
      end
_ENDEXTRACT
      subroutine comsym(ir1,ir2,ipopcn,ijk256,ijkln,maplab,nun)
c
c...  symmetry ir1 and ir2 are identical really
c...  change the highest to the lowest one and adapt ijkln combi's
c
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(convex,i8drct)
      integer *8 ijkln,ijkl,ijkl2
      integer *8 irtl,irfl
_ELSE
      integer  ijkln,ijkl,ijkl2
_ENDIF
      integer popcnt
_IF(extpopcnt)
      external popcnt
_ENDIF
      logical ijk256
c
      dimension ipopcn(256),ijk256(256),ijkln(*),maplab(*)
c
INCLUDE(common/sizes)
INCLUDE(common/helpr)
INCLUDE(common/orb)
INCLUDE(common/table)
INCLUDE(common/cic)
      logical lsymd,lsadap
      common/symchk/crtsym,excit,nrep,lsymd,lsadap
      irf = max(ir1,ir2)
      irt = min(ir1,ir2)
      nrep=nrep-1
c..      one can come here again to get numbers consequtive
10    do 20 ii=1,ninnex
20    if (iro(ii).eq.irf) iro(ii) = irt
c
c...  adapt ijkln if necessary
c
      if (nun.eq.0) go to 100
      irfl = maplab(irf)
      irtl = maplab(irt)
      lul=0
_IF(cray)
      call setsto(256,.false.,ijk256)
_ELSE
      call setstl(256,.false.,ijk256)
_ENDIF
      do 60 loop=1,nun
         ijkl=ijkln(loop)
_IF(cray,t3e,convex,t3d,t3e)
         ijkl2 = xor(ijkl,irfl)
_ELSEIF(i8,i8drct)
         ijkl2 = IXOR64(ijkl,irfl)
_ELSE
         ijkl2 = IXOR32(ijkl,irfl)
_ENDIF
c...     if irf was there there should be 1 one less so ..
_IF(cray,t3e,convex,i8,i8drct)
         if (popcnt(ijkl2).lt.popcnt(ijkl)) then
_ELSEIF(bits8)
         call pad8(ijkl2,tmp1)
         call pad8(ijkl,tmp2)
         if(popcnt(tmp1).lt.popcnt(tmp2)) then
_ELSE
         if (popcnt(pad(ijkl2)).lt.popcnt(pad(ijkl))) then
_ENDIF
c...   add  the irt  (using xor !!!) so  ....
c..    if the rep was already there it is annihilated leaving
c..    two reps , which then must be identical
_IF(cray,t3e,convex,t3d,t3e)
            ijkl = xor(ijkl2,irtl)
_ELSEIF(i8,i8drct)
            ijkl = IXOR64(ijkl2,irtl)
_ELSE
            ijkl = IXOR32(ijkl2,irtl)
_ENDIF
         end if
c...       add ijkl to the list  / check on # reps outside
      if (nrep.gt.8) then
_IF(cray,t3e,convex,i8,i8drct)
         if (popcnt(ijkl).lt.2) go to 60
_ELSEIF(bits8)
         call pad8(ijkl,tmp8)
         if (popcnt(tmp8).lt.2) go to 60
_ELSE
         if (popcnt(pad(ijkl)).lt.2) go to 60
_ENDIF
         do 59 ii=1,lul
59       if (ijkln(ii).eq.ijkl) go to 60
      else
         if (ijk256(ijkl+1).or.ipopcn(ijkl+1).lt.2) go to 60
         ijk256(ijkl+1)=.true.
      end if
      lul=lul+1
      ijkln(lul)=ijkl
60    continue
      nun=lul
c
c..      see if all reps are consecutive , if so return
c
100       if (irf.gt.nrep) go to 200
          irt = irf
          irf = irf + 1
          go to 10
c
200   return
      end
      subroutine fillin(ijkln,maplab,nun)
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(convex,i8drct)
      integer *8 ijkl,ijkln,itemp8
_ELSE
      integer  ijkln,ijkl
_ENDIF
      logical logica
_IF(extleadz)
      external leadz
_ENDIF
      dimension ijkln(nun),maplab(*)
c...
c... complete as much of multiplication table as possible from  ijkln
c...
      if (nun.eq.0) return
c...
1     logica=.false.
      do 2 loop=1,nun
      ijkl=ijkln(loop)
_IF(cray,t3e,convex,i8,i8drct)
      i = 64 - leadz(ijkl)
_IF(cray,t3e)
      ijkl = xor(ijkl,maplab(i))
_ELSE
      itemp8 = maplab(i)
      ijkl = ieor(ijkl,itemp8)
_ENDIF
      j = 64 - leadz(ijkl)
_IF(cray,t3e)
      ijkl = xor(ijkl,maplab(j))
_ELSE
      itemp8 = maplab(j)
      ijkl = ieor(ijkl,itemp8)
_ENDIF
      k = 64 - leadz(ijkl)
_IF(cray,t3e)
      ijkl = xor(ijkl,maplab(k))
_ELSE
      itemp8 = maplab(k)
      ijkl = ieor(ijkl,itemp8)
_ENDIF
      l = 64 - leadz(ijkl)
_ELSEIF(bits8)
      call pad8(ijkl,tmp8)
      i = 64 - leadz(tmp8)
      ijkl = IXOR32(ijkl,maplab(i))
      call pad8(ijkl,tmp8)
      j = 64 - leadz(tmp8)
      ijkl = IXOR32(ijkl,maplab(j))
      call pad8(ijkl,tmp8)
      k = 64 - leadz(tmp8)
      ijkl = IXOR32(ijkl,maplab(k))
      call pad8(ijkl,tmp8)
      l = 64 - leadz(tmp8)
_ELSE
      i = 64 - leadz(pad(ijkl))
_IF(cray,t3e,convex,t3d,t3e,i8drct)
      ijkl = xor(ijkl,maplab(i))
      j = 64 - leadz(pad(ijkl))
      ijkl = xor(ijkl,maplab(j))
      k = 64 - leadz(pad(ijkl))
      ijkl = xor(ijkl,maplab(k))
      l = 64 - leadz(pad(ijkl))
_ELSEIF(i8)
      ijkl = IXOR64(ijkl,maplab(i))
      j = 64 - leadz(pad(ijkl))
      ijkl = IXOR64(ijkl,maplab(j))
      k = 64 - leadz(pad(ijkl))
      ijkl = IXOR64(ijkl,maplab(k))
      l = 64 - leadz(pad(ijkl))
_ELSE
      ijkl = IXOR32(ijkl,maplab(i))
      j = 64 - leadz(pad(ijkl))
      ijkl = IXOR32(ijkl,maplab(j))
      k = 64 - leadz(pad(ijkl))
      ijkl = IXOR32(ijkl,maplab(k))
      l = 64 - leadz(pad(ijkl))
_ENDIF
_ENDIF
      call mulfix(i,j,k,l,logica)
      call mulfix(i,k,j,l,logica)
  2   call mulfix(i,l,j,k,logica)
      if(logica)goto 1
c...
      return
      end
      subroutine mulfix(i,j,k,l,logica)
c
      implicit REAL  (a-h,o-z),integer (i-n)
      logical logica
INCLUDE(common/table)
c
      if (mult(i,j).eq.0) then
         if (mult(k,l).eq.0) then
c..          both 0
            return
         else
c...       ij  0 , kl ne 0 => ij = kl
            mult(i,j) = mult(k,l)
            mult(j,i) = mult(k,l)
         end if
      else
         if (mult(k,l).eq.0) then
c..         ij ne 0  kl eq 0  => kl = ij
            mult(k,l) = mult(i,j)
            mult(l,k) = mult(i,j)
         else
c..          both ne 0
            if (mult(i,j).ne.mult(k,l))
     * call caserr('error in automulfix symmetry assigment')
            return
         end if
      end if
c
      logica=.true.
c
      return
      end
      subroutine symsap(cr)
c
      implicit REAL  (a-h,o-z),integer  (i-n)
c...
c...   ** gamess ** pick up symmetries from gamess adapt/orbitals
c...
INCLUDE(common/sizes)
       parameter (mxorb1=maxorb+1)
      dimension cr(*)
c...    common/craypk/ used to receive section 490
        common/craypk/ nira,mula(8,8),isymao(maxorb),isymmo(maxorb)
c...
INCLUDE(common/table)
INCLUDE(common/cic)
        logical lsymd,lsadap
        common/symchk/ crtsym,excit,nrep,lsymd,lsadap
INCLUDE(common/machin)
INCLUDE(common/helpr)
INCLUDE(common/orb)
INCLUDE(common/mapp)
INCLUDE(common/file12)
INCLUDE(common/iofile)
INCLUDE(common/prints)
c..     next common blocks for reading of vectors
      logical otrad
      common/lsort/ilifx(maxorb),ilifs(maxorb),ilifq(maxorb),
     * lelim(maxorb),
     *comman(maxorb),pottt(2),ncolo,nbas,newb,ncore,
     *mapcie(maxorb),ilifc(maxorb),nsact,mapaie(maxorb),ilifa(maxorb)
     *,iqsec
      common/scra  /ilifd(maxorb),ntrad(maxorb),itran(mxorb3),
     * ctran(mxorb3),otrad,itrad,deig(maxorb),
     *docc(mxorb1),nbasis,newbas,ncolq,jeig,jocc,jtemp,
     *potnuc,dx,dy,dz,sone(508)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
c
c...     get vector info
c
      call secget(isec3,1004,jblkk)
      call rdedx(comman,mach(15),jblkk,idaf)
c...
c...  now get eigen vectors from dumpfile
c...
      len1=lensec(mach(8))
      nbsq=nbas*nbas
      call secget(iqsec,3,jblkk)
      nav = lenwrd()
      call readi(ilifd,mach(9)*nav,jblkk+len1+1,idaf)
      iblkrl = iposun(idaf)
      call rdedx(cr,nbsq,iblkrl,idaf)
c
c...     if vectors were not symmetry adapted switch to symdet
c
      if (otrad) then
         lsadap = .false.
         lsymd = .true.
         write(iwr,560)
560      format(/' ** no adapt info .. switch to symdet **')
         return
      end if
c..
        call secget(isect(490),51,iblock)
        call readi(nira,mach(13)*nav,iblock,idaf)
c..
        do 5 i=1,8
        do 5 j=1,8
5       mult(i,j) = mula(i,j)
        nirr = nira
cjvl     if nirr not 1,2,4,8, make it so 
        if (nirr.eq.3) nirr = 4
        if (nirr.eq.5.or.nirr.eq.6.or.nirr.eq.7) nirr = 8
cjvl
        nrep = nirr
c..   assign labels
        sname = anumt(nirr)
        do 6 i=1,nirr
6       rname(i) = anumt(i)
c..
c..        build iro  (reordered)
c..
      aax = 0.0d0
      iix = 0
c
      do 20 loop=1,ninnex
         jc=mapaie(mapie(loop))
         ilif = (jc-1)*nbas
c...     find position of absolute max. element
         ix = 1
         ax = dabs(cr(ilif+1))
         do 10 i=2,newb
            if (dabs(cr(ilif+i)).gt.ax) then
               ax = dabs(cr(ilif+i))
               ix = i
            end if
10       continue
c
c...     symmetry of mo is symmetry of ao with largest coef
c
         iro(loop) = isymao(ix)
c
c..      check if orbital pure symmetry (i.e. find worst case)
c
         do 15 i=1,newb
            if (isymao(i).ne.iro(loop)) then
               if (dabs(cr(ilif+i)).gt.aax) then
                  aax = dabs(cr(ilif+i))
                  iix = loop
               end if
            end if
15       continue
20    continue
c
c..
      if(.not.oprint(28)) then
       write(iwr,559)
559    format(/'   ** symmetry taken from adapt ** ',//
     1       ,7x,'multiplication table'/,
     2        7x,'====================')
       do 558 j=1,nirr
558    write(iwr,557) (mult(i,j),i=1,nirr)
557    format(5x,8i3)
      endif
      if (nirr.ne.nira) write(iwr,556) nira,nirr
556   format(/' ** # representations extended from',i6,' to',i6,' ** ')
c..
       if (aax.gt.1.0d-5) then
         write(iwr,562) aax,iix
562   format(/' *** largest deviation from pure symmetry of ',e12.5,
     *       ' detected in orbital ',i4,' ***')
      end if
c..
      return
      end
      subroutine inreor(conf,nword,nref,iwr)
c
c...  issue a internal reorder to get the symmetries aligned
c...  if a detsym directive is executed
c
c...  reorder the occupations in the reference 'conf'igurations
c...  and the orbital restrictions 'mestri' and the 'eigen'-values
c...  as well
c
c...  common/expanc/(iper,dum,dumm) is used for scratch
c
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer conf
_ELSEIF(convex,i8drct)
      integer *8 conf
_ELSE
      REAL conf
_ENDIF
INCLUDE(common/hold)
      dimension nestri(4)
      equivalence (nestri(1),nesvac)
INCLUDE(common/sizes)
INCLUDE(common/helpr)
INCLUDE(common/orb)
INCLUDE(common/symci)
INCLUDE(common/mapp)
INCLUDE(common/table)
INCLUDE(common/cic)
      common/expanc/ idum(nd200),iper(nd200),idumm(100)
INCLUDE(common/prints)
c
      dimension conf(nword,nref)
c
      if(nirr.eq.0) go to 9999
      modei = .false.
c
c...  reorder the reordering / keep track of permutation in iper
c
      mcount = 0
      do  50 i=1,nirr
      mcx=mcount
         do 40 loop=1,nint
            if (iro(loop).ne.i) go to 40
               mcount = mcount + 1
               iper(loop) = mcount
40       continue
50    norbi(i)=mcount-mcx
c
      do 100 i=1,nirr
      mcx=mcount
      if(nintp1.gt.ninnex)go to 100
         do 90 loop=nintp1,ninnex
            if (iro(loop).ne.i) go to 90
               mcount = mcount + 1
               iper(loop) = mcount
90       continue
100   norbe(i)  = mcount-mcx
c
c...  now use iper to reorder mestri and conf
c
c...    mestri
      do 370 i=1,3
         nn = nestri(i)
         if (nn.le.0) go to 370
         do 340 loop=1,nn
            call permco(mestri(1,loop,i),nint,iper,idum,idumm)
340      continue
370   continue
c
c...    conf
c
      do 450 i=1,nref
         call permco(conf(1,i),nint,iper,idum,idumm)
450   continue
c
c...  put mapei, mapie and iro right
c
      call setsto(nd200,nd200,mapei)
_IF(cray)
      call scatter(ninnex,mapei,mapie,iper)
_ELSE
      call isctr(ninnex,iper,mapie,mapei)
_ENDIF
      call ibasgn(nd200,1,1,idum)
_IF(cray)
      call scatter(nd200,mapie,mapei,idum)
      call scatter(ninnex,idum,iper,iro)
_ELSE
      call isctr(nd200,idum,mapei,mapie)
      call isctr(ninnex,iro,iper,idum)
_ENDIF
_IF(cray,t3e)
      call fmove(idum,iro,ninnex)
_ELSE
      call icopy(ninnex,idum,1,iro,1)
_ENDIF
c
c...  print
c
      if(.not.oprint(28)) then
       write(iwr,4456)
4456   format(//5x,' by detsym reordered orbital sequence'/5x,38('=')/
     * //4('     function  e/i label')/
     *   5x,91('-'))
       write(iwr,4457) (mapie(loop),loop,loop=1,ninnex)
4457   format(4(i13,i11))
      endif
c
9999  return
      end
_IF(hpux11)
c HP compiler bug JAGae56123
c$HP$ OPTIMIZE LEVEL2
_ENDIF
      subroutine permco(conf,nint,ip,iocc,newocc)
c
c     permute the orbitals in a occupation-scheme ('conf') according
c     to 'ip'
c     ** see cpack , for packing convention **
c
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer conf,z3,z1,z0
_ELSEIF(convex)
      integer *8 conf,z3,z1,z0
_ELSEIF(i8drct)
      integer *8 conf,z3,z1,z0,shift
      external shift
_ELSE
      REAL conf,z3,z1,z0
_ENDIF
      dimension ip(nint),iocc(nint),newocc(nint),conf(2)
c
INCLUDE(common/cdcryi)
c
_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
_IF(cray,t3e,convex,i8,i8drct)
      data z0,z1,z3/0,1,3/
_ELSEIF(bits8)
      call pad8(0,z0)
      call pad8(1,z1)
      call pad8(3,z3)
_ELSE
      z0 = pad(0)
      z1 = pad(1)
      z3 = pad(3)
_ENDIF
      call upack(conf,nint,iocc)
c
_IF1()      do 10 i=1,nint
_IF1()   10 newocc(ip(i)) = iocc(i)
_IF(cray)
      call scatter(nint,newocc,ip,iocc)
_ELSE
      call isctr(nint,iocc,ip,newocc)
_ENDIF
c
      ns = n64
      iw = 1
_IF(cray,t3e,i8)
      conf(1) = 0
_ELSE
      conf(1) = z0
_ENDIF
      do 20 i=1,nint
         if (ns.gt.0) go to 15
c...  new word
         iw = iw + 1
         ns = n64
         conf(iw) = z0
15       ns = ns - 2
_IF(cray,t3e)
         if (newocc(i).eq.2) conf(iw) = or(conf(iw),shift(z3,ns))
         if (newocc(i).eq.1) conf(iw) = or(conf(iw),shift(z1,ns))
_ELSEIF(convex,i8)
         if(newocc(i).eq.2)conf(iw)=ior(conf(iw),ishftc(z3,ns,64))
         if(newocc(i).eq.1)conf(iw)=ior(conf(iw),ishftc(z1,ns,64))
_ELSEIF(i8drct)
         if(newocc(i).eq.2)conf(iw)=ior(conf(iw),shift(z3,ns))
         if(newocc(i).eq.1)conf(iw)=ior(conf(iw),shift(z1,ns))
_ELSEIF(bits8)
         if (newocc(i).eq.2) then
          call shift8(z3,ns,tmp8)
          call dor8(conf(iw),tmp8,conf(iw))
         endif
         if (newocc(i).eq.1) then
          call shift8(z1,ns,tmp8)
          call dor8(conf(iw),tmp8,conf(iw))
         endif
_ELSE
         if (newocc(i).eq.2) conf(iw) = dor(conf(iw),SHIFT(z3,ns))
         if (newocc(i).eq.1) conf(iw) = dor(conf(iw),SHIFT(z1,ns))
_ENDIF
20    continue
c
      return
      end
_IF(hpux11)
c$HP$ OPTIMIZE LEVEL3
_ENDIF
      subroutine purmit(iar,iscr,ip,nn)
c
c     permute the array iar  according to ip  / iscr is scratch
c
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension iar(nn),iscr(nn),ip(nn)
c
      do 10 i=1,nn
10    iscr(i) = iar(i)
c
      do 20 i=1,nn
20    iar(i) = iscr(ip(i))
c
      return
      end
      subroutine inreo2(conf,nword,nref,iwr)
c
c...  issue a internal reorder to get the symmetries alligned
c...  for extra orbitals added due to mrdci-option
c
c...  reorder the occupations in the reference 'conf'igurations
c...  orbital restrictions or 'eigen' values are not cared for
c
c...  common/expanc/(iper,dum,idumm) is used for scratch
c
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer conf
_ELSEIF(convex,i8drct)
      integer *8 conf
_ELSE
      REAL conf
_ENDIF
INCLUDE(common/sizes)
INCLUDE(common/cic)
INCLUDE(common/helpr)
INCLUDE(common/orb)
INCLUDE(common/symci)
INCLUDE(common/mapp)
INCLUDE(common/table)
c
      common/expanc/  dum(nd200),iper(nd200),idumm(100)
c
      dimension conf(nword,nref)
c
c...    first sort out new symmetry partitioning
c
      call setsto(nirr,0,norbi)
      call setsto(nirr,0,norbe)
      do 4  i=1,nint
4        norbi(iro(i)) = norbi(iro(i)) + 1
      if (nint.eq.ninnex) go to 6
      nint1 = nint + 1
      do 5  i=nint1,ninnex
         norbe(iro(i)) = norbe(iro(i)) + 1
5     continue
c
c...  new  symmetry result printing
c
6     write(iwr,601) sname
      write(iwr,602) (norbi(i),i=1,nirr)
      write(iwr,603) (norbe(i),i=1,nirr)
c
601   format(//,'-----  new symmetry assignment ------',/,
     1          '         group  ',a4)
602   format(// ' division of the internal orbitals ',8i4)
603   format(   ' division of the external orbitals ',8i4)
c
      modei = .false.
c
      do 10 i=1,nd200
      iper(i) = nd200
10    mapei(i) = nd200
c
c...  reorder the reordering / keep track of permutation in iper
c
      mcount = 0
      do  50 i=1,nirr
         do 40 loop=1,nint
            if (iro(loop).ne.i) go to 40
               mcount = mcount + 1
               mapei(mapie(loop)) = mcount
               iper(loop) = mcount
40       continue
50    continue
c
      loopb = nint + 1
      if (loopb.gt.ninnex) go to  200
      do 100 i=1,nirr
         do 90 loop=loopb,ninnex
            if (iro(loop).ne.i) go to 90
               mcount = mcount + 1
               mapei(mapie(loop)) = mcount
               iper(loop) = mcount
90       continue
100   continue
c
c...  now use iper to reorder  conf
c
200   if (nref.le.0) go to 500
      do 450 i=1,nref
         call permco(conf(1,i),nint,iper,dum,idumm)
450   continue
c
c...  put mapei and iro right
c
500   do 510 loop=1,nd200
510   mapie(mapei(loop)) = loop
      mcount = 0
      do 550 i=1,nirr
         nn = norbi(i)
         if (nn.le.0) go to 550
         do 530 loop=1,nn
            mcount = mcount + 1
            iro(mcount) = i
530      continue
550   continue
      do 600 i=1,nirr
         nn = norbe(i)
         if (nn.le.0) go to 600
         do 580 loop=1,nn
            mcount = mcount + 1
            iro(mcount) = i
580      continue
600   continue
c
c...  print
c
      write(iwr,4456)
4456  format(/5x,'reordered orbital sequence',/5x,26('=')/
     *       //4('     function  e/i label')/
     *  5x,91('-'))
      write(iwr,4457) (mapie(loop),loop,loop=1,ninnex)
4457  format(4(i13,i11))
c
      return
      end
      subroutine clrefc(iconf,nword,nlast,nint)
c
c...    clear refererence configurations from rubbish and
c...    occupied external orbitals
c...     called by rjbw
c
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf,mask,iiii,z0
_ELSEIF(convex,i8drct)
      integer *8 iconf,mask,iiii,z0
_ELSE
      REAL iconf,mask,iiii,z0
_ENDIF
      common/maskc/mask(64)
INCLUDE(common/cdcryi)
      dimension iconf(5,nlast)
      dimension iiii(4)
c
_IF(cray,t3e,convex,i8,i8drct)
      data z0/0/
_ELSEIF(bits8)
      call pad8(0,z0)
_ELSE
       z0 = pad(0)
_ENDIF
      iiii(1) = mask((min(nint,n32))*2)
      iiii(2)=z0
      if (nint.gt.n32) iiii(2) = mask((nint-n32)*2)
      iiii(3) = z0
      iiii(4) = z0
c
      nwm1 = nword - 1
      if (nwm1.gt.4) call errors(99999)
      do 20 j=1,nwm1
         do 10 i=1,nlast

_IF(cray,t3e,convex,i8,i8drct)
         iconf(j,i) = iand(iconf(j,i),iiii(j))
_ELSEIF(bits8)
         call dand8(iconf(j,i),iiii(j),iconf(j,i))
_ELSE
         iconf(j,i) = dand(iconf(j,i),iiii(j))
_ENDIF
10       continue
20    continue
c
      return
      end
      subroutine sccor(cr,conf,ngot)
c
c...  control calculation of size-consistency corrections
c...  for MR-SDCI wavefunctions
c...  prepared for various definitions of PSI-0
c...  this version takes : PSI(0) = PSI(REF)
c...  use davidv for diagonalisations in vacuum space
c...  now all dependy on external core-partitioning is removed
c...     cr(ngot) is now the usable area
c
      implicit REAL  (a-h,o-z)
_IF(cray,t3e,i8)
      integer conf
_ELSEIF(convex,i8drct)
      integer *8 conf
_ELSEIF(bits8)
      REAL conf
_ELSE
_ENDIF
c
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/auxh)
INCLUDE(common/dir_erg)
INCLUDE(common/cntrl)
INCLUDE(common/ccntl)
INCLUDE(common/disktl)
c
      dimension cr(ngot),conf(5,*)
_IF(parallel)
c     print *,' parallel size-consistency disabled '
c     return
_ENDIF
c
c...  count # csf's for proper correction => nrfcsf
c
      nrfcsf = 0
      do 10 i=1,nref
         call upack2(conf(3,i),jj,ii)
         nrfcsf = nrfcsf + ii
10    continue
c
      if (nrfcsf.le.0) return
      if (nrfcsf.gt.nval) then
c
c ... dont kill but warn ..
c    +   call caserr('ERROR ** Warn Joop, nrfcsf>nval')
c
         write(iwr,603) 
603      format(/
     +     ' *********************************************'/
     +     ' *** sorry nrfcsf > nval in SC corrections ***'/
     +     ' ***                warn JVL               ***'/
     +     ' *********************************************'//)
         return
      endif
c
c...  set up special addresses for mini-davidson
c
      kc = 1
      kz = kc + nval
      kbb = kz + nval
      kend = kbb + maxpa(4)**2
c
c...   check if core suffices
c
      navail = ngot - kend
      maxdav = navail/(2*nrfcsf)
      if (maxdav.lt.4) then
         write(iwr,601) nrfcsf,navail
601      format(' *** sorry not enough core for SC-corrections ',
     *          'with',i6,' vacuum states and',i6,' available ***')
         return
      end if
c
      efinal = h11
c
      if (nrfcsf.gt.1) then
        write(iwr,600)
600   format(///1x,52('@'),/
     *       '    MULTI-REFERENCE SIZE CONSISTENCY CORRECTIONS',
     *       /,1x,52('@'))
         write(iwr,604)
604      format(/'        === PSI-0 = PSI-ref ===')
      else
        write(iwr,'(//1x)')
      end if
c
c...  model D : PSI-O = reference  CI
c
      mmax = 100
      if (nrfcsf.eq.1) mmax = 0
_IF(ga)
      call rdedx(cr(kc+1),nval,index(26),numcz)
_ELSE
      call rdedx(cr(kc+1),nval,index(26),numci)
_ENDIF
      if (nval.gt.nrfcsf) 
     +    call vclr(cr(kc+nrfcsf+1),1,nval-nrfcsf)
c
      call davidv(cr(kc+1),cr(kz+1),cr(kbb+1),nrfcsf,e0,mmax,
     *            thresh,eshift,cr(kend+1),maxdav,conf,iwr,.true.)
c
      ecorr = efinal - e0
_IF(ga)
      call rdedx(cr(kz+1),nval,index(26),numcz)
_ELSE
      call rdedx(cr(kz+1),nval,index(26),numci)
_ENDIF
      c0 = ddot(nrfcsf,cr(kc+1),1,cr(kz+1),1)
      call sccors(nrfcsf,ecorr,c0,efinal,nelec,iwr)
c
      return
      end
      subroutine sccors(nref,ecorr,c0,ef,nelec,iwr)
c
c ********************************************************
c *                                                      *
c * size-consistency corrections for (multi-reference) ci*
c * c0 corr, ef  are determined by calling routine       *
c *                                                      *
c *references :                                          *
c * s.r. langhoff, e.r. davidson                         *
c * int.j. quant.chem. 8 ,61 , (1974)                    *
c *                                                      *
c * j.a. pople,r. seeger,r.krishnan                      *
c * int.j.quant.chem. s11, 149 ,(1977)                   *
c *                                                      *
c * per.e.m.siegbahn                                     *
c * chem.phys.lett. 55, 386, (1978)                      *
c *                                                      *
c * r.alrichs private communication                      *
c *                                                      *
c * r.j. vos & j.h. van lenthe                           *
c *                                                      *
c ********************************************************
c
      implicit REAL  (a-h,o-z)
c
c *****davidson corrections****************
c
      edav = ecorr*(1 - c0**2)
c
c ****siegbahn corrections***************
c
      esieg = ecorr*(1-c0**2)/c0**2
c
c***************pople corrections *******************
c
      thetha = dacos(c0)
      th2 = 2*thetha
      tn =  dtan(th2)
      cs2 = dcos(th2)
      tnk = tn**2
      csr = 1.0d0/cs2
      roemer = 2*(csr-1.0d0)
      wrt = (nelec*nelec) + (2*nelec*tnk)
      wrt = dsqrt(wrt)
      teller = wrt - nelec
      if (dabs(roemer).lt.1.0d-14) then
         epop = 0.0d0
      else
         epop = ( (teller/roemer) - 1) *ecorr
      end if
c
c ****************total energies *****************************
c
      edvt = edav +ef
      ept = epop + ef
      esiegt = esieg + ef
c
c ************************output**********************************
c
      write(iwr,101) nref,ecorr,c0
      write(iwr,102)
      write(iwr,103) edav,edvt,esieg,esiegt,epop,ept
101   format(' ****************************************************',
     1       /,'      size-consistency corrected (mr)-ci results ',
     2        /,'      nrefcsf',i3,' ecorr',f14.9,' c0 ',f8.6)
102   format(/12x,'   sc-correction  ', 7x, ' energy (total)' )
103   format(/' davidson ',f19.10,5x,f19.10,
     1       /' siegbahn ',f19.10,5x,f19.10,
     2       /' pople    ',f19.10,5x,f19.10)
      write(iwr,104)
104   format(' ****************************************************')
c
      return
      end
      subroutine davidv(c,z,w,ntotav,h11,maxcyc,thresh,eshift,cz,
     *                  maxdav,conf,iwr,oprdvd)
c
      implicit REAL  (a-h,o-z)
c
_IF(cray,t3e,i8)
      integer conf(*)
_ELSEIF(convex,i8drct)
      integer *8 conf(*)
_ELSE
      REAL conf(*)
_ENDIF
      logical oprdvd
c
      dimension c(*),z(*),w(*),cz(ntotav,2,maxdav)
c...
c...  davidson  diagonalizer with lock mode and no fancies
c...  in store mode / vacuum only
c
c...  *** c,z should be should be at least nval long
c...  *** remember that the read cannot stop at exactly ntotav elements
c...  *** all dependency on external core-partitioning is gone =>vrijkl
c...  Since rdedx can't read partial "records", reads nval now
c
c...  used for SC-correction determination ; so no CEPA needed
c...  the start-vector is supplied in c , usually through
c...  a read(c(khbc+1),nval,index(26),8), or from a previous davidv
c...  results are returned in c
c...
c...  the routine is derived from the Cyber205 davidson
c...  it's assumed that blank common is enough to store the
c...  rather short expansion vectors // otherwise I'll have to
c...  make it write to file and upset the disk
c...
c...
      logical jeqi
INCLUDE(common/sizes)
INCLUDE(common/helpr)
INCLUDE(common/cntrl)
INCLUDE(common/disktl)
      common/dialoc/floc(1275),f2loc(1275),xn(50),xns(50),ilifq(50)
      common/junk/igenz(4088),eig(50),g(16863)
      dimension q(2500)
      equivalence (q(1),g(1276))
INCLUDE(common/timez)
INCLUDE(common/restar)
c
c...  see how may expansion-vectors we can handle
c
      ispace =   min(maxdav,50)
      if (ispace.lt.3) call caserr('not enough space in davidv')
c
      icyc = 0
      if (maxcyc.gt.0.and.oprdvd) write(iwr,2000) ntotav
2000  format(3x,' vacuum davidson for',i7,'-dim. subspace',/
     *3x,46('-')/'   cycle  ',
     *5x,'energy(au)',6x,'tester',2x,'variance'/3x,46('-'))
c
c.... get start c-vector, normalise
c
      if (ntotav.gt.nval) then
         write(iwr,3001) ntotav, nval
3001     format(/' ntotav = ',i6/
     *           ' nval   = ',i6/)
         call caserr('davidv called for non-vacuum')
      endif
c
10    continue
      old = cpulft(1)
      call vclr(c(ntotav+1),1,nval-ntotav)
      top = ddot(ntotav,c,1,c,1)
      cnorm = 1.0d0/dsqrt(top)
      call dscal(ntotav,cnorm,c,1)
c
c...  calculate Hc and starting H11
c
      call vclr(z,1,ntotav)
c
c...   calculate Z = HC in the vacuum space
c...   scratch-space for vvijkl (maxpa(4)**2)
c
      call vrijkl(index(1),ntotav,z,c,w,conf)
_IF(parallel)
      call pg_dgop(4712,z,ntotav,'+')
_ENDIF
      h11 = ddot(ntotav,c,1,z,1)
c
      if (maxcyc.le.0) return
c
c...  prepare for iterations
c
      i199=index(199)
      call ibasgn(ispace,0,50,ilifq)
      xn(1)=1.0d0
      xns(1)=1.0d0
      floc(1)=h11
       var=ddot(ntotav,z,1,z,1)
      f2loc(1)=var
      call dcopy(ntotav,c,1,cz(1,1,1),1)
      call dcopy(ntotav,z,1,cz(1,2,1),1)
c
      do 6666 joop=2,ispace
c
c...  save current c-vector in z-place in case of convergence
c
      call dcopy(ntotav,c,1,cz(1,2,joop),1)
c
         icyc = icyc + 1
c
         jeqi=joop.eq.ispace
         ioop=joop-1
         var=var-h11*h11
c... generate pert vector
_IF(cyber205)
         call gtrian(ntotav,h11,c,z,c)
_IF(ga)
         call getvacdiag(z)
_ELSE
         call rdedx(z,nval,i199,numci)
_ENDIF
         call vsav(ntotav,eshift-h11,z,z)
         call mrgle(ntotav,0.05d0,z)
         call vvdv(ntotav,c,c,z)
         call mrgge(ntotav,0.2d0,c)
         tester=absmax(ntotav,0.0d0,c)
_ENDIF(cyber205)
_IF(cray)
         call daxpy(ntotav,-h11,c,1,z,1)
_IF(ga)
         call getvacdiag(z)
_ELSE
         call rdedx(c,nval,i199,numci)
_ENDIF
         call szero(c(ntotav+1),nval-ntotav)
         top=h11-eshift
         do 1 loop=1,ntotav
         temp=top-c(loop)
         temp=z(loop)/cvmgp(temp,dsign(0.05d0,temp),dabs(temp)-0.05d0)
         c(loop)=cvmgp(dsign(0.2d0,temp),temp,dabs(temp)-0.2d0)
    1    continue
         loop=idamax(ntotav,c,1)
         tester=dabs(c(loop))
_ENDIF(cray)
_IF(ibm,vax)
         call gtrian(ntotav,h11,c,z,c)
_IF(ga)
         call getvacdiag(z)
_ELSE
         call rdedx(z,nval,i199,numci)
_ENDIF
         call vsav(ntotav,eshift-h11,z,z)
         call mrgle(ntotav,0.05d0,z)
         call vvdv(ntotav,c,c,z)
         call mrgge(ntotav,0.2d0,c)
_IF(ibm)
         loop=idamax(ntotav,c,1)
         tester=dabs(c(loop))
_ELSE
         tester=0.0d0
         do 2 loop=1,ntotav
         temp=dabs(c(loop))
         if(temp.gt.tester)tester=temp
    2    continue
_ENDIF
_ENDIF(ibm,vax)
_IFN(cray,cyber205,ibm,vax)
         call vsmsb(c,1,h11,z,1,c,1,ntotav)
_IF(ga)
         call getvacdiag(z)
_ELSE
         call rdedx(z,nval,i199,numci)
_ENDIF
         call vsadd(z,1,eshift-h11,z,1,ntotav)
         call mrgle(ntotav,0.05d0,z)
_IFN1(f)         call vdiv(c,1,z,1,c,1,ntotav)
         call mrgge(ntotav,0.2d0,c)
         call maxmgv(c,1,tester,loop,ntotav)
_ENDIF(not,cray,cyber205,ibm,vax)
         top=cpulft(1)
         if (oprdvd.or.(maxcyc-icyc).le.10) write(iwr,111)icyc,
     +       h11,tester,var
111      format(i5,e23.14,e10.3,e11.4)
         if(tester.lt.thresh) goto 777
         if(icyc.ge.maxcyc)goto 440
c... orthogonalize pert vector
         x=ddot(ntotav,c,1,c,1)
9944     y=x
         do 609 koop=1,ioop
         x=ddot(ntotav,c,1,cz(1,1,koop),1)*xns(koop)
         call daxpy(ntotav,-x,cz(1,1,koop),1,c,1)
609      continue
         x=ddot(ntotav,c,1,c,1)
         if(x.lt.(0.05d0*y))goto 9944
         x=1.0d0/x
         xns(joop)=x
         y=dsqrt(x)
         xn(joop)=y
         call dcopy(ntotav,c,1,cz(1,1,joop),1)
c
c...   calculate Z = HC in the vacuum space
c...   C is in cr(khbc+1), Z in cr(khbz+1)
c...   scratch-space for vvijkl (maxpa(4)**2) in khbb + 1
c
      call vclr(z,1,ntotav)
      call vrijkl(index(1),ntotav,z,c,w,conf)
_IF(parallel)
         call pg_dgop(4713,z,ntotav,'+')
_ENDIF
         temp = ddot(ntotav,c,1,z,1)
c
c...     save z-vector and setup mini H-matrix
c
         call dcopy(ntotav,z,1,cz(1,2,joop),1)
         m=iky(joop)
         mm=m+joop
c... compute h and h**2 elements
         floc(mm)=temp*x
         f2loc(mm)=ddot(ntotav,z,1,z,1)*x
         do 5097 koop=1,ioop
            top=ddot(ntotav,cz(1,1,joop),1,cz(1,2,koop),1)
            bot=ddot(ntotav,cz(1,2,joop),1,cz(1,2,koop),1)
            temp=y*xn(koop)
            floc(m+koop)=temp*top
            f2loc(m+koop)=temp*bot
5097     continue
c
         call dcopy(mm,floc,1,g,1)
c
c... determine best solution (lock mode only)
c
         call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
         do 6669 loop=1,joop
6669     eig(loop)=-dabs(q(ilifq(loop)+1))
c
         ivrs=idmin(joop,eig,1)
         ivrs=ilifq(ivrs)
         var=expect(f2loc,q(ivrs+1),joop)
         h11=expect(floc,q(ivrs+1),joop)
c
c... generate best vector     and     h * best vector
c
      call vmul(q(ivrs+1),1,xn,1,eig,1,joop)
         call dscal(ntotav,eig(joop),c,1)
         do 621 koop=1,ioop
         call daxpy(ntotav,eig(koop),cz(1,1,koop),1,c,1)
621      continue
c
         if (jeqi) goto 43210
         call dscal(ntotav,eig(joop),z,1)
         do 622 koop=1,ioop
         call daxpy(ntotav,eig(koop),cz(1,2,koop),1,z,1)
622      continue
c
         top=cpulft(1)
         tlefti=timlim-top-2.0d0
         if(tlefti.lt.1.5d0*(top-old)) go to 440
         old=top
c
6666  continue
c
c...  fold
c
43210 continue
      temp=ddot(ntotav,c,1,c,1)
      call dscal(ntotav,dsqrt(1.0d0/temp),c,1)
       top = cpulft(1)
       if (oprdvd.or.(maxcyc-icyc).le.10) write(iwr,1001)icyc,top
1001   format(i8,f12.3,' space contracted')
cjvl
cjvl    it is (much) better to fold to the n lowest eigenvectors
cjvl    of the davidson matrix, instead of 1
cjvl
      goto 10
c
440   write(iwr,102)
102   format(///20x,14('*')//20x,'davidv no convergence'//20x,14('*'))
      irest = 9
c
777   continue
c
      call dcopy(ntotav,cz(1,2,joop),1,c,1)

c
      return
      end
      subroutine dbgciv( text,array,ncol,nrow,mcolr,mrowr)
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
      character *(*) text
      dimension array(*)
      dimension f(4),ia(4),ja(4)
      write(iwr,100)text
100   format(/1x,'***** ', a6 ,' *****')
      iseq = 0
      n = 0
      do i=1,ncol
      mapi = (i-1)*mcolr+1
      do j=1,nrow
       mapij = mapi + (j-1)*mrowr
       n = n + 1
       f(n)  = array(mapij)
       ia(n) = i
       ja(n) = j
       if(n.eq.4) then
       iseq = iseq + 1
       write(iwr,200)iseq,(ia(k),ja(k),f(k),k=1,n)
       n = 0
       endif
      enddo
      enddo
      if(n.ne.0) then
       iseq=iseq+1
       write(iwr,200)iseq,(ia(k),ja(k),f(k),k=1,n)
      endif
200   format(i6,2x,4(2i5,f12.5))
      return
      end
_IFN(parallel)
      subroutine genc2(cr,iconf)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
c... to generate perturbation  c-vector
      dimension cr(*),iconf(5,*)
INCLUDE(common/sizes)
INCLUDE(common/diactl)
      common/moco/isi,ic
INCLUDE(common/dir_erg)
INCLUDE(common/ccntl)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct
INCLUDE(common/disktl)
INCLUDE(common/cntrl)
INCLUDE(common/symci)
INCLUDE(common/outctl)
      z12=0.0d0
      z22=0.0d0
      z11=0.0d0
      h12=0.0d0
      s12=0.0d0
      s22=0.0d0
      tester=0.0d0
      if(instor)goto 888
c----------- out of store version --------------------------------------
      i26=index(26)
      i27=index(27)
      i28=index(28)
      if(nval)101,100,101
101   call rdedx(cr(kumpc+1),nval,i26,numci)
      call rdedx(cr(khbill+1),nval,i27,numci)
      call vvc2(cr,cr(kumpc+1),cr(khbill+1))
      call wrt3(cr(kumpc+1),nval,i28,numci)
100   modev=s22.lt.(1.0d-18)
      if(ndoub.eq.0)go to 200
_IF(cray,t3e,convex,i8,i8drct)
      iconff=iconf(2,1)
_ELSE
      iconff=ipad(iconf(2,1))
_ENDIF
      call rdedx(cr(kumpc+1),ndoub,i26+iconff,numci)
      call rdedx(cr(khbill+1),ndoub,i27+iconff,numci)
      call ddc2(cr,iconf,cr(kumpc+1),cr(khbill+1))
      call wrt3(cr(kumpc+1),ndoub,i28+iconff,numci)
200   if(nmin2)301,999,301
301   do 300 ic=nstrt2,nlast2
_IF(cray,t3e,convex,i8,i8drct)
      isi=iconf(4,ic)
_ELSE
      isi=ipad(iconf(4,ic))
_ENDIF
      call upack2(iconf(3,ic),lamt,lams)
      call upack2(iconf(2,ic),kc,jc)
      ns=lams*norbsi(isi)
      nt=lamt*norbtr(isi)
      nsnt=ns+nt
      call rdedx(cr(kumpc+1),nsnt,i26+jc,numci)
      call rdedx(cr(khbill+1),nsnt,i27+jc,numci)
      if(ns)401,400,401
401   call ssc2(cr,cr(kumpc+1),cr(khbill+1),0,lams,jc,ns)
      call wrt3(cr(kumpc+1),ns,i28+jc,numci)
      if(nt)400,300,400
400   nsnt=kumpc+ns
         call ssc2(cr,cr(nsnt+1),cr(khbill+ns+1),1,lamt,kc,nt)
      call wrt3(cr(nsnt+1),nt,i28+kc,numci)
300   continue
      goto 999
c------------ in store version ----------------
888   if(nnst.ne.0)call vvc2(cr,cr(khbc+1),cr(khbz+1))
       modev=s22.lt.(1.0d-18)
         if(nmin1.ne.0)call ddc2(cr,iconf,cr(khbcd+1),cr(khbzd+1))
      if(nmin2)1,999,1
1     do 2 ic=nstrt2,nlast2
      call upack2(iconf(1,ic),kc,jc)
       call upack2(iconf(2,ic),kcc,jcc)
      call upack2(iconf(3,ic),lamt,lams)
_IF(cray,t3e,convex,i8,i8drct)
      isi=iconf(4,ic)
_ELSE
      isi=ipad(iconf(4,ic))
_ENDIF
c... singlet
         call ssc2(cr,cr(khbcs+jc+1),cr(khbzs+jc+1),0,lams,jcc,
     *lams*norbsi(isi))
c... triplet
         call ssc2(cr,cr(khbct+kc+1),cr(khbzt+kc+1),1,lamt,kcc,
     *lamt*norbtr(isi))
2     continue
c... following obeyed by both versions
999   moded=nmin1.eq.0
      modest=nmin2.eq.0
      return
      end
      subroutine c22(ee,zz,cc,number)
      implicit REAL  (a-h,o-z),integer  (i-n)
       dimension zz(*),cc(*),ee(*)
INCLUDE(common/sizes)
INCLUDE(common/outctl)
INCLUDE(common/dir_erg)
       do 999 loop=1,number
        z=zz(loop)
        c=cc(loop)
      x=z-h11*c
      z11=z11+x*x
       if( dabs(c).le.scream)goto 998
      top=0.0d0
        goto 999
998   csq=1.0d0/(1.0d0-c*c)
      croot= dsqrt(csq)
      g12=x*croot
      g22=(ee(kumpe+loop)+eshift-c*(x+z))*csq
      top=b22(g12,g22)
      top=top*croot
      top=top/(1.0d0-top*c)
      atop= dabs(top)
      if(tester.lt.atop)tester=atop
      if(atop.gt.(1.0d0))top=top/atop
      s12=c*top+s12
       s22=top*top+s22
      h12=z*top+h12
999   cc(loop)=top
       return
      end
      subroutine vvc2(ee,c,z)
c... pert c vector for vacuum states
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension c(*),ee(*),z(*)
INCLUDE(common/outctl)
INCLUDE(common/disktl)
INCLUDE(common/cntrl)
INCLUDE(common/ccepa)
       call rdedx(ee(kumpe+1),nval,index(199),numci)
ccepa
      if (ifcepa.gt.0.and.icepa.ge.10) call 
     1 cepadvv(ee(kumpe+1),z,c,ee)
ccepa
      call c22(ee,z,c,nval)
      return
      end
      subroutine ddc2(ee,iconf,c,z)
c... pert c vector for doublet states
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension c(*),ee(*),z(*),iconf(5,*)
INCLUDE(common/disktl)
INCLUDE(common/cntrl)
INCLUDE(common/outctl)
INCLUDE(common/ccepa)
       call rdedx(ee(kumpe+1),ndoub,index(198),numci)
      if (ifcepa.gt.0) call cepddd(iconf,ee(kumpe+1),z,c)
       call c22(ee,z,c,ndoub)
      return
      end
      subroutine ssc2(ee,c,z,it,lam,jckc,num)
c... pert c vector for  n-2  states
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension z(*),ee(*),c(*)
INCLUDE(common/outctl)
INCLUDE(common/disktl)
INCLUDE(common/ccepa)
_IF(linux)
      external fget
_ENDIF
      if(lam)777,999,777
777     call rdedx(ee(kumpe+1),num,index(199)+jckc,numci)
      if (ifcepa.gt.0) 
     1   call cepadst(ee(kumpe+1),z,c,it,num,ee)
       call c22(ee,z,c,num)
999   return
       end
      subroutine gc3(b,a,n,ibl)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension a(*),b(*)
INCLUDE(common/sizes)
INCLUDE(common/dir_erg)
INCLUDE(common/disktl)
      if(n.eq.0)go to 999
      call rdedx(a,n,ibl,numci)
      call dscal(n,bbb,b,1)
      call daxpy(n,aaa,a,1,b,1)
999   return
      end
      subroutine genc3(cr,iconf)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
INCLUDE(common/sizes)
c... to revise c and z vectors
      dimension cr(*),iconf(5,*)
INCLUDE(common/disktl)
INCLUDE(common/diactl)
INCLUDE(common/symci)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct
INCLUDE(common/ccntl)
INCLUDE(common/cntrl)
      common/scra  /jndex(6500)
INCLUDE(common/outctl)
_IF(linux)
      external fget
_ENDIF
      i26=index(26)
      if(instor)goto 888
c----------- out of store version --------------------------------------
      i28=index(28)
      do 51 looper=1,2
      nt=-ntotal
114   if(nt)511,555,555
511   call search(i28,numci)
      mus=khbill
      ncount=0
112   call fget(cr(mus+1),m,numci)
      i28=i28+1
      ncount=ncount+1
      nt=nt+m
      mus=mus+m
      jndex(ncount)=m
      if((mus+511).le.khgot2.and.nt.lt.0.and.ncount.lt.6500)goto 112
      call gc3(cr(khbill+1),cr(mus+1),mus-khbill,i26)
      call search(i26,numci)
      mus=khbill
      do 117 loop=1,ncount
      m=jndex(loop)
      call put(cr(mus+1),m,numci)
117   mus=mus+m
      i26=iposun(numci)
      goto 114
555   i26=index(27)
51    i28=index(29)
      return
c------------- in store version ----------------------------------------
888   lh=khbc
      do 1 looper=1,2
c... vacuum
      call gc3(cr(lh+1),cr(khbb+1),nval,i26)
c... doublet
      kh=lh+nval
_IF(cray,t3e,convex,i8,i8drct)
       klaus=i26+iconf(2,1)
       call gc3(cr(kh+1),cr(khbb+1),ndoub,klaus)
_ELSE
      call gc3(cr(kh+1),cr(khbb+1),ndoub,i26+ipad(iconf(2,1)))
_ENDIF
      kh=kh+ndoub
c... n-2 states
      if(nmin2)2,11,2
2      ll=kh+nsingl
      do 3 loop=nstrt2,nlast2
      call upack2(iconf(2,loop),it,is)
      call upack2(iconf(3,loop),nt,ns)
_IF(cray,t3e,convex,i8,i8drct)
      isi=iconf(4,loop)
_ELSE
      isi=ipad(iconf(4,loop))
_ENDIF
      ns=norbsi(isi)*ns
      nt=norbtr(isi)*nt
c... singlet
      call gc3(cr(kh+1),cr(khbb+1),ns,i26+is)
      kh=kh+ns
c... triplet
      call gc3(cr(ll+1),cr(khbb+1),nt,i26+it)
3     ll=ll+nt
   11      call dumpcz(iconf,cr(lh+1),i26)
        i26=index(27)
1       lh=khbz
       return
      end
      subroutine genz3(cr,iconf)
c... to compute z12 and z22
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension cr(*),iconf(5,*)
INCLUDE(common/cntrl)
INCLUDE(common/symci)
INCLUDE(common/disktl)
INCLUDE(common/diactl)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt
INCLUDE(common/ccntl)
INCLUDE(common/outctl)
_IF(linux)
      external fget
_ENDIF
      i27=index(27)
      if(instor)goto 888
c-------- out of store version----------
      i29=index(29)
      nt=-ntotal
114   if(nt)111,11,11
111   call search(i29,numci)
      mus=khbill
      i29x=i29
112   call fget(cr(mus+1),m,numci)
      i29=i29+1
      nt=nt+m
      mus=mus+m
      if((mus+511).le.khgot2.and.nt.lt.0)goto 112
      call gz3(cr(khbill+1),cr(mus+1),mus-khbill,i27)
      i27=i27+i29-i29x
      goto 114
c-------- in store   version -----------
c... vacuum
888   call gz3(cr(khbz+1),cr(khbb+1),nval,i27)
c... doublet
_IF(cray,t3e,convex,i8,i8drct)
      klaus=i27+iconf(2,1)
      call gz3(cr(khbzd+1),cr(khbb+1),ndoub,klaus)
_ELSE
      call gz3(cr(khbzd+1),cr(khbb+1),ndoub,i27+ipad(iconf(2,1)))
_ENDIF
      if(nmin2)2,11,2
2     do 3 loop=nstrt2,nlast2
      call upack2(iconf(2,loop),it,is)
      call upack2(iconf(3,loop),nt,ns)
_IF(cray,t3e,convex,i8,i8drct)
      isi=iconf(4,loop)
_ELSE
      isi=ipad(iconf(4,loop))
_ENDIF
      call upack2(iconf(1,loop),mut,mus)
c... singlets
      call gz3(cr(khbzs+mus+1),cr(khbb+1),ns*norbsi(isi),i27+is)
c... triplets
3     call gz3(cr(khbzt+mut+1),cr(khbb+1),nt*norbtr(isi),i27+it)
11    return
      end
      subroutine gz3(z2,z1,n,ibl)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension z1(*),z2(*)
INCLUDE(common/sizes)
INCLUDE(common/dir_erg)
INCLUDE(common/disktl)
      if(n)1,999,1
1     call rdedx(z1,n,ibl,numci)
      z12=ddot(n,z1,1,z2,1)+z12
      z22=ddot(n,z2,1,z2,1)+z22
999   return
      end
      subroutine david3(cr,iconf,c,z,iwr)
      implicit REAL  (a-h,o-z), integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      dimension cr(*),iconf(5,*),c(*),z(*),q(2500)
c...
c...  davidson  diagonalizer -- instore mode
c...
      logical jeqi
INCLUDE(common/helpr)
INCLUDE(common/dir_erg)
INCLUDE(common/diactl)
INCLUDE(common/cntrl)
INCLUDE(common/disktl)
      common/stoctl/khbill,khbz,khbc,khbb
INCLUDE(common/auxh)
      common/dialoc/floc(1275),f2loc(1275),xn(50),xns(50),ilifq(50)
      common/junk/igenz(4088),eig(50),g(16863)
INCLUDE(common/restar)
INCLUDE(common/timez)
INCLUDE(common/ccepa)
INCLUDE(common/mr_cepctl)
INCLUDE(common/mr_pairs)
c     equivalence (q(1),g(1276))
      dimension fsloc(1275),f2sloc(1275)
      character *8 text,txcepa
      logical odumpz
      data text/'       '/,txcepa/'  cepa  '/
c
      if (icepa.ge.0)lvar=.true.
      fact=2.0d0
      fact2=2.0d0
      var=0.0d0
      i26=index(26)
      i27=index(27)
      i28=index(28)
      i199=index(199)
      call ibasgn(lcall+lcall,i28,ispsma,index(28))
      call ibasgn(lcall,0,lcall,ilifq)
      i29=index(29)
      xn(1)=1.0d0
      xns(1)=1.0d0
      ntot2=ntotal+ntotal
9999  floc(1)=h11
      old=cpulft(1)
       if(lvar)var=ddot(ntotal,z,1,z,1)
       f2loc(1)=var
      call wrt3(c,ntotal,i28,numci)
      if(odebug(40)) then
       call dbgvecv('david3-c',c,ntotal)
      endif
      call wrt3(z,ntotal,i29,numci)
      if(odebug(40)) then
       call dbgvecv('david3-z',z,ntotal)
      endif
ccepa
      if (icepa.ge.0) then
         if (moddia.gt.0) then
            call cntrci(iconf,c,z,0,0,cr(kcpcc+1),1,2,temp)
         else
            call vclr(cr(kcpcc+1),1,ncpcc*2)
            call cntrci(iconf,c,c,0,0,cr(kcpcc+1),0,1,temp)
         end if
         call dumpcc(cr(kcpcc+1),1,1)
c...  in case we come here after a space contraction ....
         if (ifcepa.eq.1) call hh2cep
     *     (cr(kcpcc+1),ncpcc,iconf,h11,h11,var,var,xn,1)
      end if
      fsloc(1) = floc(1)
      f2sloc(1) = f2loc(1)
      ncycep = 0
      h11o = h11
ccepa
      do 6666 joop=2,lcall
      jeqi=joop.eq.lcall
      odumpz=.true.
      ioop=joop-1
      h11sq=h11*h11
      if(lvar)var=var-h11sq
c... generate pert vector
      if (ifcepa.gt.0) then
         call cepvds(z,z,c,.false.,.true.,iconf)
      endif
ccepa  recompute var for cpea + emin or lock - modes
      if (ifcepa.gt.0.and.moddia.le.0) var = ddot(ntotal,z,1,z,1)
     *   -h11*h11
ccepa
_IF(cray)
      call daxpy(ntotal,-h11,c,1,z,1)
      call getcz(iconf,c,i199)
      if (ifcepa.gt.0) call cepvds(c,z,c,.true.,.false.,iconf)
      top=h11-eshift
      do 1 loop=1,ntotal
      temp=top-c(loop)
      temp=z(loop)/cvmgp(temp,dsign(0.05d0,temp),dabs(temp)-0.05d0)
      c(loop)=cvmgp(dsign(0.2d0,temp),temp,dabs(temp)-0.2d0)
    1 continue
      loop=idamax(ntotal,c,1)
      tester=dabs(c(loop))
_ELSEIF(ibm,vax)
      call gtrian(ntotal,h11,c,z,c)
      call getcz(iconf,z,i199)
      if (ifcepa.gt.0) call cepvds(z,z,c,.true.,.false.,iconf)
      call vsav(ntotal,eshift-h11,z,z)
      call mrgle(ntotal,0.05d0,z)
      call vvdv(ntotal,c,c,z)
      call mrgge(ntotal,0.2d0,c)
_IF(ibm)
      loop=idamax(ntotal,c,1)
      tester=dabs(c(loop))
_ELSE
      tester=0.0d0
      do 2 loop=1,ntotal
      temp=dabs(c(loop))
      if(temp.gt.tester)tester=temp
    2 continue
_ENDIF
_ELSE
      call vsmsb(c,1,h11,z,1,c,1,ntotal)
      call getcz(iconf,z,i199)
      if (ifcepa.gt.0) call cepvds(z,z,c,.true.,.false.,iconf)
      call vsadd(z,1,eshift-h11,z,1,ntotal)
      call mrgle(ntotal,0.05d0,z)
      call vdiv(c,1,z,1,c,1,ntotal)
      call mrgge(ntotal,0.2d0,c)
      call maxmgv(c,1,tester,loop,ntotal)
_ENDIF
      top=cpulft(1)
      if(odebug(40)) then
       call dbgvecv('david3-c',c,ntotal)
       call dbgvecv('david3-z',z,ntotal)
      endif
      write(iwr,111)icyc,top,h11,tester,var,text
111   format(1x,i7,f12.3,3e20.11,a8)
      if(tester.lt.thresh)goto 777
      if(icyc.ge.maxcyc)goto 440
c...
c...  possible jacobi davidson preconditioning on vacuum space
c...  here khbb is the first-1 word available (khbz is waisted)
c
      call gdvdmr(index(1),index(199),i26,i27,
     +            khbb,nval,maxpa(4),
     +            c,z,h11,cr,iconf,numci,ncycg)

c
      if(malt)eshift=-eshift
c... orthogonalize pert vector
      x=ddot(ntotal,c,1,c,1)
9944  y=x
      joke=28
      do 609 koop=1,ioop
      call rdedx(z,ntotal,index(joke),numci)
      x=ddot(ntotal,c,1,z,1)*xns(koop)
      call daxpy(ntotal,-x,z,1,c,1)
609   joke=joke+2
      x=ddot(ntotal,c,1,c,1)
      if(x.lt.(0.05d0*y))goto 9944
      x=1.0d0/x
      xns(joop)=x
      y=dsqrt(x)
      xn(joop)=y
      call genz(cr,iconf,index,temp)
      m=iky(joop)
      mm=m+joop
c... compute h and h**2 elements
      floc(mm)=temp*x
      if(lvar)f2loc(mm)=ddot(ntotal,z,1,z,1)*x
ccepa
      if (icepa.ge.0) then
c...     c.c and c.z for current c/z vectors
         if (moddia.gt.0) then
            call cntrci(iconf,c,z,0,0,cr(kcpcc+1),1,2,temp)
         else
            call vclr(cr(kcpcc+1),1,ncpcc*2)
            call cntrci(iconf,c,c,0,0,cr(kcpcc+1),0,1,temp)
         end if
         call dumpcc(cr(kcpcc+1),joop,joop)
      end if
ccepa
      joke=29
      do 5097 koop=1,ioop
      call search(index(joke),numci)
      top=0.0d0
      bot=0.0d0
      loop=0
      if (icepa.ge.0) goto 5095
6098  moop=min(ntotal-loop,16863)
      call reads(g,moop,numci)
      top=ddot(moop,c(loop+1),1,g,1)+top
      if(lvar)bot=ddot(moop,z(loop+1),1,g,1)+bot
      loop=loop+16863
      if(loop.lt.ntotal)goto 6098
      goto 5096
5095  call rdedx(cr(kdavcz+1),ntotal,index(joke),numci)
      top=ddot(ntotal,c,1,cr(kdavcz+1),1)
      if(lvar)bot=ddot(ntotal,z,1,cr(kdavcz+1),1)
ccepa
      if (icepa.ge.0) then
          call rdedx(cr(kdavcz+1),ntotal,index(joke-1),numci)
          call cntrci(iconf,c,cr(kdavcz+1),0,0,cr(kcpcc+1),0,1,temp)
          if (moddia.gt.0) then
            call rdedx(cr(kdavcz+1),ntotal,index(joke),numci)
            call cntrci(iconf,c,cr(kdavcz+1),0,0,cr(kcpcc+1),0,2,temp)
            call rdedx(cr(kdavcz+1),ntotal,index(joke-1),numci)
            call cntrci(iconf,z,cr(kdavcz+1),0,0,cr(kcpcc+1),0,3,top)
          else
            call vclr(cr(kcpcc+ncpcc+1),1,ncpcc*2)
          end if
          call dumpcc(cr(kcpcc+1),joop,koop)
        end if
ccepa
5096  temp=y*xn(koop)
      floc(m+koop)=temp*top
      f2loc(m+koop)=temp*bot
5097  joke=joke+2
ccepa
      call fmove(floc(m+1),fsloc(m+1),joop)
      call fmove(f2loc(m+1),f2sloc(m+1),joop)
      if (ifcepa.eq.1) go to 6000
ccepa
5098  if(moddia)6660,6661,6662
c... lock mode
6661  continue
      call dcopy(mm,floc,1,g,1)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
      do 6669 loop=1,joop
6669  eig(loop)=-dabs(q(ilifq(loop)+1))
      goto 6668
c... emin mode
6660  continue
      call dcopy(mm,floc,1,g,1)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
      goto 6668
c... vmin mode
6662  continue
      call vsma(floc,1,-(h11+h11),f2loc,1,g,1,mm)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
_IF(cray,ibm,vax)
      call vsav(joop,h11sq,eig,eig)
_ELSE
      call vsadd(eig,1,h11sq,eig,1,joop)
_ENDIF
6668  continue
      ivrs=idmin(joop,eig,1)
      ivrs=ilifq(ivrs)
      if(jeqi)goto 4010
      if (odumpz) then
         indc=index(joke-1)
         call wrt3(c,ntotal,index(joke-1),numci)
      endif
      if (odumpz) call wrt3(z,ntotal,index(joke),numci)
      if(lvar)var=expect(f2loc,q(ivrs+1),joop)
      h11=expect(floc,q(ivrs+1),joop)
c... generate best vector     and     h * best vector
4010  continue
      call vmul(q(ivrs+1),1,xn,1,eig,1,joop)
c....   next bit is strang ... Remco ?? was lt.10
      if (icepa.lt.0) then
        call dscal(ntot2,eig(joop),z,1)
        joke=28
        do 621 koop=1,ioop
        call search(index(joke),numci)
        loop=0
        temp=eig(koop)
6088    moop=min(ntotal-loop,16863)
        call reads(g,moop,numci)
        call daxpy(moop,temp,g,1,c(loop+1),1)
         loop=loop+16863
         if(loop.lt.ntotal)goto 6088
         if(jeqi)goto 621
        call search(index(joke+1),numci)
        loop=0
6087    moop=min(ntotal-loop,16863)
        call reads(g,moop,numci)
        call daxpy(moop,temp,g,1,z(loop+1),1)
         loop=loop+16863
         if(loop.lt.ntotal)goto 6087
621      joke=joke+2
         if(jeqi)goto 43210
         call dumpcz(iconf,c,i26)
         call dumpcz(iconf,z,i27)
       else
         call dumpcz(iconf,z,i27)
         call dscal(ntotal,eig(joop),c,1)
         joke=28
         do 622 koop=1,ioop
         call search(index(joke),numci)
         loop=0
         temp=eig(koop)
6089     moop=min(ntotal-loop,16863)
         call reads(g,moop,numci)
         call daxpy(moop,temp,g,1,c(loop+1),1)
         loop=loop+16863
         if(loop.lt.ntotal)goto 6089
622      joke=joke+2
         call dumpcz(iconf,c,i26)
      endif
ccepa
ccepa see if cepa proces may start / if so modify h and h**2
ccepa and start the micro-iteration-process first
ccepa
      if (icepa.ge.0.and.tester.le.critcep.and.joop.ge.itcepa) ifcepa=1
      if (ifcepa.ne.1) go to 6001
      text = txcepa
c... control/print micro-iterations
         top=cpulft(1)
         bot=0.0d0
c        call givtim(top,bot)
         temp = h11-h11o
         if (lprcep.gt.0.and.ncycep.gt.0.and.
     *       dabs(temp).gt.tester*crycep) write(6,114)
     *                        ncycep,h11o,temp
114   format(' cepa micro it.',i3,
     *       ' last e ',e20.13,' delta ',e10.3,'  cepa ')
         temp = dabs(temp)
         h11o = h11
         if (ncycep.gt.0.and.
     *     (ncycep.ge.mcycep.or.temp.le.tester*crycep)) go to 6001
c...  compute pair energies
         call cepair(cr,iconf,i26)
c...  what was the ci-energy again ?
         h11ci = expect(fsloc,q(ivrs+1),joop)
c...  compute cepa shifts
         call cepac(c,iconf)
c restore c and z
         call rdedx(c,ntotal,indc,numci)
         call getcz(iconf,z,i27)
         odumpz=.false.
c...  modify h (floc) and h**2 (f2loc) matrices
         ncycep = ncycep + 1
6000     call hh2cep(cr(kcpcc+1),ncpcc,iconf,floc,fsloc,f2loc,
     *    f2sloc,xn,joop)
         go to 5098
ccepa
6001  ncycep = 0
      h11o = h11
ccepa
      if(jeqi)goto 43210
      if (icepa.ge.0) then
      call getcz(iconf,z,i27)
      call dscal(ntotal,eig(joop),z,1)
      joke=29
      do 623 koop=1,ioop
       call search(index(joke),numci)
       loop=0
       temp=eig(koop)
6090   moop=min(ntotal-loop,16863)
       call reads(g,moop,numci)
       call daxpy(moop,temp,g,1,z(loop+1),1)
       loop=loop+16863
       if(loop.lt.ntotal)goto 6090
623   joke=joke+2
      call dumpcz(iconf,z,i27)
      endif
       top=cpulft(1)
       tlefti=timlim-top-fact2
       if(tlefti.lt.fact*(top-old)) go to 440
       old=top
6666   continue
43210 continue
      temp=ddot(ntotal,c,1,c,1)
      call dscal(ntotal,dsqrt(1.0d0/temp),c,1)
      call dumpcz(iconf,c,i26)
      call genz(cr,iconf,index(26),h11)
      call dumpcz(iconf,z,i27)
       top=cpulft(1)
       write(iwr,1001)icyc,top
1001   format(i8,2f12.3,' space contracted')
      goto 9999
440   write(iwr,102)
102   format(//20x,14('*')/20x,'no convergence'/20x,14('*'))
      irest=9
      goto 3333
777   write(iwr,101)icyc
101   format(/20x,26('*')/20x,'convergence at cycle',i6/20x,26('*'))
      irest=0
3333  return
      end
      subroutine david2(cr,iconf,c,z,iwr)
      implicit REAL  (a-h,o-z), integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension cr(*),iconf(5,*),c(*),z(*),q(2500)
c...
c...  davidson  diagonalizer -- 2 segment outstore mode
c...
      logical jeqi,odumpz
INCLUDE(common/sizes)
INCLUDE(common/helpr)
INCLUDE(common/dir_erg)
INCLUDE(common/diactl)
INCLUDE(common/auxh)
INCLUDE(common/cntrl)
INCLUDE(common/disktl)
      common/stoctl/khbill,khbz,khbc,khbb
      common/dialoc/floc(1275),f2loc(1275),xn(50),xns(50),ilifq(50)
      common/junk/igenz(4088),eig(50),g(16863)
INCLUDE(common/restar)
INCLUDE(common/timez)
INCLUDE(common/ccepa)
INCLUDE(common/mr_cepctl)
INCLUDE(common/mr_pairs)
      dimension fsloc(1275),f2sloc(1275)
      character *8 text,txcepa
      data text/'       '/,txcepa/'  cepa  '/
c     equivalence (q(1),g(1276))
      logical onorm
c we need some extra memory for kcpcc (cc) and kdavcz (cz)
      onorm=.false.
      if (icepa.ge.0)lvar=.true.
      if (onorm) then
      call setstr(50,1.0d0,xn)
      call setstr(50,1.0d0,xns)
      endif
      fact=2.0d0
      fact2=2.0d0
      var=0.0d0
      i26=index(26)
      i27=index(27)
      call getcz(iconf,c,i26)
      call getcz(iconf,z,i27)
      i28=index(28)
      i29=i28+ispsma
      index(29)=i29
      i199=index(199)
      call ibasgn(lcall,0,lcall,ilifq)
      xn(1)=1.0d0
      xns(1)=1.0d0
      ntot2=ntotal+ntotal
9999  floc(1)=h11
      old=cpulft(1)
      if(lvar)var=ddot(ntotal,z,1,z,1)
      f2loc(1)=var
      call wrt3(z,ntotal,i29,numci)
      call wrt3(c,ntotal,i28,numci)
ccepa
      if (icepa.ge.0) then
         if (moddia.gt.0) then
            call cntrci(iconf,c,z,0,0,cr(kcpcc+1),1,2,temp)
         else
            call vclr(cr(kcpcc+1),1,ncpcc*2)
            call cntrci(iconf,c,c,0,0,cr(kcpcc+1),0,1,temp)
         end if
         call dumpcc(cr(kcpcc+1),1,1)
c...  in case we come here after a space contraction ....
         if (ifcepa.eq.1) call hh2cep
     *     (cr(kcpcc+1),ncpcc,iconf,h11,h11,var,var,xn,1)
      end if
      fsloc(1) = floc(1)
      f2sloc(1) = f2loc(1)
      ncycep = 0
      h11o = h11
ccepa
      do 6666 joop=2,lcall
      jeqi=joop.eq.lcall
      odumpz=.true.
      ioop=joop-1
      h11sq=h11*h11
      if(lvar)var=var-h11sq
c... generate pert vector
      if (ifcepa.gt.0) then
         call cepvds(z,z,c,.false.,.true.,iconf)
      endif
ccepa  recompute var for cpea + emin or lock - modes
      if (ifcepa.gt.0.and.moddia.le.0) var = ddot(ntotal,z,1,z,1)
     *   -h11*h11
ccepa
_IF(cray)
      call daxpy(ntotal,-h11,c,1,z,1)
      call getcz(iconf,c,i199)
      if (ifcepa.gt.0) call cepvds(c,z,c,.true.,.false.,iconf)
      top=h11-eshift
      do 1 loop=1,ntotal
      temp=top-c(loop)
      temp=z(loop)/cvmgp(temp,dsign(0.05d0,temp),dabs(temp)-0.05d0)
      c(loop)=cvmgp(dsign(0.2d0,temp),temp,dabs(temp)-0.2d0)
    1 continue
      loop=idamax(ntotal,c,1)
      tester=dabs(c(loop))
_ELSEIF(ibm,vax)
      call gtrian(ntotal,h11,c,z,c)
      call getcz(iconf,z,i199)
      if (ifcepa.gt.0) call cepvds(z,z,c,.true.,.false.,iconf)
      call vsav(ntotal,eshift-h11,z,z)
      call mrgle(ntotal,0.05d0,z)
      call vvdv(ntotal,c,c,z)
      call mrgge(ntotal,0.2d0,c)
_IF(ibm)
      loop=idamax(ntotal,c,1)
      tester=dabs(c(loop))
_ELSE
      tester=0.0d0
      do 2 loop=1,ntotal
      temp=dabs(c(loop))
      if(temp.gt.tester)tester=temp
    2 continue
_ENDIF
_ELSE
      call daxpy(ntotal,-h11,c,1,z,1)
      call getcz(iconf,c,i199)
      if (ifcepa.gt.0) call cepvds(c,z,c,.true.,.false.,iconf)
      top=h11-eshift
      do 1 loop=1,ntotal
      temp=top-c(loop)
      if( dabs(temp).lt.(0.05d0))temp= dsign(0.05d0,temp)
      temp=z(loop)/temp
      if( dabs(temp).ge.(0.2d0))temp= dsign(0.2d0,temp)
    1 c(loop)=temp
      call maxmgv(c,1,tester,loop,ntotal)
_ENDIF
      top=cpulft(1)
      write(iwr,111)icyc,top,h11,tester,var,text
111   format(1x,i7,f12.3,3e20.11,a8)
      if(tester.lt.thresh)goto 777
      if(icyc.ge.maxcyc)goto 440
c...
c...  possible jacobi davidson preconditioning on vacuum space
c...  here khbb is the first-1 word available (khbz is waisted)
c
      call gdvdmr(index(1),index(199),i26,i27,
     +            khbb,nval,maxpa(4),
     +            c,z,h11,cr,iconf,numci,ncycg)
c
      if(malt)eshift=-eshift
c... orthogonalize pert vector
      x=ddot(ntotal,c,1,c,1)
9944  y=x
      joke=28
      do 609 koop=1,ioop
      call rdedx(z,ntotal,index(joke),numci)
      x=ddot(ntotal,c,1,z,1)*xns(koop)
      call daxpy(ntotal,-x,z,1,c,1)
609   joke=joke+2
      x=ddot(ntotal,c,1,c,1)
      if(x.lt.(0.05d0*y))goto 9944
      x=1.0d0/x
      y=dsqrt(x)
      if (onorm) then
         call dscal(ntotal,y,c,1)
         xn(joop)=1.0d0
         xns(joop)=1.0d0
      else
         xn(joop)=y
         xns(joop)=x
      endif
      if (joke.eq.30) then
         indc=index(joke)
      else
         indc=index(joke-1)+ispsma
      endif
      indz=indc+ispbig
      call dumpcz(iconf,c,indc)
      index(joke)=indc
      index(joke+1)=indz
      call genz(cr,iconf,index(joke),temp)
      call getcz(iconf,c,indc)
      call getcz(iconf,z,indz)
      m=iky(joop)
      mm=m+joop
c... compute h and h**2 elements
      if (onorm) then
        xold=x
        x=1.0d0
      endif
      floc(mm)=temp*x
      if(lvar)f2loc(mm)=ddot(ntotal,z,1,z,1)*x
      if (onorm)x=xold
ccepa
      if (icepa.ge.0) then
c...     c.c and c.z for current c/z vectors
         if (moddia.gt.0) then
            call cntrci(iconf,c,z,0,0,cr(kcpcc+1),1,2,temp)
         else
            call vclr(cr(kcpcc+1),1,ncpcc*2)
            call cntrci(iconf,c,c,0,0,cr(kcpcc+1),0,1,temp)
         end if
         call dumpcc(cr(kcpcc+1),joop,joop)
      end if
ccepa
      joke=29
      do 5097 koop=1,ioop
      call search(index(joke),numci)
      top=0.0d0
      bot=0.0d0
      loop=0
      if (icepa.ge.0) goto 5095
6098  moop=min(ntotal-loop,16863)
      call reads(g,moop,numci)
      top=ddot(moop,c(loop+1),1,g,1)+top
      if(lvar)bot=ddot(moop,z(loop+1),1,g,1)+bot
      loop=loop+16863
      if(loop.lt.ntotal)goto 6098
      goto 5096
5095  call rdedx(z,ntotal,index(joke),numci)
      top=ddot(ntotal,c,1,z,1)
      if(lvar) then
         call getcz(iconf,c,indz)
         bot=ddot(ntotal,c,1,z,1)
      endif
ccepa 
      if (icepa.ge.0) then
          call getcz(iconf,c,indc)
          call rdedx(z,ntotal,index(joke-1),numci)
          call cntrci(iconf,c,z,0,0,cr(kcpcc+1),0,1,temp)
          if (moddia.gt.0) then
            call rdedx(z,ntotal,index(joke),numci)
            call cntrci(iconf,c,z,0,0,cr(kcpcc+1),0,2,temp)
            call rdedx(c,ntotal,index(joke-1),numci)
            call getcz(iconf,z,indz)
            call cntrci(iconf,z,c,0,0,cr(kcpcc+1),0,3,top)
          else
            call vclr(cr(kcpcc+ncpcc+1),1,ncpcc*2)
          end if
          call dumpcc(cr(kcpcc+1),joop,koop)
      end if
      call getcz(iconf,c,indc)
      call getcz(iconf,z,indz)
ccepa
5096  temp=y*xn(koop)
      if (onorm)temp=1.0d0
      floc(m+koop)=temp*top
      f2loc(m+koop)=temp*bot
5097  joke=joke+2
ccepa
      call fmove(floc(m+1),fsloc(m+1),joop)
      call fmove(f2loc(m+1),f2sloc(m+1),joop)
      if (ifcepa.eq.1) go to 6000
ccepa
5098  if(moddia)6660,6661,6662
c... lock mode
6661  continue
      call dcopy(mm,floc,1,g,1)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
      do 6669 loop=1,joop
6669  eig(loop)=-dabs(q(ilifq(loop)+1))
      goto 6668
c... emin mode
6660  continue
      call dcopy(mm,floc,1,g,1)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
      goto 6668
c... vmin mode
 6662 continue
      call vsma(floc,1,-(h11+h11),f2loc,1,g,1,mm)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
_IF(cray,ibm,vax)
      call vsav(joop,h11sq,eig,eig)
_ELSE
      call vsadd(eig,1,h11sq,eig,1,joop)
_ENDIF
6668  continue
      ivrs=idmin(joop,eig,1)
      ivrs=ilifq(ivrs)
      if(jeqi)goto 4010
      if (odumpz) call wrt3(c,ntotal,indc,numci)
      if (odumpz) indz=indc+ispsma
      if (odumpz) call wrt3(z,ntotal,indz,numci)
      if (odumpz) index(joke)=indz
      if(lvar)var=expect(f2loc,q(ivrs+1),joop)
      h11=expect(floc,q(ivrs+1),joop)
c... generate best vector     and     h * best vector
4010  continue
      call vmul(q(ivrs+1),1,xn,1,eig,1,joop)
      if (icepa.lt.0) then
         call dscal(ntot2,eig(joop),z,1)
         joke=28
         do 621 koop=1,ioop
         call search(index(joke),numci)
         loop=0
         temp=eig(koop)
6088     moop=min(ntotal-loop,16863)
         call reads(g,moop,numci)
         call daxpy(moop,temp,g,1,c(loop+1),1)
         loop=loop+16863
         if(loop.lt.ntotal)goto 6088
         if(jeqi)goto 621
         call search(index(joke+1),numci)
         loop=0
6087     moop=min(ntotal-loop,16863)
         call reads(g,moop,numci)
         call daxpy(moop,temp,g,1,z(loop+1),1)
         loop=loop+16863
         if(loop.lt.ntotal)goto 6087
621      joke=joke+2
          call dumpcz(iconf,c,i26)
          call dumpcz(iconf,z,i27)
      else
         call dumpcz(iconf,z,i27)
         call dscal(ntotal,eig(joop),c,1)
         joke=28
         do 622 koop=1,ioop
         call search(index(joke),numci)
         loop=0
         temp=eig(koop)
6089     moop=min(ntotal-loop,16863)
         call reads(g,moop,numci)
         call daxpy(moop,temp,g,1,c(loop+1),1)
         loop=loop+16863
         if(loop.lt.ntotal)goto 6089
622      joke=joke+2
         call dumpcz(iconf,c,i26)
      endif
ccepa
ccepa see if cepa proces may start / if so modify h and h**2
ccepa and start the micro-iteration-process first
ccepa
      if (icepa.ge.0.and.tester.le.critcep.and.joop.ge.itcepa) ifcepa=1
      if (ifcepa.ne.1) go to 6001
      text = txcepa
c... control/print micro-iterations
         top=cpulft(1)
         bot=0.0d0
c        call givtim(top,bot)
         temp = h11-h11o
         if (lprcep.gt.0.and.ncycep.gt.0.and.
     *       dabs(temp).gt.tester*crycep) write(6,114)
     *                        ncycep,h11o,temp
114   format(' cepa micro it.',i3,
     *       ' last e ',e20.13,' delta ',e10.3,'  cepa ')
         temp = dabs(temp)
         h11o = h11
         if (ncycep.gt.0.and.
     *     (ncycep.ge.mcycep.or.temp.le.tester*crycep)) go to 6001
c...  compute pair energies
         call cepair(cr,iconf,i26)
c...  what was the ci-energy again ?
         h11ci = expect(fsloc,q(ivrs+1),joop)
c...  compute cepa shifts
         call cepac(c,iconf)
c restore c and z
         call rdedx(c,ntotal,indc,numci)
         call getcz(iconf,z,i27)
         odumpz=.false.
c...  modify h (floc) and h**2 (f2loc) matrices
         ncycep = ncycep + 1
6000     call hh2cep(cr(kcpcc+1),ncpcc,iconf,floc,fsloc,f2loc,
     *    f2sloc,xn,joop)
         go to 5098
ccepa
6001  ncycep = 0
      h11o = h11
ccepa
      if(jeqi)goto 43210
      if (icepa.ge.0) then
      call getcz(iconf,z,i27)
      call dscal(ntotal,eig(joop),z,1)
      joke=29
      do 623 koop=1,ioop
       call search(index(joke),numci)
       loop=0
       temp=eig(koop)
6090   moop=min(ntotal-loop,16863)
       call reads(g,moop,numci)
       call daxpy(moop,temp,g,1,z(loop+1),1)
       loop=loop+16863
       if(loop.lt.ntotal)goto 6090
623   joke=joke+2
      call dumpcz(iconf,z,i27)
      endif
       top=cpulft(1)
       tlefti=timlim-top-fact2
       if(tlefti.lt.fact*(top-old)) go to 440
       old=top
6666   continue
43210 continue
      temp=ddot(ntotal,c,1,c,1)
      call dscal(ntotal,dsqrt(1.0d0/temp),c,1)
      call dumpcz(iconf,c,i26)
      call genz(cr,iconf,index(26),h11)
       top=cpulft(1)
       write(iwr,1001)icyc,top
1001   format(i8,2f12.3,' space contracted')
      goto 9999
440   write(iwr,102)
102   format(//20x,14('*')/20x,'no convergence'/20x,14('*'))
      irest=9
      goto 3333
777   write(iwr,101)icyc
101   format(/20x,26('*')/20x,'convergence at cycle',i6/20x,26('*'))
      irest=0
3333  return
      end
      subroutine david1(cr,iconf,c,iwr)
      implicit REAL  (a-h,o-z), integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension cr(*),iconf(5,*),c(*),q(2500)
c...
c...  davidson  diagonalizer -- 1-segment out store mode
c...
      logical jeqi,nlvar
INCLUDE(common/sizes)
INCLUDE(common/helpr)
INCLUDE(common/dir_erg)
INCLUDE(common/cntrl)
INCLUDE(common/disktl)
INCLUDE(common/auxh)
      common/dialoc/floc(1275),f2loc(1275),xn(50),xns(50),ilifq(50)
INCLUDE(common/diactl)
      common/junk/igenz(4088),eig(50),g(16863)
      common/stoctl/khbill,khbz,khbc,khbb
INCLUDE(common/restar)
INCLUDE(common/timez)
INCLUDE(common/prnprn)
INCLUDE(common/ccepa)
INCLUDE(common/mr_cepctl)
INCLUDE(common/mr_pairs)
      dimension fsloc(1275),f2sloc(1275)
      character *8 text,txcepa
      data text/'       '/,txcepa/'  cepa  '/
      logical odumpz
c     equivalence (q(1),g(1276))
      if (icepa.ge.0)lvar=.true.
      fact=2.0d0
      fact2=2.0d0
      var=0.0d0
      i26=index(26)
      i27=index(27)
      i29=index(28)
      i28=i29+ispsma
      i199=index(199)
      index(29)=i29
      index(28)=i28
      call ibasgn(lcall,0,lcall,ilifq)
      xn(1)=1.0d0
      xns(1)=1.0d0
      call rdedx(c,ntotal,i27,numci)
      call wrt3(c,ntotal,i29,numci)
      if(lvar)var=ddot(ntotal,c,1,c,1)
      nlvar=.not.lvar
      call rdedx(c,ntotal,i26,numci)
      call wrt3(c,ntotal,i28,numci)
      if(odebug(40)) then
       call dbgvecv('david1-c',c,ntotal)
      endif
9999  floc(1)=h11
      old=cpulft(1)
      f2loc(1)=var
      j27=i29
c... c in buffer
ccepa
      if (icepa.ge.0) then
         if (moddia.gt.0) then
            call cntrci(iconf,c,cr(kcepc+1),i26,i27,cr(kcpcc+1),
     1       1,2,temp)
         else
            call vclr(cr(kcpcc+1),1,ncpcc*2)
            call cntrci(iconf,c,cr(kcepc+1),i26,i26,cr(kcpcc+1),
     1       0,1,temp)
         end if
         call dumpcc(cr(kcpcc+1),1,1)
c...  in case we come here after a space contraction ....
         if (ifcepa.eq.1) call hh2cep
     *     (cr(kcpcc+1),ncpcc,iconf,h11,h11,var,var,xn,1)
         call rdedx(c,ntotal,i26,numci)
      end if
      fsloc(1) = floc(1)
      f2sloc(1) = f2loc(1)
      ncycep = 0
      h11o = h11
ccepa
      do 6666 joop=2,lcall
      jeqi=joop.eq.lcall
      ioop=joop-1
      odumpz=.true.
      if(nlvar)goto 7009
      h11sq=h11*h11
      var=var-h11sq
      if (ifcepa.gt.0.and.moddia.le.0) var= -h11*h11
c... generate pert vector
      if (ifcepa.gt.0) 
     1 call cepvds_out(g,g,c,.false.,.false.,iconf,.true.,0)
7009  mb=0
      call search(j27,numci)
4000  call readsx(g,min(ntotal-mb,16353),nseg,numci)
      if (ifcepa.gt.0)
     1 call cepvds_out(g,g,c,.false.,.true.,iconf,.false.,
     2 nseg)
ccepa  recompute var for cpea + emin or lock - modes
      if (ifcepa.gt.0.and.moddia.le.0) var = 
     1  ddot(nseg,g,1,g,1) + var
ccepa
_IF(ibm,vax)
      call gtrian(nseg,h11,c(mb+1),g,c(mb+1))
_ELSE
      call vsmsb(c(mb+1),1,h11,g,1,c(mb+1),1,nseg)
_ENDIF
      mb=mb+nseg
      if(mb.lt.ntotal)goto 4000
      temp=eshift-h11
      if (ifcepa.gt.0) 
     1 call cepvds_out(g,g,c,.false.,.false.,iconf,.true.,0)
      mb=0
      call search(i199,numci)
4001  call readsx(g,min(ntotal-mb,16353),nseg,numci)
      if (ifcepa.gt.0)
     1 call cepvds_out(g,g,c,.true.,.false.,iconf,.false.,
     2 nseg)
_IF(ibm,vax)
      call vsav(nseg,temp,g,g)
_ELSE
      call vsadd(g,1,temp,g,1,nseg)
_ENDIF
      call mrgle(nseg,0.05d0,g)
_IF(ibm,vax)
      call vvdv(nseg,c(mb+1),c(mb+1),g)
_ELSE
      call vdiv(c(mb+1),1,g,1,c(mb+1),1,nseg)
_ENDIF
      mb=mb+nseg
      if(mb.lt.ntotal)goto 4001
      call mrgge(ntotal,0.2d0,c)
      loop=idamax(ntotal,c,1)
      tester=dabs(c(loop))
      top=cpulft(1)
      if(odebug(40)) then
       write(6,*)' david1: 111'
       call dbgvecv('david1-c',c,ntotal)
      endif
      write(iwr,111)icyc,top,h11,tester,var,text
111   format(1x,i7,f12.3,3e20.11,a8)
      if(tester.lt.thresh)goto 777
      if(icyc.ge.maxcyc)goto 440
c...
c...  possible jacobi davidson preconditioning on vacuum space
c...  here khbill is the first-1 word available and c starts there
c...  (and is ntotal long) 
c
      call gdvdmr(index(1),index(199),i26,i27,
     +            khbill+ntotal+nval,nval,maxpa(4),
     +            c,c(ntotal+1),h11,cr,iconf,numci,ncycg)
c
      if(odebug(40)) then
       write(6,*)' david1: after gdvdmr'
       call dbgvecv('david1-c',c,ntotal)
      endif
c
      if(malt)eshift=-eshift
c... orthogonalize pert vector
      x=ddot(ntotal,c,1,c,1)
9944  y=x
      joke=28
      do 609 koop=1,ioop
      x=0.0d0
      mb=0
      indj=index(joke)
      call search(indj,numci)
4002  nseg=min(ntotal-mb,16863)
      call reads(g,nseg,numci)
      x=ddot(nseg,c(mb+1),1,g,1)+x
      mb=mb+16863
      if(mb.lt.ntotal)goto 4002
      x=-x*xns(koop)
      mb=0
      call search(indj,numci)
4003  nseg=min(ntotal-mb,16863)
      call reads(g,nseg,numci)
      call daxpy(nseg,x,g,1,c(mb+1),1)
      mb=mb+16863
      if(mb.lt.ntotal)goto 4003
609   joke=joke+2
      if(odebug(40)) then
       write(6,*)' david1: 609'
       call dbgvecv('david1-c',c,ntotal)
      endif
      x=ddot(ntotal,c,1,c,1)
      if(x.lt.(0.05d0*y))goto 9944
      x=1.0d0/x
      xns(joop)=x
      y=dsqrt(x)
      xn(joop)=y
      indz=indj+ispsma
      indc=indz+ispbig
      index(joke+1)=indz
      index(joke)=indc
      call putcz(iconf,c,indc)
      call genz(cr,iconf,index(joke),temp)
      call rdedx(c,ntotal,indz,numci)
      if(odebug(40)) then
       write(6,*)' david1: indz'
       call dbgvecv('david1-c',c,ntotal)
      endif
c... z in buffer
      m=iky(joop)
      mm=m+joop
c... compute h and h**2 elements
      floc(mm)=temp*x
      if(lvar)f2loc(mm)=ddot(ntotal,c,1,c,1)*x
ccepa
      if (icepa.ge.0) then
c...     c.c and c.z for current c/z vectors
         if (moddia.gt.0) then
            call cntrci(iconf,c,cr(kcepc+1),indc,indz,cr(kcpcc+1)
     1 ,1,2,temp)
         else
            call vclr(cr(kcpcc+1),1,ncpcc*2)
            call cntrci(iconf,c,cr(kcepc+1),indc,indc,cr(kcpcc+1),
     1  0,1,temp)
         end if
         call dumpcc(cr(kcpcc+1),joop,joop)
         call rdedx(c,ntotal,indz,numci)
      end if
ccepa
      joke=29
      do 5097 koop=1,ioop
      temp=xn(koop)*y
      bot=0.0d0
      if(nlvar)goto 7000
      top=0.0d0
      call search(index(joke),numci)
      mb=0
6098  nseg=min(ntotal-mb,16863)
      call reads(g,nseg,numci)
      top=ddot(nseg,c(mb+1),1,g,1)+top
      mb=mb+16863
      if(mb.lt.ntotal)goto 6098
      f2loc(m+koop)=top*temp
      if(odebug(40)) then
       write(6,*)' david1: f2loc', m, koop, f2loc(m+koop)
      endif
7000  call search(index(joke-1),numci)
      mb=0
6099  nseg=min(ntotal-mb,16863)
      call reads(g,nseg,numci)
      bot=ddot(nseg,c(mb+1),1,g,1)+bot
      mb=mb+16863
      if(mb.lt.ntotal)goto 6099
      floc(m+koop)=bot*temp
      if(odebug(40)) then
       write(6,*)' david1: floc', m, koop, floc(m+koop)
      endif
ccepa
      if (icepa.ge.0) then
          call cntrci(iconf,c,cr(kcepc+1),indc,index(joke-1),
     1 cr(kcpcc+1),0,1,temp)
          if (moddia.gt.0) then
            call cntrci(iconf,c,cr(kcepc+1),indc,index(joke),
     1 cr(kcpcc+1),0,2,temp)
            call cntrci(iconf,c,cr(kcepc+1),indz,index(joke-1),
     1 cr(kcpcc+1),0,3,top)
          else
            call vclr(cr(kcpcc+ncpcc+1),1,ncpcc*2)
          end if
          call dumpcc(cr(kcpcc+1),joop,koop)
          call rdedx(c,ntotal,indz,numci)
      end if
ccepa
5097  joke=joke+2
ccepa
      call fmove(floc(m+1),fsloc(m+1),joop)
      call fmove(f2loc(m+1),f2sloc(m+1),joop)
      if (ifcepa.eq.1) go to 6000
ccepa
5098  if(moddia)6660,6661,6662
c... lock mode
6661  continue
      call dcopy(mm,floc,1,g,1)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
      do 6669 loop=1,joop
6669  eig(loop)=-dabs(q(ilifq(loop)+1))
      goto 6668
c... emin mode
6660  continue
      call dcopy(mm,floc,1,g,1)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
      goto 6668
c... vmin mode
 6662  continue
      call vsma(floc,1,-(h11+h11),f2loc,1,g,1,mm)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
_IF(cray,ibm,vax)
      call vsav(joop,h11sq,eig,eig)
_ELSE
      call vsadd(eig,1,h11sq,eig,1,joop)
_ENDIF
6668  continue
      ivrs=idmin(joop,eig,1)
      ivrs=ilifq(ivrs)
      call vmul(q(ivrs+1),1,xn,1,eig,1,joop)
      x=eig(joop)
      if(jeqi)goto 4010
      if (odumpz) call wrt3(c,ntotal,indz,numci)
      if(lvar)var=expect(f2loc,q(ivrs+1),joop)
      h11=expect(floc,q(ivrs+1),joop)
c... form h * best vector
      if (icepa.ge.0) goto 4010
      call dscal(ntotal,x,c,1)
      joke=29
      do 621 koop=1,ioop
      call search(index(joke),numci)
      mb=0
      temp=eig(koop)
6088  nseg=min(ntotal-mb,16863)
      call reads(g,nseg,numci)
      call daxpy(nseg,temp,g,1,c(mb+1),1)
      mb=mb+16863
      if(mb.lt.ntotal)goto 6088
621    joke=joke+2
      if(odebug(40)) then
       write(6,*)' david1: 621'
       call dbgvecv('david1-c',c,ntotal)
      endif
      call putcz(iconf,c,i27)
4010  call rdedx(c,ntotal,indc,numci)
      if(jeqi)goto 4020
      if (odumpz) indc=indz+ispsma
      if (odumpz) call wrt3(c,ntotal,indc,numci)
c... form best vector
4020  continue
      call dscal(ntotal,x,c,1)
      joke=28
      do 622 koop=1,ioop
      call search(index(joke),numci)
      mb=0
      temp=eig(koop)
6089  nseg=min(ntotal-mb,16863)
      call reads(g,nseg,numci)
      call daxpy(nseg,temp,g,1,c(mb+1),1)
      mb=mb+16863
      if(mb.lt.ntotal)goto 6089
622   joke=joke+2
      if(jeqi)goto 43210
      call putcz(iconf,c,i26)
      index(joke)=indc
      j27=i27
ccepa
ccepa see if cepa proces may start / if so modify h and h**2
ccepa and start the micro-iteration-process first
ccepa
      if (icepa.ge.0.and.tester.le.critcep.and.joop.ge.itcepa) ifcepa=1
      if (ifcepa.ne.1) go to 6001
      text = txcepa
c... control/print micro-iterations
         top=cpulft(1)
         bot=0.0d0
c        call givtim(top,bot)
         temp = h11-h11o
         if (lprcep.gt.0.and.ncycep.gt.0.and.
     *       dabs(temp).gt.tester*crycep) write(6,114)
     *                        ncycep,h11o,temp
114   format(' cepa micro it.',i3,
     *       ' last e ',e20.13,' delta ',e10.3,'  cepa ')
         temp = dabs(temp)
         h11o = h11
         if (ncycep.gt.0.and.
     *     (ncycep.ge.mcycep.or.temp.le.tester*crycep)) go to 6001
c...  compute pair energies
         call cepair(cr,iconf,i26)
c...  what was the ci-energy again ?
         h11ci = expect(fsloc,q(ivrs+1),joop)
c...  compute cepa shifts
         call cepac(c,iconf)
c restore c 
         call rdedx(c,ntotal,indc,numci)
         odumpz=.false.
c...  modify h (floc) and h**2 (f2loc) matrices
         ncycep = ncycep + 1
6000     call hh2cep(cr(kcpcc+1),ncpcc,iconf,floc,fsloc,f2loc,
     *    f2sloc,xn,joop)
         go to 5098
ccepa
ccepa
6001  ncycep = 0
      h11o = h11
ccepa
      if (jeqi)goto 43210
      if (icepa.ge.0) then
      call rdedx(c,ntotal,indz,numci)
      call dscal(ntotal,x,c,1)
      joke=29
      do 626 koop=1,ioop
      call search(index(joke),numci)
      mb=0
      temp=eig(koop)
6090  nseg=min(ntotal-mb,16863)
      call reads(g,nseg,numci)
      call daxpy(nseg,temp,g,1,c(mb+1),1)
      mb=mb+16863
      if(mb.lt.ntotal)goto 6090
626    joke=joke+2
      if(odebug(40)) then
       write(6,*)' david1: 621b'
       call dbgvecv('david1-c',c,ntotal)
      endif
      call putcz(iconf,c,i27)
      call rdedx(c,ntotal,i26,numci)
      endif
       top=cpulft(1)
       tlefti=timlim-top-fact2
       if(tlefti.lt.fact*(top-old)) go to 440
       old=top
6666   continue
43210 continue
      if(odebug(40)) then
       write(6,*)' david1: 43210'
       call dbgvecv('david1-c',c,ntotal)
      endif
      temp=ddot(ntotal,c,1,c,1)
      call dscal(ntotal,dsqrt(1.0d0/temp),c,1)
      call putcz(iconf,c,i26)
      call wrt3(c,ntotal,i28,numci)
      call genz(cr,iconf,index(26),h11)
      call rdedx(c,ntotal,i27,numci)
      call wrt3(c,ntotal,i29,numci)
      if(lvar)var=ddot(ntotal,c,1,c,1)
      call rdedx(c,ntotal,i28,numci)
      top=cpulft(1)
      write(iwr,1001)icyc,top
1001   format(i8,2f12.3,' space contracted')
      goto 9999
440   write(iwr,102)
102   format(//20x,14('*')/20x,'no convergence'/20x,14('*'))
      irest=9
      goto 3333
777   write(iwr,101)icyc
101   format(/20x,26('*')/20x,'convergence at cycle',i6/20x,26('*'))
      irest=0
3333  return
      end
      subroutine cciaib(cr,cs,b,epairs,iconf)
      implicit REAL  (a-h,o-z),integer (i-n)
      dimension cr(*),b(*),cs(*)
c...   singlet / vacuum (ia/ib) interactions   for cepa   **out-store**
c...         ** adapted to mr cepa only **
INCLUDE(common/symci)
INCLUDE(common/ccepa)
INCLUDE(common/ccntl)
INCLUDE(common/mr_cepctl)
INCLUDE(common/diactl)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh
_IF(cray,t3e,i8)
      integer iconf(5,*)
_ELSEIF(i8drct)
      integer*8 iconf(5,*)
_ELSE
      REAL iconf(5,*)
_ENDIF
      common/moco/ns,nt,lams,lamt
      common/loco/scr(202)
      common/spew/rmode(4)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      equivalence (rmode(4),ninti)
c
crwah mrdcepa and (jvl) cepa
c
      icl=ic
_IF(cray,i8,t3e,i8drct)
65    lamv=iconf(3,ir)
      itop = iconf(1,ir)
_ELSE
65    lamv=ipad(iconf(3,ir))
      itop = ipad(iconf(1,ir))
_ENDIF
      call getsb(b,lamv*lams)
      if (ir.gt.nref) go to 70
c... k(ic)
      call mxmd(cs,norbsi,1,cr(ninti+kh+1),1,1,scr,1,1,lams,norbsi,1)
c... produce singlet pair energy   contributions in z-array
      call mxmd(b,1,lamv,scr,1,1,cr(kcpz+1),1,1,lamv,lams,1)
      epairs = epairs + ddot(lamv,cr(kcpz+1),1,cr(kcpv+itop+1),1)
70    call getsb(rmode(1),1)
      if(model.ne.12) return
      call getsb(rmode(2),3)
      if(ic.eq.icl)goto 65
c
      return
      end
      subroutine cciajb(cr,cs,ct,b,epairs,epairt,iconf)
      implicit REAL  (a-h,o-z),integer (i-n)
      dimension cr(*),b(*),cs(*),ct(*)
c...   n-2 / vacuum (ia/jb) interactions  for closed and mr cepa
INCLUDE(common/ccepa)
INCLUDE(common/mr_cepctl)
INCLUDE(common/ccntl)
INCLUDE(common/diactl)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88
      common/moco/ns,nt,lams,lamt
      common/loco/scr(202)
      common/spew/rmode(4)
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      equivalence (rmode(4),ninti)
_IF(cray,t3e,i8)
      integer iconf(5,*)
_ELSEIF(i8drct)
      integer*8 iconf(5,*)
_ELSE
      REAL iconf(5,*)
_ENDIF
c
c... update z vacuum  => singlet/triplet pair - energies
c
      icl=ic
      lamst = lams + lamt
_IF(cray,i8,t3e,i8drct)
65    lamv=iconf(3,ir)
      itop = iconf(1,ir)
_ELSE
65    lamv=ipad(iconf(3,ir))
      itop = ipad(iconf(1,ir))
_ENDIF
      call getsb(b,lamv*lamst)
c...  skip if no reference config involved
      if (ir.gt.nref) go to 70
c
      ipbas=ninti+kh
      iqbas=ipbas+ns
      if (lams.le.0) go to 62
c****         ---  singlets  ---
c...  k(ic)
      call mxmd(cs,ns,1,cr(ipbas+1),1,1,scr,1,1,lams,ns,1)
c... update z vacuum  => singlet pair - energies
c     if (instor) kcpz=kcpz+itop
      call mxmd(b,1,lamv,scr,1,1,cr(kcpz+1),1,1,lamv,lams,1)
      epairs = epairs + ddot(lamv,cr(kcpz+1),1,cr(kcpv+itop+1),1)
62    if (lamt.le.0) go to 70
c****         ---  triplets  ---
c...  k(ic)
c     if (instor) kcpz=kcpz+itop
      call mxmd(ct,nt,1,cr(iqbas+1),1,1,scr,1,1,lamt,nt,1)
      call mxmd(b(lams*lamv+1),1,lamv,scr,1,1,cr(kcpz+1),1,1,
     *          lamv,lamt,1)
      epairt = epairt + ddot(lamv,cr(kcpz+1),1,cr(kcpv+itop+1),1)
c****
70    call getsb(rmode(1),1)
      if(model.ne.11) return
      call getsb(rmode(2),3)
      if(ic.eq.icl)goto 65
      return
      end
      subroutine cepinf (iconf,nword)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf,idc,mask,icon,pack2
_ELSEIF(convex,i8drct)
      integer *8 iconf,idc,mask,icon,ior81,ior84
      integer *8 pack2, shiftr
      external shiftr
_ELSEIF(bits8)
      REAL iconf,idc,mask,icon
_ELSE
      REAL iconf,idc,mask,icon,pack2
_ENDIF
c
c...  determine i,j identification of (n-2) states    (upack2-format)
c...   needed for closed shell cepa
c...   needed for    **mrd**   cepa  i,j is excitation class index
c...   needed for    mrcepa(1)  i,j is excitation class index+p/h identification
c...  conf(5,..) is used to store the information
ccepa **note** for both pair-energies and shifts goes ..
ccepa          i,j : singlet   / j,i  : triplet
ccepa          i,j is used for pair-energies and shifts (space maxpair*maxpair)
ccepa          since singlet-triplet distinction is sometimes needed i,j and j,i have meaning
ccepa          but not for MRCEPA's anymore
c
c     ** called from end of congen **
c
INCLUDE(common/ccepa)
INCLUDE(common/orb)
INCLUDE(common/cdcryz)
INCLUDE(common/ccntl)
      common/symcon/ nornor(26),ndimax
INCLUDE(common/cdcryi)
INCLUDE(common/iofile)
      common/maskc/mask(1)
INCLUDE(common/mr_mp2l)
INCLUDE(common/mr_pairs)
c
      dimension iconf(nword,*)
      dimension iorr(4),isneak(4)
_IF(cray,t3e,i8)
      integer msk
_ELSEIF(convex,i8drct)
      integer *8 msk
_ELSE
      REAL msk
_ENDIF
      dimension msk(3)
      data msk/3*0/
      integer popcnt
_IF(extpopcnt)
      external popcnt
_ENDIF
_IF(i8)
      pack2 (i1,i2)  = ior(ishft(i1,32),i2)
_ENDIF
_IF1()      integer shift
_IF1()      shift(ival,i) = ishftc(ival,i,-64)
c 
      if (icepa.ge.10) go to 100
c
c**** check if everything alright
c****  changed from 16/120 to 40/800
c
      if (nint.gt.maxpair.or.nnst.ne.1)
     * call caserr('unacceptable conditions for running CEPA012')
c
      icon = mask(nint*2)
c...  n-1 states
      do i=nstrt1,nlast1
      if (ndif(iconf(1,i),icon,iorr).ne.2) 
     1   call caserr('doublet state with ndif <> 2')
_IF(bits8)
         call pack2(0,iorr(1),iconf(5,i))
_ELSE
         iconf(5,i) = pack2(0,iorr(1))
_ENDIF
      end do
c...  n-2 states
      ndimax = 4
      do i=nstrt2,nlast2
         if (ndif(iconf(1,i),icon,iorr).ne.4)
     1   call caserr('n-2 state with ndif <> 4')
_IF(bits8)
         call pack2(iorr(4),iorr(1),iconf(5,i))
_ELSE
         iconf(5,i) = pack2(iorr(4),iorr(1))
_ENDIF
      end do
c
      return
c
c
c      for  all MRCEPA's the distinction between singlet/triplet is gone (also in cepair)
c
c      **mrd**  cepa(0)
c     i,j are :   1 / # holes + 1 in doc for vacuum
c                 2 / # holes + 1 in doc for nmin1
c                 3 / # holes + 1 in doc for nmin2  (singlet/triplet)
c      **mrd** cepa-1 
c     np = # partcles ; nh = # holes ; ib first doc hole (or 0) ; ia second doc hole (or 0)
c     mrdoc  = # doubly occuped orbitals (in reference function)
c     i,j are :   ip / jp 
c                 all                 jp = 1 + ib
c                 np=0,1,2,nh=0,1     ip = np+1
c                 np=0,1,2,nh=2       ip = 3+mrdoc*np+ia 
c     ndoc numbering is relative (done by routine docfind)
c
100   continue
c
c...  get double detection mask
c
      call set_docfind(iconf,msk,nref,mrdoc,'getmask')
c
      if (ndoc.gt.0.and.mrdoc.ne.ndoc) write(iwr,601) ndoc,mrdoc
601   format(/' ********  mrd cepa(0)  ******** ',/,
     1        '    ndoc given   ',i6,' ******** ',/,
     2        '   ndoc detected ',i6,' ******** ')
      if (ndoc.le.0) write(iwr,612)mrdoc
612   format(//,'  ** mrdcepa    # doubly occupied orbitals ',i3)
      if (lprcep.gt.0) write(iwr,613)
613   format('  ** double occupation pattern  ')
      if (lprcep.gt.0) call prcon(msk,0,0)
      ndoc = mrdoc
c
c     set code - info
c
      do i=1,nref
_IF(bits8)
         call pack2(1,1,iconf(5,i))
_ELSE
         iconf(5,i) = pack2(1,1)
_ENDIF
      end do
cmp    1,1 is better for the sneaky ones
c
c...  mr-cepa classes
c
      nref1 = nref + 1
c...   vacuum  **check on secret refs**
      nsec = 0
      do loop = nref1,nnst
         call docfind(iconf(1,loop),nnd,iorr)
         if (nnd.eq.0) then
            nsec = nsec + 1
            isneak(nsec) = loop
         end if
         if (nsec.eq.4) then 
            write(iwr,602) isneak
602         format(/' *** mrdcepa sneaky reference configs ',4i7)
            nsec = 0
         end if
c
         ip1 = 1
         ip2 = nnd+1
         if (icepa.ge.20.and.icepa.le.30) then
            ip2 = 1 + iorr(1)
            if (nnd.le.1) then 
               ip1 = 1
            else
               ip1 = 3+iorr(2)
            end if
         end if
_IF(bits8)
         call pack2(ip1,ip2,iconf(5,loop))
_ELSE
         iconf(5,loop) = pack2(ip1,ip2)
_ENDIF
      end do
      if (nsec.gt.0) write(iwr,602) (isneak(i),i=1,nsec)
c...   doublet
      do loop = nstrt1,nlast1
         call docfind(iconf(1,loop),nnd,iorr)
         ip1 = 2
         ip2 = nnd+1
         if (icepa.ge.20.and.icepa.le.30) then
            ip2 = 1 + iorr(1)
            if (nnd.le.1) then 
               ip1 = 2
            else
               ip1 = 3+mrdoc+iorr(2)
            end if
         end if
_IF(bits8)
         call pack2(ip1,ip2,iconf(5,loop))
_ELSE
         iconf(5,loop) = pack2(ip1,ip2)
_ENDIF
      end do
c...   singlet/triplet
      do loop = nstrt2,nlast2
         call docfind(iconf(1,loop),nnd,iorr)
         ip1 = 3
         ip2 = nnd+1
         if (icepa.ge.20.and.icepa.le.30) then
            ip2 = 1 + iorr(1)
            if (nnd.le.1) then 
               ip1 = 3
            else
               ip1 = 3+mrdoc*2+iorr(2)
            end if
         end if
_IF(bits8)
         call pack2(ip1,ip2,iconf(5,loop))
_ELSE
         iconf(5,loop) = pack2(ip1,ip2)
_ENDIF
      end do
c
      if ((maxpair.lt.3+mrdoc*3.and.icepa.ge.20.and.icepa.le.30)
     1    .or.lprcep.ge.1) then
         write(iwr,650) 
650      format(' MRCEPA pairs indices ',
     1        /,'     conf   - indices ')
         do loop=1,nlast2
            call upack2(iconf(5,loop),i,j)
            write(iwr,651) loop,i,j
651         format(1x,i10,2i5)
         end do
         if (maxpair.lt.3+mrdoc*3.and.icepa.ge.20.and.icepa.le.30)
     1      call caserr('maxpair in mr_pairs too small')
      end if
c
      return
      end
      subroutine docfind(conf,ndh,iorr)
c
c...  subroutine to find holes in doubly occupied space (at most 2)
c
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/natorb)
INCLUDE(common/orb)
INCLUDE(common/cdcryz)
       common/symcon/norb,norb1,norb16,nword,jspppp(22),ndimax  
_IF(cray,t3e,i8)
      integer conf,mask,msk,dub
_ELSEIF(convex,i8drct)
      integer *8 conf,mask,msk,dub
      integer *8 shiftr
      external shiftr
_ELSEIF(bits8)
      REAL conf,mask,msk,dub
_ELSE
      REAL conf,mask,msk,dub
_ENDIF
      dimension conf(5,*),msk(5),dub(5),mask(5),iorr(2),iorrr(4)
      character*(*) text
      integer popcnt
_IF(extpopcnt)
      external popcnt
_ENDIF
      save msk,doclist,ndoc
      integer doclist(62)
c
      dub(2)=0.0d0
      do i=1,nwm1
_IF(cray,t3e,i8,i8drct)
_IF(cray,t3e)
         dub(i) = and(conf(i,1),msk(i))
_ELSE
         dub(i) = iand(conf(i,1),msk(i))
_ENDIF
_ELSEIF(bits8)
         call dand8(conf(i,1),msk(i),dub(i))
_ELSE
         dub(i) = dand(conf(i,1),msk(i))
_ENDIF
      end do
c
      do i=1,4
         iorrr(i) = -1
      end do
      ndh =  ndif(msk,dub,iorrr) 
c
      if (ndh.eq.0) then
         iorr(1) = 0
         iorr(2) = 0
      else if(ndh.eq.1) then
         iorr(1) = locat1(doclist,ndoc,iorrr(1))
         if (iorr(1).eq.0) call caserr('docfind nd-1 error')
         iorr(2) = 0
      else if (ndh.eq.2) then
         ior1 = locat1(doclist,ndoc,iorrr(1))
         ior2 = locat1(doclist,ndoc,iorrr(4))
         if (ior1.eq.0.or.ior2.eq.0) call caserr('docfind nd-2 error')
         iorr(1) = min(ior1,ior2)
         iorr(2) = max(ior1,ior2)
      else
         call caserr('unexpected docfind usage')
      end if
c
      return
c
      entry set_docfind(conf,mask,nref,ndh,text)
c
      if (text.eq.'setmask') then
c
c...  set up # double detection mask
c
         ndimax_save = ndimax
         ndimax = 62
         nword = nwm1
         ndoc=0
         do i=1,nwm1
            msk(i) = n32m
            do j=1,nref
_IF(cray,t3e,i8,convex,i8drct)
               msk(i) = iand(msk(i),conf(i,j))
_ELSEIF(bits8)
               call dand8(msk(i),conf(i,j),msk(i))
_ELSE
               msk(i) = dand(msk(i),conf(i,j))
_ENDIF
            end do
            ndoc = ndoc + popcnt(msk(i))
         end do
c
         do i=1,nwm1
            mask(i) = 0
         end do
         ndh = ndif(msk,mask,doclist)
         if (ndh.ne.ndoc) call caserr('docfind foulup')
         if (ndh.gt.62) call caserr('docfind 62 foulup')
c
         do i=1,nwm1
_IF(cray,t3e,i8,convex,i8drct)
_IF(i8)
            msk(i) = ior(msk(i),ishft(msk(i),-1))
_ELSEIF(i8drct)
            msk(i) = ior(msk(i),shiftr(msk(i),1))
_ELSE
            msk(i) = or(msk(i),shiftr(msk(i),1))
_ENDIF
_ELSEIF(bits8)
            call shiftr8(msk(i),1,tmp8)
            call dor8(msk(i),tmp8,msk(i))
_ELSE
            msk(i) = dor(msk(i),shiftr(msk(i),1))
_ENDIF
         end do
         ndimax = ndimax_save
      end if
c
      do i=1,nwm1
         mask(i) = msk(i)
      end do
      ndh = ndoc
c
      return
      end
      subroutine cepddd(iconf,diag,z,c)
c
c...  adjust diagonal + z vector for doublets (ddc2)
c
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension diag(*),z(*),c(*),iconf(5,*)
c
INCLUDE(common/ccntl)
INCLUDE(common/symci)
c
INCLUDE(common/ccepa)
INCLUDE(common/mr_pairs)
c
      if (lsinsh.le.0) return
c
      kk = 0
      do 10 loop = nstrt1,nlast1
_IF(cray,t3e,convex,i8,i8drct)
         isi = iconf(4,loop)
_ELSE
         isi = ipad(iconf(4,loop))
_ENDIF
         nn = norbe(isi)
         if (icepa.lt.10) then
            bilbo = eshifd(loop-nstrt1+1)
         else
            call upack2(iconf(5,loop),ip1,ip2)
            bilbo = eshifc(ip1,ip2)
         end if
          do 5 i=1,nn
             diag(kk+i) = diag(kk+i) + bilbo
             z(kk+i) = z(kk+i) + c(kk+i)*bilbo
    5     continue
10    kk = kk + nn
c
      return
      end
      subroutine prcepa(iwr)
c
c...  print final cepa-results
c
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(common/ccepa)
INCLUDE(common/mr_pairs)
INCLUDE(common/orb)
c
      REAL ph_pair(3,3)
c
      ecorr = 0.0d0
      if (icepa.lt.10) then
        ndim=nint
      else
        ndim=3
      endif
      if (icepa.ge.20.and.icepa.lt.30) then
c...     translate to mrdcepa(0) pair energies
         do i=1,3
            do j=1,3
               ph_pair(j,i) = 0.0d0
            end do
         end do
         do i=1,3
           ph_pair(i,1) = epair(i,1)
           do j=1,mrdoc
             ph_pair(i,2) = ph_pair(i,2) + epair(i,j+1)
             do k=1,mrdoc
               ph_pair(i,3) = ph_pair(i,3) + epair(3+mrdoc*(i-1)+k,j+1)
             end do
           end do
           ecorr = ecorr + ph_pair(i,1) + ph_pair(i,2) + ph_pair(i,3)
         end do
      else
c
         do i=1,ndim
            do j=1,ndim
               ecorr = ecorr + epair(j,i)
            end do
         end do
      end if
c
      if (icepa.lt.10) write(iwr,600) icepa
      if (icepa.eq.10) write(iwr,608) 'final mrd-cepa   results '
      if (icepa.eq.11) write(iwr,608) 'final mr-cepa(1)-doc results '
      if (icepa.eq.12) write(iwr,608) 'final mr-cepa(1)   results '
      if (icepa.eq.20) write(iwr,608) 'final mr-cepa-1   results '
      if (icepa.eq.100) write(iwr,608) 'final mrcepa00  results ' 
      if (icepa.eq.101) write(iwr,608) 'final mr-aqcc    results ' 
      if (icepa.eq.102) write(iwr,608) 'final mr-acpf    results ' 
      if (lsinsh.gt.0.and.(icepa.lt.20.or.icepa.ge.30)) write(iwr,601)
      if (lcepaul.gt.0) write(iwr,602)
      write(6,610)
      write(iwr,607) 'from pair energies   ',ecorr
      if (icepa.ge.100) then
         if (icorr_mr.eq.1) write(iwr,607) 
     1      'from projection of ci',ecorr_mr
         if (icorr_mr.eq.2) write(iwr,607) 
     1      'from variation of ref',ecorr_mr
      end if
c
      if (icepa.ge.10) go to 5
      write(iwr,603)
      do 2 i=1,nint
2     write(iwr,604) (i,j,epair(i,j),j=1,i)
c
      write(iwr,605)
      do 3 i=2,nint
         ie = i-1
3     write(iwr,604) (i,j,epair(j,i),j=1,ie)
c
      goto 9
c
5     if (icepa.ge.20.and.icepa.lt.30) then
c      **mrd** cepa-1 
c     np = # particles ; nh = # holes ; ib first doc hole (or 0) ; ia second doc hole (or 0)
c     mrdoc  = # doubly occuped orbitals (in reference function)
c     i,j are :   ip / jp 
c                 all                 jp = 1 + ib
c                 np=0,1,2,nh=0,1     ip = np+1
c                 np=0,1,2,nh=2       ip = 3+mrdoc*np+ia 
c     ndoc numbering is relative (done by routine docfind)
c
         write(iwr,609) (((i-1),ph_pair(j,i),i=1,3),j=1,3)
         write(iwr,620)
         do ip=0,2
           do ih=0,2
             if (ih.eq.0) then
                write(iwr,621) ip,ih,' '
                write(iwr,622) epair(ip+1,1)
             else if (ih.eq.1) then
                write(iwr,621) ip,ih,'   --  holes =>  '
                write(iwr,623) (i,epair(ip+1,i+1),i=1,mrdoc)
             else if (ih.eq.2) then
               write(iwr,621) ip,ih,'   --  hole-hole '
               do j=1,mrdoc
                write(iwr,625) (j,i,epair(3+ip*mrdoc+j,i+1),i=1,j)
               end do
             end if
           end do
         end do
      else
c     *simple* mrcepa
         write(iwr,609) (((i-1),epair(j,i),i=1,3),j=1,3)
      end if
c
9     write(iwr,606)
c
600   format(///1x,80('*')//5x,'final cepa',i2,' results ')
601   format(5x,'(including single-shift)')
602   format(5x,'(paul s recipe)')
603   format(///5x,' singlet pair energies ')
604   format((2x,7(2i3,f12.8)))
605   format(//,5x,' triplet pair energies ')
606   format(/,1x,80('*'))
607   format(1x,4x,'total correlation energy ',a,f15.8)
608   format(///,1x,80(1h*),//,5x,a)
609   format(/  ,5x,'         pair energies  => holes ',/,
     1             5x,' vacuum  ',3(i3,f13.9),/,
     2             5x,'   n-1   ',3(i3,f13.9),/,
     3             5x,'   n-2   ',3(i3,f13.9) )
610   format(1x)
620   format(/ ,5x,'==========  detailed    pair energies ===========')
621   format(   2x,i3,' particle  ',i3,' hole ',a)
622   format(   5x,6x,f13.9)
623   format(  (5x,6(i3,3x,f13.9)))
625   format(  (5x,6(2i3,f13.9)))
c
      return
      end
      subroutine bucntl(iconf,cr,again,iwr)
c
c     control routine for reference configuration selection
c
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      logical again
      dimension iconf(5,*),cr(*)
INCLUDE(common/comrjb)
c
c     ifbuen : logical  /.true. do it
c     cradd(1) (0.1)  add config with excitation pattern iex1 (137b)
c     cradd(2) (0.04) add config with excitation pattern iex2 ( 13b)
c     cradd(3) (0.1)  keep reference config with old exc. pattern
c     iblkcc  : block-number for configs on ed7 (set in write after sele
c     iblkcr  : block-number for reference configs (after endinp)
c     iblkcv  : block-number for ci-start (initialise to 0) set in buenk
c     ncv     : number of start coef (initialize to 0 ) set in rjbw
c     nbuenk  : iteration count for rjb (initialise to 0)
c     mbuenk  : max. number of rjbs (e.g.initialise to 2)
c     weighb  : desired weigh of ref conf. (e.g.initialise .95)
c
INCLUDE(common/corctl)
INCLUDE(common/ccntl)
INCLUDE(common/cntrl)
INCLUDE(common/auxh)
INCLUDE(common/symci)
INCLUDE(common/table)
c
      again = .false.
      if (.not.ifbuen) return
      again = .true.
      nbuenk = nbuenk + 1
c
      maxs = 0
      maxt = 0
      do 10 i=1,nirr
      maxs = max(maxs,norbsi(i),norbe(i))
      maxt = max(maxt,norbtr(i))
10    continue
c
      ks = nlast2*5 + 1
      kcs = ks + maxs
      kct = kcs + maxt
      kselec = max(kcs+maxs*maxpa(1),kct+maxt*maxpa(2))
      kselec = max(kselec,ks+nval,kcs+ndoub)
      nuse = kselec + (nref+2)*2
      call corfai('buenker')
c
      call rjb(iconf,cr(ks),cr(kcs),cr(kcs),cr(kct),cr(kselec),
     1             newref,newcsf,rweig)
c
c...  move the selection data
c
      newre2 = newref*2
      do 20 i=1,newre2
20    cr(ks+i-1) = cr(kselec+i-1)
      kvc = ks + newre2
      kcs = kvc + newcsf
      kct = kcs + maxs*maxpa(1)
      nuse = kct + maxt*maxpa(2)
      call corfai('buenker')
c
      call rjbw(iconf,cr(ks),cr(kcs),cr(kct),cr(kvc),cr(kcs),
     1             newref,newcsf,again,rweig,iwr)
c
      return
      end
      subroutine rjb(iconf,s,cs,t,ct,selec,newref,newcsf,
     +               rweig)
c
c...  select most important configs to expand reference space
c
      implicit REAL  (a-h,o-z),integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension iconf(5,*),s(*),cs(*),t(*),ct(*),selec(2,*)
c
INCLUDE(common/sizes)
INCLUDE(common/comrjb)
INCLUDE(common/cntrl)
INCLUDE(common/ccntl)
INCLUDE(common/orb)
INCLUDE(common/disktl)
INCLUDE(common/table)
INCLUDE(common/symci)
INCLUDE(common/helpr)
c
      newref = 0
      newcsf = 0
      i26 = index(26)
      call search(i26,numci)
      rweig = 0.0d0
c
c...  vacuum
c
      if (nnst.eq.0) go to 50
      call reads(s,nval,numci)
      koop = 0
      do 40 loop=1,nnst
_IF(cray,t3e,convex,i8,i8drct)
         lams = iconf(3,loop)
_ELSE
         lams = ipad(iconf(3,loop))
_ENDIF
         hobbit = 0.0d0
         do 10 noop=1,lams
            koop = koop + 1
            hobbit = hobbit + s(koop)**2
10       continue
         if (loop.le.nref) rweig = rweig + hobbit
         if (hobbit.lt.cradd(2)) go to 40
c...  found a shining reference config
         newref = newref + 1
         newcsf = newcsf + lams
_IF(bits8)
         call pack3(loop,0,0,selec(1,newref))
_ELSE
         selec(1,newref) = pack3(loop,0,0)
_ENDIF
         selec(2,newref) = hobbit
40    continue
c
c...  doublets
c
50    if (nmin1.le.0) go to 100
      call reads(cs,ndoub,numci)
      koop = 1
      do 90 loop=nstrt1,nlast1
_IF(cray,t3e,convex,i8,i8drct)
         isi = iconf(4,loop)
_ELSE
         isi = ipad(iconf(4,loop))
_ENDIF
         nj = norbe(isi)
         mj = nbase(isi)
_IF(cray,t3e,convex,i8,i8drct)
         lams = iconf(3,loop)
_ELSE
         lams = ipad(iconf(3,loop))
_ENDIF
c...  sum squares of coefficients over can-set
         call avsq(s,cs(koop),nj,lams)
c...  scan the coefficient for each external orbital
         do 70 moop=1,nj
            if (s(moop).lt.cradd(2)) go to 70
            newref = newref + 1
            newcsf = newcsf + lams
_IF(bits8)
            call pack3(loop,0,moop+mj,selec(1,newref))
_ELSE
            selec(1,newref) = pack3(loop,0,moop+mj)
_ENDIF
            selec(2,newref) = s(moop)
70       continue
c
90    koop = koop + lams*nj
c
c...  singlets,triplets (all together)
c
100   if (nmin2.le.0) go to 200
c
      do 190 loop=nstrt2,nlast2
_IF(cray,t3e,convex,i8,i8drct)
         isi = iconf(4,loop)
_ELSE
         isi = ipad(iconf(4,loop))
_ENDIF
         call upack2(iconf(3,loop),lamt,lams)
         call reads(cs,lams*norbsi(isi),numci)
         call avsq(s,cs,norbsi(isi),lams)
         call reads(ct,lamt*norbtr(isi),numci)
         call avsq(t,ct,norbtr(isi),lamt)
c
c...  all summed now loop over irreps
c
         koos = 0
         koot = 0
         do 170 jsi=1,nirr
            nj = norbe(jsi)
            if (nj.eq.0) go to 170
            ksi = mult(jsi,isi)
            mj = nbase(jsi)
            if (ksi-jsi) 170,110,140
c...  c vectors totally symmetric (a1)
110         do 130 moop=1,nj
               jc = moop + mj
               do 120 joop=1,moop
                  koos = koos + 1
                  if (moop.eq.joop) go to 115
                  koot = koot + 1
                  s(koos) = s(koos) + t(koot)
115               if (s(koos).lt.cradd(2)) go to 120
c...  found a candidate
                  newref = newref + 1
                  newcsf = newcsf + lams
                  if (moop.ne.joop) newcsf = newcsf + lamt
_IF(bits8)
                  call pack3(loop,joop+mj,jc,selec(1,newref))
_ELSE
                  selec(1,newref) = pack3(loop,joop+mj,jc)
_ENDIF
                  selec(2,newref) = s(koos)
120            continue
130         continue
            go to 170
c...  c vectors not symmetric / straight
140         nk = norbe(ksi)
            if (nk.eq.0) go to 170
            mk = nbase(ksi)
            do 160 joop=1,nk
               kc = joop + mk
               do 150 moop=1,nj
                  koos = koos + 1
                  s(koos) = s(koos) + t(koos)
                  if (s(koos).lt.cradd(2)) go to 150
c...  found one
                  newref = newref + 1
                  newcsf = newcsf + lams + lamt
_IF(bits8)
                  call pack3(loop,moop+mj,kc,selec(1,newref))
_ELSE
                  selec(1,newref) = pack3(loop,moop+mj,kc)
_ENDIF
                  selec(2,newref) = s(koos)
150            continue
160         continue
170      continue
c
190   continue
c
200   continue
c
      return
      end
      subroutine avsq(r,a,nd,n)
c
c...  summ n squared vectors of dimension nd
c
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension r(*),a(*)
      call vclr(r,1,nd)
      if (n*nd.le.0) return
c
      ii=0
      do 20 i=1,n
         do 10 j=1,nd
10       r(j) = r(j) + a(j+ii)**2
20    ii=ii+nd
c
      return
      end
      subroutine rjbw(iconf,selec,cs,ct,vc,iconn,
     1                   newref,newcsf,again,rweig,iwr)
c
c     routine to write the new reference configurations and
c     the start for the next ci to ed7
c
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf,iconn
_ELSEIF(convex,i8drct)
      integer *8 iconf,iconn
_ELSE
      REAL iconf,iconn
      logical eq
_ENDIF
      logical again
      dimension iconf(5,*),iconn(5,*),selec(2,*),cs(*),ct(*),vc(*)
c
INCLUDE(common/comrjb)
c
INCLUDE(common/sizes)
INCLUDE(common/cntrl)
INCLUDE(common/ccntl)
INCLUDE(common/orb)
INCLUDE(common/cdcryi)
INCLUDE(common/disktl)
INCLUDE(common/symci)
INCLUDE(common/table)
INCLUDE(common/helpr)
INCLUDE(common/mapp)
      common /expanc/ ip(nd200),ipb(nd200),idum(nd200)
c
c...  check convergence
c
      nav=lenwrd()
      iwhy = 0
      kke = 0
      if (nbuenk.gt.mbuenk) iwhy = 2
      if (rweig.gt.weighb) iwhy = 1
c...  we will check if the reference space changed at all later
      if (iwhy.gt.0) go to 400
c
c
c...  set up the array of starting - ci vectors
c...  and set up the new reorder array to get all important
c...  external orbitals first
c
      i26 = index(26)
_IF(cray)
      call szero(ip,ninnex)
_ELSE
      call izero(ninnex,ip,1)
_ENDIF
      kkv = 0
      kke = nint
      jcon = 0
      do 200 i=1,newref
         call upack3(selec(1,i),icon,ia,ib)
c...  check if ia,ib should be added
         if (ia.le.1.or.locat1(ip,kke,ia).ne.0) go to 20
         kke = kke + 1
         ip(kke) = ia
20       if (ib.le.1.or.locat1(ip,kke,ib).ne.0) go to 30
         kke = kke + 1
         ip(kke) = ib
c...  adding of ci-vector elements to start-ci-vector
30       if (ib.eq.0) go to 50
         if (ia.eq.0) go to 100
         go to 150
c
c...  vacuum
c
50       if (jcon.gt.0) go to 60
c...  read vacuum coefficients
         call rdedx(cs,nval,i26,numci)
         jcon = icon
60     continue
_IF(cray,t3e,convex,i8,i8drct)
         kk = iconf(1,icon)
         lams = iconf(3,icon)
_ELSE
         kk = ipad(iconf(1,icon))
         lams = ipad(iconf(3,icon))
_ENDIF
         do 70 j=1,lams
            kkv = kkv + 1
            vc(kkv) = cs(kk + j)
70       continue
         go to 200
c...  doublets
100      if (jcon.ge.icon) go to 110
c...  read doublet coefficients
_IF(cray,t3e,convex,i8,i8drct)
         klaus=i26+iconf(2,1)
         call rdedx(cs,ndoub,klaus,numci)
_ELSE
         call rdedx(cs,ndoub,i26+ipad(iconf(2,1)),numci)
_ENDIF
         jcon = nstrt2 - 1
110      continue
_IF(cray,t3e,convex,i8,i8drct)
         isi = iconf(4,icon)
         kk = iconf(1,icon) + ib - nbase(isi)
_ELSE
         isi = ipad(iconf(4,icon))
         kk = ipad(iconf(1,icon)) + ib - nbase(isi)
_ENDIF
         nj = norbe(isi)
_IF(cray,t3e,convex,i8,i8drct)
         lams = iconf(3,icon)
_ELSE
         lams = ipad(iconf(3,icon))
_ENDIF
         do 120 j=1,lams
            kkv = kkv + 1
            vc(kkv) = cs(kk)
120      kk = kk + nj
         go to 200
c...  n-2 's  (singlet/triplet)
150      if (jcon.ge.icon) go to 160
c...  a new c
         jcon = icon
_IF(cray,t3e,convex,i8,i8drct)
         isi = iconf(4,icon)
_ELSE
         isi = ipad(iconf(4,icon))
_ENDIF
         call upack2(iconf(3,icon),lamt,lams)
         call upack2(iconf(2,icon),lbl,kbl)
         call rdedx(cs,lams*norbsi(isi),i26+kbl,numci)
         call rdedx(ct,lamt*norbtr(isi),i26+lbl,numci)
160      iroa = iro(ia)
         irob = iro(ib)
         ia = ia - nbase(iroa)
         ib = ib - nbase(irob)
         if (lams*norbsi(isi).le.0) go to 180
c     singlets
c...  symmetric case
         if (iroa.eq.irob)
     1   kk = ibassi(iroa,isi) + ib*(ib-1)/2  + ia
c...  asymmetric case
         if (iroa.ne.irob)
     1   kk = ibassi(iroa,isi) + (ib-1)*norbe(iroa) + ia
c...  get singlet coefficients
         do 170 j=1,lams
            kkv = kkv + 1
            vc(kkv) = cs(kk)
            kk = kk + norbsi(isi)
170      continue
180      if (lamt*norbtr(isi).le.0.or.(iroa.eq.irob.and.ia.eq.ib))
     *   go to 200
c     triplets
c...  symmetric case
         if (iroa.eq.irob)
     1   kk = ibastr(iroa,isi) + (ib-1)*(ib-2)/2  + ia
c...  asymmetric case
         if (iroa.ne.irob)
     1   kk = ibastr(iroa,isi) + (ib-1)*norbe(iroa) + ia
c...  get triplet coefficients
         do 190 j=1,lamt
            kkv = kkv + 1
            vc(kkv) = ct(kk)
            kk = kk + norbtr(isi)
190      continue
c
200   continue
c
c...  end of cv-vector generation
c
c...  reference configuration generation
c
_IF(i8drct)
      call readi8(iconf,nlast2*5,iblkcc,numci)
_ELSE
      call rdedx(iconf,nlast2*5,iblkcc,numci)
_ENDIF
c...  read reference confs + excitation patterns in
_IF(i8drct)
      call readi8(iconf,nref*5,iblkcr,numci)
_ELSE
      call rdedx(iconf,nref*5,iblkcr,numci)
_ENDIF
c...   clear not necessary parts of conf array  also externals
      call clrefc(iconf,5,nlast2,nint)
c
c...   now get new nint/next etc setup
c
      kke = kke - nint
      nint = nint + kke
      next = next - kke
      nintwo=nint+nint
      nint21=nintwo+1
      nintp1=nint+1
      nwm1=(nint-1)/n32+1
      nwp1=nwm1+1
c
      call coperm(ip,idum,ninnex,iwr)
      call pinver(ip,ipb,ninnex)
      call pinver(ip,ipb,ninnex)
c...  permutation in ip and invers one in ipb
c
      je = nw - 1
      do 300 i=1,newref
         call upack3(selec(1,i),icon,ia,ib)
c...  excitation pattern
         iconn(5,i) = iex2
         if (selec(2,i).gt.cradd(1)) iconn(5,i) = iex1
         if (selec(2,i).gt.cradd(3).and.icon.le.nref)
     1       iconn(5,i) = iconf(5,icon)
         iconn(1,i) = iconf(1,icon)
         iconn(2,i) = iconf(2,icon)
         iconn(3,i) = iconf(3,icon)
         iconn(4,i) = iconf(4,icon)
c....     note the following creations *must* succeed
         if (ia.gt.1) then
            if (icre(iconn(1,i),ipb(ia),iconn(1,i)).eq.0)
     1       call caserr('icre-1 error in rjb - zucht')
         end if
         if (ib.gt.1) then
            if (icre(iconn(1,i),ipb(ib),iconn(1,i)).eq.0)
     1       call caserr('icre-2 error in rjbw - zucht')
         end if
_IF(cray,t3e)
         call setsto(je-nwp1+1,0,iconn(nwp1,i))
_ELSE
         call setstz(je-nwp1+1,0,iconn(nwp1,i))
_ENDIF
300   continue
c
c...  check if anything happened
c
      if (kke.gt.0.or.newref.ne.nref) go to 330
c...  compare old and new reference set more carefully
      do 320 i=1,nref
         do 310 j=1,5
_IF(cray,t3e,convex,i8,i8drct)
            if(iconn(j,i).eq.iconf(j,i))go to 330
_ELSE
            if(eq(iconn(j,i),iconf(j,i)))go to 330
_ENDIF
310      continue
320   continue
c
      iwhy = 3
c
c...  cleanup
c
330   continue
c...  now get all new orbitals properly alligned
      call purmit(mapie,idum,ip,ninnex)
      call purmit(iro,idum,ip,ninnex)
c...  adapt mapei as well
      do 340 i=1,ninnex
340   mapei(mapie(i)) = i
c
c...  print
c
400   if (iwhy.gt.0) again = .false.
      sumsq = 0.0d0
      do 5 i=1,newref
   5  sumsq = sumsq + selec(2,i)
c
      write(iwr,601) nbuenk,nref,rweig,newref,sumsq,newcsf
601   format(//1x,54(1h*),
     1 /,' *      ----   reference generator ',i4,'  ----         *',
     1 /,' * old   ',i4,' configurations / weight ',f9.5,'       *',
     2 /,' * new   ',i4,' configurations / weight ',f9.5,'       *',
     3 /,' *       ',i4,' configuration-state-functions           *')
      ib = nint - kke + 1
      if (kke.gt.0) write(iwr,602) kke,(ip(i),i=ib,nint)
      if (kke.gt.0) write(iwr,608) (mapie(i),i=ib,nint)
602   format(' *',52x,'*',
     1 /,' *  the following ',i4,' external orbitals are added    *',
     2 /,(' *',8x,10i4,4x,'*'))
608   format(' *  in the original order they are known as           *',
     2 /,(' *',8x,10i4,4x,'*'))
c
      if (.not.again) write(iwr,604)
604   format(1x,54(1h*),
     1 /,' *    the ci-calculation is not reentered because     *')
      if (iwhy.eq.1) write(iwr,605)
      if (iwhy.eq.2) write(iwr,606)
      if (iwhy.eq.3) write(iwr,607)
605   format(
     1   ' *    the weight of the old set is sufficient         *',/,
     1         1x,54(1h*))
606   format(
     1   ' *    the maximum number of iterations is exceeded    *',/,
     2         1x,54(1h*))
607   format(
     1   ' *  the new reference-set is identical to the old one *',/,
     2         1x,54(1h*))
      write(iwr,603)
603   format(1x,54(1h*))
c
c...   no eigen or restrict allowed at this moment
c
      if (kke.gt.0) call inreo2(iconn,5,newref,iwr)
c
      if (.not.again) return
c
c...  write reference configs out
c
_IF(i8drct)
      call wrt3i8(iconn,5*newref,iblkcr,numci)
_ELSE
      call wrt3(iconn,5*newref,iblkcr,numci)
_ENDIF
c...  write start vector out
      iblkcv = iblkcr + lensec(5*newref)
      ncv = newcsf
      call wrt3(vc,ncv,iblkcv,numci)
c
      nxblk = iblkcv + lensec(ncv)
c
      nref = newref
c
      return
      end
      subroutine cepair(cr,iconf,iblkc)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf,idc
_ELSEIF(convex,i8drct)
      integer *8 iconf,idc
_ELSE
      REAL iconf,idc
_ENDIF
c
c... to obtain pair energies needed for cepa   for vector c1
c...   ** out store only (the rest of program can be instore)  **
c
      dimension cr(*),iconf(5,*)
      common/spcoef/nwwwww(9),iwwwww(9)
INCLUDE(common/sizes)
INCLUDE(common/timanb)
INCLUDE(common/disktl)
      common/stoctl/khbill,khbz,khbc,khbb,khbzd,khbzs,khbzt,khbcd,
     *khbcs,khbct,kh,khb1,khb2,khb3,khb4,khb5,khb6,khb7,khb8
     *,khb44,khb55,khb66,khb88
      common/spew/rmode(5)
INCLUDE(common/cdcryz)
      common/moco/ns,nt,lamsc,lamtc,lamsr,lamtr,lamss,lamtt,
     *khcsc,khzsc,khctc,khztc,isc,khcsr,khzsr,khctr,khztr,isr,
     *norsir,nortrr,lamstc
      common/symcon/root2,root2i
INCLUDE(common/symci)
INCLUDE(common/auxh)
INCLUDE(common/helpr)
INCLUDE(common/outctl)
INCLUDE(common/ccntl)
INCLUDE(common/cntrl)
INCLUDE(common/mr_cepctl)
INCLUDE(common/orb)
INCLUDE(common/ccepa)
INCLUDE(common/diactl)
c***
INCLUDE(common/mr_pairs)
      equivalence (nspnsp,nspips(1)),(nornor,norbsi(1))
      equivalence (rmode(1),model),(rmode(2),ic),(rmode(3),ir)
      equivalence (rmode(4),kkkbas),(rmode(5),icfock)
c***
c
c==================================================================
c
c...   clear relevant part of epair
c
      nepair = nint
      if (icepa.ge.10) nepair = 3
      if (icepa.ge.20.and.icepa.lt.30) nepair = maxpair
      do i=1,nepair
       do j=1,nepair
        epair(j,i) = 0.0d0
       end do
      end do
c
      nodel = 15
      tmlast = cpulft(1)
      call brtsb(index(1))
      call possb(nwwwww(9),iwwwww(9))
      call getsb(rmode(1),1)
c
      if (instor) goto 949
      kh = 0
c=================out of store version=====only=========================
92121 continue
      if (model.eq.11) go to 90011
      if (model.eq.12) go to 90012
      go to 99000
c...
c...  n-2 / vacuum (ia/jb)
c...
90011   call getsb(rmode(2),3)
        call rdedx(cr(khbill+1),m2exch,index(5),numci)
        kcpv=lump17
        kcpz=lump16
        call rdedx(cr(lump17+1),nval,iblkc,numci)
91011 continue
_IF(cray,t3e,convex,i8,i8drct)
      itop=iconf(4,ic)
_ELSE
      itop=ipad(iconf(4,ic))
_ENDIF
      ns=norbsi(itop)
      nt=norbtr(itop)
      call upack2(iconf(3,ic),lamtc,lamsc)
      call upack2(iconf(2,ic),kc,jc)
      nss=ns*lamsc
      ntt=nt*lamtc
      mus=nss+ntt
      call rdedx(cr(lump15+1),mus,iblkc+jc,numci)
      if(itop.eq.1)call renorm(cr(lump15+1),root2i,lamsc)
      call upack2(iconf(5,ic),ior1,ior2)
      if (icepa.lt.10) then
         call cciajb(cr,cr(lump15+1),cr(lump15+nss+1),cr(lump18+1),
     *               epair(ior1,ior2),epair(ior2,ior1),iconf)
      else
         esin = 0.0d0
         etri = 0.0d0
         call cciajb(cr,cr(lump15+1),cr(lump15+nss+1),cr(lump18+1),
     *               esin,etri,iconf)
         epair(ior1,ior2) = epair(ior1,ior2) + esin + etri
      end if
        if(model.eq.11)goto 91011
        goto 92121
c...
c...  n-2  /  vacuum  (ia/ib)
c...
90012 call getsb(rmode(2),3)
      call rdedx(cr(khbill+1),m2fock,index(4),numci)
      kcpv=lump22
      kcpz=lump21
      call rdedx(cr(lump22+1),nval,iblkc,numci)
91012 continue
_IF(cray,t3e,convex,i8,i8drct)
      itop=iconf(4,ic)
      lamsc=iand(iconf(3,ic),mms)
_ELSEIF(bits8)
      itop=ipad(iconf(4,ic))
      call dand8(iconf(3,ic),mms,tmp8)
      lamsc = ipad(tmp8)
_ELSE
      itop=ipad(iconf(4,ic))
      lamsc=ipad(dand(iconf(3,ic),mms))
_ENDIF
      nss=norbsi(itop)*lamsc
_IF(cray,t3e,convex,i8,i8drct)
      jc=iand(iconf(2,ic),mms)
_ELSEIF(bits8)
      call dand8(iconf(2,ic),mms,tmp8)
      jc=ipad(tmp8)
_ELSE
      jc=ipad(dand(iconf(2,ic),mms))
_ENDIF
      call rdedx(cr(lump20+1),nss,iblkc+jc,numci)
      if(itop.eq.1)call renorm(cr(lump20+1),root2i,lamsc)
      call upack2(iconf(5,ic),ior1,ior2)
      call cciaib(cr,cr(lump20+1),cr(lump23+1),epair(ior1,ior2),iconf)
      if(model.eq.12)goto 91012
      goto 92121
c
99000 if (icepa.lt.10) goto 99999
c...
c... vacuum / vacuum (ij/kl)        **mrdcepa**
c...
      kcpill=lump2
      kcpv=lump1
      kcpz=khbill
      call rdedx(cr(lump1+1),nval,iblkc,numci)
      call brtsb(indcep)
      call getsb(rmode(1),1)
      call ccijkl(iconf,cr)
c...
c... doublet / vacuum (ij/ka)       **mrdcepa**
c...
      if(ndoub.gt.0) then
        call possb(nwwwww(3),iwwwww(3))
        call getsb(rmode(1),1)
        kcpd=lump8
        kcpz=lump5
        kcpv=lump7
        kcpill=khbill
        khbb=khbill
        kcpcd=lump6
        call rdedx(cr(lump7+1),nval,iblkc,numci)
_IF(cray,i8,t3e,i8drct)
        call rdedx(cr(kcpd+1),ndoub,iblkc+iconf(2,1),numci)
_ELSE
        call rdedx(cr(kcpd+1),ndoub,iblkc+ipad(iconf(2,1)),numci)
_ENDIF
        call ccijka(iconf,cr)
      endif
c
      kcpv=lump17
      call rdedx(cr(lump17+1),nval,iblkc,numci)
c
      go to 99999
c
c=================mrd cepa in store version==============================
c
949   call renorm(cr(khbcs+1),root2i,nspnsp)
950   continue
      kh = khbb - khbill
      if (model.eq.11) go to 9011
      if (model.eq.12) go to 9012
      go to 9900
c...
c...  n-2 / vacuum (ia/jb)
c...
9011   call getsb(rmode(2),3)
       call rdedx(cr(khbb+1),m2exch,index(5),numci)
9111  continue
_IF(cray,t3e,convex,i8,i8drct)
      itop=iconf(4,ic)
_ELSE
      itop=ipad(iconf(4,ic))
_ENDIF
      ns=norbsi(itop)
      nt=norbtr(itop)
      call upack2(iconf(3,ic),lamtc,lamsc)
      call upack2(iconf(1,ic),kc,jc)
      nss=ns*lamsc
      ntt=nt*lamtc
      mus=nss+ntt
      call upack2(iconf(5,ic),ior1,ior2)
      kcpz=khbz
      kcpv=khbc
      if (icepa.lt.10) then
         call cciajb(cr,cr(khbcs+jc+1),cr(khbct+kc+1),cr(khb1+1),
     *               epair(ior1,ior2),epair(ior2,ior1),iconf)
      else
         esin = 0.0d0
         etri = 0.0d0
         call cciajb(cr,cr(khbcs+jc+1),cr(khbct+kc+1),cr(khb1+1),
     *               esin,etri,iconf)
         epair(ior1,ior2) = epair(ior1,ior2) + esin + etri
      end if
      if(model.eq.11)goto 9111
      goto 950
c...
c...  n-2  /  vacuum  (ia/ib)
c...
9012  call getsb(rmode(2),3)
      call rdedx(cr(khbb+1),m2fock,index(4),numci)
9112  continue
_IF(cray,t3e,convex,i8,i8drct)
      itop=iconf(4,ic)
_IF(cray,t3e)
      lamsc=and(iconf(3,ic),mms)
_ELSE
      lamsc=iand(iconf(3,ic),mms)
_ENDIF
_ELSEIF(bits8)
      itop=ipad(iconf(4,ic))
      call dand8(iconf(3,ic),mms,tmp8)
      lamsc = ipad(tmp8)
_ELSE
      itop=ipad(iconf(4,ic))
      lamsc=ipad(dand(iconf(3,ic),mms))
_ENDIF
      nss=norbsi(itop)*lamsc
_IF(cray,t3e,convex,i8,i8drct)
_IF(cray,t3e)
      jc=and(iconf(1,ic),mms)
_ELSE
      jc=iand(iconf(1,ic),mms)
_ENDIF
_ELSEIF(bits8)
      call dand8(iconf(1,ic),mms,tmp8)
      jc=ipad(tmp8)
_ELSE
      jc=ipad(dand(iconf(1,ic),mms))
_ENDIF
         call upack2(iconf(5,ic),ior1,ior2)
        kcpz=khbz
        kcpv=khbc
         call cciaib(cr,cr(khbcs+jc+1),cr(khb2+1),epair(ior1,ior2),
     *               iconf)
      if(model.eq.12)goto 9112
      goto 950
c
9900  if (icepa.lt.10) goto 99998
c...
c... vacuum / vacuum (ij/kl)        **mrdcepa**
c...
      call brtsb(indcep)
      call getsb(rmode(1),1)
      kcpill=khbb
      kcpv=khbc
      kcpz=khbz
      call ccijkl(iconf,cr)
c...
c... doublet / vacuum (ij/ka)       **mrdcepa**
c...
      if(ndoub.gt.0) then
        call possb(nwwwww(3),iwwwww(3))
        call getsb(rmode(1),1)
        kcpill=khbb
        kcpd=khbcd
        kcpv=khbc
        call ccijka(iconf,cr)
      endif
c
99998 call renorm(cr(khbcs+1),root2,nspnsp)
      kcpv=khbc
c
99999 c0 = ddot(nrfcsf,cr(kcpv+1),1,cr(kcpv+1),1)
      scale = 1.0d0/c0
c...    multiply by 1/c0**2 to compensate for mxmd's with ref-vector
c...    thus we scale back to c0 = 1
      do i=1,nepair
       do j=1,nepair
        epair(j,i) = epair(j,i)*scale
       end do
      end do
      c0 = dsqrt(c0)
c
      return
      end
_ENDIF(not,parallel)
_IF(parallel)
      subroutine besips(joop,eig,c,z,g,nbuf,icurc,icurz,iwr)
c
c...  generate best c and z vector for davidi
c
      implicit REAL  (a-h,o-z)
c
      dimension eig(joop)
      dimension c(nbuf),z(nbuf),g(nbuf)
c MPP
INCLUDE(common/prnprn)
INCLUDE(common/disktl)
      common/noddav/nwips,iblock1,iblock2,nwords,nread
c MPP
      if(odebug(40)) write(iwr,*)' entered besips'
      if(nwips.gt.0) then
      ntime = (nwips-1)/nbuf + 1
      nblock = nbuf/511
      nnn = nbuf
      do 100 itime = 1,ntime
         ioff = (itime-1)*nblock
         if (itime.eq.ntime) nnn = nwips - ioff*511
         call rdedx(c,nnn,index(26+2*joop)+ioff,numcip)
         call rdedx(z,nnn,index(27+2*joop)+ioff,numcip)
         call dscal(nnn,eig(joop),c,1)
         call dscal(nnn,eig(joop),z,1)
         do 10 i=1,joop-1
           call rdedx(g,nnn,index(26+2*i)+ioff,numcip)
           call daxpy(nnn,eig(i),g,1,c,1)
           call rdedx(g,nnn,index(27+2*i)+ioff,numcip)
           call daxpy(nnn,eig(i),g,1,z,1)
10       continue
         call wrt3(c,nnn,icurc+ioff,numcip)
         call wrt3(z,nnn,icurz+ioff,numcip)
100   continue
      endif
c
       if(odebug(40)) write(iwr,*)' quit besips'
       return
       end
_IF(ga)
      subroutine clrcz(iconf,iblock,iwr,numcin)
_ELSE
      subroutine clrcz(iconf,iblock,iwr)
_ENDIF
c
c...  clear a normal c/z vector starting at iblock
c      kind of dumpcz action
c
      implicit REAL  (a-h,o-z)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension iconf(5,*)
INCLUDE(common/prnprn)
       parameter (nbuf=511*10)
      common /scra/ buf(nbuf)
INCLUDE(common/disktl)
c
INCLUDE(common/cntrl)
INCLUDE(common/ccntl)
INCLUDE(common/symci)
INCLUDE(common/sizes)
INCLUDE(common/disc)
c
      if(odebug(40)) write(iwr,*)' entered clrcz'
      nwrit = 0
      loop = nstrt2-1
      istage = 2
_IF(ga)
      call switch_oswed(numcin,.true.)
      call search(iblock,numcin)
_ELSE
      call search(iblock,numcip)
_ENDIF
      naim = nval
      call vclr(buf,1,nbuf)
c
c...  dump zero's to numcip in nice dumpcz manner
c
c
c..      N-2  (singlet/triplet vectors nicely alligned)
c...     istage = 1 : busy on singlet
c...              2 : busy on triplet
c...     doublets/valence  (one long vector) are treated as 1 triplet
c...     naim = current aim of # elements to write
c
10    if (nwrit.lt.naim) then
         nw = min(nbuf,naim-nwrit)
_IF(ga)
         oadd(numcz)=.false.
         call wrt3s(buf,nw,numcin)
         oadd(numcz)=.true.
_ELSE
         call wrt3s(buf,nw,numcip)
_ENDIF
         nwrit = nwrit + nw
         go to 10
      end if
c
c...  check for finish
c
      if (nwrit.eq.ntotal) then
          if (loop.ne.nlast2.or.istage.ne.2) call caserr('clrcz1')
          if(odebug(40)) write(iwr,*)' quit clrcz'
          return
      else if (nwrit.gt.ntotal) then
         call caserr('clrcz2?')
      end if
c
c..   prepare for next vector
c
      if (nwrit.lt.nval+ndoub) then
         naim = nval+ndoub
      else if (istage.eq.2) then
         loop = loop + 1
_IF(cray,t3e,convex,i8,i8drct)
         isi=iconf(4,loop)
_ELSE
         isi=ipad(iconf(4,loop))
_ENDIF
         call upack2(iconf(3,loop),lamt,lams)
         lams=lams*norbsi(isi)
         lamt=lamt*norbtr(isi)
         naim = naim + lams
         istage = 1
         if (lamt.eq.0) istage=2
      else if (istage.eq.1) then
         naim = naim + lamt
         istage = 2
      end if
c
      go to 10
c
      call switch_oswed(numcin,.false.)
      return
      end
      subroutine davidi(cr,iconf,c,z,gg,nbuf,creal,zreal,
     +                     inipsc,iwr)
c MPP
c MPP    this one will work out store distributed for ipsc
c MPP    all dot-product etc are changed
c MPP    so not everything is flagged
c MPP    derived from david3
c MPP    since for the davidson the c/z vectors are divided 
c MPP    among the nodes  / this c/z is a lot smaller
c MPP    if the c/z vectors fit inipsc = .true.
c MPP    creal,zreal are the real (instore) c/z vectors
c MPP
      implicit REAL  (a-h,o-z), integer (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
      dimension cr(*),iconf(5,*),creal(*),zreal(*)
      dimension c(nbuf),z(nbuf),gg(nbuf),q(2500)
      logical inipsc
INCLUDE(common/prnprn)
c...
c...  davidson  diagonalizer 
c...
      logical jeqi
INCLUDE(common/sizes)
INCLUDE(common/helpr)
INCLUDE(common/dir_erg)
INCLUDE(common/diactl)
INCLUDE(common/cntrl)
INCLUDE(common/disktl)
      common/dialoc/floc(1275),f2loc(1275),xn(50),xns(50),ilifq(50)
      common/junk/igenz(4088),eig(50),g(16863)
INCLUDE(common/restar)
INCLUDE(common/timez)
      common/remc_ed19/nmnmnmn,ipipip
      equivalence (q(1),g(1276))
      if(odebug(40)) write(iwr,*)' entered davidi'
      fact=2.0d0
      fact2=2.0d0
      var=0.0d0
         call ibasgn(lcall,0,lcall,ilifq)
      i26=index(26)
      i27=index(27)
      i28=index(28)
      i29=index(29)
      i199=index(199)
          idump = index(28+lcall+lcall)
          icurrc = index(28+lcall+lcall+1)
          icurrz = index(28+lcall+lcall+2)
      xn(1)=1.0d0
      xns(1)=1.0d0
9999     continue
c MPP
        icblok=0
        call gczips(index(26),i28,.true.,iwr,cr,icblok)
        call gczips(index(27),i29,.true.,iwr,cr,icblok)
        call copy_bloks(i28,icurrc,icblok)
        call copy_bloks(i29,icurrz,icblok)
c       call gczips(index(26),icurrc,.false.,iwr,iconf)
c       call gczips(index(27),icurrz,.true.,iwr,iconf)
        call pg_synch(1111)
        h11 = dotips(i28,i29,iwr)
        if(odebug(40)) write(iwr,*)' dotips h11 ',h11
c MPP
      floc(1)=h11
      old=cpulft(1)
       if (lvar) var = dotips(i29,i29,iwr)
       f2loc(1)=var
      do 6666 joop=2,lcall
      jeqi=joop.eq.lcall
      ioop=joop-1
      h11sq=h11*h11
      if(lvar)var=var-h11sq
c... generate pert vector
      call perips(c,z,nbuf,icurrc,icurrz,idump,
     + eshift,h11,tester,x,iwr)
c MPP
      top=cpulft(1)
      write(iwr,111)icyc,top,h11,tester,var
111   format(1x,i7,f12.3,3e20.11)
      call flush(iwr)
      if(tester.lt.thresh)goto 777
      if(icyc.ge.maxcyc)goto 440
      if(malt)eshift=-eshift
c
c... orthogonalize pert vector
c
       call ortips(c,z,nbuf,inipsc,ioop,idump,x,gg,xns,iwr)
c
      x=1.0d0/x
      xns(joop)=x
      y=dsqrt(x)
      xn(joop)=y
c MPP
c MPP   dump c-vector found to index(26)
c MPP
      call dczips(cr,index(26),idump,iwr)
c MPP    clear z-vector
      call clrcz(iconf,index(27),iwr,numcz)
      if (instor) then
         call getcz(iconf,creal,i26)
         call getcz(iconf,zreal,i27)
      end if
c MPP
      call switch_oswed(numcz,.false.)
      call genz(cr,iconf,index(26))
      call pg_synch(1112)
      call switch_oswed(numcz,.true.)
c MPP
      if (instor) then
         call dumpcz(iconf,creal,i26)
         call dumpcz(iconf,zreal,i27)
      end if
c MPP
c MPP   dump vectors
c MPP
        call gczips(index(26),index(26+joop*2),.false.,iwr,cr,icbbb)
        call gczips(index(27),index(27+joop*2),.false.,iwr,cr,icbbb)
c       call gczips(index(26),index(26+joop*2),.false.,iwr,iconf)
c       call gczips(index(27),index(27+joop*2),.true.,iwr,iconf)
c MPP
c MPP   calc h1n,h2n,h1n, etc and z1n,z2n, ... as efficient as possible
c MPP
c      if (joop.eq.3) call debugj(xn,joop)
      m=iky(joop)
      mm=m+joop
c
       call hhhips(joop,lvar,floc(m+1),f2loc(m+1),xn,c,z,gg,nbuf,iwr)
       if(odebug(40)) then
        write(iwr,*)' floc ',(floc(i),i=1,joop*(joop+1)/2)
       endif
c
c..   next bit might be done by 1 node ..
c
      if(moddia)6660,6661,6662
c... lock mode
6661  continue
      call dcopy(mm,floc,1,g,1)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
      do 6669 loop=1,joop
6669  eig(loop)=-dabs(q(ilifq(loop)+1))
      goto 6668
c... emin mode
6660  continue
      call dcopy(mm,floc,1,g,1)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
      goto 6668
c... vmin mode
6662  continue
       call vsma(floc,1,-(h11+h11),f2loc,1,g,1,mm)
      call jacobi(g,iky,joop,q,ilifq,joop,eig,2,1,1.0d-9)
      call vsadd(eig,1,h11sq,eig,1,joop)
6668  continue
      ivrs=idmin(joop,eig,1)
      ivrs=ilifq(ivrs)
      if(jeqi)goto 4010
      if(lvar)var=expect(f2loc,q(ivrs+1),joop)
      h11=expect(floc,q(ivrs+1),joop)
c... generate best vector     and     h * best vector
4010  continue
      call vmul(q(ivrs+1),1,xn,1,eig,1,joop)
      call besips(joop,eig,c,z,gg,nbuf,icurrc,icurrz,iwr)
c MPP
       if(jeqi)goto 43210
       top=cpulft(1)
       tlefti=timlim-top-fact2
       if(tlefti.lt.fact*(top-old)) go to 440
       old=top
6666   continue
43210 continue
c MPP
c MPP   folding (keeping 2 vectors would be better)
      temp=dotips(icurrc,icurrc,iwr)
      call sclips(c,nbuf,dsqrt(1.0d0/temp),icurrc,iwr)
      call dczips(cr,i26,icurrc,iwr)
      call clrcz(iconf,i27,iwr,numcz)
      if (instor) then
         call getcz(iconf,creal,i26)
         call getcz(iconf,zreal,i27)
      end if
c MPP
      call switch_oswed(numcz,.false.)
      call genz(cr,iconf,index(26))
      call pg_synch(789)
      call switch_oswed(numcz,.true.)
c MPP
      if (instor) then
         call dumpcz(iconf,creal,i26)
         call dumpcz(iconf,zreal,i27)
      end if
c MPP
       top=cpulft(1)
       write(iwr,1001)icyc,top
1001   format(i8,2f12.3,' space contracted')
      goto 9999
440   write(iwr,102)
102   format(//20x,14('*')/20x,'no convergence'/20x,14('*'))
      irest=9
      goto 3333
777   write(iwr,101)icyc
101   format(/20x,26('*')/20x,'convergence at cycle',i6/20x,26('*'))
      call walltime(ddtime)
      write (iwr,6010)ddtime
 6010 format (/1x,40('+')/1x,'walltime in seconds ',F20.2,
     +/1x,40('+'))
      irest=0
3333  continue
c MPP
c     if (opg_root()) then
         call dczips(iconf,i26,icurrc,iwr)
         call dczips(iconf,i27,icurrz,iwr)
         call clrcz(iconf,i27,iwr,numcz)
         if (instor) then
            call getcz(iconf,creal,i26)
            call getcz(iconf,zreal,i27)
         end if
c     end if
c MPP
      return
      end
      subroutine davips(cr,icr,creal,zreal,iwr)
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
      dimension cr(*),icr(*),creal(*),zreal(*)
c MPP
c MPP.    sort out proper usage of private numci depending on davidson used
c MPP     all relies on 5-record chunks
c MPP     this redoes part of the diabas work
c MPP     then call davidi
c MPP
INCLUDE(common/prnprn)
INCLUDE(common/disktl)
      common/nodinf/mpid,minode,mihost,ldim,nnodes,nodscr,
     +  maxnod(maxlfn)
      common/noddav/nwips,iblock1,iblock2,nwords,nread
INCLUDE(common/cntrl)
INCLUDE(common/diactl)
      common/stoctl/ntotw,khbz,khbc,khbb
INCLUDE(common/corctl)
       logical inipsc
       common/remc_ed19/niblok,ipoint,idiagoff
c
c     index 26/27 are for complete vectors (size ispbig)
c     ispsma = # blocks for vectors as written by davidsons
c     calculate new ispsma estimate (we do it in 5-record chunks) 
c
      nwips = 0
      ispsma = ispsma/nnodes
      ispsma = (ispsma/5)*5 + 5
      if(odebug(40)) then
       write(iwr,*)' entered davips'
       write(iwr,*) ' ispsma ',ispsma
      endif
c
c        i28=index(28)
c        call ibasgn(lcall+lcall+3,i28,ispsma,index(28))
c
c...   gather diagonal and establish nwips
c
        call gczips_diag(index(199),index(199),.true.,iwr,cr(ipoint))
        ispsma=(nwips+510)/511
         i28=index(28)
         call ibasgn(lcall+lcall+3,i28,ispsma,index(28))
c
c...     core for davidi
c
      nbuf = (ngot-ntotw)/3
      if (nbuf.ge.nwips) then
         nbuf = nwips
         inipsc = .true.
      else
         nbuf = (nbuf/511)*511
         inipsc = .false.
      end if
c
      call davidi(cr,icr,cr(ntotw+1),cr(ntotw+nbuf+1),
     #            cr(ntotw+2*nbuf+1),nbuf,creal,zreal,inipsc,iwr)
c
      return
      end
      function dotips(ibl1,ibl2,iwr)
c
c...   calculate dot-product between 2 compacted distributed vectors
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
       parameter (nbb=5)
       parameter (nbuf=511*nbb)
      common /scra/ buf(nbuf),work(nbuf)
INCLUDE(common/disktl)
      common/nodinf/mpid,minode,mihost,ldim,nnodes,nodscr,maxnod(16)
      common/noddav/nwips,iblock1,iblock2,nwords,nread
c
      if(odebug(40)) write(iwr,*)' entered dotips'
      dd = 0.0d0
      nr = 0
      ib1 = ibl1
      ib2 = ibl2
10    nn = min(nbuf,nwips-nr)
      call rdedx(buf,nn,ib1,numcip)
      call rdedx(work,nn,ib2,numcip)
      dd = dd + ddot(nn,buf,1,work,1)
      nr = nr + nn
      if (nr.eq.nwips) then
_IF(parallel)
         call pg_dgop(10000,dd,1,'+')
_ENDIF
      else
         ib1 = ib1 + nbb
         ib2 = ib2 + nbb
         go to 10
      end if
c
      dotips = dd
      if(odebug(40)) write(iwr,*)' quit dotips'
c
      return
      end
      subroutine fgetmm(ar,max,iunit,iconf,ioff)
c...    read exactly max words and keep the rest
c...    (differs from getsb by not assuming full records)
c...    keeps track of block !!
c...    a call with max<0 initialises to block -max
      implicit REAL  (a-h,o-z)
      REAL iconf    
      dimension iconf(5,*)
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      common/bufc/ buf(511),nrest,krest,iblock
INCLUDE(common/iofile)
      dimension ar(*)
_IF(linux)
      external fget
_ENDIF
c
      if(odebug(40)) write(iwr,*)' entered fgetmm'
      if (max.lt.0) then
         nrest = 0
         iblock = -max
         if(odebug(40)) write(iwr,*)' quit fgetmm'
         return
      end if
c
      if(iposun(iunit).ne.iblock)call search(iblock,iunit)
c
      if (nrest.ne.0) then
         nw = min(max,nrest)
         call vmov(buf(krest),1,ar,1,nw)
         krest = krest + nw
         nrest = nrest - nw
      else
         nw = 0
      end if
      go to 20
c getnn(buf,nn,iunit,iconf,ioff)     
c 10    call fget(ar(nw+1),nn,iunit)
10    call getnn(ar(nw+1),nn,iunit,iconf,ioff)
      nw = nw + nn
20    if (nw+511.le.max) go to 10
c
      if (nw.ne.max) then
c 30       call fget(buf,nn,iunit)
30       call getnn(buf,nn,iunit,iconf,ioff)
         nm = min(max-nw,nn)
         call vmov(buf,1,ar(nw+1),1,nm)
         nw = nw + nm
         if (nw.ne.max) go to 30
         nrest = nn-nm
         krest = nm+1
      end if
c
      iblock = iposun(iunit)
      if(odebug(40)) write(iwr,*)' quit fgetmm'
c
      return
      end
      subroutine fgetmm_diag(ar,max,iunit)
c...    read exactly max words and keep the rest
c...    (differs from getsb by not assuming full records)
c...    keeps track of block !!
c...    a call with max<0 initialises to block -max
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      common/bufc/ buf(511),nrest,krest,iblock
INCLUDE(common/iofile)
      dimension ar(*)
_IF(linux)
      external fget
_ENDIF
c
      if(odebug(40)) write(iwr,*)' entered fgetmm'
      if (max.lt.0) then
         nrest = 0
         iblock = -max
         if(odebug(40)) write(iwr,*)' quit fgetmm'
         return
      end if
c
      if(iposun(iunit).ne.iblock)call search(iblock,iunit)
c
      if (nrest.ne.0) then
         nw = min(max,nrest)
         call vmov(buf(krest),1,ar,1,nw)
         krest = krest + nw
         nrest = nrest - nw
      else
         nw = 0
      end if
      go to 20
c
10    call fget(ar(nw+1),nn,iunit)
      nw = nw + nn
20    if (nw+511.le.max) go to 10
c
      if (nw.ne.max) then
30       call fget(buf,nn,iunit)
         nm = min(max-nw,nn)
         call vmov(buf,1,ar(nw+1),1,nm)
         nw = nw + nm
         if (nw.ne.max) go to 30
         nrest = nn-nm
         krest = nm+1
      end if
c
      iblock = iposun(iunit)
      if(odebug(40)) write(iwr,*)' quit fgetmm'
c
      return
      end
c
      subroutine copy_bloks(ibl1,ibl2,nbloks)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/disktl)
      parameter (nbuf=511*5)
      common /scra/ buf(nbuf),work(nbuf)
_IF(linux)
      external fget
_ENDIF
      ntot=0
      ifrom=ibl1
      ito=ibl2
      do i=1,nbloks
         call search(ifrom,numcip)
         call fget(buf,mmm,numcip)
         ntot=ntot+mmm
         ifrom=iposun(numcip)
         call search(ito,numcip)
         call wrt3(buf,mmm,ito,numcip)
         ito=iposun(numcip)
      enddo
      return
      end

      subroutine distrcz(a,ld,il,ih,jl,jh,ibl,myminb,mymaxb,
     1 iblokindex,omode,ibl2,oc,icblok)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/disktl)
INCLUDE(common/ga_file)
INCLUDE(common/ccntl)
INCLUDE(common/symci)
INCLUDE(common/cntrl)
      common/noddav/nwips,iblock1,iblock2,nwords,nread
      dimension a(il:il+ld-1,jl:jh)
      logical omode,oc
      dimension iblokindex(*)
c
c ...
c
      if (omode) then
         idj = 1
      else
         idj = 2
      endif
      iworkblok=1
      ibl3=ibl2
c
      icval=1
      ia=0
      ibuf=0
      itest=0
      ioff = 0
      iend=(myminb-1)*ga_file_block
      do iworkblok=myminb,mymaxb
         nnq=iblokindex(iworkblok)
         ibuf=ibuf+nnq
         itest=itest+nnq
         if (nnq.lt.ga_file_block) then
            istart=(iworkblok-1)*ga_file_block + 1
            iend=istart+ga_file_block-1
            do iij=ibuf,1,-1
               a(iend,idj)=a(istart-1-ioff+iij,idj)
               iend = iend - 1
            enddo
         endif
         if (ibuf.ge.511) then
            ia=ia+1
            call wrt3(a(iend+1,idj),511,ibl3,numcip)
            ibl3 = iposun(numcip)
            ibuf=ibuf-511
            iend=iend+511
         endif
         ioff = ibuf
      enddo
  60  call wrt3(a(iend+1,idj),ibuf,ibl3,numcip)
      if (nwips .ne. itest) call caserr('distrcz error')
      ibl3 = iposun(numcip)
      icblok=(nwips+510)/511
      return
      end

      subroutine hhhips(joop,lvar,floc,f2loc,xn,c,z,g,nbuf,iwr)
c
c...    calculate ipsc h and hh elements without too much IO
c... compute h and h**2 elements
c
      implicit REAL  (a-h,o-z)
c
INCLUDE(common/prnprn)
      dimension floc(joop),f2loc(joop),xn(joop),scra(50)
      logical lvar
c
c MPP
      dimension c(nbuf),z(nbuf),g(nbuf)
c MPP
INCLUDE(common/disktl)
      common/noddav/nwips,iblock1,iblock2,nwords,nread
c MPP
      if(odebug(40)) write(iwr,*)' entered hhhips'
      call vclr(floc(1),1,joop)
      call vclr(f2loc(1),1,joop)
c MPP
      if(nwips.gt.0) then
      ntime = (nwips-1)/nbuf + 1
      nblock = nbuf/511
      nnn = nbuf
      do 100 itime = 1,ntime
         ioff = (itime-1)*nblock
         if (itime.eq.ntime) nnn = nwips - ioff*511
         call rdedx(c,nnn,index(26+2*joop)+ioff,numcip)
         call rdedx(z,nnn,index(27+2*joop)+ioff,numcip)
         floc(joop) = floc(joop) + ddot(nnn,c,1,z,1)
         if (lvar) f2loc(joop) = f2loc(joop) + ddot(nnn,z,1,z,1)
         do 10 i=1,joop-1
            call rdedx(g,nnn,index(27+2*i)+ioff,numcip)
            floc(i) = floc(i) + ddot(nnn,c,1,g,1)
            if (lvar) f2loc(i) = f2loc(i) + ddot(nnn,z,1,g,1)
10       continue
100   continue
      endif
c
c...   normalise
c
      do 200 i=1,joop
         floc(i) = floc(i)*xn(i)*xn(joop)
         f2loc(i) = f2loc(i)*xn(i)*xn(joop)
200   continue
c MPP
c MPP    combine results
c MPP
_IF(parallel)
      call pg_dgop(10000,floc,joop,'+')
      if (lvar) call pg_dgop(10001,f2loc,joop,'+')
_ENDIF
      if(odebug(40)) write(iwr,*)' quit hhhips'
c
       return
       end
       subroutine ortips(c,z,nbuf,inipsc,ioop,idump,x,dd,xns,iwr)
c
c... orthogonalize pert vector  for davidi
c
      implicit REAL  (a-h,o-z), integer (i-n)
c
      dimension c(nbuf),z(nbuf),dd(ioop),xns(ioop)
INCLUDE(common/prnprn)
      logical inipsc
c
INCLUDE(common/disktl)
      common/noddav/nwips,iblock1,iblock2,nwords,nread
c
      if(odebug(40)) write(iwr,*)' entered ortips',inipsc
      if (inipsc) then
c...         c-vectors fit in core
10      y = x
        do 20 koop=1,ioop
           call rdedx(z,nbuf,index(26+2*koop),numcip)
           x = ddot(nbuf,c,1,z,1)*xns(koop)
_IF(parallel)
           call pg_dgop(10002,x,1,'+')
_ENDIF
           call daxpy(nbuf,-x,z,1,c,1)
20      continue
        x = ddot(nbuf,c,1,c,1)
_IF(parallel)
        call pg_dgop(10003,x,1,'+')
_ENDIF
        if (x.lt.(0.05d0*y)) goto 10
        if(odebug(40))  then
c        write(iwr,*) ' orthog pert vector ',(c(i),i=1,nbuf)
         call dbgvecv('orthog pert vector',c,nbuf)
        endif
        call wrt3(c,nbuf,idump,numcip)
      else
c...         c-vectors do'nt fit in core
c...         first do all dot-products
100     y = x
        call vclr(dd,1,ioop)
        ntime = (nwips-1)/nbuf + 1
        nblock = nbuf/511
        nnn = nbuf
        do 120 itime = 1,ntime
           ioff = (itime-1)*nblock
           if (itime.eq.ntime) nnn = nwips - ioff*511
           call rdedx(c,nnn,idump+ioff,numcip)
           do 110 i=1,ioop
            call rdedx(z,nnn,index(26+2*i)+ioff,numcip)
            dd(i) = dd(i) + ddot(nnn,c,1,z,1)
110        continue
120     continue
c
_IF(parallel)
      call pg_dgop(10000,dd,ioop,'+')
_ENDIF
      do 125 i=1,ioop
125   dd(i) = dd(i)*xns(i)
c
c..    now calculated orthog. vector (saxpy's)
c
        nnn = nbuf
        x = 0.0d0
        do 140 itime = 1,ntime
           ioff = (itime-1)*nblock
           if (itime.eq.ntime) nnn = nwips - ioff*511
           call rdedx(c,nnn,idump+ioff,numcip)
           do 130 i=1,ioop
            call rdedx(z,nnn,index(26+2*i)+ioff,numcip)
            call daxpy(nnn,-dd(i),z,1,c,1)
130        continue
           x = x + ddot(nnn,c,1,c,1)
           call wrt3(c,nnn,idump+ioff,numcip)
140     continue
_IF(parallel)
        call pg_dgop(10000,x,1,'+')
_ENDIF
        if (x.lt.(0.05d0*y)) goto 100
      end if
c
      if(odebug(40)) write(iwr,*)' quit ortips'
      return
      end
      subroutine perips(c,z,nbuf,icurc,icurz,idump,
     + eshift,h11,tester,x,iwr)
c
c... generate pert vector
c
      implicit REAL  (a-h,o-z)
      dimension c(nbuf),z(nbuf)
INCLUDE(common/prnprn)
c
c MPP
INCLUDE(common/disktl)
      common/noddav/nwips,iblock1,iblock2,nwords,nread
c MPP
      if(odebug(40)) write(iwr,*)' entered perips'
      i199 = index(199)
      tester = 0.0d0
      x = 0.0d0
c
      if(nwips.gt.0) then
      ntime = (nwips-1)/nbuf + 1
      nblock = nbuf/511
      nnn = nbuf
      if(odebug(40)) then
       write(iwr,*)' ntime,nnn,nwips ',ntime,nnn,nwips
       write(iwr,*)' ncurc,ncurz,i199,h11',icurc,icurz,i199,h11
      endif
      call pg_synch(90)
      do 100 itime = 1,ntime
         ioff = (itime-1)*nblock
         if (itime.eq.ntime) nnn = nwips - ioff*511
         call rdedx(c,nnn,icurc+ioff,numcip)
         if(odebug(40)) then
c         write(iwr,*)' c-vector',(c(i),i=1,nnn)
          call dbgvecv('c-vector',c,nnn)
         endif
         call rdedx(z,nnn,icurz+ioff,numcip)
         if(odebug(40)) then
          call dbgvecv('z-vector',z,nnn)
c         write(iwr,*)' z-vector',(z(i),i=1,nnn)
         endif
         call vsmsb(c,1,h11,z,1,c,1,nnn)
         call rdedx(z,nnn,i199+ioff,numcip)
         if(odebug(40)) then
c         write(iwr,*)' diagonal-vector',(z(i),i=1,nnn)
           call dbgvecv('diag',z,nnn)
         endif
         call vsadd(z,1,eshift-h11,z,1,nnn)
         call mrgle(nnn,0.05d0,z)
         call vdiv(c,1,z,1,c,1,nnn)
         call mrgge(nnn,0.2d0,c)
         call maxmgv(c,1,tttt,loop,nnn)
         tester = dmax1(tttt,tester)
         x=ddot(nnn,c,1,c,1) + x
         call wrt3(c,nnn,idump+ioff,numcip)
         if(odebug(40)) then
          call dbgvecv('pert vector',c,nnn)
c         write(iwr,*) ' pert vector ',(c(i),i=1,nnn)
         endif
100   continue
      endif
c MPP
_IF(parallel)
c synch() has been removed from here with tools
      call pg_dgop(10000,x,1,'+')
      call pg_dgop(10001,tester,1,'absmax')
_ENDIF
c MPP
      if(odebug(40)) write(iwr,*)' quit perips'
      return
      end
      subroutine sclips(c,nbuf,scal,icurc,iwr)
c
c...  scale a vector for davidi
c
      implicit REAL  (a-h,o-z)
      dimension c(nbuf)
c
c MPP
INCLUDE(common/disktl)
      common/noddav/nwips,iblock1,iblock2,nwords,nread
INCLUDE(common/prnprn)
c MPP
      if(odebug(40)) write(iwr,*)' entered sclips'
      if(nwips.gt.0) then
      ntime = (nwips-1)/nbuf + 1
      nblock = nbuf/511
      nnn = nbuf
      do 100 itime = 1,ntime
         ioff = (itime-1)*nblock
         if (itime.eq.ntime) nnn = nwips - ioff*511
         call rdedx(c,nnn,icurc+ioff,numcip)
         call dscal(nnn,scal,c,1)
         call wrt3(c,nnn,icurc+ioff,numcip)
100   continue
      endif
c
      if(odebug(40)) write(iwr,*)' quit sclips'
      return
      end
_ENDIF(parallel)
c
c...   Jacobi Davidson  and DIIS convergence acceleration
c...   Also for the TABLE CI (MRDCI7)
c...   Huub van Dam / Joop van Lenthe
c...   Utrecht 1995
c
      subroutine prjajd(ci,iconf,iwr)
c
c...  print jacobi-davidson control info
c...  also write out if diis 3*3 is activated
c
      implicit REAL (a-h,o-z), integer(i-n)
      character*(*) ci
_IF(cray,t3e,i8)
      integer   iconf
_ELSEIF(convex,i8drct)
      integer *8   iconf
_ELSE
      REAL  iconf
_ENDIF
      dimension iconf(5,*)
c
INCLUDE(common/sizes)
INCLUDE(common/ccntl)
INCLUDE(common/cntrl)
INCLUDE(common/diactl)
c
INCLUDE(common/dir_erg)
c
INCLUDE(common/cjdavid)
INCLUDE(common/prints)
c
      logical oprind
c
      if (maxsjd.eq.0) ifjajd = 0
      if (ifjajd.le.0) return
c
      oprind = .not.oprint(28).and..not.oprint(29)
      if (ci.eq.'direc') then
c
c...   by default do not do jacdav and vmin (for now)
c
         if (ifjajd.eq.10) then
            if (moddia.eq.1) then
               ifjajd = 0
               return
            end if
         end if
c
c...   determine the real maxsjd 
c...   -1 : ref ; -2 vacuum ; -10 : vacuum unless it is > ntotal/150
c
         maxsjd = min(nval,maxsjd)
         if (maxsjd.eq.-10) then
            maxsjd = -2
            if (ntotal/nval.lt.150) maxsjd = -1
         end if
         if (maxsjd.eq.-2) then
            maxsjd = nval
         else if (maxsjd.eq.-1) then
c...  count # csf's 
            nrfcsf = 0
            do 10 i=1,nref
               call upack2(iconf(3,i),jj,ii)
               nrfcsf = nrfcsf + ii
10          continue
            maxsjd = nrfcsf
         else if (maxsjd.ne.nval) then
c...        maxsjd should round to include whole configurations
            nn = 0
            do 20 i=1,nnst
               call upack2(iconf(3,i),jj,ii)
               nn = nn + ii
               if (nn.ge.maxsjd) go to 30
20          continue
30          maxsjd = nn
         end if
c         
         if (eshijd.ne.-1000) then 
            if(oprind) write(iwr,601) eshijd,threjd,maxcjd
601   format(' Jacobi-Davidson with ... ',/
     *       ' shift = ',e10.3,' threshold = ',e10.3,' and maximal ',
     *       i6,' cycles')
         else 
            if(oprind) write(iwr,602) threjd,maxcjd
602   format(' Jacobi-Davidson with ... ',/
     *       ' automatic shift, threshold = ',e10.3,' and maximal ',
     *       i6,' cycles')
         end if
         if(oprind) write(iwr,603) maxsjd
603   format(' considering',i7,' (vacuum) states '/)
      else if (ci.eq.'table') then
          if (eshijd.ne.-1000) then 
            if(oprind) write(iwr,611) eshijd,threjd,maxcjd,maxsjd,crijd
611   format(' Jacobi-Davidson with ... ',/
     *       ' shift = ',e10.3,' threshold = ',e10.3,' and maximal ',
     *       i6,' cycles',/,
     *       ' selecting maximal ',i6,' states with criterion ',e10.3/)
         else 
            if(oprind) write(iwr,612) threjd,maxcjd,maxsjd,crijd
612   format(' Jacobi-Davidson with ... ',/
     *       ' automatic shift, threshold = ',e10.3,' and maximal ',
     *       i6,' cycles',/,
     *       ' selecting maximal ',i6,' states with criterion ',e10.3/)
         end if
      else
         call caserr(' prjajd called wrongly ')
      end if
      if (lcall.eq.3) then
         if(oprind) write(iwr,616)
616      format(' * DIIS 3*3 acceleration of 2*2 is activated ',
     1          '(this requires 8 vectors on disk) *',/)
      end if
c
      return
      end
      subroutine gdvdmr(indcc,i199,i26,i27,
     +                  ibase,nvac,maxpa,
     +                  c,d,energy,cr,iconf,num8,ncyc)
      implicit REAL (a-h,o-z), integer(i-n)
_IF(cray,t3e,i8)
      integer iconf(5,*)
_ELSEIF(convex,i8drct)
      integer *8 iconf(5,*)
_ELSE
      REAL iconf(5,*)
_ENDIF
c
      dimension    c(nvac), d(nvac)
      dimension    cr(*)
c
INCLUDE(common/corctl)
INCLUDE(common/modeci)
INCLUDE(common/cjdavid)
c
c
c     Subroutine Generalized Jacobi Davidson for MRDCI (GJD)
c     ======================================================
c        VERSION for DIRECT   (March 1995)
c     Handle vacuum space only /  called from all davidsons 
c     ======================================================
c
c     Description:
c     ------------
c
c       This subroutine calculates the GJD update vector adapted
c       for Ci type problems. 
c
c       The GJD method is described in:
c
c       [1] "A Generalized Jacobi-Davidson iteration method for
c            linear eigenvalue problems"
c
c           Gerard L.G. Sleijpen, Henk A. van der Vorst
c           Department of Mathematics
c           Utrecht University
c
c           Preprint nr. 856
c           June 1994.
c
c       The method was adapted to Ci type problems in:
c
c       [2] "An improvement of Davidson's iteration method
c            Applications to MRCI and MRCEPA calculations"
c
c           H.J.J. van Dam, J.H. van Lenthe
c           Theoretical Chemistry Group
c           Debye Institute
c
c           Gerard L.G. Sleijpen, Henk A. van der Vorst
c           Department of Mathematics
c           Utrecht University
c
c           Submitted to: Journal of Computational Chemistry.
c           November 1994.
c
c     Input:
c     ------
c
c     <indcc>  : The dumpfile index for the coupling coefficients.
c     <i199>   : The dumpfile index for the diagonal of the Ci-matrix.
c     <i26>    : The dumpfile index for the current C-vector.
c     <i27>    : The dumpfile index for the current Z-vector.
c     <iunit>  : The dumpfile unit number.
c     <nvec>   : The size of the Ci-vector.
c     <nvac>   : The size of the vacuum space if appropriate.
c     <maxpa>  : The maximum number of spin coupling coefficients
c                per state. This defines the size of the workspace
c                in the subroutine vrijkl.
c     <c>      : The current approximation of the update-vector.
c     <d>      : Space (1 vector)
c     <energy> : The current approximation of the eigenvalue.
c
c     NOTE: <r> and <d> contain data other than their names suggest.
c           Reason for this is the need to be compatible with the
c           core partitioning in the MRDCI program. The data will be
c           moved to more appropriate locations in this subroutine.
c
c     Output:
c     -------
c
c     <ncyc>   : The number of iterations needed to solve linear
c                system.
c     <c>      : The new perturbation vector.
c
c     Workspace:
c     ----------
c
c     <cr>     : The blank common.
c     <iconf>  : The blank common.
c
      parameter (mspace=3)
c
      ncyc = 0
      if (ifjajd.eq.0) return
c
c
c...  build local core partitioning
c...  NOTE: <ibase> is the last word reserved in CI
c           <itop> is the last occupied word in this subroutine.
c
      itop = ibase
      ku   = itop + 1
      itop = itop + nvac
      kz   = itop + 1
      itop = itop + nvac
      kr   = itop + 1
      itop = itop + maxsjd
      ksz  = itop + 1
      itop = itop + maxsjd
      kdpr = itop + 1
      itop = itop + maxsjd
      ksp  = itop + 1
      itop = itop + maxsjd*mspace
      ksq  = itop + 1
      itop = itop + maxsjd*mspace
      kw   = itop + 1
      itop = itop + maxpa*maxpa
c
      if (itop.gt.ngot) then
          nuse = itop
          if (.not.odiag_1) then
             odiag_1 = .true.
             return
          end if
          call corfai('jacobi-')
      end if
c
c...  initialise datastructures
c
      call rdedx(cr(ku),nvac,i26,num8)
      call rdedx(cr(kz),nvac,i27,num8)
      call rdedx(d,nvac,i199,num8)
c
_IF(cray,ibm,vax)
      call gtrian(maxsjd,energy,cr(kr),cr(kz),cr(ku))
_ELSE
      call vsmsb(cr(ku),1,energy,cr(kz),1,cr(kr),1,maxsjd)
_ENDIF
c
c...  Currently the meaning of the datastructures is
c...  <c>  : the current update vector.
c...  <d>  : the diagonal of the Ci-matrix.
c...  <cr(kr)> : the current davidson residu-vector
c...  <cr(ku)> : the current eigenvector.
c...  <cr(kz)> : the Ci-matrix times current eigenvector.
c
c
c...  compute update vector
c
      call dirdav(indcc,maxsjd,maxpa,ncyc,maxcjd,
     +            cr(kr),d,cr(ku),cr(kz),cr(kdpr),
     +            energy,threjd,eshijd,iprjd,c,cr(ksz),cr(ksp),
     +            cr(ksq),cr(kw),iconf)
c
      if (iprjd.gt.0) write(6,601) ncyc
 601   format('+',t90,'JD-CG',i5)
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine gdvdtb(ndim,nvec,nfb,nfab,r,d,u,c,
     +                  energy,ncyc,tester,core,lword)
      implicit REAL (a-h,o-z), integer(i-n)
c
      dimension    r(nvec), d(nvec), u(nvec), c(ndim)
      dimension    core(lword)
c
INCLUDE(common/cjdavid)
c
c
c     Subroutine Generalized Jacobi Davidson for Table-CI (GJD)
c     =========================================================
c      called from drmdd
c     =========================================================
c
c     Description:
c     ------------
c
c       This subroutine calculates the GJD update vector adapted
c       for Ci type problems. 
c
c       The GJD method is described in:
c
c       [1] "A Generalized Jacobi-Davidson iteration method for
c            linear eigenvalue problems"
c
c           Gerard L.G. Sleijpen, Henk A. van der Vorst
c           Department of Mathematics
c           Utrecht University
c
c           Preprint nr. 856
c           June 1994.
c
c       The method was adapted to Ci type problems in:
c
c       [2] "An improvement of Davidson's iteration method
c            Applications to MRCI and MRCEPA calculations"
c
c           H.J.J. van Dam, J.H. van Lenthe
c           Theoretical Chemistry Group
c           Debye Institute
c
c           Gerard L.G. Sleijpen, Henk A. van der Vorst
c           Department of Mathematics
c           Utrecht University
c
c           Submitted to: Journal of Computational Chemistry.
c           November 1994.
c
c     Input:
c     ------
c
c     <ndim>   : The dimension of the Davidson subspace.
c     <nkojd>  : The number of reference configurations (in common).
c     <nvec>   : The size of the Ci-vector.
c     <nfb>    : The unit number for the subspace vectors.
c     <nfab>   : The unit number for the Ci-matrix time subspace
c                vectors.
c     <nfmx>   : The unit number for the Hamilton matrix.
c     <d>      : The diagonal of the Ci-matrix.
c     <energy> : The current approximation of the eigenvalue.
c     <lword>  : The size of the free core space.
c
c     Output:
c     -------
c
c     <r>      : The new perturbation vector.
c     <d>      : The diagonal of the Ci-matrix.
c     <nsel>   : The number of selected elements.
c     <ncyc>   : The number of iterations needed to solve linear
c                system.
c     <tester> : The current convergence tester.
c
c     Workspace:
c     ----------
c
c     <core>   : The free core space.
c     <u>      : Storage for the current eigenvector.
c
c     Variables:
c     ----------
c
c     <itop>   : The last occupied word on the <core> space.
c
      character*80 errmsg
c
      parameter (mspace=3)
c
c...  arrange some core for the z-vector
c
      ncyc = 0
      if (ifjajd.eq.0) return
c
c...  build core pre-partitioning
c
      itop  = 0
      klist = itop + 1
      itop  = itop + lenint(maxsjd)
      kz    = itop + 1
      itop  = itop + nvec
      if (itop.gt.lword) then
          write(6,601) itop,lword
601       format(' you want to use',i8,' you have ',i8)
          call caserr(' core failure in gdvdtb ')
      end if
c
c...  initialise datastructures
c
      call vclr(core(kz),1,nvec)
      do 20 i = 1, ndim
        read(nfab) (core(itop+j),j=1,nvec)
        do 10 j = 1, nvec
          core(kz+j-1) = core(kz+j-1) + c(i)*core(itop+j)
10      continue
20    continue
c
      call vclr(u,1,nvec)
      do 40 i = 1, ndim
        read(nfb) (core(itop+j),j=1,nvec)
        do 30 j = 1, nvec
          u(j) = u(j) + c(i)*core(itop+j)
30      continue
40    continue
      call rewftn(nfb)
c
      do 50 i = 1, nvec
        r(i) = core(kz+i-1) - energy*u(i)
50    continue
_IF()
c
c     for cyber205 alternative       tester=vecsum(r,r,nvec)
c     for cyber205 alternative       tester=dsqrt(tester)
_ENDIF
      tester = dnrm2(nvec,r,1)
c
c...  build the list of selected elements
c
      call selstt(nvec,u,nkojd,maxsjd,core(klist),nsel,crijd)
c
c...  now the amount of memory actually required is known.
c...  therefore, rebuild the core partitioning.
c...  i.e.: - truncate the list of selected states.
c...        - put the selected elements of the z-vector
c...          immediately after the list.
c...        - discard the original z-vector.
c
      itop  = 0
      klist = itop + 1
      itop  = itop + lenint(nsel)
      ksz   = itop + 1
      itop  = itop + nsel
c
c...  save the selected part of the z-vector
c
      call dgthr(nsel,core(kz),core(ksz),core(klist))
c
c...  now we can continue savely.
c...  build local core partitioning.
c
      ksr  = itop + 1
      itop = itop + nsel
      ksu  = itop + 1
      itop = itop + nsel
      ksd  = itop + 1
      itop = itop + nsel
      kmat = itop + 1
      itop = itop + nsel*(nsel+1)/2
      kcc  = itop + 1
      itop = itop + nsel
      kcz  = itop + 1
      itop = itop + nsel
      kpc  = itop + 1
      itop = itop + nsel*mspace
      kpz  = itop + 1
      itop = itop + nsel*mspace
c
c...  checking core requirements
c 
      if(itop.gt.lword) then 
        write(errmsg,'(a)')'Out of memory in : GDVDTB !!!'
        call caserr(errmsg)
      endif
c
c...  some more initialisations
c...  gather the required vectors for the slected states
c
      call getstt(core(kmat),nsel,core(klist),nvec,d)
      call dgthr(nsel,u,core(ksu),core(klist))
      call dgthr(nsel,d,core(ksd),core(klist))
      call dgthr(nsel,r,core(ksr),core(klist))
c
c...  initialise guess for solution vector
c
      do 60 i = 1, nsel
        core(kcc+i-1)=core(ksr+i-1)/(core(ksd+i-1)-energy)
60    continue
c
c...  compute update vector
c
      call tabdav(nsel,maxcjd,ncyc,energy,threjd,eshijd,iprjd,
     1            core(ksu),core(ksz),core(ksd),
     1            core(ksr),core(kcc),core(kcz),core(kpc),core(kpz),
     1            core(kmat))
c
c...  finish preconditioning
c
      do 70 i = 1, nvec
70    r(i) = r(i)/(d(i)-energy)
      if (ncyc.lt.maxcjd+5) 
     +  call dsctr(nsel,core(kcc),core(klist),r)
c
      if (iprjd.gt.0) write (6,602) nsel,ncyc
602   format(' Jacobi-Davidson Preconditioner selected ',i6,' cycles',
     1       i6)
c     print*,'*** OUT GDVDTB'
      return
      end
_EXTRACT(tabdav,mips4)
c
c-----------------------------------------------------------------------
c
      subroutine tabdav(nsel,maxcjd,ncyc,energy,thresh,eshift,
     1                  iprjd,su,sv,sd,sr,sc,sz,pc,pz,hmat)
      implicit REAL (a-h,o-z), integer(i-n)
      parameter (mspace=3)
c
      dimension  su(nsel), sv(nsel), sd(nsel), sr(nsel)
      dimension  sc(nsel), sz(nsel), pc(nsel,mspace), pz(nsel,mspace)
      dimension  hmat((nsel*(nsel+1))/2)
c
c
c      **   called from gdvdtb  **
c
c     Input:
c     ------
c
c     <nsel>   : The number of selected elements.
c     <maxcjd> : The maximum number of cycles for the the linear
c                system solver.
c     <energy> : The current approximation of the eigenvalue.
c     <thresh> : The Davidson convergence threshold.
c     <eshift> : The energy shift for convergence improvement.
c     <hmat>   : The selected part of the Hamilton matrix.
c     <su>     : The selected part of the current approximation of the
c                eigenvector.
c     <sv>     : The selected part of the Ci-matrix times the current
c                approximation of the eigenvector.
c
c     Output:
c     -------
c
c     <ncyc>   : The number of iterations performed.
c     <r>      : The generalized Davidson update vector.
c
c     Workspace:
c     ----------
c
c     <sd>     : Storage for selected part of Ci-matrix diagonal.
c     <sr>     : Storage for selected part of current residual vector.
c     <sc>     : Storage for selected part of the new update vector.
c     <sz>     : Storage for selected part of <hmat> times <sc>.
c     <pc>     : Storage for the subspace vectors <c>.
c     <pz>     : Storage for <hmat> times the subspace vectors <c>.
c                |z> = H|c>
c
c     Comment on Algorithm:
c     ---------------------
c
      ncyc = 0
      call dinicg
      eshifts = eshift
      if(eshifts.le.-1000.0d0) eshift = 0.0d0
      dnormm = 1000.0d0
c
c...  calculate eu = energy of selected vector
c...  use diagonal of (1-uuT)H(1-uuT) as preconditioner
c
      eu = ddot(nsel,su,1,sv,1)
      do 10 i = 1, nsel
        sd(i) = sd(i) - 2.0d0*su(i)*sv(i) + su(i)*su(i)*eu
10    continue
c
c***  Begin iterate
c
30    ncyc = ncyc + 1
c
c...    Compute sz = (B-eI)sc (see [2] equation (14))
c
        call vclr(sz,1,nsel)
        call msmatv(nsel,sz,hmat,sc)
        dotu = ddot(nsel,su,1,sc,1)
        dotv = ddot(nsel,sv,1,sc,1)
        do 40 i = 1, nsel
          sz(i) = sz(i) + su(i)*(dotu*eu-dotv) 
     1                  - sv(i)*dotu - sc(i)*energy
40      continue
c
c...    Compute the residue vector
c
        do 50 i = 1, nsel
          sz(i) = sr(i) - sz(i)
50      continue
c
c...    Compute residue norm and check convergence
c
        dnormn = dnrm2(nsel,sz,1)
c       write(*,*)'           GENDAV', ncyc, dnorm
        if(dnormn.lt.thresh.or.ncyc.ge.maxcjd) then
           goto 999
        else if (dnormn.gt.1.0d0.and.ncyc.gt.10) then
c...       signal no convergence very clearly 
           ncyc = maxcjd + 11
           goto 999
        endif
c
c...    Calling DIIS to prop up residue and solution vector.
c
        call diiscg(sc,sz,pc,pz,nsel)
        dnorm = dnrm2(nsel,sz,1)
c
c...    Dynamic shifting
c
        dnormm = dmin1(dnormm,dnorm)
        if (eshifts.le.-1000.d0.and.dnorm.gt.dnormm.and. 
     +    dnormm.gt.0.0d0) then
          eshift = eshift + 0.25d0*(dnorm/dnormm - 1.0d0)
        endif
c
        if (iprjd.gt.10) then
           write(6,601) ncyc,dnormn,dnorm,eshift
601        format(' CG ',i6,' norm ',d10.3,' diis-norm ',d10.3,
     1            ' shift ',f10.2)
        end if
c
c...    precondition 
c
        do 60 i = 1, nsel
          sz(i) = sz(i) / (sd(i) - energy + eshift)
60      continue
c
c...    update current solution vector
c
        do 70 i = 1, nsel
          sc(i) = sc(i) + sz(i)
70      continue
c
      goto 30
c
c***  End Iterate
c
999   eshift = eshifts
c
      return
      end
_ENDEXTRACT
c
c-----------------------------------------------------------------------
c
      subroutine dirdav(index,nvac,maxpa,ncyc,maxcjd,
     +                  r,d,u,z,dpr,e,
     +                  thresh,eshift,iprjd,
     +                  sc,sz,pc,pz,w,iconf)
      implicit REAL (a-h,o-z), integer(i-n)
_IF(cray,t3e,i8)
      integer iconf(5,*)
_ELSEIF(convex,i8drct)
      integer *8 iconf(5,*)
_ELSE
      REAL iconf(5,*)
_ENDIF
      parameter (mspace = 3)
c
      dimension  r(nvac), d(nvac), u(nvac), z(nvac)
      dimension  sc(nvac), sz(nvac), pc(nvac,mspace), pz(nvac,mspace)
      dimension  w(maxpa*maxpa), dpr(nvac)
c
c
c     Subroutine Direct-Ci Generalized Jacobi Davidson (GJD)
c     ======================================================
c     Handle the vacuum space only
c      ** called from gdvdmr **
c
c     Input:
c     ------
c
c     <index>  : A dumpfile index.
c     <nvac>   : The size of the vacuum space.
c     <maxpa>  : The maximum number of spin coupling coefficients
c                per vacuum state.
c     <r>      : The current Davidson residue vector.
c     <d>      : The diagonal of the Ci-matrix.
c     <u>      : The current approximation of the eigenvector.
c     <z>      : The current matrix-vector product of the approximation
c                of the eigenvector and the Ci-matrix.
c     <sc>     : The Davidson update-vector for the vacuum part.
c     <e>      : The current approximation of the eigenvalue.
c     <eshift> : The energy shift.
c
c     Output:
c     -------
c
c     <sc>     : The new update vector.
c     <ncyc>   : The number of iterations needed to solve linear
c                system.
c
c     Workspace:
c     ----------
c
c     <sc>     : Storage for vacuum part of the new update vector.
c     <sz>     : Storage mostly for matrix times <sc>.
c     <dpr>    : Storage for a backup of the original Davidson 
c                update-vector.
c     <pc>     : Storage for the <c>-vector in DIIS.
c     <pz>     : Storage for the residue vector in DIIS.
c     <w>      : Storage for the subroutine vrijkl.
c     <iconf>  : Storage for the Ci-matrix table.
c
      ncyc = 0
      call dinicg
      eshifts = eshift
      if(eshifts.le.-1000.0d0) eshift = 0.0d0
      dnormm = 1000.0d0
c
c...  save original update-vector
c
      do 10 i = 1, nvac
        dpr(i) = sc(i)
10    continue
c
c...  calculate z-vacuum (in z) and eu (= <zv|u>)
c
      call vclr(z,1,nvac)
      call vrijkl(index,nvac,z,u,w,iconf)
      eu = ddot(nvac,z,1,u,1)
c
c...  use diagonal of (1-uuT)H(1-uuT) as preconditioner
c
      do 20 i = 1, nvac
        d(i) = d(i) - 2.0d0*u(i)*z(i) + u(i)*u(i)*eu
20    continue
c
c***  Begin Iterate
c 
30    ncyc = ncyc + 1
c
c...    Compute sz = (B-eI)c (see [2] Equation (14))
c
        call vclr(sz,1,nvac)
        call vrijkl(index,nvac,sz,sc,w,iconf)
c
        dotu = ddot(nvac,u,1,sc,1)
        dotz = ddot(nvac,z,1,sc,1)
c
        do 40 i = 1, nvac
          sz(i) = sz(i) + u(i)*(dotu*eu-dotz) - z(i)*dotu - sc(i)*e
40      continue
c
c...    Compute the residue vector
c
        do 50 i = 1, nvac
          sz(i) = r(i) - sz(i)
50      continue
c
c...    Compute residue norm and check convergence
c
c       cyber205 alternative       dnorm = absmax(nvac,0.0d0,r)
c       cyber205 alternative       call maxmgv(r,1,dnorm,loop,nvac)
        loop=idamax(nvac,sz,1)
        dnorm=dabs(sz(loop))
        if(dnorm.lt.thresh.or.ncyc.gt.maxcjd) then
          goto 100
        else if (dnorm.gt.1.0d0.and.ncyc.gt.10) then
c...      signal no convergence and restore orignal update-vector
          ncyc = maxcjd
          do 60 i = 1, nvac
            sc(i) = dpr(i)
60        continue
          goto 100
       endif
c
c...   Calling DIIS to prop up residue and solution vector.
c
       dnormn = dnorm
       call diiscg(sc,sz,pc,pz,nvac)
       loop  = idamax(nvac,sz,1)
       dnorm = dabs(sz(loop))
c
c...   Dynamic shifting
c
       dnormm = dmin1(dnormm,dnorm)
       if(eshifts.le.-1000.d0.and.dnorm.gt.dnormm.and.
     +    dnormm.gt.0.0d0) then
         eshift = eshift + 0.25d0*(dnorm/dnormm - 1.0d0)
       endif
c
        if (iprjd.gt.10) then
           write(6,601) ncyc,dnormn,dnorm,eshift
601        format(' CG ',i6,' norm ',d10.3,' diis-norm ',d10.3,
     1            ' shift ',f10.2)
        end if
c
c...   Precondition residue vector
c
       do 70 i = 1, nvac
         sz(i) = sz(i) / (d(i)-e+eshift)
70     continue
c
c...   Update current solution vector
c
       do 80 i = 1, nvac
         sc(i) = sc(i) + sz(i)
80     continue
c
      goto 30
c
c***  End Iterate
c
100   continue
      eshift = eshifts
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine selstt(nvec,v,nref,maxsel,ilist,nsel,crit)
      implicit REAL (a-h,o-z), integer(i-n)
c
      dimension ilist(maxsel)
      dimension v(nvec)
c
c
c     Subroutine Select States
c     ========================
c
c     Description:
c     ------------
c
c     The subroutine selects the states from a Ci-vector that have 
c     a coefficient larger than <crit> in absolute value.
c
c     At most <maxsel> states are selected.
c
c      ** called from davtab ** used in TABLE-CI **
c
c     Input:
c     ------
c
c     <nvec>   : The size of the vector Ci-vector.
c     <v>      : The Ci-vector.
c     <nref>   : The number of reference configurations.
c     <maxsel> : The maximum number of coefficients that may
c                be select.
c     <crit>   : The selection criterion.
c
c     Output:
c     -------
c
c     <ilist>  : The list of seleted coefficients.
c     <nsel>   : The number of selected coefficients.
c
c
c...  initialise the number of selected coefficients
c
c     print*,'*** IN  SELSTT'
c     print*,'*** nref = ',nref
      nsel = 0
c
c...  Select coefficients
c
      do 10 i = 1, min(maxsel,nref)
        ilist(i) = i
10    continue
      nsel = min(maxsel,nref)
      do 20 i = min(maxsel,nref)+1, nvec
        if (dabs(v(i)).gt.dabs(crit)) then
          if (nsel.ge.maxsel) then
            return
          else
            nsel = nsel + 1
            ilist(nsel) = i
          endif
        endif
20    continue
c
c     print*,'*** OUT SELSTT'
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine getstt(hmat,nsel,ilist,nvec,dia)
      implicit REAL (a-h,o-z), integer(i-n)
      dimension ilist(nsel)
      dimension hmat(nsel*(nsel+1)/2), dia(nvec)
c
c     Subroutine Get States
c     =====================
c
c     Description:
c     ------------
c
c     Collects the selected part of the Ci matrix.
c
c      ** called from davtab ** used in TABLE-CI **
c
c     Algorithm:
c     ----------
c
c     The global outline of the algorithm is:
c
c     1.       Read a part of the matrix.
c     2.       If there are selected rows in the current part then
c     2.1        Count the number of elements in each row.
c     2.2        Start checking rows.
c     2.3        If this row is a selected row then
c     2.3.1        Start checking colums.
c     2.3.2        If this column is a selected column then
c     2.3.2.1        Store the matrix element.
c                  else
c     2.3.2.2        Skip to next column.
c     2.3.2.3        Go to 2.3.2
c                  endif
c                else
c     2.3.3        Skip to the next row.
c     2.3.4        Go to 2.3
c                endif
c              else
c     2.4        Skip current part.
c     2.5        Go to 1
c              endif
c
c     Remarks on the actual implementation:
c
c       It is assumed that the elements in <ilist> are sorted from 
c       small to large.
c
c       It is assumed that the indices (i,j) of the matrix elements
c       <hmx> are sorted in increasing order where <j> runs fastest.
c
c       It is assumed that <j> is less than <i> for all indices (i,j).
c
c     Input:
c     ------
c
c     <nsel>  : the number of selected states.
c     <ilist> : the indices of the selected states.
c     <nvec>  : the dimension of the Ci vector
c     <dia>   : the diagonal of the Ci matrix
c
c     Output:
c     -------
c
c     <hmat>  : the selected part of the interaction matrix.
c
c     Variables:
c     ----------
c
c     COMMON /FTAPE/ : File-IO data, mostly unit numbers.
c
      common /ftape/ nfb,nfab,nf11,nfrsut,nfsc,nfmx,
     +               nspacf(6),nfmain,nfci,idum4(6)
c
c     Workspace to store the Ci-matrix elements.
c     (See also the GAMESS subroutine DRMDD)
c
      parameter(ndmax = 1600)
      dimension  hmx(ndmax), ind(ndmax), jnd(ndmax)
c
c...  initialise matrix
c
      call vclr(hmat, 1, nsel*(nsel+1)/2)
c
      if(nsel.le.0) return
c
c...  Read part of the matrix.
c
   10 kcurr = 1
      isel  = 1
        read (nfmx,end=50,err=50) nd,hmx,ind,jnd
        if (nd.eq.-999) go to 50
   20   if(ind(1).gt.ilist(isel)) then
          if(isel.ge.nsel) then
c
c...        There is no selected row in this block.
c...        Skip to next block.
c
            go to 10
          else
c
c...        The next selected row may be in this block.
c
            isel = isel + 1
            go to 20
          endif
        else if(ind(nd).lt.ilist(isel)) then
c
c...      There is no selected row in this block.
c...      Skip to next block.
c
          go to 10
        else
c
c...      There may be a selected row in this block.
c
          nrow = ind(nd)-ind(1)+1
          if(nrow.eq.1) then
            jsel = 1
            do 100 kcurr = 1, nd
 22           if(jnd(kcurr).gt.ilist(jsel)) then
c
c...            The next selected columns may be in this block.
c
                jsel = jsel + 1
                go to 22
              else if(jnd(kcurr).eq.ilist(jsel)) then
c
c...            This is a selected element.
c...            Store it.
c
                hmat(isel*(isel-1)/2+jsel) = hmx(kcurr)
                jsel = jsel + 1
              else
c
c...            This is not an selected element.
c...            Skip to next element.
c
              endif
100        continue
          else
            indold = -1
            do 110 kcurr = 1, nd
 25           if(ind(kcurr).gt.ilist(isel)) then
                if(isel.ge.nsel) then
c
c...              There is no selected row left in this block.
c...              Skip to the next block.
c
                  go to 10
                else
c
c...              The next selected row may be in this block.
c
                  isel = isel + 1
                  go to 25
                endif
              else if(ind(kcurr).eq.ilist(isel)) then
                if(ind(kcurr).ne.indold) then
                  indold = ind(kcurr)
                  jsel   = 1
                endif
 27             if(jnd(kcurr).gt.ilist(jsel)) then
c
c...              The next selected columns may be in this block.
c
                  jsel = jsel + 1
                  go to 27
                else if(jnd(kcurr).eq.ilist(jsel)) then
c
c...              This is a selected element.
c...              Store it.
c
                  hmat(isel*(isel-1)/2+jsel) = hmx(kcurr)
                  jsel = jsel + 1
                else
c
c...              This is not an selected element.
c...              Skip to next element.
c
                endif
              else
c
c...            We don't want this element
c     
              endif
110         continue
            go to 10
          endif
        endif
c
   50 call rewftn(nfmx)
c
      do 120 i = 1, nsel
        hmat(i*(i+1)/2) = dia(ilist(i))
120   continue
c  
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine msmatv(nvec,vres,smat,vin)
      implicit REAL (a-h,o-z), integer(i-n)
      dimension vres(nvec), smat(nvec*(nvec+1)/2), vin(nvec)
c
c     Subroutine Multiply Symmetric Matrix with Vector
c     =======================================================
c
c     Description:
c     ------------
c
c     A symmetric matrix is multiplied by a vector and added to
c     the result vector.
c
c     It is assumed that the matrix is stored as a lower triangular
c     matrix with diagonal.
c
c     Input:
c     ------
c
c     <nvec>   : The size of the vectors.
c     <vres>   : The result vector.
c     <smat>   : The symmetric matrix.
c     <vin>    : The input vector.
c
c     Output:
c     -------
c
c     <vres>   : The original result vector plus the matrix-vector
c                product.
c
      k = 1
      do 20 i = 1, nvec
        do 10 j = 1, i-1
          vres(i) = vres(i) + smat(k)*vin(j)
          vres(j) = vres(j) + smat(k)*vin(i)
10      k = k + 1
        vres(i) = vres(i) + smat(k)*vin(i)
20    k = k + 1
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine vrijkl(index,nvac,z,c,w,conf)
c
      implicit none
c
      integer index,nvac
_IF(cray,t3e,i8)
      integer conf(5,*)
_ELSEIF(convex,i8drct)
      integer *8 conf(5,*)
_ELSE
      REAL conf(5,*)
_ENDIF
      REAL z(nvac), c(nvac), w(*)
c
c     Subroutine Multiply Vacuum Matrix with Vector
c     =============================================
c     Stolen from DIRECT CI (vvijkl)
c
c     Description:
c     ------------
c
c     Adds the matrix-vector product of a vector and the vacuum
c     Ci-matrix to the result vector: 
c
c       z' = z + A.c
c
c     where A is the vacuum part of the Ci-matrix.
c
c     Input:
c     ------
c
c     <index>  : A dumpfile index.
c     <nvac>   : The size of the vacuum space.
c     <c>      : The input vector.
c     <z>      : The result vector.
c
c     Output:
c     -------
c
c     <z>      : The original result vector plus the matrix-vector
c                product.
c
c     Workspace:
c     ----------
c
c     <w>      : Storage for elements of the vacuum Ci-matrix.
c     <conf>  : The blank common.
c
c
c...  Equivalencing to get addresses right.(see also subroutine vvijkl).
c
      REAL rmodel
      common/spew/ rmodel(3)
      equivalence (rmodel(1),model),(rmodel(2),ic),(rmodel(3),ir)
      integer model,ic,ir
      integer nc,nr,icbas,irbas
_IFN(cray,convex,t3e,i8,i8drct)
      integer ipad
_ENDIF
c
      call brtsb(index)
      call getsb(rmodel(1),1)
      if(model.ne.2) then
_IF(parallel)
        return
_ENDIF
        write(*,*)'*** model is ',model,' !!!'
        call caserr('*** no vacuum data !!!')
      endif
999   call getsb(rmodel(2),2)
c
_IF(cray,convex,t3e,i8,i8drct)
      nc=conf(3,ic)
      nr=conf(3,ir)
_ELSE
      nc=ipad(conf(3,ic))
      nr=ipad(conf(3,ir))
_ENDIF
      call getsb(w(1),nc*nr)
_IF(cray,convex,t3e,i8,i8drct)
      icbas=conf(1,ic)
      irbas=conf(1,ir)
_ELSE
      icbas=ipad(conf(1,ic))
      irbas=ipad(conf(1,ir))
_ENDIF
c
      if (icbas+nc.gt.nvac) return
c
      call mxmb(w,1,nr,c(icbas+1),1,1,z(irbas+1),1,1,nr,nc,1)
      if(ic.ne.ir)
     1 call mxmb(w,nr,1,c(irbas+1),1,1,z(icbas+1),1,1,nc,nr,1)
c
      call getsb(rmodel(1),1)
      if(model.eq.2)goto 999
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine srthmx(nd,hmx,ind,jnd)
      implicit REAL (a-h,o-z), integer(i-n)
      dimension hmx(nd), ind(nd), jnd(nd)
c
c     SORT NFMX
c
c     This subroutine sorts the hamilton matrix file such that
c     for all entries 
c
c        h(i), ind(i), jnd(i)
c
c     and
c
c        h(i+1), ind(i+1), jnd(i+1)
c
c     holds
c
c          (ind(i+1).gt.ind(i)).or.
c        (((ind(i+1).eq.ind(i)).and.(jnd(i+1).ge.jnd(i)))
c
c     VARIABLES:
c     ----------
c
c     <ielm>  : the index of the element currently being sorted.
c     <ifree> : the index of the free position in the arrays.
c     <velmh> : the value of the current element in the h-matrix.
c     <ielmi> : the value of the current i-index.
c     <ielmj> : the value of the current j-index.
c
c
      do 20 ielm = 2, nd
        velmh = hmx(ielm)
        ielmi = ind(ielm)
        ielmj = jnd(ielm)
        ifree = ielm
10      if (ifree.le.1) then
          hmx(1) = velmh
          ind(1) = ielmi
          jnd(1) = ielmj
        else
          if (ielmi.gt.ind(ifree-1).or.
     1      (ielmi.eq.ind(ifree-1).and.ielmj.ge.jnd(ifree-1))) then
            hmx(ifree) = velmh
            ind(ifree) = ielmi
            jnd(ifree) = ielmj
          else
            hmx(ifree) = hmx(ifree-1)
            jnd(ifree) = jnd(ifree-1)
            ind(ifree) = ind(ifree-1)
            ifree      = ifree - 1
            go to 10
          endif
        endif
20    continue
c
      return
      end
      subroutine diiscg(c,z,ci,zi,ntotal)
c
      implicit REAL (a-h,o-z), integer(i-n)
c
      dimension c(ntotal),z(ntotal),ci(ntotal,3),zi(ntotal,3),index(3)
CMR   replaced entry dinicg+save var with subroutine dinicg+save common
      common /diniini/index,ncyc
      save /diniini/
c
c...  perform 3*3 diis convergence acceleration for Conjugate-Gradient
c...  for jacobi davidson ... linear equation accellerator
c...  ** copied from the original one that worked from disk for mp2/3
c...  All in core version
c
c...   c,z   :  work-space 
c...   ci,zi :  c/z storage 
c...   index :  core-addresses for vector-storage (1/2/3) : iteration i/i-1/i-2
c...   on entrance and exit : z = current residuum 
c...                          c = current solution
c...   // peter pulay, Fayetteville Arkansas 1989 \\
c...   // Joop v.L Lenthe - Huub v. Dam  Utrecht 1995 \\
c 
      dimension diism(3,3),lll(3),mmm(3)
c
      ncyc = ncyc + 1
c
c...     swap residu and c addresse
c  
         ind = index(3)
         index(3) = index(2)
         index(2) = index(1)
         index(1) = ind
c
c...     dump current c/z vectors (z-vector is in core)
c
         call dcopy(ntotal,z,1,zi(1,index(1)),1)
         call dcopy(ntotal,c,1,ci(1,index(1)),1)
c
      if (ncyc.gt.2) then
c
c...    dot product of current residuum with  c's
c    
c        call rdedx(c,ntotal,index(1,1),8)  is already in core
c
         diism(1,1) = ddot(ntotal,ci(1,index(1)),1,zi(1,index(1)),1)
         diism(2,1) = ddot(ntotal,ci(1,index(2)),1,zi(1,index(1)),1)
         diism(3,1) = ddot(ntotal,ci(1,index(3)),1,zi(1,index(1)),1)
c
c...   dot product of previous residuum with c's   
c
         diism(1,2) = ddot(ntotal,ci(1,index(1)),1,zi(1,index(2)),1)
         diism(2,2) = ddot(ntotal,ci(1,index(2)),1,zi(1,index(2)),1)
         diism(3,2) = ddot(ntotal,ci(1,index(3)),1,zi(1,index(2)),1)
c
c...   dot product of oldest residuum with c's   
c
         diism(1,3) = ddot(ntotal,ci(1,index(1)),1,zi(1,index(3)),1)
         diism(2,3) = ddot(ntotal,ci(1,index(2)),1,zi(1,index(3)),1)
         diism(3,3) = ddot(ntotal,ci(1,index(3)),1,zi(1,index(3)),1)
c
c...  Now the matrix diism(3,3) contains <c(i)|z(j)> where i,j=1
c...  corresponds to the current value, 2 to the previous one, 
c...  and 3 to the value two before. Then row 2 is subtracted from
c...  row 1, and row 3 from row 2. This is equivalent of orthogonalizing
c...  the interpolated residue to (C1-C2) and (C2-C3). The third
c...  row of the diism matrix is all 1's: the sum of the three
c...  coefficients is 1. This arrangement minimizes numerical errors
c
c...     according to Hestenes' theory, the interpolated Z must
c...     be orthogonal to delta(present-C) and delta(previous-C)
c...     a good summary of this is given by Wormer,Paldus IJQC ???
c..      see also : ????
c        
        do 635 iii=1,3
          diism(1,iii) = diism(1,iii)-diism(2,iii)
          diism(2,iii) = diism(2,iii)-diism(3,iii)
          diism(3,iii) = 1.0d0
 635    continue
c          write(*,987) diism(1,1),diism(1,2),diism(1,3),diism(2,1),
c     1      diism(2,2),diism(2,3),diism(3,1),diism(3,2),diism(3,3)
c
c...    normalize the equations
c
        s1 = diism(1,1)
        if (abs(s1).lt.1.0d-15) s1=1.0d-15
        s2 = diism(2,2)
        if (abs(s2).lt.1.0d-15) s2=1.0d-15
        do 695 iii=1,3
           diism(1,iii)=diism(1,iii)/s1
           diism(2,iii)=diism(2,iii)/s2
 695    continue
c          write(*,987) diism(1,1),diism(1,2),diism(1,3),diism(2,1),
c     1      diism(2,2),diism(2,3),diism(3,1),diism(3,2),diism(3,3)
 987      format(' diis',5x, 3e13.6,/,10x,3e13.6,/,10x,3e13.6)
c
c...    now invert the diism matrix. Its last column will give the coefficients 
c
        call osinv1(diism,3,determ,1.0d-14 ,lll,mmm)       
c        write(6,988) determ,(diism(kk,3),kk=1,3)
c 988    format(11x,'DIIS det ',e11.4,' coeff',3F13.6) 
c
c...    calculate extrapolated c as linear combination of the c's
c...    cex = c(1)*diism(1,3) + c(2)*diism(2,3) + c(3)*diism(3,3)
c...    replace c by it's extrapolated counterpart 
c
        call dscal(ntotal,diism(1,3),c,1)
        call daxpy(ntotal,diism(2,3),ci(1,index(2)),1,c,1)
        call daxpy(ntotal,diism(3,3),ci(1,index(3)),1,c,1)
c
c...    calculate extrapolated residue as linear combination of the z's
c...    zex = z(1)*diism(1,3) + z(2)*diism(2,3) + z(3)*diism(3,3)   
c...    exactly like c's / abd leave result in z (to return)
c...    replace z by it's extrapolated friend 
c
        call dscal(ntotal,diism(1,3),z,1)
        call daxpy(ntotal,diism(2,3),zi(1,index(2)),1,z,1)
        call daxpy(ntotal,diism(3,3),zi(1,index(3)),1,z,1)
c
c...   check - calculate the scalar products with the c's
c        call getcz(cr(kmpz+1),index(30))
c        s1=ddot(ntotal,cr(kmpc+1),1,cr(kmpz+1),1)
c        call getcz(cr(kmpz+1),index(31))
c        s2=ddot(ntotal,cr(kmpc+1),1,cr(kmpz+1),1)
c        write(*,*) 's1,s2', s1,s2 
c
      endif
c
      return
      end
c
c...  initialise
c
      subroutine dinicg
c
      implicit integer(i-n)
      dimension index(3)
      common /diniini/index,ncyc
      save /diniini/
      ncyc = 0
      index(1) = 1
      index(2) = 2
      index(3) = 3
c
      return
      end

      subroutine osinv1 (a,n,d,tol,l,m)
      implicit REAL (a-h,o-z)
c
c     parameters:  a - input matrix , destroyed in computation and repla
c                      by resultant inverse (must be a general matrix)
c                  n - order of matrix a
c                  d - resultant determinant
c            l and m - work vectors of lenght n
c                tol - if pivot element is less than this parameter the
c                      matrix is taken for singular (usually = 1.0e-8)
c     a determinant of zero indicates that the matrix is singular
c
      dimension a(*), m(*), l(*)
      d=1.0d0
      nk=-n
      do 180 k=1,n
         nk=nk+n
         l(k)=k
         m(k)=k
         kk=nk+k
         biga=a(kk)
         do 20 j=k,n
            iz=n*(j-1)
         do 20 i=k,n
            ij=iz+i
c
c     10 follows
c
            if (dabs(biga)-dabs(a(ij))) 10,20,20
   10       biga=a(ij)
            l(k)=i
            m(k)=j
   20    continue
         j=l(k)
         if (j-k) 50,50,30
   30    ki=k-n
         do 40 i=1,n
            ki=ki+n
            holo=-a(ki)
            ji=ki-k+j
            a(ki)=a(ji)
   40    a(ji)=holo
   50    i=m(k)
         if (i-k) 80,80,60
   60    jp=n*(i-1)
         do 70 j=1,n
            jk=nk+j
            ji=jp+j
            holo=-a(jk)
            a(jk)=a(ji)
   70    a(ji)=holo
   80    if (dabs(biga)-tol) 90,100,100
   90    d=0.0d0
         return
  100    do 120 i=1,n
            if (i-k) 110,120,110
  110       ik=nk+i
            a(ik)=a(ik)/(-biga)
  120    continue
         do 150 i=1,n
            ik=nk+i
            ij=i-n
         do 150 j=1,n
            ij=ij+n
            if (i-k) 130,150,130
  130       if (j-k) 140,150,140
  140       kj=ij-i+k
            a(ij)=a(ik)*a(kj)+a(ij)
  150    continue
         kj=k-n
         do 170 j=1,n
            kj=kj+n
            if (j-k) 160,170,160
  160       a(kj)=a(kj)/biga
  170    continue
         d=d*biga
         a(kk)=1.0d0/biga
  180 continue
      k=n
  190 k=k-1
      if (k) 260,260,200
  200 i=l(k)
      if (i-k) 230,230,210
  210 jq=n*(k-1)
      jr=n*(i-1)
      do 220 j=1,n
         jk=jq+j
         holo=a(jk)
         ji=jr+j
         a(jk)=-a(ji)
  220 a(ji)=holo
  230 j=m(k)
      if (j-k) 190,190,240
  240 ki=k-n
      do 250 i=1,n
         ki=ki+n
         holo=a(ki)
         ji=ki+j-k
         a(ki)=-a(ji)
  250 a(ji)=holo
      go to 190
  260 return
c
      end
      subroutine diisci(c,z,index,ntotal,indc,indz,energy,
c    *                  conf,ifcepa,nref,num8,iwr)
     *                  conf,nref,num8,iwr)
c
      implicit REAL (a-h,o-z), integer(i-n)
c
      dimension c(ntotal),z(ntotal),index(3,2),conf(5,*)
c
c...  perform 3*3 CONJUGATE-GRADIENT convergence acceleration for ci
C...  i.e. overlap between gradient and C is minimised
c...   c,z   :  work-space 
c...            2 vectors in core are needed and 6 on disk
c...   index :  diskaddresses on ed7 for vector-storage
c...            (1/2/3,1) : c-vectors of iteration i/i-1/i-2
c...            (1/2/3,2) : z-vectors of iteration i/i-1/i-2 
c...   in DIRECT index is probably index(31-36)
c...   indc/indz : address for c/z vector (dumpcz/getcz format)
c...   energy : current energy (for residue calc)
c...   conf/ifcepa : info needed for cepa adjustment to z-vector 
c...   nref : # reference config's (keep overlap of reference parts
c...          of c-vectors >0 (avoid phase problems) 
c...   on entrance and exit : c/z undefined 
c...   // Peter pulay, J.H van Lenthe Fayetteville Arkansas 1989 \\
c
INCLUDE(common/cjdavid)
c
      dimension diism(3,3),lll(3),mmm(3)
      save ncyc
      data ncyc/0/
c
      ncyc = ncyc + 1
c
c...     swap residu and c addresses
c  
         do 10 i=1,2
            ind = index(3,i)
            index(3,i) = index(2,i)
            index(2,i) = index(1,i)
            index(1,i) = ind
10       continue
c
c...     dump current c/z vectors 
c 
            call getcz(conf,z,indz)
            call wrt3(z,ntotal,index(1,2),num8)
            call getcz(conf,c,indc)
            call wrt3(c,ntotal,index(1,1),num8)
c
      if (ncyc.gt.2) then
c
c...   sort out phases  (for security)
c
        nr = max(nref,1)
        call rdedx(c,nr,index(1,1),num8)
        call rdedx(z,nr,index(2,1),num8)
        fase2 = ddot(nr,c,1,z,1)
        if (fase2.eq.0.0d0) call caserr('diisci fase failure2')
        fase2 = dsign(1.0d0,fase2)
        call rdedx(z,nr,index(3,1),num8)
        fase3 = ddot(nr,c,1,z,1)
        if (fase3.eq.0.0d0) call caserr('diisci fase failure3')
        fase3 = dsign(1.0d0,fase3)

c
c...    calculate current residuum
c  
        call rdedx(c,ntotal,index(1,1),num8)
        call rdedx(z,ntotal,index(1,2),num8)
c       if (ifcepa.gt.0) call cepvds(conf,z,z,c,.false.,.true.)
        call daxpy(ntotal,-energy,c,1,z,1) 
c        
c...    dot product of current residuum with  c's
c
         diism(1,1) = ddot(ntotal,c,1,z,1)
         call rdedx(c,ntotal,index(3,1),num8)
         diism(3,1) = ddot(ntotal,c,1,z,1) * fase3
         call rdedx(c,ntotal,index(2,1),num8)
         diism(2,1) = ddot(ntotal,c,1,z,1) * fase2
c
c...    calculate previous residuum (with new energy and shifts)
c 
        call rdedx(z,ntotal,index(2,2),num8)
c       if (ifcepa.gt.0) call cepvds(conf,z,z,c,.false.,.true.)
        call daxpy(ntotal,-energy,c,1,z,1) 

c...   dot product of previous residuum with c's   
c
         diism(2,2) = ddot(ntotal,c,1,z,1)
         call rdedx(c,ntotal,index(1,1),num8)
         diism(1,2) = ddot(ntotal,c,1,z,1) * fase2
         call rdedx(c,ntotal,index(3,1),num8)
         diism(3,2) = ddot(ntotal,c,1,z,1) * fase2*fase3
c
c...    calculate oldest residuum (with new energy and shifts))
c
        call rdedx(z,ntotal,index(3,2),num8)
c       if (ifcepa.gt.0) call cepvds(conf,z,z,c,.false.,.true.)
        call daxpy(ntotal,-energy,c,1,z,1) 
c
c...   dot product of oldest residuum with c's   
c
         diism(3,3) = ddot(ntotal,c,1,z,1) 
         call rdedx(c,ntotal,index(2,1),num8)
         diism(2,3) = ddot(ntotal,c,1,z,1) * fase2*fase3
         call rdedx(c,ntotal,index(1,1),num8)
         diism(1,3) = ddot(ntotal,c,1,z,1) * fase3
c          write(*,987) diism(1,1),diism(1,2),diism(1,3),diism(2,1),
c     1      diism(2,2),diism(2,3),diism(3,1),diism(3,2),diism(3,3)
c
c...  Now the matrix diism(3,3) contains <c(i)|z(j)> where i,j=1
c...  corresponds to the current value, 2 to the previous one, 
c...  and 3 to the value two before. Then row 2 is subtracted from
c...  row 1, and row 3 from row 2. This is equivalent of orthogonalizing
c...  the interpolated residue to (C1-C2) and (C2-C3). The third
c...  row of the diism matrix is all 1's: the sum of the three
c...  coefficients is 1. This arrangement minimizes numerical errors
c
c...     according to Hestenes' theory, the interpolated Z must
c...     be orthogonal to delta(present-C) and delta(previous-C)
c...     a good summary of this is given by Wormer,Paldus IJQC ???
c..      see: P.E.S. Wormer,F.Visser,J.Paldus J.Comput.Chem.48 (1982) 23
c        
        do 635 iii=1,3
          diism(1,iii) = diism(1,iii)-diism(2,iii)
          diism(2,iii) = diism(2,iii)-diism(3,iii)
          diism(3,iii) = 1.0d0
 635    continue
c          write(*,987) diism(1,1),diism(1,2),diism(1,3),diism(2,1),
c     1      diism(2,2),diism(2,3),diism(3,1),diism(3,2),diism(3,3)
c
c...    normalize the equations
c
        s1 = diism(1,1)
        if (abs(s1).lt.1.0d-15) s1=1.0d-15
        s2 = diism(2,2)
        if (abs(s2).lt.1.0d-15) s2=1.0d-15
        do 695 iii=1,3
           diism(1,iii)=diism(1,iii)/s1
           diism(2,iii)=diism(2,iii)/s2
 695    continue
c          write(*,987) diism(1,1),diism(1,2),diism(1,3),diism(2,1),
c     1      diism(2,2),diism(2,3),diism(3,1),diism(3,2),diism(3,3)
c987      format(' diis',5x, 3e13.6,/,10x,3e13.6,/,10x,3e13.6)
c
c...    now invert the diism matrix. Its last column will give the coefficients 
c
        call osinv1(diism,3,determ,1.0d-14 ,lll,mmm)       
        if (iprjd.gt.0) write(iwr,988) determ,(diism(kk,3),kk=1,3)
 988    format(11x,'DIIS det ',e11.4,' coeff',3e13.6) 
c
c...    calculate extrapolated c as linear combination of the c's
c...    cex = c(1)*diism(1,3) + c(2)*diism(2,3) + c(3)*diism(3,3)
c...    including phases
c
        call rdedx(c,ntotal,index(3,1),num8)
        call dscal(ntotal,diism(3,3)*fase3,c,1)
        call rdedx(z,ntotal,index(2,1),num8)
        call daxpy(ntotal,diism(2,3)*fase2,z,1,c,1)
        call rdedx(z,ntotal,index(1,1),num8)
        call daxpy(ntotal,diism(1,3),z,1,c,1)
c                                              
c....   normalise c (keep normalisation for z)
c
        cnorm = ddot(ntotal,c,1,c,1)
        cnorm = 1.0d0/dsqrt(cnorm)
        call dscal(ntotal,cnorm,c,1)
c
c...    replace c by it's extrapolated counterpart 
c...    also dump it (for whoever)
c
        call wrt3(c,ntotal,index(1,1),num8)
        call dumpcz(conf,c,indc) 
c
c...    calculate extrapolated z as linear combination of the z's
c...    zex = z(1)*diism(1,3) + z(2)*diism(2,3) + z(3)*diism(3,3)   
c...    exactly like c's 
c
        call rdedx(z,ntotal,index(3,2),num8)
        call dscal(ntotal,diism(3,3)*fase3,z,1)
        call rdedx(c,ntotal,index(2,2),num8)
        call daxpy(ntotal,diism(2,3)*fase2,c,1,z,1) 
        call rdedx(c,ntotal,index(1,2),num8)
        call daxpy(ntotal,diism(1,3),c,1,z,1) 
        call dscal(ntotal,cnorm,z,1)
c
c...    replace z by it's extrapolated friend 
c
        call wrt3(z,ntotal,index(1,2),num8)
        call dumpcz(conf,z,indz)
        call rdedx(c,ntotal,index(1,1),num8) 
      end if
c
c...    diis energy
c
c       if (ifcepa.gt.0) call cepvds(conf,z,z,c,.false.,.true.)
        energo = energy
        energy = ddot(ntotal,c,1,z,1)
        if (iprjd.gt.0) write(iwr,999) energo,energy
999     format(6x,'pre-diis energy ',f16.10,' after diis ',f16.10)
c
      return
      end
      subroutine ver_dirctb(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/dirctb.m,v $
     +     "/
      data revision /"$Revision: 6195 $"/
      data date /"$Date: 2010-09-17 16:11:09 +0200 (Fri, 17 Sep 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
c
_IF(ga)
      subroutine gczips_diag(iblkf,iblkt,oc,iwr,iblokind)
c works only for ga-file ...
c      add up a c/z-vector and scatter it over the nodes
c      partial vector at iblkf, compact result at iblkt
c       if (.not.oc) do not add just compact on each node (for c)
c      in nw return # words written
c      kind of getcz action
c      *note* the n-2 c/z vectors are singlet/triplet adjacent on node
c
      implicit REAL  (a-h,o-z)
      dimension iblokind(*)
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      parameter (nbuf=511*5)
      common /scra/ buf(nbuf),work(nbuf)
INCLUDE(common/disktl)
      common/nodinf/mpid,minode,mihost,ldim,nnodes,nodscr,maxnod(16)
      common/noddav/nwips,iblock1,iblock2,nwords,nread
INCLUDE(common/cntrl)
INCLUDE(common/ccntl)
INCLUDE(common/ga_file)
INCLUDE(common/symci)
      common/remc_ed19/nnblok,ipoint,idiagoff
c
       logical oipsci,oc,omode
c
      if(odebug(40)) write(iwr,*)' entered gczips'
      ibl = iblkt
c
c
      iga_to_use = ga_file_unit_to_ga( numcz )
      call pg_distribution(ga_file_handle(iga_to_use),ipg_nodeid(),
     1 il,ih,jl,jh)
      myminb = (il-1)/ga_file_block + 1
      mymaxb = ih/ga_file_block
      nwips=0
      do iii=myminb,mymaxb
         nwips=nwips+iblokind(iii)
      enddo
      call fgetmm_diag(buf,-iblkf,numcip)
      np = 0
      nw = 0
      nnproc = ipg_nnodes()
      do iik = 0, nnproc - 1
      nr = 0
      ibuf=nwips
_IF(i8,cray,convex,t3e,i8drct)
      call pg_brdcst(2211,ibuf,8,iik)
_ELSE
      call pg_brdcst(2211,ibuf,4,iik)
_ENDIF
  70    nn = min(nbuf,ibuf-nr)
      call fgetmm_diag(buf,nn,numcip)
_IF(parallel)
      if (oc) call pg_dgop(12340,buf,nn,'+')
_ENDIF
      if (ipg_nodeid().eq.iik) then
            call wrt3(buf,nn,ibl,numcip)
            ibl = iposun(numcip)
      endif
      nr = nr + nn
      nw = nw + nn
      if (nr.ne.ibuf) go to 70
      enddo
      if (nw.ne.ntotal) then
         call caserr('gczips error')
      endif
      return
      end

      subroutine getvacdiag(diag)
      implicit REAL (a-h,o-z),integer (i-n)
INCLUDE(common/sizes)
INCLUDE(common/cntrl)
INCLUDE(common/disktl)
      parameter (nbuf=511*5)
      common /scra/ buf(nbuf),work(nbuf)
      common/noddav/nwips,iblock1,iblock2,nwords,nread
      dimension diag(nval)
_IF(linux)
      external fget
_ENDIF
      logical opg_root
      nnodes=ipg_nnodes()
      call search(index(199),numcip)
      ioff=1
      nstill=nval
      do iik=0,nnodes-1
       nonthis=nwips
       ioff_old=ioff
       iread=0
       mread=min(nonthis,nstill)
       if (ipg_nodeid().eq.iik) then
 10     call fget(buf,nn,numci)
        mm=min(mread,nn)
        iread=iread+mm
        mread=mread-mm
        nstill=nstill-mm
        call dcopy(mm,buf,1,diag(ioff),1)
        ioff=ioff+mm
        if (mread) 20,20,10
 20     continue
       endif
_IF(i8,cray,convex,t3e,i8drct)
      call pg_brdcst(2211,nstill,8,iik)
      call pg_brdcst(2210,iread,8,iik)
      call pg_brdcst(2213,ioff,8,iik)
_ELSE
      call pg_brdcst(2211,nstill,4,iik)
      call pg_brdcst(2210,iread,4,iik)
      call pg_brdcst(2213,ioff,4,iik)
_ENDIF
      call pg_brdcst(2212,diag(ioff_old),iread*8,iik)
      if ((nstill.lt.0).or.(mread.lt.0)) call caserr('negative error')
      if (nstill.eq.0) goto 30
      enddo
  30  return
      end

      subroutine getnn(buf,nn,iunit,iconf,ioff)
      implicit REAL  (a-h,o-z),integer  (i-n)  
_IF(cray,t3e,i8)
      integer   iconf
_ELSEIF(convex,i8drct)
      integer *8   iconf
_ELSE
      REAL iconf
_ENDIF
      dimension iconf(5,*),buf(*)
INCLUDE(common/sizes)
INCLUDE(common/cntrl)
INCLUDE(common/ccntl)
INCLUDE(common/cdcryz)
INCLUDE(common/symci)
INCLUDE(common/disktl)
INCLUDE(common/disc)
      save icount
c  
      icval=ioff
      call find(iunit)
      iblok=ipos(isel)
      iaddr=iblok-ioff
      if (iaddr) 20,10,20
c
c  iaddr = 0 start vacuum
c
  10  nn=nval
      icount=1
      goto 30
c
c  iaddr .ne. 0 --> doublet or n-2 if icount.eq. 1 then doublet
c  if icount.eq. 2 n-2
c
_IF(cray,convex,t3e,i8,i8drct)
  20  icdoub=iconf(2,1)
_ELSE
  20  icdoub=ipad(iconf(2,1))
_ENDIF
c
      if (icount.eq.1.and.iaddr.eq.icdoub) then
         nn=ndoub
         icount=nstrt2
         goto 30
      endif
  35  loop=icount
      call upack2(iconf(2,loop),kc,jc)
      call upack2(iconf(3,loop),lamt,lams)
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,loop)
_ELSE
      isi=ipad(iconf(4,loop))
_ENDIF
      nss=lams*norbsi(isi)
      nst=lamt*norbtr(isi)
      if (iaddr.eq.jc) then
         nn=nss
         if (nn.eq.0) goto 40
         goto 30
      endif
      if (iaddr.eq.kc) then
         nn=nst
         if (nn.eq.0) goto 40
         goto 30
      endif
  40  icount=icount+1
      if (icount.gt.nlast2) call caserr('icount getnn')
      goto 35
  30  call rdedx(buf,nn,iblok,iunit)
      return
      end

      subroutine czdumpga(a,ld,il,ih,jl,jh,iblkf,myminb,
     1 mymaxb,iblokind,omode,ibl)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/disktl)
INCLUDE(common/ga_file)
INCLUDE(common/ccntl)
INCLUDE(common/symci)
INCLUDE(common/cntrl)
      common/noddav/nwips,iblock1,iblock2,nwords,nread
      dimension a(il:il+ld-1,jl:jh)
      logical omode,oc
      dimension iblokind(*)
_IF(linux)
      external fget
_ENDIF
c
c ...
c
      if (omode) then
         idj = 1
      else
         idj = 2
      endif
      iworkblok=1
      icval=1
      ibuf=0
      ioff = 0
      jjj = 0
      kkk=1
      call search(ibl,numcip)
      itest=nwips
      do iworkblok=myminb,mymaxb
         nnq=iblokind(iworkblok)
         ibuf = ibuf + nnq
         istart=(iworkblok-1)*ga_file_block + 1
         kkk=nnq
         if (jjj.gt.0) then
             nbb=min(nnq,jjj)
             itest=itest-nbb
             if (istart.lt.igetend) then
c we have to move the block
                do iijk=jjj,1,-1
                   a(istart+iijk-1,idj)=a(istart2-1+iijk,idj) 
                enddo
                igetend=istart+jjj
                istart2=istart+nbb
                jjj=jjj-nbb
                kkk=nnq-nbb
                istart=istart+nbb
             else
                do iijk=1,nbb
                    a(istart-1+iijk,idj)=a(istart2-1+iijk,idj)
                enddo
                jjj=jjj-nbb
                kkk=nnq-nbb
                istart=istart+nbb
                if (jjj.ne.0) then
                   istart2=istart2+nbb
                endif
             endif
         endif
         if (jjj.eq.0.and.itest.gt.0) then
            call fget(a(istart,idj),mmm,numcip)
            jjj=mmm-kkk
            if (jjj.ne.0) istart2=istart+kkk
            igetend=istart+mmm
            itest=itest-(mmm-jjj)
         endif
      enddo
  60  if (nwips .ne. ibuf) call caserr('czdumpga error')
      return
      end

      subroutine create_ed19(nblok,iwr)
      implicit REAL  (a-h,o-z),integer  (i-n)        
INCLUDE(common/sizes)
INCLUDE(common/disktl)
INCLUDE(common/ga_file)
INCLUDE(common/corctl)
      common/remc_ed19/iblokk,ipoint,idiagoff,ibasp,dfrac_core,ncoretot
      numcz = 20
      ipoint = igmem_alloc(nblok)
      ipoint=ipoint-ibasp+1
      iga_to_use = ga_file_unit_to_ga( numcz )
      iblokk=nblok
      ga_file_size(iga_to_use) = 2*nblok*ga_file_block
      nga_use=(2*nblok*ga_file_block)/ipg_nnodes()+nblok
      dbest=(ncoretot-nga_use)*1.0d0/ncoretot*1.0d0
      write (iwr,6009)
 6009 format (/1x,42('*')/1x,'create global array for c/z-vectors (ed19)
     +'/1x,42('*'))
      write (iwr,6011) dfrac_core,dbest
 6011 format (/1x,42('*')/1x,'part used for Direct: ',f4.2,' best: ',
     +f4.2,/1x,42('*'))
      call search(1,numcz)
      call walltime(ddtime)
      write (iwr,6010)ddtime
 6010 format (/1x,40('+')/1x,'walltime in seconds ',F20.2,
     +/1x,40('+'))
      call flush(6)
c     call pg_distribution(ga_file_handle(iga_to_use),ipg_nodeid(),
c    1 il,ih,jl,jh)
c     if (2*nblok.gt.ga_file_size(iga_to_use)) then
c        write(iwr,'(A,3I10)')'Allocate more blocks for ed19',
c    1   2*nblok,ga_file_size(iga_to_use),iga_to_use
c        call caserr('Not enough space in global ed19 file')
c     endif
      return
      end

      subroutine create_index(q,iconf)
      implicit REAL  (a-h,o-z),integer  (i-n)        
INCLUDE(common/sizes)
INCLUDE(common/ccntl)
INCLUDE(common/symci)
INCLUDE(common/cntrl)
INCLUDE(common/ga_file)
INCLUDE(common/disktl)
      common/remc_ed19/iblokk,ipoint,idiagoff,ibasp
      REAL iconf
      dimension iconf(5,*),q(*)
      iga_to_use = ga_file_unit_to_ga( numcz )
      call pg_distribution(ga_file_handle(iga_to_use),ipg_nodeid(),
     1 il,ih,jl,jh)
      myminb = (il-1)/ga_file_block + 1
      call create_int(q(ipoint),iconf,myminb)
      return
      end

      subroutine create_int(iblok_index,iconf,myminb)
      implicit REAL  (a-h,o-z),integer  (i-n)        
INCLUDE(common/sizes)
INCLUDE(common/ccntl)
INCLUDE(common/symci)
INCLUDE(common/cntrl)
      common/remc_ed19/iblokk,ipoint,idiagoff
_IF(cray,t3e,i8)
      integer   iconf
_ELSEIF(convex,i8drct)
      integer *8   iconf
_ELSE
      REAL iconf
_ENDIF
      dimension iconf(5,*),iblok_index(*)
      iworkblok=1
      icval=1
      inval=icval
      isum=0
      idiagoff = 0
_IF(cray,convex,t3e,i8,i8drct)
      icdoub=icval+iconf(2,1)
_ELSE
      icdoub=icval+ipad(iconf(2,1))
_ENDIF
      if (nval.eq.0) inval=0
      if (ndoub.eq.0) icdoub=0
      icount=nstrt2
  10  if (iworkblok.eq.inval) then
         nncz=nval
      elseif (iworkblok.eq.icdoub) then
         nncz=ndoub
      else
  35     loop=icount
         if (icount.gt.nlast2) goto 60
         call upack2(iconf(2,loop),kc,jc)
         call upack2(iconf(3,loop),lamt,lams)
_IF(cray,convex,t3e,i8,i8drct)
         isi=iconf(4,loop)
_ELSE
         isi=ipad(iconf(4,loop))
_ENDIF
         nss=lams*norbsi(isi)
         nst=lamt*norbtr(isi)
         jc=icval+jc
         kc=icval+kc
         if (iworkblok.eq.jc) then
            nncz=nss
            if (nncz.eq.0) goto 37
            goto 45
         endif
  37     if (iworkblok.eq.kc) then
            nncz=nst
            if (nncz.eq.0) goto 40
            goto 45
         endif
  40     icount=icount+1
         if (icount.gt.nlast2) goto 60
         goto 35
      endif
  45  nnq=nncz
  50  nnq=min(nnq,511)
      iblok_index(iworkblok)=nnq
      if (iworkblok.lt.myminb) idiagoff=idiagoff+nnq
      isum=isum+nnq
      iworkblok = iworkblok+1
      nncz=nncz-511
      if (nncz .le. 0) then
         goto 10
      endif
      goto 45
  60  continue
      return
      end

      subroutine switch_oswed(iunit,omode)
      implicit REAL  (a-h,o-z),integer  (i-n)        
INCLUDE(common/sizes)
INCLUDE(common/nodeio)
      logical omode
      oswed3(iunit)=omode
      return
      end

      subroutine timer(text)
      implicit REAL (a-h,o-z), integer (i-n)
      character*(*)text
      iwr=6
      if (ipg_nodeid().ne.0)return
      call walltime(ddtime)
      write (iwr,6010)text,ipg_nodeid(),ddtime
 6010 format (/1x,40('+')/1x,'place ',A,' -- node:',I2,' -- walltime',
     +F10.2,/1x,40('+'))
      call flush(6)
      return
      end
      function ipg_whichnodeid()
      implicit REAL  (a-h,o-z)
_IF(parallel)
INCLUDE(common/parallel)
       ipg_whichnodeid = mod(icount_slb,ipg_nnodes())
_ELSE
       ipg_whichnodeid=0
_ENDIF
      return
      end

      subroutine remabcd(ic,iconf,cr,lump,lumq,lumr,oijkl,islast)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(cray,t3e,i8)
      integer iconf
_ELSEIF(convex,i8drct)
      integer *8 iconf
_ELSE
      REAL iconf
_ENDIF
c...
c... to build diagonal of hamiltonian
c... this routine used only if  =eigen= directive not presented
c...
c... used for both instore and outstore modes
c...
       logical luck,load,oijkl
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      dimension cr(*),iconf(5,*)
INCLUDE(common/disktl)
      common/spew/rmode(5)
      common/stoctl/khbill
INCLUDE(common/cntrl)
INCLUDE(common/table)
      common/moco/n128
      common/lsort/irebas(8),multt(8),multtt(8,8),scr(1000)
INCLUDE(common/orb)
INCLUDE(common/auxh)
      common/spcoef/nwwwww(9),iwwwww(9)
INCLUDE(common/symci)
INCLUDE(common/helpr)
INCLUDE(common/ccntl)
      equivalence (nornor,norbsi(1))
      equivalence (ibldub,index(198))
c
      call upack2(iconf(3,ic),lamt,lams)
_IF(cray,convex,t3e,i8,i8drct)
      isi=iconf(4,ic)
_ELSE
      isi=ipad(iconf(4,ic))
_ENDIF
      call upack2(iconf(2,ic),kc,jc)
      norsi=norbsi(isi)
      nortr=norbtr(isi)
      if(isi.eq.islast)goto 34
      islast=isi
c... retrieve diagonal of (ab/cd)  p  supermatrix
_IF(parallel)
c     call comed7
      call pried7
_ENDIF(parallel)
      call vclr(cr(lump+1),1,norsi)
      call vclr(cr(lumr+1),1,norsi*lams)
      call diagpq(cr(lump+1),cr(lumr+1),index(isi+9),norsi)
_IF(parallel)
      call pried7
_ENDIF(parallel)
      if(isi.ne.1)goto 204
c... normalize diagonal of totally symmetric block
      do 304 moop=1,nirr
      ne=norbe(moop)
      if(ne)404,304,404
404   idbas=lump+ibassi(moop,1)
      do 1 loop=1,ne
      joop=iky(loop+1)+idbas
    1 cr(joop)=cr(joop)*0.5d0
304   continue
c... retrieve diagonal of  (ab/cd)  q  supermatrix
204   continue
_IF(parallel)
c     call comed7
      call pried7
_ENDIF(parallel)
      call vclr(cr(lumq+1),1,nortr)
      call vclr(cr(lumr+1),1,nortr*lamt)
      call diagpq(cr(lumq+1),cr(lumr+1),index(isi+17),nortr)
_IF(parallel)
      call pried7
_ENDIF(parallel)
34    if(lams)54,44,54
54    lams1=lams+1
      if (oijkl) goto 41
      idbas=lumr
       icbas=khbill
      frodo=0.0d0
      do 64 loop=1,lams
_IF(cray,ibm,vax)
_IF1(c)cdir$ ivdep
      do 6 joop=1,norsi
    6 cr(icbas+joop)=cr(lump+joop)+frodo
_ELSE
      call vsadd(cr(lump+1),1,frodo,cr(icbas+1),1,norsi)
_ENDIF
      idbas=idbas+lams1
64    icbas=icbas+norsi
      call wrt3(cr(khbill+1),lams*norsi,index(199)+jc,numci)
41    if(odebug(40)) then
       kb=khbill
       if(oijkl)kb=lump
       call dbgvec('n-2sijkl',cr(kb+1),lams*norsi)
      endif
      if(lamt)44,84,44
44    lams1=lamt+1
      if (oijkl) goto 83
      idbas=lumr
      icbas=khbill
      frodo=0.0d0
      do 104 loop=1,lamt
_IF(cray,ibm,vax)
_IF1(c)cdir$ ivdep
      do 7 joop=1,nortr
    7 cr(icbas+joop)=cr(lumq+joop)+frodo
_ELSE
      call vsadd(cr(lumq+1),1,frodo,cr(icbas+1),1,nortr)
_ENDIF
      idbas=idbas+lams1
104   icbas=icbas+nortr
      call wrt3(cr(khbill+1),lamt*nortr,index(199)+kc,numci)
83    if(odebug(40)) then
       kb=khbill
       if (oijkl)kb=lumq
       call dbgvec('n-2tijkl',cr(kb+1),lamt*nortr)
      endif
84    continue
      return 
      end

      function ipg_stdub(ner,ner2,isr,ll)
      implicit REAL (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/orb)
INCLUDE(common/symci)
      common/rem_iabc/myext(62,nd200),myibass3(8,nd200),myniky(nd200,8),
     *mynorbe(nd200,8)
      ipg_stdub=0
      if (ner.eq.0) return
      istart=1
      do i=1,isr-1
        istart=istart+norbe(i)
      enddo
      iii=0
      do i=istart,istart+ner2-1
      iii=iii+1
         if (myext(ll,nintp1-1+i).eq.1) then
             ipg_stdub=iii
             return
         endif
      enddo
      if (ipg_stdub.eq.0) call caserr('0 in ipg_stdub')
      return
      end
_ENDIF
