_IF(win95)
*$pragma aux openc "!_" parm (reference,value,reference,reference,reference)
_ENDIF
      subroutine bckchn(mchain,itap61,buf,ibuf,nszbf,maxbkt,
     .                  maxval,rints,jout,
     1           nxbkt,newlen,npercl,itap91,lpercl,irew)
      implicit REAL (a-h,o-z)
INCLUDE(common/prnprn)
      logical irew
      integer ibuf(nszbf*2),mchain(maxbkt)
      REAL buf(nszbf),rints(newlen*npercl)
      if(irew) call srew(itap61)
      if(odebug(33)) write(6,*) ' in bckchn, nxbkt = ',nxbkt
c
      do 103 kl=1,nxbkt
      call vclr(rints,1,newlen*npercl)
      isctr=mchain(kl)
      if(odebug(33)) write(*,*) '  in bckch reading from sector ',isctr
 1111 call rread(itap91,buf,intowp(nszbf),isctr)
c     write (jout,*) ' isctr=',isctr
c     write(*,*) ' ibuf ',ibuf(1),ibuf(2)
c     write(*,*) ' buf ',buf(1)
c     write (jout,*) ' read tape 77'
_IF1()ctjl  isctr=ibuf(1)
      call setmbf(isctr,ibuf(1))
      if(odebug(33)) write (jout,*) ' isctr=',isctr
      call setmbf(mbuf,ibuf(2))
      if(odebug(33)) write (jout,*) ' mbuf=',mbuf,ibuf(2)
c     write (jout,*) ' maxval=',maxval
      ivoff=(maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
_IF1()ctjl  ivoff=(maxval+3)/intowp(1)
c     iboff=itri*(ibkt-1)
      if(odebug(33)) write (jout,*) ' ivoff=',ivoff
      do 101 i=1,mbuf
c        write (jout,*) ' i=',i
_IF1()cj       call setmbf(i3,ibuf(2+i))
_IF1()ctjl     i3=ibuf(i+2)
      if(odebug(33)) write(6,*) ' i,ibuf(i+2),buf(ivoff+i) ',i,
     +   ibuf(i+2),buf(ivoff+i)
         rints(ibuf(i+2))=buf(ivoff+i)
c
c     if(i1.eq.i2) then
c     write(*,*) '   in bckchn;  ij,kl,val = ',i1,i2,rints(ind)
c     end if
c
_IF1()ctjl
c     if(ind.ge.1.or.ind.le.10) then
c     write(*,*) ' in bckch:ij,kl,val=',ij,kl,rints(ind)
c     end if
_IF1()ctjl
c        write (jout,*) i1,i2,ind,ind-iboff,rints(ind-iboff)
  101 continue
      call vclr(buf,1,nszbf)
      if (isctr.ne.0) goto 1111
      ipercl=npercl
      if(kl.eq.nxbkt) ipercl=lpercl
      do 105 ip=1,ipercl
       test = 0.0d0
      do 106 jki=1,newlen
106   test=test+dabs(rints((ip-1)*newlen+jki))
      if(odebug(33)) write(6,*) ' test , npercl, newlen ',test,
     +               npercl,newlen
105   call swrit(itap61,rints((ip-1)*newlen+1),intowp(newlen))
103   continue
c
c     write (jout,*) ' exiting backchain'
      return
      end
      subroutine bktdmp(bkt,length,idisk)
c
      implicit REAL (a-h,o-z)
      dimension bkt(length)
      call swrit(idisk,bkt,intowp(length))
      return
      end
      subroutine bktdmx(bkt,length,idisk,ichan)
c
      implicit REAL (a-h,o-z)
INCLUDE(common/prnprn)
      dimension bkt(length)
      if(odebug(33)) write(*,*) ' just enetered bktdmx '
      call rgetsa(idisk,ichan)
      call swrit(idisk,bkt,intowp(length))
      if(odebug(33)) write(*,*) ' leaving bktdmx  '
      return
      end
      subroutine fmatx(buf,ibuf,lenbf,ivoff,h,ntri,eig,nbf,escf,maxval,
     1                 itap60,itap62,itap63,ni,nj,nk,nl,no,enuc)
      implicit integer (a-z)
      REAL buf(lenbf),h(ntri),eig(nbf),escf,enuc,a0
      integer ibuf(*),ni(maxval),nj(maxval),nk(maxval),nl(maxval)
c
      data itemp,a0 /255,0.0d0/
c
c  first we must read in the c integrals (file itap62)
c
_IF1(ct)cdir$ block
      call vclr(eig,1,nbf)
      call srew(itap62)
 10   call tit_sread(itap62,ibuf,intowp(lenbf))
      num = ibuf(2)
c
      do 20 ind = 1,num
_IF(cray)
      ni(ind) = shiftr(ibuf(2+ind),24)
      nj(ind) = and(itemp,shiftr(ibuf(2+ind),16))
      nk(ind) = and(itemp,shiftr(ibuf(2+ind),8))
      nl(ind) = and(itemp,ibuf(2+ind))
_ELSE
      ni(ind) = ishft(ibuf(2+ind),-24)
      nj(ind) = IAND32(itemp,ishft(ibuf(2+ind),-16))
      nk(ind) = IAND32(itemp,ishft(ibuf(2+ind),-8))
      nl(ind) = IAND32(itemp,ibuf(2+ind))
_ENDIF
 20   continue
c
      do 30 ind = 1,num
      if(ni(ind).eq.nj(ind).and.nk(ind).eq.nl(ind)) 
     1 eig(ni(ind)) = eig(ni(ind)) + buf(ivoff+ind) + buf(ivoff+ind)
 30   continue
c
      if(ibuf(1).eq.0) go to 10
c
c  now read in the d integrals and add their contribution
c
      call srew(itap63)
 40   call tit_sread(itap63,ibuf,intowp(lenbf))
      num = ibuf(2)
      do 50 ind = 1,num
_IF(cray)
      ni(ind) = shiftr(ibuf(2+ind),24)
      nj(ind) = and(itemp,shiftr(ibuf(2+ind),16))
      nk(ind) = and(itemp,shiftr(ibuf(2+ind),8))
      nl(ind) = and(itemp,ibuf(2+ind))
_ELSE
      ni(ind) = ishft(ibuf(2+ind),-24)
      nj(ind) = IAND32(itemp,ishft(ibuf(2+ind),-16))
      nk(ind) = IAND32(itemp,ishft(ibuf(2+ind),-8))
      nl(ind) = IAND32(itemp,ibuf(2+ind))
_ENDIF
 50   continue
c
      do 60 ind = 1,num
      if(ni(ind).eq.nk(ind).and.nj(ind).eq.nl(ind)) 
     1 eig(ni(ind)) = eig(ni(ind)) - buf(ivoff+ind)
 60   continue
c
      if(ibuf(1).eq.0) go to 40
c
c  now read in the a integrals - the most complicated set
c
      call srew(itap60)
 70   call tit_sread(itap60,ibuf,intowp(lenbf))
      num = ibuf(2)
      do 80 ind = 1,num
_IF(cray)
      ni(ind) = shiftr(ibuf(2+ind),24)
      nj(ind) = and(itemp,shiftr(ibuf(2+ind),16))
      nk(ind) = and(itemp,shiftr(ibuf(2+ind),8))
      nl(ind) = and(itemp,ibuf(2+ind))
_ELSE
      ni(ind) = ishft(ibuf(2+ind),-24)
      nj(ind) = IAND32(itemp,ishft(ibuf(2+ind),-16))
      nk(ind) = IAND32(itemp,ishft(ibuf(2+ind),-8))
      nl(ind) = IAND32(itemp,ibuf(2+ind))
_ENDIF
 80   continue
c
      do 90 ind = 1,num
      if(ni(ind).eq.nj(ind).and.nk(ind).eq.nl(ind)) then
      eig(ni(ind)) = eig(ni(ind)) + buf(ivoff+ind) + buf(ivoff+ind)
      eig(nk(ind)) = eig(nk(ind)) + buf(ivoff+ind) + buf(ivoff+ind)
      end if
 90   continue
c
      do 100 ind = 1,num
      if(ni(ind).eq.nk(ind).and.nj(ind).eq.nl(ind)) then
      eig(ni(ind)) = eig(ni(ind)) - buf(ivoff+ind)
      eig(nj(ind)) = eig(nj(ind)) - buf(ivoff+ind)
      end if
 100  continue
c
      do 110 ind = 1,num
      if(ni(ind).eq.nj(ind).and.nk(ind).eq.nl(ind).and.ni(ind).eq.
     1   nk(ind)) eig(ni(ind)) = eig(ni(ind)) - buf(ivoff+ind)
 110  continue
c
      if(ibuf(1).eq.0) go to 70
c
_IF1()ctjl      write(*,*) ' one electron integrals ' 
_IF1()ctjl      call print(h,ntri,nbf,6)
      escf = enuc
      do 130 i = 1,no
      escf = escf - eig(i)
 130  continue
c
      do 120 i = 1,nbf
      ii = (i+1)*i/2 
      eig(i) = eig(i) + h(ii)
 120  continue
c
      do 140 i = 1,no
      escf = escf + eig(i) + eig(i)
 140  continue
c
      write(*,435)escf,(eig(i),i=1,nbf)
 435  format(' scf energy ',e20.12,/,
     &       ' orbital energies ',/,(5e15.5))
c
_IF1(ct)cdir$ vector
      return
      end
      integer FUNCTION IADTWP(N)
_IF(cray,i8)
      IADTWP=N
_ELSE
      IADTWP=(N+3)/2
_ENDIF
      RETURN
      END
      subroutine initt2(buf,ibuf,nszbf,length,nbf,
     .     itap60,jout,indxvl,ntriv,eigval,
     1     norbo,iofosy,norsq,nirr,no,indxl,
     1    maxval,orbsym,t2,it69,ndimt1,ndimt2,ioff)
      implicit REAL (a-h,o-z)
INCLUDE(common/prnprn)
      integer orbsym
      dimension buf(length),eigval(nbf),ioff(nbf),
     .          ibuf(length*2),orbsym(nbf),t2(2)
      dimension indxvl(ntriv),norbo(nirr),indxl(no*no),
     1iofosy(nirr+1),norsq(nirr)
      data itemp /255/
_IF1()   11 format ('                                                ')
_IF1()   12 format ('************************************************')
_IF1()   19 format (4i3,f20.12)
      if(odebug(33)) write(*,*) '   in initt2 ',nirr
      if(odebug(33)) write(6,*) ' iofosy ',iofosy
      call vclr(t2,1,ndimt2)
      call srew(it69)
      call swrit(it69,t2,intowp(ndimt1))
c
      intlen=(intowp(nszbf)-2)/intowp(1)
      maxval=intowp(intlen)/(1+intowp(1))
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
      ibflen=(intowp(length)-2)/intowp(1)
      maxbuf=intowp(ibflen)/(1+intowp(1))
      iboff = (maxbuf+2)/intowp(1)
      if(intowp(iboff).ne.(maxbuf+2)) iboff = iboff + 1
c      write (jout,*) ' max values read in , maxbuf = ',maxbuf
      call srew(itap60)
  111 call vclr(buf,1,length)
      call tit_sread(itap60,buf,intowp(length))
c
      call setmbf(iflg,ibuf(1))
c     write (jout,*) ' iflg=',iflg
      call setmbf(mbuf,ibuf(2))
c     write (jout,*) ' mbuf=',mbuf
      do 101 ii=1,mbuf
         call setmbf(ijkl,ibuf(2+ii))
         i=ishft(ijkl,-24)
         j=IAND32(itemp,ishft(ijkl,-16))
         k=IAND32(itemp,ishft(ijkl,-8))
         l=IAND32(itemp,ijkl)
         rint=buf(iboff+ii)
         i=i-no
         k=k-no
_IF1()c        write (jout,11)
_IF1()c        write (jout,12)
_IF1()c        write (jout,*) ' *********  ints as read into sortj *******'
_IF1()c        write (jout,11)
_IF1()c        write (jout,12)
_IF1()c        write (jout,19) i,j,k,l,rint
c
c
c
         ik=ioff(i)+k
         jl=(j-1)*no+l
         iposik=indxvl(ik)
         iposjl=indxl(jl)
         jlsym=IXOR32(orbsym(j),orbsym(l))+1
_IF1()cj         iksym=IXOR32(orbsym(i+no),orbsym(k+no))+1
         iadr=(iposik-1)*norsq(jlsym)+iposjl+iofosy(jlsym)
_IF1()cj       t2(iadr)=rint/(eigval(j)+eigval(l)-eigval(k+no)-eigval(i+no))
         t2(iadr)=rint
         if(i.eq.k) then
            jl=(l-1)*no+j
            iposjl=indxl(jl)
            iadr=(iposik-1)*norsq(jlsym)+iposjl+iofosy(jlsym)
            t2(iadr)=rint
         end if
c
  101 continue
      if (iflg.eq.0) goto 111
      if(odebug(33)) then
       write(6,*) ' t2 finished ? '
        do 1010 i=1,iofosy(nirr+1)
        write(6,*) ' i,t2(i) ',i,t2(i)
 1010   continue
      endif
_IF1()cj      call swrit(it69,t2,intowp(ndimt2))
c
      return
      end
      integer FUNCTION INTOWP(N)
_IF(cray,i8)
      INTOWP=N
_ELSE
      INTOWP=2*N
_ENDIF
      RETURN
      END

      block data t_flname_data
INCLUDE(common/t_flname)
      DATA FNAME/'ft01','ft02','ft03','ft04','ft05','ft06','ft07',
     *  'ft08','ft09',
     1  'ft10','ft11','ft12','ft13','ft14','ft15','ft16','ft17','ft18',
     2  'ft19','ft20','ft21','ft22','ft23','ft24','ft25','ft26','ft27',
     3  'ft28','ft29','ft30','ft31','ft32','ft33','ft34','ft35','ft36',
     4  'ft37','ft38','ft39','ft40','ft41','ft42','ft43','ft44','ft45',
     5  'ft46','ft47','ft48','ft49','ft50','ft51','ft52','ft53','ft54',
     6  'ft55','ft56','ft57','ft58','ft59','ft60','ft61','ft62','ft63',
     7  'ft64','ft65','ft66','ft67','ft68','ft69','ft70','ft71','ft72',
     8  'ft73','ft74','ft75','ft76','ft77','ft78','ft79','ft80','ft81',
     9  'ft82','ft83','ft84','ft85','ft86','ft87','ft88','ft89','ft90',
     1  'ft91','ft92','ft93','ft94','ft95','ft96','ft97','ft98','ft99',
     1  'help'/
      end
      SUBROUTINE RFILE(ITAPE)
C
C USE THIS FOR ALL UNFORMATTED FILES
C
      implicit integer (A-Z)
      character*8 FNAM
C
INCLUDE(common/t_pointr)
INCLUDE(common/t_sect)
INCLUDE(common/t_flname)
C
INCLUDE(common/utilc)
c
C
C
_IF(cray)
      SECTOR = 1024
      BLOCKS = 20
      ISTATS = 0
_ELSE
      SECTOR = 512
_ENDIF
      STATUS = 1
      SIZE = 0
      FNAM = FNAME(ITAPE)
C
      if(ooprnt) WRITE(*,*) ' OPENING FILE ',FNAM,'  ON UNIT ',ITAPE 
_IF(cray)
      CALL WOPEN(ITAPE,20,ISTATS,IERR)
      IF(IERR.NE.0) THEN
        WRITE(*,*) ' ERROR ENCOUNTERED IN OPENING FILE',ITAPE
        WRITE(*,*) ' IERR = ',IERR
        call caserr('open error in rfile')
      END IF
_ELSE
      CALL OPENC(ITAPE,FNAM,SIZE,STATUS,4)
_ENDIF
      if(ooprnt) WRITE(*,*) ' OPENED FILE ',FNAM,'  ON UNIT ',ITAPE
C
C PLACE THE POINTER OF ITAPE AT THE BEGINNING
C
      POS(ITAPE) = 1
C
C RETURN
C
      RETURN
      END
      SUBROUTINE WWRITW(ITAPE,ARRAY,NLEN,FWORD,LWORD)
C
      implicit integer (A-Z)
C
INCLUDE(common/t_pointr)
INCLUDE(common/t_sect)
INCLUDE(common/utilc)
C
      DIMENSION ARRAY(NLEN)
C
      if(ooprnt)
     +    WRITE(*,*) ' IN WWRITW, FWORD,NLEN,ITAPE = ',FWORD,NLEN,ITAPE
_IF(cray)
      CALL PUTWA(ITAPE,ARRAY,FWORD,NLEN,IERR)
      IF(IERR.NE.0) THEN
        WRITE(*,*) ' ERROR IN WWRITW FOR FILE',ITAPE
        WRITE(*,*) ' IERR = ',IERR
        call caserr('write error in wwritw')
      END IF
_ELSE
      CALL WRABSF(ITAPE,ARRAY,NLEN,FWORD)
_ENDIF
C
C UPDATE LWORD AND THEN RETURN
C
      LWORD = FWORD + NLEN
      POS(ITAPE) = LWORD
C
      RETURN
      END
      SUBROUTINE WREADW(ITAPE,ARRAY,NLEN,FWORD,LWORD)
C
      implicit integer (A-Z)
C
INCLUDE(common/t_pointr)
INCLUDE(common/t_sect)
INCLUDE(common/t_files)
INCLUDE(common/utilc)
C
      DIMENSION ARRAY(NLEN)
C
      if(ooprnt)
     + WRITE(*,*) ' IN WREADW, FWORD,NLEN,ITAPE = ',FWORD,NLEN,ITAPE
_IF(cray)
      CALL GETWA(ITAPE,ARRAY,FWORD,NLEN,IERR)
      IF(IERR.NE.0) THEN
        WRITE(*,*) ' ERROR IN WREADW FOR FILE',ITAPE
        WRITE(*,*) ' IERR = ',IERR
        call caserr('read error in wreadw')
      END IF
_ELSE
      CALL RDABSF(ITAPE,ARRAY,NLEN,FWORD)
_ENDIF
C
C UPDATE LWORD AND THEN RETURN
C
      LWORD = FWORD + NLEN
      POS(ITAPE) = LWORD
      if(ooprnt)WRITE(*,*) ' LWORD ',LWORD
C
      RETURN
      END
      SUBROUTINE RGETSA(ITAPE,IADR)
C
C THIS ROUTINE RETURNS THE CURRENT SECTOR ADDRESS FOR FILES USING
C THE SREAD/SWRIT AND/OR RREAD/RWRIT IO ROUTINES.
C
      implicit integer (A-Z)
C
INCLUDE(common/t_pointr)
INCLUDE(common/t_sect)
C
      IPOS = POS(ITAPE)
      IREC = (IPOS-1)/SECTOR + 1
      TEST = SECTOR*(IREC-1) + 1
      ierr=0
      IF(IPOS.NE.TEST.OR.IERR.NE.0) THEN
         WRITE(*,*) ' ERROR ENCOUNTERED IN RGETSA FOR FILE ',ITAPE
         WRITE(*,*) ' IERR,IPOS,TEST = ',IERR,IPOS,TEST
         call caserr('RGETSA data error')
      END IF
C
      IADR=IREC
C
      RETURN
      END
      SUBROUTINE RSETSA(ITAPE,IADR)
C
C THIS ROUTINE SETS THE POINTER OF A UTILITIES FILE TO THE DESIRED
C SECTOR ADDRESS.
C
      implicit integer (A-Z)
C
INCLUDE(common/t_pointr)
INCLUDE(common/t_sect)
C
      IPOS = SEC2I(IADR-1) + 1
      POS(ITAPE) = IPOS
_IF1()      ierr=0
_IF1()      IF(IERR.NE.0) THEN
_IF1()         WRITE(*,*) ' ERROR ENCOUNTERED IN RSETSA FOR FILE ',ITAPE
_IF1()         WRITE(*,*) ' IERR,IPOS,IADR = ',IERR,IPOS,IADR
_IF1()         call caserr('RSETSA data error')
_IF1()      STOP
_IF1()      END IF
C
      RETURN
      END
      SUBROUTINE RREAD(ITAPE,ARRAY,NLEN,IREC)
C
C THIS ROUTINE READS NLEN INTEGER WORDS INTO ARRAY STARTING AT
C SECTOR IREC.
C
      implicit integer (A-Z)
C
INCLUDE(common/t_pointr)
INCLUDE(common/t_sect)
INCLUDE(common/t_iobuff)
INCLUDE(common/utilc)
C
      DIMENSION ARRAY(NLEN)
C
      FWORD = (IREC-1)*SECTOR + 1
      if(ooprnt) WRITE(*,*) ' IN RREAD FOR FILE',ITAPE
      if(ooprnt) WRITE(*,*) ' IREC,FWORD,NLEN = ',IREC,FWORD,NLEN
_IF(cray)
       CALL GETWA(ITAPE,ARRAY,FWORD,NLEN,IERR)
       IF(IERR.NE.0) THEN
          WRITE(*,*) ' ERROR ENCOUNTERED IN RREAD FOR FILE ',ITAPE
          WRITE(*,*) ' IERR = ',IERR
          call caserr('RREAD data error')
       END IF
_ELSE
       CALL RDABSF(ITAPE,ARRAY,NLEN,FWORD)
_ENDIF
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      TEST = NLEN/SECTOR
      IF(SECTOR*TEST .NE. NLEN) TEST = TEST + 1
      TESTL = SECTOR*TEST
      POS(ITAPE) = FWORD + TESTL
C
      RETURN
      END
      SUBROUTINE RWRIT(ITAPE,ARRAY,NLEN,IREC)
C
C THIS ROUTINE WRITES NLEN INTEGER WORDS FROM ARRAY TO FILE ITAPE
C STARTING AT SECTOR IREC.
C
      implicit integer (A-Z)
C
INCLUDE(common/t_pointr)
INCLUDE(common/t_sect)
INCLUDE(common/t_iobuff)
INCLUDE(common/utilc)
C
      DIMENSION ARRAY(NLEN)
C
      FWORD = (IREC-1)*SECTOR + 1
      if(ooprnt)WRITE(*,*) ' IN RWRIT FOR FILE',ITAPE
      if(ooprnt)WRITE(*,*) ' IREC,NLEN = ',IREC,NLEN
_IF(cray)
      CALL PUTWA(ITAPE,ARRAY,FWORD,NLEN,IERR)
      IF(IERR.NE.0) THEN
       WRITE(*,*) ' ERROR ENCOUNTERED IN RWRIT FOR FILE ',ITAPE
       WRITE(*,*) ' IERR = ',IERR
       call caserr('write error in rwrit')
      END IF
_ELSE
      CALL WRABSF(ITAPE,ARRAY,NLEN,FWORD)
_ENDIF
C
C WE MUST FILL THE ENTIRE BUFFER WITH ZEROS
C
      EXTRA = MOD(NLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      FWORD2 = FWORD + NLEN
      NW = SECTOR - EXTRA
_IF(cray)
_IF1()      CALL PUTWA(ITAPE,BUFF,FWORD2,NW,IERR)
_IF1()      IF(IERR.NE.0) THEN
_IF1()       WRITE(*,*) ' ERROR ENCOUNTERED WHILE POSITIONING FILE POINTER'
_IF1()       WRITE(*,*) ' FOR FILE ',ITAPE,'   IN RWRIT,  IERR = ',IERR
_IF1()      call caserr('file positioning error in RWRIT')
_IF1()      END IF
_ELSE
      CALL WRABSF(ITAPE,BUFF,NW,FWORD2)
_ENDIF
      END IF
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      TEST = NLEN/SECTOR
      IF(SECTOR*TEST .NE. NLEN) TEST = TEST + 1
      TESTL = SECTOR*TEST
      POS(ITAPE) = FWORD + TESTL
C
      RETURN
      END
      SUBROUTINE TIT_SREAD(ITAPE,ARRAY,NLEN)
C
C THIS ROUTINE READS NLEN INTEGER WORDS FROM ITAPE INTO ARRAY
C STARTING AT THE CURRENT POINTER LOCATION.
C
      implicit integer (A-Z)
C
INCLUDE(common/t_pointr)
INCLUDE(common/t_iobuff)
INCLUDE(common/t_sect)
INCLUDE(common/utilc)
C
      DIMENSION ARRAY(NLEN)
C
      FWORD = POS(ITAPE)
      if(ooprnt)
     +   WRITE(*,*) ' IN SREAD, ITAPE,NLEN,FWORD = ',ITAPE,NLEN,FWORD
_IF(cray)
      CALL GETWA(ITAPE,ARRAY,FWORD,NLEN,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED IN SREAD FOR FILE ',ITAPE
      WRITE(*,*) ' IERR = ',IERR
      call caserr('sread error encountered')
      END IF
_ELSE
      CALL RDABSF(ITAPE,ARRAY,NLEN,FWORD)
_ENDIF
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      TEST = NLEN/SECTOR
      IF(SECTOR*TEST .NE. NLEN) TEST = TEST + 1
      TESTL = SECTOR*TEST
      POS(ITAPE) = FWORD + TESTL
C
      RETURN
      END
      SUBROUTINE SWRIT(ITAPE,ARRAY,NLEN)
C
C THIS ROUTINE WRITES NLEN INTEGER WORDS FROM ARRAY TO FILE ITAPE
C STARTING AT THE CURRENT POINTER LOCATION.
C
      implicit integer (A-Z)
C
INCLUDE(common/t_pointr)
INCLUDE(common/t_iobuff)
INCLUDE(common/t_sect)
INCLUDE(common/utilc)
C
      DIMENSION ARRAY(NLEN)
C
      FWORD = POS(ITAPE)
      if(ooprnt)
     +  WRITE(*,*) ' IN SWRIT: ITAPE,NLEN,FWORD = ',ITAPE,NLEN,FWORD
_IF(cray)
      CALL PUTWA(ITAPE,ARRAY,FWORD,NLEN,IERR)
      IF(IERR.NE.0) THEN
       WRITE(*,*) ' ERROR ENCOUNTERED IN SWRIT FOR FILE ',ITAPE
       WRITE(*,*) ' IERR = ',IERR
       call caserr('error encountered in routine swrit')
      END IF
_ELSE
      CALL WRABSF(ITAPE,ARRAY,NLEN,FWORD)
_ENDIF
C
C WE MUST FILL THE ENTIRE BUFFER WITH ZEROS
C
      EXTRA = MOD(NLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
CTJL  WRITE(*,*) ' IN SWRIT: extra ne to 0; extra = ',extra
      FWORD2 = FWORD + NLEN
      NW = SECTOR - EXTRA
_IF(cray)
      CALL PUTWA(ITAPE,BUFF,FWORD2,NW,IERR)
      IF(IERR.NE.0) THEN
       WRITE(*,*) ' ERROR ENCOUNTERED WHILE POSITIONING FILE POINTER '
       WRITE(*,*) ' FOR FILE ',ITAPE,'   IN RWRIT,  IERR = ',IERR
       call caserr('search error in swrit')
      END IF
_ELSE
      CALL WRABSF(ITAPE,BUFF,NW,FWORD2)
_ENDIF
      END IF
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      TEST = NLEN/SECTOR
      IF(SECTOR*TEST .NE. NLEN) TEST = TEST + 1
      TESTL = SECTOR*TEST
      POS(ITAPE) = FWORD + TESTL
C
      RETURN
      END
      SUBROUTINE SREW(ITAPE)
C
      implicit integer (A-Z)
C
INCLUDE(common/t_pointr)
C
C PLACE THE POINTER OF ITAPE AT THE BEGINNING
C
C
C UPDATE THE POSITION ARRAY
C
      POS(ITAPE) = 1
C
C RETURN
C
      RETURN
      END
      SUBROUTINE RCLOSE(ITAPE,JCODE)
C
      implicit integer (A-Z)
C
_IF1(c)       DIMENSION DN(100)
_IF1(c)       CHARACTER*20 FNAM
_IF1()      CHARACTER*20 F61,F77,F78
C
INCLUDE(common/t_flname)
INCLUDE(common/utilc)
C
      if(ooprnt)WRITE(*,*) ' CLOSING FILE',ITAPE,'  JCODE = ',JCODE
      IF(JCODE.NE.3.AND.JCODE.NE.4) THEN
        WRITE(*,*) ' INVALID JCODE IN RCLOSE,  JCODE = ',JCODE
        WRITE(*,*) ' FILE ',ITAPE,'  CLOSED AND SAVED.'
        JCODE = 3
      END IF
      IF(ITAPE.EQ.6) STOP ' YOU CANNOT CLOSE A FILE ON UNIT 6'
      if (jcode.eq.3)then
_IF(cray)
      CALL WCLOSE(ITAPE,IERR)
      IF(IERR.NE.0) THEN
        WRITE(*,*) ' ERROR ENCOUNTERED CLOSING FILE',ITAPE
        WRITE(*,*) ' IERR,JCODE = ',IERR,JCODE
        CALL MABORT
      END IF
_ELSE
      CALL CLOSEC(ITAPE)
_ENDIF
C
      else if(JCODE .EQ. 4) THEN
_IF(cray)
_IF1()      CALL RELEASE(STAT,'DN'L,DN(ITAPE))
      STAT=0
      CALL WCLOSE(ITAPE,STAT)
      FNAM = FNAME(ITAPE)
_IF1()CTJL  IF(ITAPE.EQ.61) FNAM = F61
_IF1()CTJL  IF(ITAPE.EQ.77) FNAM = F77
_IF1()CTJL  IF(ITAPE.EQ.78) FNAM = F78
      IF(STAT.NE.0) THEN
        WRITE(6,*) '  PROBLEMS DELETING FILE ',ITAPE,'  STAT = ',STAT
        WRITE(6,*) ' FNAM = ',FNAM
        CALL MABORT
      END IF
_ELSE
      CALL DFILEC(ITAPE)
_ENDIF
      END IF
C
      RETURN
      END
      subroutine tit_izo(a,n)
c
      integer a(n)
c
      do 1 i=1,n
           a(i)=0
    1 continue
c
      return
      end
      SUBROUTINE MABORT
C
C     IABORT = 999
C     REWIND IABORT
C
      CALL CASERR('ABORT IN CCSD SORT MODULE')
C
      RETURN
      END
      subroutine nasai(doc,uoc,flov,orbsym,ptocc,nirred,nbf,no,nv,
     1                 idu,ktri,ltri,ksq,lsq,ptsq,bkof1,bkof2,tfile,
     2                 ntris,nsqs,lenbf,ivoff,itype)
      implicit integer (a-z)
      integer doc(nirred),uoc(nirred),flov(nirred,4),orbsym(nbf)
      integer ptocc(nbf),idu(nbf),ktri(ntris),ltri(ntris),ksq(nsqs)
      integer lsq(nsqs),ptsq(nirred,nirred),bkof1(16),bkof2(16)
      integer tfile(16),itype(nbf)
INCLUDE(common/prnprn)
c
_IF1(ct)cdir$ block
      do 10 irr = 1,nirred
      flov(irr,1) = 0
      flov(irr,2) = -1
      flov(irr,3) = 0
      flov(irr,4) = -1
 10   continue
c
      icnt = 0
      do 20 irr = 1,nirred
      noirr = doc(irr)
      isym = irr - 1
      icnth = icnt + 1
      do 30 i = 1,noirr
      icnt = icnt + 1
      orbsym(icnt) = isym
 30   continue
      if(icnth-1.ne.icnt) then
      flov(irr,1) = icnth
      flov(irr,2) = icnt
      end if
 20   continue
c
      icnt = no
      do 40 irr = 1,nirred
      uoirr = uoc(irr)
      isym = irr - 1
      icnth = icnt + 1
      do 50 i = 1,uoirr
      icnt = icnt + 1
      orbsym(icnt) = isym
 50   continue
      if(icnth-1.ne.icnt) then
      flov(irr,3) = icnth
      flov(irr,4) = icnt
      end if
 40   continue
       if(odebug(33)) then
        write(*,*) ' in nasai; doc,uoc ',doc,uoc
        write(*,*) ' flov ',flov
       endif
c
c  now construct ptocc
c
      icnt = 0
      do irr=1,nirred
       do itst=1,no
        if (itype(itst).eq.irr)then
         icnt=icnt+1
         ptocc(itst)=icnt
        endif
       enddo
      enddo
      do irr=1,nirred
       do itst=no+1,no+nv
        if (itype(itst).eq.irr)then
         icnt=icnt+1
         ptocc(itst)=icnt
        endif
       enddo
      enddo
c
c  set up idu array - in cc order
c
      call tit_izo(idu,nbf)
      fuoc = no + 1
      do 100 i = fuoc,nbf
      idu(i) = 1
 100  continue
c
c  set up ktri and ltri arrays - input pitzer kl and get out cc k or l
c
      icnt = 0
      fmoir = 1
      do 120 irr = 1,nirred
      nmoir = doc(irr) + uoc(irr)
      lmoir = fmoir + nmoir - 1
      do 140 i = fmoir,lmoir
      do 150 j = fmoir,i
      icnt = icnt + 1
      ktri(icnt) = max(ptocc(i),ptocc(j))
      ltri(icnt) = min(ptocc(i),ptocc(j))
 150  continue
 140  continue
      fmoir = lmoir + 1
 120  continue
c
c  set up ksq and lsq - given a pitzer kl get back a cc k and l
c  also set up ptsq
c
      call tit_izo(ptsq,nirred*nirred)
      icnt = 0
      fmoir = 1
      do 160 irr = 1,nirred
      nmoir = doc(irr) + uoc(irr)
      if(nmoir.eq.0) go to 161
      lmoir = fmoir + nmoir - 1
      fmojr = 1
      do 170 jrr = 1,irr-1
      nmojr = doc(jrr) + uoc(jrr)
      lmojr = fmojr + nmojr - 1
_IF1()ctjl      if(irr.ne.jrr) then
      ptsq(irr,jrr) = icnt + 1
      ptsq(jrr,irr) = icnt + 1
      do 190 i = fmoir,lmoir
      do 200 j = fmojr,lmojr
      icnt = icnt + 1
      ksq(icnt) = max(ptocc(i),ptocc(j))
      lsq(icnt) = min(ptocc(i),ptocc(j))
 200  continue
 190  continue
_IF1()ctjl      end if
      fmojr = lmojr + 1
 170  continue
      fmoir = lmoir + 1
 161  continue
 160  continue
c
      if(odebug(33)) then
         write(*,*) ' ptocc  ',ptocc
         write(*,*) ' orbsym ',orbsym
         write(*,*) ' ktri ',ktri
         write(*,*) ' ltri ',ltri
         write(*,*) ' ksq ',ksq
         write(*,*) ' lsq ',lsq
         write(*,*) ' ptsq ',ptsq
      endif
c
c  now set up bkof1, bkof2 and tfile
c
      call tit_izo(bkof1,16)
      call tit_izo(bkof2,16)
      call tit_izo(tfile,16)
c
      bkof1(1)  = ivoff
      bkof1(16) = ivoff +   lenbf
      bkof1(13) = ivoff + 2*lenbf
      bkof1(4)  = ivoff + 2*lenbf
      bkof1(11) = ivoff + 3*lenbf
      bkof1(9)  = ivoff + 4*lenbf
      bkof1(3)  = ivoff + 4*lenbf
      bkof1(12) = ivoff + 5*lenbf
      bkof1(15) = ivoff + 5*lenbf
c
      bkof2(1)  = 0
      bkof2(16) =   intowp(lenbf)
      bkof2(13) = 2*intowp(lenbf)
      bkof2(4)  = 2*intowp(lenbf)
      bkof2(11) = 3*intowp(lenbf)
      bkof2(9)  = 4*intowp(lenbf)
      bkof2(3)  = 4*intowp(lenbf)
      bkof2(12) = 5*intowp(lenbf)
      bkof2(15) = 5*intowp(lenbf)
c
      tfile(1)  = 60
      tfile(16) = 61
      tfile(13) = 62
      tfile(4)  = 62
      tfile(11) = 63
      tfile(9)  = 64
      tfile(3)  = 64
      tfile(12) = 65
      tfile(15) = 65 
c
c  auxilary arrays set up - read in nasa mo integrals
c
_IF1(ct)cdir$ vector
      return
      end
      subroutine newlis(flov,nbfo,noo,nbf,no,cororb,virorb,imap,
     1      eigval,wtemp,oorbsy,orbsym,nirred)
      implicit integer(a-z)
      REAL eigval,wtemp
      dimension flov(nirred,4),cororb(nirred),virorb(nirred),
     1   imap(nbfo),eigval(nbfo),wtemp(nbfo),oorbsy(nbfo),orbsym(nbf)
c
      nbf = nbfo
      no = noo
      do 10 irr = 1,nirred
      nbf = nbf - cororb(irr) - virorb(irr)
      no = no - cororb(irr)
10    continue
      write(6,*) ' nbf = ',nbf
      write(6,*) ' no = ',no
      ncore = noo - no
      if(nbf.eq.nbfo) then
      do 20 i=1,nbf
      orbsym(i)=oorbsy(i)
      imap(i)=i
   20 continue
      return
      end if
c
      do 40 irr = 1,nirred
      write(6,*) ' irr,flov ',irr,(flov(irr,jki),jki=1,4)
40    continue
c
      call tit_izo(imap,nbfo)
      coret = 0
      virt = 0
      itoto = 0
      itotv = no
      do 22 irr = 1,nirred
      if(flov(irr,1).eq.0) go to 23
      nosym = flov(irr,2) - flov(irr,1) + 1 - cororb(irr)
c
      coret = coret + cororb(irr)
      do 25 icnt = flov(irr,1)+cororb(irr),flov(irr,2)
      imap(icnt) = icnt - coret
25    continue
c
      flov(irr,1) = itoto + 1
      flov(irr,2) = flov(irr,1) + nosym - 1
c
      itoto = itoto + nosym
23    if(flov(irr,3).eq.0) go to 22
c
      nvsym = flov(irr,4) - flov(irr,3) + 1 -virorb(irr)
c
      do 27 icnt = flov(irr,3),flov(irr,4)-virorb(irr)
      imap(icnt) = icnt - virt - ncore
27    continue
      virt = virt + virorb(irr)
c
      flov(irr,3) = itotv + 1
      flov(irr,4) = flov(irr,3) + nvsym - 1
c
      itotv = itotv + nvsym
c
22    continue
c
      write(6,*) ' imap ',imap
c
      do 30 irr = 1,nirred
      write(6,*) ' irr,flov ',irr,(flov(irr,jki),jki=1,4)
30    continue
c
      do 50 i=1,nbfo
      if(imap(i).ne.0) then
      orbsym(imap(i)) = oorbsy(i)
      wtemp(imap(i)) = eigval(i)
      end if
50    continue
      call vclr(eigval,1,nbfo)
      do 60 i=1,nbf
60    eigval(i)=wtemp(i)
c
      return
      end
      subroutine prcss(buf2,bkt,ibkt,ni,nj,nk,nl,icnt,ibktsp,lnbuf,
     1                 idu,bkof1,bkof2,tfile,nbf,lenbf,maxval)
      implicit integer (a-z)
      REAL buf2(lnbuf),bkt(ibktsp)
      integer ibkt(*),ni(lnbuf),nj(lnbuf),nk(lnbuf),nl(lnbuf)
      integer idu(nbf),bkof1(16),bkof2(16),tfile(16)
c
c  first pack the doc,uoc flag for k and l into nj - which already
c  contains the doc,uoc flag for i and j - also, add 1
c
_IF1(ct)cdir$ block
      do 10 i = 1,icnt
_IF(cray)
      nj(i)=or(idu(nl(i)),shiftl(or(idu(nk(i)),shiftl(nj(i),1)),1))+1
_ELSE
      nj(i)=IOR32(idu(nl(i)),ishft(IOR32(idu(nk(i)),ishft(nj(i),1)),1))+1
_ENDIF
 10   continue
c
c  pack k and l into nk
c
      do 20 i = 1,icnt
_IF(cray)
      nk(i) = or(nl(i),shiftl(nk(i),8))
_ELSE
      nk(i) = IOR32(nl(i),ishft(nk(i),8))
_ENDIF
 20   continue
c
c  pack ijkl into nl
c
      do 30 i =1,icnt
_IF(cray)
      nl(i) = or(min(ni(i),nk(i)),shiftl(max(ni(i),nk(i)),16))
_ELSE
      nl(i) = IOR32(min(ni(i),nk(i)),ishft(max(ni(i),nk(i)),16))
_ENDIF
 30   continue
c
c  determine bkt offset for real address => ni
c  determine bkt offset for label address => nk
c
      do 40 i = 1,icnt
      ni(i) = bkof1(nj(i)) 
      nk(i) = bkof2(nj(i))
 40   continue
c
c  determine file to be written to if buffer full => nj
c
      do 50 i = 1,icnt
      nj(i) = tfile(nj(i))
 50   continue
c
c  put values in bkts and write out if necessary
c
      do 60 i = 1,icnt
      loff = nk(i)
      roff = ni(i)
      num = ibkt(loff+2) + 1
c
      if(num.gt.maxval) then
      ibkt(loff+1) = 0
      call swrit(nj(i),ibkt(loff+1),intowp(lenbf))
      num = 1
      end if
c
      bkt(roff+num) = buf2(i)
      ibkt(loff+2+num) = nl(i)
      ibkt(loff+2) = num
 60   continue
c
c  buf2, ni etc.... processed!
c
_IF1(ct)cdir$ vector
      return
      end
      subroutine rdnai(buf,ibuf,buf2,ibuf2,bkt,ibkt,ni,nj,nk,nl,h,
     2           f,doc,uoc,ptocc,idu,ktri,ltri,ksq,lsq,ptsq,bkof1,
     5           bkof2,tfile,jout,nirred,nbf,no,nv,ntri,
     6           ntris,nsqs,lnbuf,lenbf,ibktsp,maxval,enuc,orbsym,
     &           itap30,itap60,itap61,itap62,itap63,itap64,
     &           itap65,itap66,itap67,itap69,itap78,itap91,itap57)
      implicit integer (a-z)
      REAL f(ntris),h(ntri),buf(lnbuf),buf2(lnbuf),bkt(ibktsp)
      REAL enuc
      integer orbsym(nbf)
      integer ni(lnbuf),nj(lnbuf),nk(lnbuf),nl(lnbuf),doc(nirred)
      integer uoc(nirred),ptocc(nbf),idu(nbf),ktri(ntris),ltri(ntris)
      integer ksq(nsqs),lsq(nsqs),ptsq(nirred,nirred),bkof1(16)
      integer bkof2(16),tfile(16),ibuf(*),ibuf2(*),ibkt(*)
      integer florb(8,2)
INCLUDE(common/titan)
c
_IF1()      ibtsym(i,j) = IXOR32(i-1,j-1)+1
      icnt=0
      do isym=1,nirred
      florb(isym,1)=icnt+1
      icnt=icnt+doc(isym)+uoc(isym)
      florb(isym,2)=icnt
      enddo
c      print *,'orbsym ',orbsym
c
      call g_rd1el(ptocc,orbsym,h,nbf)
c      do i=1,no+nv
c      do j=1,i
c      print *,'i,j,h(ij) ',i,j,h(i*(i-1)/2+j)
c      enddo
c      enddo
c
c  one electron ints in h - now read in and sort two electron ints
c  first initialize the ibkt array
c
      do 11 i = 1,6
      ind1 = (i-1)*intowp(lenbf) + 1
      ibkt(ind1) = 0
      ibkt(ind1+1) = 0
 11   continue
      call g_rd2el(ptocc,orbsym,buf,bkt,ibkt,ni,nj,nk,nl,
     &             ibktsp,lnbuf,idu,bkof1,bkof2,tfile,nbf,
     &             lenbf,maxval)
C
      ind = 1
      ibkt(ind) = -1
      call swrit(itap60,ibkt(ind),intowp(lenbf))
      ind = ind + intowp(lenbf)
      ibkt(ind) = -1
      call swrit(itap61,ibkt(ind),intowp(lenbf))
      ind = ind + intowp(lenbf)
      ibkt(ind) = -1
      call swrit(itap62,ibkt(ind),intowp(lenbf))
      ind = ind + intowp(lenbf)
      ibkt(ind) = -1
      call swrit(itap63,ibkt(ind),intowp(lenbf))
      ind = ind + intowp(lenbf)
      ibkt(ind) = -1
      call swrit(itap64,ibkt(ind),intowp(lenbf))
      ind = ind + intowp(lenbf)
      ibkt(ind) = -1
      call swrit(itap65,ibkt(ind),intowp(lenbf))
      ibkt(1) = -1
      ibkt(2) = 0
      call swrit(itap66,ibkt,intowp(lenbf))
C
      return
      end
      subroutine read30(itap30,mpoint,mconst,mcalcs,ncalcs,nbf,
     .                  nsymhf,mxcoef,nlamda,itemp,jout,nc,no,ityp,
     .                  nirred,flov,orbsym,eigval,ptocc,wtemp)
c
      implicit REAL (a-h,o-z)
c
cray  character*8 ityp(nsymhf),icsym
      character*4 ityp(nsymhf),icsym
      integer nlamda(nsymhf),itemp(mcalcs),nc(nsymhf),
     .        flov(nirred,4),orbsym(nbf),symorb,ptocc(nbf)
INCLUDE(common/prnprn)
      dimension eigval(nbf),wtemp(nbf)
c
      jpoint = 101 + mconst + mpoint
c
c
      call wreadw (itap30,itemp,mcalcs,jpoint,jpoint)
      loccal = itemp(ncalcs)
      jpoint = loccal + 60
      call wreadw (itap30,locvec,1,jpoint,jpoint)
      locvec=locvec+intowp(mxcoef)
      call wreadw(itap30,eigval,intowp(nbf),locvec,locvec)
      call wreadw (itap30,ityp,nsymhf,locvec,locvec)
c     write(*,*) ' locvec,mxcoef = ',locvec,mxcoef
      call wreadw (itap30,nlamda,nsymhf,locvec,locvec)
      write(*,*) '  nlamda ',nlamda
      call wreadw (itap30,nc,nsymhf,locvec,locvec)
      no=0
      do 555 i=1,nsymhf
         no=no+nc(i)
         write (jout,*) ' nc=',nc(i)
  555 continue
c
      do 805 isym = 1,nirred
      flov(isym,1) = 0
      flov(isym,2) = -1
      flov(isym,3) = 0
      flov(isym,4) = -1
  805 continue
c
      icnt=0
      do 705 isym=1,nsymhf
      icsym=ityp(isym)
      icnth=icnt+1
      noi=nc(isym)
      do 700 i=1,noi
      if(odebug(33)) write(6,*) ' icsym,icnt ',icnt,icsym
      icnt=icnt+1
      if(icsym.eq.'a1  ')orbsym(icnt)=0
      if(icsym.eq.'a2  ')orbsym(icnt)=1
      if(icsym.eq.'b1  ')orbsym(icnt)=2
      if(icsym.eq.'b2  ')orbsym(icnt)=3
      if(icsym.eq.'a   ')orbsym(icnt)=0
      if(icsym.eq.'b   ')orbsym(icnt)=1
      if(icsym.eq.'a''  ')orbsym(icnt)=0
      if(icsym.eq.'a"  ')orbsym(icnt)=1
      if(icsym.eq.'ag  ')orbsym(icnt)=0
      if(icsym.eq.'b1g ')orbsym(icnt)=1
      if(icsym.eq.'b2g ')orbsym(icnt)=2
      if(icsym.eq.'b3g ')orbsym(icnt)=3
      if(icsym.eq.'au  ')orbsym(icnt)=4
      if(icsym.eq.'b1u ')orbsym(icnt)=5
      if(icsym.eq.'b2u ')orbsym(icnt)=6
      if(icsym.eq.'b3u ')orbsym(icnt)=7
      if(nirred.gt.4)go to 690
      if(icsym.eq.'ag  ')orbsym(icnt)=0
      if(icsym.eq.'bg  ')orbsym(icnt)=1
      if(icsym.eq.'au  ')orbsym(icnt)=2
      if(icsym.eq.'bu  ')orbsym(icnt)=3
  690 continue
  700 continue
      if(icnth-1.ne.icnt) then
      symorb = orbsym(icnt) + 1
      flov(symorb,1)=icnth
      flov(symorb,2)=icnt
      end if
c      write(*,*)'isym ',isym,'f1 ',flov(isym,1)
c      write(*,*)'isym ',isym,'f2 ',flov(isym,2)
  705 continue
      do 715 isym=1,nsymhf
      icsym=ityp(isym)
      icnth=icnt+1
      nvi=nlamda(isym)-nc(isym)
      do 710 i=1,nvi
      icnt=icnt+1
      if(icsym.eq.'a1  ')orbsym(icnt)=0
      if(icsym.eq.'a2  ')orbsym(icnt)=1
      if(icsym.eq.'b1  ')orbsym(icnt)=2
      if(icsym.eq.'b2  ')orbsym(icnt)=3
      if(icsym.eq.'a   ')orbsym(icnt)=0
      if(icsym.eq.'b   ')orbsym(icnt)=1
      if(icsym.eq.'a''  ')orbsym(icnt)=0
      if(icsym.eq.'a"  ')orbsym(icnt)=1
      if(icsym.eq.'ag  ')orbsym(icnt)=0
      if(icsym.eq.'b1g ')orbsym(icnt)=1
      if(icsym.eq.'b2g ')orbsym(icnt)=2
      if(icsym.eq.'b3g ')orbsym(icnt)=3
      if(icsym.eq.'au  ')orbsym(icnt)=4
      if(icsym.eq.'b1u ')orbsym(icnt)=5
      if(icsym.eq.'b2u ')orbsym(icnt)=6
      if(icsym.eq.'b3u ')orbsym(icnt)=7
      if(nirred.gt.4)go to 691
      if(icsym.eq.'ag  ')orbsym(icnt)=0
      if(icsym.eq.'bg  ')orbsym(icnt)=1
      if(icsym.eq.'au  ')orbsym(icnt)=2
      if(icsym.eq.'bu  ')orbsym(icnt)=3
  691 continue
  710 continue
      if(icnth-1.ne.icnt) then
      symorb = orbsym(icnt) + 1
      flov(symorb,3)=icnth
      flov(symorb,4)=icnt
      end if
c       write(*,*) 'isym ',isym,' f3 ',flov(isym,3)
c       write(*,*) 'isym ',isym,' f4 ',flov(isym,4)
  715 continue
      do 4041 jki=1,nirred
      if(odebug(33)) then
        write(6,*) ' jki,flov(jki,1) ',jki,flov(jki,1)
        write(6,*) ' jki,flov(jki,2) ',jki,flov(jki,2)
        write(6,*) ' jki,flov(jki,3) ',jki,flov(jki,3)
      endif
4041  write(6,*) ' jki,flov(jki,4) ',jki,flov(jki,4)
      icnt=0
      iof = 0
      do 558 i=1,nsymhf
      if(i.ne.1) iof = iof + nlamda(i-1)
      noi=nc(i)
      do 557 j=1,noi
      icnt=icnt+1
      ipt=iof+j
c     write(6,*)'ipt=',ipt,icnt
  557 ptocc(ipt)=icnt
  558 continue
c     write(*,*) '  virtuals'
      iof = 0
      do 658 i=1,nsymhf
      if(i.ne.1) iof = iof + nlamda(i-1)
      nvi=nlamda(i)-nc(i)
      do 657 j=1,nvi
      icnt=icnt+1
      ipt=iof+nc(i) + j
c     write(6,*) iof,nc(i),j
c     write(6,*)'ipt=',ipt,icnt
  657 ptocc(ipt)=icnt
  658 continue
      do 559 i=1,nbf
      ipt=ptocc(i)
  559 wtemp(ipt)=eigval(i)
      do 569 i=1,nbf
      eigval(i)=wtemp(i)
      if(odebug(33)) write(6,*)'i=',i,'e(i)=',eigval(i)
  569 continue
      write(*,*) 'orbsym ',orbsym
      return
      end
      subroutine rwtp(unit,addr)
      implicit integer (a-z)
      addr=0
      return
      end
      integer FUNCTION SEC2I(N)
_IF(cray)
      SEC2I=1024*N
_ELSE
      SEC2I=512*N
_ENDIF
      RETURN
      END
      subroutine setd3(flov,norbt,orbsym,nprb,nprsq,iofsy,indx,indxl,
     1   ioff,iofpr,fpbkt,prbkt,nv,ntriv,nbf,nirr,maxbkt,space1,
     1    nbkt,nprgs,npr,norbv,ntrio,apbkt)
      implicit integer(a-z)
INCLUDE(common/prnprn)
      dimension flov(nirr,4),indx(ntrio),indxl(nv*nv),nprb(nirr)
     1 , nprsq(nirr),iofsy(nirr+1),orbsym(nbf),norbt(nirr),norbv(nirr)
      dimension npr(nirr,3,2),apbkt(maxbkt)
      dimension nprgs(nirr,nirr)
      dimension ioff(nbf),iofpr(nirr+1),fpbkt(maxbkt),prbkt(2)
c
c       flov(nirr,4)  -  1      first occupied in symmetry
c                       2      last occupied in symmetry
c                       3      first virtual in symmetry
c                       4      last virtual in symmetry
c          norbt(nirr)    no of orbitals in each symmetry type
c          orbsym(nbf)     symmetry type of each orbital (0-nirr-1)
c          npr(nirr)      no of triangle pairs of each sym type(virt)
c          nprsq(nirr)     no of square pairs of each sym type (virt)
c         nprgs(irr,jrr)  gives no of orbitals which with orbital of
c                         symtyp irr give a pair of symm jrr
c         nprgsq(irr,jrr)  equivalent for squares
c          iofsy(nirr+1)   offset to beginning of each block(ij|kl)
c          iofpr(nirr) offset to beginning of each pair (ik)
c          indx(ij) index of ij triangle in symmetry
c          indxl(ij) index of ij square in symmetry
c
      if(odebug(33)) then
       write(6,*) ' orbsym ',orbsym
       write(6,*) ' flov ',flov
      endif
      no=nbf-nv
      do 5 i=1,nirr
      norbv(i)=flov(i,4)-flov(i,3)+1
 5    norbt(i)=flov(i,2)-flov(i,1)+1
      if(odebug(33)) write(6,*) ' norbt ',norbt
      call tit_izo(nprb,nirr)
      call tit_izo(nprgs,nirr*nirr)
      do 10 i=1,no
      do 20 j=1,i
      totsym=IXOR32(orbsym(i),orbsym(j))
      nprgs(orbsym(i)+1,totsym+1)=
     1      nprgs(orbsym(i)+1,totsym+1)+1
      nprb(totsym+1)=nprb(totsym+1)+1
20    continue
10    continue
      if(odebug(33)) then
       write(6,*) ' nprb ',nprb
       write(6,*) ' nprgs ',nprgs
      endif
      call tit_izo(nprsq,nirr)
      do 12 i=1,nv
      do 22 j=1,nv
      totsym=IXOR32(orbsym(i+no),orbsym(j+no))
      nprsq(totsym+1)=nprsq(totsym+1)+1
22    continue
12     continue
      if(odebug(33)) write(6,*) ' nprsq ',nprsq
      call tit_izo(iofsy,nirr+1)
      do 45 ir=1,nirr
45    iofsy(ir+1)=nprb(ir)*nprsq(ir)+iofsy(ir)
      call tit_izo(iofpr,nirr+1)
      do 46 ir=1,nirr
46    iofpr(ir+1)=iofpr(ir)+nprb(ir)
      if(odebug(33)) write(6,*) ' iofsy ',iofsy
c
      iofst2=0
      if(odebug(33)) write(6,*) ' nirr ',nirr
      do 30 k=1,nirr
      if(k.gt.1) iofst2=iofst2+nprb(k-1)
      iofst1=0
      do 40 l=1,nirr
      if(l.gt.1) iofst1=iofst1+nprgs(l-1,k)
      if(odebug(33)) write(6,*) ' l,norbt(l) ',l,norbt(l)
      do 50 i=1,norbt(l)
      nsym=IXOR32(l-1,k-1)+1
      if(nsym.gt.l) go to 50
      endl=norbt(nsym)
      if(k.eq.1) endl=i
      if(odebug(33)) write(6,*) 'k,endl ',k,endl
      do 60 j=1,endl
      ixj=i+flov(l,1)-1
      jxj=j+flov(nsym,1)-1
      ij=ioff(ixj)+jxj
      if(odebug(33)) then
       write(6,*) ' k,l,nsym ',k,l,nsym
       write(6,*) ' i,j,ixj,jxj,ij ',i,j,ixj,jxj,ij
      endif
      if(k.eq.1) indx(ij)=ioff(i)+j+iofst1
      if(k.gt.1) indx(ij)=(i-1)*norbt(nsym)+j+iofst1
      if(odebug(33)) write(6,*) ' indx ,ij ',ij,indx(ij)
60    continue
50    continue
40    continue
30    continue
      if(odebug(33)) then
       write(6,*) ' indx ',indx
       write(6,*) ' indxl ',indxl
       write(6,*) ' nprb ',nprb
       write(6,*) ' nprsq ',nprsq
       write(6,*) ' indx ',indx
       write(6,*) ' indxl ',indxl
      endif
c
      len=0
      nbkt=1
      fpbkt(1)=0
      if(odebug(33)) write(6,*) ' space1 ',space1
      call tit_izo(apbkt,maxbkt)
      do 401 irr=1,nirr
      do 402 ik=1,nprb(irr)
      if(odebug(33)) write(6,*) ' irr,ik ',irr,ik,iofpr(irr)
      len=len+nprsq(irr)
      if(odebug(33)) write(6,*) ' len,space1 ',len,space1
      if(len.gt.space1) then
         nbkt=nbkt+1
         len=len-nprsq(irr)
         apbkt(nbkt)=len+apbkt(nbkt-1)
         len=nprsq(irr)
         fpbkt(nbkt)=ik-1+iofpr(irr)
      end if
      prbkt(ik+iofpr(irr))=nbkt
402   continue
401   continue
      fpbkt(nbkt+1)=iofpr(nirr+1)
      if(odebug(33)) then
       write(6,*) ' nbkt at end of setup ',nbkt
       write(6,*) ' fpbkt ',fpbkt
       write(6,*) ' prbkt ',prbkt
       write(6,*) ' apbkt ',apbkt
      endif
      do 109 irr=1,nirr
      npr(irr,1,1)=nprb(irr)
109   continue
      return
      end
      subroutine sete1(flov,norbo,orbsym,norsq,iofosy,indxol,
     1iofosq,fsobkt,psobkt,no,nbf,nirr,maxbkt,space1,
     1    nbkt,norgsq,npr,nv,norbv,indxl,nvrsq,apbkt)
      implicit integer(a-z)
INCLUDE(common/prnprn)
      dimension flov(nirr,4),indxol(no*nv),norbv(nirr)
     1 , norsq(nirr),iofosy(nirr+1),orbsym(nbf),norbo(nirr)
      dimension npr(nirr,3,2),indxl(no*no),nvrsq(nirr)
      dimension norgsq(nirr,nirr),apbkt(maxbkt)
      dimension iofosq(nirr+1),fsobkt(maxbkt),psobkt(nv*no)
c
c       flov(nirr,4)  -  1      first occupied in symmetry
c                       2      last occupied in symmetry
c                       3      first virtual in symmetry
c                       4      last virtual in symmetry
c          norbo(nirr)    no of occ orbitals in each symmetry type
c          orbsym(nbf)     symmetry type of each orbital (0-nirr-1)
c          norsq(nirr)     no of square pairs of each sym type (occ)
c         norgsq(irr,jrr)  equivalent for squares
c          iofosy(nirr+1)   offset to beginning of each block(ij|kl)
c          iofosq(nirr) offset to beginning of each square (ik)
c          indx(ij) index of ij triangle in symmetry
c          indxl(ij) index of ij square in symmetry
c
      if(odebug(33))  then
       write(6,*) ' orbsym ',orbsym
       write(6,*) ' nbf,nirr,no in setupo ',nbf,nirr,no
       write(6,*) ' flov ',flov
      endif
      do 5 i=1,nirr
      norbv(i)=flov(i,4)-flov(i,3)+1
 5    norbo(i)=flov(i,2)-flov(i,1)+1
      if(odebug(33)) write(6,*) ' norbo ',norbo
      call tit_izo(norsq,nirr)
      call tit_izo(norgsq,nirr*nirr)
      do 12 i=1,no
      do 22 j=1,nv
      totsym=IXOR32(orbsym(i),orbsym(j+no))
      norgsq(orbsym(i)+1,totsym+1)=
     1 norgsq(orbsym(i)+1,totsym+1)+1
      norsq(totsym+1)=norsq(totsym+1)+1
22    continue
12     continue
      call tit_izo(nvrsq,nirr)
      do 112 i=1,no
      do 122 j=1,no
      totsym=IXOR32(orbsym(i),orbsym(j))
      nvrsq(totsym+1)=nvrsq(totsym+1)+1
122   continue
112   continue
      if(odebug(33)) write(6,*) ' nvrsq ',nvrsq
      call tit_izo(iofosy,nirr+1)
      do 45 ir=1,nirr
45    iofosy(ir+1)=norsq(ir)*nvrsq(ir)+iofosy(ir)
      call tit_izo(iofosq,nirr+1)
      do 46 ir=1,nirr
46    iofosq(ir+1)=iofosq(ir)+norsq(ir)
      if(odebug(33)) write(6,*) ' iofosy ',iofosy
c
      iofst2=0
      do 300 k=1,nirr
      if(k.gt.1) iofst2=iofst2+norsq(k-1)
      iofst1=0
      do 400 l=1,nirr
      if(l.gt.1) iofst1=iofst1+norgsq(l-1,k)
      do 500 i=1,norbo(l)
      nsym=IXOR32(l-1,k-1)+1
      endl=norbv(nsym)
      do 600 j=1,endl
      ij=(i+flov(l,1)-2)*nv+flov(nsym,3)+j-no-1
      indxol(ij)=(i-1)*norbv(nsym)+j+iofst1
600   continue
500   continue
400   continue
300   continue
      if(odebug(33)) then
       write(6,*) ' indxol ',indxol
       write(6,*) ' norsq ',norsq
      endif
c
      len=0
      nbkt=1
      bgbkt=0
      fsobkt(1)=0
      call tit_izo(fsobkt,maxbkt)
      call tit_izo(apbkt,maxbkt)
      if(odebug(33)) write(6,*) ' space1 ',space1
      do 401 irr=1,nirr
      do 402 ik=1,norsq(irr)
      if(odebug(33)) write(6,*) ' irr,ik ',irr,ik,iofosq(irr)
      len=len+nvrsq(irr)
      if(len.gt.space1) then
         nbkt=nbkt+1
         len=len-nvrsq(irr)
         apbkt(nbkt)=len+apbkt(nbkt-1)
         len=nvrsq(irr)
         fsobkt(nbkt)=ik-1+iofosq(irr)
      end if
      psobkt(ik+iofosq(irr))=nbkt
402   continue
401   continue
      fsobkt(nbkt+1)=iofosq(nirr+1)
      if(odebug(33)) then
       write(6,*) ' nbkt ',nbkt
       write(6,*) ' fsobkt ',fsobkt
       write(6,*) ' psobkt ',psobkt
       write(6,*) ' apbkt ',(apbkt(jki),jki=1,nbkt+1)
      endif
      do 109 i=1,nirr
      npr(i,3,1)=norsq(i)
109   continue
      return
      end
      subroutine setf(orbsym,norsq,iofosy,
     1iofosq,fsobkt,psobkt,no,nbf,nirr,maxbkt,space1,
     1    nbkt,nv,nvrsq,apbkt)
      implicit integer(a-z)
INCLUDE(common/prnprn)
      dimension norsq(nirr),iofosy(nirr+1),orbsym(nbf)
      dimension nvrsq(nirr)
      dimension apbkt(maxbkt)
      dimension iofosq(nirr+1),fsobkt(maxbkt),psobkt(no*nv)
c
c          norbo(nirr)    no of occ orbitals in each symmetry type
c          orbsym(nbf)     symmetry type of each orbital (0-nirr-1)
c          norsq(nirr)     no of square pairs of each sym type (occ)
c         norgsq(irr,jrr)  equivalent for squares
c          iofosy(nirr+1)   offset to beginning of each block(ij|kl)
c          iofosq(nirr) offset to beginning of each square (ik)
c          indx(ij) index of ij triangle in symmetry
c          indxl(ij) index of ij square in symmetry
c
      if(odebug(33)) then
       write(6,*) ' orbsym ',orbsym
       write(6,*) ' nbf,nirr,no in setupo ',nbf,nirr,no
      endif
      call tit_izo(norsq,nirr)
      do 12 i=1,nv
      do 22 j=1,no
      totsym=IXOR32(orbsym(i+no),orbsym(j))
      norsq(totsym+1)=norsq(totsym+1)+1
22    continue
12     continue
      call tit_izo(nvrsq,nirr)
      do 112 i=1,nv
      do 122 j=1,nv
      totsym=IXOR32(orbsym(i+no),orbsym(j+no))
      nvrsq(totsym+1)=nvrsq(totsym+1)+1
122   continue
112   continue
      if(odebug(33)) write(6,*) ' nvrsq ',nvrsq
      call tit_izo(iofosy,nirr+1)
      do 45 ir=1,nirr
45    iofosy(ir+1)=norsq(ir)*nvrsq(ir)+iofosy(ir)
      call tit_izo(iofosq,nirr+1)
      do 46 ir=1,nirr
46    iofosq(ir+1)=iofosq(ir)+norsq(ir)
      if(odebug(33)) then
       write(6,*) ' iofosy ',iofosy
       write(6,*) ' iofosq ',iofosq
       write(6,*) ' norsq ',norsq
      endif
c
      len=0
      nbkt=1
      bgbkt=0
      fsobkt(1)=0
      call tit_izo(fsobkt,maxbkt)
      call tit_izo(apbkt,maxbkt)
      if(odebug(33)) write(6,*) ' space1 ',space1
      do 401 irr=1,nirr
      do 402 ik=1,norsq(irr)
      if(odebug(33)) write(6,*) ' irr,ik ',irr,ik,iofosq(irr),nvrsq(irr)
      len=len+nvrsq(irr)
      if(odebug(33)) write(6,*) ' len in setf ',len
      if(len.gt.space1) then
         nbkt=nbkt+1
         len=len-nvrsq(irr)
         apbkt(nbkt)=len+apbkt(nbkt-1)
         len=nvrsq(irr)
         if(odebug(33)) 
     +   write(6,*) ' updating fsobkt ',nbkt,fsobkt(nbkt),ik-1
         fsobkt(nbkt)=ik-1+iofosq(irr)
      end if
      psobkt(ik+iofosq(irr))=nbkt
402   continue
401   continue
      fsobkt(nbkt+1)=iofosq(nirr+1)
      if(odebug(33)) then
       write(6,*) ' nbkt ',nbkt
       write(6,*) ' fsobkt ',fsobkt
       write(6,*) ' psobkt ',psobkt
       write(6,*) ' apbkt ',(apbkt(jki),jki=1,nbkt+1)
      endif
      return
      end
_IF(hpux11)
c HP compiler bug JAGae51337
c$HP$ OPTIMIZE LEVEL2
c$HP$ OPTIMIZE ASSUME_NO_PARAMETERS_OVERLAPS OFF
_ENDIF
      subroutine setmbf(mbuf,ibuf)
      implicit REAL (a-h,o-z)
      mbuf=ibuf
      return
      end
_IF(hpux11)
c$HP$ OPTIMIZE ASSUME_NO_PARAMETERS_OVERLAPS ON
_ENDIF
      subroutine sett2(flov,norbv,norbo,orbsym,nprb,nprsq,iofsy,indxl,
     1   ioff,iofpr,nv,ntriv,nbf,nirr,space1,indxvl,
     1    nprgs,nprgsq,ndimt1,ndimt2,npr,no)
      implicit integer(a-z)
INCLUDE(common/prnprn)
      dimension flov(nirr,4),indxl(no*no),indxvl(ntriv),nprb(nirr)
     1 , nprsq(nirr),iofsy(nirr+1),orbsym(nbf),norbv(nirr)
      dimension nprgs(nirr,nirr),nprgsq(nirr,nirr),norbo(nirr)
      dimension ioff(nbf),iofpr(nirr+1),npr(nirr,3,2)
c
c
      if(odebug(33)) then
       write(6,*) ' orbsym ',orbsym
       write(6,*) ' flov ',flov
      endif
      do 3 i=1,nbf
3      ioff(i)=i*(i-1)/2
      do 5 i=1,nirr
      norbo(i)=flov(i,2)-flov(i,1)+1
 5    norbv(i)=flov(i,4)-flov(i,3)+1
      if(odebug(33)) write(6,*) ' norbv ',norbv
      ndimt1=0
      do 109 i=1,nirr
109   ndimt1=ndimt1+norbo(i)*norbv(i)
      if(odebug(33)) write(6,*) ' ndimt1 ',ndimt1
      call tit_izo(nprb,nirr)
      call tit_izo(nprgs,nirr*nirr)
      do 10 i=1,no
      do 20 j=1,no
      totsym=IXOR32(orbsym(i),orbsym(j))
      nprgs(orbsym(i)+1,totsym+1)=
     1      nprgs(orbsym(i)+1,totsym+1)+1
      nprb(totsym+1)=nprb(totsym+1)+1
20    continue
10    continue
      if(odebug(33)) then
       write(6,*) ' nprb ',nprb
       write(6,*) ' nprgs ',nprgs
      endif
      call tit_izo(nprsq,nirr)
      call tit_izo(nprgsq,nirr*nirr)
      do 12 i=1,nv
      do 22 j=1,i
      totsym=IXOR32(orbsym(i+no),orbsym(j+no))
      nprgsq(orbsym(i+no)+1,totsym+1)=
     1 nprgsq(orbsym(i+no)+1,totsym+1)+1
      nprsq(totsym+1)=nprsq(totsym+1)+1
22    continue
12     continue
      if(odebug(33)) then
       write(6,*) ' nprsq ',nprsq
       write(6,*) ' nprgsq ',nprgsq
      endif
      call tit_izo(iofsy,nirr+1)
      do 45 ir=1,nirr
45    iofsy(ir+1)=nprb(ir)*nprsq(ir)+iofsy(ir)
      if(space1.lt.iofsy(nirr+1)) then
          write(6,*) ' not enough space to make t2 '
          write(6,*) ' need ',iofsy(nirr+1),' words : have ',
     1        space1,' words '
           call mabort
      else
          ndimt2=iofsy(nirr+1)
      end if
      if(odebug(33)) then
       write(6,*) ' iofsy ',iofsy
       write(6,*) ' ndimt2 = ',ndimt2
      endif
      call tit_izo(iofpr,nirr+1)
      do 46 ir=1,nirr
46    iofpr(ir+1)=iofpr(ir)+nprb(ir)
      if(odebug(33)) write(6,*) ' iofsy ',iofsy
c
      iofst2=0
      if(odebug(33)) write(6,*) ' nirr ',nirr
      do 30 k=1,nirr
_IF1()cj      if(k.gt.1) iofst2=iofst2+nprb(k-1)
      iofst1=0
      do 40 l=1,nirr
      if(l.gt.1) iofst1=iofst1+nprgs(l-1,k)
      do 50 i=1,norbo(l)
      nsym=IXOR32(l-1,k-1)+1
      endl=norbo(nsym)
      do 60 j=1,endl
      ij=(i+flov(l,1)-2)*no+flov(nsym,1)+j-1
      indxl(ij)=(i-1)*norbo(nsym)+j+iofst1
      if(odebug(33)) write(6,*) ' indxl ,ij ',ij,indxl(ij)
60    continue
50    continue
40    continue
30    continue
      iofst2=0
      do 300 k=1,nirr
      if(k.gt.1) iofst2=iofst2+nprsq(k-1)
      iofst1=0
      do 400 l=1,nirr
      if(l.gt.1) iofst1=iofst1+nprgsq(l-1,k)
      do 500 i=1,norbv(l)
      nsym=IXOR32(l-1,k-1)+1
      if(nsym.gt.l) go to 500
      endl=norbv(nsym)
      if(k.eq.1) endl=i
      do 600 j=1,endl
      ixj=i+flov(l,3)-no-1
      jkj=j+flov(nsym,3)-no-1
      ij=ioff(ixj)+jkj
      if(k.eq.1) indxvl(ij)=ioff(i)+j+iofst1
      if(k.gt.1) indxvl(ij)=(i-1)*norbv(nsym)+j+iofst1
600   continue
500   continue
400   continue
300   continue
      if(odebug(33)) then
       write(6,*) ' indxvl ',indxvl
       write(6,*) ' nprb ',nprb
       write(6,*) ' nprsq ',nprsq
       write(6,*) ' indxl ',indxl
      endif
      do 101 ji=1,nirr
      npr(ji,2,1)=nprsq(ji)
101   npr(ji,1,2)=nprb(ji)
c
      return
      end
      subroutine setup(flov,norbt,orbsym,nprb,nprsq,iofsy,indx,indxl,
     1   ioff,iofpr,fpbkt,prbkt,nv,ntriv,nbf,nirr,maxbkt,space1,
     1    nbkt,nprgs,nprgsq,npr,no,apbkt)
      implicit integer(a-z)
INCLUDE(common/prnprn)
      dimension flov(nirr,4),indx(ntriv),indxl(nv*nv),nprb(nirr)
     1 , nprsq(nirr),iofsy(nirr+1),orbsym(nbf),norbt(nirr)
      dimension npr(nirr,3,2),apbkt(maxbkt)
      dimension nprgs(nirr,nirr),nprgsq(nirr,nirr)
      dimension ioff(nbf),iofpr(nirr+1),fpbkt(maxbkt),prbkt(ntriv)
c
c       flov(nirr,4)  -  1      first occupied in symmetry
c                       2      last occupied in symmetry
c                       3      first virtual in symmetry
c                       4      last virtual in symmetry
c          norbt(nirr)    no of orbitals in each symmetry type
c          orbsym(nbf)     symmetry type of each orbital (0-nirr-1)
c          npr(nirr)      no of triangle pairs of each sym type(virt)
c          nprsq(nirr)     no of square pairs of each sym type (virt)
c         nprgs(irr,jrr)  gives no of orbitals which with orbital of
c                         symtyp irr give a pair of symm jrr
c         nprgsq(irr,jrr)  equivalent for squares
c          iofsy(nirr+1)   offset to beginning of each block(ij|kl)
c          iofpr(nirr) offset to beginning of each pair (ik)
c          indx(ij) index of ij triangle in symmetry
c          indxl(ij) index of ij square in symmetry
c
c      write(6,*) ' orbsym ',orbsym
c      write(6,*) ' flov ',flov
      no=nbf-nv
      do 5 i=1,nirr
 5    norbt(i)=flov(i,4)-flov(i,3)+1
      if(odebug(33)) write(6,*) ' norbt ',norbt
      call tit_izo(nprb,nirr)
      call tit_izo(nprgs,nirr*nirr)
      do 10 i=1,nv
      do 20 j=1,i
      totsym=IXOR32(orbsym(i+no),orbsym(j+no))
      nprgs(orbsym(i+no)+1,totsym+1)=
     1      nprgs(orbsym(i+no)+1,totsym+1)+1
      nprb(totsym+1)=nprb(totsym+1)+1
20    continue
10    continue
      if(odebug(33)) then
       write(6,*) ' nprb ',nprb
       write(6,*) ' nprgs ',nprgs
      endif
      call tit_izo(nprsq,nirr)
      call tit_izo(nprgsq,nirr*nirr)
      do 12 i=1,nv
      do 22 j=1,nv
      totsym=IXOR32(orbsym(i+no),orbsym(j+no))
      nprgsq(orbsym(i+no)+1,totsym+1)=
     1 nprgsq(orbsym(i+no)+1,totsym+1)+1
      nprsq(totsym+1)=nprsq(totsym+1)+1
22    continue
12     continue
      if(odebug(33)) then
       write(6,*) ' nprsq ',nprsq
       write(6,*) ' nprgsq ',nprgsq
      endif
      call tit_izo(iofsy,nirr+1)
      do 45 ir=1,nirr
45    iofsy(ir+1)=nprb(ir)*nprsq(ir)+iofsy(ir)
      call tit_izo(iofpr,nirr+1)
      do 46 ir=1,nirr
46    iofpr(ir+1)=iofpr(ir)+nprb(ir)
      if(odebug(33)) write(6,*) ' iofsy ',iofsy
c
      iofst2=0
      if(odebug(33)) write(6,*) ' nirr ',nirr
      do 30 k=1,nirr
      if(k.gt.1) iofst2=iofst2+nprb(k-1)
      iofst1=0
      do 40 l=1,nirr
      if(l.gt.1) iofst1=iofst1+nprgs(l-1,k)
      if(odebug(33)) write(6,*) ' l,norbt(l) ',l,norbt(l)
      do 50 i=1,norbt(l)
      nsym=IXOR32(l-1,k-1)+1
      if(nsym.gt.l) go to 50
      endl=norbt(nsym)
      if(k.eq.1) endl=i
      if(odebug(33)) write(6,*) 'k,endl ',k,endl
      do 60 j=1,endl
      ixj=i+flov(l,3)-no-1
      jxj=j+flov(nsym,3)-no-1
      ij=ioff(ixj)+jxj
      if(odebug(33)) then
       write(6,*) ' k,l,nsym ',k,l,nsym
       write(6,*) ' i,j,ixj,jxj,ij ',i,j,ixj,jxj,ij
      endif
      if(k.eq.1) indx(ij)=ioff(i)+j+iofst1
      if(k.gt.1) indx(ij)=(i-1)*norbt(nsym)+j+iofst1
      if(odebug(33)) write(6,*) ' indx ,ij ',ij,indx(ij)
60    continue
50    continue
40    continue
30    continue
      if(odebug(33)) write(6,*) ' indx ',indx
      iofst2=0
      do 300 k=1,nirr
      if(k.gt.1) iofst2=iofst2+nprsq(k-1)
      iofst1=0
      do 400 l=1,nirr
      if(l.gt.1) iofst1=iofst1+nprgsq(l-1,k)
      do 500 i=1,norbt(l)
      nsym=IXOR32(l-1,k-1)+1
      endl=norbt(nsym)
      do 600 j=1,endl
      ij=(i+flov(l,3)-no-2)*nv+j+flov(nsym,3)-no-1
      indxl(ij)=(i-1)*norbt(nsym)+j+iofst1
600   continue
500   continue
400   continue
300   continue
      if(odebug(33)) then
       write(6,*) ' indxl ',indxl
       write(6,*) ' nprb ',nprb
       write(6,*) ' nprsq ',nprsq
       write(6,*) ' indx ',indx
       write(6,*) ' indxl ',indxl
      endif
c
      len=0
      nbkt=1
      bgbkt=0
      fpbkt(1)=0
      if(odebug(33)) write(6,*) ' space1 ',space1
      call tit_izo(fpbkt,maxbkt)
      call tit_izo(apbkt,maxbkt)
      do 401 irr=1,nirr
      do 402 ik=1,nprb(irr)
      if(odebug(33)) write(6,*) ' irr,ik ',irr,ik,iofpr(irr)
      len=len+nprsq(irr)
      if(odebug(33)) write(6,*) ' len,space1 ',len,space1
      if(len.gt.space1) then
         nbkt=nbkt+1
         if (nbkt.gt.maxbkt) call caserr('maxbkt overflow; cf tsort')
         len=len-nprsq(irr)
         apbkt(nbkt)=len+apbkt(nbkt-1)
         len=nprsq(irr)
         fpbkt(nbkt)=ik-1+iofpr(irr)
      end if
      prbkt(ik+iofpr(irr))=nbkt
402   continue
401   continue
      fpbkt(nbkt+1)=iofpr(nirr+1)
      if(odebug(33)) then
       write(6,*) ' nbkt at end of setup ',nbkt
       write(6,*) ' fpbkt ',fpbkt
       write(6,*) ' prbkt ',prbkt
       write(6,*) ' apbkt ',apbkt
      endif
      do 109 irr=1,nirr
      npr(irr,2,2)=nprsq(irr)
      npr(irr,2,1)=nprb(irr)
109   continue
      return
      end
      subroutine setupo(flov,norbo,orbsym,norsq,iofosy,indxol,
     1iofosq,fsobkt,psobkt,no,nbf,nirr,maxbkt,space1,
     1    nbkt,npr,apbkt)
      implicit integer(a-z)
INCLUDE(common/prnprn)
      dimension flov(nirr,4),indxol(no*no)
     1 , norsq(nirr),iofosy(nirr+1),orbsym(nbf),norbo(nirr)
      dimension npr(nirr,3,2),apbkt(maxbkt)
      dimension iofosq(nirr+1),fsobkt(maxbkt),psobkt(no*no)
c
c       flov(nirr,4)  -  1      first occupied in symmetry
c                       2      last occupied in symmetry
c                       3      first virtual in symmetry
c                       4      last virtual in symmetry
c          norbo(nirr)    no of occ orbitals in each symmetry type
c          orbsym(nbf)     symmetry type of each orbital (0-nirr-1)
c          norsq(nirr)     no of square pairs of each sym type (occ)
c         norgsq(irr,jrr)  equivalent for squares
c          iofosy(nirr+1)   offset to beginning of each block(ij|kl)
c          iofosq(nirr) offset to beginning of each square (ik)
c          indx(ij) index of ij triangle in symmetry
c          indxl(ij) index of ij square in symmetry
c
      if(odebug(33)) then
       write(6,*) ' orbsym ',orbsym
       write(6,*) ' nbf,nirr,no in setupo ',nbf,nirr,no
       write(6,*) ' flov ',flov
      endif
      do 5 i=1,nirr
 5    norbo(i)=flov(i,2)-flov(i,1)+1
      if(odebug(33)) write(6,*) ' norbo ',norbo
      call tit_izo(norsq,nirr)
      do 12 i=1,no
      do 22 j=1,no
      totsym=IXOR32(orbsym(i),orbsym(j))
      norsq(totsym+1)=norsq(totsym+1)+1
22    continue
12     continue
      if(odebug(33)) write(6,*) ' norsq ',norsq
      call tit_izo(iofosy,nirr+1)
      do 45 ir=1,nirr
45    iofosy(ir+1)=norsq(ir)*norsq(ir)+iofosy(ir)
      call tit_izo(iofosq,nirr+1)
      do 46 ir=1,nirr
46    iofosq(ir+1)=iofosq(ir)+norsq(ir)
      if(odebug(33)) write(6,*) ' iofosy ',iofosy
c
c
      len=0
      nbkt=1
      bgbkt=0
      fsobkt(1)=0
      call tit_izo(fsobkt,maxbkt)
      call tit_izo(apbkt,maxbkt)
      if(odebug(33)) write(6,*) ' space1 ',space1
      do 401 irr=1,nirr
      do 402 ik=1,norsq(irr)
      if(odebug(33)) write(6,*) ' irr,ik ',irr,ik,iofosq(irr)
      len=len+norsq(irr)
      if(len.gt.space1) then
         nbkt=nbkt+1
         len=len-norsq(irr)
         apbkt(nbkt)=len+apbkt(nbkt-1)
         len=norsq(irr)
         fsobkt(nbkt)=ik-1+iofosq(irr)
      end if
      psobkt(ik+iofosq(irr))=nbkt
402   continue
401   continue
      fsobkt(nbkt+1)=iofosq(nirr+1)
      if(odebug(33)) then
       write(6,*) ' indxol ',indxol
       write(6,*) ' nbkt ',nbkt
       write(6,*) ' fsobkt ',fsobkt
       write(6,*) ' psobkt ',psobkt
       write(6,*) ' apbkt ',apbkt
      endif
      return
      end
      subroutine setvo1(flov,norbo,orbsym,norsq,iofosy,indxol,
     1iofosq,fsobkt,psobkt,no,nbf,nirr,maxbkt,space1,
     1    nbkt,norgsq,npr,nv,norbv,apbkt)
      implicit integer(a-z)
INCLUDE(common/prnprn)
      dimension flov(nirr,4),indxol(no*no),norbv(nirr)
     1 , norsq(nirr),iofosy(nirr+1),orbsym(nbf),norbo(nirr)
      dimension npr(nirr,3,2),apbkt(maxbkt)
      dimension norgsq(nirr,nirr)
      dimension iofosq(nirr+1),fsobkt(maxbkt),psobkt(no*no)
c
c       flov(nirr,4)  -  1      first occupied in symmetry
c                       2      last occupied in symmetry
c                       3      first virtual in symmetry
c                       4      last virtual in symmetry
c          norbo(nirr)    no of occ orbitals in each symmetry type
c          orbsym(nbf)     symmetry type of each orbital (0-nirr-1)
c          norsq(nirr)     no of square pairs of each sym type (occ)
c         norgsq(irr,jrr)  equivalent for squares
c          iofosy(nirr+1)   offset to beginning of each block(ij|kl)
c          iofosq(nirr) offset to beginning of each square (ik)
c          indx(ij) index of ij triangle in symmetry
c          indxl(ij) index of ij square in symmetry
c
      if(odebug(33)) write(6,*) ' flov ',flov
      do 5 i=1,nirr
      norbv(i)=flov(i,4)-flov(i,3)+1
 5    norbo(i)=flov(i,2)-flov(i,1)+1
c      write(6,*) ' norbo ',norbo
      call tit_izo(norsq,nirr)
      call tit_izo(norgsq,nirr*nirr)
      do 12 i=1,nv
      do 22 j=1,no
      totsym=IXOR32(orbsym(i+no),orbsym(j))
      norgsq(orbsym(i+no)+1,totsym+1)=
     1 norgsq(orbsym(i+no)+1,totsym+1)+1
      norsq(totsym+1)=norsq(totsym+1)+1
22    continue
12     continue
c      write(6,*) ' norsq ',norsq
      if(odebug(33)) write(6,*) ' norgsq ',norgsq
      call tit_izo(iofosy,nirr+1)
      do 45 ir=1,nirr
45    iofosy(ir+1)=norsq(ir)*norsq(ir)+iofosy(ir)
      call tit_izo(iofosq,nirr+1)
      do 46 ir=1,nirr
46    iofosq(ir+1)=iofosq(ir)+norsq(ir)
c      write(6,*) ' iofosy ',iofosy
c
      iofst2=0
      do 300 k=1,nirr
      if(k.gt.1) iofst2=iofst2+norsq(k-1)
      iofst1=0
      do 400 l=1,nirr
      if(l.gt.1) iofst1=iofst1+norgsq(l-1,k)
      do 500 i=1,norbv(l)
      nsym=IXOR32(l-1,k-1)+1
      endl=norbo(nsym)
      do 600 j=1,endl
      ij=(i+flov(l,3)-no-2)*no+flov(nsym,1)+j-1
      indxol(ij)=(i-1)*norbo(nsym)+j+iofst1
600   continue
500   continue
400   continue
300   continue
c      write(6,*) ' indxol ',indxol
      if(odebug(33)) write(6,*) ' norsq ',norsq
c
      len=0
      nbkt=1
      bgbkt=0
      fsobkt(1)=0
      call tit_izo(apbkt,maxbkt)
      call tit_izo(fsobkt,maxbkt)
c      write(6,*) ' space1 ',space1
      do 401 irr=1,nirr
      do 402 ik=1,norsq(irr)
      if(odebug(33)) write(6,*) ' irr,ik ',irr,ik,iofosq(irr)
      len=len+norsq(irr)
      if(len.gt.space1) then
         nbkt=nbkt+1
         len=len-norsq(irr)
         apbkt(nbkt)=len+apbkt(nbkt-1)
         len=norsq(irr)
         fsobkt(nbkt)=ik-1+iofosq(irr)
      end if
      psobkt(ik+iofosq(irr))=nbkt
402   continue
401   continue
      fsobkt(nbkt+1)=iofosq(nirr+1)
      if(odebug(33)) then
       write(6,*) ' fsobkt ',fsobkt
       write(6,*) ' psobkt ',psobkt
       write(6,*) ' apbkt ',(apbkt(jki),jki=1,nbkt+1)
      endif
      do 109 i=1,nirr
      npr(i,3,2)=norsq(i)
      npr(i,3,1)=0
109   continue
      return
      end
      subroutine sord1(bkt,ibkt,buf,ibuf,nszbf,nxbkt,length,nbf,
     .      npercl,itap60,jout,itap,indxol,newlen,
     1     norbo,iofosy,norsq,psobkt,iofosq,fsobkt,nirr,no,
     1    mchain,maxval,orbsym,nv,apbkt)
      implicit REAL (a-h,o-z)
INCLUDE(common/prnprn)
      integer orbsym,psobkt,fsobkt,apbkt
      dimension bkt(2),buf(length),ibkt(2),
     .          ibuf(length*2),mchain(nxbkt),orbsym(nbf)
      dimension indxol(no*no),norbo(nirr),apbkt(2),
     1iofosy(nirr+1),norsq(nirr),psobkt(no*nv),iofosq(nirr+1),fsobkt(2)
      data itemp /255/
_IF1()   11 format ('                                                ')
_IF1()   12 format ('************************************************')
   19 format (4i3,f20.12)
      if(odebug(33)) then
       write(6,*) '   in sord1 ',nirr
       write(6,*) ' iofosq ',iofosq
      endif
c
      intlen=(intowp(nszbf)-2)/intowp(1)
      maxval=intowp(intlen)/(1+intowp(1))
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
      if(odebug(33)) then
       write (jout,*) ' nszbf',nszbf
       write (jout,*) ' intlen',intlen
       write (jout,*) ' maxvalues in nszbf ',maxval
       write (jout,*) ' ivoff',ivoff
      endif
c
      ibflen=(intowp(length)-2)/intowp(1)
      maxbuf=intowp(ibflen)/(1+intowp(1))
      iboff = (maxbuf+2)/intowp(1)
      if(intowp(iboff).ne.(maxbuf+2)) iboff = iboff + 1
      if(odebug(33)) write (jout,*) 
     +          ' max values read in , maxbuf = ',maxbuf
      call srew(itap60)
  111 call vclr(buf,1,length)
      call tit_sread(itap60,buf,intowp(length))
c
      call setmbf(iflg,ibuf(1))
c     write (jout,*) ' iflg=',iflg
      call setmbf(mbuf,ibuf(2))
c     write (jout,*) ' mbuf=',mbuf
      do 101 ii=1,mbuf
         call setmbf(ijkl,ibuf(2+ii))
         i=ishft(ijkl,-24)
         j=IAND32(itemp,ishft(ijkl,-16))
         k=IAND32(itemp,ishft(ijkl,-8))
         l=IAND32(itemp,ijkl)
         rint=buf(iboff+ii)
         i=i-no
         k=k-no
c        write (jout,11)
c        write (jout,12)
c        write (jout,*) ' *********  ints as read into sortj *******'
c        write (jout,11)
c        write (jout,12)
         if(odebug(33)) write (jout,19) i,j,k,l,rint
c
c      k([il],k,j)
c
      if(odebug(33)) write(6,*) ' contribution 1 '
      ij=(i-1)*no+j
      kl=(k-1)*no+l
      iposij=indxol(ij)
      iposkl=indxol(kl)
      ijsym=IXOR32(orbsym(no+i),orbsym(j))+1
      iadr=(iposkl-1)*norsq(ijsym)+iposij+iofosy(ijsym)
      ibkt1=psobkt(iposkl+iofosq(ijsym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) write(6,*) ' ij,kl,iadrn,ifill ',ij,kl,iadrn,
     +    ifill ,indxol(ij),indxol(kl)
c
c
c
      iadr=(iposij-1)*norsq(ijsym)+iposkl+iofosy(ijsym)
      ibkt1=psobkt(iposij+iofosq(ijsym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) write(6,*) ' jk/il or kj/li ',kl,ij,
     +   iadrn,ifill,indxol(kl),indxol(ij)
c
c
  101 continue
      if (iflg.eq.0) goto 111
      if(odebug(33)) write(6,*) ' about to write to itape 91 '
      do 22 ibkt1=1,nxbkt
         noff=(ibkt1-1)*nszbf
_IF1()cj       ibkt(intowp(noff)+1)=-1
         call bktdmx(bkt(noff+1),nszbf,itap,ichan)
         mchain(ibkt1)=ichan
   22 continue
c
      return
      end
      subroutine sord2(bkt,ibkt,buf,ibuf,nszbf,nxbkt,length,nbf,
     .      npercl,itap60,jout,itap,indxol,newlen,
     1     norbo,iofosy,norsq,psobkt,iofosq,fsobkt,nirr,no,
     1    mchain,maxval,orbsym,apbkt)
      implicit REAL (a-h,o-z)
      integer orbsym,psobkt,fsobkt,apbkt
INCLUDE(common/prnprn)
      dimension bkt(2),buf(length),ibkt(2),
     .          ibuf(length*2),mchain(nxbkt),orbsym(nbf)
      dimension indxol(no*no),norbo(nirr),apbkt(2),
     1iofosy(nirr+1),norsq(nirr),psobkt(no*no),iofosq(nirr+1),fsobkt(2)
      data itemp /255/
_IF1)   11 format ('                                                ')
_IF1)   12 format ('************************************************')
   19 format (4i3,f20.12)
c      write(*,*) '   in sord2 ',nirr
      if(odebug(33))write(6,*) ' iofosq ',iofosq
c
      intlen=(intowp(nszbf)-2)/intowp(1)
      maxval=intowp(intlen)/(1+intowp(1))
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
      if(odebug(33)) then
       write (jout,*) ' nszbf',nszbf
       write (jout,*) ' intlen',intlen
       write (jout,*) ' maxvalues in nszbf ',maxval
       write (jout,*) ' ivoff',ivoff
      endif
c
      ibflen=(intowp(length)-2)/intowp(1)
      maxbuf=intowp(ibflen)/(1+intowp(1))
      iboff = (maxbuf+2)/intowp(1)
      if(intowp(iboff).ne.(maxbuf+2)) iboff = iboff + 1
      if(odebug(33))
     +  write (jout,*) ' max values read in , maxbuf = ',maxbuf
      call srew(itap60)
  111 call vclr(buf,1,length)
      call tit_sread(itap60,buf,intowp(length))
c
      call setmbf(iflg,ibuf(1))
c     write (jout,*) ' iflg=',iflg
      call setmbf(mbuf,ibuf(2))
c     write (jout,*) ' mbuf=',mbuf
      do 101 ii=1,mbuf
         call setmbf(ijkl,ibuf(2+ii))
         i=ishft(ijkl,-24)
         j=IAND32(itemp,ishft(ijkl,-16))
         k=IAND32(itemp,ishft(ijkl,-8))
         l=IAND32(itemp,ijkl)
         rint=buf(iboff+ii)
         i=i-no
         k=k-no
c        write (jout,11)
c        write (jout,12)
c        write (jout,*) ' *********  ints as read into sortj *******'
c        write (jout,11)
c        write (jout,12)
         if(odebug(33))write (jout,19) i,j,k,l,rint
c
c
c
c
c
c
c
c      k([il],k,j)
c
      if(odebug(33))write(6,*) ' contribution 3 '
      il=(i-1)*no+l
      jk=(k-1)*no+j
      iposil=indxol(il)
      iposjk=indxol(jk)
      jksym=IXOR32(orbsym(j),orbsym(k+no))+1
      ilsym=IXOR32(orbsym(no+i),orbsym(l))+1
_IF1()cj      write(6,*) ' jksym ',jksym,' ilsym ',ilsym
      iadr=(iposil-1)*norsq(jksym)+iposjk+iofosy(jksym)
      ibkt1=psobkt(iposil+iofosq(jksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33))write(6,*) 
     +  ' il,jk,iadrn,ifill ',il,jk,iadrn,ifill,indxol(il),indxol(jk)
c
c
c
c
c     k([kj],i,l)
c
      if(odebug(33))write(6,*) ' contribution 6 '
      jk=(k-1)*no +j
      il=(i-1)*no+l
      iposjk=indxol(jk)
      iposil=indxol(il)
      jksym=IXOR32(orbsym(j),orbsym(k+no))+1
      ilsym=IXOR32(orbsym(no+i),orbsym(l))+1
_IF1()cj      write(6,*) ' jksym ',jksym,' ilsym ',ilsym
      iadr=(iposjk-1)*norsq(ilsym)+iposil+iofosy(ilsym)
      ibkt1=psobkt(iposjk+iofosq(ilsym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
         if(odebug(33))write(6,*) ' jk/il or kj/li ',jk,il,iadrn,ifill
     1    ,indxol(jk),indxol(il)
c
c
  101 continue
      if (iflg.eq.0) goto 111
      if(odebug(33))write(6,*) ' about to write to itape 91 '
      do 22 ibkt1=1,nxbkt
         noff=(ibkt1-1)*nszbf
_IF1()cj       ibkt(intowp(noff)+1)=-1
         call bktdmx(bkt(noff+1),nszbf,itap,ichan)
         mchain(ibkt1)=ichan
   22 continue
c
      return
      end
      subroutine sortd3(bkt,ibkt,buf,ibuf,nszbf,nxbkt,length,nbf,
     .      npercl,itap61,jout,itap,ioff,indx,indxl,newlen,
     1     norbt,iofsy,nprsq,prbkt,iofpr,fpbkt,nirr,nv,ntriv,no,
     1    mchain,maxval,orbsym,apbkt)
      implicit REAL (a-h,o-z)
      integer orbsym,prbkt,fpbkt,apbkt
INCLUDE(common/prnprn)
      dimension bkt(2),buf(length),ibkt(2),apbkt(2),
     .          ibuf(length*2),mchain(nxbkt),orbsym(nbf)
      dimension ioff(nbf),indx(ntriv),indxl(nv*nv),norbt(nirr),
     1   iofsy(nirr),nprsq(nirr),prbkt(ntriv),iofpr(nirr),fpbkt(2)
      data itemp /255/
_IF1()   11 format ('                                                ')
_IF1()   12 format ('************************************************')
   19 format (4i3,f20.12)
c      write(*,*) '   in sortd3 '
c
      intlen=(intowp(nszbf)-2)/intowp(1)
      maxval=intowp(intlen)/(1+intowp(1))
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
      if(odebug(33))then
       write (jout,*) ' nszbf',nszbf
       write (jout,*) ' intlen',intlen
       write (jout,*) ' maxvalues in nszbf ',maxval
       write (jout,*) ' ivoff',ivoff
      endif
c
      ibflen=(intowp(length)-2)/intowp(1)
      maxbuf=intowp(ibflen)/(1+intowp(1))
      iboff = (maxbuf+2)/intowp(1)
      if(intowp(iboff).ne.(maxbuf+2)) iboff = iboff + 1
      if(odebug(33))write (jout,*) 
     +     ' max values read in , maxbuf = ',maxbuf
      call srew(itap61)
  111 call vclr(buf,1,length)
      call tit_sread(itap61,buf,intowp(length))
c
      call setmbf(iflg,ibuf(1))
c     write (jout,*) ' iflg=',iflg
      call setmbf(mbuf,ibuf(2))
c     write (jout,*) ' mbuf=',mbuf
      do 101 ii=1,mbuf
         call setmbf(ijkl,ibuf(2+ii))
         i=ishft(ijkl,-24)
         j=IAND32(itemp,ishft(ijkl,-16))
         k=IAND32(itemp,ishft(ijkl,-8))
         l=IAND32(itemp,ijkl)
         i=i-no
         k=k-no
         rint=buf(iboff+ii)
c        write (jout,11)
c        write (jout,12)
c        write (jout,*) ' *********  ints as read into sortj *******'
c        write (jout,11)
c        write (jout,12)
         if(odebug(33))write (jout,19) i,j,k,l,rint
c
c
c
c
c
c      k([jl],i,k) or k([jl],k,i)
c
c
      if(odebug(33))write(6,*) ' contribution 4 '
      if(j.ge.l) then
      jl=ioff(j)+l
      ik=(i-1)*nv+k
      else
      jl=ioff(l)+j
      ik=(k-1)*nv+i
      end if
      iposjl=indx(jl)
      iposik=indxl(ik)
      jlsym=IXOR32(orbsym(j),orbsym(l))+1
      iksym=IXOR32(orbsym(no+i),orbsym(no+k))+1
_IF1()cj      write(6,*) ' jlsym ',jlsym,' iksym ',iksym
      iadr=(iposjl-1)*nprsq(iksym)+iposik+iofsy(iksym)
      ibkt1=prbkt(iposjl+iofpr(iksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33))write(6,*) ' jl/ik or lj/ki ',jl,ik,iadrn,ifill
     1    ,indx(jl),indxl(ik)
c
c
c
c    if(j.eq.l) then do k([jl],k,i)
c      if(j.ne.l) then  it is included when (il/kj) is read
c
c
      if(j.eq.l.and.i.ne.k) then
      if(odebug(33))write(6,*) ' contribution 5 '
      jl=ioff(j)+l
      ik=(k-1)*nv+i
      iposjl=indx(jl)
      iposik=indxl(ik)
      jlsym=IXOR32(orbsym(j),orbsym(l))+1
      iksym=IXOR32(orbsym(no+i),orbsym(no+k))+1
_IF1()cj      write(6,*) ' jlsym ',jlsym,' iksym ',iksym
      iadr=(iposjl-1)*nprsq(iksym)+iposik+iofsy(iksym)
      ibkt1=prbkt(iposjl+iofpr(iksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33))write(6,*) ' jl,ki,iadrn,ifill ',jl,ik,iadrn,ifill
     1    ,indx(jl),indxl(ik)
      end if
c
c

  101 continue
      if (iflg.eq.0) goto 111
      if(odebug(33))write(6,*) ' about to write to itape 91 '
      do 22 ibkt1=1,nxbkt
         noff=(ibkt1-1)*nszbf
_IF1()cj       ibkt(intowp(noff)+1)=-1
         call bktdmx(bkt(noff+1),nszbf,itap,ichan)
         mchain(ibkt1)=ichan
   22 continue
c
      return
      end
      subroutine sorte1(bkt,ibkt,buf,ibuf,nszbf,nxbkt,length,nbf,
     .      npercl,itap60,jout,itap,indxol,newlen,
     1     iofosy,norsq,psobkt,iofosq,fsobkt,nirr,no,
     1    mchain,maxval,orbsym,isqoo,apbkt,nv)
      implicit REAL (a-h,o-z)
      integer orbsym,psobkt,fsobkt,apbkt
INCLUDE(common/prnprn)
      dimension bkt(2),buf(length),ibkt(2),apbkt(2),
     .          ibuf(length*2),mchain(nxbkt),orbsym(nbf)
      dimension indxol(nv*no),isqoo(no*no),
     1iofosy(nirr+1),norsq(nirr),psobkt(nv*no),iofosq(nirr+1),fsobkt(2)
      data itemp /255/
_IF1()   11 format ('                                                ')
_IF1()   12 format ('************************************************')
   19 format (4i3,f20.12)
      if(odebug(33)) then
       write(6,*) '   in sorte1 ',nirr
       write(6,*) ' iofosq ',iofosq
       write(6,*) ' isqov ',indxol
      endif
c
      intlen=(intowp(nszbf)-2)/intowp(1)
      maxval=intowp(intlen)/(1+intowp(1))
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
      if(odebug(33))then
       write (jout,*) ' nszbf',nszbf
       write (jout,*) ' intlen',intlen
       write (jout,*) ' maxvalues in nszbf ',maxval
       write (jout,*) ' ivoff',ivoff
      endif
c
      ibflen=(intowp(length)-2)/intowp(1)
      maxbuf=intowp(ibflen)/(1+intowp(1))
      iboff = (maxbuf+2)/intowp(1)
      if(intowp(iboff).ne.(maxbuf+2)) iboff = iboff + 1
      if(odebug(33))
     +  write (jout,*) ' max values read in , maxbuf = ',maxbuf
      call srew(itap60)
  111 call vclr(buf,1,length)
      call tit_sread(itap60,buf,intowp(length))
c
      call setmbf(iflg,ibuf(1))
c     write (jout,*) ' iflg=',iflg
      call setmbf(mbuf,ibuf(2))
c     write (jout,*) ' mbuf=',mbuf
      do 101 ii=1,mbuf
         call setmbf(ijkl,ibuf(2+ii))
         i=ishft(ijkl,-24)
         j=IAND32(itemp,ishft(ijkl,-16))
         k=IAND32(itemp,ishft(ijkl,-8))
         l=IAND32(itemp,ijkl)
         rint=buf(iboff+ii)
         i=i-no
c        write (jout,11)
c        write (jout,12)
c        write (jout,*) ' *********  ints as read into sortj *******'
c        write (jout,11)
c        write (jout,12)
         if(odebug(33))write (jout,19) i,j,k,l,rint
c
c
c
c
c       contribution to k ([ik],j,l)
c
      if(odebug(33)) write(6,*) ' contribution 1'
      ij=(i-1)*no+j
      kl=(k-1)*no+l
      iposij=indxol(ij)
      iposkl=isqoo(kl)
      klsym=IXOR32(orbsym(k),orbsym(l))+1
      iadr=(iposij-1)*norsq(klsym)+iposkl+iofosy(klsym)
      if(odebug(33)) write(6,*) ' iposkl,iofosq ',iposkl,iofosq(klsym)
      ibkt1=psobkt(iposij+iofosq(klsym))
      if(odebug(33)) write(6,*) ' ibkt1 ',ibkt1,npercl
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
_IF1()cj    write(6,*) ' iposik,jlsym,iofpr,ibkt1,ixbkt,npercl,ibit ',
_IF1()cj   1  iposik,jlsym,iofpr(jlsym),ibkt1,ixbkt,npercl,ibit
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            mchain(ixbkt)=ichan
            call vclr(bkt(noff+1),1,nszbf)
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) then
       write(6,*) ' klsym,iadr,ibkt1,iadr,iadrx,iadrn ',
     1  klsym,iadr,ibkt1,iadr,iadrx,iadrn
       write(6,*) ' ij,kl,iadrn,ifill ',ij,kl,iadrn,ifill
     1   ,indxol(ij),isqoo(kl)
      endif
c
c
      if(k.ne.l) then
      if(odebug(33)) write(6,*) ' contribution 3 '
      ij=(i-1)*no+j
      kl=(l-1)*no+k
      iposij=indxol(ij)
      iposkl=isqoo(kl)
      klsym=IXOR32(orbsym(k),orbsym(l))+1
      iadr=(iposij-1)*norsq(klsym)+iposkl+iofosy(klsym)
      ibkt1=psobkt(iposij+iofosq(klsym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
         if(odebug(33)) write(6,*) ' ij,kl,iadrn,ifill ',
     +   ij,kl,iadrn,ifill,indxol(ij),isqoo(kl)
      end if
c
  101 continue
      if (iflg.eq.0) goto 111
      if(odebug(33)) write(6,*) ' about to write to itape 91 '
      do 22 ibkt1=1,nxbkt
         noff=(ibkt1-1)*nszbf
_IF1()cj       ibkt(intowp(noff)+1)=-1
         call bktdmx(bkt(noff+1),nszbf,itap,ichan)
         mchain(ibkt1)=ichan
   22 continue
c
      return
      end
      subroutine sorte2(bkt,ibkt,buf,ibuf,nszbf,nxbkt,length,nbf,
     .      npercl,itap60,jout,itap,indxol,newlen,
     1     iofosy,norsq,psobkt,iofosq,fsobkt,nirr,no,
     1    mchain,maxval,orbsym,isqoo,apbkt)
      implicit REAL (a-h,o-z)
      integer orbsym,psobkt,fsobkt,apbkt
INCLUDE(common/prnprn)
      dimension bkt(2),buf(length),ibkt(2),apbkt(2),
     .          ibuf(length*2),mchain(nxbkt),orbsym(nbf)
      dimension indxol(no*no),isqoo(no*no),
     1iofosy(nirr+1),norsq(nirr),psobkt(no*no),iofosq(nirr+1),fsobkt(2)
      data itemp /255/
_IF1()   11 format ('                                                ')
_IF1()   12 format ('************************************************')
   19 format (4i3,f20.12)
      if(odebug(33)) then
       write(*,*) '   in sorte2 ',nirr
       write(6,*) ' iofosq ',iofosq
      endif
c
      intlen=(intowp(nszbf)-2)/intowp(1)
      maxval=intowp(intlen)/(1+intowp(1))
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
      write (jout,*) ' nszbf',nszbf
      write (jout,*) ' intlen',intlen
      write (jout,*) ' maxvalues in nszbf ',maxval
      write (jout,*) ' ivoff',ivoff
c
      ibflen=(intowp(length)-2)/intowp(1)
      maxbuf=intowp(ibflen)/(1+intowp(1))
      iboff = (maxbuf+2)/intowp(1)
      if(intowp(iboff).ne.(maxbuf+2)) iboff = iboff + 1
      if(odebug(33)) write (jout,*) ' max values read in , maxbuf = ',
     +               maxbuf
      call srew(itap60)
  111 call vclr(buf,1,length)
      call tit_sread(itap60,buf,intowp(length))
c
      call setmbf(iflg,ibuf(1))
c     write (jout,*) ' iflg=',iflg
      call setmbf(mbuf,ibuf(2))
c     write (jout,*) ' mbuf=',mbuf
      do 101 ii=1,mbuf
         call setmbf(ijkl,ibuf(2+ii))
         i=ishft(ijkl,-24)
         j=IAND32(itemp,ishft(ijkl,-16))
         k=IAND32(itemp,ishft(ijkl,-8))
         l=IAND32(itemp,ijkl)
         rint=buf(iboff+ii)
         i=i-no
c        write (jout,11)
c        write (jout,12)
c        write (jout,*) ' *********  ints as read into sortj *******'
c        write (jout,11)
c        write (jout,12)
         if(odebug(33))write (jout,19) i,j,k,l,rint
c
c
c
c
c       contribution to k ([ik],j,l)
c
      if(odebug(33)) write(6,*) ' contribution 1'
      ik=(i-1)*no+k
      jl=(j-1)*no+l
      iposik=indxol(ik)
      iposjl=isqoo(jl)
      jlsym=IXOR32(orbsym(j),orbsym(l))+1
      iadr=(iposjl-1)*norsq(jlsym)+iposik+iofosy(jlsym)
      if(odebug(33)) write(6,*) ' iposjl,iofosq ',iposjl,iofosq(jlsym)
      ibkt1=psobkt(iposjl+iofosq(jlsym))
      if(odebug(33)) write(6,*) ' ibkt1 ',ibkt1,npercl
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      if(odebug(33)) write(6,*) 
     +   ' iposik,jlsym,ibkt1,ixbkt,npercl,ibit ',
     +     iposik,jlsym,ibkt1,ixbkt,npercl,ibit
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            mchain(ixbkt)=ichan
            call vclr(bkt(noff+1),1,nszbf)
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) then
        write(6,*) ' jlsym,iadr,ibkt1,iadr,iadrx,iadrn ',
     +  jlsym,iadr,ibkt1,iadr,iadrx,iadrn
        write(6,*) ' ik,jl,iadrn,ifill ',ik,jl,iadrn,ifill
     +   ,indxol(ik),isqoo(jl)
      endif
c
c
      if(k.ne.l) then
      if(odebug(33)) write(6,*) ' contribution 3 '
      il=(i-1)*no+l
      jk=(j-1)*no+k
      iposil=indxol(il)
      iposjk=isqoo(jk)
      jksym=IXOR32(orbsym(k),orbsym(j))+1
      iadr=(iposjk-1)*norsq(jksym)+iposil+iofosy(jksym)
      ibkt1=psobkt(iposjk+iofosq(jksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) write(6,*) ' il,jk,iadrn,ifill ',il,jk,iadrn,ifill
     1    ,indxol(il),isqoo(jk)
      end if
c
  101 continue
      if (iflg.eq.0) goto 111
c      write(6,*) ' about to write to itape 91 '
      do 22 ibkt1=1,nxbkt
         noff=(ibkt1-1)*nszbf
_IF1()cj       ibkt(intowp(noff)+1)=-1
         call bktdmx(bkt(noff+1),nszbf,itap,ichan)
         mchain(ixbkt)=ichan
   22 continue
c
      return
      end
      subroutine sortf(bkt,ibkt,buf,ibuf,nszbf,nxbkt,length,nbf,
     .      npercl,itap60,jout,itap,indxol,newlen,
     1     iofosy,nvrsq,psobkt,iofosq,fsobkt,nirr,no,
     1    mchain,maxval,orbsym,isqvv,apbkt,nv)
      implicit REAL (a-h,o-z)
      integer orbsym,psobkt,fsobkt,apbkt
INCLUDE(common/prnprn)
      dimension bkt(2),buf(length),ibkt(2),apbkt(2),
     .          ibuf(length*2),mchain(nxbkt),orbsym(nbf)
      dimension indxol(nv*no),isqvv(nv*nv),
     1iofosy(nirr+1),nvrsq(nirr),psobkt(nv*no),iofosq(nirr+1),fsobkt(2)
      data itemp /255/
_IF1()   11 format ('                                                ')
_IF1()   12 format ('************************************************')
   19 format (4i3,f20.12)
      if(odebug(33)) then
       write(*,*) '   in sortf ',nirr
       write(6,*) ' iofosq ',iofosq
      endif
c
      intlen=(intowp(nszbf)-2)/intowp(1)
      maxval=intowp(intlen)/(1+intowp(1))
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
      if(odebug(33)) then
       write (jout,*) ' nszbf',nszbf
       write (jout,*) ' intlen',intlen
       write (jout,*) ' maxvalues in nszbf ',maxval
       write (jout,*) ' ivoff',ivoff
      endif
c
      ibflen=(intowp(length)-2)/intowp(1)
      maxbuf=intowp(ibflen)/(1+intowp(1))
      iboff = (maxbuf+2)/intowp(1)
      if(intowp(iboff).ne.(maxbuf+2)) iboff = iboff + 1
      if(odebug(33)) write (jout,*) ' max values read in , maxbuf = '
     +               ,maxbuf
      call srew(itap60)
  111 call vclr(buf,1,length)
      call tit_sread(itap60,buf,intowp(length))
c
      call setmbf(iflg,ibuf(1))
c     write (jout,*) ' iflg=',iflg
      call setmbf(mbuf,ibuf(2))
c     write (jout,*) ' mbuf=',mbuf
      do 101 ii=1,mbuf
         call setmbf(ijkl,ibuf(2+ii))
         i=ishft(ijkl,-24)
         j=IAND32(itemp,ishft(ijkl,-16))
         k=IAND32(itemp,ishft(ijkl,-8))
         l=IAND32(itemp,ijkl)
         rint=buf(iboff+ii)
         k=k-no
         i=i-no
         if(j.le.no) then
            iswp=j
            j=l-no
            l=iswp
            iswp=k
            k=i
            i=iswp
         else
            j=j-no
         end if
c        write (jout,11)
c        write (jout,12)
c        write (jout,*) ' *********  ints as read into sortj *******'
c        write (jout,11)
c        write (jout,12)
         if(odebug(33))write (jout,19) i,j,k,l,rint
c
c
c
c
c       contribution to k ([ik],j,l)
c
      if(odebug(33)) write(6,*) ' contribution 1'
      ik=(i-1)*nv+k
      jl=(j-1)*no+l
      iposik=isqvv(ik)
      iposjl=indxol(jl)
      jlsym=IXOR32(orbsym(i+no),orbsym(k+no))+1
      iadr=(iposjl-1)*nvrsq(jlsym)+iposik+iofosy(jlsym)
      if(odebug(33)) write(6,*) ' iposjl,iofosq ',iposjl,iofosq(jlsym)
      ibkt1=psobkt(iposjl+iofosq(jlsym))
_IF1()cj    write(6,*) ' ibkt1 ',ibkt1,npercl
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
_IF1()cj    write(6,*) ' iposik,jlsym,iofpr,ibkt1,ixbkt,npercl,ibit ',
_IF1()cj   1  iposik,jlsym,iofpr(jlsym),ibkt1,ixbkt,npercl,ibit
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            mchain(ixbkt)=ichan
            call vclr(bkt(noff+1),1,nszbf)
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) then
       write(6,*) ' jlsym,iadr,ibkt1,iadr,iadrx,iadrn ',
     1  jlsym,iadr,ibkt1,iadr,iadrx,iadrn
       write(6,*) ' ik,jl,iadrn,ifill ',ik,jl,iadrn,ifill
     1   ,isqvv(ik),indxol(jl)
      endif
c
c
      if(i.ne.j) then
      if(odebug(33)) write(6,*) ' contribution 3 '
      il=(i-1)*no+l
      jk=(j-1)*nv+k
      iposil=indxol(il)
      iposjk=isqvv(jk)
      jksym=IXOR32(orbsym(k+no),orbsym(j+no))+1
      iadr=(iposil-1)*nvrsq(jksym)+iposjk+iofosy(jksym)
      ibkt1=psobkt(iposil+iofosq(jksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) write(6,*) ' il,jk,iadrn,ifill ',il,jk,iadrn,ifill
     1    ,indxol(il),isqvv(jk)
      end if
c
  101 continue
      if (iflg.eq.0) goto 111
      if(odebug(33)) write(6,*) ' about to write to itape 91 '
      do 22 ibkt1=1,nxbkt
         noff=(ibkt1-1)*nszbf
_IF1()cj       ibkt(intowp(noff)+1)=-1
         call bktdmx(bkt(noff+1),nszbf,itap,ichan)
         mchain(ibkt1)=ichan
   22 continue
c
      return
      end
      subroutine sorti(bkt,ibkt,buf,ibuf,length,ltyp,ibktsp,intbuf,nbf,
     .                 no,itap60,itap61,itap62,itap63,itap64,itap65,
     .                 itap66,itap78,jout,imap,nbfo)
      implicit REAL (a-h,o-z)
INCLUDE(common/prnprn)
      dimension bkt(ibktsp),buf(intbuf),ibkt(ibktsp*2),
     .          ibuf(intbuf*2),imap(nbfo)
      data itemp /255/
_IF1()   11 format ('                                                ')
_IF1()   12 format ('************************************************')
   19 format (4i3,f20.12)
      call vclr(bkt,1,ibktsp)
c
      intlen=(intowp(length)-2)/intowp(1)
      maxval=intowp(intlen)/(1+intowp(1))
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
      if(odebug(33)) then
       write (jout,*) '   in sorti '
       write (jout,*) ' length ',length
       write (jout,*) ' intlen ',intlen
       write (jout,*) ' maxval ',maxval
       write (jout,*) ' ivoff ',ivoff
      endif
c
      ibflen=(intowp(intbuf)-2)/intowp(1)
      maxbuf=intowp(ibflen)/(1+intowp(1))
      iboff = (maxbuf+2)/intowp(1)
      if(intowp(iboff).ne.(maxbuf+2)) iboff = iboff + 1
c     write (jout,*) ' maxbuf=',maxbuf
      call srew(itap78)
  111 call vclr(buf,1,intbuf)
      call tit_sread(itap78,buf,intowp(intbuf))
c
      call setmbf(iflg,ibuf(1))
c     write (jout,*) ' iflg=',iflg
      call setmbf(mbuf,ibuf(2))
c     write (jout,*) ' mbuf=',mbuf
      do 101 ii=1,mbuf
         call setmbf(ijkl,ibuf(2+ii))
         i=ishft(ijkl,-24)
         j=IAND32(itemp,ishft(ijkl,-16))
         k=IAND32(itemp,ishft(ijkl,-8))
         l=IAND32(itemp,ijkl)
         rint=buf(iboff+ii)
c        write (jout,11)
c        write (jout,12)
c        write (jout,*) ' *********  ints as read into sorti *******'
c        write (jout,11)
c        write (jout,12)
         if(odebug(33))write (jout,19) i,j,k,l,rint
         i=imap(i)
         j=imap(j)
         k=imap(k)
         l=imap(l)
         ijkl2 = IOR32(j,ishft(i,8))
         ijkl2 = IOR32(k,ishft(ijkl2,8))
         ijkl2 = IOR32(l,ishft(ijkl2,8))
         if(i.eq.0.or.j.eq.0.or.k.eq.0.or.l.eq.0) go to 101
         if (l.gt.no) then
            if (j.gt.no) then
               ntyp=1
            else
_IF1()cj             ntyp=6
               ntyp=5
            end if
         else if (k.gt.no) then
            if (j.gt.no) then
               ntyp=5
            else
               ntyp=3
            end if
         else if (i.le.no) then
            ntyp=0
         else if (j.gt.no) then
            ntyp=2
         else
            ntyp=4
         end if
c        write (jout,*) ' ntyp',ntyp
         itap=60+ntyp
         noff=ntyp*length
c        write (jout,*) ' noff',noff
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' ntyp',ntyp
c           write (jout,*) ' length',length
            call bktdmp(bkt(noff+1),length,itap)
            call vclr(bkt(noff+1),1,length)
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=ijkl2
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
         if(odebug(33))then
           if (ntyp.eq.1) then
              write (jout,19) i,j,k,l,bkt(ivoff+noff+ifill)
           end if
         endif
c
  101 continue
      if (iflg.eq.0) goto 111
      do 22 jtyp=1,ltyp
         itap=59+jtyp
         noff=(jtyp-1)*length
         ibkt(intowp(noff)+1)=-1
         call bktdmp(bkt(noff+1),length,itap)
   22 continue
c
      return
      end
      subroutine sortj(bkt,ibkt,buf,ibuf,nszbf,nxbkt,length,nbf,
     .      npercl,itap61,jout,itap,ioff,indx,indxl,newlen,
     1     norbt,iofsy,nprsq,prbkt,iofpr,fpbkt,nirr,nv,ntriv,no,
     1    mchain,maxval,orbsym,apbkt)
      implicit REAL (a-h,o-z)
INCLUDE(common/prnprn)
      integer orbsym,prbkt,fpbkt,apbkt
      dimension bkt(2),buf(length),ibkt(2),apbkt(2),
     .          ibuf(length*2),mchain(nxbkt),orbsym(nbf)
      dimension ioff(nbf),indx(ntriv),indxl(nv*nv),norbt(nirr),
     1   iofsy(nirr),nprsq(nirr),prbkt(ntriv),iofpr(nirr),fpbkt(2)
      data itemp /255/
_IF1()   11 format ('                                                ')
_IF1()   12 format ('************************************************')
   19 format (4i3,f20.12)
       if (odebug(33)) write(*,*) '   in sortj '
c      write(6,*) ' orbsym ',orbsym
c
c      write(*,*) ' ioff ',ioff
c      write(*,*) ' indx ',indx
c      write(*,*) ' indxl ',indxl
c      write(*,*) ' norbt ',norbt
c      write(*,*) ' prbkt ',prbkt
c      write(*,*) ' iofpr ',iofpr
c      write(6,*) ' nprsq ',nprsq
c      write(*,*) ' apbkt ',apbkt
c      write(*,*) ' fpbkt ',fpbkt
c      write(6,*) ' iofsy ',iofsy
      intlen=(intowp(nszbf)-2)/intowp(1)
      maxval=intowp(intlen)/(1+intowp(1))
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
c    write (*,*)'nszbf,intlen,maxval,ivoff ',nszbf,intlen,maxval,ivoff
      if(odebug(33)) then
       write (jout,*) ' intlen',intlen
       write (jout,*) ' maxvalues in nszbf ',maxval
       write (jout,*) ' ivoff',ivoff
      endif
c
      ibflen=(intowp(length)-2)/intowp(1)
      maxbuf=intowp(ibflen)/(1+intowp(1))
      iboff = (maxbuf+2)/intowp(1)
      if(intowp(iboff).ne.(maxbuf+2)) iboff = iboff + 1
c      write (*,*)' length,maxbuf,iboff ',length,maxbuf,iboff
      if(odebug(33)) then
        write (jout,*) ' max values read in , maxbuf = ',maxbuf
      endif
      call srew(itap61)
  111 call vclr(buf,1,length)
      call tit_sread(itap61,buf,intowp(length))
c
      call setmbf(iflg,ibuf(1))
c     write (jout,*) ' iflg=',iflg
      call setmbf(mbuf,ibuf(2))
c     write (jout,*) ' mbuf=',mbuf
      do 101 ii=1,mbuf
         call setmbf(ijkl,ibuf(2+ii))
         i=ishft(ijkl,-24)
         j=IAND32(itemp,ishft(ijkl,-16))
         k=IAND32(itemp,ishft(ijkl,-8))
         l=IAND32(itemp,ijkl)
         i=i-no
         j=j-no
         k=k-no
         l=l-no
         rint=buf(iboff+ii)
c        write (jout,11)
c        write (jout,12)
c        write (jout,*) ' *********  ints as read into sortj *******'
c        write (jout,11)
c        write (jout,12)
         if(odebug(33))write (jout,19) i,j,k,l,rint
c
c        contribution to k ([ik],j,l)
c
_IF1()ctjl      write(6,*) ' contribution 1'
         ik=ioff(i)+k
         jl=(j-1)*nv+l
         iposik=indx(ik)
         iposjl=indxl(jl)
         jlsym=IXOR32(orbsym(no+j),orbsym(no+l))+1
         iksym=IXOR32(orbsym(no+i),orbsym(no+k))+1
         if(odebug(33)) write(6,*) ' jlsym ',jlsym,' iksym ',iksym
         iadr=(iposik-1)*nprsq(jlsym)+iposjl+iofsy(jlsym)
         ibkt1=prbkt(iposik+iofpr(jlsym))
         if(odebug(33)) write(6,*) ' ibkt1 ',ibkt1,npercl
         ixbkt=(ibkt1-1)/npercl+1
         ibit=mod(ibkt1-1,npercl)
         if(odebug(33)) 
     +     write(6,*) ' iposik,jlsym,iofpr,ibkt1,ixbkt,npercl,ibit ',
     +     iposik,jlsym,iofpr(jlsym),ibkt1,ixbkt,npercl,ibit
         iofst=ibit*newlen
         iadrx=iadr-apbkt(ibkt1)
         iadrn=iadrx+iofst
c
c        write out value in ixbkt core load in ibkt buffer
c        address in core load = iadrn
c        address in buffer = iadrx
c        original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            mchain(ixbkt)=ichan
            call vclr(bkt(noff+1),1,nszbf)
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
         if(odebug(33)) write(6,*) 
     +     ' ik,jl,iadrn,ifill ',ik,jl,iadrn,ifill
     1      ,indx(ik),indxl(jl),iadr
c
c
c        if i.eq.k contribution to k([ik],l,j) - if i.ne.k this is done
c        when (il/kj) is read
c
         if(i.eq.k.and.j.ne.l) then
_IF1()ctjl      write(6,*) ' contribution 2'
            ik=ioff(i)+k
            jl=(l-1)*nv+j
            iposik=indx(ik)
            iposjl=indxl(jl)
            jlsym=IXOR32(orbsym(no+j),orbsym(no+l))+1
            iksym=IXOR32(orbsym(no+i),orbsym(no+k))+1
_IF1()cj      write(6,*) ' jlsym ',jlsym,' iksym ',iksym
            iadr=(iposik-1)*nprsq(jlsym)+iposjl+iofsy(jlsym)
            ibkt1=prbkt(iposik+iofpr(jlsym))
            ixbkt=(ibkt1-1)/npercl+1
            ibit=mod(ibkt1-1,npercl)
            iofst=ibit*newlen
            iadrx=iadr-apbkt(ibkt1)
            iadrn=iadrx+iofst
c
c           write out value in ixbkt core load in ibkt buffer
c           address in core load = iadrn
c           address in buffer = iadrx
c           original address = iadr
c
            noff=(ixbkt-1)*nszbf
            ifill=ibkt(intowp(noff)+2)+1
c           write (jout,*) ' ifill',ifill
            if (ifill.gt.maxval) then
c              write (jout,*) ' bucket dumped'
c              write (jout,*) ' ifill',ifill
c              write (jout,*) ' nszbf',nszbf
               call bktdmx(bkt(noff+1),nszbf,itap,ichan)
               call vclr(bkt(noff+1),1,nszbf)
               mchain(ixbkt)=ichan
               ibkt(intowp(noff)+1)=ichan
               ifill=1
            end if
            ibkt(intowp(noff)+2+ifill)=iadrn
            bkt(ivoff+noff+ifill)=rint
            ibkt(intowp(noff)+2)=ifill
            if(odebug(33)) write(6,*) ' ik,lj,iadrn,ifill ',ik,jl,iadrn,
     &         ifill,indx(ik),indxl(jl)
         end if
c
c        k([il],j,k) - if k.eq.l then already included in k([ik],j,l)
c
         if(k.ne.l) then
_IF1()ctjl      write(6,*) ' contribution 3 '
            il=ioff(i)+l
            jk=(j-1)*nv+k
            iposil=indx(il)
            iposjk=indxl(jk)
            jksym=IXOR32(orbsym(no+j),orbsym(no+k))+1
            ilsym=IXOR32(orbsym(no+i),orbsym(no+l))+1
_IF1()cj      write(6,*) ' jksym ',jksym,' ilsym ',ilsym
            iadr=(iposil-1)*nprsq(jksym)+iposjk+iofsy(jksym)
            ibkt1=prbkt(iposil+iofpr(jksym))
            ixbkt=(ibkt1-1)/npercl+1
            ibit=mod(ibkt1-1,npercl)
            iofst=ibit*newlen
            iadrx=iadr-apbkt(ibkt1)
            iadrn=iadrx+iofst
c
c           write out value in ixbkt core load in ibkt buffer
c           address in core load = iadrn
c           address in buffer = iadrx
c           original address = iadr
c
            noff=(ixbkt-1)*nszbf
            ifill=ibkt(intowp(noff)+2)+1
c           write (jout,*) ' ifill',ifill
            if (ifill.gt.maxval) then
c              write (jout,*) ' bucket dumped'
c              write (jout,*) ' ifill',ifill
c              write (jout,*) ' nszbf',nszbf
               call bktdmx(bkt(noff+1),nszbf,itap,ichan)
               call vclr(bkt(noff+1),1,nszbf)
               mchain(ixbkt)=ichan
               ibkt(intowp(noff)+1)=ichan
               ifill=1
            end if
            ibkt(intowp(noff)+2+ifill)=iadrn
            bkt(ivoff+noff+ifill)=rint
            ibkt(intowp(noff)+2)=ifill
            if(odebug(33)) write(6,*) ' il,jk,iadrn,ifill ',il,jk,iadrn,
     1         ifill,indx(il),indxl(jk)
         end if
c
c
c        k([jl],i,k) or k([jl],k,i)
c
c
_IF1()ctjl      write(6,*) ' contribution 4 '
         if(j.ge.l) then
            jl=ioff(j)+l
            ik=(i-1)*nv+k
         else
            jl=ioff(l)+j
            ik=(k-1)*nv+i
         end if
         iposjl=indx(jl)
         iposik=indxl(ik)
         jlsym=IXOR32(orbsym(no+j),orbsym(no+l))+1
         iksym=IXOR32(orbsym(no+i),orbsym(no+k))+1
_IF1()cj      write(6,*) ' jlsym ',jlsym,' iksym ',iksym
         iadr=(iposjl-1)*nprsq(iksym)+iposik+iofsy(iksym)
         ibkt1=prbkt(iposjl+iofpr(iksym))
         ixbkt=(ibkt1-1)/npercl+1
         ibit=mod(ibkt1-1,npercl)
         iofst=ibit*newlen
         iadrx=iadr-apbkt(ibkt1)
         iadrn=iadrx+iofst
c
c        write out value in ixbkt core load in ibkt buffer
c        address in core load = iadrn
c        address in buffer = iadrx
c        original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
_IF1()ctjl         write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
         if(odebug(33)) write(6,*) ' jl/ik or lj/ki ',jl,ik,iadrn,ifill
     +    ,indx(jl),indxl(ik)
c
c
c        if(j.eq.l) then do k([jl],k,i)
c        if(j.ne.l) then  it is included when (il/kj) is read
c
c
         if(j.eq.l.and.i.ne.k) then
_IF1()ctjl      write(6,*) ' contribution 5 '
            jl=ioff(j)+l
            ik=(k-1)*nv+i
            iposjl=indx(jl)
            iposik=indxl(ik)
            jlsym=IXOR32(orbsym(no+j),orbsym(no+l))+1
            iksym=IXOR32(orbsym(no+i),orbsym(no+k))+1
_IF1()cj      write(6,*) ' jlsym ',jlsym,' iksym ',iksym
            iadr=(iposjl-1)*nprsq(iksym)+iposik+iofsy(iksym)
            ibkt1=prbkt(iposjl+iofpr(iksym))
            ixbkt=(ibkt1-1)/npercl+1
            ibit=mod(ibkt1-1,npercl)
            iofst=ibit*newlen
            iadrx=iadr-apbkt(ibkt1)
            iadrn=iadrx+iofst
c
c           write out value in ixbkt core load in ibkt buffer
c           address in core load = iadrn
c           address in buffer = iadrx
c           original address = iadr
c
            noff=(ixbkt-1)*nszbf
            ifill=ibkt(intowp(noff)+2)+1
c           write (jout,*) ' ifill',ifill
            if (ifill.gt.maxval) then
c              write (jout,*) ' bucket dumped'
c              write (jout,*) ' ifill',ifill
c              write (jout,*) ' nszbf',nszbf
               call bktdmx(bkt(noff+1),nszbf,itap,ichan)
               call vclr(bkt(noff+1),1,nszbf)
               mchain(ixbkt)=ichan
               ibkt(intowp(noff)+1)=ichan
               ifill=1
            end if
            ibkt(intowp(noff)+2+ifill)=iadrn
            bkt(ivoff+noff+ifill)=rint
            ibkt(intowp(noff)+2)=ifill
            if(odebug(33)) write(6,*) ' jl,ki,iadrn,ifill ',
     +         jl,ik,iadrn,ifill,indx(jl),indxl(ik)
         end if
c
c
c
c        k([jk],i,l) or k([jk],l,i)
c        if(i.eq.j) then contribution already covered
c
c
         if(i.ne.j) then
_IF1()ctjl      write(6,*) ' contribution 6 '
            if(j.ge.k) then
               jk=ioff(j)+k
               il=(i-1)*nv+l
            else
               jk=ioff(k)+j
               il=(l-1)*nv+i
            end if
            iposjk=indx(jk)
            iposil=indxl(il)
            jksym=IXOR32(orbsym(no+j),orbsym(no+k))+1
            ilsym=IXOR32(orbsym(no+i),orbsym(no+l))+1
_IF1()cj      write(6,*) ' jksym ',jksym,' ilsym ',ilsym
            iadr=(iposjk-1)*nprsq(ilsym)+iposil+iofsy(ilsym)
            ibkt1=prbkt(iposjk+iofpr(ilsym))
            ixbkt=(ibkt1-1)/npercl+1
            ibit=mod(ibkt1-1,npercl)
            iofst=ibit*newlen
            iadrx=iadr-apbkt(ibkt1)
            iadrn=iadrx+iofst
c
c           write out value in ixbkt core load in ibkt buffer
c           address in core load = iadrn
c           address in buffer = iadrx
c           original address = iadr
c
            noff=(ixbkt-1)*nszbf
            ifill=ibkt(intowp(noff)+2)+1
c           write (jout,*) ' ifill',ifill
            if (ifill.gt.maxval) then
c              write (jout,*) ' bucket dumped'
c              write (jout,*) ' ifill',ifill
c              write (jout,*) ' nszbf',nszbf
               call bktdmx(bkt(noff+1),nszbf,itap,ichan)
               call vclr(bkt(noff+1),1,nszbf)
               mchain(ixbkt)=ichan
               ibkt(intowp(noff)+1)=ichan
               ifill=1
            end if
            ibkt(intowp(noff)+2+ifill)=iadrn
            bkt(ivoff+noff+ifill)=rint
            ibkt(intowp(noff)+2)=ifill
            if(odebug(33)) write(6,*) ' jk/il or kj/li ',jk,il,iadrn,
     +         ifill,indx(jk),indxl(il)
         end if
c
c
c        if(j.eq.k) then do k([jk],l,i)
c        if(j.ne.k) then included when (ik/jl) is read
c
c
c
         if(j.eq.k.and.i.ne.l) then
_IF1()ctjl      write(6,*) ' contribution 7 '
            jk=ioff(j)+k
            il=(l-1)*nv+i
            iposjk=indx(jk)
            iposil=indxl(il)
            jksym=IXOR32(orbsym(no+j),orbsym(no+k))+1
            ilsym=IXOR32(orbsym(no+i),orbsym(no+l))+1
_IF1()cj    write(6,*) ' jksym ',jksym,' ilsym ',ilsym
            iadr=(iposjk-1)*nprsq(ilsym)+iposil+iofsy(ilsym)
            ibkt1=prbkt(iposjk+iofpr(ilsym))
            ixbkt=(ibkt1-1)/npercl+1
            ibit=mod(ibkt1-1,npercl)
            iofst=ibit*newlen
            iadrx=iadr-apbkt(ibkt1)
            iadrn=iadrx+iofst
c
c           write out value in ixbkt core load in ibkt buffer
c           address in core load = iadrn
c           address in buffer = iadrx
c           original address = iadr
c
            noff=(ixbkt-1)*nszbf
            ifill=ibkt(intowp(noff)+2)+1
c           write (jout,*) ' ifill',ifill
            if (ifill.gt.maxval) then
c              write (jout,*) ' bucket dumped'
c              write (jout,*) ' ifill',ifill
c              write (jout,*) ' nszbf',nszbf
               call bktdmx(bkt(noff+1),nszbf,itap,ichan)
               call vclr(bkt(noff+1),1,nszbf)
               mchain(ixbkt)=ichan
               ibkt(intowp(noff)+1)=ichan
               ifill=1
            end if
            ibkt(intowp(noff)+2+ifill)=iadrn
            bkt(ivoff+noff+ifill)=rint
            ibkt(intowp(noff)+2)=ifill
            if(odebug(33)) write(6,*) ' jk,il,iadrn,ifill ',
     +         jk,il,iadrn,ifill,indx(jk),indxl(il)
         end if
c
  101 continue
      if (iflg.eq.0) goto 111
      do 22 ibkt1=1,nxbkt
         noff=(ibkt1-1)*nszbf
_IF1()cj       ibkt(intowp(noff)+1)=-1
         call bktdmx(bkt(noff+1),nszbf,itap,ichan)
         mchain(ibkt1)=ichan
   22 continue
c
      return
      end
      subroutine sorto(bkt,ibkt,buf,ibuf,nszbf,nxbkt,length,nbf,
     .      npercl,itap60,jout,itap,indxol,newlen,
     1     norbo,iofosy,norsq,psobkt,iofosq,fsobkt,nirr,no,
     1    mchain,maxval,orbsym,apbkt)
      implicit REAL (a-h,o-z)
      integer orbsym,psobkt,fsobkt,apbkt
INCLUDE(common/prnprn)
      dimension bkt(2),buf(length),ibkt(2),apbkt(5),
     .          ibuf(length*2),mchain(nxbkt),orbsym(nbf)
      dimension indxol(no*no),norbo(nirr),
     1iofosy(nirr+1),norsq(nirr),psobkt(no*no),iofosq(nirr+1),fsobkt(2)
      data itemp /255/
_IF1()   11 format ('                                                ')
_IF1()   12 format ('************************************************')
   19 format (4i3,f20.12)
      if(odebug(33)) then
       write(6,*) '   in sorto ',nirr
       write(6,*) ' iofosq ',iofosq
       write(6,*) ' apbkt ',apbkt
      endif
c
      intlen=(intowp(nszbf)-2)/intowp(1)
      maxval=intowp(intlen)/(1+intowp(1))
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
      if(odebug(33)) then
       write (jout,*) ' nszbf',nszbf
       write (jout,*) ' intlen',intlen
       write (jout,*) ' maxvalues in nszbf ',maxval
       write (jout,*) ' ivoff',ivoff
      endif
c
      ibflen=(intowp(length)-2)/intowp(1)
      maxbuf=intowp(ibflen)/(1+intowp(1))
      iboff = (maxbuf+2)/intowp(1)
      if(intowp(iboff).ne.(maxbuf+2)) iboff = iboff + 1
      if(odebug(33))  write (jout,*) ' max values read in , maxbuf = '
     +                ,maxbuf
      call srew(itap60)
  111 call vclr(buf,1,length)
      call tit_sread(itap60,buf,intowp(length))
c
      call setmbf(iflg,ibuf(1))
c     write (jout,*) ' iflg=',iflg
      call setmbf(mbuf,ibuf(2))
c     write (jout,*) ' mbuf=',mbuf
      do 101 ii=1,mbuf
         call setmbf(ijkl,ibuf(2+ii))
         i=ishft(ijkl,-24)
         j=IAND32(itemp,ishft(ijkl,-16))
         k=IAND32(itemp,ishft(ijkl,-8))
         l=IAND32(itemp,ijkl)
         rint=buf(iboff+ii)
c        write (jout,11)
c        write (jout,12)
c        write (jout,*) ' *********  ints as read into sortj *******'
c        write (jout,11)
c        write (jout,12)
         if(odebug(33))write (jout,19) i,j,k,l,rint
c
c       contribution to k ([ik],j,l)
c
      if(odebug(33)) write(6,*) ' contribution 1'
      ik=(k-1)*no+i
      jl=(l-1)*no+j
      iposik=indxol(ik)
      iposjl=indxol(jl)
      jlsym=IXOR32(orbsym(j),orbsym(l))+1
      iksym=IXOR32(orbsym(i),orbsym(k))+1
_IF1()cj      write(6,*) ' jlsym ',jlsym,' iksym ',iksym
      iadr=(iposik-1)*norsq(jlsym)+iposjl+iofosy(jlsym)
      if(odebug(33)) write(6,*) ' iposik,iofosy ',
     +               iposik,iofosy(jlsym),norsq(jlsym)
      ibkt1=psobkt(iposik+iofosq(jlsym))
_IF1()cj    write(6,*) ' ibkt1 ',ibkt1,npercl
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
_IF1()cj    write(6,*) ' iposik,jlsym,iofpr,ibkt1,ixbkt,npercl,ibit ',
_IF1()cj   1  iposik,jlsym,iofpr(jlsym),ibkt1,ixbkt,npercl,ibit
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            mchain(ixbkt)=ichan
            call vclr(bkt(noff+1),1,nszbf)
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) then
         write(6,*) ' jlsym,iksym,iadr,ibkt1,iadr,iadrx,iadrn ',
     1  jlsym,iksym,iadr,ibkt1,iadr,iadrx,iadrn
         write(6,*) ' ik,jl,iadrn,ifill ',ik,jl,iadrn,ifill
     1   ,indxol(ik),indxol(jl)
      endif
c
_IF1()cj    if(i.ne.k) then
      ik=(i-1)*no+k
      jl=(j-1)*no+l
      iposik=indxol(ik)
      iposjl=indxol(jl)
      jlsym=IXOR32(orbsym(j),orbsym(l))+1
      iksym=IXOR32(orbsym(i),orbsym(k))+1
_IF1()cj      write(6,*) ' jlsym ',jlsym,' iksym ',iksym
      iadr=(iposik-1)*norsq(jlsym)+iposjl+iofosy(jlsym)
      ibkt1=psobkt(iposik+iofosq(jlsym))
_IF1()cj    write(6,*) ' ibkt1 ',ibkt1,npercl
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
_IF1()cj    write(6,*) ' iposik,jlsym,iofpr,ibkt1,ixbkt,npercl,ibit ',
_IF1()cj   1  iposik,jlsym,iofpr(jlsym),ibkt1,ixbkt,npercl,ibit
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            mchain(ixbkt)=ichan
            call vclr(bkt(noff+1),1,nszbf)
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
         if(odebug(33)) write(6,*) ' ik,jl,iadrn,ifill ',
     +         ik,jl,iadrn,ifill,indxol(ik),indxol(jl)
c
c      k([il],j,k) - if k.eq.l then already included in k([ik],j,l)
c
      if(k.ne.l) then
      if(odebug(33)) write(6,*) ' contribution 3 '
      il=(l-1)*no+i
      jk=(k-1)*no+j
      iposil=indxol(il)
      iposjk=indxol(jk)
      jksym=IXOR32(orbsym(j),orbsym(k))+1
      ilsym=IXOR32(orbsym(i),orbsym(l))+1
_IF1()cj      write(6,*) ' jksym ',jksym,' ilsym ',ilsym
      iadr=(iposil-1)*norsq(jksym)+iposjk+iofosy(jksym)
      ibkt1=psobkt(iposil+iofosq(jksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) write(6,*) ' il,jk,iadrn,ifill ',il,jk,iadrn,ifill
     1    ,indxol(il),indxol(jk)
      end if
c
      if(k.ne.l) then
      if(odebug(33)) write(6,*) ' contribution 3 '
      il=(i-1)*no+l
      jk=(j-1)*no+k
      iposil=indxol(il)
      iposjk=indxol(jk)
      jksym=IXOR32(orbsym(j),orbsym(k))+1
      ilsym=IXOR32(orbsym(i),orbsym(l))+1
_IF1()cj      write(6,*) ' jksym ',jksym,' ilsym ',ilsym
      iadr=(iposil-1)*norsq(jksym)+iposjk+iofosy(jksym)
      ibkt1=psobkt(iposil+iofosq(jksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) write(6,*) ' il,jk,iadrn,ifill ',il,jk,iadrn,ifill
     1    ,indxol(il),indxol(jk)
      end if
c
      if(k.ne.l) then
      if(odebug(33)) write(6,*) ' contribution 3 '
      il=(i-1)*no+l
      jk=(j-1)*no+k
      iposil=indxol(il)
      iposjk=indxol(jk)
      jksym=IXOR32(orbsym(j),orbsym(k))+1
      ilsym=IXOR32(orbsym(i),orbsym(l))+1
_IF1()cj      write(6,*) ' jksym ',jksym,' ilsym ',ilsym
      iadr=(iposjk-1)*norsq(jksym)+iposil+iofosy(jksym)
      ibkt1=psobkt(iposjk+iofosq(jksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) write(6,*) ' il,jk,iadrn,ifill ',il,jk,iadrn,ifill
     1    ,indxol(il),indxol(jk)
      end if
c
c
c
c
c      k([jl],i,k) or k([jl],k,i)
c
c
      if(odebug(33)) write(6,*) ' contribution 4 '
      jl=(l-1)*no+j
      ik=(k-1)*no+i
      iposjl=indxol(jl)
      iposik=indxol(ik)
      jlsym=IXOR32(orbsym(j),orbsym(l))+1
      iksym=IXOR32(orbsym(i),orbsym(k))+1
_IF1()cj      write(6,*) ' jlsym ',jlsym,' iksym ',iksym
      iadr=(iposjl-1)*norsq(iksym)+iposik+iofosy(iksym)
      ibkt1=psobkt(iposjl+iofosq(iksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) write(6,*) ' jl/ik or lj/ki ',jl,ik,iadrn,ifill
     1    ,indxol(jl),indxol(ik)
c
c
      if(odebug(33)) write(6,*) ' contribution 4 '
      jl=(j-1)*no+l
      ik=(i-1)*no+k
      iposjl=indxol(jl)
      iposik=indxol(ik)
      jlsym=IXOR32(orbsym(j),orbsym(l))+1
      iksym=IXOR32(orbsym(i),orbsym(k))+1
_IF1()cj      write(6,*) ' jlsym ',jlsym,' iksym ',iksym
      iadr=(iposjl-1)*norsq(iksym)+iposik+iofosy(iksym)
      ibkt1=psobkt(iposjl+iofosq(iksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) write(6,*) ' jl/ik or lj/ki ',jl,ik,iadrn,ifill
     1    ,indxol(jl),indxol(ik)
c
c
c
c
      if(i.ne.j) then
      if(odebug(33)) write(6,*) ' contribution 6 '
      jk=(k-1)*no +j
      il=(l-1)*no+i
      iposjk=indxol(jk)
      iposil=indxol(il)
      jksym=IXOR32(orbsym(j),orbsym(k))+1
      ilsym=IXOR32(orbsym(i),orbsym(l))+1
_IF1()cj      write(6,*) ' jksym ',jksym,' ilsym ',ilsym
      iadr=(iposjk-1)*norsq(ilsym)+iposil+iofosy(ilsym)
      ibkt1=psobkt(iposjk+iofosq(ilsym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) write(6,*) ' jk/il or kj/li ',jk,il,iadrn,ifill
     1    ,indxol(jk),indxol(il)
      end if
c
      if(i.ne.j) then
      if(odebug(33))write(6,*) ' contribution 6 '
      jk=(j-1)*no +k
      il=(i-1)*no+l
      iposjk=indxol(jk)
      iposil=indxol(il)
      jksym=IXOR32(orbsym(j),orbsym(k))+1
      ilsym=IXOR32(orbsym(i),orbsym(l))+1
_IF1()cj      write(6,*) ' jksym ',jksym,' ilsym ',ilsym
      iadr=(iposjk-1)*norsq(ilsym)+iposil+iofosy(ilsym)
      ibkt1=psobkt(iposjk+iofosq(ilsym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33)) write(6,*) ' jk/il or kj/li ',jk,il,iadrn,ifill
     1    ,indxol(jk),indxol(il)
      end if
c
  101 continue
      if (iflg.eq.0) goto 111
      if(odebug(33)) write(6,*) ' about to write to itape 91 '
      do 22 ibkt1=1,nxbkt
         noff=(ibkt1-1)*nszbf
_IF1()cj       ibkt(intowp(noff)+1)=-1
         call bktdmx(bkt(noff+1),nszbf,itap,ichan)
         mchain(ibkt1)=ichan
   22 continue
c
      return
      end
      subroutine sorvo1(bkt,ibkt,buf,ibuf,nszbf,nxbkt,length,nbf,
     .      npercl,itap60,jout,itap,indxol,newlen,
     1     norbo,iofosy,norsq,psobkt,iofosq,fsobkt,nirr,no,
     1    mchain,maxval,orbsym,apbkt)
      implicit REAL (a-h,o-z)
      integer orbsym,psobkt,fsobkt,apbkt
INCLUDE(common/prnprn)
      dimension bkt(2),buf(length),ibkt(2),apbkt(2),
     .          ibuf(length*2),mchain(nxbkt),orbsym(nbf)
      dimension indxol(no*no),norbo(nirr),
     1iofosy(nirr+1),norsq(nirr),psobkt(no*no),iofosq(nirr+1),fsobkt(2)
      data itemp /255/
_IF1()   11 format ('                                                ')
_IF1()   12 format ('************************************************')
   19 format (4i3,f20.12)
      if(odebug(33)) then
       write(6,*) ' iofosq ',iofosq
       write(6,*) '   in sorvo1 ',nirr
      endif
c
      intlen=(intowp(nszbf)-2)/intowp(1)
      maxval=intowp(intlen)/(1+intowp(1))
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
      if(odebug(33)) then
       write (jout,*) ' nszbf',nszbf
       write (jout,*) ' intlen',intlen
       write (jout,*) ' maxvalues in nszbf ',maxval
       write (jout,*) ' ivoff',ivoff
      endif
c
      ibflen=(intowp(length)-2)/intowp(1)
      maxbuf=intowp(ibflen)/(1+intowp(1))
      iboff = (maxbuf+2)/intowp(1)
      if(intowp(iboff).ne.(maxbuf+2)) iboff = iboff + 1
      if(odebug(33)) write (jout,*) ' max values read in , maxbuf = '
     +              ,maxbuf
      call srew(itap60)
  111 call vclr(buf,1,length)
      call tit_sread(itap60,buf,intowp(length))
c
      call setmbf(iflg,ibuf(1))
c     write (jout,*) ' iflg=',iflg
      call setmbf(mbuf,ibuf(2))
c     write (jout,*) ' mbuf=',mbuf
      do 101 ii=1,mbuf
         call setmbf(ijkl,ibuf(2+ii))
         i=ishft(ijkl,-24)
         j=IAND32(itemp,ishft(ijkl,-16))
         k=IAND32(itemp,ishft(ijkl,-8))
         l=IAND32(itemp,ijkl)
         rint=buf(iboff+ii)
         i=i-no
         j=j-no
c        write (jout,11)
c        write (jout,12)
c        write (jout,*) ' *********  ints as read into sortj *******'
c        write (jout,11)
c        write (jout,12)
         if(odebug(33))write (jout,19) i,j,k,l,rint
c
c
c
c
c       contribution to k ([ik],j,l)
c
      if(odebug(33))write(6,*) ' contribution 1'
      ik=(i-1)*no+k
      jl=(j-1)*no+l
      iposik=indxol(ik)
      iposjl=indxol(jl)
      jlsym=IXOR32(orbsym(j+no),orbsym(l))+1
      iksym=IXOR32(orbsym(i+no),orbsym(k))+1
_IF1()cj      write(6,*) ' jlsym ',jlsym,' iksym ',iksym
      iadr=(iposik-1)*norsq(jlsym)+iposjl+iofosy(jlsym)
      if(odebug(33))write(6,*) ' iposik,iofosq ',iposik,iofosq(jlsym)
      ibkt1=psobkt(iposik+iofosq(jlsym))
      if(odebug(33))write(6,*) ' ibkt1 ',ibkt1,npercl
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
_IF1()cj    write(6,*) ' iposik,jlsym,iofpr,ibkt1,ixbkt,npercl,ibit ',
_IF1()cj   1  iposik,jlsym,iofpr(jlsym),ibkt1,ixbkt,npercl,ibit
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            mchain(ixbkt)=ichan
            call vclr(bkt(noff+1),1,nszbf)
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
         if(odebug(33)) then
          write(6,*) ' jlsym,iksym,iadr,ibkt1,iadr,iadrx,iadrn ',
     1                 jlsym,iksym,iadr,ibkt1,iadr,iadrx,iadrn
          write(6,*) ' ik,jl,iadrn,ifill ',ik,jl,iadrn,ifill
     1                ,indxol(ik),indxol(jl)
         endif
c
c
c
c
c      k([il],j,k) - if k.eq.l then already included in k([ik],j,l)
c
      if(k.ne.l) then
      if(odebug(33))write(6,*) ' contribution 3 '
      il=(i-1)*no+l
      jk=(j-1)*no+k
      iposil=indxol(il)
      iposjk=indxol(jk)
      jksym=IXOR32(orbsym(no+j),orbsym(k))+1
      ilsym=IXOR32(orbsym(no+i),orbsym(l))+1
_IF1()cj      write(6,*) ' jksym ',jksym,' ilsym ',ilsym
      iadr=(iposil-1)*norsq(jksym)+iposjk+iofosy(jksym)
      ibkt1=psobkt(iposil+iofosq(jksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33))write(6,*) ' il,jk,iadrn,ifill ',il,jk,iadrn,ifill
     1    ,indxol(il),indxol(jk)
      end if
c
c
c
c
c
c      k([jl],i,k) or k([jl],k,i)
c
c
      if(i.ne.j) then
      if(odebug(33))write(6,*) ' contribution 4 '
      jl=(j-1)*no+l
      ik=(i-1)*no+k
      iposjl=indxol(jl)
      iposik=indxol(ik)
      jlsym=IXOR32(orbsym(no+j),orbsym(l))+1
      iksym=IXOR32(orbsym(no+i),orbsym(k))+1
_IF1()cj      write(6,*) ' jlsym ',jlsym,' iksym ',iksym
      iadr=(iposjl-1)*norsq(iksym)+iposik+iofosy(iksym)
      ibkt1=psobkt(iposjl+iofosq(iksym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33))write(6,*) ' jl/ik or lj/ki ',jl,ik,iadrn,ifill
     1    ,indxol(jl),indxol(ik)
      end if
c
c
c     k([jk],i,l) or k([jk],l,i)
c    if(i.eq.j) then contribution already covered
c
c
      if(i.ne.j.and.k.ne.l) then
      if(odebug(33))write(6,*) ' contribution 6 '
      jk=(j-1)*no +k
      il=(i-1)*no+l
      iposjk=indxol(jk)
      iposil=indxol(il)
      jksym=IXOR32(orbsym(no+j),orbsym(k))+1
      ilsym=IXOR32(orbsym(no+i),orbsym(l))+1
_IF1()cj      write(6,*) ' jksym ',jksym,' ilsym ',ilsym
      iadr=(iposjk-1)*norsq(ilsym)+iposil+iofosy(ilsym)
      ibkt1=psobkt(iposjk+iofosq(ilsym))
      ixbkt=(ibkt1-1)/npercl+1
      ibit=mod(ibkt1-1,npercl)
      iofst=ibit*newlen
      iadrx=iadr-apbkt(ibkt1)
      iadrn=iadrx+iofst
c
c    write out value in ixbkt core load in ibkt buffer
c   address in core load = iadrn
c   address in buffer = iadrx
c    original address = iadr
c
         noff=(ixbkt-1)*nszbf
         ifill=ibkt(intowp(noff)+2)+1
c        write (jout,*) ' ifill',ifill
         if (ifill.gt.maxval) then
c           write (jout,*) ' bucket dumped'
c           write (jout,*) ' ifill',ifill
c           write (jout,*) ' nszbf',nszbf
            call bktdmx(bkt(noff+1),nszbf,itap,ichan)
            call vclr(bkt(noff+1),1,nszbf)
            mchain(ixbkt)=ichan
            ibkt(intowp(noff)+1)=ichan
            ifill=1
         end if
         ibkt(intowp(noff)+2+ifill)=iadrn
         bkt(ivoff+noff+ifill)=rint
         ibkt(intowp(noff)+2)=ifill
      if(odebug(33))write(6,*) ' jk/il or kj/li ',jk,il,iadrn,ifill
     1    ,indxol(jk),indxol(il)
      end if
c
c
  101 continue
      if (iflg.eq.0) goto 111
      if(odebug(33))write(6,*) ' about to write to itape 91 '
      do 22 ibkt1=1,nxbkt
         noff=(ibkt1-1)*nszbf
_IF1()cj       ibkt(intowp(noff)+1)=-1
         call bktdmx(bkt(noff+1),nszbf,itap,ichan)
         mchain(ibkt1)=ichan
   22 continue
c
      return
      end
      subroutine tit_mp2(ioff,itriv,isqoo,no,npr,flov,nirred,ntriv,nbf,
     1  buf,ifd2,e,it69,ndimt2,escf)
      implicit integer(a-z)
INCLUDE(common/prnprn)
      integer flov(nirred,4),npr(nirred,3,2),ifd2(nirred),
     1  ioff(nbf),itriv(ntriv),isqoo(no*no)
      REAL e(nbf),buf(2),d,emp2,empx,escf
      data a2 /2.0d0/
c
      if(odebug(33))write(6,*) ' eigval ',e
      emp2=0.0d0
      off2=npr(1,1,2)
      if(odebug(33))write(6,*) ' off2 ',off2
      do 100 asym=1,nirred
      fa=flov(asym,3)-no
      la=flov(asym,4)-no
      do 200 a=fa,la
      do 300 b=fa,a-1
      ab=ioff(a)+b
      aboff=(itriv(ab)-1)*off2
      do 400 isym=1,nirred
      fi=flov(isym,1)
      li=flov(isym,2)
      do 500 i=fi,li
      do 600 j=fi,li
      ji=(j-1)*no+i
      ij=(i-1)*no+j
      ijadd=isqoo(ij)
      jiadd=isqoo(ji)
      ijab=aboff+ijadd
      jiab=aboff+jiadd
      d=e(i)+e(j)-e(a+no)-e(b+no)
      empx=(a2*buf(ijab)*buf(ijab)+a2*buf(jiab)*buf(jiab)
     1   -a2*buf(ijab)*buf(jiab))/d
      if(odebug(33))write(6,*) ' a,b,i,j,empx,d ',a,b,i,j,empx,d
      emp2=emp2+(a2*buf(ijab)*buf(ijab)+a2*buf(jiab)*buf(jiab)
     1   -a2*buf(ijab)*buf(jiab))/d
600   continue
500   continue
400   continue
300   continue
      aa=ioff(a)+a
      aaoff=(itriv(aa)-1)*off2
      do 700 isym=1,nirred
      fi=flov(isym,1)
      li=flov(isym,2)
      do 800 i=fi,li
      do 900 j=fi,li
      ij=(i-1)*no+j
      ijadd=isqoo(ij)
      ijaa=aaoff+ijadd
      d=e(i)+e(j)-e(a+no)-e(a+no)
      empx=(buf(ijaa)*buf(ijaa))/d
      if(odebug(33))write(6,*) ' a,i,j,empx,d ',a,i,j,empx,d
      emp2=emp2+(buf(ijaa)*buf(ijaa))/d
900   continue
800   continue
700   continue
200   continue
100   continue
c
c      write(6,*) ' emp2 so far ',emp2
      do 10 absym=2,nirred
      off2=npr(absym,1,2)
      if(odebug(33))write(6,*) ' off2,absym ',off2,absym
      do 20 asym=1,nirred
      bsym=IXOR32(absym-1,asym-1)+1
      if(asym.gt.bsym) then
      fa=flov(asym,3)-no
      la=flov(asym,4)-no
      fb=flov(bsym,3)-no
      lb=flov(bsym,4)-no
      if(odebug(33))write(6,*) ' asym,bsym,fa,la,fb,lb ',
     +              asym,bsym,fa,la,fb,lb
      do 30 a=fa,la
      do 40 b=fb,lb
      ab=ioff(a)+b
      aboff=(itriv(ab)-1)*off2+ifd2(absym)
      if(odebug(33))write(6,*) ' a,b,aboff,ifd2 ',a,b,aboff,ifd2(absym)
      do 50 isym=1,nirred
      jsym=IXOR32(isym-1,absym-1)+1
      fi=flov(isym,1)
      li=flov(isym,2)
      fj=flov(jsym,1)
      lj=flov(jsym,2)
      do 60 i=fi,li
      do 70 j=fj,lj
      ij=(i-1)*no+j
      ji=(j-1)*no+i
      ijadd=isqoo(ij)
      jiadd=isqoo(ji)
      if(odebug(33))write(6,*) ' i,j,ij,ijadd ',i,j,ij,ijadd
      ijab=aboff+ijadd
      jiab=aboff+jiadd
      d=e(i)+e(j)-e(a+no)-e(b+no)
      empx=(a2*buf(ijab)*buf(ijab)+a2*buf(jiab)*buf(jiab)
     1   -a2*buf(ijab)*buf(jiab))/d
      if(odebug(33))write(6,*) ' a,b,i,j,empx,d ',a,b,i,j,empx,d
      emp2=emp2+(a2*buf(ijab)*buf(ijab)+a2*buf(jiab)*buf(jiab)
     1   -a2*buf(ijab)*buf(jiab))/d
70    continue
60    continue
50    continue
40    continue
30    continue
      end if
20    continue
10    continue
c
c   make t2
c
      icnt=0
      do 1001 asym=1,nirred
      fa=flov(asym,3)-no
      la=flov(asym,4)-no
      do 2001 a=fa,la
      do 3001 b=fa,a
      do 4001 isym=1,nirred
      fi=flov(isym,1)
      li=flov(isym,2)
      do 5001 i=fi,li
      do 6001 j=fi,li
      d=e(i)+e(j)-e(a+no)-e(b+no)
      if(odebug(33))write(6,*) ' i,j,a,b,ijab,buf ',
     +                           i,j,a,b,ijab,buf(icnt)
      icnt=icnt+1
      buf(icnt)=buf(icnt)/d
6001  continue
5001  continue
4001  continue
3001  continue
2001  continue
1001  continue
      do 101 absym=2,nirred
      do 201 asym=1,nirred
      bsym=IXOR32(absym-1,asym-1)+1
      if(asym.gt.bsym) then
      fa=flov(asym,3)-no
      la=flov(asym,4)-no
      fb=flov(bsym,3)-no
      lb=flov(bsym,4)-no
      do 301 a=fa,la
      do 401 b=fb,lb
      do 501 isym=1,nirred
      jsym=IXOR32(isym-1,absym-1)+1
      fi=flov(isym,1)
      li=flov(isym,2)
      fj=flov(jsym,1)
      lj=flov(jsym,2)
      do 601 i=fi,li
      do 701 j=fj,lj
      d=e(i)+e(j)-e(a+no)-e(b+no)
      icnt=icnt+1
      buf(icnt)=buf(icnt)/d
701   continue
601   continue
501   continue
401   continue
301   continue
      end if
201   continue
101   continue
c
      write(6,342)emp2,emp2+escf
 342  format(' mp2 correlation energy ',e20.12,/,
     &       ' scf+mp2 energy         ',e20.12)
      call swrit(it69,buf,intowp(ndimt2))
      return
      end
      subroutine tsort(core,icore,lscr,nsymx,noccx,nvirx,potnucx,itype)
      implicit integer(a-z)
      REAL core(lscr)
      integer icore(2*lscr)
      integer nsymx,noccx(nsymx),nvirx(nsymx),itype(*)
      REAL potnucx,cpulft
c
      logical irew
      REAL enuc,escf      
      integer wpadti,sec2i
INCLUDE(common/prnprn)
INCLUDE(common/titan)
      character*10 charwall
c
      maxcor = lscr
cjvl      maxbkt=3000
      maxbkt=18000
      nlist=9
      in=5
      iw=6
_IF1()    1 format ('                                                ')
_IF1()    2 format ('************************************************')
c
      write(iw,99)cpulft(1) ,charwall()
 99   format(/1x,'commence titan ccsd integral sort at', f9.2,
     +           ' seconds',a10,' wall')
c
      irew=.true.
      input=5
      jout=6
      lu1el=22
      lu2el=22
      itap30=30
      itap60=60
      itap61=61
      itap62=62
      itap63=63
      itap64=64
      itap65=65
      itap66=66
      itap67=67
      itap69=69
      itap78=78
      itap91=91
      itap57=57
      intbuf=sec2i(5)/intowp(1)
c
      call rfile(itap30)
      call rfile(itap60)
      call rfile(itap61)
      call rfile(itap62)
      call rfile(itap63)
      call rfile(itap64)
      call rfile(itap65)
      call rfile(itap66)
_IF1()cj    call rfile(itap67)
      call rfile(itap69)
      call rfile(itap78)
      call rfile(itap91)
c
      ltyp=7
cjvl      niobf=200
cjvl      niobf=225
      niobf=255
c      write(6,*) ' niobf  = ',niobf
      length=niobf*(sec2i(1)/intowp(1))
      nszbf=niobf*(sec2i(1)/intowp(1))
      maxval = (intowp(length)-2)/(1+intowp(1))
      ibktsp=length*ltyp
      ivoff = (maxval+2)/intowp(1)
      if(intowp(ivoff).ne.(maxval+2)) ivoff = ivoff + 1
_IF1()c
_IF1()c  read in input for ccsd calculation
_IF1()c
_IF1()c      open(unit=5,file='input.data')
_IF1()c --  read in specifications for scratch directory on cfs
_IF1()c
_IF1()c      call daname(lu1el,'ORD1EL',0)
_IF1()c      call rwtp(lu1el,ad1el)
_IF1()c      call rlist(lu1el,ad1el,4,nsym,1,nocc,8,nvir,8,potnuc,2)
      nsym=nsymx
      do i=1,nsym
       nocc(i)=noccx(i)
       nvir(i)=nvirx(i)
      enddo
      potnuc=potnucx
c
      enuc=potnuc
      nirred=nsym
      write(6,845)nirred
 845  format(/,' number of symmetries ',i5)
      write(6,846)(nocc(irr),irr=1,nirred)
 846  format(' occupied ',8i5)
      write(6,847)(nvir(irr),irr=1,nirred)
 847  format(' virtual  ',8i5)
      doc = 1
      uoc = doc + nirred
      nbfo = 0
      noo = 0
      nvo = 0
      do 1012 irr = 1,nirred
c      read(5,*) icore(doc+irr-1),icore(uoc+irr-1)
      icore(doc+irr-1)=nocc(irr)
      icore(uoc+irr-1)=nvir(irr)
      noo = noo + icore(doc+irr-1)
      nvo = nvo + icore(uoc+irr-1) 
 1012 continue
      nbfo = noo + nvo
      nbf = nbfo
      no = noo
      nv = nvo
c
      ntris = 0
      nsqs = 0
      do 1112 irr = 1,nirred
      nmoir = icore(doc+irr-1) + icore(uoc+irr-1)
      ntris = ntris + (nmoir+1)*nmoir/2
      do 1113 jrr = 1,irr-1
      nmojr = icore(doc+jrr-1) + icore(uoc+jrr-1)
      nsqs = nsqs + nmoir*nmojr
 1113 continue
 1112 continue
c      print *,'ntris ',ntris
c
      idu = uoc + nirred
      ktri = idu + nbf
      ltri = ktri + ntris
      ksq = ltri + ntris
      lsq = ksq + nsqs
      ptsq = lsq + nsqs
      bkof1 = ptsq + nirred*nirred
      bkof2 = bkof1 + 16
      tfile = bkof2 + 16
      ibeg = tfile + 16 
c
      ntri=(nbfo*(nbfo+1))/2
c
      ic=0
      flov = ibeg
      imap = flov + nirred*4
      ptocc = imap + nbfo
      oorbsy = ptocc + nbfo
      orbsym = oorbsy + nbfo
      eigval = iadtwp(orbsym+nbfo)
      ik3 = wpadti(eigval+nbfo)
c
c  if using tjl's codes, read file30 for orbsym & eigval and 
c  construct flov & ptocc - also reorder if cor's and vir's present.
c  if using nasa codes then form orbsym, flov, ptocc (which is now
c  nasa to cc ordering).
c
      call nasai(icore(doc),icore(uoc),icore(flov),icore(orbsym),
     1           icore(ptocc),nirred,nbf,no,nv,icore(idu),icore(ktri),
     2           icore(ltri),icore(ksq),icore(lsq),icore(ptsq),
     3           icore(bkof1),icore(bkof2),icore(tfile),ntris,nsqs,
     4           length,ivoff,itype)
c
      ntriv=nv*(nv+1)/2
      ntrio=no*(no+1)/2
_IF1()ctjl      write(6,*) ' ntriv,nv,no,ntrio ',ntriv,nv,no,ntrio
_IF1()ctjl      write(*,*)  '  nv = ',nv
      npr=ik3
      ioff=npr+nirred*6
c
      isqoo=ioff+nbfo
      isqvv=isqoo+no*no
      itriv=isqvv+nv*nv
      itrio=itriv+ntriv
      isqov=itrio+ntrio
      isqvo=isqov+no*nv
      ifa=isqvo+no*nv
      ifa2=ifa+nirred+1
      ifb=ifa2+nirred+1
      ifc=ifb+nirred+1
      ifd1=ifc+nirred+1
      ifd2=ifd1+nirred+1
      ifd3=ifd2+nirred+1
      ifd4=ifd3+nirred+1
      ife1=ifd4+nirred+1
      iff=ife1+nirred+1
      fpbka=iff+nirred+1
      fpbkb=fpbka+maxbkt
      fpbkc=fpbkb+maxbkt
      fpbkd1=fpbkc+maxbkt
      fpbkd2=fpbkd1+maxbkt
      fpbkd3=fpbkd2+maxbkt
      fpbke1=fpbkd3+maxbkt
      fpbkf=fpbke1+maxbkt
      fsec=fpbkf+maxbkt
      do 1018 jki=1,nlist
1018  icore(fsec+jki-1)=1
      apbkt=fsec+nlist
      ik3=apbkt+maxbkt
c     write (jout,1)
c     write (jout,2)
c     write (jout,*) ' ************  nlamda  *************'
c     write (jout,2)
c     write (jout,1)
c     write (jout,*) icore (421)
c
      rtop=iadtwp(ik3)
      call vclr(core(rtop),1,maxcor-rtop)
c
c      lnbuf = 14336
      lnBuf=128*256
c      lnbuf = 128
      hone = iadtwp(ik3)
      ni = wpadti(hone+ntri)
      nj = ni + lnbuf
      nk = nj + lnbuf
      nl = nk + lnbuf
      fone = iadtwp(nl+lnbuf)
      buf = fone + ntri
      ibuf = wpadti(buf)
      buf2 = buf + lnbuf
      ibuf2 = wpadti(buf2)
      bkt = buf2 + lnbuf
      ibkt = wpadti(bkt)
      ttop = bkt + ibktsp
c      write(*,*) ttop,'  real words required for initial sort '
      if(ttop.gt.maxcor) then
       write(jout,*) '  not enough memory before nasa sort '
       write(jout,*) '  ttop, maxcor = ',ttop,maxcor
       call caserr('insufficient memory for integral t-sort')
      end if
c
c  core allocated - call rdnai
c
c      print *,'ntris ',ntris
      call rdnai(core(buf),icore(ibuf),core(buf2),icore(ibuf2),
     1          core(bkt),icore(ibkt),icore(ni),icore(nj),icore(nk),
     2          icore(nl),core(hone),core(fone),icore(doc),icore(uoc),
     3          icore(ptocc),icore(idu),icore(ktri),icore(ltri),
     4          icore(ksq),icore(lsq),icore(ptsq),icore(bkof1),
     5          icore(bkof2),icore(tfile),jout,nirred,nbf,no,nv,ntri,
     6          ntris,nsqs,lnbuf,length,ibktsp,maxval,enuc,
     &          icore(orbsym),itap30,itap60,itap61,itap62,itap63,
     &          itap64,
     &          itap65,itap66,itap67,itap69,itap78,itap91,itap57)
c  read in nasa integrals - now evaluate the orbital eigenvalues and
c  the scf energy.
c
      buf = hone + ntri
      ibuf = wpadti(buf)
      ni = wpadti(buf+length)
      nj = ni + maxval
      nk = nj + maxval
      nl = nk + maxval
      ttop = nl + maxval
c      write(*,*) '   ttop 2nd time around = ',ttop
c
      call fmatx(core(buf),icore(ibuf),length,ivoff,core(hone),ntri,
     1           core(eigval),nbf,escf,maxval,itap60,itap62,itap63,
     2           icore(ni),icore(nj),icore(nk),icore(nl),no,enuc)
c
      call rclose(itap30,3)
      call rclose(itap78,3)
_IF1()cj    call rclose(itap61,3)
_IF1()cj    call rclose(itap62,3)
_IF1()cj    call rclose(itap63,3)
_IF1()cj    call rclose(itap64,3)
_IF1()cj    call rclose(itap65,3)
_IF1()cj    call rclose(itap67,3)
_IF1()cj      call rclose(itap66,3)
c      write (jout,*) ' called rclose tape 78'
      nprsq=ik3
      nprb=nprsq+nirred
      iofpr=nprb+nirred
      norbv=iofpr+nirred+1
      norbo=norbv+nirred
      nprgs=norbo+nirred
      nprgsq=nprgs+nirred*nirred
      itop=nprgsq+nirred*nirred
      buf=iadtwp(itop)
      ibuf=wpadti(buf)
      t2=buf+length
      newlen=niobf*(sec2i(1)/intowp(1))
      if(nv*nv.gt.newlen) then
        write(6,*)
     +  'not enough buffer space to write 1 block of integrals'
        write(6,*) ' buffer space = ',newlen,' nv*nv = ',nv*nv
        write(6,*) ' possible cause is niobf,which is set to ',niobf
      call caserr('insufficient buffer space for integral sort')
      end if
      space1=maxcor-length-iadtwp(itop)+1
      call sett2(icore(flov),icore(norbv)
     1 ,icore(norbo),icore(orbsym),icore(nprb),icore(nprsq),
     1  icore(ifd2),icore(isqoo),
     1   icore(ioff),icore(iofpr),nv,ntriv,nbf,nirred,space1,
     1   icore(itriv),
     1 icore(nprgs),icore(nprgsq),ndimt1,ndimt2,icore(npr),no)
c
      call initt2(core(buf),icore(ibuf),nszbf,length,nbf,
     .     itap63,jout,icore(itriv),ntriv,core(eigval),
     1 icore(norbo),icore(ifd2),icore(nprb),nirred,no,icore(isqoo),
     1 maxval,icore(orbsym),core(t2),itap69,ndimt1,ndimt2,icore(ioff))
c
      call tit_mp2(icore(ioff),icore(itriv)
     1 ,icore(isqoo),no,icore(npr),icore(flov),nirred,ntriv,nbf,
     1  core(t2),icore(ifd2),core(eigval),itap69,ndimt2,escf)
c
      call rclose(itap69,3)
      ijer=1
_IF1()cj      if(ijer.eq.1) stop
      nprsq=ik3
      nprb=nprsq+nv*nv
      prbkt=nprb+nirred
      mchain=prbkt+ntriv
      iofpr=mchain+maxbkt
      norbt=iofpr+nirred+1
      nprgs=norbt+nirred
      nprgsq=nprgs+nirred*nirred
      itop=nprgsq+nirred*nirred
      space=maxcor-iadtwp(itop)+1
c
c    space1 is number of values that can be held in newlen
c
      space1=newlen
_IF1()cj    intlen=(intowp(nszbf)-2)/intowp(1)
_IF1()cj    space1=intowp(intlen)/(1+intowp(1))
c
c    no of newlen buffers in each core load in sort
c
      npercl=(space-nszbf-1)/newlen
c      write(6,*) ' npercl ',npercl
c
      if(odebug(33)) then
       write(6,*) ' flov bf setup ',(icore(flov+jki-1),jki=1,16)
        write(6,*) ' flov = ',flov
      endif
      call setup(icore(flov),icore(norbt),icore(orbsym),icore(nprb),
     1   icore(nprsq),icore(ifb),icore(itriv),icore(isqvv),icore(ioff),
     2   icore(iofpr),icore(fpbkb),icore(prbkt),nv,ntriv,nbf,nirred,
     3maxbkt,space1,nbktb,icore(nprgs),icore(nprgsq),icore(npr),no,
     4icore(apbkt))
c
c    nbkt = total no of newlen buffers needed
c    nxbkt = total no of core loads required
c
      nxbkt=(nbktb-1)/npercl+1
c      write(6,*) ' nxbkt ',nxbkt
      lpercl=mod(nbktb,npercl)
      if(lpercl.eq.0) lpercl=npercl
c      write(6,*) ' lpercl,nbktb ',lpercl,nbktb
c
      ntotsz=nxbkt*nszbf+length+iadtwp(itop)
c      write(*,*) ' ntotsz,maxcor,itop = ',ntotsz,maxcor,iadtwp(itop)
      if(ntotsz.gt.maxcor) then
        write(6,*) ' reducing sort file buffer size '
        nszbf = (maxcor - length - iadtwp(itop))/nxbkt
        nszbf = nszbf/(sec2i(1)/intowp(1))
        if(nszbf.lt.1) then
        nreq = nxbkt*sec2i(1)/intowp(1) + length + iadtwp(itop)
        write(6,*) ' not enough space for buffers in sort '
        write(6,*) ' require ',nreq,'  real words of memory'
        write(6,*) ' have    ',maxcor,'  real words of memory'
        call mabort
        end if
        nszbf = nszbf*(sec2i(1)/intowp(1))
        write(6,*) ' new buffer size = ',nszbf
      end if
c
      intbuf=iadtwp(itop)
      jintbf=wpadti(intbuf)
      bufx=intbuf+length
      ibufx=wpadti(bufx)
      ktop=bufx+nxbkt*nszbf
      if(ktop.gt.maxcor) then
      write(*,*)  '   not enough core :  ktop,maxcor = ',ktop,maxcor
      call mabort
      end if
c
      call vclr(core(bufx),1,nszbf*nxbkt)
      call sortj(core(bufx),icore(ibufx),core(intbuf),icore(jintbf),
     1    nszbf,nxbkt,length,nbf,
     1npercl,itap61,jout,itap91,icore(ioff),icore(itriv),icore(isqvv),
     2    newlen,icore(norbt),icore(ifb),icore(nprsq),icore(prbkt),
     3    icore(iofpr),icore(fpbkb),nirred,nv,ntriv,no,icore(mchain),
     4    maxval,icore(orbsym),icore(apbkt))
c
      sortbf=iadtwp(itop)
      isorbf=wpadti(sortbf)
      rints=sortbf+nszbf
      ktop2=rints+newlen*npercl
      if(ktop2.gt.maxcor) then
      write(*,*)  '   not enough core :  ktop2,maxcor = ',ktop2,maxcor
      call mabort
      end if
c
c
      call bckchn(icore(mchain),itap61,core(sortbf),icore(isorbf),
     1     nszbf,maxbkt,maxval,core(rints),jout,nxbkt,newlen,npercl,
     1     itap91,lpercl,irew)
c
      call rclose(itap61,3)
      call rclose(itap91,4)
      call rfile(itap91)
c
      ijer = 1
_IF1()cj      if(ijer.eq.1) stop
      norbo=ik3
      norsq=norbo+nirred
      iofosq=norsq+nirred
      mchain=iofosq+nirred+1
      psobkt=mchain+maxbkt
      itop=psobkt+no*no
      space=maxcor-iadtwp(itop)+1
      npercl=(space-nszbf-1)/newlen
      if(odebug(33))
     +   write(6,*) ' no,nbf,nirred before setupo ',no,nbf,nirred
      call setupo(icore(flov),icore(norbo),icore(orbsym),icore(norsq),
     1    icore(ifa),
     1   icore(isqoo),icore(iofosq),icore(fpbka),
     1  icore(psobkt),no,nbf,nirred,maxbkt,space1,nbkta,
     2   icore(npr),icore(apbkt))
c
      nxbkt=(nbkta-1)/npercl +1
c      write(6,*) ' nxbkt for sorto ',nxbkt
      lpercl=mod(nbkta,npercl)
      if(lpercl.eq.0) lpercl=npercl
c      write(6,*) ' lpercl,nbkta ',lpercl,nbkta
      ntotsz=nxbkt*nszbf+length+iadtwp(itop)
      if(ntotsz.gt.maxcor) then
         write(6,*) ' not enough space for buffers in sorto '
         write(6,*) ' require ',ntotsz,' words; have ',maxcor,
     1    ' words '
         call caserr('insufficient space for buffers')
      end if
      intbuf=iadtwp(itop)
      jintbf=wpadti(intbuf)
       bufx=intbuf+length
      ibufx=wpadti(bufx)
      ktop=bufx+nxbkt*nszbf
      if(ktop.gt.maxcor) then
       write(*,*)  '   not enough core :  ktop,maxcor = ',ktop,maxcor
      call mabort
       end if
c
       call vclr(core(bufx),1,nszbf*nxbkt)
      call sorto(core(bufx),icore(ibufx),core(intbuf),icore(jintbf),
     1    nszbf,nxbkt,length,nbf,
     1npercl,itap60,jout,itap91,icore(isqoo),
     2  newlen,icore(norbo),icore(ifa),icore(norsq),icore(psobkt),
     3    icore(iofosq),icore(fpbka),nirred,no,icore(mchain),
     4    maxval,icore(orbsym),icore(apbkt))
c
      sortbf=iadtwp(itop)
      isorbf=wpadti(sortbf)
      rints=sortbf+nszbf
      ktop2=rints+newlen*npercl
      if(ktop2.gt.maxcor) then
      write(*,*)  '   not enough core :  ktop2,maxcor = ',ktop2,maxcor
      call mabort
      end if
c
c
      call bckchn(icore(mchain),itap60,core(sortbf),icore(isorbf),
     1     nszbf,maxbkt,maxval,core(rints),jout,nxbkt,newlen,npercl,
     1     itap91,lpercl,irew)
c
       call rclose(itap60,3)
      call rclose(itap91,4)
c      write(*,*) ' called rclose for file 91 '
      call rfile(itap91)
      norbo=ik3
      norsq=norbo+nirred
      iofosq=norsq+nirred
      mchain=iofosq+nirred+1
      psobkt=mchain+maxbkt
      norgsq=psobkt+no*nv
      norbv=norgsq+nirred*nirred
      itop=norbv+nirred
      space=maxcor-iadtwp(itop)+1
      npercl=(space-nszbf-1)/newlen
c      write(*,*) ' calling setvo1 '
      call setvo1(icore(flov),icore(norbo),icore(orbsym),icore(norsq),
     1    icore(ifc),
     1   icore(isqov),icore(iofosq),icore(fpbkc),
     1  icore(psobkt),no,nbf,nirred,maxbkt,space1,nbktc,icore(norgsq),
     2   icore(npr),nv,icore(norbv),icore(apbkt))
      nxbkt=(nbktc-1)/npercl +1
      nbktd2=nbktc
      nbktd1=nbktc
      lpercl=mod(nbktc,npercl)
      if(lpercl.eq.0) lpercl=npercl
c      write(6,*) ' lpercl,nbktc ',lpercl,nbktc
      do 114 jki=1,nirred
114   icore(ifd1+jki-1)=icore(ifc+jki-1)
      do 109 jki=1,nbktd2+1
      icore(fpbkd1+jki-1)=icore(fpbkc+jki-1)
109   icore(fpbkd2+jki-1)=icore(fpbkc+jki-1)
c      write(6,*) ' nxbkt for sorto ',nxbkt
      ntotsz=nxbkt*nszbf+length+iadtwp(itop)
      if(ntotsz.gt.maxcor) then
         write(6,*) ' not enough space for buffers in sorto '
      write(6,*) ' require ',ntotsz,' words; have ',maxcor,
     1    ' words '
         call mabort
      end if
      intbuf=iadtwp(itop)
      jintbf=wpadti(intbuf)
       bufx=intbuf+length
      ibufx=wpadti(bufx)
      ktop=bufx+nxbkt*nszbf
      if(ktop.gt.maxcor) then
       write(*,*)  '   not enough core :  ktop,maxcor = ',ktop,maxcor
      call mabort
       end if
c
       call vclr(core(bufx),1,nszbf*nxbkt)
      call sorvo1(core(bufx),icore(ibufx),core(intbuf),icore(jintbf),
     1    nszbf,nxbkt,length,nbf,
     1npercl,itap62,jout,itap91,icore(isqov),
     2  newlen,icore(norbo),icore(ifc),icore(norsq),icore(psobkt),
     3    icore(iofosq),icore(fpbkc),nirred,no,icore(mchain),
     4    maxval,icore(orbsym),icore(apbkt))
c
      sortbf=iadtwp(itop)
      isorbf=wpadti(sortbf)
      rints=sortbf+nszbf
      ktop2=rints+newlen*npercl
      if(ktop2.gt.maxcor) then
      write(*,*)  '   not enough core :  ktop2,maxcor = ',ktop2,maxcor
      call mabort
      end if
c
c
      call bckchn(icore(mchain),itap62,core(sortbf),icore(isorbf),
     1     nszbf,maxbkt,maxval,core(rints),jout,nxbkt,newlen,npercl,
     1     itap91,lpercl,irew)
c
       call rclose(itap62,3)
      call rclose(itap91,4)
      call rfile(itap91)
c
      intbuf=iadtwp(itop)
      jintbf=wpadti(intbuf)
       bufx=intbuf+length
      ibufx=wpadti(bufx)
      ktop=bufx+nxbkt*nszbf
      if(ktop.gt.maxcor) then
       write(*,*)  '   not enough core :  ktop,maxcor = ',ktop,maxcor
      call mabort
       end if
c
       call vclr(core(bufx),1,nszbf*nxbkt)
      call sord1(core(bufx),icore(ibufx),core(intbuf),icore(jintbf),
     1    nszbf,nxbkt,length,nbf,
     1npercl,itap63,jout,itap91,icore(isqov),
     2  newlen,icore(norbo),icore(ifd1),icore(norsq),icore(psobkt),
     3    icore(iofosq),icore(fpbkd1),nirred,no,icore(mchain),
     4    maxval,icore(orbsym),nv,icore(apbkt))
c
      sortbf=iadtwp(itop)
      isorbf=wpadti(sortbf)
      rints=sortbf+nszbf
      ktop2=rints+newlen*npercl
      if(ktop2.gt.maxcor) then
      write(*,*)  '   not enough core :  ktop2,maxcor = ',ktop2,maxcor
      call mabort
      end if
c
c
      call bckchn(icore(mchain),itap66,core(sortbf),icore(isorbf),
     1     nszbf,maxbkt,maxval,core(rints),jout,nxbkt,newlen,npercl,
     1     itap91,lpercl,irew)
c
      call rgetsa(itap66,icore(fsec+4))
c      write(6,*) ' fsec(5) ' ,icore(fsec+4)
      call rclose(itap91,4)
c
      call rfile(itap91)
      intbuf=iadtwp(itop)
      jintbf=wpadti(intbuf)
       bufx=intbuf+length
      ibufx=wpadti(bufx)
      ktop=bufx+nxbkt*nszbf
      if(ktop.gt.maxcor) then
       write(*,*)  '   not enough core :  ktop,maxcor = ',ktop,maxcor
      call mabort
       end if
c
       call vclr(core(bufx),1,nszbf*nxbkt)
      call sord2(core(bufx),icore(ibufx),core(intbuf),icore(jintbf),
     1    nszbf,nxbkt,length,nbf,
     1npercl,itap63,jout,itap91,icore(isqov),
     2  newlen,icore(norbo),icore(ifc),icore(norsq),icore(psobkt),
     3    icore(iofosq),icore(fpbkd2),nirred,no,icore(mchain),
     4    maxval,icore(orbsym),icore(apbkt))
c
      sortbf=iadtwp(itop)
      isorbf=wpadti(sortbf)
      rints=sortbf+nszbf
      ktop2=rints+newlen*npercl
      if(ktop2.gt.maxcor) then
      write(*,*)  '   not enough core :  ktop2,maxcor = ',ktop2,maxcor
      call mabort
      end if
c
c
      irew=.false.
      call bckchn(icore(mchain),itap66,core(sortbf),icore(isorbf),
     1     nszbf,maxbkt,maxval,core(rints),jout,nxbkt,newlen,npercl,
     1     itap91,lpercl,irew)
c
      call rgetsa(itap66,icore(fsec+5))
c      write(6,*) ' fsec(6) ',icore(fsec+5)
      call rclose(itap91,4)
c
      call rfile(itap91)
      nprsq=ik3
      nprb=nprsq+nv*nv
      prbkt=nprb+nirred
      mchain=prbkt+ntrio
      iofpr=mchain+maxbkt
      norbt=iofpr+nirred+1
      norbv=norbt+nirred
      nprgs=norbv+nirred
      itop=nprgs+nirred*nirred
      newlen=niobf*(sec2i(1)/intowp(1))
c      write(6,*) ' newlen = ',newlen
      nszbf=niobf*(sec2i(1)/intowp(1))
      space=maxcor-iadtwp(itop)+1
c
c    space1 is number of values that can be held in newlen
c
      space1=newlen
_IF1()cj    intlen=(intowp(nszbf)-2)/intowp(1)
_IF1()cj    space1=intowp(intlen)/(1+intowp(1))
c
c    no of newlen buffers in each core load in sort
c
      npercl=(space-nszbf-1)/newlen
c      write(6,*) ' npercl ',npercl
c
      if(odebug(33)) then
       write(6,*) ' flov bf setup ',(icore(flov+jki-1),jki=1,16)
       write(6,*) ' flov = ',flov
      endif
      call setd3(icore(flov),icore(norbt),icore(orbsym),icore(nprb),
     1   icore(nprsq),icore(ifd3),icore(itrio),icore(isqvv),icore(ioff),
     2   icore(iofpr),icore(fpbkd3),icore(prbkt),nv,ntriv,nbf,nirred,
     3    maxbkt,space1,nbktd3,icore(nprgs),icore(npr),
     4    icore(norbv),ntrio,icore(apbkt))
c
c    nbkt = total no of newlen buffers needed
c    nxbkt = total no of core loads required
c
      nxbkt=(nbktd3-1)/npercl+1
      lpercl=mod(nbktd3,npercl)
      if(lpercl.eq.0) lpercl=npercl
c      write(6,*) ' nxbkt ',nxbkt
c
      ntotsz=nxbkt*nszbf+length+iadtwp(itop)
      if(ntotsz.gt.maxcor) then
         write(6,*) ' not enough space for buffers in sort '
         write(6,*) ' require ',ntotsz,' words; have ',maxcor,
     1    ' words '
         call caserr('insufficient space for buffers in sort')
      end if
c
      intbuf=iadtwp(itop)
      jintbf=wpadti(intbuf)
      bufx=intbuf+length
      ibufx=wpadti(bufx)
      ktop=bufx+nxbkt*nszbf
      if(ktop.gt.maxcor) then
      write(*,*)  '   not enough core :  ktop,maxcor = ',ktop,maxcor
      call mabort
      end if
c
      call vclr(core(bufx),1,nszbf*nxbkt)
      call sortd3(core(bufx),icore(ibufx),core(intbuf),icore(jintbf),
     1    nszbf,nxbkt,length,nbf,
     1npercl,itap63,jout,itap91,icore(ioff),icore(itrio),icore(isqvv),
     2    newlen,icore(norbt),icore(ifd3),icore(nprsq),icore(prbkt),
     3    icore(iofpr),icore(fpbkd3),nirred,nv,ntriv,no,icore(mchain),
     4    maxval,icore(orbsym),icore(apbkt))
c
      sortbf=iadtwp(itop)
      isorbf=wpadti(sortbf)
      rints=sortbf+nszbf
      ktop2=rints+newlen*npercl
      if(ktop2.gt.maxcor) then
      write(*,*)  '   not enough core :  ktop2,maxcor = ',ktop2,maxcor
      call mabort
      end if
c
c
      call bckchn(icore(mchain),itap66,core(sortbf),icore(isorbf),
     1     nszbf,maxbkt,maxval,core(rints),jout,nxbkt,newlen,npercl,
     1     itap91,lpercl,irew)
c
      call rclose(itap63,3)
      call rclose(itap66,3)
      call rclose(itap91,4)
c
      call rfile(itap91)
      norbo=ik3
      norsq=norbo+nirred
      iofosq=norsq+nirred
      mchain=iofosq+nirred+1
      psobkt=mchain+maxbkt
      norgsq=psobkt+no*nv
      norbv=norgsq+nirred*nirred
      nvrsq=norbv+nirred
      itop=nvrsq+nirred
      space=maxcor-iadtwp(itop)+1
      npercl=(space-nszbf-1)/newlen
      call sete1(icore(flov),icore(norbo),icore(orbsym),icore(norsq),
     1    icore(ife1),
     1   icore(isqvo),icore(iofosq),icore(fpbke1),
     1  icore(psobkt),no,nbf,nirred,maxbkt,space1,nbkte1,icore(norgsq),
     2icore(npr),nv,icore(norbv),icore(isqoo),icore(nvrsq),icore(apbkt))
      nxbkt=(nbkte1-1)/npercl +1
      lpercl=mod(nbkte1,npercl)
      if(lpercl.eq.0) lpercl=npercl
c      write(6,*) ' lpercl,nbkte1 ',lpercl,nbkte1
c      write(6,*) ' nxbkt for sorte1 ',nxbkt
      ntotsz=nxbkt*nszbf+length+iadtwp(itop)
      if(ntotsz.gt.maxcor) then
         write(6,*) ' not enough space for buffers in sorte1 '
      write(6,*) ' require ',ntotsz,' words; have ',maxcor,
     1    ' words '
         call mabort
      end if
      intbuf=iadtwp(itop)
      jintbf=wpadti(intbuf)
       bufx=intbuf+length
      ibufx=wpadti(bufx)
      ktop=bufx+nxbkt*nszbf
      if(ktop.gt.maxcor) then
       write(*,*)  '   not enough core :  ktop,maxcor = ',ktop,maxcor
      call mabort
       end if
c
       call vclr(core(bufx),1,nszbf*nxbkt)
      call sorte1(core(bufx),icore(ibufx),core(intbuf),icore(jintbf),
     1    nszbf,nxbkt,length,nbf,
     1npercl,itap64,jout,itap91,icore(isqov),
     2  newlen,icore(ife1),icore(nvrsq),icore(psobkt),
     3    icore(iofosq),icore(fpbke1),nirred,no,icore(mchain),
     4    maxval,icore(orbsym),icore(isqoo),icore(apbkt),nv)
c
      sortbf=iadtwp(itop)
      isorbf=wpadti(sortbf)
      rints=sortbf+nszbf
      ktop2=rints+newlen*npercl
      if(ktop2.gt.maxcor) then
      write(*,*)  '   not enough core :  ktop2,maxcor = ',ktop2,maxcor
      call mabort
      end if
c
c
      irew=.true.
      call bckchn(icore(mchain),itap64,core(sortbf),icore(isorbf),
     1     nszbf,maxbkt,maxval,core(rints),jout,nxbkt,newlen,npercl,
     1     itap91,lpercl,irew)
c
      call rclose(itap91,4)
c
      call rclose(itap64,3)
c
      call rfile(itap91)
      norbo=ik3
      norsq=norbo+nirred
      iofosq=norsq+nirred
      mchain=iofosq+nirred+1
      psobkt=mchain+maxbkt
      norgsq=psobkt+nv*no
      norbv=norgsq+nirred*nirred
      nvrsq=norbv+nirred
      itop=nvrsq+nirred
      space=maxcor-iadtwp(itop)+1
      npercl=(space-nszbf-1)/newlen
      call setf(icore(orbsym),icore(norsq),
     1    icore(iff),
     1   icore(iofosq),icore(fpbkf),
     1  icore(psobkt),no,nbf,nirred,maxbkt,space1,nbktf,
     2nv,icore(nvrsq),icore(apbkt))
      nxbkt=(nbktf-1)/npercl +1
      lpercl=mod(nbktf,npercl)
      if(lpercl.eq.0) lpercl=npercl
c      write(6,*) ' lpercl,nbktf ',lpercl,nbktf
c      write(6,*) ' nxbkt for sortf ',nxbkt
      ntotsz=nxbkt*nszbf+length+iadtwp(itop)
      if(ntotsz.gt.maxcor) then
         write(6,*) ' not enough space for buffers in sorte1 '
      write(6,*) ' require ',ntotsz,' words; have ',maxcor,
     1    ' words '
         call mabort
      end if
      intbuf=iadtwp(itop)
      jintbf=wpadti(intbuf)
       bufx=intbuf+length
      ibufx=wpadti(bufx)
      ktop=bufx+nxbkt*nszbf
      if(ktop.gt.maxcor) then
       write(*,*)  '   not enough core :  ktop,maxcor = ',ktop,maxcor
      call mabort
       end if
c
      call vclr(core(bufx),1,nszbf*nxbkt)
      call sortf(core(bufx),icore(ibufx),core(intbuf),icore(jintbf),
     1    nszbf,nxbkt,length,nbf,
     1npercl,itap65,jout,itap91,icore(isqov),
     2  newlen,icore(iff),icore(nvrsq),icore(psobkt),
     3    icore(iofosq),icore(fpbkf),nirred,no,icore(mchain),
     4    maxval,icore(orbsym),icore(isqvv),icore(apbkt),nv)
c
      sortbf=iadtwp(itop)
      isorbf=wpadti(sortbf)
      rints=sortbf+nszbf
      ktop2=rints+newlen*npercl
      if(ktop2.gt.maxcor) then
      write(*,*)  '   not enough core :  ktop2,maxcor = ',ktop2,maxcor
      call mabort
      end if
c
c
      irew=.true.
      call bckchn(icore(mchain),itap65,core(sortbf),icore(isorbf),
     1     nszbf,maxbkt,maxval,core(rints),jout,nxbkt,newlen,npercl,
     1     itap91,lpercl,irew)
c
      call rclose(itap65,3)
      call rclose(itap91,4)
c
      call wrt57(itap57,nirred,nv,no,icore(flov),newlen,nlist,nbkta,
     1 nbktb,nbktc,nbktd1,nbktd2,nbktd3,nbkte1,nbktf,
     1 icore(npr),
     2 icore(fpbka),icore(fpbkb),icore(fpbkc),icore(fpbkd1),
     3 icore(fpbkd2),icore(fpbkd3),icore(fpbke1),
     4 icore(fpbkf),escf,enuc,core(eigval),icore(orbsym),nbf,
     5 ndimt1,ndimt2,icore(isqvo),
     4 icore(fsec),icore(itriv),icore(itrio),icore(isqvv),icore(isqoo),
     4 icore(isqov),icore(ioff),icore(ifa),icore(ifb),icore(ifc),
     5 core(ifd1),icore(ifd2),icore(ifd3),icore(ife1),
     5 icore(iff),icore(ptocc),ntriv,ntrio,icore(ifd4),icore(ifa2))
      write(iw,999)cpulft(1) ,charwall()
 999  format(/' end of integral sort at ',f9.2,' seconds',a10,' wall')
c
c
c      call tstop(6)
c
      return
      end
      integer FUNCTION WPADTI(N)
_IF(cray,i8)
      WPADTI=N
_ELSE
      WPADTI=N*2-1
_ENDIF
      RETURN
      END
      subroutine wrt57(itap57,nirred,nv,no,flov,lnbkt,nlist,nbka,
     1 nbkb,nbkc,nbkd1,nbkd2,nbkd3,nbke1,nbkf,npr,
     2    fpbka,fpbkb,fpbkc,fpbkd1,fpbkd2,fpbkd3,fpbke1,
     3    fpbkf,escf,enuc,eigval,orbsym,nbf,ndimt1,ndimt2,isqvo,
     4    fsec,itriv,itrio,isqvv,isqoo,isqov,ioff,ifa,ifb,ifc,
     5ifd1,ifd2,ifd3,ife1,iff,ptocc,ntriv,ntrio,ifd4,ifa2)
      implicit integer(a-z)
      parameter (numint=24)
      REAL escf,enuc,eigval(nbf)
      dimension flov(nirred,4),npr(nirred,3,2),fpbka(1),
     1   fpbkb(1),fpbkc(1),fpbkd1(1),fpbkd2(2),fpbkd3(1),
     2     fpbke1(1),fpbkf(1),orbsym(nbf),fsec(nlist)
      dimension ic(numint),ptocc(nbf),iff(nirred),
     1  ife1(nirred),ifd3(nirred),ifd2(nirred),ifd1(nirred),
     2  ifc(nirred),ifb(nirred),ifa(nirred),ioff(nbf),isqov(no*nv),
     3isqvv(nv*nv),isqoo(no*no),itrio(ntrio),itriv(ntriv),ifd4(nirred+1)
     4 ,isqvo(no*nv),ifa2(nirred+1)
      call rfile(itap57)
      nsqvv=nv*nv
      nsqov=no*nv
      nsqoo=no*no
      ic(1)=nirred
      ic(2)=nbf
      ic(3)=no
      ic(4)=nv
      ic(5)=ntriv
      ic(6)=ntrio
      ic(7)=nsqvv
      ic(8)=nsqoo
      ic(9)=nsqov
      ic(10)=lnbkt
      ic(11)=nlist
      ndimwx=0
      ndimw=0
      do 509 jki=1,nirred
      ndimwx=npr(jki,1,2)*npr(jki,1,2)+ndimwx
      ndimw=npr(jki,3,2)*npr(jki,3,2)+ndimw
509   continue
      if(ndimwx.gt.ndimw) ndimw=ndimwx
      call tit_izo(ifd4,nirred+1)
      call tit_izo(ifa2,nirred+1)
      do 510 jki=1,nirred
      ifa2(jki+1)=ifa2(jki)+npr(jki,1,2)*npr(jki,1,1)
      ifd4(jki+1)=ifd4(jki)+npr(jki,1,2)*npr(jki,2,2)
510   continue
      if(nbka.eq.0) nbka=1
      if(nbkb.eq.0) nbkb=1
      if(nbkc.eq.0) nbkc=1
      if(nbkd1.eq.0) nbkd1=1
      if(nbkd2.eq.0) nbkd2=1
      if(nbkd3.eq.0) nbkd3=1
      if(nbke1.eq.0) nbke1=1
_IF1()ctjl  if(nbke2.eq.0) nbke2=1
      if(nbkf.eq.0) nbkf=1
      do 101 i=1,nbka+1
      fpbka(i)=fpbka(i)+1
101   continue
      do 102 i=1,nbkb+1
      fpbkb(i)=fpbkb(i)+1
102   continue
      do 103 i=1,nbkc+1
      fpbkc(i)=fpbkc(i)+1
103   continue
      do 104 i=1,nbkd1+1
      fpbkd1(i)=fpbkd1(i)+1
104   continue
      do 105 i=1,nbkd2+1
      fpbkd2(i)=fpbkd2(i)+1
105   continue
      do 106 i=1,nbkd3+1
      fpbkd3(i)=fpbkd3(i)+1
106   continue
      do 107 i=1,nbke1+1
      fpbke1(i)=fpbke1(i)+1
107   continue
      do 109 i=1,nbkf+1
      fpbkf(i)=fpbkf(i)+1
109   continue
      ic(12)=nbka
      ic(13)=nbkb
      ic(14)=nbkc
      ic(15)=nbkd1
      ic(16)=nbkd2
      ic(17)=nbkd3
      ic(18)=nbke1
      ic(19)=nbkf
      ic(20)=ndimt1
      ic(21)=ndimt2
      norbs=nbf
      ic(22)=norbs
      ic(23)=ndimw
      call wwritw(itap57,ic,numint,1,junk)
      call wwritw(itap57,escf,intowp(1),junk,junk)
      call wwritw(itap57,enuc,intowp(1),junk,junk)
      call wwritw(itap57,flov,nirred*4,junk,junk)
      call wwritw(itap57,npr,nirred*6,junk,junk)
      call wwritw(itap57,fpbka,nbka+1,junk,junk)
      call wwritw(itap57,fpbkb,nbkb+1,junk,junk)
      call wwritw(itap57,fpbkc,nbkc+1,junk,junk)
      call wwritw(itap57,fpbkd1,nbkd1+1,junk,junk)
      call wwritw(itap57,fpbkd2,nbkd2+1,junk,junk)
      call wwritw(itap57,fpbkd3,nbkd3+1,junk,junk)
      call wwritw(itap57,fpbke1,nbke1+1,junk,junk)
      call wwritw(itap57,fpbkf,nbkf+1,junk,junk)
      call wwritw(itap57,fsec,nlist,junk,junk)
      call wwritw(itap57,orbsym,nbf,junk,junk)
      call wwritw(itap57,itriv,ntriv,junk,junk)
      call wwritw(itap57,itrio,ntrio,junk,junk)
      call wwritw(itap57,isqvv,nv*nv,junk,junk)
      call wwritw(itap57,isqoo,no*no,junk,junk)
      call wwritw(itap57,isqov,no*nv,junk,junk)
      call wwritw(itap57,isqvo,no*nv,junk,junk)
      call wwritw(itap57,ioff,nbf,junk,junk)
      call wwritw(itap57,ifa,nirred,junk,junk)
      call wwritw(itap57,ifa2,nirred,junk,junk)
      call wwritw(itap57,ifb,nirred,junk,junk)
      call wwritw(itap57,ifc,nirred,junk,junk)
      call wwritw(itap57,ifd1,nirred,junk,junk)
      call wwritw(itap57,ifd2,nirred,junk,junk)
      call wwritw(itap57,ifd3,nirred,junk,junk)
      call wwritw(itap57,ifd4,nirred,junk,junk)
      call wwritw(itap57,ife1,nirred,junk,junk)
      call wwritw(itap57,iff,nirred,junk,junk)
      call wwritw(itap57,ptocc,nbf,junk,junk)
      call wwritw(itap57,eigval,intowp(nbf),junk,junk)
      call rclose(itap57,3)
      iprtxx=0
      if (iprtxx.ne.0)then
      write(6,*) '  in subroutine wrt57 '
      write(6,*) ' nirred ',nirred
      write(6,*) ' ntriv ',ntriv
      write(6,*) ' ntrio ',ntrio
      write(6,*) ' nsqvv ',nsqvv
      write(6,*) ' nsqoo ',nsqoo
      write(6,*) ' nsqov ',nsqov
      write(6,*) ' lnbkt ',lnbkt
      write(6,*) ' nlist ',nlist
      write(6,*) ' nbka ',nbka
      write(6,*) ' nbkb ',nbkb
      write(6,*) ' nbkc ',nbkc
      write(6,*) ' nbkd1 ',nbkd1
      write(6,*) ' nbkd2 ',nbkd2
      write(6,*) ' nbkd3 ',nbkd3
      write(6,*) ' nbke1 ',nbke1
_IF1()ctjl  write(6,*) ' nbke2 ',nbke2
      write(6,*) ' nbkf ',nbkf
      write(6,*) ' ndimt1 ',ndimt1
      write(6,*) ' ndimt2 ',ndimt2
      write(6,*) ' norbs ',norbs
      write(6,*) ' ndimw ',ndimw
_IFN(win95,linux)
      write(6,*) ' fpbka  ',fpbka
      write(6,*) ' fpbkb  ',fpbkb
      write(6,*) ' fpbkc  ',fpbkc
      write(6,*) ' fpbkd1 ',fpbkd1
      write(6,*) ' fpbkd2 ',fpbkd2
      write(6,*) ' fpbkd3 ',fpbkd3
      write(6,*) ' fpbke1 ',fpbke1
      write(6,*) ' fpbkf  ',fpbkf
      write(6,*) ' isqoo ',isqoo
      write(6,*) ' itrio ',itrio
      write(6,*) ' itriv ',itriv
      write(6,*) ' isqvv ',isqvv
      write(6,*) ' isqov ',isqov
      write(6,*) ' ifa2 ',ifa2
      write(6,*) ' ifc ',ifc
      write(6,*) ' fsec ',fsec
      write(6,*) ' npr ',npr
      write(6,*) ' ifd1 ',ifd1
      write(6,*) ' ifd3 ',ifd3
_ENDIF
      endif
      return
      end
      subroutine ver_tsort(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/tsort.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
