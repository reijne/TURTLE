      subroutine read_aivb_input(readnr)
c...  reads AiVB input, makes necessary checks and sets essential variables
      implicit none
c
INCLUDE(../m4/common/sizes)
c
INCLUDE(common/turtleparam)
INCLUDE(common/aivb)
c
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/runlab)
INCLUDE(../m4/common/iofile)
      character(len=1)  label5
      character(len=8)  label1,label2,label3,label4
      character(len=10) ctmp1
      integer :: readnr,cntr,itemp,i,j,k,idbg,ii,jrec,jump,ij,jj,ictmp1
      integer :: ist,iend,iatnel
c
      idbg = aimo_debugpr(2)
      cntr = 0
      itemp = 1
      aimo_atcntr(readnr) = 0
      if ( idbg .eq. 1 ) then
        print *,'     ==> READ AiVB INPUT subroutine - begin'
        print *,' readnr: ',readnr
        print *,' aivb_optorb: ',aivb_optorb
      end if
      if ( readnr .eq. 1 ) then
        aimo_trconf(readnr,1,1) = 0
        ii = 0
        call inpa(label1)
        if ( label1(1:2) .ne. '  ' ) then
          call cinput(jrec,jump,-1)
          call inpi(aivb_inpspp(readnr))
          write(iwr,'(a,I3)') ' Picking spin path ',aivb_inpspp(readnr)
        else
          aivb_inpspp(readnr) = 1
        end if
9004    call inpa(label1)
        if ( label1(1:2) .ne. '  ' ) then
          call cinput(jrec,jump,-1)
          ii = ii + 1
          call inpi(aimo_trconf(readnr,2,ii))
          aimo_trconf(readnr,1,1) = aimo_trconf(readnr,1,1)+ 1
          goto 9004
        end if
        if ( idbg .eq. 1 ) then
          write(iwr,'(1x,a,I3)') '_trconf 1: ',aimo_trconf(readnr,1,1)
          write(iwr,'(1x,a,20I3)') '_trconf 2: ',
     1          (aimo_trconf(readnr,2,ij),
     2          ij=1,aimo_trconf(readnr,1,1))
        end if

        call input
      else 
        call cinput(jrec,jump,-1)
      end if
9001  call inpa(label1)
      if ( idbg .eq. 1 ) print *,'label1 read aivb inp: ',label1(1:5)
      call cinput(jrec,jump,-1)
c...  we should have at the beginning something else than END or another AiVB
      if ( label1(1:3) .eq. 'end' ) then
           if ( cntr .eq. 0 ) then
             call vberr(' atomic state definitions expected first ')
           end if
      end if
c...  if everything besides END, we might have atomic state definitions, read it
      if ( label1(1:3) .ne. 'end' .or. label1(1:4) .ne. 'conf' ) then
        call inpa(label1)
        if ( idbg .eq. 1 ) print *,'label1 1: ',label1(1:5)

        cntr = cntr + 1
        if ( label1(1:8) .eq. zaname(cntr) ) then 
          aimo_atcntr(readnr) = aimo_atcntr(readnr) + 1
          aivb_atnmrdin(cntr,readnr) = label1
          aimo_atindex(aimo_atcntr(readnr),1) = cntr
        else
         call vberr(' wrong atom symbol or wrong AiVB definition - check
     & Your input ')
        end if

        j = 0

        aivb_zchr(cntr) = ' '
        call inpa(label2)
        if ( idbg .eq. 1 ) print *,'label2 1: ',label2(1:5)

        aimo_atindex(cntr,4) = 0
        if ( label2(1:2) .eq. '+ ' ) then
          aivb_zchr(cntr) = '- '
          aimo_atindex(cntr,4) = aimo_atindex(cntr,4) + 1
          call inpa4(label2)
          if ( idbg .eq. 1 ) print *,'label2 2: ',label2(1:5)
        else if ( label2(1:2) .eq. '++' ) then
          aivb_zchr(cntr) = '2-'
          aimo_atindex(cntr,4) = aimo_atindex(cntr,4) + 2
          call inpa4(label2)
          if ( idbg .eq. 1 ) print *,'label2 3: ',label2(1:5)
        else if ( label2(1:2) .eq. '- ' ) then
          aivb_zchr(cntr) = '+ '
          aimo_atindex(cntr,4) = aimo_atindex(cntr,4) - 1
          call inpa4(label2)
          if ( idbg .eq. 1 ) print *,'label2 4: ',label2(1:5)
        else if ( label2(1:2) .eq. '--' ) then
          aivb_zchr(cntr) = '2+'
          aimo_atindex(cntr,4) = aimo_atindex(cntr,4) - 2
          call inpa4(label2)
          if ( idbg .eq. 1 ) print *,'label2 5: ',label2(1:5)
        end if
        iatnel = czan(aimo_atindex(cntr,1)) + aimo_atindex(cntr,4)
        if ( idbg .eq. 1 ) print *,'n el. in at.: ',iatnel

        j = 0
        k = 0
c...  now, the same goes with the multiplicity names (sing, doub, trip, quar, quin)
        do j=1,6
          if ( label2(1:1) .ne. aivb_multname(j) ) k = k+1
          if ( label2(1:1) .eq. aivb_multname(j) ) then 
            aivb_multrdin(cntr,readnr) = label2
          end if
        end do
        if ( k .eq. 6 ) 
     &  call vberr(' wrong multiplicity label or wrong AiVB definition -
     & check Your input ')

        call inpa(label3)
        if ( idbg .eq. 1 ) print *,'label3 3: ',label3,':'
        j = 0
        k = 0
c...  and with the atomic state names (s, p, d)
        do j=1,4
          if ( label3(1:1) .ne. aivb_atstname(j) ) k = k+1
          if ( label3(1:1) .eq. aivb_atstname(j) ) then
            aivb_atstrdin(cntr,readnr) = label3
          end if
        end do
        if ( k .eq. 4 )
     &  call vberr(' wrong atomic state symbol or wrong AiVB defnition -
     & check Your input ')

        i = cntr
c
c...  row 1, 2 and 3 S-block doublets
c
          if ( iatnel .eq. 1 .or. 
     &         iatnel .eq. 3 .or. 
     &         iatnel .eq. 11 ) then
            if ( iatnel .eq. 1 ) then
              aimo_coreorb(i) = 0
              aimo_atindex(i,2) = 7
              aimo_atindex(i,3) = 1
            else if ( iatnel .eq. 3 ) then
              aimo_coreorb(i) = 1
              aimo_atindex(i,2) = 7
              aimo_atindex(i,3) = 1
            else if ( iatnel .eq. 11 ) then 
              aimo_coreorb(i) = 5
              aimo_atindex(i,2) = 7
              aimo_atindex(i,3) = 1
            end if
            label4 = label2(1:1)//label3(1:1)
            itemp = itemp + 1
            if ( label2 .ne. '2s' ) then
              call vberr(' incorrect atomic state for S-block atom - che
     &ck Your input ')
            end if
c
c...  row 1,2 and 3 S-block singlets
c
          else  if ( iatnel .eq. 2 .or. 
     &               iatnel .eq. 4 .or. 
     &               iatnel .eq. 12 ) then
            if ( iatnel .eq. 2 ) then
              aimo_coreorb(i) = 0
              aimo_atindex(i,2) = 8
              aimo_atindex(i,3) = 1
            else if ( iatnel .eq. 4 ) then
              aimo_coreorb(i) = 1
              aimo_atindex(i,2) = 8
              aimo_atindex(i,3) = 1
            else if ( iatnel .eq. 12 ) then 
              aimo_coreorb(i) = 5
              aimo_atindex(i,2) = 8
              aimo_atindex(i,3) = 1
            end if
            label4 = label2(1:1)//label3(1:1)
            itemp = itemp + 1
            if ( label2 .ne. '1s' ) then
              call vberr(' incorrect atomic state for S-block atom - che
     &ck Your input ')
            end if
c
c...  row 2
c
          else if ( iatnel .le. 10 .and. 
     &              iatnel .ge. 3 ) then
c...  row 2, P-block
            if ( iatnel .le. 10 .and. 
     &           iatnel .ge. 5 ) then
              aimo_coreorb(i) = 1
              aimo_atindex(i,2) = iatnel - 
     &                            aimo_coreorb(i)*2 - 2
              jj = 0
              label4 = label2(1:1)//label3(1:1)
              itemp = itemp + 1
              do ii=1,aimo_atstcom(aimo_atindex(i,2),1)
                if ( label2 .ne. aivb_allatst(aimo_atindex(i,2),ii) )
     & then
                  jj = jj + 1
                else
                  aimo_atindex(i,3) = ii
                end if
              end do
              if ( jj .eq. aimo_atstcom(aimo_atindex(i,2),1) ) then
                call vberr(' incorrect atomic state for P-block atom -
     & check Your input ')
              end if
            end if
c
c...  row 3
c
          else if ( iatnel .le. 18 .and. 
     &              iatnel .ge. 11 ) then
c...  row 3, P-block
            if ( iatnel .le. 18 .and. 
     &           iatnel .ge. 13 ) then
              aimo_coreorb(i) = 5
              aimo_atindex(i,2) = iatnel - 
     &                            aimo_coreorb(i)*2 - 2
              jj = 0
              label4 = label2(1:1)//label3(1:1)
              itemp = itemp + 1
                do ii=1,aimo_atstcom(aimo_atindex(i,2),1)
                  if ( label2 .ne. aivb_allatst(aimo_atindex(i,2),
     & ii) ) then
                    jj = jj + 1
                  else
                    aimo_atindex(i,3) = ii
                  end if
                end do
                if ( jj .eq. aimo_atstcom(aimo_atindex(i,2),1) ) then
                  call vberr(' incorrect atomic state for P-block atom -
     & check Your input ')
                end if
            end if
c
c...  0 s electrons
c
          else if ( iatnel .eq. 0 ) then
            aimo_coreorb(i) = 0
            aimo_atindex(i,2) = 9
            aimo_atindex(i,3) = 1
c
c...  given atom is not supported (rows 1, 2 and 3 are)
c
          else
            call vberr(' given atom not yet supported by AiVB - try usin
     &g STRUC ')
          end if
c        end do


        call input
        call inpa(label1)
        if ( idbg .eq. 1 ) print *,'label1 4: ',label1(1:5)
        if ( label1(1:4) .eq. 'conf' .or. label1(1:3) .eq. 'end' ) then
          if ( label1(1:4) .eq. 'conf' ) then
            aimo_trconf(readnr+1,1,1) = 0
            ii = 0
            call inpa(label1)
            if ( label1(1:2) .ne. '  ' ) then
              call cinput(jrec,jump,-1)
              call inpi(aivb_inpspp(readnr+1))
              write(iwr,'(a,I3)') ' Picking spin path ',
     1                            aivb_inpspp(readnr+1)
            else
              aivb_inpspp(readnr+1) = 1
            end if
9005        call inpa(label1)
            if ( label1(1:2) .ne. '  ' ) then
              call cinput(jrec,jump,-1)
              ii = ii + 1
              call inpi(aimo_trconf(readnr+1,2,ii))
              aimo_trconf(readnr+1,1,1) = aimo_trconf(readnr+1,1,1)+ 1
              goto 9005
            end if
            if ( idbg .eq. 1 ) then
              write(iwr,'(a,I3)') ' _trconf 1: ',
     1              aimo_trconf(readnr+1,1,1)
              write(iwr,'(a,20I3)') ' _trconf 2: ',
     1        (aimo_trconf(readnr+1,2,ij),
     2        ij=1,aimo_trconf(readnr+1,1,1))
            end if
          end if
          goto 9002
        end if
        if ( cntr .gt. nat ) 
     &  call vberr(' too many definitions, more than atoms in molecule -
     & check Your input')
        call cinput(jrec,jump,-1)
        goto 9001
      end if

9002  continue
      
      if ( aimo_atcntr(readnr) .ne. nat ) 
     &  call vberr(' missing some atoms in declarations - check Your inp
     &ut ')

      call cinput(jrec,jump,-1)
      if ( idbg .eq. 1 ) then
        write(iwr,'(a,30I3)') ' coreorb   : ',(aimo_coreorb(i), i=1, 
     &                           aimo_atcntr(readnr))
        write(iwr,'(a,30I3)') ' atindex,1 : ',(aimo_atindex(i,1),i=1,
     &                           aimo_atcntr(readnr))
        write(iwr,'(a,30I3)') ' atindex,2 : ',(aimo_atindex(i,2),i=1,
     &                           aimo_atcntr(readnr))
        write(iwr,'(a,30I3)') ' atindex,3 : ',(aimo_atindex(i,3),i=1,
     &                           aimo_atcntr(readnr))
        write(iwr,'(a,30I3)') ' atindex,4 : ',(aimo_atindex(i,4),i=1,
     &                           aimo_atcntr(readnr))
        write(iwr,'(a,30a,a)') ' atnmrdin  : ',(aivb_atnmrdin(i,readnr),
     &                         ' | ',i=1,cntr)
        write(iwr,'(a,30a,a)') ' multrdin  : ',(aivb_multrdin(i,readnr),
     &                         ' | ',i=1,cntr)
        write(iwr,'(a,30a,a)') ' atstrdin  : ',(aivb_atstrdin(i,readnr),
     &                         ' | ',i=1,cntr)
        print *,'     ==> READ AIVB INPUT subroutine - end'
      end if

      return
      end

c*****************************************************************************************
c
c...  generate all determinants and unique configurations 
c...  by knowing in what state given atom is (say: P triplet)
c
      subroutine generate_aivb_drv(coefff,deter,readnr,detnr,confnr,
     &           numb,maxx,nelec,ncore_i,mdet,scr)
c
      implicit none
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
c
INCLUDE(common/turtleparam)
INCLUDE(common/aivb)
c
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/runlab)
c
      common/hans/ msocc,mstrco,mdetco,mdetst,mperm,mstruc
c
      common /aimotmp/ aimo_lterms,adbg9,adbg10,vbcrreadnr
      integer aimo_lterms(maxcon,maxact),adbg9,adbg10,vbcrreadnr
c
      integer :: msocc,mstrco,mdetco,mdetst,mperm,mstruc,msp,mpair
      integer :: oldmsocc,oldmstrco,oldmdetco,oldmdetst,oldmperm
      integer :: oldmstruc,oldmsp
c
      integer :: nnd,nnc,ndall,ncall,nbothsp
      integer :: iconfnalpha,n2s,n2s1,n2s2,new1,new2,newtont
      REAL :: iconfspin,iconfmult
c
c...  'nnd'      : number of singly occupied orbitals, for which all
c...               possible, unique couplings will be generated
c...  'nnc'      : number of all possible atomic state's variations (sum over atoms)
c...  'ndall'    : max. number of possible couplings between singly
c...               occupied orbitals
c...  'ncall'    : max. number of possible couplings between all atomic
c...               state's variations of all atoms (product over atoms)
c...  'confnr'   : # of configuration in given readnr, at the output
c...  'detnr'    : # of determinants in given readnr, at the output
c...  'confnrold': contains number of configurations from previous input read
c...  'detnrold' : same as confnrold (number of determinants)
c...  'nconfv'   : # of elconfvaria configurations
c
      integer :: nelec,readnr,detnr,confnr,numb,deter(nelec,*),ij
      integer :: ncore_i,icore,mdet,detnrold,ifact,fact,igrp,nog
      integer :: confnrold,i,j,k,ii,jhi,jadd,jaddend,iposs,ielsign
      integer :: idbg3,idbg4,idbg7,status,itmp1,itmp2,idet,iconf
      integer :: inoncore,ihilimit,jx,jy,jz,jj,ising,itmp3,nconfv
      integer :: ll,l,inc,jatom,nposs,kind,nposssum,jtr
      integer :: nsingat(nat),nms,isign,idbg5,idbg6,maxx,itmp4
      integer :: kxind,kyind,kzind,kuind,klimit,kadd,kgrp,ksumdets
      integer :: scr(*),mscr,igmem_max_memory,nstruc,ndetst,itmp22
      integer :: anspath,aispath,ncc,annwdet,anconf,ndetp,ndetpc
      integer :: inconfs,indeco,nalphapconf,idbg9,idbg10,idbg11
      REAL :: coefff(mdet),anorm
      logical :: allout,projec,arumer
c
      integer, allocatable, dimension(:) :: kiconfindx,ilast,kidetps
      integer, allocatable, dimension(:) :: klogic,krumer,kdetps
      integer, allocatable, dimension(:) :: koccn,korder
      integer, allocatable, dimension(:,:) :: kconf,kdeter,kspinp
      integer, allocatable, dimension(:,:) :: elconfvaria,kconfdef
      integer, allocatable, dimension(:,:) :: kstruc,kindx,kspdet
      integer, allocatable, dimension(:,:) :: ksptmp,kcind
      REAL, allocatable, dimension(:) :: kscr,kcoefff
      REAL, allocatable, dimension(:,:) :: kcoeff
      integer, allocatable, dimension(:,:,:) :: kconfdetdef,knconf
c
c...  _debugpr explained in common/aivb
      idbg3 = aimo_debugpr(3)
      idbg4 = aimo_debugpr(4)
      idbg5 = aimo_debugpr(5)
      idbg6 = aimo_debugpr(6)
      idbg7 = aimo_debugpr(7)
      idbg9 = aimo_debugpr(9)
      idbg10 = aimo_debugpr(10)
      idbg11 = aimo_debugpr(11)
      confnrold = confnr
      detnrold = detnr
      detnr = 0
      confnr = 0
      vbcrreadnr = readnr

      if ( idbg3 .eq. 1 ) 
     1     print *,'     ==> GENERATE AIVB DRIVER subroutine - begin'
c...  settin up nn value - number of singly occupied orbitals
c...  from which all possible unique couplings will be made
c
c...  first, setting up total number of core orbitals
      icore = 0
      do i=1,aimo_atcntr(readnr),1
        icore = icore + aimo_coreorb(i)
      end do
c...  icore is a local number, setting up global core number
c...  IF CORE in CONF hasn't been specified 
      if ( idbg3 .eq. 1 ) then
        write(iwr,'(a,I3,a,I7)') ' readnr: ',readnr,'  mdet: ',mdet
        write(iwr,'(a,I3)') ' maxx: ',maxx
        write(iwr,'(a,I3)') ' icore: ',icore
        write(iwr,'(a,I3,a,I3,a,I3,a,I3,a,I3,a,I3)') ' nelec: ',nelec,
     1  '  ncore_i: ',ncore_i,'  detnr: ',detnr,'  detnrold: ',detnrold
     2  ,'  confnr: ',confnr,'  confnrold: ',confnrold
      end if

c...  setting temp ncore_i for orbital optimization
      if ( aivb_optorb ) then
        ncore_i = icore
      end if

c...  beginning of the atomic state's variations coupling part <==

c...  setting up number of all variations (sum; nnc), and
c...  number of all possible couplings between them (product; ncall)
      nnc = 0
      ncall = 1
      do i=1,aimo_atcntr(readnr),1
       nnc = nnc + aimo_atstcom(aimo_atindex(i,2),aimo_atindex(i,3)+1)
       ncall = ncall*aimo_atstcom(aimo_atindex(i,2),aimo_atindex(i,3)+1)
      end do

      if ( idbg3 .eq. 1 ) then
        write(iwr,'(a,I3,a,I3)') ' nnc: ',nnc,'    ncall: ',ncall
      end if

      ALLOCATE(kiconfindx(nnc), STAT=status)
      call allocmem('kiconfindx',nnc,'gen_aivb_drv','vbaivb.m',status)

      ALLOCATE(kconf(ncall,aimo_atcntr(readnr)), STAT=status)
      call allocmem('kconf',nnc*aimo_atcntr(readnr),
     1'gen_aivb_drv','vbaivb.m',status)

      igrp = 1
      ii = 1
      do j=1,aimo_atcntr(readnr),1
        aimo_conford(j) = aimo_atstcom(aimo_atindex(j,2),
     &                    aimo_atindex(j,3)+1)
        do k=1,aimo_conford(j)
          kiconfindx(ii) = k
          aimo_confgrps(ii) = igrp
          ii = ii + 1
        end do
        igrp = igrp + 1
      end do
c...  nog - number of groups which is equal to number of atoms
      nog = aimo_atcntr(readnr)
      if ( idbg3 .eq. 1 ) then
        write(iwr,'(1X,a,30I3)') 'iconfindx : ',
     1        (kiconfindx(i), i=1,nnc)
        write(iwr,'(1X,a,30I3)') '_conford  : ',
     1        (aimo_conford(i), i=1,nog)
        write(iwr,'(1X,a,30I3)') '_confgrps : ',
     1        (aimo_confgrps(i), i=1,nnc)
      end if

c...
c...  generating atomic states couplings
c...
      call gen_aivbconf_drv(kconf,nnc,nog,ncall,nconfv,aimo_confgrps)

      if ( idbg3 .eq. 1 ) then
        write(iwr,1010) 'mdet: ',mdet,'  nnc: ',nnc,'  nog: ',nog,
     1          '  ncall: ',ncall,'  nconfv: ',nconfv
1010    format(1X,5(a,I6))
      end if

      if ( idbg4 .eq. 1 ) then
        do i=1,ncall
          write(iwr,'(1X,a,30I3)') 'kconf: ',
     1         (kconf(i,j), j=1,nog)
        end do
      end if

c...  transforming 'kconf' to 'kconf' corresponding to the '_statedef' table
      do i=1,ncall
        do j=1,nog
          itmp1 = kconf(i,j)
          itmp2 = kiconfindx(itmp1)
          if ( idbg4 .eq. 1 ) then
            write(iwr,'(1X,a,I3)') 'kconf b: ',
     1      kconf(i,j)
            write(iwr,'(1X,a,I3)') 'itmp1  : ',itmp1
            write(iwr,'(1X,a,I3)') 'itmp2  : ',itmp2
          end if
          kconf(i,j) = itmp2
          if ( idbg4 .eq. 1 ) then
            write(iwr,'(1X,a,I3)') 'kconf a: ',
     1      kconf(i,j)
          end if
        end do
      end do

      DEALLOCATE(kiconfindx, STAT=status)
      call deallocmem('kiconfidx','gen_aivb_drv','vbaivb.m',status)

      if ( idbg4 .eq. 1 ) then
        do i=1,ncall
          write(iwr,'(1X,a,30I3)') 'kconf after: ',
     1          (kconf(i,j), j=1,nog)
        end do
      end if

c...  define electronic configurations in terms of determinant numbers
      ALLOCATE(kconfdetdef(2,ncall,nog*3*3), STAT=status)
      call allocmem('kconfdetdef',2*ncall*nog*3*3,'gen_aivb_drv',
     1'vbaivb.m',status)

      ALLOCATE(knconf(3,ncall,nog*3), STAT=status)
      call allocmem('knconf',3*ncall*nog*3,'gen_aivb_drv',
     1'vbaivb.m',status)

      nposssum = 0
      do inc=1,ncall
        nposs = 1
        kind = 1
        kadd = 1
        ksumdets = 0
        do jatom=1,nog
          kxind = aimo_atindex(jatom,2)
          kyind = aimo_atindex(jatom,3)
          kzind = kconf(inc,jatom)
          klimit = aimo_statedef(kxind,kyind,kzind,1)
          nposs = nposs*klimit
          ksumdets = ksumdets + klimit
          if ( idbg4 .eq. 1 ) then
            write(iwr,'(a,I3)') ' kxind   : ',kxind
            write(iwr,'(a,I3)') ' kyind   : ',kyind
            write(iwr,'(a,I3)') ' kzind   : ',kzind
            write(iwr,'(a,I3)') ' klimit  : ',klimit
            write(iwr,'(a,I3)') ' nposs   : ',nposs
          end if
          kuind = 2
          do k=1,klimit
            kconfdetdef(1,inc,kind  ) = kxind
            kconfdetdef(1,inc,kind+1) = aimo_statedef(kxind,kyind,kzind,
     &                                  kuind  )
            kconfdetdef(1,inc,kind+2) = aimo_statedef(kxind,kyind,kzind,
     &                                  kuind+1)
            kuind = kuind + 2
            kind = kind + 3
            knconf(2,inc,kadd) = jatom
            knconf(3,inc,kadd) = kadd
            if ( idbg4 .eq. 1 ) then
              write(iwr,'(a,I3)') ' grp: ',knconf(2,inc,kadd)
            end if
            kadd = kadd + 1
          end do
          kconfdetdef(2,inc,jatom) = klimit
          knconf(1,inc,2) = ksumdets
          if ( idbg4 .eq. 1 ) then
            write(iwr,'(a,I3)') ' ksumdets: ',ksumdets
            write(iwr,2228)
          end if
        end do
        if ( idbg4 .eq. 1 )
     1  write(iwr,2229)
2229    format(72('='))
2228    format(50('-'))
        knconf(1,inc,1) = nposs
        nposssum = nposssum + nposs
      end do

      if ( idbg4 .eq. 1 ) then
        do inc=1,ncall
          do jatom=1,nog
            kadd = (jatom-1)*3
            write(iwr,'(2(I3,a))') inc,' conf; ',jatom,' atom '
            write(iwr,'(10x,20(3I3,a))') (
     1      kconfdetdef(1,inc,kind+kadd),kconfdetdef(1,inc,kind+1+kadd),
     2      kconfdetdef(1,inc,kind+2+kadd),'  ', kind=1,
     4      kconfdetdef(2,inc,jatom)*3,3)
          end do
          write(iwr,'(a,I3)') '  # of possible coupling: ',
     1                            knconf(1,inc,1)
          write(iwr,'(a,60I3)') '  determinant groups for given conf: ',
     1    (knconf(2,inc,k), k=1,knconf(1,inc,2))
          write(iwr,'(a,60I3)') '  how many: ',
     1    (knconf(3,inc,k), k=1,knconf(1,inc,2))
          print *,'==========> next conf ================> '
        end do
        write(iwr,'(a,I3)') '  nposs sum: ',nposssum
      end if

      ALLOCATE(kcind(ncall,aimo_atcntr(readnr)), STAT=status)
      call allocmem('kcind',nnc*aimo_atcntr(readnr),
     1'gen_aivb_drv','vbaivb.m',status)
      do i=1,ncall
        do j=1,nog
          kcind(i,j) = kconf(i,j)
        end do
      end do
      if ( idbg4 .eq. 1 ) then
        do i=1,ncall
          write(iwr,'(a,30I3)') ' kcind: ',
     1          (kcind(i,j), j=1,nog)
        end do
      end if

      DEALLOCATE(kconf, STAT=status)
      call deallocmem('kconf','gen_aivb_drv','vbaivb.m',status)

      ALLOCATE(kconf(nposssum,2*nog), STAT=status)
      call allocmem('kconf',nposssum*2*nog,
     1'gen_aivb_drv','vbaivb.m',status)

      ALLOCATE(elconfvaria(maxx,nposssum), STAT=status)
      call allocmem('elconfvaria',maxx*nposssum,
     1'gen_aivb_drv','vbaivb.m',status)

      ALLOCATE(kconfdef(nposssum,2*nog+2), STAT=status)
      call allocmem('kconfdef',nposssum*2*nog+2,
     1'gen_aivb_drv','vbaivb.m',status)

      kxind = 1
      do inc=1,ncall

        if ( idbg5 .eq. 1 ) then
          write(iwr,*) ' '
          write(iwr,2230)
          write(iwr,'(a,I3)') ' generating variations for conf. # ',inc
        end if

c...  copy the group numbers of determinants to a single vector
        do k=1,knconf(1,inc,2)
          aimo_confgrps(k) = knconf(2,inc,k)
        end do

c...  generate the possible couplings between different determinants describing given
c...  atomic state withing the same 'inc'
        call gen_aivbconf_drv(kconf(kxind,1),knconf(1,inc,2),nog,
     &                        nposssum,nconfv,aimo_confgrps)

c...  loop over the # of possible couplings withing given 'inc'
       
        do kind=kxind,kxind+knconf(1,inc,1)-1
          if ( idbg5 .eq. 1 ) then 
              write(iwr,*) ' '
              write(iwr,'(a,I3)') ' kind: ',kind
          end if
c...  loop over the # of atoms
          kconfdef(kind,1) = inc
          isign = 1
          do jatom=1,nog
            idet = kconf(kind,jatom)
            kuind = (idet-1)*3 + 1
            itmp1 = kconfdetdef(1,inc,kuind)
            isign = isign*kconfdetdef(1,inc,kuind+1)
            itmp2 = kconfdetdef(1,inc,kuind+2)
            kzind = (jatom-1)*2 + 1 + 2
            kconfdef(kind,kzind) = itmp1
            kconfdef(kind,kzind+1) = itmp2
            if ( idbg5 .eq. 1 ) then
              write(iwr,'(a,I3)') ' jatom: ',jatom
              write(iwr,'(a,I3)') ' idet : ',idet
              write(iwr,'(a,I3)') ' kuind: ',kuind
              write(iwr,'(a,I3)') ' itmp1: ',itmp1
              write(iwr,'(a,I3)') ' itmp2: ',itmp2
              write(iwr,'(a,I3)') ' kzind: ',kzind
              write(iwr,'(3(a,I3))') ' kconfdef(',
     1        kind,',',kzind,'): ',kconfdef(kind,kzind)
              write(iwr,'(3(a,I3))') ' kconfdef(',
     1        kind,',',kzind+1,'): ',kconfdef(kind,kzind+1)
            end if
            if ( idbg5 .eq. 1 )
     1      print *,'          -------- next jatom ---------- '
          end do
          kconfdef(kind,2) = isign
          if ( idbg5 .eq. 1 )
     1    print *,' -------- next kind ---------- '
        end do

        kxind = kxind + knconf(1,inc,1)

        if ( idbg5 .eq. 1 ) then
          write(iwr,'(a,I3)') ' next kxind: ',kxind
          do i=1,kxind-1
            write(iwr,'(a,30I3)') ' tmp kconf: ',
     1            (kconf(i,j), j=1,nog)
          end do
        end if

        if ( idbg5 .eq. 1 ) then
          write(iwr,*) ' '
          write(iwr,2230)
        end if
2230    format(45(';;'))

      end do

      DEALLOCATE(kconfdetdef, STAT=status)
      call deallocmem('kconfdetdef','gen_aivb_drv','vbaivb.m',status)

c...  'nconfv' effectively becomes now equal to 'nposssum'
c...  which is the number of all the possible couplings between the determinants involved
c...  in the description of all the atomic states of atoms
      nconfv = nposssum

c...  'kconf' stores first the couplings between the atomic states, and their projections, of different atoms
c...  'kconf NEW' stores then the couplings between all the determinantas involved in the description of all the atomic states
c...  'kconfdef' stores the pairs of numbers x1y1 x2y2 ... describing the given electronic configuration:
c...  <x1> is the x-index in _confdef and <y1> is the y-index in _confdef (for atom 1) and so on for atom 2 (x2 y2)
      if ( idbg3 .eq. 1 ) then
        write(iwr,*) ' '
        do i=1,nconfv
          write(iwr,'(a,30I3)') ' ==> kconf NEW: ',
     1          (kconf(i,j), j=1,nog)
        end do
        write(iwr,*) ' '
        do i=1,nconfv
          write(iwr,'(a,I3,a,30I3)') ' ==> kconfdef ',i,': ',
     1          (kconfdef(i,j), j=1,nog*2+2)
        end do
      end if

      DEALLOCATE(kconf, STAT=status)
      call deallocmem('kconf','gen_aivb_drv','vbaivb.m',status)

      ALLOCATE(ilast(nconfv*2), STAT=status)
      call allocmem('ilast',nconfv*2,'gen_aivb_drv','vbaivb.m',status)

c...  jhi numerates index, thus keeping last, highes index number
c...  jadd keeps last, highest possible orbital number to add for next loop
c...  generate configuration atom by atom
c...  first, the loop over possible couplings between determinants of atomic states 'iconf'
      do iconf=1,nconfv
        jhi = 1
c...  first core (doubly occu.) orbitals - simple
        do k=1,ncore_i*2
          elconfvaria(jhi,iconf) = (k+1)/2
          jhi = jhi + 1
        end do
        jadd = ncore_i
c...  loop over atoms 'jatom'
        do jatom=1,nog
          jaddend = 0
          kind = (jatom-1)*2+1 + 2
          jx = kconfdef(iconf,kind)
          jy = kconfdef(iconf,kind+1)
          ihilimit = aimo_confdef(jx,jy,1)
          if ( idbg6 .eq. 1 ) then
            write(iwr,2232)
            write(iwr,'(3(a,I3))') ' kind: ',kind,'  jx: ',jx,
     1            '  jy: ',jy
            write(iwr,'(a,I3)') ' ihilimit: ',ihilimit
            write(iwr,'(2(a,I3))') ' jhi: ',jhi,'   jadd: ',jadd
            write(iwr,2232)
          end if
2232      format(25('--'))
c...  then loop over the rest, non-core electrons which come from
c...  given atomic state's determinant '_confdef'
          do k=1,ihilimit
c...  if the orbital number from '_confdef' multiplied by 100,
c...  then it's an empty orbital
            if ( aimo_confdef(jx,jy,k+1) .ge. 100 ) then
              if ( idbg6 .eq. 1 ) then
                write(iwr,'(a,I3,a,I3)') 'skipped orbital ',k,
     1              '   ',aimo_confdef(jx,jy,k+1)
              end if
c...  if it's an empty orbital, increase the orbital number by one
c...  every time You encounter a pair of empty orbitals
              if ( aimo_confdef(jx,jy,k+1) .eq. aimo_confdef(jx,jy,k) )
     1        jaddend = jaddend + 1
              if ( idbg6 .eq. 1 ) then
                write(iwr,'(a,I3)') ' >100 jadd: ',jadd
                write(iwr,'(a,I3)') ' >100 jaddend: ',jaddend
              end if
              goto 1122
            end if
c...  this part is executed only if an orbital is occupied (positivie or negative value in '_confdef')
c...  establish wether it's an alpha or beta orbital from '_confdef'
            jaddend = 0
            ielsign = 
     &      aimo_confdef(jx,jy,k+1)/ABS(aimo_confdef(jx,jy,k+1))
            if ( idbg6 .eq. 1 ) then
              write(iwr,'(5(a,I3))') 'jx: ',jx,'  jy: ',jy,
     1            '  aimo_confdef(jx,jy,k+1): ',aimo_confdef(jx,jy,k+1),
     2            '  k: ',k
              write(iwr,'(a,I3)') ' jadd: ',jadd
              write(iwr,'(a,I3)') ' sign: ',ielsign
            end if
c...  put an occupied orbital number on the next 'jhi' position
c...  multiply the orbital number and each addition to it, by a sign 'ielsign'
            elconfvaria(jhi,iconf) = 
     &        ielsign*ABS(aimo_confdef(jx,jy,k+1)) + 
     &        ielsign*jadd
            if ( idbg6 .eq. 1 ) then
              write(iwr,'(4X,a,I3,a,I3)') 'CONF. ELE.: ',
     1        aimo_confdef(jx,jy,k+1),'   ',elconfvaria(jhi,iconf)
            end if
c...  set the 'ilast' value for the given 'iconf' - has to be an absolute value
            ilast(iconf) = ABS(elconfvaria(jhi,iconf))
c...  increase the index in elconfvaria, for the given 'iconf'
            jhi = jhi + 1
1122        continue
          end do
c...  set the 'jadd' value which will be added to the next orbital numbers for the next 'jatom'
          jadd = ABS(elconfvaria(jhi-1,iconf)) + jaddend
          if ( idbg6 .eq. 1 ) then
            write(iwr,'(a,30I3)') ' CONF: ',
     1            (elconfvaria(jj,iconf), jj=1,jhi-1)
            write(iwr,'(a,I3,a,I3)') ' 1 jhi = ',jhi,'  1 jadd = ',jadd
            write(iwr,'(a,I3,a,I3)') ' ilast(',iconf,') = ',ilast(iconf)
            write(iwr,2231)
          end if
2231      format(35('--'))
        end do
        if ( idbg6 .eq. 1 ) then
          write(iwr,'(a,I3,a,I3)') ' 2 jhi = ',jhi,'  2 jadd = ',jadd
          write(iwr,'(a,10I3)') ' CONF FINAL: ',
     1          (elconfvaria(j,iconf), j=1,maxx)
        end if
        if ( idbg6 .eq. 1 ) 
     1    write(iwr,2233)
2233    format('<<',45('=='),'>>')
      end do
      if ( idbg6 .eq. 1 ) then
        do i=1,nconfv
          write(iwr,'(a,30I3)') ' Electronic Configurations: ',
     1          (elconfvaria(j,i), j=1,maxx)
        end do
      end if

      if ( aivb_optorb ) ncore_i = 0

c...  setting up numb value, which reffers to total number of electrons
c...  to check whether we have proper number of electrons
      do i=1,nconfv
        do j=i+1,nconfv
          if ( ilast(i) .ne. ilast(j) ) then
            write(iwr,'(a)') 'use different ilast(i) values'
            goto 1123
          end if
        end do
      end do
1123  continue

      do i=1,nconfv
        ilast(nconfv+i) = 0
        do j=1,maxx
          if ( ABS(elconfvaria(j,i)) .le. ilast(i) ) then
            ilast(nconfv+i) = ilast(nconfv+i) + 1
          end if
        end do
      end do

      do i=1,nconfv
        do j=i+1,nconfv
          if ( ilast(nconfv+i) .ne. ilast(nconfv+j) ) then
            call vberr('number of electrons in electronic configurations
     1 doesnt match - vbaivb code error. contact the developers')
          end if
        end do
      end do
      numb = ilast(nconfv+1)

      if ( idbg3 .eq. 1 ) then
        write(iwr,'(a,10I3)') ' ilast: ',(ilast(i), i=1,nconfv)
        write(iwr,'(a,I3)') ' numb = ',numb
      end if

      if ( numb .ne. maxx ) then
        call vberr(' unexpected number of electrons in configuration ')
      end if

      if ( aivb_diff .eqv. .true. ) then
        if ( readnr .gt. 1 ) then
          aivb_diffadd = aivb_diffadd + ilast(nconfv) - ncore_i
          if ( idbg3 .eq. 1 ) then
            write(iwr,'(a,I3)') ' ====> aivb_diff different orbitals: ',
     &                          aivb_diffadd
          end if
        else
          aivb_diffadd = 0
          if ( idbg3 .eq. 1 ) then
            write(iwr,'(a,I3)') ' ====> aivb_diff different orbitals: ',
     &                          aivb_diffadd
          end if
        end if
      else
        aivb_diffadd = 0
      end if

      DEALLOCATE(ilast, STAT=status)
      call deallocmem('ilast','gen_aivb_drv','vbaivb.m',status)

      ALLOCATE(kspinp(nconfv,numb), STAT=status)
      call allocmem('kspinp',nconfv*numb,'gen_aivb_drv',
     1'vbaivb.m',status)

      ALLOCATE(korder(mdet), STAT=status)
      call allocmem('korder',mdet,'gen_aivb_drv',
     1'vbaivb.m',status)


c
c...  *****************beginning of the determinant part******************
c...  *******************each configuration separately********************
c
c...  loop over the number of 'elconfvaria' configurations
c...  ..................................................................
c
      aimo_confpr(readnr) = 0
      aimo_nconfs = 0
      do iconf=1,nconfv

        if ( iconf .gt. 1 ) then
          if ( kconfdef(iconf,1) .ne. kconfdef(iconf-1,1) ) then
            ndetpc = 0
            aimo_confnr(readnr) = aimo_confnr(readnr) + 1
            if ( knconf(1,iconf,1) .gt. 1 ) then
              aimo_nconfs = 1
            else
              aimo_nconfs = 0
            end if
          end if
        else
          if ( knconf(1,iconf,1) .gt. 1 ) then
            aimo_nconfs = 1
          else
            aimo_nconfs = 0
          end if
          ndetpc = 0
          aimo_confnr(readnr) = 1
        end if

        if ( aimo_trconf(readnr,1,1) .gt. 0 ) then
          if ( idbg7 .eq. 1 ) then
            write(iwr,'(a,I3)') ' <=> truncation: ',
     1            aimo_trconf(readnr,1,1)
          end if
          do jtr=1,aimo_trconf(readnr,1,1)
            if ( idbg7 .eq. 1 ) then
              write(iwr,'(a,I3)') ' <===> jtr: ',jtr
              write(iwr,'(a,I3)') ' <===> aimo_trconf(readnr,2,jtr): '
     1                ,aimo_trconf(readnr,2,jtr)
              write(iwr,'(a,I3)') ' <===> aimo_confnr: ',
     1                            aimo_confnr(readnr)
            end if
            if ( aimo_trconf(readnr,2,jtr) .eq. 
     &           aimo_confnr(readnr) ) then
              if ( idbg7 .eq. 1 ) write(iwr,'(a)') ' <=====> true '
              if ( iconf .gt. 1 ) then
                if ( kconfdef(iconf,1) .ne. kconfdef(iconf-1,1) ) then
                  if ( idbg7 .eq. 1 ) then
                    write(iwr,*) ' CASE 1a '
                  end if
                  confnr = confnr + 1
c...  clear the 'korder' vector
                  call vclr(korder,1,mdet)
                  aimo_cpats = 1
                  aimo_lmc = 1
                  aimo_confpr(readnr) = aimo_confpr(readnr) + 1
c
c...  generate description
c
                  call gen_aivb_descr(readnr,confnr+confnrold,
     1            aimo_confnr(readnr),kcind,ncall,confnr)
c
c...
c
                else
                  if ( idbg7 .eq. 1 ) then
                    write(iwr,*) ' CASE 1b '
                  end if
                  aimo_cpats = aimo_cpats + 1
                end if
                goto 1212
              else
                if ( idbg7 .eq. 1 ) then
                  write(iwr,*) ' CASE 2 '
                end if
c...  clear the 'korder' vector
                call vclr(korder,1,mdet)
                confnr = confnr + 1
                aimo_cpats = 1
c
c...  generate description
c
                call gen_aivb_descr(readnr,confnr+confnrold,
     1          aimo_confnr(readnr),kcind,ncall,confnr)
c
c...
c
                aimo_lmc = 1
                aimo_confpr(readnr) = aimo_confpr(readnr) + 1
                goto 1212
              end if
            end if
          end do
          if ( idbg7 .eq. 1 ) then
            write(iwr,'(a)') ' <=====> false '
            write(iwr,'(a)') '                ---------------------'
            write(iwr,'(a)') '  '
          end if
          goto 1213
        else
          if ( iconf .gt. 1 ) then
            if ( kconfdef(iconf,1) .ne. kconfdef(iconf-1,1) ) then
              if ( idbg7 .eq. 1 ) then
                write(iwr,*) ' CASE 3a '
              end if
c...  clear the 'korder' vector
              call vclr(korder,1,mdet)
              confnr = confnr + 1
              aimo_cpats = 1
              aimo_lmc = 1
              aimo_confpr(readnr) = aimo_confpr(readnr) + 1
c
c...  generate description
c
              call gen_aivb_descr(readnr,confnr+confnrold,
     1        aimo_confnr(readnr),kcind,ncall,confnr)
c
c...
c
            else
              if ( idbg7 .eq. 1 ) then
                write(iwr,*) ' CASE 3b '
              end if
              aimo_cpats = aimo_cpats + 1
            end if
          else 
            if ( idbg7 .eq. 1 ) then
              write(iwr,*) ' CASE 4 '
            end if
c...  clear the 'korder' vector
            call vclr(korder,1,mdet)
            confnr = confnr + 1
            aimo_cpats =  1
c
c...  generate description
c
            call gen_aivb_descr(readnr,confnr+confnrold,
     1      aimo_confnr(readnr),kcind,ncall,confnr)
c
c...
c
            aimo_lmc = 1
            aimo_confpr(readnr) = aimo_confpr(readnr) + 1
          end if
        end if
1212    continue

        if ( idbg7 .eq. 1 ) then
          write(iwr,'(5(a,I3))') 
     1    ' confnr: ',confnr,
     2    '  aimo_confnr: ',aimo_confnr(readnr),
     3    '  _confpr: ',aimo_confpr(readnr),
     4    '  _cpats: ',aimo_cpats,
     5    '  _lmc: ',aimo_lmc
        end if

c
c...  ==============> generating determinants <==============
c
        if ( idbg7 .eq. 1 ) then
          write(iwr,2234)
2234      format(45('[]'),/,45('[]'))
          write(iwr,'(a,a,I3)') ' generating determinants for ',
     1    'configuration # ',aimo_confnr(readnr)
        end if

        if ( idbg7 .eq. 1 ) then
          write(iwr,'(a,I3)') 
     1    ' > 1 confs making up the at. state : ',aimo_nconfs
        end if

c...  counting up singly occ. orbitals, and
c...  setting up number of electrons belonging to each group
c...  nnd is the number of singly occ. electrons to be permuted
        do i=1,nat
          nsingat(i) = 0
        end do

        nnd = 0
        nalphapconf = 0
        nbothsp = 0
        do jatom=1,aimo_atcntr(readnr)
          aivb_nms(jatom) = 0.0d0
          kind = (jatom-1)*2+1 + 2
          kxind = kconfdef(iconf,kind)
          kyind = kconfdef(iconf,kind+1)
          klimit = aimo_confdef(kxind,kyind,1)
          if ( idbg7 .eq. 1 ) then
            print *,'-----------------------------------------------'
            write(iwr,'(a,I3)') ' jatom : ',jatom
            write(iwr,'(a,I3)') ' kind  : ',kind
            write(iwr,'(a,I3)') ' kxind : ',kxind
            write(iwr,'(a,I3)') ' kyind : ',kyind
            write(iwr,'(a,I3)') ' klimit: ',klimit
          end if

          do i=1,klimit
            itmp1 = aimo_confdef(kxind,kyind,(i+1))
            if ( idbg7 .eq. 1 ) then
              print *,'----------------------------------------'
              write(iwr,'(a,I3)') ' i    : ',i
              write(iwr,'(a,I3)') ' itmp1: ',itmp1
            end if
            if ( itmp1 .ge. 100 ) goto 1125
            if ( itmp1 .gt. 0 ) then
              aivb_nms(jatom) = aivb_nms(jatom) + 0.50d0
            else if ( itmp1 .lt. 0 ) then
              aivb_nms(jatom) = aivb_nms(jatom) - 0.50d0
            end if
1125        continue
          end do

          aimo_nnd(jatom) = 0
          do i=1,klimit
            itmp1 = ABS(aimo_confdef(kxind,kyind,(i+1)-1))
            itmp2 = ABS(aimo_confdef(kxind,kyind,(i+1)))
            itmp22 = aimo_confdef(kxind,kyind,(i+1))
            itmp3 = ABS(aimo_confdef(kxind,kyind,(i+1)+1))
            if ( itmp2 .ge. 100 ) goto 1126
            if ( itmp22 .gt. 0 ) then
              itmp22 = 1
            else
              itmp22 = 0
            end if
            if ( idbg7 .eq. 1 ) then
              print *,'----------------------------------------'
              write(iwr,'(a,I3)') ' i    : ',i
              write(iwr,'(a,I3)') ' itmp1: ',itmp1
              write(iwr,'(a,I3)') ' itmp2: ',itmp2
              write(iwr,'(a,I3)') ' itmp22: ',itmp22
              write(iwr,'(a,I3)') ' itmp3: ',itmp3
            end if
            if ( i .eq. 1 ) then
              if ( itmp2 .ne. itmp3 ) then
                nnd = nnd + 1
                aimo_nnd(jatom) = aimo_nnd(jatom) + 1
                aimo_elgrps(nnd,1) = jatom
                kspinp(iconf,nnd) = itmp22
                aimo_elgrps(nnd,2) = itmp22
                if ( itmp22 .eq. 1 )
     1          nalphapconf = nalphapconf + 1
              end if
            else if ( i .eq. klimit ) then
              if ( itmp2 .ne. itmp1 ) then
                aimo_nnd(jatom) = aimo_nnd(jatom) + 1
                nnd = nnd + 1
                aimo_elgrps(nnd,1) = jatom
                kspinp(iconf,nnd) = itmp22
                aimo_elgrps(nnd,2) = itmp22
                if ( itmp22 .eq. 1 )
     1          nalphapconf = nalphapconf + 1
              end if
            else
              if ( itmp2 .ne. itmp3 ) then
                if ( itmp2 .ne. itmp1 ) then
                  aimo_nnd(jatom) = aimo_nnd(jatom) + 1
                  nnd = nnd + 1
                  aimo_elgrps(nnd,1) = jatom
                  kspinp(iconf,nnd) = itmp22
                  aimo_elgrps(nnd,2) = itmp22
                  if ( itmp22 .eq. 1 )
     1            nalphapconf = nalphapconf + 1
                end if
              end if
            end if
1126        continue
          end do
          if ( idbg7 .eq. 1 ) then
            write(iwr,'(a,I3,a,F4.1)') ' _nms atom ',jatom,': ',
     1                             aivb_nms(jatom)
            write(iwr,'(2(a,i3))') ' nnd after atom ',jatom,': ',nnd
            write(iwr,'(a,I3)') ' nalphapconf: ',nalphapconf
          end if
c
c...  if there's only one singly occ. orbital in the atom
c...  it can has either an alpha or beta spin
c...  'kspinp = 2' indicates free choice for 'kterm' in subr. 'spinef'
c
          if ( aimo_nnd(jatom) .eq. 1 ) then
            kspinp(iconf,nnd) = 2
            aimo_elgrps(nnd,2) = 2
            nbothsp = nbothsp + 1
          end if

        end do
        if ( idbg7 .eq. 1 ) then
          write(iwr,'(a,I3)') ' nbothsp: ',nbothsp
          write(iwr,'(a,I3)') ' before nalphapconf: ',nalphapconf
        end if

c
c...  the number of alpha spins have to be adapted, if the number of the singly
c...  occupied orbitals equal the old number of alpha spins 
c...  (if it does, it means that we don't care about the spin projection)
c...  very rare situation
c
        if ( nbothsp .gt. 2 ) then
          nalphapconf = nalphapconf - nbothsp/2
        else
          nalphapconf = nalphapconf - (nbothsp+1)/2
        end if

        if ( idbg7 .eq. 1 ) then
          write(iwr,'(a,I3)') ' nalphapconf: ',nalphapconf
          write(iwr,'(a,30F4.1)') ' _nms : ',(aivb_nms(i), i=1,
     1    aimo_atcntr(readnr))
          write(iwr,'(a,30I3)') ' _nnd: ',(aimo_nnd(i), i=1,
     1    aimo_atcntr(readnr))
        end if

        ALLOCATE(kindx(nnd,2), STAT=status)
        call allocmem('kindx',nnd*2,'gen_aivb_drv','vbaivb.m',
     1  status)

        do i=1,nelec
          aimo_atord(i) = ABS(elconfvaria(i+ncore_i*2,iconf)) - ncore_i
        end do

c...  'kindx' contains indexes of singly occ. orbitals, from 'elconfvaria', for each configuration
        ising = 0
        do j=1,nelec
          ij = 0
          itmp1 = aimo_atord(j-1)
          itmp2 = aimo_atord(j)
          itmp3 = aimo_atord(j+1)
          if ( j .eq. 1 ) then
            if ( itmp2 .eq. itmp3 ) then
              ij = ij + 1
            end if
          else if ( j .eq. nelec ) then
            if ( itmp2 .eq. itmp1 ) then
              ij = ij + 1
            end if
          else
            if ( itmp2 .ne. itmp3 ) then
              ij = ij + 1
            end if
            if ( itmp2 .ne. itmp1 ) then
              ij = ij + 1
            end if
          end if
          if ( ij .eq. 2 .or. ij .eq. 0 ) then
            ising = ising + 1
            kindx(ising,1) = j
            kindx(ising,2) = aimo_atord(j)
          end if
        end do

c...  'aimo_atord' is now in a mixed form
        if ( idbg7 .eq. 1 ) then
          write(iwr,'(a,100I3)') ' _atord : ',
     1         (aimo_atord(j),j=1,nelec)
          write(iwr,'(a,10I3)') ' spin pattern : ',
     1        (kspinp(iconf,j),j=1,nnd)
          write(iwr,'(a,30I3)')  ' _elgrps 2: ',(aimo_elgrps(i,2), 
     1                            i=1,nnd)
          write(iwr,'(a,100I3)') ' indexes: ',
     1         (kindx(j,1),j=1,nnd)
          write(iwr,'(a,100I3)') ' values : ',
     1         (kindx(j,2),j=1,nnd)
        end if

c...  setting up 'aimo_atord' (given configurtion) in an <singly occ><doubly occ> form
        do i=1,nelec
          if ( idbg7 .eq. 1 ) print *,'i loop: ',i
          do j=1,nnd
            itmp1 = kindx(j,1)
            if ( idbg7 .eq. 1 ) then
              print *,' j loop: ',j
              print *,' itmp1 : ',itmp1
            end if
            if ( i .eq. itmp1 ) then
              if ( idbg7 .eq. 1 ) then
                print *,'  i: ',i,'  itmp1: ',itmp1
                print *,'  j: ',j,'  itmp1-1: ',itmp1-1
              end if
              do k=itmp1,j+1,-1
                if ( idbg7 .eq. 1 ) then
                  print *,'    k: ',k
                  print *,'    k-1: ',k-1
                  print *,'    k-1 val: ',aimo_atord(k-1)
                end if
                itmp2 = aimo_atord(k-1)
                aimo_atord(k) = itmp2
                if ( idbg7 .eq. 1 )
     1          print *,'    k val: ',aimo_atord(k)
              end do
              aimo_atord(j) = kindx(j,2)
              if ( idbg7 .eq. 1 )
     1        print *,'  j val: ',aimo_atord(j)
            end if
          end do
        end do
        aivb_ndocc = (nelec - nnd)/2

        DEALLOCATE(kindx, STAT=status)
        call deallocmem('kindx','gen_aivb_drv','vbaivb.m',status)

        if ( idbg7 .eq. 1 ) then
          write(iwr,'(a,I3)') ' _detnr: ',aimo_detnr
          write(iwr,'(a,I3)') ' aivb_ndocc: ',aivb_ndocc
        end if

c...
c...  checking the existing determinants in the configuration, their spin patterns
c...
        if ( aimo_cpats .gt. 1 .and. aimo_cpats .le. aimo_detnr ) then
          if ( idbg7 .eq. 1 ) then
            write(iwr,'(a,a)') ' checking the available determinants ',
     1      'in the configurations '
          end if

c          do i=1,aimo_detnr
c...  checking if we have the same orbitals or not
            itmp1 = 0
            do j=1,maxx-ncore_i*2
              do k=1,maxx-ncore_i*2
                if ( aimo_atord(j) .eq. 
     1               deter(k,detnrold-aimo_detnr+aimo_lmc) ) then
                  itmp1 = itmp1 + 1
                  goto 2022
                end if
              end do
2022          continue
            end do
c...  we do not - need new determinants
            if ( itmp1 .ne. maxx-ncore_i*2 ) then
              if ( idbg7 .eq. 1 ) then
                write(iwr,'(a)') 
     1          ' the given configuration have some different orbitals.'
                write(iwr,'(a)') 
     2          ' need to generate new determinants.'
              end if
              goto 2021
            end if
c          end do

          ALLOCATE(kspdet(aimo_detnr,nnd), STAT=status)
          call allocmem('kspdet',aimo_detnr*nnd,'gen_aivb_drv',
     1    'vbaivb.m',status)

          ALLOCATE(koccn(maxx-ncore_i*2), STAT=status)
          call allocmem('koccn',maxx-ncore_i*2,'gen_aivb_drv',
     1    'vbaivb.m',status)
          ALLOCATE(ksptmp(aimo_detnr,nnd), STAT=status)
          call allocmem('ksptmp',aimo_detnr*nnd,'gen_aivb_drv',
     1    'vbaivb.m',status)

          print *,' //=======================================\\ '
          print *,'   ITMP1 in DETERMINANT CHECK needs fixing   '
          print *,' \\=======================================// '
c...  the itmp1 needs fixing. wrong use of nalphaconf in here. 
c...  I guess the msocc/iconfalpha could be re-calculated here as well for this use only
          itmp1 = nalphapconf + aivb_ndocc - nbothsp/2

c...  setting up the spin patterns of the available determinants
          if ( idbg7 .eq. 1 ) then
            print *,' itmp1: ',itmp1
            print *,' maxx: ',maxx
            print *,' ncore_i: ',ncore_i
            print *,' maxx-ncore_i*2: ',maxx-ncore_i*2
          end if
          do j=aimo_lmc,aimo_detnr
            call izero(maxx-ncore_i*2,koccn,1) 
            do i=1,maxx-ncore_i*2
              itmp2 = deter(i,j+detnrold-aimo_detnr)
              if ( idbg7 .eq. 1 )
     1        print *,' itmp2: ',itmp2
              if ( i .gt. itmp1 ) then
                koccn(itmp2) =  koccn(itmp2) + 2
              else if ( i .le. itmp1 ) then
                koccn(itmp2) =  koccn(itmp2) + 1
              end if
            end do

            itmp3 = 0
            do i=1,maxx-ncore_i*2
              if ( koccn(i) .eq. 1 ) then
                itmp3 = itmp3 + 1
                if ( kspinp(iconf,itmp3) .eq. 2 ) then
                  kspdet(j,itmp3) = 2
                else
                  kspdet(j,itmp3) = 1
                end if
              else if ( koccn(i) .eq. 2 ) then
                itmp3 = itmp3 + 1
                if ( kspinp(iconf,itmp3) .eq. 2 ) then
                  kspdet(j,itmp3) = 2
                else
                  kspdet(j,itmp3) = 0
                end if
              end if
            end do

            if ( idbg7 .eq. 1 ) then
              write(iwr,'(a,20I3)') ' koccn: ',
     2        ( koccn(i), i=1,maxx-ncore_i*2 )
              write(iwr,'(a,20I3)') ' kspdet: ',
     2        ( kspdet(j,i), i=1,nnd )
              print *,'------------------------------------------------'
            end if
          end do

          do j=aimo_lmc,aimo_detnr
            itmp1 = 0
            do i=1,nnd
              if ( idbg7 .eq. 1 ) then
                write(iwr,'(2(a,I3))') 
     1          ' kspinp: ',kspinp(iconf,i),'     kspdet: ',kspdet(j,i)
              end if
              if ( kspinp(iconf,i) .eq. kspdet(j,i) ) itmp1 = itmp1 + 1
            end do
            if ( idbg7 .eq. 1 ) then
              write(iwr,'(5(a,I3))')
     1        ' itmp1: ',itmp1,
     2        ' nnd: ',nnd,
     3        ' _cpats: ',aimo_cpats,
     4        ' _lmc: ',aimo_lmc,
     5        ' j: ',j
            end if
            if ( itmp1 .eq. nnd ) then
              if ( idbg7 .eq. 1 ) then
                write(iwr,'(a,I3,a,I3)') ' determinant nr ',j,
     1          ' already defines configuration ',aimo_cpats
              end if

              DEALLOCATE(kspdet, STAT=status)
              call deallocmem('kspdet','gen_aivb_drv','vbaivb.m',status)

              DEALLOCATE(koccn, STAT=status)
              call deallocmem('koccn','gen_aivb_drv','vbaivb.m',status)

              DEALLOCATE(ksptmp, STAT=status)
              call deallocmem('ksptmp','gen_aivb_drv','vbaivb.m',status)

              goto 2030
            end if
            if ( idbg7 .eq. 1 ) 
     1      print *,'--------------------------------------------------'
          end do

          DEALLOCATE(kspdet, STAT=status)
          call deallocmem('kspdet','gen_aivb_drv','vbaivb.m',status)

          DEALLOCATE(koccn, STAT=status)
          call deallocmem('koccn','gen_aivb_drv','vbaivb.m',status)

          DEALLOCATE(ksptmp, STAT=status)
          call deallocmem('ksptmp','gen_aivb_drv','vbaivb.m',status)

2021      continue
        end if

        aimo_lmc = aimo_cpats

c...  saving old parameters, if they were ever set
        oldmsocc  = msocc
        oldmstrco = mstrco
        oldmdetco = mdetco
        oldmdetst = mdetst
        oldmperm  = mperm
        oldmstruc = mstruc
        oldmsp    = msp

c... calculating parameters required for 'spinef' subroutine
        msocc = nnd
        mpair = msocc/2
        msp = 0

        if ( idbg7 .eq. 1 ) then
          print *,' 0 msp ',msp
          print *,' 1 msocc ',msocc
          print *,' 2 pair ',mpair
        end if

        iconfspin = (msocc - mpair*2)/2.0d0
        if ( idbg7 .eq. 1 ) 
     1  print *,' 3 spin ',iconfspin
        iconfmult = 2.0d0*iconfspin + 1.0d0
        if ( idbg7 .eq. 1 ) 
     1  print *,' 4 mult ',iconfmult
        iconfnalpha = (msocc + iconfmult - 1)/2
        if ( idbg7 .eq. 1 ) 
     1  print *,' 5 nalpha ',iconfnalpha

        n2s = msocc/2.0d0 - iconfspin
        if ( idbg7 .eq. 1 ) 
     1  print *,' 6 n2s ',n2s
        if ( n2s .eq. 0 ) then
          mstrco = 1
        else
          new1 = newtont(msocc,n2s)
          n2s1 = n2s - 1
          new2 = newtont(msocc,n2s1)
          mstrco = new1 - new2
        end if
        if ( idbg7 .eq. 1 ) 
     1  print *,' 7 mstrco ',mstrco

        mdetst = 2**((msocc-iconfmult+1)/2)
        if ( idbg7 .eq. 1 ) 
     1  print *,' 8 mdetst ',mdetst

        mdetco = mstrco*mdetst
        if ( idbg7 .eq. 1 ) 
     1  print *,' 9 mdetco ',mdetco

        n2s2 = n2s/2
        if ( idbg7 .eq. 1 ) 
     1  print *,' 10 n2s2 ',n2s2
        if ( n2s2 .gt. 0 ) then
          mperm = newtont(n2s,n2s2)
        else
          mperm = 1
        end if
        if ( idbg7 .eq. 1 ) 
     1  print *,' 11 mperm ',mperm

        mstruc = mstrco
        if ( idbg7 .eq. 1 ) 
     1  print *,' 12 mstruc ',mstruc

        if ( idbg7 .eq. 1 ) then
          write(iwr,'(6(a,I4))') 
     1          ' msocc: ',msocc,'  mstrco: ',mstrco,
     2          '  mdetco: ',mdetco,'  mdetst: ',mdetst,
     3          '  mperm: ',mperm,'  mstruc: ',mstruc
          write(iwr,'(a,30I3)')  ' _atord : ',(aimo_atord(i), i=1,nelec)
          write(iwr,'(a,30I3)')  ' _elgrps 1: ',(aimo_elgrps(i,1), 
     1                            i=1,nnd)
          write(iwr,'(a,30I3)')  ' spinp : ',(kspinp(iconf,i), 
     1                            i=1,nnd)
        end if

        allout = .false.
        adbg9 = idbg9
        adbg10 = idbg10
        if ( idbg11 .eq. 1 ) allout = .true.
        projec = .false.
        
c...  setting up required memory for 'spinef' and 'schmic' subroutines
        mscr = igmem_max_memory()
        mscr = mscr/3
        mscr = max(mscr,25000)
        mscr = min(mscr,2000000)

        annwdet = 0
        arumer = .false.
        anconf = 1

c...  allocate memory
        ALLOCATE(kdeter(nelec,mdetco), STAT=status)
        call allocmem('kdeter',nelec*mdetco,'gen_aivb_drv','vbaivb.m',
     1  status)

        ALLOCATE(klogic(mstrco*mdetco), STAT=status)
        call allocmem('klogic',mstrco*mdetco,'gen_aivb_drv','vbaivb.m',
     1  status)

        ALLOCATE(kstruc(mdetst,mstrco), STAT=status)
        call allocmem('kstruc',mdetst*mstrco,'gen_aivb_drv','vbaivb.m',
     1  status)

        ALLOCATE(kcoeff(mdetst,mstrco), STAT=status)
        call allocmem('kcoeff',mdetst*mstrco,'gen_aivb_drv','vbaivb.m',
     1  status)

        ALLOCATE(kscr(mscr), STAT=status)
        call allocmem('kscr',mscr,'gen_aivb_drv','vbaivb.m',
     1  status)

        ALLOCATE(krumer(maxspvb*msocc), STAT=status)
        call allocmem('krumer',maxspvb*msocc,'gen_aivb_drv','vbaivb.m',
     1  status)

        ALLOCATE(kdetps(mstruc), STAT=status)
        call allocmem('kdetps',mstruc,'gen_aivb_drv','vbaivb.m',
     1  status)

        ALLOCATE(kcoefff(mstruc*mdetco), STAT=status)
        call allocmem('kcoefff',mstruc*mdetco,'gen_aivb_drv','vbaivb.m',
     1  status)

        ALLOCATE(kidetps(mstruc*mdetco), STAT=status)
        call allocmem('kidetps',mstruc*mdetco,'gen_aivb_drv','vbaivb.m',
     1  status)

        if ( idbg7 .eq. 1 )
     1  write(iwr,'(a,I10)') ' memory for spinef: ',mscr

c...  generate leading terms and all the determinants
        call spinef(aimo_atord,anconf,kdeter,klogic,
     1              annwdet,kstruc,
     2              kcoeff,ndetst,nstruc,allout,projec,
     3              kscr(1),mscr,mdetst,nelec,msp,krumer)
        
        if ( idbg7 .eq. 1 ) 
     1  print *,'leading term of inteterest: ',aimo_lterm(aivb_nlterm)

        if ( idbg7 .eq. 1 ) then
          do i=1,nstruc
            write(iwr,'(a,20I3)') ' _lterms: ',
     1      (aimo_lterms(i,j), j=1,nnd)
          end do
          write(iwr,'(a,20F8.4)') 
     1          ' kcoeff: ',(kcoeff(j,i),j=1,ndetst)
        end if

        anspath = 1
        aispath = aimo_lterm(aivb_nlterm)

c...  Schmidt-orthogonalize the rumer diagrams to get BD
        if ( nnd .gt. 0 ) then
        call schmic(kcoeff,mdetst,ndetst,nstruc,
     1              kstruc,kscr,annwdet,
     2              kcoefff,kdetps,
     3              kidetps,ncc,
     4              anspath,aispath,
     5              arumer,kdeter,nelec,klogic,
     6              krumer)
        end if

        if ( nnd .eq. 0 ) then
          annwdet = 1
          kcoefff(1) = 1.0d0
          kdetps(1) = 1
          kidetps(1) = 1
        end if

        if ( idbg7 .eq. 1 )  then
          write(iwr,'(2(a,I3))')
     1      ' nnwdet/kntdet: ',annwdet,
     2      '  nstruc/knstru: ',nstruc
          write(iwr,'(a,20I3)') ' kdetps: ',(kdetps(i), i=1,nstruc)
          write(iwr,'(a,1000I3)') 
     1          ' kidetps: ',(kidetps(i), i=1,ncc)
          write(iwr,'(a,20F8.4)') 
     1          ' kcoefff: ',(kcoefff(i), i=1,ncc)
        end if

        detnr = kdetps(1)
        do i=1,kdetps(1)
          itmp1 = kidetps(i)
          do j=1,nelec
            deter(j,detnrold+i) = kdeter(j,itmp1) + aivb_diffadd
          end do
          coefff(detnrold+i) = kcoefff(i)
          korder(aimo_detstruc(readnr,2,confnr)+i) = aimo_cpats

          if ( idbg7 .eq. 1 ) then
           write(iwr,'(a,I3,I3,I3)')
     1     'korder: ',aimo_detstruc(readnr,2,confnr),i,
     2     korder(aimo_detstruc(readnr,2,confnr)+i)
          end if
        end do
        aimo_detstruc(readnr,2,confnr) = aimo_detstruc(readnr,2,confnr)
     &                                 + detnr
        aimo_detnr = aimo_detstruc(readnr,2,confnr)

        if ( detnr .eq. 1 ) ndetpc = ndetpc + 1

        if ( idbg7 .eq. 1 )
     &  write(iwr,'(a,I3)') ' ndetpc: ',ndetpc

        if ( idbg7 .eq. 1 ) then
          do i=1,detnr
            itmp1 = aimo_detstruc(readnr,2,confnr) - detnr
            write(iwr,'(a,I3,a,I3,a,F8.4,a,30I3)') 
     1      ' det. ',i,'  grp. ',korder(itmp1+i),': ',
     2      coefff(detnrold+i),'  ',(deter(j,i+detnrold), 
     3      j=1,maxx-ncore_i*2)
          end do
        end if

        detnrold = detnr + detnrold

        if ( aimo_cpats .eq. knconf(1,aimo_confnr(readnr),1) 
     1       .and. aimo_lmc .gt. 1 ) then
          if ( idbg7 .eq. 1 ) then
            write(iwr,'(a,I3,I3,I3)') 
     1      ' _cpats,knconf,_confnr ',aimo_cpats,
     2      knconf(1,aimo_confnr(readnr),1),aimo_confnr(readnr)
          end if
          anorm = 0.0d0
          itmp1 = detnrold - aimo_detstruc(readnr,2,confnr)
          if ( idbg7 .eq. 1 )
     1    write(iwr,'(a,I3)') 'itmp1 ',itmp1
          do i=1,aimo_detstruc(readnr,2,confnr)
            itmp2 = (iconf - knconf(1,aimo_confnr(readnr),1)) +korder(i)
            if ( idbg7 .eq. 1 ) then
              write(iwr,'(a,I3)') '  itmp2    ',itmp2
              write(iwr,'(a,I3)') '  korder   ',korder(i)
              write(iwr,'(a,I3)') '  kconfdef ',kconfdef(itmp2,2)
              write(iwr,'(a,I3)') '  coefff   ',coefff(itmp1+i)
            end if
            anorm = anorm + kconfdef(itmp2,2)*kconfdef(itmp2,2)*
     1      coefff(itmp1+i)*coefff(itmp1+i)
            if ( idbg7 .eq. 1 )
     1        print *,'----------------------------------'
          end do

          if ( idbg7 .eq. 1 )
     1      print *,' anorm 1 ',anorm
          anorm = 1.0d0/dsqrt(anorm)
          if ( idbg7 .eq. 1 )
     1      print *,' anorm 2 ',anorm

          do i=1,aimo_detstruc(readnr,2,confnr)
            itmp2 = (iconf - knconf(1,aimo_confnr(readnr),1)) +korder(i)
            coefff(itmp1+i) = ABS(coefff(itmp1+i))*anorm*
     2      kconfdef(itmp2,2)
          end do
        end if

        DEALLOCATE(kdeter, STAT=status)
        call deallocmem('kdeter','gen_aivb_drv','vbaivb.m',status)

        DEALLOCATE(klogic, STAT=status)
        call deallocmem('klogic','gen_aivb_drv','vbaivb.m',status)

        DEALLOCATE(kstruc, STAT=status)
        call deallocmem('kstruc','gen_aivb_drv','vbaivb.m',status)

        DEALLOCATE(kcoeff, STAT=status)
        call deallocmem('kcoeff','gen_aivb_drv','vbaivb.m',status)

        DEALLOCATE(kscr, STAT=status)
        call deallocmem('kscr','gen_aivb_drv','vbaivb.m',status)

        DEALLOCATE(krumer, STAT=status)
        call deallocmem('krumer','gen_aivb_drv','vbaivb.m',status)

        DEALLOCATE(kdetps, STAT=status)
        call deallocmem('kdetps','gen_aivb_drv','vbaivb.m',status)

        DEALLOCATE(kcoefff, STAT=status)
        call deallocmem('kcoefff','gen_aivb_drv','vbaivb.m',status)

        DEALLOCATE(kidetps, STAT=status)
        call deallocmem('kidetps','gen_aivb_drv','vbaivb.m',status)

2030    continue

        if ( idbg7 .eq. 1 ) then
          write(iwr,'(a,I3)') ' confnr: ',confnr
          write(iwr,'(a,I3)') ' detnr : ',detnr

          write(iwr,'(a,I3)') ' confnrold: ',confnrold
          write(iwr,'(a,I3)') ' detnrold : ',detnrold

          do i=1,confnr
            write(iwr,'(I3,a,I3)') 
     1      i,' _detstruc: ',aimo_detstruc(readnr,2,i)
          end do
        end if

1213    continue
      end do
c...  ..................................................................

      DEALLOCATE(korder, STAT=status)
      call deallocmem('korder','gen_aivb_drv','vbaivb.m',status)

      DEALLOCATE(kspinp, STAT=status)
      call deallocmem('kspinp','gen_aivb_drv','vbaivb.m',status)

      DEALLOCATE(kconfdef, STAT=status)
      call deallocmem('kconfdef','gen_aivb_drv','vbaivb.m',status)

      DEALLOCATE(elconfvaria, STAT=status)
      call deallocmem('elconfvaria','gen_aivb_drv','vbaivb.m',status)

      DEALLOCATE(knconf, STAT=status)
      call deallocmem('knconf','gen_aivb_drv','vbaivb.m',status)

      DEALLOCATE(kcind, STAT=status)
      call deallocmem('kcind','gen_aivb_drv','vbaivb.m',status)

      confnr = confnrold + confnr
      detnr = detnrold

      if ( idbg3 .eq. 1 ) then
        print *,' '
        write(iwr,'(a,I3)') ' confnr: ',confnr
        write(iwr,'(a,I3)') ' detnr : ',detnr
        write(iwr,'(a,I3)') ' conf. per readnr: ',aimo_confpr(readnr)
        do i=1,detnr
          write(iwr,'(a,I3,a,F8.4,a,30I3)') 
     1    ' det. ',i,': ',
     2    coefff(i),'  ',(deter(j,i), j=1,maxx-ncore_i*2)
        end do
      end if

c...  bringing back the original values, whatever they were
      msocc  = oldmsocc
      mstrco = oldmstrco
      mdetco = oldmdetco
      mdetst = oldmdetst
      mperm  = oldmperm
      mstruc = oldmstruc
      msp    = oldmsp

      if (idbg3 .eq. 1) 
     1 print *,'     ==> GENERATAE AIVB DRIVER subroutine - end'
cmarcin      stop 'end of gen aivb drv'
      
      return
      end

c*****************************************************************************************
c...  prepare text descriptions for generated configurations
c
      subroutine gen_aivb_descr(readnr,confnr,nconf,kcind,ncall,iconf)
c
      implicit none
c
INCLUDE(common/turtleparam)
INCLUDE(common/aivb)
INCLUDE(../m4/common/iofile)
c
      character(len=32) aimotmpdescr,aimotmpdescr2,aimotmpdescr3
      character(len=32) atmpd(2)
      character(len=10) ctmp1
c
      integer :: ist1,iend1,ist2,iend2,jind,ic,jjx,jjy,jjz
      integer :: jx,jy,jz,icnt,igrp,ictmp1,jj,iconf,ii
      integer :: ist3,iend3,ist4,iend4,ist5,iend5
      integer :: jatom,idbg13,icind,ncall,status,i,j
      integer :: readnr,nconf,confnr,kcind(ncall,*)
c
      integer, allocatable, dimension(:,:) :: kind
c

      idbg13 = aimo_debugpr(13)

      ALLOCATE(kind(aimo_atcntr(readnr),3),
     1         STAT=status)
      call allocmem('kind',aimo_atcntr(readnr)*3,'gen_aivb_descr',
     1 'vbaivb.m',status)

      if ( idbg13 .eq. 1 ) then
        write(iwr,'(a,I3)') ' readnr: ',readnr
        write(iwr,'(a,I3)') ' confnr: ',confnr
        write(iwr,'(a,I3)') ' nconf: ',nconf
        write(iwr,'(a,I3)') ' ncall: ',ncall
      end if
      call izero(aimo_atcntr(readnr)*3,kind,1)
c
c...  first, gather an informations about the atoms
c
      igrp = 1
      icnt = 1
      do jatom=1,aimo_atcntr(readnr)
        if ( kind(jatom,3) .eq. 1 ) goto 1006
        kind(jatom,1) = igrp
        kind(igrp,2) = 1
        kind(jatom,3) = 1
c...  x index for _conflab
        jjx = aimo_atindex(jatom,2)
c...  y index for _conflab
        jjy = aimo_atindex(jatom,3)
c...  z index for _conflab
        jjz = kcind(nconf,jatom)
        do i=jatom+1,aimo_atcntr(readnr)
c...  x index for _conflab
          jx = aimo_atindex(i,2)
c...  y index for _conflab
          jy = aimo_atindex(i,3)
c...  z index for _conflab
          jz = kcind(nconf,i)
          if ( idbg13 .eq. 1 ) then
            write(iwr,'(2(a,I3))')
     1      ' jatom: ',jatom,' vs  iatom: ',i
            write(iwr,'(2(a,I3))')
     1      ' jjx: ',jjx,'  jx: ',jx
            write(iwr,'(2(a,I3))')
     1      ' jjy: ',jjy,'  jy: ',jy
            write(iwr,'(2(a,I3))')
     1      ' jjz: ',jjz,'  jz: ',jz
          end if
          if ( jjx .eq. jx .and. jjy .eq. jy .and. jjz .eq. jz ) then
            kind(i,1) = igrp
            kind(igrp,2) = kind(igrp,2) + 1
            kind(i,3) = 1
            icnt = icnt + 1
            if ( idbg13 .eq. 1 ) then
              write(iwr,'(2(a,I3))')
     1        ' kind,1: ',kind(i,1),'  iatom: ',i
              write(iwr,'(a,I3,a,I3)')
     1        ' kind,2: ',kind(i,2),'   ',i
              write(iwr,'(a,I3,a,I3)')
     1        ' kind,3: ',kind(i,3),'   ',i
              write(iwr,'(a,I3)') ' icnt: ',icnt
            end if
          end if
1009      continue
          if ( idbg13 .eq. 1 ) then
            write(iwr,'(a,30I4)') ' kind,1: ',
     1      (kind(ii,1), ii=1,aimo_atcntr(readnr))

            write(iwr,'(a,30I4)') ' kind,2: ',
     1      (kind(ii,2), ii=1,aimo_atcntr(readnr))

            write(iwr,'(a,30I4)') ' kind,3: ',
     1      (kind(ii,3), ii=1,aimo_atcntr(readnr))
            write(iwr,1007)
          end if
        end do
        if ( icnt .eq. aimo_atcntr(readnr) ) then
          goto 1010
        end if
        igrp = igrp + 1
1006    continue
        icnt = icnt + 1
        if ( idbg13 .eq. 1 ) then
          write(iwr,'(a,I3)') ' icnt: ',icnt
          write(iwr,'(a,I3)') ' igrp: ',igrp
          write(iwr,1008)
        end if
      end do
1007  format(30('-'))
1008  format(50('-'))
1010  continue

      if ( idbg13 .eq. 1 ) then
        write(iwr,'(a,30I4)') ' kind,1: ',
     1  (kind(i,1), i=1,aimo_atcntr(readnr))

        write(iwr,'(a,30I4)') ' kind,2: ',
     1  (kind(i,2), i=1,igrp)

        write(iwr,'(a,30I4)') ' kind,3: ',
     1  (kind(i,3), i=1,igrp)
      end if
      
c...  loop over number of group atoms
c      do jatom = 1,aimo_atcntr(readnr)
      jj = 1
      do jatom = 1,igrp
        i = kind(jatom,1)
        if ( idbg13 .eq. 1 ) then
          write(iwr,'(a,I3)') ' jatom: ',jatom
        end if

        if ( iconf .eq. 1 ) then
          ctmp1 = aivb_atnmrdin(jatom,readnr)
          ictmp1 = ichar(ctmp1(1:1)) - 32
          call stripblanks(ist1,iend1,ctmp1)
          ctmp1 = char(ictmp1)
          aivb_atnmrdin(jatom,readnr) = 
     &    ctmp1(ist1:iend1)//aivb_zchr(jatom)//'('

          do j=1,kind(jatom,2)
            call stripblanks(ist2,iend2,aivb_atnmrdin(jatom,readnr))
            aivb_atnmrdin(jatom,readnr) = 
     &      aivb_atnmrdin(jatom,readnr)(ist2:iend2)//
     &      char(jj+48)
            if ( j .lt. kind(jatom,2) ) then
              call stripblanks(ist3,iend3,aivb_atnmrdin(jatom,readnr))
              aivb_atnmrdin(jatom,readnr) = 
     &        aivb_atnmrdin(jatom,readnr)(ist3:iend3)
     &        //','
            end if
            jj = jj + 1
          end do

          call stripblanks(ist2,iend2,aivb_atnmrdin(jatom,readnr))
          aivb_atnmrdin(jatom,readnr) = 
     &    aivb_atnmrdin(jatom,readnr)(ist2:iend2)//')'

          if ( idbg13 .eq. 1 ) then
            write(iwr,'(a,30a,a)') ' atnmrdin  : ',
     &      aivb_atnmrdin(jatom,readnr)
          end if
        end if


        jind = (i-1)*2+1 + 2
c...  x index for _conflab
        jjx = aimo_atindex(jatom,2)
c...  y index for _conflab
        jjy = aimo_atindex(jatom,3)
c...  z index for _conflab
        jjz = kcind(nconf,jatom)

        if ( idbg13 .eq. 1 ) then 
          write(iwr,'(a,I3)') ' 2 atom: ',jatom
          write(iwr,'(a,I3)') ' 3 jjx : ',jjx
          write(iwr,'(a,I3)') ' 4 jjy : ',jjy
          write(iwr,'(a,I3)') ' 44 jjz: ',jjz
          write(iwr,'(a,a,a)') ' aivb_conflab  :',
     1          aivb_conflab(jjx,jjy,jjz)(1:32),'|'
          write(iwr,'(a,a,a)') ' aivb_statelab  :',
     1          aivb_statelab(jjx,jjy,1)(1:20),'|'
          write(iwr,'(a,a,a)') ' aivb_statelab+1:',
     1          aivb_statelab(jjx,jjy,jjz+1)(1:20),'|'
        end if
        call stripblanks(ist1,iend1,aivb_conflab(jjx,jjy,jjz))
        call stripblanks(ist3,iend3,aivb_statelab(jjx,jjy,1))
        call stripblanks(ist4,iend4,aivb_statelab(jjx,jjy,jjz+1))
        if ( idbg13 .eq. 1 ) then
          write(iwr,'(a,I3,a,I3)') ' 6 ist1:',ist1,' iend1: ',iend1
        end if
        aimotmpdescr = aivb_conflab(jjx,jjy,jjz)
        aimotmpdescr2 = aivb_statelab(jjx,jjy,1)
        aimotmpdescr3 = aivb_statelab(jjx,jjy,jjz+1)
        if ( idbg13 .eq. 1 ) then
          write(iwr,'(a)') '7'
          write(iwr,'(a,a,a)') 'TMPDESCR     : ',aimotmpdescr,'END '
          write(iwr,'(a,a,a)') 'TMPDESCR TRIM: ',
     1                         aimotmpdescr(ist1:iend1),'END '
          write(iwr,'(a,a,a)') 'TMPDESCR 2   : ',aimotmpdescr2,'END '
          write(iwr,'(a,a,a)') 'TMPDESCR 2 TRIM: ',
     1                         aimotmpdescr2(ist3:iend3),'END '
          write(iwr,'(a,a,a)') 'TMPDESCR 3   : ',aimotmpdescr3,'END '
          write(iwr,'(a,a,a)') 'TMPDESCR 3 TRIM: ',
     1                         aimotmpdescr3(ist4:iend4),'END '
          write(iwr,'(a)') '8'
        end if
c...  first atom, no addidtion of things before
        if ( jatom .eq. 1 ) then
          call stripblanks(ist2,iend2,aivb_atnmrdin(jatom,readnr))
          aivb_confdescr(confnr,1) =
     &    aivb_atnmrdin(jatom,readnr)(ist2:iend2)//': '//
     &    aimotmpdescr

          aivb_confdescr(confnr,3) = 
     &    aivb_atnmrdin(jatom,readnr)(ist2:iend2)//': '//
     &    aimotmpdescr2(ist3:iend3)//' '//
     &    aimotmpdescr3(ist4:iend4)

c...  next atom, adding things before
        else
          call stripblanks(ist1,iend1,aivb_confdescr(confnr,1))
          call stripblanks(ist5,iend5,aivb_confdescr(confnr,3))
          if ( idbg13 .eq. 1 ) then 
            write(iwr,'(2(a,I3))') 'ist1: ',ist1,'  iend1: ',iend1
          end if
          call stripblanks(ist2,iend2,aivb_atnmrdin(jatom,readnr))
          aivb_confdescr(confnr,1) = 
     &    aivb_confdescr(confnr,1)(ist1:iend1)
     &    //'   '//
     &    aivb_atnmrdin(jatom,readnr)(ist2:iend2)//': '//
     &    aimotmpdescr

          aivb_confdescr(confnr,3) = 
     &    aivb_confdescr(confnr,3)(ist5:iend5)
     &    //'   '//
     &    aivb_atnmrdin(jatom,readnr)(ist2:iend2)//': '//
     &    aimotmpdescr2(ist3:iend3)//' '//
     &    aimotmpdescr3(ist4:iend4)
        end if
        if ( idbg13 .eq. 1 ) then
          write(iwr,'(a)' ) '====> '
          write(iwr,'(a,a,a)') 'aivb_confdescr(,1): ',
     1          aivb_confdescr(confnr,1),' END '

          write(iwr,'(a,a,a)') 'aivb_confdescr(,3): ',
     1          aivb_confdescr(confnr,3),' END '
          write(iwr,'(a)' ) '<==== '
          write(iwr,'(a)') '---------------------------------------'
        end if
      end do

      DEALLOCATE(kind, STAT=status)
      call deallocmem('kind','gen_aivb_descr','vbaivb.m',status)

      return
      end

c*****************************************************************************************
c...  main driver for generating all possible couplings between atomic state's
c...  variations of each atom
c
      subroutine gen_aivbconf_drv(confs,nn,nog,nall,confnr,grps)
c
      implicit none
c
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/runlab)
c
INCLUDE(common/turtleparam)
INCLUDE(common/aivb)
c
      integer :: confs(nall,nog),idbg8,i,j,grps(nn)
      integer :: nn,nall,confnr,nog,ifault
      character(len=8) :: mode
c
      idbg8 = aimo_debugpr(8)
      mode = 'aimo_c  '
      if ( idbg8 .gt. 0 ) mode = 'aimo_c_d'
c
      if ( idbg8 .gt. 0 ) then
        print *,'     ==> GEN AIVBCONF DRV subroutine - begin'
        write(iwr,'(1X,a,I3,a,I3,a,I3,a,I3,a,a)') 
     1        '# diff. states: ',nn,'  # atoms: ',nog,
     2        '  # all coupling between diff. states: ',nall,
     3        '  confnr: ',confnr,
     4        '  mode: ',mode
        write(iwr,'(a,20I3)') ' grps: ',(grps(i), i=1,nn)
      end if

      call gencomb(nn,nog,confs,confnr,ifault,nall,nog,mode,
     &             grps)

      if ( idbg8 .gt. 1 ) then
        write(iwr,'(a)') ' possible couplings between atomic states '
        write(iwr,'(a,I3)') ' # possible couplings: ',confnr
        do i=1,confnr
          write(*,'(I3,a,30I3)') i,'   ',(confs(i,j), j=1,nog)
        end do
      end if

      if ( idbg8 .gt. 0 ) 
     &  print *,'     ==> GEN AIVBCONF DRV subroutine - end'

      return
      end

c*****************************************************************************************
c...  calculate x!
c
      integer function ifact(x)
      implicit none
c
      integer :: x,temp,i,j

      temp = 1
      do i=2,x
        temp = temp*i
      end do

      ifact = temp
      return
      end

c*****************************************************************************************
      block data aivbdata
      implicit none
c
INCLUDE(common/turtleparam)
INCLUDE(common/aivb)
c
      integer :: i
c
      data (aivb_multname(i), i=1,6) /'1','2','3','4',
     &'5','    '/
      data (aivb_atstname(i), i=1,4) /'s','p','d',' '/
c
      data  aivb_allatst(1,1)         /'2p'/
      data (aivb_allatst(2,i), i=1,4) /'3p','1d',
     &                                 '1s','5s'/
      data (aivb_allatst(3,i), i=1,3) /'4s','2d','2p'/
      data (aivb_allatst(4,i), i=1,3) /'3p','1d','1s'/
      data  aivb_allatst(5,1)         /'2p'/
      data  aivb_allatst(6,1)         /'1s'/
      data  aivb_allatst(7,1)         /'2s'/
      data  aivb_allatst(8,1)         /'1s'/
      data  aivb_allatst(9,1)         /'      '/
c
c     limits tables
      data (aimo_atstcom(i,1), i=1,9)   /1,4,3,3,1,1,1,1,1/ ! number of atomic states per 's' or 'p' # of electrons

cc    1 p electrons
      data  aimo_atstcom(1,2)           /3/

cc    2 p electrons
      data (aimo_atstcom(2,i+1), i=1,4) /3,5,1,1/

cc    3 p eletrons
      data (aimo_atstcom(3,i+1), i=1,3) /1,5,3/

cc    4 p electrons
      data (aimo_atstcom(4,i+1), i=1,3) /3,5,1/

cc    5 p electrons
      data  aimo_atstcom(5,2)           /3/

cc    6 p electrons
      data  aimo_atstcom(6,2)           /1/

cc    1 s electron
      data  aimo_atstcom(7,2)           /1/

cc    2 s electrons
      data  aimo_atstcom(8,2)           /1/

cc    0 s electrons
      data  aimo_atstcom(9,2)           /1/

c...
c...  1 p electron
      data (aimo_confdef(1,1,i), i=1,8) /7,1,-1,2,300,300,400,400/
      data (aimo_confdef(1,2,i), i=1,8) /7,1,-1,200,200,3,400,400/
      data (aimo_confdef(1,3,i), i=1,8) /7,1,-1,200,200,300,300,4/

cc    doub P
      data  aivb_statelab(1,1,1)           /'doub P         '/

      data (aimo_statedef(1,1,1,i), i=1,3) /1,1,1/
      data  aivb_statelab(1,1,2)           /'|L=1 M=1, cos> '/
      data   aivb_conflab(1,1,1)  /'s2px '/

      data (aimo_statedef(1,1,2,i), i=1,3) /1,1,2/
      data  aivb_statelab(1,1,3)           /'|L=1 M=1, sin> '/
      data   aivb_conflab(1,1,2)  /'s2py '/

      data (aimo_statedef(1,1,3,i), i=1,3) /1,1,3/
      data  aivb_statelab(1,1,4)           /'|L=1 M=0, cos> '/
      data   aivb_conflab(1,1,3)  /'s2pz '/

c...
c...  2 p electrons
      data (aimo_confdef(2,1,i), i=1,7) /6,1,-1,2,300,300,4/
      data (aimo_confdef(2,2,i), i=1,7) /6,1,-1,200,200,3,4/
      data (aimo_confdef(2,3,i), i=1,7) /6,1,-1,2,3,400,400/
      data (aimo_confdef(2,4,i), i=1,9) /8,1,-1,2,-2,300,300,400,400/
      data (aimo_confdef(2,5,i), i=1,9) /8,1,-1,200,200,3,-3,400,400/
      data (aimo_confdef(2,6,i), i=1,7) /6,1,-1,2,-3,400,400/
      data (aimo_confdef(2,7,i), i=1,7) /6,1,-1,-2,3,400,400/
      data (aimo_confdef(2,8,i), i=1,7) /6,1,-1,2,300,300,-4/
      data (aimo_confdef(2,9,i), i=1,7) /6,1,-1,-2,300,300,4/
      data (aimo_confdef(2,10,i), i=1,7) /6,1,-1,200,200,3,-4/
      data (aimo_confdef(2,11,i), i=1,7) /6,1,-1,200,200,-3,4/
      data (aimo_confdef(2,12,i), i=1,9) /8,1,-1,200,200,300,300,4,-4/
      data (aimo_confdef(2,13,i), i=1,5) /4,1,2,3,4/

cc    trip P
      data  aivb_statelab(2,1,1)           /'trip P         '/

      data (aimo_statedef(2,1,1,i), i=1,3) /1,1,1/
      data  aivb_statelab(2,1,2)           /'|L=1 M=1, cos> '/
      data   aivb_conflab(2,1,1)  /'s2pxpz '/

      data (aimo_statedef(2,1,2,i), i=1,3) /1,1,2/
      data  aivb_statelab(2,1,3)           /'|L=1 M=1, sin> '/
      data   aivb_conflab(2,1,2)  /'s2pypz '/

      data (aimo_statedef(2,1,3,i), i=1,3) /1,1,3/
      data  aivb_statelab(2,1,4)           /'|L=1 M=0, cos> '/
      data   aivb_conflab(2,1,3)  /'s2pxpy '/

cc    sing D
      data  aivb_statelab(2,2,1)           /'sing D         '/

      data (aimo_statedef(2,2,1,i), i=1,5) /2,1,4,-1,5/
      data  aivb_statelab(2,2,2)           /'|L=2 M=2, cos> '/
      data   aivb_conflab(2,2,1)  /'s2(px2-py2) '/

      data (aimo_statedef(2,2,2,i), i=1,5) /2,1,6,-1,7/
      data  aivb_statelab(2,2,3)           /'|L=2 M=2, sin> '/
      data   aivb_conflab(2,2,2)  /'s2(pxpy|-px|py) '/

      data (aimo_statedef(2,2,3,i), i=1,5) /2,1,8,-1,9/
      data  aivb_statelab(2,2,4)           /'|L=2 M=1, cos> '/
      data   aivb_conflab(2,2,3)  /'s2(pxpz|-px|pz) '/

      data (aimo_statedef(2,2,4,i), i=1,5) /2,1,10,-1,11/
      data  aivb_statelab(2,2,5)           /'|L=2 M=1, sin> '/
      data   aivb_conflab(2,2,4)  /'s2(pypz|-py|pz) '/

      data (aimo_statedef(2,2,5,i), i=1,7) /3,1,4,1,5,-2,12/
      data  aivb_statelab(2,2,6)           /'|L=2 M=0, cos> '/
      data   aivb_conflab(2,2,5)  /'s2(px2+py2)-2(s2pz2) '/

cc    sing S
      data  aivb_statelab(2,3,1)           /'sing S         '/

      data (aimo_statedef(2,3,1,i), i=1,7) /3,1,4,1,5,1,12/
      data  aivb_statelab(2,3,2)           /'|L=0 M=0, cos> '/
      data   aivb_conflab(2,3,1)  /'s2pz2 '/

cc    quin S
      data  aivb_statelab(2,4,1)           /'quin S         '/

      data (aimo_statedef(2,4,1,i), i=1,3) /1,1,13/
      data  aivb_statelab(2,4,2)           /'|L=0 M=0, cos> '/
      data   aivb_conflab(2,4,1)  /'spxpypz '/

c...
c...  3 p electrons
      data (aimo_confdef(3,1,i), i=1,6) /5,1,-1,2,3,4/
      data (aimo_confdef(3,2,i), i=1,8) /7,1,-1,2,-2,300,300,4/
      data (aimo_confdef(3,3,i), i=1,8) /7,1,-1,200,200,3,-3,4/
      data (aimo_confdef(3,4,i), i=1,6) /5,1,-1,2,-3,4/
      data (aimo_confdef(3,5,i), i=1,6) /5,1,-1,-2,3,4/
      data (aimo_confdef(3,6,i), i=1,8) /7,1,-1,2,3,-3,400,400/
      data (aimo_confdef(3,7,i), i=1,8) /7,1,-1,2,300,300,4,-4/
      data (aimo_confdef(3,8,i), i=1,8) /7,1,-1,2,-2,3,400,400/
      data (aimo_confdef(3,9,i), i=1,8) /7,1,-1,200,200,3,4,-4/
      data (aimo_confdef(3,10,i), i=1,6) /5,1,-1,2,3,-4/

cc    quar S
      data  aivb_statelab(3,1,1)           /'quar S         '/

      data (aimo_statedef(3,1,1,i), i=1,3) /1,1,1/
      data  aivb_statelab(3,1,2)           /'|L=0 M=0, cos> '/
      data   aivb_conflab(3,1,1)  /'s2pxpypz '/

cc    doub D
      data  aivb_statelab(3,2,1)           /'doub D         '/

      data (aimo_statedef(3,2,1,i), i=1,5) /2,1,2,-1,3/
      data  aivb_statelab(3,2,2)           /'|L=2 M=2, cos> '/
      data   aivb_conflab(3,2,1)  /'s2(px2-py2)pz '/

      data (aimo_statedef(3,2,2,i), i=1,5) /2,1,4,-1,5/
      data  aivb_statelab(3,2,3)           /'|L=2 M=2, sin> '/
      data   aivb_conflab(3,2,2)  /'s2(pxpy|-px|py)pz '/

      data (aimo_statedef(3,2,3,i), i=1,5) /2,1,6,-1,7/
      data  aivb_statelab(3,2,4)           /'|L=2 M=1, cos> '/
      data   aivb_conflab(3,2,3)  /'s2px(py2-pz2) '/

      data (aimo_statedef(3,2,4,i), i=1,5) /2,1,8,-1,9/
      data  aivb_statelab(3,2,5)           /'|L=2 M=1, sin> '/
      data   aivb_conflab(3,2,4)  /'s2(px2py-pypz2) '/

      data (aimo_statedef(3,2,5,i), i=1,7) /3,1,4,1,5,-2,10/
      data  aivb_statelab(3,2,6)           /'|L=2 M=0, cos> '/
      data   aivb_conflab(3,2,5)  
     1       /'s2(pxpy|+px|py)pz-2(s2pxpypz|) '/

cc    doub P
      data  aivb_statelab(3,3,1)           /'doub P         '/

      data (aimo_statedef(3,3,1,i), i=1,5) /2,1,7,1,6/
      data  aivb_statelab(3,3,2)           /'|L=1 M=1, cos> '/
      data   aivb_conflab(3,3,1)  /'s2px(pz2+py2)'/
c      data   aivb_conflab(3,3,1)  /'s2 px pz2 '/

      data (aimo_statedef(3,3,2,i), i=1,5) /2,1,9,1,8/
      data  aivb_statelab(3,3,3)           /'|L=1 M=1, sin> '/
      data   aivb_conflab(3,3,2)  /'s2py(pz2+px2) '/
c      data   aivb_conflab(3,3,2)  /'s2 py pz2 '/

      data (aimo_statedef(3,3,3,i), i=1,5) /2,1,2,1,3/
c      data (aimo_statedef(3,3,3,i), i=1,11) /5,1,2,1,3,2,10,-1,4,-1,5/
      data  aivb_statelab(3,3,4)           /'|L=1 M=0, cos> '/
      data   aivb_conflab(3,3,3)
     1 /'s2(px2+py2)pz '/
cmarcin     1 /'s2(px2+py2)pz+[2(s2pxpypz|)-s2(pxpy|+px|py)pz] '/

c...
c...  4 p electrons
      data (aimo_confdef(4,1,i), i=1,7) /6,1,-1,2,3,-3,4/
      data (aimo_confdef(4,2,i), i=1,7) /6,1,-1,2,-2,3,4/
      data (aimo_confdef(4,3,i), i=1,7) /6,1,-1,2,3,4,-4/
      data (aimo_confdef(4,4,i), i=1,9) /8,1,-1,2,-2,300,300,4,-4/
      data (aimo_confdef(4,5,i), i=1,9) /8,1,-1,200,200,3,-3,4,-4/
      data (aimo_confdef(4,6,i), i=1,7) /6,1,-1,2,-3,4,-4/
      data (aimo_confdef(4,7,i), i=1,7) /6,1,-1,-2,3,4,-4/
      data (aimo_confdef(4,8,i), i=1,7) /6,1,-1,2,3,-3,-4/
      data (aimo_confdef(4,9,i), i=1,7) /6,1,-1,-2,3,-3,4/
      data (aimo_confdef(4,10,i), i=1,7) /6,1,-1,2,-2,3,-4/
      data (aimo_confdef(4,11,i), i=1,7) /6,1,-1,2,-2,-3,4/
      data (aimo_confdef(4,12,i), i=1,9) /8,1,-1,2,-2,3,-3,400,400/

cc    trip P
      data  aivb_statelab(4,1,1)           /'trip P         '/

      data (aimo_statedef(4,1,1,i), i=1,3) /1,1,1/
      data  aivb_statelab(4,1,2)           /'|L=1 M=1, cos> '/
      data   aivb_conflab(4,1,1)           /'s2pxpy2pz '/

      data (aimo_statedef(4,1,2,i), i=1,3) /1,1,2/
      data  aivb_statelab(4,1,3)           /'|L=1 M=1, sin> '/
      data   aivb_conflab(4,1,2)           /'s2px2pypz '/

      data (aimo_statedef(4,1,3,i), i=1,3) /1,1,3/
      data  aivb_statelab(4,1,4)           /'|L=1 M=0, cos> '/
      data   aivb_conflab(4,1,3)           /'s2pxpypz2 '/

cc    sing D
      data  aivb_statelab(4,2,1)           /'sing D         '/

      data (aimo_statedef(4,2,1,i), i=1,5) /2,1,4,-1,5/
      data  aivb_statelab(4,2,2)           /'|L=2 M=2, cos> '/
      data   aivb_conflab(4,2,1)  /'s2(px2-py2)pz2 '/

      data (aimo_statedef(4,2,2,i), i=1,5) /2,1,6,-1,7/
      data  aivb_statelab(4,2,3)           /'|L=2 M=2, sin> '/
      data   aivb_conflab(4,2,2)  /'s2(pxpy|-px|py)pz2 '/

      data (aimo_statedef(4,2,3,i), i=1,5) /2,1,8,-1,9/
      data  aivb_statelab(4,2,4)           /'|L=2 M=1, cos> '/
      data   aivb_conflab(4,2,3)  /'s2(pxpy2pz|-px|py2pz) '/

      data (aimo_statedef(4,2,4,i), i=1,5) /2,1,10,-1,11/
      data  aivb_statelab(4,2,5)           /'|L=2 M=1, sin> '/
      data   aivb_conflab(4,2,4)  /'s2px2(pypz|-py|pz) '/

      data (aimo_statedef(4,2,5,i), i=1,7) /3,1,4,1,5,-2,12/
      data  aivb_statelab(4,2,6)           /'|L=2 M=0, cos> '/
      data   aivb_conflab(4,2,5)  
     1       /'s2(px2+py2)pz2-2(s2px2py2) '/

cc    sing S
      data  aivb_statelab(4,3,1)           /'sing S         '/
      data (aimo_statedef(4,3,1,i), i=1,7) /3,1,12,1,4,1,5/
      data  aivb_statelab(4,3,2)           /'|L=0 M=0, cos> '/
      data   aivb_conflab(4,3,1)  /'s2px2py2 '/

c...
c...  5 p electrons
      data (aimo_confdef(5,1,i), i=1,8) /7,1,-1,2,3,-3,4,-4/
      data (aimo_confdef(5,2,i), i=1,8) /7,1,-1,2,-2,3,4,-4/
      data (aimo_confdef(5,3,i), i=1,8) /7,1,-1,2,-2,3,-3,4/

cc    doub P
      data  aivb_statelab(5,1,1)           /'doub P         '/

      data (aimo_statedef(5,1,1,i), i=1,3) /1,1,1/
      data  aivb_statelab(5,1,2)           /'|L=1 M=1, cos> '/
      data   aivb_conflab(5,1,1)  /'s2pxpy2pz2 '/

      data (aimo_statedef(5,1,2,i), i=1,3) /1,1,2/
      data  aivb_statelab(5,1,3)           /'|L=1 M=1, sin> '/
      data   aivb_conflab(5,1,2)  /'s2px2pypz2 '/

      data (aimo_statedef(5,1,3,i), i=1,3) /1,1,3/
      data  aivb_statelab(5,1,4)           /'|L=1 M=0, cos> '/
      data   aivb_conflab(5,1,3)  /'s2px2py2pz '/

c...
c...  6 p electrons
      data (aimo_confdef(6,1,i), i=1,9) /8,1,-1,2,-2,3,-3,4,-4/

cc    sing S
      data  aivb_statelab(6,1,1)           /'sing S         '/

      data (aimo_statedef(6,1,1,i), i=1,3) /1,1,1/
      data  aivb_statelab(6,1,2)           /'|L=0 M=0, cos> '/
      data   aivb_conflab(6,1,1)  /'s2px2py2pz2 '/

c...
c...  1 s electron
      data (aimo_confdef(7,1,i), i=1,2)  /1,1/

cc    doub S
      data  aivb_statelab(7,1,1)         /'doub S         '/

      data (aimo_statedef(7,1,1,i), i=1,3) /1,1,1/
      data  aivb_statelab(7,1,2)         /'|L=0 M=0, cos> '/
      data   aivb_conflab(7,1,1) /'s '/

c...
c...  2 s electrons
      data (aimo_confdef(8,1,i), i=1,3) /2,1,-1/

cc    sing S
      data  aivb_statelab(8,1,1)           /'sing S         '/

      data (aimo_statedef(8,1,1,i), i=1,3) /1,1,1/
      data  aivb_statelab(8,1,2)           /'|L=0 M=0, cos> '/
      data   aivb_conflab(8,1,1)  /'s2 '/

c...
c...  0 s electrons
      data (aimo_confdef(9,1,i), i=1,3) /2,100,100/

cc    no state
      data  aivb_statelab(9,1,1)           /'no el. '/

      data (aimo_statedef(9,1,1,i), i=1,3) /1,1,1/
      data  aivb_statelab(9,1,2)           /' '/
      data   aivb_conflab(9,1,1)  /'no el. '/

      end

      subroutine allocmem(name,size,where,file,status)
      implicit none
c
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/iofile)
INCLUDE(common/aivb)
c
      character(len=*) name,where,file
      integer :: size,status,idbg12
c
      idbg12 = aimo_debugpr(12)

      if ( status == 0 ) then
        if ( idbg12 .eq. 1 ) then
          write(iwr,11) 
     1    ' allocating ',size,' words for ',name,
     2    ' in subr. ',where,' at ',file
11        format(' >> ',5('==='),a,I8,6(a))
        end if
      else 
        write(iwr,12) 
     1  ' bad allocation of ',size,' words for ',name,
     2  ' in subr. ',where,' at ',file
12      format(' >> ',5('==='),a,I8,6(a))
        call caserr(' memory alloc error ')
      end if

      end

      subroutine deallocmem(name,where,file,status)
      implicit none
c
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/iofile)
INCLUDE(common/aivb)
c
      character(len=*) name,where,file
      integer :: status,idbg12
c
      idbg12 = aimo_debugpr(12)

      if ( status == 0 ) then
        if ( idbg12 .eq. 1 ) then
          write(iwr,11) 
     1    ' deallocating ',name,
     2    ' in subr. ',where,' at ',file
11        format(' << ',5('=-='),6(a))
        end if
      else
        write(iwr,12) 
     1  ' bad deallocation of ',name,
     2  ' in subr. ',where,' at ',file
12      format(' << ',5('=-='),6(a))
        call caserr(' memory alloc error ')
      end if

      end
