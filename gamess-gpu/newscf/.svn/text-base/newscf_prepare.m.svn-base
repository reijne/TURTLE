      Subroutine newscf_prepare( core, uhf, direct, enuclear, ehf, etot, 
     +     sz, s2, ek, vir )

      ! Allocate arrays that need GAMESS-UK memory alloc scheme

      Use newscf_numbers
      Use newscf_modules
      Use f90_io_data
      
      Implicit None
      
      Real( wp ), Dimension( * ), Intent( InOut ) :: core
      Logical                   , Intent( In    ) :: uhf
      Logical                   , Intent( In    ) :: direct
      Real( wp )                , Intent(   Out ) :: enuclear
      Real( wp )                , Intent( InOut ) :: ehf
      Real( wp )                , Intent(   Out ) :: etot
      Real( wp )                , Intent(   Out ) :: sz
      Real( wp )                , Intent(   Out ) :: s2
      Real( wp )                , Intent(   Out ) :: ek
      Real( wp )                , Intent(   Out ) :: vir

      Integer, External :: igmem_alloc_inf

      Integer :: l1, l2, l22, nshtri
      Integer :: ifock,  idmat
      Integer :: ifockb, idmatb
      Integer :: iprefa, irdmat
      Integer :: iblk25
      Character( Len = 16) :: fnm
      Character( Len = 14) :: snm

      fnm = "newscf_prepare.m"
      snm = "newscf_prepare"

      l1 = num
      l2 = l1 *( l1 + 1 ) / 2

      l22 = 2 * l2

      nshtri = ikyp( nshell )

      If( direct ) Then
         irdmat = igmem_alloc_inf( 2 * nshtri, fnm, snm, 'irdmat', 
     +                             IGMEM_NORMAL )
         iprefa = irdmat + nshtri
      End If

      If( uhf ) Then
         idmat  = igmem_alloc_inf( l22, fnm, snm, 'idmat', 
     +                             IGMEM_NORMAL )
         idmatb = idmat + l2
         ifock  = igmem_alloc_inf( l22, fnm, snm, 'ifock', 
     +                             IGMEM_NORMAL )
         ifockb = ifock + l2
      Else
         idmat  = igmem_alloc_inf( l2, fnm, snm, 'idmat', 
     +                             IGMEM_NORMAL )
         idmatb = idmat
         ifock  = igmem_alloc_inf( l2, fnm, snm, 'ifock', 
     +                             IGMEM_NORMAL )
         ifockb = ifock
      End If

      Call newscf_prepare_2( core, uhf, direct, enuclear, ehf, etot, sz, 
     +     s2, ek, vir, nshtri, Q( irdmat ), Q( iprefa ), 
     +     l2, Q( ifock ), Q( ifockb ), 
     +     Q( idmat ), Q( idmatb ) )

      Call gmem_free_inf( ifock, fnm, snm, 'ifock' )
      Call gmem_free_inf( idmat, fnm, snm, 'idmat' )
      If( direct ) Then
         Call gmem_free_inf( irdmat, fnm, snm, 'irdmat' )
      End If

      End Subroutine newscf_prepare

      Subroutine newscf_prepare_2( core, uhf, direct, enuclear, ehf, 
     +     etot, sz, s2, ek, vir,
     +     nrd, rdmat_work, prefac_work,
     +     nfk, focka_work, fockb_work, densa_work, densb_work )

      Use newscf_numbers
      Use newscf_modules
      Use newscf_routines
      
      Implicit None
      
      Real( wp ), Dimension( *     ), Intent( InOut ) :: core
      Logical                       , Intent( In    ) :: uhf
      Logical                       , Intent( In    ) :: direct
      Real( wp )                    , Intent(   Out ) :: enuclear
      Real( wp )                    , Intent( InOut ) :: ehf
      Real( wp )                    , Intent(   Out ) :: etot
      Real( wp )                    , Intent(   Out ) :: sz
      Real( wp )                    , Intent(   Out ) :: s2
      Real( wp )                    , Intent(   Out ) :: ek
      Real( wp )                    , Intent(   Out ) :: vir
      Integer                       , Intent( In    ) :: nrd
      Real( wp ), Dimension( 1:nrd ), Intent(   Out ) :: rdmat_work
      Real( wp ), Dimension( 1:nrd ), Intent(   Out ) :: prefac_work
      Integer                       , Intent( In    ) :: nfk
      Real( wp ), Dimension( 1:nfk ), Intent(   Out ) :: focka_work
      Real( wp ), Dimension( 1:nfk ), Intent(   Out ) :: fockb_work
      Real( wp ), Dimension( 1:nfk ), Intent(   Out ) :: densa_work
      Real( wp ), Dimension( 1:nfk ), Intent(   Out ) :: densb_work

      Real( wp ) :: accdi1

      character *8 title,guess
      common/restrz/title(12),guess

      if (guess.eq.'atdens'.or.guess.eq.'atoms' .or.
     +    (zguess.eq.'mosaved'. and. (guess.eq.'atoms'.or.
     +      guess.eq.'atdens') ) ) then
c
c         do one cycle scf with the density matrix from denat
c
         call initial_evecs( uhf, uhfatom, direct, rdmat_work, 
     +        prefac_work, 
     +        focka_work, fockb_work, densa_work, densb_work, core )

c$$$         if (irest.ne.0) go to 160
c
c        set zguess to 'anything' so denscf will be called but once
c
         zguess = 'anything'
         guess  = 'anything'
      end if
      ! HACK HACK HACK HACK
      !--------------------
      Call get_accdi1( accdi1 )

      Call newscf( core, uhf, direct, enuclear, ehf, etot, 
     +     sz, s2, ek, vir,
     +     rdmat_work, prefac_work,
     +     focka_work, fockb_work, densa_work, densb_work,
     +     accdi1 )

      End Subroutine newscf_prepare_2

      Subroutine dlc_error( s1, s2 )

      Implicit None

INCLUDE(../m4/common/errcodes)

      Character( Len = * ), Intent( In ) :: s1
      Character( Len = * ), Intent( In ) :: s2

      Call gamerr( s1, ERR_NO_CODE, ERR_INTERNAL, 
     +     ERR_ASYNC, ERR_NO_SYS )

      End Subroutine dlc_error
 
      subroutine print_conv(uhf)

      implicit none

      logical uhf

INCLUDE(../m4/common/newscf)
INCLUDE(../m4/common/iofile)

      logical opg_root

      integer i,j,linfo

      write(iwr,*)
      write(iwr,*)'Convergence procedure'
      write(iwr,*)'*********************'

      do i =1,nphase
         write(iwr,*)
         write(iwr,*)'Phase ',i 
         if(uhf)then
            write(iwr,*)'  Level ',
     &           shift(i,1),shift(i,2)
         else
            write(iwr,*)'  Level ',shift(i,1)
         endif
         if(diis(i))write(iwr,*)'  DIIS'
         if(lock_vec(i))write(iwr,*)'  Lock'
         if(restore_vec(i))write(iwr,*)'  Restore'
         if(new_diis(i))write(iwr,*)'  NewDIIS'
         if(extrap(i))write(iwr,*)'  Extrapolate',
     &        extrap_tol(i),extrap_coef(i)
         if(smear(i).ne.SMEAR_OFF) then
            if(smear(i).eq.SMEAR_ENERGY) then
               write(iwr,*) 'Smear: ',
     +              '  Initial Energy = ',
     +              new_esmear_start(i),
     +              '  Final Energy = ',
     +              new_esmear_final(i)
            else
               write(iwr,*) 'Smear: ',
     +              '  Energy is ', egap_scale( i ),
     +              ' * the HOMO/LUMO gap'
            endif
         endif
         do j = 1,nnext(i)

            if(nextphase(i,j).eq.-1)then
               write(iwr,*)'#  Abort calculation'
            else if(nextphase(i,j).eq.0)then
               write(iwr,*)'#  Converge calculation'
            else
               write(iwr,*)'#  Switch to phase ',nextphase(i,j)
               write(iwr,*)'   next ',nextphase(i,j)
            endif

            call strtrm( phase_info(i,j) , linfo )
            if ( linfo .gt. 0 ) then
               write(iwr,*)'# Info: ',phase_info(i,j)
            endif

            if(tester_chk(i,j) .eq. CONV_BELOW)then
               write(iwr,*)'     Tester below ',tester_val(i,j)
            else if (tester_chk(i,j) .eq. CONV_ABOVE)then
               write(iwr,*)'     Tester above ',tester_val(i,j)
            endif

            if(dele_chk(i,j) .eq. CONV_BELOW)then
               write(iwr,*)'#    Energy change'
               write(iwr,*)' dE below ',
     &              dele_val(i,j)
            else if (dele_chk(i,j) .eq. CONV_ABOVE)then
               write(iwr,*)'#     Energy change '
               write(iwr,*)' dE above ',
     &              dele_val(i,j)
            endif

            if(abs_dele_chk(i,j) .eq. CONV_BELOW)then
               write(iwr,*)'#     Absolute energy change'
               write(iwr,*)'    dEabs below',
     &              abs_dele_val(i,j)
            else if (abs_dele_chk(i,j) .eq. CONV_ABOVE)then
               write(iwr,*)'#     Absolute energy change'
               write(iwr,*)'    dEabs above',
     &              abs_dele_val(i,j)
            endif

            if( ncyc_chk(i,j) .eq. CONV_BELOW)then
               write(iwr,*)'#     Phase cycle count'
               write(iwr,*)'      ncyc below ',ncyc_val(i,j)
            else if ( ncyc_chk(i,j) .eq. CONV_ABOVE)then
               write(iwr,*)'#     Phase cycle count'
               write(iwr,*)'      ncyc above ',ncyc_val(i,j)
            endif

cjens
            if( ntotcyc_chk(i,j) .eq. CONV_BELOW)then
               write(iwr,*)'#     Total cycle count'
               write(iwr,*)'      ncyc below ',ntotcyc_val(i,j)
            else if (ntotcyc_chk(i,j) .eq. CONV_ABOVE)then
               write(iwr,*)'#     Total cycle count'
               write(iwr,*)'      ncyc above ',ntotcyc_val(i,j)
            endif


         enddo
      enddo

      end

c
c Check convergence
c
      integer function check_conv(uhf,iphase,tester,dele,niter_phase,
     &     niter,o1,o2,o3,o4,o5)

      implicit none

      logical uhf
      integer iphase
      REAL tester
      REAL dele
      integer niter_phase
      integer niter
      integer linfo

INCLUDE(../m4/common/newscf)      
INCLUDE(../m4/common/iofile)

      logical o1(maxphase), o2(maxphase), o3(maxphase), o4(maxphase),
     &        o5(maxphase)
      integer i
c
c  default is no change yet
c
      check_conv = iphase

      do i = 1, nnext(iphase)
c
c perform convergence/abort tests
c if more than one test is passed, the last one takes
c effect
c
         o1(i) = (tester_chk(iphase,i) .eq. CONV_INACTIVE)
         if(tester_chk(iphase,i) .eq. CONV_ABOVE .and. 
     &        tester .gt. tester_val(iphase,i))o1(i)=.true.
         if(tester_chk(iphase,i) .eq. CONV_BELOW .and. 
     &        tester .lt. tester_val(iphase,i))o1(i)=.true.

         o2(i) = (dele_chk(iphase,i) .eq. CONV_INACTIVE)
         if(dele_chk(iphase,i) .eq. CONV_ABOVE .and. 
     &        dele .gt. dele_val(iphase,i))o2(i)=.true.
         if(dele_chk(iphase,i) .eq. CONV_BELOW .and. 
     &        dele .lt. dele_val(iphase,i))o2(i)=.true.

         o3(i) = (abs_dele_chk(iphase,i) .eq. CONV_INACTIVE)
         if(abs_dele_chk(iphase,i) .eq. CONV_ABOVE .and. 
     &        abs(dele) .gt. abs_dele_val(iphase,i))o3(i)=.true.
         if(abs_dele_chk(iphase,i) .eq. CONV_BELOW .and. 
     &        abs(dele) .lt. abs_dele_val(iphase,i))o3(i)=.true.

         o4(i) = (ncyc_chk(iphase,i) .eq. CONV_INACTIVE)
         if(ncyc_chk(iphase,i) .eq. CONV_ABOVE .and. 
     &        niter_phase .gt. ncyc_val(iphase,i))o4(i)=.true.
         if(ncyc_chk(iphase,i) .eq. CONV_BELOW .and. 
     &        niter_phase .lt. ncyc_val(iphase,i))o4(i)=.true.


         o5(i) = (ntotcyc_chk(iphase,i) .eq. CONV_INACTIVE)
         if(ntotcyc_chk(iphase,i) .eq. CONV_ABOVE .and. 
     &        niter .gt. ntotcyc_val(iphase,i))o5(i)=.true.
         if(ntotcyc_chk(iphase,i) .eq. CONV_BELOW .and. 
     &        niter .lt. ntotcyc_val(iphase,i))o5(i)=.true.


         if( o1(i) .and. o2(i) .and. o3(i) .and. o4(i) .and. o5(i) )then
            check_conv = nextphase(iphase,i)

cjmht       see if we have owt to say about this phase change
            call strtrm( phase_info(iphase,i) , linfo )
            if ( linfo .gt. 1 ) then
               pchange_info = phase_info(iphase,i)
            else
               pchange_info = ' '
            endif

         endif

      enddo

      end

      subroutine default_conv(uhf)

      implicit none

      logical uhf

INCLUDE(../m4/common/newscf)
INCLUDE(../m4/common/iofile)

      integer iphase

      if(nphase .ne. 0)then
         write(iwr,*)'Using user-defined convergence scheme'
         return
      endif

      nphase = 2

      iphase=1
      lock_vec(iphase) = lock_vec(0)
      diis(iphase) = diis(0)
      shift(iphase,1) = shift(0,1)
      shift(iphase,2) = shift(0,2)
      nnext(iphase) = 1
      tester_chk(iphase,1) = CONV_BELOW
      tester_val(iphase,1) = 0.1d0
      nextphase(iphase,1) = 2
c
c  As phase 1, but turn off level shifter
c
      iphase=2
      lock_vec(iphase) = lock_vec(0)
      diis(iphase) = diis(0)
      shift(iphase,1) = 0.0d0
      shift(iphase,2) = 0.0d0
      nnext(iphase) = 1
      tester_chk(iphase,1) = CONV_BELOW
      tester_val(iphase,1) = tester_val(0,1)
      abs_dele_chk(iphase,1) = CONV_BELOW
      abs_dele_val(iphase,1) = abs_dele_val(0,1)
      nextphase(iphase,1) = 0

      end

      subroutine get_accdi1( value )

      Use newscf_numbers

INCLUDE(../m4/common/scfopt)

      Real( wp ), Intent( Out ) :: value

      value = accdi1

      End Subroutine get_accdi1

      subroutine scfsav_new(q,p,e,pop,ndaf,l1,l2,iblkp,iblke)

      Use distributed_matrices

      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/runlab)
INCLUDE(../m4/common/harmon)
c$$$      dimension q(*),p(*),e(*),pop(*)
      Type( matrix ) :: q
      Type( matrix ) :: p
      dimension e(*),pop(*)
c
      data m1/1/
      l0 = newbas0
      call putq_new(zcom,ztitle,e,pop,l1,l1,l0,
     *m1,m1,q,ndaf,iblk)
c$$$      call wrt3(p,l2,iblkp,idaf)
      Call matrix_write_triangle( p, iblkp, idaf )
      call wrt3(e,l1,iblke,idaf)
      return
      end


      subroutine putq_new(zcomm,ztit,eig,def,norb,norbn,ncolu,ieig,
     *ideff,q,mpos,iblkq)
c... standard e.vector outputting routine(+ header blocks)
c... allow VB/servec weird dimensioning
      Use distributed_matrices
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
c$$$      dimension q(*),zcomm(*),ztit(*),eig(*),def(*)
      Type( matrix ) :: q
      dimension zcomm(*),ztit(*),eig(*),def(*)
INCLUDE(../m4/common/sector)
      common/junkc/zcom(19),ztitle(10)
      common/blkin/value(maxorb),pop(maxorb),escf,
     *nbasis,newbas,ncol,ivalue,ipop
INCLUDE(../m4/common/machin)
INCLUDE(../m4/common/scfopt)
INCLUDE(../m4/common/tran)
      data m29/29/
       if(mpos.gt.508)call caserr(
     *'invalid section specified for vector output')
      len1=lensec(mach(8))
      nbsq=norb*norb
      m3 = 3
c
      if (zcomm(5).eq.'vb'.or.zcomm(5).eq.'servec') then
         nbsq = norb*ncolu
         m3 = 33
      end if
c
      lenv=lensec(nbsq)
      len2=lensec(mach(9))
      j=len2+lenv+len1+1
      call secput(mpos,m3,j,iblk)
      iblkq=iblk+len2+len1+1
      escf=etot
      do 100 i=1,10
 100  ztitle(i)=ztit(i)
      do 200 i=1,19
 200  zcom(i)=zcomm(i)
      call dcopy(ncolu,eig,1,value,1)
      call dcopy(ncolu,def,1,pop,1)
      nbasis=norb
      newbas=norbn
      ncol=ncolu
      ivalue=ieig
      ipop=ideff
      call wrtc(zcom,m29,iblk,numdu)
c...    clear non-existent vectors abd populations and make e.v. large
      if(ncol.ne.norb.and.zcomm(5).ne.'vb'.and.zcomm(5).ne.'servec')then
c$$$         call vclr(q(ncol*norb+1),1,(norb-ncol)*norb)
         do i=ncol+1,norb
            value(i) = 9999900.0d0 + i*0.1d0
            pop(i) = 0.0d0
         end do
      end if
c
      call wrt3s(value(1),mach(8),numdu)
      nav = lenwrd()
      call wrt3is(ilifc(1),mach(9)*nav,numdu)
c$$$      call wrt3(q(1),nbsq,iblkq,numdu)
      Call matrix_write_rectangle( q, iblkq, numdu )
      call clredx
      call revind
      return
      end
