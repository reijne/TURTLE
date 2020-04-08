      subroutine JDriver(ao_tag,adens,bdens,kma,kmb)
c     
c     Anything in capitals is taken from crystal
c     
      implicit none
INCLUDE(common/claddd_crys.hf77)
INCLUDE(common/dfacomm_crys.hf77)
INCLUDE(common/tcomm_crys.hf77)
INCLUDE(common/basato.hf77)
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_basis)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
c
INCLUDE(../m4/common/timeperiods)
c
c scratch space, for the moment hardwired in but is a good
c candidate for dynamic memory
c
      integer shells_max,shellp_max
      parameter(shells_max=600)
      parameter(shellp_max=(shells_max*(shells_max+1))/2) 
      REAL dvec_scr(100)
c     integer ishlp_num,ishlp_list(shellp_max,2)
c     integer qshlp_num,qshlp_list(shellp_max),qatm_list(200)
      integer INZILA(shells_max),mono_shls
      integer bielec_list(shells_max),bi_shls
      REAL s(1225)
      REAL CJ(1225)
      REAL Sd(shells_max,49)
c
      integer ao_tag
      REAL adens(*),bdens(*)
      REAL kma(*),kmb(*)

      integer lshp,lxyz,lmon,lbie,lshl,lbf,lpol,i,lbf2
      integer p1,p2,MOCFN,lcl1,lcl2
      integer q1,q2
      integer bfp_num,bfq_num,pole_num
      integer ipos
      integer mono_total,biel_total,scre_total
      integer IDIPO,itol1,itol2
      
      REAL X2(3),X1(3)
      REAL X22OLD,Y22OLD,Z22OLD
      REAL enel,fact
      REAL ovtest,test

      call start_time_period(TP_DFT_JMULT)
      call start_time_period(TP_DFT_JMULT_INTRO)

      mono_total=0
      biel_total=0
      scre_total=0
C     
C     Initialise various common blocks used by cryspak. Also calculate
C     various basis set quantaties used by other routines
      if(totshl(ao_tag).gt.shells_max) then
        write(6,*) 'Increase size of shells_max'
        stop
      endif
      IDIPO = poleexp_num
      itol1 = over_tol
      itol2 = pener_tol
      ACCFAC=0.1d0**itol1
      test  = log(0.1d0)*dble(itol1)
      TTTTT2= log(0.1d0)*dble(itol2)-1.5d0*log(2.d0)
      write(6,*) 'Gaussian Overlap Screen Factor:',test
      write(6,*) 'Coulomb Penetration Factor:    ',TTTTT2
      pole_num=(IDIPO+1)**2
      call congen
      call buxtab(13,0.016d0,16.d0)
      call basato_fill(ao_tag)

c     call Gauss_Rys_init
C     
C     Find all interacting shell pairs
c     call basis_intshlp(ao_tag,ishlp_num,ishlp_list)
      do lcl1=1,totshl(ao_tag)
        do lcl2=1,49
          Sd(lcl1,lcl2)=0.0d0
        enddo
      enddo

      call end_time_period(TP_DFT_JMULT_INTRO)
      call start_time_period(TP_DFT_JMULT_SB)

C     Form multipole moments for interacting shell pairs
      do q1=1,totshl(ao_tag)
         X1OLD  =XL(1,q1)
         Y1OLD  =XL(2,q1)
         Z1OLD  =XL(3,q1) 
         do q2=1,totshl(ao_tag)
            do lxyz=1,3
               X1(lxyz)=XL(lxyz,q1)
               X2(lxyz)=XL(lxyz,q2)-X1(lxyz)
            enddo
            BA(1)=X2(1)
            BA(2)=X2(2)
            BA(3)=X2(3)
            RSQUAR=BA(1)**2+BA(2)**2+BA(3)**2
            if(ovtest(q1,q2,rsquar).ge.test) then
              X22OLD =XL(1,q2)-X1OLD
              Y22OLD =XL(2,q2)-Y1OLD
              Z22OLD =XL(3,q2)-Z1OLD
              EX2OLD=EXAD(q2)
              EX1OLD=EXAD(q1)
              DDDDDD=EX1OLD+EX2OLD
              FA2OLD=EX2OLD/DDDDDD
              DDDDDD=DDDDDD*0.5d0
              FATOLD=FA2OLD*EX1OLD*0.5d0
              X21OLD=X22OLD
              Y21OLD=Y22OLD
              Z21OLD=Z22OLD
              call POLIPA(q1,q2,IDIPO,S,0.0d0,0.0d0,0.0d0)

              bfq_num=LATAO(q1)*LATAO(q2)
c     

              call Pls_fill(ao_tag,q1,q2,adens,dvec_scr)
c     
c     pls_fill assumes that second index corresponds to fastest
c     varying index
c     
c           write(6,*)'patch of density',q1,q2
c           call matprt(dvec_scr,LATAO(q1),LATAO(q2),latao(q2),1,6)

            ipos=0
            do lpol=1,pole_num
               do lbf=1,bfq_num
                  ipos=ipos+1
                  Sd(q1,lpol)=Sd(q1,lpol)+S(ipos)*dvec_scr(lbf)
ccc   write(6,*) 'S:',lpol,lbf,S(ipos),dvec_scr(lbf)
               enddo
            enddo
           endif
         enddo                  ! end q2
      enddo                     ! end q1

c      write(6,*)'SB'
c      call matprt(sd,pole_num,totshl(ao_tag),30,1,6)

      call end_time_period(TP_DFT_JMULT_SB)

C     
C     Charge normalisation
c
      enel=0.d0
      do lshl=1,totshl(ao_tag)
         enel=enel+Sd(lshl,1) 
      enddo
      fact=dble(nelectrons)/enel
      write(6,*) 'fact:',fact,nelectrons,enel
      do lshl=1,totshl(ao_tag)
         do lpol=1,pole_num
            Sd(lshl,lpol)=Sd(lshl,lpol)*fact
         enddo
      enddo
C     

      call start_time_period(TP_DFT_JMULT_FOCK)

C     Loop over all interacting shell pairs
      do p1=1,totshl(ao_tag)
         do p2=1,p1
            do lxyz=1,3
               X1(lxyz)=XL(lxyz,p1)
               X2(lxyz)=XL(lxyz,p2)-X1(lxyz)
            enddo
            BA(1)=X2(1)
            BA(2)=X2(2)
            BA(3)=X2(3)
            RSQUAR=BA(1)**2+BA(2)**2+BA(3)**2
            if(ovtest(p1,p2,rsquar).ge.test) then
              bfp_num=LATAO(p1)*LATAO(p2)
              X1OLD  =XL(1,p1)
              Y1OLD  =XL(2,p1)
              Z1OLD  =XL(3,p1)
              X22OLD =XL(1,p2)-X1OLD
              Y22OLD =XL(2,p2)-Y1OLD
              Z22OLD =XL(3,p2)-Z1OLD
              EX2OLD=EXAD(p2)
              EX1OLD=EXAD(p1)
              DDDDDD=EX1OLD+EX2OLD
              FA2OLD=EX2OLD/DDDDDD
              DDDDDD=DDDDDD*0.5d0
              FATOLD=FA2OLD*EX1OLD*0.5d0
              X21OLD=X22OLD
              Y21OLD=Y22OLD
              Z21OLD=Z22OLD
C     
C     Classify shell pair with other shells to see if interaction is in
C     bielectronic zone

              call CLASSS(INZILA,mono_shls,bielec_list,bi_shls)
              mono_total=mono_total+mono_shls
              biel_total=biel_total+bi_shls
C
C     Form bielectronic integrals for the current shell pair
c           write(6,*) 'Going to Jbielec'
c           call Jbielec(ao_tag,p1,p2,bielec_list,bi_shls,
c    &                   adens,kma)
C     
C     Form monoelectronic integrals for the current shell pair
              do lcl1=1,bfp_num
                dvec_scr(lcl1)=0.d0
              enddo
              do lmon=1,mono_shls
C     Find interacting shells with shell in mono zone
                q1=INZILA(lmon)
                X1OLD  =XL(1,q1)
                Y1OLD  =XL(2,q1)
                Z1OLD  =XL(3,q1)
                do lcl1=1,1225
                  CJ(lcl1)=0.0d0
                enddo

c     write(6,*) 'Shells:',p1,p2,q1

                call start_time_period(TP_DFT_JMULT_CJAT0)
                call cjat0(p1,p2,IDIPO,bfp_num,q1,
     &               X1OLD,Y1OLD,Z1OLD,X1,CJ)
                call end_time_period(TP_DFT_JMULT_CJAT0)
               
c              write(6,*)'CJ for q1=',q1
c              call matprt(cj,pole_num,bfp_num,bfp_num,1,6)

                i=0
 
                do lpol=1,pole_num
                  do lbf=1,bfp_num
                     i=i+1
                     dvec_scr(lbf)=dvec_scr(lbf)+CJ(i)*Sd(q1,lpol)
                  enddo
               enddo
            enddo

c           write(6,*)'patch of fock',p1,p2
c           call matprt(dvec_scr,LATAO(p1),LATAO(p2),latao(p2),1,6)


            call Fmn_fill(ao_tag,p1,p2,dvec_scr,kma,kmb)
            else
              scre_total=scre_total+1
            endif
         enddo                  ! end p2

      enddo                     ! end p1
c
c
      write(6,*) '************************************'
      write(6,*) 'Number of Monopole Interactions:    ',mono_total
      write(6,*) 'Number of Bielectronic Interaction: ',biel_total
      write(6,*) 'Number of screens on p1-p2:         ',scre_total
      write(6,*) 'Charge Normalisation Factor:        ',fact
c     write(6,*) 'FINISHED IN JDRIVER'
c     stop


      call end_time_period(TP_DFT_JMULT_FOCK)
      call end_time_period(TP_DFT_JMULT)

      return
      end
      subroutine Fmn_fill(tag,ii,jj,dvec_scr,kma,kmb)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/nshel)
      integer tag,ii,jj
      REAL dvec_scr(*)
      REAL kma(*),kmb(*),fac
      integer ibas_num,jbas_num
      integer lbai,lbaj,inum,jnum
      integer dpos,ipos,jpos,maxa,mina,ijbas_num
      dpos=0
      ibas_num=0
      jbas_num=0
      ipos=kloc(ii)
      jpos=kloc(jj)
      inum=kmax(ii)-kmin(ii)
      jnum=kmax(jj)-kmin(jj)
      do lbai=0,inum
         ibas_num=ipos+lbai
         do lbaj=0,jnum
            dpos=dpos+1
            jbas_num=jpos+lbaj
            maxa=max(ibas_num,jbas_num)
            mina=min(ibas_num,jbas_num)
            ijbas_num=(((maxa*maxa)-maxa)/2)+mina
            kma(ijbas_num)=dvec_scr(dpos)
c     write(6,*) 'kma:',ijbas_num,kma(ijbas_num)
         enddo
      enddo
      return
      end 
      subroutine Pls_fill(tag,ii,jj,adens,dvec_scr)
      REAL adens(*)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/nshel)
      integer tag,ii,jj
      integer ibas_num,jbas_num
      REAL dvec_scr(*),fac
      integer lbai,lbaj,ijbas_num,inum,jnum
      integer dpos,ipos,jpos,maxa,mina
      dpos=0
      ibas_num=0
      jbas_num=0
      ipos=kloc(ii)
      jpos=kloc(jj)
      inum=(kmax(ii)-kmin(ii))
      jnum=(kmax(jj)-kmin(jj))


      do lbai=0,inum
         ibas_num=ipos+(lbai)

         do lbaj=0,jnum
            jbas_num=jpos+(lbaj)

            dpos=dpos+1
c     fac=2.d0
c     if(ibas_num.eq.jbas_num) fac=1.d0
            maxa=max(ibas_num,jbas_num)
            mina=min(ibas_num,jbas_num)
            ijbas_num=(((maxa*maxa)-maxa)/2)+mina
c     dvec_scr(dpos)=adens(ijbas_num)*fac
            dvec_scr(dpos)=adens(ijbas_num)
         enddo
      enddo
      return
      end

      subroutine matprt(a,ncol,nrow,mcol,mrow,iwr)
      implicit none
      REAL a
      integer ncol,nrow,mcol,mrow,iwr
      integer istart, ifin, irow, icol

      dimension a(*)
      integer width
      parameter(width=10)

      istart = 1

      do while(istart .le. ncol)

         ifin=min(istart+width-1,ncol)

         write(iwr,100)(icol,icol=istart,ifin)
         do irow=1,nrow
            write(iwr,101)irow,(a(1+(irow-1)*mrow + (icol-1)*mcol),
     &           icol=istart,ifin)
         enddo
         istart = ifin + 1
         write(iwr,'(/)')
      enddo

 100  format(5x,10(2x,i3,3x))
 101  format(i5,10f8.4)

      end



      subroutine basis_intshlp(tag,ishlp_num,ishlp_list)
C
C Routine to calculate set of interacting shell pairs. At the moment
C no screening is used.
C
      implicit none
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_basis)
INCLUDE(common/dft_mol_info)
      integer ishlp_num
      integer ishlp_list(200,2)
      integer tag 
  
      integer latm,lshl,jatm,jshl
      integer ltype,jtype,ishls,jshls
      ishls=1
      jshls=1
      ishlp_num=0
      do latm=1,natoms
        ltype=atom_tag(tag,latm)
        do lshl=1,num_shl(tag,ltype)
          jshls=1
c         do jatm=1,latm
          do jatm=1,natoms
            jtype=atom_tag(tag,jatm)
            do jshl=1,num_shl(tag,jtype)
c           do jshl=1,lshl
c
              ishlp_num=ishlp_num+1
              ishlp_list(ishlp_num,1)=ishls
              ishlp_list(ishlp_num,2)=jshls
              jshls=jshls+1
c
            enddo
          enddo
          ishls=ishls+1
        enddo
      enddo   
c     write(6,*) 'NUMBER OF SHELL PAIRS:',ishlp_num
c     do lshl=1,ishlp_num
c       write(6,*) 'PAIR:',lshl,ishlp_list(lshl,1),ishlp_list(lshl,2)
c     enddo
      return
      end 

      subroutine basis_intshlq(tag,shl,interact_num,interact_list,
     & atom_list)
 
      implicit none

INCLUDE(common/dft_parameters)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_basis)
      integer tag,shl
      integer interact_num,interact_list(*),atom_list(*)

      integer latm,lshl
      integer atom_type,shl_count
      interact_num=0
      shl_count=0
      do latm=1,natoms
        atom_type=atom_tag(tag,latm)
        do lshl=1,num_shl(tag,atom_type) 
          shl_count=shl_count+1
c         if(shl_count.le.shl) then
            interact_num=interact_num+1
            interact_list(interact_num)=shl_count
            atom_list(interact_num)=latm
c         endif
        enddo  
      enddo
      return
      end
c
C Fills basato 
c
      subroutine basato_fill(tag)
      implicit none
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_basis)
INCLUDE(common/dft_mol_info)
      integer maxnor
      parameter(maxnor=2000)
      REAL cunorm_list
      common/phillbas/cunorm_list(maxnor,2)

INCLUDE(common/basato.hf77)
INCLUDE(common/parinf.hf77)
INCLUDE(common/servi_inpbas.hf77)
      integer tag

      integer latm,lshl,lprm,i
      integer atom_type,nshl,crys,ccpd
      integer l,lh,ngauss,blook
      integer bf_crystal(4),bf_charge(4) 
      integer nelec,nelec_shl
    
      logical normal_sw
      integer nprim_tmp

      data bf_crystal/1,4,3,5/
      data bf_charge/2,8,6,10/

      normal_sw=.true.
      nshl = 0
      NSHPRI(1) = 1
      NDQ(1) = 0
      LAA(1) = 1
      crys = 1
c
      IOUT=6
      INF(3)=4
      INF(20)=totshl(tag)
      INF(24)=natoms
      if(totshl(tag).gt.LIM015) then
        write(6,*) 'Increase size of LIM015'
        stop
      endif
      if(natoms.gt.LIM016) then
        write(6,*) 'Increase size of LIM016'
        stop
      endif
c
      do latm=1,natoms
        nelec=ian(latm)
        atom_type=atom_tag(tag,latm)
        XA(1,latm)=atom_c(latm,1)
        XA(2,latm)=atom_c(latm,2)
        XA(3,latm)=atom_c(latm,3)
        AZNUC(latm)=dble(ian(latm))
        NAT(latm)=ian(latm)
        NSHPRI(latm+1)=NSHPRI(latm)+num_shl(tag,atom_type)
        IPSEUD(latm)=0 
        nprim_tmp=1
        do lshl=1,num_shl(tag,atom_type)
          nshl=nshl+1
c
          XL(1,nshl)=atom_c(latm,1)
          XL(2,nshl)=atom_c(latm,2)
          XL(3,nshl)=atom_c(latm,3)
c
          l      = angmom(tag,atom_type,lshl)
          lh     = hybrid(tag,atom_type,lshl)
          ngauss = nprim(tag,atom_type,lshl)
          blook  = pstart(tag,atom_type,lshl)
c
          if(l.eq.lh) then
            if(l.eq.1) LAT(nshl)=0
            if(l.eq.2) LAT(nshl)=l
            if(l.eq.3) LAT(nshl)=l
          else
            LAT(nshl)=1
          endif
c
          nelec_shl=bf_charge(LAT(nshl)+1)
          if(nelec.gt.nelec_shl) then
            CHE(nshl)=dble(nelec_shl)
            nelec=nelec-nelec_shl
          else
            CHE(nshl)=dble(nelec)
            nelec=0
          endif
c
          LAN(nshl)=ngauss
          LATAO(nshl)=bf_crystal(LAT(nshl)+1)
          LATOAT(nshl)=latm
          SCALE(nshl)=1.0d0
          ccpd=blook 
          do lprm=1,ngauss
            C1(crys)  = 0.0d0
            C2(crys)  = 0.0d0
            C3(crys)  = 0.0d0
            EXX(crys) = alpha(tag,atom_type,ccpd)
            if(LAT(nshl).eq.0) then
              C1(crys)=cont_coeff(tag,atom_type,ccpd,1)
cc
              if(normal_sw) then
                nprim_tmp=nprim_tmp+1
                C1(crys)=cunorm_list(crys,1)
c               write(6,*) 'C1 un-normalised:',crys,C1(crys)
              endif
cc
            endif
            if(LAT(nshl).eq.1) then
              C1(crys)=cont_coeff(tag,atom_type,ccpd,1)
              C2(crys)=cont_coeff(tag,atom_type,ccpd,2)
c
c
              if(normal_sw) then
                nprim_tmp=nprim_tmp+1
                C1(crys)=cunorm_list(crys,1)
                C2(crys)=cunorm_list(crys,2)
c               write(6,*) 'C1 un-normalised:',crys,C1(crys),C2(crys)
              endif
c
c
            endif
            if(LAT(nshl).eq.2) C2(crys)=cont_coeff(tag,atom_type,ccpd,2)
            if(LAT(nshl).eq.3) C3(crys)=cont_coeff(tag,atom_type,ccpd,3)
            ccpd=ccpd+1
            crys=crys+1
          enddo
          NDQ(nshl+1)=NDQ(nshl)+LATAO(nshl)
          LAA(nshl+1)=crys
        enddo
      enddo
      call GAUNOV
C
C debug print
c     call BASPRT
c     WRITE(IOUT,102)(I,SCALE(I),EXAD(I),CHE(I),I=1,nshl)
c     write(6,*) 'Leaving basato_fill'

102   FORMAT(1X,79('*')/3('   LA SCALE   EXAD    CHE')/
     *((3(I5,F6.2,1PE10.2,0PF4.0))))
      return
      end

      REAL function ovtest(shl1,shl2,rsquar)
      implicit none
INCLUDE(common/basato.hf77)
      integer shl1,shl2
      REAL rsquar
      REAL exad1,exad2
      REAL abi,p

      exad1=exad(shl1)
      exad2=exad(shl2)
     
      abi    = 1.0d0/(exad1+exad2)
      p      = exad1*exad2*abi
      ovtest = log(p*abi*4.0d0)*0.75d0
      ovtest = ovtest-rsquar*p
      return
      end
      subroutine ver_dft_mpole2(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/mpole2.m,v $
     +     "/
      data revision /
     +     "$Revision: 5774 $"
     +      /
      data date /
     +     "$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
