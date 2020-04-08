      PROGRAM FORM31
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'arrays.inc'
      INCLUDE 'common.inc'
      INCLUDE 'init.inc'
*
      parameter (mxorbi=400)
      parameter (mxo=mxorbi+1)
*
* -- Files form MRDCI
      integer*4 idk,nform
      dimension idk(mxo)
      integer*4 oiir,lj2
      dimension oiir(MXORBI),ibias(9),lj2(8)
c     COMMON /I2/ OIJ(8),OXX(8),OISYM(8),OJSYM(36),OIIR(MXORB),                 
c    *OMJ(8),OKJ(8),ONJ(8),OLJ(8),OIBAL(8),                                     
c    *ONTIL(8),ONBAL(9),ODUM1(3),ONCOMP(100),                                   
c    *OITIL(8),OMCOMP(100)                                                      
* 
      NAMELIST/F31/ iprint,EZERO,nform
*
*----------------------------------------------------------------------*
*   Routine zur Erstellung eines STONEY files aus den MOLCAS files     *
*    TRAONE  ubd TRAINT.                                               *
*                                                                      *
*   Input files:                                                       *
*                                                                      *
*   - TRAONE                                                           *
*   - TRAINT                                                           *
*  (- INPORB)                                                          *)
*  (- TRANSF)                                                          *      
*                                                                      *     
*   Output files:                                                      *
*   - STONEY                                                           * 
*  (- TRANSF)                                                          *
*                                                                      *
*----------------------------------------------------------------------*
      Luorb  =12
      Lutra1 =3     
      Lutra2 =4     
      Lutra3 =14     
      Lusto  =31     
      IN=5
      iout=6
      trafil1='TRAONE'
      trafil2='TRAINT'
      trafil3='TRANSF'
      orbfile='INPORB'
      stofile='STONEY'
      iprint=1
      EZERO=0.0d0
      Trans=.false.
*     
*----------------------------------------------------------------------*
*     Konstruiere IDK                                                  *
*----------------------------------------------------------------------*
*
      idk(1) = 0
      Do i=2,mxo
         idk(i) = idk(i-1) + i - 1
      End Do 
*     
*----------------------------------------------------------------------*
*     Read NAMELIST input F31                                          *
*----------------------------------------------------------------------*
*     
      WRITE(IOUT,*) '======================================'
      WRITE(IOUT,*) '===  FORM31 =========================='
      WRITE(IOUT,*) '======================================'
      WRITE(IOUT,*) 
      Read(IN,F31)
      write(6,*)
      write(6,*) 'nform=',nform
      write(6,*) 'EZERO=',Ezero
*     
*----------------------------------------------------------------------*
*     Initialize arrays                                                *
*----------------------------------------------------------------------*
*     
      Do iSym=1,mxSym
        nFro(iSym)=0
        nDel(iSym)=0
        nOrb(iSym)=0
      End Do
      Blank(1:72)=' '
*     
*----------------------------------------------------------------------*
*     Read the input stream line by line and identify key command      *
*----------------------------------------------------------------------*
*     
100   Read(*,'(A)') Line
      If ( Line(1:72).eq.Blank(1:72) .or. Line(1:1).eq.'*' ) Goto 100
      Call UpCase(Line)
110   jCmd=0
      Do iCmd=1,nCmd
        If ( Line(1:lCmd).eq.CmdTab(iCmd)(1:lCmd) ) jCmd=iCmd
      End Do
      If ( jCmd.eq.0 ) Goto 992
*     
*----------------------------------------------------------------------*
*     Branch to the processing of the command sections                 *
*----------------------------------------------------------------------*
*     
      Goto (1,2,3),jCmd
*---  Process the "TITLe" command -------------------------------------*
1     mTit=0
15    Read(*,'(A)') Line
      If ( Line(1:72).eq.Blank(1:72) .or. Line(1:1).eq.'*' ) Goto 15
      Call UpCase(Line)
      Do iCmd=1,nCmd
        If ( Line(1:lCmd).eq.CmdTab(iCmd)(1:lCmd) ) Goto 110
      End Do
      mTit=mTit+1
      If ( mTit.gt.MxTitL ) Goto 993
      Title(mTit)=Line
      Goto 15
*---  Process the "TRANSform" command ---------------------------------*
2     Continue
      okay=.false.
      Inquire(File=trafil3,exist=okay)
      If(.not.okay) then
       Write(*,*) ' Option TRANSform needs file',trafil3
       stop 20
      Endif
*
      okay=.false.
      Inquire(File=orbfile,exist=okay)
      If(.not.okay) then
       Write(*,*) ' Option TRANSform needs file',orbfile
       stop 20
      Endif
      TRANS=.true. 
      Goto 100
*---  Process the "END of input" command ------------------------------*
    3 Continue
*     
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
*     
*
*----------------------------------------------------------------------*
*     OPEN 1EL INTEGRAL FILE (TRAONE) AND READ INPUT INFORMATION
*----------------------------------------------------------------------*
*     
*
      Call DaName(lutra1,trafil1)
      Call GsList(IOlist,8,iToc,64,Ecore,2,nSym,1,nBas,8,nOrb,8,nFro,8,
     * nDel,8,NameDu,800)
      iDisk=0
      Call dafile(lutra1,4,IOlist,iDum,iDisk)
      PotNuc = Ecore
      write(6,*) 'POTNUC=',potnuc
*
       write(*,*) ' nSym=',nSym
       write(*,*) ' nBas'
       write(*,999) nBas
       write(*,*) ' nOrb'
       write(*,999) nOrb
       write(*,*) ' nFro'
       write(*,999) nFro
       write(*,*) ' nDel'
       write(*,999) nDel
*     If(iPrint.ge.99) then
       write(*,*) ' ECor(1)',ECor(1)
       write(*,*) ' ECor(2)',ECor(2)
*     Endif
*
      Do iSym=1,nsym+1
        NTIL(iSym)=0
      End Do
*
*----------------------------------------------------------------------*
*     CALCULATE VALUES FOR 1 ST OUTPUT-RECORD, PACK THEM TO OCCUPIED
*     IRREPS    
*----------------------------------------------------------------------*
* 
      NOCC=0
      NBOX=0
      NCO=0
      KSUM=0
      NTIL(1)=0
      IORBS=0
      DO   IS=1,NSYM
        IORBS=IORBS+NBAS(IS)
        NCV=NBAS(IS)-NDEL(IS)
        IF(NCV.EQ.0) GOTO 134
        NOCC=NOCC+1
*
*-----# MO'S / OCCUPIED SYMMETRY --------------------------------------*
* 
        MJ(NOCC)=NORB(IS)+NDEL(IS)+NFRO(IS)
        KSUM=KSUM+NORB(IS)+NDEL(IS)+NFRO(IS)
*
*-----# FROZEN MO'S  / OCCUPIED SYMMETRY-------------------------------*
*
        KJ(NOCC)=NFRO(IS)
        NCO=NCO+NFRO(IS)
*
*----# VALENCE+VIRTUAL MO'S  / OCCUPIED SYMMETRY-----------------------*
* 
        NO=NORB(IS)
        LJ(NOCC)=NO
        NBOX=NBOX+NO
*
*----# SO'S / OCCUPIED SYMMETRY----------------------------------------* 
* 
        NJ(NOCC)=NBAS(IS)
*
*---- SYMMETRY-PACKING LABELS------------------------------------------* 
* 
        IOSYM(NOCC)=IS
*
*---- LAST MO IN EACH SYMMETRY-----------------------------------------*
*
134     IF(IS.NE.NSYM) NTIL(IS+1)=NTIL(IS)+NBAS(IS)
*
      End Do
*
*---- Calculate NIT(667)-----------------------------------------------*
*
      write(6,*) "IOSYM:",IOSYM
      Call MkNIT(NSYM,NIT,LJ,IOSYM,JAB)
*
*---- Calculate LG-----------------------------------------------------*
*
      IX=0                                                                      
      JX=0                                                                      
      KX=0                                                                      
      DO  I=1,NOCC                                                            
        LG(I)  =IX                                                               
        ITIL(I)=JX                                                              
        IBAL(I)=KX                                                              
        IX=IX+LJ(I)*(LJ(I)+1)/2                                               
        JX=JX+MJ(I)*NJ(I)                                                     
        KX=KX+MJ(I)                                                              
      End Do
*
*----------------------------------------------------------------------*
*     OPEN  STONEY file and write first record
*----------------------------------------------------------------------*
*     
      OPEN(UNIT=LUSTO,file=stofile,FORM='UNFORMATTED')
*
*     Set some parameters not necessarly used
*
      nod = 1000
      jod = 2000
*** Potnuc aus DEFAULT GESETZT
      VNUC =200.0d0
*
      if (nform.eq.0) then
        Write(LUSTO) NOCC,NOD,JOD,KSUM,NBOX,MJ,KJ,LJ,NJ,NSYM,NTIL,
     *    NBAL,IOSYM,JAB,IORBS,KNU,LSYM,NCOMP,Eo,co,VNUC,EZERO
      Else
       Write(LUSTO) NOCC,NOD,JOD,KSUM,NBOX,MJ,KJ,LJ,NJ,NSYM,NTIL,
     *    NBAL,IOSYM,JAB,IORBS,KNU,LSYM,NCOMP,co,VNUC,EZERO
        write(6,*) 'NBOx=',nbox
      End If
*
*----------------------------------------------------------------------*
*     Mache oiir                              
*----------------------------------------------------------------------*
*    
      ival = 0
      ibias(1) = 0
      Do i=1,nsym  
        nni = lj(i)      
        Do k=1,nni        
          ival = ival+1
          oiir(ival)=i
        End Do
        ibias(i+1) = ibias(i) + lj(i)
      End Do 
      DO i=1,8
        lj2(i) = (lj(i))
      End Do
*---- Write second  record of STONEY file------------------------------*
*
      if (nform.eq.0) then
        Write(LUSTO) NIT,LG,CFo,ICFo,IBAL,ITIL,MCOMP
      Else
        Write(LUSTO) NIT,LG,ICFo,IBAL,ITIL,MCOMP
      End If
*     
*---- DUMP FILES        -----------------------------------------------*  
*    
      Call Dumper(nit,oiir,idk,ibias,lj2)  
*     
*---- Close  STONEY file-----------------------------------------------*  
*     
      Close(LUSTO)
*
      If(Trans) then
*
*----------------------------------------------------------------------*
*     OPEN  TRANSF file 
*----------------------------------------------------------------------*
*     
      OPEN(UNIT=LUTRA3,file=trafil3,FORM='UNFORMATTED')
*
      Read(LUTRA3) iSOx,iPrimc,iTest
*
      If (iSOx.gt.MxAO) then
       write(*,*) ' iSOx.gt.MxAO: Increase MxAO'
       stop 20
      Endif
*
      If (iPrimc.gt.MaxPrim) then
       write(*,*) ' iPrimc.gt.MaxPrim: Increase MaxPrim'
       stop 20
      Endif
*
      If (iTest.eq.1) then
       write(*,*) ' File TRANSF contains already a transformation
     *matrix from primitives to MOs '
       stop 20
      Endif
*
      Read(LUTRA3) ((Ccord(iPrim,isel),iPrim=1,iPrimc),isel=1,6)
      Read(LUTRA3) ((Nlcomp(iPrim,isel),iPrim=1,iPrimc),isel=1,3)
      Read(LUTRA3) ((TraPriSO(iSO,iPrim),iSO=1,iSOx),iPrim=1,iPrimc)
*
      if(iPrint.gt.119) then
       call RecPrt1('TraPriso',trapriso,1,iSOx,1,iPrimc,MxAO,MaxPrim)
       call RecPrt1('Ccord',ccord,1,iPrimc,1,6,MaxPrim,6)
       call RecPrt2('Nlcomp',Nlcomp,1,iPrimc,1,3,MaxPrim,3)
      endif
*
*----------------------------------------------------------------------*
*     OPEN MO-FILE (INPORB) AND READ MO-COEFFICIENTS to CMO
*----------------------------------------------------------------------*
* 
      Call RdCMO 
*     
*---- Copy CMO to ROC--------------------------------------------------*
*
      FMT='(4E18.12)'
      KCMO  = 0
      NDIV  = 4
      icSO  =0
      iorbc =0
      icbas =0
*    
      DO 600 ISYM=1,NSYM
         DO 610 IORB=1,NORB(ISYM)-NFRO(ISYM)-NDEL(ISYM)
            iorbc=iorbc+1 
            DO 611 IBAS=1,NBAS(ISYM)
             ROC(iorbc,icBas+IBAS)= CMO(IBAS+KCMO)
611         CONTINUE
            KCMO=KCMO+NBAS(ISYM)
610      CONTINUE
         icSO=icSO+NBAS(ISYM)
         icBas=icBas+NBAS(ISYM)
600   CONTINUE
      if(iPrint.gt.119) then
       Call Recprt1(' Roc',Roc,1,iSOx,1,iSOx,MxAO,MxAO)
      Endif
*     
*---- Transformation from SO's to MO's---------------------------------*
*     
      If(Nbox.ne.iorbc) then
       write(*,*) ' Error in Form31: NBox.ne.iorbc '
       stop 20
      Endif 
      If(iSOx.ne.icSO) then
       write(*,*) ' Error in Form31: iSOx.ne.icSO'
       stop 20
      Endif 
*
      Do i=1,Nbox
       Do k=1,iPrimc
        tmp=0.0e0
        Do j=1,iSOx
         tmp=tmp+TraPriSO(j,k)*ROC(i,j)
        End Do
        TraPriMO(i,k)=tmp
       End Do
      End Do
      if(iPrint.gt.119) then
       Call Recprt1(' TRaPriMo',TRaPrimo,1,nbox,1,iPrimc,
     * MxAO,MaxPrim)
      Endif
*     
*---- Write Transformation matrix back to TRANSF-----------------------*
*     
      Rewind(Lutra3)
*
      Write(Lutra3) iSOx,iPrimc,1
      Write(Lutra3) ((Nlcomp(iPrim,isel),iPrim=1,iPrimc),isel=1,3),
     *((Ccord(iPrim,isel),iPrim=1,iPrimc),isel=1,6) 
      Write(LUTRA3) ((TraPriSO(iSO,iPrim),iSO=1,iSOx),iPrim=1,iPrimc)
      Write(LUTRA3) ((TraPriMO(iSO,iPrim),iSO=1,iSOx),iPrim=1,iPrimc)
*
      Close(Lutra3)
*
      Endif
*
      WRITE(*,*) 'Normal Ending Of Form31'
      Stop 'End Of Form31'
*
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
991   Write(*,*) ' '
      Write(*,'(6X,A)') '***  Error in subroutine FORM31 ***'
      Write(*,'(6X,A)') 'Failed in locating namelist'
      Stop 20
992   Write(*,*) ' '
      Write(*,'(6X,A)') '***  Error in subroutine FORM31 ***'
      Write(*,'(6X,A)') 'Unknwon command'
      Stop 20
993   Write(*,*) ' '
      Write(*,'(6X,A)') '***  Error in subroutine FORM31 ***'
      Write(*,'(6X,A)') 'Failed in title processing'
      Stop 20
999   format (8(2x,i3))
      End
