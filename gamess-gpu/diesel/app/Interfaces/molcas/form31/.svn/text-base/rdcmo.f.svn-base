      Subroutine RdCmo
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Read molecular orbitals from input                               *
*     WORD='MOLO': read from unit 5 in format 6F12.8 (not active)      *
*     WORD='LUMO': read from unit INPORB in the same format            *
*     WORD='JOBI': read from an old casscf interface unit JOBIPH.      *
*                                                                      *
************************************************************************
*     
      Implicit Real*8 (A-H,O-Z)
      Include  'arrays.inc'
      Include  'common.inc'
      okay=.false.
*----------------------------------------------------------------------*
*     Read MO coefficients from a formatted vector file                *
*----------------------------------------------------------------------*
        Inquire (File=orbfile,exist=okay)
        If ( okay ) Then
          Open (unit=Luorb,
     *          file=orbfile,
     *          status='old',
     *          form='formatted')
************************************************************************
*   Mod. K.Pfingst 3/6/93
          Call RdVec(Luorb,nSym,nBas,nOrb,Cmo,Occno,iAutoCut,VecTit)
************************************************************************
          Close (unit=Luorb)
        Else
          Goto 992
        End If
************************************************************************
      return
992   Write(*,*) ' '
      Write(*,'(6X,A)') '***  Error in subroutine RdCmo ***'
      Write(*,'(6X,A)') 'Failed to locate vector file'
      Stop 20
      End
