      Subroutine GetH0(H,Ecore)
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Load one electron integrals                                      *
*                                                                      *
*     Calling parameters:                                              *
*     H   : core Hamiltonian matrix                                    *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*     
      Implicit Real*8 (A-H,O-Z)
*     
      Parameter( LuOne = 14)
*     
      Dimension H(*)
      Dimension IOList(20),iToc(64),nBas(8),nOrb(8),nFro(8),nDel(8)
      Character*4 Name(2,400)
*----------------------------------------------------------------------*
*     Start procedure:                                                 *
*     open the transformed one-electron integral file                  *
*----------------------------------------------------------------------*
      Call DaName(LuOne,'TRAONE')
*----------------------------------------------------------------------*
*     Set up the scatter/gather list and read table of contents        *
*----------------------------------------------------------------------*
      iDisk=0
      Call GSList(IOList,8,iToc,64,Ecore,2,nSym,1,
     *            nBas,8,nOrb,8,nFro,8,nDel,8,Name,800)
      Call DaFile(LuOne,4,IOList,iDum,iDisk)
*----------------------------------------------------------------------*
*     Determine the number of integrals (symmetry blocked)             *
*----------------------------------------------------------------------*
      NorbTT=0
      Do iSym=1,nSym
        NorbTT=NorbTT+Norb(iSym)*(Norb(iSym)+1)/2
      End Do
cdebug
      write(6,*) NorbTT,' symmetry blocked integrals read.'
*----------------------------------------------------------------------*
*     Load the core Hamiltonian matrix                                 *
*----------------------------------------------------------------------*
      iDisk=iToc(2)
      Call DaFile(LuOne,2,H,2*NorbTT,iDisk)
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
      Return
      End
