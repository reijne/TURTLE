************************************************************************
*                                                                      *
*     The present file includes the work space definition used by      *
*     the memory manager. Two versions are available. The code         *
*     encapsulated by comment cards corresponds to the MOLCAS-2        *
*     standard using a fixed size memory and codes need to be          *
*     recompiled to chance the max size. In contrast, MOLCAS-3         *
*     queries the available work space from the system by means        *
*     of the malloc procedure (C_language) and returns the size.       *
*     The Common/WrkSpc/ is then used as an ancor for relative.        *
*                                                                      *
************************************************************************
*----------------------------------------------------------------------*
*     MOLCAS-2                                                         *
*----------------------------------------------------------------------*
*     Parameter ( mxMem = 4*1024*1024 )
*     Real*8    Work(mxMem)
*     Integer   iWork(2*mxMem)
*     Equivalence (Work,iWork)
*     Common /WrkSpc/Work
*----------------------------------------------------------------------*
*     MOLCAS-3.5                                                       *
*----------------------------------------------------------------------*
      Real*8    Work(8)
      Integer   iWork(8)
      Equivalence (Work,iWork)
      Common /WrkSpc/Work,mxMem
