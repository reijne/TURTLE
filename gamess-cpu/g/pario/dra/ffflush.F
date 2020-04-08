c $Id: ffflush.F,v 1.1.1.2 2002-12-23 14:11:40 psh Exp $
      subroutine ffflush(unit)
      integer unit
c
#ifdef CRAY
*     if(unit.eq.6)then
*       call flush(101)
*     else
*       call flush(unit)
*     endif
#elif defined(NX)
      call forflush(unit)
#elif !defined(AIX) && !defined(HPUX) && !defined(IFCLINUX)
      call flush(unit)
#endif
c
      end
