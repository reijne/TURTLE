*
* $Id: xpress.com,v 1.2 2000-10-26 15:38:01 psh Exp $
*
C	This is a common block provided by Express for certain message passing 
C	communication configurations. 
c
      integer nocare
      integer norder
      integer nonode
      integer ihost
      integer ialnod
      integer ialprc
      common/xpress/ nocare,norder,nonode,ihost,ialnod,ialprc
