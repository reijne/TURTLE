c
c     $Header: /c/qcg/cvs/psh/GAMESS-UK/g/tcgmsg/ipcv4.0/hpuxargs.f,v 1.1.1.5 2007-10-30 10:14:14 jmht Exp $
c
      integer function hpargc()
      hpargc = iargc() + 1
      end
      integer function hpargv(index, arg, maxlen)
      character*256 arg
      hpargv = igetarg(index,arg,maxlen)
      end
