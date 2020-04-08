c     $Header: /c/qcg/cvs/psh/GAMESS-UK/g/tcgmsg/ipcv4.0/testpf.f,v 1.1.1.5 2007-10-30 10:14:16 jmht Exp $
      character*60 fname
      call pbeginf
      fname = ' '
      write(fname,'(a,i3.3)') '/tmp/pfcopy.test',nodeid()
      call pfcopy(5, 0, fname)
      call pend
      end
