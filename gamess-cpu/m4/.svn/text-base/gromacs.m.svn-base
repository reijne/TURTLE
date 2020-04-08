c
c     Getfilenames for GROMACS. 
c     -------------------------
c
c     Both the taskfarm interface and the Chemshell interface expect
c     the getfilenames routine to exist to specify the names of the
c     input and output files. 
c
c     Here we follow the chemshell integration so we have to provide 
c     the corresponding routine. 
c
      subroutine getfilenames(infile,outfile,length)
      implicit none
      character*(*) infile, outfile
      integer length
c
      character*9 inname
      character*6 outname
      parameter (inname = 'gamess.in')
      parameter (outname = 'stdout')
c
      if (len(inname).le.length) then
        infile = inname
      else
        call caserr2('getfilenames: failed, infile too short')
      endif
c
      if (len(outname).le.length) then
        outfile = outname
      else
        call caserr2('getfilenames: failed, outfile too short')
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
c     For GROMACS we also need to provide the function default_heap
c     as that is what Chemshell provides instead of icoresize.
c     As we generally follow the Chemshell integration we have to
c     introduce this function, but we will provide it as a wrapper
c     around icoresize.
c
      integer function default_heap
      implicit none
c
      integer  icoresize
      external icoresize
c
      default_heap = icoresize()
c
      return
      end
