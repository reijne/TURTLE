      subroutine PropIntR1(nirrep)
      implicit real*8 (a-h,o-z)
      rewind(32)
      read(32)nirrep
      return
      end
      subroutine PropIntR2(nirrep,inirrep)
      implicit real*8 (a-h,o-z)
      dimension inirrep(nirrep)
      read(32)(inirrep(i),i=1,nirrep)
      return
      end
      subroutine PropIntR3(n1e,ncomp,iopx)
      implicit real*8 (a-h,o-z)
      rewind(32)
      read(32)
      read(32)
      do i=1,ncomp-1
        read(32)
        read(32)
      enddo
      read(32)n1e,iopx
      return
      end
      subroutine PropIntR4(mult,nirrep)
      implicit real*8 (a-h,o-z)
      dimension mult(nirrep,nirrep)
      do i=1,nirrep
         read(32)(mult(i,j),j=1,nirrep)
      enddo
      return
      end
      subroutine PropIntR5(ncomp,n1e,bob)
      implicit real*8 (a-h,o-z)
      dimension bob(n1e)
c     rewind(32)
c     read(32)
c     read(32)
c     do i=1,ncomp-1
c        read(32)
c        read(32)
c     enddo
      read(32)bob
      return
      end
