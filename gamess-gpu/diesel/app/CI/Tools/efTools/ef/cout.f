      subroutine cout(nvarn,grad,en,done)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param'
      integer en,endone,done
      dimension grad(nvarn),buf(numatm,numatm),endone(numatm)

      open(12,file='hesscinfo')
      rewind 12

      read(12,*) cstepp
      do k=1,done
        read(12,*) endone(k)
        read(12,'(5f14.8)') (buf(k,j),j=1,nvarn)
c       write(*,'(5f14.8)') (buf(k,j),j=1,nvarn)
      enddo
c     close(12)
c     open(12,file='hesscinfo')
      rewind 12
      write(12,*) cstepp
      do k=1,done
        write(12,*) endone(k)
        write(12,'(5f14.8)') (buf(k,j),j=1,nvarn)
c       write(*,'(5f14.8)') (buf(k,j),j=1,nvarn)
      enddo

      write(12,*) en
      write(12,'(5f14.8)') (grad(j),j=1,nvarn)

      close(12)

      return
      end
