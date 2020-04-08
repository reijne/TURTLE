      subroutine iout(grad,nvarn,en,done)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param'
      integer en,endone,done
      dimension grad(nvarn),buf(numatm,numatm),endone(numatm)

      open(13,file='hessinfo')
      rewind 13

      do k=1,done
        read(13,*) endone(k)
        read(13,'(5f14.8)') (buf(k,j),j=1,nvarn)
c       write(*,'(5f14.8)') (buf(k,j),j=1,nvarn)
      enddo
      rewind 13
      do k=1,done
        write(13,*) endone(k)
        write(13,'(5f14.8)') (buf(k,j),j=1,nvarn)
c       write(*,'(5f14.8)') (buf(k,j),j=1,nvarn)
      enddo

      write(13,*) en
      write(13,'(5f14.8)') (grad(j),j=1,nvarn)

      close(13)

      return
      end
