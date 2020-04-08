*Process Opt(0) NoVector
*Deck RecPrt
      SUBROUTINE RECPRT1(LINE,A,I,il,J,jl,m,n)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(m,n)
      CHARACTER LINE*(*), FORMAT*20
      WRITE (*,*)
      WRITE (*,*) LINE
      WRITE (*,'(1X,A,I3,A,I3,A)') ' Size ( ',Il,' X ',Jl,' )'
      WRITE (*,*)
      L = Min(jl,6)
      WRITE (FORMAT,'(A,I2,A)') '(1X,',L,'(E12.6,1X))'
      DO 10 K = 1 ,il
         WRITE (*,FORMAT) (A(I+K-1,J+L-1),L=1,jl)
10    CONTINUE
      WRITE (*,*)
      RETURN
      END
