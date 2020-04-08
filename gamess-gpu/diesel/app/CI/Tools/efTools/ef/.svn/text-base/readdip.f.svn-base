      subroutine readdip (dipx,dipy,dipz)
C     ghost
C     . reads the cartesian components 
C       of dipole momenta from file no. 3 (grad.out)
      implicit none
     
      integer*4 maxSpalte
      character zeile*80,pattern*9,a*128
      real*8 dip(10)
      real*8 dipx,dipy,dipz

      open(unit=3,file='grad.out')                                              
      rewind (3)
      pattern = 'total(MP2'
 10   read (3,'(A80)',end=999) zeile
      if (index(zeile,pattern).ne.0) then
         read(3,'(A80)') zeile
         read(3,'(A80)') zeile
         call readl(128,zeile,dip,maxSpalte)
         dipx=dip(5)
         read(3,'(A80)') zeile
         call readl(128,zeile,dip,maxSpalte)
         dipy=dip(5)
         read(3,'(A80)') zeile
         call readl(128,zeile,dip,maxSpalte)
         dipz=dip(5)
         close(3)
C         write(*,*)dipx,dipy,dipz
         return
      endif
      goto10
 999  write(*,*) '****************************************************'
      write(*,*) '* WARNING!!! No dipole info found in file grad.out *'
      write(*,*) '*         I am proceeding anyway - beware!         *'
      write(*,*) '****************************************************'
      close(3)
      return
      end
      
