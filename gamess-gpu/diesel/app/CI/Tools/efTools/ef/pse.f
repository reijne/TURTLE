      program pse

C     print out the periodic table 
c     based on standard masses from block data!
c     ghost  3.98

      implicit double precision (a-h,o-z)
      include 'param'
      common /mass/ ams(107),dummy(107)
      common /elemts/ elemnt(107)
      dimension  amass(107)
      character*2 elemnt

      print*
      print*, 'THESE ARE THE MASSES OF THE MOST ABUNDANT ISOTOPES'
      print*
      print*, ' MIND:'
      print*, ' Sometimes, that is only an excess of some percent!'
      print*, ' Ref.: Handbook of Chemistry and Physics'
      print*, '       63rd edition 1982-1983'
      print*
      print*


      do i = 1,107
         if (mod(i,2).eq.0) then
            write(*,50) i,elemnt(i),ams(i),
     $           i+1,elemnt(i+1),ams(i+1)
         elseif (i.eq.107.or.i.eq.1) then
            write(*,60) i,elemnt(i),ams(i)
         endif
      enddo
 50   format (i3,' ',a3,':',f16.8,i3,'   ',a3,':',f16.8)
 60   format (i3,' ',a3,':',f16.8)
      

      stop
      end


