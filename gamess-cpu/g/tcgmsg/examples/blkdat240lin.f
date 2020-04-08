      block data
C$Id: blkdat240lin.f,v 1.1.1.1 2000-10-26 16:29:40 psh Exp $
      implicit double precision (a-h, o-z)
      include 'cscf.h'
c
c     initalize data in common ... clumsy but avoids code to read in data
c
c     line of sixteen be atoms 4.0 a.u. apart with 240 orbitals
c
c     have 9s functions on each center and simulate p's by having
c     s function at +- 1 in each of x, y, z
c
      data ax /    0.0d0,  4.0d0,  8.0d0, 12.0d0,
     $            16.0d0, 20.0d0, 24.0d0, 28.0d0,
     $            32.0d0, 36.0d0, 40.0d0, 44.0d0,
     $            48.0d0, 52.0d0, 56.0d0, 60.0d0/
      data ay /   16*0.0d0/
      data az /   16*0.0d0/
      data  q /16*4.0d0/, enrep/152.3666555456/
c
      data x /9*0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     $        9*4.0d0, 5.6d0,  2.4d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0,
     $        9*8.0d0, 9.6d0,  6.4d0, 8.0d0, 8.0d0, 8.0d0, 8.0d0,
     $       9*12.0d0,13.6d0, 10.4d0,12.0d0,12.0d0,12.0d0,12.0d0,
     $       9*16.0d0,17.6d0, 14.4d0,16.0d0,16.0d0,16.0d0,16.0d0,
     $       9*20.0d0,21.6d0, 18.4d0,20.0d0,20.0d0,20.0d0,20.0d0,
     $       9*24.0d0,25.6d0, 22.4d0,24.0d0,24.0d0,24.0d0,24.0d0,
     $       9*28.0d0,29.6d0, 26.4d0,28.0d0,28.0d0,28.0d0,28.0d0,
     $       9*32.0d0,33.6d0, 30.4d0,32.0d0,32.0d0,32.0d0,32.0d0,
     $       9*36.0d0,37.6d0, 34.4d0,36.0d0,36.0d0,36.0d0,36.0d0,
     $       9*40.0d0,41.6d0, 38.4d0,40.0d0,40.0d0,40.0d0,40.0d0,
     $       9*44.0d0,45.6d0, 42.4d0,44.0d0,44.0d0,44.0d0,44.0d0,
     $       9*48.0d0,49.6d0, 46.4d0,48.0d0,48.0d0,48.0d0,48.0d0,
     $       9*52.0d0,53.6d0, 50.4d0,52.0d0,52.0d0,52.0d0,52.0d0,
     $       9*56.0d0,57.6d0, 54.4d0,56.0d0,56.0d0,56.0d0,56.0d0,
     $       9*60.0d0,61.6d0, 58.4d0,60.0d0,60.0d0,60.0d0,60.0d0/
      data y /9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0/
      data z /9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0/
      data expnt /1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0/
      end
