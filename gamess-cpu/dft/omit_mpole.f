      subroutine JDriver(ao_tag,adens,bdens,kma,kmb)
      implicit none
      integer ao_tag
      real*8 adens(*),bdens(*)
      real*8 kma(*),kmb(*)
c
      call caserr('no multipole dft options available')
      end

      real*8 function ovtest(shl1,shl2,rsquar)
      implicit none
      integer shl1,shl2
      real*8 rsquar
      call caserr('no multipole dft options available')
      ovtest = 0.0d0
      end

      subroutine BASATO_FILL(itag)
      integer itag
c - no action
      end

      subroutine exfit(ao_tag,natyp,kf_tag,nff,
     &                 apts,awpt,prpt,prwt,
     &                 fitmat,CKfit,indx,
     &                 bfn_val,abfn_val,
     &                 adens,bdens,kma,kmb,
     &                 memory_int,memory_fp,extwr_sw)
      implicit none
      integer ao_tag,natyp,kf_tag,nff
      logical extwr_sw
      real*8 adens(*),bdens(*)
C *Out variables
      real*8 kma(*),kmb(*)
C *Scratch space and pointers
      integer memory_int(*)
      real*8 memory_fp(*)
      real*8 apts(3,*),awpt(*)
      real*8 prpt(natyp,*)
      real*8 prwt(natyp,*)
      real*8 fitmat(nff,*)
      real*8 CKfit(*)
      integer indx(*)
      real*8 bfn_val(*),abfn_val(*)
c
      call caserr('no XF fit options available')
      end
