c   non existent subroutines
      subroutine qppope
       stop 'call from dppop.f to qppope, missing'
      end
c ---
c     subroutine trial(a,b,c,d,e,f)
c      print *,a,b,c,d,e,f
c      stop 'call from scfnew.f to trial, missing'
c     end
c ---
      subroutine solqd(a,b,c,d,e,f)
       print *,a,b,c,d,e,f
       stop 'call from drfsub.f to solqd, missing'
      end
c ---
      subroutine solqde(a,b,c)
       print *,a,b,c
       stop 'call from drfsub.f to solqde, missing'
      end
c ---
      subroutine solext(a,b)
       print *,a,b
       stop 'call from drfsub.f to solext, missing'
      end
c ---
      subroutine solexp
       stop 'call from drfsub.f to solexp, missing'
      end
c ---
      subroutine itrf
       stop 'call from drfsub.f to itrf, missing'
      end
c ---
      subroutine omgait(a,b,c,d,e,f,g,h,i,j,k,l)
       print *,a,b,c,d,e,f,g,h,i,j,k,l
       stop 'call from drfctl.f to omgait, missing'
      end
