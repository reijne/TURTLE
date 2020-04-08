      double precision function timer()
c
c     return the time since the last call to timer.
c
c     must be initialized by calling once and throwing away the
c     value 
c     ... use cpu time on multi-user machines
c     ... use elapsed time on dedicated or single user machines.
c
*mdc*if unix
*      real*4 dtime, tt(2)
*      timer = dble(dtime(tt))
*mdc*elseif tcgmsg
C$Id: timer.F,v 1.1.1.1 2006-10-03 12:33:30 jmht Exp $
      save mlast
      data mlast/0/
      m = mtime()
      timer = dble(m - mlast) * 0.01d0
      mlast = m
*mdc*endif
c
      end
