c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/plot.m,v $
c  $State: Exp $
c  
      subroutine aline(x1,y1,x2,y2)
c
c a line drawing routine using plot10 routines
c
      common/ploty/shiftx,shifty,xmill,chsze,iplot
      common/plotz/xmax,xmin,ymax,ymin,xscale,yscale
c
c convert to screen millimetres
      xr=(x1-xmin)*xscale
      yr=(y1-ymin)*yscale
      call movea(xr+shiftx,yr+shifty)
      xr=(x2-xmin)*xscale
      yr=(y2-ymin)*yscale
      call drawa(xr+shiftx,yr+shifty)
c     call alfmod
      return
      end
      subroutine andbit(un,pl,n)
c-------------------------------------------------------------
      logical un,pl
      dimension un(*),pl(*)
      common/junk2/ilifg(1)
      n1 = n - 1
      n2 = n - 2
      do 10 i = 1,n1
      ii=ilifg(i)
      do 10 j = 1,n2
      ij=ii+j
      un(ij)=un(ij).and.pl(ij)
 10   continue
      return
      end
      subroutine andlog(rmesh,plate,unused,h,n)
      logical plate,unused
c-------------------------------------------------------------
       dimension rmesh(*)
      dimension plate(*),unused(*)
      common/junk2/ilifg(1)
      n1 = n-1
      n2 = n1-1
      do 50 j=1,n2
      do 50 i=1,n1
      ix = n*j + i
50    unused(j+ilifg(i)) = (rmesh(ix).lt.h).and.(rmesh(ix+1).ge.h)
      call andbit(unused,plate, n)
      return
      end
      subroutine arch3d(grid8,grid,title,ngrid,
     *scamax,scamin,facmax,facmin,anga,angb,dist3d)
      REAL scamax,scamin,facmax,facmin,anga,angb
      REAL dist3d,grid8
      character*8 title(*)
      dimension grid(*),grid8(*)
      common/junk2/ilifg(200)
      common/ploty/xshift,yshift,xmill
      logical iswmx,iswmn
      data zero/0.0/
      call erase1
c     call tekden(0)
c
      gmax = scamax
      gmin = scamin
      fmax = facmax
      fmin = facmin
      ax = anga
      az = angb
      dist = dist3d
      iswmx=.true.
      iswmn=.true.
      ng=ngrid
      mm=ng-1
      ng2=ng/2
      ng1=ng2-1
      isq=ng*ng
      xng=ng
c iresol = 400 originally
      iresol=500
      iisq=(isq+1)/2
      do 406 loop=1,isq
 406  grid(loop) = grid8(loop)
c ... find local max and min and compare with(scamax,scamin)
      vmax=-1.0e37
      vmin=1.0e37
      do 3004 j=2,mm
      j4=(j-1)*ng
      k4=j4+ng
      l4=j4-ng
      do 3004 i=2,mm
      ii=i+1
      nt1=i-1
      pow=grid(j4+i)
      pow1=grid(j4+ii)
      pow2=grid(j4+nt1)
      pow3=grid(k4+i)
      pow4=grid(l4+i)
       if(pow.lt.pow1.or.pow.lt.pow2.or.pow.lt.pow3.
     *or.pow.lt.pow4)go to 3005
c ... maxima
      if(pow.gt.gmax.or.pow.lt.vmax)go to 3004
      iswmx=.false.
      vmax=pow
      go to 3004
 3005 if(pow.gt.pow1.or.pow.gt.pow2.or.pow.gt.pow3.
     *or.pow.gt.pow4)go to 3004
c ... minima
      if(pow.lt.gmin.or.pow.gt.vmin)go to 3004
      iswmn=.false.
      vmin=pow
 3004 continue
      if(iswmn)vmin=gmin
      if(iswmx)vmax=gmax
      vmax=vmax*fmax
      vmin=vmin*fmin
c      write(6,3006)vmax,vmin
c 3006 format(/20x,'grid maxima fixed at ',f10.4,' units'//
c     *20x,'grid minima fixed at ',f10.4,' units')
      do 3002 i=1,isq
      pow=grid(i)
      if(pow.lt.vmax)go to 3003
      grid(i)=vmax
      go to 3002
 3003 if(pow.ge.vmin)go to 3002
      grid(i)=vmin
 3002 continue
      trunc=xng*1.5
      factor=trunc/abs(vmax-vmin)*0.833333
      do 30 i=1,isq
 30   grid(i)=grid(i)*factor
      xresol=iresol
      xllim=-10.0*xresol/512.0
      call limits(xllim,xllim,xresol,xresol,1.0,1.0)
      call movea(xshift+45.0,yshift+10.0)
      call wrtchr(title,8,10,0,0)
c
c     call ghost 'surplt' for perspective plot
c
      call surplt(grid,1,ng,ng,1,ng,ng)
      call picnow
      call frame
c     call tekden(1)
c     call scrscl(1,23)
      return
      end
      subroutine ark(xv,yv,m,n,clopen,gridno,thron,throff)
c-------------------------------------------------------------
      real*4 newsin,newcos
      logical contin
      integer clopen,fml,gridno,times
      common/iofile/iread,iwrit
      dimension xv(*),yv(*)
c  xv,yv are vectors of x- and y-coordinates resp.
c  of points on a contour
c  m,n point to lowest and highest relevant elements  resp.
c  of xv,yv.
c  clopen is 1 for closed contour and 2 for open contour
c  gridno is thhe size of grid in use
c  error checks follow
      if(n-m-1) 800,900,910
  800    write(iwrit,801)
  801    format('less than two points specified on contour')
         goto 802
  900    clopen=2
         goto 920
  910    if(clopen.ne.1)clopen=2
c
  920 z=12.0-gridno*0.1
      if(gridno.gt.100)z=2.0
      times =ifix(z)
      contin=.true.
      k     =m
      goto(1,2),clopen
    1   fml=2
        goto 7
    2   fml=1
c trigs provides initial values for oldsin and oldcos
    7 call trigs(xv,yv,m,n,k,clopen,fml,oldsin,oldcos)
      x=xv(m)
      y=yv(m)
c***
      xjm=x
      yjm=y
      dist=0.0
      iplot=1
c start of loop
   10 k=k+1
      if(k-n)100,150,200
  100    goto(250,102),clopen
  102    fml=2
         goto 250
  200    deltax=xv(m)-xv(n)
         deltay=yv(m)-yv(n)
         k=m
         d=xv(n)
         s=yv(n)
         contin=.false.
         goto 300
  150    goto(250,152),clopen
  152    fml=3
         contin=.false.
  250 deltax=xv(k)-xv(k-1)
      deltay=yv(k)-yv(k-1)
      d=xv(k-1)
      s=yv(k-1)
c  trigs provides new values for newsin and newcos
  300 call trigs(xv,yv,m,n,k,clopen,fml,newsin,newcos)
      e     =7-oldcos*newcos-oldsin*newsin
      sumcos=oldcos+newcos
      sumsin=oldsin+newsin
      f     =deltax*sumcos+deltay*sumsin
      g     =deltax*deltax+deltay*deltay
      tk    =(3/e)*((f*f+2*e*g)**0.5-f)
      vk    =1/tk
      gamma =3*vk*vk
      eta   =2*vk
      bk    =gamma*(eta*deltax-sumcos)
      ck    =gamma*(eta*deltay-sumsin)
      alpha =0.5*vk
      beta  =0.5*tk
      a     =-0.333333333*bk
      b     =(newcos-oldcos)*alpha+bk*beta
      p     =-0.333333333*ck
      q     =(newsin-oldsin)*alpha+ck*beta
      deltat=tk/times
      t     =0
c  a loop to generate an arc from 1 point of definition to
c  the next follows
      do 500 i=1,times
         t=t+deltat
c      nested mult. to save time
         x=(((a*t+b)*t+oldcos)*t+d)
         y=(((p*t+q)*t+oldsin)*t+s)
c***
         dist=dist+sqrt((x-xjm)**2 + (y-yjm)**2)
         if(iplot)502,502,504
 504        if(dist.lt.thron)go to 503
            iplot=-1
            dist=0.0
            go to 503
 502        if(dist.lt.throff)go to 501
            iplot=1
            dist=0.0
            go to 501
 503     call aline(xjm,yjm,x,y)
 501     xjm=x
         yjm=y
 500     continue
      oldsin=newsin
      oldcos=newcos
      if(contin) go to 10
 802  continue
      return
      end
      subroutine bdry(ms,pts,trues,n)
      integer pts,trues
      dimension pts(*),trues(*)
      dimension ms(*)
c get single lines list
      max = ms(1) +1
      maxbd=8*n-5
      now = 2
      ixa = ms(3)
      iya = ms(4)
      do 100 j = 2,max
      i = j
      if(i .eq. max) i=1
      ix = 2*i + 1
      ixb = ms(ix)
      iyb = ms(ix+1)
10    inowx = 2*now -1
      if(inowx.gt.maxbd)then
         call gsperr('s/r bdry: inowx is greater than maxbd')
         call fatal
         endif
      pts(inowx) = ixa
      pts(inowx+1) = iya
      now = now + 1
      if (ixb - ixa) 1,2,3
    2 if(iyb.eq.iya) go to 100
      goto 7
1     ixa = ixa - 1
      goto 7
3     ixa = ixa + 1
7     if (iyb - iya) 4,10,6
4     iya = iya - 1
      goto 10
6     iya = iya + 1
      goto 10
100   now = now -1
      pts(1) = now - 1
c  set logical array to -1
      lim = n*n
      do 200 i=1,lim
200   trues(i) = -1
c  set boundary points to positive number corr. to order
      max = pts(1)*2 +1
      do 300 i=3,max,2
      ixa = pts(i)
      iya = pts(i+1)
      ix = n*(iya-1) + ixa
c   first point a special case
300   trues(ix) =(i-1)/2
      trues(ix) = trues(ix) +lim
      return
      end
      subroutine box
c
c
      common/ploty/shiftx,shifty,xmill
c
      call colour(2)
      xp=shiftx+xmill+1.5
      yp=shifty+xmill-14.5
      call movea(xp,yp)
      call drawr(0.0,-20.0)
      call drawr(35.0,0.0)
      call drawr(0.0,20.0)
      call drawr(-35.0,0.0)
      call colour(1)
      return
      end
      subroutine cntour(rmesh8,work,title,ngrid,rho,npoint,ktype,nucl
     +,oexp,nuccol,sizenu,sizel2,pltlab,pltcro)
      REAL rmesh8,rho
      character *8 title
      character *16 zgrid
      integer q,pts
      logical*4 pltlab,pltcro
      logical nucl,oexp
      dimension m(10),zgrid(8)
      dimension rmesh8(*),work(*),title(*),rho(*)
      common/ploty/shiftx,shifty,xmill,chsze,iplot
      common/junk2/ilifg(200),height(50),
     *xv(1000),yv(1000),q(200),pts(1596)
c
c  the routine produces a contour map of size ngrid by ngrid
c  it assumes that the co-ordinate of the bottom left hand
c  corner is (1.0,1.0), and that the next pt along the bottom
c  line is at (2.0,1.0) etc.
c  hence the bottom right hand corner is (ngrid,1.0)
c  hence the top   right hand corner is (ngrid,ngrid)
c  hence the top   left hand corner is (1.0,ngrid)
c
c  the real array rmesh contains the heights in the following
c  order
c  (1,1) , (2,1) ,   ...     , (n,1) ,
c  (1,2) , (2,2) ,   ...     , (n,2) ,
c  (1,n) , (2,n) ,   ...     , (n,n)
c
c  i.e. the values along the bottom    line from left to right
c  i.e. the values along the 2nd bottomline from left to right
c  i.e. the values along the top      line from left to right
c
c
c
c  npoint gives the no. of contours to be drawn (max. of 50)
c  heights supplied in rho
c***********************************************************************
c   the only limitation is that ngrid is less than or equal to 200
c***********************************************************************
      data m(1),m(2),m(3),m(4),m(6),m(9)
     1 /4,0,1,1,1,1/
      data zgrid/'electron density',
     +           'amplitude       ',
     +           'atom-difference ',
     +           'potential       ',
     +           'mol.- difference',
     +           'unknown         ',
     +           'unknown         ',
     +           'transition dens.'
     +/
c  define square grid as boundary
c  now for GHOST pre-amble
      call picnow
      call erase1
c     call tekden(0)
c write the title and grid type
      call movea(shiftx+1.8,shifty-4.0)
      call wrtchr(title,8,10,0,0)
      call movea(shiftx+2.0,shifty+xmill+2.0)
      call wrtchr(zgrid(ktype),16,1,0,0)
      if(nucl.and.oexp)call box
      n=ngrid
      nsq=n*n
      nocont=npoint
      mkx=1
      mxline=1000
      mxopen=50
      mkr=-1
c convert grid from real*8 to real*4
      do 406 loop=1,nsq
 406  work(loop) = rmesh8(loop)
      do 409 i=1,nocont
 409     height(i)=rho(i)
      n2=n/2
       do 22 i=1,n2
      ii=n-i+1
      ia=(i-1)*n
      ib=(ii-1)*n
      do 22 j=1,n
      temp=work(ib+j)
      work(ib+j)=work(ia+j)
 22   work(ia+j)=temp
      m(5) = n
      m(7) = n
      m(8) = n
      m(10) = n
      n1=n-1
      n2=n-2
      do 23 i=1,n1
 23   ilifg(i)=(i-1)*n2
      irmesh=1
      itrues=irmesh+nsq
      iplate=itrues+nsq
      iunu  =iplate + n1*n2
      itop  =iunu   + n1*n2
c ...
c ..  rmesh    trues        plate       unused
c ..   n*n      n*n        (n-2)*(n-1)   (n-2)*(n-1)
c ..   r*4      i*4           l*4          l*4
c
      call bdry(m,pts,work(itrues),n)
      call tmplt(work(itrues),pts,work(iplate),kount,n)
      call mainco(work(irmesh),pts,work(itrues),work(iplate),
     *work(iunu),xv,yv,q,n,height,nocont,mkr,mkx,mxline,mxopen)
c
      if(nucl) call shonuc(oexp,nuccol,sizenu,sizel2,pltlab,pltcro)
      call picnow
c     call tekden(1)
c     call scrscl(1,23)
      call frame
      return
      end
      subroutine colour(icol)
c
c --- sets line colour index to icol
c     but doesn't influence colours of the texts through wrtchr
c
      dimension ichar(4)
      integer*4 igr
      common/grdata/igr
      data ichar/27,77,76,48/
c               esc  m  l  0
      ichar(4)=icol+48
c     call toutst(4,ichar)
      if(igr.ne.0) write(igr,*)'0 0 c ',icol
      return
      end
_IF1()      subroutine corchk(lword,iparam,itype)
_IF1()      character*8 type(4)
_IF1()      common/iofile/iread,iwrit
_IF1()      data type/'basic','grid','atomscf','plot'/
_IF1()      lwor=lword-iparam
_IF1()      if(lwor)1,2,2
_IF1() 1    lwor=-lwor
_IF1()      lwor=(lwor*8)/1024+1
_IF1()      write(iwrit,3)type(itype),lwor
_IF1() 3    format(//10x,'**** insufficient core assigned in ',a8,' mode'//
_IF1()     *10x,'**** additional ',i3,' kbytes required')
_IF1()      call fatal
_IF1() 2    lwor=(lwor*8)/1024
_IF1()c      write(iwrit,4)type(itype),lwor
_IF1()c 4    format(//10x,'main core not used in ',a8,' mode =',i6,
_IF1()c     *' kbytes')
_IF1()      return
_IF1()      end
      subroutine cross(x,y,sz)
c
c  x,y are coordinates in mm for the cross position
c  s   is the size, in mm, of the cross
c
      s=sz*0.5
      call movea(x-s,y)
      call drawa(x+s,y)
      call movea(x,y+s)
      call drawa(x,y-s)
c     call alfmod
      return
      end
      subroutine crprod(a,b,c)
c
c       calculates cross product:   a x b = c
c                                   -   -   -
      REAL a(3),b(3),c(3)
c
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end
      subroutine disco(cont,ncon,nlins)
      REAL  cont
c
       common/iofile/iread,iwrit
      dimension cont(*)
c
      if(ncon.le.0)then
         call gsperr('cont: there are no contours to display')
         goto 9999
         endif
      write(iwrit,10)
   10 format(/'current contour values :'/)
      ncols= nint(float(ncon)/float(nlins)+.5)
      do 12 j=1,ncols
         write(iwrit,13)
   13    format('     contour  value')
   12    continue
      do 20 i=1,nlins
         if(i.gt.ncon)goto9998
         write(iwrit,15)i,cont(i)
   15    format(/5x,i3,2x,f13.8)
         if(ncon.le.nlins)goto20
         ipn=i+nlins
         if(ipn.gt.ncon)goto20
         write(iwrit,25)ipn,cont(ipn)
   25    format(5x,i3,2x,f13.8)
         if(ncon.le.2*nlins)goto20
         ipn=i+2*nlins
         if(ipn.gt.ncon)goto20
         write(iwrit,25)ipn,cont(ipn)
   20    continue
 9998 write(iwrit,'(/)')
 9999 return
      end
      subroutine drawa(x,y)
c
c****** drawa(x,y) ***************************************
c
      integer*4 igr,ipi
      common/grdata/igr,ipi,curx,cury
      call join(x,y)
      if(igr.ne.0)then
      curx=x
      cury=y
      write(igr,10)x,y
10    format(2f11.4)
      endif
      return
      end
      subroutine drawr(dx,dy)
c
c****** drawr(dx,dy) *************************************
c
      integer*4 igr,ipi
      common/grdata/igr,ipi,curx,cury
      call line(dx,dy)
      if(igr.ne.0)then
      curx=curx+dx
      cury=cury+dy
      write(igr,10)curx,cury
10    format(2f11.4)
      endif
      return
      end
      subroutine drwbdy(pts)
c-------------------------------------------------------------
       integer pts
       dimension pts(*)
c   draw boundary
c   only one call of join for each straight line
      max = 2*pts(1) + 1
       ax=pts(3)
       ay=pts(4)
c***
      xjm=ax
      yjm=ay
       bx=pts(5)
       by=pts(6)
      dx = bx-ax
      dy = by-ay
      do 100 i=7,max,2
       cx=pts(i)
       cy=pts(i+1)
      if(((cx-bx).eq.dx).and.((cy-by).eq.dy)) go to 50
c***
      call aline(xjm,yjm,bx,by)
      xjm=bx
      yjm=by
      dx = cx-bx
      dy = cy-by
50    bx = cx
100   by = cy
      if(((ax-bx).eq.dx).and.((cy-by).eq.dy)) go to 60
c***
      call aline(xjm,yjm,bx,by)
      xjm=bx
      yjm=by
c***
  60  call aline(xjm,yjm,ax,ay)
      return
      end
      subroutine endtst(u,i,j,n)
c-------------------------------------------------------------
       logical u
       dimension u(*)
       common/junk2/ilifg(1)
      n1 = n - 1
      n2 = n - 2
      do 10 i = 1,n1
      ii=ilifg(i)
      do 10 j = 1,n2
      ij=ii+j
      if(u(ij))go to 30
 10   continue
      i = 0
      j = 0
 30   return
      end
      subroutine erase1
c
c****** erase1 *******************************************
c
      call erase
      return
      end
_IF1()      subroutine errmsg(i,j)
_IF1()      common/iofile/iread,iwrit
_IF1()      call gsperr('surface drawing routines call to errmsg')
_IF1()      write(iwrit,900) i,j
_IF1() 900  format(' error number',i5,'  value in error =',i5)
_IF1()      call fatal
_IF1()      end
      subroutine fatal
c
c    a routine to handle fatal errors during graphics
c    and tidy up afterwards so that the terminal will still
c    talk to you properly
c
      call gsperr('sorry but this is fatal .... goodbye ....')
c     call anmode
c     call tekden(1)
c     call finitt(0,0)
      call stoplt
 _IFN1(h)     stop
 _IF1(h)     call exit
      end
_IF(convex)
_IF(cio)
      subroutine filech
c
c --- check files
c
      implicit REAL (a-h,o-z)
      character*4 ydd,yed3
      character*132 filnam
      character*8 ztitf(11)
INCLUDE(common/sizes)
INCLUDE(common/discc)
      common/disc/isel,iselr,iselw,irep,ichek,ipos(maxlfn),
     1 iblksz(maxlfn)
      common/iofile/ ir,iw,ipun(18),nav
      common/ddnam/yed3,ydd(maxlfn),filnam(maxlfn)
      common/setfil/ndd
c
      do 50 isel=1,maxlfn
         call rdtpnm(ztitf)
         if(irep.lt.0) then
            call gsperr('problem in rdtpnm - isel = '//char(isel))
            return
            endif
         if(ipos(isel).le.0) goto 50
         ndd=ndd+1
         ydd(ndd)=yed(isel)
         j=1
         do 45 k=1,10
            filnam(ndd)(j:j+7)=ztitf(k)
            j=j+8
 45         continue
 50      continue
c
      return
      end
_ENDIF
_ENDIF
      subroutine flcont(rmesh,unused,q,trues,h,xv,yv,n,mxline)
c-------------------------------------------------------------
      integer trues,q
      logical unused
      common/iofile/iread,iwrit
      common/junk2/ilifg(200)
       dimension trues(*),q(*)
       dimension unused(*)
       dimension rmesh(*),xv(*),yv(*)
      line = 2
c  test if q empty
      if(q(3).eq.0) go to 500
c  no
      yv(1) = +2.0
      idiag = 1
      ng2=n*n
      mb = q(2)
      q(2) = mb - 1
      q(3) = q(3) -1
      mx = (mb-1) * 4+1
      i= q(mx+2)
      j= q(mx+3)
      l= q(mx)
      m= q(mx+1)
      ia= l - i
      ja= m - j
      if (ia.eq.-1) then
      unused(m-1+ilifg(l)) = .false.
      endif
c   inverse lin. interpolation
50    idi = n*(j-1) + i
      idl = n*(m-1) + l
      z = rmesh(idi)
      za= rmesh(idl)
      t = (z-h)/(z-za)
      xv(line)= float(i) +t* float(ia)
      yv(line)= float(j) +t* float(ja)
c   swing tail of diag line to right
      line = line + 1
      if(line.gt.mxline)go to 5001
      if(line.gt.3)go to 10
      if(i-l) 1,2,3
1     if(j-m) 4,2,6
3     if(j-m)6,2,4
4     l = i
      goto 7
6     m = j
7     ia = l-i
      ja = m-j
      goto 10
c   is bdry diag across 1st square
2     in = i+ja
      jn = j-ia
      idn = n*(jn-1)+in
      idl = n*(m -1)+l
      if(trues(idn).eq.(trues(idl)-1)) go to 28
      in = in+ia
      jn = jn+ja
      idn = n*(jn-1)+in
      idi = n*(j-1)+i
      if(trues(idn).eq.(trues(idi)+1)) go to 27
10    l = i+ja
      m = j-ia
      idl = n*(m-1)+l
      if(rmesh(idl).lt.h) go to 20
      l=l+ia
      m=m+ja
      idl = n*(m-1)+l
      if(rmesh(idl).lt.h) go to 30
      i = l
      l=l-ja
      j = m
      m = m+ia
      goto 20
30    i = i+ja
      j = j-ia
20    ia = l-i
      ja =m-j
      if (ia.eq.-1) then
      unused(m-1+ilifg(l)) = .false.
      endif
      idl = n*(m-1)+l
      idi = n*(j-1)+i
      left = trues(idl)
      irght= trues(idi)
c   test for last point
      if(irght.gt.0) go to 22
      if(left.gt.0) go to 23
      goto 50
22    if (left.gt.ng2) left=left-ng2
      if(left.eq.(irght+1)) go to 29
      if((left.eq.2).and.(irght.gt.ng2)) go to 29
c   test for diag from i,j
      in = i+ja+ia
      jn = j-ia+ja
      idiag = 1
      idn = n*(jn-1)+in
      left= trues(idn)
      if(left.gt.ng2) left = left -ng2
      if(left.eq.(irght+1)) go to 26
      if((left.eq.2).and.(irght.gt.ng2)) go to 26
      goto 50
c   test for diag from l,m
23    in = i+ja
      jn = j-ia
      idiag = 2
      idn = n*(jn-1)+in
      irght= trues(idn)
      if(left.gt.ng2) left = left -ng2
      if(left.eq.(irght+1)) go to 26
      if((left.eq.2).and.(irght.gt.ng2)) go to 26
      goto 50
c   next line may be diag boundary
26    idi = n*(j-1)+i
      idl = n*(m-1)+l
      z = rmesh(idi)
      za = rmesh(idl)
      t = (z-h)/(z-za)
      xv(line)= float(i) +t*float(ia)
      yv(line)= float(j) +t*float(ja)
      line = line + 1
      if(line.gt.mxline)go to 5001
      goto (27,28),idiag
27    idn = n*(jn-1)+in
      if(rmesh(idn).ge.h) go to 275
      l = in
      m = jn
      ia = l-i
      ja = m-j
      goto 29
275   i = in
      j = jn
      goto 20
285   l = in
      m = jn
      goto 20
28    idn = n*(jn-1)+in
      if(rmesh(idn).lt.h) go to 285
      i = in
      j=jn
      ia = l-i
      ja = m-j
29    idi = n*(j-1)+i
      idl = n*(m-1)+l
      z = rmesh(idi)
      za = rmesh(idl)
      t=(z-h)/(z-za)
      xv(line)= float(i) +t*float(ia)
      yv(line)= float(j) +t*float(ja)
      xv(1)= line
      go to 10000
c  test if finished
500   call endtst(unused,in,jn,n)
      if(in.gt.0) go to 65
      yv(1)= 0
      xv(1)= 1
      return
c   closed curve
65    unused(jn+ilifg(in))  = .false.
      yv(1)= 1.0
      i = in + 1
      l = in
      j = jn + 1
      m = j
      ia = -1
      ja = 0
68    idi = n*(j-1)+i
      idl = n*(m-1)+l
      z = rmesh(idi)
      za = rmesh(idl)
      t = (z-h) / (z-za)
      xv(line)= float(i) + t*float(ia)
      yv(line)= float(j) + t*float(ja)
      l = i + ja
      m = j-ia
      idl = n*(m-1)+l
      if(rmesh(idl).lt.h) go to 70
      l = l+ia
      m = m+ja
      idl = n*(m-1)+l
      if(rmesh(idl).lt.h) go to 80
      i = l
      l = l-ja
      j = m
      m = m + ia
      goto 70
80    i = i+ja
      j = j-ia
70    ia = l-i
      ja = m-j
      if(ia.eq.-1) then
      if(.not.unused(m-1+ilifg(l)))go to 75
      unused(m-1+ilifg(l)) = .false.
      endif
      line = line + 1
      if(line.gt.mxline)go to 5001
      goto 68
75    xv(1)= line
10000 return
 5001 write(iwrit,5002)h
 5002 format(//' error detected in constructing ',
     +f10.4,' contour'//)
      call fatal
      return
      end
      subroutine gamplt(work,acgrid,rmesh)
      implicit real*4 (a-h,o-z)
      parameter(mxcont=61)
      parameter(defsizenu=2.0,deftextsize=0.038)
      REAL re8,rmesh,work,acgrid
      REAL scamax,scamin,facmax,facmin
      REAL popocc,plane,anga,angb,dist3d,scalpt,cont8
      REAL czan,c,cont,savec,csaved
      character*132 filenm
      character*16 zgrid(8)
      character*8 zcom,ztitle(11),ztext,atlab,atlsav
      character*4 acomm,ch4,yed3,ydd,mode,ytrunc
      character*1 zzz
      logical*4 fail,pltlab,onlyp,pltcro
      logical sttyp,nucl,oexp
      integer ncon, isavec
c     integer system
INCLUDE(common/sizes)
c
      parameter(igraph=42,igradump=43)
c
      dimension atlsav(maxat)
      dimension csaved(3,maxat)
      dimension cont(mxcont),savec(mxcont)
      dimension rmesh(*),work(*),acgrid(*)
      common/startyp/sttyp
      common/ploty/shiftx,shifty,xmill,chsze,iplot
      common/iofile/ir,iw
      common/ddnam/yed3,ydd(maxlfn),filenm(maxlfn)
      common/junk/popocc(maxorb),plane(26)
      common/infoa/nat,ich(7),czan(maxat),c(3,maxat)
      common/junkc/zcom(24),atlab(maxat)
      common/setfil/ndd
      common/grdata/igfile,ipict,xdummy(2),icolor,textsize
c
      parameter(ncomm=50)
      character*4 comm(ncomm)
      data scalpt/627.53d0/
      data comm/'quit','rest','setd','file','cont','dcon','econ','gcon',
c                1      2      3      4      5      6      7      8
     +          'icon','rcon','help','plot','pl  ','scal','dist','view',
c                9      10     11     12     13     14     15     16
     +          'info','vrot','hrot','fmax','fmin','smax','smin','clea',
c                17     18     19     20     21     22     23     24
     +          'add ','sub ','curr','titl','coor','exit','nucl','labe',
c                25     26     27     28     29     30     31     32
     +          'inde','nato','mode','gecp','csav','cloa','vdef','acon',
c                33     34     35     36     37     38     39     40
     +          'mul ','mcon','mgri','agri','glib','gdum','atom','asav',
c                41     42     43     44     45     46     47     48
     +          'aloa','aswa'
c                49     50
     +         /
c
      data zgrid/'electron density',
     +           'amplitude       ',
     +           'atom-difference ',
     +           'potential       ',
     +           'mol.- difference',
     +           'unknown         ',
     +           'unknown         ',
     +           'transition dens.'/
c
      icolor=0
      sttyp=.true.
      textsize=deftextsize
      sizel2=deftextsize
      sizenu=defsizenu
      nuccol=2
      maxcon= mxcont
      ncon=0
      isavec=0
      ngrid=0
      nucl=.true.
      pltcro=.true.
      oexp=.false.
      pltlab=.false.
      yed3='    '
      ndd=0
c
c  don't use picture graphics by default
c  if it is used, then
c  the fort.42 contains several picture inputs, for each use of plot c
c  directive one. They can be separated by divide.c program - see help
c
       igfile=0
c
c
      ipict=0
      idump=0
c     now 0 pictures and 0 dumped grids
c
c give initial values to scale and view variables (3d-surf.)
      scamax=0.7d0
      scamin=-0.7d0
      facmax=1.2d0
      facmin=1.2d0
      angb=30.0d0
      anga=30.0d0
      dist3d=1.3d0
c
c     open up gridfile
c
      call iniplt
c
c --- write headings
c
      write(iw,1)
    1 format(///28x,'gamsplot'/28x,'========')
      write(iw,2)
    2 format(6x,'GAMESS contour and 3d-surface',
     +          ' plotting module.'//)
c
c      write(iw,3)
c    3 format(
c     +       )
c
      write(iw,4)
    4 format(//12x,' *** gamsplot ready for command input ***'
     +      //5x,'(for help type "help",',
     +        ' to leave the program type "quit".)'/)
_IF(convex)
_IF(cio)
c
c --- open up all files available
c
c----------------------------------------
       call filech
c----------------------------------------
_ENDIF
_ENDIF
c
c --- command read and process loop:
c
   10 continue
      write(iw,15)
   15 format('gmpt>> ')
      call input
      call inpa(ztext)
      acomm=ytrunc(ztext)
c     call inpa4(acomm)
      if(acomm(1:1).eq.' ') goto 10
      do 20 icomm=1,ncomm
         if(comm(icomm).eq.acomm)goto 30
   20    continue
      call gsperr('command not recognised')
      goto 10
c
c --- on icomm jump to appropriate command code.
c
   30 continue
      goto(100,200,300,400,500,600,700,800,900,1000,
     +     1100,1200,1200,1400,1500,1600,1700,1800,1900,2000,
     +     2100,2200,2300,2400,2500,2600,2700,2800,2900,100,3100,3200,
     +     3300 ,3400,3500,3601,3700,3800,3900,4000,4100,4200,4300,4400
     +     ,4500,4600,4700,4800,4900,5000)icomm
c
c --- 1. quit
c --- 30. exit
c
  100 continue
      write(iw,105)
  105 format(/20x,' *** quitting gamsplot ***'/)
      call stoplt
      go to 3600
c
c --- 2. restore (restore a grid from section isec of
c                 the current dumpfile)
c
  200 continue
      call inpi(isec)
      fail=.false.
      call grid(rmesh,ngrid,ztitle,ktype,isec,
     +          fail)
      if(fail) then
       write(iw,3003)
3003     format('rest: not successfull in restoring grid')
         goto 10
         endif
      nsq=ngrid*ngrid
c transition densities will have ktype=8 because 5,6,7 
c was something different
      if(ktype.eq.6) ktype=8
c
c
c if this is an electrostatic potential grid, it must be
c multiplied by 627.53 .
      if (ktype.eq.4) then
      call dscal(nsq,scalpt,rmesh,1)
c      write(iw,*)'electrostatic potential grid multiplied by 627.53'
          endif
c generate a default set of contours according to ktype
      call gerco(ktype,rmesh,ngrid,cont,ncon,maxcon,0)
      write(iw,*)
     * 'default contours according to grid type generated'
      goto 10
c
c --- 3. setd (set current dumpfile)
c
  300 continue
      ch4='    '
      write(iw,*)'set new dumpfile name'
      call inpa(ch4)
      call setd(ch4)
      goto 10
c
c --- 4. file (display file information)
c
  400 continue
      write(iw,410)
  410 format(/2x,'files available:'/5x,'lfn',4x,'filename')
      do 430 i=1,ndd
         write(iw,420)ydd(i),filenm(i)
  420    format(5x,a4,3x,a80)
  430    continue
      if(yed3.ne.'    ')then
         write(iw,440)yed3
  440    format(/' current dumpfile has lfn ',a4)
         endif
      write(iw,*)
      goto 10
c
c --- 5. cont (display contour values)
c
  500 continue
      write(iw,*)'displaying list of contour values'
      call disco(cont,ncon,maxcon)
      goto 10
c
c --- 6. dcon (delete a contour from list of contour values)
c
  600 continue
         if(ngrid.le.0)then
         call gsperr('dcon: a grid has not been restored')
         goto 10
         endif
      call inpi(icon)
      if(icon.gt.ncon.or.icon.le.0)then
         call gsperr('dcon: contour '//char(icon)//' does not exist')
         goto 10
         endif
      cont8=cont(icon)
      do 630 i=icon,ncon
         cont(i)=cont(i+1)
  630    continue
      ncon=ncon-1
      write(iw,635)icon,cont8
  635 format(/'dcon: contour deleted:',i3,f10.5/)
      goto 10
c
c --- 7. econ (enter a new set of contours by hand)
c
  700 continue
      if(ngrid.le.0)then
         call gsperr('econ: a grid has not been restored')
         goto 10
         endif
      write(iw,705)
  705 format(//'enter contour values. enter  e  to finish'//
     +         'contour',5x,'value')
      icon=1
  710 write(iw,715)icon
  715 format(5x,i2,'   >>')
      read(ir,*,err=799)cont(icon)
      icon=icon+1
      if(icon.le.maxcon) goto 710
  797 read(ir,*,err=799) dummy
      write(iw,796)
  796 format('econ: too many contours'/)
      goto 797
  799 ncon=icon-1
      write(iw,3004) ncon
3004  format(i4,' contours will be used'/)
      goto 10
c
c --- 8. gcon (generate default set of contours)
c
  800 continue
      if(ngrid.le.0)then
         call gsperr('gcon: a grid has not been restored')
         goto 10
         endif
      call gerco(ktype,rmesh,ngrid,cont,ncon,maxcon,0)
      write(iw,*) 'default contours generated'
      call disco(cont,ncon,maxcon)
      goto 10
c
c --- 9. icon (insert a new contour after the one specified)
c
  900 continue
      if(ngrid.le.0)then
         call gsperr('icon: a grid has not been restored')
         goto 10
         endif
      if(ncon.eq.maxcon)then
         call gsperr('icon: too many contours')
         goto 10
         endif
      call inpi(icon)
      if(icon.gt.ncon.or.icon.lt.0)then
         call gsperr('icon: contour '//char(icon)//' does not exist')
         goto 10
         endif
      icon=icon+1
      do 910 i=ncon,icon,-1
         cont(i+1)=cont(i)
  910    continue
      ncon=ncon+1
      call inpf(cont(icon))
      write(iw,915)icon,cont(icon)
  915 format('new contour:',i3,3x,f10.5/)
      goto 10
c
c --- 10. rcon (replace the contour specified )
c
 1000 continue
      if(ngrid.le.0)then
         call gsperr('rcon: a grid has not been restored')
         goto 10
         endif
      call inpi(icon)
      if(icon.gt.ncon.or.icon.le.0)then
         call gsperr('contour '//char(icon)//' does not exist')
         goto 10
         endif
      call inpf(cont(icon))
      write(iw,1010)icon,cont(icon)
 1010 format('contour',i3,' replaced by value ',f10.5/)
      goto 10
c
c --- 11. help |
c
 1100 continue
      write(iw,*) ' ******************************************************************************'
      write(iw,*) '                Description of commands of Gamess utility plot'
      write(iw,*) ' =============================================================================='
      write(iw,*) ' '
      write(iw,*) '  1) quit: leaving of utility plot'
      write(iw,*) ' '
      write(iw,*) '  2) rest section_number: restore section from file with name set in setd'
      write(iw,*) '         command; data will be writen to "current register"; maximal section_number is 190'
      write(iw,*) '         This command has side effect - generates default contours according'
      write(iw,*) '         to type of restored grid. If you have developed your own contours, you'
      write(iw,*) '         should save them before rest and then load - see csav,cloa,acon,gcon,gecp'
      write(iw,*) ' '
      write(iw,*) '  3) setd filename: file with this name (4 characters) will be used for input'
      write(iw,*) ' '
      write(iw,*) '  4) file: display information about current dumpfile'
      write(iw,*) ' '
      write(iw,*) '  5) cont: print values of contours'
      write(iw,*) ' '
      write(iw,*) '  6) dcon contour_number: delete contour with this number'
      write(iw,*) ' '
      write(iw,*) '  7) econ <LF>value<LF> ... e<LF>: values for contours from this list will be '
      write(iw,*) '       used instead of default. <LF> means new line, e ends this list'
      write(iw,*) '       see also gcon,gecp,csav,cloa,acon'
      write(iw,*) ' '
      write(iw,*) '  8) gcon: sets default values for contours according to type of grid, which was'
      write(iw,*) '           restored from dumpfile. This procedure is called automatically whenever you restore grid'
      write(iw,*) '       see also rest,gecp,csav,cloa,acon'
      write(iw,*) ' '
      write(iw,*) '  9) icon contour_number value: insert contour with value after contour no...'
      write(iw,*) ' '
      write(iw,*) ' 10) rcon contour_number value: change value of contour no. ...'
      write(iw,*) ' '
      write(iw,*) ' 11) help: prints this file'
      write(iw,*) ' '
      write(iw,*) ' 12) plot picture_type: c=cont => plot contours; s=surf => plot 3-dimensional'
      write(iw,*) '        picture; grid information in "current register" is used - see also'
      write(iw,*) '        add,sub,clear,current directives'
      write(iw,*) ' '
      write(iw,*) ' 13) pl:   the same as plot, but usage of number 13 is not recomended'
      write(iw,*) ' '
      write(iw,*) ' 14) scale: prints functions cut-off values for surface plot'
      write(iw,*) ' '
      write(iw,*) ' 15) dist value: new value for view distance in surface plotting'
      write(iw,*) ' '
      write(iw,*) ' 16) view: prints information about point of view for surface plot'
      write(iw,*) ' '
      write(iw,*) ' 17) info: display information about grid'
      write(iw,*) ' '
      write(iw,*) ' 18) vrot angle: enter new value for x-rotation of 3d-surf.'
      write(iw,*) ' '
      write(iw,*) ' 19) hrot angle: enter new value for z-rotation of 3d-surf.'
      write(iw,*) ' '
      write(iw,*) ' 20) fmax value: new value for facmax'
      write(iw,*) ' '
      write(iw,*) ' 21) fmin value: new value for facmin'
      write(iw,*) ' '
      write(iw,*) ' 22) smax value: new value for scamax'
      write(iw,*) ' '
      write(iw,*) ' 23) smin value: new value for scamin'
      write(iw,*) ' '
      write(iw,*) ' 24) clear: clear the "accumulator register"'
      write(iw,*) ' '
      write(iw,*) ' 25) add: adds "current reg." to accumulator'
      write(iw,*) ' '
      write(iw,*) ' 26) sub: subtracts "current reg." from accumulator'
      write(iw,*) ' '
      write(iw,*) ' 27) current <LF> new_title: copy accumulator to current register and give new'
      write(iw,*) '         title to this grid'
      write(iw,*) '     This command has a side effect. It is supposed to be used for plotting'
      write(iw,*) '     of differential density maps and so it generates new contours with'
      write(iw,*) '     minus values of now existing contours, if old contours have no minus '
      write(iw,*) '     values else it lets the old contours. If you want to plot total densi-'
      write(iw,*) '     ty map later with your own contours, you should save them before using'
      write(iw,*) '     of this directive - see also csav,cloa,rest'
      write(iw,*) ' '
      write(iw,*) '     '
      write(iw,*) ' 28) title <LF> new_title: give new title to grid in "current register"'
      write(iw,*) ' '
      write(iw,*) ' 29) coord: print values of atomic coordinates'
      write(iw,*) ' '
      write(iw,*) ' 30) exit: the same as quit'
      write(iw,*) ' '
      write(iw,*) ' 31) nuclei [star|cros|nocros] on - off - explain - set colour size:'
      write(iw,*) '                      do you want to plot nuclei?'
      write(iw,*) '                      this command plots positions of nuclei in the plane.'
      write(iw,*) '                  Also colour and relative size can be setted.'
      write(iw,*) '                  By default is on. Explain will draw frame with following'
      write(iw,*) '                  explanation to each picture:'
      write(iw,*) '                  nuclei:'
      write(iw,*) '                  * in plane (or ''x in plane'' if cros was specified)'
      write(iw,*) '                  + projected'
      write(iw,*) ' '
      write(iw,*) ' 32) label on - off - set rel_size: switch atom labels (default off);'
      write(iw,*) '                             relative size can be setted (default 1.0)'
      write(iw,*) ' '
      write(iw,*) ' 33) index: print index of grids on current dumpfile'
      write(iw,*) ' '
      write(iw,*) ' 34) natom integer_number: set number of atoms'
      write(iw,*) ' '
      write(iw,*) ' 35) mode format_a4: set if ghost or internal surf routine should be used'
      write(iw,*) ' '
      write(iw,*) ' 36) gecp: generation of default contours for ecp calculation'
      write(iw,*) ' '
      write(iw,*) ' 37) csav: save current values of contours - see also rest,current'
      write(iw,*) ' '
      write(iw,*) ' 38) cloa: load current values of contours from saved ones - see also rest,curr'
      write(iw,*) ' '
      write(iw,*) ' 39) vdef: generate default view for surface plot'
      write(iw,*) ' '
      write(iw,*) ' 40) acon number: generate automatical contour values for given number of con-'
      write(iw,*) '                  tours with positive values. If there are negative values in'
      write(iw,*) '              the grid, it will generate also contours for negative values'
      write(iw,*) '              and the number of all contours will be 2*number+1.'
      write(iw,*) '              This procedure was tested only on several total and differen-'
      write(iw,*) '              tial electron density maps, but it should be suitable for'
      write(iw,*) '              every grid. It is not called automatically, because in this '
      write(iw,*) '              procedure error can occur for very "wild" grid.'
      write(iw,*) '              This procedure takes into account only positive values in grid'
      write(iw,*) '              and supposes that part of grid with negative values (if exists)'
      write(iw,*) '              is similar to the positive one. You can try to multiply the'
      write(iw,*) '              grid by -1 using mgri command. You can also try to multiply'
      write(iw,*) '              the grid by any constant, because the region of +/- (0.005-0.1)'
      write(iw,*) '              is supposed to be most important and is emphasized by smaller'
      write(iw,*) '              step of contour values.'
      write(iw,*) '                   Recomended numbers of contours:'
      write(iw,*) '                      total density - 10-17'
      write(iw,*) '                      difer density -  6-13'
      write(iw,*) '              See also: gcon,gecp,econ,mgri,agri'
      write(iw,*) ' '
      write(iw,*) ' 41) mul: multiplies accumulating grid by current grid'
      write(iw,*) ' '
      write(iw,*) ' 42) mcon value: multiply all contours by a given value'
      write(iw,*) ' '
      write(iw,*) ' 43) mgri value: multiply grid by constant - see also acon'
      write(iw,*) ' '
      write(iw,*) ' 44) agri value: add constant value to grid'
      write(iw,*) ' '
      write(iw,*) ' 45) glib [set colour size]:'
      write(iw,*) '             will create input for program picture.'
      write(iw,*) '             Colour of texts and their relative size can be setted.'
      write(iw,*) '             It works only for contour plots, not for surface plots!'
      write(iw,*) '             This input will be fort.42 . Don''t use "file ftn042 ftn042 keep"'
      write(iw,*) '             because in this way specified files must not be used for'
      write(iw,*) '             formatted output - they are OPENed for unformatted IO.'
      write(iw,*) '             To get inputs for each picture use command ''divide fort.42'''
      write(iw,*) '             Then the input file will be picture1 ,... This files are'
      write(iw,*) '             readable using vi, for example. To put them on screen use'
      write(iw,*) '             script /usr/local/dopict and .../toscreen or toplotter:'
      write(iw,*) '             example follows ...'
      write(iw,*) '             divide fort.42'
      write(iw,*) '             rm fort.42'
      write(iw,*) '             dopict 3'
      write(iw,*) '             toscreen'
      write(iw,*) '             dopict 1'
      write(iw,*) '             toplotter'
      write(iw,*) '             ...'
      write(iw,*) '             Default is off and pictures are plotted using ''ghost library'''
      write(iw,*) '             to file GAMESS_GRID.'
      write(iw,*) '             The advantage of this new possibilities is, that you can put'
      write(iw,*) '             any new text in the pictures and modify them using program'
      write(iw,*) '             picture. Finaly, in the title after graphics directive you can'
      write(iw,*) '             use super&sub scripts and greek letters (because it goes'
      write(iw,*) '             through program picture) in TeX like format. Use command'
      write(iw,*) '             ''picture -H'' or ''man picture'' to get more information.'
      write(iw,*) '             (The program "picture" is local at chemie at FU Berlin. Contact'
      write(iw,*) '             its author, Jiri Pittner - jiri@hpsiepsi.chemie.fu-berlin.de,'
      write(iw,*) '             if you are interested.)'
      write(iw,*) ' '
      write(iw,*) ' 46) gdump   writes in fort.43 ascii representation of current grid'
      write(iw,*) '             so that you can use other plotting programs to process it'
      write(iw,*) '             use command ''divide fort.43'' to split the file in several'
      write(iw,*) '             numbered pieces, each of them corresponds to one successfull'
      write(iw,*) '             use of gdump directive. The data format is directly useable'
      write(iw,*) '             for public domain "gnuplot" program, which can convert it'
      write(iw,*) '             to PostScript, X-window display or plenty of old fashioned'
      write(iw,*) '             devices.'
      write(iw,*) '             (The utility divide is trivial C program (should be'
      write(iw,*) '             distributed with Gamess), because the fort.43 is very simple:'
      write(iw,*) ' '
      write(iw,*) '             filename1'
      write(iw,*) '             1st line of file1'
      write(iw,*) '             ...'
      write(iw,*) '             last line of file1'
      write(iw,*) '             ASCII(1)'
      write(iw,*) '             filename2'
      write(iw,*) '             ...'
      write(iw,*) ' '
      write(iw,*) ' '
      write(iw,*) ' '
      write(iw,*) ' 47) atoms   on next natom lines follows table: label x y z'
      write(iw,*) '             sets coordinates and label for all atoms'
      write(iw,*) '             can be used if there is a bug in restoring atomic coordinates'
      write(iw,*) '             from dumpfile. To change number of atoms use command natom'
      write(iw,*) ' '
      write(iw,*) ' 48) asave   saves defined coordinates so that the restore does not destrou them'
      write(iw,*) ' '
      write(iw,*) ' 49) aload   loads saved coordinats in the buffer used by plot'
      write(iw,*) ' '
      write(iw,*) ' '
      write(iw,*) ' 50) aswap ix : swaps orientation of ix axis for atomic coordinates'
      write(iw,*) '                ix=1,2,3'
      write(iw,*) ' '
      write(iw,*) ' '
      write(iw,*) ' '
      write(iw,*) ' *******************************************************************************'
      write(iw,*) '                              example of input'
      write(iw,*) ' ==============================================================================='
      write(iw,*) ' '
      write(iw,*) ' cd $TMPDIR'
      write(iw,*) ' rm *       because there can be old GAMESS_GRID and new pictures would be added'
      write(iw,*) '        at the end of it'
      write(iw,*) ' cp /work1/store/eda/user/saved.ed3 ed3'
      write(iw,*) ' gamess << !'
      write(iw,*) ' file ed3 ed3 keep'
      write(iw,*) ' time 10'
      write(iw,*) ' utility'
      write(iw,*) ' plot'
      write(iw,*) ' help'
      write(iw,*) ' setd ed3'
      write(iw,*) ' glib'
      write(iw,*) ' rest 100'
      write(iw,*) ' gdump'
      write(iw,*) ' clear'
      write(iw,*) ' add'
      write(iw,*) ' acon 13'
      write(iw,*) ' csav'
      write(iw,*) ' plot c'
      write(iw,*) ' rest 110'
      write(iw,*) ' sub'
      write(iw,*) ' cloa'
      write(iw,*) ' plot c              plot with the same contours as previous sect. for comparisn'
      write(iw,*) ' curr'
      write(iw,*) '  Differential density 1^1A_{1g} - 2^1A_{1g}'
      write(iw,*) ' acon 9'
      write(iw,*) ' plot c'
      write(iw,*) ' plot s'
      write(iw,*) ' quit'
      write(iw,*) ' exit'
      write(iw,*) ' !'
      write(iw,*) ' cp GAMESS_GRID /mnt/user/'
      write(iw,*) ' cp fort.42 /mnt/user/'
      write(iw,*) ' cp fort.43 /mnt/user/'
      goto 10
c
c --- 12. plot (start plotting contours or 3d-surface)
c
 1200 continue
      ch4='    '
      if(ngrid.le.0)then
         call gsperr('plot:a grid has not been restored')
         goto 10
         endif
      call inpa(ch4)
      if(ch4.eq.'    ')then
         call gsperr('plot: missing plot argument')
         goto 10
         endif
       ipict=ipict+1
       if(igfile.ne.0)then
       if(ipict.lt.10)then
        write(igfile,810)ipict
       else
        if(ipict.lt.100)then
         write(igfile,8100)ipict
        else
         write(igfile,81000)ipict
        endif
       endif
810    format('picture',i1)
8100   format('picture',i2)
81000  format('picture',i3)
       write(igfile,*)'0 0 m'
       write(igfile,*)'240 180 m'
      endif
      if(ch4.eq.'cont')goto 1210
      if(ch4.eq.'c   ')goto 1210
      if(ch4.eq.'surf')goto 1300
      if(ch4.eq.'s   ')goto 1300
      call gsperr('plot: argument '//ch4//' not recognised')
      ipict=ipict-1
      goto 10
 1210 continue
c plot contours :
      if(ncon.le.0)then
         call gsperr('plot: contour values have not been set')
         ipict=ipict-1
         goto 10
         endif
      write(iw,3000) ipict
3000  format('plotting contours to picture no.',i3,
     +       ' in GAMESS_GRID (appending to it)'/)
      call cntour(rmesh,work,ztitle,ngrid,
     + cont,ncon,ktype,nucl,oexp,nuccol,sizenu,sizel2,pltlab,pltcro)
c
c  end picture when my graphics is used
c
      if(igfile.ne.0)then
      zzz='b'-'a'
      write(igfile,8888)zzz
8888  format(a1)
      endif
      go to 10
 1300 continue
c 13. plot 3d surface
      write(iw,3001) ipict
3001  format('plotting surface to picture no.',i3,
     * ' in GAMESS_GRID (appending to the file)'/)
c save a copy of the grid in sgri - not needed now, when cntour and 
c arch3d use work for temporary storage
      call arch3d(rmesh,work,ztitle,ngrid,scamax,scamin,
     +            facmax,facmin,anga,angb,dist3d)
      goto 10
c
c --- 14. scale (set function cut-off values for surface plot)
c
 1400 continue
      if(ngrid.le.0)then
         call gsperr('scale: a grid has not been restored')
         goto 10
         endif
      write(iw,1405)scamax,facmax,scamin,facmin
 1405 format(/'function scaling :'
     +       /'   scamax =',f10.3
     +       /'   facmax =',f10.3
     +       /'   scamin =',f10.3
     +       /'   facmin =',f10.3/)
      goto 10
c
c --- 15. dist (enter new value for view distance - 3d-surf.)
c
 1500 continue
      call inpf(dist3d)
      goto 10
c
c --- 16. view (display surface plot view data)
c
 1600 continue
      if(ngrid.le.0)then
         call gsperr('view: a grid has not been restored')
         goto 10
         endif
      write(iw,1615)anga,angb,dist3d
 1615 format(/'view :'
     +       /'    vrot =',f10.3
     +       /'    hrot =',f10.3
     +       /'    dist =',f10.3/)
      goto 10
c
c --- 17. info (display grid information)
c
 1700 continue
      if(ngrid.le.0)then
         call gsperr('info: a grid has not been restored')
         goto 10
         endif
      ng1=ngrid-2
      side=plane(26)
      write(iw,1705)
 1705 format(//15x,'grid information :'/
     +         15x,'------------------')
      write(iw,1710)(ztitle(i),i=1,10),ng1,ng1,isec,zgrid(ktype),
     +              side,side
 1710 format(/' title       :  ',10a8/
     +        ' grid size   :',i4,' *',i4, ' points'/
     +        ' section no. :',i4/
     +        ' type        :  ',a16/
     +        ' dimensions  :',f8.3,' *',f8.3,' bohr')
      write(iw,1715)
 1715 format(' plane definition :')
      k=0
      do 1719 n=1,7,3
         k=k+1
 1719    write(iw,1720)k,(plane(i),i=n,n+2)
 1720    format(8x,'point',i2,' =',3f9.3)
c if this is a mol. difference plot, skip the orb.occ.
      if (ktype.eq.5) goto 1790
c      write(iw,1725)
c 1725 format(' orbital occupation :')
c      write(iw,1730)(popocc(i),i=1,32)
c 1730 format(8f9.3)
 1790 write(iw,'(//)')
      goto 10
c
c --- 18. vrot (enter a new value for x rotation of 3d-surf.)
c
 1800 continue
      call inpf(anga)
      goto 10
c
c --- 19. hrot (enter a new value for z rotation of 3d-surf.)
c
 1900 continue
      call inpf(angb)
      goto 10
c
c --- 20. fmax (set new value for facmax.)
c
 2000 continue
      call inpf(facmax)
      goto 10
c
c --- 21. fmin (set a new value for facmin.)
 2100 continue
      call inpf(facmin)
      goto 10
c
c --- 22. smax (set a new value for scamax.)
c
 2200 continue
      call inpf(scamax)
      goto 10
c
c --- 23. smin (set a new value for scamin.)
c
 2300 continue
      call inpf(scamin)
      goto 10
c
c --- 24. clear (clear the accumulating grid.) *mdp*
c
 2400 continue
      write(iw,2405)
 2405 format('clearing the accumulating grid.')
      call vclr(acgrid,1,nsq)
      goto 10
c
c --- 25. add (add the current grid to the accum.grid) *mdp*
c
 2500 continue
      if(ngrid.le.0)then
         call gsperr('add: a grid has not been restored')
         goto 10
         endif
      write(iw,2505)
 2505 format('adding the current grid to the accumulating grid.')
      call vadd(acgrid,1,rmesh,1,acgrid,1,nsq)
      goto 10
c
c --- 26. sub (subtract the current grid from the accum.grid) *mdp*
c
 2600 continue
      if(ngrid.le.0)then
         call gsperr('sub: a grid has not been restored')
         goto 10
         endif
      write(iw,2605)
 2605 format('subtracting the current grid from the accumulating grid.')
       call vsub(acgrid,1,rmesh,1,acgrid,1,nsq)
      goto 10
c
c --- 27. current (copy accum.grid to current grid) *mdp*
c
 2700 continue
      if(ngrid.le.0)then
         call gsperr('curr: a grid has not been restored')
         goto 10
         endif
      write(iw,2705)
 2705 format('copying the accumulating grid to the current grid.')
      call dcopy(nsq,acgrid,1,rmesh,1)
      write(iw,2720)
 2720 format('enter a title for the new current grid:'/'>>')
      read(ir,'(10a8)')(ztitle(j),j=1,10)
      ztitle(11)=' '
c set ktype to molecular difference plot.
      ktype=5
c this should blow info's mind .... no number of section to be printed
      isec=123456
c generate the new contours for the difference plot if
c  old contours were only positive
      onlyp=.true.
      do 2730, i=1,ncon
      if (cont(ncon).lt.0.0) onlyp=.false.
2730  continue
      if(onlyp.or.ncon.eq.0) then
      call gerco(ktype,rmesh,ngrid,cont,ncon,maxcon,0)
      write(iw,*) 'new contours for difference plot were generated '
      write(iw,*) '(as opposite values to old contours)'
      endif
      goto 10
c
c --- 28. title  (enter a new title for the current grid)
c
 2800 continue
      write(iw,2810)
 2810 format('enter new title for the current grid:'/'>>')
      read(ir,'(10a8)')(ztitle(j),j=1,10)
      ztitle(11)=' '
      goto 10
c
c --- 29. coords. (print atomic coordinates)
c
 2900 continue
      if(ngrid.le.0)then
         call gsperr('coords: a grid has not been restored')
         goto 10
         endif
      write(iw,2910)
 2910 format(/15x,'atomic coordinates (bohr) : '
     +       /5x,'atom',10x,'x',14x,'y',14x,'z')
      do 2920 i=1,nat
 2920    write(iw,2930)i,atlab(i),(c(j,i),j=1,3)
 2930    format(1x,i5,2x,a4,1x,3(f12.7,3x))
      write(iw,'(/)')
      goto 10
c
c --- 30. exit  -> quit
c
c --- 31. nuclei (reset plot-nuclei flag)
c
 3100 continue
      ch4='    '
      call inpa4(ch4)
      if(ch4.eq.'star')then
        sttyp=.true.
        call inpa4(ch4)
        endif
      if(ch4.eq.'cros')then
        sttyp=.false.
        call inpa4(ch4)
        endif
      if(ch4.eq.'nocr')then
        pltcro=.false.
        call inpa4(ch4)
        endif
      nucl = (ch4.ne.'off ')
      oexp = (ch4.eq.'expl')
      if(ch4.eq.'set ')then
        call inpi(nuccol)
        call inpf(re8)
        sizenu=defsizenu*re8
        write(iw,9898)nuccol,sizenu
9898    format('Setting colour',i3,' and size ',f10.5,
     +         ' for marks of positions of nuclei.'/)
      endif
      goto 10
c
c --- 32. label (switch atom labels on/off)
c
 3200 continue
      pltlab=.true.
      ch4='    '
      call inpa(ch4)
      if(ch4.eq.'off ')then
         pltlab=.false.
         write(iw,*)' Atoms will not be labeled'
         goto 10
         endif
      if(ch4.eq.'set ')then
      call inpf(re8)
      sizel2=deftextsize*re8
      write(iw,*)' Atoms will be labeled with size ',sizel2
      endif
      goto 10
c
c --- 33. index (index of grids on dumpfile)
c
 3300 continue
      if(yed3.eq.'    ')then
         call gsperr('index: a dumpfile has not been set')
         goto 10
         endif
      write(iw,*) 'print index of grids on dumpfile'
      call index1
      goto 10
c
c --- 34. natom (set the number of atoms. added sept 1986)
c
 3400 call inpi(nat)
      write(iw,*) 'number of atoms set to ',nat
      goto 10
c
c --- 35. mode (set if ghost or internal surf routine is used. sept 86)
c
 3500 call inpa4(mode)
      write(iw,*) 'mode setted to ',mode
      goto 10
c
c --- 36. gecp (generate default set of contours suitable for
c               ecp-electron density maps)
c
3601  continue
      if(ngrid.le.0)then
         call gsperr('gecp: a grid has not been restored')
         goto 10
         endif
      call gerco(6,rmesh,ngrid,cont,ncon,maxcon,0)
      write(iw,*) 'default contours suitable for ecp-calc. generated'
      call disco(cont,ncon,maxcon)
      goto 10
c
c --- 37. csave (save current contours to savec array)
c
3700  continue
      if(ncon.le.0) then
      call gsperr('csav: no contour values are defined')
      goto 10
      endif
      isavec=ncon
      do 3701, i=1,ncon
3701  savec(i)=cont(i)
      write(iw,3702)
3702  format('contours were saved')
      goto 10
c
c --- 38. cload (load current contours from savec array)
c
3800  continue
      if(isavec.le.0) then
      call gsperr('cloa: no contour values were saved')
      goto 10
      endif
      ncon=isavec
      do 3801, i=1,ncon
3801  cont(i)=savec(i)
      write(iw,3802)
3802  format('contours were loaded')
      goto 10
c
c --- vdef give initial values to scale and view variables (3d-surf.)
c
3900  continue
      write(iw,*)
     * 'vdef: default values set for scale and view variables'
      scamax=0.7d0
      scamin=-0.7d0
      facmax=1.2d0
      facmin=1.2d0
      angb=30.0d0
      anga=30.0d0
      dist3d=1.3d0
      goto 10
c
c --- acon  automatical generation of contours
c
4000  continue
      call inpi(ngenc)
      if (ngenc.lt.5) ngenc=5
      if (ngenc.gt.mxcont) ngenc=mxcont
      if(ngrid.le.0)then
         call gsperr('acon: a grid has not been restored')
         goto 10
         endif
      call gerco(7,rmesh,ngrid,cont,ncon,maxcon,ngenc)
      write(iw,*) 'automatic contours generated'
      call disco(cont,ncon,maxcon)
      goto 10
c
c --- mul  multiplication of two grids
c
4100  continue
      if(ngrid.le.0) then
       call gsperr('mul: a grid has not been restored')
       goto 10
       endif
      do 4101, i=1,nsq
4101  acgrid(i)=acgrid(i)*rmesh(i)
      write(iw,*)'accumulating grid was multiplied by current grid'
      goto 10
c
c --- mcon  multiply contour values by a given constant
c
4200  continue
      call inpf(contm)
      if(ncon.le.0) then
       call gsperr('mcon: no contours exist')
       goto 10
       endif
      if(contm.eq.0.0 .or. contm.eq.1.0) goto 10
      do 4201, i=1,ncon
4201  cont(i)=cont(i)*contm
      write(iw,4202)contm
4202  format('contours were multiplied by ',f15.9)
      goto 10
c
c --- mgrid  multiplication of grid by constant
c
4300  continue
      call inpf(const)
      if(ngrid.le.0) then
       call gsperr('mgrid: a grid has not been restored')
       goto 10
       endif
      if(const.eq.0.0 .or. const.eq.1.0) goto 10
      do 4301, i=1,nsq
4301  rmesh(i)=rmesh(i)*const
      write(iw,4302)const
4302  format('current grid was multiplied by constant value ',f15.9)
      goto 10
c
c --- agrid adition of constant to grid
c
4400  continue
      call inpf(const)
      if(ngrid.le.0) then
       call gsperr('agrid: a grid has not been restored')
       goto 10
       endif
      if(const.ne.0.0) then
       do 4401, i=1,nsq
4401   rmesh(i)=rmesh(i)+const
       write(iw,4402)const
4402   format('constant value ',f15.9,' was added to current grid')
      endif
      goto 10
c
c --- glibrary sets flag for use of my plotting programs
c     the generated file is very simple and can be easily converted
c     for any other program - following example clears it
c  this is the example (without the fortran's c, of course)
c    123. 125. m
c    456.1 457.8
c    0 0 c 45
c    200. 200. p 12 3 0.023 10 20 tREST OF THE LINE IS TEXT
c    the 4 lines above mean: move to [x,y]; draw to [x,y]; 
c    change pen to number[n] (in this case [x,y] is ignored)
c    draw some text on given [x,y] with color 12, 
c    position 3 ([x,y]=right upper corner of the text)
c    size 0.023.size of viewing window, cursive angle 10 deg., 
c    angle to x axis 20 degrees
c
c     you are interested in program which converts the file above
c     to postscript, x11 terminal, tektronix terminal, ...
c     contact jiri@hpsiepsi.chemie.fu-berlin.de
c
c
4500  igfile=igraph
      ch4='    '
      call inpa4(ch4)
      if (ch4.eq.'set ')then
        call inpi(icolor)
        call inpf(re8)
        textsize=deftextsize*re8
      endif
      write(iw,9897)icolor,textsize
9897  format(/'Input for program picture will be created,',
     +       ' text colour and size will be:',i3,f10.7/)
      goto 10
c
c --- gdump saves grid in ascii form into igradump fortran stream
c     is assumed to be used if you want to plot surface or contours 
c     by some other graphical software - like public domain 
c     GNUPLOT, for example
c
4600  if(ngrid.le.0) then
         call gsperr('gdump: a grid has not been restored')
         goto 10
         endif
      idump=idump+1
      if(idump.lt.10)then
        write(igradump,811)idump
      else
        if(idump.lt.100)then
          write(igradump,8101)idump
        else
          write(igradump,81001)idump
        endif
      endif
811   format('dumpgrid',i1)
8101  format('dumpgrid',i2)
81001 format('dumpgrid',i3)
c the dump itsself
      write(igradump,'(128a)')
     + '#this is dump of grid calculated by GAMESS graphical analysis',
     + ' section'
      write(igradump,822)(ztitle(i),i=1,10)
822   format('#title of this grid:@',11a8)
      write(igradump,823)zgrid(ktype)
823   format('#type of this grid:',a16)
      write(igradump,824)plane(26)
824   format('#side of the window in a.u.:',f8.4)
      write(igradump,899)ngrid
899   format('#size of this grid in points:',i5)
      write(igradump,'(128a)')
     + '#definition of the plane of the grid by xyz coords of its',
     + ' three points in a.u.'
      i=1 
      write(igradump,826)(plane(j),j=i,i+2)
      i=7
      write(igradump,826)(plane(j),j=i,i+2)
      i=4
      write(igradump,826)(plane(j),j=i,i+2)
826   format('# ',3(f12.7,3x))
      write(igradump,'(128a)')
     + '#geometry of the calculated system in cartesian coords in a.u.'
      write(igradump,890)nat
890   format('#number of atoms is:',i5)
      do 827,i=1,nat
827   write(igradump,828)atlab(i),(c(j,i),j=1,3)
828   format('# ',a8,4x,3(f12.7,3x))
      write(igradump,'(128a)')
     + '# the data of the grid follow in the form directly acceptable',
     + ' by GNUPLOT program'
c
      do 820,i=1,ngrid
      do 821,j=1,ngrid
821   write(igradump,830)rmesh(j+(ngrid-i+1)*ngrid)
830   format(e12.6)
      if(i.ne.ngrid)write(igradump,*)
820   continue
c closing mark
      zzz='b'-'a'
      write(igradump,8888)zzz
      goto 10
c
c
c      47. atoms - defines positions and labels of all atoms
c
4700  continue
      do 4701,iat=1,nat
      read(ir,*)atlab(iat),c(1,iat),c(2,iat),c(3,iat)
4701  continue
      write(iw,*)' atoms  defined'
      goto 10
c
c      48. asave - save atomic coordinates
c
4800  continue
      natsav=nat
      do 4810,i=1,nat
      atlsav(i)=atlab(i)
      csaved(1,i)=c(1,i)
      csaved(2,i)=c(2,i)
      csaved(3,i)=c(3,i)
4810  continue
      write(iw,*)' atomic coordinates have been saved'
      goto 10
c
c      49. aload - loads atomic coordinates
c
4900  continue
      nat=natsav
      do 4910,i=1,nat
      atlab(i)=atlsav(i)
      c(1,i)=csaved(1,i)
      c(2,i)=csaved(2,i)
      c(3,i)=csaved(3,i)
4910  continue
      write(iw,*)' atomic coordinates have been restored'
      goto 10
c
c      50. aswap - swaps atomic coordinates
c
5000  continue
      call inpi(iswap)
      if(iswap.lt.1 .or. iswap.gt.3) 
     +     call gsperr('wrong copordinate in aswap')
      do 5001,iat=1,nat
5001  c(iswap,iat)=-c(iswap,iat)
      write(iw,*)' atomic coordinates swapped in axis ',iswap
      goto 10
c
 3600 return
      end
      subroutine gerco(ktype,grid,ngrid,val,ncon,maxcon,nautoc)
      implicit REAL (a-h,o-z)
c  empirical parameters for automatic generation of contours
      parameter (mxcont=61,qmin=6.0,qmax=1.8,nstep=10,podmx=0.0043,
     *pdmin=0.15,rsmin0=2.5d-5,rsmax0=120.0,
     *dlimit=12.0,qratio=0.15,pnlog0=0.3,
     *qnlog0=0.34,snlog=0.75,r0=0.05)
      logical*4 exstmi,wasswap
      character*8 zname
      dimension grid(*),val(*)
      dimension difer(2*mxcont+1),difer2(2*mxcont+1)
      dimension val1(2*mxcont+1)
      dimension rhotra(12),rhopot(12),rhotot(16),rhomo(12),
     +          rhodif(12),ecptot(16)
      common/iofile/ir,iw
      data ecptot/
     *2.0d-2,1.0d-2,8.0d-3,6.0d-3,4.0d-3,2.0d-3,1.0d-3,8.0d-4,6.0d-4,
     *4.0d-4,2.0d-4,1.0d-4,8.0d-5,6.0d-5,4.0d-5,2.0d-5/
      data rhopot/
     *210.0d0,180.0d0,150.0d0,120.0d0,90.0d0,75.0d0,
     * 60.0d0, 40.0d0, 20.0d0, 10.0d0, 5.0d0, 2.0d0/
      data rhotot/
     *64.7837d0,16.1959d0,4.0490d0,1.0122d0,0.5061d0,
     * 0.2531d0, 0.1265d0,0.0633d0,0.0316d0,0.0158d0,
     * 0.0079d0, 0.0040d0,0.0020d0,0.0010d0,0.0005d0,0.0002d0/
      data rhomo/
     *1.0d0,0.5d0,0.25d0,0.125d0,0.0625d0,
     *0.03125d0,0.01562d0,0.00781d0,0.00391d0,
     *0.00195d0,0.00098d0,0.00049d0/
      data rhodif/
     *0.86910d0,0.43455d0,0.21727d0,0.10864d0,
     *0.05432d0,0.02716d0,0.01358d0,0.00697d0,
     *0.00339d0,0.00170d0,0.00085d0,0.00042d0/
      data rhotra/
     *2.0d-3,1.0d-3,8.0d-4,4.0d-4,3.0d-4,2.0d-4,1.0d-4,
     +8.0d-5,4.0d-5,2.0d-5,1.0d-5,5.0d-6/
      data dzero,m4/
     *0.0d0,4/
c
c
c
      nsq=ngrid*ngrid
      go to (1,2,3,4,55,66,77,88,2003),ktype
c        9 and so on = unknown types
c
c --- electron density
 1    nplot=16
      do 5 j=1,nplot
 5       val(j)=rhotot(j)
      do 20 itest=1,nsq
         if(grid(itest).lt.dzero)go to 3
 20      continue
      go to 2002
c --- amplitude
 2    nplot=25
      do 7 i=1,12
         top=rhomo(i)
         val(i)=top
 7       val(26-i)=-top
      val(13)=dzero
      go to 2002
c --- atom-difference or molecular-difference default
 3    nplot=25
      do 8 i=1,12
         top=rhodif(i)
         val(i)=top
 8       val(26-i)=-top
      val(13)=dzero
      go to 2002
c --- transition density default
 88   nplot=25
      do 888, i=1,12
      top=rhotra(i)
      val(i)=top
888   val(26-i)=-top
      val(13)=dzero
      go to 2002
c --- electrostatic potential
 4    nplot=25
      do 9 i=1,12
         top=rhopot(i)
         val(i)=top
 9       val(26-i)=-top
      val(13)=dzero
      go to 2002
c --- atom-difference or molecular-difference
c     when values for total density  were set before
55    if (ncon.le.0) go to 3
c      if nothing was set, use default
      exstmi=.false.
      do 56, i=1,ncon
      if (val(i).lt. 0.0) exstmi=.true.
56    continue
      if (exstmi) goto 3
c     when there were minus values, you probably wanted default
58    if (2*ncon+1.gt.maxcon) then
      nplot=maxcon
      if (nplot.eq.(nplot/2)*2) nplot=nplot-1
      ncnold=ncon
        ncon=(nplot-1)/2
        do 11, i=1,ncon
11        val(i)=val(i+ncnold-ncon)
      else
      nplot=2*ncon +1
      endif
      do 10 i=1,ncon
 10      val(nplot+1-i)=-val(i)
      val(ncon+1)=dzero
      goto 2002
c --- total density map for ecp calculation
c
 66   nplot=16
      do 866 i=1,nplot
 866     val(i)= ecptot(i)
      go to 2002
c
c --- automatic generation of contours
c
 77   continue
      write(iw,111) nautoc
111   format(20x,'Automatic generation of',i3,' contours'/)
      if(nsq.le.1) goto 55
c
c  if absol. minimum is in abs. value greater then abs. maximum,
c  mult. grid by -1 and after generation return it back
c
      tmin=0d0
      tmax=0d0
      do 1001,i=1,nsq
      g=grid(i)
      if(g.gt.tmax)tmax=g
      if(g.lt.tmin)tmin=g
1001  continue
      if(tmax.eq.tmin.and.tmax.eq.0d0) 
     + call caserr('gerco: cannot generate contours, zero grid')
      if(tmax.eq.0d0) write(iw,*)' grid has only negative values'
      if(tmin.eq.0d0) write(iw,*)' grid has only positive values'
      write(iw,1002) tmax,tmin
1002  format('Absolute maximum and minimum of the grid is: ',2f14.8)
      tmin=dabs(tmin)
      if(tmax.ge.tmin) then
        wasswap=.false.
        zname='positive'
      else
        wasswap=.true.
        zname='negative'
        do 1004,i=1,nsq
1004    grid(i)=-grid(i)
      endif
      write(iw,1003)zname
1003  format('Contour generation will be based on ',a8,
     +       ' grid points.'//)
c
c  main algorithm for generation of contours
c
      nminus=0
      nreson=0
      amax=-1.0
      amin=1d9
      aver=0.0
      averm=0.0
      avera=0.0
      raver=0.0
      do 78, i=1,nsq
      a=grid(i)
      if(a.le.0.0) then
      nminus=nminus+1
        else
        if(a.gt.amax) amax=a
      if(a.lt.amin) amin=a
      b=log10(a)
      aver=aver+b
      averm=averm+a
      avera=avera+a*a
      endif
78    continue
      resmin=rsmin0
      resmax=rsmax0
      if(nminus.gt.0) then
        resmin=3.0*resmin      
      resmax=0.5*resmax
      endif
      do 781, i=1,nsq
        a=grid(i)
      if(a.ge.resmin .and. a.le.resmax) then
        nreson=nreson+1
          raver=raver+a
        endif
781   continue
      totmax=amax
      totmin=amin
      if (amax.gt.resmax) amax=resmax
      if (amin.lt.resmin) amin=resmin
      absmax=amax
      absmin=amin
      amax=log10(amax)+0.1+0.003*nautoc**2
      amin=log10(amin)-0.0005*nautoc**2
      nplus=nsq-nminus
      aver=aver/float(nplus)
      averm=averm/float(nplus)
      avera=log10(avera/float(nplus))/2.0
      pdmina=pdmin+(avera-aver)/8.0
      if(pdmina.gt.0.35) pdmina=0.35
      pdmaxa=podmx+(avera-aver)/500.0
      if(pdmaxa.gt.0.015) pdmaxa=0.015
      write(iw,112)nsq,zname,nplus,zname,totmax,totmin,10**aver,
     + averm,10**avera
112   format(1x,'Some information about grid:'/
     *       1x,'number of points: ',i8/1x,
     *'number of points with ',a8,' value: ',i8/1x,
     *'  --- following information valid for absolute values of ',a8,
     *' points ---'/1x,
     *'absolute maximum:    ',f15.9/1x,
     *'absolute minimum:    ',f15.9/1x,
     *'geometrical average: ',f15.9/1x,
     *'arithmetical average:',f15.9/1x,
     *'kvadratical average: ',f15.9/
     *)
      aver=log10(averm)
c   if it is possible, use averages from reasonable values only
      if(nreson/float(nplus).ge. 0.15) then
         raver=raver/float(nreson)
c write(iw,*)'RAVER= ',raver
         aver=log10(raver)
      endif
c  information extracted, now cut very high/low peeks
      stepup=(amax-aver)/float(nstep)
      stpdwn=(amin-aver)/float(nstep)
      alim=aver
      do 79, i=1,nstep
      alim=alim+stepup
      elim=10**alim
      nhiold=nhi
      nhi=0
      do 791, j=1,nsq
      if(grid(j).ge.elim) nhi=nhi+1
791   continue
      if(float(nhi)/float(nplus).le.pdmaxa) goto 792
      if(i.ge.5 .and.(abs(nhiold-nhi)/float(nplus)).lt.
     *  (pdmaxa/5.0)) goto 792
79    continue
792   amax=alim
      alim=aver
      do 793, i=1,nstep
      alim=alim+stpdwn
      elim=10**alim
      nlo=0
      do 794, j=1,nsq
      if(grid(j).gt.0.0 .and. grid(j).le.elim) nlo=nlo+1
794   continue
      if(float(nlo)/float(nplus).le.pdmina) goto 795
793   continue
795   amin=alim
      a=log10(absmax)-0.02
      if(amax.gt.a) amax=a
      a=log10(absmin)+0.02
      if(amin.lt.a) amin=a
c   now calculate contour decrement
c      write(iw,*)'amin=,amax=',amin,amax
      p=amax-amin
      if(p.lt.qmax) then
      amax=amax+(qmax-p)/2.0
      p=qmax
      endif
      if(p.gt.qmin) then
      amax=amax-(p-qmin)/2.0
      p=qmin
      endif
      p=p/float(nautoc-1)
c   and make positive contours
      q=10**(-p/2)
      do 800, i=1,2*nautoc+1
      if(i.eq.1) then
      elim=10**amax/q
      else
      elim=elim*q
      endif
      val1(i)=elim
      nhiold=nhi
      nhi=0
      do 801, j=1,nsq
      if(grid(j).ge.elim) nhi=nhi+1
801   continue
      if(i.eq.1) nhiold=nhi
      difer(i)=(nhi-nhiold)/float(nplus)
800   continue
      df2max=0.0
      do 802, i=1,2*nautoc
      difer2(i)=(difer(i+1)-difer(i))/(1.0+nsq/1000.0)
      if (difer2(i).gt.df2max) df2max=difer2(i)
802   continue
      dif2av=0.0
      ndifpl=0
      do 806, i=1,2*nautoc
      difer2(i)=difer2(i)/df2max
      if(difer2(i).gt.0.0) then
      dif2av=dif2av+difer2(i)
      ndifpl=ndifpl+1
      endif
806   continue
      dif2av=dif2av/float(ndifpl)
      q=10**(-p)
      val(1)=10**amax
      rrr=log10(r0*averm)/2.0
c      write(iw,*)'rrr= ',rrr,' q=',q
      do 803, i=2,nautoc
      val(i)=val(i-1)*q
      a=log10(val(i))-rrr
      do 804, j=2*nautoc+1,1,-1
      if(val(i-1).le.val1(j)) goto 805
804   continue
805   j=j-1
      if(j.le.1) j=2
      if(j.gt.2*nautoc) goto 808
c    calculate perturbation of next contour
      if(difer2(j)*difer2(j-1).lt.0) then
      if((difer2(j-1)/dif2av).ge.dlimit) then
      val(i)=val(i)*10**(p*qratio)
      endif
      endif
c     normal perturbation
      pnlog=pnlog0
      qnlog=qnlog0
      if (nminus.ge.0) then
      pnlog=0.9*pnlog
      qnlog=0.95*qnlog
      endif
      if(i.gt.1) val(i)=val(i)*10**(p*pnlog*difer2(j-1))
808   if(i.gt.1) then
      val(i)=val(i)*10**(p*qnlog*(2*exp(-snlog*a*a)-1.0))
      if(val(i).ge. 0.95*val(i-1)) val(i)=0.95*val(i-1)
      endif
803   continue
c
c  put grid back to the original state
c
      if(wasswap)then
      do 1005,i=1,nsq
1005  grid(i)=-grid(i)
      endif
c
c  return value of contours
c
      ncon=nautoc
      nplot=nautoc
      if(nminus.ge.1.or.wasswap) goto 58
c
 2002 ncon=nplot
 2003 return
      end
      subroutine grid(rmesh,ngpts,ztitle,ktype,isec,
     +                ofail)
      implicit REAL (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
c this routine is adapted for new analb graphics code
c
c grid definition parameters
c
      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &             igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &             igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &             ngrid
c
c data calculation parameters
c
      common/dfcalc/cdata(5,mxcalc),ictype(mxcalc),icsect(mxcalc),
     &            icgrid(mxcalc),icgsec(mxcalc),ndata(mxcalc),
     &            icstat(mxcalc),icdata(5,mxcalc),ncalc
     &            iresec(mxrest),nrest
c
c labels and titles
c
      common/cplot/zgridt(10,mxgrid),zgrid(mxgrid),
     &             zcalct(10,mxcalc),zcalc(mxcalc),
     &             zplott(10,mxplot),zplot(mxplot)
c
      logical oerr
      dimension rmesh(*)
      dimension ztitle(*)
      common/iofile/ ir,iw,ipun(18),nav
      common/ddnam/yed3
      common/junk/popocc(maxorb),plane(26)
      common/junk2/kstart(7,mxshel),nshell,non,num(2)
      common/junkc/zhead(4),ztit(10)
     *, zcom(10),ztag(maxat)
c
c here was error - the last array in infoa should be cc, not c
c
      common/infoa/nat,ich(7),czan(maxat),cc(3,maxat)
INCLUDE(common/machin)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
INCLUDE(common/sector)
      data m14,m26,m255,itgrid/14,26,255,51/
      data m1/1/
      m110=10+maxat
      m1420=mxprim*5+maxat
      ngpts = 0
      call secfet(isect(491),m1,iblk,ofail)
      if(ofail)return
      iblk=iblk+lensec(m1420)
      call rdchr(zcom,m110,iblk,numdu)
      call reads(kstart,mach(2),numdu)
      nat=non
c
      if(yed3.eq.'    ')then
         call gsperr('dumpfile has not been defined')
         ofail=.true.
         goto 999
      endif
      if(isec.le.1.or.isec.gt.188)then
         call gsperr('invalid section number')
         ofail=.true.
         goto 999
      endif
      ierr = 1
      ibld1 = 1
      call rcasec(isec,icalc,numdu,ibld1,rmesh,.true.,ierr)
c there was memory error - core dumped, than disapeared when this 'write'
c was added - very strange. When it was removed, the problem reappeared!
c 
c
      write(6,*)' grid: restoring grid data'
      if(ierr.eq.1)then
         call gsperr('error retreiving data from dumpfile')
      else if(ierr.eq.2)then
         call gsperr('error retreiving grid from dumpfile')
      else if(igtype(icgrid(icalc)).ne.1)then
         call gsperr('grid must be 2D')
      else if(npt(1,icalc).ne.npt(2,icalc))then
         call gsperr
     + ('grid must have equal numbers of points in x and y directions')
      else
         igrid = icgrid(icalc)
         sidex = dist(geom(1,igrid),geom(4,igrid))
         sidey = dist(geom(1,igrid),geom(7,igrid))
         if(dabs(sidex-sidey).gt.1.0d-5)then
            call gsperr('grid must be square')
            goto 999
         endif
c
c convert grid parameters
c
         do 100 j = 1,10
            ztitle(j)=zcalct(j,icalc)
 100     continue
         ngpts = npt(1,icgrid(icalc))
         ktype = ictype(icalc)
         a1 = ngpts
         a2 = ngpts-2
         plane(26) = sidex*a2/a1
c
c  recalculate geom to plane for back compatibility
c  is used at least by gamplt and shonuc
c
c geom(10,...) does not seem to be the 4th corner !!!
c in one case it was equal to geom(1,...)
c
      do 110, j=1,3
      plane(j)=(geom(j+6,igrid)+geom(j+3,igrid))/2.0
c      plane(3+j)=geom(6+j,igrid)-geom(j,igrid)
c swap it on HP
      plane(3+j)=-geom(6+j,igrid)+geom(j,igrid)
      plane(6+j)=geom(3+j,igrid)-geom(j,igrid)
110   continue
c
c --- read atomic coordinates from dumpfile
c
         call rdrec2(idum,idum,idum,idum,cc,cc,dum,dum,dum,cc,cc,oerr)
         if(oerr)then
            call gsperr('problem in restoring coordinates')
            ofail=.true.
         endif
      endif
 999  return
      end
      subroutine gsperr(mesg)
c
      character*(*) mesg
      common/iofile/iread,iw
c
      write(iw,20)mesg
      return
  20  format(/' error:  ',a/)
c
      end
      subroutine index1
INCLUDE(common/sizes)
c
c generates a list of the grids available on the current dumpfile
c
      implicit REAL (a-h,o-z)
      logical orevis,err
      character*132 filnam
      character*16 zgrid(6)
      character*8 zhead,ztitle
      character*4 yed3,ydd
c
INCLUDE(common/sector)
      common/junk/space(maxorb+26),mgrid,ktype
      common/junkc/zhead(4),ztitle(10)
      common/iofile/ir,iwrit,ispp(18),nav
      common/ddnam/yed3,ydd(16),filnam(16)
      common/setfil/ndd,ncdd
c
c here the name was added for transition dens.
c files getgri.f and fetgri.f do not exist already
c
      data zgrid/'electron density',
     +           'amplitude',
     +           'atom-difference',
     +           'potential',
     +           'mol.-difference',
     +           'transition dens.'
     +/
      mlen=maxorb+26+4/nav
c
      write(iwrit,120)yed3,filnam(ncdd)
 120  format(/1x,'dumpfile ( lfn ',a4,') : ',a80)
      write(iwrit,110)
 110  format(1x,'section',6x,'type',9x,'created',13x,'title'/)
      do 20 i=1,508
         call upack3(apos(i),ipos,iclass,ilen)
         if(ipos.le.0) goto 20
         if(iclass.eq.50)then
            call secfet(i,iclass,iblk,err)
            if(err)goto 20
            call rdedx(space,mlen,iblk,numdu)
            call rdchrs(zhead,14,numdu)
            write(iwrit,100)i,zgrid(ktype),zhead(2),zhead(3),
     +                  (ztitle(j),j=1,4)
 100        format(2x,i3,4x,a16,2x,a8,2x,a7,3x,4a8)
            endif
 20      continue
      write(iwrit,*)
      return
      end
_IF1()      subroutine infree(imxfre,imxbuk)
_IF1()c check on permitted maximum values and initialize free storage
_IF1()      common/maxmum/maxfre,maxres,maxint
_IF1()      logical ibad,contrs,shoshr
_IF1()      maxfre=imxfre
_IF1()      maxres=imxbuk
_IF1()      maxint=31
_IF1()c check on the size of maxfre
_IF1()      if(maxfre.lt.100) go to 1
_IF1()c check on resolution
_IF1()      if(maxres.lt.2.or.maxres.gt.1024) go to 3
_IF1()      return
_IF1()c type out error message and go home
_IF1()    1 call errmsg(1,maxfre)
_IF1()      return
_IF1()    3 call errmsg(3,maxres)
_IF1()      return
_IF1()      end
      subroutine iniplt
      common/ploty/shiftx,shifty,xmill,chsze,iplot
      common/plotz/xmax,xmin,ymax,ymin,xscale,yscale
c
c ----- initialise
c
c xmill is the dimension (in mm) of the plot
      xmill=150.0
c shiftx and shifty determine the position (in mm)
c of the plot origin on the screen
      shiftx=40.0
      shifty=15.0
c
      xscale=1.0
      yscale=1.0
      xmin=0.0
      ymin=0.0
      xmax=xmill
      ymax=xmill
c*      npiby2=0
      call paper(1)
      ichan=0
c  define screen size in millimetres (240mm*180mm).
c     call vwindo(0.0,240.0,0.0,180.0)
      call map(0.0,240.0,0.0,180.0)
      call erase1
c     call scrscl(0,23)
      call filnam('GAMESS_GRID')
      return
      end
      subroutine label(x,y,s,lab)
      common/grdata/igr,ipi,curx,cury,icol,size
c
      character*4 lab
c
c        x=x+s
c        y=y+s
c
      sizex=size
      size=s
      call movea(x,y)
      call wrtchx(lab,4,1,0,0)
      size=sizex
      return
      end
      subroutine mainco(rmesh,pts,trues,plate,unused,
     *xv,yv,q,n,height,nocont,mkr,mkx,mxline,mxopen)
      integer trues,pts,q
      integer clopen
      logical plate,unused
      common/iofile/iread,iwrit
      dimension rmesh(*),height(*),xv(*),yv(*)
      dimension trues(*),pts(*),q(*)
      dimension plate(*),unused(*)
      data tonp,tonz,tonn,toffp,toffz,toffn/
     *100000.0,10.0,1.0,
     *0.0,2.0,1.0/
c   set up frame
c***
      t=n
      call scalep(0.0,0.0,t,t)
c   draw boundary
      call drwbdy(pts)
      do 10 ih=1,nocont
         h = height(ih)
         if(h) 9004,9005,9006
 9006    thron=tonp
         throff=toffp
         goto 9007
 9005    thron=tonz
         throff=toffz
         goto 9007
 9004    thron=tonn
         throff=toffn
 9007    continue
c
c      write(iwrit,9000)ih,h
c 9000 format(////' contour ',i3,1x,f10.4/)
c
c      prepare unused at each height
         call andlog(rmesh,plate,unused,h, n)
c      find open contour ends
         call qrite(pts,rmesh,q,h,n,mxopen)
c      follow contour
   9     call flcont(rmesh,unused,q,trues,h,xv,yv,n,mxline)
c
c     write(iwrit,9001)xv(1),yv(1)
c 9001 format(5x,2f10.5)
c
c      remove multiple points
         call redunp(xv,yv)
c
c      write(iwrit,9001)xv(1),yv(1)
c      line=xv(1)
c      do 9002 i=1,line
c 9002 write(iwrit,9003)i,xv(i),yv(i)
c 9003 format(i5,2f10.3)
c
c      test if all lines already plotted
         if(yv(1).eq.0.0) go to 10
c      draw smooth lines
         call ark(xv,yv,2,ifix(xv(1)),ifix(yv(1)),n,thron,throff)
c
c   annotate lines
c     if(iabs(mkr).lt.2)go to 9
c     call plotht(xv,yv,2,ifix(xv(1)),ih)
c
         goto 9
c        call alfmod
10       continue
c
c*** spot indication
c  if(mkr.gt.0)c1all op1(n)
c     write(iwrit,100)
c 100 format(1h0,' end contour')
c      if(iabs(mkr).lt.2)return
c  call frames(1)
c  call formit('(4x,$contours$,//,4x,$no$,2x,$height$,/,
c    110(4x,i2,f9.3,/))')
c  do 11 i=1,nocont
c11  call outvar(i,height(i))
c  call jctext(dum,dum,0,1,'unpacked.')
c  call resetf
c
      return
      end
_IF1()c
_IF1()c     ****************************************************
_IF1()c     *                                                  *
_IF1()c     * extra routines to convert to <ghost> standard    *
_IF1()c     *                                                  *
_IF1()c     ****************************************************
      subroutine movea(x,y)
c
c****** movea(x,y) ***************************************
c
      integer*4 igr,ipi
      common/grdata/igr,ipi,curx,cury
      call positn(x,y)
      if(igr.ne.0)then
      if(abs(curx-x).gt.1e-6.or.abs(cury-y).gt.1e-6) write(igr,10)x,y
      curx=x
      cury=y
10    format(2f11.4,' m')
      endif
      return
      end
      subroutine mover(dx,dy)
c
c****** mover(dx,dy) *************************************
c
      integer*4 igr,ipi
      common/grdata/igr,ipi,curx,cury
      call move(dx,dy)
      if(igr.ne.0)then
      curx=curx+dx
      cury=cury+dy
      if(abs(dx).gt.1e-6.or.abs(dy).gt.1e-6) write(igr,10)curx,cury
10    format(2f11.4,' m')
      endif
      return
      end
      subroutine qrite(pts,rmesh,q,h,n,mxopen)
c-------------------------------------------------------------
c   note one end of each open contour by hunting round boundary
       integer pts,q
       common/iofile/iread,iwrit
       dimension rmesh(*),pts(*),q(*)
      mb = 1
      idx = n*(pts(4)-1) + pts(3)
      hnew = rmesh(idx)
      max = 2*pts(1) + 1
      do 100 i=3,max,2
      hold =hnew
      idx = n*(pts(i+1)-1) + pts(i)
      hnew = rmesh(idx)
      if(hold.ge.hnew) go to 100
      if(hold.ge.h.and.hnew.gt.h) go to 100
      if(hold.le.h.and.hnew.lt.h) go to 100
      mb = mb + 1
       if(mb.le.mxopen)go to 102
      write(iwrit,5002)h
 5002 format(/' error detected in constructing ',
     +f10.4,' contour'//)
      call fatal
 102  mx = 4*(mb-1)+1
      q(mx) = pts(i-2)
      q(mx+1) = pts(i-1)
      q(mx+2) = pts(i)
      q(mx+3) = pts(i+1)
100   continue
       mx=4*(mb-1)+1
      q(1) = 1
      q(2) = mb
      q(3) = mb-1
       mxpp=mx+3
c     write(6,101) (q(k),k=1,mxpp)
c 101 format(1h0,' ends of open contours',/,4(2x,i6),/,13(10(2x,i6),/))
      return
      end
      subroutine rdrec2(n1,n2,n3,n4,b,c,d1,d2,d3,e,f,oerr)
      implicit REAL (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
      common/blkin/h(2)
      common/iofile/ ir,iw,ipu(18),nav
INCLUDE(common/sector)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      dimension b(*),c(*),e(*),f(*),ih(2)
      equivalence (ih(1),h(1))
      data m15/15/
      mach=12*nat+3+4/nav
      mxcen3 = nat * 3
      maxtot = mxcen3*4+3
      call secfet(isect(493),m15,iblk,oerr)
      if(oerr)return
      call rdedx(h,mach,iblk,numd)
      call dcopy(mxcen3,h,1,b,1)
      call dcopy(mxcen3,h(  mxcen3+1),1,c,1)
      call dcopy(mxcen3,h(2*mxcen3+1),1,e,1)
      call dcopy(mxcen3,h(3*mxcen3+1),1,f,1)
      itemp=4*mxcen3+1
      d1 = h(itemp  )
      d2 = h(itemp+1)
      d3 = h(itemp+2)
      imax = nav*maxtot
      n1 = ih(imax+1)
      n2 = ih(imax+2)
      n3 = ih(imax+3)
      n4 = ih(imax+4)
      return
      end
      subroutine redunp(xv,yv)
c-------------------------------------------------------------
c   remove identical copies of consecutive points
       dimension xv(*),yv(*)
      if(yv(1).eq.0.0) go to 12
      max = ifix(xv(1))
      imax = max
      icur = 3
      xold = xv(2)
      yold = yv(2)
      do 10 i=3,imax
      xnew = xv(i)
      ynew = yv(i)
      if((xold.eq.xnew).and.(yold.eq.ynew)) go to 5
      xold = xnew
      yold = ynew
      xv(icur) = xold
      yv(icur) = yold
      icur = icur + 1
      goto 10
5     max = max - 1
10    continue
c   first and last points compared
      if((xv(2).eq.xv(max)).and.(yv(2).eq.yv(max)))max=max-1
c   do not plot if only one point
      if(max.le.2) yv(1) = 0.0
      xv(1) = float(max)
12    return
      end
      subroutine scalep(xst,yst,xfin,yfin)
      common/ploty/shiftx,shifty,xmill,chsze,iplot
      common/plotz/xmax,xmin,ymax,ymin,xscale,yscale
      common/iofile/iread,iwrit
      goto 50
      entry limits(xst,yst,xfin,yfin,sx,sy)
 50   continue
      xd=xfin-xst
      yd=yfin-yst
      if(xd.lt.0.0 .or. yd.lt.0.0) goto 100
c
c --- convert to millimetres
c
      xscale=xmill/xd
      yscale=xmill/yd
c
c add a correction factor into xmin and ymin so the nuclei and
c contours agree on screen coordinates.
c
      xmin=xst+0.5
      ymin=yst+0.5
      xmax=xfin
      ymax=yfin
      return
c
c --- error message
c
 100  call gsperr('in scale')
      write(iwrit,110)xst,yst,xfin,yfin
 110  format(1x,'xst,yst,xfin,yfin=',4d14.6)
      call fatal
      end
      subroutine secfet(mpos,mtype,iblock,oerr)
c
c -- a modified secget, edited especially for gamsplot
c
      implicit REAL (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sector)
      common/iofile/iread,iwrit
c
      oerr=.false.
      if(mpos.lt.1) then
         call gsperr('attempting to use invalid dumpfile section')
         oerr=.true.
         return
         endif
      call upack3(apos(mpos),ipos,iclass,ilen)
      if(ipos.ge.1)go to 20
      call gsperr('attempt to retreive undefined dumpfile section')
      oerr=.true.
      return
 20   if(mtype.eq.0)mtype=iclass
      if(iclass.ne.mtype) then
         call gsperr('retreived dumpfile section is of wrong type')
         oerr=.true.
         return
         endif
      iblock=ipos+iblkdu
      call search(iblock,num3)
      return
      end
_IF(convex,hp700,hpux11)
_IF(cio)
      subroutine setd(yin)
c
c --- set dumpfile ddname
c
      implicit REAL (a-h,o-z)
      character*4 ydd,yin,yed3
INCLUDE(common/sizes)
      character *132 filnam
      common/ddnam/yed3,ydd(maxlfn),filnam(maxlfn)
      common/setfil/ndd,ncdd
      common/iofile/iread,iwrit
INCLUDE(common/discc)
c
      if(yin.eq.'    ')then
         call gsperr('setd: missing lfn for dumpfile')
         goto 1030
         endif
      do 1020 i=1,ndd
         if(yin.eq.ydd(i))goto 1040
 1020    continue
         ndd = ndd + 1
         ydd(ndd) = yin
         i = ndd
 1040 yed3=yin
      ncdd=i
      write(iwrit,1035)yed3
 1035 format(/'new dumpfile = ',a4/)
c
c initialise dumpfile i/o
c
      ibld1=1
      do 1210  i=1,maxlfn
         numd1=i
         if(yed3.eq.yed(i)) goto 1220
 1210    continue
         ndd = ndd - 1
         call gsperr('setd: ipos(numd1) .le. 0  -- s/r setd')
      return
 1220 continue
      filnam(ndd) = zedfil(numd1)
      call secini(ibld1,numd1)
 1030 return
      end
_ENDIF
_IF(fortio)
      subroutine setd(yin)
c
c --- set dumpfile ddname
c
      implicit REAL (a-h,o-z)
      character*4 ydd,yin,yed3
      character *132 filnam
INCLUDE(common/sizes)
      common/ddnam/yed3,ydd(maxlfn),filnam(maxlfn)
      common/setfil/ndd,ncdd
INCLUDE(common/discc)
      common/iofile/iread,iwrit
c
      if(yin.eq.'    ')then
         call gsperr('setd: missing lfn for dumpfile')
         goto 1030
         endif
      do 1020 i=1,ndd
         if(yin.eq.ydd(i))goto 1040
 1020    continue
         ndd=ndd+1
         ydd(ndd)=yin
         i=ndd
 1040 yed3=yin
      ncdd=i
      write(iwrit,1035)yed3
 1035 format(/'new dumpfile = ',a4/)
c
c initialise dumpfile i/o
c
      ibld1=1
      do 1210  i=1,maxlfn
         numd1=i
         if(yed3.eq.yed(i)) goto 1220
 1210    continue
         ndd=ndd-1
         call gsperr('setd: invalid lfn for dumpfile')
      return
 1220 filnam(ndd)=zedfil(numd1)
      call secini(ibld1,numd1)
 1030 return
      end
_ENDIF
_ENDIF
      subroutine shonuc(oexp,nuccol,size,sizel2,pltlab,pltcro)
INCLUDE(common/sizes)
c
c --- routine to show the positions of the nuclei projected onto
c     the grid plane.
c
c  don't forget to declare `ddot' as double prec!!!
      REAL `ddot',popocc,plane,c,czan
      REAL d(3),e(3),f(3),g(3),orig(3),x(3)
      REAL ex,fx,gx,hside
      logical*4 pltlab,pltcro
      logical oexp
      character*4 atlab4
      character *8 zcom,atlab
      common/infoa/nat,ich(7),czan(maxat),c(3,maxat)
      common/junk/popocc(maxorb),plane(26)
      common/ploty/shiftx,shifty,xmill
      common/junkc/zcom(24),atlab(maxat)
      common/iofile/iread,iw
      common/grdata/igr,ipi,curx,cury,icol,sz
c
      scale=xmill/sngl(plane(26))
      hside=plane(26)/2.0d0
c change colour
      call colour(nuccol)
c
c set up axes of plane
c
c points defining plane:
c    a is plane(1 to 3)
c    b is plane(4 to 6)
c    c is plane(7 to 9)
      do 10 i=1,3
c  d = c - a
         d(i)=plane(6+i)-plane(i)
c  e = b - a
         e(i)=plane(3+i)-plane(i)
c  a is the origin of the plane axes
         orig(i)=plane(i)
  10     continue
c calculate orthogonal set of axis
      call crprod(d,e,f)
c g defines horizontal axis of plot
      call crprod(e,f,g)
c make set of axes e,f,g unit vectors
c  e,g are in plane, f is perpendicular
      call univec(e)
      call univec(f)
      call univec(g)
c
c loop for each nuclear coordinate
c
      write(iw,*)' -------------------------------------'
      write(iw,*)' positions of nuclei from the dumpfile'
      write(iw,*)' -------------------------------------'
      do 20 i=1,nat
         do 30 j=1,3
c  x is current nucleus coord. translated to new origin
  30        x(j)=c(j,i)-orig(j)
         write(iw,22) i,(c(j,i),j=1,3)
22       format(i4,2x,3f14.5)
c find components of x in new coord. system
         ex=ddot(3,e,1,x,1)
         fx=ddot(3,f,1,x,1)
         gx=ddot(3,g,1,x,1)
c    plot position of nucleus:
c translate since origin is centre of plot
         gx=gx+hside
         ex=ex+hside
c convert coords to screen mm and shift onto screen origin
         xp=sngl(gx)*scale+shiftx
         yp=sngl(ex)*scale+shifty
c if nucleus is in plane, plot a star, if it has been projected
c onto the plane, plot a cross.
         if(pltcro) then
         if(fx.lt.1.0d-5.and.fx.gt.-1.0d-5) then
            call sstar(xp,yp,size)
          else
            call cross(xp,yp,size)
          endif
          endif
          if(pltlab) then
          atlab4=atlab(i)(1:4)
          icolsave=icol
          icol=nuccol
          call label(xp,yp,sizel2,atlab4)
          icol=icolsave
          endif
   20     continue
      write(iw,*)' -------------------------------------'
c write key to nuclei symbols
      if(oexp) then
      xp=shiftx+xmill+3.0
      yp=shifty+xmill-20.0
      call movea(xp,yp)
      call wrtchr('nuclei:  ',9,1,0,0)
      xp=xp+1.5
      yp=yp-4.0
      call sstar(xp,yp,3.0)
      call mover(4.0,0.0)
      call wrtchr('in plane ',9,1,0,0)
      yp=yp-5.0
      call cross(xp,yp,3.0)
      call mover(4.0,0.0)
      call wrtchr('projected',9,1,0,0)
c change colour back
      endif
      call colour(1)
      return
      end
      subroutine sstar(x,y,sz)
c
c  x,y are coordinates in mm for the star position
c  s   is the size, in mm, of the star
c
      logical typ
      common/startyp/typ
      s=sz*0.5
      d=s*0.8
      if(typ)then
      call movea(x-s,y)
      call drawa(x+s,y)
      endif
      call movea(x-d,y+d)
      call drawa(x+d,y-d)
      call movea(x+d,y+d)
      call drawa(x-d,y-d)
      if(typ)then
      call movea(x,y+s)
      call drawa(x,y-s)
      endif
c     call alfmod
      return
      end
_IF1()c
_IF1()c****** finitt(i,j) ***************************************
_IF1()c
_IF1()c     subroutine finitt(i,j)
      subroutine stoplt
      call grend
      return
      end
      subroutine tmplt(trues,pts,plate,kount,n)
c-------------------------------------------------------------
      integer trues,pts
      logical plate
      logical flip,l
      common/junk2/ilifg(200)
      dimension trues(*),pts(*)
      dimension plate(*)
      kount=0
      ng2=n*n
      n1 = n-1
      n2 = n1-1
c   set up as false at start of each line
      do 100 j=1,n2
      l=.false.
      flip = .false.
      do 100 i=1,n1
      ix = n*(j) +i
      if(trues(ix).gt.0) go to 50
      if(flip) go to 95
      goto 100
c   point on boundary
   50 if(l) go to 70
      idx = trues(ix)
      if(idx.gt.ng2) idx = idx - ng2
      idx = idx*2+2
c   previous boundary point below,alongside,above present one
      if (pts(idx-2) - pts(idx)) 55,56,57
55    ivar = 2
      idif = 1
      goto 70
57    ivar = 2
      idif =-1
      goto 70
56    idx = trues(ix)
      if (idx.gt.ng2) idx = 1
      idx = idx*2+2
c   next boundary point below,alongside,above present one
      if(pts(idx+2)-pts(idx)) 61,100,63
61    ivar =-2
      idif = 1
      goto 70
63    ivar =-2
      idif =-1
70    l = .true.
c   l is true while boundary crossing is incomplete
      idx = trues(ix)
      if(ivar) 71,73,72
71    if (idx.gt.ng2) idx = idx - ng2
      goto 73
72    if (idx.gt.ng2) idx=1
73    idx = idx*2+2
      k =idx + ivar
      if((pts(k)-pts(idx)).eq.idif) l=.false.
      if(l) go to 78
      flip = .not. flip
   78 if(flip) go to 95
      goto 100
   95 if(trues(ix+1).gt.0) go to 100
      kount = kount + 1
      trues(ix+1) = 0
100   plate(j+ilifg(i))  = flip
      return
      end
      subroutine toutst(j,ichar)
c
c****** toutst(j,ichar) **********************************
c
      return
      end
      subroutine trigs(xv,yv,m,n,i,ipen,itype,sina,cosa)
c-------------------------------------------------------------
      dimension xv(*),yv(*)
      iuse=i
      if (i.lt.m.or.i.gt.n) goto 20
      goto  (14,4,16),itype
14    isave = i
      if (ipen.eq.2) goto 5
      j= i+1
      if(j.gt.n) j=m
      dxold =xv(j) - xv(i)
      dyold = yv(j) - yv(i)
      i = j
      j=  i+1
      if (j.gt.n) j=m
      dxnew = xv(j) - xv(i)
      dynew = yv(j) - yv(i)
      dolds = dxold*dxold + dyold*dyold
      dnews = dxnew*dxnew + dynew*dynew
      alpha = dnews + 2*sqrt(dnews*dolds)
      beta  = -dolds
1     c         = dxold *alpha + dxnew*beta
      s         = dyold*alpha  + dynew*beta
      an    = sqrt(c*c + s*s)
      sina  = s/an
      cosa  = c/an
      i=iuse
      return
5     isave = i
      j         = i + 1
      if (j.gt.n) j=m
      dxnew = xv(j) -xv(i)
      dynew = yv(j) - yv(i)
      j         = i
      i         = j-1
      if (i.lt.m) i=n
      dxold = xv(j) - xv(i)
      dyold = yv(j) - yv(i)
      dolds = dxold*dxold + dyold*dyold
2     dnews = dxnew*dxnew + dynew*dynew
      c         = dxold*dnews + dxnew*dolds
      s         = dyold*dnews + dynew*dolds
      an    = sqrt(c*c + s*s)
      sina  = s/an
      cosa  = c/an
      i=iuse
      return
4     if (i-isave-1) 5,6,5
6     isave = i
      dxold = dxnew
      dyold = dynew
      dolds = dnews
      j         = i + 1
      if (j.gt.n) j=m
      dxnew = xv(j) - xv(i)
      dynew = yv(j) - yv(i)
      goto 2
16    if (ipen.eq.1) goto 4
      if (i-isave-1) 9,10,9
10    dxold = dxnew
      dyold = dynew
      dolds = dnews
      j         = i
      i         = i - 1
      if ( i.lt.m) i=n
      dxnew = xv(j) - xv(i)
      dynew = yv(j) - yv(i)
      goto 11
9     j         = i
      i         = i - 1
      if (i.lt.m) i=n
      dxnew = xv(j) - xv(i)
      dynew = yv(j) - yv(i)
      j         = i
      i         = i - 1
      if ( i.lt.m) i=n
      dxold = xv(j) - xv(i)
      dyold = yv(j) - yv(i)
      dolds = dxold*dxold + dyold*dyold
11    dnews = dxnew*dxnew + dynew*dynew
      alpha = -dnews
      beta  = dolds + 2*sqrt(dolds*dnews)
      isave = -1
      goto 1
20    sina  = 10
      cosa  =10
      i=iuse
      return
      end
      subroutine univec(a)
c
c converts vector a into the corresponding unit vector
c
      REAL a(3),t,l
c
      t=0.0d0
      do 10 i=1,3
  10     t=t+a(i)*a(i)
      l=dsqrt(t)
      do 20 i=1,3
  20     a(i)=a(i)/l
      return
      end
      subroutine wrtchx(istr,idim,jdim,x,y)
c     special version for atom labels
      integer*4 igr,ipi,icol
      character*(*) istr(jdim)
      character*1 z
      dimension z(256)
      common/grdata/igr,ipi,curx,cury,icol,size
      do 20 j=1,jdim
20        call typecs(istr(j))
      if(igr.ne.0)then
      jmax=idim*jdim
      jndex=-idim
      do 30,j=0,jdim-1
      jndex=jndex+idim
      do 30,i=1,idim
      z(i+jndex)=istr(j+1)(i:i)
30    continue
      j=jmax
40    continue
      if(j.gt.1.and.z(j).le.' ')then
        j=j-1
        goto 40
      endif
      jmax=j
      write(igr,10)curx,cury,icol,size,(z(j),j=1,jmax)
10    format(2f11.4,' p ',i2,' 9 ',f10.7,' 0 0 t',256a1)
      endif
      return
      end
      subroutine wrtchr(istr,idim,jdim,x,y)
      integer*4 igr,ipi,icol
      character*(*) istr(jdim)
      character*1 z
      dimension z(256)
      common/grdata/igr,ipi,curx,cury,icol,size
      do 20 j=1,jdim
20        call typecs(istr(j))
      if(igr.ne.0)then
      jmax=idim*jdim
      jndex=-idim
      do 30,j=0,jdim-1
      jndex=jndex+idim
      do 30,i=1,idim
      z(i+jndex)=istr(j+1)(i:i)
30    continue
      j=jmax
40    continue
      if(j.gt.1.and.z(j).le.' ')then
        j=j-1
        goto 40
      endif
      jmax=j
      if(cury.lt.100.0)then
      write(igr,10)curx,cury,icol,size,(z(j),j=1,jmax)
      else
      write(igr,50)curx,cury,icol,size,(z(j),j=1,jmax)
      endif
10    format(2f11.4,' p ',i2,' 1 ',f10.7,' 0 0 t',256a1)
50    format(2f11.4,' p ',i2,' 2 ',f10.7,' 0 0 t',256a1)
      endif
      return
      end
      if(j.gt.1.and.z(j).le.' ')then
        j=j-1
        goto 40
      endif
      jmax=j
      if(cury.lt.100.0)then
      write(igr,10)curx,cury,icol,size,(z(j),j=1,jmax)
      else
      write(igr,50)curx,cury,icol,size,(z(j),j=1,jmax)
      endif
10    format(2f11.4,' p ',i2,' 1 ',f10.7,' 0 0 t',256a1)
50    format(2f11.4,' p ',i2,' 2 ',f10.7,' 0 0 t',256a1)
      endif
      return
      end

      subroutine ver_plot(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/plot.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
