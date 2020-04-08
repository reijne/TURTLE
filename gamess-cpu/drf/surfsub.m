      subroutine drfsurf(ieps,xscm)
c-----------------------------------------------------------------------
c   this is the driver subroutine for the construction of an enveloping
c   molecular surface.
c
c   the structure of this routine is as follows:
c
c   - calculation of raw surface points if ibem = 1
c
c   - construction of the initial surface sphere with radius 1
c     and 60, 240 or 960 triangles (32, 122 or 482 vertices)
c
c   - construction of enveloping molecular surface from:      ibem:
c     * raw surface points                                      1
c     * connolly surface points (to be read in from unit 42)    2
c     * quantum mechanical and classical atom positions         3
c     * sphere of radius spherad, with origin in (0,0,0)        4
c     * connolly surface points, areas and normal vectors used
c       directly from unit42                                    5
c     * n faces of cube circumscribing sphere of radius sherad -n
c
c   - calculation of midpoints, normal vectors and areas of the
c     triangles that make up the surface
c
c   - checks on surface
c
c   - (debug) printing
c
c     the midpoint vectors (rm), normal components at midpoints (nm)
c     and triangle areas (area) are written to the drf da-file with
c     labels 31, 32 and 33, respectively
c
c-----------------------------------------------------------------------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
INCLUDE(comdrf/sizescon)
INCLUDE(../m4/common/connolly)
c
c------  maxinv and maxcon needed for declarations
c
cmw      parameter (maxinv=10,maxnrc=10000)
c
c------  work array
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/max)
INCLUDE(comdrf/mem)
c
c nu als argument
c
c      include 'gen/scm'
c
INCLUDE(comdrf/dafil)
c
c     common/drf_in/dstgrp,igroup,iunit
INCLUDE(comdrf/drfin)

INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/extinf)
INCLUDE(comdrf/drfdaf)
c
INCLUDE(comdrf/rad)
INCLUDE(comdrf/drfbem)
c
INCLUDE(comdrf/runpar)
INCLUDE(comdrf/mcinf)
c
      parameter (nrsl=2000)
c
      REAL xscm
      dimension xscm(*)
c
c------  hondo an  drf declarations
c
      dimension p(3)
c
c------  local variables
c
      dimension isurf1(5,12)
      integer isurf1
c
      integer maxnrs,maxntr,nrs,ntr
c
      integer ixrs,ixsurf,ixinv,ixrm,ixnm,ixarea,ixnorm,ixedge,
     1 ixraw,ixcon,ixouts,ixshell,ixclose,ixdclos,ixrad2,last
c
      integer nrspts,nrc
c
      REAL maxedge,minedge,pdens,scale,delta,sdis
      REAL pi,one
c
      character*8 forw,back
c
      logical lcni,lpla,lmid,luni,lmidout
c
      integer itel
      dimension cgeom(3)
c
c-----  dummy variables
c
      dimension xsurf(3,6),xnorm(3,6),area(6)
c
      data zero,one,two,four /0.0d00,1.0d00,2.0d00,4.0d00/
      data forw,back /'forward ','backward'/
c      data scale /3.0d00/
      data scale /1.0d00/
c
      data bohr /0.529177249d00/
c
c-----  begin
c
      write(iwr,11)
   11 format(/1x,104('-')/)
c
      pi = four*atan(one)
      epsilon = eps1
      isur = isur1
      ilxs = 81
      ilxn = 82
      ilar = 83
      if ((ieps .eq. 1) .and. (itwosur .eq. 1)) then
        epsilon = eps2
        isur = isur2
        open (unit=isur2, form='formatted', file='csurfop',
     +        status='unknown')
        ilxs = 84
        ilxn = 85
        ilar = 86
      endif
      epsfact = one/(two*pi*(one+epsilon))
c
      if (iunit .eq. 1) then
        if (insphr .eq. 1) spherad = spherad/bohr
        if (insp .eq. 1) d = d*bohr*bohr
        if (inpro .eq. 1)  rprobe = rprobe/bohr
        if (inrp .eq. 1)  rp = rp/bohr
        if (inprj .eq. 1)  rprobej = rprobej/bohr
        if (insw .eq. 1) swidth = swidth/bohr
        if (insd .eq. 1) sdist = sdist/bohr
        if (insr .eq. 1) solrad = solrad/bohr
      endif
c
      if (ibem .gt. 0) goto 1
c
c-----  define 6 faces of cube circumscribing sphere of radius spherad
c
      r = spherad
      nbem = -ibem
c
      xsurf(1,1) = r
      xsurf(2,1) = 0.0d0
      xsurf(3,1) = 0.0d0
      xnorm(1,1) = 1.0d0
      xnorm(2,1) = 0.0d0
      xnorm(3,1) = 0.0d0
c
      xsurf(1,2) = -r
      xsurf(2,2) = 0.0d0
      xsurf(3,2) = 0.0d0
      xnorm(1,2) = -1.0d0
      xnorm(2,2) = 0.0d0
      xnorm(3,2) = 0.0d0
c
      xsurf(1,3) = 0.0d0
      xsurf(2,3) = r
      xsurf(3,3) = 0.0d0
      xnorm(1,3) = 0.0d0
      xnorm(2,3) = 1.0d0
      xnorm(3,3) = 0.0d0
c
      xsurf(1,4) = 0.0d0
      xsurf(2,4) = -r
      xsurf(3,4) = 0.0d0
      xnorm(1,4) = 0.0d0
      xnorm(2,4) = -1.0d0
      xnorm(3,4) = 0.0d0
c
      xsurf(1,5) = 0.0d0
      xsurf(2,5) = 0.0d0
      xsurf(3,5) = r
      xnorm(1,5) = 0.0d0
      xnorm(2,5) = 0.0d0
      xnorm(3,5) = 1.0d0
c
      xsurf(1,6) = 0.0d0
      xsurf(2,6) = 0.0d0
      xsurf(3,6) = -r
      xnorm(1,6) = 0.0d0
      xnorm(2,6) = 0.0d0
      xnorm(3,6) = -1.0d0
c
      area(1) = 0.666666666666667d0*pi*r**2*epsfact
      area(2) = 0.666666666666667d0*pi*r**2*epsfact
      area(3) = 0.666666666666667d0*pi*r**2*epsfact
      area(4) = 0.666666666666667d0*pi*r**2*epsfact
      area(5) = 0.666666666666667d0*pi*r**2*epsfact
      area(6) = 0.666666666666667d0*pi*r**2*epsfact
c
      ntr = 6
c
      call dawrit(idafdrf,iodadrf,xsurf,3*ntr,ilxs,navdrf)
      call dawrit(idafdrf,iodadrf,xnorm,3*ntr,ilxn,navdrf)
      call dawrit(idafdrf,iodadrf,area,   ntr,ilar,navdrf)
c
      goto 300
c
    1 continue
c
c     call setscm(i10)
c     call cmem(loadcm)
c     loc10 = loccm(xscm(i10))
      lwor = igmem_max_memory()
      i10 = igmem_alloc(lwor)
c
cxxx  if ((ibem .eq. 2) .or. (ibem .eq. 5)) then
c 1-----
c  -----  read connolly's vdwaals surface
c         first count # boundary elements on connolly's surface
c
c        rewind (isur)
c        ncon = 0
c 1000   format(e12.6)
c 123    read(isur,1000,end=125)  dum
c        ncon = ncon + 1
c        goto 123
c 125    continue
c        rewind (isur)
        if (ibem .eq. 5) then
c   2-----
          write(iwr,1005)
 1005     format(/,' Connolly surface construction')
          nbem = 1
c         ixrm = i10
c         ixnm = ixrm + 3*nbem
c         ixarea = ixnm + 3*nbem
c         ixlast = ixarea + nbem + 1
c         need = loc10 + ixlast - i10
c     need = ixlast - i10
c
          natm = natom
          ixrtp = i10
          ixrad = ixrtp + natom + 1
          ixua = ixrad + natom + 1
          ixeva = ixua + 3*maxsph*natom + 1
          ixpy = ixeva + 3*maxsph*natom + 1
          ixay = ixpy + 3*maxyon + 1
          ixnua = ixay + 3*maxyon + 1
          ixity = ixnua + natom + 1 
          ixias = ixity + natom + 2
          ixiat = ixias + natom + 2
          ixico = ixiat + natom + 2
          ixicu = ixico + 3*maxyon + 1
          ixmol = ixicu + maxyon + 1
          ixsrs = ixmol + maxyon + 1
          last1 = ixsrs + natom + 1
c
          ixrm = last1
          ixnm = ixrm + 3*nbem
          ixarea = ixnm + 3*nbem
          ixlast = ixarea + nbem + 1
          need = ixlast - i10
c
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for drfsurf'
	call caserr('Insufficient memory for drfsurf')
      endif
c         call setc(need)
          call mscon(xscm(ixrm),xscm(ixnm),xscm(ixarea),
     1    xscm(ixrtp),xscm(ixrad),xscm(ixua),xscm(ixeva),
     2    xscm(ixpy),xscm(ixay),xscm(ixnua),xscm(ixity),
     3    xscm(ixias),xscm(ixiat),xscm(ixico),xscm(ixicu),
     4    xscm(ixmol),xscm(ixsrs),
     5    natm,nbem,ncon,1)
c
c   find out how many surface elements will be generated
c

          ntr = ncon
          nbem = ncon
          ixrm = last1
          ixnm = ixrm + 3*nbem
          ixarea = ixnm + 3*nbem
          ixlast = ixarea + nbem + 1
c         need = loc10 + ixlast - i10
      need = ixlast - i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for drfsurf'
	call caserr('Insufficient memory for drfsurf')
      endif
c         call setc(need)
c
c  Define surface elements
c
          call mscon(xscm(ixrm),xscm(ixnm),xscm(ixarea),
     1    xscm(ixrtp),xscm(ixrad),xscm(ixua),xscm(ixeva),
     2    xscm(ixpy),xscm(ixay),xscm(ixnua),xscm(ixity),
     3    xscm(ixias),xscm(ixiat),xscm(ixico),xscm(ixicu),
     4    xscm(ixmol),xscm(ixsrs),
     1    natm,nbem,ncon,2)
          if (ibemout .lt. 2) then
            write(iwr,1010) nbem
 1010       format(/,' Generated ',i4, ' surface elements')
          endif
          goto 200
c   2-----
        endif
c        close (unit=isur,status='keep')
c 1-----
cxxx  endif
c
      delta = zero
      sdis = zero
      nav = lenwrd()
c
      lcni = .false.
      lpla = .false.
      lmid = .false.
      luni = .false.
      lmidout = .false.
c
c
c-----  enveloping surface constructed from raw surface points
c       if ibem = 1: set up pointers for xscm and call appropriate
c       subroutine bemraws
c
      if (ibem .eq. 1) then
c 1-----
c  -----  set maximum dimensions for local sphere
c         maxnrs = # of surface points; maxntr = # of surface triangles
c
c  -----  recursion formula, starting at option leveli=0, maxnrs=32
c         each following level: maxnrs(i+1) = 4*maxnrs(i) - 6
c         allowed values: leveli=0,1,2 (checked in rfin)
c
        maxnrs = 32
        do 100, i = 1, leveli
          maxnrs = 4 * maxnrs - 6
  100   continue
        maxntr = 2 * maxnrs - 4
c
c  -----  calculate pointers:
c          -ixrs    unit sphere coordinate array
c          -ixsurf  surface location array
c          -ixrad2  atom radii**2 array
c
c
        ixrs = i10
        ixsurf = ixrs + 3 * maxnrs
        ixrad2 = ixsurf + 3 * maxntr
        ixnew  = ixrad2 + natmax + maxpol
        last   = ixnew + (nrsl*nrsl)/nav + 1
c
c  -----  check if xscm is sufficiently large
c
      need = last - i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for drfsurf'
	call caserr('Insufficient memory for drfsurf')
      endif
c       call setc(need)
c
        call bemraws(xscm(ixrs),xscm(ixsurf),
     2               maxnrs,maxntr,xscm(ixrad2),
     3               xscm(ixnew),nrsl)
c
c  -----  calculate an approximate point-density for the raw surface:
c         because the point-density is not homogeneous, but depends loca
c         on the radius of the atoms, a minimum point-density is obtaine
c         from pdens=number of surface_points / 4*pi*(radmax**2)
c         (4pir**2 is the area of a sphere with radius r)
c         in which radmax is the biggest atom radius, determined in
c         subroutine rfin
c
        pdens = maxnrs / (4 * pi * (radmax**2))
c 1-----
      endif
c
c-----  construct initial unit sphere with number of
c       triangles determined by nbem (via levelo in rfin)
c            nbem = 4**levelo * 60
c
c-----  set dimensions
c
      maxntr = nbem
      maxnrs = (nbem + 4) / 2
c
c-----  set up pointers for construction of unit sphere
c
c       all previously used work-space is overwritten
c
c-----  ixrs   : begin of rs(3,maxnrs) surface points array
c       ixsurf : begin of surf(3,maxntr) array (integer!!)
c       ixinv  : begin of involv(maxinv,maxnrs) array (integer)
c       ixrm   : begin of triangle midpoint rm(3,maxntr) array
c       ixnm   : begin of triangle midpoint normal component
c                tnorm(3,maxntr) array
c       ixarea : begin of triangle area area(maxntr) array
c
c-----  these 6 arrays are the absolute necessity, for storage of
c       initial (final) surface point coordinates (rs), connection
c       between surface points (surf) ,list of triangles in which
c       surface points are involved (involv), triangle midpoints (rm),
c       normal components at triangle midpoints (tnorm) and triangle
c       areas (area).
c
c-----  only rs and surf are permanently needed
c
c-----  remaining pointers:
c
c        ixnorm : begin of vertex normal vector component array
c                 norms(3,maxnrs)
c        ixedge : begin of triangle edge length array edge(3,maxntr)
c
      ixrs = i10
      ixsurf = ixrs + 3 * maxnrs
c
      ixinv = ixsurf + 3 * maxntr
      ixrm = ixinv + maxinv
      ixnm = ixrm + 3
      ixarea = ixnm + 3
      ixnorm = ixarea + 1
      ixedge = ixnorm + 3
      ixnew  = ixedge + 3
c
c-----  check if xscm is sufficiently large
c
      last   = ixnew + (nrsl*nrsl)/nav + 1
c     need = loc10 + last - i10
      need = last - i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for drfsurf'
	call caserr('Insufficient memory for drfsurf')
      endif
c     call setc(need)
c
c-----  calculate surface points on unit sphere
c
      call setups(one,maxnrs,xscm(ixrs))
      call initsu(isurf1)
      call level1(one,isurf1,maxnrs,nrs,xscm(ixrs),
     1            maxntr,ntr,xscm(ixsurf))
c
      if (levelo .gt. 0)
     1 call levelh(one,levelo,maxnrs,nrs,xscm(ixrs),maxntr,
     2 ntr,xscm(ixsurf),xscm(ixnew),nrsl)
c
c-----  the unit sphere has been constructed.
c
c-----  set flag for property calculation of unit sphere
c
      luni = .true.
c
c-----  if extra printing is required, some extra pointers need
c       to be calculated, and some extra space reserved.
c
      if (iuniout .ge. 10) then
c 1-----
c  -----  vertex coordinates, normal vectors and involvements
c         of the unit sphere are stored and printed on output
c
        lcni = .true.
c
c  -----  adjust remaining pointer(s)
c
c  -----  for involv(maxinv,maxnrs)
c
        ixrm = ixinv + maxinv * maxnrs
c
        ixnm = ixrm + 3
        ixarea = ixnm + 3
        ixnorm = ixarea + 1
c
c  -----  for norms(3,maxnrs)
c
        ixedge = ixnorm + 3 * maxnrs
        last = ixedge + 3 + 1
c
        if (iuniout .ge. 100) then
c   2-----
c    -----  triangle pointers, length of edges and areas are stored
c           and printed.
c
          lpla = .true.
c
c    -----  adjust pointer(s)
c
c    -----  for area(maxntr)
c
          ixnorm = ixarea + maxntr
          ixedge = ixnorm + 3 * maxnrs
c
c    -----  for edge(3,maxntr)
c
          last = ixedge + 3 * maxntr + 1
c   2-----
        endif
c
c       need = loc10 + last - i10
      need = last - i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for drfsurf'
	call caserr('Insufficient memory for drfsurf')
      endif
c       call setc(need)
c
        call compgen(lcni,lpla,lmid,luni,minedge,maxedge,maxinv,nrs,ntr,
     1   xscm(ixrs),xscm(ixsurf),xscm(ixinv),xscm(ixrm),
     2   xscm(ixnm),xscm(ixarea),tarea,xscm(ixnorm),xscm(ixedge))
c
c  -----  write results to output and to file according
c         to value of iuniout
c         if iuniout = 0 writing and computation is skipped
c
        call writgen(lcni,lpla,lmid,luni,lmidout,nrs,ntr,minedge,
     1   maxedge,maxinv,xscm(ixrs),xscm(ixsurf),xscm(ixinv),xscm(ixrm),
     2   xscm(ixnm),xscm(ixarea),tarea,xscm(ixnorm),xscm(ixedge))
c 1-----
      endif
c
c-----  calculate center of nuclear charge of complete system
c
      ztotal = zero
      call clear(cgrav,3)
      iupd = 1
      do 20, n = 1, nxtpts
c 1-----
        if (mcupdt .and.
     2     ((n .ge. istrt) .and. (n .lt. istrt+nnpts))) then
c   2-----
          p(1) = xnpts(1,iupd)
          p(2) = xnpts(2,iupd)
          p(3) = xnpts(3,iupd)
          iupd = iupd + 1
c   2-----
        else
c   2-----
          p(1) = xpts(1,n)
          p(2) = xpts(2,n)
          p(3) = xpts(3,n)
c   2-----
        endif
c
        if ((nxcent(n) .ne. 'qq') .and.
     1      (nxcent(n)(:2) .ne. 'gr')) then
c   2-----
          call drfnval(nxcent(n),nval,znuc)
          ztotal = ztotal + znuc
          cgrav(1) = cgrav(1) + znuc*p(1)
          cgrav(2) = cgrav(2) + znuc*p(2)
          cgrav(3) = cgrav(3) + znuc*p(3)
c   2-----
        endif
c 1-----
   20 continue
c
      do 30, n = 1, nat
c 1-----
        ztotal = ztotal + czan(n)
        cgrav(1) = cgrav(1) + czan(n)*c(1,n)
        cgrav(2) = cgrav(2) + czan(n)*c(2,n)
        cgrav(3) = cgrav(3) + czan(n)*c(3,n)
c 1-----
   30 continue
c
      if (ztotal .gt. 1.0d-5) then
      do 40, n = 1, 3
        cgrav(n) = cgrav(n)/ztotal
   40 continue
       endif
c
      call hatout(cgrav,3,1,2,'cgrav')
c
c-----  reset flags
c
      lcni = .false.
      lpla = .false.
      luni = .false.
c
c-----  construction of the enveloping surface:
c       if from raw surface points, or connolly points, these
c       are to be read in, after the appropriate space has been
c       reserved.
c
      if ((ibem .eq. 1) .or. (ibem .eq. 2)) then
c 1-----
        if (ibem .eq. 1) then
c   2-----
c    -----  raw surface points define molecular surface
c
          if (ibemout .gt. 0) write(iwr,101)
  101     format(// 1x,'the raw surface points determine the ',
     1           'enveloping surface')
c
c    -----  introduce pointer ixraw : marks beginnig of raw surface-poin
c           array - everything from involv on is overwritten
c
          ixraw = ixsurf + 3 * maxntr
          last = ixraw + 3 * nraw + 1
c
c         need = loc10 + last - i10
      need = last - i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for drfsurf'
	call caserr('Insufficient memory for drfsurf')
      endif
c         call setc(need)
c
c    -----  the raw surface points are read in from unit32, which
c           is an unformatted file
c
          open (unit=32,form='unformatted',file='rawsurf',
     +        status='unknown')
          rewind (32)
c
c     -----  initialize help pointer ixrawt
c
          ixrawt = ixraw
c
c    -----  calculate the number of fully filled records on unit32
c
          nrec = nraw / 170
c
c    -----  calculate the significant number of rest-points on record
c           nrec + 1
c
          nrest = nraw - nrec * 170
c
c    -----  loop over filled records
c
          do 10, i = 1, nrec
            read(32) (xscm(j), j = ixrawt, ixrawt + 509)
            ixrawt = ixrawt + 510
  10      continue
c
c    -----  read in rest-points off record nrec + 1
c
          read(32) (xscm(j), j = ixrawt, ixrawt + nrest * 3-1)
c
          close (unit=32,status='keep')
c
          nrspts = nraw
c   2-----
        else
c   2-----
c    -----  connolly points define molecular surface
c
          if (ibemout .gt. 0) write(iwr,102)
  102     format(// 1x,'the connolly points determine the ',
     1           'enveloping surface')
c
c    -----  introduce pointer ixcon : mark beginning of connolly point
c           array. everything from involv on will be overwritten
c
          ixcon = ixsurf + 3 * maxntr
          last = ixcon + 3 * ncon + 1
c
c         need = loc10 + last - i10
      need = last - i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for drfsurf'
	call caserr('Insufficient memory for drfsurf')
      endif
c         call setc(need)
c
          call conread(isur,ibem,pdens,nrc,xscm(ixcon),dum1,dum2)
c
          nrspts = nrc
c
c    -----  calculate geometric centre  of complete system
c
          itel = 0
          call clear(cgeom,3)
          iupd = 1
          do 220, n = 1, nxtpts
c     3-----
            if (nxcent(n) .ne. 'qq') then
c       4-----
              itel = itel + 1
              if (mcupdt .and.
     2     ((n .ge. istrt) .and. (n .lt. istrt+nnpts))) then
c         5-----
                p(1) = xnpts(1,iupd)
                p(2) = xnpts(2,iupd)
                p(3) = xnpts(3,iupd)
                iupd = iupd + 1
c         5-----
              else
c         5-----
                p(1) = xpts(1,n)
                p(2) = xpts(2,n)
                p(3) = xpts(3,n)
c         5-----
              endif
c
              cgeom(1) = cgeom(1) + p(1)
              cgeom(2) = cgeom(2) + p(2)
              cgeom(3) = cgeom(3) + p(3)
c       4-----
            endif
c     3-----
  220     continue
c
          do 230, n = 1, nat
c     3-----
            itel = itel + 1
            cgeom(1) = cgeom(1) + c(1,n)
            cgeom(2) = cgeom(2) + c(2,n)
            cgeom(3) = cgeom(3) + c(3,n)
c     3-----
  230     continue
c
          do 240, n = 1, 3
            cgeom(n) = cgeom(n)/itel
  240     continue
c
          call hatout(cgeom,3,1,2,'cgeom')
c
c    -----  the connolly points are transformed backwards with the
c           centre of geometry to reflect again the input geometry
c
          call cortran(back,cgeom,nrspts,xscm(ixcon))
c   2-----
        endif
c
c  -----  calculate position vectors of the enveloping surface points
c         from molecular surface points spts
c
c  -----  adjust pointers
c
c  -----  ixspts : begin of surface (raw or connolly) points array
c                  rc(3,nrspts)
c         ixused : begin of integer array that determines whether a surf
c                  point has been allocated to a spoke (iused(nrspts))
c
        ixspts = ixsurf + 3 * maxntr
        ixused = ixspts + 3 * nrspts
c
        last = ixused + nrspts + 1
c
c       need = loc10 + last - i10
      need = last - i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for drfsurf'
	call caserr('Insufficient memory for drfsurf')
      endif
c       call setc(need)
c
c  -----  transform surface-point coordinates to have centre of gravity
c         origin (because initial sphere has been defined w.r.t. (0,0,0)
c
        call cortran(forw,cgrav,nrspts,xscm(ixspts))
c
c  -----  calculate the cylinder radius around a spoke from the surface
c         point density
c  -----  the average distance between two raw or connolly points is
c         proportional to 1/sqrt(point-density), easily verified
c         considering a uniform grid on a two-dimensional surface
c         in order to find a surface point in the neighbourhood of
c         a spoke, a safety margin is taken into account for the
c         cylinder radius delta.
c         presently the scale factor scale is set to 3.0, which works
c         well.
c         too large a delta can create artifacts on the enveloping
c         surface, especially if the molecule inside is not
c         more or less spherical.
c
        if (swidth .ge. 0.0d0) then
          delta = swidth
        else
          delta = scale / sqrt(pdens)
        endif
c
        call spokec(delta,nrs,xscm(ixrs),nrspts,
     1    xscm(ixused),xscm(ixspts))
c
c  -----  the molecular surface has been generated from the surface poin
c 1-----
      else
c 1-----
        if (ibem .eq. 3) then
c   2-----
c    -----  the atompositions and radii determine the molecular surface
c
          write(iwr,103)
  103     format(// 1x,'the atomic (qm and/or classical) positions ',
     1   'determine the surface through Juffer''s procedure') 
c
c    -----  the cylinder radius (delta) and the minimum atom-spoke dista
c           (sdis) are set
c
          if (swidth .ge. 0.0d0) then
            delta = swidth
          else
            delta = 2.0d0* radmax/(levelo/2.0d0+1) + rprobej
          endif
          if (sdist .lt. 0.0d0) then
            sdis = radmax + rprobej
          else 
            sdis = sdist
          endif
c
c    -----  transform atom coordinates to have centre of gravity as orig
c           (because unit sphere has been defined w.r.t. (0,0,0)
c
          call cortran(forw,cgrav,nat,c)
          call cortran(forw,cgrav,nxtpts,xpts)
c
          call spoke(delta,sdis,maxnrs,nrs,xscm(ixrs),natmax,nat,
     1     c,maxpol,nxtpts,xpts)
c
c    -----  back-transform atom coordinates to original input values
c           the surface point coordinates can only be transformed
c           backwards after the calculation of the normal components
c           in compgen
c
          call cortran(back,cgrav,nat,c)
          call cortran(back,cgrav,nxtpts,xpts)
c   2-----
        else if (ibem .eq.4) then
c   2-----
c    -----  the surface is a single triangulized sphere of radius
c           -spherad- constructed from the unit sphere by scaling
c           the vertex coordinates
c
          write(iwr,113) spherad
  113     format(// 1x,' Triangulized spherical surface of ', f12.6, 
     1   ' bohr radius around (0,0,0)') 
          do 60, i = 1, 3*nrs
            xscm(ixrs-1+i) = xscm(ixrs-1+i)*spherad
   60     continue
c   2-----
        endif
c 1-----
      endif
c
c-----  the molecular surface-points have been determined (rs)
c       next, depending on the output-flag ibemout, a number of
c       properties of the surface are calculated.
c       in any case, the triangle-midpoints, the normal components
c       at the triangle midpoints and the triangle areas are computed.
c
c       other properties that may be computed are:
c       - unit normal vectors at the vertices  (norms)
c
      lmid = .true.
c
c-----  set pointers
c
c-----  space is already reserved for rs and surf; in fact, these
c       have been calculated and are used further on.
c
      ixinv = ixsurf + 3 * maxntr
c
c-----  reserve space for involv(maxinv,maxnrs)
c
      ixrm = ixinv + maxinv * maxnrs
c
c-----  reserve space for rm(3,maxntr)
c
      ixnm = ixrm + 3 * maxntr
c
c-----  reserve space for nm(3,maxntr)
c
      ixarea = ixnm + 3 * maxntr
c
c-----  reserve space for area(maxntr)
c
      ixnorm = ixarea + maxntr
c
      ixedge = ixnorm + 3
c
      last = ixedge + 3 + 1
c
c-----  determine output requests
c
      if (ibemout .ge. 5) lmidout = .true.
c
c-----  if (lmidout) midpoint vectors,
c       normal vectors at midpoints
c       and triangle areas will be written on output
c
      if (ibemout .ge. 10) then
c 1-----
c  -----  vertex coordinates, normal vectors and involvements
c         of the enveloping surface are printed on output
c
        lcni = .true.
c
c  -----  adjust pointer(s)
c
c  -----  reserve space for norms(3,maxnrs)
c
        ixedge = ixnorm + 3 * maxnrs
c
        last = ixedge + 3 + 1
c
        if (ibemout .ge. 100) then
c   2-----
c    -----  triangle pointers, length of edges and areas are printed
c
          lpla = .true.
c
c    -----  adjust pointer(s)
c
c    -----  reserve space for edge(3,maxntr)
c
          last = ixedge + 3 * maxntr + 1
c   2-----
        endif
c 1-----
      endif
c
c     need = loc10 + last - i10
      need = last - i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for drfsurf'
	call caserr('Insufficient memory for drfsurf')
      endif
c     call setc(need)
c
      call compgen(lcni,lpla,lmid,luni,minedge,maxedge,maxinv,nrs,ntr,
     1 xscm(ixrs),xscm(ixsurf),xscm(ixinv),xscm(ixrm),
     2 xscm(ixnm),xscm(ixarea),tarea,xscm(ixnorm),xscm(ixedge))
c
c-----  back-transform midpoint and vertex vectors to atom-input
c       coordinate system
c
      call cortran(back,cgrav,ntr,xscm(ixrm))
      call cortran(back,cgrav,nrs,xscm(ixrs))
c
c-----  write calculated data for molecular surface to output
c
      call writgen(lcni,lpla,lmid,luni,lmidout,nrs,ntr,minedge,
     1     maxedge,maxinv,
     2     xscm(ixrs),xscm(ixsurf),xscm(ixinv),xscm(ixrm),
     3     xscm(ixnm),xscm(ixarea),tarea,xscm(ixnorm),xscm(ixedge))
c
c-----  check surface: atoms within?, distance atoms-surface, etc..
c
c-----  adjust pointers; introduce new pointers
c
c        ixouts : begin of iouts(nat+nxtpts) array; determines which
c                 atoms are outside enveloping surface
c        ixshell: begin of ishell(nat+nxtpts) array; determines which
c                 atoms are within maxedge distance from surface
c        ixclose: begin of iclose(nat+nxtpts) array; determines the
c                 vertex closest to each atom
c        ixdclos: begin of dclose(nat+nxtpts) array; the distance
c                 between each atom and its closest vertex
c
c-----  reserve space for all new arrays, everything from area on
c       will be overwritten
c
      ixouts = ixarea + 3 * maxntr
      ixshell = ixouts + nat + nxtpts
      ixclose = ixshell + nat + nxtpts
      ixdclos = ixclose + nat + nxtpts
      last = ixdclos + nat + nxtpts + 1
c
c     need = loc10 + last - i10
      need = last - i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for drfsurf'
	call caserr('Insufficient memory for drfsurf')
      endif
c     call setc(need)
c
      call chk(nouts,nshell,maxedge,maxinv,nrs,natmax,nat,c,
     1 maxpol,nxtpts,xpts,xscm(ixrs),xscm(ixsurf),xscm(ixinv),
     2 xscm(ixrm),xscm(ixnm),xscm(ixouts),xscm(ixshell),xscm(ixclose),
     3 xscm(ixdclos))
c
c-----  write general output for surface construction to standard
c       output according to ibemout
c       if ibemout .ge.2, information about the position of the
c       atoms w.r.t. the surface is written, otherwise suppressed
c
      call writgms(ibemout,nouts,nshell,minedge,maxedge,delta,
     $             sdis,nrs,ntr,tarea,nat,nxtpts,xscm(ixouts),
     $             xscm(ixshell),xscm(ixclose),xscm(ixdclos))
c
c-----  write essential data to drf da-file
c       these are: midpoint coordinates (rm)
c                  midpoint normal components (nm)
c                  triangle areas (area)
c
  200 call dawrit(idafdrf,iodadrf,xscm(ixrm),3*ntr,ilxs,navdrf)
      call dawrit(idafdrf,iodadrf,xscm(ixnm),3*ntr,ilxn,navdrf)
      call dawrit(idafdrf,iodadrf,xscm(ixarea),ntr,ilar,navdrf)
c
c-----  reset core memory
c
c     call setc(loadcm)
      call gmem_free(i10)
c
c-----  set dimension of boundary problem
c
  300 if (ieps .eq. 0) then
        nbem1 = nbem
      else
        nbem2 = nbem
      endif
c
      if (ibemout .gt. 0) write(iwr,301)
  301 format(/1x,104('-')/)
c
      return
      end
      subroutine conread(isur,ibem,pdens,nrc,rc,tnorm,area)
c------
c
c       this subroutine takes care of the reading in and conversion
c       to a.u. of the connolly molecular surface points. it reads
c       the point-density from unit99 and the molecular surface
c       points from unit 42 or 43.
c
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  declaration of dummy variables
c
      integer nrc,ibem
      REAL pdens
      dimension rc(3,*),tnorm(3,*),area(*)
c
INCLUDE(comdrf/iofil)
c
c-----  declaration of local variables
c
      integer n,i,j
c
c-----  read the point-density from the connolly control file
c       unit99. ultimately, the point-density pdens determines
c       the cylinder radius around a spoke (in drfsurf)
c
      open (unit=99,form='formatted',file='coninp',
     +        status='unknown')
      rewind (99)
      read (99,fmt=*) pdens
      close (99,status='keep')
c
c-----  taken from juffer's subroutine rdgms
c
      n = 0
  140 n = n + 1
c     if (n .gt. maxnrc) goto 150
      if (ibem .eq. 5) then
c
c  -----  long output format is assumed
c
        read(isur,51,end=148)
     1                        rc(1,n),rc(2,n),rc(3,n),area(n),
     1                        tnorm(1,n),tnorm(2,n),tnorm(3,n)
      else
c
c  -----  short output format is assumed (only coordinates)
c
        read(isur,51,end=148)
     1                        rc(1,n),rc(2,n),rc(3,n)
      endif
  51  format(7e12.6)
      goto 140
  148 continue
      nrc = n - 1
      goto 152
c 150 continue
c     nrc = n
  152 continue
      write (iwr,57)
  57  format(1x,'       nrc' /)
      write (iwr,15) nrc
  15  format(1x,i10,2x,i10)
c
      return
      end
      subroutine bemraws(rs,isurf,maxnrs,maxntr,rad2,
     +                   new,nrsl)
c------
c     --------- defines a raw set of molecular surface points
c               by putting points on spheres around atoms
c     --------- the raw surface points are written unformatted
c               to unit 32 in batches of 170 points (3*170)
c               coordinates
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy variables
c
      dimension rs(3,maxnrs),isurf1(5,12),isurf(3,maxntr),craw(3,170)
      dimension rad2(*), new(nrsl,nrsl)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/extinf)
c
INCLUDE(comdrf/rad)
INCLUDE(comdrf/drfbem)
c
c-----  local variables
c
      dimension rss(3)
c
      REAL one, radius
c
      data radius/1.d00/
      data one /1.d00/
c
c-----  begin
c-----  initialize nrs and ntr
c
      nrs = 0
      ntr = 0
c
c-----  open unit 32 for writing raw surface points
c
      open (unit=32,form='unformatted',file='rawsurf',
     +        status='unknown')
c
c-----  prepare a list of all atom radii**2
c
      npoints = nat + nxtpts
      do 50, n = 1, npoints
        if (n .le. nat) then
          rad2(n) = radat(n)**2
        else
          rad2(n) = radext(n-nat)**2
        endif
   50 continue
c
c-----  first put some(level) points on a unit sphere
c
c               level=0: 32
c               level=1: 122
c               level=2: 482
c
      call setups(one,maxnrs,rs)
      call initsu(isurf1)
      call level1(radius,isurf1,maxnrs,nrs,rs,maxntr,ntr,isurf)
      if(leveli.gt.0)
     1  call levelh(radius,leveli,maxnrs,nrs,rs,maxntr,ntr,isurf,
     +              new,nrsl)
c
c-----  scan all atoms
c
      nraw = 0
      nrw = 0
      do 2000, np1 = 1, npoints
        if (np1 .eq. npoints) then
           nlast = npoints-1
        else
           nlast = npoints
        endif
c
c  -----  only "real" atoms
c
        if (np1 .le. nat) then
c
c    -----  internal atoms
c
           xp1 = c(1,np1)
           yp1 = c(2,np1)
           zp1 = c(3,np1)
           rad = radat(np1)
        else
c
c    -----  external atoms
c
           if (nxcent(np1-nat) .eq. 'qq') goto 2000
           xp1 = xpts(1,np1-nat)
           yp1 = xpts(2,np1-nat)
           zp1 = xpts(3,np1-nat)
           rad = radext(np1-nat)
        endif
c
c  -----  scale rs to wanted radius ands shift to atom center
c
        do 1000, n = 1, nrs
          rss(1) = rs(1,n)*rad + xp1
          rss(2) = rs(2,n)*rad + yp1
          rss(3) = rs(3,n)*rad + zp1
c
c    -----  scan atoms again
c
          do 500, np2 = 1, npoints
c
c      -----  exclude atom np2=np1 (senseless check)
c
            if ((np2 .eq. np1) .and. (npoints.ne.1)) goto 500
c
c      -----  calculate distance between atom coordinate
c             and surface-point-in-spe (rss)
c
            dis2 = 0.0d0
            do 100, k = 1, 3
c
              if (np2 .le. nat) then
                dis2 = dis2 + (rss(k) - c(k,np2))**2
              else
                dis2 = dis2 + (rss(k) - xpts(k,np2-nat))**2
              endif
  100       continue
c
c      -----  get rid of points in the molecular volume
c
            if (dis2 .lt. rad2(np2)) then
c
c        -----  point is not on the surface; get next point rss
c
              goto 1000
            else
c
c        -----  point may be on the surface, but second loop over
c               atoms must be finished
c
              if (np2 .ge. nlast) then
c
c          -----  last atom in second loop has been reached:
c                 point is indeed on the surface: count and store
c
                nrw = nrw + 1
                craw(1,nrw) = rss(1)
                craw(2,nrw) = rss(2)
                craw(3,nrw) = rss(3)
                if (nrw .eq. 170) then
                  write(32) craw
                  nraw = nraw + nrw
                  nrw = 0
                endif
c
c          -----  get next surface point
c
                goto 1000
c
              endif
            endif
  500     continue
 1000   continue
 2000 continue
      write (32) craw
      nraw = nraw+nrw
c
      close (unit=32, status='keep')
c
      return
      end
      subroutine setups(radius,maxnrs,rs)
c------
c subroutine setups determines the coordinates of the edge-points rs
c of the pentagons of a dodecahedron.
c the calculation is based on the coordinates of a cube.
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  declarations of dummy arguments.
c
      integer maxnrs
      dimension rs(3,maxnrs)
c
c-----  begin
c
c-----  first compute the length of a side of a pentagon (a) and a cube
c
      a=2.0d0*radius/sqrt(3.0d0)
      b=4.0d0*radius/((1.0d0+sqrt(5.0d0))*sqrt(3.0d0))
c
c-----  d is just for convenience.
c
      d=b*sqrt((3.0d0-sqrt(5.0d0))/8.0d0)
c
c-----  now compute the coordinates.
c
      rs(1,1)=0.5d0*a-d
      rs(2,1)=0.0d0
      rs(3,1)=0.5d0*a+0.5d0*b
      rs(1,2)=0.5d0*a
      rs(2,2)=0.5d0*a
      rs(3,2)=0.5d0*a
      rs(1,3)=0.0d0
      rs(2,3)=0.5d0*a+0.5d0*b
      rs(3,3)=0.5d0*a-d
      rs(1,4)=-0.5d0*a
      rs(2,4)=0.5d0*a
      rs(3,4)=0.5d0*a
      rs(1,5)=-0.5d0*a+d
      rs(2,5)=0.0d0
      rs(3,5)=0.5d0*a+0.5d0*b
      rs(1,6)=0.0d0
      rs(2,6)=0.5d0*a+0.5d0*b
      rs(3,6)=-0.5d0*a+d
      rs(1,7)=0.5d0*a
      rs(2,7)=0.5d0*a
      rs(3,7)=-0.5d0*a
      rs(1,8)=0.5d0*a-d
      rs(2,8)=0.0d0
      rs(3,8)=-0.5d0*a-0.5d0*b
      rs(1,9)=-0.5d0*a+d
      rs(2,9)=0.0d0
      rs(3,9)=-0.5d0*a-0.5d0*b
      rs(1,10)=-0.5d0*a
      rs(2,10)=0.5d0*a
      rs(3,10)=-0.5d0*a
      rs(1,11)=0.5d0*a
      rs(2,11)=-0.5d0*a
      rs(3,11)=-0.5d0*a
      rs(1,12)=-0.5d0*a
      rs(2,12)=-0.5d0*a
      rs(3,12)=-0.5d0*a
      rs(1,13)=0.0d0
      rs(2,13)=-0.5d0*a-0.5d0*b
      rs(3,13)=-0.5d0*a+d
      rs(1,14)=0.5d0*a
      rs(2,14)=-0.5d0*a
      rs(3,14)=0.5d0*a
      rs(1,15)=0.0d0
      rs(2,15)=-0.5d0*a-0.5d0*b
      rs(3,15)=0.5d0*a-d
      rs(1,16)=-0.5d0*a
      rs(2,16)=-0.5d0*a
      rs(3,16)=0.5d0*a
      rs(1,17)=0.5d0*a+0.5d0*b
      rs(2,17)=-0.5d0*a+d
      rs(3,17)=0.0d0
      rs(1,18)=0.5d0*a+0.5d0*b
      rs(2,18)=0.5d0*a-d
      rs(3,18)=0.0d0
      rs(1,19)=-0.5d0*a-0.5d0*b
      rs(2,19)=-0.5d0*a+d
      rs(3,19)=0.0d0
      rs(1,20)=-0.5d0*a-0.5d0*b
      rs(2,20)=0.5d0*a-d
      rs(3,20)=0.0d0
c
      return
      end
      subroutine initsu(surf1)
c------
c subroutine initsu sets up the pointer surf1 for the pentagons. it
c points to five vertices, whose components are in rs, which were
c calculated in subroutine setups.
c note that surf1(j,i)=k means that edge-point j of pentagon i is
c vertex k. note also that rs(1,k) is the x-component of vertex k.
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  declarations of dummy arguments.
      integer surf1(5,12)
c
c-----  begin
c
      surf1(1,1)=1
      surf1(2,1)=2
      surf1(3,1)=3
      surf1(4,1)=4
      surf1(5,1)=5
      surf1(1,2)=5
      surf1(2,2)=4
      surf1(3,2)=20
      surf1(4,2)=19
      surf1(5,2)=16
      surf1(1,3)=5
      surf1(2,3)=16
      surf1(3,3)=15
      surf1(4,3)=14
      surf1(5,3)=1
      surf1(1,4)=1
      surf1(2,4)=14
      surf1(3,4)=17
      surf1(4,4)=18
      surf1(5,4)=2
      surf1(1,5)=8
      surf1(2,5)=7
      surf1(3,5)=6
      surf1(4,5)=10
      surf1(5,5)=9
      surf1(1,6)=9
      surf1(2,6)=10
      surf1(3,6)=20
      surf1(4,6)=19
      surf1(5,6)=12
      surf1(1,7)=9
      surf1(2,7)=12
      surf1(3,7)=13
      surf1(4,7)=11
      surf1(5,7)=8
      surf1(1,8)=8
      surf1(2,8)=11
      surf1(3,8)=17
      surf1(4,8)=18
      surf1(5,8)=7
      surf1(1,9)=18
      surf1(2,9)=7
      surf1(3,9)=6
      surf1(4,9)=3
      surf1(5,9)=2
      surf1(1,10)=6
      surf1(2,10)=10
      surf1(3,10)=20
      surf1(4,10)=4
      surf1(5,10)=3
      surf1(1,11)=17
      surf1(2,11)=14
      surf1(3,11)=15
      surf1(4,11)=13
      surf1(5,11)=11
      surf1(1,12)=13
      surf1(2,12)=12
      surf1(3,12)=19
      surf1(4,12)=16
      surf1(5,12)=15
c
      return
      end
      subroutine level1(radius,surf1,maxnrs,nrs,rs,maxntr,ntr,surf)
c------
c subroutine level1 sets up the pointer for the triangles. this pointer
c is constructed from surf1. therefore, each pentagon is split up in
c five triangles. this is done by 'connecting' the midpoint of the
c pentagon, whose radial vector is corrected to the value of radius,
c with the edge-points of the pentagons.
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  declaration of dummy arguments.
c
      integer surf1(5,12),surf(3,maxntr)
      dimension rs(3,maxnrs)
c
c----- declarations of local variables.
c
      dimension rm(3)
c
c-----  begin
c
c-----  for all pentagons i do:
c
      do 10,l=1,12
c
c  -----  calculate position-vector of midpoint of pentagon.
c
        do 20,k=1,3
          rm(k)=0
20      continue
        do 30,i=1,5
          j=surf1(i,l)
          do 40,k=1,3
            rm(k)=rm(k)+rs(k,j)
40        continue
30      continue
        do 50,k=1,3
          rm(k)=rm(k)/5.0d0
50      continue
c
c  -----  the length of the radial vector of rm has to be corrected to t
c         radius of the sphere.
c
        call corrm(radius,rm)
c
c  -----  copy new coordinates into rs. set up the new triangles. set up
c         pointer surf.
c
        j=20+l
        do 60,k=1,3
          rs(k,j)=rm(k)
60      continue
        i5=5*l-5
        surf(1,i5+1)=surf1(1,l)
        surf(2,i5+1)=surf1(2,l)
        surf(3,i5+1)=j
        surf(1,i5+2)=surf1(2,l)
        surf(2,i5+2)=surf1(3,l)
        surf(3,i5+2)=j
        surf(1,i5+3)=surf1(3,l)
        surf(2,i5+3)=surf1(4,l)
        surf(3,i5+3)=j
        surf(1,i5+4)=surf1(4,l)
        surf(2,i5+4)=surf1(5,l)
        surf(3,i5+4)=j
        surf(1,i5+5)=surf1(5,l)
        surf(2,i5+5)=surf1(1,l)
        surf(3,i5+5)=j
c
10    continue
c
c-----  set initial value for nrs and ntr.
c
      nrs=32
      ntr=60
c
      return
      end
      subroutine corrm(radius,rm)
c------
c subroutine corrm corrects the components of an position-vector rm
c such that the length of rm cooresponds to the radius of a sphere.
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  declarations of dummy arguments.
c
      dimension rm(3)
c
c-----  declarations of local variables.
c
      REAL lrm
c
c-----  begin
c
c-----  calculate length lrm of rm and scalingfactor scalef.
c
      lrm=sqrt(rm(1)**2+rm(2)**2+rm(3)**2)
      scalef=radius/lrm
c
c-----  scale the length of rm to the radius of the sphere by multiplica
c       by scalef of its components.
c
      do 10,k=1,3
        rm(k)=rm(k)*scalef
10    continue
c
      return
      end
      subroutine levelh(radius,option,maxnrs,nrs,rs,maxntr,ntr,surf,
     +                  new,nrsl)
c------
c subroutine levelh generates from the triangles constructed in
c subroutine level1 more triangles. this is done by splitting each edge
c of each triangle in two equal parts. so, in each 'old' triangle four
c 'new' triangles are formed by 'connecting' the midpoints of the
c edges with the edge-points of the 'old' triangle.
c
c list of important local variables:
c   new(nrsl,nrsl) (integer): new(i,j)=l means: the midpoint of the side
c                             between vertex i and j is point l. point l
c                             is one of the points in rs also.
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  declarations of dummy arguments.
c
      integer option
      integer surf(3,maxntr)
****  surely surft is also integer ?
      integer surft
      dimension rs(3,maxnrs)
c
c-----  declarations of local variables.
c
      integer h
      integer nrsl
*     parameter (nrsl=2000)
      dimension new(nrsl,nrsl)
      parameter (ntrl=4000)
      dimension surft(3,ntrl)
      dimension rt(3,3),rme(3,3),rm(3)
c
c-----  begin
c
c-----  generate triangles until value of option is reached.
c
      do 10,h=1,option
c
c  -----  calculate current number of triangles ntrt.
c
        ntrt=(4**(h-1))*60
c
c  -----  set elements of new to zero.
c
        do 20,l=1,nrsl
          do 30,m=1,nrsl
            new(m,l)=0
  30      continue
  20    continue
c
c  -----  save original surf in surft.
c
        do 40,i=1,ntrt
          do 50,k=1,3
            surft(k,i)=surf(k,i)
  50      continue
  40    continue
c
c  ----- now construct new triangles.
c        for all triangles i do:
c
        do 60,i=1,ntrt
c
c    -----  get coordinates of vertices of triangle.
c
          do 70,l=1,3
            j=surft(l,i)
            do 80,k=1,3
              rt(k,l)=rs(k,j)
  80        continue
  70      continue
c
c    -----  split each side of the triangle in two parts.
c
          do 90,k=1,3
            rme(k,1)=rt(k,1)+0.5d0*(rt(k,2)-rt(k,1))
            rme(k,2)=rt(k,2)+0.5d0*(rt(k,3)-rt(k,2))
            rme(k,3)=rt(k,3)+0.5d0*(rt(k,1)-rt(k,3))
  90      continue
c
c    -----  correct the lengths of the radial vectors of the position-ve
c           of the midpoints of the edges of the triangle to the radius
c           sphere.
c
          do 100,l=1,3
            do 110,k=1,3
              rm(k)=rme(k,l)
  110       continue
            call corrm(radius,rm)
            do 120,k=1,3
              rme(k,l)=rm(k)
  120       continue
  100     continue
c
c    -----  check if some of the midpoints rme are already in rs.
c           if not, add to rs.
c
          if (new(surft(1,i),surft(2,i)) .eq. 0) then
            nrs=nrs+1
            new(surft(1,i),surft(2,i))=nrs
            new(surft(2,i),surft(1,i))=nrs
            do 130,k=1,3
              rs(k,nrs)=rme(k,1)
  130       continue
          endif
          if (new(surft(2,i),surft(3,i)) .eq. 0) then
            nrs=nrs+1
            new(surft(2,i),surft(3,i))=nrs
            new(surft(3,i),surft(2,i))=nrs
            do 140,k=1,3
              rs(k,nrs)=rme(k,2)
  140       continue
          endif
          if (new(surft(3,i),surft(1,i)) .eq. 0) then
            nrs=nrs+1
            new(surft(3,i),surft(1,i))=nrs
            new(surft(1,i),surft(3,i))=nrs
            do 150,k=1,3
              rs(k,nrs)=rme(k,3)
  150       continue
          endif
c
c    -----  now set up the new triangles. adjust the pointer surf to the
c           situation.
c
          i3=4*i-4
          surf(1,i3+1)=surft(1,i)
          surf(2,i3+1)=new(surft(1,i),surft(2,i))
          surf(3,i3+1)=new(surft(3,i),surft(1,i))
          surf(1,i3+2)=surft(2,i)
          surf(2,i3+2)=new(surft(2,i),surft(3,i))
          surf(3,i3+2)=new(surft(1,i),surft(2,i))
          surf(1,i3+3)=surft(3,i)
          surf(2,i3+3)=new(surft(3,i),surft(1,i))
          surf(3,i3+3)=new(surft(2,i),surft(3,i))
          surf(1,i3+4)=new(surft(1,i),surft(2,i))
          surf(2,i3+4)=new(surft(2,i),surft(3,i))
          surf(3,i3+4)=new(surft(3,i),surft(1,i))
c
  60    continue
  10  continue
c
c-----  calculate final number of triangles.
c
      ntr=(4**option)*60
c
      return
      end
      subroutine compgen(lcni,lpla,lmid,luni,minedge,maxedge,maxinv,
     $                   nrs,ntr,rs,surf,involv,rm,tnorm,area,tarea,
     $                   norms,edge)
c------
c subroutine compgen calculates unit normal vectors norms at the vertice
c on the surface, determines a list of triangles involv in which a
c vertex is involved, determines the area of the triangles area (for
c that purpose the triangles are considered to be flat) and calculates
c the length of the edges of each triangle edge.
c
c furthermore, it calculates the positions of the midpoints of the
c triangles (rm) and the normal vectors of the triangles (tnorm).
c these are the
c relevant quantities (together with the areas) for the constant element
c approach.
c
c involvements and edge-lengths are always calculated, but only stored
c and printed if lcni, resp. lpla are set .true.
c
c
c the stored quantities are governed by logicals
c     lcni: vertex coordinates, normal vectors and involvments are
c           calculated, stored and printed
c     lpla: triangle pointers, edge-lengths and areas are calculated
c           stored and printed
c     lmid: triangle midpoints, areas and normal components are calculat
c           and stored
c     luni: in case the unit sphere vertex normals are computed
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesahsurfsub)
c
c-----  declarations of dummy parameters.
c
      REAL minedge,maxedge
      REAL norms
      integer surf(3,*),involv(maxinv,*)
      dimension rs(3,*),norms(3,*),area(*),edge(3,*),
     $       tnorm(3,*),rm(3,*)
      logical lcni,lpla,lmid,luni
c
c-----  declarations of local variables.
c
      common/iofile/ir,iwr,ipnch
c
      dimension jvertex(3)
      REAL normlen
      REAL tedge(3),nm(3),rmt(3),rt1(3),rt2(3),rt3(3)
      REAL bmina(3),cmina(3),bminc(3)
c
      character *8 errmsg(3)
c
      data errmsg /'program', 'stop in', '-compgen'/
c
c-----  initialize
c
      minedge=1.0d+30
      maxedge=0.0d0
      tarea=0.0d0
c
c-----  for all triangles i do:
c
      do 10,i=1,ntr
c
c  -----  get coordinates of vertices of triangle.
c
c            rt1(3) : vertex a
c            rt2(3) : vertex b       of triangle abc
c            rt3(3) : vertex c
c
        do 20,l=1,3
          jvertex(l) = surf(l,i)
  20    continue
        do 30,k=1,3
          rt1(k)=rs(k,jvertex(1))
          rt2(k)=rs(k,jvertex(2))
          rt3(k)=rs(k,jvertex(3))
  30    continue
c
c  -----  calculate the difference vectors b-a, c-a and b-c
c         put these in vectors bmina,cmina and bminc
c
        call vecdiff(rt2,rt1,bmina,3)
        call vecdiff(rt3,rt1,cmina,3)
        call vecdiff(rt2,rt3,bminc,3)
c
c  -----  calculate the lengths of the triangle edges from the
c         difference vectors, these are returned in tedge
c
        call vlength(bmina,tedge(1),3)
        call vlength(cmina,tedge(2),3)
        call vlength(bminc,tedge(3),3)
c
c  -----  determine maximum and minimum edge-lengths
c
c  -----  if output is requested, store edge-lengths in edge
c
c
        do 40,l=1,3
          if (tedge(l) .gt. maxedge) maxedge=tedge(l)
          if (tedge(l) .lt. minedge) minedge=tedge(l)
c
          if (lpla) then
            edge(l,i)=tedge(l)
          endif
c
  40    continue
c
        if (lmid .or. lpla) then
c
c    -----  calculate area of element. calculate vector nm perpendicular
c           to plane containing the three vertices of triangle i,
c           and the triangle's midpoint rmt.
c
c    -----  this is done either because they are needed for the constant
c           element approach (lmid) or because the area and edge-lengths
c           are requested as output (lpla)
c
c    -----  calculate area of triangle via vector-product of two
c           edge-vectors.
c           the triangle area equals half the length of the
c           triangle's normal vector (the vector product of the edge-vec
c           the normalized normal vector is the traingle's requested
c           normal vector, stored in tnorm
c
          call vecprod(bmina,cmina,nm)
c
c    -----  calculate the length of the normal vector
c
          call vlength(nm,normlen,3)
c
c    -----  calculate and store the triangle area
c
          area(i) = 0.5d0 * normlen
c
c    -----  calculate total area
c
          tarea = tarea + area(i)
c
c    -----  only if requested, calculate and store triangle normal
c           components in tnorm
c    -----  as well as triangle midpoint vector components in rm
c
          if (lmid) then
c
c      -----  normalize the triangle normal vector
c
            call nlizv(nm,nm,normlen,3)
c
c      -----  calculate triangle midpoints from vertex coordinates
c
            call midcalc(rt1,rt2,rt3,rmt)
c
c      -----  check orientation of normal vector via dotproduct with
c             triangle midpoints
c
c           call dotprod(rmt,nm,dot,3)
            dot = ddot(3,rmt,1,nm,1)
            if (dot .lt. 0.0d0) then
              do 45, k = 1, 3
                nm(k) = -nm(k)
  45          continue
            endif
c
c      -----  store the normal vector components nm in tnorm
c
c      -----  store midpoint vector coordinates rmt in rm
c
            do 50, k = 1, 3
              tnorm(k,i) = nm(k)
              rm(k,i) = rmt(k)
  50        continue
c
c      -----  end of (lmid) if statement
c
          endif
c
c    -----  end of (lmid) .or. (lpla) if statement
c
        endif
c
c  -----  end of i do loop (over triangles)
c
  10  continue
c
c-----  determine list of triangles in which vertex j is involved.
c       calculate unit normal vector in vertex j.
c
c-----  for all vertices j do:
c
      do 60,j=1,nrs
c
c  -----  initialize vertex normal vector component array norms
c         if requested on output
c
        if (lcni) then
          do 70,k=1,3
            norms(k,j)=0.0d0
 70       continue
        endif
c
c  -----  initialize involv.
c
        involv(1,j)=0
        m=2
c
c  -----  for all triangles i do:
c
        do 80,i=1,ntr
c
c    -----  adjust involv if necessary.
c
          do 90,l=1,3
            if (j .eq. surf(l,i)) then
              if (m .gt. maxinv-1) then
                write(iwr,11) j,maxinv-1
 11             format(// 1x,'error: number of triangles in which',
     $                     ' vertex ',i4,' is involved > maximum',
     $                     ' (= ',i2,')' /)
                write(iwr,13)
 13             format(// 1x,'subroutine compgen aborted ',
     $   'with error status. increase maxinv in drfsurf. ')
                call hnderr(3,errmsg)
              else
c
c          -----  calculate normal vector in vertex j if requested
c
                involv(1,j)=involv(1,j)+1
                involv(m,j)=i
                m=m+1
                if (lcni) then
                  do 100,k=1,3
c
c              -----  in case of sphere, normal vector components are
c                     proportional to vertex vectors
c
                    if (luni) then
                      norms(k,j) = rs(k,j)
                    endif
c
c              -----  in case of general surface, an averaging procedure
c                     is followed to calculate vertex normal vector
c                     the vertex normal vector is the average of all
c                     triangle normal vectors in which the vertex is inv
c
                    if (lmid) then
                      norms(k,j)=norms(k,j)+tnorm(k,i)
                    endif
 100              continue
c
c            -----  end of (lcni) if statement
c
                endif
c
c          -----  end of (m .gt. maxinv-1) if statement
c
              endif
c
c        -----  end of (j .eq. surf(l,i)) if statement
c
            endif
c
c      -----  end of l do loop (over triangle vertex involvements)
c
 90       continue
c
c    -----  end of i do loop (over triangles)
c
 80     continue
c
c  ----- the length of the normal vector in vertex j must be one.
c        normalize the normal vector in vertex j
c
        if (lcni) then
          call normliz(norms(1,j),norms(1,j),3)
        endif
c
c  -----  end of j do loop (over vertices)
c
 60   continue
c
      return
      end
      subroutine writgen(lcni,lpla,lmid,luni,lmidout,nrs,ntr,
     1           minedge,maxedge,maxinv,
     2           rs,surf,involv,rm,tnorm,area,tarea,norms,edge)
c------
c subroutine writgen writes results of surface computations
c to output, according to printflags
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(comdrf/sizesahsurfsub)
c
c-----  declarations of dummy arguments.
c
      integer surf(3,*),involv(maxinv,*)
      REAL minedge,maxedge
      REAL rs(3,*),rm(3,*),tnorm(3,*),area(*),norms(3,*),edge(3,*)
      logical lcni,lpla,lmid,luni,lmidout
c
c-----  declarations of local variables.
c
      common/iofile/ir,iwr,ipnch
c
      integer i,j,k,l
c
c-----  write data surface.
c-----  write results formatted to standard output file.
c
      write(iwr,11)
  11  format(// 1x,'subroutine drfsurf results:' /)
      if (luni) then
        write(iwr,12)
  12    format(1x,'for unit sphere' /)
      endif
      if (lmid) then
        write(iwr,13)
  13    format(1x,'for final enveloping surface' /)
      endif
c
c-----  write number of vertices nrs and number of triangles ntr.
c-----  write maximum and minimum edge lengths
c
c-----  for final surface this is done in subroutine writgms
c
      if (luni) then
        write(iwr,14)
  14    format(1x,'       nrs         ntr' /)
        write(iwr,15) nrs,ntr
  15    format(1x,i10,2x,i10)
        write(iwr,33)
  33    format(// 1x,'     minedge       maxedge' /)
        write(iwr,35) minedge,maxedge
  35    format(1x,e12.5,2x,e12.5)
c
      endif
c
      if (lcni) then
c
c  -----  write vector components of vertices, unit normal vectors and n
c         of triangles in which vertex is involved.
c
        write(iwr,17)
  17    format(/ 1x,'vertex',2x,7x,'vector components vertex',7x,2x,1x,
     $              'vector components unit normal vector')
        write(iwr,19)
  19    format(1x,6x,2x,1x,'x-component',2x,1x,'y-component',2x,1x,
     $         'z-component',2x,1x,'x-component',2x,1x,'y-component',
     $         2x,1x,'z-component' /)
        do 100,j=1,nrs
          write(iwr,21) j,rs(1,j),rs(2,j),rs(3,j),norms(1,j),norms(2,j),
     $             norms(3,j)
100     continue
  21    format(1x,i6,6 (2x,e12.5))
c
c  -----  write total number of triangles in which certain vertex is inv
c         give list of triangles in which vertex is involved.
c
        write(iwr,23)
  23    format(/ 1x,'vertex',2x,'number of triangles',4x,'list of',
     $   ' triangles in which a vertex is involved. ' )
c
        do 110,j=1,nrs
          write(iwr,27) j,involv(1,j),(involv(l,j),l=2,involv(1,j)+1)
  110   continue
  27    format(1x,i6,2x,i19,4x,50 (i4,1x))
      endif
c
      if (lpla) then
c
c  -----  write pointer, lengths of edges and area of triangle.
c
        write(iwr,29)
  29    format(/ 1x,'triangle',2x,4x,'vertices',4x,2x,11x,
     $  'lengths of edges',11x,2x,8x,'area' /)
        do 120,i=1,ntr
          write(iwr,31) i,surf(1,i),surf(2,i),surf(3,i),edge(1,i),
     $             edge(2,i),edge(3,i),area(i)
  120   continue
  31    format(1x,i8,3 (2x,i4),4 (2x,e12.5))
c
      endif
c
      if (lmidout) then
c
c  -----  write midpoint vectors, normal vectors and triangle areas
c
        write(iwr,39)
  39    format(/ 1x,'triangle',7x,'vector components midpoint',7x,
     $           1x,'vector components unit normal vector',7x,1x,
     $           'triangle area')
        write(iwr,19)
        do 200,j=1,ntr
          write(iwr,41) j,rm(1,j),rm(2,j),rm(3,j),tnorm(1,j),
     $                   tnorm(2,j),tnorm(3,j),area(j)
  200   continue
  41    format(1x,i6,7 (2x,e12.5))
      endif
c
      return
      end
      subroutine cortran(way,corcent,number,coord)
c------
c ------ this routine shifts the coordinates of coord forward or
c        backward (way) with a vector corcent
c
c        forward: the coordinates to be transformed are defined with
c                 respect to (0,0,0) and become defined w.r.t. corcent
c        backward:the coordinates to be transformed are defined w.r.t.
c                 corcent and become defined w.r.t. (0,0,0)
c
c        if corcent should equal (0,0,0) then this transformation is
c        a waste of time and should be avoided.
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      character*8 way
c
      dimension corcent(3)
      dimension coord(3,number)
c
      common/iofile/ir,iwr,ipnch
c
c
c-----  check if transformation is needed
c
      if ((corcent(1) .eq. 0.0d0) .and. (corcent(2) .eq. 0.0d0)
     1   .and. (corcent(3) .eq. 0.0d0)) then
        write (iwr,1001)
 1001   format (//' coordinate shift skipped, because corcent ',
     $         'coincides with (0,0,0) ')
        return
      else
c
c  -----  determine the direction of coordinate shift
c
        if (way .eq. 'forward ') then
          factor = -1.0d0
        else
          factor = 1.0d0
        endif
c
c  -----  perform coordinate transformation
c
        do 100, i = 1, number
          do 200, k = 1, 3
            coord(k,i) = coord(k,i) + factor * corcent(k)
 200      continue
 100    continue
c
      endif
      return
      end
      subroutine spokec(delta,nrs,rs,nrc,iused,rc)
c------
c subroutine spokec uses the surface points from the connolly
c program or bemraws to define the macromolecule surface.
c
c------  the surface points are in array rc, the unit sphere surface
c        coordinates are in rs, and are replaced by the enveloping
c        surface point vertex coordinates.
c
c------  the integer array iused keeps track of the surface points
c        allocated to a spoke (i.e. copied into rs) to avoid
c        coincidence of enveloping surface vertices.
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  declarations of dummy parameters.
c
      dimension iused(*)
      dimension rs(3,*),rc(3,*)
c
c-----  declarations of local variables.
c
      common/iofile/ir,iwr,ip
c
      REAL lmax,lrs,imp,lrc,lt
      logical found
c
      character *8 errmsg(3)
c
      data errmsg /'program','stop in','-spokec-'/
c
c----  initialize. no surface-point has been used yet.
c
      do 10,n=1,nrc
        iused(n)=0
10    continue
c
c-----  calculate length spoke.
c-----  as the initial surface is a sphere, the length of all
c        spokes is the same
c
      call vlength(rs(1,1),lrs,3)
c
c-----  for all spokes j do:
c
      do 20,j=1,nrs
c
c  -----  initialize.
c
        found=.false.
        lmax=0
c
c  -----  for all surface-points n do:
c
        do 40,n=1,nrc
c
c    -----  if not used yet.
c
          if (iused(n) .eq. 0) then
c
c      -----  alternative search method
c        (1) check quadrant of atom position
c        (2) calculate distance to spoke
c        (3) determine furthest away
c
c      -----  determine if atom n is within cylinder around spoke with
c             radius delta.
c             this is determined via the dotproduct of the atom and spok
c             vectors, and the length of the atom vector
c
c         if (rs(1,j) .ne. 0.0) then
c           if (sign(1.0,rc(1,n)) .ne. sign(1.0,rs(1,j))) goto 40
c         endif
c         if (rs(1,j) .ne. 0.0) then
c           if (sign(1.0,rc(2,n)) .ne. sign(1.0,rs(2,j))) goto 40
c         endif
c         if (rs(1,j) .ne. 0.0) then
c           if (sign(1.0,rc(3,n)) .ne. sign(1.0,rs(3,j))) goto 40
c         endif
c
c------  the atom is in the right quadrant
c
c      -----  determine if surface-point n is within cylinder
c             around spoke with radius delta.
c             this is done through the dotproduct of the surface point
c             and spoke vectors and their lengths
c
c           call dotprod(rs(1,j),rc(1,n),imp,3)
            imp = ddot(3,rs(1,j),1,rc(1,n),1)
c
            if (imp .gt. 0) then
c
c        -----  surface-point n is at the right side with respect to spo
c
              call vlength(rc(1,n),lrc,3)
              cosinus=imp/(lrs*lrc)
c
c        -----  lt is the length of the projection of rc on rs
c
              lt=lrc*cosinus
c
c        -----  dis is the distance of rc to the spoke through rs
c
              dis=lrc**2-lt**2
              if (dis .lt. delta**2) then
c
c          -----  surface-point n is within cylinder.
c
                if (lrc .ge. lmax) then
c
c            -----  surface-point is furthest away.
c                   replace coordinates spoke by surface-point,
c                   if not used already by other spoke.
c
                  lmax=lrc
                  found=.true.
                  nf=n
                endif
              endif
            endif
          endif
c
40      continue
c
        if (.not. (found)) then
          write(iwr,11) j
  11      format(// 1x,'error: no surface-point within range ',
     $          'delta for spoke ', i3 /)
          write(iwr,13)
  13      format(1x,'       delta' /)
          write(iwr,15) delta
  15      format(1x,e12.5)
          write(iwr,17)
  17      format(/,1x,'probably too few atoms for definition ',
     $ 'of Juffer surface; you may setting cylwdth ',
     $ '(dielectric subdirective juffer) to increase hit-rate')
          call hnderr(3,errmsg)
        else
          do 60,k=1,3
            rs(k,j)=rc(k,nf)
60        continue
          iused(nf)=1
        endif
c
c  -----  the choice has been made.
c
20      continue
c
      return
      end
      subroutine spoke(delta,sdis,maxnrs,nrs,rs,natmax,nat,c,
     1                 maxpol,nxtpts,xpts)
c------
c subroutine spoke uses atoms to define the
c macromolecule surface.
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  declarations of dummy parameters.
c
      common/iofile/ir,iwr,ipnch
c
      dimension rs(3,maxnrs),c(3,natmax)
      dimension xpts(3,maxpol)
c
c-----  declarations of local variables.
c
      REAL lmax,lrs,imp,lx,lt2,lt3,lt
      logical found
c
      character *8 errmsg(3)
c
      data errmsg /'program','stop in','-spoke-'/
c
c-----  begin
c
      npoints = nat + nxtpts
c
c-----  calculate length spoke.
c       as the initial surface is a sphere, the length of all
c       spokes is the same, so length is only calculated once
c
      call vlength(rs(1,1),lrs,3)
c
c-----  for all spokes j do:
c
      do 10,j=1,nrs
        lmax=0
        found=.false.
c
c  -----  for all atoms n (qm and classical)  do:
c
        do 30,n=1,npoints
c
c    -----  determine if atom n is within cylinder around spoke with
c           radius delta.
c           this is determined via the dotproduct of the atom and spoke
c           vectors, and the length of the atom vector
c
          if (n .le. nat) then
c
c           call dotprod(rs(1,j),c(1,n),imp,3)
            imp = ddot(3,rs(1,j),1,c(1,n),1)
            call vlength(c(1,n),lx,3)
          else
c
c           call dotprod(rs(1,j),xpts(1,n-nat),imp,3)
            imp = ddot(3,rs(1,j),1,xpts(1,n-nat),1)
            call vlength(xpts(1,n-nat),lx,3)
          endif
c
          cosinus=imp/(lrs*lx)
          if (cosinus .gt. 0) then
c
c      -----  atom n is at the right side with respect to spoke j.
c
            lt2=lx*cosinus
            dist=sqrt(lx**2-lt2**2)
            if (dist .lt. delta) then
c
c        -----  atom n is within cylinder.
c
              if (lx .gt. lmax) then
                lt=lt2
                dis=dist
                lmax=lx
                found=.true.
              endif
            endif
          endif
c
30        continue
c
        if (.not. (found)) then
          write(iwr,11) delta,j
  11      format(// 1x,'error: no atom within range delta=',
     $           f10.5,' for ',
     $           'spoke ',i3,/)
          write(iwr,13)
  13      format(// 1x,'subroutine spoke aborted with error status.')
          write(iwr,17)
  17      format(/,1x,'probably too few atoms for definition ',
     $ 'of Juffer surface; you may try ',/,1x,
     $ '1. setting cylwdth >< delta ',
     $ '(dielectric subdirective juffer) or',/,1x,
     $ '2. enlarging rprobe ',
     $ '(dielectric subdirective juffer) or',/,1x,
     $ '3. decreasing bemlev (dielectric subdirective) ',/,
     $ 1x,'to increase hit-rate')
          call hnderr(3,errmsg)
          call hnderr(3,errmsg)
        endif
c
c  -----  the choice has been made. scale rs such that no atom within
c         cylinder is closer to vertex than sdis
c         (=rwater+ratom(maximum value)).
c         at present sdis = radmax + rprobej: can be changed in subroutin
c
        if (dis .gt. sdis) then
          lt3=0
        else
          lt3=sqrt(sdis**2-dis**2)
        endif
        lt=lt+lt3
        scalef=lt/lrs
        do 50,k=1,3
          rs(k,j)=scalef*rs(k,j)
50        continue
c
10      continue
c
      return
      end
      subroutine chk(nouts,nshell,maxedge,maxinv,nrs,natmax,nat,c,
     $               maxpol,nxtpts,xpts,rs,surf,involv,rm,tnorm,
     $               iouts,ishell,iclose,dclose)
c------
c subroutine chk determines if atoms are outside the generated surface.
c the triangle are considered to be flat.
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(comdrf/sizesahsurfsub)
c
c-----  declarations of dummy arguments.
c
      integer surf(3,*)
      dimension involv(maxinv,*),iouts(*),iclose(*),ishell(*)
      REAL maxedge
      dimension rs(3,*),rm(3,*),tnorm(3,*),dclose(*)
      dimension xpts(3,maxpol),c(3,natmax)
c
c-----  declarations of local variables.
c
      REAL imp
      dimension disv(3),rshort(3)
c
c-----  initialize.
c
      nouts=0
      nshell=0
      npoints = nat + nxtpts
c
      do 10,n=1,npoints
        iouts(n)=0
        ishell(n)=0
10    continue
c
c-----  determine shortest distance to vertices.
c-----  for all atoms n do:
c
      do 20,n=1,npoints
c
c  -----  just a startvalue for shortest distance to surface.
c
        dclose(n)=1.0d30
c
c  -----  for all vertices j do:
c
        do 30,j=1,nrs
c
c    -----  calculate distance from vertex j to atom n.
c
          dis=0.0d0
          do 40,k=1,3
            if (n .le. nat) then
              dis=dis+(rs(k,j)-c(k,n))**2
            else
              dis=dis+(rs(k,j)-xpts(k,n-nat))**2
            endif
40        continue
c
c    -----  is distance smaller then shortest distance till now?
c
          if (dis .lt. dclose(n)) then
c
c      -----  save data of vertex.
c
            dclose(n)=dis
            iclose(n)=j
          endif
c
c    -----  next vertex.
c
30      continue
c
c  -----  is atom n within shell with depth maxedge?
c
        dclose(n) = sqrt(dclose(n))
        if (dclose(n) .lt. maxedge) then
          ishell(n)=1
          nshell=nshell+1
        endif
c
c  -----  next atom.
c
20    continue
c
c-----  which atoms are outside?
c-----  for all atoms n do:
c
      do 50,n=1,npoints
c
c  -----  initialize.
c
        dism=1.0d+30
c
c  -----  get closest vertex.
c
        j=iclose(n)
c
c  -----  for all triangles of which vertex j is an edge-point do:
c
        do 60,m=2,involv(1,j)+1
c
c    -----  get triangle.
c
          i=involv(m,j)
c
c    -----  calculate distance to midpoint.
c
          dis=0.0d0
          do 110,k=1,3
            if (n .le. nat) then
              disv(k)=c(k,n)-rm(k,i)
            else
              disv(k)=xpts(k,n-nat)-rm(k,i)
            endif
            dis=dis+disv(k)**2
110       continue
c
          if (dis .lt. dism) then
            is=i
            dism=dis
            do 120,k=1,3
              rshort(k)=disv(k)
120         continue
          endif
c
60      continue
c
c  -----  calculate improduct of difference vector with normal of
c         nearest triangle (is).
c
c       call dotprod(rshort,tnorm(1,is),imp,3)
        imp = ddot(3,rshort,1,tnorm(1,is),1)
c
        if (imp .gt. 0.0) then
c
c    -----  atom n is outside with respect to triangle l. correct ishell
c           since atom is outside.
c
          iouts(n)=1
          ishell(n)=0
          nouts=nouts+1
          if (nshell .gt. 0) nshell=nshell-1
        endif
c
50    continue
c
      return
      end
      subroutine writgms(ibemout,nouts,nshell,minedge,maxedge,delta,
     $                 sdis,nrs,ntr,tarea,nat,nxtpts,iouts,ishell,
     $                 iclose,dclose)
c------
c subroutine writgms write general results of subroutine drfsurf
c and, if ibemout .ge.2, specific information about the position of the
c atoms w.r.t. the surface (as calculated in subroutine chk)
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  declarations of dummy parameters.
c
      dimension iouts(*),ishell(*),iclose(*)
      REAL minedge,maxedge
      dimension dclose(*)
c
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/drfdaf)
c
      character *8 errmsg(3)
c
      data errmsg /'program','stop in','-writgms'/
c
c
c-----  write results formatted to standard output file.
c
      write(iwr,11)
  11  format(// 1x,'subroutine drfsurf results:' /)
c
c-----  write number of vertices nrs and number of triangles ntr.
c
      write(iwr,13)
  13  format(1x,'       nrs         ntr' /)
      write(iwr,15) nrs,ntr
  15  format(1x,i10,2x,i10)
c
      write(iwr,33)
  33  format(// 1x,'     minedge       maxedge         delta',
     $             '          sdis         tarea' /)
      write(iwr,35) minedge,maxedge,delta,sdis,tarea
  35  format(1x,e12.5,4(2x,e12.5))
c
c-----  write information about location of atoms with respect to surfac
c
c-----  write number of atoms outside surface, number of atoms
c       within shell and inside surface.
c
      write(iwr,37)
  37  format(// 1x,'     nouts      nshell' /)
      write(iwr,15) nouts,nshell
c
      if ((ibemout .ge. 2) .or. (nouts .gt. 0)) then
        write(iwr,39)
 39     format(/, ' internal atom')
        do 130,n=1,nat
          if (iouts(n) .eq. 1) then
            write(iwr,41) n,iclose(n),dclose(n)
  41        format(1x,i5,'  outside surface: nearest vertex ',i5,
     $              ' distance ',e12.5)
          else
            if (ishell(n) .eq. 1) then
              write(iwr,43) n,iclose(n),dclose(n)
  43          format(1x,i5,'  inside shell: nearest vertex ',i5,
     $                ' distance ',e12.5)
            else
              write(iwr,45) n,iclose(n),dclose(n)
  45          format(1x,i5,'                nearest vertex ',i5,
     $                ' distance ',e12.5)
            endif
          endif
130     continue
c
        write (iwr,139)
 139    format (/,' classical atom')
        do 140,n=1,nxtpts
          if (iouts(n+nat) .eq. 1) then
            write(iwr,41) n,iclose(n+nat),dclose(n+nat)
          else
            if (ishell(n+nat) .eq. 1) then
              write(iwr,43) n,iclose(n+nat),dclose(n+nat)
            else
              write(iwr,45) n,iclose(n+nat),dclose(n+nat)
            endif
          endif
140     continue
      endif
c
c-----  if there are atoms outside the surface, then the calculation
c       can't be right: get out of the program altogether and
c       make some suggestions
c
      if (nouts .gt. 0) then
        write(iwr,47) nouts
  47    format (// 1x,i4,' atoms found outside the surface:',/,
     1  1x,' the calculation cannot be right. improve the description',
     2  ' of the surface.')
        write(iwr,49)
  49    format (/ 1x,' some suggestions:',/,
cxxx 1  1x,' (1) decrease the cylinder radius around a spoke in',
cxxx 2  ' subroutine drfsurf - data scale',/,
cxxx 3  1x,' (2) increase the point density on the raw/connolly',
cxxx 4  ' surface: leveli in $bem/ in connolly program.')
     3  1x,' (1) set cylwdth >< delta ',
     4  '(dielectric subdirective juffer)',/,
     5  1x,' (2) set surfdist > sdis ',
     4  '(dielectric subdirective juffer)',/,
     4  1x,' (3) increase number of surface points',
     4  '(dielectric subdirective bemlev)')
c
        call hnderr(3,errmsg)
c
      else
c
c  -----  write final message.
c
        write(iwr,51)
  51    format(/// 1x,'subroutine drfsurf ran succesfully without',
     $              ' errors.' ///)
      endif
c
      return
      end
      subroutine vlength(a,alen,veclen)
c------
c       subroutine calculates the length of a vector a of
c       dimension veclen
c       and returns it in alen
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy variables
c
      integer veclen
      REAL a, alen
c
      dimension a(veclen)
c
c-----  local variables
c
      integer i
c
c-----  initialize alen
c
      alen = 0.0d0
c
c-----  loop over vector components
c
      do 10, i = 1, veclen
        alen = alen + a(i)**2
  10  continue
c
c-----  take square root for length
c
      alen = sqrt(alen)
c
      return
      end
_IF(notused)
      subroutine vecsum(a,b,c,veclen)
c------
c       subroutine calculates the sum
c       of two vectors a and b, both assumed to be of
c       dimension veclen and returns the difference in c
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy variables
c
      integer veclen
c
      dimension a(veclen), b(veclen), c(veclen)
c
c-----  local variables
c
      integer i
c
c-----  loop over vector components
c
      do 10, i = 1, veclen
c
c  -----  calculate sum of a and b
c
        c(i) = a(i) + b(i)
  10  continue
c
      return
      end
_ENDIF
      subroutine vecdiff(a,b,c,veclen)
c------
c       subroutine calculates the difference
c       of two vectors a and b, both assumed to be of
c       dimension veclen and returns the difference in c
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy variables
c
      integer veclen
c
      dimension a(veclen), b(veclen), c(veclen)
c
c-----  local variables
c
      integer i
c
c-----  loop over vector components
c
      do 10, i = 1, veclen
c
c  -----  calculate sum of a and b
c
        c(i) = a(i) - b(i)
  10  continue
c
      return
      end
c      subroutine dotprod(a,b,dot,veclen)
cc------
cc       subroutine calculates the dot product of two
cc       vectors a and b, both assumed to be of length veclen
cc       dimension veclen and return it in dot
cc------
c      implicit REAL  (a-h,o-z),integer  (i-n)
cINCLUDE(../m4/common/sizes)
cINCLUDE(comdrf/sizesrf)
cc
cc-----  dummy variables
cc
c      integer veclen
cc
c      REAL dot
cc
c      dimension a(veclen),b(veclen)
cc
cc-----  local variables
cc
c      integer i
cc
cc-----  initialize dot
cc
c      dot = 0.0d0
cc
cc-----  loop over vector components
cc
c      do 10, i = 1, veclen
cc
cc  -----  calculate dot-product
cc
c        dot = dot + a(i) * b(i)
c  10  continue
c
c     return
c     end
      subroutine vecprod(a,b,c)
c------
c       subroutine calculates the vector (or outer) product
c       of vectors a and b, both assumed to be of
c       dimension 3 and returns it in vector c
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy variables
c
      dimension a(3), b(3), c(3)
c
c-----  calculate components of vector product , via
c       determinant method
c
      c(1) = a(2) * b(3) - a(3) * b(2)
c
      c(2) = a(3) * b(1) - a(1) * b(3)
c
      c(3) = a(1) * b(2) - a(2) * b(1)
c
      return
      end
      subroutine normliz(a,anorm,veclen)
c------
c       subroutine normalizes a vector a of
c       dimension veclen and returns the normalized
c       vector in anorm
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy variables
c
      integer veclen
c
      dimension a(veclen), anorm(veclen)
c
c
c-----  determine length of vector a
c
      call vlength(a,alen,veclen)
c
c-----  determine norm to scale with, equals 1/length
c
      scale = 1.0d0 / alen
c
c-----  loop over vector components
c
      do 10, i = 1, veclen
c
c  -----  scale vector components with 1/length
c
        anorm(i) = scale * a(i)
  10  continue
c
      return
      end
      subroutine nlizv(a,anorm,alen,veclen)
c------
c       subroutine normalizes a vector a of dimension
c       veclen, given that its length is alen
c       and returns the normalized vector in anorm
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy variables
c
      integer veclen
c
      REAL alen
c
      dimension a(veclen), anorm(veclen)
c
c-----  determine norm to scale with, equals 1/length
c
      scale = 1.0d0 / alen
c
c-----  loop over vector components
c
      do 10, i = 1, veclen
c
c  -----  scale vector components with 1/length
c
        anorm(i) = scale * a(i)
  10  continue
c
      return
      end
      subroutine midcalc(a,b,c,midvec)
c------
c       subroutine calculates the position vector that is
c       the geometric mean of vectors a,b and c,
c       supposed to be of dimension 3
c       and returns it in vector midvec
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy variables
c
      REAL midvec
c
      dimension a(3), b(3), c(3), midvec(3)
c
c
      data third/.33333333333333333333333333d00/
c
c-----  loop over vector components
c
      do 10, i = 1, 3
c
c  -----  add vector components of a,b and c
c
        midvec(i) = a(i) + b(i) + c(i)
c
c  -----  average vector components
c
        midvec(i) = third * midvec(i)
c
  10  continue
c
      return
      end
