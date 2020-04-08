      subroutine mscon(rc,tnorm,areas,
     1 rtppls,rad,ua,eva,py,ay,
     2 nua,itype,ias,iat,ico,icuptr,
     3 molyon,srs,
     4 nxatm,nelm,ncon,iturn)
c------
c molecular surface program
c ms
c
c december 16, 1983
c
c copyright c 1983
c by michael connolly
c
c written by michael connolly
c
c mailing address:
c
c michael connolly
c molecular biology department
c scripps clinic and research foundation
c 10666 north torrey pines rd.
c la jolla, ca 92037, u.s.a.
c
c phone number: (619) 455-9100  ext. 2838
c
c
c references:
c
c m.l. connolly, "solvent-accessible surfaces of proteins
c and nucleic acids", science, volume 221, pages 709-713 (1983).
c
c m.l. connolly, "analytical molecular surface calculation",
c journal of applied crystallography,
c volume 16, pages 548-558 (1983).
c
c
c the program and its input and output are briefly
c described below.  further information may be found
c in the runoff documents:
c
c ms: user's guide
c ms: maintenance manual
c
c this program may be freely distributed to anyone.
c
c it is written in fortran 77.
c
c ms calculates the molecular surface of a molecule
c given the coordinates of its atoms.  van der waals radii for
c the atoms and the probe radius must also be specified.
c
c the term molecular surface was introduced by f.m. richards
c (annual reviews of biophysics and bioengineering, 1977,
c pages 151-176) with a specific meaning.  the surface
c richards defined consists of two parts:  (1) the contact
c surface and (2) the reentrant surface.  he defines the
c contact surface to be that part of the van der waals
c surface of each atom which is accessible to a probe sphere
c of a given radius.  he defines the reentrant surface to be
c the inward-facing part of the probe sphere when it is
c simultaneously in contact with more than one atom.
c
c in implementing this definition i have found that there are
c two kinds of reentrant surface:  (1) concave reentrant
c surface, which is generated when the probe sphere
c simultaneously touches three atoms and (2) saddle-shaped
c reentrant surface, which is generated as the probe sphere
c rolls along the crevice between two atoms.
c i have also found that reentrant surface belonging to one
c probe may be contained in the interior volume of an
c overlapping probe and so must be removed.
c
c the input to this program consists of three files.
c
c the first file contains one record specifying the requested
c density of surface points (the number of surface points
c per unit area), the probe radius, the buried surface
c flag and the ascii/binary long/short output flag.
c normally the buried surface flag will be blank or zero.
c if this flag is equal to 1, then the
c only surface calculated will be the surface of each
c molecule that is buried by the other molecules.
c if this flag is equal to 2, both buried and unburied
c surface are calculated, but the buried surface points
c are flagged by a 1 in the buried output field to the
c right of the surface normal.
c the ascii/binary long/short output flag
c may have one of these four values:
c
c 0      ascii       long
c 1      binary      long
c 2      ascii       short
c 3      binary      short
c
c the various output formats are discussed below.
c
c the format of the parameter record is: (2f10.5,2i5)
c
c the second file contains atomic radius records.
c each record has an integer that is the atom type
c and a van der waals radius for that type.
c the atom types need not be contiguous or in
c increasing order. the format is: (i5,f10.5).
c
c the third file contains the atomic coordinate records.
c each atomic coordinate record has the x, y and z coordinates
c of the atoms, followed by the atom type, a surface

c the surface request number may be 0, 1 or 2.
c 0 means that the atom is to be ignored,
c except for occupying a place in the sequence of atom numbers.
c 1 means that no surface is to be generated for this atom,
c but probe positions that collide with this atom will still
c be disallowed.  2 means that we are requested to generate
c surface for this atom.  in most cases 2 should be specified
c for all atoms. the molecule number is an integer that is
c used to divide the atoms into groups so that each group
c may be given its own surface. these groups need not
c correspond to actual molecules, as the program knows
c nothing about bonding. the characters to the right of
c these six fields are not read and may contain anything.
c
c the output from this program consists of two files,
c called 'contact' and 'reentrant', containing the contact
c and reentrant surface points.  all lines in both files have
c the same format.  the first three fields are the atom
c numbers of the atoms the probe was touching when it
c generated the given surface point.  the first number is
c the atom whose van der waals surface the point is closest to.
c the fourth field is the number of atoms the probe was
c touching. if this number is less then three, one or both
c of the second and third fields will be zero. the fifth, sixth
c and seventh fields are the coordinates of the surface point.
c the eighth field is the molecular area associated with
c the surface point.  the ninth, tenth and eleventh fields
c are a unit vector perpendicular to the surface pointing in
c the direction of the probe center.
c if the buried surface flag equals 2, there is a '1' written
c at the end of the record for buried surface points, and
c a '0' for accessible points. the output format is:
c (3i5,i2,3f9.3,4f7.3,i2)
c there is also a short output format: (3i5,i2,3f9.3,i2).
c the corresponding binary records are:
c (4i*2,7r*4,i*2) and (4i*2,3r*4,i*2)
c
c the contact file is generated in atom number order.
c the reentrant file should be sorted
c into atom number order and then merged with the contact file
c using the sort/merge utility programs available at your
c installation.  then all the surface points belonging to a
c given atom will be in a contiguous series of records.
c
c example command file for vax-11/780:
c
c
c the program has the ability to calculate the
c van der waals surface
c of a molecule simply by specifying a probe radius of zero.
c the reentrant code is bypassed completely and there are no
c before and reentrant files. for a probe radius of zero
c the van der waals surface and the contact surface
c are equivalent.
c
c the flow of the program may be described in general terms
c as follows. first all the input is read. then the contact
c and reentrant surface is generated. the contact surface
c is in its final form, but the reentrant surface is
c written to a temporary file, called 'before'.
c each reentrant probe position is written to this file,
c followed by all its surface points.
c after all the contact and reentrant surface has been generated,
c the 'before' file is read and a final reentrant surface file is
c written which contains all the reentrant surface points
c not lying within any reentrant probe.
c
c
c sometimes it is desirable to calculate the surface of only
c part of a molecule. one cannot simply remove the remaining
c atoms from the input file as this will generate false
c surfaces. one could calculate a surface for the entire
c molecule and then edit out that part of the surface belonging
c to the atoms of interest, but this is a needless
c waste of computer time.  the proper way to accomplish this
c task is to place the probe next to the atoms whose surface
c is requested and to use the remaining atoms only for
c collision checks so that false surfaces will not be
c generated. when the probe is placed next to two or three
c atoms, only one of these need be an atom whose surface
c is requested. only that part of the arc or spherical
c triangle that belongs to the atoms whose surface is
c requested will be written as output.
c
c this program allows the calculation of individual surfaces
c for several interacting molecules. this may be done with
c or without the buried surface option.
c
c this program should be compiled with integer*4 as the default
c it uses no include files or libraries
c
c
c
cxxx 5/91
c      changes marked thus made by alex de vries (month/year)
c
c      general: logical*1 to logical
cxxx
      implicit REAL (a-h,o-z)
c-----
INCLUDE(comdrf/sizescon)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/connolly)
INCLUDE(comdrf/drfbem)
INCLUDE(comdrf/iofil)

      dimension rc(3,nelm),tnorm(3,nelm),areas(nelm)
c     dimension rtppls(maxatm)
      dimension rtppls(nxatm)
c
c atom type arrays
c
c itype     atom type number
c rtype     van der waals radius
c nua       number of unit vectors on sphere
c
c dimensioned one more in case input file is too long
c
      dimension itype(nxatm+1)
      dimension nua(nxatm)
c
c arrays for all atoms
c
c co        atomic coordinates
c ias       surface request number
c iat       atom itype
c molnum    molecule number
c rad       radius
c srs       some reentrant surface
c
c dimension arrays 1 more in case input file is too long
c
cxxx defined in common/connolly: co and molnum
cxxx      dimension co(3,maxatm+1)
cxxx  integer*2 ias(maxatm+1),molnum(maxatm+1)
      integer*2 ias(nxatm+1)
c     integer*2 ias(maxatm+1)
      integer*2 iat(nxatm+1)
c     integer*2 iat(maxatm+1)
      dimension rad(nxatm)
c     dimension rad(maxatm)
      logical srs(nxatm)
c     logical srs(maxatm)
c
c cube arrays
c
c ico       integer cube coordinates
c icuptr    pointer to next atom in cube
c comin     minimum atomic coordinates (cube corner)
c icube     pointer to first atom in list for cube
c scube     some atom with srn=2 in cube
c sscube    some atom with srn=2 in cube or adjoining cubes
c
      integer*2 ico(3,maxyon),icuptr(maxyon)
c     integer*2 ico(3,nxatm),icuptr(nxatm)
c     integer*2 ico(3,maxatm),icuptr(maxatm)
      dimension comin(3)
      integer*2 icube
      logical scube, sscube
      common /scrtch/ scube(maxcub,maxcub,maxcub), 
     +                sscube(maxcub,maxcub,maxcub),
     +                icube(maxcub,maxcub,maxcub)
c
c neighbor arrays
c
c inbr      atom number
c cnbr      coordinates
c rnbr      radius
c snbr      true if srn = 2 for neighbor
c mnbr      mutual neighbor of iatom and jatom
c molnbr    molecule number
c ernbr     expanded radius (rnbr + rp)
c disnbr    distance from neighbor to iatom
c lknbr     link to next farthest out neighbor
c itnl      temporary neighbor list (before sort)
c
      integer*2 inbr(maxnbr)
      dimension cnbr(3,maxnbr),rnbr(maxnbr)
      logical snbr(maxnbr),mnbr(maxnbr)
      integer*2 molnbr(maxnbr)
      dimension ernbr(maxnbr)
      dimension disnbr(maxnbr)
      integer*2 lknbr(maxnbr)
      integer*2 itnl(maxnbr)
c
c circle and sphere unit vector arrays
c
c up        unit vectors for probe
c ua        unit vectors for atom
c eva       extended vectors for atom
c circle    points on a circle
c
      common/scra/up(3,maxsph),circle(3,maxcir)
      dimension  ua(3,maxsph,nxatm),
     +           eva(3,maxsph,nxatm)
c
c ci        coordinates of atom i
c cj        coordinates of atom j
c ck        coordinates of atom k
c si        srn = 2 for atom i
c sj        srn = 2 for atom j
c sk        srn = 2 for atom k
c
      dimension ci(3), cj(3), ck(3)
      logical si,sj,sk,sns
c
c geometric construction vectors
c
c vij       vector from atom i to atom j
c uij       unit vector of vij
c q,t       two perpendicular vectors in the saddle plane
c cijk      center of circle of intersection of expanded sphere
c           of atom k with the saddle plane of atoms i and j
c vijk      vector from torus center to cijk
c uijk      unit vector of vijk
c bij       torus center
c aij       starting altitude vector from torus center
c bijk      base point for atoms i, j and k
c aijk      altitude vector for atoms i, j and k
c aijp      altitude vector to probe (rotated)
c a         altitude vector, general, used in orsr
c p         probe coordinates (general, used in orsr)
c pijp      center of probe placed tangent to atoms i and j
c pij       starting center of probe tangent to atoms i and j
c pijk      probe placed tangent to atoms i, j and k
c pipt      probe placed tangent to atom i
c vpi       vector from probe center to contact point with atom i
c vpj       vector from probe center to contact point with atom j
c vpk       vector from probe center to contact point with atom k
c vps0      starting arc points relative to probe center
c vbs0      starting arc points relative to torus center
c vbs       rotated arc point relative to torus center
c arca      area of arc point
c ayon      arc point on yon side of symmetry element
c vector    temporary vector storage
c
c
      dimension vij(3),uij(3),q(3),t(3),cijk(3),vijk(3),uijk(3)
      dimension bij(3),aij(3),bijk(3),aijk(3,2),aijp(3,2),a(3)
      dimension p(3),pij(3),pijp(3,2),pijk(3,2),pipt(3)
      dimension vpi(3),vpj(3),vpk(3)
      dimension vps0(3,maxarc),vbs0(3,maxarc,2),vbs(3),arca(maxarc)
      dimension vector(3)
c
c g         uij, q, t frame for torus
c h         rotation about x-axis
c ghgt      rotation about uij axis
c pow       powers of ghgt
c
c rotation matrices
c
       dimension g(3,3),h(3,3),ghgt(3,3),pow(3,3)
c
c s         surface points for probe
c torus     surface points for torus
c n1        atom surface point is on or closest to
c n2,n3     other atoms probe is touching
c yon       whether point lies on yon side of symmetry element
c
c reentrant probe record
c
      dimension s(3,maxppp),torus(maxppp)
      integer*2 n1(maxppp),n2(maxppp),n3(maxppp)
      logical yon(maxppp)
      integer*2 imol,np
c
c both      both probe positions free of collisions
c pair      this member of pair free from collisions
c yonprb    probe crosses symmetry element
c found     search flag
c
c logical variables
c
      logical both,pair(2),ayon(maxarc)
      logical found
      logical yonprb
c
c ascii     ascii output records
c long      long output records
c
      logical ascii,long
c
c orsr for non-symmetry-related probes
c the factor of three for the victim arrays
c is based upon experience
c
c py        center of yon probe
c ay        altitude vector of yon probe
c pv        center of victim probe
c av        altitude vector of victim probe
c ivic      list of probe numbers of victims
c ivicp     pointer to next victim with same hash
c molyon    molecule number of yon probe
c molvic    molecule number of victim probe
c eat       coordinates of eaters of probe
c
c     dimension py(3,nxatm),ay(3,nxatm)
      dimension py(3,maxyon),ay(3,maxyon)
      dimension pv(3,maxvic),av(3,maxvic)
      integer*2 ivic(maxvic),ivicp(maxvic)
c     integer*2 molyon(nxatm),molvic(maxvic)
      integer*2 molyon(maxyon),molvic(maxvic)
      dimension eat(3,maxeat)
c
c multiple-molecule variables
c
c bury      probe position is buried
c
      logical bury
c
c counters
c
c nias      number of atoms with given srn
c nshape    number of surface points with given shape
c nlost     number of surface points lost in non-symmetry orsr
c
      dimension nias(3)
      integer nshape(3),nlost(3)
c
c output variables
c
c iatnum   atom numbers
c ishape   surface shape number
c ib       buried surface flag
c
      integer*2 iatnum(3),ishape,ib
      REAL outco(3),outvec(3)
c
c logical functions
c
      logical collid
      logical buried
c
c big loop arrays
c
      parameter (maxit=50)
      dimension dd(maxit)
      logical oupper
      REAL nn(maxit)
      REAL nmid
      logical exist
      data small /1.0d-03/
      data nrepeat /5/
      data exist /.false./
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c initialize value of pi
c
      pi = acos(-1.0d0)
c
      ascii = .false.
      long =  .false.
c
caleko
c
c     If the user is working with surface density, this is
c     used, otherwise if the user is working with a maximum
c     number of surface points, a loop over the whole surface 
c     calculation is performed to get the
c     right surface density.
c
      do 10 i=1,maxit
        dd(i)=0.0d0
        nn(i)=0.0d0
10    continue

      if (isurdens) then 
        goto 555
      else
        nlow = 0
        nhigh = 0
        oupper = .false.
        iround=1  
        nmncon=0.8d0*nmaxcon
        if (iminpoint) nmncon=nmincon
cxxx
        nn(iround) = 0.0d0
        dd(iround) = 0.0d0
        imin = 1
        imax = 2
        iround = iround + 1
cxxx
        nmid=(nmaxcon+nmncon)/2.0d0
        dd(iround)=1.0d0/natom
        d=dd(iround)
        goto 333
 111    nn(iround-1)=ncon 
        dd(iround)= dd(iround-1)*nmid/nn(iround-1)
        d=dd(iround)
        goto 333
 222    nn(iround-1)=ncon
cxxx
c-----  guard against consecutive passes that 
c       deliver an equal number of points
cxxx    dd(iround)=
cxxx +    (abs(dd(iround-1)-dd(iround-2))+0.001)*nmid
cxxx +              /abs(nn(iround-1)-nn(iround-2)+1) 
          if (nn(iround-1) .eq. nn(imin)) nlow = nlow + 1
          if (nn(iround-1) .eq. nn(imax)) nhigh= nhigh+ 1
           if (.not. oupper) then
             if (nn(iround-1) .gt. nn(imax)) then
               if (nn(imax) .lt. nmncon) then
                  imin = imax
               endif
               imax = iround - 1
             endif
             if (nn(iround-1) .gt. nmaxcon) then
               oupper = .true.
             endif
           else  
             if (nn(iround-1) .gt. nmaxcon) then
               if (nn(iround-1) .le. nn(imax)) then
                 imax = iround - 1
               else if (dd(iround-1) .lt. dd(imax)) then
                 imax = iround - 1
               endif
             endif
           endif
           if (nn(iround-1) .lt. nmncon) then
             if (nn(iround-1) .ge. nn(imin)) then
               imin = iround - 1
             else if (dd(iround-1) .gt. dd(imin)) then
               imin = iround - 1
             endif
           endif
c       if (abs(nn(iround-1)-nn(iround-2)) .lt. small) then
c         if (nn(iround-1) .gt. nmaxcon) then
c           dd(iround) = dd(iround-1)*0.99
c         else if(nn(iround-1) .lt. nmncon) then
c           dd(iround) = dd(iround-1)*1.01
          if ((nlow .gt. nrepeat) .or.
     1        (nhigh .gt. nrepeat)) then
            write(6,*) ' Connolly surface not refined further'
            if (nn(imax) .lt. nmaxcon*1.1d0) then
              d=dd(imax)
              ncon=nn(imax)
            else
              d=dd(imin)
              ncon=nn(imin)
              if (d .lt. 1.d-10) then
                d = dd(imax)
                ncon=nn(imax)
              endif
            endif
            goto 448
          endif
c       else
          dd(iround) = (dd(imax) + dd(imin)) /2.0d0
          if (nn(imax) .le. nmncon) 
     1     dd(iround) = dd(iround-1)*2.0d0
c       endif
        d=dd(iround)   
 333    continue
      endif
 555  continue

cxxx
      ncon=0
c
c input value checking
c check for negative probe radius
c
      if (rp .lt. 0.0d0) call error(120,0,rp)
c
c initialize srn counters and coordinate minima
c
      do 350 k = 1,3
         nias(k) = 0
350   continue
      do 400 k = 1,3
         comin(k) = 1000000.0d0
400   continue

      radmax = 0.0d0
c
c     checking loop over atoms 
c
      do 750 n=1,natom 
c
caleko
c
c     ibury is standard set to 0 and ias to 2
c
      ibury=0
      ias(n)=2
c
c atom radii must be zero or positive
c
      if (rtype(n) .lt. 0.0d0) call error(140,n,rtype(n))
c
c check for new maximum radius
c
c     if (rtype(n) .gt. radmax) radmax = rtype(n)
      if (rtype(n)+solrad .gt. radmax) 
     1   radmax = rtype(n)+solrad
cxxx
      rad(n) = rtype(n) + solrad
cxxx
c
caleko
c
c Increment atom radii with the solvent radius to get the 
c surface at a reasonable distance
c
      rtppls(n)=rtype(n)+solrad

c
c number of unit vectors depends on sphere area and input density
c
      nua(n) = (4.0d0*pi*rtppls(n)**2) * d
c
c decrease to array size if too large
c
      if (nua(n) .gt. maxsph) nua(n) = maxsph
      if (nua(n) .lt. 1) nua(n) = 1
c
c create unit vector arrays
c
      call genun(ua(1,1,n),nua(n))
c
c compute extended vectors for later probe placement
c
      do 250 isph = 1,nua(n)
         do 200 k = 1,3
            eva(k,isph,n) = (rtppls(n) + rp) * ua(k,isph,n)
200      continue
250   continue

c
c calculate width of cube from maximum atom radius and probe radius
c
      width = 2.0d0 * (radmax + rp)
c
c check for atom overflow
c
      if (n .gt. maxatm) call error(150,n,0.0d0)
c
c check surface request number
c
      if (ias(n) .ne. 2) call error(160,n,0.0d0)
c
c check for new coordinate minima
c
      do 550 k = 1,3
         if (co(k,n) .lt. comin(k)) comin(k) = co(k,n)
550   continue
c
c increment counters for each srn type
c
      do 600 k = 1,3
         if (ias(n) .eq. k-1) nias(k) = nias(k) + 1
600   continue
c
c   End of input loop
c

750   continue

c
c     set up cube arrays
c     first the integer coordinate arrays
c
      do 850 i = 1,natom
         do 800 k = 1,3
            ico(k,i) = (co(k,i)-comin(k))/width + 1
            if (ico(k,i) .lt. 1) call caserr
     +         ('cube coordinate too small')
            if (ico(k,i) .gt. maxcub)call caserr
     +         ('cube coordinate too large')
 800      continue
 850   continue
c
c initialize head pointer and srn=2 arrays
c
      do 1000 k = 1,maxcub
         do 950 j = 1,maxcub
            do 900 i = 1,maxcub
               icube(i,j,k) = 0
               scube(i,j,k) = .false.
               sscube(i,j,k) = .false.
 900         continue
 950      continue
 1000  continue
c
c initialize linked list pointers
c
      do 1050 i = 1,natom
         icuptr(i) = 0
 1050  continue
c
c set up head and later pointers for each atom
c
      do 1250 iatom = 1,natom
c
c skip atoms with surface request numbers of zero
c
         if (ias(iatom) .eq. 0) go to 1250
         i = ico(1,iatom)
         j = ico(2,iatom)
         k = ico(3,iatom)
         if (icube(i,j,k) .le. 0) then
c
c     first atom in this cube
c
            icube(i,j,k) = iatom
         else
c
c     add to end of linked list
c
            iptr = icube(i,j,k)
1100        continue
c
c check for duplicate coordinates
c
            d2 = dist2(co(1,iatom),co(1,iptr))
            if (d2 .le. 0.0d0) then
               ias(iatom) = 0
cxxx           write (6,1150) iatom,iptr
               if (ibemout .ge. 2) write (6,1150) iatom
cxxx           format(1x,'atom',i5,' dropped (same co as ',i5,')')
1150           format(1x,'atom',i5,
     1       ' dropped (overlapping coordinates)')
               go to 1250
            end if
            if (icuptr(iptr) .le. 0) go to 1200
c
c move on down the list
c
            iptr = icuptr(iptr)
            go to 1100
1200        continue
c
c store atom number
c
            icuptr(iptr) = iatom
         end if
c
c check for surfaced atom
c
         if (ias(iatom) .eq. 2) scube(i,j,k) = .true.
1250  continue
c
c check for 3 x 3 x 3 with some srn = 2
c
      do 1550 k = 1,maxcub
         do 1500 j = 1,maxcub
            do 1450 i = 1,maxcub
               if (icube(i,j,k) .eq. 0) go to 1450
c
c check whether this cube or any adjacent cube has srn = 2
c
               do 1400 k1 = k-1,k+1
                  if (k1 .lt. 1 .or. k1 .gt. maxcub) go to 1400
                  do 1350 j1 = j-1,j+1
                     if (j1 .lt. 1 .or. j1 .gt. maxcub) go to 1350
                     do 1300 i1 = i-1,i+1
                        if (i1 .lt. 1 .or. i1 .gt. maxcub) go to 1300
                        if (scube(i1,j1,k1)) sscube(i,j,k) = .true.
1300                 continue
1350              continue
1400           continue
1450        continue
1500     continue
1550  continue
c
c initialization
c maximum number of neighbors any atom has
c
      maxnb = 0
c
c numbers of surface points
c
      do 1600 k = 1,3
         nshape(k) = 0
         nlost(k) = 0
1600  continue
c
c number of yon probes
c
      ny = 0
c
c contact and reentrant areas
c
      areac = 0.0d0
      arear = 0.0d0
c
c write out messages
c
      if (ibemout .ge. 2) then
      write (6,1650) natom
1650  format(1x,i5,1x,'atoms')
      write (6,1700) nias(1)
1700  format(1x,i5,1x,'omitted')
      write (6,1750) nias(2)
1750  format(1x,i5,1x,'collision only')
      write (6,1800) nias(3)
1800  format(1x,i5,1x,'surface')
      write (6,1850) d,rp
1850  format(1x,'surface point density = ',f10.5,5x,
     1 'probe radius = ',f10.5)
      if (ibury .eq. 1) write (6,1900)
1900  format(1x,'buried surface only')
      if (ibury .eq. 2) write (6,1950)
1950  format(1x,'buried surface flagged')
      if (.not. ascii) write (6,2000)
2000  format(1x,'binary contact and reentrant files')
      if (.not. long) write (6,2050)
2050  format(1x,'short output records')
      endif
c
c stop if density is not positive
c
      if (d .le. 0.0d0) call caserr
     +   ('non-positive density')
c
c skip probe and circle setup if van der waals surface
c
      if (rp .eq. 0.0d0) go to 2150
c
c set up probe sphere and circle
c
      nup = (4.0d0 * pi * rp ** 2) * d
      if (nup .lt. 1) nup = 1
       if (nup .gt. maxsph) nup = maxsph
      call genun(up,nup)
      ncirc = (2.0d0 * pi * rp) * sqrt(d)
      if (ncirc .lt. 1) ncirc = 1
      if (ncirc .gt. maxcir) ncirc = maxcir
      do 2100 i = 1, ncirc
         fi = (2.0d0 * pi * (i-1))/ncirc
         circle(1,i) = rp * cos(fi)
         circle(2,i) = rp * sin(fi)
         circle(3,i) = 0.0d0
2100  continue
c
      inquire(file='con4',exist=exist)
      if (exist) then
        open(4,form='unformatted',file='con4',status='unknown')
        do i = 1, maxyon
          backspace(4)
        enddo
      else
        open(4,form='unformatted',file='con4',status='new')
      endif
c open before file for writing
cxxx 5/91
c     open(4,form='unformatted',status='new')
c     open(4,form='unformatted',file='con4')
c     open(4,form='unformatted',file='con4',status='delete')
c     open(4,form='unformatted',file='con4',status='new',err=2125)
c     goto 2126
cxxx
c2125 continue
c2125 write(6,*) " Remove the file con4 from the working directory"
c     stop
c2126 continue
      rewind(4)
c
c skip to here if no reentrant surface will be calculated
c
2150  continue
c
c open contact file for writing
c
caleko
c  Commented out
c
c      if (ascii) then
cxxx 5/91
c        open(9,status='new')
c         open(9,form='formatted',file='csurf')
c      else
c
c        open(9,form='unformatted',status='new')
c
c         open(9,form='unformatted',file='csurf')
cxxx
c      end if
c      open(7,form='unformatted',file='con7')
c      rewind(7)
c
c initialize some reentrant surface to false for each atom
c
      do 2200 iatom = 1,natom
         srs(iatom) = .false.
2200  continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c big loop for each atom
c
      do 5650 iatom = 1, natom
c
c skip ignored atoms
c
         if (ias(iatom) .eq. 0) go to 5650
c
         ici = ico(1,iatom)
         icj = ico(2,iatom)
         ick = ico(3,iatom)
c
c skip iatom if its cube and adjoining cubes contain only blockers
c
         if (.not. sscube(ici,icj,ick)) go to 5650
c
c transfer values from large arrays to iatom variables
c
         ri = rtype(iatom)+solrad
         si = ias(iatom) .eq. 2
         do 2250 k = 1,3
            ci(k) = co(k,iatom)
2250     continue
         imol = molnum(iatom)
c
c gather the neighboring atoms of iatom
c initialize number of neighbors, and number of neighbors in the
c same molecule as atom i
c
         nnbr = 0
         nimol = 0
c
c initialize srn = 2 for some neighbor to false
c
         sns = .false.
c
c save a little time for distance check
c
         sumi = 2.0d0 * rp + ri
c
c check iatom cube and adjacent cubes for neighboring atoms
c
         do 2550 jck = ick-1,ick+1
            if (jck .lt. 1 .or. jck .gt. maxcub) go to 2550
            do 2500 jcj = icj-1,icj+1
              if (jcj .lt. 1 .or. jcj .gt. maxcub) go to 2500
                do 2450 jci = ici-1,ici+1
                  if (jci .lt. 1 .or. jci .gt. maxcub) go to 2450
                  jatom = icube(jci,jcj,jck)
2300              continue
c
c check for end of linked list for this cube
c
                  if (jatom .le. 0) go to 2400
c
c distance check
c
                  sum = sumi + rad(jatom)
                  vect1 = abs(co(1,jatom) - ci(1))
                  if (vect1 .ge. sum) go to 2350
                  vect2 = abs(co(2,jatom) - ci(2))
                  if (vect2 .ge. sum) go to 2350
                  vect3 = abs(co(3,jatom) - ci(3))
                  if (vect3 .ge. sum) go to 2350
                  d2 = vect1 ** 2 + vect2 ** 2 + vect3 ** 2
                  if (d2 .ge. sum ** 2) go to 2350
c
c iatom is not its own neighbor
c
                  if (iatom .eq. jatom) go to 2350
c
c we have a new neighbor
c
                  nnbr = nnbr + 1
c
c check for neighbor overflow
c
                 if (nnbr .gt. maxnbr) call error(210,nnbr,0.0d0)
c
c save atom number in temporary array

                  itnl(nnbr) = jatom
c
c check whether surfaced neighbor in same molecule
c
                  if (ias(jatom) .eq. 2 .and. molnum(jatom) .eq. imol)
     1  sns = .true.
c
c count the number of atoms in the same molecule as iatom
c
                  if (imol .eq. molnum(jatom)) nimol = nimol + 1
2350              continue
c
c get number of next atom in cube
c
                  jatom = icuptr(jatom)
                  go to 2300
2400              continue
2450           continue
2500        continue
2550     continue
c
c keep track of maximum number of neighbors
c for array-dimensioning purposes
c
         if (nnbr .gt. maxnb) maxnb = nnbr
c
c no surface for atom i if buried only flag set and
c there are no neighbors from the other molecules
c
         if (ibury .eq. 1 .and. nimol .eq. nnbr) go to 5650
c
c no surface if iatom and all neighbor atoms
c in the same molecule have surface request numbers < 2
c
         if (.not. si .and. .not. sns) go to 5650
c
c set up neighbors arrays with jatom in increasing order

c initialize minimum neighbor atom number
c
         jmold = 0
         do 2700 iuse = 1,nnbr
            jmin = natom + 1
            do 2600 jnbr = 1,nnbr
c
c don't use ones already sorted
c
               if (itnl(jnbr) .le. jmold) go to 2600
               if (itnl(jnbr) .lt. jmin) then
                  jmin = itnl(jnbr)
                  jminbr = jnbr
               end if
2600        continue
            jmold = jmin
            jnbr = jminbr
            jatom = itnl(jnbr)
c
c transfer atom number, coordinates, radius, surface request number,
c molecule number, expanded radius, distance from iatom
c
            inbr(iuse) = jatom
            do 2650 k = 1,3
               cnbr(k,iuse) = co(k,jatom)
2650        continue
            rnbr(iuse) = rad(jatom)
            snbr(iuse) = ias(jatom) .eq. 2
            molnbr(iuse) = molnum(jatom)
            ernbr(iuse) = rnbr(iuse) + rp
            disnbr(iuse) = dist2(ci,cnbr(1,iuse))
c
c initialize link to next farthest out neighbor
c
            lknbr(iuse) = 0
2700     continue
c
c set up a linked list of neighbors in order of
c increasing distance from iatom
c initialize pointer to first neighbor to 0
c
         lkf = 0
c
c look for neighbor in same molecule
c we want only atoms in same molecule for collision check
c
         do 2750 l = 1,nnbr
            if (imol .ne. molnbr(l)) go to 2750
            lkf = l
            go to 2800
2750     continue
         if (lkf .eq. 0) go to 3000
2800     continue
c
c put remaining neighbors in linked list at proper position
c
         do 2950 l = lkf+1,nnbr
            if (imol .ne. molnbr(l)) go to 2950
            l1 = 0
            l2 = lkf
2850        continue
            if (disnbr(l) .lt. disnbr(l2)) go to 2900
            l1 = l2
            l2 = lknbr(l2)
            if (l2 .ne. 0) go to 2850
2900        continue
c
c add to list
c
             if (l1 .eq. 0) then
               lkf = l
               lknbr(l) = l2
            else
               lknbr(l1) = l
               lknbr(l) = l2
             end if
2950     continue
3000     continue

c no reentrant surface will be calculated if we are
c calculating the van der waals surface
c instead of the molecular surface
c
        if (rp .eq. 0.0d0) go to 5200

c no reentrant surface if iatom has no neighbors
c
         if (nimol .le. 0) go to 5200
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c medium loop for each neighbor of iatom
c
         do 5150 jnbr = 1,nnbr
            jatom = inbr(jnbr)
c
c each pair of atoms is considered only once
c
            if (jatom .le. iatom) go to 5150
c
c each molecule gets a separate surface
c
            if (imol .ne. molnbr(jnbr)) go to 5150
c
c tranfer from neighbor arrays to jatom variables
c
            rj = rnbr(jnbr)
            sj = snbr(jnbr)
            do 3050 k = 1,3
               cj(k) = cnbr(k,jnbr)
3050        continue
c
c here follow geometric calculations of points, vectors and
c distances used for probe placement in both saddle and
c concave reentrant surface generation
c
c calculate the intersection
c of the expanded spheres of iatom and jatom
c this circle is called the saddle circle
c the plane it lies in is called the saddle plane
c
            do 3100 k = 1,3
                vij(k) = cj(k) - ci(k)
3100        continue
c
c create an orthonormal frame
c with uij pointing along the inter-atomic axis
c and q and t defining the saddle plane
c
            if (anrm(vij) .le. 0.0d0) then
               write (6,3150) iatom,jatom
3150           format(1x,'atoms',2i5,' have the same center')
               go to 5150
            end if
            call vnorm(vij,uij)
            call vperp(uij,q)
            call croscon(uij,q,t)
c
c calculate the saddle circle center and radius
c
            dij = anrm(vij)
            f = 0.5d0*(1.0d0+((ri+rp)**2-(rj+rp)**2)/dij**2)
c
c base point
c
            do 3200 k = 1,3
               bij(k) = ci(k) + f * vij(k)
3200        continue
            f1 = (ri+rj+2*rp)**2 - dij**2
c
c skip to bottom of middle loop if atoms are too far apart
c
            if (f1 .le. 0.0d0) go to 5150
            f2 = dij**2 - (ri-rj)**2
c
c skip to bottom of middle loop if one atom inside the other
c
            if (f2 .le. 0.0d0) go to 5150
c
c height (radius of saddle circle)
c
            hij = sqrt(f1*f2) / (2*dij)
c
c a starting altitude
c
            do 3250 k = 1,3
               aij(k) = hij * q(k)
3250        continue
c
c concave reentrant surface
c
c gather mutual neighbors of iatom and jatom
c
            mutual = 0
            do 3300 knbr = 1, nnbr
               d2 = dist2(cj,cnbr(1,knbr))
               mnbr(knbr) = (d2 .lt. (2.d0*rp+rj+rnbr(knbr))**2)
     1 .and. (knbr .ne. jnbr)
               if (mnbr(knbr)) mutual=mutual+1
3300        continue
c
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c inner loop for each mutual neighbor of iatom and jatom
c
            ishape = 3
            do 4200 knbr = 1,nnbr
               if (.not. mnbr(knbr)) go to 4200
               katom = inbr(knbr)
c
c iatom < jatom < katom
c
               if (katom .le. jatom) go to 4200
               sk = snbr(knbr)
c
c skip neighbor if all three atom not marked to be surfaced
c
               if (.not. (si .or. sj .or. sk)) go to 4200
c
c each molecule gets a separate surface
c
               if (imol .ne. molnbr(knbr)) go to 4200
c
c tranfer from neighbor array to katom variables
c
               rk = rnbr(knbr)
               do 3350 k = 1,3
                  ck(k) = cnbr(k,knbr)
3350           continue
c
c calculate intersection of expanded sphere of katom
c with saddle plane. we will call this the katom circle.
c
c projection of vector,
c from katom to a point on the saddle plane,
c onto iatom-jatom axis,
c in order to get distance katom is from saddle plane
c
               dk = uij(1) * (bij(1)-ck(1)) + 
     +              uij(2) * (bij(2)-ck(2)) +
     +              uij(3) * (bij(3)-ck(3))
c
c calculate radius of katom circle
c
               rijk = (rk+rp) ** 2 - dk ** 2
c
c skip concave calculation if no intersection
c
               if (rijk .le. 0.0d0) go to 4200
               rijk = sqrt(rijk)
c
c calculate center of katom circle
c
               do 3400 k = 1,3
                  cijk(k) = ck(k) + dk * uij(k)
3400           continue
c
c calculate intersection of the katom circle with the saddle circle
c
               do 3450 k = 1,3
                  vijk(k) = cijk(k) - bij(k)
3450           continue
               dijk = anrm(vijk)
               if (dijk .le. 0.0d0) then
                  write (6,3500) iatom,jatom,katom
3500              format(1x,'atoms',3i5,' have concentric circles')
                  go to 4200
               end if
               f = 0.5d0 * (1.0d0+(hij**2-rijk**2)/dijk**2)
c
c base point bijk is on symmetry plane and saddle plane
c
               do 3550 k = 1,3
                  bijk(k) = bij(k) + f * vijk(k)
3550           continue
               f1 = (hij+rijk)**2 - dijk**2
c
c skip to bottom of inner loop if katom too far away
c
               if (f1 .le. 0.0d0) go to 4200
               f2 = dijk**2 - (hij-rijk)**2
c
c skip to bottom of inner loop if katom circle inside saddle circle
c or vice-versa
c
               if (f2 .le. 0.0d0) go to 4200
               hijk = sqrt(f1*f2) / (2*dijk)
               call vnorm(vijk,uijk)
c
c uij and uijk lie in the symmetry plane passing through the atoms
c so their cross product is perpendicular to this plane
c
               call croscon(uij,uijk,aijk)
c
c two altitudes
c
               do 3600 k = 1,3
                  aijk(k,1) = hijk * aijk(k,1)
                  aijk(k,2) = - aijk(k,1)
3600           continue
c
c probe placement at ends of altitude vectors
c
               do 3700 ip = 1,2
                  do 3650 k = 1,3
                     pijk(k,ip) = bijk(k) + aijk(k,ip)
3650              continue
c
c collision check with mutual neighbors
c
                  pair(ip) = .not. collid(pijk(1,ip),rp,cnbr,ernbr,mnbr,
     1  nnbr,maxnbr,ishape,jnbr,knbr,molnbr,imol,lkf,lknbr)
3700           continue
c
c if neither probe position is allowed, skip to bottom of inner loop
c
               if (.not. pair(1) .and. .not. pair(2)) go to 4200
               both = pair(1) .and. pair(2)
c
c some reentrant surface for all three atoms
c
               srs(iatom) = .true.
               srs(jatom) = .true.
               srs(katom) = .true.
c
c generate surface points
c
               area = (4d0 * 3.14159d0 * rp ** 2)/nup
               do 4150 ip = 1,2
                  if (.not. pair(ip)) go to 4150
c
c give it some kind of value, in case we don't call buried
c
                  bury = .false.
c
c only call buried if we care what the answer is
c
                  if (ibury .gt. 0) bury = buried(pijk(1,ip),rp,cnbr,rnb
     1 r,mnbr,nnbr,maxnbr,ishape,jnbr,knbr,molnbr,imol)
c
c skip if not buried and buried surface only flag set
c
                  if (ibury .eq. 1 .and. .not. bury) go to 4150
c
c determine whether probe has surface on far side of plane
c
                  yonprb = hijk .lt. rp .and. .not. both
c
c calculate vectors defining spherical triangle
c the vectors are given the probe radius as a length
c only for the purpose of making the geometry more clear
c
                  do 3750 k = 1,3
                     vpi(k) = (ci(k) - pijk(k,ip)) * rp / (ri + rp)
                     vpj(k) = (cj(k) - pijk(k,ip)) * rp / (rj + rp)
                     vpk(k) = (ck(k) - pijk(k,ip)) * rp / (rk + rp)
3750              continue
                  sign = det(vpi,vpj,vpk)
c
c initialize number of surface points written
c
                  np = 1
c
c gather points on probe sphere lying within triangle
c
                  do 4000 i = 1,nup
c
c if the unit vector is pointing away from the symmetry plane
c the surface point cannot lie within the inward-facing triangle
c
                     if (ddot(3,up(1,i),1,aijk(1,ip),1) .gt. 0.0d0) 
     +                    go to 4000
                     if (sign * det(up(1,i),vpj,vpk) .lt. 0.0d0)
     1                 go to 4000
                     if (sign * det(vpi,up(1,i),vpk) .lt. 0.0d0)
     1                 go to 4000
                     if (sign * det(vpi,vpj,up(1,i)) .lt. 0.0d0)
     1                 go to 4000
                     if (np .gt. maxppp) call error(320,np,0.0d0)
c
c calculated whether point is on yon side of plane
c
                     yon(np) = aijk(1,ip) * (aijk(1,ip) + up(1,i)) +
     1 aijk(2,ip) * (aijk(2,ip) + up(2,i)) +
     2 aijk(3,ip) * (aijk(3,ip) + up(3,i)) .lt. 0.0d0
c
c overlapping reentrant surface removal
c for symmetry-related probe positions
c
                     if (yon(np) .and. both) go to 4000
c
c calculate coordinates of surface point
c
                     do 3800 k = 1,3
                        s(k,np) = pijk(k,ip) + up(k,i) * rp
3800                 continue
c
c find the closest atom and put the three atom numbers
c in the proper order
c n1 is closest, n2 < n3
c
                     dsi = dist3(s(1,np),ci) - ri
                     dsj = dist3(s(1,np),cj) - rj
                     dsk = dist3(s(1,np),ck) - rk
                     if (dsi .le. dsj .and. dsi .le. dsk) go to 3850
                     if (dsj .le. dsi .and. dsj .le. dsk) go to 3900
                     if (.not. sk) go to 4000
                     n1(np) = katom
                     n2(np) = iatom
                     n3(np) = jatom
                     go to 3950
3850                 continue
                     if (.not. si) go to 4000
                     n1(np) = iatom
                     n2(np) = jatom
                     n3(np) = katom
                     go to 3950
3900                 continue
                     if (.not. sj) go to 4000
                     n1(np) = jatom
                     n2(np) = iatom
                     n3(np) = katom
3950                 continue
                     np = np + 1
c
c end of nup loop
c
4000              continue
                  np = np - 1
c
c skip the write if no points
c
                  if (np .le. 0) go to 4150
c
c write the molecule number, shape, number of points,
c probe position and
c the vector from the base to the probe center
c
                  write (4) imol,ishape,np,(pijk(k,ip),k=1,3),
     1 (aijk(k,ip),k=1,3),yonprb,bury
c
c save probe in yon probe arrays
c
                  if (yonprb) then
c
c check for overflow
c
                     if (ny .ge. maxyon) call error(720,ny,0.0d0)
                     ny = ny + 1
                     molyon(ny) = imol
                     do 4050 k = 1,3
                        py(k,ny) = pijk(k,ip)
                        ay(k,ny) = aijk(k,ip)
4050                 continue
                  end if
c
c write surface points for this probe position
c
                  do 4100 i = 1,np
                     write (4) n1(i),n2(i),n3(i),(s(k,i),k=1,3),
     1               area,yon(i)
4100              continue
c
c end of ip loop
c
4150           continue
c
c end of concave reentrant loop
c
4200        continue
c
c saddle-shaped reentrant
c
            ishape = 2
c
c check for neither atom to be surfaces
c
            if (.not. (si .or. sj)) go to 5150
c
c special check for buried tori
c
c if both atoms are marked to be surface,
c but neither atom has any reentrant surface so far
c (after triangles with all katoms have been checked)
c and if there is some mutual neighbor in the same molecule
c close enough so that the torus cannot be free,
c then we know that this must be a buried torus
c
            if (si .and. sj .and. .not. srs(iatom) .and.
     1        .not. srs(jatom) .and. mutual .gt. 0) then
               do 4250 knbr = 1,nnbr
                  if (.not. mnbr(knbr)) go to 4250
                  if (imol .ne. molnbr(knbr)) go to 4250
                  d2 = dist2(bij,cnbr(1,knbr))
                  rk2 = ernbr(knbr) ** 2 - hij ** 2
                  if (d2 .lt. rk2) go to 5150
4250           continue
            end if
c
c calculate number of rotations of probe pair,
c rotation angle and rotation matrix
c
            rij = ri/(ri + rp) + rj/(rj + rp)
            avh = (abs(hij-rp) + hij * rij)/3.0d0
            nrot = sqrt(d) * pi * avh
            if (nrot .lt. 1) nrot = 1
            angle = 3.14159d0/nrot
c
c set up rotation matrix around x-axis
c
            call imatx(h)
            h(2,2) = cos(angle)
            h(3,3) = h(2,2)
            h(3,2) = sin(angle)
            h(2,3) = - h(3,2)
c
c calculate matrix to rotate x-axis onto iatom-jatom axis
c
            do 4300 k = 1,3
               g(k,1) = uij(k)
               g(k,2) = q(k)
               g(k,3) = t(k)
4300        continue
c
c make the probe pair rotation matrix be about the iatom-jatom axis
c
            call conj(h,g,ghgt)
c
c arc generation
c
            do 4350 k = 1,3
               pij(k) = bij(k) + aij(k)
               vpi(k) = (ci(k) - pij(k)) * rp / (ri + rp)
               vpj(k) = (cj(k) - pij(k)) * rp / (rj + rp)
4350        continue
c
c rotate circle onto iatom-jatom-probe plane
c and select points between probe-iatom and
c probe-jatom vector to form the arc
c
            narc = 1
            do 4500 i = 1,ncirc
              if (narc .gt. maxarc) call error(440,narc,0.0d0)
c
c rotation
c
               call multv(circle(1,i),g,vps0(1,narc))
c
c if the vector is pointing away from the symmetry line
c the surface point cannot lie on the inward-facing arc
c
               if (ddot(3,vps0(1,narc),1,aij,1) .gt. 0.0d0) go to 4500
               call croscon(vpi,vps0(1,narc),vector)
               if (ddot(3,g(1,3),1,vector,1) .lt. 0.0d0) go to 4500
               call croscon(vps0(1,narc),vpj,vector)
               if (ddot(3,g(1,3),1,vector,1) .lt. 0.0d0) go to 4500
c
c make arc point vectors originate with saddle circle center bij
c rather than probe center because they will be
c rotated around the iatom-jatom axis
c
               do 4400 k = 1,3
                  vbs0(k,narc,1) = vps0(k,narc) + aij(k)
4400           continue
c
c invert arc through line of symmetry
c
               duij = ddot(3,uij,1,vbs0(1,narc,1),1)
               do 4450 k = 1,3
                  vbs0(k,narc,2) = - vbs0(k,narc,1) + 2 * duij * uij(k)
4450           continue
c
c check whether the arc point crosses the iatom-jatom axis
c and calculate the area associated with the point
c
               ht = ddot(3,aij,1,vbs0(1,narc,1),1)/hij
               ayon(narc) = ht .lt. 0.d00
               arca(narc) = (2d0*3.14159d0**2*rp*abs(ht))/(ncirc*nrot)
               narc = narc + 1
4500        continue
            narc = narc - 1
c
c initialize power matrix to identity
c
            call imatx(pow)
c
c set knbr to zero for collision and buried checks
c
            knbr = 0
c
c rotate the probe pair around the pair of atoms
c
            do 5100 irot = 1,nrot
c
c multiply altitude vector by power matrix
c
               call multv(aij,pow,aijp(1,1))
c
c set up opposing altitude
c
               do 4550 k = 1,3
                  aijp(k,2) = - aijp(k,1)
4550           continue
c
c set up probe sphere positions
c
               do 4650 ip = 1,2
                  do 4600 k = 1,3
                     pijp(k,ip) = bij(k) + aijp(k,ip)
4600              continue
c
c check for collisions with neighboring atoms
c
                  pair(ip) = .not. collid(pijp(1,ip),rp,cnbr,ernbr,mnbr,
     1 nnbr,maxnbr,ishape,jnbr,knbr,molnbr,imol,lkf,lknbr)
4650           continue
c
c no surface generation if neither probe position is allowed
c
               if (.not. pair(1) .and. .not. pair(2)) go to 5050
               both = pair(1) .and. pair(2)
c
c some reentrant surface for both atoms
c
               srs(iatom) = .true.
               srs(jatom) = .true.
c
c skip to bottom of middle loop if iatom and jatom
c are close enough and the surface point density is
c low enough so that the arc has no points
c
               if (narc .le. 0) go to 5050
c
c surface generation
c
               do 5000 ip = 1,2
                  if (.not. pair(ip)) go to 5000
c
c set default value for bury
c
                  bury = .false.
c
c don't check for probe collisions against other molecules
c unless we need to
c
                  if (ibury .gt. 0) bury = buried(pijp(1,ip),rp,cnbr,
     1  rnbr,mnbr,nnbr,maxnbr,ishape,jnbr,knbr,molnbr,imol)
c
c skip if not buried and buried surface only flag set
c
                  if (ibury .eq. 1 .and. .not. bury) go to 5000
c
c determine whether probe has surface on far side of line
c
                  yonprb = hij .lt. rp .and. .not. both
                  np = 1
c
c the saddle-shaped reentrant surface points come from the arc
c
                  do 4850 i = 1,narc
c
c overlapping reentrant surface removal
c for symmetry-related probe positions
c
                     if (both .and. ayon(i)) go to 4850
                     if (np .gt. maxppp) call error(480,np,0.0d0)
c
c rotate the arc from the xy plane onto the iatom-jatom-probe plane
c
                     call multv(vbs0(1,i,ip),pow,vbs)
c
c make coordinates relative to origin
c
                     do 4700 k = 1,3
                        s(k,np) = bij(k) + vbs(k)
4700                 continue
c
c find the closest atom and set up the atom numbers for the point
c
                     dsi = dist3(s(1,np),ci) - ri
                     dsj = dist3(s(1,np),cj) - rj
                     if (dsi .le. dsj) go to 4750
                     if (.not. sj) go to 4850
                     n1(np) = jatom
                     n2(np) = iatom
                     n3(np) = 0
                     go to 4800
4750                 continue
                     if (.not. si) go to 4850
                     n1(np) = iatom
                     n2(np) = jatom
                     n3(np) = 0
4800                 continue
c
c we've got a surface point
c
                     yon(np) = ayon(i)
                     torus(np) = arca(i)
                     np = np + 1
c
c end of arc point loop
c
4850              continue
                  np = np - 1
                  if (np .le. 0) go to 5000
c
c write the molecule number, shape,number of points,
c probe position and the vector from the base to the probe center
c
                  write (4) imol,ishape,np,(pijp(k,ip),k=1,3),
     1 (aijp(k,ip),k=1,3),yonprb,bury
                  if (yonprb) then
c
c save probe in yon probe arrays
c check for overflow
c
                     if (ny .ge. maxyon) call error(720,ny,0.0d0)
                     ny = ny + 1
                     molyon(ny) = imol
                     do 4900 k = 1,3
                        py(k,ny) = pijp(k,ip)
                        ay(k,ny) = aijp(k,ip)
4900                 continue
                  end if
c
c write surface points for this probe position
c
                  do 4950 i = 1,np
                     write (4) n1(i),n2(i),n3(i),(s(k,i),k=1,3),
     1 torus(i),yon(i)
c
c end of arc point loop
c
4950              continue
c
c end of probe pair loop
c
5000           continue
c
c skip to here if both probe positions disallowed or no arc points
c
5050           continue
c
c calculate new power matrix
c
               call cat(pow,ghgt)
c
c end of rotation loop
c
5100        continue
c
c end of neighbor loop
c
5150     continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c skip to here if van der waals surface calculation
c
5200     continue
c
c contact surface
c
         ishape = 1
c
c skip atom i if marked no surface requested
c
         if (.not. si) go to 5650
c
c if we are not calculating buried surface
c and the probe radius is greater than zero
c and iatom has at least one neighbor, but no reentrant surface,
c then iatom must be completely inaccessible to the probe
c
         if (.not. bury .and. rp .gt. 0.0d0.and.
     1 nimol .gt. 0 .and. .not. srs(iatom)) go to 5650
c
c find the index into the atom type arrays for iatom
c
caleko
c  Commented out, try to use only atom arrays
c
c        do 5250 idx = 1,ntype
c            if (itype(idx) .eq. iat(iatom)) go to 5300
c5250     continue
c         stop 'logic error in ms regarding atom types'
c5300     continue
         area = (4.0d0 * pi * ri ** 2) / nua(iatom)
c
c set jnbr, knbr to zero for collision, buried checks
c
         jnbr = 0
         knbr = 0
c
c contact probe placement loop
c
         do 5600 i = 1,nua(iatom)
c
c set up probe coordinates
c
            do 5350 k = 1,3
               pipt(k) = ci(k) + eva(k,i,iatom)
5350        continue
c
c check for collision with neighboring atoms
c
            if (collid(pipt,rp,cnbr,ernbr,mnbr,nnbr,maxnbr,ishape,
     1 jnbr,knbr,molnbr,imol,lkf,lknbr)) go to 5600
c
c go write it out if we don't care about buried surface
c
            if (ibury .eq. 0) then
               ib = 0
               go to 5400
            end if
            bury = buried(pipt,rp,cnbr,rnbr,mnbr,nnbr,maxnbr,ishape,
     1 jnbr,knbr,molnbr,imol)
            if (ibury .eq. 1 .and. .not. bury) go to 5600
            if (bury) then
               ib = 1
            else
               ib = 0
            end if
c
5400        continue
c
c increment surface point counter for convex surface
c
            nshape(1) = nshape(1) + 1
c
c add surface point area to contact area
c
            areac = areac + area
            iatnum(1) = iatom
            iatnum(2) = 0
            iatnum(3) = 0
            do 5450 k = 1,3
               outco(k) = ci(k) + ri * ua(k,i,iatom)
               outvec(k) = ua(k,i,iatom)
5450        continue
c
cxxx   
c
c      The point coordinates, area and normal vectors
c      are stored in common/conout.
c
cxxx   if (i.eq.1) then 
cxxx     ncon=j
cxxx   else
         ncon=ncon+1
cxxx   endif

       if (iturn.gt.1) then
         do 2221 k=1,3
           rc(k,ncon)=outco(k) 
           tnorm(k,ncon)=outvec(k)
2221     continue
         areas(ncon)=area
       endif
c
c four different output formats
caleko
c            if (ascii .and. long) then
c              write (9,5500) iatnum,ishape,outco,area,outvec,ib
c               write (9,5500) outco,area,outvec
cxxx
c5500           format(7e12.6)
c5500          format(3i5,i2,3f9.3,4f7.3,i2)
cxxx
c            else if (.not. ascii .and. long) then
c               write (9) iatnum,ishape,outco,area,outvec,ib
c            else if (ascii .and. .not. long) then
c               write (9,5550) iatnum,ishape,outco,ib
c5550           format(3i5,i2,3f9.3,i2)
c            else if (.not. ascii .and. .not. long) then
c               write (9) iatnum,ishape,outco,ib
c            end if
c
c end of nua loop
c
5600     continue
c
c end of iatom loop
c
5650  continue
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c write out messages
c
c for array dimensioning
c
      if (ibemout .ge. 2) write(6,5700) maxnb
5700  format(1x,i5,' neighbors maximum')
c
c close contact and before files
c
c      close(7)
c
c if van der waals surface we are finished
c
      if (rp .eq. 0.0d0) go to 8000
      close(4)
c
c orsr    orsr    orsr    orsr    orsr    orsr    orsr    orsr    orsr
c
c
c overlapping reentrant surface removal
c for non-symmetry-related probes
c probe diameter
c
      dp = 2 * rp
c
c diameter squared
c
      dp2 = dp ** 2
c
c radius squared
c
      rp2 = rp ** 2
c
c width for cubing algorithm
c     width = dp
c
c     set up cube arrays
c     first the integer coordinate arrays
c
      do 5850 iy = 1,ny
         do 5800 k = 1,3
            ico(k,iy) = (py(k,iy)-comin(k)-radmax-rp)/width + 1
            if (ico(k,iy) .lt. 1) ico(k,iy) = 1
            if (ico(k,iy) .gt. maxcub) call caserr
     +         ('cube coordinate too large')
5800     continue
5850  continue
c
c initialize head pointer array
c
       do 6000 k = 1,maxcub
         do 5950 j = 1,maxcub
            do 5900 i = 1,maxcub
               icube(i,j,k) = 0
5900        continue
5950     continue
6000  continue
c
c initialize linked list pointers
c
      do 6050 iy = 1,maxyon
         icuptr(iy) = 0
6050  continue
c
c set up head and later pointers for each yon probe
c
      do 6200 iy = 1,ny
c
c skip atoms with surface request numbers of zero
c
         i = ico(1,iy)
         j = ico(2,iy)
         k = ico(3,iy)
         if (icube(i,j,k) .le. 0) then
c
c     first atom in this cube
c
            icube(i,j,k) = iy
         else
c
c     add to end of linked list
c
            iptr = icube(i,j,k)
6100        continue
            if (icuptr(iptr) .le. 0) go to 6150
            iptr = icuptr(iptr)
            go to 6100
6150        continue
            icuptr(iptr) = iy
         end if
6200  continue
c
c reopen before file for reading
cxxx 5/91
c     open(4,form='unformatted',status='old')
      open(4,form='unformatted',file='con4',status='old')
c
c first pass
c gather victim probes
c
      nv = 1
c
c no victim probes if no yon probes
c
      if (ny .le. 0) go to 6950
      rewind(4)
c
c initialize victim hashing array
c
      do 6250 j = 1,maxvic
        ivic(j) = 0
6250  continue
c
c initialize index of free slot last used
c
      ifrlst = 0
c
c initialize probe record number
c
      i = 1
6300  continue
c
c check for victim overflow
c
      if (nv .gt. maxvic) call error(760,nv,0.0d0)
c
c read reentrant probe and points
c
      read (4,end=6950) molvic(nv),ishape,np,
     1 (pv(k,nv),k=1,3),(av(k,nv),k=1,3),yonprb,bury
      do 6350 j = 1,np
         read (4) n1(1),n2(1),n3(1),(s(k,1),k=1,3),area,yon(1)
6350  continue
      if (yonprb) go to 6900
c
c check if probe too far from symmetry element for possible overlap
c
      if (anrm(av(1,nv)) .gt. dp) go to 6900
c
c look for overlap with any yon probe in the same molecule
c use cubing algorithm to save time
c
c calculate which cube this probe lies in
c
      ici = (pv(1,nv)-comin(1)-radmax-rp)/width + 1
      if (ici .lt. 1) ici = 1
      if (ici .gt. maxcub) call caserr
     +        ('cube coordinate too large')
      icj = (pv(2,nv)-comin(2)-radmax-rp)/width + 1
      if (icj .lt. 1) icj = 1
      if (icj .gt. maxcub) call caserr
     +         ('cube coordinate too large')
      ick = (pv(3,nv)-comin(3)-radmax-rp)/width + 1
      if (ick .lt. 1) ick = 1
      if (ick .gt. maxcub) call caserr
     +          ('cube coordinate too large')
c
c check for overlap with probes in adjoining cubes
c
      do 6850 jck = ick-1,ick+1
         if (jck .lt. 1 .or. jck .gt. maxcub) go to 6850
         do 6800 jcj = icj-1,icj+1
            if (jcj .lt. 1 .or. jcj .gt. maxcub) go to 6800
            do 6750 jci = ici-1,ici+1
               if (jci .lt. 1 .or. jci .gt. maxcub) go to 6750
               jp = icube(jci,jcj,jck)
c
6400           continue
               if (jp .le. 0) go to 6700
               if (molyon(jp) .ne. molvic(nv)) go to 6650
               x = abs(py(1,jp) - pv(1,nv))
               if (x .ge. dp) go to 6650
               y = abs(py(2,jp) - pv(2,nv))
               if (y .ge. dp) go to 6650
               z = abs(py(3,jp) - pv(3,nv))
               if (z .ge. dp) go to 6650
               d2 = x ** 2 + y ** 2 + z ** 2
               if (d2 .ge. dp2) go to 6650
c
c check that probes face each other
c
               if (ddot(3,ay(1,jp),1,av(1,nv),1) .ge. 0.0d0) go to 6650
c
c new victim probe
c put into hashing table
c
               ihash = mod(i,maxvic) + 1
               if (ivic(ihash) .eq. 0) then
c
c empty slot
c
                  ivic(ihash) = i
                  ivicp(ihash) = 0
               else
                  iprev = ihash
                  iptr = ivicp(ihash)
6450              continue
c
c check for end of linked list
c
                  if (iptr .eq. 0) go to 6500
                   iprev = iptr
                  iptr = ivicp(iptr)
                  go to 6450
6500              continue
c
c look for a free slot
c
                  do 6550 ifree = ifrlst+1,maxvic
                     if (ivic(ifree) .eq. 0) go to 6600
6550              continue
                  call caserr('victim oveflow')
6600              continue
c
c store record number in free slot
c
                  ivic(ifree) = i
                  ivicp(iprev) = ifree
                  ivicp(ifree) = 0
c
c new index to last free slot used
c
                  ifrlst = ifree
               end if
               nv = nv + 1
c
c one overlap makes this probe a victim
c we don't need to check any more
c
               go to 6900
6650           continue
               jp = icuptr(jp)
               go to 6400
6700           continue
6750        continue
6800     continue
6850  continue
c
c end of yon probe loop
c skip to here if finished with hunt for overlapping probes
c
6900  continue
      i = i + 1
      go to 6300
c
c skip to here if there are no yon probes and hence no victims
c
6950  continue
      nv = nv - 1
c
c open reentrant file for writing
c
c      if (ascii) then
cxxx 5/91
caleko
c  Commented out
c
c        open(8,status='new')
c         open(8,form='formatted',file='rentry',status='unknown')
c      else
c        open(8,form='unformatted',status='new')
c         open(8,form='unformatted',file='rentry',status='unknown')
cxxx
c      end if
c      rewind(8)
c
c second pass
c read, check and write surface points
c
      rewind(4)
      i = 1
7000  continue
      read (4,end=7850) imol,ishape,np,
     1 (p(k),k=1,3),(a(k),k=1,3),yonprb,bury
c
c no points can be eaten if this probe is neither yon nor a victim
c
      neat = 0
      nyeat = 0
      if (ny .le. 0) go to 7450
      ipt = 0
c
c determine if probe is a yon or victim probe
c
      if (.not. yonprb) go to 7050
c
c we've got a yon probe here
c
      ipt = 2
      go to 7200
7050  continue
      if (nv .le. 0) go to 7450
c
c hash into table of victim probes
c
      iptr = mod(i,maxvic) + 1
7100  continue
      if (iptr .eq. 0) go to 7200
      if (ivic(iptr) .eq. 0) go to 7200
      if (ivic(iptr) .eq. i) go to 7150
      iptr = ivicp(iptr)
      go to 7100
7150  continue
c
c we've got a victim
c
      ipt = 1
7200  continue
c
      if (ipt .le. 0) go to 7450
c
c check this victim or yon probe against all yon probes
c
      do 7300 j = 1,ny
         if (imol .ne. molyon(j)) go to 7300
         if (dist2(p,py(1,j)) .ge. dp2) go to 7300
         if (ddot(3,a,1,ay(1,j),1) .ge. 0.0d0) go to 7300
c
c this yon probe could eat some of the probe's points
c
         neat = neat + 1
         nyeat = nyeat + 1
         if (neat .gt. maxeat) call error(830,neat,0.0d0)
         do 7250 k = 1,3
            eat(k,neat) = py(k,j)
7250     continue
c
c end of yon probe loop
c
7300  continue
c
c only yon probes can have their points eaten by victims
c
      if (ipt .le. 1) go to 7450
c
c check this yon probe against all victim probes
c
      do 7400 j = 1,nv
         if (imol .ne. molvic(j)) go to 7400
         if (dist2(p,pv(1,j)) .ge. dp2) go to 7400
         if (ddot(3,a,1,av(1,j),1) .ge. 0.0d0) go to 7400
c
c this victim probe could eat some of the probe's points
c
         neat = neat + 1
         if (neat .gt. maxeat) call error(850,neat,0.0d0)
         do 7350 k = 1,3
            eat(k,neat) = pv(k,j)
7350     continue
c
c end of victim probe loop
c
7400  continue
c
c skip to here if victim or both probe overlap checks omitted
c
7450  continue
c
c read the surface points belonging to the probe
c
       do 7750 j = 1,np
         read (4) n1(1),n2(1),n3(1),(s(k,1),k=1,3),area,yon(1)
         if (neat .le. 0) go to 7550
c
c check surface point against all eaters of this probe
c
         do 7500 k = 1,neat
c
c victim probes cannot eat non-yon points of yon probes
c
            if (yonprb .and. .not. yon(1) .and. k .gt. nyeat) go to 7500
            if (dist2(eat(1,k),s(1,1)) .lt. rp2) go to 7700
7500     continue
c
c skip to here if no overlapping probes could eat this point
c
7550     continue
         do 7600 k = 1,3
            outvec(k) = (p(k) - s(k,1))/rp
7600     continue
c
c reentrant surface point
c
         nshape(ishape) = nshape(ishape) + 1
         arear = arear + area
c
c mark whether buried
c
         ib = 0
         if (ibury .gt. 0 .and. bury) ib = 1
c
c four possible output formats
c
         iatnum(1) = n1(1)
         iatnum(2) = n2(1)
         iatnum(3) = n3(1)
         do 7650 k = 1,3
            outco(k) = s(k,1)
7650     continue
c
caleko
c
c      The point coordinates, area and normal vectors
c      are stored in common/conout.
c
cxxx   if (i.eq.1) then 
cxxx     ncon=j
cxxx   else
         ncon=ncon+1
cxxx   endif

       if (iturn.gt.1) then
         do 2222 k=1,3
           rc(k,ncon)=outco(k) 
           tnorm(k,ncon)=outvec(k)
2222     continue
         areas(ncon)=area
       endif
c
caleko
c  Commented out

c         if (ascii .and. long) then
c            write (8,5500) iatnum,ishape,outco,area,outvec,ib
c            write (8,5500) outco,area,outvec
c         else if (.not. ascii .and. long) then
c            write (8) iatnum,ishape,outco,area,outvec,ib
c         else if (ascii .and. .not. long) then
c            write (8,5550) iatnum,ishape,outco,ib
c         else if (.not. ascii .and. .not. long) then
c            write (8) iatnum,ishape,outco,ib
c         end if
c         go to 7750
7700     continue
         nlost(ishape) = nlost(ishape) + 1
c

c end of np loop
c
7750  continue
c
c end of i loop
c
7800  continue
      i = i + 1
      go to 7000
7850  continue
c
c close files
c
      close(4)
c      close(8)
caleko

c
c messages
c
      if (ibemout .ge. 2) then
      write (6,7900) ny,nv
7900  format(1x,i5,' yon and ',i5,' victim probes')
      write(6,7950) nlost(2),nlost(3)
7950  format(1x,i5,' saddle and ',i5,
     1 ' concave surface points removed during non-symmetry orsr')
      endif
8000  continue
c
c write out how many points
c
      if (ibemout .ge. 2) then
      write (6,8050) nshape(1),nshape(2),nshape(3)
8050  format(1x,i5,' contact and ',i5,' saddle and ',
     1 i5,' concave surface points')
      write (6,8100) nshape(1) + nshape(2) + nshape(3)
8100  format(1x,i8,' total surface points')
      write (6,8150) areac,arear,areac+arear
c
8150  format(1x,'contact area:',f10.3,2x,'reentrant area:',f10.3,
     1 2x,'total area:',f10.3)
      endif
caleko
c     End of loop for defining the right surface density
c
      if (.not. isurdens) then
        if ((ncon.le.nmaxcon).and.(ncon.ge.nmncon)) goto 444
          if (iround.eq.1) then
            iround=iround+1
            goto 111
          else
        if(iround.ge.maxit) call error(880,0,0.0d0)
            iround=iround+1
            goto 222
          endif
 444    continue
      endif
 448  continue

      if ((iturn .gt. 1) .and. (ibemout .ge. 5)) then
         write(iwr,8887)
 8887    format(/,' Connolly surface: coordinates, ',
     1   'area, and normal vectors')
         do 8888 k=1,ncon
            write(iwr,8889) k,(rc(i,k),i=1,3), areas(k),
     1        (tnorm(i,k),i=1,3)
 8889       format(i4,1x,7(f12.6,1x))
 8888    continue
      endif
cxxx
      isurdens = .true.
      return
      end
c
c subroutines and functions
c
c
c general vector and matrix routines
c
      function dist3(a,b)
      implicit REAL (a-h,o-z)
c------
c distance between a and b
c------
      dimension a(3), b(3)
      dist3 = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)
      return
      end
      function dist2(a,b)
      implicit REAL (a-h,o-z)
c------
c distance between a and b squared
c------
      dimension a(3), b(3)
      dist2 = (a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2
      return
      end
      function anrm(a)
      implicit REAL (a-h,o-z)
c------
c norm of a
c------
      dimension a(3)
      anrm = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
      return
      end
c      function dot(a,b)
c      implicit REAL (a-h,o-z)
cc------
cc dot product
cc------
c      dimension a(3), b(3)
c      dot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
c      return
c      end
      subroutine croscon(a,b,c)
      implicit REAL (a-h,o-z)
c------
c cross product
c------
      dimension a(3), b(3), c(3)
      c(1) = a(2) * b(3) - a(3) * b(2)
      c(2) = a(3) * b(1) - a(1) * b(3)
      c(3) = a(1) * b(2) - a(2) * b(1)
      return
      end
      subroutine multv(v,a,w)
      implicit REAL (a-h,o-z)
c------
c multiply v by a giving w
c------
      dimension a(3,3)
      dimension v(3), w(3)
      do 50 i = 1, 3
         w(i) = a(i,1)*v(1) + a(i,2)*v(2) + a(i,3)*v(3)
50    continue
      return
      end
      subroutine vnorm(a,b)
      implicit REAL (a-h,o-z)
c------
c normalize a giving b
c------
      dimension a(3), b(3)
      v = anrm(a)
      do 50 k = 1,3
         b(k) = a(k) / v
50    continue
      return
      end
      subroutine vperp(a,b)
      implicit REAL (a-h,o-z)
c------
c return b perpendicular to a
c------
      dimension a(3), b(3), p(3)
c
c find smallest component
c
      small = 10000.0d0
      m = 0
      do 50 k = 1,3
         if (abs(a(k)) .ge. small) go to 50
         small = abs(a(k))
         m = k
50    continue
      do 100 k = 1,3
         b(k) = 0.0d0
         if (k .eq. m) b(k) = 1.0d0
100   continue
c
c take projection along a
c
      dt = a(m) / (a(1)**2 + a(2)**2 + a(3)**2)
      do 150 k = 1, 3
         p(k) = dt * a(k)
c
c subtract projection from b

          b(k) = b(k) - p(k)
150   continue
c
c renormalize b
c
      call vnorm(b,b)
      return
      end
      subroutine cat(a,b)
      implicit REAL (a-h,o-z)
c------
c concatenate matrix b into matrix a
c------
      dimension a(3,3), b(3,3)
      dimension temp(3,3)
      do 100 i = 1,3
         do 50 j = 1,3
            temp(i,j) = a(i,1)*b(1,j) + a(i,2)*b(2,j) + a(i,3)*b(3,j)
50       continue
100   continue
      do 200 i = 1,3
         do 150 j = 1,3
            a(i,j) = temp(i,j)
150      continue
200   continue
      return
      end
      subroutine conj(h,g,ghgt)
      implicit REAL (a-h,o-z)
c------
c conjugate matrix g with matrix h giving ghgt
c------
      dimension h(3,3), g(3,3), ghgt(3,3)
      dimension gt(3,3)
c
c initialize ghgt matrix to identity
c concatenate g h gt
c
      call imatx(ghgt)
      call cat(ghgt,g)
      call cat(ghgt,h)
c
c calculate gt
c
      do 100 k = 1,3
         do 50 l = 1,3
            gt(k,l) = g(l,k)
50       continue
100   continue
      call cat(ghgt,gt)
      return
      end
      function det(a,b,c)
      implicit REAL (a-h,o-z)
c------
c return triple product of the three vectors
c------
      dimension a(3), b(3), c(3)
      dimension ab(3)
      call croscon(a,b,ab)
      det = ddot(3,ab,1,c,1)
      return
      end
c
c geometric routines
c
c
      logical function collid(p,rp,cnbr,ernbr,mnbr,nnbr,maxnbr,ishape,
     1 jnbr,knbr,molnbr,imol,lkf,lknbr)
      implicit REAL (a-h,o-z)
c------
c collision check of probe with neighboring atoms
c belonging to the same molecule
c------
      dimension p(3), ernbr(maxnbr)
      dimension cnbr(3,maxnbr)
      logical mnbr(maxnbr)
      integer*2 molnbr(maxnbr)
      integer*2 lknbr(maxnbr)
      integer*2 imol,ishape
c
c check whether probe is too close to any neighbor
c
      i = lkf
      go to 100
50    continue
      i = lknbr(i)
100   continue
      if (i .eq. 0) go to 150
      vect1 = abs(p(1) - cnbr(1,i))
      if (vect1 .ge. ernbr(i)) go to 50
      vect2 = abs(p(2) - cnbr(2,i))
      if (vect2 .ge. ernbr(i)) go to 50
      vect3 = abs(p(3) - cnbr(3,i))
      if (vect3 .ge. ernbr(i)) go to 50
      if (i .eq. jnbr .or. i .eq. knbr) go to 50
      sr2 = ernbr(i) ** 2
      dd2 = vect1 ** 2 + vect2 ** 2 + vect3 ** 2
      if (dd2 .ge. sr2) go to 50
      collid = .true.
      return
150   continue
      collid = .false.
      return
      end
      logical function buried(p,rp,cnbr,rnbr,mnbr,nnbr,maxnbr,ishape,
     1 jnbr,knbr,molnbr,imol)
      implicit REAL (a-h,o-z)
c------
c collision check of probe with neighboring atoms
c belonging to a different molecule
c------
      dimension p(3), rnbr(maxnbr)
      dimension cnbr(3,maxnbr)
      logical mnbr(maxnbr)
      integer*2 molnbr(maxnbr)
      integer*2 imol,ishape
c
      if (nnbr .le. 0) go to 100
c
c check whether probe is too close to any neighbor
c
      do 50 i = 1, nnbr
         if (imol .eq. molnbr(i)) go to 50
         if (ishape .gt. 1 .and. i .eq. jnbr) go to 50
         if (ishape .eq. 3 .and. (i .eq. knbr .or. .not. mnbr(i)))
     1 go to 50
         sumrad = rp + rnbr(i)
         vect1 = abs(p(1) - cnbr(1,i))
         if (vect1 .ge. sumrad) go to 50
         vect2 = abs(p(2) - cnbr(2,i))
         if (vect2 .ge. sumrad) go to 50
         vect3 = abs(p(3) - cnbr(3,i))
         if (vect3 .ge. sumrad) go to 50
         sr2 = sumrad ** 2
         dd2 = vect1 ** 2 + vect2 ** 2 + vect3 ** 2
         if (dd2 .lt. sr2) go to 150
50    continue
100   continue
      buried = .false.
      go to 200
150   continue
      buried = .true.
200   continue
      return
      end
      subroutine genun(u,n)
      implicit REAL (a-h,o-z)
c------
c generate unit vectors over sphere
c------
      dimension u(3,n)
c
      nequat = sqrt(n * 3.14159d0)
      nvert = 0.5d0 * nequat
      if (nvert .lt. 1) nvert = 1
      nu = 0
      do 100 i = 0,nvert
         fi = (3.14159d0 * i) / nvert
         z = cos(fi)
         xy = sin(fi)
         nhor = nequat * xy
         if (nhor .lt. 1) nhor = 1
         do 50 j = 0,nhor-1
            fj = (2d0 * 3.14159d0 * j) / nhor
            x = cos(fj) * xy
            y = sin(fj) * xy
            if (nu .ge. n) go to 150
            nu = nu + 1
            u(1,nu) = x
            u(2,nu) = y
            u(3,nu) = z
50       continue
100   continue
150   continue
      n = nu
      return
      end
c
c error message subroutine
c
      subroutine error(number,int,float)
      implicit REAL (a-h,o-z)
c------
c error message subroutine
c------
      dimension list(16)
      data list/120,127,130,140,150,160,170,
     1 210,320,440,480,720,760,830,850,880/
c
      do 50 i = 1,16
         if (list(i) .eq. number) go to 150
50    continue
      write (6,100)
100   format(1x,'error of unidentifiable type')
      call caserr('error of unidentifiable type')
150   continue
c
      go to (200,300,400,500,600,700,800,900,1000,1100,1200,
     1 1300,1400,1500,1600,1700) i
c
200   write (6,250) number,float
250   format(1x,'error',i5,2x,'negative probe radius:',f10.5)
      stop
300   write (6,350) number,int
350   format(1x,'error',i5,2x,'bad buried surface flag:',i5)
      stop
400   write (6,450) number,int
450   format(1x,'error',i5,2x,'too few or too many atom types:',i5)
      stop
500   write (6,550) number,float,int
550   format(1x,'error',i5,2x,'negative atom radius:',
     1 f10.5,' atom',i5)
      stop
600   write (6,650) number,int
650   format(1x,'error',i5,2x,'too many atoms:',i5)
      stop
700   write (6,750) number,int
750   format(1x,'error',i5,2x,
     1 'invalid surface request number for atom:',i5)
      stop
800   write (6,850) number,int
850   format(1x,'error',i5,2x,'invalid atom type for atom:',i5)
      stop
900   write (6,950) number,int
950   format(1x,'error',i5,2x,'too many neighbors:',i5)
      stop
1000  write (6,1050) number,int
1050  format(1x,'error',i5,2x,'too many points for reentrant probe:',
     1 i5)
      stop
1100  write (6,1150) number,int
1150  format(1x,'error',i5,2x,'too many points for arc:',i5)
      stop
1200  write (6,1250) number,int
1250  format(1x,'error',i5,2x,'too many points for reentrant probe:',
     1 i5)
      stop
1300  write (6,1350) number,int
1350  format(1x,'error',i5,2x,'too many yon probes:',i5)
      stop
1400  write (6,1450) number,int
1450  format(1x,'error',i5,2x,'too many victim probes:',i5)
      stop
1500  write (6,1550) number,int
1550  format(1x,'error',i5,2x,'too many eaters:',i5)
      stop
1600  write (6,1650) number,int
1650  format(1x,'error',i5,2x,'too many eaters:',i5)
      stop
1700  write (6,1750) number
1750  format(1x,'error',i5,2x,'point generation does not converge, try rerunning
     1 with a different number of points. DIELOUT SOME in the input conveys
     2 information about the convergence of the number of points')
      stop
      end
