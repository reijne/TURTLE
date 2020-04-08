      function radius(name,iradex,charge,afct)
c=======================================================================
c     this function assigns a radius to an atom
c     value is in bohr according to :
c
c     bondi:j.phys.chem.(1964) 68 441
c     table i:  rw
c
cmas  changed function feb 7, 1996
c=======================================================================
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c * * * this is for communication with hondo8
c
      character *8 errmsg(3)
c
INCLUDE(comdrf/iofil)
c
      parameter (ndif=106)
      character*2 names(ndif),name*4
      dimension rad(ndif)
      dimension afct(2)
c
      parameter ( nfrec = 18 )
      dimension rfrec0(nfrec),rfrec1(nfrec),rfrec2(nfrec),rfrec3(nfrec)
c
      data errmsg/'program','stop in','-radius-'/
      data names/'h ','he','li','be','b ','c ','n ','o ','f ','ne',
     *           'na','mg','al','si','p ','s ','cl','ar','k ','ca',
     *           'sc','ti','v ','cr','mn','fe','co','ni','cu','zn',
     *           'ga','ge','as','se','br','kr','rb','sr','y ','zr',
     *           'nb','mo','tc','ru','rh','pd','ag','cd','in','sn',
     *           'sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     *           'pm','sm','eu','gd','tb','dy','ho','er','tm','yb',
     *           'lu','hf','ta','w ','re','os','ir','pt','au','hg',
     *           'tl','pb','bi','po','at','rn','fr','ra','ac','th',
     *           'pa','u ','np','pu','am','cm','bk','cf','es','fm',
     *           'md','no','lr','e','qq','xx'/
c
c    -----  note: e refers to a center without nuclear charge (bsse)
c
      data rad  /
c         h        he       li       be       b        c
     1 0.2267d1,0.2646d1,0.0000d1,0.0000d1,0.0000d1,0.3213d1,
c         n        o        f        ne       na       mg
     1 0.2929d1,0.2872d1,0.2778d1,0.2910d1,0.0000d1,0.0000d1,
c         al       si       p        s        cl       ar       k
     1 0.0000d1,0.3968d1,0.3401d1,0.3402d1,0.3307d1,0.3553d1,0.0000d1,
c         ca       sc       ti       v        cr       mn
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         fe       co       ni       cu       zn       ga
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         ge       as       se       br       kr       rb
     1 0.0000d1,0.3500d1,0.3590d1,0.3500d1,0.3817d1,0.0000d1,
c         sr       y        zr       nb       mo       tc
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         ru       rh       pd       ag       cd       in
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         sn       sb       te       i        xe       cs
     1 0.0000d1,0.0000d1,0.3892d1,0.3741d1,0.3817d1,0.0000d1,
c         ba       la       ce       pr       nd       pm
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         sm       eu       gd       tb       dy       ho
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         er       tm       yb       lu       hf       ta
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         w        re       os       ir       pt       au
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         hg       tl       pb       bi       po       at
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         rn       fr       ra       ac       th       pa
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         u        np       pu       am       cm       bk
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         cf       es       fm       md       no       lr
     1 0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,0.0000d1,
c         e        qq       xx
     1 0.0000d1,0.0000d1,0.1000d1/
c
c=======================================================================
c     data obtained by fitting of polynoma (3rd degree) to rvdw-values
c     for ions with charge ranging from q=-2 to q=+2
c     by monte carlo sampling of ab initio charge densities
c
c                                        2         3
c     rvdw(q) = rvdw(0) + r1 * q + r2 * q  + r3 * q
c
c     values are in a (angstroms) and e (electron charge)
c
c     fitting done by vladimir frecer
c=======================================================================
c-----------------------------------------------------------------------
c     h
c-----------------------------------------------------------------------
      data rfrec0(1),rfrec1(1),rfrec2(1),rfrec3(1)/
     .       1.520d0, -1.054d0, -0.293d0,  0.000d0/
c-----------------------------------------------------------------------
c     he
c-----------------------------------------------------------------------
      data rfrec0(2),rfrec1(2),rfrec2(2),rfrec3(2)/
     .       1.330d0, -0.300d0, -0.129d0, -0.037d0/
c-----------------------------------------------------------------------
c     li
c-----------------------------------------------------------------------
      data rfrec0(3),rfrec1(3),rfrec2(3),rfrec3(3)/
     .       2.210d0, -0.996d0,  0.003d0,  0.096d0/
c-----------------------------------------------------------------------
c     be
c-----------------------------------------------------------------------
      data rfrec0(4),rfrec1(4),rfrec2(4),rfrec3(4)/
     .       2.180d0, -0.335d0, -0.103d0, -0.051d0/
c-----------------------------------------------------------------------
c     b
c-----------------------------------------------------------------------
      data rfrec0(5),rfrec1(5),rfrec2(5),rfrec3(5)/
     .       2.060d0, -0.298d0, -0.007d0,  0.009d0/
c-----------------------------------------------------------------------
c     c
c-----------------------------------------------------------------------
      data rfrec0(6),rfrec1(6),rfrec2(6),rfrec3(6)/
     .       1.930d0, -0.234d0, -0.013d0,  0.008d0/
c-----------------------------------------------------------------------
c     n
c-----------------------------------------------------------------------
      data rfrec0(7),rfrec1(7),rfrec2(7),rfrec3(7)/
     .       1.800d0, -0.170d0, -0.015d0,  0.003d0/
c-----------------------------------------------------------------------
c     o
c-----------------------------------------------------------------------
      data rfrec0(8),rfrec1(8),rfrec2(8),rfrec3(8)/
     .       1.700d0, -0.140d0, -0.014d0,  0.002d0/
c-----------------------------------------------------------------------
c     f
c-----------------------------------------------------------------------
      data rfrec0(9),rfrec1(9),rfrec2(9),rfrec3(9)/
     .       1.610d0, -0.096d0, -0.012d0,  0.004d0/
c-----------------------------------------------------------------------
c     ne
c-----------------------------------------------------------------------
      data rfrec0(10),rfrec1(10),rfrec2(10),rfrec3(10)/
     .       1.510d0, -0.098d0, -0.005d0,  0.003d0/
c-----------------------------------------------------------------------
c     na
c-----------------------------------------------------------------------
      data rfrec0(11),rfrec1(11),rfrec2(11),rfrec3(11)/
     .       2.240d0, -0.835d0, -0.009d0,  0.094d0/
c-----------------------------------------------------------------------
c     mg
c-----------------------------------------------------------------------
      data rfrec0(12),rfrec1(12),rfrec2(12),rfrec3(12)/
     .       2.420d0, -0.375d0, -0.062d0, -0.039d0/
c-----------------------------------------------------------------------
c     al
c-----------------------------------------------------------------------
      data rfrec0(13),rfrec1(13),rfrec2(13),rfrec3(13)/
     .       2.410d0, -0.365d0,  0.021d0,  0.005d0/
c-----------------------------------------------------------------------
c     si
c-----------------------------------------------------------------------
      data rfrec0(14),rfrec1(14),rfrec2(14),rfrec3(14)/
     .       2.330d0, -0.291d0,  0.007d0,  0.006d0/
c-----------------------------------------------------------------------
c     p
c-----------------------------------------------------------------------
      data rfrec0(15),rfrec1(15),rfrec2(15),rfrec3(15)/
     .       2.260d0, -0.275d0,  0.015d0,  0.006d0/
c-----------------------------------------------------------------------
c     s
c-----------------------------------------------------------------------
      data rfrec0(16),rfrec1(16),rfrec2(16),rfrec3(16)/
     .       2.150d0, -0.240d0,  0.020d0,  0.004d0/
c-----------------------------------------------------------------------
c     cl
c-----------------------------------------------------------------------
      data rfrec0(17),rfrec1(17),rfrec2(17),rfrec3(17)/
     .       2.080d0, -0.208d0,  0.017d0,  0.001d0/
c-----------------------------------------------------------------------
c     ar
c-----------------------------------------------------------------------
      data rfrec0(18),rfrec1(18),rfrec2(18),rfrec3(18)/
     .       1.970d0, -0.251d0,  0.010d0,  0.017d0/
      data zero,one,two,three/0.0d0,1.0d0,2.0d0,3.0d0/
      data bohr/0.529177249d0/
c-----------------------------------------------------------------------
c     begin
c-----------------------------------------------------------------------
cafc
cafc  gamess-names for dummy-atoms
      if (name(1:2) .eq. 'bq') then
        radius = 0.0d0
        return
      endif
      do i=1,ndif
        if (name(:2).eq.names(i))  goto 20
      enddo
c-----------------------------------------------------------------------
c     error : symbol not found
c-----------------------------------------------------------------------
      write(iwr,*) ' no radius found for ',name
      call hnderr(3,errmsg)
c
   20 if (iradex.eq.4.and.i.le.nfrec) then
c-----------------------------------------------------------------------
c     use data obtained by vladimir frecer
c     remember that radius is then in angstroms
c-----------------------------------------------------------------------
         radius = rfrec0(i) + rfrec1(i)*charge + rfrec2(i)*charge**two
     .                      + rfrec3(i)*charge**three
         radius = radius / bohr
      elseif (iradex.eq.4) then
c-----------------------------------------------------------------------
c     try to use frecer's model for an atom which was not fitted
c-----------------------------------------------------------------------
         write (iwr,*) 'no frecer values for ' ,name
         call hnderr(3,errmsg)
      elseif (iradex.eq.0) then
c-----------------------------------------------------------------------
c     use value from table
c-----------------------------------------------------------------------
         radius = rad(i)
      else
c-----------------------------------------------------------------------
c     compute radius from polarizability if iradex = 1 or 2
c-----------------------------------------------------------------------
         radius = afct(iradex)*alfa(name,iradex,ier)**(one/three)
      endif
      return
      end
      subroutine drfzfa(xexp)
c------
c      calculates and updates the (expanded) field due to
c      external charges at the expansion centra of the
c      quantum mechanical system
c
c      the interaction energy between the external charges
c      and the nuclei (extnuc) is also calculated
c
c      penetration (overlap) effects
c      resulting from close contacts
c      can be accounted for (modxza flag)
c
c      --------  p.th. van duijnen, ibm-kingston 1985, and
c               groningen, dec. 1991.
c
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy variables
c
      dimension xexp(3,nexp)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/mollab)
c
cmw      include 'drf/dimpar'
INCLUDE(../m4/common/drfopt)
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/drfzfa)
INCLUDE(comdrf/rfene)
INCLUDE(comdrf/extinf)
INCLUDE(comdrf/grinf)
INCLUDE(comdrf/drfamb)
c
INCLUDE(comdrf/drfbem)
INCLUDE(comdrf/rad)
c
INCLUDE(comdrf/runpar)
INCLUDE(comdrf/mcener)
c
c-----  local variables
c
      dimension p(3), q(3), pq(3), w(3)
      dimension b(3,3)
c
c
      character*80 card
      character*16 namj, nami
      logical group
      logical ipol
c
      data third/.33333333333333333333333333333d00/
      data sixth/.16666666666666666666666666667d00/
      data two,three,four,twelve/2.0d00,3.0d00,4.0d00,12.0d00/
      data zero,one,pt5,pt75/0.0d00,1.0d00,0.5d00,0.75d00/
      data small /1.0d-03/
c
c-----  begin
c
      ier = 0
c
c-----  read expansion centers.
c       = # internal  atoms
c       (+ center of internal nuclear charge)
c
      call daread(idafdrf,iodadrf,xexp,3*nexp,1)
c
      if (mcupdt) then
c 1-----
c  -----  if update of zfa is required, read zfa from disk
c
        call daread(idafdrf,iodadrf,zfa,4*nzfa,3)
        call daread(idafdrf,iodadrf,zpa,nzfa,4)
        call daread(idafdrf,iodadrf,dzfa,12*nzfa,5)
        call daread(idafdrf,iodadrf,dpseu,3*nzfa,6)
c
c  -----  update field at expansion centra according to
c         monte carlo move
c
        call expupd(xexp,extnuc(iactst),repmod(iactst),
     1              extnucg,repmodg,iactst)
c
        call dawrit(idafdrf,iodadrf,zfa,4*nzfa,76,navdrf)
        call dawrit(idafdrf,iodadrf,zpa,nzfa,77,navdrf)
        call dawrit(idafdrf,iodadrf,dzfa,12*nzfa,78,navdrf)
        call dawrit(idafdrf,iodadrf,dpseu,3*nzfa,79,navdrf)
c 1-----
      else
c 1-----
c  -----  calculate zfa, dzfa (derivative) and model repulsion
c
        pseufac = pt75
c
c  -----  initialize pseudo repulsion and nuclear-external charge
c         interaction
c
        epseu = zero
        extnuc(iactst) = zero
        repmod(iactst) = zero
        do 400, igr = 1, ngran
          extnucg(igr,iactst) = zero
          repmodg(igr,iactst) = zero
  400   continue
c
c  -----  nspec counts the external points close to internal centra,
c         for which a model repulsion is calculated, and for which
c         the "reverse" dispersion is calculated, if required
c
        nspec = 0
c
        call clear(zfa,4*nzfa)
        call clear(zpa,nzfa)
        call clear(dzfa,12*nzfa)
        call clear(dpseu,3*nzfa)
c
c  -----  loop over the expansion centra
c
        do 500, j = 1, nexp
c   2-----
c    -----  exclude ambiguous atom: its interaction with the
c           discrete classical system is treated classically
c           in subr. drfzfp
c
          if (nambpt(j) .ne. 0) goto 500
c
          if (j .le. nat) then
c     3-----
c      -----  current expansion centre is at a nucleus
c
            za = czan(j)
            namj = anam(j)
            alfj = alfat(j)
c           alfj = alfa(namj,0,ier)
c
c      -----  # of  valence electrons
c
            call drfnval(namj,nvalj,znuc)
c
c      -----  scaling factor for slater-kirkwood estimate
c             of dispersion and repulsion energy
c
            if (alfj .ne. zero) then
              aovernj= sqrt(alfj/nvalj)
            endif
c     3-----
          else
c     3-----
c      -----  non-nuclei are given a polarizability of 1.
c
            alfj = one
c     3-----
          endif
c
c    -----  pointers to arrays w.r.t. expansion center -j-
c
          jf = (j-1)*4
          nzf = jf*3
c
c    -----  position vector -j- into -q-
c
          q(1) = xexp(1,j)
          q(2) = xexp(2,j)
          q(3) = xexp(3,j)
c
c    -----  pointer to external points in -xtpts-
c           with special properties
c           imp: polarizable external points
c
          imp = 1
c
c    -----  loop over external points
c
          do 300, ii = 1, nxtpts
c     3-----
c      -----  check whether external point is polarizable
c
            ipol = ii .eq. mpol(imp)
c
c      -----  check whether external point is an ambiguous
c             atom. the interaction with the qm system is taken care of
c             by inclusion into the qm system (whereby it is treated exa
c
            if (ncutpt(ii) .ne. 0) goto 290
c
c      -----  position vector of external point -ii- in -p-
c
            p(1) = xpts(1,ii)
            p(2) = xpts(2,ii)
            p(3) = xpts(3,ii)
c
            call distab(q,p,pq,dist)
c
            if (dist .lt. small) then
              write(iwr,95) anam(j), nxcent(ii), dist
  95          format(/,1x,'WARNING: interaction ',
     1   'between ', a16, ' and ', a16, ' at', e15.8,
     2   ' bohr distance: skipped')
              goto 300
            endif
c
            dmind1 = one/dist
            dmin3 = dmind1*(dmind1**2)
            factp = one
            factf = one
            v = one
c
            alfi = alfext(ii)
            if (alfi .eq. zero) alfi = one
c
c           if (chrg(ii) .eq. zero) goto 250
c
c      -----  get index of analysis group of external point -ii-
c
            igr = igranl(ii)
c
c      -----  account for penetration effects (optional)
c
            if (modxza .ne. 0) then
c       4-----
c        -----  set distance criterium for penetration effects
c               depending on polarizabilities of -j- and -ii-
c
              s = (alfj*alfi)**sixth
              v = dist/s
c
c        -----  check distance criterium
c               modify interactions (thole's cone model)
c
              if (ithole .eq. 1) then
c         5-----
                if (v .le. afact) then
                  av = v/afact
                  factp = av**4 - two*av**3 + two*av
                  factf = four*av**3 - three*av**4
                endif
              else if (ithole .eq. 2) then
                au = afact*v
                factp = (one - (pt5*au + one)*exp(-au))
                factf = (one - (pt5*au**2 + au + one)*exp(-au))
c         5-----
              endif
c       4-----
            endif
            factp = factp*dmind1
            factf = factf*dmin3
c
c      -----  update coulomb interaction with internal nuclei
c
            if ((field(:4) .ne. ' ') .and. (j .le. nat)) then
              extnuca = za*chrg(ii)*factp
              extnucg(igr,iactst) = extnucg(igr,iactst) + extnuca
              extnuc(iactst) = extnuc(iactst) + extnuca
            endif
c
c      -----  update -zfa-
c
c             -zfa- : sum(p) {f(p;q) + v(p;q) - f(p;q).q} * q(p)
c
c             the field operator and potential of all external charges
c             in p (-p-) at expansion centre q (-j- = -q-)
c
            sum = zero
            do 210, k = 1, 3
              zfa(jf+k,igr) = zfa(jf+k,igr) + chrg(ii)*(pq(k)*factf)
              sum = sum - q(k)*(pq(k)*factf)
  210       continue
c
c      -----  the expanded term f(p;q).q is not to be added if
c             the field of the external charges is treated
c             as a perturbation and the qm charge distribution
c             is modelled by point charges (and dipoles)
c
            zfa(jf+4,igr) = zfa(jf+4,igr) + chrg(ii)*(sum+factp)
c
c      -----  update -zpa- : the potential of the external charges
c             at the expansion centra
c
            zpa(j,igr) = zpa(j,igr) + chrg(ii)*factp
c
c      -----  classical field gradient tensor into -b-
c
            call clear(b,9)
c
c      -----  b(3,3) = t(p;q) : field gradient of charge in -p- at -q-
c             the gradient is scaled according to v
c
            call drftpq(pq,dmind1,dmin3,ithole,afact,v,b)
c
c      -----  multiply position vector of -j- with -b- to give -w-
c
            call matvec(b,q,w,3,.false.)
c
c      -----  update derivative of -zfa-
c
            do 214, k = 1, 3
              do 213, l = 1, 3
                dzfa(nzf+(k-1)*4+l,igr) = dzfa(nzf+(k-1)*4+l,igr)
     1                               - b(k,l)*chrg(ii)
  213         continue
              dzfa(nzf+(k-1)*4+4,igr) = dzfa(nzf+(k-1)*4+4,igr)
     1                              + w(k)*chrg(ii)
  214       continue
c
  250       continue
c
c      -----  non bonding repulsion taken from charmm:
c             j.comp.chem. 4 (1983) 213
c
c      -----  this is called the pseudo- or model repulsion throughout
c             the program
c
            if ((field(:4) .ne. ' ') .and.
     1          (j .le. nat) .and.
     1          (iqmclr .eq. 1) .and. 
     1          (namj(:2) .ne. 'bq') .and.
     1          (nxcent(ii)(:2) .ne. 'qq') .and.
     2          (nxcent(ii)(:2) .ne. 'e ') .and.
     3          (nxcent(ii)(:2) .ne. '  ')) then
c       4-----
c        -----  this might be a close atom: apply distance
c               criterium for pauli repulsion
c
              if (dist .lt. dstmin) then
c         5-----
c          -----  store external point counter in array relating to
c                 atom, for calculation of exact integrals
c
c          -----  calculate pauli repulsion and
c                 derivative, using some model potential
c
c               if (irepopt .eq. 0) then
c
c                 the necessary ingredients:
c
c            -----  # valence electrons
c
                  nami=nxcent(ii)
c
c            -----  check if the point represents a group
c
                  if (nami(:5) .eq. 'group') then
c           6-----
c              -----  non-bonded interactions evaluated only between
c                     internal atoms and individual members
c                     of a group:
c                     skip the group-representing point
c
                    goto 290
                  else
c
c              -----  get the number of valence electrons of the
c                     external point
c
                    call drfnval(nxcent(ii),nvali,znuc)
c
c              -----  check if it was already found
c                     close to another internal atom
c
c                   do 220, nsp = 1, nspec
c                     if (nspecl(nsp) .eq. ii) goto 230
c 220               continue
c
c              -----  this external point has not been found
c                     close to an internal atom:
c                     add it to the list
c
c                   nspec = nspec + 1
c                   nvalel(nspec) = nvali
c                   polars(nspec) = alfa(nxcent(ii),0,ier)
c                   polars(nspec) = alfext(ii)
c                   nspecl(nspec) = ii
c 230               continue
c
c              -----  1/r12 term for selected atoms
c
c              -----  get polarizability of external atom
c
c                   if (ipol) then
c                     alfi = polar(imp)
c                   else
c                     alfi = alfa(nami,0,ier)
c                   endif
                    alfi = alfext(ii)
                                                
                    if (alfi .ne. zero) then
                      aoverni = sqrt(alfi/nvali)
                    endif
c
c              -----  calculate equilibrium lj distance from
c                     "size" of internal and external atom
c
c              -----  size of external atom: ri
c
                    ri = radext(ii)*rfact
c
c              -----  size of external atom: rjx
c
                    rjx = radat(j)*rfact
c
c              -----  in case of a h-bond, the sizes need to be
c                     modified
c
c              -----  check h-bond criterium
c
                    if ((ihbond .eq. 1) 
     1                  .and. (dist .le. hbondl)) then
                      if(namj(:2).eq.'h'.and.
     1                  (nami(:2).eq.'n'.or.
     2                   nami(:2).eq.'o'.or.
     3                   nami(:2).eq.'f')) rjx=hbondr
                      if(nami(:2).eq.'h'.and.
     1                  (namj(:2).eq.'n'.or.
     2                   namj(:2).eq.'o'.or.
     3                   namj(:2).eq.'f')) ri=hbondr
                    endif
cmas  changed feb 7, 1996
c
c              -----  calculate the 1/r12 lj term for the close contact
c
                    if (alfj .ne. zero .and. alfi .ne. zero) then
                      pseu = pseufac*alfj*alfi/(aovernj+aoverni)
                      pseu1 = (ri+rjx)**6
                      repmodg(igr,iactst) = repmodg(igr,iactst) +
     1                                      pseu*pseu1*(dmind1**12)
                      epseu = epseu + pseu*pseu1*(dmind1**12)
c
c                -----  and derivative
c
                      do 270, k = 1, 3
                        dpseu(k,j,igr) = dpseu(k,j,igr) +
     1                         twelve*pq(k)*pseu*pseu1*(dmind1**14)
  270                 continue
                    endif
c          6-----
                  endif
c               else if (irepopt .eq. 1) then
c
c            -----  set some necessary blabla for integrals
c
c                  call drfpseu
c
c               endif
c         5-----
              endif
c       4-----
            endif
  290       if (ipol) imp = imp + 1
c     3-----
  300     continue
c   2-----
  500   continue
c
c  -----  write zfa, zpa, dzfa and dpseu to disk.
c
        call dawrit(idafdrf,iodadrf,zfa,4*nzfa,3,navdrf)
        call dawrit(idafdrf,iodadrf,zpa,nzfa,4,navdrf)
        call dawrit(idafdrf,iodadrf,dzfa,12*nzfa,5,navdrf)
        call dawrit(idafdrf,iodadrf,dpseu,3*nzfa,6,navdrf)
c
        repmod(iactst) = epseu
c 1-----
      endif
c
      if (idrfout .eq. 5) then
        do 600, igr = 1, ngran
          call hatout(dzfa(1,igr),4,3*nexp,5,'dzfa')
  600   continue
      endif
c
      if (idrfout .eq. 3) then
        do 700, igr = 1, ngran
          call hatout(zfa(1,igr),4,nexp,5,'zfa')
          call hatout(zpa(1,igr),1,nexp,5,'zpa')
  700   continue
      endif
c
      return
      end
