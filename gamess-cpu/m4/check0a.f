cvp configuration comparison for dk=0
c   establishing the interacting configurations through the use of 
c   creation and annihilation operators
      subroutine check0a(twoe,nteint,iottr,iotnew,iotm,
     +   maindf,maxr,jconb,nomax2,
     +   ioint,isk,jsk,ibs,nspiel,mspiel,jnum1
     &  ,nspp1,nsppa,nsppb,nspp4)
cvp
      implicit real*8 (a-h,p-z)
cvp parameters and global common-blocks are in the include-file
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
      real*8 twoe
      integer nteint
      dimension twoe(nteint)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr,jconb,nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
cvp
c special for dk=0
c jmerk loops directly controlled by the loop counter
c because of this the access succeeds faster
c changed version: only the number of interacting configurations to one
c with same sk succeeds the comparison in relation to the lower 
c triangular matrix
c (vp 13.3.1992)
c  persistent configurations get established
c  passed on are: iswh=number of sk
c                 isk=left sk
c                 jsk=right sk
c                 ibs=number of configuration from the sk isk
c  results are:   nspiel=number of ww configurations.
c                 ispiel=vector, containing the number of ww 
c                                configurations
c
      ncisk=nconf(isk)
      inytl  = nytl(isk)
      indub  = ndub(isk)
      inopen = inytl - indub
      istar=niot(isk)+1
cvp
      iloc=istar+(ibs-1)*inytl
c do loop over the configurations with respect to the one that should
c be tested
c nullifying the jcon field
        do 215 k=1,imo
        jcon(k)=-1
        jcon1(k)=1
        jcon2(k)=2
215     continue
cvp storing the test configuration on itest
         do 205 k=1,inytl
  205    itest(k)=iottr(istar+(k-1)*ncisk+ibs-1)
         do 220 k=1,inopen
         jposo(itest(k))=k
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 1
         jcon(itest(k)) = 0
  220    continue
c
         do 221 k=inopen+1,inytl
c first occupying jcon1, such that starting from single excitations and
c finding double excitations do not get counted
         jposc(itest(k)) = k-inopen
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 0
         jcon (itest(k)) = 1
  221    continue
c
         if(jsk.lt.1) go to 250
         jnytl=nytl(jsk)
         jndub=ndub(jsk)
         jnopen=jnytl - jndub
         ncjsk = nconf (jsk)
c building the field ispiel (number of players) removing all that were
c identified already as double excitations
cvp comparison with respect to the lower triangular matrix
cvp nspiel is here the maximum number of configurations to be tested,
cvp  that fit onto the vectors
cvp  if all fit onto the vectors then nspiel=ibs-1
         nspiel=mspiel
c ispiel : number of configurations that play with sk
c nspiel : number of players within sk jsk
c
c establishing the interacting configurations with the help of 
c  creation and annihilation operators
*     call chkww(iotnew,iotm,maindf,maxr,jconb,nomax2,
*    +           isk,jsk,ibs,nspiel,jnum1,ncisk,ncjsk)
c aufbau von jcona for vergleich auf basis von erzeugern
c build jcona for comparison with respect to generators
         do k=1,imo
c  for refcon bzw. refcon1
           jcona(k)=jcon1(k)
           jcona(k+momax)=jcon2(k)-jcon1(k)
         enddo
c
c aufbau von idref, d.h. vergleich der testkonfiguration ibs mit
c allen mains
c build idref, i.e. compare the test configuration ibs to all mains
c    berechnung der startadresse iostt for das iotnew-feld
c    calculate the start address iostt for the iotnew field
       iostt=iotnst(isk)
c    vorbelegung von idref
c    initialisation of idref
       mainr=iotnew(iostt+ibs)
       do k=1,nko
         idref(k)=maindf(k,mainr)
       enddo
c    korrektur von idref um die beitraege von jcona
c    correct idref with the contributions of jcona
       moai=iotnew(iostt+  ncisk+ibs)
       moaj=iotnew(iostt+2*ncisk+ibs)
       moak=iotnew(iostt+3*ncisk+ibs)
       moal=iotnew(iostt+4*ncisk+ibs)
       do k=1,nko
         idref(k)=idref(k)+jconb(k,moai)+jconb(k,moaj)
     &                    -jconb(k,moak)-jconb(k,moal)
       enddo
c
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nspiel : number of players within the sk jsk
c
c
c berechnung von idiff
c calculation of idiff
c    berechnung der startadresse iosttj for das iotnew-feld
c    calculate the start address iosttj for the iotnew field
         iosttj =iotnst(jsk)+jnum1
         iostj1=iosttj + ncjsk
         iostj2=iostj1 + ncjsk
         iostj3=iostj2 + ncjsk
         iostj4=iostj3 + ncjsk
         do kdoo=1,nspiel
           idiff(kdoo)=idref( iotnew(iosttj+kdoo) )
     &                -jcona( iotnew(iostj3+kdoo) )
     &                -jcona( iotnew(iostj4+kdoo) )
         enddo
         do kdoo=1,nspiel
          if (idiff(kdoo).le.2) then
           idiff(kdoo)=idiff(kdoo)
     &                 +jcona( iotnew(iostj1+kdoo) )
     &                 +jcona( iotnew(iostj2+kdoo) )
          endif
         enddo
c
cbs
cbs now compress
cbs
          nspie1=0
*vocl loop,novrec
          do 385 kdoo = 1,nspiel
           if(idiff(kdoo).le.2) then
             nspie1 = nspie1 + 1
             ispiel(nspie1) = kdoo + jnum1
             idiff(nspie1)=idiff(kdoo)
           endif
 385      continue
cvp jstar is set back to its old value
         jstar=niot(jsk)
cvp
cvp  end of finding the interacting configurations
cvp
          nspiel=nspie1
cbs   and now the whole thing once again
cvp   initialising
           do 261 jbs=1,nspiel
 261       jmerko(jbs)=0
cvp loop over all open shells for the cause of p-case separation
cvp
c test of open shells
c as with the closed shells start at the last open shell
         kstar=jstar+jnopen*ncjsk
         do 3050 jdoo = 1,jnopen
          kstar=kstar-ncjsk
c establish the difference to jcon
*vocl loop,novrec
          do 3040 kdoo = 1,nspiel
            jmerko(kdoo)=jmerko(kdoo)+abs( jcon( iottr(kstar+
     &                ispiel(kdoo)) ) )
 3040     continue
 3050    continue
cvp establishing the p-case using jmerko and idiff
c notation after buenker, i.e. npfal=2 --> p=1 etc.
          nspp1=0
          nspp2=0
          nspp3=0
          nspp4=0
      do 3100 kdoo=1,nspiel
       if (jmerko(kdoo).eq.2) then
         nspp1=nspp1+1
         nwwmo(nspp1,1)=ispiel(kdoo)
cvp      npfal(kdoo)=2
cvp    endif
        else if (jmerko(kdoo).eq.1.and.idiff(kdoo).eq.2) then
         nspp2=nspp2+1
         nwwmo(nspp2,2)=ispiel(kdoo)
cvp      npfal(kdoo)=3
cvp    endif
        else if (idiff(kdoo).eq.1) then
         nspp3=nspp3+1
         nwwmo(nspp3,3)=ispiel(kdoo)
cvp      npfal(kdoo)=4
cvp    endif
cvp    if (jmerko(kdoo).eq.0) then
        else
         nspp4=nspp4+1
         nwwmo(nspp4,4)=ispiel(kdoo)
cvp      npfal(kdoo)=5
       endif
 3100 continue
cvp sorted configuration numbers according to p-cases --> npfal
      do 3210 kdoo=1,nspp1
        npfal(kdoo)=nwwmo(kdoo,1)
        jmerko(kdoo)=2
        jmerkd(kdoo)=0
        nqrfal(kdoo)=ibinom(inopen,2)
 3210 continue
      do 3220 kdoo=1,nspp2
        npfal(kdoo+nspp1)=nwwmo(kdoo,2)
        jmerkd(kdoo+nspp1)=0
        nqrfal(kdoo+nspp1)=ideks(inopen+1)
        nwwmo(kdoo+nspp1,1)=ideks(indub+1)
 3220 continue
      nsppa=nspp1+nspp2
      do 3230 kdoo=1,nspp3
        npfal(kdoo+nsppa)=nwwmo(kdoo,3)
        jmerko(kdoo+nsppa)=0
        jmerkd(kdoo+nsppa)=0
        nqrfal(kdoo+nsppa)=ideks(inopen+1)
        sumint(kdoo+nsppa)=0.0d0
 3230 continue
      nsppb=nsppa+nspp3
      do 3240 kdoo=1,nspp4
        npfal(kdoo+nsppb)=nwwmo(kdoo,4)
        moafal(kdoo+nsppb)=ideks(indub+1)
 3240 continue
cvp loop over all doubly occupied mo's of phi_l
c        write(6,*) 'nspiel = ',nspiel
c start at the back, because there are the most differences
         kstar=jstar+jnytl*ncjsk
         do 1260 jdoo = 1,jndub
          kstar=kstar-ncjsk
cvp loop for p=1
cvp jmerkd counts the occupation differences with respect to the doubly
c   occupied mo's of \phi_l
*vocl loop,novrec
          do 281 kdoo=1,nspp1
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              jmerkd(kdoo)=jmerkd(kdoo)+1
              npos(kdoo,1)=jposo(iottr(kstar+npfal(kdoo)))
            endif
 281      continue
cvp loop for p=2
*vocl loop,novrec
          do 282 kdoo=nspp1+1,nsppa
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              nwwmo(kdoo,2)=iottr(kstar+npfal(kdoo))
cvp jmerkd(*)=1: r=1
ct            if (jdiff(kdoo).eq.-2) jmerkd(kdoo)=1
              if (jdiff(kdoo).eq.2) jmerkd(kdoo)=1
             else
              nx=jposc(iottr(kstar+npfal(kdoo)))
              nwwmo(kdoo,1)=nwwmo(kdoo,1)-nx
            endif
 282      continue
cvp loop for p=3
*vocl loop,novrec
          do 283 kdoo=nsppa+1,nsppb
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              jmerkd(kdoo)=jmerkd(kdoo)+1
            endif
 283      continue
cvp loop for p=4
cvp jmerkd counts the occupation differences with respect to the doubly
c   occupied mo's of \phi_l
*vocl loop,novrec
          do 284 kdoo=nsppb+1,nsppb+nspp4
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              mobfal(kdoo)=iottr(kstar+npfal(kdoo))
             else
              nx=jposc(iottr(kstar+npfal(kdoo)))
              moafal(kdoo)=moafal(kdoo)-nx
          endif
 284      continue
1260     continue
cvp
cvp loop over the open shells starting from the back
cvp
         kstar=jstar+jnopen*ncjsk
         do 1360 jdoo = 1,jnopen
          kstar=kstar-ncjsk
c         write(6,*) 'open shell nr. ',jdoo,
c    *'jstar  ',jstar
c establish the difference to jcon
cvp jmerko counts the occupation differences with respect to the singly
c   occupied mo's of \phi_l, is counted down to zero again here.
cvp p=1-fall
*vocl loop,novrec
          do 481 kdoo=1,nspp1
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
          if(jdiff(kdoo).ne.0) then
            jmerko(kdoo)=jmerko(kdoo)-1
            npos(kdoo,2)=jdiff(kdoo)
ct       nasty: reconsider r-case determination 
ct        or define jcon = -1 for unoccupied mo's and include
ct        abs in establishing jmerko ?
ct          npos(kdoo,2)=1-jcon2( iottr(kstar+npfal(kdoo)) )
            nwwmo(kdoo,4-jmerko(kdoo))=inopen+1-jdoo
           else
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-
     &       ibinom(nx-1,inopen+1-jdoo-jmerko(kdoo))
          endif
481       continue
cvp p=2-fall
*vocl loop,novrec
          do 482 kdoo=nspp1+1,nsppa
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
          if(jdiff(kdoo).ne.0) then
            nqlfal(kdoo)=inopen-jdoo+1
           else
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-nx
          endif
482       continue
cvp p=3-fall
*vocl loop,novrec
          do 483 kdoo=nsppa+1,nsppb
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
          if(jdiff(kdoo).ne.0) then
            nqlfal(kdoo)=inopen-jdoo+1
            mobfal(kdoo)=iottr(kstar+npfal(kdoo))
           else
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-nx
          endif
483       continue
 1360   continue
cvp determination of r-case for p=1
      do 491 kdoo=1,nspp1
          if (jmerkd(kdoo).eq.1) then
cvp        dopppelt besetztes mo von phi_l trifft auf das hoehere
cvp         einfach besetzte von phi_r
            if (npos(kdoo,1).eq.imoq(nqrfal(kdoo),1)) then
                npos(kdoo,1)=-1
              else
                npos(kdoo,1)=1
            endif
cvp          npos(*,1) is now the relative occupation of the last open
cvp           mo's starting from phi_r
            if (npos(kdoo,1).eq.npos(kdoo,2)) then
              nrfal(kdoo)=2
             else
              nrfal(kdoo)=3
            endif
           else
            nrfal(kdoo)=1
          endif
  491 continue
cvp
cvp p=1
cvp determination of ql-case from the ww mo's of \phi_l
      do 581 kdoo=1,nspp1
        nqlfal(kdoo)=iqmo(nwwmo(kdoo,3),nwwmo(kdoo,4))
  581 continue
cvp p=3
*vocl loop,novrec
      do 493 kdoo=nsppa+1,nsppb
cvp determination of the irreducable representation and storing the
c   higher and lower mo-numbers within the irreducable representation.
cvp
cvp idiff: position of the open mo's of \phi_r --> needed for
cvp cb- and ex-integrals 
c       idiff(kdoo)=moafal(kdoo)
        idiff(kdoo)=nqrfal(kdoo)
cvp
        nrfal(kdoo)=nirred(mobfal(kdoo))
        mobfal(kdoo)=mopos(mobfal(kdoo))
        moafal(kdoo)=mopos(itest(nqrfal(kdoo)))
cvp moafal contains the bigger mo
        if (moafal(kdoo).lt.mobfal(kdoo)) then
          nx=moafal(kdoo)
          moafal(kdoo)=mobfal(kdoo)
          mobfal(kdoo)=nx
        endif
cvp coding the r-case for p=3 with the sign of qr
        if (jmerkd(kdoo).eq.1) nqrfal(kdoo)=-1*nqrfal(kdoo)
 493  continue
cvp
cvp
cvp determining the absolute mo-numbers for p=4
      do 494 kdoo=nsppb+1,nsppb+nspp4
          moafal(kdoo)=itest(inopen+moafal(kdoo))
c         mobfal(kdoo)=iot(iloc+inopen-1+mobfal(kdoo))
 494  continue
cvp
cvp  end of label determination
cvp
cvp no calculation of the integral addresses in case the integrals
cvp  are presorted already, i.e. ioint = 1
      if (ioint.eq.1) goto 251
cvp
cvp determining the ww mo's of \phi_r from qr,
cvp calculating and sorting the absolute mo-numbers for the integral
cvp calculation
         kstar = niot(jsk) - ncjsk
*vocl loop,novrec
      do 591 kdoo=1,nspp1
cvp consider the Coulomb-integral in the chemists notation (d,c | b,a),
cvp then the classification applies:
cvp                            d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp consider the exchange-integral in the chemists notation (d,c | b,a),
cvp then the classification applies:
cvp                            d = nwwmo(*,1)
cvp                            c = idiff(*)
cvp                            b = jdiff(*)
cvp                            a = jmerko(*)
cvp  idiff, jdiff and jmerko are freely applicable here
          nwwmo(kdoo,3)=iottr(kstar+npfal(kdoo)+nwwmo(kdoo,3)*ncjsk)
          nwwmo(kdoo,4)=iottr(kstar+npfal(kdoo)+nwwmo(kdoo,4)*ncjsk)
          nwwmo(kdoo,1)=itest(imoq(nqrfal(kdoo),1))
          nwwmo(kdoo,2)=itest(imoq(nqrfal(kdoo),2))
c       write(6,*) ' r-case, ww mo: ',nrfal(kdoo),(nwwmo(kdoo,k),k=1,4)
          if (nwwmo(kdoo,1).lt.nwwmo(kdoo,3)) then
            nx1=nwwmo(kdoo,1)
            ny1=nwwmo(kdoo,2)
            nwwmo(kdoo,1)=nwwmo(kdoo,3)
            nwwmo(kdoo,2)=nwwmo(kdoo,4)
            nwwmo(kdoo,3)=nx1
            nwwmo(kdoo,4)=ny1
          endif
cvp sorting the exchange-integrals
          if (nrfal(kdoo).eq.2) then
            idiff(kdoo)=nwwmo(kdoo,3)
            if (nwwmo(kdoo,2).gt.nwwmo(kdoo,4)) then
              jdiff(kdoo)=nwwmo(kdoo,2)
              jmerko(kdoo)=nwwmo(kdoo,4)
             else
              jdiff(kdoo)=nwwmo(kdoo,4)
              jmerko(kdoo)=nwwmo(kdoo,2)
            endif
           else
            idiff(kdoo)=nwwmo(kdoo,4)
            if (nwwmo(kdoo,2).gt.nwwmo(kdoo,3)) then
              jdiff(kdoo)=nwwmo(kdoo,2)
              jmerko(kdoo)=nwwmo(kdoo,3)
             else
              jdiff(kdoo)=nwwmo(kdoo,3)
              jmerko(kdoo)=nwwmo(kdoo,2)
            endif
          endif
cvp sorting the Coulomb-integrals
          if (jmerkd(kdoo).ne.1) then
            nxx=nwwmo(kdoo,2)
            nwwmo(kdoo,2)=nwwmo(kdoo,3)
            nwwmo(kdoo,3)=nxx
            if (nwwmo(kdoo,3).lt.nwwmo(kdoo,4)) then
              nyy=nwwmo(kdoo,4)
              nwwmo(kdoo,4)=nwwmo(kdoo,3)
              nwwmo(kdoo,3)=nyy
            endif
          endif
c       write(6,*) ' r-case, ww mo: ',nrfal(kdoo),(nwwmo(kdoo,k),k=1,4)
c     write(6,*)
  591 continue
cvp p=2 case
cvp calculation of the absolute mo-numbers and sorting
*vocl loop,novrec
      do 492 kdoo=nspp1+1,nsppa
        nwwmo(kdoo,1)=itest(inopen+nwwmo(kdoo,1))
        if (jmerkd(kdoo).eq.1) then
cvp r=1
          nwwmo(kdoo,3)=itest(nqrfal(kdoo))
          nx3=nwwmo(kdoo,1)
          nwwmo(kdoo,1)=nwwmo(kdoo,2)
          nwwmo(kdoo,2)=nx3
         else
cvp r=2
          nwwmo(kdoo,3)=iottr(kstar+npfal(kdoo)+nqlfal(kdoo)*ncjsk)
        endif
cvp mo's are sorted as follows:
cvp  nwwmo(*,1) = doubly occupied mo which is empty in the other 
cvp               configuration
cvp  nwwmo(*,2) = doubly occupied mo which is singly occupied in the 
cvp               other configuration
cvp  nwwmo(*,3) = singly occupied mo which is empty in the other 
cvp               configuration
  492 continue
cvp
c     do 3500 kdoo=nspp1+1,nsppa
c       write(6,*) ' phi_l: ',(iot(k),k=jstar+npfal(kdoo)*jnytl,
c    &   jstar+npfal(kdoo)*jnytl+jnytl-1)
c       write(6,*) ' ww mo: ',(nwwmo(kdoo,k),k=1,3)
c3500 continue
cvp
*vocl loop,novrec
      do 592 kdoo=nspp1+1,nsppa
        if (nwwmo(kdoo,2).lt.nwwmo(kdoo,3)) then
          nx4=nwwmo(kdoo,2)
          nwwmo(kdoo,2)=nwwmo(kdoo,3)
          nwwmo(kdoo,3)=nx4
        endif
  592 continue
cvp
cvp
cvp  calculation of the integral addresses
cvp
c     iloc1=iloc-1
*vocl loop,temp(ndum1,ndum2,ndum3,ndum4,idum1,idum2,idum3,idum4,idum)
      do 3201 kdoo=1,nspp1
cvp consider the Coulomb-integral in the chemists notation (d,c | b,a),
cvp then the classification applies:
cvp                            d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp access to the nit-field
        ndum1=nirred(nwwmo(kdoo,1))
        ndum2=nirred(nwwmo(kdoo,2))
        ndum3=nirred(nwwmo(kdoo,3))
        ndum4=nirred(nwwmo(kdoo,4))
cvp switch to mo-number with the corresponding irreducable 
cvp representation
        idum1=mopos(nwwmo(kdoo,1))
        idum2=mopos(nwwmo(kdoo,2))
        idum3=mopos(nwwmo(kdoo,3))
        idum4=mopos(nwwmo(kdoo,4))
        idum=ideks(ndum1)
cvp determination of the symmetry case (Coulomb-integral)
*vocl stmt,if(80)
        if (ndum1.eq.ndum2) then
c       if (ndum2.eq.ndum3) then
*vocl stmt,if(60)
          if (ndum1.eq.ndum3) then
cvp      case 1
c           write(6,*) ' case 1'
            intcb(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &       + idum4 +
     &       ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)
            intex(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &       + mopos(jmerko(kdoo)) +
     &       ideks(ideks(idum1)+
     &        mopos(idiff(kdoo)))
     &        +ideks(mopos(jdiff(kdoo)))
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
c          else if (ndum1.eq.ndum2) then
cvp      case 2
c           write(6,*) ' case 2'
            intcb(kdoo)=nit( ideks(ideks(ndum1+1))
     &       +ideks(ndum3+1) ) + idum4 +
     &       ideks(ncimo(ndum3)+1)*(ideks(idum1)+idum2-1)
     &       +ideks(idum3)
            intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo))+1) )
     &       + mopos(jmerko(kdoo)) +
     &       ideks( ncimo(nirred(idiff(kdoo)))*(idum1-1)
     &       +mopos(idiff(kdoo)) )
     &       +ncimo(nirred(idiff(kdoo)))
     &       *(mopos(jdiff(kdoo))-1)
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
         else
*vocl stmt,if(80)
          if (ndum1.eq.ndum3) then
c         else if (ndum1.eq.ndum3) then
cvp      case 3
c           write(6,*) ' case 3'
            intcb(kdoo)=nit( ideks(idum+ndum2+1) )
     &       +idum4 +
     &       ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1)
            if (nrfal(kdoo).ne.2) then
              intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo))+1) )
     &         + mopos(jmerko(kdoo)) +
     &         ideks( ncimo(nirred(idiff(kdoo)))*(idum1-1)
     &         +mopos(idiff(kdoo)) )
     &         +ncimo(nirred(idiff(kdoo)))
     &         *(mopos(jdiff(kdoo))-1)
             else
              intex(kdoo)=nit( ideks(ideks(ndum1+1))
     &         +ideks(nirred(jdiff(kdoo))+1) )
     &         + mopos(jmerko(kdoo)) +
     &         ideks(ncimo(nirred(jdiff(kdoo)))+1)
     &         *(ideks(idum1)+mopos(idiff(kdoo))-1)
     &         +ideks(mopos(jdiff(kdoo)))
            endif
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
cvp      case 4
c           write(6,*) ' case 4'
cvp only i-1 shows up in the equations
            idum1=idum1-1
            intcb(kdoo)=nit( ideks(idum+ndum2)
     &       +ideks(ndum3)+ndum4 ) + idum4 +
     &       ncimo(ndum4)*(idum3-1+ncimo(ndum3)*(idum2-1+
     &       ncimo(ndum2)*idum1))
            intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo)))
     &       +ideks(nirred(jdiff(kdoo)))+nirred(jmerko(kdoo)) )
     &       + mopos(jmerko(kdoo)) +
     &       ncimo(nirred(jmerko(kdoo)))
     &       *(mopos(jdiff(kdoo))-1
     &       +ncimo(nirred(jdiff(kdoo)))
     &       *(mopos(idiff(kdoo))-1
     &       +ncimo(nirred(idiff(kdoo)))*idum1))
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
        endif
 3201 continue
cvp calculation of the integral addresses for p=2
*vocl loop,novrec
      do 3202 kdoo=nspp1+1,nsppa
cvp das integral sei (ab|ac) oder (ac|ab) usw.
cvp dann ist nwwmo(*,1)=a
cvp          nwwmo(*,2)=max(b,c)
cvp          nwwmo(*,3)=min(b,c)
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp if the integral is (ab|ac) or (ac|ab) etc.
cvp then     nwwmo(*,1)=a
cvp          nwwmo(*,2)=max(b,c)
cvp          nwwmo(*,3)=min(b,c)
cvp however it should still be sorted according to size
cvp access to nit field
        ndum1=nirred(max(nwwmo(kdoo,1),nwwmo(kdoo,2)))
        ndum2=nirred(min(nwwmo(kdoo,1),nwwmo(kdoo,2)))
c       ndum3=nirred(max(nwwmo(kdoo,1),nwwmo(kdoo,3)))
c       ndum4=nirred(min(nwwmo(kdoo,1),nwwmo(kdoo,3)))
cvp uebergang zu mo-nummern innerhalb der jeweiligen irred. darst.
cvp change over to mo-numbers within the current irrep.
        idum1=mopos(max(nwwmo(kdoo,1),nwwmo(kdoo,2)))
        idum2=mopos(min(nwwmo(kdoo,1),nwwmo(kdoo,2)))
        idum3=mopos(max(nwwmo(kdoo,1),nwwmo(kdoo,3)))
        idum4=mopos(min(nwwmo(kdoo,1),nwwmo(kdoo,3)))
c       idum=ideks(ndum1)
c       intcb(kdoo)=nit( ideks(idum+ndum2)
c    &    +ideks(ndum3)+ndum4 ) + idum4
cvp bestimmung des symmetrie-falles
cvp determining the symmetry case
        if (ndum1.eq.ndum2) then
cvp        fall 1
cvp        case 1
c           write(6,*) ' case 1'
        intcb(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &     + idum4 +
c           intcb(kdoo)=intcb(kdoo)+
     &       ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)
          else
cvp        fall 3
cvp        case 3
c           write(6,*) ' case 3'
        intcb(kdoo)=nit( ideks(ideks(ndum1)+ndum2+1) )
     &     + idum4 +
c           intcb(kdoo)=intcb(kdoo)+
     &       ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1)
        endif
 3202 continue
cvp
cvp  p=3
cvp
cvp  berechnung der summe uber cb- und ex-integrale sowie
cvp  der adressen der ex-integrale bei gemeinsamen offenen
cvp  schalen
cvp  calculate the sum over cb- and ex-integrals and also 
cvp  the addresses of the ex-integrals in case of common open shells
cvp
cvp  sei d gemeinsames offenes mo: betrachtung der integrale
cvp     (ab|dd) und (ad|bd), es gilt a > b
cvp  if d is a common open mo: considering the integrals
cvp     (ab|dd) and (ad|bd), where a > b
cvp
      do 593 jdoo=1,inopen-1
        do 603 kdoo=nsppa+1,nsppb
          mogd=jdoo+jmerko(kdoo)
          if (mogd.eq.idiff(kdoo)) jmerko(kdoo)=1
  603   continue
        do 613 kdoo=nsppa+1,nsppb
          mogd=itest(jdoo+jmerko(kdoo))
cvp das ex-integral sei (ad|bd), a>b
cvp dann ist a=moafal(*)
cvp          b=mobfal(*)
cvp          d=mogd
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp if the ex-integral is (ad|bd), a>b
cvp then     a=moafal(*)
cvp          b=mobfal(*)
cvp          d=mogd
cvp it still has to be sorted according to size
cvp access to nit field
        ndum1=nrfal(kdoo)
        ndum2=nirred(mogd)
cvp bestimmung des symmetrie-falles
cvp determining the symmetrie case
        if (ndum1.eq.ndum2) then
           idum1=max(moafal(kdoo),mopos(mogd))
           idum2=min(moafal(kdoo),mopos(mogd))
           idum3=max(mobfal(kdoo),mopos(mogd))
           idum4=min(mobfal(kdoo),mopos(mogd))
cvp        fall 1
cvp        case 1
c           write(6,*) ' case 1'
            nwwmo(kdoo,jdoo)=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)+idum4
            ix=ideks(mopos(mogd)+1)
            iy=ideks(moafal(kdoo))+mobfal(kdoo)
            icb=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks( max(ix,iy) )
     &       +       min(ix,iy)
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)
          else if (ndum1.gt.ndum2) then
           idum1=moafal(kdoo)
           idum2=mopos(mogd)
           idum3=mobfal(kdoo)
cvp        case 3
c           write(6,*) ' case 3'
            nwwmo(kdoo,jdoo)=nit( ideks(ideks(ndum1)+ndum2+1) )
     &       +ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1) + idum2
cvp        fall 2
cvp        case 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum1+1))+ideks(ndum2+1) )
     &       +ideks(ncimo(ndum2)+1)
     &       *(ideks(idum1)+idum3-1)+ideks(idum2+1)
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)
          else if (ndum1.lt.ndum2) then
           idum2=moafal(kdoo)
           idum1=mopos(mogd)
           idum4=mobfal(kdoo)
cvp        fall 3
cvp        case 3
c           write(6,*) ' case 3'
            nwwmo(kdoo,jdoo)=nit( ideks(ideks(ndum2)+ndum1+1) )
     &       +ideks(ncimo(ndum1)*(idum1-1)+idum2)
     &       +ncimo(ndum1)*(idum1-1) + idum4
cvp        fall 2
cvp        case 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum2+1))+ideks(ndum1+1) )
     &       +ideks(ncimo(ndum1)+1)
     &       *(ideks(idum1+1)-1)+ideks(idum2)+idum4
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)
        endif
  613   continue
  593 continue
cvp schleife ueber die (gemeinsamen) geschlossenen schalen
cvp loop over the (common) closed shells
      do 533 jdoo=1,indub
        do 623 kdoo=nsppa+1,nsppb
          mogc=itest(inopen+jdoo)
cvp das ex-integral sei (ac|bc), a>b
cvp dann ist a=moafal(*)
cvp          b=mobfal(*)
cvp          c=mogc
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp if the ex-integral is (ac|bc), a>b
cvp then     a=moafal(*)
cvp          b=mobfal(*)
cvp          c=mogc
cvp it should still be sorted according to size
cvp access to nit field
        ndum1=nrfal(kdoo)
        ndum2=nirred(mogc)
cvp bestimmung des symmetrie-falles
cvp determining the symmetry case
        if (ndum1.eq.ndum2) then
           idum1=max(moafal(kdoo),mopos(mogc))
           idum2=min(moafal(kdoo),mopos(mogc))
           idum3=max(mobfal(kdoo),mopos(mogc))
           idum4=min(mobfal(kdoo),mopos(mogc))
cvp        case 1
c           write(6,*) ' case 1'
            iex=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)+idum4
            ix=ideks(mopos(mogc)+1)
            iy=ideks(moafal(kdoo))+mobfal(kdoo)
            icb=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks( max(ix,iy) )
     &       +       min(ix,iy)
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)+twoe(icb)-twoe(iex)
          else if (ndum1.gt.ndum2) then
           idum1=moafal(kdoo)
           idum2=mopos(mogc)
           idum3=mobfal(kdoo)
cvp        case 3
c           write(6,*) ' case 3'
            iex=nit( ideks(ideks(ndum1)+ndum2+1) )
     &       +ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1) + idum2
cvp        case 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum1+1))+ideks(ndum2+1) )
     &       +ideks(ncimo(ndum2)+1)
     &       *(ideks(idum1)+idum3-1)+ideks(idum2+1)
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)+twoe(icb)-twoe(iex)
          else if (ndum1.lt.ndum2) then
           idum2=moafal(kdoo)
           idum1=mopos(mogc)
           idum4=mobfal(kdoo)
cvp        fall 3
cvp        case 3
c           write(6,*) ' case 3'
            iex=nit( ideks(ideks(ndum2)+ndum1+1) )
     &       +ideks(ncimo(ndum1)*(idum1-1)+idum2)
     &       +ncimo(ndum1)*(idum1-1) + idum4
cvp        fall 2
cvp        case 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum2+1))+ideks(ndum1+1) )
     &       +ideks(ncimo(ndum1)+1)
     &       *(ideks(idum1+1)-1)+ideks(idum2)+idum4
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)+twoe(icb)-twoe(iex)
        endif
  623   continue
  533 continue
cvp korrektur der integralsumme for r=2
cvp correction of the integral sum for r=2
        do 633 kdoo=nsppa+1,nsppb
          if (nqrfal(kdoo).lt.0) then
cvp zugriff auf nit-feld
cvp access to nit field
            ndum1=nrfal(kdoo)
            idum1=itest(idiff(kdoo))
            ix=ideks(mopos(idum1)+1)
            iy=ideks(moafal(kdoo))+mobfal(kdoo)
            icb=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks( max(ix,iy) )
     &       +       min(ix,iy)
            sumint(kdoo)=sumint(kdoo)+twoe(icb)
          endif
  633   continue
cvp
cvp  ende berechnung der integraladressen
cvp  finished calculating the integral addresses
cvp
 251  continue
cvp
c     if (ibs.eq.4.or.ibs.eq.42) then
c     do 3500 kdoo=1,nspp1
c       write(6,*) ' phi_l: ',(iot(k),k=jstar+npfal(kdoo)*jnytl,
c    &   jstar+npfal(kdoo)*jnytl+jnytl-1)
c       write(6,*) ' ww mo: ',(nwwmo(kdoo,k),k=1,4)
c       write(6,*) ' intcb= ',intcb(kdoo)
c3500 continue
c     endif
cvp
cvp ruecksortierung wird uebergangen !!!!!!!!!!!!!!!!!!!!!!!!!!
cvp  ispiel enthaelt jetzt die adressen der wechselwirkenden
cvp  konfigurationen: zuerst die for p=1, dann for p=2 usw.
cvp  auf npfal wird vorlaeufig der p-fall abgespeichert
cvp  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cvp back-sort is skipped !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cvp  ispiel now holds the addresses of the interacting 
cvp  configurations: first those for p=1, then for p=2 usw.
cvp  provisionaly the p-case is stored on npfal
cvp  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do 5210 kdoo=1,nspp1
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=2
 5210 continue
      do 5220 kdoo=nspp1+1,nsppa
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=3
 5220 continue
      do 5230 kdoo=nsppa+1,nsppb
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=4
 5230 continue
      do 5240 kdoo=nsppb+1,nsppb+nspp4
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=5
 5240 continue
cvp
cvp p-faelle werden in urspruengliche reihenfolge gebracht
cvp p-cases being put in original order
cvp
cre    lspp1=1
cre    lspp2=1
cre    lspp3=1
cre    lspp4=1
cvp
c      write(6,*) ' jbs=,nspp3=,ispp3= ',jbs,nspp3,(ispp3(kkkk),
c    &   kkkk=1,nspp3)
cvp
cre    do 5100 kdoo=1,nspiel
cre      if (ispiel(kdoo).eq.npfal(lspp1).and.lspp1.le.nspp1) then
cre        ndfal1(kdoo)=2
cre        ndfal2(kdoo)=nrfal(lspp1)
cre        ndfal3(kdoo)=nqlfal(lspp1)
cre        ndfal4(kdoo)=nqrfal(lspp1)
cre        ndfal7(kdoo)=intcb(lspp1)
cre        ndfal8(kdoo)=intex(lspp1)
cre        lspp1=lspp1+1
cre      endif
cre      if (ispiel(kdoo).eq.npfal(lspp2+nspp1).and.lspp2.le.nspp2)
cre  &  then
cre        ndfal1(kdoo)=3
cre        ndfal3(kdoo)=nqlfal(lspp2+nspp1)
cre        ndfal4(kdoo)=nqrfal(lspp2+nspp1)
cre        ndfal7(kdoo)=intcb(lspp2+nspp1)
cre        lspp2=lspp2+1
cre      endif
cre      if (ispiel(kdoo).eq.npfal(lspp3+nsppa).and.lspp3.le.nspp3)
cre  &  then
c          write(6,*) ' kdoo=,ispiel=,ispp3= ',kdoo,ispiel(kdoo)
c    &       ,ispp3(lspp3)
cre        ndfal1(kdoo)=4
cre        ndfal2(kdoo)=nrfal(lspp3+nsppa)
cre        ndfal3(kdoo)=nqlfal(lspp3+nsppa)
cre        ndfal4(kdoo)=nqrfal(lspp3+nsppa)
cre        ndfal5(kdoo)=moafal(lspp3+nsppa)
cre        ndfal6(kdoo)=mobfal(lspp3+nsppa)
cre        ndfeld(kdoo,1)=nwwmo(lspp3+nsppa,1)
cre        ndfeld(kdoo,2)=nwwmo(lspp3+nsppa,2)
cre        ndfeld(kdoo,3)=nwwmo(lspp3+nsppa,3)
cre        ndfeld(kdoo,4)=nwwmo(lspp3+nsppa,4)
cre        ndfeld(kdoo,5)=nwwmo(lspp3+nsppa,5)
cre        ndfeld(kdoo,6)=nwwmo(lspp3+nsppa,6)
cre        ndfeld(kdoo,7)=nwwmo(lspp3+nsppa,7)
cre        ndfeld(kdoo,8)=nwwmo(lspp3+nsppa,8)
cre        ndfeld(kdoo,9)=nwwmo(lspp3+nsppa,9)
cre        rdum(kdoo)=sumint(lspp3+nsppa)
cre        lspp3=lspp3+1
cre      endif
cre      if (ispiel(kdoo).eq.npfal(lspp4+nsppb).and.lspp4.le.nspp4)
cre  &  then
cre        ndfal1(kdoo)=5
cre        ndfal5(kdoo)=moafal(lspp4+nsppb)
cre        ndfal6(kdoo)=mobfal(lspp4+nsppb)
cre        lspp4=lspp4+1
cre      endif
c5100  continue
 250   continue
      return
      end
