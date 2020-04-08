c 
c  $Author: jmht $
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/newmrd4.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   refcon = newmrd4    =
c ******************************************************
c ******************************************************
c
cvp
cvp building fields that should facilitate comparing configurations
      subroutine chfeldr(iottr,iot,iotm)
cvp
cvp parameter and global common-blocks are in the include file
      implicit real*8 (a-h,p-z)
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c NB 
c sneak() field for sneaky references
c nsneak: # of sneaky references
      integer sneak(1000),nsneak
      common /cepa_sn1/sneak,nsneak
c
c modoc is a doubly occupied orbital in the ref. space
c docnum is the number doc. MO in the references
c doccnt is counter for doubly occ. MOs in the ref. space
c doctest is current doc. MO in the first ref.
c jjopen is the numer of open shells due to the state multiplicity
c
      integer doccnt,modoc(maxshl),docnum,doctest,jjopen 
      integer movac(256*maxshl),nvacuum
      common/cepa_sn2/doccnt,modoc,docnum,doctest,movac,nvacuum,jjopen
c
c cepai: if cepa calculation is requested
c icepa: specifies the cepa variant
c ipcepa: print level
c 
      logical cepai
      integer icepa,ipcepa
      common/cepa_mrd1/cepai,icepa,ipcepa
c
c NB
      integer iot, iottr, iotm
      dimension iot(iotm), iottr(iotm)
cvp
cvp  iot is a dummy field here used to transpose the configuration 
cvp  field iottr and next copy it back onto iottr
cvp
cvp
      do 10 i=1,4
        do 10 j=1,3
   10     idk2r(i,j)=0
      do 20 i=2,4
        do 20 j=1,i-1
   20     idk2r(i,j)=i-j
      idk2r(3,2)=3
cvp
      nra(1)=4
      nra(2)=4
      nra(3)=2
      nrb(1)=2
      nrb(2)=3
      nrb(3)=3
      nrc(1)=3
      nrc(2)=2
      nrc(3)=4
cvp
cvp building of the ibinom field
cvp ibinom contains the binomial coefficients
cvp (i uber 0) = 1
      do 110 i=0,nopmax+1
  110 ibinom(i,0)=1
      do 115 i=1,nopmax+1
      do 115 j=1,nopmax+1
  115 ibinom(i,j)=0
cvp (0 uber 1) = 0
      do 120 i=0,nopmax+1
  120 ibinom(i,1)=i
      do 130 i=2,nopmax+1
      do 130 j=2,i
  130 ibinom(i,j)=ibinom(i-1,j)+ibinom(i-1,j-1)
cvp
cvp storage of the q-case classification <--> mo-numbers for dk=0,p=1
cvp iqmo(i,j): i=higher mo, j=lower mo; From the mo-numbers follows
cvp the q-case:
cvp imoq(q-fall,1): higher mo
cvp imoq(q-fall,2): lower mo
cvp
      do 140 i=2,nopmax
      do 140 j=1,i-1
      iqmo(i,j)=ibinom(i,2)-i+j+1
      imoq(iqmo(i,j),1)=i
  140 imoq(iqmo(i,j),2)=j
cvp
cvp storage of the q-case classification <--> mo-numbers for dk=1,p=1
cvp imo3q(q-fall,1): highest mo
cvp imo3q(q-fall,2): middle mo
cvp imo3q(q-fall,3): lowest mo
      do 160 i=3,nopmax
      do 160 j=2,i-1
      do 160 k=1,j-1
      ndum=ibinom(i-1,3)+ibinom(j-1,2)+k
      imo3q(ndum,1)=i
      imo3q(ndum,2)=j
  160 imo3q(ndum,3)=k
cvp
cvp storage of the q-case classification <--> mo-numbers fuer dk=2
cvp imo4q(q-fall,1): highest mo
cvp imo4q(q-fall,2): second highest mo
cvp imo4q(q-fall,3): third highest mo
cvp imo4q(q-fall,4): lowest mo
      do 170 i=4,nopmax
      do 170 j=3,i-1
      do 170 k=2,j-1
      do 170 l=1,k-1
      ndum=ibinom(i-1,4)+ibinom(j-1,3)+ibinom(k-1,2)+l
      imo4q(ndum,1)=i
      imo4q(ndum,2)=j
      imo4q(ndum,3)=k
  170 imo4q(ndum,4)=l
c NB1
      if (cepai) then
         kcount = 0
         nsneak = 0
c
         if (ipcepa.ge.1) then
c
      write(iwr,'(//" symmetry and spin adapted",
     & " CAS reference configurations",/
     & " surviving the configuration selection",/
     & " per super-category",/)')
c
         endif
      end if
c NB1
cvp
cvp  transposising the iot-field
cvp
cx    write(iwr,*) ' printing out all configurations '
      do 200 i=1,iswhm
        istar=niot(i)
        do 210 j=1,nytl(i)
          do 220 k=1,nconf(i)
            iot(istar+k+(j-1)*nconf(i))=
     &           iottr(istar+j+(k-1)*nytl(i))
  220     continue
  210   continue
c
c NB2
c
      if(cepai) then
c jj runs over conf. # in super-category i
c kdoc runs over doc.MOs
c kconf runs over conf.MOs
c ii label for MO # in iottr()
c kcount sums all ref.# up to sk i
c
      if (ipcepa.ge.1) then
         write(iwr,'(/"super-category: ",i2,/)')i
      endif
c
         do 211 jj=1,nconf(i)
            doccnt = 0
            do kdoc=1,docnum
               kk = 2*(i-1)+1+jjopen
               do 225 kconf=kk,nytl(i)
                  ii =istar+kconf+(jj-1)*nytl(i)
                  if(modoc(kdoc).ne.iottr(ii)) go to 225
                  doccnt = doccnt + 1
 225           continue
            enddo
c
c kvac runs over vac.MOs 
c
         if(doccnt.eq.docnum) then 
c           write(iwr,*)'all doc.MOs found' 
            doccnt = 0 
            do kconf=1,nytl(i) 
               do 226 kvac=1,nvacuum
                  ii =istar+kconf+(jj-1)*nytl(i)
                  if(movac(kvac).eq.iottr(ii)) doccnt=doccnt+1
 226           continue
            enddo
            if(doccnt.eq.nytl(i)) then
              nsneak = nsneak + 1
              sneak(nsneak) = jj + kcount
c             write(iwr,*)'sneaky ref.',nsneak,sneak(nsneak)
c
              write(iwr,'(i4,1x,24i4)')sneak(nsneak),
     &               (iottr(istar+ii+(jj-1)*nytl(i))
     &                                 ,ii=1,nytl(i) )
c
            endif
         end if
c        do 211 jj=1,nconf(i)
c         write(iwr,'("konf. nr.",i4,1x,24i4)')jj,
c    &               (iottr(istar+ii+(jj-1)*nytl(i))
c    &                                 ,ii=1,nytl(i) )
  211    continue
      kcount = kcount + nconf(i)
c     write(iwr,'("# of conf in sk",i2,"=",i8)')i,kcount
c
      end if
c
c NB2 
c
  200 continue
c
c NB3
      if (ipcepa.ge.1) then
      write(iwr,'("  number of survived CAS reference configurations: "
     & ,i4)')nsneak
      end if
c NB3
cx    write(iwr,*) ' end of print '
cvp restore onto iottr
      do 230 i=1,iotm
  230 iottr(i)=iot(i)
cvp
      return
      end
cvp
      subroutine compar(iotnew,iottr,iotm,
     +    mcsum,jsk,jref,nspiel,mspiel)
cvp
      implicit real*8 (a-h,p-z)
      parameter (iwod=1000)
cvp parameters and global common-blocks are in include-file
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer iotnew, iottr, iotm
      dimension iotnew(iotm), iottr(iotm)
cvp
c new fields: inop = # open shells per sk
c             iskop = sk for a given # open shells
c             iref: iref(0,*) = sk of the reference, followed by the
c                   singly and doubly occupied mo's
      integer inop,iskop
      dimension inop(iswhm),iskop(0:nopmax)
      common/mains/inop,iskop
c
c iotnew: new iot-field, contains the number of the generating reference
c         and the creators and annihilators.
c         max. size: 5 * (# konfigurationen)
      integer ipoint(ncmax)
      common/iotn/ipoint
c
c konflg: configurations in the long form to establish the annihilators
c rwork is 'scratch-common'
      integer konflg(momax*ncmax)
      common/rwork/konflg
c
c  establishing the permanent configurations
c  passed down: iswh=number of sk
c               isk=left sk
c               jsk=right sk
c               ibs=number of configurations from the sk isk
c  results are: nspiel=number of ww konf.
c               ispiel=vector, holding the number of ww configurations.
cvp
      isk=iref(0,jref)
cvp
        inytl  = nytl(isk)
        indub  = ndub(isk)
        inopen = inytl - indub
c nullifying the jcon field
c   imo = number of electrons, read from ft31, therefore here
        imo=momax
        do 215 k=1,imo
        jcon(k)=1
        jcon2(k)=2
        jcon1(k)=1
215     continue
cvp storing the test configuration on itest
         do 205 k=1,inytl
  205    itest(k)=iref(k,jref)
cx       write(iwr,*) ' itest= ',(itest(ii),ii=1,inytl)
         do 220 k=1,inopen
         jposo(itest(k))=k
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 1
         jcon(itest(k)) = 0
  220    continue
c
         do 221 k=inopen+1,inytl
c first occupy jcon1 such that finding doubles starting from singles
c doesn't get counted.
         jposc(itest(k)) = k-inopen
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 0
  221    continue
cvp
         jnytl=nytl(jsk)
         jndub=ndub(jsk)
         jnopen=jnytl - jndub
         ncjsk = nconf (jsk)
c building field ispiel (# of members) removing all configurations
c that were already identified as doubly excited.
cvp compare against lower triangle matrix
cvp nspiel is here the max. # configurations to be tested that fit on
cvp the vectors
         nspiel=mspiel
c ispiel : number of configurations that play within sk
c nspiel : number of players within sk jsk
c nullifying idiff as neccesary
         do 255 ialex = 1,nspiel
         ispiel(ialex) = 0
255      idiff(ialex) = 0
c do loop over open shells
c start at the end because that is were the most differences are.
         jstar=niot(jsk)
c
         kstar=jstar+jnopen*ncjsk
         do 360 jdoo = 1,jnopen
          kstar=kstar-ncjsk
*vocl loop,novrec
          do 380 kdoo = 1,nspiel
          idiff(kdoo) = idiff(kdoo)
     &                        + jcon1( iottr(kstar+ipoint(kdoo)) )
380       continue
360      continue
c do loop over closed shells
c start at the end because that is were the most differences are.
         kstar=jstar+jnytl*ncjsk
         do 260 jdoo = 1,jndub
          kstar=kstar-ncjsk
c sum the differences in relation to the test configurations on idiff
          do 270 kdoo = 1,nspiel
          idiff(kdoo) = idiff(kdoo)
     &                        + jcon2( iottr(kstar+ipoint(kdoo)) )
270       continue
260      continue
cbs
cbs now compress
cbs
          kspiel=0
*vocl loop,novrec
          do 385 kdoo = 1,nspiel
           if(idiff(kdoo).le.2) then
             iotnew(5*mcsum+ipoint(kdoo))=jref
             kspiel=kspiel+1
             ispiel(kspiel)=ipoint(kdoo)
           endif
 385      continue
c
c  on ipoint all not yet classified configurations get projected.
          nspie1=0
*vocl loop,novrec
          do 386 kdoo = 1,nspiel
             if(idiff(kdoo).gt.2) then
             nspie1 = nspie1 + 1
             ipoint(nspie1)=ipoint(kdoo)
             endif
 386      continue
c
      mspiel=nspie1
c
c
          do 387 kdoo=1,kspiel
 387      idiff(kdoo)=0
c
c  establishing the creators and annihilators in relation to the 
c  reference configuration
c    ispiel contains the addresses of the kspiel configurations that
c    need to be tested.
c
c   establishing the creators
c
c    do loop over the singly occupied mo's of the configurations
         kstar=jstar-ncjsk
         do 460 jdoo = 1,jnopen
          kstar=kstar+ncjsk
          do 410 kdoo=1,kspiel
            jdiff(kdoo)=jcon1( iottr(kstar+ispiel(kdoo)) )
            if (jdiff(kdoo).ne.0) then
             idiff(kdoo)=idiff(kdoo)+1
c  set iotnew setzen, ni respectively nj already set to 1
             iotnew( 5*mcsum+idiff(kdoo)*ncjsk
     &              +ispiel(kdoo) ) = iottr(kstar+ispiel(kdoo))
            endif
 410      continue
 460     continue
c    do loop over the doubly occupied mo's of the configurations
         do 480 jdoo = 1,jndub
          kstar=kstar+ncjsk
*vocl loop,novrec
          do 470 kdoo = 1,kspiel
            jdiff(kdoo)=jcon2( iottr(kstar+ispiel(kdoo)) )
            if (jdiff(kdoo).ne.0) then
             idiff(kdoo)=idiff(kdoo)+1
c  set iotnew, coding the doubly occupied by the shift of the mo-numbers
             iotnew( 5*mcsum+idiff(kdoo)*ncjsk
     &              +ispiel(kdoo) ) = iottr(kstar+ispiel(kdoo))
     &                               +momax
            endif
c   special case detection for jdiff = 2 !
            if (jdiff(kdoo).eq.2) then
             idiff(kdoo)=idiff(kdoo)+1
             iotnew( 5*mcsum+idiff(kdoo)*ncjsk
     &              +ispiel(kdoo) ) = iottr(kstar+ispiel(kdoo))
            endif
  470     continue
  480    continue
cx       if (jsk.eq.4) then
cx         write(iwr,*)
cx         write(iwr,*) ' jref=,mspiel=',jref,mspiel
cx         write(iwr,*) ' ipoint=',(ipoint(ii),ii=1,mspiel)
cx       endif
c
c   establishing the annihilators
c
      do 505 kdoo=1,kspiel
        idiff(kdoo)=0
  505 continue
c    initialisation of konflg
      do 510 kdoo=1,imo*kspiel
        konflg(kdoo)=0
  510 continue
      kstar=jstar-ncjsk
c    setting the singly occupied mo's
      do 520 jdoo=1,jnopen
        kstar=kstar+ncjsk
        do 530 kdoo=1,kspiel
c        monumb=number of singly occupied mo's
          monumb=iottr(kstar+ispiel(kdoo))
          konflg( (monumb-1)*kspiel+kdoo ) = 1
  530   continue
  520 continue
c    setting the doubly occupied mo's
      do 540 jdoo=1,jndub
        kstar=kstar+ncjsk
        do 550 kdoo=1,kspiel
c        monumb=number of doubly occupied mo's
          monumb=iottr(kstar+ispiel(kdoo))
          konflg( (monumb-1)*kspiel+kdoo ) = 2
  550   continue
  540 continue
c
c   compare the reference against configurations 
c    with singly occupied mo's of the reference
      do 560 jdoo=1,inopen
        do 570 kdoo=1,kspiel
          jdiff(kdoo)=konflg( (itest(jdoo)-1)*kspiel+kdoo )
c  set iotnew
           if (jdiff(kdoo).eq.0) then
             idiff(kdoo)=idiff(kdoo)+1
             iotnew( 5*mcsum+(idiff(kdoo)+2)*ncjsk
     &              +ispiel(kdoo) ) = itest(jdoo)
           endif
  570   continue
  560 continue
c    doubly occupied mo's of the reference
      do 580 jdoo=inopen+1,inytl
*vocl loop,novrec
        do 590 kdoo=1,kspiel
c     special case detection for singly respectively doubly 
c     occupied mo's
          jdiff(kdoo)=konflg( (itest(jdoo)-1)*kspiel+kdoo )
           if (jdiff(kdoo).lt.2) then
             idiff(kdoo)=idiff(kdoo)+1
c  coding the excitation from a doubly occupied mo according to the
c   shift of the mo-number
             iotnew( 5*mcsum+(idiff(kdoo)+2)*ncjsk
     &              +ispiel(kdoo) ) = itest(jdoo) + momax
c   special case detection for double-excitations !
            if (jdiff(kdoo).eq.0) then
             iotnew( 5*mcsum+(idiff(kdoo)+3)*ncjsk
     &              +ispiel(kdoo) ) = itest(jdoo)
            endif
           endif
  590   continue
  580 continue
c
      return
      end
c
c establishing the excitations of configurations from the mains
      subroutine excit(iottr,iotnew,iotm,mtape)
cvp
      implicit real*8 (a-h,p-z)
      parameter (iwod=1000)
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer iotnew, iottr, iotm
      dimension iotnew(iotm), iottr(iotm)
c new fields: inop = # open shells per sk
c             iskop = sk for given # open shells
      integer inop,iskop
      dimension inop(iswhm),iskop(0:nopmax)
      common/mains/inop,iskop
c
c iotnew: new iot-field, contains the number of the generating 
c         reference as well as the creators and annihilators
c         max. size: 5 * (# configurations)
      integer ipoint(ncmax)
      common/iotn/ipoint
c
c
c loop over the sk
      mcsum=0
      do 1000 jsk=1,iswhm
        mc=nconf(jsk)
        if (mc.eq.0) goto 1000
c      nullifying the number of generating references
        do 1010 ii=1,mc
          iotnew(5*mcsum+ii)=0
 1010   continue
c      initialisation of the remaining entries with 1
          do 1020 ii=1,4*mc
            iotnew(5*mcsum+mc+ii)=1
 1020     continue
c   calculation of the number nstrip of cycles over the configurations
c    with smaller sk, just in case that number exceeds ncmax
          nstrip=(mc-1)/ncmax+1
c
c   splitting the loop over all mc configuration into loops with 
c    maximal ncmax cycles
        do 3028 istrip=1,nstrip
c   mspiel = number of configurations to be tested
          mspiel=ncmax
c   jnum1 = number of the first configuration to be tested -1
          jnum1=(istrip-1)*ncmax
          if (istrip.eq.nstrip) then
             mspiel=mc-(nstrip-1)*ncmax
          endif
c   setting the address vector
        do 1030 ii=1,mspiel
          ipoint(ii)=ii + jnum1
 1030   continue
c   loop over all mains
        do 1100 jref=1,nko
          if (abs(iref(0,jref)-jsk).gt.2) goto 1100
             call compar(iotnew,iottr,iotm,
     +       mcsum,jsk,jref,nspiel,mspiel)
c       write(iwr,*) ' jsk=,jref=,mspiel= ',jsk,jref,mspiel
c
 1100   continue
        if (mspiel.ne.0) then
           write(iwr,*) ' super category ',jsk,' ; ',
     &      mspiel,' configurations were not classified! '
           call caserr('problem with super category')
        endif
c
 3028   continue
*       do 123 ii=1,mc
*         write(iwr,*) ' konf=',ii,' generated from: ',
*    &    iotnew(5*mcsum+ii),
*    &   ' erz. mo=',iotnew(5*mcsum+ii+mc),iotnew(5*mcsum+ii+2*mc),
*    &   ' ver. mo=',iotnew(5*mcsum+ii+3*mc),iotnew(5*mcsum+ii+4*mc)
* 123   continue
c
        mcsum=mcsum+mc
 1000 continue
c   writing the second record of mtape
      write(mtape) iotnew
c
      return
      end
      subroutine refcon(q,odebug)
c
      implicit real*8 (a-h,o-z)
c
      logical odebug
      dimension q(*)
c
      real*8 hstori, rtoli, cptoli, cptolmi
      logical bypass, debugd, debugtab, debugci
      common /adlin/ hstori, rtoli, cptoli, cptolmi, cptolcci,
     +               nrooti, ntchi,
     +               kprini, iseleci,
     .               ndecii, icodei, konfi, kepsi,
     .               iointi, norhpi,
     .               izusi, isski,
     .               ifirsti, ilasti, istarti, ndavi,
     .               iggeyi, nformi,
     +               bypass(24),debugd, debugtab, debugci,
     +               nteintin,mdiin,iotmin,nedimin
c
c iotm   : field length for the selected configurations
c          note this may now be ovewritten by data input (mrdcin2)
c mx     : max reference space
c
      integer iotm, mx, momax, maxref
      parameter (momax = 256)
      parameter (maxref = 256)
      parameter (mxroot = 50)
      parameter (mx = 256)
c
      nav = lenwrd()
      mvect = mxroot * maxref
      iotm = (iotmin-1)/nav+1
      iotm = iotm * nav
c
      ivect = 1
      iottr = ivect + mvect
      iotnew = iottr + iotm / nav
      iot = iotnew + iotm / nav
      maindf = iot + iotm / nav
      jconb =  maindf + mx * mx / nav
      momax2 = momax + momax
      need =   jconb  + mx * momax2 / nav
c
c     allocate required memory
c
      ivect = igmem_alloc(need)
      iottr = ivect + mvect
      iotnew = iottr + iotm / nav
      iot    = iotnew + iotm / nav
      maindf = iot + iotm / nav
      jconb = maindf + mx * mx / nav
c
      call refcon2(q(ivect),mvect,q(iottr), q(iotnew), 
     +             q(iot), iotm, 
     +             q(maindf), mx, q(jconb), momax2,
     +             odebug)
c
c      now free memory
c
       call gmem_free(ivect)
c
      return
      end
      subroutine refcon2(vect1,mvect,iottr,iotnew,iot,iotm,
     +                   maindf,mrx,jconb,momax2,odebug)
c
cvp  program to compare all configurations with the reference 
cvp  configurations. Establishes from which reference configuration
cvp  a configuration was generated and which are the corresponding 
cvp  creation and annihilation operators. The special case detection
cvp  with respect to singly and doubly occupied mo's gets coded in the
cvp  mo numbers.
c
      implicit real*8 (a-h,p-z)
      real*8 vect1
      integer mvect
      dimension vect1(mvect)
      integer iottr, iotnew, iot, iotm
      dimension iottr(iotm), iotnew(iotm), iot(iotm)
      integer maindf, mrx, jconb, momax2
      dimension maindf(mrx,mrx), jconb(mrx,momax2)
      logical odebug
c
      parameter (iwod=1000)
cvp parameter und globale common-bloecke sind im include-file
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
cvp
      real*8 ew
      common/junk/ew(mxroot)
c
      integer jkan
      dimension jkan(iwod)
c
      integer mconf
      integer nplu
      dimension mconf(iswhm),nplu(iswhm)
c
c itest2: mains 
c
      integer itest2
      dimension itest2(maxshl)
c
c new fields: inop = # open shells per sk
c             iskop = sk for given # open shells
      integer inop,iskop
      dimension inop(iswhm),iskop(0:nopmax)
      common/mains/inop,iskop
c
      integer  ntape,mtap32,mdisk
      integer  ideli,ltype,linf,ntab,kfile,kclab
      integer  ntype, mstvt, nf01, nf62, nhead
      integer  nf99, nf78, nf11
      common /ftap/ ntape, mtap32, mdisk,
     .               ideli, ltype, linf,  ntab,  kfile,
     .               kclab, ntype, mstvt, nf01, nf62,
     +               nhead, nf99,  nf78, nf11
      character*10 charwall
c
c iotnew: new iot-field, contains the number of the generating reference
c         and the creators and annihilators.
c         max. size: 5 * (# konfigurationen)
c
c     namelist /refin/ nele,nko
cvp
c
c ft99 should contain the reference configurations
  11   format(/1x,104('*'))
  12   format(/5x,'***  start of refcon at ',f10.2,' seconds',
     +         a10,' wall'/)
       write (iwr,11)
       write (iwr,12) cpulft(1) ,charwall()
c
      ktape=nf99
c ft33 contains the configurations
      ltape=linf
c ft78 contains the output of the program
      mtape=nf78
      call rewftn(ktape)
      call rewftn(ltape)
      call rewftn(mtape)
c
c reading the configurations
      read (ltape) jsec,nrootx,nytl,nplu,ndub,vect1,ew,mconf
c     write(iwr,*) ' mconf= ',mconf
c
      ny=0
      do 7 i=1,iswhm
      niot(i)=ny
      nod(i)=nytl(i)-ndub(i)
      if (mconf(i).eq.0) go to 7
      ig=1
      nx=nytl(i)
      jg=0
      read (ltape)
      read (ltape)
      read (ltape) jkan
   8  if (jkan(ig).eq.0) go to 9
      jg=jg+1
      do 10 j=1,nx
      ny=ny+1
      if(ny.le.iotm) iottr(ny)=jkan(ig)
      if (ig.lt.iwod) go to 10
      read (ltape) jkan
      ig=0
   10 ig=ig+1
      go to 8
    9 nconf(i)=jg
    7 continue
      if (ny.gt.iotm) then
       write(iwr,51234) ny, iotm
51234  format(1x,
     +      'insufficient memory for configurations in refcon'/
     +      1x,'required (ny)    = ', i12/
     +      1x,'available (iotm) = ', i12/)
         call caserr(
     + 'insufficient memory for configurations - increase iotm')
      endif
      call rewftn(ltape)
c
      write(iwr,61200) nconf
61200 format(/1x,
     +   'number of selected configurations in each super category'/
     +    5x,12i7)
c abfrage, ob iotnew zu klein
      mcsum=0
      do 30 i=1,iswhm
        mcsum=mcsum+nconf(i)
   30 continue
      if (5*mcsum.gt.iotm) then
       write(iwr,51235) 5*mcsum, iotm
51235  format(1x,
     +     'insufficient memory for configurations in refcon'/
     +      1x,'required (5*msum)= ', i12/
     +      1x,'available (iotm) = ', i12/)
         call caserr(
     + 'insufficient memory for configurations - increase iotm')
      endif
cvp
cvp building fields that should facilitate comparing configurations:
cvp here in particular the transposition of the iot-field
      call chfeldr(iottr,iot,iotm)
c
      do 5 i=1,iswhm
      inop(i)=nytl(i)-ndub(i)
      iskop(inop(i))=i
    5 continue
c
      nref=0
c     read(ktape) nr
c     write(iwr,*) ' number of ??????????: ',nr
c     read(ktape) (itest(ll),ll=1,maxshl)
c     nref=nref+1
c     write(iwr,*) ' reference nr ',nref,':',(itest(ll),ll=1,maxshl)
c input vom ft05
cread read(5,refin,end=20)
cre20 continue
cread do 15 i=1,nko
cread read(5,*) nr
c lesen von ft99
      read(ktape,*) 
      read(ktape,*)
      read(ktape,*) nko,nele
      write(iwr,*)
      write(iwr,*) ' number of reference configurations: ',nko
      write(iwr,*) ' number of electrons:                ',nele
      write(iwr,*)
      do 15 i=1,nko
      read(ktape,*) nr
      norb=(nele+nr)/2
      read(ktape,*) (itest2(ii),ii=1,norb)
      do 23 j=1,norb
	itest(j)=itest2(j)
   23 continue
cread read(5,*) (itest(ll),ll=1,norb)
c   establishing the sk of the reference and storing it on iref
      isk=iskop(nr)
      iref(0,i)=isk
      do 17 j=1,norb
        iref(j,i)=itest(j)
   17 continue
      if (odebug) then
       write(iwr,*) ' number of singly occupied shells ',nr
       write(iwr,*) ' reference from the sk: ',iref(0,i)
       write(iwr,*) ' reference ',(itest(ll),ll=1,norb)
       write(iwr,*)
      endif
   15 continue
c
c comparing the mains among eachother
      call mcomp(maindf,mrx,jconb,momax2,mtape)
c
c establishing the excitation of the configurations from the mains
      call excit(iottr,iotnew,iotm,mtape)
c
      call rewftn(ktape)
      call rewftn(ltape)
      call rewftn(mtape)
c
 230  format(/104('-'))
      write(iwr,22)cpulft(1) ,charwall()
  22  format(/5x,'***  end of refcon at ',f10.2,' seconds',a10,' wall')
      write(iwr,230)
      return
      end
c
c comparing the main among eachother
      subroutine mcomp(maindf,mrx,jconb,momax2,mtape)
cvp
      implicit real*8 (a-h,p-z)
      parameter (iwod=1000)
cvp parameters and globale common-bloecks are in include-file
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer maindf, mrx, jconb, momax2
      dimension maindf(mrx,mrx),jconb(mrx,momax2)
c neue felder: inop = # offene schalen pro sk
c              iskop = sk fuer geg. # offener schalen
c              iref: iref(0,*) = sk der referenz, dahinter die
c                    einfach und doppelt bes. mo's
      integer inop,iskop
      dimension inop(iswhm),iskop(0:nopmax)
      common/mains/inop,iskop
c
c iotnew: new iot-field, contains the number of the generating reference
c         and the creators and annihilators.
c         max. size: 5 * (# configurationen)
      integer ipoint(ncmax)
      common/iotn/ipoint
c
c maindf: idiff of the mains among eachother
c jconb : long jcon-field for all mains
c rwork is 'scratch-common'
c
c zeroing of maindf
      do 90 i=1,maxref
       do 95 j=1,maxref
        maindf(j,i)=0
   95  continue
   90 continue
c loop over all configuration paires
      do 100 lref=1,nko
c      initialisation of the test configuration
        isk=iref(0,lref)
c
        inytl  = nytl(isk)
        indub  = ndub(isk)
        inopen = inytl - indub
c nullifying the jcon field
c   imo = number of electrons, is read from ft31, therefore here
        imo=momax
        do 215 k=1,imo
        jcon2(k)=2
        jcon1(k)=1
215     continue
c   storing the test configuration on itest
         do 205 k=1,inytl
  205    itest(k)=iref(k,lref)
c        write(iwr,*) ' itest= ',(itest(ii),ii=1,inytl)
         do 220 k=1,inopen
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 1
  220    continue
         do 221 k=inopen+1,inytl
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 0
  221    continue
c
c   building jconb: contains the special case detection, the excitation
c    needed with respect to singly or doubly occupied mo's
         do 225 k=1,imo
            jconb(lref,k)=jcon1(k)
            jconb(lref,k+momax)=jcon2(k)-jcon1(k)
  225    continue
c
         do 200 jref=lref+1,nko
           jsk=iref(0,jref)
           jnytl  = nytl(jsk)
           jndub  = ndub(jsk)
           jnopen = jnytl - jndub
c
c         compare for each pair seperately
           ldiff=0
           do 230 jdoo=1,jnopen
             ldiff=ldiff+jcon1( iref(jdoo,jref) )
  230      continue
           do 240 jdoo=jnopen+1,jnytl
             ldiff=ldiff+jcon2( iref(jdoo,jref) )
  240      continue
c        maindf is used as a square matrix
            maindf(lref,jref)=ldiff
            maindf(jref,lref)=ldiff
  200    continue
  100 continue
c
*      write(iwr,*)' excitation levels of the mains among each other '
*      do 310 i=1,nko
*       write(iwr,*)
*       write(iwr,*) (maindf(i,j),j=1,nko)
*  310 continue
c  writing the first record of mtape
      write(mtape) nko,maindf,iref,jconb
c
      return
      end
      subroutine ver_newmrd4(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/newmrd4.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
