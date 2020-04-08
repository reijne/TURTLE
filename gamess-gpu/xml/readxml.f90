
! A simple Fortran 90 example
! Read in molecular structures, atomic basis sets, basis set assignments and molecular orbitals
! Print data to the screen

SUBROUTINE readxml_init(file)
  
  USE fagentxi
  
  IMPLICIT NONE

  character file*(*)  

  INTEGER::retval=0

  ! initialise
  retval=axParserStart()

  ! retval is used to collect error codes.
  ! an error has occured if this is < 0
  
  ! set the data format to XML
  
  retval=axFormat(1)
  
  ! specify the number of evaluations to cache */
  
  retval=axCache(3000)
  
  ! open the data documents
  
  write(6,*)'xml init on file',file
  retval=axDataGetUri(file)
  
  ! open supporting documents
  
  ! ontologies (logical model)
  
  retval=axGetUri("ontology.owl")
  
  ! mappings (logical-physical mappings)

  ! load Molpro AgentX extensions
  
  retval=axGetUri("mapext.rdf")

  ! load standard mappings
  
  retval=axGetUri("map.rdf")
  
  ! All concepts and properties are identified by unique URIs. The following 
  ! function call sets the prefix for these URIs.
  
  retval=axBaseUri("http://www.grids.ac.uk/eccp/owl-ontologies#")
  
  ! Concepts used in function calls are prefixed by the base URI.
  ! For example, axSelect("Molecule") is transformed to axSelect
  ! ("http://www.grids.ac.uk/eccp/owl-ontologies#Molecule").

END SUBROUTINE readxml_init

subroutine readxml_getcoords(natoms,coord,symbols)

  ! first pass through loading coordinates and figuring out 
  ! how many unique basis sets there are for subsequent loading

  USE fagentxi
  USE xml_modules
  USE known_basis

  INTEGER natoms
  REAL(KIND=8):: coord(3,*)
  CHARACTER*8 symbols(*)
  INTEGER, parameter :: skind=KIND(1.0), dkind=KIND(1D0)
  INTEGER::i=0,j=0,k=0,l=0,noBasisSetAssignment=0,noMolecule=0
  INTEGER::noprop=0
  INTEGER::noAtom=0,noAtomicBasisSet=0,noBasisGroup=0
  INTEGER::noMolecularOrbital=0,noMolecularOrbitalSet=0,noMOCoeffs=0
  CHARACTER(LEN=30)::id,elementtype,prop
  REAL(kind=dkind)::xcoordinate,ycoordinate,zcoordinate
  REAL(kind=dkind)::coeff,expon,dprop
  INTEGER::minl=0,maxl=0,iprop

  logical debug
  common/dbg/debug

  real*8 pi, twopi, root3, rt3by2, bohr, degree
  common/cnstnt/pi,twopi,root3,rt3by2,bohr,degree


  logical known_bas

  integer ibas


  integer npos
  character*2 elmnt(maxat)

  integer ix1
  !  integer icharg(maxat)

  integer nprim
  integer nda

  !
  character*2 sbf

  ! mapping atomic numbers and symbols
  !integer ich
  !external ich

  integer  iq, if, lq, nq, ml, nctr
  real*8 anorm

  logical basis_match

  noMolecule=axSelect("Molecule")

  write(6,*)'noMolecule',noMolecule
  
  ! This call locates all datasets that relate to the 
  ! concept 'http://www.grids.ac.uk/eccp/owl-ontologies#Molecule'.
  ! The number of data sets found is returned

  num_atomic_basis_sets = 0

  do j=1,noMolecule
     
     noprop=axSelect("identifier")

     ! This call locates all data elements that relate to the
     ! property 'http://www.grids.ac.uk/eccp/owl-ontologies#identifier'.
     ! The number of data elements found is returned

     if(noprop.gt.0)then
        call axValue(id)

        ! This call gets the value of the previously selected data element

        print*,"serialisation of Molecule :",id
        retval=axDeselect()

        ! This call drops the data element selection

     endif

     noAtom=axSelect("Atom")

     natoms = noAtom
     write(6,*)'noAtom assigned',noatom, natoms
     
     do k=1,noAtom
        
        noprop=axSelect("identifier")
        if(noprop.gt.0)then
           call axValue(id)
           print*,"serialisation of Atom :",id
           retval=axDeselect()
        endif
        
        noprop=axSelect("elementType")
        if(noprop.gt.0)then
           call axValue(elementtype)
           print*,"element = ",elementtype
           retval=axDeselect()
        endif

        symbols(k)=elementtype
!!       icharg(k)=ich(elmnt(k))

        noprop=axSelect("xCoordinate")
        if(noprop.gt.0)then
           retval=axValueConvert(xcoordinate)

           ! axValueConvert will set the argument to the value of the 
           ! selected data element.
           ! It is a generic function - types are determined
           print*,"xCoordinate = ",xcoordinate
           retval=axDeselect()
        endif
        
        noprop=axSelect("yCoordinate")
        if(noprop.gt.0)then
           retval=axValueConvert(ycoordinate)
           print*,"yCoordinate = ",ycoordinate
           retval=axDeselect()
        endif
        
        noprop=axSelect("zCoordinate")
        if(noprop.gt.0)then
           retval=axValueConvert(zcoordinate)
           print*,"zCoordinate = ",zcoordinate
           retval=axDeselect()
        endif
        
        coord(1,k) = xcoordinate / 0.529177249
        coord(2,k) = ycoordinate / 0.529177249
        coord(3,k) = zcoordinate / 0.529177249

        ! locate first identifier
 
        noBasisGroup=axSelect("BasisGroup")
        noprop=axSelect("identifier")
        call axValue(id)
        write(6,*)'first id of basis group',id
        retval=axDeselect()
        known_bas = .false.
        do ibas = 1, num_atomic_basis_sets
           if(id .eq. basis_ids(ibas,1) ) then
              known_bas = .true.
              atom_basis(k) = ibas
           endif
        enddo
        ! deselect the group 
        retval=axDeselect()
        ! WHY NOT

        if (.not. known_bas) then
           noBasisGroup=axSelect("BasisGroup")
           num_atomic_basis_sets =  num_atomic_basis_sets + 1  
           len_atomic_basis_set(num_atomic_basis_sets) =  noBasisGroup
           atom_basis(k) = num_atomic_basis_sets
           do l=1,noBasisGroup
              noprop=axSelect("identifier")
              if(noprop.gt.0)then
                 call axValue(id)
                 print*,"basis group id = ",id
                 retval=axDeselect()
                 basis_ids(num_atomic_basis_sets,l) = id
              endif
              retval=axSelectNext()
              ! this axSelectNext call selects the next data set located by the 
              ! last axSelect call.           
           enddo
           if (noBasisGroup .gt. 0) then
              retval=axDeselect()
           endif
        endif

        retval=axSelectNext()
        
     end do
     
     if (noAtom .gt. 0) then
        retval=axDeselect()
     endif
     retval=axSelectNext()
     
  end do
  
  if (noMolecule .gt. 0) then
     retval=axDeselect()
  endif

  write(6,*)'Basis sets detected',num_atomic_basis_sets
  do i = 1,num_atomic_basis_sets
     write(6,*)(basis_ids(i,j)(1:4),j=1,len_atomic_basis_set(i))
  enddo

  END SUBROUTINE readxml_getcoords

  subroutine readxml_loadbasis(csinpi,cpinpi,cdinpi)


  USE fagentxi
  USE xml_modules
  USE known_basis

  IMPLICIT NONE

  REAL*8 csinpi(*)
  REAL*8 cpinpi(*)
  REAL*8 cdinpi(*)

  INTEGER, parameter :: skind=KIND(1.0), dkind=KIND(1D0)
  INTEGER::i=0,j=0,k=0,l=0,noBasisSetAssignment=0,noMolecule=0
  INTEGER::noprop=0
  INTEGER::noAtom=0,noAtomicBasisSet=0,noBasisGroup=0
  INTEGER::noMolecularOrbital=0,noMolecularOrbitalSet=0,noMOCoeffs=0
  CHARACTER(LEN=30)::id,elementtype,prop
  REAL(kind=dkind)::xcoordinate,ycoordinate,zcoordinate
  REAL(kind=dkind)::coeff,expon,dprop
  INTEGER::minl=0,maxl=0,iprop

  INTEGER retval, ibas, ix, nprim, natoms
  logical known_bas

  logical readxml_bas_atom_check


  logical debug
  common/dbg/debug

   INTEGER shell, atom

   REAL*8 trmat, ptr, dtr, ftr, gtr
   INTEGER intyp, ns, ks, intypi, nsi,ksi
   REAL*8 exi
   REAL*8 csi, cpi, cdi, cfi, cgi 
   INTEGER kstari, katomi, ktypi, kngi, kloci, kmini, kmaxi
   INTEGER nshi 
   common/junk/trmat(432), &
        ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),  &
        intyp(mxshel),ns(maxat),ks(maxat), &
        intypi(mxshel),nsi(maxat),ksi(maxat), &
        exi(mxprim), &
        csi(mxprim),cpi(mxprim),cdi(mxprim),cfi(mxprim),cgi(mxprim), &
        kstari(mxshel),katomi(mxshel),ktypi(mxshel),kngi(mxshel), &
        kloci(mxshel),kmini(mxshel),kmaxi(mxshel),nshi

   INTEGER itmp

  noAtomicBasisSet=axSelect("AtomicBasisSet")
  
  nshi = 0
  ! should just be one instance here
  do j=1,noAtomicBasisSet
     
     noprop=axSelect("identifier")
     if(noprop.gt.0)then
        call axValue(id)
        print*,"serialisation of AtomicBasisSet :",id
        retval=axDeselect()
     endif
     
     noprop=axSelect("length")
     if(noprop.gt.0)then
        retval=axValueConvert(iprop)
        print*,"length :",iprop
        retval=axDeselect()
     endif
     

     noprop=axSelect("angular")
     if(noprop.gt.0)then
        call axValue(prop)
        print*,"angular :",prop
        retval=axDeselect()
     endif
     
     noprop=axSelect("type")
     if(noprop.gt.0)then
        call axValue(prop)
        print*,"type :",prop
        retval=axDeselect()
     endif
     
     noprop=axSelect("groups")
     if(noprop.gt.0)then
        retval=axValueConvert(iprop)
        print*,"groups :",iprop
        retval=axDeselect()
     endif
     
     noprop=axSelect("cartesianLength")
     if(noprop.gt.0)then
        retval=axValueConvert(iprop)
        print*,"cartesianLength :",iprop
        retval=axDeselect()
     endif
     
     noprop=axSelect("noPrimitives")
     if(noprop.gt.0)then
        retval=axValueConvert(iprop)
        print*,"noPrimitives :",iprop
        retval=axDeselect()
     endif
     
     noBasisGroup=axSelect("BasisGroup")

     if (j .eq. 2)then
        stop 'too many basis sets'
     endif

     ! counter for primitives
     ix = 0
     nprim = 0


     do k=1,noBasisGroup

        nshi = nshi + 1
        
        noprop=axSelect("identifier")
        if(noprop.gt.0)then
           call axValue(id)
           print*,"serialisation of BasisGroup: ",id
           retval=axDeselect()
        endif
        
        noprop=axSelect("minL")
        if(noprop.gt.0)then
           retval=axValueConvert(minl)
           print*,"minL = ",minL
           retval=axDeselect()
        endif
        
        noprop=axSelect("maxL")
        if(noprop.gt.0)then
           retval=axValueConvert(maxl)
           print*,"maxL = ",maxL
           retval=axDeselect()
        endif
        
        noprop=axSelect("noPrimitives")
        if(noprop.gt.0)then
           retval=axValueConvert(iprop)
           print*,"noPrimitives :",iprop
           retval=axDeselect()
        endif
        
        noprop=axSelect("angular")
        if(noprop.gt.0)then
           call axValue(prop)
           print*,"angular :",prop
           retval=axDeselect()
        endif
        
        noprop=axSelect("noContractions")
        if(noprop.gt.0)then
           retval=axValueConvert(iprop)
           print*,"noContractions :",iprop
           retval=axDeselect()
        endif
        
        noprop=axSelect("coefficient")
        print*,"coefficients: "
        do i = 1,noprop
           retval=axValueConvert(coeff)
           print*,coeff
           ! coef1(ix + i)=coeff

           if ( (minL .eq. 0) .and. (maxL .eq. 0)) then
              csi(nprim + i) = coeff
              csinpi(nprim + i) = coeff
           else if ( (minL .eq. 1) .and. (maxL .eq. 1)) then
              cpi(nprim + i) = coeff
              cpinpi(nprim + i) = coeff
           else if ( (minL .eq. 2) .and. (maxL .eq. 2)) then
              cdi(nprim + i) = coeff
              cdinpi(nprim + i) = coeff
           else
              stop 'no fs yet'
           endif

           retval=axSelectNext()
        end do
        if (noprop .gt. 0) then
           retval=axDeselect()
        endif
        
        noprop=axSelect("exponent")
        print*,"exponents: "
        do i = 1,noprop
           retval=axValueConvert(expon)
           print*,expon
           ! expo1(ix+i)=expon

           exi(nprim + i) = expon

           retval=axSelectNext()
        end do
        if (noprop .gt. 0) then
           retval=axDeselect()
        endif

        if ( (minL .eq. 0) .and. (maxL .eq. 0)) then
           intypi(nshi) = 16
           kmini(nshi) = 1
           kmaxi(nshi) = 1
           ktypi(nshi) = 1
        else if ( (minL .eq. 1) .and. (maxL .eq. 1)) then
           intypi(nshi) = 17
           kmini(nshi) = 2
           kmaxi(nshi) = 4
           ktypi(nshi) = 2
        else if ( (minL .eq. 2) .and. (maxL .eq. 2)) then
           intypi(nshi) = 18
           kmini(nshi) = 5
           kmaxi(nshi) = 10
           ktypi(nshi) = 3
        else
           stop 'no fs yet'
        endif

        itmp = 0
        do atom = 1,nat
           if ( (itmp .eq. 0) .and.  readxml_bas_atom_check(nshi,atom) ) then
              itmp = atom
           endif
        enddo

        katomi(nshi) = itmp
        kstari(nshi) = nprim + 1
        kngi(nshi) = noprop

        nprim = nprim + noprop

        retval=axSelectNext()
        
     end do
     
     write(6,*)'total number or primitives processed',nprim

     nshi = noBasisGroup

     if (noBasisGroup .gt. 0) then
        retval=axDeselect()
     endif
     retval=axSelectNext()
     
  end do
  
  if (noAtomicBasisSet .gt. 0) then      
     retval=axDeselect()
  endif

1230 format(2x,i2,4x,a4,2x,2h( ,f7.4,3h , , f7.4,3h , ,f7.4,2h ),2i4)

  if(nprim.gt.mxprim)stop 'maximum number of primitive basis fns exceeded'
  !if(ntbf.gt.maxorb) stop 'maximum number of contracted basis fns exceeded'

   if(debug)then
      write(6,*)'primitive functions'
      do i=1,nprim
!         write(nout,901)coef(i),expo(i),icnt(i),ityp(i),ictr(i)
      enddo
   endif
901 format(1x,2f12.5,3i5)


END SUBROUTINE readxml_loadbasis

SUBROUTINE readxml_getorbitals(mocoeff,moener,numao)

  USE fagentxi
  USE xml_modules
  USE known_basis

  INTEGER::noprop=0
  INTEGER::noAtom=0,noAtomicBasisSet=0,noBasisGroup=0
  INTEGER::noMolecularOrbital=0,noMolecularOrbitalSet=0,noMOCoeffs=0
  CHARACTER(LEN=30)::id,elementtype,prop

  INTEGER numao
  REAL*8 mocoeff(numao,numao)
  REAL*8 moener(numao)
  

  logical debug
  common/dbg/debug


  noMolecularOrbitalSet=axSelect("MolecularOrbitalSet");
  
  do i=1,noMolecularOrbitalSet
     
     print*,"serialisation of MolecularOrbitalSet :",i
     
     noprop=axSelect("basisSetReference")
     if(noprop.gt.0)then
        call axValue(prop)
        print*,"basisSetReference :",prop
        retval=axDeselect()
     endif
     
     noprop=axSelect("spin")
     if(noprop.gt.0)then
        call axValue(prop)
        print*,"spin :",prop
        retval=axDeselect()
     endif
     
     noprop=axSelect("method")
     if(noprop.gt.0)then
        call axValue(prop)
        print*,"method:",prop
        retval=axDeselect()
     endif
     
     noMolecularOrbital=axSelect("MolecularOrbital");
     
     do j=1,noMolecularOrbital
        
        print*,"serialisation of MolecularOrbital :",j
        
        noprop=axSelect("basisSetReference")
        if(noprop.gt.0)then
           call axValue(prop)
           print*,"basisSetReference :",prop
           retval=axDeselect()
        endif
        
        noprop=axSelect("spin")
        if(noprop.gt.0)then
           call axValue(prop)
           print*,"spin :",prop
           retval=axDeselect()
        endif
        
        noprop=axSelect("method")
        if(noprop.gt.0)then
           call axValue(prop)
           print*,"method:",prop
           retval=axDeselect()
        endif
        
        noprop=axSelect("occupation")
        if(noprop.gt.0)then
           retval=axValueConvert(iprop)
           print*,"occupation:",iprop
           retval=axDeselect()
        endif
        
        noprop=axSelect("energy")
        if(noprop.gt.0)then
           retval=axValueConvert(dprop)
           print*,"energy:",dprop
           moener(j) = dprop
           retval=axDeselect()
        endif
        
        noprop=axSelect("symmetryId")
        if(noprop.gt.0)then
           retval=axValueConvert(iprop)
           print*,"symmetryId:",iprop
           retval=axDeselect()
        endif
        
        noprop=axSelect("mOCoefficient")
        print*,"coefficients: "
        do k = 1,noprop
           retval=axValueConvert(coeff)
           print*,j,k,coeff
           mocoeff(k,j) = coeff
           retval=axSelectNext()
        end do
        if (noprop .gt. 0) then
           retval=axDeselect()
        endif
        
        retval=axSelectNext()

     end do
     
     if (noMolecularOrbital .gt. 0) then
        retval=axDeselect()
     endif
     
  end do
  
  if (noMolecularOrbitalSet .gt. 0) then
     retval=axDeselect()
  endif
  
  retval=axParserFinish()
  
  ! the final call to ParserFinish cleans up
  

!   write(6,*)'final npos',npos
!   write(6,*)'final elmnt',elmnt

END SUBROUTINE readxml_getorbitals


LOGICAL FUNCTION  readxml_bas_atom_check(shell,atom)

  USE xml_modules
  USE known_basis

  INTEGER shell
  INTEGER atom

  INTEGER ibas, jbas, kbas

  INTEGER itmp

  ! first identify which group the shell belongs to
  do ibas = 1, num_atomic_basis_sets
     do jbas = 1, len_atomic_basis_set(ibas)

        read(basis_ids(ibas,jbas),'(I1)'),itmp
        
        if(shell .eq. itmp ) then
           kbas = ibas
        endif
     enddo
  enddo

  if (atom_basis(atom) .eq. kbas) then
     readxml_bas_atom_check =  .true.
  else
     readxml_bas_atom_check =  .false.
  endif

END FUNCTION  readxml_bas_atom_check
