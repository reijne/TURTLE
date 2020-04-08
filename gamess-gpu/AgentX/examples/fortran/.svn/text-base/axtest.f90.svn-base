
! A simple Fortran 90 example
! Read in molecular structures, atomic basis sets, basis set assignments and molecular orbitals
! Print data to the screen

PROGRAM simplefcode
  
  USE fagentxi
  
  IMPLICIT NONE
  
  INTEGER, parameter :: skind=KIND(1.0), dkind=KIND(1D0)
  INTEGER::retval=0
  INTEGER::i=0,j=0,k=0,l=0,noBasisSetAssignment=0,noMolecule=0
  INTEGER::noprop=0
  INTEGER::noAtom=0,noAtomicBasisSet=0,noBasisGroup=0
  INTEGER::noMolecularOrbital=0,noMolecularOrbitalSet=0,noMOCoeffs=0
  CHARACTER(LEN=30)::id,elementtype,prop
  REAL(kind=dkind)::xcoordinate,ycoordinate,zcoordinate
  REAL(kind=dkind)::coeff,expon,dprop
  INTEGER::minl=0,maxl=0,iprop
  
  ! initialise
  
  retval=axParserStart()

  ! retval is used to collect error codes.
  ! an error has occured if this is < 0
  
  ! set the data format to XML
  
  retval=axFormat(1)
  
  ! specify the number of evaluations to cache */
  
  retval=axCache(3000)
  
  ! open the data documents
  
  retval=axDataGetUri("../xml/molpro_corrected.xml")
  
  ! open supporting documents
  
  ! ontologies (logical model)
  
  retval=axGetUri("../../ontology/ontology.owl")
  
  ! mappings (logical-physical mappings)

  ! load Molpro AgentX extensions
  
  retval=axGetUri("../../map/molpro.rdf")

  ! load standard mappings
  
  retval=axGetUri("../../map/map.rdf")
  
  ! All concepts and properties are identified by unique URIs. The following 
  ! function call sets the prefix for these URIs.
  
  retval=axBaseUri("http://www.grids.ac.uk/eccp/owl-ontologies#")
  
  ! Concepts used in function calls are prefixed by the base URI.
  ! For example, axSelect("Molecule") is transformed to axSelect
  ! ("http://www.grids.ac.uk/eccp/owl-ontologies#Molecule").
  
  noMolecule=axSelect("Molecule")
  
  ! This call locates all datasets that relate to the 
  ! concept 'http://www.grids.ac.uk/eccp/owl-ontologies#Molecule'.
  ! The number of data sets found is returned
  
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
        
        noprop = axSelect("xCoordinate")
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
        
        noBasisGroup=axSelect("BasisGroup")
        
        do l=1,noBasisGroup
           
           noprop=axSelect("identifier")
           if(noprop.gt.0)then
              call axValue(id)
              print*,"basis group id = ",id
              retval=axDeselect()
           endif
           
           retval=axSelectNext()

        ! this axSelectNext call selects the next data set located by the 
        ! last axSelect call.           
 
       enddo
        
        if (noBasisGroup .gt. 0) then
           retval=axDeselect()
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
  
  noAtomicBasisSet=axSelect("AtomicBasisSet")
  
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
     
     do k=1,noBasisGroup
        
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
           retval=axSelectNext()
        end do
        if (noprop .gt. 0) then
           retval=axDeselect()
        endif
        
        retval=axSelectNext()
        
     end do
     
     if (noBasisGroup .gt. 0) then
        retval=axDeselect()
     endif
     retval=axSelectNext()
     
  end do
  
  if (noAtomicBasisSet .gt. 0) then      
     retval=axDeselect()
  endif
  
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
           print*,coeff
           retval=axSelectNext()
        end do
        if (noprop .gt. 0) then
           retval=axDeselect()
        endif
        
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
  
END PROGRAM simplefcode

