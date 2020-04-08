  subroutine xmlStopWriter()
    
    use writexmli

    logical opg_root
    if ( .not. opg_root() ) return
!    if( xmlStatus .eq. 0 ) call xmlStartWriter()
    if( xmlStatus .eq. 0 ) return

    if ( xmlCML ) call xmlEndCML()

    call xmlDecreaseIndent()
    write( xmlFileId, * )""
    call xmlIndenter()
    write( xmlFileId, * )"</eccp>"
    xmlStatus = 2
    close(xmlFileId)
    
  end subroutine xmlStopWriter
  
  subroutine xmlWriteCoords( ztitle, atom_tag, charges, coords, natoms, index )
    
    use writexmli

    integer, parameter :: lentitle = 10
    integer :: natoms
    character(len=*), dimension( lentitle ) :: ztitle
    character(len=8*lentitle) :: moltitle
    character(len=*), dimension( 1:natoms ) :: atom_tag
    character(len=2) :: element, charge2element
    real( kind=xmlDP ), dimension( 1:3, 1:natoms ) :: coords
    real( kind=xmlDP ), dimension( 1:natoms ) :: charges
    real( kind=xmlDP ) :: x, y, z
!   see m4/common/convert
    real( kind=xmlDP ), parameter :: bohr_to_angs = 0.529177249
    integer :: index
    logical :: opg_root ! For parallel code

    character(len=lentitle) :: molname, strindex
    integer start,end

    ! Only root should do this in parallel
    if (.not. opg_root() ) return
    
    if( xmlStatus.eq.0 ) call xmlStartWriter()
!    call xmlStartCML()

    call xmlNumberToString( index, strindex )
    molname = "mol" // strindex

    write(moltitle,'(10a8)')ztitle
    call xmlStartMolecule( molname , moltitle )

    do i = 1, natoms
       element = charge2element( charges(i) )
       x = coords( 1, i ) * bohr_to_angs
       y = coords( 2, i ) * bohr_to_angs
       z = coords( 3, i ) * bohr_to_angs
       call xmlAddAtom( atom_tag( i ), element, x, y, z )
    end do
    call xmlEndMolecule
!    call xmlEndCML()
    
  end subroutine xmlWriteCoords

  subroutine xmlWriteInitialMetadata
  ! Write out any helpful metadata
    use writexmli
    
    Implicit None
    character(len=8) :: date,time,year,accno,jobname
    character(len=256) :: executable
    character(len=300) :: m4keys
    character(len=44) :: hostname
    integer :: jobsecs
    !Hack required for itanium
    integer(4),parameter ::zero = 0
!   For calling blkrun (to get the title)
    character(len=80):: title
    character(len=25):: runt
    character(len=51):: scftype

    logical :: opg_root ! for parallel case
    if (.not. opg_root() ) return

    if( xmlStatus.eq.0 ) call xmlStartWriter()

    call xmlAddMetadataList()

    !Static stuff
    call xmlAddMetadata( "Program", "GAMESS-UK" )
    call xmlAddMetadata( "Version", "7.0 Development Branch" )
    call getm4keys(m4keys)
    call xmlAddMetadata( "m4keys", trim(m4keys) )

    call tidajt(date,time,year,accno,jobname,jobsecs)
    call xmlAddMetadata( "date", trim(date)//' '//trim(year)//' '//trim(time) )
    call xmlAddMetadata( "uid",trim(accno) )
    call xmlAddMetadata( "user",trim(jobname) )

    !call get_command_argument(0,executable)
    call getarg(zero,executable)
    call xmlAddMetadata( "executable",trim(executable) )

    call gethost(hostname)
    call xmlAddMetadata( "hostname",trim(hostname) )

    call blkrun(.false.,runt,scftype,title)
    call xmlAddMetadata( "title",trim(title) )

    call xmlEndMetadataList()

  end subroutine xmlWriteInitialMetadata    


  subroutine xmlAddMetadata( name, content )
    !
    ! Add metadata to a metadatalist
    !
    use writexmli
    
    character(len=*), intent(in) :: name, content
    logical :: opg_root ! for parallel case
    if (.not. opg_root() ) return

    if( xmlStatus.eq.0 ) call xmlStartWriter()

    call xmlIndenter()
    write( xmlFileId, * ) "<metadata name='",trim(name),"' content='",trim(content),"'></metadata>"

  end subroutine xmlAddMetadata
  
  subroutine xmlAddExtraMetadata( name, value )
    !
    ! Subroutine to write additional metadata that may not fit elsewhere
    ! This starts/stops the cml block if necessary
    !
    use writexmli
    
    character(len=*), intent(in) :: name, value
    logical :: opg_root ! for parallel case
    logical :: openedCML = .False. ! Flag to see if we need to close the CML block

    if (.not. opg_root() ) return

    if( xmlStatus.eq.0 ) call xmlStartWriter()

    ! Start CML if necessary
!    if (.not. xmlCML ) then
!       call xmlStartCML()
!       openedCML=.True.
!    endif

    call xmlIndenter()
    write( xmlFileId, * ) "<metadata dictRef='",trim(name),"' content='",trim(value),"'></metadata>"

!   Use the xmlns field to specify that this is part of the cml namespace
!    write( xmlFileId, * ) "<metadata xmlns='http://www.xml-cml.org/schema' dictRef='",trim(name),"' content='",trim(value),"'></metadata>"

!    if ( openedCML ) call xmlEndCML()

  end subroutine xmlAddExtraMetadata

  subroutine xmlAddParameterList( title )
    use writexmli

    CHARACTER(LEN=*), INTENT(IN) :: title
    
    logical opg_root
    if ( .not. opg_root() ) return

    call xmlIndenter()
    write( xmlFileId, * ) '<parameterList title="',trim(title),'">'
    call xmlIncreaseIndent()

  end subroutine xmlAddParameterList

  subroutine xmlEndParameterList
    use writexmli
    logical opg_root
    if ( .not. opg_root() ) return
    call xmlDecreaseIndent()
    call xmlIndenter()
    write( xmlFileId, * )"</parameterList>"
  end subroutine xmlEndParameterList

  subroutine xmlAddParameterStr( name, value )
    ! Add a parameter as a name-value, where value is a string
    use writexmli
    CHARACTER(LEN=*), INTENT(IN) :: name, value
    
    logical opg_root
    if ( .not. opg_root() ) return

    call xmlIndenter()
    write( xmlFileId, * ) '<parameter name="',trim(name),'" value="',trim(value),'"/>'
  end subroutine xmlAddParameterStr

  subroutine xmlAddParameterInt( name, value )
    ! Add a parameter as a name-value, where value is an integer
    use writexmli
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: value
    CHARACTER(LEN=30) :: strvalue
    
    logical opg_root
    if ( .not. opg_root() ) return

    call xmlIndenter()
    call xmlNumberToString( value, strvalue )
    write( xmlFileId, * ) '<parameter name="',trim(name),'" value="',trim(strvalue),'"/>'
  end subroutine xmlAddParameterInt


!  subroutine xmlEndParameter
!    use writexmli
!    call xmlDecreaseIndent()
!    call xmlIndenter()
!    write( xmlFileId, * )"</parameter>"
!  end subroutine xmlEndParameter

  subroutine xmlWriteModuleSCF( ConvThresh, ConvIndex, Maxcyc, DmpCutOff, &
       &                         DiisOn, DiisOff )
  ! Write out all of the information on an scf convergenve

    use writexmli
    
    INTEGER, INTENT(IN) :: ConvThresh, ConvIndex, MaxCyc
    REAL(kind=xmlDP), INTENT(IN) :: DmpCutOff, DiisOn, DiisOff

    logical :: opg_root ! for parallel case
    if (.not. opg_root() ) return

    if( xmlStatus.eq.0 ) call xmlStartWriter()

    call xmlStartSCF()
    call xmlAddPropertyList1( "SCF parameters" )

    call xmlAddProperty1( "gamessuk:SCFConvThresh" )
    call xmlAddDPScalar( 10.0d0**(-ConvThresh) )
    call xmlEndProperty()

    call xmlAddProperty1( "gamessuk:SCFConvIndex" )
    call xmlAddDPScalar( 10.0d0**(-ConvIndex) )
    call xmlEndProperty()

    call xmlAddProperty1( "gamessuk:SCFMaxCyc" )
    call xmlAddInteger( MaxCyc )
    call xmlEndProperty()

    call xmlAddProperty1( "gamessuk:SCFDampingCutOff" )
    call xmlAddDPScalar( DmpCutOff )
    call xmlEndProperty()

    call xmlAddProperty1( "gamessuk:SCFDiisOnset" )
    call xmlAddDPScalar( DiisOn )
    call xmlEndProperty()

    call xmlAddProperty1( "gamessuk:SCFDiisOff" )
    call xmlAddDPScalar( DiisOff )
    call xmlEndProperty()

    call xmlEndPropertyList()
!    call xmlEndModule()

  end subroutine xmlWriteModuleSCF    

  subroutine xmlWriteModuleSCFCycle( Cycle, TotalEnergy, ElectronicEnergy, Econv, Tester )
  ! Write out all of the information on an scf convergenve

    use writexmli
    
    INTEGER, INTENT(IN) :: Cycle
    REAL(kind=xmlDP), INTENT(IN) :: TotalEnergy, ElectronicEnergy, Econv, Tester

    logical :: opg_root ! for parallel case
    if (.not. opg_root() ) return

    if( xmlStatus.eq.0 ) call xmlStartWriter()

    call xmlAddModuleSCFCycle( Cycle )
    call xmlAddPropertyList1( "SCF Cycle" )

    call xmlAddProperty1( "gamessuk:SCFCycTotalEnergy" )
    call xmlAddDPScalar( TotalEnergy )
    call xmlEndProperty()

    call xmlAddProperty1( "gamessuk:SCFCycElectronicEnergy" )
    call xmlAddDPScalar( ElectronicEnergy )
    call xmlEndProperty()

    call xmlAddProperty1( "gamessuk:SCFCycEconv" )
    call xmlAddDPScalar( Econv )
    call xmlEndProperty()

    call xmlAddProperty1( "gamessuk:SCFCycTester" )
    call xmlAddDPScalar( Tester )
    call xmlEndProperty()

    call xmlEndPropertyList()
    call xmlEndModule()

  end subroutine xmlWriteModuleSCFCycle

  subroutine xmlSCFFinalEnergies( ElectronicEnergy, NuclearEnergy, TotalEnergy)

    ! Write out the final SCF values and close the SCF module
    use writexmli

    REAL(kind=xmlDP), INTENT(IN) :: ElectronicEnergy, NuclearEnergy, TotalEnergy

    logical :: opg_root ! for parallel case
    if (.not. opg_root() ) return

    if (.not.xmlModuleSCF) call xmlStartSCF()

    call xmlAddPropertyList1( "SCF Final Properties" )

    call xmlAddProperty1( "gamessuk:SCFElectronicEnergy" )
    call xmlAddDPScalar( ElectronicEnergy )
    call xmlEndProperty()

    call xmlAddProperty1( "gamessuk:SCFNuclearEnergy" )
    call xmlAddDPScalar( NuclearEnergy )
    call xmlEndProperty()

    call xmlAddProperty1( "gamessuk:SCFTotalEnergy" )
    call xmlAddDPScalar( TotalEnergy )
    call xmlEndProperty()

    call xmlEndPropertyList()
    call xmlEndSCF()

  end subroutine xmlSCFFinalEnergies


  subroutine xmlAddPropertyList1( title )
    use writexmli
    CHARACTER(LEN=*), INTENT(IN) :: title
    
    logical opg_root
    if ( .not. opg_root() ) return

    call xmlIndenter()
    write( xmlFileId, * ) '<propertyList title="',trim(title),'">'
    call xmlIncreaseIndent()
  end subroutine xmlAddPropertyList1

  subroutine xmlEndPropertyList
    use writexmli

    logical opg_root
    if ( .not. opg_root() ) return

    call xmlDecreaseIndent()
    call xmlIndenter()
    write( xmlFileId, * )"</propertyList>"
  end subroutine xmlEndPropertyList

  subroutine xmlAddProperty1( dictRef )
    use writexmli
    CHARACTER(LEN=*), INTENT(IN) :: dictRef
    
    logical opg_root
    if ( .not. opg_root() ) return

    call xmlIndenter()
    write( xmlFileId, * ) '<property dictRef="',trim(dictRef),'">'
    call xmlIncreaseIndent()
  end subroutine xmlAddProperty1

  subroutine xmlEndProperty
    use writexmli

    logical opg_root
    if ( .not. opg_root() ) return

    call xmlDecreaseIndent()
    call xmlIndenter()
    write( xmlFileId, * )"</property>"
  end subroutine xmlEndProperty

  subroutine xmlAddInteger( value )

    use writexmli

    INTEGER, INTENT(in) :: value
    character(len=30) :: valueStr

    logical :: opg_root ! for parallel case
    if (.not. opg_root() ) return

    if( xmlStatus.eq.0 ) call xmlStartWriter()

    call xmlIndenter()
    call xmlNumberToString( value, valueStr )
    write( xmlFileId, * ) '<scalar datatype="xsd:integer">',trim(valueStr),'</scalar>'
    
  end subroutine xmlAddInteger

  subroutine xmlAddDPScalar( value )

    use writexmli

    REAL(kind=xmlDP), INTENT(in) :: value
    character(len=30) :: valueStr

    logical :: opg_root ! for parallel case
    if (.not. opg_root() ) return

    if( xmlStatus.eq.0 ) call xmlStartWriter()

    call xmlIndenter()
    call xmlNumberToString( value, valueStr )
    write( xmlFileId, * ) '<scalar datatype="xsd:double">',trim(valueStr),'</scalar>'
    
  end subroutine xmlAddDPScalar

  subroutine xmlAddDPScalarDR( dictRef, value )

    use writexmli

    character(len=*), intent(in) :: dictRef
    REAL(kind=xmlDP), intent(in) :: value
    character(len=30) :: valueStr

    logical :: opg_root ! for parallel case
    if (.not. opg_root() ) return

    if( xmlStatus.eq.0 ) call xmlStartWriter()

!    call xmlStartCML()
    call xmlIndenter()

    call xmlNumberToString( value, valueStr )
    
    write( xmlFileId, * ) "<scalar dictRef='",trim(dictRef),"'>",trim(valueStr),"</scalar>"
!    call xmlEndCML()
    
  end subroutine xmlAddDPScalarDR
  
  subroutine xmlAddError( desc, code )
    
    use writexmli

    character(len=*), intent(in) :: desc
    character(len=30) :: codeStr
    integer :: code

    logical :: opg_root ! for parallel case
    if (.not. opg_root() ) return
    
    call xmlNumberToString( code, codeStr )

    call xmlAddExtraMetadata( "gamess-uk:errorDesc", desc )
    call xmlAddExtraMetadata( "gamess-uk:errorCode", trim(codeStr) )
    
  end subroutine xmlAddError
