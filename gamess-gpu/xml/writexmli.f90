MODULE writexmli
  
  INTEGER, parameter :: xmlSP = KIND( 1.0 ), xmlDP = KIND( 1D0 )
  INTEGER, save :: xmlIndent, xmlFileId, xmlStatus, xmlAtomId
  INTEGER, save :: xmlSCFCount ! To number the scf cycles
  LOGICAL, save :: xmlCML, xmlModuleSCF ! Track CML namespace & SCF modules
  data xmlStatus,xmlFileId,xmlIndent / 0,0,0/
  data xmlSCFCount,xmlModuleSCF,xmlCML / 0,.False.,.False. /

  INTERFACE xmlNumberToString
     
     module procedure xmlINumberToString
     module procedure xmlDPNumberToString

  END INTERFACE

CONTAINS
  
  subroutine xmlStartWriter()

    logical opg_root
    if ( .not. opg_root() ) return
    if( xmlStatus .ne. 0 ) return
    if( xmlFileId .eq. 0 ) xmlFileId = 59
    xmlIndent = 0
    xmlStatus = 1
    
    open( xmlFileId, file="gamessout.xml", err=10, status="replace", action="write" )
    write( xmlFileId, FMT="(A38)" ) "<?xml version='1.0' encoding='UTF-8'?>"
    write( xmlFileId, * ) ""
    call xmlIncreaseIndent()
    call xmlIndenter()
    write( xmlFileId, * ) "<eccp xmlns='http://www.grids.ac.uk/eccp/ns#' >"
    call xmlIncreaseIndent()
    write( xmlFileId, * ) ""

    ! For the time being assume everything is in CML space
    call xmlStartCML()
    return
10  write( *, * )"An error occured opening the file: out.xml"
    
  end subroutine xmlStartWriter
  
  subroutine xmlStartCML()
    
    call xmlIndenter()
    write( xmlFileId, * )"<cml xmlns='http://www.xml-cml.org/schema' >"
    call xmlIncreaseIndent()
    ! Flag we can query to check if we are in CML namespace
    xmlCML=.True.
    
  end subroutine xmlStartCML
  
  subroutine xmlEndCML()
    
    call xmlDecreaseIndent()
    call xmlIndenter()
    write( xmlFileId, * )"</cml>"
    xmlCML=.False.
    
  end subroutine xmlEndCML
  
  subroutine xmlStartMolecule( moleculeID, moleculeTitle )
    
    CHARACTER(LEN=*), INTENT(IN) :: moleculeID, moleculeTitle
    
!   Need to set xmlAtomId to 0 as otherwise each consecutive write
!   has differently identified atoms
    xmlAtomId = 0
    
    call xmlIndenter()
    write( xmlFileId, * ) "<molecule id='",trim( moleculeID ),"' title='",trim( moleculeTitle ),"' >"
    call xmlIncreaseIndent()
    
  end subroutine xmlStartMolecule
  
  subroutine xmlEndMolecule()
    call xmlDecreaseIndent()
    call xmlIndenter()
    write( xmlFileId, * )"</molecule>"
    
  end subroutine xmlEndMolecule
  
  subroutine xmlAddAtom( atomID, elementType, xCoord, yCoord, zCoord )
    
    CHARACTER(LEN=*), INTENT(IN) :: atomID
    REAL(kind=xmlDP), INTENT(IN) :: xCoord, yCoord, zCoord
    character(len=20) :: xmlAtomIdStr ! Used for numbering atoms (see xmlStartMolecule )
    character(len=20) :: xmlIdStr
    character(len=2)  :: elementType
    character(len=30) :: xCoordStr, yCoordStr, zCoordStr

    call xmlIndenter()

    call xmlNumberToString( xmlAtomId, xmlAtomIdStr )
    call xmlNumberToString( xCoord, xCoordStr )
    call xmlNumberToString( yCoord, yCoordStr ) 
    call xmlNumberToString( zCoord, zCoordStr )

    xmlIdStr = trim(atomID)//'_'//trim(xmlAtomIdStr)
    write( xmlFileId, * ) "<atom id='",trim(xmlIdStr),"' elementType='", &
         trim(elementType),"' x3='",trim(xCoordStr),"' y3='",trim(yCoordStr),"' z3='",trim(zCoordStr),"' />"
    
    xmlAtomId = xmlAtomId + 1

  end subroutine xmlAddAtom
  
  subroutine xmlStartSCF()
    
!   Adds the header for the SCF module and indents the writer
    character(len=5) :: strCount
    logical opg_root

    if (.not.opg_root()) return
    
    if ( xmlModuleSCF) then
       write(*,*) "*** xmlStartSCF adding module when already open! ***"
    endif

    ! Increment couter of SCF modules
    xmlSCFCount = xmlSCFCount + 1

    call xmlIndenter()
    call xmlNumberToString( xmlSCFCount, strCount )
    write( xmlFileId, * ) '<module dictRef="gamessuk:SCF" role="step" serial="',trim(strCount),'">'
    call xmlIncreaseIndent()

    ! Flag we are in an xml Module
    xmlModuleSCF = .True.
    
  end subroutine xmlStartSCF

  subroutine xmlEndSCF( )
    
    logical opg_root
    if (.not.opg_root()) return
    
    if ( .not. xmlModuleSCF) then
       write(*,*) "*** xmlEndSCF closing unopened module! ***"
    endif

    ! Flag we are in an xml Module
    xmlModuleSCF = .False.

    call xmlEndModule()
    
  end subroutine xmlEndSCF


  subroutine xmlAddModuleSCFCycle( serial )
    
!   Adds the header for the SCF module and indents the writer
    INTEGER, INTENT(IN) :: serial
    character(len=5) :: serialStr
    
    call xmlIndenter()
    call xmlNumberToString( serial, serialStr )
    write( xmlFileId, * ) '<module dictRef="gamessuk:SCF Cycle" role="step" serial="',trim(serialStr),'">'
    call xmlIncreaseIndent()
    
  end subroutine xmlAddModuleSCFCycle

  subroutine xmlEndModule()
    call xmlDecreaseIndent()
    call xmlIndenter()
    write( xmlFileId, * )"</module>"
  end subroutine xmlEndModule

  subroutine xmlAddMetadataList
    call xmlIndenter()
    write( xmlFileId, * ) '<metadataList>'
    call xmlIncreaseIndent()
  end subroutine xmlAddMetadataList

  subroutine xmlEndMetadataList
    call xmlDecreaseIndent()
    call xmlIndenter()
    write( xmlFileId, * )"</metadataList>"
  end subroutine xmlEndMetadataList


! Utility routines
  subroutine xmlDPNumberToString( value, strvalue )

    real(kind=xmlDP), intent(in) :: value
    character(len=30) :: tvalue
    character(len=*),intent(out) :: strvalue
    integer :: i,j

    write( tvalue, * ) value
    j = 1
    do i = 1,30
       if( tvalue(i:i) .ne. " " ) then
          strvalue(j:j) = tvalue(i:i)
          j = j + 1
       end if
    end do
   
    do i = j,30
       strvalue(i:i)=" "
    end do

  end subroutine xmlDPNumberToString

  subroutine xmlINumberToString( value, strvalue )

    integer, intent(in) :: value
    character(len=30) :: tvalue
    character(len=*),intent(out) :: strvalue
    integer :: i,j

    write( tvalue, * ) value
    j = 1
    do i = 1,len(tvalue)
       if( tvalue(i:i) .ne. " " ) then
          strvalue(j:j) = tvalue(i:i)
          j = j + 1
       end if
    end do
    
    do i = j,len(strvalue)
       strvalue(i:i)=" "
    end do

  end subroutine xmlINumberToString
         
  subroutine xmlIncreaseIndent()
    
    xmlIndent = xmlIndent + 2
    
  end subroutine xmlIncreaseIndent
  
  subroutine xmlDecreaseIndent()
    
    xmlIndent = xmlIndent - 2
    
  end subroutine xmlDecreaseIndent
  
  subroutine xmlIndenter()
    
    integer :: i
    
    i = 0
    do while ( i .le. xmlIndent )
       write( xmlFileId, FMT='(A1)', advance='no' )" "
       i = i + 1
    end do
    
  end subroutine xmlIndenter
  
end module writexmli
