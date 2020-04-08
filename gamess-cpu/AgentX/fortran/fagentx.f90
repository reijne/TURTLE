
integer function axvalueconvertint(value)
  INTEGER, INTENT(INOUT) :: value
  CHARACTER(LEN=30) vstring
  call axvalue(vstring)
  read(vstring,*)value
  axvalueconvertint=0
  return
end function axvalueconvertint

integer function axvalueconvertfloat(value)
  INTEGER, parameter :: prec = KIND(1.0)
  REAL(KIND=prec), INTENT(INOUT) :: value
  CHARACTER(LEN=30) vstring
  call axvalue(vstring)
  read(vstring,*)value
  axvalueconvertfloat=0
  return
end function axvalueconvertfloat

integer function axvalueconvertdouble(value)
  INTEGER, parameter :: prec = KIND(1D0)
  REAL(KIND=prec), INTENT(INOUT) :: value
  CHARACTER(LEN=30) vstring
  call axvalue(vstring)
  read(vstring,*)value
  axvalueconvertdouble=0
  return 
end function axvalueconvertdouble

INTEGER FUNCTION axproxy (str)
  IMPLICIT NONE
  INTEGER :: strlen = 0, axwrapproxy
  CHARACTER(LEN=*), INTENT(IN) :: str
  strlen = len(str)
  axproxy = axwrapproxy(strlen,str)
  return
END FUNCTION axproxy

INTEGER FUNCTION axgeturi (str)
  IMPLICIT NONE
  INTEGER :: strlen = 0, axwrapgeturi
  CHARACTER(LEN=*), INTENT(IN) :: str
  strlen = len(str)
  axgeturi = axwrapgeturi(strlen,str)
  return
END FUNCTION axgeturi

INTEGER FUNCTION axdatageturi (str)
  IMPLICIT NONE
  INTEGER :: strlen = 0, axwrapdatageturi
  CHARACTER(LEN=*), INTENT(IN) :: str
  strlen = len(str)
  axdatageturi = axwrapdatageturi(strlen,str)
  return
END FUNCTION axdatageturi

INTEGER FUNCTION axselect (str)
  IMPLICIT NONE
  INTEGER :: strlen = 0, axwrapselect
  CHARACTER(LEN=*), INTENT(IN) :: str
  strlen = len(str)
  axselect = axwrapselect(strlen,str)
  return
END FUNCTION axselect

INTEGER FUNCTION axselectnext()
  INTEGER :: axwrapselectnext
  axselectnext = axwrapselectnext()
  return
END FUNCTION axselectnext

INTEGER FUNCTION axdeselect()
  INTEGER :: axwrapdeselect
  axdeselect = axwrapdeselect()
  return
END FUNCTION axdeselect

INTEGER FUNCTION axbaseuri(str)
  INTEGER :: strlen = 0, axwrapbaseuri
  CHARACTER(LEN=*), INTENT(IN) :: str
  strlen = len(str)
  axbaseuri = axwrapbaseuri(strlen,str)
  return
END FUNCTION axbaseuri

INTEGER FUNCTION axparserfinish()
  INTEGER :: axwrapparserfinish
  axparserfinish = axwrapparserfinish()
  return
END FUNCTION axparserfinish

INTEGER FUNCTION axparserstart()
  INTEGER :: axwrapparserstart
  axparserstart = axwrapparserstart()
  return
END FUNCTION axparserstart

INTEGER FUNCTION axrefine(str)
  IMPLICIT NONE
  INTEGER :: strlen = 0, axwraprefine
  CHARACTER(LEN=*), INTENT(IN) :: str
  strlen = len(str)
  axrefine = axwraprefine(strlen,str)
  return
END FUNCTION axrefine

INTEGER FUNCTION axformat(dataformat)
  IMPLICIT NONE
  INTEGER, INTENT(IN)::dataformat
  INTEGER :: axwrapformat
  axformat = axwrapformat(dataformat)
  return
END FUNCTION axformat

INTEGER FUNCTION axcurrent()
  INTEGER :: axwrapcurrent
  axcurrent = axwrapcurrent()
  return
END FUNCTION axcurrent

INTEGER FUNCTION axlen()
  INTEGER :: axwraplen
  axlen = axwraplen()
  return
END FUNCTION axlen

SUBROUTINE axvalue(str)
  IMPLICIT NONE
  INTEGER :: strlen = 0
  CHARACTER(LEN=*), INTENT(INOUT) :: str
  strlen = len(str)
  call axwrapvalue(strlen,str)
  return
END SUBROUTINE axvalue

INTEGER FUNCTION axbuffer(size)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: size
  INTEGER :: axwrapbuffer
  axbuffer = axwrapbuffer(size)
  return
END FUNCTION axbuffer

INTEGER FUNCTION axnamespace(str)
  IMPLICIT NONE
  INTEGER :: strlen = 0, axwrapnamespace
  CHARACTER(LEN=*), INTENT(IN) :: str
  strlen = len(str)
  axnamespace = axwrapnamespace(strlen,str)
  return
END FUNCTION axnamespace
      
INTEGER FUNCTION axcache(size)
  IMPLICIT NONE
  INTEGER, INTENT(IN)::size
  INTEGER :: axwrapcache
  axcache = axwrapcache(size)
  return
END FUNCTION axcache

INTEGER FUNCTION axselectno(position)
  IMPLICIT NONE
  INTEGER, INTENT(IN)::position
  INTEGER :: axwrapselectno
  axselectno = axwrapselectno(position)
  return
END FUNCTION axselectno
