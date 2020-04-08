MODULE fagentxi

  INTERFACE

     INTEGER FUNCTION axproxy (proxystr)
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN)::proxystr
     END FUNCTION axproxy

     INTEGER FUNCTION axgetUri (uri)
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN)::uri
     END FUNCTION axgeturi

     INTEGER FUNCTION axdatageturi (uri)
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN)::uri
     END FUNCTION axdatageturi

     INTEGER FUNCTION axselect (entity)
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN)::entity
     END FUNCTION axselect

     INTEGER FUNCTION axselecttype(entity,type)
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN)::entity
       INTEGER, INTENT(IN)::type
     END FUNCTION axselecttype

     INTEGER FUNCTION axselectnext()
     END FUNCTION axselectnext
     
     INTEGER FUNCTION axdeselect()
     END FUNCTION axdeselect

     INTEGER FUNCTION axbaseuri(baseuristring)
       CHARACTER(LEN=*), INTENT(IN)::baseuristring
     END FUNCTION axbaseuri
  
     INTEGER FUNCTION axparserfinish()
     END FUNCTION axparserfinish

     INTEGER FUNCTION axparserstart()
     END FUNCTION axparserstart
     
     INTEGER FUNCTION axrefine(expression)
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN)::expression
     END FUNCTION axrefine

     INTEGER FUNCTION axformat(dataformat)
       IMPLICIT NONE
       INTEGER, INTENT(IN)::dataformat
     END FUNCTION axformat

     INTEGER FUNCTION axcurrent()
     END FUNCTION axcurrent

     INTEGER FUNCTION axlen()
     END FUNCTION axlen

     SUBROUTINE axvalue(strvar)
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(OUT)::strvar
     END SUBROUTINE axvalue
       
     INTEGER FUNCTION axbuffer(size)
       IMPLICIT NONE
       INTEGER, INTENT(IN)::size
     END FUNCTION axbuffer

     INTEGER FUNCTION axnamespace(namespaceuri)
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN)::namespaceuri
     END FUNCTION axnamespace

     INTEGER FUNCTION axcache(size)
       IMPLICIT NONE
       INTEGER, INTENT(IN)::size
     END FUNCTION axcache

     INTEGER FUNCTION axselectno(position)
       IMPLICIT NONE
       INTEGER, INTENT(IN)::position
     END FUNCTION axselectno

  END INTERFACE

  INTERFACE axvalueconvert
     
     module procedure axvalueconvertfloatdum
     module procedure axvalueconvertdoubledum
     module procedure axvalueconvertintdum
 
  END INTERFACE

CONTAINS

  integer function axvalueconvertintdum(value)
    INTEGER, INTENT(IN) :: value
    axvalueconvertintdum = axvalueconvertint(value)
    return
  end function axvalueconvertintdum
  
  integer function axvalueconvertfloatdum(value)
    INTEGER, parameter :: prec = KIND(1.0)
    REAL(KIND=prec), INTENT(IN) :: value
    axvalueconvertfloatdum = axvalueconvertfloat(value)
    return
  end function axvalueconvertfloatdum
  
  integer function axvalueconvertdoubledum(value)
    INTEGER, parameter :: prec = KIND(1D0)
    REAL(KIND=prec), INTENT(IN) :: value
    axvalueconvertdoubledum = axvalueconvertdouble(value)
    return 
  end function axvalueconvertdoubledum
  
END MODULE
