
_IF(alpha)
alpha
_ENDIF

_IF(beta)
beta
_ENDIF

_IF(gamma)
gamma
_ENDIF

_IF(_NOT(gamma,beta))
not (gamma or beta) (new)
_ENDIF

_IFN(gamma,beta)
not (gamma or beta) (old)
_ENDIF

_IF(_NOT(_AND(gamma,beta)))
(not gamma) or (not beta) (new)
_ENDIF

_IFN(gamma)
not gamma (old)
_ENDIF

_IF(beta,_NOT(alpha))
alpha or beta
_ENDIF


_IF(_AND(alpha,beta,gamma))
alpha and beta and gamma
_ENDIF

_IF(_AND(alpha,beta))
alpha and beta
_ENDIF

_IF(_AND(gamma,beta))
gamma and beta
_ENDIF

_IF(TRUE)
simple true
_ENDIF

_IF(alpha)
_MACRO(QM,qq(ivoff+`$1'))dnl
_ENDIF

      QM(1) = 1.0
