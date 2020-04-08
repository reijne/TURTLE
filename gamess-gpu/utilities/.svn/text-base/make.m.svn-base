#
_IF(hp700)
LDLIBS=-lm
_ELSE
LDLIBS=
_ENDIF
_IF(rs6000)
CFLAGS=-g
_ELSE
CFLAGS=
_ENDIF
all: wrap72 validate
wrap72: wrap72.o
	cc -o wrap72 $(CFLAGS) wrap72.o
validate:	validate.o
	cc -o validate $(CFLAGS) validate.o $(LDLIBS)
