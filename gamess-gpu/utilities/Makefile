#
CC=cc
#
LDLIBS=-lm
CFLAGS=-g -D_POSIX_C_SOURCE=200112L -D_SYSTYPE_SVR4

all: wrap72 validate ccnvrt quote timer

wrap72: wrap72.o
	$(CC) -o wrap72 $(CFLAGS) wrap72.o $(LDLIBS)
validate:	validate.o
	$(CC) -o validate $(CFLAGS) validate.o $(LDLIBS)
ccnvrt:	ccnvrt.o
	$(CC) -o ccnvrt $(CFLAGS) ccnvrt.o $(LDLIBS)
quote:	quote.o
	$(CC) -o quote $(CFLAGS) quote.o $(LDLIBS)
timer:	timer.o
	$(CC) -o timer $(CFLAGS) timer.o ${LDLIBS}

clean:
	touch x.o
	rm -rf wrap72 validate ccnvrt quote timer *.o *.exe

MACH=rs6000.m
M4=m4
M4OPTS=-DM4_OPTIONS=beta
#M4OPTS=-Dalpha
test:
	cat ../machines/$(MACH) ../utilities/gener.m test.m | $(M4) $(M4OPTS) - 

# DO NOT DELETE
