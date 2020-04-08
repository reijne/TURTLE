-include ../../../../Makefile.conf
LIB=../../libQM.a

k := $(shell awk '$$1=="\#define" && $$2=="SIZEOF_VOID_P" { printf("        integer*%d pUMat, pUMat1, pEMat, php5dar\n",$$3); }' ../../../../config.h >pointers.inc)

#FCOMPILER = f77
#FFLAGS = -O
#CC = gcc


SourcesF =	bum03.f       bum12.f       calcbum.f     cycle2.f\
			bum01.f       bum04.f       bum13.f       calcdarhp5.f\
			bum02.f       bum11.f       bum21.f       cycle1.f      cycle.f \
			buminit.f\
			sga.f branch.f weight.f eltran.f

OBJECTS =	${SourcesF:.f=.o} 


main:   $(LIB) 

  
$(LIB):  $(OBJECTS) loc.o
	ar $(ARFLAGS) $(LIB) $(OBJECTS) loc.o ;\


loc.o:	loc.c
	$(CC) $(CFLAGS) -c loc.c -o loc.o


.f.o: global sga.inc sgapars.inc $<
	$(FCOMPILER) $(FFLAGS) -c $<

$(OBJECTS): global sga.inc sgapars.inc

clean:
	rm -f *.o *~ *% lib*.a pointers.inc

realclean:
	rm -f *.o *~ *% lib*.a pointers.inc
