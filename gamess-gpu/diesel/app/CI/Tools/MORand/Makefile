QMLIBS = \
	../../../../lib/QM/libQM.a \
	../../../../lib/Container/libContainer.a \
	../../../../lib/Math/libMath.a \
	../../../../../m4/libdiesel.a

APP = morand

-include ../../../../Makefile.conf
include ../../../../Makefile.common

main.o:	main.cc $(QMLIB) Compiled.h ../../../../Configured.h

$(APP): 	$(QMLIBS) $(OBJECTS) main.o 
	$(LD) $(LDFLAGS) -o $@ main.o \
	$(OBJECTS) $(QMLIBS) $(STDLIBS)
