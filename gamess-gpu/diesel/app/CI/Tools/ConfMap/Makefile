QMLIBS = \
	../../../../lib/QM/libQM.a \
	../../../../lib/Container/libContainer.a \
	../../../../lib/Math/libMath.a \
	../../../../../m4/libdiesel.a

APP = confMap

-include ../../../../Makefile.conf
include ../../../../Makefile.common

main.o:	main.cc $(QMLIBS)

$(APP): 	$(QMLIBS) $(OBJECTS) main.o
	$(LD) $(LDFLAGS) -o $@ main.o  $(OBJECTS) $(QMLIBS) $(STDLIBS)
