QMLIBS = \
	../../../../lib/QM/libQM.a \
	../../../../lib/Container/libContainer.a \
	../../../../lib/Math/libMath.a \
	../../../common/libCommon.a \
	../../../../../m4/libdiesel.a

STDLIBS = -lfl

APP = diagExt

-include ../../../../Makefile.conf
include ../../../../Makefile.common

main.o:	main.cc $(QMLIBS) Compiled.h Configured.h

$(APP): 	$(QMLIBS) $(OBJECTS) main.o
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS) main.o \
	../Input.o ../DiagInput.o \
	$(QMLIBS) $(STDLIBS)




