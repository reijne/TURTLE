QMLIBS = \
	../../../../lib/QM/libQM.a \
	../../../../lib/Container/libContainer.a \
	../../../../lib/Math/libMath.a \
	../../../common/libCommon.a \
	../../../../../m4/libdiesel.a

STDLIBS = -lfl

INCS += -I/usr/include

APP = diesel

-include ../../../../Makefile.conf
include ../../../../Makefile.common

main.o:	main.cc $(QMLIBS) *.h Compiled.h ../../../../Configured.h

$(APP): 	$(QMLIBS) $(OBJECTS) main.o
	$(LD) $(LDFLAGS) -o $@ main.o $(OBJECTS) $(QMLIBS) $(STDLIBS)


