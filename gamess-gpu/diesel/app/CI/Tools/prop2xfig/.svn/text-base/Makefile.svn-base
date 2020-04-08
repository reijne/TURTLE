QMLIBS = \
	../../../../lib/QM/libQM.a \
	../../../../lib/Container/libContainer.a \
	../../../../lib/Math/libMath.a \
	../../../../../m4/libdiesel.a

-include ../../../../Makefile.conf
include ../../../../Makefile.common

APP =	${SOURCESCC:.cc=}
all:	$(APP)


%:	%.o
	$(LD) $(LDFLAGS) -o $@ $@.o $(QMLIBS) $(STDLIBS)
