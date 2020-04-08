QMLIBS = \
        ../../../../lib/QM/libQM.a \
        ../../../../lib/Container/libContainer.a \
        ../../../../lib/Math/libMath.a \
	../../../../../m4/libdiesel.a

APP = propints

-include ../../../../Makefile.conf
include ../../../../Makefile.common


%.o: %.f
	g77 -c $<


$(APP):  $(OBJECTS) read.o
	gcc $(LDFLAGS) -o $@ \
	$(OBJECTS) $(STDLIBS) -lf2c




