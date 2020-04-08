-include ../../../Makefile.conf

LIB = libStream.a


sourcesCC =	   binstream.cc

OBJECTS =	${sourcesCC:.cc=.o} ${sourcesC:.c=.o}


all: $(LIB)


main: $(LIB)
	g++ $(LIB) -o main

$(LIB):  $(OBJECTS) 
	ar qv $(LIB) ;\
	@set -e; for i in $(SUBDIRS); do \
		ar r $(LIB) $$i/*.o; \
	done
	ar r $(LIB) $(OBJECTS) ;\



ifneq ($(CLEAN),clean)
-include $(patsubst %.cc,.%.d,$(sourcesCC))
endif



.%.d: %.cc
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed '\''s/$*.o/& $@/g'\'' > $@'


clean:
	@set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i clean; done
	rm -f *.o .*.d *.bak *~ main $(LIB)
