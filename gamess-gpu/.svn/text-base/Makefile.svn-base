# A simple Makefile to call the main Makefile in m4

all: check_config
	@echo " The build will take place in the m4 directory.";echo
	@cd m4; make

check_config:
	@(if test -f m4/Makefile; then \
	echo ; echo " Found a Makefile in m4, using this to build GAMESS-UK." ; echo; \
	else \
	echo; \
	echo " ###################################################################" ; \
	echo " #  ERROR! No Makefile found in the m4 directory!                  #" ; \
	echo " #  You need to run the configure script in this directory first.  #" ; \
	echo " #  For information on configuring and building GAMESS-UK, please  #" ; \
	echo " #  read the INSTALL documentation in the config directory.        #" ; \
	echo " ###################################################################" ; \
	echo; exit -1 ; fi)

clean:  check_config
	@(cd m4; make clean)

realclean: check_config
	@(cd m4; make realclean)
	(cd utilities; make clean; make)
	(cd m4; make depend)

squeakyclean: check_config
	@(cd m4; make clean)
	find . -name "*.o" -exec /bin/rm {} \;
	find . -name "*.a" -exec /bin/rm {} \;
	find . -name "dependencies.txt" -exec /bin/rm {} \;
	(cd utilities; make clean; make)
	(cd m4; make depend)
