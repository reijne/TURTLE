
SUBDIRS = Container Math QM

all:
	@set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i; done

clean:
	for i in $(SUBDIRS); do $(MAKE) -C $$i clean CLEAN=clean; done
