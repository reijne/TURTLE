
SUBDIRS = common CI Interfaces

all:
	@set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i; done

clean:
	@set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i clean CLEAN=clean; done
