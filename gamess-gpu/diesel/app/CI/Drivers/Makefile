SUBDIRS = $(dir $(wildcard */Makefile))

all:
	@set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i; done

clean:
	@for i in $(SUBDIRS); do $(MAKE) -C $$i clean CLEAN=clean; done
