SUBDIRS = $(dir $(wildcard */Makefile))

#SUBDIRS = \
#ConfMap \
#GeomOpt \
#MORand \
#PrintDens \
#RefSym \
#SelSym \
#SetOps \
#cmpPT \
#confStat \
#diffConf \
#efTools \
#excited \
#lstConfs \
#lstInts \
#orb2xfig \
#prop2xfig \
#symTree


all:
	@set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i; done

clean:
	@set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i clean CLEAN=clean; done
