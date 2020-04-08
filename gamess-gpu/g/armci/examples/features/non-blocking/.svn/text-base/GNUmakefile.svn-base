#$Id: GNUmakefile,v 1.1.1.1 2007-10-30 10:13:45 jmht Exp $
all: simple overlap

.PHONY : simple
simple:
	@echo '>>> Making simple only <<<'
	(cd simple;make)

.PHONY : overlap
overlap:
	@echo '>>> Making overlap only <<<'
	(cd overlap;make)

clean:
	(cd simple;make clean)
	(cd overlap;make clean)
