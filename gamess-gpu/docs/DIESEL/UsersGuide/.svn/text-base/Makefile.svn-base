SUBDIRS = $(dir $(wildcard */Makefile))


TEXDEPFILES = toc lof lot bbl


MAINDEP = $(wildcard *.tex */*.bib */*.tex */*/*.tex */*.pstex.ps */*.pstex */*/*.pstex.ps */*/*.pstex */*/*.ps */*/*/*.tex */*/*/*.pstex.ps */*/*/*.pstex */*/*/*.ps)

SOURCESFIG =	$(wildcard [^~]*.fig)
SOURCESGNUPLOT =	$(wildcard [^~]*.gnuplot)
SOURCESPS =	$(wildcard [^~]*.ps.source)
SOURCESTIFF =	$(wildcard [^~]*.tif)

OBJECTS =	\
	${SOURCESGNUPLOT:.gnuplot=.ps} \
	${SOURCESFIG:.fig=.pstex} ${SOURCESFIG:.fig=.pstex.ps} \
	${SOURCESTIFF:.tif=.ps} \
	${SOURCESPS:.ps.source=.ps} \


DEPENDENCIES = SubDirMake $(OBJECTS)

MAIN = $(wildcard main.tex)

	
ifneq ($(MAIN),)
	DEPENDENCIES += main.dvi
endif


all: $(DEPENDENCIES)


main.dvi: $(MAINDEP) flag
	echo "\centerline{"Version: `date +%H`":"`date +%M`", "`date +%d`"."`date +%m`"."`date +%Y`"}" > Date.tex
	@for i in $(TEXDEPFILES); do touch main.$$i; cp main.$$i main.$$i.old; done
	@if ! latex main; then rm main.dvi; exit; fi
#	bibtex main
	@sleep 1
	@for i in $(TEXDEPFILES); do if ! diff -q main.$$i main.$$i.old >/dev/null; then touch flag; fi; done
	@rm main.*.old
	@if grep "Rerun to" main.log >/dev/null; then touch flag; fi;
	@make main.dvi


main.ps: main.dvi 
	dvips main.dvi -o main.ps

main.1200dpi.ps: main.dvi 
	dvips -Plex24 main.dvi -o main.1200dpi.ps

main.pcl: main.ps
	gs -sDEVICE=ljet3 -r600 -sPAPERSIZE=a4 -dNOPAUSE -dSAFER -sOutputFile=main.pcl -q main.ps <<"quit"


%.ps: %.gnuplot %.data Makefile
	cd $(dir $<);gnuplot $< >$@

%.gnuplot:
	
%.ps: %.ps.source Makefile
	ln -sf $(notdir $*).ps.source $(notdir $*).ps


%.data:
	
%.ps: %.tif Makefile
	tifftopnm $< | pnmtops -scale 1.0 >$@


%.pstex.ps: %.fig %.gnuplot Makefile
	fig2dev -Lpstex $< >$@

%.pstex: %.fig %.gnuplot Makefile
	fig2dev -Lpstex_t -p `pwd`/$@.ps $< >$@


#-------------------------------------------------------------------------------

SubDirMake:
	@set -e; if [ "$(SUBDIRS)" != "" ]; then make _SubDirMake; fi
	
	
_SubDirMake:
	@for i in $(SUBDIRS); do	$(MAKE) -C $$i; done
	

clean:
	@set -e; if [ "$(SUBDIRS)" != "" ]; then make _SubDirClean; fi
	@rm -f *.pcl *.ps *.pstex* *.dvi *.lof *.log *.lot *.toc *.aux *~ *.bck *.bak *.iso main.bbl main.blg main.lop
	@find -name "*.bck" -exec rm {} \;;

_SubDirClean:
	@for i in $(SUBDIRS); do	$(MAKE) clean -C $$i; done



tar:
	make clean
	tar -czvf ../Diss.tgz *
	cp ../Diss.tgz ../Diss.Versions/Diss.`hostname`.`date +%y%m%d.%H%M`.tgz
