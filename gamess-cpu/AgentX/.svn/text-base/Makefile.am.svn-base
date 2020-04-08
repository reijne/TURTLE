
SUBDIRS  = @AXLIBXML2@ @AXCORE@ @AXFORTRAN@ @AXPYTHON@ @AXPERL@ @AXAXT@ @AXCEXA@ @AXFEXA@ @AXJAVA@

maintainer-clean-local:
			@find $(top_builddir) -type f -name "*~" -exec rm -f {} \;
			@find $(top_builddir) -type f -name "#*" -exec rm -f {} \;	
			@rm -f ontology/*.dot*
			@rm -f examples/python/*.pyc
