
lib_LTLIBRARIES 	= libagentx.la
include_HEADERS		= $(top_builddir)/include/agentx.h

bin_SCRIPTS		= $(top_builddir)/agentx-config $(top_builddir)/axcshenv $(top_builddir)/axbashenv
libagentx_la_SOURCES 	= axPath.c agentx.c eval.c iniparser.c dictionary.c strlib.c axError.c axArray.c $(top_builddir)/include/dictionary.h $(top_builddir)/include/iniparser.h $(top_builddir)/include/strlib.h $(top_builddir)/include/agentxInternals.h $(top_builddir)/include/agentx.h $(top_builddir)/include/axArray.h $(top_builddir)/include/config.h $(top_builddir)/include/axArray.h $(top_builddir)/include/axError.h
libagentx_la_LDFLAGS	= -release @PACKAGE_VERSION@ @CYGWIN_EXTRA_LDFLAGS@
libagentx_la_LIBADD	= @libxml2_libs@


install-data-local:
			$(mkinstalldirs) "$(DESTDIR)$(prefix)/share/AgentX/ontology"
			$(INSTALL_DATA) "$(top_builddir)/ontology/ontology.owl" "$(DESTDIR)$(prefix)/share/AgentX/ontology"
			$(INSTALL_DATA) "$(top_builddir)/ontology/ontology.pprj" "$(DESTDIR)$(prefix)/share/AgentX/ontology"
			$(mkinstalldirs) "$(DESTDIR)$(prefix)/share/AgentX/map"
			$(INSTALL_DATA) "$(top_builddir)/map/map.rdf" "$(DESTDIR)$(prefix)/share/AgentX/map"

uninstall-local:
			rm -rf "$(DESTDIR)$(prefix)/share/AgentX"
