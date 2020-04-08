
INCLUDES 		= -I@PERL5EXT@
lib_LTLIBRARIES 	= libplagentx.la

libplagentx_la_SOURCES	= plagentx.c
libplagentx_la_LDFLAGS	= -module -avoid-version @CYGWIN_EXTRA_LDFLAGS@ @CYGWIN_EXTRA_PERL_LDFLAGS@
libplagentx_la_LIBADD	= $(top_builddir)/src/libagentx.la @CYGWIN_EXTRA_PERL_LIBADD@

install-data-local:
			$(mkinstalldirs) "$(DESTDIR)$(prefix)/share/AgentX/perl"
			$(INSTALL_DATA) libplagentx.pm "$(DESTDIR)$(prefix)/share/AgentX/perl"
			ln -fs "$(prefix)/lib/libplagentx.so" "$(DESTDIR)$(prefix)/lib/libplagentx.bundle"

