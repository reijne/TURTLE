
INCLUDES 		= @PYINCLUDE@
lib_LTLIBRARIES 	= _libpyagentx.la

_libpyagentx_la_SOURCES	= pyagentx.c
_libpyagentx_la_LDFLAGS	= -module -avoid-version @CYGWIN_EXTRA_LDFLAGS@ @CYGWIN_EXTRA_PYTHON_LDFLAGS@
_libpyagentx_la_LIBADD	= $(top_builddir)/src/libagentx.la @CYGWIN_EXTRA_PYTHON_LIBADD@


install-data-local:
			$(mkinstalldirs) "$(DESTDIR)$(prefix)/share/AgentX/python"
			$(INSTALL_DATA) libpyagentx.py "$(DESTDIR)$(prefix)/share/AgentX/python"
