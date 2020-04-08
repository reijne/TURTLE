
INCLUDES 		= @JINCLUDE@
lib_LTLIBRARIES 	= libjagentx.la

libjagentx_la_SOURCES	= jagentx.c
libjagentx_la_LDFLAGS	= -module -avoid-version @CYGWIN_EXTRA_LDFLAGS@ @CYGWIN_EXTRA_JAVA_LDFLAGS@
libjagentx_la_LIBADD	= $(top_builddir)/src/libagentx.la @CYGWIN_EXTRA_JAVA_LIBADD@


install-data-local:
			$(mkinstalldirs) "$(DESTDIR)$(prefix)/share/AgentX/java"
			$(INSTALL_DATA) libjagentx.java  "$(DESTDIR)$(prefix)/share/AgentX/java"
			$(INSTALL_DATA) libjagentxJNI.java "$(DESTDIR)$(prefix)/share/AgentX/java"
			$(INSTALL_DATA) SWIGTYPE_p_FILE.java "$(DESTDIR)$(prefix)/share/AgentX/java"
			$(INSTALL_DATA) SWIGTYPE_p_p_char.java "$(DESTDIR)$(prefix)/share/AgentX/java"
			$(INSTALL_DATA) SWIGTYPE_p_p_p_char.java "$(DESTDIR)$(prefix)/share/AgentX/java"