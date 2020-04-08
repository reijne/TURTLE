
SUFFIXES = .o .f90
.f90.o:
	@FC@ @FCFLAGS@ @FCFLAGS_free_f90@ @DEFAULT_INCLUDES@ -c $<

lib_LTLIBRARIES 	= libfagentx.la

libfagentx_la_SOURCES	= fagentxi.f90 fagentxc.c fagentx.f90 $(top_builddir)/include/fagentxInternals.h $(top_builddir)/include/agentx.h $(top_builddir)/include/config.h 
libfagentx_la_LDFLAGS	= -static -release @PACKAGE_VERSION@

clean-local:
			rm -f @FMOD_NAME@

maintainer-clean-local:
			rm -f @FMOD_NAME@

install-data-local:
			$(mkinstalldirs) "$(DESTDIR)$(prefix)/share/AgentX/fortran"
			$(INSTALL_DATA) @FMOD_NAME@ "$(DESTDIR)$(prefix)/share/AgentX/fortran/"
