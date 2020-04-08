
SUFFIXES = .o .f90

.f90.o:
	@FC@ @FCFLAGS@ @FCFLAGS_free_f90@ ${DEFAULT_INCLUDES} @FC_MOD_FLAG@${top_builddir}/fortran -c $<

noinst_PROGRAMS	= axtest

axtest_SOURCES 	= axtest.f90 
axtest_LDADD	= $(top_builddir)/fortran/libfagentx.la $(top_builddir)/src/libagentx.la
