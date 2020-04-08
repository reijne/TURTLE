## Process this file with automake to produce Makefile.in


INCLUDES = -I$(top_builddir)/include -I@srcdir@/include @THREAD_CFLAGS@ @Z_CFLAGS@

lib_LTLIBRARIES = libaxxml2.la
libaxxml2_la_LIBADD = @THREAD_LIBS@ @Z_LIBS@ $(ICONV_LIBS) @M_LIBS@ @WIN32_EXTRA_LIBADD@

libaxxml2_la_LDFLAGS = @CYGWIN_EXTRA_LDFLAGS@ @WIN32_EXTRA_LDFLAGS@ -version-info @LIBXML_VERSION_INFO@ @MODULE_PLATFORM_LIBS@

if WITH_TRIO_SOURCES
libaxxml2_la_SOURCES = entities.c encoding.c error.c parserInternals.c  \
		parser.c hash.c list.c xmlIO.c xmlmemory.c uri.c  \
		xpath.c  \
		xpointer.c nanohttp.c \
		globals.c xmlstring.c \
		xmlunicode.c \
		triostr.c trio.c dict.c SAX2.c \
		chvalid.c tree.c threads.c SAX.c xmlmodule.c xmlregexp.c legacy.c valid.c
else
libaxxml2_la_SOURCES  = entities.c encoding.c error.c parserInternals.c  \
		parser.c hash.c list.c xmlIO.c xmlmemory.c uri.c  \
		xpath.c  \
		xpointer.c nanohttp.c \
		globals.c xmlstring.c \
		xmlunicode.c \
		dict.c SAX2.c \
		chvalid.c tree.c threads.c SAX.c xmlmodule.c xmlregexp.c legacy.c valid.c

endif

maintainer-clean-local:
			@rm -f *~
			@rm -f */*~
			@rm -f */*/*~
			@rm -f */*/*/*~
			@rm -f #*
			@rm -f */#*
			@rm -f */*/#*
			@rm -f */*/*/#*

