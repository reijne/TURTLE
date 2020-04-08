#
# $Id: GNUmakefile,v 1.2 2000-10-26 15:38:00 psh Exp $
#
#include ../../config/makefile.h
include ../DEFS

MXOBJ	=	mxbrod.o \
		mxcombv1.o \
		mxsubs.o \
		sum.o \
		maxx.o \
		peigs_dgop.o

all:		$(MXOBJ)

.EXPORT_ALL_VARIABLES:

.PHONY :	clean

#
# Note that under Cygwin the rm commands are case-insensitive, so the
# commands below intended to remove intermediate .f files actually remove
# the original .F files!  Since intermediate files are not generated for
# Cygwin, this is avoided by commenting out the rm statements.
#
clean:		
		-rm -f *.o
ifneq ($(peigs_TARGET),$(findstring $(peigs_TARGET),CYGNUS X86OSX))
		-rm -f mxbrod.f
		-rm -f mxcombv1.f
		-rm -f mxsubs.f
		-rm -f sum.f peigs_dgop.f
endif

realclean:
		clean

.SUFFIXES:	.o .s .F .f .c

.c.o:;	$(peigs_CC) $(peigs_COPT) -c $<
.f.o:;	$(peigs_FC) $(peigs_FOPT) -c $<
.F.o:;	$(peigs_FC) $(peigs_FOPT) -c $<


