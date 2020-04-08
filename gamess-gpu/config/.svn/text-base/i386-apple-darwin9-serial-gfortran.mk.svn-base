#doc  This is a the machine dependant file for the serial build of GAMESS-UK
#doc  on Mac OSX on Intel i686
#doc
#doc  This is an almost supported build - the file is provided for developers
#doc  and others who which to experiment with this compiler
#doc  gfortran > 4.2 is recommended
#doc 
#doc  Options:
#doc     universal - build a fat binary that supports both ppc & i686 architectures
#doc     i8 - build Fortran code with 64-bit integers. NOTE!!! This has only been tested
#doc          with gfortran 4.4.0, and was broken as the addresses of integers passed
#doc          between C and Fortran was getting mangled - not sure how to fix yet.
#doc
#doc     for static mv all dylib files out of the way (in /usr/local/lib) - a bit nasty
#
# DEFAULT AND OPTIONAL OPTIONS
#dopt zora vb drf mrdci nbo vdw mopac sysmo dl-find newscf O2 
#opt demo debug O0 O1 universal i8
#
# ================ M4 Processing options
#
MACHINE_KEY=G
# Need macosx for getlog in machscf.m & printing banner
MACHOPT=linux,pclinux,macosx,littleendian,cio,unix,upck-equiv,GFS

IAND32=iand($$1,$$2)
IOR32=ior($$1,$$2)
IXOR32=ieor($$1,$$2)
ISHFT32=ishft($$1,$$2)

FC = gfortran
LD = gfortran
FC90 = gfortran
CC = gcc
CXX = g++
FFLAGSTMP = -c
RANLIB = ranlib
GMAKE=make
#
# ================ Default compilation Flags
#
#--#if i8#
FFLAGSI8 = -fdefault-integer-8
CFLAGSI8 = -DLONG_LONG_INTEGER
I8_M4_OPT=,i8
#--#else#
FFLAGSI8 =
CFLAGSI8 =
I8_M4_OPT=,i8drct
#--#endif i8#

#--#if universal#
FFLAGSUNI  = -arch i386 -arch ppc
CFLAGSUNI  = -arch i386 -arch ppc
LDFLAGSUNI = -arch i386 -arch ppc
#--#endif#

FFLAGSTMP = -c -DGFS -fno-second-underscore $(FFLAGSI8) $(FFLAGSUNI) -fd-lines-as-comments
CFLAGSTMP = -c -DLINUX -D_FILE_OFFSET_BITS=64 -DGFS $(CFLAGSI8) ${CFLAGSUNI}
LDFLAGSTMP= $(LDFLAGSUNI)

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g -DCDBG="    "
FFLAGSS = ${FFLAGSTMP} -g -DCDBG="    "
FFLAGSN = ${FFLAGSTMP} -g -DCDBG="    "
FFLAGSN0 = ${FFLAGSTMP} -g -DCDBG="    " 
CFLAGS = ${CFLAGSTMP} -g 
LDFLAGS = ${LDFLAGSTMP} -g
#--#elseif O0#
FFLAGSV = ${FFLAGSTMP}  -O0
FFLAGSS = ${FFLAGSTMP}  -O0
FFLAGSN = ${FFLAGSTMP}  -O0
FFLAGSN0 = ${FFLAGSTMP} -O0
CFLAGS = ${CFLAGSTMP}
LDFLAGSTMP =  ${LDFLAGSTMP} -O0
#--#elseif O1#
FFLAGSV = ${FFLAGSTMP}  -O1
FFLAGSS = ${FFLAGSTMP}  -O1
FFLAGSN = ${FFLAGSTMP}  -O1
FFLAGSN0 = ${FFLAGSTMP} -O0
CFLAGS = ${CFLAGSTMP} 
LDFLAGSTMP =  ${LDFLAGSTMP} -O1
#--#else# 
#### Beware: gfortran optimization level>1 is not optimal
#
# Need to disable -malign-double as it causes segfaults with the
# latest version of gfortran. It seems that this is because the
# gfortran libraries themselves have not been compiled with the
# -malign-double flag
# However, this change causes chap2/c2015_e to give an error
# indicating that the casscf code may be broken in this version
FFLAGSV = ${FFLAGSTMP} -O2
FFLAGSS = ${FFLAGSTMP} -O1
FFLAGSN = ${FFLAGSTMP} -O1
FFLAGSN0 = ${FFLAGSTMP} -O0
CFLAGS = ${CFLAGSTMP}
CXXFLAGS = -c -O1
LDFLAGS = ${LDFLAGSTMP} -O2
#--#endif#

## M4 Options ##
OPTIONS=$(MACHOPT)$(I8_M4_OPT)$(BLASOPT)

CHMHOST=gnu

# Diesel build options 
LD_DIESEL = ${CXX}
DIESEL_LIBS = -lg2c
FLEX = /usr/bin/flex
__UNDERBAR = 1
SIZEOF_VOID_P = 4
LONG_LONG_INT = long long int
LONG_INT = long int
INT = int
SHORT_INT = short int


#
# ===============  Additional Files 
#
EXTRA=gethes.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o

#
# ===============  Compiler Exceptions
#

# Exceptions for the problem also dealt with by the GFS m4 macro (see Marc's comments
# on BOZ constants in (e.g.) machscf.m
mainci.o:	mainci.m
	cat  ../utilities/gener.m mainci.m | $(M4) $(M4OPTS) > mainci.f
	$(FC) $(FFLAGSV) -fno-range-check $(OBJNAME) $*.f
util5.o:	util5.m
	cat  ../utilities/gener.m util5.m | $(M4) $(M4OPTS) > util5.f
	$(FC) $(FFLAGSV) -fno-range-check $(OBJNAME) $*.f
# End GFS files

gethes.o:	casb.m
	cat  ../utilities/gener.m casb.m | $(M4) -DGEN_EXTRACTFILE=gethes $(M4OPTS) > gethes.f
	$(FC) $(FFLAGSN) gethes.f

mkmakw.o:	drvmp.m
	cat  ../utilities/gener.m drvmp.m | $(M4) -DGEN_EXTRACTFILE=mkmakw $(M4OPTS) > mkmakw.f
	$(FC) $(FFLAGSS) mkmakw.f

mpmakw.o:	secmp2.m
	cat  ../utilities/gener.m secmp2.m | $(M4) -DGEN_EXTRACTFILE=mpmakw $(M4OPTS) > mpmakw.f
	$(FC) $(FFLAGSS) mpmakw.f

umpe3a.o:	mp3.m
	cat  ../utilities/gener.m mp3.m | $(M4) -DGEN_EXTRACTFILE=umpe3a $(M4OPTS) > umpe3a.f
	$(FC) $(FFLAGSS) umpe3a.f

umpe3b.o:	mp3.m
	cat  ../utilities/gener.m mp3.m | $(M4) -DGEN_EXTRACTFILE=umpe3b $(M4OPTS) > umpe3b.f
	$(FC) $(FFLAGSS) umpe3b.f

check0a.o:	newmrd5.m
	cat ../utilities/gener.m newmrd5.m | $(M4) -DGEN_EXTRACTFILE=check0a $(M4OPTS) > check0a.f
	$(FC) $(FFLAGSN) check0a.f

tst1s.o:	direct.m
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f 

dircta.o:	dircta.f
	$(FC) $(FFLAGSN0) $*.f

rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

mrdci1.o:	mrdci1.f
	$(FC) $(FFLAGSN0) $*.f
