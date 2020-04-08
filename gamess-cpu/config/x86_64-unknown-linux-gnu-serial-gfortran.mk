#doc  This is a the machine dependant file for the serial build of GAMESS-UK
#doc  on 64-bit Linux with the gfortran compiler
#doc
#doc               *  *  *     B E W A R E   *  *  *
#doc
#doc  This is not a supported build - the file is provided for developers
#doc  and others who which to experiment with this compiler
#doc  Even at the default optimisation level (-O0) some functionality
#doc  will not work at present.
#doc 
#doc NB: dl-find is not supported with gfortran < 4.2 as it uses allocatable
#doc     arrays in type declarations which is an extension to the F90 standard:
#doc     http://coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2004-12/0452.html
#doc     newscf is also not supported as it uses complex character array constructors
#doc 
#doc  Options:
#doc  static - create a statically linked binary
##
# Default options:
#dopt vb mrdci zora drf nbo vdw sysmo O0 static 
#opt i8  mopac
#
# ================ M4 Processing options

#
#
#  acml   - link against PGI's ACML blas library 
#           (default location:  ~/acml2.7.0/pathscale64/lib/libacml.a)

#
MACHINE_KEY=G
MACHOPT=opteron,linux,littleendian,cio,unix,upck-equiv,GFS,glibc

# Bitwise operators that differ from the defaults.
IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)
# 
#
RANLIB=ar -s
CPP=/lib/cpp
GA_F77_DEFS = -traditional
#
# ===============  Compiler Options
#
#ensure 64bit compilers invoked
FC = gfortran
FC90 = gfortran
LD = gfortran
CC = gcc

#Default compilation Flags
#--#if i8#
FFLAGSI8 = -fdefault-integer-8
CFLAGSI8 = -DLONG_INTEGER
I8_M4_OPT = ,i8
#--#else#
I8_M4_OPT = ,i8drct,64bitpointers
#--#endif i8#

#--#if coverage#
FLAGSCOV  = -fprofile-arcs -ftest-coverage
#--#endif#

#--#if profile#
FFLAGSPG  = -pg
CFLAGSPG  = -pg
LDFLAGSPG = -pg
#--#endif#

FFLAGSD = -c -fno-second-underscore ${FFLAGSI8} -fd-lines-as-comments ${FFLAGSPG} ${FLAGSCOV}
CFLAGSD = -c ${CFLAGSI8} ${CFLAGSPG} ${FLAGSCOV}

#--#if debug#
FFLAGSV  = ${FFLAGSD} -g -ffpe-trap=invalid,zero,overflow
FFLAGSS  = ${FFLAGSD} -g -ffpe-trap=invalid,zero,overflow
FFLAGSN  = ${FFLAGSD} -g -ffpe-trap=invalid,zero,overflow
FFLAGSN0 = ${FFLAGSD} -g -ffpe-trap=invalid,zero,overflow
CFLAGS   = ${CFLAGSD} -g
LDTMP    = -g
#--#elseif O0#
FFLAGSV  = ${FFLAGSD} -O0 
FFLAGSS  = ${FFLAGSD} -O0 
FFLAGSN  = ${FFLAGSD} -O0 
FFLAGSN0 = ${FFLAGSD} -O0
CFLAGS   = ${CFLAGSD} 
LDTMP    = -g
#--#else# 
FFLAGSV  = ${FFLAGSD}  -O2 
FFLAGSS  = ${FFLAGSD}  -O
FFLAGSN  = ${FFLAGSD}  -O1
FFLAGSN0 = ${FFLAGSD} 
CFLAGS   = -D_REENTRANT ${CFLAGSD}
LDTMP    = 
#--#endif debug#
#
#--#if static#
LDFLAGS=${LDTMP} ${LDFLAGSPG} ${FLAGSCOV} -static
#--#else#
LDFLAGS=${LDTMP} ${LDFLAGSPG} ${FLAGSCOV}
#--#endif static#
#
#
#--#if acml#
BLASOPT=,blas
#acml needs to be linked statically for it to work.
LIBBLAS= ~/acml2.7.0/pathscale64/lib/libacml.a
#--#endif acml#
#
OPTIONS=${MACHOPT}${I8_M4_OPT}${BLASOPT}
BL_LIB=${LIBBLAS}
#

# ===============  Additional Files 
#
EXTRA_MP2= check0a.o
EXTRA_DFT = jkint_dft.o jkder_dft.o
#
# ===============  Compiler Exceptions 
#
rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

eispack.o:	eispack.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

check0a.o:	newmrd5.m
	cat ../utilities/gener.m newmrd5.m | $(M4) -DGEN_EXTRACTFILE=check0a $(M4OPTS) > check0a.f
	$(FC) $(FFLAGSN) check0a.f
#
#  ========== DFT Exceptions (PGF) ===============
#
jkint_dft.o:	integ2e.m
	cat ../utilities/gener.m  integ2e.m | $(M4) -DGEN_EXTRACTFILE=jkint_dft $(M4OPTS) > jkint_dft.f
	$(FC) $(FFLAGSN) jkint_dft.f
#
jkder_dft.o:	deriv2e.m
	cat ../utilities/gener.m  deriv2e.m | $(M4) -DGEN_EXTRACTFILE=jkder_dft $(M4OPTS) > jkder_dft.f
	$(FC) $(FFLAGSN) jkder_dft.f

