#doc  This is a the machine dependant file for the serial build of GAMESS-UK
#doc  on 64-bit Linux with the gfortran versions > 4.2
#doc
#doc 
#doc  Options:
#doc  static - create a statically linked binary
#doc  i8 - build with 64bit FORTRAN integers (integer*8)
#doc  xml - enable module to read/write XML(CML) data files.
##
# Default options:
#dopt vb mrdci zora drf nbo vdw sysmo dl-find newscf 
#opt i8 static 00 xml mopac
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=opteron,linux,pclinux,littleendian,cio,unix,upck-equiv,GFS,glibc

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

FFLAGSD = -fno-second-underscore ${FFLAGSI8} -fd-lines-as-comments
CFLAGSD =  ${CFLAGSI8}

#--#if debug#
FFLAGSV = -c ${FFLAGSD} -g 
FFLAGSS = -c ${FFLAGSD} -g 
FFLAGSN = -c ${FFLAGSD} -g 
FFLAGSN0 = -c ${FFLAGSD} -g
CFLAGS =  -c ${CFLAGSD} -g
LDTMP  = -g
#--#elseif O0#
FFLAGSV = -c ${FFLAGSD} -O0 
FFLAGSS = -c ${FFLAGSD} -O0 
FFLAGSN = -c ${FFLAGSD} -O0 
FFLAGSN0 = -c ${FFLAGSD} -O0
CFLAGS =  -c ${CFLAGSD} 
LDTMP  = -g
#--#else# 
FFLAGSV = -c ${FFLAGSD}  -O2 
FFLAGSS = -c ${FFLAGSD}  -O
FFLAGSN = -c ${FFLAGSD}  -O1
FFLAGSN0 = -c ${FFLAGSD} 
CFLAGS  = -c -D_REENTRANT ${CFLAGSD}
LDTMP = 
#--#endif debug#
#
#--#if static#
LDFLAGS=${LDTMP} -static
#--#else#
LDFLAGS=${LDTMP}
#--#endif static#
#
OPTIONS=${MACHOPT}${I8_M4_OPT}${BLASOPT}
BL_LIB=${LIBBLAS}
#

# ===============  Additional Files 
#
EXTRA=gethes.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o
EXTRA_DFT = jkint_dft.o jkder_dft.o
#
# ===============  Compiler Exceptions 
#
# Exceptions for the problem also dealt with by the GFS m4 macro (see Marc's comments
# on BOZ constants in (e.g.) machscf.m
mainci.o:       mainci.m
	cat  ../utilities/gener.m mainci.m | $(M4) $(M4OPTS) > mainci.f
	$(FC) $(FFLAGSV) -fno-range-check $(OBJNAME) $*.f
util5.o:        util5.m
	cat  ../utilities/gener.m util5.m | $(M4) $(M4OPTS) > util5.f
	$(FC) $(FFLAGSV) -fno-range-check $(OBJNAME) $*.f
# End GFS files
anala.o:	anala.f
	$(FC) $(FFLAGSN) $*.f

rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

eispack.o:	eispack.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

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

