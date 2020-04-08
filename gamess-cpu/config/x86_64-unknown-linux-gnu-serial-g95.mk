#doc This is the machine-specific file for the serial version of GAMESS-UK
#doc on a Xeon E5320 or Opteron with the g95 compiler
#doc NB: nbo should be fixed
#doc NOTE! - this build _sort of_ works. A number of
#doc         examples fail due to what appear to be
#doc         bugs in g95, so please use with care.
#doc
#doc Options:
#doc datain - force GAMESS-UK to read it's input from 
#doc          a file called datain instead of standard input
#doc mopac  - include mopac code in the build
#doc debug  - include debugging information in the objects (no optimsation)
#doc static - build a static binary
#
# Default options:
#dopt zora vb drf mrdci ma nbo vdw i8 
#opt debug static datain
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=opteron,linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,glibc

# Bitwise operators that differ from the defaults.
IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)
IAND32=iand($$1,$$2)
IOR32=iand($$1,$$2)
IAND64=iand($$1,$$2)
IOR64=iand($$1,$$2)
#
#
RANLIB=ranlib
CPP=/lib/cpp
#
# ===============  Compiler Options
#

FC = g95
FC90 = g95
LD = g95
CC = gcc

#Default compilation Flags
#--#if ma#
GA_F77_DEFS=-traditional
GA_TARGET=LINUX64
#MA_OPTIONS=F2C_TWO_UNDERSCORES=1
#--#endif ma#

#--#if i8#
FFLAGSI8 = -i8
CFLAGSI8 = -DLONG_INTEGER
I8_M4_OPT=,i8
#--#else#
FFLAGSI8 = 
CFLAGSI8 = -DSTD_INT
I8_M4_OPT=,64bitpointers
#--#endif#

# Default flags for all compilations (including debug)
FFLAGSD = -c  -fsloppy-char -fno-second-underscore ${FFLAGSI8}
CFLAGSD= -c  -DLINUX -D_FILE_OFFSET_BITS=64 ${CFLAGSI8}

#--#if static#
LDFLAGSD= -D_FILE_OFFSET_BITS=64 -static -static-libgcc 
#--#else#
LDFLAGSD= -D_FILE_OFFSET_BITS=64
#--#endif static#
#

#--#if debug#
FFLAGSV  = ${FFLAGSD}  -g -DCDBG="    "
FFLAGSS  = ${FFLAGSD}  -g -DCDBG="    "
FFLAGSN  = ${FFLAGSD}  -g -DCDBG="    "
FFLAGSN0 = ${FFLAGSD}  -g
CFLAGS   = ${CFLAGSD}  -g
LDFLAGS  = ${LDFLAGSD} -g 
#--#else#
FFLAGSV  =  ${FFLAGSD}  -O2
FFLAGSS  =  ${FFLAGSD}  -O
FFLAGSN  =  ${FFLAGSD}  -O1
FFLAGSN0 =  ${FFLAGSD}
CFLAGS   =  ${CFLAGSD}
LDFLAGS  =  ${LDFLAGSD}
#--#endif debug#
#
#

# Bring all the options together
OPTIONS=${MACHOPT}${BLASOPT}${I8_M4_OPT}

# Add any extra libraries here
BL_LIB =

#
# ===============  Additional Files 
#
EXTRA=gethes.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o
EXTRA_DFT = jkint_dft.o jkder_dft.o
#
# ===============  Compiler Exceptions
#
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
	cat ../machines/$(MACH) ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f 

rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

eispack.o:	eispack.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

#wrong results for c2029_a
tsort.o:	tsort.f
	$(FC) $(FFLAGSN) $*.f
#
#  ========== DFT Exceptions (PGF) ===============
#
jkint_dft.o:	integ2e.m
	cat ../machines/$(MACH) ../utilities/gener.m  integ2e.m | $(M4) -DGEN_EXTRACTFILE=jkint_dft $(M4OPTS) $(SNGL) > jkint_dft.f
	$(FC) $(FFLAGSN) jkint_dft.f
#
jkder_dft.o:	deriv2e.m
	cat ../machines/$(MACH) ../utilities/gener.m  deriv2e.m | $(M4) -DGEN_EXTRACTFILE=jkder_dft $(M4OPTS) $(SNGL) > jkder_dft.f
	$(FC) $(FFLAGSN) jkder_dft.f
