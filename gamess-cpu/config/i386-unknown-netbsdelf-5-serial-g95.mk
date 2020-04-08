#doc Machine-dependent file for GAMESS-UK serial build on Linux/Pentium
#doc with g95 compiler G95 (GCC 4.0.3 (g95 0.90!) Jul 27 2006)
#doc NB: Mopac is not available with this build as Real and double precision
#doc DO loop index variables are not implemented in g95.
#
# DEFAULT AND OPTIONAL OPTIONS
#dopt zora vb drf mrdci vdw nbo sysmo dl-find 
#opt demo debug static mopac
#
# ================ M4 Processing options
#
MACHINE_KEY=G
#jmht - removed bits8 & replaced with i8drct
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,i8drct

IAND32=iand($$1,$$2)
IOR32=ior($$1,$$2)
IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)

#jmht ISHFT32 doesn't appear to be needed
#ISHFT32=ishft($$1,$$2)

FC = g95
LD = g95
FC90 = g95
CC = gcc
CXX = g++
# -Wno-globals to not complain about inconsistent types between subroutines
FFLAGSTMP = -c -fsloppy-char -Wno-globals
RANLIB = ranlib
BL_LIB=-lcompat
#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g -DCDBG="    "
FFLAGSS = ${FFLAGSTMP} -g -DCDBG="    "
FFLAGSN = ${FFLAGSTMP} -g -DCDBG="    "
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS = -g -c -DLINUX -D_FILE_OFFSET_BITS=64 -DLINUXF2C
CXXFLAGS = -g -c 
LDFLAGSTMP = -D_FILE_OFFSET_BITS=64 -g
#--#else# 
#FFLAGSV = ${FFLAGSTMP} -O3 -malign-double
#FFLAGSS = ${FFLAGSTMP} -O -malign-double
#FFLAGSN = ${FFLAGSTMP} -O1 -malign-double
FFLAGSV = ${FFLAGSTMP} 
FFLAGSS = ${FFLAGSTMP} 
FFLAGSN = ${FFLAGSTMP} 
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS  = -c -DLINUX -D_FILE_OFFSET_BITS=64 -DLINUXF2C
CXXFLAGS = -c -O2
LDFLAGSTMP = -D_FILE_OFFSET_BITS=64
#--#endif#
# Need -fno-globals for drf/drfsub.f or else g77 crashes without warning

#--#if static#
LDFLAGS  = ${LDFLAGSTMP} -static -static-libgcc 
#--#else#
LDFLAGS  = ${LDFLAGSTMP}
#--#endif#

OPTIONS=${MACHOPT}

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

# Hack to add -DMAF_INTERNAL to compilation flags for glue_uk_us.o
# The file mafdecls.fh is included in glue_uk_us.m. This file declares
# constants for the different data types as well as the function names
# for the MA library. However, glue_uk_us.m  NEVER actually calls any
# of the MA functions. A g95 quirk is that if a function is declared
# the corresponding symbol must exist at link time even if it is never 
# used. As a result the link will fail on an irrelevant technically. 
# In mafdecls.fh there is a way to suppress the function type definitions
# for when the file is included in the MA source code itself by specifying 
# the -DMAF_INTERNAL on the preprocessor line. We hijack GA_F77_DEFS for this.
GA_F77_DEFS=-DMAF_INTERNAL

#
# ===============  Additional Files 
#
EXTRA=gethes.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o

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

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f
#
# Information on the machine this build was configured on
#
# Output of uname -a :
# Linux csevig6 2.4.24_usb #1 SMP Mon Feb 16 17:03:30 GMT 2004 i686 i686 i386 GNU/Linux
# 
# Output of cat /proc/cpuinfo :
# 
# processor	: 0
# vendor_id	: GenuineIntel
# cpu family	: 15
# model		: 2
# model name	: Intel(R) Pentium(R) 4 CPU 2.80GHz
# stepping	: 9
# cpu MHz		: 2793.036
# cache size	: 512 KB
# fdiv_bug	: no
# hlt_bug		: no
# f00f_bug	: no
# coma_bug	: no
# fpu		: yes
# fpu_exception	: yes
# cpuid level	: 2
# wp		: yes
# flags		: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe cid
# bogomips	: 5570.56
