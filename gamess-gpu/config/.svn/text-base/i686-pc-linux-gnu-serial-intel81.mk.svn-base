#doc  This is the machine-specific file for the serial based version of GAMESS-UK
#doc  using Intel compiler 8.1 running on 2GHz Xeon CPU's
#doc
#doc  Options:
#doc  mkl:     build against Intel's mkl library (expected in /opt/intel/mkl701)
#doc  goto:    build against goto blas library (expected in /usr/lib)
#doc  static:  make a static binary
#doc  gcc:     use GNU GCC compiler if Intel C-compiler unavailable (vdw may not work with gcc<4.0)
#doc  i8: Use 64 bit integers as default (Don't use - is bust).
#
#dopt mrdci zora vb drf nbo vdw sysmo dl-find
#opt  mkl goto static gcc i8 mopac
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,upck-equiv

#
# ===============  Compiler Options
#
## INTEL COMPILER VERSION 8.1 ##
FC = ifort
LD = ifort
FC90 = ifort
#--#if gcc#
CC=  gcc
CXX = g++
#--#else#
CC=  icc
CXX = icc
#--#endif#

RANLIB=ar -s


#--#if profile#
FFLAGSPROFILE  = -pg
CFLAGSPROFILE  = -pg
LDFLAGSPROFILE = -pg
#--#endif profile#

#--#if i8#
FFLAGSI8=-i8
CFLAGSI8= -DLONG_INTEGER
I8_M4_OPT=,i8
#--#else#
FFLAGSI8=
CFLAGSI8=
I8_M4_OPT=,i8drct
#--#endif i8#


FFLAGSTMP  = -c ${FFLAGSPROFILE} ${FFLAGSI8}
CFLAGSTMP  = -c ${CFLAGSPROFILE} ${CFLAGSI8}
LDFLAGSTMP =    ${LDFLAGSPROFILE} 



#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g -fpe0
FFLAGSS = ${FFLAGSTMP} -g -fpe0
FFLAGSN = ${FFLAGSTMP} -g -fpe0
FFLAGSN0 = ${FFLAGSTMP} -g -O0 -fpe0
CFLAGS = ${CFLAGSTMP} -g -c -DLINUX -D_FILE_OFFSET_BITS=64 
CXXFLAGS = ${CFLAGSTMP} -g -c
#--#else# 
FFLAGSV = ${FFLAGSTMP} -ip -O3
FFLAGSS = ${FFLAGSTMP} -ip -O2 
FFLAGSN = ${FFLAGSTMP} -ip -O1
FFLAGSN0 = ${FFLAGSTMP} -O0
CFLAGS  = ${CFLAGSTMP} -DLINUX -D_FILE_OFFSET_BITS=64
CXXFLAGS = ${CFLAGSTMP} -O2
#--#endif#

#--#if static#
#---#if gcc#
LDFLAGS  = ${LDFLAGSTMP} -Wl,-Map load.map -D_FILE_OFFSET_BITS=64 -static
#---#else#
LDFLAGS  = ${LDFLAGSTMP} -Wl,-Map load.map -D_FILE_OFFSET_BITS=64 -static -static-libcxa
#---#endif gcc#
#--#else#
LDFLAGS  = ${LDFLAGSTMP} -Wl,-Map load.map -D_FILE_OFFSET_BITS=64
#--#endif#

#--#if mkl#
#---#if static#
lBLAS= /opt/intel/mkl701/lib/32/
LIBBLAS= ${IBLAS}libmkl_lapack.a ${IBLAS}libmkl_ia32.a ${IBLAS}libguide.a -lpthread
# NOTE (jmht): mfg needed the -i-static flag to get this to build but it seems to be redundant.
# ALSO using mkl 7.01, the static binary runs fine under bash but seg faults under csh
#---#else#
LIBBLAS= -L/opt/intel/mkl701/lib/32  -lguide -lmkl -lmkl_lapack64
#---#endif static#
BLASOPT=,blas
#--#elseif goto#
LIBBLAS=-L/usr/lib -lgotoblas
BLASOPT=,blas
#--#endif#

# Bring all the options together
OPTIONS=${MACHOPT}${BLASOPT}${I8_M4_OPT}
BL_LIB = ${LIBBLAS}

# Diesel build options
LD_DIESEL = ${CXX}
#--#if mkl#
DIESEL_LIBS = ${LIBBLAS} -lifcore
#--#elseif goto#
DIESEL_LIBS = ${LIBBLAS} -lifcore
#--#else#
DIESEL_LIBS = -lifcore
#--#endif#
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

mrdci1.o:	mrdci1.f
	$(FC) $(FFLAGSN0) $*.f

dircta.o:	dircta.f
	$(FC) $(FFLAGSN0) $*.f

#
# Information on the machine this build was configured on
# Curtailed output of /proc/cpuinfo
# vendor_id       : GenuineIntel
# cpu family      : 15
# model           : 2
# model name      : Intel(R) XEON(TM) CPU 2.00GHz
# stepping        : 4
# cpu MHz         : 1999.816
# cache size      : 512 KB
# 
# Output of uname -a
# Linux enterprise 2.4.20-28.7smp #1 SMP Thu Dec 18 11:18:31 EST 2003 i686 unknown
# 
# Using Intel Compilers, version 8.1
# Output of: /opt/intel_fc_80/bin/ifc -V
# Intel(R) Fortran Compiler for 32-bit applications, Version 8.1    Build 20041019Z Package ID: l_fc_pu_8.1.021
# 
# Output of: /opt/intel_cc_80/bin/icc -V
# Intel(R) C++ Compiler for 32-bit applications, Version 8.1    Build 20041019Z Package ID: l_cc_pu_8.1.024
# 
