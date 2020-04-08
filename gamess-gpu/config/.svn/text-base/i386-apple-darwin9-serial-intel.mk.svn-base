#doc  This is the machine-specific file for the serial based version of GAMESS-UK
#doc  Working example using Intel compilers 11.1.089 on a i686 platform
#doc  Using Xcode 3.2.3
#doc  Should work with Intel 10.x and 11.x Compilers and MKL 10.x
#doc
#doc  Has been tested on Mac OS X 10.6.4. Chap2 and VB checkes out OK
#doc  
#doc  Compiler/Library/Include paths are assumed and set by these lines 
#doc   source /opt/intel/Compiler/<current version>/bin/intel64/iccvars_intel64.csh
#doc   source /opt/intel/Compiler/<current version>/bin/intel64/ifortvars_intel64.csh
#doc   source /opt/intel/Compiler/<current version>/Frameworks/mkl/tools/environment/mklvarsem64t.csh
#doc  Set these sources in Your .cshrc file.
#doc
#doc  Options:
#doc  default build assumes 64bit and i8
#
#dopt mrdci zora vb drf nbo vdw mkl i8
#opt  debug mopac sysmo dl-find
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,macosx,pclinux,littleendian,cio,unix,upck-equiv

IAND32=iand($$1,$$2)
IOR32=ior($$1,$$2)
IXOR32=ieor($$1,$$2)
ISHFT32=ishft($$1,$$2)
IAND64=iand($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)

MKLDIR=${MKLROOT}/lib/em64t

#
# ===============  Compiler Options
#
FC = ifort
FC90 = ifort
CC= icc
CXX = icc
LD = ifort

RANLIB=ranlib

#--#if profile#
FFLAGSPROFILE = -pg 
LDFLAGSPROFILE =  -pg
#--#else #
FFLAGSPROFILE =
LDFLAGSPROFILE =
#--#endif#

#--#if static#
LDFLAGSSTATIC  = -static -static-libcxa
#--#endif#

#--#if i8#
FFLAGSI8 = -i8
CFLAGSI8 = -DLONG_INTEGER
I8_M4_OPT=,i8
#--#else#
FFLAGSI8=
CFLAGSI8=
I8_M4_OPT=,64bitpointers
#--#endif i8#

# Set the default flags
FFLAGSTMP= -c -I${FPATH} $(FFLAGSPROFILE) ${FFLAGSI8}
CFLAGSTMP= -c -no-multibyte-chars -DLINUX -D_FILE_OFFSET_BITS=64 -I${CPATH} ${CFLAGSI8}
CXXFLAGSTMP= $(CFLAGSTMP)
LDFLAGSTMP= $(LDFLAGSSTATIC) $(LDFLAGSPROFILE) 

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g 
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g -O0 
CFLAGS = ${CFLAGSTMP} -g 
CXXFLAGS = ${CXXFLAGSTMP} 
LDFLAGS =  -g $(LDFLAGSTMP)
#--#else# 
FFLAGSV = ${FFLAGSTMP} -ip -O1
FFLAGSS = ${FFLAGSTMP} -ip -O1 
FFLAGSN = ${FFLAGSTMP} -ip -O1
FFLAGSN0 = ${FFLAGSTMP} -O0
CFLAGS  = ${CFLAGSTMP}
CXXFLAGS = -g ${CXXFLAGSTMP}
LDFLAGS = $(LDFLAGSTMP)
#--#endif#


#--#if mkl#
#---#if static#
lBLAS= -L${MKLDIR}
LIBBLAS= ${IBLAS}libmkl_intel_ilp64.a ${IBLAS}libmkl_intel_sp2dp.a ${IBLAS}libmkl_core.a ${IBLAS}libmkl_mc.a
#---#else#
LIBBLAS= -lmkl_intel_ilp64 -lmkl_intel_sp2dp -lmkl_core -lmkl_mc
#---#endif static#
BLASOPT=,blas
#--#elseif goto#
LIBBLAS=-L/usr/lib -lgotoblas
BLASOPT=,blas
#--#else#
LIBBLAS=
BLASOPT=
#--#endif#

# Bring all the options together
OPTIONS=${MACHOPT}${I8_M4_OPT}${BLASOPT}
BL_LIB = ${LIBBLAS}

# for CHARMM, determines location of library
CHMHOST=gnu

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
# Output of /usr/sbin/system_profiler SPHardwareDataType 
#  Hardware Overview:
#
#  Model Name: Mac Pro
#  Model Identifier: MacPro3,1
#  Processor Name: Quad-Core Intel Xeon
#  Processor Speed: 2,8 GHz
#  Number Of Processors: 1
#  Total Number Of Cores: 4
#  L2 Cache: 12 MB
#  Memory: 8 GB
#  Bus Speed: 1,6 GHz
#  Boot ROM Version: MP31.006C.B05
#  SMC Version (system): 1.25f4
#  Serial Number (system): CK8100MSXYL
#  Hardware UUID: 75AD3FD6-6641-5CB4-AD4B-3C4A55607B8A
# 
# Output of uname -a (it is 10.6.4 Mac OS X)
#Darwin mac100.chem.uu.nl 10.4.0 Darwin Kernel Version 10.4.0: Fri Apr 23 18:28:53 PDT 2010; root:xnu-1504.7.4~1/RELEASE_I386 i386
# 
# Using Intel Compilers, version 11.1
# Output of: /opt/intel/Compiler/11.1/089/bin/intel64/ifort -V
# Intel(R) Fortran Intel(R)64 Compiler Professional for applications running on Intel(R)64, Version 11.1 Build 20100806 Package ID: m_cprof_p_11.1.089
# Copyright (C) 1985-2010 Intel Corporation.  All rights reserved.
#
# Output of: /opt/intel/Compiler/11.1/089/bin/intel64/icc -V
# Intel(R) C Intel(R) 64 Compiler Professional for applications running on Intel(R) 64, Version 11.1    Build 20100806 Package ID: m_cproc_p_11.1.089
# Copyright (C) 1985-2010 Intel Corporation.  All rights reserved.
# 
