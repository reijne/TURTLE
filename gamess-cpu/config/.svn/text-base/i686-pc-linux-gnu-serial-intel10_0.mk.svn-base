#doc  This is the machine-specific file for the serial based version of GAMESS-UK
#doc  Test version using Intel compilers 10.0.023 on a i686 platform
#doc
#doc  Trial version only; the default build is 'debug' for a reason
#doc  Beware: This version becomes untested after 27-10-2007 (when my trial license expires)
#doc  
#doc  Compiler/Library/Include paths are hardcoded for the Intel demo install tree
#doc
#doc  Options:
#doc  
#
#dopt mrdci zora vb drf nbo vdw sysmo dl-find 
#opt  mopac debug
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,upck-equiv,bits8


IAND32=iand($$1,$$2)
IOR32=ior($$1,$$2)
IXOR32=ieor($$1,$$2)
ISHFT32=ishft($$1,$$2)
IAND64=iand($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)

#
# ===============  Compiler Options
#
## INTEL COMPILER VERSION 10.0.023
FC = /opt/intel/fc/10.0.023/bin/ifort
LD = /opt/intel/fc/10.0.023/bin/ifort
FC90 = /opt/intel/fc/10.0.023/bin/ifort
CC=  /opt/intel/cc/10.0.023/bin/icc
CXX = /opt/intel/cc/10.0.023/bin/icc

LDFLAGSTMPL =  -L /opt/intel/cc/10.0.023/lib -L/opt/intel/fc/10.0.023/lib
#--#if profile#
FFLAGSTMP = -c -pg -I/opt/intel/fc/10.0.023/include
CFLAGSTMP = -c -no-multibyte-chars -pg -I/opt/intel/cc/10.0.023/include
LDFLAGSTMP =  -pg ${LDFLAGSTMPL}
#--#else #
FFLAGSTMP = -c -I/opt/intel/fc/10.0.023/include
CFLAGSTMP = -c -no-multibyte-chars -I/opt/intel/cc/10.0.023/include
LDFLAGSTMP =  ${LDFLAGGSTMPL}
#--#endif#
RANLIB=ranlib

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g 
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g -O0 
CFLAGS = ${CFLAGSTMP} -g -c -DLINUX -D_FILE_OFFSET_BITS=64 
CXXFLAGS = ${CFLAGSTMP} -g -c
LDFLAGSTMP =  ${LDFLAGSTMPL} -g 
#--#else# 
FFLAGSV = ${FFLAGSTMP} -ip -O2
FFLAGSS = ${FFLAGSTMP} -ip -O2 
FFLAGSN = ${FFLAGSTMP} -ip -O1
FFLAGSN0 = ${FFLAGSTMP} -O0
CFLAGS  = ${CFLAGSTMP} -DLINUX -D_FILE_OFFSET_BITS=64
CXXFLAGS = ${CFLAGSTMP} -c -O2
LDFLAGSTMP =  ${LDFLAGSTMPL} -O2
#--#endif#

#--#if static#
#---#if gcc#
LDFLAGS  = ${LDFLAGSTMP} -Wl,-Map load.map -D_FILE_OFFSET_BITS=64 -static
#---#else#
LDFLAGS  = ${LDFLAGSTMP} -Wl,-Map load.map -D_FILE_OFFSET_BITS=64 -static -static-libcxa
#---#endif gcc#
#--#else#
LDFLAGS  = ${LDFLAGSTMP} -D_FILE_OFFSET_BITS=64
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
OPTIONS=${MACHOPT}${BLASOPT}
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
#cat /proc/cpuinfo
#processor       : 0
#vendor_id       : GenuineIntel
#cpu family      : 6
#model           : 15
#model name      : Intel(R) Core(TM)2 Duo CPU     T7300  @ 2.00GHz
#stepping        : 10
#cpu MHz         : 2001.000
#cache size      : 4096 KB
#physical id     : 0
#siblings        : 2
#core id         : 0
#cpu cores       : 2
#fdiv_bug        : no
#hlt_bug         : no
#f00f_bug        : no
#coma_bug        : no
#fpu             : yes
#fpu_exception   : yes
#cpuid level     : 10
#wp              : yes
#flags           : fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe lm constant_tsc pni monitor ds_cpl vmx est tm2 ssse3 cx16 xtpr lahf_lm
#bogomips        : 3993.17
#clflush size    : 64
#
#processor       : 1
#vendor_id       : GenuineIntel
#cpu family      : 6
#model           : 15
#model name      : Intel(R) Core(TM)2 Duo CPU     T7300  @ 2.00GHz
#stepping        : 10
#cpu MHz         : 2001.000
#cache size      : 4096 KB
#physical id     : 0
#siblings        : 2
#core id         : 1
#cpu cores       : 2
#fdiv_bug        : no
#hlt_bug         : no
#f00f_bug        : no
#coma_bug        : no
#fpu             : yes
#fpu_exception   : yes
#cpuid level     : 10
#wp              : yes
#flags           : fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe lm constant_tsc pni monitor ds_cpl vmx est tm2 ssse3 cx16 xtpr lahf_lm
#bogomips        : 3990.13
#clflush size    : 64
#
# 
# Output of uname -a
#Linux lapmarc 2.6.22.9 #2 SMP Sat Oct 6 23:22:23 CEST 2007 i686 i686 i386 GNU/Linux
# 
# Using Intel Compilers, version 10.0
# Output of: /opt/intel/fc/10.0.023/bin/ifort -V
#Intel(R) Fortran Compiler for applications running on IA-32, Version 10.0    Build 20070426 Package ID: l_fc_p_10.0.023
#Copyright (C) 1985-2007 Intel Corporation.  All rights reserved.

# Output of: /opt/intel/cc/10.0.023/bin/icc -V
# Intel(R) C Compiler for applications running on IA-32, Version 10.0    Build 20070426 Package ID: l_cc_p_10.0.023
#Copyright (C) 1985-2007 Intel Corporation.  All rights reserved.
# 
