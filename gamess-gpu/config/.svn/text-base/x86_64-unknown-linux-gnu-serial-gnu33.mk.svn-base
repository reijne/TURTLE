#doc  This is the machine dependant file for the serial build of GAMESS-UK
#doc  on an Opteron processor with g77 version 3.3.5
#doc
#doc  Options:
#doc  static - produce a statically linked binary
##
# Default options:
#dopt vb nbo mrdci zora drf sysmo f77
#opt static mopac
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=opteron,linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,64bitpointers,fio_access_append,glibc

# Bitwise operators that differ from the defaults.
IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)
# 
#
RANLIB=ranlib
CPP=/lib/cpp
GA_F77_DEFS = -traditional
#
# ===============  Compiler Options
#
#ensure 64bit compilers invoked?
FC = g77
LD = g77
FC90 = g77
CC = gcc
FFLAGSTMP = -c
#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS = -g -c -DLINUX  -DLINUXF2C  -D_FILE_OFFSET_BITS=64
#--#else# 
# rem removed -align-doubles
FFLAGSV = ${FFLAGSTMP} -O3 -fno-globals
FFLAGSS = ${FFLAGSTMP} -O -fno-globals
FFLAGSN = ${FFLAGSTMP} -O1 -fno-globals
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS = -c -DLINUX  -DLINUXF2C  -D_FILE_OFFSET_BITS=64
#--#endif#
# Need -fno-globals for drf/drfsub.f or else g77 crashes without warning
#
#--#if static#
LDFLAGS  = -static -static-libgcc -D_FILE_OFFSET_BITS=64
#--#else#
LDFLAGS  = -D_FILE_OFFSET_BITS=64
#--#endif#
#
#
#
OPTIONS=${MACHOPT}${BLASOPT}
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
	cat ../machines/$(MACH) ../utilities/gener.m  integ2e.m | $(M4) -DGEN_EXTRACTFILE=jkint_dft $(M4OPTS) $(SNGL) > jkint_dft.f
	$(FC) $(FFLAGSN) jkint_dft.f
#
jkder_dft.o:	deriv2e.m
	cat ../machines/$(MACH) ../utilities/gener.m  deriv2e.m | $(M4) -DGEN_EXTRACTFILE=jkder_dft $(M4OPTS) $(SNGL) > jkder_dft.f
	$(FC) $(FFLAGSN) jkder_dft.f
#
# Details of the machine this build was configured on (date=Fri Sep  2 15:43:04 BST 2005):
#
# output of uname -a:
# Linux ladon 2.6.11.4-21.7-smp #1 SMP Thu Jun 2 14:23:14 UTC 2005 x86_64 x86_64 x86_64 GNU/Linux
# 
# output of cat /etc/SuSE-release:
# SuSE Linux 9.3 (x86-64)
# VERSION = 9.3
# 
# output of cat /proc/cpuinfo:
# processor	: 0
# vendor_id	: AuthenticAMD
# cpu family	: 15
# model		: 5
# model name	: AMD Opteron(tm) Processor 246
# stepping	: 10
# cpu MHz		: 1994.377
# cache size	: 1024 KB
# physical id	: 255
# siblings	: 1
# fpu		: yes
# fpu_exception	: yes
# cpuid level	: 1
# wp		: yes
# flags		: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 pni syscall nx mmxext lm 3dnowext 3dnow
# bogomips	: 3923.96
# TLB size	: 1024 4K pages
# clflush size	: 64
# cache_alignment	: 64
# address sizes	: 40 bits physical, 48 bits virtual
# power management: ts fid vid ttp
# 
# output of gcc --version:
# gcc (GCC) 3.3.5 20050117 (prerelease) (SUSE Linux)
# 
# output of gcc  -print-libgcc-file-name
# /usr/lib64/gcc-lib/x86_64-suse-linux/3.3.5/libgcc.a
# 
# output of gcc -dumpmachine
# x86_64-suse-linux
