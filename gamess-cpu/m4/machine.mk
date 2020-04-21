#doc  MK file for the serial build of GAMESS-UK on an Opteron processor with
#doc  PGI compilers (tested with 5.2, 6.0 and 6.1)
#doc 
#doc  Options:
#doc  blas - build against blas (expected in: /usr/pgi/linux86-64/6.1/lib/libacml.a)
#doc  i8   - build with 8-byte integers.
##
# Default options:
#dopt vb nbo mrdci zora drf vdw sysmo dl-find
#opt blas static i8 mopac
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=opteron,linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv

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
FC = pgf90
FC90 = pgf90
LD = pgf90
CC = pgcc

#--#if i8#
FFLAGSI8=-i8
CFLAGSI8=-DLONG_INTEGER
M4_I8_OPT=,i8
#--#else #
FFLAGSI8=
CFLAGSI8=
M4_I8_OPT=,64bitpointers
#--#endif#

# Default compilation flags
FFLAGSD = -c -pg ${FFLAGSI8}
CFLAGSD= -c -pg ${CFLAGSI8}

#--#if debug#
FFLAGSV = ${FFLAGSD} -g
FFLAGSS = ${FFLAGSD} -g
FFLAGSN = ${FFLAGSD} -g
FFLAGSN0 = ${FFLAGSD} -g
CFLAGS = ${CFLAGSD} -g
LDTMP  = -g
#--#else# 
FFLAGSV =  ${FFLAGSD} -tp=k8-64 -Mcache_align -O2
FFLAGSS =  ${FFLAGSD} -tp=k8-64 -Mcache_align -O
FFLAGSN =  ${FFLAGSD} -tp=k8-64 -Mcache_align -O1
FFLAGSN0 = ${FFLAGSD}
CFLAGS  = ${CFLAGSD} -tp=k8-64 -D_REENTRANT
LDTMP = -tp=k8-64
#--#endif debug#
#
#--#if static#
LDFLAGS=${LDTMP} -Bstatic
#--#else#
LDFLAGS=${LDTMP}
#--#endif static#
#
#
#--#if blas#
BLASOPT=,blas
#acml needs to be linked statically for it to work.
LIBBLAS=/usr/pgi/linux86-64/6.1/lib/libacml.a
#--#endif blas#
#
OPTIONS=${MACHOPT}${BLASOPT}${M4_I8_OPT}
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

sp0011.o:	integs.m
	cat  ../utilities/gener.m integs.m | $(M4) -DGEN_EXTRACTFILE=sp0011 $(M4OPTS) > sp0011.f
	$(FC) $(FFLAGSN) sp0011.f

chfeq.o:	cphf.m
	cat  ../utilities/gener.m cphf.m | $(M4) -DGEN_EXTRACTFILE=chfeq $(M4OPTS) > chfeq.f
	$(FC) $(FFLAGSN) chfeq.f

bmove.o:	mclr.m
	cat  ../utilities/gener.m mclr.m | $(M4) -DGEN_EXTRACTFILE=bmove $(M4OPTS) > bmove.f
	$(FC) $(FFLAGSN) bmove.f

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

## Build Documentation
#
# Output of uname -a
# Linux scarf 2.4.21-27.0.2.ELsmp #1 SMP Wed Jan 12 23:25:44 EST 2005 x86_64 x86_64 x86_64 GNU/Linux
# 
# Output of cat /proc/cpuinfo
# processor	: 0
# vendor_id	: AuthenticAMD
# cpu family	: 15
# model		: 5
# model name	: AMD Opteron(tm) Processor 248
# physical id	: 0
# siblings	: 1
# stepping	: 10
# cpu MHz		: 2193.497
# cache size	: 1024 KB
# fpu		: yes
# fpu_exception	: yes
# cpuid level	: 1
# wp		: yes
# flags		: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 syscall nx mmxext lm 3dnowext 3dnow
# bogomips	: 4377.80
# TLB size	: 1088 4K pages
# clflush size	: 64
# address sizes	: 40 bits physical, 48 bits virtual
# power management: ts fid vid ttp
# 
# Output of cat /etc/issue
# Red Hat Enterprise Linux AS release 3 (Taroon Update 4)
# 
# Output of /usr/pgi/linux86-64/5.2/bin/pgf77 -V
# pgf77 5.2-4
# Copyright 1989-2000, The Portland Group, Inc.  All Rights Reserved.
# Copyright 2000-2004, STMicroelectronics, Inc.  All Rights Reserved.
# 
# Output of /usr/pgi/linux86-64/5.2/bin/pgfcc -V
# pgcc 5.2-4
# Copyright 1989-2000, The Portland Group, Inc.  All Rights Reserved.
# Copyright 2000-2004, STMicroelectronics, Inc.  All Rights Reserved.
#
