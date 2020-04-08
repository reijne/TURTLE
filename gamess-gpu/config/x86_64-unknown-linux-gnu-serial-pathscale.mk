#doc  This is the machine dependant file for the serial build of GAMESS-UK
#doc  on Redhat Enterprise Linux AS3 running on an Opteron processor with
#doc  the pathscale compilers version 2.2.
#doc
#doc  Options:
#doc  static - create a statically linked binary
#doc  acml   - link against PGI's ACML blas library 
#doc           (default location:  ~/acml2.7.0/pathscale64/lib/libacml.a)
##
# Default options:
#dopt vb nbo mrdci zora drf vdw sysmo dl-find static 
#opt i8 acml mopac
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=opteron,linux,pclinux,littleendian,cio,unix,upck-equiv,64bitpointers,i8drct

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
#ensure 64bit compilers invoked
FC = /opt/pathscale/bin/pathf90
FC90 = /opt/pathscale/bin/pathf90
LD = /opt/pathscale/bin/pathf90
CC = /opt/pathscale/bin/pathcc

#Default compilation Flags
#--#if i8#
FFLAGSI8 = -i8
CFLAGSI8 = -DLONG_INTEGER
#--#endif i8#

FFLAGSD = -fno-second-underscore ${FFLAGSI8}
CFLAGSD =  ${CFLAGSI8}

#--#if debug#
FFLAGSV = -c -g
FFLAGSS = -c -g
FFLAGSN = -c -g
CFLAGS = -g -c
LDTMP  = -g
#--#else# 
FFLAGSV = -c ${FFLAGSD}  -O2
FFLAGSS = -c ${FFLAGSD}  -O
FFLAGSN = -c ${FFLAGSD}  -O1
FFLAGSN0 = -c
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
#
#--#if acml#
BLASOPT=,blas
#acml needs to be linked statically for it to work.
LIBBLAS= ~/acml2.7.0/pathscale64/lib/libacml.a
#--#endif acml#
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
	cat ../utilities/gener.m  integ2e.m | $(M4) -DGEN_EXTRACTFILE=jkint_dft $(M4OPTS) $(SNGL) > jkint_dft.f
	$(FC) $(FFLAGSN) jkint_dft.f
#
jkder_dft.o:	deriv2e.m
	cat ../utilities/gener.m  deriv2e.m | $(M4) -DGEN_EXTRACTFILE=jkder_dft $(M4OPTS) $(SNGL) > jkder_dft.f
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
