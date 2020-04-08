#doc Machine-dependant file for GAMESS-UK serial build on Linux/Pentium
#doc with PGI compilers version 6.0-2. Compilers and libararies are
#doc expected to be found in: /opt/pgi/linux86/6.0
#doc
#doc Options:
#doc acml   - use pgi blas library in /opt/pgi/linux86/6.2/lib
#doc mopac  - include mopac code in the build
#doc debug  - include debugging information in the objects (no optimsation)
#doc static - build a static binary
#doc xml    - build AgentX to allow the code to read/write xml (Experimental)
#
# DEFAULT AND OPTIONAL OPTIONS
#dopt zora vb drf nbo mrdci vdw sysmo dl-find
#opt acml  debug static xml mopac
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,i8drct

#--#if charmm#
# currently CHARMM builds with pgf77 so GAMESS-UK follows suit (this may change)
FC = pgf77
LD = pgf77
#FC90 = pgf90
#--#else#
FC = pgf90
LD = pgf90
FC90 = pgf90
#--#endif#
CC=pgcc

FFLAGSTMP = -c
CFLAGSTMP = -c -DLINUX -D_FILE_OFFSET_BITS=64 
LDFLAGSTMP=-Wl,-Map,link.map

RANLIB=ranlib

#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g -Ktrap=fp
FFLAGSS  = ${FFLAGSTMP} -g -Ktrap=fp
FFLAGSN  = ${FFLAGSTMP} -g -Ktrap=fp
FFLAGSN0 = ${FFLAGSTMP} -g -Ktrap=fp
CFLAGS   = ${CFLAGSTMP} -g
LDDEBUG=-g -Ktrap=fp 
#--#else# 
FFLAGSV  = ${FFLAGSTMP} -fastsse
FFLAGSS  = ${FFLAGSTMP} -O
FFLAGSN  = ${FFLAGSTMP} -O1
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS   = ${CFLAGSTMP}
LDDEBUG=
#--#endif#

#--#if static#
LDFLAGSSTATIC= -Bstatic -g77libs
#--#endif#
LDFLAGS  =  $(LDFLAGSTMP) $(LDDEBUG) $(LDFLAGSSTATIC)

#--#if acml#
LIBBLAS=-L/opt/pgi/linux86/6.2/lib -lacml
# Below causes segmentation faults
#lBLAS= -lacml -lpgsse2
#--#endif acml#

OPTIONS=${MACHOPT}

# Bring all the options together
BL_LIB = ${LIBBLAS}

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

# dlamch.o:	dlapack.m
# 	cat ../utilities/gener.m  dlapack.m | $(M4)  -DGEN_EXTRACTFILE=dlamch  $(M4OPTS)  > dlamch.f
# 	$(FC) $(FFLAGSN0) dlamch.f

tst1s.o:	direct.m
	cat ../machines/$(MACH) ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f 

rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

#wrong results for c2029_a
tsort.o:	tsort.f
	$(FC) $(FFLAGSN) $*.f

#
# Information on the machine this build was configured on
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
# 
# Output of cat /etc/issue
# Welcome to SuSE Linux 9.0 (i586)
# 
# Output of /opt/pgi/linux86/6.0/bin/pgf90 -V
# pgf90 6.0-2 32-bit target on x86 Linux
# Copyright 1989-2000, The Portland Group, Inc.  All Rights Reserved.
# Copyright 2000-2005, STMicroelectronics, Inc.  All Rights Reserved.
