#doc Machine-dependant file for GAMESS-UK MPI build on Linux/Pentium
#doc with PGI compilers version 6.0-2 and mpich compiler wrappers.
#doc
#doc Options:
#doc blas - link against PGI blas library (default location: /opt/pgi/linux86/6.0/lib)
#
# DEFAULT AND OPTIONAL OPTIONS
#dopt base dynamic_lb
#opt blas debug static
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,i8drct

FC = mpif77
LD = mpif77
FC90 = mpif77
CC=mpicc
FFLAGSTMP = -c
RANLIB=ranlib

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS = -g -c -DLINUX -D_FILE_OFFSET_BITS=64 
#--#else# 
FFLAGSV = ${FFLAGSTMP} -fast
FFLAGSS = ${FFLAGSTMP} -O
FFLAGSN = ${FFLAGSTMP} -O1
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS  = -c -DLINUX -D_FILE_OFFSET_BITS=64
#--#endif#

#--#if static#
LDFLAGS  = -Bstatic -g77libs -D_FILE_OFFSET_BITS=64
#--#else#
LDFLAGS  = -D_FILE_OFFSET_BITS=64
#--#endif#

#--#if blas#
LIBBLAS=-L/opt/pgi/linux86/6.0/lib  -llapack -lblas
BLASOPT=,blas
#--#endif#

OPTIONS=${MACHOPT}${BLASOPT}

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
