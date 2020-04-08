#doc Machine-dependant file for GAMESS-UK serial build on Linux/Pentium
#doc with GNU compilers version 3.3.1.
#doc
#doc Options:
#doc static - create a statically linked binary
#doc diesel - include the C++ DIESEL CI programm
#
# DEFAULT AND OPTIONAL OPTIONS
#dopt zora vb drf nbo mrdci sysmo f77 
#opt demo debug static diesel mopac
#
# ================ M4 Processing options
#
#
#last option added for gcc 322 version of g77 build with charmm
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,bits8,fio_access_append,glibc

FC = g77
LD = g77
FC90 = g77
CC = gcc
#--#if diesel#
CXX = g++
#--#endif diesel#
#--#if profile#
FFLAGSTMP  = -pg -c -fno-globals -Wno-globals
CFLAGSTMP  = -pg
LDFLAGSTMP = -pg
#--#else#
FFLAGSTMP  = -c -fno-globals -Wno-globals
CFLAGSTMP  = 
LDFLAGSTMP =
#--#endif#
RANLIB = ranlib
#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g
FFLAGSS  = ${FFLAGSTMP} -g
FFLAGSN  = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS   = ${CFLAGSTMP} -g -c -DLINUX -D_FILE_OFFSET_BITS=64 -DLINUXF2C
CXXFLAGS = ${CFLAGSTMP} -g -c 
#--#else# 
FFLAGSV  = ${FFLAGSTMP} -O3 -malign-double
FFLAGSS  = ${FFLAGSTMP} -O  -malign-double
FFLAGSN  = ${FFLAGSTMP} -O1 -malign-double
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS   = ${CFLAGSTMP} -c -DLINUX -D_FILE_OFFSET_BITS=64 -DLINUXF2C
CXXFLAGS = ${CFLAGSTMP} -c -O2
#--#endif#
# Need -fno-globals for drf/drfsub.f or else g77 crashes without warning

#--#if static#
LDFLAGS  = ${LDFLAGSTMP} -static -static-libgcc -D_FILE_OFFSET_BITS=64
#--#else#
LDFLAGS  = ${LDFLAGSTMP} -D_FILE_OFFSET_BITS=64
#--#endif#

OPTIONS=${MACHOPT}

#--#if diesel#
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
#--#endif diesel#

#
# ===============  Additional Files 
#
EXTRA=gethes.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o

#
# ===============  Compiler Exceptions
#
##
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
# 
# Output of cat /etc/issue
# Welcome to SuSE Linux 9.0 (i586)
#
# Output of /usr/bin/g77 --version
# GNU Fortran (GCC) 3.3.1 (SuSE Linux)
# Copyright (C) 2002 Free Software Foundation, Inc.
# 
# Output of /usr/bin/gcc --version
# gcc (GCC) 3.3.1 (SuSE Linux)
# Copyright (C) 2003 Free Software Foundation, Inc.
