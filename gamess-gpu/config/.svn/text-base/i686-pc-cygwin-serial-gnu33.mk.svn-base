#doc This is the machine-specific file for the serial Windows build 
#doc on a Pentium 4 with the g77 under cygwin
#
# 
#dopt mrdci drf vb zora nbo mopac sysmo f77 
#opt debug static
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,bits8

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
CFLAGS = -g -c -DLINUX -D_FILE_OFFSET_BITS=64 -DLINUXF2C
#--#else# 
FFLAGSV = ${FFLAGSTMP} -O3 -fno-globals -malign-double
FFLAGSS = ${FFLAGSTMP} -O -fno-globals -malign-double
FFLAGSN = ${FFLAGSTMP} -O1 -fno-globals -malign-double
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS  = -c -DLINUX -D_FILE_OFFSET_BITS=64 -DLINUXF2C
#--#endif#
# Need -fno-globals for drf/drfsub.f or else g77 crashes without warning

#--#if static#
LDFLAGS  = -static -static-libgcc -D_FILE_OFFSET_BITS=64
#--#else#
LDFLAGS  = -D_FILE_OFFSET_BITS=64
#--#endif#

OPTIONS=${MACHOPT}

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
# Information on the machine this build was configured on:
#
# System Information
# OS Name	Microsoft Windows 2000 Professional
# Version	5.0.2195 Service Pack 4 Build 2195
# OS Manufacturer	Microsoft Corporation
# System Manufacturer	Viglen                         
# System Model	GENIE                          
# System Type	X86-based PC
# Processor	x86 Family 15 Model 2 Stepping 4 GenuineIntel ~1999 Mhz
# BIOS Version	BIOS Date: 07/03/02 14:09:02  Ver: 08.00.00
#
# Output of uname -a:
# CYGWIN_NT-5.0 MIAVIG4 1.5.16(0.128/4/2) 2005-04-25 20:26 i686 unknown unknown Cygwin
# 
# Output of: g77 -dumpversion
# GNU Fortran (GCC) 3.3.3 (cygwin special)
# Copyright (C) 2002 Free Software Foundation, Inc.
# 
# Output of: g77 -dumpmachine
# i686-pc-cygwin
#  
# Output of: g77 -print-libgcc-file-name
# /usr/lib/gcc-lib/i686-pc-cygwin/3.3.3/libgcc.a
#  
