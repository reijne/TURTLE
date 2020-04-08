#doc This is the machine-specific file for the serial Windows build 
#doc on a Pentium 4 with g95 under cygwin
# 
#dopt mrdci drf vb zora nbo vdw mopac sysmo dl-find
#opt debug static
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,bits8,GFS

IAND32=iand($$1,$$2)
IOR32=ior($$1,$$2)
IXOR32=ieor($$1,$$2)
ISHFT32=ishft($$1,$$2)

FC = gfortran
LD = gfortran
FC90 = gfortran
CC = gcc
FFLAGSTMP = -c -fno-range-check -fno-second-underscore -DGFS
#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS = -g -c -DLINUX -D_FILE_OFFSET_BITS=64 
#--#else# 
FFLAGSV = ${FFLAGSTMP} 
FFLAGSS = ${FFLAGSTMP} 
FFLAGSN = ${FFLAGSTMP} 
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS  = -c -DLINUX -D_FILE_OFFSET_BITS=64 
#--#endif#

#--#if dl-find#
# Need to set or else we can't make the dl-find library
GMAKE=/usr/bin/make
#--#endif dl-find#

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
# OS Name	Microsoft Windows XP Professional
# Version  5.1, build 2600.xpsp_sp2_rtm.040803-2158: SP2 
# OS Manufacturer	Microsoft Corporation
# System Manufacturer	Toshiba
# System Model	Tecra
# Processor	Pentium M 1.73 Ghz
#
# Output of uname -a:
# CYGWIN_NT-5.1 ps96port1 1.5.19(0.150/4/2) 2006-01-20 13:28 i686 Cygwin
# 
# Output of: g95 -dumpversion
#
# G95 (GCC 4.0.3 (g95!) Jul  2 2006)
# Copyright (C) 2002-2005 Free Software Foundation, Inc.
#
# G95 comes with NO WARRANTY, to the extent permitted by law.
# You may redistribute copies of G95
# under the terms of the GNU General Public License.
# For more information about these matters, see the file named COPYING


