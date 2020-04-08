#doc This is the machine-specific file for the serial Windows build 
#doc on a Pentium 4 with the Intel compiler version 7.1 and optional mkl version 6.0
#doc the compiler libraries are expected to be found in E:\Program Files\Intel\Compiler70\IA32\
#doc and the MKL libraries in E:\Program Files\Intel\MKL60\ia32\lib
# 
#dopt mrdci zora vb drf nbo vdw sysmo mopac dl-find
#opt blas


# Rename library extensions and add the flag so we can use .o and not .obj
OBJNAME= -Fo$*.o
# LD doesn't accept -o flag so...
LDNAME=-Fe

# Machine-specific 
MACHINE_KEY=G
IAND32=iand($$1,$$2)
IOR32=ior($$1,$$2)
IXOR32=ieor($$1,$$2)
SHIFT=ishift($$1,$$2)

# M4 Options
MACHOPT=linux,win95,ifc,littleendian,cio,unix,doublebackslash,upck-equiv,i8drct

# Compiler options
FC = ifl
LD = ifl
CC = icl
FFLAGSV  = -c -w90 -w95 -Qip -O3
FFLAGSS  = -c -w90 -w95 -Qip -O2
FFLAGSN  = -c -w90 -w95 -Qip -O1
FFLAGSN0 = -c -w90 -w95 -Qip -O1
CFLAGS   = -c -Qip -O3 -I"E:\Program Files\Microsoft Visual Studio .NET\Vc7\include"

LDBASE= -Qip -Qoption,link,-map,-LIBPATH:"E:\Program Files\Intel\Compiler70\IA32\Lib",-LIBPATH="E:\Program Files\Microsoft Visual Studio .NET\Vc7\lib",-LIBPATH="E:\Program Files\Microsoft Visual Studio .NET\Vc7\PlatformSDK\lib"

#--#if blas#
LDBLAS = -Qoption,link,-LIBPATH:"E:\Program Files\Intel\MKL60\ia32\lib"
LIBBLAS =  mkl_c_dll.lib
BLASOPT=,blas
#--#endif blas#

LDFLAGS= ${LDBASE} ${LDBLAS}

OPTIONS=${MACHOPT}$(BLASOPT)
BL_LIB = ${LIBBLAS}

EXTRA_BASE= bitmanip.o
EXTRA_BENCH= bitmanip.o
# Need to check about the below
#EXTRA=gethes.o
#EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o
EXTRA_MP2= check0a.o

#
#
#  ================ Exceptions for linux (g77, pgf77, absoft, efc) =======
#
#
tst1s.o:	direct.m
	cat  ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) -Fotst1s.o tst1s.f 
#
rpa.o:	rpa.f
	$(FC) $(FFLAGSN) -Fo$*.o $*.f
mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) -Fo$*.o $*.f

gethes.o:	casb.m
	cat  ../utilities/gener.m casb.m | $(M4) -DGEN_EXTRACTFILE=gethes $(M4OPTS) > gethes.f
	$(FC) $(FFLAGSN) -Fogethes.o gethes.f

mkmakw.o:	drvmp.m
	cat  ../utilities/gener.m drvmp.m | $(M4) -DGEN_EXTRACTFILE=mkmakw $(M4OPTS) > mkmakw.f
	$(FC) $(FFLAGSS) -Fomkmakw.o mkmakw.f

mpmakw.o:	secmp2.m
	cat  ../utilities/gener.m secmp2.m | $(M4) -DGEN_EXTRACTFILE=mpmakw $(M4OPTS) > mpmakw.f
	$(FC) $(FFLAGSS) -Fompmakw.o mpmakw.f

umpe3a.o:	mp3.m
	cat  ../utilities/gener.m mp3.m | $(M4) -DGEN_EXTRACTFILE=umpe3a $(M4OPTS) > umpe3a.f
	$(FC) $(FFLAGSS) -Foumpe3a.o umpe3a.f

umpe3b.o:	mp3.m
	cat  ../utilities/gener.m mp3.m | $(M4) -DGEN_EXTRACTFILE=umpe3b $(M4OPTS) > umpe3b.f
	$(FC) $(FFLAGSS) -Foumpe3b.o umpe3b.f

check0a.o:	newmrd5.m
	cat ../utilities/gener.m newmrd5.m | $(M4) -DGEN_EXTRACTFILE=check0a $(M4OPTS) > check0a.f
	$(FC) $(FFLAGSN) -Focheck0a.o check0a.f
#
#
#
#
#
#  ====================== Exceptions for Windows95 =======================
#
#parallel.f:	parallel.m
#	rm -f parallel.f
#	cat ../utilities/gener.m parallel.m | $(M4) $(M4OPTS) $(SNGL) | \
#	../utilities/quote >> parallel.f
#
bitmanip.f:	machscf.m
	cat ../utilities/gener.m machscf.m | $(M4) -DGEN_EXTRACTFILE=bitmanip \
	$(M4OPTS) $(SNGL) | ../utilities/quote > bitmanip.f
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
# Information on the machine follows
# Output of ifl -v
# Intel(R) Fortran Compiler for 32-bit applications, Version 7.1   Build 20030307Z
# Copyright (C) 1985-2003 Intel Corporation.  All rights reserved.
#
# Output of icl -V
# Intel(R) C++ Compiler for 32-bit applications, Version 7.1   Build 20030307Z
# Copyright (C) 1985-2003 Intel Corporation.  All rights reserved.
#
