#doc Makefile settings for serial build on Itanium 2 running Linux 
#doc using the Intel compiler 7.0 (named efc)
#doc
#doc Options:
#doc mkl - use the Intel MKL numerical library (default location: /opt/intel/mkl701/lib/64)
#doc i8  -  build with integer*8 as opposed to the default integer*4
#
#dopt drf mrdci zora nbo vb vdw mopac sysmo dl-find i8
#opt mkl 
#
#
# ================ M4 Processing options
#

# Machine-specific options 
MACHINE_KEY=G
MACHOPT=itanium,linux,pclinux,littleendian,cio,unix,upck-equiv,64bitpointers,doublebackslash

IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)

#--#if mkl#
# BLAS options #
LIBBLAS = -L/opt/intel/mkl701/lib/64  -lmkl_ipf -lguide
BLASOPT=,blas
#--#endif#


#
# ===============  Compiler Options
#
FC = efc
LD = efc
CC = ecc -no-gcc
FC90 = ifort
CXX = icc
RANLIB = ranlib

#--#if i8#
FFLAGI8= -i8
CFLAGI8= -DLONG_INTEGER
I8_M4_OPT=,i8
#--#else#
FFLAGI8=
CFLAGI8=
I8_M4_OPT=,i8drct
#--#endif i8#
FFLAGSTMP = -c -ip -ftz ${FFLAGI8}

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS = -g -c ${CFLAGI8}
LDFLAGS  = -g -Vaxlib
CXXFLAGS = -g -c
#--#else# 
FFLAGSV = ${FFLAGSTMP} -O2 -cm -w90 -w95
FFLAGSS = ${FFLAGSTMP} -O2 -cm -w90 -w95
FFLAGSN = ${FFLAGSTMP} -O1 -cm -w90 -w95
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS  = -c   ${CFLAGI8}
LDFLAGS =   -Vaxlib
CXXFLAGS = -O2 -c
#--#endif#

# Diesel build options
LD_DIESEL = ${CXX}
#--#if mkl#
DIESEL_LIBS = ${LBLAS} ${lBLAS} -lifcore -lcxa -lunwind -lstdc++
#--#else#
DIESEL_LIBS = -lifcore -lcxa -lunwind -lstdc++
#--#endif#
FLEX = /usr/bin/flex
__UNDERBAR = 1
SIZEOF_VOID_P = 8
#--#if i8#
LONG_LONG_INT = long int
LONG_INT = long int
INT = long int
SHORT_INT = short int
#--#else#
LONG_LONG_INT = long long int
LONG_INT = long int
INT = int
SHORT_INT = short int
#--#endif i8#

# Put it all together
BL_LIB= ${MPI_LIBS} ${LIBBLAS}

## M4 Options ##
OPTIONS=${MACHOPT}${BLASOPT}${I8_M4_OPT}

#
# ===============  Additional Files 
#
EXTRA=gethes.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o mcdab_ga.o tst1s.o check0a.o
#
# ===============  Compiler Exceptions
#
tst1s.o:	direct.m
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f 
#
rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f
mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

gethes.o:	casb.m
	cat ../utilities/gener.m casb.m | $(M4) -DGEN_EXTRACTFILE=gethes $(M4OPTS) > gethes.f
	$(FC) $(FFLAGSN) gethes.f

mkmakw.o:	drvmp.m
	cat ../utilities/gener.m drvmp.m | $(M4) -DGEN_EXTRACTFILE=mkmakw $(M4OPTS) > mkmakw.f
	$(FC) $(FFLAGSS) mkmakw.f

mpmakw.o:	secmp2.m
	cat ../utilities/gener.m secmp2.m | $(M4) -DGEN_EXTRACTFILE=mpmakw $(M4OPTS) > mpmakw.f
	$(FC) $(FFLAGSS) mpmakw.f

umpe3a.o:	mp3.m
	cat ../utilities/gener.m mp3.m | $(M4) -DGEN_EXTRACTFILE=umpe3a $(M4OPTS) > umpe3a.f
	$(FC) $(FFLAGSS) umpe3a.f

umpe3b.o:	mp3.m
	cat ../utilities/gener.m mp3.m | $(M4) -DGEN_EXTRACTFILE=umpe3b $(M4OPTS) > umpe3b.f
	$(FC) $(FFLAGSS) umpe3b.f

check0a.o:	newmrd5.m
	cat ../utilities/gener.m newmrd5.m | $(M4) -DGEN_EXTRACTFILE=check0a $(M4OPTS) > check0a.f
	$(FC) $(FFLAGSN) check0a.f
aprq34d.o:	mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=aprq34d $(M4OPTS) > aprq34d.f
	$(FC) $(FFLAGSS) aprq34d.f

mcdab_ga.o:	mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=mcdab_ga $(M4OPTS) > mcdab_ga.f
	$(FC) $(FFLAGSS) mcdab_ga.f

#
# Information on the machine this build was configured on:
#
# Output of uname -a
# Linux tc2.chem.uu.nl 2.6.7-B64k.1.5 #1 SMP Sun Dec 12 17:33:17 CET 2004 ia64 ia64 ia64 GNU/Linux
#
# Output of cat /etc/issue
# Bull Linux Advanced Server release 3AS (Bull V1)
#
# Output of cat /proc/cpuinfo
# processor  : 0
# vendor     : GenuineIntel
# arch       : IA-64
# family     : Itanium 2
# model      : 2
# revision   : 1
# archrev    : 0
# features   : branchlong
# cpu number : 0
# cpu regs   : 4
# cpu MHz    : 1495.948995
# itc MHz    : 1495.948995
# BogoMIPS   : 2239.75
#
# Output of /opt/intel_fc_80/bin/ifort -V
# Intel(R) Fortran Itanium(R) Compiler for Itanium(R)-based applications
# Version 8.1    Build 20050203 Package ID: l_fc_pc_8.1.024
# Copyright (C) 1985-2005 Intel Corporation.  All rights reserved.
#
