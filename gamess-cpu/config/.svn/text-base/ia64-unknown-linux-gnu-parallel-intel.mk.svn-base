#doc Makefile settings for the GAMESS-UK parallel (GA) build on Itanium 2 running Linux
#doc using the Intel compiler as specified in the module command
#doc 
#doc Options:
#doc i8  - specify a build with integer*8 as opposed to the default integer*4
#doc mpiwrap - use mpif90, mpicc etc. wrappers
#doc           (otherwise MPICH is expected to be located in /usr/local/mpi_bull2-0.2-mono/intel_v8)
#doc mkl - use the Intel MKL numerical library (default location: /opt/intel/mkl/10.0.010/lib/64 
#doc scalapack - link in the ScaLAPACK diagonaliser

#
#dopt ga mpi peigs drf ci zora nbo vb vdw dl-find i8 mkl
#opt mpiwrap mopac newscf scalapack debug 
#
#
# ================ M4 Processing options
#

# Machine-specific options 
MACHINE_KEY=G
MACHOPT=itanium,linux,pclinux,littleendian,cio,unix,upck-equiv,64bitpointers,altix

IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)

# MPI Libraries
#--#if mpiwrap#
MPI_INCLUDE =
MPI_LIB =
LIBMPI=
MPI_LIBS     =
#--#else#
MPI_INCLUDE = /usr/local/mpi_bull2-0.2-mono/intel_v8/include
MPI_LIB = /usr/local/mpi_bull2-0.2-mono/intel_v8/lib/
LIBMPI= -lmpich -lmpio -lcpuset -lzerocopy
MPI_LIBS     = -L ${MPI_LIB} ${LIBMPI}
#--#endif mpiwrap#

#GA Stuff
GA_TARGET=LINUX64
GA_F77_DEFS = -traditional
GA_VERSION_PAR= GA_VERSION=SHMEM  MPI_INCLUDE=${MPI_INCLUDE}  MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "

#--#if mkl#
#---#if i8#
# The below is a shameful hack of which we should be embarrassed.
# If we are in i8 mode and linking to an i4 mkl library, then as
# this is a littleendian architecture the correct bits get passed
# through when we pass an i8 integer to the i4 library. This trick
# only works for cases where we have single integer arguments.
# dgthr, dsctr and ddoti take integer array arguments and so
# we can't use the i4 mkl versions and need to select our own 
# i8-compiled ones.
#LIBBLAS = ../linalg/dsctr.o ../linalg/dgthr.o ../linalg/ddoti.o  -L/opt/intel/mkl/10.0.010/lib/64  -lmkl_ipf -lguide
LIBBLAS = -lmkl_intel_ilp64 -lmkl_ipf -lmkl_core -lguide
#---#else#
LIBBLAS = -L/opt/intel/mkl/10.0.010/lib/64 -lmkl_ipf -lguide
#---#endif i8#
BLASOPT=,blas
#--#endif#

#--#if scalapack#
#LIB_SCALAPACK= /opt/intel/mkl70cluster/lib/64/libmkl_scalapack.a /apps/libs/BLACS/1.1/LIB/blacsF77init_MPI-LINUX-0.a /apps/libs/BLACS/1.1/LIB/blacs_MPI-LINUX-0.a
LIB_SCALAPACK= /opt/intel/mkl70cluster/lib/64/libmkl_scalapack.a /opt/intel/mkl70cluster/lib/64/libmkl_blacsF77init.a /opt/intel/mkl70cluster/lib/64/libmkl_blacs.a
#--#endif scalapack#


# Put it all together
BL_LIB= ${MPI_LIBS} ${LIBBLAS}

## M4 Options ##
#--#if i8#
OPTIONS=${MACHOPT}${BLASOPT},i8
#--#else#
OPTIONS=${MACHOPT}${BLASOPT},i8drct,64bitpointers
#--#endif i8#

#
# ===============  Compiler Options
#
#--#if mpiwrap#
FC = mpif90
LD = mpif90
CC = mpicc -no-gcc
FC90 = mpif90
#--#else#
FC = ifort
LD = ifort
CC = icc -no-gcc
FC90 = ifort
FINC = -I${MPI_INCLUDE}
#--#endif mpiwrap#
CXX = icc
RANLIB = ar -s

#--#if i8#
FFLAGSI8= -i8
CFLAGSI8= -DLONG_INTEGER -DEXT_INT
I8_M4_OPT=,i8
#--#else#
FFLAGSI8=
CFLAGSI8= -DSTD_INT
#--#endif i8#
FFLAGSTMP = -c -ftz ${FFLAGSI8}
LDFLAGSTMP= -Wl,-Map,load.map

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS = -g -c ${CFLAGSI8}
LDFLAGS  = -g ${LDFLAGSTMP}
CXXFLAGS = -g -c
#--#else# 
FFLAGSV = ${FFLAGSTMP} -O2 -ip -cm -w90 -w95
FFLAGSS = ${FFLAGSTMP} -O2 -ip -cm -w90 -w95
FFLAGSN = ${FFLAGSTMP} -O1 -ip -cm -w90 -w95
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS  = -c   ${CFLAGSI8}
LDFLAGS = ${LDFLAGSTMP}
CXXFLAGS = -O2 -c
#--#endif#


# Diesel build options
LD_DIESEL = ${CXX}
#--#if mkl#
DIESEL_LIBS = ${LBLAS} ${lBLAS} -L/opt/intel/fc/10.1.008/lib -lifcore -lcxa -lunwind -lstdc++
#--#else#
DIESEL_LIBS = -L/opt/intel/fc/x10.1.008/lib -lifcore -lcxa -lunwind -lstdc++
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

#
# ===============  Additional Files 
#
EXTRA=gethes.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o aprq34d.o mcdab_ga.o tst1s.o check0a.o
#
# ===============  Compiler Exceptions
#
#
# rpa.o:  rpa.f
#	   $(FC) $(FFLAGSN) $*.f
# mrdci5.o:	 mrdci5.f
#	   $(FC) $(FFLAGSN) $*.f

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

mcdab_ga.o:     mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=mcdab_ga $(M4OPTS) > mcdab_ga.f
	$(FC) $(FFLAGSS) mcdab_ga.f

tst1s.o:	direct.m
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f 
#
vbin.o:	vbin.f
	$(FC) $(FFLAGSN) $*.f
#
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
# Intel(R) Fortran IA-64 Compiler for applications running on IA-64, Version 10.1    
# Build 20070913 Package ID: l_fc_p_10.1.008
# Copyright (C) 1985-2007 Intel Corporation.  All rights reserved.
