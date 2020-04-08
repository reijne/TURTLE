#doc Makefile settings for Global Array based build on Itanium 2 running Linux
#doc using the Intel compiler 8.1
#doc NB: If see the error: undefined reference to _ReadULEB
#doc when linking then uncomment the INTELBUG variable in this file
#doc
#doc Options:
#doc mpiwrap - use mpif90, mpicc etc. wrappers
#doc           (otherwise MPICH is expected to be located in /usr/local/mpi_bull2-0.2-mono/intel_v8)
#doc mkl      - build with Intel MKL ( expected in /opt/intel/mkl70/lib/64 )
#doc scalapack - link in the ScaLAPACK diagonaliser
#
#dopt ga mpi ga peigs ci zora nbo vb vdw dl-find i8
#opt debug mkl mpiwrap mopac newscf scalapack
#
#
# ================ M4 Processing options
#

# Machine-specific options 
MACHINE_KEY=G
MACHOPT=itanium,linux,pclinux,littleendian,cio,unix,upck-equiv,altix

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
#
# The below is a shameful hack of which we should be embarrassed.
# We are in i8 mode but the mkl library is i4. For arguments that take
# single integer arguments this works, but for the routines that take integer
# array arguments, we need to use our own i8 routines.
#
LIBBLAS = ../linalg/dsctr.o ../linalg/dgthr.o ../linalg/ddoti.o -L/opt/intel/mkl70/lib/64 -lmkl_ipf -lguide
#---#else#
LIBBLAS = -L/opt/intel/mkl70/lib/64 -lmkl_ipf -lguide
#---#endif#
BLASOPT=,blas
#--#endif#

#--#if scalapack#
#LIB_SCALAPACK= /opt/intel/mkl70cluster/lib/64/libmkl_scalapack.a /apps/libs/BLACS/1.1/LIB/blacsF77init_MPI-LINUX-0.a /apps/libs/BLACS/1.1/LIB/blacs_MPI-LINUX-0.a
LIB_SCALAPACK= /opt/intel/mkl70cluster/lib/64/libmkl_scalapack.a /opt/intel/mkl70cluster/lib/64/libmkl_blacsF77init.a /opt/intel/mkl70cluster/lib/64/libmkl_blacs.a
#--#endif scalapack#


# If you see an error similar to the following when linking:
# /opt/intel/ifc_80/lib/libcxa.so.6: undefined reference to `_ReadULEB'
# /opt/intel/ifc_80/lib/libcxa.so.6: undefined reference to `_ReadSLEB'
# Then uncomment the line below, changing the path to match the directory
# where the Intel library files can be found 
#INTELBUG= /opt/intel/ifc_80/lib/libcxa.a /opt/intel/ifc_80/lib/libunwind.a

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
RANLIB = ranlib

#--#if i8#
FFLAGSI8= -i8
CFLAGSI8= -DLONG_INTEGER -DEXT_INT
I8_M4_OPT=,i8
#--#else#
FFLAGSI8=
CFLAGSI8= -DSTD_INT
I8_M4_OPT=,64bitpointers
#--#endif#

FFLAGSTMP = -c -ftz ${FFLAGSI8}
CFLAGSTMP = -c ${CFLAGSI8}

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g ${FINC}
FFLAGSS = ${FFLAGSTMP} -g ${FINC}
FFLAGSN = ${FFLAGSTMP} -g ${FINC}
FFLAGSN0 = ${FFLAGSTMP} -g 
CFLAGS = -g ${CFLAGSTMP}
LDFLAGS  = -g
CXXFLAGS = -g -c
#--#else# 
FFLAGSV = ${FFLAGSTMP} -O2 -ip -cm -w90 -w95 ${FINC}
FFLAGSS = ${FFLAGSTMP} -O2 -ip -cm -w90 -w95 ${FINC}
FFLAGSN = ${FFLAGSTMP} -O1 -ip -cm -w90 -w95 ${FINC}
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS  = ${CFLAGSTMP}
LDFLAGS =  
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
LONG_LONG_INT = long long int
LONG_INT = long int
INT = int
SHORT_INT = short int
#--#else#
LONG_LONG_INT = long int
LONG_INT = long int
INT = long int
SHORT_INT = short int
#--#endif#


# Put the link-line options together 
BL_LIB= ${LIB_SCALAPACK} ${MPI_LIBS} ${LIBBLAS} ${INTELBUG}

## M4 Options ##
OPTIONS=${MACHOPT}${BLASOPT}${I8_M4_OPT}

#
# ===============  Additional Files 
#
EXTRA=gethes.o 
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o aprq34d.o mcdab_ga.o tst1s.o
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
