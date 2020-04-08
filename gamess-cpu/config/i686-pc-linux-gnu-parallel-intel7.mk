#doc Makefile settings for the parallel Global Array-based build
#doc on Pentium-class processors with Intel 7 compilers on an SCORE system
#doc The compiler libraries are expected to be found in: 
#doc /opt/intel/compiler70/ia32/lib 
#doc
#doc  Options:
#doc  mkl - use the Intel mkl libraries, expected /opt/intel/mkl61/lib/32
#doc  score - build for an score system using the score mpi wrappers
#doc  mpiwrap - use the mpich mpixx wrappers, other mpich is expected to
#doc            be found in /opt/mpich-1.2.6/
#
#dopt ga mpi i8 peigs mp2 zora vb vdw
#opt mkl score mpiwrap 
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,i8drct
#
# ===============  Compiler Options
#
## INTEL COMPILER VERSION 7 ##
#--#if score#
FC = mpif90 -compiler intel7
LD = mpif90 -compiler intel7
FC90 = mpif90 -compiler intel7
CC= mpicc -compiler intel7
#--#if mpiwrap
FC = mpif90
LD = mpif90
FC90 = mpif90
CC= mpicc
#--#elseif
FC = ifc 
LD = ifc
FC90 = ifc
CC= icc
#--#endif#

CXX= icc
FFLAGSTMP = -c
RANLIB=ranlib

#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g
FFLAGSS  = ${FFLAGSTMP} -g
FFLAGSN  = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS = -g -c -DLINUX -D_FILE_OFFSET_BITS=64 
CXXFLAGS = -g -c
#--#else# 
FFLAGSV = ${FFLAGSTMP} -ip -O3
FFLAGSS = ${FFLAGSTMP} -ip -O2 
FFLAGSN = ${FFLAGSTMP} -ip -O1
FFLAGSN0 = ${FFLAGSTMP} -O0
CFLAGS  = -c -DLINUX -D_FILE_OFFSET_BITS=64
CXXFLAGS = -c -O2
#--#endif#


#--#if mkl#
LIBBLAS= -L/opt/intel/mkl61/lib/32  -lguide -lmkl -lmkl_lapack64
BLASOPT=,blas
#--#endif#

LDFLAGS  = -D_FILE_OFFSET_BITS=64

LIBDIR=/opt/intel/compiler70/ia32/lib 
lLIB= -lPEPCF90


# GA and MPI variables
GA_F77_DEFS = -traditional
GA_TARGET=LINUX

#--if score mpiwrap#
MPI_INCLUDE = 
MPI_LIB =
MPI_LIBS =
GA_VERSION_PAR= GA_VERSION=SHMEM USE_MPI=YES

#--#else#
MPI_INCLUDE = /opt/mpich-1.2.6/include 
MPI_LIB =/opt/mpich-1.2.6/lib
LIBMPI =  -lmpich
MPI_LIBS =  -L${MPI_LIB} ${LIBMPI}
GA_VERSION_PAR= GA_VERSION=SHMEM USE_MPI=YES MPI_INCLUDE=${MPI_INCLUDE}  MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#--#endif#

OPTIONS=${MACHOPT}${BLASOPT}
BL_LIB = ${LLIB} ${lLIB}  ${LIBBLAS} ${MPI_LIBS}

# Diesel build options
LD_DIESEL = ${CXX}
#--#if mkl#
DIESEL_LIBS = -lifcore
#--#else#
DIESEL_LIBS = -llapack -lblas -lifcore
#--#endif mkl#
FLEX = /usr/bin/flex
__UNDERBAR = 1
SIZEOF_VOID_P = 4
LONG_LONG_INT = long long int
LONG_INT = long int
INT = int
SHORT_INT = short int

#
# ===============  Additional Files 
#
EXTRA=gethes.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o aprq34d.o mcdab_ga.o 

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

aprq34d.o:	mp2_parallel.m
	cat ../machines/$(MACH) ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=aprq34d $(M4OPTS) > aprq34d.f
	$(FC) $(FFLAGSS) aprq34d.f
