#doc  This is the machine-specific file for the MPI parallel version
#doc  of GAMESS-UK using Intel compiler 8.1  running on 2GHz Xeon CPU's
#doc  with MPICH 1.2.6 compiled with ifc & icc
#doc
#doc  Options:
#doc  mkl        - use Intel mkl libraries ( expected in /opt/intel/mkl61/lib/32 )
#doc  goto       - use goto blas ( expected in /usr/lib )
#doc  score      - build on an SCORE system using the score wrappers
#doc  mpiwrap    - build using the mpixxx wrappers for the compiler - otherwise
#doc               mpich is expected in: /usr/local/mpich-1.2.5
#doc  newscf     - build the parallel ScaLAPACK SCF/DFT driver. BLACS & ScaLAPACK
#doc               are required and expected in: /usr/local/lib/mpich (see SCALIB variable)
#
#
#dopt base newscf
#opt mkl score mpiwrap goto vdw
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,upck-equiv,i8drct

#
# ===============  Compiler Options
#
## INTEL COMPILER VERSION 8.1 ##
#--#if score#
FC = mpif90 -compiler intel8
LD = mpif90 -compiler intel8
FC90 = mpif90 -compiler intel8
CC=  mpicc -compiler intel8
CXX = mpicc -compiler intel8
# LD flag -shared does not work with SCORE 5.8.2 so we need to use -nostatic or else
# the gamess-uk binary is not recognised as an SCORE application
LDSCORE = -nostatic
#--#elseif mpiwrap#
FC = mpif90
LD = mpif90
FC90 = mpif90
CC = mpicc
CXX = mpicc
#--#else#
FC = ifort 
LD = ifort 
FC90 = ifort
CC=  icc
CXX = icc
# MPI locations - Please set these for your system
MPI_INCLUDE = -I/usr/local/mpich-1.2.5/include
LIBMPI= -L/usr/local/mpich-1.2.5/lib -lmpich
#--#endif#

FFLAGSTMP = -c ${MPI_INCLUDE}
RANLIB=ranlib

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
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

#--#if static#
# This should work but doesn't...
LDFLAGS  = ${LDSCORE} -static -D_FILE_OFFSET_BITS=64
#--#else#
LDFLAGS  = ${LDSCORE} -Wl,-Map load.map
#--#endif#

#--#if mkl#
LIBBLAS=  -L/opt/intel/mkl61/lib/32  -lguide -lmkl -lmkl_lapack64
BLASOPT=,blas
#--#elseif goto#
LIBBLAS=-L/usr/lib -lgotoblas
BLASOPT=,blas
#--#endif#

#--#if newscf#
SCADIR= /usr/local/lib/mpich
SCALIB= ${SCADIR}/blacsF77init_MPI-LINUX-0.a ${SCADIR}/blacs_MPI-LINUX-0.a ${SCADIR}/libscalapack.a
#--#endif newscf#

# Bring all the options together
OPTIONS=${MACHOPT}${BLASOPT}
BL_LIB = ${SCALIB} ${LIBBLAS}  ${SCALIB} ${LIBMPI}

# Diesel build options
LD_DIESEL = ${CXX}
#--#if mkl#
DIESEL_LIBS = ${LIBBLAS} -lifcore
#--#elseif goto#
DIESEL_LIBS = ${LIBBLAS} -lifcore
#--#else#
DIESEL_LIBS = -lifcore
#--#endif#
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
#
# Information on the machine this build was configured on
# Curtailed output of /proc/cpuinfo
# vendor_id       : GenuineIntel
# cpu family      : 15
# model           : 2
# model name      : Intel(R) XEON(TM) CPU 2.00GHz
# stepping        : 4
# cpu MHz         : 1999.816
# cache size      : 512 KB
# 
# Output of uname -a
# Linux enterprise 2.4.20-28.7smp #1 SMP Thu Dec 18 11:18:31 EST 2003 i686 unknown
# 
# Using Intel Compilers, version 8.1
# Output of: ifc -V
# Intel(R) Fortran Compiler for 32-bit applications, Version 8.1    Build 20041019Z Package ID: l_fc_pu_8.1.021
# 
# Output of: icc -V
# Intel(R) C++ Compiler for 32-bit applications, Version 8.1    Build 20041019Z Package ID: l_cc_pu_8.1.024
# 
# Output of: /home/mpich.intel/bin/mpichversion
# MPICH Version:          1.2.6
# MPICH Release date:     $Date: 2007-08-17 15:46:11 $
# MPICH Patches applied:  none
# MPICH configure:        --prefix=/home/scott/mpich-1.2.6 -cc=icc -fc=ifort
# MPICH Device:           ch_p4
# 
