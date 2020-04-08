#doc Makefile settings for the GAMESS-UK serial build on Itanium 2 running Linux
#doc using the Intel compiler  as determined by module
#doc 
#doc Options:
#doc mkl - use the Intel MKL numerical library (default location: /opt/intel/mkl70/lib/64)
#doc i4  - specify a build with integer*4 as opposed to the default integer*8
#
#dopt drf mrdci zora nbo vb vdw sysmo mopac dl-find mkl
#opt  i4 debug 
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

#--#if mkl#
#---#if i4#
LIBBLAS =  -lmkl_ipf -lguide
#---#else#
# The below is a shameful hack of which we should be embarrassed.
# If we are in i8 mode and linking to an i4 mkl library, then as
# this is a littleendian architecture the correct bits get passed
# through when we pass an i8 integer to the i4 library. This trick
# only works for cases where we have single integer arguments.
# dgthr, dsctr and ddoti take integer array arguments and so
# we can't use the i4 mkl versions and need to select our own 
# i8-compiled ones.
# LIBBLAS = ../linalg/dsctr.o ../linalg/dgthr.o ../linalg/ddoti.o -lmkl_ipf -lguide
#  this uses the i*8 arguments
LIBBLAS = -lmkl_intel_ilp64 -lmkl_ipf -lmkl_core -lguide
#---#endif#
BLASOPT=,blas
#--#endif#
LIBMPI=

# Put it all together
LIB= ${MPI_LIBS} ${LIBBLAS}

## M4 Options ##
#--#if i4#
OPTIONS=${MACHOPT}${BLASOPT},i8drct,64bitpointers
#--#else#
OPTIONS=${MACHOPT}${BLASOPT},i8
#--#endif i4#

#
# ===============  Compiler Options
#
FC = ifort
LD = ifort
CC = icc -no-gcc
FC90 = ifort
CXX = icc
RANLIB = ar -s

#--#if i4#
FFLAGI8=
CFLAGI8=
#--#else#
FFLAGI8= -i8
CFLAGI8= -DLONG_INTEGER
#--#endif i4#
FFLAGSTMP = -c -ftz ${FFLAGI8}
LDFLAGSTMP= -Wl,-Map,load.map

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS = -g -c ${CFLAGI8}
LDFLAGS  = -g ${LDFLAGSTMP}
CXXFLAGS = -g -c
#--#else# 
FFLAGSV = ${FFLAGSTMP} -O2 -ip -cm -w90 -w95
FFLAGSS = ${FFLAGSTMP} -O2 -ip -cm -w90 -w95
FFLAGSN = ${FFLAGSTMP} -O1 -ip -cm -w90 -w95
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS  = -c   ${CFLAGI8}
LDFLAGS = ${LDFLAGSTMP}
CXXFLAGS = -O2 -c
#--#endif#

# Diesel build options
LD_DIESEL = ${CXX}
#--#if mkl#
DIESEL_LIBS = ${LBLAS} ${lBLAS} -L/opt/intel/fc/9.0/lib -lifcore -lcxa -lunwind -lstdc++
#--#else#
DIESEL_LIBS = -L/opt/intel/fc/9.0/lib -lifcore -lcxa -lunwind -lstdc++
#--#endif#
FLEX = /usr/bin/flex
__UNDERBAR = 1
SIZEOF_VOID_P = 8
#--#if i4#
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

LIBMPI = 
# Put the link line together
BL_LIB = ${SCALAPACK_LIB} ${LIBBLAS} ${LIBMPI} ${LIB_EXTRA}


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

vbin.o:	vbin.f
	$(FC) $(FFLAGSN) $*.f

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

# Compiler crashed on 9.0 (jvl oct05)

mclr.o:	mclr.m
	cat ../utilities/gener.m  mclr.m | $(M4) $(M4OPTS) > mclr.f
	$(FC)   $(FFLAGSS) -g mclr.f

# for intel 11 (jvl oct10):x
newmrd2.o:	newmrd2.f
	$(FC) $(FFLAGSN) $*.f


#
# Information on the machine this build was configured on:
#
# Output of uname -a
# Linux tc2.chem.uu.nl 2.6.7-B64k.1.5 #1 SMP Sun Dec 12 17:33:17 CET 2004 ia64 ia64 ia64 GNU/Linux
#
# Output of cat /etc/issue
# Bull Linux Advanced Server release 3AS (Bull V1) 
# Kernel \r on an \m
#
# processor  : 0 (of 8)
# vendor     : GenuineIntel
# arch       : IA-64
# family     : 32
# model      : 0
# revision   : 7
# archrev    : 0
# features   : branchlong, 16-byte atomic ops
# cpu number : 0
# cpu regs   : 4
# cpu MHz    : 1595.672000
# itc MHz    : 398.918000
# BogoMIPS   : 3187.67
#

