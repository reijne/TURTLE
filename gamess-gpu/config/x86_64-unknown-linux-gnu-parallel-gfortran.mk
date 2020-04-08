#doc  Makefile settings for Global Array build on commodity clusters
#doc  with Optern/Intel EM64T processors using the gfortran compiler
#doc
#doc  Options:
#doc  mpiwrap - use mpif90 etc to compile and do not manually set path to libs & includes
#doc  i8 - build with integer*8 (implies tcgmsg on top of mpi)
#doc  vdw - use vanderwaals terms
#doc  debug - add debugging information to the executable
#doc
#doc NB: dl-find is not supported with gfortran < 4.2 as it uses allocatable
#doc     arrays in type declarations which is an extension to the F90 standard:
#doc     http://coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2004-12/0452.html
#doc     newscf is also not supported as it uses complex character array constructors
#doc     Since version 4.4 leadz is a gfortran intrinsic: hence the extleadz flag...
#
# DEFAULT OPTIONS
#dopt ci vdw masscf peigs ga mpi mpiwrap i8 vb zora
#opt debug infinband openib scalapack 
#
# ================ M4 Processing options
#

# Machine-specific options 
MACHINE_KEY=G
MACHOPT=linux,littleendian,cio,unix,upck-equiv,opteron,glibc,extleadz,GIGA_DUMP

IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)


# Link in ScaLAPACK
#--#if scalapack#
SPCKDIR=/usr/local/lib64
NEWSCF90_LIB=${SPCKDIR}/libscalapack.a ${SPCKDIR}/blacsF77init_MPI-LINUX-0.a ${SPCKDIR}/blacs_MPI-LINUX-0.a
#--#endif scalapack#


#--#if mpiwrap#
MPI_INCLUDE=
MPI_LIB=
LIBMPI =
#--#else#
MPI_INCLUDE=/software/software/OpenMPI/1.8.6-gcccuda-2.7.11/include
MPI_LIB=/software/software/OpenMPI/1.8.6-gcccuda-2.7.11/lib
LIBMPI = -lmpi -lmpi_mpifh
#--#endif mpiwrap#


# IC_INCLUDE is the includes required
# to build against the interconnect
IC_INCLUDE=
LIB_EXTRA=
#--#if myrinet#
# MPI and GM Locations for MYRINET switch
GM_INCLUDE=/opt/gm/2.0.19/include
IC_INCLUDE=-I${GM_INCLUDE}
GM_LIB=/opt/gm/2.0.19/lib64
LIB_EXTRA=-L$(GM_LIB) -lgm
#
#--#elseif infiniband openib#
#
# MPI Locations for INFINIBAND switch
IB_INCLUDE=/opt/ibgd/driver/infinihost/include
IC_INCLUDE=-I${IB_INCLUDE}
IB_LIB=/opt/ibgd/driver/infinihost/lib64
LIB_EXTRA= -L${IB_LIB} -lvapi
#--#endif#


# PEIGS definitions
PEIGS_TARGET_CPU_PAR=x86_64

GA_TARGET=LINUX64
#--#if ga#
# GA definitions & MPI and Switch library locations
GA_F77_DEFS = -traditional
# Change below to compile the GA's with a different c-compiler
#GA_CC = gcc
#---#if myrinet#
GA_VERSION_PAR= GA_VERSION=SHMEM ARMCI_NETWORK=GM GM_INCLUDE=${GM_INCLUDE} GM_LIB=${GM_LIB} MPI_INCLUDE=${MPI_INCLUDE} MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#
#---#elseif infiniband openib#
#
GA_VERSION_PAR= GA_VERSION=SHMEM ARMCI_NETWORK=OPENIB IB_INCLUDE=${IB_INCLUDE} MPI_INCLUDE=${MPI_INCLUDE} MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#---#else#
#
# Default sockets version
GA_VERSION_PAR= GA_VERSION=SHMEM MPI_INCLUDE=${MPI_INCLUDE} MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
LIB_EXTRA=-lpthread
#---#endif#
#--#endif ga#


#
# ===============  Compiler Options
#
CPP = cpp -traditional
#--#if mpiwrap#
FC = mpif90
LD = mpif90
FC90 = mpif90
CC = mpicc
CXX= mpiCC
FINC =
#--#else#
FC = gfortran 
LD = gfortran 
FC90 = gfortran
CC=  gcc
CXX = gcc
FINC = -I${MPI_INCLUDE} ${IC_INCLUDE}
#--#endif#


#--#if i8#
FFLAGSI8 = -fdefault-integer-8
CFLAGSI8 = -DLONG_INTEGER -DEXT_INT
# Added linux64 to get tst1s extracted from direct.m
I8_M4_OPT=,i8,linux64
#--#else#
FFLAGSI8=
CFLAGSI8= -DSTD_INT
I8_M4_OPT=,64bitpointers
#--#endif mpi#

#--#if coverage#
FLAGSCOV = -fprofile-arcs -ftest-coverage
#--#endif#

#--#if profile#
FLAGSPG = -pg
#--#endif#

FFLAGSTMP = -c ${FFLAGSI8} ${FLAGSPG} ${FLAGSCOV} ${FINC} -fno-range-check -mcmodel=large
CFLAGSTMP = -c ${CFLAGSI8} ${FLAGSPG} ${FLAGSCOV} ${FINC} -mcmodel=large
LDFLAGSTMP= ${FLAGSPG} ${FLAGSCOV} -mcmodel=large


#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g
FFLAGSS  = ${FFLAGSTMP} -g
FFLAGSN  = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g 
CFLAGS   = ${CFLAGSTMP}  -g
LDFLAGS  = ${LDFLAGSTMP} -g 
#--#else# 
FFLAGSV  = ${FFLAGSTMP} -O2 -w
FFLAGSS  = ${FFLAGSTMP} -O2 -w
FFLAGSN  = ${FFLAGSTMP} -O1 -w
FFLAGSN0 = ${FFLAGSTMP} -O0
CFLAGS   = ${CFLAGSTMP}
LDFLAGS  = ${LDFLAGSTMP}
#--#endif#


# Put the link line together
BL_LIB = ${NEWSCF90_LIB} ${LIBMPI} ${LIB_EXTRA}

# M4 Options 
OPTIONS=${MACHOPT}${I8_M4_OPT}

#
# ===============  Additional Files 
#
#--#if i8#
EXTRAI8= tst1s.o
#--#else#
EXTRAI8=
#--#endif#

EXTRA=$(EXTRAI8)
EXTRA_MP2= aprq34d.o mcdab_ga.o
EXTRA_DFT = jkint_dft.o jkder_dft.o

#
# ===============  Compiler Exceptions
#
tst1s.o:	direct.m
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f 
#
# anala - -O1 needed for gfortran 4.3.2 was OK with version 4.1.2
anala.o:	anala.f
	$(FC) $(FFLAGSN) $*.f
analb.o:	analb.f
	$(FC) $(FFLAGSN0) $*.f
rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f
mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

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
#  ========== DFT Exceptions ===============
#
jkint_dft.o:	integ2e.m
	cat ../utilities/gener.m  integ2e.m | $(M4) -DGEN_EXTRACTFILE=jkint_dft $(M4OPTS) $(SNGL) > jkint_dft.f
	$(FC) $(FFLAGSN) jkint_dft.f
#
jkder_dft.o:	deriv2e.m
	cat ../utilities/gener.m  deriv2e.m | $(M4) -DGEN_EXTRACTFILE=jkder_dft $(M4OPTS) $(SNGL) > jkder_dft.f
	$(FC) $(FFLAGSN) jkder_dft.f


