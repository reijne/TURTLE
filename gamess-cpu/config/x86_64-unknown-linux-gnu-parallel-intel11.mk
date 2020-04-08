#doc  Makefile settings for Global Array build on commodity clusters
#doc  with Optern/Intel EM64T processors using the Intel compiler 11 (11.1/072
#doc
#doc adds exception for helfey
#doc
#doc  Options:
#doc  mpiwrap - use mpif90 etc to compile and do not manually set path to libs & includes
#doc  myrinet - build for Myrinet switch. Default library locations:
#doc            libgm.a:    /opt/gm/2.0.19/lib64
#doc            lib*mpi.a:  /opt/MPI/intel8/mpich-gm/1.2.6..14/lib 
#doc
#doc  infiniband - build for INFINIBAND switch. Default library locations:
#doc              libvapi.a:     /opt/ibgd/driver/infinihost/lib64
#doc              lib*mpi.a:     /opt/MPI/intel8/mvapich/0.9.5/lib 
#doc
#doc  mkl - build against Intel MKL library - NOT AVAILABLE WITH I8
#doc  scalapack - link to external scalapack library - NOT AVAILABLE WITH I8
#doc  mkl_scalapack - build with Intel MKL BLAS and ScaLAPACK - NOT AVAILABLE WITH I8
#doc  score  - build on an SCORE system
#doc  newscf - build distributed scf driver - NOT AVAILABLE WITH i8
#doc  vdw - use vanderwaals terms
#doc  datain - read input from a file named "datain" as opposed to stdin
#doc  debug - add debugging information to the executable
#
# DEFAULT OPTIONS
#dopt ga mpi ci peigs vdw mpiwrap dl-find masscf vb zora i8 drf
#opt debug score myrinet infiniband openib mkl scalapack newscf i8 datain
#
# ================ M4 Processing options
#

# Machine-specific options 
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,upck-equiv,helfey

IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)

#--#if mkl_scalapack#
# 
# Using both the MKL BLAS snd BLACS/ScaLAPACK
SCALAPACK_HOME=/opt/intel/mkl/10.0.1.014/lib/em64t
#SCALAPACK_LIB=-Wl,--start-group $(SCALAPACK_HOME)/libmkl_scalapack_lp64.a $(SCALAPACK_HOME)/libmkl_blacs_intelmpi_lp64.a $(SCALAPACK_HOME)/libmkl_intel_lp64.a $(SCALAPACK_HOME)/libmkl_intel_thread.a  $(SCALAPACK_HOME)/libmkl_core.a $(SCALAPACK_HOME)/libguide.a -lpthread -Wl,--end-group
BLASOPT=,blas
#
# #--#else#
#
#---#if mkl#
# -lguide is used to pick the threading library for dynamic linking
LIBBLAS=-L/software/software/ifort/2015.2.164/composerxe/mkl/lib/intel64 -lmkl_core -lmkl_intel_ilp64 -lpthread -lmkl -lmkl_scalapack_ilp64 -lmkl_sequential -lmkl_blacs_intelmpi_ilp64 -lmkl_core 
BLASOPT=,blas
#---#endif mkl#
#
#---#if scalapack#
# include ScaLAPACK in the build
SPCKDIR=/software/software/ScaLAPACK/2.0.2-gompic-2.7.11-OpenBLAS-0.2.13-LAPACK-3.5.0/lib
SCALAPACK_LIB=${SPCKDIR}/libscalapack.a 
#---#endif scalapack#
#--#endif mkl_scalapack#

#--#if mpiwrap#
MPI_INCLUDE=
MPI_LIB=-L/software/software/OpenMPI/1.8.4-iccifort-2015.2.164/lib
LIBMPI = -lmpi -lmpi_ilp64 -lmpifort
# IC_INCLUDE is the include required
# to build against the interconnect
IC_INCLUDE=
LIB_EXTRA=
#--#else#
MPI_INCLUDE=/cm/shared/apps/openmpi/open64/64/1.10.1/include
MPI_LIB=/cm/shared/apps/openmpi/open64/64/1.10.1/lib64
LIBMPI =  -L$(MPI_LIB) -lmpi -lmpi_mpifh
#---#if myrinet#
# MPI and GM Locations for MYRINET switch
GM_INCLUDE=/opt/gm/2.0.19/include
IC_INCLUDE=-I${GM_INCLUDE}
GM_LIB=/opt/gm/2.0.19/lib64
LIB_EXTRA=-L$(GM_LIB) -lgm
#
#---#elseif infiniband openib#
#
# MPI Locations for INFINIBAND switch
IB_INCLUDE=/opt/ibgd/driver/infinihost/include
IC_INCLUDE=-I${IB_INCLUDE}
IB_LIB=/opt/ibgd/driver/infinihost/lib64
#LIB_EXTRA= -L${IB_LIB} -lvapi
#---#endif#
#--#endif mpiwrap#


# The following two GA-flags are actually used in other (non-GA) places
# as well. So just set the bastards.
GA_TARGET=LINUX64
GA_F77_DEFS = -traditional
# PEIGS definitions
PEIGS_TARGET_CPU_PAR=x86_64

#--#if ga#
# GA definitions & MPI and Switch library locations
# Change below to compile the GA's with a different c-compiler
#GA_CC = gcc
#
#---#if myrinet#
#
GA_VERSION_PAR= GA_VERSION=SHMEM ARMCI_NETWORK=GM GM_INCLUDE=${GM_INCLUDE} GM_LIB=${GM_LIB} MPI_INCLUDE=${MPI_INCLUDE} MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#
#---#elseif infiniband openib#
#
GA_VERSION_PAR= GA_VERSION=SHMEM ARMCI_NETWORK=OPENIB IB_INCLUDE=${IB_INCLUDE} MPI_INCLUDE=${MPI_INCLUDE} MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#
#---#else#
#
# Default sockets version
GA_VERSION_PAR= GA_VERSION=SHMEM MPI_INCLUDE=${MPI_INCLUDE} MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
LIB_EXTRA=-lpthread
#
#---#endif#
#--#endif ga#


#
# ===============  Compiler Options
#
#--#if mpiwrap#
FC = mpif90 -mcmodel=large -shared-intel
LD =  mpif90 -mcmodel=large -shared-intel
FC90 = mpif90 -mcmodel=large -shared-intel
#CC = mpicc -no-gcc
CC = mpicc -cc=icc -mcmodel=large -shared-intel
CXX= mpiCC
FINC =
#--#elseif score#
# Compiler locations for SCORE
FC = mpif90 -compiler intel9
LD = mpif90 -compiler intel9
FC90 = mpif90 -compiler intel9
CC= mpicc -compiler intel9
CXX= mpicc -compiler intel9
FINC =
LDSCORE=-notstatic
#--#else#
FC = ifort 
LD = ifort 
FC90 = ifort
CC=  icc -no-gcc
CXX = icc
FINC = -I${MPI_INCLUDE} ${IC_INCLUDE}
#--#endif#

#--#if i8#
FFLAGSI8 = -i8
CFLAGSI8 = -DLONG_INTEGER -DEXT_INT
# Added linux64 to get tst1s extracted from direct.m
# jmht - don't think this is actually needed, but have left it.
I8_M4_OPT=,i8,linux64
#--#else#
FFLAGSI8=
CFLAGSI8= -DSTD_INT
I8_M4_OPT=,i8drct,64bitpointers
#--#endif i8#

FFLAGSTMP = -c ${FFLAGSI8} ${FINC} 
# -std=c99 - might need this if int64_t can't be found
CFLAGSTMP = -c ${CFLAGSI8}  ${FINC} 

# Taken the load map generation out. With the intel 11.x compilers if the file
# load.map does not exist the loader skips to do something fancy to the next file on
# the line. Typically that means that mains.o gets corrupted after which the link
# stage fatally fails. If you touch load.map before attempting linking everything
# is fine!
#LDFLAGSTMP= -Wl,-Map load.map ${LDSCORE}
LDFLAGSTMP= ${LDSCORE} 

#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g
FFLAGSS  = ${FFLAGSTMP} -g
FFLAGSN  = ${FFLAGSTMP} -g
FFLAGSN0 =  ${FFLAGSTMP} -g 
CFLAGS   = ${CFLAGSTMP}  -g
LDFLAGS  = ${LDFLAGSTMP} -g 
#--#else# 
FFLAGSV  = ${FFLAGSTMP} -O3 -ip -ftz
FFLAGSS  = ${FFLAGSTMP} -O2 -ip -ftz
FFLAGSN  = ${FFLAGSTMP} -O1 -ip -ftz
FFLAGSN0 = ${FFLAGSTMP} -ip -ftz 
CFLAGS   = ${CFLAGSTMP}
LDFLAGS  = ${LDFLAGSTMP}
#--#endif#


# Put the link line together
BL_LIB = ${SCALAPACK_LIB} ${LIBBLAS} ${LIBMPI} ${LIB_EXTRA}

# M4 Options 
OPTIONS=${MACHOPT}${I8_M4_OPT}${BLASOPT}

#
# ===============  Additional Files 
#
#--#if i8#
# jmht - don't think this is actually needed, but have left it.
EXTRAI8= tst1s.o helfey.o gethes.o
EXTRA= $(EXTRAI8)
#--#else#
EXTRAI8=
EXTRA= helfey.o
#--#endif#

EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o aprq34d.o mcdab_ga.o
EXTRA_DFT = jkint_dft.o jkder_dft.o

#
# ===============  Compiler Exceptions
#
# jmht - don't think this is actually needed, but have left it.
tst1s.o:	direct.m
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f 

gethes.o:	casb.m
	cat ../utilities/gener.m casb.m | $(M4) -DGEN_EXTRACTFILE=gethes $(M4OPTS) > gethes.f
	$(FC) $(FFLAGSN) gethes.f

#
rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

helfey.o:       drv1e.m
	cat ../utilities/gener.m drv1e.m | $(M4) -DGEN_EXTRACTFILE=helfey $(M4OPTS) > helfey.f
	$(FC) $(FFLAGSN) helfey.f

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
#  ========== DFT Exceptions ===============
#
jkint_dft.o:	integ2e.m
	cat ../utilities/gener.m  integ2e.m | $(M4) -DGEN_EXTRACTFILE=jkint_dft $(M4OPTS) $(SNGL) > jkint_dft.f
	$(FC) $(FFLAGSN) jkint_dft.f
#
jkder_dft.o:	deriv2e.m
	cat ../utilities/gener.m  deriv2e.m | $(M4) -DGEN_EXTRACTFILE=jkder_dft $(M4OPTS) $(SNGL) > jkder_dft.f
	$(FC) $(FFLAGSN) jkder_dft.f

