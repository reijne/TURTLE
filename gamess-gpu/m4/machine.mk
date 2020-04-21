#doc This is the machine-specific file for the GLOBAL ARRAY based parallel version
#doc of GAMESS-UK on Linux running on Opteron processors with PGI compilers (tested with 5.2.)
#doc
#doc Options:
#doc datain - force GAMESS-UK to read it's input from a file called datain instead of standard input
#doc score  - build on an SCORE system (WARNING - this build only works on SCORE 5.8.4 & above)
#doc mpiwrap - use mpi wrappers scripts to locate mpi libraries and headers
#doc newscf - include distributed data scf module
#doc scalapack - include scalapack interface if the library is available
#doc blas   - inlclude external blas libraries.
#doc myrinet - build against MX drivers on a Myrinet switch (GM not supported yet)
#doc static_lb - statically load-balanced version (no processor set aside as nxtval server)
#
# Default options:
#dopt ga mpi ci peigs newscf dl-find vdw masscf mpiwrap
#opt score blas myrinet datain i8 scalapack static_lb vb zora
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=opteron,linux,littleendian,cio,unix,doublebackslash,upck-equiv,extpopcnt,extleadz

# Bitwise operators that differ from the defaults.
IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)
#
#
RANLIB=ranlib
CPP=/lib/cpp

#--#if mpiwrap#
MPI_INCLUDE  =
IMPI_INCLUDE  =
MPI_LIB      = 
LIBMPI       = 
#--#else#
MPI_INCLUDE  = /opt/pgi/linux86-64/2019/mpi/openmpi-3.1.3/include/
IMPI_INCLUDE= -I${MPI_INCLUDE}
MPI_LIB      = /opt/pgi/linux86-64/2019/mpi/openmpi-3.1.3/lib/
LIBMPI =  -L$(MPI_LIB) -lmpi -lmpi_mpifh
#--#endif mpiwrap#

#--#if myrinet#
# Building against GM for a Myrinet switch
GM_INCLUDE = /opt/gm/include
IGM_INCLUDE= -I${GM_INCLUDE}
GM_LIB = /opt/gm/lib64
LIBGM =  -L$(GM_LIB) -lgm -lpthread
#--#else#
GM_INCLUDE =
IGM_INCLUDE =
GM_LIB =
LIBGM =
#--#endif myrinet#

#
# ===============  Compiler Options
#

#--#if score#
FC = mpif90 -compiler pgi
FC90 = mpif90 -compiler pgi
LD = mpif90 -compiler pgi
CC = mpicc -compiler pgi
#--#elseif mpiwrap#
FC =  mpif90
FC90 =  mpif90
LD =  mpif90
CC =  mpicc
#--#else#
FC = pgf90
FC90 = pgf90
LD = pgf90
CC = pgcc
#--#endif#

#Default compilation Flags

#--#if i8#
#---#if score#
# The Mlfs is a hack discovered by Cliff Addisson at Liverpool in
# his attempts to get VASP to work with 64bit/i8 and SCORE
# also add -Mdalign as it is used by the Global Arrays
# (which were also hacked to use -Mlfs)
FFLAGSI8 = -i8 -Mlfs -Mdalign
#---#else#
FFLAGSI8 = -i8
#---#endif score#
CFLAGSI8 = -DLONG_INTEGER -DEXT_INT
I8_M4_OPT=,i8
#--#else#
FFLAGSI8 = 
CFLAGSI8 = -DSTD_INT
#---#if !ga#
# If we're not running with GA (& therefore ma) we need to 
# use 64bitpointers
I8_M4_OPT=,64bitpointers
#---#else#
I8_M4_OPT=
#---#endif#
#--#endif#

FFLAGSD = -c -tp=k8-64 -Mcache_align ${FFLAGSI8} ${IMPI_INCLUDE} ${IGM_INCLUDE}
CFLAGSD= -c -tp=k8-64 ${CFLAGSI8}  ${IMPI_INCLUDE} ${IGM_INCLUDE}

#--#if static#
LDFLAGSD= -tp=k8-64 -Bstatic
#--#elseif score#
LDFLAGSD= -tp=k8-64 -no-static
#--#else#
LDFLAGSD= -tp=k8-64 
#--#endif static#
#

#--#if debug#
FFLAGSV = ${FFLAGSD} -g
FFLAGSS = ${FFLAGSD} -g
FFLAGSN = ${FFLAGSD} -g
FFLAGSN0 = ${FFLAGSD} -g
CFLAGS = ${CFLAGSD} -g
LDFLAGS  = ${LDFLAGSD} -g 
#--#else#
FFLAGSV =  ${FFLAGSD}  -O2
FFLAGSS =  ${FFLAGSD}  -O
FFLAGSN =  ${FFLAGSD}  -O1
FFLAGSN0 = ${FFLAGSD}
CFLAGS  = -D_REENTRANT ${CFLAGSD}
LDFLAGS = ${LDFLAGSD}
#--#endif debug#
#
#

#--#if blas#
#---#if i8#
#
BLAS NOT SUPPORTED FOR I8 build!!
#
#---#else#
#----#if score#
BLASOPT=,blas
# Below is o.k. but sources lapack from the pgi directory
LIBBLAS= /usr/local/lib64/liblapack_common.a /usr/local/lib64/liblapack_pgi.a /usr/local/lib64/liblapack_common.a -lg2c /usr/local/lib64/libblas.so
#----#else#
BLASOPT=,blas
LIBBLAS=-L/apps/libs/acml/3.6.0/pgi-64bit/pgi64/lib -lacml
#----#endif score#
#---#endif i8#
#--#endif blas#

#--#if scalapack#
#---#if score#
#LIBF90=-L/usr/local/lib64  -lscalapack  -lblacsF77init -lblacs
LIB_SCALAPACK= /usr/local/lib64/libscalapack.a  /usr/local/lib64/libblacsF77init.a /usr/local/lib64/libblacs.a
#---#else#
#LIBF90= -L/usr/local/lib/SCALAPACK -lscalapack -L/usr/local/lib/BLACS/LIB -lblacsF77init_MPI-LINUX64-1 -lblacs_MPI-LINUX64-1
LIB_SCALAPACK= /apps/libs/ScaLAPACK/scalapack-1.8.0/lib/libscalapack.a /apps/libs/BLACS/1.1/LIB/blacsF77init_MPI-LINUX-0.a /apps/libs/BLACS/1.1/LIB/blacs_MPI-LINUX-0.a
#---#endif score#
#--#endif scalapack#


# PEIGS definitions
PEIGS_TARGET_CPU_PAR=x86_64

#--#if ga#
# Global Array definitions
GA_F77_DEFS = -traditional
GA_TARGET=LINUX64
# May need to use GCC for compiling the GAs
#GA_CC = gcc
#---#if myrinet#
GA_VERSION_PAR= GA_VERSION=SHMEM ARMCI_NETWORK=GM GM_INCLUDE=${GM_INCLUDE} GM_LIB=${GM_LIB} MPI_INCLUDE=${MPI_INCLUDE} MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI}"
#---#else#
# Default MPI/TCP/IP socket code
GA_VERSION_PAR= GA_VERSION=SHMEM MPI_INCLUDE=${MPI_INCLUDE} MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#---#endif myrinet#
#--#else#
#--#endif ga#


# Bring all the options together
OPTIONS=${MACHOPT}${BLASOPT}${I8_M4_OPT}

BL_LIB = ${LIB_SCALAPACK} ${LIBBLAS} ${LIBMPI} ${LIBGM}

#
# ===============  Additional Files
#
EXTRA=
EXTRA_MP2= check0a.o aprq34d.o mcdab_ga.o
EXTRA_DFT = jkint_dft.o jkder_dft.o
#
# ===============  Compiler Exceptions
#

check0a.o: newmrd5.m
	cat ../utilities/gener.m newmrd5.m | $(M4) -DGEN_EXTRACTFILE=check0a $(M4OPTS) > check0a.f
	$(FC) $(FFLAGSN) check0a.f

# The below are extracted by the linux keyword and are probably not required for this compiler
#
aprq34d.o: mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=aprq34d $(M4OPTS) > aprq34d.f
	$(FC) $(FFLAGSS) aprq34d.f

mcdab_ga.o: mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=mcdab_ga $(M4OPTS) > mcdab_ga.f
	$(FC) $(FFLAGSS) mcdab_ga.f

#
#  ========== DFT Exceptions (PGF) ===============
#
# The below are extracted by the opteron adn em64t keywords and are probably not required for this compiler
#
jkint_dft.o:    integ2e.m
	cat ../utilities/gener.m  integ2e.m | $(M4) -DGEN_EXTRACTFILE=jkint_dft $(M4OPTS) $(SNGL) > jkint_dft.f
	$(FC) $(FFLAGSN) jkint_dft.f
#
jkder_dft.o: deriv2e.m
	cat ../utilities/gener.m  deriv2e.m | $(M4) -DGEN_EXTRACTFILE=jkder_dft $(M4OPTS) $(SNGL) > jkder_dft.f
	$(FC) $(FFLAGSN) jkder_dft.f

