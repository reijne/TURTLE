#doc This is the machine-specific file for the HECToR (http://www.hector.ac.uk)
#doc GLOBAL ARRAY based parallel version of GAMESS-UK on Linux running on 
#doc Opteron processors with PGI compilers (tested with 7.1.4.)
#doc
#doc Options:
#doc datain    - force GAMESS-UK to read it's input from a file called datain
#doc             instead of standard input
#doc score     - build on an SCORE system (WARNING - this build only works on
#doc             SCORE 5.8.4 & above)
#doc mpiwrap   - use mpi wrappers scripts to locate mpi libraries and headers
#doc newscf    - include distributed data scf module
#doc scalapack - include scalapack interface if the library is available
#doc blas      - include external blas libraries.
#doc static_lb - statically load-balanced version (no processor set aside as
#doc             nxtval server)
#
# Default options:
#dopt ga mpi ci peigs dl-find vdw masscf mpiwrap i8 datain
#opt blas scalapack static_lb vb zora newscf
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
MPI_INCLUDE  = /opt/mpich/PGI/include
IMPI_INCLUDE= -I${MPI_INCLUDE}
MPI_LIB      = /opt/mpich/PGI/lib
LIBMPI =  -L$(MPI_LIB) -lmpich -lfmpich
#--#endif mpiwrap#


#
# ===============  Compiler Options
#

FC = ftn
FC90 = ftn
LD = ftn
CC = cc

#Default compilation Flags

#--#if i8#
FFLAGSI8 = -i8
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
CFLAGSD = -c -tp=k8-64 ${CFLAGSI8}  ${IMPI_INCLUDE} ${IGM_INCLUDE}

#--#if static#
LDFLAGSD= -tp=k8-64 -Bstatic
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
BLASOPT=,blas
LIBBLAS=/opt/xt-libsci/10.4.6/pgi/lib/libsci.a
#---#endif i8#
#--#endif blas#

#--#if scalapack#
LIB_SCALAPACK= /apps/libs/ScaLAPACK/scalapack-1.8.0/lib/libscalapack.a /apps/libs/BLACS/1.1/LIB/blacsF77init_MPI-LINUX-0.a /apps/libs/BLACS/1.1/LIB/blacs_MPI-LINUX-0.a
#--#endif scalapack#


# PEIGS definitions
PEIGS_TARGET_CPU_PAR=x86_64

#--#if ga#
# Global Array definitions
GA_F77_DEFS = -traditional
GA_TARGET=LINUX64
# May need to use GCC for compiling the GAs
#GA_CC = gcc
# Default MPI/TCP/IP socket code
GA_VERSION_PAR=MA_USE_ARMCI_MEM=1 XT_SYMMETRIC_HEAP_SIZE=128M ARMCI_NETWORK=PORTALS FC=ftn CC=cc FOPT="-fastsse" COPT="-fastsse"
# Need to overwrite GA_MPI
#---#if i8#
GA_MPI=USE_MPI=yes
#---#else#
GA_MPI=USE_INTEGER4=yes USE_MPI=yes
#---#endif i8#
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

