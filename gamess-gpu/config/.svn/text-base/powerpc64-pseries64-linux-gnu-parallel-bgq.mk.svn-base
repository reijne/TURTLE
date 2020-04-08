#  Machine-dependent makefile settings for serial linux build 
#  on IBM BlueGene/Q machines. This is a 64-bit build with integer*4.
# 
#  This file was created on Blue Joule at Daresbury Laboratory. 
#  TWK 31 Jan 2013
#
#doc For IBM BlueGene/Q using XLF compilers
#doc blas option uses IBM ESSL

#dopt base mpi datain blas newscf qmmm scalapack
#opt nbo drf zora debug dynamic_lb
#
#
# ================ M4 Processing options
#
MACHINE_KEY=r
MACHOPT=bluegene,rs6000,cio,unix,doublebackslash,upck-equiv,64bitpointers

## Machine specific options from: rs6000.m ##
IAND64 = iand($$1,$$2)
IOR64 = ior($$1,$$2)
IXOR64 = ieor($$1,$$2)
EXTRA_M4_DEFINITIONS = -Dleadz='leadzz($$1)'
## end of rs6000.m  ##

#
# ===============  Compiler Options
#
FC = /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlf77_r
LD = /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlf90_r
FC90 = /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlf90_r
CC=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlc_r

FFLAGSTMP = -c
RANLIB=ranlib
ARCHIVE  =  ar rcv
ARCH1=-qarch=qp -qtune=qp
LARGEFILES = -D_LARGE_FILES
EXTNAME=-qEXTNAME
LDFLAGS0  = 

#--#if 64bit#
ARCHIVE = ar -X 64 rcv
FLAGS64 = -q64
#--#endif#

#--#if debug#
FFLAGSV = ${FFLAGSTMP} $(EXTNAME) -g
FFLAGSS = ${FFLAGSTMP} $(EXTNAME) -g
FFLAGSN = ${FFLAGSTMP} $(EXTNAME) -g
CFLAGS  = $(ARCH1) ${CFLAGSI8}  -O -c ${NOCU} ${FLAGS64} ${LARGEFILES}
#--#else# 

FFLAGSV = -g -O3 -c -qnosave $(ARCH1) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
FFLAGSS = -g -O  -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGS1 = -g -O  -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGSN = -g  -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGS0 = -g -O0 -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
CFLAGS  = $(ARCH1) ${CFLAGSI8}  -O -c ${NOCU} ${FLAGS64} ${LARGEFILES}

#--#endif#

#--#if blas#
LIBBLAS= -L/bgsys/ibm_essl/prod/opt/ibmmath/essl/5.1/lib64 -lesslbg
BLASOPT=,blas
#--#endif#

# MPI  locations
#MPI_LIB      =
#MPI_INCLUDE  = 
#LIBMPI =

LDFLAGS = $(LDFLAGS0) $(FLAGS64) $(RENAMES) 

OPTIONS=${MACHOPT}${BLASOPT}

#--#if scalapack#
LIBSCALAPACK = -L/gpfs/packages/ibm/scalapack/2.0.2/lib -lscalapack_extname
#--#endif#
NSSLIB=
BL_LIB = ${LIBSCALAPACK} ${LIBBLAS} ${LIBMPI} ${NSSLIB}
#
# ===============  Additional Files 
#
EXTRA_BASE= iterate.o
EXTRA=
EXTRA_MP2= 
#
#  ========== Exceptions for IBM r6000 power series ===============
# 
iterate.o: morokuma.m
	cat ../utilities/gener.m  morokuma.m | $(M4) -DGEN_EXTRACTFILE=iterate $(M4OPTS)  > iterate.f
	$(FC) $(FFLAGS0) iterate.f
#
