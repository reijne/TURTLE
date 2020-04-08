#doc MPI AIX build on IBM Machines. This is a 64-bit build with integer*4 and without -qEXTNAME
#doc The build currently does not link PeIGS (but probably should)
#doc
#doc Options:
#doc pSeries_mpi - triggers poe-based rungamess environment
#doc datain - force GAMESS-UK to read input from a file called datain instead of standard input
#doc essl   - link against the ESSL numeric library
#doc dynamic_lb - build the dynamic loadbalancing MPI version (incompatible with newscf)
#doc newscf - build the new distributed data SCF/DFT driver (requires BLACS and ScaLAPACK)
#doc              newscf assumes scalapack lapack and blacs etc are in /usr/local/lib
#
#
#dopt base pSeries_mpi datain
#opt blas nbo drf zora newscf debug dynamic_lb vampir
#
#
# ================ M4 Processing options
#
MACHINE_KEY=r
MACHOPT=rs6000,cio,unix,doublebackslash,upck-equiv,rs6000_noextname

## Machine specific options from: rs6000.m ##
IAND64 = iand($$1,$$2)
IOR64 = ior($$1,$$2)
IXOR64 = ieor($$1,$$2)
EXTRA_M4_DEFINITIONS = -Dleadz='leadzz($$1)'
## end of rs600.m  ##

#
# ===============  Compiler Options
#
FC = mpxlf
LD = mpxlf
FC90 = mpxlf90
CC= mpcc

FFLAGSTMP = -c
RANLIB=ranlib
ARCHIVE  =  ar -X 64 rcv
ARCH1=-qarch=auto
ARCH2=-qarch=auto -qtune=auto
LARGEFILES = -D_LARGE_FILES
FLAGS64 = -q64
NOCU=-DNOC_
EXTNAME=
#--#if debug#
FFLAGSV = ${FFLAGSTMP} $(EXTNAME) ${FLAGS64} -g -qxlf77=leadzero
FFLAGSS = ${FFLAGSTMP} $(EXTNAME) ${FLAGS64} -g -qxlf77=leadzero
FFLAGSN = ${FFLAGSTMP} $(EXTNAME) ${FLAGS64}-g -qxlf77=leadzero
FFLAGS0 = ${FFLAGSTMP} $(EXTNAME) -g ${FLAGS64} -qxlf77=leadzero
CFLAGS  = $(ARCH1) ${CFLAGSI8}  -O -c ${NOCU} ${FLAGS64} ${LARGEFILES}
FFLAGSV90 = -C -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSV)
FFLAGSS90 = -C -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSS)
#--#else# 
###FFLAGSV = -O3 -c $(ARCH2) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
#FFLAGSV = -O2 -c $(ARCH2) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
FFLAGSV = -O3 -qmaxmem=-1 -c $(ARCH2) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64} -qxlf77=leadzero
FFLAGSS = -O -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64} -qxlf77=leadzero
FFLAGSN = -g -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64} -qxlf77=leadzero
FFLAGS0 =    -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64} -qxlf77=leadzero
CFLAGS  = $(ARCH1) ${CFLAGSI8}  -O -c ${NOCU} ${FLAGS64} ${LARGEFILES}
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSV)
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSS)
#--#endif#

LDFLAGS  = -q64 -b loadmap:load.map -b bigtoc

#--#if essl#
LIBBLAS= -lessl
BLASOPT=,blas
#--#endif essl#

#--#if vampir#
LIBVT = -L/hpcx/usrlocal/packages/vampir/lib64 -lVT
VTOPT=,vampir
#--#endif#

# no need to set the mpi stuff if we use mpxlf etc
MPI_INCLUDE = .
MPI_LIB =
LIBMPI = 

# seems to be a problem with 64bitpointers
###OPTIONS=${MACHOPT}${BLASOPT},64bitpointers
OPTIONS=${MACHOPT}${BLASOPT}

#--#if newscf#
PAR_M4_OPTS=parallel,mpi${VTOPT}
BL_LIB= -L /usr/local/lib -lblacsF77init -lblacs -lscalapack -ltools -llapack -lessl ${LIBVT}
#--#elseif dynamic_lb#
PAR_M4_OPTS=parallel,mpi,dynamic${VTOPT}
BL_LIB = ${LIBBLAS} ${LIBVT}
#--#else#
PAR_M4_OPTS=parallel,mpi${VTOPT}
BL_LIB = ${LIBBLAS} ${LIBVT}
#--#endif#


# no need to set with mpxlf
MPI_LIBS = 

#
# ===============  Additional Files 
#
EXTRA_BASE= iterate.o
EXTRA=
EXTRA_MP2= 
#
#  ========== Exceptions for IBM r6000 power1, 2 and 3 (and SP2) ===============
#
iterate.o: morokuma.m
	cat ../utilities/gener.m  morokuma.m | $(M4) -DGEN_EXTRACTFILE=iterate $(M4OPTS)  > iterate.f
	$(FC) $(FFLAGS0) iterate.f

# DFT module
integ2e.o:	integ2e.m
	cat  ../utilities/gener.m integ2e.m | $(M4) $(M4OPTS) > integ2e.f
	$(FC) $(FFLAGSS) integ2e.f

xc.o:	xc.m
	cat  ../utilities/gener.m xc.m | $(M4) $(M4OPTS) > xc.f
	$(FC) $(FFLAGSS) xc.f
