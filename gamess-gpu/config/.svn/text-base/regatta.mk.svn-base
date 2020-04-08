#doc  This is the machine-specific file for the GLOBAL ARRAY based parallel 
#doc  version on a single Regatta node (does not use LAPI)
#doc  Integer*8 and 64 bit options are enabled
#doc
#doc Options:
#doc pSeries_mpi - triggers poe-based rungamess environment
#doc essl        - link against the essl libaray
#
#dopt ga tcgmsg-mpi mp2 peigs pSeries_mpi
#opt essl 
#
# ================ M4 Processing options
#
MACHINE_KEY=r
MACHOPT=rs6000,cio,unix,doublebackslash,upck-equiv

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
CPP = /usr/lib/cpp

CFLAGSI8 =  -DLONG_INTEGER
FFLAGSI8 = -qintsize=8
I8_M4_OPT = ,i8

# this build is based on -qEXTNAME
#####NOCU=-DNOC_
EXTNAME=-qEXTNAME

#--#if debug#
FFLAGSV = ${FFLAGSTMP} ${FFLAGSI8} -g
FFLAGSS = ${FFLAGSTMP} ${FFLAGSI8} -g
FFLAGSN = ${FFLAGSTMP} ${FFLAGSI8} -g
CFLAGS  = $(ARCH1) ${CFLAGSI8}  -O -c ${NOCU} ${FLAGS64} ${LARGEFILES}
#
#--#else# 
#
###FFLAGSV = -O3 -c $(ARCH2) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
#FFLAGSV = -O2 -c $(ARCH2) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
FFLAGSV = -O3 -qmaxmem=-1 -c $(ARCH2) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
FFLAGSS = -O -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGSN = -g -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGS0 =    -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
CFLAGS  = $(ARCH1) ${CFLAGSI8}  -O -c ${NOCU} ${FLAGS64} ${LARGEFILES}
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSV)
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSS)
#--#endif#

#--#if static#
echo static unimplemented
exit
#--#else#
LDFLAGS  = -q64 -b loadmap:load.map -b bigtoc
#--#endif#

#--#if essl#
lBLAS= -lessl
BLASOPT=,blas
#--#essl#

MPI_INCLUDE = .
MPI_LIB =
LIBMPI = 

#GA stuff
GA_F77_DEFS = -traditional
GA_TARGET=IBM64
# seems to be no longer needed...
GA_VERSION_PAR="USE_MPI=YES"

# problem with 64bitpointers
##OPTIONS=${MACHOPT}${BLASOPT},64bitpointers
OPTIONS=${MACHOPT}${BLASOPT}${I8_M4_OPT}

BL_LIB = ${LIBBLAS}
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
#
# DFT module
#
integ2e.o:	integ2e.m
	cat  ../utilities/gener.m integ2e.m | $(M4) $(M4OPTS) > integ2e.f
	$(FC) $(FFLAGSS) integ2e.f

xc.o:	xc.m
	cat  ../utilities/gener.m xc.m | $(M4) $(M4OPTS) > xc.f
	$(FC) $(FFLAGSS) xc.f

