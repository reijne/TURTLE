#doc This is the machine-specific file for the GLOBAL ARRAY based parallel version on HUYGENS
#doc derived from the one on HPCx
#doc it includes i8 and 64bit options for compatibility with GA
#doc currently  tcgmsg are  bad idea
#doc nbo is serial and is there for copatibility
#doc
#doc Keywords:
#doc edgafs - put the ed files in a Global Array by default (added following problems on Huygens)
#doc
#
#dopt ga mpi i8 ci peigs vb zora dl-find vdw newscf nbo
#opt essl mp3 edgafs
#
# ================ M4 Processing options
#
MACHINE_KEY=r
MACHOPT=rs6000,cio,unix,doublebackslash,upck-equiv,i8

## Machine specific options from: rs6000.m ##
IAND64 = iand($$1,$$2)
IOR64 = ior($$1,$$2)
IXOR64 = ieor($$1,$$2)
EXTRA_M4_DEFINITIONS = -Dleadz='leadzz($$1)'
## end of rs600.m  ##

#
# ===============  Compiler Options
#

FC = mpfort -compiler xlf90_r -qfixed
LD = mpfort -compiler xlf90_r
FC90 = mpfort -compiler xlf90_r
CC= mpcc

FFLAGSTMP = -c
RANLIB=ranlib
ARCHIVE  =  ar rcv
ARCH1=-qarch=auto
ARCH2=-qarch=auto -qtune=auto
LARGEFILES = -D_LARGE_FILES
CFLAGS64 = -q64
FLAGS64 = -q64
CPP = /usr/bin/cpp
CFLAGSI8 =  -DLONG_INTEGER -DEXT_INT
FFLAGSI8 = -qintsize=8
I8_M4_OPTS = ,i8

# this build is based on -qEXTNAME
#####NOCU=-DNOC_
EXTNAME=-qEXTNAME

#--#if debug#
FFLAGSV = ${FFLAGSTMP} ${FFLAGSI8} ${EXTNAME} -g
FFLAGSS = ${FFLAGSTMP} ${FFLAGSI8} ${EXTNAME} -g
FFLAGSN = ${FFLAGSTMP} ${FFLAGSI8} ${EXTNAME} -g
CFLAGS  = $(ARCH1) ${CFLAGSI8}  -O -c ${NOCU} ${CFLAGS64} ${LARGEFILES}
#
#--#else# 
#
###FFLAGSV = -O3 -c $(ARCH2) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
#FFLAGSV = -O2 -c $(ARCH2) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
#FFLAGSV = -O3 -qmaxmem=-1 -c $(ARCH2) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
FFLAGSV = -O2 -qmaxmem=-1 -c $(ARCH2) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
FFLAGSS = -O -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGSN = -g -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGS0 =    -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
CFLAGS  = ${CFLAGSI8}  -O -c ${NOCU} ${CFLAGS64} ${LARGEFILES}
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 -qfree=f90 $(FFLAGSV)
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 -qfree=f90 $(FFLAGSS)
#--#endif#

#--#if static#
echo static unimplemented
exit
#--#else#
#####-b loadmap:load.map -b bigtoc
LDFLAGS  = -q64 
#--#endif#

#--#if i8#
#--#else#
echo i4 not provided
#--#endif#

#--#if essl#
LIBBLAS= -lessl6464
BLASOPT=,blas
#--#endif essl#

# MPI variables
MPI_INCLUDE = .
MPI_LIB =
LIBMPI = 

#Global Array stuff
GA_F77_DEFS = -traditional
GA_TARGET=LINUX64
# seems to be no longer needed...
GA_VERSION_PAR=

# problem with 64bitpointers
##OPTIONS=${MACHOPT}${BLASOPT},64bitpointers
#--#if edgafs#
# Put ed files in GA by default
EDGAFSM4=,edgafs
#--#endif edgafs#

OPTIONS=${MACHOPT}${BLASOPT}${EDGAFSM4}

BL_LIB = ${LIBBLAS}
#MPI_LIBS = -L${MPI_LIB} ${LIBMPI}
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

