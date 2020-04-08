#  Machine-dependent makefile settings for serial AIX build 
#  on IBM Machines. This is a 64-bit build with integer*4.
#  This file was created on the BlueGene machine at EPCC.
#
#doc For IBM AIX Power440 using XLF compilers
#doc blas option uses IBM ESSL

#dopt base mpi datain blas newscf qmmm
#opt nbo drf zora debug dynamic_lb
#
#
# ================ M4 Processing options
#
MACHINE_KEY=r
MACHOPT=bluegene,rs6000,cio,unix,doublebackslash,upck-equiv,rs6000_noextname

## Machine specific options from: rs6000.m ##
IAND64 = iand($$1,$$2)
IOR64 = ior($$1,$$2)
IXOR64 = ieor($$1,$$2)
EXTRA_M4_DEFINITIONS = -Dleadz='leadzz($$1)'
## end of rs6000.m  ##

#
# ===============  Compiler Options
#
#FC = blrts_xlf
#LD = blrts_xlf90
#FC90 = blrts_xlf90
#CC= blrts_xlc
FC = mpixlf
LD = mpixlf90
FC90 = mpixlf90
CC= mpixlc

FFLAGSTMP = -c
RANLIB=ranlib
ARCHIVE  =  ar rcv
ARCH1=-qarch=440d
ARCH1=-qarch=440d -qflttrap=overflow:zerodivide:invalid
ARCH2=-qarch=440d -qtune=440
LARGEFILES = -D_LARGE_FILES
#EXTNAME=-qEXTNAME
EXTNAME=
#LDFLAGS0  = -b loadmap:load.map -b bigtoc
LDFLAGS0  = 
####NOCU=-DNOC_

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

#FFLAGSV = -g -O3 -c -qnosave $(ARCH1) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
#FFLAGSS = -g -O  -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
#FFLAGS1 = -g -O  -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
#FFLAGSN = -g  -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
#FFLAGS0 = -g -O0 -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
#CFLAGS  = $(ARCH1) ${CFLAGSI8}  -O -c ${NOCU} ${FLAGS64} ${LARGEFILES}
#
# First go back to basics making sure the code is "right"
#
FFLAGSV = -g -O0 -c -qnosave $(ARCH1) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
FFLAGSS = -g -O0 -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGS1 = -g -O0 -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGSN = -g -O0 -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGS0 = -g -O0 -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
CFLAGS  = $(ARCH1) ${CFLAGSI8} -O0 -c ${NOCU} ${FLAGS64} ${LARGEFILES}
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSV)
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSS)

#--#endif#

#--#if blas#
LIBBLAS= -L/bgl/local/lib -lessln
BLASOPT=,blas
#--#endif#

# MPI  locations
#MPI_LIB      = /bgl/BlueLight/ppcfloor/bglsys/lib
#MPI_INCLUDE  = /bgl/BlueLight/ppcfloor/bglsys/include
#LIBMPI =  -L$(MPI_LIB) -lmpich -lfmpich


LDFLAGS = $(LDFLAGS0) $(FLAGS64) $(RENAMES) 

###OPTIONS=${MACHOPT}${BLASOPT},64bitpointers
OPTIONS=${MACHOPT}${BLASOPT}

#--#if newscf#
LIBSCALAPACK = -L/bgl/local/lib -lscalapack -lblacsF77init_MPI-BGENE-1 -lblacs_MPI-BGENE-1 -llapack440
#--#endif#
NSSLIB = -lc -lnss_files -lnss_dns -lresolv
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
#  ==================  Exceptions for IBM r6000 440d  ====================
#
intega.o:	intega.f
	$(FC) $(FFLAGS0) $*.f
guess.o:	guess.f
	$(FC) $(FFLAGSS) $*.f
analb.o:	analb.f
	$(FC) $(FFLAGSS) $*.f
util1.o:	util1.f
	$(FC) $(FFLAGSS) $*.f
master.o:	master.f
	$(FC) $(FFLAGSS) $*.f
analc.o:	analc.f
	$(FC) $(FFLAGSS) $*.f
server.o:	server.f
	$(FC) $(FFLAGSS) $*.f
direct.o:	direct.f
	$(FC) $(FFLAGSS) $*.f
sec2e.o:	sec2e.f
	$(FC) $(FFLAGSS) $*.f
util6.o:	util6.f
	$(FC) $(FFLAGSS) $*.f
dirrpa.o:	dirrpa.f
	$(FC) $(FFLAGSS) $*.f
util7.o:	util7.f
	$(FC) $(FFLAGSS) $*.f
cphf.o:	cphf.f
	$(FC) $(FFLAGSS) $*.f
secmp2.o:	secmp2.f
	$(FC) $(FFLAGS1) $*.f
#
#  ========== Exceptions for IBM r6000 power1,2,3 and 4 (and SP2) ===============
#
integ2e.o:	integ2e.m
	cat ../utilities/gener.m integ2e.m | $(M4) $(M4OPTS) > integ2e.f
	$(FC) $(FFLAGSS) integ2e.f

xc.o:	xc.m
	cat ../utilities/gener.m xc.m | $(M4) $(M4OPTS) > xc.f
	$(FC) $(FFLAGSS) xc.f

