#  Machine-dependent makefile settings for parallel AIX GA build 
#  on IBM BlueGene Machines. This is a 32-bit build with integer*4.
#
#  This file was created on the BlueGene machine at EPCC.
#
#doc For IBM AIX Power440 using XLF compilers
#doc the Fortran is compiled with -qextname this means that if ScaLAPACK
#doc is used it must be compiled as follows:
#doc BLAS: edit make.inc and set
#doc       FORTRAN        = blrts_xlf
#doc       OPTS           = -qarch=440d -qtune=440 -qextname -O3 -qhot
#doc       NOOPT          = -qarch=440d -qtune=440 -qextname
#doc LAPACK: edit make.inc and set
#doc       FORTRAN        = blrts_xlf
#doc       OPTS           = -qarch=440d -qtune=440 -qextname -O3 -qstrict
#doc       NOOPT          = -qarch=440d -qtune=440 -qextname -qnoopt
#doc       TIMER          = EXT_ETIME
#doc BLACS: edit Bmake.inc and set
#doc       COMMLIB        = MPI
#doc       MPILIB         = -L$(MPILIBdir) -lfmpich_.rts -lmpich.rts -lfmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts
#doc       INTFACE        = -DAdd_
#doc       F77            = blrts_xlf
#doc       F77NO_OPTFLAGS = -qarch=440d -qtune=440 -qextname -qnoopt
#doc       F77FLAGS       = -qarch=440d -qtune=440 -qextname -O3 -qstrict
#doc       CC             = powerpc-bgl-blrts-gnu-gcc
#doc       CCFLAGS        = -O3
#doc       ARCH           = powerpc-bgl-blrts-gnu-ar
#doc       RANLIB         = powerpc-bgl-blrts-gnu-ranlib
#doc ScaLAPACK: edit SLmake.inc and set
#doc       USEMPI         = -DUsingMpiBlacs
#doc       SMPLIB         = -L/bgl/BlueLight/ppcfloor/bglsys/lib -lfmpich_.rts -lmpich.rts -lfmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts
#doc       F77            = blrts_xlf
#doc       CC             = powerpc-bgl-blrts-gnu-gcc
#doc       NOOPT          = -qarch=440d -qtune=440 -qextname -qnoopt
#doc       F77FLAGS       = -qarch=440d -qtune=440 -qextname -O3 -qstrict
#doc       CCFLAGS        = -O3
#doc       CDEFS          = -DAdd_ -DNO_IEEE $(USEMPI)
#doc       ARCH           = powerpc-bgl-blrts-gnu-ar
#doc       RANLIB         = powerpc-bgl-blrts-gnu-ranlib
#doc To build LAPACK you will also have to change the Makefile and
#doc TESTING/Makefile and insert mpirun -np 1 to run the tests. 
#doc The LAPACK build will fail if the tests do.


#dopt ga mpi peigs mp2 datain blas newscf
#opt nbo drf zora debug dynamic_lb
#
#
# ================ M4 Processing options
#
MACHINE_KEY=r
MACHOPT=bluegene,rs6000,cio,unix,doublebackslash,upck-equiv,idamin,idmax

## Machine specific options from: rs6000.m ##
IAND32 = iand($$1,$$2)
IOR32 = ior($$1,$$2)
IXOR32 = ieor($$1,$$2)
IAND64 = iand($$1,$$2)
IOR64 = ior($$1,$$2)
IXOR64 = ieor($$1,$$2)
EXTRA_M4_DEFINITIONS = -Dleadz='leadzz($$1)'
## end of rs6000.m  ##

#
# ===============  Compiler Options
#
#--#if mpiwrap#
FC      = mpixlf
LD      = mpixlf90
FC90    = mpixlf90
CC      = mpixlc
CPP     = cpp -traditional
RANLIB  = ranlib
ARCHIVE = ar rcv
#--#else#
#FC      = /opt/ibmcmp/xlf/bg/10.1/bin/blrts_xlf
#LD      = /opt/ibmcmp/xlf/bg/10.1/bin/blrts_xlf
FC      = /opt/ibmcmp/xlf/bg/10.1/bin/blrts_xlf90
LD      = /opt/ibmcmp/xlf/bg/10.1/bin/blrts_xlf90
FC90    = /opt/ibmcmp/xlf/bg/10.1/bin/blrts_xlf90
CC      = /bgl/BlueLight/ppcfloor/blrts-gnu/bin/powerpc-bgl-blrts-gnu-gcc 
CPP     = /bgl/BlueLight/ppcfloor/blrts-gnu/bin/powerpc-bgl-blrts-gnu-cpp -traditional
RANLIB  = /bgl/BlueLight/ppcfloor/blrts-gnu/bin/powerpc-bgl-blrts-gnu-ranlib
ARCHIVE = /bgl/BlueLight/ppcfloor/blrts-gnu/bin/powerpc-bgl-blrts-gnu-ar rcv
#--#endif#

GA_TARGET = BGL
PEIGS_TARGET = BGL
ARCH1=-qarch=440d -qtune=440 -qflttrap=overflow:zerodivide:invalid
ARCH2=-qarch=440d -qtune=440 -qflttrap=overflow:zerodivide:invalid
LARGEFILES = -D_LARGE_FILES
GA_F77_DEFS=-DBGL -DBLRTS -DBGML
#--#if scalapack#
# implies using mpi directly and => i4
GA_VERSION_PAR=ARMCI_NETWORK=BGMLMPI USE_SCALAPACK=yes
GA_MPI=MSG_COMMS=BGMLMPI USE_INTEGER4=yes
#--#else#
# no scalapack => tcgmsg-mpi & i8
GA_VERSION_PAR=ARMCI_NETWORK=BGMLMPI 
GA_MPI=MSG_COMMS=BGMLMPI USE_MPI=yes
#--#endif#
EXTNAME=-qextname

FFLAGSTMP = -c 

#--#if 64bit#
ARCHIVE = ar -X 64 rcv
FLAGS64 = -q64
#--#endif#

CFLAGSI8 = -DSTD_INT
FFLAGSI8 =
I8_M4_OPTS =
# CFLAGSI8 =  -DLONG_INTEGER -DEXT_INT
# FFLAGSI8 = -qintsize=8
# I8_M4_OPTS = ,i8

# MPI  locations
#--#if mpiwrap#
MPI_LIBS =
#--#else#
MPI_LIB     = /bgl/BlueLight/ppcfloor/bglsys/lib
MPI_INCLUDE = /bgl/BlueLight/ppcfloor/bglsys/include
LIBMPI      = -L${MPI_LIB} -lfmpich_.rts -lmpich.rts -lfmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts
MPI_LIBS    = ${LIBMPI}
FFLAGSMPI   = -I${MPI_INCLUDE}
#--#endif#

#--#if debug#
FFLAGSV  = ${FFLAGSTMP} $(EXTNAME) -g -O0 ${ARCH1} ${FFLAGSMPI} -qfixed
FFLAGSS  = ${FFLAGSTMP} $(EXTNAME) -g -O0 ${ARCH1} ${FFLAGSMPI} -qfixed
FFLAGSN  = ${FFLAGSTMP} $(EXTNAME) -g -O0 ${ARCH1} ${FFLAGSMPI} -qfixed
FFLAGS0  = ${FFLAGSTMP} $(EXTNAME) -g -O0 ${ARCH1} ${FFLAGSMPI} -qfixed
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 ${FFLAGSTMP} $(EXTNAME) -g -O0 ${ARCH1} ${FFLAGSMPI}
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 ${FFLAGSTMP} $(EXTNAME) -g -O0 ${ARCH1} ${FFLAGSMPI}
CFLAGS   = ${CFLAGSI8}  -g -O0 -c ${NOCU} ${FLAGS64} ${LARGEFILES}
LDFLAGS0 = -g 
#--#else# 
FFLAGSV = -O3 -c $(ARCH1) ${FFLAGSTMP} ${FFLAGSI8} ${EXTNAME} ${FLAGS64} ${FFLAGSMPI} -qfixed
FFLAGSS = -O3 -c $(ARCH1) ${FFLAGSTMP} ${FFLAGSI8} ${EXTNAME} ${FLAGS64} ${FFLAGSMPI} -qfixed
FFLAGS1 = -O3 -c $(ARCH1) ${FFLAGSTMP} ${FFLAGSI8} ${EXTNAME} ${FLAGS64} ${FFLAGSMPI} -qfixed
FFLAGSN = -O3 -c $(ARCH1) ${FFLAGSTMP} ${FFLAGSI8} ${EXTNAME} ${FLAGS64} ${FFLAGSMPI} -qfixed
FFLAGS0 = -O3 -c $(ARCH1) ${FFLAGSTMP} ${FFLAGSI8} ${EXTNAME} ${FLAGS64} ${FFLAGSMPI} -qfixed
CFLAGS  = $(ARCH1) ${CFLAGSI8} -O3 -c ${NOCU} ${FLAGS64} ${LARGEFILES}
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 ${FFLAGSTMP} $(EXTNAME) -O3 ${ARCH1} ${FFLAGSMPI}
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 ${FFLAGSTMP} $(EXTNAME) -O3 ${ARCH1} ${FFLAGSMPI}
LDFLAGS0  = 
#--#endif#

#--#if scalapack#
# This implies we are picking up ScaLAPACK and BLAS from somewhere (if only
# the reference implementation included in ScaLAPACK for the latter).
LIBSCALAPACK = -L${HOME}/ScaLAPACK/lib -lscalapack_BGL.440 -lblacsF77init_MPI-BGL.440-0 -lblacs_MPI-BGL.440-0 -llapack_BGL.440
LIBBLAS= -lblas_BGL.440
BLASOPT=,blas,scalapack
#--#elseif blas#
LIBBLAS= -lblas_BGL.440
BLASOPT=,blas
#--#endif#

# MPI  locations
#--#if mpiwrap#
MPI_LIBS =
#--#else#
MPI_LIB     = /bgl/BlueLight/ppcfloor/bglsys/lib
MPI_INCLUDE = /bgl/BlueLight/ppcfloor/bglsys/include
LIBMPI      = -L${MPI_LIB} -lmpich.rts -lfmpich.rts -lmsglayer.rts -ldevices.rts -lrts.rts
MPI_LIBS    = ${LIBMPI}
#--#endif#

LDFLAGS = $(LDFLAGS0) $(FLAGS64) $(RENAMES) 

OPTIONS=${MACHOPT}${BLASOPT}

#--#if newscf#
LIBSCALAPACK = -L${HOME}/ScaLAPACK/lib -lscalapack_BGL.440 -lblacsF77init_MPI-BGL.440-0 -lblacs_MPI-BGL.440-0 -llapack_BGL.440
#--#endif#
NSSLIB = -lc -lnss_files -lnss_dns -lresolv
#LIB = ${LIBSCALAPACK} ${LIBBLAS} ${LIBMPI} ${NSSLIB}
#LIB = ${LIBBLAS} ${LIBSCALAPACK} ${LIBMPI} ${NSSLIB}
BL_LIB = ${LIBBLAS} ${LIBSCALAPACK} ${LIBBLAS} ${LIBMPI} ${NSSLIB}
#
# ===============  Additional Files 
#
EXTRA_BASE= iterate.o cmpi.o
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

