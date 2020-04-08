#  Machine-dependent makefile settings on IBM Machines. 
#  This is a 64-bit build 
#  This file was created on the Huygens machine at SARA.
#
#doc For IBM pseries64 + Suse using XLF compilers
#doc blas option uses IBM ESSL

#dopt base mpi blas 64bit
#opt nbo drf zora debug dynamic_lb newscf datain qmmm
#
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
## end of rs6000.m  ##

#
# ===============  Compiler Options
#
FC = mpfort
LD = mpfort
FC90 = xlf90
CC= mpcc

FFLAGSTMP = -c
RANLIB=ranlib
ARCHIVE  =  ar rcv
ARCH1=-qarch=auto
ARCH1=-qarch=auto -qflttrap=overflow:zerodivide:invalid
ARCH2=-qarch=auto -qtune=auto
LARGEFILES = -D_LARGE_FILES
EXTNAME=-qEXTNAME
#EXTNAME=
LDFLAGS0  = 
####NOCU=-DNOC_

#--#if 64bit#
ARCHIVE = ar rcv
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
CFLAGS  = $(ARCH1) ${CFLAGSI8} -g -O0 -c ${NOCU} ${FLAGS64} ${LARGEFILES}

#--#endif#

#--#if blas#
#---#if 64bit#
LIBBLAS= -L/usr/lib64 -lessl
#---#else#
LIBBLAS= -L/usr/lib -lessl
#---#endif#
BLASOPT=,blas
#--#endif#


LDFLAGS = $(LDFLAGS0) $(FLAGS64) $(RENAMES) 
#--#if 64bits#
MPI_LIB      = /opt/ibmhpc/ppe.poe/lib/libmpi64
MPI_INCLUDE  = /opt/ibmhpc/ppe.poe/include/thread64
LIBMPI =  -L$(MPI_LIB) -lmpi_ibm -lpoe
#--#else#
MPI_LIB      = /opt/ibmhpc/ppe.poe/lib/libmpi
MPI_INCLUDE  = /opt/ibmhpc/ppe.poe/include/thread
LIBMPI =  -L$(MPI_LIB) -lmpi_ibm -lpoe
#--#endif#
FFLAGSV90 = $(FFLAGSV) -I$(MPI_INCLUDE) $(LIBMPI)
FFLAGSS90 = $(FFLAGSS) -I$(MPI_INCLUDE) $(LIBMPI)

###OPTIONS=${MACHOPT}${BLASOPT},64bitpointers
OPTIONS=${MACHOPT}${BLASOPT}

#--#if newscf#
LIBSCALAPACK = -L/bgl/local/lib -lscalapack -lblacsF77init_MPI-BGENE-1 -lblacs_MPI-BGENE-1 -llapack
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

