#doc Makefile settings for Global Array build on MPI Ultra Sparc Sun-Fire-15000 
#doc running Solaris 2.9 with  Forte Developer 7 Fortran 95 7.0
#doc Default location for MPI libraries is: /opt/SUNWhpc/
#doc Please set MPI_INC & MPI_LIBS if your libraries are elsewhere
#doc
#doc Options
#doc datain - force GAMESS-UK to read input from a file called datain instead of from standard input
#
#dopt ga mpi i8 peigs ci datain vb zora
#opt 
#
# ================ M4 Processing options
#
MACHINE_KEY=s
MACHOPT=sun,ultra,cio,unix,signals,doublebackslash,upck-equiv,64bit

FC = mpf77
LD = mpf77
FC90 = mpf77
CC= cc

# See Makefile.vb.in to see where this needed.
WHOAMI=/usr/ucb/whoami

#MPI flags
MPI_INC  = /opt/SUNWhpc/include
MPI_L =/opt/SUNWhpc/lib
LIBMPI= -lmpi
MPI_LIBS= -L${MPI_L} ${LIBMPI}

LIBBLAS = -lmvec

#GA Stuff
GA_F77_DEFS = -traditional
GA_TARGET=SOLARIS
GA_VERSION_PAR= GA_VERSION=SHMEM USE_MPI=YES MPI_INCLUDE=$(MPI_INC)  MPI_LIB=$(MPI_L)
GA_TARGET_CPU_PAR= 

# Put the link line together
BL_LIB= ${LLIB} ${lLIB} ${MPI_LIBS} ${LIBBLAS} ${lBLAS}

## M4 Options ##
OPTIONS=${MACHOPT}


# ===============  Compiler Options
## real problems here between serial and parallel ARCH parameters
## and use of SOLARIS64 .. fall back to SOLARIS
## Using SOLARIS64 gives compiler error in GA 3.3 unless
## USE_INTEGER4=y
# ======= check for xarch specification and XO5 vs O5  ========================
TARGET  = ultra3cu
ARCH    = v8plusb
FFLAGSV  = -c -fast -xtarget=$(TARGET) -xarch=$(ARCH) -xO5 -xsafe=mem -fsimple=2
FFLAGSV4 = -c -fnonstd -xtarget=$(TARGET) -xarch=$(ARCH) -O4 -libmil -dalign
FFLAGSS  = -c -fnonstd -xtarget=$(TARGET) -xarch=$(ARCH) -O2 -libmil -dalign
CFLAGS   = ${CFLAGSI8} -xtarget=$(TARGET) -xarch=$(ARCH) -O -c
LDFLAGS= -xlic_lib=sunperf -lsocket -lnsl -xtarget=$(TARGET) -xarch=$(ARCH)

RANLIB=ranlib


EXTRA_BASE= sp0111.o sp1111.o tvder.o lagrng.o graph.o
EXTRA=wtedz0.o 
EXTRA_CI=wtedz0.o

#
#  ======================  Exceptions for UltraSPARC/Solaris ==============
#
sp0111.o:	integs.m
	cat ../utilities/gener.m integs.m | $(M4) -DGEN_EXTRACTFILE=sp0111 $(M4OPTS) > sp0111.f
	$(FC) $(FFLAGSS) sp0111.f
sp1111.o:	integs.m
	cat ../utilities/gener.m integs.m | $(M4) -DGEN_EXTRACTFILE=sp1111 $(M4OPTS) > sp1111.f
	$(FC) $(FFLAGSS) sp1111.f
#
tvder.o:	drv1e.m
	cat ../utilities/gener.m drv1e.m | $(M4) -DGEN_EXTRACTFILE=tvder $(M4OPTS) > tvder.f
	$(FC) $(FFLAGSV4) tvder.f
#
lagrng.o:	scf.m
	cat ../utilities/gener.m scf.m | $(M4) -DGEN_EXTRACTFILE=lagrng  $(M4OPTS) > lagrng.f
	$(FC) $(FFLAGSS) lagrng.f
#
graph.o:	analb.m
	cat ../utilities/gener.m analb.m | $(M4) -DGEN_EXTRACTFILE=graph  $(M4OPTS) > graph.f
	$(FC) $(FFLAGSS) graph.f
#
wtedz0.o:	util7.m
	cat ../utilities/gener.m util7.m | $(M4) -DGEN_EXTRACTFILE=wtedz0  $(M4OPTS) > wtedz0.f
	$(FC) $(FFLAGSS) wtedz0.f
dircta.o:	dircta.m
	cat ../utilities/gener.m dircta.m | $(M4) $(M4OPTS) > dircta.f
	$(FC) $(FFLAGSV4) dircta.f
dirctb.o:	dirctb.m
	cat ../utilities/gener.m dirctb.m | $(M4) $(M4OPTS) > dirctb.f
	$(FC) $(FFLAGSV4) dirctb.f
#
rpa.o:	rpa.m
	cat ../utilities/gener.m rpa.m | $(M4) $(M4OPTS) > rpa.f
	$(FC) $(FFLAGSS) rpa.f
#
optim.o:	optim.m
	cat ../utilities/gener.m optim.m | $(M4) $(M4OPTS) > optim.f
	$(FC) $(FFLAGSS) optim.f
#
guess.o:	guess.f
	$(FC) $(FFLAGSS) $*.f
#
intege.o:	intege.f
	$(FC) $(FFLAGSS) $*.f

#
# Build Information
#
#  uname -a = SunOS frontend 5.9 Generic_112233-12 sun4u sparc SUNW,Sun-Fire-15000
#  f90 -V = f90: Forte Developer 7 Fortran 95 7.0 2002/03/09
