#
#doc This is the machine-specific file for the GLOBAL ARRAY based parallel version
#doc This version is tested on Huygens.
#doc it implicitly includes i8 and 64bit options for compatibility with GA
#doc syslog option is for writing a start/stop message in /var/log/messages
#
#dopt ga mpi ci peigs zora vb dl-find vdw i8 openib
#opt vampir newscf essl datain sockets debug masscf scalapack mp2 base syslog
#
# ================ M4 Processing options
#
MACHINE_KEY=r
MACHOPT=rs6000,cio,unix,doublebackslash,upck-equiv,64bitpointers,glibc

## Machine specific options from: rs6000.m ##
IAND64 = iand($$1,$$2)
IOR64 = ior($$1,$$2)
IXOR64 = ieor($$1,$$2)
EXTRA_M4_DEFINITIONS = -Dleadz='leadzz($$1)'
## end of rs600.m  ##

#
# ===============  Compiler Options
#

FC = mpfort -compiler xlf95_r -qfixed 
LD = mpfort -compiler xlf95_r
FC90 = mpfort -compiler xlf95_r 
CC= mpcc -compiler xlc_r 

# Can't use -D preprocessing for DL-FIND so need to use -WF
# See dl-find/Makefile.in.gamess
DLF-PPFLAGS=-WF,-DGAMESS,-DOLDALLOC

RANLIB=ranlib
ARCHIVE  =  ar rcv
ARCH1=-qarch=auto
ARCH2=-qarch=auto -qtune=auto
LARGEFILES = -D_LARGE_FILES
FLAGS64 = -q64
##CPP = /usr/lib/cpp ## on aix
#CPP = /usr/bin/cpp  ## sara 4.1.x version
#CPP = /sara/sw/gcc/4.2.1/bin/cpp ## new sara version
CPP=cpp
#####NOCU=-DNOC_
EXTNAME=-qEXTNAME

#--#if i8#
CFLAGSI8 =  -DLONG_INTEGER -DEXT_INT -DXLCLINUX
FFLAGSI8 = -qintsize=8 
I8_M4_OPTS=,i8
#--#else#
CFLAGSI8 = -DSTD_INT
FFLAGSI8 =
I8_M4_OPTS=
#--#endif i8#

FFLAGSTMP= -c $(FFLAGSI8) $(EXTNAME) $(FLAGS64)
CFLAGSTMP= -c $(CFLAGSI8) $(FLAGS64) $(NOCU) $(LARGEFILES)

#--#if scalapack#
############ USING SCALAPACK (implies external BLAS) ####
#---#if i8#
Cannnot have i8 and SCALAPACK!!!
#---#endif i8#
LIBBLAS= -lessl
BLASOPT=,blas
SCALAPACK_LIB = -L/usr/local/lib -lblacsF77init -lblacs -lscalapack -ltools -llapack
############ END OF SCALAPACK OPTIONS ###################
#--#elseif essl#
############ USING EXTERNAL BLAS ########################
#---#if i8#
LIBBLAS=  -lessl6464
#---#else#
LIBBLAS= -L/usr/local/lib -ltools -llapack -lessl
#---#endif i8#
BLASOPT=,blas
############ END OF EXTERNAL BLAS OPTIONS ###############
#--#else#
LIBBLAS=
BLASOPT=
############ END OF NON-EXTERNAL BLAS ###################
#--#endif#


#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g -qxlf77=leadzero
FFLAGSS = ${FFLAGSTMP} -g -qxlf77=leadzero
FFLAGSN = ${FFLAGSTMP} -g -qxlf77=leadzero
FFLAGS0 = ${FFLAGSTMP} -g -qxlf77=leadzero
CFLAGS  = $(CFLAGSTMP) $(ARCH1) -g
LDOPT   = -g
#
#--#else# 
#
FFLAGSV = $(FFLAGSTMP) -O3 -qmaxmem=-1 -c $(ARCH2) -qxlf77=leadzero
FFLAGSS = $(FFLAGSTMP) -O $(ARCH1) -qxlf77=leadzero
FFLAGSN = $(FFLAGSTMP) -g $(ARCH1) -qxlf77=leadzero
FFLAGS0 = $(FFLAGSTMP) $(ARCH1) -qxlf77=leadzero
CFLAGS  = $(CFLAGSTMP) $(ARCH1) -O0
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSV)
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSS)
LDOPT     = -O3
#--#endif#

#--#if static#
echo static unimplemented
exit
#--#else#
#--#if openib#
LDFLAGS  = ${LDOPT} -q64 -libverbs -lpthread -lrt
#--#else#
LDFLAGS  = ${LDOPT} -q64 
#--#endif#
#--#endif#

# MPI variables Using blanks for mpfort/mpcc ...
MPI_INCLUDE = .
MPI_LIB = ""
LIBMPI = ""

#Global Array stuff
GA_F77_DEFS = -traditional
GA_TARGET=LINUX64
PEIGS_TARGET=LINUX64

#--#if openib#
	IB_INCLUDE=/usr/include/infiniband
	GA_VERSION_PAR= GA_VERSION=SHMEM ARMCI_NETWORK=OPENIB IB_INCLUDE=${IB_INCLUDE} LIBMPI=""
#--#else#
	#--#if lapi64#
		GA_VERSION_PAR= GA_VERSION=SHMEM ARMCI_NETWORK=LAPI64 LIBMPI=""
	#--#else#
		#--#if vampir#
			GA_VERSION_PAR=GA_USE_VAMPIR=yes LIBMPI=""
			LIBVT = -L/hpcx/usrlocal/packages/vampir/lib64 -lVT
			INCVT = /hpcx/usrlocal/packages/vampir/include
			VTOPT=,vampir
		#--#else#
#			Sockets
			GA_VERSION_PAR= GA_VERSION=SHMEM ARMCI_NETWORK=SOCKETS LIBMPI=""
		#--#endif#
	#--#endif#
#--#endif#


OPTIONS=${MACHOPT}$(I8_M4_OPTS)${BLASOPT}${VTOPT}
BL_LIB = ${SCALAPACK_LIB} ${LIBBLAS} ${LIBVT} 
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

