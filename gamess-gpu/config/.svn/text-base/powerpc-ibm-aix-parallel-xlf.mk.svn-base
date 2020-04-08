#
#doc This is the machine-specific file for the GLOBAL ARRAY based parallel version
#doc on HPCx
#doc it implicitly includes i8 and 64bit options for compatibility with GA
#
#dopt ga mpi mp2 peigs newscf dl-find vdw masscf essl scalapack datain 
#opt i8 vampir zora vb
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

FC = mpxlf_r
LD = mpxlf_r
FC90 = mpxlf90_r
CC= mpcc_r

# Can't use -D preprocessing for DL-FIND so need to use -WF
# See dl-find/Makefile.in.gamess
DLF-PPFLAGS=-WF,-DGAMESS,-DOLDALLOC

RANLIB=ranlib
ARCHIVE  =  ar -X 64 rcv
ARCH1=-qarch=auto
ARCH2=-qarch=auto -qtune=auto
LARGEFILES = -D_LARGE_FILES
FLAGS64 = -q64
CPP = /usr/lib/cpp
#####NOCU=-DNOC_
EXTNAME=-qEXTNAME

#--#if i8#
CFLAGSI8 =  -DLONG_INTEGER -DEXT_INT
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
Cannot currently use i8 with essl
#---#else#
LIBBLAS= -lessl
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
CFLAGS  = $(CFLAGSTMP) $(ARCH1) -O
#
#--#else# 
#
FFLAGSV = $(FFLAGSTMP) -O3 -qmaxmem=-1 -c $(ARCH2) -qxlf77=leadzero
FFLAGSS = $(FFLAGSTMP) -O $(ARCH1) -qxlf77=leadzero
FFLAGSN = $(FFLAGSTMP) -g $(ARCH1) -qxlf77=leadzero
FFLAGS0 = $(FFLAGSTMP) $(ARCH1) -qxlf77=leadzero
CFLAGS  = $(CFLAGSTMP) $(ARCH1) -O
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSV)
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSS)
#--#endif#

#--#if static#
echo static unimplemented
exit
#--#else#
LDFLAGS  = -q64 \
#--#if vampir#
	   -b rename:.vtclassdef_,.vtclassdef \
	   -b rename:.vtfuncdef_,.vtfuncdef \
	   -b rename:.vttraceon_,.vttraceon \
	   -b rename:.vttraceoff_,.vttraceoff \
	   -b rename:.vtbegin_,.vtbegin \
	   -b rename:.vtend_,.vtend \
#--#endif#
#--#if scalapack#
#----#if newscf#
	   -b rename:.indxg2l_,.indxg2l \
	   -b rename:.blacs_gridexit_,.blacs_gridexit \
	   -b rename:.blacs_gridmap_,.blacs_gridmap \
	   -b rename:.pdelget_,.pdelget \
	   -b rename:.pdgemv_,.pdgemv \
	   -b rename:.pdgetri_,.pdgetri \
	   -b rename:.pzgemm_,.pzgemm \
	   -b rename:.pdgemm_,.pdgemm \
	   -b rename:.pdgemr2d_,.pdgemr2d \
#----#endif newscf#
	   -b rename:.blacs_get_,.blacs_get \
	   -b rename:.blacs_pinfo_,.blacs_pinfo \
	   -b rename:.blacs_gridinit_,.blacs_gridinit \
	   -b rename:.blacs_gridinfo_,.blacs_gridinfo \
	   -b rename:.numroc_,.numroc \
	   -b rename:.descinit_,.descinit \
	   -b rename:.pdgetrf_,.pdgetrf \
	   -b rename:.pdpotrf_,.pdpotrf \
	   -b rename:.pdpotri_,.pdpotri \
	   -b rename:.pdpotrs_,.pdpotrs \
	   -b rename:.pdgetrs_,.pdgetrs \
	   -b rename:.iceil_,.iceil \
	   -b rename:.pdlamch_,.pdlamch \
	   -b rename:.pdsyev_,.pdsyev \
	   -b rename:.pdsyevx_,.pdsyevx \
	   -b rename:.pdsyevd_,.pdsyevd \
	   -b rename:.indxg2p_,.indxg2p \
#--#elseif newscf#
	   -b rename:.blacs_get_,.blacs_get \
	   -b rename:.indxg2p_,.indxg2p \
	   -b rename:.indxg2l_,.indxg2l \
	   -b rename:.blacs_gridexit_,.blacs_gridexit \
	   -b rename:.blacs_gridmap_,.blacs_gridmap \
	   -b rename:.blacs_gridinfo_,.blacs_gridinfo \
	   -b rename:.numroc_,.numroc \
	   -b rename:.descinit_,.descinit \
	   -b rename:.pdelget_,.pdelget \
	   -b rename:.pdgemv_,.pdgemv \
	   -b rename:.pdgetrf_,.pdgetrf \
	   -b rename:.pdgetri_,.pdgetri \
	   -b rename:.pdsyevd_,.pdsyevd \
	   -b rename:.pdsyev_,.pdsyev \
	   -b rename:.pdgemm_,.pdgemm \
	   -b rename:.pzgemm_,.pzgemm \
	   -b rename:.pdgemr2d_,.pdgemr2d \
#--#endif newscf#
#--#if vb#
	   -b rename:.mp_stdout_mode_,.mp_stdout_mode \
#--#endif vb#
	   -b loadmap:load.map -b bigtoc
#--#endif#

# MPI variables
MPI_INCLUDE = .
MPI_LIB =
LIBMPI = 

#Global Array stuff
GA_F77_DEFS = -traditional
GA_TARGET=LAPI64
PEIGS_TARGET=LAPI64

# Mess around with the GA and PeIGS version based on Vampir and MPI
#---#if vampir#
GA_VERSION_PAR=GA_USE_VAMPIR=yes
LIBVT = -L/hpcx/usrlocal/packages/vampir/lib64 -lVT
INCVT = /hpcx/usrlocal/packages/vampir/include
VTOPT=,vampir
#---#else#
GA_VERSION_PAR=
#---#endif vampir#

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

