#doc Makefile settings serial IBM AIX Power440 using XLF compilers
#doc Default build is 64-bit with integer*4 and with -qEXTNAME
#doc 
#doc If using 32-bit build on huygens set environment : *** OBJECT_MODE = 32 ***
#doc
#doc If using i8 ensure you remove the essl option unless you have an i8 essl library! - Huygens does
#doc
#doc Options:
#doc essl - use IBM ESSL (currently only availble with i4)
#doc datain - force GAMESS-UK to read input from a file called datain instead of from standard input
#doc i8 - build with 64-bit integers (requries 64bit option to be set too)

#dopt mrdci nbo drf zora vb vdw sysmo dl-find 64bit essl  i8 mopac
#opt debug datain mp2 scf ci
#
#
# ================ M4 Processing options
#
MACHINE_KEY=r
MACHOPT=rs6000,cio,unix,doublebackslash,upck-equiv,rs6000_extname,64bitpointers

IAND64 = iand($$1,$$2)
IOR64 = ior($$1,$$2)
IXOR64 = ieor($$1,$$2)
EXTRA_M4_DEFINITIONS = -Dleadz='leadzz($$1)'

#
# ===============  Compiler Options
#
#--#if 64bit i8#
FC = xlf_r
#--#else#
FC = xlf
## xlf pthread kludge...
#--#endif#
FC90 = xlf90
CC= xlc
CPP=cpp
LD = xlf_r
RANLIB=ranlib
ARCHIVE  =  ar rcv

#--#if dl-find#
# Can't use -D preprocessing for DL-FIND so need to use -WF
# See dl-find/Makefile.in.gamess
DLF-PPFLAGS=-WF,-DGAMESS,-DOLDALLOC
#--#endif dl-find#

# Sort out flags for 64bit and i8 builds
#--#if 64bit i8#
FLAGS64= -q64
#--#else#
FLAGS64=
#--#endif#

#--#if i8#
CFLAGSI8   =  -DLONG_INTEGER
FFLAGSI8   = -qintsize=8
I8_M4_OPT =,i8
#--#else#
CFLAGSI8   =
FFLAGSI8   =
I8_M4_OPT =,64bitpointers
#--#endif#


# Default compiler/linker flags for all compililations
FFLAGSTMP  = -c -qEXTNAME ${FLAGS64} ${FFLAGSI8}
CFLAGSTMP  = -c -D_LARGE_FILES ${FLAGS64} ${CFLAGSI8}
LDFLAGS    = ${FLAGS64} 


#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g -fullpath
FFLAGSS = ${FFLAGSTMP} -g -fullpath
FFLAGS1 = ${FFLAGSTMP} -g -fullpath
FFLAGSN = ${FFLAGSTMP} -g -fullpath
FFLAGS0 = ${FFLAGSTMP} -g -fullpath
CFLAGS  = ${CFLAGSTMP} -g
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSV)
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSS)
#--#elseif 64bit#
FFLAGSV = ${FFLAGSTMP} -O3 -qnosave 
FFLAGSS = ${FFLAGSTMP} -O
FFLAGS1 = ${FFLAGSTMP} -O
FFLAGSN = ${FFLAGSTMP}
FFLAGS0 = ${FFLAGSTMP} -O0 
CFLAGS  = ${CFLAGSTMP}  -O 
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSV)
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSS)
#--#else# 
# Same optimisation flags for 32 bit and i8 build
FFLAGSV = ${FFLAGSTMP} -O2 -qnosave -qarch=auto -qtune=auto
FFLAGSS = ${FFLAGSTMP} -O2 -qnosave -qarch=auto -qtune=auto
FFLAGS1 = ${FFLAGSTMP} -O1 -qarch=auto -qtune=auto
FFLAGSN = ${FFLAGSTMP} -O0 -qarch=auto -qtune=auto
FFLAGS0 = ${FFLAGSTMP} -O0 -qarch=auto -qtune=auto
CFLAGS  = ${CFLAGSTMP} -O
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSV)
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSS)
#--#endif#

#--#if essl#
# For some unfathomable reas n, the mclr case c2028_e fails when we use
# the essl version of dcopy. Someone who understands the code should really
# try and work out what is going on there...
#---#if i8#
LIBBLAS= ../linalg/dcopy.o -lessl6464
#---#else#
LIBBLAS= ../linalg/dcopy.o  -lessl
#---#endif i8#
BLASOPT=,blas
#--#endif blas#

OPTIONS=${MACHOPT}${BLASOPT}${I8_M4_OPT}
BL_LIB = ${LIBBLAS}
#
# ===============  Additional Files 
#
EXTRA_BASE= iterate.o
EXTRA=
#----#if mrdci#
EXTRA=skiny.o foxy.o outp.o
#----#endif mrdci#
EXTRA_MP2= 
#
#  ========== Exceptions for IBM r6000 power series ===============
#
iterate.o: morokuma.m
	cat ../utilities/gener.m  morokuma.m | $(M4) -DGEN_EXTRACTFILE=iterate $(M4OPTS)  > iterate.f
	$(FC) $(FFLAGS0) iterate.f
#
#--#if mrdci#
skiny.o:	newmrd5.m
	cat ../utilities/gener.m  newmrd5.m | $(M4) -DGEN_EXTRACTFILE=skiny $(M4OPTS) $(SNGL) > skiny.f
	$(FC) $(FFLAGSS) skiny.f
foxy.o:	newmrd5.m
	cat ../utilities/gener.m  newmrd5.m | $(M4) -DGEN_EXTRACTFILE=foxy $(M4OPTS) $(SNGL) > foxy.f
	$(FC) $(FFLAGSS) foxy.f
outp.o:	newmrd5.m
	cat ../utilities/gener.m  newmrd5.m | $(M4) -DGEN_EXTRACTFILE=outp $(M4OPTS) $(SNGL) > outp.f
	$(FC) $(FFLAGSS) outp.f             
#--#endif mrdci#
#
#  ==================  Exceptions for IBM r6000 440d  ====================
#
#--#if 64bit i8#
intega.o:	intega.f
	$(FC) $(FFLAGS0) $*.f
#--#else#
# 32 bit exceptions:
#
# -O4 xlf-fails
dircta.o:	dircta.f
	$(FC) $(FFLAGSS) $*.f
fullci.o:	fullci.f
	$(FC) $(FFLAGS1) $*.f
## -O4 incorrect Chap2
mcscfa.o:	mcscfa.f
	$(FC) $(FFLAGS1) $*.f
xc_lib.o: xc_lib.m
	cat ../utilities/gener.m xc_lib.m | $(M4) $(M4OPTS) > xc_lib.f
	$(FC) $(FFLAGSS) xc_lib.f
#--#endif#
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
casb.o:		casb.f
	$(FC) $(FFLAGSS) $*.f
newmrd2.o:	newmrd2.f
	$(FC) $(FFLAGSS) $*.f
newmrd3.o:	newmrd3.f
	$(FC) $(FFLAGSS) $*.f
newmrd4.o:	newmrd4.f
	$(FC) $(FFLAGSS) $*.f
newmrd5.o:	newmrd5.f
	$(FC) $(FFLAGSS) $*.f
newmrd6.o:	newmrd6.f
	$(FC) $(FFLAGSS) $*.f
mrdci2.o:	mrdci2.f
	$(FC) $(FFLAGS1) $*.f
mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGS1) $*.f
mrdci7.o:	mrdci7.f
	$(FC) $(FFLAGS1) $*.f
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

