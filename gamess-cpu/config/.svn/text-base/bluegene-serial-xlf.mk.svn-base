#doc Makefile settings serial IBM AIX Power440 using XLF compilers
#doc This is a 64-bit build with integer*4 and with -qEXTNAME
#doc
#doc Options:
#doc essl - use IBM ESSL
#doc datain - force GAMESS-UK to read input from a file called datain instead of from standard input

#dopt mrdci nbo drf zora vb vdw sysmo mopac dl-find
#opt debug 64bit datain bluegene blas
#
#
# ================ M4 Processing options
#
MACHINE_KEY=r
MACHOPT=rs6000,cio,unix,doublebackslash,upck-equiv,rs6000_extname

## Machine specific options from: rs6000.m ##
IAND64 = iand($$1,$$2)
IOR64 = ior($$1,$$2)
IXOR64 = ieor($$1,$$2)
EXTRA_M4_DEFINITIONS = -Dleadz='leadzz($$1)'
## end of rs600.m  ##

#
# ===============  Compiler Options
#
#--#if bluegene#
FC = blrts_xlf
LD = blrts_xlf90
FC90 = blrts_xlf90
CC= blrts_xlc
#--#else#
FC = xlf
LD = xlf90
FC90 = xlf90
CC= xlc
#--#endif#

FFLAGSTMP = -c
RANLIB=ranlib
ARCHIVE  =  ar rcv
#--#if bluegene#
#ARCH1=-qarch=440d
#ARCH1=-qarch=440d -qflttrap=overflow:zerodivide:invalid
#ARCH2=-qarch=440d -qtune=440
#--#else#
ARCH1=
#--#endif#
LARGEFILES = -D_LARGE_FILES
EXTNAME=-qEXTNAME
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

#--#if bluegene#
#
# First go back to basics making sure the code is "right"
#
#FFLAGSV = -g -O0 -c -qnosave $(ARCH1) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
#FFLAGSS = -g -O0 -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
#FFLAGS1 = -g -O0 -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
#FFLAGSN = -g -O0 -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
#FFLAGS0 = -g -O0 -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
#CFLAGS  = $(ARCH1) ${CFLAGSI8} -O0 -c ${NOCU} ${FLAGS64} ${LARGEFILES}
#--#else#
FFLAGSV = -g -O3 -c -qnosave $(ARCH1) ${FFLAGSI8} ${EXTNAME}  ${FLAGS64}
FFLAGSS = -g -O  -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGS1 = -g -O  -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGSN = -g  -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
FFLAGS0 = -g -O0 -c $(ARCH1) ${FFLAGSI8} ${EXTNAME} ${FLAGS64}
CFLAGS  = $(ARCH1) ${CFLAGSI8}  -O -c ${NOCU} ${FLAGS64} ${LARGEFILES}
#--#endif#
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSV)
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSS)

#--#endif#

#--#if blas#
LIBBLAS= -L/home/y01/fiona/lib -lessln
BLASOPT=,blas
#--#endif blas#

LDFLAGS = $(LDFLAGS0) $(FLAGS64) $(RENAMES) 
###OPTIONS=${OPTIONS0}${BLASOPT},64bitpointers
OPTIONS=${MACHOPT}${BLASOPT}

#--#if bluegene#
NSSLIB = -lc -lnss_files -lnss_dns -lresolv
#--#endif bluegene#

BL_LIB = ${LIBBLAS} ${NSSLIB}
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

