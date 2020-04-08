#doc Machine-dependent settings for serial IBM AIX Power4 using XLF compilers.
#doc This is a 32-bit build with integer*4 and with -qEXTNAME
#doc
#doc Options:
#doc essl  - link against IBM's ESSL blas libarary
#doc 64bit - create a 64-bit binary
#
#dopt mrdci nbo drf zora vb vdw sysmo dl-find essl 64bit
#opt debug static mopac
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
FC = xlf_r
LD = xlf90_r
FC90 = xlf90_r
CC= xlc_r

# Can't use -D preprocessing for DL-FIND so need to use -WF
# See dl-find/Makefile.in.gamess
DLF-PPFLAGS=-WF,-DGAMESS,-DOLDALLOC

#--#if static#
# Static is horrible to get working. I needed to include libcrypt as below
# and then create a file called "crypt.exp" in m4 which contained the text
# between <start> and <end> (without the first hash & space)
# <start>
# #!
# __crypt_r
# __encrypt_r
# __setkey_r
# <end>
# 
# This was culled together from:
# http://publib.boulder.ibm.com/infocenter/comphelp/v9v111/index.jsp?topic=/com.ibm.xlf111.aix.doc/xlfcr/linkingpgms.htm
#
# and
#
# http://unix.derkeiler.com/Newsgroups/comp.unix.aix/2004-11/0553.html
FFLAGS_STATIC = -bstatic
CFLAGS_STATIC = -bstatic
LDFLAGS_STATIC = -bstatic -bnso  -bI:/usr/lib/syscalls.exp /usr/lib/libcrypt.a -bI:crypt.exp
#--#else#
FFLAGS_STATIC=
CFLAGS_STATIC=
LDFLAGS_STATIC = 
#--#endif static#

#--#if 64bit#
ARCHIVE = ar -X64 rcv
FLAGS64 = -q64
#--#else#
ARCHIVE  =  ar rcv
#--#endif#

RANLIB = ranlib
#ARCH1 =-qarch=pwr4
ARCH1 =-qarch=auto -qtune=auto
ARCH1 =-qarch=com
LDFLAGS0 = -b loadmap:load.map -b bigtoc
####NOCU=-DNOC_

FFLAGSTMP = -c $(FFLAGS_STATIC) $(FFLAGSI8) $(FLAGS64) -qEXTNAME
CFLAGSTMP = -c $(CFLAGS_STATIC) $(CFLAGSI8) $(FLAGS64) -D_LARGE_FILES

#--#if debug#
FFLAGSV = $(FFLAGSTMP) -g -qxlf77=leadzero
FFLAGSS = $(FFLAGSTMP) -g -qxlf77=leadzero
FFLAGSN = $(FFLAGSTMP) -g -qxlf77=leadzero
FFLAGS0 = $(FFLAGSTMP) -g -qxlf77=leadzero
FFLAGS1 = $(FFLAGSTMP) -g -qxlf77=leadzero
CFLAGS  = $(CFLAGSTMP) -g $(ARCH1) -O  $(NOCU) 
#--#else# 
FFLAGSV = $(FFLAGSTMP) -O3  -qnosave $(ARCH1) -qxlf77=leadzero
FFLAGSS = $(FFLAGSTMP) -O  $(ARCH1) -qxlf77=leadzero
FFLAGS1 =$(FFLAGSTMP)  -O  $(ARCH1) -qxlf77=leadzero
FFLAGSN = $(FFLAGSTMP) -g  $(ARCH1) -qxlf77=leadzero
FFLAGS0 = $(FFLAGSTMP)  $(ARCH1) -qxlf77=leadzero
FFLAGSV90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSV)
FFLAGSS90 = -qlanglvl=95pure -qsuffix=f=f90 $(FFLAGSS)
CFLAGS =  $(CFLAGSTMP) -O
#--#endif#

#--#if essl#
LIBBLAS= -lessl
BLASOPT=,blas
# The current version of essl exports names with an without an underscore
# so the below is not needed - the below is kept just in case it becomes
# useful again at some point in the future.
RENAMES=\
#-b rename:.ddot_,.ddot \
#-b rename:.daxpy_,.daxpy \
#-b rename:.dswap_,.dswap \
#-b rename:.idamax_,.idamax \
#-b rename:.dnrm2_,.dnrm2 \
#-b rename:.dgemm_,.dgemm \
#-b rename:.dgemv_,.dgemv \
#-b rename:.dcopy_,.dcopy \
#-b rename:.dgthr_,.dgthr \
#-b rename:.dscal_,.dscal \
#-b rename:.icopy_,.icopy \
#-b rename:.dsum_,.dsum \
#-b rename:.idmax_,.idmax \
#-b rename:.igthr_,.igthr \
#-b rename:.idamin_,.idamin \
#-b rename:.dasum_,.dasum \
#-b rename:.drot_,.drot  \
#-b rename:.xerbla_,.xerbla  \
#-b rename:.lsame_,.lsame  \
#-b rename:.dsctr_,.dsctr
#--#endif essl#

LDFLAGS = $(LDFLAGS0) $(FLAGS64) $(RENAMES) $(LDFLAGS_STATIC)

###OPTIONS=${MACHOPT}${BLASOPT},64bitpointers
OPTIONS=${MACHOPT}${BLASOPT}

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
#  ==================  Exceptions for IBM r6000 power2,3,4  ====================
#
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
# If basis is compiled with any optimisation c2024_a & b 
# give incorrect results
basis.o:	basis.f
	$(FC) $(FFLAGSN) $*.f
#
#  ========== Exceptions for IBM r6000 power1,2,3 and 4 (and SP2) ===============
#
integ2e.o:	integ2e.m
	cat ../utilities/gener.m integ2e.m | $(M4) $(M4OPTS) > integ2e.f
	$(FC) $(FFLAGSS) integ2e.f

xc.o:	xc.m
	cat ../utilities/gener.m xc.m | $(M4) $(M4OPTS) > xc.f
	$(FC) $(FFLAGSS) xc.f

