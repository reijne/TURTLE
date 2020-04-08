#doc Makefile settings for serial build on Sgi Origin with MPISpro compiler
#doc By default this build expects to link against the scs numerical library.
#
#doc Options:
#doc No options for this build - blas selected by default
#
#dopt zora vb drf nbo mrdci mopac sysmo vdw dl-find blas
#opt  i8
#
# ================ M4 Processing options
#
MACHINE_KEY=g
MACHOPT=sgi,cio,blas,signals,mips4,doublebackslash,upck-equiv,r10000

## Machine dependant stuff from sgitfp.m ##
REAL=real*8
COMPLEX=complex*16
DFLOAT=dble($$1)
IAND32=iand($$1,$$2)
IOR32=ior($$1,$$2)
IXOR32=ieor($$1,$$2)
IAND64=iand($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)
SHIFT=shiftc($$1,$$2)

#
# ===============  Compiler Options
#

FC = f77
LD = f77
FC90 = f90

GMAKE  = gmake
ARCHIVE  =  ar rcv
RANLIB   =  ar ts

#--#if i8#
n32 = -n32
FFLAGSI8=
CFLAGSI8=
M4_I8_OPT=
#--#else#
n32 = -n64
FFLAGSI8= -i8
CFLAGSI8= -DLONG_INTEGER
M4_I8_OPT=,i8
#--#endif i8#

FFLAGSTMP = -c ${FFLAGSI8}
CFLAGSTMP = -c ${CFLAGSI8}

#--#if debug#
#to do...
#--#else# 
KAP   = -WK,-so=1,-o=1
NOKAP = -WK,-so=0,-o=0
ROF	 = -OPT:roundoff=3:IEEE_arithmetic=3 -TENV:X=3
FOLD	 = -OPT:fold_arith_limit=4000
ARCH	 = -r10000

FFLAGSV  = -O3 -G 0 -mips4 $(ARCH) -c $(KAP) $(ROF) -r8 $(n32)
FFLAGSV3 = -O3 -G 0 -mips3 $(ARCH) -c $(KAP) $(ROF) -r8 $(n32)
FFLAGS2  = -O2 -G 0 -mips4 $(ARCH) -c $(KAP) $(ROF) -r8 $(n32)
FFLAGSS  = -O1 -G 0 -mips4 $(ARCH) -c $(KAP) $(ROF) -r8 $(n32)
FFLAGS0  = -O0 -G 0 -mips4 $(ARCH) -c $(KAP) $(ROF) -r8 $(n32)
FFLAGSS0 = -O1 -G 0 -mips4 $(ARCH) -c $(NOKAP) -r8 $(n32)
FFLAGS1  = -O3 -G 0 -mips4 $(ARCH) -c $(NOKAP) -r8 $(n32)
FFLAGSVN = -O3 -G 0 -mips4 $(ARCH) -c $(NOKAP) $(ROF) $(n32) -r8
CFLAGS 	  = ${CFLAGSTMP} -O3 -G 0 $(n32) -mips4 $(ARCH) $(KAP) $(ROF) 

CPP	= /usr/lib/cpp
RANLIB	= echo no need for ranlib for library
LDFLAGS	= $(n32) -mips4 $(ARCH) ${FFLAGSI8}
#--#endif#

#
# ===============  Additional Files 
#

#--#if blas#
#---#if i8#
LIBBLAS = -lscs_i8
#---#else#
LIBBLAS = -lscs
#---#endif i8#
BLASOPT=,blas
#--#endif blas#

OPTIONS=${MACHOPT}${BLASOPT}${I8_M4_OPT}
BL_LIB=${LIBBLAS}

#--#if vb#
VBEXTRA = poporb.o
#--#endif#

#--#if mrdci#
EXTRA_MRDCI=aftci.o
#--#endif#

EXTRA_BASE=stvint.o hamd1.o cteisg.o n21g.o eigen.o ffun.o \
	sfun2.o sfun4.o adapt.o iterate.o $(EXTRA_MRDCI)
EXTRA_BENCH=stvint.o hamd1.o n21g.o eigen.o ffun.o \
	sfun2.o sfun4.o adapt.o 
EXTRA_MP2=rhsgvb.o drcfun.o fockb.o
EXTRA=outtda.o calc1.o trasor.o smrd0.o \
	dmrd80.o tql2.o 
EXTRA_CI=tabdav.o gabcd.o genlim.o


#
# ===============  Compiler Exceptions
#

smrd0.o:	mrdci3.m
	cat ../utilities/gener.m mrdci3.m | $(M4) -DGEN_EXTRACTFILE=smrd0 $(M4OPTS) > smrd0.f
	$(FC) $(FFLAGSS) smrd0.f

#
# intega. bench_3
#
stvint.o:	intega.m
	cat ../utilities/gener.m intega.m | $(M4) -DGEN_EXTRACTFILE=stvint $(M4OPTS)  > stvint.f
	$(FC) $(FFLAGS0) stvint.f
#
# drv1e. bench_10
#
hamd1.o:	drv1e.m
	cat ../utilities/gener.m drv1e.m | $(M4) -DGEN_EXTRACTFILE=hamd1 $(M4OPTS)  > hamd1.f
	$(FC) $(FFLAGS2) hamd1.f
#
# analb. c2010_a
#
cteisg.o:	analb.m
	cat ../utilities/gener.m analb.m | $(M4) -DGEN_EXTRACTFILE=cteisg $(M4OPTS)   > cteisg.f
	$(FC) $(FFLAGS2) cteisg.f
#
# basis. c2024_a
#
n21g.o:	basis.m
	cat ../utilities/gener.m basis.m | $(M4) -DGEN_EXTRACTFILE=n21g $(M4OPTS)  > n21g.f
	$(FC) $(FFLAGS1) n21g.f
#
# optim. c3009
#
eigen.o:	optim.m
	cat ../utilities/gener.m optim.m | $(M4) -DGEN_EXTRACTFILE=eigen $(M4OPTS)  > eigen.f
	$(FC) $(FFLAGSS) eigen.f
 
#   mips4 fouls up (see ex c2ttt)
integs.o:	integs.m
	cat ../utilities/gener.m integs.m | $(M4) $(M4OPTS)  > integs.f
	$(FC) $(FFLAGSV3) integs.f
#
# Compile time
outtda.o : tdaf.m
	cat ../utilities/gener.m tdaf.m | $(M4) -DGEN_EXTRACTFILE=outtda $(M4OPTS)  > outtda.f
	$(FC) $(FFLAGS1) outtda.f
#
# bench_5 , casa.f
calc1.o:	casa.m
	cat ../utilities/gener.m casa.m | $(M4) -DGEN_EXTRACTFILE=calc1 $(M4OPTS)  > calc1.f
	$(FC) $(FFLAGS1) calc1.f
#
# bench_6 , mcscfb.f
trasor.o:	mcscfb.m
	cat ../utilities/gener.m mcscfb.m | $(M4) -DGEN_EXTRACTFILE=trasor $(M4OPTS)  > trasor.f
	$(FC) $(FFLAGSS) trasor.f
#
# bench_7 , dircta.f genlim + (r1000) gabcd
genlim.o:	dircta.m
	cat ../utilities/gener.m dircta.m | $(M4) -DGEN_EXTRACTFILE=genlim $(M4OPTS)  > genlim.f
	$(FC) $(FFLAGS1) genlim.f
#
mrdci4.o:	mrdci4.m
	cat ../utilities/gener.m mrdci4.m | $(M4) $(M4OPTS)  > mrdci4.f
	$(FC) $(FFLAGS2) mrdci4.f
#
mrdci5.o:	mrdci5.m
	cat ../utilities/gener.m mrdci5.m | $(M4) $(M4OPTS)  > mrdci5.f
	$(FC) $(FFLAGSS0) mrdci5.f
#
# dmrd80 (c11019)
#
dmrd80.o:	mrdci7.m
	cat ../utilities/gener.m  mrdci7.m | $(M4) -DGEN_EXTRACTFILE=dmrd80 $(M4OPTS)  > dmrd80.f
	$(FC) $(FFLAGSS) dmrd80.f
#
mrdci6.o:	mrdci6.m
	cat ../utilities/gener.m mrdci6.m | $(M4) $(M4OPTS)  > mrdci6.f
	$(FC) $(FFLAGSS) mrdci6.f
#
dirrpa.o:	dirrpa.m
	cat ../utilities/gener.m dirrpa.m | $(M4) $(M4OPTS)  > dirrpa.f
	$(FC) $(FFLAGS2) dirrpa.f
#
tabdav.o:	dirctb.m
	cat ../utilities/gener.m dirctb.m | $(M4) -DGEN_EXTRACTFILE=tabdav $(M4OPTS)  > tabdav.f
	$(FC) $(FFLAGS0) tabdav.f
#
# direct, bench_12 
fockb.o:	direct.m
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=fockb $(M4OPTS)   > fockb.f
	$(FC) $(FFLAGS2) fockb.f
#
rhsgvb.o:	cphf.m
	cat ../utilities/gener.m  cphf.m | $(M4) -DGEN_EXTRACTFILE=rhsgvb $(M4OPTS)  > rhsgvb.f
	$(FC) $(FFLAGS2) rhsgvb.f
#
mp2_parallel.o:	mp2_parallel.m
	cat ../utilities/gener.m  mp2_parallel.m| $(M4) $(M4OPTS)  > mp2_parallel.f
	$(FC) $(FFLAGS2) mp2_parallel.f
#
iterate.o: morokuma.m
	cat ../utilities/gener.m  morokuma.m | $(M4) -DGEN_EXTRACTFILE=iterate $(M4OPTS)   > iterate.f
	$(FC) $(FFLAGSS) iterate.f
#
gabcd.o:	dircta.m
	cat ../utilities/gener.m  dircta.m | $(M4) -DGEN_EXTRACTFILE=gabcd $(M4OPTS)  > gabcd.f
	$(FC) $(FFLAGS2) gabcd.f
#
drcfun.o: direct.m
	cat ../utilities/gener.m  direct.m | $(M4) -DGEN_EXTRACTFILE=drcfun $(M4OPTS)  > drcfun.f
	$(FC) $(FFLAGSS) drcfun.f

#
analf.o:	analf.m
	cat ../utilities/gener.m analf.m | $(M4) $(M4OPTS)  > analf.f
	$(FC) $(FFLAGS2) analf.f
#
analg.o:	analg.m
	cat ../utilities/gener.m analg.m | $(M4) $(M4OPTS)  > analg.f
	$(FC) $(FFLAGS2) analg.f
#
mcscfa.o:	mcscfa.m
	cat ../utilities/gener.m mcscfa.m | $(M4) $(M4OPTS)  > mcscfa.f
	$(FC) $(FFLAGS2) mcscfa.f
#
mcscfb.o:	mcscfb.m
	cat ../utilities/gener.m mcscfb.m | $(M4) $(M4OPTS)  > mcscfb.f
	$(FC) $(FFLAGS2) mcscfb.f
#
mcscfc.o:	mcscfc.m
	cat ../utilities/gener.m mcscfc.m | $(M4) $(M4OPTS)  > mcscfc.f
	$(FC) $(FFLAGS2) mcscfc.f
#
nbo.o:	nbo.m
	cat ../utilities/gener.m nbo.m | $(M4) $(M4OPTS)  > nbo.f
	$(FC) $(FFLAGS2) nbo.f
#
# chap2
mclr.o:	mclr.m
	cat ../utilities/gener.m mclr.m | $(M4) $(M4OPTS)  > mclr.f
	$(FC) $(FFLAGS2) mclr.f
#
adapt.o:	guess.m
	cat ../utilities/gener.m guess.m | $(M4) -DGEN_EXTRACTFILE=adapt $(M4OPTS)  > adapt.f
	$(FC) $(FFLAGS2) adapt.f
#
ffun.o: intega.m
	cat ../utilities/gener.m  intega.m | $(M4) -DGEN_EXTRACTFILE=ffun $(M4OPTS)  > ffun.f
	$(FC) $(FFLAGS2) ffun.f
#
sfun2.o: optim.m
	cat ../utilities/gener.m  optim.m | $(M4) -DGEN_EXTRACTFILE=sfun2 $(M4OPTS)  > sfun2.f
	$(FC) $(FFLAGS2) sfun2.f
#
sfun4.o: optim.m
	cat ../utilities/gener.m  optim.m | $(M4) -DGEN_EXTRACTFILE=sfun4 $(M4OPTS)  > sfun4.f
	$(FC) $(FFLAGS2) sfun4.f
# ...     iky gets overwritten (??? jvl)
zora.o:	zora.m
	cat ../utilities/gener.m zora.m | $(M4) $(M4OPTS)  > zora.f
	$(FC) $(FFLAGS2) zora.f
#
tql2.o: eispack.m
	cat ../utilities/gener.m  eispack.m | $(M4) -DGEN_EXTRACTFILE=tql2 $(M4OPTS)  > tql2.f
	$(FC) $(FFLAGSS) tql2.f

newmrd2.o:	newmrd2.f
	$(FC) $(FFLAGSVN) $*.f

newmrd3.o:	newmrd3.f
	$(FC) $(FFLAGSVN) $*.f
#

# MRDCI exceptions
#--#if mrdci#
aftci.o: newmrd5.m
	cat ../utilities/gener.m  newmrd5.m | $(M4) -DGEN_EXTRACTFILE=aftci $(M4OPTS)   > aftci.f
	$(FC) $(FFLAGS2) aftci.f
#--#endif#
#
#  DFT Exceptions for SGI R10K, R12K, R5k ====
#
xc.o:	xc.m
	cat  ../utilities/gener.m xc.m | $(M4) $(M4OPTS) > xc.f
	$(FC) $(FFLAGS2) xc.f

exp_dksm_hess.o:	exp_dksm_hess.m
	cat  ../utilities/gener.m exp_dksm_hess.m | $(M4) $(M4OPTS) > exp_dksm_hess.f
	$(FC) $(FFLAGS2) exp_dksm_hess.f
#
# VB Exceptions
#--#if vb#
vbdens.o:       vbdens.m
	cat ../machines/$(MACH) ../utilities/gener.m vbdens.m | $(M4) $(M4OPTS) > vbdens.f
	$(FC) $(FFLAGS2) vbdens.f

poporb.o:       vbscf.m
	cat ../machines/$(MACH) ../utilities/gener.m  vbscf.m | $(M4) -DGEN_EXTRACTFILE=poporb $(M4OPTS) $(SNGL) > poporb.f
	$(FC) $(FFLAGS2) poporb.f
#--#endif#



