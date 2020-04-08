#doc Include file for the Global Array build on Sgi Origin with MPISpro compiler
#doc This uses the scs numerical library by default.
#doc
#doc Options:
#doc blas - selected by default
#doc vb   - include the Valence Bond module
#doc zora - include Zeroth Order Relativistic Approximation module
#
#dopt ga mpi ci i8 peigs blas vb zora
#opt  vb zora
#
# ================ M4 Processing options
#
MACHINE_KEY=g
MACHOPT=sgi,cio,blas,signals,mips4,doublebackslash,upck-equiv,r10000,64bit,i8

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
FFLAGSTMP = -c

CFLAGSI8= -DLONG_INTEGER -DEXT_INT
FFLAGSI8= -i8

GMAKE  = gmake
ARCHIVE  =  ar rcv
RANLIB   =  ar ts

n32 = -64

KAP   = -WK,-so=1,-o=1
NOKAP = -WK,-so=0,-o=0
ROF	 = -OPT:roundoff=3:IEEE_arithmetic=3 -TENV:X=3
FOLD	 = -OPT:fold_arith_limit=4000
ARCH	 = -r10000


FFLAGSV  = -O3 -G 0 -mips4 $(ARCH) -c $(KAP) $(ROF) -r8 $(n32)  ${FFLAGSI8}
FFLAGSV3 = -O3 -G 0 -mips3 $(ARCH) -c $(KAP) $(ROF) -r8 $(n32)  ${FFLAGSI8}
FFLAGS2  = -O2 -G 0 -mips4 $(ARCH) -c $(KAP) $(ROF) -r8 $(n32) ${FFLAGSI8}
FFLAGSS  = -O1 -G 0 -mips4 $(ARCH) -c $(KAP) $(ROF) -r8 $(n32) ${FFLAGSI8}
FFLAGS0  = -O0 -G 0 -mips4 $(ARCH) -c $(KAP) $(ROF) -r8 $(n32) ${FFLAGSI8}
FFLAGSS0 = -O1 -G 0 -mips4 $(ARCH) -c $(NOKAP) -r8 $(n32) ${FFLAGSI8}
FFLAGS1  = -O3 -G 0 -mips4 $(ARCH) -c $(NOKAP) -r8 $(n32) ${FFLAGSI8}
FFLAGSVN = -O3 -G 0 -mips4 $(ARCH) -c $(NOKAP) $(ROF) $(n32) -r8 ${FFLAGSI8}
CFLAGS 	  = ${CFLAGSI8} -O3 -G 0 $(n32) -mips4 $(ARCH) -c $(KAP) $(ROF) 

CPP	= /usr/lib/cpp
RANLIB	= echo no need for ranlib for library
LDFLAGS	= $(n32) -mips4 $(ARCH) ${FFLAGSI8}


# MPI libraries
MPI_LIBS = -lmpi


# Global Array variables
GA_TARGET=SGITFP
GA_INTEGER_DECL=integer*8
GA_LOGICAL_DECL=logical*8
GA_VERSION_PAR= GA_VERSION=SHMEM USE_MPI=YES 
GA_TARGET_CPU_PAR= GA_TARGET_CPU=R10000

# Bring it all together
OPTIONS=${MACHOPT}
BL_LIB= ${MPI_LIBS} -lscs_i8
#
# ===============  Extra Files
#
EXTRA_BASE=stvint.o hamd1.o cteisg.o n21g.o eigen.o ffun.o \
	sfun2.o sfun4.o adapt.o iterate.o
EXTRA_BENCH=stvint.o hamd1.o n21g.o eigen.o ffun.o \
	sfun2.o sfun4.o adapt.o
EXTRA_MP2=rhsgvb.o drcfun.o fockb.o
EXTRA_CI=tabdav.o gabcd.o genlim.o
#--#if vb#
VBEXTRA = poporb.o
EXTRA = calc1.o
#--#endif#

#
# ===============  Compiler Exceptions
#

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
#
# bench_7 , dircta.f genlim + (r1000) gabcd
genlim.o:	dircta.m
	cat ../utilities/gener.m dircta.m | $(M4) -DGEN_EXTRACTFILE=genlim $(M4OPTS)  > genlim.f
	$(FC) $(FFLAGS1) genlim.f
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
nbo.o:	nbo.m
	cat ../utilities/gener.m nbo.m | $(M4) $(M4OPTS)  > nbo.f
	$(FC) $(FFLAGS2) nbo.f
#
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

#
#  ==========  Exceptions for SGI R10K, R12K, R5k ====
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

calc1.o:        casa.m  
	cat ../utilities/gener.m casa.m | $(M4) -DGEN_EXTRACTFILE=calc1 $(M4OPTS)  > calc1.f 
	$(FC) $(FFLAGS1) calc1.f 

#--#endif#

