#doc Makefile settings for MPI Ultra Sparc Sun-Fire-15000 running Solaris 2.9
#doc with  Forte Developer 7 Fortran 95 7.0
#doc
#doc Options:
#doc mopac - include MOPAC in the GAMESS-UK build
#
#dopt drf zora mrdci nbo vb vdw mopac sysmo
#opt 
#
#
# ================ M4 Processing options
#
MACHINE_KEY=s
MACHOPT=sun,ultra,cio,unix,signals,doublebackslash,upck-equiv

# ======= check for xarch specification and XO5 vs O5  ========================
TARGET  = ultra3cu
ARCH    = v9b
FFLAGSV  = -c -fast -xtarget=$(TARGET) -xarch=$(ARCH) -xO5 -xsafe=mem -fsimple=2
FFLAGSV4 = -c -fnonstd -xtarget=$(TARGET) -xarch=$(ARCH) -O4 -libmil -dalign
FFLAGSS  = -c -fnonstd -xtarget=$(TARGET) -xarch=$(ARCH) -O2 -libmil -dalign
CFLAGS   = ${CFLAGSI8} -xtarget=$(TARGET) -xarch=$(ARCH) -O -c
LDFLAGS= -xlic_lib=sunperf -xtarget=$(TARGET) -xarch=$(ARCH)

RANLIB=ranlib
BL_LIB = -lmvec

# See Makefile.vb.in to see where this needed.
WHOAMI=/usr/ucb/whoami


EXTRA_BASE= sp0111.o sp1111.o tvder.o lagrng.o graph.o
EXTRA=wtedz0.o smrd0.o

#
#  ==========  Exceptions for SuperSPARC/Solaris and SGI R5, R8 and R10K ====
#
smrd0.o:	mrdci3.m
	cat ../utilities/gener.m mrdci3.m | $(M4) -DGEN_EXTRACTFILE=smrd0 $(M4OPTS) > smrd0.f
	$(FC) $(FFLAGSS) smrd0.f

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
mrdci5.o:	mrdci5.m
	cat ../utilities/gener.m mrdci5.m | $(M4) $(M4OPTS) > mrdci5.f
	$(FC) $(FFLAGSS) mrdci5.f
#
rpa.o:	rpa.m
	cat ../utilities/gener.m rpa.m | $(M4) $(M4OPTS) > rpa.f
	$(FC) $(FFLAGSS) rpa.f
#
tdaf.o:	tdaf.m
	cat ../utilities/gener.m tdaf.m | $(M4) $(M4OPTS) > tdaf.f
	$(FC) $(FFLAGSV4) tdaf.f
#
optim.o:	optim.m
	cat ../utilities/gener.m optim.m | $(M4) $(M4OPTS) > optim.f
	$(FC) $(FFLAGSS) optim.f
#
guess.o:	guess.f
	$(FC) $(FFLAGSS) $*.f
#
casb.o:	casb.f
	$(FC) $(FFLAGSS) $*.f
#
mcscfc.o:	mcscfc.f
	$(FC) $(FFLAGSS) $*.f
#
intege.o:	intege.f
	$(FC) $(FFLAGSS) $*.f

newmrd2.o:	newmrd2.f
	$(FC) $(FFLAGSS) $*.f
#
newmrd5.o:	newmrd5.f
	$(FC) $(FFLAGSS) $*.f
#
newmrd6.o:	newmrd6.f
	$(FC) $(FFLAGSS) $*.f
#
# Build Information
#
#  uname -a = SunOS frontend 5.9 Generic_112233-12 sun4u sparc SUNW,Sun-Fire-15000
#  f90 -V = f90: Forte Developer 7 Fortran 95 7.0 2002/03/09
