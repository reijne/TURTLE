#doc  This is the machine-specific file for the serial version of GAMESS-UK
#doc  using  the Intel compiler 8.1 with Suse 9.3 running on Opteron processors
#doc
#doc Options:
#doc i8 - build the code with 64bit integers
#doc goto - link against the GOTO blas library 
#doc mkl - link against the Intel MKL library 
#
#dopt vb nbo mrdci zora drf vdw sysmo dl-find i8
#opt goto mkl static mopac
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,upck-equiv,opteron

# Bitwise operators that differ from the defaults.
IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)

#
# ===============  Compiler Options
#

## INTEL COMPILER VERSION 8.1 ##
FC = ifort
LD = ifort
FC90 = ifort
CC= icc
RANLIB=ranlib

#--#if i8#
FFLAGSI8=-i8
CFLAGSI8=-DLONG_INTEGER
# Added linux64 to get tst1s extracted from direct.m
I8_M4_OPT=,i8,linux64
#--#else#
FFLAGSI8=
CFLAGSI8=-D_FILE_OFFSET_BITS=64
I8_M4_OPT=,64bitpointers
#--#endif#

FFLAGSTMP = -c -ftz ${FFLAGSI8}
CFLAGSTMP= -c -no-gcc  -DLINUX ${CFLAGSI8}

#--#if static#
LDSTATIC  = -i-static
#--#endif#
LDFLAGSTMP  = -Wl,-Map load.map ${LDSTATIC}

#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g
FFLAGSS  = ${FFLAGSTMP} -g
FFLAGSN  = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS   = ${CFLAGSTMP} -g 
LDFLAGS  = ${LDFLAGSTMP} -g
#--#else# 
#  note that -O3 for FFLAGSV gives a no. of errors
#  when running chap3 - this should be OK based on
#  Xeon 32-bit, but not here .. maybe I8 related
FFLAGSV  = ${FFLAGSTMP}  -ip -cm -w90 -w95 -O2 -vec_report0
FFLAGSS  = ${FFLAGSTMP}  -ip -cm -w90 -w95 -O2 -vec_report0
FFLAGSN  = ${FFLAGSTMP}  -ip -cm -w90 -w95 -O1 -vec_report0
FFLAGSN0 = ${FFLAGSTMP} -cm -w90 -w95 -O0
CFLAGS   = ${CFLAGSTMP}
LDFLAGS  = ${LDFLAGSTMP} 
#--#endif#


### NOTE ####
# This version of GAMESS-UK is built with i8 by default.
# It is only fortuitous that linking against numerical 
# libraries works with i8 (the significant bits are at the correct end).
# dgthr and dsctr take integer array arguments and so  the library 
# versions crash when used with i8. So, when using blas with i8, 
# we need to pull in our own versions of dsctr, dgthr & ddoti from m4/util3.m

#--#if goto#
LIBBLAS = -L/usr/local/lib64 -lblas
BLASOPT=,blas
#--#elseif mkl#
LIBBLAS=-L/opt/intel/compiler81/fc/lib -lmkl_em64t -lguide -L/usr/lib64 -lpthread
BLASOPT=,blas
#--#endif goto#

# Bring all the options together
OPTIONS=${MACHOPT}${I8_M4_OPT}${BLASOPT}

BL_LIB = ${LIBBLAS} 

#
# ===============  Additional Files 
#
#--#if i8#
EXTRAI8=tst1s.o
#--#endif#
EXTRA=gethes.o ${EXTRAI8}
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o
EXTRA_DFT = jkint_dft.o jkder_dft.o

#
# ===============  Compiler Exceptions
#
tst1s.o:	direct.m
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f 

gethes.o:	casb.m
	cat  ../utilities/gener.m casb.m | $(M4) -DGEN_EXTRACTFILE=gethes $(M4OPTS) > gethes.f
	$(FC) $(FFLAGSN) gethes.f

mkmakw.o:	drvmp.m
	cat  ../utilities/gener.m drvmp.m | $(M4) -DGEN_EXTRACTFILE=mkmakw $(M4OPTS) > mkmakw.f
	$(FC) $(FFLAGSS) mkmakw.f

mpmakw.o:	secmp2.m
	cat  ../utilities/gener.m secmp2.m | $(M4) -DGEN_EXTRACTFILE=mpmakw $(M4OPTS) > mpmakw.f
	$(FC) $(FFLAGSS) mpmakw.f

umpe3a.o:	mp3.m
	cat  ../utilities/gener.m mp3.m | $(M4) -DGEN_EXTRACTFILE=umpe3a $(M4OPTS) > umpe3a.f
	$(FC) $(FFLAGSS) umpe3a.f

umpe3b.o:	mp3.m
	cat  ../utilities/gener.m mp3.m | $(M4) -DGEN_EXTRACTFILE=umpe3b $(M4OPTS) > umpe3b.f
	$(FC) $(FFLAGSS) umpe3b.f

check0a.o:	newmrd5.m
	cat ../utilities/gener.m newmrd5.m | $(M4) -DGEN_EXTRACTFILE=check0a $(M4OPTS) > check0a.f
	$(FC) $(FFLAGSN) check0a.f

rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

eispack.o:	eispack.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

# optim.o included as a number of the optimisation jobs on pauling failed or
# required an excessive number of optimisation steps to converge
optim.o: optim.f
	$(FC) $(FFLAGSN) $*.f
#$(FC) $(FFLAGSS) $*.f

# If compiled with O3 mclr spends it's entire time swapping in and out of memory
# and eventually dies with:  Fatal compilation error: Out of memory asking for 69632.
mclr.o: mclr.f
	$(FC) $(FFLAGSS) $*.f

# With O3, a large number of the hessian jobs fail in serial ( chap2/c2016* )
cphf.o: cphf.f
	$(FC) $(FFLAGSS) $*.f

# With O3, all the mrdci test jobs c2021a-i all fail to validate
mrdci6.o: mrdci6.f
	$(FC) $(FFLAGSS) $*.f

# With O3 c2021_f test 1 and c2021_g and h fail to validate
mrdci3.o: mrdci3.f
	$(FC) $(FFLAGSS) $*.f

# With O3 c2021_j,l and o fail to validate
newmrd3.o: newmrd3.f
	$(FC) $(FFLAGSS) $*.f

# With O3 all the ccsd jobs (c2029*) fail to validate
tsort.o: tsort.f
	$(FC) $(FFLAGSS) $*.f

# With O3 direct mp2 job c2009_f (2nd job with scftype direct mp2) fails with 
# the scf generating lots of NaN
direct.o: direct.f
	$(FC) $(FFLAGSS) $*.f

#
#  ========== DFT Exceptions (PGF) ===============
#
jkint_dft.o:    integ2e.m
	cat  ../utilities/gener.m  integ2e.m | $(M4) -DGEN_EXTRACTFILE=jkint_dft $(M4OPTS) $(SNGL) > jkint_dft.f
	$(FC) $(FFLAGSN) jkint_dft.f
#
jkder_dft.o: deriv2e.m
	cat ../utilities/gener.m deriv2e.m | $(M4) -DGEN_EXTRACTFILE=jkder_dft $(M4OPTS) $(SNGL) > jkder_dft.f
	$(FC) $(FFLAGSN) jkder_dft.f
#

