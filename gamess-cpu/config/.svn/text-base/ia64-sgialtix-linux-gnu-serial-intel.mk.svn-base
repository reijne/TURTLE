#doc Makefile settings for the GAMESS-UK serial build on the SGI Altix
#doc running Linux on Itanium 2 processors.
#doc
#doc This file has been tested with Intel compiler versions 9.1 and 10.1
#doc It has not been tested with 10.0, but we have found this version to
#doc be a disaster for GAMESS-UK so use it at your peril...
#doc
#doc NB!!! Make sure that you are linking to the Intel version of libimf.so
#doc       and not the SGI one from /usr/lib/sgi/libimf.so.6 as this will
#doc       cause the code to seg fault.
#doc
#doc Options:
#doc i8  - use integer*8 (64-bit fortran integers)
#doc mkl - use Intel MKL library (default location: /opt/modules/intel/mkl72/lib/64)
#doc scs - use SGI scs numerical library
#
#dopt drf mrdci zora nbo vb  sysmo vdw dl-find scs i8
#opt  mkl mopac
#
#
# ================ M4 Processing options
#

# Machine-specific options 
MACHINE_KEY=G
MACHOPT=itanium,linux,littleendian,cio,unix,upck-equiv,altix

IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)

#--#if mkl#
LIBBLAS = -L/opt/modules/intel/mkl721/lib/64  -lmkl_ipf -lguide
BLASOPT=,blas
#--#elseif scs#
#---#if i8#
LIBBLAS=-lscs_i8
#---#else#
LIBBLAS=-lscs
#---#endif i8#
BLASOPT=,blas
#--#endif mkl scs#

# Put it all together
BL_LIB= ${MPI_LIBS} ${LIBBLAS}

## M4 Options ##
#--#if i8#
OPTIONS=${MACHOPT}${BLASOPT},i8
#--#else#
OPTIONS=${MACHOPT}${BLASOPT},i8drct,64bitpointers
#--#endif i8#

#
# ===============  Compiler Options
#
FC = ifort
LD = ifort
CC = icc -no-gcc
FC90 = ifort
CXX = icc
RANLIB = ranlib

#--#if i8#
FFLAGI8= -i8
CFLAGI8= -DLONG_INTEGER
#--#else#
FFLAGI8=
CFLAGI8=
#--#endif i8#
FFLAGSTMP = -c  ${FFLAGI8}

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS = -g -c ${CFLAGI8}
LDFLAGS  = -g
CXXFLAGS = -g -c
#--#else# 
FFLAGSV = ${FFLAGSTMP} -O2 -cm -w90 -w95 -ip -ftz
FFLAGSS = ${FFLAGSTMP} -O2 -cm -w90 -w95  -ftz
FFLAGSN = ${FFLAGSTMP} -O1 -cm -w90 -w95  -ftz
FFLAGSN0 = ${FFLAGSTMP}  -O0
CFLAGS  = -c   ${CFLAGI8}
LDFLAGS =  
CXXFLAGS = -O2 -c
#--#endif#

# Diesel build options
LD_DIESEL = icpc
#--#if mkl scs#
DIESEL_LIBS = ${LBLAS} ${lBLAS} -lifcore
#--#else#
DIESEL_LIBS = -lifcore
#--#endif#
FLEX = /usr/bin/flex
__UNDERBAR = 1
SIZEOF_VOID_P = 8
#--#if i8#
LONG_LONG_INT = long int
LONG_INT = long int
INT = long int
SHORT_INT = short int
#--#else#
LONG_LONG_INT = long long int
LONG_INT = long int
INT = int
SHORT_INT = short int
#--#endif i8#

#
# ===============  Additional Files 
#
#EXTRA=gethes.o
EXTRA_MP2= tst1s.o check0a.o
#
# ===============  Compiler Exceptions
#
# check0a and tst1s are extracted by the "linux" m4 option so it's not clear if these are
# actually needed or not
#
check0a.o:	newmrd5.m
	cat ../utilities/gener.m newmrd5.m | $(M4) -DGEN_EXTRACTFILE=check0a $(M4OPTS) > check0a.f
	$(FC) $(FFLAGSN) check0a.f

tst1s.o:	direct.m
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f

### Below exceptions required for Intel 10.1 
# Compiling with o2 with intel 10.1 causes numerous errors
integs.o: integs.f
	$(FC) $(FFLAGSN) integs.f

# With the -ip flag included (e.g.) c2021_a fails with
# "forrtl: severe (71): integer divide by zero" under intel 10.1
mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSS) $*.f

# With the -ip flag included c2031 fails with intel 10.1
newmrd2.o:	newmrd2.f
	$(FC) $(FFLAGSS) $*.f

# Without -O0, c2024_a gives incorrect answers with intel 10.1
basis.o:	basis.f
	$(FC) $(FFLAGSN0) $*.f
### End Intel 10.1

# The below are a list of all the exceptions we've ever included for the intel compilers
# but which aren't needed currently (for 9.1 or 10.1). They are kept in case a future 
# release of the intel compiler breaks them again.
#
# With O3, a large number of the hessian jobs fail in serial ( chap2/c2016* )
#cphf.o: cphf.f
#	$(FC) $(FFLAGSS) $*.f

#dircta.o:	dircta.f
#	$(FC) $(FFLAGSN0) $*.f

# With O3 direct mp2 job c2009_f (2nd job with scftype direct mp2) fails with 
# the scf generating lots of NaN
#direct.o: direct.f
#	$(FC) $(FFLAGSS) $*.f

#eispack.o:	eispack.f
#	$(FC) $(FFLAGSN) $*.f


#gethes.o:	casb.m
#	cat  ../utilities/gener.m casb.m | $(M4) -DGEN_EXTRACTFILE=gethes $(M4OPTS) > gethes.f
#	$(FC) $(FFLAGSN) gethes.f

# If compiled with O3 mclr spends it's entire time swapping in and out of memory
# and eventually dies with:  Fatal compilation error: Out of memory asking for 69632.
#mclr.o: mclr.f
#	$(FC) $(FFLAGSS) $*.f

#mkmakw.o:	drvmp.m
#	cat  ../utilities/gener.m drvmp.m | $(M4) -DGEN_EXTRACTFILE=mkmakw $(M4OPTS) > mkmakw.f
#	$(FC) $(FFLAGSS) mkmakw.f

#mpmakw.o:	secmp2.m
#	cat  ../utilities/gener.m secmp2.m | $(M4) -DGEN_EXTRACTFILE=mpmakw $(M4OPTS) > mpmakw.f
#	$(FC) $(FFLAGSS) mpmakw.f

#mrdci1.o:	mrdci1.f
#	$(FC) $(FFLAGSN0) $*.f

# With O3 c2021_f test 1 and c2021_g and h fail to validate
#mrdci3.o: mrdci3.f
#	$(FC) $(FFLAGSS) $*.f

#mrdci5.o:	mrdci5.f
#	$(FC) $(FFLAGSN) $*.f

# With O3, all the mrdci test jobs c2021a-i all fail to validate
#mrdci6.o: mrdci6.f
#	$(FC) $(FFLAGSS) $*.f

# With O3 c2021_j,l and o fail to validate
#newmrd3.o: newmrd3.f
#	$(FC) $(FFLAGSS) $*.f

# optim.o included as a number of the optimisation jobs on pauling failed or
# required an excessive number of optimisation steps to converge
# optim.o: optim.f
#	$(FC) $(FFLAGSN) $*.f


#rpa.o:	rpa.f
#	$(FC) $(FFLAGSN) $*.f

# With O3 all the ccsd jobs (c2029*) fail to validate
#tsort.o: tsort.f
#	$(FC) $(FFLAGSS) $*.f


#umpe3a.o:	mp3.m
#	cat  ../utilities/gener.m mp3.m | $(M4) -DGEN_EXTRACTFILE=umpe3a $(M4OPTS) > umpe3a.f
#	$(FC) $(FFLAGSS) umpe3a.f

#umpe3b.o:	mp3.m
#	cat  ../utilities/gener.m mp3.m | $(M4) -DGEN_EXTRACTFILE=umpe3b $(M4OPTS) > umpe3b.f
#	$(FC) $(FFLAGSS) umpe3b.f

#vbin.o:	vbin.f
#	$(FC) $(FFLAGSN) $*.f
