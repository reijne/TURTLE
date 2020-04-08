#doc Makefile settings for Serial Pentium/Linux build with Intel 7.1 compiler
#doc Default location for compiler libraries is: /usr/local/Cluster-Apps/compiler70/ia32/lib
#doc Please set the INTELLIB variable if your libraries are located elsewhere
#doc
#doc Options:
#doc mkl - use the Intel MKL blas libraries (default location: /usr/local/Cluster-Apps/mkl61/lib/32)
#
#dopt zora vb drf nbo mrdci vdw dl-find sysmo
#opt mkl debug mopac
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,i8drct

FC = ifc
FC90= ifc
CC= icc
CXX= icc
LD = ifc
FFLAGSTMP = -c
RANLIB=ranlib

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS = -g -c -DLINUX -D_FILE_OFFSET_BITS=64 
CXXFLAGS = -g -c
#--#else# 
FFLAGSV = ${FFLAGSTMP} -ip -O3
FFLAGSS = ${FFLAGSTMP} -ip -O2 
FFLAGSN = ${FFLAGSTMP} -ip -O1
FFLAGSN0 = ${FFLAGSTMP} -O0
CFLAGS  = -c -DLINUX -D_FILE_OFFSET_BITS=64
CXXFLAGS = -c -O2
#--#endif#

#--#if static#
# This should work but doesn't...
LDFLAGS  = -static -D_FILE_OFFSET_BITS=64
#--#else#
LDFLAGS  = -D_FILE_OFFSET_BITS=64
#--#endif#

#--#if mkl#
LIBBLAS=-L/usr/local/Cluster-Apps/mkl61/lib/32 -lmkl_p4 -lmkl_lapack64 -lguide
BLASOPT=,blas
#--#endif mkl#

# Set this to point at your Intel compiler libraries
INTELIB=-L/usr/local/Cluster-Apps/compiler70/ia32/lib -lPEPCF90

# Diesel build options 
LD_DIESEL = ${CXX}
#--#if mkl#
DIESEL_LIBS = -lifcore
#--#else#
DIESEL_LIBS = -llapack -lblas -lifcore
#--#endif mkl#
FLEX = /usr/bin/flex
__UNDERBAR = 1
SIZEOF_VOID_P = 4
LONG_LONG_INT = long long int
LONG_INT = long int
INT = int
SHORT_INT = short int

# Bring all the options together
OPTIONS=${MACHOPT}${BLASOPT}
BL_LIB = ${INTELIB} ${LIBBLAS}

#
# ===============  Additional Files 
#
EXTRA=gethes.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o
#
# ===============  Compiler Exceptions
#
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

tst1s.o:	direct.m
	cat ../machines/$(MACH) ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f 

rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

