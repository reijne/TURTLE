#doc Machine-dependent file for GAMESS-UK Global Array based build on
#doc Pentium/Xeon processors with GNU compilers
#doc
#doc If on linking you see errors due to the ma libraries, please edit this 
#doc file to set GA_VERSION_PAR to include FOPT_REN=-fno-second-underscore
#doc
#doc Options:
#doc mpiwrap - use mpi wrappers ( e.g. mpif77 )to link in the mpi libraries
#doc           build currently untested without this option
#doc
#doc i8  - build with integer*8. This is only available with the GA build and
#doc       implies using the tcgmsg-mpi Global Array library

#
# DEFAULT AND OPTIONAL OPTIONS
#dopt ga mpi i8 mp2 zora vb peigs mpiwrap f77
#opt debug static
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,bits8,glibc

#--#if mpiwrap#
FC   = mpif77
LD   = mpif77
FC90 = mpif90
CC   = mpicc
CXX  = mpicxx
#--#else#
FC   = g77
LD   = g77
FC90 = g77
CC   = gcc
CXX  = g++
#--#endif mpiwrap#

#--#if profile#
FFLAGSTMP = -pg -c -fno-second-underscore -fno-globals
#--#else#
FFLAGSTMP = -c -fno-second-underscore -fno-globals
#--#endif#

#--#if mpi#
#----#if mpiwrap#
#----#else#
FINCLUDE = -I${MPI_INCLUDE}
CINCLUDE = -I${MPI_INCLUDE}
#----#endif mpiwrap#
#--#endif mpi#

#--#if debug#
FFLAGSV  = ${FFLAGSTMP} ${FINCLUDE} -g
FFLAGSS  = ${FFLAGSTMP} ${FINCLUDE} -g
FFLAGSN  = ${FFLAGSTMP} ${FINCLUDE} -g
FFLAGSN0 = ${FFLAGSTMP} ${FINCLUDE} -g
CFLAGS   = ${CINCLUDE} -g -c -DLINUX -D_FILE_OFFSET_BITS=64
CXXFLAGS = ${CINCLUDE} -g -c 
#--#else# 
FFLAGSV  = ${FFLAGSTMP} ${FINCLUDE} -O3 -malign-double
FFLAGSS  = ${FFLAGSTMP} ${FINCLUDE} -O  -malign-double
FFLAGSN  = ${FFLAGSTMP} ${FINCLUDE} -O1 -malign-double
FFLAGSN0 = ${FFLAGSTMP} ${FINCLUDE}  
CFLAGS   = ${CINCLUDE} -c -DLINUX -D_FILE_OFFSET_BITS=64
CXXFLAGS = ${CINCLUDE} -c -O2
#--#endif#
# Need -fno-globals for drf/drfsub.f or else g77 crashes without warning

#--#if profile#
#----#if static#
LDFLAGS  = -pg -static -static-libgcc -D_FILE_OFFSET_BITS=64
#----#else#
LDFLAGS  = -pg -D_FILE_OFFSET_BITS=64
#----#endif#
#--#else#
#----#if static#
LDFLAGS  = -static -static-libgcc -D_FILE_OFFSET_BITS=64
#----#else#
LDFLAGS  = -D_FILE_OFFSET_BITS=64
#----#endif#
#--#endif#

# Global Array-based & MPI variables
GA_F77_DEFS = -traditional
GA_TARGET=LINUX
PEIGS_TARGET=LINUX

#--#if mpiwrap#
# We don't need to pass through any mpi variables to the GA's as the mpiwrappers take care of it
LIBMPI=
# Below may be required if there are link errors due to the ma tools etc being compiled with 2 underscores
#GA_VERSION_PAR= GA_VERSION=SHMEM FOPT_REN=-fno-second-underscore
#----#if scalapack#
GA_VERSION_PAR= GA_VERSION=SHMEM USE_SCALAPACK=yes
#----#else#
GA_VERSION_PAR= GA_VERSION=SHMEM 
#----#endif#
#--#else#
MPI_INCLUDE = /usr/local/include
MPI_LIB     = /usr/local/lib
LIBMPI      = -L$(MPI_LIB) -lmpich -lfmpich -lpthread
#----#if scalapack#
GA_VERSION_PAR= GA_VERSION=SHMEM USE_SCALAPACK=yes MPI_INCLUDE=${MPI_INCLUDE}  MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#----#else#
GA_VERSION_PAR= GA_VERSION=SHMEM MPI_INCLUDE=${MPI_INCLUDE}  MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#----#endif#
#--#endif mpiwrap#

#--#if scalapack#
LIBBLAS=-L/usr/local/lib  -lscalapack -lblacsF77init -lblacsCinit -lblacs -llapack -lblas
BLASOPT=,blas,scalapack
#--#endif#
OPTIONS=${MACHOPT}${BLASOPT}
BL_LIB= ${LIBBLAS} ${LIBMPI}

# Diesel build options 
LD_DIESEL = ${CXX}
DIESEL_LIBS = -lg2c
FLEX = /usr/bin/flex
__UNDERBAR = 1
SIZEOF_VOID_P = 4
LONG_LONG_INT = long long int
LONG_INT = long int
INT = int
SHORT_INT = short int

#
# ===============  Additional Files 
#
EXTRA=gethes.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o aprq34d.o mcdab_ga.o 

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

aprq34d.o:	mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=aprq34d $(M4OPTS) > aprq34d.f
	$(FC) $(FFLAGSS) aprq34d.f

# Information on the machine this build was configured on
#
# uname -a :
# Linux ccp1.dl.ac.uk 2.4.20-24.7smp #1 SMP Mon Dec 1 13:18:03 EST 2003 i686 unknown
# 
# cat /proc/cpuinfo :
# processor	: 0
# vendor_id	: GenuineIntel
# cpu family	: 15
# model		: 2
# model name	: Intel(R) Xeon(TM) CPU 2.40GHz
# stepping	: 9
# cpu MHz		: 2399.365
# cache size	: 512 KB
# physical id	: 0
# siblings	: 1
# fdiv_bug	: no
# hlt_bug		: no
# f00f_bug	: no
# coma_bug	: no
# fpu		: yes
# fpu_exception	: yes
# cpuid level	: 2
# wp		: yes
# flags		: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm
# bogomips	: 4784.12
# 
# mpif77 -compiler gnu -show :
# /opt/score/bin/scoref77 -compiler=gnu -compiler=gnu -L/opt/score/mpi/mpich-1.2.5/i386-redhat7-linux2_4_gnu/lib -lmpich
# 
# mpif77 -compiler gnu --version :
# GNU Fortran 0.5.26 20000731 (Red Hat Linux 7.3 2.96-113)
# Copyright (C) 1997 Free Software Foundation, Inc.
# 
# mpicc -compiler gnu -show :
# /opt/score/bin/scorecc -compiler=gnu -DUSE_STDARG -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 -DUSE_STDARG=1 -DMALLOC_RET_VOID=1 -compiler=gnu -L/opt/score/mpi/mpich-1.2.5/i386-redhat7-linux2_4_gnu/lib -lmpich
# 
# mpicc -compiler gnu --version :
# 2.96
# 
# Score system was 5.8.2
