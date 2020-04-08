#doc Machine-dependant file for GAMESS-UK GA build on Linux/Pentium
#doc with g95 compiler G95 (GCC 4.0.3 (g95 0.90!) Jul 27 2006)
#doc NB: Mopac is not available with this build as Real and double precision
#doc DO loop index variables are not implemented in g95.
#doc
#doc Options:
#doc mpiwrap    - use mpiwrapper scripts (e.g. mpif90)
#
# DEFAULT AND OPTIONAL OPTIONS
#dopt ga mpi i8 mp2 peigs zora vb
#opt debug static mpiwrap
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,i8drct,glibc

IAND32=iand($$1,$$2)
IOR32=ior($$1,$$2)
IXOR32=ieor($$1,$$2)
ISHFT32=ishft($$1,$$2)

#--#if mpiwrap#
# Use MPI wrapper scripts
FC = mpif90
LD = mpif90
FC90 = mpif90
CC = mpicc
CXX = mpicxx
#--#else#
FC = g95
LD = g95
FC90 = g95
CC = gcc
CXX = g++
#--#endif mpiwrap#

FFLAGSTMP = -c -fsloppy-char -fno-second-underscore
RANLIB = ranlib

# Set where to find the MPI include files
#--#if mpiwrap#
MPI_INCLUDE = 
FINCLUDE =
CINCLUDE =
#--#else#
MPI_INCLUDE = /usr/local/include
FINCLUDE = -I${MPI_INCLUDE}
CINCLUDE = -I${MPI_INCLUDE}
#--#endif mpiwrap#

# Compiler flags
#--#if debug#
FFLAGSV = ${FFLAGSTMP} ${FINCLUDE} -g 
FFLAGSS = ${FFLAGSTMP} ${FINCLUDE} -g 
FFLAGSN = ${FFLAGSTMP} ${FINCLUDE} -g 
FFLAGSN0 = ${FFLAGSTMP} ${FINCLUDE} -g
CFLAGS = ${CINCLUDE} -g -c -DLINUX -D_FILE_OFFSET_BITS=64
CXXFLAGS = ${CINCLUDE} -g -c 
LDFLAGSTMP = -D_FILE_OFFSET_BITS=64 -g
DBG_M4_OPT=,debug
#--#else# 
FFLAGSV = ${FFLAGSTMP} ${FINCLUDE} -O3
FFLAGSS = ${FFLAGSTMP} ${FINCLUDE} -O 
FFLAGSN = ${FFLAGSTMP} ${FINCLUDE} -O1
FFLAGSN0 = ${FFLAGSTMP} ${FINCLUDE}
CFLAGS  = ${CINCLUDE} -c -DLINUX -D_FILE_OFFSET_BITS=64
CXXFLAGS = ${CINCLUDE} -c -O2
LDFLAGSTMP = -D_FILE_OFFSET_BITS=64
#--#endif#
# Need -fno-globals for drf/drfsub.f or else g77 crashes without warning

#--#if static#
LDFLAGS  = ${LDFLAGSTMP} -static -static-libgcc 
#--#else#
LDFLAGS  = ${LDFLAGSTMP}
#--#endif#

# Global Array-based & MPI variables
GA_F77_DEFS = -traditional
GA_TARGET=LINUX

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
# MPI_INCLUDE SET ABOVE
MPI_LIB     = /usr/local/lib
LIBMPI      = -L$(MPI_LIB) -lmpich -lfmpich -lpthread
#----#if scalapack#
GA_VERSION_PAR= GA_VERSION=SHMEM USE_SCALAPACK=yes MPI_INCLUDE=${MPI_INCLUDE}  MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#----#else#
GA_VERSION_PAR= GA_VERSION=SHMEM MPI_INCLUDE=${MPI_INCLUDE}  MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#----#endif#
#--#endif mpiwrap#

# Peigs definitions
PEIGS_TARGET=LINUX

#--#if scalapack#
LIBBLAS=-L/usr/local/lib  -lscalapack -lblacsF77init -lblacsCinit -lblacs -llapack -lblas
BLASOPT=,blas,scalapack
#--#endif scalapack#

OPTIONS=${MACHOPT}${DBG_M4_OPT}
BL_LIB= ${LIBBLAS}${LIBMPI}

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

#
# Information on the machine this build was configured on
#
# Output of uname -a :
# Linux csevig7 2.6.11.4-21.7-smp #1 SMP Thu Jun 2 14:23:14 UTC 2005 i686 i686 i386 GNU/Linux
# 
# Output of cat /proc/cpuinfo :
# 
# processor       : 0
# vendor_id       : GenuineIntel
# cpu family      : 15
# model           : 2
# model name      : Intel(R) Pentium(R) 4 CPU 2.80GHz
# stepping        : 9
# cpu MHz         : 2793.927
# cache size      : 512 KB
# physical id     : 0
# siblings        : 2
# fdiv_bug        : no
# hlt_bug         : no
# f00f_bug        : no
# coma_bug        : no
# fpu             : yes
# fpu_exception   : yes
# cpuid level     : 2
# wp              : yes
# flags           : fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe cid xtpr
# bogomips        : 5537.79
# 
# Output of cat /etc/issue
# Welcome to SuSE Linux 9.3 (i586) - Kernel
#
# Output of /usr/bin/g95 --version
# G95 (GCC 4.0.0 20050129 (experimental) (g95!) Jun 13 2005)
# Copyright (C) 2002 Free Software Foundation, Inc.
# 
# Output of /usr/bin/gcc --version
# gcc (GCC) 3.3.5 20050117 (prerelease) (SUSE Linux)
# Copyright (C) 2003 Free Software Foundation, Inc.
