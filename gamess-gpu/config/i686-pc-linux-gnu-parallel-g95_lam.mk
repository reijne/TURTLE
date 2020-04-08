#doc Machine-dependent file for GAMESS-UK Global Array based build on
#doc Pentium/Xeon processors with GNU compilers (gfortran or g95 & gcc 4.X)
#doc
#doc Options:
#doc mpiwrap - use mpi wrappers ( e.g. mpif77 ) to link in the mpi libraries
#doc         - mpi without wrappers is for default LAM-mpi on gfortran installation
#doc i8  - build with integer*8. This is only available with the GA build and
#doc       implies using the tcgmsg-mpi Global Array library
#doc
#doc static option has been disabled (for static executable rebuild lam-mpi with:
#doc configure --without-memory-manager --enable-mpi-threads --enable-progress-threads 
#doc           --enable-static --disable-shared --enable-mca-static --with-devel-headers 
#doc           and relink the gamess o's and a's with -static )
#doc
# DEFAULT AND OPTIONAL OPTIONS
#dopt ga mpi i8 mp2 zora vb peigs mpiwrap nbo vdw
#opt mopac debug 
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,bits8

#--#if mpiwrap#
FC   = mpif77
LD   = mpif77
FC90 = mpif77
CC   = mpicc
CXX  = mpicxx
#--#else#
FC   = gfortran
LD   = gfortran
FC90 = gfortran
CC   = gcc
CXX  = g++
#--#endif mpiwrap#

## Gfortran >4.2.0 has now deliberately implemented the flawed F95/F2003 standard for BOZ (accidentally
## this standard demands an intermediate conversion to a strictly positive signed integer)
## As a result z'7*' is now the highest bitpatten (7 bits in integer*1, 15 bits in integer*2, 31 bits in integer*4 etc.).
## As a result b'0*' is now the highest bitpatten (0 bits in a single bit etc.).
## The -fno-range-check flag will disable all range checking. 
FFLAGSTMP = -c -fno-second-underscore  -DGFS

#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g
FFLAGSS  = ${FFLAGSTMP} -g
FFLAGSN  = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS   = -g -c -DLINUX -D_FILE_OFFSET_BITS=64 -DGFS
CXXFLAGS = -g -c 
#--#else# 
FFLAGSV  = ${FFLAGSTMP} -O2 -malign-double
FFLAGSS  = ${FFLAGSTMP} -O  -malign-double
FFLAGSN  = ${FFLAGSTMP} -O1 -malign-double
FFLAGSN0 = ${FFLAGSTMP} -O0 -malign-double
CFLAGS   = -c -DLINUX -D_FILE_OFFSET_BITS=64 -DGFS
CXXFLAGS = -c -O2
#--#endif#
# Need -fno-globals for drf/drfsub.f or else g77 crashes without warning

#--#if static#
LDFLAGSTMP  = -static -static-libgcc -D_FILE_OFFSET_BITS=64
#--#else#
LDFLAGSTMP  = -D_FILE_OFFSET_BITS=64 
#--#endif#

#--#if debug#
LDFLAGS  = ${LDFLAGSTMP} -g
#--#else
LDFLAGS  = ${LDFLAGSTMP} -O2 -malign-double
#--#endif#


# Global Array-based & MPI variables
GA_F77_DEFS = -traditional
GA_TARGET=LINUX

#--#if mpiwrap#
# We don't need to pass through any mpi variables to the GA's as the mpiwrappers take care of it
LIBMPI=
# Below may be required if there are link errors due to the ma tools etc being compiled with 2 underscores
#GA_VERSION_PAR= GA_VERSION=SHMEM FOPT_REN=-fno-second-underscore
GA_VERSION_PAR= GA_VERSION=SHMEM
#--#else#
MPI_INCLUDE=/usr/include/lam
MPI_LIB=/usr/lib/lam
LIBMPI= -L$(MPI_LIB) -lmpi
GA_VERSION_PAR= GA_VERSION=SHMEM MPI_INCLUDE=${MPI_INCLUDE}  MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#--#endif mpiwrap#

NUMERIC_OBJ=dblas.o iblas.o dlapack.o dlamch.o mxm_noblas3.o eispack.o linpack.o
OPTIONS=${MACHOPT}
BL_LIB= ${LIBMPI}

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

dircta.o:	dircta.f
	$(FC) $(FFLAGSN0) $*.f 

rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

aprq34d.o:	mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=aprq34d $(M4OPTS) > aprq34d.f
	$(FC) $(FFLAGSS) aprq34d.f

# Information on the machine this build was configured on (Marc's Notebook)
#
# uname -a :
# Linux pctje.xs4all.nl 2.6.18-8.1.6.el5 #1 SMP Thu Jun 14 11:57:17 EDT 2007 i686 i686 i386 GNU/Linux
# 
# cat /proc/cpuinfo :
#processor       : 0
#vendor_id       : GenuineIntel
#cpu family      : 6
#model           : 15
#model name      : Intel(R) Core(TM)2 Duo CPU     T7300  @ 2.00GHz
#stepping        : 10
#cpu MHz         : 2001.000
#cache size      : 4096 KB
#physical id     : 0
#siblings        : 2
#core id         : 0
#cpu cores       : 2
#fdiv_bug        : no
#hlt_bug         : no
#f00f_bug        : no
#coma_bug        : no
#fpu             : yes
#fpu_exception   : yes
#cpuid level     : 10
#wp              : yes
#flags           : fpu vme de pse tsc msr pae mce cx8 apic mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe lm constant_tsc pni monitor ds_cpl vmx est tm2 cx16 xtpr lahf_lm
#bogomips        : 3993.08
#
#processor       : 1
#vendor_id       : GenuineIntel
#cpu family      : 6
#model           : 15
#model name      : Intel(R) Core(TM)2 Duo CPU     T7300  @ 2.00GHz
#stepping        : 10
#cpu MHz         : 2001.000
#cache size      : 4096 KB
#physical id     : 0
#siblings        : 2
#core id         : 1
#cpu cores       : 2
#fdiv_bug        : no
#hlt_bug         : no
#f00f_bug        : no
#coma_bug        : no
#fpu             : yes
#fpu_exception   : yes
#cpuid level     : 10
#wp              : yes
#flags           : fpu vme de pse tsc msr pae mce cx8 apic mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe lm constant_tsc pni monitor ds_cpl vmx est tm2 cx16 xtpr lahf_lm
#bogomips        : 4149.61
# 
# mpif77 -compiler gnu -show :
# f95 -I/usr/include/lam -I/usr/include/lam/32 -m32 -pthread -compiler gnu -L/usr/lib/lam -llammpio -llamf77mpi -lmpi -llam -laio -laio -lutil -ldl
# 
# mpif77 -compiler gnu --version :
# GNU Fortran 95 (GCC) 4.1.1 20070105 (Red Hat 4.1.1-52)
# Note: On this platform f95 is a symbolic link to gfortran
#
# mpicc -compiler gnu -show :
# gcc -I/usr/include/lam -I/usr/include/lam/32 -m32 -pthread -compiler gnu -L/usr/lib/lam -llammpio -llamf77mpi -lmpi -llam -laio -laio -lutil -ldl
# 
# mpicc -compiler gnu --version :
# gcc (GCC) 4.1.1 20070105 (Red Hat 4.1.1-52)
# 
