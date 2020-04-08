#doc Machine-dependent file for GAMESS-UK Global Array based build on
#doc Pentium/Xeon processors with GNU compilers (gfortran or g95 & gcc 4.X)
#doc
#doc Options:
#doc mpiwrap - use mpi wrappers ( e.g. mpif77 ) to link in the mpi libraries
#doc         - mpi without wrappers is for default LAM-mpi on gfortran installation
#doc i8  - build with integer*8. This is only available with the GA build and 
#doc       implies using the tcgmsg-mpi Global Array library
#doc
#     static option has been disabled (for static executable rebuild lam-mpi with:
#     configure --without-memory-manager --enable-mpi-threads --enable-progress-threads 
#                  --enable-static --disable-shared --enable-mca-static --with-devel-headers 
#                   and relink the gamess o's and a's with -static )
#doc
# DEFAULT AND OPTIONAL OPTIONS
#dopt ga mpi mp2 nbo vb vdw zora newscf dl-find mopac masscf nomalign peigs 
#opt debug i8 scalapack mpiwrap
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,linux64,littleendian,cio,unix,doublebackslash,upck-equiv,macosx

#--#if mpiwrap#
FC   = mpif90
LD   = mpif90
FC90 = mpif90
CC   = mpicc
CXX  = mpicxx
#--#else#
FC   = gfortran -I/usr/local/include
LD   = gfortran
FC90 = gfortran
# CC   = gcc -I/usr/local/include -D_REENTRANT -L/usr/local/lib -llammpio -llammpi++ -llamf77mpi -lmpi -llam -lutil -ldl 
CC   = gcc -I/usr/local/include -D_REENTRANT  
CXX  = g++
#--#endif mpiwrap#

#--#if i8#
FFLAGSI8=-fdefault-integer-8
CFLAGSI8=-DLONG_LONG_INTEGER -DEXT_INT
I8_M4_OPTS=,i8
#--#else#
FFLAGSI8=
CFLAGSI8=-DSTD_INT
I8_M4_OPTS=,bits8
#--#endif#

## Gfortran >4.2.0 has now deliberately implemented the flawed F95/F2003 standard for BOZ (accidentally
## this standard demands an intermediate conversion to a strictly positive signed integer)
## As a result z'7*' is now the highest bitpatten (7 bits in integer*1, 15 bits in integer*2, 31 bits in integer*4 etc.).
## As a result b'0*' is now the highest bitpatten (0 bits in a single bit etc.).
## The -fno-range-check flag will disable all range checking. 
FFLAGSTMP = -c -fno-second-underscore -fno-range-check  -DGFS $(FFLAGSI8)
CFLAGSTMP = -c -DLINUX -D_FILE_OFFSET_BITS=64 -DGFS $(CFLAGSI8)

#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g
FFLAGSS  = ${FFLAGSTMP} -g
FFLAGSN  = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS   = -g ${CFLAGSTMP}
CXXFLAGS = -g -c  
#--#elseif nomalign# 
FFLAGSV  = ${FFLAGSTMP} -O2 
FFLAGSS  = ${FFLAGSTMP} -O 
FFLAGSN  = ${FFLAGSTMP} -O1
FFLAGSN0 = ${FFLAGSTMP} -O0
CFLAGS   =  ${CFLAGSTMP}
CXXFLAGS = -c -O2 
#--#else#
FFLAGSV  = ${FFLAGSTMP} -O2 -malign-double
FFLAGSS  = ${FFLAGSTMP} -O  -malign-double
FFLAGSN  = ${FFLAGSTMP} -O1 -malign-double
FFLAGSN0 = ${FFLAGSTMP} -O0 -malign-double
CFLAGS   = -c ${CFLAGSTMP}
CXXFLAGS = -c -O2 
#--#endif#

#--#if static#
LDFLAGSTMP  = -static -static-libgcc 
#--#elseif nomalign#
LDFLAGSTMP  =
#--#else#
LDFLAGSTMP  =  -malign-double
#--#endif#

#--#if debug#
LDFLAGS  = ${LDFLAGSTMP} -g
#--#else
LDFLAGS  = ${LDFLAGSTMP} -O2 
#--#endif#


# Global Array-based & MPI variables
GMAKE=/usr/bin/make
GA_F77_DEFS = -traditional
GA_TARGET=MACX
#PEIGS_TARGET=LINUX
PEIGS_TARGET=X86OSX

#--#if mpiwrap#
# We don't need to pass through any mpi variables to the GA's as the mpiwrappers take care of it
GA_VERSION_PAR= GA_VERSION=SHMEM 
#--#else#
MPI_INCLUDE=/usr/local/include
MPI_LIB=/usr/local/lib
LIBMPI= -L$(MPI_LIB) -lmpi_f90 -lmpi_f77 -lmpi
GA_VERSION_PAR= GA_VERSION=SHMEM MPI_INCLUDE=${MPI_INCLUDE}  MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
#--#endif mpiwrap#

#--#if scalapack#
LIB_SCALAPACK=/Users/jmht/work/codes/scalapack_installer_0.91/lib/libscalapack.a \
/Users/jmht/work/codes/scalapack_installer_0.91/lib/blacsF77.a \
/Users/jmht/work/codes/scalapack_installer_0.91/lib/blacs.a 
#--#endif scalapack#

# The m4 pre-processing directives
OPTIONS=${MACHOPT}$(I8_M4_OPTS)
BL_LIB= ${LIBMPI} $(LIB_SCALAPACK)

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
EXTRA = gethes.o 
EXTRA_MP2 = mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o aprq34d.o mcdab_ga.o tst1s.o
EXTRA_CI = 

#
# ===============  Compiler Exceptions
#

util5.o:	util5.f
	$(FC) $(FFLAGSN) -fno-range-check $*.f 

mainci.o:	mainci.f
	$(FC) $(FFLAGSN) -fno-range-check $*.f 

machscf.o:	machscf.f
	$(FC) $(FFLAGSN0) -fno-range-check $*.f 

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
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
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

