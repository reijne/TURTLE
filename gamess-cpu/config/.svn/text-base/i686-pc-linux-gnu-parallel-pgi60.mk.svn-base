#doc This is the machine-specific file for the GLOBAL ARRAY based parallel version
#doc of GAMESS-UK running on 2GHz Xeon CPU's, compiled with PGI compilers (versions 5 & 6 )
#doc and using the mpich parallel libraries.
#doc
#doc Options:
#doc mpiwrap - use the mpixx (e.g. mpif90) wrappers
#doc           otherwise mpich is expected in: /opt/mpi/mpich-1.2.4/
#doc blas    - use the PGI blas libraries
#doc score   - build for an SCORE system
#
#dopt ga mpi peigs mp2 vb zora vdw newscf dl-find mpiwrap
#opt  blas score
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv,i8drct

#
# ===============  Compiler Options
#
## PGI COMPILER ( tested with versions 5.2 and 6.0 ) ##
#--#if mpiwrap#
FC   = mpif90
FC90 = mpif90
LD   = mpif90
CC   = mpicc
#--#elseif score#
FC   = /opt/score/bin/mpif90 -compiler pgi
FC90 = /opt/score/bin/mpif90 -compiler pgi
LD   = /opt/score/bin/mpif90 -compiler pgi
CC   = /opt/score/bin/mpicc  -compiler pgi
#--#else#
FC   = pgf90 
FC90 = pgf90
LD   = pgf90
CC   = pgcc
#--#endif#

FFLAGSTMP = -c
CFLAGSTMP =  -c -DLINUX -D_FILE_OFFSET_BITS=64
RANLIB=ranlib

#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g
FFLAGSS  = ${FFLAGSTMP} -g
FFLAGSN  = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS   = -g ${CFLAGSTMP}
#--#else# 
FFLAGSV  = ${FFLAGSTMP} -fast
FFLAGSS  = ${FFLAGSTMP} -O
FFLAGSN  = ${FFLAGSTMP} -O1
FFLAGSN0 = ${FFLAGSTMP}
CFLAGS   = -c ${CFLAGSTMP}
#--#endif#

#--#if static#
LDFLAGS  = -Bstatic -g77libs -D_FILE_OFFSET_BITS=64
#--#else#
LDFLAGS  = -D_FILE_OFFSET_BITS=64
#--#endif#

#--#if blas#
LIBBLAS=-L/opt/pgi/linux86/5.2/lib -llapack -lblas
BLASOPT=,blas
#--#endif#

#MPI and GA variables
GA_F77_DEFS = -traditional
GA_TARGET=LINUX
#--#if mpiwrap#
MPI_INCLUDE =
MPI_LIB =
LIBMPI=
#--#else#
MPI_INCLUDE = /usr/local/mpich/mpich-1.2.5.2/pgi/include
MPI_LIB = /usr/local/mpich/mpich-1.2.5.2/pgi/lib
LIBMPI=  -L${MPI_LIB} -lmpich
#--#endif mpiwrap#
GA_VERSION_PAR= GA_VERSION=SHMEM MPI_INCLUDE=${MPI_INCLUDE}  MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "

# Bring all the options together
OPTIONS=${MACHOPT}${BLASOPT}
BL_LIB = ${LIBBLAS}  ${LIBMPI}

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
# if mp2
aprq34d.o:      mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=aprq34d $(M4OPTS) > aprq34d.f
	$(FC) $(FFLAGSS) aprq34d.f

mcdab_ga.o:     mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=mcdab_ga $(M4OPTS) > mcdab_ga.f
	$(FC) $(FFLAGSS) mcdab_ga.f
# endif
