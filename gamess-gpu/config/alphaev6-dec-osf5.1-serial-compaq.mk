#doc Makefile settings for serial Compaq AlphaServer DS20E/500 under Tru64
#doc with Compaq Fortran compiler
#
#dopt zora vb drf nbo mrdci mopac vdw sysmo dl-find ksh
#opt  

#MACH is axpev6.m
MACHINE_KEY=d
MACHOPT=dec,osf,alpha,unix,blas,cio,signals,ev6,littleendian,upck-equiv,64bitpointers

# Bit-wise operators that differ from the defaults
IAND32=iand($$1,$$2)
IOR32=ior($$1,$$2)
IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)


# Compiler options
FFLAGSI8=

FC = f90
CC = cc
FFLAGSV = -fast -arch ev6 ${FFLAGSI8}  -c
FFLAGSS=  -O2 -c  ${FFLAGSI8}
FFLAGSN =   ${FFLAGSI8}
CFLAGS = ${CFLAGSI8} -O -arch ev6 -c 

# Loader
LD = $(FC) 
LDFLAGS = -fast  ${FFLAGSI8}

GMAKE = make


#GA_F77_DEFS=-DDECOSF -traditional
#GA_TARGET=DECOSF
#GA_VERSION_PAR= GA_VERSION=SHMEM
#MACH_M4_OPTS=,64bitpointers

# Objects
BL_LIB=-ldxml

#
#  ======================  Exceptions for axpev5 and axpev6 ===============
#
mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSS) $*.f

tsort.o:	tsort.f
	$(FC) $(FFLAGSV) -assume dummy_aliases $*.f
