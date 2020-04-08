#doc Machine-dependant file for GAMESS-UK serial build on Linux/Pentium
#doc with PGI compilers version 7.1.6
#doc
#doc Options:
#doc acml   - use pgi blas library in /opt/pgi/linux86/6.2/lib
#doc mopac  - include mopac code in the build
#doc debug  - include debugging information in the objects (no optimsation)
#doc static - build a static binary
#doc xml    - build AgentX to allow the code to read/write xml (Experimental)
#
# DEFAULT AND OPTIONAL OPTIONS
#dopt zora vb drf nbo mrdci vdw sysmo dl-find
#opt acml debug static xml mopac
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,littleendian,cio,unix,doublebackslash,upck-equiv,i8drct,extpopcnt,extleadz

#--#if charmm#
# currently CHARMM builds with pgf77 so GAMESS-UK follows suit (this may change)
FC = pgf77
LD = pgf77
#FC90 = pgf90
#--#else#
FC = pgf90
LD = pgf90
FC90 = pgf90
#--#endif#
CC=pgcc

FFLAGSTMP = -c
CFLAGSTMP = -c -DLINUX -D_FILE_OFFSET_BITS=64 
LDFLAGSTMP=-Wl,-Map,link.map

RANLIB=ar -s

#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g -Ktrap=fp
FFLAGSS  = ${FFLAGSTMP} -g -Ktrap=fp
FFLAGSN  = ${FFLAGSTMP} -g -Ktrap=fp
FFLAGSN0 = ${FFLAGSTMP} -g -Ktrap=fp
CFLAGS   = ${CFLAGSTMP} -g
LDDEBUG=-g -Ktrap=fp 
#--#else# 
FFLAGSV  = ${FFLAGSTMP} -fastsse
FFLAGSS  = ${FFLAGSTMP} -O
FFLAGSN  = ${FFLAGSTMP} -O1
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS   = ${CFLAGSTMP}
LDDEBUG=
#--#endif#

#--#if static#
LDFLAGSSTATIC= -Bstatic -g77libs
#--#endif#
LDFLAGS  =  $(LDFLAGSTMP) $(LDDEBUG) $(LDFLAGSSTATIC)

#--#if acml#
LIBBLAS=-L/opt/pgi/linux86/6.2/lib -lacml
# Below causes segmentation faults
#lBLAS= -lacml -lpgsse2
#--#endif acml#

OPTIONS=${MACHOPT}

# Bring all the options together
BL_LIB = ${LIBBLAS}

#
# ===============  Additional Files 
#
EXTRA_MP2=  check0a.o
#
# ===============  Compiler Exceptions
#
check0a.o:	newmrd5.m
	cat ../utilities/gener.m newmrd5.m | $(M4) -DGEN_EXTRACTFILE=check0a $(M4OPTS) > check0a.f
	$(FC) $(FFLAGSN) check0a.f
