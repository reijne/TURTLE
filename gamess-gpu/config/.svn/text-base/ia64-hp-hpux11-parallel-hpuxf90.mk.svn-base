#doc  Makefile settings for the GAMESS-UK Global Array build on HP hardware running HPUX
#doc  on Itanium 2 processors and using the HP F90 v2.8.7 compiler.
#doc  This build is uses integer*8 and builds with the system veclib by default
#doc
#doc  Options:
#doc  veclib - use the HPs veclib numerical library
#doc  base   - build the minimal functionality base code.
#
#dopt ga mpi i8 peigs veclib mp2 vdw vb zora
#opt base ci
#
#
# ================ M4 Processing options
#
# Machine-specific options 
MACHINE_KEY=b
MACHOPT=hpux11,itanium,cio,unix,upck-equiv

IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IOR32=ior($$1,$$2)
IAND32=iand($$1,$$2)
IXOR64=ieor($$1,$$2)


#--#if veclib#
# Using HPs veclib numerical library
BLASOPT=,blas,veclib
LIBBLAS = /opt/mlib/lib/hpux64/libveclib8.a
#--#else#
## Needs doing so set errchk to warn the poor punter
errchk:
	@echo Your build options do not include a veclib option!
	@echo The current build requires that you include the veclib
	@echo option in your build options at configure time.
	@exit 1
#--#endif veclib#



# --------------------------- MPI/GA Flags  ----------------------------
#
MPI_L           = /opt/mpi/lib/hpux64
MPI_INC         = /opt/mpi/include
LIBMPI          = -lmpi
MPI_LIBS        = -L${MPI_L} ${LIBMPI}

# Need gmake to build the GA's
GMAKE = /usr/local/bin/gmake
GA_TARGET       = HPUX64
GA_VERSION_PAR  = GA_VERSION=SHMEM USE_MPI=YES MPI_INCLUDE=$(MPI_INC) MPI_LIB=$(MPI_L)
GA_F77_DEFS     = -D${GA_TARGET}

# Put it all together
BL_LIB= ${MPI_LIBS} ${LIBBLAS}

## M4 Options ##
#--#if i8#
OPTIONS=${MACHOPT}${BLASOPT},i8
#--#else#
OPTIONS=${MACHOPT}${BLASOPT},i8drct
#--#endif i8#


#
# ===============  Compiler Options
#
FC = f90
LD = f90
CC = cc 
FC90 = f90
CPP=/usr/ccs/lbin/cpp


# Build with i8 by default (i4 build hasn't been attempted yet).
#--#if i8#
FFLAGI8= +i8
CFLAGI8= -DLONG_INTEGER
#--#else#
FFLAGI8=
CFLAGI8=
#--#endif i8#

FFLAGSTMP = -c +ppu +DD64 +DSitanium2 ${FFLAGI8}

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g
FFLAGSS = ${FFLAGSTMP} -g
FFLAGS2 = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP}
CFLAGS = -g -c -Ae +DD64 +DSitanium2 ${CFLAGI8}
LDFLAGS  = -g
#--#else# 
FFLAGSV = ${FFLAGSTMP} +Ofast +O2 +Ofltacc=default
FFLAGSS = ${FFLAGSTMP} +O1 +Ofltacc=default +Oparmsoverlap
FFLAGS2 = ${FFLAGSTMP} +O2 +Ofltacc=default 
FFLAGSN = ${FFLAGSTMP} +O0
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS  = -c -Ae +DD64 +DSitanium2 -O ${CFLAGI8}
LDFLAGS =  +DD64 +DSitanium2 -lm +Ofast
#--#endif#

#
# ===============  Additional Files 
#
EXTRA=
EXTRA_MP2= tst1s.o
#
# ===============  Compiler Exceptions for HP 9000 (hpux-HP-UX.11.xx) ========
#
tst1s.o:	direct.m
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSS) tst1s.f
#
### Exceptions from the old Makfile ###

#
#  ==================  Exceptions for HP 9000 (hpux-HP-UX.11.xx) ========
#
# Compiler bug JAGae55172
casb.o:	casb.m
	cat ../utilities/gener.m  casb.m | $(M4) $(M4OPTS) > casb.f
	$(FC)	$(FFLAGSS) casb.f

# Compiler bug JAGae55178
newmrd1.o:	newmrd1.m
	cat ../utilities/gener.m  newmrd1.m | $(M4) $(M4OPTS) > newmrd1.f
	$(FC)	$(FFLAGSS) newmrd1.f

# Compiler bug JAGae51271
scf.o:	scf.m
	cat ../utilities/gener.m  scf.m | $(M4) $(M4OPTS) > scf.f
	$(FC)	$(FFLAGSS) scf.f

# Compiler bug JAGae51337
surfsub.o:	surfsub.m
	cat ../utilities/gener.m  surfsub.m | $(M4) $(M4OPTS) > surfsub.f
	$(FC)	$(FFLAGSS) surfsub.f

# Compiler bug JAGae55370
mopac.o:	mopac.m
	cat ../utilities/gener.m  mopac.m | $(M4) $(M4OPTS) > mopac.f
	$(FC)	$(FFLAGSS) mopac.f

# Compiler bug JAGae56538
nbo.o:	nbo.m
	cat ../utilities/gener.m  nbo.m | $(M4) $(M4OPTS) > nbo.f
	$(FC)	$(FFLAGSS) nbo.f

#### Information on the machine this build was configured on ####
# Output of uname -a:
# HP-UX ops59 B.11.23 U ia64 4200161395 unlimited-user license
#
# Output of model:
# ia64 hp superdome server SD64A
#
# Output of f90 +version:
# HP F90 v2.9.2
