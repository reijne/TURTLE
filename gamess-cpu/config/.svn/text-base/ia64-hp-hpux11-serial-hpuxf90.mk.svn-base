#doc  Makefile settings for the GAMESS-UK serial build on HP hardware running HPUX
#doc  on Itanium 2 processors and using the HP F90 v2.8.7 compiler.
#doc  The veclib ibrary is selected by default. This supports i8.
#doc
#doc  Options
#doc  i8 - switch from the default integer*4 build to integer*8
#
#dopt drf mrdci zora nbo vb mopac vdw sysmo dl-find veclib i8
#opt 
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


#--#if i8#
FFLAGI8= +i8
CFLAGI8= -DLONG_INTEGER
#--#else#
FFLAGI8=
CFLAGI8=
#--#endif i8#

# Using HPs veclib numerical library
BLASOPT=,blas
LIBBLAS = /opt/mlib/lib/hpux64/libveclib8.a

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
# Compiler bug JAGae55172
casb.o:	casb.m
	cat ../utilities/gener.m  casb.m | $(M4) $(M4OPTS) > casb.f
	$(FC)   $(FFLAGSS) casb.f

mrdci1.o:      mrdci1.m
	cat ../utilities/gener.m  mrdci1.m | $(M4) $(M4OPTS) > mrdci1.f
	$(FC)   $(FFLAGS2) mrdci1.f

newmrd1.o:      newmrd1.m
	cat ../utilities/gener.m  newmrd1.m | $(M4) $(M4OPTS) > newmrd1.f
	$(FC)   $(FFLAGSV) newmrd1.f

newmrd2.o:      newmrd2.m
	cat ../utilities/gener.m  newmrd2.m | $(M4) $(M4OPTS) > newmrd2.f
	$(FC)   $(FFLAGSV) newmrd2.f

newmrd3.o:      newmrd3.m
	cat ../utilities/gener.m  newmrd3.m | $(M4) $(M4OPTS) > newmrd3.f
	$(FC)   $(FFLAGSV) newmrd3.f

newmrd4.o:      newmrd4.m
	cat ../utilities/gener.m  newmrd4.m | $(M4) $(M4OPTS) > newmrd4.f
		$(FC)   $(FFLAGSV) newmrd4.f

newmrd5.o:      newmrd5.m
	cat ../utilities/gener.m  newmrd5.m | $(M4) $(M4OPTS) > newmrd5.f
	$(FC)   $(FFLAGSV) newmrd5.f

newmrd6.o:      newmrd6.m
	cat ../utilities/gener.m  newmrd6.m | $(M4) $(M4OPTS) > newmrd6.f
	$(FC)   $(FFLAGSV) newmrd6.f

# Compiler bug JAGae51271
scf.o:  scf.m
	cat ../utilities/gener.m  scf.m | $(M4) $(M4OPTS) > scf.f
	$(FC)   $(FFLAGSS) scf.f

# Compiler bug JAGae51337
surfsub.o:      surfsub.m
	cat ../utilities/gener.m surfsub.m | $(M4) $(M4OPTS) > surfsub.f
	$(FC)   $(FFLAGSS) surfsub.f

# Compiler bug JAGae55370
mopac.o:        mopac.m
	cat ../utilities/gener.m  mopac.m | $(M4) $(M4OPTS) > mopac.f
	$(FC)   $(FFLAGSS) mopac.f

# Compiler bug JAGae56538
nbo.o:	nbo.m
	cat ../utilities/gener.m  nbo.m | $(M4) $(M4OPTS) > nbo.f
	$(FC)   $(FFLAGSS) nbo.f
