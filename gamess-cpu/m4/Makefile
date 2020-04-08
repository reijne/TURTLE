#
#  $Author: jmht $
#  $Date: 2008-01-21 15:20:37 $
#  $Revision: 1.316 $
#  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/Makefile.in,v $
#  $State: Exp $
#
#   makefile for GAMESS Version 8.0 direct from M4 
#
#startdoc
#
#   issue ./configure and then follow the instructions. This will generate
#   a file configure.dat containing a set of m4 keywords which control
#   pre-processing of the M4 source, listed below.
#   Please be aware that there is *not* a direct correspondence between 
#   configure keywords and the resulting M4 keywords.
#
# restricted builds: base         - mimimal SCF/DFT 
#                    mp2          - SCF, DFT, MP2, SCF Hessian (not qm/mm)
#                    ci           - SCF, DFT, MP2, SCF Hessian (not qm/mm), direct-CI, GF
#
# compile keywords:  debug
#                    profile
#                    protect-source   chmod o-w on fortran/C files
#                    timings  Detailed timing breakdown
#                    f77      Build without f90 stuff
#                    f90test  Use F90 Test features (developmental)
#                    i8       Request 64 bit integer variables
#                    blas     use blas libraries
#                    veclib   use veclib libraries on HP (advised)
#
# optional modules:  diesel   include diesel CI code 
#                             [requires code in ../diesel]
#                    drf      include DRF code
#                             [requires code in ../drf]
#                    mopac
#                 or mopac7   include mopac version 7 code
#                             [requires code in ../mopac]
#                    mopac6   include mopac version 6 code
#                    zora     include developmental ZORA code
#                    newscf   include developmental SCF code (F90 version)
#                             [requires code in ../newscf]
#                    dl-find  include DL-FIND interface
#                    vb       include developmental VB code
#                    vdw      include developmental Van der Waals code
#                    sysmo    include developmental sysmo code
#                    masscf   include developmental parallel mcscf code
#
# parallel keywords: tcgmsg   parallel code (no GA)
#                             (available for base build only)
#                    mpi      parallel code (no GA) using MPI
#                             (available for base (+VB) build only)
#                    ga       link GA tools
#                             (available for base + mp2 + ci builds only)
#                    peigs    PeIGS Parallel diagonaliser 
#                    scalapack ScaLAPACK Parallel diagonalisers
#
# keywords determining the system size limits:
#
#                    demo     200  basis functions static (macosx)
#                    small    1024 basis functions
#                    medium   2048 basis functions
#                    large    4096 basis functions
#                    huge     8192 basis functions   2000 atoms
#                    giga     16384 basis functions  4000 atoms
#                    tera     19900 basis functions  5000 atoms
#                    qmmm     3000 nuclear centres 
#                             also suppresses some large (eg NAT*NAT)
#                             memory arrays
#
# Additional options
#
#                    nbo      Activate NBO code
#                    mrdci    Activate new MRDCI code
#                    old-dft  Activate old post-HF DFT code
#                    rpagrad  RPA excited state gradients
#
#enddoc
#
#   WARNING - make clean also deletes fortran and C source 
#             leaving only .m files
#
#  Most processing is controlled by the machine file in directory 
#  ../machines, supplemented by extra M4 control flags  (see PAR_M4_OPTS
#  PEIGS_M4_OPTS, DFT_M4_OPTS, DRF_M4_OPTS, MOPAC_M4_OPTS, NEWSCF_M4_OPTS
#  settings below) which  are passed to m4 on the command line
#
# ----------------------- CONVENTIONS ------------------------------
#
# some of the current variable usage:
#
#  GA_F77_DEFS   cpp defines and flags for GA codes
#  GA_F77_INC    cpp include directory for GA codes (F77)
#  GA_C_INC      cpp includes for GA codes (C)
#  PAR_LIBS      non-mpi libraries (GA, TCGMSG etc) (-L and -l flags)
#  LIBMPI        -l settings for MPI libraries
#  MPI_LIBS      holds -L and -l flags for MPI libraries
#  MPI_LIB       directory for MPI libraries
#
# ----------------------- PATHNAMES ------------------------------
#
# target names for executables
#
EXE= ../bin/gamess
#
# ----------------------- DEFAULT MAKE TARGET ---------------------
#
GAMESS_TARGET=exe
#
# ----------------------- DEFAULT DEFINITIONS ---------------------
#
.SUFFIXES:
.SUFFIXES:	.m .o .s .c .f .f90
SHELL=/bin/sh
M4 = m4
#
FC = f77
FC90 = f90
CC = cc
CC = cc
CXX = cc
AS = as
LD = f77
#Need seperate LDNAME as windows doesn't take -o flag
LDNAME = -o 
CPP = cpp
FFLAGSV = -c
FFLAGSS = -c
FFLAGSN = -c
FFLAGSV90 = $(FFLAGSV)
FFLAGSS90 = $(FFLAGSS)
FFLAGSN90 = $(FFLAGSN)
FFLAGS_CPP =
CFLAGS = -c
CXXFLAGS = -c
GMAKE  = gmake
ARCHIVE  =  ar rcv
ARCHIVE_NOREPLACE  =  ar qcv
UNARCHIVE = ar xo
RANLIB   =  ar ts
# should be set to the same flags as used to build fortran code

# whoami used in vb Makefile, but not found in path on some machines (e.g. sparc at lomond)
WHOAMI=whoami

# M4 Stuff
INCLUDE=_INCLUDE($$1)
MACHINE_KEY=G
IAND32=and($$1,$$2)
IOR32=or($$1,$$2)
IXOR32=xor($$1,$$2)
IAND64=iand($$1,$$2)
IOR64=or($$1,$$2)
IXOR64=xor($$1,$$2)
SHIFT=shift($$1,$$2)
COMPLEX=complex*16
REAL=real*8
MPIINT=integer*4
MPILOG=logical*4
DFLOAT=dble($$1)
EXTRA_M4_DEFINITIONS=
# --------------------------- LINEAR ALGEBRA LIBRARY --------------------------
# required for all builds
MAKE_LINALG=check_linalg
LINALG_LIB =   ../linalg/liblinalg.a
# ---------------------------- GLOBAL ARRAY OPTIONS ---------------------------
# set to make_ga if GA tools need to be built
GA_INTEGER_DECL=integer
GA_LOGICAL_DECL=logical
#
# -------------------------- PEIGS definitions  ------------------
#
#
#
# ---------------------------- F90 OPTIONS ---------------------------
#
# Use the f90test keyword to switch on the use of F90 test code in analc
#         f77 to switch off F90 production code
#
#
# ----------------------- CCPDFT OPTIONS ------------------------------------
MAKE_DFT=check_dft
DFT_M4_OPTS=,ccpdft
DFT_LIB =   ../dft/libdft.a
#
# ----------------------- DIESEL OPTIONS -----------------------------------
# ----------------------- DRF OPTIONS --------------------------------------
HONDO_OBJ=gamess_hondo.o
MAKE_DRF=check_drf
DRF_M4_OPTS=,drf
DRF_LIB =   ../drf/libdrf.a
# ----------------------- MOPAC and MOPAC7 OPTIONS --------------------------
#
#------------------------ ZORA OPTIONS --------------------------------------
MAKE_ZORA=check_zora
ZORA_M4_OPTS=,zora
ZORA_LIB = ./zora.o
#
#------------------------ NMR OPTIONS --------------------------------------
#
#------------------------ VB OPTIONS ----------------------------------------
MAKE_VB=check_vb
VB_M4_OPTS=,vb
VB_LIB = ../vb/libvb.a 
# Specify a default for MPI_INCLUDE as it is needed in vb/Makefile for sendrecv
MPI_INCLUDE=../g/include
VB_EXTRA=
# ----------------------- Van der Waals Options -----------------------------
VDW_M4_OPTS=,vdw
VDW_EXTRA=vdwaals_corr.o inter_vdwaals.o
# ----------------------- Sysmo Options -----------------------------
SYSMO_M4_OPTS=,sysmo
SYSMO_EXTRA=sysmo.o
# ----------------------- Parallel MSCSF (MASSCF) Options -----------------------------
MASSCF_M4_OPTS=
MASSCF_EXTRA=
# ----------------------- NEWSCF OPTIONS ----------------------
#
#
#--------------------------- XML Options ----------------------
#
#
# ----------------------- DL_FIND OPTIONS ----------------------
#
MAKE_DLFIND=check_dlf
DLFIND_M4_OPTS=,dl-find
DLFIND_LIB =   ../dl-find/libdlf.a
#
#
#
# ----------------------- TASKFARM OPTIONS ------------------------------------
#
#
#
#
#--------------------------- PARALLEL M4 Options -----------------------------
#
# These options define the M4 options that determine which parallel 
# functionality is included in the various builds.
#
#
#
# Default case
PAR_M4_OPTS=serial
#
#
# ----------------------- SYSTEM SIZE OPTIONS --------------------------------
# defaults ( currently the same as small )
#
#
MAXORB	=	1024
MXSHEL	=	600
MXPRIM	=	2000
MXGRPS	=	560
MXGAUS	=	5600
MAXAT   =       750
MAXNZ   =       700
MAXVAR  =       2000
# VB turtleparam
MXORBVB =       255
MAXEXVB =       20000

# modifications requested by configure options



MAXORB	=	4096
MXORBVB =    1000
MXSHEL	=	2048
MXPRIM	=	8192
MXGAUS	=	11600
#
#
# ----------------------- QM/MM INTERFACE OPTIONS --------------------
#
# Add qmmm here to restrict memory use for large number of qm centres
#
#
# ------------------ DEVELOPMENT OPTIONS --------------------------
D1=,nbo
D2=,zora
D4=,mrdci
D6=,vb
D11=,vdw
D12=,sysmo
DEV_M4_OPTS=$(D1)$(D2)$(D3)$(D4)$(D5)$(D6)$(D7)$(D8)$(D9)$(D10)$(D11)$(D12)
#
#--------------------------- Debugging & Profiling -------------------
#
#
#
#--------------------------- Target ----------------------------------
#
all:	errchk utilities $(GAMESS_TARGET)

utilities:
	(cd ../utilities; make)
#
# Dummy target to force an object to be rebuilt each time
FORCE:
#
# ======================= Machine options =============================
#
#doc  MK file for the serial build of GAMESS-UK on an Opteron processor with
#doc  PGI compilers (tested with 5.2, 6.0 and 6.1)
#doc 
#doc  Options:
#doc  blas - build against blas (expected in: /usr/pgi/linux86-64/6.1/lib/libacml.a)
#doc  i8   - build with 8-byte integers.
##
# Default options:
#dopt vb nbo mrdci zora drf vdw sysmo dl-find
#opt blas static i8 mopac
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=opteron,linux,pclinux,littleendian,cio,unix,doublebackslash,upck-equiv

# Bitwise operators that differ from the defaults.
IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)
# 
#
RANLIB=ar -s
CPP=/lib/cpp
GA_F77_DEFS = -traditional
#
# ===============  Compiler Options
#
#ensure 64bit compilers invoked
FC = pgf90
FC90 = pgf90
LD = pgf90
CC = pgcc

FFLAGSI8=-i8
CFLAGSI8=-DLONG_INTEGER
M4_I8_OPT=,i8

# Default compilation flags
FFLAGSD = -c ${FFLAGSI8}
CFLAGSD= -c ${CFLAGSI8}

FFLAGSV =  ${FFLAGSD} -tp=k8-64 -Mcache_align -O2
FFLAGSS =  ${FFLAGSD} -tp=k8-64 -Mcache_align -O
FFLAGSN =  ${FFLAGSD} -tp=k8-64 -Mcache_align -O1
FFLAGSN0 = ${FFLAGSD}
CFLAGS  = ${CFLAGSD} -tp=k8-64 -D_REENTRANT
LDTMP = -tp=k8-64
#
LDFLAGS=${LDTMP}
#
#
#
OPTIONS=${MACHOPT}${BLASOPT}${M4_I8_OPT}
BL_LIB=${LIBBLAS}
#
# ===============  Additional Files 
#
EXTRA=gethes.o sp0011.o chfeq.o bmove.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o
EXTRA_DFT = jkint_dft.o jkder_dft.o
#
# ===============  Compiler Exceptions 
#
rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

eispack.o:	eispack.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

gethes.o:	casb.m
	cat  ../utilities/gener.m casb.m | $(M4) -DGEN_EXTRACTFILE=gethes $(M4OPTS) > gethes.f
	$(FC) $(FFLAGSN) gethes.f

sp0011.o:	integs.m
	cat  ../utilities/gener.m integs.m | $(M4) -DGEN_EXTRACTFILE=sp0011 $(M4OPTS) > sp0011.f
	$(FC) $(FFLAGSN) sp0011.f

chfeq.o:	cphf.m
	cat  ../utilities/gener.m cphf.m | $(M4) -DGEN_EXTRACTFILE=chfeq $(M4OPTS) > chfeq.f
	$(FC) $(FFLAGSN) chfeq.f

bmove.o:	mclr.m
	cat  ../utilities/gener.m mclr.m | $(M4) -DGEN_EXTRACTFILE=bmove $(M4OPTS) > bmove.f
	$(FC) $(FFLAGSN) bmove.f

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
#
#  ========== DFT Exceptions (PGF) ===============
#
jkint_dft.o:	integ2e.m
	cat ../machines/$(MACH) ../utilities/gener.m  integ2e.m | $(M4) -DGEN_EXTRACTFILE=jkint_dft $(M4OPTS) $(SNGL) > jkint_dft.f
	$(FC) $(FFLAGSN) jkint_dft.f
#
jkder_dft.o:	deriv2e.m
	cat ../machines/$(MACH) ../utilities/gener.m  deriv2e.m | $(M4) -DGEN_EXTRACTFILE=jkder_dft $(M4OPTS) $(SNGL) > jkder_dft.f
	$(FC) $(FFLAGSN) jkder_dft.f

## Build Documentation
#
# Output of uname -a
# Linux scarf 2.4.21-27.0.2.ELsmp #1 SMP Wed Jan 12 23:25:44 EST 2005 x86_64 x86_64 x86_64 GNU/Linux
# 
# Output of cat /proc/cpuinfo
# processor	: 0
# vendor_id	: AuthenticAMD
# cpu family	: 15
# model		: 5
# model name	: AMD Opteron(tm) Processor 248
# physical id	: 0
# siblings	: 1
# stepping	: 10
# cpu MHz		: 2193.497
# cache size	: 1024 KB
# fpu		: yes
# fpu_exception	: yes
# cpuid level	: 1
# wp		: yes
# flags		: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 syscall nx mmxext lm 3dnowext 3dnow
# bogomips	: 4377.80
# TLB size	: 1088 4K pages
# clflush size	: 64
# address sizes	: 40 bits physical, 48 bits virtual
# power management: ts fid vid ttp
# 
# Output of cat /etc/issue
# Red Hat Enterprise Linux AS release 3 (Taroon Update 4)
# 
# Output of /usr/pgi/linux86-64/5.2/bin/pgf77 -V
# pgf77 5.2-4
# Copyright 1989-2000, The Portland Group, Inc.  All Rights Reserved.
# Copyright 2000-2004, STMicroelectronics, Inc.  All Rights Reserved.
# 
# Output of /usr/pgi/linux86-64/5.2/bin/pgfcc -V
# pgcc 5.2-4
# Copyright 1989-2000, The Portland Group, Inc.  All Rights Reserved.
# Copyright 2000-2004, STMicroelectronics, Inc.  All Rights Reserved.
#
#
# ======================= End of machine options =============================
#

check_linalg:
	(cd ../linalg; make linalg)
$(LINALG_LIB):
	make check_linalg

check_dft:
	(cd ../dft; make dftlib )
$(DFT_LIB):
	make check_dft


check_drf:
	(cd ../drf; make drflib )
$(DRF_LIB):
	make check_drf


check_zora:
	(cd ../m4; make zora.o )


check_vb:
	(cd ../vb; make vb )
$(VB_LIB):
	make check_vb



check_dlf:
	(cd ../dl-find; $(GMAKE) libdlf.a )
$(DLFINDLIB):
	make check_dlf

#
MODOBJ = $(MODF90) $(MODF90TEST) 
#
# for qmmm work, additionally remove optef 
# delete common/small from mains (scratch in util1 - but check)
#
OPTEF = optef.o

#
# combine library requirements
#
LDLIBS	= 	$(DFT_LIB) $(DRF_LIB) $(MOPAC_LIB) $(ZORA_LIB) $(VB_LIB) $(NEWSCF_LIB) $(NEWSCF_LIB_F90) $(DLFIND_LIB) $(NMR_LIB) $(XML_LIB) $(AGENTX_LIB) $(GA_LIBS) $(PEIGS_LIB) $(MA_LIB) $(BL_LIB) $(LINALG_LIB)

MISC	=	$(PEIGS_OBJ) $(TASKFARM_OBJ) $(GROMACS_OBJ)
#

#
# ------------- Processing of build selection --------------
#
#
# This target is overridden by a fatal one if the build
# cannot proceed because of incompatible choices
#
errchk:
	@echo Building full code

OBJ_VERSION = gms_version.o

OBJ_BASE = mains.o anala.o analb.o analc.o analf.o analg.o basis.o basis1.o basis2.o basis3.o \
	dft.o input.o intega.o \
	machscf.o  drv2e.o drv80.o integs.o scf.o server.o optim.o ${OPTEF} \
	integb.o integb_nl.o integb_lib.o integc.o guess.o master.o  util1.o util2.o util3.o \
	parallel.o  mass.o drv1e.o c.o  $(MODOBJ) \
	morokuma.o nbo.o $(MISC) omit.o ${EXTRA_BASE} $(VB_EXTRA) $(VDW_EXTRA) $(SYSMO_EXTRA) $(HONDO_OBJ) $(EXTRA) \
	gdiesel.o glue_uk_us.o $(OBJ_VERSION)

# Additional files for GA case
OBJ_GA_SCF = ga.o

# These files appear in all MP2+second derivative builds
OBJ_MP2_COMMON = \
	mainci.o anald.o anale.o \
	cphf.o derdrv.o \
	drvmp.o inpci.o  \
	integd.o intege.o machci.o \
	tran4.o index4.o \
	util4.o util5.o util6.o util7.o util8.o \
	sec1e.o sec2e.o ${EXTRA_MP2}

# files to get MP2 and SCF second derivatives (serial case)
OBJ_MP2 = $(OBJ_MP2_COMMON) mp2a.o mp2b.o mp2c.o mp2d.o mp3.o \
	secmp2.o

OBJ_GA_MP2 = $(OBJ_MP2_COMMON) mp2_parallel.o \
	secd_parallel.o secdchf.o dirrpa.o rpa.o direct.o 

# files for Direct-CI, GFF (additional parallel CI modules here)
OBJ_CI = dircta.o dirctb.o dirctc.o dirctd.o gff.o ${EXTRA_CI}

# Remaining object files (post-HF etc)
# Supported for serial code only
OBJ_SERIAL =   tdaf.o model.o \
	newmrd1.o \
	newmrd2.o newmrd3.o newmrd4.o newmrd5.o newmrd6.o \
	mcscfa.o mcscfb.o mcscfc.o \
        mrdci1.o mrdci2.o mrdci3.o mrdci4.o mrdci5.o mrdci6.o mrdci7.o fullci.o \
        direct.o casa.o casb.o \
        ccsd.o dirrpa.o mclr.o nvccsd.o rpa.o rpagrad.o tsort.o \
	tsortc.o 

# Construct list of object files

#
#   Complete code (only valid for the serial case
#
TOBJ = ${OBJ_BASE} ${OBJ_MP2} ${OBJ_SERIAL} ${OBJ_CI}
OMIT = ,full_build


exe:	$(EXE)

$(EXE):	 $(MAKE_LINALG) $(MAKE_GA) $(MAKE_PEIGS) $(MAKE_MA) $(MAKE_DFT) $(MAKE_MOPAC) $(MAKE_DIESEL) $(MAKE_DRF) $(MAKE_NEWSCF) $(MAKE_DLFIND) $(MAKE_VB) $(MAKE_ZORA) $(MAKE_NMR) $(MAKE_XML) $(MAKE_DIESEL) $(TOBJ) 
	$(LD) $(LDNAME)$(EXE) $(LDFLAGS) $(TOBJ) $(ASOBJ) $(LDLIBS) 
	@echo ----- made $(EXE) -----




# Combine M4 processing options together

M4STUFF= 

M4OPTS	= -DGEN_OPTIONS=${OPTIONS} -DGEN_MACHINE=${MACHINE_KEY} \
		-DREAL='${REAL}' -DCOMPLEX='${COMPLEX}' -DINCLUDE='${INCLUDE}' \
		-Ddfloat='${DFLOAT}' -DMPIINT='${MPIINT}' -DMPILOG='${MPILOG}' \
		-DIAND32='${IAND32}' -DIOR32='${IOR32}' -DIXOR32='${IXOR32}' \
		-DIAND64='${IAND64}' -DIOR64='${IOR64}' -DIXOR64='${IXOR64}' \
		-DSHIFT='${SHIFT}' -DQ='qq(ivoff+$$1)' ${EXTRA_M4_DEFINITIONS} \
		-DM4_OPTIONS=$(PAR_M4_OPTS)$(MACH_M4_OPTS)$(PEIGS_M4_OPTS)$(PDIAG_M4_OPTS)$(DFT_M4_OPTS)$(DIESEL_M4_OPTS)$(DRF_M4_OPTS)$(MOPAC_M4_OPTS)$(VDW_M4_OPTS)$(SYSMO_M4_OPTS)$(MASSCF_M4_OPTS)$(DLFIND_M4_OPTS)$(XML_M4_OPTS)$(QMMM_M4_OPTS)$(NEWSCF_M4_OPTS)$(TASKFARM_M4_OPTS)$(DEV_M4_OPTS)$(F90_M4_OPTS)${OMIT}$(DATAIN_M4_OPTS)$(PROF_M4_OPTS)$(DEBUG_M4_OPTS)$(PLATFORM_M4_OPTS)${MA_M4_OPTS}${SCALAPACK_M4_OPTS} \
		-DMA_LOGICAL=${GA_LOGICAL_DECL} \
		-DMA_INTEGER=${GA_INTEGER_DECL} \
		-DMAXAT=$(MAXAT)   -DMAXORB=$(MAXORB) \
		-DMAXVAR=$(MAXVAR) -DMAXNZ=$(MAXNZ)   \
		-DMXSHEL=$(MXSHEL) -DMXPRIM=$(MXPRIM) \
		-DMXGRPS=$(MXGRPS) -DMXGAUS=$(MXGAUS) \
		-DMXORBVB=$(MXORBVB) -DMAXEXVB=$(MAXEXVB) -
#
# ======================== Rules to build object code =========================
#
.m.o:
	rm -f $*.f
	cat ../utilities/gener.m $*.m | $(M4) $(M4OPTS)  >> $*.f
	$(FC) $(FFLAGSV) $(OBJNAME) $*.f
#
#
#  default rules
#
.f.o:
	$(FC) $(FFLAGSV) $(OBJNAME) $*.f
.f90.o:
	$(FC90) $(FFLAGSV90) $(OBJNAME) $*.f90
#
#  explicit rules so that .m->.o conversion goes by .c (not .f)
#
c.o:	c.c
	$(CC) $(CFLAGS) $(OBJNAME) c.c
#
cmpi.o:	cmpi.c
	$(CC) $(CFLAGS) $(OBJNAME) cmpi.c
#
tsortc.o:	tsortc.c
	$(CC) $(CFLAGS) $(OBJNAME) tsortc.c
#

#
#  parallel.m is the first source file that needs CPP processing
#  for GA tools etc
#
parallel.o:	parallel.m
	rm -f parallel.GUK parallel.f
	cat  ../utilities/gener.m parallel.m | $(M4) $(M4OPTS)  >> parallel.f
	$(FC) $(FFLAGSS) $(OBJNAME) parallel.f

#  glue_uk_us.m is the second...
glue_uk_us.o:	glue_uk_us.m
	rm -f glue_uk_us.GUK glue_uk_us.f
	$(CPP) -traditional $(GA_F77_INC) $(GA_F77_DEFS)  < glue_uk_us.m | sed '/^#/D' >> glue_uk_us.GUK
	cat  ../utilities/gener.m glue_uk_us.GUK | $(M4) $(M4OPTS)  > glue_uk_us.f
	rm -f glue_uk_us.GUK
	$(FC) $(FFLAGSS) $(OBJNAME) glue_uk_us.f

#
# ======================== Rules to build source code =========================
#
.m.f:	
	rm -f $*.f
	cat ../utilities/gener.m $*.m | $(M4) $(M4OPTS)  >> $*.f
#
# special rule to include CPP
#
parallel.f:	parallel.m
	rm -f parallel.GUK parallel.f
	cat ../utilities/gener.m parallel.m | $(M4) $(M4OPTS)  >> parallel.f
#
.m.c:
	rm -f $*.c
	cat ../utilities/gener.m $*.m | $(M4) $(M4OPTS) >> $*.c
.m.s:
	cat ../utilities/gener.m $*.m | $(M4) $(M4OPTS)> $*.s
#
#
# ============================= Common Block Dependencies  ===========================
#
# After adding new code, this list can be regenerated using the 
# makedeps script in ../utilities
#
anala.o:		common/atmol3
anala.o:		common/canon_nat
anala.o:		common/ccpdft.hf77
anala.o:		common/chmgms
anala.o:		common/cigrad
anala.o:		common/cndx41
anala.o:		common/common
anala.o:		common/cslosc
anala.o:		common/discc
anala.o:		common/dm
anala.o:		common/dump3
anala.o:		common/gaanal
anala.o:		common/gmempara
anala.o:		common/groups
anala.o:		common/harmon
anala.o:		common/infoa
anala.o:		common/iofile
anala.o:		common/limy
anala.o:		common/machin
anala.o:		common/mapper
anala.o:		common/newscf
anala.o:		common/nodeio
anala.o:		common/nshel
anala.o:		common/parcntl
anala.o:		common/prints
anala.o:		common/restar
anala.o:		common/restri
anala.o:		common/restrj
anala.o:		common/runlab
anala.o:		common/scfopt
anala.o:		common/scfwfn
anala.o:		common/scra7
anala.o:		common/sector
anala.o:		common/segm
anala.o:		common/sizes
anala.o:		common/statis
anala.o:		common/symtry
anala.o:		common/timez
anala.o:		common/tran
anala.o:		common/unocas
anala.o:		common/xfield
anala.o:		common/zorac
analb.o:		common/atmol3
analb.o:		common/datplt
analb.o:		common/datprp
analb.o:		common/discc
analb.o:		common/dump3
analb.o:		common/infoa
analb.o:		common/iofile
analb.o:		common/limy
analb.o:		common/linkan
analb.o:		common/machin
analb.o:		common/mapper
analb.o:		common/nshel
analb.o:		common/overrule
analb.o:		common/prints
analb.o:		common/restar
analb.o:		common/restrj
analb.o:		common/root
analb.o:		common/runlab
analb.o:		common/scra7
analb.o:		common/sector
analb.o:		common/sizes
analb.o:		common/timez
analb.o:		common/tran
analb.o:		common/vdwrad
analc.o:		common/atmol3
analc.o:		common/coval
analc.o:		common/direc
analc.o:		common/discc
analc.o:		common/hermit
analc.o:		common/infoa
analc.o:		common/infob
analc.o:		common/iofile
analc.o:		common/limy
analc.o:		common/machin
analc.o:		common/mapper
analc.o:		common/nshel
analc.o:		common/restar
analc.o:		common/restrj
analc.o:		common/runlab
analc.o:		common/scfwfn
analc.o:		common/scra7
analc.o:		common/segm
analc.o:		common/sizes
analc.o:		common/symtry
analc.o:		common/timez
analc.o:		common/tran
analc.o:		common/transf
analc.o:		common/unocas
analc.o:		common/wermit
analc.o:		common/work
analc.o:		junk_dma.o
anald.o:		common/atmblk
anald.o:		common/cndx40
anald.o:		common/cndx41
anald.o:		common/common
anald.o:		common/cslosc
anald.o:		common/drive_dft
anald.o:		common/ijlab
anald.o:		common/infoa
anald.o:		common/iofile
anald.o:		common/machin
anald.o:		common/mapper
anald.o:		common/nshel
anald.o:		common/prnprn
anald.o:		common/restar
anald.o:		common/restri
anald.o:		common/restrj
anald.o:		common/restrl
anald.o:		common/runlab
anald.o:		common/segm
anald.o:		common/sizes
anald.o:		common/statis
anald.o:		common/symtry
anald.o:		common/tdhf
anald.o:		common/tdhfx
anald.o:		common/timez
anald.o:		common/vcore
anald.o:		common/vectrn
anale.o:		common/common
anale.o:		common/funcon
anale.o:		common/gmempara
anale.o:		common/infoa
anale.o:		common/iofile
anale.o:		common/mapper
anale.o:		common/prnprn
anale.o:		common/runlab
anale.o:		common/runopt
anale.o:		common/sizes
analf.o:		common/atmol3
analf.o:		common/cslosc
analf.o:		common/direc
analf.o:		common/dump3
analf.o:		common/g80nb
analf.o:		common/gmempara
analf.o:		common/infoa
analf.o:		common/iofile
analf.o:		common/limy
analf.o:		common/machin
analf.o:		common/mapper
analf.o:		common/modj
analf.o:		common/nshel
analf.o:		common/parallel
analf.o:		common/phycon
analf.o:		common/prints
analf.o:		common/psscrf
analf.o:		common/restar
analf.o:		common/root
analf.o:		common/runlab
analf.o:		common/scra7
analf.o:		common/segm
analf.o:		common/sizes
analf.o:		common/statis
analf.o:		common/symtry
analf.o:		common/timez
analf.o:		common/tran
analf.o:		common/vdwrad
analg.o:		common/atmol3
analg.o:		common/cndx41
analg.o:		common/datprp
analg.o:		common/discc
analg.o:		common/dm
analg.o:		common/dump3
analg.o:		common/gvalue
analg.o:		common/hermit
analg.o:		common/infoa
analg.o:		common/iofile
analg.o:		common/limy
analg.o:		common/mapper
analg.o:		common/natorb
analg.o:		common/nshel
analg.o:		common/restar
analg.o:		common/restrj
analg.o:		common/runlab
analg.o:		common/scfwfn
analg.o:		common/scra7
analg.o:		common/segm
analg.o:		common/sizes
analg.o:		common/timez
analg.o:		common/tran
analg.o:		common/unocas
analg.o:		common/wermit
basis.o:		common/atmol3
basis.o:		common/chmgms
basis.o:		common/datbas
basis.o:		common/fsymas
basis.o:		common/ghost
basis.o:		common/harmon
basis.o:		common/infoa
basis.o:		common/infoa
basis.o:		common/infob
basis.o:		common/iofile
basis.o:		common/iofile
basis.o:		common/machin
basis.o:		common/nshel
basis.o:		common/nshel
basis.o:		common/prints
basis.o:		common/restar
basis.o:		common/restrj
basis.o:		common/runlab
basis.o:		common/runlab
basis.o:		common/scra7
basis.o:		common/sizes
basis.o:		common/symtry
basis.o:		common/tol
basis.o:		common/transf
basis.o:		common/work
basis.o:		common/xmlin
basis.o:		common/zorac
basis1.o:		common/restar
basis1.o:		common/sizes
basis2.o:		common/iofile
basis2.o:		common/restar
basis2.o:		common/sizes
basis3.o:		common/restar
basis3.o:		common/sizes
casa.o:		common/atmblk
casa.o:		common/avstat
casa.o:		common/blksiz
casa.o:		common/caspar
casa.o:		common/ciconv
casa.o:		common/cndx40
casa.o:		common/ctrl
casa.o:		common/discc
casa.o:		common/dm
casa.o:		common/dump3
casa.o:		common/exc
casa.o:		common/filel
casa.o:		common/files
casa.o:		common/finish
casa.o:		common/gjs
casa.o:		common/harmon
casa.o:		common/iofile
casa.o:		common/machin
casa.o:		common/mapper
casa.o:		common/popnos
casa.o:		common/prints
casa.o:		common/qice
casa.o:		common/restar
casa.o:		common/runlab
casa.o:		common/segm
casa.o:		common/simul
casa.o:		common/sizes
casa.o:		common/stak
casa.o:		common/trntim
casb.o:		common/atmblk
casb.o:		common/blksiz
casb.o:		common/ciconv
casb.o:		common/cndx40
casb.o:		common/ctrl
casb.o:		common/dm
casb.o:		common/dump3
casb.o:		common/exc
casb.o:		common/filel
casb.o:		common/files
casb.o:		common/finish
casb.o:		common/gjs
casb.o:		common/harmon
casb.o:		common/infoa
casb.o:		common/iofile
casb.o:		common/machin
casb.o:		common/mapper
casb.o:		common/nshel
casb.o:		common/popnos
casb.o:		common/prints
casb.o:		common/prnprn
casb.o:		common/qice
casb.o:		common/restar
casb.o:		common/runlab
casb.o:		common/scra7
casb.o:		common/simul
casb.o:		common/sizes
casb.o:		common/stak
casb.o:		common/symtry
casb.o:		common/timez
ccsd.o:		common/atmblk
ccsd.o:		common/direc
ccsd.o:		common/discc
ccsd.o:		common/filel
ccsd.o:		common/iofile
ccsd.o:		common/machin
ccsd.o:		common/restar
ccsd.o:		common/runlab
ccsd.o:		common/savem
ccsd.o:		common/segm
ccsd.o:		common/sizes
ccsd.o:		common/statis
ccsd.o:		common/work
charmm_omit.o:		common/blur
charmm_omit.o:		common/chmgms
charmm_omit.o:		common/funct
charmm_omit.o:		common/phycon
charmm_omit.o:		common/sizes
cphf.o:		common/atmblk
cphf.o:		common/ccpdft.hf77
cphf.o:		common/cigrad
cphf.o:		common/cndx40
cphf.o:		common/cndx41
cphf.o:		common/common
cphf.o:		common/drive_dft
cphf.o:		common/ghfblk
cphf.o:		common/gmempara
cphf.o:		common/infoa
cphf.o:		common/iofile
cphf.o:		common/mapper
cphf.o:		common/mp2grd
cphf.o:		common/nshel
cphf.o:		common/prnprn
cphf.o:		common/rpadcom
cphf.o:		common/sizes
cphf.o:		common/statis
cphf.o:		common/symtry
cphf.o:		common/tdhf
cphf.o:		common/tdhfx
cphf.o:		common/timeperiods
cphf.o:		common/timez
cphf.o:		common/vcore
derdrv.o:		common/atmblk
derdrv.o:		common/atmol3
derdrv.o:		common/ccpdft.hf77
derdrv.o:		common/cigrad
derdrv.o:		common/cndx40
derdrv.o:		common/cndx41
derdrv.o:		common/common
derdrv.o:		common/cslosc
derdrv.o:		common/debug
derdrv.o:		common/disp
derdrv.o:		common/drive_dft
derdrv.o:		common/dump3
derdrv.o:		common/funct
derdrv.o:		common/ghfblk
derdrv.o:		common/gmempara
derdrv.o:		common/grad2
derdrv.o:		common/harmon
derdrv.o:		common/ijlab
derdrv.o:		common/infoa
derdrv.o:		common/iofile
derdrv.o:		common/machin
derdrv.o:		common/mapper
derdrv.o:		common/nshel
derdrv.o:		common/prnprn
derdrv.o:		common/restar
derdrv.o:		common/restrj
derdrv.o:		common/restrl
derdrv.o:		common/runlab
derdrv.o:		common/segm
derdrv.o:		common/sizes
derdrv.o:		common/statis
derdrv.o:		common/symtry
derdrv.o:		common/tdhf
derdrv.o:		common/tdhfx
derdrv.o:		common/timeperiods
derdrv.o:		common/timez
derdrv.o:		common/vectrn
dft.o:		common/angxyz
dft.o:		common/dnfnw
dft.o:		common/field
dft.o:		common/infoa
dft.o:		common/infob
dft.o:		common/iofile
dft.o:		common/nshel
dft.o:		common/phycon
dft.o:		common/radxxx
dft.o:		common/radyyy
dft.o:		common/scfopt
dft.o:		common/scfwfn
dft.o:		common/sizes
dircta.o:		common/atmblk
dircta.o:		common/auxh
dircta.o:		common/blksiz
dircta.o:		common/blksiz2
dircta.o:		common/ccepa
dircta.o:		common/ccntl
dircta.o:		common/cdcryi
dircta.o:		common/cdcryz
dircta.o:		common/cic
dircta.o:		common/cjdavid
dircta.o:		common/cntrl
dircta.o:		common/comrjb
dircta.o:		common/corctl
dircta.o:		common/diactl
dircta.o:		common/dir_erg
dircta.o:		common/direc
dircta.o:		common/disc
dircta.o:		common/discc
dircta.o:		common/disktl
dircta.o:		common/disktl
dircta.o:		common/file12
dircta.o:		common/fpinfo
dircta.o:		common/helpr
dircta.o:		common/hold
dircta.o:		common/icon
dircta.o:		common/infoa
dircta.o:		common/iofile
dircta.o:		common/machin
dircta.o:		common/mapp
dircta.o:		common/modeci
dircta.o:		common/mr_cepctl
dircta.o:		common/mr_comci
dircta.o:		common/mr_mp2
dircta.o:		common/mr_pairs
dircta.o:		common/natorb
dircta.o:		common/orb
dircta.o:		common/presrt
dircta.o:		common/prints
dircta.o:		common/prnprn
dircta.o:		common/restar
dircta.o:		common/restri
dircta.o:		common/scfwfn
dircta.o:		common/segm
dircta.o:		common/sizes
dircta.o:		common/stak
dircta.o:		common/statis
dircta.o:		common/symchk
dircta.o:		common/symci
dircta.o:		common/table
dircta.o:		common/timanb
dircta.o:		common/timez
dircta.o:		common/trial
dircta.o:		common/work
dirctb.o:		../drf/comdrf/dafil
dirctb.o:		../drf/comdrf/darw
dirctb.o:		../drf/comdrf/drfpar
dirctb.o:		../drf/comdrf/sizesrf
dirctb.o:		common/atmblk
dirctb.o:		common/auxh
dirctb.o:		common/ccepa
dirctb.o:		common/ccntl
dirctb.o:		common/cdcryi
dirctb.o:		common/cdcryz
dirctb.o:		common/cic
dirctb.o:		common/cjdavid
dirctb.o:		common/cntrl
dirctb.o:		common/comrjb
dirctb.o:		common/corctl
dirctb.o:		common/diactl
dirctb.o:		common/dir_erg
dirctb.o:		common/disc
dirctb.o:		common/discc
dirctb.o:		common/disktl
dirctb.o:		common/disktl
dirctb.o:		common/drfopt
dirctb.o:		common/dump3
dirctb.o:		common/file12
dirctb.o:		common/ga_file
dirctb.o:		common/helpr
dirctb.o:		common/hold
dirctb.o:		common/iofile
dirctb.o:		common/machin
dirctb.o:		common/mapp
dirctb.o:		common/modeci
dirctb.o:		common/mr_cepctl
dirctb.o:		common/mr_mp2
dirctb.o:		common/mr_mp2l
dirctb.o:		common/mr_pairs
dirctb.o:		common/natorb
dirctb.o:		common/nodeio
dirctb.o:		common/orb
dirctb.o:		common/outctl
dirctb.o:		common/parallel
dirctb.o:		common/prints
dirctb.o:		common/prnprn
dirctb.o:		common/restar
dirctb.o:		common/sizes
dirctb.o:		common/symci
dirctb.o:		common/table
dirctb.o:		common/timanb
dirctb.o:		common/timez
dirctb.o:		common/trial
dirctc.o:		common/atmblk
dirctc.o:		common/auxh
dirctc.o:		common/auxh
dirctc.o:		common/blksiz
dirctc.o:		common/ccepa
dirctc.o:		common/ccepa
dirctc.o:		common/ccntl
dirctc.o:		common/ccntl
dirctc.o:		common/ccntl
dirctc.o:		common/cdcryi
dirctc.o:		common/cdcryz
dirctc.o:		common/cntrl
dirctc.o:		common/corctl
dirctc.o:		common/diactl
dirctc.o:		common/dir_erg
dirctc.o:		common/discc
dirctc.o:		common/disktl
dirctc.o:		common/disktl
dirctc.o:		common/file12
dirctc.o:		common/helpmp3
dirctc.o:		common/helpr
dirctc.o:		common/iofile
dirctc.o:		common/machin
dirctc.o:		common/machin
dirctc.o:		common/mapp
dirctc.o:		common/mr_cepctl
dirctc.o:		common/mr_comci
dirctc.o:		common/mr_mp2
dirctc.o:		common/mr_mp2con
dirctc.o:		common/mr_mp2l
dirctc.o:		common/mr_pairs
dirctc.o:		common/natorb
dirctc.o:		common/orb
dirctc.o:		common/presrt
dirctc.o:		common/prints
dirctc.o:		common/prnprn
dirctc.o:		common/sizes
dirctc.o:		common/symci
dirctc.o:		common/symci
dirctc.o:		common/table
dirctc.o:		common/trial
dirctc.o:		common/trial
dirctc.o:		common/work
dirctd.o:		common/ccepa
dirctd.o:		common/ccntl
dirctd.o:		common/cdcryz
dirctd.o:		common/cntrl
dirctd.o:		common/diactl
dirctd.o:		common/disktl
dirctd.o:		common/helpr
dirctd.o:		common/mr_cepctl
dirctd.o:		common/mr_pairs
dirctd.o:		common/sizes
dirctd.o:		common/symci
dirctd.o:		common/work
direct.o:		common/atmol3
direct.o:		common/cndx40
direct.o:		common/cndx41
direct.o:		common/dump3
direct.o:		common/funct
direct.o:		common/gmempara
direct.o:		common/harmon
direct.o:		common/infoa
direct.o:		common/iofile
direct.o:		common/machin
direct.o:		common/mapper
direct.o:		common/nshel
direct.o:		common/parallel
direct.o:		common/prints
direct.o:		common/restar
direct.o:		common/runlab
direct.o:		common/scra7
direct.o:		common/sizes
direct.o:		common/statis
direct.o:		common/symtry
direct.o:		common/tran
direct.o:		common/wrtd
dirrpa.o:		common/gmempara
dirrpa.o:		common/iofile
dirrpa.o:		common/machin
dirrpa.o:		common/nshel
dirrpa.o:		common/parallel
dirrpa.o:		common/polars
dirrpa.o:		common/prints
dirrpa.o:		common/response_timep
dirrpa.o:		common/restrj
dirrpa.o:		common/rpfile
dirrpa.o:		common/runlab
dirrpa.o:		common/sizes
dirrpa.o:		common/statis
dirrpa.o:		common/symtry
dirrpa.o:		common/timeperiods
drv1e.o:		../drf/comdrf/drfbem
drv1e.o:		../drf/comdrf/drfpar
drv1e.o:		../drf/comdrf/rfene
drv1e.o:		../drf/comdrf/rfgrad
drv1e.o:		../drf/comdrf/runpar
drv1e.o:		../drf/comdrf/sizesrf
drv1e.o:		common/atmblk
drv1e.o:		common/blur
drv1e.o:		common/ccpdft.hf77
drv1e.o:		common/chgcc
drv1e.o:		common/chmgms
drv1e.o:		common/cigrad
drv1e.o:		common/cndx40
drv1e.o:		common/cndx41
drv1e.o:		common/common
drv1e.o:		common/cslosc
drv1e.o:		common/dm
drv1e.o:		common/drfopt
drv1e.o:		common/dump3
drv1e.o:		common/funct
drv1e.o:		common/g80nb
drv1e.o:		common/ghfblk
drv1e.o:		common/gmempara
drv1e.o:		common/grad2
drv1e.o:		common/hermit
drv1e.o:		common/ijlab
drv1e.o:		common/infoa
drv1e.o:		common/iofile
drv1e.o:		common/machin
drv1e.o:		common/mapper
drv1e.o:		common/misc
drv1e.o:		common/modj
drv1e.o:		common/mp2grad_pointers
drv1e.o:		common/nodinf
drv1e.o:		common/nshel
drv1e.o:		common/prints
drv1e.o:		common/prnprn
drv1e.o:		common/repp
drv1e.o:		common/restar
drv1e.o:		common/restri
drv1e.o:		common/restrl
drv1e.o:		common/root
drv1e.o:		common/runlab
drv1e.o:		common/scfwfn
drv1e.o:		common/scra7
drv1e.o:		common/segm
drv1e.o:		common/sizes
drv1e.o:		common/specal
drv1e.o:		common/statis
drv1e.o:		common/symtry
drv1e.o:		common/timeperiods
drv1e.o:		common/timez
drv1e.o:		common/wermit
drv1e.o:		common/zorac
drv2e.o:		common/atmblk
drv2e.o:		common/ccpdft.hf77
drv2e.o:		common/chmgms
drv2e.o:		common/cigrad
drv2e.o:		common/cndx41
drv2e.o:		common/cslosc
drv2e.o:		common/dmisc
drv2e.o:		common/dshlnf
drv2e.o:		common/dshlno
drv2e.o:		common/dump3
drv2e.o:		common/funct
drv2e.o:		common/ghfblk
drv2e.o:		common/gmempara
drv2e.o:		common/grad2
drv2e.o:		common/ijlab
drv2e.o:		common/incrd
drv2e.o:		common/infoa
drv2e.o:		common/iofile
drv2e.o:		common/mapper
drv2e.o:		common/nshel
drv2e.o:		common/parallel
drv2e.o:		common/picon
drv2e.o:		common/prints
drv2e.o:		common/prnprn
drv2e.o:		common/restar
drv2e.o:		common/restri
drv2e.o:		common/restrl
drv2e.o:		common/root
drv2e.o:		common/rpadcom
drv2e.o:		common/runlab
drv2e.o:		common/runopt
drv2e.o:		common/scfwfn
drv2e.o:		common/sector
drv2e.o:		common/segm
drv2e.o:		common/sizes
drv2e.o:		common/specal
drv2e.o:		common/symtry
drv2e.o:		common/timez
drv2e.o:		common/vcore
drv80.o:		common/auxvar
drv80.o:		common/ccpdft.hf77
drv80.o:		common/cslosc
drv80.o:		common/dindex
drv80.o:		common/dmisc
drv80.o:		common/ffq
drv80.o:		common/funct
drv80.o:		common/grad2
drv80.o:		common/ijlab
drv80.o:		common/infoa
drv80.o:		common/iofile
drv80.o:		common/mapper
drv80.o:		common/nshel
drv80.o:		common/parallel
drv80.o:		common/prints
drv80.o:		common/restar
drv80.o:		common/scfwfn
drv80.o:		common/sizes
drv80.o:		common/specal
drv80.o:		common/symtry
drv80.o:		common/timez
drvmp.o:		common/atmblk
drvmp.o:		common/blksiz
drvmp.o:		common/cigrad
drvmp.o:		common/cndx41
drvmp.o:		common/common
drvmp.o:		common/infoa
drvmp.o:		common/iofile
drvmp.o:		common/mapper
drvmp.o:		common/nshel
drvmp.o:		common/sizes
drvmp.o:		common/stak
drvmp.o:		common/symtry
drvmp.o:		common/three
fullci.o:		common/atmblk
fullci.o:		common/atmol3
fullci.o:		common/blkorbs
fullci.o:		common/cciab
fullci.o:		common/ci
fullci.o:		common/cone
fullci.o:		common/cstrng
fullci.o:		common/direc
fullci.o:		common/discc
fullci.o:		common/filel
fullci.o:		common/harmon
fullci.o:		common/infoa
fullci.o:		common/iofile
fullci.o:		common/machin
fullci.o:		common/mapper
fullci.o:		common/prnprn
fullci.o:		common/restar
fullci.o:		common/restri
fullci.o:		common/runlab
fullci.o:		common/savem
fullci.o:		common/scfwfn
fullci.o:		common/sector
fullci.o:		common/sizes
fullci.o:		common/statis
fullci.o:		common/timez
fullci.o:		common/tran
fullci.o:		common/work
ga.o:		common/disc
ga.o:		common/discc
ga.o:		common/errcodes
ga.o:		common/ga_file
ga.o:		common/gaanal
ga.o:		common/gadiis
ga.o:		common/gmempara
ga.o:		common/harmon
ga.o:		common/infoa
ga.o:		common/iofile
ga.o:		common/mapper
ga.o:		common/nodeio
ga.o:		common/parcntl
ga.o:		common/prints
ga.o:		common/scfopt
ga.o:		common/scra7
ga.o:		common/sizes
ga.o:		common/timeperiods
ga.o:		common/vcore
gamess_hondo.o:		../drf/comdrf/darw
gamess_hondo.o:		../drf/comdrf/sizesrf
gamess_hondo.o:		../m4/common/disc
gamess_hondo.o:		../m4/common/discc
gamess_hondo.o:		../m4/common/sizes
gamess_hondo.o:		../m4/common/utilc
gamess_hondo.o:		common/restar
gdiesel.o:		common/diesel
gdiesel.o:		common/iofile
gff.o:		common/atmblk
gff.o:		common/discc
gff.o:		common/filel
gff.o:		common/iofile
gff.o:		common/machin
gff.o:		common/prints
gff.o:		common/restar
gff.o:		common/restri
gff.o:		common/runlab
gff.o:		common/sizes
gff.o:		common/statis
gff.o:		common/work
glue_uk_us.o:		common/global
glue_uk_us.o:		common/iofile
glue_uk_us.o:		common/molsym
glue_uk_us.o:		common/scfopt
glue_uk_us.o:		common/sizes
glue_uk_us.o:		common/symtry
guess.o:		../drf/comdrf/dafil
guess.o:		../drf/comdrf/drfpar
guess.o:		common/atmol3
guess.o:		common/blockc
guess.o:		common/ccpdft.hf77
guess.o:		common/cslosc
guess.o:		common/datgue
guess.o:		common/datgue
guess.o:		common/discc
guess.o:		common/drfopt
guess.o:		common/dump3
guess.o:		common/gadiis
guess.o:		common/gjs
guess.o:		common/gmempara
guess.o:		common/harmon
guess.o:		common/hermit
guess.o:		common/infoa
guess.o:		common/infob
guess.o:		common/iofile
guess.o:		common/iofile
guess.o:		common/llvmo
guess.o:		common/machin
guess.o:		common/mapper
guess.o:		common/multic
guess.o:		common/nodeio
guess.o:		common/nshel
guess.o:		common/parcntl
guess.o:		common/phycon
guess.o:		common/prints
guess.o:		common/restar
guess.o:		common/restri
guess.o:		common/restrj
guess.o:		common/restrl
guess.o:		common/runlab
guess.o:		common/scfopt
guess.o:		common/scfwfn
guess.o:		common/scra7
guess.o:		common/sector
guess.o:		common/segm
guess.o:		common/sizes
guess.o:		common/sizes
guess.o:		common/statis
guess.o:		common/syminf
guess.o:		common/symtry
guess.o:		common/timeperiods
guess.o:		common/timez
guess.o:		common/tran
guess.o:		common/tranao
guess.o:		common/wermit
guess.o:		common/zorac
index4.o:		common/atmblk
index4.o:		common/blksiz
index4.o:		common/cigrad
index4.o:		common/cndx40
index4.o:		common/cndx41
index4.o:		common/common
index4.o:		common/ghfblk
index4.o:		common/harmon
index4.o:		common/infoa
index4.o:		common/iofile
index4.o:		common/machin
index4.o:		common/mapper
index4.o:		common/nshel
index4.o:		common/prints
index4.o:		common/prnprn
index4.o:		common/restar
index4.o:		common/restrj
index4.o:		common/restrl
index4.o:		common/runlab
index4.o:		common/sizes
index4.o:		common/stak
index4.o:		common/symtry
index4.o:		common/three
index4.o:		common/timez
index4.o:		common/tran
index4.o:		common/uhfspn
index4.o:		common/vectrn
index4.o:		common/zorac
inpci.o:		common/avstat
inpci.o:		common/blksiz
inpci.o:		common/cepa_mrd1
inpci.o:		common/cepa_mrd2
inpci.o:		common/ciconv
inpci.o:		common/cigrad
inpci.o:		common/cjdavid
inpci.o:		common/cndx40
inpci.o:		common/comrjb2
inpci.o:		common/ctrl
inpci.o:		common/direc
inpci.o:		common/discc
inpci.o:		common/disktl
inpci.o:		common/dm
inpci.o:		common/dump3
inpci.o:		common/exc
inpci.o:		common/files
inpci.o:		common/finish
inpci.o:		common/fpinfo
inpci.o:		common/gjs
inpci.o:		common/harmon
inpci.o:		common/infoa
inpci.o:		common/iofile
inpci.o:		common/jobopt
inpci.o:		common/limy
inpci.o:		common/machin
inpci.o:		common/mapper
inpci.o:		common/mrdcid
inpci.o:		common/multic
inpci.o:		common/natorb
inpci.o:		common/newmrd_sort
inpci.o:		common/prints
inpci.o:		common/restar
inpci.o:		common/restrj
inpci.o:		common/runlab
inpci.o:		common/simul
inpci.o:		common/sizes
inpci.o:		common/statis
inpci.o:		common/syminf
inpci.o:		common/timez
inpci.o:		common/work
inpci.o:		common/workc
inpci.o:		common/zorac
input.o:		../drf/comdrf/trnsf
input.o:		common/atmol3
input.o:		common/blur
input.o:		common/canon_nat
input.o:		common/ccpdft.hf77
input.o:		common/chmgms
input.o:		common/cndx41
input.o:		common/cntl1
input.o:		common/common
input.o:		common/coval
input.o:		common/crnams
input.o:		common/crnamx
input.o:		common/csubch
input.o:		common/csubst
input.o:		common/czmat
input.o:		common/datbas
input.o:		common/datgue
input.o:		common/direc
input.o:		common/discc
input.o:		common/dump3
input.o:		common/errcodes
input.o:		common/fiddle_hesc
input.o:		common/foropt
input.o:		common/fpinfo
input.o:		common/frame
input.o:		common/fsymas
input.o:		common/ghfblk
input.o:		common/gmempara
input.o:		common/harmon
input.o:		common/infoa
input.o:		common/infob
input.o:		common/iofile
input.o:		common/machin
input.o:		common/mapper
input.o:		common/mod2d
input.o:		common/moda
input.o:		common/modf
input.o:		common/modj
input.o:		common/modmd
input.o:		common/modmin
input.o:		common/molsym
input.o:		common/mp2grd
input.o:		common/nbterm
input.o:		common/nshel
input.o:		common/phycon
input.o:		common/prints
input.o:		common/prnprn
input.o:		common/qercmx
input.o:		common/qercom
input.o:		common/restar
input.o:		common/restri
input.o:		common/restrj
input.o:		common/restrl
input.o:		common/restrz
input.o:		common/rtdata
input.o:		common/runlab
input.o:		common/runopt
input.o:		common/scfopt
input.o:		common/scfwfn
input.o:		common/scra7
input.o:		common/secd_pointers
input.o:		common/sector
input.o:		common/sizes
input.o:		common/structchk
input.o:		common/symtry
input.o:		common/tdhf
input.o:		common/tdhfx
input.o:		common/tol
input.o:		common/transf
input.o:		common/work
input.o:		common/workc
input.o:		common/xmlin
input.o:		common/zorac
intega.o:		../drf/comdrf/dafil
intega.o:		../drf/comdrf/darw
intega.o:		../drf/comdrf/drfbem
intega.o:		../drf/comdrf/drfdaf
intega.o:		../drf/comdrf/drfexp
intega.o:		../drf/comdrf/drfint
intega.o:		../drf/comdrf/drfpar
intega.o:		../drf/comdrf/ijpair
intega.o:		../drf/comdrf/int
intega.o:		../drf/comdrf/iofil
intega.o:		../drf/comdrf/mem
intega.o:		../drf/comdrf/opt
intega.o:		../drf/comdrf/runpar
intega.o:		../drf/comdrf/scfopt
intega.o:		../drf/comdrf/scmprt
intega.o:		../drf/comdrf/sizesrf
intega.o:		common/atmblk
intega.o:		common/bshell
intega.o:		common/ccpdft.hf77
intega.o:		common/chgcc
intega.o:		common/cndx41
intega.o:		common/common
intega.o:		common/cslosc
intega.o:		common/denss
intega.o:		common/discc
intega.o:		common/drfopt
intega.o:		common/dump3
intega.o:		common/files
intega.o:		common/flip70
intega.o:		common/flips
intega.o:		common/g80nb
intega.o:		common/gmempara
intega.o:		common/hermit
intega.o:		common/ijlab
intega.o:		common/imisc
intega.o:		common/indez
intega.o:		common/infoa
intega.o:		common/inxblk
intega.o:		common/iofile
intega.o:		common/machin
intega.o:		common/mapper
intega.o:		common/maxlen
intega.o:		common/misc
intega.o:		common/modj
intega.o:		common/morokuma
intega.o:		common/newscf
intega.o:		common/nshel
intega.o:		common/parallel
intega.o:		common/parcntl
intega.o:		common/picon
intega.o:		common/pkfil
intega.o:		common/prints
intega.o:		common/prnprn
intega.o:		common/restar
intega.o:		common/restri
intega.o:		common/restrl
intega.o:		common/root
intega.o:		common/rtdata
intega.o:		common/runlab
intega.o:		common/scra7
intega.o:		common/segm
intega.o:		common/setint
intega.o:		common/shlg70
intega.o:		common/shlinf
intega.o:		common/shlnos
intega.o:		common/shlt
intega.o:		common/sizes
intega.o:		common/sortp
intega.o:		common/statis
intega.o:		common/symtry
intega.o:		common/timez
intega.o:		common/utilc
intega.o:		common/wermit
intega.o:		common/xfield
intega.o:		common/zorac
integb.o:		common/cigrad
integb.o:		common/common
integb.o:		common/dawfcm
integb.o:		common/derfcm
integb.o:		common/ecp1
integb.o:		common/ecp3
integb.o:		common/ecp4
integb.o:		common/ecpdim
integb.o:		common/gbase
integb.o:		common/infoa
integb.o:		common/iofile
integb.o:		common/mapper
integb.o:		common/maxlen
integb.o:		common/misc
integb.o:		common/nshel
integb.o:		common/picon
integb.o:		common/prnprn
integb.o:		common/restar
integb.o:		common/restri
integb.o:		common/runlab
integb.o:		common/sizes
integb.o:		common/symtry
integb_lib.o:		common/direc
integb_lib.o:		common/discc
integb_lib.o:		common/ecpdim
integb_lib.o:		common/infoa
integb_lib.o:		common/iofile
integb_lib.o:		common/machin
integb_lib.o:		common/nshel
integb_lib.o:		common/repp
integb_lib.o:		common/restri
integb_lib.o:		common/runlab
integb_lib.o:		common/sizes
integb_lib.o:		common/work
integb_lib.o:		common/workc
integb_nl.o:		common/infoa
integb_nl.o:		common/iofile
integb_nl.o:		common/mapper
integb_nl.o:		common/nshel
integb_nl.o:		common/rcombd
integb_nl.o:		common/restar
integb_nl.o:		common/restrl
integb_nl.o:		common/runlab
integb_nl.o:		common/sizes
integb_nl.o:		common/symtry
integb_nl.o:		common/timez
integc.o:		common/atmblk
integc.o:		common/bshell
integc.o:		common/const
integc.o:		common/flip70
integc.o:		common/geom
integc.o:		common/ginf
integc.o:		common/gmempara
integc.o:		common/ijlab
integc.o:		common/infoa
integc.o:		common/iofile
integc.o:		common/mapper
integc.o:		common/maxc
integc.o:		common/miscg
integc.o:		common/nshel
integc.o:		common/parallel
integc.o:		common/pgeom
integc.o:		common/picon
integc.o:		common/pkfil
integc.o:		common/qgeom
integc.o:		common/restar
integc.o:		common/restri
integc.o:		common/shlg70
integc.o:		common/shllfo
integc.o:		common/shlnos
integc.o:		common/shlt
integc.o:		common/sizes
integc.o:		common/sortp
integc.o:		common/symtry
integc.o:		common/timez
integd.o:		common/atmblk
integd.o:		common/cigrad
integd.o:		common/cndx41
integd.o:		common/common
integd.o:		common/infoa
integd.o:		common/iofile
integd.o:		common/mapper
integd.o:		common/nshel
integd.o:		common/prnprn
integd.o:		common/runlab
integd.o:		common/sizes
integd.o:		common/symtry
integd.o:		common/timez
intege.o:		common/atmblk
intege.o:		common/blksiz
intege.o:		common/mapper
intege.o:		common/sizes
intege.o:		common/three
integs.o:		common/astore
integs.o:		common/auxvar
integs.o:		common/bshell
integs.o:		common/ccpdft.hf77
integs.o:		common/confus
integs.o:		common/const
integs.o:		common/cslosc
integs.o:		common/flip70
integs.o:		common/geom
integs.o:		common/ginf
integs.o:		common/gmempara
integs.o:		common/ijlab
integs.o:		common/indez
integs.o:		common/infoa
integs.o:		common/iofile
integs.o:		common/machin
integs.o:		common/mapper
integs.o:		common/maxc
integs.o:		common/misc
integs.o:		common/miscg
integs.o:		common/morokuma
integs.o:		common/newscf
integs.o:		common/nshel
integs.o:		common/parallel
integs.o:		common/pgeom
integs.o:		common/picon
integs.o:		common/pkfil
integs.o:		common/pqgeom
integs.o:		common/qgeom
integs.o:		common/restar
integs.o:		common/restri
integs.o:		common/runlab
integs.o:		common/scfopt
integs.o:		common/segm
integs.o:		common/shlg70
integs.o:		common/shllfo
integs.o:		common/shlnos
integs.o:		common/shlt
integs.o:		common/sizes
integs.o:		common/sortp
integs.o:		common/symtry
integs.o:		common/timez
integs.o:		common/utilc
inter_vdwaals.o:		common/infoa
inter_vdwaals.o:		common/iofile
inter_vdwaals.o:		common/runlab
inter_vdwaals.o:		common/sizes
inter_vdwaals.o:		vdwaals_corr.o
machci.o:		common/atmblk
machci.o:		common/blksiz
machci.o:		common/blksiz2
machci.o:		common/cfsinfo
machci.o:		common/iofile
machci.o:		common/nodinf
machci.o:		common/sizes
machci.o:		common/stak
machci.o:		common/utilc
machscf.o:		common/atmblk
machscf.o:		common/blksiz
machscf.o:		common/blocks
machscf.o:		common/cfsinfo
machscf.o:		common/chemshell
machscf.o:		common/direc
machscf.o:		common/disc
machscf.o:		common/discc
machscf.o:		common/discc1
machscf.o:		common/errcodes
machscf.o:		common/files
machscf.o:		common/fock_critic
machscf.o:		common/ga_file
machscf.o:		common/gmempara
machscf.o:		common/infoa
machscf.o:		common/iofile
machscf.o:		common/ipsc
machscf.o:		common/jinfo
machscf.o:		common/machin
machscf.o:		common/maxlen
machscf.o:		common/nodeio
machscf.o:		common/nodinf
machscf.o:		common/parcntl
machscf.o:		common/prints
machscf.o:		common/restar
machscf.o:		common/segm
machscf.o:		common/sizes
machscf.o:		common/sizes
machscf.o:		common/symtry
machscf.o:		common/syscod
machscf.o:		common/taskfarm
machscf.o:		common/timanb
machscf.o:		common/timeperiods
machscf.o:		common/timez
machscf.o:		common/utilc
machscf.o:		common/vectrn
machscf.o:		common/work
machscf.o:		common/workc
machscf.o:		common/xmlin
machscf.o:		common/xmlout
mainci.o:		common/auxh
mainci.o:		common/avstat
mainci.o:		common/caspar
mainci.o:		common/ccepa
mainci.o:		common/cdcryi
mainci.o:		common/cic
mainci.o:		common/ciconv
mainci.o:		common/comrjb
mainci.o:		common/corctl
mainci.o:		common/ctrl
mainci.o:		common/diactl
mainci.o:		common/dm
mainci.o:		common/finish
mainci.o:		common/helpr
mainci.o:		common/hold
mainci.o:		common/jobopt
mainci.o:		common/mcff
mainci.o:		common/multic
mainci.o:		common/natorb
mainci.o:		common/outctl
mainci.o:		common/presrt
mainci.o:		common/qice
mainci.o:		common/simul
mainci.o:		common/sizes
mainci.o:		common/symci
mainci.o:		common/syminf
mainci.o:		common/timanb
mainci.o:		common/trial
mains.o:		common/astore
mains.o:		common/atmblk
mains.o:		common/atmol3
mains.o:		common/auxvar
mains.o:		common/blockc
mains.o:		common/blocks
mains.o:		common/bshell
mains.o:		common/canon_nat
mains.o:		common/chemshell
mains.o:		common/chgcc
mains.o:		common/chmgms
mains.o:		common/cndx40
mains.o:		common/cndx41
mains.o:		common/cntl1
mains.o:		common/cntl2
mains.o:		common/const
mains.o:		common/copt
mains.o:		common/coval
mains.o:		common/crnams
mains.o:		common/crnamx
mains.o:		common/csubch
mains.o:		common/csubst
mains.o:		common/datbas
mains.o:		common/datgue
mains.o:		common/datplt
mains.o:		common/datprp
mains.o:		common/defunk
mains.o:		common/denss
mains.o:		common/disc
mains.o:		common/discc
mains.o:		common/dmisc
mains.o:		common/drfopt
mains.o:		common/errcodes
mains.o:		common/extsav
mains.o:		common/ffq
mains.o:		common/fiddle_hesc
mains.o:		common/field
mains.o:		common/filel
mains.o:		common/files
mains.o:		common/fock_critic
mains.o:		common/foropt
mains.o:		common/fpinfo
mains.o:		common/frame
mains.o:		common/fsymas
mains.o:		common/funct
mains.o:		common/ghfblk
mains.o:		common/ginf
mains.o:		common/gjs
mains.o:		common/gmempara
mains.o:		common/grad2
mains.o:		common/gvalue
mains.o:		common/harmon
mains.o:		common/hermit
mains.o:		common/ijlab
mains.o:		common/imisc
mains.o:		common/indez
mains.o:		common/indsyx
mains.o:		common/indxsv
mains.o:		common/infoa
mains.o:		common/infob
mains.o:		common/inxblk
mains.o:		common/iofile
mains.o:		common/iotrp
mains.o:		common/jinfo
mains.o:		common/jspec
mains.o:		common/lattic
mains.o:		common/linkan
mains.o:		common/machin
mains.o:		common/mapper
mains.o:		common/maxlen
mains.o:		common/misc
mains.o:		common/miscg
mains.o:		common/modfov
mains.o:		common/modj
mains.o:		common/molsym
mains.o:		common/mpidata
mains.o:		common/nshel
mains.o:		common/overrule
mains.o:		common/periodic
mains.o:		common/phycon
mains.o:		common/picon
mains.o:		common/prints
mains.o:		common/qercmx
mains.o:		common/qercom
mains.o:		common/rcombd
mains.o:		common/restar
mains.o:		common/restri
mains.o:		common/restrj
mains.o:		common/rtdata
mains.o:		common/runlab
mains.o:		common/scra7
mains.o:		common/scrf
mains.o:		common/sector
mains.o:		common/seerch
mains.o:		common/segm
mains.o:		common/shlg70
mains.o:		common/shlinf
mains.o:		common/shlt
mains.o:		common/sizes
mains.o:		common/statis
mains.o:		common/taskfarm
mains.o:		common/tdhf
mains.o:		common/tdhfx
mains.o:		common/timez
mains.o:		common/tol
mains.o:		common/tranao
mains.o:		common/transf
mains.o:		common/trntim
mains.o:		common/unocas
mains.o:		common/utilc
mains.o:		common/vcore
mains.o:		common/vdwrad
mains.o:		common/vibrtn
mains.o:		common/wermit
mains.o:		common/work
mains.o:		common/workc
mains.o:		common/worksp
mains.o:		common/xmlin
mains.o:		common/xmlout
mains.o:		common/zjorgs
mains.o:		common/zlabs
mains.o:		common/zorac
mass.o:		common/errcodes
mass.o:		common/infoa
mass.o:		common/infob
mass.o:		common/iofile
mass.o:		common/masses
mass.o:		common/periodic
mass.o:		common/runlab
mass.o:		common/sizes
mass.o:		common/work
masscf.o:		common/atmol3
masscf.o:		common/cigrad
masscf.o:		common/cndx41
masscf.o:		common/cntl1
masscf.o:		common/cntl2
masscf.o:		common/common
masscf.o:		common/cslosc
masscf.o:		common/direc
masscf.o:		common/disp
masscf.o:		common/dump3
masscf.o:		common/fccwfn
masscf.o:		common/fsymas
masscf.o:		common/global
masscf.o:		common/gmempara
masscf.o:		common/harmon
masscf.o:		common/ijlab
masscf.o:		common/infoa
masscf.o:		common/iofile
masscf.o:		common/machin
masscf.o:		common/mapper
masscf.o:		common/modj
masscf.o:		common/mp2grad_pointers
masscf.o:		common/mp2grd
masscf.o:		common/nshel
masscf.o:		common/parcntl
masscf.o:		common/prints
masscf.o:		common/prnprn
masscf.o:		common/psscrf
masscf.o:		common/restar
masscf.o:		common/restrj
masscf.o:		common/runlab
masscf.o:		common/scfopt
masscf.o:		common/scra7
masscf.o:		common/segm
masscf.o:		common/sizes
masscf.o:		common/sortp
masscf.o:		common/symtry
masscf.o:		common/timeperiods
masscf.o:		common/timez
masscf.o:		common/tran
masscf.o:		common/vcore
masscf.o:		common/vectrn
masscf.o:		common/work
master.o:		../drf/comdrf/ciopt
master.o:		../drf/comdrf/dafil
master.o:		../drf/comdrf/darw
master.o:		../drf/comdrf/dppop
master.o:		../drf/comdrf/drfbem
master.o:		../drf/comdrf/drfdaf
master.o:		../drf/comdrf/drfpar
master.o:		../drf/comdrf/drfzfa
master.o:		../drf/comdrf/enrgci
master.o:		../drf/comdrf/extinf
master.o:		../drf/comdrf/rfene
master.o:		../drf/comdrf/runpar
master.o:		../drf/comdrf/sizesrf
master.o:		../drf/comdrf/trnsf
master.o:		common/atmblk
master.o:		common/atmol3
master.o:		common/auxvar
master.o:		common/blksiz
master.o:		common/blocks
master.o:		common/blur
master.o:		common/ccpdft.hf77
master.o:		common/chgcc
master.o:		common/chmgms
master.o:		common/ciconv
master.o:		common/cigrad
master.o:		common/cndx40
master.o:		common/cndx41
master.o:		common/cntl1
master.o:		common/common
master.o:		common/convert
master.o:		common/copt
master.o:		common/coptbas
master.o:		common/crnams
master.o:		common/crnamx
master.o:		common/cslosc
master.o:		common/csubst
master.o:		common/ctrl
master.o:		common/czmat
master.o:		common/datgue
master.o:		common/datprp
master.o:		common/direc
master.o:		common/disc
master.o:		common/discc
master.o:		common/disp
master.o:		common/dlfind
master.o:		common/dm
master.o:		common/dnfnw
master.o:		common/drfopt
master.o:		common/drive_dft
master.o:		common/dump3
master.o:		common/errcodes
master.o:		common/fccwfn
master.o:		common/fermidirac
master.o:		common/field
master.o:		common/files
master.o:		common/finish
master.o:		common/foropt
master.o:		common/fpinfo
master.o:		common/fsymas
master.o:		common/funcon
master.o:		common/funct
master.o:		common/g80nb
master.o:		common/gen3
master.o:		common/ghfblk
master.o:		common/ghost
master.o:		common/gjs
master.o:		common/gmempara
master.o:		common/groups
master.o:		common/gvalue
master.o:		common/harmon
master.o:		common/harmon
master.o:		common/indsyx
master.o:		common/indxsv
master.o:		common/infoa
master.o:		common/infob
master.o:		common/iofile
master.o:		common/jinfo
master.o:		common/jspec
master.o:		common/junkc
master.o:		common/lattic
master.o:		common/limy
master.o:		common/machin
master.o:		common/mapper
master.o:		common/maxlen
master.o:		common/mod2d
master.o:		common/moda
master.o:		common/modf
master.o:		common/modfov
master.o:		common/modj
master.o:		common/modmd
master.o:		common/modmin
master.o:		common/molsym
master.o:		common/morokuma
master.o:		common/mp2grd
master.o:		common/mpcut
master.o:		common/natorb
master.o:		common/nbterm
master.o:		common/newmrd_sort
master.o:		common/newscf
master.o:		common/nshel
master.o:		common/occnr_z
master.o:		common/orb
master.o:		common/overrule
master.o:		common/parcntl
master.o:		common/periodic
master.o:		common/phycon
master.o:		common/picon
master.o:		common/prints
master.o:		common/prnprn
master.o:		common/psscrf
master.o:		common/qercmx
master.o:		common/qercom
master.o:		common/repp
master.o:		common/restar
master.o:		common/restri
master.o:		common/restrj
master.o:		common/restrl
master.o:		common/root
master.o:		common/runlab
master.o:		common/runopt
master.o:		common/scfopt
master.o:		common/scfwfn
master.o:		common/scra7
master.o:		common/scrf
master.o:		common/sector
master.o:		common/seerch
master.o:		common/segm
master.o:		common/simul
master.o:		common/sizes
master.o:		common/sortp
master.o:		common/statis
master.o:		common/structchk
master.o:		common/symtry
master.o:		common/taskfarm
master.o:		common/tdhf
master.o:		common/tdhfx
master.o:		common/timeperiods
master.o:		common/timez
master.o:		common/tol
master.o:		common/tran
master.o:		common/trntim
master.o:		common/unocas
master.o:		common/utilc
master.o:		common/vectrn
master.o:		common/work
master.o:		common/workc
master.o:		common/wrtd
master.o:		common/xfield
master.o:		common/xmlin
master.o:		common/xmlout
master.o:		common/zjorgs
master.o:		common/zorac
mcdab_ga.o:		common/cigrad
mcdab_ga.o:		common/dmisc
mcdab_ga.o:		common/incrd
mcdab_ga.o:		common/nshel
mcdab_ga.o:		common/sizes
mclr.o:		common/blocko
mclr.o:		common/debug
mclr.o:		common/iofile
mclr.o:		common/lradrs
mclr.o:		common/lrfile
mclr.o:		common/machin
mclr.o:		common/mcff
mclr.o:		common/mclr_blocks
mclr.o:		common/mcrun
mclr.o:		common/optima
mclr.o:		common/polars
mclr.o:		common/restar
mclr.o:		common/restri
mclr.o:		common/sector
mclr.o:		common/sizes
mclr.o:		common/tran
mclr.o:		common/work
mcscfa.o:		common/atmblk
mcscfa.o:		common/count
mcscfa.o:		common/detcic
mcscfa.o:		common/direc
mcscfa.o:		common/discc
mcscfa.o:		common/dump3
mcscfa.o:		common/gjs
mcscfa.o:		common/harmon
mcscfa.o:		common/iofile
mcscfa.o:		common/jinfo
mcscfa.o:		common/jobopt
mcscfa.o:		common/machin
mcscfa.o:		common/mapper
mcscfa.o:		common/mcff
mcscfa.o:		common/mcscra
mcscfa.o:		common/multic
mcscfa.o:		common/prnprn
mcscfa.o:		common/runlab
mcscfa.o:		common/sizes
mcscfa.o:		common/stak
mcscfa.o:		common/syminf
mcscfb.o:		common/atmblk
mcscfb.o:		common/atmol3
mcscfb.o:		common/count
mcscfb.o:		common/detcic
mcscfb.o:		common/harmon
mcscfb.o:		common/infoa
mcscfb.o:		common/iofile
mcscfb.o:		common/jinfo
mcscfb.o:		common/jobopt
mcscfb.o:		common/machin
mcscfb.o:		common/mapper
mcscfb.o:		common/mcaddr
mcscfb.o:		common/mcff
mcscfb.o:		common/mcscra
mcscfb.o:		common/multic
mcscfb.o:		common/natorb
mcscfb.o:		common/prints
mcscfb.o:		common/prnprn
mcscfb.o:		common/restri
mcscfb.o:		common/runlab
mcscfb.o:		common/segm
mcscfb.o:		common/sizes
mcscfb.o:		common/syminf
mcscfc.o:		common/jobopt
mcscfc.o:		common/multic
mcscfc.o:		common/sizes
mcscfc.o:		common/syminf
model.o:		common/cntl1
model.o:		common/csubst
model.o:		common/finish
model.o:		common/funct
model.o:		common/g80nb
model.o:		common/gmempara
model.o:		common/infoa
model.o:		common/infob
model.o:		common/iofile
model.o:		common/mod2d
model.o:		common/moda
model.o:		common/modf
model.o:		common/modfov
model.o:		common/modj
model.o:		common/modmd
model.o:		common/modmin
model.o:		common/molsym
model.o:		common/nbterm
model.o:		common/phycon
model.o:		common/prints
model.o:		common/restar
model.o:		common/runhed
model.o:		common/runlab
model.o:		common/segm
model.o:		common/sizes
model.o:		common/statis
morokuma.o:		common/atmblk
morokuma.o:		common/atmol3
morokuma.o:		common/ccpdft.hf77
morokuma.o:		common/cslosc
morokuma.o:		common/dump3
morokuma.o:		common/extsav
morokuma.o:		common/filel
morokuma.o:		common/fsymas
morokuma.o:		common/funct
morokuma.o:		common/gadiis
morokuma.o:		common/gmempara
morokuma.o:		common/indez
morokuma.o:		common/infoa
morokuma.o:		common/iofile
morokuma.o:		common/mapper
morokuma.o:		common/maxlen
morokuma.o:		common/misc
morokuma.o:		common/morokuma
morokuma.o:		common/nodeio
morokuma.o:		common/nshel
morokuma.o:		common/parcntl
morokuma.o:		common/phycon
morokuma.o:		common/prints
morokuma.o:		common/prnprn
morokuma.o:		common/restar
morokuma.o:		common/restri
morokuma.o:		common/restrl
morokuma.o:		common/runlab
morokuma.o:		common/scfopt
morokuma.o:		common/scra7
morokuma.o:		common/shlg70
morokuma.o:		common/shlnos
morokuma.o:		common/shlt
morokuma.o:		common/sizes
morokuma.o:		common/sortp
morokuma.o:		common/statis
morokuma.o:		common/symtry
morokuma.o:		common/timeperiods
morokuma.o:		common/vcore
morokuma.o:		common/vcore
morokuma.o:		common/xfield
mp2_parallel.o:		common/atmblk
mp2_parallel.o:		common/atmol3
mp2_parallel.o:		common/cigrad
mp2_parallel.o:		common/cndx40
mp2_parallel.o:		common/cndx41
mp2_parallel.o:		common/cntl1
mp2_parallel.o:		common/cntl2
mp2_parallel.o:		common/common
mp2_parallel.o:		common/cslosc
mp2_parallel.o:		common/direc
mp2_parallel.o:		common/disp
mp2_parallel.o:		common/dmisc
mp2_parallel.o:		common/dshlno
mp2_parallel.o:		common/dump3
mp2_parallel.o:		common/errcodes
mp2_parallel.o:		common/flips
mp2_parallel.o:		common/fsymas
mp2_parallel.o:		common/ghfblk
mp2_parallel.o:		common/global
mp2_parallel.o:		common/gmempara
mp2_parallel.o:		common/harmon
mp2_parallel.o:		common/ijlab
mp2_parallel.o:		common/incrd
mp2_parallel.o:		common/indez
mp2_parallel.o:		common/infoa
mp2_parallel.o:		common/iofile
mp2_parallel.o:		common/machin
mp2_parallel.o:		common/mapper
mp2_parallel.o:		common/misc
mp2_parallel.o:		common/modj
mp2_parallel.o:		common/mp2grad_pointers
mp2_parallel.o:		common/mp2grd
mp2_parallel.o:		common/nshel
mp2_parallel.o:		common/parallel
mp2_parallel.o:		common/parcntl
mp2_parallel.o:		common/picon
mp2_parallel.o:		common/pkfil
mp2_parallel.o:		common/prints
mp2_parallel.o:		common/prnprn
mp2_parallel.o:		common/psscrf
mp2_parallel.o:		common/restar
mp2_parallel.o:		common/restri
mp2_parallel.o:		common/restrj
mp2_parallel.o:		common/restrl
mp2_parallel.o:		common/root
mp2_parallel.o:		common/runlab
mp2_parallel.o:		common/scfopt
mp2_parallel.o:		common/scra7
mp2_parallel.o:		common/segm
mp2_parallel.o:		common/shlg70
mp2_parallel.o:		common/shlinf
mp2_parallel.o:		common/shlnos
mp2_parallel.o:		common/shlt
mp2_parallel.o:		common/sizes
mp2_parallel.o:		common/sortp
mp2_parallel.o:		common/specal
mp2_parallel.o:		common/statis
mp2_parallel.o:		common/symtry
mp2_parallel.o:		common/timeperiods
mp2_parallel.o:		common/timez
mp2_parallel.o:		common/tran
mp2_parallel.o:		common/utilc
mp2_parallel.o:		common/vcore
mp2_parallel.o:		common/vectrn
mp2a.o:		common/atmblk
mp2a.o:		common/blksiz
mp2a.o:		common/cigrad
mp2a.o:		common/cndx40
mp2a.o:		common/cndx41
mp2a.o:		common/common
mp2a.o:		common/disp
mp2a.o:		common/infoa
mp2a.o:		common/iofile
mp2a.o:		common/mapper
mp2a.o:		common/nshel
mp2a.o:		common/prnprn
mp2a.o:		common/restrj
mp2a.o:		common/sizes
mp2a.o:		common/stak
mp2a.o:		common/statis
mp2a.o:		common/three
mp2a.o:		common/timeperiods
mp2a.o:		common/vectrn
mp2b.o:		common/atmblk
mp2b.o:		common/blksiz
mp2b.o:		common/cigrad
mp2b.o:		common/cndx40
mp2b.o:		common/cndx41
mp2b.o:		common/common
mp2b.o:		common/gmempara
mp2b.o:		common/infoa
mp2b.o:		common/iofile
mp2b.o:		common/mapper
mp2b.o:		common/nshel
mp2b.o:		common/prnprn
mp2b.o:		common/sizes
mp2b.o:		common/symtry
mp2b.o:		common/symtry
mp2b.o:		common/uhfspn
mp2c.o:		common/blksiz
mp2c.o:		common/iofile
mp2c.o:		common/prnprn
mp2c.o:		common/stak
mp2c.o:		common/three
mp2d.o:		common/atmblk
mp2d.o:		common/cigrad
mp2d.o:		common/cndx40
mp2d.o:		common/cndx41
mp2d.o:		common/common
mp2d.o:		common/infoa
mp2d.o:		common/iofile
mp2d.o:		common/mapper
mp2d.o:		common/nshel
mp2d.o:		common/sizes
mp2d.o:		common/symtry
mp2d.o:		common/tdhf
mp2d.o:		common/tdhfx
mp3.o:		common/atmblk
mp3.o:		common/cigrad
mp3.o:		common/cndx40
mp3.o:		common/cndx41
mp3.o:		common/common
mp3.o:		common/harmon
mp3.o:		common/infoa
mp3.o:		common/iofile
mp3.o:		common/mapper
mp3.o:		common/nshel
mp3.o:		common/prints
mp3.o:		common/restar
mp3.o:		common/sizes
mp3.o:		common/symtry
mp3.o:		common/uhfspn
mrdci1.o:		common/atmblk
mrdci1.o:		common/blksiz
mrdci1.o:		common/blockc
mrdci1.o:		common/discc
mrdci1.o:		common/dump3
mrdci1.o:		common/filel
mrdci1.o:		common/files
mrdci1.o:		common/ftape
mrdci1.o:		common/infoa
mrdci1.o:		common/iofile
mrdci1.o:		common/machin
mrdci1.o:		common/mapper
mrdci1.o:		common/prints
mrdci1.o:		common/restar
mrdci1.o:		common/sizes
mrdci1.o:		common/zorac
mrdci2.o:		common/atmol3
mrdci2.o:		common/blksiz
mrdci2.o:		common/blockc
mrdci2.o:		common/discc
mrdci2.o:		common/ftape
mrdci2.o:		common/infoa
mrdci2.o:		common/iofile
mrdci2.o:		common/machin
mrdci2.o:		common/mapper
mrdci2.o:		common/prints
mrdci2.o:		common/sizes
mrdci3.o:		common/iofile
mrdci3.o:		common/prints
mrdci3.o:		common/sizes
mrdci4.o:		common/discc
mrdci4.o:		common/infoa
mrdci4.o:		common/iofile
mrdci4.o:		common/jinfo
mrdci4.o:		common/machin
mrdci4.o:		common/mapper
mrdci4.o:		common/nshel
mrdci4.o:		common/prints
mrdci4.o:		common/runlab
mrdci4.o:		common/sizes
mrdci5.o:		common/iofile
mrdci5.o:		common/prints
mrdci6.o:		common/iofile
mrdci6.o:		common/prints
mrdci6.o:		common/sizes
mrdci7.o:		common/iofile
mrdci7.o:		common/prints
mrdci7.o:		common/runlab
mrdci7.o:		common/sizes
nbo.o:		common/atmol3
nbo.o:		common/common
nbo.o:		common/cslosc
nbo.o:		common/dump3
nbo.o:		common/infoa
nbo.o:		common/iofile
nbo.o:		common/limy
nbo.o:		common/mapper
nbo.o:		common/nshel
nbo.o:		common/phycon
nbo.o:		common/restar
nbo.o:		common/restri
nbo.o:		common/runlab
nbo.o:		common/scfopt
nbo.o:		common/scra7
nbo.o:		common/sector
nbo.o:		common/sizes
nbo.o:		common/tran
newmrd1.o:		common/atmblk
newmrd1.o:		common/cic
newmrd1.o:		common/cndx40
newmrd1.o:		common/discc
newmrd1.o:		common/file12
newmrd1.o:		common/filel
newmrd1.o:		common/harmon
newmrd1.o:		common/infoa
newmrd1.o:		common/iofile
newmrd1.o:		common/machin
newmrd1.o:		common/mapper
newmrd1.o:		common/mrdcid
newmrd1.o:		common/newmrd_parinc
newmrd1.o:		common/newmrd_sort
newmrd1.o:		common/orb
newmrd1.o:		common/prints
newmrd1.o:		common/restar
newmrd1.o:		common/runlab
newmrd1.o:		common/scfwfn
newmrd1.o:		common/sizes
newmrd1.o:		common/table
newmrd2.o:		common/iofile
newmrd2.o:		common/prints
newmrd3.o:		common/cepa_mrd1
newmrd3.o:		common/cepa_sn2
newmrd3.o:		common/comrjb2
newmrd3.o:		common/discc
newmrd3.o:		common/infoa
newmrd3.o:		common/iofile
newmrd3.o:		common/newmrd_parinc
newmrd3.o:		common/newmrd_sort
newmrd3.o:		common/prints
newmrd3.o:		common/sizes
newmrd3.o:		common/utilc
newmrd4.o:		common/cepa_mrd1
newmrd4.o:		common/cepa_sn1
newmrd4.o:		common/cepa_sn2
newmrd4.o:		common/iofile
newmrd4.o:		common/newmrd_rglob
newmrd5.o:		common/cepa_mrd1
newmrd5.o:		common/cepa_mrd2
newmrd5.o:		common/cepa_sn1
newmrd5.o:		common/comrjb2
newmrd5.o:		common/infoa
newmrd5.o:		common/iofile
newmrd5.o:		common/newmrd_adlinc
newmrd5.o:		common/newmrd_rglob
newmrd5.o:		common/prints
newmrd5.o:		common/sizes
newmrd6.o:		common/atmol3
newmrd6.o:		common/blockc
newmrd6.o:		common/discc
newmrd6.o:		common/dump3
newmrd6.o:		common/filel
newmrd6.o:		common/files
newmrd6.o:		common/ftape
newmrd6.o:		common/harmon
newmrd6.o:		common/infoa
newmrd6.o:		common/iofile
newmrd6.o:		common/jinfo
newmrd6.o:		common/machin
newmrd6.o:		common/mapper
newmrd6.o:		common/prints
newmrd6.o:		common/restar
newmrd6.o:		common/runlab
newmrd6.o:		common/scra7
newmrd6.o:		common/sizes
newmrd6.o:		common/tran
newmrd6.o:		common/zorac
nvccsd.o:		common/prnprn
nvccsd.o:		common/t_adata
nvccsd.o:		common/t_files
nvccsd.o:		common/t_int
nvccsd.o:		common/t_iter
nvccsd.o:		common/t_parm
omit.o:		common/errcodes
optef.o:		common/cntl1
optef.o:		common/cntl2
optef.o:		common/csubch
optef.o:		common/csubst
optef.o:		common/czmat
optef.o:		common/dump3
optef.o:		common/funct
optef.o:		common/infoa
optef.o:		common/infob
optef.o:		common/iofile
optef.o:		common/machin
optef.o:		common/overrule
optef.o:		common/prints
optef.o:		common/restar
optef.o:		common/restri
optef.o:		common/runlab
optef.o:		common/scfopt
optef.o:		common/sector
optef.o:		common/segm
optef.o:		common/sizes
optef.o:		common/statis
optef.o:		common/timez
optim.o:		../drf/comdrf/drfbem
optim.o:		../drf/comdrf/drfpar
optim.o:		../drf/comdrf/rfene
optim.o:		../drf/comdrf/runpar
optim.o:		../drf/comdrf/sizesrf
optim.o:		common/atmol3
optim.o:		common/ccpdft.hf77
optim.o:		common/chmgms
optim.o:		common/cigrad
optim.o:		common/cndx41
optim.o:		common/cntl1
optim.o:		common/cntl2
optim.o:		common/common
optim.o:		common/copt
optim.o:		common/cslosc
optim.o:		common/csubch
optim.o:		common/csubst
optim.o:		common/czmat
optim.o:		common/datgue
optim.o:		common/discc
optim.o:		common/drfopt
optim.o:		common/dump3
optim.o:		common/fccwfn
optim.o:		common/fiddle_hesc
optim.o:		common/foropt
optim.o:		common/fpinfo
optim.o:		common/funct
optim.o:		common/gmempara
optim.o:		common/harmon
optim.o:		common/infoa
optim.o:		common/infob
optim.o:		common/iofile
optim.o:		common/machin
optim.o:		common/mapper
optim.o:		common/modj
optim.o:		common/molsym
optim.o:		common/natorb
optim.o:		common/nshel
optim.o:		common/phycon
optim.o:		common/prints
optim.o:		common/prnprn
optim.o:		common/restar
optim.o:		common/restri
optim.o:		common/restrj
optim.o:		common/restrl
optim.o:		common/rpadcom
optim.o:		common/runlab
optim.o:		common/runopt
optim.o:		common/scfwfn
optim.o:		common/scra7
optim.o:		common/seerch
optim.o:		common/segm
optim.o:		common/sizes
optim.o:		common/structchk
optim.o:		common/symtry
optim.o:		common/taskfarm
optim.o:		common/timez
optim.o:		common/tran
optim.o:		common/vibrtn
optim.o:		common/work
optim.o:		common/zorac
ormas1.o:		common/fccwfn
ormas2.o:		common/atmblk
ormas2.o:		common/restar
ormas2.o:		common/restri
parallel.o:		../m4/common/mpidata
parallel.o:		../m4/common/mpistatus
parallel.o:		../m4/common/timeperiods
parallel.o:		common/cntrl
parallel.o:		common/disktl
parallel.o:		common/ga_file
parallel.o:		common/gmemcount
parallel.o:		common/gmemdata
parallel.o:		common/gmemnull
parallel.o:		common/gmempara
parallel.o:		common/gmempara
parallel.o:		common/iofile
parallel.o:		common/mpidata
parallel.o:		common/mpistatus
parallel.o:		common/nodinf
parallel.o:		common/parallel
parallel.o:		common/parcntl
parallel.o:		common/parcntl
parallel.o:		common/prnprn
parallel.o:		common/segm
parallel.o:		common/sizes
parallel.o:		common/taskinfo
parallel.o:		common/timeperiods
parallel.o:		common/vcore
parallel.o:		common/vcore
pdiag.o:		common/timeperiods
pdiag.o:		common/vcore
plot.o:		common/discc
plot.o:		common/infoa
plot.o:		common/machin
plot.o:		common/sector
plot.o:		common/sizes
ptest.o:		common/gadiis
ptest.o:		common/iofile
ptest.o:		common/parcntl
ptest.o:		common/vcore
rpa.o:		common/atmblk
rpa.o:		common/atmol3
rpa.o:		common/blocko
rpa.o:		common/direc
rpa.o:		common/dump3
rpa.o:		common/infoa
rpa.o:		common/iofile
rpa.o:		common/machin
rpa.o:		common/nshel
rpa.o:		common/polars
rpa.o:		common/restrj
rpa.o:		common/rpacom
rpa.o:		common/rpfile
rpa.o:		common/sector
rpa.o:		common/sizes
rpa.o:		common/timeperiods
rpa.o:		common/tran
rpa.o:		common/work
rpa.o:		common/workc
rpagrad.o:		common/atmblk
rpagrad.o:		common/atmol3
rpagrad.o:		common/cndx40
rpagrad.o:		common/cndx41
rpagrad.o:		common/cslosc
rpagrad.o:		common/funct
rpagrad.o:		common/gjs
rpagrad.o:		common/gmempara
rpagrad.o:		common/grad2
rpagrad.o:		common/infoa
rpagrad.o:		common/iofile
rpagrad.o:		common/machin
rpagrad.o:		common/mapper
rpagrad.o:		common/maxlen
rpagrad.o:		common/nshel
rpagrad.o:		common/restar
rpagrad.o:		common/restri
rpagrad.o:		common/restrj
rpagrad.o:		common/rpacom
rpagrad.o:		common/rpadcom
rpagrad.o:		common/runlab
rpagrad.o:		common/sector
rpagrad.o:		common/sizes
rpagrad.o:		common/symtry
rpagrad.o:		common/timez
rpagrad.o:		common/vectrn
rpagrad.o:		common/work
scf.o:		../drf/comdrf/drfbem
scf.o:		../drf/comdrf/drfdaf
scf.o:		../drf/comdrf/drfpar
scf.o:		../drf/comdrf/gsspar
scf.o:		../drf/comdrf/opt
scf.o:		../drf/comdrf/runpar
scf.o:		../drf/comdrf/sizesrf
scf.o:		common/atmblk
scf.o:		common/atmol3
scf.o:		common/blur
scf.o:		common/ccpdft.hf77
scf.o:		common/chmgms
scf.o:		common/cslosc
scf.o:		common/datgue
scf.o:		common/disc
scf.o:		common/dm
scf.o:		common/dnfnw
scf.o:		common/drfopt
scf.o:		common/drive_dft
scf.o:		common/dump3
scf.o:		common/extsav
scf.o:		common/fermidirac
scf.o:		common/field
scf.o:		common/filel
scf.o:		common/fock_critic
scf.o:		common/fpinfo
scf.o:		common/fsymas
scf.o:		common/funct
scf.o:		common/gadiis
scf.o:		common/gjs
scf.o:		common/gmempara
scf.o:		common/gvalue
scf.o:		common/harmon
scf.o:		common/ijlab
scf.o:		common/infoa
scf.o:		common/iofile
scf.o:		common/machin
scf.o:		common/mapper
scf.o:		common/maxlen
scf.o:		common/modj
scf.o:		common/nodeio
scf.o:		common/nshel
scf.o:		common/parcntl
scf.o:		common/phycon
scf.o:		common/prints
scf.o:		common/prnprn
scf.o:		common/psscrf
scf.o:		common/restar
scf.o:		common/restrj
scf.o:		common/restrl
scf.o:		common/runlab
scf.o:		common/scfopt
scf.o:		common/scfwfn
scf.o:		common/scra7
scf.o:		common/scrf
scf.o:		common/segm
scf.o:		common/sizes
scf.o:		common/sortp
scf.o:		common/statis
scf.o:		common/symtry
scf.o:		common/timeperiods
scf.o:		common/timez
scf.o:		common/tran
scf.o:		common/vcore
scf.o:		common/xfield
scf.o:		common/zorac
sec1e.o:		common/cigrad
sec1e.o:		common/common
sec1e.o:		common/infoa
sec1e.o:		common/iofile
sec1e.o:		common/lattic
sec1e.o:		common/mapper
sec1e.o:		common/maxlen
sec1e.o:		common/nshel
sec1e.o:		common/picon
sec1e.o:		common/prnprn
sec1e.o:		common/repp
sec1e.o:		common/root
sec1e.o:		common/sizes
sec1e.o:		common/symtry
sec2e.o:		common/atmblk
sec2e.o:		common/ccpdft.hf77
sec2e.o:		common/cigrad
sec2e.o:		common/cndx41
sec2e.o:		common/common
sec2e.o:		common/ddmisc
sec2e.o:		common/ddmisc
sec2e.o:		common/ddshli
sec2e.o:		common/ddshln
sec2e.o:		common/drive_dft
sec2e.o:		common/ghfblk
sec2e.o:		common/gmempara
sec2e.o:		common/grad2
sec2e.o:		common/incrdd
sec2e.o:		common/infoa
sec2e.o:		common/iofile
sec2e.o:		common/mapper
sec2e.o:		common/maxlen
sec2e.o:		common/nshel
sec2e.o:		common/parallel
sec2e.o:		common/phycon
sec2e.o:		common/prnprn
sec2e.o:		common/restar
sec2e.o:		common/restri
sec2e.o:		common/root
sec2e.o:		common/runlab
sec2e.o:		common/runopt
sec2e.o:		common/setdd
sec2e.o:		common/sizes
sec2e.o:		common/symtry
sec2e.o:		common/timez
secd_parallel.o:		common/atmblk
secd_parallel.o:		common/ccpdft.hf77
secd_parallel.o:		common/cigrad
secd_parallel.o:		common/cndx40
secd_parallel.o:		common/cndx41
secd_parallel.o:		common/common
secd_parallel.o:		common/debug
secd_parallel.o:		common/drive_dft
secd_parallel.o:		common/dump3
secd_parallel.o:		common/ghfblk
secd_parallel.o:		common/global
secd_parallel.o:		common/gmempara
secd_parallel.o:		common/infoa
secd_parallel.o:		common/iofile
secd_parallel.o:		common/mapper
secd_parallel.o:		common/mp2grad_pointers
secd_parallel.o:		common/mp2grd
secd_parallel.o:		common/nshel
secd_parallel.o:		common/prnprn
secd_parallel.o:		common/restrj
secd_parallel.o:		common/runopt
secd_parallel.o:		common/secd_pointers
secd_parallel.o:		common/segm
secd_parallel.o:		common/sizes
secd_parallel.o:		common/sortp
secd_parallel.o:		common/statis
secd_parallel.o:		common/symtry
secd_parallel.o:		common/timeperiods
secd_parallel.o:		common/timez
secd_parallel.o:		common/vcore
secd_parallel.o:		common/vectrn
secdchf.o:		common/cndx41
secdchf.o:		common/common
secdchf.o:		common/global
secdchf.o:		common/gmempara
secdchf.o:		common/iofile
secdchf.o:		common/prnprn
secdchf.o:		common/sizes
secdchf.o:		common/timeperiods
secdchf.o:		common/timez
secdchf.o:		common/vcore
secmp2.o:		common/atmblk
secmp2.o:		common/blksiz
secmp2.o:		common/cigrad
secmp2.o:		common/cndx41
secmp2.o:		common/common
secmp2.o:		common/infoa
secmp2.o:		common/iofile
secmp2.o:		common/mapper
secmp2.o:		common/nshel
secmp2.o:		common/prnprn
secmp2.o:		common/restri
secmp2.o:		common/sizes
secmp2.o:		common/stak
secmp2.o:		common/symtry
secmp2.o:		common/three
server.o:		../drf/comdrf/drfpar
server.o:		../drf/comdrf/enrgci
server.o:		../drf/comdrf/extinf
server.o:		../drf/comdrf/rad
server.o:		../drf/comdrf/rfene
server.o:		../drf/comdrf/sizesrf
server.o:		common/aimpac
server.o:		common/atmblk
server.o:		common/atmol3
server.o:		common/blocks
server.o:		common/blur
server.o:		common/ccpdft.hf77
server.o:		common/cndx40
server.o:		common/cndx41
server.o:		common/cntl1
server.o:		common/cntl2
server.o:		common/common
server.o:		common/coptbas
server.o:		common/coptbm
server.o:		common/coval
server.o:		common/csubch
server.o:		common/csubst
server.o:		common/czmat
server.o:		common/direc
server.o:		common/disc
server.o:		common/discc
server.o:		common/drfopt
server.o:		common/drfopt
server.o:		common/dump3
server.o:		common/errcodes
server.o:		common/files
server.o:		common/harmon
server.o:		common/indxsv
server.o:		common/infoa
server.o:		common/infob
server.o:		common/iofile
server.o:		common/machin
server.o:		common/mapper
server.o:		common/maxlen
server.o:		common/molsym
server.o:		common/nodeio
server.o:		common/nodinf
server.o:		common/nshel
server.o:		common/periodic
server.o:		common/phycon
server.o:		common/prints
server.o:		common/restar
server.o:		common/restri
server.o:		common/restrl
server.o:		common/runlab
server.o:		common/savem
server.o:		common/scfopt
server.o:		common/scfwfn
server.o:		common/sector
server.o:		common/segm
server.o:		common/sizes
server.o:		common/sortbeta
server.o:		common/statis
server.o:		common/symtry
server.o:		common/sysmo
server.o:		common/tran
server.o:		common/transf
server.o:		common/unocas
server.o:		common/work
server.o:		common/workc
server.o:		common/xmlout
sysmo.o:		../dft/common/ccpdft.hf77
sysmo.o:		common/atmblk
sysmo.o:		common/atmol3
sysmo.o:		common/cndx40
sysmo.o:		common/cndx41
sysmo.o:		common/common
sysmo.o:		common/dump3
sysmo.o:		common/file12
sysmo.o:		common/filel
sysmo.o:		common/files
sysmo.o:		common/gjs
sysmo.o:		common/harmon
sysmo.o:		common/infoa
sysmo.o:		common/iofile
sysmo.o:		common/machin
sysmo.o:		common/mapper
sysmo.o:		common/nshel
sysmo.o:		common/periodic
sysmo.o:		common/restar
sysmo.o:		common/restri
sysmo.o:		common/scra7
sysmo.o:		common/sizes
sysmo.o:		common/syminf
sysmo.o:		common/symtry
sysmo.o:		common/sysmo
sysmo.o:		common/tran
sysmo.o:		common/vectrn
sysmo.o:		common/zorac
taskfarm.o:		common/errcodes
taskfarm.o:		common/iofile
taskfarm.o:		common/maxlen
taskfarm.o:		common/sizes
taskfarm.o:		common/taskfarm
taskfarm.o:		common/taskmain
taskfarm.o:		common/work
taskfarm.o:		common/workc
tdaf.o:		common/atmblk
tdaf.o:		common/discc
tdaf.o:		common/filel
tdaf.o:		common/iofile
tdaf.o:		common/machin
tdaf.o:		common/restar
tdaf.o:		common/restri
tdaf.o:		common/runlab
tdaf.o:		common/segm
tdaf.o:		common/sizes
tdaf.o:		common/statis
tran4.o:		common/atmblk
tran4.o:		common/blksiz
tran4.o:		common/blocko
tran4.o:		common/cndx40
tran4.o:		common/discc
tran4.o:		common/dump3
tran4.o:		common/field
tran4.o:		common/filel
tran4.o:		common/files
tran4.o:		common/fpinfo
tran4.o:		common/gmempara
tran4.o:		common/harmon
tran4.o:		common/infoa
tran4.o:		common/iofile
tran4.o:		common/machin
tran4.o:		common/mlngth
tran4.o:		common/newmrd_sort
tran4.o:		common/prints
tran4.o:		common/prnprn
tran4.o:		common/restar
tran4.o:		common/restrj
tran4.o:		common/runlab
tran4.o:		common/scra7
tran4.o:		common/sizes
tran4.o:		common/statis
tran4.o:		common/timez
tran4.o:		common/tran
tran4.o:		common/trntim
tran4.o:		common/zorac
tsort.o:		common/prnprn
tsort.o:		common/t_files
tsort.o:		common/t_flname
tsort.o:		common/t_iobuff
tsort.o:		common/t_pointr
tsort.o:		common/t_sect
tsort.o:		common/titan
tsort.o:		common/utilc
util1.o:		../drf/comdrf/clas
util1.o:		../drf/comdrf/drfamb
util1.o:		../drf/comdrf/drfbem
util1.o:		../drf/comdrf/drfin
util1.o:		../drf/comdrf/drfpar
util1.o:		../drf/comdrf/extinf
util1.o:		../drf/comdrf/mcener
util1.o:		../drf/comdrf/mcinp
util1.o:		../drf/comdrf/rfene
util1.o:		../drf/comdrf/rota
util1.o:		../drf/comdrf/runpar
util1.o:		../drf/comdrf/sizescon
util1.o:		../drf/comdrf/sizesrf
util1.o:		common/atmblk
util1.o:		common/chgcc
util1.o:		common/cigrad
util1.o:		common/cndx40
util1.o:		common/cndx41
util1.o:		common/common
util1.o:		common/connolly
util1.o:		common/coval
util1.o:		common/csubch
util1.o:		common/czmat
util1.o:		common/disc
util1.o:		common/discc
util1.o:		common/dm
util1.o:		common/drfopt
util1.o:		common/dump3
util1.o:		common/errcodes
util1.o:		common/gmempara
util1.o:		common/harmon
util1.o:		common/infoa
util1.o:		common/infob
util1.o:		common/iofile
util1.o:		common/junkcc
util1.o:		common/machin
util1.o:		common/mapper
util1.o:		common/maxlen
util1.o:		common/nshel
util1.o:		common/overrule
util1.o:		common/parcntl
util1.o:		common/phycon
util1.o:		common/prints
util1.o:		common/repp
util1.o:		common/restar
util1.o:		common/restri
util1.o:		common/restrj
util1.o:		common/restrl
util1.o:		common/runlab
util1.o:		common/scfopt
util1.o:		common/scfwfn
util1.o:		common/scra7
util1.o:		common/sector
util1.o:		common/segm
util1.o:		common/sizes
util1.o:		common/sizes
util1.o:		common/structchk
util1.o:		common/timeperiods
util1.o:		common/timez
util1.o:		common/tran
util1.o:		common/utilc
util1.o:		common/work
util1.o:		common/workc
util1.o:		common/zorac
util2.o:		common/atmblk
util2.o:		common/ccpdft.hf77
util2.o:		common/cslosc
util2.o:		common/disc
util2.o:		common/discc
util2.o:		common/dm
util2.o:		common/dmisc
util2.o:		common/filel
util2.o:		common/files
util2.o:		common/flips
util2.o:		common/grad2
util2.o:		common/ijlab
util2.o:		common/incrd
util2.o:		common/indez
util2.o:		common/infoa
util2.o:		common/iofile
util2.o:		common/junkcc
util2.o:		common/mapper
util2.o:		common/maxlen
util2.o:		common/misc
util2.o:		common/nshel
util2.o:		common/pkfil
util2.o:		common/restar
util2.o:		common/shlnos
util2.o:		common/shlt
util2.o:		common/sizes
util2.o:		common/sortp
util3.o:		common/atmblk
util3.o:		common/ccpdft.hf77
util3.o:		common/chmgms
util3.o:		common/cigrad
util3.o:		common/cndx41
util3.o:		common/common
util3.o:		common/czmat
util3.o:		common/discc
util3.o:		common/dump3
util3.o:		common/files
util3.o:		common/funct
util3.o:		common/infoa
util3.o:		common/infoa
util3.o:		common/iofile
util3.o:		common/iofile
util3.o:		common/mapper
util3.o:		common/nshel
util3.o:		common/parallel
util3.o:		common/parcntl
util3.o:		common/runlab
util3.o:		common/sizes
util3.o:		common/statis
util3.o:		common/timeperiods
util3.o:		common/timez
util3.o:		common/work
util4.o:		common/iofile
util4.o:		common/jobopt
util4.o:		common/mcff
util4.o:		common/multic
util4.o:		common/prints
util4.o:		common/prnprn
util4.o:		common/sizes
util4.o:		common/syminf
util5.o:		common/atmblk
util5.o:		common/blksiz
util5.o:		common/cndx40
util5.o:		common/disktl
util5.o:		common/helpr
util5.o:		common/iofile
util5.o:		common/mapper
util5.o:		common/mlngth
util5.o:		common/prnprn
util5.o:		common/sizes
util5.o:		common/stak
util5.o:		common/three
util6.o:		common/atmblk
util6.o:		common/cigrad
util6.o:		common/cndx41
util6.o:		common/common
util6.o:		common/ghfblk
util6.o:		common/infoa
util6.o:		common/iofile
util6.o:		common/prnprn
util6.o:		common/sizes
util7.o:		common/atmblk
util7.o:		common/blksiz
util7.o:		common/cndx41
util7.o:		common/common
util7.o:		common/gmempara
util7.o:		common/iofile
util7.o:		common/mapper
util7.o:		common/mpcut
util7.o:		common/sizes
util7.o:		common/stak
util7.o:		common/three
util7.o:		common/uhfspn
util8.o:		common/sizes
zora.o:		common/atmblk
zora.o:		common/atmol3
zora.o:		common/common
zora.o:		common/cslosc
zora.o:		common/czmat
zora.o:		common/datgue
zora.o:		common/dump3
zora.o:		common/funct
zora.o:		common/harmon
zora.o:		common/infoa
zora.o:		common/infoa
zora.o:		common/infob
zora.o:		common/iofile
zora.o:		common/machin
zora.o:		common/mapper
zora.o:		common/maxlen
zora.o:		common/nshel
zora.o:		common/occnr_z
zora.o:		common/restar
zora.o:		common/restri
zora.o:		common/restrl
zora.o:		common/runlab
zora.o:		common/scfopt
zora.o:		common/scfwfn
zora.o:		common/scra7
zora.o:		common/sizes
zora.o:		common/sortp
zora.o:		common/symtry
zora.o:		common/tran
zora.o:		common/transf
zora.o:		common/unocas
zora.o:		common/zorac
#
# ============================ NON-COMPILED TARGETS ===========================
#
clean:
	touch x.c x.o x~ x.f \#x
	rm -f *.c *.o *.f *.mod *~ \#* $(EXE) *.exe *.lib
	rm -f gms_version.m
	rm -f notes.stamp
	rm -f ../libs/TABLE ../libs/ecplib
	@(if test -d ../dft; then \
	cd ../dft; make clean; fi)
	@(if test -d ../drf; then \
	cd ../drf; make clean; fi)
	@(if test -d ../linalg; then \
	cd ../linalg; make clean; fi)
	@(if test -d ../dl-find; then \
	cd ../dl-find; $(GMAKE) clean; fi)
	@(if test -d ../vb; then \
	cd ../vb; make clean; fi)
	@(cd ../utilities; make clean)
#
realclean:
	make clean






depend:
	../utilities/makedeps
	@(if test -d ../dft; then \
	cd ../dft; make depend; fi)
	@(if test -d ../drf; then \
	cd ../drf; make depend; fi)
	@(if test -d ../mopac; then \
	cd ../mopac; make depend; fi)
	@(if test -d ../newscf; then \
	cd ../newscf; make depend; fi)
	@(if test -d ../nmr; then \
	cd ../nmr; make depend; fi)
	@(if test -d ../vb; then \
	cd ../vb; make depend; fi)
#
qa:
	@cat ../utilities/gener.m ../qa/testdirs.txt.in | $(M4) $(M4OPTS) > ../qa/testdirs.txt
keys:
	@echo KEYS GEN_OPTIONS M4_OPTIONS | ${M4} $(M4OPTS)|  sed 's/,/ /g'
#
scfsrc:	
	@make ` echo $(SCFOBJ) | sed 's/\.o/\.f/g' `
	@make c.c
#
src:
	@$(MAKE) ` echo $(OBJF) | sed 's/\.o/\.f/g' `
	@$(MAKE) ` echo $(OBJC) | sed 's/\.o/\.c/g' `
	@$(MAKE) makefile_f
#
makefile_f:
	cat ../utilities/gener.m  makefile.m | $(M4) $(M4OPTS) > makefile_f
#
#
#

gms_version.m:  FORCE
		which svnversion > ./m4_svnversion; chmod +x ./m4_svnversion
		echo "      subroutine gms_version(cdate,ctime,cname,cversion)" > gms_version.m
		echo "      character*10 cdate,cname,cversion" >> gms_version.m
		echo "      character*5 ctime" >> gms_version.m
		echo "      cdate=\""`date '+%d-%m-%Y'`"\"" >> gms_version.m
		echo "      ctime=\""`date '+%H:%M'`"\"" >> gms_version.m
		echo "      cname=\""`${WHOAMI}`"\"" >> gms_version.m
		echo "      cversion=\""`./m4_svnversion ..`"\"" >> gms_version.m
		echo "      return" >> gms_version.m
		echo "      end" >> gms_version.m
		echo "      subroutine getm4keys(m4keys)" >> gms_version.m
		echo "      implicit none" >> gms_version.m
		echo "      character*(*) m4keys" >> gms_version.m
		make keys | grep KEYS | sed 's/KEYS//' | fmt -w50 |awk 'BEGIN {i=0} {if (i==0) {printf "       m4keys=\"%s\"\n",$$0} else {printf "     &//\"%s \"\n", $$0} i++}' >> gms_version.m
		echo "      return" >> gms_version.m
		echo "      end" >> gms_version.m
		echo "      subroutine wrtkeys(iwr)" >> gms_version.m
		make keys | grep KEYS | sed 's/KEYS//' | fmt -w50 | awk '{printf "      write(iwr,*)\n     &\"M4-Keys: %s \"\n", $$0}'  >> gms_version.m
		echo "      return" >> gms_version.m
		echo "      end" >> gms_version.m

#
config:
	echo there is no config target - Use ./configure
