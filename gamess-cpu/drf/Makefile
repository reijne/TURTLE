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
all:	
	@(cd ../m4;make)

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
FFLAGSD = -c -pg ${FFLAGSI8}
CFLAGSD= -c -pg ${CFLAGSI8}

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
EXTRA=gethes.o
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
dppop.o:		../m4/common/chgcc
dppop.o:		../m4/common/drfopt
dppop.o:		../m4/common/dump3
dppop.o:		../m4/common/g80nb
dppop.o:		../m4/common/infoa
dppop.o:		../m4/common/iofile
dppop.o:		../m4/common/limy
dppop.o:		../m4/common/modj
dppop.o:		../m4/common/nshel
dppop.o:		../m4/common/prints
dppop.o:		../m4/common/prnprn
dppop.o:		../m4/common/restar
dppop.o:		../m4/common/restri
dppop.o:		../m4/common/root
dppop.o:		../m4/common/runlab
dppop.o:		../m4/common/segm
dppop.o:		../m4/common/sizes
dppop.o:		../m4/common/statis
dppop.o:		../m4/common/symtry
dppop.o:		../m4/common/timez
dppop.o:		../m4/common/xfield
dppop.o:		comdrf/auxdrf
dppop.o:		comdrf/bas
dppop.o:		comdrf/cenmas
dppop.o:		comdrf/dafil
dppop.o:		comdrf/darw
dppop.o:		comdrf/drfamb
dppop.o:		comdrf/drfdaf
dppop.o:		comdrf/elpinf
dppop.o:		comdrf/expan
dppop.o:		comdrf/ihelp
dppop.o:		comdrf/ijpair
dppop.o:		comdrf/iofil
dppop.o:		comdrf/mollab
dppop.o:		comdrf/nmorb
dppop.o:		comdrf/opt
dppop.o:		comdrf/qpanl
dppop.o:		comdrf/runpar
dppop.o:		comdrf/rys
dppop.o:		comdrf/scfopt
dppop.o:		comdrf/scm
dppop.o:		comdrf/sizesrf
dppop.o:		comdrf/ssgg
dppop.o:		comdrf/stv
dppop.o:		comdrf/sym
dppop.o:		comdrf/tim
drfctl.o:		../m4/common/drfopt
drfctl.o:		../m4/common/infoa
drfctl.o:		../m4/common/sizes
drfctl.o:		comdrf/clas
drfctl.o:		comdrf/dafil
drfctl.o:		comdrf/drfbem
drfctl.o:		comdrf/drfdaf
drfctl.o:		comdrf/drfpar
drfctl.o:		comdrf/extinf
drfctl.o:		comdrf/free
drfctl.o:		comdrf/grinf
drfctl.o:		comdrf/ijpair
drfctl.o:		comdrf/inddaf
drfctl.o:		comdrf/iofil
drfctl.o:		comdrf/mcener
drfctl.o:		comdrf/mcmod
drfctl.o:		comdrf/mollab
drfctl.o:		comdrf/neqpar
drfctl.o:		comdrf/neqtex
drfctl.o:		comdrf/rfene
drfctl.o:		comdrf/runpar
drfctl.o:		comdrf/sizesrf
drfder.o:		../m4/common/cndx41
drfder.o:		../m4/common/common
drfder.o:		../m4/common/hermit
drfder.o:		../m4/common/infoa
drfder.o:		../m4/common/inxblk
drfder.o:		../m4/common/iofile
drfder.o:		../m4/common/mapper
drfder.o:		../m4/common/misc
drfder.o:		../m4/common/nodinf
drfder.o:		../m4/common/nshel
drfder.o:		../m4/common/prnprn
drfder.o:		../m4/common/restar
drfder.o:		../m4/common/runlab
drfder.o:		../m4/common/segm
drfder.o:		../m4/common/sizes
drfder.o:		../m4/common/symtry
drfder.o:		../m4/common/timez
drfder.o:		../m4/common/wermit
drfextra.o:		../m4/common/cslosc
drfextra.o:		../m4/common/infoa
drfextra.o:		../m4/common/nshel
drfextra.o:		../m4/common/sizes
drfextra.o:		comdrf/alpbet
drfextra.o:		comdrf/bas
drfextra.o:		comdrf/dafil
drfextra.o:		comdrf/darw
drfextra.o:		comdrf/drfpar
drfextra.o:		comdrf/drfzfa
drfextra.o:		comdrf/ijpair
drfextra.o:		comdrf/iofil
drfextra.o:		comdrf/mem
drfextra.o:		comdrf/mollab
drfextra.o:		comdrf/opt
drfextra.o:		comdrf/runpar
drfextra.o:		comdrf/scfopt
drfextra.o:		comdrf/sizesrf
drfgrad.o:		../m4/common/drfopt
drfgrad.o:		../m4/common/hermit
drfgrad.o:		../m4/common/infoa
drfgrad.o:		../m4/common/sizes
drfgrad.o:		../m4/common/wermit
drfgrad.o:		comdrf/alpbet
drfgrad.o:		comdrf/auxdrf
drfgrad.o:		comdrf/commonrf
drfgrad.o:		comdrf/dafil
drfgrad.o:		comdrf/drfamb
drfgrad.o:		comdrf/drfbem
drfgrad.o:		comdrf/drfdaf
drfgrad.o:		comdrf/drfexp
drfgrad.o:		comdrf/drfind
drfgrad.o:		comdrf/drfint
drfgrad.o:		comdrf/drfpar
drfgrad.o:		comdrf/drfzfa
drfgrad.o:		comdrf/extinf
drfgrad.o:		comdrf/grinf
drfgrad.o:		comdrf/ihelp
drfgrad.o:		comdrf/ijpair
drfgrad.o:		comdrf/iofil
drfgrad.o:		comdrf/mollab
drfgrad.o:		comdrf/neqpar
drfgrad.o:		comdrf/rad
drfgrad.o:		comdrf/rfgrad
drfgrad.o:		comdrf/runpar
drfgrad.o:		comdrf/scfopt
drfgrad.o:		comdrf/scfpar
drfgrad.o:		comdrf/sizesrf
drfscf.o:		../m4/common/infoa
drfscf.o:		../m4/common/sizes
drfscf.o:		comdrf/drfbem
drfscf.o:		comdrf/drfdaf
drfscf.o:		comdrf/drfexp
drfscf.o:		comdrf/drfind
drfscf.o:		comdrf/drfpar
drfscf.o:		comdrf/ijpair
drfscf.o:		comdrf/iofil
drfscf.o:		comdrf/opt
drfscf.o:		comdrf/shlint
drfscf.o:		comdrf/sizesrf
drfsolv.o:		../m4/common/infoa
drfsolv.o:		../m4/common/sizes
drfsolv.o:		comdrf/auxdrf
drfsolv.o:		comdrf/drfamb
drfsolv.o:		comdrf/drfbem
drfsolv.o:		comdrf/drfdaf
drfsolv.o:		comdrf/drfexp
drfsolv.o:		comdrf/drfpar
drfsolv.o:		comdrf/extinf
drfsolv.o:		comdrf/grinf
drfsolv.o:		comdrf/ijpair
drfsolv.o:		comdrf/iofil
drfsolv.o:		comdrf/mollab
drfsolv.o:		comdrf/neqpar
drfsolv.o:		comdrf/rad
drfsolv.o:		comdrf/runpar
drfsolv.o:		comdrf/sizesrf
drfsub.o:		../m4/common/connolly
drfsub.o:		../m4/common/drfopt
drfsub.o:		../m4/common/dump3
drfsub.o:		../m4/common/gvalue
drfsub.o:		../m4/common/infoa
drfsub.o:		../m4/common/iofile
drfsub.o:		../m4/common/nshel
drfsub.o:		../m4/common/runlab
drfsub.o:		../m4/common/scfopt
drfsub.o:		../m4/common/sector
drfsub.o:		../m4/common/sizes
drfsub.o:		../m4/common/symtry
drfsub.o:		comdrf/alpbet
drfsub.o:		comdrf/assign
drfsub.o:		comdrf/auxdrf
drfsub.o:		comdrf/bas
drfsub.o:		comdrf/ciopt
drfsub.o:		comdrf/clas
drfsub.o:		comdrf/dafil
drfsub.o:		comdrf/darw
drfsub.o:		comdrf/defpar
drfsub.o:		comdrf/drfamb
drfsub.o:		comdrf/drfbem
drfsub.o:		comdrf/drfdaf
drfsub.o:		comdrf/drfexp
drfsub.o:		comdrf/drfin
drfsub.o:		comdrf/drfind
drfsub.o:		comdrf/drfint
drfsub.o:		comdrf/drfpar
drfsub.o:		comdrf/drfzfa
drfsub.o:		comdrf/ecp
drfsub.o:		comdrf/elpinf
drfsub.o:		comdrf/ene
drfsub.o:		comdrf/enrgci
drfsub.o:		comdrf/extinf
drfsub.o:		comdrf/free
drfsub.o:		comdrf/grinf
drfsub.o:		comdrf/ihelp
drfsub.o:		comdrf/ijpair
drfsub.o:		comdrf/inddaf
drfsub.o:		comdrf/iofil
drfsub.o:		comdrf/max
drfsub.o:		comdrf/mcinf
drfsub.o:		comdrf/mollab
drfsub.o:		comdrf/neqpar
drfsub.o:		comdrf/neqtex
drfsub.o:		comdrf/opt
drfsub.o:		comdrf/rad
drfsub.o:		comdrf/rfene
drfsub.o:		comdrf/runpar
drfsub.o:		comdrf/rys
drfsub.o:		comdrf/scfopt
drfsub.o:		comdrf/scfpar
drfsub.o:		comdrf/scm
drfsub.o:		comdrf/shlint
drfsub.o:		comdrf/sizescon
drfsub.o:		comdrf/sizesrf
drfsub.o:		comdrf/ssgg
drfsub.o:		comdrf/stv
drfsub.o:		comdrf/sym
drfsub.o:		comdrf/tim
extgrad.o:		../m4/common/drfopt
extgrad.o:		../m4/common/infoa
extgrad.o:		../m4/common/iofile
extgrad.o:		../m4/common/sizes
extgrad.o:		comdrf/alpbet
extgrad.o:		comdrf/bas
extgrad.o:		comdrf/commonrf
extgrad.o:		comdrf/drfamb
extgrad.o:		comdrf/drfbem
extgrad.o:		comdrf/drfexp
extgrad.o:		comdrf/drfind
extgrad.o:		comdrf/drfpar
extgrad.o:		comdrf/extinf
extgrad.o:		comdrf/grinf
extgrad.o:		comdrf/ihelp
extgrad.o:		comdrf/ijpair
extgrad.o:		comdrf/iofil
extgrad.o:		comdrf/mollab
extgrad.o:		comdrf/neqpar
extgrad.o:		comdrf/opt
extgrad.o:		comdrf/rad
extgrad.o:		comdrf/rfgrad
extgrad.o:		comdrf/scfopt
extgrad.o:		comdrf/scfpar
extgrad.o:		comdrf/shlint
extgrad.o:		comdrf/sizesrf
fld.o:		../m4/common/infoa
fld.o:		../m4/common/nshel
fld.o:		../m4/common/sizes
fld.o:		comdrf/bas
fld.o:		comdrf/dafil
fld.o:		comdrf/ijpair
fld.o:		comdrf/mem
fld.o:		comdrf/mollab
fld.o:		comdrf/opt
fld.o:		comdrf/scfopt
fld.o:		comdrf/sizesrf
hondext.o:		../m4/common/infoa
hondext.o:		../m4/common/sizes
hondext.o:		comdrf/alpbet
hondext.o:		comdrf/cenmas
hondext.o:		comdrf/dafil
hondext.o:		comdrf/drfpar
hondext.o:		comdrf/expan
hondext.o:		comdrf/extinf
hondext.o:		comdrf/ijpair
hondext.o:		comdrf/iofil
hondext.o:		comdrf/mollab
hondext.o:		comdrf/nmorb
hondext.o:		comdrf/opt
hondext.o:		comdrf/runpar
hondext.o:		comdrf/scfopt
hondext.o:		comdrf/sizesrf
mcsub.o:		../m4/common/drfopt
mcsub.o:		../m4/common/infoa
mcsub.o:		../m4/common/restar
mcsub.o:		../m4/common/sizes
mcsub.o:		comdrf/auxdrf
mcsub.o:		comdrf/clas
mcsub.o:		comdrf/dafil
mcsub.o:		comdrf/drfamb
mcsub.o:		comdrf/drfbem
mcsub.o:		comdrf/drfdaf
mcsub.o:		comdrf/drfpar
mcsub.o:		comdrf/drfzfa
mcsub.o:		comdrf/enrgci
mcsub.o:		comdrf/extinf
mcsub.o:		comdrf/grinf
mcsub.o:		comdrf/ijpair
mcsub.o:		comdrf/inddaf
mcsub.o:		comdrf/iofil
mcsub.o:		comdrf/mcener
mcsub.o:		comdrf/mcgeom
mcsub.o:		comdrf/mcinf
mcsub.o:		comdrf/mcinp
mcsub.o:		comdrf/mcmod
mcsub.o:		comdrf/mollab
mcsub.o:		comdrf/neqpar
mcsub.o:		comdrf/neqtex
mcsub.o:		comdrf/opt
mcsub.o:		comdrf/rad
mcsub.o:		comdrf/rfene
mcsub.o:		comdrf/runpar
mcsub.o:		comdrf/scfopt
mcsub.o:		comdrf/sizesrf
mscon.o:		../m4/common/connolly
mscon.o:		../m4/common/sizes
mscon.o:		comdrf/drfbem
mscon.o:		comdrf/iofil
mscon.o:		comdrf/sizescon
neqpot.o:		../m4/common/infoa
neqpot.o:		../m4/common/sizes
neqpot.o:		comdrf/auxdrf
neqpot.o:		comdrf/dafil
neqpot.o:		comdrf/drfdaf
neqpot.o:		comdrf/drfpar
neqpot.o:		comdrf/drfzfa
neqpot.o:		comdrf/ijpair
neqpot.o:		comdrf/inddaf
neqpot.o:		comdrf/iofil
neqpot.o:		comdrf/mollab
neqpot.o:		comdrf/neqpar
neqpot.o:		comdrf/neqtex
neqpot.o:		comdrf/rfene
neqpot.o:		comdrf/runpar
neqpot.o:		comdrf/sizesahnqcost
neqpot.o:		comdrf/sizesrf
radius.o:		../m4/common/drfopt
radius.o:		../m4/common/infoa
radius.o:		../m4/common/sizes
radius.o:		comdrf/drfamb
radius.o:		comdrf/drfbem
radius.o:		comdrf/drfdaf
radius.o:		comdrf/drfpar
radius.o:		comdrf/drfzfa
radius.o:		comdrf/extinf
radius.o:		comdrf/grinf
radius.o:		comdrf/iofil
radius.o:		comdrf/mcener
radius.o:		comdrf/mollab
radius.o:		comdrf/rad
radius.o:		comdrf/rfene
radius.o:		comdrf/runpar
radius.o:		comdrf/sizesrf
rfin.o:		../m4/common/connolly
rfin.o:		../m4/common/infob
rfin.o:		../m4/common/sizes
rfin.o:		comdrf/assign
rfin.o:		comdrf/ciopt
rfin.o:		comdrf/commonrf
rfin.o:		comdrf/drfin
rfin.o:		comdrf/scfopt
rfin.o:		comdrf/sizescon
rfin.o:		comdrf/sizesrf
rfin.o:		comdrf/trnsf
surfsub.o:		../m4/common/connolly
surfsub.o:		../m4/common/infoa
surfsub.o:		../m4/common/sizes
surfsub.o:		comdrf/dafil
surfsub.o:		comdrf/drfbem
surfsub.o:		comdrf/drfdaf
surfsub.o:		comdrf/drfin
surfsub.o:		comdrf/drfpar
surfsub.o:		comdrf/extinf
surfsub.o:		comdrf/iofil
surfsub.o:		comdrf/max
surfsub.o:		comdrf/mcinf
surfsub.o:		comdrf/mem
surfsub.o:		comdrf/rad
surfsub.o:		comdrf/runpar
surfsub.o:		comdrf/sizesahsurfsub
surfsub.o:		comdrf/sizescon
surfsub.o:		comdrf/sizesrf
#
# ============================ NON-COMPILED TARGETS ===========================
#
#






depend:
	../utilities/makedeps
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
#
# This is an extra section designed to be
# appended to a copy the GAMESS-UK Makefile
# 
# Use the configure script in ../m4
#
#
DRFOBJ = drfctl.o \
	drfsub.o \
	neqpot.o \
	radius.o \
	drfsolv.o \
	drfscf.o \
	surfsub.o \
	hondext.o \
	dppop.o \
	fld.o \
	missing.o \
	mcsub.o \
	drfextra.o \
	rfin.o \
	mscon.o \
	drfgrad.o \
	extgrad.o \
	drfder.o 

DRFLIB = libdrf.a

drflib:	$(DRFLIB)

$(DRFLIB):	$(DRFOBJ)
		$(ARCHIVE) $(DRFLIB)  $(DRFOBJ) 
		$(RANLIB) $(DRFLIB)

clean:
	touch x.o x~ x.f \#x
	rm -f *.o *.f *~ \#* $(DRFLIB) 

.DEFAULT:
	(cd ../m4; make $@)

#
