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
chf.o:		../m4/common/mapper
chf.o:		../m4/common/sizes
chf.o:		common/dft_dfder
chf.o:		common/dft_parameters
coulomb.o:		../m4/common/errcodes
coulomb.o:		../m4/common/gmempara
coulomb.o:		../m4/common/iofile
coulomb.o:		../m4/common/restar
coulomb.o:		../m4/common/symtry
coulomb.o:		../m4/common/timeperiods
coulomb.o:		common/dft_api
coulomb.o:		common/dft_basis
coulomb.o:		common/dft_basis
coulomb.o:		common/dft_basis_api
coulomb.o:		common/dft_basis_cont_inf
coulomb.o:		common/dft_dunlap
coulomb.o:		common/dft_dunlap_ga
coulomb.o:		common/dft_mbasis
coulomb.o:		common/dft_mbasis
coulomb.o:		common/dft_memory_info
coulomb.o:		common/dft_module_comm
coulomb.o:		common/dft_mol_info
coulomb.o:		common/dft_parameters
deriv2e.o:		../m4/common/cigrad
deriv2e.o:		../m4/common/cslosc
deriv2e.o:		../m4/common/ghfblk
deriv2e.o:		../m4/common/ijlab
deriv2e.o:		../m4/common/infoa
deriv2e.o:		../m4/common/parallel
deriv2e.o:		../m4/common/prints
deriv2e.o:		../m4/common/restar
deriv2e.o:		../m4/common/restrl
deriv2e.o:		../m4/common/segm
deriv2e.o:		../m4/common/sizes
deriv2e.o:		../m4/common/sizes
deriv2e.o:		../m4/common/specal
deriv2e.o:		../m4/common/timez
deriv2e.o:		common/ccpdft.hf77
deriv2e.o:		common/dft_d2escr
deriv2e.o:		common/dft_dfgfac
deriv2e.o:		common/dft_dmisc
deriv2e.o:		common/dft_dshlnf
deriv2e.o:		common/dft_dshlno
deriv2e.o:		common/dft_grad2
deriv2e.o:		common/dft_incrd
deriv2e.o:		common/dft_intctl
deriv2e.o:		common/dft_iofile
deriv2e.o:		common/dft_mapper
deriv2e.o:		common/dft_mbasis
deriv2e.o:		common/dft_module_comm
deriv2e.o:		common/dft_module_comm
deriv2e.o:		common/dft_mol_info
deriv2e.o:		common/dft_parameters
deriv2e.o:		common/dft_picon
deriv2e.o:		common/dft_root
dftbasis.o:		../m4/common/sizes
dftbasis.o:		common/dft_basis
dftbasis.o:		common/dft_basis
dftbasis.o:		common/dft_basis_api
dftbasis.o:		common/dft_basis_api
dftbasis.o:		common/dft_mbasis
dftbasis.o:		common/dft_memory_info
dftbasis.o:		common/dft_module_comm
dftbasis.o:		common/dft_mol_info
dftbasis.o:		common/dft_order_info
dftbasis.o:		common/dft_parameters
dftbasis.o:		common/dft_parameters
exp_dksm_hess.o:		../m4/common/mapper
exp_dksm_hess.o:		../m4/common/sizes
exp_dksm_hess.o:		common/dft_basder
exp_dksm_hess.o:		common/dft_basis
exp_dksm_hess.o:		common/dft_dfder
exp_dksm_hess.o:		common/dft_parameters
gadft.o:		../m4/common/gmempara
gadft.o:		../m4/common/iofile
gadft.o:		../m4/common/parcntl
gadft.o:		../m4/common/vcore
gamess.o:		../m4/common/gmempara
gamess.o:		../m4/common/infoa
gamess.o:		../m4/common/infob
gamess.o:		../m4/common/nshel
gamess.o:		../m4/common/runlab
gamess.o:		../m4/common/sizes
gamess.o:		common/dft_basis
gamess.o:		common/dft_basis_api
gamess.o:		common/dft_iofile
gamess.o:		common/dft_module_comm
gamess.o:		common/dft_mol_info
gamess.o:		common/dft_parameters
gden.o:		../m4/common/blur
gden.o:		../m4/common/iofile
gden.o:		../m4/common/iofile
gden.o:		../m4/common/prints
gden.o:		../m4/common/restar
gden.o:		../m4/common/runlab
gden.o:		../m4/common/sizes
gden.o:		../m4/common/symtry
gden.o:		common/ccpdft.hf77
gden.o:		common/dft_api
gden.o:		common/dft_basis_api
gden.o:		common/dft_intctl
gden.o:		common/dft_iofile
gden.o:		common/dft_memory_info
gden.o:		common/dft_module_comm
gden.o:		common/dft_mol_info
gden.o:		common/dft_parameters
gden.o:		common/dft_parameters
gden.o:		common/mden
geometry.o:		common/dft_module_comm
geometry.o:		common/dft_mol_info
geometry.o:		common/dft_parameters
global.o:		../m4/common/statis
global.o:		common/dft_api
global.o:		common/dft_basis
global.o:		common/dft_basis_api
global.o:		common/dft_module_comm
global.o:		common/dft_mol_info
global.o:		common/dft_parameters
global.o:		common/dft_parameters
integ2e.o:		../m4/common/cslosc
integ2e.o:		../m4/common/infoa
integ2e.o:		../m4/common/parallel
integ2e.o:		../m4/common/parcntl
integ2e.o:		../m4/common/restar
integ2e.o:		../m4/common/sizes
integ2e.o:		../m4/common/statis
integ2e.o:		../m4/common/timez
integ2e.o:		common/dft_auxg
integ2e.o:		common/dft_basis
integ2e.o:		common/dft_basis_api
integ2e.o:		common/dft_d2escr
integ2e.o:		common/dft_denss
integ2e.o:		common/dft_dmisc
integ2e.o:		common/dft_dshlnf
integ2e.o:		common/dft_flips
integ2e.o:		common/dft_ijlab
integ2e.o:		common/dft_incrd
integ2e.o:		common/dft_indez
integ2e.o:		common/dft_intctl
integ2e.o:		common/dft_inxblk
integ2e.o:		common/dft_iofile
integ2e.o:		common/dft_mapper
integ2e.o:		common/dft_mbasis
integ2e.o:		common/dft_misc
integ2e.o:		common/dft_module_comm
integ2e.o:		common/dft_parameters
integ2e.o:		common/dft_picon
integ2e.o:		common/dft_root
integ2e.o:		common/dft_rtdata
integ2e.o:		common/dft_setint
integ2e.o:		common/dft_shlinf
integ2e.o:		common/dft_shlnos
integ2e.o:		common/dft_shlt
integ_data.o:		common/dft_inxblk
integ_data.o:		common/dft_rtdata
integ_te2c_norm.o:		../m4/common/cslosc
integ_te2c_norm.o:		../m4/common/infoa
integ_te2c_norm.o:		../m4/common/parallel
integ_te2c_norm.o:		../m4/common/parcntl
integ_te2c_norm.o:		../m4/common/restar
integ_te2c_norm.o:		../m4/common/sizes
integ_te2c_norm.o:		../m4/common/statis
integ_te2c_norm.o:		../m4/common/timez
integ_te2c_norm.o:		common/dft_flips
integ_te2c_norm.o:		common/dft_ijlab
integ_te2c_norm.o:		common/dft_indez
integ_te2c_norm.o:		common/dft_iofile
integ_te2c_norm.o:		common/dft_mapper
integ_te2c_norm.o:		common/dft_mbasis
integ_te2c_norm.o:		common/dft_misc
integ_te2c_norm.o:		common/dft_picon
integ_te2c_norm.o:		common/dft_root
integ_te2c_norm.o:		common/dft_shlg70
integ_te2c_norm.o:		common/dft_shlnos
integ_te2c_norm.o:		common/dft_shlt
integ_te3c_rep.o:		../m4/common/cslosc
integ_te3c_rep.o:		../m4/common/files
integ_te3c_rep.o:		../m4/common/infoa
integ_te3c_rep.o:		../m4/common/nshel
integ_te3c_rep.o:		../m4/common/parallel
integ_te3c_rep.o:		../m4/common/parcntl
integ_te3c_rep.o:		../m4/common/prints
integ_te3c_rep.o:		../m4/common/restar
integ_te3c_rep.o:		../m4/common/sizes
integ_te3c_rep.o:		../m4/common/statis
integ_te3c_rep.o:		../m4/common/timez
integ_te3c_rep.o:		common/dft_auxg
integ_te3c_rep.o:		common/dft_bshell
integ_te3c_rep.o:		common/dft_denss
integ_te3c_rep.o:		common/dft_flips
integ_te3c_rep.o:		common/dft_ijlab
integ_te3c_rep.o:		common/dft_indez
integ_te3c_rep.o:		common/dft_intctl
integ_te3c_rep.o:		common/dft_inxblk
integ_te3c_rep.o:		common/dft_iofile
integ_te3c_rep.o:		common/dft_mapper
integ_te3c_rep.o:		common/dft_mbasis
integ_te3c_rep.o:		common/dft_misc
integ_te3c_rep.o:		common/dft_picon
integ_te3c_rep.o:		common/dft_root
integ_te3c_rep.o:		common/dft_setint
integ_te3c_rep.o:		common/dft_shlg70
integ_te3c_rep.o:		common/dft_shlinf
integ_te3c_rep.o:		common/dft_shlnos
integ_te3c_rep.o:		common/dft_shlt
interface.o:		../m4/common/errcodes
interface.o:		../m4/common/nshel
interface.o:		../m4/common/sizes
interface.o:		common/basato.hf77
interface.o:		common/ccpdft.hf77
interface.o:		common/dft_api
interface.o:		common/dft_basis
interface.o:		common/dft_basis_api
interface.o:		common/dft_dunlap
interface.o:		common/dft_intctl
interface.o:		common/dft_iofile
interface.o:		common/dft_mbasis
interface.o:		common/dft_module_comm
interface.o:		common/dft_module_comm
interface.o:		common/dft_mol_info
interface.o:		common/dft_parameters
interface.o:		common/dft_xc
interface.o:		common/dft_xc
interface.o:		common/test4
intpack.o:		common/dft_basis
intpack.o:		common/dft_memory_info
intpack.o:		common/dft_module_comm
intpack.o:		common/dft_mol_info
intpack.o:		common/dft_numbers
intpack.o:		common/dft_order_info
intpack.o:		common/dft_parameters
intpack.o:		common/dft_physical_constants
jfitg.o:		../m4/common/errcodes
jfitg.o:		../m4/common/iofile
jfitg.o:		../m4/common/symtry
jfitg.o:		../m4/common/timeperiods
jfitg.o:		common/dft_api
jfitg.o:		common/dft_basis
jfitg.o:		common/dft_basis_api
jfitg.o:		common/dft_basis_cont_inf
jfitg.o:		common/dft_dunlap
jfitg.o:		common/dft_dunlap_ga
jfitg.o:		common/dft_mbasis
jfitg.o:		common/dft_memory_info
jfitg.o:		common/dft_module_comm
jfitg.o:		common/dft_mol_info
jfitg.o:		common/dft_parameters
matform.o:		../m4/common/nshel
matform.o:		../m4/common/sizes
matform.o:		common/dft_dunlap_ga
matform.o:		common/dft_indez
matform.o:		common/dft_mbasis
matform.o:		common/dft_shlnos
matpack.o:		common/dft_module_comm
matpack.o:		common/dft_parameters
mpole1.o:		common/basato.hf77
mpole1.o:		common/crys_params.hf77
mpole1.o:		common/crys_params.hf77
mpole2.o:		../m4/common/nshel
mpole2.o:		../m4/common/sizes
mpole2.o:		../m4/common/timeperiods
mpole2.o:		common/basato.hf77
mpole2.o:		common/claddd_crys.hf77
mpole2.o:		common/dfacomm_crys.hf77
mpole2.o:		common/dft_basis
mpole2.o:		common/dft_module_comm
mpole2.o:		common/dft_mol_info
mpole2.o:		common/dft_parameters
mpole2.o:		common/parinf.hf77
mpole2.o:		common/servi_inpbas.hf77
mpole2.o:		common/tcomm_crys.hf77
readbasis.o:		../m4/common/datbas
readbasis.o:		../m4/common/infoa
readbasis.o:		../m4/common/infoa
readbasis.o:		../m4/common/infob
readbasis.o:		../m4/common/iofile
readbasis.o:		../m4/common/iofile
readbasis.o:		../m4/common/runlab
readbasis.o:		../m4/common/runlab
readbasis.o:		../m4/common/sizes
readbasis.o:		../m4/common/work
readbasis.o:		common/readbasis
weights.o:		../m4/common/mapper
weights.o:		../m4/common/sizes
weights.o:		common/dft_basder
weights.o:		common/dft_dij
weights.o:		common/dft_module_comm
weights.o:		common/dft_mol_info
weights.o:		common/dft_parameters
xc.o:		../m4/common/errcodes
xc.o:		../m4/common/mapper
xc.o:		../m4/common/parallel
xc.o:		../m4/common/parcntl
xc.o:		../m4/common/sizes
xc.o:		../m4/common/timeperiods
xc.o:		common/dft_api
xc.o:		common/dft_basder
xc.o:		common/dft_basis
xc.o:		common/dft_basis_api
xc.o:		common/dft_dfder
xc.o:		common/dft_dij
xc.o:		common/dft_module_comm
xc.o:		common/dft_mol_info
xc.o:		common/dft_numbers
xc.o:		common/dft_order_info
xc.o:		common/dft_parameters
xc.o:		common/dft_parameters
xc.o:		common/dft_physical_constants
xc.o:		common/dft_xc
xc_fit.o:		common/dft_basis
xc_fit.o:		common/dft_basis_api
xc_fit.o:		common/dft_basis_cont_inf
xc_fit.o:		common/dft_module_comm
xc_fit.o:		common/dft_mol_info
xc_fit.o:		common/dft_numbers
xc_fit.o:		common/dft_order_info
xc_fit.o:		common/dft_parameters
xc_lib.o:		common/citation
xc_lib.o:		common/dft_dfder
xc_lib.o:		common/dft_xcfunc
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
MPOLE=	omit_mpole.o

#
# currently exluded
#xc_fit.o \
#geometry.o \
#

DFTOBJ = dftbasis.o \
	coulomb.o \
	gamess.o \
	global.o \
	interface.o \
	intpack.o \
	matpack.o \
	xc.o \
        xc_lib.o \
        weights.o \
        exp_dksm_hess.o \
        chf.o \
        Lebedev-Laikov.o \
        integ_data.o \
        integ_te2c_norm.o \
	readbasis.o \
	jfitg.o \
	integ2e.o \
	deriv2e.o \
	matform.o \
	gden.o $(GADFT) $(EXTRA_DFT) \
	${MPOLE}

DFTLIB = libdft.a

dftlib:	$(DFTLIB)

$(DFTLIB):	$(DFTOBJ)
		$(ARCHIVE) $(DFTLIB) $(DFTOBJ)
		$(RANLIB) $(DFTLIB)


.DEFAULT:
	(cd ../m4; make $@)

clean:
	touch x.c x.o x~ x.f \#x
	rm -f *.c *.o *.f *~ \#* $(DFTLIB)


gadft.o:	gadft.m
	rm -f gadft.GUK gadft.f
	cat  ../utilities/gener.m gadft.m | $(M4) $(M4OPTS) $(SNGL) >> gadft.f
	$(FC) $(FFLAGSS) gadft.f
