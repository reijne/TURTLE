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
EXE= ../bin/gamess-uk
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
MAKE_GA=make_ga
GA_DIR=../g
GA_F77_DEFS =
GA_CC = ${CC} ${CFLAGS}
GA_F77_INC= -I$(GA_DIR)/include 
GA_C_INC= -I$(GA_DIR)/include 
GA_TARGET=You_must_provide_a_GA_TARGET_setting
GA_TARGET_CPU_PAR=
GA_VERSION_PAR=
GA_SCALAPACK=
#
# NB - assumption is here is that we are always running with
# mpi, we can add tcgmsg into the mix if we need to run with i8
#
# GA build using MPI directly for communications
GA_MPI=MSG_COMMS=MPI
GA_LIBS=\
	${GA_DIR}/lib/${GA_TARGET}/libglobal.a \
        ${GA_DIR}/lib/${GA_TARGET}/libarmci.a \
        ${GA_DIR}/lib/${GA_TARGET}/libma.a
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
VB_EXTRA=casa.o casb.o  
# ----------------------- Van der Waals Options -----------------------------
VDW_M4_OPTS=,vdw
VDW_EXTRA=vdwaals_corr.o inter_vdwaals.o
# ----------------------- Sysmo Options -----------------------------
SYSMO_M4_OPTS=
SYSMO_EXTRA=
# ----------------------- Parallel MSCSF (MASSCF) Options -----------------------------
MASSCF_M4_OPTS=,masscf
MASSCF_EXTRA=masscf.o ormas1.o ormas2.o
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
SCALAPACK_M4_OPTS=
#
#
PDIAG_M4_OPTS=
#
# Sort out the functionality we are including to set PAR_M4_BUILD
#
PAR_M4_BUILD=,mp2_parallel,secd_parallel
#
#
# Now determine the parallel library options
PAR_M4_OPTS=parallel,ma,ga,mpi${PAR_M4_BUILD}
# End of GA options
#
# End of various parallel options - now just serial
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
D2=,zora
D6=,vb
D11=,vdw
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
#doc This is the machine-specific file for the GLOBAL ARRAY based parallel version
#doc of GAMESS-UK on Linux running on Opteron processors with PGI compilers (tested with 5.2.)
#doc
#doc Options:
#doc datain - force GAMESS-UK to read it's input from a file called datain instead of standard input
#doc score  - build on an SCORE system (WARNING - this build only works on SCORE 5.8.4 & above)
#doc mpiwrap - use mpi wrappers scripts to locate mpi libraries and headers
#doc newscf - include distributed data scf module
#doc scalapack - include scalapack interface if the library is available
#doc blas   - inlclude external blas libraries.
#doc myrinet - build against MX drivers on a Myrinet switch (GM not supported yet)
#doc static_lb - statically load-balanced version (no processor set aside as nxtval server)
#
# Default options:
#dopt ga mpi ci peigs newscf dl-find vdw masscf mpiwrap
#opt score blas myrinet datain i8 scalapack static_lb vb zora
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=opteron,linux,littleendian,cio,unix,doublebackslash,upck-equiv,extpopcnt,extleadz

# Bitwise operators that differ from the defaults.
IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)
#
#
RANLIB=ranlib
CPP=/lib/cpp

MPI_INCLUDE  = /opt/pgi/linux86-64/2019/mpi/openmpi-3.1.3/include/
IMPI_INCLUDE= -I${MPI_INCLUDE}
MPI_LIB      = /opt/pgi/linux86-64/2019/mpi/openmpi-3.1.3/lib/
LIBMPI =  -L$(MPI_LIB) -lmpi -lmpi_mpifh

GM_INCLUDE =
IGM_INCLUDE =
GM_LIB =
LIBGM =

#
# ===============  Compiler Options
#

FC = pgf90
FC90 = pgf90
LD = pgf90
CC = pgcc

#Default compilation Flags

FFLAGSI8 = -i8
CFLAGSI8 = -DLONG_INTEGER -DEXT_INT
I8_M4_OPT=,i8

FFLAGSD = -c -tp=k8-64 -Mcache_align ${FFLAGSI8} ${IMPI_INCLUDE} ${IGM_INCLUDE}
CFLAGSD= -c -tp=k8-64 ${CFLAGSI8}  ${IMPI_INCLUDE} ${IGM_INCLUDE}

LDFLAGSD= -tp=k8-64 
#

FFLAGSV =  ${FFLAGSD}  -O2
FFLAGSS =  ${FFLAGSD}  -O
FFLAGSN =  ${FFLAGSD}  -O1
FFLAGSN0 = ${FFLAGSD}
CFLAGS  = -D_REENTRANT ${CFLAGSD}
LDFLAGS = ${LDFLAGSD}
#
#




# PEIGS definitions
PEIGS_TARGET_CPU_PAR=x86_64

# Global Array definitions
GA_F77_DEFS = -traditional
GA_TARGET=LINUX64
# May need to use GCC for compiling the GAs
#GA_CC = gcc
# Default MPI/TCP/IP socket code
GA_VERSION_PAR= GA_VERSION=SHMEM MPI_INCLUDE=${MPI_INCLUDE} MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "


# Bring all the options together
OPTIONS=${MACHOPT}${BLASOPT}${I8_M4_OPT}

BL_LIB = ${LIB_SCALAPACK} ${LIBBLAS} ${LIBMPI} ${LIBGM}

#
# ===============  Additional Files
#
EXTRA=
EXTRA_MP2= check0a.o aprq34d.o mcdab_ga.o
EXTRA_DFT = jkint_dft.o jkder_dft.o
#
# ===============  Compiler Exceptions
#

check0a.o: newmrd5.m
	cat ../utilities/gener.m newmrd5.m | $(M4) -DGEN_EXTRACTFILE=check0a $(M4OPTS) > check0a.f
	$(FC) $(FFLAGSN) check0a.f

# The below are extracted by the linux keyword and are probably not required for this compiler
#
aprq34d.o: mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=aprq34d $(M4OPTS) > aprq34d.f
	$(FC) $(FFLAGSS) aprq34d.f

mcdab_ga.o: mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=mcdab_ga $(M4OPTS) > mcdab_ga.f
	$(FC) $(FFLAGSS) mcdab_ga.f

#
#  ========== DFT Exceptions (PGF) ===============
#
# The below are extracted by the opteron adn em64t keywords and are probably not required for this compiler
#
jkint_dft.o:    integ2e.m
	cat ../utilities/gener.m  integ2e.m | $(M4) -DGEN_EXTRACTFILE=jkint_dft $(M4OPTS) $(SNGL) > jkint_dft.f
	$(FC) $(FFLAGSN) jkint_dft.f
#
jkder_dft.o: deriv2e.m
	cat ../utilities/gener.m  deriv2e.m | $(M4) -DGEN_EXTRACTFILE=jkder_dft $(M4OPTS) $(SNGL) > jkder_dft.f
	$(FC) $(FFLAGSN) jkder_dft.f

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
	@echo Building CI code

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
	mcscfa.o mcscfb.o mcscfc.o \
        mrdci1.o mrdci2.o mrdci3.o mrdci4.o mrdci5.o mrdci6.o mrdci7.o fullci.o \
        direct.o casa.o casb.o \
        ccsd.o dirrpa.o mclr.o nvccsd.o rpa.o rpagrad.o tsort.o \
	tsortc.o 

# Construct list of object files

#
#  MP2 + SCF Seconds + direct-CI + GF
#
TOBJ = ${OBJ_BASE} ${OBJ_GA_MP2} ${OBJ_GA_SCF} ${OBJ_CI} $(MASSCF_EXTRA)
OMIT = ,ga_ci_build


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
	cat  ../utilities/gener.m parallel.m | $(M4) $(M4OPTS)  > parallel.GUK
	$(CPP) $(GA_F77_INC) $(GA_F77_DEFS)  < parallel.GUK | sed '/^#/D' >> parallel.f
	rm -f parallel.GUK
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
	cat ../utilities/gener.m parallel.m | $(M4) $(M4OPTS)  > parallel.GUK
	$(CPP) $(GA_F77_INC) $(GA_F77_DEFS)  < parallel.GUK | sed '/^#/D' >> parallel.f
	rm -f parallel.GUK
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
gwmix.o:		common/tractlt
servec.o:		../m4/common/blkorbs
servec.o:		../m4/common/discc
servec.o:		../m4/common/dump3
servec.o:		../m4/common/infoa
servec.o:		../m4/common/iofile
servec.o:		../m4/common/mapper
servec.o:		../m4/common/runlab
servec.o:		../m4/common/sector
servec.o:		../m4/common/sizes
servec.o:		../m4/common/tran
servec.o:		../m4/common/work
servec.o:		../m4/common/workc
vbaivb.o:		../m4/common/infoa
vbaivb.o:		../m4/common/iofile
vbaivb.o:		../m4/common/runlab
vbaivb.o:		../m4/common/sizes
vbaivb.o:		common/aivb
vbaivb.o:		common/turtleparam
vbci.o:		../m4/common/gmemdata
vbci.o:		../m4/common/gmempara
vbci.o:		../m4/common/iofile
vbci.o:		../m4/common/restri
vbci.o:		../m4/common/sizes
vbci.o:		../m4/common/timeperiods
vbci.o:		../m4/common/vcore
vbci.o:		common/brill
vbci.o:		common/c8_16vb
vbci.o:		common/ffile
vbci.o:		common/gainfo
vbci.o:		common/gainfo
vbci.o:		common/hsinfo
vbci.o:		common/scftvb
vbci.o:		common/splice
vbci.o:		common/tractlt
vbci.o:		common/turtleparam
vbci.o:		common/twice
vbci.o:		common/vbcri
vbci.o:		common/vbparr
vbci.o:		common/vbpert
vbci.o:		common/vbtimess
vbci.o:		common/vbvirt
vbcrestr.o:		../m4/common/infoa
vbcrestr.o:		../m4/common/iofile
vbcrestr.o:		../m4/common/nshel
vbcrestr.o:		../m4/common/sizes
vbcrestr.o:		../m4/common/timeperiods
vbcrestr.o:		../m4/common/work
vbcrestr.o:		common/aivb
vbcrestr.o:		common/allconf
vbcrestr.o:		common/c8_16vb
vbcrestr.o:		common/c8_16vb
vbcrestr.o:		common/first_vb
vbcrestr.o:		common/logics
vbcrestr.o:		common/scra7vb
vbcrestr.o:		common/splice
vbcrestr.o:		common/tractlt
vbcrestr.o:		common/turtleparam
vbcrestr.o:		common/vbcri
vbdebug.o:		common/turtleparam
vbdens.o:		../m4/common/atmblk
vbdens.o:		../m4/common/atmol3
vbdens.o:		../m4/common/cndx40
vbdens.o:		../m4/common/dump3
vbdens.o:		../m4/common/files
vbdens.o:		../m4/common/funct
vbdens.o:		../m4/common/gjs
vbdens.o:		../m4/common/gmemdata
vbdens.o:		../m4/common/gmempara
vbdens.o:		../m4/common/gmempara
vbdens.o:		../m4/common/infoa
vbdens.o:		../m4/common/iofile
vbdens.o:		../m4/common/mapper
vbdens.o:		../m4/common/qice
vbdens.o:		../m4/common/restar
vbdens.o:		../m4/common/restri
vbdens.o:		../m4/common/runlab
vbdens.o:		../m4/common/sizes
vbdens.o:		../m4/common/timeperiods
vbdens.o:		../m4/common/vcore
vbdens.o:		../m4/common/vcore
vbdens.o:		common/brill
vbdens.o:		common/ffile
vbdens.o:		common/infato
vbdens.o:		common/scftvb
vbdens.o:		common/splice
vbdens.o:		common/tractlt
vbdens.o:		common/turtleparam
vbdens.o:		common/vbcri
vbdens.o:		common/vblimit
vbdens.o:		common/vbtimess
vbgamess.o:		../m4/common/atmol3
vbgamess.o:		../m4/common/blkorbs
vbgamess.o:		../m4/common/blksiz
vbgamess.o:		../m4/common/cndx40
vbgamess.o:		../m4/common/discc
vbgamess.o:		../m4/common/dm
vbgamess.o:		../m4/common/dump3
vbgamess.o:		../m4/common/files
vbgamess.o:		../m4/common/gjs
vbgamess.o:		../m4/common/infoa
vbgamess.o:		../m4/common/iofile
vbgamess.o:		../m4/common/machin
vbgamess.o:		../m4/common/mapper
vbgamess.o:		../m4/common/prnprn
vbgamess.o:		../m4/common/qice
vbgamess.o:		../m4/common/restar
vbgamess.o:		../m4/common/runlab
vbgamess.o:		../m4/common/scfopt
vbgamess.o:		../m4/common/scra7
vbgamess.o:		../m4/common/sector
vbgamess.o:		../m4/common/sizes
vbgamess.o:		../m4/common/statis
vbgamess.o:		../m4/common/tran
vbgamess.o:		../m4/common/trntim
vbgamess.o:		../m4/common/work
vbgamess.o:		common/scftvb
vbgamess.o:		common/tractlt
vbin.o:		../m4/common/atmol3
vbin.o:		../m4/common/common
vbin.o:		../m4/common/dm
vbin.o:		../m4/common/dump3
vbin.o:		../m4/common/gmemdata
vbin.o:		../m4/common/gmempara
vbin.o:		../m4/common/harmon
vbin.o:		../m4/common/infoa
vbin.o:		../m4/common/iofile
vbin.o:		../m4/common/mapper
vbin.o:		../m4/common/prints
vbin.o:		../m4/common/restar
vbin.o:		../m4/common/restri
vbin.o:		../m4/common/restrj
vbin.o:		../m4/common/runlab
vbin.o:		../m4/common/scra7
vbin.o:		../m4/common/sector
vbin.o:		../m4/common/sizes
vbin.o:		../m4/common/sizes
vbin.o:		../m4/common/timeperiods
vbin.o:		../m4/common/tran
vbin.o:		../m4/common/vcore
vbin.o:		common/aivb
vbin.o:		common/basisvb
vbin.o:		common/brill
vbin.o:		common/c8_16vb
vbin.o:		common/dumpvbmo
vbin.o:		common/ffile
vbin.o:		common/first_vb
vbin.o:		common/hsinfo
vbin.o:		common/hsmattype
vbin.o:		common/infato
vbin.o:		common/scftvb
vbin.o:		common/scra7vb
vbin.o:		common/splice
vbin.o:		common/tractlt
vbin.o:		common/turtleparam
vbin.o:		common/twice
vbin.o:		common/vbbonds
vbin.o:		common/vbcri
vbin.o:		common/vbdiist
vbin.o:		common/vbequiv
vbin.o:		common/vbparr
vbin.o:		common/vbpert
vbin.o:		common/vbproper
vbin.o:		common/vbqc
vbin.o:		common/vbtimess
vbin.o:		common/vbvirt
vbmatre.1.o:		../m4/common/gmemdata
vbmatre.1.o:		../m4/common/gmempara
vbmatre.1.o:		../m4/common/iofile
vbmatre.1.o:		../m4/common/iofile
vbmatre.1.o:		../m4/common/qice
vbmatre.1.o:		../m4/common/restri
vbmatre.1.o:		../m4/common/sizes
vbmatre.1.o:		common/brill
vbmatre.1.o:		common/c8_16vb
vbmatre.1.o:		common/ffile
vbmatre.1.o:		common/gainfo
vbmatre.1.o:		common/hsinfo
vbmatre.1.o:		common/parinf
vbmatre.1.o:		common/scftvb
vbmatre.1.o:		common/splice
vbmatre.1.o:		common/tractlt
vbmatre.1.o:		common/turtleparam
vbmatre.1.o:		common/twice
vbmatre.1.o:		common/vbcri
vbmatre.1.o:		common/vblimit
vbmatre.1.o:		common/vbpert
vbmatre.gpu.o:		../m4/common/gmemdata
vbmatre.gpu.o:		../m4/common/gmempara
vbmatre.gpu.o:		../m4/common/iofile
vbmatre.gpu.o:		../m4/common/iofile
vbmatre.gpu.o:		../m4/common/qice
vbmatre.gpu.o:		../m4/common/restri
vbmatre.gpu.o:		../m4/common/sizes
vbmatre.gpu.o:		common/brill
vbmatre.gpu.o:		common/c8_16vb
vbmatre.gpu.o:		common/ffile
vbmatre.gpu.o:		common/gainfo
vbmatre.gpu.o:		common/hsinfo
vbmatre.gpu.o:		common/parinf
vbmatre.gpu.o:		common/scftvb
vbmatre.gpu.o:		common/splice
vbmatre.gpu.o:		common/tractlt
vbmatre.gpu.o:		common/turtleparam
vbmatre.gpu.o:		common/twice
vbmatre.gpu.o:		common/vbcri
vbmatre.gpu.o:		common/vblimit
vbmatre.gpu.o:		common/vbpert
vbmatre.o:		../m4/common/gmemdata
vbmatre.o:		../m4/common/gmempara
vbmatre.o:		../m4/common/iofile
vbmatre.o:		../m4/common/iofile
vbmatre.o:		../m4/common/qice
vbmatre.o:		../m4/common/restri
vbmatre.o:		../m4/common/sizes
vbmatre.o:		common/brill
vbmatre.o:		common/c8_16vb
vbmatre.o:		common/ffile
vbmatre.o:		common/gainfo
vbmatre.o:		common/hsinfo
vbmatre.o:		common/parinf
vbmatre.o:		common/scftvb
vbmatre.o:		common/splice
vbmatre.o:		common/tractlt
vbmatre.o:		common/turtleparam
vbmatre.o:		common/twice
vbmatre.o:		common/vbcri
vbmatre.o:		common/vblimit
vbmatre.o:		common/vbpert
vbmatre.orig.o:		../m4/common/gmemdata
vbmatre.orig.o:		../m4/common/gmempara
vbmatre.orig.o:		../m4/common/iofile
vbmatre.orig.o:		../m4/common/iofile
vbmatre.orig.o:		../m4/common/qice
vbmatre.orig.o:		../m4/common/restri
vbmatre.orig.o:		../m4/common/sizes
vbmatre.orig.o:		common/brill
vbmatre.orig.o:		common/c8_16vb
vbmatre.orig.o:		common/ffile
vbmatre.orig.o:		common/gainfo
vbmatre.orig.o:		common/hsinfo
vbmatre.orig.o:		common/parinf
vbmatre.orig.o:		common/scftvb
vbmatre.orig.o:		common/splice
vbmatre.orig.o:		common/tractlt
vbmatre.orig.o:		common/turtleparam
vbmatre.orig.o:		common/twice
vbmatre.orig.o:		common/vbcri
vbmatre.orig.o:		common/vblimit
vbmatre.orig.o:		common/vbpert
vbmatre.v1.o:		../m4/common/gmemdata
vbmatre.v1.o:		../m4/common/gmempara
vbmatre.v1.o:		../m4/common/iofile
vbmatre.v1.o:		../m4/common/iofile
vbmatre.v1.o:		../m4/common/qice
vbmatre.v1.o:		../m4/common/restri
vbmatre.v1.o:		../m4/common/sizes
vbmatre.v1.o:		common/brill
vbmatre.v1.o:		common/c8_16vb
vbmatre.v1.o:		common/ffile
vbmatre.v1.o:		common/gainfo
vbmatre.v1.o:		common/hsinfo
vbmatre.v1.o:		common/parinf
vbmatre.v1.o:		common/scftvb
vbmatre.v1.o:		common/splice
vbmatre.v1.o:		common/tractlt
vbmatre.v1.o:		common/turtleparam
vbmatre.v1.o:		common/twice
vbmatre.v1.o:		common/vbcri
vbmatre.v1.o:		common/vblimit
vbmatre.v1.o:		common/vbpert
vbmatre.v2.o:		../m4/common/gmemdata
vbmatre.v2.o:		../m4/common/gmempara
vbmatre.v2.o:		../m4/common/iofile
vbmatre.v2.o:		../m4/common/iofile
vbmatre.v2.o:		../m4/common/qice
vbmatre.v2.o:		../m4/common/restri
vbmatre.v2.o:		../m4/common/sizes
vbmatre.v2.o:		common/brill
vbmatre.v2.o:		common/c8_16vb
vbmatre.v2.o:		common/ffile
vbmatre.v2.o:		common/gainfo
vbmatre.v2.o:		common/hsinfo
vbmatre.v2.o:		common/parinf
vbmatre.v2.o:		common/scftvb
vbmatre.v2.o:		common/splice
vbmatre.v2.o:		common/tractlt
vbmatre.v2.o:		common/turtleparam
vbmatre.v2.o:		common/twice
vbmatre.v2.o:		common/vbcri
vbmatre.v2.o:		common/vblimit
vbmatre.v2.o:		common/vbpert
vbprop.o:		../m4/common/iofile
vbprop.o:		../m4/common/iofile
vbprop.o:		../m4/common/mapper
vbprop.o:		../m4/common/restri
vbprop.o:		../m4/common/sizes
vbprop.o:		common/brill
vbprop.o:		common/c8_16vb
vbprop.o:		common/c8_16vb
vbprop.o:		common/ffile
vbprop.o:		common/infato
vbprop.o:		common/parinf
vbprop.o:		common/scftvb
vbprop.o:		common/splice
vbprop.o:		common/tractlt
vbprop.o:		common/turtleparam
vbprop.o:		common/twice
vbprop.o:		common/vblimit
vbprop.o:		common/vbproper
vbqcscf.o:		../m4/common/gmemdata
vbqcscf.o:		../m4/common/gmempara
vbqcscf.o:		../m4/common/gmempara
vbqcscf.o:		../m4/common/iofile
vbqcscf.o:		../m4/common/mapper
vbqcscf.o:		../m4/common/prints
vbqcscf.o:		../m4/common/restri
vbqcscf.o:		../m4/common/runlab
vbqcscf.o:		../m4/common/sizes
vbqcscf.o:		../m4/common/sizes
vbqcscf.o:		../m4/common/timeperiods
vbqcscf.o:		common/aivb
vbqcscf.o:		common/brill
vbqcscf.o:		common/c8_16vb
vbqcscf.o:		common/ffile
vbqcscf.o:		common/first_vb
vbqcscf.o:		common/hsinfo
vbqcscf.o:		common/hsmattype
vbqcscf.o:		common/infato
vbqcscf.o:		common/parinf
vbqcscf.o:		common/scftvb
vbqcscf.o:		common/splice
vbqcscf.o:		common/tractlt
vbqcscf.o:		common/turtleparam
vbqcscf.o:		common/twice
vbqcscf.o:		common/vbcri
vbqcscf.o:		common/vbdiist
vbqcscf.o:		common/vbequiv
vbqcscf.o:		common/vbpert
vbqcscf.o:		common/vbproper
vbqcscf.o:		common/vbqc
vbqcscf.o:		common/vbtimess
vbscf.o:		../m4/common/atmol3
vbscf.o:		../m4/common/files
vbscf.o:		../m4/common/gmemdata
vbscf.o:		../m4/common/gmempara
vbscf.o:		../m4/common/gmempara
vbscf.o:		../m4/common/harmon
vbscf.o:		../m4/common/infoa
vbscf.o:		../m4/common/iofile
vbscf.o:		../m4/common/mapper
vbscf.o:		../m4/common/restri
vbscf.o:		../m4/common/runlab
vbscf.o:		../m4/common/sizes
vbscf.o:		../m4/common/timeperiods
vbscf.o:		../m4/common/tran
vbscf.o:		../m4/common/vcore
vbscf.o:		../m4/common/vcore
vbscf.o:		common/aivb
vbscf.o:		common/basisvb
vbscf.o:		common/brill
vbscf.o:		common/c8_16vb
vbscf.o:		common/c8_16vb
vbscf.o:		common/ffile
vbscf.o:		common/first_vb
vbscf.o:		common/gainfo
vbscf.o:		common/hsinfo
vbscf.o:		common/hsmattype
vbscf.o:		common/infato
vbscf.o:		common/scftvb
vbscf.o:		common/scra7vb
vbscf.o:		common/splice
vbscf.o:		common/tractlt
vbscf.o:		common/turtleparam
vbscf.o:		common/turtleparam
vbscf.o:		common/twice
vbscf.o:		common/vbcri
vbscf.o:		common/vbdiist
vbscf.o:		common/vbequiv
vbscf.o:		common/vbpert
vbscf.o:		common/vbproper
vbscf.o:		common/vbtimess
vbscf.o:		common/vbvirt
vbtran.o:		../m4/common/atmblk
vbtran.o:		../m4/common/blksiz
vbtran.o:		../m4/common/cndx40
vbtran.o:		../m4/common/files
vbtran.o:		../m4/common/gmemdata
vbtran.o:		../m4/common/gmempara
vbtran.o:		../m4/common/infoa
vbtran.o:		../m4/common/iofile
vbtran.o:		../m4/common/mapper
vbtran.o:		../m4/common/psscrf
vbtran.o:		../m4/common/restri
vbtran.o:		../m4/common/sizes
vbtran.o:		../m4/common/sizes
vbtran.o:		../m4/common/vcore
vbtran.o:		common/basisvb
vbtran.o:		common/c8_16vb
vbtran.o:		common/c8_16vb
vbtran.o:		common/ffile
vbtran.o:		common/scftvb
vbtran.o:		common/splice
vbtran.o:		common/tractlt
vbtran.o:		common/turtleparam
vbtran.o:		common/vbcri
vbtran.o:		common/vbparr
vbtran.o:		ikyg.o
vbutil.o:		../m4/common/iofile
vbutil.o:		../m4/common/vcore
vbutil.o:		common/ffile
vbutil.o:		common/gainfo
vbutil.o:		common/logics
vbutil.o:		common/scra7vb
vbutil.o:		common/turtleparam
vbutil.o:		common/vbcri
#
# ============================ NON-COMPILED TARGETS ===========================
#
#




ga_clean:
	@(if test -d ${GA_DIR}; then \
	cd ${GA_DIR}; $(GMAKE) TARGET=$(GA_TARGET) clean; fi)
	\rm -f ../bin/parallel


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
make_gatest:	check_ga make_ga $(GA_DIR)/global/testing/test.x
	ls -l $(GA_DIR)/global/testing/test.x

make_ga:	check_ga $(GA_LIBS)
	@(if test -f $(GA_DIR)/tcgmsg/ipcv4.0/parallel; then \
	cp $(GA_DIR)/tcgmsg/ipcv4.0/parallel ../bin ; fi )

GA_DIAG_FLAG=

$(GA_LIBS):
	@(echo Making Global Array Tools)
	(cd $(GA_DIR); $(GMAKE) TARGET=$(GA_TARGET) FC="$(FC)" CC="$(GA_CC)" DIAG=$(GA_DIAG_FLAG) $(GA_TARGET_CPU_PAR) $(GA_VERSION_PAR) $(GA_MPI) $(GA_SCALAPACK))

$(GA_DIR)/global/testing/test.x:
	@(echo Making Global Array Test)
	@(cd $(GA_DIR); $(GMAKE) TARGET=$(GA_TARGET) FC="$(FC)" CC="$(GA_CC)" DIAG=$(GA_DIAG_FLAG) $(GA_TARGET_CPU_PAR) $(GA_VERSION_PAR))

#
check_ga:	check_gmake
	@(if test -d $(GA_DIR); then \
	echo Checking for GA Tools - directory $(GA_DIR) seems OK .. ;\
	else \
	echo '      ================================================='; \
	echo '      === GA TOOLS ARE NEEDED TO BUILD THIS VERSION ==='; \
	echo '      ===             OF GAMESS-UK                  ==='; \
	echo '      ================================================='; \
	exit -1 ;fi)
check_gmake:
	@(if $(GMAKE) -v > /dev/null; then \
	echo Checking for gmake - seems to be installed OK .. ;\
	else \
	echo '      ==============================================================='; \
	echo '      === ERROR:  gmake IS NEEDED TO BUILD THE GA TOOLS and PEIGS ==='; \
	echo '      ==============================================================='; \
	exit -1; fi)

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
VBOBJ = vbaivb.o\
	vbci.o\
	vbdebug.o\
	vbin.o\
	vbscf.o\
	vbqcscf.o\
	vbutil.o\
	vbcrestr.o\
	vbgamess.o\
	vbmatre.o\
	vbtran.o\
	vbdens.o\
	vbprop.o\
	vbversion.o\
	servec.o\
	sendrecv.o

VBLIB = libvb.a

vb:	$(VBLIB)
	

$(VBLIB):	$(VBOBJ) $(VBEXTRA)
		$(ARCHIVE) $(VBLIB)  $(VBOBJ)  $(VBEXTRA)
		$(RANLIB) $(VBLIB)


vbversion.f:	vbaivb.m vbci.m vbdebug.m vbin.m vbutil.m vbcrestr.m vbgamess.m vbmatre.m vbtran.m vbdens.m vbqcscf.m servec.m
		echo "      subroutine vbversion(cdate,ctime,cname)" > vbversion.f
		echo "      character*10 cdate,cname" >> vbversion.f
		echo "      character*5 ctime" >> vbversion.f
		echo "      cdate=\""`date '+%d-%m-%Y'`"\"" >> vbversion.f
		echo "      ctime=\""`date '+%H:%M'`"\"" >> vbversion.f
		echo "      cname=\""`${WHOAMI}`"\"" >> vbversion.f
		echo "      return" >> vbversion.f
		echo "      end" >> vbversion.f

# Only include MPI_INCLUDE dir if is has been set
sendrecv.o:	sendrecv.c
	if  [ "${MPI_INCLUDE}" ] ; \
	then $(CC) $(CFLAGS)  -I${MPI_INCLUDE} -I../g/tcgmsg-mpi $(OBJNAME) sendrecv.c ; \
	else $(CC) $(CFLAGS)  -I../g/tcgmsg-mpi $(OBJNAME) sendrecv.c ; \
	fi

clean:
	touch x.o x~ x.f \#x
	rm -f *.o *.f *.c *~ \#* $(VBLIB) 

.DEFAULT: vb
	
#	(cd ../m4; make $@)
#
