#doc  Makefile settings for Global Array build on the SGI Altix (Itanium2)
#doc  with the Intel compiler 8.1  and MPI libraries in /usr/lib
#doc  do not setenv USE_MPI with tcgmsg-mpi
#doc
#doc  Options:
#doc  mkl - use Intel MKL numerical library (default location: /opt/modules/intel/mkl721/lib/64)
#doc        The use of MKL requires the  path to the Intel compiler libraries to be set,
#doc        with the default being: /opt/Modules/cmplrs/intel/8.1_new/lib
#doc  scs - use the SGI scs numerical library instead of mkl
#doc  suse - use sara propack4 (suse) setup
#
# DEFAULT OPTIONS
#dopt ga mpi i8 ci vb zora peigs mkl vdw dl-find 
#opt scs suse newscf
#
# ================ M4 Processing options
#

# Machine-specific options 
MACHINE_KEY=G
MACHOPT=itanium,linux,littleendian,cio,unix,upck-equiv,altix

IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)

#General Libraries
# lsma in /usr/lib and is for the shared memory stuff
LLIB= -lsma

#MPI flags
MPI_INCLUDE =/usr/include
MPI_LIB =/usr/lib
LIBMPI= -lmpi
MPI_LIBS= -L${MPI_LIB} ${LIBMPI}

#--#if suse#
MPI_INCLUDE = /opt/Modules/mpt-ia64/sgi-mpt-1.9-1/include
# MPI_LIB = /opt/Modules/mpt-ia64/sgi-mpt-1.9-1/lib
MPI_LIB = 
#--#endif#


#--#if mkl#
# only really need -lthread for static linking
# -lguide is used to pick the threading library for dynamic linking
Lmkl=-L/opt/modules/intel/mkl721/lib/64/
lmkl=-lmkl_ipf -lpthread
Lintel=-L/opt/modules/cmplrs/intel/8.1.023/lib
lintel=-lguide
LIBBLAS=${Lintel} ${lintel} ${Lmkl} ${lmkl}
BLASOPT=,blas
#--#elseif scs#
#----#if i8#
LIBBLAS=-lscs_i8
#----#else#
LIBBLAS=-lscs
#----#endif#
BLASOPT=,blas
#--#endif mkl scs#

#GA Stuff
GA_F77_DEFS = -traditional
GA_TARGET=LINUX64
GA_VERSION_PAR= GA_VERSION=SHMEM MPI_INCLUDE=${MPI_INCLUDE}  MPI_LIB=${MPI_LIB} LIBMPI=" ${LIBMPI} "
GA_TARGET_CPU_PAR= GA_TARGET_CPU=SGIALTIX

# Put the link line together
BL_LIB= ${LLIB} ${lLIB} ${MPI_LIBS} ${LIBBLAS}

## M4 Options ##
#--#if i8#
OPTIONS=${MACHOPT}${BLASOPT},i8
#--#else#
OPTIONS=${MACHOPT}${BLASOPT},i8drct,64bitpointers
#--#endif i8#

#
# ===============  Compiler Options
#
FC = ifort
LD = ifort
CC = icc -no-gcc
FC90 = ifort
CXX = icc
RANLIB = ranlib

#--#if i8#
FFLAGI8= -i8
CFLAGI8= -DLONG_INTEGER
#--#else#
FFLAGI8=
CFLAGI8= -DSTD_INT
#--#endif i8#
FFLAGSTMP = -c ${FFLAGI8}
FINC = -I${MPI_INCLUDE}

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g ${FINC}
FFLAGSS = ${FFLAGSTMP} -g ${FINC}
FFLAGSN = ${FFLAGSTMP} -g ${FINC}
FFLAGSN0 =  ${FFLAGSTMP} -g 
CFLAGS = -g -c ${CFLAGI8}
LDFLAGS  = -g
CXXFLAGS = -g -c
#--#else# 
FFLAGSV = ${FFLAGSTMP} -O2 -cm -w90 -w95 -ip -ftz ${FINC}
FFLAGSS = ${FFLAGSTMP} -O2 -cm -w90 -w95 -ftz ${FINC}
FFLAGSN = ${FFLAGSTMP} -O1 -cm -w90 -w95 -ftz ${FINC}
FFLAGSN0 = ${FFLAGSTMP} -O0
CFLAGS  = -c ${CFLAGI8}
LDFLAGS =  
CXXFLAGS = -O2 -c
#--#endif#

# Diesel build options
LD_DIESEL = icpc
#--#if mkl scs#
DIESEL_LIBS = ${LBLAS} ${lBLAS} -lifcore
#--#else#
DIESEL_LIBS = -lifcore
#--#endif#
FLEX = /usr/bin/flex
__UNDERBAR = 1
SIZEOF_VOID_P = 8
#--#if i8#
LONG_LONG_INT = long int
LONG_INT = long int
INT = long int
SHORT_INT = short int
#--#else#
LONG_LONG_INT = long long int
LONG_INT = long int
INT = int
SHORT_INT = short int
#--#endif i8#

#
# ===============  Additional Files 
#
EXTRA=
EXTRA_MP2= aprq34d.o mcdab_ga.o tst1s.o

#
# ===============  Compiler Exceptions
#
tst1s.o:	direct.m
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f

### Below exceptions required for Intel 10.1 
# Compiling with o2 with intel 10.1 causes numerous errors
integs.o: integs.f
	$(FC) $(FFLAGSN) integs.f

rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSS) $*.f

check0a.o:	newmrd5.m
	cat ../utilities/gener.m newmrd5.m | $(M4) -DGEN_EXTRACTFILE=check0a $(M4OPTS) > check0a.f
	$(FC) $(FFLAGSN) check0a.f
aprq34d.o:	mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=aprq34d $(M4OPTS) > aprq34d.f
	$(FC) $(FFLAGSS) aprq34d.f

mcdab_ga.o:	mp2_parallel.m
	cat ../utilities/gener.m mp2_parallel.m | $(M4) -DGEN_EXTRACTFILE=mcdab_ga $(M4OPTS) > mcdab_ga.f
	$(FC) $(FFLAGSS) mcdab_ga.f
#
# Information on the machine this build was configured on.
# Output of uname -a
# Linux newtd 2.4.21-sgi304rp05031014_10149 #1 SMP Thu Mar 10 14:43:36 PST 2005 ia64 ia64 ia64 GNU/Linux
# 
# Output of cat /proc/cpuinfo (for first processor):
# 
# processor  : 0
# vendor     : GenuineIntel
# arch       : IA-64
# family     : Itanium 2
# model      : 1
# revision   : 5
# archrev    : 0
# features   : branchlong
# cpu number : 0
# cpu regs   : 4
# cpu MHz    : 1300.000000
# itc MHz    : 1300.000000
# BogoMIPS   : 1071.64
# 
# Output of hinv:
# 
# NOTICE: non root user will produce less information
# Serial number: N0000767
# 1 IX-Brick
# 63 C-Brick
# 40 R-Brick
# 
# 128 1500 MHz Itanium 2 Rev. 5 Processor
# 124 1300 MHz Itanium 2 Rev. 5 Processor
# Main memory size: 479.39 Gb
# Co-processor: Silicon Graphics, Inc. IOC4 I/O controller (rev 4f) on pci01.01.0
#    QLogic 12160 Dual Channel Ultra3 SCSI (Rev 06) on pci01.03.0 on pci01.03.0
#       Disk SGI ST336753LC 36 GB
#       Disk SGI ST336753LC 36 GB
#    Ethernet controller: NetXtreme BCM5701 Gigabit Ethernet (rev 15) on pci01.04.0
# Ethernet controller: Broadcom Corporation NetXtreme BCM5704 Gigabit Ethernet (rev 03) on pci04.01.0
# Ethernet controller: Broadcom Corporation NetXtreme BCM5704 Gigabit Ethernet (rev 03) on pci04.01.1
# SCSI Controller: QLogic 12160 Dual Channel Ultra3 SCSI (Rev 06) on pci03.01.0
#   Disk Drive: unit   1 lun  0 on SCSI controller pci03.01.0-1  0
# <cut>
# SCSI Controller: Fibre Channel QLA2342 (Rev 02) on pci05.01 port 1
#   Fabric Disk: node 20000080e511a226 port 3 lun 0   on SCSI controller pci05.01.0  0
# <cut>
# SCSI Controller: Fibre Channel QLA2342 (Rev 02) on pci05.01 port 2
#   Fabric Disk: node 20000080e511a226 port 3 lun 0   on SCSI controller pci05.01.1  0
# <cut>
# 
# ifort installed in:
# /opt/modules/cmplrs/intel/8.1.023/bin/ifort
# 
# ifort -V:
# Intel(R) Fortran Itanium(R) Compiler for Itanium(R)-based applications
# Version 8.1    Build 20041123 Package ID: l_fc_pc_8.1.023
# Copyright (C) 1985-2004 Intel Corporation.  All rights reserved.
# 
# icc installed in:
# /opt/modules/cmplrs/intel/8.1.023/bin/icc
# 
# icc -V:
# Intel(R) C++ Itanium(R) Compiler for Itanium(R)-based applications
# Version 8.1    Build 20041021 Package ID: l_cc_pu_8.1.024
# Copyright (C) 1985-2004 Intel Corporation.  All rights reserved.
#
