#doc  Machine-dependent makefile settings for MPI Linux/ build on
#doc  the SGI Altix (Itanium2) with Intel compilers version 8.1
#doc  MPI is expected to be found in /opt/modules/mpt-ia64/sgi-mpt-1.11/
#doc
#doc  Options:
#doc  mkl - use Intel MKL numerical library (default location: /opt/modules/intel/mkl721/lib/64)
#doc  scs - use the SGI scs numerical library instead of mkl
#doc  newscf - build the new distributed data SCF/DFT driver (requires BLACS and ScaLAPACK)
#doc  taskfarm   - build the GAMESS-UK taskfarming binary
#
# DEFAULT AND AVAILABLE OPTIONS
#dopt base newscf mkl
#opt debug scs taskfarm vdw
#
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,upck-equiv,itanium,altix

#MPI flags for efc compiler
MPI_INCLUDE = /opt/modules/mpt-ia64/sgi-mpt-1.11/include
LMPI = -L/opt/modules/mpt-ia64/sgi-mpt-1.11/lib
lMPI= -lmpi

#
# ===============  Compiler Options
#
FC = ifort
LD = ifort
FC90 = ifort
CC= icc
FFLAGSTMP = -c -ftz -ip -I${MPI_INCLUDE}
#--#if debug#
FFLAGSV = ${FFLAGSTMP} -g
FFLAGSS = ${FFLAGSTMP} -g
FFLAGSN = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS = -g -c -DLINUX 
LDFLAGS  = -g
CXXFLAGS = -g -c
#--#else# 
FFLAGSV = ${FFLAGSTMP} -O2
FFLAGSS = ${FFLAGSTMP} -O
FFLAGSN = ${FFLAGSTMP} -O1
FFLAGSN0 = ${FFLAGSTMP} 
CFLAGS  = -c -DLINUX
LDFLAGS =  
CXXFLAGS = -O2 -c
#--#endif#


OPTIONS=${MACHOPT}

#--#if mkl#
LIBBLAS=-L/opt/modules/intel/mkl721/lib/64 -lmkl
#--#elseif scs#
LIBBLAS=-lscs
#--#endif#

#--#if newscf#
# Blacs and ScaLAPACK
MPI90_LIB=-lsdsm
#--#endif#

BL_LIB= $(LLIB) $(lLIB) $(LMPI) $(lMPI) $(LIBBLAS) $(MPI90_LIB) 

# Diesel build options
LD_DIESEL = icpc
#--#if mkl scs#
DIESEL_LIBS = ${LIBBLAS} -lifcore
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
EXTRA=gethes.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o tst1s.o


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
tst1s.o:	direct.m
	cat ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f 

# Information on the machine this build was configured on
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
