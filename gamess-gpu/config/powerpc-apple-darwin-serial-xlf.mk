#doc Makefile settings for serial PowerPC ppc970/MacOSX Darwin Version 7.9.0
#doc build with xlf compiler version 8.1
#doc Options:
#doc blas   - use the framework veclib
#doc tiger  - build for MacOSX Tiger
#doc static - build a statically linked binary
#
# DEFAULT AND OPTIONAL OPTIONS
#dopt zora vb drf mrdci nbo mopac vdw sysmo dl-find blas 
#opt debug static tiger diesel
#
#
# ================ M4 Processing options
#
MACHINE_KEY=r
MACHOPT=rs6000,cio,unix,doublebackslash,upck-equiv,macosx,xlf


IOR64 = ior($$1,$$2)
IXOR64 = ieor($$1,$$2)
EXTRA_M4_DEFINITIONS = -Dleadz='leadzz($$1)'

## Compiler Options ##

#--#if static#
# Links to static libraries must be made and xlf.cfg extended
# see ABSOFT mailinglist mac
FC = xlf -F:xlf_s
LD = xlf -F:xlf_s
#--#else#
FC = xlf
LD = xlf
#--#endif#

FC90 = xlf90
CC = cc
CPP = /usr/lib/cpp
RANLIB = ar ts
FFLAGSTMP = -c 
FFLAGSV90 = -c -qsuffix=f=f90
FFLAGSS90 = -c -qsuffix=f=f90
FFLAGSN90 = -c -qsuffix=f=f90
CXX = c++

#--#if debug#
FFLAGSV = ${FFLAGSTMP} -qEXTNAME ${FFLAGSI8} -g
FFLAGSS = ${FFLAGSTMP} -qEXTNAME ${FFLAGSI8} -g
FFLAGSN = ${FFLAGSTMP} -qEXTNAME ${FFLAGSI8} -g
FFLAGS0 = ${FFLAGSTMP} -qEXTNAME ${FFLAGSI8} -g
FFLAGS1 = ${FFLAGSTMP} -qEXTNAME ${FFLAGSI8} -g
CFLAGS = -g ${FFLAGSI8} -c 
CXXFLAGS = -g -c
#--#else# 
FFLAGSV = -O  -c  ${FFLAGSI8} -qEXTNAME -qmaxmem=-1
FFLAGSS=  -O  -c  ${FFLAGSI8} -qEXTNAME -qmaxmem=-1
FFLAGSN = -g  -c  ${FFLAGSI8} -qEXTNAME -qmaxmem=-1
FFLAGS0 = -c   ${FFLAGSI8} -qEXTNAME -qmaxmem=-1
FFLAGS1 = -c   ${FFLAGSI8} -qEXTNAME -qmaxmem=-1
CFLAGS = ${CFLAGSI8}  -O  -c
CXXFLAGS = -O3 -c
#--#endif#

#--#if blas#
LDFLAGS =  -qEXTNAME -framework veclib
#---#if tiger#
LDFLAGS =  -qEXTNAME -framework veclib -lSystemStubs
#---#endif#
#--#else#
LDFLAGS  = -qEXTNAME
#---#if tiger#
LDFLAGS =  -qEXTNAME lSystemStubs
#---#endif#
#--#endif#


OPTIONS=${MACHOPT}

# Diesel build options
LD_DIESEL = ${CXX}
#--#if blas#
DIESEL_LIBS = -framework veclib -lstdc++ -L/opt/ibmcmp/xlf/8.1/lib -lxlf90 -lxlfmath -lxl
#--#else#
DIESEL_LIBS = -lstdc++ -L/opt/ibmcmp/xlf/8.1/lib -lxlf90 -lxlfmath -lxl
#--#endif#
FLEX = /usr/bin/flex
__UNDERBAR = 0
SIZEOF_VOID_P = 4
LONG_LONG_INT = long long int
LONG_INT = long int
INT = int
SHORT_INT = short int

#
# ===============  Additional Files 
#

EXTRA_BASE= iterate.o
#--#if mrdci#
EXTRA=skiny.o foxy.o outp.o
#--#endif#

#
# ===============  Compiler Exceptions
#
#  extra exceptions for macosx (xlf)
#
ccsd.o: ccsd.m
	cat ../utilities/gener.m  ccsd.m | $(M4)  $(M4OPTS)  > ccsd.f
	$(FC) $(FFLAGS0) ccsd.f
nvccsd.o: nvccsd.m
	cat ../utilities/gener.m  nvccsd.m | $(M4)  $(M4OPTS)  > nvccsd.f
	$(FC) $(FFLAGS0) nvccsd.f
#
#
#  ========== Exceptions for IBM r6000 power1, 2 and 3 (and SP2) ===============
#

iterate.o: morokuma.m
	cat ../utilities/gener.m  morokuma.m | $(M4) -DGEN_EXTRACTFILE=iterate $(M4OPTS)  > iterate.f
	$(FC) $(FFLAGS0) iterate.f

skiny.o:	newmrd5.m
	cat ../utilities/gener.m  newmrd5.m | $(M4) -DGEN_EXTRACTFILE=skiny $(M4OPTS) $(SNGL) > skiny.f
	$(FC) $(FFLAGS1) skiny.f
foxy.o:	newmrd5.m
	cat ../utilities/gener.m  newmrd5.m | $(M4) -DGEN_EXTRACTFILE=foxy $(M4OPTS) $(SNGL) > foxy.f
	$(FC) $(FFLAGS1) foxy.f
outp.o:	newmrd5.m
	cat ../utilities/gener.m  newmrd5.m | $(M4) -DGEN_EXTRACTFILE=outp $(M4OPTS) $(SNGL) > outp.f
	$(FC) $(FFLAGS1) outp.f             

# Build Information
#
# Output of uname -a :
# Darwin tc7 7.9.0 Darwin Kernel Version 7.9.0: Wed Mar 30 20:11:17 PST 2005; root:xnu/xnu-517.12.7.obj~1/RELEASE_PPC  Power Macintosh powerpc
# 
# Output of hostinfo:
# hostinfo
# Mach kernel version:
#          Darwin Kernel Version 7.9.0:
# Wed Mar 30 20:11:17 PST 2005; root:xnu/xnu-517.12.7.obj~1/RELEASE_PPC
#  
# Kernel configured for up to 2 processors.
# 2 processors are physically available.
# Processor type: ppc970 (PowerPC 970)
# Processors active: 0 1
# Primary memory available: 2048.00 megabytes.
# Default processor set: 85 tasks, 199 threads, 2 processors
# Load average: 1.52, Mach factor: 0.59
# 
# Curtailed output of: system_profiler -detailLevel -1
# Hardware:
#     Hardware Overview:
#       Machine Model: Power Mac G5
#       CPU Type: PowerPC 970  (2.2)
#       Number Of CPUs: 2
#       CPU Speed: 2 GHz
#       L2 Cache (per CPU): 512 KB
#       Memory: 2.5 GB
#       Bus Speed: 1 GHz
#       Boot ROM Version: 5.1.4f0
#       Serial Number: CK3460CYQE8
#     
# Software:
#     System Software Overview:
#       System Version: Mac OS X 10.3.9 (7W98)
#       Kernel Version: Darwin 7.9.0
#       Boot Volume: G5
# 
# Network:
#     Internal Modem:
#       Interface: modem
#       Type: PPP (PPPSerial)
#     Built-in Ethernet:
#       Interface: en0
#       Type: Ethernet
#     VPN (PPTP):
#       Type: PPP (PPTP)
# 
# Memory:
#     DIMM0/J11:
#       Size: 256 MB
#       Type: DDR SDRAM
#       Speed: PC3200U-30330
#     DIMM1/J12:
#       Size: 256 MB
#       Type: DDR SDRAM
#       Speed: PC3200U-30330
#     DIMM2/J13:
#       Size: 512 MB
#       Type: DDR SDRAM
#       Speed: PC3200U-30330
#     DIMM3/J14:
#       Size: 512 MB
#       Type: DDR SDRAM
#       Speed: PC3200U-30330
#     DIMM4/J41:
#       Size: 512 MB
#       Type: DDR SDRAM
#       Speed: PC3200U-30330
#     DIMM5/J42:
#       Size: 512 MB
#       Type: DDR SDRAM
#       Speed: PC3200U-30330
#     DIMM6/J43:
#       Size: Empty
#       Type: Empty
#       Speed: Empty
#     DIMM7/J44:
#       Size: Empty
#       Type: Empty
#       Speed: Empty
# 
# Information on: /opt/ibmcmp/xlf/8.1/bin/xlf
# IBM XL Fortran Advanced Edition Version 8.1 for Mac OS X
# C) Copyright IBM Corp. 1990, 2003.   All Rights Reserved.
