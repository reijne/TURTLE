#doc  This is the machine-specific file for the serial version of GAMESS-UK
#doc  using the Intel compiler as determined by module (now 11)
#doc
#doc  Options:
#doc  gcc  - option to compile all C-code with gcc (if you have no Intel C-compiler).
#doc  mkl - use the Intel MKL numerical library 
#doc  goto - link against the GOTO blas library 
#doc  vdw  - include vanderwaals terms
#
#dopt vb nbo mrdci zora drf vdw sysmo dl-find i8  mkl
#opt gcc goto debug qmmm newscf mopac
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,upck-equiv,opteron,helfey

# Bitwise operators that differ from the defaults.
IXOR32=ieor($$1,$$2)
IOR64=ior($$1,$$2)
IXOR64=ieor($$1,$$2)

#
# ===============  Compiler Options
#

## INTEL COMPILER##
FC = ifort
LD = ifort
FC90 = ifort
#--#if gcc#
CC= gcc
CFLAGSGCC= 
#--#else#
CC= icc
CFLAGSGCC= -no-gcc
#--#endif gcc#
RANLIB=ranlib

#--#if i8#
FFLAGSI8 = -i8
CFLAGSI8 = -DLONG_INTEGER
I8_M4_OPT=,i8
#--#else#
FFLAGSI8=
CFLAGSI8=
I8_M4_OPT=,64bitpointers
#--#endif i8#

#--#if ma#
# If we are building the memory allocator library we
# need to set some settings for the GAs
GA_F77_DEFS = -traditional
GA_TARGET=LINUX64
#--#endif ma#

FFLAGSTMP = -c -ftz ${FFLAGSI8}
CFLAGSTMP = -c ${CFLAGSGCC} ${CFLAGSI8} -DLINUX -D_FILE_OFFSET_BITS=64
LDFLAGSTMP  = -Wl,-Map load.map


#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g
FFLAGSS  = ${FFLAGSTMP} -g
FFLAGSN  = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS   = -g ${CFLAGSTMP}
LDFLAGS  = ${LDFLAGSTMP} -g
#--#else# 
FFLAGSV  = ${FFLAGSTMP} -ip -O3  -vec-report0
FFLAGSS  = ${FFLAGSTMP} -ip -O2  -vec-report0
FFLAGSN  = ${FFLAGSTMP} -ip -O1 
FFLAGSN0 = ${FFLAGSTMP} -O0 
CFLAGS   = ${CFLAGSTMP} 
LDFLAGS  = ${LDFLAGSTMP} -g
#--#endif#

#--#if mkl#
#---#if i4#
LIBBLAS =  -lmkl_ipf -lguide 
#---#else#
# The below is a shameful hack of which we should be embarrassed.
# If we are in i8 mode and linking to an i4 mkl library, then as 
# this is a littleendian architecture the correct bits get passed
# through when we pass an i8 integer to the i4 library. This trick
# only works for cases where we have single integer arguments.
# dgthr, dsctr and ddoti take integer array arguments and so
# we can't use the i4 mkl versions and need to select our own 
# i8-compiled ones.
#LIBBLAS = ../linalg/dsctr.o ../linalg/dgthr.o ../linalg/ddoti.o -lmkl_ipf -lguide 
#LIBBLAS= -lmkl_intel_ilp64 -lmkl_ipf -lmkl_core  -lmkl -lguide
LIBBLAS= -lmkl_intel_ilp64 -lmkl_intel_sp2dp -lmkl_core  -lmkl_def -lmkl_mc3 -lmkl -lguide
#---#endif#
BLASOPT=,blas
#--#endif#

#--#if goto#
LIBBLAS = -L/usr/local/lib64 -lblas
BLASOPT=,blas
### NOTE ####
# This version of GAMESS-UK is built with i8 by default. It is only fortuitous that linking
# against numerical libraries works with i8 (the significant bits are at the correct end).
# dgthr and dsctr take integer array arguments and so  the library versions crash when used
# with i8. So, when using blas with i8, we need to pull in our own versions of dsctr, dgthr
# & ddoti from m4/util3.m
#--#endif goto#

# Bring all the options together
OPTIONS=${MACHOPT}${I8_M4_OPT}${BLASOPT}

BL_LIB = ${LBLAS} ${LIBBLAS} 

#
# ===============  Additional Files 
#
EXTRA=gethes.o helfey.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o
EXTRA_DFT = jkint_dft.o jkder_dft.o

#
# ===============  Compiler Exceptions
#
gethes.o:	casb.m
	cat  ../utilities/gener.m casb.m | $(M4) -DGEN_EXTRACTFILE=gethes $(M4OPTS) > gethes.f
	$(FC) $(FFLAGSN) gethes.f

helfey.o:	drv1e.m
	cat ../utilities/gener.m drv1e.m | $(M4) -DGEN_EXTRACTFILE=helfey $(M4OPTS) > helfey.f
	$(FC) $(FFLAGSN) helfey.f

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

tst1s.o:	direct.m
	cat ../machines/$(MACH) ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) tst1s.f 

rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $*.f

eispack.o:	eispack.f
	$(FC) $(FFLAGSN) $*.f

mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $*.f

newmrd2.o:	newmrd2.f
	$(FC) $(FFLAGSN) $*.f
#
#  ========== DFT Exceptions (PGF) ===============
#
jkint_dft.o:    integ2e.m
	cat ../machines/$(MACH) ../utilities/gener.m  integ2e.m | $(M4) -DGEN_EXTRACTFILE=jkint_dft $(M4OPTS) $(SNGL) > jkint_dft.f
	$(FC) $(FFLAGSN) jkint_dft.f
#
jkder_dft.o:    deriv2e.m
	cat ../machines/$(MACH) ../utilities/gener.m  deriv2e.m | $(M4) -DGEN_EXTRACTFILE=jkder_dft $(M4OPTS) $(SNGL) > jkder_dft.f
	$(FC) $(FFLAGSN) jkder_dft.f 

# uname -a:
# Linux ladon 2.6.11.4-21.7-smp #1 SMP Thu Jun 2 14:23:14 UTC 2005 x86_64 x86_64 x86_64 GNU/Linux
# 
# cat /proc/cpuinfo:
# processor	: 0
# vendor_id	: AuthenticAMD
# cpu family	: 15
# model		: 5
# model name	: AMD Opteron(tm) Processor 246
# stepping	: 10
# cpu MHz		: 1994.377
# cache size	: 1024 KB
# physical id	: 255
# siblings	: 1
# fpu		: yes
# fpu_exception	: yes
# cpuid level	: 1
# wp		: yes
# flags		: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 pni syscall nx mmxext lm 3dnowext 3dnow
# bogomips	: 3923.96
# TLB size	: 1024 4K pages
# clflush size	: 64
# cache_alignment	: 64
# address sizes	: 40 bits physical, 48 bits virtual
# power management: ts fid vid ttp
# 
# mpif90 -compiler intel8 -V
# Intel(R) Fortran Compiler for Intel(R) EM64T-based applications, Version 9.0    Build 20050809 
# Copyright (C) 1985-2005 Intel Corporation.  All rights reserved.
# 
# mpif90 -compiler intel8 -show
# ln -s /opt/score/mpi/mpich-1.2.5/x86_64-suse-linux2_6_intel8/include/mpif.h mpif.h
# /opt/score/bin/scoref90 -compiler=intel8 -compiler=intel8 -L/opt/score/mpi/mpich-1.2.5/x86_64-suse-linux2_6_intel8/lib -lmpichf90 -lpmpich -lmpich -lpmpich -lmpich
# rm mpif.h
# 
# mpicc -compiler gnu --version
# gcc (GCC) 3.3.5 20050117 (prerelease) (SUSE Linux)
# Copyright (C) 2003 Free Software Foundation, Inc.
# 
# 
