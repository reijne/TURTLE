#doc  This is the machine-specific file for the serial version of GAMESS-UK
#doc  using the Intel compiler 15.0 running on x86-64 processors.
#doc
#doc  Options:
#doc  gcc  - option to compile all C-code with gcc (if you have no Intel C-compiler).
#doc  goto - link against the GOTO blas library 
#doc  vdw  - include vanderwaals terms
#
#dopt vb nbo mrdci zora drf vdw sysmo dl-find i8 
#opt gcc goto debug qmmm newscf mopac serial static 
#
# ================ M4 Processing options
#
MACHINE_KEY=G
MACHOPT=linux,pclinux,littleendian,cio,unix,upck-equiv,xeon,helfey

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
# this causes error in ifort 11.1
# LDFLAGSTMP  = -Wl,-Map load.map
LDFLAGSTMP  = 


#--#if debug#
FFLAGSV  = ${FFLAGSTMP} -g
FFLAGSS  = ${FFLAGSTMP} -g
FFLAGSN  = ${FFLAGSTMP} -g
FFLAGSN0 = ${FFLAGSTMP} -g
CFLAGS   = -g ${CFLAGSTMP}
#--#else# 
FFLAGSV  = ${FFLAGSTMP} -ip -O3  -qopt-report-phase=vec
FFLAGSS  = ${FFLAGSTMP} -ip -O2  -qopt-report-phase=vec
FFLAGSN  = ${FFLAGSTMP} -ip -O1 
FFLAGSN0 = ${FFLAGSTMP} -O0 
CFLAGS   = ${CFLAGSTMP} 
#--#endif#

#--#if static#
LDFLAGS=${LDFLAGSTMP} -static -g
#--#else#
LDFLAGS=${LDFLAGSTMP} -g
#--#endif static#
#


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

BL_LIB = ${LBLAS} ${lBLAS} 

#
# ===============  Additional Files 
#
EXTRA=gethes.o helfey.o sp0011.o
EXTRA_MP2= mkmakw.o mpmakw.o umpe3a.o umpe3b.o check0a.o bmove.o chfeq.o
EXTRA_DFT = jkint_dft.o jkder_dft.o

#
# ===============  Compiler Exceptions
#
gethes.o:	casb.m
	cat  ../utilities/gener.m casb.m | $(M4) -DGEN_EXTRACTFILE=gethes $(M4OPTS) > gethes.f
	$(FC) $(FFLAGSN) gethes.f

helfey.o:     drv1e.m
	cat ../utilities/gener.m drv1e.m | $(M4) -DGEN_EXTRACTFILE=helfey $(M4OPTS) > helfey.f
	$(FC) $(FFLAGSV) helfey.f

sp0011.o:     integs.m
	cat ../utilities/gener.m integs.m | $(M4) -DGEN_EXTRACTFILE=sp0011 $(M4OPTS) > sp0011.f
	$(FC) $(FFLAGSN) sp0011.f

mkmakw.o:	drvmp.m
	cat  ../utilities/gener.m drvmp.m | $(M4) -DGEN_EXTRACTFILE=mkmakw $(M4OPTS) > mkmakw.f
	$(FC) $(FFLAGSS) mkmakw.f

chfeq.o:	cphf.m
	cat ../utilities/gener.m cphf.m | $(M4) -DGEN_EXTRACTFILE=chfeq $(M4OPTS) > chfeq.f
	$(FC) $(FFLAGSN) chfeq.f

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

bmove.o:	mclr.m
	cat ../utilities/gener.m mclr.m | $(M4) -DGEN_EXTRACTFILE=bmove $(M4OPTS) > bmove.f
	$(FC) $(FFLAGSN) bmove.f

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
# Linux raven13 2.6.32-504.el6.x86_64 #1 SMP Tue Sep 16 01:56:35 EDT 2014 x86_64 x86_64 x86_64 GNU/Linux
# 
# cat /proc/cpuinfo:
# processor       : 0
# vendor_id       : GenuineIntel
# cpu family      : 6
# model           : 45
# model name      : Intel(R) Xeon(R) CPU E5-2650 0 @ 2.00GHz
# stepping        : 7
# microcode       : 1803
# cpu MHz         : 2000.000
# cache size      : 20480 KB
# physical id     : 0
# siblings        : 8
# core id         : 0
# cpu cores       : 8
# apicid          : 0
# initial apicid  : 0
# fpu             : yes
# fpu_exception   : yes
# cpuid level     : 13
# wp              : yes
# flags           : fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc arch_perfmon pebs bts rep_good xtopology nonstop_tsc aperfmperf pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 cx16 xtpr pdcm pcid dca sse4_1 sse4_2 x2apic popcnt tsc_deadline_timer aes xsave avx lahf_lm ida arat epb xsaveopt pln pts dts tpr_shadow vnmi flexpriority ept vpid
# bogomips        : 4000.32
# clflush size    : 64
# cache_alignment : 64
# address sizes   : 46 bits physical, 48 bits virtual
# power management:
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
