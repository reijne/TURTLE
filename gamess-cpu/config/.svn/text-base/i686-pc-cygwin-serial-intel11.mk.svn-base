#doc This is the machine-specific file for the serial Windows build under Cygwin
#doc on a Xeon processor with the Intel compiler version 11.0
#doc
#doc The build environment that this file was created on was not ideal so paths
#doc have had to be hard-wired in. The compilers are expected to be found in:
#doc C:\Program Files (x86)\Intel\Compiler\11.0\074
#doc
#doc The M$ Visual Stuiod installation is expected to be found in:
#doc C:\Program Files (x86)\Microsoft Visual Studio 9.0
#doc 
#doc If the MKL libraries are to be used they are expected to be found in:
#doc E:\Program Files\Intel\MKL60\ia32\lib
#doc 
#doc Currently we need to run:
#doc export PATH=/cygdrive/c/Program\ Files\ \(x86\)/Microsoft\ Visual\ Studio\ 9.0/Common7/IDE:/cygdrive/c/Program\ Files\ \(x86\)/Microsoft\ Visual\ Studio\ 9.0/VC/BIN:$PATH
#doc
#doc before building to get the c-compilation to work.
#doc See the notes at the end of this file for a further explanation
# 
#dopt mrdci zora vb drf nbo vdw sysmo mopac dl-find newscf
#opt mkl

# Rename library extensions and add the flag so we can use .o and not .obj
OBJNAME=-Fo:$*.o

# LD doesn't accept -o flag so...
LDNAME=-Fe

# Machine-specific 
MACHINE_KEY=G
IAND32=iand($$1,$$2)
IOR32=ior($$1,$$2)
IXOR32=ieor($$1,$$2)
SHIFT=ishift($$1,$$2)

# M4 Options
MACHOPT=linux,win95,ifc,littleendian,cio,unix,upck-equiv,i8drct

# Compiler options
FC=/cygdrive/c/Program\ Files\ \(x86\)/Intel/Compiler/11.0/074/fortran/Bin/IA32/ifort.exe
FC90=/cygdrive/c/Program\ Files\ \(x86\)/Intel/Compiler/11.0/074/fortran/Bin/IA32/ifort.exe
CC=/cygdrive/c/Program\ Files\ \(x86\)/Intel/Compiler/11.0/074/cpp/bin/ia32/icl.exe

# -Qlocation,link - says where to find the linker
LD=/cygdrive/c/Program\ Files\ \(x86\)/Intel/Compiler/11.0/074/fortran/Bin/IA32/ifort.exe -Qlocation,link,"C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC/bin"

FFLAGSBASE = -c -w90 -w95 -Qvec-report0
FFLAGSV  = $(FFLAGSBASE) -Qip -O3
FFLAGSS  = $(FFLAGSBASE) -Qip -O2
FFLAGSN  = $(FFLAGSBASE) -Qip -O1
FFLAGSN0 = $(FFLAGSBASE) -O1
CFLAGS   = -c -Qip -O3 -I"C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\include"

#--#if dl-find#
# Need to pass the -fpp flag to ifort
# See dl-find/Makefile.in.gamess
DLF-PPFLAGS= -fpp -DGAMESS -DOLDALLOC

# Need to set or else we can't make the dl-find library
GMAKE=/usr/bin/make
#--#endif dl-find#

LDBASE=-Qip -Qoption,link,/LIBPATH:"C:\Program Files (x86)\Intel\Compiler\11.0\074\cpp\Lib\ia32",/LIBPATH:"C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\LIB",/LIBPATH:"C:\Program Files\Microsoft SDKs\Windows\v6.0A\lib",/LIBPATH:"C:\Program Files (x86)\Intel\Compiler\11.0\074\fortran\lib\ia32"

#--#if mkl#
LDBLAS = -Qoption,link,-LIBPATH:"C:\Program Files (x86)\Intel\Compiler\11.0\074\fortran\mkl\ia32\lib"
LIBBLAS =  mkl_c_dll.lib
BLASOPT=,blas
#--#endif mkl#

LDFLAGS= ${LDBASE} ${LDBLAS}

OPTIONS=${MACHOPT}$(BLASOPT)
BL_LIB = ${LIBBLAS}

EXTRA_BASE= bitmanip.o
EXTRA_BENCH= bitmanip.o
EXTRA_MP2= check0a.o

#
#
#  ================ Exceptions for linux (g77, pgf77, absoft, efc) =======
#
#
tst1s.o:	direct.m
	cat  ../utilities/gener.m direct.m | $(M4) -DGEN_EXTRACTFILE=tst1s $(M4OPTS) > tst1s.f
	$(FC) $(FFLAGSN) -Fotst1s.o tst1s.f 
#
rpa.o:	rpa.f
	$(FC) $(FFLAGSN) $(OBJNAME) $*.f
	
mrdci5.o:	mrdci5.f
	$(FC) $(FFLAGSN) $(OBJNAME) $*.f

check0a.o:	newmrd5.m
	cat ../utilities/gener.m newmrd5.m | $(M4) -DGEN_EXTRACTFILE=check0a $(M4OPTS) > check0a.f
	$(FC) $(FFLAGSN) -object:check0a.o check0a.f

# With O1 c2021_j-q and c2031 fail to validate
# NB - only tested with FFLAGSN0 - might work with higher optimisation.
newmrd1.o: newmrd1.f
	$(FC) $(FFLAGSN0) $(OBJNAME) $*.f

#
#  ====================== Exceptions for Windows95 =======================
#
bitmanip.f:	machscf.m
	cat ../utilities/gener.m machscf.m | $(M4) -DGEN_EXTRACTFILE=bitmanip \
	$(M4OPTS) $(SNGL) | ../utilities/quote > bitmanip.f
	

# NOTES
#
# This build was configured at Daresbury Laborato1ry on Igor Kozin's csehtms.csems.dl.ac.uk machine
# 
# Unfortunately, it seems that the Intel compilers are not correctly integrated with Cygwin so a lot of this stuff is about
# fixing that and probably wouldn't be required on other machines.
# 
# Running the script:
# 
# "C:\Program Files (x86)\Intel\Compiler\11.0\074\cpp\bin\iclvars.bat" ia32
# 
# Fixed the problems in the DOS shell, but weren't propagaged to the Cygwin shell
# 
# The initial problem was with getting simple c-programs to compile.
# 
# Calling icl directly with the include path set to:
# 
# /I"C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\INCLUDE"
# 
# Caused a stream of errors with the stdio.h include.
# 
# Using the /showIncludes flag to ICL showed that it was only including files from:
# "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\INCLUDE" directory.
# 
# It turned out it's setting the PATH as below that fixes the problem:
# export PATH=/cygdrive/c/Program\ Files\ \(x86\)/Microsoft\ Visual\ Studio\ 9.0/Common7/IDE:\
# /cygdrive/c/Program\ Files\ \(x86\)/Microsoft\ Visual\ Studio\ 9.0/VC/BIN:$PATH
# 
# The link step then failed as several dll's couldn't be found.
# 
# These are described below, together with their location:
# 
# kernel32.lib
# set WindowsSdkDir=C:\Program Files\Microsoft SDKs\Windows\v6.0A\
# set LIB=%WindowsSdkDir%\lib
# 
# libcmt.lib
# C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\ATLMFC\LIB;C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\LIB
# 
# libmmt.lib
# SET ICPP_COMPILER11=C:\Program Files (x86)\Intel\Compiler\11.0\074\cpp\Lib\ia32
# SET LIB=%ICPP_COMPILER11%\Lib\ia32
# 

