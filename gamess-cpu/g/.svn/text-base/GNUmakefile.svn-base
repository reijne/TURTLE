# $Id: GNUmakefile,v 1.13 2007-10-30 11:48:30 jmht Exp $
#
# This makefile has to be processed by GNU make
#
EXPERIMENTAL=SGI64 LAPI64 CYGWIN HITACHI CATAMOUNT
SUPPORTED_TARGETS = SOLARIS SOLARIS64 SGI SGI_N32 SGITFP HPUX HPUX64 IBM IBM64 CRAY-SV1 cray-sv2
SUPPORTED_TARGETS+= LAPI DECOSF LINUX LINUX64 CRAY-T3E CRAY-YMP CYGNUS CYGWIN INTERIX
SUPPORTED_TARGETS+= FUJITSU-VPP FUJITSU-VPP64 NEC MACX BGL $(EXPERIMENTAL)
ALL_SUPPORTED_TARGETS = $(SUPPORTED_TARGETS) WIN32
ifndef TARGET
 TARGET=UNDEFINED
endif

# Makesure to avoid periods (i.e. dot) in GA_VERSION as that might(!) be problematic 
# in Windows
 GA_VERSION = 4-1-1B
###     TOPDIR = ga-$(GA_VERSION)
     TOPDIR = g
         RM = rm
      MKDIR = mkdir
MAKESUBDIRS = for dir in $(SUBDIRS); do $(MAKE)  -C $$dir $@ || exit 1 ; done
#    SUBDIRS = ma global tcgmsg-mpi LinAlg/lapack+blas tcgmsg armci/src pario ga++
    SUBDIRS = ma global tcgmsg-mpi tcgmsg armci/src pario ga++
  USE_ARMCI = yes

ifdef USE_MPI
     MP_VER = MPI
endif
ifeq ($(MSG_COMMS),MPI)
ifdef USE_MPI
error: 
	@echo Cannot define MSG_COMMS and USE_MPI together
	exit 1
endif
     MP_VER = MPI
endif

ifdef GA_C_CORE
ifdef USE_SCALAPACK
error:
	@echo Cannot define GA_C_CORE and SCALAPACK together
endif
ifeq ($(DIAG),PAR)
error:
	@echo Cannot assign DIAG=PAR, when GA_C_CORE is defined
endif
endif

ifndef MP_VER
  EXTRA_INC = "TCG_INCLUDE=../../include"
      EXTRA = "MSG_COMMS=TCGMSG"
     MP_VER = TCG
endif

ifdef GA_C_CORE
VER="-CC"
endif

STAMP = $(TARGET)-$(MP_VER)$(VER).stamp

MAKEDIRS = lib/$(TARGET)
ifdef GA_SHLIB
   SHLIB_DIR = lib/$(TARGET)/shared
   MAKEDIRS += $(SHLIB_DIR)
endif

#all: libs test
all: libs

libs:   info $(STAMP) $(MAKEDIRS) include
	@echo '>>> making Linear Algebra utility libraries <<< '
	(cd ./LinAlg/lapack+blas; $(MAKE))
	@echo '>>> making Memory Allocator (MA) library <<< '
	(cd ./ma; $(MAKE))
        ifdef USE_MPI
	      @echo '>>> making TCGMSG library on top of MPI <<< '
	      (cd ./tcgmsg-mpi; $(MAKE) MACHINE=$(TARGET) )
        else
           ifneq ($(MSG_COMMS),MPI)
	      @echo '>>> making TCGMSG library <<<'
	      (cd ./tcgmsg; $(MAKE) MACHINE=$(TARGET))
           endif
        endif
	@echo '>>> making ARMCI remote memory copy library <<< '
	(cd ./armci/src; \
	$(MAKE) LIB_DISTRIB=../../lib $(EXTRA) $(EXTRA_INC) )
	@echo '>>> making Global Arrays library <<< '
	(cd ./global; $(MAKE))
	(cd ./pario; $(MAKE))
#	(cd ./pario; $(MAKE) LARGE_FILES=yes)
#	(cd ./ga++; $(MAKE))

include ./armci/config/makecoms.h
include ./armci/config/makemp.h
include ./config/makefile.h
ifndef LIBDIR
   LIB_DISTRIB = $(PWD)/lib
endif

test:
	@echo ' ============================================================ '
	@echo ' >>> GA Libraries compiled. Making an executable test <<< '
	(cd ./global/testing; $(MAKE) test.x)
	@echo ' ------------------------------------------------------------ '
	@echo ' An executable test program for GA is ./global/testing/test.x'
	@echo ' There are also other test programs in that directory.'
	@echo ' ------------------------------------------------------------ '
	@echo ' Also, to test your GA programs, suggested compiler/linker '
	@echo ' options are as follows.'
	@echo ' GA libraries are built in' `pwd`/lib/$(TARGET)
	@echo ' INCLUDES = ' -I`pwd`/include
	@echo ' '
	@echo ' For Fortran Programs: '
	@echo '    FLAGS = $(NOPT) $(FLDOPT)'
	@echo '    LIBS  = $(LIBS) $(FLIBS)'
	@echo ' '
	@echo ' For C Programs: '
#	@echo '    FLAGS = $(NOPT) $(CLDOPT)'
	@echo '    LIBS  = $(LIBS) $(CLIBS)'
	@echo ' ============================================================ '

info:
	@echo ''
	@echo '*** This is the version $(GA_VERSION) of the Global Arrays package ****'
	@echo ''
ifneq ($(TARGET),$(findstring $(TARGET),$(SUPPORTED_TARGETS)))
	@echo "*** The following target platforms are supported: $(ALL_SUPPORTED_TARGETS)"
	@echo '*** where CRAY-YMP works for J90/SV1, LAPI works for IBM SP,'
	@echo '*** and LINUX also supports FreeBSD/Intel'
	@echo ''
	@echo "*** You must define TARGET platforms, for example"
	@echo "***     setenv TARGET IBM"
	@echo "*** or "
	@echo "***     make TARGET=IBM"
	@echo "*** See `pwd`/README file for details"
	@echo ''
	exit 1
endif

ma:	info $(STAMP) lib/$(TARGET) include
	@echo '*** Building Memory Allocator Library ****'
	(cd ./ma; $(MAKE))


lib/$(TARGET):
	$(MKDIR) -p $@

ifdef GA_SHLIB
$(SHLIB_DIR):
	$(MKDIR) -p $@
endif

include:
	$(MKDIR) -p $@

clean:
	$(MAKESUBDIRS) $^
	$(RM) -rf include/* *.stamp
	(cd ./armci/src; $(MAKE) LIB_DISTRIB=../../lib clean )

realclean: clean

$(STAMP):
	$(MAKE) clean
	echo "" > $(STAMP)

cleanstamp:
	$(MAKESUBDIRS) $^

stampit: 
	echo "" > $(STAMP)

HAS_TAR_EXCLUDE = SUN LINUX LINUX64 SOLARIS SOLARIS64 CYGNUS CYGWIN

# we might have GNU tar available on any platform (requires setting GNUTAR)
ifdef GNUTAR
  HAS_TAR_EXCLUDE += $(TARGET)
  TAR = gtar
else
  TAR = tar
endif


%.tar: 
ifeq ($(TARGET),$(findstring $(TARGET),$(HAS_TAR_EXCLUDE)))
	(cd ..;\
	 rm -f exclude.GA;\
	 echo '' > exclude.GA;\
	 find ./$(TOPDIR) -name "*" -type d -print |grep CVS | \
		sed  -e 's/\.\///' >> exclude.GA;\
	 find ./$(TOPDIR) -name "*.[oxpa]" -print|sed -e 's/\.\///' >>exclude.GA;\
	 find ./$(TOPDIR) -name "*core" -print|sed -e 's/\.\///' >>exclude.GA;\
	 find ./$(TOPDIR) -name "lib*" -print|sed -e 's/\.\///' >>exclude.GA;\
	 find ./$(TOPDIR) -name "include*" -print|sed -e 's/\.\///' >>exclude.GA;\
	 find ./$(TOPDIR) -name "*obj" -print|sed -e 's/\.\///' >>exclude.GA;\
	 find ./$(TOPDIR) -name "*exe" -print|sed -e 's/\.\///' >>exclude.GA;\
	 find ./$(TOPDIR) -name "*stamp" -print|sed -e 's/\.\///' >>exclude.GA;\
	 find ./$(TOPDIR) -name "include/?*" -print|sed -e 's/\.\///' >>exclude.GA;\
	 $(TAR) cvfX $@ exclude.GA $(TOPDIR))
else
	(cd ..; $(TAR) cvf $@ $(TOPDIR))
endif

.PHONY: cleanF
cleanF:
	@echo Not done here
source:

depend:
	@echo No dependency analysis here

