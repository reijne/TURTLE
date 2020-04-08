# Include file for makefiles in individual subdirectories 


# LIBDIR defines where the libraries should be built.
# LIB_DISTRIB defines default location for the GA distribution package
# when built as a standalone package.
#

ifndef LIBDIR
   ifdef LIB_DISTRIB 
     LIBDIR = $(LIB_DISTRIB)/$(TARGET)
   else
     LIBDIR = .
   endif
endif

ifdef GA_SHLIB
	SHLIB_DIR    = $(LIBDIR)/shared
	OBJ_DIR      = ./obj
endif
# makefile in each directory might define error message for undefined 
# symbols etc.
# When error message is defined, it should be displayed and then
# make processing aborted
#

define print_error
        @echo $(ERRMSG)
        exit 1
        @echo
endef

ifdef ERRMSG
error:
	$(print_error)
endif

# Make sure that nothing gets compiled in case of error 
#
ifdef ERRMSG
      CC = $(print_error)
      FC = $(print_error)
endif

ifdef LIBRARY_STAMP
      SYMBOL_STAMP = $(TARGET).$(LIBRARY_STAMP)
else
      SYMBOL_STAMP = $(TARGET)
endif

FULL_LIBRARY_PATH := $(LIBDIR)/$(LIBRARY)
FULL_LIBRARY_PATH_SHARED := $(SHLIB_DIR)/$(LIBRARY_SHARED)

OBJECTS := $(OBJ) $(OBJ_FRAGILE)
SHLIB_OBJECTS := $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

 LIBOBJ := $(patsubst %,$(FULL_LIBRARY_PATH)(%),$(OBJECTS))
 LIBOBJ_NOPT := $(patsubst %,$(FULL_LIBRARY_PATH)(%),$(OBJ_FRAGILE))

$(FULL_LIBRARY_PATH): $(LIBDIR) $(SYMBOL_STAMP).stamp  $(LIBOBJ) $(LIBOBJ_NOPT) $(FULL_LIBRARY_PATH_SHARED)
	$(RANLIB) $@

$(FULL_LIBRARY_PATH_SHARED):  $(SHLIB_DIR) $(OBJ_DIR) $(LIBOBJ) $(LIBOBJ_NOPT)
ifdef GA_SHLIB
	@for i in $(OBJECTS); do \
		if [ -f $$i ]; then \
			cp $$i $(OBJ_DIR); \
		fi; \
	done
	$(CC) $(SHLIB_LDFLAGS) -o $(SHLIB_DIR)/$(LIBRARY_SHARED) $(SHLIB_OBJECTS);
endif

$(SYMBOL_STAMP).stamp:
	if [ -f *.stamp ]; then\
		 $(MAKE) cleanstamp;\
	fi
	echo "" > $(SYMBOL_STAMP).stamp

ifndef GA_SHLIB
ifdef OBJ_FRAGILE
ifndef FRAGILE
.PHONY: $(LIBOBJ_NOPT)
$(LIBOBJ_NOPT):
	@$(MAKE) FRAGILE="Yes"
endif
endif
endif


MAKESUBDIRS = for dir in $(SUBDIRS); do $(MAKE)  -C $$dir $@ || exit 1 ; done


ifdef SUBDIRS

$(FULL_LIBRARY_PATH):        subdirs

.PHONY: subdirs
subdirs:
	@for dir in $(SUBDIRS); do \
		cd $$dir ; \
		$(MAKE) || exit 1  ; \
		cd .. ;\
	done
endif
#		$(MAKE)  -C $$dir || exit 1 ;  \

.PHONY: clean
clean:
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	-$(RM) -f *.o *.p *core *stamp mputil.mp* *trace *.x *.exe obj/* *events* $(LIB_TARGETS)
	-$(RM) -rf ./obj
ifdef HARDCLEAN 
	-$(RM) -f $(FULL_LIBRARY_PATH)
else
	if [ -f $(FULL_LIBRARY_PATH) ] ; then \
		$(AR) d $(FULL_LIBRARY_PATH) $(OBJ) $(OBJ_FRAGILE) ; \
		if [ `$(AR) t $(FULL_LIBRARY_PATH) | wc | awk ' {print $$1;}'` -eq 0 ] ; then \
			$(RM) -f $(FULL_LIBRARY_PATH) ; \
		fi ; \
	fi ;
endif
ifdef GA_SHLIB
	@if [ -f $(SHLIB_DIR)/$(LIBRARY_SHARED) ]; then\
		$(RM) -f $(SHLIB_DIR)/$(LIBRARY_SHARED); \
		$(RM) -rf $(OBJ_DIR); \
	fi;
endif


.PHONY: realclean
realclean:      clean
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	-$(RM) -rf *.o?* *~ \#*\#

.PHONY: cleanstamp
cleanstamp: clean
ifdef SUBDIRS
	$(MAKESUBDIRS)
endif
	-$(RM) -rf *.stamp

$(LIBDIR):
	$(MKDIR) -p $@

ifdef GA_SHLIB
$(SHLIB_DIR):
	$(MKDIR) -p $@
$(OBJ_DIR):
	$(MKDIR) -p $@
endif

