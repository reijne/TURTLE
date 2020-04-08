##################################
# This is a gnu make (gmake)     #
# Makefile. You will need        #
# to have gmake on your machine! #
#                                #
#     see last lines!            #
#                                #
#                                #
#               ghost. 3.98      #
#                                #
##################################

-include ../../../../../Makefile.conf

LINKER		= $(FCOMPILER)
LINTFLAGS	=
MAKEFILE	= Makefile
PRINT		= pr
PROGRAM		= ./ef.x
ISOVIB		= ./isovib     
PSE	        = ./pse     
EXTHDRS		= 
HDRS		= 
EFOBJS= reada.o readl.o readmo.o main.o second.o xyzint.o gmetry.o wrttxt.o nuchar.o getgeo.o \
      getsym.o symtry.o geout.o haddon.o depvar.o getgeg.o gettxt.o axis.o rsp.o \
      block.o out.o compfg.o compfgc.o egrad.o egradint.o matpak.o jcarin.o dot.o  \
      optout.o optoutint.o hqrii.o hess.o hessc.o prmat.o hessout.o perm.o irc.o esthess.o \
      cout.o iout.o readdip.o calcint.o transformmu.o ef.o

ISOVIBOBJS= reada.o readl.o readmo.o second.o xyzint.o gmetry.o wrttxt.o nuchar.o getgeo.o \
      getsym.o symtry.o geout.o haddon.o depvar.o getgeg.o gettxt.o axis.o rsp.o \
      block.o out.o compfg.o compfgc.o egrad.o egradint.o matpak.o jcarin.o dot.o  \
      optout.o hqrii.o hess.o hessc.o prmat.o hessout.o perm.o irc.o esthess.o \
      cout.o  optoutint.o iout.o readdip.o calcint.o transformmu.o isovib.o

PSEOBJS= pse.o block.o
      

all:		$(PROGRAM) $(ISOVIB) $(PSE)

$(PROGRAM):     $(EFOBJS) $(SPECIAL) $(LIBS)
		@echo -n "Loading $(PROGRAM) ... "
		@$(LINKER) $(LDFLAGS) $(EFOBJS) $(SPECIAL) $(LIBS) -o $(PROGRAM)
		@echo "ef.x done"

$(ISOVIB):      $(ISOVIBOBJS) $(SPECIAL) $(LIBS)
		@echo -n "Loading $(ISOVIB) ... "
		@$(LINKER) $(LDFLAGS) $(ISOVIBOBJS) $(SPECIAL) $(LIBS) -o $(ISOVIB)
		@echo "isovib done"

$(PSE):         $(PSEOBJS) $(SPECIAL) $(LIBS)
		@echo -n "Loading $(PSE) ... "
		@$(LINKER) $(LDFLAGS) $(PSEOBJS) $(SPECIAL) $(LIBS) -o $(PSE)
		@echo "pse done"

depend:;	@mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

dist:;		@rdist

index:;		@ctags -wx $(HDRS) $(SRCS)

clean:;		rm -f $(PROGRAM) $(PSE) $(ISOVIB) $(EFOBJS) $(PSEOBJS) $(SPECIAL) isovib.o

lint:;		lint ${LINTFLAGS} ${SRCS}

print:;		@$(PRINT) $(HDRS) $(SRCS)

program:        $(PROGRAM)


tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)


%.o	:%.f
	$(FCOMPILER) $(FFLAGS) -c $<

###
