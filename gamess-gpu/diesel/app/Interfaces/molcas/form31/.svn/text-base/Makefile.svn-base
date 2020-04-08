#


F77 = g77  
FFLAGS  = -O
#FFLAGS  = -K -g
#FFLAGS  = 
LOADER  = g77        
LOPT    = #-lblas
libs    = $(MOLCAS)/lib/g77/libmolcas.a 
#libs    = $(MOLCAS)/lib/g77/libaixio.a \
#          $(MOLCAS)/lib/g77/libutil1.a \
#          $(MOLCAS)/lib/g77/libutil2.a \
#          $(MOLCAS)/lib/g77/libutil3.a \
#          $(MOLCAS)/lib/g77/libessl.a \
#          $(LOPT)


#################################################
FORM31OBS=\
form31.o upcase.o dumper.o kotz.o csort.o mknit.o \
 geth0.o recprt1.o recprt2.o rs4ind.o rdcmo.o rdvec.o blkdata.o

.f.o:
	$(F77) $(FFLAGS) -c $<

form31.x: $(FORM31OBS)
	$(F77) -o form31.x $(FORM31OBS) $(libs)
#################################################
