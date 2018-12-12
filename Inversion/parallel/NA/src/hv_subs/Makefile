#
#	Make discription file for Receiver function example subroutines
#

#FFLAGS = -O4
#F77 = /opt/SUNWspro/bin/f77

DEMO = ../..

include $(DEMO)/src/makena.macros   # Edit this file to change location or
                        # type of compilers (e.g. to gcc, g77)

DIR = .
HV = $(DIR)/read_seed.o \
      $(DIR)/forward_modelling.o $(DIR)/calcmisfit.o \
      $(DIR)/lofw.o  \
      $(DIR)/read_vmodelNA.o  \
      $(DIR)/ran3.o $(DIR)/indexx.o 

# make hv_subs 

hv_subs: $(HV) 
	#cd ./rfi_plot; make -f Makefile

.f.o:
	$(F77) $(FFLAGS) -c $*.f
.c.o:
	$(CC) $(CFLAGS) -c $*.c
