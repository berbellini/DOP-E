#
#	Make discription file for NA
#

#FFLAGS = -O
#F77 = /opt/local/SUNWspro/bin/f77
DEMO = ../../../..

# Edit this file to change location or
# type of compilers (e.g. to gcc, g77)

include $(DEMO)/src/makena.macros

all: nad2rfi rfi2nad nad2nad nad2asc asc2nad

nad2rfi: nad2rfi.o
	$(F77) $(FFLAGS) nad2rfi.o -o $(DEMO)/bin/nad2rfi 

rfi2nad: rfi2nad.o
	$(F77) $(FFLAGS) rfi2nad.o -o $(DEMO)/bin/rfi2nad 
	cd $(DEMO)/data/demo_results;\
	\rm -f na.nad_demo ;\
	ln -s ../../bin/rfi2nad . ;\
	ln -s ../../bin/readnad . ;\
	./rfi2nad rfi_models_500_20_2 na.nad_demo ;\
	\rm -f na.nad_500_20_2 ;\
	ln -s na.nad_demo na.nad_500_20_2 ;\
	\rm -f junk

nad2nad: nad2nad.o
	$(F77) $(FFLAGS) nad2nad.o -o $(DEMO)/bin/nad2nad 

nad2asc: nad2asc.o
	$(F77) $(FFLAGS) nad2asc.o -o $(DEMO)/bin/nad2asc 

asc2nad: asc2nad.o
	$(F77) $(FFLAGS) asc2nad.o -o $(DEMO)/bin/asc2nad 

nad_sun2dec: nad_sun2dec.o
	$(F77) $(FFLAGS) nad_sun2dec.o -o $(DEMO)/bin/nad_sun2dec 
	cd $(DEMO)/data/demo_results;\
	\rm -f na.nad_500_20_2 ;\
	../../bin/nad_sun2dec na.nad_demo_sun na.nad_demo_dec ;\
	ln -s na.nad_demo_dec na.nad_500_20_2;\


.f.o    :
	$(F77) $(FFLAGS) -O -c $*.f

