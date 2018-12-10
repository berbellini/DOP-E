#
#	Make discription file for NA
#

#FFLAGS = -O
#F77 = /opt/local/SUNWspro/bin/f77

DEMO = ../../..

# Edit this file to change location or
# type of compilers (e.g. to gcc, g77)

include $(DEMO)/src/makena.macros

.f.o    :
	$(F77) $(FFLAGS) -c $*.f

readnad: readnad.o
	$(F77) $(FFLAGS) readnad.o -o $(DEMO)/bin/readnad 
