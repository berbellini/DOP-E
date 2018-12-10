#
#	Makefile for program splot
#

# root directory for Demo program

DEMO    = ../../..

# Edit this file to change location or
# type of compilers (e.g. to gcc, g77)

include $(DEMO)/src/makena.macros 

LDIR    = -L./Lib

all:	$(X)

yes:	libs splot-x splot-xb splot-p

no:	libs splot-p

splot-x : 
	$(F77) $(FFLAGS) $(LDIR) $(F77_COMPAT) splot.F -lxpakw $(XFLAGS) \
         -o $(DEMO)/bin/splot-x

splot-xb : 
	$(F77) $(FFLAGS) $(LDIR) $(F77_COMPAT) splot.F -lxpak $(XFLAGS) \
         -o $(DEMO)/bin/splot-xb

splot-p : 
	$(F77) $(FFLAGS) $(LDIR) $(F77_COMPAT) splot.F -lpspak \
         -o $(DEMO)/bin/splot-p

libs : 
	cd subs;make -f Makefile
