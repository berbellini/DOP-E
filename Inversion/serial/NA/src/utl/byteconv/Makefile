# Make for sacsun2dec, 1998-08-08
# T. Shibutani, RCEP, DPRI, Kyoto Univ.

DEMO = ../../..

include $(DEMO)/src/makena.macros   # Edit this file to change location or
       	                            # type of compilers (e.g. to gcc, g77)

#CCOPT=-O4

.c.o:
	$(CC) $(CCOPT) -I. -c $*.c

sacsun2dec : sacsun2dec.o sac.h
	$(CC) $(CCOPT) -o sacsun2dec sacsun2dec.o

