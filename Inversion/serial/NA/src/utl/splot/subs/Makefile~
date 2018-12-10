#
#	Makefile for program splot
#

# root directory for Demo program

DEMO    = ../../../..

# Edit this file to change location or
# type of compilers (e.g. to gcc, g77)

include $(DEMO)/src/makena.macros 


SUBS    = ps_comm.o ps_colr.o cfopen.o
SUBSX   = xpakw.o xpakb.o xtra.o

all:	$(X) 

yes:	libs libsx
no:	libs

.f.o	:
	$(F77) $*.f -O -c

.c.o	:
	$(CC) $*.c -O -c

libs : $(SUBS)
	$(AR) rv ../Lib/libpspak.a ps_comm.o ps_colr.o

libsx : $(SUBSX)
	$(AR) rv ../Lib/libxpakw.a xpakw.o xtra.o cfopen.o
	$(AR) rv ../Lib/libxpak.a xpakb.o xtra.o cfopen.o
