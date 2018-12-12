/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GPHTXT                                                c
c                                                                     c
c      COPYRIGHT (C)  1997 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      3507 Laclede Avenue                                            c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
*/

#include "sysinit.h"
#include "callib.h"
#include "calplot.h"

extern struct Gcplot Gcplt;
void gottxt(int cnt, char *s);
void gintxt(int cnt, char *s);



/*
c-----
c	routines to implement string output to terminals
c	on some terminals this will use the terminal character generator, e.g.,
c	tektronix, on others it will use a bit map, e.g., PC, and others
c	will use a vector draw
c-----
*/

void gwrtxt(float xx,float yy,char *text,int flgbln)
{
	int lstr;
	plot(xx,yy,3);
	lstr = strlen(text);
	newpen(2000 + flgbln);
	gottxt(lstr,text);
}

void grdtxt(char *text,int lstr)
{
	gintxt(lstr,text);
}
