/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GWIDTH                                                c
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

#include "callib.h"
#include "dosubs.h"
#include "calplot.h"

extern struct Gcplot Gcplt;

void gwidth(float width)
{
	INT ix;
	float xx;
	Gcplt.owidth = width;
	xx = 1000.0*Gcplt.xstp*Gcplt.owidth;
	if(xx > 1000000000.0)
		ix = 1000000000;
	else
		ix = xx;
	if(ix < 0)ix = 0;
	(*do_gwid)(ix);
}
