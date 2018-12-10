/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PEND                                                  c
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
c-----
c      TERMINATE PLOT
c-----
*/
#include "calplot.h"

void clospl(int mode);


void pend(void )
{
	/* reset to default width and pen */
	gwidth(0.0);
	newpen(1);
	clospl(0);
}

void gend(int mode )
{
	/* reset to default width and pen */
	gwidth(0.0);
	newpen(1);
	if(mode >1 || mode < 0)
		clospl(0);
	else
		clospl(mode);
}
