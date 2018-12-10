/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: NEWPEN                                                c
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

void newpen(int j)
{
/*
c-----
c       j is number of pen to use
c	limited in the device specific file
c	specific device limitations are defined
c	later down the line
c-----
*/
        INT i;
	i=j;
	if(i < 0)i=0;
	(*do_pen)(i);
}
