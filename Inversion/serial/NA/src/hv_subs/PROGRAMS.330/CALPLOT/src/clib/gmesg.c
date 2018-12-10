/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GMESG                                                 c
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

void gmesg(char *mesg);
void gomesg(int cnt, char *mesg);

/* routine to but message into space in menu bar */

void gmesg(char *mesg)
{
	int lstr;
	lstr = strlen(mesg);
	gomesg(lstr,mesg);
}

