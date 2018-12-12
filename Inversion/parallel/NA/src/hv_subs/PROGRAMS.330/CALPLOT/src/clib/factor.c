/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: FACTOR                                                c
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
/*
c-----
c       multiply all pen moves by fact
c       useful in reducing size of plots, just
c       do a call factor immediately after
c       the call plots or call pinit
c----
*/

#include "callib.h"
#include "calplot.h"

extern struct Gcplot Gcplt;

void factor (float fact) 
{
	float factnw;
	int oneooone = 1001;
	if(fact <= 0.0)
		factnw=1.0;
	else
		factnw = fact;
	if(Gcplt.iunit == 1)
		factnw = factnw * 0.3937;
	plot(factnw,factnw,oneooone);
}
