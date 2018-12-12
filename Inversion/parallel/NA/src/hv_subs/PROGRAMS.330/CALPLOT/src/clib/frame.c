/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: FRAME                                                 c
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

void frame(void )
/*
c-----
c       start new plot on new page
c       move current position to lower left corner
c       define new origin
c       do not redefine factor scaling, thougn
c-----
*/
{
	float xfac, yfac;
	int mthree = -3;
	xfac = -Gcplt.xold/Gcplt.xstp;
	yfac = -Gcplt.yold/Gcplt.ystp;
        plot(xfac,yfac,mthree);
        Gcplt.xcur=0.0;
        Gcplt.ycur=0.0;
        (*do_erase)(0);
}

void gframe(int mode )
/*
c-----
c       start new plot on new page
c       move current position to lower left corner
c       define new origin
c       do not redefine factor scaling, thougn
c-----
*/
{
	float xfac, yfac;
	INT pode;
	int mthree = -3;
	xfac = -Gcplt.xold/Gcplt.xstp;
	yfac = -Gcplt.yold/Gcplt.ystp;
        plot(xfac,yfac,mthree);
        Gcplt.xcur=0.0;
        Gcplt.ycur=0.0;
	pode = mode;
        (*do_erase)(pode);
}
