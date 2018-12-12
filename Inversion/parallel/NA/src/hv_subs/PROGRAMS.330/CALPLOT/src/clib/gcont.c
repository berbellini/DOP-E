/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GCONT                                                 c
c                                                                     c
c      COPYRIGHT (C)  2004 R. B. Herrmann                             c
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

void gcont(float x,float y)
{
/*
c-----
c       basic routine to move point on plotter
c
c       x       abscissa in inches
c       y       ordinate in inches
c-----
c-----
c      xold and yold are absolute plotter coordinates
c      xcur and ycur are coordinates of point relative
c            to current plotter position in current
c            software units
c-----
*/
        INT ix,iy;
	float xn, yn, xx, yy;
/*
c-----
c	x and y are page coordinates in units given by
c	last call gunit()
c-----
*/
	Gcplt.xcur=x;
	Gcplt.ycur=y;
	xn=Gcplt.xcur*Gcplt.xstp;
	yn=Gcplt.ycur*Gcplt.ystp;
/*
c-----
c	current position is in inches, ala calcomp
c	but to use unix plot filters we need integer*2
c	the factor 1000 is really 11000 counts / 11.0 inches
c-----
*/
       	xx = 1000. * (xn + Gcplt.xold);
       	yy = 1000. * (yn + Gcplt.yold);
	if(xx > 1000000000.0)
		ix = 1000000000;
	else if(xx < -1000000000.0)
		ix = -1000000000;
	else
		ix = xx ;
	if(yy > 1000000000.0)
		iy = 1000000000;
	else if(yy < -1000000000.0)
		iy = -1000000000;
	else
		iy = yy ;
       	(*do_cont)(ix,iy);
}
