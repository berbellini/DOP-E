/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SHDSEI                                                c
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


extern struct Gcplot Gcplt;

void shdsei(float x1,float y1,int ixy,int istnd,int iplmn)
{
/*
c-----
c	permit shading of seismic traces
c-----
c	x1	- x coordinate of first point local zero
c	y1	- y coordinate of first point local zero
c	ixy	- time axis is parallel to x-axis (0)
c		                           y-axis (1)
c	istnd   - 0 turn off shading  
c		  1 turn on  shading
c	iplmn	- 0 shade positive amplitudes
c		- 1 shade negative amplitudes
c-----
c	For example, if time is parallel to x-axis
c	and iplmn = 0, then values of y > y1 will be shaded from y to y1
c-----
c-----
c      xold and yold are absolute plotter coordinates
c      xcur and ycur are coordinates of point relative
c            to current plotter position in current
c            software units
c-----
*/
	INT jx1, jy1, jxy, jstnd, jplmn;
	INT ix, iy;
	float xn, yn, xx, yy;
	Gcplt.xcur = x1;
	Gcplt.ycur = y1;
	xn=Gcplt.xcur*Gcplt.xstp;
	yn=Gcplt.ycur*Gcplt.ystp;
/*
c-----
c current position is in inches, ala calcomp
c but to use unix plot filters we need integer*2
c the factor 1000 is really 11000 counts / 11.0 inches
c-----
*/
        xx = 1000. * (xn + Gcplt.xold);
        yy = 1000. * (yn + Gcplt.yold);
	if(xx > 100000000.0)
		ix = 100000000;
	else if(xx < -100000000.0)
		ix = -100000000;
	else
		ix = xx;
	if(yy > 100000000.0)
		iy = 100000000;
	else if(yy < -100000000.0)
		iy = -100000000;
	else
		iy = yy;
	jx1 = ix;
	jy1 = iy;
	jxy = ixy;
	jstnd = istnd;
	jplmn = iplmn;
	(*do_fills)(jx1,jy1,jxy,jstnd,jplmn);
}
