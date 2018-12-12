/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLOTD                                                 c
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
#include <math.h>
#include "calplot.h"

static void segplt(float *s,float sx,float sy,float xlen,float *xs,
	float *ys,int ipat,int *iret);

extern struct Gcplot Gcplt;

/*
c-----
c	is	- index to bit in pattern ipat
c	xsv	- offset within current bit in xlen units
c	ipatld	- previous pattern. If pattern changes reinitialize
c-----
*/


void plotd(float xx,float yy,int ipat,float xlen)
{
/*
c-----
c	draw a line between the current point
c	and the coordinates (xx,yy) using
c	the pattern ipat. ipat is considered
c	to be a binary pattern with 5 bits
c	set to either 0 or 1. The length of
c	each bit is xlen. The 0 bit indicates
c	a dark vector and the 1 bit a plotted vector
c-----
*/
	float x0, y0, fct, s, sx, sy, xs, ys;
	int iret;

	if(ipat != Gcplt.ipatld){
		Gcplt.is = 0;
		Gcplt.xsv = 0.0;
		if(ipat > 31 || ipat < 0)ipat = 0;
		Gcplt.ipatld = ipat;
	}
/*
c-----
c	test for solid line
c-----
*/
	if(ipat == 0){
		plot(xx,yy,2);
		return;
	}
/*
c-----
c	get current pen position
c-----
*/
	where(&x0,&y0,&fct);
/*
c-----
c	determine length of segment to be plotted, is zero return
c-----
*/
	s = sqrt( (xx-x0)*(xx-x0) + (yy-y0)*(yy-y0) );
	if(s <= 0.001)
		return;
/*
c-----
c	determine direction cosines for parametric representation of line
c	also remember that the ultimate precision is 0.001 units in plot space
c-----
*/
	if(fabs(xx-x0) < 0.001)
		sx = 0.0;
	else
		sx = (xx-x0)/s;
	if(fabs(yy-y0) < 0.001)
		sy = 0.0;
	else
		sy = (yy-y0)/s;
	xs = x0;
	ys = y0;
/*
c-----
c	plot segments
c-----
*/
	iret = 1;
	do {
		segplt(&s,sx,sy,xlen,&xs,&ys,ipat,&iret);
	} while(iret == 1);
/*
c-----
c	exit subroutine
c-----
*/
	return;
}

static void segplt(float *s,float sx,float sy,float xlen,float *xs,
	float *ys,int ipat,int *iret)
{
/*
c-----
c	plot bit pattern segments, either partial or complete
c-----
c	s	- length of segment yet to be plotted
c	sx	- direction cosine in x-direction
c	sy	- direction cosine in y-direction
c	xlen	- length of bit segment
c	xsv	- fraction of current segment drawn
c	is	- index to bit in pattern ipat
c	xs	- current x position
c	ys	- current y position
c	ipat	- bit pattern for plot
c	iret	- -1 partial segment drawn
c	iret	- +1 complete segment drawn
c-----
*/
	float dsmx, dsx, dsy;
	int ipen, i;
	dsmx = (1.0- (Gcplt.xsv))*xlen;
	dsx = dsmx*(sx);
	dsy = dsmx*(sy);
	if(*s < dsmx){
		(*xs) = (*xs) + (*s)*(sx);
		(*ys) = (*ys) + (*s)*(sy);
		(Gcplt.xsv) = (Gcplt.xsv) + (*s)/xlen ;
		if(Gcplt.xsv >= 1.0){
			Gcplt.xsv=0.0;
			Gcplt.is = Gcplt.is + 1;
			if(Gcplt.is > 5)Gcplt.is = 0;
		}
		*iret = -1;
	} else {
		*xs = *xs + dsx;
		*ys = *ys + dsy;
		Gcplt.xsv = 0.0;
		*iret = 1;
		if(dsmx > 0.0){
			*s = *s -dsmx;
		} else {
			*s = *s - xlen;
		}
		if(*s < 0.0)*s = 0.0;
	}
	/*
	if( ((ipat) >> (Gcplt.is) ) && 01)
	*/
	i = ipat /(1 << (int) Gcplt.is);
	i%=2;
	if(i == 1)
		ipen = 2;
	else
		ipen = 3;
	plot(*xs,*ys,ipen);
	if(*iret == 1){
		Gcplt.is = Gcplt.is + 1;
		if(Gcplt.is > 5)
			Gcplt.is=0;
	}
}
