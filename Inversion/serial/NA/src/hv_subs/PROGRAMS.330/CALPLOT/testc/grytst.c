/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GRYTST                                                c
c                                                                     c
c      COPYRIGHT (C)  1986 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c-----
c	this program tests an improvement to the
c	CALPLOT device drivers to permit gray scale plots under
c	certain circumstances. This will be useful for plotting
c	seismic traces together with a trace attribute modulating
c	either the shading or the color of the background
c
c	To do this we invoke solid rectangle shading with ipatx = ipaty = 0
c	but with the color in newpen(ipen) > 1000. The present gray
c	shading will fill solid with ipen=1000, and with an
c	areal density decreasing with increasing ipen
c-----
*/

#include <math.h>
#include	"calplot.h"

main()
{
	float xwid, ywid, dx, dy, ampmin, ampmax, x, y, amp, xfac, yfac;
	int nx, ny, ishdrv, jplot, ix, iy;
	float tri();
	pinitf("GRYTST.PLT");
	gunit("in");
/*
c-----
c	program to test gray shading
c-----
*/
	xwid = 6.0;
	ywid = 6.0;
	nx = 40;
	ny = 40;
	dx = xwid/nx;
	dy = ywid/ny;
	ampmin = -1.0;
	ampmax = 1.0;
	ishdrv = 1;
/*
c-----
c	jplot = 1 multiply sin functionc
c	jplot = 2 multiply triangular functions
c-----
*/
	for(jplot=1;jplot<=1;jplot++){
		plot(1.0,1.0,-3);
		for(ix=1;ix<=nx;ix++){
			x = (ix-1)*dx;
			if(jplot == 1){
				xfac = sin(3.1415927*x/xwid);
			} else {
				xfac = tri(x,xwid);
			}
			for(iy=1;iy<=ny;iy++){
				y = (iy-1)*dy;
				if(jplot == 1){
					yfac=sin(6.2831853*y/ywid);
				} else {
					yfac = tri(y+y,ywid);
				}
				amp = xfac*yfac;
				putshd(x,x+dx,y,y+dy,amp,ishdrv,ampmin,ampmax);
			}
		}
		box(0.0,0.0,xwid,ywid);
/*
c-----
c		put out a nice gray scale according to ipen
c-----
*/
		graysc(ampmin,ampmax,ishdrv);
		plot(-1.0,-1.0,-3);
		if(jplot == 1)frame();
	}
	pend();
}

float tri(x,xl)
float x, xl;
{
/*
c-----
c	function to produce triangular function with a period
c	of 2*xl
c-----
c-----
c	first reduce periodicity to variable xval between [0,1]
c-----
*/
	float xval; 
	int ival;
	xval = x/(xl+xl);
	ival = xval;
	xval = xval - (float)ival;
	if(xval <= 0.25){
		return(4.0*xval);
	} else if(xval > 0.25 && xval <= 0.75){
		return( 4.0*(0.50 - xval));
	} else if(xval > 0.75 && xval <= 1.00){
		return( 4.0*(xval - 1.0));
	}
}

	
box(x0,y0,x1,y1)
float x0, y0, x1, y1;
{
	plot(x0,y0,3);
	plot(x0,y1,2);
	plot(x1,y1,2);
	plot(x1,y0,2);
	plot(x0,y0,2);
	plot(x0,y0,2);
}
putshd(xl,xh,yl,yh,amp,ishdrv,ampmin,ampmax)
float xl, xh, yl, yh, amp, ampmin, ampmax;
int ishdrv;
{
	float xcol, xlen, ylen;
	int ipen, ipatx, ipaty;
	xcol = (amp - ampmin)/(ampmax-ampmin);
	xcol = 100.0*xcol  + 0.5;
	ipen = 1000 + xcol;
	if(ishdrv  < 0)ipen = 1100 - xcol;
	if(ipen  <  1000)ipen = 1000;
	if(ipen >  1100)ipen=1100;
	ipatx = 0;
	ipaty = 0;
	xlen = 0.02;
	ylen = 0.02;
/*
c-----
c	ipen = 1000 is red, 1100 = blue or 1000 = dark, 1100 = light halftone
c-----
*/
	newpen(ipen);
	shader(xl,yl,xh,yh,ipatx,ipaty,xlen,ylen);
}

graysc(ampmin,ampmax,ishdrv)
float ampmin, ampmax;
int ishdrv;
{
	int nbox, i;
	float damp, xl, xh, dy, amp, yl, yh;
	nbox=21;
	damp = (ampmax-ampmin)/(nbox-1);
	xl=7.0;
	xh=7.5;
	dy = 6.0/nbox;
	for(i=1;i<=nbox;i++){
		amp = ampmin + (i-1)*damp;
		yl=(i-1)*dy;
		yh=yl+dy;
		putshd(xl,xh,yl,yh,amp,ishdrv,ampmin,ampmax);
		newpen(1);
		number(xh+0.1,yh-dy/2.,0.07,amp,0.0,2);
	}
	box(7.0,0.0,7.5,6.0);
}
