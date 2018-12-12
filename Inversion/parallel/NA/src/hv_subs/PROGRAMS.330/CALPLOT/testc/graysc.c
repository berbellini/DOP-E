/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GRAYSC                                                c
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

#include "calplot.h"
void box(float x0,float y0,float x1,float y1);
void graysc(float x0,float y0,float x1,float y1,int ipen0,int ipen1,int ipinc);

int main()
{
	int i;
	float xx,yy;
	pinitf("GRAYSC.PLT");
	gunit("in");
/*
c-----
c	draw some colored lines
c-----
*/
	for(i=0;i<20;i++){
		newpen(i);
		xx = 8.0;
		yy = 7.0 - 0.3*i;
		gwidth(0.10);
		plot(xx+0.5,yy,3);
		plot(xx+1.00,yy,2);
		newpen(1);
		gwidth(0.01);
		number(xx+1.1,yy,0.10,(float)i,0.0,-1);
	}
	gwidth(0.01);
/*
c-----
c	put out a nice gray scale according to ipen
c-----
*/
	graysc(1.0,1.0,1.5,7.0,1000,1020,1);
	graysc(2.5,1.0,3.0,7.0,1020,1040,1);
	graysc(4.0,1.0,4.5,7.0,1040,1060,1);
	graysc(5.5,1.0,6.0,7.0,1060,1080,1);
	graysc(7.0,1.0,7.5,7.0,1080,1100,1);
	pend();
}

void box(float x0,float y0,float x1,float y1)
{
	plot(x0,y0,3);
	plot(x0,y1,2);
	plot(x1,y1,2);
	plot(x1,y0,2);
	plot(x0,y0,2);
}

void graysc(float x0,float y0,float x1,float y1,int ipen0,int ipen1,int ipinc)
{
	int ipnmn, ipnmx, ipnd, ipen, ipatx, ipaty;
	float xlen, ylen, dy, xl, xh, yl, yh;
	box(x0,y0,x1,y1);
	ipnmn=ipen0;
	ipnmx=ipen1;
	ipnd=ipnmx-ipnmn;
	ylen = y1 - y0;
	dy =  ylen/(ipnd+1);
	for(ipen=ipnmn;ipen<=ipnmx;ipen+=ipinc){
		xl=x0;
		xh=x1;
		yh = y1 - (ipen-ipnmn)*dy;
		yl = yh - dy*ipinc;
		if(yl < y0)yl=y0;
		ipatx=0;
		ipaty=0;
		xlen=0.01;
		ylen=0.01;
		newpen(ipen);
		shader(xl,yl,xh,yh,ipatx,ipaty,xlen,ylen);
		newpen(1);
		number(xh+0.1,yh-dy/2.,0.10,(float)ipen,0.0,-1);
	}
}
