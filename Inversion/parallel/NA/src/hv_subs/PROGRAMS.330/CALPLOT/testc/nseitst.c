/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: NSEITST                                               c
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
#include "calplot.h"
#define NPTS 81

main()
{
	float x[NPTS],y[NPTS],amp[NPTS];
	float xl, xh, ampav, xcol;
	int i, ishdrv, ipen;
	float	ampmin =  1.0e+38,
		ampmax = -1.0e+38,
		tramin =  1.0e+38,
		tramax = -1.0e+38;
	pinitf("NSEITST.PLT");
	gunit("in");
/*
c-----
c	program to test seismic trace shading
c-----
c	fill up trace array and the attribute array
c-----
*/
	for(i=0;i< NPTS; i++){
		x[i] = i*0.05;
		amp[i] = fabs(0.25 * cos(3.1415927*x[i]/2.));
		if(ampmin > amp[i])ampmin=amp[i];
		if(ampmax < amp[i])ampmax=amp[i];
		y[i] = amp[i]*sin(3.1415927*x[i]*5./2.);
		if(tramin > fabs(y[i]) )tramin=fabs(y[i]);
		if(tramax < fabs(y[i]) )tramax=fabs(y[i]);
	}
	puttrc(1.0,1.0,x,y,NPTS);
	symbol(5.2,1.0,0.1,"TRACE",0.0,5);
	ishdrv = 1;
	putshd(1.0,1.5,x,amp,NPTS,ishdrv,ampmin,ampmax);
	symbol(5.2,1.5,0.1,"SHADE",0.0,5);
	plot(1.0,2.0,-3);
	for(i=1;i<NPTS;i++){
		ampav = (amp[i-1]+amp[i])/2;
		xl = x[i-1];
		xh = x[i];
		xcol = (ampav - ampmin)/(ampmax -ampmin);
		if(i == 1){
			ipen = 3;
		} else {
			ipen = 2;
		}
		plot(xl,xcol,ipen);
		plot(xh,xcol,2);
	}
	plot(-1.0,-2.0,-3);
	symbol(5.2,2.0,0.10,"ATTRIBUTE",0.0,9);
	ishdrv = 1;
	putshd(1.0,3.5,x,amp,NPTS,ishdrv,ampmin,ampmax);
	puttrc(1.0,3.5,x,y,NPTS);
	symbol(5.2,3.5,0.10,"ZERO ATTRIB=DARK",0.0,16);
	symbol(5.2,3.3,0.10,"ZERO ATTRIB=RED ",0.0,16);
	ishdrv = -1;
	putshd(1.0,4.5,x,amp,NPTS,ishdrv,ampmin,ampmax);
	puttrc(1.0,4.5,x,y,NPTS);
	symbol(5.2,4.5,0.10,"MAX  ATTRIB=DARK",0.0,16);
	symbol(5.2,4.3,0.10,"MAX  ATTRIB=RED ",0.0,16);
/*
c-----
c	put out a nice gray scale according to ipen
c-----
*/
	graysc();
/*
c-----
c       put up a trace, but have the shading keyed according to the
c       peak amplitude
c       This is tricky, since we must first search for local maxima
c-----
*/
        puttra(1.0,5.3,x,y,NPTS,tramax,tramin,0);
        symbol(5.2,5.3,0.10,"POSITIVE  TRACE SHADE",0.0,21);
        puttra(1.0,5.8,x,y,NPTS,tramax,tramin,1);
        symbol(5.2,5.8,0.10,"NEGATIVE  TRACE SHADE",0.0,21);
        puttra(1.0,6.3,x,y,NPTS,tramax,tramin,2);
        symbol(5.2,6.3,0.10,"TWO-SIDED TRACE SHADE",0.0,21);
	pend();
}

puttrc(x0,y0,x,y,n)
float x0, y0, x[], y[];
int n;
{
	int i;
	float xx, yy;
	plot(x0,y0,-3);
	shdsei(0.0,0.0,0,1,0);
	for(i=0;i<n;i++){
		xx = x[i];
		yy = y[i];
		if(i == 0){
			plot(xx,yy,3);
		} else {
			plot(xx,yy,2);
		}
	}
	shdsei(0.0,0.0,0,0,0);
	plot(-x0,-y0,-3);
}
 
#define SIGN(X) ( ( (X)> 0) ? 1.0 : -1.0 )

puttra(x0,y0,x,y,n,ampmx,ampmn,iplmn)
float x0, y0, x[], y[], ampmx, ampmn;
int n, iplmn;
{
	int kolor[NPTS];
	float ampmax, x2, color,fac,xx,yy;
	int i, j, kolr, iold;
/*
c-----
c	The idea here is to follow the trace output, saving the
c	current maximum. At a zero crossing, we output the
c	color corresponding to this maximum, and reset
c	we must have the color attribute defined beforehand
c-----
*/
	ampmax = fabs(y[0]);
	iold = 0;
	for(i=1;i<n;i++){
		if(fabs(y[i]) > ampmax)ampmax = fabs(y[i]);
		if(SIGN(y[i]) != SIGN(y[i-1]) ) {
			color =  1100.0 - 100.0*(ampmax-ampmn)/(ampmx - ampmn);	
			kolr = color;	
			if(kolr < 1000)kolr = 1000;	
			if(kolr > 1100)kolr = 1100;	
			for(j=iold; j< i-1;j++)
				kolor[j] = kolr;	
			iold = i;	
			ampmax = 0.0;	
		}
	}
	kolor[n-1] = kolor[n-2];
	plot(x0,y0,-3);	
	shdsei(0.0,0.0,0,1,iplmn);	
	for(i=0;i<n;i++){
		xx = x[i];
		yy = y[i];
		if(i == 0) {
			newpen(kolor[0]);
			plot(xx,yy,3);
		} else {
			if(SIGN(y[i]) != SIGN(y[i-1]) ){
				fac = (y[i] - y[i-1])/(x[i] - x[i-1]);
				x2 = x[i-1] - y[i-1]/fac;
				plot(x2,0.0,2);
				newpen(kolor[i-1]);
/*
c-----
c	force plotting of segment in this color
c-----
*/
				shdsei(0.0,0.0,0,0,0);
/*
c-----
c	initialize plotting of segment for next color
c-----
*/
				shdsei(0.0,0.0,0,1,iplmn);	
				plot(x2,0.0,3);
				newpen(kolor[i]);
				plot(xx,yy,2);
			} else {
				plot(xx,yy,2);
			}
		}
	}
	shdsei(0.0,0.0,0,0,0);
/*
c-----
c	now plot the trace in black so it can be seen
c-----
*/
	newpen(1);
	for(i=0;i<n;i++){
		xx = x[i];
		yy = y[i];
		if(i == 0){
			plot(xx,yy,3);
		} else {
			plot(xx,yy,2);
		}
	}
	plot(-x0,-y0,-3);
}

putshd(x0,y0,x,amp,n,ishdrv,ampmin,ampmax)
float x0, y0, x[], amp[], ampmin, ampmax;
int n, ishdrv;
{
	int i, ipen, ipatx, ipaty;
	float xl, xh, yl, yh, ampav, xcol, xlen, ylen;
	plot(x0,y0,-3);
	for(i=1;i<n;i++){
		xl = x[i-1];
		xh = x[i];
		ampav = (amp[i-1] + amp[i])/2;
		xcol = (ampav - ampmin)/(ampmax-ampmin);
		xcol = 100.0*xcol  + 0.5;
		ipen = 1000 + xcol;
		if(ishdrv < 0)ipen = 1100 - xcol;
		if(ipen  < 1000)ipen = 1000;
		if(ipen > 1100)ipen=1100;
		yl = -0.25;
		yh =  0.25;
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
/*
c-----
c	undo special shading
c-----
*/
	newpen(1);
	plot(-x0,-y0,-3);
}

graysc()
{
	int ipnmn, ipnmx, ipnd, ipatx, ipaty, ipen, ipinc;
	float dy, xlen, ylen, yl, yh, xl, xh;
	plot(7.0,5.0,-3);
	plot(0.0,-4.0,2);
	plot(0.5,-4.0,2);
	plot(0.5,0.0,2);
	plot(0.0,0.0,2);
	ipnmn=1000;
	ipnmx=1100;
	ipnd=ipnmx-ipnmn;
	ipinc=5;
	dy = -4.0/ipnd;
	for(ipen=ipnmn;ipen < ipnmx;ipen+=ipinc){
		xl=0.0;
		xh=0.5;
		yl=(ipen-ipnmn)*dy;
		yh=yl+dy*ipinc;
		ipatx=0;
		ipaty=0;
		xlen=0.01;
		ylen=0.01;
		newpen(ipen);
		shader(xl,yl,xh,yh,ipatx,ipaty,xlen,ylen);
		newpen(1);
		number(xh+0.1,yh-dy/2.,0.07,(float)ipen,0.0,-1);
	}
	plot(-7.0,-5.0,-3);
}
