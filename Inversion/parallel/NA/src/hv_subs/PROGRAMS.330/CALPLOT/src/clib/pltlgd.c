/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLTLGD                                                c
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

#include <math.h>
#include "calplot.h"

void pltlgd(float *x,float *y,int n,float x1,float y1,
		float deltax,float deltay,int lintyp,int inteq,
		float ht,int nocx,int nocy,int ipat,float xlen)
{
/*
c-----
c	plot y(x) in linear, log or semilog plots. Check
c	for bounds of plot
c	ipat	integer providing bit pattern for dashed lines
c	xlen	length of each on bit in pattern ipat
c-----
*/
	float hht, xmin, xmax, xx, yy;
	float ymin, ymax;
	int i, ipen;
	char c[4];

	if(ht <= 0.0){
		hht = 0.07;
	} else {
		hht = ht;
	}
	if(nocx > 0){
		xmin = pow(10.0,x1);
		xmax = pow(10.0,(x1+(float)nocx));
	}
	if(nocy > 0){
		ymin = pow(10.0,y1);
		ymax = pow(10.0,(y1+(float)nocy));
	}
/*
c-----
c	plot, checking that bounds are not exceeded in log plots
c-----
*/
	ipen = 3;
	c[0] = (char)inteq;
	for(i=0;i<n;i++){
		if(nocx > 0){
			if(x[i] >= xmin  &&  x[i] <= xmax){
				xx = (log10(x[i])-x1)*deltax;
			} else if(x[i] < xmin){
				xx = 0.0;
			} else {
				xx =  nocx * deltax;
			}
		} else {
			xx = (x[i]-x1)/deltax;
		}
		if(nocy > 0){
			if(y[i] >= ymin  &&  y[i] <= ymax){
				yy = (log10(y[i])-y1)*deltay;
			} else if(y[i] < ymin){
				yy = 0.0;
			} else {
				yy =   deltay *nocy;
			}
		} else {
			yy = (y[i]-y1)/deltay;
		}
		if(ipat == 0 || ipen == 3){
			plot(xx,yy,ipen);
		} else {
			plotd(xx,yy,ipat,xlen);
		}
		if(lintyp != 0)symbol(xx,yy,hht,c,0.0,-1);
		if(lintyp >= 0)ipen=2;
	}
	return;
}
