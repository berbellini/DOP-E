/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLTLOG                                                c
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
#include <math.h>
#include "calplot.h"

extern struct Gcplot Gcplt;


void pltlog(float x[],float y[],int n,float x1,float y1,
	float deltax,float deltay,int lintyp,int inteq,
	float ht,int nocx,int nocy)
{
/*
c-----
c	plot y(x) in linear, log or semilog plots. Check
c	for bounds of plot
c-----
*/
	float hht, xmin, xmax, ymin, ymax, xx, yy;
	int ipen, i;
	float zero = 0.0;
	int mone = -1;
	char is[2];

	is[0] = (char)inteq;
	if(ht <= 0.0)
		hht = 0.07;
	else
		hht = ht;
	if(nocx > 0){
		xmin = pow(10.0, x1);
		xmax = pow(10.0, x1 + (float)nocx);
	}
	if(nocy > 0){
		ymin = pow(10.0, y1);
		ymax = pow(10.0, y1 + (float)nocy);
	}
/*
c-----
c	plot, checking that bounds are not exceeded in log plots
c-----
*/
	ipen = 3;
	for(i=0;i<n;i++){
		if(nocx > 0){
			if(x[i] >= xmin && x[i] <= xmax){
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
			if(y[i] >= ymin && y[i] <= ymax){
				yy = (log10(y[i])-y1)*deltay;
			} else if(y[i] < ymin){
				yy = 0.0;
			}else{
				yy =   deltay *nocy;
			}
		} else {
			yy = (y[i]-y1)/deltay;
		}
		plot(xx,yy,ipen);
		if(lintyp != 0)symbol(xx,yy,hht,is,zero,mone);
		if(lintyp >= 0)ipen=2;
	}
}
