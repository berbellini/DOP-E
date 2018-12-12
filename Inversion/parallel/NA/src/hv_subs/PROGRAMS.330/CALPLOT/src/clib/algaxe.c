/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: ALGAXE                                                c
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

void algaxe(float xaxlen,float yaxlen,int nocx,int nocy,char *ttlx,char *ttly,int mtx,int mty,float x1,float y1,float deltax,float deltay)
{
/*
c-----
c	ALGAXE
c		PLOTS AXES
c		xaxlen	-	length of x axis in inches
c				(if <= 0.0 do not plot )
c		yaxlen	- 	length of y axis in inches
c				(if <= 0.0 do not plot )
c		nocx	- 	number of logarithmic cycles
c				on x axis if log plot
c				= 0 if linear plot
c		nocy	-	number of logarithmic cycles
c				= 0 if linear plot
c		ttlx	-	x axis title a character string
c		ttly	-	y axis title a character string
c		mtx	-	number of characters in x axis
c				title
c		mty	-	number of characters in y ayis
c				title
c		x1	-	first value in x-axis
c		y1	-	first value in y-axis
c				if log plot these are the powers of 10
c		deltax	-	if linear number of x-units per inch
c				if log length of cycle
c		deltay	-	yaxis
c				on y axis if log plot
c-----
*/
	int i, j;
	float slt, sst, sp, ss, ssp;
	float ttlp, sttl, xnum, xpo, ypo;
	float yl, yu, x, xtl, y, ytl;
	float factx, facty;

	slt = 0.02*yaxlen ;
	sst = 0.01 * yaxlen ;
	sp = -0.06*yaxlen ;
	ss = 0.035*yaxlen ;
	ssp = sp + ss - 0.06 ;
	ttlp = -0.11*yaxlen - 0.1 ;
	sttl = 0.035*yaxlen ;
	xnum = 1 ;
	yl = y1 ;
	yu = y1 + (float)nocy ;
	if(fabs(yl) >= 10.  || fabs(yu) >= 10. )xnum = xnum + 1. ;
	if(fabs(yl) >= 100. || fabs(yu) >= 100.)xnum = xnum + 1. ;
	if(y1 < 0) xnum = xnum + 1.0 ;
	xpo = x1 ;
	ypo = y1 ;
/*
c-----
c	draw x-axis
c-----
*/
	if(xaxlen > 0.0){
		if(nocx == 0) {
			axis(0.0,0.0,ttlx,-mtx,xaxlen,0.0,x1,deltax);
		} else {
			plot(0.0,-slt,3);
			plot(0.0,0.0,2);
			factx = xaxlen/(float)nocx;
			symbol(-.6*ss,sp,ss,"10",0.0,2);
			number(999.0,ssp,0.5*ss,x1,0.0,-1);
			plot(0.0,0.0,3);
			for(j=1;j<=nocx;j++){
				for(i=1;i<=10;i++){
					x = i;
					x = log10(x) *factx + (j-1)*factx;
					if(i != 1){
						plot(x,0.0,2);
						plot(x,-sst,2);
					}
					plot(x,0.0,3);
				}
				plot(x,-slt,2);
				symbol(x-.6*ss,sp,ss,"10",0.0,2);
				xpo = xpo + 1.0;
				number(999.0,ssp,0.5*ss,xpo,0.0,-1);
				plot(x,0.0,3);
			}
			xtl = mtx;
			xtl = (xaxlen-xtl*sttl)/2.0;
			symbol(xtl,ttlp,sttl,ttlx,0.0,mtx);
		}
		plot(0.0,0.0,3);
	}
/*
c-----
c	draw y-axis
c-----
*/
	if(yaxlen > 0.0){
		if(nocy == 0) {
			axis(0.0,0.0,ttly,mty,yaxlen,90.,y1,deltay);
		} else {
			plot(-slt,0.0,3);
			plot(0.0,0.0,2);
			sp = sp - (xnum - 1.5) * 0.5 * ss;
			ttlp = ttlp - (xnum-1.)*0.5*ss;
			facty = yaxlen/(float)nocy;
			symbol(sp-0.4,-0.5*ss,ss,"10",0.0,2);
			number(999.0,.5*ss-.06,.5*ss,y1,0.0,-1);
			plot(0.0,0.0,3);
			for(j=1;j<=nocy;j++){
				for(i=1;i<=10;i++){
					y = i;
					y = log10(y) * facty + (j-1)*facty;
					if(i != 1){
						plot(0.0,y,2);
						plot(-sst,y,2);
					}
					plot(0.0,y,3);
				}
				plot(-slt,y,2);
				symbol(sp-.4,y-.5*ss,ss,"10",0.0,2);
				ypo = ypo + 1;
				number(999.0,y+.5*ss-.06,.5*ss,ypo,0.0,-1);
				plot(0.0,y,3);
			}
			ytl=mty;
			ytl = (yaxlen-ytl*sttl)/2.0;
			symbol(ttlp-.2,ytl,sttl,ttly,90.,mty);
		}
	}
}
