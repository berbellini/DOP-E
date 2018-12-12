/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: TESTL                                                 c
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
c----
c	program to test the scale and line routines
c----
*/
#include <math.h>
#include "calplot.h"

main()
{
	float	xaxlen = 5.0,
		yaxlen = 1.0;
	float x[120],c[120],s[120], cc[120],ss[120];
	float xx, yy, ang;
	int i, j;
	xx = xaxlen;
	yy = yaxlen;
	for(i=0;i<100;i++){
		ang=0.01*6.2831853*(i);
		x[i]=ang;
		c[i]=1234.0*cos(ang);
		cc[i]=0.015 * cos(ang);
		ss[i]=0.015 * cos(ang);
		s[i]=1234.0*sin(ang);
	}
	pinitf("TESTL.PLT");
	gunit("in");
	plot(0.0,0.0,-3);
	for(j=1;j<=3;j++){
		if(j == 1){
			plot(1.0,5.0,-3);
			gscale(x,xaxlen,100,1);
			gscale(c,yaxlen,100,1);
			gscale(s,yaxlen,100,1);
			axis(0.0,0.0,"X-AXIS",-6,xaxlen,0.0,x[100],x[101]);
			axis(0.0,0.0,"Y-AXIS",6,yaxlen,90.0,s[100],s[101]);
			line(x,c,100,1,0,0);
			line(x,s,100,1,-10,1);
		} else if(j == 2){
			plot(0.0,-2.0,-3);
			gscale(x,xaxlen,50,-2);
			gscale(c,yaxlen,50,2);
			gscale(s,yaxlen,50,2);
			axis(0.0,0.0,"X-AXIS",-6,xaxlen,0.0,x[100],x[102]);
			axis(0.0,0.0,"Y-AXIS",6,yaxlen,90.0,s[100],s[102]);
			line(x,c,50,2,-2,2);
			line(x,s,50,2,0,0);
		} else if(j == 3){
			plot(0.0,-2.0,-3);
			gscale(x,xaxlen,100,1);
			gscale(cc,yaxlen,100,1);
			gscale(ss,yaxlen,100,1);
			axis(0.0,0.0,"X-AXIS",-6,xaxlen,0.0,x[100],x[101]);
			axis(0.0,0.0,"Y-AXIS",6,yaxlen,90.0,ss[100],ss[101]);
			line(x,cc,100,1,0,0);
			line(x,ss,100,1,-10,1);
		}
	}
	pend();
}
