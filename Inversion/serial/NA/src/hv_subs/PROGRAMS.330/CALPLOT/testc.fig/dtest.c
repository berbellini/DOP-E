/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: DTEST                                                 c
c                                                                     c
c      COPYRIGHT (C) 1987 R. B. Herrmann                              c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c-----
c	program to test dashed line routines
c-----
*/

#include <math.h>
#include "calplot.h"

main()
{
	int i, ipage, ipat, k, inc;
	float xlen, ylen, xx, yy, xi, x0, y0, rad, ang;
	pinitf("DTEST.PLT");
	gunit("in");
	plot(0.0,0.0,-3);
/*
c-----
c	test line algorithms
c-----
*/
	for(ipage =1 ; ipage <=3; ipage++){
		if(ipage == 1){
			xlen = 0.05;
		} else if(ipage == 2){
			xlen = 0.1;
		} else {
			xlen = 0.2;
		}
		symbol(1.4,8.6,0.14,"TEST OF PLOTD: XLEN= ",0.0,21);
		number(999.,999.,0.14,xlen,0.0,2);
		for(ipat=0;ipat <=31;ipat++){
			yy = 8.5 - 0.2*ipat;
			xi = ipat;
			number(1.0,yy-0.05,0.10,xi,0.0,-1);
			plot(1.4,yy,3);
			plotd(7.5,yy,ipat,xlen);
		}
/*
c-----
c		test ability to get correct pattern for a line
c		composed of many small segments
c-----
*/
		for(i=0;i<=61;i++){
			yy = 2.0;
			xx = 1.4 +  i*0.1;
			if(i == 0){
				plot(xx,yy,3);
			} else {
				plotd(xx,yy,31,xlen);
			}
		}
		axis(1.4,1.6,"LENGTH",-6,6.0,0.0,0.0,1.0);
		if(ipage != 3)
			frame();

	}
/*
c-----
c	draw a sequence of polygons to test algorithm
c	ability to correctly dash non-vertical and non-horizontal
c	line segments
c-----
*/
	frame();
	for(k=0;k<=31;k+=3){
		rad =  (k+1)*0.1;
		x0 = 4.0;
		y0 = 5.0;
		xlen = 0.05;
		inc = 15;
		ipat = k;
		for(i=0;i<=360;i+=inc){
			ang = i*3.1415927/180.0;
			xx = x0 + rad*cos(ang);
			yy = y0 + rad*sin(ang);
			if(i == 0)plot(xx,yy,3);
			plotd(xx,yy,ipat,xlen);
		}
	}
	pend();
}
