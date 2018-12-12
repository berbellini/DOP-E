/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLTTST                                                c
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
*/

#include <math.h>
#include "calplot.h"

main()
{
	float five = 5.0;
	float half = 0.5;
	float height = 0.25;
	float zero = 0.0;
	float eight = 8.0;
	float ht = 0.10;
	int one = 1;
	int two = 2;
	int three = 3;
	int i, icolor, j, ipen;
	float rmax, xx, yy, ang, xmax, ymax;
	
	pinitf("PLTTST.PLT");
	gunit("in");
	symbol(five,half,height,"X",zero,one);
	symbol(half,five,height,"Y",zero,one);
	for(i=1;i<=40;i++){
		icolor = i;
		newpen(icolor);
		rmax= i*0.5;
		number(height,rmax,ht,rmax,zero,one);
		number(rmax,height,ht,rmax,zero,one);
		for(j=0;j<=90;j++){
			ang = 3.1415927*(float)j/180.0;
			xx = rmax*cos(ang);
			yy= rmax*sin(ang);
			if(j > 0){
				ipen=2;
			} else {
				ipen=3;
			}
			plot(xx,yy,ipen);
		}
	}
	plot(zero,zero,three);
	plot(eight,zero,two);
	plot(eight,eight,two);
	plot(zero,eight,two);
	plot(zero,zero,two);
	frame();
	newpen(one);
	for(i=1;i<=40;i++){
		xmax= i*0.5;
		ymax=xmax;
		number(height,ymax,ht,ymax,zero,one);
		number(xmax,height,ht,xmax,zero,one);
		plot(zero,zero,3);
		plot(xmax,zero,2);
		plot(xmax,ymax,2);
		plot(zero,ymax,2);
		plot(zero,zero,2);
	}
	symbol(five,half,height,"X",zero,one);
	symbol(half,five,height,"Y",zero,one);
	pend();
}
