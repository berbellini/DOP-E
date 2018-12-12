/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SHADER                                                c
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
#include "dosubs.h"
#include "calplot.h"


extern struct Gcplot Gcplt;

void shader(float x1,float y1,float x2,float y2,int ipatx,
	int ipaty,float xlen,float ylen)
{
/*
c-----
c	shade a rectangular region
c-----
c	x1,y1	- coordinate of one corner of triangle
c	x2,y2	- coordinate of opposite corner of rectangle
c	ipatx	- pattern for shading used in drawing lines 
c		  parallel to x-axis
c	ipaty	- pattern for shading used in drawing lines
c		  parallel to y-axis
c		  (patterns are integer from 0 to 31)
c	xlen	- length in inches of one bit of pattern
c	ylen	- length in inches of one bit of pattern for y-axis
c		  (see discussion of plotd, etc)
c-----
c-----
c      xold and yold are absolute plotter coordinates
c      xcur and ycur are coordinates of point relative
c            to current plotter position in current
c            software units
c-----
*/
	INT jx1,jy1,jx2,jy2,jpatx,jpaty,jlenx,jleny;
	float arr[4];
	INT iarr[4];
	INT jarr[4];
	int i;
	arr[0] = x1;
	arr[1] = y1;
	arr[2] = x2;
	arr[3] = y2;
/*
c-----
c current position is in inches, ala calcomp
c but to use unix plot filters we need integer*2
c the factor 1000 is really 11000 counts / 11.0 inches
c-----
*/
	for(i=0;i<4;i++){
		if(i == 0 || i ==2 ){
			arr[i] = arr[i] * Gcplt.xstp;
			arr[i] = 1000.0*(arr[i] + Gcplt.xold);
			iarr[i] = arr[i];
		} else {
			arr[i] = arr[i] * Gcplt.ystp;
			arr[i] = 1000.0*(arr[i] + Gcplt.yold);
			iarr[i] = arr[i];
		}
		if(iarr[i] > 1000000000)iarr[i]=1000000000;
		if(iarr[i] < -1000000000)iarr[i]=-1000000000;
		jarr[i] = iarr[i];
	}
	jx1 = jarr[0];
	jy1 = jarr[1];
	jx2 = jarr[2];
	jy2 = jarr[3];
	jpatx = ipatx;
	jpaty = ipaty;
	jlenx = ( 1000.0 * Gcplt.xstp * xlen);
	jleny = ( 1000.0 * Gcplt.ystp * ylen);
	if(jlenx <= 0)jlenx = 1;
	if(jleny <= 0)jleny = 1;
	(*do_fillr)(jx1,jy1,jx2,jy2,jpatx,jpaty,jlenx,jleny);
}
