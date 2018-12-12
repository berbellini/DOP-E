/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SHADET                                                c
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

void shadet(float x1,float y1,float x2,float y2,
	float x3,float y3,int ipatx,int ipaty,float xlen,float ylen)
{
/*
c-----
c	shade a triangular region
c-----
c	x1,y1	- coordinate of one corner of triangle
c	x2,y2	- coordinate of second corner of triangle
c	x3,y3	- coordinate of third corner of trianlge
c		  (one or more corners may be identical, e.g.,
c		   degenerates to a line or a point )
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
	INT jx1,jy1,jx2,jy2,jx3,jy3,jpatx,jpaty,jlenx,jleny;
	float arr[6];
	INT iarr[6];
	INT jarr[6];
	int i;
	arr[0] = x1;
	arr[1] = y1;
	arr[2] = x2;
	arr[3] = y2;
	arr[4] = x3;
	arr[5] = y3;
/*
c-----
c current position is in inches, ala calcomp
c but to use unix plot filters we need integer*2
c the factor 1000 is really 11000 counts / 11.0 inches
c-----
*/
	for(i=0;i<6;i++){
		if(i == 0 || i == 2 || i == 4){
			arr[i] = arr[i] * Gcplt.xstp;
			arr[i] = 1000.0*(arr[i] + Gcplt.xold);
			iarr[i] = arr[i];
		} else {
			arr[i] = arr[i] * Gcplt.ystp;
			arr[i] = 1000.0*(arr[i] + Gcplt.yold);
			iarr[i] = arr[i];
		}
		if(iarr[i] > 100000000)iarr[i]=100000000;
		if(iarr[i] < -100000000)iarr[i]=-100000000;
		jarr[i] = iarr[i];
	}
	jx1 = jarr[0];
	jy1 = jarr[1];
	jx2 = jarr[2];
	jy2 = jarr[3];
	jx3 = jarr[4];
	jy3 = jarr[5];
	jpatx = ipatx;
	jpaty = ipaty;
	jlenx = ( 1000.0 * Gcplt.xstp * xlen);
	jleny = ( 1000.0 * Gcplt.ystp * ylen);
/*
c-----
c	protect against later divide errors
c-----
*/
	if(jlenx <= 0)jlenx = 1;
	if(jleny <= 0)jleny = 1;
	(*do_fillt)(jx1,jy1,jx2,jy2,jx3,jy3,jpatx,jpaty,jlenx,jleny);
}
