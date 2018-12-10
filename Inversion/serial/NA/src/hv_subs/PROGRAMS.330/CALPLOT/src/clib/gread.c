/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      ROUTINE: GREAD                                               c
c                                                                     c
c      COPYRIGHT (C)  1997 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
*/
#include	<stdio.h>
#include	<string.h>
#include "callib.h"
#include "dosubs.h"
#include "calplot.h"
extern struct Gcplot Gcplt;

static void doconv(float xn, float yn, INT *ix, INT *iy);
void dv_gread(int lf,char *fname, INT NumX, INT NumY, INT LowX, INT LowY, INT HighX, INT HighY, INT Num, INT Sclx, INT Scly);

/*
	Read the PageNum'th page from the file fname.
	Clip the plot page about the limits (Xlow,Ylow) -> (Xhigh,Yhigh)
	The origin is redefined to be the lower left corner of the
	clip region. This corner is then shifted to
	position (X0, Y0)
*/
void gread(char *fname, float X0, float Y0, float Xlow, float Ylow, 
	float Xhigh, float Yhigh, int PageNum, float tsclx, float tscly)
{
	INT NumX, NumY, LowX, LowY, HighX, HighY, Num;
	int lf;
	INT Sclx;
	INT Scly;
	Num = (INT)PageNum;
	Sclx = (INT)(1000.0*tsclx);
	Scly = (INT)(1000.0*tscly);
	doconv(X0-Xlow*tsclx, Y0-Ylow*tscly, &NumX, &NumY);
	doconv(Xlow, Ylow, &LowX, &LowY);
	doconv(Xhigh, Yhigh, &HighX, &HighY);
	lf = strlen(fname);
	dv_gread(lf,fname, (INT)NumX, (INT)NumY, (INT)LowX, (INT)LowY, 
		(INT)HighX, (INT)HighY, (INT)Num, (INT)Sclx, (INT)Scly);
}

static void doconv(float xn, float yn, INT *ix, INT *iy);
static void doconv(float xn, float yn, INT *ix, INT *iy){
		float xx, yy;
        	xx = 1000. * (xn + Gcplt.xold);
        	yy = 1000. * (yn + Gcplt.yold);
		if(xx > 1000000000.0)
			*ix = 1000000000;
		else if(xx < -1000000000.0)
			*ix = -1000000000;
		else
			*ix = xx ;
		if(yy > 1000000000.0)
			*iy = 1000000000;
		else if(yy < -1000000000.0)
			*iy = -1000000000;
		else
			*iy = yy ;
}
