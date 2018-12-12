/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLOT                                                  c
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

void clospl(int mode);
extern struct Gcplot Gcplt;

void plot(float x,float y,int ipen)
{
/*
c-----
c       basic routine to move point on plotter
c
c       x       abscissa in inches
c       y       ordinate in inches
c       ipen    -3 pen up define a new origin
c               -2 pen down define a new origin after
c                  pen movement
c                2 pen down during movement
c                3 pen up during movement
c              999 terminate plotting after pen movement
c             1001 reset magnification - subroutine factor
c-----
c-----
c      xold and yold are absolute plotter coordinates
c      xcur and ycur are coordinates of point relative
c            to current plotter position in current
c            software units
c-----
*/
        INT ix,iy;
	float xn, yn, xx, yy;
/*
c-----
c	x and y are page coordinates in units given by
c	last call gunit()
c-----
*/
	if(ipen == 1001){
/*
c-----
c	implement factor, also reissue last width command
c-----
*/
		Gcplt.xstp=x;
		Gcplt.ystp=y;
		gwidth(Gcplt.owidth);
	} else {
		Gcplt.xcur=x;
		Gcplt.ycur=y;
		xn=Gcplt.xcur*Gcplt.xstp;
		yn=Gcplt.ycur*Gcplt.ystp;
/*
c-----
c	current position is in inches, ala calcomp
c	but to use unix plot filters we need integer*2
c	the factor 1000 is really 11000 counts / 11.0 inches
c-----
*/;
        	xx = 1000. * (xn + Gcplt.xold);
        	yy = 1000. * (yn + Gcplt.yold);
		if(xx > 1000000000.0)
			ix = 1000000000;
		else if(xx < -1000000000.0)
			ix = -1000000000;
		else
			ix = xx ;
		if(yy > 1000000000.0)
			iy = 1000000000;
		else if(yy < -1000000000.0)
			iy = -1000000000;
		else
			iy = yy ;
        	if(ipen == 3)(*do_move)(ix,iy);
        	if(ipen == 2)(*do_cont)(ix,iy);
        	if(ipen == -3)(*do_move)(ix,iy);
        	if(ipen == -2)(*do_cont)(ix,iy);
		if(ipen == 999)clospl(0);
		if(ipen < 0){
			Gcplt.xold=xn+Gcplt.xold;
			Gcplt.yold=yn+Gcplt.yold;
			Gcplt.xcur=0.0;
			Gcplt.ycur=0.0;
		}
	}
}
