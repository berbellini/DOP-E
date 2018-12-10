/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SYMBOL                                                c
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
#include <string.h>
#include "callib.h"
#include "dosubs.h"
#include "calplot.h"

extern struct Gcplot Gcplt;

void symbol(float xloc,float yloc,float height,
	char *inbuf,float angle,int nocar) 
{

/*
c-----
c       produces either 1) one symbol of choice or 2) character string
c
c       xloc    x-coordinate lower left corner of symbol
c               if nocar < 0 and symbol is centered, center x point
c       yloc    y-coordinate lower left corner of symbol or string
c       height  height in inches of symbol
c       inbuf   nocar > 0 inbuf is ascii string of characters
c               nocar long
c               nocar < 0 inbuf is a single character of an
c                         integer equivalent, e.g., char(inteq)
c
c                         of the symbol to be plotted
c       angle   angle of rotation of symbol in degrees
c       nocar   > 0 number of characters in inbuf string
c               -1 pen is up during move after which single symbol
c                  is plotted
c               -2 pen is down, a line is drawn from present position
c                  to point where symbol is drawn
c-----
*/
	INT ix, iy, iht, iang, nchar;
	float xn, yn, xx, yy, xht, ang;

		if(xloc  <  999.0) Gcplt.x0 = xloc;
		if(yloc  <  999.0) Gcplt.y0 = yloc;
		if(nocar  >  0  &&  strlen(inbuf) < nocar){
			nchar = strlen(inbuf);
		} else {
 			nchar = nocar;
		}
		Gcplt.xcur=Gcplt.x0;
		Gcplt.ycur=Gcplt.y0;
/*
c-----
c	convert into inches if necessary
c-----
*/
		xn=Gcplt.xcur*Gcplt.xstp;
		yn=Gcplt.ycur*Gcplt.ystp;
/*
c-----
c	current position is in inches, ala calcomp
c	but to use unix plot filters we need integer*2
c	the factor 1000 is really 11000 counts / 11.0 inches
c-----
*/
		xx = 1000.0 * (xn + Gcplt.xold);
		yy = 1000.0 * (yn + Gcplt.yold);
		if(xx > 1000000000.0)
			ix = 1000000000;
		else if(xx < -1000000000.0)
			ix = -1000000000;
		else 
			ix = xx;
		
		if(yy > 1000000000.0)
			iy = 1000000000;
		else if(yy < -1000000000.0)
			iy = -1000000000;
		else 
			iy = yy;
		
/*
c-----
c	convert height to proper units
c-----
*/
		xht = (1000.0 * height * Gcplt.xstp);
		if(xht < 0.0){
			iht = 0;
		} else if(xht > 1000000000.0){
			iht = 1000000000;
		} else {
			iht = xht;
		}
		iang = angle;
		(*do_gsymb)(ix,iy,iht,iang,nchar,inbuf);
/*
c-----
c	adjust current position
c-----
*/
		if(nocar  >  -1){
			ang = angle * 3.1415927/180.0;
			Gcplt.x0 = Gcplt.x0 + height * cos(ang) * nocar;
			Gcplt.y0 = Gcplt.y0 + height * sin(ang) * nocar;
			Gcplt.xcur=Gcplt.x0;
			Gcplt.ycur=Gcplt.y0;
		}
}
