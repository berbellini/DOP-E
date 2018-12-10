/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: LINED                                                 c
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

#define SYMSIZ 0.07
#define MAX(a,b) ( (a) > (b) ? (a):(b) )

void lined(float x[],float y[],int n,int inc,int lintyp,
		int inteq,int ipat,float xlen)
{
/*
c-----
c       x       array of abscissa
c       y       array of ordinates
c       n       number of points to be plotted
c       inc     plot every inc'th point
c               (there are n*inc + inc + 1 points in array)
c       lintyp  != 0 plot symbol inteq at each abs(lintyp) point
c               =0   plot line only
c               <0   plot only symbols, no connecting line
c               >0   line connects symbol
c       inteq   if lintyp != 0, inteq symbol plotted at point
c	ipat	integer providing bit pattern for dashed lines
c	xlen	length of each on bit in pattern ipat
c-----
*/
	int dj, j1, j2, l1, l2, l3, ipen, icode, ia, m;
	float cx, cy, cf, dx, dy, x1, y1, xn, yn, dl, dr;
	float sym = SYMSIZ;
	char cstr[2];
	cstr[1] = '\0';
/*
c-----
c       start plotting from end of line closest to current pen position
c-----
*/
	l1 = n*inc ;
	l2 = l1 + inc;
	l3 = l1 - inc;
	x1 = x[l1];
	dx = x[l2];
	y1 = y[l1];
	dy = y[l2];
/*
c-----
c       determine current pen position
c-----
*/
	where(&cx,&cy,&cf);
/*
c-----
c       determine which end of plot array is closest to current pen
c       position
c-----
*/
	dl = MAX(fabs((x[0] -x1)/dx-cx), fabs((y[0] -y1)/dy-cy));
	dr = MAX(fabs((x[l3]-x1)/dx-cx), fabs((y[l3]-y1)/dy-cy));
        if(dr < dl){
                dj = -inc;
                j1 = l3;
                j2 = 0;
        } else {
                dj = inc;
                j1 = 0;
                j2 = l3;
        }
/*
c-----
c       do the plotting
c-----
*/
	ia=lintyp;
	if(ia < 0)
		ia = - ia;
	ipen = 3;
	icode = -1;
	cstr[0] = (char)inteq;
	m=j1 -dj;
	while((m<j2 && dj > 0) || (m>j2 && dj < 0)){ 
		m+=dj;
		xn = (x[m]-x1)/dx;
		yn = (y[m]-y1)/dy;
		if( ipat == 0 || ipen == 3){
			plot(xn,yn,ipen);
		} else {
			plotd(xn,yn,ipat,xlen);
		}
		if(ia > 0){
			if((m%ia)==0 || ia==1)symbol(xn,yn,sym,cstr,0.0,icode);
                }
		if(lintyp >= 0)ipen=2;
	}
/*
c-----
c       penup
c-----
*/
        plot(xn,yn,3);
}
