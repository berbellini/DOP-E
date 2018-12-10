/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SCALE                                                 c
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
#include <stdlib.h>
#define LN10 2.302585093

static float saver[7] = {
		1.0,2.0,4.0,5.0,8.0,10.0,20.0};

void gscale(float x[],float xaxlen,int n,int inc)
{
/*
c-----
c       purpose of scale is to find extreme values of array x
c       and set parameters for a nice plot by a combined use
c       of axis and line.
c-----
c       x       array to be scaled dimension n*inc+inc +1
c       xaxlen  axis length in inches into which data must fit
c       n       total number of points to be fit, spaced inc apart
c       inc     in the original array, each |inc|'th point is
c               used for scaling. Others are ignored. If array is
c               not well bounded, parameter will not be the best
c               but will suffice for the actual points used
c               The idea is to use the line routine to get a fast
c               plot
c               inc < 1, the first value is assumed to be a maximum
c-----
c       The object of scale is to define two parameters firstx and 
c       deltax
c       inc > 0 firstx is a min and deltax is positive
c       inc < 0 firstx is a max and deltax is negative
c-----
*/
	int k, nn, i, is;
	float xmin, xmax, tmax, firstv, deltav, 
		xscal, xlog, p, range, xmid, xx;
/*
c-----
c       for a nice plot we want scale to start at some nice number
c       raised to some power
c-----
*/
        k = abs(inc);
        nn=n*k;
/*
c-----
c       get largest and smallest numbers in the array
c-----
*/
        xmin = x[1];
        xmax=xmin;
	for(i=0;i<nn;i+=k){
                xx=x[i];
                if(xx > xmax)xmax=xx;
                if(xx < xmin)xmin=xx;
	}
        firstv = xmin;
        deltav = (xmax-xmin)/xaxlen;
/*
c-----
c       make deltav a multiple of one of the saver(i)
c-----
*/
        if(deltav <= 0.0){
                deltav=2.0*firstv;
                if(deltav <= 0.0)deltav=1.0;
		tmax = firstv + deltav;
        } else {
                xlog=log(deltav)/LN10;
                i = xlog + 500;
/*
c-----
c assume computer never has floating point number outside the
c range > 10**500 or < 10**(-500)
c-----
*/
                p=pow(10.0,(double)(i-500));
                deltav=deltav/p - 0.01;
                is=0;
		for(i=0;i<6;i++){
                        if(saver[i] < deltav)is=i;
		}
/*
c-----
c       construct final deltav
c-----
*/
s1000:
		is = is + 1;
		if(is > 6)goto s1001;
                deltav = saver[is]*p;
                range = (int)(xaxlen+0.49)*deltav;
                xmid= 0.5*(xmin+xmax);
		if(xmid < 0.0)xmid = xmid - deltav/2.0;
/*
c-----
c       find scaled element nearest the center
c-----
*/
                xscal = deltav*(int)(xmid/(0.99*deltav));
 		firstv= xscal - 0.5*deltav*(int)(xaxlen+0.49);
/*
c-----
c       test total range
c-----
*/
                tmax = firstv + range;
		if(tmax  <  1.01*xmax || firstv > 0.99*xmin) goto s1000;
                if(xmin*firstv <= 0.0)firstv=0.0;
	}
s1001:
/*
c-----
c       invoke sign of inc
c-----
*/
        if(inc > 0){
                x[nn]=firstv;
                x[nn+k]=deltav;
        } else {
                x[nn]=tmax;
                x[nn+k]= -deltav;
        }
}
