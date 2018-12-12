/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLTSCL                                                c
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
#include <math.h>
#include "calplot.h"

extern struct Gcplot Gcplt;

#define LN10 2.302585093

void pltscl(float x[],float axlen,int n,float *x1,float *deltax,int nocx)
{
/*
c-----
c	PLTSCL
c		GIVEN ARRAY, ESTABLISHES SCALING PARAMETERS
c		FOR LOG scale, sets according to largest value
c	
c		x	-	array of values
c		axlen	-	length of axis in inches
c		n	-	number of points
c		x1	-	value of first point of axis
c				returned value
c		deltax 	-	inches per cycle if log
c				units per inch is linear
c		nocx	-	positive - number of cycles on
c				log scale
c
c				zero or negative - linear scale
c-----
*/
	float xmax;
	int i, lmax;
	if(nocx > 0){
		xmax = -1.0e+38;
		for(i=0;i< n;i++)
			if(x[i] > xmax)xmax = x[i];
		if(xmax > 0.0){
			xmax = log(xmax)/LN10;
			if(xmax > 0.0)
				lmax = xmax +1;
			else
				lmax = xmax;
			*x1 = lmax;
		} else {
			*x1 = 0.0;
		}
		*deltax = axlen/nocx;
		*x1 = *x1 - nocx;
	} else {
		gscale(x,axlen,n,1);
		*x1 = x[n];
		*deltax = x[n+1];
	}
}
