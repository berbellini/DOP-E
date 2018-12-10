/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SHADEP                                                c
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

#include <stdlib.h>
#include <stddef.h>
#include "callib.h"
#include "dosubs.h"
#include "calplot.h"
extern struct Gcplot Gcplt;


void shadep(int narr, float *xarr, float *yarr)
{
	float xx, yy;
	INT i;
	INT n;
	INT ilw, jlw;
	INT *x, *y;
	if((x = (INT *)calloc(narr, sizeof(INT))) == NULL)
		return;
	if((y = (INT *)calloc(narr, sizeof(INT))) == NULL)
		return;
	n = narr;
	for(i = 0 ; i < n ; i++){
       		xx = 1000. * (xarr[i]*Gcplt.xstp + Gcplt.xold);
       		yy = 1000. * (yarr[i]*Gcplt.ystp + Gcplt.yold);
		if(xx > 1000000000.0)
			ilw = 1000000000;
		else if(xx < -1000000000.0)
			ilw = -1000000000;
		else
			ilw = xx ;
		if(yy > 1000000000.0)
			jlw = 1000000000;
		else if(yy < -1000000000.0)
			jlw = -1000000000;
		else
			jlw = yy ;
		x[i] = ilw;
		y[i] = jlw;
	}
	(*do_fillp)(n,x,y);
	free(x);
	free(y);
}
