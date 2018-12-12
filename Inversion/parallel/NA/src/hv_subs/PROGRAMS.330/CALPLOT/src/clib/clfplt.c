/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: CLFPLT                                                c
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

#include "sysinit.h"
#include "callib.h"
#include "calplot.h"

extern struct Gcplot Gcplt;

/*
c-----
c	dummy subroutines for cursor interaction
c	returns (0,0) always
c-----
*/
void curixy(int *ix,int *iy,char *ic)
{
/*
c-----
c	returns cursor position in absolute integer coordinates
c-----
*/
	cross(ix,iy,ic);
}

void curaxy(float *xx,float *yy,char *ic)
{
/*
c-----
c	returns coordinates in inches of cursor position
c	with respect to the absolute origin
c-----
*/
	int ix, iy;
	curixy(&ix,&iy,ic);
		*xx = ix/1000.0;
		*yy = iy/1000.0;
		if(Gcplt.iunit == 1){
			*xx = (*xx) *2.54;
			*yy = (*yy) *2.54;
		}
}

void currxy(float *xx,float *yy,char *ic)
{
/*
c-----
c	returns cursor position relative to current origin,
c	taking into account any user scaling ,e.g., call factor
c-----
*/
	float ax, ay;
	curaxy(&ax,&ay,ic);
	if(Gcplt.iunit == 1){
		ax = ax * 0.3937;
		ay = ay * 0.3937;
	}
	*xx = (ax-Gcplt.xold)/Gcplt.xstp;
	*yy = (ay-Gcplt.yold)/Gcplt.xstp;
}

