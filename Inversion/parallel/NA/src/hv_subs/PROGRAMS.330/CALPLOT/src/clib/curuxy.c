/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: CURUXY                                                c
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


void curuxy(float *xx, float *yy, float x1, float y1, 
	float deltax, float deltay,int nocx,int nocy,char *ic)
{
/*
c-----
c	subroutine to return user coordinates from cursor pick
c	this must be used in conjunction with algaxe.f, pltlog.f,
c	pltscl.f for proper mapping
c
c	xx	- user x-coordinate of cursor
c	yy	- user y-coordinate of cursor
c	ic	- character returned with cursor coordinates
c
c	x1	- user value of first point plotted at
c		  coordinate x = 0
c	y1	- user value of first y point plotted at
c		  coordinate y = 0
c	nocx	- >0 number of log cycles on x-axis
c		  <=0 linear x-axis
c	nocy	- >0 number of log cycles on y-axis
c		- <=0 linear y-axis
c	deltax  - units per inch linear plot on x-axis
c		- inches per cycle on logarithmic x-axis
c	deltay  - units per inch linear plot on y-axis
c		- inches per cycle on logarithmic y-axis
c-----
c
c	get coordinates relative to current origin and with current
c	scaling in effect through factor() call
c
c-----
*/
		currxy(xx,yy,ic);
/*
c-----
c	convert to user coordinates from screen coordinates
c-----
c	x - coordinate
c-----
*/
		if(nocx > 0)
			*xx = pow(10.0,(x1 + (*xx)/deltax));
		else
			*xx = x1 + (*xx)*deltax;
/*
c-----
c	y - coordinate
c-----
*/
		if(nocy > 0)
			*yy = pow(10.0,(y1 + (*yy)/deltay));
		else
			*yy = y1 + (*yy)*deltay;

}
