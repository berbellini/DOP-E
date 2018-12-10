/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: LINECLIP                                              c
c                                                                     c
c      COPYRIGHT (C)  1986, 1989 R. B. Herrmann                       c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      SaINT Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
*/

	
/* line clipping algorithm for a rectangular area 			*/
/* Principles of Interactive Computer Graphics, W. Newman and		*/
/*	R. Sproull, 1979, McGraw-Hill pp 66-67				*/
/*									*/
/* Computer Graphics with Pascal, M. Berger, 1986,			*/
/*	Benjamin/Cummings pp 74-77 (Cohen-Sutherland Algorithm)		*/
/* returns 1 if the line is not plottable, otherwise returns corrected	*/
/*	coordinates							*/
/* assume coordinate axes are such that left < x < right		*/
/*				      bottom < y < top			*/
#include <stdio.h>
#ifdef MSDOS
#define INT long
#else
#define INT int
#endif

static INT	Xmin = 01;
static INT	Xmax = 02;
static INT	Ymax = 010;
static INT	Ymin = 04;
static INT	None = 00;

static void linecode(INT x,INT y,INT MaxY,INT MinX,INT MinY,INT MaxX,INT *c);

INT dv_lineclip(INT *x0,INT *z0,INT *x1,INT *z1,
	INT MinX,INT MinY,INT MaxX,INT MaxY)
{
	INT xx0, yz0, xx1, yz1;
	INT c, c0, c1;
	INT x,y;
	double slope, slopeinv;
	xx0 = *x0;
	xx1 = *x1;
	yz0 = *z0;
	yz1 = *z1;
	linecode(xx0,yz0,MaxY,MinX,MinY,MaxX,&c0);
	linecode(xx1,yz1,MaxY,MinX,MinY,MaxX,&c1);
	if( xx0 != xx1)
		slope = (double)(yz1-yz0)/(double)(xx1-xx0);
	if( yz0 != yz1)
		slopeinv = (double)(xx1-xx0)/(double)(yz1-yz0);
	while( (c0 != None) || (c1 != None) ){
		if( (c0&c1) != None)
			return (1);
		if(c0 == c1)
			return (1);
		if(c0 == None){
			c = c1;
		} else {
			c = c0;
		}
		if(c & Xmin){
			y = yz0 + (INT)(slope*(double)(MinX-xx0));
			x = MinX;
		}
		if( c & Xmax){
			y = yz0 + (INT)(slope*(double)(MaxX-xx0));
			x = MaxX;
		}
		if( c & Ymax){
			x = (INT)(slopeinv*(double)(MaxY-yz0)) + xx0;
			y = MaxY;
			}
		if( c & Ymin){
			x = (INT)(slopeinv*(double)(MinY-yz0)) + xx0;
			y = MinY;
		}	
		if(c == c0){
			xx0 = x;
			yz0 = y;
			linecode(xx0,yz0,MaxY,MinX,MinY,MaxX,&c0);
		} else {
			xx1 = x;
			yz1 = y;
			linecode(xx1,yz1,MaxY,MinX,MinY,MaxX,&c1);
		}
	}
	*x0 = xx0;
	*x1 = xx1;
	*z0 = yz0;
	*z1 = yz1;
	return (0);
}

static void linecode(INT x,INT y,INT MaxY,INT MinX,INT MinY,INT MaxX,INT *c)
{
	*c = None;
	if( x < MinX)
		*c = Xmin;
	else if(x > MaxX)
		*c = Xmax;
	if( y < MinY)
		*c |= Ymin;
	else if( y > MaxY)
		*c |= Ymax;
}
