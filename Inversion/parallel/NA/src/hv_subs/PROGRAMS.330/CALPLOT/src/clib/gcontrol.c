/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GCONTROL                                              c
c                                                                     c
c      COPYRIGHT (C)  2005 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      3507 Laclede Avenue                                            c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c-----
c      send information to low level graphics
c-----
*/
#include "calplot.h"
void 	gocontrol(int type, int i1, int i2, int i3, int i4);

void gcontrol(int type, float p1, float p2, float p3, float p4)
{
	int i1, i2, i3, i4;
	i1 = (int)(p1 *10000.0) ;
	i2 = (int)(p2 *10000.0) ;
	i3 = (int)(p3 *10000.0) ;
	i4 = (int)(p4 *10000.0) ;
	gocontrol( type,  i1, i2, i3, i4);
}

