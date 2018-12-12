/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: TRITST                                                c
c                                                                     c
c      COPYRIGHT (C)  1986 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
*/
#include "calplot.h"
main()
{
	int ipatx, ipaty;
	float xlen, ylen;
	pinitf("TRITST.PLT");
	gunit("in");
	ipatx = 4;
	ipaty = 4;
	xlen = 0.05;
	ylen = 0.05;
	shadet(0.0,0.0,0.0,1.0,1.0,1.0,ipatx,ipaty,xlen,ylen);
	plot(1.0,1.0,-3);
	factor(0.5);
	shadet(0.0,0.0,0.0,1.0,1.0,1.0,ipatx,ipaty,xlen,ylen);
	factor(1.0);
	plot(1.0,0.0,-3);
	ipatx = 0;
	ipaty = 7;
	shadet(0.0,0.0,0.0,1.0,1.0,1.0,ipatx,ipaty,xlen,ylen);
	plot(-2.0,-1.0,-3);
	plot(0.0,2.0,-3);
	factor(2.0);
	ipatx = 0;
	shadet(0.0,0.0,0.0,1.0,1.0,1.0,ipatx,ipaty,xlen,ylen);
	pend();
}
