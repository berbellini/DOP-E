/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: TSTCRS                                                c
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

#include <stdio.h>
#include "calplot.h"

main()
{
	char iarr[10];
	char outstr[40];
	float xx, yy, xxx, yyy;
	int i, ix, iy;

	ginitf("TSTCRS.PLT","TSTCRS");
	gunit("in");
	plot(1.0,1.0,-3);
	plot(2.0,2.0,2);
	plot(2.0,0.0,2);
	plot(0.0,2.0,2);
	plot(0.0,0.0,2);
	xx = 5.0;
	yy = 5.0;
	plus(xx,yy);
	for(i=1;i<=3;i++){
	if(i == 1){
		curixy(&ix,&iy,iarr);
		xx = ix/1000.0;
		yy = iy/1000.0;
	} else if(i == 2){
		curaxy(&xx,&yy,iarr);
	} else if(i == 3){
		factor(0.5);
		currxy(&xx,&yy,iarr);
	}
		sprintf(outstr,"%11.4e%11.4e%c",xx,yy,iarr[0]);
		xxx = 1.0;
		yyy = 4.0 - 0.5*i;
		gwrtxt(xxx,yyy,outstr,1);
		number(xx+0.40,yy+0.25,0.20,xx,0.0,3);
		number(xx+0.40,yy+0.00,0.20,yy,0.0,3);
		symbol(xx+0.40,yy-0.25,0.20,iarr,0.0,1);
	}
	pend();
}

plus(xx,yy)
float xx, yy;
{
		plot(xx,yy,3);
		plot(xx+0.2,yy,2);
		plot(xx-0.2,yy,2);
		plot(xx,yy+0.2,3);
		plot(xx,yy-0.2,2);
		plot(xx,yy,3);
}
