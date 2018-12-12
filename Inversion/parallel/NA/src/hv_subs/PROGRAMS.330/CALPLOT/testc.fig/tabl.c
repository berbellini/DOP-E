/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: TABL                                                  c
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
c
c     this routine produces a symbol table, which shows the characters
c     available in the symbol routine.
c-----
*/
#include "calplot.h"
float znum[4] = {
	10293.84756,0.019375,-204.86,-12345.6789};

void rect(float llx, float lly, float urx, float ury);

main()
{
	float z, xs, ys, x, y;
	char cc[2];
	int m, ia, ib, ndec, mdec;
	cc[1]='\0';
		pinitf("TABL.PLT");
		gunit("in");
		rect(0.4,0.4,9.6,7.5);
		rect(0.6,0.6,9.4,7.1);
		gfont(1);
		symbol(2.15,7.25,0.15, 
			"characters available in symbol routine" ,
			0.0,38);
		z=0.0;
		m=0;
		x=0.6;
		y=6.57;
		xs=0.75;
		ys=0.25;
		gfont(0);
		for(ia=0;ia <8;ia ++){
			for(ib=0;ib<14;ib++){
				number(x+.1,y+.18,.10,z,0.0,-1);
				cc[0] = (char)m;
				symbol(x+xs,y+ys ,.30 ,cc,0.0,-1);
				z=z+1.0;
				m=m+1;
				y=y-0.45;
			}
			x=x+1.1;
			plot(x,0.6,3);
			plot(x,7.1,2);
			y=6.57;
			xs=.55;
			ys=.05;
		}
		gfont(1);
		symbol(2.05,0.45,.10,
		"integer for use in symbol call shown to left of each symbol"
				,0.0,59);
/* begin a new page */
		frame();
/*
c-----
c	the following tests the number subroutine for precision
c-----
*/
		symbol(1.0,2.5,.14,
			"example of number subroutine",90.0,28 );
			symbol(2.62,7.5,.14,
				"call number(xx,yy,ht,fpn,ang,ndec)",
				0.0,34);
		symbol( 1.5,7.2,.10,"ndec",0.0,4);
		symbol( 2.0,7.2,.10,"number",0.0,6);
		symbol( 4.5,7.2,.10,"ndec",0.0,4);
		symbol( 5.0,7.2,.10,"number",0.0,6);
		symbol( 7.5,7.2,.10,"ndec",0.0,4);
		symbol( 8.0,7.2,.10,"number",0.0,6);
		y=7.0;
		for(ia=0;ia<4;ia++){
			for(ib=1;ib<=11;ib++){
				ndec = ib-6;
				number(1.5,y,.07,(float)ndec,0.0,-1);
				number(2.0,y,.07,znum[ia],0.0,ndec);

				ndec = 999 + ib;
				mdec = ndec + 1000;
				number(4.5,y,.07,(float)ndec,0.0,-1);
				number(5.0,y,.07,znum[ia],0.0,ndec);
				number(7.5,y,.07,(float)mdec,0.0,-1);
				number(8.0,y,.07,znum[ia],0.0,mdec);
				y=y-0.14;
			}
			y=y-0.2;
		}
	pend();
}

void rect(float llx, float lly, float urx, float ury)
{
	plot(llx,lly,3);
	plot(urx,lly,2);
	plot(urx,ury,2);
	plot(llx,ury,2);
	plot(llx,lly,2);
	plot(llx,lly,3);
}
