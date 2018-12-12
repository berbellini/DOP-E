/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: TSTCUR                                                c
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

#include <math.h>
#include "calplot.h"


main()
{
	float x[22],y[22];
	char	ic[2], is[2];
	char	ttlx[10],ttly[10];
	char 	outstr[40];
	float fct, xaxlen, yaxlen, xx, yy, x1, y1, deltax, deltay, ht;
	int i, n;
	int nocx, nocy, inteq, lintyp, mtx, mty;

	ginitf("INTER","TSTCUR");
	gunit("in");
	n = 20;
	ht = 0.07;
	for(i=1;i<=n;i++){
		x[i-1] = i;
		y[i-1] = i*i*i;
	}
	for(i=1;i<=4;i++){
		plot(2.0,1.0,-3);
		if(i == 1){
			nocx=2;
			nocy=3;
			fct = 1.0;
			inteq = 0;
			lintyp = 0;
		} else if(i == 2){
			nocx=0;
			nocy=2;
			fct = 0.8;
			inteq = i;
			lintyp = -1;
		} else if(i == 3){
			nocy=0;
			nocx=2;
			fct = 1.0;
			inteq = 1;
			lintyp =  1;
		} else if(i == 4){
			nocx=0;
			nocy=0;
			fct = 0.8;
			inteq = i;
			lintyp = 0;
		}
		is[0] = (char)inteq;
		factor(fct);
		yaxlen=6.0;
		xaxlen=5.0;
		pltscl(x,xaxlen,n,&x1,&deltax,nocx);
		pltscl(y,yaxlen,n,&y1,&deltay,nocy);
		strcpy(ttlx , "X-AXIS");
		mtx = 6;
		strcpy(ttly , "Y-AXIS");
		mty = 6;
		algaxe(xaxlen,yaxlen,nocx,nocy,ttlx,ttly,mtx,mty,
			x1,y1,deltax,deltay);
		if(i <= 2){
			pltlog(x,y,n,x1,y1,deltax,deltay,lintyp,inteq,ht,
				nocx,nocy);
		} else {
			pltlgd(x,y,n,x1,y1,deltax,deltay,lintyp,inteq,ht,
				nocx,nocy,8*i -1 , 0.10);
		}
		curuxy(&xx,&yy,x1,y1,deltax,deltay,nocx,nocy,ic);
		sprintf(outstr,"%10.3f%10.3f %s",xx,yy,ic);
		gwrtxt(0.0,5.0,outstr,0);
		factor(1.0);
		plot(-2.0,-1.0,-3);
		frame();
	}
	pend();
}
