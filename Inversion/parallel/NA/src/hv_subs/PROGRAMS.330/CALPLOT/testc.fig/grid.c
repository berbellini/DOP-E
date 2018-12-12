#include <stdio.h>
#include <math.h>
#include "calplot.h"

#define NARR 16

float xarr[NARR] ;
float yarr[NARR] ;

void box(float xl, float yl, float xh, float yh);

char strout[100];

main()
{
	
	int i,j,ipen;
	float x1,y1,x2,y2,x3,y3,rad,ang;
	pinitf("GRID.PLT");
	/* put out a sequence of tics */
	newpen(1);
	for(i=0;i <= 100 ; i+=5){
		if((i%10)==0){
			plot((float)i*0.1, 0.0 ,3);
			plot((float)i*0.1, 8.0 ,2);
	}
	}
	for(i=0;i <= 100 ; i+=5){
		if((i%10)==0){
			plot( 0.0 ,(float)i*.1,3);
			plot(10.00,(float)i*.1,2);
	}
	}

	/* put in little squares */
	/* draw a box for orientation of X and Y axes */

	ipen = 1000;
	for(j=0;j<8;j+=2){
		for(i=0;i<10;i+=2){
			newpen(ipen+=5);
			shader((float)i,(float)j,(float)(i+1),
				(float)(j+1), 0, 0, 0.01,0.01);
		}
	}
	newpen(0);
	plot( 9.0,7.0,3);
	plot( 8.0,6.0,2);


	ipen = 1000;
	for(j=0;j<8;j+=2){
		y1 = j + 1.0;
		y2 = y1 + 1.0;
		y3 = y1;
		for(i=0;i<10;i+=2){
			x1 = i + 1.0;
			x2 = x1 + 0.5;
			x3 = x1 + 1.0;
			newpen(ipen+=5);
			shadet(x1,y1,x2,y2,x3,y3, 0, 0, 0.01,0.01);
		}
	}
	

	newpen(2);
	for(i=0; i < NARR; i++){
		ang = 6.2831853*(float)i/(float)(NARR-1);
		xarr[i] = 4.0 + 1.0*cos(ang);
		yarr[i] = 4.0 + 1.0*sin(ang);
	}
		shadep(NARR,xarr,yarr);
		


pend();

	
}
	
	
void box(float xl, float yl, float xh, float yh)
{
	plot(xl,yl,3);
	plot(xh,yl,2);
	plot(xh,yh,2);
	plot(xl,yh,2);
	plot(xl,yl,2);
	plot(xl,yl,3);
}
