
#include <stdio.h>
#include <math.h>
#include "calplot.h"

float xarr[] = { .900, 1.500, 1.700, 2.100, 1.700, 1.500, 1.400};
float yarr[] = {1.400,  .900, 1.350, 1.300, 1.900, 1.600, 1.900};

float nxarr[] = { 2.000, 2.500, 2.500, 2.000};
float nyarr[] = { 2.000, 2.000, 2.500, 2.500};
#define NPTS 161


main()
{
	float x[NPTS],y[NPTS],amp[NPTS];
	
	int narr, nnarr;	
	int i;
	int j;
	float xl, xh, yl, yh;

	ginitf("INTER","CLIP");
	narr = sizeof(xarr)/sizeof(float);
	nnarr = sizeof(nxarr)/sizeof(float);

	for(i=0;i< NPTS; i++){
		x[i] = i*0.01;
		amp[i] = fabs(0.25 * cos(3.1415927*x[i]/2.));
		y[i] = amp[i]*sin(3.1415927*x[i]*5./2.);
	}
	
	for(j=0; j<=1300;j+=100){
		gmesg("Click mouse button to change figure");
		gwidth(0.0);
		newpen(2);
		xl = 0.5 + 0.001*(float)(500 -700);
		yl = 0.5 + 0.001*(float)(500 -700);
		xh = 0.5 + 0.001*(float)(500 -700 +j+j);
		yh = 0.5 + 0.001*(float)(500 -700 +j+j);
		plot(xl,yl,3);
		plot(xh,yl,2);
		plot(xh,yh,2);
		plot(xl,yh,2);
		plot(xl,yl,2);

/*
*/

		newpen(3);
		plot(xarr[0],yarr[0],3);
		for (i=1;i<narr;i++)
			plot(xarr[i],yarr[i],2);
		plot(xarr[0],yarr[0],2);
		gclip("on", xl, yl, xh, yh);
		newpen(4);
		shadep( narr, xarr,  yarr);
/*
*/
	newpen(1);
	puttrc(1.0,1.0,x,y,NPTS);
		gwidth(0.0);
		newpen(3);
		plot(nxarr[0],nyarr[0],3);
		for (i=1;i<nnarr;i++)
			plot(nxarr[i],nyarr[i],2);
		plot(nxarr[0],nyarr[0],2);
		shadep(nnarr,nxarr, nyarr);
/*
*/
		gclip("off", xl, yl, xh, yh);
		frame();
	}
	pend();

}
	
puttrc(x0,y0,x,y,n)
float x0, y0, x[], y[];
int n;
{
	int i;
	float xx, yy;
	plot(x0,y0,-3);
	shdsei(0.0,0.0,0,1,0);
	for(i=0;i<n;i++){
		xx = x[i];
		yy = y[i];
		if(i == 0){
			plot(xx,yy,3);
		} else {
			plot(xx,yy,2);
		}
	}
	shdsei(0.0,0.0,0,0,0);
	plot(-x0,-y0,-3);
}
