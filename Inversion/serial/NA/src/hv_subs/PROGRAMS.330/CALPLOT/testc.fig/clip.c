
#include <stdio.h>
#include <math.h>
#include "calplot.h"

float xarr[] = { .900, 1.500, 1.700, 2.100, 1.700, 1.500, 1.400};
float yarr[] = {1.400,  .900, 1.350, 1.300, 1.900, 1.600, 1.900};

float nxarr[] = { 2.000, 2.500, 2.500, 2.000};
float nyarr[] = { 2.000, 2.000, 2.500, 2.500};


main()
{
	
	int narr, nnarr;	
	int i;
	int j;
	float xl, xh, yl, yh;

	ginitf("INTER","CLIP");
	narr = sizeof(xarr)/sizeof(float);
	nnarr = sizeof(nxarr)/sizeof(float);
	
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
	
