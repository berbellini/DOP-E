#include <math.h>
#include <stdlib.h>
#include "calplot.h"

#ifndef MSDOS
float My_Random(int ival)
{
/* rand() returns a value between 0 and RAND_MAX from stdlib */
	return ((float)ival*(float)rand()/(float)RAND_MAX);
}
#else
#define My_Random random
#endif

/* prototypes */
void dowheel();
void dostar(float x0,float y0,float rad,float dang );

/* produce a color wheel, that invokes xoring */

#define NUMCOLOR 24
#define NPOLY 10
float xarr[NPOLY], yarr[NPOLY];

main()
{
	pinitf("STARS.PLT");
	dowheel();
	pend();
}

void dowheel()
{
	int i;
	int color;
	float dang, ang;
	float rad;
	float x0, y0;
	/* define the center */
	plot(5.0,4.0,-3);
	gclip("on",-4.0,-2.0,4.0,2.0);
	x0 = 0.0; y0 = 0.0; dang = 0.0; rad = 5.0;
	for(i=0 ; i < 10000; i++){
		x0 = 0.0 + 8.0*(1.0-0.02*My_Random(100));
		y0 = 0.0 + 8.0*(1.0-0.02*My_Random(100));
		rad = 0.5 * 0.01*My_Random(100);
		dang = My_Random(30);
		newpen(1000 + My_Random(100));
		
		dostar(x0,y0,rad,dang );
	}
}
void dostar(float x0,float y0,float rad,float dang )
{

	int i;
	float degrad = 3.1415927/180.0;
	float  ang;
	float radius;
	for (i = 0 ; i <NPOLY; i++){
		ang = dang + (i-1)*36;
		ang = ang * degrad;
		if(i%2)	
			radius = rad;
		else
			radius = 0.3819*rad;
		xarr[i] = x0 + radius*cos(ang);
		yarr[i] = y0 + radius*sin(ang);
	}
	shadep(NPOLY, xarr, yarr);
}

