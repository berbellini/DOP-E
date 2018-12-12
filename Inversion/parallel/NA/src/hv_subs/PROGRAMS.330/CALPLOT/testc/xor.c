#include <math.h>
#include "calplot.h"

/* prototypes */
void dosector(float rad, float ang1, float ang2);
int getpen(int i);

/* produce a color wheel, that invokes xoring */

#define NUMCOLOR 36
#define NUMANG 1
float xarr[NUMANG+2], yarr[NUMANG+2];

main()
{
	pinitf("INTER");
	dowheel();
	pend();
}

dowheel()
{
	int i, j ;
	int irad;
	float dang, ang;
	float rad, drad, radmin;
	/* define the center */
	/* first put out a color wheel of maximum radius */
	dang = 6.2831853/NUMCOLOR;
	drad = 0.250;
	frame();
	for(j=0 ; j < 2 ; j++){
		if(j == 0)
			gmesg("Using Exclusive Or");
		else
			gmesg("Using Background Overwrite");
		plot(5.0,4.0,-3);
		newpen(3000);
		rad = 4.0;
		for(i=0 ; i < NUMCOLOR; i++){
			newpen(getpen(i));
			dosector(rad,i*dang, (i+1)*dang);
		}
		for(rad=4.0 ; rad > 0.5 ; rad = rad - drad){
			for(i=0 ; i < NUMCOLOR; i++){
				if(j == 0){
					newpen(3001);
					newpen(getpen(i));
				} else {
					newpen(0);
				}
				dosector(rad,i*dang, (i+1)*dang);
				newpen(3000);
				newpen(getpen(i));
				dosector(rad-drad,i*dang, (i+1)*dang);
			}
			radmin = rad;
		}
		for(rad=radmin ; rad <= 4.0 ; rad = rad + drad){
			for(i=0 ; i < NUMCOLOR; i++){
				if(j == 0){
					newpen(3001);
					newpen(getpen(i));
				} else {
					newpen(0);
				}
				dosector(rad,i*dang, (i+1)*dang);
				newpen(3000);
				newpen(getpen(i));
				dosector(rad-drad,i*dang, (i+1)*dang);
			}
		}
		plot(-5.0,-4.0,-3);
		if(j == 0)frame();
	}
}
void dosector(float rad, float ang1, float ang2)
{
	int narr;
	int i;
	float dang, ang;
	narr = NUMANG + 2;
	dang = ( ang2 - ang1)/NUMANG;
	for (i = 0 ; i <=NUMANG; i++){
		ang = ang1 + i*dang;
		xarr[i] = rad*cos(ang);
		yarr[i] = rad*sin(ang);
	}
	xarr[NUMANG+1] = 0.0;
	yarr[NUMANG+1] = 0.0;
	shadep(narr, xarr, yarr);
}

int getpen(int i)
{
	int dcolor;
	float p;
	/* first cycle from 0 to 7 then 1000 to 1100 */
	if(i <= 7){
		return (i);
	} else {
		/* map uniformly between 1000 and 1100 */
		if(NUMCOLOR > 1){
			p = (float)(i-8)/(float)(NUMCOLOR-8);
			dcolor = (int) ( (1-p)*1000.0 + (p*1100.0));
		} else
			dcolor = 100;
		return((int)dcolor);
	} 
	
}
