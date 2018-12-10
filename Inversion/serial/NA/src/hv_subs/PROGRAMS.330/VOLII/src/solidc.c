#include <math.h>
#include "calplot.h"

#define DEGRAD 0.017453293

	

void gsolid(float x, float y, float h, int symb)
{
	/* solid fill with 
		symb	symbol
		0	square
		1	circle
		2 	triangle
		3	diamond
		4	octagon
		5	inverted triangle
	*/
	float dang, ang;
	int nside, i;
	float rad;
	float xp[17];
	float yp[17];
	symb %=  6;
	dang = 0.0 ;
	nside = 4 ;
	rad = 0.0 ;
	switch(symb){
		case 0:
			/* square */
			rad = 1.414 * h;
			nside = 4;
			dang = 45.0;
			break;
		case 1:
			/* circle */
			rad =  h;
			nside = 16;
			dang = 0.0;
			break;
		case 2:
			/* triangle */
			rad =  1.5*h;
			nside = 3;
			dang = 90.0;
			break;
		case 3:
			/* diamond */
			rad =  1.5*h;
			nside = 4;
			dang = 0.0;
			break;
		case 4:
			/* octagon */
			rad = 1.0823 * h;
			nside = 8;
			dang = 22.5;
			break;
		case 5:
			/* inverted triangle */
			rad =  1.5*h;
			nside = 3;
			dang = 270.0;
			break;
		default:
			break;
	}
	for(i=0; i< nside; i++){
		ang = dang + (i-1)*360.0/(float)nside;
		ang = ang * DEGRAD;
		xp[i] = x + rad * cos(ang);
		yp[i] = y + rad * sin(ang);
	}
	shadep(nside,xp,yp);
}
