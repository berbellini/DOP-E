#include "calplot.h"
main()
{
	int ifont, i, ix, iy;
	float xx, yy, ang;
	float zero = 0.0;
	float height = 0.5;
	float width = 0.05;
	char str[2];
	int mone = -1;
	str[1] = '\0';
	pinitf("WIDTH.PLT");
	gunit("in");
	for(ifont=1;ifont<=4;ifont++){
		gfont(ifont);
		for(i=16;i<=111;i++){
			ix = (i-1)%16;
			iy = (i-1)/16;
			xx = 1.0 + ix*0.5;
			yy = 7.0 - (iy-1)*1.0;
			ang = 0.0;
			gwidth(width);
			str[0] = i;
			str[1] = '\0';
			symbol(xx,yy,height,str,ang,mone);
			gwidth(zero);
		}
		frame();
	}
	pend();
}
