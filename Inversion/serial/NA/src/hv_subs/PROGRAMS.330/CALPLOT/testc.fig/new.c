#include <stdio.h>
#include <math.h>
#include "calplot.h"

float xarr[] = {1.0, 2.0, 3.0, 2.8, 2.5, 2.0, 1.0};
float yarr[] = {1.0, 1.0, 2.0, 2.5, 2.0, 3.0, 2.0};

void box(float xl, float yl, float xh, float yh);

char strout[100];

main()
{
	
	int narr;	
	int i;
	int ix, iy;
	char txt[20];
	char txt1[10];
	char tetx[80];

	char ch[2];
float xv, yv;
	int HasMouse; 
	float XminDev, YminDev, 
	XmaxDev, YmaxDev, XminClip, 
	YminClip, XmaxClip, YmaxClip;
	int Color;

	narr = sizeof(xarr)/sizeof(float);
	ginitf("NEW.PLT","NEW");
	/* put out a sequence of tics */
	newpen(1);
	for(i=0;i <= 100 ; i+=5){
		if(i%10){
			plot((float)i*0.1,0.0 ,3);
			plot((float)i*0.1,0.15,2);
		} else {
			plot((float)i*0.1,0.0 ,3);
			plot((float)i*0.1,0.3 ,2);
			number((float)i*0.1-0.07,0.35,0.14,(float)(0.1*i),0.0,-1);
		}
	}
	for(i=0;i <= 100 ; i+=5){
		if(i%10){
			plot(0.0 ,(float)i*.1,3);
			plot(0.15,(float)i*.1,2);
		} else {
			plot(0.0 ,(float)i*.1,3);
			plot(0.3 ,(float)i*.1,2);
			number(0.35,(float)i*0.1-0.07,0.14,(float)(0.1*i),0.0,-1);
		}
	}
	/* draw a box for orientation of X and Y axes */
	newpen(1);
	plot(0.0,0.0,3); 
	plot(1.0,0.0,2); 
	plot(1.0,1.0,2); 
	plot(0.0,1.0,2);
	plot(0.0,0.0,2);
	symbol(0.5,0.10,0.10,"X",0.0,1);
	symbol(0.10,0.5,0.10,"Y",0.0,1);
	newpen(2);
	plot(0.0,0.0,3); 
	plot(2.0,0.0,2); 
	plot(2.0,1.5,2); 
	plot(0.0,1.5,2);
	
	newpen(1000);


	newpen(4);
	/* plot a border around the polygon - note this is
		not a polygon clip */
	gwidth(0.05);
	for(i=0;i< narr ; i++){
		if(i == 0 )
			plot(xarr[i],yarr[i],3);
		else
			plot(xarr[i],yarr[i],2);
	}
	plot(xarr[0],yarr[0],2);

	/* draw clip box */
	plot (2.0,2.0,3);
	plotd(3.0,2.0,9,0.05);
	plotd(3.0,4.0,9,0.05);
	plotd(2.0,4.0,9,0.05);
	plotd(2.0,2.0,9,0.05);
	/* turn on clip box */
	gclip("on",2.0,2.0,3.0,4.0);

	/* now do a polygon fill in green */
	newpen(1050);
	shadep(narr,xarr,yarr);
	newpen(4);

	gwidth(0.001);
	/* crudely demostrate  clip region */
	for(i=250;i<10000;i+=250){
		yv = (float)i/1000.0;
		plot(0.0,yv,3); plot(10.0,yv,2);
	}

	/* shade a triangle - this should not be seen */
	shadet(1.0,1.0,2.0,5.0,1.0,0.5, 0,0, 0.05,0.05);

	
	/* draw a straight line to see how it is clipped */
	newpen(5);
	gwidth(0.10);
	plot(0.5,1.5,3);plot(4.0,4.0,2);
	gwidth(0.001);
	gclip("on",0.0,0.0,5.0,5.0);
	plot(-10.0,-10.0,3);
	plot(5.0,10.0,2);
	plot(40.0,40.0,3);
	gclip("off",-1.0,-1.0,-1.0,-1.0);
	box(5.0,4.0,8.0,6.0);
	gclip("on",5.0,4.0,8.0,6.0);
	gmesg("Crosshair Cursor in Box Arrow Outside");
	gcursor("Cross");
	cross(&ix,&iy,ch);
sprintf(strout,"ix: %d iy: %d c:%o",ix,iy,ch[0]);
	gwrtxt(5.5,5.0,strout,1);

	ginfo(&HasMouse, &XminDev, &YminDev, 
	&XmaxDev, &YmaxDev, &XminClip, 
	&YminClip, &XmaxClip, &YmaxClip,&Color);
	sprintf(strout,"HasMouse = %d",HasMouse);
	gwrtxt(0.5,7.5,strout,1);
	sprintf(strout,"Device Limits: (%f,%f) - (%f,%f)",XminDev,YminDev,
		XmaxDev,YmaxDev);
	gwrtxt(0.5,7.0,strout,1);
	sprintf(strout,"Clip Limits  : (%f,%f) - (%f,%f)",XminClip,YminClip,
		XmaxClip,YmaxClip);
	gwrtxt(0.5,6.5,strout,1);
	sprintf(strout,"Color = %d",Color);
	gwrtxt(0.5,6.0,strout,1);
	/* turn off the clip region to see if all is recovered */

	gclip("off",5.0,4.0,8.0,6.0);
	box(5.0,6.9,8.0,7.9);
	gmesg("Arrow Cursor");
	gcursor("Arrow");
		cross(&ix,&iy,ch);
		sprintf(strout,"ix: %d iy: %d c:%o",ix,iy,ch[0]);
		newpen(2);
		gwrtxt(5.5,7.5,strout,1);
	gmesg("Crosshair Cursor");
	gcursor("Cross");
		cross(&ix,&iy,ch);
		sprintf(strout,"ix: %d iy: %d c:%o",ix,iy,ch[0]);
		newpen(3);
		gwrtxt(5.5,7.0,strout,1);
	gmesg("Plus Cursor");
	gcursor("Plus");
		cross(&ix,&iy,ch);
		sprintf(strout,"ix: %d iy: %d c:%o",ix,iy,ch[0]);
		newpen(4);
		gwrtxt(5.5,6.5,strout,1);
	gmesg("XOR Arrow Cursor");
	gcursor("Xorarrow");
		cross(&ix,&iy,ch);
		sprintf(strout,"ix: %d iy: %d c:%o",ix,iy,ch[0]);
		newpen(4);
		gwrtxt(5.5,6.0,strout,1);
	gwrtxt(3.0,2.0,"HELLO WORLD   ",0);
	gwrtxt(3.0,1.5,"HELLO WORLD   ",1);
	gwrtxt(3.0,1.0,"HELLO WORLD   ",2);
	gwrtxt(3.0,0.5,"ENTER STRING  ",2);
	gwrtxt(3.0,0.5,"              ",2);
	plot(3.0,0.5,3);
	gcursor("Plus");
	grdtxt(txt1,10);
	gwrtxt(3.0,0.2,txt1,2);

	frame();
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
