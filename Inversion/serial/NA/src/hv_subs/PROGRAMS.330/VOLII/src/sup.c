#include "calplot.h"
#include "nmenu.h"
#include <stdio.h>
#include <string.h>

void mgwrtxt(float x0, float y0, char *str, int cmd, int color);
void clearregion(float xl, float yl, float xh, float yh);
void show_menu (float x0, float y0, struct menu *m, int size, int *nm);
int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);
void draw_button(float xl, float yl, char *str, float *xlw, 
	float *ylw, float *xup, float *yup, int butrev, int lstrmax) ;
void gwrstr(float x0, float y0, int color, char *fmt, char *str);
void gwrflt(float x0, float y0, int color, char *fmt, float fval);
void gwrint(float x0, float y0, int color, char *fmt, int   ival);

/* GLOBAL VARIABLE DECLARATIONS */
	/* display information */
extern int HasMouse; 
extern float XminDev, YminDev, 
	XmaxDev, YmaxDev, XminClip, 
	YminClip, XmaxClip, YmaxClip;
extern int Color;

/* DEFINES */
#define ON      1
#define OFF     0
#define XL      0.0
#define XH      10.00
#define YL      1.5
#define YH      8.0

static float xoff=0.0, yoff=0.0;

extern int  dc_top, dc_right;
void mgwrtxt(float x0, float y0, char *str, int cmd, int color)
{
	float ypix_per_in;
	float xpix_per_in;
	int usepixel;
	/* get pixel space for device */
	float DevWid, DevHgt;
	float height;
	DevWid = XmaxDev - XminDev;
	DevHgt = YmaxDev - YminDev;
	ypix_per_in = DevHgt/8.0;
	xpix_per_in = DevWid/10.0;
	/* if window is smaller than 640x480 VGA, do not use gwrtxt */
	if(DevWid < 620.0 || DevHgt < 440){
		usepixel = OFF;
		height = 0.125;
	} else {
		usepixel = ON;
		height = 0.105;
	}
	if(usepixel == ON){
		if(cmd > 0){
			if(cmd == 1)
				newpen(0);
			else
				newpen(1);
			shader(x0,y0-0.6*height,x0+strlen(str)*height,y0+2.0*height,0,0,0.01,0.01);
		}
		if(cmd == 2)
			newpen(0);
		else
			newpen(color);
		gwrtxt(x0,  y0 , str ,0);
		newpen(1);
	} else {
/*
		gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
*/
		if(cmd > 0){
			newpen(0);
			shader(x0,y0-0.6*height,x0+strlen(str)*height,y0+1.6*height,0,0,0.01,0.01);
			newpen(color);
		}
		newpen(color);
		symbol(x0,y0,height,str,0.0,strlen(str));
		newpen(1);
/*
		gclip("on" , pxl[1],pyl[1],pxh[1],pyh[1]);
fprintf(stderr,"(%f,%f) %f %s %d\n",x0,y0,height,str,strlen(str));
*/
	}
}


extern int black, kolor;

float border;
float title;
int Button_Color;
int Button_Color_Light;
int Button_Color_Dark;
int Button_Color_Fore;
int Button_Color_Back;
void show_menu (float x0, float y0, struct menu *m, int size, int *nm)
{
	int i;
	float xl, yl, xh, yh;
	int lstr, l;
	float but_wid, but_hgt;
	float xpos, ypos;
	*nm = size/sizeof(struct menu);
	lstr = 0;
	/* get the maximum string length in this menu category */
	for(i=0 ; i < *nm ; i++){
		l = strlen(m[i].str);
		if(l > lstr)lstr = l;
	}
	xpos = x0 - 1.125;
	ypos = y0;

	for(i=0 ; i < *nm ; i++){
		m[i].lstrmx = lstr;
		if(m[i].visible){
			if(m[i].xl < 0.0){
				xpos = xpos + 1.125;
				if((xpos +1.125) > XH){
					xpos = x0;
					ypos = ypos - 0.35;
				} 
				draw_button(xpos, ypos ,m[i].str,
					&xl,&yl,&xh,&yh,OFF,m[i].lstrmx);
				m[i].xl = xl;
				m[i].yl = yl;
				m[i].xh = xh;
				m[i].yh = yh;
				but_wid = xh - xl;
				but_hgt = yh - yl;
			} else {
				draw_button(m[i].xl, m[i].yl ,m[i].str,
					&xl,&yl,&xh,&yh,OFF,m[i].lstrmx);
				m[i].yl = yl;
				m[i].yh = yh;
			}
		}
	}
}
void draw_button(float xl, float yl, char *str, float *xlw, float *ylw,
float *xup, float *yup, int butrev, int lstrmax)
{
	/* draw the button. The button width is controlled by lstrmax
		but the string, smaller is centered in this window */
	int lstr;
	int tmp;
	float bdr = 0.05;
	float title ;
	float x[6], y[6];
	float ypix_per_in;
	float xpix_per_in;
	float wid, wids;
	float DevWid, DevHgt;
	float height;

	/* the button is defined by the bounding box
		(xl,yl) -> (xl+wid, yl+title)
		The actual region for the string is
		(xl+bdr,yl+bdr) -> (xl+wid-bdr,yl+title-bdr)
	*/
		
	/* get pixel space for device */
	DevWid = XmaxDev - XminDev;
	DevHgt = YmaxDev - YminDev;
	ypix_per_in = DevHgt/8.0;
	xpix_per_in = DevWid/10.0;
	
	lstr = strlen(str);
	lstrmax++;
	if(lstrmax <= 1)
		lstrmax = lstr;
	/* if window is smaller than 640x480 VGA, do not use gwrtxt */
	if(DevWid < 620.0 || DevHgt < 440){
		height = 0.125;
		wid = (float)lstrmax*height +bdr + bdr;
		wids= (float)lstr   *height;
		title = 1.3*height + bdr + bdr;
		yoff =  0.5*(title ) - height/2.0;
		xoff =  0.5*(wid - wids);
	} else {
		wid = (float)lstrmax*8.0/xpix_per_in +bdr + bdr;
		wids= (float)lstr   *8.0/xpix_per_in;
		title = 15.0/ypix_per_in + bdr + bdr;
		yoff =  0.5*(title ) - 6.5/ypix_per_in;
		xoff =  0.5*(wid - wids);
	}


	/* add a space */
	gclip("off", XL, YL, XH, YH);
	/* define colors for button */
	if(kolor == 1){
			/* gray */
		                Button_Color_Dark  = 1050;
                Button_Color       = 1040;
                Button_Color_Light = 1010;
                if(black){
                        Button_Color_Fore  =    1;
                } else {
                        Button_Color_Fore  =    0;
                }
        } else if(kolor == 2 || kolor == 0){
                Button_Color_Dark  = 1100;
                Button_Color       = 1085;
                Button_Color_Light = 1070;
                /* color */
                if(black){
                        Button_Color_Fore  =    0;
                } else {
                        Button_Color_Fore  =    1;
                }
        }


	if(butrev == ON){
		tmp = Button_Color_Dark;
		Button_Color_Dark = Button_Color_Light;
		Button_Color_Light = tmp;
	}

	newpen(Button_Color_Back);
	x[0] = xl ;
	y[0] = yl ;
	x[1] = xl + wid ;
	y[1] = yl ;
	x[2] = xl + wid ;
	y[2] = yl + title ;
	x[3] = xl ;
	y[3] = yl + title ;
	shadep(4,x,y);

	newpen(Button_Color_Dark);
	x[0] = xl;
	y[0] = yl;
	x[1] = xl + wid;
	y[1] = yl;
	x[2] = xl + wid;
	y[2] = yl + title ;
	x[3] = xl + wid - bdr;
	y[3] = yl + title - bdr ;
	x[4] = xl + wid - bdr;
	y[4] = yl + bdr ;
	x[5] = xl + bdr;
	y[5] = yl + bdr ;
	shadep(6,x,y);

	newpen(Button_Color_Light);
	x[0] = xl;
	y[0] = yl;
	x[1] = xl + bdr;
	y[1] = yl + bdr ;
	x[2] = xl + bdr;
	y[2] = yl + title - bdr ;
	x[3] = xl + wid - bdr;
	y[3] = yl + title - bdr ;
	x[4] = xl + wid ;
	y[4] = yl + title  ;
	x[5] = xl  ;
	y[5] = yl + title  ;
	shadep(6,x,y);

	newpen(Button_Color);
	x[0] = xl + bdr;
	y[0] = yl + bdr;
	x[1] = xl + wid - bdr;
	y[1] = yl + bdr;
	x[2] = xl + wid - bdr;
	y[2] = yl + title - bdr;
	x[3] = xl + bdr;
	y[3] = yl + title - bdr;
	shadep(4,x,y);
/*
	newpen(Button_Color_Fore);
*/
	/* this positioning should change according to the
			window size, with information from ginfo */
	mgwrtxt(xl+xoff,yl+yoff,str,0, Button_Color_Fore);
	*xlw = xl;
	*ylw = yl;
	*xup = xl +  bdr +bdr + wid;
	*yup = yl + title;
}

int inside(float xv, float yv, float xlb, 
float ylb, float xhb, float yhb)
{
	/* determine if the point is inside the bounding box */
	if(xv >= xlb && xv <= xhb && yv >= ylb && yv <= yhb)
		return 1;
		else
		return 0;
}

void clearregion(float xl, float yl, float xh, float yh)
{
	float x[4], y[4];
	x[0] = xl;
	y[0] = yl;
	x[1] = xh;
	y[1] = yl;
	x[2] = xh;
	y[2] = yh;
	x[3] = xl;
	y[3] = yh;
	newpen(0);
	shadep(4, x, y);
}

char outstr[80];
void gwrstr(float x0, float y0, int color, char *fmt, char *str){
	sprintf(outstr,fmt,str);
	mgwrtxt(x0,y0+yoff,outstr,1,color);
}
void gwrflt(float x0, float y0, int color, char *fmt, float fval){
	sprintf(outstr,fmt,fval);
	mgwrtxt(x0,y0+yoff,outstr,1,color);
}
void gwrint(float x0, float y0, int color, char *fmt, int   ival){
	sprintf(outstr,fmt,ival);
	mgwrtxt(x0,y0+yoff,outstr,1,color);
}
