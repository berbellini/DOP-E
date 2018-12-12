#include	<stdio.h>
#include "calplot.h"
#include <math.h>

float border;
float title;
int Button_Color;
int Button_Color_Light;
int Button_Color_Dark;
int Button_Color_Fore;
int kolor;
int black;
void box(float xl, float yl, float xh, float yh);
void draw_button(float xl, float yl, char *str, float *xlw, float *ylw,
		float *xup, float *yup);
static int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);
struct menu { float xl; float yl; float xh 
	; float yh ; char *str ; int action;  } ;
static struct menu m[] = { 
	{  -1.0, -1.0, -1.0, -1.0, "  Mode  \0" , 1},
	{  -1.0, -1.0, -1.0, -1.0, " Select \0" , 2},
	{  -1.0, -1.0, -1.0, -1.0, " Delete \0" , 3},
	{  -1.0, -1.0, -1.0, -1.0, "  Copy  \0" , 4},
	{  -1.0, -1.0, -1.0, -1.0, "  Exit  \0" , 5}
} ;

FILE *buterr;

 main()
{
	float xl,yl,xh,yh;
	int nm, i;
	int cmd;
	float xv, yv;
	char c[2];
	buterr= fopen("button.err","w");
	ginitf("INTER","BUTTON");
	border = 8.0/60.0;
	kolor = 1;
	black = 1;

	drawit();
	nm = sizeof(m)/sizeof(struct menu);
	for(i=0 ; i < nm ; i++){
		draw_button(8.5,3.0 -i*0.4,m[i].str,&xl,&yl,&xh,&yh);
		m[i].xl = xl;
		m[i].yl = yl;
		m[i].xh = xh;
		m[i].yh = yh;
	}
	/* loop for values */
	cmd = -1;
	for(; ;){
		gcursor("Arrow");
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nm ; i++){
			if(inside(xv,yv,m[i].xl,m[i].yl,m[i].xh,m[i].yh))
			{
				cmd = m[i].action;
				gmesg(m[i].str);
			}
		}
		if(cmd == 5)
			goto end;
		else if(cmd == 4){
			pinitf("BUTTON.PLT");
			drawit();
			pend();
			ginitf("INTER","BUTTON");
			gcursor("Arrow");
		}
	}
end:
	gmesg("Closing Session - Any Key to End");
	pend();
fclose(buterr);
}

static int inside(float xv, float yv, float xlb, 
	float ylb, float xhb, float yhb)
{
	if(xv >= xlb && xv <= xhb && yv >= ylb && yv <= yhb)
		return 1;
	else
		return 0;
}

void draw_button(float xl, float yl, char *str, float *xlw, float *ylw,
		float *xup, float *yup)
{
	int lstr;
	int nc;
	float bdr = 0.04;
	float title = 20.0/60.0;
	float x[6], y[6];
	lstr = strlen(str);
	/* add a space */
	nc = lstr*8;
		Button_Color_Dark  = 1100;
		Button_Color       = 1090;
		Button_Color_Light = 1070;
		/* color */
		if(black){
			Button_Color_Fore  =    0;
		} else {
			Button_Color_Fore  =    1;
		}
	newpen(Button_Color_Dark);
	x[0] = xl;
	y[0] = yl;
	x[1] = xl + lstr*0.125;
	y[1] = yl;
	x[2] = xl + lstr*0.125;
	y[2] = yl + title ;
	x[3] = xl + lstr*0.125 - bdr;
	y[3] = yl + title - bdr ;
	x[4] = xl + lstr*0.125 - bdr;
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
	x[3] = xl + lstr*0.125 - bdr;
	y[3] = yl + title - bdr ;
	x[4] = xl + lstr*0.125 ;
	y[4] = yl + title  ;
	x[5] = xl  ;
	y[5] = yl + title  ;
	shadep(6,x,y);

	newpen(Button_Color);
	x[0] = xl + bdr;
	y[0] = yl + bdr;
	x[1] = xl + lstr*0.125 - bdr;
	y[1] = yl + bdr;
	x[2] = xl + lstr*0.125 - bdr;
	y[2] = yl + title - bdr;
	x[3] = xl + bdr;
	y[3] = yl + title - bdr;
	shadep(4,x,y);
	newpen(Button_Color_Fore);
	gwrtxt(xl+0.0625,yl+0.015,str,0);
	*xlw = xl;
	*ylw = yl;
	*xup = xl + lstr*0.125;
	*yup = yl + title;
}

drawit()
{
	int ix, iy;
	char ch[2];
	int i, ipen;
	float xx, yy;
	box(3.0,4.0,6.0,6.0);
	gclip("on",3.0,4.0,6.0,6.0);
	ipen = 3;
	newpen(4);
	for(i=0;i<= 100 ; i++){
		xx = 3.0 + (i)*(6.0-3.0)/100.0;
		yy = 5.0 + 2.0*sin(10.*xx)+cos(10.*xx);
		plot(xx,yy,ipen);
		ipen = 2;
	}
	gclip("off",3.0,4.0,6.0,6.0);
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
