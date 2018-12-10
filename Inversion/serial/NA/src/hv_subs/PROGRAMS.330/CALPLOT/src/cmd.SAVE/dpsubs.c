#include	<stdio.h> 
#include <stdio.h> 
#define INT int

extern FILE *grfp;
static void order (INT *x, INT *y);

/***** everything below here is generic output */

extern int is_long;
static INT devClipXmn = 0;
static INT devClipYmn = 0;
static INT devClipXmx = 1000000000;
static INT devClipYmx = 1000000000;
static INT curClipXmn = 0;
static INT curClipYmn = 0;
static INT curClipXmx = 1000000000;
static INT curClipYmx = 1000000000;

int debug;

void putpi(INT a)
{
	if(debug)fprintf(stderr," PUTPI: %d\n",a);
	if(is_long == 0){
		putc((char)a&0377,grfp);
		putc((char)(a>>8)&0377,grfp);
	if(debug)fprintf(stderr," short: %8x %4x %4x \n",a,(char)a&0377,(char)(a>>8)&0377);
	} else {
		putc((char) a&0377,grfp);
		putc((char)(a>>8)&0377,grfp);
		putc((char)(a>>16)&0377,grfp);
		putc((char)(a>>24)&0377,grfp);
	if(debug)fprintf(stderr," long : %8x %4x %4x %4x %4x \n",a,(char)a&0377,(char)(a>>8)&0377,(char)(a>>16)&0377,(char)(a>>24)&0377);
	}
}

#define I_L(x) ( ((x)>32767 || (x)<-32767) ? 1 : 0 )

void dp_arc(INT Xi,INT Yi,INT X0,INT Y0,INT X1,INT Y1)
{
	is_long = I_L(Xi) || I_L(Yi) || I_L(X0)||I_L(Y0) || I_L(X1) || I_L(Y1) ;
	if(is_long == 0){
		putc('a',grfp);
	} else {
		putc('A',grfp);
	}
	putpi(Xi);
	putpi(Yi);
	putpi(X0);
	putpi(Y0);
	putpi(X1);
	putpi(Y1);
}

void dp_font(INT Xi)
{
	putc('B',grfp);
	is_long = 0;
	putpi(Xi);
}

void dp_circle(INT Xi,INT Yi,INT r)
{
	is_long = I_L(Xi) || I_L(Yi) || I_L(r);
	if(is_long == 0){
		putc('c',grfp);
	} else {
		putc('C',grfp);
	}
	putpi(Xi);
	putpi(Yi);
	putpi(r);
}

void  dp_info(INT *HasMouse, INT *XminDev, INT *YminDev, 
	INT *XmaxDev, INT *YmaxDev, INT *XminClip, 
	INT *YminClip, INT *XmaxClip, INT *YmaxClip, INT *Color)
{
	/* this is for the plot file, thus nothing is interactive */
	*HasMouse = 0;
	*XminDev = devClipXmn ;
	*YminDev = devClipYmn ;
	*XmaxDev = devClipXmx ;
	*YmaxDev = devClipYmx ;
	*XminClip = curClipXmn ;
	*YminClip = curClipYmn ;
	*XmaxClip = curClipXmx ;
	*YmaxClip = curClipYmx ;
	*Color    = 2;
	
}

void dp_cursor(INT curstyp)
{
	putc('g',grfp);
	is_long = 0;
	putpi(curstyp);
}

void dp_gsymb(INT X0,INT Y0,INT X1,INT Y1,INT n,char *s)
{
	int j;
	is_long = I_L(X0) || I_L(Y0) || I_L(X1) || I_L(Y1) || I_L(n) ;
	if(is_long == 0){
		putc('D',grfp);
	} else {
		putc('d',grfp);
	}
	putpi(X0);
	putpi(Y0);
	putpi(X1);
	putpi(Y1);
	putpi(n);
	if(n < 0)
		putc(s[0],grfp);
	else
		for(j=0;j< n;j++)
			putc(s[j],grfp);
	putc('\n',grfp);
}
void dp_erase(INT mode)
{
	if(mode == 0)
		putc('e',grfp);
	else
		putc('E',grfp);
}
void dp_linemod(char *s)
{
	putc('f',grfp);
	fputs(s,grfp);
}
void dp_pen(INT Xi)
{
	putc('h',grfp);
	is_long = 0;
	putpi(Xi);
}
void dp_clip(INT cmd, INT X0,INT Y0,INT X1,INT Y1)
{
	order(&X0,&X1);
	order(&Y0,&Y1);
	is_long = I_L(X0) || I_L(Y0) || I_L(X1) || I_L(Y1) ;
	if(is_long == 0){
		putc('k',grfp);
	} else {
		putc('K',grfp);
	}
	putpi(cmd);
	putpi(X0);
	putpi(Y0);
	putpi(X1);
	putpi(Y1);
	if(cmd == 1){
		curClipXmn = X0;
		curClipXmx = Y0;
		curClipYmn = X1;
		curClipYmx = Y1;
	} else {
		/* reset clipping */
		curClipXmn = devClipXmn;
		curClipXmx = devClipXmx;
		curClipYmn = devClipYmn;
		curClipYmx = devClipYmx;
	}
}
void dp_linec(INT X0,INT Y0,INT X1,INT Y1)
{
	is_long = I_L(X0) || I_L(Y0) || I_L(X1) || I_L(Y1) ;
	if(is_long == 0){
		putc('l',grfp);
	} else {
		putc('L',grfp);
	}
	putpi(X0);
	putpi(Y0);
	putpi(X1);
	putpi(Y1);
}
void dp_move(INT X0,INT Y0)
{
	is_long = I_L(X0) || I_L(Y0) ;
	if(is_long == 0){
		putc('m',grfp);
	} else {
		putc('M',grfp);
	}
	putpi(X0);
	putpi(Y0);
}
void dp_cont(INT X0,INT Y0)
{
	is_long = I_L(X0) || I_L(Y0) ;
	if(is_long == 0){
		putc('n',grfp);
	} else {
		putc('N',grfp);
	}
	putpi(X0);
	putpi(Y0);
}
void dp_gottxt(char *s)
{
	char *c;
	/*dp_move(xX0,yy0);*/
	putc('o',grfp);
	for(c=s;*c != '\0';c++)
		putc(*c,grfp);
	putc('\n',grfp);
}
void dp_point(INT X0,INT Y0)
{
	is_long = I_L(X0) || I_L(Y0) ;
	if(is_long == 0){
		putc('p',grfp);
	} else {
		putc('P',grfp);
	}
	putpi(X0);
	putpi(Y0);
}
void dp_space(INT X0,INT Y0,INT X1,INT Y1)
{
	order(&X0,&X1);
	order(&Y0,&Y1);
	is_long = I_L(X0) || I_L(Y0) || I_L(X1) || I_L(Y1) ;
	if(is_long == 0){
		putc('s',grfp);
	} else {
		putc('S',grfp);
	}
	putpi(X0);
	putpi(Y0);
	putpi(X1);
	putpi(Y1);
}
void dp_label(char *s)
{
	putc('t',grfp);
	fputs(s,grfp);
}
void dp_fillp(INT n,INT *x,INT *y)
{
	int i, ilong;
	INT xx, yy;
	/* now check the size of numbers */
	ilong = 0;
	for(i=0 ; i < n ; i++){
		ilong |= I_L(x[i]);
		ilong |= I_L(y[i]);
	}
	if(ilong == 0){
		putc('v',grfp);
	} else {
		putc('V',grfp);
	}
	is_long = ilong;
	putpi(n);
	for(i=0; i< n ; i++){
		putpi(x[i]);
		putpi(y[i]);
	}
}
void dp_gwid(INT wid)
{
	is_long = I_L(wid);
	if(is_long == 0){
		putc('W',grfp);
	} else {
		putc('w',grfp);
	}
	putpi(wid);
}
void dp_fillt(INT X0,INT Y0,INT X1,INT Y1,INT X2,INT Y2,
	INT patx,INT paty,INT lenx,INT leny)
{
	is_long = I_L(X0) || I_L(Y0) || I_L(X1) || I_L(Y1) || 
		I_L(X2) || I_L(Y2) || I_L(patx) || I_L(paty) 
		|| I_L(lenx) || I_L(leny) ;
	if(is_long == 0){
		putc('x',grfp);
	} else {
		putc('X',grfp);
	}
	putpi(X0);
	putpi(Y0);
	putpi(X1);
	putpi(Y1);
	putpi(X2);
	putpi(Y2);
	putpi(patx);
	putpi(paty);
	putpi(lenx);
	putpi(leny);
}
void dp_fillr(INT X0,INT Y0,INT X1,INT Y1,
	INT patx,INT paty,INT lenx,INT leny)
{
	is_long = I_L(X0) || I_L(Y0) || I_L(X1) || I_L(Y1) || 
		I_L(patx) || I_L(paty) 
		|| I_L(lenx) || I_L(leny) ;
	if(is_long == 0){
		putc('y',grfp);
	} else {
		putc('Y',grfp);
	}
	putpi(X0);
	putpi(Y0);
	putpi(X1);
	putpi(Y1);
	putpi(patx);
	putpi(paty);
	putpi(lenx);
	putpi(leny);
}
void dp_fills(INT X0,INT Y0,INT ixy,INT istnd,INT iplmn)
{
	is_long = I_L(X0) || I_L(Y0) || I_L(ixy) || I_L(istnd) || 
		I_L(iplmn) ;
	if(is_long == 0){
		putc('z',grfp);
	} else {
		putc('Z',grfp);
	}
	putpi(X0);
	putpi(Y0);
	putpi(ixy);
	putpi(istnd);
	putpi(iplmn);
}

	
static void order(INT *x,INT *y)
{
	INT tmp;
	if(*y < *x){
		tmp = *x;
		*x = *y;
		*y = tmp;
	}
}

/* get text string from terminal */
void dp_gintxt(int  cnt, char *s)
{
	s[0]='\0';
}

void dp_cross(INT *ix, INT *iy, char *c)
{
	*ix = 0;
	*iy = 0;
	c[0] = ' ';
	c[1] = '\0';
}

void dp_control(int type, int i1, int i2, int i3, int i4)
{
}
