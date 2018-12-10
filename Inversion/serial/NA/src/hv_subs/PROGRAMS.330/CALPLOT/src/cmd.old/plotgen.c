/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLOTGEN                                               c
c                                                                     c
c      COPYRIGHT (C)  1986, 1989 R. B. Herrmann                       c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
CHANGES
	01 01 2001 - changed definition of dc_scalex dc_xoff to
		be compatible with new CALPLOT which permits
		x,y to be negative -- this is invoked
		with the CALPLOT standard space(0,0,7620,7620)
		Otherwise the plot psace is mapped into the
		physical space

		added extern INT dc_CLipLeft_sv since this is
		required by the corrected EPS routines

	10 Jan 2001 - removed unused calls to dv_filltri and dv_fillbox
	02 APR 2004 - inplemented gend(mode) at higher level which
		required change in dv_closepl
	20 JUL 2004 - changes cast of strlen on line 310 to
		di_gsymb (userx,usery+10,100,0,(INT)strlen(s),s );
*/
/* general interface between plot driver and hardware specific driver 	*/
/* this provides service for the routines				*/
/*	linemod, arc, circle, linec, point, space, label, fillt		*/
/*	fills and fillr move cont					*/



/* set up pointers to graphics functions 
	that do the work */


/* hardware  versions */
/* Prototype for Device Dependent Functions */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#ifndef INT
#ifdef MSDOS
#define INT long
#else
#define INT int
#endif
#endif


void di_arc(INT Xi,INT Yi,INT X0,INT Y0,INT X1,INT Y1);
void di_circle(INT Xi,INT Yi,INT r);
void di_clip(INT cmd, INT X0,INT Y0,INT X1,INT Y1);
void di_closepl(int mode);
void di_cont(INT X0,INT Y0);
void di_cursor(INT curstyp);
void di_erase(INT mode);
void di_fillp(INT n, INT *x, INT *y);
void dtv_fillp(INT n, INT *x, INT *y); /* interface to dv_fillp that forces 
					  polygon clipping check */
void di_fillr(INT X0,INT Y0,INT X1,INT Y1,
	INT patx,INT paty,INT lenx,INT leny);
void di_fills(INT X0,INT Y0,INT ixy,INT istnd,INT iplmn);
void di_fillt(INT X0,INT Y0,INT X1,INT Y1,INT X2,INT Y2,
	INT patx,INT paty,INT lenx,INT leny);
void di_font(INT Xi);
void di_gottxt(char *s);
void di_gwid(INT wid);
void di_info(INT *HasMouse, INT *XminDev, INT *YminDev, 
	INT *XmaxDev, INT *YmaxDev, INT *XminClip, 
	INT *YminClip, INT *XmaxClip, INT *YmaxClip, INT *Color);
void di_label(char *s);
void di_linec(INT X0,INT Y0,INT X1,INT Y1);
void di_linemod(char *s);
void di_move(INT X0,INT Y0);
void di_openpl(int showmenu);
void di_pen(INT Xi);
void di_point(INT X0,INT Y0);
void di_space(INT X0,INT Y0,INT X1,INT Y1);
static int farea(float x0,float z0,float x1,float z1,float x2,float z2);
static int clipSuthHodge(INT n, INT *xin, INT *yin);

static void order(INT *x,INT *y);

/* prototype definitions */
extern void dv_rline(INT X0,INT Y0,INT X1,INT Y1);
extern void dv_closepl(int mode);
extern void dv_erase(INT mode);
extern void dv_font(INT Xi);
extern void di_gsymb(INT x, INT y, INT ht, INT ang, INT n, char *s);
extern void dv_cursor(INT curstyp);
extern void dv_openpl(int showmenu);
extern void dv_pen(INT xi);
extern void dv_fillp(INT n, INT *x, INT *y);
extern void dv_clip();
extern void dv_rliner(INT x0,INT z0,INT x1,INT z1);

/* external variables */
extern double dc_scalefac;
extern int dc_hasmouse;

extern INT dc_shadeoff;		/* if non-zero turn off shading		*/
extern INT dc_color;		/* is hardware capable of fast fill shading */
extern INT dc_right, dc_left, dc_top, dc_bottom;
extern INT dc_ClipRight, dc_ClipTop, dc_ClipBottom, dc_ClipLeft;
extern INT dc_ClipRight_sv, dc_ClipTop_sv, dc_ClipBottom_sv, dc_ClipLeft_sv;
extern INT dc_iseps;		/* EPS needs special scaling */

extern INT dc_ColorInfo;	/* current color map
	White background	0 = mono, 1 = Gray, 2 = Color 
	Black brakground	4 = inverted mono, 5 = inverted gray, 6 = inverted color
				*/

/* the following flags control seismic trace shading in which the
   trace is plotted parallel to either the x or y axes */
INT	dc_shdxy	= 0 ;	/* trace time parallel to x-axis (0) y (1) */
INT	dc_shdpm	= 0 ;	/* shade positive amplitudes (0) negative (1) */
INT	dc_shdse	= 0 ;	/* end trace shading (0) start (1)	      */
INT	dc_shdX0	= 0 ;	/* baseline definition for shading	      */
INT	dc_shdY0	= 0 ;	/* baseline coordinate for shading	      */
INT	dc_shdon	= 0 ;	/* flag set to specify whether in process of
				shading */
INT	dc_shdcur = 0 ;	/* current shading status 		      */
INT	dc_shdstart = 0;	/* has shading just started */
INT	dc_shdpts = 0;
void	dc_shd_flush(void);
void	dc_shd_fill(INT ix, INT iy);
#define NUMSHD 100
INT	dc_shd_x[NUMSHD], dc_shd_y[NUMSHD];




/* coordinate information					      */
INT	dc_oldx	= 0 ;	/* previous value of x coordinate */
INT	dc_oldy	= 0 ;	/* previous value of y coordinate */

double dc_xscale = 1.0;
double dc_yscale = 1.0;
double dc_xoff   = 1.0;
double dc_yoff   = 1.0;

static INT userx, usery, ouserx, ousery;

#define IX(n)	( dc_xscale *  dc_scalefac * ( n + dc_xoff ) + 0.49)
#define IY(n)   ( dc_yscale *  dc_scalefac * ( n + dc_yoff ) + 0.49)
#define LX(n)	( dc_xscale *  dc_scalefac * ( n  ) + 0.49)
#define LY(n)   ( dc_yscale *  dc_scalefac * ( n  ) + 0.49)
#define INVX(n) ( (n - 0.5)/(dc_xscale*dc_scalefac) - dc_xoff)
#define INVY(n) ( (n - 0.5)/(dc_yscale*dc_scalefac) - dc_yoff)

extern INT	dc_hardwarefill ;	/* flag to have all filling done by hardware */
				/* = 0 do shading here
				   = 1 dc_hardwarefill including seismic
				   = 2 dc_hardwarefill but not for seismic -
					which use low level triangles
					 */

extern double dc_Nxpi;
extern double dc_Nypi;
extern INT dc_rotate;
static void conts(INT xi,INT yi);
static void shade(INT u0,INT v0,INT u1,INT v1,INT dc_shdv0,INT xory,INT dc_shdpm);
static void shadex(INT X0,INT Y0,INT X1,INT Y1);
static void shadey(INT X0,INT Y0,INT X1,INT Y1);
static void hfill(INT i,INT Y1,INT Y2,INT patx,INT paty,INT lenx,INT leny) ;
static void vfill(INT j,INT X1,INT X2,INT patx,INT paty,INT lenx,INT leny) ;
static void trifil(INT *x,INT *y,INT patx,INT paty,INT lenx,INT leny);
static void ssort(INT *x,INT *y,INT n);

#define SGN(x)  ( (x) > 0 ? (1) : (-1))
#define	ABS(x)	( (x) > 0 ? (x) : (-x))
#define	MAX(a,b) ( (a) > (b) ? (a):(b) )
#define	MIN(a,b) ( (b) > (a) ? (a):(b) )




void di_arc(INT xi,INT yi,INT X0,INT Y0,INT X1,INT Y1)
{
}

void di_circle(INT xi,INT yi,INT r)
{
}

void di_move(INT xi,INT yi)
{
	ouserx = userx;
	ousery = usery;
	userx = xi;
	usery = yi;
	xi = IX(xi);
	yi = IY(yi);
	dc_oldx = xi;
	dc_oldy = yi;
}

void di_cont(INT xi,INT yi)
{
	INT xz,yz;
	ouserx = userx;
	ousery = usery;
	userx = xi;
	usery = yi;
	xi = IX(xi);
	yi = IY(yi);
	if(dc_shdse){
		/* dc_shdse != 0 sahded area seismic trace plot */
		if(dc_hardwarefill > 1){
			/* set flags for seismic shading */
			/* compute zero crossing */
			if(dc_shdxy == 0){	/* y(x)	*/
				/* 0 time axis is parallel to x-axis */
				if((SGN(yi-dc_shdY0) == SGN(dc_oldy-dc_shdY0) &&
					yi != dc_shdY0 && dc_oldy != dc_shdY0 )
						|| (yi == dc_oldy))
					conts(xi,yi);
				else {
					xz=dc_oldx+((double)(dc_shdY0-dc_oldy)*
						(double)(xi-dc_oldx)/
						(double)(yi-dc_oldy));
					conts(xz,dc_shdY0);
					conts(xi,yi);
				}
			} else {	/* x(y)	*/
				/* 0 time axis is parallel to y-axis */
				if((SGN(xi-dc_shdX0) == SGN(dc_oldx-dc_shdX0) && xi != dc_shdX0 && dc_oldx != dc_shdX0 ) || (xi == dc_oldx) ){
					conts(xi,yi);
				} else {
					yz=dc_oldy+((double)(dc_shdX0-dc_oldx)*
						(double)(yi-dc_oldy)/
						(double)(xi-dc_oldx));
					conts(dc_shdX0,yz);
					conts(xi,yi);
				}
			}
		} else {
			if(dc_shdxy)
				shadex(dc_oldx,dc_oldy,xi,yi);
			else
				shadey(dc_oldx,dc_oldy,xi,yi);
			dv_rline(dc_oldx,dc_oldy,xi,yi);
		}
	} else {
		dv_rline(dc_oldx,dc_oldy,xi,yi);
	}
	dc_oldx = xi;
	dc_oldy = yi;
}

static void conts(INT xi,INT yi)
{
	/* coming in here the  path is guaranteed one sided */
	if(dc_oldx == xi && dc_oldy == yi)return;
		if(dc_shdxy) {
			/* dc_shdxy != 0  time is y-axis for seismic */
			if((dc_shdpm==0||dc_shdpm==2) && xi >= dc_shdX0 && dc_oldx >= dc_shdX0)
				/* dc_shdpm == 0 positive amplitudes == 2 both positive and negative  */
				dc_shdon = 1;
			else if((dc_shdpm==1||dc_shdpm==2)&&xi<= dc_shdX0 && dc_oldx <= dc_shdX0)
				/* dc_shdpm == 0 positive amplitudes == 2 both positive and negative  */
				dc_shdon = 1;
			else
				dc_shdon = 0;
			if(dc_shdon && dc_shdstart){
/*
				dv_rline(dc_shdX0,dc_oldy,dc_shdX0,dc_oldy);
				dv_rline(dc_shdX0,dc_oldy,dc_oldx,dc_oldy);
*/
				dc_shd_fill( dc_shdX0,dc_oldy);
				dc_shd_fill( dc_oldx,dc_oldy);
				dc_shdstart = 0;
			}
		} else {
			/* dc_shdxy == 0 time is x-axis for seismic */
			if((dc_shdpm==0||dc_shdpm==2) && yi >= dc_shdY0 && dc_oldy >= dc_shdY0 )
				dc_shdon = 1;
			else if((dc_shdpm==1||dc_shdpm==2) && yi <= dc_shdY0 && dc_oldy <= dc_shdY0 )
				dc_shdon = 1;
			else
				dc_shdon = 0;
			if(dc_shdon && dc_shdstart){
/*
				dv_rline(dc_oldx,dc_shdY0,dc_oldx,dc_shdY0);
				dv_rline(dc_oldx,dc_shdY0,dc_oldx,dc_oldy);
*/
				dc_shd_fill( dc_oldx,dc_shdY0);
				dc_shd_fill( dc_oldx,dc_oldy);
				dc_shdstart = 0;
			}
		}
fprintf(stderr,"dc_shdstart %d dc_shdon %d dc_shdcur %d old [ %d,%d] new [%d,%d] \n",dc_shdstart,dc_shdon,dc_shdcur,dc_oldx,dc_oldy,xi,yi);
/*
	dv_rline(dc_oldx,dc_oldy,xi,yi);
*/
dc_shd_fill(dc_oldx,  dc_oldy);
dc_shd_fill(xi,  yi);
	dc_oldx = xi;
	dc_oldy = yi;
}

void di_linec(INT X0,INT Y0,INT X1,INT Y1)
{
	di_move(X0,Y0);
	di_cont(X1,Y1);
}

void di_linemod(char *s)
{
}

void di_point(INT xi,INT yi)
{
	ouserx = userx;
	ousery = usery;
	userx = xi;
	usery = yi;
	xi = IX(xi);
	yi = IY(yi);
	dv_rline(xi,yi,xi,yi);
}


void di_label(char* s)
{
	INT ox, oy;
	ox = userx;
	oy = usery;
/*
	di_gsymb (userx-30,usery-50,100,0,(INTbh)strlen(s),s );
*/
	di_gsymb (userx,usery+10,100,0,strlen(s),s );
	/* since dv_symvec calls move() cont() which change userx,
	modify userx, usery for proper recovery. Else
	two calls to label will not plot on same line 
	The 100 means a height of 0.10 inches. The character
	height is 100 and its actual width, neglecting
	intercharacter gaps is 0.6 of its height. Thus
	the choice of the units 50 and 30 for repositioning */
	userx = ox  + strlen(s)*100;
	usery = oy ;
}


void di_fillt(INT X0,INT Y0,INT X1,INT Y1,INT X2,INT Y2,
		INT patx,INT paty,INT lenx,INT leny)	
		/* fill  a triangular region */
{
	INT x[3],y[3];
	if(dc_shadeoff)
		return;
	/* temporary array set if solid fill */
	x[0] = X0; x[1] = X1 ; x[2] = X2;
	y[0] = Y0; y[1] = Y1 ; y[2] = Y2;
	if(dc_hardwarefill > 1 || (patx == 0 && paty == 0 && dc_color > 0)){
		di_fillp(3, x, y);
		return;
	}

	X0 = IX(X0);
	X1 = IX(X1);
	X2 = IX(X2);
	Y0 = IY(Y0);
	Y1 = IY(Y1);
	Y2 = IY(Y2);
	lenx = LX(lenx);
	leny = LY(leny);
	if(lenx <=0)lenx=1;
	if(leny <=0)leny=1;
	if( (X0==X1) && (X1==X2)){
		hfill(X0,Y0,Y1,patx,paty,lenx,leny);
		hfill(X0,Y0,Y2,patx,paty,lenx,leny);
	} else if( (Y0==Y1) && (Y1==Y2)){
		vfill(Y0,X0,X1,patx,paty,lenx,leny);
		vfill(Y0,X0,X2,patx,paty,lenx,leny);
	} else {
		x[0] = X0;
		x[1] = X1;
		x[2] = X2;
		y[0] = Y0;
		y[1] = Y1;
		y[2] = Y2;
		trifil(x,y,patx,paty,lenx,leny);
	}
}

static void trifil(INT *x,INT *y,INT patx,INT paty,INT lenx,INT leny)
{
	INT X0,X1,X2,Y0,Y1,Y2;
	INT i;
	INT j1,j2;
	float yinc1 = 0.0;
	float yinc2 = 0.0;
	float ycur1, ycur2;
	ssort(x,y,3);
	X0 = x[0];
	X1 = x[1];
	X2 = x[2];
	Y0 = y[0];
	Y1 = y[1];
	Y2 = y[2];
	if(X0 == X1 && X1 == X2){
		hfill(X0,Y0,Y1,patx,paty,lenx,leny);
		hfill(X0,Y0,Y2,patx,paty,lenx,leny);
		hfill(X0,Y1,Y2,patx,paty,lenx,leny);
		return;
	}
	if(X0 == X1){
		ycur1 = (float)Y0;
		ycur2 = (float)Y1;
		yinc2 = (float)(Y2-Y0)/(float)(X2-X0);
	} else {
		ycur1 = (float)Y0 + 0.49;
		ycur2 = (float)Y0 + 0.49;
		yinc2 = (float)(Y2-Y0)/(float)(X2-X0);
		yinc1 = (float)(Y1-Y0)/(float)(X1-X0);
	}
	j1 = (INT)ycur1 + 0.49;
	j2 = (INT)ycur2 + 0.49;
	hfill(X0,j1,j2,patx,paty,lenx,leny);
	for(i=X0+1;i<X1+1;i++){
		ycur1 = ycur1 + yinc1;
		ycur2 = ycur2 + yinc2;
		j1 = (INT)ycur1 + 0.49;
		j2 = (INT)ycur2 + 0.49;
		hfill(i,j1,j2,patx,paty,lenx,leny);
		}
	if(X1 != X2){
		yinc1 = (float)(Y2-Y1)/(float)(X2-X1);
	}
	else {
		yinc1 = 1.0;
	}
	for (i= X1+1;i < X2 +1 ; i++){
		if(X0 == X1){
			ycur1 = ycur1 +  yinc2;
			ycur2 = ycur2 +  yinc1;
		}
		else {
			ycur1 = ycur1 + yinc1;
			ycur2 = ycur2 + yinc2;
		}
		j1 = (INT)ycur1 + 0.49;
		j2 = (INT)ycur2 + 0.49;
		hfill(i,j1,j2,patx,paty,lenx,leny);
	}
}


/* fill horizontally */
static void hfill(INT i,INT Y1,INT Y2,INT patx,INT paty,INT lenx,INT leny) 
{
	INT yseg, ynext,y;
	if(Y2 < Y1){
		y = Y1;
		Y1 = Y2;
		Y2 = y;
	}
	if( patx == 0 || patx & 01<<( i/lenx)%6) {
		if(paty==0)
			dv_rliner(i,Y1,i,Y2);
		else {
			yseg = Y1/leny;
			ynext = (yseg+1)*leny;
			if(ynext < Y2){
				if(paty & 01<<(yseg%6))
					dv_rliner(i,Y1,i,ynext);
				for(y=ynext;y<= Y2-leny;y+=leny){
					if(paty & 01<<(y/leny)%6)
						dv_rliner(i,y,i,y+leny);
				}
			} else {
				y = Y1;
			}
			if(paty & 01<<((Y2)/leny)%6){
				dv_rliner(i,y,i,Y2);
			}
		}
	}
}

/* fill vertically */
static void vfill(INT j,INT X1,INT X2,INT patx,INT paty,INT lenx,INT leny) 
{
	if( paty == 0 || paty & 01<<( j/leny)%6) {
		dv_rliner(X1,j,X2,j);
	}
}

/* order coordinate pairs in increasing x */
static void ssort(INT *x,INT *y,INT n)		
{			/* shell sort The C Programming Language */
			/* Kernighan and Ritchie, 1978 Prentice Hall */
			/* p 58					*/
	INT gap, i , j, temp;
	for (gap = n/2; gap > 0 ; gap /= 2)
		for(i = gap; i < n; i++)
			for(j=i-gap;j>=0 && x[j]>x[j+gap];j-=gap){
				temp = x[j];
				x[j] = x[j+gap];
				x[j+gap] = temp;
				temp = y[j];
				y[j] = y[j+gap];
				y[j+gap] = temp;
			}
}

void di_fillr(INT X0,INT Y0,INT X1,INT Y1,INT patx,INT paty,INT lenx,INT leny)
{
	INT i;
	INT x[4],y[4];
	if(dc_shadeoff)
		return;
	/* temporary array set if solid fill */
	x[0] = X0; x[1] = X1 ; x[2] = X1; x[3] = X0;
	y[0] = Y0; y[1] = Y0 ; y[2] = Y1; y[3] = Y1;
	if(dc_hardwarefill > 1 || (patx == 0 && paty == 0 && dc_color > 0)){
		di_fillp(4, x, y);
		return;
	} 

	X0 = IX(X0);
	Y0 = IY(Y0);
	X1 = IX(X1);
	Y1 = IY(Y1);
	lenx = LX(lenx);
	leny = LY(leny);
	if(lenx <=0)lenx=1;
	if(leny <=0)leny=1;
	if(X1 < X0){
		i = X0;
		X0 = X1;
		X1 = i;
	}
	for(i=X0;i<X1+1;i++)
		hfill(i,Y0,Y1,patx,paty,lenx,leny);
}

/* seismic fill initiation
	X0, Y0 coordinates of start of trace
	ixy =	0 -> time axis = x-axis
			Y0 used to define +- trace region
		1 -> time axis = y-axis
			X0 used to define +- trace region
	istnd=	0 turn off trace shading
		1 turn on  trace shading
	iplmn=	0 fill in positive amplitudes
			y > Y0 for ixy = 0 or
			x > X0 for ixy = 1
		1 fill in negative amplitudes
			y < Y0 for ixy = 0 or
			x < X0 for ixy = 1
*/

void di_fills(INT X0,INT Y0,INT ixy,INT istnd,INT iplmn)
{
	if(dc_shadeoff)
		istnd = 0;
	dc_shdse = istnd;
	/* do proper termination when turning off before resetting all flags */
	if(!dc_shdse && dc_hardwarefill > 1){
		dc_shd_flush();
		if(dc_shdon){
			if(dc_shdxy){
dc_shd_fill(dc_oldx, dc_oldy);
dc_shd_fill(dc_shdX0, dc_oldy);
/*
				dv_rline(dc_oldx,dc_oldy,dc_shdX0,dc_oldy);
				*/
				dc_shdon = dc_shdse;
				/*
				dv_rline(dc_shdX0,dc_oldy,dc_shdX0,dc_oldy);
				*/
dc_shd_fill(dc_shdX0, dc_oldy);
dc_shd_fill(dc_shdX0, dc_oldy);
			} else {
				/*
				dv_rline(dc_oldx,dc_oldy,dc_oldx,dc_shdY0);
				*/
				dc_shd_fill(dc_oldx,dc_oldy);dc_shd_fill(dc_oldx,dc_shdY0);
				dc_shdon = dc_shdse;
				/*
				dv_rline(dc_oldx,dc_shdY0,dc_oldx,dc_shdY0);
				*/
				dc_shd_fill(dc_oldx,dc_shdY0);dc_shd_fill(dc_oldx,dc_shdY0);
				dc_shdon = dc_shdse;
			}
		}
			
	}
	dc_shdxy = ixy;
	dc_shdX0 = IX(X0);
	dc_shdY0 = IY(Y0);
	dc_shdpm = iplmn;
	dc_shdon = istnd;
	if(dc_shdon)
		dc_shdstart = 1;
}

INT dc_sdx, dc_sdy;	/* an increment in device x,y space used by dv_symvec to
		implement pseudobold */
		
void di_space(INT X0,INT Y0,INT X1,INT Y1)
{ /* force rectangular region defined by (X0,Y0) (X1,Y1) to map into */
/* specific plot space - in particular (0,0) (7620,7620) should		*/
/* force all coordinates to map 1000 input units into 1.00 inches	*/
/* for systems that are not so preset, we map (X0,X1) into		*/
/*	dc_left and some upper magic number to guarantee that 0,7620	*/
/*	maps into 0.0 to 7.62 inches					*/
	dc_xoff = 0.0;
	dc_yoff = 0.0;
/*
	 this is complicated 
	The original purpose is to map user space into the device space
	However this conflicts with CALPLOT Encapsulated PostScript
	and CALPLOT in general where the user space is the entire
	x-y plane and since CALPLOT forces the users to use 
	ABSOLUTE spatial coordinates instead of view window.

	So first test for (0,0,7620,7620)

	If this is seen then the saclefac = 1.0 and offset = 0
	
	Else do the mapping into device space.
*/
	if(dc_iseps){
		dc_xscale = (double)(dc_Nypi*7.620)/(double)(X1-X0);
		dc_yscale = (double)(dc_Nxpi*7.620)/(double)(Y1-Y0);
		dc_xoff = 0.0;
		dc_yoff = 0.0;
		dc_sdx = (INT)( 1.0/(dc_scalefac * dc_xscale) );
		dc_sdy = (INT)( 1.0/(dc_scalefac * dc_yscale) );
	} else {
		if(dc_rotate){
			dc_xscale = (double)(dc_Nypi*7.620)/(double)(X1-X0);
			if(dc_xscale != 0.0)
				dc_xoff = (double)dc_bottom/dc_xscale - (double)X0;
			dc_yscale = (double)(dc_Nxpi*7.620)/(double)(Y1-Y0);
			if(dc_yscale != 0.0)
				dc_yoff = (double)dc_left/dc_xscale - (double)Y0;
		} else {
			dc_xscale = (double)(dc_Nxpi*7.620)/(double)(X1-X0);
			if(dc_xscale != 0.0)
				dc_xoff = (double)dc_left/dc_xscale - (double)X0;
			dc_yscale = (double)(dc_Nypi*7.620)/(double)(Y1-Y0);
			if(dc_yscale != 0.0)
				dc_yoff = (double)dc_bottom/dc_xscale - (double)Y0;
		}
		if(dc_scalefac != 0.0){
			if(dc_xscale != 0.0)
				dc_sdx = (INT)( 1.0/(dc_scalefac * dc_xscale) );
			else
				dc_sdx = 0;
			if(dc_yscale != 0.0)
				dc_sdy = (INT)( 1.0/(dc_scalefac * dc_yscale) );
			else
				dc_sdy = 0;
		} else {
			dc_sdx = 0;
			dc_sdy = 0;
		}
	}
/*
fprintf(stderr,"scale x y (%f,%f) off x y (%f,%f)\n",dc_xscale,dc_yscale,dc_xoff,dc_yoff);
*/
}

static void shadey(INT X0,INT Y0,INT X1,INT Y1)
{
	shade(X0,Y0,X1,Y1,dc_shdY0,0,dc_shdpm);
}

static void shadex(INT X0,INT Y0,INT X1,INT Y1)
{
	shade(Y0,X0,Y1,X1,dc_shdX0,1,dc_shdpm);
}

static void shade(INT u0,INT v0,INT u1,INT v1,INT dc_shdv0,INT xory,INT dc_shdpm)
{
	/* general purpose routine for shading. Assume trace is aligned
	   along x-axis, e.g., y(x) and xory == 0 , and do a
	   call shade(X0,Y0,X1,Y1,dc_shdY0,0,dc_shdpm); This will give
	   shading when
		y(x) > dc_shdY0 if dc_shdpm == 0
		or
		y(x) < dc_shdY0 if dc_shdpm == 1
		or
		y(x) < dc_shdY0 and y(x) > dc_shdY0 if dc_shdpm == 2

	  To use this for x(y) to do the same thing, we set xory = 1, and do
	  a call  shadey(Y0,X0,Y1,X1,dc_shdX0,1,dc_shdpm); This will give
	  shading when
		x(y) > dc_shdX0 && dc_shdpm == 0
		or
		x(y) < dc_shdX0 && dc_shdpm == 1
		or
		x(y) < dc_shdX0 and x(y) > dc_shdX0  for dc_shdpm == 2
	  The reason for the xory flag is to get the correct arguments in
	  the external routines dv_rliner(4,) and dv_fillp()
	*/
	double slope;
	INT iu[4], iv[4];
	INT ulw, uup, vlw;
	INT u,v;
	if(u0 == u1) {
		slope = 0.0; /* a vertical step */
		ulw = u0;
		uup = u0;
		vlw = v1;
	}
	else {
		slope = (double)(v1-v0)/(double)(u1-u0);
		if( u1 > u0){
			ulw = u0;
			uup = u1;
			vlw = v0;
		} else {
			ulw = u1;
			uup = u0;
			vlw = v1;
		}
	}
	if(dc_hardwarefill == 0){
		for(u=ulw;u<uup;u++){
			v = vlw + (INT)(slope*(double)(u-ulw));
			if(v > dc_shdv0 && (dc_shdpm == 0 || dc_shdpm ==2 )){
				if(xory != 0){
					dv_rliner(dc_shdv0,u,v,u);
				} else {
					dv_rliner(u,dc_shdv0,u,v);
				}
			} else if(v < dc_shdv0 && (dc_shdpm == 1 || dc_shdpm ==2)){
				if(xory != 0){
					dv_rliner(dc_shdv0,u,v,u);
				} else {
					dv_rliner(u,dc_shdv0,u,v);
				}
			}
		}
	} else {
		if(v0 >= dc_shdv0 && v1 >= dc_shdv0 && (dc_shdpm == 0||dc_shdpm == 2)){
			iu[0] = u0; iv[0] = dc_shdv0;
			iu[1] = u0; iv[1] = v0;
			iu[2] = u1; iv[2] = v1;
			iu[3] = u1; iv[3] = dc_shdv0;
			if(xory !=0)
				dtv_fillp(4,iv,iu);
			else
				dtv_fillp(4,iu,iv);
		} else if(v0 <= dc_shdv0 && v1 <= dc_shdv0&&(dc_shdpm==1||dc_shdpm==2)){
			iu[0] = u0; iv[0] = dc_shdv0;
			iu[1] = u0; iv[1] = v0;
			iu[2] = u1; iv[2] = v1;
			iu[3] = u1; iv[3] = dc_shdv0;
			if(xory !=0) {
				dtv_fillp(4,iv,iu);
			} else {
				dtv_fillp(4,iu,iv);
			}
		} else {
			/* here define a triangle as trace crosses dc_shdv0 */
			if(u0 != u1 && v0 != v1){
				/* get zero u(dc_shdv0) */
				u = u0 + (INT)((double)(dc_shdv0-v0)/slope);
				if (((dc_shdpm==2|| dc_shdpm==0) && v1 > dc_shdv0)||(
					(dc_shdpm==2||dc_shdpm==1) && v1 < dc_shdv0)){
					iu[0] = u; iv[0] = dc_shdv0;
					iu[1] = u1; iv[1] = v1;
					iu[2] = u1; iv[2] = dc_shdv0;
					if(xory !=0)
						dtv_fillp(3,iv,iu);
					else
						dtv_fillp(3,iu,iv);
				} else if(((dc_shdpm==0||dc_shdpm==2)&& v0>dc_shdv0)||
					((dc_shdpm==1||dc_shdpm==2) && v0 < dc_shdv0)){
					iu[0] = u; iv[0] = dc_shdv0;
					iu[1] = u0; iv[1] = v0;
					iu[2] = u0; iv[2] = dc_shdv0;
					if(xory !=0)
						dtv_fillp(3,iv,iu);
					else
						dtv_fillp(3,iu,iv);
				}
			}	/* note we ignore the tiny vertical segment
				 when u0 = u1 */
		}
	}
}

extern INT dc_xlinewidth;
extern INT dc_ylinewidth;
extern INT dc_linewidth;
extern INT dc_newlinewidth;
void di_gwid(INT width)
{
	if(dc_rotate){
		dc_xlinewidth = (INT)((double)width*dc_Nypi*dc_scalefac/1000.0);
		dc_ylinewidth = (INT)((double)width*dc_Nxpi*dc_scalefac/1000.0);
	} else {
		dc_xlinewidth = (INT)((double)width*dc_Nxpi*dc_scalefac/1000.0);
		dc_ylinewidth = (INT)((double)width*dc_Nypi*dc_scalefac/1000.0);
	}
	dc_linewidth = (INT)((double)width*dc_scalefac);
	dc_newlinewidth = 1;
}


void di_clip(INT reset, INT xl, INT yl, INT xh, INT yh)
{
	/* CALPLOT does its own clipping of lines, points
		and polygons. Hardware clipping is only
		used to define the plot region 
		Thus there is no dv_clip() 04/18/97 */
	INT tmp;
	order(&xl,&xh);
	order(&yl,&yh);
	/* procede with clipping */
	xl = IX(xl);
	yl = IY(yl);
	xh = IX(xh);
	yh = IY(yh);
        if(dc_rotate){
                tmp = xl;
                xl = yl ;
                yl = dc_top - tmp;
                tmp = xh;
                xh = yh ;
                yh = dc_top - tmp;
        }
	order(&xl,&xh);
	order(&yl,&yh);
        if(reset == 0){
                dc_ClipLeft = dc_ClipLeft_sv;
                dc_ClipBottom = dc_ClipBottom_sv;
                dc_ClipRight = dc_ClipRight_sv;
                dc_ClipTop = dc_ClipTop_sv;
        } else {
                dc_ClipTop    = MIN(dc_ClipTop_sv,yh);
                dc_ClipRight  = MIN(dc_ClipRight_sv,xh);
                dc_ClipLeft   = MAX(xl,dc_ClipLeft_sv) ;
                dc_ClipBottom = MAX(yl,dc_ClipBottom_sv) ;
        }
	dv_clip();

}

void  di_info(INT *HasMouse, INT *XminDev, INT *YminDev, 
	INT *XmaxDev, INT *YmaxDev, INT *XminClip, 
	INT *YminClip, INT *XmaxClip, INT *YmaxClip,INT *Color)
{
	/* this is for the plot file, thus nothing is interactive */
	*HasMouse = dc_hasmouse;
	*XminDev = dc_left ;
	*YminDev = dc_bottom ;
	*XmaxDev = dc_right ;
	*YmaxDev = dc_top ;
	*XminClip = INVX(dc_ClipLeft) ;
	*YminClip = INVY(dc_ClipBottom) ;
	*XmaxClip = INVX(dc_ClipRight) ;
	*YmaxClip = INVY(dc_ClipTop) ;
	*Color = dc_ColorInfo;
}

void di_pen(INT xi)
{
	dv_pen(xi);
}

void di_closepl(int mode)
{
	dv_closepl(mode);
}
void di_cursor(INT cursor)
{
	dv_cursor(cursor);
}
void di_erase(INT mode)
{
	dv_erase(mode);
}
void di_font(INT font)
{
	dv_font(font);
}
void di_openpl(int showmenu)
{
	dv_openpl(showmenu);
}

/* first cut for poly clip */
static void draw_poly(INT n,INT *x,INT *y);

static INT *xarro, *yarro;
static int poly_clip(INT n, INT *x, INT *y);


void di_fillp(INT n, INT *x, INT *y)
{
	INT i;
	INT tmp;
	for(i=0;i < n ; i ++){
		x[i] = IX(x[i]);
		y[i] = IY(y[i]);
/* commenting this out fixes all
*/
		if(dc_rotate){
			tmp = x[i];
			x[i] = y[i] ;
			y[i] = dc_top - tmp;
		}
	}
	dtv_fillp(n, x, y);
}

/* force clipping into a polygon */
void dtv_fillp(INT n, INT *x, INT *y)
{
	INT nout;
	/* set aside temporary storage for output array */
	if((xarro = (INT *)calloc(n+n, sizeof(INT))) == NULL)
		return;
	if((yarro = (INT *)calloc(n+n, sizeof(INT))) == NULL){
		free(xarro);
		return;
	}
	/* do the polygon clipping */ 
	nout = poly_clip(n,x,y); 
	/* code change 2011/10/04 */
	if(nout > 2){
		if(dc_hardwarefill >= 1)
			dv_fillp(nout,xarro,yarro);
		else {
		/* OK here turn into triangles and do triangle fill */
			draw_poly(nout,xarro,yarro);
		}
	}
	free(xarro);
	free(yarro);
}
/* ------------------------------------------------------------ */

#ifndef INT
#ifdef MSDOS
#define INT long
#else
#define INT int
#endif
#endif

static int clipSuthHodge();

static int poly_clip(INT n, INT *x, INT *y)
{
	int iout;
	iout = clipSuthHodge(n, x, y);
	return (iout);
}


/*	
	clipSuthhodge
	Hearn, D., and M. Pauline Baker (1994).
	Computer Graphics, Prentice Hall, Englewood Cliffs 652pp

	pages 239-242
	converted from Pascal to C
*/

#define True 1
#define False 0
#define Left	0
#define Right	1
#define Bottom	2
#define Top	3

struct point {float x; float y ;};
typedef struct point wcPt2;

#define Boundary 4
static wcPt2 firstPoint[Boundary];
static int newBoundary[Boundary];
static wcPt2 s[Boundary];

static int inside( wcPt2 , int );
static wcPt2 intersect(wcPt2, wcPt2 , int );
static int lcross( wcPt2, wcPt2, int );
static void clipPoint (wcPt2, int);
static void closeClip();
static int succ(int);

static int succ(int b)
{
	return (b+1);

}

	

static int ptCnt = 0;


static int clipSuthHodge(INT n, INT *xin, INT *yin)
{
	int b;
	int i;
	wcPt2 ptsIn;
	ptCnt = 0;
	for(b = Left ; b <= Top; b++)newBoundary[b] = True;
	for (i=0; i< n; i++){
		ptsIn.x = xin[i];
		ptsIn.y = yin[i];
		clipPoint (ptsIn, Left);
	}
	closeClip();
	return ( ptCnt);
}

static int inside(wcPt2 p, int b)
{
	int Inside;
	Inside = True;
	switch (b)
	{
		case Left : if(p.x < dc_ClipLeft ) Inside = False ; break ;
		case Right : if(p.x > dc_ClipRight ) Inside = False ; break ;
		case Bottom : if(p.y < dc_ClipBottom ) Inside = False ; break ;
		case Top : if(p.y >dc_ClipTop ) Inside = False ; break ;
	}
	return (Inside);	
}
static int lcross( wcPt2 p1, wcPt2 p2, int b )
{
	int Cross;
	if ( inside(p1, b) == inside(p2, b) )
		Cross = False;
	else
		Cross = True;
        return (Cross);
}

static wcPt2 intersect(wcPt2 p1, wcPt2 p2, int b)
{
	float m;
	wcPt2 i;
	float n,d;
	n = (p1.y - p2.y);
	d = (p1.x - p2.x);
	if(n == 0.0)n = 0.001;
	if(d == 0.0)d = 0.001;
	m = n/d;
	switch (b)
	{
		case Left :
			i.x = dc_ClipLeft ;
			if(p1.y == p2.y)
				i.y = p1.y  ;
			else
				i.y = p2.y + (dc_ClipLeft  - p2.x )*m;
			break;
		case Right :
			i.x = dc_ClipRight ;
			if(p1.y == p2.y)
				i.y = p1.y ;
			else
				i.y = p2.y + (dc_ClipRight  - p2.x )*m;
			break;
		case Bottom:
			i.y = dc_ClipBottom ;
			if(p1.x == p2.x)
				i.x = p1.x ;
			else
				i.x = p2.x + (dc_ClipBottom  - p2.y )/m;
			break;
		case Top:
			i.y =dc_ClipTop ;
			if(p1.x == p2.x)
				i.x = p1.x ;
			else
				i.x = p2.x + (dc_ClipTop  - p2.y )/m;
			break;
	}
	return (i);
}

static void clipPoint(wcPt2 p, int b)
{
	wcPt2 i;
	if(newBoundary[b]){
		firstPoint[b] = p;
		newBoundary[b] = False;
	} else {
		if(lcross(p, s[b], b) == True){
			i = intersect(p, s[b], b);
			if(b < Top){
				clipPoint(i,succ(b));
			} else {
				xarro[ptCnt] = i.x;
				yarro[ptCnt] = i.y;
				ptCnt++ ;
			}
		}
	}
	s[b] = p;
	if(inside(p,b)){
		if(b < Top){
			clipPoint(p,succ(b));
		} else {
			xarro[ptCnt] = p.x;
			yarro[ptCnt] = p.y;
			ptCnt++ ;
		}
	}
}

static void closeClip()
{
	wcPt2 i;
	int b;
	for(b = Left ; b <= Top ; b++){
		if(lcross(s[b],firstPoint[b],b)){
			i = intersect(s[b],firstPoint[b], b);
			if(b < Top)
				clipPoint (i, succ(b));
			else {
				xarro[ptCnt] = i.x;
				yarro[ptCnt] = i.y;
				ptCnt++ ;
			}
		}
	}
}

/* ------------------------------------------------------------------ */

/*
 * poly_tri.c
 *
 * Program to take a polygon definition and convert it into triangles
 * that may then be rendered by the standard triangle rendering
 * algorithms.  This assumes all transformations have been performed
 * already and cuts them up into optimal triangles based on their
 * screen-space representation.
 *
 *	Copyright (c) 1988, Evans & Sutherland Computer Corporation
 *
 * Permission to use all or part of this program without fee is
 * granted provided that it is not used or distributed for direct
 * commercial gain, the above copyright notice appears, and
 * notice is given that use is by permission of Evans & Sutherland
 * Computer Corporation.
 *
 *	Written by Reid Judd and Scott R. Nelson at
 *	Evans & Sutherland Computer Corporation (January, 1988)
 *
 * To use this program, either write your own "draw_triangle" routine
 * that can draw triangles from the definitions below, or modify the
 * code to call your own triangle or polygon rendering code.  Call
 * "draw_poly" from your main program.
 */





#define V_MAX 1000	/* Maximum number of vertices allowed (arbitrary) */

#define BIG 1.0e30	/* A number bigger than we expect to find here */

#define COUNTER_CLOCKWISE 0
#define CLOCKWISE 1



/*
 * orientation
 *
 * Return either clockwise or counter_clockwise for the orientation
 * of the polygon.
 */

static int orientation( n, x, y )
    INT n;			/* Number of vertices */
    INT x[], y[];		/* The vertex list */
{
    float area;
    int i;

    /* Do the wrap-around first */
    area = x[n-1] * y[0] - x[0] * y[n-1];

    /* Compute the area (times 2) of the polygon */
    for (i = 0; i < n-1; i++)
	area += x[i] * y[i+1] - x[i+1] * y[i];

    if (area >= 0.0)
	return COUNTER_CLOCKWISE;
    else
	return CLOCKWISE;
} /* End of orientation */



/*
 * determinant
 *
 * Computes the determinant of the three points.
 * Returns whether the triangle is clockwise or counter-clockwise.
 */

static int determinant( p1, p2, p3, x, y )
    int p1, p2, p3;		/* The vertices to consider */
    INT x[], y[];			/* The vertex list */
{
    float x1, x2, x3, y1, y2, y3;
    float determ;

    x1 = x[p1];
    y1 = y[p1];
    x2 = x[p2];
    y2 = y[p2];
    x3 = x[p3];
    y3 = y[p3];

    determ = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    if (determ >= 0.0)
	return COUNTER_CLOCKWISE;
    else
	return CLOCKWISE;
} /* End of determinant */



/*
 * distance2
 *
 * Returns the square of the distance between the two points
 */

static float distance2( x1, y1, x2, y2 )
    float x1, y1, x2, y2;
{
    float xd, yd;		/* The distances in X and Y */
    float dist2;		/* The square of the actual distance */

    xd = x1 - x2;
    yd = y1 - y2;
    dist2 = xd * xd + yd * yd;

    return dist2;
} /* End of distance2 */



/*
 * no_interior
 *
 * Returns 1 if no other point in the vertex list is inside
 * the triangle specified by the three points.  Returns
 * 0 otherwise.
 */

static int no_interior( p1, p2, p3, x, y, vp, n, poly_or )
    int p1, p2, p3;		/* The vertices to consider */
    INT x[], y[];			/* The vertex list */
    int vp[];			/* The vertex pointers (which are left) */
    int n;			/* Number of vertices */
    int poly_or;		/* Polygon orientation */
{
    int i;			/* Iterative counter */
    int p;			/* The test point */

    for (i = 0; i < n; i++) {
	p = vp[i];		/* The point to test */
	if ((p == p1) || (p == p2) || (p == p3))
	    continue;		/* Don't bother checking against yourself */
	if (   (determinant( p2, p1, p, x, y ) == poly_or)
	    || (determinant( p1, p3, p, x, y ) == poly_or)
	    || (determinant( p3, p2, p, x, y ) == poly_or) ) {
	    continue;		/* This point is outside */
	} else {
	    return 0;		/* The point is inside */
	}
    }
    return 1;			/* No points inside this triangle */
} /* End of no_interior */



/*
 * draw_poly
 *
 * Call this procedure with a polygon, this divides it into triangles
 * and calls the triangle routine once for each triangle.
 *
 * Note that this does not work for polygons with holes or self
 * penetrations.
 */

static void draw_poly(INT n,INT *x,INT *y)
{
    int prev, cur, next;	/* Three points currently being considered */
    int vp[V_MAX];		/* Pointers to vertices still left */
    int count;			/* How many vertices left */
    int min_vert;		/* Vertex with minimum distance */
    int i;			/* Iterative counter */
    float dist;			/* Distance alcross this one */
    float min_dist;		/* Minimum distance so far */
    int poly_orientation;	/* Polygon orientation */
    float x0,x1,x2,z0,z1,z2;	/* Triangle structure */
	INT X0[3], Y0[3];
	INT patx, paty;
	INT xlen, ylen;
	patx = 0;
	paty = 0;
	xlen = 50.0;
	ylen = 50.0;

    if (n > V_MAX) {
	fprintf( stderr, "Error, more than %d vertices.\n", V_MAX);
	return;
    }

    poly_orientation = orientation( n, x, y );

    for (i = 0; i < n; i++)
	vp[i] = i;		/* Put vertices in order to begin */

/* Slice off clean triangles until nothing remains */

    count = n;
    while (count > 2) {
	min_dist = BIG;		/* A real big number */
	min_vert = 0;		/* Just in case we don't find one... */
	for (cur = 0; cur < count; cur++) {
		prev = ( cur==0       ? count-1 : cur-1);
		next = ( cur==count-1 ? 0       : cur+1);
	    /* Pick out shortest distance that forms a good triangle */
	    if (   (determinant( vp[prev], vp[cur], vp[next], x, y ) == poly_orientation)
		    /* Same orientation as polygon */
		&& no_interior( vp[prev], vp[cur], vp[next], x, y, vp, count, poly_orientation )
		    /* No points inside */
		&& ((dist = distance2( (float)x[vp[prev]], (float)y[vp[prev]],
				       (float)x[vp[next]], (float)y[vp[next]] )) < min_dist) )
		    /* Better than any so far */
	    {
		min_dist = dist;
		min_vert = cur;
	    }
	} /* End of for each vertex (cur) */

	/* The following error should "never happen". */
	if (min_dist == BIG)
	    fprintf( stderr, "Error: Didn't find a triangle.\n" );

	prev = ( min_vert==0       ? count-1 : min_vert-1);
	next = ( min_vert==count-1 ? 0       : min_vert+1);

/* Output this triangle */

	x0 = x[vp[prev]];
	z0 = y[vp[prev]];
	x1 = x[vp[min_vert]];
	z1 = y[vp[min_vert]];
	x2 = x[vp[next]];
	z2 = y[vp[next]];

	if(farea(x0,z0,x1,z1,x2,z2) == 1){
		X0[0] = x0;
		Y0[0] = z0;
		X0[1] = x1;
		Y0[1] = z1;
		X0[2] = x2;
		Y0[2] = z2;
		trifil(X0,Y0,patx, paty, xlen, ylen);
	}

/* Remove the triangle from the polygon */

	count -= 1;
	for (i = min_vert; i < count; i++)
	    vp[i] = vp[i+1];
    }


} /* End of draw_poly */

static int farea(float x0,float z0,float x1,float z1,float x2,float z2)
{
	if( ((x1-x0)*(z2-z0) - (x2-x0)*(z1-z0)) != 0.0)
		return (1);
	else
		return (0);
}


/* End of poly_tri.c */
static void order(INT *x,INT *y)
{
	INT tmp;
	if(*y < *x){
		tmp = *x;
		*x = *y;
		*y = tmp;
	}
}
void dc_shd_flush(void)
{
int i;
	if ( dc_shdpts > 1 ){
for(i=0;i<dc_shdpts;i++){
	fprintf(stderr,"flush %d [%d,%d]\n",i,dc_shd_x[i],dc_shd_y[i]);
}
		dtv_fillp(dc_shdpts, dc_shd_x, dc_shd_y);
	}
	dc_shdpts = 0;
}

void	dc_shd_fill(INT ix, INT iy)
{
	if(dc_shdon){
		if( dc_shdpts == NUMSHD -1 ){
			dc_shd_flush();
			dc_shd_x[0] = dc_shd_x[dc_shdpts];
			dc_shd_y[0] = dc_shd_y[dc_shdpts];
			dc_shdpts = 1 ;
		}
		dc_shd_x[dc_shdpts] = ix ;
		dc_shd_y[dc_shdpts] = iy ;
		dc_shdpts++;
	}
}
