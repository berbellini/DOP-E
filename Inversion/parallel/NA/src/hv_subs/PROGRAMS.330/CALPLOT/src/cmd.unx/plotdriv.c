/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLOTDRIVER                                            c
c                                                                     c
c      COPYRIGHT (C)  1997 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
*/
/* BASIC PLOT FILTER DRIVER FOR VECTOR DEVICES
	( graphics terminals, x-y plotters )
 */
#include <stdio.h>
#include <signal.h>
#include <stdlib.h>
#include <stddef.h>

	static int sigs[]= { SIGHUP, SIGINT , SIGQUIT, SIGBUS, SIGKILL,
		SIGABRT, SIGTERM, SIGSEGV, SIGFPE, SIGILL };

/* set up pointers to graphics functions 
	that do the work */

#ifndef INT
#ifdef MSDOS
#define INT long
#else
#define INT int
#endif
#endif
void (*do_arc)(INT Xi,INT Yi,INT X0,INT Y0,INT X1,INT Y1)=0;
void (*do_circle)(INT Xi,INT Yi,INT r)=0;
void (*do_clip)(INT cmd, INT X0,INT Y0,INT X1,INT Y1)=0;
void (*do_cont)(INT X0,INT Y0)=0;
void (*do_cursor)(INT curstyp)=0;
void (*do_erase)(INT mode)=0;
void (*do_fillp)(INT n,INT *x,INT *y)=0;
void (*do_fillr)(INT X0,INT Y0,INT X1,INT Y1,
	INT patx,INT paty,INT lenx,INT leny)=0;
void (*do_fills)(INT X0,INT Y0,INT ixy,INT istnd,INT iplmn)=0;
void (*do_fillt)(INT X0,INT Y0,INT X1,INT Y1,INT X2,INT Y2,
	INT patx,INT paty,INT lenx,INT leny)=0;
void (*do_font)(INT Xi)=0;
void (*do_gottxt)(char *s)=0;
void (*do_gsymb)(INT X0,INT Y0,INT X1,INT Y1,INT n,char *s)=0;
void (*do_gwid)(INT wid)=0;
void  (*do_info)(INT *HasMouse, INT *XminDev, INT *YminDev, 
	INT *XmaxDev, INT *YmaxDev, INT *XminClip, 
	INT *YminClip, INT *XmaxClip, INT *YmaxClip, INT *Color)=0;
void (*do_label)(char *s)=0;
void (*do_linec)(INT X0,INT Y0,INT X1,INT Y1)=0;
void (*do_linemod)(char *s)=0;
void (*do_move)(INT X0,INT Y0)=0;
void (*do_pen)(INT Xi)=0;
void (*do_point)(INT X0,INT Y0)=0;
void (*do_space)(INT X0,INT Y0,INT X1,INT Y1)=0;

extern	void di_gcmdln( int, char ** );

#include "disubs.h"

/* internal function definitions */
void	openpl( int showmenu );
void	clospl( int mode );

void get_arc(INT *Xi,INT *Yi,INT *X0,INT *Y0,INT *X1,INT *Y1);
void get_circle(INT *Xi,INT *Yi,INT *r);
void get_clip(INT *reset, INT *X0,INT *Y0,INT *X1,INT *Y1);
void get_line(INT *X0,INT *Y0,INT *X1,INT *Y1);
void get_move(INT *Xi,INT *Yi);
void get_cont(INT *Xi,INT *Yi);
void get_point(INT *Xi,INT *Yi);
void get_space(INT *X0,INT *Y0,INT *X1,INT *Y1);
void get_poly(INT *n,INT **x,INT **y);
void get_wid(INT *Xi);
void get_tri(INT *X0,INT *Y0,INT *X1,INT *Y1,
		INT *X2,INT *Y2,INT *patx,INT *paty,INT *lenx,INT *leny);
void get_rect(INT *X0,INT *Y0,INT *X1,INT *Y1,
		INT *patx,INT *paty,INT *lenx,INT *leny);
void get_sei(INT *X0,INT *Y0,INT *ixy,INT *istnd,INT *iplmn);
void get_sym(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *n,char *s);
static void order (INT *x, INT *y);
static INT getpi(FILE *);
static void getsf(char *c, FILE *f);
static void ngetsf(char *c, FILE *f, INT n);

static void fplt(FILE *);

static int	gphopen = 0 ;/* truth value for openpl used to delay open until
			   there is actual input useful in pipelining
			   where program does not start plot stream
			   until after interactive mode and we do not
			   want screen messed up beforehand */
int   argcc;
char **argvv;

int dc_sleeptime = 2;

/* external global variables */
extern INT     NumX     ;
extern INT     NumY     ;

#define NAMLEN  80
char dv_icon_name[NAMLEN];

extern void onintr(int arg);



int main(int argc, char **argv)
{
	int i;
	argcc = argc; /* save command arguments for other routines */
	argvv = argv;
	dv_icon_name[0] = '\0';
	for ( i=0; i< sizeof(sigs)/sizeof(int) ; ++i)
		if(signal( sigs[i],SIG_IGN) != SIG_IGN )
			signal ( sigs[i], onintr);

	di_gcmdln(argc,argv);
	fplt(stdin);
	return(0);
}

int is_long;
FILE *fin;

static void fplt(dfin)
FILE *dfin;
{
	int c;
	char s[256];
	INT Xi,Yi,X0,Y0,X1,Y1,X2,Y2,patx,paty,ixy,istnd,iplmn,r,n,i;
	INT reset;
	INT lenx, leny;
	INT *x, *y;

	fin = dfin;
#ifdef	MSDOS
	setmode(fileno(fin), O_BINARY);
#endif
	while((c=getc(fin)) != EOF){
		if(gphopen == 0){
			gphopen = 1;
			openpl(1);
		}
		switch(c){
		case 'a':			/* draw arc - first */
			is_long = 0;
			get_arc(&Xi,&Yi,&X0,&Y0,&X1,&Y1);
			(*do_arc)(Xi+NumX,Yi+NumY,X0,Y0,X1,Y1);
			break;
		case 'A':
			is_long = 1;
			get_arc(&Xi,&Yi,&X0,&Y0,&X1,&Y1);
			(*do_arc)(Xi+NumX,Yi+NumY,X0,Y0,X1,Y1);
			break;
		case 'B':			/* define font integer */
			is_long = 0;
			Xi = getpi(fin);
			(*do_font)(Xi);
			break;
		case 'c':			/* circle	- first */
			is_long = 0;
			get_circle(&Xi,&Yi,&r);
			(*do_circle)(Xi+NumX,Yi+NumY,r);
			break;
		case 'C':			/* 4 byte CIRCLE  */
			is_long = 1;
			get_circle(&Xi,&Yi,&r);
			(*do_circle)(Xi+NumX,Yi+NumY,r);
			break;
		case 'D':			/* implement symbol */
			is_long = 0;
			get_sym(&X0,&Y0,&X1,&Y1,&n,s);
			(*do_gsymb)(X0+NumX,Y0+NumY,X1,Y1,n,s);
			break;
		case 'd':			/* implement symbol */
			is_long = 1;
			get_sym(&X0,&Y0,&X1,&Y1,&n,s);
			(*do_gsymb)(X0+NumX,Y0+NumY,X1,Y1,n,s);
			break;
		case 'e':			/* erase plot frame */
			(*do_erase)(0);
			break;
		case 'E':			/* erase plot frame */
			(*do_erase)(1);
			break;
		case 'f':			/* read linemod string */
			getsf(s,fin);
			(*do_linemod)(s);
			break;
		case 'g':			/* set cursor */
			is_long = 0;
			Xi = getpi(fin);
			(*do_cursor)(Xi);
			break;
		case 'h':			/* define color integer */
			is_long = 0;
			Xi = getpi(fin);
			(*do_pen)(Xi);
			break;
		case 'k':			/* draw line from first */
			is_long = 0;
			get_clip(&reset,&X0,&Y0,&X1,&Y1);
			order(&X0,&X1);
			order(&Y0,&Y1);
			(*do_clip)(reset,X0+NumX,Y0+NumY,X1+NumX,Y1+NumY);
			break;
		case 'K':			/* draw line from first */
			is_long = 1;
			get_clip(&reset,&X0,&Y0,&X1,&Y1);
			order(&X0,&X1);
			order(&Y0,&Y1);
			(*do_clip)(reset,X0+NumX,Y0+NumY,X1+NumX,Y1+NumY);
			break;
		case 'l':			/* draw line from first */
			is_long = 0;
			get_line(&X0,&Y0,&X1,&Y1);
			(*do_linec)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY);
			break;
		case 'L':			/* 4 byte LINE */
			is_long = 1;
			get_line(&X0,&Y0,&X1,&Y1);
			(*do_linec)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY);
			break;
		case 'm':			/* move - two integers */
			is_long = 0;
			get_move(&Xi,&Yi);
			(*do_move)(Xi+NumX,Yi+NumY);	/* pen up	*/
			break;
		case 'M':			/* 4 byte MOVE */
			is_long = 1;
			get_move(&Xi,&Yi);
			(*do_move)(Xi+NumX,Yi+NumY);		/* pen up	*/
			break;
		case 'n':			/* continue - two integers */
			is_long = 0;
			get_cont(&Xi,&Yi);
			(*do_cont)(Xi+NumX,Yi+NumY);	/* new becomes current */
			break;
		case 'N':			/* 4 byte CONT */
			is_long = 1;
			get_cont(&Xi,&Yi);
			(*do_cont)(Xi+NumX,Yi+NumY);	/* new becomes current */
			break;
		case 'o':			/* output text to terminal */
			getsf(s,fin); 		/* at current position */
			(*do_gottxt)(s);
			break;
		case 'p':			/* draw a point at the two */
			is_long = 0;
			get_point(&Xi,&Yi);
			(*do_point)(Xi+NumX,Yi+NumY);
			break;
		case 'P':			/* 4 byte POINT */
			is_long = 1;
			get_point(&Xi,&Yi);
			(*do_point)(Xi+NumX,Yi+NumY);
			break;
		case 's':			/* space - define plot space */
			is_long = 0;
			get_space(&X0,&Y0,&X1,&Y1);
			order(&X0,&X1);
			order(&Y0,&Y1);
			(*do_space)(X0,Y0,X1,Y1);
			break;
		case 'S':			/* 4 byte SPACE */
			is_long = 1;
			get_space(&X0,&Y0,&X1,&Y1);
			order(&X0,&X1);
			order(&Y0,&Y1);
			(*do_space)(X0,Y0,X1,Y1);
			break;
		case 't':			/* text - place string so */
			getsf(s,fin);		/* first corner lies on   */
			(*do_label)(s);		/* current point	  */
			break;
		case 'v':			/* polygon fill */
			is_long = 0;
			get_poly(&n,&x,&y);
			for(i=0;i < n ; i++){
				x[i] = x[i] + NumX;
				y[i] = y[i] + NumY;
			}
			(*do_fillp)(n,x,y);
			free(x);
			free(y);
			break;
		case 'V':			/* polygon fill */
			is_long = 1;
			get_poly(&n,&x,&y);
			for(i=0;i < n ; i++){
				x[i] = x[i] + NumX;
				y[i] = y[i] + NumY;
			}
			(*do_fillp)(n,x,y);
			free(x);
			free(y);
			break;
		case 'W':			/* define line width */
			is_long = 0;
			get_wid(&Xi);
			(*do_gwid)(Xi);
			break;
		case 'w':			/* define line width */
			is_long = 1;
			get_wid(&Xi);
			(*do_gwid)(Xi);
			break;
		case 'x':			/* fill triangular region */
			is_long = 0;
			get_tri(&X0,&Y0,&X1,&Y1,&X2,&Y2,&patx,&paty,&lenx,&leny);
			(*do_fillt)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,X2+NumX,Y2+NumY,patx,paty,lenx,leny);
			break;
		case 'X':			/* fill triangular region */
			is_long = 1;
			get_tri(&X0,&Y0,&X1,&Y1,&X2,&Y2,&patx,&paty,&lenx,&leny);
			(*do_fillt)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,X2+NumX,Y2+NumY,patx,paty,lenx,leny);
			break;
		case 'y':			/* shade rectangular region*/
			is_long = 0;
			get_rect(&X0,&Y0,&X1,&Y1,&patx,&paty,&lenx,&leny);
			(*do_fillr)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,patx,paty,lenx,leny);
			break;
		case 'Y':			/* shade rectangular region*/
			is_long = 1;
			get_rect(&X0,&Y0,&X1,&Y1,&patx,&paty,&lenx,&leny);
			(*do_fillr)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,patx,paty,lenx,leny);
			break;
		case 'z':			/* seismic line shading */
			is_long = 0;
			get_sei(&X0,&Y0,&ixy,&istnd,&iplmn);
			(*do_fills)(X0+NumX,Y0+NumY,ixy,istnd,iplmn);
			break;
		case 'Z':			/* seismic line shading */
			is_long = 1;
			get_sei(&X0,&Y0,&ixy,&istnd,&iplmn);
			(*do_fills)(X0+NumX,Y0+NumY,ixy,istnd,iplmn);
			break;
		}
	}
	if(gphopen == 1)clospl(0);
}
static INT getpi(FILE *fin)  /* get an integer stored in 2/4 ascii bytes */
{
	short a, b;
	long la, lb, lc, ld;
	if(is_long == 0){
		if((a = getc(fin)) == EOF)
			return(EOF);
		if((b = getc(fin)) == EOF)
			return(EOF);
		b = b<<8;
		return(b|a);
	} else {
		if((la = getc(fin)) == EOF) return(EOF);
		if((lb = getc(fin)) == EOF) return(EOF);
		if((lc = getc(fin)) == EOF) return(EOF);
		if((ld = getc(fin)) == EOF) return(EOF);
		ld =ld << 24;
		lc =lc << 16;
		lb =lb <<  8;
		return(ld|lc|lb|la);
	}
}

static void getsf(char *s,FILE *fin) 		/* read in a string 
						terminated by a newline */
{
	for( ; *s = (char)getc(fin); s++)
		if(*s == '\n')
			break;
	*s = '\0';
	return;
}


static void ngetsf(char *s,FILE *fin, INT n)		/* read in a string 
						terminated by a newline */
				/* number of characters - 
					we only flag negative with a newline */
{
	int i;
	char *sp;
	sp = s;
	*s = 0;
	if(n < 0){	/* single character followed by newline */
		*s++ = (char)getc(fin);
		*s = (char)getc(fin); /* here get newline */
	} else {
		for( ; *s = (char)getc(fin); s++)
			if(*s == '\n')
				break;
		/* check for printability */
		for(i=0;i<n;i++)
			if(sp[i] < 32 || sp[i] >= 127)sp[i] = ' ';
	}
	*s = '\0';
	return;
}


void get_arc(INT *Xi,INT *Yi,INT *X0,INT *Y0,INT *X1,INT *Y1)
{
			*Xi = getpi(fin);	/* integers = center */
			*Yi = getpi(fin);	/* second two = start */
			*X0 = getpi(fin);	/* last two = end */
			*Y0 = getpi(fin);
			*X1 = getpi(fin);
			*Y1 = getpi(fin);
}

void get_circle(INT *Xi,INT *Yi,INT *r)
{
			*Xi = getpi(fin);	/* two integers give center */
			*Yi = getpi(fin);	/* third gives radius	*/
			*r  = getpi(fin);
}

void get_line(INT *X0,INT *Y0,INT *X1,INT *Y1)
{
			*X0 = getpi(fin);	/* two integers to second */
			*Y0 = getpi(fin);	/* two integers		*/
			*X1 = getpi(fin);
			*Y1 = getpi(fin);
}
void get_move(INT *Xi,INT *Yi)
{
			*Xi = getpi(fin);	/* give new current point */
			*Yi = getpi(fin);	/* this is a silent move, */
}
void get_cont(INT *Xi,INT *Yi)
{
			*Xi = getpi(fin);	/* give new current point */
			*Yi = getpi(fin);	/* this is a silent move, */
}
void get_point(INT *Xi,INT *Yi)
{
			*Xi = getpi(fin);	/* give new current point */
			*Yi = getpi(fin);	/* this is a silent move, */
}
void get_space(INT *X0,INT *Y0,INT *X1,INT *Y1)
{
			*X0 = getpi(fin);	/* this changes aspect ratio */
			*Y0 = getpi(fin);	/* two integers give lower */
			*X1 = getpi(fin);	/* left corner, next two */
			*Y1 = getpi(fin);	/* right corner	*/
}
void get_tri(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *X2,
		INT *Y2,INT *patx,INT *paty,INT *lenx,INT *leny)
{
			*X0 = getpi(fin);	/* 3 x,y pairs of trianlge */
			*Y0 = getpi(fin);
			*X1 = getpi(fin);
			*Y1 = getpi(fin);
			*X2 = getpi(fin);
			*Y2 = getpi(fin);
			*patx = getpi(fin);	/* fill pattern x-aXis */
			*paty = getpi(fin);	/* fill pattern y-aXis */
			*lenx = getpi(fin);
			*leny = getpi(fin);
}
void get_rect(INT *X0,INT *Y0,INT *X1,INT *Y1,
			INT *patx,INT *paty,INT *lenx,INT *leny)
{
			*X0 = getpi(fin);	/* x,y lower left corner */
			*Y0 = getpi(fin);
			*X1 = getpi(fin);
			*Y1 = getpi(fin);
			*patx = getpi(fin);	/* fill pattern x-aXis */
			*paty = getpi(fin);	/* fill pattern y-aXis */
			*lenx = getpi(fin);	/* fill pattern y-aXis */
			*leny = getpi(fin);	/* fill pattern y-aXis */
}
void get_sei(INT *X0, INT *Y0, INT *ixy, INT *istnd, INT *iplmn)
{
			*X0 = getpi(fin);	/* x,y coord first point */
			*Y0 = getpi(fin);
			*ixy= getpi(fin);	/* shade direction */
			*istnd=getpi(fin);	/* start , end	*/
			*iplmn=getpi(fin);	/* plus, neg amplitudes */
}
void get_sym(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *n,char *s)
{
			*X0 = getpi(fin);	/* x position */
			*Y0 = getpi(fin);	/* y position */
			*X1 = getpi(fin);	/* height */
			*Y1 = getpi(fin);	/* angle */
			*n  = getpi(fin);	/* number of characters */
			ngetsf(s,fin,*n); 		/* character string */
}
void get_wid(INT *Xi)
{
			*Xi = getpi(fin);
}

void get_poly(INT *n,INT **x,INT **y)
{
	int i;
	INT xx, yy;
	int ok = 1;
	*n = getpi(fin);
	if((*x = (INT *) calloc(*n, sizeof(INT))) ==NULL)
		ok = 0;
	if((*y = (INT *) calloc(*n, sizeof(INT))) ==NULL)
		ok = 0;
	for(i=0; i< *n ; i++){
		xx = getpi(fin);
		yy = getpi(fin);
		if(ok){
			(*x)[i] = xx;
			(*y)[i] = yy;
		}
	}
	if(!ok)*n = 0;
}
	
void get_clip(INT *reset, INT *X0,INT *Y0,INT *X1,INT *Y1)
{
			*reset = getpi(fin);	
			*X0 = getpi(fin);	/* one corner */
			*Y0 = getpi(fin);	
			*X1 = getpi(fin);	/* other corner */
			*Y1 = getpi(fin);	
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


void openpl(int showmenu)
{
	/* set up function pointers for plot_dbg */
	do_arc 		= di_arc ;
	do_circle 	= di_circle ;
	do_clip 	= di_clip ;
	do_cont 	= di_cont ;
	do_erase 	= di_erase ;
	do_fillp 	= di_fillp ;
	do_fillr 	= di_fillr ;
	do_fills 	= di_fills ;
	do_fillt 	= di_fillt ;
	do_font 	= di_font ;
	do_gottxt 	= di_gottxt ;
	do_cursor 	= di_cursor ;
	do_gsymb 	= di_gsymb ;
	do_gwid 	= di_gwid ;
	do_label 	= di_label ;
	do_linec 	= di_linec ;
	do_linemod 	= di_linemod ;
	do_move 	= di_move ;
	do_pen	 	= di_pen ;
	do_point 	= di_point ;
	do_space 	= di_space ;

	di_openpl(1);
}


void clospl(int mode)
{
	di_closepl(mode);
}

