/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      ROUTINE: GREAD                                               c
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
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef	MSDOS
#include <fcntl.h>
#include <io.h>
#define GETPID rand
#else
#define GETPID getpid
#include <unistd.h>
#endif
#include "dosubs.h"

#ifdef	MSDOS
#define INT long
#else
#define INT int
#endif

extern void (*do_arc)(INT Xi,INT Yi,INT X0,INT Y0,INT X1,INT Y1);
extern void (*do_circle)(INT Xi,INT Yi,INT r);
extern void (*do_clip)(INT cmd, INT X0,INT Y0,INT X1,INT Y1);
extern void (*do_cont)(INT X0,INT Y0);
extern void (*do_cross)(int *X0,int *Y0, char *c);
extern void (*do_cursor)(INT curstyp);
extern void (*do_erase)(INT mode);
extern void (*do_fillp)(INT n,INT *x,INT *y);
extern void (*do_fillr)(INT X0,INT Y0,INT X1,INT Y1,
	INT patx,INT paty,INT lenx,INT leny);
extern void (*do_fills)(INT X0,INT Y0,INT ixy,INT istnd,INT iplmn);
extern void (*do_fillt)(INT X0,INT Y0,INT X1,INT Y1,INT X2,INT Y2,
	INT patx,INT paty,INT lenx,INT leny);
extern void (*do_font)(INT Xi);
extern void (*do_gintxt)(int cnt, char *s);
extern void (*do_gottxt)(char *s);
extern void (*do_gsymb)(INT X0,INT Y0,INT X1,INT Y1,INT n,char *s);
extern void (*do_gwid)(INT wid);
extern void (*do_info)(INT *HasMouse, INT *XminDev, INT *YminDev, 
	INT *XmaxDev, INT *YmaxDev, INT *XminClip, 
	INT *YminClip, INT *XmaxClip, INT *YmaxClip, INT *Color);
extern void (*do_label)(char *s);
extern void (*do_linec)(INT X0,INT Y0,INT X1,INT Y1);
extern void (*do_linemod)(char *s);
extern void (*do_move)(INT X0,INT Y0);
extern void (*do_pen)(INT Xi);
extern void (*do_point)(INT X0,INT Y0);
extern void (*do_space)(INT X0,INT Y0,INT X1,INT Y1);


/* internal function definitions */
static void get_arc(INT *Xi,INT *Yi,INT *X0,INT *Y0,INT *X1,INT *Y1);
static void get_circle(INT *Xi,INT *Yi,INT *r);
static void get_clip(INT *reset, INT *X0,INT *Y0,INT *X1,INT *Y1);
static void get_line(INT *X0,INT *Y0,INT *X1,INT *Y1);
static void get_move(INT *Xi,INT *Yi);
static void get_cont(INT *Xi,INT *Yi);
static void get_point(INT *Xi,INT *Yi);
static void get_space(INT *X0,INT *Y0,INT *X1,INT *Y1);
static void get_poly(INT *n,INT **x,INT **y);
static void get_wid(INT *Xi);
static void get_tri(INT *X0,INT *Y0,INT *X1,INT *Y1,
		INT *X2,INT *y2,INT *patx,INT *paty,INT *lenx,INT *leny);
static void get_rect(INT *X0,INT *Y0,INT *X1,INT *Y1,
		INT *patx,INT *paty,INT *lenx,INT *leny);
static void get_sei(INT *X0,INT *Y0,INT *ixy,INT *istnd,INT *iplmn);
static void get_sym(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *n,char *s);
static void order (INT *x, INT *y);
static INT getpi(FILE *);
static void getsf(char *c, FILE *f);
static void ngetsf(char *c, FILE *f, INT n);
static void ClipRegion( INT ill1, INT jll1, INT iur1, INT jur1,
 		INT ill2, INT jll2, INT iur2, INT jur2,
 		INT *ill3, INT *jll3, INT *iur3, INT *jur3);

static void fplt(FILE *);

/* define extreme limits of plot space - the numbers are not 2^31
	but 1,000,000 inches is a lot */
static INT	MaxX = 1000000000;
static INT	MaxY = 1000000000;
static INT	MINY = 0;
static INT	MINX = 0;

static int	count	= 1 ;	/* current page number */
static int	kount	= 0 ;	/* count of characters in input file.
			   used only to skip an initial sequence of
			   erase characters */
static int	gphopen	= 0 ;	/* flag to guarantee that output file is
			   only opened once. This is useful for
			   plot merging	*/
char *fname;
static int	Num	=1 ;
static INT	NumX	=0 ;
static INT	NumY	=0 ;
static INT	HighX	=  1000000000 ;
static INT	HighY	=  1000000000 ;
static INT	LowX	= -1000000000 ;
static INT	LowY	= -1000000000 ;
static INT	defHighX	=  1000000000 ;
static INT	defHighY	=  1000000000 ;
static INT	defLowX		= -1000000000 ;
static INT	defLowY		= -1000000000 ;
static int	PlotStdout = 1 ;	/* 1 -> output is plot file, 0 -> stdout */
static int	std;
static int	fpltout = 0; 	/* counter on merged files so only  one space */
static INT	ill3, jll3, iur3, jur3;
static float	Sclx	= 1.0;
static float	Scly	= 1.0;



static FILE	*ffin	;	/* file for input data stream 	*/

static int   argcc;
static char **argvv;

extern INT dp_ClipRightSave, dp_ClipLeftSave, dp_ClipTopSave, dp_ClipBottomSave;

#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)

void dv_gread(int lf,char *fname, INT NumX, INT NumY, INT LowX, INT LowY, 
	INT HighX, INT HighY, INT Num, INT tSclx, INT tScly);
/*
	This routine maintains the functionality of reframe(I)
	except that no file merging is permitted
	and clipping uses current local clippng bound rather than
	device bound
*/
#define NAMSTR  81
static char name[NAMSTR];

void dv_gread(int lf,char *fname, INT tNumX, INT tNumY, INT tLowX, INT tLowY, 
	INT tHighX, INT tHighY, INT tNum, INT tSclx, INT tScly)
{
	int j;
	INT Tmp;
	NumX = tNumX;
	NumY = tNumY;
	LowX = tLowX;
	HighX = tHighX;
	LowY = tLowY;
	HighY = tHighY;
	Num = tNum;
	Sclx = (float)tSclx/1000.0;
	Scly = (float)tScly/1000.0;
	/* copy the string in */
	j = lf;
	if(j >= NAMSTR)j=NAMSTR-1;
	strncpy(name,fname,j);
	name[j] = '\0';

	/* get correct ordering for clipping */
	defHighX = dp_ClipRightSave;
	defHighY = dp_ClipTopSave;
	defLowX = dp_ClipLeftSave;
	defLowY = dp_ClipBottomSave;
	order(&LowX,&HighX);
	order(&LowY,&HighY);
	Tmp = MIN(defHighY,(HighY*Scly+NumY));
	defHighY = Tmp;
	Tmp = MIN(defHighX,(HighX*Sclx+NumX));
	defHighX = Tmp;
	Tmp = MAX(defLowY,(LowY*Scly+NumY));
	defLowY = Tmp;
	Tmp = MAX(defLowX,(LowX*Sclx+NumX));
	defLowX = Tmp;
	if(strlen(name)){
		if((ffin = fopen(name, "r")) == NULL) {
			return ;
		}
		fplt(ffin);
		fclose(ffin);
	}
	/* reset clipping */
	defHighX = dp_ClipRightSave;
	defHighY = dp_ClipTopSave;
	defLowX = dp_ClipLeftSave;
	defLowY = dp_ClipBottomSave;
	(*do_clip)((INT)1,defLowX, defLowY, defHighX, defHighY);
	
}

static int is_long;
FILE *fin;

static void fplt(dfin)
FILE *dfin;
{
	int c;
	int i;
	char s[256];
	INT Xi,Yi,X0,Y0,X1,Y1,X2,y2,patx,paty,ixy,istnd,iplmn,r,n;
	INT reset;
	INT lenx, leny;
	INT *x, *y;

	fin = dfin;
#ifdef	MSDOS
	setmode(fileno(fin), O_BINARY);
#endif
	count = 1;
	while((c=getc(fin)) != EOF){
		if(gphopen == 0){
			gphopen = 1;
		}
		switch(c){
		case 'a':			/* draw arc - first */
			is_long = 0;
			get_arc(&Xi,&Yi,&X0,&Y0,&X1,&Y1);
			if(count==Num)(*do_arc)(Xi+NumX,Yi+NumY,X0+NumX,Y0+NumY,X1+NumX,Y1);
			break;
		case 'A':
			is_long = 1;
			get_arc(&Xi,&Yi,&X0,&Y0,&X1,&Y1);
			if(count==Num)(*do_arc)(Xi+NumX,Yi+NumY,X0+NumX,Y0+NumY,X1+NumX,Y1);
			break;
		case 'B':			/* define font integer */
			is_long = 0;
			Xi = getpi(fin);
			if(count==Num)(*do_font)(Xi);
			break;
		case 'c':			/* circle	- first */
			is_long = 0;
			get_circle(&Xi,&Yi,&r);
			if(count==Num)(*do_circle)(Xi+NumX,Yi+NumY,r);
			break;
		case 'C':			/* 4 byte CIRCLE  */
			is_long = 1;
			get_circle(&Xi,&Yi,&r);
			if(count==Num)(*do_circle)(Xi+NumX,Yi+NumY,r);
			break;
		case 'D':			/* implement symbol */
			is_long = 0;
			get_sym(&X0,&Y0,&X1,&Y1,&n,s);
			if(count==Num)(*do_gsymb)(X0+NumX,Y0+NumY,X1,Y1,n,s);
			break;
		case 'd':			/* implement symbol */
			is_long = 1;
			get_sym(&X0,&Y0,&X1,&Y1,&n,s);
			if(count==Num)(*do_gsymb)(X0+NumX,Y0+NumY,X1,Y1,n,s);
			break;
		case 'e':			/* erase plot frame */
			if(kount!=1)count++;
			if(count > Num)goto cnt;
			break;
		case 'E':			/* erase plot frame */
			if(kount!=1)count++;
			if(count > Num)goto cnt;
			break;
		case 'f':			/* read linemod string */
			getsf(s,fin);
			if(count==Num)(*do_linemod)(s);
			break;
		case 'g':			/* set cursor */
			is_long = 0;
			Xi = getpi(fin);
			(*do_cursor)(Xi);
			break;
		case 'h':			/* define color integer */
			is_long = 0;
			Xi = getpi(fin);
			if(count==Num)(*do_pen)(Xi);
			break;
		case 'k':			/* clipping */
			is_long = 0;
			get_clip(&reset,&X0,&Y0,&X1,&Y1);
			order(&X0,&X1);
			order(&Y0,&Y1);
			if(count==Num){
				if(reset == 0){
				(*do_clip)((INT)1,defLowX, defLowY,
					defHighX, defHighY);
				} else {
				ClipRegion( X0+NumX , Y0+NumY, X1+NumX, Y1+NumY,
					defLowX, defLowY, defHighX, defHighY,
 					&ill3, &jll3, &iur3, &jur3);
				(*do_clip)(reset,ill3,jll3,iur3,jur3);
				}
			}
			break;
		case 'K':			/* clipping */
			is_long = 1;
			get_clip(&reset,&X0,&Y0,&X1,&Y1);
			order(&X0,&X1);
			order(&Y0,&Y1);
			if(count==Num){
				if(reset == 0){
				(*do_clip)((INT)1,defLowX, defLowY,
					defHighX, defHighY);
				} else {
				ClipRegion( X0+NumX , Y0+NumY, X1+NumX, Y1+NumY,
					defLowX, defLowY, defHighX, defHighY,
 					&ill3, &jll3, &iur3, &jur3);
				(*do_clip)(reset,ill3,jll3,iur3,jur3);
				}
			}
			break;
		case 'l':			/* draw line from first */
			is_long = 0;
			get_line(&X0,&Y0,&X1,&Y1);
			if(count==Num)(*do_linec)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY);
			break;
		case 'L':			/* 4 byte LINE */
			is_long = 1;
			get_line(&X0,&Y0,&X1,&Y1);
			if(count==Num)(*do_linec)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY);
			break;
		case 'm':			/* move - two integers */
			is_long = 0;
			get_move(&Xi,&Yi);
			if(count==Num)(*do_move)(Xi+NumX,Yi+NumY);	/* pen up */ 		
			break;
		case 'M':			/* 4 byte MOVE */
			is_long = 1;
			get_move(&Xi,&Yi);
			if(count==Num)(*do_move)(Xi+NumX,Yi+NumY);	/* pen up */
			break;
		case 'n':			/* continue - two integers */
			is_long = 0;
			get_cont(&Xi,&Yi);
			if(count==Num)(*do_cont)(Xi+NumX,Yi+NumY);	
				/* new becomes current */
			break;
		case 'N':			/* 4 byte CONT */
			is_long = 1;
			get_cont(&Xi,&Yi);
			if(count==Num)(*do_cont)(Xi+NumX,Yi+NumY);	
				/* new becomes current */
			break;
		case 'o':			/* output text to terminal */
			getsf(s,fin); 		/* at current position */
			if(count==Num)(*do_gottxt)(s);
			break;
		case 'p':			/* draw a point at the two */
			is_long = 0;
			get_point(&Xi,&Yi);
			if(count==Num)(*do_point)(Xi+NumX,Yi+NumY);
			break;
		case 'P':			/* 4 byte POINT */
			is_long = 1;
			get_point(&Xi,&Yi);
			if(count==Num)(*do_point)(Xi+NumX,Yi+NumY);
			break;
		case 's':			/* space - define plot space */
			is_long = 0;
			get_space(&X0,&Y0,&X1,&Y1);
			order(&X0,&X1);
			order(&Y0,&Y1);
			if(fpltout == 0 ){
				(*do_space)(X0,Y0,X1,Y1);
				/* IF THIS IS THE INITIAL TIME SET CLIPPING */
				reset = 1;
				(*do_clip)(reset,defLowX,defLowY,
					defHighX,defHighY);
			}
			break;
		case 'S':			/* 4 byte SPACE */
			is_long = 1;
			get_space(&X0,&Y0,&X1,&Y1);
			order(&X0,&X1);
			order(&Y0,&Y1);
			if(fpltout == 0 ){
				(*do_space)(X0,Y0,X1,Y1);
				/* IF THIS IS THE INITIAL TIME SET CLIPPING */
				reset = 1;
				(*do_clip)(reset,defLowX,defLowY,
					defHighX,defHighY);
			}
			break;
		case 't':			/* text - place string so */
			getsf(s,fin);		/* first corner lies on   */
			if(count==Num)(*do_label)(s);		
					/* current point	  */
			break;
		case 'v':			/* polygon fill */
			is_long = 0;
			get_poly(&n,&x,&y);
			for(i=0;i < n ; i++){
				x[i] = x[i] + NumX;
				y[i] = y[i] + NumY;
			}
			if(count==Num)(*do_fillp)(n,x,y);
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
			if(count==Num)(*do_fillp)(n,x,y);
			free(x);
			free(y);
			break;
		case 'W':			/* define line width */
			is_long = 0;
			get_wid(&Xi);
			if(count==Num)(*do_gwid)(Xi);
			break;
		case 'w':			/* define line width */
			is_long = 1;
			get_wid(&Xi);
			if(count==Num)(*do_gwid)(Xi);
			break;
		case 'x':			/* fill triangular region */
			is_long = 0;
			get_tri(&X0,&Y0,&X1,&Y1,&X2,&y2,&patx,&paty,&lenx,&leny);
			if(count==Num)(*do_fillt)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,
				X2+NumX,y2+NumY,patx,paty,lenx,leny);
			break;
		case 'X':			/* fill triangular region */
			is_long = 1;
			get_tri(&X0,&Y0,&X1,&Y1,&X2,&y2,&patx,&paty,&lenx,&leny);
			if(count==Num)(*do_fillt)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,
				X2+NumX,y2+NumY,patx,paty,lenx,leny);
			break;
		case 'y':			/* shade rectangular region*/
			is_long = 0;
			get_rect(&X0,&Y0,&X1,&Y1,&patx,&paty,&lenx,&leny);
			if(count==Num)(*do_fillr)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,
				patx,paty,lenx,leny);
			break;
		case 'Y':			/* shade rectangular region*/
			is_long = 1;
			get_rect(&X0,&Y0,&X1,&Y1,&patx,&paty,&lenx,&leny);
			if(count==Num)(*do_fillr)(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,
				patx,paty,lenx,leny);
			break;
		case 'z':			/* seismic line shading */
			is_long = 0;
			get_sei(&X0,&Y0,&ixy,&istnd,&iplmn);
			if(count==Num)(*do_fills)(X0+NumX,Y0+NumY,ixy,istnd,iplmn);
			break;
		case 'Z':			/* seismic line shading */
			is_long = 1;
			get_sei(&X0,&Y0,&ixy,&istnd,&iplmn);
			if(count==Num)(*do_fills)(X0+NumX,Y0+NumY,ixy,istnd,iplmn);
			break;
		}
	}
cnt: 
	/* reset clipping */
	reset = 0;
	(*do_clip)(reset,defLowX,defLowY,defHighX,defHighY);
}
static INT getpi(fin)  /* get an integer stored in 2 ascii bytes */
FILE *fin;	  /* this is inefficient but hardware independent */
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
	for( ; *s = (char )getc(fin); s++)
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
	int nread;
	if(n < 0){	/* single character followed by newline */
		*s++ = (char )getc(fin);
		*s = (char )getc(fin); /* here get newline */
		*s = '\0';
	} else {
		for(nread=0 ; s[nread] = (char )getc(fin);  nread++, nread < n){
			if(s[nread] == '\n'){
				break;
			}
		}
		s[nread] = '\0';
		/* check for printability */
		for(i=0;i<nread;i++)
			if(s[i] < 32 || s[i] >= 127)s[i] = ' ';
	}
	return;
}


static void get_arc(INT *Xi,INT *Yi,INT *X0,INT *Y0,INT *X1,INT *Y1)
{
	*Xi = (INT)(Sclx * getpi(fin));	/* integers = center */
	*Yi = (INT)(Scly * getpi(fin));	/* second two = start */
	*X0 = (INT)(Sclx * getpi(fin));	/* last two = end */
	*Y0 = (INT)(Scly * getpi(fin));
	*X1 = (INT)(Sclx * getpi(fin));
	*Y1 = (INT)(Scly * getpi(fin));
}

static void get_circle(INT *Xi,INT *Yi,INT *r)
{
	*Xi = (INT)(Sclx * getpi(fin));	/* two integers give center */
	*Yi = (INT)(Scly * getpi(fin));	/* third gives radius	*/
	*r  = (INT)(Sclx * getpi(fin));
}

static void get_line(INT *X0,INT *Y0,INT *X1,INT *Y1)
{
	*X0 = (INT)(Sclx * getpi(fin));	/* two integers to second */
	*Y0 = (INT)(Scly * getpi(fin));	/* two integers		*/
	*X1 = (INT)(Sclx * getpi(fin));
	*Y1 = (INT)(Scly * getpi(fin));
}
static void get_move(INT *Xi,INT *Yi)
{
	*Xi = (INT)(Sclx * getpi(fin));	/* give new current point */
	*Yi = (INT)(Scly * getpi(fin));	/* this is a silent move, */
}
static void get_cont(INT *Xi,INT *Yi)
{
	*Xi = (INT)(Sclx * getpi(fin));	/* give new current point */
	*Yi = (INT)(Scly * getpi(fin));	/* this is a silent move, */
}
static void get_point(INT *Xi,INT *Yi)
{
	*Xi = (INT)(Sclx * getpi(fin));	/* give new current point */
	*Yi = (INT)(Scly * getpi(fin));	/* this is a silent move, */
}
static void get_space(INT *X0,INT *Y0,INT *X1,INT *Y1)
{
	*X0 = getpi(fin);	/* this changes aspect ratio */
	*Y0 = getpi(fin);	/* two integers give lower */
	*X1 = getpi(fin);	/* left corner, next two */
	*Y1 = getpi(fin);	/* right corner	*/
}
static void get_tri(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *X2,
		INT *Y2,INT *patx,INT *paty,INT *lenx,INT *leny)
{
	*X0 = (INT)(Sclx * getpi(fin));	/* 3 x,y pairs of trianlge */
	*Y0 = (INT)(Scly * getpi(fin));
	*X1 = (INT)(Sclx * getpi(fin));
	*Y1 = (INT)(Scly * getpi(fin));
	*X2 = (INT)(Sclx * getpi(fin));
	*Y2 = (INT)(Scly * getpi(fin));
	*patx = getpi(fin);	/* fill pattern x-aXis */
	*paty = getpi(fin);	/* fill pattern y-aXis */
	*lenx = (INT)(Sclx * getpi(fin));
	*leny = (INT)(Scly * getpi(fin));
}
static void get_rect(INT *X0,INT *Y0,INT *X1,INT *Y1,
			INT *patx,INT *paty,INT *lenx,INT *leny)
{
	*X0 = (INT)(Sclx * getpi(fin));	/* x,y lower left corner */
	*Y0 = (INT)(Scly * getpi(fin));
	*X1 = (INT)(Sclx * getpi(fin));
	*Y1 = (INT)(Scly * getpi(fin));
	*patx = getpi(fin);	/* fill pattern x-aXis */
	*paty = getpi(fin);	/* fill pattern y-aXis */
	*lenx = (INT)(Sclx * getpi(fin));	/* fill pattern y-aXis */
	*leny = (INT)(Scly * getpi(fin));	/* fill pattern y-aXis */
}
static void get_sei(INT *X0, INT *Y0, INT *ixy, INT *istnd, INT *iplmn)
{
	*X0 = (INT)(Sclx * getpi(fin));	/* x,y coord first point */
	*Y0 = (INT)(Scly * getpi(fin));
	*ixy= getpi(fin);	/* shade direction */
	*istnd=getpi(fin);	/* start , end	*/
	*iplmn=getpi(fin);	/* plus, neg amplitudes */
}
static void get_sym(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *n,char *s)
{
	*X0 = (INT)(Sclx * getpi(fin));	/* x position */
	*Y0 = (INT)(Scly * getpi(fin));	/* y position */
	*X1 = (INT)(Sclx * getpi(fin));	/* height */
	*Y1 = getpi(fin);	/* angle */
	*n  = getpi(fin);	/* number of characters */
	ngetsf(s,fin,*n); 		/* character string */
}
static void get_wid(INT *Xi)
{
	*Xi = (INT)(Sclx * getpi(fin));
}

static void get_poly(INT *n,INT **x,INT **y)
{
	int i;
	INT xx, yy;
	int ok = 1;
	*n =  getpi(fin);
	if((*x = (INT *) calloc(*n, sizeof(INT))) ==NULL)
		ok = 0;
	if((*y = (INT *) calloc(*n, sizeof(INT))) ==NULL)
		ok = 0;
	for(i=0; i< *n ; i++){
		xx = (INT)(Sclx * getpi(fin));
		yy = (INT)(Scly * getpi(fin));
		if(ok){
			(*x)[i] = xx;
			(*y)[i] = yy;
		}
	}
	if(!ok)*n = 0;
}
	
static void get_clip(INT *reset,INT *X0,INT *Y0,INT *X1,INT *Y1)
{
	*reset = getpi(fin);	
	*X0 = (INT)(Sclx * getpi(fin));	/* one corner */
	*Y0 = (INT)(Scly * getpi(fin));	
	*X1 = (INT)(Sclx * getpi(fin));	/* other corner */
	*Y1 = (INT)(Scly * getpi(fin));	
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




static void ClipRegion( INT ill1, INT jll1, INT iur1, INT jur1,
 		INT ill2, INT jll2, INT iur2, INT jur2,
 		INT *ill3, INT *jll3, INT *iur3, INT *jur3)
{
	/* determine the smallest intersection of two
		rectangular regions */
	/* for safety ensure that the rectangles are in fact defined
		by LOWER LEFT and UPPER RIGHT */
	order(&ill1,&iur1);
	order(&jll1,&jur1);
	order(&ill2,&iur2);
	order(&jll2,&jur2);
	
	*ill3 = MAX(ill1, ill2);
	*jll3 = MAX(jll1, jll2);
	*iur3 = MIN(iur1, iur2);
	*jur3 = MIN(jur1, jur2);
	/* if the rectangles overlap, then these must be the
		lower left and upper right */
	if( *ill3 > *iur3 || *jll3 > *jur3 ){
		/* set a useless default */
		*ill3 = -1;
		*jll3 = -1;
		*iur3 = -1;
		*jur3 = -1;
	}
}
