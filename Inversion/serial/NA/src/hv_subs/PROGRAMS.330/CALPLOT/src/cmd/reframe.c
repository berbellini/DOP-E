/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: REFRAME                                               c
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
#include <string.h>
#include <signal.h>
#include <stdlib.h>
#include <stddef.h>
#define GETPID getpid
#include <unistd.h>
#include "dpsubs.h"

	static int sigs[]= { SIGHUP, SIGINT , SIGQUIT, SIGTERM };
#define INT int

extern void 	dp_arc(INT X1, INT Y1, INT X2, INT Y2, INT x3, INT z3),
	dp_circle(),
	dp_clip(),
	dp_cont(),
	dp_erase(),
	dp_fillp(),
	dp_fillr(),
	dp_fills(),
	dp_fillt(),
	dp_font(),
	dp_gottxt(),
	dp_gsymb(),
	dp_gwid(),
	dp_label(),
	dp_linec(),
	dp_linemod(),
	dp_move(),
	dp_pen(),
	dp_point(),
	dp_space();

void	openpl( int showmenu ),
		di_gcmdln( int, char ** ),
		clospl( int mode );

/* internal function definitions */
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
		INT *X2,INT *y2,INT *patx,INT *paty,INT *lenx,INT *leny);
void get_rect(INT *X0,INT *Y0,INT *X1,INT *Y1,
		INT *patx,INT *paty,INT *lenx,INT *leny);
void get_sei(INT *X0,INT *Y0,INT *ixy,INT *istnd,INT *iplmn);
void get_sym(INT *X0,INT *Y0,INT *X1,INT *Y1,INT *n,char *s);
void order (INT *x, INT *y);
static INT getpi(FILE *);
static void getsf(char *c, FILE *f);
static void ngetsf(char *c, FILE *f, INT n);
static void ClipRegion( INT ill1, INT jll1, INT iur1, INT jur1,
 		INT ill2, INT jll2, INT iur2, INT jur2,
 		INT *ill3, INT *jll3, INT *iur3, INT *jur3);

static void fplt(FILE *);

/* define extreme limits of plot space - the numbers are not 2^31
	but 1,000,000 inches is a lot */
INT	MaxX = 1000000000;
INT	MaxY = 1000000000;
INT	MINY = 0;
INT	MINX = 0;

int	count	= 1 ;	/* current page number */
int	kount	= 0 ;	/* count of characters in input file.
			   used only to skip an initial sequence of
			   erase characters */
int	gphopen	= 0 ;	/* flag to guarantee that output file is
			   only opened once. This is useful for
			   plot merging	*/
char	Merge[100];
int	Num	=1 ;
INT	NumX	=0 ;
INT	NumY	=0 ;
INT	HighX	=  1000000000 ;
INT	HighY	=  1000000000 ;
INT	LowX	= -1000000000 ;
INT	LowY	= -1000000000 ;
INT	defHighX	=  1000000000 ;
INT	defHighY	=  1000000000 ;
INT	defLowX		= -1000000000 ;
INT	defLowY		= -1000000000 ;
int	PlotStdout = 1 ;	/* 1 -> output is plot file, 0 -> stdout */
int	std;
int	fpltout = 0; 	/* counter on merged files so only  one space */
static INT      ill3, jll3, iur3, jur3;

FILE	*ffin	;	/* file for input data stream 		*/
FILE	*grfp;			/* file for output plot stream		*/


static int   argcc;
static char **argvv;


#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)



int main(int argc,char **argv)
{
	int i;
	INT Tmp;
	void onintr();
	argcc = argc; /* save command arguments for other routines */
	argvv = argv;
	ffin = stdin ;
	for ( i=0; i< sizeof(sigs)/sizeof(int) ; ++i)
		if(signal( sigs[i],onintr) == SIG_IGN )
			signal ( sigs[i], SIG_IGN);

	di_gcmdln(argc,argv);
	/* get correct ordering for clipping */
	order(&LowX,&HighX);
	order(&LowY,&HighY);
	Tmp = MIN(defHighY,HighY+NumY);
	defHighY = Tmp;
	Tmp = MIN(defHighX,HighX+NumX);
	defHighX = Tmp;
	Tmp = MAX(defLowY,LowY+NumY);
	defLowY = Tmp;
	Tmp = MAX(defLowX,LowX+NumX);
	defLowX = Tmp;
	fplt(ffin);
	fpltout++;
	if(strlen(Merge)){
		if((ffin = fopen(Merge, "r")) == NULL) {
			fprintf(stderr,"can't open %s\n",Merge);
			exit(1);
		}
		fplt(ffin);
		fpltout++;
	}
	return(0);
}

int is_long;
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
	while((c=getc(fin)) != EOF){
		if(gphopen == 0){
			gphopen = 1;
			openpl(1);
		}
		switch(c){
		case 'a':			/* draw arc - first */
			is_long = 0;
			get_arc(&Xi,&Yi,&X0,&Y0,&X1,&Y1);
			if(count==Num)dp_arc(Xi+NumX,Yi+NumY,X0+NumX,Y0+NumY,X1+NumX,Y1);
			break;
		case 'A':
			is_long = 1;
			get_arc(&Xi,&Yi,&X0,&Y0,&X1,&Y1);
			if(count==Num)dp_arc(Xi+NumX,Yi+NumY,X0+NumX,Y0+NumY,X1+NumX,Y1);
			break;
		case 'B':			/* define font integer */
			is_long = 0;
			Xi = getpi(fin);
			if(count==Num)dp_font(Xi);
			break;
		case 'c':			/* circle	- first */
			is_long = 0;
			get_circle(&Xi,&Yi,&r);
			if(count==Num)dp_circle(Xi+NumX,Yi+NumY,r);
			break;
		case 'C':			/* 4 byte CIRCLE  */
			is_long = 1;
			get_circle(&Xi,&Yi,&r);
			if(count==Num)dp_circle(Xi+NumX,Yi+NumY,r);
			break;
		case 'D':			/* implement symbol */
			is_long = 0;
			get_sym(&X0,&Y0,&X1,&Y1,&n,s);
			if(count==Num)dp_gsymb(X0+NumX,Y0+NumY,X1,Y1,n,s);
			break;
		case 'd':			/* implement symbol */
			is_long = 1;
			get_sym(&X0,&Y0,&X1,&Y1,&n,s);
			if(count==Num)dp_gsymb(X0+NumX,Y0+NumY,X1,Y1,n,s);
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
			if(count==Num)dp_linemod(s);
			break;
		case 'g':			/* set cursor */
			is_long = 0;
			Xi = getpi(fin);
			dp_cursor(Xi);
			break;
		case 'h':			/* define color integer */
			is_long = 0;
			Xi = getpi(fin);
			if(count==Num)dp_pen(Xi);
			break;
		case 'k':			/* clipping */
			is_long = 0;
			get_clip(&reset,&X0,&Y0,&X1,&Y1);
			order(&X0,&X1);
			order(&Y0,&Y1);
			if(count==Num){
				if(reset == 0){
				dp_clip((INT)1,defLowX, defLowY,
					defHighX, defHighY);
				} else {
				ClipRegion( X0+NumX , Y0+NumY, X1+NumX, Y1+NumY,
					defLowX, defLowY, defHighX, defHighY,
 					&ill3, &jll3, &iur3, &jur3);
				dp_clip(reset,ill3,jll3,iur3,jur3);
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
				dp_clip((INT)1,defLowX, defLowY,
					defHighX, defHighY);
				} else {
				ClipRegion( X0+NumX , Y0+NumY, X1+NumX, Y1+NumY,
					defLowX, defLowY, defHighX, defHighY,
 					&ill3, &jll3, &iur3, &jur3);
				dp_clip(reset,ill3,jll3,iur3,jur3);
				}
			}
			break;
		case 'l':			/* draw line from first */
			is_long = 0;
			get_line(&X0,&Y0,&X1,&Y1);
			if(count==Num)dp_linec(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY);
			break;
		case 'L':			/* 4 byte LINE */
			is_long = 1;
			get_line(&X0,&Y0,&X1,&Y1);
			if(count==Num)dp_linec(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY);
			break;
		case 'm':			/* move - two integers */
			is_long = 0;
			get_move(&Xi,&Yi);
			if(count==Num)dp_move(Xi+NumX,Yi+NumY);		/* pen up	*/
			break;
		case 'M':			/* 4 byte MOVE */
			is_long = 1;
			get_move(&Xi,&Yi);
			if(count==Num)dp_move(Xi+NumX,Yi+NumY);		/* pen up	*/
			break;
		case 'n':			/* continue - two integers */
			is_long = 0;
			get_cont(&Xi,&Yi);
			if(count==Num)dp_cont(Xi+NumX,Yi+NumY);	
				/* new becomes current */
			break;
		case 'N':			/* 4 byte CONT */
			is_long = 1;
			get_cont(&Xi,&Yi);
			if(count==Num)dp_cont(Xi+NumX,Yi+NumY);	/* new becomes current */
			break;
		case 'o':			/* output text to terminal */
			getsf(s,fin); 		/* at current position */
			if(count==Num)dp_gottxt(s);
			break;
		case 'p':			/* draw a point at the two */
			is_long = 0;
			get_point(&Xi,&Yi);
			if(count==Num)dp_point(Xi+NumX,Yi+NumY);
			break;
		case 'P':			/* 4 byte POINT */
			is_long = 1;
			get_point(&Xi,&Yi);
			if(count==Num)dp_point(Xi+NumX,Yi+NumY);
			break;
		case 's':			/* space - define plot space */
			is_long = 0;
			get_space(&X0,&Y0,&X1,&Y1);
			order(&X0,&X1);
			order(&Y0,&Y1);
			if(fpltout == 0  ){
				dp_space(X0,Y0,X1,Y1);
				/* IF THIS IS THE INITIAL TIME SET CLIPPING */
				reset = 1;
				dp_clip(reset,defLowX,defLowY,
					defHighX,defHighY);
			}
			break;
		case 'S':			/* 4 byte SPACE */
			is_long = 1;
			get_space(&X0,&Y0,&X1,&Y1);
			order(&X0,&X1);
			order(&Y0,&Y1);
			if(fpltout == 0  ){
				dp_space(X0,Y0,X1,Y1);
				/* IF THIS IS THE INITIAL TIME SET CLIPPING */
				reset = 1;
				dp_clip(reset,defLowX,defLowY,
					defHighX,defHighY);
			}
			break;
		case 't':			/* text - place string so */
			getsf(s,fin);		/* first corner lies on   */
			if(count==Num)dp_label(s);		
					/* current point	  */
			break;
		case 'v':			/* polygon fill */
			is_long = 0;
			get_poly(&n,&x,&y);
			for(i=0;i < n ; i++){
				x[i] = x[i] + NumX;
				y[i] = y[i] + NumY;
			}
			if(count==Num)dp_fillp(n,x,y);
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
			if(count==Num)dp_fillp(n,x,y);
			free(x);
			free(y);
			break;
		case 'W':			/* define line width */
			is_long = 0;
			get_wid(&Xi);
			if(count==Num)dp_gwid(Xi);
			break;
		case 'w':			/* define line width */
			is_long = 1;
			get_wid(&Xi);
			if(count==Num)dp_gwid(Xi);
			break;
		case 'x':			/* fill triangular region */
			is_long = 0;
			get_tri(&X0,&Y0,&X1,&Y1,&X2,&y2,&patx,&paty,&lenx,&leny);
			if(count==Num)dp_fillt(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,X2+NumX,y2+NumY,patx,paty,lenx,leny);
			break;
		case 'X':			/* fill triangular region */
			is_long = 1;
			get_tri(&X0,&Y0,&X1,&Y1,&X2,&y2,&patx,&paty,&lenx,&leny);
			if(count==Num)dp_fillt(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,X2+NumX,y2+NumY,patx,paty,lenx,leny);
			break;
		case 'y':			/* shade rectangular region*/
			is_long = 0;
			get_rect(&X0,&Y0,&X1,&Y1,&patx,&paty,&lenx,&leny);
			if(count==Num)dp_fillr(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,patx,paty,lenx,leny);
			break;
		case 'Y':			/* shade rectangular region*/
			is_long = 1;
			get_rect(&X0,&Y0,&X1,&Y1,&patx,&paty,&lenx,&leny);
			if(count==Num)dp_fillr(X0+NumX,Y0+NumY,X1+NumX,Y1+NumY,patx,paty,lenx,leny);
			break;
		case 'z':			/* seismic line shading */
			is_long = 0;
			get_sei(&X0,&Y0,&ixy,&istnd,&iplmn);
			if(count==Num)dp_fills(X0+NumX,Y0+NumY,ixy,istnd,iplmn);
			break;
		case 'Z':			/* seismic line shading */
			is_long = 1;
			get_sei(&X0,&Y0,&ixy,&istnd,&iplmn);
			if(count==Num)dp_fills(X0+NumX,Y0+NumY,ixy,istnd,iplmn);
			break;
		}
	}
cnt:
	/* reset clipping */
	reset = 0;
	dp_clip(reset,defLowX,defLowY,defHighX,defHighY);
	if(gphopen == 1)clospl(0);
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
	if(n < 0){	/* single character followed by newline */
		*s++ = (char )getc(fin);
		*s = (char )getc(fin); /* here get newline */
	} else {
		for( ; *s = (char )getc(fin); s++)
			if(*s == '\n')
				break;
		/* check for printability */
		for(i=0;i<n;i++)
			if(s[i] < 32 || s[i] >= 127)s[i] = ' ';
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
		INT *y2,INT *patx,INT *paty,INT *lenx,INT *leny)
{
			*X0 = getpi(fin);	/* 3 x,y pairs of trianlge */
			*Y0 = getpi(fin);
			*X1 = getpi(fin);
			*Y1 = getpi(fin);
			*X2 = getpi(fin);
			*y2 = getpi(fin);
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
	
void get_clip(INT *reset,INT *X0,INT *Y0,INT *X1,INT *Y1)
{
			*reset = getpi(fin);	
			*X0 = getpi(fin);	/* one corner */
			*Y0 = getpi(fin);	
			*X1 = getpi(fin);	/* other corner */
			*Y1 = getpi(fin);	
}

void di_gcmdln(int argc,char **argv)
{
	int verbos = 0;
	char *cp, *tp;
	while(argc-- > 1 ) {
		if(*argv[1] == '-'){
			switch(argv[1][1]){
			case 'V':
				fprintf(stderr,"CALPLOT (2.1) COPYRIGHT (C) 1989 Saint Louis University\n");
				verbos = 1;
				break;
			case 'M':
				cp = argv[1];
				cp++;
				cp++;
				strncpy(Merge,cp,strlen(cp));
				break;
			case 'N':
				Num = atoi(&argv[1][2]);
				break;
			case 'X':
				tp = &argv[1][3];
				if(argv[1][2] == '0'){
					NumX = (INT)atoi(tp);
				}
				else if(argv[1][2] == 'L'){
					LowX = (INT)atoi(tp);
				}
				else if(argv[1][2] == 'H'){
					HighX = (INT)atoi(tp);
				}
				order(&LowX,&HighX);
				break;
			case 'Y':
				tp = &argv[1][3];
				if(argv[1][2] == '0'){
					NumY = (INT)atoi(tp);
				}
				else if(argv[1][2] == 'L'){
					LowY = (INT)atoi(tp);
				}
				else if(argv[1][2] == 'H'){
					HighY = (INT)atoi(tp);
				}
				order(&LowY,&HighY);
				break;
			case 'P':
				PlotStdout = 1;
				break;
			case 'O':
				PlotStdout = 0;
				break;
			default:
				break;
			}
		}
		else {
			if((ffin = fopen(argv[1], "r")) == NULL) {
				fprintf(stderr,"can't open %s\n",argv[1]);
				exit(1);
			}
		}
		argv++;
	}
	if(verbos){
		fprintf(stderr,"HighX %d\n",HighX);
		fprintf(stderr,"HighY %d\n",HighY);
		fprintf(stderr,"LowX %d\n",LowX);
		fprintf(stderr,"LowY %d\n",LowY);
		fprintf(stderr,"NumX %d\n",NumX);
		fprintf(stderr,"NumY %d\n",NumY);
		fprintf(stderr,"Num  %d\n",Num );
	}
}

void order(INT *x,INT *y)
{
	INT tmp;
	if(*y < *x){
		tmp = *x;
		*x = *y;
		*y = tmp;
	}
}

void openpl( int showmenu )
{
	char	name[256];
	if(PlotStdout) {
		sprintf(name,"plot%d",GETPID() );
		grfp = fopen(name,"w");
	}
	else 
		grfp = stdout;
#ifdef MSDOS
	setmode(fileno(grfp),O_BINARY);
#endif
}

void clospl( int mode ){
	fflush(grfp);
}

void onintr()
{
	exit( 1 );
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
