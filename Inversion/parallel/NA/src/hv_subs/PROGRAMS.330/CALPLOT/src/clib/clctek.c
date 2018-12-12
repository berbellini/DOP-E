/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: CLCPLT                                                c
c                                                                     c
c      COPYRIGHT (C)  1997 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      3507 Laclede Avenue                                            c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
FILE *grfp;	/* this must be seen outside */
int is_long;
#include "dpsubs.h"
#include "disubs.h"
#ifdef MSDOS
#include <io.h>
#include <fcntl.h>
#define GETPID rand
#else
#define GETPID getpid
#endif
#include "calplot.h"

/* signal handling */
#include <signal.h>
#ifdef	MSDOS
	static int sigs[]= { SIGINT, SIGTERM };
#else
	static int sigs[]= { SIGHUP, SIGINT , SIGQUIT, SIGBUS, SIGKILL,
		SIGABRT, SIGTERM, SIGSEGV, SIGFPE, SIGILL };
#endif
extern void onintr(int arg);
static int sig_set = 0;

int dc_sleeptime = -1;

/* get function definitions */
/* set up pointers to graphics functions 
	that do the work */

void 	(*do_arc)()=0,
	(*do_circle)()=0,
	(*do_clip)()=0,
	(*do_closepl)()=0,
	(*do_cont)()=0,
	(*do_control)()=0,
	(*do_cross)()=0,
	(*do_erase)()=0,
	(*do_fillp)()=0,
	(*do_fillr)()=0,
	(*do_fills)()=0,
	(*do_fillt)()=0,
	(*do_font)()=0,
	(*do_gintxt)()=0,
	(*do_gottxt)()=0,
	(*do_cursor)()=0,
	(*do_info)()=0,
	(*do_gsymb)()=0,
	(*do_gwid)()=0,
	(*do_label)()=0,
	(*do_linec)()=0,
	(*do_linemod)()=0,
	(*do_move)()=0,
	(*do_openpl)()=0,
	(*do_pen)()=0,
	(*do_point)()=0,
	(*do_space)()=0;

#define YES 1
#define NO  0

static int plotfile = NO;
static int openplot = NO;
static int openinter = NO;
INT dp_ClipRightSave, dp_ClipLeftSave, dp_ClipTopSave, dp_ClipBottomSave;

static void dgclip(INT icmd, INT ilw, INT jlw, INT iup, INT jup);
static void order(INT *x,INT *y);

void clospl(int mode){
	if(plotfile == YES){
        	fflush(grfp);
		fclose(grfp);
		openplot = NO;
	} else {
		(*do_closepl)(mode);
		openinter = NO;
	}
	if(plotfile == YES)plotfile = NO;
}


#define NAMLEN  80
char dv_icon_name[NAMLEN];
void openpl(int ls,char *s, int lis, char *iconstr)
{
        FILE *fopen();
        char name[NAMLEN];
	int i;
	/* set up signals for termination gracefully */
	if(sig_set == 0){
		for ( i=0; i< sizeof(sigs)/sizeof(int) ; ++i)
			if(signal( sigs[i],SIG_IGN) != SIG_IGN )
				signal ( sigs[i], onintr);
		sig_set = 1;
	}
	/* set up global for clipping */
	dp_ClipRightSave	=  100000000;
	dp_ClipLeftSave		= -100000000;
	dp_ClipTopSave		=  100000000;
	dp_ClipBottomSave	= -100000000;
	/* parse the name of the file string */
	for(i=0 ; i< lis && *iconstr != ' ' && *iconstr != '\0' ; i++)
		dv_icon_name[i] = *iconstr++ ;
	dv_icon_name[i] = '\0' ;
	if(strncmp(s,"INTE",4) == 0 ){
		if(openplot == YES && plotfile == YES)clospl(0);
		plotfile = NO;
	} else {
		
		if(strncmp(s,"plot",4) == 0 ){	
#ifdef MSDOS
			int stime;
			long ltime;
			ltime = time(NULL);
			stime = (unsigned int)ltime/2;
			srand(stime);
			stime = rand()%9999;
			srand(stime);
        		sprintf(name,"plot%d",rand()%9999 );
#else
        		sprintf(name,"plot%d",GETPID() );
#endif
			if(openplot == NO){
#ifdef MSDOS
        		grfp = fopen(name,"wb");
#else
        		grfp = fopen(name,"w");
#endif
			openplot = YES;
			}
		} else if(strncmp(s,"stdout",6) == 0){
			grfp = stdout;
#ifdef MSDOS
			setmode(fileno(stdout),O_BINARY);
#endif
		} else {
			/* strip off blanks at the end of the string s */
			for(i=0 ; i< ls && *s != ' ' && *s != '\0' ; i++)
				name[i] = *s++ ;
			name[i] = '\0' ;
			if(openplot == NO){
#ifdef MSDOS
        		grfp = fopen(name,"wb");
#else
        		grfp = fopen(name,"w");
#endif
			openplot = YES;
			}
		}
		plotfile = YES;
	}



	/* set up the function definitions 
		for an interactive/hardcopy
		this is where the switch would occur */
	if(plotfile == YES){
		do_arc 		= dp_arc ;
		do_circle 	= dp_circle ;
		do_clip 	= dgclip ;
		do_cont 	= dp_cont ;
		do_control 	= dp_control ;
		do_cross 	= dp_cross ;
		do_erase 	= dp_erase ;
		do_fillp 	= dp_fillp ;
		do_fillr 	= dp_fillr ;
		do_fills 	= dp_fills ;
		do_fillt 	= dp_fillt ;
		do_font 	= dp_font ;
		do_gintxt 	= dp_gintxt ;
		do_gottxt 	= dp_gottxt ;
		do_cursor 	= dp_cursor ;
		do_info 	= dp_info ;
		do_gsymb 	= dp_gsymb ;
		do_gwid 	= dp_gwid ;
		do_label 	= dp_label ;
		do_linec 	= dp_linec ;
		do_linemod 	= dp_linemod ;
		do_move 	= dp_move ;
		do_pen	 	= dp_pen ;
		do_point 	= dp_point ;
		do_space 	= dp_space ;
	} else {
		do_arc 		= di_arc ;
		do_circle 	= di_circle ;
		do_clip 	= dgclip ;
		do_closepl 	= di_closepl ;
		do_cont 	= di_cont ;
		do_control 	= di_control ;
		do_cross 	= di_cross ;
		do_erase 	= di_erase ;
		do_fillp 	= di_fillp ;
		do_fillr 	= di_fillr ;
		do_fills 	= di_fills ;
		do_fillt 	= di_fillt ;
		do_font 	= di_font ;
		do_gintxt 	= di_gintxt ;
		do_gottxt 	= di_gottxt ;
		do_cursor 	= di_cursor ;
		do_info 	= di_info ;
		do_gsymb 	= di_gsymb ;
		do_gwid 	= di_gwid ;
		do_label 	= di_label ;
		do_linec 	= di_linec ;
		do_linemod 	= di_linemod ;
		do_move 	= di_move ;
		do_openpl 	= di_openpl ;
		do_pen	 	= di_pen ;
		do_point 	= di_point ;
		do_space 	= di_space ;

		/* now invoke the device open plot but
		only open it once to avoid screen clear */
		if(openinter == NO ){
			if(strncmp(s,"INTER",5) == 0 )
				(*do_openpl)(1);
			else if(strncmp(s,"INTEM",5) == 0 )
				(*do_openpl)(0);
			openinter = YES;
		}
	}
}

#define NAMSTR  81
static char name[NAMSTR];
/* output text string to graphics terminal at current position */
void gottxt(int cnt,char *s)
{
	int i,j;
	j = cnt;
	if(j >= NAMSTR)j=NAMSTR-1;
	strncpy(name,s,j);
	name[j] = '\0';
	(*do_gottxt)(name);
}
/* get text string from terminal */
void gintxt(int cnt, char *s)
{
	(*do_gintxt)( cnt, s);
}

void gomesg(int cnt,char *s)
{
	int i,j;
	if(plotfile == NO){
		j = cnt;
		if(j >= NAMSTR)j=NAMSTR-1;
		strncpy(name,s,j);
		name[j] = '\0';
		dv_mesg(name);
	}
}

void gocontrol(int type, int i1, int i2, int i3, int i4)
{
	(do_control)(type, i1, i2, i3 ,i4);
}

static void dgclip(INT icmd, INT ilw, INT jlw, INT iup, INT jup)
{
	if(icmd == 1){
		order(&ilw, &iup);
		order(&jlw, &jup);
	}
	dp_ClipRightSave	= iup;
	dp_ClipLeftSave		= ilw;
	dp_ClipTopSave		= jup;
	dp_ClipBottomSave	= jlw;
	if(plotfile == YES)
		dp_clip (icmd, ilw, jlw, iup, jup);
	else
		di_clip (icmd, ilw, jlw, iup, jup);
	if(icmd == 0){
		dp_ClipRightSave	=  100000000;
		dp_ClipLeftSave		= -100000000;
		dp_ClipTopSave		=  100000000;
		dp_ClipBottomSave	= -100000000;
	}
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
