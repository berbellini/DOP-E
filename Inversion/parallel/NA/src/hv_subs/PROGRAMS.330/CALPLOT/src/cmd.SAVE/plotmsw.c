/*
Command line priorities on geometry (highest first)
	In .Xdefaults
	plotxvig.calxvig.plotxvig.geometry: 1000x800+100+50 for standalone
		program e.g.,plotxvig
	plotxvig.calxvig.MFTMENU96.geometry: 1000x800+10+10 for
		library based program using
		ginitf("INTER","MFTMENU96")

	Then 
		if standalone
			COMMAND LINE
			ENVIRONMENT
		if library
			ENVIRONMENT
	
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLOTXVIG                                              c
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
CHANGES
	10 Jan 2001 - removed unused calls to dv_filltri and dv_fillbox
	13 Jan 2001 - introduced -p flag for positioning grid
	02 APR 2004 - inplemented gend(mode) at higher level which
		required change in dv_closepl
*/
/* device dependent plot program interface 			*/
/* X11 using XviG library					*/
/* It is assumed that the plotgen program clips and requires	*/
/*	the very low level primitives				*/
/*	supports the calls 					*/
/*		dv_clip						*/
/*		dv_closepl					*/
/*		dv_erase					*/
/*		dv_fillp					*/
/*		dv_font						*/
/*		di_gcmdln					*/
/*		di_gsymb					*/
/*		onintr						*/
/*		dv_openpl					*/
/*		dv_pen						*/
/*		zpoint						*/
/*		dv_rliner						*/

#ifndef INT
#ifdef MSDOS
#define INT long
#else
#define INT int
#endif
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<ctype.h>
#include	<signal.h>

#define IX(n)	( dc_xscale *  dc_scalefac * ( n + dc_xoff ) + 0.5)
#define IY(n)   ( dc_yscale *  dc_scalefac * ( n + dc_yoff ) + 0.5)
#define LX(n)	( dc_xscale *  dc_scalefac * ( n  ) + 0.5)
#define LY(n)   ( dc_yscale *  dc_scalefac * ( n  ) + 0.5)
#define INVX(n) ( (n - 0.5)/(dc_xscale*dc_scalefac) - dc_xoff)
#define INVY(n) ( (n - 0.5)/(dc_yscale*dc_scalefac) - dc_yoff)
#define MAX(a,b) ( (a) > (b) ? (a):(b) )
#define MIN(a,b) ( (b) > (a) ? (a):(b) )
#define ON      1
#define OFF     0
#define ACTION_NEXT     0
#define ACTION_INVERT   1
#define ACTION_GRAY     2
#define ACTION_QUIT     3


/* Global Function Prototypes	*/
void dv_clip();
void di_cross(int *ix,int *iy,char *s);
void dv_closepl(int mode);
void dv_concur(INT x, INT y);
void dv_erase(INT mode);
void dv_fillp(INT n, INT *x,INT *y);
void di_gottxt(char *s);
void di_gsymb(INT x,INT y,INT ht,INT ang,INT nchar,char *s);
extern INT dv_lineclip(INT *x0,INT *z0,INT *x1,INT *z1,
	INT MinX,INT MinY,INT MaxX,INT MaxY);
void dv_mesg(char *mesg);
void dv_movcur(INT x, INT y);
void dv_openpl(int showmenu);
void dv_pen(INT p);
void dv_rliner(INT x0,INT z0,INT x1,INT z1);
extern void dv_symvec(INT xloc,INT yloc,INT ht,char *s,INT ang,INT nochar);
void  dv_zpoint(int x,int y);
void dv_zzpoint(int x,int y);

/* Global Variables		*/

extern double dc_xscale ;
extern double dc_yscale ;
extern double dc_xoff   ;
extern double dc_yoff   ;
INT dc_left	= 0;	/* device plotting limits */
INT dc_right	= 600;
INT dc_bottom	= 0;
INT dc_top	= 480;
double dc_Nxpi	= 409.0;
double dc_Nypi	= 409.0;
INT dc_oldx1	= -1;
INT dc_oldy1 	= -1;
int dc_curpen = 1 ;	/* current pen value */
extern int	dc_sleeptime ;
int	dc_mapcurpen;	/* mapping into dither for dv_zzpoint */
int	scaledplot = OFF;
INT	dc_rotate = OFF;
INT	dc_shadeoff = OFF ;	/* permit shading unless turned off 
				in command line */
INT	dc_color = ON ;		/* hardware fast color shading */
INT	dc_ColorInfo = 2;
INT	dc_xlinewidth = 0;
INT	dc_ylinewidth = 0;
INT	dc_linewidth = 0;
INT	dc_herelinewidth = 1; /* set linewidth here */
INT	dc_newlinewidth = 0;	/* flag to indicate a line width change */
INT	dc_hardwarefill = 1;	/* flag to have all filling done by hardware */
				/* 0 raster shade here
                                   1 hardware fill
                                   2 hardware fill and hardware seismic fill
				*/
INT	dc_curfont = 0;
int	defaultfont = 0;

double	dc_scalefac = 1.0;
int dc_hasmouse = 0;
extern INT dc_oldx;	/* current plot point position */
extern INT dc_oldy;	/* needed by gottxt() and gintxt() */
INT dc_ClipRight, dc_ClipTop, dc_ClipBottom, dc_ClipLeft;
INT dc_ClipRight_sv, dc_ClipTop_sv, dc_ClipBottom_sv, dc_ClipLeft_sv;
INT	dc_iseps = 0;	/* invoke special scaling for EPS */
#define NAMLEN  80
extern char dv_icon_name[NAMLEN];

/* Local  Variables		*/

static short xpos,ypos,xposold,yposold;
       int Ldc_top, Ldc_right, Ldc_bottom, Ldc_left; /* for internal clipping */
static int savpen = 1;	/* current color */
static int dc_dim2 = 255;	/* dimensions of dither matrix dv_zzpoint */
static int xoring = OFF;	/* 0 regular plot 1 XOR for cursor */
static int black = 1;	/* back ground is white */
static int kolor = 1;   /* do not use halftone to emulate color */
static int Kolor = 1;
static int maxcolors;	/* either 2 for monochrome or 16
				this is the number of different lines */
static double true_dc_scalefac = 1.0;
static int revvideo = 0;/* 0 gottxt straight write out strips trailing blank */
			/* 1 gottxt writes blanks then characters */
			/* 2 reversed video with trailing blanks */
static int cross_inc = 10;
static int curstype;
static int border;
static int title;
static int lastbuttonx;         /* X-value of last button */
static int Button_Color_Dark  = 1100;
static int Button_Color       = 1090;
static int Button_Color_Light = 1070;
static int Button_Color_Fore  =    1;
static int local_pen = 1;
static int dogrid = 0;
static int dogrid_inc = 1;

#define NAMSTR  256
static char name[NAMSTR];
static char omesg[NAMSTR];

static INT mx, my;

#define MAXCOL 35
int	Maxcol = MAXCOL;	/* 35 unique colors */
float   red_array[MAXCOL];
float green_array[MAXCOL];
float  blue_array[MAXCOL];
int	CBLK = 0;	/* map position for preset colors */
int	CWHIT= 1;
int	CRED = 1;
int	CORNG = 1;
int	CYEL = 1;
int	CGRN = 1;
int	CBLGR= 1;
int	CBLU = 1;
static int colmap[MAXCOL][3] ;

static int user_curstype;

#define EVENT_NONE			0
#define EVENT_MOUSE_MOTION		1
#define EVENT_BUTTON_LEFT		2
#define EVENT_BUTTON_RIGHT		3
#define EVENT_BUTTON_CENTER		4
#define EVENT_KEY_PRESS			5
#define EVENT_ARROW_LEFT		6
#define EVENT_ARROW_RIGHT		7
#define EVENT_ARROW_UP			8
#define EVENT_ARROW_DOWN		9
#define EVENT_PAGE_UP			10
#define EVENT_PAGE_DOWN			11

static int char_pressed;
static int getevent();
static int gcmdln_read = OFF;   /* flag to determine if di_gcmdln called */
static int gcmdln_change = OFF;
static  int eargc;
static  char    **eargv = NULL;
static int *coords;
static int Width=-1, Height=-1, Xpos=-1, Ypos=-1;

/* Local  Function Prototypes	*/

static void draw_border(int xl, int yl, int xh, int yh);
static void draw_grid();
static void draw_dash(int k, INT ix0, INT iy0,INT ix1, INT iy1);
static int check_menu(int xpos, int ypos);
static void show_menu(void);

static void clearscr();
static void prgerr();
static void lsleep(int n);
static void cblank();
static void gottxtp(char *s);
static void draw_cursor(int curstyp, int xc, int yc, int xl, int yl, int xh, int yh);
static int inclipregion(int xpos, int ypos);
static void (*xline)(int y,int xs,int xe);
static void (*yline)(int x,int ys,int ye);
static void (*point)(int x,int y);
static void xvigxline(int y,int xs,int xe);
static void xvigyline(int x,int ys,int ye);
static void xvigpoint(int x,int y);
static void terminate();
static void gtenv();

static void setcolor(int cmd);
static void cal_coord(float *red,float *green,float *blue,int index);
static float value(float n1, float n2, float hue);
static void hls_to_rgb(float *r,float *g,float *b,float h,float l,float s);
static void svcolorindices( int ncolor,float red_array[],
		float green_array[],float blue_array[],int cmd);
static void set_line_index(int index);
static void fillquad(int x0,int z0,int x1,int z1,
	int x2, int z2, int x3, int z3);
static void do_parsegeom(char *str);
static int ReadInteger(char *string, char **NextString);
static void usage(void);

#include "xvig.h"



void onintr(int arg)
{
	/* go directly to the XviG_Exit() */
	XviG_Exit();
	/* give some information about the signal */
	switch (arg) {
		case SIGINT :
 			fprintf(stderr,"control C clean up\n");
			break;
		case SIGILL:
			fprintf(stderr,"SIGILL exit\n");
			break;
		case SIGABRT:
			fprintf(stderr,"SIGABRT exit\n");
			break;
		case SIGFPE:
			fprintf(stderr,"SIGFPE exit\n");
			break;
		case SIGSEGV:
			fprintf(stderr,"SIGSEGV exit\n");
			break;
		case SIGTERM:
			fprintf(stderr,"SIGTERM exit\n");
			break;
		default:
			fprintf(stderr,"Unknown signal=%d \n",arg);
	}
	exit(arg);
}





#include <time.h>
static void lsleep(int n)
{
	time_t t, ot;
	time(&ot);
	t = ot + n;
	while(ot < t)
		time(&ot);
}


void dv_zpoint(int x,int y)
{
	(*point)(x+border,Ldc_top-y-border);
}




/* routine to parse environment with fields separated by colons ':'
	The routine returns eargc which is the number of arguments 
	and eargv which is a char **. This is an array of strings
	for the various arguments. 

	To make this compatible with the C language argc and argv used in
	main(int argc ; char **argv)
	eargv is incremented by 1, and argv[0] = NULL
*/
static void gtenv()
{
	char *cp;
	int i,lstr;
	int ibeg, iend, iarg;
	char *eptr;
	char *tptr;
	int	*iarr;
	int icnt;

#ifdef MSW
	eptr = getenv("PLOTMSW");
#else
	eptr = getenv("PLOTXVIG");
#endif
	if(eptr != (char *)NULL){
		/* first ensure that the desired structure of
			beginning and ending with a ':' is true */
		lstr = strlen(eptr);
		tptr = calloc(lstr+3,sizeof(char));
		if(eptr[0] != ':')
			strcpy(tptr,":");
		strcat(tptr,eptr);
		if(eptr[lstr-1] != ':')
			strcat(tptr,":");
	

		/* get number of arguments */
		eargc = 0;
		lstr = strlen(tptr);
		for(i=0 ; i < lstr; i++){
			if(tptr[i] == ':'){
					eargc++;
			}
		}
		/* now allocate a temporary array for colon pointers */
		iarr = ( int *) calloc (eargc, sizeof(int ) );
		lstr = strlen(tptr);
		for(i=0,icnt=0;i< lstr;i++){
			if(tptr[i] == ':')
				iarr[icnt++] = i;
		}
				
		/* now allocate the number of strings */
		eargv = ( char **) calloc (eargc, sizeof(char * ) );
		/* now systematically go through to get the
			beginning and end for a copy */
		/* for the eargv[0] pointer */
		eargv[0] = calloc(1,sizeof(char));
		eargv[0][0] = '\0';
		for(i=1;i < icnt; i++){
			lstr = iarr[i] - iarr[i-1];
			eargv[i] = calloc(lstr,sizeof(char));
			strncpy(eargv[i],&tptr[iarr[i-1]+1],lstr-1);
			eargv[i][lstr] = '\0';
		}
		free(tptr);
	}
}


void di_gcmdln(int argc,char *argv[])
{
	char *cp;
	int int_val;
	gcmdln_read = ON;
	while(argc-- > 1 ) {
		cp = argv[1];
		if(*argv[1] == '-'){
			switch(argv[1][1]){
			case 'V':
				fprintf(stderr,"CALPLOT (3.0) COPYRIGHT (C) 1989 Saint Louis University\n");
				break;
			case 'S':
				scaledplot = ON;
				cp = argv[1];
				cp++;
				cp++;
				true_dc_scalefac = atof(cp);
				if(true_dc_scalefac <= 0.0)
					true_dc_scalefac = 1.0;
				dc_scalefac = true_dc_scalefac;
				gcmdln_change = ON;
				break;
			case 'W':
				cp = argv[1];
				cp++;
				cp++;
				dc_sleeptime = atoi(cp);
				gcmdln_change = ON;
				break;
			case 'R':
				dc_rotate = 1;
				gcmdln_change = ON;
				break;
			case 'N':
				kolor = 0;
				dc_shadeoff = ON;
				dc_hardwarefill = 0;
				gcmdln_change = ON;
				break;
			case 'G':
				/* VGA color */
				kolor = 2;
				gcmdln_change = ON;
				break;
			case 'K':
				kolor = 1;
				switch(argv[1][2]){
				case 'R':
					Kolor = 2;
					break;
				case 'B':
					Kolor = 3;
					break;
				default:
					break;
				}
				gcmdln_change = ON;
				break;
			case 'I' :
				/* make background black by inverting */
				black = 0;
				gcmdln_change = ON;
				break ;
			case 'F':	/* default font value */
				cp = argv[1];
				cp++;
				cp++;
				defaultfont = atoi(cp);
				if(defaultfont < 0)defaultfont = 0;
				dc_curfont = defaultfont;
				gcmdln_change = ON;
				break;
			case 'g':	/* parse a geometry e.g.,
				-geometry widthxheight{+-}xoff{+-}yoff
				we only need the -g */
				argv++;
				argc--;
				cp = argv[1];
				do_parsegeom(cp);
				gcmdln_change = ON;
				break;
			case 'p':	/* put up a positioning grid */
				dogrid = 1;
				cp = argv[1];
				cp++;
				cp++;
				dogrid_inc = 1;
				if(strlen(cp)> 0){
					int_val = atoi(cp);
					if(int_val < 0){
						dogrid_inc = 1;
					} else if(int_val > 10){
						dogrid_inc = 10;
					} else if(int_val < 10 && int_val > 4){
						dogrid_inc = 5;
					} else if(int_val == 3){
						dogrid_inc = 4;
					} else {
						dogrid_inc = int_val;
					}
				}
				break;
			case 'h':
			case '?':
				usage();
				break;
			default:
				break;
			}
		}
		argv++;
	}
	/* error check */
}


static void prgerr()
{
}




void dv_font(INT fontvalue)
{
	if(fontvalue > 0)
		dc_curfont = (fontvalue-1)%4 +1;
	else
		dc_curfont = defaultfont;
}


void di_gsymb(INT x,INT y,INT ht,INT ang,INT nchar,char *s)
{
	if(nchar > 0){
		dv_symvec(x,y,ht,s,ang,nchar);
	} else {
		dv_symvec(x,y,ht,s,ang,nchar);
	}
}




static int dim = 16;
static int bit[16][16] = {
/* 4x4 magic square dithered int 16x16 */
{  1,225, 49,209, 15,239, 63,223,  4,228, 52,212, 14,238, 62,222},
{177, 81,129, 97,191, 95,143,111,180, 84,132,100,190, 94,142,110},
{193, 33,241, 17,207, 47,255, 31,196, 36,244, 20,206, 46,254, 30},
{113,145, 65,161,127,159, 79,175,116,148, 68,164,126,158, 78,174},
{ 12,236, 60,220,  6,230, 54,214,  9,233, 57,217,  7,231, 55,215},
{188, 92,140,108,182, 86,134,102,185, 89,137,105,183, 87,135,103},
{204, 44,252, 28,198, 38,246, 22,201, 41,249, 25,199, 39,247, 23},
{124,156, 76,172,118,150, 70,166,121,153, 73,169,119,151, 71,167},
{ 13,237, 61,221,  3,227, 51,211, 16,240, 64,224,  2,226, 50,210},
{189, 93,141,109,179, 83,131, 99,192, 96,144,112,178, 82,130, 98},
{205, 45,253, 29,195, 35,243, 19,208, 48,256, 32,194, 34,242, 18},
{125,157, 77,173,115,147, 67,163,128,160, 80,176,114,146, 66,162},
{  8,232, 56,216, 10,234, 58,218,  5,229, 53,213, 11,235, 59,219},
{184, 88,136,104,186, 90,138,106,181, 85,133,101,187, 91,139,107},
{200, 40,248, 24,202, 42,250, 26,197, 37,245, 21,203, 43,251, 27},
{120,152, 72,168,122,154, 74,170,117,149, 69,165,123,155, 75,171}
};



void dv_zzpoint(int x,int y)
{
	static int ipen, jpen;
	INT i, ix, iy, ixm, iym;
	if(dc_curpen < 1000)
		dv_zpoint(x,y);
	else {
		if(bit[x%dim][y%dim] < dc_mapcurpen)
			dv_zpoint(x,y);
	}
}

void dv_cursor(INT curstyp)
{
	/* establish the cursor to use for interactive plot */
	/* curstyp	0	arrow
	        	1	XOR arrow
			2	XOR cross hair
			3	XOR plus
	*/
	user_curstype = (int) curstyp;
	if(user_curstype < 0)
		user_curstype = 0;
	else if(user_curstype > 6)
		user_curstype = 0;
}








/************************* XviG Stuff ********************************/
/*
#include "xvig.h"
*/

/* get text string from terminal */
void di_gintxt(int cnt,char *s)
{

int xx, yy;
	char ch[2], us[2];
	int key;
	int in_char, in_ch;
	int last=0;

	strcpy(name,"");
	ch[0] = '\000';
	ch[1] = '\000';
	us[0] = '\137';
	us[1] = '\000';
	/* force a wait for a RETURN perhaps permit ENTER LF? */
	if( last < cnt ){
		cblank();
		gottxtp(us);
		dc_oldx-=8;
	}
	while(1){
		/* grab new character */
		in_ch = XviG_GetChar();
		ch[0] = (in_ch ) &0xff;
		if(ch[0] == 13 || ch[0] == 12)goto out;
		/* add character if room */
		if( last < cnt && isprint(ch[0])){
			/* clear bit map and then put out an underscore */
			name[last] = ch[0];
			last++;
			cblank();	/* clear bits of past dregs */
			gottxtp(ch);
			/* clean up the keyboard buffer to avoid problems with
				other keyboard operations */
		}
		if( (ch[0] == 8) && (strlen(name) > 0) ){ /* error correct */
			if(last < cnt)cblank();
			if(last > 0 )dc_oldx-=8;
			cblank();
			name[last--] = 0;
			if(last <0)last=0;
		}
		if( last < cnt ){
			cblank();
			gottxtp(us);
			dc_oldx-=8;
		}
	}
out:
	/* get rid of last underscore */
	cblank();
	name[last] = '\0';
	if(last > cnt)last = cnt ;
	strncpy(s,name,last);
	s[last] = '\0';
}


/* set up color tables for this device */

static void setcolor(int cmd)
{
	float red, green, blue;
	int  index;
	for(index=0;index<Maxcol;index++){
		cal_coord(&red,&green,&blue,index);
		red_array[index] = red;
		green_array[index] = green;
		blue_array[index] = blue;
	}
	if(Maxcol > 2){
		CRED = 0;
		CBLU = Maxcol - 3;
		CGRN = (CBLU+CRED)/2;
		CYEL = (CGRN+CRED)/2;
		/* note linear interpolation does not work
			we need to do trig interpolation */
		CYEL = (int)( (float)(0.707*CGRN + CRED) + 0.5);
		CORNG= (int)( (float)(0.43*CGRN + CRED) + 0.5);
		CBLGR= (int)( (float)(CBLU - 0.707*CGRN) + 0.5);

		CBLK = 1;
		CWHIT = 0;
		/* correct for offset in colors */
		CRED+=2;
		CBLU+=2;
		CGRN+=2;
		CYEL+=2;
		CORNG+=2;
		CBLGR+=2;
/*
fprintf(stderr,"CRED   : %d\n",CRED);
fprintf(stderr,"CBLU   : %d\n",CBLU);
fprintf(stderr,"CGRN   : %d\n",CGRN);
fprintf(stderr,"CYEL   : %d\n",CYEL);
fprintf(stderr,"CORNG  : %d\n",CORNG);
fprintf(stderr,"CBLGR  : %d\n",CBLGR);
*/
/* EXACT FIX FOR MAXCOL 35 using X11 definitions
	orange  255 165   0
	yellow  255 255   0
        green     0 255   0
	red     255   0   0
	blue      0   0 255
        cyan      0 255 255
*/
CRED=2;
CYEL=12;
CGRN=18;
CBLU=34;
CORNG=7;
CBLGR=24;
	} else {
		CRED = 1; CBLU = 1; CGRN = 1; CYEL = 1; CORNG = 1;  CBLGR = 1;
		CBLK = 1; CWHIT = 0;
	}
/*
fprintf(stderr,"CRED   : %d\n",CRED);
fprintf(stderr,"CBLU   : %d\n",CBLU);
fprintf(stderr,"CGRN   : %d\n",CGRN);
fprintf(stderr,"CYEL   : %d\n",CYEL);
fprintf(stderr,"CORNG  : %d\n",CORNG);
fprintf(stderr,"CBLGR  : %d\n",CBLGR);
*/
	svcolorindices( Maxcol, red_array, green_array, blue_array,cmd);
}

static void cal_coord(float *red,float *green,float *blue,int index)
{
	float rred, rgrn, rblu;
	float l = 0.5;	/* lumination	*/
	float s = 1.0;	/* saturation	*/
	float h = 0.0 ; /* hue		*/
	float fac ;	/* map index to gray non linearly */
	int jindex;
	if(index == 0){
		rred = 0.0;
		rgrn = 0.0;
		rblu = 0.0;
	}
	else if(index == 1){
		rred = 1.0;
		rgrn = 1.0;
		rblu = 1.0;
	}
	else {
		jindex = index;
		if(kolor < 2){	/* actual color */
			if(Kolor == 1){
				/* construct scale Red -> Green -> Blue */
				l = 0.5;
				s = 1.0 ;
/*
				h = 240.0*(float)(jindex-2)/(float)(Maxcol-3);
*/
				h = (float)(jindex-2)/(float)(Maxcol-3);
				if(h <= 0.5)
					fac = 2.0*h*h;
				else
					fac = 1.0 - 2.0*(1.0-h)*(1.0-h);
				fac = (fac + h )/2.0;
				h = fac * 240.0;

				hls_to_rgb(&rred,&rgrn,&rblu,h,l,s);
			} else {
				/* construct scale Red -> Green -> Blue */
				fac = (float)(jindex-2)/(float)(Maxcol-3) ;
				if(Kolor == 3)fac = 1.0 - fac;
					/*Blue -> White -> Red */
				if(fac <=0.5){
					rred = 1.0;
					rgrn = 2.0*fac;
					rblu = 2.0*fac;
				} else {
					rred = 2.0*(1.0 - fac);
					rgrn = 2.0*(1.0 - fac);
					rblu = 1.0;
				}
			}
		} else {		/* grayscale */
			fac = (float)(jindex-2)/(float)(Maxcol-3) ;
			/* for PC combatibility */
			if(black==0)
				fac = 0.3 + 0.5*fac;
			else
				fac = 0.3 + 0.5*fac;
			rred = 1.0-fac;
			rgrn = 1.0-fac;
			rblu = 1.0-fac;
		}
	}
	*red = rred;
	*green= rgrn;
	*blue= rblu;
}

static float value(float n1, float n2, float hue)
{
	if(hue > 360.0)hue-=360.;
	if(hue < 0.0)hue+=360.;
	if(hue < 60.0){
		return( n1 + (n2-n1)*hue/60. );
	}
	else if(hue >= 60.0 && hue < 180.0){
		return( n2);
	}
	else if(hue >= 180.0 && hue < 240.0){
		return( n1 + (n2-n1)*(240.0-hue)/60.0 );
	}
	else if(hue >= 240.0){
		return (n1);
	}
}

/* Mapping of h(hue), l (lumination), s (saturation) to
   strength of red, green, blue
   Foley and Van Dam, Fundamentals of Interactive Computer Graphics
   1982, Addison-Wesley, page 619 */
/* s = 0 -> gray, s = 1.0 -> full strengh color
   l = 0 -> black, l=0.5 -> full color, l=1.0 -> white
   h = 0 -> red, h=30 -> orange, h=60 -> yellow, h=90-> yellow green
   h = 120->green, h=180->blue green, h=240->blue, h=300->magenta
   we only use 0 <= h <= 240 */

static void hls_to_rgb(float *r,float *g,float *b,float h,float l,float s)
{
	float m1, m2;
	if(l <= 0.5){
		m2 = l*(1. + s);
	}
	else
		m2 = l + s - l*s;
	m1 = 2*l - m2;
	if(s == 0){
		*r = l;
		*g = l;
		*b = l;
	} else {
		*r = value(m1,m2,h+120.0);
		*g = value(m1,m2,h);
		*b = value(m1,m2,h-120.0);
	}
}


static void svcolorindices( int ncolor,float red_array[],
		float green_array[],float blue_array[],int cmd)
{
	int i;
	for(i=2;i<ncolor;i++){
		colmap[i-1][0] = (255.0  * red_array[i]);
		colmap[i-1][1] = (255.0 * green_array[i]);
		colmap[i-1][2] = (255.0 * blue_array[i]);
	}
	/* foreground */
		colmap[ncolor-1][0] = 255*(1-black);
		colmap[ncolor-1][1] = 255*(1-black);
		colmap[ncolor-1][2] = 255*(1-black);
	/* background */
		colmap[0][0] = 255*black;
		colmap[0][1] = 255*black;
		colmap[0][2] = 255*black;
/* debug 
for(i=0;i < ncolor; i++)
	fprintf(stderr,"%5d %5d %5d %5d\n",i,colmap[i][0],
		colmap[i][1], colmap[i][2]);
*/

}

void dv_closepl(int mode)
{
	dv_erase(mode);
	terminate();
}

void dv_openpl(int showmenu)
{
	int numx, numy, i, j;
	int c;
	int savesleep;
	unsigned int win_width=600, win_height=480;
	unsigned int old_width=600, old_height=480;
	unsigned int new_win_width, new_win_height;
	int win_xpos, win_ypos;
	int depth;
	/* only get environment if command line has not been invoked */
        if(gcmdln_read == OFF || gcmdln_change == OFF ){
                gtenv();
                di_gcmdln(eargc,eargv);
        }
	/* always force interactive */
	dc_sleeptime = -1;
	/* initialize color map */
	setcolor(1);

	border = 8;
	title = 16;
	if(dv_icon_name[0] == '\0')
#ifndef MSW
		strcpy(dv_icon_name,"plotxvig");
#else
		strcpy(dv_icon_name,"plotmsw");
#endif
	if (!(depth=XviG_Init(dv_icon_name, colmap, MAXCOL))){
		exit(1);
	}
	if(depth == 1){
		maxcolors = 2;
		dc_hardwarefill = 1;
	} else {
		maxcolors = 8;
		dc_hardwarefill = 1;
	}
	/* dc_colorInfo 
		0	monochrome on white background
		1	gray on white background
		2	color on white background
		4	monochrome on black background
		5	gray on black background
		6	color on black background

	kolor	0
		1	Color     
		2	Grayscale
	Kolor	1	Red -> Green -> Blue
		2	Red -> White -> Blue
	*/

	if(maxcolors == 2){
		if(black == 0)
			dc_ColorInfo = 4;
		else
			dc_ColorInfo = 0;
	} else {
		if(black == 0){		/* background is black */
			dc_ColorInfo = 7 - kolor;
		} else {
			dc_ColorInfo = 3 - kolor; 
		}
	}
	/* define pointers to functions */
	point = xvigpoint;
	xline = xvigxline;
	yline = xvigyline;
	/* define window coordinates if not the default 
		from command line or environment override */
	if(Width > 0 && Height > 0){
		win_width = Width;
		old_width = win_width;
		win_height = Height;
		old_height = win_height;
	}


	XviG_OpenWindow(dv_icon_name,Xpos,Ypos,&win_width, &win_height);
	sleep(1);
	XviG_WindowSize(&new_win_width, &new_win_height);


#ifndef MSW
	XviG_WindowPosition(&win_xpos, &win_ypos);
        if(old_width != new_win_width && old_height != new_win_height){
                /* adjust for CALPLOT ASPECT RATIO */
                new_win_height = ( 8. * new_win_width ) / 10.;
                win_height = new_win_height + border + border + title;
                win_width = new_win_width + border + border;
                XviG_WindowPosition(&win_xpos, &win_ypos);
                XviG_CloseWindow(dv_icon_name);
                XviG_OpenWindow(dv_icon_name,win_xpos,win_ypos,
			&win_width, &win_height);
        }
#endif
        dc_top = win_height -1 ;
        dc_right = win_width -1 ;
	dc_bottom = 0;
	dc_left = 0 ;


	/* save the actual device limits */
	Ldc_top = (int)dc_top;
	Ldc_right = (int)dc_right;
	Ldc_left = (int)dc_left;
	Ldc_bottom = (int)dc_bottom;




	dc_right = dc_right - border -border;
	dc_top = dc_top - border - border - title;

	dc_Nxpi = (double)dc_right/10.00;
	dc_Nypi = (double)dc_top  / 8.00;

	dc_ClipLeft = dc_left;
	dc_ClipRight = dc_right;
	dc_ClipTop = dc_top;
	dc_ClipBottom = dc_bottom;
	dc_ClipLeft_sv = dc_left;
	dc_ClipRight_sv = dc_right;
	dc_ClipTop_sv = dc_top;
	dc_ClipBottom_sv = dc_bottom;
	/* check for mouse */
	dc_hasmouse = 1;

	/* Draw a border */
	xpos = (dc_right + dc_left)/2;
	ypos = (dc_top + dc_bottom)/2;
	omesg[0] = '\0';
	dv_cursor((INT) 0);
	
	savesleep = dc_sleeptime;
	/* force an automatic erase */
	dc_sleeptime = 0;
	dv_erase(0);
	dc_sleeptime = savesleep;

	/* send bounds */
	XviG_SendMessage( 1, border, title, 0, 0);
	/* send clip */
	XviG_SendMessage( 2, dc_ClipLeft, dc_ClipBottom, dc_ClipRight, dc_ClipTop);
	user_curstype = 0;


}

void dv_erase(INT mode)
{
	int chr, xx, yy;
	/* mode
		1 - automatically erase without waiting and
			without user intervention
	*/
	/* erase automatically if dc_sleeptime >= 0
		else wait for a key to be hit or left mouse
		except when cursor is in menu region */
	/* xpos, ypos measured from lower left border
		dc_left <= xpos <= dc_right
		dc_bottom <= ypos <= dc_top + title
	*/
	int menu_res;
	if(mode == 1 )goto junk;
	if(dogrid)draw_grid();
	if(dc_sleeptime >= 0){
		lsleep(dc_sleeptime);
	} else {
		XviG_Flush();
		XviG_SetCursor(XviG_CURSOR_ARROW);
		while(1){
			chr = XviG_GetCursor(XviG_BUTTON, &xx, &yy);
			/* convert from absolute screen coordinates to
				user coordinates
			user               screen
			(0,0)              (border,Ldc_top-border)
			(dc_right,dc_top)  (border+dc_right,border+title)
			(xpos,ypos)        (xx,yy)
			*/
			xpos = xx - border;
			ypos = Ldc_top - border - yy;
			if(ypos > dc_top){
				menu_res = check_menu(xpos,ypos);
				if(menu_res == ACTION_NEXT){
					goto enderase;
				} else if(menu_res == ACTION_QUIT){
					terminate();
				}
			} else {
				goto enderase;
			}
		}
	}
/*
*/
enderase:
junk:
	clearscr();
	dv_pen((INT)1);
	XviG_SetCursor(XviG_CURSOR_ARROW);
	/* Draw a border */
	draw_border(Ldc_left,Ldc_bottom,Ldc_right,Ldc_top);
}

void dv_pen(INT jpen)
{
	/* let data structure tell number of colors */
	int index;
	int ipen;
	ipen = jpen;
	switch(ipen){
		case 2000:
			revvideo = 0;
			return;
		case 2001:
			revvideo = 1;
			return;
		case 2002:
			revvideo = 2;
			return;
		case 3000:
			xoring = 0;
			XviG_SetGC(0);
			return;
		case 3001:
			xoring = 1;
			XviG_SetGC(1);
			return;
		case 3002:
			xoring = 1 - xoring;
			if(xoring == 1)
				XviG_SetGC(1);
			else
				XviG_SetGC(0);
			return;
		default:
			break;
	}
	dc_curpen = ipen;
	if(dc_curpen >= 1000 && dc_curpen <=1100){
		/* map from user 1000-1100 to device 0-255 */
		dc_mapcurpen = (int)( (float)(dc_dim2)*(dc_curpen-1000)/100.0);
		if(dc_mapcurpen>dc_dim2)dc_mapcurpen=dc_dim2;
		/* non linear map since human eye beholds intensity rather than
		 dots, e.g., density versus linear change in points plotted */
		dc_mapcurpen = (int) ( (float)dc_mapcurpen *
			(float)(dc_dim2+dc_dim2-dc_mapcurpen) / (float)dc_dim2 ) ;
	}
	if(dc_curpen >= 1000 && kolor == 0)return;
	if(ipen == 0 )
		(set_line_index)(CWHIT); 
        else if(ipen < 0){
                ipen = 1;
                savpen = ipen;
        } else if(ipen < 1000){
			ipen = ipen%(maxcolors-1);
			if(ipen == 0)ipen = (maxcolors -1);
                savpen = ipen;
		if(ipen == 0)
                	set_line_index(CWHIT);
		else if(ipen == 1)
                	set_line_index(CBLK);
		else if(ipen == 2)
                	set_line_index(CRED);
		else if(ipen == 3)
                	set_line_index(CGRN);
		else if(ipen == 4)
                	set_line_index(CBLU);
		else if(ipen == 5)
                	set_line_index(CORNG);
		else if(ipen == 6)
                	set_line_index(CBLGR);
		else if(ipen == 7)
                	set_line_index(CYEL);
        } else if(ipen >= 1000){
		if(ipen > 1100)ipen=1100;
                index = (int)( (Maxcol-3) * (float)((ipen-1000)/100.0)+0.5) +2;
                if(index >= Maxcol)index = Maxcol - 1;
                set_line_index(index);
        }

}

static void set_line_index(int index)
{
	int linecolor;
        if(index == 0)
                linecolor = index;
        else if(index == 1)
                linecolor = Maxcol -1;
        else
                linecolor = index -1;
	XviG_SetColor(linecolor);
	local_pen = linecolor;
}


void di_cross(int* x, int* y, char *s)
{
	long cursor;
	int chr;
	int xx, yy, buttonpress;
	int oldpen, savepen;
	int menu_res;

	oldpen = dc_curpen;
	savepen = 1;
	dv_pen((INT)savepen);
	XviG_SendMessage( 2, dc_ClipLeft, dc_ClipBottom, 
		dc_ClipRight, dc_ClipTop);
	XviG_Flush();
	if(user_curstype == 0)
		XviG_SetCursor(XviG_CURSOR_ARROW);
	else if(user_curstype == 1)
		XviG_SetCursor(XviG_CURSOR_XORARROW);
	else if(user_curstype == 2)
		XviG_SetCursor(XviG_CURSOR_XHAIR);
	else if(user_curstype == 3)
		XviG_SetCursor(XviG_CURSOR_PLUS);
	else if(user_curstype == 4)
		XviG_SetCursor(XviG_CURSOR_BOX);
	else if(user_curstype == 5)
		XviG_SetCursor(XviG_CURSOR_RUBBER);
	else if(user_curstype == 6)
		XviG_SetCursor(XviG_CURSOR_OFF);
	/* begin loop */
	while(1){
			chr = XviG_GetCursor(XviG_KEY_BUTTON, &xx, &yy);
		/* convert from absolute screen coordinates to
			user coordinates
		user               screen
		(0,0)              (border,Ldc_top-border)
		(dc_right,dc_top)  (border+dc_right,border+title)
		(xpos,ypos)        (xx,yy)
		*/
		xpos = xx - border;
		ypos = Ldc_top - border - yy;
		if(chr == XviG_BUTTON1)		/* left */
			chr = 1;
		else if(chr == XviG_BUTTON2)	/* center */
			chr = 3;
		else if(chr == XviG_BUTTON3)	/* right */
			chr = 2;
		if(ypos > dc_top){
			menu_res = check_menu(xpos,ypos);
			if(menu_res == ACTION_NEXT){
				continue;
			} else if(menu_res == ACTION_QUIT){
				terminate();
			}
		} else {
			goto enderase;
		}
	}
enderase:
	XviG_CloseCursor();
	XviG_SetCursor(XviG_CURSOR_ARROW);
	/* safety */
	if(xpos < dc_left)xpos = dc_left;
	if(xpos > dc_right)xpos = dc_right;
	if(ypos < dc_bottom)ypos = dc_bottom;
	if(ypos > dc_top)ypos = dc_top;
	*x=INVX(xpos);
	*y=INVY(ypos);

	if(*x < 0 )*x = 0;
	if(*y < 0 )*y = 0;
	if(!isprint(chr)){
		if(chr < 0 )
			chr= 040;
		else if(chr > 3)
			chr = 040;
	}
	s[0] = chr;
	s[1] = '\0';
	dv_pen((INT)oldpen);
}

struct menu { int xl; int xh ; char *str ; int action; int cev; } ;
static struct menu m[] = { 
	{  -1, -1, "  Next  \0" , ACTION_NEXT, 0},
	{  -1, -1, "  Quit  \0" , ACTION_QUIT, 0}
} ;

static void draw_border(int xl, int yl, int xh, int yh)
{
	int oldx_old, oldy_old, oldpen;
	int i,j;
	int ii,jj;

	oldx_old = dc_oldx;
	oldy_old = dc_oldy;
	oldpen = dc_curpen;
	if(kolor == 2){
		/* gray */
		Button_Color_Dark  = 1050;
		Button_Color       = 1040;
		Button_Color_Light = 1010;
		if(black){
			Button_Color_Fore  =    1;
		} else {
			Button_Color_Fore  =    0;
		}
	} else if(kolor == 1 || kolor == 0){
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
	for(i = 0; i< border; i++){
		if(i == 0 || i == 1)
			dv_pen((INT)Button_Color_Light);
		else if(i==(border-1) || i == (border-2))
			dv_pen((INT)Button_Color_Dark);
		else
			dv_pen((INT)Button_Color);
	(*xline)(Ldc_top -yl-i,xl+i,xh-i);
	(*xline)(Ldc_top -yh+i,xl+i,xh-i);
	(*yline)(xl+i,Ldc_top -yl-i,Ldc_top -yh+i);
	(*yline)(xh-i,Ldc_top -yl-i,Ldc_top -yh+i);
	}

	/* reset all position parameters */
	show_menu();
	dc_oldx = oldx_old;
	dc_oldy = oldy_old;
	dc_curpen = oldpen;
	dv_pen((INT)dc_curpen);

}

static void show_menu()
{
	int nm, i, dx, nc;
	int revvideo_old;
	int savepen;
	int j;
	INT oldcll, oldclr, oldclt, oldclb;
	/* put out the command options at top */
	nm = sizeof(m)/sizeof(struct menu);
	nc = strlen(m[0].str)*8;
	/* get absolute locations */
	dx = (dc_right - dc_left - nc)/(nm-1);
	dx = nc;
	dc_oldy =  Ldc_top - title -border - border + 2 ;
	savepen = dc_curpen;
	/* set background */
		dv_pen((INT)Button_Color);
/*
		for(j=border;j<title+border;j++)
			(*xline)(j,(int)dc_left+border, (int)dc_right+border);
*/
	fillquad( 	(int)dc_left+border , border, 
			(int)dc_right+border, border,
		 	(int)dc_right+border, title+border-1,  
			(int)dc_left+border , title+border-1);

	dv_pen((INT)Button_Color_Fore);
	revvideo_old = revvideo;
	for(i=0,j=0;i<nm;i++){
		m[i].xl = j*dx;
		m[i].xh = j*dx + nc -1;
		dc_oldx = m[i].xl ; 
		lastbuttonx = m[i].xh+8;
		j++;
		/* draw shadow box */
		dv_pen((INT)Button_Color_Dark);
/*
			(*xline)(border+title-1,
				m[i].xl  +border,m[i].xh+border);
			(*xline)(border+title-2,
				m[i].xl+1+border,m[i].xh+border);
			(*xline)(border+title-3,
				m[i].xl+1+border,m[i].xh+border);
*/
		fillquad( 	m[i].xl+border, title+border-1, 
				m[i].xh+border, title+border-1,
	 			m[i].xh+border, title+border-3,  
				m[i].xl+border, title+border-3);

/*
			(*yline)(m[i].xh  +border,
				border,border+title-2);
			(*yline)(m[i].xh-1+border,
				border,border+title-2);
			(*yline)(m[i].xh-2+border,
				border,border+title-2);
*/
		fillquad( 	m[i].xh  +border, border, 
				m[i].xh-2+border, border,
	 			m[i].xh-2+border, border+title-2,  
				m[i].xh  +border, border+title-2);


		dv_pen((INT)Button_Color_Light);
/*
			(*xline)(border  ,
				m[i].xl+border,m[i].xh-1+border);
			(*xline)(border+1,
				m[i].xl+border,m[i].xh-2+border);
*/
		fillquad( 	m[i].xl+  border, border, 
				m[i].xh-1+border, border,
	 			m[i].xh-2+border, border+1,  
				m[i].xl+  border, border+1);
/*
			(*yline)(m[i].xl  +border,
				border,border+title-2);
			(*yline)(m[i].xl+1+border,
				border,border+title-3);
			(*yline)(m[i].xl+2+border,
				border,border+title-4);
*/
		fillquad( 	m[i].xl  +border, border, 
				m[i].xl+2+border, border,
	 			m[i].xl+2+border, border+title-4,  
				m[i].xl  +border, border+title-2);
		revvideo = 0;
		dv_pen((INT)Button_Color_Fore);
		di_gottxt(m[i].str);
		oldcll = dc_ClipLeft;
		oldclr = dc_ClipRight;
		oldclt = dc_ClipTop;
		oldclb = dc_ClipBottom;
		dc_ClipLeft = Ldc_left;
		dc_ClipRight = Ldc_right;
		dc_ClipTop = Ldc_top;
		dc_ClipBottom = Ldc_bottom;
		dc_ClipLeft = oldcll;
		dc_ClipRight = oldclr;
		dc_ClipTop = oldclt;
		dc_ClipBottom = oldclb;
	}
	dv_pen((INT)savepen);
	/* put a box around */
	
	revvideo = revvideo_old;
}


static int check_menu(int xpos, int ypos)
{
	int nm, i, dx, nc;
	nm = sizeof(m)/sizeof(struct menu);
	nc = strlen(m[0].str)*8;
	if(ypos <= dc_top ) return(-1);
	for(i=0;i<nm;i++){
		if(xpos >= (m[i].xl) && xpos <= (m[i].xh)) {
			return(m[i].action);
		}
	}
	return(-1);
}

static void clearscr()
{
	int j;
	dv_pen((INT)0);
/* this works but is slow use the fill primitive
	for(j=border+title;j <= Ldc_top-border ; j++)
		(*xline)(j,Ldc_left+border,Ldc_right-border);
*/
	XviG_FillRectangle(Ldc_left+border,border+title,
		Ldc_right-border,Ldc_top - border);
}


char egaglyph[] = {
'\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000', /*   0 */
'\377','\377','\377','\377','\377','\377','\377','\377','\377','\377','\377','\377','\377','\377', /*   1 */
'\000','\000','\176','\377','\333','\377','\377','\303','\347','\377','\176','\000','\000','\000', /*   2 */
'\000','\000','\000','\154','\376','\376','\376','\376','\174','\070','\020','\000','\000','\000', /*   3 */
'\000','\000','\000','\020','\070','\174','\376','\174','\070','\020','\000','\000','\000','\000', /*   4 */
'\000','\000','\030','\074','\074','\347','\347','\347','\030','\030','\074','\000','\000','\000', /*   5 */
'\000','\000','\030','\074','\176','\377','\377','\176','\030','\030','\074','\000','\000','\000', /*   6 */
'\000','\000','\000','\000','\000','\030','\074','\074','\030','\000','\000','\000','\000','\000', /*   7 */
'\377','\377','\377','\377','\377','\347','\303','\303','\347','\377','\377','\377','\377','\377', /*   8 */
'\000','\000','\000','\000','\074','\146','\102','\102','\146','\074','\000','\000','\000','\000', /*   9 */
'\377','\377','\377','\377','\303','\231','\275','\275','\231','\303','\377','\377','\377','\377', /*  10 */
'\000','\000','\036','\016','\032','\062','\170','\314','\314','\314','\170','\000','\000','\000', /*  11 */
'\000','\000','\074','\146','\146','\146','\074','\030','\176','\030','\030','\000','\000','\000', /*  12 */
'\000','\000','\077','\063','\077','\060','\060','\060','\160','\360','\340','\000','\000','\000', /*  13 */
'\000','\000','\177','\143','\177','\143','\143','\143','\147','\347','\346','\300','\000','\000', /*  14 */
'\000','\000','\030','\030','\333','\074','\347','\074','\333','\030','\030','\000','\000','\000', /*  15 */
'\000','\000','\200','\300','\340','\370','\376','\370','\340','\300','\200','\000','\000','\000', /*  16 */
'\000','\000','\002','\006','\016','\076','\376','\076','\016','\006','\002','\000','\000','\000', /*  17 */
'\000','\000','\030','\074','\176','\030','\030','\030','\176','\074','\030','\000','\000','\000', /*  18 */
'\000','\000','\146','\146','\146','\146','\146','\146','\000','\146','\146','\000','\000','\000', /*  19 */
'\000','\000','\177','\333','\333','\333','\173','\033','\033','\033','\033','\000','\000','\000', /*  20 */
'\000','\174','\306','\140','\070','\154','\306','\306','\154','\070','\014','\306','\174','\000', /*  21 */
'\000','\000','\000','\000','\000','\000','\000','\000','\376','\376','\376','\000','\000','\000', /*  22 */
'\000','\000','\030','\074','\176','\030','\030','\030','\176','\074','\030','\176','\000','\000', /*  23 */
'\000','\000','\030','\074','\176','\030','\030','\030','\030','\030','\030','\000','\000','\000', /*  24 */
'\000','\000','\030','\030','\030','\030','\030','\030','\176','\074','\030','\000','\000','\000', /*  25 */
'\000','\000','\000','\000','\030','\014','\376','\014','\030','\000','\000','\000','\000','\000', /*  26 */
'\000','\000','\000','\000','\060','\140','\376','\140','\060','\000','\000','\000','\000','\000', /*  27 */
'\000','\000','\000','\000','\000','\300','\300','\300','\376','\000','\000','\000','\000','\000', /*  28 */
'\000','\000','\000','\000','\050','\154','\376','\154','\050','\000','\000','\000','\000','\000', /*  29 */
'\000','\000','\000','\020','\070','\070','\174','\174','\376','\376','\000','\000','\000','\000', /*  30 */
'\000','\000','\000','\376','\376','\174','\174','\070','\070','\020','\000','\000','\000','\000', /*  31 */
'\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000', /*  32 */
'\000','\000','\030','\074','\074','\074','\030','\030','\000','\030','\030','\000','\000','\000', /*  33 */
'\000','\146','\146','\146','\044','\000','\000','\000','\000','\000','\000','\000','\000','\000', /*  34 */
'\000','\000','\154','\154','\376','\154','\154','\154','\376','\154','\154','\000','\000','\000', /*  35 */
'\030','\030','\174','\306','\302','\300','\174','\006','\206','\306','\174','\030','\030','\000', /*  36 */
'\000','\000','\000','\000','\302','\306','\014','\030','\060','\146','\306','\000','\000','\000', /*  37 */
'\000','\000','\070','\154','\154','\070','\166','\334','\314','\314','\166','\000','\000','\000', /*  38 */
'\000','\060','\060','\060','\140','\000','\000','\000','\000','\000','\000','\000','\000','\000', /*  39 */
'\000','\000','\014','\030','\060','\060','\060','\060','\060','\030','\014','\000','\000','\000', /*  40 */
'\000','\000','\060','\030','\014','\014','\014','\014','\014','\030','\060','\000','\000','\000', /*  41 */
'\000','\000','\000','\000','\146','\074','\377','\074','\146','\000','\000','\000','\000','\000', /*  42 */
'\000','\000','\000','\000','\030','\030','\176','\030','\030','\000','\000','\000','\000','\000', /*  43 */
'\000','\000','\000','\000','\000','\000','\000','\000','\030','\030','\030','\060','\000','\000', /*  44 */
'\000','\000','\000','\000','\000','\000','\376','\000','\000','\000','\000','\000','\000','\000', /*  45 */
'\000','\000','\000','\000','\000','\000','\000','\000','\000','\030','\030','\000','\000','\000', /*  46 */
'\000','\000','\002','\006','\014','\030','\060','\140','\300','\200','\000','\000','\000','\000', /*  47 */
'\000','\000','\174','\306','\316','\336','\366','\346','\306','\306','\174','\000','\000','\000', /*  48 */
'\000','\000','\030','\070','\170','\030','\030','\030','\030','\030','\176','\000','\000','\000', /*  49 */
'\000','\000','\174','\306','\006','\014','\030','\060','\140','\306','\376','\000','\000','\000', /*  50 */
'\000','\000','\174','\306','\006','\006','\074','\006','\006','\306','\174','\000','\000','\000', /*  51 */
'\000','\000','\014','\034','\074','\154','\314','\376','\014','\014','\036','\000','\000','\000', /*  52 */
'\000','\000','\376','\300','\300','\300','\374','\006','\006','\306','\174','\000','\000','\000', /*  53 */
'\000','\000','\070','\140','\300','\300','\374','\306','\306','\306','\174','\000','\000','\000', /*  54 */
'\000','\000','\376','\306','\006','\014','\030','\060','\060','\060','\060','\000','\000','\000', /*  55 */
'\000','\000','\174','\306','\306','\306','\174','\306','\306','\306','\174','\000','\000','\000', /*  56 */
'\000','\000','\174','\306','\306','\306','\176','\006','\006','\014','\170','\000','\000','\000', /*  57 */
'\000','\000','\000','\030','\030','\000','\000','\000','\030','\030','\000','\000','\000','\000', /*  58 */
'\000','\000','\000','\030','\030','\000','\000','\000','\030','\030','\060','\000','\000','\000', /*  59 */
'\000','\000','\006','\014','\030','\060','\140','\060','\030','\014','\006','\000','\000','\000', /*  60 */
'\000','\000','\000','\000','\000','\176','\000','\000','\176','\000','\000','\000','\000','\000', /*  61 */
'\000','\000','\140','\060','\030','\014','\006','\014','\030','\060','\140','\000','\000','\000', /*  62 */
'\000','\000','\174','\306','\306','\014','\030','\030','\000','\030','\030','\000','\000','\000', /*  63 */
'\000','\000','\174','\306','\306','\336','\336','\336','\334','\300','\174','\000','\000','\000', /*  64 */
'\000','\000','\020','\070','\154','\306','\306','\376','\306','\306','\306','\000','\000','\000', /*  65 */
'\000','\000','\374','\146','\146','\146','\174','\146','\146','\146','\374','\000','\000','\000', /*  66 */
'\000','\000','\074','\146','\302','\300','\300','\300','\302','\146','\074','\000','\000','\000', /*  67 */
'\000','\000','\370','\154','\146','\146','\146','\146','\146','\154','\370','\000','\000','\000', /*  68 */
'\000','\000','\376','\146','\142','\150','\170','\150','\142','\146','\376','\000','\000','\000', /*  69 */
'\000','\000','\376','\146','\142','\150','\170','\150','\140','\140','\360','\000','\000','\000', /*  70 */
'\000','\000','\074','\146','\302','\300','\300','\336','\306','\146','\072','\000','\000','\000', /*  71 */
'\000','\000','\306','\306','\306','\306','\376','\306','\306','\306','\306','\000','\000','\000', /*  72 */
'\000','\000','\074','\030','\030','\030','\030','\030','\030','\030','\074','\000','\000','\000', /*  73 */
'\000','\000','\036','\014','\014','\014','\014','\014','\314','\314','\170','\000','\000','\000', /*  74 */
'\000','\000','\346','\146','\154','\154','\170','\154','\154','\146','\346','\000','\000','\000', /*  75 */
'\000','\000','\360','\140','\140','\140','\140','\140','\142','\146','\376','\000','\000','\000', /*  76 */
'\000','\000','\306','\356','\376','\376','\326','\306','\306','\306','\306','\000','\000','\000', /*  77 */
'\000','\000','\306','\346','\366','\376','\336','\316','\306','\306','\306','\000','\000','\000', /*  78 */
'\000','\000','\070','\154','\306','\306','\306','\306','\306','\154','\070','\000','\000','\000', /*  79 */
'\000','\000','\374','\146','\146','\146','\174','\140','\140','\140','\360','\000','\000','\000', /*  80 */
'\000','\000','\174','\306','\306','\306','\306','\326','\336','\174','\014','\016','\000','\000', /*  81 */
'\000','\000','\374','\146','\146','\146','\174','\154','\146','\146','\346','\000','\000','\000', /*  82 */
'\000','\000','\174','\306','\306','\140','\070','\014','\306','\306','\174','\000','\000','\000', /*  83 */
'\000','\000','\176','\176','\132','\030','\030','\030','\030','\030','\074','\000','\000','\000', /*  84 */
'\000','\000','\306','\306','\306','\306','\306','\306','\306','\306','\174','\000','\000','\000', /*  85 */
'\000','\000','\306','\306','\306','\306','\306','\306','\154','\070','\020','\000','\000','\000', /*  86 */
'\000','\000','\306','\306','\306','\306','\326','\326','\376','\174','\154','\000','\000','\000', /*  87 */
'\000','\000','\306','\306','\154','\070','\070','\070','\154','\306','\306','\000','\000','\000', /*  88 */
'\000','\000','\146','\146','\146','\146','\074','\030','\030','\030','\074','\000','\000','\000', /*  89 */
'\000','\000','\376','\306','\214','\030','\060','\140','\302','\306','\376','\000','\000','\000', /*  90 */
'\000','\000','\074','\060','\060','\060','\060','\060','\060','\060','\074','\000','\000','\000', /*  91 */
'\000','\000','\200','\300','\340','\160','\070','\034','\016','\006','\002','\000','\000','\000', /*  92 */
'\000','\000','\074','\014','\014','\014','\014','\014','\014','\014','\074','\000','\000','\000', /*  93 */
'\020','\070','\154','\306','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000', /*  94 */
'\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\377','\000', /*  95 */
'\060','\060','\030','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000', /*  96 */
'\000','\000','\000','\000','\000','\170','\014','\174','\314','\314','\166','\000','\000','\000', /*  97 */
'\000','\000','\340','\140','\140','\170','\154','\146','\146','\146','\174','\000','\000','\000', /*  98 */
'\000','\000','\000','\000','\000','\174','\306','\300','\300','\306','\174','\000','\000','\000', /*  99 */
'\000','\000','\034','\014','\014','\074','\154','\314','\314','\314','\166','\000','\000','\000', /* 100 */
'\000','\000','\000','\000','\000','\174','\306','\376','\300','\306','\174','\000','\000','\000', /* 101 */
'\000','\000','\070','\154','\144','\140','\360','\140','\140','\140','\360','\000','\000','\000', /* 102 */
'\000','\000','\000','\000','\000','\166','\314','\314','\314','\174','\014','\314','\170','\000', /* 103 */
'\000','\000','\340','\140','\140','\154','\166','\146','\146','\146','\346','\000','\000','\000', /* 104 */
'\000','\000','\030','\030','\000','\070','\030','\030','\030','\030','\074','\000','\000','\000', /* 105 */
'\000','\000','\006','\006','\000','\016','\006','\006','\006','\006','\146','\146','\074','\000', /* 106 */
'\000','\000','\340','\140','\140','\146','\154','\170','\154','\146','\346','\000','\000','\000', /* 107 */
'\000','\000','\070','\030','\030','\030','\030','\030','\030','\030','\074','\000','\000','\000', /* 108 */
'\000','\000','\000','\000','\000','\354','\376','\326','\326','\326','\306','\000','\000','\000', /* 109 */
'\000','\000','\000','\000','\000','\334','\146','\146','\146','\146','\146','\000','\000','\000', /* 110 */
'\000','\000','\000','\000','\000','\174','\306','\306','\306','\306','\174','\000','\000','\000', /* 111 */
'\000','\000','\000','\000','\000','\334','\146','\146','\146','\174','\140','\140','\360','\000', /* 112 */
'\000','\000','\000','\000','\000','\166','\314','\314','\314','\174','\014','\014','\036','\000', /* 113 */
'\000','\000','\000','\000','\000','\334','\166','\146','\140','\140','\360','\000','\000','\000', /* 114 */
'\000','\000','\000','\000','\000','\174','\306','\160','\034','\306','\174','\000','\000','\000', /* 115 */
'\000','\000','\020','\060','\060','\374','\060','\060','\060','\066','\034','\000','\000','\000', /* 116 */
'\000','\000','\000','\000','\000','\314','\314','\314','\314','\314','\166','\000','\000','\000', /* 117 */
'\000','\000','\000','\000','\000','\146','\146','\146','\146','\074','\030','\000','\000','\000', /* 118 */
'\000','\000','\000','\000','\000','\306','\306','\326','\326','\376','\154','\000','\000','\000', /* 119 */
'\000','\000','\000','\000','\000','\306','\154','\070','\070','\154','\306','\000','\000','\000', /* 120 */
'\000','\000','\000','\000','\000','\306','\306','\306','\306','\176','\006','\014','\370','\000', /* 121 */
'\000','\000','\000','\000','\000','\376','\314','\030','\060','\146','\376','\000','\000','\000', /* 122 */
'\000','\000','\016','\030','\030','\030','\160','\030','\030','\030','\016','\000','\000','\000', /* 123 */
'\000','\000','\030','\030','\030','\030','\000','\030','\030','\030','\030','\000','\000','\000', /* 124 */
'\000','\000','\160','\030','\030','\030','\016','\030','\030','\030','\160','\000','\000','\000', /* 125 */
'\000','\000','\166','\334','\000','\000','\000','\000','\000','\000','\000','\000','\000','\000', /* 126 */
'\000','\000','\000','\000','\020','\070','\154','\306','\306','\376','\000','\000','\000','\000', /* 127 */
} ;


static void gottxtp(char *s)	/* actual output for PC */
{
	int i,j,k,kmax;
	char *glyph;
	unsigned int mask, l;
	int m,x,y;
	kmax = 14;
	glyph = &egaglyph[0];
	while((i = *s) != '\0'){
		j = i*14;
		y = dc_oldy+kmax-1;
		for(k=0;k<kmax;k++){
			l = glyph[j+k];
			mask = 0200;
			x = dc_oldx;
			for (m=0;m<8;m++){
				if( (l&mask)!=0 && x >= 0 && 
					x <= dc_right && y >=border 
					&& y <=Ldc_top)
					dv_zpoint(x,y);
				mask = (mask>>1);
				x++;
			}
			y--;
		}
		dc_oldx+=8;
		s++;
	}
}

/* output text string to graphics terminal at current position */
void di_gottxt( char *s)
{
	int oldpen;
	int i,j;
	char  ch[2];
	oldpen = savpen;
	ch[1] = '\000';
	/* strip off blanks at the end of the string s */
	j = strlen(s);
	if(j >= NAMSTR)j=NAMSTR - 1;
	strncpy(name,s,j);
	/* drop trailing blanks */
	if(revvideo == 0){
		for(i=j-1;i>=0;i--)
			if(name[i] != ' ' && name[i] != '\0')break;
		name[++i] = '\0';
	} else
		i = j;
	for(j=0;j<i;j++){
		if(revvideo )cblank();
		ch[0]=name[j];
		if(revvideo == 2)
			XviG_SetColor(0);
		gottxtp(ch);
	}
	/* reset reversed video to normal */
	if(revvideo == 2)
		dv_pen((INT)oldpen);
	revvideo = 0;
}

static void cblank()
{
	int oldpen;
	char ch[2];
	/* use special character set to output a blank */
	/* and then return cursor to original position */
	oldpen=savpen;
	if(revvideo == 2)
		XviG_SetColor(Maxcol -1);
	else
		XviG_SetColor(0);
	ch[0] = '\001';
	ch[1] = '\000';
	gottxtp(ch);
	dv_pen((INT)oldpen);
	dc_oldx-=8;
}
void dv_movcur (INT isx,INT isy)
{
        mx = isx;
        my = isy;
}

void dv_concur(INT isx,INT isy)
{
        XviG_DrawLine(mx+border,Ldc_top-my-border,isx+border,Ldc_top-isy-border);
        mx = isx;
        my = isy;
}
static void xvigxline(int y,int xs,int xe)
{
	XviG_DrawLine(xs, y, xe, y);
}

static void xvigyline(int x,int ys,int ye)
{
	XviG_DrawLine(x, ys, x, ye);
}

static void xvigpoint(int x,int y)
{
	XviG_DrawPoint(x, y);
}


void dv_mesg(char *mesg)
{
	int ls;
	int revvideo_old;
	int xoring_old;
	int curpen_old;
	ls = strlen(mesg);
	dc_oldy = Ldc_top - title -border - border + 2 ;
	dc_oldx = lastbuttonx;
	revvideo_old = revvideo;
	xoring_old = xoring;
	if(xoring == ON)XviG_SetGC(0);
	xoring = OFF;
	revvideo = 0;
	curpen_old = dc_curpen;
	/* clear out the old mesg */
	if(strlen(omesg)){
		dv_pen((INT)Button_Color);
		di_gottxt(omesg);
	}
	/* output string but do not exceed the limits of the window */
	dc_oldx = lastbuttonx;
	dv_pen((INT)Button_Color_Fore);
	di_gottxt(mesg);
	revvideo = revvideo_old;
	xoring = xoring_old;
	if(xoring == ON)XviG_SetGC(1);
	dc_curpen = curpen_old;
	dv_pen((INT)dc_curpen);
	strcpy(omesg,mesg);
	XviG_Flush();
}


void dv_fillp(INT n, INT *x,INT *y)
{
	int i,j;
	/*RBH_SetDither(local_pen);*/
	if((coords = (int *)calloc(n+n, sizeof(int))) == NULL)
		return;
	for(i=0,j=0;i<n; i++){
		coords[j++] = x[i]+border;
		coords[j++] = Ldc_top - y[i] - border;
	}
	XviG_FillPolygon(coords, n);
	free(coords);
}


static void terminate()
{
	XviG_Exit();
	exit(0);
}


static void fillquad(int x0,int z0,int x1,int z1,
	int x2, int z2, int x3, int z3)
{
	/*RBH_SetDither(local_pen);*/
	/* allocate space for 4 corners */
	if((coords = (int *)calloc(8, sizeof(int))) == NULL)
		return;
	coords[0] = x0 ;
	coords[1] = z0 ;
	coords[2] = x1 ;
	coords[3] = z1 ;
	coords[4] = x2 ;
	coords[5] = z2 ;
	coords[6] = x3 ;
	coords[7] = z3 ;
	XviG_FillPolygon(coords, 4);
	free(coords);
}

/* $XConsortium: ParseGeom.c,v 11.19 94/04/17 20:20:23 rws Exp $ */

/*

Copyright (c) 1985, 1986, 1987  X Consortium

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE X CONSORTIUM BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

Except as contained in this notice, the name of the X Consortium shall
not be used in advertising or otherwise to promote the sale, use or
other dealings in this Software without prior written authorization
from the X Consortium.

*/


/*
 *    XParseGeometry parses strings of the form
 *   "=<width>x<height>{+-}<xoffset>{+-}<yoffset>", where
 *   width, height, xoffset, and yoffset are unsigned integers.
 *   Example:  "=80x24+300-49"
 *   The equal sign is optional.
 *   It returns a bitmask that indicates which of the four values
 *   were actually found in the string.  For each value found,
 *   the corresponding argument is updated;  for each value
 *   not found, the corresponding argument is left unchanged. 
 */

static void do_parsegeom(char *str)
{
	/* modified from Xconsortium ParseGeom.c */
	/* parse the X11 geometry string
		=<width>x<height>{+-}<xoffset>{+-}<yoffset>
	*/
	char *strind;
	int width, height, xoff, yoff;
	int tempWidth, tempHeight, tempX, tempY;
	char *nextCharacter;
	tempWidth = -1;
	tempHeight = -1;
	tempX = 9999;
	tempY = 9999;

	
	if( (str == NULL) || (*str == '\0'))
		return;
	if(*str == '=')
		str++;	/* ignore leading = */
	strind = (char *)str;
	if(*strind == '+' || *strind == '-' || *strind == 'x')
			return;
	/* now try to get geometry */
	
	/* get the width */
	tempWidth = ReadInteger(strind, &nextCharacter);
	if (strind == nextCharacter) return ;
	strind = nextCharacter;

	if (*strind == 'x' || *strind == 'X') {	
		strind++;
		tempHeight = ReadInteger(strind, &nextCharacter);
		if (strind == nextCharacter) return ;
		strind = nextCharacter;
	}

	if ((*strind == '+') || (*strind == '-')) {
		if (*strind == '-') {
  			strind++;
			tempX = -ReadInteger(strind, &nextCharacter);
			if (strind == nextCharacter)
			    goto CONT ;
			strind = nextCharacter;

		}
		else
		{	strind++;
			tempX = ReadInteger(strind, &nextCharacter);
			if (strind == nextCharacter)
			    goto CONT;
			strind = nextCharacter;
		}
		if ((*strind == '+') || (*strind == '-')) {
			if (*strind == '-') {
				strind++;
				tempY = -ReadInteger(strind, &nextCharacter);
				if (strind == nextCharacter)
			    	    goto CONT;
				strind = nextCharacter;

			}
			else
			{
				strind++;
				tempY = ReadInteger(strind, &nextCharacter);
				if (strind == nextCharacter)
			    	    goto CONT;
				strind = nextCharacter;
			}
		}
	}
CONT:
	if(tempWidth > 0 && tempHeight > 0){
		Width = tempWidth ;
		Height = tempHeight;
	}
	if(tempX != 9999 && tempY != 9999){
		Xpos = tempX;
		Ypos = tempY;
	}
	
}

static int ReadInteger(char *string, char **NextString)
{
    register int Result = 0;
    int Sign = 1;
    
    if (*string == '+')
	string++;
    else if (*string == '-')
    {
	string++;
	Sign = -1;
    }
    for (; (*string >= '0') && (*string <= '9'); string++)
    {
	Result = (Result * 10) + (*string - '0');
    }
    *NextString = string;
    if (Sign >= 0)
	return (Result);
    else
	return (-Result);
}

static void usage(void)
{
fprintf(stderr,"plotxvig [options]\n");
fprintf(stderr,"-V                         Program Version \n" );
fprintf(stderr,"-Sscalefac  (default=1.0)  Plot magnifier  \n" );
fprintf(stderr,"-R          (default off)  Rotate plot 90 degrees \n" );
fprintf(stderr,"-N          (default off)  Turn off shading \n" );
fprintf(stderr,"-I         (default off)   Invert background (e.g., black)\n");
fprintf(stderr,"-Ffont      (default 0)    Default font \n" );
fprintf(stderr,"                             0 Roman \n" );
fprintf(stderr,"                             1 Roman \n" );
fprintf(stderr,"                             2 Italic \n" );
fprintf(stderr,"                             3 Bold \n" );
fprintf(stderr,"                             4 Symbol (Greek) \n" );
fprintf(stderr,"-K          (default color) Color output Red->Green->Blue \n" );
fprintf(stderr,"-KR         (default -K   ) Color output Red->White->Blue \n" );
fprintf(stderr,"-KB         (default -K   ) Color output Blue->White->Red \n" );
fprintf(stderr,"-G          (default -K   ) Gray output \n" );
fprintf(stderr,"-W          (default 0)    Line width in units of 0.001 in or 0.0025cm) \n" );
fprintf(stderr,"-geometry   (default off)  Parse X11 geometry string\n");
fprintf(stderr,"       -geometry WIDTHxHEIGHT{+-}XOFF{+-}YOFF\n");
fprintf(stderr,"-p                         Put up positioning grid every 1000 CALPLOT units, same as -p1\n" );
fprintf(stderr,"-p2 every 500, -p4 every 250, -p10 every 100\n");
fprintf(stderr,"-h                         Do not execute, show options \n" );
fprintf(stderr,"-?                         Do not execute, show options \n" );
	exit (0);
}
void dv_clip()
{
}

static void draw_grid()
{
	int i,j;
	/* new 01/12/2001 draw reference grid */
	INT DX, UX;
	INT DY, UY;
	DX = (INT)1000/dogrid_inc;
	DY = (INT)1000/dogrid_inc;
	UX = (INT)(10000.0/dc_scalefac);
	UY = (INT)( 8000.0/dc_scalefac);
	dv_pen((INT)1);
	if(dogrid){
		for(i=0 ; i <= UX ; i+=DX){
			if(dc_rotate){
				draw_dash(i, i, 0, i, UY);
			} else {
				draw_dash(i, i, 0, i, UY);
			}
		}
		for(j=0 ; j <= UY ; j+=DY){
			if(dc_rotate){
				draw_dash(j, 0, j, UX, j);
			} else {
				draw_dash(j, 0, j, UX, j);
			}
		}
	}
}

static void draw_dash(int k, INT ix0, INT iy0,INT ix1, INT iy1)
{
	INT i;
	int din, di2;
	INT jx0,jx1,jy0,jy1;
	/* first order everything */
	if(iy0 > iy1){
		i = iy0; iy0 = iy1; iy1 = i;
	}
	if(ix0 > ix1){
		i = ix0; ix0 = ix1; ix1 = i;
	}
	if(dogrid_inc <=2)
		din = 100;
	else if(dogrid_inc == 4)
		din = 50;
	else if(dogrid_inc == 5)
		din = 40;
	else if(dogrid_inc == 10)
		din = 20;
	di2 = 2*din;
	if( (k % 1000) ==  0){
		/* no dashes */
		if(dc_rotate){
			dv_movcur(IY(ix0), IX(iy0));
			dv_concur(IY(ix1), IX(iy1));
		} else {
			dv_movcur(IX(ix0), IY(iy0));
			dv_concur(IX(ix1), IY(iy1));
		}
	} else {
		if(ix0 == ix1){
			/* dash in y direction */
			for(i=iy0 ; i < iy1 ; i+=di2){
				if(dc_rotate){
					jx0 = IY(ix0);
					jx1 = IY(ix1);
					jy0 = IX(i    );
					jy1 = IX(i+din);
				} else {
					jx0 = IX(ix0);
					jx1 = IX(ix1);
					jy0 = IY(i    );
					jy1 = IY(i+din);
				}
				dv_movcur(jx0,jy0);dv_concur(jx1,jy1);
			}
		} else if(iy0 = iy1){
			/* dash in x direction */
			for(i=ix0 ; i < ix1 ; i+=di2){
				if(dc_rotate){
					jx0 = IY(i    );
					jx1 = IY(i+din);
					jy0 = IX(iy0);
					jy1 = IX(iy1);
				} else {
					jx0 = IX(i    );
					jx1 = IX(i+din);
					jy0 = IY(iy0);
					jy1 = IY(iy1);
				}
				dv_movcur(jx0,jy0);dv_concur(jx1,jy1);
			}
		} else {
			/* safety diagonal line */
			if(dc_rotate){
				dv_movcur(IY(ix0), IX(iy0));
				dv_concur(IY(ix1), IX(iy1));
			} else {
				dv_movcur(IX(ix0), IY(iy0));
				dv_concur(IX(ix1), IY(iy1));
			}
		}
	}
}
