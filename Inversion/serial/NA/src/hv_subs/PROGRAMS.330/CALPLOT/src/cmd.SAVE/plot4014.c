/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLOT4014                                              c
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
	02 APR 2004 - inplemented gend(mode) at higher level which
		required change in dv_closepl
*/
/* device dependent plot program interface 			*/
/* Tektronix 4014 with extensions for Kermit 3.1 and TeraTerm   */
/* It is assumed that the plotgen program clips and requires	*/
/*	the very low level primitives				*/
/*		dv_clip						*/
/*		dv_closepl					*/
/*		dv_concur					*/
/*		dv_erase					*/
/*		dv_fillp					*/
/*		dv_font						*/
/*		di_gsymb					*/
/*		di_gcmdln					*/
/*		dv_movcur					*/
/*		onintr						*/
/*		dv_openpl					*/
/*		dv_pen						*/
/*		zpoint						*/

#include	<stdio.h>
#include	<ctype.h>
#include	<math.h>
#include	<string.h>
#include	<stdlib.h>
#include	<stddef.h>
#ifdef MSDOS
#include	<fcntl.h>
#include	<io.h>
#define INT 	long
#else
#define	INT	int
#endif



#define INVX(n) ( (n - 0.5)/(dc_xscale*dc_scalefac) - dc_xoff)
#define INVY(n) ( (n - 0.5)/(dc_yscale*dc_scalefac) - dc_yoff)




/* Global Function Prototypes	*/
void dv_clip();
void dv_closepl(int mode);
void di_cross(int *x, int *y, char *c);
void dv_concur(INT isx,INT isy);
void dv_fillp(INT n, INT *x,INT *y);
void di_gintxt(int cnt,char *s);
void dv_openpl(int showmenu);
void dv_movcur(INT isx,INT isy);
void dv_pen(INT p);
void dv_zpoint(INT x,INT y);

/* Global Variables		*/
extern INT	dc_oldx	;	/* previous value of x coordinate */
extern INT	dc_oldy	;	/* previous value of y coordinate */
INT 		dc_ClipRight, dc_ClipTop, 
		dc_ClipBottom, dc_ClipLeft; /* current clip region*/
INT dc_ClipRight_sv, dc_ClipTop_sv, dc_ClipBottom_sv, dc_ClipLeft_sv;
INT	dc_iseps = 0;	/* invoke special scaling for EPS */
INT		dc_xlinewidth = 0;
INT		dc_ylinewidth = 0;
INT		dc_linewidth = 0;
INT		dc_minlinewidth = 0;
INT		dc_herelinewidth = 0;	/*set linewidth here in this routine */
INT		dc_newlinewidth = 1;	/*flag to indicate a line width change*/
INT		dc_hardwarefill = 2; 	/* flag to have all filling 								done by hardware */
INT		dc_curfont = 0;
double		dc_scalefac = 1.0;
int		dc_hasmouse = 0;
	/* flags and values set by command line argument		*/
static INT	scaledplot = 0;
INT		dc_rotate = 0;
extern int	dc_sleeptime ;
INT		dc_shadeoff = 0 ;	/* permit shading unless turned 
					off in command line */
int		dc_curpen = 1 ;	/* current pen value */
extern INT	dc_dim2;	/* dimensions of dither matrix dv_zzpoint */
int		dc_mapcurpen;	/* mapping into dither for dv_zzpoint */
INT		dc_color = 0 ;	/* no hardware color shading */
INT		dc_ColorInfo = 0; /* no color shading */
extern double	dc_xscale ;
extern double	dc_yscale ;
extern double	dc_xoff   ;
extern double	dc_yoff   ;
/* device clipping region */
INT 		dc_left = 0, 	/* device limits for plot window */
		dc_right=4095, 
		dc_top=3119,  
		dc_bottom = 0;
double dc_Nxpi	= 409.5;
double dc_Nypi	= 389.9;
INT dc_oldx1	= -1;
INT dc_oldy1 	= -1;


/* Local  Function Prototypes	*/
static INT maxval(INT x,INT y);
static INT minval(INT x,INT y);
static void do_plaid(INT xmn,INT ymn,INT xmx,INT ymx,INT patx,INT paty,INT lenx,INT leny);
static void symps(INT x,INT y,INT ht,char *s,INT ang,INT nchar);
static float value(float n1, float n2, float hue);
static void eps_rotate_box(INT x,INT y,INT ht,INT angle,INT nchar);
static void eps_bounding_box( INT x, INT y);
static void eps_box_rotate(INT *tx,INT *ty,float ct,float st,INT x0,INT y0,INT x,INT y);
static void show_clip(INT lx,INT ly,INT ux,INT uy);
static void psinit(void );
static void hls_to_rgb(float *r,float *g,float *b,float h,float l,float s);
static void coord(float *red,float *green,float *blue,INT index);
static void setcolor();
static void con (INT isx,INT isy);
static void mov (INT isx,INT isy);
static void endpage();
static void newpage();
static void gvect(int ix,  int iy);
static void outcolor(int c);
static char gin(int *ix,int *iy);
static void galpha();
static void ggraph();
static void usage(void);

/* Local  Variables		*/

#define NFONTS 4
static int reduc = 1;
static int oi1 = -1;
static int oi2 = -1;
static int oi3 = -1;
static int oi4 = -1;
static int oi5 = -1;

static INT	defaultfont = 0;
static INT revvideo = 0;
static	int is4025 = 0;		/* 0 vanilla 4014 */
				/* 1 Tektronix 4025	*/
				/* 2 MS-DOS Kermit 3.00 color support */
static INT	black = 0;	/* paint onto white paper */
static INT	white = 1;
#define BLACK	0
#define RED	1
#define GREEN	2
#define YELLOW	3
#define BLUE	4
#define PURPLE	5 
#define CYAN	6
#define WHITE	7
static int  colormap[] = {
	BLACK,
	WHITE,
	RED,
	GREEN,
	BLUE,
	PURPLE,
	CYAN,
	YELLOW
	} ;






void onintr(int arg)
{
	if(is4025==1) {
		fputs("!MON H !COM 31\n", stdout);
	}
	exit(0);
}

void dv_closepl(int mode)
{
	fflush(stdout);
	fputc(31,stdout);
	if(is4025==1) {
		/* fputs("!WORK 0 !MON H !COM 31\n", stdout); */
		fputs("!MON H !COM 31\n", stdout);
	}
}

void dv_openpl(int showmenu)
{
	if(is4025==1) {
		fputs("\037WORK 33 \037GRAPHIC 1,33 \037SHR \n", stdout);
		fputs("\037COM 33\n!ERA G!WORK H\037\035\037\n", stdout);
	} else if (is4025==2) {
		printf("\033\014");
	} else if (is4025==3) {
		printf("\033\014");
	} else  {
		printf("\033\014");
	}
	dc_ClipLeft = dc_left;
	dc_ClipRight = dc_right;
	dc_ClipTop = dc_top;
	dc_ClipBottom = dc_bottom;
	dc_ClipLeft_sv = dc_left;
	dc_ClipRight_sv = dc_right;
	dc_ClipTop_sv = dc_top;
	dc_ClipBottom_sv = dc_bottom;
}

void dv_erase(INT mode)
{
	fflush(stdout);
	if(dc_sleeptime)
		sleep(dc_sleeptime);
	if(is4025==1) {
		fputs("\037\035\037!ERA G\n", stdout);
	} else {
		fputc(27,stdout);
		fputc(12,stdout);
	}
	dc_oldx1 = -1;
	dc_oldy1 = -1;
	oi1 = -1 ; oi2 = -1; oi3 = -1; oi4 = -1; oi5 = -1;
}

void dv_zpoint(int x,int y)
{
	dv_movcur(x,y);
	dv_concur(x,y);
}

void dv_movcur (INT isx,INT isy)
{
	fputc(29,stdout);
	gvect(reduc*isx,reduc*isy);
}

void dv_concur(INT isx,INT isy)
{
	if(dc_curpen)
		gvect(reduc*isx,reduc*isy);
}

void dv_pen(INT jpen)
{
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
	if(is4025==2){
		ipen = (( ipen-1)%7 + 1 );
		outcolor(colormap[ipen]);
	} else if (is4025==3) {
		fputc('\033',stdout);
		fputc('M',stdout);
		fputc('L',stdout);
		if(ipen > 0)
			ipen = (ipen-1)%15 + 1 ;
		fputc(ipen+48  ,stdout);
	}

}

static void outcolor(int c)
{
	printf("\033[0;3%1dm",c);
}

static void gvect(int ix,int iy)
{
	int i1,i2,i3,i4,i5;
	i1 = 32 + (iy/128)%32;		/* high Y */
	i2 = 96 + (4*iy)%13 + ix%4;	/* LSBXY */
	i3 = 96 + (iy/4)%32;		/* Low Y */
	i4 = 32 + (ix/128)%32;		/* high X */
	i5 = 64 + (ix/4)%32;		/* Low X */
		if(i1 != oi1){
			fputc(i1,stdout);
		}
		if(i2 != oi2){
			fputc(i2,stdout);
			fputc(i3,stdout);
		} else if(i4 != oi4 || i3 != oi3) {
			fputc(i3,stdout);
		}
		if (i4 != oi4){
			fputc(i4,stdout);
		}
		fputc(i5,stdout);
		oi1 = i1;
		oi2 = i2;
		oi3 = i3;
		oi4 = i4;
		oi5 = i5;
}

void di_gcmdln(int argc,char *argv[])
{
	char *cp;
	double atof();
	while(argc-- > 1 ) {
		if(*argv[1] == '-'){
			switch(argv[1][1]){
			case 'V':
				fprintf(stderr,"CALPLOT (3.0) COPYRIGHT (C) 1997 Saint Louis University\n");
				break;
			case 'S':
				scaledplot = 1;
				cp = argv[1];
				cp++;
				cp++;
				dc_scalefac = atof(cp);
				if(dc_scalefac <= 0.0)
					dc_scalefac = 1.0;
				break;
			case 'R':
				dc_rotate = 1;
				break;
			case 'N':
				dc_shadeoff = 1;
				break;
			case 'D':
				cp = argv[1];
				cp++;
				cp++;
				reduc = atoi(cp);
				if(reduc < 1)reduc = 1;
				break;
			case 'F':	/* default font value */
				cp = argv[1];
				cp++;
				cp++;
				defaultfont = atoi(cp);
				if(defaultfont < 0)defaultfont = 0;
				if(defaultfont>0)
					defaultfont = (defaultfont-1)%4 + 1;
				dc_curfont = defaultfont;
				break;
			case 'W':
				cp = argv[1];
				cp++;
				cp++;
				dc_sleeptime = atoi(cp);
				break;
			case 'T':
				cp = &argv[1][1];
				/* TEKTRONIX 4025 */
				if(strncmp(cp,"T4025",5)==0)
					is4025 = 1;
				/* TERATERM V12  */
				if(strncmp(cp,"TT",2)==0)
					is4025 = 3;
				break;
			case 'K':
				is4025 = 2;
				/* MS-DOS KERMIT 3.00 color tek support */
				break;
			case '?':
			case 'h':
				usage();
				break;
			default:
				break;
			}
		}
		argv++;
	}
	/* some protection 	*/
	dc_scalefac = dc_scalefac/(double)reduc;
	dc_right = (int) ( (double)dc_right/(double)reduc );
	dc_top = (int) ( (double)dc_top/(double)reduc );
	dc_bottom = 0;
	dc_left = 0;
}


void dv_font(INT fontvalue)
{
	if(fontvalue > 0)
		dc_curfont = (fontvalue-1)%NFONTS + 1;
	else
		dc_curfont = defaultfont;
}


	
static INT maxval(INT x,INT y)
{
	if(x > y)
		return x;
	else
		return y;
}

static INT minval(INT x,INT y)
{
	if(x < y)
		return x;
	else
		return y;
}

void dv_cursor(INT curstyp)
{
}



void dv_fillp(INT n, INT x[],INT y[])
{
}
void di_gsymb(INT x,INT y,INT ht,INT ang,INT nchar,char *s)
{
	int des, lst ;
	int savpen;
	if(nchar > 0){
		if(revvideo){
			savpen = dc_curpen;
			lst = nchar * ht;
			des =  ht/2;
			if(revvideo == 1)
				dc_curpen = 0;
			else
				dc_curpen = savpen;
			di_fillr(x -10,y-des,x+lst,y+ht+des/2,0,0,0,0);
			if(revvideo == 1)
				dc_curpen = savpen;
			else
				dc_curpen = 0;
		}
		dv_symvec(x,y,ht,s,ang,nchar);
		if(revvideo){
			dc_curpen = savpen;
			revvideo = 0;
		}
	} else {
		dv_symvec(x,y,ht,s,ang,nchar);
	}
}



/* C- level routine to turn on cross hair and to receive characters
   The TEKTRONICS sends 4 bytes coordinate, 1 byte character and a
   CARRIAGE_RETURN. The carriage return is converted to a NEWLINE byu
   the UNIX terminal handler. The routine eats up characters from
   the terminal line buffer until it sees one terminated by a 
   NEWLINE. If exactly six characters are input, then we have
   a correct input sequence */
#include <termios.h>


void di_cross(int *x, int *y, char *c)
{
/* the lincoln lab version did not work so used Carl Johnson's gin */
	char gin();
	int ix,iy;
	*c=gin(&ix,&iy);
	*x=(int)INVX(ix);
	*y=(int)INVY(iy);
	c[1] = '\0';
}

static char gin(int *ix,int *iy)
{
	int ic;
	int i,cc[100];
	int nc;
	struct termios old,new;
	tcgetattr(0,&old);
	tcgetattr(0,&new);
	new.c_lflag &= ~ECHO;
	/* reset terminal line characteristics */
	tcsetattr(0,TCSANOW,&new);
	putchar('\027');
	putchar('\033');
	putchar('\032');
	fflush(stdout);
	nc = 0;
	for(i=0; i<100; i++) {
		cc[nc++] = (  getc(stdin) & 0177 );
		if(nc > 1 && ( cc[nc-1] == '\n' ) )
			break;
	}
	tcsetattr(0,TCSANOW,&old);
	cc[nc]='\0';
	printf("\n");
	if(nc != 6) {
		printf("nc = %d\n",nc);
		return('?');
	}
	*ix=4*(((cc[1]&037)<<5)+(cc[2]&037));
	*iy=4*(((cc[3]&037)<<5)+(cc[4]&037));
	return(cc[0]);
}

/* turn off vector plot mode and enable text mode */

static void galpha()
{
	putc(31,stdout);
}
/* enter into graphics mode */
static void ggraph()
{
	putc(29,stdout);
}

/* output text string to graphics terminal at current position */

#define NAMSTR  81
static char name[NAMSTR];

/* get text string from terminal */
void di_gintxt(int cnt,char *s)
{
	char ch[2];
	int last=0;
	strcpy(name,"");
	ch[0] = 0;
	ch[1] = 0;
	galpha();
	while(ch[0] != 13 && ch[0] != 10  && last < NAMSTR){/* keep doing this until we get  CR */
		ch[0] = getc(stdin);
		if(isprint(ch[0])){
			name[last] = ch[0];
			last++;
		}
	}
	if(last >= cnt)last = cnt-1 ;
	name[last] = '\0';
	strcpy(s,name);
}


void dv_mesg(char *mesg)
{
}

static void usage(void)
{
fprintf(stderr,"plot4014 [options]\n");
fprintf(stderr,"-V                         Program Version \n" );
fprintf(stderr,"-Sscalefac  (default=1.0)  Plot magnifier  \n" );
fprintf(stderr,"-R          (default off)  Rotate plot 90 degrees \n" );
fprintf(stderr,"-N          (default off)  Turn off shading \n" );
fprintf(stderr,"-Ddec       (default 1)    Decimate plot resolution  \n" );
fprintf(stderr,"            This reduces transmission time on modem, e.g., -D4\n");
fprintf(stderr,"-Ffont      (default 0)    Default font \n" );
fprintf(stderr,"                             0 Roman \n" );
fprintf(stderr,"                             1 Roman \n" );
fprintf(stderr,"                             2 Italic \n" );
fprintf(stderr,"                             3 Bold \n" );
fprintf(stderr,"                             4 Symbol (Greek) \n" );
fprintf(stderr,"-Wsleeptime  (default 0)   Delay between screen erase\n" );
fprintf(stderr,"-T4025      (default off)  Initialize Tektronix 4025\n");
fprintf(stderr,"-TT         (default off)  Initialize Teraterm Emulator\n");
fprintf(stderr,"-K          (default off)  Initialize Kermit 3.00 Emulator\n");
fprintf(stderr,"-h                         Do not execute, show options \n" );
fprintf(stderr,"-?                         Do not execute, show options \n" );
	exit(0);
}

void dv_clip()
{
}
