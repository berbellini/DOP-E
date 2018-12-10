/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLOTPNG                                               c
c                                                                     c
c      COPYRIGHT (C)  2003 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
CHANGES
	01 Jan 2003 - created based on plotgif
	02 APR 2004 - inplemented gend(mode) at higher level which
		required change in dv_closepl
*/
/* device dependent plot program interface 			*/
/* Plot rasterizer for fast plot on PC				*/
/* It is assumed that the plotgen program clips and requires	*/
/*	the very low level primitives				*/
/*	supports the calls 					*/
/*		dv_clip						*/
/*		dv_erase					*/
/*		dv_fillp					*/
/*		dv_font						*/
/*		di_gcmdln					*/
/*		di_gsymb					*/
/*		onintr						*/
/*		dv_openpl					*/
/*		dv_pen						*/
/*		dv_zpoint						*/
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
#include	<ctype.h>
#ifdef MSDOS
#include	<dos.h>
#include	<malloc.h>
#include	<stdlib.h>
#define FAR	far 
#include	<io.h>
#endif
#include	<fcntl.h>

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

/* Global Function Prototypes	*/
void dv_clip();
void		dv_closepl(int mode);
void		dv_erase(INT mode);
void		dv_fillp(INT n, INT *x,INT *y);
extern void	di_fillr(INT X0,INT Y0,INT X1,INT Y1,
			INT patx,INT paty,INT lenx,INT leny);
void		dv_font(INT fontvalue);
void		di_gsymb(INT x,INT y,INT ht,INT ang,INT nchar,char *s);
extern INT	dv_lineclip(INT *x0,INT *z0,INT *x1,INT *z1,
			INT MinX,INT MinY,INT MaxX,INT MaxY);
void		dv_openpl(int showmenu);
void		dv_pen(INT p);
extern void	dv_symvec(INT xloc,INT yloc,INT ht,char *s,INT ang,INT nochar);
void		dv_zpoint(int x,int y);
void		dv_zzpoint(int x,int y);

/* Global Variables		*/
INT dc_left	= 0;	/* device plotting limits */
INT dc_right	= 319;
INT dc_bottom	= 0;
INT dc_top	=199;
double dc_Nxpi	= 409.0;
double dc_Nypi	= 409.0;
INT dc_oldx1	= -1;
INT dc_oldy1 	= -1;

int dc_curpen = 1 ;	/* current pen value */
int	dc_mapcurpen;	/* mapping into dither for dv_zzpoint */
	/* flags and values set by command line argument		*/
int	scaledplot = 0;
INT	dc_rotate = 0;
INT	dc_shadeoff = 0 ;	/* permit shading unless turned 
					off in command line */
INT	dc_color = 0 ;	/* if = 0 no hardware fast color shading */
INT	dc_ColorInfo = 2;
INT	dc_xlinewidth = 0;
INT	dc_ylinewidth = 0;
INT	dc_linewidth = 0;
INT	dc_herelinewidth = 1; /* set linewidth here */
INT	dc_newlinewidth = 0;	/* flag to indicate a line width change */
INT	dc_hardwarefill = 0;	/* flag to have all filling done by hardware */
INT	dc_curfont = 0;
int	defaultfont = 0;
extern INT dc_oldx;	/* current plot point position */
extern INT dc_oldy;	/* needed by gottxt() and gintxt() */
INT dc_ClipRight, dc_ClipTop, dc_ClipBottom, dc_ClipLeft;
INT dc_ClipRight_sv, dc_ClipTop_sv, dc_ClipBottom_sv, dc_ClipLeft_sv;
INT	dc_iseps = 0;	/* invoke special scaling for EPS */

double	dc_scalefac = 1.0;
int dc_hasmouse = 0;

static int Ldc_top, Ldc_right, Ldc_bottom, Ldc_left; /* for internal clipping */

/* Local  Variables		*/

	/* graphics structure */
struct videoconfig {
	int numxpixels;	/* Number of pixels in x axis	*/
	int numypixels;	/* Number of pixels in y axis	*/
	int numtextcols;	/* Number of text columns available	*/
	int numtextrows;	/* Number of text rows available */
	int numcolors;	/* Number of actual colors	*/
	int bitsperpixel;	/* Number of bits representing a pixel */
	int numvideopages;	/* Number of available video pages */
	} ;
static struct videoconfig config;

static int savpen = 1;	/* current color */
static int dc_dim2 =  255;	/* dimensions of dither matrix dv_zzpoint */
static int maxcolors = 1;
static int black = 1;
static int kolor = 1;
static int Kolor = 1;
static int Whiten = 0;		/* lighten blues */
static int	Maxcol = 128;	/* 128 unique colors */
static int	CBLK = 0;	/* map position for preset colors */
static int	CWHIT= 1;
static int	CRED = 1;
static int	CORNG = 1;
static int	CYEL = 1;
static int	CGRN = 1;
static int	CBLGR= 1;
static int	CBLU = 1;
static int revvideo = 0;/* 0 gottxt straight write out strips trailing blank */
			/* 1 gottxt writes blanks then characters */
			/* 2 reversed video with trailing blanks */

#define MAXCOLORS 256
static int Red[MAXCOLORS], Green[MAXCOLORS], Blue[MAXCOLORS];
static int Getpixel();

static int	gphshift = 0,
	gphmaskshift = 0,
	gphbitmask = 0,
	gphcharmask = 0,
	gphbits = 0;
static int	gphcolor;


/*	arr is a pointer to an array of Nx pointers to 		*/
/*		raster buffers					*/
/*	hw is the number of buffers created so far (for very 	*/
/*		lone plots					*/

static int  **arr = NULL;
static int	hw = 0;
static int	NxB;
static int	Ny;
static int	ccolor;

static int pagedirty = 0;	/* flag to indicate that page is not blank */
int	Num	= 1 ;		/* output only this page for multipage */
int	count	= 1 ;		/* current page number */

/* Local  Function Prototypes	*/
static void	(*clearscreen)()=0,
		(*getvideoconfig)(struct videoconfig *config)=0,
		(*imagesize)()=0,
		(*point)(int x, int y)=0,
		(*xline)(int y, int xs, int xe)=0,
		(*yline)(int x, int ys, int ye)=0,
		(*setcolor)(int color, int pen)=0,
		(*putimage)()=0,
		(*getimage)()=0;
	/* define specific functions to accomplish this all */
static void rassetpalette() ,
	raspoint(int x, int y),
	rasxline(int y, int xs, int xe),
	rasyline(int x, int ys, int ye),
	setcolor2(),
	setcolor4(),
	setcolor16(),
	setcolor256();
	/* raster output */
static void pcputimage();
static void prgerr();
static void (*pointz)(int x, int y) = 0;

static void setcolorrgb(int Maxcol);
static void coord(int *red,int *green,int *blue,int index,int Maxcol);
static float value(float n1, float n2, float hue);
static void hls_to_rgb(float *r,float *g,float *b,float h,float l,float s);
static void setcolor2(int color,int ipen);
static void setcolor4(int color,int ipen);
static void setcolor16();
static void setcolor256(int color,int ipen);
static void usage(void);

/* PNG definitions */
#define PROGNAME  "plotpng"
#define VERSION   "1.00 of 01 January 2003"
#define APPNAME   "Simple CALPLOT to PNG Converter"
#include <string.h>
#include <setjmp.h>     /* for jmpbuf declaration in writepng.h */
#include <time.h>
#include "writepng.h"
int rc;
ulg rowbytes;
int error = 0;
int text = FALSE;


/* local prototypes */

static int  wpng_isvalid_latin1(uch *p, int len);
static void wpng_cleanup(void);
static mainprog_info wpng_info;   /* lone global */






void onintr(int arg)
{
	dv_closepl(0);
	exit(0);
}

void dv_closepl(int mode)
{
	if(pagedirty)dv_erase(0);
}


void dv_openpl(int showmenu)
{
	int numx, numy, i, j;
	int c;
	char *getenv(), *sptr;
	/* check environment for terminal definition */
	/* note di_gcmdln is called first in plotdriver.c */
	/* allocate data storage for graphics arrays  */
#ifdef MSDOS
	setmode(fileno(stdout), O_BINARY);
	setmode(fileno(stdin ), O_BINARY);
#endif
		NxB = config.numxpixels;
		point = raspoint;
		xline = rasxline;
		yline = rasyline;
		setcolor = setcolor256;
		gphshift = 0;
		gphmaskshift = 0;
		gphbitmask = 0x00;
		gphcharmask = 0x0;
		gphbits = 8;
	Ny = config.numypixels;
	arr = ( int **) calloc (config.numypixels, sizeof(int * ) );
	while ( hw< config.numypixels) {	/* get some rasters */
		if( (arr[hw++] = (int *)calloc(NxB,sizeof(int)))==NULL){
			fprintf(stderr,"lpdriver:alloc failed at %d\n",
				--hw);
			exit( 3 );
		}
	}
	maxcolors = (int)config.numcolors;
	/* maxcolors must be set before defining the palette */
	rassetpalette();
	dc_right = (int)config.numxpixels - 1;
	dc_top   = (int)config.numypixels - 1;
	dc_Nxpi = (double)dc_right/10.00;
	dc_Nypi = (double)dc_top  / 8.00;
	dv_pen((INT)1);
	dc_ClipLeft   = dc_left;
	dc_ClipRight  = dc_right;
	dc_ClipTop    = dc_top;
	dc_ClipBottom = dc_bottom;
	dc_ClipLeft_sv   = dc_left;
	dc_ClipRight_sv  = dc_right;
	dc_ClipTop_sv    = dc_top;
	dc_ClipBottom_sv = dc_bottom;
	Ldc_top = (int)dc_top;
	Ldc_right = (int)dc_right;
	Ldc_left = (int)dc_left;
	Ldc_bottom = (int)dc_bottom;
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

}

static void rassetpalette()
{
	int i;
	if(maxcolors == 2){
		/* BLACK */
		Blue[0] = Green[0] = Red[0] = 0;
		/* WHITE */
		Blue[1] = Green[1] = Red[1] = 255;
	} else if(maxcolors == 4){
		/* BLACK */
		Blue[0] = Green[0] = Red[0] = 0;
		/* RED */
		Blue[1] = 0;
		Green[1] = 0;
		Red[1] = 255;
		/* BLUE */
		Blue[2] = 255;
		Green[2] = 0;
		Red[2] = 0;
		/* WHITE */
		Blue[3] = Green[3] = Red[3] = 255;
	} else if(maxcolors == 16){
		setcolorrgb(16);
	} else if(maxcolors == 256){
		setcolorrgb(256);
	}
}


void dv_erase(INT mode)
{
	int x, y;
	if(pagedirty && count==Num )
		pcputimage();
	if(pagedirty)count++;
	/* erase the page */
	for (x=0;x< NxB;x++){
		for(y=0;y<Ny;y++){
			arr[y][x] = 0;
		}
	}
	pagedirty = 0;
	dv_pen((INT)1);
}


void dv_zpoint(int x,int y)
{
	(*point)(x,Ldc_top-y);
	pagedirty = 1;
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
		dc_mapcurpen = dc_dim2 - dc_mapcurpen;
	}
	if(dc_curpen >= 1000 && kolor == 0)return;
	if(ipen == 0 )
		(*setcolor)((int)0,(int)0); 
	else if(ipen < 0){
		ipen = 1;
		savpen = ipen;
	} else if(ipen < 1000){
		ipen = ipen%(maxcolors-1) ;
		if(ipen == 0)ipen = (maxcolors -1);
		savpen = ipen;
		(*setcolor)((int)ipen,(int)ipen);
	} else if(ipen >= 1000){
		index = (int)( (maxcolors-3) * (float)(ipen-1000)/100.0) +2;
		if(index >= maxcolors)index = maxcolors - 1;
		(*setcolor)((int)index,(int)ipen);
	}
}


/* set up color tables for this device */
static void setcolorrgb(int Maxcol)
{
	int red, green, blue;
	int  index;

	for(index=0;index<Maxcol;index++){
		coord(&red,&green,&blue,index,Maxcol);
		Red[index] = red;
		Green[index] = green;
		Blue[index] = blue;
	}
	if(Maxcol > 1){
		CBLK = 0;
		CWHIT = 1;
		CRED = 2;
		CBLU = Maxcol - 1;
		CGRN = (int)( (float)(CBLU+CRED)*0.5 + 0.5);
		CYEL = (int)( (float)(CBLU + 3.*CRED)*0.25 + 0.5);
		CORNG= (int)( (float)(CBLU + 7.*CRED)*0.125 + 0.5);
		CBLGR= (int)( (float)(3.*CBLU + CRED)*0.25 + 0.5);
	}
}

static void coord(int *red,int *green,int *blue,int index,int Maxcol)
{
	float rred, rgrn, rblu;
	float l = 0.5;	/* lumination	*/
	float s = 1.0;	/* saturation	*/
	float h = 0.0 ; /* hue		*/
	float fac ;	/* map index to gray non linearly */
	if(index == 0){
		rred = 0.0;
		rgrn = 0.0;
		rblu = 0.0;
	} else if(index == 1){
		rred = 1.0;
		rgrn = 1.0;
		rblu = 1.0;
	} else {
		if(kolor < 2){
			if(Kolor == 1){
				l = 0.5;
				s = 1.0 ;
				h = 240.0*(float)(index-2)/(float)(Maxcol-3);
				hls_to_rgb(&rred,&rgrn,&rblu,h,l,s);
			} else {
				/* construct scale Red -> Green -> Blue */
				fac = (float)(index-2)/(float)(Maxcol-3) ;
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
			fac = (float)(index-2)/(float)(Maxcol-3) ;
			if(black==0)
				fac = 0.2 + 0.7*fac;
			else
				fac = 0.2 + 0.7*fac;
			rred = 1.0-fac;
			rgrn = 1.0-fac;
			rblu = 1.0-fac;
			}
	}
	*red = (int)(255.0*rred);
	if(*red < 0)
		*red = 0;
	else if(*red > 255)
		*red = 255;
	*green= (int)(255.0*rgrn);
	if(*green < 0)
		*green = 0;
	else if(*green > 255)
		*green = 255;
	*blue= (int)(255.0*rblu);
	if(*blue < 0)
		*blue = 0;
	else if(*blue > 255)
		*blue = 255;
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




static void setcolor2(int color,int ipen)
{
	static int pt[] =
		{ 0,	/* black */
		1,	/* white */
		0 };
	if(ipen == 0)
		gphcolor = pt[0];
	else
		gphcolor = pt[1];
}

static void setcolor4(int color,int ipen)
{
	static int pt[] =
		{ 0,	/* black */
		3,	/* white */
		1,	/* red */
		2,	/* blue */
		0 };
	if(ipen == 0)
		gphcolor = 0;
	else {
		color = (color%(maxcolors - 1)) + 1;
		gphcolor = pt[color];
	}
}
static void setcolor16()
{
}
static void setcolor256(int color,int ipen)
{
	if(ipen != color)
		gphcolor = color;
	else {
		if(ipen == 0)
			gphcolor = CBLK;
		else {
			ipen = (ipen-1)%7 + 1;
			if(ipen == 1){
				gphcolor=CWHIT;
			} 
			/* line black on white background, 
				line white on black background */
			else if(ipen == 2){
				gphcolor=CRED;
			} else if(ipen == 3) {
				gphcolor=CGRN;
			} else if(ipen == 4) {
				gphcolor=CBLU;	/* blue */
			} else if(ipen == 5) {
				gphcolor=CORNG;	/* orange */
			} else if(ipen == 6) {
				gphcolor=CBLGR;	/* blue green */
			} else if(ipen == 7) {
				gphcolor=CYEL;	/* yellow */
			}
		}
	}
}

static void raspoint(int x,int y)
{
	char *lfcp;
	char c;
	unsigned char color;
	int bitpos,j,mask;
        unsigned char temp;

	if(x < Ldc_left)return;
	if(x > Ldc_right)return;
	if(y < Ldc_bottom)return;
	if(y > Ldc_top)return;
	arr[y][x] = gphcolor;
}

static void rasxline(int y,int xs,int xe)
{
	int x;
	if(xe < xs){
		x  = xs;
		xs = xe;
		xe = x ;
	}
	for(x=xs;x<=xe;x++)
		raspoint(x,y);
}

static void rasyline(int x,int ys,int ye)
{
	int y;
	if(ye < ys){
		y  = ys;
		ys = ye;
		ye = y ;
	}
	for(y=ys;y<=ye;y++)
		raspoint(x,y);
}
	



void di_gcmdln(int argc,char *argv[])
{
	char *cp;
	double atof();
	/* DEFAULTS */
		config.numcolors = 256;
		config.numxpixels = 640;
		config.numypixels = 480;
	while(argc-- > 1 ) {
		if(*argv[1] == '-'){
			switch(argv[1][1]){
			case 'V':
				fprintf(stderr,"CALPLOT (3.25) COPYRIGHT (C) 2003 Saint Louis University\n");
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
				Num = atoi(&argv[1][2]);
				break;
			case 'I':
				black = 0;
				break;
			case 'K':
				kolor = 1;
				switch(argv[1][2]){
				case 'W':
					Whiten = 1;
					break;
				case 'R':
					Whiten = 0;
					Kolor = 2;
					break;
				case 'B':
					Whiten = 0;
					Kolor = 3;
					break;
				default:
					Whiten = 0;
					break;
				}
				break;
			case 'G':
				kolor = 2;
				break;
			case 'F':	/* default font value */
				cp = argv[1];
				cp++;
				cp++;
				defaultfont = atoi(cp);
				if(defaultfont < 0)defaultfont = 0;
				dc_curfont = defaultfont;
				break;
			case 'C':
				cp = argv[1];
				cp++;
				cp++;
				config.numcolors = 2;
				if(strncmp(cp,"256",3)==0)
					config.numcolors = 256;
				else if(strncmp(cp,"16",2)==0)
					config.numcolors = 16;
				else if(strncmp(cp,"4",1)==0)
					config.numcolors = 4;
				else if(strncmp(cp,"2",1)==0)
					config.numcolors = 2;
				break;
			case 'X':
				cp = argv[1];
				cp++;
				cp++;
				config.numxpixels = 640;
				if(strncmp(cp,"800",3)==0)
					config.numxpixels = 800;
				else if(strncmp(cp,"400",3)==0)
					config.numxpixels = 400;
				else if(strncmp(cp,"1000",3)==0)
					config.numxpixels = 1000;
				else if(strncmp(cp,"2000",3)==0)
					config.numxpixels = 2000;
				break;
			case 'Y':
				cp = argv[1];
				cp++;
				cp++;
				config.numypixels = 480;
				if(strncmp(cp,"600",3)==0)
					config.numypixels = 600;
				else if(strncmp(cp,"320",2)==0)
					config.numypixels = 320;
				else if(strncmp(cp,"800",2)==0)
					config.numypixels = 800;
				else if(strncmp(cp,"1600",2)==0)
					config.numypixels = 1600;
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
}

struct errstr {
	char *item;
	} ;

struct errstr err[] = {
	" ",
	(char *)0 
	} ;

static void prgerr()
{
	struct errstr *eptr;
	for(eptr = err ; eptr->item != (char *)NULL ; eptr++)
		fprintf(stderr,"%s\n",eptr->item);
	exit(0);
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
	INT des, lst ;
	INT savpen;
	if(nchar > 0){
		if(revvideo){
			savpen = dc_curpen;
			lst = nchar * ht;
			des =  ht/2;
			if(revvideo == 1)
				dv_pen((INT)0);
			else
				dv_pen((INT)1);
			di_fillr(x -10,y-des,x+lst,y+ht+des/2,(INT)0,(INT)0,(INT)0,(INT)0);
			if(revvideo == 1)
				dv_pen(savpen);
			else
				dv_pen((INT)0);
		}
		dv_symvec(x,y,ht,s,ang,nchar);
		if(revvideo){
			dc_curpen = savpen;
			revvideo = 0;
			dv_pen((INT)dc_curpen);
		}
	} else {
		dv_symvec(x,y,ht,s,ang,nchar);
	}
}


void dv_fillp(INT n, INT *x,INT *y)
{
}


void dv_rliner(INT x0,INT z0,INT x1,INT z1)
{
	int x,y;
	int count,error;
	int ypx, ymx;
	int deltax,deltay;
	INT tmp;
	/* invoke rotation if required */
	if(dc_rotate){
		tmp = z0;
		z0 = x0;
		x0 = dc_right - tmp;
		tmp = z1;
		z1 = x1 ;
		x1 = dc_right - tmp;
	}
	if(dv_lineclip(&x0,&z0,&x1,&z1,dc_ClipLeft,dc_ClipBottom,dc_ClipRight,dc_ClipTop))
		return;
	if(dc_curpen < 1000)
		pointz=dv_zpoint;
	else if(dc_curpen >= 1000 && kolor >= 1){
		pointz=dv_zpoint;
		}
	else {
		pointz=dv_zzpoint;
	}
	error = 0;
	deltax = (int)(x1 - x0);
	deltay = (int)(z1 - z0);
	if(deltax == 0 &&  dc_curpen < 1000 )
		(*yline)((int)x1,(int)(dc_top-z0),(int)(dc_top-z1));
	else if( deltax == 0 && kolor == 1 && dc_curpen >= 1000 )
		(*yline)((int)(x1),(int)(dc_top-z0),(int)(dc_top-z1));
	else if(deltay == 0 && dc_curpen < 1000 )
		(*xline)((int)(dc_top-z0),(int)(x0),(int)(x1));
	else if ( deltay == 0 && kolor == 1 && dc_curpen >= 1000 )
		(*xline)((int)(dc_top-z0),(int)(x0),(int)(x1));
	else {
		if(deltay < 0)
		{
			tmp = x0;
			x0 = x1;
			x1 = tmp;
			tmp = z0;
			z0 = z1;
			z1 = tmp;
			deltax = -deltax;
			deltay = -deltay;
		}
		ypx = deltax+deltay;
		ymx = deltay-deltax;
		x = (int)x0;
		y = (int)z0;
		(*pointz) (x,y);
		if(deltax >= 0){			/* positive slope */
			if(deltax > deltay){		/* 0 < slope < 1 */
				for(count =1;count < deltax;count++){
					if(error <=  0){
						x++;
						(*pointz)(x,y);
						error+= deltay;
					}
					else{
						x++;
						y++;
						(*pointz)(x,y);
						error+=ymx;
					}
				}
			}
			else if(deltax < deltay) {	/* slope > 1 */
				for(count=1;count < deltay;count++){
					if(error <  0){
						x++;
						y++;
						(*pointz)(x,y);
						error+= ymx;
					}
					else {
						y++;
						(*pointz)(x,y);
						error-=deltax;
					}
				}
			}
			else {				/* slope = 1 */
				x++;
				y++;
				while(x<(int)x1){
					(*pointz)(x++,y++);
				}
			}
		}
		else {					/* negative slope */
			if(-deltax > deltay){		/* -1 < slope < 0 */
				for(count=1;count < -deltax;count++){
					if(error <=  0){
						x--;
						(*pointz)(x,y);
						error+=deltay;
					}
					else {
						x--;
						y++;
						(*pointz)(x,y);
						error+= ypx;
					}
				}
			}
			else if (-deltax < deltay) {	/* slope < -1 */
				for(count=1;count<deltay;count++){
					if(error <  0){
						x--;
						y++;
						(*pointz)(x,y);
						error+= ypx;
					}
					else {
						y++;
						(*pointz)(x,y);
						error+=deltax;
					}
				}
			}
			else {
				x--;
				y++;
				while(x>(int)x1){
					(*pointz)(x--,y++);
				}
			}
		}
		(*pointz)((int)x1,(int)z1);
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
		if(bit[x%dim][y%dim] >= dc_mapcurpen)
			dv_zpoint(x,y);
	}
}

void dv_cursor(INT curstyp)
{
}



int Bitsperpixel;

static void pcputimage()
{
	int i ;
	int k;
    int cols, rows;
	/* reverse the sense of black and white if the image is to
		have reversed background
		ASSUME THAT POSITION 0 is black in current color
		map and that Maxcolors-1 is white  */
	if(black == 0){
		Red[0] = 0;
		Green[0] = 0;
		Blue[0] = 0;
		Red[ 1] = 255;
		Green[1] = 255;
		Blue[1] = 255;
	} else {
		Red[0] = 255;
		Green[0] = 255;
		Blue[0] = 255;
		Red[ 1] = 0;
		Green[1] = 0;
		Blue[1] = 0;
	}

    Bitsperpixel = gphbits;
    cols = config.numxpixels;
    rows = config.numypixels;
	wpng_info.infile = NULL;
	wpng_info.outfile = stdout;
	wpng_info.image_data = NULL;
	wpng_info.row_pointers = NULL;
	wpng_info.filter = FALSE;
	wpng_info.interlaced = FALSE;
	wpng_info.have_bg = FALSE;
	wpng_info.have_time = FALSE;
	wpng_info.have_text = 0;
	wpng_info.gamma = 0.0;
	wpng_info.modtime = time(NULL);
	wpng_info.width = cols;
	wpng_info.height = rows;
	wpng_info.sample_depth = 8;
	/* force RGB */
	wpng_info.pnmtype = 6 ;
/* allocate libpng stuff, initialize transformations, write pre-IDAT data */
	if ((rc = writepng_init(&wpng_info)) != 0) {
		switch (rc) {
			case 2:
				fprintf(stderr, PROGNAME ":  libpng initialization problem (longjmp)\n");
			break;
			case 4:
				fprintf(stderr, PROGNAME ":  insufficient memory\n");
			break;
			case 11:
				fprintf(stderr, PROGNAME ":  internal logic error (unexpected PNM type)\n");
			break;
			default:
				fprintf(stderr, PROGNAME ":  unknown writepng_init() error\n");
			break;
		}
		exit(rc);
	}
	/* calculate rowbytes on basis of image type; note that this becomes much
	*      * more complicated if we choose to support PBM type, ASCII PNM types, or
	*           * 16-bit-per-sample binary data [currently not an official NetPBM type] */
	if (wpng_info.pnmtype == 5)
		rowbytes = wpng_info.width;
	else if (wpng_info.pnmtype == 6)
		rowbytes = wpng_info.width * 3;
	else /* if (wpng_info.pnmtype == 8) */
		rowbytes = wpng_info.width * 4;


	/* not interlaced:  write progressively (row by row) */ 
        long j;
        ulg bytes;

        wpng_info.image_data = (uch *)calloc(rowbytes, sizeof(uch ));
        if (wpng_info.image_data == NULL) {
            fprintf(stderr, PROGNAME ":  insufficient memory for row data\n");
            writepng_cleanup(&wpng_info);
            wpng_cleanup();
            exit(5);
        }
        error = 0;
        for (j=0 ; j < wpng_info.height;   j++) {

		bytes=0L;
		for(k=0;k<cols;k++){
			wpng_info.image_data[bytes++] = (uch)(Red  [arr[j][k]]);
			wpng_info.image_data[bytes++] = (uch)(Green[arr[j][k]]);
			wpng_info.image_data[bytes++] = (uch)(Blue [arr[j][k]]);
		}
            if (bytes != rowbytes) {
                fprintf(stderr, PROGNAME
                  ":  expected %lu bytes, got %lu bytes (row %ld)\n", rowbytes,
                  bytes, j);
                ++error;
                break;
            }
            if (writepng_encode_row(&wpng_info) != 0) {
                fprintf(stderr, PROGNAME
                  ":  libpng problem (longjmp) while writing row %ld\n",j);
                ++error;
                break;
            }
        }
        if (error) {
            writepng_cleanup(&wpng_info);
            wpng_cleanup();
            exit(2);
        }
        if (writepng_encode_finish(&wpng_info) != 0) {
            fprintf(stderr, PROGNAME ":  error on final libpng call\n");
            writepng_cleanup(&wpng_info);
            wpng_cleanup();
            exit(2);
        }




	writepng_cleanup(&wpng_info);
	wpng_cleanup();

    exit( 0 );
}


static int wpng_isvalid_latin1(uch *p, int len)
{
    int i, result = -1;

    for (i = 0;  i < len;  ++i) {
        if (p[i] == 10 || (p[i] > 31 && p[i] < 127) || p[i] > 160)
            continue;           /* character is completely OK */
        if (result < 0 || (p[result] != 27 && p[i] == 27))
            result = i;         /* mark location of first questionable one */
    }                           /*  or of first escape character (bad) */

    return result;
}





static void wpng_cleanup(void)
{
    if (wpng_info.outfile) {
        fclose(wpng_info.outfile);
        wpng_info.outfile = NULL;
    }

    if (wpng_info.infile) {
        fclose(wpng_info.infile);
        wpng_info.infile = NULL;
    }

    if (wpng_info.image_data) {
        free(wpng_info.image_data);
        wpng_info.image_data = NULL;
    }

    if (wpng_info.row_pointers) {
        free(wpng_info.row_pointers);
        wpng_info.row_pointers = NULL;
    }
}




/* assume that x is to the right */
static int Getpixel( x, y )
int x, y;
    {
    int color;
	unsigned int mask;
	if(Bitsperpixel == 1){
		mask = (0x0080) >> (x%8) ;
		color=(int)(((arr[y][x/8]&(unsigned char)mask)>>(7-(x%8)))&01);
	} else if(Bitsperpixel == 2){
		mask = (0x00C0) >> (2*(x%4)) ;
		color=(int)(((arr[y][x/4]&(unsigned char)mask )>>(6-(2*(x%4))))&03);
	} else if(Bitsperpixel == 4){
		if( x%2 == 0)
			color = (int)(( arr[y][x/2] & 0xF0 ) >> 4 );
		else
			color = (int)( arr[y][x/2] & 0xF );
	} else if(Bitsperpixel == 8)
		color = (int)(arr[y][x] &0xFF);
    return color;
}




/* get text string from terminal */
void di_gintxt(int cnt,char *s)
{
	s[0]='\0';
}
void di_cross(INT *ix, INT *iy, char *c)
{
	*ix = 0;
	*iy = 0;
	c[0] = '\0';
}

void dv_mesg(char *mesg)
{
}

static void usage(void)
{
fprintf(stderr,"plotpng [options]\n");
fprintf(stderr,"-V                         Program Version \n" );
fprintf(stderr,"-Sscalefac (default=1.0)   Plot magnifier  \n" );
fprintf(stderr,"-R         (default off)   Rotate plot 90 degrees \n" );
fprintf(stderr,"-Nnum     (default 1)      Convery page num \n" );
fprintf(stderr,"-I         (default off)   Invert background (e.g., black)\n");
fprintf(stderr,"-Ffont     (default 0)     Default font \n" );
fprintf(stderr,"                             0 Roman \n" );
fprintf(stderr,"                             1 Roman \n" );
fprintf(stderr,"                             2 Italic \n" );
fprintf(stderr,"                             3 Bold \n" );
fprintf(stderr,"                             4 Symbol (Greek) \n" );
fprintf(stderr,"-Xnumx     (default 640)   X-pixels one of 640,800,400,1000,2000\n");
fprintf(stderr,"-Ynumy     (default 480)   Y-pixels one of 480,600,320,800,1600 \n");
fprintf(stderr,"-K         (default color) Color output \n" );
fprintf(stderr,"-KW        (default -K   ) Color output, whitened spectrum\n" );
fprintf(stderr,"-KR         (default -K  ) Color output Red->White->Blue \n" );
fprintf(stderr,"-KB         (default -K  ) Color output Blue->White->Red \n" );
fprintf(stderr,"-G          (default -K  ) Gray output \n" );
fprintf(stderr,"-Ccolors   (default 2)     Size of Colormap 2, 4, 16 or 256\n");
fprintf(stderr,"-h                         Do not execute, show options \n" );
fprintf(stderr,"-?                         Do not execute, show options \n" );
	exit (0);
}

void dv_clip()
{
}
