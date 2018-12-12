/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLOTGIF                                               c
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
	10 Jan 2001 - removed unused calls to dv_filltri and dv_fillbox
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

#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)

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
	short numxpixels;	/* Number of pixels in x axis	*/
	short numypixels;	/* Number of pixels in y axis	*/
	short numtextcols;	/* Number of text columns available	*/
	short numtextrows;	/* Number of text rows available */
	short numcolors;	/* Number of actual colors	*/
	short bitsperpixel;	/* Number of bits representing a pixel */
	short numvideopages;	/* Number of available video pages */
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

static char 	**arr = NULL;
static int	hw = 0;
static int	NxB;
static int	Ny;
static int	ccolor;

static int pagedirty = 0;	/* flag to indicate that page is not blank */

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
static void gifsetpalette() ,
	gifpoint(int x, int y),
	gifxline(int y, int xs, int xe),
	gifyline(int x, int ys, int ye),
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
static void putsi(short a,FILE *grfp);
static void usage(void);

/* GIF Stuff */
/*
 * Pointer to function returning an int
 */

#ifdef SIGNED_COMPARE_SLOW
typedef unsigned long int count_int;
typedef unsigned short int count_short;
#else /*SIGNED_COMPARE_SLOW*/
typedef long int          count_int;
#endif /*SIGNED_COMPARE_SLOW*/
/*
 * a code_int must be able to hold 2**BITS values of type int, and also -1
 */
typedef int             code_int;

typedef int (* ifunptr)();
void compress( int init_bits, FILE *outfile, ifunptr ReadValue );
void cl_block (void);
void GIFEncode(FILE *fp, int GWidth, int GHeight, int GInterlace, 
	int Background, int Bitsperpixel, 
	int *Red, int *Green, int *Blue, ifunptr Getpixel );
void Putword(int w, FILE *fp );
void writeerr(void);
void cl_hash(count_int hsize);
void output(code_int code );
void char_init( void);
void char_out( int c );
void flush_char(void );

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
	if(config.numcolors == 2){
		NxB = config.numxpixels/8;
		point = gifpoint;
		xline = gifxline;
		yline = gifyline;
		setcolor = setcolor2;
		gphshift = 3;
		gphmaskshift = 7;
		gphbitmask = 0x07;
		gphcharmask = 0xFF7F;
		gphbits = 1;
		kolor = 0;
	} else if(config.numcolors == 4){
		NxB = config.numxpixels/4;
		point = gifpoint;
		xline = gifxline;
		yline = gifyline;
		setcolor = setcolor4;
		gphshift = 2;
		gphmaskshift = 6;
		gphbitmask = 0x03;
		gphcharmask = 0xFF3F;
		gphbits = 2;
		kolor = 0;
	} else if(config.numcolors == 16){
		NxB = config.numxpixels/2;
		point = gifpoint;
		xline = gifxline;
		yline = gifyline;
		setcolor = setcolor256;
		gphshift = 1;
		gphmaskshift = 4;
		gphbitmask = 0x01;
		gphcharmask = 0xFF0F;
		gphbits = 4;
	} else if(config.numcolors == 256){
		NxB = config.numxpixels;
		point = gifpoint;
		xline = gifxline;
		yline = gifyline;
		setcolor = setcolor256;
		gphshift = 0;
		gphmaskshift = 0;
		gphbitmask = 0x00;
		gphcharmask = 0x0;
		gphbits = 8;
	}
	Ny = config.numypixels;
	arr = ( char **) calloc (config.numypixels, sizeof(char * ) );
	while ( hw< config.numypixels) {	/* get some rasters */
		if( (arr[hw++] = (char *)calloc(NxB,sizeof(char)))==NULL){
			fprintf(stderr,"lpdriver:alloc failed at %d\n",
				--hw);
			exit( 3 );
		}
	}
	maxcolors = (int)config.numcolors;
	/* maxcolors must be set before defining the palette */
	gifsetpalette();
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

static void gifsetpalette()
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
	if(pagedirty)
		pcputimage();
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

static void gifpoint(int x,int y)
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
	mask = gphcharmask;
	color = (char)gphcolor << gphmaskshift;
	lfcp = &arr[y][x>>gphshift];
	bitpos = x & gphbitmask;
	/* default for 8 bit */
	for(j=0;j<bitpos;j++){
		mask  >>= gphbits;
		color >>= gphbits;
	}
	*lfcp = (*lfcp & (char)mask) | color;
}

static void gifxline(int y,int xs,int xe)
{
	int x;
	if(xe < xs){
		x  = xs;
		xs = xe;
		xe = x ;
	}
	for(x=xs;x<=xe;x++)
		gifpoint(x,y);
}

static void gifyline(int x,int ys,int ye)
{
	int y;
	if(ye < ys){
		y  = ys;
		ys = ye;
		ye = y ;
	}
	for(y=ys;y<=ye;y++)
		gifpoint(x,y);
}
	


char *gifown = "The Graphics Interchange Format(c) is the Copyright property of CompuServe Incorporated. GIF(sm) is a Service Mark property of CompuServe Incorporated." ;

void di_gcmdln(int argc,char *argv[])
{
	char *cp;
	double atof();
	/* DEFAULTS */
		config.numcolors = 2;
		config.numxpixels = 640;
		config.numypixels = 480;
	while(argc-- > 1 ) {
		if(*argv[1] == '-'){
			switch(argv[1][1]){
			case 'V':
				fprintf(stderr,"CALPLOT (3.0) COPYRIGHT (C) 1989 Saint Louis University\n");
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
				else if(strncmp(cp,"320",3)==0)
					config.numxpixels = 320;
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
				if(strncmp(cp,"350",3)==0)
					config.numypixels = 350;
				else if(strncmp(cp,"600",3)==0)
					config.numypixels = 600;
				else if(strncmp(cp,"200",2)==0)
					config.numypixels = 200;
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


static void putsi(short a,FILE *grfp)
{
        putc((char)a,grfp);
        putc((char)(a>>8),grfp);
}

int Bitsperpixel;

static void pcputimage()
{
	int i ;
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

    /* All set, let's do it. */
    GIFEncode(
	stdout, cols, rows, 0, 0, Bitsperpixel, Red, Green, Blue, Getpixel );

    exit( 0 );
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


/* PUBLIC */
/*
**
** Based on GIFENCOD by David Rowley <mgardi@watdscu.waterloo.edu>.A
** Lempel-Zim compression based on "compress".
**/
/*****************************************************************************
 *
 * GIFENCODE.C    - GIF Image compression interface
 *
 * GIFEncode( FName, GHeight, GWidth, GInterlace, Background,
 *            Bitsperpixel, Red, Green, Blue, Getpixel )
 *
 *****************************************************************************/


#define TRUE 1
#define FALSE 0

int Width, Height;
int curx, cury;
long CountDown;
int Pass = 0;
int Interlace;

/*
 * Bump the 'curx' and 'cury' to point to the next pixel
 */
void BumpPixel()
{
        /*
         * Bump the current X position
         */
        curx++;

        /*
         * If we are at the end of a scan line, set curx back to the beginning
         * If we are interlaced, bump the cury to the appropriate spot,
         * otherwise, just increment it.
         */
        if( curx == Width ) {
                curx = 0;

                if( !Interlace )
                        cury++;
                else {
                     switch( Pass ) {

                       case 0:
                          cury += 8;
                          if( cury >= Height ) {
                                Pass++;
                                cury = 4;
                          }
                          break;

                       case 1:
                          cury += 8;
                          if( cury >= Height ) {
                                Pass++;
                                cury = 2;
                          }
                          break;

                       case 2:
                          cury += 4;
                          if( cury >= Height ) {
                             Pass++;
                             cury = 1;
                          }
                          break;

                       case 3:
                          cury += 2;
                          break;
                        }
                }
        }
}

/*
 * Return the next pixel from the image
 */
GIFNextPixel( getpixel )
ifunptr getpixel;
{
        int r;

        if( CountDown == 0 )
                return EOF;

        CountDown--;

        r = ( * getpixel )( curx, cury );

        BumpPixel();

        return r;
}

/* public */

void GIFEncode(FILE *fp, int GWidth, int GHeight, int GInterlace, 
	int Background, int Bitsperpixel, 
	int *Red, int *Green, int *Blue, ifunptr Getpixel )
{
        int B;
        int RWidth, RHeight;
        int LeftOfs, TopOfs;
        int Resolution;
        int ColorMapSize;
        int InitCodeSize;
        int i;

        Interlace = GInterlace;

        ColorMapSize = 1 << Bitsperpixel;

        RWidth = Width = GWidth;
        RHeight = Height = GHeight;
        LeftOfs = TopOfs = 0;

        Resolution = Bitsperpixel;

        /*
         * Calculate number of bits we are expecting
         */
        CountDown = (long)Width * (long)Height;

        /*
         * Indicate which pass we are on (if interlace)
         */
        Pass = 0;

        /*
         * The initial code size
         */
        if( Bitsperpixel <= 1 )
                InitCodeSize = 2;
        else
                InitCodeSize = Bitsperpixel;

        /*
         * Set up the current x and y position
         */
        curx = cury = 0;

        /*
         * Write the Magic header
         */
        fwrite( "GIF87a", 1, 6, fp );

        /*
         * Write out the screen width and height
         */
        Putword( RWidth, fp );
        Putword( RHeight, fp );

        /*
         * Indicate that there is a global colour map
         */
        B = 0x80;       /* Yes, there is a color map */

        /*
         * OR in the resolution
         */
        B |= (Resolution - 1) << 5;

        /*
         * OR in the Bits per Pixel
         */
        B |= (Bitsperpixel - 1);

        /*
         * Write it out
         */
        fputc( B, fp );

        /*
         * Write out the Background colour
         */
        fputc( Background, fp );

        /*
         * Byte of 0's (future expansion)
         */
        fputc( 0, fp );

        /*
         * Write out the Global Colour Map
         */
        for( i=0; i<ColorMapSize; i++ ) {
                fputc( Red[i], fp );
                fputc( Green[i], fp );
                fputc( Blue[i], fp );
        }

        /*
         * Write an Image separator
         */
        fputc( ',', fp );

        /*
         * Write the Image header
         */

        Putword( LeftOfs, fp );
        Putword( TopOfs, fp );
        Putword( Width, fp );
        Putword( Height, fp );

        /*
         * Write out whether or not the image is interlaced
         */
        if( Interlace )
                fputc( 0x40, fp );
        else
                fputc( 0x00, fp );

        /*
         * Write out the initial code size
         */
        fputc( InitCodeSize, fp );

        /*
         * Go and actually compress the data
         */
        compress( InitCodeSize+1, fp, Getpixel );

        /*
         * Write out a Zero-length packet (to end the series)
         */
        fputc( 0, fp );

        /*
         * Write the GIF file terminator
         */
        fputc( ';', fp );

        /*
         * And close the file
         */
        fclose( fp );

}

/*
 * Write out a word to the GIF file
 */
void Putword(int w, FILE *fp )
{
        fputc( w & 0xff, fp );
        fputc( (w / 256) & 0xff, fp );
}


/***************************************************************************
 *
 *  GIFCOMPR.C       - GIF Image compression routines
 *
 *  Lempel-Ziv compression based on 'compress'.  GIF modifications by
 *  David Rowley (mgardi@watdcsu.waterloo.edu)
 *
 ***************************************************************************/

/*
 * General DEFINEs
 */

#define BITS    12

#define HSIZE  5003            /* 80% occupancy */

#ifdef NO_UCHAR
 typedef char   char_type;
#else /*NO_UCHAR*/
 typedef        unsigned char   char_type;
#endif /*NO_UCHAR*/

/*
 *
 * GIF Image compression - modified 'compress'
 *
 * Based on: compress.c - File compression ala IEEE Computer, June 1984.
 *
 * By Authors:  Spencer W. Thomas       (decvax!harpo!utah-cs!utah-gr!thomas)
 *              Jim McKie               (decvax!mcvax!jim)
 *              Steve Davies            (decvax!vax135!petsd!peora!srd)
 *              Ken Turkowski           (decvax!decwrl!turtlevax!ken)
 *              James A. Woods          (decvax!ihnp4!ames!jaw)
 *              Joe Orost               (decvax!vax135!petsd!joe)
 *
 */
#include <ctype.h>
/* #include <signal.h> */

#define ARGVAL() (*++(*argv) || (--argc && *++argv))

int n_bits;                        /* number of bits/code */
int maxbits = BITS;                /* user settable max # bits/code */
code_int maxcode;                  /* maximum code, given n_bits */
code_int maxmaxcode = (code_int)1 << BITS; /* should NEVER generate this
code */
#ifdef COMPATIBLE               /* But wrong! */
# define MAXCODE(n_bits)        ((code_int) 1 << (n_bits) - 1)
#else /*COMPATIBLE*/
# define MAXCODE(n_bits)        (((code_int) 1 << (n_bits)) - 1)
#endif /*COMPATIBLE*/

count_int htab [HSIZE];
unsigned short codetab [HSIZE];
#define HashTabOf(i)       htab[i]
#define CodeTabOf(i)    codetab[i]

code_int hsize = HSIZE;                 /* for dynamic table sizing */

/*
 * To save much memory, we overlay the table used by compress() with those
 * used by decompress().  The tab_prefix table is the same size and type
 * as the codetab.  The tab_suffix table needs 2**BITS characters.  We
 * get this from the beginning of htab.  The output stack uses the rest
 * of htab, and contains characters.  There is plenty of room for any
 * possible stack (stack used to be 8000 characters).
 */

#define tab_prefixof(i) CodeTabOf(i)
#define tab_suffixof(i)        ((char_type *)(htab))[i]
#define de_stack               ((char_type *)&tab_suffixof((code_int)1<<BITS))

code_int free_ent = 0;                  /* first unused entry */

/*
 * block compression parameters -- after all codes are used up,
 * and compression rate changes, start over.
 */
int clear_flg = 0;

int offset;
long int in_count = 1;            /* length of input */
long int out_count = 0;           /* # of codes output (for debugging) */

/*
 * compress stdin to stdout
 *
 * Algorithm:  use open addressing double hashing (no chaining) on the
 * prefix code / next character combination.  We do a variant of Knuth's
 * algorithm D (vol. 3, sec. 6.4) along with G. Knott's relatively-prime
 * secondary probe.  Here, the modular division first probe is gives way
 * to a faster exclusive-or manipulation.  Also do block compression with
 * an adaptive reset, whereby the code table is cleared when the compression
 * ratio decreases, but after the table fills.  The variable-length output
 * codes are re-sized at this point, and a special CLEAR code is generated
 * for the decompressor.  Late addition:  construct the table according to
 * file size for noticeable speed improvement on small files.  Please direct
 * questions about this implementation to ames!jaw.
 */

int g_init_bits;
FILE *g_outfile;

int ClearCode;
int EOFCode;

void compress( int init_bits, FILE *outfile, ifunptr ReadValue )
{
    register long fcode;
    register code_int i = 0;
    register int c;
    register code_int ent;
    register code_int disp;
    register code_int hsize_reg;
    register int hshift;

    /*
     * Set up the globals:  g_init_bits - initial number of bits
     *                      g_outfile   - pointer to output file
     */
    g_init_bits = init_bits;
    g_outfile = outfile;

    /*
     * Set up the necessary values
     */
    offset = 0;
    out_count = 0;
    clear_flg = 0;
    in_count = 1;
    maxcode = MAXCODE(n_bits = g_init_bits);

    ClearCode = (1 << (init_bits - 1));
    EOFCode = ClearCode + 1;
    free_ent = ClearCode + 2;

    char_init();

    ent = GIFNextPixel( ReadValue );

    hshift = 0;
    for ( fcode = (long) hsize;  fcode < 65536L; fcode *= 2L )
        hshift++;
    hshift = 8 - hshift;                /* set hash code range bound */

    hsize_reg = hsize;
    cl_hash( (count_int) hsize_reg);            /* clear hash table */

    output( (code_int)ClearCode );

#ifdef SIGNED_COMPARE_SLOW
    while ( (c = GIFNextPixel( ReadValue )) != (unsigned) EOF ) {
#else /*SIGNED_COMPARE_SLOW*/
    while ( (c = GIFNextPixel( ReadValue )) != EOF ) {
#endif /*SIGNED_COMPARE_SLOW*/

        in_count++;

        fcode = (long) (((long) c << maxbits) + ent);
        i = (((code_int)c << hshift) ^ ent);    /* xor hashing */

        if ( HashTabOf (i) == fcode ) {
            ent = CodeTabOf (i);
            continue;
        } else if ( (long)HashTabOf (i) < 0 )      /* empty slot */
            goto nomatch;
        disp = hsize_reg - i;           /* secondary hash (after G. Knott) */
        if ( i == 0 )
            disp = 1;
probe:
        if ( (i -= disp) < 0 )
            i += hsize_reg;

        if ( HashTabOf (i) == fcode ) {
            ent = CodeTabOf (i);
            continue;
        }
        if ( (long)HashTabOf (i) > 0 )
            goto probe;
nomatch:
        output ( (code_int) ent );
        out_count++;
        ent = c;
#ifdef SIGNED_COMPARE_SLOW
        if ( (unsigned) free_ent < (unsigned) maxmaxcode) {
#else /*SIGNED_COMPARE_SLOW*/
        if ( free_ent < maxmaxcode ) {
#endif /*SIGNED_COMPARE_SLOW*/
            CodeTabOf (i) = free_ent++; /* code -> hashtable */
            HashTabOf (i) = fcode;
        } else
                cl_block();
    }
    /*
     * Put out the final code.
     */
    output( (code_int)ent );
    out_count++;
    output( (code_int) EOFCode );

    return;
}

/*****************************************************************
 * TAG( output )
 *
 * Output the given code.
 * Inputs:
 *      code:   A n_bits-bit integer.  If == -1, then EOF.  This assumes
 *              that n_bits =< (long)wordsize - 1.
 * Outputs:
 *      Outputs code to the file.
 * Assumptions:
 *      Chars are 8 bits long.
 * Algorithm:
 *      Maintain a BITS character long buffer (so that 8 codes will
 * fit in it exactly).  Use the VAX insv instruction to insert each
 * code in turn.  When the buffer fills up empty it and start over.
 */

unsigned long cur_accum = 0;
int cur_bits = 0;

unsigned long masks[] = { 0x0000, 0x0001, 0x0003, 0x0007, 0x000F,
                                  0x001F, 0x003F, 0x007F, 0x00FF,
                                  0x01FF, 0x03FF, 0x07FF, 0x0FFF,
                                  0x1FFF, 0x3FFF, 0x7FFF, 0xFFFF };

void output(code_int code )
{
    cur_accum &= masks[ cur_bits ];

    if( cur_bits > 0 )
        cur_accum |= ((long)code << cur_bits);
    else
        cur_accum = code;

    cur_bits += n_bits;

    while( cur_bits >= 8 ) {
        char_out( (unsigned int)(cur_accum & 0xff) );
        cur_accum >>= 8;
        cur_bits -= 8;
    }

    /*
     * If the next entry is going to be too big for the code size,
     * then increase it, if possible.
     */
   if ( free_ent > maxcode || clear_flg ) {

            if( clear_flg ) {

                maxcode = MAXCODE (n_bits = g_init_bits);
                clear_flg = 0;

            } else {

                n_bits++;
                if ( n_bits == maxbits )
                    maxcode = maxmaxcode;
                else
                    maxcode = MAXCODE(n_bits);
            }
        }

    if( code == EOFCode ) {
        /*
         * At EOF, write the rest of the buffer.
         */
        while( cur_bits > 0 ) {
                char_out( (unsigned int)(cur_accum & 0xff) );
                cur_accum >>= 8;
                cur_bits -= 8;
        }

        flush_char();

        fflush( g_outfile );

        if( ferror( g_outfile ) )
                writeerr();
    }
}

/*
 * Clear out the hash table
 */
void cl_block (void)             /* table clear for block compress */
{

        cl_hash ( (count_int) hsize );
        free_ent = ClearCode + 2;
        clear_flg = 1;

        output( (code_int)ClearCode );
}

void cl_hash(hsize)          /* reset code table */
register count_int hsize;
{

        register count_int *htab_p = htab+hsize;

        register long i;
        register long m1 = -1;

        i = hsize - 16;
        do {                            /* might use Sys V memset(3) here */
                *(htab_p-16) = m1;
                *(htab_p-15) = m1;
                *(htab_p-14) = m1;
                *(htab_p-13) = m1;
                *(htab_p-12) = m1;
                *(htab_p-11) = m1;
                *(htab_p-10) = m1;
                *(htab_p-9) = m1;
                *(htab_p-8) = m1;
                *(htab_p-7) = m1;
                *(htab_p-6) = m1;
                *(htab_p-5) = m1;
                *(htab_p-4) = m1;
                *(htab_p-3) = m1;
                *(htab_p-2) = m1;
                *(htab_p-1) = m1;
                htab_p -= 16;
        } while ((i -= 16) >= 0);

        for ( i += 16; i > 0; i-- )
                *--htab_p = m1;
}

void writeerr(void)
{
        printf( "error writing output file\n" );
        exit(1);
}

/******************************************************************************
 *
 * GIF Specific routines
 *
 ******************************************************************************/

/*
 * Number of characters so far in this 'packet'
 */
int a_count;

/*
 * Set up the 'byte output' routine
 */
void  char_init(void  )
{
        a_count = 0;
}

/*
 * Define the storage for the packet accumulator
 */
char accum[ 256 ];

/*
 * Add a character to the end of the current packet, and if it is 254
 * characters, flush the packet to disk.
 */
void char_out( int c )
{
        accum[ a_count++ ] = c;
        if( a_count >= 254 )
                flush_char();
}

/*
 * Flush the packet to disk, and reset the accumulator
 */
void flush_char(void )
{
        if( a_count > 0 ) {
                fputc( a_count, g_outfile );
                fwrite( accum, 1, a_count, g_outfile );
                a_count = 0;
        }
}

/* The End */

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
fprintf(stderr,"plotgif [options]\n");
fprintf(stderr,"-V                         Program Version \n" );
fprintf(stderr,"-Sscalefac (default=1.0)   Plot magnifier  \n" );
fprintf(stderr,"-R         (default off)   Rotate plot 90 degrees \n" );
fprintf(stderr,"-N         (default off)   Turn off shading \n" );
fprintf(stderr,"-I         (default off)   Invert background (e.g., black)\n");
fprintf(stderr,"-Ffont     (default 0)     Default font \n" );
fprintf(stderr,"                             0 Roman \n" );
fprintf(stderr,"                             1 Roman \n" );
fprintf(stderr,"                             2 Italic \n" );
fprintf(stderr,"                             3 Bold \n" );
fprintf(stderr,"                             4 Symbol (Greek) \n" );
fprintf(stderr,"-Xnumx     (default 640)   X-pixels one of 640,800,320,1000,2000\n");
fprintf(stderr,"-Ynumy     (default 480)   Y-pixels one of 480,600,200,800,1600 \n");
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
