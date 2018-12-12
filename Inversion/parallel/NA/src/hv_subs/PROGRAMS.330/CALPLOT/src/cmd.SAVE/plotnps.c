/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLOTNPS                                               c
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
	21 AUG 2004 - removed MSDOS #ifdef changed %d to %d
*/
/* device dependent plot program interface 			*/
/* PostScript Native Language					*/
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
#include	<math.h>
#include	<stdlib.h>
#include	<stddef.h>
#include	<string.h>
static char	pstream[100];
FILE	*popen();
#define	INT	int

#define	MAXLINE 1000
#define	OUTSTR(x)	(void)fprintf(stream,x );
#define	MAX(a,b) ( (a) > (b) ? (a):(b) )
#define	MIN(a,b) ( (b) > (a) ? (a):(b) )

/* Global Function Prototypes	*/
void dv_clip();
void dv_closepl(int mode);
void dv_concur(INT isx,INT isy);
void dv_fillp(INT n, INT *x,INT *y);
void di_gsymb(INT x,INT y,INT ht,INT ang,INT nchar,char *s);
void dv_openpl(int showmenu);
void dv_movcur(INT isx,INT isy);
void dv_pen(INT p);
void dv_pen(INT p);
void dv_zpoint(INT x,INT y);
void di_cont(INT x, INT y);
void di_control(int type, int i1, int i2, int i3, int i4);

/* Global Variables		*/
	/* device clipping region */
INT dc_left = 0, dc_right, dc_top,  dc_bottom = 0;
double dc_Nxpi	= 1000.0;
double dc_Nypi	= 1000.0;
INT dc_oldx1	= -1;
INT dc_oldy1 	= -1;
	/* flags and values set by command line argument		*/
INT	dc_rotate = 0;
INT	dc_shadeoff = 0 ;	/* permit shading unless turned off */
INT	dc_color = 3 ;	/* hardware color shading implemented by gray scale */
INT	dc_ColorInfo = 2;

int dc_curpen = 1 ;	/* current pen value */
extern INT dc_dim2;	/* dimensions of dither matrix dv_zzpoint */
int	dc_mapcurpen;	/* mapping into dither for dv_zzpoint */

INT	dc_xlinewidth = 0;
INT	dc_ylinewidth = 0;
INT	dc_linewidth = 0;
INT	dc_minlinewidth = 0;
INT	dc_herelinewidth = 0; /* set linewidth here in this routine */
INT	dc_newlinewidth = 1;	/* flag to indicate a line width change */
INT	dc_hardwarefill = 2; 	/* flag to have all filling done by hardware */
INT	dc_curfont = 0;
double	dc_scalefac = 1.0;
int dc_hasmouse = 0;
extern INT dc_shdon ;
extern INT dc_shdcur;
extern INT dc_shdse;


INT dc_ClipRight, dc_ClipTop, dc_ClipBottom, dc_ClipLeft;
INT dc_ClipRight_sv, dc_ClipTop_sv, dc_ClipBottom_sv, dc_ClipLeft_sv;
INT	dc_iseps = 0;	/* invoke special scaling for EPS */
static INT cxl, cyl, cxh, cyh;


/* Local  Function Prototypes	*/
static void con (INT isx,INT isy);
static void coord(float *red,float *green,float *blue,INT index);
static void do_plaid(INT xmn,INT ymn,INT xmx,
		INT ymx,INT patx,INT paty,INT lenx,INT leny);
static void endpage(void);
static void eps_rotate_box(INT x,INT y,INT ht,INT angle,INT nchar);
static void eps_bounding_box( INT x, INT y);
static void eps_box_rotate(INT *tx,INT *ty,
		float ct,float st,INT x0,INT y0,INT x,INT y);
static void hls_to_rgb(float *r,float *g,float *b,float h,float l,float s);
static void mov (INT isx,INT isy);
static void newpage(void);
static void psinit(void );
static void rotateout(INT *x,INT *y);
static void show_clip(INT lx,INT ly,INT ux,INT uy);
static void setcolor(void);
static void set_fill_index(INT index);
static void set_line_index(INT index);
static void symps(INT x,INT y,INT ht,char *s,INT ang,INT nchar);
static void usage(void);

 
/* Local  Variables		*/
static INT	scaledplot = 0;
static INT	pltstream = 0;
static INT	black = 0;	/* paint onto white paper */
static INT	white = 1;
static double grayshade = 0.0;
static INT grayline = 0.0;

static INT	linecolor = 0;
static INT	shadecolor = 0;
static INT	defaultfont = 0;
static INT	half30 = 0;		/* set halftone screen to 30/inch */
static INT Whiten = 0;		/* lighten blues */
static INT	Maxcol = 128;	/* 128 unique colors */
static float   red_array[128];
static float green_array[128];
static float  blue_array[128];
static INT	Kolor =1;
static INT	CBLK = 0;	/* map position for preset colors */
static INT	CWHIT= 1;
static INT	CRED = 1;
static INT	CORNG = 1;
static INT	CYEL = 1;
static INT	CGRN = 1;
static INT	CBLGR= 1;
static INT	CBLU = 1;
static INT	revvideo = 0;
			/* to specify if page is dirty */
static INT	pagedirty = 0;		
			/* to specify if within a line draw */
static INT	Inline	= 0;		
			/* page number of output */
static INT	Page_Number = 0;		
static INT	criticalio = 0;
static INT	pleasequit = 0;
static INT iseps = 0;		/* flag for EPS document */
static INT size_paper = 0;	/* flag for 11x14 paper document 
				0 8.5 x 11
				1 11  x 14
				2 8.5  x 14
				3 A3 842x1191 11.694x16.542
				4 A4 595x 842  8.264x11.694
				 */
static INT pagenumber=0;	/* page number of PS output */
static INT epsrotate=0;		/* rotation for EPS */
static INT bblx, bbly, bbux, bbuy;
static INT dc_top_eps, dc_right_eps, dc_bottom_eps=0, dc_left_eps=0;

static long linecount = 0;	/* to guarantee limitcheck 
				not exceeded on printer*/

static	FILE	*stream = NULL;
#define NOSTR 256
static char ostr[NOSTR];		/* for title bar at bottom of page */



static float value(float n1, float n2, float hue);
/* define PostScript fonts and their maximum number */
#define NFONTS 12
static char *fonts[NFONTS+1] = {
	"/Times-Roman",
	"/Times-Roman",
	"/Times-Italic",
	"/Times-Bold",
	"/Symbol",
	"/Helvetica",
	"/Helvetica-Oblique",
	"/Helvetica-Bold",
	"/Symbol",
	"/Courier",
	"/Courier-Oblique",
	"/Courier-Bold",
	"/Symbol"
	};
static INT fonts_used[NFONTS+1] =
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

void onintr(int arg)
{
	if(criticalio){
		pleasequit = 1;
	} else {
		exit (1);
	}
}

void dv_closepl(int mode)
{
	INT num_fonts_used, i,first;
	/* count the number of fonts used */
	for(i=1, num_fonts_used=0;i <= NFONTS; i++)
		if(fonts_used[i] > 0)num_fonts_used++;

	endpage();

	OUTSTR("%%%%Trailer\n")
	(void)fprintf(stream,"%%%%Pages: %d\n",Page_Number);
	if(num_fonts_used > 0){
		OUTSTR("%%%%DocumentNeededResources:")
		for(i=1,first=0;i<=NFONTS;i++){
			if(fonts_used[i] > 0){
				if(first == 0){
					fprintf(stream," font %s\n",fonts[i]);
					first++;
				} else {
					first = 1;
					fprintf(stream,"%%%%+ font %s\n",fonts[i]);
				}
			}
		}
	} else {
		OUTSTR("%%%%DocumentNeededResources: \n")
	}
	if(iseps) {
		if(bbuy > dc_top_eps  )bbuy = dc_top_eps;
		if(bbux > dc_right_eps)bbux = dc_right_eps;
		if(bbly > dc_top_eps  )bbly = dc_top_eps;
		if(bblx > dc_right_eps)bblx = dc_right_eps;
		
		if(bbuy < dc_bottom_eps  )bbuy = dc_bottom_eps;
		if(bbux < dc_left_eps    )bbux = dc_left_eps;
		if(bbly < dc_bottom_eps  )bbly = dc_bottom_eps;
		if(bblx < dc_left_eps    )bblx = dc_left;

		/* safety check */
		if(bbuy <= bbly)
			bbly = bbuy;
		if(bbux <= bblx)
			bblx = bbux;
		/* put in a margin of 0.125" = 3.175 mm = 9 PostScript units */
		bbux += 125; 
		bbuy += 125;
		bblx -= 125;
		bbly -= 125;
		
		bblx = (INT)( (float)bblx * 0.072);
		bbly = (INT)( (float)bbly * 0.072);
		bbux = (INT)( (float)bbux * 0.072);
		bbuy = (INT)( (float)bbuy * 0.072);
		(void)fprintf(stream,"%%%%BoundingBox: %d %d %d %d\n",
			bblx,bbly,bbux,bbuy);
		OUTSTR("%%%%EOF\n")
	} else {
		OUTSTR("%%%%EOF\n")
	}
#ifndef MSDOS
	if(stream != NULL && stream != stdout)
		fclose(stream);
#endif
	stream = NULL;
}

void dv_openpl(int showmenu)
{
	if(stream == NULL){
#ifndef MSDOS
		if(pltstream){
			stream = popen(pstream,"w");
			if(stream == NULL){
				fprintf(stderr,"plotps:cannnot open pipe\n");
				exit (1);
			}
		} else 
#else
		/* setmode(fileno(stdout), O_BINARY); */
#endif
			stream = stdout;
		psinit();
		newpage();
		dc_iseps = iseps;
		if(iseps == 1 ){
			dc_top = 0000;
			dc_right = 1000000;
			dc_bottom = -1000000;
			dc_left = -1000000;
			dc_top_eps = 1000000;
			dc_right_eps = 1000000;
			dc_bottom_eps = -1000000;
			dc_left_eps = -1000000;
			show_clip(dc_bottom_eps,dc_left_eps,dc_top_eps,dc_right_eps);
		} else if(!iseps) {
			OUTSTR(" 0 0 translate\n") 
			dc_bottom = 0;
			dc_left = 0;
			if(size_paper=0){	/* Letter */
				dc_right = 11000 -125 ;
				dc_left = 125;
				dc_top   = 8500 -125 ;
				dc_bottom = 125;
			} else if(size_paper == 1){ /* Big */
				dc_right = 14000 -125 ;
				dc_left = 125;
				dc_top   = 11000 -125 ;
				dc_bottom = 125;
			} else if(size_paper == 2){ /* Big */
				dc_right = 14000 -125 ;
				dc_left = 125;
				dc_top   =  8500 -125 ;
				dc_bottom = 125;
			} else if(size_paper == 3){	/* A3 */
				dc_right = 16542 -125 ;
				dc_left = 125;
				dc_top   = 11694 -125 ;
				dc_bottom = 125;
			} else if(size_paper == 4){	/* A4 */
				dc_right = 11694 -125 ;
				dc_left = 125;
				dc_top   =  8264 -125 ;
				dc_bottom = 125;
			} else {
				dc_right = 11000 -125;
				dc_left = 125;
				dc_top = 8500 -125;
				dc_bottom = 125;
			}
			show_clip(dc_bottom,dc_left,dc_top,dc_right);
		}
	OUTSTR("newpath\n")
	setcolor();	
	dv_pen((INT)1);
	dc_shdon = 0;
	dc_shdcur = 0;
	}
/*
	dc_ClipTop =dc_top;
	dc_ClipBottom = dc_bottom;
	dc_ClipLeft = dc_left;
	dc_ClipRight = dc_right;
	dc_ClipTop_sv =dc_top;
	dc_ClipBottom_sv = dc_bottom;
	dc_ClipLeft_sv = dc_left;
	dc_ClipRight_sv = dc_right;
	cxl = dc_top - dc_ClipBottom;
	cyl = dc_ClipLeft;
	cxh = dc_top - dc_ClipTop;
	cyh = dc_ClipRight;
*/
	if(iseps == 1){
		dc_ClipTop =dc_top_eps;
		dc_ClipBottom = dc_bottom_eps;
		dc_ClipLeft = dc_left_eps;
		dc_ClipRight = dc_right_eps;
		dc_ClipTop_sv =dc_top_eps;
		dc_ClipBottom_sv = dc_bottom_eps;
		dc_ClipLeft_sv = dc_left_eps;
		dc_ClipRight_sv = dc_right_eps;
		cxl = dc_ClipBottom;
		cyl = dc_ClipLeft;
		cxh = dc_ClipTop;
		cyh = dc_ClipRight;
	} else {
		dc_ClipTop =dc_top;
		dc_ClipBottom = dc_bottom;
		dc_ClipLeft = dc_left;
		dc_ClipRight = dc_right;
		dc_ClipTop_sv =dc_top;
		dc_ClipBottom_sv = dc_bottom;
		dc_ClipLeft_sv = dc_left;
		dc_ClipRight_sv = dc_right;
		cxl = dc_top - dc_ClipBottom;
		cyl = dc_ClipLeft;
		cxh = dc_top - dc_ClipTop;
		cyh = dc_ClipRight;
	}


	/* dc_colorInfo 
		0	monochrome on white background
		1	gray on white background
		2	color on white background
		4	monochrome on black background
		5	gray on black background
		6	color on black background

	dc_color	0
			1	Color
			2	Grayscale
			3	Grayscale but lines are solid black
	*/

	
	if(dc_color == 0)
		dc_ColorInfo = 0;
	else if(dc_color == 1)
		dc_ColorInfo = 2;
	else
		dc_ColorInfo = 1;
	
}

void dv_erase(INT mode)
{
	if(pagedirty){
		endpage();
		newpage();
		/* START OUTPUT */
		if(half30==1)
		OUTSTR("30 45 { dup mul exch dup mul add 1 exch sub } setscreen\n")
		if(!iseps)
			OUTSTR("af\n")
		if(!iseps){
			OUTSTR(" 0 0 translate\n") 
			show_clip(dc_bottom,dc_left,dc_top,dc_right);
		} else if(iseps == 1){
			show_clip(dc_bottom_eps,dc_left_eps,dc_top_eps,dc_right_eps);
		}
		if(dc_linewidth < dc_minlinewidth)
			(void)fprintf(stream,"%d unit setlinewidth\n",dc_minlinewidth);
		else
			(void)fprintf(stream,"%d unit setlinewidth\n",dc_linewidth);
		OUTSTR("1 setlinecap 1 setlinejoin newpath\n")
	}
	dc_shdon = 0;
	dc_shdcur = 0;
	pagedirty = 0;
}

static void psinit()
{
	if(iseps){
		OUTSTR("%%!PS-Adobe-3.0 EPSF-3.0\n")
		OUTSTR("%%%%BoundingBox: (atend)\n")
		/* set absolute limits for BoundingBox */
		bblx = 100000000; bbly = 100000000; bbux = -100000000; bbuy = -100000000;
	} else {
		OUTSTR("%%!PS-Adobe-3.0\n")
		if(dc_rotate )
			OUTSTR("%%%%Orientation: Portrait\n")
		else
			OUTSTR("%%%%Orientation: Portrait\n")
	}
	OUTSTR("%%%%Creator: plotnps\n")
	OUTSTR("%%%%DocumentNeededResources: (atend)\n")
	OUTSTR("%%%%Pages: (atend)\n")
	OUTSTR("%%%%PageOrder: Ascend\n")
	OUTSTR("%%%%EndComments\n")

	OUTSTR("%%%%BeginProlog\n")
	if(!iseps){
		OUTSTR("initmatrix\n")
		OUTSTR("/ps { print flush } def\n")
		OUTSTR("/home { newpath 0 pgtop moveto } def\n")
		OUTSTR("/mf { statusdict /manualfeed true put\n")
		OUTSTR(" } def\n")
		OUTSTR("/af { statusdict /manualfeed false put } def\n")
		OUTSTR("/jobname (-) def\n")
		OUTSTR("userdict /jobname jobname put\n")
		OUTSTR("clippath pathbbox pop pop exch pop 0 exch translate\n")
		OUTSTR("clippath pathbbox /pgtop exch def pop pop pop\n")
	}
	OUTSTR("/ShowEqui 	%% HEIGHT NCHAR STRING ShowEqui\n")
	OUTSTR("	{\n")
	OUTSTR("		%% stack is HEIGHT NCHAR STRING\n")
	OUTSTR("		dup\n")
	OUTSTR("		%% stack is HEIGHT NCHAR STRING STRING \n")
	OUTSTR("		stringwidth pop exch\n")
	OUTSTR("		%% stack is HEIGHT NCHAR XWIDTH STRING\n")
	OUTSTR("		4 1 roll\n")
	OUTSTR("		%% stack is STRING HEIGHT NCHAR XWIDTH\n")
	OUTSTR("		exch div   \n")
	OUTSTR("		%% stack is STRING HEIGHT EXTRA_PER_CHARACTER\n")
	OUTSTR("		sub\n")
	OUTSTR("		%% stack is STRING EXTRAX\n")
	OUTSTR("		exch  0 exch  ashow\n")
	OUTSTR("	}\n")
	OUTSTR("	def\n")
	OUTSTR("/plaid {		\n")
	OUTSTR("		unit /xmin exch def 	\n")
	OUTSTR("		unit /dx exch def 	\n")
	OUTSTR("		unit /xmax exch def	\n")
	OUTSTR("		unit /ymin exch def 	\n")
	OUTSTR("		unit /dy exch def 	\n")
	OUTSTR("		unit /ymax exch def	\n")
	OUTSTR("		/strn exch def	\n")
	OUTSTR("	\n")
	OUTSTR("		xmin dx xmax 			%% loop over x	\n")
	OUTSTR("			{ ymin dy ymax 		%% loop over y	\n")
	OUTSTR("				{ exch dup 3 2 roll dx dy 4 2 roll	\n")
	OUTSTR("					maskit   	\n")
	OUTSTR("				} for 	\n")
	OUTSTR("				pop 	\n")
	OUTSTR("			} for	\n")
	OUTSTR("} def	\n")
	OUTSTR("	\n")
	OUTSTR("	\n")
	OUTSTR("/maskit {	\n")
	OUTSTR("		1 setgray	\n")
	OUTSTR("		gsave	\n")
	OUTSTR("		4 copy	\n")
	OUTSTR("		translate	\n")
	OUTSTR("		scale	\n")
	OUTSTR("		6 6 	\n")
	OUTSTR("		false	\n")
	OUTSTR("		[6 0 0 6 0 0]	\n")
	OUTSTR("		strn	\n")
	OUTSTR("		imagemask	\n")
	OUTSTR("		%% return to initial position	\n")
	OUTSTR("		4 2 roll	\n")
	OUTSTR("		1 exch div exch  1 exch div exch  	\n")
	OUTSTR("		scale	\n")
	OUTSTR("		neg exch neg exch 	\n")
	OUTSTR("		translate	\n")
	OUTSTR("		grestore	\n")
	OUTSTR("} def	\n")

	/* START OUTPUT */
	if(half30==1)
		OUTSTR("30 45 { dup mul exch dup mul add 1 exch sub } setscreen\n")
	OUTSTR("/unit { 0.072 mul } def \n")
	OUTSTR("/m { unit exch unit exch moveto } def \n")
	OUTSTR("/n { unit exch unit exch rlineto } def \n")
	OUTSTR("/l { unit exch unit exch lineto } def \n")
	OUTSTR("/sn { stroke newpath } def\n")
	OUTSTR("/cgf { closepath gsave  0 setgray fill grestore stroke } def\n")
	if(!iseps)
		OUTSTR("af\n")
	OUTSTR("1 setlinecap 1 setlinejoin newpath\n")
	OUTSTR("%%%%EndPrologue\n")
	pagedirty = 0;
}



static void endpage()
{
	if(Inline)
		OUTSTR("stroke\n")
	if(strlen(ostr)> 0 ){
	dv_pen((INT)1);
	di_gsymb((INT) 500,(INT) 100,(INT) 50,(INT) 0,(INT) strlen(ostr),ostr);
	}
	if(!iseps)OUTSTR(" 0 0 translate\n")
	OUTSTR("showpage\n")
}

static void newpage()
{
	Page_Number++;
	(void)fprintf(stream,"%%%%Page: %d %d\n", Page_Number,Page_Number);
	OUTSTR("%%%%BeginPageSetup\n")
	OUTSTR("%%%%EndPageSetup\n")
}
	


void dv_zpoint(INT x,INT y)
{
	dv_movcur(x,y);
	dv_concur(x,y);
}

static void rotateout(INT *x,INT *y)
{
	/* rotate x, y, coordinates for laserwriter */
	INT tmp;
	tmp = *y;
	*y = *x;
	*x = dc_top - tmp;
}

static void mov (INT isx,INT isy)
{
	pagedirty = 1;
	(void)fprintf(stream,"%d %d m\n",isx,isy);
}

static INT mx, my;	/* temporary move save */
void dv_movcur(INT isx,INT isy)
{
	linecount = 0L;
	rotateout(&isx, &isy);
	mx = isx;
	my = isy;
	if(dc_shdon != dc_shdcur ){
		if(dc_shdon == 1){
			/* OUTSTR("%% movcur toggle dc_shdon == 1 \n") */
			OUTSTR("sn\n")
		} else {
			/* OUTSTR("%% movcur toggle dc_shdon == 0 \n")  */
			if(dc_color > 0 ){
				fprintf(stream,
					"closepath gsave %5.3f %5.3f %5.3f setrgbcolor fill grestore stroke\n",
					red_array[shadecolor],green_array[shadecolor],
					blue_array[shadecolor]);
			} else
				OUTSTR("/cgf\n")
		}
		dc_shdcur = dc_shdon;
		Inline = 0;
	} else if(Inline){		/* terminate current line */
		OUTSTR("stroke\n");
		Inline = 0;
	}
}

static void con (INT isx,INT isy)
{
	(void)fprintf(stream,"%d %d n\n",isx,isy);
}

void dv_concur(isx,isy)
INT isx,isy;
{
	INT sisx, sisy;
	if(dc_shdon != dc_shdcur ){
		if(dc_shdon == 1){
			/* OUTSTR("%% concur toggle dc_shdon == 1 \n") */
			OUTSTR("stroke\n")
		} else {
			/* OUTSTR("%% concur toggle dc_shdon == 0 \n")  */
			if(dc_color > 0 ){
				fprintf(stream,
					"closepath gsave %5.3f %5.3f %5.3f setrgbcolor fill grestore stroke\n",
					red_array[shadecolor],green_array[shadecolor],
					blue_array[shadecolor]);
			} else
				OUTSTR("cgf\n")
			
		}
		dc_shdcur = dc_shdon;
		Inline = 0;
	}
	if(!Inline){
		OUTSTR("newpath\n")
		mov(mx,my);
		Inline = 1;
		if( dc_newlinewidth == 1){
			if(dc_linewidth < dc_minlinewidth)
				(void)fprintf(stream,"%d unit setlinewidth\n",dc_minlinewidth);
			else
				(void)fprintf(stream,"%d unit setlinewidth\n",dc_linewidth);
			dc_newlinewidth = 0;
		}
	}
	sisx = isx;
	sisy = isy;
	rotateout(&isx, &isy);
	con(isx-mx,isy-my);
	if(iseps){
		eps_bounding_box( mx, my);
		eps_bounding_box( isx, isy);
	}
	mx = isx;
	my = isy;
	linecount = linecount + 1L;
	if(linecount > MAXLINE)dv_movcur(sisx,sisy);
}

void dv_pen(INT jpen)
{
	INT index;
	int ipen;
	ipen = jpen;
	if(Inline && !dc_shdse){		/* terminate current line */
		OUTSTR("stroke\n");
		Inline = 0;
	}
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
	if(dc_curpen >= 1000 && dc_color == 0)return;
	if(dc_color == 3){
		if(dc_curpen >= 1000 && dc_curpen <= 1100){
			index =(INT)((Maxcol-2)*(float)(ipen-1000)/100.0 )+2;
			if(index >= Maxcol)index = Maxcol-1;
			set_fill_index(index); 
		} else {
			if(dc_curpen == 0){
				set_fill_index(1);
				set_line_index(1);
			} else {
				set_fill_index(0);
				set_line_index(0);
			}
		}
	} else if(dc_color >0 && dc_color < 3){
		if(black == 1){
			if(ipen == 0){
				set_fill_index(0);
				set_line_index(0);
			}
		} else {
			if(ipen == 0){
				set_fill_index(1);
				set_line_index(1);
			}
		}
		if(ipen == 0)return;
		if(ipen < 0)ipen = 1;
		if(ipen < 1000){
			ipen = (ipen-1)%7 + 1;
			if(ipen == 1){
				/* line black on white background */
				set_line_index(black);
				set_fill_index(black);
			} else if(ipen == 2){
				set_line_index(CRED);
				set_fill_index(CRED);
			} else if(ipen == 3) {
				set_line_index(CGRN);
				set_fill_index(CGRN);
			} else if(ipen == 4) {
				set_line_index(CBLU);
				set_fill_index(CBLU);
			} else if(ipen == 5) {
				set_line_index(CORNG);
				set_fill_index(CORNG);
			} else if(ipen == 6) {
				set_line_index(CBLGR);
				set_fill_index(CBLGR);
			} else if(ipen == 7) {
				set_line_index(CYEL);
				set_fill_index(CYEL);
			}
		} else if(ipen >= 1000){
			if(ipen > 1100)ipen=1100;
			if(Whiten)		/* whiten blues */
				ipen = 1000 + (INT)(0.75 * (float)(ipen-1000));
			index =(INT)((Maxcol-3)*(float)(ipen-1000)/100.0 )+2;
			set_fill_index(index); 
			set_line_index(index); 
		}
	}
}

static void set_fill_index(index)
INT index;
{
	shadecolor = index;
}

static void set_line_index(index)
INT index;
{
	linecolor = index;
	if(!dc_shdse)
	fprintf(stream,"%5.3f %5.3f %5.3f setrgbcolor\n",red_array[linecolor],
		green_array[linecolor],blue_array[linecolor]);
}

/* set up color tables for this device */
static void setcolor()
{
	float red, green, blue;
	INT  index;

	for(index=0;index<Maxcol;index++){
		coord(&red,&green,&blue,index);
		red_array[index] = red;
		green_array[index] = green;
		blue_array[index] = blue;
	}
	if(Maxcol > 1){
		CRED = 2;
		CBLU = Maxcol - 1;
		CGRN = (INT)( (float)(CBLU+CRED)*0.5 + 0.5);
		CYEL = (INT)( (float)(CBLU + 3.*CRED)*0.25 + 0.5);
		CORNG= (INT)( (float)(CBLU + 7.*CRED)*0.125 + 0.5);
		CBLGR= (INT)( (float)(3.*CBLU + CRED)*0.25 + 0.5);
	}
}

static void coord(float *red,float *green,float *blue,INT index)
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
		if(dc_color < 2){
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
	} else if(hue >= 60.0 && hue < 180.0){
		return( n2);
	} else if(hue >= 180.0 && hue < 240.0){
		return( n1 + (n2-n1)*(240.0-hue)/60.0 );
	} else if(hue >= 240.0){
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


void di_gcmdln(argc,argv)
int argc;
char *argv[];
{
	char *cp, *tp;
	double atof();
	while(argc-- > 1 ) {
		if(*argv[1] == '-'){
			switch(argv[1][1]){
			case 'V':
				fprintf(stderr,"CALPLOT (3.0) COPYRIGHT (C) 1997 R. B. Herrmann\n");
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
#ifndef MSDOS
			case 'P':
				pltstream = 1;
				cp = argv[1];
				cp++;
				cp++;
				strcpy(pstream,cp);
				break;
#endif
			case 'R':
				dc_rotate = 1;
				break;
			case 'N':
				dc_shadeoff = 1;
				break;
			case 'K':
				dc_color = 1;
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
			case 'B':
				size_paper = 1;
				break;
			case 'L':
				size_paper = 2;
				break;
			case 'A':
				cp = &argv[1][1];
				if(strncmp(cp,"A3",3)==0)
					size_paper = 3;
				else if(strncmp(cp,"A4",3)==0)
					size_paper = 4;
				break;
			case 'G':
				dc_color = 2;
				break;
			case 'F':	/* default font value */
				cp = argv[1];
				cp++;
				cp++;
				defaultfont = atoi(cp);
				if(defaultfont < 0)defaultfont = 0;
				if(defaultfont > 0)
					defaultfont = (defaultfont-1)%NFONTS+1;
				dc_curfont = defaultfont;
				break;
			case 'H':
				cp = &argv[1][1];
				if(strncmp(cp,"H30",3)==0)
					half30 = 1;
				else if(strncmp(cp,"H60",3)==0)
					half30 = 0;
				break;
			case 'W':
				cp = argv[1];
				cp++;
				cp++;
				dc_minlinewidth = atoi(cp);
				if(dc_minlinewidth < 0) dc_minlinewidth = 0;
			case 'E':
				cp = &argv[1][1];
				if(strncmp(cp,"EPS",3)==0)
					iseps = 1;
				break;
			case 'T':	/* title at bottom of page  */
				cp = argv[1];
				cp++;
				cp++;
				if(strlen(cp) < NOSTR)
					strcpy(ostr,cp);
				else
					strncpy(ostr,cp,NOSTR-1);
					ostr[NOSTR -1 ] = '\0';
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
	/* for EPS, we change the rotate flag - this is because
		we usually consider  the x-axis horizontal.
		The plotnps default  is  that x is in the long direction
		unless the -R flag is used.
		For the printed page, we make the x-horizontal,
		which is the short direction
	*/
	if(iseps==1){
		if(dc_rotate != 0){
			dc_rotate = 0;
		} else {
			dc_rotate = 1;
		}
	}

}

void dv_font(INT fontvalue)
{
	if(fontvalue > 0)
		dc_curfont = (fontvalue-1)%NFONTS + 1;
	else
		dc_curfont = defaultfont;
}



void di_gsymb(INT x,INT y,INT ht,INT ang,INT nchar,char *s)
{
	if(nchar > 0){
		symps(x,y,ht,s,ang,nchar);
	} else if(nchar < 0 ) {
		if((int)s[0] < 16)
			dv_symvec(x,y,ht,s,ang,nchar);
		else{
			if(nchar < -1)
				di_cont(x,y);
			s[0] = (char)((int)s[0]+16);
			symps(x,y,ht,s,ang,1);
		}
			
	}
}


void dv_fillp(INT n, INT x[],INT y[])
{
	INT i;
	if(n < 3)return;
	/* beware we must rotate coordinates */
		for(i = 0;i < n ; i++)
			rotateout(&x[i],&y[i]);
	if(Inline){
		OUTSTR("stroke\n")
	}
	Inline = 0;
	OUTSTR("newpath\n")
	mov(x[0],y[0]);
	if(iseps)eps_bounding_box( x[0], y[0]);
	for(i=1;i<n;i++){
		if(iseps)eps_bounding_box( x[i], y[i]);
		con(x[i]-x[i-1],y[i]-y[i-1]);
	}
	OUTSTR("closepath\n")
	Inline = 0;
	if(dc_color > 0 ){
		fprintf(stream,"%5.3f %5.3f %5.3f setrgbcolor\n",
			red_array[shadecolor],green_array[shadecolor],
			blue_array[shadecolor]);
		OUTSTR("fill \n")
		if(linecolor != shadecolor){
		fprintf(stream,"%5.3f %5.3f %5.3f setrgbcolor\n",
			red_array[linecolor],green_array[linecolor],
			blue_array[linecolor]);
		}
	}
}

extern double dc_xscale ;
extern double dc_yscale ;
extern double dc_xoff   ;
extern double dc_yoff   ;

#define IX(n)	( dc_xscale *  dc_scalefac * ( n + dc_xoff ) + 0.49)
#define IY(n)   ( dc_yscale *  dc_scalefac * ( n + dc_yoff ) + 0.49)

static void symps(INT x,INT y,INT ht,char *s,INT ang,INT nchar)
{
	INT  angle;
	int i, j;
	float height, fontsize;
	char *outstr;
	int olddc_curpen;
	INT lst,des;
	/* beware we must rotate coordinates */
	outstr = calloc(2*nchar+2,sizeof(char));
	olddc_curpen = dc_curpen;
	if(revvideo){
		lst = nchar * ht;
		des =  ht/2;
		if(revvideo == 1)
			dv_pen((INT)0);
		else if(revvideo == 2)
			dv_pen((INT)1);
		di_fillr(x -10,y-des,x+lst,y+ht+des/2,0,0,0,0);
	}
	x = IX(x);
	y = IY(y);
	if(!dc_rotate){
		rotateout(&x,&y);
		angle = ang + 90;
	} else {
		angle = ang;
	}
	if(Inline){
		OUTSTR("stroke\n")
	}
	Inline = 0;
	pagedirty = 1;
	OUTSTR("newpath\n")
	OUTSTR("gsave\n")
/*
*/
		show_clip(cxl,cyl,cxh,cyh);
	if(iseps)eps_rotate_box(x,y,ht,angle,nchar);
		if(revvideo == 2)
			dv_pen((INT)0);
		else
			dv_pen((INT)olddc_curpen);
	fprintf(stream,"%d %d m\n",x,y);
	fprintf(stream,"%d rotate\n",angle);
	height = dc_scalefac * (float)ht * 72.0 / 1000.0 ;
	fontsize = 1.5 * height ;
	if(dc_curfont < 0 || dc_curfont > NFONTS){
		OUTSTR("/Times-Roman findfont\n")
		fonts_used[1] = 1;
	} else {
		fprintf(stream,"%s findfont\n",fonts[dc_curfont]);
		fonts_used[dc_curfont] = 1;
	}
	fprintf(stream,"%f scalefont setfont\n",fontsize);
	for(i=0,j=0;i<nchar;i++){
		if(s[i] == '(')
			outstr[j++] = '\\';
		else if(s[i] == ')')
			outstr[j++] = '\\';
		else if(s[i] == '\\' )
			outstr[j++] = '\\';
		else if(s[i] == '\177' )
			s[i] = ' ';
		outstr[j++] = s[i];
	}
	outstr[j] = '\0';
	fprintf(stream,"%f %d (%s) ShowEqui\n",height,nchar,outstr);
	revvideo = 0;
	OUTSTR("grestore\n")
	dv_pen((INT)olddc_curpen);
	free((void *)outstr);
}

static void eps_bounding_box( INT x, INT y)
{
	/* the purpose of this is to find the smallest bounding box
		within the plot space */
	if(x < cxl ) x = cxl;
	if(x > cxh ) x = cxh;
	if(y < cyl ) y = cyl;
	if(y > cyh ) y = cyh;
	if( x >= dc_left_eps && 
	    x <= dc_right_eps && 
	    y >= dc_bottom_eps && 
	    y <= dc_top_eps ) {
		if( x < bblx ) bblx = x;
		if( y < bbly ) bbly = y;
		if( x > bbux ) bbux = x;
		if( y > bbuy ) bbuy = y;
/*
fprintf(stderr,"EBB (%d,%d) [%d %d %d %d]\n",x,y,cxl,cyl,cxh,cyh);
fprintf(stderr," BB: [%d %d %d %d]\n",bblx,bbly,bbux,bbuy);
*/

	}
}

static void eps_rotate_box(INT x,INT y,INT ht,INT angle,INT nchar)
{
	/* get the corners of the box defined by the
		symbol string */
	float degrad, ct, st;
	INT width;
	INT tx,ty;
	degrad = 3.1415927/180.0;
	ct = cos(angle*degrad);
	st = sin(angle*degrad);
	width = nchar*ht;
	eps_box_rotate(&tx,&ty,ct,st,x,y,    0,-ht/2);
	eps_bounding_box( tx, ty);
	eps_box_rotate(&tx,&ty,ct,st,x,y,width,-ht/2);
	eps_bounding_box( tx, ty);
	eps_box_rotate(&tx,&ty,ct,st,x,y,    0,ht);
	eps_bounding_box( tx, ty);
	eps_box_rotate(&tx,&ty,ct,st,x,y,width,ht);
	eps_bounding_box( tx, ty);
}

static void eps_box_rotate(INT *tx,INT *ty,float ct,float st,INT x0,INT y0,INT x,INT y)
{
	*tx = (INT)(ct*x -st*y + x0);
	*ty = (INT)(st*x +ct*y + y0);
}
	
static void show_clip(INT lx,INT ly,INT ux,INT uy)
{
	fprintf(stream,"%d %d m %d %d l %d %d l %d %d l closepath clip\n",
		lx,ly,ux,ly,ux,uy,lx,uy);
}

static unsigned int patmask[6] = { 0x80, 0x40, 0x20, 0x10, 0x8, 0x4 };
static char x_intmask[6] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20 };
static char y_intmask[6] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20 };
static void do_plaid(INT xmn,INT ymn,INT xmx,INT ymx,INT patx,INT paty,INT lenx,INT leny)
{
	INT j,k;
	unsigned INT pat[12] ;
	unsigned INT mask ;
	/* adjust pattern length to be compatible with mapping
		in PostScript code 
	*/
	lenx*= 6;
	leny*= 6;
	/* determine bounding box for paqttern */
	xmn /= lenx; xmn *= lenx;
	ymn /= leny; ymn *= leny;
	xmx /= lenx ; xmx++; xmx *= lenx;
	ymx /= leny ; ymx++; ymx *= leny;
	if(patx == 0)
		mask = 0xFC;
	else {
		for(j=0, mask = 0x0; j<6; j++){	
			if( (patx & 0x3F) & x_intmask[j] ){
		 	mask |=   patmask[j] ;
			}
		}
	}
	for(j=0,k=0; j<6 ; j++)
		if( ((paty ) & y_intmask[j]) || paty == 0 ) {
			pat[k++] = (mask >> 4) & 0xF ;
			pat[k++] = (mask & 0xF);
		} else {
			pat[k++] = 0x0;
			pat[k++] = 0x0;
		}
	/*
	PATTERN YMAX DY YMIN XMAX DX XMIN plaid
	*/
	if(!dc_rotate){
	/* for synchronization we must refer to dc_top */
	k = dc_top - xmx;
	k /= lenx;
	k = k * lenx;
	k = dc_top - k;
	fprintf(stream,"< %1X%1X %1X%1X %1X%1X %1X%1X %1X%1X %1X%1X > %d %d %d %d %d %d plaid \n",
		pat[0],pat[1],pat[2],pat[3],pat[4],pat[5],
		pat[6],pat[7],pat[8],pat[9],pat[10],pat[11],
		ymx,leny,ymn,xmn,-lenx,k );
	} else {
	fprintf(stream,"< %1X%1X %1X%1X %1X%1X %1X%1X %1X%1X %1X%1X > %d %d %d %d %d %d plaid \n",
		pat[0],pat[1],pat[2],pat[3],pat[4],pat[5],
		pat[6],pat[7],pat[8],pat[9],pat[10],pat[11],
		ymx,leny,ymn,xmx,lenx,xmn );
	}
}


void dv_cursor(INT curstyp)
{
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



void usage(void)
{
fprintf(stderr,"plotnps [options]\n");
fprintf(stderr,"-V                           Program Version \n" );
fprintf(stderr,"-Sscalefac  (default=1.0)    Plot magnifier  \n" );
#ifndef MSDOS
fprintf(stderr,"-Ppipe      (default stdout) Pipe output to process pipe \n" );
#endif
fprintf(stderr,"-R          (default off)    Rotate plot 90 degrees \n" );
fprintf(stderr,"-N          (default off)    Turn off shading \n" );
fprintf(stderr,"-Ffont      (default 0)      Default font \n" );
fprintf(stderr,"                             0 Times-Roman \n" );
fprintf(stderr,"                             1 Times-Roman \n" );
fprintf(stderr,"                             2 Times-Italic \n" );
fprintf(stderr,"                             3 Times-Bold \n" );
fprintf(stderr,"                             4 Symbol (Greek) \n" );
fprintf(stderr,"                             5 Helvetica \n" );
fprintf(stderr,"                             6 Helvetica-Oblique \n" );
fprintf(stderr,"                             7 Helvetica-Bold \n" );
fprintf(stderr,"                             8 Symbol (Greek) \n" );
fprintf(stderr,"                             9 Courier \n" );
fprintf(stderr,"                            10 Courier-Oblique \n" );
fprintf(stderr,"                            11 Courier-Bold \n" );
fprintf(stderr,"                            12 Symbol (Greek) \n" );
fprintf(stderr,"-H30        (default H60)    Halftone for shading larger dots \n" );
fprintf(stderr,"-H60        (default H60)    Halftone for shading \n" );
fprintf(stderr,"-K          (default gray-)  Color PostScript output \n" );
fprintf(stderr,"-KW         (default gray-)  Color PostScript output, whitened spectrum \n" );
fprintf(stderr,"-KR         (default gray-)  Color output Red->White->Blue\n");
fprintf(stderr,"-KB         (default gray-)  Color output Blue->White->Red\n");
fprintf(stderr,"-G          (default gray-)  Gray PostScript output \n" );
fprintf(stderr,"            (default is gray shading but all black colored lines) \n" );
fprintf(stderr,"-B          (default 8.5x11) Paper is 11 x 14 \n" );
fprintf(stderr,"-L          (default 8.5x11) Paper is 8.5x14 \n" );
fprintf(stderr,"-A3         (default 8.5x11) Paper is A3 \n" );
fprintf(stderr,"-A4         (default 8.5x11) Paper is A4 \n" );
fprintf(stderr,"-W          (default 0)      Line width in units of 0.001 in or 0.0025cm) \n" );
fprintf(stderr,"-EPS                         EPS output \n" );
fprintf(stderr,"-Ttitle     (default off)    Title at bottom left of page \n" );
fprintf(stderr,"-h                           Do not execute, show options \n" );
fprintf(stderr,"-?                           Do not execute, show options \n" );
	exit (0);
}

void dv_clip()
{
/*
fprintf(stderr,"dv_clip [%d %d %d %d] [%d %d %d %d]\n",
	dc_left,dc_bottom,dc_right,dc_top,
	dc_ClipLeft,dc_ClipBottom,dc_ClipRight,dc_ClipTop);
*/
		if(iseps == 1){
			if(dc_rotate){

				cxl = dc_top - dc_ClipTop;
				cyl = dc_ClipLeft;
				cxh = dc_top - dc_ClipBottom;
				cyh = dc_ClipRight;
			} else {
				cxl = dc_top - dc_ClipTop;
				cyl = dc_ClipLeft;
				cxh = dc_top - dc_ClipBottom;
				cyh = dc_ClipRight;
			}
		} else {
			if(dc_rotate){
				cxl = dc_top - dc_ClipBottom;
				cyl = dc_ClipLeft;
				cxh = dc_top - dc_ClipTop;
				cyh = dc_ClipRight;
			} else {
				cxl = dc_top - dc_ClipBottom;
				cyl = dc_ClipLeft;
				cxh = dc_top - dc_ClipTop;
				cyh = dc_ClipRight;
			}
		}
/*
		show_clip(cxl,cyl,cxh,cyh);
*/
}
void di_control(int type, int i1, int i2, int i3, int i4)
{
}
