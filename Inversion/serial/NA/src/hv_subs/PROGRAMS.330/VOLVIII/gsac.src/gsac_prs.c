/* CHANGES
 * 16 SEP 2004 - If O is specified plot travel time for the time axis
 * instead of time from the first sample
 *
 * 21 DEC 2004 - modified Landscape display so that the display of
 * positive amplitude is always in the direction ot the left of increasing 
 * time
 *
 * 27 JAN 2005 - eliminated any trace plot or trace annotation outside of the 
 		plot window
 * 15 JUN 2005 - corrected default menu option to actually reset all
 * 15 AUG 2005 - implemented interactive picking of P and S arrivals
 * 	TODO: make this work for all projections
 * 	if the output exists plot even in refr off mode
 * 08 FEB 2008 - corrected error for clipping when annotating
 * 14 AUG 2010 - added prs mag to the commands
 * 09 JUN 2011 - do not annotate trace if KSTNM is -12345
 * 23 JUL 2011 - add LAST kolor FIRST kolor to permit annotating a shaded trace
 * 25 MAR 2012 - minor change on axis title - also implement TITLE annotation
 * */
#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include	"csstim.h"

extern struct sacfile_ *sacdata;
extern void dogrid(void);
extern int *sortptr;

extern void gsac_exec_hold(void);

#define PRS_LS	-1
#define	PRS_PO 	-2
#define	PRS_SE 	-3
#define	PRS_RE 	-4
#define	PRS_ABS	-5
#define	PRS_REL	-6
#define	PRS_TIT	-7
#define PRS_VMIN -8
#define PRS_VMAX -9
#define PRS_DTDX -10
#define PRS_DTDD -11
#define PRS_TLIM -12
#define PRS_VLIM -13
#define PRS_DEF -14
#define PRS_AMP -15
#define PRS_SHD -16
#define PRS_COL -17
#define PRS_ANN -18
#define LANDSCAPE 0
#define PORTRAIT 1
#define SEASCAPE 2
#define REVERSE 3
#define ABSOLUTE 4
#define RELATIVE 5
#define PRS_SCALE_RELATIVE -19
#define PRS_SCALE_ABSOLUTE -20
#define PRS_SHOW_PICK -21
#define PRS_COL_FIRST -40
#define PRS_COL_LAST -41


static int prs_typeaxis = 0 ;	/* 0 LANDSCAPE, 1 PORTRAIT, 2 SEASCAPE */
static int prs_whichaxis = 50 ;	/* header variable.  */
static int prs_absolute = YES;
static char prs_str[1000];
static char prs_titstr[1000];
static char prs_annstr[100];
static float prs_amp = 0.5 ;		/* amplitude in inches */
static int prs_exttitle = NO; /* use external title */
static int prs_ann = NO; /* annotate trace with STA etc */
static int prs_doptaux = NO;	/* apply p-tau for dt/dx */
static int prs_doptaud = NO;	/* apply p-tau for dt/dd */
static float prs_dtdx = 0.0;
static float prs_dtdd = 0.0;
static float prs_tmin = -1.0e+30;
static float prs_tmax =  1.0e+30;
static float prs_vmin = 0;
static float prs_vmax = 0;
static int   dotlim   = NO;
static int   dovlim   = NO;
static int prs_doshd = NO;
static int prs_doshdplmn = 0;
static int prs_shdcolor = 1;
static int prs_shdcolor_first = -1;
static int prs_shdcolor_last  = -1;
static int prs_showpick = NO;
static float *py = (float *)NULL;
static float *px = (float *)NULL;
static int prs_scale_relative = YES;
static float prs_absolute_power;

static float prs_abs_max_amp;
static float prs_abs_min_amp;

FILE *prsctlfid ;


/* title from header value title[bhdr - 33] */
char *title[] = {
	"Receiver depth ",
	"",
	"",
	"",
	"",
	"Depth (km)",
	"Magnitude",
	"User0",
	"User1",
	"User2",
	"User3",
	"User4",
	"User5",
	"User6",
	"User7",
	"User8",
	"User9",
	"Distance (km)",
	"Azimuth (deg)",
	"Back Azimuth (deg)",
	"Distance (deg)",
};

struct arghdr prsarg[] = {
	{PRS_LS    , "LANDSCAPE", NHDR, 0, 0, NO, "Landscape", 1},
	{PRS_PO    , "PORTRAIT"	, NHDR, 0, 0, NO, "PORtrait ", 3},
	{PRS_SE    , "SEASCAPE"	, NHDR, 0, 0, NO, "SEascape ", 2},
	{PRS_SE    , "S"	, NHDR, 0, 0, NO, "SEASCAPE ", -1},
	{PRS_RE    , "REVERSE"	, NHDR, 0, 0, NO, "REVerse ",3},
	{PRS_REL   , "RELATIVE"	, NHDR, 0, 0, NO, "RELative ",3},
	{PRS_REL   , "R"	, NHDR, 0, 0, NO, "RELATIVE ",-1},
	{PRS_ABS   , "ABSOLUTE"	, NHDR, 0, 0, NO, "ABsolute ", 2},
	{PRS_ABS   , "A"	, NHDR, 0, 0, NO, "ABSOLUTE ",-1},
	{PRS_TIT   , "TITLE"	, CHDR, 0, 1, NO, "TItle string ", 2},
	{PRS_ANN   , "ANNOTATE"	, CHDR, 0, 1, NO, "ANnotate string ", 2},
	{PRS_DTDX  , "P"	, RHDR, 0, 1, NO, "DTDX p ",-1},
	{PRS_DTDX  , "PX"	, RHDR, 0, 1, NO, "DTDX px ",-1},
	{PRS_DTDX  , "DTDX"	, RHDR, 0, 1, NO, "DTDX dtdx ",-1},
	{PRS_DTDD  , "DTDD"	, RHDR, 0, 1, NO, "DTDD dtdd ",-1},
	{PRS_DTDD  , "PDEL"	, RHDR, 0, 1, NO, "PDel dtdd ", 2},
	{PRS_TLIM  , "TLIMIT"	, RHDR, 0, 2, NO, "TLimit tmin tmax", 2},
	{PRS_VLIM  , "VLIMIT"	, RHDR, 0, 2, NO, "VLimit vmin vmax", 2},
	{PRS_DEF   , "DEFAULT"	, NHDR, 0, 0, NO, "DEfault ", 2},
	{PRS_AMP   , "AMP"	, RHDR, 0, 1, NO, "AMp amp ", 2},
	{PRS_SHD   , "SHADE"	, CHDR, 0, 1, NO, "SHade [POS|NEG|ALL|OFF] ", 2},
	{PRS_COL   , "COLOR", IHDR, 0, 1, NO, "Color color ", 1},
	{PRS_COL_FIRST   , "KF", IHDR, 0, 1, NO, "KF first_trace_color ", -1},
	{PRS_COL_LAST    , "KL", IHDR, 0, 1, NO, "KL last_trace_color ", -1},
	{PRS_SCALE_RELATIVE  , "SCALERELATIVE", NHDR, 0, 0, NO, "ScaleRelative", -1},
	{PRS_SCALE_RELATIVE  , "SR", NHDR, 0, 0, NO, "ScaleRelative", -1},
	{PRS_SCALE_ABSOLUTE  , "SCALEABSOLUTE"	, RHDR, 0, 1, NO, "ScaleAbsolute power ",-1},
	{PRS_SCALE_ABSOLUTE  , "SA"	, RHDR, 0, 1, NO, "ScaleAbsolute power ",-1},
	{PRS_SHOW_PICK   , "PICK"	, YHDR, NO, 1, NO, "PIck [ON|OFF]  ", 2},
	{ H_EVDP,  "EVDP"    , RHDR, 0, 0, NO, "" ,-1},
	{ H_USER0, "USER0"   , RHDR, 0, 0, NO, "" ,-1},
	{ H_USER1, "USER1"   , RHDR, 0, 0, NO, "" ,-1},
	{ H_USER2, "USER2"   , RHDR, 0, 0, NO, "" ,-1},
	{ H_USER3, "USER3"   , RHDR, 0, 0, NO, "" ,-1},
	{ H_USER4, "USER4"   , RHDR, 0, 0, NO, "" ,-1},
	{ H_USER5, "USER5"   , RHDR, 0, 0, NO, "" ,-1},
	{ H_USER6, "USER6"   , RHDR, 0, 0, NO, "" ,-1},
	{ H_USER7, "USER7"   , RHDR, 0, 0, NO, "" ,-1},
	{ H_USER8, "USER8"   , RHDR, 0, 0, NO, "" ,-1},
	{ H_USER9, "USER9"   , RHDR, 0, 0, NO, "" ,-1},
	{ H_DIST,  "DIST"    , RHDR, 0, 0, NO, "" , 2},
	{ H_AZ,    "AZ"      , RHDR, 0, 0, NO, "" ,-1},
	{ H_BAZ,   "BAZ"     , RHDR, 0, 0, NO, "" ,-1},
	{ H_GCARC, "GCARC"   , RHDR, 0, 0, NO, "" , 1},
	{ H_STEL,  "STEL"    , RHDR, 0, 0, NO, "" ,-1},
	{ H_MAG ,  "MAG"     , RHDR, 0, 0, NO, "" ,-1},
	{  0, ""	, IHDR,  0, 0, NO, "",-1}
};

#define TITLE_LOC_TOP    0
#define TITLE_LOC_BOTTOM 1
#define TITLE_LOC_LEFT   2
#define TITLE_LOC_RIGHT  3

#define TITLE_SIZE_TINY   0
#define TITLE_SIZE_SMALL  1
#define TITLE_SIZE_MEDIUM 2
#define TITLE_SIZE_LARGE  3

extern int title_on ;
extern int title_size ;
extern int title_location;
extern char title_text[] ;

float tit_siz;


/* these are temporary variables only used here */
static float prs_real[10];
static int   prs_int [10];
static int   prs_yn;
static char  prs_tstr[100];
static int   prs_tstrl;
static char  *chdr_default = "-12345  ";

extern void XviG_Flush();
void prs_get_time_window(double *ptes,double *ptss,double *ptwin,int rs_absolute,float rs_tmin,float rs_tmax,int *prs_oset, int rs_doptaux, int rs_doptaud, float rs_dtdx, float rs_dtdd, int ntrc);
void gsac_prs_plot( int ntrc, float hvmax, float hvmin, double ts, double te, double twin, float x0 , float y0, float xlen, float ylen,float rs_amp, int rs_typeaxis, int rs_doshd , int rs_doshdplmn, int rs_ann, char *rs_annstr, int rs_whichaxis, int rs_absolute, double tss, float rs_dtdx, int rs_doptaux, float rs_dtdd, int rs_doptaud, int rs_doclip, float rs_cliplev, int rs_shdcolor, int rs_shdcolor_first, int rs_shdcolor_last);

void gsac_prs_plot_amp( int ntrc, float hvmax, float hvmin, double ts, double te, double twin, float x0 , float y0, float xlen, float ylen,float rs_amp, int rs_typeaxis, int rs_doshd , int rs_doshdplmn, int rs_ann, char *rs_annstr, int rs_whichaxis, int rs_absolute, double tss, float rs_dtdx, int rs_doptaux, float rs_dtdd, int rs_doptaud, int rs_doclip);


void gsac_set_param_prs(int ncmd, char **cmdstr)
{
	int i ;
int HasMouse; 
float XminDev, YminDev, 
	XmaxDev, YmaxDev, XminClip, 
	YminClip, XmaxClip, YmaxClip;
int Color;
	/*
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	*/
	if(gsac_control.prs > 0){
		if(gsac_control.prshist == NULL)
                        gsac_control.prshist = fopen("prshist.tmp","w+");
		for(i=0;i<ncmd;i++)
			fprintf(gsac_control.prshist,"%s ",cmdstr[i]);
		fprintf(gsac_control.prshist,"\n");
		fflush(gsac_control.prshist);
	}
	/* initialize graphics */
	if(gsac_control.plotinit == NO){
		if(gsac_control.plotdevice==WIN){
			ginitf("INTEM","GSAC");
			printf("Initializing Interactive Graphics\n");
			gmesg("Initializing Interactive Graphics");
			gsac_control.everinteractive = YES;
			gsac_control.plotinit = YES;
			gsac_control.plotchange = NO;
			ginfo(&HasMouse, &XminDev, &YminDev, 
				&XmaxDev, &YmaxDev, &XminClip, 
				&YminClip, &XmaxClip, &YmaxClip,&Color);
			gsac_control.XmaxDev = XmaxDev;
			gsac_control.YmaxDev = YmaxDev;
			if(Color >= 4)
				gsac_control.black = 0;
			else
				gsac_control.black = 1;
			gsac_control.kolor = Color%4;

		}
	}
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, prsarg, NO, YES))
	       return	;
	for(i=0 ; prsarg[i].key[0] != '\0' ; i++){
		if(prsarg[i].used > 0){
			if(prsarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, prsarg[i].key, 
					prsarg[i].mfit,prsarg[i].narg, prs_real);
			} else if(prsarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, prsarg[i].key, 
					prsarg[i].mfit,prsarg[i].narg, prs_int );
			} else if(prsarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, prsarg[i].key, 
					prsarg[i].mfit,prsarg[i].narg, prs_str );
			}
			
			if(prsarg[i].id >=  0 )	{
				/* get what item to sort */
				if(prs_whichaxis != prsarg[i].id){
					prs_whichaxis = prsarg[i].id;
					/* 25 MAR 2012
					prs_exttitle = NO;
					*/
				}
			} else {
				switch(prsarg[i].id){
					/* for this to work we must 
					 * */
					case PRS_DEF:
						prs_amp = 0.5 ;
						prs_exttitle = NO;
						prs_ann = NO;
						prs_doptaux = NO;
						prs_doptaud = NO;
						prs_dtdx = 0.0;
						prs_dtdd = 0.0;
						prs_tmin = -1.0e+30;
						prs_tmax =  1.0e+30;
						prs_vmin = 0;
						prs_vmax = 0;
						dotlim   = NO;
						dovlim   = NO;
						prs_doshd = NO;
						prs_doshdplmn = 0;
						prs_shdcolor = 1;
						prs_shdcolor_first = -1;
						prs_shdcolor_last = -1;
						prs_absolute = YES ;
						prs_scale_relative=YES;
						prs_exttitle = NO;
						break;
					case PRS_LS:
						prs_typeaxis = LANDSCAPE ;
						break;
					case PRS_PO:
						prs_typeaxis = PORTRAIT ;
						break;
					case PRS_SE:
						prs_typeaxis = SEASCAPE ;
						break;
					case PRS_RE:
						prs_typeaxis = REVERSE ;
						break;
					case PRS_REL:
						prs_absolute = NO ;
						break;
					case PRS_ABS:
						prs_absolute = YES ;
						break;
					case PRS_TIT:
						strcpy(prs_titstr,prs_str);
						prs_exttitle = YES;
						break;
					case PRS_ANN:
						gsac_strupr(prs_str);
						if(strcmp(prs_str,"OFF")==0){
							prs_ann = NO;
						} else {
							strcpy(prs_annstr,prs_str);
							prs_ann = YES;
						}
						break;
					case PRS_DTDX:
						prs_doptaux = YES;
						prs_doptaud = NO ;
						prs_dtdx = prs_real[0];
						break;
					case PRS_DTDD:
						prs_doptaud = YES;
						prs_doptaux = NO ;
						prs_dtdd = prs_real[0];
						break;
					case PRS_TLIM:
						dotlim = YES;
						prs_tmin = prs_real[0];
						prs_tmax = prs_real[1];
						if(prs_tmax == prs_tmin){
							printf("Time limits must differ\n");
							dotlim = NO;
						}
						break;
					case PRS_VLIM:
						dovlim = YES;
						prs_vmin = prs_real[0];
						prs_vmax = prs_real[1];
						if(prs_vmax == prs_vmin){
							printf("Axis limits must differ\n");
							dovlim = NO;
						}
						break;
					case PRS_AMP:
						prs_amp = prs_real[0] ;
						break;
					case PRS_COL:
						prs_shdcolor = prs_int[0] ;
						break;
					case PRS_COL_LAST:
						prs_shdcolor_last = prs_int[0] ;
						break;
					case PRS_COL_FIRST:
						prs_shdcolor_first = prs_int[0] ;
						break;
					case PRS_SHOW_PICK:
                                                prs_showpick = prs_yn;
						break;
					case PRS_SHD:
						gsac_strupr(prs_str);
						if(strncmp(prs_str,"P",1)==0){
							prs_doshdplmn = 0;
							prs_doshd = YES;
						} else if(strncmp(prs_str,"N",1)==0){
							prs_doshdplmn = 1;
							prs_doshd = YES;
						} else if(strncmp(prs_str,"A",1)==0){
							prs_doshdplmn = 2;
							prs_doshd = YES;
						} else if(strncmp(prs_str,"O",1)==0){
							prs_doshd = NO;
						}
						break;
					case PRS_SCALE_RELATIVE:
						prs_scale_relative = YES ;
						break;
					case PRS_SCALE_ABSOLUTE:
						/* enforce the power choice */
						if(prs_real[0] >= 0.0 
							&& prs_real[0] <= 2.0){

							prs_scale_relative=NO;
						prs_absolute_power=prs_real[0];
						} else {
							prs_scale_relative=YES;
						}
						break;

				}
			}
		}
	}
}

void  gsac_prs(float rs_dtdx,float rs_dtdd,int rs_absolute,int rs_doptaux,
	int rs_doptaud,float rs_tmin,float rs_tmax,float rs_vmin,float rs_vmax,
	float x0, float y0, float xlen, float ylen, float rs_amp, int rs_typeaxis,
	int rs_doshd,int rs_doshdplmn,int rs_ann,char *rs_annstr,int rs_exttitle,
	char *rs_titstr, char  *rs_tstr, int rs_tstrl, int rs_whichaxis, int rs_shdcolor,
	int rs_shdcolor_first, int rs_shdcolor_last,
	double *tss,double *tes,float *hvmin,float *hvmax,int dotlim,int dovlim,int labx,
	int rs_doclip,float rs_cliplev);

void gsac_exec_prs(void)
{
	int ntrc;
	char pltname[12];
	char fname[12];
	float x0, y0, xlen, ylen;
	double tss,tes;
	float hvmin, hvmax;


	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	gsac_exec_hold();
	if(ntrc < 1)
		return;
	/* initialize */
	if(gsac_control.plotdevice==WIN){
		if(gsac_control.hold == NO){
			gframe(2);
		} else if(gsac_control.hold == 1){
			gframe(2);
			gsac_control.hold++ ;
		}
	}
	if(gsac_control.plotdevice==PLT){
		if(gsac_control.hold == NO || 
			(gsac_control.hold != NO && 
			 gsac_control.inpltmode == NO)){
		sprintf(pltname,"PRS%3.3d.PLT",gsac_control.plotcount_prs);
		printf("Initializing %s\n",pltname);
		ginitf(pltname,"GSAC");
		gsac_control.inpltmode = YES ;


		}
	} else {
		if(gsac_control.hold == NO){
			gframe(2);
		}
	}
	xlen = gsac_control.xlen ;
	ylen = gsac_control.ylen ;
	x0   = gsac_control.x0 ;
	y0   = gsac_control.y0 ;

	gsac_prs(prs_dtdx,prs_dtdd,prs_absolute,prs_doptaux,prs_doptaud,
		prs_tmin, prs_tmax, prs_vmin, prs_vmax, x0, y0, xlen, ylen, prs_amp,
		prs_typeaxis,prs_doshd,prs_doshdplmn,prs_ann,prs_annstr,
		prs_exttitle, prs_titstr, prs_tstr, prs_tstrl, prs_whichaxis, prs_shdcolor,
		prs_shdcolor_first, prs_shdcolor_last,
		&tss,&tes,&hvmin,&hvmax,dotlim,dovlim,YES,NO,1.0);
	/* output the CTL file but only after the plot has been made which defines boundaries */
	if(gsac_control.plotdevice==PLT){
		if(gsac_control.hold == NO || 
			(gsac_control.hold != NO && 
			 gsac_control.inpltmode == NO)){
		/* create PRS%3.3d.CTL  eventually check for */
		sprintf(fname,"PRS%3.3d.CTL",gsac_control.plotcount_prs );
		prsctlfid=fopen(fname,"w+");
		fprintf(prsctlfid,
			"refmod96 -XMIN %f -XMAX %f -TMIN %f -TMAX %f -XLEN %f -YLEN %f -X0 %f -Y0 %f -KF 2 -KR 4 -NO -NM 1 -P -SH -M model\n",
			hvmin, hvmax, tss, tes,
			gsac_control.xlen, gsac_control.ylen, 
			gsac_control.x0, gsac_control.y0
		);
		/* close PRS%3.3d.CTL" */
		fclose(prsctlfid);
		/* now update the counter since we need this for the CTL file */
		gsac_control.plotcount_prs++;
		}
	}
}

void  gsac_prs(float rs_dtdx,float rs_dtdd,int rs_absolute,int rs_doptaux,
	int rs_doptaud,float rs_tmin,float rs_tmax,float rs_vmin,float rs_vmax,
	float x0, float y0, float xlen, float ylen, float rs_amp, int rs_typeaxis,
	int rs_doshd,int rs_doshdplmn,int rs_ann,char *rs_annstr,int rs_exttitle,
	char *rs_titstr, char  *rs_tstr, int rs_tstrl, int rs_whichaxis, int rs_shdcolor,
	int rs_shdcolor_first, int rs_shdcolor_last,
	double *ttss,double *ttes,float *hhvmin,float *hhvmax,int dotlim,int dovlim,int labx,
	int rs_doclip,float rs_cliplev){
	/* float rs_dtdx - ray parameter sec/km
	 * float rs_dtdd - ray parameter sec/deg
	 * int rs_absolute 
	 * int rs_doptaux - plot t - px
	 * int rs_doptaud - plot t - p DELTA
	 * float rs_tmin - minimum time
	 * float rs_tmax - maximum time
	 * float rs_vmin - minimum value, dypically distance
	 * float rs_vmax - maximum value, dypically distance
	 * float x0      - lower left corner
	 * float y0      - lower left corner
	 * float xlen    - horizontal axis length
	 * float ylen    - vertical axis length
	 * float rs_amp  - amplitude factor
	 * int rs_typeaxis - PORTRAIT, LANDSCAPE, SEASCAPE, REVERSE
	 * int rs_doshd    - do shaded area plot
	 * int rs_doshdplmn - 0 pos,1 neg,2 all
	 * int rs_ann       - annotate with
	 * char *rs_annstr  -  this string, e.g., STA or NAME (file)
	 * int rs_exttitle  - YES NO
	 * char *rs_titstr  - the title
	 * char  *rs_tstr   - time axis string, e.g., T - p D
	 * int rs_tstrl	    - length of time axis string
	 * int rs_whichaxis - header variable for axis, e.g., 
	 * 		DIST GCARC, ray parameter
	 * int rs_shdcolor  - CALPLOT color for shading
	 * int rs_shdcolor_first  - CALPLOT color for shading for first trace read 
	 * int rs_shdcolor_last   - CALPLOT color for shading for last trace read
	 * double ttss       - starting time
	 * double ttes	    - ending time
	 * float hhvmin      - minimum axis value
	 * float hhvmax      - maximum axis value
	 * int labx          - if YES label Dist (e.g.) axis and annotate with Origin, Ray parameter
	 * int rs_doclip     - YES clip trace according to xlen/ntrc multiple 
	 * float rs_cliplev  - default 2
	 * */

	int ntrc;
	float dy;
	float hvmin, hvmax;
	float v ;
	int k, ls ;
	double ts,te,twin;
	double tss,tes;
	float dclip;
	float av;	/* plot amplitude of trace in page units */

	char timstr[32];
	int rs_oset ;
	char ostr[100];
	char lstr[100];
	float vlen ;  /* length of vaxis , e.g., not time axis */

	ntrc = gsac_control.number_itraces;

	prs_get_time_window(&tes,&tss,&twin,rs_absolute,rs_tmin,rs_tmax,
		&rs_oset, rs_doptaux, rs_doptaud, rs_dtdx, rs_dtdd, ntrc);

	*ttss = tss;
	*ttes = tes;

	/* we must define the time window for the plot */


	dy = ylen ;
	if(rs_typeaxis == LANDSCAPE || rs_typeaxis == SEASCAPE){
		vlen = xlen;
		dclip = xlen/(float)ntrc;
	} else {
		vlen = ylen;
		dclip = ylen/(float)ntrc;
	}
	rs_tstrl = strlen(rs_tstr);
	gclip("off", x0, y0, x0+xlen, y0+dy);
	if(gsac_control.hold == NO)gframe(1);
	/* determine the bounds on the real header variable in the data set */
	if(dovlim == YES){
		hvmin = rs_vmin;
		hvmax = rs_vmax;
	} else {
		hvmin =  1.0e+37;
		hvmax = -1.0e+37;
		for ( k=0 ; k < ntrc  ; k++){
			v = sacdata[k].sachdr.rhdr[rs_whichaxis];
			if(hvmin > v ) hvmin = v;
			if(hvmax < v ) hvmax = v;
		}
		/* since not exact clipping, adjust for the trace amplitude */
		av = hvmax - hvmin;
		if(rs_doclip){
			av = SIGN(av)*MIN(rs_cliplev*dclip,ABS(av));
		}
		hvmin -= av*rs_amp/vlen;
		hvmax += av*rs_amp/vlen;
	}
	*hhvmin = hvmin;
	*hhvmax = hvmax;
	/* safety */
	if(hvmax == hvmin){
		if(hvmax == 0.0){
			hvmax = 1.0;
		}
		hvmin = - hvmax ; /* or perhaps 0 -> 1 instead of -1 -> 1 */
	}

	/* to make a nice plot we will adjust the hvmin and hvmax
	 * to accomodate the peak amplitude 
	 * normally hvmax - hvmin -> axis length (inches), 
	 * */
	/* now begin to plot . We initially do things dimensionlessly
	 * and then actually imnplement the PORTRAIT, LANDSCAPE and SEASCAPE */
	if(gsac_control.grid)
		dogrid();

	if(gsac_control.background == YES &&
		gsac_control.background_color >= 0){
		newpen(gsac_control.background_color);
		shader(x0,y0,x0+xlen,y0+dy,0,0,0.01,0.01);
		newpen(1);
	}
	gbox(x0,y0,x0+xlen,y0+ylen);
	/* put up the axes - note that this relies on the specific
	 * header identification */
	if(rs_exttitle == NO)
		strcpy(lstr, title[rs_whichaxis-33]);
	else
		strcpy(lstr, rs_titstr);
	ls = strlen(lstr);
	if(rs_doptaux == YES){
		strcpy(rs_tstr,"T -p X (s)");
	} else if(rs_doptaud == YES){
		strcpy(rs_tstr,"T -p D (s)");
	} else {
		strcpy(rs_tstr,"Time (s)");
	}
	rs_tstrl = strlen(rs_tstr);
	if(rs_typeaxis == LANDSCAPE){
		if(gsac_control.ygrid == YES) dolnygrid(x0,x0+xlen,y0,ylen,
			tss,tes,0.10, YES, gsac_control.ygrid_color, 
			gsac_control.ygrid_type,gsac_control.ygrid_minor);
		doliny(x0,y0,ylen,tss,tes,0.10,YES,YES,YES,rs_tstrl,rs_tstr);
		doliny(x0+xlen,y0,ylen,tss,tes,0.10,NO,NO,NO,rs_tstrl,rs_tstr);
		if(gsac_control.xgrid == YES)
		dolnxgrid(x0,y0,y0+ylen,xlen,hvmin,hvmax,0.10,
			YES, gsac_control.xgrid_color, 
			gsac_control.xgrid_type,gsac_control.xgrid_minor);
		dolinx(x0,y0,xlen,hvmin,hvmax,0.10,NO,NO,YES,ls,lstr);
		dolinx(x0,y0+ylen,xlen,hvmin,hvmax,0.10,YES,NO,NO,ls,lstr);
	} else if(rs_typeaxis == SEASCAPE){
		if(gsac_control.ygrid == YES) dolnygrid(x0,x0+xlen,y0,ylen,
			-tss,-tes,0.10, NO, gsac_control.ygrid_color, 
			gsac_control.ygrid_type,gsac_control.ygrid_minor);
		dnliny(x0,y0,ylen,-tss,-tes,0.10,YES,YES,YES,rs_tstrl,rs_tstr);
		dnliny(x0+xlen,y0,ylen,-tss,-tes,0.10, NO, NO, NO, rs_tstrl, rs_tstr);
		if(gsac_control.xgrid == YES)
		dolnxgrid(x0,y0,y0+ylen,xlen,hvmin,hvmax,0.10,
			YES, gsac_control.xgrid_color, 
			gsac_control.xgrid_type,gsac_control.xgrid_minor);
		dolinx(x0,y0,xlen,hvmin,hvmax,0.10,NO,NO,NO,ls,lstr);
		dolinx(x0,y0+ylen,xlen,hvmin,hvmax,0.10,YES,YES,YES,ls,lstr);
	} else if(rs_typeaxis == PORTRAIT){
		if(gsac_control.xgrid == YES)
		dolnxgrid(x0,y0,y0+ylen,xlen,tss,tes,0.10,
			YES, gsac_control.xgrid_color, 
			gsac_control.xgrid_type,gsac_control.xgrid_minor);
		dolinx(x0,y0,xlen,tss,tes,0.10,NO,NO,YES,rs_tstrl,rs_tstr);
		dolinx(x0,y0+ylen,xlen,tss,tes,0.10,YES,NO,NO,rs_tstrl,rs_tstr);
		if(gsac_control.ygrid == YES) dolnygrid(x0,x0+xlen,y0,ylen,
			hvmin,hvmax,0.10, YES, gsac_control.ygrid_color, 
			gsac_control.ygrid_type,gsac_control.ygrid_minor);
		doliny(x0,y0,ylen,hvmin,hvmax,0.10,YES,YES,YES,ls,lstr);
		doliny(x0+xlen,y0,ylen,hvmin,hvmax,0.10,NO,NO,NO,ls,lstr);
	} else if(rs_typeaxis == REVERSE){
		if(gsac_control.xgrid == YES)
		dolnxgrid(x0,y0,y0+ylen,xlen,tss,tes,0.10,
			YES, gsac_control.xgrid_color, 
			gsac_control.xgrid_type,gsac_control.xgrid_minor);
		dolinx(x0,y0,xlen,tss,tes,0.10,YES,NO,NO,rs_tstrl,rs_tstr);
		dolinx(x0,y0+ylen,xlen,tss,tes,0.10,NO,YES,YES,rs_tstrl,rs_tstr);
		if(gsac_control.ygrid == YES) dolnygrid(x0,x0+xlen,y0,ylen,
			-hvmin,-hvmax,0.10, NO, gsac_control.ygrid_color, 
			gsac_control.ygrid_type,gsac_control.ygrid_minor);
		dnliny(x0,y0,ylen,-hvmin,-hvmax,0.10, YES, YES, YES, ls, lstr);
		dnliny(x0+xlen,y0,ylen,-hvmin,-hvmax,0.10, NO, NO, NO, ls, lstr);
	}
		if(gsac_control.plotdevice==WIN)
			XviG_Flush();
		/* annotate with the plot titles
			added 25 MAR 2012 */
		if(title_on == YES){
			switch(title_size){
				case TITLE_SIZE_TINY   :
					tit_siz = 0.05 ;
					break;
				case TITLE_SIZE_SMALL  :
					tit_siz = 0.1 ;
					break;
				case TITLE_SIZE_MEDIUM :
					tit_siz = 0.15 ;
					break;
				case TITLE_SIZE_LARGE  :
					tit_siz = 0.2 ;
					break;
			}
			if(title_location == TITLE_LOC_TOP){
				gcent(x0+0.5*xlen,y0+ylen+0.2,tit_siz,title_text,0.0);
			} else if(title_location == TITLE_LOC_BOTTOM){
				gcent(x0+0.5*xlen,y0-0.7,tit_siz,title_text,0.0);
			} else if(title_location == TITLE_LOC_LEFT){
				gcent(x0-0.8,y0+0.5*ylen,tit_siz,title_text, 90.0);
			} else if(title_location == TITLE_LOC_RIGHT){
				gcent(x0+xlen+0.2,y0+0.5*ylen,tit_siz,title_text,-90.0);
			}

		}
	/* annotate with reference time if absolute and O not set */
	if(labx){
		if(rs_absolute){
			if(rs_oset == YES){
				gleft(x0,y0 - 0.1*ylen,0.01*ylen,"Origin set",0.0);
			} else {
				printtimestr(gsac_control.begmin,timstr);
				gleft(x0,y0 - 0.1*ylen,0.01*ylen,timstr,0.0);
			}
		} else {
			gleft(x0,y0 - 0.1*ylen,0.01*ylen,"Relative",0.0);
		}
		if(rs_doptaux){
			sprintf(ostr,"p = %5.3f s/km",rs_dtdx);
			gleft(x0,y0 - 0.12*ylen,0.01*ylen,ostr,0.0);
		}
		if(rs_doptaud){
			sprintf(ostr,"p = %5.3f s/deg",rs_dtdd);
			gleft(x0,y0 - 0.12*ylen,0.01*ylen,ostr,0.0);
		}
	}
	/* initial cut at this change later */
		gclip("on", x0, y0, x0+xlen, y0+dy);

	/* define the maximum and minimum amplitude of each trace
		within the plot window */

	
	gsac_prs_plot_amp(  ntrc,  hvmax,  hvmin,  ts,  te,  twin,  x0 ,  y0,  xlen,  ylen, rs_amp,  rs_typeaxis,  rs_doshd ,  rs_doshdplmn,  rs_ann, rs_annstr, rs_whichaxis, rs_absolute, tss,  rs_dtdx,  rs_doptaux,  rs_dtdd, rs_doptaud,  rs_doclip);

	/* now do the plot */

	gsac_prs_plot(  ntrc,  hvmax,  hvmin,  ts,  te,  twin,  x0 ,  y0,  xlen,  ylen, rs_amp,  rs_typeaxis,  rs_doshd ,  rs_doshdplmn,  rs_ann, rs_annstr, rs_whichaxis, rs_absolute, tss,  rs_dtdx,  rs_doptaux,  rs_dtdd, rs_doptaud,  rs_doclip,rs_cliplev, rs_shdcolor,  rs_shdcolor_first,  rs_shdcolor_last);
		gclip("off", x0, y0, x0+xlen, y0+dy);
	/* flush the buffer if in interactive mode */
	if(gsac_control.plotdevice==WIN)
		XviG_Flush();
	/* clean up */
	if(gsac_control.plotdevice == PLT && gsac_control.hold == NO){
			/* force new Pnnnn.PLT on next call */
		gsac_control.plotchange = NO; 
		gsac_control.inpltmode = NO;
		gend(0);
	}
	if(gsac_control.hold == NO){
		if(gsac_control.everinteractive == YES){
			ginitf("INTEM","GSAC");
		}
	}

}
void prs_get_time_window(double *ptes,double *ptss,double *ptwin,int rs_absolute,float rs_tmin,float rs_tmax,int *prs_oset,int rs_doptaux, int rs_doptaud, float rs_dtdx, float rs_dtdd, int ntrc)
{
	int k, kk ;
	double otss, otes;
	double ts, te;
	double tss, tes, twin;
	int rs_oset;
	ntrc = gsac_control.number_itraces;
	/* define the timing window for the plot */
	rs_oset = NO;
	if(rs_absolute == YES) {
		/* instead of using preset begin and end we must recalcuate
		 * because of the possibility of a p-tau plot - note
		 * we only do a p-tau for traces which have a defined distance 
		 * and we only plot those traces */
		/* if O is set define tss and tes */
		otss =  1.0e+30;
		otes = -1.0e+30;
		twin = 0.0;
		for ( k=0 ; k < ntrc ; k++){
			if(sacdata[k].sachdr.rhdr[H_O] != -12345.){
				tss =   sacdata[k].tzbeg
					- ( sacdata[k].tzref + 
					sacdata[k].sachdr.rhdr[H_O]);
				tes =   tss +
					sacdata[k].sachdr.rhdr[H_E] -
					sacdata[k].sachdr.rhdr[H_B];
				/* make correction for p-tau */
				if(rs_doptaux){
				tss = tss - rs_dtdx*sacdata[k].sachdr.rhdr[H_DIST] ;
				tes = tes - rs_dtdx*sacdata[k].sachdr.rhdr[H_DIST] ;
				} else if(rs_doptaud){
				tss = tss - rs_dtdd*sacdata[k].sachdr.rhdr[H_GCARC] ;
				tes = tes - rs_dtdd*sacdata[k].sachdr.rhdr[H_GCARC] ;
				}
				tss = MIN(otss,tss);
				otss = tss;
				tes = MAX(otes,tes);
				otes = tes;
				twin = MAX(twin,tes - tss);
			}
		}
		/* now that we have examined all the traces consider case the
		 * O is not set */
		if( twin == 0.0){
			twin = (gsac_control.endmax - gsac_control.begmin);
			tss = 0.0;
			tes = twin;
		} else {
			rs_oset = YES;
		}
	} else {
		twin = 0.0;
		tss = 0.0;
		for ( kk=0 ; kk < ntrc  ; kk++){
			k = sortptr[kk];
			ts = sacdata[k].sachdr.rhdr[H_B];
			te = sacdata[k].sachdr.rhdr[H_E];
			if((te-ts) > twin )
				twin = te - ts ;
			/* note that the LAST instance of ts te is preserved
			 * for the time scale */
			if(kk==0)
				tss = ts;
			tes = tss + twin;

		}

	}
	if(dotlim==YES ){
		tss = MIN(rs_tmin,rs_tmax);
		tes = MAX(rs_tmin,rs_tmax);
		twin = tes - tss;
	}
	*prs_oset = rs_oset;
	*ptss = tss;
	*ptes = tes;
	*ptwin = twin;
}

void gsac_prs_plot( int ntrc, float hvmax, float hvmin, double ts, double te, double twin, float x0 , float y0, float xlen, float ylen,float rs_amp, int rs_typeaxis, int rs_doshd , int rs_doshdplmn, int rs_ann, char *rs_annstr, int rs_whichaxis, int rs_absolute, double tss, float rs_dtdx, int rs_doptaux, float rs_dtdd, int rs_doptaud, int rs_doclip, float rs_cliplev, int rs_shdcolor, int rs_shdcolor_first, int rs_shdcolor_last)
{
	int k, i, ls, kk, nii;
	float v, uu, vv, xx, yy;
	float dy, dv;
	float trmax;    /* trace maximum amplitude */
	float tv;       /* normalized trace value [-1,1] */
	float tb, to;
	float dclip;
	float av;       /* plot amplitude of trace in page units */
	char rs_strann[100];
	int npts;
	float tamp;

	/* safety feature for shading */
	if(rs_doshd == YES){
		if(rs_shdcolor_last < 0) rs_shdcolor_last = rs_shdcolor;
		if(rs_shdcolor_first < 0) rs_shdcolor_first = rs_shdcolor;
	}
	
	for ( kk=0 ; kk < ntrc  ; kk++){
		k = sortptr[kk];
		if(prs_scale_relative == YES){
			/* use the extreme of the trace within the
				plot window */
			trmax = MAX(ABS(sacdata[k].winmax),
				ABS(sacdata[k].winmin));
		} else {
			trmax = MAX(ABS(prs_abs_max_amp),
				ABS(prs_abs_min_amp));
		}
		/* safety */
		if(trmax == 0.0)
			trmax = 1.0;
		v = sacdata[k].sachdr.rhdr[rs_whichaxis] ;
		tamp = rs_amp;
		if(prs_scale_relative == NO){
			if(prs_absolute_power == 0.0)
				tamp *= 1.0;
			if(prs_absolute_power == 1.0)
				tamp *= v;
			else
				tamp *= pow(v,prs_absolute_power);
		}
		uu = (v - hvmin)/(hvmax - hvmin);
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		tb = sacdata[k].sachdr.rhdr[H_B];
		te = sacdata[k].sachdr.rhdr[H_E];
		to = sacdata[k].sachdr.rhdr[H_O];
		dv = sacdata[k].sachdr.rhdr[H_DELTA]/twin;
		if(kk == 0)
			ts = tb;
		if(rs_absolute == YES) {
			vv = (float)(sacdata[k].tzref + tb - gsac_control.begmin -tss)/twin;
			if(to != -12345.){
				vv = (float)(tb - to - tss)/twin ;
			}
			if(rs_doptaux){
				vv = vv - rs_dtdx*sacdata[k].sachdr.rhdr[H_DIST]/twin ;
			} else if(rs_doptaud){
				vv = vv - rs_dtdd*sacdata[k].sachdr.rhdr[H_GCARC]/twin ;
			}
		} else {
/* changed 18 MAR 2006 was tb - tss and did not do relative */
			vv = (ts-tss)/twin;
			if(rs_doptaux){
				vv = vv - rs_dtdx*sacdata[k].sachdr.rhdr[H_DIST]/twin ;
			} else if(rs_doptaud){
				vv = vv - rs_dtdd*sacdata[k].sachdr.rhdr[H_GCARC]/twin ;
			}
		}
		/* allocate array for the trace */
		if(px == (float *)NULL)
			px = (float *)calloc(npts,sizeof(float));
		else
			px = (float *)realloc(px,npts*sizeof(float));
		if(py == (float *)NULL)
			py = (float *)calloc(npts,sizeof(float));
		else
			py = (float *)realloc(py,npts*sizeof(float));
		/* nii is a counter of the actual points plotted */
		for(i = 0 , nii = 0; i < npts ; i++){
			tv = sacdata[k].sac_data[i]/trmax;
			av = tamp*tv;
			if(rs_doclip){
				av = SIGN(av)*MIN(rs_cliplev*dclip,ABS(av));
			}
			if(rs_typeaxis == LANDSCAPE){
				xx = x0 + uu*xlen - av;
				yy = y0 + vv*ylen;
			} else if(rs_typeaxis == SEASCAPE){
				xx = x0 + uu*xlen + av;
				yy = y0 + ylen - vv*ylen ;
			} else if(rs_typeaxis == PORTRAIT){
				xx = x0 + vv*xlen;
				yy = y0 + uu*ylen + av;
			} else if(rs_typeaxis == REVERSE){
				xx = x0 + vv*xlen;
				yy = y0 + ylen - uu*ylen + av;
			}
			if(uu >= 0.0 && uu <= 1.0 && vv >= 0.0 && vv <= 1.0){
				px[nii] = xx ;
				py[nii] = yy ;
				nii++;
			}
			vv += dv;
		}
		/* if required, plot shaded trace first */
				
		if(rs_doshd == YES){
			if(rs_typeaxis == LANDSCAPE){
				if(rs_doshdplmn == 0)
					shdsei(x0+uu*xlen,y0,YES ,YES,1);
				else if(rs_doshdplmn == 1)
					shdsei(x0+uu*xlen,y0,YES ,YES,0);
				else if(rs_doshdplmn == 2)
					shdsei(x0+uu*xlen,y0,YES ,YES,2);
			} else if(rs_typeaxis == SEASCAPE){
				shdsei(x0+uu*xlen,y0,YES ,YES,rs_doshdplmn);
			} else if(rs_typeaxis == PORTRAIT){
				shdsei(x0,y0+uu*ylen,NO ,YES,rs_doshdplmn);
			} else if(rs_typeaxis == REVERSE){
				shdsei(x0,y0+(1.0-uu)*ylen,NO ,YES,rs_doshdplmn);
			}
			if(kk == 0)
				newpen(rs_shdcolor_first);
			else if ( kk == ntrc -1 )
				newpen(rs_shdcolor_last);
			else
				newpen(rs_shdcolor);
			if(nii > 0){
				plot(px[0],py[0],3);
/* never put out a very long trace - just use 500 sample chunks 
 	THIS IS A NICE STATEMENT BUT NOT IMPLEMENTED */
				for(i = 0 ; i < nii ; i++){
					plot(px[i],py[i],2);
				}
				plot(px[nii-1],py[nii-1],3);
			}
			shdsei(x0,y0,NO,NO,0);
		}
		
		gsac_setcolor(YES, k, ntrc);
		if(nii > 0){
			plot(px[0],py[0],3);
			for(i = 0 ; i < nii ; i++){
				plot(px[i],py[i],2);
			}
			plot(px[nii-1],py[nii-1],3);
		}
		gsac_setcolor(NO , k, ntrc);
		/* annotate the trace with picks */
		/* annotate the trace if desired, but only if the trace 
			is plotted  */
		if(rs_ann == YES && uu >= 0.0 && uu <= 1.0){
			/* we will put this outside the axis - we currently
			 * do this only if the annotation is STA we can later
			 * permit date  or filename - we use rs_strann
			 * so that we can later use other things here
			 * */
			if(strncmp(rs_annstr,"STA",3)==0){
				/* annotate */
                                if(strncmp(sacdata[k].schdr[H_KSTNM],chdr_default,8)!=0)
				   strcpy(rs_strann,sacdata[k].schdr[H_KSTNM]);
                                else
				   strcpy(rs_strann,"");
			} else if(strncmp(rs_annstr,"NAME",4)==0){
				/* annotate */
				ls = strlen(sacdata[k].sac_ifile_name);
				if(ls > sizeof(rs_strann)-1)ls--;
				strncpy(rs_strann,sacdata[k].sac_ifile_name,ls);
				rs_strann[ls]='\0';
			}
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			if(rs_typeaxis == LANDSCAPE){
				xx = x0 + uu*xlen ;
				yy = y0 + ylen + 0.2;
				gleft(xx,yy,0.1,rs_strann,90.0);
			} else if(rs_typeaxis == SEASCAPE){
				xx = x0 + uu*xlen ;
				yy = y0 + ylen - ylen -0.2 ;
				gright(xx,yy,0.1,rs_strann,90.0);
			} else if(rs_typeaxis == PORTRAIT){
				xx = x0 + xlen + 0.2;
				yy = y0 + uu*ylen ;
				gleft(xx,yy,0.1,rs_strann,0.0);
			} else if(rs_typeaxis == REVERSE){
				xx = x0 + xlen + 0.2;
				yy = y0 + ylen - uu*ylen ;
				gleft(xx,yy,0.1,rs_strann,0.0);
			}
			/* now determine where to put the string */
			gmesg(" ");
			gclip("on", x0, y0, x0+xlen, y0+ylen);
		}
	}
}

/* this duplicates much of gsac_prs_plot but is separate for clarity of code
   its only purpose to to examine the trace within the plot window and
   determine the trace amplitude extremes
*/
void gsac_prs_plot_amp( int ntrc, float hvmax, float hvmin, double ts, double te, double twin, float x0 , float y0, float xlen, float ylen,float rs_amp, int rs_typeaxis, int rs_doshd , int rs_doshdplmn, int rs_ann, char *rs_annstr, int rs_whichaxis, int rs_absolute, double tss, float rs_dtdx, int rs_doptaux, float rs_dtdd, int rs_doptaud, int rs_doclip)
{
	int k, i, kk, nii;
	float v, uu, vv;
	float dv;
	float trmax;    /* trace maximum amplitude */
	float tv;       /* normalized trace value [-1,1] */
	float tb, to;
	int npts;
	
	prs_abs_max_amp = -1.0e+38 ;
	prs_abs_min_amp =  1.0e+38 ;
	for ( kk=0 ; kk < ntrc  ; kk++){
		k = sortptr[kk];
		sacdata[k].winmax = -1.0e+38;
		sacdata[k].winmin =  1.0e+38;
		trmax = MAX(ABS(sacdata[k].sachdr.rhdr[H_DEPMIN]),
			ABS(sacdata[k].sachdr.rhdr[H_DEPMAX]));
		/* safety */
		if(trmax == 0.0)
			trmax = 1.0;
		v = sacdata[k].sachdr.rhdr[rs_whichaxis] ;
		uu = (v - hvmin)/(hvmax - hvmin);
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		tb = sacdata[k].sachdr.rhdr[H_B];
		te = sacdata[k].sachdr.rhdr[H_E];
		to = sacdata[k].sachdr.rhdr[H_O];
		dv = sacdata[k].sachdr.rhdr[H_DELTA]/twin;
		if(kk == 0)
			ts = tb;
		if(rs_absolute == YES) {
			vv = (float)(sacdata[k].tzref + tb - gsac_control.begmin -tss)/twin;
			if(to != -12345.){
				vv = (float)(tb - to - tss)/twin ;
			}
			if(rs_doptaux){
				vv = vv - rs_dtdx*sacdata[k].sachdr.rhdr[H_DIST]/twin ;
			} else if(rs_doptaud){
				vv = vv - rs_dtdd*sacdata[k].sachdr.rhdr[H_GCARC]/twin ;
			}
		} else {
/* changed 18 MAR 2006 was tb - tss and did not do relative */
			vv = (ts-tss)/twin;
			if(rs_doptaux){
				vv = vv - rs_dtdx*sacdata[k].sachdr.rhdr[H_DIST]/twin ;
			} else if(rs_doptaud){
				vv = vv - rs_dtdd*sacdata[k].sachdr.rhdr[H_GCARC]/twin ;
			}
		}
		for(i = 0 , nii = 0; i < npts ; i++){
			if(uu >= 0.0 && uu <= 1.0 && vv >= 0.0 && vv <= 1.0){
				tv = sacdata[k].sac_data[i];
				if(prs_scale_relative == NO){
					if(prs_absolute_power == 0.0)
						tv *= 1.0;
					else if(prs_absolute_power == 1.0)
						tv *= v;
					else
					tv *= pow(v,prs_absolute_power);
				}
				if(tv > sacdata[k].winmax)
					sacdata[k].winmax = tv;
				if(tv < sacdata[k].winmin)
					sacdata[k].winmin = tv;
			}
			vv += dv;
		}
	
		if(sacdata[k].winmax > prs_abs_max_amp)
			prs_abs_max_amp = sacdata[k].winmax;
		if(sacdata[k].winmin < prs_abs_min_amp)
			prs_abs_min_amp = sacdata[k].winmax;
		}
}
