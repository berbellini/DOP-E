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
 * 22 JUN 2007 - corrected the y-offset in the REFR???.CTL to permit refrmod96 overlay
 * 11 JUL 2010 - added the ~ command which then creates a screen dump 
		with a file name DUMPxxx.PLT  The purpose is to get a hardcopy
		of the complete interactive menu for documentation
 * 14 MAR 2011 - clean up for actual use
 * 22 MAR 2011 - the filters are now 2-pole 1-pass, thus causal!!
 * */
#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include <unistd.h>
#include	"csstim.h"

/* external variables */
extern struct sacfile_ *sacdata;
extern int *sortptr;
extern float agc_window_pct  ;
extern int agc_do_pct ;
extern int sort_which;
extern int sort_type;
extern int sort_reverse;
extern int sort_do;
/* external routines invoked */
extern void gsac_exec_read(void);
extern void gsac_exec_dagc(void);
extern void gsac_exec_wh(void);
extern void gsac_exec_sort(void);
extern void dogrid(void);
extern void gsac_exec_hold(void);


#define REFR_LS	-1
#define	REFR_PO 	-2
#define	REFR_SE 	-3
#define	REFR_RE 	-4
#define	REFR_ABS	-5
#define	REFR_REL	-6
#define REFR_COL -17
#define REFR_ANN -18
#define LANDSCAPE 0
#define PORTRAIT 1
#define SEASCAPE 2
#define REVERSE 3
#define ABSOLUTE 4
#define RELATIVE 5

#define REFR_TYPE_E	-19
#define REFR_TYPE_R	-20
#define REFR_TYPE_T	-21
#define	EXPLORATION 1
#define REGIONAL 2
#define TELESEISM 3

/* default filter parameters */
#define BUTTERWORTH 0
#define HIGHPASS 2
#define LOWPASS 3
#define NUMPASS 1
#define NPOLE 2


static int refr_typeaxis = 0 ;	/* 0 LANDSCAPE, 1 PORTRAIT, 2 SEASCAPE */
static int refr_whichaxis = 50 ;	/* header variable.  */
static int refr_absolute = YES;
static char refr_str[1000];
static char refr_titstr[1000];
static char refr_annstr[100];
static char refr_strann[100];
static float refr_amp = 0.5 ;		/* amplitude in inches */
static int refr_exttitle = NO; /* use external title */
static int refr_ann = NO; /* annotate trace with STA etc */
static int refr_doptaux = NO;	/* apply p-tau for dt/dx */
static int refr_doptaud = NO;	/* apply p-tau for dt/dd */
static float refr_dtdx = 0.0;
static float refr_dtdd = 0.0;
static float refr_tmin = -1.0e+30;
static float refr_tmax =  1.0e+30;
static float refr_vmin = 0;
static float refr_vmax = 0;
static int   dotlim   = NO;
static int   dovlim   = NO;
static int refr_doshd = NO;
static int refr_doshdplmn = 0;
static int refr_shdcolor = 1;
static float *py = (float *)NULL;
static float *px = (float *)NULL;
static int refr_doclip = NO;
static float refr_cliplev = 2.0 ;
static int refr_dorefr = YES ;	/* in refraction processing mode
				   else in reflection processing mode */
static int refr_type = EXPLORATION ; /*used for filter settings */
static float refr_lpc;
static float refr_hpc;

static	double tss,tes;
static	float hvmin, hvmax;
static	float refrs_x0, refrs_y0, refrs_xlen, refrs_ylen, dy, dv;
static	char refr_pltname[12];

#define MXGATE 1000
int do_getoffset(int k,int ilw,int iup,int ngate,int npts,float *s,float *x);
static float s[MXGATE];
void do_outstr(float x0,float y0,float *xx,float *yy,char *str);



void  gsac_refrrs(float rs_dtdx,float rs_dtdd,int rs_absolute,int rs_doptaux,
	int rs_doptaud,float rs_tmin,float rs_tmax,float rs_vmin,float rs_vmax,
	float x0, float y0, float xlen, float ylen, float rs_amp, int rs_typeaxis,
	int rs_doshd,int rs_doshdplmn,int rs_ann,char *rs_annstr,int rs_exttitle,
	char *rs_titstr, char  *rs_tstr, int rs_tstrl, int rs_whichaxis, int rs_shdcolor,
	double *tss,double *tes,float *hvmin,float *hvmax,int dotlim,int dovlim,int labx,
	int rs_doclip,float rs_cliplev,int rs_shwpick);
int dorefrinteractive(float x0,float y0,float xlen,float ylen,
	int refr_typeaxis, 
	double tss,double tes,float hvmin,float hvmax);

/* title from header value title[bhdr - 40] */
static char *title[] = {
	"Receiver depth ",
	"",
	"",
	"",
	"",
	"Depth (km)",
	"",
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

struct arghdr refrarg[] = {
	{REFR_ANN   , "ANNOTATE"	, CHDR, 0, 1, NO, "ANNOTATE string ", 2},
	{REFR_COL   , "COLOR", IHDR, 0, 1, NO, "Color color ", 1},
	{REFR_TYPE_T  , "TEL", IHDR, 0, 0, NO, "Tel", 1},
	{REFR_TYPE_R  , "REG", IHDR, 0, 0, NO, "Reg", 1},
	{REFR_TYPE_E  , "EXP", IHDR, 0, 0, NO, "Exp", 1},
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
	{  0, ""	, IHDR,  0, 0, NO, "",-1}
};


/* these are temporary variables only used here */
float refr_real[10];
int   refr_int [10];
char  refr_tstr[100];
char  refr_vstr[100];
int   refr_tstrl;

extern void XviG_Flush();



void gsac_set_param_refr(int ncmd, char **cmdstr)
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
	if(testarg(ncmd, cmdstr, refrarg, NO, YES))
	       return	;
	for(i=0 ; refrarg[i].key[0] != '\0' ; i++){
		if(refrarg[i].used > 0){
			if(refrarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, refrarg[i].key, 
					refrarg[i].mfit,refrarg[i].narg, refr_real);
			} else if(refrarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, refrarg[i].key, 
					refrarg[i].mfit,refrarg[i].narg, refr_int );
			} else if(refrarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, refrarg[i].key, 
					refrarg[i].mfit,refrarg[i].narg, refr_str );
			}
			
			if(refrarg[i].id >=  0 )	{
				/* get what item to sort */
				if(refr_whichaxis != refrarg[i].id){
					refr_whichaxis = refrarg[i].id;
					refr_exttitle = NO;
				}
			} else {
				switch(refrarg[i].id){
					/* for this to work we must 
					 * */
					case REFR_ANN:
						gsac_strupr(refr_str);
						if(strcmp(refr_str,"OFF")==0){
							refr_ann = NO;
						} else {
							strcpy(refr_annstr,refr_str);
							refr_ann = YES;
						}
						break;
					case REFR_TYPE_E:
						refr_type = EXPLORATION ; 
						break;
					case REFR_TYPE_T:
						refr_type = TELESEISM ; 
						break;
					case REFR_TYPE_R:
						refr_type = REGIONAL ; 
						break;
					case REFR_COL:
						refr_shdcolor = refr_int[0] ;
						break;

				}
			}
		}
	}
}

void gsac_exec_refr(void)
{
	int ntrc;
	float tv;	/* normalized trace value [-1,1] */
	int iret;

	gsac_control.prs = 1;
       /* this routine only sets processing states
	*         printf("gsac_control.prs %d\n",gsac_control.prs);
	*                 */
	if(gsac_control.refrpick == NULL)
		gsac_control.refrpick = fopen("refrpick.tmp","w+");


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
	} else {
		printf("Must be in interactive mode: bg x!\n");
		return;
	}

	refrs_xlen = gsac_control.xlen ;
	refrs_ylen = gsac_control.ylen ;
	refrs_x0   = gsac_control.x0 ;
	refrs_y0   = gsac_control.y0 + 0.8 ;
	refr_hpc = -1. ;
	refr_lpc = -1. ;
	agc_window_pct = -1. ;
	

	if(gsac_control.hold == NO)gframe(2);

	/* force s sort by distance */
	sort_do = YES;
	sort_reverse = NO;
	sort_type = RHDR;
	sort_which = H_DIST;
	gsac_exec_sort();

	/* put up initial display */
	iret = YES;
	while(iret == YES){
		gsac_refrrs(refr_dtdx,refr_dtdd,refr_absolute,refr_doptaux,
			refr_doptaud,refr_tmin,refr_tmax,refr_vmin,
			refr_vmax, refrs_x0, refrs_y0, refrs_xlen, refrs_ylen, refr_amp,
			refr_typeaxis,refr_doshd,refr_doshdplmn,refr_ann,
			refr_annstr, YES, " ", refr_tstr, refr_tstrl, 
			refr_whichaxis, refr_shdcolor, &tss, &tes, 
			&hvmin,&hvmax,dotlim,dovlim,NO,refr_doclip,
			refr_cliplev, YES);
		iret = dorefrinteractive(refrs_x0,refrs_y0,refrs_xlen,refrs_ylen,refr_typeaxis,
			tss,tes,hvmin,hvmax);
	}

}

/* menu routines */
#include "nmenu.h"

void clearregion(float xl, float yl, float xh, float yh);
void show_menu (float x0, float y0, struct menu *m, int size, int *nm);
int inside(float xv, float yv, float xlb, 
	float ylb, float xhb, float yhb);
int proc_menu(float cx, float cy, int nm, struct menu *p);
int menu_sel(float x0, float y0, struct menu *men, int sizeofmen, char *mesg);
int dorefrpick(char *phase,float x0,float y0,float xlen,float ylen,
	int refr_typeaxis, float refr_dtdx,
	double tss,double tes,float hvmin,float hvmax,float *t0,float* p,
	float *x1, float *t1, float *x2, float *t2);
int doreflpick(char *phase,float x0,float y0,float xlen,float ylen,
	int refr_typeaxis, float refr_dtdx,
	double tss,double tes,float hvmin,float hvmax,float *t0,float* p,
	float *x1, float *t1, float *x2, float *t2);
void refrshwpick(float x0,float y0,float xlen,float ylen,
	double tss,double tes,float hvmin,float hvmax,float refr_dtdx);
int do_zoomm(float x0,float y0,float xlen,float ylen,double tss,double tes,float hvmin,float hvmax);
float pinterp(float x,float x1,float x2);
void do_prefine(float p,float t0,float xx1,float tt1,float xx2,float tt2,float gate,float ttlw,float tthg,float tvlw,float tvhg,int isp);
int dops_refr(float x0,float y0,float xlen,float ylen,double tss,double tes,float hvmin,float hvmax,char *msg,char *phase);
int dops_refl(float x0,float y0,float xlen,float ylen,double tss,double tes,float hvmin,float hvmax,char *msg,char *phase);
void do_setat0(float p,float t0,float xx1,float xx2,char *phase);
void do_hardcopy_refr(void);



#define MENU_REFR_Y  YES
#define MENU_REFR_N  NO
static struct menu refryn[] = {
		{  -1.0, -1.0, -1.0, -1.0, "Yes\0" , MENU_REFR_Y , -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "No \0" , MENU_REFR_N , -1, 1, 1}
};


#define MENU_REFR_PSY  YES
#define MENU_REFR_PSN  NO
#define MENU_REFR_PSREFINE 2
static struct menu refrpsyn[] = {
	{  -1.0, -1.0, -1.0, -1.0, "Yes\0" , MENU_REFR_PSY , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "No \0" , MENU_REFR_PSN , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Refine\0" , MENU_REFR_PSREFINE , -1, 1, 1}
};

#define MENU_REFR_REF1  1
#define MENU_REFR_REF2  2
#define MENU_REFR_REF3  3
#define MENU_REFR_REF4  4
#define MENU_REFR_REF5  5
#define MENU_REFR_REF6  6
#define MENU_REFR_REF7  7
#define MENU_REFR_REF8  8
static struct menu prref[] = {
	{  -1.0, -1.0, -1.0, -1.0, "1\0", MENU_REFR_REF1 , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "2\0", MENU_REFR_REF2 , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "3\0", MENU_REFR_REF3 , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "4\0", MENU_REFR_REF4 , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "5\0", MENU_REFR_REF5 , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "6\0", MENU_REFR_REF6 , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "7\0", MENU_REFR_REF7 , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "8\0", MENU_REFR_REF8 , -1, 1, 1}
};

struct refr_tim {
	int set;
	float t0;
	float p;
	float x1;
	float t1;
	float x2;
	float t2;
};

#define NUM_PICK 8

static struct refr_tim gsac_prefr[8] = {
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 }
};

static struct refr_tim gsac_srefr[8] = {
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
	{ -1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 }
};



#define MENU_REFR_GOREFL	1
#define MENU_REFR_DOREFR	2
#define MENU_REFR_EXIT		3
#define MENU_REFR_RESET		4
#define MENU_REFR_ZOOM		5
#define MENU_REFR_UNZOOM	6
#define MENU_REFR_AMPDEC	7
#define MENU_REFR_AMPINC	8
#define MENU_REFR_SHADE		9
#define MENU_REFR_VEL		10
#define MENU_REFR_CLIP		11
#define MENU_REFR_READ		12
#define MENU_REFR_HP		13
#define MENU_REFR_LP		14
#define MENU_REFR_AGC		15
#define MENU_REFR_HARDCOPY	16
#define MENU_REFR_GOREFR	17
#define MENU_REFR_DOREFL	18
static struct menu refrph[] = {
	{  -1.0, -1.0, -1.0, -1.0, "DoRefr\0"    , MENU_REFR_DOREFR    , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "GoRefl\0"    , MENU_REFR_GOREFL    , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Zoom\0" , MENU_REFR_ZOOM , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Unzoom\0" , MENU_REFR_UNZOOM , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "-Amp\0" , MENU_REFR_AMPDEC , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "+Amp\0" , MENU_REFR_AMPINC , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Shd\0" , MENU_REFR_SHADE , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Vel\0" , MENU_REFR_VEL , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Clip\0" , MENU_REFR_CLIP , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "HP\0" , MENU_REFR_HP , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "LP\0" , MENU_REFR_LP , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "AGC\0" , MENU_REFR_AGC , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "ReRead\0" , MENU_REFR_READ , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Reset\0" , MENU_REFR_RESET , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Hardcopy\0" , MENU_REFR_HARDCOPY , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Exit\0" , MENU_REFR_EXIT , -1, 1, 1}
};

static struct menu reflph[] = {
	{  -1.0, -1.0, -1.0, -1.0, "DoRefl\0"    , MENU_REFR_DOREFL    , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "GoRefr\0"    , MENU_REFR_GOREFR    , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Zoom\0" , MENU_REFR_ZOOM , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Unzoom\0" , MENU_REFR_UNZOOM , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "-Amp\0" , MENU_REFR_AMPDEC , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "+Amp\0" , MENU_REFR_AMPINC , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Shd\0" , MENU_REFR_SHADE , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Clip\0" , MENU_REFR_CLIP , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "HP\0" , MENU_REFR_HP , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "LP\0" , MENU_REFR_LP , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "AGC\0" , MENU_REFR_AGC , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "ReRead\0" , MENU_REFR_READ , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Reset\0" , MENU_REFR_RESET , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Hardcopy\0" , MENU_REFR_HARDCOPY , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Exit\0" , MENU_REFR_EXIT , -1, 1, 1}
};


#define MENU_SHD_POS  1
#define MENU_SHD_NEG  2
#define MENU_SHD_ALL  3
#define MENU_SHD_OFF  4
static struct menu refrshd[] = {
	{  -1.0, -1.0, -1.0, -1.0, "Pos\0"    , MENU_SHD_POS    , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Neg\0"    , MENU_SHD_NEG    , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Both\0"   , MENU_SHD_ALL    , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "Off\0"    , MENU_SHD_OFF    , -1, 1, 1}
};

#define MENU_PS_P 1
#define MENU_PS_S 2

static struct menu refrps[] = {
	{  -1.0, -1.0, -1.0, -1.0, "P\0"    , MENU_PS_P    , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "S\0"    , MENU_PS_S    , -1, 1, 1},
};

#define MENU_VEL_OFF  0
#define MENU_VEL_01  1
#define MENU_VEL_02  2
#define MENU_VEL_03  3
#define MENU_VEL_04  4
#define MENU_VEL_05  5
#define MENU_VEL_06  6
#define MENU_VEL_07  7
#define MENU_VEL_08  8
#define MENU_VEL_09  9
#define MENU_VEL_10  10
#define MENU_VEL_15  15
#define MENU_VEL_20  20
#define MENU_VEL_25  25
#define MENU_VEL_30  30
#define MENU_VEL_35  35
#define MENU_VEL_40  40
#define MENU_VEL_50  50
#define MENU_VEL_60  60
#define MENU_VEL_70  70
#define MENU_VEL_80  80
#define MENU_VEL_90  90
static struct menu refrvel[] = {
	{  -1.0, -1.0, -1.0, -1.0, "OFF\0"    , MENU_VEL_OFF  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.1\0"    , MENU_VEL_01  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.2\0"    , MENU_VEL_02  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.3\0"    , MENU_VEL_03  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.4\0"    , MENU_VEL_04  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.5\0"    , MENU_VEL_05  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.6\0"    , MENU_VEL_06  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.7\0"    , MENU_VEL_07  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.8\0"    , MENU_VEL_08  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.9\0"    , MENU_VEL_09  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "1.0\0"    , MENU_VEL_10  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "1.5\0"    , MENU_VEL_15  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "2.0\0"    , MENU_VEL_20  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "2.5\0"    , MENU_VEL_25  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "3.0\0"    , MENU_VEL_30  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "3.5\0"    , MENU_VEL_35  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "4.0\0"    , MENU_VEL_40  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "5.0\0"    , MENU_VEL_50  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "6.0\0"    , MENU_VEL_60  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "7.0\0"    , MENU_VEL_70  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "8.0\0"    , MENU_VEL_80  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "9.0\0"    , MENU_VEL_90  , -1, 1, 1}
};

#define MENU_CLIP_OFF  0
#define MENU_CLIP_01   5
#define MENU_CLIP_02  10
#define MENU_CLIP_03  20
#define MENU_CLIP_04  30
#define MENU_CLIP_05  40
#define MENU_CLIP_06  50
#define MENU_CLIP_07  100
static struct menu refrclip[] = {
	{  -1.0, -1.0, -1.0, -1.0, "OFF\0"    , MENU_CLIP_OFF  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.5\0"    , MENU_CLIP_01  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 1\0"    , MENU_CLIP_02  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 2\0"    , MENU_CLIP_03  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 3\0"    , MENU_CLIP_04  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 4\0"    , MENU_CLIP_05  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 5\0"    , MENU_CLIP_06  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "10\0"    , MENU_CLIP_07  , -1, 1, 1}
};

#define MENU_HP_OFF  0
#define MENU_HP_00   5
#define MENU_HP_01  10
#define MENU_HP_02  20
#define MENU_HP_03  30
#define MENU_HP_04  40
#define MENU_HP_05  50
#define MENU_HP_06  60
#define MENU_HP_07  70
#define MENU_HP_08  80
#define MENU_HP_09  90
#define MENU_HP_10  100
#define MENU_HP_11  150
#define MENU_HP_12  200
#define MENU_HP_13  250
#define MENU_HP_14  500
static struct menu refrhpe[] = {
	{  -1.0, -1.0, -1.0, -1.0, "OFF\0"    , MENU_HP_OFF  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "  5\0"    , MENU_HP_00  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 10\0"    , MENU_HP_01  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 20\0"    , MENU_HP_02  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 30\0"    , MENU_HP_03  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 40\0"    , MENU_HP_04  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 50\0"    , MENU_HP_05  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 60\0"    , MENU_HP_06  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 70\0"    , MENU_HP_07  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 80\0"    , MENU_HP_08  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 90\0"    , MENU_HP_09  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "100\0"    , MENU_HP_10  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "150\0"    , MENU_HP_11  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "200\0"    , MENU_HP_12  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "250\0"    , MENU_HP_13  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "500\0"    , MENU_HP_14  , -1, 1, 1}
};
static struct menu refrhpr[] = {
	{  -1.0, -1.0, -1.0, -1.0, "OFF\0"   , MENU_HP_OFF  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".5\0"    , MENU_HP_00  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 1\0"    , MENU_HP_01  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 2\0"    , MENU_HP_02  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 3\0"    , MENU_HP_03  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 4\0"    , MENU_HP_04  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 5\0"    , MENU_HP_05  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 6\0"    , MENU_HP_06  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 7\0"    , MENU_HP_07  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 8\0"    , MENU_HP_08  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 9\0"    , MENU_HP_09  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "10\0"    , MENU_HP_10  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "15\0"    , MENU_HP_11  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "20\0"    , MENU_HP_12  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "25\0"    , MENU_HP_13  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "50\0"    , MENU_HP_14  , -1, 1, 1}
};
static struct menu refrhpt[] = {
	{  -1.0, -1.0, -1.0, -1.0, "OFF\0"   , MENU_HP_OFF  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".01\0"    , MENU_HP_01  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".02\0"    , MENU_HP_02  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".03\0"    , MENU_HP_03  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".04\0"    , MENU_HP_04  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".05\0"    , MENU_HP_05  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".06\0"    , MENU_HP_06  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".07\0"    , MENU_HP_07  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".08\0"    , MENU_HP_08  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".09\0"    , MENU_HP_09  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".10\0"    , MENU_HP_10  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".15\0"    , MENU_HP_11  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".20\0"    , MENU_HP_12  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".25\0"    , MENU_HP_13  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".50\0"    , MENU_HP_14  , -1, 1, 1}
};


#define MENU_LP_OFF  0
#define MENU_LP_00   5
#define MENU_LP_01  10
#define MENU_LP_02  20
#define MENU_LP_03  30
#define MENU_LP_04  40
#define MENU_LP_05  50
#define MENU_LP_06  60
#define MENU_LP_07  70
#define MENU_LP_08  80
#define MENU_LP_09  90
#define MENU_LP_10  100
#define MENU_LP_11  150
#define MENU_LP_12  200
#define MENU_LP_13  250
#define MENU_LP_14  500
static struct menu refrlpe[] = {
	{  -1.0, -1.0, -1.0, -1.0, "OFF\0"    , MENU_LP_OFF  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "  5\0"    , MENU_LP_00  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 10\0"    , MENU_LP_01  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 20\0"    , MENU_LP_02  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 30\0"    , MENU_LP_03  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 40\0"    , MENU_LP_04  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 50\0"    , MENU_LP_05  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 60\0"    , MENU_LP_06  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 70\0"    , MENU_LP_07  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 80\0"    , MENU_LP_08  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 90\0"    , MENU_LP_09  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "100\0"    , MENU_LP_10  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "150\0"    , MENU_LP_11  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "200\0"    , MENU_LP_12  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "250\0"    , MENU_LP_13  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "500\0"    , MENU_LP_14  , -1, 1, 1}
};
static struct menu refrlpr[] = {
	{  -1.0, -1.0, -1.0, -1.0, "OFF\0"   , MENU_LP_OFF  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".5\0"    , MENU_LP_00  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 1\0"    , MENU_LP_01  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 2\0"    , MENU_LP_02  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 3\0"    , MENU_LP_03  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 4\0"    , MENU_LP_04  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 5\0"    , MENU_LP_05  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 6\0"    , MENU_LP_06  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 7\0"    , MENU_LP_07  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 8\0"    , MENU_LP_08  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, " 9\0"    , MENU_LP_09  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "10\0"    , MENU_LP_10  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "15\0"    , MENU_LP_11  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "20\0"    , MENU_LP_12  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "25\0"    , MENU_LP_13  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "50\0"    , MENU_LP_14  , -1, 1, 1}
};
static struct menu refrlpt[] = {
	{  -1.0, -1.0, -1.0, -1.0, "OFF\0"    , MENU_LP_OFF  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".01\0"    , MENU_LP_01  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".02\0"    , MENU_LP_02  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".03\0"    , MENU_LP_03  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".04\0"    , MENU_LP_04  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".05\0"    , MENU_LP_05  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".06\0"    , MENU_LP_06  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".07\0"    , MENU_LP_07  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".08\0"    , MENU_LP_08  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".09\0"    , MENU_LP_09  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".10\0"    , MENU_LP_10  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".15\0"    , MENU_LP_11  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".20\0"    , MENU_LP_12  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".25\0"    , MENU_LP_13  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, ".50\0"    , MENU_LP_14  , -1, 1, 1}
};

#define MENU_AGC_OFF  0
#define MENU_AGC_01  5
#define MENU_AGC_02  10
#define MENU_AGC_03  15
#define MENU_AGC_04  20
#define MENU_AGC_05  25
#define MENU_AGC_06  30
#define MENU_AGC_07  35
#define MENU_AGC_08  40
#define MENU_AGC_09  45
#define MENU_AGC_10  50
static struct menu refragc[] = {
	{  -1.0, -1.0, -1.0, -1.0, "OFF\0"    , MENU_AGC_OFF  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.05\0"    , MENU_AGC_01  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.10\0"    , MENU_AGC_02  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.15\0"    , MENU_AGC_03  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.20\0"    , MENU_AGC_04  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.25\0"    , MENU_AGC_05  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.30\0"    , MENU_AGC_06  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.35\0"    , MENU_AGC_07  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.40\0"    , MENU_AGC_08  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.45\0"    , MENU_AGC_09  , -1, 1, 1},
	{  -1.0, -1.0, -1.0, -1.0, "0.50\0"    , MENU_AGC_10  , -1, 1, 1},
};

int dorefrinteractive(float x0,float y0,float xlen,float ylen,
	int refr_typeaxis, 
	double tss,double tes,float hvmin,float hvmax)
{
int i, cmd;
int j;
float lx, ly, hx, hy;
int iret, rret, psret;
int nmph;	/* number of entires in refrph menu */
float t0, p;
float gate;
float xx1, tt1, xx2, tt2 ; /* limits used to define the line */
float tvmin, tvmax, ttmin, ttmax;
float xx, yy;
float vel;

	if(refr_dorefr == YES){
		show_menu(0.25, 1.2, refrph, sizeof(refrph),&nmph);
	} else {
		show_menu(0.25, 1.2, reflph, sizeof(reflph),&nmph);
	}
	if(gsac_control.plotdevice == WIN) XviG_Flush();
	/* now go interactive to get menu choice */
	/* set clip region to the actual plot region not the entire box */
	/* cross hair cursors within the trace plot region , default
	 * 	arrow outside this region 
	 * if acceptable then display this line - also if 
	 * interactive and phases exist display
	 * */
	/* get dtdx from global */
	/* infinite loop */
	for(;;){
		/* display pick lines here */
		/*
		refrshwpick(x0,y0,xlen,ylen,tss,tes,hvmin,hvmax,refr_dtdx);
		*/
		sleep(1);
		if(gsac_control.plotdevice == WIN) XviG_Flush();
		if(refr_dorefr == YES){
			iret = menu_sel(0.5,0.5,refrph, sizeof(refrph),"Select ");
		} else {
			iret = menu_sel(0.5,0.5,reflph, sizeof(reflph),"Select ");
		}
		if(iret == MENU_REFR_EXIT ){
			rret = menu_sel(0.5,0.5,refryn,sizeof(refryn),
				"Write SAC Headers? ");
			if(rret == YES){
				printf("Arrival time picks written to SAC header\n");
				gsac_exec_wh();
			} else {
				printf("Arrival time picks not written to SAC header\n");
			}
			return NO;
		} else if(iret == MENU_REFR_GOREFR ){
			refr_dorefr = YES;
			return YES;
		} else if(iret == MENU_REFR_GOREFL ){
			refr_dorefr = NO;
			refr_doptaux = NO;
			refr_doptaud = NO ;
			refr_dtdx = 0.0 ;
			strcpy(refr_vstr, "");
			return YES;
		} else if(iret == MENU_REFR_DOREFR ){
			/* Ask P or S , then get Refraction */
			psret = menu_sel(0.5,0.8,refrps,sizeof(refrps),
				"Select P or S");
			switch(psret){
				case MENU_PS_P:
				return dops_refr(x0,y0,xlen,ylen,tss,tes,hvmin,hvmax,"Select P Refraction","P");
				break;
				case MENU_PS_S:
				return dops_refr(x0,y0,xlen,ylen,tss,tes,hvmin,hvmax,"Select S Refraction","S");
				break;
			}
		} else if(iret == MENU_REFR_DOREFL ){
			/* Ask P or S, then get Refraction, then nultiple */
			psret = menu_sel(0.5,0.8,refrps,sizeof(refrps),
				"Select P or S");
			switch(psret){
				case MENU_PS_P:
				return dops_refl(x0,y0,xlen,ylen,tss,tes,hvmin,hvmax,"Select P Reflection","P");
				break;
				case MENU_PS_S:
				return dops_refl(x0,y0,xlen,ylen,tss,tes,hvmin,hvmax,"Select S Reflection","S");
				break;
			}
		} else if(iret == MENU_REFR_AMPDEC ){
					refr_amp /= 2.0 ;
			return YES;
		} else if(iret == MENU_REFR_AMPINC ){
					refr_amp *= 2.0 ;
			return YES;
		} else if(iret == MENU_REFR_VEL ){
			rret = menu_sel(0.5,0.8,refrvel,sizeof(refrvel),
				"Select reduction velocity");
			switch(rret){
				case MENU_VEL_OFF:
					refr_doptaux = NO;
					refr_doptaud = NO ;
					refr_dtdx = 0.0 ;
					strcpy(refr_vstr, "");
					break;
				default:
					refr_doptaux = YES;
					refr_doptaud = NO ;
					/* to get the correct velocity
					 * I will assume nothing about the
					 * naming convention and instead
					 * use the contents of the
					 * string in the data menu structure */
					for(i = 0 ; 
					i < sizeof(refrvel)/sizeof(struct menu); i++){
					if(refrvel[i].action == rret){
						strcpy(refr_vstr, refrvel[i].str);
						vel=strtod(refr_vstr,(char **)NULL);
					refr_dtdx = 1.0/vel;
					}
}
					break;
			}
			return YES;
		} else if(iret == MENU_REFR_CLIP ){
			rret = menu_sel(0.5,0.5,refrclip,sizeof(refrclip),
				"Select clip level");
			switch(rret){
				case MENU_CLIP_OFF:
					refr_doclip = NO;
					break;
				default:
					refr_doclip = YES;
					refr_cliplev = 0.1*(float)rret;
					break;
			}
			return YES;
		} else if(iret == MENU_REFR_SHADE ){
			rret = menu_sel(0.5,0.5,refrshd,sizeof(refrshd),
				"Select Trace Shading");
			switch(rret){
				case MENU_SHD_POS:
					refr_doshdplmn = 0;
					refr_doshd = YES;
					break;
				case MENU_SHD_NEG:
					refr_doshdplmn = 1;
					refr_doshd = YES;
					break;
				case MENU_SHD_ALL:
					refr_doshdplmn = 2;
					refr_doshd = YES;
					break;
				case MENU_SHD_OFF:
					refr_doshd = NO;
					break;
			}
			return YES;
		} else if(iret == MENU_REFR_READ ){
			gsac_exec_read();
			return YES;
		} else if(iret == MENU_REFR_HP ){
			if(refr_type == EXPLORATION){
				rret = menu_sel(0.5,0.5,refrhpe,sizeof(refrhpe),
					"Define High Pass Corner");
			} else if(refr_type == REGIONAL){
				rret = menu_sel(0.5,0.5,refrhpr,sizeof(refrhpr),
					"Define High Pass Corner");
			} else if(refr_type == TELESEISM){
				rret = menu_sel(0.5,0.5,refrhpt,sizeof(refrhpt),
					"Define High Pass Corner");
			}
			switch(rret){
				case MENU_HP_OFF:
					gsac_exec_read();
					refr_hpc = -1. ;
					break;
				default:
					if(refr_type == EXPLORATION){
						refr_hpc = (float)rret;
					} else if(refr_type == REGIONAL){
						refr_hpc = 0.1*(float)rret;
					} else if(refr_type == TELESEISM){
						refr_hpc = 0.001*(float)rret;
					}
					gsac_filt(refr_hpc, 1.0e+10, NPOLE, NUMPASS, HIGHPASS, BUTTERWORTH, 1.0);
					break;
			}
			return YES;
		} else if(iret == MENU_REFR_LP ){
			if(refr_type == EXPLORATION){
				rret = menu_sel(0.5,0.5,refrlpe,sizeof(refrlpe),
					"Define Low Pass Corner");
			} else if(refr_type == REGIONAL){
				rret = menu_sel(0.5,0.5,refrlpr,sizeof(refrlpr),
					"Define Low Pass Corner");
			} else if(refr_type == TELESEISM){
				rret = menu_sel(0.5,0.5,refrlpt,sizeof(refrlpt),
					"Define Low Pass Corner");
			}
			switch(rret){
				case MENU_LP_OFF:
					gsac_exec_read();
					refr_lpc = -1. ;
					break;
				default:
					if(refr_type == EXPLORATION){
						refr_lpc = (float)rret;
					} else if(refr_type == REGIONAL){
						refr_lpc = 0.1*(float)rret;
					} else if(refr_type == TELESEISM){
						refr_lpc = 0.001*(float)rret;
					}
					gsac_filt(0.0, refr_lpc,  NPOLE, NUMPASS, LOWPASS, BUTTERWORTH, 1.0);
					break;
			}
			return YES;
		} else if(iret == MENU_REFR_AGC ){
			rret = menu_sel(0.5,0.5,refragc,sizeof(refragc),
				"Define AGC Window Percent ");
			switch(rret){
				case MENU_AGC_OFF:
					gsac_exec_read();
					break;
				default:
					agc_window_pct = (float)rret;
					agc_do_pct = YES;
					gsac_exec_dagc();
					break;
			}
			return YES;
		} else if(iret == MENU_REFR_RESET ){
			refr_amp = 0.5 ;
			refr_exttitle = NO;
			refr_ann = NO;
			refr_doptaux = NO;
			refr_doptaud = NO;
			refr_dtdx = 0.0;
			refr_dtdd = 0.0;
			refr_tmin = -1.0e+30;
			refr_tmax =  1.0e+30;
			refr_vmin = 0;
			refr_vmax = 0;
			dotlim   = NO;
			dovlim   = NO;
			refr_doshd = NO;
			refr_doshdplmn = 0;
			refr_shdcolor = 1;
			refr_absolute = YES ;
			refr_doclip = NO;
			refr_cliplev = 2.0;
			return YES;
		} else if(iret == MENU_REFR_ZOOM ){
			/* this is from do_mft */
			gcursor("Arrow");
			rret = do_zoomm(x0,y0,xlen,ylen,tss,tes,hvmin,hvmax);
			gcursor("Arrow");
			if(rret == YES){
				dotlim = YES;
				dovlim = YES;
				return YES;
			}
		} else if(iret == MENU_REFR_UNZOOM ){
			/* this is from do_mft */
			refr_tmin = -1.0e+30;
			refr_tmax =  1.0e+30;
			refr_vmin = 0;
			refr_vmax = 0;
			dotlim = NO;
			dovlim = NO;
			return YES;
		} else if(iret == MENU_REFR_HARDCOPY ){
			do_hardcopy_refr();
		}
	}
	return NO;
}

int do_zoomm(float x0,float y0,float xlen,float ylen,double tss,double tes,float hvmin,float hvmax)
{
	/* the return indicates whether to replot everything */
	char c[2];
	float xv1, xv2, yv1, yv2;
	float x1, y1, x2, y2;
	float xfrac, yfrac;
	float px1,px2,py1,py2;
	float tx1,tx2,ty1,ty2;
	/* turn off XOR */
	newpen(3000);
	curaxy(&xv1, &yv1, c);
	if( !inside(xv1,yv1, x0,y0,x0+xlen,y0+ylen))
		return NO;
	gclip("on",x0,y0,x0+xlen,y0+ylen);
	gcursor("Box");
	curaxy(&xv2, &yv2, c);
	gclip("off",x0,y0,x0+xlen,y0+ylen);
	gcursor("Arrow");
	if(xv1==xv2 && yv1 == yv2) return NO;  /* need two points for a line */
	/* 
	We work with three coordinate systems, each consisting of two
	opposite corners of a view rectangle

	(pxl,pyl) -> (pxh,pyl)   is display screen which does not change
	(opxl,opyl) -> (opxh,opyh) is the mapped portion of the original plot
	(apxl,apyl) -> (apxh,apyl) are user coordinates of the original space

	To do the mapping,use a simple interpolation

	V(p) = (1-p)V1 + pV2   where  0 <= p <= 1
	*/
	if( inside(xv1,yv1, x0,y0,x0+xlen,y0+ylen)
	    && inside(xv2,yv2, x0,y0,x0+xlen,y0+ylen)) {
		/* these are screen coordinates */
		x1 = MIN(xv1,xv2) ;
		x2 = MAX(xv1,xv2) ;
		y1 = MIN(yv1,yv2) ;
		y2 = MAX(yv1,yv2) ;
		px1 = pinterp(x1,x0,x0+xlen);
		px2 = pinterp(x2,x0,x0+xlen);
		py1 = pinterp(y1,y0,y0+ylen);
		py2 = pinterp(y2,y0,y0+ylen);
		/* now estimate the coordinates in the
					original plot space */
		tx1 = (1.0 - px1) * hvmin + px1 * hvmax;
		tx2 = (1.0 - px2) * hvmin + px2 * hvmax;
		ty1 = (1.0 - py1) * tss + py1 * tes;
		ty2 = (1.0 - py2) * tss + py2 * tes;
		/* update our perception of the view of original space */

		refr_vmin = tx1;
		refr_vmax = tx2;
		refr_tmin = ty1;
		refr_tmax = ty2;
		return YES;
	} else {
		return NO;
	}
}

int menu_sel(float x0, float y0, struct menu *men, int sizeofmen, char *mesg)
{
float cx, cy;
char c[2];
int iret;
int i;
float xl, yl, xh, yh;
int nmen;
int switchchar;
	show_menu(x0, y0, men, sizeofmen,&nmen);
	xl = men[0].xl;
	yl = men[0].yl;
	xh = men[0].xh;
	yh = men[0].yh;
	for(i = 0 ; i < nmen; i++){
		if(men[i].xl < xl)xl = men[i].xl;
		if(men[i].yl < yl)yl = men[i].yl;
		if(men[i].xh > xh)xh = men[i].xh;
		if(men[i].yh > yh)yh = men[i].yh;
	}
	if(gsac_control.plotdevice == WIN) XviG_Flush();
	gmesg(mesg);
	gcursor("Arrow");
	for( ; ; ){
		curaxy(&cx, &cy, c);
		switchchar = c[0];
		switch(switchchar){
			case '~':
				sprintf(refr_pltname,"DUMP%3.3d.PLT",
					gsac_control.plotcount_dump);
				gsac_control.plotcount_dump++;
				ginitf(refr_pltname,"GSAC");
				/* start plot mode */
        			gsac_control.hold = YES;
        			gsac_control.plotdevice = PLT;
        			gsac_control.everinteractive = NO ;
				gsac_refrrs(refr_dtdx,refr_dtdd,refr_absolute,refr_doptaux,
					refr_doptaud,refr_tmin,refr_tmax,refr_vmin,
					refr_vmax, refrs_x0, refrs_y0, refrs_xlen, refrs_ylen, refr_amp,
					refr_typeaxis,refr_doshd,refr_doshdplmn,refr_ann,
					refr_annstr, YES, " ", refr_tstr, refr_tstrl, 
					refr_whichaxis, refr_shdcolor, &tss, &tes, 
					&hvmin,&hvmax,dotlim,dovlim,NO,refr_doclip,
					refr_cliplev, YES);
				show_menu(x0, y0, men, sizeofmen,&nmen);


				gsac_control.plotdevice = WIN;
				gsac_control.hold = NO;
				gsac_control.everinteractive = YES ;
				ginitf("INTEM","GSAC");

				break;
			default:
				iret=proc_menu(cx,cy,nmen,men);
				if(iret >= 0){
					clearregion(xl,yl,xh,yh);
					if(gsac_control.plotdevice == WIN) XviG_Flush();
				return iret;
			}
		}
	}
}

int proc_menu(float cx, float cy, int nm, struct menu *p)
{
	/* return the index into the penu structure action field */
int i;
int cmd;
float lx, ly, hx, hy;
	for(i=0 ; i < nm ; i++){
		lx = p[i].xl;
		ly = p[i].yl;
		hx = p[i].xh;
		hy = p[i].yh;
		cmd = p[i].action;
		if(inside(cx,cy,lx,ly,hx,hy))
			return cmd;
	}
	return -1;
}

void refrshwpick(float x0,float y0,float xlen,float ylen,
	double tss,double tes,float hvmin,float hvmax,float refr_dtdx)
{
	/* plot the picks - use a solid line for the actual range - dashed for
	 * extrapolation
	 * */
int kolor;
float x1, y1, x2, y2;
int i,j;
float xx1,tt1,xx2,tt2;
int set;
	kolor = 2;
	gclip("on" ,x0,y0,x0+xlen,y0+ylen);
	for(j=0;j < 2 ; j ++){
		for(i=0;i < NUM_PICK ; i++){
			if(j == 0){
				set = gsac_prefr[i].set ;
				xx1 = gsac_prefr[i].x1;
				tt1 = gsac_prefr[i].t1;
				xx2 = gsac_prefr[i].x2;
				tt2 = gsac_prefr[i].t2;
			} else if (j==1){
				set = gsac_srefr[i].set ;
				xx1 = gsac_srefr[i].x1;
				tt1 = gsac_srefr[i].t1;
				xx2 = gsac_srefr[i].x2;
				tt2 = gsac_srefr[i].t2;
			}
			if(set > 0) {
				tt1 = tt1 - refr_dtdx * ABS(xx1);
				tt2 = tt2 - refr_dtdx * ABS(xx2);
				x1 = x0 + xlen*(xx1 - hvmin)/(hvmax - hvmin);
				x2 = x0 + xlen*(xx2 - hvmin)/(hvmax - hvmin);
				y1 = y0 + ylen*(tt1 - tss)/(tes - tss);
				y2 = y0 + ylen*(tt2 - tss)/(tes - tss);
				newpen(2);
				plot(x1,y1,3); plot(x2,y2,2);plot(x2,y2,3);
				newpen(1);
			}
		}
	}
	gclip("off",x0,y0,x0+xlen,y0+ylen);
}

float pinterp(float x,float x1,float x2)
{
	/* linear interpolation parameter */
	/* get the value of p in   x = (1-p)*x1 + p x2 */
	if(x1 == x2){
		return (0.0);
	} else {
		return (  (x - x1) / ( x2 - x1) );
	}
}


void do_prefine(float p,float t0,float xx1,float tt1,float xx2,float tt2,float gate,float ttlw,float tthg,float tvlw,float tvhg,int isp)
{
	/* p	ray parameter
	 * t0	intercept time
	 * (xx1,tt1) - (xx2,tt2) limits of line for fit
	 * (tvlw,ttlw) - (tvhg,tthg) - window to define stacking function
	 * isp	P or S
	 * */
	int i,j,k,l, ilw, iup, ngate, ntrc, nstk;
	int npts;
	int offset;
	float delta, dist;
	double tpred;
	float peak, tmp;
	double tss,tes;
	/* use the linear estimate to start a correlation procedure to
	 * determine arrival times
	 * This is a two step procedure. First stack along the line, and
	 * then use this stacked waveform to get bettern arrival times through
	 * a cross-correlation. A pick will only be made if the peak
	 * cross-correlation is within given factor of the auto-correlation
	 * */

	ntrc = gsac_control.number_itraces;
	/* get the delta and the number of points from the
	 * first trace - this should all be
	 * the same */
	if(ntrc < 1 )
		return;

	/* initialize the stack */
	nstk = 0;
	for(i=0 ; i < MXGATE; i++)
		s[i] = 0.0;

	npts = sacdata[0].sachdr.ihdr[H_NPTS];
	delta = sacdata[0].sachdr.rhdr[H_DELTA];
	for ( k=0 ; k < ntrc ; k++){
		npts  = sacdata[k].sachdr.ihdr[H_NPTS];
		delta = sacdata[k].sachdr.rhdr[H_DELTA];
		dist  = sacdata[k].sachdr.rhdr[H_DIST ];
		tss = sacdata[k].sachdr.rhdr[H_B] -
			sacdata[k].sachdr.rhdr[H_O];
		ngate = gate/delta + 1;
		/* now make the stack only if the traces are within
		 * the correct distance range */
		if(dist >= MIN(tvlw,tvhg) && dist <= MAX(tvlw,tvhg)){
			/* assume that Origin is set and that we
			 * give the offset with respect to the 
			 * */
			tpred = t0 + p*dist;
			ilw = (tpred - tss)/delta;
			iup = ilw + ngate;

		/* for this to work should somehow normalize by peak amplitude in the window */
			peak = 0.0;
			for(j=0; j < ngate; j++){
				l = j+ilw;
				if(l >= 0 && l < npts)
					tmp = ABS(sacdata[k].sac_data[l]);
					if(tmp > peak )
						peak = tmp;
			}
			if(peak > 0.0){
				for(j=0; j < ngate; j++){
					l = j+ilw;
					if(l >= 0 && l < npts)
						s[j] += sacdata[k].sac_data[l]/peak;
				}
			}
			nstk++;
		}
	}
	/* now systematically get the cross-correlation and look 
	 * for the maximum - this will be done by brute force since the 
	 * series are small 
	 * */
	for ( k=0 ; k < ntrc ; k++){
		npts  = sacdata[k].sachdr.ihdr[H_NPTS];
		delta = sacdata[k].sachdr.rhdr[H_DELTA];
		dist  = sacdata[k].sachdr.rhdr[H_DIST ];
		if(dist >= MIN(xx1,xx2) && dist <= MAX(xx1,xx2)){
			/* we know that O is = 0 = reference time */
			/* get offset with respect to the reference trace */
			tpred = t0 + p*dist;
			ilw = (tpred - tss)/delta;
			iup = ilw + ngate;
			/* offset is shift of trace relative to the stack in range +- ngate samples */
			offset = do_getoffset(k,ilw,iup,ngate,npts,s,sacdata[k].sac_data);
			if(isp == YES)
				sacdata[k].sachdr.rhdr[H_A]  = tpred + offset * delta ;
			else
				sacdata[k].sachdr.rhdr[H_T0] = tpred + offset * delta ;
		}
	}
}

int do_getoffset(int k,int ilw,int iup,int ngate,int npts,float *s,float *x)
{
	float sum, summax;
	int i, j, l;
	int offset ;
	summax = 0.0;
	offset = 0;
	/* time shift */
	for(i= -ngate/2; i < ngate/2; i++){
		sum = 0.0;
		for(j=0; j < ngate; j++){
			l = j + ilw +i;
			if(l >= 0 && l < npts)
				sum = s[j] * x[l];
			if(sum > summax){
				offset = i;
				summax = sum;
			}
		}
	}
	return (offset);
}

int dops_refr(float x0,float y0,float xlen,float ylen,double tss,double tes,float hvmin,float hvmax,char *msg,char *phase)
{
int rret;
float p, t0, xx1, xx2, tt1, tt2;
float tvmin, tvmax, ttmin, ttmax;
double ttss,ttes;
float thvmin, thvmax;
float xx, yy;
float gate;
int refrret;

	ttss = tss;
	ttes = tes;
	refrret = menu_sel(0.5,0.5,prref,sizeof(prref), msg);
	dorefrpick(phase,x0,y0,xlen,ylen,
        	refr_typeaxis, refr_dtdx,
        	tss,tes,hvmin,hvmax, &t0, &p, 
		&xx1, &tt1, &xx2, &tt2);
	rret = menu_sel(0.5,0.5,refrpsyn,sizeof(refrpsyn),
		"Accept Y/N Refine?"); 
	if(rret == MENU_REFR_PSY){
		printf("%s: Refr t0 %f p %f (sec/km) Vel %f (km/sec) Refractor %d\n",phase,t0,p,1.0/p,refrret);
		fprintf(gsac_control.refrpick,
			"%s: Refr %d t0 %f p %f %f %f %f %f\n",phase,refrret,t0,p,xx1,tt1,xx2,tt2);
		fflush(gsac_control.refrpick);
		gsac_srefr[rret].set = rret;
		gsac_srefr[rret].p = p;
		gsac_srefr[rret].t0 = t0;
		gsac_srefr[rret].x1 = xx1;
		gsac_srefr[rret].t1 = tt1;
		gsac_srefr[rret].x2 = xx2;
		gsac_srefr[rret].t2 = tt2;
		do_setat0(p,t0,xx1,xx2,phase);
		return YES;
	} else if (rret == MENU_REFR_PSREFINE){
		printf("%s: t0 %f p %f\n",phase,t0,p);
	gsac_refrrs(p,refr_dtdd,refr_absolute,YES,refr_doptaud,
		refr_tmin, refr_tmax, refr_vmin, refr_vmax, x0, y0, xlen, ylen, refr_amp,
		refr_typeaxis,refr_doshd,refr_doshdplmn,refr_ann,refr_annstr,
		YES, " ", refr_tstr, refr_tstrl, refr_whichaxis, refr_shdcolor,
		&ttss,&ttes,&thvmin,&thvmax,YES,dovlim,NO,refr_doclip,refr_cliplev,YES);
		/* force a line to indicate the reference line */
		xx = x0 ; yy = y0 + ylen*(t0+p*thvmin-ttss)/(ttes - ttss);
		plot (xx,yy,3);
		newpen(3);
		xx = x0 + xlen ; 
		plot (xx,yy,2);
		newpen(1);
		/* save previous values */
		tvmin = refr_vmin ;
		tvmax = refr_vmax ;
		ttmin = refr_tmin ;
		ttmax = refr_tmax ;
		rret = NO;
		gmesg("Select correlation gate");
		while(rret == NO)
			rret = do_zoomm(x0,y0,xlen,ylen,ttss,ttes,thvmin,thvmax);
		/* restore previous values */
		gate = ABS(refr_tmax - refr_tmin);
		if(strcmp(phase,"S") == 0){
		/* S */
			do_prefine(p,t0,xx1,tt1,xx2,tt2,gate,
				refr_tmin,refr_tmax,refr_vmin,refr_vmax,NO);
		} else if(strcmp(phase,"P") == 0) {
		/* P */
			do_prefine(p,t0,xx1,tt1,xx2,tt2,gate,
				refr_tmin,refr_tmax,refr_vmin,refr_vmax,YES);
		}
		refr_vmin = tvmin ;
		refr_vmax = tvmax ;
		refr_tmin = ttmin ;
		refr_tmax = ttmax ;
		return YES;
	} else {
		return YES;
	}
}

int dops_refl(float x0,float y0,float xlen,float ylen,double tss,double tes,float hvmin,float hvmax,char *msg,char *phase)
{
int rret;
float p, t0, xx1, xx2, tt1, tt2;
float tvmin, tvmax, ttmin, ttmax;
double ttss,ttes;
float thvmin, thvmax;
float xx, yy;
float gate;
int refrret;
int reflmult, refllyr;
	/* define the reflector */
	refllyr = menu_sel(0.5,0.5,prref,sizeof(prref), "Select Reflection");
	reflmult = menu_sel(0.5,0.5,prref,sizeof(prref), "Select Multiple");
	doreflpick(phase,x0,y0,xlen,ylen,
        	refr_typeaxis, refr_dtdx,
        	tss,tes,hvmin,hvmax, &t0, &p, 
		&xx1, &tt1, &xx2, &tt2);
	rret = menu_sel(0.5,0.5,refrpsyn,sizeof(refrpsyn),
		"Accept Y/N ?"); 
	if(rret == MENU_REFR_PSY){
		printf("%s: Refl t0 %f p %f (sec/km) Vrms %f (km/sec) Reflector %d Multiple %d\n",phase,t0,p,1.0/p,refllyr,reflmult);
		fprintf(gsac_control.refrpick,
			"%s: Refl %d t0 %f p %f Reflector %d Multiple %d %f %f %f %f\n",phase,refrret,t0,p,refllyr,reflmult,xx1,tt1,xx2,tt2);
		/*
			*/
		return YES;
	} else {
		return YES;
	}
}

void do_setat0(float p,float t0,float xx1,float xx2,char *phase)
{
/* when the straight line technique is used, set the arrival time values */
int k, ntrc;
float dist;
float tpred;

	ntrc = gsac_control.number_itraces;
	/* get the delta and the number of points from the
	 * first trace - this should all be
	 * the same */
	if(ntrc < 1 )
		return;
	for ( k=0 ; k < ntrc ; k++){
		dist  = sacdata[k].sachdr.rhdr[H_DIST ];
		if(dist >= MIN(xx1,xx2) && dist <= MAX(xx1,xx2)){
			/* we know that O is = 0 = reference time */
			/* get offset with respect to the reference trace */
			tpred = t0 + p*dist;
			if(strcmp(phase,"P") == 0)
				sacdata[k].sachdr.rhdr[H_A]  = tpred ;
			else if(strcmp(phase,"S") == 0)
				sacdata[k].sachdr.rhdr[H_T0] = tpred ;
		}
	}
}

static int timelist[] = { H_A, H_T0, -1};
void  gsac_refrrs(float rs_dtdx,float rs_dtdd,int rs_absolute,int rs_doptaux,
	int rs_doptaud,float rs_tmin,float rs_tmax,float rs_vmin,float rs_vmax,
	float x0, float y0, float xlen, float ylen, float rs_amp, int rs_typeaxis,
	int rs_doshd,int rs_doshdplmn,int rs_ann,char *rs_annstr,int rs_exttitle,
	char *rs_titstr, char  *rs_tstr, int rs_tstrl, int rs_whichaxis, int rs_shdcolor,
	double *ttss,double *ttes,float *hhvmin,float *hhvmax,int dotlim,int dovlim,int labx,
	int rs_doclip,float rs_cliplev,int rs_shwpick){
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
	 * double ttss       - starting time
	 * double ttes	    - ending time
	 * float hhvmin      - minimum axis value
	 * float hhvmax      - maximum axis value
	 * int labx          - if YES label Dist (e.g.) axis and annotate with Origin, Ray parameter
	 * int rs_doclip     - YES clip trace according to xlen/ntrc multiple 
	 * float rs_cliplev  - default 2
	 * int rs_shwpick    - YES show A t0 picks only
	 * */

	int ntrc;
	float dy, dv;
	float hvmin, hvmax;
	float v, uu, vv, xx, yy;
	int npts;
	int k, i, j, ls, kk, nii;
	float trmax;	/* trace maximum amplitude */
	float tv;	/* normalized trace value [-1,1] */
	double ts,te,twin;
	double tss,tes;
	float tb, to;
	float dclip;
	float av;	/* plot amplitude of trace in page units */

	char timstr[32];
	char rs_strann[100];
	int rs_oset ;
	char ostr[100];
	char lstr[100];
	double otss,otes;
	float vlen ;  /* length of vaxis , e.g., not time axis */

	ntrc = gsac_control.number_itraces;
	/* define the timing window for the plot */
	if(rs_absolute == YES) {
		/* instead of using preset begin and end we must recalcuate
		 * because of the possibility of a p-tau plot - note
		 * we only do a p-tau for traces which have a defined distance 
		 * and we only plot those traces */
		/* if O is set define tss and tes */
		otss =  1.0e+30;
		otes = -1.0e+30;
		twin = 0.0;
		rs_oset = NO;
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
		for ( k=0 ; k < ntrc ; k++){
			ts = sacdata[k].sachdr.rhdr[H_B];
			te = sacdata[k].sachdr.rhdr[H_E];
			if((te-ts) > twin )
				twin = te - ts ;
			/* note that the LAST instance of ts te is preserved
			 * for the time scale */
			if(k==0)
				tss = ts;
			tes = tss + twin;

		}

	}
	if(dotlim==YES ){
		tss = MIN(rs_tmin,rs_tmax);
		tes = MAX(rs_tmin,rs_tmax);
		twin = tes - tss;
	}
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
		strcpy(rs_tstr,"T - X/");
		strcat(rs_tstr,refr_vstr);
		strcat(rs_tstr, " (s)");
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
	}
		if(gsac_control.plotdevice == WIN) XviG_Flush();
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

	for ( kk=0 ; kk < ntrc  ; kk++){
		k = sortptr[kk];
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
		if(k == 0)
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
			av = rs_amp*tv;
			if(rs_doclip){
				av = SIGN(av)*MIN(rs_cliplev*dclip,ABS(av));
			}
			if(rs_typeaxis == LANDSCAPE){
				xx = x0 + uu*xlen - av;
				yy = y0 + vv*ylen;
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
			}
			newpen(rs_shdcolor);
			if(nii > 0){
				plot(px[0],py[0],3);
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
			for(i = 0 ; i < nii ; i ++){
				plot(px[i],py[i],2);
			}
			plot(px[nii-1],py[nii-1],3);
		}
		gsac_setcolor(NO , k, ntrc);
		/* annotate the trace with picks */
		if(rs_shwpick == YES){
			newpen(2);
			for(j=0 ; timelist[j] >= 0 ; j++){
				for ( k=0 ; k < ntrc  ; k++){
					v = sacdata[k].sachdr.rhdr[rs_whichaxis] ;
					uu = (v - hvmin)/(hvmax - hvmin);
					tb = sacdata[k].sachdr.rhdr[H_B];
					te = sacdata[k].sachdr.rhdr[H_E];
					to = sacdata[k].sachdr.rhdr[H_O];
					dv = sacdata[k].sachdr.rhdr[H_DELTA]/twin;
					tv = sacdata[k].sachdr.rhdr[timelist[j]];
					if(tv != -12345) {
					if(k == 0)
						ts = tb;
					if(rs_absolute == YES) {
						vv = (float)(sacdata[k].tzref + tv - gsac_control.begmin -tss)/twin;
						if(tv != -12345.){
							vv = (float)(tv - to - tss)/twin ;
						}
						if(rs_doptaux){
							vv = vv - rs_dtdx*sacdata[k].sachdr.rhdr[H_DIST]/twin ;
						} else if(rs_doptaud){
							vv = vv - rs_dtdd*sacdata[k].sachdr.rhdr[H_GCARC]/twin ;
						}
					}
					if(rs_typeaxis == LANDSCAPE){
						xx = x0 + uu*xlen ;
						yy = y0 + vv*ylen;
					}
					switch(timelist[j]){
						case H_A:
							newpen(2);
							plot(xx-0.05,yy,3);plot(xx+0.05,yy,2);
							break;
						case H_T0:
							newpen(4);
							plot(xx-0.05,yy,3);plot(xx+0.05,yy,2);
							break;
					}
					}
				}
			}
			newpen(1);
		}
		/* annotate the trace if desired, but only if the trace 
			is plotted  */
		if(rs_ann == YES && uu >= 0.0 && uu <= 1.0){
			/* we will put this outside the axis - we currently
			 * do this only is the annotation is STA we can later
			 * permit date  or filename - we use rs_strann
			 * so that we can later use other things here
			 * */
			if(strncmp(rs_annstr,"STA",3)==0){
				/* annotate */
				strcpy(rs_strann,sacdata[k].schdr[H_KSTNM]);
			} else if(strncmp(rs_annstr,"NAME",4)==0){
				/* annotate */
				ls = strlen(sacdata[k].sac_ifile_name);
				if(ls > sizeof(rs_strann)-1)ls--;
				strncpy(rs_strann,sacdata[k].sac_ifile_name,ls);
				rs_strann[ls]='\0';
			}
			gclip("off", x0, y0, x0+xlen, y0+dy);
			if(rs_typeaxis == LANDSCAPE){
				xx = x0 + uu*xlen ;
				yy = y0 + ylen + 0.2;
				gleft(xx,yy,0.1,rs_strann,90.0);
			}
			/* now determine where to put the string */
			gmesg(" ");
			gclip("on", x0, y0, x0+xlen, y0+dy);
		}
	}
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

int doreflpick(char *phase,float x0,float y0,float xlen,float ylen,
	int refr_typeaxis, float refr_dtdx,
	double tss,double tes,float hvmin,float hvmax,float *t0,float *p,
	float *xx1, float *tt1, float *xx2, float *tt2)
{
float ax, ay;
char c[2];
float px, py;
int HasMouse; 
float XminDev, YminDev, 
	XmaxDev, YmaxDev, XminClip, 
	YminClip, XmaxClip, YmaxClip;
int Color;
float vel;
	gmesg("Define reflection hyperbola");
	ginfo(&HasMouse, &XminDev, &YminDev, 
		&XmaxDev, &YmaxDev, &XminClip, 
		&YminClip, &XmaxClip, &YmaxClip,&Color);
	px =  (0.0 - hvmin)/(hvmax - hvmin) ;
	py =  (0.0 - tss  )/(tes - tss) ;
	gcontrol(8,px,py, 0.0, 0.0);

	gclip("on",x0,y0,x0+xlen,y0+ylen);
	gcursor("Hyperbola");
	curaxy(&ax, &ay, c);
	/* map to user coordinates, then plot */
	px = (ay - y0)/(ylen);
	*t0 = px*(tes - tss);
	if(ay > 0.0){
		px = ( ax - x0)/(xlen);
		vel = px*(hvmax - hvmin) / (*t0);
		*p = 1.0/vel ;
	} else {
		vel = 1.0;
		*p = 1.0/vel ;
	}
	gclip("off",x0,y0,x0+xlen,y0+ylen);
}

int dorefrpick(char *phase,float x0,float y0,float xlen,float ylen,
	int refr_typeaxis, float refr_dtdx,
	double tss,double tes,float hvmin,float hvmax,float *t0,float* p,
	float *xx1, float *tt1, float *xx2, float *tt2)
{
/* t = t0 + p x */
float cx1, cy1,  cx2, cy2;	/* screen coordinates */
float x1, t1, x2, t2;		/* user x-t coordinates */
float px1,py1,px2,py2;			/* fraction of window space */
char c[2];
int cmd;
	/* Put up menu to define the refraction arrival 1 = direct 
	 * pick end points, remap to distance time
	 * define intercept slope, adjust slope for reduction velocity
	 * if acceptable , write to gsac.control.refrpick and then plot the
	 * phase
	 * */
	gmesg("Select two points for the refraction");
	gclip("on",x0,y0,x0+xlen,y0+ylen);
	gcursor("Cross");
	for(; ;){
		curaxy(&cx1, &cy1, c);
		cmd = -1;
		if(inside(cx1,cy1,x0,y0,x0+xlen,y0+ylen)){
			px1 = (cx1 - x0)/xlen;
			py1 = (cy1 - y0)/ylen;
			x1 = hvmin + (hvmax - hvmin)*px1;
			t1 = tss + (tes - tss)*py1 + refr_dtdx*ABS(x1);
			gcursor("Rubber");
			while(cmd == -1){
				curaxy(&cx2, &cy2, c);
				if(inside(cx2,cy2,x0,y0,x0+xlen,y0+ylen)){
				if(cx1 != cx2 && cy1 != cy2){
					px2 = (cx2 - x0)/xlen;
					py2 = (cy2 - y0)/ylen;
					x2 = hvmin + (hvmax - hvmin)*px2;
					t2 = tss + (tes - tss)*py2 +refr_dtdx*ABS(x2);
					*p = (t2 - t1 )/(x2 - x1) ;
					*t0 = t1 - (*p)*x1;
					*xx1 = x1;
					*tt1 = t1;
					*xx2 = x2;
					*tt2 = t2;
/*
printf("x-t: (%f,%f) (%f,%f)\n",x1,t1,x2,t2);
*/
					gcursor("Arrow");
				return;
				}
				}
			}
			return;
		}
	}
	gclip("off",x0,y0,x0+xlen,y0+ylen);
}

void do_hardcopy_refr(void)
{
	/* make a hard copy by switching from
		interactive graphics to plot
		mode. Then when completed, switch
		back to interactive mode
	*/
	char fname[12];
	char ostr[100];
	FILE *refrctlfid ;
	float x0, y0, xlen, ylen;
	double tss,tes;
	float hvmin, hvmax;
	float xx, yy;

        xlen = gsac_control.xlen ;
        ylen = gsac_control.ylen ;
        x0   = gsac_control.x0 ;
        y0   = gsac_control.y0 + 0.8 ;


	sprintf(refr_pltname,"REFR%3.3d.PLT",gsac_control.plotcount_refr);
	sprintf(fname  ,"REFR%3.3d.CTL",gsac_control.plotcount_refr);
	gsac_control.plotcount_refr++;
	/* start plot mode */
	ginitf(refr_pltname,"GSAC");
	/* annotate before the plot to preserve the state */
	xx = x0; yy = y0 - 0.8;
	gsac_control.hold = YES;
	gsac_control.plotdevice = PLT;
	gsac_control.everinteractive = NO ;
	/* plot the traces */
	gsac_refrrs(refr_dtdx,refr_dtdd,refr_absolute,refr_doptaux,
		refr_doptaud,refr_tmin,refr_tmax,refr_vmin,
		refr_vmax, refrs_x0, refrs_y0, refrs_xlen, refrs_ylen, refr_amp,
		refr_typeaxis,refr_doshd,refr_doshdplmn,refr_ann,
		refr_annstr, YES, " ", refr_tstr, refr_tstrl, 
		refr_whichaxis, refr_shdcolor, &tss, &tes, 
		&hvmin,&hvmax,dotlim,dovlim,NO,refr_doclip,
		refr_cliplev, YES);
	/* annotate plot with processing commands */
	gclip("off", x0, y0, x0+xlen, y0+ylen);
	newpen(1);
	if(refr_hpc > 0.0){
		sprintf(ostr,"hp c %f n 2 p 2",refr_hpc);
		do_outstr(x0,y0,&xx,&yy,ostr);
	}

	if(refr_lpc > 0.0){
		sprintf(ostr,"lp c %f n 2 p 2",refr_lpc);
		do_outstr(x0,y0,&xx,&yy,ostr);
	}
	if(refr_doshd ==YES){
		if(refr_doshdplmn == 0)
			do_outstr(x0,y0,&xx,&yy,"SHADE POSITIVE");
		else if(refr_doshdplmn == 1)
			do_outstr(x0,y0,&xx,&yy,"SHADE NEGATIVE");
		else if(refr_doshdplmn == 2)
			do_outstr(x0,y0,&xx,&yy,"SHADE ALL");
	}
	if(refr_doclip){
		sprintf(ostr,"CLIPLEVEL %5.1f",refr_cliplev);
		do_outstr(x0,y0,&xx,&yy,ostr);
	}
	if(agc_window_pct > 0.0){
		sprintf(ostr,"AGCWINDOW %5.0f (%%)",agc_window_pct);
		do_outstr(x0,y0,&xx,&yy,ostr);
	}

	/* save processing commands on control file */
	refrctlfid=fopen(fname,"w+");
	fprintf(refrctlfid,
		"refmod96 -XMIN %f -XMAX %f -TMIN %f -TMAX %f -XLEN %f -YLEN %f -X0 %f -Y0 %f -KF 2 -KR 4 -NO -NM 1 -P -SH  -M model",
		hvmin, hvmax, tss, tes,
		gsac_control.xlen, gsac_control.ylen, 
		gsac_control.x0, gsac_control.y0 +0.8
		);
		if(refr_doptaux)
			fprintf(refrctlfid," -VRED %f",1.0/refr_dtdx);
		fprintf(refrctlfid,"\n");
	fclose(refrctlfid);
	
	/* return to interactive graphics mode */
	gsac_control.plotdevice = WIN;
	gsac_control.hold = NO;
	gsac_control.everinteractive = YES ;
	ginitf("INTEM","GSAC");
}	

void do_outstr(float x0,float y0,float *xx,float *yy,char *str)
{
	float x,y;
	x = *xx;
	y = *yy;
	gleft(x,y,0.1,str,0.0) ;
	/* get new position */
	y -= 0.20;
	if(y < 0.25){
		x += 2.0; y = y0 - 0.8;
	}
	*xx = x;
	*yy = y;
}
