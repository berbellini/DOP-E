/* NOTES
 * before ploting traces in a group determine the extremes of the absolute
 * B and E so that the time scale can be defined
 *
 * When doing ppk with perplot - the overscale scale will be determined
 * by perplot. Note however that the clip region may be smaller if
 * on the last page
 */
/* Changes:
	26 JAN 2004 - update phase list for teleseisms - also care taken that phase ID
	differs between Regional and Teleseism
	16 JUN 2005 - corrected Pp type to pP
	30 AUG 2006 - line 627 To get the ABS() to work correctly I had to
		rewrite the definition of ABS in gsac.h
        01 JUL 2008 - added the F marker to ppk and changed gsac.h
	30 SEP 2009 - added an Undo to the regional/teleseismic phase pick menu
			which erases the last pick
	11 JUL 2010 - added the ~ command which then creates a screen dump 
		with a file name DUMPxxx.PLT  The purpose is to get a hardcopy
		of the complete interactive menu for documentation
        20 SEP 2012 - the use of the 'm' key has the same effect as the '*' key
*/
#include	<stdio.h>
#include	"gsac.h"

#include "gsac_plot.h"
#include "gsac_sac.h"
#include "gsac_arg.h"
#include "csstim.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

extern void gsac_exec_doxlim(void);

#define	PK_PERPLOT	0
#define	PK_PERSTA	1
#define	PK_ABSOLUTE	2
#define	PK_RELATIVE	3
#define	PK_MARKALL	4
#define	PK_MARKALLOFF	5
#define	PK_DEFAULT	6
#define	PK_REGIONAL	7
#define	PK_TELESEISM	8
#define	PK_QUALITY	9
#define	PK_PQUALITY	10


struct arghdr pkarg[] = {
	{PK_PERPLOT,  "PERPLOT"	, NHDR, 0, 1, NO, "PERPLOT [n|OFF]",-1},
	{PK_PERSTA,   "PERSTA"	, YHDR, 0, 1, NO, "PERSTA [ON|OFF] ",-1},
	{PK_ABSOLUTE, "ABSOLUTE", RHDR, 0, 0, NO, "ABSOLUTE", 1},
	{PK_RELATIVE, "RELATIVE", RHDR, 0, 0, NO, "RELATIVE", 3},
	{PK_ABSOLUTE, "A", RHDR, 0, 0, NO, "ABSOLUTE",-1},
	{PK_RELATIVE, "R", RHDR, 0, 0, NO, "RELATIVE",-1},
	{PK_MARKALL, "MARKALL", RHDR, 0, 0, NO, "MARKALL",-1},
	{PK_MARKALLOFF, "MARKALLOFF", RHDR, 0, 0, NO, "MARKALLOFF",-1},
	{PK_DEFAULT, "DEFAULT", RHDR, 0, 0, NO, "DEFAULT", 1},
	{PK_REGIONAL, "REGIONAL", RHDR, 0, 0, NO, "REGIONAL",3},
	{PK_TELESEISM, "TELESEISM", RHDR, 0, 0, NO, "TELESEISM", 1},
	{PK_QUALITY, "QUALITY", RHDR, 0, 0, NO, "QUALITY", 1},
	{PK_PQUALITY, "PQUALITY", RHDR, 0, 0, NO, "PQUALITY", 2},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};


/* these are temporary variables only used here */
float pk_real[10];
int   pk_int [10];
int   pk_yn;
int   pk_num;

static int pkperplot = -1;
static int pkpersta = NO;
static int pkabsolute = YES;
static int pkmarkall = NO;
static int pkphasemenu = PK_DEFAULT ;
static int pickundo = NO;  /* this will be a binary flag for a single use */
static int last_h = -1 ; /* parameters for undo */
static int last_kh = -1;
static int last_k = -1;

/* commands for interactive picking */
#define NEXT 0
#define QUIT 1
#define CONT 2
#define PREVIOUS 3
extern struct plt_ctl *pmap ;

extern void gsac_show_plotpk(float x0, float y0, float xlen, float ylen, int ns, int ne, float ts, float te, float dy, int ntrc, double twin,int numperframe, int pabsolute, float *uy, float *ly, int setup, int ylimctrl, float ylimlow , float ylimhigh, int overlay);
int gsac_show_intpk(float x0, float y0, float xlen, float ylen, int ns, int ne, float ts, float te, float dy, int ntrc, double twin,int numperframe, float uy, float ly);
void traceinfo(float cx,float cy,float x0,float y0, float dy,float xlen,float ylen,int ns, int ne, int ntrc, int *k, float *ux, float *uy , int numperframe);
int inside(float cx, float cy, float xl, float yl, float xh, float yh);
void updatehdr(int k,int switchchar,double tp,int ns,int ne,int numperframe);
void updthd(int k,int tn, int ktn, char *kstr ,double tp,int ns,int ne,int numperframe);
void showdec(float xl, float yl, float xh, float yh, int inc);
void upd_hdr_time(int k, float tph, int H, int HK, char *str);
void upd_hdr_ihdr(int k, int val,  int H );
void do_check(int val, int k);
void show_check(float x0, float y0, float xlen, float ylen, int ns, int ne, float dy, int ntrc, int numperframe, float yh, float yl);
static char dump_pltname[12];

/* menu routines for regional and teleseism */
#include "nmenu.h"
void clearregion(float xl, float yl, float xh, float yh);
void show_menu (float x0, float y0, struct menu *m, int size, int *nm);
int inside(float xv, float yv, float xlb, 
	float ylb, float xhb, float yhb);

int pk_parse_key(float x0, float y0, float xlen, float ylen, int ns, int ne, float ts, float te, float dy, int ntrc, double twin,int numperframe, float yh, float yl, float cx,float cy,int switchchar,int *lpos);

/* define the phase menus for regional and teleseismic phases */
/* NOTE UNIQUE NUMBERS ARE USED FOR THE REGIONAL AND TELESEISMIC PHASES */

#define MENU_P_UNDO 200
#define MENU_PREG_P  1
#define MENU_PREG_S  2
#define MENU_PREG_Lg 3
#define MENU_PREG_Pg 4
static struct menu preg[] = {
		{  -1.0, -1.0, -1.0, -1.0, "P\0" , MENU_PREG_P , -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "S\0" , MENU_PREG_S , -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Lg\0" , MENU_PREG_Lg, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Pg\0" , MENU_PREG_Pg, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Undo\0" , MENU_P_UNDO, -1, 1, 1},
};

#define MENU_PTEL_P  101
#define MENU_PTEL_S  102
#define MENU_PTEL_PKP  103
#define MENU_PTEL_PcP  104
#define MENU_PTEL_PP  105
#define MENU_PTEL_pP  106
static struct menu ptel[] = {
		{  -1.0, -1.0, -1.0, -1.0, "P\0" , MENU_PTEL_P , -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "S\0" , MENU_PTEL_S , -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "PKP\0" , MENU_PTEL_PKP, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "PcP\0" , MENU_PTEL_PcP, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "PP\0" , MENU_PTEL_PP, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "pP\0" , MENU_PTEL_pP, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Undo\0" , MENU_P_UNDO, -1, 1, 1},
};

/* define the phase menu for phase quality */

#define MENU_PQUAL_I 1
#define MENU_PQUAL_E 2
static struct menu pqual[] = {
		{  -1.0, -1.0, -1.0, -1.0, "i\0" , MENU_PQUAL_I, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "e\0" , MENU_PQUAL_E, -1, 1, 1},
};


/* define the phase menu for phase polarity */
#define MENU_PPOL_C 1
#define MENU_PPOL_p 2
#define MENU_PPOL_D 3
#define MENU_PPOL_m 4
#define MENU_PPOL_X 5
static struct menu ppol[] = {
		{  -1.0, -1.0, -1.0, -1.0, "C\0" , MENU_PPOL_C, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "+\0" , MENU_PPOL_p, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "D\0" , MENU_PPOL_D, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "-\0" , MENU_PPOL_m, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "X\0" , MENU_PPOL_X, -1, 1, 1}
};

/* Qualioy menu items to set  state */
#define MENU_QUALITY_ACCEPT  201
#define MENU_QUALITY_REJECT  202
static struct menu pquality[] = {
		{  -1.0, -1.0, -1.0, -1.0, "Accept\0" , MENU_QUALITY_ACCEPT, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Reject\0" , MENU_QUALITY_REJECT, -1, 1, 1},
};
static int pkquality ;

void show_polmenu(void);
void show_qualmenu(void);
void proc_phmenu(float cx,float cy,float x0,float y0,float dy,
	float xlen,float ylen,int ns,int ne, float yh, float yl,
	int ntrc, float ts, float te, double twin, int numperframe, 
	struct menu *p);
void proc_qlmenu(float cx,float cy,float x0,float y0,float dy,
	float xlen,float ylen,int ns,int ne, float yh, float yl,
	int ntrc, float ts, float te, double twin, int numperframe, 
	struct menu *p);
void proc_polmenu(char *ostr,float x0,float y0,float xlen,float ylen);
void proc_qualmenu(char *ostr,float x0,float y0,float xlen,float ylen);



extern void XviG_Flush();

char ostr1[80];
char ostr2[80];

int nmp;	/* number of phase entries for regional/teleseismic */
int npl;	/* number of entries in polarity meny */
int nqa;	/* number of entries in quality meny */


void gsac_set_param_plotpk(int ncmd, char **cmdstr)
{
	/* the only parameter to be set is MORE 
	 *
	 */
int i;
int HasMouse; 
float XminDev, YminDev, 
	XmaxDev, YmaxDev, XminClip, 
	YminClip, XmaxClip, YmaxClip;
int Color;
	/* initialize interactive  graphics */
	if(gsac_control.plotinit == NO){
		if(gsac_control.plotdevice==WIN){
			ginitf("INTEM","GSAC");
			gmesg("Initializing Interactive Graphics");
			gsac_control.everinteractive = YES;
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
		gsac_control.plotinit = YES;
		gsac_control.plotchange = NO;
	} else {
		gframe(2);
		if(gsac_control.plotchange == YES){
			if(gsac_control.plotdevice==WIN){
				ginitf("INTEM","GSAC");
				gmesg("Initializing Interactive Graphics");
				gsac_control.everinteractive = YES;
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
	}
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, pkarg, NO, YES))
		return;
	/* set defaults */
	pkquality = MENU_QUALITY_ACCEPT ;
	/* parse commands */
	for(i=0 ; pkarg[i].key[0] != '\0' ; i++){
		if(pkarg[i].used > 0){
			if(pkarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, pkarg[i].key, 
					pkarg[i].mfit,pkarg[i].narg, pk_real);
			} else if(pkarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, pkarg[i].key, 
					pkarg[i].mfit,pkarg[i].narg, pk_int );
			} else if(pkarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, pkarg[i].key, 
					pkarg[i].mfit,pkarg[i].narg, &pk_yn );
			} else if(pkarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, pkarg[i].key, 
					pkarg[i].mfit,pkarg[i].narg, &pk_num );
			}
			switch(pkarg[i].id){
				case PK_PERPLOT:
					pkperplot = pk_num;
					break;
				case PK_ABSOLUTE:
					pkabsolute = YES;
					break;
				case PK_RELATIVE:
					pkabsolute = NO;
					break;
				case PK_PERSTA:
					if(pk_yn == NO)
						pkpersta = NO;
					else if(pk_yn == YES)
						pkpersta = YES;
					break;
				case PK_MARKALL:
					pkmarkall = YES;
					break;
				case PK_MARKALLOFF:
					pkmarkall = NO;
					break;
				case PK_DEFAULT:
					pkmarkall = NO;
					pkpersta = NO;
					pkphasemenu = PK_DEFAULT;
					break;
				case PK_REGIONAL:
					pkphasemenu = PK_REGIONAL;
					break;
				case PK_TELESEISM:
					pkphasemenu = PK_TELESEISM;
					break;
				case PK_QUALITY:
					pkphasemenu = PK_QUALITY;
					break;
				case PK_PQUALITY:
					pkphasemenu = PK_PQUALITY;
					break;

			}
		}
	}

}

void gsac_exec_plotpk(void)
{
	int iret,k, kkk;
	int ntrc;
	float x0, y0, xlen, ylen, dy;

	int numperframe;

	double ts, te;
	double twin;

	float uy, ly;

	/* if there are no traces return */
	ntrc = gsac_control.number_otraces;
	if(ntrc < 1)
		return;
	gsac_exec_doxlim();
	/* must have interactive graphics */
	if(gsac_control.plotdevice!=WIN){
		printf("PLOTPK requires interactive graphics\n");
		return;
	}
	if(pmap == (struct plt_ctl *)NULL)
		pmap = (struct plt_ctl *)calloc(ntrc, sizeof(struct plt_ctl));
	else
		pmap = (struct plt_ctl *) realloc(pmap,ntrc*sizeof(struct plt_ctl));
	gsac_control.uxmul = 1.0;
	gsac_control.uxcen = 0.5;
	/* beginning of trace plot 
	 * For absolute time plot this is just 
	 * 	gsac_control.tzend - gsac_control.tzbeg
	 * For relative, we need the maximum time window of
	 * 	a trace
	 * Note then when we build in plot limits this will change
	 * But here we introduce the concept time for the X-axis */
	if(pkabsolute == YES) {
		twin = (gsac_control.endmaxx - gsac_control.begminx);
		ts = 0.0;
		te = twin;
	} else {
		twin = 0.0;
		for ( k=0 ; k < ntrc ; k++){
			ts = sacdata[k].tzbegx - sacdata[k].tzref ;
			te = sacdata[k].tzendx - sacdata[k].tzref ;
			if((te-ts) > twin )
				twin = te - ts ;
		}

	}

	if(pkperplot > 0)
		if(pkperplot > ntrc)
			numperframe = ntrc;
		else
			numperframe = pkperplot;
	else
		numperframe = ntrc;

	/* we must define the time window for the plot */

	xlen = 8.0;
	ylen = 6.0;
	x0 = 1.5;
	y0 = 1.5;
	dy = ylen / numperframe;
	gclip("off", x0, y0, x0+xlen, y0+ylen);

	/* do the bottom axis */

	iret = NEXT;
	for ( kkk=0 ; kkk < ntrc && iret == NEXT ; kkk+=numperframe){
		/* put up the traces */
		gsac_control.uxmul = 1.0;
		gsac_control.uxcen = 0.5;
		gframe(2);
		gsac_show_plotpk( x0,  y0,  xlen,  ylen,  kkk,  
				MIN(kkk+numperframe,ntrc),  ts,  te,  dy, ntrc, 
				twin, numperframe, pkabsolute, &uy, &ly, 
				YES, YLIM_OFF, -1, +1, NO);
		/* do the interactive picking */
		iret = gsac_show_intpk(x0, y0, xlen, ylen, kkk, 
				MIN(kkk+numperframe,ntrc), ts, te, dy, ntrc, 
				twin, numperframe,uy,ly);
		/* special case to handle previous */
		if(iret == PREVIOUS){
			kkk = MAX(-numperframe, kkk-numperframe-numperframe);
			iret = NEXT ;
		}
	}
	/* clean up */

	gclip("off", x0, y0, x0+xlen, y0+ylen);
	gcursor("Arrow");
}

int gsac_show_intpk(float x0, float y0, float xlen, float ylen, int ns, int ne, float ts, float te, float dy, int ntrc, double twin,int numperframe, float yh, float yl)
{
float cx, cy;
char ch[2];
int switchchar;
int doloop;
int iret;
int lpos;

/* if have special phase menus put them up now */
	if(pkphasemenu == PK_REGIONAL) {
		show_menu(1.0, 0.7, preg, sizeof(preg),&nmp);
	} else if(pkphasemenu == PK_TELESEISM) {
		show_menu(1.0, 0.7, ptel, sizeof(ptel),&nmp);
	} else if(pkphasemenu == PK_QUALITY  ){
		/* display the quality picks in memory */
		show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
		show_menu(1.0, 0.7, pquality, sizeof(pquality),&nmp);
		if(pkquality == MENU_QUALITY_ACCEPT)
			gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
		"Mouse Click to Accept (also a - Accept, r - Reject, n - Next, b - Back)",0.0);
		else if(pkquality == MENU_QUALITY_REJECT)
			gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
		"Mouse Click to Reject (also a - Accept, r - Reject, n - Next, b - Back)",0.0);
	} else if(pkphasemenu == PK_PQUALITY  ){
		/* display the quality picks in memory */
		show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
		show_menu(1.0, 0.7, pquality, sizeof(pquality),&nmp);
		if(pkquality == MENU_QUALITY_ACCEPT)
			gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
		"Mouse Click to Accept/Set P (also a - Accept, r - Reject, P - pick P, n - Next, b - Back)",0.0);
		else if(pkquality == MENU_QUALITY_REJECT)
			gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
		"Mouse Click to Reject (also a - Accept, r - Reject, n - Next, b - Back)",0.0);
	}
/* work interactively with the current traces */
	/* set clip region to the actual plot region not the entire box */
	/* cross hair cursors within the trace plot region , default
	 * 	arrow outside this region */
	gcursor("Cross");
	gclip("on", x0, yl , x0+xlen, yh);
	doloop = YES;
	/* infinite loop */
	ch[1] = '\0';
	lpos = 0;
	while(doloop){
		gsac_control.ylim_rsc = YES;
		curaxy(&cx,&cy,ch);
	/* we must decide if the cursor was in the button region
	 * or if the cursor was in the trace region */
		switchchar = ch[0];
		switchchar = toupper(switchchar);
		if(inside(cx,cy,x0,y0,x0+xlen,y0+ylen)){
			iret = pk_parse_key(x0, y0, xlen, ylen, ns, ne, 
				ts, te, dy, ntrc, twin,numperframe, yh, yl,
				cx,cy,switchchar,&lpos);
			if(iret == NEXT)
				return (NEXT);
			else if(iret == PREVIOUS)
				return (PREVIOUS);
			else if(iret == QUIT)
				return (QUIT);
		} else {
			if(pkphasemenu == PK_REGIONAL)
				proc_phmenu(cx,cy,x0,y0,dy,xlen,ylen,ns,ne,
					yh, yl, ntrc, ts, te, twin, numperframe,
					preg);
			if(pkphasemenu == PK_TELESEISM)
				proc_phmenu(cx,cy,x0,y0,dy,xlen,ylen,ns,ne,
					yh, yl, ntrc, ts, te, twin, numperframe,
					ptel);
			if(pkphasemenu == PK_QUALITY )
				proc_qlmenu(cx,cy,x0,y0,dy,xlen,ylen,ns,ne,
					yh, yl, ntrc, ts, te, twin, numperframe,
					pquality);
			if(pkphasemenu == PK_PQUALITY )
				proc_qlmenu(cx,cy,x0,y0,dy,xlen,ylen,ns,ne,
					yh, yl, ntrc, ts, te, twin, numperframe,
					pquality);

		}
	}
	gsac_control.ylim_rsc = NO;
	return(CONT);
}

int pk_parse_key(float x0, float y0, float xlen, float ylen, int ns, int ne, float ts, float te, float dy, int ntrc, double twin,int numperframe, float yh, float yl, float cx,float cy,int switchchar,int *lpos) 
{
/* remember that switchchar is upper case from the line before the call */
int k, kk;
char ch[2];
float bx, by;
float ux, uy;
float uv;
double tp;
double ttime;
int isxtime;
float amp;
float cx_X, twin_X;	/* for use with the X X to define window */
	traceinfo(cx,cy,x0,y0,dy,xlen,ylen,ns,ne,ntrc,&k,&ux,&uy ,numperframe);
	cx_X = -1;
	if(k >= 0)
	switch(switchchar){
		case '~':
			sprintf(dump_pltname,"DUMP%3.3d.PLT",
				gsac_control.plotcount_dump);
			gsac_control.plotcount_dump++;
			ginitf(dump_pltname,"GSAC");
        		gsac_control.hold = YES;
        		gsac_control.plotdevice = PLT;
        		gsac_control.everinteractive = NO ;
			gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,ntrc, 
				twin, numperframe, pkabsolute,&yh,&yl,NO, 
				YLIM_OFF, -1, +1, NO);
 			if(pkphasemenu == PK_QUALITY || pkphasemenu == PK_PQUALITY ){
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			}
			/* if have special phase menus put them up now */
			if(pkphasemenu == PK_REGIONAL) {
				show_menu(1.0, 0.7, preg, sizeof(preg),&nmp);
			} else if(pkphasemenu == PK_TELESEISM) {
				show_menu(1.0, 0.7, ptel, sizeof(ptel),&nmp);
			} else if(pkphasemenu == PK_QUALITY  ){
				/* display the quality picks in memory */
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
				show_menu(1.0, 0.7, pquality, sizeof(pquality),&nmp);
				if(pkquality == MENU_QUALITY_ACCEPT)
					gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
				"Mouse Click to Accept (also a - Accept, r - Reject, n - Next, b - Back)",0.0);
				else if(pkquality == MENU_QUALITY_REJECT)
					gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
				"Mouse Click to Reject (also a - Accept, r - Reject, n - Next, b - Back)",0.0);
			} else if(pkphasemenu == PK_PQUALITY  ){
				/* display the quality picks in memory */
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
				show_menu(1.0, 0.7, pquality, sizeof(pquality),&nmp);
				if(pkquality == MENU_QUALITY_ACCEPT)
					gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
				"Mouse Click to Accept/Set P (also a - Accept, r - Reject, P - pick P, n - Next, b - Back)",0.0);
				else if(pkquality == MENU_QUALITY_REJECT)
					gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
				"Mouse Click to Reject (also a - Accept, r - Reject, n - Next, b - Back)",0.0);
			}
			gsac_control.plotdevice = WIN;
			gsac_control.hold = NO;
			gsac_control.everinteractive = YES ;
			ginitf("INTEM","GSAC");
			break;
		case '-':
		case '_':
			/* compress time scale */
			gsac_control.uxcen += (ux - 0.5)/gsac_control.uxmul;
			gsac_control.uxmul /= 2.0 ;
			gsac_control.uxmul = MAX(1.0, gsac_control.uxmul/2.0);
			if(gsac_control.uxmul == 1.0)
				gsac_control.uxcen = 0.5;
			*lpos = 0;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			clearregion(0.0,1.0,10.0,8.0);
			gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,ntrc, 
				twin, numperframe, pkabsolute,&yh,&yl,NO, 
				YLIM_OFF, -1, +1, NO);
 			if(pkphasemenu == PK_QUALITY || pkphasemenu == PK_PQUALITY ){
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			}
			gclip("on", x0, yl , x0+xlen, yh);
			break;
		case '+':
		case '=':
			/* expand time scale */
			gsac_control.uxcen += (ux - 0.5)/gsac_control.uxmul;
			gsac_control.uxmul *= 2.0 ;
			*lpos = 0;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			clearregion(0.0,1.0,10.0,8.0);
			gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,ntrc, 
				twin, numperframe, pkabsolute,&yh,&yl,NO, 
				YLIM_OFF, -1, +1, NO);
 			if(pkphasemenu == PK_QUALITY || pkphasemenu == PK_PQUALITY ){
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			}
			gclip("on", x0, yl , x0+xlen, yh);
			break;
		case ' ':
			/* center cursor point */
			gsac_control.uxcen += (ux - 0.5)/gsac_control.uxmul;
			*lpos = 0;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			clearregion(0.0,1.0,10.0,8.0);
			gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,ntrc, 
				twin, numperframe, pkabsolute,&yh,&yl,NO, 
				YLIM_OFF, -1, +1, NO);
 			if(pkphasemenu == PK_QUALITY || pkphasemenu == PK_PQUALITY ){
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			}
			gclip("on", x0, yl , x0+xlen, yh);
			break;
		case '*':
		case 'M':
			/* amplify trace */
			pmap[k].uymul *= 2.0 ;
			*lpos = 0;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			clearregion(0.0,1.0,10.0,8.0);
			gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,ntrc, 
				twin, numperframe, pkabsolute,&yh,&yl,NO, 
				YLIM_OFF, -1, +1, NO);
 			if(pkphasemenu == PK_QUALITY || pkphasemenu == PK_PQUALITY ){
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			}
			gclip("on", x0, yl , x0+xlen, yh);
			break;
		case '/':
			/* reduce trace amplitude */
			pmap[k].uymul /= 2.0;
			*lpos = 0;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			clearregion(0.0,1.0,10.0,8.0);
			gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,ntrc, 
				twin, numperframe, pkabsolute,&yh,&yl,NO, 
				YLIM_OFF, -1, +1, NO);
 			if(pkphasemenu == PK_QUALITY || pkphasemenu == PK_PQUALITY ){
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			}
			gclip("on", x0, yl , x0+xlen, yh);
			break;
		case 'B': /* display previous page  */
			  return PREVIOUS;
			  break;
		case 'D': break; /* set phase Down */
		case 'U': break; /* set phase up */
		case 'E': break; /* set phase onset emergent */
		case 'I': break; /* set phase onset impulsive */
		case 'L':        /* give cursor location and amplitude
		 			note use a small unclip window */
			     uv = gsac_control.uxcen + (ux - 0.5)/gsac_control.uxmul;
			     if(pkabsolute == YES) {
    				ttime = uv*twin  + gsac_control.begminx;
			     } else {
    				ttime = uv*twin + sacdata[pmap[k].k].tzbegx;
			     }
			     if(sacdata[pmap[k].k].sachdr.ihdr[15] == 1){
                		     isxtime = YES ;
			     	     printtimestr(ttime, ostr1);
			     } else {
			     	     isxtime = NO;
				    sprintf(ostr1,"Freq %f Hz",ttime);
			     }
			     /* get the sample  reference to the B time */
			  kk = (int)(0.5 +(float)
				(ttime - sacdata[pmap[k].k].tzref 
				- sacdata[pmap[k].k].sachdr.rhdr[H_B])/
				pmap[k].delta);
			     amp = sacdata[pmap[k].k].sac_data[kk];
			     sprintf(ostr2,"  %12.4g",amp);
			     strcat(ostr1,ostr2);
				/* change 17 JAN 2008 - output Epoch time */
			     sprintf(ostr2,"  %20.4g",ttime);
			     strcat(ostr1,ostr2);
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			  gleft(x0,7.8-(*lpos)*0.15,0.10,ostr1,0.0);
 			if(pkphasemenu == PK_QUALITY || pkphasemenu == PK_PQUALITY ){
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			}
			gclip("on", x0, yl , x0+xlen, yh);
			  *lpos = *lpos + 1;
			  break;
		case 'N':  /* next plot of PERPLOT is set  */
			  return NEXT;
			  break;
		case 'O': /* pervious plot window  */
			  /* this requires a history of windows
			   * and a replot which requires turning of
			   * clipping while a redraw is done */
			*lpos = 0;
			  gsac_control.uxmul = 1.0;
			  gsac_control.uxcen = 0.5;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			clearregion(0.0,1.0,10.0,8.0);
			gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,ntrc, 
				twin, numperframe, pkabsolute,&yh,&yl,YES, 
				YLIM_OFF, -1, +1, NO);
 			if(pkphasemenu == PK_QUALITY || pkphasemenu == PK_PQUALITY ){
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			}
			gclip("on", x0, yl , x0+xlen, yh);
			break;
		case 'P': 
		case 'S': 
		case 'F': 
			/* P time   or S time
                           F is end of useful signal time
			 * get coordinate, get time - problem here
			 * is the ppk relative vs absolute since
			 * we could have the ts and te be trace 
			 * trace dependent - also worry about inside
			 * bounds 
			 * Give the absolute arrival time so that we 
			 * do not have to synchronize */
			     uv = gsac_control.uxcen + (ux - 0.5)/gsac_control.uxmul;
			     if(pkabsolute == YES) {

			     tp = uv*twin + gsac_control.begminx;
			     } else {
			     tp = uv*twin + sacdata[pmap[k].k].tzbegx;
			     }
			     updatehdr(k,switchchar,tp,ns,ne,numperframe);
			*lpos = 0;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			clearregion(0.0,1.0,10.0,8.0);
			gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,ntrc, 
				twin, numperframe, pkabsolute,&yh,&yl,NO, 
				YLIM_OFF, -1, +1, NO);
 			if(pkphasemenu == PK_QUALITY || pkphasemenu == PK_PQUALITY ){
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			}
			gclip("on", x0, yl , x0+xlen, yh);
			     /* change the header and redisplay */
				

			break;
		case 'Q': /* end interactive  */
			return QUIT;
		case 'T': 
			/* user time - then enter n value */
			curaxy(&bx,&by,ch);
			switchchar = ch[0];
			if(isdigit(ch[0])){
			     uv = gsac_control.uxcen + (ux - 0.5)/gsac_control.uxmul;
			     if(pkabsolute == YES) {

			     tp = uv*twin + gsac_control.begminx;
			     } else {
		    	     tp = uv*twin + sacdata[pmap[k].k].tzbegx;
			     }
			     updatehdr(k,switchchar,tp,ns,ne,numperframe);
			*lpos = 0;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			clearregion(0.0,1.0,10.0,8.0);
			gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,ntrc, 
				twin, numperframe, pkabsolute,&yh,&yl,NO, 
				YLIM_OFF, -1, +1, NO);
 			if(pkphasemenu == PK_QUALITY || pkphasemenu == PK_PQUALITY ){
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			}
			gclip("on", x0, yl , x0+xlen, yh);
			}
			break;
		case 'X': 
			cx_X = gsac_control.uxcen + (ux - 0.5)/gsac_control.uxmul;
			/* now invoke the cursor again and only use
			      * of we see a second X */
			curaxy(&cx,&cy,ch);
			switchchar = ch[0];
			switchchar = toupper(switchchar);
			if(inside(cx,cy,x0,y0,x0+xlen,y0+ylen)){
				traceinfo(cx,cy,x0,y0,dy,xlen,ylen,ns,ne,
					ntrc,&k,&ux,&uy ,numperframe);
			if(switchchar == 'X'){
				uv = gsac_control.uxcen + (ux - 0.5)/gsac_control.uxmul;
				/* re map */
				gsac_control.uxcen = (uv + cx_X)/2.0;
				/* redisplay */
				twin_X=ABS(uv - cx_X)*twin;
				gsac_control.uxmul = ABS(twin/twin_X);
				*lpos = 0;
				gclip("off", x0, y0, x0+xlen, y0+ylen);
				clearregion(0.0,1.0,10.0,8.0);
				gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,
					dy,ntrc, twin, numperframe, pkabsolute,
					&yh,&yl,NO, YLIM_OFF, -1, +1, NO);
 			if(pkphasemenu == PK_QUALITY || pkphasemenu == PK_PQUALITY ){
				show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			}
				gclip("on", x0, yl , x0+xlen, yh);
			}
			}
			  break; /* invoke twice to define new window */
		case 'n': 
			break; /* quality value n = 0 1 2 3 4 */
		/*  process quality menu */
		case 'A': 
			if(pkphasemenu == PK_QUALITY ||
				pkphasemenu == PK_PQUALITY ){
				upd_hdr_ihdr( k, YES ,  H_IHDR20 ) ;
				do_check(YES ,k );
			}
			break; 
		case 'R': 
			if(pkphasemenu == PK_QUALITY ||
				pkphasemenu == PK_PQUALITY ){
				upd_hdr_ihdr( k, NO ,  H_IHDR20 ) ;
				do_check(NO ,k );
			}
			break; 
		case 1:	/* left mouse button */
		case 2:	/* right mouse button */
		case 3:	/* center mouse button */
			if(pkphasemenu == PK_QUALITY){
				upd_hdr_ihdr( k, (pkquality == MENU_QUALITY_ACCEPT) ,  H_IHDR20 ) ;
				do_check((pkquality == MENU_QUALITY_ACCEPT) ,k );
			} else if(pkphasemenu == PK_PQUALITY){
			/* P time   or S time'
			* get coordinate, get time - problem here
			 * is the ppk relative vs absolute since
			 * we could have the ts and te be trace 
			 * trace dependent - also worry about inside
			 * bounds 
			 * Give the absolute arrival time so that we 
			 * do not have to synchronize */
			     uv = gsac_control.uxcen + (ux - 0.5)/gsac_control.uxmul;
			     if(pkabsolute == YES) {

			     tp = uv*twin + gsac_control.begminx;
			     } else {
			     tp = uv*twin + sacdata[pmap[k].k].tzbegx;
			     }
			     updatehdr(k,'P',tp,ns,ne,numperframe);
			*lpos = 0;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			clearregion(0.0,1.0,10.0,8.0);
			gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,ntrc, 
				twin, numperframe, pkabsolute,&yh,&yl,NO, 
				YLIM_OFF, -1, +1, NO);
			     /* change the header and redisplay */
				upd_hdr_ihdr( k, (pkquality == MENU_QUALITY_ACCEPT) ,  H_IHDR20 ) ;
				do_check((pkquality == MENU_QUALITY_ACCEPT) ,k );
		show_check( x0, y0, xlen, ylen, ns, ne, dy, ntrc, numperframe, yh, yl);
			gclip("on", x0, yl , x0+xlen, yh);
			}
			break;

	}
	return (-1);
}

void traceinfo(float cx,float cy,float x0,float y0, float dy, float xlen,float ylen,int ns, int ne, int ntrc, int *k, float *ux, float *uy , int numperframe)
{
	float uxv, uyv;
	int kv;
	int kk;
	/* first determine the location of the */
	uxv = (cx - x0)/xlen;
	uyv = (cy - y0)/ylen;
	kv = (ne - ns) - ((uyv)*numperframe);
	*ux = uxv; 
	*uy = uxv; 
	*k = kv;
	*k = -1;
	for(kk = ns%numperframe ; kk <= (ne-1)%numperframe ; kk++){
		if(pmap[kk].xl <= cx && cx <= pmap[kk].xh
				&& pmap[kk].yl <= cy && cy <= pmap[kk].yh){
			*k = kk;
			return;
		}
	}
}

void upd_hdr_time(int k, float tph, int H, int HK, char *str){
	sacdata[pmap[k].k].sachdr.rhdr[H] = tph;
	strncpy(sacdata[pmap[k].k].sachdr.chdr[HK],str,8);
	strcpy (sacdata[pmap[k].k].schdr[HK]      ,str  );
}

void upd_hdr_ihdr(int k, int val,  int H )
{
	sacdata[pmap[k].k].sachdr.ihdr[H] = val;
}

void updatehdr(int k,int switchchar,double tp,int ns,int ne,int numperframe)
{
	int i;
	float tph;	/* offset with respect to trace reference time */
	tph = (float)(tp - sacdata[pmap[k].k].tzref);
/*
printf("swithchar %c H_F %d H_KF %d tph %f\n",switchchar,H_F,H_KF,tph);
*/
	switch(switchchar){
		case 'P': upd_hdr_time(k,tph,H_A ,H_KA ,"iP      "); break;
		case 'S': upd_hdr_time(k,tph,H_T0,H_KT0,"iS      "); break;
		case '0': upd_hdr_time(k,tph,H_T0,H_KT0,"T0      "); break;
		case '1': upd_hdr_time(k,tph,H_T1,H_KT1,"T1      "); break;
		case '2': upd_hdr_time(k,tph,H_T2,H_KT2,"T2      "); break;
		case '3': upd_hdr_time(k,tph,H_T3,H_KT3,"T3      "); break;
		case '4': upd_hdr_time(k,tph,H_T4,H_KT4,"T4      "); break;
		case '5': upd_hdr_time(k,tph,H_T5,H_KT5,"T5      "); break;
		case '6': upd_hdr_time(k,tph,H_T6,H_KT6,"T6      "); break;
		case '7': upd_hdr_time(k,tph,H_T7,H_KT7,"T7      "); break;
		case '8': upd_hdr_time(k,tph,H_T8,H_KT8,"T8      "); break;
		case '9': upd_hdr_time(k,tph,H_T9,H_KT9,"T9      "); break;
		case 'F': upd_hdr_time(k,tph,H_F ,H_KF ,"F       "); break;
	}
	if(pkmarkall){
		for(i = ns%numperframe ; i <= (ne-1)%numperframe ; i++){
			tph = (float)(tp - sacdata[pmap[i].k].tzref);
		switch(switchchar){
			case 'P': upd_hdr_time(i,tph,H_A ,H_KA ,"P       "); break;
			case 'S': upd_hdr_time(i,tph,H_T0,H_KT0,"S       "); break;
			case '0': upd_hdr_time(i,tph,H_T0,H_KT0,"T0      "); break;
			case '1': upd_hdr_time(i,tph,H_T1,H_KT1,"T1      "); break;
			case '2': upd_hdr_time(i,tph,H_T2,H_KT2,"T2      "); break;
			case '3': upd_hdr_time(i,tph,H_T3,H_KT3,"T3      "); break;
			case '4': upd_hdr_time(i,tph,H_T4,H_KT4,"T4      "); break;
			case '5': upd_hdr_time(i,tph,H_T5,H_KT5,"T5      "); break;
			case '6': upd_hdr_time(i,tph,H_T6,H_KT6,"T6      "); break;
			case '7': upd_hdr_time(i,tph,H_T7,H_KT7,"T7      "); break;
			case '8': upd_hdr_time(i,tph,H_T8,H_KT8,"T8      "); break;
			case '9': upd_hdr_time(i,tph,H_T9,H_KT9,"T9      "); break;
		}
		}
	}
}

void updthd(int k,int tn, int ktn, char *kstr ,double tp,int ns,int ne,int numperframe)
{
	int i;
	float tph;	/* offset with respect to trace reference time */
	/* clean up the string by replacing nulls by blanks */
	for(i = 0 ; i < 8 ; i++)
		if(kstr[i] == '\0')
			kstr[i] = ' ' ;
	if(pickundo == YES){
		upd_hdr_time(k,-12345.0,tn ,ktn , kstr);
	} else {
		tph = (float)(tp - sacdata[pmap[k].k].tzref);
		upd_hdr_time(k,tph,tn ,ktn , kstr);
	}
	if(pkmarkall){
		for(i = ns%numperframe ; i <= (ne-1)%numperframe ; i++){
			if(pickundo == YES){
				upd_hdr_time(i,-12345.0,tn ,ktn , kstr);
			} else {
				tph = (float)(tp - sacdata[pmap[i].k].tzref);
				upd_hdr_time(i,tph,tn ,ktn , kstr);
			}
		}
	}
}



/* phase menu codes */

/* initially put up the menu horizontally 
* if a phase is picked then put up quality
* if P is picked then put up polarity
* */


void show_polmenu(void)
{
	show_menu(1.0, 0.3, ppol, sizeof(ppol),&npl);
}

void show_qualmenu(void)
{
	show_menu(1.0, 0.3, pqual, sizeof(pqual),&nqa);
}

/* process phase menu 
 * if arrow is in the region do something, else nothing - out in a loop
 * if any phase force trave pick THEN
 * put up quality of phase - then erase
 * then if P put up polarity then erase
 * Update the phase header values, the string and also redraw at the end
 * */
void proc_phmenu(float cx,float cy,float x0,float y0,float dy,
	float xlen,float ylen,int ns,int ne, float yh, float yl,
	int ntrc, float ts, float te, double twin, int numperframe,
	struct menu *p)
{
float lx, ly, hx, hy;
int i, j, cmd;
char ch[2];
float vx,vy;
int k;
float ux, uy, uv;
double tp;
int doloop;
char ostr[9];
	/* initialize */
	for(i=0 ; i < nmp ; i++){
		lx = p[i].xl;
		ly = p[i].yl;
		hx = p[i].xh;
		hy = p[i].yh;
		cmd = p[i].action;
                if(cmd ==  MENU_P_UNDO ){
			if(last_h > 0){
				pickundo = YES;
				updthd(last_k,last_h ,last_kh ,"-12345   ",tp,ns,ne,numperframe);
				last_h = -1;
				pickundo = NO;
				gmesg("                     ");
				gclip("off", x0, y0, x0+xlen, y0+ylen);
				clearregion(0.0,1.0,10.0,8.0);
				gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,
					ntrc, twin, numperframe, pkabsolute,&yh,&yl,
					NO, YLIM_OFF, -1, +1, NO);
				gclip("on", x0, yl , x0+xlen, yh);
			}
			return;
		}
		if(inside(cx,cy,lx,ly,hx,hy)){
			/* force trace pick */
			gmesg("Pick the phase");
			doloop = YES;
			while(doloop){
			/* loop until we get a trace pick */
			gclip("on", x0, yl , x0+xlen, yh);
			curaxy(&vx,&vy,ch);
			if(inside(vx,vy,x0,y0,x0+xlen,y0+ylen)){
				doloop = NO;
				traceinfo(vx,vy,x0,y0,dy,xlen,ylen,ns,ne,
					ntrc,&k,&ux,&uy , numperframe);
			     	uv = gsac_control.uxcen + (ux - 0.5)/gsac_control.uxmul;
			     	if(pkabsolute == YES) {

			     	tp = uv*twin + gsac_control.begminx;
			     	} else {
			     	tp = uv*twin + sacdata[pmap[k].k].tzbegx;
			     	}
				for(j=0;j < 8; j++)
					ostr[j] = ' ' ;
				ostr[8] = '\0';
					/* get the quality */
					proc_qualmenu(ostr,x0,y0,xlen,ylen);
					/* put in the phase name */
					strcat(ostr,p[i].str);
					if(cmd == MENU_PREG_P ){
						proc_polmenu(ostr,x0,y0,xlen,ylen);
					} else  if(cmd == MENU_PTEL_P ){
						proc_polmenu(ostr,x0,y0,xlen,ylen);
					}
				switch(cmd){
				/* regional phases */
				case MENU_PREG_P  :
					last_h  = H_A;
					last_kh = H_KA;
					last_k = k;
					updthd(k,H_A ,H_KA ,ostr,tp,ns,ne,numperframe);
					break;
				case MENU_PREG_S  :
					last_h  = H_T0;
					last_kh = H_KT0;
					last_k = k;
					updthd(k,H_T0,H_KT0,ostr,tp,ns,ne,numperframe);
					break;
				case MENU_PREG_Lg :
					last_h  = H_T2;
					last_kh = H_KT2;
					last_k = k;
					updthd(k,H_T2,H_KT2,ostr,tp,ns,ne,numperframe);
					break;
				case MENU_PREG_Pg :
					last_h  = H_T3;
					last_kh = H_KT3;
					last_k = k;
					updthd(k,H_T3,H_KT3,ostr,tp,ns,ne,numperframe);
					break;
				/* teleseismic phases */
				case MENU_PTEL_P :
					last_h  = H_A;
					last_kh = H_KA;
					last_k = k;
					updthd(k,H_A ,H_KA ,ostr,tp,ns,ne,numperframe);
					break;
				case MENU_PTEL_S :
					last_h  = H_T0;
					last_kh = H_KT0;
					last_k = k;
					updthd(k,H_T0,H_KT0,ostr,tp,ns,ne,numperframe);
					break;
				case MENU_PTEL_PKP :
					last_h  = H_T1;
					last_kh = H_KT1;
					last_k = k;
					updthd(k,H_T1,H_KT1,ostr,tp,ns,ne,numperframe);
					break;
				case MENU_PTEL_PcP :
					last_h  = H_T2;
					last_kh = H_KT2;
					last_k = k;
					updthd(k,H_T2,H_KT2,ostr,tp,ns,ne,numperframe);
					break;
				case MENU_PTEL_PP :
					last_h  = H_T3;
					last_kh = H_KT3;
					last_k = k;
					updthd(k,H_T3,H_KT3,ostr,tp,ns,ne,numperframe);
					break;
				case MENU_PTEL_pP :
					last_h  = H_T4;
					last_kh = H_KT4;
					last_k = k;
					updthd(k,H_T4,H_KT4,ostr,tp,ns,ne,numperframe);
					break;
				case MENU_QUALITY_ACCEPT :
					pkquality = MENU_QUALITY_ACCEPT ;
fprintf(stderr,"MENU_QUALITY_ACCEPT\n");
					gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
						"Mouse Click to Accept (also a - Accept, r - Reject)",0.0);
					break;
				case MENU_QUALITY_REJECT :
fprintf(stderr,"MENU_QUALITY_REJECT\n");
					pkquality = MENU_QUALITY_REJECT ;
					gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
						"Mouse Click to Reject (also a - Accept, r - Reject)",0.0);
					break;
			}
			gmesg("                     ");
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			clearregion(0.0,1.0,10.0,8.0);
			gsac_show_plotpk(x0,y0,xlen,ylen,ns,ne,ts,te,dy,
				ntrc, twin, numperframe, pkabsolute,&yh,&yl,
				NO, YLIM_OFF, -1, +1, NO);
			gclip("on", x0, yl , x0+xlen, yh);
			}
			}
			return;
		}
	}
}

/* process quality menu 
 * if arrow is in the region do something, else nothing - out in a loop
 * Update the  header values, the string and also redraw at the end
 * */
void proc_qlmenu(float cx,float cy,float x0,float y0,float dy,
	float xlen,float ylen,int ns,int ne, float yh, float yl,
	int ntrc, float ts, float te, double twin, int numperframe,
	struct menu *p)
{
float lx, ly, hx, hy;
int i, cmd;
	for(i=0 ; i < nmp ; i++){
		lx = p[i].xl;
		ly = p[i].yl;
		hx = p[i].xh;
		hy = p[i].yh;
		cmd = p[i].action;
		if(inside(cx,cy,lx,ly,hx,hy)){
				switch(cmd){
				case MENU_QUALITY_ACCEPT :
					pkquality = MENU_QUALITY_ACCEPT ;
					gclip("off", x0, y0, x0+xlen, y0+ylen);
					clearregion(x0,y0+ylen+0.10,x0+xlen,y0+ylen+0.30);
					newpen(1);
					gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
						"Mouse Click to Accept (also a - Accept, r - Reject)",0.0);
					gclip("on", x0, yl , x0+xlen, yh);
					break;
				case MENU_QUALITY_REJECT :
					pkquality = MENU_QUALITY_REJECT ;
					gclip("off", x0, y0, x0+xlen, y0+ylen);
					clearregion(x0,y0+ylen+0.10,x0+xlen,y0+ylen+0.30);
					newpen(1);
					gcent(x0+xlen/2.,y0+ylen+2.*0.10,0.10,
						"Mouse Click to Reject (also a - Accept, r - Reject)",0.0);
					gclip("on", x0, yl , x0+xlen, yh);
					break;
			}
		}
	}
	return;
}


void proc_polmenu(char *ostr, float x0, float y0, float xlen, float ylen)
{
float xl, yl, xh, yh;
int i, cmd;
char ch[2];
float cx, cy;
int doloop = YES;
	show_polmenu();
	gclip("on", x0, y0, x0+xlen, y0+ylen);
	gmesg("Select P polarity");
	while(doloop){
		curaxy(&cx,&cy,ch);
		for(i=0 ; i < npl ; i++){
			xl = ppol[i].xl;
			yl = ppol[i].yl;
			xh = ppol[i].xh;
			yh = ppol[i].yh;
			cmd = ppol[i].action;
			if(inside(cx,cy,xl,yl,xh,yh)){
				strcat(ostr,"_");
				strcat(ostr,ppol[i].str);
				doloop = NO;
			}
		}
	}
	gclip("off", x0, y0, x0+xlen, y0+ylen);
	clearregion(0.0,0.0,10.0,0.6);
	gclip("on", x0, y0, x0+xlen, y0+ylen);
}

void proc_qualmenu(char *ostr,float x0, float y0, float xlen, float ylen)
{
float xl, yl, xh, yh;
int i, cmd;
char ch[2];
float cx, cy;

int doloop = YES;
	show_qualmenu();
	gclip("on", x0, y0, x0+xlen, y0+ylen);
	gmesg("Select phase quality");
	while(doloop){
		curaxy(&cx,&cy,ch);
		for(i=0 ; i < nqa ; i++){
			xl = pqual[i].xl;
			yl = pqual[i].yl;
			xh = pqual[i].xh;
			yh = pqual[i].yh;
			cmd = pqual[i].action;
			if(inside(cx,cy,xl,yl,xh,yh)){
				strcpy(ostr,pqual[i].str);
				doloop = NO;
			}
		}
	}
	gclip("off", x0, y0, x0+xlen, y0+ylen);
	clearregion(0.0,0.0,10.0,0.6);
	gclip("on", x0, y0, x0+xlen, y0+ylen);
}

void do_check(int val, int k)
{
	float xl, xh, yl, yh;
	float dx, dy, ht;
	xl = pmap[k].xl ;
	xh = pmap[k].xh ;
	yl = pmap[k].yl ;
	yh = pmap[k].yh ;
	/* put in the X or the check here */
	dx = ABS(xh-xl);
	dy = ABS(yh-yl);
	ht = 0.5*dy ;
	clearregion(xl+0.25*ht,yl+0.25*ht,xl+0.75*ht,yl+0.75*ht);
	if(val == YES){
		newpen(2);
		gcent(xl+0.5*ht, yl+0.25*ht,0.5*ht,"+",0.0);
		newpen(1);
	} else if (val == NO){
		newpen(4);
		gcent(xl+0.5*ht, yl+0.25*ht,0.25*ht,"X",0.0);
		newpen(1);
	}
}	

void show_check(float x0, float y0, float xlen, float ylen, int ns, int ne, float dy, int ntrc, int numperframe, float yh, float yl)
{
	int kkk, k;
	for(kkk = ns ; kkk < ne ; kkk++){
		k = sortptr[kkk];
		do_check(sacdata[k].sachdr.ihdr[H_IHDR20], kkk%numperframe);
	}

}
