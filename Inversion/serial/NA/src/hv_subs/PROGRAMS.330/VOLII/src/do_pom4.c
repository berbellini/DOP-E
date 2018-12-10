/* to do
	for MATCH
	Add, Delete
		Add is a fictitious that must be removed on return
		Delete can delete a previous pick, but also a fictittious
	Mode should highlight the current selection - e.g., XOR
	The zoom box should always plot the proper period at the bottom,
		or should just zoom on the period only not on velocity
	If zooom stat set when tomath reset

	for zoom and interp
		use zoom to define the period limit
	go through the list and select, sort observations then interpolate
	use pointer? since many more -- now for match only have
	one mode to work with so can easily define simpel output array



  03 16 00 - mode select implemented - now for do pick we should
	use the mode value and not the original amplitude value
	also select on from already picked instead of from 
        nearest neighbor - also note in the matchprep we use
        instantaneous period not actual period

  10 10 00 - new output STAcomp.ods

  12 27 00 - minor changes to colors
  05 MAR 2002 - use the flag 'X' for phase velocity output in SURF96 format
  02 NOV 2003 - cleaned up WIN32/MSDOS interface
*/

#include	<stdio.h>
#include	<ctype.h>
#include	<stdlib.h>
#include	<string.h>
#include	"calplot.h"
#include	<math.h>
#include	"grphsubc.h"
#include	"nmenu.h"



FILE *cerr;


/* DEFINES */
#define ON      1
#define OFF     0
#define YES     1
#define NO      0
/* plot window for POM96 graphs 
	   The plot will be placed between (XL,YL) and (XH,YH) 
	   The POM96.PLT will be shifted by adding (XOFF,YOFF) */
#define XL      0.0
#define XH      10.00
#define YL      1.5
#define YH      8.0
#define XOFF	0.0
#define YOFF	-0.5;
#define MIN(a,b) ( (b) > (a) ? (a):(b) )
#define MAX(a,b) ( (a) > (b) ? (a):(b) )
#define ABS(a)   ( (a) > (0) ? (a): -(a) )
#define LIN 0
#define LOG 1

/* STRUCTURE DEFINITIONS */

struct disp {
	char 	type[7]	;	/* identification of type 	*/
	char	lvry[2]	;	/* A undef L Love R Rayleigh	*/
	char	cug[2]	;	/* U = group, C= phase, G= gamma */
	int	mode	;	/* -1 not defined, 0 = Fund, 1 = first */
	float	per	;	/* period seconds */
	float	vel	;	/* velocity value	*/
	float	evel	;	/* error velocity value	*/
	float	amp	;	/* spectral amplitude */
	float	dist	;	/* distance */
	float	lat1	;	/* event latitude */
	float	lon1	;	/* event longitude */
	float	lat2	;	/* station latitude */
	float	lon2	;	/* station longitude */
	int	ictl	;	/* control flag */
	int	symb	;	/* symbol used in plot */
	float	instper	;	/* instantaneous period */
	char	comment[10];	/* comment field */
	char	kstnm[9];	/* station name */
	char	kcmpnm[9];	/* componet name */
	int	nzyear	;	/* year */
	int	nzjday	;	/* julian day */
	int	nzhour	;	/* hour */
	int	nzmin	;	/* minute */
};

struct pos {
	float xl;
	float xh;
	float yl;
	float yh;
	float axl;
	float axh;
	float ayl;
	float ayh;
	char xlnlg[4];
	char ylnlg[4];
	char fname[20];
} pompos[3];


#define	ZOOM	1
#define	UNZOOM	2
#define	AUTO	3
#define	MODE	4
#define	PICK	5
#define	RESTART	6
#define	MATCH	7
#define	EXIT	8
#define INTERP	9
#define ADD	10
#define RETURN	11
#define DELETE	12

static struct menu m[] = {
		{  -1.0, -1.0, -1.0, -1.0, "Zoom\0"   , ZOOM   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "UnZoom\0" , UNZOOM ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Mode\0"   , MODE   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Auto\0"   , AUTO   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Pick\0"   , PICK   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Restart\0", RESTART,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Exit\0"   , EXIT   ,-1, 1, 1},
};

static struct menu mm[] = {
		{  -1.0, -1.0, -1.0, -1.0, "Zoom\0"   , ZOOM   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "UnZoom\0" , UNZOOM ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Mode\0"   , MODE   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Interp\0" , INTERP ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Add\0"    , ADD    ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Delete\0" , DELETE ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Return\0" , RETURN, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Exit\0"   , EXIT   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Match\0"  , MATCH  ,-1, 1, 1},
};


static struct menu md[] = {
		{  -1.0, -1.0, -1.0, -1.0, "Fund\0" , 1, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "1 st\0" , 2, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "2 nd\0" , 3, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "3 rd\0" , 4, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "4 th\0" , 5, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "5 th\0" , 6, -1, 1, 1},
};

static struct menu my[] = {
		{  -1.0, -1.0, -1.0, -1.0, "Yes\0" , 1, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "No \0" , 2, -1, 1, 1},
};


/* GLOBAL VARIABLE DECLARATIONS */
	/* display information */
extern int HasMouse; 
extern float XminDev, YminDev, 
	XmaxDev, YmaxDev, XminClip, 
	YminClip, XmaxClip, YmaxClip;
extern int Color;
extern int doshade;

extern int Units;
extern float border;
extern float title;
extern int black, kolor;
extern int Button_Color;
extern int Button_Color_Light;
extern int Button_Color_Dark;
extern int Button_Color_Fore;
extern int Button_Color_Back;

int Mode  = -1;
int iMode = -1;	/* pointer into mode array */
int InPick = OFF;
int InLine = OFF;
static int XAxisPeriod = NO;
static int ModeSelect  = -1;	/* for use in interpolation selection */
char ostr[512];
struct disp *pomdsp;
int num_disp;
/* arrays for amplitude and dispersion plots */
/* we will place 3 graphs in the page :
	[0]  spectra, [1] group velocity, [2] trace
	The positions of the original plot file will be mapped
	on to the following boxes:
*/
/* screenplacement */
float  	pxl[3]= { 0.5, 2.1, 9.4},
	pxh[3]= { 1.5, 7.9, 9.9},
	pyl[3]= { 2.0, 1.5, 2.0},
	pyh[3]= { 7.0, 7.3, 7.0};
float apxl[3], apyl[3], apxh[3], apyh[3];	/* user values of corners */

float opxl[3], opyl[3], opxh[3], opyh[3];	/* original cut corners */
float sclx[3], scly[3], scl[3]		;	/* scale factors */

int xlnlg[3], ylnlg[3];

float *period;
int num_period;
int Pick_Color;


/* FUNCTION DECLARATIONS */
extern void gsolid(float x, float y, float h, int symb);
extern char *my_pathfind( char *path, const char *name, const char *mode);
void do_gpv_pick(void);
void get_pom_output(void );
void draw_button(float xl, float yl, char *str, float *xlw, float *ylw,
float *xup, float *yup, int butrev, int lstrmax) ;
void do_scale(float per,float vel,float *xval,float *yval,int i);
void do_invscale(float *per,float *vel,float xval,float yval,int i);
int get_action(int nm, char *fname);
int get_action1(int nm, char *fname);
void do_pick( float xv, float yv, int lev);
void do_line( float xv, float yv);
void do_zoom( void);
void do_zoomm( void);
void do_clear( int fig);
void do_unzoom( void);
void do_unzoomm( void);
void do_resetzoom(void);
int do_mode( void);
void do_showpick1(void);
void do_showmode1(int mode);
void do_showpk();
void do_showmd(int mode);
void do_showperscl(int i);
void do_restart(int nm);
void do_box(float xl,float yl,float xh,float yh);
void show_menu (float x0, float y0, struct menu *men, int size, int *nm);
void redisplay_pick(void);
void do_start(int *nm, int menuon);
	float pinterp(float x,float x1,float x2);
int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);
int do_matchprep(char *fname);
void doout(int *ret);
void doout1(int mode);

void do_axes(int i);
void do_clearscale(int i);
void mgwrtxt(float x0, float y0, char *str, int cmd, int color);

#define		RED	2
#define		WHITE	0
#define		BLACK	1



char *wavetype;
int wave;


int do_page4(char *type, int dotype)
{
	int nm, ret;
	if(doshade)
		Pick_Color = WHITE;
	else
		Pick_Color = RED;
	do_start(&nm, ON);
	wavetype = type;
	wave = dotype;
	ret = get_action(nm," ");
	newpen(1);
	gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
	gframe(1);
	free(pomdsp);
	return (ret);
}

int get_action(int nm ,char *fname)
{
	int cmd, i, ii;
	int ret;
	int j;
	float xv, yv;
	float xv1, yv1;
	float xl, yl, xh, yh;
	float per, vel, xval, yval;
	char c[2];
	/* annotate top */
/*
	mgwrtxt(2.0,7.4,"None\0",2, 1);
*/
	do_showmd(Mode);
	do_showpk();
	mgwrtxt(0.5,7.7,wavetype,1,2);
	
	/* prohibit writing at the bottom */
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	gcursor("Arrow");
	for(; ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
		/* loop on commands */
		for(i=0 ; i < nm ; i++){
			if(inside(xv,yv,m[i].xl,m[i].yl,m[i].xh,m[i].yh))
			{
				cmd = m[i].action;
				gmesg(m[i].str);
				ii = i;
				break;
			}
		}
		if(cmd > 0){
			switch(cmd){
			case ZOOM:	/* Zoom */
				draw_button(0.5+ii*1.125,0.6 ,
				    m[ii].str,&xl,&yl,&xh,&yh,ON ,m[ii].lstrmx);
				gcursor("Arrow");
				do_zoom();
				gcursor("Arrow");
				draw_button(0.5+ii*1.125,0.6 ,
				    m[ii].str,&xl,&yl,&xh,&yh,OFF,m[ii].lstrmx);
				break;
			case UNZOOM:	/* UnZoom */
				InPick = OFF;
				InLine = OFF;
				do_showpk();
				do_unzoom();
				break;
			case AUTO:	/* Automatic Pick from nearest line */
				if(Mode <0){
					gmesg("Choose a Mode First");
					Mode = do_mode();
					do_showmd(Mode);
				}
				InPick = OFF;
				InLine = ON;
				do_showpk();
				curaxy(&xv1, &yv1, c);
				do_line(xv1, yv1);
				break;
			case MODE:	/* Mode */
				Mode = do_mode();
				do_showmd(Mode);
				break;
			case PICK:	/* Pick */
				if(Mode <0){
					gmesg("Choose a Mode First");
					Mode = do_mode();
					do_showmd(Mode);
				}
				if(InPick == ON) {
					InPick = OFF ;
				} else {
					InPick = ON;
				}
				InLine = OFF;
				do_showpk();
				break;
			case RESTART:	/* Restart */
				do_restart(nm);
				InPick = OFF;
				InLine = OFF;
				do_showpk();
				break;
			case MATCH:	/* Match */
				gframe(1);
				gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
				/* turn off XOR */
				newpen(3000);
				/* reset the zoom parameters */
				do_resetzoom();
				ret = do_matchprep(fname);
				gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
				if(ret > 0){
					goto end;
				} else {
					gmesg(fname);
					show_menu(0.5, 0.6, m, sizeof(m), &nm);
					InPick = OFF;
					InLine = ON;
					Mode = -1;
					/* redisplay the plots */
					do_resetzoom();
					for(j=1;j<2;j++){
						gread(pompos[j].fname,pxl[j],pyl[j],
	    						opxl[j],opyl[j],opxh[j],opyh[j],
	    						1, sclx[j],scly[j]);
					}
					do_axes(1);
					/* redisplay the current picked 
						values of both 
						amplitude and velocity in  XOR */
					do_showmd(Mode);
					do_showpk();
					mgwrtxt(0.5,7.7,wavetype,1,2);
					/* turn on XOR */
					newpen(3001);
					for(j=0;j<num_disp;j++){
						if(pomdsp[j].mode >= 0){
							per = pomdsp[j].per;
							vel = pomdsp[j].vel;
							/* display velocity value */
							do_scale(per,vel,&xval,&yval,1);
	newpen(Pick_Color);
							gsolid(xval,yval,0.02*scl[1],pomdsp[j].symb);
							/* display spectral amplitude value */
						}
					}
				/* turn off XOR */
				newpen(3000);
				gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
				}
				break;
			case EXIT:
				ret = 1;
				goto end;
			default:
				break;

			}
		} else {
			if(InPick == ON){
				if(inside(xv,yv,pxl[1]-0.1,pyl[1],pxh[1]+0.1,pyh[1]))
				{
					do_pick(xv, yv, 0);
				}
			}
			if(InLine == ON){
				if(inside(xv,yv,pxl[1]-0.1,pyl[1],pxh[1]+0.1,pyh[1]))
				{
					do_line(xv, yv);
				}
			}
		}
	}
end:
	/* turn off XOR */
	newpen(3000);
	gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
	gmesg("Saving Dispersion Files");
	/* interactively ask if the selected dispersion should be saved */
	doout(&ret);
	/* free period array */
	free(period);
	return(ret);
} 

void do_pick( float xv1, float yv1, int lev)
{
	/* from the screen coordinates, (xv1, yv1) 
		obtain the nearest point in the dispersion table 
	   lev	= 0 do with original list
		= 1 work with the already picked list and inst period
	*/
	float xval, yval;
	float fxval;
	int i;
	float sum_min;
	float sum;
	int ii = -1;
	float per, vel;
	float xs = 0.0, ys = 0.0;
	float ps = 0.0, vs = 0.0;
	/* convert screen coordinates to the (period,velocity) pair */
	sum_min = 1.0e+38; 

	/* scan through dispersion file and line nearest to line */
	/* find the closet x-value first, then the closest y-value */
	for(i=0; i< num_disp ; i++){
		if(lev == 0) {
			per = pomdsp[i].per;
			vel = pomdsp[i].vel;
		} else {
                        if(!XAxisPeriod)
                                per = 1.0/pomdsp[i].instper;
                        else
                                per = pomdsp[i].instper;
			vel = pomdsp[i].vel;
			if( pomdsp[i].mode != ModeSelect)
				per = -1.0;
		}
		if(per >= 0.0){
			do_scale(per,vel,&xval,&yval,1);
			do_scale(pomdsp[i].instper,vel,&fxval,&yval,1);
			sum = (xval-xv1)*(xval -xv1)+(yval-yv1)*(yval-yv1); 
			if(sum < sum_min){
				sum_min = sum ;
				ii = i;
				xs = xval; 
				ys = yval;
				ps = per; 
				vs = vel; 
			}
		}

	}
	if(ii >= 0 && ii < num_disp){
		/* plot a new XOR symbol there */
		newpen(3001);
		newpen(Pick_Color);
		/* put up the velocity  point */
		gsolid(xs,ys,0.02*scl[1],pomdsp[ii].symb);
		gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
		newpen(4);
		/* turn off XOR */
		newpen(3000);
		/* toggle the assignment */
		if(pomdsp[ii].mode > -1)
			pomdsp[ii].mode = -1;
		else
			pomdsp[ii].mode = Mode;
		gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	}

}


void do_line(float xv1, float yv1)
{
	char c[2];
	float xv2, yv2;
	float x1,x2,y1,y2;
	float ypred;
	float per, vel;
	float dataper;
	float xval, yval;
	float xxval;
	int i,j;
	int ii;
	float sum_min_x, sum_min_y;
	float sum_x, sum_y;
	float xs, ys, ps, vs;

	gcursor("Rubber");
	curaxy(&xv2, &yv2, c);
	gcursor("Arrow");
	x1 = MIN(xv1,xv2) ;
	x2 = MAX(xv1,xv2) ;
	y1 = MIN(yv1,yv2) ;
	y2 = MAX(yv1,yv2) ;
	if(x1==x2 && y1 == y2) return;  /* need two points for a line */
	if( inside(x1,y1,pxl[1]-0.1,pyl[1],
	    pxh[1]+0.1, pyh[1])
	    && inside(x2,y2,pxl[1]+0.1,pyl[1],
	    pxh[1]+0.1, pyh[1]) ) {
		newpen(3001);
		/* scan through dispersion file and line nearest to line */
		/* since we are plotting according to the filter periods */
		for(j = 0 ; j < num_period ; j++){
			/* rest goodness of fit */
			sum_min_x = 1.0e+38; 
			sum_min_y = 1.0e+38; 
			if(!XAxisPeriod)
				per = 1.0/period[j];
			else
				per = period[j];
			vel = 1.0;
			do_scale(per,vel,&xxval,&yval,1);
			if(xxval >= x1 && xxval <= x2){
				ii = -1;
				/* this period is within the cursor range */
				/* find the data period closest to this value */
				for(i=0; i< num_disp ; i++){
					sum_x = ABS(per - pomdsp[i].per);
					if(sum_x < sum_min_x){
						sum_min_x = sum_x;
						dataper = pomdsp[i].per;
					}
				}

				/* now make a line prediction */
				/* beware of divide by zero */
				ypred = yv1 + (xxval -xv1)*(yv2-yv1)/(xv2-xv1); 
				/*we now find the point closest to this period*/
				for(i=0; i< num_disp ; i++){
					per = pomdsp[i].per;
					vel = pomdsp[i].vel;
					if(ABS(per - dataper) < 0.001*dataper){
					do_scale(per,vel,&xval,&yval,1);
					if(xval >= x1 && xval <= x2){
						sum_y = ABS(yval-ypred); 
						if(sum_y < sum_min_y){
							sum_min_y = sum_y ;
							ii = i;
							xs = xval; 
							ys = yval;
							ps = per; 
							vs = vel; 
						}
					}
					}
/*
*/
				}
				/* now plot only if we have a hit */
				if(ii >= 0 && ii < num_disp && sum_min_y < 0.10){
					if(pomdsp[ii].mode < 0){
					gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
					newpen(Pick_Color);
					gsolid(xs,ys,0.02*scl[1],pomdsp[ii].symb);
					gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
					gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
					}
					pomdsp[ii].mode = Mode;
				}
			}
		}
		/* turn off XOR */
		newpen(3000);
	}
}

void do_zoom(void)
{
	char c[2];
	float xv1, xv2, yv1, yv2;
	float xfrac, yfrac;
	float x1,x2,y1,y2;
	float tx1,tx2,ty1,ty2;
	float per1, vel1;
	float per2, vel2;
	float px1,px2,py1,py2;
	/* turn off XOR */
	newpen(3000);
	curaxy(&xv1, &yv1, c);
	if( !inside(xv1,yv1, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1]))return;
	gcursor("Box");
	curaxy(&xv2, &yv2, c);
	gcursor("Arrow");
	if(xv1==xv2 && yv1 == yv2) return;  /* need two points for a line */
	/* 
	We work with three coordinate systems, each consisting of two
	opposite corners of a view rectangle

	(pxl,pyl) -> (pxh,pyl)   is display screen which does not change
	(opxl,opyl) -> (opxh,opyh) is the mapped portion of the original plot
	(apxl,apyl) -> (apxh,apyl) are user coordinates of the original space

	To do the mapping, use a simple interpolation

	V(p) = (1-p)V1 + pV2   where  0 <= p <= 1
	*/
	if( inside(xv1,yv1, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1])
	    && inside(xv2,yv2, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1]) ) {
		/* display all current picks to turn them off */
		/* since the point may extend across clip region ? */
		redisplay_pick();
		do_clear(1);
		/* these are screen coordinates */
		x1 = MIN(xv1,xv2) ;
		x2 = MAX(xv1,xv2) ;
		y1 = MIN(yv1,yv2) ;
		y2 = MAX(yv1,yv2) ;
		px1 = pinterp(x1,pxl[1],pxh[1]);
		px2 = pinterp(x2,pxl[1],pxh[1]);
		py1 = pinterp(y1,pyl[1],pyh[1]);
		py2 = pinterp(y2,pyl[1],pyh[1]);
		/* now estimate the coordinates in the
					original plot space */
		tx1 = (1.0 - px1) * opxl[1] + px1 * opxh[1];
		tx2 = (1.0 - px2) * opxl[1] + px2 * opxh[1];
		ty1 = (1.0 - py1) * opyl[1] + py1 * opyh[1];
		ty2 = (1.0 - py2) * opyl[1] + py2 * opyh[1];
		/* update our perception of the view of original space */
		opxl[1] = tx1;
		opxh[1] = tx2;
		opyl[1] = ty1;
		opyh[1] = ty2;
		/* use this perception to get scale factor */
		xfrac = (opxh[1] - opxl[1])/(pxh[1] - pxl[1]);
		yfrac = (opyh[1] - opyl[1])/(pyh[1] - pyl[1]);
		/* now read in the figure */
		sclx[1] = 1.0/xfrac;
		scly[1] = 1.0/yfrac;
		scl[1] = MIN(sclx[1],scly[1]);
		do_clearscale(1 );
		gread(pompos[1].fname,pxl[1],pyl[1],
		    opxl[1],opyl[1],opxh[1],opyh[1],
		    1, sclx[1],scly[1]);
		do_box( pxl[1], pyl[1], pxh[1], pyh[1]);
		/* reset the plot values at the new corners */
		do_invscale(&per1,&vel1,x1,y1,1);
		do_invscale(&per2,&vel2,x2,y2,1);
		apxl[1] = per1;
		apyl[1] = vel1;
		apxh[1] = per2;
		apyh[1] = vel2;
		redisplay_pick();
	}
}

void do_unzoom(void)
{
	float xfrac, yfrac;
	/* turn off XOR */
	newpen(3000);
	do_clear(1);
	/* reset the pxl and pyl */
	/*
		pxl[1] = pompos[1].xl ;
		pxh[1] = pompos[1].xh ;
		pyl[1] = pompos[1].yl ;
		pyh[1] = pompos[1].yh ;
		*/
	apxl[1] = pompos[1].axl ;
	apxh[1] = pompos[1].axh ;
	apyl[1] = pompos[1].ayl ;
	apyh[1] = pompos[1].ayh ;
	opxl[1] = pompos[1].xl ;
	opxh[1] = pompos[1].xh ;
	opyl[1] = pompos[1].yl ;
	opyh[1] = pompos[1].yh ;
	xfrac = (opxh[1] - opxl[1])/(pxh[1] - pxl[1]);
	yfrac = (opyh[1] - opyl[1])/(pyh[1] - pyl[1]);
	sclx[1] = 1.0/xfrac;
	scly[1] = 1.0/yfrac;
	scl[1] = MIN(sclx[1],scly[1]);
	do_clearscale(1 );
	gread(pompos[1].fname,pxl[1],pyl[1],
	    opxl[1],opyl[1],opxh[1],opyh[1],
	    1, sclx[1],scly[1]);
	do_box( pxl[1], pyl[1], pxh[1], pyh[1]);
	/* redisplay picks */
	redisplay_pick();
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
}

void do_clear(int i )
{
	/* clear the graphic for the i'th figure */
	gclip("off", pxl[i],pyl[i],pxh[i],pyh[i]);
	newpen(0);
	shader(pxl[i], pyl[i], pxh[i], pyh[i], 0, 0, 0.01, 0.01);
	newpen(1);
}


void get_pom_output(void )
{
	/* get the output of pom96	*/
	FILE *pom96ctl;
	FILE *pom96dsp;
	float per;
	float fval;
	int i, j;
	char instr[512];


	/* GET PLOT CONTROL INFORMATION */
	if(( pom96ctl = fopen("pom96.ctl","r")) == NULL) return;
	fscanf(pom96ctl,"%d",&num_disp);
	fscanf(pom96ctl,"%s",instr);
	if(strncmp("XAXIS-FREQUENCY",instr,15) == 0)
		XAxisPeriod = NO ;
	else
		XAxisPeriod = YES;

	for (i=0;i<1;i++){
		fscanf(pom96ctl,"%f %f %f %f %e %e %e %e %s %s %s",
		    &pompos[i].xl, &pompos[i].yl,
		    &pompos[i].xh, &pompos[i].yh,
		    &pompos[i].axl, &pompos[i].ayl,
		    &pompos[i].axh, &pompos[i].ayh,
		    &pompos[i].xlnlg, &pompos[i].ylnlg,
		    &pompos[i].fname);
	}
	pompos[1].xl = pompos[0].xl;
	pompos[1].xh = pompos[0].xh;
	pompos[1].yl = pompos[0].yl;
	pompos[1].yh = pompos[0].yh;
	pompos[1].axl = pompos[0].axl;
	pompos[1].axh = pompos[0].axh;
	pompos[1].ayl = pompos[0].ayl;
	pompos[1].ayh = pompos[0].ayh;
	strcpy(pompos[1].xlnlg , pompos[0].xlnlg);
	strcpy(pompos[1].ylnlg , pompos[0].ylnlg);
	strcpy(pompos[1].fname,pompos[0].fname);
	/* get the number of filter periods */
	fscanf(pom96ctl,"%d",&num_period);
	/* get the filter periods */
	period = (float *)calloc(num_period,sizeof(float));
	for(i=0; i < num_period; i++){
		fscanf(pom96ctl,"%f",&fval);
		period[i] = fval;
	}
	fclose(pom96ctl);
	/* GET DISPERSION INFORMATION */
	if(( pom96dsp = fopen("pom96.dsp","r")) == NULL) return;
	/* allocate the array to hold the dispersion information */
	pomdsp = (struct disp *)calloc(num_disp+20,sizeof(struct disp));
	for(j=0;j<num_disp;j++){
		fscanf(pom96dsp,"%s %s %s %d %f %f %f %f %d",
		    &pomdsp[j].type,
		    &pomdsp[j].lvry, &pomdsp[j].cug,
		    &pomdsp[j].mode,
		             &per, &pomdsp[j].vel, 
		    &pomdsp[j].evel,
		    &pomdsp[j].amp, 
		    &pomdsp[j].symb);
	if(!XAxisPeriod)
		pomdsp[j].per = 1.0/per;
	else
		pomdsp[j].per = per;

	}
	fclose(pom96dsp);
}


void do_scale(float per,float vel,float *xval,float *yval,int i)
{
/*
	per	FLOAT	x-axis value (either period or frequency)
	vel	FLOAT	y-axis velocity value or spectral amplitude value
	*xval	FLOAT	returned screen x-coordinate
	*yval	FLOAT	returned screen y-coordinate
	i	INT	figure index, e.g., 0 for amp vs per, 1 for vel vs per
*/
	/* do the mapping between user coordinates and screen coordinates */
	float xfac, yfac;
	if(xlnlg[i] == LIN){
		xfac = (pxh[i] - pxl[i])
		    / (apxh[i] - apxl[i]);
		*xval = pxl[i] + (per - apxl[i])
		    * xfac;
	} else{
		xfac = (pxh[i] - pxl[i])
		    / log(apxh[i] / apxl[i]);
		*xval = pxl[i] + log(per / apxl[i])
		    * xfac;
	}
	if(ylnlg[i] == LIN){
		yfac = (pyh[i] - pyl[i])
		    / (apyh[i] - apyl[i]);
		*yval = pyl[i] + (vel - apyl[i])
		    * yfac ;
	} else{
		yfac = (pyh[1] - pyl[1])
		    / log(apyh[i] / apyl[i]);
		*yval = pyl[i] + log(vel / apyl[i])
		    * yfac;
	}
}

void do_invscale(float *per,float *vel,float xval,float yval,int i)
{
	/* do the mapping from screen coordinates to user coordinates */
	float xfac, yfac;
	if(xlnlg[i] == LIN){
		xfac = (pxh[i] - pxl[i])
		    / (apxh[i] - apxl[i]);
		*per = apxl[i] + (xval - pxl[i])/xfac ;
	} else{
		xfac = (pxh[i] - pxl[i])
		    / log(apxh[i] / apxl[i]);
		*per = apxl[i] * 
		    exp((xval - pxl[i])/xfac);
	}
	if(ylnlg[i] == LIN){
		yfac = (pyh[i] - pyl[i])
		    / (apyh[i] - apyl[i]);
		*vel = apyl[i] + (yval - pyl[i])/yfac ;
	} else{
		yfac = (pyh[i] - pyl[i])
		    / log(apyh[i] / apyl[i]);
		*vel = apyl[i] * 
		    exp((yval - pyl[i])/yfac);
	}
}


int do_mode(void)
{
	int mode;
	int i, cmd, nmd, ii ;
	char c[2];
	float xv, yv;
	/* output a mode menu for mode selection */
	show_menu(1.0, 0.25, md, sizeof(md), &nmd);
	/* place current mode at top of page */
	mode = -1 ;
	cmd = -1;
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nmd ; i++){
			if(inside(xv,yv,md[i].xl,md[i].yl,md[i].xh,md[i].yh))
			{
				cmd = md[i].action;
				gmesg(md[i].str);
				ii = i;
				iMode = ii;
				break;
			}
		}
		if(cmd > 0){
			mode = cmd -1;
		}
	}
	/* clear the window */
	gclip("off", XL, YL, XH, YH);
	newpen(0);
	shader(md[0].xl,md[0].yl,
	    md[nmd-1].xh,md[nmd-1].yh, 0, 0, 0.01, 0.01);
	/* put the mode title at the top */
	do_showmd(mode);
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	return (mode);
}

void redisplay_pick()
{
	float per, vel;
	float xval, yval;
	int i;
	/* redisplay picks */
	/* first show the velocities */
	/* clear the scale */
	do_axes(1);
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	/* turn on XOR */
	newpen(3001);
	newpen(Pick_Color);
	for(i=0; i< num_disp ; i++){
		per = pomdsp[i].per;
		vel = pomdsp[i].vel;
		if(pomdsp[i].mode >= 0){
			do_scale(per,vel,&xval,&yval,1);
			if(inside(xval,yval,pxl[1],pyl[1],pxh[1],pyh[1])){
				/* put up the velocity  point */
				gsolid(xval,yval,0.02*scl[1],pomdsp[i].symb);
			}
		}
	}
	/* turn off XOR */
	newpen(3000);

}

void do_start(int *nm, int menuon)
{
	int i;
	float xfrac, yfrac;
	get_pom_output();
	for(i=0;i<3;i++){
		if(strncmp(pompos[i].xlnlg,"lin",3) == 0)
			xlnlg[i] = LIN;
			else
			xlnlg[i] = LOG;
		if(strncmp(pompos[i].ylnlg,"lin",3) == 0)
			ylnlg[i] = LIN;
			else
			ylnlg[i] = LOG;
		/*
				pxl[i] = pompos[i].xl ;
				pyl[i] = pompos[i].yl ;
				pxh[i] = pompos[i].xh ;
				pyh[i] = pompos[i].yh ;
				*/
		apxl[i] = pompos[i].axl;
		apyl[i] = pompos[i].ayl;
		apxh[i] = pompos[i].axh;
		apyh[i] = pompos[i].ayh;
		opxl[i] = pompos[i].xl ;
		opyl[i] = pompos[i].yl ;
		opxh[i] = pompos[i].xh ;
		opyh[i] = pompos[i].yh ;
		/* use this perception to get scale factor */
		xfrac = (opxh[i] - opxl[i])/(pxh[i] - pxl[i]);
		yfrac = (opyh[i] - opyl[i])/(pyh[i] - pyl[i]);
		/* now read in the figure */
		sclx[i] = 1.0/xfrac;
		scly[i] = 1.0/yfrac;
		scl[i] = MIN(sclx[i],scly[i]);
	}
	/* initialize plot scale parameters */
	/*
		printf("%d\n",num_disp);
		for(i=0;i<3;i++)
		printf("%f %f %f %f %e %e %e %e %s %s %s\n",
			pompos[i].xl, pompos[i].yl,
			pompos[i].xh, pompos[i].yh,
			pompos[i].axl, pompos[i].ayl,
			pompos[i].axh, pompos[i].ayh,
			pompos[i].xlnlg, pompos[i].ylnlg,
			pompos[i].fname);
	
		for(i=0;i<num_disp;i++)
		printf("%s %s %s %d %f %f %f %f %e %f %f %f %f %d %d %f %s %s %s %d %d %d %d\n",
			pomdsp[i].type, 
			pomdsp[i].lvry, pomdsp[i].cug,
			pomdsp[i].mode,
			pomdsp[i].per, pomdsp[i].vel, 
			pomdsp[i].evel, 
			pomdsp[i].dist, pomdsp[i].amp, 
			pomdsp[i].lat1, pomdsp[i].lon1, 
			pomdsp[i].lat2, pomdsp[i].lon2,
			pomdsp[i].ictl, pomdsp[i].symb,
			pomdsp[i].instper,
			pomdsp[j].comment,
			pomdsp[j].kstnm,
			pomdsp[j].kcmpnm,
			pomdsp[j].nzyear,
			pomdsp[j].nzjday,
			pomdsp[j].nzhour,
			pomdsp[j].nzmin);
	*/

	/* put up the plots */
	for(i=1;i<2;i++){
		gread(pompos[i].fname,pxl[i],pyl[i],
		    pompos[i].xl,pompos[i].yl,
		    pompos[i].xh,pompos[i].yh,
		    1, sclx[i],scly[i]);
		do_box( pxl[i], pyl[i], pxh[i], pyh[i]);
	}
	/* now put up the plot scales */
	do_axes(1);
	if(menuon)show_menu(0.5, 0.6, m, sizeof(m), nm);
	Mode = -1;
}

void do_restart(int nm)
{
	/* clear the graphics area and do not redisplay the menu */
	gclip("off", XL, YL, XH, YH);
	/* turn off XOR and erase entire display area */
	newpen(3000);
	newpen(0);
	shader(XL, YL, XH, YH, 0, 0, 0.01, 0.01);
	Mode = -1;
	iMode = -1;
	InPick = OFF;
	InLine = OFF;
	do_showmd(Mode);
	do_showpk();
	mgwrtxt(0.5,7.7,wavetype,1,2);
	newpen(1);
	do_start(&nm, OFF);
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	/* should mean reset all picking */
}


void do_box(float xl,float yl,float xh,float yh)
{
	/* draw a box */
	plot(xl,yl,3);
	plot(xh,yl,2);
	plot(xh,yh,2);
	plot(xl,yh,2);
	plot(xl,yl,2);
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

void do_axes(int i )
{
	float xx, yy, xlen, ylen, xlow, xhgh, ylow, yhgh;
	xx = pxl[i];
	yy = pyl[i];
	xlen = pxh[i] - pxl[i];
	ylen = pyh[i] - pyl[i];
	xlow = apxl[i];
	ylow = apyl[i];
	xhgh = apxh[i];
	yhgh = apyh[i];
	newpen(1);
	if(xlnlg[i] == LIN){
		if(XAxisPeriod)
		dolinx(xx,yy,xlen,xhgh,xlow,0.10,OFF,OFF,ON,10,"Period (s)");
		else
		dolinx(xx,yy,xlen,xhgh,xlow,0.10,OFF,OFF,ON,14,"Frequency (Hz)");
	} else {
		if(XAxisPeriod)
		dologx(xx,yy,xlen,xhgh,xlow,0.10,OFF,OFF,ON,10,"Period (s)");
		else
		dologx(xx,yy,xlen,xhgh,xlow,0.10,OFF,OFF,ON,14,"Frequency (Hz)");
	}
	if(ylnlg[i] == LIN){
		doliny(xx , yy ,ylen,yhgh,ylow,0.07,ON,ON,ON,0," ");
	} else {
		dology(xx , yy ,ylen,yhgh,ylow,0.07,ON,ON,ON,0," ");
	}
	if ( i == 1){
		gright(pxl[i]+xlen,pyh[i]+0.1,0.15,"C(km/sec)",0.0);
	}
}

void do_clearscale(int i )
{
	float xl, yl, xh, yh;
	/* clear y -axis */
	if(i == 0)
		xl = 0.0;
		else
		xl = pxh[i-1];
	xl = xl + 0.05;
	yl = pyl[i] - 0.10;
	xh = pxl[i];
	yh = pyh[i] + 0.10;
	newpen(0);
	shader(xl, yl, xh, yh, 0, 0, 0.01, 0.01);
	/* clear y-axis on right */
	xl = pxh[i];
	xl = xl + 0.05;
	yl = pyl[i] - 0.10;
	xh = pxl[i] + 0.10;
	yh = pyh[i] + 0.10;
	newpen(0);
	shader(xl, yl, xh, yh, 0, 0, 0.01, 0.01);
	/* clear x-axis */
	xl = pxl[i];
	xh = pxh[i] + 0.05;
	yl = pyl[i] - 0.40;
	yh = pyl[i];
	shader(xl, yl, xh, yh, 0, 0, 0.01, 0.01);
	newpen(1);
}

void do_showpk()
{
	gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
	mgwrtxt(3.0,7.7,"   Action   ",1,2);
	if(InPick == ON)
		mgwrtxt(3.0,7.4," Picking On ",1,1);
	else if(InLine == ON)
		mgwrtxt(3.0,7.4,"Auto Picking",1,1);
	else
		mgwrtxt(3.0,7.4," Picking Off",1,1);
}
void do_showmd(int mode)
{
	mgwrtxt(2.0,7.7,"Mode\0",1,2);
	if(mode < 0)
		mgwrtxt(2.0,7.4,"None\0",1,1);
	else
		mgwrtxt(2.0,7.4,md[iMode].str,1,1);
}


void doout(int *ret)
{
	int j, myret;
	int i, ii;
	float xv, yv;
	char c[2];
	int cmd, nmy;
	float per;
	FILE *pom96odsp;
	char ofname[20];
	if(wave == 1)
		strcpy(ofname,"love.dsp");
	else if(wave == 2)
		strcpy(ofname,"rayl.dsp");
	else
		strcpy(ofname,"unkn.dsp");
	fprintf(stderr,"%s\n",ofname);
	/* ask whether to also save the file as STACOMP.dsp */
	newpen(3000);
	mgwrtxt(0.5,0.25,"Save as",1,1);
	mgwrtxt(1.5,0.25,ofname,1,1);
	show_menu(3.5, 0.15, my, sizeof(my), &nmy);
	/* place current mode at top of page */
	myret = -1 ;
	cmd = -1;
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nmy ; i++){
			if(inside(xv,yv,my[i].xl,my[i].yl,my[i].xh,my[i].yh))
			{
				cmd = my[i].action;
				gmesg(my[i].str);
				ii = i;
				iMode = ii;
				break;
			}
		}
		if(cmd > 0){
			if(cmd == 1)
				myret = 1;
			else
				myret = -1;
		}
	}
	/* if myret == 1 create the new file else use the pom96.ods */
	if(myret > 0){
		if(( pom96odsp = fopen(ofname     ,"w")) == NULL) return;
		if(*ret != 2)
			*ret = myret;
	} else {
		if(( pom96odsp = fopen("pom96.ods","w")) == NULL) return;
	}
	for(j=0;j<num_disp;j++){
		if(pomdsp[j].mode >= 0){
			if(!XAxisPeriod)
				per = 1.0/pomdsp[j].per;
			else
				per = pomdsp[j].per;
			fprintf(pom96odsp,"SURF96 %1s %1s X %2d %11.4g %11.5f %11.5f %10.4f \n",
				pomdsp[j].lvry, pomdsp[j].cug,
				pomdsp[j].mode,
				per, 
				pomdsp[j].vel, 
				pomdsp[j].evel,
				pomdsp[j].amp);
		}
	}
	fclose(pom96odsp);
}

int do_matchprep(char *fname)
{
	/* interpolate the instantaneous period and the chosen spectra, according to
		the period array. Also select the mode for output
	*/
	int i ;
	int nm;
	gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
	ModeSelect = -1;
	newpen(1);
	do_axes(1);
	/* now put up a menu and then an event loop */
	for(i=0;i<2;i++){
		do_box( pxl[i], pyl[i], pxh[i], pyh[i]);
	}
	newpen(1);
	do_showpick1();
	/* now show the desired periods */
	do_showperscl(0);
	do_showperscl(1);
	show_menu(0.5, 1.1, mm, sizeof(mm), &nm);
	return(get_action1(nm,fname));
}

char *pathname;
int get_action1(int nm, char *fname)
{
	int cmd, i, ii;
	int ret;
	float xv, yv;
	float xl, yl, xh, yh;
	char c[2];
	int InDel = OFF;
	int InAdd = OFF;
	int oModeSelect ;
	/* annotate top */
	oModeSelect = 1;
	mgwrtxt(2.0,7.4,"None\0",2, 1);
	Mode = -1;
	do_showmd(ModeSelect);
	mgwrtxt(0.5,7.7,wavetype,1,2);
	/* prohibit writing at the bottom */
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	gcursor("Arrow");
	for(; ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
		/* loop on commands */
		for(i=0 ; i < nm ; i++){
			if(inside(xv,yv,mm[i].xl,mm[i].yl,mm[i].xh,mm[i].yh))
			{
				cmd = mm[i].action;
				gmesg(mm[i].str);
				ii = i;
				break;
			}
		}
		if(cmd > 0){
printf( "cmd: %d\n",cmd);
			switch(cmd){
			case ZOOM:	/* Zoom */
				draw_button(0.5+ii*1.125,0.6 ,
				    mm[ii].str,&xl,&yl,&xh,&yh,ON ,mm[ii].lstrmx);
				gcursor("Arrow");
				do_zoomm();
				gcursor("Arrow");
				draw_button(0.5+ii*1.125,0.6 ,
				    mm[ii].str,&xl,&yl,&xh,&yh,OFF,mm[ii].lstrmx);
				do_clear(1);
				do_clearscale(1);
				do_box( pxl[1], pyl[1], pxh[1], pyh[1]);
				do_axes(1);
				do_showpick1();
				do_showmode1(ModeSelect);
				do_showperscl(1);
				break;
			case UNZOOM:	/* UnZoom */
				InPick = OFF;
				InLine = OFF;
				do_unzoomm();
				do_clear(1);
				do_clearscale(1);
				do_box( pxl[1], pyl[1], pxh[1], pyh[1]);
				do_axes(1);
				do_showpick1();
				do_showmode1(ModeSelect);
				do_showperscl(1);
				break;
			case MODE:	/* Mode */
				oModeSelect = ModeSelect;
				ModeSelect = do_mode();
				do_showmode1(oModeSelect);
				do_showmode1(ModeSelect);
				do_showmd(ModeSelect);
				break;
			case INTERP:
				printf("Interpolate Periods within current box\n");
				break;
			case ADD:
				if(InAdd == ON) {
					InAdd = OFF ;
					InDel = OFF;
					draw_button(0.5+ii*1.125,1.1 ,
				    		mm[ii].str,&xl,&yl,&xh,&yh,OFF,mm[ii].lstrmx);
				} else {
					InAdd = ON;
					InDel = OFF;
					draw_button(0.5+ii*1.125,1.1 ,
				    		mm[ii].str,&xl,&yl,&xh,&yh,ON ,mm[ii].lstrmx);
				}
				printf("Add a fictitious point for Match \n");
				break;
			case DELETE:
				if(InDel == ON) {
					InDel = OFF ;
					InAdd = OFF ;
					InPick = OFF;
					draw_button(0.5+ii*1.125,1.1 ,
				    		mm[ii].str,&xl,&yl,&xh,&yh,OFF,mm[ii].lstrmx);
				} else {
					draw_button(0.5+ii*1.125,1.1 ,
				    		mm[ii].str,&xl,&yl,&xh,&yh,ON ,mm[ii].lstrmx);
					InDel = ON;
					InAdd = OFF ;
					InPick = ON;
				}
				printf("Delete a real/fictitious dispersion for Match \n");
				break;
			case RETURN:	/* Return */
				printf("Return\n");
				gframe(1);
				return(0);
			case MATCH:	/* Match */
				while(ModeSelect < 0 ){
					gmesg("No mode chosen");
					oModeSelect = ModeSelect;
					ModeSelect = do_mode();
					do_showmode1(oModeSelect);
					do_showmode1(ModeSelect);
					do_showmd(ModeSelect);
				}
				doout1(ModeSelect);

#ifdef MSDOS
sprintf(ostr,"sacpom96.exe %s -F %s -D disp.d -AUTO",pathname, fname);
#else
fprintf(stderr,"ModeSelect: %d pathname: %s\n",ModeSelect, pathname);
        pathname = (char *)my_pathfind(getenv("PATH"), "sacpom96", "rx");
fprintf(stderr,"ModeSelect: %d pathname: %s\n",ModeSelect, pathname);
sprintf(ostr,"%s -F %s -D disp.d -AUTO",pathname, fname);
#endif
fprintf(stderr,"ostr: %s\n",ostr);
				printf("run sacpom96 Mode: %d\n",Mode);
					ret = system(ostr);
					if(ret >= 0 ){
						/* success */
						return(2);
					}
				break;
			case EXIT:
				printf("Exit\n");
				gframe(1);
				return(-1);
			default:
				break;

			}
		} else {
			if(InPick == ON){
				if(inside(xv,yv,pxl[1]-0.1,pyl[1],pxh[1]+0.1,pyh[1]))
				{
					do_pick(xv, yv, 1);
				}
			}
		}
	}
} 

void do_zoomm(void)
{
	char c[2];
	float xv1, xv2, yv1, yv2;
	float xfrac, yfrac;
	float x1,x2,y1,y2;
	float tx1,tx2,ty1,ty2;
	float per1, vel1;
	float per2, vel2;
	float px1,px2,py1,py2;
	/* turn off XOR */
	newpen(3000);
	curaxy(&xv1, &yv1, c);
	if( !inside(xv1,yv1, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1]))return;
	gcursor("Box");
	curaxy(&xv2, &yv2, c);
	gcursor("Arrow");
	if(xv1==xv2 && yv1 == yv2) return;  /* need two points for a line */
	/* 
	We work with three coordinate systems, each consisting of two
	opposite corners of a view rectangle

	(pxl,pyl) -> (pxh,pyl)   is display screen which does not change
	(opxl,opyl) -> (opxh,opyh) is the mapped portion of the original plot
	(apxl,apyl) -> (apxh,apyl) are user coordinates of the original space

	To do the mapping, use a simple interpolation

	V(p) = (1-p)V1 + pV2   where  0 <= p <= 1
	*/
	if( inside(xv1,yv1, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1])
	    && inside(xv2,yv2, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1]) ) {
		/* display all current picks to turn them off */
		/* since the point may extend across clip region ? */
/*
		redisplay_pick();
*/
		do_clear(1);
		/* these are screen coordinates */
		x1 = MIN(xv1,xv2) ;
		x2 = MAX(xv1,xv2) ;
		y1 = MIN(yv1,yv2) ;
		y2 = MAX(yv1,yv2) ;
		px1 = pinterp(x1,pxl[1],pxh[1]);
		px2 = pinterp(x2,pxl[1],pxh[1]);
		py1 = pinterp(y1,pyl[1],pyh[1]);
		py2 = pinterp(y2,pyl[1],pyh[1]);
		/* now estimate the coordinates in the
					original plot space */
		tx1 = (1.0 - px1) * opxl[1] + px1 * opxh[1];
		tx2 = (1.0 - px2) * opxl[1] + px2 * opxh[1];
		ty1 = (1.0 - py1) * opyl[1] + py1 * opyh[1];
		ty2 = (1.0 - py2) * opyl[1] + py2 * opyh[1];
		/* update our perception of the view of original space */
		opxl[1] = tx1;
		opxh[1] = tx2;
		opyl[1] = ty1;
		opyh[1] = ty2;
		/* use this perception to get scale factor */
		xfrac = (opxh[1] - opxl[1])/(pxh[1] - pxl[1]);
		yfrac = (opyh[1] - opyl[1])/(pyh[1] - pyl[1]);
		sclx[1] = 1.0/xfrac;
		scly[1] = 1.0/yfrac;
		scl[1] = MIN(sclx[1],scly[1]);
		/* reset the plot values at the new corners */
		do_invscale(&per1,&vel1,x1,y1,1);
		do_invscale(&per2,&vel2,x2,y2,1);
		apxl[1] = per1;
		apyl[1] = vel1;
		apxh[1] = per2;
		apyh[1] = vel2;
	}
}

void do_unzoomm(void)
{
	float xfrac, yfrac;
	/* turn off XOR */
	newpen(3000);
	do_clear(1);
	/* reset the pxl and pyl */
	/*
		pxl[1] = pompos[1].xl ;
		pxh[1] = pompos[1].xh ;
		pyl[1] = pompos[1].yl ;
		pyh[1] = pompos[1].yh ;
		*/
	apxl[1] = pompos[1].axl ;
	apxh[1] = pompos[1].axh ;
	apyl[1] = pompos[1].ayl ;
	apyh[1] = pompos[1].ayh ;
	opxl[1] = pompos[1].xl ;
	opxh[1] = pompos[1].xh ;
	opyl[1] = pompos[1].yl ;
	opyh[1] = pompos[1].yh ;
	xfrac = (opxh[1] - opxl[1])/(pxh[1] - pxl[1]);
	yfrac = (opyh[1] - opyl[1])/(pyh[1] - pyl[1]);
	sclx[1] = 1.0/xfrac;
	scly[1] = 1.0/yfrac;
	scl[1] = MIN(sclx[1],scly[1]);
	/* redisplay picks */
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
}

/* display a particular mode - this will permit XOR of
   switching the entire mode on/off and also the interactive picking
*/
void do_showmode1(int mode){
	int j;
	float xval, yval, vel, per;
	if(mode < 0)return;
	for(j=0;j<num_disp;j++){
		if(pomdsp[j].mode == mode  ){
			if(!XAxisPeriod)
				per = 1.0/pomdsp[j].instper;
			else
				per = pomdsp[j].instper;

			vel = pomdsp[j].vel;
			/* display velocity value */
			/* turn on XOR */
			newpen(3001);
			newpen(1);
			do_scale(per,vel,&xval,&yval,1);
			if(inside(xval,yval,pxl[1],pyl[1],pxh[1],pyh[1])){
				gsolid(xval,yval,0.035*scl[1],pomdsp[j].mode);
			}
			/* turn off XOR */
			newpen(3000);
		}
	}
}

/* display picked values */
void do_showpick1(void)
{
	int j;
	float xval, yval, vel, per;
	for(j=0;j<num_disp;j++){
		if(pomdsp[j].mode >= 0){
			if(!XAxisPeriod)
				per = 1.0/pomdsp[j].instper;
			else
				per = pomdsp[j].instper;

			vel = pomdsp[j].vel;
			/* display velocity value */
			do_scale(per,vel,&xval,&yval,1);
			if(inside(xval,yval,pxl[1],pyl[1],pxh[1],pyh[1])){
				gsolid(xval,yval,0.02*scl[1],pomdsp[j].mode);
			}
		}
	}
}

void do_showperscl(int i)
{
	int j;
	float per, vel, xval,yval;
	newpen(2);
	for(j=0 ; j < num_period; j++){
		if(!XAxisPeriod)
			per = 1.0/period[j];
		else
			per = period[j];
		vel = pomdsp[0].vel;
		do_scale(per,vel,&xval,&yval,i);
		if(xval >= pxl[i] && xval <= pxh[i]){
			plot(xval,pyl[i]    ,3);
			plot(xval,pyl[i]+0.2,2);
			plot(xval,pyh[i]    ,3);
			plot(xval,pyh[i]-0.2,2);
		}
	}
	newpen(1);
}

void do_resetzoom(void)
{
float xfrac, yfrac;
	apxl[1] = pompos[1].axl ;
	apxh[1] = pompos[1].axh ;
	apyl[1] = pompos[1].ayl ;
	apyh[1] = pompos[1].ayh ;
	opxl[1] = pompos[1].xl ;
	opxh[1] = pompos[1].xh ;
	opyl[1] = pompos[1].yl ;
	opyh[1] = pompos[1].yh ;
	xfrac = (opxh[1] - opxl[1])/(pxh[1] - pxl[1]);
	yfrac = (opyh[1] - opyl[1])/(pyh[1] - pyl[1]);
	sclx[1] = 1.0/xfrac;
	scly[1] = 1.0/yfrac;
	scl[1] = MIN(sclx[1],scly[1]);
}

void doout1(int mode)
{
	int j;
	FILE *pom96mat;
	FILE *pom96matp;
	if(( pom96mat = fopen("disp.d","w")) == NULL) return;
	if(( pom96matp = fopen("disp.dp","w")) == NULL) return;
	for(j=0;j<num_disp;j++){
		if(pomdsp[j].mode == mode){
		fprintf(pom96mat,"SURF96 %1s %1s X %2d %11.4g %11.4g %11.4g \n",
			pomdsp[j].lvry, 
			pomdsp[j].cug,
			pomdsp[j].mode,
			pomdsp[j].instper, 
			pomdsp[j].vel,
			pomdsp[j].evel); 
		fprintf(pom96matp,"SURF96 %1s %1s X %2d %11.4g %11.4g %11.4g \n",
			pomdsp[j].lvry, 
			pomdsp[j].cug,
			pomdsp[j].mode,
			pomdsp[j].per, 
			pomdsp[j].vel,
			pomdsp[j].evel); 

		}
	}
	fclose(pom96mat);
	fclose(pom96matp);
}
