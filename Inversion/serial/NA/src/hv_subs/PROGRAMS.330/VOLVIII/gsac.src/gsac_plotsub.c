/* NOTES
 * before ploting traces in a group determine the extremes of the absolute
 * B and E so that the time scale can be defined
 *
 * When doing ppk with perplot - the overscale scale will be determined
 * by perplot. Note however that the clip region may be smaller if
 * on the last page

 * changes: 06 JAN 2005 - changed x-positioning to doubles because of
 *                        problems in plotting very long time series -
 *                        rethink the loops to be smarter and not to clip
 *                        entire trace
 *
 *                        Added Mark time
 *
 *                        Introduced ilow, ihgh so that entire trace is
 *                        now plotted with pen up for a window - this
 *                        makes the plot files smaller
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
extern int markt_on ;
extern int markt_doo ;
extern int markt_dod ;
extern int markt_numvel ;
extern float markt_dist ;
extern double markt_o ;
#define NUMMARKT 11
extern float markt_vel[NUMMARKT];
extern void gsac_plot_fileid(float x0,float y0,float xlen,float ylen, int k);


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


static char *title[] = {
	"Frequency (Hz)" , "Time (s)"  
	} ;



struct plt_ctl *pmap = (struct plt_ctl *)NULL; 

void gsac_show_pick(float x0, float y0, float xlen, float ylen, int k, float yy0, float dy, float twin, int pabsolute);
void gsac_show_markt(float x0, float y0, float xlen, float ylen, int k, float yy0, float dy, float twin, int pabsolute);
void gsac_show_plotpk(float x0, float y0, float xlen, float ylen, int ns, int ne, float ts, float te, float dy, int ntrc, double twin,int numperframe, int pabsolute, float *uy, float *ly, int setup, int ylimctrl, float ylimlow , float ylimhigh, int overlay);
void showdec(float xl, float yl, float xh, float yh, int inc);
void dogrid(void);

extern void XviG_Flush();

char timstr[32];


void gsac_show_plotpk(float x0, float y0, float xlen, float ylen, int ns, int ne, float tb, float te, float dy, int ntrc, double twin,int numperframe, int pabsolute_in, float *uy, float *ly, int setup, int ylimctrl, float ylimlow, float ylimhigh, int overlay)
{
	int i,k, kkk, kk,inc;
	float yy0, yyy0;
	double dx;
	double xx;
	float  yy;
	double txx;
	float tdepmax, tdepmin;
	int indmax, indmin;
	float ht;
	int npts;
	float delta, depmax, depmin;
	float u;
	double twinl, twinh;
	double twinll, twinhh;
	int nwind;			/* number of data points in  window */

	float floor;			/* lower bound for log y plot */
	float ceiling;

	float uxcen, uxmul;
	int mqdp;
	int ilow, ihgh;			/* rough indices of visible plot 
						added 14 JAN 2005 for 
						faster/smaller plot*/
	double tzs, tze;
	int isxtime;			/* is the axis time or frequency */
	int pabsolute;			/* local version - overrides user
					default if trace is IXY == freq */

	uxcen = gsac_control.uxcen;
	uxmul = gsac_control.uxmul;
	/* beginning of trace plot */

	*ly =  1.0e+38;
	*uy = -1.0e+38;

	ht = 0.10;
	
	if(gsac_control.grid)
		dogrid();

	newpen(1);
	gclip("off", x0, y0, x0+xlen, y0+dy);
/*
	gbox(x0, y0, x0+xlen, y0+ylen);
*/

	/* do the bottom axis */
	/* put trace plot separate after the setup  to use the
	 * struct parameters */
	/* for absolute plot put in the B b +TWIN e.g, */
	twinl = twin*(uxcen - 0.5/uxmul);
	twinh = twin*(uxcen + 0.5/uxmul);
	/* first trace defines Time or Frequency
		XLOG is permitted only for IXY 
		also IXY implies pabsolute, no markt, eventually reduced ppk */
	if(sacdata[0].sachdr.ihdr[15] == 1){
		isxtime = YES ;
		pabsolute = pabsolute_in ;
	} else {
		isxtime = NO ;
		pabsolute = YES ;
	}
	/* last trace define the scale */
	if(pabsolute == YES) {
		twinll = twinl + sacdata[ne-1].tzbegx - sacdata[ne-1].tzbeg ;
		twinhh = twinh + sacdata[ne-1].tzbegx - sacdata[ne-1].tzbeg ;
		tzs = sacdata[ne-1].tzbeg + twinll;
		tze = sacdata[ne-1].tzbeg + twinhh - twinll;
		/* by definition tzs, tze are > 0 for YEAR > 1970 1 00 etc */
		printtimestr(sacdata[ne-1].tzbeg,timstr);
		timstr[31]='\0';
	} else {
		twinll = twinl;
		twinhh = twinh;
		if(isxtime == YES)
			sprintf(timstr,"Time relative to %s %f\n", 
				gsac_control.xlimkey[0],
				gsac_control.xlimoff[0]);
	}

	/* YLIM - if ylimctrl == YLIM_ALL then we must search through all
	 * traces to find the extreme limits to get largest depmax and least
	 * depmin . Recall that
	 * YLIM_OFF means automatic control
	 * YLIM_ALL means all traces have same scaling
	 * YLIM_USR means that user has specified the scaling limits to be
	 * applied to all traces
	 * THESE VALUES ARE FIXED
	 * */
	if(ylimctrl == YLIM_ALL){
		depmax = -1.0e+37;
		depmin =  1.0e+37;
		for ( kkk=ns ; kkk < ne ; kkk++){
			k = sortptr[kkk];
			depmin = MIN(sacdata[k].sachdr.rhdr[H_DEPMIN],depmin);
 			depmax = MAX(sacdata[k].sachdr.rhdr[H_DEPMAX],depmax);
		}
	} else if (ylimctrl == YLIM_USR){
			depmin = MIN(ylimlow, ylimhigh);
			depmax = MAX(ylimlow, ylimhigh);
	}
	/* set the background */
	if(gsac_control.background == YES &&
		gsac_control.background_color >= 0){
		newpen(gsac_control.background_color);
		if(overlay == YES){
		shader(x0,y0,x0+xlen,y0+ylen,0,0,0.01,0.01);
		} else {
		shader(x0,y0+ylen-(ne-ns)*dy,x0+xlen,y0+ylen,0,0,0.01,0.01);
		}
		newpen(1);
	}
	/* NOW PROCESS TRACES */
	for ( kkk=ns ; kkk < ne ; kkk++){
		k = sortptr[kkk];
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		delta  = sacdata[k].sachdr.rhdr[H_DELTA];
		tb = sacdata[k].sachdr.rhdr[H_B];
		te = sacdata[k].sachdr.rhdr[H_E];
		dx = xlen*delta/twin;


		/* use depmax depmin */
		kk = kkk%numperframe;
		if(overlay == YES){
			yy0 = y0 ;
			yyy0 = yy0 - kk*1.5*ht ;
		} else {
			yy0 = y0 + (numperframe -1  - kk )*dy;
			yyy0 = yy0  ;
		}

		nwind = 1 + (twinh - twinl)/delta;
		if(gsac_control.plotdevice==WIN){
			if(gsac_control.qdp < 0)
				inc = 1;
			else if(gsac_control.qdp > 0)
				inc = gsac_control.qdp;
			else {
				/* automatic determination */
				mqdp = MIN(gsac_control.XmaxDev,gsac_control.YmaxDev);
				inc = MAX(1, npts/mqdp);
			}
		} else {
			/* never do a QDP for an CALPLOT external plot */
			inc = 1;
		}

		/* adjust the decimation */
		if(nwind <= 4000)
			inc = 1;

		pmap[kk].abstime = pabsolute;
		pmap[kk].xl = x0;
		pmap[kk].xh = x0 +xlen;;
		pmap[kk].yl = yy0;
		pmap[kk].yh = yy0 +dy;
		pmap[kk].k  = k;
		pmap[kk].n  = kk;
		pmap[kk].tfirst = sacdata[k].tzbegx;
		pmap[kk].tlast  = sacdata[k].tzendx;
		pmap[kk].ifirst = 0;
		pmap[kk].ilast  = npts -1;
		pmap[kk].npts  = npts;
		pmap[kk].xlen = xlen;
		pmap[kk].delta = delta;
		/* not used but also float = double */
		pmap[kk].ymult = 1.0;
		if(setup == YES){
			if(ylimctrl == YLIM_OFF){
				/* define the depmax depmin to use if
				 * not already defined by YLIM_ALL or 
				 * YLIM_USR */
				depmin = sacdata[k].sachdr.rhdr[H_DEPMIN];
				depmax = sacdata[k].sachdr.rhdr[H_DEPMAX];
				depmax += 0.05*fabs(depmax);
				depmin -= 0.05*fabs(depmin);
				pmap[kk].depmax = depmax;
				pmap[kk].depmin = depmin;
				/* to implement the automatic trace scaling
				 * for the current window, we modify the
				 * pmap[kk].uymul; Note that since this is 
				 * based on depmax, depmin, uymul >= 1 
				 * */
				pmap[kk].uymul = 1.0;
			} else {
				/* use previously set values */
				pmap[kk].depmax = depmax;
				pmap[kk].depmin = depmin;
				pmap[kk].uymul = 1.0;
			}
		}

		if(*uy < pmap[kk].yh)
			*uy = pmap[kk].yh ;
		if(*ly > pmap[kk].yl)
			*ly = pmap[kk].yl ;
		if(pabsolute == YES) {
			u = (float)(sacdata[k].tzref + tb  
				- gsac_control.begminx)/twin;
			xx = x0 + xlen*(uxmul*(u-uxcen) +0.5 );
		} else {
			u = -(float)(sacdata[k].tzbegx - sacdata[k].tzbeg)/twin;
			xx = x0 + xlen*(uxmul*(u-uxcen) +0.5 );
		}
		/* define range of data that fills within the clip space
			This is done to reduce the number of non-plotable
			move calls, which significantly reduces the
			size of the PLOT file and also speeds up the
			on-line plot */
		ilow = MAX(0 ,(x0-xx)/(uxmul*dx) -1 );
		ihgh = MIN(npts , (x0+xlen-xx)/(uxmul*dx) +1 );
		/* to implement the autoscaling within the YLIM_OFF
		 * regime I must look at the points within the current
		 * plot window and then readjust the scaling
		 * */
		if(ylimctrl == YLIM_OFF ){
			txx = xx + (ilow)*uxmul*dx ;
			tdepmax = -1.0e+37;
			tdepmin =  1.0e+37;
			for(i=ilow ; i < ihgh ; i+= inc ){
				if(txx >= x0 && txx <= x0+xlen){
				tdepmin = MIN(sacdata[k].sac_data[i],tdepmin);
 				tdepmax = MAX(sacdata[k].sac_data[i],tdepmax);
				}
				txx += uxmul*dx*inc;
			}
			if (tdepmax == tdepmin) {
				if(tdepmax == 0){
					tdepmax = 1.0;
				} else {
					tdepmin = -tdepmax;
				}
			}
			depmax = tdepmax + 0.05*fabs(tdepmax);
			depmin = tdepmin - 0.05*fabs(tdepmin);
		}
		if(gsac_control.ylim_rsc == YES){
			depmax /= pmap[kk].uymul;
			depmin /= pmap[kk].uymul;
		}
		/* up to here xx is the absolute position of first point
		 * with clipping we only plot in the range [ 0 , xlen ]
		 * so with inverse we can map the window back to the sample
		 * or absolute time */
		/* never divide by 0 */
		if(gsac_control.xgrid == YES)
		dolnxgrid(x0,yy0,yy0+dy,xlen,twinll,twinhh,0.10, 
			YES, gsac_control.xgrid_color, gsac_control.xgrid_type,gsac_control.xgrid_minor);
		dolinx(x0,yy0+dy,xlen,twinll,twinhh,0.10, NO, NO, NO, 1, " ");
		if(kkk == ne -1){
			dolinx(x0,yy0   ,xlen,twinll,twinhh,0.10, 
				YES, NO, YES, strlen(title[isxtime]), title[isxtime]);
			if(isxtime == YES)
				gleft(x0,yy0-0.30,0.07,timstr,0.0);
		}
		
		if(gsac_control.plotliny){
		if(gsac_control.ygrid == YES)
			dolnygrid(x0,x0+xlen,yy0,dy,depmax,depmin,0.10, 
				YES, gsac_control.ygrid_color, 
				gsac_control.ygrid_type,gsac_control.ygrid_minor);
			doliny(x0+xlen,yy0,dy,depmax,depmin,0.10, YES, NO, NO, 1, " ");
			doliny(x0,yy0,dy,depmax,depmin,0.10, NO, YES, YES, 1, " ");
		} else {
			if(depmin > 0)
				floor = depmin;
			else
				floor = 1.0e-10;
			if(depmax > 0)
				ceiling = depmax;
			else
				ceiling = 1.0e-9;
			if(gsac_control.ygrid == YES){
				dologygrid(x0,x0+xlen,yy0,dy,ceiling,floor,
					0.10,gsac_control.ygrid_color, 
					gsac_control.ygrid_type,
					gsac_control.ygrid_minor);
			}
			dology(x0+xlen,yy0,dy,ceiling,floor,0.10, YES, NO, NO, 1, " ");
			dology(x0,yy0,dy,ceiling,floor,0.10, NO, YES, YES, 1, " ");
		}


		/* end annotate plot titles */
		gbox(x0, yy0, x0+xlen, yy0+dy);
		gclip("on", x0, yy0, x0+xlen, yy0+dy);
		gsac_setcolor(YES, kkk, ntrc);
		gsac_plot_fileid(x0,yyy0,xlen,dy, k);
		txx = xx + (ilow)*uxmul*dx ;
		for(i=ilow ; i < ihgh ; i+= inc ){
			if(gsac_control.plotliny){
			yy = yy0 +dy*(sacdata[k].sac_data[i] - depmin)
				/(depmax - depmin);
			} else {
				/* never take a log of a negative or zero
				 * put in the floor of the depmin from the
				 * ylim */
				if(sacdata[k].sac_data[i] > 0)
					yy = yy0 +dy*log(sacdata[k].sac_data[i] 
					/ floor) /log(ceiling / floor);
				else
					yy = yy0;
			}
			if(i == ilow )
				plot((float)txx,yy,3);
			else
				plot((float)txx,yy,2);
			txx += uxmul*dx*inc;
		}
		/* if decimate, so indicate in trace box */
		if(inc > 1)
			showdec(x0, yy0, x0+xlen, yy0+dy, inc);
		/* show pick */
		if(isxtime == YES)
		gsac_show_pick(x0, y0, xlen, ylen, k, yy0, dy, twin, pabsolute);
		/* show marktime */
		if(markt_on == YES && isxtime == YES)
		gsac_show_markt(x0, y0, xlen, ylen, k, yy0, dy, twin, pabsolute);

		gclip("off", x0, y0, x0+xlen, y0+ylen);
		gsac_setcolor(NO , kkk, ntrc);
		if(gsac_control.plotdevice==WIN)
			XviG_Flush();
	}
		/* annotate with the plot titles
			added 02 OCT 2008 */
/*
fprintf(stderr,"ON %d LOC %d TEXXT %s\n",title_on,title_location,title_text);
fprintf(stderr,"X0 %f XLEN %f Y0 %f YLEN %f\n",x0,xlen,y0,ylen);
*/
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
	gmesg(" ");

}


/* a shorthand way of accessing the SAC header time values 
 * for O A T0 .. T9 */
static int timelist[] = {
  H_O, H_A, H_T0, H_T1, H_T2, H_T3, H_T4, H_T5, H_T6, H_T7, H_T8, H_T9, -1};
static int timechar[] = {
  H_KO, H_KA,  H_KT0, H_KT1, H_KT2, H_KT3, H_KT4, H_KT5, H_KT6, H_KT7, H_KT8, H_KT9, -1 };
void gsac_show_pick(float x0, float y0, float xlen, float ylen, int k, float yy0, float dy, float twin, int pabsolute)
{
	/* search through basic timing values */
	int i;
	float xx, tb, tv, uv, ux;
	float size;
	float uxcen, uxmul;
	uxcen = gsac_control.uxcen;
	uxmul = gsac_control.uxmul;

	for(i=0; timelist[i] >= 0 ; i++){
		tv = sacdata[k].sachdr.rhdr[timelist[i]];
		if(tv != -12345.){
		if(pabsolute == YES) {
			uv = (float)(sacdata[k].tzref + tv - gsac_control.begminx)/twin;
		} else {
			tb = sacdata[k].tzbegx - sacdata[k].tzref ;
			uv = (float)(tv -  tb )/twin;
		}
		/* now get the actual position */
		ux = uxmul*(uv - uxcen) + 0.5;
		xx = x0 + xlen*ux ;
		newpen(1);
		plot(xx,yy0+0.15*dy,3);
		newpen(2);
		plot(xx,yy0+0.8*dy,2);
		plot(xx,yy0+0.8*dy,3);
		newpen(1);
		size=MIN(0.1, 0.1*dy);
		if(timelist[i] == 7)
		gleft(xx,yy0+0.85*dy,size,"O",0.0);
		else
		if(strncmp(sacdata[k].schdr[timechar[i]],"-12345",6)!=0)
		gleft(xx,yy0+0.85*dy,size,sacdata[k].schdr[timechar[i]],0.0);
		}
	}
}

void gsac_show_markt(float x0, float y0, float xlen, float ylen, int k, float yy0, float dy, float twin, int pabsolute)
{
	/* search through basic timing values */
	int i;
	float xx, tb, tv, uv, ux;
	float size;
	float uxcen, uxmul;
	float dist,o;
	char ostr[4];
	uxcen = gsac_control.uxcen;
	uxmul = gsac_control.uxmul;

	for(i=0; i < markt_numvel  ; i++){
		dist = sacdata[k].sachdr.rhdr[H_DIST];
		o    = sacdata[k].sachdr.rhdr[H_O];
		if(dist != -12345.){
		tv = o + dist/markt_vel[i];
		if(pabsolute == YES) {
			uv = (float)(sacdata[k].tzref + tv - gsac_control.begminx)/twin;
		} else {
			tb = sacdata[k].tzbegx - sacdata[k].tzref ;
			uv = (float)(tv -  tb )/twin;
		}
		/* now get the actual position */
		ux = uxmul*(uv - uxcen) + 0.5;
		xx = x0 + xlen*ux ;
		newpen(1);
		plot(xx,yy0+0.3*dy,3);
		newpen(4);
		plot(xx,yy0+0.7*dy,2);
		plot(xx,yy0+0.7*dy,3);
		newpen(1);
		size=MIN(0.07, 0.07*dy);
		if(sprintf(ostr,"%3.1f",markt_vel[i]) < 4){
		gcent(xx,yy0+0.07*dy,size,ostr,0.0);
		}
		}
	}
}
/* put decimation value in an inset in trace window if possible */


void showdec(float xl, float yl, float xh, float yh, int inc){
	/*dedicate lower right corner, and also with a height of only about
		0.2 ylen at most */
	int ndec;
	float ht;
	float xln, yln;
	xln = ABS(xh-xl);
	yln = ABS(yh-yl);
	if(inc > 999)
		ndec = 4;
	else if(inc > 99)
		ndec = 3;
	else if(inc > 9)
		ndec = 2;
	else
		ndec = 1;
	ht = MAX(0.05*yln,0.1);
	gbox(xh-ndec*ht,yl,xh,yl+ht);
	number(xh-ndec*ht+0.2*ht,yl+0.15*ht,0.7*ht,(float)inc,0.0,-1);
}

void dogrid()
{
	float len = 0.2;
	int pat = 21; /* 010101  in binary */
	int i;
	/* draw x- y-axes in blue */
	newpen(4);
	for(i=1 ; i < 10 ; i++){
		plot(i,0.0,3);
		plotd(i,8.0,pat,len);
	}
	for(i=1 ; i < 10 ; i++){
		plot(0.0,i,3);
		plotd(10.0,i,pat,len);
	}
	newpen(1);
}


