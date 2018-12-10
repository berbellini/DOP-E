/* NOTES
 * before ploting traces in a grounp determine the extremes of the absolute
 * B and E so that the time scale can be defined
 */
/* CHANGES
 09 JAN 2005 - gsac_plotsp - yy0 was not defined before first call to
	glip - actualy makes not difference but initialized baker@usgs.gov
*/
#include	<stdio.h>
#include	"gsac.h"
#include	"gsac_plotsp.h"
#include	"gsac_sac.h"
#include	"gsac_arg.h"

extern struct sacfile_ *sacdata;
extern int  sfgetline(FILE *fp, char s[], int lim);



#define PSP_AM 0
#define PSP_PH 1
#define PSP_OVERLAY 2
#define PSP_PERPLOT 3
#define PSP_XLOG 4
#define PSP_XLIN 5
#define PSP_YLOG 6
#define PSP_YLIN 7
#define PSP_FMIN 8
#define PSP_FMAX 9
#define PSP_DEFAULT 10
#define PSP_AMIN 11
#define PSP_AMAX 12

#define PLOT_AM 0
#define PLOT_PH 1

static int psp_doamph = PLOT_AM;

struct arghdr psparg[] = {
	{PSP_AM, "AMPLITUDE"	, IHDR, 0, 0, NO, "AMPLITUDE ",3},
	{PSP_AM, "AM"	, IHDR, 0, 0, NO, "AMPLITUDE ",-1},
	{PSP_PH, "PHASE"	, IHDR, 0, 0, NO, "PHASE", 2},
	{PSP_OVERLAY, "OVERLAY" , YHDR, 0, 1, NO, "OVERLAY [ON|OFF] ", 1},
	{PSP_PERPLOT, "PERPLOT" , NHDR, 0, 1, NO, "PERLOT [n|OFF]", 2},
	{PSP_XLIN, "XLIN", IHDR, 0, 0, NO, "XLIN", 3},
	{PSP_XLOG, "XLOG", IHDR, 0, 0, NO, "XLOG", 3},
	{PSP_YLIN, "YLIN", IHDR, 0, 0, NO, "YLIN", 3},
	{PSP_YLOG, "YLOG", IHDR, 0, 0, NO, "YLOG", 3},
	{PSP_FMIN, "FMIN", RHDR, 0, 1, NO, "FMIN", 3},
	{PSP_FMAX, "FMAX", RHDR, 0, 1, NO, "FMAX", 3},
	{PSP_AMIN, "AMIN", RHDR, 0, 1, NO, "AMIN", 3},
	{PSP_AMAX, "AMAX", RHDR, 0, 1, NO, "AMAX", 3},
	{PSP_DEFAULT, "DEFAULT", RHDR, 0, 0, NO, "DEFAULT", 1},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};

static int pspperplot = -1;
static int psp_overlay = NO;
static int psp_yn;
static int psp_num;
static int psp_ylin = NO;
static int psp_xlin = NO;
static float psp_fmin = -1;
static float psp_fmax = 1.0e+38;
static float psp_amin = 0.0;
static float psp_amax = 1.0e+38;

static float psp_real[10];

float *y = (float *)NULL;
extern int *sortptr;

extern void XviG_Flush();
extern void dogrid(void);

extern void gsac_plot_fileid(float x0,float y0,float xlen,float ylen, int k);
extern void gsac_exec_hold(void);

void gsac_plotsp(float x0,float y0,float xlen,float ylen,float fe,float fs,int ns,int ne,int ntrc, float depmax, float depmin, int numpspframe, float dy);

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


void gsac_set_param_plotsp(int ncmd, char **cmdstr)
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

	float tmpmx, tmpmn;

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
			if(Color >= 4)
				gsac_control.black = 0;
			else
				gsac_control.black = 1;
			gsac_control.kolor = Color%4;


		}
	}
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, psparg, YES, YES))
	       	return	;
	/* now determine which was set */
	for(i=0 ; psparg[i].key[0] != '\0' ; i++){
		if(psparg[i].used > 0){	
			if(psparg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, psparg[i].key, 
					psparg[i].mfit,psparg[i].narg, &psp_yn );
			} else if(psparg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, psparg[i].key, 
					psparg[i].mfit,psparg[i].narg, &psp_num );
			} else if(psparg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, psparg[i].key, 
					psparg[i].mfit,psparg[i].narg, psp_real );
			}
			psparg[i].used = 0;
			switch(psparg[i].id){
				case PSP_PERPLOT:
					pspperplot = psp_num;
					break;
				case PSP_AM:
					psp_doamph = PLOT_AM;
					break;
				case PSP_PH:
					psp_doamph = PLOT_PH;
					break;
				case PSP_OVERLAY:
					if(psp_yn == NO)
						psp_overlay = NO;
					else if(psp_yn == YES)
						psp_overlay = YES;
					break;
				case PSP_XLOG:
					psp_xlin = NO;
					break;
				case PSP_XLIN:
					psp_xlin = YES;
					break;
				case PSP_YLOG:
					psp_ylin = NO;
					break;
				case PSP_YLIN:
					psp_ylin = YES;
					break;
				case PSP_FMIN:
					psp_fmin = psp_real[0];;
					break;
				case PSP_FMAX:
					psp_fmax = psp_real[0];;
					break;
				case PSP_AMIN:
					psp_amin = psp_real[0];;
					break;
				case PSP_AMAX:
					psp_amax = psp_real[0];;
					break;
				case PSP_DEFAULT:
					pspperplot = -1;
					psp_overlay = NO;
					psp_ylin = NO;
					psp_xlin = NO;
					psp_fmin = -1;
					psp_fmax = 1.0e+38;
					psp_amin = 0.0;
					psp_amax = 1.0e+38;
					break;
			}
		}
	}
	/* safety check on fmax */
	tmpmx = MAX(psp_fmax, psp_fmin);
	tmpmn = MIN(psp_fmax, psp_fmin);
	psp_fmax =tmpmx;
	psp_fmin =tmpmn;
	tmpmx = MAX(psp_amax, psp_amin);
	tmpmn = MIN(psp_amax, psp_amin);
	psp_amax =tmpmx;
	psp_amin =tmpmn;

}


void gsac_exec_plotsp(void)
{
int i,j,k,kkk;
int ntrc;
float x0, y0, xlen, ylen, dy;
float depmax, depmin;
int indmax, indmin;
char pltname[10];
char instr[100];
int numpspframe;
int n2;
float temp;
float fs, fe;
float tr, ti;


	/* if there are no traces return */
	ntrc = gsac_control.number_otraces;
	gsac_exec_hold();
	if(ntrc < 1)
		return;

	if(gsac_control.fft == NO){
		printf("Execute FFT first before plot\n");
		return;
	}
	/* initialize */
	if(gsac_control.plotdevice==WIN){
		if(gsac_control.hold == NO){
			gframe(2);
		} else if(gsac_control.hold == 1){
			gframe(2);
			gsac_control.hold++ ;
		}
	}

	if(pspperplot > 0)
		if(pspperplot > ntrc)
			numpspframe = ntrc;
		else
			numpspframe = pspperplot;
	else
		numpspframe = ntrc;

	xlen = gsac_control.xlen ;
	ylen = gsac_control.ylen ;
	x0   = gsac_control.x0 ;
	y0   = gsac_control.y0 ;

	if(psp_overlay == YES){
		dy = ylen;
	} else {
		dy = ylen / numpspframe;
	}

	/* temporary */
	/* since all are plotted to the same frequency scale we must
	 * define the lower and upper limits. As a kludge, the minimum
	 * frequency plotted with log x-axis is the DF of first trace 
	 * in memory */
	/* get global fmax */
	fe = MIN(gsac_control.fmax, psp_fmax);
	if(psp_xlin == YES){
		fs = MAX (psp_fmin, 0.0);
	} else {
		fs = MAX (psp_fmin, sacdata[0].df);
	}


	/* get extreme amplitudes for the case of an overlay */
	if(psp_overlay == YES ){
		/* never overlay phase spectra */
		if(psp_doamph == PLOT_AM){
			depmax = 0.0;
			for ( kkk=0 ; kkk < ntrc ; kkk++){
				k = sortptr[kkk];
				n2 = sacdata[k].npow2/2.0  ;
				for(i=0 , j= 0; i <= n2 ; i++){
					tr = sacdata[k].sac_spectra[j++];
					ti = sacdata[k].sac_spectra[j++];
					temp = sqrt(tr*tr + ti*ti);
					if(temp > depmax)
						depmax = temp;
				}
			}
			/* adjust amplitude for AMIN AMAX */
			if(psp_ylin == YES){
				depmax = MIN(depmax * 1.2, psp_amax);
				depmin = psp_amin ;
			} else {
				if(psp_amax < 1.0e+37)
					depmax = psp_amax;
				else
					depmax *= 1.5;
				if(psp_amin > 0.0)
					depmin = psp_amin;
				else
					depmin = depmax/10000.0 ;
			}
		}
	}

	for ( kkk=0 ; kkk < ntrc ; kkk+=numpspframe){
		if(gsac_control.plotdevice==PLT){
			if(gsac_control.hold == NO || 
				(gsac_control.hold != NO && 
				 gsac_control.inpltmode == NO)){
			sprintf(pltname,"P%3.3d.PLT",gsac_control.plotcount++);
			printf("Initializing %s\n",pltname);
			ginitf(pltname,"GSAC");
			gsac_control.inpltmode = YES ;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			if(gsac_control.grid)
				dogrid();
			}
		} else {
			if(gsac_control.hold == NO){
				gframe(2);
				gclip("off", x0, y0, x0+xlen, y0+ylen);
				if(gsac_control.grid)
					dogrid();
			}
		}
		newpen(1);
		gsac_plotsp(x0,y0,xlen,ylen,fe,fs,kkk,MIN(kkk+numpspframe,ntrc),
				ntrc,depmax,depmin,numpspframe, dy);
		if(gsac_control.plotdevice!=PLT){
			if( pspperplot > 0){
				printf("More? y/n\n");
				sfgetline(stdin, instr, 100);
				if(instr[0] == 'n' || instr[0] == 'N')
					goto jump;
			}
		}
		if(gsac_control.plotdevice==PLT && gsac_control.hold == NO){
			/* force new Pnnnn.PLT on next call */
			gsac_control.plotchange = NO; 
			gend(0);
			gsac_control.inpltmode = NO; 
		}
	}

	gclip("off", x0, y0, x0+xlen, y0+dy);
jump:
	/* clean up */
	if(gsac_control.hold == NO){
		if(gsac_control.everinteractive == YES){
			ginitf("INTEM","GSAC");
		}
	}

}

void gsac_plotsp(float x0,float y0,float xlen,float ylen,float fe,float fs,int ns,int ne,int ntrc, float depmax, float depmin, int numpspframe, float dy)
{
int kk, kkk, n2;
int i, j, k, is;
float tr, ti;
float depmen;
	int indmax, indmin;
float df;
float yy0, yyy0;
float xx, yy;
float ht;
float freq;
	/* do the bottom axis */
	ht = 0.10;
	yy0 = 1.0;
	gclip("off", x0, yy0, x0+xlen, yy0+dy);
	/* set the background */
	if(gsac_control.background == YES &&
		gsac_control.background_color >= 0){
		newpen(gsac_control.background_color);
		if(psp_overlay == YES){
		shader(x0,y0,x0+xlen,y0+ylen,0,0,0.01,0.01);
		} else {
		shader(x0,y0+ylen-(ne-ns)*dy,x0+xlen,y0+ylen,0,0,0.01,0.01);
		}
		newpen(1);
	}
	for(kkk=ns;kkk < ne;kkk++){
		k = sortptr[kkk];
		kk = kkk%numpspframe;

		n2 = sacdata[k].npow2/2.0  ;
		y = (float *)realloc(y,(n2+1)*sizeof(float));
		/* now fill array with amplitude or phase spectra */
		for(i=0 , j= 0; i <= n2 ; i++){
			tr = sacdata[k].sac_spectra[j++];
			ti = sacdata[k].sac_spectra[j++];
			if(psp_doamph == PLOT_AM){
				y[i] = sqrt(tr*tr + ti*ti);
			} else if (psp_doamph == PLOT_PH){
				y[i] = atan2(ti,tr);
			}
		}
		if(psp_doamph == PLOT_AM && psp_overlay == NO){
			getmxmn(y, n2+1,&depmax, &depmin, &depmen,&indmax,&indmin);
			/* adjust amplitude for AMIN AMAX */
			if(psp_ylin == YES){
				depmax = MIN(depmax * 1.1, psp_amax);
				depmin = psp_amin ;
			} else {
				if(psp_amax < 1.0e+37)
					depmax = psp_amax;
				else
					depmax *= 1.5;
				if(psp_amin > 0.0)
					depmin = psp_amin;
				else
					depmin = depmax/10000.0 ;
			}
		} else if (psp_doamph == PLOT_PH){
			depmax =  3.1415927;
			depmin = -3.1415927;
			depmen = 0.0;
		}



		df  = sacdata[k].df;

		/* use depmax depmin */
		if(psp_overlay == YES){
			yy0 = y0;
			yyy0 = yy0 -kk*1.5*ht ;
		} else {
			yy0 = y0 + (numpspframe -1  - kk )*dy;
			yyy0 = yy0;
		}
		if(psp_xlin == YES){
			if(gsac_control.xgrid == YES)
				dolnxgrid(x0,yy0,yy0+dy,xlen,fe,fs,0.10, YES, 
					gsac_control.xgrid_color, 
					gsac_control.xgrid_type,
					gsac_control.xgrid_minor);
			dolinx(x0,yy0+dy,xlen,fe,fs,0.10, NO, NO, NO, 0, " ");
		} else {
			if(gsac_control.xgrid == YES){
				dologxgrid(x0,yy0,yy0+dy,xlen,fe,fs,0.10,
					gsac_control.xgrid_color,
					gsac_control.xgrid_type,
					gsac_control.xgrid_minor);
			}
			dologx(x0,yy0+dy,xlen,fe,fs,0.10, NO, NO, NO, 0, " ");
		}
		gbox(x0, yy0, x0+xlen, yy0+dy);
		if(kkk == ne -1 ){
			if(psp_xlin == YES) {
				dolinx(x0,yy0,xlen,fe,fs,0.10, YES, NO, YES, 14, "Frequency (Hz)");
			} else {
				dologx(x0,yy0,xlen,fe,fs,0.10,YES,NO,YES,14,"Frequency (Hz)");
			}
		}
		if(psp_doamph == PLOT_AM){
			if(psp_ylin == YES){
				if(gsac_control.ygrid == YES)
					dolnygrid(x0,x0+xlen,yy0,dy,depmax,
						depmin, 0.10, YES, 
						gsac_control.ygrid_color, 
						gsac_control.ygrid_type,
						gsac_control.ygrid_minor);
				doliny(x0,yy0,dy,depmax,depmin,0.10, NO, YES, YES, 1, " ");
				doliny(x0+xlen,yy0,dy,depmax,depmin,0.10, YES, NO, NO, 0, " ");
			} else {
				if(gsac_control.ygrid == YES){
					dologygrid(x0,x0+xlen,yy0,dy,depmax,
						depmin,0.10,
						gsac_control.ygrid_color,
						gsac_control.ygrid_type,
						gsac_control.ygrid_minor);
				}
				dology(x0,yy0,dy,depmax,depmin,0.10,NO,YES,YES,1," ");
				dology(x0+xlen,yy0,dy,depmax,depmin,0.10,YES,NO,NO,0," ");
			}
		} else {
			if(gsac_control.ygrid == YES){
				dolnygrid(x0,x0+xlen,yy0,dy,depmax,depmin,0.10, 
					YES, gsac_control.ygrid_color, 
					gsac_control.ygrid_type,
					gsac_control.ygrid_minor);
			}
			doliny(x0,yy0,dy,depmax,depmin,0.10, NO, YES, YES, 1, " ");
			doliny(x0+xlen,yy0,dy,depmax,depmin,0.10, YES, NO, NO, 0, " ");
		}
		gclip("on", x0, yy0, x0+xlen, yy0+dy);
		gsac_setcolor(YES, kkk, ntrc);
		gsac_plot_fileid(x0,yyy0,xlen,dy, k);
		/* plot the spectra */
		for(i=0, is=0 ; i <= n2 ; i++){
				freq = i * df;
				if(freq >= fs && freq <= fe){
				if(psp_doamph == PLOT_AM){
					if(psp_ylin == YES){
				yy = yy0 +dy*(y[i] - depmin)/(depmax - depmin);
					} else {
				yy = yy0 + dy*log10(y[i] / depmin)/log10(depmax / depmin);
					}
				} else {
				yy = yy0 +dy*(y[i] - depmin)/(depmax - depmin);
				}
				if(psp_xlin == YES){
				xx = x0 + xlen*(freq-fs)/(fe-fs);
				} else {
				xx = x0 + xlen*log10(freq/fs)/log10(fe/fs);
				}
	
		
				if(is == 0 ){
					plot(xx,yy,3);
					is = 1;
				} else {
					plot(xx,yy,2);
				}
			}
		}
		gsac_setcolor(NO, kkk, ntrc);
		gclip("off", x0, yy0, x0+xlen, yy0+dy);
	}
		/* annotate with the plot titles
			added 29 MAY 2009 */
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
