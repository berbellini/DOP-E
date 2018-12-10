#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	TAPER_DFLT	0
#define	TAPER_COS	1
#define	TAPER_HAM	2
#define	TAPER_HAN	3
#define	TAPER_WID	4


struct arghdr taperarg[] = {
	{TAPER_DFLT, "DEFAULT", IHDR, 0, 0, NO, "", -1},
	{TAPER_COS , "COSINE" , IHDR, 0, 0, NO, "COSINE", 1},
	{TAPER_HAM , "HAMMING", IHDR, 0, 0, NO, "HAMMING",  2},
	{TAPER_HAM , "HANNING", IHDR, 0, 0, NO, "HANNING",  2},
	{TAPER_WID , "WIDTH",   RHDR, 0, 1, NO, "Width w ", 1},
	{0,	""		, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float taper_real[10];
int   taper_int [10];
int   taper_yn;
int   taper_num;

static int taper_type = TAPER_HAN ;
static float taper_w    = 0.05;
static float gsac_taper(int i, int npts, int nt);
static float gsac_taperfunc(float x);

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_taper(int ncmd, char **cmdstr)
{
	int i;

	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, taperarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; taperarg[i].key[0] != '\0' ; i++){
		if(taperarg[i].used > 0){
			if(taperarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, taperarg[i].key, 
					taperarg[i].mfit,taperarg[i].narg, taper_real);
			} else if(taperarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, taperarg[i].key, 
					taperarg[i].mfit,taperarg[i].narg, taper_int );
			} else if(taperarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, taperarg[i].key, 
					taperarg[i].mfit,taperarg[i].narg, &taper_yn );
			} else if(taperarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, taperarg[i].key, 
					taperarg[i].mfit,taperarg[i].narg, &taper_num );
			}
			switch(taperarg[i].id){
				case TAPER_COS:
					taper_type = TAPER_COS;
					break;
				case TAPER_HAM:
					taper_type = TAPER_HAM;
					break;
				case TAPER_HAN:
					taper_type = TAPER_HAN;
					break;
				case TAPER_WID:
					taper_w= taper_real[0];
					if(taper_w > 0.5)
						taper_w = 0.5;
					if(taper_w < 0.0)
						taper_w = 0.0;
					break;

			}
		}
	}
			
		
}

void gsac_exec_taper(void)
{
	/* to the tpaer here */
	int i, k, ntrc, npts, nt;
	float depmax, depmin, depmen;
	int indmax, indmin;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
	if(taper_w <= 0.0)
		return;
		
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		nt = (int)(taper_w * npts);
		if(nt > 0){
			/* perform the taper - since the npts can vary
				cannot save time by precomputing */
			for(i=0; i < npts ; i++){
				sacdata[k].sac_data[i] *= gsac_taper(i,npts,nt);
			}
			getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		}
	}
}

static float gsac_taper(int i, int npts, int nt){
	int nl, nh;
	nl = nt;
	nh = npts - nt;
	if(i < nl)
		return( gsac_taperfunc((float)i/(float)nt) );
	else if(i > nh)
		return( gsac_taperfunc((float)(npts-1-i)/(float)nt) );
	else
		return (1.0);
}
static float gsac_taperfunc(float x){
	if(taper_type == TAPER_COS)
		return(sin(3.1415927*x/2.0));
	else if(taper_type == TAPER_HAN)
		return(0.5 - 0.5*cos(3.1415927*x));
	else if(taper_type == TAPER_HAM)
		return(0.54 - 0.46*cos(3.1415927*x));
	else
		return(sin(3.1415927*x/2.0));	/* this is the default */
		
}
