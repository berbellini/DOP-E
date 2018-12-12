#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

/* Changes:
	14 APR 2007 - ensure that if HALFWIDTH < DT, nothing is done */

extern struct sacfile_ *sacdata;

#define	TRAPEZOID_DFLT	0
#define	TRAPEZOID_WIDTH		1


struct arghdr trapezoidarg[] = {
	{TRAPEZOID_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{TRAPEZOID_WIDTH  , "WIDTH"  , RHDR, NO, 3, NO, "Width t1 t2 t3 x0", 1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
static float trapezoid_real[10];
static int   trapezoid_int [10];
static int   trapezoid_yn;
static int   trapezoid_num;
static float   trapezoid_t1;
static float   trapezoid_t2;
static float   trapezoid_t3;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_trapezoid(int ncmd, char **cmdstr)
{
	int i;
	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, trapezoidarg, NO, YES))
		return;
	/* parse commands */
	trapezoid_t1 = -1.0 ;
	trapezoid_t2 = -1.0 ;
	trapezoid_t3 = -1.0 ;
	for(i=0 ; trapezoidarg[i].key[0] != '\0' ; i++){
		if(trapezoidarg[i].used > 0){
			if(trapezoidarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, trapezoidarg[i].key, 
					trapezoidarg[i].mfit,trapezoidarg[i].narg, trapezoid_real);
			} else if(trapezoidarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, trapezoidarg[i].key, 
					trapezoidarg[i].mfit,trapezoidarg[i].narg, trapezoid_int );
			} else if(trapezoidarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, trapezoidarg[i].key, 
					trapezoidarg[i].mfit,trapezoidarg[i].narg, &trapezoid_yn );
			} else if(trapezoidarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, trapezoidarg[i].key, 
					trapezoidarg[i].mfit,trapezoidarg[i].narg, &trapezoid_num );
			}
			switch(trapezoidarg[i].id){
				case TRAPEZOID_WIDTH:
					trapezoid_t1 = trapezoid_real[0];;
					trapezoid_t2 = trapezoid_real[1];;
					trapezoid_t3 = trapezoid_real[2];;
					break;

			}
		}
	}
			
		
}

void gsac_exec_trapezoid(void)
{
	/* test for bad arguments */
	if(trapezoid_t1 >= 0.0 && trapezoid_t2 >= 0.0 && trapezoid_t3 >= 0)
		if(trapezoid_t1 > 0.0 || trapezoid_t2 > 0.0 || trapezoid_t3 > 0)
			gsac_pulconv(trapezoid_t1, trapezoid_t2, trapezoid_t3);
}

static double *f = (double *)NULL;
static double *t = (double *)NULL;

void gsac_pulconv(float t1, float t2, float t3)
{
	float ht;
	int i, j, k, m, ntrc, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float dt, tmp;
	float time, time1, time2, time3;
	int np;
	float p;
	float areaadj;

	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
		
	ht = 2.0/(t1 + t2 + t2 + t3);
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		dt = sacdata[k].sachdr.rhdr[H_DELTA];
		if(npts > 0 && dt > 0.0 && (t1 + t2 + t2 + t3) > dt){
		/* create a temporary array for the initial file
			and for the filter */
			np = 1 + (int)(t1 + t2 + t3)/dt ;
			if(t == (double *)NULL)
				t = (double *)calloc(npts,sizeof(double));
			else
				t = (double *)realloc(t,npts*sizeof(double));
			if(f == (double *)NULL)
				f = (double *)calloc(np,sizeof(double));
			else
				f = (double *)realloc(f,np*sizeof(double));
			for(i=0; i < npts ; i++)
				t[i] = sacdata[k].sac_data[i] ;
			time = 0.0 ;
			time1 = t1;
			time2 = t1 + t2;
			time3 = t1 + t2 + t3;
			areaadj = 0.0 ;
			for(j=0; j < np ; j++){
				p = 0.0;
				if(time <= time1 && time1 > 0.0) {
					p = (time - 0.0)/(time1 - 0.0);
				} else if(time > time1 && time <= time2 ) {
					p = 1.0;
				} else if(time > time2 && time <= time3 && time3 > time2 ) {
					p = 1.- (time - time2)/(time3 - time2);
				} else {
					/* should never get here */
					p = 0.0;
				}
				f[j] = p*ht;
				areaadj += f[j] ;
				time +=  dt;
			}
			areaadj *=dt ;
			if(areaadj > 0.0){
			/* now convolve */
			for(i=0;i < npts ; i++){
				tmp = 0.0 ;
				for (j=0 ; j < np ; j++){
					m = i - j ;	
					if( m >= 0 && m < npts)
					tmp += f[j]*sacdata[k].sac_data[m];
				}
				t[i] = tmp*dt/areaadj;
			}

			for(i=0; i < npts ; i++)
				sacdata[k].sac_data[i] =t[i];
			getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
			}
	
		}
	}
}

