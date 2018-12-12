#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	RICKER_DFLT	0
#define	RICKER_F	1
#define PI 3.1415927


struct arghdr rickerarg[] = {
	{RICKER_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", 1},
	{RICKER_F  , "F"  , RHDR, NO, 1, NO, "F frequency", 1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float ricker_real[10];
static float ricker_freq = 25.0;
/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_ricker(int ncmd, char **cmdstr)
{
	int i;
	/* default */

	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, rickerarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; rickerarg[i].key[0] != '\0' ; i++){
		if(rickerarg[i].used > 0){
			if(rickerarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, rickerarg[i].key, 
					rickerarg[i].mfit,rickerarg[i].narg, ricker_real);
			}
			switch(rickerarg[i].id){
				case RICKER_F:
					ricker_freq = ricker_real[0];
					break;
				case RICKER_DFLT:
					ricker_freq = 25.0;
					break;

			}
		}
	}
			
		
}

static double *f = (double *)NULL;
static double *t = (double *)NULL;

void gsac_exec_ricker(void)
{
	fprintf(stderr,"Ricker frequency %f\n",ricker_freq);
	int i, j, k, m, ntrc, npts;
	float time, toffset;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float dt, tmp;
	int np;
	float p;
	int w;

	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
	/* at time = 0 the Ricker = 1 */
		
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		dt = sacdata[k].sachdr.rhdr[H_DELTA];
		if(npts > 0 && dt > 0.0 ){
		/* create a temporary array for the initial file
			and for the filter */
			w = 1.0/(ricker_freq*dt);
			np = 2*w + 1;
			if(t == (double *)NULL)
				t = (double *)calloc(2*npts,sizeof(double));
			else
				t = (double *)realloc(t,2*npts*sizeof(double));
			if(f == (double *)NULL)
				f = (double *)calloc(np,sizeof(double));
			else
				f = (double *)realloc(f,np*sizeof(double));
			for(i=0; i < npts ; i++)
				t[i] = sacdata[k].sac_data[i] ;
			toffset = 0.0;
			for(j=0; j < np ; j++){
				time = (-w + j-1)*dt;
				tmp = PI*ricker_freq*time ;
				tmp = (1.0-2.0*PI*PI*ricker_freq*ricker_freq*time*time)*exp(-tmp*tmp);
					f[j] = tmp;

			}
			/* now convolve */
			for(i=0;i < npts + w ; i++){
				tmp = 0.0 ;
				for (j=0 ; j < np ; j++){
					m = i - j ;	
					if( m >= 0 && m < npts)
					tmp += f[j]*sacdata[k].sac_data[m];
				}
				if( i >= w)
					t[i-w] = tmp*dt;
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
