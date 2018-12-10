#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	SMTH_DFLT	0
#define	SMTH_MEAN	1
#define	SMTH_MEDIAN	2
#define	SMTH_HALF	3
#define	SMTH_PASS	4

#define SMTH_HALF_MAX 128


struct arghdr smtharg[] = {
	{SMTH_DFLT, "DEFAULT", IHDR, 0, 0, NO, "", 1},
	{SMTH_MEAN  , "MEAN"  , IHDR, 0, 0, NO, "", 3},
	{SMTH_MEDIAN, "MEDIAN"  , IHDR, 0, 0, NO, "", 3},
	{SMTH_HALF, "HALFWIDTH", IHDR, 0, 1, NO, "HALFWIDTH n ", 1},
	{SMTH_PASS, "PASS", IHDR, 0, 1, NO, "PASS p ", 1},
	{0,	""		, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
static int   smth_int [10];

/* control parameters */
static int smth_type = SMTH_MEAN;
static int smth_half = 1;
static int smth_pass = 1;
static float *tx = (float *)NULL;
void smth_mean(float *x, int npts, int pass, int half);
void smth_median(float *x, int npts, int pass, int half);
static float smth_med_ar[SMTH_HALF_MAX];
static int smth_med_key[SMTH_HALF_MAX];


/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_smth(int ncmd, char **cmdstr)
{
	int i;
	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, smtharg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; smtharg[i].key[0] != '\0' ; i++){
		if(smtharg[i].used > 0){
			if(smtharg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, smtharg[i].key, 
					smtharg[i].mfit,smtharg[i].narg, smth_int );
			}
			switch(smtharg[i].id){
				case SMTH_MEAN:
					smth_type = SMTH_MEAN;
					break;
				case SMTH_MEDIAN:
					smth_type = SMTH_MEDIAN;
					break;
				case SMTH_HALF:
					if(smth_int[0] > SMTH_HALF_MAX)
						smth_int[0] = SMTH_HALF_MAX;
					if(smth_int[0] > 0)
						smth_half = smth_int[0];
					else
						printf("SMOOTH HALF n must have n > 0\n");
					break;
				case SMTH_PASS:
					if(smth_int[0] >= 0)
						smth_pass = smth_int[0];
					else
						printf("SMOOTH PASS n must have n >= 0\n");
					break;
				case SMTH_DFLT:
					smth_half = 1;
					smth_type = SMTH_MEAN;
					smth_pass = 1;
			}
		}
	}
			
		
}

void gsac_exec_smth(void)
{

	int k, ntrc, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	/* if there are no traces return */
printf("smooth type %d pass %d half %d\n",smth_type, smth_pass, smth_half);
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
		
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];

		if(npts > 0){
			if(smth_type == SMTH_MEAN)
				smth_mean(sacdata[k].sac_data, npts, smth_pass, smth_half);
			else if(smth_type == SMTH_MEDIAN)
				smth_median(sacdata[k].sac_data, npts, smth_pass, smth_half);
			getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		}
	}
}

void smth_mean(float *x, int npts, int pass, int half)
{
	int np, nwid;
	int i;
	float sum;
	nwid = 2*half + 1;
	/* safety */
	if(nwid > npts)
		return;
	/* set up a temporary array */
	if(tx == (float *)NULL)
		tx = (float *)calloc(npts,sizeof(float));
	else
		tx = (float *)realloc(tx,npts*sizeof(float));
	for(np = 0 ; np < pass ; np++){
		/* start running average */
		for(i=0, sum=0.0; i < nwid ; i++)
			sum += x[i];
		/* get running average - later fill in first and last points */
		for(i=half+1;i < npts-half; i++){
			sum += x[i+half] - x[i-half-1];
			tx[i] = sum/nwid;
		}
		/* get the beginning */
		for(i=0 ; i < half+1; i++)
			tx[i] = tx[half];
		/* get the end */
		for(i=npts-half ; i < npts; i++)
			tx[i] = tx[npts-half-1];
		/* restore the array */
		for(i=0;i < npts; i++)
			x[i] = tx[i];
	}
}

void gnomesort(int n, float ar[], int key[] ) ;
void smth_median(float *x, int npts, int pass, int half)
{
	int np, nwid;
	int i, j, k;
	nwid = 2*half + 1;
	/* safety */
	if(nwid > npts)
		return;
	/* set up a temporary array */
	if(tx == (float *)NULL)
		tx = (float *)calloc(npts,sizeof(float));
	else
		tx = (float *)realloc(tx,npts*sizeof(float));
	for(np = 0 ; np < pass ; np++){
		/* get running average - later fill in first and last points */
		for(i=half;i < npts-half; i++){
			/* perhaps this can be made smarter to eliminate
			 * 	this for loop
			 */
			for(j=i-half,k=0;j < i+half+1;j++, k++){
				smth_med_ar[k]=x[j];
				smth_med_key[k]=k;
			}
			gnomesort(nwid, smth_med_ar, smth_med_key ) ;
			tx[i] = smth_med_ar[half];
		}
		/* get the beginning */
		for(i=0 ; i < half; i++)
			tx[i] = tx[half];
		/* get the end */
		for(i=npts-half ; i < npts; i++)
			tx[i] = tx[npts-half-1];
		/* restore the array */
		for(i=0;i < npts; i++)
			x[i] = tx[i];
	}
}

