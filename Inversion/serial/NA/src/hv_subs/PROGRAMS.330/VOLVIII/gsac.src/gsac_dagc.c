#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

#define AGCWINDOW	1

struct arghdr agcarg[] = {
{AGCWINDOW, "WINDOW" , RHDR, 0, 1, NO, "W window ", 1},
{0      , ""      , IHDR, 0, 0, NO, "",-1}
};

int  agc_do_pct = NO;
float agc_window = 0.0;
float agc_window_pct = 0.0;
static float agc_real[1];

void gsac_set_param_dagc(int ncmd, char **cmdstr)
{
	int i;
	if(ncmd == 1)
		return;
	/* is the command syntax correct ? Also reset */
	if(testarg(ncmd, cmdstr, agcarg, NO, YES))
		return;
	for(i=0 ; agcarg[i].key[0] != '\0' ; i++){
		/* check for special commands */
		if(agcarg[i].used > 0){
			if(agcarg[i].id == AGCWINDOW){
				getargr(ncmd, cmdstr, agcarg[i].key, agcarg[i].mfit, agcarg[i].narg, agc_real);
				agc_window  = agc_real[0];
			}
		}
	}
	if(gsac_control.prs > 0){
		if(gsac_control.prshist == NULL)
			gsac_control.prshist = fopen("prshist.tmp","w+");
		fprintf(gsac_control.prshist,"agc w %g\n",agc_window);
		fflush(gsac_control.prshist);
	}
}

void gsac_doagc(float *x, int n, float dt, float win );

void gsac_exec_dagc(void)
{

	int k, ntrc, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float delta;
	float gate;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
		
	for ( k=0 ; k < ntrc ; k ++){
		npts  = sacdata[k].sachdr.ihdr[H_NPTS];
		delta = sacdata[k].sachdr.rhdr[H_DELTA];

		if(npts > 0){
			if(agc_do_pct == YES){
				gate = agc_window_pct * 0.01 * npts * delta;
			} else {
				gate = agc_window;
			}
			gsac_doagc(sacdata[k].sac_data, npts, delta,  gate );
			getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		}
	}
}

void gsac_doagc(float *x,  int n, float dt, float win )
{
	/* x - array from which to form the AGC 
	 * g - gain factor from trace
	 * n - number of points in x
	 * dt - sampling interval of trace
	 * win - smoothing window for agc
	 * */
	/*
	 * m - length of smoothing window, odd number
	 * g - gain factor from trace
	 * */
	int m, m2, i;
	float sum;
	float S1;
	float lower_bound;
	float max;
	float *g;
	/* safety */
	if(n < 2)
		return;
	/* create the gain array */
	if((g = (float *)calloc(n, sizeof(float))) == (float *)NULL){
		printf("Cannot allocate memory to create gain array in agc\n");
		return;
	}
	/* define the smoothing width in samples taking care to
	 * address extreme values. Also make the width an odd number
	 * */
	m = win/dt;
	if( m > n)
		m = n/2;
	if(m < 3 )
		m = 3;
	m += (m%2 -1 );
	m2 = m / 2;
	/* we define the operator as a running mean over the trace
	 * for efficiency we note that for the next value we just have to add
	 * a new endpoint and remove the first point
	 *
	 * For simplicity the first m points are set to the same value as are
	 * the last m points
	 * */
	for ( i= 0, S1=0.0 ; i < m ; i++)
		S1 += ABS(x[i]);
	/* get gain in the first section */
	for(i = 0 ; i <= m2  ; i++){
		g[i] =  S1 / m ;
	}
	/* now get everything in the middle */
	sum = S1;
	for(i = m2 +1  ; i < n - m2  ; i++ ){
		sum += ABS( x[i+m2]) - ABS( x[i-m2-1]);
		g[i] = sum / m;
	}
	/* now fill in the tail */
	for(i = n -m2 ; i < n ; i++)
		g[i] = sum / m;
	/* as of now the g[i] represents an envelope of the trace
	 * We will now eliminate any zeros, and then invert it to
	 * get the gain factor
	 * */
	max = 0.0;
	for(i=0 ; i < n ; i++){
		if(g[i] > max) max = g[i];
	}
	/* safety */
	if(max == 0.0){
		lower_bound = 1.0 ;
	} else {
		lower_bound = 0.001 * max;
	}
	for(i = 0 ; i < n ; i++)
			g[i] = 1.0/ MAX(g[i], lower_bound);
	/* now adjust the amplitude */
	for(i = 0 ; i < n ; i++)
			x[i] *= g[i];
	/* clean up */
	free(g);
	return;
}
