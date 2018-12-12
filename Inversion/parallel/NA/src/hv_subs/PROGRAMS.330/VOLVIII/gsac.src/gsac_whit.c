/* Created 10 JAN 2005
   purpose: to whiten spectral according to ground noise model
   BAsically we want to reduce the effect fo the miroseisms
   Eventually we will use the signal spectrum to define the
   whitening filter
*/

#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

static float *y = (float *)NULL;
static float *amp = (float *)NULL;

#define	WHIT_DFLT	0
#define WHIT_FREQ	1
#define WHIT_ABS	2


void gsac_dowhit(float **x, float dt, int npts);
float gsac_taper3(float xl, float xh, float x);

struct arghdr whitarg[] = {
	{WHIT_DFLT, "DEFAULT", IHDR, 0, 0, NO, "", 1},
	{WHIT_FREQ, "FREQLIMITS" , RHDR, 0, 4,YES, "FREQLIMITS F1 F2 F3 F4",4},
	{WHIT_ABS , "ABSOLUTE",IHDR, 0, 0, NO, "ABSOLUTE", 1},
	{0,	""  	     , IHDR, 0, 0, NO, "", -1}
};

/* frequencies for freqlimit control important for deconvolution */
static float whit_f1 = -10.;
static float whit_f2 =  -5.;
static float whit_f3 =  1.0e6;
static float whit_f4 =  1.0e7;
static int whit_dotaper = NO;
static int whit_doabs = NO;

/* these are temporary variables only used here */
float whit_real[10];
int   whit_int [10];
int   whit_yn;
int   whit_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_whit(int ncmd, char **cmdstr)
{
	int i;
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, whitarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; whitarg[i].key[0] != '\0' ; i++){
		if(whitarg[i].used > 0){
			if(whitarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, whitarg[i].key, 
					whitarg[i].mfit,whitarg[i].narg, whit_real);
			} else if(whitarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, whitarg[i].key, 
					whitarg[i].mfit,whitarg[i].narg, whit_int );
			} else if(whitarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, whitarg[i].key, 
					whitarg[i].mfit,whitarg[i].narg, &whit_yn );
			} else if(whitarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, whitarg[i].key, 
					whitarg[i].mfit,whitarg[i].narg, &whit_num );
			}
			switch(whitarg[i].id){
				case WHIT_DFLT:
					whit_f1 = -10.;
					whit_f2 =  -5.;
					whit_f3 =  1.0e6;
					whit_f4 =  1.0e7;
					whit_dotaper = NO;
					whit_doabs = NO;
					break;
				case WHIT_FREQ:
					whit_f1 = whit_real[0];
					whit_f2 = whit_real[1];
					whit_f3 = whit_real[2];
					whit_f4 = whit_real[3];
					whit_dotaper = YES;
					break;
				case WHIT_ABS:
					whit_doabs = YES;
					break;
			}
		}
	}
}

void gsac_exec_whit(void)
{
	int k, ntrc, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float dt;
	float permin, permax;

	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
printf("WHITEN ABSOLUTE %s FREQLIMIT %s ",
		whit_doabs == YES ? "YES" : "NO",
		whit_dotaper == YES ? "YES" : "NO");
if(whit_dotaper == YES)
	printf( "%f %f %f %f\n",whit_f1,whit_f2,whit_f3,whit_f4);
else
	printf("\n");

	/* process the traces */
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		dt = sacdata[k].sachdr.rhdr[H_DELTA];
		/* filter the trace */
		gsac_dowhit(&sacdata[k].sac_data, dt, npts);
		/* update the header values */
		getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
		sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
		sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
		sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		/* set USER1 = permin USER2 = permax */
		if(whit_dotaper){
			/* if filter bound is set at default, finally use */
			if(whit_f3 > 0.0){
				permin = 1.0/whit_f3;
				/* if filter bound is set at default, finally use */
				if(sacdata[k].permin == -12345.)
					sacdata[k].permin = permin;
				/* never expand range of periods */
				if(sacdata[k].permin < permin)
					sacdata[k].permin = permin;
			}
			if(whit_f2 > 0.0){
				permax = 1.0/whit_f2;
				/* if filter bound is set at default, finally use */
				if(sacdata[k].permax == -12345.)
					sacdata[k].permax = permax;
				/* never expand range of periods */
				if(sacdata[k].permax > permax)
					sacdata[k].permax = permax;
			}
		}
	}
}


void gsac_dowhit(float **x, float dt, int npts)
{
	int i, j;
	int n2, n22;
	int jr, ji, kr, ki;
	float df, freq, f;
	float xr, xi;
	float taper;
	float ampmx, ampmn;
	int nsmooth;
	float ta;

	n2 = npow2(npts);
	/* define Nyquist index */
	n22 = n2 / 2;
	/* ensure that the temporary array is of the proper size 
	 * The y array will initially represent the complex
	 * time series, and then the complex spectra */

	if(y == (float *)NULL)
		y = (float *)calloc(2*n2,sizeof(float));
	else
		y = (float *)realloc(y,2*n2*sizeof(float));
	if(amp == (float *)NULL)
		amp = (float *)calloc(n2,sizeof(float));
	else
		amp = (float *)realloc(amp,n2*sizeof(float));
	for(i=0, j=0; i < n2 ; i++){
		/* real part of time series */
		if(i < npts )
			y[j] = (*x)[i];
		else
			y[j] = 0.0;
		/* imaginary part of time series */
		j++;
		y[j++] = 0.0;
	}
	/* get the DFT  with the proper dimensions */
	df = 0.0;
	four(y,n2,-1,&dt,&df);

	if(whit_doabs == NO){
		/* get amplitude spectrum */
		ampmx = 0.0;
		ampmn = 1.0e+37;
		for(i=0, j=0; i < n22 ; i++ , j+=2){
			jr = j   ;
			ji = j+1 ;
			amp[i] = sqrt(y[jr]*y[jr] + y[ji]*y[ji] );
		}
		/* smooth like mad */
		for(nsmooth = 0 ; nsmooth < 10 ; nsmooth++){
			for(i=0; i < n22  ; i++){
				if(i == 0)
					y[i] = (amp[i] + amp[i+1])/2.;
				else if(i == n22 -1)
					y[i] = (amp[i] + amp[i-1])/2.;
				else
					y[i] = (amp[i-1]+2.*amp[i]+amp[i+1])/4. ;
			}
			for(i=0; i < n22  ; i++){
				amp[i] = y[i];
			}
		}
		/* get the extrema */
		for(i = 0 ; i < n22 ; i++){
			if(amp[i] > ampmx) ampmx = amp[i];
			if(amp[i] < ampmn) ampmn = amp[i];
		}
		/* water level */
		if(ampmn < 0.0001 * ampmx)ampmn = 0.0001*ampmx;
	}
		
	/* recompute the Fourier transform from the original data */
	for(i=0, j=0; i < n2 ; i++){
		/* real part of time series */
		if(i < npts )
			y[j] = (*x)[i];
		else
			y[j] = 0.0;
		/* imaginary part of time series */
		j++;
		y[j++] = 0.0;
	}
	/* get the DFT  with the proper dimensions */
	df = 0.0;
	four(y,n2,-1,&dt,&df);


	for(i=0, j = 0; i <= n22 ; i++){
		freq = i*df;
		jr = j   ;
		ji = j+1 ;
		xr = y[jr];
		xi = y[ji];
		/* safety for inverse */
		if(freq == 0.0)
			f = 0.01 * df;
		else
			f = freq ;
		/* taper */
		if(whit_dotaper){
			if(freq < whit_f2)
				taper = gsac_taper3(whit_f1, whit_f2, freq);
			else if(freq > whit_f3)
				taper = gsac_taper3(whit_f4, whit_f3, freq);
			else
				taper = 1.0;
		} else {
			taper = 1.0;
		}
		if(whit_doabs == YES){
			ta = sqrt(xr*xr + xi*xi);
			if(ta > 1.0e-35 && ta < 1.0e+35){
				y[jr] = taper*xr/ta;
				y[ji] = taper*xi/ta;
			} else {
				y[jr] = taper;
				y[ji] =   0.0;
			}
		} else {
			y[jr] = taper*xr/MAX(amp[i],ampmn);
			y[ji] = taper*xi/MAX(amp[i],ampmn);
		}

		/* ensure a real time series */
		if(i > 0){
			kr = 2*(n2 - i );
			ki = kr + 1;
			y[kr] =   y[jr] ;
			y[ki] = - y[ji] ;
		}
		j+=2 ;
	}
	/* take care of Nyquist frequency n2 is the real part at the Nyquist
	 * n2+1 is the imaginary part */

	y[n2+1] = 0.0;

	/* now obtain the inverse DFT */
	four(y,n2,+1,&dt,&df);
	/* now update the time series */
	for(i=0 , j=0 ; i < npts ; i++){
		(*x)[i] = y[j];
		j+=2 ;
	}
}

