#include	<stdio.h>
#include	"gsac.h"
#include	"gsac_sac.h"
#include "gsac_docommand.h"
#include "gsac_arg.h"

extern struct sacfile_ *sacdata;

void gsac_set_param_int(int ncmd, char **cmdstr)
{
	int i;
	for(i=0; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
}

void gsac_exec_int(void)
{

	int i, k, ntrc, npts;
	double sum;
	float y0, dt;
	float depmax, depmin, depmen;
	int indmax, indmin;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
		
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		dt = sacdata[k].sachdr.rhdr[H_DELTA];
		if(npts > 0){
			y0 = sacdata[k].sac_data[0];
			sum = 0.0;
			y0 = 0.0;
			for(i=0; i < npts ; i++){
				sum += (sacdata[k].sac_data[i]-y0)*dt;
				sacdata[k].sac_data[i] = sum;
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
