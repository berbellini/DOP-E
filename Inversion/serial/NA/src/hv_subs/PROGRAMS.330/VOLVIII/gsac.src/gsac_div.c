#include	<stdio.h>
#include	"gsac.h"
#include	"gsac_sac.h"
#include "gsac_docommand.h"
#include "gsac_arg.h"

#define	DIV_DIV 1

struct arghdr divarg[] = {
	{DIV_DIV, "DIV", RHDR, 0, 1, NO, "DIV value",-1},
	{0,     ""              , IHDR, 0, 0, NO, "",-1}
};

/* these are temporary variables only used here */
float div_real[1];
float div_value = 1.0;

extern struct sacfile_ *sacdata;

void gsac_set_param_div(int ncmd, char **cmdstr)
{
	float tmp;
	if (ncmd == 1)
		return;
	if(isargr(cmdstr[1], &tmp) == 1){
		if(tmp == 0.0)
			return;
		div_value = tmp;
	}
}

void gsac_exec_div(void)
{

	int i, k, ntrc, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
		
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		if(npts > 0){
			for(i=0; i < npts ; i++)
				sacdata[k].sac_data[i] /= div_value;
			getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		}
	}
}
