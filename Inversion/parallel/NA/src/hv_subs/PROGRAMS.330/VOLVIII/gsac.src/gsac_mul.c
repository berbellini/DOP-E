#include	<stdio.h>
#include	"gsac.h"
#include	"gsac_sac.h"
#include "gsac_docommand.h"
#include "gsac_arg.h"

#define	MUL_MUL 1

struct arghdr mularg[] = {
	{MUL_MUL, "MUL", RHDR, 0, 1, NO, "MUL value",-1},
	{0,     ""              , IHDR, 0, 0, NO, "",-1}
};

/* these are temporary variables only used here */
float mul_real[1];
float mul_value = 1.0;

extern struct sacfile_ *sacdata;

void gsac_set_param_mul(int ncmd, char **cmdstr)
{
	float tmp;
	if (ncmd == 1)
		return;
	if(isargr(cmdstr[1], &tmp) == 1)
		mul_value = tmp;
}

void gsac_exec_mul(void)
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
				sacdata[k].sac_data[i] *= mul_value;
			getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		}
	}
}
