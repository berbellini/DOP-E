#include	<stdio.h>
#include	"gsac.h"
#include	"gsac_sac.h"
#include "gsac_docommand.h"
#include "gsac_arg.h"

/* changes
	19 OCT 2007 error in the array indexing caught by Harley Benz.
	Change
		sacdata[k].sac_data[i] = dif;
	To
	sacdata[k].sac_data[i-1] = dif;
*/


extern struct sacfile_ *sacdata;

void gsac_set_param_dif(int ncmd, char **cmdstr)
{
	int i;
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
}

void gsac_exec_dif(void)
{

	int i, k, ntrc, npts;
	double dif;
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
			dif = 0.0;
			for(i=1; i < npts ; i++){
				dif = (sacdata[k].sac_data[i]-y0)/dt;
				y0 = sacdata[k].sac_data[i];
				sacdata[k].sac_data[i-1] = dif;
			}
			getmxmn(sacdata[k].sac_data, npts -1,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
			/* increment B by 0.5 DT, decrement E by 0.5 DT
			 * and decrease npts by 1 */
			sacdata[k].sachdr.rhdr[H_B] += 0.5*dt ;
			sacdata[k].sachdr.rhdr[H_E] -= 0.5*dt ;
			sacdata[k].sachdr.ihdr[H_NPTS] -=  1;
		}
	}
}
