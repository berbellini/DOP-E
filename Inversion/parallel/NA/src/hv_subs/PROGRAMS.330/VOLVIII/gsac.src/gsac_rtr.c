/* Changes:
        16 JUL 2007 - correct problems with numbers by
		using x = i/npts instead of i for the
		time variable
	`	Thus sumxx += i*i;  is replaced by sumxx += x*x
*/

#include	<stdio.h>
#include	"gsac.h"
#include	"gsac_sac.h"
#include	 "gsac_docommand.h"

extern struct sacfile_ *sacdata;

void gsac_set_param_rtr(int ncmd, char **cmdstr)
{
	int i;
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	if(gsac_control.prs > 0){
		if(gsac_control.prshist == NULL)
			gsac_control.prshist = fopen("prshist.tmp","w+");
		fprintf(gsac_control.prshist,"rtr \n");
		fflush(gsac_control.prshist);
	}
}

void gsac_exec_rtr(void)
{

	int i, k, ntrc, npts;
	double sumn, sumx, sumxx, sumy, sumxy;
	float depmax, depmin, depmen;
	int indmax, indmin;
	double yval;
	double det, a, b;
	double x;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
		
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		if(npts > 0){
			sumn  = 0.0;
			sumx  = 0.0;
			sumy  = 0.0;
			sumxx = 0.0;
			sumxy = 0.0;
			for(i=0; i < npts ; i++){
				x = (double)i/(double)npts;
				yval = sacdata[k].sac_data[i];
				sumy  += yval;
				sumxy += x*yval;
				sumn  += 1.0;
				sumx  += x;
				sumxx += x*x;
			}
			det = sumn*sumxx - sumx *sumx;
			if(det == 0.0) {
				a = 0.0;
				b = 0.0;
			} else {
				a = ( sumxx * sumy -  sumx*sumxy)/det;
				b = (-sumx  * sumy +  sumn*sumxy)/det;
			}

			for(i=0; i < npts ; i++){
				x = (double)i/(double)npts;
				sacdata[k].sac_data[i] -= (a + b*x);
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
