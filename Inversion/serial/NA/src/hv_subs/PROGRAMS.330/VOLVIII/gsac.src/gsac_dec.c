#include        <stdio.h>
#include        "gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

/* implement INTERPOLATE 
	CHANGES: routine inter re-written 26 May 2006 INGV Roma
	because logic was bad. Also cubic interpolation works at
	the ends instead of linear.
*/

extern struct sacfile_ *sacdata;

#define	IN_NDEC 0

int dec_ndec = 0;



struct arghdr decarg[] = {
	{IN_NDEC, "NDEC", IHDR, 0, 1, YES, "NDEC ndec", 2},
	{-10, ""        , CHDR, 0, 0, YES, "" ,-1}
};

/* these are temporary variables only used here */
int dec_int[10];
int do_dec;


void gsac_set_param_dec(int ncmd, char **cmdstr)
{
	int itmp;
	dec_ndec = 0;
	do_dec = NO;
	if (ncmd == 1)
		return;
	if(isargi(cmdstr[1], &itmp) == 1){
		dec_ndec = itmp;;
		do_dec = YES;
	}

}

void gsac_exec_dec(void)
{
	int i, k, ntrc, npts, mpts;
	int j;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float dt;

	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
	if(do_dec != YES)
		return;

	gsac_control.begmin= 1.0e+30;
	gsac_control.endmax=-1.0e+30;

	/* process the traces */
	/* note we do not have to reallocate since the number of points
		always decreases*/
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		mpts = npts/dec_ndec;
		if(mpts < 1)
			mpts = 1;
		for(i=0,j=0; i < mpts ; i++,j+=dec_ndec)
			sacdata[k].sac_data[i] = sacdata[k].sac_data[j];

		/* update the header values */
		getmxmn(sacdata[k].sac_data, mpts,&depmax, &depmin, &depmen,&indmax,&indmin);
		sacdata[k].sachdr.ihdr[H_NPTS] = mpts;
		dt = dec_ndec*sacdata[k].sachdr.rhdr[H_DELTA];
		sacdata[k].sachdr.rhdr[H_DELTA] = dt;
		sacdata[k].sachdr.rhdr[H_E] = sacdata[k].sachdr.rhdr[H_B] +
				(mpts -1 )*dt;
		sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
		sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
		sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;

		sacdata[k].tzbeg=sacdata[k].tzref+sacdata[k].sachdr.rhdr[H_B];
		sacdata[k].tzend=sacdata[k].tzref+sacdata[k].sachdr.rhdr[H_E];
		/* get bounds for absolute plotting */
		sacdata[k].tzbegx = sacdata[k].tzbeg;
		sacdata[k].tzendx = sacdata[k].tzend;
		if(sacdata[k].tzbeg < gsac_control.begmin)
			gsac_control.begmin = sacdata[k].tzbeg;
		if(sacdata[k].tzend > gsac_control.endmax)
				gsac_control.endmax = sacdata[k].tzend;
	}

}

