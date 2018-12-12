/* PURPOSE:
	Sychronize will do the following:
	default) select the earliest absolute B time
		set this as the reference time
		adjust all time offsets, e.g., A O Tn so that
	 	   the absolute time of these markers does not change
	O)	set the reference time to the origin time
		adjust all time offsets, e.g., A O Tn so that
	 	   the absolute time of these markers does not change
		As a result ot this the O value is 0.0
	A)	set the reference time to the P  time
		adjust all time offsets, e.g., A O Tn so that
	 	   the absolute time of these markers does not change
		As a result ot this the A value is 0.0
   CHANGES:
      22 JUL 2004 - coorected synchronize - used tzb instead of tzmin for
  	dtime
      13 JUL 2010 - corrected the logic which gave bad results
      21 NOV 2011 - fixed line 207 - it is not necessary to call getmxmn since we do not change the trace
 */
#include	<stdio.h>
#include        <string.h>
#include "gsac.h"
#include "gsac_docommand.h"
#include "gsac_sac.h"
#include "gsac_sachdr.h"
#include "gsac_arg.h"
#include "csstim.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


static float  fhdr_default = -12345.;
#define SYNC_O	1
#define SYNC_A	2

struct arghdr syncarg[] = {
	{SYNC_O, "O",  IHDR, 0, 0, NO, "O", 1},
	{SYNC_A, "A",  IHDR, 0, 0, NO, "A", 1},
	{0,     ""              , IHDR, 0, 0, NO, "",-1}
};

/* these are temporary variables only used here */
int   sync_int [10];

static int sync_o = NO;
static int sync_a = NO;


void gsac_set_param_sync(int ncmd, char **cmdstr)
{

	/* initialization consists only of examining for the 
		O - origin time flag */
	int i;
	/* force the default case */
	sync_o = NO;
	sync_a = NO;
	if(ncmd == 1){
		return;
	}
	if(testarg(ncmd, cmdstr, syncarg, NO, YES)){
		return   ;
	}
	/* parse command line to see if there is an O flag */
	
	for(i=0 ; syncarg[i].key[0] != '\0' ; i++){
		if(syncarg[i].used > 0){
			switch(syncarg[i].id){
				case SYNC_O:
					sync_o = YES;
					break;
				case SYNC_A:
					sync_a = YES;
					break;
			}
		}
	}
}

/*                  B  E  O  A  T0  T1 ............................ T9 TIMMIN TIMMAX*/
static int tm[] = { 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 64,    65 };

void gsac_exec_sync(void)
{

	/* this is a two step process
	 * 1) look at all files to determine the smallest absolute B time
	 * 2) use to time to reset the reference times of all other traces
	      and to reset the other time markers, B, E, O, A, Tn if they are 
	      NOT equal to the -12345. (unknown) value 
	      */
	int k, ntrchdr, i, npts;
	double tzmin, tzb, tztmp;
	int month, day;
	double dtime;
	int otimeset;
	int atimeset;
	float otime;
	float atime;
	float depmax, depmin, depmen;
	int indmax, indmin;

	otimeset  = NO;
	atimeset  = NO;
	otime = 0.0;
	atime = 0.0;
	/* if there are no traces return */
	ntrchdr = gsac_control.number_iheaders;
	if(ntrchdr < 1)
		return;
	gsac_control.begmin= 1.0e+30;
	gsac_control.endmax=-1.0e+30;
	/* since synchronize will not change the limits of the
	   earliest and latest sample of points in epoch time,
	   we will not change the gsac_control.begmin or
	   gsac_control.endmax here */
	/* phase 1 - get the minimum absolute B time */
	for ( k=0 ; k < ntrchdr ; k ++){
		tzb = sacdata[k].tzref + sacdata[k].sachdr.rhdr[H_B];
		if(k == 0)
		 	tzmin = tzb;
		else {
			if(tzb < tzmin)tzmin = tzb;
		}
	}




	/* phase 2 - correct B E A 0 Tn if they are set and then reset
	 * reference time of the trace */
	for ( k=0 ; k < ntrchdr ; k ++){
		tztmp = sacdata[k].tzref ;
		dtime = (double)(tztmp - tzmin);
		for(i= 0 ; i < 16 ; i ++){
			/* B and E always reset */
			if ( i < 2 )
				sacdata[k].sachdr.rhdr[tm[i]] += dtime;
			else
				if(sacdata[k].sachdr.rhdr[tm[i]] != fhdr_default)
					sacdata[k].sachdr.rhdr[tm[i]] += dtime;
				

		}
                /* update with new values */
		sacdata[k].tzref = tzmin;
		sacdata[k].tzbeg=sacdata[k].tzref+sacdata[k].sachdr.rhdr[H_B];
		sacdata[k].tzend=sacdata[k].tzref+sacdata[k].sachdr.rhdr[H_E];
		sacdata[k].tzbegx = sacdata[k].tzbeg;
		sacdata[k].tzendx = sacdata[k].tzend;
		/* get bounds for absolute plotting */
		if(sacdata[k].tzbeg < gsac_control.begmin)
			gsac_control.begmin = sacdata[k].tzbeg;
		if(sacdata[k].tzend > gsac_control.endmax)
			gsac_control.endmax = sacdata[k].tzend;

		/* NOW RESET THE OTHER TIME FIELDS */
		etoh(tzmin, &sacdata[k].sachdr.ihdr[H_NZYEAR], 
			&sacdata[k].sachdr.ihdr[H_NZJDAY], &month, &day,
			&sacdata[k].sachdr.ihdr[H_NZHOUR], 
			&sacdata[k].sachdr.ihdr[H_NZMIN],
			&sacdata[k].sachdr.ihdr[H_NZSEC], 
			&sacdata[k].sachdr.ihdr[H_NZMSEC]);

		/* special case for the origin time which is position 3 
			in the array - we just save one */
		if(otimeset == NO ){
			if(sacdata[k].sachdr.rhdr[H_O] != fhdr_default){
				otimeset = YES;
				otime = sacdata[k].sachdr.rhdr[H_O] ;
			}
		} 
		if(atimeset == NO ){
			if(sacdata[k].sachdr.rhdr[H_A] != fhdr_default){
				atimeset = YES;
				atime = sacdata[k].sachdr.rhdr[H_A] ;
			}
		}
		sacdata[k].sachdr.ihdr[H_IZTYPE] = -12345 ;
	}
	/* phase 3 - if sync_o == YES or sync_a == YES then
		just adjust the headers and the NZYEAR ... NZMSEC fields
		only 
	*/
	if(sync_o == YES || sync_a == YES){
		/* for sync_o everything changes since there is 
			assumed only one origin time */
		for ( k=0 ; k < ntrchdr ; k ++){
			if(sync_o == YES){
				sacdata[k].sachdr.ihdr[H_IZTYPE] = ENUM_IO ;
				dtime = sacdata[k].sachdr.rhdr[H_O];
			} else if(sync_a == YES){
				sacdata[k].sachdr.ihdr[H_IZTYPE] = ENUM_IA ;
				dtime = sacdata[k].sachdr.rhdr[H_A];
			}
			for(i= 0 ; i < 16 ; i ++)
			  /* B and E always reset */
			  if ( i < 2 ){
			    sacdata[k].sachdr.rhdr[tm[i]] -= dtime;
			  } else {
			    if(sacdata[k].sachdr.rhdr[tm[i]] != fhdr_default)
                              sacdata[k].sachdr.rhdr[tm[i]] -= dtime;
			  }
				/* timmin and timmax are always reset */
			etoh((tzmin+dtime), &sacdata[k].sachdr.ihdr[H_NZYEAR], 
				&sacdata[k].sachdr.ihdr[H_NZJDAY], &month, &day,
				&sacdata[k].sachdr.ihdr[H_NZHOUR], 
				&sacdata[k].sachdr.ihdr[H_NZMIN],
				&sacdata[k].sachdr.ihdr[H_NZSEC], 
				&sacdata[k].sachdr.ihdr[H_NZMSEC]);
			sacdata[k].tzref += dtime;
		}
	}
}
