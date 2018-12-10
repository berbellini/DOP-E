#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include	"csstim.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

#define ROT_TO 1
#define ROT_SUF 2
static float rot_ang;
static int   rot_gc ;
static int   rot_do ;
static char  rot_suffix[80];
static int   rot_dosuffix;

struct arghdr rotarg[] = {
	{ROT_TO, "TO"	, CHDR, 0, 1, YES, "",-1},
	{ROT_SUF, "SUFFIX", CHDR, 0, 1, YES, "",1},
	{-10, ""        , CHDR, 0, 0, YES, "" ,-1}
};


void gsac_set_param_rot(int ncmd, char **cmdstr)
{
	int i;
	char instr[80];
	/* note when the testrg routine is used, if the argument is
		NO then you must use internal variables to define the 
		state of the operation - if you use YES, then things are
		not changed until the input is proven correct. An exmple of
		this concept with YES is the following:
		Assume we wish aa LP filter with fc 1 np 2 p 1 
		If we enter  fc 2 np2   there is a syntax error and we
		should not chnge the fc since the np2 is wrong. One way to
		do this in the code would be to do two calls

			if(testarg(ncmd, cmdstr, cmdargs, YES) is OK
			then
				testarc,ncmd, cmdstr, cmdargs, NO)
		*/
	if(testarg(ncmd, cmdstr, rotarg, NO, YES))
		return;
	rot_dosuffix = NO;
	for(i=0 ; rotarg[i].key[0] != '\0' ; i++){
		if(rotarg[i].used > 0){
			if(rotarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, rotarg[i].key, 
					rotarg[i].mfit, rotarg[i].narg, instr );
				if(rotarg[i].id == ROT_TO){
					if(strncmp(instr,"GC",2)==0 ||strncmp(instr,"gc",2)==0){
						rot_gc = YES;
						rot_do = YES;
					} else {
						if(isargr(instr,&rot_ang)==YES){
							rot_gc = NO;
							rot_do = YES;
						} else {
							rot_do = NO;
						}
					}
				} else if(rotarg[i].id == ROT_SUF){
						if(strlen(instr) < 80){
						strcat(rot_suffix,instr);
						rot_dosuffix = YES;
						}
				}
			}
		}
	}
}

void gsac_exec_rot(void)
{
	int i, k, ntrc, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float cmpaz[2], az[2], baz[2], delta[2];
	float o[2], a[2], t0[2];
	char cmpstr[2][100];
	double tzbeg[2], tzend[2], tzref[2];
	double ts, te;
	float CAZ0, CAZ1, SAZ0, SAZ1, C, S;
	int ns[2];	/* index for starting fo each file */
	int month, day;
	float x0, x1, y0, y1;
	cmpstr[0][0] = '\0';
	cmpstr[1][0] = '\0';
	/* is the command to rotate ? */
	if(rot_do == NO)
		return;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if( ntrc < 1 )
		return;
	if(ntrc != 2){
		printf("Only two traces are permitted in memory\n");
		return;
	}
	gsac_control.begmin= 1.0e+30;
	gsac_control.endmax=-1.0e+30;
	/* steps
	 * get cmpinc of each trace,
	 * get ts and te of each
	 * get delta of each
	 * find common time window
	 * change ts te
	 * change component name in both headers
	 * change file name
	 * reset depmax and depmin
	 * */
	/* set parameters for checking */
	cmpaz[0] = sacdata[0].sachdr.rhdr[H_CMPAZ];
	   az[0] = sacdata[0].sachdr.rhdr[H_AZ];
	  baz[0] = sacdata[0].sachdr.rhdr[H_BAZ];
	delta[0] = sacdata[0].sachdr.rhdr[H_DELTA];
	tzbeg[0] = sacdata[0].tzbeg;
	tzend[0] = sacdata[0].tzend;
	tzref[0] = sacdata[0].tzref;
	    o[0] = sacdata[0].sachdr.rhdr[H_O];
	    a[0] = sacdata[0].sachdr.rhdr[H_A];
	   t0[0] = sacdata[0].sachdr.rhdr[H_T0];
	cmpaz[1] = sacdata[1].sachdr.rhdr[H_CMPAZ];
	delta[1] = sacdata[1].sachdr.rhdr[H_DELTA];
	   az[1] = sacdata[1].sachdr.rhdr[H_AZ];
	  baz[1] = sacdata[1].sachdr.rhdr[H_BAZ];
	tzbeg[1] = sacdata[1].tzbeg;
	tzend[1] = sacdata[1].tzend;
	tzref[1] = sacdata[1].tzref;
	    o[1] = sacdata[1].sachdr.rhdr[H_O];
	    a[1] = sacdata[1].sachdr.rhdr[H_A];
	   t0[1] = sacdata[1].sachdr.rhdr[H_T0];

	if(delta[0] != delta[1]){
		printf("Sample rates DELTA different\n");
		return;
	}
	if(ABS (az[0] - az[1]) > 0.01 ){
		printf("Azimuths AZ different\n");
		return;
	}
	if(ABS( baz[0] - baz[1] ) > 0.01 ){
		printf("Back azimuths BAZ different\n");
		return;
	}
	if(cmpaz[0] == cmpaz[1]){
		printf("Component azimuths CMPAZ are same\n");
		return;
	}
	ts = MAX(tzbeg[0], tzbeg[1]);
	te = MIN(tzend[0], tzend[1]);

	if(ts > te){
		printf("Cannot rotate - no trace overlap\n");
		return;
	}
	ns[0] = (int)((ts - tzbeg[0])/delta[0] + 0.49);
	ns[1] = (int)((ts - tzbeg[1])/delta[0] + 0.49);
	npts  = (int)((te -ts)/delta[0] + 0.49) + 1;
	
	if(rot_gc){
		strcat(cmpstr[0],"R\0");
		strcat(cmpstr[1],"T\0");
		rot_ang = baz[0] + 180.0;
		printf("Rotating to great circle ");
		printf("to form %s and %s\n",cmpstr[0],cmpstr[1]);
	} else {
		sprintf(cmpstr[0],"%3.3d",(int)(rot_ang   ));
		sprintf(cmpstr[1],"%3.3d",(int)(rot_ang+90));
		printf("Rotating to angle %f ",rot_ang);
		printf("to form %s and %s\n",cmpstr[0],cmpstr[1]);
	}
	/* to do the rotation we will use vector projection
	 * to get the vector decomponsed into [N, E]
	 * x0 = [ t0 cos(cmpaz[0]), t0 sin(cmpaz[0]) ]
	 * x1 = [ t1 cos(cmpaz[1]), t1 sin(cmpaz[1]) ]
	 * 	where t0 and t1 are the two time series and cmpaz[] is the
	 * 	direction of positive displacement.
	 * Define C = cos(rot_ang), S = sin(rot_ang), then the 
	 * vector in the rot_ang direction is [ C, S ] and the vector in the
	 * rot_ang + 90 direction is [ -S, C ] and this the
	 * two rotated traces are
	 * y0 = [ C, S] dot [ t0 cos(cmpaz[0])+t1 cos(cmpaz[1]) ,
	 * 			t0 sin(cmpaz[0])+t1 sin(cmpaz[1]) ]
	 * y1 = [-S, C] dot [ t0 cos(cmpaz[0])+t1 cos(cmpaz[1]) ,
	 * 			t0 sin(cmpaz[0])+t1 sin(cmpaz[1]) ]
	 * */
	cmpaz[0] *= 3.1415927/180.0 ;
	cmpaz[1] *= 3.1415927/180.0 ;
	CAZ0 = cos(cmpaz[0]);
	CAZ1 = cos(cmpaz[1]);
	SAZ0 = sin(cmpaz[0]);
	SAZ1 = sin(cmpaz[1]);
	C = cos(rot_ang * 3.1415927/180.0);
	S = sin(rot_ang * 3.1415927/180.0);
	/* now do the rotation */
	for(i=0 ; i < npts ; i++){
		x0 = sacdata[0].sac_data[i + ns[0]]; 
		x1 = sacdata[1].sac_data[i + ns[1]]; 
		y0 =  C*(x0*CAZ0+x1*CAZ1) + S*(x0*SAZ0+x1*SAZ1);
		y1 = -S*(x0*CAZ0+x1*CAZ1) + C*(x0*SAZ0+x1*SAZ1);
		sacdata[0].sac_data[i] = y0;
		sacdata[1].sac_data[i] = y1;
	}


	/* update headers 
	 * k=0 is for new radial or for new az
	 * k=1 is for new transverse or for new az+90 */
		sacdata[0].sachdr.rhdr[H_CMPAZ] = fmod(rot_ang,       360.);
		sacdata[1].sachdr.rhdr[H_CMPAZ] = fmod(rot_ang + 90., 360.);
		for ( k=0 ; k < 2 ; k ++){
		if(npts > 0){
			getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.ihdr[H_NPTS] = npts;
			sacdata[k].sachdr.rhdr[H_B] = ts - tzref[0];
			sacdata[k].sachdr.rhdr[H_E] = te - tzref[0];
			sacdata[k].sachdr.rhdr[H_O] = o[0];
			sacdata[k].sachdr.rhdr[H_A] = a[0];
			sacdata[k].sachdr.rhdr[H_T0] = t0[0];
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
			/* update the headers */
		etoh(tzref[0], &sacdata[k].sachdr.ihdr[H_NZYEAR], 
			&sacdata[k].sachdr.ihdr[H_NZJDAY], &month, &day,
			&sacdata[k].sachdr.ihdr[H_NZHOUR], 
			&sacdata[k].sachdr.ihdr[H_NZMIN],
			&sacdata[k].sachdr.ihdr[H_NZSEC], 
			&sacdata[k].sachdr.ihdr[H_NZMSEC]);
		}
		/* redefine component names and default file name */
		strcpy(sacdata[k].sac_ofile_name, sacdata[k].sac_ifile_name);
		strcpy(sacdata[k].schdr[H_KCMPNM], sacdata[k].ocmpnm);
		chcmpnm(sacdata[k].schdr[H_KSTNM],sacdata[k].schdr[H_KCMPNM], 
				cmpstr[k],sacdata[k].sac_ofile_name,rot_dosuffix,rot_suffix);
		strncpy(sacdata[k].sachdr.chdr[H_KCMPNM], sacdata[k].schdr[H_KCMPNM],8);

		/* set epoch times of first and last sample */
	       	sacdata[k].tzref = tzref[0] ;
		sacdata[k].tzbeg=sacdata[k].tzref+sacdata[k].sachdr.rhdr[H_B];
		sacdata[k].tzend=sacdata[k].tzref+sacdata[k].sachdr.rhdr[H_E];
		sacdata[k].tzbegx = sacdata[k].tzbeg;
		sacdata[k].tzendx = sacdata[k].tzend;
		/* get bounds for absolute plotting */
		/*
		if(sacdata[k].tzbeg < gsac_control.begmin)
		*/
			gsac_control.begmin = sacdata[k].tzbeg;
		/*
		if(sacdata[k].tzend > gsac_control.endmax)
		*/
			gsac_control.endmax = sacdata[k].tzend;


	}
}

