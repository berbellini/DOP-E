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

#define ROT3_TO 1
#define ROT3_SUF 2
#define ROT3_TYPE_GC 0
#define ROT3_TYPE_ANG 1
#define ROT3_TYPE_UVWSTS2 2
#define ROT3_TYPE_UVWTRIL 3
#define ROT3_TYPE_ZNE 4
static int   rot3_type;
static float rot3_ang;
static int   rot3_do ;
static char  rot3_suffix[80];
static int   rot3_dosuffix ;

struct arghdr rot3arg[] = {
	{ROT3_TO, "TO"	, CHDR, 0, 1, YES, "",-1},
	{ROT3_SUF, "SUFFIX"	, CHDR, 0, 1, YES, "",1},
	{-10, ""        , CHDR, 0, 0, YES, "" ,-1}
};


void gsac_set_param_rot3(int ncmd, char **cmdstr)
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
		should not change the fc since the np2 is wrong. One way to
		do this in the code would be to do two calls

			if(testarg(ncmd, cmdstr, cmdargs, YES) is OK
			then
				testarc,ncmd, cmdstr, cmdargs, NO)
		*/
printf("Executing rot3\n");
	if(testarg(ncmd, cmdstr, rot3arg, NO, YES))
		return;
	rot3_dosuffix = NO;
	for(i=0 ; rot3arg[i].key[0] != '\0' ; i++){
		if(rot3arg[i].used > 0){
			if(rot3arg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, rot3arg[i].key, 
					rot3arg[i].mfit, rot3arg[i].narg, instr );
				if(rot3arg[i].id == ROT3_TO){
					if(strncmp(instr,"GC",2)==0 ||strncmp(instr,"gc",2)==0){
						rot3_do = YES;
						rot3_type = ROT3_TYPE_GC ;
					} else if(strncmp(instr,"UVWSTS2",7)==0 ||strncmp(instr,"uvwsts2",7)==0){
						rot3_do = YES;
						rot3_type = ROT3_TYPE_UVWSTS2 ;
					} else if(strncmp(instr,"UVWTRIL",7)==0 ||strncmp(instr,"uvwtril",7)==0){
						rot3_do = YES;
						rot3_type = ROT3_TYPE_UVWTRIL ;
					} else if(strncmp(instr,"ZNE",3)==0 ||strncmp(instr,"zne",3)==0){
						rot3_do = YES;
						rot3_type = ROT3_TYPE_ZNE ;
					} else {
						if(isargr(instr,&rot3_ang)==YES){
							rot3_do = YES;
							rot3_type = ROT3_TYPE_ANG ;
						} else {
							rot3_do = NO;
						}
					}
				} else if(rot3arg[i].id == ROT3_SUF){
						if(strlen(instr) < 80){
						strcat(rot3_suffix,instr);
						rot3_dosuffix = YES;
						}
				}
			}
		}
	}
}

void gsac_exec_rot3(void)
{
	int i, k, ntrc, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float cmpaz[3], cmpinc[3], az[3], baz[3], delta[3];
	float o[3], a[3], t0[3];
	char cmpstr[3][100];
	double tzbeg[3], tzend[3], tzref[3];
	double ts, te;
	float CAZ0, CAZ1, CAZ2, SAZ0, SAZ1, SAZ2, C, S;
	float CIN0, CIN1, CIN2, SIN0, SIN1, SIN2;
	int ns[3];	/* index for starting fo each file */
	int month, day;
	float x0, x1, x2, y0, y1, y2;
        float a11, a12, a13, a21, a22, a23, a31, a32, a33;
	cmpstr[0][0] = '\0';
	cmpstr[1][0] = '\0';
	cmpstr[2][0] = '\0';
	/* is the command to rot3ate ? */
	if(rot3_do == NO)
		return;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1 )
		return;
	if(ntrc != 3){
		printf("Only three traces are permitted in memory\n");
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
	cmpaz[0]  = sacdata[0].sachdr.rhdr[H_CMPAZ];
	cmpinc[0] = sacdata[0].sachdr.rhdr[H_CMPINC];
	   az[0]  = sacdata[0].sachdr.rhdr[H_AZ];
	  baz[0]  = sacdata[0].sachdr.rhdr[H_BAZ];
	delta[0]  = sacdata[0].sachdr.rhdr[H_DELTA];
	tzbeg[0]  = sacdata[0].tzbeg;
	tzend[0]  = sacdata[0].tzend;
	tzref[0]  = sacdata[0].tzref;
	    o[0]  = sacdata[0].sachdr.rhdr[H_O];
	    a[0]  = sacdata[0].sachdr.rhdr[H_A];
	   t0[0]  = sacdata[0].sachdr.rhdr[H_T0];

	cmpaz[1]  = sacdata[1].sachdr.rhdr[H_CMPAZ];
	cmpinc[1] = sacdata[1].sachdr.rhdr[H_CMPINC];
	delta[1]  = sacdata[1].sachdr.rhdr[H_DELTA];
	   az[1]  = sacdata[1].sachdr.rhdr[H_AZ];
	  baz[1]  = sacdata[1].sachdr.rhdr[H_BAZ];
	tzbeg[1]  = sacdata[1].tzbeg;
	tzend[1]  = sacdata[1].tzend;
	tzref[1]  = sacdata[1].tzref;
	    o[1]  = sacdata[1].sachdr.rhdr[H_O];
	    a[1]  = sacdata[1].sachdr.rhdr[H_A];
	   t0[1]  = sacdata[1].sachdr.rhdr[H_T0];

	cmpaz[2]  = sacdata[2].sachdr.rhdr[H_CMPAZ];
	cmpinc[2] = sacdata[2].sachdr.rhdr[H_CMPINC];
	delta[2]  = sacdata[2].sachdr.rhdr[H_DELTA];
	   az[2]  = sacdata[2].sachdr.rhdr[H_AZ];
	  baz[2]  = sacdata[2].sachdr.rhdr[H_BAZ];
	tzbeg[2]  = sacdata[2].tzbeg;
	tzend[2]  = sacdata[2].tzend;
	tzref[2]  = sacdata[2].tzref;
	    o[2]  = sacdata[2].sachdr.rhdr[H_O];
	    a[2]  = sacdata[2].sachdr.rhdr[H_A];
	   t0[2]  = sacdata[2].sachdr.rhdr[H_T0];


	if(delta[0] != delta[1] || delta[0] != delta[2]){
		printf("Sample rates DELTA different\n");
		return;
	}
	if(ABS( az[0] - az[1]) > 0.01 || ABS( az[0] - az[2]) > 0.01){
		printf("Azimuths AZ different\n");
		return;
	}
	if(ABS( baz[0] - baz[1]) > 0.01 || ABS( baz[0] - baz[2]) > 0.01){
		printf("Back azimuths BAZ different\n");
		return;
	}
	ts = MAX(tzbeg[2],MAX(tzbeg[0], tzbeg[1]));
	te = MIN(tzend[2],MIN(tzend[0], tzend[1]));

	if(ts > te){
		printf("Cannot rotate - no trace overlap\n");
		return;
	}
	ns[0] = (int)((ts - tzbeg[0])/delta[0] + 0.49);
	ns[1] = (int)((ts - tzbeg[1])/delta[0] + 0.49);
	ns[2] = (int)((ts - tzbeg[2])/delta[0] + 0.49);
	npts  = (int)((te -ts)/delta[0] + 0.49) + 1;

	/* define the naming convention of the output file */
	/* The SAC file definition of CMPINC and CMPAZ is that these are the
	 * directions of positive motion. CMPINC = 0 corresponds to up, = 90
	 * corresponds to horizontal, and 180 = down. CMPAZ is measured from north
	 * with these definitions, 
	 * Z (pos up): CMPINC 0 CMPAZ not_used
	 * N (pos North): CMPINC 90 CMPAZ 0
	 * E (pos East ): CMPINC 90 CMPAZ 90
	 * For an arbitrary CMPINC and CMPAZ the positive amplitude corresponds
	 * to the following vector in (N, E, Z) vector space
	 * (cos Az sin Inc, sin Az, sin Inc, cos Inc)
	 *
	 * The net N, E, Z is a summation of the three trace projections onto N, E, Z
	 * To rotate to an arbitrary horizontal azimuth RANG (rot3_ang) , we must just do 
	 * another projection of the horizontals
	 *
	 * H(RANG)     =  C N + S E + 0 Z
	 * H(RANG+90)  = -S N + C E + 0 Z
	 * Z           =              1 Z
         * [30 MAR 2009] to generalize this to permit computation of UVW we write this as
           y0           = a12 N + a11 E + a13 Z
           y1           = a22 N + a21 E + a23 Z
           y2           = a32 N + a31 E + a33 Z
           where the a matrix elements are 
              a12 =  C, a11 = S, a13 = 0
              a22 = -S, a21 = C, a23 = 0
              a32 =  0, a31 = 0, a23 = 1
           for rotation to  R(y0), T(y1) Z(y2) in 'rotate to gc' 
           for rotation to  ANG(y0), (ANG+90)(y1) Z(y2) in 'rotate to ANG' 
           and
              a11 = 2 a12 = 0 a13=sqrt(2)
	      a21 = -1, a22 = sqrt(3) a23=sqrt(2)
              a31 = -1, a32 = -sqrt(3), a33=sqrt(2)
           for rotation to  U(y0), V(y1) W(y2) in 'rotate to uvw' 
	 *
	 * Note that for ROTATE3 TO GC, rot3_ang is BAZ + 180
	 * */
	
	switch(rot3_type){
		case ROT3_TYPE_GC:
			strcat(cmpstr[0],"R\0");
			strcat(cmpstr[1],"T\0");
			strcat(cmpstr[2],"Z\0");
			rot3_ang = baz[0] + 180.0;
			printf("Rotating to great circle ");
			printf("to form %s, %s and %s in this order\n",
				cmpstr[0],cmpstr[1],cmpstr[2]);
			C = cos(rot3_ang * 3.1415927/180.0);
			S = sin(rot3_ang * 3.1415927/180.0);
                        a11 = S;
			a12 = C;
			a13 = 0.0;
			a21 = a12;
			a22 = -S;
			a23 = 0.0;
			a31 = a13;
			a32 = a23;
			a33 = 1.0;
			sacdata[0].sachdr.rhdr[H_CMPAZ] = 
				fmod(rot3_ang,       360.);
			sacdata[1].sachdr.rhdr[H_CMPAZ] = 
				fmod(rot3_ang + 90., 360.);
			sacdata[2].sachdr.rhdr[H_CMPAZ] = 0.0;
			sacdata[0].sachdr.rhdr[H_CMPINC] = 90.0;
			sacdata[1].sachdr.rhdr[H_CMPINC] = 90.0;
			sacdata[2].sachdr.rhdr[H_CMPINC] =  0.0;
			break;
		case ROT3_TYPE_ANG:
			sprintf(cmpstr[0],"%3.3d",(int)(rot3_ang   ));
			sprintf(cmpstr[1],"%3.3d",(int)(rot3_ang+90));
			strcat(cmpstr[2],"Z\0");
			printf("Rotating to angle %f ",rot3_ang);
			printf("to form %s, %s and %s\n",
				cmpstr[0],cmpstr[1],cmpstr[2]);
			C = cos(rot3_ang * 3.1415927/180.0);
			S = sin(rot3_ang * 3.1415927/180.0);
                        a11 = S;
			a12 = C;
			a13 = 0.0;
			a21 = a12;
			a22 = -S;
			a23 = 0.0;
			a31 = a13;
			a32 = a23;
			a33 = 1.0;
			sacdata[0].sachdr.rhdr[H_CMPAZ] = 
				fmod(rot3_ang,       360.);
			sacdata[1].sachdr.rhdr[H_CMPAZ] = 
				fmod(rot3_ang + 90., 360.);
			sacdata[2].sachdr.rhdr[H_CMPAZ] = 0.0;
			sacdata[0].sachdr.rhdr[H_CMPINC] = 90.0;
			sacdata[1].sachdr.rhdr[H_CMPINC] = 90.0;
			sacdata[2].sachdr.rhdr[H_CMPINC] =  0.0;
			break;
		case ROT3_TYPE_UVWTRIL:
			strcat(cmpstr[0],"U\0");
			strcat(cmpstr[1],"V\0");
			strcat(cmpstr[2],"W\0");
			printf("Rotating to UVW(Trillium)  ");
			printf("to form %s, %s and %s in this order\n",
				cmpstr[0],cmpstr[1],cmpstr[2]);
			a11 = 2.0 / sqrt(6.0) ;
			a12 = 0.0;
			a13 = sqrt(2.0) / sqrt(6.0);
			a21 = -1.0 / sqrt(6.0) ;
			a22 = sqrt(3.0) / sqrt(6.0);
			a23 = a13;
			a31 = a21;
			a32 = -a22;
			a33 = a13;
			/* U */
			sacdata[0].sachdr.rhdr[H_CMPAZ]  =  90.0 ;
			sacdata[0].sachdr.rhdr[H_CMPINC] =  54.7;
			/* V */
			sacdata[1].sachdr.rhdr[H_CMPAZ]  = 330.0;
			sacdata[1].sachdr.rhdr[H_CMPINC] =  54.7;
			/* W */
			sacdata[2].sachdr.rhdr[H_CMPAZ]  = 210.0;
			sacdata[2].sachdr.rhdr[H_CMPINC] =  54.7;
			break;
		case ROT3_TYPE_UVWSTS2:
			strcat(cmpstr[0],"U\0");
			strcat(cmpstr[1],"V\0");
			strcat(cmpstr[2],"W\0");
			printf("Rotating to UVW(STS2)  ");
			printf("to form %s, %s and %s in this order\n",
				cmpstr[0],cmpstr[1],cmpstr[2]);
			a11 = -2.0 / sqrt(6.0) ;
			a12 = 0.0;
			a13 = sqrt(2.0) / sqrt(6.0);
			a21 =  1.0 / sqrt(6.0) ;
			a22 = sqrt(3.0) / sqrt(6.0);
			a23 = a13;
			a31 = a21;
			a32 = -a22;
			a33 = a13;
			/* U */
			sacdata[0].sachdr.rhdr[H_CMPAZ]  = 270.0 ;
			sacdata[0].sachdr.rhdr[H_CMPINC] =  54.7;
			/* V */
			sacdata[1].sachdr.rhdr[H_CMPAZ]  =  30.0;
			sacdata[1].sachdr.rhdr[H_CMPINC] =  54.7;
			/* W */
			sacdata[2].sachdr.rhdr[H_CMPAZ]  = 150.0;
			sacdata[2].sachdr.rhdr[H_CMPINC] =  54.7;
			break;
	}
	cmpaz[0] *= 3.1415927/180.0 ;
	cmpaz[1] *= 3.1415927/180.0 ;
	cmpaz[2] *= 3.1415927/180.0 ;
	cmpinc[0] *= 3.1415927/180.0 ;
	cmpinc[1] *= 3.1415927/180.0 ;
	cmpinc[2] *= 3.1415927/180.0 ;
	CAZ0 = cos(cmpaz[0]);
	CAZ1 = cos(cmpaz[1]);
	CAZ2 = cos(cmpaz[2]);
	SAZ0 = sin(cmpaz[0]);
	SAZ1 = sin(cmpaz[1]);
	SAZ2 = sin(cmpaz[2]);
	CIN0 = cos(cmpinc[0]);
	CIN1 = cos(cmpinc[1]);
	CIN2 = cos(cmpinc[2]);
	SIN0 = sin(cmpinc[0]);
	SIN1 = sin(cmpinc[1]);
	SIN2 = sin(cmpinc[2]);
	/* define the transformation matrix from x(E) y(N) z(U) to
		desired result
	*/
	/* now do the rotation */
	for(i=0 ; i < npts ; i++){
		/* get the data sample */
		x0 = sacdata[0].sac_data[i + ns[0]]; 
		x1 = sacdata[1].sac_data[i + ns[1]]; 
		x2 = sacdata[2].sac_data[i + ns[2]]; 
		y0 =  a12*(x0*CAZ0*SIN0+x1*CAZ1*SIN1+x2*CAZ2*SIN2) 
		    + a11*(x0*SAZ0*SIN0+x1*SAZ1*SIN1+x2*SAZ2*SIN2)
		    + a13*(x0*CIN0 + x1*CIN1 + x2*CIN2);
		y1 =  a22*(x0*CAZ0*SIN0+x1*CAZ1*SIN1+x2*CAZ2*SIN2) 
		    + a21*(x0*SAZ0*SIN0+x1*SAZ1*SIN1+x2*SAZ2*SIN2)
		    + a23*(x0*CIN0 + x1*CIN1 + x2*CIN2);
		y2 =  a32*(x0*CAZ0*SIN0+x1*CAZ1*SIN1+x2*CAZ2*SIN2) 
		    + a31*(x0*SAZ0*SIN0+x1*SAZ1*SIN1+x2*SAZ2*SIN2)
		    + a33*(x0*CIN0 + x1*CIN1 + x2*CIN2);
		sacdata[0].sac_data[i] = y0;
		sacdata[1].sac_data[i] = y1;
		sacdata[2].sac_data[i] = y2;
	}


	/* update headers 
	 * k=0 is for new radial or for new az
	 * k=1 is for new transverse or for new az+90 */
		for ( k=0 ; k < 3 ; k ++){
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
				cmpstr[k],sacdata[k].sac_ofile_name,rot3_dosuffix,rot3_suffix);
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

