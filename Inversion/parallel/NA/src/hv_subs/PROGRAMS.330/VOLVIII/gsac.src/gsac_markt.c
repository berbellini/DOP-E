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


#define	MARKT_DFLT	0
#define	MARKT_VEL	1
#define	MARKT_DIST	2
#define	MARKT_OGMT	3
#define	MARKT_OCAL	4
#define MARKT_ON        5
#define MARKT_OFF       6

#define NUMMARKT 11

struct arghdr marktarg[] = {
	{MARKT_DFLT,"DEFAULT", IHDR, 0, 0, NO, "", 2},
	{MARKT_VEL ,"X0"  , RHDR, 0, 1, NO, "X0 x0", -1},
	{MARKT_DIST,"DISTANCE"  , RHDR, 0, 1, NO, "Y0 y0", 1},
	{MARKT_DIST,"D"  , RHDR, 0, 1, NO, "Y0 y0", -1},
	{MARKT_OGMT,"OGMT", IHDR, 0, 6, YES, "OGMT YY DOY HH MM SS MSEC",-1},
	{MARKT_OCAL,"OCAL", IHDR, 0, 7, YES, "OCAL YY MM DD HH MM SS MSEC",-1},
	{MARKT_ON  ,"ON", IHDR, 0, 0, NO, "", -1},
	{MARKT_OFF ,"OFF", IHDR, 0, 0, NO, "", -1},
	{0,	""	, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float markt_real[NUMMARKT];
int   markt_int [NUMMARKT];
int   markt_yn;
int   markt_num;

/* these are prototypes for global variables to be used by the routine */
int markt_on ;
int markt_doo = NO;
int markt_dod = NO;
int markt_numvel = 8;
float markt_dist ;
double markt_o ;
float markt_vel[NUMMARKT] = {1., 2., 3., 3.5, 4., 5., 6., 7., 8., -1., -1.};


void gsac_set_param_markt(int ncmd, char **cmdstr)
{
	int i;
	int numgmt, numcal;
	numgmt = gsac_countgmt(ncmd, cmdstr,"GMT") ;
	if(numgmt > 0){
		for(i=0;i<numgmt;i++){
			gsac_mergegmt(&ncmd,cmdstr,"GMT");
		}
	}
	numcal = gsac_countgmt(ncmd, cmdstr,"CAL") ;
	if(numcal > 0){
		for(i=0;i<numgmt;i++){
			gsac_mergegmt(&ncmd,cmdstr,"CAL");
		}
	}
	if(numgmt > 0 || numcal > 0){
		printf("Executing: ");
		for(i=0;i<ncmd;i++)
			printf("%s ",cmdstr[i]);
		printf("\n");
	}
	/* initialize */
	markt_on = YES;
	
	
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, marktarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; marktarg[i].key[0] != '\0' ; i++){
		if(marktarg[i].used > 0){
			if(marktarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, marktarg[i].key, 
					marktarg[i].mfit,marktarg[i].narg, markt_real);
			} else if(marktarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, marktarg[i].key, 
					marktarg[i].mfit,marktarg[i].narg, markt_int );
			} else if(marktarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, marktarg[i].key, 
					marktarg[i].mfit,marktarg[i].narg, &markt_yn );
			} else if(marktarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, marktarg[i].key, 
					marktarg[i].mfit,marktarg[i].narg, &markt_num );
			}
			switch(marktarg[i].id){
				case MARKT_DFLT:
					markt_doo = NO;
					markt_dod = NO;
					markt_numvel =  9 ;
					markt_vel[0] =  1.;
					markt_vel[1] =  2.;
					markt_vel[2] =  3.;
					markt_vel[3] =  3.5;
					markt_vel[4] =  4.;
					markt_vel[5] =  5.;
					markt_vel[6] =  6.;
					markt_vel[7] =  7.;
					markt_vel[8] =  8.;
					markt_vel[9] = -1.;
					markt_vel[10] = -1.;
					markt_on = YES;
					break;
				case MARKT_ON:
					markt_on = YES;
					break;
				case MARKT_OFF:
					markt_on = NO;
					break;
				case MARKT_DIST:
					markt_dist = markt_real[0];
					break;
				case MARKT_VEL:
					break;
				case MARKT_OGMT:
					htoe1(  markt_int[0], /* year */
						markt_int[1], /* doy  */
						markt_int[2], /* hour */
						markt_int[3], /* minute */
						markt_int[4], /* second */
						markt_int[5], /* msec */
						&markt_o);
					break;
				case MARKT_OCAL:
					htoe2(  markt_int[0], /* year */
						markt_int[1], /* month  */
						markt_int[2], /* day  */
						markt_int[3], /* hour */
						markt_int[4], /* minute */
						markt_int[5], /* second */
						markt_int[6], /* msec */
						&markt_o);
					break;
			}
		}
	}
			
		
}

void gsac_exec_markt(void)
{
	/* nothing is done since the values are passed to the gsac_plotsub.c routine */
}
