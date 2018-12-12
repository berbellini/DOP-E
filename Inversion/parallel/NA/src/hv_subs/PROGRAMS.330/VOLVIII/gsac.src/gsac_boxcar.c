#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

/* Changes:
	14 APR 2007 - ensure that if HALFWIDTH < DT, nothing is done */


#define	BOXCAR_DFLT	0
#define	BOXCAR_WIDTH		1


struct arghdr boxcararg[] = {
	{BOXCAR_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{BOXCAR_WIDTH  , "WIDTH"  , RHDR, NO, 1, NO, "Width width", 1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
static float boxcar_real[10];
static int   boxcar_int [10];
static float boxcar_width;
static int   boxcar_yn;
static int   boxcar_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_boxcar(int ncmd, char **cmdstr)
{
	int i;
	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, boxcararg, NO, YES))
		return;
	/* parse commands */
	boxcar_width = -1. ;
	for(i=0 ; boxcararg[i].key[0] != '\0' ; i++){
		if(boxcararg[i].used > 0){
			if(boxcararg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, boxcararg[i].key, 
					boxcararg[i].mfit,boxcararg[i].narg, boxcar_real);
			} else if(boxcararg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, boxcararg[i].key, 
					boxcararg[i].mfit,boxcararg[i].narg, boxcar_int );
			} else if(boxcararg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, boxcararg[i].key, 
					boxcararg[i].mfit,boxcararg[i].narg, &boxcar_yn );
			} else if(boxcararg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, boxcararg[i].key, 
					boxcararg[i].mfit,boxcararg[i].narg, &boxcar_num );
			}
			switch(boxcararg[i].id){
				case BOXCAR_WIDTH:
					boxcar_width = boxcar_real[0];
					break;

			}
		}
	}
			
		
}

void gsac_exec_boxcar(void)
{
	if(boxcar_width > 0.0)
		gsac_pulconv(0.0,boxcar_width,0.0);
}
