/* CHANGES
	06 JUN 2010 - created
*/

#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	SHIFT_DFLT	0
#define	SHIFT_FIXEDAMOUNT	1

static float shiftfixed;
static int doshiftfixed;


struct arghdr shiftarg[] = {
	{SHIFT_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{SHIFT_FIXEDAMOUNT, "FIXED"  , RHDR, NO, 1, NO, "Fixed amount ", 1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float shift_real[10];
int   shift_int [10];
int   shift_yn;
int   shift_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_shift(int ncmd, char **cmdstr)
{
	int i;
	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, shiftarg, NO, YES))
		return;
	/* parse commands */
        doshiftfixed  = NO;
	for(i=0 ; shiftarg[i].key[0] != '\0' ; i++){
		if(shiftarg[i].used > 0){
			if(shiftarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, shiftarg[i].key, 
					shiftarg[i].mfit,shiftarg[i].narg, shift_real);
			} else if(shiftarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, shiftarg[i].key, 
					shiftarg[i].mfit,shiftarg[i].narg, shift_int );
			} else if(shiftarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, shiftarg[i].key, 
					shiftarg[i].mfit,shiftarg[i].narg, &shift_yn );
			} else if(shiftarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, shiftarg[i].key, 
					shiftarg[i].mfit,shiftarg[i].narg, &shift_num );
			}
			switch(shiftarg[i].id){
				case SHIFT_FIXEDAMOUNT:
					shiftfixed = shift_real[0] ;
					doshiftfixed = YES;
					break;

			}
		}
	}
			
		
}

void gsac_exec_shift(void)
{
	int i, k, ntrc, npts;
	float depmax, depmin, depmen;
	int indmax, indmin;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
		
	for ( k=0 ; k < ntrc ; k ++){
                if(doshiftfixed == YES){
			sacdata[k].sachdr.rhdr[H_B] += shiftfixed;
			sacdata[k].sachdr.rhdr[H_E] += shiftfixed;
                        if(sacdata[k].sachdr.rhdr[H_A] != -12345.)
				sacdata[k].sachdr.rhdr[H_A] += shiftfixed;
                        if(sacdata[k].sachdr.rhdr[H_T0] != -12345.)
				sacdata[k].sachdr.rhdr[H_T0] += shiftfixed;
                        if(sacdata[k].sachdr.rhdr[H_T1] != -12345.)
				sacdata[k].sachdr.rhdr[H_T1] += shiftfixed;
                        if(sacdata[k].sachdr.rhdr[H_T2] != -12345.)
				sacdata[k].sachdr.rhdr[H_T2] += shiftfixed;
                        if(sacdata[k].sachdr.rhdr[H_T3] != -12345.)
				sacdata[k].sachdr.rhdr[H_T3] += shiftfixed;
                        if(sacdata[k].sachdr.rhdr[H_T4] != -12345.)
				sacdata[k].sachdr.rhdr[H_T4] += shiftfixed;
                        if(sacdata[k].sachdr.rhdr[H_T5] != -12345.)
				sacdata[k].sachdr.rhdr[H_T5] += shiftfixed;
                        if(sacdata[k].sachdr.rhdr[H_T6] != -12345.)
				sacdata[k].sachdr.rhdr[H_T6] += shiftfixed;
                        if(sacdata[k].sachdr.rhdr[H_T7] != -12345.)
				sacdata[k].sachdr.rhdr[H_T7] += shiftfixed;
                        if(sacdata[k].sachdr.rhdr[H_T8] != -12345.)
				sacdata[k].sachdr.rhdr[H_T8] += shiftfixed;
                        if(sacdata[k].sachdr.rhdr[H_T9] != -12345.)
				sacdata[k].sachdr.rhdr[H_T9] += shiftfixed;
			sacdata[k].tzbeg = 
				sacdata[k].tzref + sacdata[k].sachdr.rhdr[H_B];
			sacdata[k].tzend = 
				sacdata[k].tzref + sacdata[k].sachdr.rhdr[H_E];
			sacdata[k].tzbegx = sacdata[k].tzbeg;
			sacdata[k].tzendx = sacdata[k].tzend;
			/* get bounds for absolute plotting this is easy 
				since there is only one in memory at a time */
				gsac_control.begmin = sacdata[k].tzbeg;
				gsac_control.endmax = sacdata[k].tzend;
		}
	}
}
