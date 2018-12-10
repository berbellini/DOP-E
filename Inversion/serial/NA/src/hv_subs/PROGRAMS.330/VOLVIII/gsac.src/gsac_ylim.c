#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

#define YLIMOFF     0
#define YLIMALL     1
#define YLIMUSR     2

struct arghdr ylimarg[] = {
{YLIMOFF, "OFF"   , YHDR, 0, 0, NO, "OFF ", 1},
{YLIMALL, "ALL"   , YHDR, 0, 0, NO, "ALL  ", 1},
{YLIMUSR, "SCALE" , RHDR, 0, 2, NO, "SCALE min max  ", 1},
{0      , ""      , IHDR, 0, 0, NO, "",-1}
};
/*
{YLIMUSR, "ABSOLUTE", RHDR, 0, 2, NO, "ABSOLUTE"},
*/


float ylim_real[2];


void gsac_set_param_ylim(int ncmd, char **cmdstr)
{
	int i;
	if(ncmd == 1)
		return;
	/* is the command syntax correct ? Also reset */
	if(testarg(ncmd, cmdstr, ylimarg, NO, YES))
		return;
	for(i=0 ; ylimarg[i].key[0] != '\0' ; i++){
		/* check for special commands */
		if(ylimarg[i].used > 0){
			if(ylimarg[i].id == YLIMOFF){
				gsac_control.ylim_ctrl = YLIM_OFF;
			} else if(ylimarg[i].id == YLIMALL){
				gsac_control.ylim_ctrl = YLIM_ALL;
			} else if(ylimarg[i].id == YLIMUSR){
				getargr(ncmd, cmdstr, ylimarg[i].key, 
					ylimarg[i].mfit,ylimarg[i].narg, ylim_real);
				gsac_control.ylim_low  = ylim_real[0];
				gsac_control.ylim_high = ylim_real[1];
				gsac_control.ylim_ctrl = YLIM_USR;
			}
		}
	}
}

void gsac_exec_ylim(void)
{
	/* there is nothing to do since this just sets a flag used by
	 * gsac_plot */
	if(gsac_control.ylim_ctrl == YLIM_OFF)
		printf("YLIM is turned off\n");
	else if(gsac_control.ylim_ctrl == YLIM_ALL)
		printf("YLIM applied to all \n");
	else if(gsac_control.ylim_ctrl == YLIM_USR)
		printf("YLIM applied to all with user limits \n");

}
