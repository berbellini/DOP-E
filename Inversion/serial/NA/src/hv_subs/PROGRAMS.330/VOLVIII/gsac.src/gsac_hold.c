#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	HOLD_DFLT	0
#define	HOLD_ON		1
#define	HOLD_OFF	2


struct arghdr holdarg[] = {
	{HOLD_ON   ,"ON"     , IHDR,  0, 0, NO, "", 2},
	{HOLD_OFF  ,"OFF"    , IHDR,  0, 0, NO, "", 2},
	{0,	""           , IHDR,  0, 0, NO, "",-1}
};

/* these are temporary variables only used here */
float hold_real[10];
int   hold_int [10];
int   hold_yn;
int   hold_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_hold(int ncmd, char **cmdstr)
{
	int i;
	/* implement a binary change of state of ncmd == 1 */
	if(ncmd == 1){
		return;
	}
	if(testarg(ncmd, cmdstr, holdarg, NO, YES))
		return;
	for(i=0 ; holdarg[i].key[0] != '\0' ; i++){
		if(holdarg[i].used > 0){
			if(holdarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, holdarg[i].key, 
					holdarg[i].mfit,holdarg[i].narg, hold_real);
			} else if(holdarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, holdarg[i].key, 
					holdarg[i].mfit,holdarg[i].narg, hold_int );
			} else if(holdarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, holdarg[i].key, 
					holdarg[i].mfit,holdarg[i].narg, &hold_yn );
			} else if(holdarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, holdarg[i].key, 
					holdarg[i].mfit,holdarg[i].narg, &hold_num );
			}
			switch(holdarg[i].id){
				case HOLD_ON:
					gsac_control.hold = 1 ;
					break;
				case HOLD_OFF:
					gsac_control.hold = NO ;
					if(gsac_control.inpltmode == YES){
						gend(0);
						gsac_control.inpltmode = NO;
					}
					break;
				case HOLD_DFLT:
					gsac_control.hold = NO ;
					break;
			}
		}
	}
			
		
}

void gsac_exec_hold(void)
{
	if(gsac_control.hold == NO){
		printf("Hold is OFF\n");
	} else {
		printf("Hold is ON \n");
	}
}
