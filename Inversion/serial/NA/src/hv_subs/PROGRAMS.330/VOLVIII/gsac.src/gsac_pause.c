#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

#include <unistd.h>

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	PAUSE_DFLT	0
#define	PAUSE_PERIOD	1


struct arghdr pausearg[] = {
	{PAUSE_DFLT    , "DEFAULT", IHDR, 0, 0, NO, "", 1},
	{PAUSE_PERIOD  , "PERIOD" , IHDR,  0, 1, NO, "PERIOD delay", 1},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};

/* these are temporary variables only used here */
static float pause_real[10];
static int   pause_int [10];
static int   pause_yn;
static int   pause_num;

static int pause_delay;



/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_pause(int ncmd, char **cmdstr)
{
	int i;
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, pausearg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; pausearg[i].key[0] != '\0' ; i++){
		if(pausearg[i].used > 0){
			if(pausearg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, pausearg[i].key, 
					pausearg[i].mfit, pausearg[i].narg, pause_real);
			} else if(pausearg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, pausearg[i].key, 
					pausearg[i].mfit, pausearg[i].narg, pause_int );
			} else if(pausearg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, pausearg[i].key, 
					pausearg[i].mfit, pausearg[i].narg, &pause_yn );
			} else if(pausearg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, pausearg[i].key, 
					pausearg[i].mfit, pausearg[i].narg, &pause_num );
			}


			switch(pausearg[i].id){
				case PAUSE_DFLT:
					pause_delay = 1;
					break;
				case PAUSE_PERIOD:
					pause_delay = pause_int[0];
					break;

			}
		}
	}
			
		
}


void gsac_exec_pause(void)
{
	/*
	char c;
	printf("Pause: Enter/return to continue\n");
	scanf("%c", &c);
	*/
	sleep(pause_delay);
}
