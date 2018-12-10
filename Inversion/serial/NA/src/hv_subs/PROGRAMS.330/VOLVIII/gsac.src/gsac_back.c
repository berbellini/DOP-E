#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	BACK_DFLT	0
#define	BACK_ON		1
#define	BACK_OFF	2
#define	BACK_COLOR	5


struct arghdr backarg[] = {
	{BACK_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{BACK_ON  , "ON"  , IHDR, NO, 0, NO, "ON", -1},
	{BACK_OFF  , "OFF"  , IHDR, NO, 0, NO, "OFF", -1},
	{BACK_COLOR, "COLOR", IHDR, NO, 1, NO, "Color int_value", 1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float back_real[10];
int   back_int [10];
int   back_yn;
int   back_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_back(int ncmd, char **cmdstr)
{
	int i;
	/* initial debug */
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, backarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; backarg[i].key[0] != '\0' ; i++){
		if(backarg[i].used > 0){
			if(backarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, backarg[i].key,
					backarg[i].mfit,backarg[i].narg, back_int );
                        }

			switch(backarg[i].id){
				case BACK_COLOR:
					if(back_int[0] < 0){
					gsac_control.background_color = 0;
					} else if(back_int[0] > 1100){
					gsac_control.background_color = 0;
					} else{
					gsac_control.background_color = back_int[0];
					}
					break;
				case BACK_DFLT:
					gsac_control.background_color = 0;
					gsac_control.background = NO;
					break;
				case BACK_ON:
					gsac_control.background = YES;
					break;
				case BACK_OFF:
					gsac_control.background = NO ;
					break;
			}
		}
	}
			
		
}

void gsac_exec_back(void)
{
	/* this just sets paramters */
}
