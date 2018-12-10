#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	YGRID_DFLT	0
#define	YGRID_ON	1
#define	YGRID_OFF	2
#define	YGRID_SOLID	3
#define	YGRID_DOTTED	4
#define	YGRID_COLOR	5
#define	YGRID_MINOR	6


struct arghdr ygridarg[] = {
	{YGRID_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{YGRID_ON  , "ON"  , IHDR, NO, 0, NO, "ON ", -1},
	{YGRID_OFF  , "OFF"  , IHDR, NO, 0, NO, "OFF", -1},
	{YGRID_SOLID, "SOLID", IHDR, NO, 0, NO, "Solid", 1},
	{YGRID_DOTTED, "DOTTED", IHDR, NO, 0, NO, "Dotted", 1},
	{YGRID_COLOR, "COLOR", IHDR, NO, 1, NO, "Color int_value", 1},
	{YGRID_MINOR, "MINOR", YHDR, NO, 1, NO, "Minor [ON|OFF]", 1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float ygrid_real[10];
int   ygrid_int [10];
int   ygrid_yn;
int   ygrid_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_ygrid(int ncmd, char **cmdstr)
{
	int i;
	/* initial debug */
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, ygridarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; ygridarg[i].key[0] != '\0' ; i++){
		if(ygridarg[i].used > 0){
			if(ygridarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, ygridarg[i].key,
					ygridarg[i].mfit,ygridarg[i].narg, ygrid_int );
			} else if(ygridarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, ygridarg[i].key, 
					ygridarg[i].mfit,ygridarg[i].narg, &ygrid_yn );
                        }
			switch(ygridarg[i].id){
				case YGRID_COLOR:
					if(ygrid_int[0] < 0) {
					gsac_control.ygrid_color = 4;
					} else if(ygrid_int[0] > 1100) {
					gsac_control.ygrid_color = 1030;
					} else {
					gsac_control.ygrid_color = ygrid_int[0];
					}
					break;
				case YGRID_DFLT:
					gsac_control.ygrid = NO;
					gsac_control.ygrid_type = 2;
					gsac_control.ygrid_color = 1030;
					gsac_control.ygrid_minor = NO;
					break;
				case YGRID_ON:
					gsac_control.ygrid = YES;
					break;
				case YGRID_OFF:
					gsac_control.ygrid = NO;
					break;
				case YGRID_SOLID:
					gsac_control.ygrid_type = 1;
					break;
				case YGRID_DOTTED:
					gsac_control.ygrid_type = 2;
					break;
				case YGRID_MINOR:
					if(ygrid_yn == NO){
						gsac_control.ygrid_minor = NO;
					} if(ygrid_yn == YES){
						gsac_control.ygrid_minor = YES;
					}
					break;
			}
		}
	}
			
		
}

void gsac_exec_ygrid(void)
{
	/* this just sets paramters */
}
