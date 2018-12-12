#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	XGRID_DFLT	0
#define	XGRID_ON		1
#define	XGRID_OFF	2
#define	XGRID_SOLID	3
#define	XGRID_DOTTED	4
#define	XGRID_COLOR	5
#define	XGRID_MINOR	6


struct arghdr xgridarg[] = {
	{XGRID_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{XGRID_ON  , "ON"  , IHDR, NO, 0, NO, "ON", -1},
	{XGRID_OFF  , "OFF"  , IHDR, NO, 0, NO, "OFF", -1},
	{XGRID_SOLID, "SOLID", IHDR, NO, 0, NO, "Solid", 1},
	{XGRID_DOTTED, "DOTTED", IHDR, NO, 0, NO, "Dotted", 1},
	{XGRID_COLOR, "COLOR", IHDR, NO, 1, NO, "Color int_value", 1},
	{XGRID_MINOR, "MINOR", YHDR, NO, 1, NO, "Minor [ON|OFF]", 1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float xgrid_real[10];
int   xgrid_int [10];
int   xgrid_yn;
int   xgrid_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_xgrid(int ncmd, char **cmdstr)
{
	int i;
	/* initial debug */
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, xgridarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; xgridarg[i].key[0] != '\0' ; i++){
		if(xgridarg[i].used > 0){
			if(xgridarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, xgridarg[i].key,
					xgridarg[i].mfit,xgridarg[i].narg, xgrid_int );
			} else if(xgridarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, xgridarg[i].key, 
					xgridarg[i].mfit,xgridarg[i].narg, &xgrid_yn );
                        }
			switch(xgridarg[i].id){
				case XGRID_COLOR:
					if(xgrid_int[0] < 0){
					gsac_control.xgrid_color = 1030;
					}else if(xgrid_int[0] > 1100){
					gsac_control.xgrid_color = 1030;
					}else{
					gsac_control.xgrid_color = xgrid_int[0];
					}
					break;
				case XGRID_DFLT:
					gsac_control.xgrid = NO;
					gsac_control.xgrid_type = 2;
					gsac_control.xgrid_color = 1030;
					gsac_control.xgrid_minor = YES;
					break;
				case XGRID_ON:
					gsac_control.xgrid = YES;
					break;
				case XGRID_OFF:
					gsac_control.xgrid = NO;
					break;
				case XGRID_SOLID:
					gsac_control.xgrid_type = 1;
					break;
				case XGRID_DOTTED:
					gsac_control.xgrid_type = 2;
					break;
				case XGRID_MINOR:
					if(xgrid_yn == NO){
						gsac_control.xgrid_minor = NO;
					} if(xgrid_yn == YES){
						gsac_control.xgrid_minor = YES;
					}
					break;
			}
		}
	}
			
		
}

void gsac_exec_xgrid(void)
{
	/* this just sets paramters */
}
