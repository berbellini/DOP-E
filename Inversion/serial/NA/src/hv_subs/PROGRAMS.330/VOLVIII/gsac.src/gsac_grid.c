#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	GRID_DFLT	0
#define	GRID_ON		1
#define	GRID_OFF	2
#define	GRID_SOLID	3
#define	GRID_DOTTED	4
#define	GRID_COLOR	5
#define	GRID_MINOR	6


struct arghdr gridarg[] = {
	{GRID_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{GRID_ON  , "ON"  , IHDR, NO, 0, NO, "ON", -1},
	{GRID_OFF  , "OFF"  , IHDR, NO, 0, NO, "OFF", -1},
	{GRID_SOLID, "SOLID", IHDR, NO, 0, NO, "Solid", 1},
	{GRID_DOTTED, "DOTTED", IHDR, NO, 0, NO, "Dotted", 1},
	{GRID_COLOR, "COLOR", IHDR, NO, 1, NO, "Color int_value", 1},
	{GRID_MINOR, "MINOR", YHDR, NO, 1, NO, "Minor [ON|OFF]", 1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float grid_real[10];
int   grid_int [10];
int   grid_yn;
int   grid_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_grid(int ncmd, char **cmdstr)
{
	int i;
	/* initial debug */
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, gridarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; gridarg[i].key[0] != '\0' ; i++){
		if(gridarg[i].used > 0){
			if(gridarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, gridarg[i].key,
					gridarg[i].mfit,gridarg[i].narg, grid_int );
			} else if(gridarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, gridarg[i].key, 
					gridarg[i].mfit,gridarg[i].narg, &grid_yn );
                        }

			switch(gridarg[i].id){
				case GRID_COLOR:
					if(grid_int[0] < 0){
					gsac_control.xgrid_color = 1030;
					gsac_control.ygrid_color = 1030;
					} else if(grid_int[0] > 1100){
					gsac_control.xgrid_color = 1030;
					gsac_control.ygrid_color = 1030;
					} else{
					gsac_control.xgrid_color = grid_int[0];
					gsac_control.ygrid_color = grid_int[0];
					}
					break;
				case GRID_DFLT:
					gsac_control.xgrid = NO;
					gsac_control.ygrid = NO;
					gsac_control.xgrid_type = 2;
					gsac_control.ygrid_type = 2;
					gsac_control.xgrid_color = 1030;
					gsac_control.ygrid_color = 1030;
					gsac_control.xgrid_minor = YES;
					gsac_control.ygrid_minor = NO;
					break;
				case GRID_ON:
					gsac_control.xgrid = YES;
					gsac_control.ygrid = YES;
					break;
				case GRID_OFF:
					gsac_control.xgrid = NO;
					gsac_control.ygrid = NO;
					break;
				case GRID_SOLID:
					gsac_control.xgrid_type = 1;
					gsac_control.ygrid_type = 1;
					break;
				case GRID_DOTTED:
					gsac_control.xgrid_type = 2;
					gsac_control.ygrid_type = 2;
					break;
				case GRID_MINOR:
printf("grid_yn %d\n",grid_yn);
					if(grid_yn == NO){
						gsac_control.xgrid_minor = NO;
						gsac_control.ygrid_minor = NO;
					} if(grid_yn == YES){
						gsac_control.xgrid_minor = YES;
						gsac_control.ygrid_minor = YES;
					}
					break;
			}
		}
	}
			
		
}

void gsac_exec_grid(void)
{
	/* this just sets paramters */
}
