#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

#define	PCTL_DFLT	0
#define	PCTL_X0		1
#define	PCTL_Y0		2
#define	PCTL_XLEN	3
#define	PCTL_YLEN	4
#define	PCTL_XLAB	5
#define	PCTL_YLAB	6
#define	PCTL_GRID	7


struct arghdr pctlarg[] = {
	{PCTL_DFLT, "DEFAULT", IHDR, 0, 0, NO, "", 1},
	{PCTL_X0  , "X0"  , RHDR, 0, 1, NO, "X0 x0",2},
	{PCTL_Y0  , "Y0"  , RHDR, 0, 1, NO, "Y0 y0",2},
	{PCTL_XLEN, "XLEN", RHDR, 0, 1, NO, "XLEN xlen",3},
	{PCTL_YLEN, "YLEN", RHDR, 0, 1, NO, "YLEN ylen",3},
	{PCTL_XLAB, "XLAB", CHDR, 0, 1, NO, "XLAB xlabel",3},
	{PCTL_YLAB, "YLAB", CHDR, 0, 1, NO, "YLAB ylabel",3},
	{PCTL_GRID, "GRID", YHDR, 0, 1, NO, "GRID [ON|OFF]", 1},
	{0,	""	  , IHDR, 0, 0, NO, "",-1}
};

/* these are temporary variables only used here */
float pctl_real[10];
int   pctl_int [10];
int   pctl_yn;
int   pctl_num;

void gsac_set_param_pctl(int ncmd, char **cmdstr)
{
	int i;
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, pctlarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; pctlarg[i].key[0] != '\0' ; i++){
		if(pctlarg[i].used > 0){
			if(pctlarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, pctlarg[i].key, 
					pctlarg[i].mfit,pctlarg[i].narg, pctl_real);
			} else if(pctlarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, pctlarg[i].key, 
					pctlarg[i].mfit,pctlarg[i].narg, pctl_int );
			} else if(pctlarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, pctlarg[i].key, 
					pctlarg[i].mfit,pctlarg[i].narg, &pctl_yn );
			} else if(pctlarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, pctlarg[i].key, 
					pctlarg[i].mfit,pctlarg[i].narg, &pctl_num );
			}
			switch(pctlarg[i].id){
				case PCTL_X0:
					gsac_control.x0 = pctl_real[0];
					break;
				case PCTL_Y0:
					gsac_control.y0 = pctl_real[0];
					break;
				case PCTL_XLEN:
					gsac_control.xlen = pctl_real[0];
					break;
				case PCTL_YLEN:
					gsac_control.ylen = pctl_real[0];
					break;
				case PCTL_XLAB:
					break;
				case PCTL_YLAB:
					break;
				case PCTL_GRID:
					gsac_control.grid = pctl_yn ;
					break;
				case PCTL_DFLT:
					gsac_control.x0 = 1.25 ;
					gsac_control.y0 = 1.0 ;
					gsac_control.xlen = 8.0 ;
					gsac_control.ylen = 6.0 ;
					break;

			}
		}
	}
}

void gsac_exec_pctl(void)
{
	/* set up plot controls */
	printf("X0 %f Y0 %f XLEN %f YLEN %f\n",gsac_control.x0,
			gsac_control.y0, gsac_control.xlen, gsac_control.ylen);
}
