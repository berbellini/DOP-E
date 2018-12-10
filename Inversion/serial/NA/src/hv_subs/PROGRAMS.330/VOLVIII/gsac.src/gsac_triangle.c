#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

/* Changes:
	14 APR 2007 - ensure that if HALFWIDTH < DT, nothing is done */

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	TRIANGLE_DFLT	0
#define	TRIANGLE_WIDTH		1
#define	TRIANGLE_HALFWIDTH	2


struct arghdr trianglearg[] = {
	{TRIANGLE_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{TRIANGLE_WIDTH  , "WIDTH"  , RHDR, NO, 1, NO, "Width width", 1},
	{TRIANGLE_HALFWIDTH  , "Half"  , RHDR, NO, 1, NO, "Half half_width", 1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
static float triangle_real[10];
static int   triangle_int [10];
static float   triangle_halfwidth;
static int   triangle_yn;
static int   triangle_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_triangle(int ncmd, char **cmdstr)
{
	int i;
	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, trianglearg, NO, YES))
		return;
	/* parse commands */
	triangle_halfwidth = -1. ;
	for(i=0 ; trianglearg[i].key[0] != '\0' ; i++){
		if(trianglearg[i].used > 0){
			if(trianglearg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, trianglearg[i].key, 
					trianglearg[i].mfit,trianglearg[i].narg, triangle_real);
			} else if(trianglearg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, trianglearg[i].key, 
					trianglearg[i].mfit,trianglearg[i].narg, triangle_int );
			} else if(trianglearg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, trianglearg[i].key, 
					trianglearg[i].mfit,trianglearg[i].narg, &triangle_yn );
			} else if(trianglearg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, trianglearg[i].key, 
					trianglearg[i].mfit,trianglearg[i].narg, &triangle_num );
			}
			switch(trianglearg[i].id){
				case TRIANGLE_WIDTH:
					triangle_halfwidth = triangle_real[0]/2.0;
					break;
				case TRIANGLE_HALFWIDTH:
					triangle_halfwidth = triangle_real[0];
					break;

			}
		}
	}
			
		
}

void gsac_exec_triangle(void)
{
	if(triangle_halfwidth > 0.0)
		gsac_pulconv(triangle_halfwidth,0.0,triangle_halfwidth);
}
