#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

void gsac_set_param_qdp(int ncmd, char **cmdstr)
{
	int i;
	if(ncmd == 1)
			return;
	/* look for the magin key work of ON OFF or and integer n
	 * */

	gsac_strupr(cmdstr[1]);
	if(strcmp(cmdstr[1],"ON")==0)
		gsac_control.qdp = 0;
	else if(strcmp(cmdstr[1],"OFF")==0)
		gsac_control.qdp = -1;
	else {
		if(isargi(cmdstr[1],&i) == YES)
			if(i > 0)
				gsac_control.qdp = i;
	}
		
}

void gsac_exec_qdp(void)
{
	/* nothing to do here 
	 * we have already defined the 
	 * gac_control.qdp parameter which is used by the plot routines
	 * PLOT1 and PLOTPK but not by the PLOTSP (yet?)  */
}
