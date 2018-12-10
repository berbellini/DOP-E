#include        <stdio.h>
#include        <string.h>
#include "gsac.h"
#include "gsac_docommand.h"
#include "gsac_sac.h"
#include "gsac_sachdr.h"
#include "gsac_arg.h"
#include "csstim.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


void gsac_set_param_wh(int ncmd, char **cmdstr)
{
	int i;
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
}

void gsac_exec_wh(void)
{
	int k, ntrchdr;
	/* safety - do not write header if within a cut - force user to
	 * overwrite */
	if(gsac_control.docut){
		printf("Cannot write headers while CUT is on\n");
		printf("Note that a WRITE will replace original file\n");
	} else {
		ntrchdr = gsac_control.number_iheaders;
		if(ntrchdr < 1)
			return;
		for ( k=0 ; k < ntrchdr ; k ++){
			bwsach(sacdata[k].sac_ofile_name,sacdata[k].sachdr);
		}
	}
}
