#include	<stdio.h>
#include	"gsac.h"
#include	"gsac_sac.h"
#include	"gsac_docommand.h"

extern struct sacfile_ *sacdata;

void gsac_set_param_nop(int ncmd, char **cmdstr)
{
	int i;
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
}

void gsac_exec_nop(void)
{
}
