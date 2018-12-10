#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

#define FMT_PRESERVE	0
#define FMT_LOCAL	1

struct arghdr fmtarg[] = {
	{FMT_PRESERVE,"PRESERVE", NHDR, 0, 0, NO, "PRESERVE", 1},
	{FMT_LOCAL   ,"LOCAL"   , NHDR, 0, 0, NO, "LOCAL", 1}
};


void gsac_set_param_format(int ncmd, char **cmdstr)
{
	int i;
	if(ncmd == 1 )
		return;
	if(testarg(ncmd, cmdstr, fmtarg, NO, YES))
		return   ;
	for(i=0 ; fmtarg[i].key[0] != '\0' ; i++){
		if(fmtarg[i].used > 0){
			if(fmtarg[i].ricell == NHDR){
				switch(fmtarg[i].id){
					case FMT_PRESERVE:
						gsac_control.local = NO ;
						break;
					case FMT_LOCAL:
						gsac_control.local = YES ;
						break;
				}
			}
		}
	}
}

void gsac_exec_format(void)
{
}
