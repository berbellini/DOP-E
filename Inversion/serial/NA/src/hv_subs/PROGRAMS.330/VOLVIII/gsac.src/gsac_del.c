#include        <stdio.h>
#include        "gsac.h"
#include        "gsac_sac.h"
#include	"gsac_docommand.h"
#include	"gsac_arg.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

static int deletefile = -1;




void gsac_set_param_del(int ncmd, char **cmdstr)
{
	int tmp;
	if(ncmd == 1)
		return;
	/* for this command the second argument must be an integer */
	if(isargi(cmdstr[1], &tmp)==YES){
		if(tmp < 0 || tmp > gsac_control.max_number_traces -1){
			printf("DELETE trace_number - incorrect trace_number\n");
		} else {
			deletefile = tmp;
		}
	}
	
}

void gsac_exec_del(void)
{
	int k;
	k = deletefile;
	sacdata[k].display = NO;
}
