#include	<stdio.h>
#include	<string.h>
#include	"gsac.h"

/* this is not clean
 * gsac_command is defined in gsac_command.h but the command list in initialized there too
*/

void gsac_help(int ncmd, char **cmdstr, char **helper)
{
	int i;
		/* eventually put a help help here */
	/*
	if(ncmd == 1){
		return;
	}
	*/
		for(i = 0 ; strlen(helper[i]) > 0 ; i++ )
			printf("%s",helper[i]);
}


