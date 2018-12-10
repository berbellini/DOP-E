#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	"gsac.h"
#include	"gsac_command.h"
#include	"gsac_docommand.h"

/* this is actually getting very lengthy 
 * do not forget to pass argument strings
*/

void (*do_set_param());
void (*do_exec());

/* clean this up later */
char tmpstr[1000];

int parsecommand(char *cmdstr)
{
	int i;
	strcpy(tmpstr,cmdstr);
	gsac_strupr(tmpstr);
	for (i=0;i< sizeof(gsac_command_list)/sizeof(struct gsac_command);i++){
		if(strcmp(tmpstr,gsac_command_list[i].str) == 0){
			return(i);
		}
	}
	/* failure */
	return(-1);
}

int gsac_docommand(int ncmd, char  **cmdstr, char *input_lineptr)
{
	int cmd, cmdh, iret;
	cmd = parsecommand(cmdstr[0]);
	/* this is not clean too if then else */
	if(cmd < 0){
		/* pass this off to the system */
		printf("%s\n",input_lineptr);
		iret = system(input_lineptr);
		return(0);
	} else {
		if(gsac_command_list[cmd].cmd == EXIT) {
			return (1);
		} else if (gsac_command_list[cmd].cmd == HELP){
			if(ncmd == 1  ) {
				cmdh = 0;
				gsac_help(ncmd, cmdstr, gsac_command_list[cmdh].helper);
			} else {
				cmdh = parsecommand(cmdstr[1]);
				if(cmdh < 0  ){
					printf("%s is not a keyword\n",cmdstr[1]);
				} else {
				gsac_help(ncmd, cmdstr, gsac_command_list[cmdh].helper);
				}
			}
			return (0);

		} else {
			gsac_command_list[cmd].gsac_set_param(ncmd,cmdstr);
			gsac_command_list[cmd].gsac_exec();
			return (0);
		}
	}
}

