/* CHANGES
 * 2004 09 24 - this merges the EXP SQR SQRT and LOG routines
 * 	The nature of the command is determined by looking at
 * 	the command name, e.g., cmdstr[0]
 *
 * */
#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

void gsac_set_param_plotctl(int ncmd, char **cmdstr)
{
	/* convert command to upper case */
	gsac_strupr(cmdstr[0]);
	/* define the command - this is ugly but we only have 4 functions */
	if(strcmp(cmdstr[0],"XLIN") == 0){
		gsac_control.plotlinx = YES;
	} else if(strcmp(cmdstr[0],"YLIN") == 0){
		gsac_control.plotliny = YES;
	} else if(strcmp(cmdstr[0],"XLOG") == 0){
		gsac_control.plotlinx = NO ;
	} else if(strcmp(cmdstr[0],"YLOG") == 0){
		gsac_control.plotliny = NO ;
	} else if(strcmp(cmdstr[0],"LINLIN") == 0){
		gsac_control.plotlinx = YES;
		gsac_control.plotliny = YES;
	} else if(strcmp(cmdstr[0],"LINLOG") == 0){
		gsac_control.plotlinx = YES;
		gsac_control.plotliny = NO ;
	} else if(strcmp(cmdstr[0],"LOGLIN") == 0){
		gsac_control.plotlinx = NO ;
		gsac_control.plotliny = YES;
	} else if(strcmp(cmdstr[0],"LOGLOG") == 0){
		gsac_control.plotlinx = NO ;
		gsac_control.plotliny = NO ;
	}
}

void gsac_exec_plotctl(void)
{
if(gsac_control.plotlinx)
	printf("xlin-");
else
	printf("xlog-");
if(gsac_control.plotliny)
	printf("ylin \n");
else
	printf("ylog \n");
}
