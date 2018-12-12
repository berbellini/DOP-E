#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	ECHO_DFLT	0
#define	ECHO_X0		1
#define	ECHO_Y0		2
#define	ECHO_XLEN	3
#define	ECHO_YLEN	4
#define	ECHO_XLAB	5
#define	ECHO_YLAB	6


struct arghdr echoarg[] = {
	{ECHO_DFLT, "DEFAULT", IHDR, 0, 0, NO, "", 1},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};

/* these are temporary variables only used here */
float echo_real[10];
int   echo_int [10];
int   echo_yn;
int   echo_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_echo(int ncmd, char **cmdstr)
{
	int i;
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
}

void gsac_exec_echo(void)
{
}
