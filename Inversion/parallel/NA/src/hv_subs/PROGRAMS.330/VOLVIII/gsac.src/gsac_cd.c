#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

#include <unistd.h>
#include <errno.h>

#ifdef WIN32
#include	"glob.h"
#else
#include	<glob.h>
#endif
#include	"gsac.h"

/* safety for old cc compiler suite on SOLARIS */
#ifndef GLOB_TILDE
#define  GLOB_TILDE      0x0000  /* Expand tilde names from the passwd file. */
#endif

extern int errno;
#define BSIZ 100
char obuf[BSIZ + 1];
glob_t globdbuf;



/* special include for change directory */
#include <unistd.h>

extern struct sacfile_ *sacdata;
extern int *sortptr;


void gsac_set_param_cd(int ncmd, char **cmdstr)
{
	if(ncmd == 1)
		return;
	glob(cmdstr[1],  GLOB_TILDE, NULL, &globdbuf);
	if(globdbuf.gl_pathc > 0){
		if(chdir(globdbuf.gl_pathv[0]) < 0){
			printf("errno = %d\n",errno);
			perror(cmdstr[1]);
		}
	} else {
		printf("Cannot expand %s\n",cmdstr[1]);
	}
	printf("Current directory is %s\n",getcwd(obuf,BSIZ));
			
		
}

void gsac_exec_cd(void)
{
}
