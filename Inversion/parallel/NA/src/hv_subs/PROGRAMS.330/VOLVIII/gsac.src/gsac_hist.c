#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	HIST_DFLT	0
#define	HIST_LIST	1

static int hist_count = -1;


struct arghdr histarg[] = {
	{HIST_DFLT, "DEFAULT", IHDR, 0, 0, NO, "", 1},
	{HIST_LIST , "LIST"  , IHDR, 0, 1, NO, "LIST n", 1},
	{0,	""		, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float hist_real[10];
int   hist_int [10];
int   hist_yn;
int   hist_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_hist(int ncmd, char **cmdstr)
{
	int i;
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, histarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; histarg[i].key[0] != '\0' ; i++){
		if(histarg[i].used > 0){
			if(histarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, histarg[i].key, 
					histarg[i].mfit,histarg[i].narg, hist_real);
			} else if(histarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, histarg[i].key, 
					histarg[i].mfit,histarg[i].narg, hist_int );
			} else if(histarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, histarg[i].key, 
					histarg[i].mfit,histarg[i].narg, &hist_yn );
			} else if(histarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, histarg[i].key, 
					histarg[i].mfit,histarg[i].narg, &hist_num );
			}
			switch(histarg[i].id){
				case HIST_DFLT:
					hist_count = -1;
					break;
				case HIST_LIST:
					hist_count = hist_int[0];
					if(hist_count < 0 )
						hist_count = -1 ;
					break;

			}
		}
	}
			
		
}

void gsac_exec_hist(void)
{
#ifdef READLINE_LIBRARY
#  include "readline.h"
#  include "history.h"
	int i,is;
	int list_count;
	extern HIST_ENTRY **history_list ();
	HIST_ENTRY **list;
	list = history_list ();
	/* get the number of entries in the list */
	if (list){
		list_count = 0;
		for (i = 0; list[i]; i++)
			list_count++;
	}
	if(hist_count > 0){
		is = list_count - hist_count ;
		is = MAX(0, is);
	} else {
		is = 0;
	}
	if (list){
		for (i = is; list[i] && i < list_count ; i++)
			printf("%s\n",list[i]->line);
	}
#endif
}
