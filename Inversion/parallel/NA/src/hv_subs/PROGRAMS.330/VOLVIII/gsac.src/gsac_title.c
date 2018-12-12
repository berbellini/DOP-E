/*
	Changes:
	02 APR 2009 - to get this to work with the getargs I had to have a perfect comparison, e.g.,
		use TEXT text fully instead of T text
		This is a problem is getargs?
*/
#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	TITLE_DFLT	0
#define	TITLE_ON	1
#define	TITLE_OFF	2
#define	TITLE_LOCATION	3
#define	TITLE_SIZE	4
#define	TITLE_TEXT	5

#define TITLE_LOC_TOP    0
#define TITLE_LOC_BOTTOM 1
#define TITLE_LOC_LEFT   2
#define TITLE_LOC_RIGHT  3

#define TITLE_SIZE_TINY   0
#define TITLE_SIZE_SMALL  1
#define TITLE_SIZE_MEDIUM 2
#define TITLE_SIZE_LARGE  3



struct arghdr titlearg[] = {
	{TITLE_DFLT, "DEFAULT"    , IHDR, NO, 0, NO, "",  2},
	{TITLE_ON  , "ON"         , YHDR, NO, 0, NO, "ON ",  2},
	{TITLE_OFF  , "OFF"       , YHDR, NO, 0, NO, "OFF",  3},
	{TITLE_LOCATION,"LOCATION", CHDR, NO, 1, NO, "Location Top|Bottom|Left|Right ", 1},
	{TITLE_SIZE, "SIZE"       , CHDR, NO, 1, NO, "Size Tiny|Small|Medium|Large ", 1},
	{TITLE_TEXT, "TEXT"       , CHDR, NO, 1, NO, "TEXT text ", 4},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float title_real[10];
int   title_int [10];
int   title_yn;
int   title_num;

/* these are prototypes for global variables to be used by the routine */
int title_on ;
int title_size ;
int title_location;
char title_text[121] = "None";

void gsac_set_param_title(int ncmd, char **cmdstr)
{
	int i;
	char instr[121];
        char instr_size[80];
        char instr_location[80];
	/* initial debug */
	/*
	*/
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, titlearg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; titlearg[i].key[0] != '\0' ; i++){
		if(titlearg[i].used > 0){
			if(titlearg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, titlearg[i].key, 
					titlearg[i].mfit,titlearg[i].narg, title_real);
			} else if(titlearg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, titlearg[i].key, 
					titlearg[i].mfit,titlearg[i].narg, title_int );
			} else if(titlearg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, titlearg[i].key, 
					titlearg[i].mfit,titlearg[i].narg, &title_yn );
			} else if(titlearg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, titlearg[i].key, 
					titlearg[i].mfit,titlearg[i].narg, &title_num );
			} else if(titlearg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, titlearg[i].key, 
					titlearg[i].mfit, titlearg[i].narg, instr );
printf("%s\n",instr);
			}
			switch(titlearg[i].id){
				case TITLE_SIZE:
					switch(instr[0]){
					   case 'T':
					   case 't':
					      title_size = TITLE_SIZE_TINY;
					      break;
					   case 'S':
					   case 's':
					      title_size = TITLE_SIZE_SMALL;
					      break;
					   case 'M':
					   case 'm':
					      title_size = TITLE_SIZE_MEDIUM;
					      break;
					   case 'L':
					   case 'l':
					      title_size = TITLE_SIZE_LARGE;
					      break;
					   default:
					      title_size = TITLE_SIZE_SMALL;
					}
					break;
				case TITLE_LOCATION:
					switch(instr[0]){
					   case 'T':
					   case 't':
					      title_location = TITLE_LOC_TOP;
					      break;
					   case 'r':
					   case 'R':
					      title_location = TITLE_LOC_RIGHT;
					      break;
					   case 'B':
					   case 'b':
					      title_location = TITLE_LOC_BOTTOM;
					      break;
					   case 'L':
					   case 'l':
					      title_location = TITLE_LOC_LEFT;
					      break;
					   default:
					      title_location = TITLE_LOC_TOP;
					}
					break;
				case TITLE_TEXT:
					strcpy(title_text,instr);
printf("%s\n",instr);
printf("%s\n",title_text);
					break;
				case TITLE_ON:
					title_on = YES ;
					break;
				case TITLE_OFF:
					title_on = NO ;
					break;
				case TITLE_DFLT:
					title_on = NO ;
					title_location = TITLE_LOC_TOP;
					title_size = TITLE_SIZE_SMALL;
					break;

			}
		}
	}
			
		
}

void gsac_exec_title(void)
{
	printf("TITLE_ON        : %d\n",title_on);
	printf("TITLE_LOCATION  : %d\n",title_location);
	printf("TITLE_SIZE      : %d\n",title_size);
	printf("TITLE_TEXT      : %s\n",title_text);
}
