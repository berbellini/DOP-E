#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

#define	CE_ON		0
#define CE_OFF		1
#define CE_FILLZ	2

struct arghdr cearg[] = {

	{CE_ON, "ON"		, IHDR, 0, 0, NO, "ON",2},
	{CE_OFF, "OFF"		, IHDR, 0, 0, NO, "OFF",2},
	{CE_FILLZ, "FILLZ"	, IHDR, 0, 0, NO, "FILLZ ", 1},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};


void gsac_set_param_cuterr(int ncmd, char **cmdstr)
{
	int i;
	/* note when the testrg routine is used, if the argument is
		NO then you must use internal variables to define the 
		state of the operation - if you use YES, then things are
		not changed until the input is proven correct. An exmple of
		this concept with YES is the following:
		Assume we wish aa LP filter with fc 1 np 2 p 1 
		If we enter  fc 2 np2   there is a syntax error and we
		should not chnge the fc since the np2 is wrong. One way to
		do this in the code would be to do two calls

			if(testarc,ncmd, cmdstr, cmdargs, YES) is OK
			then
				testarc,ncmd, cmdstr, cmdargs, NO)
		*/
	if(testarg(ncmd, cmdstr, cearg, NO, YES))
		return  ;
	for(i=0 ; cearg[i].key[0] != '\0' ; i++){
		if(cearg[i].used > 0){
			switch(cearg[i].id){
				case CE_ON:
					gsac_control.docut = YES;
					break;
				case CE_OFF:
					gsac_control.docut = NO;
					break;
				case CE_FILLZ:
					break;
			}
		}
	}
			
		
}

void gsac_exec_cuterr(void)
{
}
