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

#define MATH_SQR    1
#define MATH_SQRT   2
#define MATH_EXP    3
#define MATH_LOG    4
#define MATH_EXP10  5
#define MATH_LOG10  6
#define MATH_ABS    7
static int math_which;


void gsac_set_param_math(int ncmd, char **cmdstr)
{
	/* convert command to upper case */
	gsac_strupr(cmdstr[0]);
	/* define the command - this is ugly but we only have 4 functions */
	if(strcmp(cmdstr[0],"SQR") == 0)
		math_which = MATH_SQR;
	else if(strcmp(cmdstr[0],"SQRT") == 0)
		math_which = MATH_SQRT;
	else if(strcmp(cmdstr[0],"EXP") == 0)
		math_which = MATH_EXP;
	else if(strcmp(cmdstr[0],"EXP10") == 0)
		math_which = MATH_EXP10;
	else if(strcmp(cmdstr[0],"LOG") == 0)
		math_which = MATH_LOG;
	else if(strcmp(cmdstr[0],"LOG10") == 0)
		math_which = MATH_LOG10;
	else if(strcmp(cmdstr[0],"ABS") == 0)
		math_which = MATH_ABS;
	else
		math_which = -1;		/* this should be impossible */

}

void gsac_exec_math(void)
{
	int i, k, npts, ntrc;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float x;
	int doproc;

	/* only apply the filtering if permitted */

	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;

	/* process the traces */
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		if(npts > 0){
			depmax = sacdata[k].sachdr.rhdr[H_DEPMAX];
			depmin = sacdata[k].sachdr.rhdr[H_DEPMIN];
			/* SANITY CHECK */
			doproc = YES;
			switch(math_which){
				case MATH_SQRT:
					if(depmin < 0.0 || depmax < 0.0){
						printf("Cannot take SQRT of negative numbers in %s\n",sacdata[k].sac_ifile_name);
						doproc = NO;
					}
					break;
				case MATH_LOG:
				case MATH_LOG10:
					if(depmin <= 0.0 || depmax <= 0.0){
						printf("Cannot take LOG or LOG10 of negative or zero numbers in %s\n",sacdata[k].sac_ifile_name);
						doproc = NO;
					}
					break;
				case MATH_EXP:
					if(depmax > 85){
						printf("Cannot take EXP of trace values > 85 in %s\n",sacdata[k].sac_ifile_name);
						doproc = NO;
					}
					break;
				case MATH_EXP10:
					if(depmax > 85){
						printf("Cannot take EXP10 of trace values > 37 in %s\n",sacdata[k].sac_ifile_name);
						doproc = NO;
					}
					break;
			}
			if(doproc){
			for(i=0; i < npts ; i++){
				x = sacdata[k].sac_data[i];
				/* this is not efficient */
				switch(math_which){
					case MATH_ABS:
						sacdata[k].sac_data[i] = ABS(x);
						break;
					case MATH_SQR:
						sacdata[k].sac_data[i] = x*x;
						break;
					case MATH_SQRT:
						sacdata[k].sac_data[i] = sqrt(ABS(x));
						break;
					case MATH_EXP:
						sacdata[k].sac_data[i] = exp(x);
						break;
					case MATH_EXP10:
						sacdata[k].sac_data[i] = 2.30259*exp(x);
						break;
					case MATH_LOG:
						if(x <= 0)
							x = 0.00001 * depmax;
						sacdata[k].sac_data[i] = log(x);
						break;
					case MATH_LOG10:
						if(x <= 0)
							x = 0.00001 * depmax;
						sacdata[k].sac_data[i] = log10(x);
						break;
				}
			}
			}
			/* get the DFT */
			/* redetermine the depmax depmin depmen */
			getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		}
	}
}
