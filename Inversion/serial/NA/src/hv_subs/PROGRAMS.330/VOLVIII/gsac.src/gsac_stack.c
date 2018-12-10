#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include	"csstim.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

static float *xa = (float *)NULL;

#define	STACK_DFLT	0
#define	STACK_ABS	1
#define	STACK_REL	2
#define	STACK_SUF	3
#define STACK_NORM	4

static int   stack_absolute = YES ;
static char  stack_suffix[80];
static int   stack_dosuffix = NO ;


struct arghdr stackarg[] = {
	{STACK_DFLT, "DEFAULT", IHDR, 0, 0, NO, "", -1},
	{STACK_ABS , "ABSOLUTE"  , IHDR, 0, 0, NO, "Absolute", 1},
	{STACK_REL , "RELATIVE"  , IHDR, 0, 0, NO, "Relative", 1},
	{STACK_SUF , "SUFFIX"    , CHDR, 0, 1, YES, "",1},
	{STACK_NORM,  "NORM" , YHDR, 0, 1, NO, "NORM [ON|OFF] ", 1},
	{0,	""		, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float stack_real[10];
int   stack_int [10];
int   stack_yn;
int   stack_num;

static int stacknorm = NO ;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_stack(int ncmd, char **cmdstr)
{
	int i;
	char instr[80];
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, stackarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; stackarg[i].key[0] != '\0' ; i++){
		if(stackarg[i].used > 0){
			if(stackarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, stackarg[i].key, 
					stackarg[i].mfit,stackarg[i].narg, stack_real);
			} else if(stackarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, stackarg[i].key, 
					stackarg[i].mfit,stackarg[i].narg, stack_int );
			} else if(stackarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, stackarg[i].key, 
					stackarg[i].mfit,stackarg[i].narg, &stack_yn );
			} else if(stackarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, stackarg[i].key, 
					stackarg[i].mfit,stackarg[i].narg, &stack_num );
			} else if(stackarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, stackarg[i].key, 
					stackarg[i].mfit, stackarg[i].narg, instr );
			}
			switch(stackarg[i].id){
				case STACK_DFLT:
					stack_absolute = YES;
					stack_dosuffix = NO;
				case STACK_REL:
					stack_absolute = NO;
					break;
				case STACK_ABS:
					stack_absolute = YES;
					break;
				case STACK_NORM:
					if(stack_yn == NO)
						stacknorm = NO;
					else if(stack_yn == YES)
						stacknorm = YES;
					break;
				case STACK_SUF:
					if(strlen(instr) < 80){
					strcat(stack_suffix,instr);
						stack_dosuffix = YES;
					}
					break;
			}
		}
	}
			
		
}

void gsac_exec_stack(void)
{
	double tzbeg, tzend, tzref;
	double ts, te;
	int i, j, k, ntrc, npts, newpts;
	int xl;
	float dtx, dty;
	int stack_master;
	float delta;
	float depmax, depmin, depmen;
	int indmax, indmin;
	int month, day;

	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1){
		return;
	}

	/* CHECK TO SEE THAT ALL DT ARE THE SAME ELSE RETURN */
	/* define the master trace */
	stack_master = 0;
	dtx = sacdata[stack_master].sachdr.rhdr[H_DELTA];
	for ( k=0 ; k < ntrc ; k ++){
		dty = sacdata[k].sachdr.rhdr[H_DELTA];
		if(ABS (dtx - dty) > 0.01 * ABS(dtx)){
			printf("DT not equal - stacking not done\n");
			printf("Files %d %d with %f %f, respectively\n",stack_master,k,dtx,dty);
			return;
		}
	}
	/* now look for the proper time window */
	delta = dtx;
	if(stack_absolute == YES){
		for(k=0; k < ntrc ; k++){
			if(k == 0 ){
				ts = sacdata[k].tzbeg;
				te = sacdata[k].tzend;
				tzref = sacdata[k].tzref;
			} else {
				tzbeg = sacdata[k].tzbeg;
				tzend = sacdata[k].tzend;
				ts = MAX(ts,tzbeg);
				te = MIN(te,tzend);
			}
		}
		newpts = 1 + (int)((te - ts)/delta + 0.49);
	} else {
		for(k=0; k < ntrc ; k++){
			npts = sacdata[k].sachdr.ihdr[H_NPTS];
			if(k == 0 ){
				newpts = npts;
				ts = sacdata[k].tzbeg;
				te = sacdata[k].tzend;
				tzref = sacdata[k].tzref;
			} else {
				newpts = MIN(newpts,npts);
			}
		}
	}
	printf("New time series length:  %d\n",newpts);
	/* create the new time series using calloc to ensure that it is
	 * full of zeros */
	if(xa == (float *)NULL)
		xa = (float *)calloc(newpts,sizeof(float));
	else
		xa = (float *)realloc(xa,newpts*sizeof(float));
	/* force initialization */
	for(i=0 ; i < newpts ; i++)
		xa[i] = 0.0;
	/* now stack on the basis of the ts te */
	
	/* now redefine the parameters of the first trace */
	if(stack_absolute){
		xl = 0;
		for(k=0; k < ntrc ; k++){
			tzbeg = sacdata[k].tzbeg;
			npts = sacdata[k].sachdr.ihdr[H_NPTS];
			j = (int)((ts - tzbeg )/dtx + 0.49);
			for(i = 0 ; i < newpts ; i++){
				xa[i] += sacdata[k].sac_data[i + j];
			}
		}
	} else {
		/* relative stacking */
		for(k=0; k < ntrc ; k++){
			for(i=0 ; i < newpts; i++){
				xa[i] += sacdata[k].sac_data[i];
			}
		}
	}

	sacdata[0].tzref = tzref;
	sacdata[0].sachdr.ihdr[H_NPTS] = newpts;
	sacdata[0].sac_data = realloc(sacdata[0].sac_data,newpts*sizeof(float));
	sacdata[0].sachdr.rhdr[H_B] = ts - tzref;
	sacdata[0].sachdr.rhdr[H_E] = ts - tzref + (newpts -1 ) * delta;
	/* copy array in to save */
	if(stacknorm == YES){
		printf("Stack is divided by number of traces \n");
		for(i=0 ; i < newpts ; i++){
			sacdata[0].sac_data[i] = xa[i]/(float)ntrc;
		}
	} else {
		for(i=0 ; i < newpts ; i++){
			sacdata[0].sac_data[i] = xa[i];
		}
	}
	getmxmn(sacdata[0].sac_data, newpts,&depmax, &depmin, &depmen,&indmax,&indmin);
	sacdata[0].sachdr.rhdr[H_TIMMAX] = sacdata[0].sachdr.rhdr[H_B]  + ( indmax)*sacdata[0].sachdr.rhdr[H_DELTA] ;
	sacdata[0].sachdr.rhdr[H_TIMMIN] = sacdata[0].sachdr.rhdr[H_B]  + ( indmin)*sacdata[0].sachdr.rhdr[H_DELTA] ;
	sacdata[0].sachdr.rhdr[H_DEPMIN] = depmin;
	sacdata[0].sachdr.rhdr[H_DEPMEN] = depmen;
	sacdata[0].sachdr.rhdr[H_DEPMAX] = depmax;

	sacdata[0].tzbeg = ts;
	sacdata[0].tzend = te;
	sacdata[0].tzbegx = sacdata[0].tzbeg;
	sacdata[0].tzendx = sacdata[0].tzend;
	/* get bounds for absolute plotting */
	gsac_control.begmin = sacdata[0].tzbeg;
	gsac_control.endmax = sacdata[0].tzend;

	/* NOW RESET THE OTHER TIME FIELDS */
	etoh(tzref, &sacdata[0].sachdr.ihdr[H_NZYEAR], 
		&sacdata[0].sachdr.ihdr[H_NZJDAY], &month, &day,
		&sacdata[0].sachdr.ihdr[H_NZHOUR], 
		&sacdata[0].sachdr.ihdr[H_NZMIN],
		&sacdata[0].sachdr.ihdr[H_NZSEC], 
		&sacdata[0].sachdr.ihdr[H_NZMSEC]);

		if(stack_dosuffix == YES){
			strcpy(sacdata[0].sac_ofile_name,sacdata[0].sac_ifile_name);
			strcat(sacdata[0].sac_ofile_name,stack_suffix);
		} else {
			strcpy(sacdata[0].sac_ofile_name,sacdata[0].sac_ifile_name);
			strcat(sacdata[0].sac_ofile_name,".stk");
		}
	/* redefine the number of traces to 1 */
		gsac_control.number_otraces = 1 ;
	/* save the number of traces actually stacked */
		sacdata[0].sachdr.ihdr[H_IHDR11] =  ntrc;
	printf("New default output filename: %s\n",sacdata[0].sac_ofile_name);


}
