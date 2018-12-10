/* 
	Changes:
	10 JAN 2005 - all previous traces are ignored
		only one new trace results from this operation

	27 AUG 2005 - set gsac_control.fmax correctly
*/
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


#define FG_IMPULSE	0
#define FG_DELTA	1
#define FG_NPTS		2
#define FG_TRIANGLE	3
#define FG_BOX		4
#define FG_LENGTH	5
#define FG_COMB		6


struct arghdr fgarg[] = {
	{FG_IMPULSE, "IMPULSE"	, XHDR, 0, 0, NO, "IMPULSE", 1},
	{FG_TRIANGLE, "TRIANGLE", XHDR, 0, 0, NO, "TRIANGLE", 1},
	{FG_BOX, "BOX", XHDR, 0, 0, NO, "BOX", 1},
	{FG_DELTA, "DELTA"	, RHDR, 0, 1, NO, "DELTA dt", 1},
	{FG_NPTS, "NPTS"	, IHDR, 0, 1, NO, "NPTS npts ", 1},
	{FG_LENGTH, "LENGTH"	, RHDR, 0, 1, NO, "LENGTH length", 1},
	{FG_COMB, "COMB"	, RHDR, 0, 2, NO, "COMB ncomb delay", 1},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};

/* these are temporary variables only used here */
float fg_real[10];
int   fg_int [10];

static float fg_delta = 1.0;
static float fg_length = 2.0;
static int fg_npts = 1024;
static int do_funcgen;
static int fg_type = FG_IMPULSE;
static float comb_num = 1;
static int comb_n = 1;
static float comb_delay = 0.0;
static int comb_ndelay = 0;



void gsac_set_param_fg(int ncmd, char **cmdstr)
{
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
	int i;

	do_funcgen = NO;
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, fgarg, NO, YES))
	       	return	;
	do_funcgen = YES;
	for(i=0 ; fgarg[i].key[0] != '\0' ; i++){
		if(fgarg[i].used > 0){
			if(fgarg[i].ricell == RHDR)
				getargr(ncmd, cmdstr, fgarg[i].key, 
					fgarg[i].mfit, fgarg[i].narg, fg_real);
			else if(fgarg[i].ricell == IHDR)
				getargi(ncmd, cmdstr, fgarg[i].key, 
					fgarg[i].mfit, fgarg[i].narg, fg_int );
			switch(fgarg[i].id){
				case FG_NPTS:
					if(fg_int[0] < 1)
						fg_npts = 1024 ;
					else
						fg_npts = fg_int[0];
					break;
				case FG_DELTA:
					if(fg_int[0] < 0)
						fg_delta = 1;
					else
						fg_delta = fg_real[0];
					break;
				case FG_COMB:
					comb_num = fg_real[0];
					if(comb_num <= 1.0){
						comb_n = 1;
						comb_delay = 0.0;
					} else {
						comb_n = (int)comb_num;
						comb_delay = fg_real[1];
					}
					break;
				case FG_LENGTH:
					if(fg_int[0] < 0)
						fg_length = 1;
					else
						fg_length = fg_real[0];
					break;
				case FG_IMPULSE:
					fg_type = FG_IMPULSE;
					break;
				case FG_TRIANGLE:
					fg_type = FG_TRIANGLE;
					break;
				case FG_BOX:
					fg_type = FG_BOX;
					break;
			}
		}
	}
}

			
		

void gsac_exec_fg(void)
{
	int k, ntrc, i, j, nl, coff, pos;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float fg_norm;

	if(do_funcgen == NO){
		printf("FUNCGEN incorrect syntax: FUNCGEN IMPULSE DELTA delta NPTS npts\n");
		return;
	} else {
		printf("FUNCGEN: creating synthetic with DELTA = %f and NPTS = %d\n",fg_delta,fg_npts); 
		if(ntrc > 0)
		printf("       : previous traces in memory are removed\n");
	}
	/* create a new trace structure */
	gsac_alloc_trace(0);
	gsac_control.number_itraces = 1;
	/* now initialize headers */

	ntrc = gsac_control.number_itraces;
	k = ntrc  - 1;
	/* perhaps put this as a separate routine called newhdr */
	for(i=0 ; i < 70 ; i++)
		sacdata[k].sachdr.rhdr[i] = -12345. ;
	for(i=0 ; i < 40 ; i++)
		sacdata[k].sachdr.ihdr[i] = -12345 ;
	for(i=0;i<24;i++){
		for(j=0;j<8;j++){
			strncpy(sacdata[k].sachdr.chdr[i],"-12345  ",8);
			strncpy(sacdata[k].schdr[i],"-12345  ",8);
		}
		sacdata[k].schdr[i][8] = '\0';
	}
	/* the header must now have actual values put in 
	 * NVHDR  8 9 DELTA DEPMAX DEPMIN O B E NPTS IFTYPE IZTYPE */
	/* initilize data stream to 0.0 */
	sacdata[k].sac_data =  (float *)calloc(fg_npts,sizeof(float));
	/* put in the IMPULSE HERE */
	if(comb_delay > 0.0){
		comb_ndelay = comb_delay/fg_delta;
		fg_norm = fg_delta*comb_n;
	} else {
		comb_ndelay = 0;
		fg_norm = fg_delta;
	}
printf("comb_delay %f comb_ndelay %d comb_n %d fg_norm %f\n",comb_delay,comb_ndelay,comb_n, fg_norm);
	/* this is very simple, but the comb option forces us to
		ensure that array dimensions are not exceeded -
		so a little ugly */
	for(j=0 ; j < comb_n ; j++){
		coff= j*comb_ndelay - (comb_n -1)*comb_ndelay/2;
		if(fg_type == FG_IMPULSE){
			pos = coff+fg_npts/2;
			if(pos >= 0 && pos < fg_npts)
				sacdata[k].sac_data[pos] = 1.0/fg_norm;
		} else if(fg_type == FG_TRIANGLE){
			pos = coff+fg_npts/2 -1;
			if(pos >= 0 && pos < fg_npts)
				sacdata[k].sac_data[pos] = 0.25/fg_norm;
			pos = coff+fg_npts/2;
			if(pos >= 0 && pos < fg_npts)
				sacdata[k].sac_data[pos]    = 0.50/fg_norm;
			pos = coff+fg_npts/2 +1;
			if(pos >= 0 && pos < fg_npts)
				sacdata[k].sac_data[pos] = 0.25/fg_norm;
		} else if(fg_type == FG_BOX){
			/* safety - if length < 10 *dt change it */
			nl = fg_length / fg_norm;
			if(nl < 10)nl = 10;
			for(i=0;i< nl; i++){
				pos = coff+fg_npts/2 +i;
				if(pos >= 0 && pos < fg_npts)
				sacdata[k].sac_data[pos] = 1.0/(fg_norm*nl);
			}
		} 
	}
	/* define integer header values */
	
	sacdata[k].sachdr.ihdr[H_NZYEAR] = 1970;
	sacdata[k].sachdr.ihdr[H_NZJDAY] = 1;
	sacdata[k].sachdr.ihdr[H_NZHOUR] = 0;
	sacdata[k].sachdr.ihdr[H_NZMIN] = 0;
	sacdata[k].sachdr.ihdr[H_NZSEC] = 0;
	sacdata[k].sachdr.ihdr[H_NZMSEC] = 0;
	sacdata[k].sachdr.ihdr[H_NVHDR] = 6;
	sacdata[k].sachdr.ihdr[8] = -12345;
	sacdata[k].sachdr.ihdr[H_NPTS] = fg_npts;
	/* set enumerated header values */
	sacdata[k].sachdr.ihdr[H_IFTYPE] = ENUM_ITIME;
	sacdata[k].sachdr.ihdr[H_IZTYPE] = ENUM_IO;
	/* set logical header values */
	sacdata[k].sachdr.ihdr[35] = 1;
	sacdata[k].sachdr.ihdr[H_LPSPOL] = 0;
	sacdata[k].sachdr.ihdr[H_LOVROK] = 1;
	sacdata[k].sachdr.ihdr[H_LCALDA] = 1;
	sacdata[k].sachdr.ihdr[H_LHDR5] = 0;
	/* set real header values */
	sacdata[k].sachdr.rhdr[H_DELTA] = fg_delta;
	sacdata[k].sachdr.rhdr[H_B] = -fg_delta*(fg_npts/2);
	sacdata[k].sachdr.rhdr[H_E] = -fg_delta*(fg_npts/2) + fg_delta*(fg_npts-1);
	getmxmn(sacdata[k].sac_data, fg_npts,&depmax, &depmin, &depmen,&indmax,&indmin);
	sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
	sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
	sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
	sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
	sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
	sacdata[k].display = YES;
	/* define the file name */
	gsac_control.number_otraces = 1 ;
	gsac_control.fmax = 0.5/sacdata[k].sachdr.rhdr[H_DELTA];
	if(fg_type == FG_IMPULSE){
		strcpy(sacdata[k].sac_ofile_name,"impulse.sac");
	} else if(fg_type == FG_TRIANGLE){
		strcpy(sacdata[k].sac_ofile_name,"triangle.sac");
	}
	/* set epoch time */
	htoe1(sacdata[k].sachdr.ihdr[H_NZYEAR], 
		sacdata[k].sachdr.ihdr[H_NZJDAY], 
		sacdata[k].sachdr.ihdr[H_NZHOUR],
		sacdata[k].sachdr.ihdr[H_NZMIN],
		sacdata[k].sachdr.ihdr[H_NZSEC],
		sacdata[k].sachdr.ihdr[H_NZMSEC],
		&sacdata[k].tzref);
	sacdata[k].tzbeg = sacdata[k].tzref + sacdata[k].sachdr.rhdr[H_B];
	sacdata[k].tzend = sacdata[k].tzref + sacdata[k].sachdr.rhdr[H_E];
	sacdata[k].tzbegx = sacdata[k].tzbeg;
	sacdata[k].tzendx = sacdata[k].tzend;
	/* get bounds for absolute plotting 
		this is easy since there is only one in memory at a time */
		gsac_control.begmin = sacdata[k].tzbeg;
		gsac_control.endmax = sacdata[k].tzend;
	return;
}
