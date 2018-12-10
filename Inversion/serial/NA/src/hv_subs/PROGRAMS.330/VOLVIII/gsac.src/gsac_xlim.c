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


void gsac_exec_doxlim(void);

/* note that here the important thing is the positioning
 * and not the B E etc so be careful later
 *
 * In the design, we can use either the order
 * xlim o -10 o 20, or xlim o 20 o -10. When we actually do the xlim
 * as part of do_read.c, then we will worry about the least and the greatest
 * time., as in
 * start = MIN(sacdata[k].sachdr.rhdr[XLIM_INT1] + offset1, 
 * 		sacdata[k].sachdr.rhdr[XLIM_INT2] + offset2)
 *  end  = MAX(sacdata[k].sachdr.rhdr[XLIM_INT1] + offset1, 
 * 		sacdata[k].sachdr.rhdr[XLIM_INT2] + offset2)
 *
 * */
#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

/* the numbering is crucial since the desac routine will know that if the
	flag is > 0 that the OLD sac careless syntax must be recognized. If
	the input were always of the form
		cut b 10 b 20
	instead of
		cut 10 20
		cut b 10 20
		cut 10 b 20
	we would not need the desac routine and its additional complexity because of
	the addition of the CAL GMT ON OFF flags. Note that a
		cut a 0 o 0
	is the same as
		cut a  o 
		cut a o 0
		cut a 0 o
	at least these have letter codes
	desac is shared between cut and xlim so the same numbering care and naming must be
	applied to both
*/


#define	XLIM_ON		-1
#define XLIM_OFF	-2
#define XLIM_B		5
#define XLIM_E		6
#define XLIM_O		7
#define XLIM_A		8
#define XLIM_T0		10
#define XLIM_T1		11
#define XLIM_T2		12
#define XLIM_T3		13
#define XLIM_T4		14
#define XLIM_T5		15
#define XLIM_T6		16
#define XLIM_T7		17
#define XLIM_T8		18
#define XLIM_T9		19
#define XLIM_GMT	-20
#define XLIM_CAL	-21

struct arghdr xlimarg[] = {
	{XLIM_ON, "ON"	, IHDR, 0, 0, NO, "ON",-1},
	{XLIM_OFF, "OFF", IHDR, 0, 0, NO, "OFF", 2},
	{XLIM_B, "B"	, RHDR, 0, 1, NO, "B offset ",-1},
	{XLIM_E, "E"	, RHDR, 0, 1, NO, "E offset ",-1},
	{XLIM_O, "O"	, RHDR, 0, 1, NO, "O offset ",-1},
	{XLIM_A, "A"	, RHDR, 0, 1, NO, "A offset ",-1},
	{XLIM_T0, "T0"	, RHDR, 0, 1, NO, "T0 offset ",-1},
	{XLIM_T1, "T1"	, RHDR, 0, 1, NO, "T1 offset ",-1},
	{XLIM_T2, "T2"	, RHDR, 0, 1, NO, "T2 offset ",-1},
	{XLIM_T3, "T3"	, RHDR, 0, 1, NO, "T3 offset ",-1},
	{XLIM_T4, "T4"	, RHDR, 0, 1, NO, "T4 offset ",-1},
	{XLIM_T5, "T5"	, RHDR, 0, 1, NO, "T5 offset ",-1},
	{XLIM_T6, "T6"	, RHDR, 0, 1, NO, "T6 offset ",-1},
	{XLIM_T7, "T7"	, RHDR, 0, 1, NO, "T7 offset ",-1},
	{XLIM_T8, "T8"	, RHDR, 0, 1, NO, "T8 offset ",-1},
	{XLIM_T9, "T9"	, RHDR, 0, 1, NO, "T9 offset ",-1},
	{XLIM_GMT, "GMT", IHDR, 0, 6, NO, "GMT YY DOY HH MM SS MSEC",-1},
	{XLIM_CAL, "CAL", IHDR, 0, 7, NO, "CAL YY MO DY HH MN SS MSEC",-1},
	{0,	""	, IHDR, 0, 0, NO, "",-1}
};

static float xlim_real[10];
static int xlim_int[10];
static int xlimint[2];
static double xlimepoch[2];
int desac(int mcmd, char **tcmdstr, int *ncmd, char *cmdstr[5]);
/* potentially we have the following for cut and xlim:
	an entry of
	cmd start end (3 entries)
		interpreted as B start B end
	to
	cmd CAL YR MO DY HR MN SE MS CAL YR MO DY HR MN SE MS (15 entries)
*/
	
static char *cmdstr[] = { (char *)NULL, (char *)NULL, 
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL};

static int xlimyear[2];		/* for use with xlim CAL or GMT */
static int xlimjday[2];
static int xlimmon[2];
static int xlimday[2];
static int xlimhour[2];
static int xlimmin[2];
static int xlimsec[2];
static int xlimmsec[2];

char timestr[100];

void gsac_set_param_xlim(int mcmd, char **tcmdstr)
{
	int i, cnt, ind;
	int ncmd;
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
	if(desac(mcmd,tcmdstr,&ncmd,cmdstr)==NO){
		printf("Use syntax\n  XLIM ref offset ref offset | GMT YR DY HR MN SEC MSEC | CAL YR MO DY HR MN SEC MSEC\n");
		gsac_control.doxlim = NO;
		return;
	}

	if(testarg(ncmd, &cmdstr[0], xlimarg, NO, YES))
		return  ;

	cnt = 0;
	for(i=0 ; xlimarg[i].key[0] != '\0' && cnt < 2 ; i++){
		if(xlimarg[i].used > 0){
			switch(xlimarg[i].id){
				case XLIM_ON:
					gsac_control.doxlim = YES;
					break;
				case XLIM_OFF:
					gsac_control.doxlim = NO;
					break;
				case XLIM_GMT:
					/* the use in ind get position within line */ 
					ind =getargi(ncmd, cmdstr, xlimarg[i].key, 
					xlimarg[i].mfit, xlimarg[i].narg, xlim_int);
					xlimyear[cnt]=xlim_int[0];
					xlimjday[cnt]=xlim_int[1];
					xlimhour[cnt]=xlim_int[2];
					xlimmin[cnt]=xlim_int[3];
					xlimsec[cnt]=xlim_int[4];
					xlimmsec[cnt]=xlim_int[5];
					htoe1(xlim_int[0], 
						xlim_int[1],
						xlim_int[2],
						xlim_int[3],
						xlim_int[4],
						xlim_int[5],
					&xlimepoch[cnt]);
					gsac_control.doxlim = YES;
					gsac_control.xlimint[cnt]
						=xlimarg[i].id;
					cnt++;
					if(xlimarg[i].used > 6){
					ind = getargi(ncmd-ind, &cmdstr[ind], xlimarg[i].key, 
						xlimarg[i].mfit, xlimarg[i].narg, xlim_int);
						xlimyear[cnt]=xlim_int[0];
						xlimjday[cnt]=xlim_int[1];
						xlimhour[cnt]=xlim_int[2];
						xlimmin[cnt]=xlim_int[3];
						xlimsec[cnt]=xlim_int[4];
						xlimmsec[cnt]=xlim_int[5];
						htoe1(xlim_int[0], 
							xlim_int[1],
							xlim_int[2],
							xlim_int[3],
							xlim_int[4],
							xlim_int[5],
						&xlimepoch[cnt]);
					gsac_control.xlimint[cnt]=xlimarg[i].id;
					cnt++;
					}
					break;
				case XLIM_CAL:
					ind =getargi(ncmd, cmdstr, xlimarg[i].key, 
						xlimarg[i].mfit, xlimarg[i].narg, xlim_int);
						xlimyear[cnt]=xlim_int[0];
						xlimmon[cnt]=xlim_int[1];
						xlimday[cnt]=xlim_int[2];
						xlimhour[cnt]=xlim_int[3]; 
							xlimmin[cnt]=xlim_int[4]; 
							xlimsec[cnt]=xlim_int[5];
							xlimmsec[cnt]=xlim_int[6]; 						htoe2(xlim_int[0], 
							xlim_int[1], 
							xlim_int[2], 
							xlim_int[3], 
							xlim_int[4],
							xlim_int[5],
							xlim_int[6],
						&xlimepoch[cnt]);
						gsac_control.doxlim = YES;
						gsac_control.xlimint[cnt]
							=xlimarg[i].id;
						cnt++;
					if(xlimarg[i].used > 7){
					ind = getargi(ncmd-ind, &cmdstr[ind], xlimarg[i].key, 
						xlimarg[i].mfit, xlimarg[i].narg, xlim_int);
						xlimyear[cnt]=xlim_int[0];
						xlimmon[cnt]=xlim_int[1];
						xlimday[cnt]=xlim_int[2];
						xlimhour[cnt]=xlim_int[3];
						xlimmin[cnt]=xlim_int[4];
						xlimsec[cnt]=xlim_int[5];
						xlimmsec[cnt]=xlim_int[6];
						xlimint[cnt]=xlimarg[i].id;
						htoe2(xlim_int[0], 
							xlim_int[1],
							xlim_int[2],
							xlim_int[3],
							xlim_int[4],
							xlim_int[5],
							xlim_int[6],
						&xlimepoch[cnt]);
					gsac_control.xlimint[cnt]=xlimarg[i].id;
					cnt++;
					}
					break;
				default:

					ind =getargr(ncmd, cmdstr, xlimarg[i].key, 
						xlimarg[i].mfit, xlimarg[i].narg, xlim_real);
						gsac_control.doxlim = YES;
						gsac_control.xlimint[cnt]=xlimarg[i].id;
						gsac_control.xlimoff[cnt] = xlim_real[0];
						strcpy(gsac_control.xlimkey[cnt],xlimarg[i].key);
					cnt++;
					if(xlimarg[i].used > 1){
					/* special case like xlim o -10 o 20
					 *  instead of simpler case of
					 *  xlim a -10 t0 +20
					 */
					ind = getargr(ncmd-ind, &cmdstr[ind], xlimarg[i].key, 
						xlimarg[i].mfit, xlimarg[i].narg, xlim_real);
					gsac_control.xlimoff[cnt] = xlim_real[0];
					strcpy(gsac_control.xlimkey[cnt],xlimarg[i].key);
					gsac_control.xlimint[cnt]=xlimarg[i].id;
					}
					break;
			}		
		}
	}
			
		
}

void gsac_exec_xlim(void)
{
}


void gsac_exec_doxlim(void)
{
	int ntrc ;
	int k, i;
	double tw[2];
	int cnt;
	int doxlim;
	int success ;
	ntrc = gsac_control.number_otraces;
	if(gsac_control.doxlim){
		printf("xlim ");
		for(i=0 ; i < 2 ; i++){
			if(gsac_control.xlimint[i] == XLIM_GMT){
				printf("GMT %04d %03d %02d %02d %02d %03d ",
					xlimyear[i],
					xlimjday[i],
					xlimhour[i],
					xlimmin[i],
					xlimsec[i],
					xlimmsec[i]);
			} else if(gsac_control.xlimint[i] == XLIM_CAL){     
				printf("CAL %04d %02d %02d %02d %02d %02d %03d ",
					xlimyear[i],
					xlimmon[i],
					xlimday[i],
					xlimhour[i],
					xlimmin[i],
					xlimsec[i],
					xlimmsec[i]);
			} else {
				printf("%s %f ",
					gsac_control.xlimkey[i],
					gsac_control.xlimoff[i]);
			}
		}
		printf("\n");	
/*
		printf("xlimint %d xlimoff %f\n",
			gsac_control.xlimint[0],gsac_control.xlimoff[0]);
		printf("xlimint %d xlimoff %f\n",
			gsac_control.xlimint[1],gsac_control.xlimoff[1]);
*/
		 /* now loop through all traces and set or reset the 
		  * tzbegx and tzbege */
		if(ntrc < 1)
			 return ;
		/* process all traces */

		gsac_control.begminx=   1.0e+36 ;
		gsac_control.endmaxx=  -1.0e+36 ;
		
		success = YES;

		for ( k=0 ; k < ntrc ; k++){
			/* define the xlim window */
			doxlim = 0;
			for(cnt = 0 ; cnt < 2 ; cnt ++){
			if(gsac_control.xlimint[cnt] == XLIM_GMT){
				tw[cnt] = xlimepoch[cnt];
				doxlim++;
			} else if(gsac_control.xlimint[cnt] == XLIM_CAL){
				tw[cnt] = xlimepoch[cnt];
				doxlim++;
			} else {

			if(sacdata[k].sachdr.rhdr[gsac_control.xlimint[cnt]]!= -12345.){
			/* only do this if the headers are set */
				tw[cnt] =  sacdata[k].tzref  
				+ sacdata[k].sachdr.rhdr[gsac_control.xlimint[cnt]]
				+ gsac_control.xlimoff[cnt] ;
				doxlim++;
			} else {
printf("%s is not set for trace %s\n",gsac_control.xlimkey[cnt],sacdata[k].sac_ifile_name);
success = NO;

			}
			} 
			if(doxlim == 2){
			sacdata[k].tzbegx =  MAX(sacdata[k].tzbeg,MIN(tw[0],tw[1]));
			sacdata[k].tzendx =  MIN(sacdata[k].tzend,MAX(tw[0],tw[1]));
			if(sacdata[k].tzbegx < gsac_control.begminx){
				gsac_control.begminx = sacdata[k].tzbegx;
			} 
			
			if(sacdata[k].tzendx > gsac_control.endmaxx){
				gsac_control.endmaxx = sacdata[k].tzendx;
			}

			}
		}
		}
		if(doxlim != 2 || success == NO){
			printf(" xlim not applied\n");
			/* safety */
			gsac_control.begminx= gsac_control.begmin ;
			gsac_control.endmaxx= gsac_control.endmax ;
		}
		

	} else {
		printf("XLIM is turned off\n");
		if(ntrc < 1)
			 return ;
		gsac_control.begminx= gsac_control.begmin ;
		gsac_control.endmaxx= gsac_control.endmax ;
                /* 28 MAR 2009 reset the tzbegx for each trace */
		for ( k=0 ; k < ntrc ; k++){
			sacdata[k].tzbegx =  sacdata[k].tzbeg;
		}
	}

}


