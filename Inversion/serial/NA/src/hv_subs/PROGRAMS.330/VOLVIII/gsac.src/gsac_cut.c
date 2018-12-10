/* note that here the important thing is the positioning
 * and not the B E etc so be careful later
 *
 * In the design, we can use either the order
 * cut o -10 o 20, or cut o 20 o -10. When we actually do the cut
 * as part of do_read.c, then we will worry about the least and the greatest
 * time., as in
 * start = MIN(sacdata[k].sachdr.rhdr[CUT_INT1] + offset1, 
 * 		sacdata[k].sachdr.rhdr[CUT_INT2] + offset2)
 *  end  = MAX(sacdata[k].sachdr.rhdr[CUT_INT1] + offset1, 
 * 		sacdata[k].sachdr.rhdr[CUT_INT2] + offset2)
 *
 * */
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

#define	CUT_ON		-1
#define CUT_OFF		-2
#define CUT_B		5
#define CUT_E		6
#define CUT_O		7
#define CUT_A		8
#define CUT_T0		10
#define CUT_T1		11
#define CUT_T2		12
#define CUT_T3		13
#define CUT_T4		14
#define CUT_T5		15
#define CUT_T6		16
#define CUT_T7		17
#define CUT_T8		18
#define CUT_T9		19
#define CUT_GMT		-20
#define CUT_CAL		-21

struct arghdr cutarg[] = {
	{CUT_ON, "ON"	, IHDR, 0, 0, NO, "ON",2},
	{CUT_OFF, "OFF"	, IHDR, 0, 0, NO, "OFF",2},
	{CUT_B, "B"	, RHDR, 0, 1, NO, "B offset ",-1},
	{CUT_E, "E"	, RHDR, 0, 1, NO, "E offset ",-1},
	{CUT_O, "O"	, RHDR, 0, 1, NO, "O offset ",-1},
	{CUT_A, "A"	, RHDR, 0, 1, NO, "A offset ",-1},
	{CUT_T0, "T0"	, RHDR, 0, 1, NO, "T0 offset ",-1},
	{CUT_T1, "T1"	, RHDR, 0, 1, NO, "T1 offset ",-1},
	{CUT_T2, "T2"	, RHDR, 0, 1, NO, "T2 offset ",-1},
	{CUT_T3, "T3"	, RHDR, 0, 1, NO, "T3 offset ",-1},
	{CUT_T4, "T4"	, RHDR, 0, 1, NO, "T4 offset ",-1},
	{CUT_T5, "T5"	, RHDR, 0, 1, NO, "T5 offset ",-1},
	{CUT_T6, "T6"	, RHDR, 0, 1, NO, "T6 offset ",-1},
	{CUT_T7, "T7"	, RHDR, 0, 1, NO, "T7 offset ",-1},
	{CUT_T8, "T8"	, RHDR, 0, 1, NO, "T8 offset ",-1},
	{CUT_T9, "T9"	, RHDR, 0, 1, NO, "T9 offset ",-1},
	{CUT_GMT, "GMT"  , IHDR, 0, 6, YES, "GMT YY DOY HH MM SS MSEC",-1},
	{CUT_CAL, "CAL"  , IHDR, 0, 7, YES, "CAL YY MO DD HH MN SS MSEC",-1},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};

static float cut_real[10];
static int cut_int[10];
int cutint[2];
float cutoff[2];
extern int desac(int mcmd, char **tcmdstr, int *ncmd, char *cmdstr[5]);

	
static char *cmdstr[] = { (char *)NULL, (char *)NULL, 
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL,(char *)NULL,
			(char *)NULL};


static int cutyear[2];		/* for use with cut CAL or GMT */
static int cutjday[2];
static int cutmon[2];
static int cutday[2];
static int cuthour[2];
static int cutmin[2];
static int cutsec[2];
static int cutmsec[2];

static char ostr[100];

void gsac_set_param_cut(int mcmd, char **tcmdstr)
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
for(i=0;i<mcmd;i++)printf("%s ",tcmdstr[i]);printf("\n");
	if(desac(mcmd,tcmdstr,&ncmd,cmdstr)==NO){
		printf("Use syntax\n  CUT ref offset ref offset | GMT YR DY HR MN SEC MSEC | CAL YR MO DY HR MN SEC MSEC\n");
		gsac_control.docut = NO;
		return;
	}

	/* test what was returned 
		for(i=0 ; i < ncmd; i++)printf("%s ",cmdstr[i]);printf("\n");
	*/

	if(testarg(ncmd, cmdstr, cutarg, NO, YES)){
		gsac_control.docut = NO;
		return  ;
	}

	cnt = 0;
	for(i=0 ; cutarg[i].key[0] != '\0' && cnt < 2 ; i++){
		if(cutarg[i].used > 0){
			switch(cutarg[i].id){
				case CUT_ON:
					gsac_control.docut = YES;
					break;
				case CUT_OFF:
					gsac_control.docut = NO;
					break;
				case CUT_GMT:
					/* the use in ind get position within line */ 
					ind =getargi(ncmd, cmdstr, cutarg[i].key, 
					cutarg[i].mfit, cutarg[i].narg, cut_int);
					cutyear[cnt]=cut_int[0];
					cutjday[cnt]=cut_int[1];
					cuthour[cnt]=cut_int[2];
					cutmin[cnt]=cut_int[3];
					cutsec[cnt]=cut_int[4];
					cutmsec[cnt]=cut_int[5];
					htoe1(cut_int[0], 
						cut_int[1],
						cut_int[2],
						cut_int[3],
						cut_int[4],
						cut_int[5],
					&gsac_control.cutepoch[cnt]);
					gsac_control.docut = YES;
					gsac_control.cutint[cnt]
						=cutarg[i].id;
printf("GMT %lf ",gsac_control.cutepoch[cnt]);
printtimestr(gsac_control.cutepoch[cnt], ostr); printf("%s\n",ostr);
					cnt++;
					if(cutarg[i].used > 6){
					ind = getargi(ncmd-ind, &cmdstr[ind], cutarg[i].key, 
						cutarg[i].mfit, cutarg[i].narg, cut_int);
						cutyear[cnt]=cut_int[0];
						cutjday[cnt]=cut_int[1];
						cuthour[cnt]=cut_int[2];
						cutmin[cnt]=cut_int[3];
						cutsec[cnt]=cut_int[4];
						cutmsec[cnt]=cut_int[5];
						htoe1(cut_int[0], 
							cut_int[1],
							cut_int[2],
							cut_int[3],
							cut_int[4],
							cut_int[5],
						&gsac_control.cutepoch[cnt]);
					gsac_control.cutint[cnt]=cutarg[i].id;
printf("GMT %lf ",gsac_control.cutepoch[cnt]);
printtimestr(gsac_control.cutepoch[cnt], ostr); printf("%s\n",ostr);
					cnt++;
					}
					break;
				case CUT_CAL:
					ind =getargi(ncmd, cmdstr, cutarg[i].key, 
						cutarg[i].mfit, cutarg[i].narg, cut_int);
						cutyear[cnt]=cut_int[0];
						cutmon[cnt]=cut_int[1];
						cutday[cnt]=cut_int[2];
						cuthour[cnt]=cut_int[3]; 
							cutmin[cnt]=cut_int[4]; 
							cutsec[cnt]=cut_int[5];
							cutmsec[cnt]=cut_int[6]; 						htoe2(cut_int[0], 
							cut_int[1], 
							cut_int[2], 
							cut_int[3], 
							cut_int[4],
							cut_int[5],
							cut_int[6],
						&gsac_control.cutepoch[cnt]);
printf("CAL %lf ",gsac_control.cutepoch[cnt]);
printtimestr(gsac_control.cutepoch[cnt], ostr); printf("%s\n",ostr);
						gsac_control.docut = YES;
						gsac_control.cutint[cnt]
							=cutarg[i].id;
						cnt++;
					if(cutarg[i].used > 7){
					ind = getargi(ncmd-ind, &cmdstr[ind], cutarg[i].key, 
						cutarg[i].mfit, cutarg[i].narg, cut_int);
						cutyear[cnt]=cut_int[0];
						cutmon[cnt]=cut_int[1];
						cutday[cnt]=cut_int[2];
						cuthour[cnt]=cut_int[3];
						cutmin[cnt]=cut_int[4];
						cutsec[cnt]=cut_int[5];
						cutmsec[cnt]=cut_int[6];
						cutint[cnt]=cutarg[i].id;
						htoe2(cut_int[0], 
							cut_int[1],
							cut_int[2],
							cut_int[3],
							cut_int[4],
							cut_int[5],
							cut_int[6],
						&gsac_control.cutepoch[cnt]);
					gsac_control.cutint[cnt]=cutarg[i].id;
printf("CAL %lf ",gsac_control.cutepoch[cnt]);
printtimestr(gsac_control.cutepoch[cnt], ostr); printf("%s\n",ostr);
					cnt++;
					}
					break;

				case  CUT_B  :
				case  CUT_E  :
				case  CUT_O  :
				case  CUT_A  :
				case  CUT_T0 :
				case  CUT_T1 :
				case  CUT_T2 :
				case  CUT_T3 :
				case  CUT_T4 :
				case  CUT_T5 :
				case  CUT_T6 :
				case  CUT_T7 :
				case  CUT_T8 :
				case  CUT_T9 :
					ind =getargr(ncmd, cmdstr, cutarg[i].key, 
						cutarg[i].mfit, cutarg[i].narg, cut_real);
						gsac_control.docut = YES;
						gsac_control.cutint[cnt]=cutarg[i].id;
						gsac_control.cutoff[cnt] = cut_real[0];
						strcpy(gsac_control.cutkey[cnt],cutarg[i].key);
					cnt++;
					if(cutarg[i].used > 1){
					/* special case like cut o -10 o 20
					 *  instead of simpler case of
					 *  cut a -10 t0 +20
					 */
					ind = getargr(ncmd-ind, &cmdstr[ind], cutarg[i].key, 
						cutarg[i].mfit, cutarg[i].narg, cut_real);
					gsac_control.cutoff[cnt] = cut_real[0];
					strcpy(gsac_control.cutkey[cnt],cutarg[i].key);
					gsac_control.cutint[cnt]=cutarg[i].id;
					}
					break;
			}		
		}
	}
}

void gsac_exec_cut(void)
{
	int i;
	/* note that gsac_control.cutkey is used only for nice output here */
	if(gsac_control.docut){
		for(i=0 ; i < 2 ; i++){
			if(gsac_control.cutint[i] == CUT_GMT){
				printf("GMT %04d %03d %02d %02d %02d %03d ",
					cutyear[i],
					cutjday[i],
					cuthour[i],
					cutmin[i],
					cutsec[i],
					cutmsec[i]);
			} else if(gsac_control.cutint[i] == CUT_CAL){     
				printf("CAL  %04d %02d %02d %02d %02d %02d %03d ",
					cutyear[i],
					cutmon[i],
					cutday[i],
					cuthour[i],
					cutmin[i],
					cutsec[i],
					cutmsec[i]);
			} else {
				printf("%s %f ",
					gsac_control.cutkey[i],
					gsac_control.cutoff[i]);
			}
		}
		printf("\n");	
	} else {
		printf("CUT is turned off\n");
	}

}
