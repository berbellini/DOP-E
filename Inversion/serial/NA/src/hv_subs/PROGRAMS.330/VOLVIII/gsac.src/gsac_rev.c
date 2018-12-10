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


#define	REV_DFLT	0
#define	REV_SUF		1

static char  rev_suffix[80];
static float rev_dosuffix = NO;


struct arghdr revarg[] = {
	{REV_DFLT, "DEFAULT", IHDR, 0, 0, NO, "", -1},
	{REV_SUF , "SUFFIX" , CHDR, 0, 1, NO, "SUFFIX suffix x0", 1},
	{0,	""		, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float rev_real[10];
int   rev_int [10];
int   rev_yn;
int   rev_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_rev(int ncmd, char **cmdstr)
{
	int i;
	char instr[80];
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, revarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; revarg[i].key[0] != '\0' ; i++){
		if(revarg[i].used > 0){
			if(revarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, revarg[i].key, 
					revarg[i].mfit,revarg[i].narg, rev_real);
			} else if(revarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, revarg[i].key, 
					revarg[i].mfit,revarg[i].narg, rev_int );
			} else if(revarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, revarg[i].key, 
					revarg[i].mfit,revarg[i].narg, &rev_yn );
			} else if(revarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, revarg[i].key, 
					revarg[i].mfit,revarg[i].narg, &rev_num );
			} else if(revarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, revarg[i].key, 
					revarg[i].mfit, revarg[i].narg, instr );
			}
			switch(revarg[i].id){
				case REV_SUF:
					if(strlen(instr) < 80){
						strcpy(rev_suffix, instr);
						rev_dosuffix = YES;
					}
					break;
				case REV_DFLT:
					rev_dosuffix = NO;
			}
		}
	}
		
}

void gsac_exec_rev(void)
{
	int k, ntrc, i, j, nup, npts;
	float tmp;
	double dtmp;
	double otzref;
	double ntzref;
	int month, day;
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		if( (npts%2) == 0){
			/* even */
			nup = npts/2 ;
		} else {
			/* odd  */
			nup = npts/2 ;
		}
		nup = npts/2 ;
		for(i=0, j=npts-1 ; i < nup ; i++, j--){
			tmp = sacdata[k].sac_data[i];
			sacdata[k].sac_data[i] = sacdata[k].sac_data[j];
			sacdata[k].sac_data[j] = tmp;
		}
		if(rev_dosuffix == YES){
			strcat(sacdata[k].sac_ofile_name,rev_suffix);
		} else {
			strcat(sacdata[k].sac_ofile_name,".rev");
		}
		/* systematically change all time headers */ 
		htoe1(sacdata[k].sachdr.ihdr[H_NZYEAR], 
			sacdata[k].sachdr.ihdr[H_NZJDAY], 
			sacdata[k].sachdr.ihdr[H_NZHOUR],
			sacdata[k].sachdr.ihdr[H_NZMIN],
			sacdata[k].sachdr.ihdr[H_NZSEC],
			sacdata[k].sachdr.ihdr[H_NZMSEC],
			&otzref);
		ntzref = -otzref ;
		sacdata[k].tzref =  ntzref;
		dtmp = sacdata[k].tzbeg;
		sacdata[k].tzbeg = - sacdata[k].tzend ;
		sacdata[k].tzend = -dtmp;
		sacdata[k].tzbegx = sacdata[k].tzbeg;
		sacdata[k].tzendx = sacdata[k].tzend;
		etoh(sacdata[k].tzref, 
				&sacdata[k].sachdr.ihdr[H_NZYEAR], 
				&sacdata[k].sachdr.ihdr[H_NZJDAY], &month, &day,
				&sacdata[k].sachdr.ihdr[H_NZHOUR], 
				&sacdata[k].sachdr.ihdr[H_NZMIN],
				&sacdata[k].sachdr.ihdr[H_NZSEC], 
				&sacdata[k].sachdr.ihdr[H_NZMSEC]);

		if(sacdata[k].sachdr.rhdr[H_A] != -12345.)
			sacdata[k].sachdr.rhdr[H_A] = 
				-sacdata[k].sachdr.rhdr[H_A];
		if(sacdata[k].sachdr.rhdr[H_O] != -12345.)
			sacdata[k].sachdr.rhdr[H_O] = 
				-sacdata[k].sachdr.rhdr[H_O];
		if(sacdata[k].sachdr.rhdr[H_T0] != -12345.)
			sacdata[k].sachdr.rhdr[H_T0] = 
				-sacdata[k].sachdr.rhdr[H_T0];
		if(sacdata[k].sachdr.rhdr[H_T1] != -12345.)
			sacdata[k].sachdr.rhdr[H_T1] = 
				-sacdata[k].sachdr.rhdr[H_T1];
		if(sacdata[k].sachdr.rhdr[H_T2] != -12345.)
			sacdata[k].sachdr.rhdr[H_T2] = 
				-sacdata[k].sachdr.rhdr[H_T2];
		if(sacdata[k].sachdr.rhdr[H_T3] != -12345.)
			sacdata[k].sachdr.rhdr[H_T3] = 
				-sacdata[k].sachdr.rhdr[H_T3];
		if(sacdata[k].sachdr.rhdr[H_T4] != -12345.)
			sacdata[k].sachdr.rhdr[H_T4] = 
				-sacdata[k].sachdr.rhdr[H_T4];
		if(sacdata[k].sachdr.rhdr[H_T5] != -12345.)
			sacdata[k].sachdr.rhdr[H_T5] = 
				-sacdata[k].sachdr.rhdr[H_T5];
		if(sacdata[k].sachdr.rhdr[H_T6] != -12345.)
			sacdata[k].sachdr.rhdr[H_T6] = 
				-sacdata[k].sachdr.rhdr[H_T6];
		if(sacdata[k].sachdr.rhdr[H_T7] != -12345.)
			sacdata[k].sachdr.rhdr[H_T7] = 
				-sacdata[k].sachdr.rhdr[H_T7];
		if(sacdata[k].sachdr.rhdr[H_T8] != -12345.)
			sacdata[k].sachdr.rhdr[H_T8] = 
				-sacdata[k].sachdr.rhdr[H_T8];
		if(sacdata[k].sachdr.rhdr[H_T9] != -12345.)
			sacdata[k].sachdr.rhdr[H_T9] = 
				-sacdata[k].sachdr.rhdr[H_T9];
		/* special care for B and E since - old B -> new E */
		tmp = sacdata[k].sachdr.rhdr[H_B];
		sacdata[k].sachdr.rhdr[H_B] = - sacdata[k].sachdr.rhdr[H_E];
		sacdata[k].sachdr.rhdr[H_E] = -tmp;
		
	 
	}
	/* redefine the bounds for absolute plotting */
	gsac_control.begmin= 1.0e+30;
	gsac_control.endmax=-1.0e+30;
	for ( k=0 ; k < ntrc ; k ++){
		if(sacdata[k].tzbeg < gsac_control.begmin)
			gsac_control.begmin = sacdata[k].tzbeg;
		if(sacdata[k].tzend > gsac_control.endmax)
			gsac_control.endmax = sacdata[k].tzend;
	}

}
