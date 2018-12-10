/* Changes:
	20 JUL 2007  - added the option to shift the traces 
		prior to subtraction according to the T9
		variable, which would be set by the CORRELATE
		command.  This is under developement and the
		logic must be tested for failure. The purpose
		is for comparing two synthetics from different 
		techniques to define a goodness of fit criteria
*/

#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include	"csstim.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	SUBF_DFLT	0
#define	SUBF_MAST	1
#define	SUBF_SUF	3
#define	SUBF_ABSOLUTE	4
#define	SUBF_RELATIVE	5

static int   subf_master = 0;
static int   subf_absolute = YES;
static char  subf_suffix[80];
static int   subf_dosuffix = NO;

struct arghdr subfarg[] = {
	{SUBF_DFLT, "DEFAULT", IHDR, 0, 0, NO, "",  1},
	{SUBF_SUF , "SUFFIX" , CHDR, 0, 1, NO, "SUFFIX suffix",  1},
	{SUBF_MAST, "MASTER" , IHDR, 0, 1, NO, "MASTER n",  1},
	{SUBF_ABSOLUTE, "ABSOLUTE", IHDR, 0, 0, NO, "",  1},
	{SUBF_RELATIVE, "RELATIVE", IHDR, 0, 0, NO, "",  1},
	{0,	""		, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float subf_real[10];
int   subf_int [10];
int   subf_yn;
int   subf_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_subf(int ncmd, char **cmdstr)
{
	int i;
	char instr[80];
	/* parsing code here */
	/* set up default */
	subf_dosuffix = YES;
	subf_suffix[0]='\0';
		strcpy(subf_suffix,".sub");
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, subfarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; subfarg[i].key[0] != '\0' ; i++){
		if(subfarg[i].used > 0){
			if(subfarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, subfarg[i].key, 
					subfarg[i].mfit,subfarg[i].narg, subf_real);
			} else if(subfarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, subfarg[i].key, 
					subfarg[i].mfit,subfarg[i].narg, subf_int );
			} else if(subfarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, subfarg[i].key, 
					subfarg[i].mfit,subfarg[i].narg, &subf_yn );
			} else if(subfarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, subfarg[i].key, 
					subfarg[i].mfit,subfarg[i].narg, &subf_num );
			} else if(subfarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, subfarg[i].key, 
					subfarg[i].mfit, subfarg[i].narg, instr );
			}
			switch(subfarg[i].id){
				case SUBF_SUF:
					if(strlen(instr) < 80){
						strcpy(subf_suffix, instr);
						subf_dosuffix = YES;
					}
					break;
				case SUBF_DFLT:
					subf_dosuffix = YES;
					strcpy(subf_suffix,".mul");
					subf_master = 0;
					subf_absolute = YES;
					break;
				case SUBF_MAST:
					subf_master = subf_int[0];
					break;
				case SUBF_ABSOLUTE:
					subf_absolute = YES;
					break;
				case SUBF_RELATIVE:
					subf_absolute = NO;
					break;
			}
		}
	}
			
}


void gsac_exec_subf(void)
{

	double ts, te;
	int i, k;
	int npts, ntrc;
	float dtx, dty;
	float x0, y0;
	int nptsx;
	int noff_master, noff;
	float depmax, depmin, depmen;
	int indmax, indmin;
	int month, day;
	double tzref;
	float f_t9_offset;


	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
	if(subf_master < 0 || subf_master > ntrc){
		printf("subf_master = %d not in range of [0,%d]\n",
				subf_master,ntrc-1);
		return;
	}
	printf("SUBF MASTER %d   \n", subf_master );


	/* CHECK TO SEE THAT ALL DT ARE THE SAME ELSE RETURN 
		ALSO DETERMINE LATEST START TIME AND EARLIEST END TIME */
	/* note that we will use the idea of rotate3 to look at the
		common time window - unless we set relative in which 
		case we will ignore and use the shortest length */
	/* define the master trace parameters */
	nptsx = sacdata[subf_master].sachdr.ihdr[H_NPTS];
	dtx = sacdata[subf_master].sachdr.rhdr[H_DELTA];
	ts = sacdata[subf_master].tzbeg;
	te = sacdata[subf_master].tzend;
	tzref = sacdata[subf_master].tzref;
	for ( k=0 ; k < ntrc ; k ++){
		dty = sacdata[k].sachdr.rhdr[H_DELTA];
		if(ABS (dtx - dty) > 0.01 * ABS(dtx)){
			printf("DT not equal - subf not done\n");
			printf("Files %d %d with %f %f, respectively\n",
				subf_master,k,dtx,dty);
			return;
		}
			if(sacdata[k].tzbeg > ts)ts = sacdata[k].tzbeg;
			if(sacdata[k].tzend < te)te = sacdata[k].tzend;
	}
        if(ts > te){
                printf("Cannot subtract traces - no trace overlap\n");
                return;
        }
	/* define the number of points in the output file */
	npts  = MIN(((int)((te -ts)/dtx + 0.49) + 1),nptsx);
	/* now cycle through the data set performing the operation 
		also reset the headers */

	/* compute offset for the master trace */
	noff_master = (int)((ts - sacdata[subf_master].tzbeg)/dtx + 0.49);
	/* process by output point - this is inefficient in terms of 
		memory access but saves the effort of copying the 
		master trace since everything is
		overwritten in memory by this process */
	for(i=0;i< npts ; i++){
		x0 = sacdata[subf_master].sac_data[i + noff_master];
		noff = (int)((ts - sacdata[k].tzbeg)/dtx + 0.49);
		for ( k=0 ; k < ntrc ; k ++){
			/* define offset for each trace */
			noff = (int)((ts - sacdata[k].tzbeg)/dtx + 0.49);
			y0 = sacdata[k].sac_data[i + noff];
			sacdata[k].sac_data[i] = y0-x0;
		}
	}
	/* now update headers and output file name */
	for(k = 0 ; k < ntrc ; k++){
		getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
		sacdata[k].sachdr.ihdr[H_NPTS] = npts;
		sacdata[k].sachdr.rhdr[H_B] = ts - tzref;
		sacdata[k].sachdr.rhdr[H_E] = te - tzref;
		sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
		sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
		sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		/* update the headers */
		etoh(tzref, &sacdata[k].sachdr.ihdr[H_NZYEAR], 
			&sacdata[k].sachdr.ihdr[H_NZJDAY], &month, &day,
			&sacdata[k].sachdr.ihdr[H_NZHOUR], 
			&sacdata[k].sachdr.ihdr[H_NZMIN],
			&sacdata[k].sachdr.ihdr[H_NZSEC], 
			&sacdata[k].sachdr.ihdr[H_NZMSEC]);
	
		chofname(sacdata[subf_master].schdr[H_KSTNM],
			sacdata[subf_master].schdr[H_KCMPNM],
			sacdata[k].schdr[H_KSTNM],
			sacdata[k].schdr[H_KCMPNM],
			sacdata[k].sac_ofile_name,
			subf_dosuffix, subf_suffix);

	}


	
}

