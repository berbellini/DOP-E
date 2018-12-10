#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include	"csstim.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	DIVF_DFLT	0
#define	DIVF_MAST	1
#define	DIVF_SUF	3
#define	DIVF_ABSOLUTE	4
#define	DIVF_RELATIVE	5
#define	DIVF_WATER	6

static int   divf_master = 0;
static int   divf_absolute = YES;
static char  divf_suffix[80];
static int   divf_dosuffix = NO;
static float divf_waterlevel = 0.0001;


struct arghdr divfarg[] = {
	{DIVF_DFLT, "DEFAULT", IHDR, 0, 0, NO, "",  1},
	{DIVF_SUF , "SUFFIX" , CHDR, 0, 1, NO, "SUFFIX suffix",  1},
	{DIVF_MAST, "MASTER" , IHDR, 0, 1, NO, "MASTER n",  1},
	{DIVF_ABSOLUTE, "ABSOLUTE", IHDR, 0, 0, NO, "",  1},
	{DIVF_RELATIVE, "RELATIVE", IHDR, 0, 0, NO, "",  1},
	{0,	""		, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float divf_real[10];
int   divf_int [10];
int   divf_yn;
int   divf_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_divf(int ncmd, char **cmdstr)
{
	int i;
	char instr[80];
	/* parsing code here */
	/* set up default */
	divf_dosuffix = YES;
	divf_suffix[0]='\0';
		strcpy(divf_suffix,".div");
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, divfarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; divfarg[i].key[0] != '\0' ; i++){
		if(divfarg[i].used > 0){
			if(divfarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, divfarg[i].key, 
					divfarg[i].mfit,divfarg[i].narg, divf_real);
			} else if(divfarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, divfarg[i].key, 
					divfarg[i].mfit,divfarg[i].narg, divf_int );
			} else if(divfarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, divfarg[i].key, 
					divfarg[i].mfit,divfarg[i].narg, &divf_yn );
			} else if(divfarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, divfarg[i].key, 
					divfarg[i].mfit,divfarg[i].narg, &divf_num );
			} else if(divfarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, divfarg[i].key, 
					divfarg[i].mfit, divfarg[i].narg, instr );
			}
			switch(divfarg[i].id){
				case DIVF_SUF:
					if(strlen(instr) < 80){
						strcpy(divf_suffix, instr);
						divf_dosuffix = YES;
					}
					break;
				case DIVF_DFLT:
					divf_dosuffix = YES;
					strcpy(divf_suffix,".mul");
					divf_master = 0;
					divf_absolute = YES;
					break;
				case DIVF_MAST:
					divf_master = divf_int[0];
					break;
				case DIVF_WATER:
					divf_waterlevel = divf_real[0];
					break;
				case DIVF_ABSOLUTE:
					divf_absolute = YES;
					break;
				case DIVF_RELATIVE:
					divf_absolute = NO;
					break;
			}
		}
	}
			
}


void gsac_exec_divf(void)
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
	float water;


	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
	if(divf_master < 0 || divf_master > ntrc){
		printf("divf_master = %d not in range of [0,%d]\n",
				divf_master,ntrc-1);
		return;
	}
	printf("DIVF MASTER %d Water level factor %f  \n", divf_master,divf_waterlevel );


	/* CHECK TO SEE THAT ALL DT ARE THE SAME ELSE RETURN 
		ALSO DETERMINE LATEST START TIME AND EARLIEST END TIME */
	/* note that we will use the idea of rotate3 to look at the
		common time window - unless we set relative in which 
		case we will ignore and use the shortest length */
	/* define the master trace parameters */
	nptsx = sacdata[divf_master].sachdr.ihdr[H_NPTS];
	dtx = sacdata[divf_master].sachdr.rhdr[H_DELTA];
	ts = sacdata[divf_master].tzbeg;
	te = sacdata[divf_master].tzend;
	tzref = sacdata[divf_master].tzref;
	for ( k=0 ; k < ntrc ; k ++){
		dty = sacdata[k].sachdr.rhdr[H_DELTA];
		if(ABS (dtx - dty) > 0.01 * ABS(dtx)){
			printf("DT not equal - divf not done\n");
			printf("Files %d %d with %f %f, respectively\n",
				divf_master,k,dtx,dty);
			return;
		}
		if(sacdata[k].tzbeg > ts)ts = sacdata[k].tzbeg;
		if(sacdata[k].tzend < te)te = sacdata[k].tzend;
	}
        if(ts > te){
                printf("Cannot divide traces - no trace overlap\n");
                return;
        }
	/* define the number of points in the output file */
	npts  = MIN(((int)((te -ts)/dtx + 0.49) + 1),nptsx);
	/* now cycle through the data set performing the operation 
		also reset the headers */

	/* compute offset for the master trace */
	noff_master = (int)((ts - sacdata[divf_master].tzbeg)/dtx + 0.49);
	/* get the maximum amplitude within the common region */
	water = 0.0 ;
	for(i=0 ; i < npts ; i ++)
		if(sacdata[divf_master].sac_data[i +noff_master] > water)
			water = sacdata[divf_master].sac_data[i +noff_master];
	water = water*divf_waterlevel;
	if(water <= 0.0)
		water = 1.0;
	/* process by output point - this is inefficient in terms of 
		memory access but saves the effort of copying the 
		master trace since everything is
		overwritten in memory by this process */
	for(i=0;i< npts ; i++){
		x0 = sacdata[divf_master].sac_data[i + noff_master];
		noff = (int)((ts - sacdata[k].tzbeg)/dtx + 0.49);
		for ( k=0 ; k < ntrc ; k ++){
			/* define offset for each trace */
			noff = (int)((ts - sacdata[k].tzbeg)/dtx + 0.49);
			y0 = sacdata[k].sac_data[i + noff];
			sacdata[k].sac_data[i] = y0/MAX(x0,water);
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
	
		chofname(sacdata[divf_master].schdr[H_KSTNM],
			sacdata[divf_master].schdr[H_KCMPNM],
			sacdata[k].schdr[H_KSTNM],
			sacdata[k].schdr[H_KCMPNM],
			sacdata[k].sac_ofile_name,
			divf_dosuffix, divf_suffix);

	}


	
}

