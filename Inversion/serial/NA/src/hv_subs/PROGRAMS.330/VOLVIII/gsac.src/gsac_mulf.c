#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include	"csstim.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	MULF_DFLT	0
#define	MULF_MAST	1
#define	MULF_SUF	3
#define	MULF_ABSOLUTE	4
#define	MULF_RELATIVE	5

static int   mulf_master = 0;
static int   mulf_absolute = YES;
static char  mulf_suffix[80];
static int   mulf_dosuffix = NO;

struct arghdr mulfarg[] = {
	{MULF_DFLT, "DEFAULT", IHDR, 0, 0, NO, "",  1},
	{MULF_SUF , "SUFFIX" , CHDR, 0, 1, NO, "SUFFIX suffix",  1},
	{MULF_MAST, "MASTER" , IHDR, 0, 1, NO, "MASTER n",  1},
	{MULF_ABSOLUTE, "ABSOLUTE", IHDR, 0, 0, NO, "",  1},
	{MULF_RELATIVE, "RELATIVE", IHDR, 0, 0, NO, "",  1},
	{0,	""		, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float mulf_real[10];
int   mulf_int [10];
int   mulf_yn;
int   mulf_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_mulf(int ncmd, char **cmdstr)
{
	int i;
	char instr[80];
	/* parsing code here */
	/* set up default */
	mulf_dosuffix = YES;
	mulf_suffix[0]='\0';
		strcpy(mulf_suffix,".mul");
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, mulfarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; mulfarg[i].key[0] != '\0' ; i++){
		if(mulfarg[i].used > 0){
			if(mulfarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, mulfarg[i].key, 
					mulfarg[i].mfit,mulfarg[i].narg, mulf_real);
			} else if(mulfarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, mulfarg[i].key, 
					mulfarg[i].mfit,mulfarg[i].narg, mulf_int );
			} else if(mulfarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, mulfarg[i].key, 
					mulfarg[i].mfit,mulfarg[i].narg, &mulf_yn );
			} else if(mulfarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, mulfarg[i].key, 
					mulfarg[i].mfit,mulfarg[i].narg, &mulf_num );
			} else if(mulfarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, mulfarg[i].key, 
					mulfarg[i].mfit, mulfarg[i].narg, instr );
			}
			switch(mulfarg[i].id){
				case MULF_SUF:
					if(strlen(instr) < 80){
						strcpy(mulf_suffix, instr);
						mulf_dosuffix = YES;
					}
					break;
				case MULF_DFLT:
					mulf_dosuffix = YES;
					strcpy(mulf_suffix,".mul");
					mulf_master = 0;
					mulf_absolute = YES;
					break;
				case MULF_MAST:
					mulf_master = mulf_int[0];
					break;
				case MULF_ABSOLUTE:
					mulf_absolute = YES;
					break;
				case MULF_RELATIVE:
					mulf_absolute = NO;
					break;
			}
		}
	}
			
}


void gsac_exec_mulf(void)
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


	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
	if(mulf_master < 0 || mulf_master > ntrc){
		printf("mulf_master = %d not in range of [0,%d]\n",
				mulf_master,ntrc-1);
		return;
	}
	printf("MULF MASTER %d   \n", mulf_master );


	/* CHECK TO SEE THAT ALL DT ARE THE SAME ELSE RETURN 
		ALSO DETERMINE LATEST START TIME AND EARLIEST END TIME */
	/* note that we will use the idea of rotate3 to look at the
		common time window - unless we set relative in which 
		case we will ignore and use the shortest length */
	/* define the master trace parameters */
	nptsx = sacdata[mulf_master].sachdr.ihdr[H_NPTS];
	dtx = sacdata[mulf_master].sachdr.rhdr[H_DELTA];
	ts = sacdata[mulf_master].tzbeg;
	te = sacdata[mulf_master].tzend;
	tzref = sacdata[mulf_master].tzref;
	for ( k=0 ; k < ntrc ; k ++){
		dty = sacdata[k].sachdr.rhdr[H_DELTA];
		if(ABS (dtx - dty) > 0.01 * ABS(dtx)){
			printf("DT not equal - mulf not done\n");
			printf("Files %d %d with %f %f, respectively\n",
				mulf_master,k,dtx,dty);
			return;
		}
		if(sacdata[k].tzbeg > ts)ts = sacdata[k].tzbeg;
		if(sacdata[k].tzend < te)te = sacdata[k].tzend;
	}
        if(ts > te){
                printf("Cannot multiply traces - no trace overlap\n");
                return;
        }
	/* define the number of points in the output file */
	npts  = MIN(((int)((te -ts)/dtx + 0.49) + 1),nptsx);
	/* now cycle through the data set performing the operation 
		also reset the headers */

	/* compute offset for the master trace */
	noff_master = (int)((ts - sacdata[mulf_master].tzbeg)/dtx + 0.49);
	/* process by output point - this is inefficient in terms of 
		memory access but saves the effort of copying the 
		master trace since everything is
		overwritten in memory by this process */
	for(i=0;i< npts ; i++){
		x0 = sacdata[mulf_master].sac_data[i + noff_master];
		noff = (int)((ts - sacdata[k].tzbeg)/dtx + 0.49);
		for ( k=0 ; k < ntrc ; k ++){
			/* define offset for each trace */
			noff = (int)((ts - sacdata[k].tzbeg)/dtx + 0.49);
			y0 = sacdata[k].sac_data[i + noff];
			sacdata[k].sac_data[i] = y0*x0;
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
	
		chofname(sacdata[mulf_master].schdr[H_KSTNM],
			sacdata[mulf_master].schdr[H_KCMPNM],
			sacdata[k].schdr[H_KSTNM],
			sacdata[k].schdr[H_KCMPNM],
			sacdata[k].sac_ofile_name,
			mulf_dosuffix, mulf_suffix);

	}


	
}

