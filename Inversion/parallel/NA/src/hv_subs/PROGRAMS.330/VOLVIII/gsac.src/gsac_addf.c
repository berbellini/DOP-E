#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include	"csstim.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	ADDF_DFLT	0
#define	ADDF_MAST	1
#define	ADDF_SUF	3
#define	ADDF_ABSOLUTE	4
#define	ADDF_RELATIVE	5

static int   addf_master = 0;
static int   addf_absolute = YES;
static char  addf_suffix[80];
static int   addf_dosuffix = NO;

struct arghdr addfarg[] = {
	{ADDF_DFLT, "DEFAULT", IHDR, 0, 0, NO, "",  1},
	{ADDF_SUF , "SUFFIX" , CHDR, 0, 1, NO, "SUFFIX suffix",  1},
	{ADDF_MAST, "MASTER" , IHDR, 0, 1, NO, "MASTER n",  1},
	{ADDF_ABSOLUTE, "ABSOLUTE", IHDR, 0, 0, NO, "",  1},
	{ADDF_RELATIVE, "RELATIVE", IHDR, 0, 0, NO, "",  1},
	{0,	""		, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float addf_real[10];
int   addf_int [10];
int   addf_yn;
int   addf_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_addf(int ncmd, char **cmdstr)
{
	int i;
	char instr[80];
	/* parsing code here */
	/* set up default */
	addf_dosuffix = YES;
	addf_suffix[0]='\0';
		strcpy(addf_suffix,".add");
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, addfarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; addfarg[i].key[0] != '\0' ; i++){
		if(addfarg[i].used > 0){
			if(addfarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, addfarg[i].key, 
					addfarg[i].mfit,addfarg[i].narg, addf_real);
			} else if(addfarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, addfarg[i].key, 
					addfarg[i].mfit,addfarg[i].narg, addf_int );
			} else if(addfarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, addfarg[i].key, 
					addfarg[i].mfit,addfarg[i].narg, &addf_yn );
			} else if(addfarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, addfarg[i].key, 
					addfarg[i].mfit,addfarg[i].narg, &addf_num );
			} else if(addfarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, addfarg[i].key, 
					addfarg[i].mfit, addfarg[i].narg, instr );
			}
			switch(addfarg[i].id){
				case ADDF_SUF:
					if(strlen(instr) < 80){
						strcpy(addf_suffix, instr);
						addf_dosuffix = YES;
					}
					break;
				case ADDF_DFLT:
					addf_dosuffix = YES;
					strcpy(addf_suffix,".mul");
					addf_master = 0;
					addf_absolute = YES;
					break;
				case ADDF_MAST:
					addf_master = addf_int[0];
					break;
				case ADDF_ABSOLUTE:
					addf_absolute = YES;
					break;
				case ADDF_RELATIVE:
					addf_absolute = NO;
					break;
			}
		}
	}
			
}


void gsac_exec_addf(void)
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
	if(addf_master < 0 || addf_master > ntrc){
		printf("addf_master = %d not in range of [0,%d]\n",
				addf_master,ntrc-1);
		return;
	}
	printf("ADDF MASTER %d   \n", addf_master );


	/* CHECK TO SEE THAT ALL DT ARE THE SAME ELSE RETURN 
		ALSO DETERMINE LATEST START TIME AND EARLIEST END TIME */
	/* note that we will use the idea of rotate3 to look at the
		common time window - unless we set relative in which 
		case we will ignore and use the shortest length */
	/* define the master trace parameters */
	nptsx = sacdata[addf_master].sachdr.ihdr[H_NPTS];
	dtx = sacdata[addf_master].sachdr.rhdr[H_DELTA];
	ts = sacdata[addf_master].tzbeg;
	te = sacdata[addf_master].tzend;
	tzref = sacdata[addf_master].tzref;
	for ( k=0 ; k < ntrc ; k ++){
		dty = sacdata[k].sachdr.rhdr[H_DELTA];
		if(ABS (dtx - dty) > 0.01 * ABS(dtx)){
			printf("DT not equal - addf not done\n");
			printf("Files %d %d with %f %f, respectively\n",
				addf_master,k,dtx,dty);
			return;
		}
		if(sacdata[k].tzbeg > ts)ts = sacdata[k].tzbeg;
		if(sacdata[k].tzend < te)te = sacdata[k].tzend;
	}
        if(ts > te){
                printf("Cannot add traces - no trace overlap\n");
                return;
        }
	/* define the number of points in the output file */
	npts  = MIN(((int)((te -ts)/dtx + 0.49) + 1),nptsx);
	/* now cycle through the data set performing the operation 
		also reset the headers */

	/* compute offset for the master trace */
	noff_master = (int)((ts - sacdata[addf_master].tzbeg)/dtx + 0.49);
	/* process by output point - this is inefficient in terms of 
		memory access but saves the effort of copying the 
		master trace since everything is
		overwritten in memory by this process */
	for(i=0;i< npts ; i++){
		x0 = sacdata[addf_master].sac_data[i + noff_master];
		noff = (int)((ts - sacdata[k].tzbeg)/dtx + 0.49);
		for ( k=0 ; k < ntrc ; k ++){
			/* define offset for each trace */
			noff = (int)((ts - sacdata[k].tzbeg)/dtx + 0.49);
			y0 = sacdata[k].sac_data[i + noff];
			sacdata[k].sac_data[i] = y0+x0;
		}
	}
	/* now update headers and output file name */
	for(k = 0 ; k < ntrc ; k++){
		getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
		sacdata[k].sachdr.ihdr[H_NPTS] = npts;
		sacdata[k].sachdr.rhdr[H_B] = ts - tzref;
		sacdata[k].sachdr.rhdr[H_E] = te - tzref;
		sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
		sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
		sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		/* update the headers */
		etoh(tzref, &sacdata[k].sachdr.ihdr[H_NZYEAR], 
			&sacdata[k].sachdr.ihdr[H_NZJDAY], &month, &day,
			&sacdata[k].sachdr.ihdr[H_NZHOUR], 
			&sacdata[k].sachdr.ihdr[H_NZMIN],
			&sacdata[k].sachdr.ihdr[H_NZSEC], 
			&sacdata[k].sachdr.ihdr[H_NZMSEC]);
	
		chofname(sacdata[addf_master].schdr[H_KSTNM],
			sacdata[addf_master].schdr[H_KCMPNM],
			sacdata[k].schdr[H_KSTNM],
			sacdata[k].schdr[H_KCMPNM],
			sacdata[k].sac_ofile_name,
			addf_dosuffix, addf_suffix);

	}


	
}

