#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	OUTCSV_DFLT	0
#define	OUTCSV_X0		1
#define	OUTCSV_Y0		2
#define	OUTCSV_XLEN	3
#define	OUTCSV_YLEN	4
#define	OUTCSV_XLAB	5
#define	OUTCSV_YLAB	6


struct arghdr outcsvarg[] = {
	{OUTCSV_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float outcsv_real[10];
int   outcsv_int [10];
int   outcsv_yn;
int   outcsv_num;

static int *noff = (int *)NULL;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_outcsv(int ncmd, char **cmdstr)
{
	int i;
	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, outcsvarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; outcsvarg[i].key[0] != '\0' ; i++){
		if(outcsvarg[i].used > 0){
			if(outcsvarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, outcsvarg[i].key, 
					outcsvarg[i].mfit,outcsvarg[i].narg, outcsv_real);
			} else if(outcsvarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, outcsvarg[i].key, 
					outcsvarg[i].mfit,outcsvarg[i].narg, outcsv_int );
			} else if(outcsvarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, outcsvarg[i].key, 
					outcsvarg[i].mfit,outcsvarg[i].narg, &outcsv_yn );
			} else if(outcsvarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, outcsvarg[i].key, 
					outcsvarg[i].mfit,outcsvarg[i].narg, &outcsv_num );
			}
		}
	}
}

void gsac_exec_outcsv(void)
{
	/* output all traces in memory in CSV formal */
	int npts, ntrc;
	int nptsx;
	double ts, te;
	double tzref;
	float dtx, dty;
	int i, k;
	float toff;
	float x0;
	FILE *filecsv;
	int isxtime;


	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;

	/* CHECK TO SEE THAT ALL DT ARE THE SAME ELSE RETURN 
		ALSO DETERMINE LATEST START TIME AND EARLIEST END TIME */
	/* note that we will use the idea of rotate3 to look at the
		common time window - unless we set relative in which 
		case we will ignore and use the shortest length */
	/* define the master trace parameters */
	nptsx = sacdata[0].sachdr.ihdr[H_NPTS];
	dtx = sacdata[0].sachdr.rhdr[H_DELTA];
	ts = sacdata[0].tzbeg;
	te = sacdata[0].tzend;
	tzref = sacdata[0].tzref;
	if(sacdata[0].sachdr.ihdr[15] == 1){
		isxtime = YES ;
	} else {
		isxtime = NO ;
	}

	for ( k=0 ; k < ntrc ; k ++){
		dty = sacdata[k].sachdr.rhdr[H_DELTA];
		if(ABS (dtx - dty) > 0.01 * ABS(dtx)){
			printf("DT not equal - outcsv not done\n");
			printf("Files %d %d with %f %f, respectively\n",
				1,k,dtx,dty);
			return;
		}
		if(sacdata[k].tzbeg > ts)ts = sacdata[k].tzbeg;
		if(sacdata[k].tzend < te)te = sacdata[k].tzend;
	}
        if(ts > te){
                printf("Cannot output traces - no trace overlap\n");
                return;
        }
	/* define the number of points in the output file */
	npts  = MIN(((int)((te -ts)/dtx + 0.49) + 1),nptsx);
	/* now cycle through the data set performing the operation 
		also reset the headers */

	/* define offset for each trace */
	if(noff == (int *)NULL)
		noff = (int *)calloc(ntrc,sizeof(int));
	else
		noff = (int *)realloc(noff,ntrc*sizeof(int));
	/* get the offset */
	for(k=0;k< ntrc ; k++){
		noff[k] = (int)((ts - sacdata[k].tzbeg)/dtx + 0.49);
	}
	printf("Creating f001.csv\n");
	filecsv = fopen("f001.csv","w+");
	/* output the Column Titles */
	if(isxtime == NO )
		fprintf(filecsv,"Freq (Hz)");
	else
		fprintf(filecsv,"Time (sec)");
	for(k=0;k< ntrc ; k++){
		fprintf(filecsv,",%s%s",sacdata[k].schdr[H_KSTNM],
			sacdata[k].schdr[H_KCMPNM]);
	}
	fprintf(filecsv,"\n");
	
	/* now out the points */
	toff = sacdata[0].sachdr.rhdr[H_B] + noff[0]*dtx;
	for(i=0;i< npts ; i++){
		fprintf(filecsv,"%f",toff);
		for(k=0;k< ntrc ; k++){
			x0 = sacdata[k].sac_data[i + noff[k]];
			fprintf(filecsv,",%f",x0);
		}
		fprintf(filecsv,"\n");
		toff += dtx;
	}
	fclose(filecsv);
}
