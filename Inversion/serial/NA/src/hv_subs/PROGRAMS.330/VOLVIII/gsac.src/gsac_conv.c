/* Changes:
	24 OCT 2010 - corrected the indexing which was a problem
		when the interpolated series of the filter had fewer points than
                the files in memory  of the filter. This also pointed out
                a problem with gsac_in 
*/
#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	CONV_DFLT	0
#define	CONV_ZERO 	1
#define	CONV_FILE	2
#define	CONV_SUF	3

static char conv_suffix[80];
static int conv_dosuffix;
static char *conv_filter_file;
static int conv_zero_marker;
static char conv_zero_string[5];
static int conv_do_convolve;
static float *sy = (float *)NULL;
static double *sx = (double *)NULL;
static double *x = (double *)NULL;
static float *h = (float *)NULL;
static float *y = (float *)NULL;

extern void gsac_inter(double *sx,float *sy,double *x,float *y,int npts, int *mpts);




struct arghdr convarg[] = {
	{CONV_DFLT , "DEFAULT", IHDR, 0, 0, NO, "", -1},
	{CONV_ZERO , "ZERO "  , CHDR, 0, 1, NO, "Zero [O|A|T0|...|T9]", 1},
	{CONV_FILE , "FILE"   , CHDR, 0, 1, NO, "File filter_file", 1},
	{CONV_SUF  , "SUFFIX" , CHDR, 0, 1, NO, "Suffix suffix",  1},
	{0,	""		, IHDR, 0, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float conv_real[10];
int   conv_int [10];
int   conv_yn;
int   conv_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_conv(int ncmd, char **cmdstr)
{
	int i,ls;
	char instr[1000];

	/* set up default */
	conv_dosuffix = YES;
        conv_zero_marker = H_O;	
	conv_do_convolve = NO;
	conv_filter_file = (char *)NULL;
	conv_suffix[0]= '\0';
		strcpy(conv_suffix,".con");

	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, convarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; convarg[i].key[0] != '\0' ; i++){
		if(convarg[i].used > 0){
			if(convarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, convarg[i].key, 
					convarg[i].mfit,convarg[i].narg, conv_real);
			} else if(convarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, convarg[i].key, 
					convarg[i].mfit,convarg[i].narg, conv_int );
			} else if(convarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, convarg[i].key, 
					convarg[i].mfit,convarg[i].narg, &conv_yn );
			} else if(convarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, convarg[i].key, 
					convarg[i].mfit,convarg[i].narg, &conv_num );
			} else if(convarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, convarg[i].key, 
					convarg[i].mfit,convarg[i].narg, instr );
			}
			switch(convarg[i].id){
				case CONV_SUF:
					if(strlen(instr) < 80){
						strcpy(conv_suffix, instr);
						conv_dosuffix = YES;
					}
					break;
				case CONV_ZERO:
					/* a lengthy comparison */
					gsac_strupr(instr);
					strncpy(conv_zero_string,instr,4);
					conv_zero_string[4]= '\0';
					if(strncmp("O",instr,1)==0 )
						conv_zero_marker = H_O;
					else if(strncmp("A",instr,1)==0 )
						conv_zero_marker = H_A;
					else if(strncmp("T0",instr,2)==0 )
						conv_zero_marker = H_T0;
					else if(strncmp("T1",instr,2)==0 )
						conv_zero_marker = H_T1;
					else if(strncmp("T2",instr,2)==0 )
						conv_zero_marker = H_T2;
					else if(strncmp("T3",instr,2)==0 )
						conv_zero_marker = H_T3;
					else if(strncmp("T4",instr,2)==0 )
						conv_zero_marker = H_T4;
					else if(strncmp("T5",instr,2)==0 )
						conv_zero_marker = H_T5;
					else if(strncmp("T6",instr,2)==0 )
						conv_zero_marker = H_T6;
					else if(strncmp("T7",instr,2)==0 )
						conv_zero_marker = H_T7;
					else if(strncmp("T8",instr,2)==0 )
						conv_zero_marker = H_T8;
					else if(strncmp("T9",instr,2)==0 )
						conv_zero_marker = H_T9;
					else 
						conv_zero_marker = -1;
					break;
				case CONV_FILE:
					ls = strlen(instr);
					conv_filter_file = calloc(ls+1,sizeof(char));
					strcpy(conv_filter_file,instr);
					conv_do_convolve = YES;
					break;

			}
		}
	}
}
			
		
struct sacfile_ *convdata ;

void gsac_exec_conv(void)
{
	int i, j, k, ntrc, iret;
	int hnpts;
	float hdelta;
	int npts;
	int newhpts;
	float delta;
	float zero_marker_offset;
	float depmax, depmin, depmen;
	int indmax, indmin;
	double sum;
	int zero_offset;
	float tmpx;
/*
printf("conv_do_convolve %d\nb",conv_do_convolve);
printf("suffix           %s\n",conv_suffix);
printf("conv_zero_marker %d\n",conv_zero_marker);
printf("conv_filter_file %s\n",conv_filter_file);
*/

	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
	/* begin extensive error checks */
	/* check to see if the file is a valid sac file */
	if(conv_do_convolve == NO){
		printf("CONV Fail: Filter file not specified\n");
		return;
	}
	if(gsac_valid_sacfile(conv_filter_file) == 0){
		printf("CONV Fail: Not a valid sac file: %s\n",conv_filter_file);
		free(conv_filter_file);
		return;
	}
	/* allocate a data structure for the external file */
	convdata=(struct sacfile_ *)calloc(1, sizeof(struct sacfile_));
	convdata[0].sac_data = (float *)NULL;
	/* safety check on the external file */
	iret = brsac(conv_filter_file,&convdata[0].sachdr,
		&convdata[0].sac_data);
	if(iret < 0 ){
		printf("CONV Fail: corrupt sac file \n");
                       free(convdata);
		return;
	}
	/* check to see if the zero marker is set */
	if( convdata[0].sachdr.rhdr[conv_zero_marker] == -12345.){
		printf("CONV Fail: Zero marker %s not set in file %s\n",
			 conv_zero_string,conv_filter_file);
                       free(convdata);
		return;
	}
	/* determine sampling parameters of the external file */
	hnpts = convdata[0].sachdr.ihdr[H_NPTS];
	hdelta = convdata[0].sachdr.rhdr[H_DELTA];
	/* get the position of the zero marker in the file */
	zero_marker_offset = 
			(convdata[0].sachdr.rhdr[conv_zero_marker]-convdata[0].sachdr.rhdr[H_B])/hdelta;
	/* copy the external data into a temporary array 
		for resampling using gsac_inter */
	sx = (double *)calloc(hnpts,sizeof(double));
	sy = (float *)calloc(hnpts,sizeof(float));
	/* the (sx,sy) pair corresponds to the external 
		trace in time with respect to B */
	for(i=0;i< hnpts;i++){
		sx[i] = (double)i*(double)hdelta;
		sy[i] = convdata[0].sac_data[i];
	}

	/* process the traces */

	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		delta = sacdata[k].sachdr.rhdr[H_DELTA];
		/* resample the external trace at the delta of the
			trace in memory. Also ensure that the sampling
			is with respect tot he zero marker */
		newhpts = (int)( hnpts*hdelta/delta);
		/* define the sampling for the interpolated trace */
		h = (float *)calloc(newhpts,sizeof(float));
		x = (double *)calloc(newhpts,sizeof(double));
                /* the output (x,h) pairs correspond to external trace
                   resampled with new delta 
			with respect to the zero_marker of external trace */
		tmpx = (convdata[0].sachdr.rhdr[conv_zero_marker]
				-convdata[0].sachdr.rhdr[H_B]) ; 
		zero_offset = tmpx/delta ;
		/* offset of resampled */
		for(i=0; i < newhpts; i++){
			x[i] = (double)i*(double)delta + 
				tmpx - zero_offset*delta;
		}
		/* now interpolate the external trace */
		gsac_inter(sx,sy,x,h,hnpts, &newhpts);
		/* allocate space for the convolution output and
			initialize to zero */
		y = (float *)calloc(npts,sizeof(float));
		/* now do the convolution */
			/* begin convolution 
				y[i] = SUM x[j] h[i+o-j] for 0 <=i < N
                                and x is defined for 0, ..., N-1
                                and h is defined for 0, ..., M-1
                                Here o is the offset on h
                           the j summation is constrained by
                                Limits on x:   0 <= j < N 
                                Limits on h:   0 <= i+o-j < M
                                   or
                                               j > i+o-M which is same as
                                               j>= i+o-M-1
                                   and
                                               j <= i+o which is the same as
                                               j < i+o+1
                           putting together j>= MAX(0, i+o-M+1) and
					    j < MIN(N,i+o+1)
			*/
		for(i=0;i < npts ; i++){
			sum = 0.0;
			/*
			if(i + zero_offset < newhpts){
				for(j=0 ; i >= j - zero_offset  && j < newhpts ; j++){
					sum += sacdata[k].sac_data[j] * h[i + zero_offset -j];
				}
			}
			*/
				for(j=MAX(0, i+zero_offset-newhpts-1)  ;  j < MIN(npts,i+zero_offset+1) ; j++){
					sum += sacdata[k].sac_data[j] * h[i + zero_offset -j];
				}
                        y[i] = sum;
		}
		/* multiply by dt to finish the convolution
			and overwrite the data array in momory */
		for(i=0 ; i < npts ; i++){
			sacdata[k].sac_data[i] = y[i]*delta;
		}
		/* update the header values for the k th trace */
		getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
		sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
		sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
		sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
		sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		/* update the output file name */
		strcat(sacdata[k].sac_ofile_name,conv_suffix);
		/* free temporary arrays defined for this external trace */
		free(x);
		free(h);
		free(y);

	}
	/* clean up */
	free(sx);
	free(sy);
	free(conv_filter_file);
	free(convdata[0].sac_data);
	free(convdata);
}

