#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include	"csstim.h"

/* CHANGES
	16 SEP 2009 - improper index on about line 214 to change
	sacdata[k].tzbegx = sacdata[k].tzbeg;
	sacdata[k].tzendx = sacdata[k].tzend;

	sacdata[0].tzbegx = sacdata[0].tzbeg;
	sacdata[0].tzendx = sacdata[0].tzend;
*/

extern struct sacfile_ *sacdata;
extern int *sortptr;

#define MERGE_DFLT      0
#define MERGE_GAP       1

/* these are temporary variables only used here */
static float merge_real[10];

static int   merge_gap_fill;


struct arghdr mergearg[] = {
	{MERGE_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{MERGE_GAP  , "GAP"  , RHDR, NO, 1, NO, "GAP gap_fill", 1},
	{0,	""		, IHDR, NO, 0, NO, "", -1}
};

static float *xa = (float *)NULL;
static int   *xl = (int   *)NULL;
static char mkstnm [9];
static char mkcmpnm[9];

static float  fhdr_default = -12345.;

/*                  B  E  O  A  T0  T1 ............................ T9 */
static int tm[] = { 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };

void gsac_set_param_merge(int ncmd, char **cmdstr)
{
	int i;
	/* initial debug */
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, mergearg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; mergearg[i].key[0] != '\0' ; i++){
		if(mergearg[i].used > 0){
			if(mergearg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, mergearg[i].key,
					mergearg[i].mfit,mergearg[i].narg, merge_real);
}
			switch(mergearg[i].id){
				case MERGE_DFLT:
					merge_gap_fill = 0.0;
					break;
				case MERGE_GAP:
					merge_gap_fill = merge_real[0];
			}
		}
	}
		
}

void gsac_exec_merge(void)
{
char ostr1[100];
char ostr2[100];
	/* 
	 * 1) look at all traces to determine the extreme time window
	 * 2) create merged trace array and find a flag for array set
	 * 3) fill the arrays and not exceptions at first and last points
	 * 4) create the merged array - also do a QC do denote any gaps
	      */
	int k, ntrc, i, l , j, npts;
	int j1, j2;
	double tzmin, tzmax, tzb, tze;
	int newpts;
	float dtime;
	float delta, tdelta;
	float depmax, depmin, depmen;
	int indmax, indmin;
	int month, day;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
	gsac_control.begmin= 1.0e+30;
	gsac_control.endmax=-1.0e+30;
	/* since synchronize will not change the limits of the
	 * earliest and latest sample of points in epoch time,
	 * we will not change the gsac_control.begmin or
	 * gsac_control.endmax here */
	/* phase 1 - get the minimum B time */
	for ( k=0 ; k < ntrc ; k ++){
		tzb = sacdata[k].tzref + sacdata[k].sachdr.rhdr[H_B];
		tze = sacdata[k].tzref + sacdata[k].sachdr.rhdr[H_E];
		printtimestr(tzb,ostr1);printtimestr(tze,ostr2);
		printf("%s %s (%d)\n",ostr1,ostr2,k);
		tdelta = sacdata[k].sachdr.rhdr[H_DELTA];
		if(k == 0){
		 	tzmin = tzb;
		 	tzmax = tze;
			delta = tdelta;
			strcpy(mkstnm , sacdata[k].schdr[H_KSTNM] );
			strcpy(mkcmpnm, sacdata[k].schdr[H_KCMPNM]);
		} else {
			if(tzb < tzmin)tzmin = tzb;
			if(tze > tzmax)tzmax = tze;
			if(ABS(delta - tdelta) > 0.25*delta){
				printf("Fail: DT different %f vs %f\n",delta,tdelta);
				return;
			}
			if(strcmp(mkstnm,sacdata[k].schdr[H_KSTNM]) !=0){
				printf("Fail: STATION name different %s vs %s\n",
				mkstnm,sacdata[k].schdr[H_KSTNM]);
					return;
			}
			if(strcmp(mkcmpnm,sacdata[k].schdr[H_KCMPNM]) !=0){
				printf("Fail: COMPONENT name different %s vs %s\n",
				mkstnm,sacdata[k].schdr[H_KCMPNM]);
					return;
			}
		}
	}
	printtimestr(tzmin,ostr1);printtimestr(tzmax,ostr2);
	printf("\n%s %s (Merge window)\n",ostr1,ostr2);
	/* we now have the extreme time window define the new arrays
	 * and intitialize. The xl array is 0 if point is now set */
	newpts = 1 + (int)((tzmax - tzmin)/delta + 0.49);
	printf("New time series length:  %d\n",newpts);
	if(xa == (float *)NULL)
		xa = (float *)calloc(newpts,sizeof(float));
	else
		xa = (float *)realloc(xa,newpts*sizeof(float));
	if(xl == (int *)NULL)
		xl = (int *)calloc(newpts,sizeof(int));
	else
		xl = (int *)realloc(xl,newpts*sizeof(int));
	for(i=0;i < newpts ; i++){
		xa[i] = merge_gap_fill;
		xl[i] = 0  ;
	}
	/* preserve any time picks even though why would one make these
	 * with a partial trace */


	for ( k=0 ; k < ntrc ; k ++){
		tzb = sacdata[k].tzref + sacdata[k].sachdr.rhdr[H_B];
		dtime = (float)(tzb - tzmin);
		for(i= 0 ; i < 14 ; i ++){
			/* B and E always reset */
			if ( i < 2 )
				sacdata[k].sachdr.rhdr[tm[i]] += dtime;
			else
				if(sacdata[k].sachdr.rhdr[tm[i]] != fhdr_default)
				sacdata[k].sachdr.rhdr[tm[i]] += dtime;

		}
		sacdata[k].tzref = tzmin;
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		j = (int)((tzb - tzmin)/delta + 0.49);
printf("k %d j %d - %d  npts %d\n",k,j,j+npts -1,npts);
		for(i = j, l=0 ; i < j+npts ; i++, l++){
			/* put in an exception check here */
			if(xl[i] != 0 && (i == j || i == j+npts -1 )  ){
				if (xa[i] != sacdata[k].sac_data[l] ){
					printf("(%d) warning just for 1st and last\n",k);
					if(i>0)printf("xa[%d] = %f\n",i-1,xa[i-1]);
					printf("xa[%d] = %f != sacdata[%d].sac_data[%d] = %f\n",i,xa[i],k,l,sacdata[k].sac_data[l]);
					if(i< newpts-2)printf("xa[%d] = %f\n",i+1,xa[i+1]);
				}
			}
			xa[i] = sacdata[k].sac_data[l];
			xl[i] = 1 ;
		}
	}
	/* now do quality control */
	j1 = 0;
	for(i=0 ; i < newpts ; i++){
		if(xl[i] == 0 && j1 == 0){
			printf("Data gap begins at sample %d\n",i);
			j1 = 1;
		}
		if(j1 == 1 && xl[i] == 1){
			printf("Data gap ends   at sample %d\n",i-1);
			j1 = 0;
		}
	}
		/* now redefine the parameters of the first trace */

		sacdata[0].tzref = tzmin;
		sacdata[0].sachdr.ihdr[H_NPTS] = newpts;
		sacdata[0].sac_data = realloc(sacdata[0].sac_data,newpts*sizeof(float));
		sacdata[0].sachdr.rhdr[H_B] = 0;
		sacdata[0].sachdr.rhdr[H_E] = (newpts -1 ) * delta;
		/* copy array in to save */
		for(i=0 ; i < newpts ; i++){
			sacdata[0].sac_data[i] = xa[i];
		}
		getmxmn(sacdata[0].sac_data, newpts,&depmax, &depmin, &depmen,&indmax,&indmin);
		sacdata[0].sachdr.rhdr[H_TIMMAX] = sacdata[0].sachdr.rhdr[H_B]  + ( indmax)*sacdata[0].sachdr.rhdr[H_DELTA] ;
		sacdata[0].sachdr.rhdr[H_TIMMIN] = sacdata[0].sachdr.rhdr[H_B]  + ( indmin)*sacdata[0].sachdr.rhdr[H_DELTA] ;
		sacdata[0].sachdr.rhdr[H_DEPMIN] = depmin;
		sacdata[0].sachdr.rhdr[H_DEPMEN] = depmen;
		sacdata[0].sachdr.rhdr[H_DEPMAX] = depmax;

		sacdata[0].tzbeg = tzmin;
		sacdata[0].tzend = tzmax;
		sacdata[0].tzbegx = sacdata[0].tzbeg;
		sacdata[0].tzendx = sacdata[0].tzend;
		/* get bounds for absolute plotting */
		gsac_control.begmin = sacdata[0].tzbeg;
		gsac_control.endmax = sacdata[0].tzend;

		/* NOW RESET THE OTHER TIME FIELDS */
		etoh(tzmin, &sacdata[0].sachdr.ihdr[H_NZYEAR], 
			&sacdata[0].sachdr.ihdr[H_NZJDAY], &month, &day,
			&sacdata[0].sachdr.ihdr[H_NZHOUR], 
			&sacdata[0].sachdr.ihdr[H_NZMIN],
			&sacdata[0].sachdr.ihdr[H_NZSEC], 
			&sacdata[0].sachdr.ihdr[H_NZMSEC]);
		/* redefine component names and default file name */
		j1 = findblank(sacdata[0].schdr[H_KSTNM] );
		j2 = findblank(sacdata[0].schdr[H_KCMPNM]);
		sacdata[0].sac_ofile_name[0] = '\0';
		strncpy(sacdata[0].sac_ofile_name,sacdata[0].schdr[H_KSTNM ],j1);
		sacdata[0].sac_ofile_name[j1] = '\0';
		strncat(sacdata[0].sac_ofile_name,sacdata[0].schdr[H_KCMPNM],j2);
		sacdata[0].sac_ofile_name[j1+j2] = '\0';
		/* redefine the number of traces to 1 */
		printf("New default output filename: %s\n",sacdata[0].sac_ofile_name);
		gsac_control.number_otraces = 1 ;
}
